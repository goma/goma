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
 *$Id: mm_input_bc.c,v 5.21 2010-07-21 16:39:27 hkmoffa Exp $
 */

#include <exodusII.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "el_elm.h"
#include "md_timer.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_bc.h"
#include "mm_eh.h"
#include "mm_input.h"
#include "mm_interface.h"
#include "mm_mp_const.h"
#include "mm_post_def.h"
#include "mm_species.h"
#include "rf_allo.h"
#include "rf_bc_const.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_io.h"
#include "rf_io_const.h"
#include "rf_mp.h"
#include "rf_solver.h"
#include "rf_solver_const.h"
#include "std.h"

/*
 * Hey! This is the *one* place where these are defined. All other locations
 * have a mm_mp_structs and mm_mp.h to declare what these are.
 */

extern int Num_Var_Init_Mat[MAX_NUMBER_MATLS]; /* number of variables to overwrite  *
                                                * with material-specific            *
                                                * initialization                    */

extern struct Boundary_Condition *BC_Types;

extern struct Rotation_Specs *ROT_Types;

static Spfrtn sr;

/*
 * What to look for each time...
 */

/*************** R O U T I N E S   I N   T H E   F I L E ***********************
 *
 *    NAME				TYPE		CALLED_BY
 *--------------------------------------------------------------------
 *
 *    rd_bc_specs()			void		read_input_file
 *
 */

static int BC_consistency(struct Boundary_Condition *); /* BC_Type                           */

static int detect_BC(int, int);

/*
 *	This file is a break-off from the very large file mm_input.c.
 *      See the comments found therein.
 *
 *
 *	Author:			Edward D. Wilkes, GRAM Inc.
 *	Revised:
 */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 * rd_bc_specs -- read input file for boundary condition specifications
 *
 * Comments:	This code was lifted out of the ever-growing read_input_file()
 *           	routine above. This makes things more modular, with division
 *           	into "sections".
 *
 *		Someday we need to comb through all these placeholder options
 *		to see if they are useful.
 *
 * Revised:			Fri Oct 29 10:15:45 MDT 1993 pasacki@sandia.gov
 */

void rd_bc_specs(FILE *ifp, char *input) {
  char err_msg[MAX_CHAR_IN_INPUT];
  int i;
  int j;
  //  int error;
  int ibc, var, num_const;
  char ts[MAX_LINE_LENGTH];

  char seed_string[80];
  char instruction_string[DIM][80];
  char component_string[DIM][120];
  char condition_string[DIM][80];
  char coordinate_string[DIM][80] = {"x", "y", "z"};
  char topo_string[80];
  char rot_eq_string[80];

  static const char yo[] = "rd_bc_specs";
  char echo_string[MAX_CHAR_IN_INPUT] = "\0";
  char *echo_file = Echo_Input_File;

  int overlap_bc = FALSE;
  int iread;
  int NO_SPECIES = -1;
  int k, eq_found, eqn, irc, p;
  dbl sx, sy, sz;

  struct Rotation_Specs *rot;

  /*
   * Initialize...
   */

  num_new_BC_Desc = 0;

  /*
   * Read boundary condition specifications
   */

  iread = look_for_optional(ifp, "Boundary Condition Specifications", input, '=');

  look_for(ifp, "Number of BC", input, '=');
  if (fscanf(ifp, "%d", &Num_BC) != 1) {
    fprintf(stderr, "error reading number of boundary conditions try to read list");
    Num_BC = -1;
  }
  /*
   * Count boundary conditions if Num_BC is set to -1
   */
  if (Num_BC == -1) {
    Num_BC = count_list(ifp, "BC", input, '=', "END OF BC");
  }

  /*
   *  Allocate space for the vector, BC_Types.
   */
  BC_Types = alloc_struct_1(struct Boundary_Condition, Num_BC);

  if (Debug_Flag && ProcID == 0) {
    printf("%s:\tallocated %d copies of the %lu sized\n", yo, Num_BC,
           (unsigned long int)sizeof(struct Boundary_Condition));
    printf("%s:\tBoundary_Condition structure.\n", yo);
  }

  ECHO("\n***Boundary Condition Specifications***\n", echo_file);

  SPF(echo_string, "%s = %d", "Number of BC", Num_BC);
  ECHO(echo_string, echo_file);

  /* initialize number of interface sources*/
  Num_Interface_Srcs = 0;
  for (ibc = 0; ibc < Num_BC; ibc++) {

    /*
     *   Initialize the boundary condition to its
     *   default state.
     */
    initialize_Boundary_Condition(BC_Types + ibc);

    look_for(ifp, "BC", input, '=');
    if (fscanf(ifp, "%s", ts) != 1) {
      GOMA_EH(-1, "error reading BC_Types[ibc].BC_Name");
    }

    /*
     *  Set the BC_Name based upon the string read into
     *  the character variable ts.
     *
     *  -> Also, assign the pointer to the BC_Description structure
     *     that matches the string, here.
     */
    BC_Types[ibc].desc = NULL;
    for (k = 0; k < Num_BC_Names; k++) {
      if (!strcmp(ts, BC_Desc[k].name1) || !strcmp(ts, BC_Desc[k].name2)) {
        BC_Types[ibc].BC_Name = BC_Desc[k].BC_Name;
        BC_Types[ibc].desc = &BC_Desc[k];
        BC_Types[ibc].BC_Desc_index = k;
      }
    }

    if (BC_Types[ibc].desc == NULL) {
      fprintf(stderr, "%s:\tBC %s not recognized.\n", yo, ts);
      exit(-1);
    }

    /*
     * Read the boundary condition Set_type.
     *
     * There are five valid strings,
     * N* and S* for node sets and side sets and LS for level set.
     *
     * NS, SS: No continuation
     * NC, SC: Continuation
     *
     */
    if (fscanf(ifp, "%80s", input) != 1) {
      fprintf(stderr, "%s:\tError reading BC_Types[ibc].Set_Type\n", yo);
      exit(-1);
    }

    if (!strcmp(input, "NS")) {
      (void)strncpy(BC_Types[ibc].Set_Type, "NS", 3);
    } else if (!strcmp(input, "NC")) {
      (void)strncpy(BC_Types[ibc].Set_Type, "NS", 3);
      cont->upBCID = ibc;
      printf("BC Continuation: BCID = %4d\n", cont->upBCID);
    } else if (!strcmp(input, "SS")) {
      (void)strncpy(BC_Types[ibc].Set_Type, "SS", 3);
    } else if (!strcmp(input, "SC")) {
      (void)strncpy(BC_Types[ibc].Set_Type, "SS", 3);
      cont->upBCID = ibc;
      printf("BC Continuation: BCID = %4d\n", cont->upBCID);
    } else if (!strcmp(input, "LS")) {
      (void)strncpy(BC_Types[ibc].Set_Type, "LS", 3);

      if (!strcmp(BC_Types[ibc].desc->name1, "BAAIJENS_FLUID_SOLID") ||
          !strcmp(BC_Types[ibc].desc->name1, "LS_NO_SLIP")) {
        /* bad noble hack
        GOMA_EH(GOMA_ERROR,"Cannot apply BAAIJENS_SOLID_FLUID or LS_NO_SLIP to set type LS. Use
        PF");
        */
      }
    } else if (!strcmp(input, "PF")) {
      (void)strncpy(BC_Types[ibc].Set_Type, "PF", 3);
    } else {
      fprintf(stderr, "%s:\tunrecognized Set_Type for current boundary condition %d - %s\n", yo,
              ibc, input);
      exit(-1);
    }

    /*
     * Read in the ID of the boundary condition.  This ID refers to an ID in
     * the EXODUS binary input file.  The existence of this ID will be checked
     * later...
     */

    if (fscanf(ifp, "%d", &BC_Types[ibc].BC_ID) != 1) {
      fprintf(stderr, "%s:\tError reading BC_Types[ibc].BC_ID\n", yo);
      exit(-1);
    }

    /*
     * Summarize BC's that were successfully parsed so far...(argument
     * processing still to come)...
     */
    for (k = 0; k < Num_BC_Names; k++) {
      if (!strcmp(ts, BC_Desc[k].name1) || !strcmp(ts, BC_Desc[k].name2)) {

        SPF(echo_string, "BC = %s %s %d", BC_Desc[k].name1, BC_Types[ibc].Set_Type,
            BC_Types[ibc].BC_ID);
      }
    }

    /* Check for fluid/solid interaction BC's (embedded or contact) */
    if (BC_Types[ibc].desc->method == EMBEDDED_SURF || BC_Types[ibc].desc->method == CONTACT_SURF)
      overlap_bc = TRUE;

    /*
     * Read in a second SS number for EDGE conditions
     */
    if (BC_Types[ibc].desc->method == COLLOCATE_EDGE ||
        BC_Types[ibc].desc->method == WEAK_INT_EDGE ||
        BC_Types[ibc].desc->method == STRONG_INT_EDGE) {
      if (fscanf(ifp, "%d", &BC_Types[ibc].BC_ID2) != 1) {
        sprintf(err_msg, "Expected 2nd int for SSID for EDGE on %s.", BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_ID2);
    }

    /*
     * Here's a RECIPE for adding new boundary conditions so you don't have any
     * excuses not to add new ones.  The changes sould be made in at least
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
     *  Boundary Condition  Step 2: add a description of your boundary condition
     *                              to the BC_Desc structure array in mm_names.h.
     *      mm_names.h              This structure includes information about the
     *                              type of boundary condition, which equation it
     *                              applies to, what variables it is sensitive to,
     *                              whether to rotate mesh or momentum equations,
     *                              etc.  It is very important that you fill out
     *                              this structure carefully, otherwise the code
     *                              won't know what to do.
     *  Boundary Condition  Step 3: add your BC case to the correct format listing
     *                              for reading the necessary arguments from the
     *      mm_input_bc.c           input file in mm_input_bc.c (this file).
     *
     *  Boundary Condition  Step 4: Add a function call (and a function) in the
     *                              correct routine for evaluating your boundary
     *      bc_colloc.c             condition.  This will probably in bc_colloc.c
     *      bc_integ.c              for collocated conditions or bc_integ.c for
     *                              strong or weak integrated conditions.
     *  Boundary Conditions Step 5: Add type into BC conflict lists.
     *  Boundary Condition  Step 5: use and enjoy your new boundary condition
     *
     * Step 3 should be done below:
     */
    /*
     * Read in data depending on the type of boundary condition...
     */

    switch (BC_Types[ibc].BC_Name) {
      /* Fall through for all cases which don't require any further
       * data input...
       */

    case LS_SOLID_FLUID_BC:
      overlap_bc = TRUE;
    case PSPG_BC:
    case KIN_ELECTRODEPOSITION_BC:   /*  RSL 5/27/02  */
    case VNORM_ELECTRODEPOSITION_BC: /*  RSL 5/30/02  */
    case Q_VELO_SLIP_BC:
    case LS_NO_SLIP_BC:
    case LS_CAPILLARY_BC:
    case LS_CAP_CURVE_BC:
    case LS_CAP_DIV_N_BC:
    case LS_CAP_DIV_S_N_BC:
    case H_FREE_BC:
    case QNOBC_BC:
    case APR_NOBC_BC:
    case API_NOBC_BC:
    case POTENTIAL_NOBC_BC:
    case LS_INLET_BC:
    case LS_CONT_T_BC:
    case LS_CONT_FLUX_BC:
    case LS_CONT_VEL_BC:
    case LS_CONT_TRACTION_BC:
    case SHELL_SURFACE_CHARGE_BC:
    case SHELL_SURFACE_CHARGE_SIC_BC:
    case SHELL_DIFF_KINEMATIC_BC:
    case LS_EXTV_FLUID_SIC_BC:
    case LS_EIK_KINEMATIC_BC:
    case LS_EXTV_KINEMATIC_BC:
    case LS_EXTV_KIN_LEAK_BC:
    case REP_FORCE_SHU_BC:
    case REP_FORCE_SHU_SIC_BC:
    case SHELL_GRAD_FP_NOBC_BC:
    case SHELL_FLOW_DEVELOPED_BC:
    case SHELL_GRAD_FH_NOBC_BC:
    case SHELL_GRAD_PC_NOBC_BC:
    case STRESS_DEVELOPED_BC:
    case SHELL_TFMP_FREE_LIQ_BC:
    case SHELL_TFMP_NUM_DIFF_BC:
    case SHELL_TFMP_GRAD_S_BC:
    case SHELL_TFMP_FREE_GAS_BC:
    case SHELL_LUBRICATION_OUTFLOW_BC:
    case ZERO_VELO_TANGENT_3D_BC:
    case EM_ER_FREE_BC:
    case EM_EI_FREE_BC:
    case EM_HR_FREE_BC:
    case EM_HI_FREE_BC:
    case E_ER_2D_BC:
    case E_EI_2D_BC:
    case RESTIME_NOBC_BC:
    case EM_MMS_SIDE_BC:
    case EM_MMS_SIDE_IMAG_BC:
    case EM_ABSORBING_REAL_BC:
    case EM_ABSORBING_IMAG_BC:

      break;

      /* Fall through for all cases which requires a single integer value
       * as data input
       */
    case MOVING_PLANE_ETCH_BC:

      if (fscanf(ifp, "%d", &BC_Types[ibc].BC_Data_Int[0]) != 1) {
        sr = sprintf(err_msg, "%s: Expected 1 int for %s.", yo, BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_Data_Int[0]);

      break;

      /* Fall through for all cases which require a single floating point
       * value as data input and one additional optional parameter
       */
    case SURFTANG_SCALAR_BC:
      GOMA_WH(-1, "Use CAP_ENDFORCE_SCALAR for consistent sign convention");
      /*FALLTHROUGH*/
    case QSIDE_BC:
      /*	case TNRMLSIDE_BC:  12/6/01 */
      /*	case TSHRSIDE_BC:   ..dal   */
    case CURRENT_BC:
    case CURRENT_SIC_BC:
    case KINEMATIC_BC:
    case LUB_KINEMATIC_BC:
    case KIN_LEAK_HEAT_BC:
    case KINEMATIC_PETROV_BC:
    case KINEMATIC_COLLOC_BC:
    case KINEMATIC_DISC_BC:
    case KINEMATIC_EDGE_BC:
    case T_MELT_BC:
    case VELO_NORMAL_BC:
    case VELO_NORMAL_LS_BC:
    case VELO_NORMAL_LS_COLLOC_BC:
    case VELO_NORMAL_LS_PETROV_BC:
    case VELO_TANGENT_LS_BC:
    case VELO_NORM_COLLOC_BC:
    case VELO_TANG1_COLLOC_BC:
    case VELO_TANG2_COLLOC_BC:
    case VELO_NORMAL_DISC_BC:
    case CAP_ENDFORCE_SCALAR_BC:
    case SURFTANG_SCALAR_EDGE_BC:
    case FLOW_PRESSURE_BC:
    case FLOW_PRESSURE_VAR_BC:
    case FLOW_STRESSNOBC_BC:
    case FLOW_GRADV_BC:
    case FLOW_GRADV_T_BC:
    case FLOW_GRADV_SIC_BC:
    case FILL_INLET_BC:
    case FILL_CA_BC:
    case SHARP_CA_2D_BC:
    case STRONG_FILL_CA_BC:
    case WETTING_TENSION_BC:
    case LS_CA_H_BC:
    case LS_T_BC:
    case LS_U_BC:
    case LS_V_BC:
    case LS_W_BC:
    case LS_Q_BC:
    case LS_FLOW_PRESSURE_BC:
    case LS_CAP_HYSING_BC:
    case LS_CAP_DENNER_DIFF_BC:
    case SH_FLUID_STRESS_BC:
    case SH_LUBP_SOLID_BC:
    case SH_LUBP_SOLID_RS_BC:
    case LS_ATTACH_BC:
    case SH_SLOPE_X_BC:
    case SH_SLOPE_Y_BC:
    case APR_VELOCITY_BC:
    case API_VELOCITY_BC:
    case TENSION_SHEET_BC:
    case FRICTION_BC:
    case FRICTION_RS_BC:
    case SHEAR_TO_SHELL_BC:
    case GRAD_LUB_PRESS_BC:
    case LUB_STATIC_BC:
    case SHELL_GRAD_FP_BC:
    case SHELL_GRAD_FH_BC:
    case SHELL_GRAD_PC_BC:
    case LS_WALL_ANGLE_BC:
    case SH_SDET_BC:
    case SH_MESH2_WEAK_BC:
    case RESTIME_GRADSIC_BC:
    case GRAD_LUBP_NOBC_BC:

      if (fscanf(ifp, "%lf", &BC_Types[ibc].BC_Data_Float[0]) != 1) {
        sr = sprintf(err_msg, "%s: Expected 1 flt for %s.", yo, BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      BC_Types[ibc].max_DFlt = 1;

      SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].BC_Data_Float[0]);
      if (BC_Types[ibc].BC_Name == GRAD_LUB_PRESS_BC ||
          BC_Types[ibc].BC_Name == GRAD_LUBP_NOBC_BC) {
        BC_Types[ibc].BC_Data_Float[1] = 0.;
        BC_Types[ibc].BC_Data_Float[2] = 1.;
        if (fscanf(ifp, "%lf %lf", &BC_Types[ibc].BC_Data_Float[1],
                   &BC_Types[ibc].BC_Data_Float[2]) != 2) {
        }
        BC_Types[ibc].max_DFlt = 3;
        SPF(endofstring(echo_string), " %.4g %.4g", BC_Types[ibc].BC_Data_Float[1],
            BC_Types[ibc].BC_Data_Float[2]);
      }
      if (fscanf(ifp, "%d", &BC_Types[ibc].BC_Data_Int[0]) != 1) {
        BC_Types[ibc].BC_Data_Int[0] = -1;
        /* The default for this int now becomes an added sign needed to resolve unhandled issues
         * with determining the sign automagically. Put it in [1], because [0] gets overridden.
         */
        if (BC_Types[ibc].BC_Name == CAP_ENDFORCE_SCALAR_BC ||
            BC_Types[ibc].BC_Name == SURFTANG_SCALAR_BC) {
          BC_Types[ibc].BC_Data_Int[1] = 1;
        }
      } else {
        SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_Data_Int[0]);
        if (BC_Types[ibc].BC_Name == CAP_ENDFORCE_SCALAR_BC ||
            BC_Types[ibc].BC_Name == SURFTANG_SCALAR_BC) {
          if (BC_Types[ibc].BC_Data_Int[0] != 1 && BC_Types[ibc].BC_Data_Int[0] != -1) {
            sr = sprintf(err_msg, "%s: Optional int has to have a value of 1 or -1 (sign): %d", yo,
                         BC_Types[ibc].BC_Data_Int[0]);
            GOMA_EH(GOMA_ERROR, err_msg);
          } else {
            BC_Types[ibc].BC_Data_Int[1] = BC_Types[ibc].BC_Data_Int[0];
          }
        }
      }
      if (BC_Types[ibc].BC_Name == VELO_NORMAL_LS_BC ||
          BC_Types[ibc].BC_Name == VELO_NORMAL_LS_PETROV_BC ||
          BC_Types[ibc].BC_Name == VELO_NORMAL_LS_COLLOC_BC ||
          BC_Types[ibc].BC_Name == VELO_TANGENT_LS_BC) {
        BC_Types[ibc].BC_Data_Float[1] = 0.;
        BC_Types[ibc].BC_Data_Float[2] = 0.;
        BC_Types[ibc].BC_Data_Float[3] = 135.;
        if (fscanf(ifp, "%lf %lf %lf", &BC_Types[ibc].BC_Data_Float[1],
                   &BC_Types[ibc].BC_Data_Float[2], &BC_Types[ibc].BC_Data_Float[3]) != 3) {
          sr = sprintf(err_msg, "%s: Expected 3 flts for %s.", yo, BC_Types[ibc].desc->name1);
          GOMA_EH(GOMA_ERROR, err_msg);
        }
        BC_Types[ibc].max_DFlt = 4;
        SPF(endofstring(echo_string), " %.4g %.4g %.4g", BC_Types[ibc].BC_Data_Float[1],
            BC_Types[ibc].BC_Data_Float[2], BC_Types[ibc].BC_Data_Float[3]);
      }
      break;

      /*
       * DIRICHLET CONDITIONS ONLY!!
       * Fall through for all cases which require a single floating point
       * value as data input and two optional parameters:
       *   BC_relax      -> Time constant for relaxation of bc
       *                    (specify -1.0 for instantaneous)
       *   BC_EBID_Apply -> limit application of this Dirichlet bc to
       *                    the given element block (specify -1 for
       *                    no limitation)
       */
    case DISTNG_BC:
    case DX_BC:
    case DX_RS_BC:
    case DXDISTNG_BC:
    case DY_BC:
    case DY_RS_BC:
    case DYDISTNG_BC:
    case DZ_BC:
    case DZ_RS_BC:
    case DZDISTNG_BC:
    case DX_NOTHING_BC:
    case DY_NOTHING_BC:
    case DZ_NOTHING_BC:
    case P_BC:
    case PSTAR_BC:
    case T_BC:
    case U_BC:
    case V_BC:
    case W_BC:
    case USTAR_BC:
    case VSTAR_BC:
    case WSTAR_BC:
    case PU_BC:
    case PV_BC:
    case PW_BC:
    case N1_BC:
    case N2_BC:
    case N3_BC:
    case LM1_BC:
    case LM2_BC:
    case LM3_BC:
    case S11_BC:
    case S12_BC:
    case S13_BC:
    case S22_BC:
    case S23_BC:
    case S33_BC:
    case S11_1_BC:
    case S12_1_BC:
    case S13_1_BC:
    case S22_1_BC:
    case S23_1_BC:
    case S33_1_BC:
    case S11_2_BC:
    case S12_2_BC:
    case S13_2_BC:
    case S22_2_BC:
    case S23_2_BC:
    case S33_2_BC:
    case S11_3_BC:
    case S12_3_BC:
    case S13_3_BC:
    case S22_3_BC:
    case S23_3_BC:
    case S33_3_BC:
    case S11_4_BC:
    case S12_4_BC:
    case S13_4_BC:
    case S22_4_BC:
    case S23_4_BC:
    case S33_4_BC:
    case S11_5_BC:
    case S12_5_BC:
    case S13_5_BC:
    case S22_5_BC:
    case S23_5_BC:
    case S33_5_BC:
    case S11_6_BC:
    case S12_6_BC:
    case S13_6_BC:
    case S22_6_BC:
    case S23_6_BC:
    case S33_6_BC:
    case S11_7_BC:
    case S12_7_BC:
    case S13_7_BC:
    case S22_7_BC:
    case S23_7_BC:
    case S33_7_BC:

    case G11_BC:
    case G12_BC:
    case G13_BC:
    case G21_BC:
    case G22_BC:
    case G23_BC:
    case G31_BC:
    case G32_BC:
    case G33_BC:
    case VOLT_BC:
    case QS_BC:
    case F_BC:
    case F_DIODE_BC:
    case SH_BC:
    case NN_BC:
    case EXT_V_BC:
    case H_BC:
    case CONT_TANG_VEL_BC:
    case CONT_NORM_VEL_BC:
    case VELO_NORMAL_EDGE_BC:
    case VELO_NORMAL_EDGE_INT_BC:
    case CA_EDGE_CURVE_BC:
    case CA_EDGE_CURVE_INT_BC:
    case POROUS_LIQ_PRESSURE_BC:
    case POROUS_LIQ_FLUX_CONST_BC:
    case POROUS_GAS_PRESSURE_BC:
    case POROUS_GAS_FLUX_CONST_BC:
    case POROUS_TEMP_BC:
    case POROUS_SINK_BC:
    case SH_X_BC:
    case SH_Y_BC:
    case SH_K_BC:
    case SH_TENS_BC:
    case F1_BC:
    case F2_BC:
    case F3_BC:
    case F4_BC:
    case F5_BC:
    case APR_BC:
    case INTP_BC:
    case INTM_BC:
    case INTD_BC:
    case RESTIME_BC:
    case API_BC:
    case SH_USER_BC:
    case SH_LUBP_BC:
    case LUB_PRESS_BC:
    case LUB_PRESS_2_BC:
    case SHELL_TEMP_BC:
    case SHELL_FILMP_BC:
    case SHELL_FILMH_BC:
    case SHELL_PARTC_BC:
    case SHELL_OPEN_PRESS_BC:
    case SHELL_OPEN_PRESS_2_BC:
    case SHELL_SAT_1_BC:
    case SHELL_SAT_2_BC:
    case SHELL_SAT_3_BC:
    case SH_GAMMA1_BC:
    case SHELL_TFMP_PRES_BC:
    case EM_E1R_BC:
    case EM_E1I_BC:
    case EM_E2R_BC:
    case EM_E2I_BC:
    case EM_E3R_BC:
    case EM_E3I_BC:
    case EM_H1R_BC:
    case EM_H1I_BC:
    case EM_H2R_BC:
    case EM_H2I_BC:
    case EM_H3R_BC:
    case EM_H3I_BC:
    case EM_CONT_REAL_BC:
    case EM_CONT_IMAG_BC:
    case SHELL_TFMP_SAT_BC:

      if (fscanf(ifp, "%lf", &BC_Types[ibc].BC_Data_Float[0]) != 1) {
        sprintf(err_msg, "%s: Expected 1 flt for %s.", yo, BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      BC_Types[ibc].max_DFlt = 2;

      SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].BC_Data_Float[0]);

      if (fscanf(ifp, "%lf", &BC_Types[ibc].BC_relax) != 1) {
        BC_Types[ibc].BC_relax = -1.0;
      } else {

        SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].BC_relax);

        if (fscanf(ifp, "%d", &BC_Types[ibc].BC_EBID_Apply) != 1) {
          BC_Types[ibc].BC_EBID_Apply = -1;
        } else
          SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_EBID_Apply);
      }
      break;

      /*
       * Fall through for all cases which require two floating point
       * values as data input
       */
    case QCONV_BC:
    case KIN_LEAK_BC:
    case VNORM_LEAK_BC:
    case VELO_SLIP_EK_BC:
    case LS_EIK_KIN_LEAK_BC:
    case SHEET_ENDSLOPE_BC:
    case LS_EXTV_LATENT_BC:
    case LS_LATENT_HEAT_BC:
    case LS_ACOUSTIC_SOURCE_BC:

      if (fscanf(ifp, "%lf %lf", &BC_Types[ibc].BC_Data_Float[0],
                 &BC_Types[ibc].BC_Data_Float[1]) != 2) {
        sr = sprintf(err_msg, "%s: Expected 2 flts for %s.", yo, BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      BC_Types[ibc].max_DFlt = 2;
      SPF(endofstring(echo_string), " %.4g %.4g", BC_Types[ibc].BC_Data_Float[0],
          BC_Types[ibc].BC_Data_Float[1]);
      break;

    case LS_STRESS_JUMP_BC:
      if (fscanf(ifp, "%lf %lf", &BC_Types[ibc].BC_Data_Float[0],
                 &BC_Types[ibc].BC_Data_Float[1]) != 2) {
        sr = sprintf(err_msg, "%s: Expected 2 flts for %s.", yo, BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      BC_Types[ibc].max_DFlt = 2;
      /* Try reading optional int */
      if (fscanf(ifp, "%d", &BC_Types[ibc].BC_Data_Int[0]) != 1) {
        BC_Types[ibc].BC_Data_Int[0] = 0;
      }
      if ((BC_Types[ibc].BC_Data_Int[0] != 0) && (BC_Types[ibc].BC_Data_Int[0] != 1) &&
          (BC_Types[ibc].BC_Data_Int[0] != -1)) {
        GOMA_EH(GOMA_ERROR,
                "Expected LS_STRESS_JUMP float float [-1,0,1] got LS_STRESS_JUMP %g %g %d",
                BC_Types[ibc].BC_Data_Float[0], BC_Types[ibc].BC_Data_Float[1],
                BC_Types[ibc].BC_Data_Int[0]);
      }
      SPF(endofstring(echo_string), " %.4g %.4g %d", BC_Types[ibc].BC_Data_Float[0],
          BC_Types[ibc].BC_Data_Float[1], BC_Types[ibc].BC_Data_Int[0]);
      break;

    case CAPILLARY_SHEAR_VISC_BC:
      if (fscanf(ifp, "%lf %lf", &BC_Types[ibc].BC_Data_Float[0],
                 &BC_Types[ibc].BC_Data_Float[1]) != 2) {
        sr = sprintf(err_msg, "%s: Expected 2 flts for %s.", yo, BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      /* Try reading two optional floats [2] is t_start_on , while [3] is t_fully_on */
      if (fscanf(ifp, "%lf %lf", &BC_Types[ibc].BC_Data_Float[2],
                 &BC_Types[ibc].BC_Data_Float[3]) != 2) {
      }
      BC_Types[ibc].max_DFlt = 4;
      SPF(endofstring(echo_string), " %.4g %.4g %.4g %.4g", BC_Types[ibc].BC_Data_Float[0],
          BC_Types[ibc].BC_Data_Float[1], BC_Types[ibc].BC_Data_Float[2],
          BC_Types[ibc].BC_Data_Float[3]);

      break;

      /* Fall through for all cases which requires two integer values
       * as data input
       */
    case YFLUX_ETCH_BC:

      if (fscanf(ifp, "%d %d", &BC_Types[ibc].BC_Data_Int[0], &BC_Types[ibc].BC_Data_Int[1]) != 2) {
        sr = sprintf(err_msg, "%s: Expected 2 int for %s.", yo, BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      SPF(endofstring(echo_string), " %d %d", BC_Types[ibc].BC_Data_Int[0],
          BC_Types[ibc].BC_Data_Int[1]);
      BC_Types[ibc].species_eq = BC_Types[ibc].BC_Data_Int[0];
      break;

      /*
       * Fall through for all cases which require three floating point
       * values as data input
       */
    case LS_ADC_OLD_BC:
    case LS_ADC_BC:
      srand((long)ut()); /* Seed the random number generator  when LS_ADC is used */
                         /* fall through */
    case FORCE_BC:
    case FORCE_SIC_BC:
    case FORCE_RS_BC:
    case NORM_FORCE_BC:
    case NORM_FORCE_RS_BC:
    case QSIDE_DIR_BC:
    case SLOPEX_BC:
    case SLOPEY_BC:
    case SLOPEZ_BC:
    case SLOPE_BC:
    case VELO_TANGENT_EDGE_BC:
    case VELO_TANGENT_EDGE_INT_BC:
    case SH_S11_WEAK_BC:
    case SH_S22_WEAK_BC:
    case VELO_NORMAL_LUB_BC:

      if (fscanf(ifp, "%lf %lf %lf", &BC_Types[ibc].BC_Data_Float[0],
                 &BC_Types[ibc].BC_Data_Float[1], &BC_Types[ibc].BC_Data_Float[2]) != 3) {
        sr = sprintf(err_msg, "%s: Expected 3 flts for %s.", yo, BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      SPF_DBL_VEC(endofstring(echo_string), 3, BC_Types[ibc].BC_Data_Float);
      BC_Types[ibc].max_DFlt = 3;

      break;

      /*
       * Fall through for all cases which require three floating point
       * values as data input and one optional integer.
       */
    case CAPILLARY_BC:
    case CAPILLARY_TABLE_BC:
    case PF_CAPILLARY_BC:
    case LIGHTP_TRANS_BC:
    case LIGHTM_TRANS_BC:
    case LIGHTD_TRANS_BC:

      if (fscanf(ifp, "%lf %lf %lf", &BC_Types[ibc].BC_Data_Float[0],
                 &BC_Types[ibc].BC_Data_Float[1], &BC_Types[ibc].BC_Data_Float[2]) != 3) {
        sr = sprintf(err_msg, "%s: Expected 3 flts for %s.", yo, BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      SPF_DBL_VEC(endofstring(echo_string), 3, BC_Types[ibc].BC_Data_Float);
      BC_Types[ibc].max_DFlt = 3;

      /* Try reading the optional integer. */
      if (fscanf(ifp, "%d", &BC_Types[ibc].BC_Data_Int[0]) != 1) {
        BC_Types[ibc].BC_Data_Int[0] = -1;
      } else {
        SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_Data_Int[0]);
      }

      if (BC_Types[ibc].BC_Name == CAPILLARY_TABLE_BC) {
        new_BC_Desc = ((struct BC_descriptions **)realloc(
            new_BC_Desc, (num_new_BC_Desc + 1) * sizeof(struct BC_descriptions *)));

        new_BC_Desc[num_new_BC_Desc] = alloc_BC_description(BC_Types[ibc].desc);

        BC_Types[ibc].desc = new_BC_Desc[num_new_BC_Desc];

        BC_Types[ibc].index_dad = num_new_BC_Desc++; /* This is important to Phil */

        if (num_BC_Tables == MAX_BC_TABLES) {
          GOMA_EH(GOMA_ERROR, "Maximum TABLE_BCs exceeded .");
        }

        BC_Tables[num_BC_Tables] = setup_table_BC(ifp, input, &BC_Types[ibc], echo_string);

        BC_Types[ibc].table_index = num_BC_Tables++;
      }
      break;

      /*
       * Fall through for all cases which require four floating point
       * values as data input
       */
    case SURFTANG_BC:
      GOMA_WH(-1, "Use CAP_ENDFORCE for consistent sign convention");
      /*FALLTHROUGH*/
    case PLANEX_BC:
    case PLANEY_BC:
    case PLANEZ_BC:
    case PLANE_BC:
    case CA_EDGE_BC:
    case CA_EDGE_INT_BC:
    case CAP_ENDFORCE_BC:
    case QRAD_BC:
    case LS_QRAD_BC:
    case SURFTANG_EDGE_BC:
    case FLOW_HYDROSTATIC_BC:
    case LUB_PRESS_HYDROSTATIC_BC:
    case VELO_TANGENT_3D_BC:
    case SHARP_WETLIN_VELOCITY_BC:
    case APR_PLANE_TRANS_BC:
    case API_PLANE_TRANS_BC:
    case DX_USER_NODE_BC:
    case DY_USER_NODE_BC:
    case DZ_USER_NODE_BC:
      if (fscanf(ifp, "%lf %lf %lf %lf", &BC_Types[ibc].BC_Data_Float[0],
                 &BC_Types[ibc].BC_Data_Float[1], &BC_Types[ibc].BC_Data_Float[2],
                 &BC_Types[ibc].BC_Data_Float[3]) != 4) {
        sr = sprintf(err_msg, "%s: Expected 4 flts for %s.", yo, BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      BC_Types[ibc].max_DFlt = 4;
      SPF_DBL_VEC(endofstring(echo_string), 4, BC_Types[ibc].BC_Data_Float);
      break;
      /*
       * Fall through for all cases which require four floating point
       * values as data input and an optional integer element block
       */
    case CA_BC:
    case CA_MOMENTUM_BC:
      if (fscanf(ifp, "%lf %lf %lf %lf", &BC_Types[ibc].BC_Data_Float[0],
                 &BC_Types[ibc].BC_Data_Float[1], &BC_Types[ibc].BC_Data_Float[2],
                 &BC_Types[ibc].BC_Data_Float[3]) != 4) {
        sr = sprintf(err_msg, "%s: Expected 4 flts for %s.", yo, BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      } else {
        SPF_DBL_VEC(endofstring(echo_string), 4, BC_Types[ibc].BC_Data_Float);
        /* Scan for the optional int. If not present, put a -1 in second data position */
        /* This is to ensure a nonzero entry in BC_Data_Int[2] for CA */
        /* note: optional int isn't listed in the manual. What's it for? */
        if (fscanf(ifp, "%d", &BC_Types[ibc].BC_Data_Int[2]) != 1) {
          BC_Types[ibc].BC_Data_Int[2] = -1;
        } else
          SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_Data_Int[2]);
      }
      BC_Types[ibc].max_DFlt = 4;

      break;

      /*
       * Fall through for all cases which require four floating point
       * values as data input and result in the creation of solid model
       * (CGM) entities
       */
    case SM_PLANE_BC: /* Solid Model PLANE BC */
      GOMA_EH(GOMA_ERROR, "CGM not supported, SM_PLANE_BC");
      break;

    /*
     *   modify VELO_SLIP conditions for position dependent slip
     */
    case VELO_SLIP_BC:
    case VELO_SLIP_ROT_BC:
    case VELO_SLIP_FLUID_BC:
    case VELO_SLIP_ROT_FLUID_BC:
    case VELO_SLIP_FILL_BC:
    case VELO_SLIP_ROT_FILL_BC:
    case AIR_FILM_BC:
    case AIR_FILM_ROT_BC:
      if (fscanf(ifp, "%lf %lf %lf %lf", &BC_Types[ibc].BC_Data_Float[0],
                 &BC_Types[ibc].BC_Data_Float[1], &BC_Types[ibc].BC_Data_Float[2],
                 &BC_Types[ibc].BC_Data_Float[3]) != 4) {
        sr = sprintf(err_msg, "%s: Expected 4 flts for %s.", yo, BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      BC_Types[ibc].max_DFlt = 4;
      SPF_DBL_VEC(endofstring(echo_string), 4, BC_Types[ibc].BC_Data_Float);

      if (fscanf(ifp, "%d %lf", &BC_Types[ibc].BC_Data_Int[0], &BC_Types[ibc].BC_Data_Float[4]) !=
          2) {
        BC_Types[ibc].BC_Data_Int[0] = -1;
        BC_Types[ibc].BC_Data_Float[4] = 0.;
      } else
        SPF(endofstring(echo_string), " %d %.4g", BC_Types[ibc].BC_Data_Int[0],
            BC_Types[ibc].BC_Data_Float[4]);
      if (fscanf(ifp, "%lf", &BC_Types[ibc].BC_Data_Float[5]) != 1) {
        BC_Types[ibc].BC_Data_Float[5] = 0.;
      } else
        SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].BC_Data_Float[5]);
      BC_Types[ibc].max_DFlt = 5;

      if (fscanf(ifp, "%lf", &BC_Types[ibc].BC_Data_Float[6]) != 1) {
        BC_Types[ibc].BC_Data_Float[6] = 0.;
      } else
        SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].BC_Data_Float[6]);
      BC_Types[ibc].max_DFlt = 6;

      if (fscanf(ifp, "%lf", &BC_Types[ibc].BC_Data_Float[7]) != 1) {
        BC_Types[ibc].BC_Data_Float[7] = 1;
      } else
        SPF(endofstring(echo_string), " %g", BC_Types[ibc].BC_Data_Float[7]);
      BC_Types[ibc].max_DFlt = 7;

      break;
    case VELO_SLIP_POWER_CARD_BC:
      if (fscanf(ifp, "%lf %lf %lf %lf %lf", &BC_Types[ibc].BC_Data_Float[0],
                 &BC_Types[ibc].BC_Data_Float[1], &BC_Types[ibc].BC_Data_Float[2],
                 &BC_Types[ibc].BC_Data_Float[3], &BC_Types[ibc].BC_Data_Float[4]) != 5) {
        sr = sprintf(err_msg, "%s: Expected 5 flts for %s.", yo, BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      BC_Types[ibc].max_DFlt = 5;
      SPF_DBL_VEC(endofstring(echo_string), 5, BC_Types[ibc].BC_Data_Float);
      break;
    case VELO_SLIP_POWER_BC:
      if (fscanf(ifp, "%lf %lf %lf %lf %lf", &BC_Types[ibc].BC_Data_Float[0],
                 &BC_Types[ibc].BC_Data_Float[1], &BC_Types[ibc].BC_Data_Float[2],
                 &BC_Types[ibc].BC_Data_Float[3], &BC_Types[ibc].BC_Data_Float[4]) != 5) {
        sr = sprintf(err_msg, "%s: Expected 5 flts for %s.", yo, BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      BC_Types[ibc].max_DFlt = 5;
      SPF_DBL_VEC(endofstring(echo_string), 5, BC_Types[ibc].BC_Data_Float);

      // Read in 3 tangential vector components for 3D

      if (fscanf(ifp, "%lf", &BC_Types[ibc].BC_Data_Float[5]) != 1) {
        BC_Types[ibc].BC_Data_Float[5] = 0.0;
      } else {
        SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].BC_Data_Float[5]);
        BC_Types[ibc].max_DFlt = 6;
      }

      if (fscanf(ifp, "%lf", &BC_Types[ibc].BC_Data_Float[6]) != 1) {
        BC_Types[ibc].BC_Data_Float[6] = 0.0;
      } else {
        SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].BC_Data_Float[6]);
        BC_Types[ibc].max_DFlt = 7;
      }

      if (fscanf(ifp, "%lf", &BC_Types[ibc].BC_Data_Float[7]) != 1) {
        BC_Types[ibc].BC_Data_Float[7] = 0.0;
      } else {
        SPF(endofstring(echo_string), " %g", BC_Types[ibc].BC_Data_Float[7]);
        BC_Types[ibc].max_DFlt = 8;
      }
      break;

      /*
       * Fall through for all cases which require five floating point
       * values as data input
       */

    case VELO_EK_3D_BC:
    case VELO_SLIP_LEVEL_BC:
    case VELO_SLIP_LEVEL_SIC_BC:
    case VELO_SLIP_LS_ROT_BC:
    case SHARP_HOFFMAN_VELOCITY_BC:
      if (fscanf(ifp, "%lf %lf %lf %lf %lf", &BC_Types[ibc].BC_Data_Float[0],
                 &BC_Types[ibc].BC_Data_Float[1], &BC_Types[ibc].BC_Data_Float[2],
                 &BC_Types[ibc].BC_Data_Float[3], &BC_Types[ibc].BC_Data_Float[4]) != 5) {
        sr = sprintf(err_msg, "%s: Expected 5 flts for %s.", yo, BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      /*   outside slip coeff, gas_factor,  contact_friction   */
      if (fscanf(ifp, "%lf", &BC_Types[ibc].BC_Data_Float[5]) != 1) {
        BC_Types[ibc].BC_Data_Float[5] = 1.0 / LITTLE_PENALTY;
      }
      if (fscanf(ifp, "%lf", &BC_Types[ibc].BC_Data_Float[6]) != 1) {
        BC_Types[ibc].BC_Data_Float[6] = 8.0;
      }
      if (fscanf(ifp, "%lf", &BC_Types[ibc].BC_Data_Float[7]) != 1) {
        BC_Types[ibc].BC_Data_Float[7] = 1e-6;
      }
      if (fscanf(ifp, "%lf", &BC_Types[ibc].BC_Data_Float[8]) != 1) {
        BC_Types[ibc].BC_Data_Float[8] = 0.0;
      }

      BC_Types[ibc].max_DFlt = 9;
      SPF_DBL_VEC(endofstring(echo_string), 9, BC_Types[ibc].BC_Data_Float);

      break;

    case VELO_SLIP_LS_ORIENTED_BC:
      if (fscanf(ifp, "%lf %lf %lf %lf %lf",
                 &BC_Types[ibc].BC_Data_Float[0],       // ls_width
                 &BC_Types[ibc].BC_Data_Float[1],       // beta_negative
                 &BC_Types[ibc].BC_Data_Float[2],       // beta_positive
                 &BC_Types[ibc].BC_Data_Float[3],       // gamma_negative
                 &BC_Types[ibc].BC_Data_Float[4]) != 5) // gamma_positive
      {
        GOMA_EH(GOMA_ERROR, "%s: Expected 5 flts for %s.", yo, BC_Types[ibc].desc->name1);
      }

      if (fscanf(ifp, "%lf %lf %lf",
                 &BC_Types[ibc].BC_Data_Float[5],       // v_x
                 &BC_Types[ibc].BC_Data_Float[6],       // v_y
                 &BC_Types[ibc].BC_Data_Float[7]) != 5) // v_z
      {
        BC_Types[ibc].BC_Data_Float[5] = 0.0;
        BC_Types[ibc].BC_Data_Float[6] = 0.0;
        BC_Types[ibc].BC_Data_Float[7] = 0.0;
      }

      SPF_DBL_VEC(endofstring(echo_string), 8, BC_Types[ibc].BC_Data_Float);
      break;

    case VELO_SLIP_LS_HEAVISIDE_BC:
      if (fscanf(ifp, "%lf %lf %lf %lf %lf %lf",
                 &BC_Types[ibc].BC_Data_Float[0],       // ls_width
                 &BC_Types[ibc].BC_Data_Float[1],       // beta_negative
                 &BC_Types[ibc].BC_Data_Float[2],       // beta_positive
                 &BC_Types[ibc].BC_Data_Float[3],       // v_x
                 &BC_Types[ibc].BC_Data_Float[4],       // v_y
                 &BC_Types[ibc].BC_Data_Float[5]) != 6) // v_z
      {
        sr = sprintf(err_msg, "%s: Expected 6 flts for %s.", yo, BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      SPF_DBL_VEC(endofstring(echo_string), 6, BC_Types[ibc].BC_Data_Float);
      break;
      /*
       * Fall through for all cases which require five floating point
       * values as data input plus optional parameters
       */

    case REP_FORCE_BC:
    case REP_FORCE_RS_BC:
    case ATTR_FORCE_BC:
    case ATTR_FORCE_RS_BC:
      if (fscanf(ifp, "%lf %lf %lf %lf %lf", &BC_Types[ibc].BC_Data_Float[0],
                 &BC_Types[ibc].BC_Data_Float[1], &BC_Types[ibc].BC_Data_Float[2],
                 &BC_Types[ibc].BC_Data_Float[3], &BC_Types[ibc].BC_Data_Float[4]) != 5) {
        sr = sprintf(err_msg, "%s: Expected 5 flts for %s.", yo, BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      /*  optional repulsive force exponent   */
      if (fscanf(ifp, "%lf", &BC_Types[ibc].BC_Data_Float[5]) != 1) {
        BC_Types[ibc].BC_Data_Float[5] = 4.0;
      }
      /*  optional friction coefficient (static)   */
      if (fscanf(ifp, "%lf", &BC_Types[ibc].BC_Data_Float[6]) != 1) {
        BC_Types[ibc].BC_Data_Float[6] = 0.0;
      }

      BC_Types[ibc].max_DFlt = 7;
      SPF_DBL_VEC(endofstring(echo_string), 7, BC_Types[ibc].BC_Data_Float);

      break;
      /*
       * Fall through for all cases which require 3 integers and
       * five floating point values as data input
       */
    case VL_EQUIL_BC:
      if (fscanf(ifp, "%d %d %d %lf %lf %lf %lf %lf", &BC_Types[ibc].BC_Data_Int[0],
                 &BC_Types[ibc].BC_Data_Int[1], &BC_Types[ibc].BC_Data_Int[2],
                 &BC_Types[ibc].BC_Data_Float[0], &BC_Types[ibc].BC_Data_Float[1],
                 &BC_Types[ibc].BC_Data_Float[2], &BC_Types[ibc].BC_Data_Float[3],
                 &BC_Types[ibc].BC_Data_Float[4]) != 8) {
        sprintf(err_msg, "%s: Expected 3 ints and 5 flts for %s.", yo, BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      BC_Types[ibc].species_eq = BC_Types[ibc].BC_Data_Int[0];
      assign_global_species_var_type(SPECIES_MASS_FRACTION, FALSE);

      SPF_INT_VEC(endofstring(echo_string), 3, BC_Types[ibc].BC_Data_Int);
      SPF_DBL_VEC(endofstring(echo_string), 5, BC_Types[ibc].BC_Data_Float);
      BC_Types[ibc].max_DFlt = 5;

      break;

      /*
       * Fall through for all cases which require 1 keyword, 1 integer
       *  and 3 floating point values as data input
       */

    case YFLUX_EQUIL_BC:

      if (fscanf(ifp, "%80s", input) != 1) {
        sprintf(err_msg, "%s: Expected a keyword RAOULT or FLORY for %s.", yo,
                BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      } else {
        if (!strcmp(input, "RAOULT")) {
          BC_Types[ibc].BC_Data_Int[2] = RAOULT;
        } else if (!strcmp(input, "FLORY")) {
          BC_Types[ibc].BC_Data_Int[2] = FLORY;
        } else if (!strcmp(input, "FLORY_CC")) {
          BC_Types[ibc].BC_Data_Int[2] = FLORY_CC;
        } else {
          GOMA_EH(GOMA_ERROR, "I don't recognize your YFLUX_EQUIL Keyword!");
        }
      }
      if (fscanf(ifp, "%d %lf %lf %lf", &BC_Types[ibc].BC_Data_Int[0],
                 &BC_Types[ibc].BC_Data_Float[0], &BC_Types[ibc].BC_Data_Float[1],
                 &BC_Types[ibc].BC_Data_Float[2]) != 4) {
        sr = sprintf(err_msg, "%s: Expected 1 int and 3 flt following %s for  %s.", yo, input,
                     BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      SPF(endofstring(echo_string), " %s %d", input, BC_Types[ibc].BC_Data_Int[0]);
      SPF_DBL_VEC(endofstring(echo_string), 3, BC_Types[ibc].BC_Data_Float);

      /**	extra floats for Flory/Chilton-Coburn correlation	**/

      if (BC_Types[ibc].BC_Data_Int[2] == FLORY_CC) {
        if (fscanf(ifp, "%lf %lf", &BC_Types[ibc].BC_Data_Float[3],
                   &BC_Types[ibc].BC_Data_Float[4]) != 2) {
          sr = sprintf(err_msg, "%s: Expected 2 additional flts for  %s.", yo,
                       BC_Types[ibc].desc->name1);
          GOMA_EH(GOMA_ERROR, err_msg);
        }

        SPF_DBL_VEC(endofstring(echo_string), 2, &(BC_Types[ibc].BC_Data_Float[3]));

      } else {
        BC_Types[ibc].BC_Data_Float[3] = 0.;
        BC_Types[ibc].BC_Data_Float[4] = 0.;
      }

      BC_Types[ibc].max_DFlt = 5;
      BC_Types[ibc].species_eq = BC_Types[ibc].BC_Data_Int[0];

      break;

      /*
       * Fall through for all cases which require 1 keyword, 1 integer
       *  and 10 floating point values as data input
       */

    case YFLUX_SULFIDATION_BC:

      if (fscanf(ifp, "%80s", input) != 1) {
        sprintf(err_msg, "%s: Expected a keyword SOLID_DIFFUSION or GAS_DIFFUSION for %s.", yo,
                BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      } else {
        if (!strcmp(input, "SOLID_DIFFUSION_SIMPLIFIED")) {
          BC_Types[ibc].BC_Data_Int[2] = SOLID_DIFFUSION_SIMPLIFIED;
        } else if (!strcmp(input, "SOLID_DIFFUSION_ELECTRONEUTRALITY")) {
          BC_Types[ibc].BC_Data_Int[2] = SOLID_DIFFUSION_ELECTRONEUTRALITY;
        } else if (!strcmp(input, "SOLID_DIFFUSION_ELECTRONEUTRALITY_LINEAR")) {
          BC_Types[ibc].BC_Data_Int[2] = SOLID_DIFFUSION_ELECTRONEUTRALITY_LINEAR;
        } else if (!strcmp(input, "SOLID_DIFFUSION")) {
          BC_Types[ibc].BC_Data_Int[2] = SOLID_DIFFUSION;
        } else if (!strcmp(input, "GAS_DIFFUSION")) {
          BC_Types[ibc].BC_Data_Int[2] = GAS_DIFFUSION;
        } else if (!strcmp(input, "FULL")) {
          BC_Types[ibc].BC_Data_Int[2] = FULL;
        } else if (!strcmp(input, "ANNIHILATION_ELECTRONEUTRALITY")) {
          BC_Types[ibc].BC_Data_Int[2] = ANNIHILATION_ELECTRONEUTRALITY;
        } else if (!strcmp(input, "ANNIHILATION")) {
          BC_Types[ibc].BC_Data_Int[2] = ANNIHILATION;
        } else {
          GOMA_EH(GOMA_ERROR, "I don't recognize your YFLUX_SULFIDATION Keyword!");
        }
      }
      if (fscanf(ifp, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &BC_Types[ibc].BC_Data_Int[0],
                 &BC_Types[ibc].BC_Data_Float[0], &BC_Types[ibc].BC_Data_Float[1],
                 &BC_Types[ibc].BC_Data_Float[2], &BC_Types[ibc].BC_Data_Float[3],
                 &BC_Types[ibc].BC_Data_Float[4], &BC_Types[ibc].BC_Data_Float[5],
                 &BC_Types[ibc].BC_Data_Float[6], &BC_Types[ibc].BC_Data_Float[7],
                 &BC_Types[ibc].BC_Data_Float[8], &BC_Types[ibc].BC_Data_Float[9]) != 11) {
        sr = sprintf(err_msg, "%s: Expected 1 int and 10 flt following %s for  %s.", yo, input,
                     BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      BC_Types[ibc].species_eq = BC_Types[ibc].BC_Data_Int[0];

      SPF(endofstring(echo_string), " %s %d", input, BC_Types[ibc].BC_Data_Int[0]);
      SPF_DBL_VEC(endofstring(echo_string), 10, BC_Types[ibc].BC_Data_Float);
      BC_Types[ibc].max_DFlt = 10;

      break;

      /*
       * Fall through for all cases which require 1 keyword, 4 integers
       *  and 1 floating point values as data input
       */

    case VL_POLY_BC:

      if (fscanf(ifp, "%80s", input) != 1) {
        sr = sprintf(err_msg, "%s: Expected a keyword MASS or VOLUME for %s.", yo,
                     BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      } else {
        if (!strcmp(input, "VOLUME")) {
          BC_Types[ibc].BC_Data_Int[0] = VOLUME;
        } else if (!strcmp(input, "MASS")) {
          BC_Types[ibc].BC_Data_Int[0] = MASS;
        } else {
          GOMA_EH(GOMA_ERROR, "I don't recognize your VL_POLY Keyword!");
        }
      }

      if (fscanf(ifp, "%d %d %d %lf", &BC_Types[ibc].BC_Data_Int[1], &BC_Types[ibc].BC_Data_Int[2],
                 &BC_Types[ibc].BC_Data_Int[3], &BC_Types[ibc].BC_Data_Float[0]) != 4) {
        sr = sprintf(err_msg, "%s: Expected 3 ints and 1 flt for %s.", yo,
                     BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      BC_Types[ibc].species_eq = BC_Types[ibc].BC_Data_Int[1];

      SPF(endofstring(echo_string), " %s", input);
      for (i = 0; i < 3; i++)
        SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_Data_Int[i]);
      SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].BC_Data_Float[0]);

      break;

      /*
       * Fall through for all cases which require six floating point
       * values as data input
       */

    case SHARP_BLAKE_VELOCITY_BC:
    case SHARP_COX_VELOCITY_BC:
    case SHARP_SHIK_VELOCITY_BC:
      if (fscanf(ifp, "%lf %lf %lf %lf %lf %lf ", &BC_Types[ibc].BC_Data_Float[0],
                 &BC_Types[ibc].BC_Data_Float[1], &BC_Types[ibc].BC_Data_Float[2],
                 &BC_Types[ibc].BC_Data_Float[3], &BC_Types[ibc].BC_Data_Float[4],
                 &BC_Types[ibc].BC_Data_Float[5]) != 6) {
        sr = sprintf(err_msg, "%s: Expected 6 flts for %s.", yo, BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      for (i = 0; i < 6; i++)
        SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].BC_Data_Float[i]);
      break;

      /*
       * Fall through for all cases which require seven floating point
       * values as data input
       */
    case CA_OR_FIX_BC:
    case MOVING_PLANE_BC:
    case WETTING_SPEED_LIN_BC:

      if (fscanf(ifp, "%lf %lf %lf %lf %lf %lf %lf", &BC_Types[ibc].BC_Data_Float[0],
                 &BC_Types[ibc].BC_Data_Float[1], &BC_Types[ibc].BC_Data_Float[2],
                 &BC_Types[ibc].BC_Data_Float[3], &BC_Types[ibc].BC_Data_Float[4],
                 &BC_Types[ibc].BC_Data_Float[5], &BC_Types[ibc].BC_Data_Float[6]) != 7) {
        sr = sprintf(err_msg, "%s: Expected 7 flts for %s.", yo, BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      for (i = 0; i < 7; i++)
        SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].BC_Data_Float[i]);

      /*initialize pinning flag to zero, which is BC_Data_Int[0]
        for each CA_OR_FIX Card */
      if (!strcmp(BC_Types[ibc].desc->name1, "CA_OR_FIX")) {
        BC_Types[ibc].BC_Data_Int[0] = 0;
      }

      break;

      /*
       * Fall through for all cases which require eight floating point (used to be 7)
       * values as data input and one optional integer (BC_Data_int[0]).
       */
    case CAP_REPULSE_BC:

      if (fscanf(ifp, "%lf %lf %lf %lf %lf %lf %lf %lf", &BC_Types[ibc].BC_Data_Float[0],
                 &BC_Types[ibc].BC_Data_Float[1], &BC_Types[ibc].BC_Data_Float[2],
                 &BC_Types[ibc].BC_Data_Float[3], &BC_Types[ibc].BC_Data_Float[4],
                 &BC_Types[ibc].BC_Data_Float[5], &BC_Types[ibc].BC_Data_Float[6],
                 &BC_Types[ibc].BC_Data_Float[7]) != 8) {
        sr = sprintf(err_msg, "%s: Expected 8 flts for %s.", yo, BC_Types[ibc].desc->name1);
        GOMA_EH(-1, err_msg);
      }

      for (i = 0; i < 8; i++)
        SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].BC_Data_Float[i]);

      /* Try reading the optional integer. */
      if (fscanf(ifp, "%d", &BC_Types[ibc].BC_Data_Int[0]) != 1) {
        BC_Types[ibc].BC_Data_Int[0] = -1;
      } else
        SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_Data_Int[0]);

      break;

      /*
       * Fall through for all cases which require seven floating point
       * values as data input and one optional integer.
       */
    case LS_RECOIL_PRESSURE_BC:
    case CAP_RECOIL_PRESS_BC:

      if (fscanf(ifp, "%lf %lf %lf %lf %lf %lf %lf", &BC_Types[ibc].BC_Data_Float[0],
                 &BC_Types[ibc].BC_Data_Float[1], &BC_Types[ibc].BC_Data_Float[2],
                 &BC_Types[ibc].BC_Data_Float[3], &BC_Types[ibc].BC_Data_Float[4],
                 &BC_Types[ibc].BC_Data_Float[5], &BC_Types[ibc].BC_Data_Float[6]) != 7) {
        sr = sprintf(err_msg, "%s: Expected 7 flts for %s.", yo, BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      for (i = 0; i < 7; i++)
        SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].BC_Data_Float[i]);

      /* Try reading the optional integer. */
      if (fscanf(ifp, "%d", &BC_Types[ibc].BC_Data_Int[0]) != 1) {
        BC_Types[ibc].BC_Data_Int[0] = -1;
      } else
        SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_Data_Int[0]);

      break;

      /*
       * Fall through for all cases which require eight floating point
       * values as data input
       */

    case VAR_CA_EDGE_BC:
    case VELO_THETA_TPL_BC:
    case VELO_THETA_HOFFMAN_BC:
    case VELO_THETA_COX_BC:
    case WETTING_SPEED_BLAKE_BC:
    case WETTING_SPEED_HOFFMAN_BC:
    case WETTING_SPEED_COX_BC:
    case WETTING_SPEED_SHIK_BC:
    case LINEAR_WETTING_SIC_BC:
    case BLAKE_DIRICHLET_BC:
    case HOFFMAN_DIRICHLET_BC:
    case COX_DIRICHLET_BC:
    case BLAKE_DIRICH_ROLL_BC:
    case HOFFMAN_DIRICH_ROLL_BC:
    case COX_DIRICH_ROLL_BC:
    case EM_ER_FARFIELD_DIRECT_BC:
    case EM_EI_FARFIELD_DIRECT_BC:
    case EM_HR_FARFIELD_DIRECT_BC:
    case EM_HI_FARFIELD_DIRECT_BC:
    case E_ER_FARFIELD_BC:
    case E_EI_FARFIELD_BC:
    case EM_FARFIELD_REAL_NED_BC:
    case EM_FARFIELD_IMAG_NED_BC:

      if (fscanf(ifp, "%lf %lf %lf %lf %lf %lf %lf %lf", &BC_Types[ibc].BC_Data_Float[0],
                 &BC_Types[ibc].BC_Data_Float[1], &BC_Types[ibc].BC_Data_Float[2],
                 &BC_Types[ibc].BC_Data_Float[3], &BC_Types[ibc].BC_Data_Float[4],
                 &BC_Types[ibc].BC_Data_Float[5], &BC_Types[ibc].BC_Data_Float[6],
                 &BC_Types[ibc].BC_Data_Float[7]) != 8) {
        sprintf(err_msg, "%s: Expected 8 flts for BC %s.\n", yo, BC_Types[ibc].desc->name1);
        strcat(err_msg,
               "Be advised syntax to this boundary condition has changed. Consult GOMA manual.\n");
        GOMA_EH(GOMA_ERROR, err_msg);
      } else { /* 1*/
        for (i = 0; i < 8; i++)
          SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].BC_Data_Float[i]);

        if (fscanf(ifp, "%d", &BC_Types[ibc].BC_Data_Int[2]) != 1) {
          BC_Types[ibc].BC_Data_Int[2] = -1;
        } else
          SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_Data_Int[2]);

        if (fscanf(ifp, "%lf", &BC_Types[ibc].BC_Data_Float[8]) != 1) {
          BC_Types[ibc].BC_Data_Float[8] = 1.0;
        } else
          SPF(endofstring(echo_string), " %lf", BC_Types[ibc].BC_Data_Float[8]);
        /*   outside slip coeff, gas_factor,  contact_friction   */
        if (fscanf(ifp, "%lf", &BC_Types[ibc].BC_Data_Float[9]) != 1) {
          if (!strcmp(BC_Types[ibc].desc->name1, "VELO_THETA_HOFFMAN")) {
            /* Max DCA for Hoffman condition  */
            BC_Types[ibc].BC_Data_Float[9] = 180.0;
          } else {
            BC_Types[ibc].BC_Data_Float[9] = 1.0 / sqrt(LITTLE_PENALTY * BIG_PENALTY);
          }
        } else
          SPF(endofstring(echo_string), " %lf", BC_Types[ibc].BC_Data_Float[9]);
        if (fscanf(ifp, "%lf", &BC_Types[ibc].BC_Data_Float[10]) != 1) {
          if (!strcmp(BC_Types[ibc].desc->name1, "VELO_THETA_HOFFMAN") ||
              !strcmp(BC_Types[ibc].desc->name1, "VELO_THETA_COX")) {
            /* DCL shearrate for Hoffman/Cox condition  */
            BC_Types[ibc].BC_Data_Float[10] = -1.0;
          } else {
            BC_Types[ibc].BC_Data_Float[10] = 8.0;
          }
        } else
          SPF(endofstring(echo_string), " %lf", BC_Types[ibc].BC_Data_Float[10]);
        if (fscanf(ifp, "%lf", &BC_Types[ibc].BC_Data_Float[11]) != 1) {
          BC_Types[ibc].BC_Data_Float[11] = 1e-6;
        } else
          SPF(endofstring(echo_string), " %lf", BC_Types[ibc].BC_Data_Float[11]);
      }

      BC_Types[ibc].BC_Data_Int[3] = (BC_Types[ibc].BC_Name == BLAKE_DIRICH_ROLL_BC ||
                                      BC_Types[ibc].BC_Name == HOFFMAN_DIRICH_ROLL_BC ||
                                      BC_Types[ibc].BC_Name == COX_DIRICH_ROLL_BC)
                                         ? TRUE
                                         : FALSE;

      break;
      /*
       * Fall through for all cases which require nine floating point
       * values as data input
       */

    case VELO_THETA_SHIK_BC:
    case SHIK_DIRICHLET_BC:
    case SHIK_DIRICH_ROLL_BC:
    case HYSTERESIS_WETTING_BC:
    case EM_ER_SOMMERFELD_BC:
    case EM_EI_SOMMERFELD_BC:
    case EM_HR_SOMMERFELD_BC:
    case EM_HI_SOMMERFELD_BC:
      if (fscanf(ifp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &BC_Types[ibc].BC_Data_Float[0],
                 &BC_Types[ibc].BC_Data_Float[1], &BC_Types[ibc].BC_Data_Float[2],
                 &BC_Types[ibc].BC_Data_Float[3], &BC_Types[ibc].BC_Data_Float[4],
                 &BC_Types[ibc].BC_Data_Float[5], &BC_Types[ibc].BC_Data_Float[6],
                 &BC_Types[ibc].BC_Data_Float[7], &BC_Types[ibc].BC_Data_Float[8]) != 9) {
        sprintf(err_msg, "%s: Expected 9 flts for BC %s.\n", yo, BC_Types[ibc].desc->name1);
        strcat(err_msg,
               "Be advised syntax to this boundary condition has changed. Consult GOMA manual.\n");
        GOMA_EH(GOMA_ERROR, err_msg);
      } else { /* 1*/
        for (i = 0; i < 9; i++)
          SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].BC_Data_Float[i]);

        if (fscanf(ifp, "%d", &BC_Types[ibc].BC_Data_Int[2]) != 1) {
          BC_Types[ibc].BC_Data_Int[2] = -1;
        } else
          SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_Data_Int[2]);
        /*   outside slip coeff, gas_factor,  contact_friction   */
        if (fscanf(ifp, "%lf", &BC_Types[ibc].BC_Data_Float[9]) != 1) {
          BC_Types[ibc].BC_Data_Float[9] = 1.0 / sqrt(LITTLE_PENALTY * BIG_PENALTY);
        } else
          SPF(endofstring(echo_string), " %lf", BC_Types[ibc].BC_Data_Float[9]);
        if (fscanf(ifp, "%lf", &BC_Types[ibc].BC_Data_Float[10]) != 1) {
          BC_Types[ibc].BC_Data_Float[10] = 8.0;
        } else
          SPF(endofstring(echo_string), " %lf", BC_Types[ibc].BC_Data_Float[10]);
        if (fscanf(ifp, "%lf", &BC_Types[ibc].BC_Data_Float[11]) != 1) {
          BC_Types[ibc].BC_Data_Float[11] = 1e-6;
        } else
          SPF(endofstring(echo_string), " %lf", BC_Types[ibc].BC_Data_Float[11]);
      }
      BC_Types[ibc].BC_Data_Int[3] = (BC_Types[ibc].BC_Name == SHIK_DIRICH_ROLL_BC) ? TRUE : FALSE;

      break;

      /*
       * Fall through for all cases which require ten floating point
       * values as data input
       */
    case MOVING_CA_BC:
    case REP_FORCE_ROLL_BC:
    case REP_FORCE_ROLL_RS_BC:
    case CAP_REPULSE_TABLE_BC:
      if (fscanf(ifp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &BC_Types[ibc].BC_Data_Float[0],
                 &BC_Types[ibc].BC_Data_Float[1], &BC_Types[ibc].BC_Data_Float[2],
                 &BC_Types[ibc].BC_Data_Float[3], &BC_Types[ibc].BC_Data_Float[4],
                 &BC_Types[ibc].BC_Data_Float[5], &BC_Types[ibc].BC_Data_Float[6],
                 &BC_Types[ibc].BC_Data_Float[7], &BC_Types[ibc].BC_Data_Float[8],
                 &BC_Types[ibc].BC_Data_Float[9]) != 10) {
        sr = sprintf(err_msg, "%s: Expected 10(!) flts for %s.", yo, BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      for (i = 0; i < 10; i++)
        SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].BC_Data_Float[i]);

      if (fscanf(ifp, "%d", &BC_Types[ibc].BC_Data_Int[2]) != 1) {
        BC_Types[ibc].BC_Data_Int[2] = -1;
      } else
        SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_Data_Int[2]);

      break;
      /*
       * Fall through for all cases which require twelve floating point
       * values as data input
       */
    case E_ER_PLANEWAVE_BC:
    case E_EI_PLANEWAVE_BC:
      if (fscanf(ifp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                 &BC_Types[ibc].BC_Data_Float[0], &BC_Types[ibc].BC_Data_Float[1],
                 &BC_Types[ibc].BC_Data_Float[2], &BC_Types[ibc].BC_Data_Float[3],
                 &BC_Types[ibc].BC_Data_Float[4], &BC_Types[ibc].BC_Data_Float[5],
                 &BC_Types[ibc].BC_Data_Float[6], &BC_Types[ibc].BC_Data_Float[7],
                 &BC_Types[ibc].BC_Data_Float[8], &BC_Types[ibc].BC_Data_Float[9],
                 &BC_Types[ibc].BC_Data_Float[10], &BC_Types[ibc].BC_Data_Float[11]) != 12) {
        sr = sprintf(err_msg, "%s: Expected 12(!) flts for %s.", yo, BC_Types[ibc].desc->name1);
        GOMA_EH(-1, err_msg);
      }
      for (i = 0; i < 12; i++)
        SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].BC_Data_Float[i]);

      if (fscanf(ifp, "%d", &BC_Types[ibc].BC_Data_Int[2]) != 1) {
        BC_Types[ibc].BC_Data_Int[2] = -1;
      } else
        SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_Data_Int[2]);

      break;

    case CAP_REPULSE_USER_BC:
    case QRAD_REPULSE_ROLL_BC:

      if (fscanf(ifp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                 &BC_Types[ibc].BC_Data_Float[0], &BC_Types[ibc].BC_Data_Float[1],
                 &BC_Types[ibc].BC_Data_Float[2], &BC_Types[ibc].BC_Data_Float[3],
                 &BC_Types[ibc].BC_Data_Float[4], &BC_Types[ibc].BC_Data_Float[5],
                 &BC_Types[ibc].BC_Data_Float[6], &BC_Types[ibc].BC_Data_Float[7],
                 &BC_Types[ibc].BC_Data_Float[8], &BC_Types[ibc].BC_Data_Float[9],
                 &BC_Types[ibc].BC_Data_Float[10], &BC_Types[ibc].BC_Data_Float[11],
                 &BC_Types[ibc].BC_Data_Float[12], &BC_Types[ibc].BC_Data_Float[13],
                 &BC_Types[ibc].BC_Data_Float[14]) != 15) {
        sr = sprintf(err_msg, "%s: Expected 15 flts for %s.", yo, BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      for (i = 0; i < 15; i++)
        SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].BC_Data_Float[i]);

      /* Try reading the optional integer. */
      if (fscanf(ifp, "%d", &BC_Types[ibc].BC_Data_Int[2]) != 1) {
        BC_Types[ibc].BC_Data_Int[2] = -1;
      } else
        SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_Data_Int[2]);

      break;
    case CAP_REPULSE_ROLL_BC:

      if (fscanf(ifp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                 &BC_Types[ibc].BC_Data_Float[0], &BC_Types[ibc].BC_Data_Float[1],
                 &BC_Types[ibc].BC_Data_Float[2], &BC_Types[ibc].BC_Data_Float[3],
                 &BC_Types[ibc].BC_Data_Float[4], &BC_Types[ibc].BC_Data_Float[5],
                 &BC_Types[ibc].BC_Data_Float[6], &BC_Types[ibc].BC_Data_Float[7],
                 &BC_Types[ibc].BC_Data_Float[8], &BC_Types[ibc].BC_Data_Float[9],
                 &BC_Types[ibc].BC_Data_Float[10], &BC_Types[ibc].BC_Data_Float[11]) != 12) {
        sr = sprintf(err_msg, "%s: Expected 12 required flts for %s.", yo,
                     BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      for (i = 0; i < 12; i++)
        SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].BC_Data_Float[i]);
      /* Try reading the optional integer. */
      if (fscanf(ifp, "%d", &BC_Types[ibc].BC_Data_Int[0]) != 1) {
        BC_Types[ibc].BC_Data_Int[0] = -1;
      } else
        SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_Data_Int[0]);

      if (fscanf(ifp, "%lf %lf %lf %lf", &BC_Types[ibc].BC_Data_Float[12],
                 &BC_Types[ibc].BC_Data_Float[13], &BC_Types[ibc].BC_Data_Float[14],
                 &BC_Types[ibc].BC_Data_Float[15]) != 4) {
        sr = sprintf(err_msg, "%s: Expected 12 required flts for %s.", yo,
                     BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
        BC_Types[ibc].BC_Data_Float[12] = 1.;
        BC_Types[ibc].BC_Data_Float[13] = 2.;
        BC_Types[ibc].BC_Data_Float[14] = 0.;
        BC_Types[ibc].BC_Data_Float[15] = 1.;
      }

      /* Try reading the optional integer. */
      if (fscanf(ifp, "%d", &BC_Types[ibc].BC_Data_Int[2]) != 1) {
        BC_Types[ibc].BC_Data_Int[2] = -1;
      } else
        SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_Data_Int[2]);

      break;
      /*
       * for USER DEFINED BOUNDARY CONDITIONS
       * Fall through for all cases which require an arbitrary number of constants
       */
    case FLOW_PRESS_USER_BC:
      GOMA_WH(-1, " FLOW_PRESS_USER is no longer used.  Use PRESSURE_USER instead.");
      /* Fall through */
    case UVARY_BC:
    case VVARY_BC:
    case WVARY_BC:
    case U_PARABOLA_BC:
    case V_PARABOLA_BC:
    case W_PARABOLA_BC:
    case FILLET_BC:
    case DOUBLE_RAD_BC:
    case ROLL_FLUID_BC:
    case SPLINEX_BC:
    case SPLINEY_BC:
    case SPLINEZ_BC:
    case SPLINE_BC:
    case SPLINEX_RS_BC:
    case SPLINEY_RS_BC:
    case SPLINEZ_RS_BC:
    case SPLINE_RS_BC:
    case T_USER_BC:
    case UUSER_BC:
    case VUSER_BC:
    case WUSER_BC:
    case UUSER_COLLOC_BC:
    case VUSER_COLLOC_BC:
    case WUSER_COLLOC_BC:
    case QUSER_BC:
    case DX_USER_BC:
    case DY_USER_BC:
    case DZ_USER_BC:
    case P_LIQ_USER_BC:
    case SH_P_OPEN_USER_BC:
    case Q_LASER_WELD_BC:
    case Q_VAPOR_BC:
    case LS_QLASER_BC:
    case LS_QVAPOR_BC:
    case FORCE_USER_BC:
    case FORCE_USER_SIC_BC:
    case FORCE_USER_RS_BC:
    case PRESSURE_USER_BC:
    case VAR_CA_USER_BC:
    case KIN_CHEM_BC:
    case CURRENT_USER_BC:
    case CURRENT_USER_SIC_BC:
    case VOLT_USER_BC:
    case U_VES11_PARABOLA_BC:
    case U_VES12_PARABOLA_BC:
    case U_VES22_PARABOLA_BC:
    case U_VES13_PARABOLA_BC:
    case U_VES23_PARABOLA_BC:
    case U_VES33_PARABOLA_BC:
    case U_VES11_1_PARABOLA_BC:
    case U_VES12_1_PARABOLA_BC:
    case U_VES22_1_PARABOLA_BC:
    case U_VES13_1_PARABOLA_BC:
    case U_VES23_1_PARABOLA_BC:
    case U_VES33_1_PARABOLA_BC:
    case U_VES11_2_PARABOLA_BC:
    case U_VES12_2_PARABOLA_BC:
    case U_VES22_2_PARABOLA_BC:
    case U_VES13_2_PARABOLA_BC:
    case U_VES23_2_PARABOLA_BC:
    case U_VES33_2_PARABOLA_BC:
    case U_VES11_3_PARABOLA_BC:
    case U_VES12_3_PARABOLA_BC:
    case U_VES22_3_PARABOLA_BC:
    case U_VES13_3_PARABOLA_BC:
    case U_VES23_3_PARABOLA_BC:
    case U_VES33_3_PARABOLA_BC:
    case U_VES11_4_PARABOLA_BC:
    case U_VES12_4_PARABOLA_BC:
    case U_VES22_4_PARABOLA_BC:
    case U_VES13_4_PARABOLA_BC:
    case U_VES23_4_PARABOLA_BC:
    case U_VES33_4_PARABOLA_BC:
    case U_VES11_5_PARABOLA_BC:
    case U_VES12_5_PARABOLA_BC:
    case U_VES22_5_PARABOLA_BC:
    case U_VES13_5_PARABOLA_BC:
    case U_VES23_5_PARABOLA_BC:
    case U_VES33_5_PARABOLA_BC:
    case U_VES11_6_PARABOLA_BC:
    case U_VES12_6_PARABOLA_BC:
    case U_VES22_6_PARABOLA_BC:
    case U_VES13_6_PARABOLA_BC:
    case U_VES23_6_PARABOLA_BC:
    case U_VES33_6_PARABOLA_BC:
    case U_VES11_7_PARABOLA_BC:
    case U_VES12_7_PARABOLA_BC:
    case U_VES22_7_PARABOLA_BC:
    case U_VES13_7_PARABOLA_BC:
    case U_VES23_7_PARABOLA_BC:
    case U_VES33_7_PARABOLA_BC:
      num_const = read_constants(ifp, &(BC_Types[ibc].u_BC), NO_SPECIES);

      if (num_const < 0) {
        sr = sprintf(err_msg, "? user BC consts for %s on %sID=%d\n", BC_Types[ibc].desc->name1,
                     BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      BC_Types[ibc].len_u_BC = num_const;
      BC_Types[ibc].max_DFlt = num_const;

      for (i = 0; i < num_const; i++)
        SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].u_BC[i]);
      if (BC_Types[ibc].BC_Name == ROLL_FLUID_BC)
        BC_Types[ibc].BC_Data_Int[0] = (int)BC_Types[ibc].u_BC[9];

      break;

    case FRICTION_ACOUSTIC_BC:
    case FRICTION_ACOUSTIC_RS_BC:

      if (fscanf(ifp, "%lf %d", &BC_Types[ibc].BC_Data_Float[0], &BC_Types[ibc].BC_Data_Int[0]) !=
          2) {
        sr = sprintf(err_msg, "%s: Expected 1 flt, 1 int for %s.", yo, BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].BC_Data_Float[0]);
      SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_Data_Int[0]);

      num_const = read_constants(ifp, &(BC_Types[ibc].u_BC), NO_SPECIES);
      if (num_const < 0) {
        sr = sprintf(err_msg, "? user BC consts for %s on %sID=%d\n", BC_Types[ibc].desc->name1,
                     BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      BC_Types[ibc].len_u_BC = num_const;
      BC_Types[ibc].max_DFlt = num_const + 1;

      for (i = 0; i < num_const; i++)
        SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].u_BC[i]);

      break;

    case VELO_TANGENT_USER_BC:
      if (fscanf(ifp, "%d %lf %lf %lf", &BC_Types[ibc].BC_Data_Int[0],
                 &BC_Types[ibc].BC_Data_Float[0], &BC_Types[ibc].BC_Data_Float[1],
                 &BC_Types[ibc].BC_Data_Float[2]) != 4) {
        sr = sprintf(err_msg, "Expected 1 int, 3 flts for %s on %sID=%d\n",
                     BC_Types[ibc].desc->name1, BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_Data_Int[0]);
      for (i = 0; i < 3; i++)
        SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].BC_Data_Float[i]);

      num_const = read_constants(ifp, &(BC_Types[ibc].u_BC), NO_SPECIES);
      if (num_const < 0) {
        sr = sprintf(err_msg, "? user BC consts for %s on %sID=%d\n", BC_Types[ibc].desc->name1,
                     BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      BC_Types[ibc].len_u_BC = num_const;
      BC_Types[ibc].max_DFlt = num_const + 3;

      for (i = 0; i < num_const; i++)
        SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].u_BC[i]);

      break;
      /*
       * for USER DEFINED BOUNDARY CONDITIONS
       * Fall through for all cases which require an integer plus an
       * arbitrary number of double constants
       */
    case YFLUX_USER_BC:
    case YUSER_BC:
    case YFLUX_ALLOY_BC:
    case FEATURE_ROLLON_BC:
      if (fscanf(ifp, "%d", &BC_Types[ibc].BC_Data_Int[0]) != 1) {
        sr = sprintf(err_msg, "Expected 1 int for %s on %sID=%d\n", BC_Types[ibc].desc->name1,
                     BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      BC_Types[ibc].species_eq = BC_Types[ibc].BC_Data_Int[0];
      num_const = read_constants(ifp, &(BC_Types[ibc].u_BC), NO_SPECIES);
      if (num_const < 0) {
        sr = sprintf(err_msg, "? user BC consts for %s on %sID=%d\n", BC_Types[ibc].desc->name1,
                     BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      BC_Types[ibc].len_u_BC = num_const;
      BC_Types[ibc].max_DFlt = num_const;

      SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_Data_Int[0]);
      for (i = 0; i < num_const; i++)
        SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].u_BC[i]);

      break;

      /*for MODEL options with USER-DEFINED BOUNDARY CONDITION option
       * Fall through for all cases which require a Keyword  plus an
       * arbitrary number of double constants
       */
    case CA_EDGE_OR_FIX_BC:
      if (fscanf(ifp, "%80s", input) != 1) {
        sprintf(err_msg, "%s: Expected Keyword for %s.", yo, BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      if (!strcmp(BC_Types[ibc].desc->name1, "CA_EDGE_OR_FIX")) {
        if (!strcmp(input, "CIRCLE")) {
          BC_Types[ibc].BC_Data_Int[0] = CIRCLE;
        } else if (!strcmp(input, "USER")) {
          BC_Types[ibc].BC_Data_Int[0] = B_USER;
        } else {
          GOMA_EH(GOMA_ERROR, "I don't recognize your CA_EDGE_OR_FIX Keyword!");
        }
      }
      num_const = read_constants(ifp, &(BC_Types[ibc].u_BC), NO_SPECIES);
      if (num_const < 0) {
        sprintf(err_msg, "? user BC consts for %s on %sID=%d\n", BC_Types[ibc].desc->name1,
                BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      BC_Types[ibc].len_u_BC = num_const;
      BC_Types[ibc].max_DFlt = num_const;

      SPF(endofstring(echo_string), " %s", input);
      for (i = 0; i < num_const; i++)
        SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].u_BC[i]);

      break;

      /*
       * Fall through for all cases which require one integer as data input
       */
    case POROUS_CONV_BC:
    case YFLUX_SUS_BC:
    case KIN_DISPLACEMENT_PETROV_BC:
    case KIN_DISPLACEMENT_COLLOC_BC:
    case KIN_DISPLACEMENT_BC:
    case KIN_DISPLACEMENT_RS_BC:
    case PERIODIC_BC:
      if (fscanf(ifp, "%d", &BC_Types[ibc].BC_Data_Int[0]) != 1) {
        sprintf(err_msg, "Expected 1 int for %s on %sID=%d\n", BC_Types[ibc].desc->name1,
                BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      BC_Types[ibc].species_eq = BC_Types[ibc].BC_Data_Int[0];

      if (BC_Types[ibc].BC_Name == PERIODIC_BC)
        GOMA_EH(GOMA_ERROR, "Don't you wish ...there were PERIODIC_BC");

      SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_Data_Int[0]);
      if (BC_Types[ibc].BC_Name == KIN_DISPLACEMENT_BC ||
          BC_Types[ibc].BC_Name == KIN_DISPLACEMENT_COLLOC_BC ||
          BC_Types[ibc].BC_Name == KIN_DISPLACEMENT_PETROV_BC) {
        num_const = read_constants(ifp, &(BC_Types[ibc].u_BC), NO_SPECIES);
        BC_Types[ibc].len_u_BC = num_const;

        for (i = 0; i < num_const; i++)
          SPF(endofstring(echo_string), " %g", BC_Types[ibc].u_BC[i]);
      }

      break;

      /*
       * Fall through for all cases which require one integer and one
       * float as data input
       */
    case ELEC_TRACTION_BC:
    case ELEC_TRACTION_SOLID_BC:
    case SH_GAMMA1_DERIV_SYMM_BC:
    case SH_GAMMA2_DERIV_SYMM_BC:
    case DVZDR_ZERO_BC:
      if (fscanf(ifp, "%d %lf", &BC_Types[ibc].BC_Data_Int[0], &BC_Types[ibc].BC_Data_Float[0]) !=
          2) {
        sprintf(err_msg, "Expected 1 int and 1 float for %s on %sID=%d\n",
                BC_Types[ibc].desc->name1, BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      SPF(endofstring(echo_string), " %d %.4g", BC_Types[ibc].BC_Data_Int[0],
          BC_Types[ibc].BC_Data_Float[0]);

      if (BC_Types[ibc].BC_Name == DVZDR_ZERO_BC) {
        // specify the radial direction as the wall
        BC_Types[ibc].BC_Data_Float[1] = 0.0;
        BC_Types[ibc].BC_Data_Float[2] = 1.0;
        BC_Types[ibc].BC_Data_Float[3] = 0.0;
      }
      break;

      /*
       * Stefan flow and kinematic stefan flow take one int
       * and one string
       *  -> first int is the element block of the pertinent side
       *     to apply the boundary condition on
       *  -> string is the source term to move the boundary.
       */
    case SDC_STEFANFLOW_BC:
    case SDC_KIN_SF_BC:
      if (fscanf(ifp, "%d %s", &BC_Types[ibc].BC_Data_Int[0], ts) != 2) {
        sprintf(err_msg, "Expected 1 int and a string for %s on %sID=%d\n",
                BC_Types[ibc].desc->name1, BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      BC_Types[ibc].BC_Data_Int[1] = match_interface_source_string(ts);
      BC_Types[ibc].Storage_ID = BC_Types[ibc].BC_Data_Int[1];
      SPF(endofstring(echo_string), " %d %s", BC_Types[ibc].BC_Data_Int[0], ts);
      break;

      /*
       * Fall through for all cases which require one integer and one
       * floating point values as data input, and one floating
       * point value and one integer value as optional input.
       */
    case Y_BC:
    case Y_DISCONTINUOUS_BC:
    case YFLUX_CONST_BC:
    case YTOTALFLUX_CONST_BC:
    case KINEMATIC_SPECIES_BC:
    case SURFACE_CHARGE_BC:
    case FICK_CHRGD_SURF_GRAD_BC:

      if (fscanf(ifp, "%d %lf", &BC_Types[ibc].BC_Data_Int[0], &BC_Types[ibc].BC_Data_Float[0]) !=
          2) {
        sprintf(err_msg, "Expected 1 int, 1 flt for %s on %sID=%d\n", BC_Types[ibc].desc->name1,
                BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      SPF(endofstring(echo_string), " %d %.4g", BC_Types[ibc].BC_Data_Int[0],
          BC_Types[ibc].BC_Data_Float[0]);

      if (!strcmp(BC_Types[ibc].Set_Type, "NS")) {
        if (fscanf(ifp, "%lf", &BC_Types[ibc].BC_relax) != 1) {
          BC_Types[ibc].BC_relax = -1.0;
        } else {

          SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].BC_relax);

          if (fscanf(ifp, "%d", &BC_Types[ibc].BC_EBID_Apply) != 1) {
            BC_Types[ibc].BC_EBID_Apply = -1;
          } else
            SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_EBID_Apply);
        }
      }
      BC_Types[ibc].species_eq = BC_Types[ibc].BC_Data_Int[0];
      break;

      /*
       * Fall through for all cases which require one integer and seven
       * floating point values as data input
       */
    case CURRENT_BV_BC:
      if (fscanf(ifp, "%d %lf %lf %lf %lf %lf %lf %lf", &BC_Types[ibc].BC_Data_Int[0],
                 &BC_Types[ibc].BC_Data_Float[0], &BC_Types[ibc].BC_Data_Float[1],
                 &BC_Types[ibc].BC_Data_Float[2], &BC_Types[ibc].BC_Data_Float[3],
                 &BC_Types[ibc].BC_Data_Float[4], &BC_Types[ibc].BC_Data_Float[5],
                 &BC_Types[ibc].BC_Data_Float[6]) != 8) {
        sr = sprintf(err_msg, "Expected 1 int, 7 flts for %s on %sID=%d\n",
                     BC_Types[ibc].desc->name1, BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      BC_Types[ibc].species_eq = BC_Types[ibc].BC_Data_Int[0];

      SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_Data_Int[0]);
      SPF_DBL_VEC(endofstring(echo_string), 7, BC_Types[ibc].BC_Data_Float);
      BC_Types[ibc].max_DFlt = 7;
      break;

      /*
       * Fall through for all cases which require one integer and eight
       * floating point values as data input
       */
    case CURRENT_ORR_BC:
    case YFLUX_H2O_ANODE_BC:
    case YFLUX_H2O_CATHODE_BC:
      if (fscanf(ifp, "%d %lf %lf %lf %lf %lf %lf %lf %lf", &BC_Types[ibc].BC_Data_Int[0],
                 &BC_Types[ibc].BC_Data_Float[0], &BC_Types[ibc].BC_Data_Float[1],
                 &BC_Types[ibc].BC_Data_Float[2], &BC_Types[ibc].BC_Data_Float[3],
                 &BC_Types[ibc].BC_Data_Float[4], &BC_Types[ibc].BC_Data_Float[5],
                 &BC_Types[ibc].BC_Data_Float[6], &BC_Types[ibc].BC_Data_Float[7]) != 9) {
        sr = sprintf(err_msg, "Expected 1 int, 8 flts for %s on %sID=%d\n",
                     BC_Types[ibc].desc->name1, BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      BC_Types[ibc].species_eq = BC_Types[ibc].BC_Data_Int[0];

      SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_Data_Int[0]);
      SPF_DBL_VEC(endofstring(echo_string), 8, BC_Types[ibc].BC_Data_Float);
      BC_Types[ibc].max_DFlt = 8;
      break;

      /*
       * Fall through for all cases which require one integer and nine
       * floating point values as data input
       */
    case YFLUX_ORR_BC:
    case CURRENT_HOR_BC:
    case YFLUX_BV_BC:
      if (fscanf(ifp, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf", &BC_Types[ibc].BC_Data_Int[0],
                 &BC_Types[ibc].BC_Data_Float[0], &BC_Types[ibc].BC_Data_Float[1],
                 &BC_Types[ibc].BC_Data_Float[2], &BC_Types[ibc].BC_Data_Float[3],
                 &BC_Types[ibc].BC_Data_Float[4], &BC_Types[ibc].BC_Data_Float[5],
                 &BC_Types[ibc].BC_Data_Float[6], &BC_Types[ibc].BC_Data_Float[7],
                 &BC_Types[ibc].BC_Data_Float[8]) != 10) {
        sr = sprintf(err_msg, "Expected 1 int, 9 flts for %s on %sID=%d\n",
                     BC_Types[ibc].desc->name1, BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      BC_Types[ibc].species_eq = BC_Types[ibc].BC_Data_Int[0];

      SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_Data_Int[0]);
      SPF_DBL_VEC(endofstring(echo_string), 9, BC_Types[ibc].BC_Data_Float);
      BC_Types[ibc].max_DFlt = 9;
      break;

      /*
       * Fall through for all cases which require one integer and ten
       * floating point values as data input
       */
    case YFLUX_HOR_BC:
      if (fscanf(ifp, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &BC_Types[ibc].BC_Data_Int[0],
                 &BC_Types[ibc].BC_Data_Float[0], &BC_Types[ibc].BC_Data_Float[1],
                 &BC_Types[ibc].BC_Data_Float[2], &BC_Types[ibc].BC_Data_Float[3],
                 &BC_Types[ibc].BC_Data_Float[4], &BC_Types[ibc].BC_Data_Float[5],
                 &BC_Types[ibc].BC_Data_Float[6], &BC_Types[ibc].BC_Data_Float[7],
                 &BC_Types[ibc].BC_Data_Float[8], &BC_Types[ibc].BC_Data_Float[9]) != 11) {
        sr = sprintf(err_msg, "Expected 1 int, 10 flts for %s on %sID=%d\n",
                     BC_Types[ibc].desc->name1, BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      BC_Types[ibc].species_eq = BC_Types[ibc].BC_Data_Int[0];

      SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_Data_Int[0]);
      SPF_DBL_VEC(endofstring(echo_string), 10, BC_Types[ibc].BC_Data_Float);
      BC_Types[ibc].max_DFlt = 10;
      break;

    case YFLUX_BV2_BC:
    case CURRENT_BV2_BC: /* RSL 1/15/01 */
      if (fscanf(ifp, "%d %lf %lf %lf %lf %lf %lf %lf %lf", &BC_Types[ibc].BC_Data_Int[0],
                 &BC_Types[ibc].BC_Data_Float[0], &BC_Types[ibc].BC_Data_Float[1],
                 &BC_Types[ibc].BC_Data_Float[2], &BC_Types[ibc].BC_Data_Float[3],
                 &BC_Types[ibc].BC_Data_Float[4], &BC_Types[ibc].BC_Data_Float[5],
                 &BC_Types[ibc].BC_Data_Float[6], &BC_Types[ibc].BC_Data_Float[7]) != 9) {
        sr = sprintf(err_msg, "Expected 1 int, 8 flts for %s on %sID=%d\n",
                     BC_Types[ibc].desc->name1, BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      BC_Types[ibc].species_eq = BC_Types[ibc].BC_Data_Int[0];

      SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_Data_Int[0]);
      SPF_DBL_VEC(endofstring(echo_string), 9, BC_Types[ibc].BC_Data_Float);
      BC_Types[ibc].max_DFlt = 9;

      break;

    case YFLUX_BC:
    case YFLUX_NI_BC:
    case LS_YFLUX_BC:
    case CURRENT_NI_BC: /* RSL 3/9/01 */
      if (fscanf(ifp, "%d %lf %lf", &BC_Types[ibc].BC_Data_Int[0], &BC_Types[ibc].BC_Data_Float[0],
                 &BC_Types[ibc].BC_Data_Float[1]) != 3) {
        sprintf(err_msg, "Expected 1 int, 2 flts for %s on %sID=%d\n", BC_Types[ibc].desc->name1,
                BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      BC_Types[ibc].species_eq = BC_Types[ibc].BC_Data_Int[0];

      SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_Data_Int[0]);
      SPF_DBL_VEC(endofstring(echo_string), 2, BC_Types[ibc].BC_Data_Float);
      BC_Types[ibc].max_DFlt = 2;

      break;

      /*
       * Fall through for all cases which require two integers as data input
       * Plus an optional float
       */
    case SOLID_FLUID_BC:
    case FLUID_SOLID_BC:
    case SOLID_FLUID_RS_BC:
    case FLUID_SOLID_RS_BC:
    case SOLID_FLUID_CONTACT_BC:
    case FLUID_SOLID_CONTACT_BC:
    case SOLID_LAGRANGE_MULT_BC:
    case BAAIJENS_SOLID_FLUID_BC:
    case BAAIJENS_FLUID_SOLID_BC:
    case LAGRANGE_NO_SLIP_BC:
    case POROUS_PRESSURE_BC:
    case POROUS_PRESSURE_LUB_BC:
    case DARCY_CONTINUOUS_BC:
    case DARCY_LUB_BC:
    case NO_SLIP_BC:
    case NO_SLIP_RS_BC:
    case VELO_TANGENT_SOLID_BC:
    case VELO_NORMAL_SOLID_BC:
    case SURFACE_ELECTRIC_FIELD_BC:
    case SURFACE_ACOUSTIC_VELOCITY_BC:
    case SURFACE_USER_SHELL_BC:
    case SURFACE_LUBRICATION_BC:
    case LUBP_SH_FP_MATCH_BC:
    case LUBP_SH_FP_FLUX_BC:
    case T_CONTACT_RESIS_BC:
    case T_CONTACT_RESIS_2_BC:
    case LIGHTP_JUMP_BC:
    case LIGHTM_JUMP_BC:
    case LIGHTP_JUMP_2_BC:
    case LIGHTM_JUMP_2_BC:

      if (fscanf(ifp, "%d %d", &BC_Types[ibc].BC_Data_Int[0], &BC_Types[ibc].BC_Data_Int[1]) != 2) {
        sr = sprintf(err_msg, "Expected 2 ints for %s on %sID=%d\n", BC_Types[ibc].desc->name1,
                     BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      SPF(endofstring(echo_string), " %d %d", BC_Types[ibc].BC_Data_Int[0],
          BC_Types[ibc].BC_Data_Int[1]);
      /* add optional scaling float for fluid-structure bcs */

      if (!strcmp(BC_Types[ibc].desc->name1, "VELO_TANGENT_SOLID")) {
        if (fscanf(ifp, "%d ", &BC_Types[ibc].BC_Data_Int[2]) != 1) {
          BC_Types[ibc].BC_Data_Int[2] = 0;
        } else
          SPF(endofstring(echo_string), " %d ", BC_Types[ibc].BC_Data_Int[2]);
        /* these also need to be initialized */
        BC_Types[ibc].BC_Data_Float[0] = 1.; /*beta*/
        BC_Types[ibc].BC_Data_Float[1] = 0.;
      }

      if (!strcmp(BC_Types[ibc].desc->name1, "DARCY_CONTINUOUS")) {
        if (fscanf(ifp, "%lf %lf", &BC_Types[ibc].BC_Data_Float[0],
                   &BC_Types[ibc].BC_Data_Float[1]) != 2) {
          /* Default scaling factor to 0.0 */
          SPF(endofstring(echo_string), " %s", "(DARCY_CONTINUOUS length scale defaults to 0.0)");

          BC_Types[ibc].BC_Data_Float[0] = 0.0; /*the length scale */
          BC_Types[ibc].BC_Data_Float[1] = 0.0; /*the threshold time */
        } else
          SPF(endofstring(echo_string), " %.4g %.4g", BC_Types[ibc].BC_Data_Float[0],
              BC_Types[ibc].BC_Data_Float[1]);
      }
      if (!strcmp(BC_Types[ibc].desc->name1, "T_CONTACT_RESIS") ||
          !strcmp(BC_Types[ibc].desc->name1, "T_CONTACT_RESIS_2")) {
        if (fscanf(ifp, "%lf", &BC_Types[ibc].BC_Data_Float[0]) != 1) {
          /* Default scaling factor to 0.0 */
          SPF(endofstring(echo_string), " %s", "(Warning, no thermal contact resistance)");

          BC_Types[ibc].BC_Data_Float[0] = 0.0;
        } else
          SPF(endofstring(echo_string), " %.4g ", BC_Types[ibc].BC_Data_Float[0]);
      }
      if (!strcmp(BC_Types[ibc].desc->name1, "SOLID_FLUID_RS") ||
          !strcmp(BC_Types[ibc].desc->name1, "SOLID_FLUID") ||
          !strcmp(BC_Types[ibc].desc->name1, "FLUID_SOLID") ||
          !strcmp(BC_Types[ibc].desc->name1, "FLUID_SOLID_RS")) {
        if (fscanf(ifp, "%lf", &BC_Types[ibc].BC_Data_Float[0]) != 1) {
          /* Default scaling factor to 1.0 */
          SPF(endofstring(echo_string), " %s", "(Solid-Fluid Scaling Factor defaults to 1.0)");

          BC_Types[ibc].BC_Data_Float[0] = 1.0;
        } else
          SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].BC_Data_Float[0]);
      }
      if (!strcmp(BC_Types[ibc].desc->name1, "FLUID_SOLID") && Linear_Solver == FRONT)
        GOMA_EH(GOMA_ERROR,
                "Mucho Problemo with FLUID_SOLID boundary condition and frontal solver. USE lu");
      if ((!strcmp(BC_Types[ibc].desc->name1, "SOLID_FLUID") ||
           !strcmp(BC_Types[ibc].desc->name1, "SOLID_FLUID_RS")) &&
          Linear_Solver == FRONT) {
        SPF(endofstring(echo_string), "\n  %s",
            "Warning: Poquito problemo con el SOLID_FLUID BC y frontal solver.");
        SPF(endofstring(echo_string), " %s",
            "Llame el senor Schunk @ 5058458991.  Gracias, que tenga una buena dia.");
      }

      break;

      /*
       * Fall through cases for 2 integers, then 1 float, followed by optional int. and  float
       */

    case VELO_SLIP_SOLID_BC:
      if (fscanf(ifp, "%d %d %lf", &BC_Types[ibc].BC_Data_Int[0], &BC_Types[ibc].BC_Data_Int[1],
                 &BC_Types[ibc].BC_Data_Float[0]) != 3)

      {
        sr = sprintf(err_msg, "Expected 2 ints and 1 float for %s on %sID=%d\n",
                     BC_Types[ibc].desc->name1, BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      SPF(endofstring(echo_string), " %d %d %.4g", BC_Types[ibc].BC_Data_Int[0],
          BC_Types[ibc].BC_Data_Int[1], BC_Types[ibc].BC_Data_Float[0]);

      if (fscanf(ifp, "%d %lf", &BC_Types[ibc].BC_Data_Int[2], &BC_Types[ibc].BC_Data_Float[1]) !=
          2) {
        BC_Types[ibc].BC_Data_Int[2] = -1;
        BC_Types[ibc].BC_Data_Float[1] = 0.;
      } else
        SPF(endofstring(echo_string), " %d %.4g", BC_Types[ibc].BC_Data_Int[2],
            BC_Types[ibc].BC_Data_Float[1]);

      break;

    case QSIDE_LS_BC:

      if (fscanf(ifp, "%d %d %lf %lf", &BC_Types[ibc].BC_Data_Int[0], &BC_Types[ibc].BC_Data_Int[1],
                 &BC_Types[ibc].BC_Data_Float[0], &BC_Types[ibc].BC_Data_Float[1]) != 4)

      {
        sr = sprintf(err_msg, "Expected 2 ints and 2 float for %s on %sID=%d\n",
                     BC_Types[ibc].desc->name1, BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      SPF(endofstring(echo_string), " %d %d %.4g %.4g", BC_Types[ibc].BC_Data_Int[0],
          BC_Types[ibc].BC_Data_Int[1], BC_Types[ibc].BC_Data_Float[0],
          BC_Types[ibc].BC_Data_Float[1]);
      break;

      /*
       * Fall through cases for 2 integers, then 3 floats
       */
    case POROUS_LIQ_PRESSURE_FILL_BC:
      if (fscanf(ifp, "%d %d %lf %lf %lf", &BC_Types[ibc].BC_Data_Int[0],
                 &BC_Types[ibc].BC_Data_Int[1], &BC_Types[ibc].BC_Data_Float[0],
                 &BC_Types[ibc].BC_Data_Float[1], &BC_Types[ibc].BC_Data_Float[2]) != 5)

      {
        sr = SPF(err_msg, "Expected 2 ints and 3 floats for %s on %sID=%d\n",
                 BC_Types[ibc].desc->name1, BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      SPF(endofstring(echo_string), " %d %d", BC_Types[ibc].BC_Data_Int[0],
          BC_Types[ibc].BC_Data_Int[1]);
      for (i = 0; i < 3; i++)
        SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].BC_Data_Float[i]);

      break;

      /*
       * Fall through cases for 2 integers, then 4 floats
       */
    case POR_LIQ_FLUX_FILL_BC:
      if (fscanf(ifp, "%d %d %lf %lf %lf %lf", &BC_Types[ibc].BC_Data_Int[0],
                 &BC_Types[ibc].BC_Data_Int[1], &BC_Types[ibc].BC_Data_Float[0],
                 &BC_Types[ibc].BC_Data_Float[1], &BC_Types[ibc].BC_Data_Float[2],
                 &BC_Types[ibc].BC_Data_Float[3]) != 6)

      {
        sr = SPF(err_msg, "Expected 2 ints and 4 floats for %s on %sID=%d\n",
                 BC_Types[ibc].desc->name1, BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      SPF(endofstring(echo_string), " %d %d", BC_Types[ibc].BC_Data_Int[0],
          BC_Types[ibc].BC_Data_Int[1]);
      for (i = 0; i < 4; i++)
        SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].BC_Data_Float[i]);

      break;

      /*
       * Fall through for all cases which require three integers as
       * data input
       */
    case P_EQUIL_BC:
      if (fscanf(ifp, "%d %d %d", &BC_Types[ibc].BC_Data_Int[0], &BC_Types[ibc].BC_Data_Int[1],
                 &BC_Types[ibc].BC_Data_Int[2]) != 3) {
        SPF(err_msg, "Expected 3 ints for %s on %sID=%d\n", BC_Types[ibc].desc->name1,
            BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      for (i = 0; i < 3; i++)
        SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_Data_Int[i]);

      break;

      /*
       * Fall through for all cases which require three integers
       * and one float
       */
    case VN_POROUS_BC:
    case VP_EQUIL_BC:
    case VL_EQUIL_PRXN_BC:
    case IS_EQUIL_PRXN_BC:
      if (fscanf(ifp, "%d %d %d %lf", &BC_Types[ibc].BC_Data_Int[0], &BC_Types[ibc].BC_Data_Int[1],
                 &BC_Types[ibc].BC_Data_Int[2], &BC_Types[ibc].BC_Data_Float[0]) != 4) {
        SPF(err_msg, "Expected 3 ints, 1 flt for %s on %sID=%d\n", BC_Types[ibc].desc->name1,
            BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      BC_Types[ibc].species_eq = BC_Types[ibc].BC_Data_Int[0];
      if (BC_Types[ibc].BC_Name == VL_EQUIL_PRXN_BC) {
        assign_global_species_var_type(SPECIES_MASS_FRACTION, FALSE);
      }
      if (BC_Types[ibc].BC_Name == IS_EQUIL_PRXN_BC) {
        assign_global_species_var_type(SPECIES_CONCENTRATION, FALSE);
      }
      if (BC_Types[ibc].species_eq == 0 && (BC_Types[ibc].BC_Name == IS_EQUIL_PRXN_BC ||
                                            BC_Types[ibc].BC_Name == VL_EQUIL_PRXN_BC)) {
        IntSrc_BCID[Num_Interface_Srcs] = BC_Types[ibc].BC_ID;
        Num_Interface_Srcs++;
      }

      for (i = 0; i < 3; i++)
        SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_Data_Int[i]);
      SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].BC_Data_Float[0]);
      break;

      /*
       * Fall through for all cases which require three integers
       * and two float
       */
    case POROUS_GAS_BC:
    case YFLUX_DISC_RXN_BC:
      if (fscanf(ifp, "%d %d %d %lf %lf", &BC_Types[ibc].BC_Data_Int[0],
                 &BC_Types[ibc].BC_Data_Int[1], &BC_Types[ibc].BC_Data_Int[2],
                 &BC_Types[ibc].BC_Data_Float[0], &BC_Types[ibc].BC_Data_Float[1]) != 5) {
        sprintf(err_msg, "Expected 3 ints, 2 flts for %s on %sID=%d\n", BC_Types[ibc].desc->name1,
                BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      for (i = 0; i < 3; i++)
        SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_Data_Int[i]);
      SPF(endofstring(echo_string), " %.4g %.4g", BC_Types[ibc].BC_Data_Float[0],
          BC_Types[ibc].BC_Data_Float[1]);

      break;

      /*
       * Fall through for all cases which require one integer and three
       * floating point values as data input
       */

    case VELO_TANGENT_BC:
    case VELO_STREAMING_BC:
    case LATENT_HEAT_BC:
      if (fscanf(ifp, "%d %lf %lf %lf", &BC_Types[ibc].BC_Data_Int[0],
                 &BC_Types[ibc].BC_Data_Float[0], &BC_Types[ibc].BC_Data_Float[1],
                 &BC_Types[ibc].BC_Data_Float[2]) != 4) {
        sr = sprintf(err_msg, "Expected 1 int, 3 flts for %s on %sID=%d\n",
                     BC_Types[ibc].desc->name1, BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      /* LATENT_HEAT is a leak condition, so no species designation
         is needed. */
      /* BC_Types[ibc].species_eq = BC_Types[ibc].BC_Data_Int[0]; */

      BC_Types[ibc].max_DFlt = 3;
      SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_Data_Int[0]);
      for (i = 0; i < 3; i++)
        SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].BC_Data_Float[i]);

      if (fscanf(ifp, "%lf", &BC_Types[ibc].BC_Data_Float[3]) != 1) {
        BC_Types[ibc].BC_Data_Float[3] = 0.0;
        BC_Types[ibc].max_DFlt = 4;
      } else
        SPF(endofstring(echo_string), " %lf", BC_Types[ibc].BC_Data_Float[11]);

      break;

      /* fall through case with one int and 4 floating points */
    case POROUS_FLUX_BC:
      if (fscanf(ifp, "%d %lf %lf %lf %lf", &BC_Types[ibc].BC_Data_Int[0],
                 &BC_Types[ibc].BC_Data_Float[0], &BC_Types[ibc].BC_Data_Float[1],
                 &BC_Types[ibc].BC_Data_Float[2], &BC_Types[ibc].BC_Data_Float[3]) != 5) {
        sr = sprintf(err_msg, "Expected 1 int, 4 flts for %s on %sID=%d\n",
                     BC_Types[ibc].desc->name1, BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      BC_Types[ibc].species_eq = BC_Types[ibc].BC_Data_Int[0];

      SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_Data_Int[0]);
      for (i = 0; i < 4; i++)
        SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].BC_Data_Float[i]);

      break;

      /*
       * Fall through for all cases which require one string and one integer
       * value as data input
       */
    case FIX_BC:
      /*
       *  Read the variable name to be fixed
       */
      if (fscanf(ifp, "%80s", input) != 1) {
        GOMA_EH(GOMA_ERROR, "Error reading equation type for FIX_BC");
      }

      /* loop through variable names and compare to input */

      eq_found = FALSE;
      for (k = 0; k < MAX_VARIABLE_TYPES + 1; k++) /* only check variables (including porosity) */
      {
        if (!strcmp(input, Var_Name[k].name1) || !strcmp(input, Var_Name[k].name2)) {
          BC_Types[ibc].BC_Data_Int[0] = Var_Name[k].Index;
          BC_Types[ibc].equation = EQ_Name[k].Index;
          eq_found = TRUE;
        }
      }
      if (!eq_found) {
        sr = SPF(err_msg, "? Variable name = \"%s\" for %s on %sID=%d\n", input,
                 BC_Types[ibc].desc->name1, BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      if (fscanf(ifp, "%d", &BC_Types[ibc].BC_Data_Int[1]) != 1) {
        sr = SPF(err_msg, "? 1 int species ID after %s (or %s) for %s on %sID=%d\n",
                 Var_Name[BC_Types[ibc].BC_Data_Int[0]].name1,
                 Var_Name[BC_Types[ibc].BC_Data_Int[0]].name2, BC_Types[ibc].desc->name1,
                 BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      BC_Types[ibc].species_eq = BC_Types[ibc].BC_Data_Int[1];
      BC_Types[ibc].BC_relax = -1.0;

      /*
       * Need to Allocate Space for the BC Description of this Condition.
       *
       * Now, though, rather than allocate directly into whichever
       *     BC_Types[say_what].desc
       * we'll create a more organized data structure and copy the pointer
       * into the motley place.
       */

      new_BC_Desc = ((struct BC_descriptions **)realloc(
          new_BC_Desc, (num_new_BC_Desc + 1) * sizeof(struct BC_descriptions *)));

      new_BC_Desc[num_new_BC_Desc] = alloc_BC_description(BC_Types[ibc].desc);

      BC_Types[ibc].desc = new_BC_Desc[num_new_BC_Desc];

      BC_Types[ibc].index_dad = num_new_BC_Desc;

      num_new_BC_Desc++;

      BC_Types[ibc].desc->equation = BC_Types[ibc].BC_Data_Int[0];
      BC_Types[ibc].desc->sens[BC_Types[ibc].BC_Data_Int[0]] = 1;

      SPF(endofstring(echo_string), " %s %d", input, BC_Types[ibc].BC_Data_Int[1]);
      break;

      /*
       *  Fall through case here which requires 1 keyword, 2 integers,
       *  and a float
       */

    case DISCONTINUOUS_VELO_BC:
    case LATENT_HEAT_INTERNAL_BC:

      if (fscanf(ifp, "%80s", input) != 1) {
        sprintf(err_msg, "%s: Expected Keyword for %s.", yo, BC_Types[ibc].desc->name1);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      /* loop through variable names and compare to input */

      if (!strcmp(BC_Types[ibc].desc->name1, "LATENT_HEAT_INTERNAL")) {
        if (!strcmp(input, "LIQUID_VAPOR")) {
          BC_Types[ibc].BC_Data_Int[0] = LIQUID_VAPOR;
        } else if (!strcmp(input, "SOLID_LIQUID")) {
          BC_Types[ibc].BC_Data_Int[0] = SOLID_LIQUID;
        } else {
          GOMA_EH(GOMA_ERROR, "I don't recognize your LATENT_HEAT_INTERNAL Keyword!");
        }

        if (fscanf(ifp, "%d %d", &BC_Types[ibc].BC_Data_Int[1], &BC_Types[ibc].BC_Data_Int[2]) !=
            2) {
          GOMA_EH(GOMA_ERROR, "The LATENT_HEAT_INTERNAL card requires 2 material block numbers");
        }

        if (fscanf(ifp, "%lf", &BC_Types[ibc].BC_Data_Float[0]) != 1) {
          GOMA_EH(GOMA_ERROR,
                  "The LATENT_HEAT_INTERNAL card requires one float after the block numbers");
        }
        SPF(endofstring(echo_string), " %s", input);
        SPF(endofstring(echo_string), " %d %d %.4g", BC_Types[ibc].BC_Data_Int[1],
            BC_Types[ibc].BC_Data_Int[2], BC_Types[ibc].BC_Data_Float[0]);
      }

      if (!strcmp(BC_Types[ibc].desc->name1, "DISCONTINUOUS_VELO")) {
        if (!strcmp(input, "EVAPORATION")) {
          BC_Types[ibc].BC_Data_Int[0] = EVAPORATION;
        } else if (!strcmp(input, "DISSOLUTION")) {
          BC_Types[ibc].BC_Data_Int[0] = DISSOLUTION;
        } else {
          GOMA_EH(GOMA_ERROR, "I don't recognize your DISCONTINUOUS Keyword!");
        }

        if (fscanf(ifp, "%d %d", &BC_Types[ibc].BC_Data_Int[1], &BC_Types[ibc].BC_Data_Int[2]) !=
            2) {
          GOMA_EH(GOMA_ERROR, "The DISCONTINUOUS_VELO card requires 2 Element block numbers");
        }
        SPF(endofstring(echo_string), " %s", input);
        SPF(endofstring(echo_string), " %d %d", BC_Types[ibc].BC_Data_Int[1],
            BC_Types[ibc].BC_Data_Int[2]);
      }

      break;

      /*
       * Fall through for all cases which require two ( string and integer )
       * and multiple floating point constants as data input
       */
    case GD_CONST_BC:
    case GD_LINEAR_BC:
    case GD_INVERSE_BC:
    case GD_PARAB_BC:
    case GD_PARAB_OFFSET_BC:
    case GD_POLYN_BC:
    case GD_TIME_BC:
    case GD_CIRC_BC:
    case GD_TABLE_BC:

      /*
       *  Read the Equation to be replaced by a Generalized Dirichlet Condition
       */
      if (fscanf(ifp, "%80s", input) != 1) {
        GOMA_EH(GOMA_ERROR, "Error reading equation type for GD_*_BC");
      }

      /* loop through equation names and compare to input */

      eq_found = FALSE;
      for (k = 0; !eq_found && k < Num_EQ_Names; k++) {
        if (!strcmp(input, EQ_Name[k].name1) || !strcmp(input, EQ_Name[k].name2)) {
          BC_Types[ibc].BC_Data_Int[0] = EQ_Name[k].Index;
          eq_found = TRUE;
          BC_Types[ibc].equation = EQ_Name[k].Index;
        }
      }
      if (!eq_found) {
        sr = SPF(err_msg, "Equation name = \"%s\" for %s on %sID=%d ???\n", input,
                 BC_Types[ibc].desc->name1, BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      /*
       * Read the species number
       */

      if (fscanf(ifp, "%d", &BC_Types[ibc].BC_Data_Int[1]) != 1) {
        sr = SPF(err_msg, "Expected 1 int species ID after %s for %s on %sID=%d\n", input,
                 BC_Types[ibc].desc->name1, BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      BC_Types[ibc].species_eq = BC_Types[ibc].BC_Data_Int[1];

      k--;
      SPF(endofstring(echo_string), " %s %d", EQ_Name[k].name1, BC_Types[ibc].BC_Data_Int[1]);

      /*
       *  Read the variable name which this condition uses
       */
      if (fscanf(ifp, "%80s", input) != 1) {
        GOMA_EH(GOMA_ERROR, "Error reading variable type for GD_BC");
      }

      /* loop through variable names and compare to input */

      eq_found = FALSE;
      for (k = 0; !eq_found && k < Num_Var_Names; k++) {
        if (!strcmp(input, Var_Name[k].name1) || !strcmp(input, Var_Name[k].name2)) {
          BC_Types[ibc].BC_Data_Int[2] = Var_Name[k].Index;
          eq_found = TRUE;
        }
      }

      /****     Adding a block since the GD_TIME names won't be in the
                      variable list ; constants are defined in the
                      GD block of rf_bc_const.h         ****/

      if (!eq_found) {
        if (strcmp(input, "LINEAR") == 0) {
          BC_Types[ibc].BC_Data_Int[2] = GD_TIME_LIN;
          eq_found = TRUE;
        } else if (strcmp(input, "EXPONENTIAL") == 0) {
          BC_Types[ibc].BC_Data_Int[2] = GD_TIME_EXP;
          eq_found = TRUE;
        } else if (strcmp(input, "SINUSOIDAL") == 0) {
          BC_Types[ibc].BC_Data_Int[2] = GD_TIME_SIN;
          eq_found = TRUE;
        } else if (strcmp(input, "TABLE") == 0) {
          BC_Types[ibc].BC_Data_Int[2] = GD_TIME_TABLE;
          eq_found = TRUE;
        }
      }

      /****         end of modified block   ******/

      if (!eq_found) {
        sr = SPF(err_msg, "Variable name = \"%s\" for %s on %sID=%d ???\n", input,
                 BC_Types[ibc].desc->name1, BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      SPF(endofstring(echo_string), " %s", input);

      if (fscanf(ifp, "%d", &BC_Types[ibc].BC_Data_Int[3]) != 1) {
        sr = SPF(err_msg, "Expecting 1 int for %s on %sID=%d.\n", BC_Types[ibc].desc->name1,
                 BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      SPF(endofstring(echo_string), " %d", BC_Types[ibc].BC_Data_Int[3]);

      /*
       * read the floating point values that correspond to this
       * boundary condition type
       */
      switch (BC_Types[ibc].BC_Name) {
      case GD_CONST_BC:
        if (fscanf(ifp, "%lf", &BC_Types[ibc].BC_Data_Float[0]) != 1) {
          sr = sprintf(err_msg, "Expected 1 flt for %s on %sID=%d.\n", BC_Types[ibc].desc->name1,
                       BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
          GOMA_EH(GOMA_ERROR, err_msg);
        }
        SPF_DBL_VEC(endofstring(echo_string), 1, BC_Types[ibc].BC_Data_Float);
        BC_Types[ibc].max_DFlt = 1;
        break;
      case GD_LINEAR_BC:
      case GD_INVERSE_BC:
        if (fscanf(ifp, "%lf %lf", &BC_Types[ibc].BC_Data_Float[0],
                   &BC_Types[ibc].BC_Data_Float[1]) != 2) {
          sr = sprintf(err_msg, "Expected 2 flts for %s on %sID=%d.\n", BC_Types[ibc].desc->name1,
                       BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
          GOMA_EH(GOMA_ERROR, err_msg);
        }

        BC_Types[ibc].max_DFlt = 2;
        SPF_DBL_VEC(endofstring(echo_string), 2, BC_Types[ibc].BC_Data_Float);
        break;
      case GD_PARAB_BC:
      case GD_CIRC_BC:
        if (fscanf(ifp, "%lf %lf %lf", &BC_Types[ibc].BC_Data_Float[0],
                   &BC_Types[ibc].BC_Data_Float[1], &BC_Types[ibc].BC_Data_Float[2]) != 3) {
          sr = sprintf(err_msg, "Expected 3 flts for %s on %sID=%d.\n", BC_Types[ibc].desc->name1,
                       BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
          GOMA_EH(GOMA_ERROR, err_msg);
        }

        BC_Types[ibc].max_DFlt = 3;
        SPF_DBL_VEC(endofstring(echo_string), 3, BC_Types[ibc].BC_Data_Float);

        break;
      case GD_PARAB_OFFSET_BC:
        if (fscanf(ifp, "%lf %lf %lf %lf", &BC_Types[ibc].BC_Data_Float[0],
                   &BC_Types[ibc].BC_Data_Float[1], &BC_Types[ibc].BC_Data_Float[2],
                   &BC_Types[ibc].BC_Data_Float[3]) != 4) {
          sr = sprintf(err_msg, "Expected 4 flts for %s on %sID=%d.\n", BC_Types[ibc].desc->name1,
                       BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
          GOMA_EH(GOMA_ERROR, err_msg);
        }

        BC_Types[ibc].max_DFlt = 4;
        SPF_DBL_VEC(endofstring(echo_string), 4, BC_Types[ibc].BC_Data_Float);

        break;

      case GD_POLYN_BC:
        if (fscanf(ifp, "%lf %lf %lf %lf %lf %lf %lf", &BC_Types[ibc].BC_Data_Float[0],
                   &BC_Types[ibc].BC_Data_Float[1], &BC_Types[ibc].BC_Data_Float[2],
                   &BC_Types[ibc].BC_Data_Float[3], &BC_Types[ibc].BC_Data_Float[4],
                   &BC_Types[ibc].BC_Data_Float[5], &BC_Types[ibc].BC_Data_Float[6]) != 7) {
          sr = sprintf(err_msg, "Expected 7 flts for %s on %sID=%d.\n", BC_Types[ibc].desc->name1,
                       BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
          GOMA_EH(GOMA_ERROR, err_msg);
        }

        BC_Types[ibc].max_DFlt = 7;
        SPF_DBL_VEC(endofstring(echo_string), 7, BC_Types[ibc].BC_Data_Float);
        break;

      case GD_TIME_BC:
        BC_Types[ibc].BC_Data_Int[3] = GD_TIME_BC;

        if (BC_Types[ibc].BC_Data_Int[2] != GD_TIME_TABLE) {

          double *gd_time_values = NULL;
          int lfdcount = read_constants(ifp, &gd_time_values, 0);

          if (lfdcount < 2 || lfdcount > 3) {
            sr = sprintf(err_msg, "Expected 2 or 3 flts for %s on %sID=%d.\n",
                         BC_Types[ibc].desc->name1, BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
            GOMA_EH(GOMA_ERROR, err_msg);
          }
          SPF_DBL_VEC(endofstring(echo_string), lfdcount, gd_time_values);

          for (k = 0; k < lfdcount; k++) {
            BC_Types[ibc].BC_Data_Float[k] = gd_time_values[k];
          }

          free(gd_time_values);

          BC_Types[ibc].BC_Data_Int[4] = 0;
          if (lfdcount == 3) {
            if (BC_Types[ibc].BC_Data_Float[2] < 0) {
              sr =
                  sprintf(err_msg, "Expected a positive value for maximum time on %s on %sID=%d.\n",
                          BC_Types[ibc].desc->name1, BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
              GOMA_EH(GOMA_ERROR, err_msg);
            }
            BC_Types[ibc].BC_Data_Int[4] = GD_TIME_MAX;
          }

          BC_Types[ibc].max_DFlt = lfdcount;
          SPF_DBL_VEC(endofstring(echo_string), lfdcount, BC_Types[ibc].BC_Data_Float);
        } else {
          if (num_BC_Tables == MAX_BC_TABLES) {
            GOMA_EH(GOMA_ERROR, "Maximum TABLE_BCs exceeded .");
          }

          BC_Tables[num_BC_Tables] = setup_gd_table_BC(ifp, input, &BC_Types[ibc], echo_string);
          BC_Types[ibc].table_index = num_BC_Tables++;
        }

        break;

      case GD_TABLE_BC:

        /* Read scaling factor */

        if (fscanf(ifp, "%lf", &BC_Types[ibc].BC_Data_Float[0]) != 1) {
          sr = sprintf(err_msg, "Expected 1 flts for %s on %sID=%d.\n", BC_Types[ibc].desc->name1,
                       BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
          GOMA_EH(GOMA_ERROR, err_msg);
          for (i = 0; i < 1; i++)
            SPF(endofstring(echo_string), " %.4g", BC_Types[ibc].BC_Data_Float[i]);
        }

        if (num_BC_Tables == MAX_BC_TABLES) {
          GOMA_EH(GOMA_ERROR, "Maximum TABLE_BCs exceeded .");
        }

        BC_Tables[num_BC_Tables] = setup_gd_table_BC(ifp, input, &BC_Types[ibc], echo_string);

        BC_Types[ibc].table_index = num_BC_Tables++;

        break;

      default:
        GOMA_EH(GOMA_ERROR, "Unimplemented GD type condition");
        break;
      }

      /*
       * Need to Allocate Space for the BC Description for the
       * current boundary condition
       *
       *  We do this by reallocating the pointer array to add an
       *  extra pointer onto the end of the current list. Then
       *  we malloc and initialize a new BC description structure
       *  and add it to the end of the list.
       */
      new_BC_Desc = (struct BC_descriptions **)realloc(
          new_BC_Desc, (num_new_BC_Desc + 1) * sizeof(struct BC_descriptions *));
      new_BC_Desc[num_new_BC_Desc] = alloc_BC_description(BC_Types[ibc].desc);

      /*
       *  Now fill in the new BC_Types condition and BC_description
       *  structure with information about the boundary condition
       */
      BC_Types[ibc].desc = new_BC_Desc[num_new_BC_Desc];
      BC_Types[ibc].index_dad = num_new_BC_Desc;
      BC_Types[ibc].desc->equation = BC_Types[ibc].BC_Data_Int[0];
      var = BC_Types[ibc].BC_Data_Int[2];
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
      if (var >= V_FIRST && var < V_LAST) {
        BC_Types[ibc].desc->sens[var] = 1;
      } else if (var >= MESH_POSITION1 && var <= MESH_POSITION3) {
        BC_Types[ibc].desc->sens[var - MESH_POSITION1 + MESH_DISPLACEMENT1] = 1;
      } else if (var >= SOLID_POSITION1 && var <= SOLID_POSITION3) {
        BC_Types[ibc].desc->sens[var - SOLID_POSITION1 + SOLID_DISPLACEMENT1] = 1;
      } else if (var >= D_VEL1_DT && var <= D_P_DT) {
        BC_Types[ibc].desc->sens[var - D_VEL1_DT] = 1;
      }

      /* Check to see if this condition implies rotation of mesh or momentum equations */

      eqn = BC_Types[ibc].BC_Data_Int[0];
      if ((eqn >= R_MESH_NORMAL) && (eqn <= R_MESH_TANG2)) {
        BC_Types[ibc].desc->rotate = R_MESH1;
      }
      if ((eqn >= R_SOLID_NORMAL) && (eqn <= R_SOLID_TANG2)) {
        BC_Types[ibc].desc->rotate = R_SOLID1;
      }
      if ((eqn >= R_MOM_NORMAL) && (eqn <= R_MOM_TANG2)) {
        BC_Types[ibc].desc->rotate = R_MOMENTUM1;
      }

      break;

    case TABLE_WICV_BC:
    case TABLE_WICS_BC:
    case TABLE_BC:
      /*
       * Read in boundary data as a Data_Table structure
       */
      /*
       * Because there may be more than open TABLE_BC we must allocated
       * another BC_desc and fill it with the appropriate information.
       * We do this via the new_BC_Desc route so that it will integrate
       * into the parallel code.
       */

      new_BC_Desc = ((struct BC_descriptions **)realloc(
          new_BC_Desc, (num_new_BC_Desc + 1) * sizeof(struct BC_descriptions *)));

      /* Allocated  a new BC_desc for this Table BC and copy over the old info
       * that is listed in mm_names.h under TABLE_BC. Later this will be changed to reflect
       * the input deck card
       */

      new_BC_Desc[num_new_BC_Desc] = alloc_BC_description(BC_Types[ibc].desc);

      BC_Types[ibc].desc = new_BC_Desc[num_new_BC_Desc];

      BC_Types[ibc].index_dad = num_new_BC_Desc++; /* This is important to Phil */

      if (num_BC_Tables == MAX_BC_TABLES) {
        GOMA_EH(GOMA_ERROR, "Maximum TABLE_BCs exceeded .");
      }

      /*
       * Now read through the remainder of the BC card. Fill in the new BC_desc to
       * reflect was read ( equation, species_no etc.).  Also allocate space for
       * Data_Table struct and save a pointer to it in BC_Types[ibc]->table.
       * Finally, read and sort the tabular data itself
       */

      BC_Tables[num_BC_Tables] = setup_table_BC(ifp, input, &BC_Types[ibc], echo_string);

      /*
       * Record which table in BC_Tables is associate with this BC.  Important
       * when sending info to other processors
       */

      BC_Types[ibc].table_index = num_BC_Tables++;

      break;

      /*
       *
       * Read in Lagrange multiplier boundary conditions and add
       * associated augmenting condition
       *
       */
    case LGR_FLOWRATE_BC:
      /*
       * This boundary condition enforces a flowrate requirement over
       * the application sideset.  The flowrate requirement is a global
       * integrated constraint on the normal velocity.  It resembles very
       * much an augmenting condition.  This boundary condition is enforced
       * via a single constant lagrange multiplier unknown associated with
       * the sideset.  The value of the Lagrange multiplier is the additional
       * degree of freedom that appears just like any other augmenting
       * condition, but it isn't associated with a boundary condition parameter
       * or material parameter.  In addition, the lagrange multiplier multiplies
       * a surface integral on the sideset which is added to the residual of the
       * fluid momentum equation as an additional momentum source.  This is in
       * contrast to the regular augmenting condition algorithm which in general
       * doesn't change the residual vector.  When you work out the math you discover
       * that this additional momentum source term looks exactly like the
       * the surface integral that is used to apply a constant FLOW_PRESSURE
       * boundary condition with the multiplier taking the place of the pressure.
       * Actually, its value is the negative of the pressure.
       *
       * So what we have here is a hybrid sort of thing.  On the one hand,  it
       * applies a global integrated constraint via the augmenting condition
       * algorithm so we need to set this condition up here.  And also it applies
       * a more local additional the residual and jacobian matrices, which we
       * can apply via the regular bc_integ, weak integrated formality.  So pay
       * attention and here we go.....
       */

      /* First, read in the flowrate and a starting guess for the multiplier
       * as two float parameters on the BC = FLOWRATE line
       */

      if (fscanf(ifp, "%lf", &BC_Types[ibc].BC_Data_Float[0]) != 1) {

        sr = sprintf(err_msg, "Expected flt to start for %s on %s ID=%d.\n",
                     BC_Types[ibc].desc->name1, BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      if (fscanf(ifp, "%s", input) != 1) {
        sr = sprintf(err_msg, "Expected 2nd parameter for %s on %s ID=%d.\n",
                     BC_Types[ibc].desc->name1, BC_Types[ibc].Set_Type, BC_Types[ibc].BC_ID);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      SPF(endofstring(echo_string), " %.4g %s", BC_Types[ibc].BC_Data_Float[0], input);

      BC_Types[ibc].BC_Data_Float[1] = (float)strtod(input, NULL);

      /*
       * Next see if ACs have been defined above.  If so add space for one more, if not, allocate
       * an array of one
       */

      if (augc == NULL) {
        augc = alloc_struct_1(struct AC_Information, 1);
        nAC = 1;
      } else {
        augc = realloc_struct_1(augc, struct AC_Information, nAC + 1, nAC);
        nAC++;
      }

      /* Record the index of the attached AC in the BC_Types structure */

      BC_Types[ibc].BC_Data_Int[0] = nAC - 1;

      /*
       * Next up...fill in the pieces of the AC structure with the appropriate info.
       */
      augc[nAC - 1].Type = AC_LGRM;     /* Attached AC is a Lagrange multiplier constraint */
      augc[nAC - 1].MTID = -1;          /* this integrates over all blocks sharing sideset */
      augc[nAC - 1].MFID = VOLUME_FLUX; /* This flux type used in evaluate_flux */
      augc[nAC - 1].COMPID = 0.0;
      augc[nAC - 1].SSID = BC_Types[ibc].BC_ID;               /* Sideset set id */
      augc[nAC - 1].CONSTV = -BC_Types[ibc].BC_Data_Float[0]; /* Store the flow rate value */
      augc[nAC - 1].BCID = ibc; /* record which boundary condition AC is attached to */
      augc[nAC - 1].DFID = 0;   /* record index of flow rate parameter (for continuation */
      augc[nAC - 1].iread = 0;

      /*
       * Now check the BC card again to see if we are suppose to
       * read the Lagrange multiplier from the guess file
       */

      if (!strcmp(input, "read"))
        augc[nAC - 1].iread = 1;

      /*
       * Lastly,  DataFlt[0] gets the starting LM guess from BC_Data_Float[1]
       *
       */

      augc[nAC - 1].DataFlt = alloc_dbl_1(1, -BC_Types[ibc].BC_Data_Float[1]);
      augc[nAC - 1].len_AC = 1;

      /* Oh, one more thing..reset the nAC value in all the AC's */

      for (i = 0; i < nAC; i++)
        augc[i].nAC = nAC;

      break;

      /*
       * Fall through for all unimplemented cases
       */
    default:
      fprintf(stderr, "%s:\tBC type not yet implemented.\n", yo);
      exit(-1);
      break;
    }

    BC_consistency(&BC_Types[ibc]);

    ECHO(echo_string, echo_file);

    /* check hunting BC data floats for out-of-bounds*/
    if (nHC) {
      for (i = 0; i < nHC; i++) {
        if (hunt[i].Type == 1 && hunt[i].BCID == ibc) {
          if (hunt[i].DFID >= BC_Types[ibc].max_DFlt) {
            GOMA_WH(-1, "Whoa.... hunting data float outside range\n");
            fprintf(stderr, " HC %d BC %d of %d\n", i, hunt[i].DFID, BC_Types[ibc].max_DFlt);
          }
        }
      }
    }
    /* check loca BC data floats for out-of-bounds*/
    if (nCC) {
      for (i = 0; i < nCC; i++) {
        if (cpcc[i].Type == 1 && cpcc[i].BCID == ibc) {
          if (cpcc[i].DFID >= BC_Types[ibc].max_DFlt)
            GOMA_WH(-1, "Whoa.... loca data float outside range\n");
        }
      }
    }
  }

  ECHO("END OF BC", echo_file);

  /* Check for overlap AC with no relevant boundary conditions */
  if (Do_Overlap && !overlap_bc) {
    GOMA_EH(GOMA_ERROR, "Overlap AC requested with no embedded or contact BC's!");
  }

  /*
   * Check for a pressure Datum Condition
   */
  iread = look_for_optional(ifp, "PRESSURE DATUM", input, '=');

  /*
   * Initialize, even if we don't use them.
   */

  pressure_datum_element = -1;
  pressure_datum_value = -1.0e12;

  if (iread == 1) {
    PRESSURE_DATUM = 1;
    if (fscanf(ifp, "%d %lf ", &pressure_datum_element, &pressure_datum_value) != 2) {
      SPF(echo_string, "\t(%s)", "PRESSURE DATUM =  missing data");
    } else {
      SPF(echo_string, "%s = %d %.4g", "PRESSURE DATUM", pressure_datum_element,
          pressure_datum_value);
    }
    ECHO(echo_string, echo_file);
  } else {
    PRESSURE_DATUM = 0;
  }

  /*
   * Read directions for Rotating Equations in 3D
   */
  iread = look_for_optional(ifp, "Rotation Specifications", input, '=');

  /* count number of rotation or non-rotation specifications */
  if (iread == 1) {

    SPF(echo_string, "\n%s =\n", "Rotation Specifications");
    ECHO(echo_string, echo_file);

    Num_ROT = count_list(ifp, "ROT", input, '=', "END OF ROT");
  } else {
    Num_ROT = 0;
  }

  SPF(echo_string, "%s = %d", "Number of rotation conditions", Num_ROT);
  ECHO(echo_string, echo_file);

  /*
   *  Allocate space for the vector, ROT_Types.
   */
  if (Num_ROT != 0) {
    ROT_Types = (struct Rotation_Specs *)smalloc(Num_ROT * sizeof(struct Rotation_Specs));

    /*
     * When they are "-1", these indeces mean a corresponding NULL pointer
     * for the BC_Desc ptr. If otherwise, they point into the so-indexed
     * BC_Desc[] array of structures. All necessitated by need for parallel
     * processors to know this information.
     */

#ifndef CDIM
#define CDIM 3
#endif

    for (i = 0; i < Num_ROT; i++) {
      for (j = 0; j < CDIM; j++) {
        ROT_Types[i].BC_desc_index[j] = -1;
      }
    }

    for (irc = 0; irc < Num_ROT; irc++) {

      /*
       * Handy quick pointer for *this* ROT_Type.
       */

      rot = ROT_Types + irc;

      strcpy(rot_eq_string, "?");

      look_for(ifp, "ROT", input, '=');

      /* Read equation type: MESH or MOM */

      if (fscanf(ifp, "%s", ts) != 1) {
        GOMA_EH(-1, "error reading Rotation equation name");
      }

      if (!strcmp(ts, "MESH")) {
        rot->eq_type = R_MESH1;
        strcpy(rot_eq_string, "d");
      } else if (!strcmp(ts, "MOM")) {
        rot->eq_type = R_MOMENTUM1;
        strcpy(rot_eq_string, "v");
      } else {
        sr = SPF(err_msg, "Invalid ROT[%d] eq_type \"%s\". Valid are: MOM, MESH.", irc, ts);
        GOMA_EH(-1, err_msg);
      }

      SPF(echo_string, "%s = %s", "ROT", ts);

      /* Read topology type */

      if (fscanf(ifp, "%s", ts) != 1) {
        GOMA_EH(-1, "error reading Rotation topology");
      }

      if (!strcmp(ts, "SURFACE")) {
        rot->type = FACE;

        /* Read list of SS needed to define topology */

        if (fscanf(ifp, "%d", &rot->ss_id[0]) != 1) {
          sr = sprintf(err_msg, "Error reading 1 int for ROT SURFACE ssid1 ...");
          GOMA_EH(GOMA_ERROR, err_msg);
        }

        sr = SPF(topo_string, "surf (%d)", rot->ss_id[0]);
        rot->ss_id[1] = 0;
        rot->ss_id[2] = 0;

        SPF(endofstring(echo_string), " %s %d", ts, rot->ss_id[0]);
      } else if (!strcmp(ts, "EDGE")) {
        rot->type = CURVE;
        /* Read list of SS needed to define topology */
        if (fscanf(ifp, "%d %d", &rot->ss_id[0], &rot->ss_id[1]) != 2) {
          sr = sprintf(err_msg, "Error reading 2 ints for ROT EDGE ssid1 ssid2 ...");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
        sr = sprintf(topo_string, "edge (%d & %d)", rot->ss_id[0], rot->ss_id[1]);
        rot->ss_id[2] = 0;

        SPF(endofstring(echo_string), " %s %d %d", ts, rot->ss_id[0], rot->ss_id[1]);
      } else if (!strcmp(ts, "VERTEX")) {
        rot->type = VERTEX;
        /* Read list of SS needed to define topology */
        if (fscanf(ifp, "%d %d %d", &rot->ss_id[0], &rot->ss_id[1], &rot->ss_id[2]) != 3) {
          sprintf(err_msg, "Error reading 3 ints for ROT VERTEX ssid1 ssid2 ssid3 ...");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
        sprintf(topo_string, "vrtx (%d & %d & %d)", rot->ss_id[0], rot->ss_id[1], rot->ss_id[2]);

        SPF(endofstring(echo_string), " %s %d %d %d", ts, rot->ss_id[0], rot->ss_id[1],
            rot->ss_id[2]);

      } else if (!strcmp(ts, "VOLUME")) {
        rot->type = BODY;
        rot->ss_id[0] = 0;
        rot->ss_id[1] = 0;
        rot->ss_id[2] = 0;
        sr = sprintf(topo_string, "volume");

        SPF(endofstring(echo_string), " %s %s", ts, "volume");
      } else {
        SPF(err_msg, "Invalid ROT[%d] topology \"%s\". Valid are: SURFACE, EDGE, VERTEX, VOLUME.",
            irc, ts);
        GOMA_EH(-1, err_msg);
      }

      rot->ROTATE = 0;
      rot->elems = NULL;
      /* loop over equation components and read BC list w/ SS#'s */
      for (p = 0; p < DIM; p++) { /* p */
        if (fscanf(ifp, "%s %d", ts, &rot->BC_SS[p]) != 2) {
          sprintf(err_msg, "Expected 1 string, 1 int for rotation instruction %d for ROT %s", p + 1,
                  topo_string);
          GOMA_EH(-1, err_msg);
        }

        strcpy(condition_string[p], ts);

        /* compare string to valid rotation vectors */
        if (!strcmp(ts, "NONE") || !strcmp(ts, "NA") || !strcmp(ts, "NO")) {
          rot->BC_Type[p] = ROT_NONE;
          rot->BC_desc[p] = NULL;
          rot->BC_SS[p] = -1;
          /*
           * Case: a = a	 self replacement ain't even worth the effort
           */
          strcpy(instruction_string[p], coordinate_string[p]);
        } else if (!strcmp(ts, "N")) {
          rot->BC_Type[p] = ROT_N;
          rot->BC_desc[p] = NULL;
          rot->ROTATE = 1;
          rot->BC_SS[p] = -1;
          strcpy(instruction_string[p], "n");
        } else if (!strcmp(ts, "N2")) {
          rot->BC_Type[p] = ROT_N2;
          rot->BC_desc[p] = NULL;
          rot->ROTATE = 1;
          rot->BC_SS[p] = -1;
          strcpy(instruction_string[p], "n2");
        } else if (!strcmp(ts, "N3")) {
          rot->BC_Type[p] = ROT_N3;
          rot->BC_desc[p] = NULL;
          rot->ROTATE = 1;
          rot->BC_SS[p] = -1;
          strcpy(instruction_string[p], "n3");
        } else if (!strcmp(ts, "T")) {
          rot->BC_Type[p] = ROT_T;
          rot->BC_desc[p] = NULL;
          rot->ROTATE = 1;
          rot->BC_SS[p] = -1;
          strcpy(instruction_string[p], "t");
        } else if (!strcmp(ts, "T1")) {
          rot->BC_Type[p] = ROT_T1;
          rot->BC_desc[p] = NULL;
          rot->ROTATE = 1;
          rot->BC_SS[p] = -1;
          strcpy(instruction_string[p], "t1");
        } else if (!strcmp(ts, "T2")) {
          rot->BC_Type[p] = ROT_T2;
          rot->BC_desc[p] = NULL;
          rot->ROTATE = 1;
          rot->BC_SS[p] = -1;
          strcpy(instruction_string[p], "t2");
        } else if (!strcmp(ts, "B")) {
          rot->BC_Type[p] = ROT_B;
          rot->BC_desc[p] = NULL;
          rot->ROTATE = 1;
          rot->BC_SS[p] = -1;
          strcpy(instruction_string[p], "b");
        } else if (!strcmp(ts, "S")) {
          rot->BC_Type[p] = ROT_S;
          rot->BC_desc[p] = NULL;
          rot->ROTATE = 1;
          rot->BC_SS[p] = -1;
          strcpy(instruction_string[p], "s");
        } else if (!strcmp(ts, "X")) {
          rot->BC_Type[p] = ROT_X;
          rot->BC_desc[p] = NULL;
          rot->ROTATE = 1;
          rot->BC_SS[p] = -1;
          strcpy(instruction_string[p], "ex");
        } else if (!strcmp(ts, "Y")) {
          rot->BC_Type[p] = ROT_Y;
          rot->BC_desc[p] = NULL;
          rot->ROTATE = 1;
          rot->BC_SS[p] = -1;
          strcpy(instruction_string[p], "ey");
        } else if (!strcmp(ts, "Z")) {
          rot->BC_Type[p] = ROT_Z;
          rot->BC_desc[p] = NULL;
          rot->ROTATE = 1;
          rot->BC_SS[p] = -1;
          strcpy(instruction_string[p], "ez");
        } else if (!strcmp(ts, "GD")) {
          rot->BC_Type[p] = ROT_GD;
          rot->BC_desc[p] = NULL;
          rot->BC_SS[p] = -1;
          strcpy(instruction_string[p], "GD");
        } else {
          for (k = 0; k < Num_BC_Names; k++) {
            if (!strcmp(ts, BC_Desc[k].name1) || !strcmp(ts, BC_Desc[k].name2)) {
              rot->BC_Type[p] = BC_Desc[k].BC_Name;
              rot->BC_desc[p] = &BC_Desc[k];
              rot->BC_desc_index[p] = k;
            }
          }

          if (detect_BC(rot->BC_Type[p], rot->BC_SS[p]) == 0) {
            SPF(err_msg, "ROT card: %d  \t >->->BC : %s %d not found<-<-<\n", irc, ts,
                rot->BC_SS[p]);
            GOMA_EH(GOMA_ERROR, err_msg);
          }
        }

        SPF(endofstring(echo_string), " %s %d", ts, (rot->BC_SS[p] == -1 ? 0 : rot->BC_SS[p]));
      }

      /* Read method of tangent calc and list of floats */
      strcpy(seed_string, "unseeded");
      if (fscanf(ifp, "%s", ts) != 1) {
        rot->method = ROT_BASIS;
        GOMA_WH(-1, "No Rotation method");
      }

      if (!strcmp(ts, "SEED")) {
        rot->method = ROT_SEED;
        /* Read list of SS needed to define topology */
        if (fscanf(ifp, "%lf %lf %lf", &sx, &sy, &sz) != 3)
          GOMA_EH(GOMA_ERROR, "Error reading SEED values for Rotation");
        rot->seed[0] = sx;
        rot->seed[1] = sy;
        rot->seed[2] = sz;
        sr = SPF(seed_string, "SEED [%g,%g,%g]", sx, sy, sz);

      } else if (!strcmp(ts, "NONE")) {
        rot->method = ROT_BASIS;
        strcpy(seed_string, "NONE");
      } else if (!strcmp(ts, "BASIS")) {
        rot->method = ROT_BASIS;
        strcpy(seed_string, "BASIS");
      } else if (!strcmp(ts, "BASIS_RESEED")) {
        rot->method = ROT_BASIS_RESEED;
        strcpy(seed_string, "BASIS_RESEED");
      } else if (!strcmp(ts, "BASIS_ONCE")) {
        rot->method = ROT_BASIS_ONCE;
        strcpy(seed_string, "BASIS_ONCE");
      } else {
        sr = SPF(err_msg,
                 "Invalid ROT[%d] seed_method \"%s\". Valid are: NONE, SEED, BASIS, BASIS_RESEED, "
                 "BASIS_ONCE.",
                 irc, ts);
        GOMA_EH(-1, err_msg);
      }

      SPF(endofstring(echo_string), " %s", seed_string);

      /*
       * A feeble attempt to wrought output to aid in
       * understanding what the ?!@# is going with this ROT.
       */

      /*
       * The rotation component instruction condition will be either
       * SOMEBCNAME ssid or it will be ssid=0 and some "x" or "n", for example.
       */
      for (p = 0; p < DIM; p++) {
        if (rot->BC_SS[p] == -1) {
          /*
           * Write a nice projected component like (n.R_d), (t1.R_d),
           * (t2.R_d), (b.R_d), (ex.R_d), (ey.R_d), (ez.R_d) or
           * (1.d) for "unchanged".
           */

          sr = sprintf(condition_string[p], "(%s.R_%s)", instruction_string[p], rot_eq_string);
        } else {
          /*
           * Write a nice BC constraint for this component.
           */

          sr = sprintf(condition_string[p], "(%s @ %d)", rot->BC_desc[p]->name1, rot->BC_SS[p]);
        }

        sr = snprintf(component_string[p], 119, "R_%s%s=%s", rot_eq_string, coordinate_string[p],
                      condition_string[p]);
      }

      SPF(endofstring(echo_string), "\n\t (%3d. %s %s %s %s %s)", irc + 1, topo_string,
          component_string[0], component_string[1], component_string[2], seed_string);

      ECHO(echo_string, echo_file);
    } /* end of loop over irc, rotation specs. */
    ECHO("END OF ROT", echo_file);
  } /* end of if Num_Rot != 0 */
}
/* rd_bc_specs -- read input file for boundary condition specifications */

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/
/*
 * Check the consistency of a Boundary Condition structure.
 * Initially this routine was implemented to check and make sure that
 * NS designators are not used on BC's that are applicable only to
 * sidesets, and vice-versa.  However, it is also a good place to put
 * code that checks other BC input parameters, i.e. angles are always
 * in the range (0,2Pi).
 *
 * Author: T.A. Baer
 * Date  : December 1999
 */

static int BC_consistency(struct Boundary_Condition *BC_Type) {
  char err_msg[MAX_CHAR_IN_INPUT];
  static const char yo[] = "rd_bc_specs";

  /* Check that BC is appropriate to SS */
  if (!strcmp(BC_Type->Set_Type, "SS")) {
    switch (BC_Type->desc->method) {
    case COLLOCATE_SURF:
    case WEAK_INT_SURF:
    case STRONG_INT_SURF:
    case WEAK_SHELL_GRAD:
    case STRONG_SHELL_GRAD:
    case STRONG_INT_EDGE:
    case WEAK_INT_EDGE:
    case COLLOCATE_EDGE:
    case CONTACT_SURF:
    case EMBEDDED_SURF:
    case WEAK_SHARP_INT:
    case WEAK_INT_NEDELEC:
    case STRONG_INT_NEDELEC:
      break;
    case LS_SPECIAL:
      if (BC_Type->desc->BC_Name != LS_ADC_OLD_BC) {
        sprintf(err_msg, "%s %s %d\n\t\t %s", "BC Consistency error detected. ",
                BC_Type->desc->name1, BC_Type->BC_ID, " BC is not applicable to side sets ");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      break;
    default:
      sprintf(err_msg, "%s %s %d\n\t\t %s", "BC Consistency error detected. ", BC_Type->desc->name1,
              BC_Type->BC_ID, " BC is not applicable to side sets ");
      GOMA_EH(GOMA_ERROR, err_msg);
      break;
    }
  }

  if (!strcmp(BC_Type->Set_Type, "NS")) {
    switch (BC_Type->desc->method) {
    case DIRICHLET:
    case SPECIAL:
      break;
    case LS_SPECIAL:

      if (BC_Type->desc->BC_Name != LS_ADC_BC) {
        sprintf(err_msg, "%s %s %d\n\t\t %s", "BC Consistency error detected. ",
                BC_Type->desc->name1, BC_Type->BC_ID, " BC is not applicable to node sets ");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      break;
    default:
      // This boundary condition uses a side on a shell element in 2D. This is
      // a node set consisting of one node. Therefore, it's ok.
      switch (BC_Type->desc->BC_Name) {
      case SH_GAMMA1_DERIV_SYMM_BC:
      case SH_GAMMA2_DERIV_SYMM_BC:
      case SHELL_TFMP_FREE_LIQ_BC:
      case SHELL_TFMP_FREE_GAS_BC:
      case SHELL_TFMP_GRAD_S_BC:
      case GRAD_LUB_PRESS_BC:
      case SH_SDET_BC:
      case SH_MESH2_WEAK_BC:
      case SHELL_LUBRICATION_OUTFLOW_BC:
        break;
      default:
        sprintf(err_msg, "%s %s %d\n\t\t %s", "BC Consistency error detected. ",
                BC_Type->desc->name1, BC_Type->BC_ID, " BC is not applicable to node sets ");
        GOMA_EH(GOMA_ERROR, err_msg);
        break;
      }
    }
  }

  if (!strcmp(BC_Type->desc->name1, "VELO_THETA_TPL") ||
      !strcmp(BC_Type->desc->name1, "VELO_THETA_HOFFMAN") ||
      !strcmp(BC_Type->desc->name1, "VELO_THETA_COX") ||
      !strcmp(BC_Type->desc->name1, "VELO_THETA_SHIK") ||
      !strcmp(BC_Type->desc->name1, "WETTING_SPEED_BLAKE") ||
      !strcmp(BC_Type->desc->name1, "WETTING_SPEED_HOFFMAN") ||
      !strcmp(BC_Type->desc->name1, "WETTING_SPEED_SHIK") ||
      !strcmp(BC_Type->desc->name1, "WETTING_SPEED_COX")) { /* 1 */
    /*
     * Check some parameters for sanity.
     *
     */

    /*
     * Equilibrium contact angles shall be between 1 and 179 degrees. If you can't
     * fit within that range you have no business here.
     */

    if (BC_Type->BC_Data_Float[0] < 1. || BC_Type->BC_Data_Float[0] > 179.) {
      sr = sprintf(err_msg, "BC %s angle %g not in (1,179)", BC_Type->desc->name1,
                   BC_Type->BC_Data_Float[0]);
      GOMA_EH(GOMA_ERROR, err_msg);
    }

    /*
     * Normal vector magnitudes ought to be pretty close to 1.
     */

    if (fabs(BC_Type->BC_Data_Float[1] * BC_Type->BC_Data_Float[1] +
             BC_Type->BC_Data_Float[2] * BC_Type->BC_Data_Float[2] +
             BC_Type->BC_Data_Float[3] * BC_Type->BC_Data_Float[3] - 1) > 1e-3) {
      sr = sprintf(err_msg, "BC %s ss normal (%g,%g,%g) not unit.", BC_Type->desc->name1,
                   BC_Type->BC_Data_Float[1], BC_Type->BC_Data_Float[2], BC_Type->BC_Data_Float[3]);
      GOMA_WH(-1, "Will use variable wall normal for wetting bc");
    }

    /*
     * Other parameters in our Blake model should be positive...
     */

    if (BC_Type->BC_Data_Float[4] < 0) {
      sr = sprintf(err_msg, "BC %s preexponential %g negative.", BC_Type->desc->name1,
                   BC_Type->BC_Data_Float[4]);
      GOMA_EH(GOMA_ERROR, err_msg);
    }

    if (BC_Type->BC_Data_Float[5] < 0) {
      sr = sprintf(err_msg, "BC %s thermally-scaled surface energy %g negative.",
                   BC_Type->desc->name1, BC_Type->BC_Data_Float[5]);
      GOMA_EH(GOMA_ERROR, err_msg);
    }

    /*
     * t_relax = 0 is OK; means just do it with no relaxation...
     */

    if (BC_Type->BC_Data_Float[6] < 0) {
      sr = sprintf(err_msg, "BC %s relaxation time %g negative.", BC_Type->desc->name1,
                   BC_Type->BC_Data_Float[6]);
      GOMA_EH(GOMA_ERROR, err_msg);
    }

    if (BC_Type->BC_Data_Float[7] < 0) {
      sr = sprintf(err_msg, "BC %s old tpl velocity %g negative.", BC_Type->desc->name1,
                   BC_Type->BC_Data_Float[7]);
      GOMA_EH(GOMA_ERROR, err_msg);
    }
  }

  if (!strcmp(BC_Type->desc->name1, "VAR_CA_EDGE")) {
    /*
     * Normalize substrate normal vector
     */
    dbl *norm = &BC_Type->BC_Data_Float[5];
    dbl mag = sqrt(SQUARE(norm[0]) + SQUARE(norm[1]) + SQUARE(norm[2]));
    if (mag == 0.0) {
      sprintf(err_msg, "%s: ERROR for BC %s, zero length normal vector\n", yo,
              BC_Type->desc->name1);
      GOMA_EH(GOMA_ERROR, err_msg);
    } else {
      norm[0] /= mag;
      norm[1] /= mag;
      norm[2] /= mag;
    }
  }

  return 0;
}

/* Search through the BC list to find a match for BC_name and BC_SS
 *
 * Return TRUE if match of both found on same BC
 * Return FALSE if  no such match found
 *
 * Used when reading ROT cards to assure BCs listed are present
 *
 *
 * Author : T A Baer
 * Date   : Jan 10, 2001
 */

static int detect_BC(int BC_name, int BC_ID)

{

  int k = 0, found = FALSE;

  while (k < Num_BC && !found) {
    found = (BC_Types[k].BC_Name == BC_name) && (BC_Types[k].BC_ID == BC_ID);

    k++;
  }

  return (found);
}

/* END of file mm_input_bc.c */
