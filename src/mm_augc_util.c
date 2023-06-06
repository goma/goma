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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ac_update_parameter.h"
#include "az_aztec.h"
#include "bc_contact.h"
#include "dp_types.h"
#include "dpi.h"
#include "el_elm.h"
#include "el_geom.h"
#include "exo_struct.h"
#include "loca_const.h"
#include "mm_as.h"
#include "mm_as_structs.h"
#include "mm_augc_util.h"
#include "mm_eh.h"
#include "mm_fill.h"
#include "mm_fill_ls.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_terms.h"
#include "mm_fill_util.h"
#include "mm_flux.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "mm_post_def.h"
#include "mm_unknown_map.h"
#include "mpi.h"
#include "rd_mesh.h"
#include "rf_allo.h"
#include "rf_bc.h"
#include "rf_bc_const.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_io.h"
#include "rf_io_const.h"
#include "rf_mp.h"
#include "rf_node_const.h"
#include "rf_solve.h"
#include "rf_util.h"
#include "rf_vars_const.h"
#include "sl_auxutil.h"
#include "sl_util_structs.h"
#include "std.h"
#include "user_ac.h"
#include "util/aprepro_helper.h"

#define GOMA_MM_AUGC_UTIL_C

/*
 * Prototype declarations of static functions.
 */

static int estimate_bAC(int, double[], double **, int, Comm_Ex *, MF_Args *);

#ifdef NOT_USED
static int estimate_dAC_LSvel(int, double[], double **, int, Comm_Ex *, MF_Args *);
#endif

static int estimate_dAC_ALC(int, double[], double **, int, Comm_Ex *, MF_Args *);

static double
getPositionAC(struct AC_Information *augc, double *cAC_iAC, double *soln, Exo_DB *exo);

static double getAngleAC(struct AC_Information *augc, double *cAC_iAC, double *soln, Exo_DB *exo);

/*

   AUGMENTING CONDITION UTILITY ROUTINES

   BY IAN GATES

   2/98 - 9/98



*/
//! Get the value of the unknown corresponding to the augmented condition
//! and load it into the vector of Augmented condition unknowns
/*!
 *
 */
void load_extra_unknownsAC(int iAC,     /* ID NUMBER OF AC'S */
                           double *xa,  /* VECTOR OF EXTRA UNKNOWNS */
                           Comm_Ex *cx, /* array of communications structures */
                           Exo_DB *exo, /* ptr to the finite element mesh database */
                           Dpi *dpi)    /* distributed processing information */
{

  int mn;
  int ibc, idf;
  struct Boundary_Condition *BC_Type;

  /*
   * 		BEGIN EXECUTION
   */

  if (iAC < 0)
    return;

  if (augc[iAC].Type == AC_LGRM) {
    xa[iAC] = augc[iAC].DataFlt[0];
    return;
  }

  if (augc[iAC].Type == AC_OVERLAP || augc[iAC].Type == AC_PF_CONSTRAINT) {
    xa[iAC] = augc[iAC].lm_value;
  }

  if (augc[iAC].Type == AC_PERIODIC) {
    xa[iAC] = augc[iAC].lm_value;
  }

  if (augc[iAC].Type == AC_USERBC || augc[iAC].Type == AC_VOLUME || augc[iAC].Type == AC_FLUX ||
      augc[iAC].Type == AC_LS_VEL || augc[iAC].Type == AC_POSITION || augc[iAC].Type == AC_ANGLE) {

    ibc = augc[iAC].BCID;
    idf = augc[iAC].DFID;

    /*  case for user defined parameter  */

    if (ibc == APREPRO_AC_BCID || ibc == APREPRO_LIB_AC_BCID) {
      xa[iAC] = augc[iAC].DataFlt[0];
    } else {

      switch (BC_Types[ibc].BC_Name) {
      case SPLINE_BC:
      case SPLINEX_BC:
      case SPLINEY_BC:
      case SPLINEZ_BC:
      case SPLINE_RS_BC:
      case SPLINEX_RS_BC:
      case SPLINEY_RS_BC:
      case SPLINEZ_RS_BC:
      case FILLET_BC:
      case DOUBLE_RAD_BC:
      case FEATURE_ROLLON_BC:
      case ROLL_FLUID_BC:
      case UVARY_BC:
      case VVARY_BC:
      case WVARY_BC:
      case U_PARABOLA_BC:
      case V_PARABOLA_BC:
      case W_PARABOLA_BC:
      case PRESSURE_USER_BC:
      case FLOW_PRESS_USER_BC:
      case T_USER_BC:
      case UUSER_BC:
      case VUSER_BC:
      case WUSER_BC:
      case QUSER_BC:
      case DX_USER_BC:
      case DY_USER_BC:
      case DZ_USER_BC:
      case P_LIQ_USER_BC:
      case SH_P_OPEN_USER_BC:
      case VAR_CA_USER_BC:
      case CA_EDGE_OR_FIX_BC:
      case YFLUX_USER_BC:
      case SURFACE_CHARGE_BC:
      case YUSER_BC:
      case FORCE_USER_BC:
      case FORCE_USER_RS_BC:
        xa[iAC] = BC_Types[ibc].u_BC[idf];
        break;
      case TABLE_BC:
      case TABLE_WICV_BC:
      case TABLE_WICS_BC:
        BC_Type = &BC_Types[ibc];
        xa[iAC] = BC_Type->table->f[idf];
        break;
      default:
        xa[iAC] = BC_Types[ibc].BC_Data_Float[idf];
        break;
      }
    }

  } else if (augc[iAC].Type == AC_USERMAT || augc[iAC].Type == AC_FLUX_MAT) {

    mn = map_mat_index(augc[iAC].MTID);

    switch (augc[iAC].MPID) {

      /*
       * General Model Constants
       */

    case TAGC_THERMAL_CONDUCTIVITY:
      xa[iAC] = mp_glob[mn]->thermal_conductivity;
      break;

    case TAGC_ACOUSTIC_WAVENUMBER:
      xa[iAC] = mp_glob[mn]->wave_number;
      break;

    case TAGC_ACOUSTIC_IMPEDANCE:
      xa[iAC] = mp_glob[mn]->acoustic_impedance;
      break;

    case TAGC_ACOUSTIC_ABSORPTION:
      xa[iAC] = mp_glob[mn]->acoustic_absorption;
      break;

    case TAGC_REFRACTIVE_INDEX:
      xa[iAC] = mp_glob[mn]->refractive_index;
      break;

    case TAGC_LIGHT_ABSORPTION:
      xa[iAC] = mp_glob[mn]->light_absorption;
      break;

    case TAGC_EXTINCTION_INDEX:
      xa[iAC] = mp_glob[mn]->extinction_index;
      break;

    case TAGC_ELECTRICAL_CONDUCTIVITY:
      xa[iAC] = mp_glob[mn]->electrical_conductivity;
      break;

    case TAGC_PERMITTIVITY:
      xa[iAC] = mp_glob[mn]->permittivity;
      break;

    case TAGC_VISCOSITY:
      xa[iAC] = mp_glob[mn]->viscosity;
      break;

    case TAGC_SURFACE_TENSION:
      xa[iAC] = mp_glob[mn]->surface_tension;
      break;

    case TAGC_HEAT_CAPACITY:
      xa[iAC] = mp_glob[mn]->heat_capacity;
      break;

    case TAGC_VOLUME_EXPANSION:
      xa[iAC] = mp_glob[mn]->Volume_Expansion;
      break;

    case TAGC_DENSITY:
      xa[iAC] = mp_glob[mn]->density;
      break;

    case TAGC_POROSITY:
      xa[iAC] = mp_glob[mn]->porosity;
      break;

    case TAGC_POROUS_COMPRESSIBILITY:
      xa[iAC] = mp_glob[mn]->porous_compressibility;
      break;

    case TAGC_PERMEABILITY:
      xa[iAC] = mp_glob[mn]->permeability;
      break;

    case TAGC_REL_GAS_PERM:
      xa[iAC] = mp_glob[mn]->rel_gas_perm;
      break;

    case TAGC_REL_LIQ_PERM:
      xa[iAC] = mp_glob[mn]->rel_liq_perm;
      break;

    case TAGC_SATURATION:
      xa[iAC] = mp_glob[mn]->saturation;
      break;

    case TAGC_MELTING_POINT_LIQUIDUS:
      xa[iAC] = mp_glob[mn]->melting_point_liquidus;
      break;

    case TAGC_MELTING_POINT_SOLIDUS:
      xa[iAC] = mp_glob[mn]->melting_point_solidus;
      break;

    case TAGC_FLOWINGLIQUID_VISCOSITY:
      xa[iAC] = mp_glob[mn]->FlowingLiquid_viscosity;
      break;

    case TAGC_DIFFUSIVITY_0:
      xa[iAC] = mp_glob[mn]->diffusivity[0];
      break;

    case TAGC_DIFFUSIVITY_1:
      xa[iAC] = mp_glob[mn]->diffusivity[1];
      break;

      /*
       * Generalized Newtonian Models:
       * Newtonian, Power Law, Carreau or Bingham(1,2,3)
       */

    case TAGC_MU0:
      xa[iAC] = gn_glob[mn]->mu0;
      break;
    case TAGC_NEXP:
      xa[iAC] = gn_glob[mn]->nexp;
      break;
    case TAGC_MUINF:
      xa[iAC] = gn_glob[mn]->muinf;
      break;
    case TAGC_LAM:
      xa[iAC] = gn_glob[mn]->lam;
      break;
    case TAGC_AEXP:
      xa[iAC] = gn_glob[mn]->aexp;
      break;
    case TAGC_ATEXP:
      xa[iAC] = gn_glob[mn]->atexp;
      break;

      /*  CARREAU_WLF    */

    case TAGC_WLFC2:
      xa[iAC] = gn_glob[mn]->wlfc2;
      break;

    case TAGC_REFTEMP:
      xa[iAC] = mp_glob[mn]->reference[TEMPERATURE];
      break;

      /* these are for the BINGHAM yielding material model */

    case TAGC_TAU_Y:
      xa[iAC] = gn_glob[mn]->tau_y;
      break;
    case TAGC_FEXP:
      xa[iAC] = gn_glob[mn]->fexp;
      break;

      /* these are for SUSPENSION/FILLED_EPOXY models */

    case TAGC_MAXPACK:
      xa[iAC] = gn_glob[mn]->maxpack;
      break;
    case TAGC_FICKDIFF_X:
      xa[iAC] = mp_glob[mn]->u_fdiffusivity[0][0];
      break;
    case TAGC_FICKDIFF_Y:
      xa[iAC] = mp_glob[mn]->u_fdiffusivity[0][1];
      break;

      /* these can be implemented for Ryan's Qtensor model as needed */

#if 0
      case TAGC_QTENSOR_EXTENSION_P:
	xa[iAC] = mp_glob[mn]->Qtensor_Extension_P;
	break;
      case TAGC_QTENSOR_NCT:
	xa[iAC] = mp_glob[mn]->Qtensor_Nct;
	break;
#endif

      /* these are for CURE/EPOXY/FILLED_EPOXY models */

    case TAGC_GELPOINT:
      xa[iAC] = gn_glob[mn]->gelpoint;
      break;
    case TAGC_CUREAEXP:
      xa[iAC] = gn_glob[mn]->cureaexp;
      break;
    case TAGC_CUREBEXP:
      xa[iAC] = gn_glob[mn]->curebexp;
      break;

      /*
       * Constants used in the Viscoelastic Constitutive Equations
       */

    case TAGC_TIME_CONST:
    case TAGC_TIME_CONST1:
    case TAGC_TIME_CONST2:
    case TAGC_TIME_CONST3:
    case TAGC_TIME_CONST4:
    case TAGC_TIME_CONST5:
    case TAGC_TIME_CONST6:
    case TAGC_TIME_CONST7:
      xa[iAC] = ve_glob[mn][augc[iAC].MPID - TAGC_TIME_CONST]->time_const;
      break;

    case TAGC_WT_FUNC:
      xa[iAC] = vn_glob[mn]->wt_func;
      break;

    case TAGC_ALPHA:
    case TAGC_ALPHA1:
    case TAGC_ALPHA2:
    case TAGC_ALPHA3:
    case TAGC_ALPHA4:
    case TAGC_ALPHA5:
    case TAGC_ALPHA6:
    case TAGC_ALPHA7:
      xa[iAC] = ve_glob[mn][augc[iAC].MPID - TAGC_ALPHA]->alpha;
      break;

    case TAGC_PTT_XI:
      xa[iAC] = ve_glob[mn][augc[iAC].MPID - TAGC_PTT_XI]->xi;
      break;

    case TAGC_PTT_EPS:
      xa[iAC] = ve_glob[mn][augc[iAC].MPID - TAGC_PTT_EPS]->eps;
      break;

    case TAGC_SHIFT_FUNC:
    case TAGC_SHIFT_FUNC1:
      xa[iAC] = vn_glob[mn]->shift[augc[iAC].MPID - TAGC_SHIFT_FUNC];
      break;

    case TAGC_POLYMER_YIELD_STRESS:
      xa[iAC] = ve_glob[mn][augc[iAC].MPID - TAGC_POLYMER_YIELD_STRESS]->gn->tau_y;
      break;

    case TAGC_POLYMER_YIELD_EXPONENT:
      xa[iAC] = ve_glob[mn][augc[iAC].MPID - TAGC_POLYMER_YIELD_EXPONENT]->gn->fexp;
      break;
      /*
       * Constants used in the Elasticity Constitutive Equations
       */

    case TAGC_LAME_MU:
      xa[iAC] = elc_glob[mn]->lame_mu;
      break;

    case TAGC_LAME_MU_CONTACT_LINE_G0:
      xa[iAC] = *(elc_glob[mn]->u_mu + 1);
      break;

    case TAGC_LAME_MU_CONTACT_LINE_G1:
      xa[iAC] = *(elc_glob[mn]->u_mu + 2);
      break;

    case TAGC_LAME_MU_CONTACT_LINE_R0:
      xa[iAC] = *(elc_glob[mn]->u_mu + 3);
      break;

    case TAGC_LAME_LAMBDA:
      xa[iAC] = elc_glob[mn]->lame_lambda;
      break;

    case TAGC_BEND_STIFFNESS:
      xa[iAC] = elc_glob[mn]->bend_stiffness;
      break;

    case TAGC_CONV_LAG_VELX:
      xa[iAC] = *(elc_glob[mn]->v_mesh_sfs);
      break;

    case TAGC_CONV_LAG_VELY:
      xa[iAC] = *(elc_glob[mn]->v_mesh_sfs + 1);
      break;

    case TAGC_CONV_LAG_VELZ:
      xa[iAC] = *(elc_glob[mn]->v_mesh_sfs + 2);
      break;

    case TAGC_CONV_LAG_ROTRATE:
      xa[iAC] = *(elc_glob[mn]->u_v_mesh_sfs);
      break;

    case TAGC_CONV_LAG_ROT_X0:
      xa[iAC] = *(elc_glob[mn]->u_v_mesh_sfs + 1);
      break;

    case TAGC_CONV_LAG_ROT_Y0:
      xa[iAC] = *(elc_glob[mn]->u_v_mesh_sfs + 2);
      break;

    case TAGC_CONV_LAG_ROT_Z0:
      xa[iAC] = *(elc_glob[mn]->u_v_mesh_sfs + 3);
      break;

    case TAGC_RS_LAME_MU:
      xa[iAC] = elc_rs_glob[mn]->lame_mu;
      break;

    case TAGC_RS_LAME_LAMBDA:
      xa[iAC] = elc_rs_glob[mn]->lame_lambda;
      break;

    case TAGC_RS_CONV_LAG_VELX:
      xa[iAC] = *(elc_rs_glob[mn]->v_mesh_sfs);
      break;

    case TAGC_RS_CONV_LAG_VELY:
      xa[iAC] = *(elc_rs_glob[mn]->v_mesh_sfs + 1);
      break;

    case TAGC_RS_CONV_LAG_VELZ:
      xa[iAC] = *(elc_rs_glob[mn]->v_mesh_sfs + 2);
      break;

    case TAGC_RS_CONV_LAG_ROTRATE:
      xa[iAC] = *(elc_rs_glob[mn]->u_v_mesh_sfs);
      break;

    case TAGC_RS_CONV_LAG_ROT_X0:
      xa[iAC] = *(elc_rs_glob[mn]->u_v_mesh_sfs + 1);
      break;

    case TAGC_RS_CONV_LAG_ROT_Y0:
      xa[iAC] = *(elc_rs_glob[mn]->u_v_mesh_sfs + 2);
      break;

    case TAGC_RS_CONV_LAG_ROT_Z0:
      xa[iAC] = *(elc_rs_glob[mn]->u_v_mesh_sfs + 3);
      break;

    case TAGC_POISSON:
      xa[iAC] = elc_glob[mn]->poisson;
      break;

    case TAGC_STRSS_FR_SOL_VOL_FRAC:
      xa[iAC] = elc_glob[mn]->Strss_fr_sol_vol_frac;
      break;

      /*
       * Constants used for Source Term Models
       */

    case TAGC_NSS_A0:
      xa[iAC] = mp_glob[mn]->momentum_source[0];
      break;

    case TAGC_NSS_A1:
      xa[iAC] = mp_glob[mn]->momentum_source[1];
      break;

    case TAGC_NSS_A2:
      xa[iAC] = mp_glob[mn]->momentum_source[2];
      break;

    case TAGC_NSS_A3:
      xa[iAC] = mp_glob[mn]->u_momentum_source[0];
      break;

    case TAGC_LUB_HGT_U0:
      xa[iAC] = mp_glob[mn]->u_heightU_function_constants[0];
      break;

    case TAGC_LUB_HGT_U1:
      xa[iAC] = mp_glob[mn]->u_heightU_function_constants[1];
      break;

    case TAGC_LUB_HGT_U2:
      xa[iAC] = mp_glob[mn]->u_heightU_function_constants[2];
      break;

    case TAGC_LUB_HGT_U3:
      xa[iAC] = mp_glob[mn]->u_heightU_function_constants[3];
      break;

    case TAGC_LUB_HGT_U4:
      xa[iAC] = mp_glob[mn]->u_heightU_function_constants[4];
      break;

    case TAGC_LUB_HGT_U5:
      xa[iAC] = mp_glob[mn]->u_heightU_function_constants[5];
      break;

    case TAGC_LUB_HGT_U6:
      xa[iAC] = mp_glob[mn]->u_heightU_function_constants[6];
      break;

    case TAGC_LUB_HGT_U7:
      xa[iAC] = mp_glob[mn]->u_heightU_function_constants[7];
      break;

    case TAGC_LUB_HGT_L0:
      xa[iAC] = mp_glob[mn]->u_heightL_function_constants[0];
      break;

    case TAGC_LUB_HGT_L1:
      xa[iAC] = mp_glob[mn]->u_heightL_function_constants[1];
      break;

    case TAGC_LUB_HGT_L2:
      xa[iAC] = mp_glob[mn]->u_heightL_function_constants[2];
      break;

    case TAGC_LUB_HGT_L3:
      xa[iAC] = mp_glob[mn]->u_heightL_function_constants[3];
      break;

    case TAGC_LUB_HGT_L4:
      xa[iAC] = mp_glob[mn]->u_heightL_function_constants[4];
      break;

    case TAGC_LUB_HGT_L5:
      xa[iAC] = mp_glob[mn]->u_heightL_function_constants[5];
      break;

    case TAGC_LUB_HGT_L6:
      xa[iAC] = mp_glob[mn]->u_heightL_function_constants[6];
      break;

    case TAGC_LUB_HGT_L7:
      xa[iAC] = mp_glob[mn]->u_heightL_function_constants[7];
      break;

    case TAGC_U_LUB_VELO_U0:
      xa[iAC] = mp_glob[mn]->u_veloU_function_constants[0];
      break;

    case TAGC_U_LUB_VELO_U1:
      xa[iAC] = mp_glob[mn]->u_veloU_function_constants[1];
      break;

    case TAGC_U_LUB_VELO_U2:
      xa[iAC] = mp_glob[mn]->u_veloU_function_constants[2];
      break;

    case TAGC_U_LUB_VELO_U3:
      xa[iAC] = mp_glob[mn]->u_veloU_function_constants[3];
      break;

    case TAGC_U_LUB_VELO_U4:
      xa[iAC] = mp_glob[mn]->u_veloU_function_constants[4];
      break;

    case TAGC_U_LUB_VELO_U5:
      xa[iAC] = mp_glob[mn]->u_veloU_function_constants[5];
      break;

    case TAGC_U_LUB_VELO_L0:
      xa[iAC] = mp_glob[mn]->u_veloL_function_constants[0];
      break;

    case TAGC_U_LUB_VELO_L1:
      xa[iAC] = mp_glob[mn]->u_veloL_function_constants[1];
      break;

    case TAGC_U_LUB_VELO_L2:
      xa[iAC] = mp_glob[mn]->u_veloL_function_constants[2];
      break;

    case TAGC_U_LUB_VELO_L3:
      xa[iAC] = mp_glob[mn]->u_veloL_function_constants[3];
      break;

    case TAGC_U_LUB_VELO_L4:
      xa[iAC] = mp_glob[mn]->u_veloL_function_constants[4];
      break;

    case TAGC_U_LUB_VELO_L5:
      xa[iAC] = mp_glob[mn]->u_veloL_function_constants[5];
      break;

    case TAGC_LUB_VELO_U0:
      xa[iAC] = mp_glob[mn]->veloU[0];
      break;

    case TAGC_LUB_VELO_U1:
      xa[iAC] = mp_glob[mn]->veloU[1];
      break;

    case TAGC_LUB_VELO_U2:
      xa[iAC] = mp_glob[mn]->veloU[2];
      break;

    case TAGC_LUB_VELO_L0:
      xa[iAC] = mp_glob[mn]->veloL[0];
      break;

    case TAGC_LUB_VELO_L1:
      xa[iAC] = mp_glob[mn]->veloL[1];
      break;

    case TAGC_LUB_VELO_L2:
      xa[iAC] = mp_glob[mn]->veloL[2];
      break;

    case TAGC_LUB_DCA_U0:
      xa[iAC] = mp_glob[mn]->u_dcaU_function_constants[0];
      break;

    case TAGC_LUB_DCA_U1:
      xa[iAC] = mp_glob[mn]->u_dcaU_function_constants[1];
      break;

    case TAGC_LUB_DCA_U2:
      xa[iAC] = mp_glob[mn]->u_dcaU_function_constants[2];
      break;

    case TAGC_LUB_DCA_U3:
      xa[iAC] = mp_glob[mn]->u_dcaU_function_constants[3];
      break;

    case TAGC_LUB_DCA_L0:
      xa[iAC] = mp_glob[mn]->u_dcaL_function_constants[0];
      break;

    case TAGC_LUB_DCA_L1:
      xa[iAC] = mp_glob[mn]->u_dcaL_function_constants[1];
      break;

    case TAGC_LUB_DCA_L2:
      xa[iAC] = mp_glob[mn]->u_dcaL_function_constants[2];
      break;

    case TAGC_LUB_DCA_L3:
      xa[iAC] = mp_glob[mn]->u_dcaL_function_constants[3];
      break;

    case TAGC_LUB_SOURCE_0:
      xa[iAC] = mp_glob[mn]->u_lubsource_function_constants[0];
      break;

    case TAGC_LUB_SOURCE_1:
      xa[iAC] = mp_glob[mn]->u_lubsource_function_constants[1];
      break;

    case TAGC_LUB_SOURCE_2:
      xa[iAC] = mp_glob[mn]->u_lubsource_function_constants[2];
      break;

    case TAGC_RST_FUNC_0:
      xa[iAC] = mp_glob[mn]->Rst_func;
      break;

    case TAGC_RST_FUNC_1:
      xa[iAC] = mp_glob[mn]->Rst_diffusion;
      break;

    case TAGC_RST_FUNC_2:
      xa[iAC] = mp_glob[mn]->Rst_func_supg;
      break;

    case TAGC_HEAT_SOURCE_0:
      xa[iAC] = mp_glob[mn]->u_heat_source[0];
      break;

    case TAGC_SPECIES_SOURCE_0_P0:
      xa[iAC] = mp_glob[mn]->u_species_source[0][0];
      break;

    case TAGC_SPECIES_SOURCE_0_P1:
      xa[iAC] = mp_glob[mn]->u_species_source[0][1];
      break;

    case TAGC_SPECIES_SOURCE_0_P2:
      xa[iAC] = mp_glob[mn]->u_species_source[0][2];
      break;

    case TAGC_SPECIES_SOURCE_0_P3:
      xa[iAC] = mp_glob[mn]->u_species_source[0][3];
      break;

    case TAGC_SPECIES_SOURCE_1_P0:
      xa[iAC] = mp_glob[mn]->u_species_source[1][0];
      break;

    case TAGC_SPECIES_SOURCE_1_P1:
      xa[iAC] = mp_glob[mn]->u_species_source[1][1];
      break;

    case TAGC_SPECIES_SOURCE_1_P2:
      xa[iAC] = mp_glob[mn]->u_species_source[1][2];
      break;

    case TAGC_SPECIES_SOURCE_1_P3:
      xa[iAC] = mp_glob[mn]->u_species_source[1][3];
      break;

    case TAGC_LATENT_HEAT_0:
      xa[iAC] = mp_glob[mn]->latent_heat_vap[0];
      break;

    case TAGC_LATENT_HEAT_1:
      xa[iAC] = mp_glob[mn]->latent_heat_vap[1];
      break;

    default:
      printf("\n\t Error: Invalid Material Property Tag %d\n", augc[iAC].MPID);
      exit(0);
      break;
    }
  }
  augc[iAC].tmp1 = xa[iAC];
} /* END of routine load_extra_unknownsAC  */
/*****************************************************************************/

/*

UPDATE PARAMETER FOR AUGMENTING CONDITIONS

BY IAN GATES

2/98 - 9/98

In this routine you take the value of the augmenting condition, from the aug solution
vector, and you put it back into the structures used by Goma to generate the problem.
*/

void update_parameterAC(
    int iAC,      /* ID NUMBER OF The AC */
    double *x,    /* goma SOLUTION VECTOR */
    double *xdot, /* goma SOLUTION VECTOR TIME DERIVATIVE */
    double *xa,   /* VECTOR OF EXTRA UNKNOWNS associated with augmented conditions */
    Comm_Ex *cx,  /* array of communications structures */
    Exo_DB *exo,  /* ptr to the finite element mesh database */
    Dpi *dpi)     /* distributed processing information */
{
  int mn;
  int ibc, idf;
  int iCC;
  struct Boundary_Condition *BC_Type;

  double lambda, delta, value;

  /*
   * 		BEGIN EXECUTION
   */

  /*
   * augc[iAC] is a AC_Information. We fill in information into AC_Information
   * based on the current value of xa[iAC]
   */

  lambda = xa[iAC];
  augc[iAC].tmp1 = lambda;

  if (augc[iAC].Type == AC_LGRM) {
    augc[iAC].DataFlt[0] = lambda;
    return;
  }

  if (augc[iAC].Type == AC_OVERLAP || augc[iAC].Type == AC_PF_CONSTRAINT) {
    augc[iAC].lm_value = lambda;
  }

  if (augc[iAC].Type == AC_PERIODIC) {
    augc[iAC].lm_value = lambda;
  }

  if (augc[iAC].Type == AC_ARC_LENGTH) {
    delta = lambda - cont->BegParameterValue;

    for (iCC = 0; iCC < nCC; iCC++) {
      value = cpcc[iCC].Beg_CC_Value + cpcc[iCC].ratio * delta;
      update_parameterC(iCC, value, x, xdot, xa, cont->Delta_s0, cx, exo, dpi);
    }
  }

  if (augc[iAC].Type == AC_USERBC || augc[iAC].Type == AC_VOLUME || augc[iAC].Type == AC_FLUX ||
      augc[iAC].Type == AC_LS_VEL || augc[iAC].Type == AC_POSITION || augc[iAC].Type == AC_ANGLE) {

    ibc = augc[iAC].BCID;
    idf = augc[iAC].DFID;

    /*  case for user defined parameter  */

#ifdef GOMA_ENABLE_APREPRO_LIB
    if (ibc == APREPRO_LIB_AC_BCID) {
      goma_error err = aprepro_parse_goma_augc(&augc[iAC], lambda);
      GOMA_EH(err, "Issue with aprepro augmenting condition parsing");
      augc[iAC].DataFlt[0] = lambda;
    } else if (ibc == APREPRO_AC_BCID) {
#else
    if (ibc == APREPRO_AC_BCID) {
#endif // GOMA_ENABLE_APREPRO_LIB
#ifndef tflop
      int err;
      FILE *jfp = NULL;
      char cmd_str[MAX_SYSTEM_COMMAND_LENGTH];
      double temp, lambda_user;
      int ibc_user, idf_user, count;

      if (Num_Proc != 1) {
        GOMA_EH(GOMA_ERROR, "aprepro AC condition not ready for parallel");
      }
      sprintf(cmd_str, "%s %s %s %s %s %s %.20g", "bcdiff.pl", "-p", augc[iAC].Params_File, "-i",
              Input_File, augc[iAC].AP_param, lambda);
      err = system(cmd_str);
      GOMA_EH(err, "Error could not create process for bcdiff.pl");

      augc[iAC].DataFlt[0] = lambda;

      jfp = fopen("tmp.bcdiff", "r");

      ibc_user = 0;
      count = 0;
      while ((ibc_user != -1) && (count < 50)) {
        count++;
        if (fscanf(jfp, "%d", &ibc_user) != 1) {
          GOMA_EH(GOMA_ERROR, "error reading bcdiff");
        }
        if (ibc_user != -1) {
          if (fscanf(jfp, "%d %lf %lf", &idf_user, &temp, &lambda_user) != 3) {
            GOMA_EH(GOMA_ERROR, "error reading bcdiff");
          }
          switch (BC_Types[ibc_user].BC_Name) {
          case SPLINE_BC:
          case SPLINEX_BC:
          case SPLINEY_BC:
          case SPLINEZ_BC:
          case SPLINE_RS_BC:
          case SPLINEX_RS_BC:
          case SPLINEY_RS_BC:
          case SPLINEZ_RS_BC:
          case FILLET_BC:
          case DOUBLE_RAD_BC:
          case FEATURE_ROLLON_BC:
          case ROLL_FLUID_BC:
          case UVARY_BC:
          case VVARY_BC:
          case WVARY_BC:
          case U_PARABOLA_BC:
          case V_PARABOLA_BC:
          case W_PARABOLA_BC:
          case PRESSURE_USER_BC:
          case FLOW_PRESS_USER_BC:
          case T_USER_BC:
          case UUSER_BC:
          case VUSER_BC:
          case WUSER_BC:
          case QUSER_BC:
          case DX_USER_BC:
          case DY_USER_BC:
          case DZ_USER_BC:
          case P_LIQ_USER_BC:
          case SH_P_OPEN_USER_BC:
          case VAR_CA_USER_BC:
          case CA_EDGE_OR_FIX_BC:
          case YFLUX_USER_BC:
          case SURFACE_CHARGE_BC:
          case YUSER_BC:
          case FORCE_USER_BC:
          case FORCE_USER_RS_BC:
            BC_Types[ibc_user].u_BC[idf_user] = lambda_user;
            break;
          default:
            BC_Types[ibc_user].BC_Data_Float[idf_user] = lambda_user;
            break;
          } /*  switch loop */
        }   /* if ibc_user  */
      }     /* while ibc_user  */

      /*    now do AC floats if any  */
      ibc_user = 0;
      count = 0;
      while ((ibc_user != -1) && (count < 50)) {
        count++;
        if (fscanf(jfp, "%d", &ibc_user) != 1) {
          GOMA_EH(GOMA_ERROR, "error reading bcdiff");
        }
        if (ibc_user != -1) {
          if (fscanf(jfp, "%d %lf %lf", &idf_user, &temp, &lambda_user) != 3) {
            GOMA_EH(GOMA_ERROR, "error reading bcdiff");
          }
          augc[ibc_user].DataFlt[idf_user] = lambda_user;
        } /* if ibc_user  */
      }   /* while ibc_user  */
      fclose(jfp);
#else
      GOMA_EH(GOMA_ERROR, "aprepro must be run prior to running goma on this platform.");
#endif
    } else {

      switch (BC_Types[ibc].BC_Name) {
      case SPLINE_BC:
      case SPLINEX_BC:
      case SPLINEY_BC:
      case SPLINEZ_BC:
      case SPLINE_RS_BC:
      case SPLINEX_RS_BC:
      case SPLINEY_RS_BC:
      case SPLINEZ_RS_BC:
      case FILLET_BC:
      case DOUBLE_RAD_BC:
      case FEATURE_ROLLON_BC:
      case ROLL_FLUID_BC:
      case UVARY_BC:
      case VVARY_BC:
      case WVARY_BC:
      case U_PARABOLA_BC:
      case V_PARABOLA_BC:
      case W_PARABOLA_BC:
      case PRESSURE_USER_BC:
      case FLOW_PRESS_USER_BC:
      case T_USER_BC:
      case UUSER_BC:
      case VUSER_BC:
      case WUSER_BC:
      case QUSER_BC:
      case DX_USER_BC:
      case DY_USER_BC:
      case DZ_USER_BC:
      case P_LIQ_USER_BC:
      case SH_P_OPEN_USER_BC:
      case VAR_CA_USER_BC:
      case CA_EDGE_OR_FIX_BC:
      case YFLUX_USER_BC:
      case SURFACE_CHARGE_BC:
      case YUSER_BC:
      case FORCE_USER_BC:
      case FORCE_USER_RS_BC:
        BC_Types[ibc].u_BC[idf] = lambda;
        break;
      case TABLE_BC:
      case TABLE_WICV_BC:
      case TABLE_WICS_BC:
        BC_Type = &BC_Types[ibc];
        BC_Type->table->f[idf] = lambda;
        break;
      default:
        BC_Types[ibc].BC_Data_Float[idf] = lambda;
        break;
      }
    }

  } else if (augc[iAC].Type == AC_USERMAT || augc[iAC].Type == AC_FLUX_MAT) {

    mn = map_mat_index(augc[iAC].MTID);

    switch (augc[iAC].MPID) {

      /*
       * General Model Constants
       */

    case TAGC_THERMAL_CONDUCTIVITY:
      mp_glob[mn]->thermal_conductivity = lambda;
      break;

    case TAGC_ACOUSTIC_WAVENUMBER:
      mp_glob[mn]->wave_number = lambda;
      break;

    case TAGC_ACOUSTIC_IMPEDANCE:
      mp_glob[mn]->acoustic_impedance = lambda;
      break;

    case TAGC_ACOUSTIC_ABSORPTION:
      mp_glob[mn]->acoustic_absorption = lambda;
      break;

    case TAGC_REFRACTIVE_INDEX:
      mp_glob[mn]->refractive_index = lambda;
      break;

    case TAGC_LIGHT_ABSORPTION:
      mp_glob[mn]->light_absorption = lambda;
      break;

    case TAGC_EXTINCTION_INDEX:
      mp_glob[mn]->extinction_index = lambda;
      break;

    case TAGC_ELECTRICAL_CONDUCTIVITY:
      mp_glob[mn]->electrical_conductivity = lambda;
      break;

    case TAGC_VISCOSITY:
      mp_glob[mn]->viscosity = lambda;
      break;

    case TAGC_SURFACE_TENSION:
      mp_glob[mn]->surface_tension = lambda;
      break;

    case TAGC_HEAT_CAPACITY:
      mp_glob[mn]->heat_capacity = lambda;
      break;

    case TAGC_VOLUME_EXPANSION:
      mp_glob[mn]->Volume_Expansion = lambda;
      break;

    case TAGC_DENSITY:
      mp_glob[mn]->density = lambda;
      break;

    case TAGC_POROSITY:
      mp_glob[mn]->porosity = lambda;
      break;

    case TAGC_PERMEABILITY:
      mp_glob[mn]->permeability = lambda;
      break;

    case TAGC_REL_GAS_PERM:
      mp_glob[mn]->rel_gas_perm = lambda;
      break;

    case TAGC_REL_LIQ_PERM:
      mp_glob[mn]->rel_liq_perm = lambda;
      break;

    case TAGC_SATURATION:
      mp_glob[mn]->saturation = lambda;
      break;

    case TAGC_MELTING_POINT_LIQUIDUS:
      mp_glob[mn]->melting_point_liquidus = lambda;
      break;

    case TAGC_MELTING_POINT_SOLIDUS:
      mp_glob[mn]->melting_point_solidus = lambda;
      break;

    case TAGC_FLOWINGLIQUID_VISCOSITY:
      mp_glob[mn]->FlowingLiquid_viscosity = lambda;
      break;

    case TAGC_DIFFUSIVITY_0:
      mp_glob[mn]->diffusivity[0] = lambda;
      break;

    case TAGC_DIFFUSIVITY_1:
      mp_glob[mn]->diffusivity[1] = lambda;
      break;

      /*
       * Generalized Newtonian Models:
       * Newtonian, Power Law, Carreau or Bingham(1,2,3)
       */

    case TAGC_MU0:
      gn_glob[mn]->mu0 = lambda;
      break;
    case TAGC_NEXP:
      gn_glob[mn]->nexp = lambda;
      break;
    case TAGC_MUINF:
      gn_glob[mn]->muinf = lambda;
      break;
    case TAGC_LAM:
      gn_glob[mn]->lam = lambda;
      break;
    case TAGC_AEXP:
      gn_glob[mn]->aexp = lambda;
      break;
    case TAGC_ATEXP:
      gn_glob[mn]->atexp = lambda;
      break;

      /*  CARREAU_WLF  */

    case TAGC_WLFC2:
      gn_glob[mn]->wlfc2 = lambda;
      break;

    case TAGC_REFTEMP:
      mp_glob[mn]->reference[TEMPERATURE] = lambda;
      break;

      /* these are for the BINGHAM yielding material model */

    case TAGC_TAU_Y:
      gn_glob[mn]->tau_y = lambda;
      break;
    case TAGC_FEXP:
      gn_glob[mn]->fexp = lambda;
      break;

      /* these are for SUSPENSION/FILLED_EPOXY models */

    case TAGC_MAXPACK:
      gn_glob[mn]->maxpack = lambda;
      break;
    case TAGC_FICKDIFF_X:
      mp_glob[mn]->u_fdiffusivity[0][0] = lambda;
      break;
    case TAGC_FICKDIFF_Y:
      mp_glob[mn]->u_fdiffusivity[0][1] = lambda;
      break;

      /* these are for CURE/EPOXY/FILLED_EPOXY models */

    case TAGC_GELPOINT:
      gn_glob[mn]->gelpoint = lambda;
      break;
    case TAGC_CUREAEXP:
      gn_glob[mn]->cureaexp = lambda;
      break;
    case TAGC_CUREBEXP:
      gn_glob[mn]->curebexp = lambda;
      break;

      /*
       * Constants used in the Viscoelastic Constitutive Equations
       */

    case TAGC_TIME_CONST:
    case TAGC_TIME_CONST1:
    case TAGC_TIME_CONST2:
    case TAGC_TIME_CONST3:
    case TAGC_TIME_CONST4:
    case TAGC_TIME_CONST5:
    case TAGC_TIME_CONST6:
    case TAGC_TIME_CONST7:
      ve_glob[mn][augc[iAC].MPID - TAGC_TIME_CONST]->time_const = lambda;
      break;

    case TAGC_WT_FUNC:
      vn_glob[mn]->wt_func = lambda;
      break;

    case TAGC_ALPHA:
    case TAGC_ALPHA1:
    case TAGC_ALPHA2:
    case TAGC_ALPHA3:
    case TAGC_ALPHA4:
    case TAGC_ALPHA5:
    case TAGC_ALPHA6:
    case TAGC_ALPHA7:
      ve_glob[mn][augc[iAC].MPID - TAGC_ALPHA]->alpha = lambda;
      break;

    case TAGC_PTT_XI:
      ve_glob[mn][augc[iAC].MPID - TAGC_PTT_XI]->xi = lambda;
      break;

    case TAGC_PTT_EPS:
      ve_glob[mn][augc[iAC].MPID - TAGC_PTT_EPS]->eps = lambda;
      break;

    case TAGC_SHIFT_FUNC:
    case TAGC_SHIFT_FUNC1:
      vn_glob[mn]->shift[augc[iAC].MPID - TAGC_SHIFT_FUNC] = lambda;
      break;

    case TAGC_POLYMER_YIELD_STRESS:
      ve_glob[mn][augc[iAC].MPID - TAGC_POLYMER_YIELD_STRESS]->gn->tau_y = lambda;
      break;

    case TAGC_POLYMER_YIELD_EXPONENT:
      ve_glob[mn][augc[iAC].MPID - TAGC_POLYMER_YIELD_EXPONENT]->gn->fexp = lambda;
      break;

      /*
       * Constants used in the Elasticity Constitutive Equations
       */

    case TAGC_LAME_MU:
      elc_glob[mn]->lame_mu = lambda;
      break;

    case TAGC_LAME_MU_CONTACT_LINE_G0:
      *(elc_glob[mn]->u_mu + 1) = lambda;
      break;

    case TAGC_LAME_MU_CONTACT_LINE_G1:
      *(elc_glob[mn]->u_mu + 2) = lambda;
      break;

    case TAGC_LAME_MU_CONTACT_LINE_R0:
      *(elc_glob[mn]->u_mu + 3) = lambda;
      break;

    case TAGC_LAME_LAMBDA:
      elc_glob[mn]->lame_lambda = lambda;
      break;

    case TAGC_BEND_STIFFNESS:
      elc_glob[mn]->bend_stiffness = lambda;
      break;

    case TAGC_CONV_LAG_VELX:
      *(elc_glob[mn]->v_mesh_sfs) = lambda;
      break;

    case TAGC_CONV_LAG_VELY:
      *(elc_glob[mn]->v_mesh_sfs + 1) = lambda;
      break;

    case TAGC_CONV_LAG_VELZ:
      *(elc_glob[mn]->v_mesh_sfs + 2) = lambda;
      break;

    case TAGC_CONV_LAG_ROTRATE:
      *(elc_glob[mn]->u_v_mesh_sfs) = lambda;
      break;

    case TAGC_CONV_LAG_ROT_X0:
      *(elc_glob[mn]->u_v_mesh_sfs + 1) = lambda;
      break;

    case TAGC_CONV_LAG_ROT_Y0:
      *(elc_glob[mn]->u_v_mesh_sfs + 2) = lambda;
      break;

    case TAGC_CONV_LAG_ROT_Z0:
      *(elc_glob[mn]->u_v_mesh_sfs + 3) = lambda;
      break;

    case TAGC_RS_LAME_MU:
      elc_rs_glob[mn]->lame_mu = lambda;
      break;

    case TAGC_RS_LAME_LAMBDA:
      elc_rs_glob[mn]->lame_lambda = lambda;
      break;

    case TAGC_RS_CONV_LAG_VELX:
      *(elc_rs_glob[mn]->v_mesh_sfs) = lambda;
      break;

    case TAGC_RS_CONV_LAG_VELY:
      *(elc_rs_glob[mn]->v_mesh_sfs + 1) = lambda;
      break;

    case TAGC_RS_CONV_LAG_VELZ:
      *(elc_rs_glob[mn]->v_mesh_sfs + 2) = lambda;
      break;

    case TAGC_RS_CONV_LAG_ROTRATE:
      *(elc_rs_glob[mn]->u_v_mesh_sfs) = lambda;
      break;

    case TAGC_RS_CONV_LAG_ROT_X0:
      *(elc_rs_glob[mn]->u_v_mesh_sfs + 1) = lambda;
      break;

    case TAGC_RS_CONV_LAG_ROT_Y0:
      *(elc_rs_glob[mn]->u_v_mesh_sfs + 2) = lambda;
      break;

    case TAGC_RS_CONV_LAG_ROT_Z0:
      *(elc_rs_glob[mn]->u_v_mesh_sfs + 3) = lambda;
      break;

    case TAGC_POISSON:
      elc_glob[mn]->poisson = lambda;
      break;

    case TAGC_STRSS_FR_SOL_VOL_FRAC:
      elc_glob[mn]->Strss_fr_sol_vol_frac = lambda;
      break;

      /*
       * Constants used for Source Term Models
       */

    case TAGC_NSS_A0:
      mp_glob[mn]->momentum_source[0] = lambda;
      break;

    case TAGC_NSS_A1:
      mp_glob[mn]->momentum_source[1] = lambda;
      break;

    case TAGC_NSS_A2:
      mp_glob[mn]->momentum_source[2] = lambda;
      break;

    case TAGC_NSS_A3:
      mp_glob[mn]->u_momentum_source[0] = lambda;
      break;

    case TAGC_LUB_HGT_U0:
      mp_glob[mn]->u_heightU_function_constants[0] = lambda;
      break;

    case TAGC_LUB_HGT_U1:
      mp_glob[mn]->u_heightU_function_constants[1] = lambda;
      break;

    case TAGC_LUB_HGT_U2:
      mp_glob[mn]->u_heightU_function_constants[2] = lambda;
      break;

    case TAGC_LUB_HGT_U3:
      mp_glob[mn]->u_heightU_function_constants[3] = lambda;
      break;

    case TAGC_LUB_HGT_U4:
      mp_glob[mn]->u_heightU_function_constants[4] = lambda;
      break;

    case TAGC_LUB_HGT_U5:
      mp_glob[mn]->u_heightU_function_constants[5] = lambda;
      break;

    case TAGC_LUB_HGT_U6:
      mp_glob[mn]->u_heightU_function_constants[6] = lambda;
      break;

    case TAGC_LUB_HGT_U7:
      mp_glob[mn]->u_heightU_function_constants[7] = lambda;
      break;

    case TAGC_LUB_HGT_L0:
      mp_glob[mn]->u_heightL_function_constants[0] = lambda;
      break;

    case TAGC_LUB_HGT_L1:
      mp_glob[mn]->u_heightL_function_constants[1] = lambda;
      break;

    case TAGC_LUB_HGT_L2:
      mp_glob[mn]->u_heightL_function_constants[2] = lambda;
      break;

    case TAGC_LUB_HGT_L3:
      mp_glob[mn]->u_heightL_function_constants[3] = lambda;
      break;

    case TAGC_LUB_HGT_L4:
      mp_glob[mn]->u_heightL_function_constants[4] = lambda;
      break;

    case TAGC_LUB_HGT_L5:
      mp_glob[mn]->u_heightL_function_constants[5] = lambda;
      break;

    case TAGC_LUB_HGT_L6:
      mp_glob[mn]->u_heightL_function_constants[6] = lambda;
      break;

    case TAGC_LUB_HGT_L7:
      mp_glob[mn]->u_heightL_function_constants[7] = lambda;
      break;

    case TAGC_U_LUB_VELO_U0:
      mp_glob[mn]->u_veloU_function_constants[0] = lambda;
      break;

    case TAGC_U_LUB_VELO_U1:
      mp_glob[mn]->u_veloU_function_constants[1] = lambda;
      break;

    case TAGC_U_LUB_VELO_U2:
      mp_glob[mn]->u_veloU_function_constants[2] = lambda;
      break;

    case TAGC_U_LUB_VELO_U3:
      mp_glob[mn]->u_veloU_function_constants[3] = lambda;
      break;

    case TAGC_U_LUB_VELO_U4:
      mp_glob[mn]->u_veloU_function_constants[4] = lambda;
      break;

    case TAGC_U_LUB_VELO_U5:
      mp_glob[mn]->u_veloU_function_constants[5] = lambda;
      break;

    case TAGC_U_LUB_VELO_L0:
      mp_glob[mn]->u_veloL_function_constants[0] = lambda;
      break;

    case TAGC_U_LUB_VELO_L1:
      mp_glob[mn]->u_veloL_function_constants[1] = lambda;
      break;

    case TAGC_U_LUB_VELO_L2:
      mp_glob[mn]->u_veloL_function_constants[2] = lambda;
      break;

    case TAGC_U_LUB_VELO_L3:
      mp_glob[mn]->u_veloL_function_constants[3] = lambda;
      break;

    case TAGC_U_LUB_VELO_L4:
      mp_glob[mn]->u_veloL_function_constants[4] = lambda;
      break;

    case TAGC_U_LUB_VELO_L5:
      mp_glob[mn]->u_veloL_function_constants[5] = lambda;
      break;

    case TAGC_LUB_VELO_U0:
      mp_glob[mn]->veloU[0] = lambda;
      break;

    case TAGC_LUB_VELO_U1:
      mp_glob[mn]->veloU[1] = lambda;
      break;

    case TAGC_LUB_VELO_U2:
      mp_glob[mn]->veloU[2] = lambda;
      break;

    case TAGC_LUB_VELO_L0:
      mp_glob[mn]->veloL[0] = lambda;
      break;

    case TAGC_LUB_VELO_L1:
      mp_glob[mn]->veloL[1] = lambda;
      break;

    case TAGC_LUB_VELO_L2:
      mp_glob[mn]->veloL[2] = lambda;
      break;

    case TAGC_LUB_DCA_U0:
      mp_glob[mn]->u_dcaU_function_constants[0] = lambda;
      break;

    case TAGC_LUB_DCA_U1:
      mp_glob[mn]->u_dcaU_function_constants[1] = lambda;
      break;

    case TAGC_LUB_DCA_U2:
      mp_glob[mn]->u_dcaU_function_constants[2] = lambda;
      break;

    case TAGC_LUB_DCA_U3:
      mp_glob[mn]->u_dcaU_function_constants[3] = lambda;
      break;

    case TAGC_LUB_DCA_L0:
      mp_glob[mn]->u_dcaL_function_constants[0] = lambda;
      break;

    case TAGC_LUB_DCA_L1:
      mp_glob[mn]->u_dcaL_function_constants[1] = lambda;
      break;

    case TAGC_LUB_DCA_L2:
      mp_glob[mn]->u_dcaL_function_constants[2] = lambda;
      break;

    case TAGC_LUB_DCA_L3:
      mp_glob[mn]->u_dcaL_function_constants[3] = lambda;
      break;

    case TAGC_LUB_SOURCE_0:
      mp_glob[mn]->u_lubsource_function_constants[0] = lambda;
      break;

    case TAGC_LUB_SOURCE_1:
      mp_glob[mn]->u_lubsource_function_constants[1] = lambda;
      break;

    case TAGC_LUB_SOURCE_2:
      mp_glob[mn]->u_lubsource_function_constants[2] = lambda;
      break;

    case TAGC_RST_FUNC_0:
      mp_glob[mn]->Rst_func = lambda;
      break;

    case TAGC_RST_FUNC_1:
      mp_glob[mn]->Rst_diffusion = lambda;
      break;

    case TAGC_RST_FUNC_2:
      mp_glob[mn]->Rst_func_supg = lambda;
      break;

    case TAGC_HEAT_SOURCE_0:
      mp_glob[mn]->u_heat_source[0] = lambda;
      break;

    case TAGC_SPECIES_SOURCE_0_P0:
      mp_glob[mn]->u_species_source[0][0] = lambda;
      break;

    case TAGC_SPECIES_SOURCE_0_P1:
      mp_glob[mn]->u_species_source[0][1] = lambda;
      break;

    case TAGC_SPECIES_SOURCE_0_P2:
      mp_glob[mn]->u_species_source[0][2] = lambda;
      break;

    case TAGC_SPECIES_SOURCE_0_P3:
      mp_glob[mn]->u_species_source[0][3] = lambda;
      break;

    case TAGC_SPECIES_SOURCE_1_P0:
      mp_glob[mn]->u_species_source[1][0] = lambda;
      break;

    case TAGC_SPECIES_SOURCE_1_P1:
      mp_glob[mn]->u_species_source[1][1] = lambda;
      break;

    case TAGC_SPECIES_SOURCE_1_P2:
      mp_glob[mn]->u_species_source[1][2] = lambda;
      break;

    case TAGC_SPECIES_SOURCE_1_P3:
      mp_glob[mn]->u_species_source[1][3] = lambda;
      break;

    case TAGC_LATENT_HEAT_0:
      mp_glob[mn]->latent_heat_vap[0] = lambda;
      break;

    case TAGC_LATENT_HEAT_1:
      mp_glob[mn]->latent_heat_vap[1] = lambda;
      break;

    default:
      printf("\n\t Error: Invalid Material Property Tag %d\n", augc[iAC].MPID);
      exit(0);
      break;
    }
  }

} /* END of routine update_parameterAC  */
/*****************************************************************************/

//! Routine called from sol_nonlinear() that handles setting up the expanded
//! matrix problem involved with augmented conditions
//! added by users.
int user_aug_cond(int iAC,
                  int nAC,
                  double x_AC[],
                  double **bAC,
                  double **cAC,
                  double **dAC,
                  double *gAC,
                  int numProcUnknowns,
                  double **x_sens_p,
                  Comm_Ex *cx,
                  MF_Args *mf_args) {

  /*  flags to determine if analytical AC matrices are available */

  int jAC = 0, have_bAC = FALSE, have_cAC = FALSE, have_dAC = FALSE;

  int i, err = -99;
  const int Jac_state = af->Assemble_Jacobian;

  double p_save, dp_save = 0.0;
  double fd_factor = FD_FACTOR;

  double *res_p, *res_m, xm, *x = mf_args->x;

  asdv(&res_p, numProcUnknowns);
  asdv(&res_m, numProcUnknowns);

  /*
    LOAD EXTRA UNKNOWNS INTO x_AC
  */

  load_extra_unknownsAC(iAC, x_AC, cx, mf_args->exo, mf_args->dpi);

  /*

  RESIDUALS OF THE AUGMENTING CONDITION

  */

  user_aug_cond_residuals(iAC, mf_args->x, mf_args->xdot, *(mf_args->delta_t), *(mf_args->time),
                          x_sens_p, gAC, &have_bAC, &have_cAC, &have_dAC, bAC, cAC, dAC,
                          mf_args->exo, mf_args->dpi, cx);

  /*

  B Matrix     dR/dp

  */

  if (!have_bAC) {
    /* Compute bAC with numerical difference
     *
     */

    p_save = x_AC[iAC];
    /*  */
    dp_save = fd_factor * p_save;
    dp_save = (fabs(dp_save) < fd_factor ? fd_factor : dp_save);
    /*  */
    x_AC[iAC] = p_save + dp_save;

    update_parameterAC(iAC, mf_args->x, mf_args->xdot, x_AC, cx, mf_args->exo, mf_args->dpi);

    init_vec_value(res_p, 0.0, numProcUnknowns);
    af->Assemble_Residual = TRUE;
    af->Assemble_Jacobian = FALSE;
    af->Assemble_LSA_Jacobian_Matrix = FALSE;
    af->Assemble_LSA_Mass_Matrix = FALSE;

    err = matrix_fill_full(mf_args->ams, mf_args->x, res_p, mf_args->x_old, mf_args->x_older,
                           mf_args->xdot, mf_args->xdot_old, mf_args->x_update, mf_args->delta_t,
                           mf_args->theta_, mf_args->first_elem_side_bc, mf_args->time,
                           mf_args->exo, mf_args->dpi, mf_args->num_total_nodes,
                           mf_args->h_elem_avg, mf_args->U_norm, mf_args->estifm);

    if (err == -1)
      return (err);
    /*  */
    x_AC[iAC] = p_save - dp_save;

    update_parameterAC(iAC, mf_args->x, mf_args->xdot, x_AC, cx, mf_args->exo, mf_args->dpi);

    init_vec_value(res_m, 0.0, numProcUnknowns);

    af->Assemble_Residual = TRUE;
    af->Assemble_Jacobian = FALSE;
    af->Assemble_LSA_Jacobian_Matrix = FALSE;
    af->Assemble_LSA_Mass_Matrix = FALSE;

    err = matrix_fill_full(mf_args->ams, mf_args->x, res_m, mf_args->x_old, mf_args->x_older,
                           mf_args->xdot, mf_args->xdot_old, mf_args->x_update, mf_args->delta_t,
                           mf_args->theta_, mf_args->first_elem_side_bc, mf_args->time,
                           mf_args->exo, mf_args->dpi, mf_args->num_total_nodes,
                           mf_args->h_elem_avg, mf_args->U_norm, mf_args->estifm);

    if (err == -1)
      return (err);
    /*  */
    xm = 0.5 / dp_save;
    /*  */
    v2sum(numProcUnknowns, &bAC[iAC][0], xm, &res_p[0], -xm, &res_m[0]);
    /*  */
    x_AC[iAC] = p_save;

    update_parameterAC(iAC, mf_args->x, mf_args->xdot, x_AC, cx, mf_args->exo, mf_args->dpi);

  } /*  end of if have_bAC */

  /*

  C Matrix     dN/dx

  */

  if (!have_cAC) {
    /*
     * Compute cAC with numerical differencing
     */

    for (i = 0; i < numProcUnknowns; i++) {
      p_save = x[i];
      /*  */
      x[i] = p_save + dp_save;

      user_aug_cond_residuals(iAC, mf_args->x, mf_args->xdot, *(mf_args->delta_t), *(mf_args->time),
                              x_sens_p, res_p, &have_bAC, &have_cAC, &have_dAC, bAC, cAC, dAC,
                              mf_args->exo, mf_args->dpi, cx);
      /*  */
      x[i] = p_save - dp_save;

      user_aug_cond_residuals(iAC, mf_args->x, mf_args->xdot, *(mf_args->delta_t), *(mf_args->time),
                              x_sens_p, res_m, &have_bAC, &have_cAC, &have_dAC, bAC, cAC, dAC,
                              mf_args->exo, mf_args->dpi, cx);
      /*  */
      xm = 0.5 / dp_save;
      /*  */

      cAC[iAC][i] = xm * (res_p[iAC] - res_m[iAC]);

      /*  */
      x[i] = p_save;
    }
  } /*  end of if !have_cAC */

  /*

  D Matrix     dN/dp

  */

  if (!have_dAC) {

    for (jAC = 0; jAC < nAC; jAC++) {
      p_save = x_AC[jAC];

      dp_save = fd_factor * p_save;
      dp_save = (fabs(dp_save) < fd_factor ? fd_factor : dp_save);
      /*  */
      x_AC[jAC] = p_save + dp_save;

      update_parameterAC(jAC, mf_args->x, mf_args->xdot, x_AC, cx, mf_args->exo, mf_args->dpi);

      user_aug_cond_residuals(iAC, mf_args->x, mf_args->xdot, *(mf_args->delta_t), *(mf_args->time),
                              x_sens_p, res_p, &have_bAC, &have_cAC, &have_dAC, bAC, cAC, dAC,
                              mf_args->exo, mf_args->dpi, cx);
      /*  */
      x_AC[jAC] = p_save - dp_save;

      update_parameterAC(jAC, mf_args->x, mf_args->xdot, x_AC, cx, mf_args->exo, mf_args->dpi);

      user_aug_cond_residuals(iAC, mf_args->x, mf_args->xdot, *(mf_args->delta_t), *(mf_args->time),
                              x_sens_p, res_m, &have_bAC, &have_cAC, &have_dAC, bAC, cAC, dAC,
                              mf_args->exo, mf_args->dpi, cx);
      /*  */
      xm = 0.5 / dp_save;
      /*  */
      dAC[iAC][jAC] = xm * (res_p[iAC] - res_m[iAC]);

      /*  */
      x_AC[jAC] = p_save;

      update_parameterAC(jAC, mf_args->x, mf_args->xdot, x_AC, cx, mf_args->exo, mf_args->dpi);
    }

  } /* end of if have_dAC  */

  safe_free((void *)res_p);
  safe_free((void *)res_m);
  af->Assemble_Jacobian = Jac_state;

  return (TRUE);
}

//! Routine called from sol_nonlinear() that handles calculating all
//! of the extra matrix and residual components
//! involved with augmented conditions that are "standard"
/*!
 *            AC_FLUX:
 *	      AC_VOLUME:
 *	      AC_LS_VEL:
 *	      AC_POSITION:
 *	      AC_ANGLE:
 *            AC_FLUX_MAT:
 *
 * Output Variables
 * ---------------------
 *   x_AC[iAC]   Value of the unknowns in the augmented condition
 *   gAC[iAC]    Value of the Residual for the iAC augmented condition
 *   cAC[iAC][i] Derivative of the AugResidual equation wrt to the
 *               ith unknown on this processor.
 *   bAC[iAC][i]  Derivative of the ith residual on this processor
 *                wrt to the iAC augmented condition
 *   dAC[iAC][jAC] Derivative of the iAC AugResidual equation wrt to the
 *                 jAC unknown corresponding to the jAC aug condition
 */
int std_aug_cond(int iAC,
                 int nAC,
                 double x_AC[],
                 double **bAC,
                 double **cAC,
                 double **dAC,
                 double *gAC,
                 int numProcUnknowns,
                 Comm_Ex *cx,
                 MF_Args *mf_args) {
  double inventory;
  double LSvel_inventory;
  int i, jAC;
#ifdef PARALLEL
  double global_inventory = 0.0;
  double global_LSvel_inventory = 0.0;
#endif

  /*
   *  LOAD EXTRA UNKNOWNS INTO x_AC
   *    Get the value of the unknown corresponding to the augmented condition
   *    and load it into the vector of Augmented condition unknowns, x_AC[]
   */
  load_extra_unknownsAC(iAC, x_AC, cx, mf_args->exo, mf_args->dpi);

  /*
   * dN/dx  This is the sensitivity of the flux constraint with respect all the
   *        FEM degrees of freedom ( i.e. everything except the augmenting parameter )
   */
  if (augc[iAC].Type == AC_FLUX) {
    inventory = evaluate_flux(mf_args->exo, mf_args->dpi, augc[iAC].SSID, augc[iAC].MFID, NULL,
                              augc[iAC].MTID, augc[iAC].COMPID, NULL, FALSE, mf_args->x,
                              mf_args->xdot, cAC[iAC], *(mf_args->delta_t), *(mf_args->time), 0);
    if (augc[iAC].len_AC > 2) {
      int mfid, ssid, mtid, compid;
      double inventory1, ac_factor = 1.0;

      mfid = (int)augc[iAC].DataFlt[0];
      ssid = (int)augc[iAC].DataFlt[1];
      mtid = (int)augc[iAC].DataFlt[2];
      if (augc[iAC].len_AC > 3) {
        compid = (int)augc[iAC].DataFlt[3];
      } else {
        compid = 0;
      }
      if (augc[iAC].len_AC > 4) {
        ac_factor = augc[iAC].DataFlt[4];
      }
      inventory1 = evaluate_flux(mf_args->exo, mf_args->dpi, ssid, mfid, NULL, mtid, compid, NULL,
                                 FALSE, mf_args->x, mf_args->xdot, cAC[iAC], *(mf_args->delta_t),
                                 *(mf_args->time), 0);
      inventory += ac_factor * inventory1;
    }
    /*
     * And last of all. set the diagonal sensitivities to zero (but not for Level Set Velocity AC)
     */
    for (jAC = 0; jAC < nAC; jAC++) {
      dAC[iAC][jAC] = 0.0;
    }
    // Formulate and store the residual for the augmented condition
    gAC[iAC] = (inventory - augc[iAC].CONSTV);
  } else if (augc[iAC].Type == AC_FLUX_MAT) {
    inventory = evaluate_flux(mf_args->exo, mf_args->dpi, augc[iAC].SSID, augc[iAC].MFID, NULL,
                              augc[iAC].MDID, augc[iAC].COMPID, NULL, FALSE, mf_args->x,
                              mf_args->xdot, cAC[iAC], *(mf_args->delta_t), *(mf_args->time), 0);

    /*
     * And last of all. set the diagonal sensitivities to zero (but not for Level Set Velocity AC)
     */
    for (jAC = 0; jAC < nAC; jAC++) {
      dAC[iAC][jAC] = 0.0;
    }
    // Formulate and store the residual for the augmented condition
    gAC[iAC] = (inventory - augc[iAC].CONSTV);
  } else if (augc[iAC].Type == AC_VOLUME) {
    if (augc[iAC].VOLID >= 10) {
      user_aug_cond_volume_residuals(iAC, mf_args->x, mf_args->xdot, *(mf_args->delta_t),
                                     *(mf_args->time), x_AC, gAC, cAC, dAC, numProcUnknowns,
                                     mf_args->exo, mf_args->dpi, cx);
    } else {
      // I guess we stored the volume calculation in evol on a previous
      // step. Here we combine the contributions across processors and come up
      // with a global unknown.
      inventory = augc[iAC].evol;
#ifdef PARALLEL
      if (Num_Proc > 1) {
        MPI_Allreduce(&inventory, &global_inventory, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        inventory = global_inventory;
      }
#endif
      // Formulate and store the residual for the augmented condition
      gAC[iAC] = (inventory - augc[iAC].CONSTV);

      for (i = 0; i < numProcUnknowns; i++) {
        cAC[iAC][i] = augc[iAC].d_evol_dx[i];
      }
      for (jAC = 0; jAC < nAC; jAC++) {
        dAC[iAC][jAC] = 0.0;
      }
    }
  } else if (augc[iAC].Type == AC_LS_VEL) {
    if (augc[iAC].d_lsvol_dx == 0)
      DPRINTF(stderr, "\n can't see vectors ");

    inventory = augc[iAC].lsvol;
    LSvel_inventory = augc[iAC].lsvel;

#ifdef PARALLEL
    if (Num_Proc > 1) {
      MPI_Allreduce(&inventory, &global_inventory, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&LSvel_inventory, &global_LSvel_inventory, 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD);
      inventory = global_inventory;
      LSvel_inventory = global_LSvel_inventory;
    }

#endif /* PARALLEL */

    for (i = 0; i < numProcUnknowns; i++) {
      cAC[iAC][i] = (augc[iAC].d_lsvel_dx[i] / inventory -
                     augc[iAC].d_lsvol_dx[i] * augc[iAC].lsvel / inventory / inventory);
    }

    inventory = LSvel_inventory / inventory;

    /*
     * And last of all. set the diagonal sensitivities to zero (but not for Level Set Velocity AC)
     */
    for (jAC = 0; jAC < nAC; jAC++) {
      dAC[iAC][jAC] = 0.0;
    }
    // Formulate and store the residual for the augmented condition
    gAC[iAC] = (inventory - augc[iAC].CONSTV);

  } else if (augc[iAC].Type == AC_POSITION) {
    inventory = getPositionAC(augc + iAC, cAC[iAC], mf_args->x, mf_args->exo);
#ifdef PARALLEL
    if (Num_Proc > 1) {
      MPI_Allreduce(&inventory, &global_inventory, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      inventory = global_inventory;
    }
#endif
    /*
     * And last of all. set the diagonal sensitivities to zero (but not for Level Set Velocity AC)
     */
    for (jAC = 0; jAC < nAC; jAC++) {
      dAC[iAC][jAC] = 0.0;
    }
    // Formulate and store the residual for the augmented condition
    gAC[iAC] = (inventory - augc[iAC].CONSTV);
  } else if (augc[iAC].Type == AC_ANGLE) {
    inventory = getAngleAC(augc + iAC, cAC[iAC], mf_args->x, mf_args->exo);
#ifdef PARALLEL
    if (Num_Proc > 1) {
      MPI_Allreduce(&inventory, &global_inventory, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      inventory = global_inventory;
    }
#endif
    /*
     * And last of all. set the diagonal sensitivities to zero (but not for Level Set Velocity AC)
     */
    for (jAC = 0; jAC < nAC; jAC++) {
      dAC[iAC][jAC] = 0.0;
    }
    // Formulate and store the residual for the augmented condition
    gAC[iAC] = (inventory - augc[iAC].CONSTV);
  }

  // Time dependent augmenting conditions. This small section provides an
  // example for turning on a time dependent
  // augmenting condition. Seems to work fine, and is less disruptive.
  // Cp is an adjustable constant, to set the time scale.
  //  if (TimeIntegration != 0) {
  //  double Cp = augc[iAC].LewisNum;
  //  double tt  = 0.0;
  //  gAC[iAC] = Cp * augc[iAC].tmp2 + (inventory -  augc[iAC].CONSTV);
  //  dAC[iAC][iAC] += Cp * (1. + 2.*tt) / *(mf_args->delta_t);
  //}

  /*
   * dR/dp    This is the sensitivity of the residual vector with respect to the
   *          additional augmenting parameter
   */
  estimate_bAC(iAC, x_AC, bAC, numProcUnknowns, cx, mf_args);

  return (TRUE);
}

int alc_aug_cond(int iAC,
                 int nAC,
                 double x_AC[],
                 double **bAC,
                 double **cAC,
                 double **dAC,
                 double *gAC,
                 double *equation,
                 int numProcUnknowns,
                 Comm_Ex *cx,
                 struct con_struct *con,
                 MF_Args *mf_args)
/*
 * This routine allows the arc length bordering algorithm to be
 * performed as an additional augmenting condition, thus doing
 * some of the work normally done inside LOCA.
 */
{

  /*  flags to determine if analytical AC matrices are available */

  int i;

  double alceq_s = 0.0, alceq_p = 0.0, resid = 0.0;

  double *x = mf_args->x, *x_old = mf_args->x_old;

  /* Data needed from LOCA's "con" structure */
  double lambda = con->general_info.param;
  double lambda_old = con->private_info.param_old;
  double dS = con->private_info.arc_step;
  double dp_ds = con->private_info.dp_ds;
  double *dx = con->private_info.x_tang;
  double *sv = con->private_info.scale_vec;

  /* Residual is arc length equation from LOCA */

  /* Parameter contribution: dp_ds * (p - p_old) */
  alceq_p = dp_ds * (lambda - lambda_old);

  /* Solution contribution:  dx_ds * (x - x_old) with scale factor */
  for (i = 0; i < NumUnknowns[pg->imtrx]; i++) {
    alceq_s += sv[i] * sv[i] * dx[i] * (x[i] - x_old[i]);
  }
  if (Num_Proc > 1) {
    alceq_s = AZ_gsum_double(alceq_s, mf_args->ams->proc_config);
  }

  /* Construct and load total residual */
  resid = alceq_p + alceq_s - dS;
  *equation = resid;
  gAC[iAC] = resid;

  /* Here, the extra unknown is just the continuation parameter */
  x_AC[iAC] = lambda;

  /* Get bAC by differencing residual with parameter */
  estimate_bAC(iAC, x_AC, bAC, numProcUnknowns, cx, mf_args);

  /* cAC is solution sensitivity from LOCA */
  for (i = 0; i < NumUnknowns[pg->imtrx]; i++) {
    cAC[iAC][i] = dx[i] * sv[i] * sv[i];
  }

  /* Extra right column of dAC is AC residual sensitivity to lambda */
  estimate_dAC_ALC(iAC, x_AC, bAC, numProcUnknowns, cx, mf_args);

  /* Bottom right entry of dAC is dp_ds */
  dAC[iAC][iAC] = dp_ds;

  /* Done */
  return (TRUE);
}

void overlap_aug_cond(int ija[],
                      double a[],
                      double x_AC[],
                      double *gAC,
                      double **bAC,
                      double **cAC,
                      double **dAC,
                      Comm_Ex *cx,
                      MF_Args *mf_args)
/*
 * This function is called just once per iteration to fill in all
 * augmenting condition terms associated with CONTACT_SURF and
 * EMBEDDED_SURF boundary conditions for overlapping grid
 * fluid/solid problems, where space for the terms has not been
 * allocated in the Jacobian array. It calls the functions
 * apply_embedded_bc() and apply_contact_bc() on an element by
 * element basis, which places the terms in the gAC, bAC, and
 * cAC arrays as appropriate.
 */
{
  int err, ibf = -1, ibs = -1, ib, i, jAC;
  int e_start, e_end, ielem, iside, jelem, jside;
  int solid_eb, fluid_eb, lm_eb;
  double h[DIM];
  double mu_avg = 0;
  double *x = mf_args->x;
  Exo_DB *exo = mf_args->exo;
  struct elem_side_bc_struct *elem_side_bc;
  int ac_lm = -1;

  /* Get element block indices for solid, fluid, and Lagrange surface material */
  solid_eb = augc[nAC - 1].solid_eb + 1;
  fluid_eb = augc[nAC - 1].fluid_eb + 1;
  lm_eb = augc[nAC - 1].lm_eb + 1;
  if (solid_eb == fluid_eb)
    GOMA_EH(GOMA_ERROR, "Solid and fluid must be on different blocks!");

  if (Do_Overlap && lm_eb == solid_eb)
    ac_lm = 1;
  else if (Do_Overlap && lm_eb == fluid_eb)
    ac_lm = 2;
  else
    GOMA_EH(GOMA_ERROR, "Shouldn't be here.");

  /* Just fill this in for use as argument */
  for (i = 0; i < DIM; i++) {
    h[i] = 1.0;
  }

  /* Locate fluid and solid blocks in exo structure */
  for (ib = 0; ib < exo->num_elem_blocks; ib++) {
    if (exo->eb_id[ib] == solid_eb)
      ibs = ib;
    if (exo->eb_id[ib] == fluid_eb)
      ibf = ib;
  }

  if (ac_lm == 1) /* AC on solid */
  {
    /* Fill in terms on level set (fluid element block) */
    e_start = exo->eb_ptr[ibf];
    e_end = exo->eb_ptr[ibf + 1];

    for (ielem = e_start; ielem < e_end; ielem++) {

      /* Check for level set going through element (fill unknown sign change) */
      if (elem_overlaps_interface(ielem, x, exo, ls->Length_Scale)) {
        load_elem_dofptr(ielem, exo, x_static, x_old_static, xdot_static, xdot_old_static, 0);

        /* Assemble embedded BC's into gAC, bAC, and cAC */
        if (ls->Length_Scale == 0.) {
          apply_embedded_bc(ielem, x, *mf_args->delta_t, *mf_args->theta_, *mf_args->time, NULL, 0,
                            gAC, bAC, cAC, exo);
        } else {
          apply_distributed_sources(ielem, ls->Length_Scale, x, exo, *mf_args->delta_t,
                                    *mf_args->theta_, *mf_args->time, NULL, 0, gAC, bAC, cAC);
        }
      }
    }

    /*
     * Now handle the CONTACT_SURF type boundary conditions.
     * Proceed by looping over OVERLAP type AC's, of which there
     * are at least two per relevant element side. It should be
     * more efficient to make only one call to apply_contact_bc()
     * per element side, which will fill in all AC entries for
     * that side.
     */
    jelem = -1;
    jside = -1;

    for (jAC = 0; jAC < nAC; jAC++) {
      if (augc[jAC].Type == AC_OVERLAP) {
        ielem = augc[jAC].lm_elem;
        iside = augc[jAC].lm_side;
        if (ielem != jelem || iside != jside) {

          /* This is the next element side that needs to be processed */
          /* If this is a new element, get pointers */
          if (ielem != jelem) {
            load_elem_dofptr(ielem, exo, x_static, x_old_static, xdot_static, xdot_old_static, 0);
          }

          /* Find the appropriate side and point elem_side_bc to it */
          elem_side_bc = First_Elem_Side_BC_Array[pg->imtrx][ielem];
          jside = elem_side_bc->id_side;
          if (iside != jside) {
            do {
              elem_side_bc = elem_side_bc->next_side_bc;
              jside = elem_side_bc->id_side;
            } while (iside != jside);
          }

          /* Now fill in all gAC, bAC, and cAC terms for this side */
          err = apply_contact_bc(
              x, mf_args->resid, *mf_args->delta_t, *mf_args->theta_, *mf_args->h_elem_avg, h,
              mu_avg, *mf_args->U_norm, First_Elem_Side_BC_Array[pg->imtrx], ielem,
              ei[pg->imtrx]->ielem_type, ei[pg->imtrx]->num_local_nodes, ei[pg->imtrx]->ielem_dim,
              Proc_Connect_Ptr[ielem], elem_side_bc, *mf_args->num_total_nodes, CONTACT_SURF, jAC,
              gAC, bAC, cAC, dAC, *mf_args->time, exo);
          GOMA_EH(err, " apply_contact_bc");
        } /* END of "if (ielem != jelem || iside != jside)" */

        jelem = ielem;
        jside = iside;
      }                  /* END of "if (augc[jAC].Type == AC_OVERLAP)" */
    }                    /* END of loop over augmenting conditions (jAC) */
  } else if (ac_lm == 2) /* AC on fluid */
  {
    jelem = -1;
    jside = -1;

    for (jAC = 0; jAC < nAC; jAC++) {
      if (augc[jAC].Type == AC_OVERLAP) {
        ielem = augc[jAC].lm_elem;
        iside = augc[jAC].lm_side;
        if (ielem == -1 || ielem != jelem || iside != jside) {

          /* This is the next element side that needs to be processed */
          /* If this is a new element, get pointers */
          if (ielem >= 0 && ielem != jelem) {
            load_elem_dofptr(ielem, exo, x_static, x_old_static, xdot_static, xdot_old_static, 0);
          }
          if (ielem >= 0 && current_elem_overlaps_interface(ls->Length_Scale)) {
            if (ls->Length_Scale == 0.) {
              apply_embedded_bc(ielem, x, *mf_args->delta_t, *mf_args->theta_, *mf_args->time, NULL,
                                jAC, gAC, bAC, cAC, exo);
            } else {
              apply_distributed_sources(ielem, ls->Length_Scale, x, exo, *mf_args->delta_t,
                                        *mf_args->theta_, *mf_args->time, NULL, jAC, gAC, bAC, cAC);
            }
          } else {
            /* trivialize augmenting conditions in elements that don't cross interface */
            gAC[jAC] = augc[jAC].lm_value;
            dAC[jAC][jAC] = 1.;
          }
        } /* END of "if (ielem != jelem || iside != jside)" */

        jelem = ielem;
        jside = iside;
      } /* END of "if (augc[jAC].Type == AC_OVERLAP)" */
    }   /* END of loop over augmenting conditions (jAC) */

    /* Fill in terms on solid element block */
    e_start = exo->eb_ptr[ibs];
    e_end = exo->eb_ptr[ibs + 1];

    for (ielem = e_start; ielem < e_end; ielem++) {

      if (First_Elem_Side_BC_Array[ielem] != NULL) {
        int call_contact;
        int ibc, bc_input_id, bct;
        elem_side_bc = First_Elem_Side_BC_Array[pg->imtrx][ielem];

        load_elem_dofptr(ielem, exo, x_static, x_old_static, xdot_static, xdot_old_static, 0);
        /*****************************************************************************/
        do { /* begining of do while construct */
          /* which loops over the sides of this element that have boundary
             conditions */

          call_contact = 0;
          for (ibc = 0; (bc_input_id = (int)elem_side_bc->BC_input_id[ibc]) != -1; ibc++) {
            bct = BC_Types[bc_input_id].desc->method;
            if (bct == CONTACT_SURF)
              call_contact = 1;
          }

          if (call_contact) {
            err = apply_contact_bc(
                x, mf_args->resid, *mf_args->delta_t, *mf_args->theta_, *mf_args->h_elem_avg, h,
                mu_avg, *mf_args->U_norm, First_Elem_Side_BC_Array[pg->imtrx], ielem,
                ei[pg->imtrx]->ielem_type, ei[pg->imtrx]->num_local_nodes, ei[pg->imtrx]->ielem_dim,
                Proc_Connect_Ptr[ielem], elem_side_bc, *mf_args->num_total_nodes, CONTACT_SURF, jAC,
                gAC, bAC, cAC, dAC, *mf_args->time, exo);
          }
          /****************************************************************************/
        } while ((elem_side_bc = elem_side_bc->next_side_bc) != NULL);
        /* END of do  while () construct				      */
        /******************************************************************************/
      } /* END if (First_Elem_Side_BC_Array[ielem] != NULL) 		      */
    }
  }
#if 0
   {
     int index;
     int numProcUnknowns = NumUnknowns[pg->imtrx] + NumExtUnknowns[pg->imtrx];
     for (jAC = 0; jAC < nAC; jAC++)
       {
         printf("jAC=%d, LM=%g, gAC=%g, dAC[%d][%d]=%g\n",jAC,augc[jAC].lm_value,gAC[jAC],jAC,jAC,dAC[jAC][jAC]);
         for ( index = 0; index<numProcUnknowns; index++)
           {
             if (bAC[jAC][index] != 0. ) printf(" eqntype=%d, N=%d, bAC[jAC][%d]=%g\n",idv[pg->imtrx][index][0],idv[pg->imtrx][index][2],index,bAC[jAC][index]);
             if (cAC[jAC][index] != 0. ) printf(" vartype=%d, N=%d, cAC[jAC][%d]=%g\n",idv[pg->imtrx][index][0],idv[pg->imtrx][index][2],index,cAC[jAC][index]);
           }
       }
   }
#endif
  /* Done! */
  return;
} /* END of function overlap_aug_cond() */

static int estimate_bAC(
    int iAC, double x_AC[], double **bAC, int numProcUnknowns, Comm_Ex *cx, MF_Args *mf_args) {
  double p_save, dp_save;
  double fd_factor = FD_FACTOR;

  double *res_p, *res_m, xm;
  int err = -99;

  asdv(&res_p, numProcUnknowns);
  asdv(&res_m, numProcUnknowns);

  p_save = x_AC[iAC];
  /*  */
  dp_save = fd_factor * p_save;
  dp_save = (fabs(dp_save) < fd_factor ? fd_factor : dp_save);
  /*  */
  x_AC[iAC] = p_save + dp_save;

  update_parameterAC(iAC, mf_args->x, mf_args->xdot, x_AC, cx, mf_args->exo, mf_args->dpi);

  af->Assemble_Residual = TRUE;
  af->Assemble_Jacobian = FALSE;
  af->Assemble_LSA_Jacobian_Matrix = FALSE;
  af->Assemble_LSA_Mass_Matrix = FALSE;

  if (augc[iAC].Type == AC_VOLUME) {
    init_vec_value(augc[iAC].d_evol_dx, 0.0, numProcUnknowns);
    augc[iAC].evol = 0.;
  }

  if (augc[iAC].Type == AC_LS_VEL) {
    init_vec_value(augc[iAC].d_lsvel_dx, 0.0, numProcUnknowns);
    init_vec_value(augc[iAC].d_lsvol_dx, 0.0, numProcUnknowns);
    augc[iAC].lsvel = 0.;
    augc[iAC].lsvol = 0.;
  }

  err = matrix_fill_full(mf_args->ams, mf_args->x, res_p, mf_args->x_old, mf_args->x_older,
                         mf_args->xdot, mf_args->xdot_old, mf_args->x_update, mf_args->delta_t,
                         mf_args->theta_, mf_args->first_elem_side_bc, mf_args->time, mf_args->exo,
                         mf_args->dpi, mf_args->num_total_nodes, mf_args->h_elem_avg,
                         mf_args->U_norm, mf_args->estifm);

  if (err == -1)
    return (err);

  /*  */
  x_AC[iAC] = p_save - dp_save;

  update_parameterAC(iAC, mf_args->x, mf_args->xdot, x_AC, cx, mf_args->exo, mf_args->dpi);

  af->Assemble_Residual = TRUE;
  af->Assemble_Jacobian = FALSE;
  af->Assemble_LSA_Jacobian_Matrix = FALSE;
  af->Assemble_LSA_Mass_Matrix = FALSE;

  if (augc[iAC].Type == AC_VOLUME) {
    init_vec_value(augc[iAC].d_evol_dx, 0.0, numProcUnknowns);
    augc[iAC].evol = 0.;
  }

  if (augc[iAC].Type == AC_LS_VEL) {
    init_vec_value(augc[iAC].d_lsvel_dx, 0.0, numProcUnknowns);
    init_vec_value(augc[iAC].d_lsvol_dx, 0.0, numProcUnknowns);
    augc[iAC].lsvel = 0.;
    augc[iAC].lsvol = 0.;
  }

  err = matrix_fill_full(mf_args->ams, mf_args->x, res_m, mf_args->x_old, mf_args->x_older,
                         mf_args->xdot, mf_args->xdot_old, mf_args->x_update, mf_args->delta_t,
                         mf_args->theta_, mf_args->first_elem_side_bc, mf_args->time, mf_args->exo,
                         mf_args->dpi, mf_args->num_total_nodes, mf_args->h_elem_avg,
                         mf_args->U_norm, mf_args->estifm);

  if (err == -1)
    return (err);

  xm = 0.5 / dp_save;

  v2sum(numProcUnknowns, &bAC[iAC][0], xm, &res_p[0], -xm, &res_m[0]);

  x_AC[iAC] = p_save;

  update_parameterAC(iAC, mf_args->x, mf_args->xdot, x_AC, cx, mf_args->exo, mf_args->dpi);

  af->Assemble_Jacobian = TRUE;

  safe_free((void *)res_p);
  safe_free((void *)res_m);

  return (TRUE);
} /* End of routine estimate_bAC */

#ifdef NOT_USED
static int estimate_dAC_LSvel(
    int iAC, double x_AC[], double **dAC, int numProcUnknowns, Comm_Ex *cx, MF_Args *mf_args) {

  static char *yo = "estimate_dAC";

  double p_save, dp_save;
  double fd_factor = FD_FACTOR;

  double inventory_p, inventory_m, xm, *x = mf_args->x;
  /*  int ielem, err=-99;
      int e_start = mf_args->exo->eb_ptr[0];
      int e_end   = mf_args->exo->eb_ptr[mf_args->exo->num_elem_blocks]; */

  int jAC = 0, chosen_vel = 0;

  for (jAC = 0; jAC < nAC; jAC++) {
    if (augc[jAC].Type == AC_LS_VEL && iAC == jAC) {

      chosen_vel = augc[jAC].DIR;

      p_save = x_AC[jAC];

      dp_save = fd_factor * p_save;
      dp_save = (fabs(dp_save) < fd_factor ? fd_factor : dp_save);
      /*  */

      x_AC[jAC] = p_save + dp_save;

      update_parameterAC(jAC, mf_args->x, mf_args->xdot, x_AC, cx, mf_args->exo, mf_args->dpi);

      inventory_p =
          find_LS_vel(mf_args->exo, mf_args->dpi, NULL, chosen_vel, mf_args->x, numProcUnknowns);
      inventory_p -= augc[iAC].CONSTV;

      /*  */
      x_AC[jAC] = p_save - dp_save;

      update_parameterAC(jAC, mf_args->x, mf_args->xdot, x_AC, cx, mf_args->exo, mf_args->dpi);

      inventory_m =
          find_LS_vel(mf_args->exo, mf_args->dpi, NULL, chosen_vel, mf_args->x, numProcUnknowns);
      inventory_m -= augc[iAC].CONSTV;

      /*  */
      xm = 0.5 / dp_save;
      /*  */
      dAC[iAC][jAC] = xm * (inventory_p - inventory_m);

      /*  */
      x_AC[jAC] = p_save;

      update_parameterAC(jAC, mf_args->x, mf_args->xdot, x_AC, cx, mf_args->exo, mf_args->dpi);

    } /* end if augc[jAC] == AC_LS_VEL  */
  }   /* end  for{jAC...  */
  return (TRUE);
} /* End of routine estimate_dAC */
#endif

static int estimate_dAC_ALC(
    int iAC, double x_AC[], double **dAC, int numProcUnknowns, Comm_Ex *cx, MF_Args *mf_args)
/*
 * This routine merely calculates d(gAC) / d(lambda) by centered
 * differencing and fills in the last column of dAC for arc
 * length continuation done as an additional AC.
 */
{
  double p_save, dp_save;
  double fd_factor = FD_FACTOR;

  double *res_p, *res_m, *x = mf_args->x, *xdot = mf_args->xdot;
  int jAC;
  int ibc, idf;
  int have_bAC = TRUE, have_cAC = TRUE, have_dAC = TRUE;

  double inventory;
  double LSvel_inventory;
#ifdef PARALLEL
  double global_inventory = 0.0;
  double global_LSvel_inventory = 0.0;
#endif

  /* Initialize perturbed AC residual arrays */
  asdv(&res_p, nAC);
  asdv(&res_m, nAC);

  /* Set perturbation to continuation parameter (in x_AC) */
  p_save = x_AC[iAC];
  dp_save = fd_factor * p_save;
  dp_save = (fabs(dp_save) < fd_factor ? fd_factor : dp_save);

  /* Apply positive perturbation */
  x_AC[iAC] = p_save + dp_save;
  update_parameterAC(iAC, x, xdot, x_AC, cx, mf_args->exo, mf_args->dpi);

  /* Evaluate perturbed residuals of all AC's except arc length */
  for (jAC = 0; jAC < iAC; jAC++) {
    switch (augc[jAC].Type) {
    case AC_FLUX:
      inventory = evaluate_flux(mf_args->exo, mf_args->dpi, augc[jAC].SSID, augc[jAC].MFID, NULL,
                                augc[jAC].MTID, augc[jAC].COMPID, NULL, FALSE, x, xdot, NULL,
                                *(mf_args->delta_t), *(mf_args->time), 0);
#ifdef PARALLEL
      if (Num_Proc > 1) {
        MPI_Allreduce(&inventory, &global_inventory, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        inventory = global_inventory;
      }
#endif
      res_p[jAC] = inventory - augc[jAC].CONSTV;
      break;

    case AC_FLUX_MAT:
      inventory = evaluate_flux(mf_args->exo, mf_args->dpi, augc[jAC].SSID, augc[jAC].MFID, NULL,
                                augc[jAC].MDID, augc[jAC].COMPID, NULL, FALSE, x, xdot, NULL,
                                *(mf_args->delta_t), *(mf_args->time), 0);
#ifdef PARALLEL
      if (Num_Proc > 1) {
        MPI_Allreduce(&inventory, &global_inventory, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        inventory = global_inventory;
      }
#endif
      res_p[jAC] = inventory - augc[jAC].CONSTV;
      break;

    case AC_VOLUME:
      inventory = augc[jAC].evol;
#ifdef PARALLEL
      if (Num_Proc > 1) {
        MPI_Allreduce(&inventory, &global_inventory, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        inventory = global_inventory;
      }
#endif
      res_p[jAC] = inventory - augc[jAC].CONSTV;
      break;

    case AC_LS_VEL:
      if (augc[iAC].d_lsvol_dx == 0)
        DPRINTF(stderr, "\n can't see vectors ");
      inventory = augc[jAC].lsvol;
      LSvel_inventory = augc[jAC].lsvel;
#ifdef PARALLEL
      if (Num_Proc > 1) {
        MPI_Allreduce(&inventory, &global_inventory, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&LSvel_inventory, &global_LSvel_inventory, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        inventory = global_inventory;
        LSvel_inventory = global_LSvel_inventory;
      }
#endif
      inventory = LSvel_inventory / inventory;
      res_p[jAC] = inventory - augc[jAC].CONSTV;
      break;

    case AC_LGRM:
      if (augc[jAC].MFID != VOLUME_FLUX) {
        GOMA_EH(GOMA_ERROR, "This Lagrange multiplier type is not yet supported!");
      } else {
        ibc = augc[jAC].BCID;
        idf = augc[jAC].DFID;
        inventory = evaluate_flux(mf_args->exo, mf_args->dpi, augc[jAC].SSID, augc[jAC].MFID, NULL,
                                  augc[jAC].MTID, augc[jAC].COMPID, NULL, FALSE, x, xdot, NULL,
                                  *(mf_args->delta_t), *(mf_args->time), 0);
        res_p[jAC] = inventory - (-BC_Types[ibc].BC_Data_Float[idf]);
      }
      break;

    case AC_USERBC:
    case AC_USERMAT:
      load_extra_unknownsAC(jAC, x_AC, cx, mf_args->exo, mf_args->dpi);
      user_aug_cond_residuals(jAC, x, xdot, *(mf_args->delta_t), *(mf_args->time), NULL, res_p,
                              &have_bAC, &have_cAC, &have_dAC, NULL, NULL, NULL, mf_args->exo,
                              mf_args->dpi, cx);
      break;

    case AC_ARC_LENGTH:
      break;

    default:
      GOMA_EH(-1, "Unknown AC type!");
    }
  }

  /* Apply negative perturbation */
  x_AC[iAC] = p_save - dp_save;
  update_parameterAC(iAC, x, xdot, x_AC, cx, mf_args->exo, mf_args->dpi);

  /* Evaluate perturbed residuals of all AC's except arc length */
  for (jAC = 0; jAC < iAC; jAC++) {
    switch (augc[jAC].Type) {
    case AC_FLUX:
      inventory = evaluate_flux(mf_args->exo, mf_args->dpi, augc[jAC].SSID, augc[jAC].MFID, NULL,
                                augc[jAC].MTID, augc[jAC].COMPID, NULL, FALSE, x, xdot, NULL,
                                *(mf_args->delta_t), *(mf_args->time), 0);
#ifdef PARALLEL
      if (Num_Proc > 1) {
        MPI_Allreduce(&inventory, &global_inventory, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        inventory = global_inventory;
      }
#endif
      res_m[jAC] = inventory - augc[jAC].CONSTV;
      break;

    case AC_FLUX_MAT:
      inventory = evaluate_flux(mf_args->exo, mf_args->dpi, augc[jAC].SSID, augc[jAC].MFID, NULL,
                                augc[jAC].MDID, augc[jAC].COMPID, NULL, FALSE, x, xdot, NULL,
                                *(mf_args->delta_t), *(mf_args->time), 0);
#ifdef PARALLEL
      if (Num_Proc > 1) {
        MPI_Allreduce(&inventory, &global_inventory, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        inventory = global_inventory;
      }
#endif
      res_m[jAC] = inventory - augc[jAC].CONSTV;
      break;

    case AC_VOLUME:
      inventory = augc[jAC].evol;
#ifdef PARALLEL
      if (Num_Proc > 1) {
        MPI_Allreduce(&inventory, &global_inventory, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        inventory = global_inventory;
      }
#endif
      res_m[jAC] = inventory - augc[jAC].CONSTV;
      break;

    case AC_LS_VEL:
      if (augc[iAC].d_lsvol_dx == 0)
        DPRINTF(stderr, "\n can't see vectors ");
      inventory = augc[jAC].lsvol;
      LSvel_inventory = augc[jAC].lsvel;
#ifdef PARALLEL
      if (Num_Proc > 1) {
        MPI_Allreduce(&inventory, &global_inventory, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&LSvel_inventory, &global_LSvel_inventory, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        inventory = global_inventory;
        LSvel_inventory = global_LSvel_inventory;
      }
#endif
      inventory = LSvel_inventory / inventory;
      res_m[jAC] = inventory - augc[jAC].CONSTV;
      break;

    case AC_LGRM:
      if (augc[jAC].MFID != VOLUME_FLUX) {
        GOMA_EH(GOMA_ERROR, "This Lagrange multiplier type is not yet supported!");
      } else {
        ibc = augc[jAC].BCID;
        idf = augc[jAC].DFID;
        inventory = evaluate_flux(mf_args->exo, mf_args->dpi, augc[jAC].SSID, augc[jAC].MFID, NULL,
                                  augc[jAC].MTID, augc[jAC].COMPID, NULL, FALSE, x, xdot, NULL,
                                  *(mf_args->delta_t), *(mf_args->time), 0);
        res_m[jAC] = inventory - (-BC_Types[ibc].BC_Data_Float[idf]);
      }
      break;

    case AC_USERBC:
    case AC_USERMAT:
      load_extra_unknownsAC(jAC, x_AC, cx, mf_args->exo, mf_args->dpi);
      user_aug_cond_residuals(jAC, x, xdot, *(mf_args->delta_t), *(mf_args->time), NULL, res_m,
                              &have_bAC, &have_cAC, &have_dAC, NULL, NULL, NULL, mf_args->exo,
                              mf_args->dpi, cx);
      break;

    case AC_ARC_LENGTH:
      break;

    default:
      GOMA_EH(-1, "Unknown AC type!");
    }
  }

  /* Restore original continuation parameter */
  x_AC[iAC] = p_save;
  update_parameterAC(iAC, x, xdot, x_AC, cx, mf_args->exo, mf_args->dpi);

  /* Now fill in right column of dAC */
  for (jAC = 0; jAC < iAC; jAC++) {
    dAC[jAC][iAC] = (res_p[jAC] - res_m[jAC]) * 0.5 / dp_save;
  }

  /* clean up and return */
  safe_free((void *)res_p);
  safe_free((void *)res_m);
  return (TRUE);
} /* End of estimate_dAC_ALC */

/*
 * Compute standard lagrange multiplier constraints and
 * insert into appropriate vectors
 */

int std_lgr_cond(int iAC,
                 int nAC,
                 double x_AC[],
                 double **bAC,
                 double **cAC,
                 double **dAC,
                 double *gAC,
                 double *resid_vector,
                 double *scale,
                 int numProcUnknowns,
                 Comm_Ex *cx,
                 MF_Args *mf_args)

{
  double inventory;
  int i, jAC;
  double theta_scale = 1.0;
#ifdef PARALLEL
  double global_inventory = 0.0;
#endif

  if (upd->CoordinateSystem == CYLINDRICAL || upd->CoordinateSystem == SWIRLING) {
    theta_scale =
        2.0 *
        M_PIE; /* this scaling is to remove the azimuthal factor that is applied by evaluate_flux */
  }

  /*
   *  LOAD EXTRA UNKNOWNS INTO x_AC
   *    Get the value of the unknown corresponding to the augmented condition
   *    and load it into the vector of Augmented condition unknowns, x_AC[]
   */
  load_extra_unknownsAC(iAC, x_AC, cx, mf_args->exo, mf_args->dpi);

  if (augc[iAC].Type == AC_LGRM) {

    switch (augc[iAC].MFID) {

    case VOLUME_FLUX: {
      int inode, i_offset;
      VARIABLE_DESCRIPTION_STRUCT *vd;
      int ibc = augc[iAC].BCID, idf = augc[iAC].DFID;

      inventory = evaluate_flux(mf_args->exo, mf_args->dpi, augc[iAC].SSID, augc[iAC].MFID, NULL,
                                augc[iAC].MTID, augc[iAC].COMPID, NULL, FALSE, mf_args->x,
                                mf_args->xdot, cAC[iAC], *(mf_args->delta_t), *(mf_args->time), 0);

      gAC[iAC] = inventory - (-BC_Types[ibc].BC_Data_Float[idf]);

      for (inode = 0, i = 0; i < numProcUnknowns; i++) {
        int var, isDBC = -1;
        vd = Index_Solution_Inv(i, &inode, NULL, &i_offset, NULL, pg->imtrx);

        var = vd->Variable_Type;

        if (Nodes[inode]->DBC[pg->imtrx])
          isDBC = Nodes[inode]->DBC[pg->imtrx][i_offset];

        if ((var >= VELOCITY1 &&
             var <= VELOCITY3) && /* this stuff only gets added to velocity dofs */
            isDBC == -1)          /* if a Dirichlet condition exists skip  */
        {
          bAC[iAC][i] = cAC[iAC][i] / theta_scale; /* this is true for LM type of constraints... */
        }
      }
    } break;
    default:

      GOMA_EH(GOMA_ERROR, "Cannot find Lagrange Multiplier condition.");
    }
  } else if (augc[iAC].Type == AC_PF_CONSTRAINT) {
    inventory = pfd->Constraint_Integral;
#ifdef PARALLEL
    if (Num_Proc > 1) {
      MPI_Allreduce(&inventory, &global_inventory, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      inventory = global_inventory;
    }
#endif
    for (i = 0; i < numProcUnknowns; i++) {
      cAC[iAC][i] = pfd->jac_info->d_lm_pf[i];
      bAC[iAC][i] = pfd->jac_info->d_pf_lm[i];
    }

    gAC[iAC] = (inventory - augc[iAC].CONSTV);
    /*	  gAC[iAC] = (inventory ) ;  */
  }

  /*
   * Last of all, zero out the constraint sensitivity wrt the Lagrange multiplier
   *  since it is always zero.....
   */

  for (jAC = 0; jAC < nAC; jAC++)
    dAC[iAC][jAC] = 0.0;

  return (1);
}

int arc_length_status(
    struct con_struct *con, double equation, double delta_y, double reltol, double abstol)
/*
 * This routine is called when the arc length bordering algorithm
 * has been done within an augmenting condition outside of LOCA.
 * Here, the convergence status is checked and reported only.
 */
{
  struct general_info_struct *cgi = &(con->general_info);
  struct private_info_struct *cpi = &(con->private_info);
  double param_update, scaled_resid;

  /*
   * Check whether or not the arc-length equation and continuation param
   * are converged (after first arc length step).
   */

  if (cpi->step_num == 0)
    return (TRUE);
  param_update = fabs(delta_y) / (reltol * fabs(cgi->param) + abstol);
  scaled_resid = fabs(equation) / (reltol * fabs(cpi->arc_step) + abstol);

  if (cgi->printproc > 4) {
    printf("\n\tEXTERNAL Arc Length Continuation: Convergence Criteria\n");
    printf("\tVariable	 Scaled Update (<1)  Unscaled Update  New Value\n");
    printf("\t***********************************************************\n");
    printf("\tparameter	 %e	   %e	 %g\n", param_update, -delta_y, cgi->param - delta_y);
    printf("\tarc length	 %e	   %e	 %g\n", scaled_resid, equation,
           cpi->arc_step + equation);
    printf("\t***********************************************************\n");
  }
  cgi->param -= delta_y;

  if ((param_update < 1.0) && (scaled_resid < 1.0))
    return (TRUE);
  else
    return (FALSE);
}

int periodic_bc_cond(int iAC,
                     int nAC,
                     double x_AC[],
                     double **bAC,
                     double **cAC,
                     double **dAC,
                     double *gAC,
                     double *resid_vector,
                     double *scale,
                     int numProcUnknowns,
                     Comm_Ex *cx,
                     MF_Args *mf_args)

{
  int var1, var2, jAC;

  /* DRN: probably a dumb idea, but I hate to just make the list of augmenting vars longer
   *      so I am using the fluid_eb and solid_eb to hold the equation indices that are tied
   * together by this augmenting condition
   */
  var1 = augc[iAC].fluid_eb;
  var2 = augc[iAC].solid_eb;
  if (var1 >= 0 && var2 >= 0) {
    /*fprintf(stderr,"v1=%g,v2=%g,r1=%g,r2=%g,lm=%g\n",
    mf_args->x[var1],mf_args->x[var2],resid_vector[var1],resid_vector[var2],augc[iAC].lm_value);*/

    resid_vector[var1] += augc[iAC].lm_value / scale[var1];
    resid_vector[var2] -= augc[iAC].lm_value / scale[var2];

    bAC[iAC][var1] = 1.;
    bAC[iAC][var2] = -1.;

    gAC[iAC] = mf_args->x[var1] - mf_args->x[var2];
    cAC[iAC][var1] = 1.;
    cAC[iAC][var2] = -1.;
    for (jAC = 0; jAC < nAC; jAC++)
      dAC[iAC][jAC] = 0.0;
  } else {
    gAC[iAC] = augc[iAC].lm_value;
    for (jAC = 0; jAC < nAC; jAC++)
      dAC[iAC][jAC] = 0.0;
    dAC[iAC][iAC] = 1.;
  }

  return (1);
}

int create_overlap_acs(Exo_DB *exo, int iAC)
/*
 * This routine is invoked when an OVERLAP type AC input card is specified.
 * Here, the original copy of that AC structure is replaced with a set of
 * new AC structures for each Lagrange multiplier constraint quadrature pt.
 * A call to realloc_struct_1 is made to accommodate the new records.
 *
 * INPUTS:
 *
 *      exo:            ptr to ExodusII database
 *      iAC:            number of Overlap AC (from card)
 *
 */
{
  int old_nAC, new_ACs = -1, new_nAC;
  int ac_count, i, j, k, nqp, si = 0, ss = -1;
  int sb, fb, lb;
  int e, s;
  int found = FALSE;
  int dim = pd->Num_Dim;
  int sse_first[MAX_NGV];
  int inode, ac_lm = 0;
  int num_elem = 0;
  FILE *ac_id;

  /* Save old nAC value */
  old_nAC = nAC;
  ac_count = old_nAC - 1;

  /* Read inputs from overlap AC card */
  sb = augc[iAC].solid_eb;
  fb = augc[iAC].fluid_eb;
  lb = augc[iAC].lm_eb;

  /* Check for valid Lagrange multiplier specification */
  if (lb == sb) {
    ac_lm = 1;
  } else if (lb == fb) {
    ac_lm = 2;
  } else {
    GOMA_EH(GOMA_ERROR,
            "Augmenting Condition Lagrange multiplier must live on fluid or solid block!");
  }

  /* At present, only 2D problems are supported */
  if (dim != 2)
    GOMA_EH(GOMA_ERROR, "Overlap AC algorithm is only for 2D problems!");

  /* Determine number of quadrature points per side (on solid block) */
  /* For now, only P0 interpolation of Lagrange equation will be supported */
  nqp = 1;

  if (ac_lm == 1) {
    /*
     *  Number of new AC's = (# of element sides in side set)
     *                     * (# of quadrature points per side)
     *                     * (# of problem dimensions)
     */

    /* Find the solid side set index */
    si = augc[iAC].SSID;
    for (i = 0; i < exo->num_side_sets; i++) {
      if (!found && exo->ss_id[i] == si) {
        found = TRUE;
        ss = i;
      }
    }

    /* Here, just count element sides along the side set */
    for (i = 0; i < MAX_NGV; i++) {
      sse_first[i] = TRUE;
    }
    new_ACs = exo->ss_num_sides[ss] * nqp * dim;
  } else if (ac_lm == 2) {
    /* hmmm, hardcoded to maximum */
    num_elem = (MAX_NGV - 5) / dim;
    new_ACs = num_elem * dim;
  }

  new_nAC = old_nAC + new_ACs - 1;

  /* See if MAX_NGV will allow this number of AC's */
  if (new_nAC > (MAX_NGV - 5)) {
    DPRINTF(stderr, "Must increase MAX_NGV to at least %d\n", new_nAC + 5);
    exit(-1);
  }

  /* Now, expand augc to make room for them */
  augc = (struct AC_Information *)realloc_struct_1(augc, struct AC_Information, new_nAC, old_nAC);

  if (ac_lm == 1) {
    /* Open ID output file for overlap AC data */
    ac_id = fopen("overlap_ac_id.dat", "w");

    /*
     * Now fill in data for each new AC constraint.
     * Numbering convention: Element or side, then local LM DOF, then dimension.
     */
    for (i = 0; i < exo->ss_num_sides[ss]; i++) {
      e = exo->ss_elem_list[exo->ss_elem_index[ss] + i];
      s = exo->ss_side_list[exo->ss_elem_index[ss] + i];
      /* These are to identify mid-node or last node on this side */
      k = exo->ss_node_side_index[ss][i + 1];
      inode = exo->ss_node_list[ss][k - 1];

      /*
       * sse_first array was set to identify the first side of each element
       * which lies on the boundary side set. This is to avoid duplication
       * when Lagrange multiplier equation is defined.
       */
      if (sse_first[i]) {
        for (j = 0; j < nqp; j++) {
          for (k = 0; k < pd->Num_Dim; k++) {
            augc[ac_count].Type = AC_OVERLAP;
            augc[ac_count].BCID = -99;
            augc[ac_count].VOLID = -99;
            augc[ac_count].MPID = -99;
            augc[ac_count].MFID = -99;
            augc[ac_count].MTID = -99;
            augc[ac_count].MDID = -99;
            augc[ac_count].SSID = si;
            augc[ac_count].fluid_eb = fb - 1;
            augc[ac_count].solid_eb = sb - 1;
            augc[ac_count].lm_eb = lb - 1;
            augc[ac_count].lm_elem = e;
            augc[ac_count].lm_side = s;
            augc[ac_count].lm_dim = k;
            augc[ac_count].lm_value = 0.0;
            augc[ac_count].lm_resid = 0.0;

            ac_count++;
          }

          /* Write AC ID data */
          fprintf(ac_id, "Solid element %5d side %d  <%g, %g> AC's:  %3d (X), %3d (Y)\n", e + 1, s,
                  Coor[0][inode], Coor[1][inode], ac_count - 1, ac_count);
        }
      }
    }
    fclose(ac_id);
  } else if (ac_lm == 2) {
    for (i = 0; i < num_elem; i++) {
      for (k = 0; k < pd->Num_Dim; k++) {
        augc[ac_count].Type = AC_OVERLAP;
        augc[ac_count].BCID = -99;
        augc[ac_count].VOLID = -99;
        augc[ac_count].MPID = -99;
        augc[ac_count].MFID = -99;
        augc[ac_count].MTID = -99;
        augc[ac_count].MDID = -99;
        augc[ac_count].SSID = si;
        augc[ac_count].fluid_eb = fb - 1;
        augc[ac_count].solid_eb = sb - 1;
        augc[ac_count].lm_eb = lb - 1;
        augc[ac_count].lm_elem = -1; /* no elem assigned yet */
        augc[ac_count].lm_side = -1; /* side has no meaning in this formulation */
        augc[ac_count].lm_dim = k;
        augc[ac_count].lm_value = 0.0;
        augc[ac_count].lm_resid = 0.0;

        ac_count++;
      }
    }
  }

  /* Update ALL AC's with new nAC value */
  for (i = 0; i < new_nAC; i++) {
    augc[i].nAC = new_nAC;
  }
  DPRINTF(stdout, "\tCreated %d new AC's for overlap algorithm.\n", new_ACs);

  /* Change nAC and return */
  if (ac_count == new_nAC) {
    nAC = new_nAC;
    return 0;
  } else
    return -1;
} /* END of function create_overlap_acs */

int assign_overlap_acs(double x[], Exo_DB *exo)
/*
 *
 */
{
  int iAC, k;
  int fb, lb;
  int e;
  int e_start, e_end;
  int ibf = -1, ib;

  for (iAC = 0; iAC < nAC; iAC++) {
    if (augc[iAC].Type == AC_OVERLAP)
      break;
  }
  if (augc[iAC].Type != AC_OVERLAP) {
    GOMA_EH(GOMA_ERROR, "Cannot locate first overlap AC\n");
  }

  /* Read inputs from existing overlap AC */
  fb = augc[iAC].fluid_eb + 1;
  lb = augc[iAC].lm_eb + 1;

  /* only need to reassign if AC LM lives on fluid mesh */
  if (lb != fb)
    return 0;

  /* Locate fluid and solid blocks in exo structure */
  for (ib = 0; ib < exo->num_elem_blocks; ib++) {
    if (exo->eb_id[ib] == fb)
      ibf = ib;
  }

  /* Locate AC on every element in fluid element block that needs one */
  e_start = exo->eb_ptr[ibf];
  e_end = exo->eb_ptr[ibf + 1];

  /* now populate augc with elem number */
  for (e = e_start; e < e_end; e++) {
    /* Check for level set going through element (fill unknown sign change) */
    if (elem_overlaps_interface(e, x, exo, ls->Length_Scale)) {
      for (k = 0; k < pd->Num_Dim; k++) {
        if (iAC >= nAC || (augc[iAC].Type != AC_OVERLAP)) {
          /* Sometimes we get here cuz there is a defgrad.  In this
           * case, we want the time step to fail so that we can
           * continue and not abort.
           */
          return -1;
        }
        augc[iAC].lm_elem = e;

        iAC++;
      }
    }
  }
  printf("Overlap AC: %d of the available %d ACs are currently assigned\n", iAC, nAC);
  for (; iAC < nAC; iAC++) {
    if (augc[iAC].Type == AC_OVERLAP)
      augc[iAC].lm_elem = -1;
  }

  return (-2);
} /* END of function assign_overlap_acs */

int create_periodic_acs(Exo_DB *exo)
/*
 * This routine is invoked when an PERIODIC type AC input card is specified.
 * Here, new entries in the AC structure are created for each BC.
 * A call to realloc_struct_1 is made to accommodate the new records.
 *
 * INPUTS:
 *
 *      exo:            ptr to ExodusII database
 *      iAC:            number of Periodic AC (from card)
 *
 */
{
  int old_nAC, new_ACs, new_nAC, iAC;
  int ac_count, i, j, k;
  int ssid1 = 0, ssid2 = 0, ss1 = -1, ss2 = -1, var, ie;
  int found = FALSE;
  int num_nodes_on_side, max_sides, max_nodes;
  int *node_list1, *node_list2, *node_list2match;
  int count1, count2, beg, end, n, node_num;
  int dim_count, dim = pd->Num_Dim;
  double scale;

  /* Save old nAC value */
  old_nAC = nAC;

  for (iAC = 0; iAC < old_nAC; iAC++) {
    if (augc[iAC].Type == AC_PERIODIC) {

      /* Read inputs from overlap AC card */
      var = augc[iAC].VAR;
      ssid1 = augc[iAC].SSID;
      ssid2 = augc[iAC].SSID2;

      /* set vars to -1 for this entry */
      augc[iAC].fluid_eb = -1;
      augc[iAC].solid_eb = -1;

      found = FALSE;
      for (i = 0; i < exo->num_side_sets; i++) {
        if (!found && exo->ss_id[i] == ssid1) {
          found = TRUE;
          ss1 = i;
        }
      }
      if (!found)
        GOMA_EH(GOMA_ERROR, "SSID1 not found in mesh for periodic bc");

      found = FALSE;
      for (i = 0; i < exo->num_side_sets; i++) {
        if (!found && exo->ss_id[i] == ssid2) {
          found = TRUE;
          ss2 = i;
        }
      }
      if (!found)
        GOMA_EH(GOMA_ERROR, "SSID2 not found in mesh for periodic bc");

      /* require SS's match */
      if (exo->ss_num_sides[ss1] != exo->ss_num_sides[ss2]) {
        GOMA_EH(GOMA_ERROR, "Periodic bc's only supported for matching side sets!");
      }

      num_nodes_on_side = (exo->ss_node_side_index[ss1][1] - exo->ss_node_side_index[ss1][0]);
      max_sides = exo->ss_num_sides[ss1];
      max_nodes = max_sides * num_nodes_on_side;

      /* allocate temp array to store list of nodes */
      node_list1 = (int *)smalloc(max_nodes * sizeof(int));
      node_list2 = (int *)smalloc(max_nodes * sizeof(int));

      count1 = 0;
      beg = exo->ss_node_side_index[ss1][0];
      end = exo->ss_node_side_index[ss1][exo->ss_num_sides[ss1]];
      for (n = beg; n < end; n++) {
        node_num = exo->ss_node_list[ss1][n];
        found = FALSE;
        for (i = 0; i < count1 && !found; i++) {
          if (node_num == node_list1[i])
            found = TRUE;
        }
        if (!found) {
          node_list1[count1++] = node_num;
        }
      }

      count2 = 0;
      beg = exo->ss_node_side_index[ss2][0];
      end = exo->ss_node_side_index[ss2][exo->ss_num_sides[ss2]];
      for (n = beg; n < end; n++) {
        node_num = exo->ss_node_list[ss2][n];
        found = FALSE;
        for (i = 0; i < count2 && !found; i++) {
          if (node_num == node_list2[i])
            found = TRUE;
        }
        if (!found) {
          node_list2[count2++] = node_num;
        }
      }

      if (count1 != count2) {
        GOMA_EH(GOMA_ERROR, "Periodic bc's only supported for matching side sets!");
      }

      /* now look to match nodes between lists */
      node_list2match = (int *)smalloc(max_nodes * sizeof(int));

      /* to make comparisons, we need a length scale
       * crude, but we'll get it from the distance between two points in list
       */
      scale = 0.;
      for (k = 0; k < dim; k++)
        scale += (Coor[k][node_list1[1]] - Coor[k][node_list1[0]]) *
                 (Coor[k][node_list1[1]] - Coor[k][node_list1[0]]);
      scale = sqrt(scale);

      for (i = 0; i < count1; i++) {
        found = FALSE;
        for (j = 0; j < count1 && !found; j++) {
          dim_count = 0;
          for (k = 0; k < dim; k++) {
            if (fabs(Coor[k][node_list1[i]] - Coor[k][node_list2[j]]) < 1.e-5 * scale)
              dim_count++;
          }
          if (dim_count == dim - 1) {
            found = TRUE;
            node_list2match[i] = node_list2[j];
          }
        }
        if (!found) {
          GOMA_EH(GOMA_ERROR, "Periodic bc's only supported for matching side sets!");
        }
      }

      /* reuse count2 to keep track of how many nodes are there that have the
         desired variable defined
       */
      count2 = 0;
      for (i = 0; i < count1; i++) {
        ie = Index_Solution(node_list1[i], var, 0, 0, -1, pg->imtrx);
        if (ie >= 0)
          count2++;
      }

      new_ACs = count2;
      new_nAC = nAC + new_ACs;

      /* See if MAX_NGV will allow this number of AC's */
      if (new_nAC > (MAX_NGV)) {
        DPRINTF(stderr, "Must increase MAX_NGV to at least %d\n", new_nAC);
        exit(-1);
      }

      /* Now, expand augc to make room for them */
      augc = (struct AC_Information *)realloc_struct_1(augc, struct AC_Information, new_nAC, nAC);

      /*
       * Now fill in data for each new AC constraint.
       */
      ac_count = nAC;
      for (i = 0; (i < count1) && (ac_count < new_nAC); i++) {
        augc[ac_count].Type = AC_PERIODIC;
        augc[ac_count].VAR = var;
        augc[ac_count].SSID = ssid1;
        augc[ac_count].SSID2 = ssid2;
        augc[ac_count].lm_value = 0.;
        /* DRN: probably a dumb idea, but I hate to just make the list of augmenting vars longer
         *      so I am using the fluid_eb and solid_eb to hold the equation indices that are tied
         * together by this augmenting condition
         */
        ie = Index_Solution(node_list1[i], var, 0, 0, -1, pg->imtrx);
        if (ie == -1)
          continue;
        augc[ac_count].fluid_eb = ie;
        ie = Index_Solution(node_list2match[i], var, 0, 0, -1, pg->imtrx);
        augc[ac_count].solid_eb = ie;
        /*
        fprintf(stderr,"ac_count=%d,node1=%d,node2=%d,ie1=%d,ie2=%d\n",
          ac_count,node_list1[i],node_list2match[i],augc[ac_count].fluid_eb,augc[ac_count].solid_eb);
        */

        ac_count++;
      }

      /* Update ALL AC's with new nAC value */
      for (i = 0; i < new_nAC; i++) {
        augc[i].nAC = new_nAC;
      }
      nAC = new_nAC;
      DPRINTF(stderr, "\tCreated %d new AC's for periodic bc.\n", new_ACs);

      safe_free(node_list1);
      safe_free(node_list2);
      safe_free(node_list2match);
    }
  }

  return 0;
} /* END of function create_overlap_acs */

static void RETN_COORD(int node, double *x, double *y, double *z) {
  if (node >= 0) {
    *x = Coor[0][node];
    *y = Coor[1][node];
    if (Num_Dim > 2)
      *z = Coor[2][node];
    else
      *z = 0.0;
  } else {
    *x = *y = *z = 0.0;
  }
}

#ifdef NOT_USED
static double fBal_Special(struct AC_Information *augc, double *cAC_iAC, double *x, Exo_DB *exo) {
  return 0.0;
}
#endif

//! This function return the real position of a single node in the mesh.
/*!
 *  This is used by the AC_POSITION augmenting condition to provide an
 *  additional equation to a system.
 */
double getPositionAC(struct AC_Information *augc, double *cAC_iAC, double *x, Exo_DB *exo) {
  double posNode[3];
  int i, ins, inode;
  NODE_INFO_STRUCT *node_ptr;
  double pos1D = 0.0;

  if (augc->COMPID != 0 && augc->COMPID != 1) {
    GOMA_EH(-1, "!augc->COMPID != 0,1 - don't know what to do?");
  }

  // Node set index of 1 node NS containing the position node
  int NSIndexPosition = augc->MTID;

  // Get the global node number from the structure
  int globalNodeNum = -1;

  for (ins = 0; ins < exo->num_node_sets; ins++) {
    /* Check for a match between the ID of the current node set and the node set
       ID specified by NSIndexPosition */
    if (exo->ns_id[ins] == NSIndexPosition) {
      if (exo->ns_num_nodes[ins] == 0 && Num_Proc > 1) {
        // we probably aren't this ns owner
        break;
      } else if (exo->ns_num_nodes[ins] != 1) {
        printf("NS %d, nodes should be equal to one %d\n", exo->ns_id[ins], exo->ns_num_nodes[ins]);
        exit(-1);
      }
      /* Get the 0th local node in the current node set */
      inode = exo->ns_node_list[exo->ns_node_index[ins] + 0];
      node_ptr = Nodes[inode];
      globalNodeNum = node_ptr->Global_Node_Num;
    }
  }

  // Decide whether this processor owns that node and get local node number
  int node = -1;
  for (i = 0; i < Num_Internal_Nodes + Num_Border_Nodes; i++) {
    node_ptr = Nodes[i];
    if (globalNodeNum == node_ptr->Global_Node_Num) {
      node = i;
      break;
    }
  }

  // get the real position of the node

  if (node >= 0) {
    //   First find the static position
    RETN_COORD(node, &posNode[0], &posNode[1], &posNode[2]);
    // Now add in the displacements
    //            ( we assume here no complications due to subparametric mapping)
    i = Index_Solution(node, MESH_DISPLACEMENT1, 0, 0, -1, pg->imtrx);
    if (i != -1) {
      posNode[0] += x[i];
      posNode[1] += x[i + 1];
      if (exo->num_dim > 2) {
        posNode[2] += x[i + 2];
      }
    }

    // Process a 1D idea of what we need from the position
    int coordDir = augc->VOLID;
    if (coordDir < 0 || coordDir > 2) {
      GOMA_EH(GOMA_ERROR, "shouldn't be here");
    }

    if (augc->COMPID == 1) {
      // Special section for buoyant cylinder

    } else {
      // for some reason GCC 12 will throw a warning despite the if above,
      // if that if (coorDir < 0 || coordDir > 2) is placed here in addition
      // to the one above it will not throw a warning for some reason
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
      pos1D = posNode[coordDir];
      // restore Warray-bounds
#pragma GCC diagnostic pop

      // Figure out cAC -> there are cases where the position may be an unknown
      // in the problem.
      if (i != -1) {
        cAC_iAC[i + coordDir] = 1.0;
      }
      augc->evol = pos1D;
      // Force balance (or zero)
      augc->lm_resid = pos1D - augc->CONSTV;
    }
  }
  return pos1D;
}

//! This function return the real position of a single node in the mesh.
/*!
 *  This is used by the AC_ANGLE augmenting condition to provide an
 *  additional equation to a system.
 */
double getAngleAC(struct AC_Information *augc, double *cAC_iAC, double *x, Exo_DB *exo) {
  double ordinate = 0., n1[DIM], n2[DIM], n_dot_n, n1_mag, n2_mag, xi[DIM];
  double n1_dx[DIM][DIM][MDE], n2_dx[DIM][DIM][MDE], d_ord_dx[4][DIM][DIM][MDE];
  int i, ins, inode, mat_num, face, node2, err, ielem, num_nodes_on_side;
  int p, q, k, ielem_dim;
  int local_side[2], side_nodes[3]; /* Assume quad has no more than 3 per side. */
  int elem_list[4], local_node[4], elem_ct;
  NODE_INFO_STRUCT *node_ptr;

  if (augc->COMPID != 0 && augc->COMPID != 1) {
    GOMA_EH(-1, "!augc->COMPID != 0,1 - don't know what to do?");
  }

  // Node set index of 1 node NS containing the position node
  int NSIndexPosition = augc->MTID;

  // Get the global node number from the structure
  int globalNodeNum = -1;

  for (ins = 0; ins < exo->num_node_sets; ins++) {
    /* Check for a match between the ID of the current node set and the node set
       ID specified by NSIndexPosition */
    if (exo->ns_id[ins] == NSIndexPosition) {
      if (exo->ns_num_nodes[ins] == 0 && Num_Proc > 1) {
        // we probably aren't this ns owner
        break;
      } else if (exo->ns_num_nodes[ins] != 1) {
        inode = exo->ns_node_list[exo->ns_node_index[ins] + 0];
        for (i = 1; i < exo->ns_num_nodes[ins]; i++) {
          if (exo->ns_node_list[exo->ns_node_index[ins] + i] != inode) {
            printf("NS %d, non-unique nodes in list of %d\n", exo->ns_id[ins],
                   exo->ns_num_nodes[ins]);
            printf("Proc %d Nodes %d %d\n", ProcID, inode,
                   exo->ns_node_list[exo->ns_node_index[ins] + i]);
          }
        }
      }
      /* Get the 0th local node in the current node set */
      inode = exo->ns_node_list[exo->ns_node_index[ins] + 0];
      node_ptr = Nodes[inode];
      globalNodeNum = node_ptr->Global_Node_Num;
    }
  }

  // Decide whether this processor owns that node and get local node number
  int node = -1;
  for (i = 0; i < Num_Internal_Nodes + Num_Border_Nodes; i++) {
    node_ptr = Nodes[i];
    if (globalNodeNum == node_ptr->Global_Node_Num) {
      node = i;
      break;
    }
  }

  if (node >= 0) {
    elem_list[0] = elem_list[1] = elem_list[2] = elem_list[3] = -1;
    local_node[0] = local_node[1] = local_node[2] = local_node[3] = -1;
    if (!exo->node_elem_conn_exists) {
      GOMA_EH(-1, "Cannot compute angle without node_elem_conn.");
    }

    elem_list[0] = exo->node_elem_list[exo->node_elem_pntr[node]];

    /*
     * Find out where this node appears in the elements local
     * node ordering scheme...
     */

    local_node[0] = in_list(node, exo->elem_node_pntr[elem_list[0]],
                            exo->elem_node_pntr[elem_list[0] + 1], exo->elem_node_list);

    GOMA_EH(local_node[0], "Can not find node in elem node connectivity!?! ");
    local_node[0] -= exo->elem_node_pntr[elem_list[0]];
    /* check for neighbors*/

    mat_num = map_mat_index(augc->VOLID);
    if (mat_num == find_mat_number(elem_list[0], exo)) {
      elem_ct = 1;
    } else {
      GOMA_WH(-1, "block id doesn't match first element");
    }
    load_ei(elem_list[0], exo, 0, pg->imtrx);
    for (face = 0; face < ei[pg->imtrx]->num_sides; face++) {
      ielem = exo->elem_elem_list[exo->elem_elem_pntr[elem_list[0]] + face];
      if (ielem != -1) {
        node2 = in_list(node, exo->elem_node_pntr[ielem], exo->elem_node_pntr[ielem + 1],
                        exo->elem_node_list);
        if (node2 != -1 && (mat_num == find_mat_number(ielem, exo))) {
          elem_list[elem_ct] = ielem;
          local_node[elem_ct] = node2;
          local_node[elem_ct] -= exo->elem_node_pntr[ielem];
          elem_ct++;
        }
      }
    }
    memset(d_ord_dx, 0.0, sizeof(dbl) * 4 * DIM * DIM * MDE);

    for (ielem = 0; ielem < elem_ct; ielem++) {
      if (local_node[ielem] < 0 || local_node[ielem] > 3) {
        GOMA_EH(-1, "Node out of bounds.");
      }

      /*
       * Now, determine the local name of the sides adjacent to this
       * node...this works for the exo patran convention for quads...
       *
       * Again, local_node and local_side are zero based...
       */

      local_side[0] = (local_node[ielem] + 3) % 4;
      local_side[1] = local_node[ielem];

      /*
       * With the side names, we can find the normal vector.
       * Again, assume the sides live on the same element.
       */
      load_ei(elem_list[ielem], exo, 0, pg->imtrx);

      /*
       * We abuse the argument list under the conditions that
       * we're going to do read-only operations and that
       * we're not interested in old time steps, time derivatives
       * etc.
       */
      err = load_elem_dofptr(elem_list[ielem], exo, x, x, x, x, 0);
      GOMA_EH(err, "load_elem_dofptr");

      err = bf_mp_init(pd);
      GOMA_EH(err, "bf_mp_init");

      /*
       * What are the local coordinates of the nodes in a quadrilateral?
       */

      find_nodal_stu(local_node[ielem], ei[pg->imtrx]->ielem_type, &xi[0], &xi[1], &xi[2]);

      err = load_basis_functions(xi, bfd);
      GOMA_EH(err, "problem from load_basis_functions");

      err = beer_belly();
      GOMA_EH(err, "beer_belly");

      err = load_fv();
      GOMA_EH(err, "load_fv");

      /* First, one side... */

      get_side_info(ei[pg->imtrx]->ielem_type, local_side[0] + 1, &num_nodes_on_side, side_nodes);

      surface_determinant_and_normal(elem_list[ielem], exo->elem_node_pntr[elem_list[ielem]],
                                     ei[pg->imtrx]->num_local_nodes, ei[pg->imtrx]->ielem_dim - 1,
                                     local_side[0] + 1, num_nodes_on_side, side_nodes);

      n1[0] = fv->snormal[0];
      n1[1] = fv->snormal[1];
      ielem_dim = ei[pg->imtrx]->ielem_dim;
      for (p = 0; p < ielem_dim; p++) {
        for (q = 0; q < ielem_dim; q++) {
          for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
            n1_dx[p][q][k] = fv->dsnormal_dx[p][q][k];
          }
        }
      }

      /* Second, the adjacent side of the quad... */

      get_side_info(ei[pg->imtrx]->ielem_type, local_side[1] + 1, &num_nodes_on_side, side_nodes);

      surface_determinant_and_normal(elem_list[ielem], exo->elem_node_pntr[elem_list[ielem]],
                                     ei[pg->imtrx]->num_local_nodes, ei[pg->imtrx]->ielem_dim - 1,
                                     local_side[1] + 1, num_nodes_on_side, side_nodes);

      n2[0] = fv->snormal[0];
      n2[1] = fv->snormal[1];
      for (p = 0; p < ielem_dim; p++) {
        for (q = 0; q < ielem_dim; q++) {
          for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
            n2_dx[p][q][k] = fv->dsnormal_dx[p][q][k];
          }
        }
      }

      /* cos (theta) = n1.n2 / ||n1|| ||n2|| */

      n_dot_n = n1_mag = n2_mag = 0.;
      for (p = 0; p < ielem_dim; p++) {
        n_dot_n += n1[p] * n2[p];
        n1_mag += n1[p] * n1[p];
        n2_mag += n2[p] * n2[p];
      }
      ordinate += M_PI - acos(n_dot_n / sqrt(n1_mag * n2_mag));
      for (p = 0; p < ielem_dim; p++) {
        for (q = 0; q < ielem_dim; q++) {
          for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
            d_ord_dx[ielem][p][q][k] +=
                (n1[p] * n2_dx[p][q][k] + n2[p] * n1_dx[p][q][k]) / sqrt(1. - SQUARE(n_dot_n));
          }
        }
      }

    } /*ielem loop    */

    /* For nonlocal element information, we do a
     * direct injection into a through Jac_BC.
     */

    int je;
    i = Index_Solution(node, MESH_DISPLACEMENT1, 0, 0, mat_num, pg->imtrx);
    for (ielem = 0; ielem < elem_ct; ielem++) {
      load_ei(elem_list[ielem], exo, 0, pg->imtrx);

      for (p = 0; p < ielem_dim; p++) {
        for (q = 0; q < ielem_dim; q++) {
          for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
            je = ei[pg->imtrx]->gun_list[MESH_DISPLACEMENT1][k];
            GOMA_EH(je, "Bad var index.");
            cAC_iAC[je] += d_ord_dx[ielem][p][q][k];
          }
        }
      }
    }

    // Figure out cAC -> there are cases where the position may be an unknown
    // in the problem.
    if (i != -1) {
      augc->evol = ordinate;
      augc->lm_resid = ordinate - augc->CONSTV;
    }
  }
  return ordinate;
}
/*****************************************************************************/
/* End of file mm_augc_util.c */
/*****************************************************************************/
