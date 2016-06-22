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
 

#ifdef USE_RCSID
static const char rcs_id[] = "$Id: ac_update_parameter.c,v 5.3 2008-12-19 22:54:25 rbsecor Exp $";
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "std.h"

#include "exo_struct.h"
#include "dpi.h"

#include "rf_allo.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_io_const.h"
#include "rf_io_structs.h"
#include "rf_io.h"
#include "rf_mp.h"
#include "rf_solver_const.h"
#include "rf_solver.h"

#include "rf_masks.h"

#include "el_geom.h"
#include "rf_vars_const.h"
#include "mm_mp_const.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"

#include "mm_mp.h"
#include "mm_mp_structs.h"

#include "mm_post_def.h"
#include "mm_post_proc.h"
#include "mm_eh.h"

#include "exo_struct.h"		/* defn of Exo_DB */
#include "dp_types.h"

#include "sl_util.h"		/* defines sl_init() */

#define _RF_UPDATE_PARAMETER_C
#include "user_continuation.h"
#include "goma.h"

/*

   UPDATE PARAMETER FOR CONTINUATION

   IAN GATES

   2/98 - 9/98

   Case lists consolidated - EDW 2/01

*/

void
update_parameterC(int iCC,       /* CONDITION NUMBER */
                  double lambda, /* PARAMETER VALUE */
		  double *x, 	 /* UNKNOWN VECTOR */
		  double *xdot,  /* UNKNOWN_DOT VECTOR */
		  double *x_AC,  /* x_AC VECTOR */
		  double delta_s,/* STEP */
		  Comm_Ex *cx,	 /* array of communications structures */
		  Exo_DB *exo,	 /* ptr to the finite element mesh database */
		  Dpi *dpi)	 /* distributed processing information */
{

  int ic, mn, mpr;
  int ibc, idf;
  //struct Boundary_Condition *BC_Type;

#ifdef DEBUG
  static const char yo[] = "update_parameterC";

  /*
   * 		BEGIN EXECUTION
   */

  fprintf(stderr, "%s() begins...\n", yo);
#endif

/* Calls from ac_conti.c (nCC=0) */
  if (nCC == 0) {
    if (cont->upType == 1) { /* BC */
      ibc = cont->upBCID;
      idf = cont->upDFID;
      update_BC_parameter(lambda, ibc, idf, cx, exo, dpi);
    }
    else if (cont->upType == 2) {
      mn = cont->upMTID;
      mpr = cont->upMPID;
      ic = cont->upMDID;
      update_MT_parameter(lambda, mn, mpr, ic, cx, exo, dpi);
    }
   else if (cont->upType == 3) {
      ibc = cont->upBCID;
      idf = cont->upDFID;
      update_AC_parameter(lambda, ibc, idf, cx, exo, dpi);
    }
    else if (cont->upType == 4) {
      mn = cont->upMTID;
      mpr = cont->upMPID;
      ic = cont->upMDID;
      update_UM_parameter(lambda, mn, mpr, ic, cx, exo, dpi);
    }
    else if (cont->upType == 5) {
      update_user_parameter(lambda, x, xdot, x_AC, cx, exo, dpi);
    }
 
    else EH(-1, "Bad continuation type!");
  }

/* Calls from LOCA */
 else {
    
  /*
   *    BC parameters
   */

  if (cpcc[iCC].Type == 1) {
    ibc = cpcc[iCC].BCID;
    idf = cpcc[iCC].DFID;
    update_BC_parameter(lambda, ibc, idf, cx, exo, dpi);
  }

  /*
   *    MT parameters
   */

  else if (cpcc[iCC].Type == 2) {
    mn = cpcc[iCC].MTID;
    mpr = cpcc[iCC].MPID;
    ic = cpcc[iCC].MDID;
    update_MT_parameter(lambda, mn, mpr, ic, cx, exo, dpi);
  }

  /*
   *    AC parameters
   */

  else if (cpcc[iCC].Type == 3) {
    ibc = cpcc[iCC].BCID;
    idf = cpcc[iCC].DFID;
    update_AC_parameter(lambda, ibc, idf, cx, exo, dpi);
  }

  /*
   *    UM parameters
   */

  else if (cpcc[iCC].Type == 4) {
    mn = cpcc[iCC].MTID;
    mpr = cpcc[iCC].MPID;
    ic = cpcc[iCC].MDID;
    update_UM_parameter(lambda, mn, mpr, ic, cx, exo, dpi);
  }

  /*
   *    UF parameters
   */

  else if (cpcc[iCC].Type == 5) {
    update_user_parameter(lambda, x, xdot, x_AC, cx, exo, dpi);
  }

  else {
    if (cont->upType != 6) EH(-1, "Bad continuation type!");
  }
 }

}

void
update_parameterTP(int iTC,       /* CONDITION NUMBER */
                   double lambda, /* PARAMETER VALUE */
		   double *x, 	  /* UNKNOWN VECTOR */
		   double *xdot,  /* UNKNOWN_DOT VECTOR */
		   double *x_AC,  /* x_AC VECTOR */
		   double delta_s,/* STEP */
		   Comm_Ex *cx,	  /* array of communications structures */
		   Exo_DB *exo,	  /* ptr to the finite element mesh database */
		   Dpi *dpi)	  /* distributed processing information */
{

  int ic, mn, mpr;
  int ibc, idf;
  //struct Boundary_Condition *BC_Type;

#ifdef DEBUG
  static const char yo[] = "update_parameterTP";

  /*
   * 		BEGIN EXECUTION
   */

  fprintf(stderr, "%s() begins...\n", yo);
#endif

  /*
   *    BC parameters
   */

  if (tpcc[iTC].Type == 1) {
    ibc = tpcc[iTC].BCID;
    idf = tpcc[iTC].DFID;
    update_BC_parameter(lambda, ibc, idf, cx, exo, dpi);
  }

  /*
   *    MT parameters
   */

  else if (tpcc[iTC].Type == 2) {
    mn = tpcc[iTC].MTID;
    mpr = tpcc[iTC].MPID;
    ic = tpcc[iTC].MDID;
    update_MT_parameter(lambda, mn, mpr, ic, cx, exo, dpi);
  }

  /*
   *    AC parameters
   */

  else if (tpcc[iTC].Type == 3) {
    ibc = tpcc[iTC].BCID;
    idf = tpcc[iTC].DFID;
    update_AC_parameter(lambda, ibc, idf, cx, exo, dpi);
  }

  /*
   *    UM parameters
   */

  else if (tpcc[iTC].Type == 4)
    {
      mn = tpcc[iTC].MTID;
      mpr = tpcc[iTC].MPID;
      ic = tpcc[iTC].MDID;
      update_UM_parameter(lambda, mn, mpr, ic, cx, exo, dpi);
    }

  /*
   *    UF parameters
   */

  else if (tpcc[iTC].Type == 5) {
    update_user_TP_parameter(lambda, x, xdot, x_AC, cx, exo, dpi);
  }

  else {
    if (loca_in->TPupType != 6) EH(-1, "Bad TP continuation type!");
  }
}

void
update_parameterHC(int iHC,      /* Hunting condition number */
                   double lambda, /* PARAMETER VALUE */
		   double *x, 	 /* UNKNOWN VECTOR */
		   double *xdot,  /* UNKNOWN_DOT VECTOR */
		   double *x_AC,  /* x_AC VECTOR */
		   double delta_s,/* STEP */
		   Comm_Ex *cx,	 /* array of communications structures */
		   Exo_DB *exo,	 /* ptr to the finite element mesh database */
		   Dpi *dpi)	 /* distributed processing information */
{

  int ic, mn, mpr;
  int ibc, idf;
  //struct Boundary_Condition *BC_Type;

#ifdef DEBUG
  static const char yo[] = "update_parameterHC";


  /*
   * 		BEGIN EXECUTION
   */

  fprintf(stderr, "%s() begins...\n", yo);
#endif

  /*
   *    BC parameters
   */

  if (hunt[iHC].Type == 1) {
    ibc = hunt[iHC].BCID;
    idf = hunt[iHC].DFID;
    update_BC_parameter(lambda, ibc, idf, cx, exo, dpi);
  }

  /*
   *    MT parameters
   */

  else if (hunt[iHC].Type == 2) {
    mn = hunt[iHC].MTID;
    mpr = hunt[iHC].MPID;
    ic = hunt[iHC].MDID;
    update_MT_parameter(lambda, mn, mpr, ic, cx, exo, dpi);
  }

  /*
   *    AC parameters
   */

  else if (hunt[iHC].Type == 3) {
    ibc = hunt[iHC].BCID;
    idf = hunt[iHC].DFID;
    update_AC_parameter(lambda, ibc, idf, cx, exo, dpi);
    }

  /*
   *    UM parameters
   */

  else if (hunt[iHC].Type == 4) {
    mn = hunt[iHC].MTID;
    mpr = hunt[iHC].MPID;
    ic = hunt[iHC].MDID;
    update_UM_parameter(lambda, mn, mpr, ic, cx, exo, dpi);
  }

  /*
   *    UF parameters
   */

  else if (hunt[iHC].Type == 5) {
    update_user_parameter(lambda, x, xdot, x_AC, cx, exo, dpi);
  }

  else EH(-1, "Bad HC continuation type!");

}

void
update_parameterS(double lambda, /* PARAMETER VALUE */
		  double *x, 	 /* UNKNOWN VECTOR */
		  double *xdot,  /* UNKNOWN_DOT VECTOR */
		  int sens_type,/*  type of sensitivity variable(BC or MT) */
		  int sens_id,		/*  variable id BCID or MT# */
		  int sens_flt,		/* data float id or matl prop id */
		  int sens_flt2,		/* data float id for UM */
		  Comm_Ex *cx,	 /* array of communications structures */
		  Exo_DB *exo,	 /* ptr to the finite element mesh database */
		  Dpi *dpi)	 /* distributed processing information */
{

  int ic, mn, mpr;
  int ibc, idf;
  //struct Boundary_Condition *BC_Type;

#ifdef DEBUG
  static const char yo[] = "update_parameterS";

  /*
   * 		BEGIN EXECUTION
   */

  fprintf(stderr, "%s() begins...\n", yo);
#endif

  /*
   *    BC parameters
   */

  if (sens_type == 1) {
    ibc = sens_id;
    idf = sens_flt;
    update_BC_parameter(lambda, ibc, idf, cx, exo, dpi);
  }

  /*
   *    MT parameters
   */

  else if (sens_type == 2) {
    mn = sens_id;
    mpr = sens_flt;
    ic = sens_flt;
    update_MT_parameter(lambda, mn, mpr, ic, cx, exo, dpi);
  }

  /*
   *    AC parameters
   */

  else if (sens_type == 3) {
    ibc = sens_id;
    idf = sens_flt;
    update_AC_parameter(lambda, ibc, idf, cx, exo, dpi);
  }

  /*
   *    UM parameters
   */

  else if (sens_type == 4) {
    mn = sens_id;
    mpr = sens_flt;
    ic = sens_flt2;
    update_UM_parameter(lambda, mn, mpr, ic, cx, exo, dpi);
  }

  /*
   *    UF parameters
   */

  else if (sens_type == 5) {
    update_user_parameter(lambda, x, xdot, NULL, cx, exo, dpi);
  }

  else EH(-1, "Bad sens. parameter type!");

}

  /*
   * Original parameter lists are now separated into the following
   * separate routines for BC and MT types - EDW 2/01.
   */

void
update_BC_parameter(double lambda, /* Parameter value */
                    int ibc,       /* Boundary condition index */
                    int idf,       /* Boundary condition float tag */
                    Comm_Ex *cx,   /* array of communications structures */
                    Exo_DB *exo,   /* ptr to the finite element mesh database */
                    Dpi *dpi)      /* distributed processing information */
{

  struct Boundary_Condition *BC_Type;

	switch (BC_Types[ibc].BC_Name)
	{
	case SPLINE_BC:
	case SPLINEX_BC:
	case SPLINEY_BC:
	case SPLINEZ_BC:
 	case SPLINE_RS_BC:
	case SPLINEX_RS_BC:
 	case SPLINEY_RS_BC:
 	case SPLINEZ_RS_BC:
	case FILLET_BC:
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
	case YUSER_BC:
		BC_Types[ibc].u_BC[idf] = lambda;
		break;
        case TABLE_BC:
        case TABLE_WICV_BC:
        case TABLE_WICS_BC:
                BC_Type = &BC_Types[ibc];
                BC_Type->table->f[idf] = lambda;
                break;
	case FRICTION_ACOUSTIC_BC:
		if(idf > 1)
			{ BC_Types[ibc].u_BC[idf-2] = lambda;}
		else
			{ BC_Types[ibc].BC_Data_Float[idf] = lambda;}
		break;
	default:
  	  	BC_Types[ibc].BC_Data_Float[idf] = lambda;
		break;
	}

} /* END of routine update_BC_parameter  */
/*****************************************************************************/

void
update_MT_parameter(double lambda, /* Parameter value */
                    int mn,        /* Material number index */
                    int mpr,       /* Material property index */
                    int ic,        /* Material property tag subindex */
                    Comm_Ex *cx,   /* array of communications structures */
                    Exo_DB *exo,   /* ptr to the finite element mesh database */
                    Dpi *dpi)      /* distributed processing information */
{

    switch (mpr) {

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
      
    case TAGC_ELECTRICAL_CONDUCTIVITY: 
      mp_glob[mn]->electrical_conductivity = lambda;
      break;
      
    case TAGC_PERMITTIVITY: 
      mp_glob[mn]->permittivity = lambda;
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

      /* CARREAU_WLF */

    case TAGC_WLFC2:
      gn_glob[mn]->wlfc2 = lambda;
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
      mp_glob[mn]->u_fdiffusivity[0][0]=lambda;
      break;
    case TAGC_FICKDIFF_Y:
      mp_glob[mn]->u_fdiffusivity[0][1]=lambda;
      break;

      /* these can be implemented for Ryan's Qtensor model as needed */

#if 0
    case TAGC_QTENSOR_EXTENSION_P:
      mp_glob[mn]->u_Qtensor_Extension_P[0][0]=lambda;
      break;
    case TAGC_QTENSOR_NCT:
      mp_glob[mn]->u_Qtensor_Nct[0][0]=lambda;
      break;
#endif
      
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
      ve_glob[mn][mpr-TAGC_TIME_CONST]->time_const = lambda;
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
      ve_glob[mn][mpr-TAGC_ALPHA]->alpha = lambda;
      break;

    case TAGC_PTT_XI:
      ve_glob[mn][cont->upMPID-TAGC_PTT_XI]->xi = lambda;
      break;

    case TAGC_PTT_EPS:
      ve_glob[mn][cont->upMPID-TAGC_PTT_EPS]->eps = lambda;
      break;

    case TAGC_SHIFT_FUNC:
    case TAGC_SHIFT_FUNC1:
      vn_glob[mn]->shift[mpr-TAGC_SHIFT_FUNC] = lambda;
      break;

      /* 
       * Constants used in the Elasticity Constitutive Equations
       */

    case TAGC_LAME_MU:
      elc_glob[mn]->lame_mu = lambda;
      break;

    case TAGC_LAME_MU_CONTACT_LINE_G0:
      *(elc_glob[mn]->u_mu+1) = lambda;
      break;

    case TAGC_LAME_MU_CONTACT_LINE_G1:
      *(elc_glob[mn]->u_mu+2) = lambda;
      break;

    case TAGC_LAME_MU_CONTACT_LINE_R0:
      *(elc_glob[mn]->u_mu+3) = lambda;
      break;

    case TAGC_LAME_LAMBDA:
      elc_glob[mn]->lame_lambda = lambda;
      break;

    case TAGC_BEND_STIFFNESS:
      elc_glob[mn]->bend_stiffness = lambda;
      break;

    case TAGC_CONV_LAG_VELX:
      elc_glob[mn]->v_mesh_sfs[0] = lambda;
      break;

    case TAGC_CONV_LAG_VELY:
      elc_glob[mn]->v_mesh_sfs[1] = lambda;
      break;

    case TAGC_CONV_LAG_VELZ:
      elc_glob[mn]->v_mesh_sfs[2] = lambda;
      break;

    case TAGC_CONV_LAG_ROTRATE:
      *(elc_glob[mn]->u_v_mesh_sfs) = lambda;
      break;

    case TAGC_CONV_LAG_ROT_X0:
      *(elc_glob[mn]->u_v_mesh_sfs+1) = lambda;
      break;
     
    case TAGC_CONV_LAG_ROT_Y0:
      *(elc_glob[mn]->u_v_mesh_sfs+2) = lambda;
      break;

    case TAGC_CONV_LAG_ROT_Z0:
      *(elc_glob[mn]->u_v_mesh_sfs+3) = lambda;
      break;

    case TAGC_RS_LAME_MU:
      elc_rs_glob[mn]->lame_mu = lambda;
      break;

    case TAGC_RS_LAME_LAMBDA:
      elc_rs_glob[mn]->lame_lambda = lambda;
      break;

    case TAGC_RS_CONV_LAG_VELX:
      elc_rs_glob[mn]->v_mesh_sfs[0] = lambda;
      break;

    case TAGC_RS_CONV_LAG_VELY:
      elc_rs_glob[mn]->v_mesh_sfs[1] = lambda;
      break;

    case TAGC_RS_CONV_LAG_VELZ:
      elc_rs_glob[mn]->v_mesh_sfs[2] = lambda;
      break;

    case TAGC_RS_CONV_LAG_ROTRATE:
      *(elc_rs_glob[mn]->u_v_mesh_sfs) = lambda;
      break;

    case TAGC_RS_CONV_LAG_ROT_X0:
      *(elc_rs_glob[mn]->u_v_mesh_sfs+1) = lambda;
      break;
     
    case TAGC_RS_CONV_LAG_ROT_Y0:
      *(elc_rs_glob[mn]->u_v_mesh_sfs+2) = lambda;
      break;

    case TAGC_RS_CONV_LAG_ROT_Z0:
      *(elc_rs_glob[mn]->u_v_mesh_sfs+3) = lambda;
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
      
    case TAGC_SHU_QFLOW:
      *(mp_glob[mn]->u_shell_user_par+3) = lambda;
      break;

    case TAGC_SHU_VWEB:
      *(mp_glob[mn]->u_shell_user_par+4) = lambda;
      break;

    case TAGC_SHU_ROLLRAD:
      *(mp_glob[mn]->u_shell_user_par+5) = lambda;
      break;

    case TAGC_SHU_X0:
      *(mp_glob[mn]->u_shell_user_par+6) = lambda;
      break;

    case TAGC_SHU_GAPN:
      *(mp_glob[mn]->u_shell_user_par+7) = lambda;
      break;

    case TAGC_SHU_UPS_XLOC:
      *(mp_glob[mn]->u_shell_user_par+8) = lambda;
      break;

    case TAGC_SHU_DNS_XLOC:
      *(mp_glob[mn]->u_shell_user_par+9) = lambda;
      break;

    default: 
      printf("\n\t Error: Invalid Material Property Tag %d\n", mpr);
      exit(0);
      break;
    }

}/* END of routine update_MT_parameter  */
/*******************************************************************************/

void
update_UM_parameter(double lambda, /* Parameter value */
                    int mn,        /* Material number index */
                    int mpr,       /* Material property index */
                    int ic,        /* User model float index */
                    Comm_Ex *cx,   /* array of communications structures */
                    Exo_DB *exo,   /* ptr to the finite element mesh database */                    Dpi *dpi)      /* distributed processing information */
/*
 * This function allows continuation in a user-defined material property
 * model parameter. It is invoked by setting the continuation type to "UM".
 * Many of the commonly-used models are already provided here.
 * More can be added as follows:
 *   1. Determine if the property of interest already has a user model option.
 *      If so, then the Material_Properties structure will contain at least
 *      the following members (see mm_mp_structs.h):
 *         dbl propname;
 *         dbl len_u_propname;
 *         dbl *u_propname;
 *         dbl d_propname[...]
 *         int PropnameModel;
 *      Add these if they are not already present.
 *   2. Determine if there is already a model read input for the property
 *      in mm_input_mp.c. If so, there will be a call to look_for_mat_prop()
 *      for this property; otherwise, you will need to put one in.
 *   3. Determine if there is a material property tag already defined.
 *      If so, it will appear in mm_mp_const.h as follows:
 *         #define TAGC_PROPNAME         xxxx
 *      Otherwise, add one here with a 4-digit integer xxxx which is not
 *      used by another TAGC_ define.
 *   4. Determine if there is already a function to read the model constants
 *      [i.e. usr_propname()] in user_mp.c. If not, create one using one of
 *      the existing usr_ functions as a template, and don't forget to put
 *      a prototype for it in user_mp.h.
 *
 * Now, to use this model (old or new) for continuation:
 *   1. Define the model in user_mp.c, noting which constants are which.
 *      As always, include all necessary sensitivities.
 *   2. Set the Propname Model card in the .mat file = USER and follow
 *      with the appropriate number of float constants. These will go into
 *      the mp->u_propname array in this order starting with index 0.
 *   3. In the Continuation Specifications section of the input deck,
 *      set these cards as follows:
 *         Continuation Type                 = UM
 *         Material id                       = i
 *         Material property tag             = j
 *         Material property tag subindex    = k
 *      Here, i and j identify the material number and specific property
 *      (the same as for type MT), and k is the (0-based) index of the 
 *      user constant to be used for continuation.
 */
{

  int float_length = -1;
  switch (mpr)
    {
      case TAGC_THERMAL_CONDUCTIVITY:
        float_length = mp_glob[mn]->len_u_thermal_conductivity;
        break;
      case TAGC_ELECTRICAL_CONDUCTIVITY:
        float_length = mp_glob[mn]->len_u_electrical_conductivity;
        break;
      case TAGC_VISCOSITY:
        float_length = mp_glob[mn]->len_u_viscosity;
        break;
      case TAGC_SURFACE_TENSION:
        float_length = mp_glob[mn]->len_u_surface_tension;
        break;
      case TAGC_HEAT_CAPACITY:
        float_length = mp_glob[mn]->len_u_heat_capacity;
        break;
      case TAGC_VOLUME_EXPANSION:
        float_length = mp_glob[mn]->len_u_Volume_Expansion;
        break;
      case TAGC_DENSITY:
        float_length = mp_glob[mn]->len_u_density;
        break;
      case TAGC_POROSITY:
        float_length = mp_glob[mn]->len_u_porosity;
        break;
      case TAGC_PERMEABILITY:
        float_length = mp_glob[mn]->len_u_permeability;
        break;
      case TAGC_REL_GAS_PERM:
        float_length = mp_glob[mn]->len_u_rel_gas_perm;
        break;
      case TAGC_REL_LIQ_PERM:
        float_length = mp_glob[mn]->len_u_rel_liq_perm;
        break;
      case TAGC_SATURATION:
        float_length = mp_glob[mn]->len_u_saturation;
        break;
      case TAGC_FLOWINGLIQUID_VISCOSITY:
        float_length = mp_glob[mn]->len_u_FlowingLiquid_viscosity;
        break;
      case TAGC_TAU_Y:
        float_length = gn_glob[mn]->len_u_tau_y;
        break;
      default:
        EH(-1, "No model available for that property!");
        break;
    }
        if (ic >= float_length)
          {
            EH(-1, "Float index larger than user parameter list length!");
          }


  switch (mpr)
    {
      case TAGC_THERMAL_CONDUCTIVITY:
            mp_glob[mn]->u_thermal_conductivity[ic] = lambda;
        break;
      case TAGC_ELECTRICAL_CONDUCTIVITY:
            mp_glob[mn]->u_electrical_conductivity[ic] = lambda;
        break;
      case TAGC_VISCOSITY:
            mp_glob[mn]->u_viscosity[ic] = lambda;
        break;
      case TAGC_SURFACE_TENSION:
            mp_glob[mn]->u_surface_tension[ic] = lambda;
        break;
      case TAGC_HEAT_CAPACITY:
            mp_glob[mn]->u_heat_capacity[ic] = lambda;
        break;
      case TAGC_VOLUME_EXPANSION:
            mp_glob[mn]->u_Volume_Expansion[ic] = lambda;
        break;
      case TAGC_DENSITY:
            mp_glob[mn]->u_density[ic] = lambda;
        break;
      case TAGC_POROSITY:
            mp_glob[mn]->u_porosity[ic] = lambda;
        break;
      case TAGC_PERMEABILITY:
            mp_glob[mn]->u_permeability[ic] = lambda;
        break;
      case TAGC_REL_GAS_PERM:
            mp_glob[mn]->u_rel_gas_perm[ic] = lambda;
        break;
      case TAGC_REL_LIQ_PERM:
            mp_glob[mn]->u_rel_liq_perm[ic] = lambda;
        break;
      case TAGC_SATURATION:
            mp_glob[mn]->u_saturation[ic] = lambda;
        break;
      case TAGC_TAU_Y:
            gn_glob[mn]->u_tau_y[ic] = lambda;
        break;
      case TAGC_FLOWINGLIQUID_VISCOSITY:
            mp_glob[mn]->u_FlowingLiquid_viscosity[ic] = lambda;
        break;
      default:
        EH(-1, "No model available for that property!");
        break;
    }
}/* END of routine update_UM_parameter  */
/*******************************************************************************/

void
update_AC_parameter(double lambda, /* Parameter value */
                    int ibc,       /* Boundary condition index */
                    int idf,       /* Boundary condition float tag */
                    Comm_Ex *cx,   /* array of communications structures */
                    Exo_DB *exo,   /* ptr to the finite element mesh database */
                    Dpi *dpi)      /* distributed processing information */
{

 	if ( idf == -1 || augc[ibc].len_AC == 0 )	
 		{
 		augc[ibc].CONSTV = lambda;
 		}
 	else
 		{
 		augc[ibc].DataFlt[idf] = lambda;
 		}
}/* END of routine update_AC_parameter  */
/*******************************************************************************/

void
retrieve_parameterS(double *lambda, /* PARAMETER VALUE */
		  double *x, 	 /* UNKNOWN VECTOR */
		  double *xdot,  /* UNKNOWN_DOT VECTOR */
		  int sens_type,/*  type of sensitivity variable(BC or MT) */
		  int sens_id,		/*  variable id BCID or MT# */
		  int sens_flt,		/* data float id or matl prop id */
		  int sens_flt2,	/* data float id for UM */
		  Comm_Ex *cx,	 /* array of communications structures */
		  Exo_DB *exo,	 /* ptr to the finite element mesh database */
		  Dpi *dpi)	 /* distributed processing information */
{

  int ic;
  int mn,mpr;
  int ibc, idf;
  //struct Boundary_Condition *BC_Type;
  
#ifdef DEBUG
  static const char yo[]="retrieve_parameterS";

  fprintf(stderr, "%s() begins...\n", yo);
#endif

  if (sens_type == 1) {

    ibc = sens_id;
    idf = sens_flt;

    retrieve_BC_parameter(lambda, ibc, idf, cx, exo, dpi);

  	}
  else
  if (sens_type == 2) {

    mn = sens_id;
    mpr = sens_flt;

      retrieve_MT_parameter(lambda, mn, mpr, cx, exo, dpi);
	}
  else
  if (sens_type == 3) {
 
     ibc = sens_id;
     idf = sens_flt;
 
      retrieve_AC_parameter(lambda, ibc, idf, cx, exo, dpi);
 	}
  else 
  if (sens_type == 4) {

    mn = sens_id;
    mpr = sens_flt;
    ic = sens_flt2;

      retrieve_UM_parameter(lambda, mn, mpr, ic, cx, exo, dpi);
  }

  else EH(-1, "Bad sens. parameter type!");
 
}/* END of routine retrieve_parameterS  */
/*****************************************************************************/
void
retrieve_BC_parameter(double *lambda, /* Parameter value */
                    int ibc,       /* Boundary condition index */
                    int idf,       /* Boundary condition float tag */
                    Comm_Ex *cx,   /* array of communications structures */
                    Exo_DB *exo,   /* ptr to the finite element mesh database */
                    Dpi *dpi)      /* distributed processing information */
{

  struct Boundary_Condition *BC_Type;

	switch (BC_Types[ibc].BC_Name)
	{
	case SPLINE_BC:
	case SPLINEX_BC:
	case SPLINEY_BC:
	case SPLINEZ_BC:
 	case SPLINE_RS_BC:
 	case SPLINEX_RS_BC:
 	case SPLINEY_RS_BC:
 	case SPLINEZ_RS_BC:
	case FILLET_BC:
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
	case YUSER_BC:
		*lambda = BC_Types[ibc].u_BC[idf];
		break;
        case TABLE_BC:
        case TABLE_WICV_BC:
        case TABLE_WICS_BC:
                BC_Type = &BC_Types[ibc];
                *lambda = BC_Type->table->f[idf];
                break;
	case FRICTION_ACOUSTIC_BC:
		if(idf > 1)
			{ *lambda = BC_Types[ibc].u_BC[idf-2];}
		else
			{ *lambda = BC_Types[ibc].BC_Data_Float[idf];}
		break;
	default:
  	  	*lambda = BC_Types[ibc].BC_Data_Float[idf];
		break;
	}

} /* END of routine retrieve_BC_parameter  */
/*****************************************************************************/

void
retrieve_MT_parameter(double *lambda, /* Parameter value */
                    int mn,        /* Material number index */
                    int mpr,       /* Material property index */
                    Comm_Ex *cx,   /* array of communications structures */
                    Exo_DB *exo,   /* ptr to the finite element mesh database */
                    Dpi *dpi)      /* distributed processing information */
{

    switch (mpr) {

      /* 
       * General Model Constants
       */

    case TAGC_THERMAL_CONDUCTIVITY: 
      *lambda = mp_glob[mn]->thermal_conductivity;
      break;
      
    case TAGC_ACOUSTIC_WAVENUMBER: 
      *lambda = mp_glob[mn]->wave_number;
      break;
      
    case TAGC_ACOUSTIC_IMPEDANCE: 
      *lambda = mp_glob[mn]->acoustic_impedance;
      break;
      
    case TAGC_ACOUSTIC_ABSORPTION: 
      *lambda = mp_glob[mn]->acoustic_absorption;
      break;
      
    case TAGC_REFRACTIVE_INDEX: 
      *lambda = mp_glob[mn]->refractive_index;
      break;
      
    case TAGC_LIGHT_ABSORPTION: 
      *lambda = mp_glob[mn]->light_absorption;
      break;
      
    case TAGC_ELECTRICAL_CONDUCTIVITY: 
      *lambda = mp_glob[mn]->electrical_conductivity;
      break;
      
    case TAGC_VISCOSITY: 
      *lambda = mp_glob[mn]->viscosity;
      break;         
      
    case TAGC_SURFACE_TENSION: 
      *lambda = mp_glob[mn]->surface_tension;
      break;
      
    case TAGC_HEAT_CAPACITY: 
      *lambda = mp_glob[mn]->heat_capacity;
      break;
      
    case TAGC_VOLUME_EXPANSION: 
      *lambda = mp_glob[mn]->Volume_Expansion;
      break;
      
    case TAGC_DENSITY: 
      *lambda = mp_glob[mn]->density;
      break;
      
    case TAGC_POROSITY: 
      *lambda = mp_glob[mn]->porosity;
      break;              
      
    case TAGC_PERMEABILITY: 
      *lambda = mp_glob[mn]->permeability;
      break;
      
    case TAGC_REL_GAS_PERM: 
      *lambda = mp_glob[mn]->rel_gas_perm;
      break;
      
    case TAGC_REL_LIQ_PERM: 
      *lambda = mp_glob[mn]->rel_liq_perm;
      break;
      
    case TAGC_SATURATION: 
      *lambda = mp_glob[mn]->saturation;
      break;
      
    case TAGC_MELTING_POINT_LIQUIDUS: 
      *lambda = mp_glob[mn]->melting_point_liquidus;
      break;
      
    case TAGC_MELTING_POINT_SOLIDUS: 
      *lambda = mp_glob[mn]->melting_point_solidus;
      break;
      
    case TAGC_FLOWINGLIQUID_VISCOSITY: 
      *lambda = mp_glob[mn]->FlowingLiquid_viscosity;
      break;
      
      /* 
       * Generalized Newtonian Models: 
       * Newtonian, Power Law, Carreau or Bingham(1,2,3)
       */
      
    case TAGC_MU0: 
      *lambda = gn_glob[mn]->mu0;
      break;
    case TAGC_NEXP: 
      *lambda = gn_glob[mn]->nexp;
      break;
    case TAGC_MUINF: 
      *lambda = gn_glob[mn]->muinf;
      break;
    case TAGC_LAM: 
      *lambda = gn_glob[mn]->lam;
      break;
    case TAGC_AEXP: 
      *lambda = gn_glob[mn]->aexp;
      break;
    case TAGC_ATEXP: 
      *lambda = gn_glob[mn]->atexp;
      break;

      /* CARREAU_WLF */

    case TAGC_WLFC2:
      *lambda = gn_glob[mn]->wlfc2;
      break;
      
      /* these are for the BINGHAM yielding material model */
      
    case TAGC_TAU_Y: 
      *lambda = gn_glob[mn]->tau_y;
      break;
    case TAGC_FEXP: 
      *lambda = gn_glob[mn]->fexp;
      break;
      
      /* these are for SUSPENSION/FILLED_EPOXY models */
      
    case TAGC_MAXPACK: 
      *lambda = gn_glob[mn]->maxpack;
      break;
    case TAGC_FICKDIFF_X:
      *lambda = mp_glob[mn]->u_fdiffusivity[0][0];
      break;
    case TAGC_FICKDIFF_Y:
      *lambda = mp_glob[mn]->u_fdiffusivity[0][1];
      break;
      
      /* these can be implemented for Ryan's Qtensor model as needed */

#if 0
    case TAGC_QTENSOR_EXTENSION_P:
      *lambda = mp_glob[mn]->u_Qtensor_Extension_P[0][0];
      break;
    case TAGC_QTENSOR_NCT:
      *lambda = mp_glob[mn]->u_Qtensor_Nct[0][0];
      break;
#endif

      /* these are for CURE/EPOXY/FILLED_EPOXY models */
      
    case TAGC_GELPOINT: 
      *lambda = gn_glob[mn]->gelpoint;
      break;
    case TAGC_CUREAEXP: 
      *lambda = gn_glob[mn]->cureaexp;
      break;
    case TAGC_CUREBEXP: 
      *lambda = gn_glob[mn]->curebexp;
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
      *lambda = ve_glob[mn][mpr-TAGC_TIME_CONST]->time_const;
      break;

    case TAGC_WT_FUNC:
      *lambda = vn_glob[mn]->wt_func;
      break;

    case TAGC_ALPHA:
    case TAGC_ALPHA1:
    case TAGC_ALPHA2:
    case TAGC_ALPHA3:
    case TAGC_ALPHA4:
    case TAGC_ALPHA5:
    case TAGC_ALPHA6:
    case TAGC_ALPHA7:
      *lambda = ve_glob[mn][mpr-TAGC_ALPHA]->alpha;
      break;

    case TAGC_PTT_XI:
      *lambda = ve_glob[mn][cont->upMPID-TAGC_PTT_XI]->xi;
      break;

    case TAGC_PTT_EPS:
      *lambda = ve_glob[mn][cont->upMPID-TAGC_PTT_EPS]->eps;
      break;

    case TAGC_SHIFT_FUNC:
    case TAGC_SHIFT_FUNC1:
      *lambda = vn_glob[mn]->shift[mpr-TAGC_SHIFT_FUNC];
      break;

      /* 
       * Constants used in the Elasticity Constitutive Equations
       */

    case TAGC_LAME_MU:
      *lambda = elc_glob[mn]->lame_mu;
      break;

    case TAGC_LAME_MU_CONTACT_LINE_G0:
      *lambda = *(elc_glob[mn]->u_mu+1);
      break;

    case TAGC_LAME_MU_CONTACT_LINE_G1:
      *lambda = *(elc_glob[mn]->u_mu+2);
      break;

    case TAGC_LAME_MU_CONTACT_LINE_R0:
      *lambda = *(elc_glob[mn]->u_mu+3);
      break;

    case TAGC_LAME_LAMBDA:
      *lambda = elc_glob[mn]->lame_lambda;
      break;

    case TAGC_CONV_LAG_VELX:
      *lambda = elc_glob[mn]->v_mesh_sfs[0];
      break;

    case TAGC_CONV_LAG_VELY:
      *lambda = elc_glob[mn]->v_mesh_sfs[1];
      break;

    case TAGC_CONV_LAG_VELZ:
      *lambda = elc_glob[mn]->v_mesh_sfs[2];
      break;

    case TAGC_CONV_LAG_ROTRATE:
      *lambda = *(elc_glob[mn]->u_v_mesh_sfs);
      break;

    case TAGC_CONV_LAG_ROT_X0:
      *lambda = *(elc_glob[mn]->u_v_mesh_sfs+1);
      break;
     
    case TAGC_CONV_LAG_ROT_Y0:
      *lambda = *(elc_glob[mn]->u_v_mesh_sfs+2);
      break;

    case TAGC_CONV_LAG_ROT_Z0:
      *lambda = *(elc_glob[mn]->u_v_mesh_sfs+3);
      break;

    case TAGC_RS_LAME_MU:
      *lambda = elc_rs_glob[mn]->lame_mu;
      break;

    case TAGC_RS_LAME_LAMBDA:
      *lambda = elc_rs_glob[mn]->lame_lambda;
      break;

    case TAGC_RS_CONV_LAG_VELX:
      *lambda = elc_rs_glob[mn]->v_mesh_sfs[0];
      break;

    case TAGC_RS_CONV_LAG_VELY:
      *lambda = elc_rs_glob[mn]->v_mesh_sfs[1];
      break;

    case TAGC_RS_CONV_LAG_VELZ:
      *lambda = elc_rs_glob[mn]->v_mesh_sfs[2];
      break;

    case TAGC_RS_CONV_LAG_ROTRATE:
      *lambda = *(elc_rs_glob[mn]->u_v_mesh_sfs);
      break;

    case TAGC_RS_CONV_LAG_ROT_X0:
      *lambda = *(elc_rs_glob[mn]->u_v_mesh_sfs+1);
      break;
     
    case TAGC_RS_CONV_LAG_ROT_Y0:
      *lambda = *(elc_rs_glob[mn]->u_v_mesh_sfs+2);
      break;

    case TAGC_RS_CONV_LAG_ROT_Z0:
      *lambda = *(elc_rs_glob[mn]->u_v_mesh_sfs+3);
      break;

    case TAGC_POISSON:
      *lambda = elc_glob[mn]->poisson;
      break;

    case TAGC_STRSS_FR_SOL_VOL_FRAC:
      *lambda = elc_glob[mn]->Strss_fr_sol_vol_frac;
      break;

      /* 
       * Constants used for Source Term Models
       */

    case TAGC_NSS_A0: 
      *lambda = mp_glob[mn]->momentum_source[0];
      break;
      
    case TAGC_NSS_A1: 
      *lambda = mp_glob[mn]->momentum_source[1];
      break;
      
    case TAGC_NSS_A2: 
      *lambda = mp_glob[mn]->momentum_source[2];
      break;
      
    case TAGC_NSS_A3: 
      *lambda = mp_glob[mn]->u_momentum_source[0];
      break;
      
    case TAGC_SHU_QFLOW:
      *lambda = *(mp_glob[mn]->u_shell_user_par+3);
      break;

    case TAGC_SHU_VWEB:
      *lambda = *(mp_glob[mn]->u_shell_user_par+4);
      break;

    case TAGC_SHU_ROLLRAD:
      *lambda = *(mp_glob[mn]->u_shell_user_par+5);
      break;

    case TAGC_SHU_X0:
      *lambda = *(mp_glob[mn]->u_shell_user_par+6);
      break;

    case TAGC_SHU_GAPN:
      *lambda = *(mp_glob[mn]->u_shell_user_par+7);
      break;

    case TAGC_SHU_UPS_XLOC:
      *lambda = *(mp_glob[mn]->u_shell_user_par+8);
      break;

    case TAGC_SHU_DNS_XLOC:
      *lambda = *(mp_glob[mn]->u_shell_user_par+9);
      break;

    default: 
      printf("\n\t Error: Invalid Material Property Tag %d\n", mpr);
      exit(0);
      break;
    }

}/* END of routine retrieve_MT_parameter  */
/*******************************************************************************/

void
retrieve_AC_parameter(double *lambda, /* Parameter value */
                    int ibc,       /* Boundary condition index */
                    int idf,       /* Boundary condition float tag */
                    Comm_Ex *cx,   /* array of communications structures */
                    Exo_DB *exo,   /* ptr to the finite element mesh database */
                    Dpi *dpi)      /* distributed processing information */
{

 	if ( idf == -1 || augc[ibc].len_AC == 0 )	
 		{
 		*lambda = augc[ibc].CONSTV;
 		}
 	else
 		{
 		*lambda = augc[ibc].DataFlt[idf];
 		}
}/* END of routine update_AC_parameter  */
/*******************************************************************************/

void
retrieve_UM_parameter(double *lambda, /* Parameter value */
                    int mn,        /* Material number index */
                    int mpr,       /* Material property index */
                    int ic,        /* User model float index */
                    Comm_Ex *cx,   /* array of communications structures */
                    Exo_DB *exo,   /* ptr to the finite element mesh database */                    Dpi *dpi)      /* distributed processing information */
{

  int float_length = -1;
  switch (mpr)
    {
      case TAGC_THERMAL_CONDUCTIVITY:
        float_length = mp_glob[mn]->len_u_thermal_conductivity;
        break;
      case TAGC_ELECTRICAL_CONDUCTIVITY:
        float_length = mp_glob[mn]->len_u_electrical_conductivity;
        break;
      case TAGC_VISCOSITY:
        float_length = mp_glob[mn]->len_u_viscosity;
        break;
      case TAGC_SURFACE_TENSION:
        float_length = mp_glob[mn]->len_u_surface_tension;
        break;
      case TAGC_HEAT_CAPACITY:
        float_length = mp_glob[mn]->len_u_heat_capacity;
        break;
      case TAGC_VOLUME_EXPANSION:
        float_length = mp_glob[mn]->len_u_Volume_Expansion;
        break;
      case TAGC_DENSITY:
        float_length = mp_glob[mn]->len_u_density;
        break;
      case TAGC_POROSITY:
        float_length = mp_glob[mn]->len_u_porosity;
        break;
      case TAGC_PERMEABILITY:
        float_length = mp_glob[mn]->len_u_permeability;
        break;
      case TAGC_REL_GAS_PERM:
        float_length = mp_glob[mn]->len_u_rel_gas_perm;
        break;
      case TAGC_REL_LIQ_PERM:
        float_length = mp_glob[mn]->len_u_rel_liq_perm;
        break;
      case TAGC_SATURATION:
        float_length = mp_glob[mn]->len_u_saturation;
        break;
      case TAGC_FLOWINGLIQUID_VISCOSITY:
        float_length = mp_glob[mn]->len_u_FlowingLiquid_viscosity;
        break;
      case TAGC_TAU_Y:
        float_length = gn_glob[mn]->len_u_tau_y;
        break;
      default:
        EH(-1, "No model available for that property!");
        break;
    }
        if (ic >= float_length)
          {
            EH(-1, "Float index larger than user parameter list length!");
          }


  switch (mpr)
    {
      case TAGC_THERMAL_CONDUCTIVITY:
        *lambda = mp_glob[mn]->u_thermal_conductivity[ic];
        break;
      case TAGC_ELECTRICAL_CONDUCTIVITY:
        *lambda = mp_glob[mn]->u_electrical_conductivity[ic];
        break;
      case TAGC_VISCOSITY:
        *lambda = mp_glob[mn]->u_viscosity[ic];
        break;
      case TAGC_SURFACE_TENSION:
        *lambda = mp_glob[mn]->u_surface_tension[ic];
        break;
      case TAGC_HEAT_CAPACITY:
        *lambda = mp_glob[mn]->u_heat_capacity[ic];
        break;
      case TAGC_VOLUME_EXPANSION:
        *lambda = mp_glob[mn]->u_Volume_Expansion[ic];
        break;
      case TAGC_DENSITY:
        *lambda = mp_glob[mn]->u_density[ic];
        break;
      case TAGC_POROSITY:
        *lambda = mp_glob[mn]->u_porosity[ic];
        break;
      case TAGC_PERMEABILITY:
        *lambda = mp_glob[mn]->u_permeability[ic];
        break;
      case TAGC_REL_GAS_PERM:
        *lambda = mp_glob[mn]->u_rel_gas_perm[ic];
        break;
      case TAGC_REL_LIQ_PERM:
        *lambda = mp_glob[mn]->u_rel_liq_perm[ic];
        break;
      case TAGC_SATURATION:
        *lambda = mp_glob[mn]->u_saturation[ic];
        break;
      case TAGC_FLOWINGLIQUID_VISCOSITY:
        *lambda = mp_glob[mn]->u_FlowingLiquid_viscosity[ic];
        break;
      case TAGC_TAU_Y:
        *lambda = gn_glob[mn]->u_tau_y[ic];
        break;
      default:
        EH(-1, "No model available for that property!");
        break;
    }
}/* END of routine retrieve_UM_parameter  */
/*******************************************************************************/

/* END of file ac_update_parameter.c  */
