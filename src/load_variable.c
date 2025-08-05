/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2025 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/
#include <string.h>

#include "el_elm.h"
#include "mm_as.h"
#include "mm_eh.h"

int load_variable(double *x_var,       /* variable value */
                  double *d_x_var,     /* sensitivities of variable value */
                  int jvar,            /* variable number */
                  int wspec,           /* species number */
                  double tt,           /* parameter to vary time integration from
                                          explicit (tt = 1) to implicit (tt = 0) */
                  double dt,           /* current time step size */
                  double d_vect_var[]) /* vector sensitivities  */

/******************************************************************************

   Function which calculates the value of a chosen variable at the current point
    (using the field variables); available variable names listed in mm_names.h

   Author:          Rich Cairncross (1511)
   Date:            22 MAY 1995
   Revised:
******************************************************************************/

{
  int var = -1, b;
  *x_var = 0.;
  *d_x_var = 0.;

  memset(d_vect_var, 0, DIM * sizeof(double));

  if (jvar >= D_VEL1_DT && jvar <= D_P_DT) {
    if (pd->TimeIntegration == STEADY)
      GOMA_EH(GOMA_ERROR, "Unsteady GD for Steady problem");
  }

  /* ---- Find variable value, sensitivities, and species number */

  switch (jvar) {
  case VELOCITY1:
    *x_var = fv->v[0];
    var = VELOCITY1;
    *d_x_var = 1.;
    break;
  case VELOCITY2:
    *x_var = fv->v[1];
    var = VELOCITY2;
    *d_x_var = 1.;
    break;
  case VELOCITY3:
    *x_var = fv->v[2];
    var = VELOCITY3;
    *d_x_var = 1.;
    break;
  case USTAR:
    *x_var = fv->v_star[0];
    var = USTAR;
    *d_x_var = 1.;
    break;
  case VSTAR:
    *x_var = fv->v_star[1];
    var = VSTAR;
    *d_x_var = 1.;
    break;
  case WSTAR:
    *x_var = fv->v_star[2];
    var = WSTAR;
    *d_x_var = 1.;
    break;
  case PVELOCITY1:
    *x_var = fv->pv[0];
    var = PVELOCITY1;
    *d_x_var = 1.;
    break;
  case PVELOCITY2:
    *x_var = fv->pv[1];
    var = PVELOCITY2;
    *d_x_var = 1.;
    break;
  case PVELOCITY3:
    *x_var = fv->pv[2];
    var = PVELOCITY3;
    *d_x_var = 1.;
    break;
  case TEMPERATURE:
    *x_var = fv->T;
    var = TEMPERATURE;
    *d_x_var = 1.;
    break;
  case VOLTAGE:
    *x_var = fv->V;
    var = VOLTAGE;
    *d_x_var = 1.;
    break;
  case SURF_CHARGE:
    *x_var = fv->qs;
    var = SURF_CHARGE;
    *d_x_var = 1.;
    break;
  case SHELL_CURVATURE:
    *x_var = fv->sh_K;
    var = SHELL_CURVATURE;
    *d_x_var = 1.;
    break;
  case SHELL_CURVATURE2:
    *x_var = fv->sh_K2;
    var = SHELL_CURVATURE2;
    *d_x_var = 1.;
    break;
  case SHELL_TENSION:
    *x_var = fv->sh_tens;
    var = SHELL_TENSION;
    *d_x_var = 1.;
    break;
  case SHELL_X:
    *x_var = fv->sh_x;
    var = SHELL_X;
    *d_x_var = 1.;
    break;
  case SHELL_Y:
    *x_var = fv->sh_y;
    var = SHELL_Y;
    *d_x_var = 1.;
    break;
  case SHELL_USER:
    *x_var = fv->sh_u;
    var = SHELL_USER;
    *d_x_var = 1.;
    break;
  case SHELL_ANGLE1:
    *x_var = fv->sh_ang[0];
    var = SHELL_ANGLE1;
    *d_x_var = 1.;
    break;
  case SHELL_ANGLE2:
    *x_var = fv->sh_ang[1];
    var = SHELL_ANGLE2;
    *d_x_var = 1.;
    break;
  case SHELL_SURF_DIV_V:
    *x_var = fv->div_s_v;
    var = SHELL_SURF_DIV_V;
    *d_x_var = 1.;
    break;
  case SHELL_SURF_CURV:
    *x_var = fv->curv;
    var = SHELL_SURF_CURV;
    *d_x_var = 1.;
    break;
  case N_DOT_CURL_V:
    *x_var = fv->n_dot_curl_s_v;
    var = N_DOT_CURL_V;
    *d_x_var = 1.;
    break;
  case GRAD_S_V_DOT_N1:
    *x_var = fv->grad_v_dot_n[0];
    var = GRAD_S_V_DOT_N1;
    *d_x_var = 1.;
    break;
  case GRAD_S_V_DOT_N2:
    *x_var = fv->grad_v_dot_n[1];
    var = GRAD_S_V_DOT_N2;
    *d_x_var = 1.;
    break;
  case GRAD_S_V_DOT_N3:
    *x_var = fv->grad_v_dot_n[2];
    var = GRAD_S_V_DOT_N3;
    *d_x_var = 1.;
    break;
  case SHELL_DIFF_FLUX:
    *x_var = fv->sh_J;
    var = SHELL_DIFF_FLUX;
    *d_x_var = 1.;
    break;
  case SHELL_DIFF_CURVATURE:
    *x_var = fv->sh_Kd;
    var = SHELL_DIFF_CURVATURE;
    *d_x_var = 1.;
    break;
  case SHELL_NORMAL1:
    *x_var = fv->n[0];
    var = SHELL_NORMAL1;
    *d_x_var = 1.;
    break;
  case SHELL_NORMAL2:
    *x_var = fv->n[1];
    var = SHELL_NORMAL2;
    *d_x_var = 1.;
    break;
  case SHELL_NORMAL3:
    *x_var = fv->n[2];
    var = SHELL_NORMAL3;
    *d_x_var = 1.;
    break;
  case ACOUS_PREAL:
    *x_var = fv->apr;
    var = ACOUS_PREAL;
    *d_x_var = 1.;
    break;
  case ACOUS_PIMAG:
    *x_var = fv->api;
    var = ACOUS_PIMAG;
    *d_x_var = 1.;
    break;
  case EM_CONT_REAL:
    *x_var = fv->epr;
    var = EM_CONT_REAL;
    *d_x_var = 1.;
    break;
  case EM_CONT_IMAG:
    *x_var = fv->epi;
    var = EM_CONT_IMAG;
    *d_x_var = 1.;
    break;
  case POR_SINK_MASS:
    *x_var = fv->sink_mass;
    var = POR_SINK_MASS;
    *d_x_var = 1.;
    break;
  case ACOUS_REYN_STRESS:
    *x_var = fv->ars;
    var = ACOUS_REYN_STRESS;
    *d_x_var = 1.;
    break;
  case SHELL_BDYVELO:
    *x_var = fv->sh_bv;
    var = SHELL_BDYVELO;
    *d_x_var = 1.;
    break;
  case SHELL_LUBP:
    *x_var = fv->sh_p;
    var = SHELL_LUBP;
    *d_x_var = 1.;
    break;
  case SHELL_FILMP:
    *x_var = fv->sh_fp;
    var = SHELL_FILMP;
    *d_x_var = 1.;
    break;
  case SHELL_FILMH:
    *x_var = fv->sh_fh;
    var = SHELL_FILMH;
    *d_x_var = 1.;
    break;
  case SHELL_PARTC:
    *x_var = fv->sh_pc;
    var = SHELL_PARTC;
    *d_x_var = 1.;
    break;
  case LUBP:
    *x_var = fv->lubp;
    var = LUBP;
    *d_x_var = 1.;
    break;
  case LUBP_2:
    *x_var = fv->lubp_2;
    var = LUBP_2;
    *d_x_var = 1.;
    break;
  case SHELL_SAT_CLOSED:
    *x_var = fv->sh_sat_closed;
    var = SHELL_SAT_CLOSED;
    *d_x_var = 1.;
    break;
  case SHELL_PRESS_OPEN:
    *x_var = fv->sh_p_open;
    var = SHELL_PRESS_OPEN;
    *d_x_var = 1.;
    break;
  case SHELL_PRESS_OPEN_2:
    *x_var = fv->sh_p_open_2;
    var = SHELL_PRESS_OPEN_2;
    *d_x_var = 1.;
    break;
  case SHELL_SAT_1:
    *x_var = fv->sh_sat_1;
    var = SHELL_SAT_1;
    *d_x_var = 1.;
    break;
  case SHELL_SAT_2:
    *x_var = fv->sh_sat_2;
    var = SHELL_SAT_2;
    *d_x_var = 1.;
    break;
  case SHELL_SAT_3:
    *x_var = fv->sh_sat_3;
    var = SHELL_SAT_3;
    *d_x_var = 1.;
    break;
  case SHELL_TEMPERATURE:
    *x_var = fv->sh_t;
    var = SHELL_TEMPERATURE;
    *d_x_var = 1.;
    break;
  case SHELL_DELTAH:
    *x_var = fv->sh_dh;
    var = SHELL_DELTAH;
    *d_x_var = 1.;
    break;
  case SHELL_LUB_CURV:
    *x_var = fv->sh_l_curv;
    var = SHELL_LUB_CURV;
    *d_x_var = 1.;
    break;
  case SHELL_LUB_CURV_2:
    *x_var = fv->sh_l_curv_2;
    var = SHELL_LUB_CURV_2;
    *d_x_var = 1.;
    break;
  case SHELL_SAT_GASN:
    *x_var = fv->sh_sat_gasn;
    var = SHELL_SAT_GASN;
    *d_x_var = 1.;
    break;
  case SHELL_SHEAR_TOP:
    *x_var = fv->sh_shear_top;
    var = SHELL_SHEAR_TOP;
    *d_x_var = 1.;
    break;
  case SHELL_SHEAR_BOT:
    *x_var = fv->sh_shear_bot;
    var = SHELL_SHEAR_BOT;
    *d_x_var = 1.;
    break;
  case SHELL_CROSS_SHEAR:
    *x_var = fv->sh_cross_shear;
    var = SHELL_CROSS_SHEAR;
    *d_x_var = 1.;
    break;
  case TFMP_PRES:
    *x_var = fv->tfmp_pres;
    var = TFMP_PRES;
    *d_x_var = 1;
    break;
  case TFMP_SAT:
    *x_var = fv->tfmp_sat;
    var = TFMP_SAT;
    *d_x_var = 1;
    break;
  case MAX_STRAIN:
    *x_var = fv->max_strain;
    var = MAX_STRAIN;
    *d_x_var = 1.;
    break;
  case CUR_STRAIN:
    *x_var = fv->cur_strain;
    var = CUR_STRAIN;
    *d_x_var = 1.;
    break;
  case EDDY_NU:
    *x_var = fv->eddy_nu;
    var = EDDY_NU;
    *d_x_var = 1.;
    break;
  case TURB_K:
    *x_var = fv->turb_k;
    var = TURB_K;
    *d_x_var = 1.;
    break;
  case TURB_OMEGA:
    *x_var = fv->turb_omega;
    var = TURB_OMEGA;
    *d_x_var = 1.;
    break;
  case LIGHT_INTP:
    *x_var = fv->poynt[0];
    var = LIGHT_INTP;
    *d_x_var = 1.;
    break;
  case LIGHT_INTM:
    *x_var = fv->poynt[1];
    var = LIGHT_INTM;
    *d_x_var = 1.;
    break;
  case LIGHT_INTD:
    *x_var = fv->poynt[2];
    var = LIGHT_INTD;
    *d_x_var = 1.;
    break;
  case RESTIME:
    *x_var = fv->restime;
    var = RESTIME;
    *d_x_var = 1.;
    break;
  case EM_E1_REAL:
    *x_var = fv->em_er[0];
    var = EM_E1_REAL;
    *d_x_var = 1.;
    break;
  case EM_E2_REAL:
    *x_var = fv->em_er[1];
    var = EM_E2_REAL;
    *d_x_var = 1.;
    break;
  case EM_E3_REAL:
    *x_var = fv->em_er[2];
    var = EM_E3_REAL;
    *d_x_var = 1.;
    break;
  case EM_E1_IMAG:
    *x_var = fv->em_ei[0];
    var = EM_E1_IMAG;
    *d_x_var = 1.;
    break;
  case EM_E2_IMAG:
    *x_var = fv->em_ei[1];
    var = EM_E2_IMAG;
    *d_x_var = 1.;
    break;
  case EM_E3_IMAG:
    *x_var = fv->em_ei[2];
    var = EM_E3_IMAG;
    *d_x_var = 1.;
    break;
  case EM_H1_REAL:
    *x_var = fv->em_hr[0];
    var = EM_H1_REAL;
    *d_x_var = 1.;
    break;
  case EM_H2_REAL:
    *x_var = fv->em_hr[1];
    var = EM_H2_REAL;
    *d_x_var = 1.;
    break;
  case EM_H3_REAL:
    *x_var = fv->em_hr[2];
    var = EM_H3_REAL;
    *d_x_var = 1.;
    break;
  case EM_H1_IMAG:
    *x_var = fv->em_hi[0];
    var = EM_H1_IMAG;
    *d_x_var = 1.;
    break;
  case EM_H2_IMAG:
    *x_var = fv->em_hi[1];
    var = EM_H2_IMAG;
    *d_x_var = 1.;
    break;
  case EM_H3_IMAG:
    *x_var = fv->em_hi[2];
    var = EM_H3_IMAG;
    *d_x_var = 1.;
    break;

  case MASS_FRACTION:
    *x_var = fv->c[wspec];
    var = MASS_FRACTION;
    *d_x_var = 1.;
    break;
  case MESH_DISPLACEMENT1:
    *x_var = fv->d[0];
    var = MESH_DISPLACEMENT1;
    *d_x_var = 1.;
    break;
  case MESH_DISPLACEMENT2:
    *x_var = fv->d[1];
    var = MESH_DISPLACEMENT2;
    *d_x_var = 1.;
    break;
  case MESH_DISPLACEMENT3:
    *x_var = fv->d[2];
    var = MESH_DISPLACEMENT3;
    *d_x_var = 1.;
    break;
  case SOLID_DISPLACEMENT1:
    *x_var = fv->d_rs[0];
    var = SOLID_DISPLACEMENT1;
    *d_x_var = 1.;
    break;
  case SOLID_DISPLACEMENT2:
    *x_var = fv->d_rs[1];
    var = SOLID_DISPLACEMENT2;
    *d_x_var = 1.;
    break;
  case SOLID_DISPLACEMENT3:
    *x_var = fv->d_rs[2];
    var = SOLID_DISPLACEMENT3;
    *d_x_var = 1.;
    break;
  case SURFACE:
    GOMA_EH(GOMA_ERROR, "SURFACE variables not defined yet");
    break;
  case PRESSURE:
    *x_var = fv->P;
    var = PRESSURE;
    *d_x_var = 1.;
    break;
  case SHEAR_RATE:
    *x_var = fv->SH;
    var = SHEAR_RATE;
    *d_x_var = 1.;
    break;

  case EXT_VELOCITY:
    *x_var = fv->ext_v;
    var = EXT_VELOCITY;
    *d_x_var = 1.;
    break;

  case EFIELD1:
    *x_var = fv->E_field[0];
    var = EFIELD1;
    *d_x_var = 1.;
    break;

  case EFIELD2:
    *x_var = fv->E_field[1];
    var = EFIELD2;
    *d_x_var = 1.;
    break;
  case EFIELD3:
    *x_var = fv->E_field[2];
    var = EFIELD3;
    *d_x_var = 1.;
    break;

  case ENORM:
    *x_var = fv->Enorm;
    var = ENORM;
    *d_x_var = 1.;
    break;

  case CURVATURE:
    *x_var = fv->H;
    var = CURVATURE;
    *d_x_var = 1.;
    break;

  case NORMAL1:
    *x_var = fv->n[0];
    var = NORMAL1;
    *d_x_var = 1;
    break;

  case NORMAL2:
    *x_var = fv->n[1];
    var = NORMAL2;
    *d_x_var = 1;
    break;

  case NORMAL3:
    *x_var = fv->n[2];
    var = NORMAL3;
    *d_x_var = 1;
    break;

  case POLYMER_STRESS11:
    *x_var = fv->S[0][0][0];
    var = POLYMER_STRESS11;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS12:
    *x_var = fv->S[0][0][1];
    var = POLYMER_STRESS12;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS22:
    *x_var = fv->S[0][1][1];
    var = POLYMER_STRESS22;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS13:
    *x_var = fv->S[0][0][2];
    var = POLYMER_STRESS13;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS23:
    *x_var = fv->S[0][1][2];
    var = POLYMER_STRESS23;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS33:
    *x_var = fv->S[0][2][2];
    var = POLYMER_STRESS33;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS11_1:
    *x_var = fv->S[1][0][0];
    var = POLYMER_STRESS11_1;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS12_1:
    *x_var = fv->S[1][0][1];
    var = POLYMER_STRESS12_1;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS22_1:
    *x_var = fv->S[1][1][1];
    var = POLYMER_STRESS22_1;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS13_1:
    *x_var = fv->S[1][0][2];
    var = POLYMER_STRESS13_1;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS23_1:
    *x_var = fv->S[1][1][2];
    var = POLYMER_STRESS23_1;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS33_1:
    *x_var = fv->S[1][2][2];
    var = POLYMER_STRESS33_1;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS11_2:
    *x_var = fv->S[2][0][0];
    var = POLYMER_STRESS11_2;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS12_2:
    *x_var = fv->S[2][0][1];
    var = POLYMER_STRESS12_2;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS22_2:
    *x_var = fv->S[2][1][1];
    var = POLYMER_STRESS22_2;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS13_2:
    *x_var = fv->S[2][0][2];
    var = POLYMER_STRESS13_2;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS23_2:
    *x_var = fv->S[2][1][2];
    var = POLYMER_STRESS23_2;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS33_2:
    *x_var = fv->S[2][2][2];
    var = POLYMER_STRESS33_2;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS11_3:
    *x_var = fv->S[3][0][0];
    var = POLYMER_STRESS11_3;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS12_3:
    *x_var = fv->S[3][0][1];
    var = POLYMER_STRESS12_3;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS22_3:
    *x_var = fv->S[3][1][1];
    var = POLYMER_STRESS22_3;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS13_3:
    *x_var = fv->S[3][0][2];
    var = POLYMER_STRESS13_3;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS23_3:
    *x_var = fv->S[3][1][2];
    var = POLYMER_STRESS23_3;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS33_3:
    *x_var = fv->S[3][2][2];
    var = POLYMER_STRESS33_3;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS11_4:
    *x_var = fv->S[4][0][0];
    var = POLYMER_STRESS11_4;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS12_4:
    *x_var = fv->S[4][0][1];
    var = POLYMER_STRESS12_4;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS22_4:
    *x_var = fv->S[4][1][1];
    var = POLYMER_STRESS22_4;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS13_4:
    *x_var = fv->S[4][0][2];
    var = POLYMER_STRESS13_4;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS23_4:
    *x_var = fv->S[4][1][2];
    var = POLYMER_STRESS23_4;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS33_4:
    *x_var = fv->S[4][2][2];
    var = POLYMER_STRESS33_4;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS11_5:
    *x_var = fv->S[5][0][0];
    var = POLYMER_STRESS11_5;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS12_5:
    *x_var = fv->S[5][0][1];
    var = POLYMER_STRESS12_5;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS22_5:
    *x_var = fv->S[5][1][1];
    var = POLYMER_STRESS22_5;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS13_5:
    *x_var = fv->S[5][0][2];
    var = POLYMER_STRESS13_5;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS23_5:
    *x_var = fv->S[5][1][2];
    var = POLYMER_STRESS23_5;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS33_5:
    *x_var = fv->S[5][2][2];
    var = POLYMER_STRESS33_5;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS11_6:
    *x_var = fv->S[6][0][0];
    var = POLYMER_STRESS11_6;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS12_6:
    *x_var = fv->S[6][0][1];
    var = POLYMER_STRESS12_6;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS22_6:
    *x_var = fv->S[6][1][1];
    var = POLYMER_STRESS22_6;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS13_6:
    *x_var = fv->S[6][0][2];
    var = POLYMER_STRESS13_6;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS23_6:
    *x_var = fv->S[6][1][2];
    var = POLYMER_STRESS23_6;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS33_6:
    *x_var = fv->S[6][2][2];
    var = POLYMER_STRESS33_6;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS11_7:
    *x_var = fv->S[7][0][0];
    var = POLYMER_STRESS11_7;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS12_7:
    *x_var = fv->S[7][0][1];
    var = POLYMER_STRESS12_7;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS22_7:
    *x_var = fv->S[7][1][1];
    var = POLYMER_STRESS22_7;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS13_7:
    *x_var = fv->S[7][0][2];
    var = POLYMER_STRESS13_7;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS23_7:
    *x_var = fv->S[7][1][2];
    var = POLYMER_STRESS23_7;
    *d_x_var = 1.;
    break;

  case POLYMER_STRESS33_7:
    *x_var = fv->S[7][2][2];
    var = POLYMER_STRESS33_7;
    *d_x_var = 1.;
    break;

  case VELOCITY_GRADIENT11:
    *x_var = fv->G[0][0];
    var = VELOCITY_GRADIENT11;
    *d_x_var = 1.;
    break;

  case VELOCITY_GRADIENT12:
    *x_var = fv->G[0][1];
    var = VELOCITY_GRADIENT12;
    *d_x_var = 1.;
    break;

  case VELOCITY_GRADIENT21:
    *x_var = fv->G[1][0];
    var = VELOCITY_GRADIENT21;
    *d_x_var = 1.;
    break;

  case VELOCITY_GRADIENT22:
    *x_var = fv->G[1][1];
    var = VELOCITY_GRADIENT22;
    *d_x_var = 1.;
    break;

  case VELOCITY_GRADIENT13:
    *x_var = fv->G[0][2];
    var = VELOCITY_GRADIENT13;
    *d_x_var = 1.;
    break;

  case VELOCITY_GRADIENT23:
    *x_var = fv->G[1][2];
    var = VELOCITY_GRADIENT23;
    *d_x_var = 1.;
    break;

  case VELOCITY_GRADIENT31:
    *x_var = fv->G[2][0];
    var = VELOCITY_GRADIENT31;
    *d_x_var = 1.;
    break;

  case VELOCITY_GRADIENT32:
    *x_var = fv->G[2][1];
    var = VELOCITY_GRADIENT32;
    *d_x_var = 1.;
    break;

  case VELOCITY_GRADIENT33:
    *x_var = fv->G[2][2];
    var = VELOCITY_GRADIENT33;
    *d_x_var = 1.;
    break;

  case PHASE1:
  case PHASE2:
  case PHASE3:
  case PHASE4:
  case PHASE5:

    b = jvar - PHASE1;
    *x_var = fv->pF[b];
    var = PHASE1 + b;
    *d_x_var = 1.;
    break;

    /* if variable type is mesh position **not** mesh displacement*/
  case MESH_POSITION1:
    *x_var = fv->x[0];
    var = MESH_DISPLACEMENT1;
    *d_x_var = 1.;
    break;
  case MESH_POSITION2:
    *x_var = fv->x[1];
    var = MESH_DISPLACEMENT2;
    *d_x_var = 1.;
    break;
  case MESH_POSITION3:
    *x_var = fv->x[2];
    var = MESH_DISPLACEMENT3;
    *d_x_var = 1.;
    break;
  case SOLID_POSITION1:
    *x_var = fv->x[0];
    var = SOLID_DISPLACEMENT1;
    *d_x_var = 1.;
    break;
  case SOLID_POSITION2:
    *x_var = fv->x[1];
    var = SOLID_DISPLACEMENT2;
    *d_x_var = 1.;
    break;
  case SOLID_POSITION3:
    *x_var = fv->x[2];
    var = SOLID_DISPLACEMENT3;
    *d_x_var = 1.;
    break;
  case POR_LIQ_PRES:
    *x_var = fv->p_liq;
    var = POR_LIQ_PRES;
    *d_x_var = 1.;
    break;
  case POR_GAS_PRES:
    *x_var = fv->p_gas;
    var = POR_GAS_PRES;
    *d_x_var = 1.;
    break;
  case POR_POROSITY:
    *x_var = fv->porosity;
    var = POR_GAS_PRES;
    *d_x_var = 1.;
    break;
    /* adding velocity magnitude, i.e. SPEED  */
  case SPEED:
    for (b = 0; b < pd->Num_Dim; b++) {
      *x_var += SQUARE(fv->v[b]);
    }
    if (pd->CoordinateSystem == SWIRLING || pd->CoordinateSystem == PROJECTED_CARTESIAN ||
        pd->CoordinateSystem == CARTESIAN_2pt5D) {
      *x_var += SQUARE(fv->v[pd->Num_Dim]);
    }
    *x_var = sqrt(*x_var);
    var = VELOCITY1;
    *d_x_var = 1. / (*x_var);
    for (b = 0; b < pd->Num_Dim; b++) {
      d_vect_var[b] += fv->v[b] * (*d_x_var);
    }
    if (pd->CoordinateSystem == SWIRLING || pd->CoordinateSystem == PROJECTED_CARTESIAN ||
        pd->CoordinateSystem == CARTESIAN_2pt5D) {
      d_vect_var[pd->Num_Dim] += fv->v[pd->Num_Dim] * (*d_x_var);
    }
    break;

    /* if variable type is a time derivative */
  case D_VEL1_DT:
    *x_var = fv_dot->v[0];
    var = VELOCITY1;
    *d_x_var = (1. + 2. * tt) / dt;
    break;
  case D_VEL2_DT:
    *x_var = fv_dot->v[1];
    var = VELOCITY2;
    *d_x_var = (1. + 2. * tt) / dt;
    break;
  case D_VEL3_DT:
    *x_var = fv_dot->v[2];
    var = VELOCITY3;
    *d_x_var = (1. + 2. * tt) / dt;
    break;
  case D_T_DT:
    *x_var = fv_dot->T;
    var = TEMPERATURE;
    *d_x_var = (1. + 2. * tt) / dt;
    break;
  case D_Y_DT:
    *x_var = fv_dot->c[wspec];
    var = MASS_FRACTION;
    *d_x_var = (1. + 2. * tt) / dt;
    break;
  case D_X1_DT:
    *x_var = fv_dot->d[0];
    var = MESH_DISPLACEMENT1;
    *d_x_var = (1. + 2. * tt) / dt;
    break;
  case D_X2_DT:
    *x_var = fv_dot->d[1];
    var = MESH_DISPLACEMENT2;
    *d_x_var = (1. + 2. * tt) / dt;
    break;
  case D_X3_DT:
    *x_var = fv_dot->d[2];
    var = MESH_DISPLACEMENT3;
    *d_x_var = (1. + 2. * tt) / dt;
    break;
  case D_S_DT:
    GOMA_EH(GOMA_ERROR, "SURFACE variables not defined yet");
    break;
  case D_P_DT:
    *x_var = fv_dot->P;
    var = PRESSURE;
    *d_x_var = (1. + 2. * tt) / dt;
    break;
  default:
    GOMA_EH(GOMA_ERROR, "Illegal option in load_variable");
  } /* end of switch on jvar */

  if (wspec != 0 && var != MASS_FRACTION && var != POR_GAS_PRES) {
    GOMA_EH(GOMA_ERROR, "Non-zero species number for wrong variable");
  }
  if (var == MASS_FRACTION)
    var = MAX_VARIABLE_TYPES + wspec;

  return (var);
} /* END of routine load_variable()                                          */