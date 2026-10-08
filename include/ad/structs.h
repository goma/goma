/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) Goma Developers, National Technology & Engineering        *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/
#ifndef GOMA_AD_STRUCTS_H
#define GOMA_AD_STRUCTS_H

#ifdef GOMA_ENABLE_SACADO

#ifdef __cplusplus
#include <Sacado.hpp>
extern "C" {
#include "el_elm.h"
#include "mm_mp_const.h"
}
using ADType = Sacado::Fad::DFad<double>;
struct AD_Basis {
  ADType d_phi[MDE][DIM];                /* d_phi[i][a]    = d(phi_i)/d(q_a) */
  ADType grad_phi[MDE][DIM];             /* grad_phi[i][a] = e_a . grad(phi_i) */
  ADType grad_phi_e[MDE][DIM][DIM][DIM]; /* grad_phi_e[i][a][p][q] */
                                         /* = (e_p e_q): grad(phi_i e_a) */
  ADType curl_phi_e[MDE][DIM][DIM];
};
struct AD_Field_Variables {
  AD_Field_Variables() = default;
  std::vector<AD_Basis> basis;
  ADType detJ;
  ADType J[DIM][DIM];
  ADType B[DIM][DIM];
  ADType v[DIM];
  ADType v_dot[DIM];
  ADType x[DIM];
  ADType d[DIM];
  ADType x_dot[DIM];
  ADType grad_v[DIM][DIM];
  ADType G[DIM][DIM];
  ADType grad_G[DIM][DIM][DIM];
  ADType div_G[DIM];
  ADType S[MAX_MODES][DIM][DIM];
  ADType S_dot[MAX_MODES][DIM][DIM];
  ADType grad_S[MAX_MODES][DIM][DIM][DIM];
  ADType div_S[MAX_MODES][DIM];
  ADType grad_SH[DIM];
  ADType P;
  ADType T;
  ADType SH;
  ADType grad_P[DIM];
  ADType eddy_nu;
  ADType eddy_nu_dot;
  ADType grad_eddy_nu[DIM];
  ADType turb_k;
  ADType turb_k_dot;
  ADType grad_turb_k[DIM];
  ADType turb_omega;
  ADType turb_omega_dot;
  ADType film_height;
  ADType film_height_dot;
  ADType grad_film_height[DIM];
  ADType grad_turb_omega[DIM];
  ADType sh_sat_1;
  ADType sh_sat_1_dot;
  ADType grad_sh_sat_1[DIM];
  ADType sh_sat_2;
  ADType sh_sat_2_dot;
  ADType grad_sh_sat_2[DIM];
  ADType sh_sat_3;
  ADType sh_sat_3_dot;
  ADType grad_sh_sat_3[DIM];
  int total_ad_variables;
  int ielem;
  int offset[V_LAST];
};

extern AD_Field_Variables *ad_fv;
#endif

#endif

#endif // GOMA_AD_STRUCTS_H