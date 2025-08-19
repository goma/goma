#ifdef GOMA_ENABLE_SACADO

#include <ad_turbulence.h>

#include <ad_porous.h>

extern "C" {
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_mp.h"
#include "mm_shell_util.h"
#include "mm_std_models_shell.h"
#include "rf_allo.h"
#include "rf_fem_const.h"
}

static inline ADType set_ad_or_dbl(dbl val, int eqn, int dof) {
  ADType tmp;
  if (ad_fv->total_ad_variables > 0 && af->Assemble_Jacobian == TRUE) {
    if (pd->gv[eqn]) {
      tmp = ADType(ad_fv->total_ad_variables, ad_fv->offset[eqn] + dof, val);
    } else {
      tmp = val;
    }
  } else {
    tmp = val;
  }
  return tmp;
}

ADType ad_load_cap_pres(int ipore, int ilnode, int ignode, ADType saturation)

/*************************************************************************
 * load_cap_pres -- calculate capillary pressure  in a porous
 *                    media from the unknowns used in the problem
 *
 *  input:   Assume that load_fv and load_fv_grads and load_fv_mesh_derivs
 *  -----    have all been called already, so unknowns and their
 *           sensitivies are known at this gauss point.
 *
 *  output:  calculates the saturation and its first and second
 *  -----    derivatives with respect to all the problem unknowns
 *           For transient calculations, the saturation at the
 *           old time step must be recalculated as well, in order
 *           to be used in the capacitance term.
 *
 *      mp->cap_pres
 *      mp->d_cap_pres
 *      mp->d_d_cap_pres
 *	    mp_old->cap_pres
 ************************************************************************/
{
  ADType con_a, con_b, con_c, con_d;
  ADType sat_min, sat_max;
  ADType sat_norm, sat_clip, d_sat_norm_d_saturation, d_sat_clip_d_saturation;
  ADType cap_pres = 0.0;
  ADType mexp, nexp, n_inv, m_inv;
  ADType brack_in, d_brack_in_d_sat_norm, d_d_brack_in_d_sat_norm;

  int i_ext_field;
  ADType val_ext_field;
  ADType sat_min_1, sat_max_1, sat_min_2, sat_max_2;
  ADType con_c_1, nexp_1, con_c_2, nexp_2;

  int draining_curve, draining_curve_old;
  int switch_now;
  ADType cap_pres_switch, sat_switch, sat_norm_switch;
  ADType brack_in_inv, sat_min_wet, sat_max_dry;

  /*
   *  Find which model is used for saturation and calculate:
   *    1)  The capillary pressure, mp->cap_pres,
   *    2)  The sensitivity of saturation to all concentrations, porosity
   *        and temperature
   *          i.e., the first derivative of saturation w.r.t. each variable
   *           put this in mp->d_cap_pres[var]
   *    3)  The second derivatives of saturation w.r.t. each variable,
   *        including cross-terms, put this in mp->d_d_cap_pres[var][var]
   *        The second derivative is needed for sensitivity of fluxes
   *        which depend on the first derivative of capillary pressure
   */

  /**********************************************************************
   *                   ATANH MODEL FOR CAPILLARY PRESSURE
   **********************************************************************/
  if (mp->PorousShellCapPresModel[ipore] == ATANH) {
    /*
     * FOR ATANH EQUATION
     *  mp->u_saturation[0] is the irreduceable water saturation
     *  mp->u_saturation[1] is the irreduceable air saturation
     *  mp->u_saturation[2] is shift factor
     *  mp->u_saturation[3] is multiplier in tanh function
     */
    //    sat_min = mp->u_cap_pres[0];
    //    sat_max = mp->u_cap_pres[1];
    //    con_c = mp->u_cap_pres[2];
    //    con_d = mp->u_cap_pres[3];

    /* Normalized saturation - should be between -1 to +1 */
    //    sat_norm = (2.0 * saturation - sat_max - sat_min) / (sat_max - sat_min);
    //    d_sat_norm_d_saturation = 2.0/(sat_max - sat_min);

    con_a = mp->u_PorousShellCapPres[ipore][0];
    con_b = mp->u_PorousShellCapPres[ipore][1];
    con_c = mp->u_PorousShellCapPres[ipore][2];
    con_d = mp->u_PorousShellCapPres[ipore][3];

    sat_norm = (con_a - saturation) / con_b;
    d_sat_norm_d_saturation = -1.0 / con_b;

    /* Clip if necessary */
    if (sat_norm <= -0.995) {
      sat_clip = -0.995;
    } else if (sat_norm >= 0.995) {
      sat_clip = 0.995;
    } else {
      sat_clip = sat_norm;
    }
    //    cap_pres = mp->cap_pres = con_c - con_d * atanh(sat_clip);
    cap_pres = con_d / (con_c - atanh(sat_clip));
    mp->cap_pres = cap_pres.val();

  } else if (mp->PorousShellCapPresModel[ipore] == SINH) {
    /*
     * FOR SINH EQUATION
     *  mp->u_saturation[0] is the irreducable water saturation
     *  mp->u_saturation[1] is the irreduceable air saturation
     *  mp->u_saturation[2] is shift factor
     *  mp->u_saturation[3] is multiplier in tanh function
     */

    sat_min = mp->u_PorousShellCapPres[ipore][0];
    sat_max = mp->u_PorousShellCapPres[ipore][1];
    con_c = mp->u_PorousShellCapPres[ipore][2];
    con_d = mp->u_PorousShellCapPres[ipore][3];

    /* Normalized saturation - should be between -6 to +6 */
    sat_norm = (12.0 * saturation - 6.0 * (sat_max + sat_min)) / (sat_max - sat_min);
    d_sat_norm_d_saturation = 12.0 / (sat_max - sat_min);

    /* Clip if necessary */
    //    if (saturation <= 0.0) {
    //      sat_clip = -6.0;
    //      d_sat_clip_d_saturation = 0.0;
    //    }else if (saturation >= 0.995) {
    //      sat_clip = 6.0;
    //      d_sat_clip_d_saturation = 0.0;
    //    } else {
    sat_clip = sat_norm;
    d_sat_clip_d_saturation = d_sat_norm_d_saturation;
    //    }
    cap_pres = con_c - con_d * sinh(sat_clip);

  } else if (mp->PorousShellCapPresModel[ipore] == VAN_GENUCHTEN) {
    /*
     * FOR VAN_GENUCHTEN EQUATION
     *  mp->u_saturation[0] is the irreducable water saturation
     *  mp->u_saturation[1] is the irreduceable air saturation
     *  mp->u_saturation[2] is entry pressure
     *  mp->u_saturation[3] is exponent n
     */

    sat_min = mp->u_PorousShellCapPres[ipore][0];
    sat_max = mp->u_PorousShellCapPres[ipore][1];
    con_c = mp->u_PorousShellCapPres[ipore][2];
    nexp = mp->u_PorousShellCapPres[ipore][3];

    n_inv = 1.0 / nexp;
    mexp = 1.0 - n_inv;
    m_inv = -1.0 / mexp;

    /* Normalized saturation - should be between 0 to 1 */
    sat_norm = (saturation - sat_min) / (sat_max - sat_min);
    d_sat_norm_d_saturation = 1.0 / (sat_max - sat_min);

    brack_in = pow(sat_norm, m_inv) - 1.0;
    d_brack_in_d_sat_norm = m_inv * pow(sat_norm, (m_inv - 1.0));
    d_d_brack_in_d_sat_norm = m_inv * (m_inv - 1.0) * pow(sat_norm, (m_inv - 2.0));

    cap_pres = con_c * pow(brack_in, n_inv);

  } else if (mp->PorousShellCapPresModel[ipore] == VAN_GENUCHTEN_EXTERNAL) {
    /*
     * FOR VAN_GENUCHTEN EXTERNAL EQUATION
     *  mp->u_saturation[0] is the irreducable water saturation for curve 1
     *  mp->u_saturation[1] is the irreduceable air saturation for curve 1
     *  mp->u_saturation[2] is entry pressure for curve 1
     *  mp->u_saturation[3] is exponent n for curve 1
     *
     *  mp->u_saturation[4] is the irreducable water saturation for curve 2
     *  mp->u_saturation[5] is the irreduceable air saturation for curve 2
     *  mp->u_saturation[6] is entry pressure for curve 2
     *  mp->u_saturation[7] is exponent n for curve 2
     *
     */

    i_ext_field = mp->por_shell_cap_pres_ext_field_index[ipore];
    /* Here I assume the external field value ranges from 0 to 1 */
    val_ext_field = fv->external_field[i_ext_field];

    sat_min_1 = mp->u_PorousShellCapPres[ipore][0];
    sat_max_1 = mp->u_PorousShellCapPres[ipore][1];
    con_c_1 = mp->u_PorousShellCapPres[ipore][2];
    nexp_1 = mp->u_PorousShellCapPres[ipore][3];

    sat_min_2 = mp->u_PorousShellCapPres[ipore][4];
    sat_max_2 = mp->u_PorousShellCapPres[ipore][5];
    con_c_2 = mp->u_PorousShellCapPres[ipore][6];
    nexp_2 = mp->u_PorousShellCapPres[ipore][7];

    /* Interpolate the fitting parameters with external field value */
    sat_min = sat_min_1 + val_ext_field * (sat_min_2 - sat_min_1);
    sat_max = sat_max_1 + val_ext_field * (sat_max_2 - sat_max_1);
    con_c = con_c_1 + val_ext_field * (con_c_2 - con_c_1);
    nexp = nexp_1 + val_ext_field * (nexp_2 - nexp_1);

    n_inv = 1.0 / nexp;
    mexp = 1.0 - n_inv;
    m_inv = -1.0 / mexp;

    /* Normalized saturation - should be between 0 to 1 */
    sat_norm = (saturation - sat_min) / (sat_max - sat_min);
    d_sat_norm_d_saturation = 1.0 / (sat_max - sat_min);

    brack_in = pow(sat_norm, m_inv) - 1.0;
    d_brack_in_d_sat_norm = m_inv * pow(sat_norm, (m_inv - 1.0));

    cap_pres = con_c * pow(std::max(brack_in, 1e-32), n_inv);

  } else if (mp->PorousShellCapPresModel[ipore] == VAN_GENUCHTEN_HYST) {
    /*
     * FOR VAN_GENUCHTEN_HYST EQUATION
     *  mp->u_saturation[0] is the irreducable water saturation for wetting curve
     *  mp->u_saturation[1] is the irreduceable air saturation for wetting curve
     *  mp->u_saturation[2] is entry pressure for wetting curve
     *  mp->u_saturation[3] is exponent n for wetting curve
     *
     *  mp->u_saturation[4] is the irreducable water saturation for draining curve
     *  mp->u_saturation[5] is the irreduceable air saturation for draining curve
     *  mp->u_saturation[6] is entry pressure for draining curve
     *  mp->u_saturation[7] is exponent n for draining curve
     *
     *  Other input
     *  mp->u_saturation[8] Initial saturation curve for the material, viz.
     *                      if 1.0 then on the draining curve, and 0.0 then
     *                      on the wetting curve.
     *  mp->u_saturation[9] is Liquid_inventory_rate threshold for curve switching
     */

    /* Quality check of global node number ignode*/
    if (ignode < 0) {
      if (ipore == 0) {
        ignode = ei[pg->imtrx]->gnn_list[SHELL_SAT_1][ilnode];
      } else if (ipore == 1) {
        ignode = ei[pg->imtrx]->gnn_list[SHELL_SAT_2][ilnode];
      } else if (ipore == 2) {
        ignode = ei[pg->imtrx]->gnn_list[SHELL_SAT_3][ilnode];
      } else {
        GOMA_EH(GOMA_ERROR, "Something's wrong here!");
      }
    }

    /* First, see if we are switching */
    switch_now = pmv_hyst->curve_switch[ipore][ignode];

    /* Get the minimum and maximum saturation for the scanning curves */
    sat_min_wet = pmv_hyst->sat_min_imbibe[ipore][ignode];
    sat_max_dry = pmv_hyst->sat_max_drain[ipore][ignode];

    /* Find out whether on draining or wetting curve */
    draining_curve = pmv_hyst->curve_type[ipore][ignode];
    draining_curve_old = pmv_hyst->curve_type_old[ipore][ignode];

    /* If we are not switching */
    if (switch_now == 0) {
      if (draining_curve == 1) /* Stay on main draining curve */
      {
        sat_min = mp->u_PorousShellCapPres[ipore][4];
        sat_max = sat_max_dry;
        con_c = mp->u_PorousShellCapPres[ipore][6];
        nexp = mp->u_PorousShellCapPres[ipore][7];

        n_inv = 1.0 / nexp;
        mexp = 1.0 - n_inv;
        m_inv = -1.0 / mexp;

        /* Normalized saturation - should be between 0 to 1 */
        sat_norm = (saturation - sat_min) / (sat_max - sat_min);
        d_sat_norm_d_saturation = 1.0 / (sat_max - sat_min);

        brack_in = pow(sat_norm, m_inv) - 1.0;
        d_brack_in_d_sat_norm = m_inv * pow(sat_norm, (m_inv - 1.0));

        /* Evaluate capillary pressure and its sensitivities */
        cap_pres = con_c * pow(brack_in, n_inv);

      } else if (draining_curve == 0) /* Stay on main wetting curve */
      {
        sat_min = sat_min_wet;
        sat_max = mp->u_PorousShellCapPres[ipore][1];
        con_c = mp->u_PorousShellCapPres[ipore][2];
        nexp = mp->u_PorousShellCapPres[ipore][3];

        n_inv = 1.0 / nexp;
        mexp = 1.0 - n_inv;
        m_inv = -1.0 / mexp;

        /* Normalized saturation - should be between 0 to 1 */
        sat_norm = (saturation - sat_min) / (sat_max - sat_min);
        d_sat_norm_d_saturation = 1.0 / (sat_max - sat_min);

        brack_in = pow(sat_norm, m_inv) - 1.0;
        d_brack_in_d_sat_norm = m_inv * pow(sat_norm, (m_inv - 1.0));

        /* Evaluate capillary pressure and its sensitivities */
        cap_pres = con_c * pow(brack_in, n_inv);
      } else {
        GOMA_EH(GOMA_ERROR, "Either we are on draining curve or wetting curve!");
      }
    } else if (switch_now == 1) /* We are switching */
    {
      if ((draining_curve == 1) && (draining_curve_old == 0)) /* Wetting going to draining */
      {
        /* Get all parameters from draining curve */
        sat_min = mp->u_PorousShellCapPres[ipore][4];
        sat_max = mp->u_PorousShellCapPres[ipore][5];
        con_c = mp->u_PorousShellCapPres[ipore][6];
        nexp = mp->u_PorousShellCapPres[ipore][7];

        n_inv = 1.0 / nexp;
        mexp = 1.0 - n_inv;
        m_inv = -1.0 / mexp;

        /* Get the saturation and capillary pressure at switching or reversal point */
        sat_switch = pmv_hyst->sat_switch[ipore][ignode];
        cap_pres_switch = pmv_hyst->cap_pres_switch[ipore][ignode];

        /* Calculate normalized saturation at reversal point */
        brack_in_inv = pow((cap_pres_switch / con_c), nexp) + 1.0;
        sat_norm_switch = pow(brack_in_inv, (-mexp));

        /* Quality check */
        if ((sat_norm_switch > 1.0) || (sat_norm_switch < 0.0)) {
          GOMA_EH(GOMA_ERROR, "Invalid value of normalized saturation at switching point");
        }

        /* Calculate sat_max for the scanning draining curve */
        sat_max_dry = (1.0 / sat_norm_switch) * (sat_switch - sat_min * (1.0 - sat_norm_switch));

        /* More quality check */
        if (sat_max_dry > sat_max) {
          GOMA_EH(GOMA_ERROR, "Invalid value of maximum saturation at the scanning drying curve");
        }
        /* Update sat_max_dry in pmv_hyst structure*/
        pmv_hyst->sat_max_drain[ipore][ignode] = sat_max_dry.val();

        /* Calculate normalized saturation based on newly calculated sat_max_dry */
        sat_norm = (saturation - sat_min) / (sat_max_dry - sat_min);
        d_sat_norm_d_saturation = 1.0 / (sat_max_dry - sat_min);

        brack_in = pow(sat_norm, m_inv) - 1.0;
        d_brack_in_d_sat_norm = m_inv * pow(sat_norm, (m_inv - 1.0));

        /* Evaluate capillary pressure and its sensitivities */
        cap_pres = con_c * pow(brack_in, n_inv);

      } else if ((draining_curve == 0) && (draining_curve_old == 1)) /* draining going to wetting */
      {
        /* Get all parameters from wetting curve */
        sat_min = mp->u_PorousShellCapPres[ipore][0];
        sat_max = mp->u_PorousShellCapPres[ipore][1];
        con_c = mp->u_PorousShellCapPres[ipore][2];
        nexp = mp->u_PorousShellCapPres[ipore][3];

        n_inv = 1.0 / nexp;
        mexp = 1.0 - n_inv;
        m_inv = -1.0 / mexp;

        /* Get the saturation and capillary pressure at switching or reversal point */
        sat_switch = pmv_hyst->sat_switch[ipore][ignode];
        cap_pres_switch = pmv_hyst->cap_pres_switch[ipore][ignode];

        /* Calculate normalized saturation at reversal point */
        brack_in_inv = pow((cap_pres_switch / con_c), nexp) + 1.0;
        sat_norm_switch = pow(brack_in_inv, (-mexp));

        /* Quality check */
        if ((sat_norm_switch > 1.0) || (sat_norm_switch < 0.0)) {
          GOMA_EH(GOMA_ERROR, "Invalid value of normalized saturation at switching point");
        }

        /* Calculate sat_min for the scanning wetting curve */
        sat_min_wet = (sat_switch - sat_max * sat_norm_switch) / (1.0 - sat_norm_switch);

        /* More quality check */
        if (sat_min_wet < sat_min) {
          GOMA_EH(GOMA_ERROR, "Invalid value of minimum saturation at the scanning wetting curve");
        }
        /* Update sat_min_wet in pmv_hyst structure*/
        pmv_hyst->sat_min_imbibe[ipore][ignode] = sat_min_wet.val();

        /* Calculate normalized saturation based on newly calculated sat_min_wet */
        sat_norm = (saturation - sat_min_wet) / (sat_max - sat_min_wet);
        d_sat_norm_d_saturation = 1.0 / (sat_max - sat_min_wet);

        brack_in = pow(sat_norm, m_inv) - 1.0;
        d_brack_in_d_sat_norm = m_inv * pow(sat_norm, (m_inv - 1.0));

        /* Evaluate capillary pressure and its sensitivities */
        cap_pres = con_c * pow(brack_in, n_inv);
      }
    }
  } else if (mp->PorousShellCapPresModel[ipore] == VAN_GENUCHTEN_HYST_EXT) {
    /*
     * FOR VAN_GENUCHTEN_HYST_EXT EQUATION
     *  mp->u_saturation[0] is the irreducable water saturation for wetting curve 1
     *  mp->u_saturation[1] is the irreduceable air saturation for wetting curve 1
     *  mp->u_saturation[2] is entry pressure for wetting curve 1
     *  mp->u_saturation[3] is exponent n for wetting curve 1
     *
     *  mp->u_saturation[4] is the irreducable water saturation for draining curve 1
     *  mp->u_saturation[5] is the irreduceable air saturation for draining curve 1
     *  mp->u_saturation[6] is entry pressure for draining curve 1
     *  mp->u_saturation[7] is exponent n for draining curve 1
     *
     *  mp->u_saturation[8] is the irreducable water saturation for wetting curve 2
     *  mp->u_saturation[9] is the irreduceable air saturation for wetting curve 2
     *  mp->u_saturation[10] is entry pressure for wetting curve 2
     *  mp->u_saturation[11] is exponent n for wetting curve 2
     *
     *  mp->u_saturation[12] is the irreducable water saturation for draining curve 2
     *  mp->u_saturation[13] is the irreduceable air saturation for draining curve 2
     *  mp->u_saturation[14] is entry pressure for draining curve 2
     *  mp->u_saturation[15] is exponent n for draining curve 2

     *  Other input
     *  mp->u_saturation[16] Initial saturation curve for the material, viz.
     *                       if 1.0 then on the draining curve, and 0.0 then
     *                       on the wetting curve.
     *  mp->u_saturation[17] is Liquid_inventory_rate threshold for curve switching
     */

    /* Quality check of global node number ignode*/
    if (ignode < 0) {
      if (ipore == 0) {
        ignode = ei[pg->imtrx]->gnn_list[SHELL_SAT_1][ilnode];
      } else if (ipore == 1) {
        ignode = ei[pg->imtrx]->gnn_list[SHELL_SAT_2][ilnode];
      } else if (ipore == 2) {
        ignode = ei[pg->imtrx]->gnn_list[SHELL_SAT_3][ilnode];
      } else {
        GOMA_EH(GOMA_ERROR, "Something's wrong here!");
      }
    }

    /* First, see if we are switching */
    switch_now = pmv_hyst->curve_switch[ipore][ignode];

    i_ext_field = mp->por_shell_cap_pres_ext_field_index[ipore];
    /* Here I assume the external field value ranges from 0 to 1 */
    val_ext_field = *evp->external_field[i_ext_field][ilnode];

    /* Get the minimum and maximum saturation for the scanning curves */
    sat_min_wet = pmv_hyst->sat_min_imbibe[ipore][ignode];
    sat_max_dry = pmv_hyst->sat_max_drain[ipore][ignode];

    /* Find out whether on draining or wetting curve */
    draining_curve = pmv_hyst->curve_type[ipore][ignode];
    draining_curve_old = pmv_hyst->curve_type_old[ipore][ignode];

    /* If we are not switching */
    if (switch_now == 0) {
      if (draining_curve == 1) /* Stay on main draining curve */
      {

        /* Get parameters from draining curve 1*/
        sat_min_1 = mp->u_PorousShellCapPres[ipore][4];
        con_c_1 = mp->u_PorousShellCapPres[ipore][6];
        nexp_1 = mp->u_PorousShellCapPres[ipore][7];

        /* Get parameters from draining curve 2*/
        sat_min_2 = mp->u_PorousShellCapPres[ipore][12];
        con_c_2 = mp->u_PorousShellCapPres[ipore][14];
        nexp_2 = mp->u_PorousShellCapPres[ipore][15];

        /* Interpolate the fitting parameters with external field value */
        sat_min = sat_min_1 + val_ext_field * (sat_min_2 - sat_min_1);
        con_c = con_c_1 + val_ext_field * (con_c_2 - con_c_1);
        nexp = nexp_1 + val_ext_field * (nexp_2 - nexp_1);

        sat_max = sat_max_dry;

        n_inv = 1.0 / nexp;
        mexp = 1.0 - n_inv;
        m_inv = -1.0 / mexp;

        /* Normalized saturation - should be between 0 to 1 */
        sat_norm = (saturation - sat_min) / (sat_max - sat_min);
        d_sat_norm_d_saturation = 1.0 / (sat_max - sat_min);

        brack_in = pow(sat_norm, m_inv) - 1.0;
        d_brack_in_d_sat_norm = m_inv * pow(sat_norm, (m_inv - 1.0));

        /* Evaluate capillary pressure and its sensitivities */
        cap_pres = con_c * pow(brack_in, n_inv);
      } else if (draining_curve == 0) /* Stay on main wetting curve */
      {
        sat_min = sat_min_wet;

        /* Get parameters from wetting curve 1*/
        sat_max_1 = mp->u_PorousShellCapPres[ipore][1];
        con_c_1 = mp->u_PorousShellCapPres[ipore][2];
        nexp_1 = mp->u_PorousShellCapPres[ipore][3];

        /* Get parameters from wetting curve 2*/
        sat_max_2 = mp->u_PorousShellCapPres[ipore][9];
        con_c_2 = mp->u_PorousShellCapPres[ipore][10];
        nexp_2 = mp->u_PorousShellCapPres[ipore][11];

        /* Interpolate the fitting parameters with external field value */
        sat_max = sat_max_1 + val_ext_field * (sat_max_2 - sat_max_1);
        con_c = con_c_1 + val_ext_field * (con_c_2 - con_c_1);
        nexp = nexp_1 + val_ext_field * (nexp_2 - nexp_1);

        n_inv = 1.0 / nexp;
        mexp = 1.0 - n_inv;
        m_inv = -1.0 / mexp;

        /* Normalized saturation - should be between 0 to 1 */
        sat_norm = (saturation - sat_min) / (sat_max - sat_min);
        d_sat_norm_d_saturation = 1.0 / (sat_max - sat_min);

        brack_in = pow(sat_norm, m_inv) - 1.0;
        d_brack_in_d_sat_norm = m_inv * pow(sat_norm, (m_inv - 1.0));

        /* Evaluate capillary pressure and its sensitivities */
        cap_pres = con_c * pow(brack_in, n_inv);
      } else {
        GOMA_EH(GOMA_ERROR, "Either we are on draining curve or wetting curve!");
      }
    } else if (switch_now == 1) /* We are switching */
    {
      if ((draining_curve == 1) && (draining_curve_old == 0)) /* Wetting going to draining */
      {
        /* Get all parameters from draining curve 1 */
        sat_min_1 = mp->u_PorousShellCapPres[ipore][4];
        sat_max_1 = mp->u_PorousShellCapPres[ipore][5];
        con_c_1 = mp->u_PorousShellCapPres[ipore][6];
        nexp_1 = mp->u_PorousShellCapPres[ipore][7];

        /* Get all parameters from draining curve 2 */
        sat_min_2 = mp->u_PorousShellCapPres[ipore][12];
        sat_max_2 = mp->u_PorousShellCapPres[ipore][13];
        con_c_2 = mp->u_PorousShellCapPres[ipore][14];
        nexp_2 = mp->u_PorousShellCapPres[ipore][15];

        /* Interpolate the fitting parameters with external field value */
        sat_min = sat_min_1 + val_ext_field * (sat_min_2 - sat_min_1);
        sat_max = sat_max_1 + val_ext_field * (sat_max_2 - sat_max_1);
        con_c = con_c_1 + val_ext_field * (con_c_2 - con_c_1);
        nexp = nexp_1 + val_ext_field * (nexp_2 - nexp_1);

        n_inv = 1.0 / nexp;
        mexp = 1.0 - n_inv;
        m_inv = -1.0 / mexp;

        /* Get the saturation and capillary pressure at switching or reversal point */
        sat_switch = pmv_hyst->sat_switch[ipore][ignode];
        cap_pres_switch = pmv_hyst->cap_pres_switch[ipore][ignode];

        /* Calculate normalized saturation at reversal point */
        brack_in_inv = pow((cap_pres_switch / con_c), nexp) + 1.0;
        sat_norm_switch = pow(brack_in_inv, (-mexp));

        /* Quality check */
        if ((sat_norm_switch > 1.0) || (sat_norm_switch < 0.0)) {
          GOMA_EH(GOMA_ERROR, "Invalid value of normalized saturation at switching point");
        }

        /* Calculate sat_max for the scanning draining curve */
        sat_max_dry = (1.0 / sat_norm_switch) * (sat_switch - sat_min * (1.0 - sat_norm_switch));

        /* More quality check */
        if (sat_max_dry > sat_max) {
          GOMA_EH(GOMA_ERROR, "Invalid value of maximum saturation at the scanning drying curve");
        }
        /* Update sat_max_dry in element storage structure*/
        pmv_hyst->sat_max_drain[ipore][ignode] = sat_max_dry.val();

        /* Calculate normalized saturation based on newly calculated sat_max_dry */
        sat_norm = (saturation - sat_min) / (sat_max_dry - sat_min);
        d_sat_norm_d_saturation = 1.0 / (sat_max_dry - sat_min);

        brack_in = pow(sat_norm, m_inv) - 1.0;
        d_brack_in_d_sat_norm = m_inv * pow(sat_norm, (m_inv - 1.0));

        /* Evaluate capillary pressure and its sensitivities */
        cap_pres = con_c * pow(brack_in, n_inv);
      } else if ((draining_curve == 0) && (draining_curve_old == 1)) /* draining going to wetting */
      {
        /* Get all parameters from wetting curve 1 */
        sat_min_1 = mp->u_PorousShellCapPres[ipore][0];
        sat_max_1 = mp->u_PorousShellCapPres[ipore][1];
        con_c_1 = mp->u_PorousShellCapPres[ipore][2];
        nexp_1 = mp->u_PorousShellCapPres[ipore][3];

        /* Get all parameters from wetting curve 2 */
        sat_min_2 = mp->u_PorousShellCapPres[ipore][8];
        sat_max_2 = mp->u_PorousShellCapPres[ipore][9];
        con_c_2 = mp->u_PorousShellCapPres[ipore][10];
        nexp_2 = mp->u_PorousShellCapPres[ipore][11];

        /* Interpolate the fitting parameters with external field value */
        sat_min = sat_min_1 + val_ext_field * (sat_min_2 - sat_min_1);
        sat_max = sat_max_1 + val_ext_field * (sat_max_2 - sat_max_1);
        con_c = con_c_1 + val_ext_field * (con_c_2 - con_c_1);
        nexp = nexp_1 + val_ext_field * (nexp_2 - nexp_1);

        n_inv = 1.0 / nexp;
        mexp = 1.0 - n_inv;
        m_inv = -1.0 / mexp;

        /* Get the saturation and capillary pressure at switching or reversal point */
        sat_switch = pmv_hyst->sat_switch[ipore][ignode];
        cap_pres_switch = pmv_hyst->cap_pres_switch[ipore][ignode];

        /* Calculate normalized saturation at reversal point */
        brack_in_inv = pow((cap_pres_switch / con_c), nexp) + 1.0;
        sat_norm_switch = pow(brack_in_inv, (-mexp));

        /* Quality check */
        if ((sat_norm_switch > 1.0) || (sat_norm_switch < 0.0)) {
          GOMA_EH(GOMA_ERROR, "Invalid value of normalized saturation at switching point");
        }

        /* Calculate sat_min for the scanning wetting curve */
        sat_min_wet = (sat_switch - sat_max * sat_norm_switch) / (1.0 - sat_norm_switch);

        /* More quality check */
        if (sat_min_wet < sat_min) {
          GOMA_EH(GOMA_ERROR, "Invalid value of minimum saturation at the scanning wetting curve");
        }
        /* Update sat_min_wet in element storage structure*/
        pmv_hyst->sat_min_imbibe[ipore][ignode] = sat_min_wet.val();

        /* Calculate normalized saturation based on newly calculated sat_min_wet */
        sat_norm = (saturation - sat_min_wet) / (sat_max - sat_min_wet);
        d_sat_norm_d_saturation = 1.0 / (sat_max - sat_min_wet);

        brack_in = pow(sat_norm, m_inv) - 1.0;
        d_brack_in_d_sat_norm = m_inv * pow(sat_norm, (m_inv - 1.0));

        /* Evaluate capillary pressure and its sensitivities */
        cap_pres = con_c * pow(brack_in, n_inv);
      }
    } else {
      GOMA_EH(GOMA_ERROR, "To switch or not to switch; that is the question");
    }
  } else {
    GOMA_EH(GOMA_ERROR, "No models for capillary pressure");
  }
  return cap_pres;
} /* end   load_cap_pres  */

ADType ad_por_mass_source_model(void)
/******************************************************************************
 *
 *  A function which computes mass source term in porous media transport and their
 *  Jacobian sensitivities
 *
 *  Kristianto Tjiptowidjojo (June 2015)
 * Modified for AD (2025)
 *
 *
 ******************************************************************************/
{
  ADType MassSource = 0.0;

  ADType sat_min;
  ADType width, alpha;
  ADType sat_center, sat_normalized;
  ADType Hside = 0.0;
  ADType dHside_dS = 0.0;
  ADType tau = 0.0;
  ADType sink_mass_max = 0.0, nexp = 0.0;
  ADType sink_mass_clip = fv->sink_mass;

  ADType saturation = 0.0;

  /* Get saturation based on the formulation used*/

  if (pd->v[pg->imtrx][SHELL_SAT_1]) /* Saturation formulation */
  {
    saturation = ad_fv->sh_sat_1;
  } else {
    GOMA_EH(GOMA_ERROR, "Cannot find appropriate saturation model");
  }

  /*** Calculate MassSource based on the selected models ***/

  switch (mp->PorousSinkConstantsModel) {

  case LINEAR:

    if (saturation >= mp->u_porous_sink_constants[6]) {
      tau = mp->u_porous_sink_constants[0];
    } else {
      tau = 0.;
    }

    sink_mass_max = mp->u_porous_sink_constants[1];

    MassSource = -tau * mp->u_porous_sink_constants[2] * (sink_mass_max - fv->sink_mass) *
                 saturation / sink_mass_max / fv->volume_change;

    /* Again, I don't think this term belongs, contrary to the paper from which
       it came.   Disappears with application of the Reynolds Transport theorem */

    /*
    for(b=0; b < VIM; b++)
       {
        MassSource += fv_dot->x[b] * mp->density *
                      (mp->porosity * mp->d_saturation[POR_LIQ_PRES]*fv->grad_p_liq[b] +
                       mp->saturation * fv->grad_porosity[b]);
       }
    */

    break;

  case POWER_LAW:

    /* Evaluate heaviside function based on minimum saturation */
    sat_min = mp->u_porous_sink_constants[3];
    width = mp->u_porous_sink_constants[4];
    alpha = 0.5 * width;
    sat_center = sat_min - alpha;
    sat_normalized = saturation - sat_center;

    if (saturation >= sat_min) {
      Hside = 1.0;
      dHside_dS = 0.0;
    } else if (saturation <= (sat_min - width)) {
      Hside = 0.0;
      dHside_dS = 0.0;
    } else {
      Hside = 0.5 * (1. + sat_normalized / alpha + sin(M_PIE * sat_normalized / alpha) / M_PIE);
      dHside_dS = 0.5 * (1.0 / alpha + cos(M_PIE * sat_normalized / alpha) / alpha);
    }

    tau = mp->u_porous_sink_constants[0];
    sink_mass_max = mp->u_porous_sink_constants[1];
    nexp = mp->u_porous_sink_constants[2];

    if (fv->sink_mass < sink_mass_max) {
      sink_mass_clip = fv->sink_mass;
    } else {
      sink_mass_clip = sink_mass_max;
    }

    MassSource = -tau * pow(((sink_mass_max - sink_mass_clip) / sink_mass_max), nexp) * saturation /
                 mp->density;

    MassSource *= Hside;

    break;

  default:

    GOMA_EH(GOMA_ERROR, "No valid porous sink models were found");
  }

  /****** Calculate their Jacobian sensitivities */
  return (MassSource);
}

#ifdef AD_POROUS_SHELL_NODE_SOURCE
void ad_porous_shell_open_source_model(ADType j_1_2[MDE], // Flux between porous layers 1 and 2
                                       ADType j_2_3[MDE]  // Flux between porous layers 2 and 3
                                       )
/*****************************************************************************
 * This function calculates inter-layer fluxes amongst porous shell layers.
 * As of now, each layer is assumed to be stacked on top one another.
 * i.e. Layer 3 on top of layer 2 on top of layer 1.
 *
 *
 * Kristianto Tjiptowidjojo   (tjiptowi@unm.edu)  October 2019
 * Modified for AD (2025)
 *
 *****************************************************************************/
{
  int ipore, var, j;

  int porous_shell_var[MAX_POR_SHELL];
  porous_shell_var[0] = SHELL_SAT_1;
  porous_shell_var[1] = SHELL_SAT_2;
  porous_shell_var[2] = SHELL_SAT_3;

  ADType H[MAX_POR_SHELL] = {0.0};     // Pore height (vertical)
  ADType kappa[MAX_POR_SHELL] = {0.0}; // Cross permeability
  ADType mu = mp->viscosity;           // Viscosity

  ADType sat_nodes[MAX_POR_SHELL][MDE] = {{0.0}};
  ADType cap_pres_nodes[MAX_POR_SHELL][MDE] = {{0.0}};

  ADType Hside_1[MDE] = {0.0};
  ADType Hside_1_square[MDE] = {0.0};
  ADType sat_max_1 = 0.99;
  ADType width_1 = 0.24;
  ADType alpha_1 = 0.5 * width_1;
  ADType sat_center_1 = sat_max_1 - alpha_1;
  ADType sat_normalized_1[MDE] = {0.0};

  ADType Hside_2[MDE] = {0.0};
  ADType sat_min_2 = 0.5;
  ADType width_2 = 0.05;
  ADType alpha_2 = 0.5 * width_2;
  ADType sat_center_2 = sat_min_2 - alpha_2;
  ADType sat_normalized_2[MDE] = {0.0};

  ADType Hside_3[MDE] = {0.0};
  ADType Hside_3_square[MDE] = {0.0};
  ADType sat_max_3 = 0.99;
  ADType width_3 = 0.24;
  ADType alpha_3 = 0.5 * width_3;
  ADType sat_center_3 = sat_max_3 - alpha_3;
  ADType sat_normalized_3[MDE] = {0.0};

  /* Extract all of the necessary information */
  for (ipore = 0; ipore < pd->Num_Porous_Shell_Eqn; ipore++) {
    var = porous_shell_var[ipore];
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        switch (ipore) {
        case 0:
          sat_nodes[ipore][j] = set_ad_or_dbl(*esp->sh_sat_1[j], R_SHELL_SAT_1, j);
          sat_normalized_1[j] = sat_nodes[ipore][j] - sat_center_1;
          break;
        case 1:
          sat_nodes[ipore][j] = set_ad_or_dbl(*esp->sh_sat_2[j], R_SHELL_SAT_2, j);
          sat_normalized_2[j] = sat_nodes[ipore][j] - sat_center_2;
          break;
        case 2:
          sat_nodes[ipore][j] = set_ad_or_dbl(*esp->sh_sat_3[j], R_SHELL_SAT_3, j);
          sat_normalized_3[j] = sat_nodes[ipore][j] - sat_center_3;
          break;
        }
        cap_pres_nodes[ipore][j] = ad_load_cap_pres(ipore, j, -1, sat_nodes[ipore][j]);
      }
      H[ipore] = porous_shell_height_model(ipore);
      kappa[ipore] = porous_shell_cross_perm_model(ipore);
    }
  }

  /* Apply heaviside function to deactivate flux at S_1 > 0.99*/

  for (j = 0; j < ei[pg->imtrx]->dof[SHELL_SAT_1]; j++) {
    if (sat_nodes[0][j] >= sat_max_1) {
      Hside_1[j] = 0.0;
      Hside_1_square[j] = 0.0;
    } else if (sat_nodes[0][j] <= (sat_max_1 - width_1)) {
      Hside_1[j] = 1.0;
      Hside_1_square[j] = 1.0;
    } else {
      Hside_1[j] = 1.0 - 0.5 * (1. + sat_normalized_1[j] / alpha_1 -
                                sin(M_PIE * sat_normalized_1[j] / alpha_1) / M_PIE);
      Hside_1_square[j] = Hside_1[j] * Hside_1[j];
    }
  }

  /* Apply heaviside function to activate flux at S_2 > S_min*/
  for (j = 0; j < ei[pg->imtrx]->dof[SHELL_SAT_2]; j++) {
    if (sat_nodes[1][j] >= sat_min_2) {
      Hside_2[j] = 1.0;
    } else if (sat_nodes[1][j] <= (sat_min_2 - width_2)) {
      Hside_2[j] = 0.0;
    } else {
      Hside_2[j] = 0.5 * (1. + sat_normalized_2[j] / alpha_2 +
                          sin(M_PIE * sat_normalized_2[j] / alpha_2) / M_PIE);
    }
  }

  /* Apply heaviside function to deactivate flux at S_3 > 0.99*/
  if (pd->e[pg->imtrx][R_SHELL_SAT_3]) {
    for (j = 0; j < ei[pg->imtrx]->dof[SHELL_SAT_3]; j++) {
      if (sat_nodes[2][j] >= sat_max_3) {
        Hside_3[j] = 0.0;
        Hside_3_square[j] = 0.0;
      } else if (sat_nodes[2][j] <= (sat_max_3 - width_3)) {
        Hside_3[j] = 1.0;
        Hside_3_square[j] = 1.0;
      } else {
        Hside_3[j] = 1.0 - 0.5 * (1. + sat_normalized_3[j] / alpha_3 -
                                  sin(M_PIE * sat_normalized_3[j] / alpha_3) / M_PIE);
        Hside_3_square[j] = Hside_3[j] * Hside_3[j];
      }
    }
  }

  /* Populate the interporous flux */
  for (j = 0; j < ei[pg->imtrx]->dof[SHELL_SAT_1]; j++) {
    j_1_2[j] = (kappa[0] / mu) * (cap_pres_nodes[0][j] - cap_pres_nodes[1][j]) / (2.0 * H[0]);
    j_1_2[j] += (kappa[1] / mu) * (cap_pres_nodes[0][j] - cap_pres_nodes[1][j]) / (2.0 * H[1]);
    j_1_2[j] *= Hside_1_square[j] * Hside_2[j];

    if (pd->e[pg->imtrx][R_SHELL_SAT_3]) {
      j_2_3[j] = (kappa[1] / mu) * (cap_pres_nodes[2][j] - cap_pres_nodes[1][j]) / (2.0 * H[1]);
      j_2_3[j] += (kappa[2] / mu) * (cap_pres_nodes[2][j] - cap_pres_nodes[1][j]) / (2.0 * H[2]);
      j_2_3[j] *= Hside_2[j] * Hside_3_square[j];
    }
  }

  /* Populate the interporous flux sensitivity */
  return;
}
#else
void ad_porous_shell_open_source_model(ADType j_1_2, // Flux between porous layers 1 and 2
                                       ADType j_2_3, // Flux between porous layers 2 and 3
                                       ADType cap_pres_fv[MAX_POR_SHELL])
/*****************************************************************************
 * This function calculates inter-layer fluxes amongst porous shell layers.
 * As of now, each layer is assumed to be stacked on top one another.
 * i.e. Layer 3 on top of layer 2 on top of layer 1.
 *
 *
 * Kristianto Tjiptowidjojo   (tjiptowi@unm.edu)  October 2019
 * Modified for AD (2025)
 *
 *****************************************************************************/
{
  int ipore, var;

  int porous_shell_var[MAX_POR_SHELL];
  porous_shell_var[0] = SHELL_SAT_1;
  porous_shell_var[1] = SHELL_SAT_2;
  porous_shell_var[2] = SHELL_SAT_3;

  ADType H[MAX_POR_SHELL] = {0.0};     // Pore height (vertical)
  ADType kappa[MAX_POR_SHELL] = {0.0}; // Cross permeability
  ADType mu = mp->viscosity;           // Viscosity

  ADType sat_nodes[MAX_POR_SHELL] = {{0.0}};
  ADType cap_pres_nodes[MAX_POR_SHELL] = {{0.0}};

  ADType Hside_1 = 0;
  ADType Hside_1_square = 0;
  ADType sat_max_1 = 0.99;
  ADType width_1 = 0.24;
  ADType alpha_1 = 0.5 * width_1;
  ADType sat_center_1 = sat_max_1 - alpha_1;
  ADType sat_normalized_1 = 0;

  ADType Hside_2 = 0;
  ADType sat_min_2 = 0.5;
  ADType width_2 = 0.05;
  ADType alpha_2 = 0.5 * width_2;
  ADType sat_center_2 = sat_min_2 - alpha_2;
  ADType sat_normalized_2 = 0;

  ADType Hside_3 = 0;
  ADType Hside_3_square = 0;
  ADType sat_max_3 = 0.99;
  ADType width_3 = 0.24;
  ADType alpha_3 = 0.5 * width_3;
  ADType sat_center_3 = sat_max_3 - alpha_3;
  ADType sat_normalized_3 = 0;

  /* Extract all of the necessary information */
  for (ipore = 0; ipore < pd->Num_Porous_Shell_Eqn; ipore++) {
    var = porous_shell_var[ipore];
    if (pd->v[pg->imtrx][var]) {
      switch (ipore) {
      case 0:
        sat_normalized_1 = sat_nodes[ipore] - sat_center_1;
        sat_nodes[ipore] = ad_fv->sh_sat_3;
        break;
      case 1:
        sat_nodes[ipore] = ad_fv->sh_sat_2;
        sat_normalized_2 = sat_nodes[ipore] - sat_center_2;
        break;
      case 2:
        sat_nodes[ipore] = ad_fv->sh_sat_3;
        sat_normalized_3 = sat_nodes[ipore] - sat_center_3;
        break;
      }
    }
    H[ipore] = porous_shell_height_model(ipore);
    kappa[ipore] = porous_shell_cross_perm_model(ipore);
  }

  /* Apply heaviside function to deactivate flux at S_1 > 0.99*/

  if (sat_nodes[0] >= sat_max_1) {
    Hside_1 = 0.0;
    Hside_1_square = 0.0;
  } else if (sat_nodes[0] <= (sat_max_1 - width_1)) {
    Hside_1 = 1.0;
    Hside_1_square = 1.0;
  } else {
    Hside_1 = 1.0 - 0.5 * (1. + sat_normalized_1 / alpha_1 -
                           sin(M_PIE * sat_normalized_1 / alpha_1) / M_PIE);
    Hside_1_square = Hside_1 * Hside_1;
  }

  /* Apply heaviside function to activate flux at S_2 > S_min*/
  if (sat_nodes[1] >= sat_min_2) {
    Hside_2 = 1.0;
  } else if (sat_nodes[1] <= (sat_min_2 - width_2)) {
    Hside_2 = 0.0;
  } else {
    Hside_2 =
        0.5 * (1. + sat_normalized_2 / alpha_2 + sin(M_PIE * sat_normalized_2 / alpha_2) / M_PIE);
  }

  /* Apply heaviside function to deactivate flux at S_3 > 0.99*/
  if (pd->e[pg->imtrx][R_SHELL_SAT_3]) {
    if (sat_nodes[2] >= sat_max_3) {
      Hside_3 = 0.0;
      Hside_3_square = 0.0;
    } else if (sat_nodes[2] <= (sat_max_3 - width_3)) {
      Hside_3 = 1.0;
      Hside_3_square = 1.0;
    } else {
      Hside_3 = 1.0 - 0.5 * (1. + sat_normalized_3 / alpha_3 -
                             sin(M_PIE * sat_normalized_3 / alpha_3) / M_PIE);
      Hside_3_square = Hside_3 * Hside_3;
    }
  }

  /* Populate the interporous flux */
  j_1_2 = (kappa[0] / mu) * (cap_pres_fv[0] - cap_pres_fv[1]) / (2.0 * H[0]);
  j_1_2 += (kappa[1] / mu) * (cap_pres_fv[0] - cap_pres_fv[1]) / (2.0 * H[1]);
  j_1_2 *= Hside_1_square * Hside_2;

  if (pd->e[pg->imtrx][R_SHELL_SAT_3]) {
    j_2_3 = (kappa[1] / mu) * (cap_pres_fv[2] - cap_pres_fv[1]) / (2.0 * H[1]);
    j_2_3 += (kappa[2] / mu) * (cap_pres_fv[2] - cap_pres_fv[1]) / (2.0 * H[2]);
    j_2_3 *= Hside_2 * Hside_3_square;
  }

  /* Populate the interporous flux sensitivity */
  return;
}

#endif

ADType porous_shell_rel_perm_model(int ipore, ADType saturation) {
  /******************************************************************************
   *
   *  This function computes the relative permeability of an open porous shell.
   *  Used with the function assemble_porous_shell_saturation.
   *
   *  Kristianto Tjiptowidjojo (tjiptowi@unm.edu) - October 2018
   *
   ******************************************************************************/
  ADType k_rel = 0.0;
  ADType s_eff;
  ADType d_s_eff, a1, factor, factor2;
  ADType expon2, sat_min, sat_max, viscosity, lambda;
  int i_rel_perm_ev;
  double scale;

  switch (mp->PorousShellRelPermModel[ipore]) {

  case CONSTANT:

    k_rel = mp->PorousShellRelPerm[ipore];

    break;

  case VAN_GENUCHTEN:

    /*
     *
     * FOR VAN_GENUCHTEN EQUATION
     *  mp->u_PorousShellRelPerm[ipore][0] is the irreduceable water saturation
     *  mp->u_PorousShellRelPerm[ipore][1] is the irreduceable air saturation
     *  mp->u_PorousShellRelPerm[ipore][2] is the exponent, 1 - 1/beta
     *  mp->u_PorousShellRelPerm[ipore][3] is the liquid viscosity
     *
     *  Store some temporary variables
     */

    sat_min = mp->u_PorousShellRelPerm[ipore][0];
    sat_max = 1.0 - mp->u_PorousShellRelPerm[ipore][1];
    s_eff = (saturation - sat_min) / (sat_max - sat_min);
    viscosity = mp->u_PorousShellRelPerm[ipore][3];
    lambda = mp->u_PorousShellRelPerm[ipore][2];

    /*
     *  Clip the relative permeability to zero if the effective saturation
     *  is equal to or less than zero. -> there can be no transport
     *  in a liquid phase if there is no continguous pathway in that phase.
     */
    if (s_eff < 0.0) {
      k_rel = mp->PorousShellRelPerm[ipore] = 0.0;
    }

    /*
     *  Clip the relative permeability at one -> it can never be
     *  greater than one.  Actually, note that if somehow s_eff is
     *  very close to 1.0 and fails this test, then you are dividing
     *  by zero as factor=1.0 below.    Now and then GOMA aborts due
     *  to this.
     */
    else if (s_eff >= 0.99999) {
      k_rel = 1.0 / viscosity;
      mp->PorousShellRelPerm[ipore] = k_rel.val();

    }

    /*
     *  Otherwise, apply Van Genuchten formula
     */
    else {
      expon2 = 1.0 / lambda;
      factor = pow(s_eff, expon2);
      factor2 = pow(1.0 - factor, lambda);
      a1 = 1.0 - factor2;
      k_rel = sqrt(s_eff) * a1 * a1 / viscosity;

      mp->PorousShellRelPerm[ipore] = k_rel.val();
    }

    break;

  case EXTERNAL_FIELD:

    /*
     *
     * FOR EXTERNAL FIELD
     *  mp->u_PorousShellRelPerm[ipore][0] is the scaling factor for the read-in external field
     * variable
     */

    i_rel_perm_ev = mp->por_shell_rel_perm_ext_field_index[ipore];
    scale = mp->u_rel_liq_perm[0];

    k_rel = mp->PorousShellRelPerm[ipore] = scale * fv->external_field[i_rel_perm_ev];

    break;

  case VAN_GENUCHTEN_EXTERNAL:

    /*
     *
     * FOR VAN_GENUCHTEN_EXTERNAL
     *  mp->u_PorousShellRelPerm[ipore][0] is the irreduceable water saturation
     *  mp->u_PorousShellRelPerm[ipore][1] is the irreduceable air saturation
     *  mp->u_PorousShellRelPerm[ipore][2] is the exponent, 1 - 1/beta for external field value of 0
     *  mp->u_PorousShellRelPerm[ipore][3] is the liquid viscosity
     *  mp->u_PorousShellRelPerm[ipore][4] is the exponent, 1 - 1/beta for external field value of 1
     */

    i_rel_perm_ev = mp->por_shell_rel_perm_ext_field_index[ipore];

    sat_min = mp->u_PorousShellRelPerm[ipore][0];
    sat_max = 1.0 - mp->u_PorousShellRelPerm[ipore][1];
    s_eff = (saturation - sat_min) / (sat_max - sat_min);
    viscosity = mp->u_PorousShellRelPerm[ipore][3];

    /* Here I assume that efv is bounded between 0 and 1 */
    lambda = fv->external_field[i_rel_perm_ev] *
                 (mp->u_PorousShellRelPerm[ipore][4] - mp->u_PorousShellRelPerm[ipore][2]) +
             mp->u_PorousShellRelPerm[ipore][2];

    /*
     *  Clip the relative permeability to zero if the effective saturation
     *  is equal to or less than zero. -> there can be no transport
     *  in a liquid phase if there is no continguous pathway in that phase.
     */
    if (s_eff < 0.0) {
      k_rel = mp->PorousShellRelPerm[ipore] = 0.0;
    }

    /*
     *  Clip the relative permeability at one -> it can never be
     *  greater than one.  Actually, note that if somehow s_eff is
     *  very close to 1.0 and fails this test, then you are dividing
     *  by zero as factor=1.0 below.    Now and then GOMA aborts due
     *  to this.
     */
    else if (s_eff >= 0.99999) {
      k_rel = 1.0 / viscosity;
      mp->PorousShellRelPerm[ipore] = k_rel.val();
    }

    /*
     *  Otherwise, apply Van Genuchten formula
     */
    else {
      expon2 = 1.0 / lambda;
      factor = pow(s_eff, expon2);
      factor2 = pow(1.0 - factor, lambda);
      a1 = 1.0 - factor2;
      k_rel = sqrt(s_eff) * a1 * a1 / viscosity;

      mp->PorousShellRelPerm[ipore] = k_rel.val();
    }

    break;

  default:
    GOMA_EH(GOMA_ERROR, "Unrecognized Porous Shell Relative Permeability model");
    break;
  }

  return (k_rel);
}

void Inn_ad(ADType v[DIM], // Input vector
            ADType w[DIM]  // Output rotated vector
            )
/******************************************************************************
 *
 * Inn()
 *
 * Function to rotate into shell cordinates by the following transformation:
 *       w = (I-nn)*v
 *
 * Scott A Roberts (1514) sarober@sandia.gov
 *
 ******************************************************************************/
{
  int i, j;
  for (i = 0; i < DIM; i++) {
    w[i] = 0.0;
    for (j = 0; j < DIM; j++) {
      w[i] += (v[j] * delta(i, j) - v[j] * fv->snormal[i] * fv->snormal[j]);
    }
  }
  return;
} /* End of Inn */

extern "C" int ad_assemble_porous_shell_saturation(dbl tt,           // Time integration form
                                                   dbl dt,           // Time step size
                                                   dbl xi[DIM],      // Current coordinates
                                                   const Exo_DB *exo // ExoII handle
) {

  /* --- Initialization -----------------------------------------------------*/

  // Initialize output status function
  int status = 0;

  // Bail out of function if there is nothing to do
  if ((!pd->e[pg->imtrx][R_SHELL_SAT_1]) && (!pd->e[pg->imtrx][R_SHELL_SAT_2]) &&
      (!pd->e[pg->imtrx][R_SHELL_SAT_3]))
    return (status);

  // Variable definitions
  int eqn, peqn, var;                            // Equation / variables
  int i, j, a, b, ipore;                         // Counter variables
  dbl phi_i, grad_phi_i[DIM], gradII_phi_i[DIM]; // Basis functions (i)
  dbl d_gradII_phi_i_dmesh[DIM][DIM][MDE];

  // Bookkeeping arrays
  int porous_shell_eqn[MAX_POR_SHELL];
  porous_shell_eqn[0] = R_SHELL_SAT_1;
  porous_shell_eqn[1] = R_SHELL_SAT_2;
  porous_shell_eqn[2] = R_SHELL_SAT_3;

  int porous_shell_var[MAX_POR_SHELL];
  porous_shell_var[0] = SHELL_SAT_1;
  porous_shell_var[1] = SHELL_SAT_2;
  porous_shell_var[2] = SHELL_SAT_3;

  // Group saturation values at nodes and Gauss point
  ADType sat_nodes[MAX_POR_SHELL][MDE] = {{0.0}};
  ADType sat_dot_nodes[MAX_POR_SHELL][MDE] = {{0.0}};
  ADType sat_gauss[MAX_POR_SHELL] = {0.0};
  ADType grad_sat_gauss[MAX_POR_SHELL][DIM] = {{0.0}};
  ADType grad_II_sat_gauss[MAX_POR_SHELL][DIM] = {{0.0}};
  ADType sat_dot_gauss[MAX_POR_SHELL] = {0.0};
  for (ipore = 0; ipore < pd->Num_Porous_Shell_Eqn; ipore++) {
    var = porous_shell_var[ipore];
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        switch (ipore) {
        case 0: {
          sat_nodes[ipore][j] = set_ad_or_dbl(*esp->sh_sat_1[j], R_SHELL_SAT_1, j);
          ADType udot = set_ad_or_dbl(*esp_dot->sh_sat_1[j], R_SHELL_SAT_1, j);
          if (af->Assemble_Jacobian == TRUE) {
            udot.fastAccessDx(ad_fv->offset[R_SHELL_SAT_1] + j) =
                (1. + 2. * tran->current_theta) / tran->delta_t;
          }
          sat_dot_nodes[ipore][j] = udot;
        } break;
        case 1: {
          sat_nodes[ipore][j] = set_ad_or_dbl(*esp->sh_sat_2[j], R_SHELL_SAT_2, j);
          ADType udot = set_ad_or_dbl(*esp_dot->sh_sat_2[j], R_SHELL_SAT_2, j);
          if (af->Assemble_Jacobian == TRUE) {
            udot.fastAccessDx(ad_fv->offset[R_SHELL_SAT_2] + j) =
                (1. + 2. * tran->current_theta) / tran->delta_t;
          }
          sat_dot_nodes[ipore][j] = udot;
        } break;
        case 2: {
          sat_nodes[ipore][j] = set_ad_or_dbl(*esp->sh_sat_3[j], R_SHELL_SAT_3, j);
          ADType udot = set_ad_or_dbl(*esp_dot->sh_sat_3[j], R_SHELL_SAT_3, j);
          if (af->Assemble_Jacobian == TRUE) {
            udot.fastAccessDx(ad_fv->offset[R_SHELL_SAT_3] + j) =
                (1. + 2. * tran->current_theta) / tran->delta_t;
          }
          sat_dot_nodes[ipore][j] = udot;
        } break;
        }
      }
      switch (ipore) {
      case 0:
        sat_gauss[ipore] = ad_fv->sh_sat_1;
        sat_dot_gauss[ipore] = ad_fv->sh_sat_1_dot;
        break;
      case 1:
        sat_gauss[ipore] = ad_fv->sh_sat_2;
        sat_dot_gauss[ipore] = ad_fv->sh_sat_2_dot;
        break;
      case 2:
        sat_gauss[ipore] = ad_fv->sh_sat_3;
        sat_dot_gauss[ipore] = ad_fv->sh_sat_3_dot;
        break;
      }
      for (a = 0; a < DIM; a++) {
        switch (ipore) {
        case 0:
          grad_sat_gauss[ipore][a] = ad_fv->grad_sh_sat_1[a];
          break;
        case 1:
          grad_sat_gauss[ipore][a] = ad_fv->grad_sh_sat_2[a];
          break;
        case 2:
          grad_sat_gauss[ipore][a] = ad_fv->grad_sh_sat_3[a];
          break;
        }
      }
      Inn_ad(grad_sat_gauss[ipore], grad_II_sat_gauss[ipore]);
    }
  }

  // Setup lubrication
  int *n_dof = NULL;
  int dof_map[MDE];
  dbl wt_old = fv->wt;
  n_dof = (int *)array_alloc(1, MAX_VARIABLE_TYPES, sizeof(int));
  lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

  // Unpack FEM variables from global structures
  dbl wt = fv->wt;      // Gauss point weight
  dbl h3 = fv->h3;      // Differential volume element
  dbl det_J = fv->sdet; // Jacobian of transformation
  dbl dA = det_J * wt * h3;

  /* Equation term multipliers*/
  dbl etm_mass[MAX_POR_SHELL] = {0.0}, etm_diff[MAX_POR_SHELL] = {0.0},
      etm_sour[MAX_POR_SHELL] = {0.0};
  for (ipore = 0; ipore < pd->Num_Porous_Shell_Eqn; ipore++) {
    eqn = porous_shell_eqn[ipore];
    if (pd->e[pg->imtrx][eqn]) {
      etm_mass[ipore] = pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
      etm_diff[ipore] = pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
      etm_sour[ipore] = pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
    }
  }

  // Load porous medium parameters
  ADType phi[MAX_POR_SHELL] = {0.0}; // Porosity
  ADType H[MAX_POR_SHELL] = {0.0};   // Pore height (vertical)
  for (ipore = 0; ipore < pd->Num_Porous_Shell_Eqn; ipore++) {
    var = porous_shell_var[ipore];
    if (pd->v[pg->imtrx][var]) {
      phi[ipore] = porous_shell_porosity_model(ipore);
      H[ipore] = porous_shell_height_model(ipore);
      porous_shell_permeability_model(ipore);
    }
  }

  /* --- Calculate equation components ---------------------------------------*/

  /* With saturation formulation, capillary pressure has to be evaluated from
   * capillary pressure - saturation curve. This is problematic in computing
   * capillary pressure gradient in diffusion term viz. having to evaluate
   * second derivative of capillary pressure w.r.t. saturation, i.e. d^2_cap_pres/d_S^2
   * to evaluate for Jacobian.
   *
   * I propose to do this instead: Use the curve to evaluate capillary pressure at the nodes
   * then use basis functions of the saturation DOF to compute gradient at the Gauss point
   * That way we only require first derivative to compute the Jacobian.
   *
   */
  ADType cap_pres[MAX_POR_SHELL][MDE] = {{0.0}};
  ADType cap_pres_fv[MAX_POR_SHELL] = {0.0};
  ADType grad_p[MAX_POR_SHELL][DIM] = {{0.0}};
  ADType grad_II_p[MAX_POR_SHELL][DIM] = {{0.0}};
  for (ipore = 0; ipore < pd->Num_Porous_Shell_Eqn; ipore++) {
    var = porous_shell_var[ipore];
    if (pd->v[pg->imtrx][var]) {
      /* Old method */
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        cap_pres[ipore][j] = ad_load_cap_pres(ipore, j, -1, sat_nodes[ipore][j]);
      }
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        cap_pres_fv[ipore] += cap_pres[ipore][j] * bf[var]->phi[j];
      }
      for (a = 0; a < DIM; a++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          grad_p[ipore][a] += cap_pres[ipore][j] * bf[var]->grad_phi[j][a];
        }
      }
      Inn_ad(grad_p[ipore], grad_II_p[ipore]);
    }
  }

  /* Assemble each component of the equation */

  ADType E_MASS[MAX_POR_SHELL][MDE] = {{0.0}};

  int mass_lump = 1;

  for (ipore = 0; ipore < pd->Num_Porous_Shell_Eqn; ipore++) {
    var = porous_shell_var[ipore];
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        if (mass_lump) {
          E_MASS[ipore][j] = H[ipore] * phi[ipore] * sat_dot_nodes[ipore][j];
        } else {
          E_MASS[ipore][j] = H[ipore] * phi[ipore] * sat_dot_gauss[ipore];
        }
      }
    }
  }

  // Load relative permeability as a function of saturation
  ADType rel_liq_perm[MAX_POR_SHELL] = {0.0};

  for (ipore = 0; ipore < pd->Num_Porous_Shell_Eqn; ipore++) {
    var = porous_shell_var[ipore];
    if (pd->v[pg->imtrx][var]) {
      rel_liq_perm[ipore] = porous_shell_rel_perm_model(ipore, sat_gauss[ipore]);
    }
  }

  // Calculate DIFFUSION terms
  ADType E_DIFF[MAX_POR_SHELL][DIM] = {{0.0}};

  for (ipore = 0; ipore < pd->Num_Porous_Shell_Eqn; ipore++) {
    var = porous_shell_var[ipore];
    if (pd->v[pg->imtrx][var]) {
      if (mp->PorousShellPermeabilityModel[ipore] == ORTHOTROPIC) {
        for (a = 0; a < DIM; a++) {
          for (b = 0; b < DIM; b++) {
            E_DIFF[ipore][a] += H[ipore] * mp->PorousShellPermTensor[ipore][a][b] *
                                rel_liq_perm[ipore] *
                                (grad_II_p[ipore][b] + mp->momentum_source[b]);
          }
        }
      } else {
        for (a = 0; a < DIM; a++) {
          E_DIFF[ipore][a] += H[ipore] * mp->PorousShellPermeability[ipore] * rel_liq_perm[ipore] *
                              (grad_II_p[ipore][a] + mp->momentum_source[a]);
        }
      }
    }
  }

  // Calculate SOURCE term from adjacent porous shells
  ADType E_SOUR[MAX_POR_SHELL][MDE] = {{0.0}};

  if (pd->Num_Porous_Shell_Eqn > 1) {

#ifdef AD_POROUS_SHELL_NODE_SOURCE
    ADType j_1_2[MDE] = {0.0}; // Flux between porous layers 1 and 2
    ADType j_2_3[MDE] = {0.0}; // Flux between porous layers 2 and 3

    /* Calculate the fluxes and their sensitivities */
    ad_porous_shell_open_source_model(j_1_2, j_2_3);

    /* Store the interlayer fluxes and the sensitivities */

    for (j = 0; j < ei[pg->imtrx]->dof[SHELL_SAT_1]; j++) {
      E_SOUR[0][j] = -j_1_2[j];
      E_SOUR[1][j] = j_1_2[j];

      if (pd->Num_Porous_Shell_Eqn > 2) {
        E_SOUR[1][j] += j_2_3[j];
        E_SOUR[2][j] = -j_2_3[j];
      }
    }
#else

    ADType j_1_2 = 0;
    ADType j_2_3 = 0;

    /* Calculate the fluxes and their sensitivities */
    ad_porous_shell_open_source_model(j_1_2, j_2_3, cap_pres_fv);

    /* Store the interlayer fluxes and the sensitivities */

    for (j = 0; j < ei[pg->imtrx]->dof[SHELL_SAT_1]; j++) {
      E_SOUR[0][j] = -j_1_2;
      E_SOUR[1][j] = j_1_2;

      if (pd->Num_Porous_Shell_Eqn > 2) {
        E_SOUR[1][j] += j_2_3;
        E_SOUR[2][j] = -j_2_3;
      }
    }

#endif
  }

  // Load sink terms due to adsorption and its sensitivities
  // Right now it only applies to first porous shell layer - SHELL_SAT_1
  ADType E_SINK = 0.0;

  if (pd->e[pg->imtrx][R_POR_SINK_MASS]) {
    E_SINK = ad_por_mass_source_model();
  }

  /* --- Assemble residuals --------------------------------------------------*/

  // Assemble residual contribution to this equation

  /* Loop over porous shell layer equations*/
  std::vector<std::vector<ADType>> resid(pd->Num_Porous_Shell_Eqn);

  for (ipore = 0; ipore < pd->Num_Porous_Shell_Eqn; ipore++) {
    eqn = porous_shell_eqn[ipore];
    resid[ipore].resize(ei[pg->imtrx]->dof[eqn]);
    for (int i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      resid[ipore][i] = 0;
    }
  }

  for (ipore = 0; ipore < pd->Num_Porous_Shell_Eqn; ipore++) {
    eqn = porous_shell_eqn[ipore];
    if ((af->Assemble_Residual) && pd->e[pg->imtrx][eqn]) {
      peqn = upd->ep[pg->imtrx][eqn];

      // Loop over DOF (i)
      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

        // Load basis functions
        ShellBF(eqn, i, &phi_i, grad_phi_i, gradII_phi_i, d_gradII_phi_i_dmesh,
                n_dof[MESH_DISPLACEMENT1], dof_map);

        // Assemble mass term
        ADType mass = 0.0;
        if (pd->e[pg->imtrx][eqn] & T_MASS) {
          mass += E_MASS[ipore][i] * phi_i;
        }
        mass *= dA * etm_mass[ipore];

        // Assemble diffusion term
        ADType diff = 0.0;
        if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
          for (a = 0; a < DIM; a++) {
            diff -= E_DIFF[ipore][a] * gradII_phi_i[a];
          }
        }
        diff *= dA * etm_diff[ipore];

        // Assemble source term
        ADType sour = 0.0;
        if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
          sour += E_SOUR[ipore][0] * phi_i;
          if (ipore == 0) {
            sour -= E_SINK * phi_i;
          }
        }
        sour *= dA * etm_sour[ipore];

        // Assemble full residual
        resid[ipore][i] += mass + diff + sour;
        lec->R[LEC_R_INDEX(peqn, i)] += mass.val() + diff.val() + sour.val();
      } // End of loop over DOF (i)
    } // End of residual assembly of R_SHELL_SAT_1
  } // End of loop over porous shell layers

  /* --- Assemble Jacobian --------------------------------------------------*/

  /* Loop over porous shell layer equations*/
  for (ipore = 0; ipore < pd->Num_Porous_Shell_Eqn; ipore++) {
    eqn = porous_shell_eqn[ipore];
    if ((af->Assemble_Jacobian) && (pd->e[pg->imtrx][eqn])) {
      peqn = upd->ep[pg->imtrx][eqn];

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        for (int var = V_FIRST; var < V_LAST; var++) {

          /* Sensitivity w.r.t. velocity */
          if (pd->v[pg->imtrx][var]) {
            int pvar = upd->vp[pg->imtrx][var];

            for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              // J = &(lec->J[LEC_J_INDEX(peqn, pvar, ii, 0)]);
              lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += resid[ipore][i].dx(ad_fv->offset[var] + j);

            } /* End of loop over j */
          } /* End of if the variale is active */
        }
      }
    } /* End of if assemble Jacobian */

  } // End of loop over porous shell layers

  // Finalize and exit
  fv->wt = wt_old;
  safe_free((void *)n_dof);
  return (status);

} // End of assemble_porous_shell_saturation
#endif // GOMA_ENABLE_SACADO