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
 *$Id: mm_std_models_shell.c,v 5.31 2010-07-30 20:48:38 prschun Exp $
 */

#ifdef USE_RCSID
static char rcsid[] = "$Id: mm_std_models_shell.c,v 5.31 2010-07-30 20:48:38 prschun Exp $";
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

/* GOMA include files */

#define GOMA_MM_STD_MODELS_SHELL_C
#include "el_elm.h"
#include "el_geom.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_mp_const.h"
#include "rf_bc_const.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_fill_const.h"
#include "rf_io.h"
#include "rf_io_const.h"
#include "rf_io_structs.h"
#include "rf_masks.h"
#include "rf_mp.h"
#include "rf_solver_const.h"
#include "rf_vars_const.h"
#include "std.h"

#include "mm_mp.h"
#include "mm_mp_structs.h"

#include "mm_eh.h"

#include "mm_fill_porous.h"
#include "mm_shell_util.h"
#include "mm_std_models_shell.h"

/*********** R O U T I N E S   I N   T H I S   F I L E ************************
 *
 *       NAME            TYPE            CALLED_BY
 *    ------------             ---------               --------------
 * .....
 *
 ******************************************************************************/
/*
 * This file contains most implemented models for material properties and equation sources
 * for SHELL material types.  It is an exact analogy for mm_std_models.c for continuum.
 ******************************************************************************/

/*****************************************************************************/
double height_function_model(double *H_U,
                             double *dH_U_dtime,
                             double *H_L,
                             double *dH_L_dtime,
                             double dH_U_dX[DIM],
                             double dH_L_dX[DIM],
                             double *dH_U_dp,
                             double *dH_U_ddh,
                             double time,    /* present time value           */
                             double delta_t) /* present time step             */

/******************************************************************************
 *
 *  A function which computes the height at the current time and the rate-of-change of
 *  height.  This model is used for the lubrication capability
 *
 *  P. Randall Schunk (March 2009, Somewhere over Texas)
 *
 *
 ******************************************************************************/

{
  double H = 0.;
  double H_dot, H_init, H_delta, H_low, x_0, Length, r;
  double R, origin[3], dir_angle[3], t, axis_pt[3], dist, cos_denom;

  dH_L_dX[0] = 0.0;
  dH_L_dX[1] = 0.0;
  dH_L_dX[2] = 0.0;

  if (pd->TimeIntegration == STEADY)
    time = 0.;

  if (mp->HeightUFunctionModel == CONSTANT) {
    *H_U = mp->heightU;
    *dH_U_dtime = 0.0;
    *dH_U_dp = 0.0;
    *dH_U_ddh = 0.0;
    dH_U_dX[0] = dH_U_dX[1] = dH_U_dX[2] = 0.0;
  }

  else if (mp->HeightUFunctionModel == CONSTANT_SPEED) {
    H_dot = mp->u_heightU_function_constants[0];
    H_init = mp->u_heightU_function_constants[1];
    *H_U = H_dot * time + H_init;

    /* add on external field height if there is one. The scale factor will be the third user const*/
    if (mp->heightU_ext_field_index >= 0)
      *H_U += mp->u_heightU_function_constants[2] * fv->external_field[mp->heightU_ext_field_index];

    *dH_U_dtime = H_dot;
    dH_U_dX[0] = dH_U_dX[1] = dH_U_dX[2] = 0.;
    *dH_U_dp = 0.;
    *dH_U_ddh = 0.0;
  }

  else if (mp->HeightUFunctionModel == EXTERNAL_FIELD) {
    /* Note that this isthe same model as CONSTANT_SPEED, but allows for more generality on input */
    H_dot = mp->u_heightU_function_constants[0];
    H_init = mp->u_heightU_function_constants[1];
    *H_U = H_dot * time + H_init;

    *H_U += mp->u_heightU_function_constants[2] * fv->external_field[mp->heightU_ext_field_index];

    *dH_U_dtime = H_dot;
    dH_U_dX[0] = dH_U_dX[1] = dH_U_dX[2] = 0.;
    *dH_U_dp = 0.;
    *dH_U_ddh = 0.0;
  }

  else if (mp->HeightUFunctionModel == CONSTANT_SPEED_DEFORM) {
    double E_mod, L_0, Pext;
    H_dot = mp->u_heightU_function_constants[0];
    H_init = mp->u_heightU_function_constants[1];
    E_mod = mp->u_heightU_function_constants[2];
    L_0 = mp->u_heightU_function_constants[3];
    Pext = mp->u_heightU_function_constants[4];

    *H_U = H_dot * time + H_init + (fv->lubp - Pext) / (E_mod / L_0);
    // Right now this isn't complete because we need an augmenting condition
    // to give us a fv_dot->lubp kicker.
    *dH_U_dtime = H_dot + fv_dot->lubp / (E_mod / L_0);
    *dH_U_dp = 1. / (E_mod / L_0);
    *dH_U_ddh = 0.0;
    dH_U_dX[0] = dH_U_dX[1] = dH_U_dX[2] = 0.;
  } else if (mp->HeightUFunctionModel == CONSTANT_SPEED_MELT) {
    H_dot = mp->u_heightU_function_constants[0];
    H_init = mp->u_heightU_function_constants[1];

    *H_U = H_dot * time + H_init + fv->sh_dh; // 0.01*fv->external_field[0];

    // PRS: I'll leave this in here for now
    // because this is where you would add
    // a volume expansion effect upon phase change.
    // otherwise it is 0.*fv_dot->sh_dh

    *dH_U_dtime = H_dot - 0. * fv_dot->sh_dh;
    *dH_U_dp = 0.;
    *dH_U_ddh = 1.0;
    dH_U_dX[0] = dH_U_dX[1] = dH_U_dX[2] = 0.;

    if (*H_U <= 0.00001 * H_init) {
      *H_U = 0.00001 * H_init;
    }

  } else if (mp->HeightUFunctionModel == ROLL_ON) {
    Length = mp->u_heightU_function_constants[4];
    H_dot = mp->u_heightU_function_constants[3];
    H_delta = mp->u_heightU_function_constants[2];
    x_0 = mp->u_heightU_function_constants[0];
    H_low = mp->u_heightU_function_constants[1];

    *H_U = (H_dot * time + H_delta) * ((fv->x[0] - x_0) / Length) + H_low;

    /* add on external field height if there is one. The scale factor will be the third user const*/
    if (mp->heightU_ext_field_index >= 0)
      *H_U += mp->u_heightU_function_constants[5] * fv->external_field[mp->heightU_ext_field_index];

    *dH_U_dtime = H_dot * (fv->x[0] - x_0) / Length;
    dH_U_dX[0] = (H_dot * time + H_delta) / Length;
    dH_U_dX[1] = 0.;
    dH_U_dX[2] = 0.;
    *dH_U_dp = 0.;
    *dH_U_ddh = 0.0;

  } else if (mp->HeightUFunctionModel == ROLL_ON_MELT) {
    Length = mp->u_heightU_function_constants[4];
    H_dot = mp->u_heightU_function_constants[3];
    H_delta = mp->u_heightU_function_constants[2];
    x_0 = mp->u_heightU_function_constants[0];
    H_low = mp->u_heightU_function_constants[1];

    *H_U = (H_dot * time + H_delta) * ((fv->x[0] - x_0) / Length) + H_low + fv->sh_dh;
    *dH_U_dtime = H_dot * (fv->x[0] - x_0) / Length;
    dH_U_dX[0] = (H_dot * time + H_delta) / Length;
    dH_U_dX[1] = 0.;
    dH_U_dX[2] = 0.;
    *dH_U_dp = 0.;
    *dH_U_ddh = 1.0;

  } else if (mp->HeightUFunctionModel == ROLL) {
    R = mp->u_heightU_function_constants[0];
    /*  origin and direction of rotation axis	*/
    origin[0] = mp->u_heightU_function_constants[1];
    origin[1] = mp->u_heightU_function_constants[2];
    origin[2] = mp->u_heightU_function_constants[3];
    dir_angle[0] = mp->u_heightU_function_constants[4];
    dir_angle[1] = mp->u_heightU_function_constants[5];
    dir_angle[2] = mp->u_heightU_function_constants[6];
    H_dot = mp->u_heightU_function_constants[7];
    origin[2] += H_dot * time;
    *dH_U_dp = 0.;
    *dH_U_ddh = 0.0;

    /*  find intersection of axis with normal plane - i.e., locate point on
            axis that intersects plane normal to axis that contains local point. */

    cos_denom = (SQUARE(dir_angle[0]) + SQUARE(dir_angle[1]) + SQUARE(dir_angle[2]));
    t = (dir_angle[0] * (fv->x[0] - origin[0]) + dir_angle[1] * (fv->x[1] - origin[1]) +
         dir_angle[2] * (fv->x[2] - origin[2])) /
        cos_denom;
    axis_pt[0] = origin[0] + dir_angle[0] * t;
    axis_pt[1] = origin[1] + dir_angle[1] * t;
    axis_pt[2] = origin[2] + dir_angle[2] * t;

    /*  compute radial direction	*/

    dist = sqrt(SQUARE(fv->x[0] - axis_pt[0]) + SQUARE(fv->x[1] - axis_pt[1]));
    if (dist > fabs(R)) {
      *H_U = fabs(R);
      dH_U_dX[0] = dH_U_dX[1] = dH_U_dX[2] = 0.0;
    } else {
      *H_U = SGN(R) * (axis_pt[2] - fv->x[2] - sqrt(SQUARE(R) - SQUARE(dist)));
      *dH_U_dtime = 0.; /* finish later  */
      dH_U_dX[0] = ((fv->x[0] - axis_pt[0]) * (1. - SQUARE(dir_angle[0]) / cos_denom) +
                    (fv->x[1] - axis_pt[1]) * (-dir_angle[0] * dir_angle[1] / cos_denom) +
                    (fv->x[2] - axis_pt[2]) * (-dir_angle[0] * dir_angle[2] / cos_denom)) /
                   sqrt(SQUARE(R) - SQUARE(dist));
      dH_U_dX[1] = ((fv->x[0] - axis_pt[0]) * (-dir_angle[1] * dir_angle[0] / cos_denom) +
                    (fv->x[1] - axis_pt[1]) * (1. - SQUARE(dir_angle[1]) / cos_denom) +
                    (fv->x[2] - axis_pt[2]) * (-dir_angle[1] * dir_angle[2] / cos_denom)) /
                   sqrt(SQUARE(R) - SQUARE(dist));

      dH_U_dX[2] = -1. + ((fv->x[0] - axis_pt[0]) * (-dir_angle[2] * dir_angle[0] / cos_denom) +
                          (fv->x[1] - axis_pt[1]) * (-dir_angle[2] * dir_angle[1] / cos_denom) +
                          (fv->x[2] - axis_pt[2]) * (-SQUARE(dir_angle[2]) / cos_denom)) /
                             sqrt(SQUARE(R) - SQUARE(dist));
    }
  }

  else if (mp->HeightUFunctionModel == CAP_SQUEEZE) {
    H_dot = mp->u_heightU_function_constants[0];
    H_low = mp->u_heightU_function_constants[1];
    R = mp->u_heightU_function_constants[2];
    dbl x_0 = mp->u_heightU_function_constants[3];
    dbl z_0 = mp->u_heightU_function_constants[4];

    *H_U = H_dot * time + H_low + R -
           sqrt(R * R - (fv->x[0] - x_0) * (fv->x[0] - x_0) - (fv->x[2] - z_0) * (fv->x[2] - z_0));
    *dH_U_dtime = H_dot;
    dH_U_dX[0] = (fv->x[0] - x_0) / sqrt(R * R - (fv->x[0] - x_0) * (fv->x[0] - x_0) -
                                         (fv->x[2] - z_0) * (fv->x[2] - z_0));
    dH_U_dX[1] = 0.;

    dH_U_dX[2] = (fv->x[2] - z_0) / sqrt(R * R - (fv->x[0] - x_0) * (fv->x[0] - x_0) -
                                         (fv->x[2] - z_0) * (fv->x[2] - z_0));
    *dH_U_ddh = 0.0;
  }

  else if ((mp->HeightUFunctionModel == FLAT_GRAD_FLAT) ||
           (mp->HeightUFunctionModel == FLAT_GRAD_FLAT_MELT)) {

    // Read in parameters
    dbl x1 = mp->u_heightU_function_constants[0];
    dbl h1 = mp->u_heightU_function_constants[1];
    dbl x2 = mp->u_heightU_function_constants[2];
    dbl h2 = mp->u_heightU_function_constants[3];
    dbl p = mp->u_heightU_function_constants[4];
    dbl x = fv->x[0];

    // Shortcuts
    dbl pp = pow(p, 2);
    dbl xx = pow(x, 2);
    dbl xx1 = pow(x1, 2);
    dbl xx2 = pow(x2, 2);

    // Define factors
    dbl n = (h2 - h1) / (x2 - x1);
    dbl f = h1 + n * (x - x1);
    dbl z1 = 0.5 + atan(p * (x - x1)) / PI;
    dbl z2 = 0.5 + atan(p * (x - x2)) / PI;
    dbl z1_x = p / (PI * (pp * xx - 2 * pp * x1 * x + pp * xx1 + 1));
    dbl z2_x = p / (PI * (pp * xx - 2 * pp * x2 * x + pp * xx2 + 1));

    // Assemble
    *H_U = (1 - z1) * h1 + z1 * (1 - z2) * f + z2 * h2;
    dH_U_dX[0] = -h1 * z1_x + z1_x * (1 - z2) * f - z1 * z2_x * f + z1 * (1 - z2) * n + z2_x * h2;
    dH_U_dX[1] = 0.0;
    dH_U_dX[2] = 0.0;
    *dH_U_dtime = 0.0;
    *dH_U_dp = 0.0;
    *dH_U_ddh = 0.0;

    // Add in any melting
    if (mp->HeightUFunctionModel == FLAT_GRAD_FLAT_MELT) {
      *H_U += fv->sh_dh;
    }

  }

  else if (mp->HeightUFunctionModel == POLY_TIME) {

    // Define variables and initialize
    int i;
    dbl np;
    *H_U = *dH_U_dtime = *dH_U_dp = *dH_U_ddh = 0.0;
    dH_U_dX[0] = dH_U_dX[1] = dH_U_dX[2] = 0.0;

    // Read in parameters
    np = mp->len_u_heightU_function_constants;

    // Assemble
    *H_U = mp->u_heightU_function_constants[0];
    for (i = 1; i < np; i++) {
      *H_U += mp->u_heightU_function_constants[i] * pow(time, i);
      *dH_U_dtime += mp->u_heightU_function_constants[i] * pow(time, i - 1) * i;
    }

  }

  else if (mp->HeightUFunctionModel == JOURNAL) {

    // Define variables and initialize
    dbl C, ecc;
    dbl x, y;
    dbl Ri, theta;
    *H_U = *dH_U_dtime = *dH_U_dp = *dH_U_ddh = 0.0;
    dH_U_dX[0] = dH_U_dX[1] = dH_U_dX[2] = 0.0;

    // Read in parameters
    C = mp->u_heightU_function_constants[0];
    ecc = mp->u_heightU_function_constants[1];

    // Read in current point
    x = fv->x[0];
    y = fv->x[1];

    // Calculate cylinder radius
    Ri = sqrt(x * x + y * y);

    // Calculate angle
    if (x > 0) {
      theta = acos(y / Ri);
    } else {
      theta = 2 * PI - acos(y / Ri);
    }

    // Calculate height and slopes
    *H_U = C * (1 + ecc * y / Ri);
    dH_U_dX[0] = -C * ecc / Ri * sin(theta) * cos(theta);
    dH_U_dX[1] = +C * ecc / Ri * sin(theta) * sin(theta);

  } else if (mp->HeightUFunctionModel == CIRCLE_MELT) {
    x_0 = mp->u_heightU_function_constants[0];
    r = mp->u_heightU_function_constants[1];
    H_low = mp->u_heightU_function_constants[2];

    Length = fv->x[0] - x_0;
    if (Length > 0.95 * r)
      GOMA_EH(GOMA_ERROR, "Problem in calculating height function model CIRCLE_MELT");

    *H_U = H_low + r - sqrt(r * r - Length * Length) + fv->sh_dh;
    *dH_U_dtime = Length / sqrt(r * r - Length * Length);
    *dH_U_dtime = 0.0;
    dH_U_dX[1] = 0.;
    dH_U_dX[2] = 0.;
    *dH_U_dp = 0.;
    *dH_U_ddh = 1.0;
  } else if (mp->HeightUFunctionModel == TABLE) {
    struct Data_Table *table_local;
    table_local = MP_Tables[mp->heightU_function_constants_tableid];

    if (!strcmp(table_local->t_name[0], "LINEAR_TIME")) {
      dbl time_local[1];
      time_local[0] = time;

      *H_U = interpolate_table(table_local, time_local, dH_U_dtime, NULL);
      dH_U_dX[0] = dH_U_dX[1] = dH_U_dX[2] = 0.0;
    }
  } else if (mp->HeightUFunctionModel == ROLLER) {
    // implement for bar elements in 2d space for now
    double hmin = mp->u_heightU_function_constants[0];
    double r = mp->u_heightU_function_constants[1];
    double xc = mp->u_heightU_function_constants[2];
    double external_field_multiplier = mp->u_heightU_function_constants[3];
    double x = fv->x[0];

    // we're all efv Sherman! It's likely that gap thickness
    // should be defined radially for this problem
    if (external_field_multiplier == 1.0) {
      *H_U = 0.0;
    } else {
      *H_U = hmin + r - sqrt(SQUARE(r) - SQUARE(x - xc));
    }

    /* add on external field height if there is one. The scale factor will be the fourth user
     * const*/
    if (mp->heightU_ext_field_index >= 0) {
      *H_U += mp->u_heightU_function_constants[3] * fv->external_field[mp->heightU_ext_field_index];
      if (*H_U < 0.0) {
        GOMA_WH(GOMA_ERROR, "read in a negative external field in height_function_model()");
      }
    }
    if (external_field_multiplier == 1.0) {
      dH_U_dX[0] = 0.0;
    } else {
      dH_U_dX[0] = (x - xc) / sqrt(SQUARE(r) - SQUARE(x - xc));
    }
    // dH_U_DX[0] = dH_ds for my_normal == primitive_s
    // so handle the external field gradients
    if (mp->heightU_ext_field_index >= 0) {

      // load ds_dcsi, that is det_J here
      double det_J;
      double d_det_J_dmeshkj[DIM][MDE];
      memset(d_det_J_dmeshkj, 0.0, sizeof(double) * DIM * MDE);
      detJ_2d_bar(&det_J, d_det_J_dmeshkj);

      int i;
      double dHext_ds, dHext_dcsi;
      dHext_ds = 0.0;
      dHext_dcsi = 0.0;

      // assume that the height field has the same dof as displacement
      for (i = 0; i < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; i++) {
        dHext_dcsi += mp->u_heightU_function_constants[3] *
                      *evp->external_field[mp->heightU_ext_field_index][i] *
                      bf[MESH_DISPLACEMENT1]->dphidxi[i][0];
      }

      dHext_ds = dHext_dcsi / det_J;
      dH_U_dX[0] += dHext_ds;
    } // end handling of the external field gradients
    dH_U_dX[1] = 0.0;
    dH_U_dX[2] = 0.0;

  } else {
    GOMA_EH(GOMA_ERROR, "Not a supported height-function model");
  }

  if (mp->HeightLFunctionModel == CONSTANT) {
    *H_L = mp->heightL;
    *dH_L_dtime = 0.0;
    dH_L_dX[0] = dH_L_dX[1] = dH_L_dX[2] = 0.0;
  }

  else if (mp->HeightLFunctionModel == CONSTANT_SPEED) {
    H_dot = mp->u_heightL_function_constants[0];
    H_init = mp->u_heightL_function_constants[1];

    *H_L = H_dot * time + H_init;

    /* add on external field height if there is one. The scale factor will be the third user const*/
    if (mp->heightL_ext_field_index >= 0) {
      *H_L += mp->u_heightL_function_constants[2] * fv->external_field[mp->heightL_ext_field_index];
    }

    *dH_L_dtime = H_dot;
    dH_L_dX[0] = dH_L_dX[1] = dH_L_dX[2] = 0.;
  } else if (mp->HeightLFunctionModel == EXTERNAL_FIELD) {
    H_dot = mp->u_heightL_function_constants[0];
    H_init = mp->u_heightL_function_constants[1];

    *H_L = H_dot * time + H_init;

    /* add on external field height if there is one. The scale factor will be the third user const*/
    *H_L += mp->u_heightL_function_constants[2] * fv->external_field[mp->heightL_ext_field_index];
    *dH_L_dtime = H_dot;
    dH_L_dX[0] = dH_L_dX[1] = dH_L_dX[2] = 0.;
  } else if (mp->HeightLFunctionModel == ROLL_ON) {
    Length = mp->u_heightL_function_constants[4];
    H_dot = mp->u_heightL_function_constants[3];
    H_delta = mp->u_heightL_function_constants[2];
    x_0 = mp->u_heightL_function_constants[0];
    H_low = mp->u_heightL_function_constants[1];

    *H_L = (H_dot * time + H_delta) * ((fv->x[0] - x_0) / Length) + H_low;
    *dH_L_dtime = H_dot * (fv->x[0] - x_0) / Length;
    dH_L_dX[0] = (H_dot * time + H_delta) / Length;
    dH_L_dX[1] = 0.;
    dH_U_dX[2] = 0.;
  } else if (mp->HeightLFunctionModel == ROLL) {
    R = mp->u_heightL_function_constants[0];
    /*  origin and direction of rotation axis	*/
    origin[0] = mp->u_heightL_function_constants[1];
    origin[1] = mp->u_heightL_function_constants[2];
    origin[2] = mp->u_heightL_function_constants[3];
    dir_angle[0] = mp->u_heightL_function_constants[4];
    dir_angle[1] = mp->u_heightL_function_constants[5];
    dir_angle[2] = mp->u_heightL_function_constants[6];
    H_dot = mp->u_heightL_function_constants[7];
    origin[2] += H_dot * time;

    /*  find intersection of axis with normal plane - i.e., locate point on
            axis that intersects plane normal to axis that contains local point. */

    cos_denom = (SQUARE(dir_angle[0]) + SQUARE(dir_angle[1]) + SQUARE(dir_angle[2]));
    t = (dir_angle[0] * (fv->x[0] - origin[0]) + dir_angle[1] * (fv->x[1] - origin[1]) +
         dir_angle[2] * (fv->x[2] - origin[2])) /
        cos_denom;
    axis_pt[0] = origin[0] + dir_angle[0] * t;
    axis_pt[1] = origin[1] + dir_angle[1] * t;
    axis_pt[2] = origin[2] + dir_angle[2] * t;

    /*  compute radial direction	*/

    dist = sqrt(SQUARE(fv->x[0] - axis_pt[0]) + SQUARE(fv->x[1] - axis_pt[1]));
    if (dist > fabs(R)) {
      *H_L = -fabs(R);
      dH_L_dX[0] = dH_L_dX[1] = dH_L_dX[2] = 0.0;
    } else {
      *H_L = SGN(R) * (axis_pt[2] - fv->x[2] - sqrt(SQUARE(R) - SQUARE(dist)));
      *dH_L_dtime = 0.; /* finish later  */
      dH_L_dX[0] = -((fv->x[0] - axis_pt[0]) * (1. - SQUARE(dir_angle[0]) / cos_denom) +
                     (fv->x[1] - axis_pt[1]) * (-dir_angle[0] * dir_angle[1] / cos_denom) +
                     (fv->x[2] - axis_pt[2]) * (-dir_angle[0] * dir_angle[2] / cos_denom)) /
                   sqrt(SQUARE(R) - SQUARE(dist));
      dH_L_dX[1] = -((fv->x[0] - axis_pt[0]) * (-dir_angle[1] * dir_angle[0] / cos_denom) +
                     (fv->x[1] - axis_pt[1]) * (1. - SQUARE(dir_angle[1]) / cos_denom) +
                     (fv->x[2] - axis_pt[2]) * (-dir_angle[1] * dir_angle[2] / cos_denom)) /
                   sqrt(SQUARE(R) - SQUARE(dist));
      dH_L_dX[2] = -1. - ((fv->x[0] - axis_pt[0]) * (-dir_angle[2] * dir_angle[0] / cos_denom) +
                          (fv->x[1] - axis_pt[1]) * (-dir_angle[2] * dir_angle[1] / cos_denom) +
                          (fv->x[2] - axis_pt[2]) * (-SQUARE(dir_angle[2]) / cos_denom)) /
                             sqrt(SQUARE(R) - SQUARE(dist));
    }
  } else if (mp->HeightLFunctionModel == TABLE) {
    struct Data_Table *table_local;
    table_local = MP_Tables[mp->heightL_function_constants_tableid];
    if (!strcmp(table_local->t_name[0], "LOWER_DISTANCE")) {
      /*
       * The LOWER_DISTANCE model looks at the lower velocity function to determine
       * how far the surface has traveled throughout the run.  It then
       * translates that into a position shift in the x direction.  This actual
       * function is a table look-up that has a lower height function model as a
       * function of distance.  However, that distance is shifted by the motion
       * of the lower surface.
       *
       * AUTHOR:    Scott A. Roberts, 1514
       * DATE:      May 2, 2012
       * LOCATION:  Somewhere over the US
       */

      // Calculate amount that the lower surface has shifted.
      int i;
      double np = 0, time_scale = 0, tn, disp = 0.0;
      if (mp->VeloLFunctionModel == SLIDER_POLY_TIME) {
        np = mp->len_u_veloL_function_constants;
        time_scale = mp->u_veloL_function_constants[0];
        tn = time_scale * time;
        for (i = 1; i < np; i++) {
          disp += mp->u_veloL_function_constants[i] / i * pow(tn, i);
        }
      } else {
        GOMA_EH(GOMA_ERROR, "To use LOWER_DISTANCE for Lower Height Function Model, "
                            "SLIDER_POLY_TIME must be used for the Lower Velocity Function.");
      }

      // Get height from table lookup
      double var1[1], slope;
      var1[0] = fv->x[0] - disp;
      *H_L = interpolate_table(table_local, var1, &slope, NULL);

      // Calculate spatial derivatives
      dH_L_dX[0] = -slope;
      dH_L_dX[1] = 0.0;
      dH_L_dX[2] = 0.0;

      // Calculate time derivative
      double H2;
      tn = time_scale * time * 1.0001;
      for (i = 1; i < np; i++) {
        disp += mp->u_veloL_function_constants[i] / i * pow(tn, i);
      }
      var1[0] = fv->x[0] - disp;
      H2 = interpolate_table(table_local, var1, &slope, NULL);
      *dH_L_dtime = (H2 - *H_L) / (time_scale * time * 0.0001);

    } else {
      GOMA_EH(GOMA_ERROR, "Lower Height Function does not know how to handle this type of TABLE.");
    }
  } else {
    GOMA_EH(GOMA_ERROR, "Not a supported height-function model");
  }

  H = MAX(*H_U - *H_L, DBL_SEMI_SMALL); /* Negative H would be seem to be a bad thing...*/
  return (H);
}

/*****************************************************************************/
double velocity_function_model(double veloU[DIM],
                               double veloL[DIM],
                               double time,    /* present time value           */
                               double delta_t) /* present time step             */

/******************************************************************************
 *
 *  A function which computes the height at the current time and the rate-of-change of
 *  height.  This model is used for the lubrication capability
 *
 *  P. Randall Schunk (March 2009, Somewhere over Texas)
 *
 *
 ******************************************************************************/

{
  double speed;
  double H_dot, hgt;
  double R, origin[3], dir_angle[3], t, axis_pt[3], rad_dir[3], dist;
  double omega, cos_denom, v_dir[3];

  if (mp->VeloUFunctionModel == CONSTANT) {
    veloU[0] = mp->veloU[0];
    veloU[1] = mp->veloU[1];
    veloU[2] = mp->veloU[2];
  } else if (mp->VeloUFunctionModel == LINEAR_TIME) {
    veloU[0] = mp->u_veloU_function_constants[0] + mp->u_veloU_function_constants[3] * time;
    veloU[1] = mp->u_veloU_function_constants[1] + mp->u_veloU_function_constants[4] * time;
    veloU[2] = 0.;
  } else if (mp->VeloUFunctionModel == ROLL) {
    if (mp->HeightUFunctionModel == ROLL) {
      R = mp->u_heightU_function_constants[0];
      /*  origin and direction of rotation axis	*/
      origin[0] = mp->u_heightU_function_constants[1];
      origin[1] = mp->u_heightU_function_constants[2];
      origin[2] = mp->u_heightU_function_constants[3];
      dir_angle[0] = mp->u_heightU_function_constants[4];
      dir_angle[1] = mp->u_heightU_function_constants[5];
      dir_angle[2] = mp->u_heightU_function_constants[6];
      H_dot = mp->u_heightU_function_constants[7];
      origin[2] += H_dot * time;

      /*  find intersection of axis with normal plane - i.e., locate point on
              axis that intersects plane normal to axis that contains local point. */

      cos_denom = (SQUARE(dir_angle[0]) + SQUARE(dir_angle[1]) + SQUARE(dir_angle[2]));
      t = (dir_angle[0] * (fv->x[0] - origin[0]) + dir_angle[1] * (fv->x[1] - origin[1]) +
           dir_angle[2] * (-origin[2])) /
          cos_denom;
      axis_pt[0] = origin[0] + dir_angle[0] * t;
      axis_pt[1] = origin[1] + dir_angle[1] * t;
      axis_pt[2] = origin[2] + dir_angle[2] * t;

      /*  compute radial direction	*/

      dist = sqrt(SQUARE(fv->x[0] - axis_pt[0]) + SQUARE(fv->x[1] - axis_pt[1]));
      if (dist > R) {
        veloU[0] = 0.0;
        veloU[1] = 0.0;
        veloU[2] = 0.0;
      } else {
        hgt = axis_pt[2] - sqrt(SQUARE(R) - SQUARE(dist));
        t = (dir_angle[0] * (fv->x[0] - origin[0]) + dir_angle[1] * (fv->x[1] - origin[1]) +
             dir_angle[2] * (hgt - origin[2])) /
            cos_denom;
        axis_pt[0] = origin[0] + dir_angle[0] * t;
        axis_pt[1] = origin[1] + dir_angle[1] * t;
        axis_pt[2] = origin[2] + dir_angle[2] * t;

        /*  compute radius and radial direction	*/

        rad_dir[0] = (fv->x[0] - axis_pt[0]) / R;
        rad_dir[1] = (fv->x[1] - axis_pt[1]) / R;
        rad_dir[2] = (hgt - axis_pt[2]) / R;

        /* compute velocity direction as perpendicular to both axis and radial
                direction.  Positive direction is determined by right hand rule */

        v_dir[0] = dir_angle[1] * rad_dir[2] - dir_angle[2] * rad_dir[1];
        v_dir[1] = dir_angle[2] * rad_dir[0] - dir_angle[0] * rad_dir[2];
        v_dir[2] = dir_angle[0] * rad_dir[1] - dir_angle[1] * rad_dir[0];

        omega = mp->u_veloU_function_constants[0];
        veloU[0] = omega * R * v_dir[0];
        veloU[1] = omega * R * v_dir[1];
        veloU[2] = omega * R * v_dir[2];
      }
    } else {
      GOMA_WH(GOMA_ERROR, "VelocityU and HeightU ROLL functions don't match.");
    }
  } else if (mp->VeloUFunctionModel == TANGENTIAL_ROTATE) {
    dbl n[DIM], v[DIM], t1[DIM], t2[DIM];
    dbl vel1, vel2, t1mag, t2mag;
    int i, j;
    /* dbl H; */
    dbl H_U, H_L, dH_U_dtime, dH_L_dtime, dH_U_dp, dH_U_ddh;
    dbl dH_U_dX[DIM], dH_L_dX[DIM];
    dbl thetax, thetay;
    dbl n2[DIM] = {0.0};
    dbl R[DIM][DIM];

    /* Import material parameters */
    for (i = 0; i < DIM; i++)
      v[i] = mp->u_veloU_function_constants[i];
    for (i = 0; i < DIM; i++)
      n[i] = fv->snormal[i];
    vel1 = mp->u_veloU_function_constants[3];
    vel2 = mp->u_veloU_function_constants[4];

    /* Rotate normal vector according to height function model slope */
    /* NOTE:  This functionality is not quite complete yet.  I'm
     * only confidient that it works when the shell is oriented with
     * the normal pointing in the z or -z direction, i.e. the shell
     * is the x-y plane.  Will have to see if it works in other
     * situations.  Regardless, if dH_U_dX is zero, this function will work
     * in any situation.  If you need dH_U_dX on a curved surface
     * we really need to think about what dH_U_dX really means. --SAR */

    /* H = height_function_model(&H_U, &dH_U_dtime, &H_L, &dH_L_dtime, dH_U_dX, dH_L_dX, &dH_U_dp,
     * &dH_U_ddh, time, delta_t); */

    height_function_model(&H_U, &dH_U_dtime, &H_L, &dH_L_dtime, dH_U_dX, dH_L_dX, &dH_U_dp,
                          &dH_U_ddh, time, delta_t);
    thetax = atan(dH_U_dX[1]);
    thetay = -atan(dH_U_dX[0]);
    R[0][0] = cos(thetay);
    R[0][1] = 0.0;
    R[0][2] = sin(thetay);
    R[1][0] = sin(thetax) * sin(thetay);
    R[1][1] = cos(thetax);
    R[1][2] = -sin(thetax) * cos(thetay);
    R[2][0] = -cos(thetax) * sin(thetay);
    R[2][1] = sin(thetax);
    R[2][2] = cos(thetax) * cos(thetay);
    for (i = 0; i < DIM; i++) {
      for (j = 0; j < DIM; j++) {
        n2[i] += R[i][j] * n[j];
      }
    }
    for (i = 0; i < DIM; i++)
      n[i] = n2[i];

    /* Calculate first tangent vector */
    t1[0] = n[1] * v[2] - n[2] * v[1];
    t1[1] = n[2] * v[0] - n[0] * v[2];
    t1[2] = n[0] * v[1] - n[1] * v[0];
    t1mag = sqrt(pow(t1[0], 2) + pow(t1[1], 2) + pow(t1[2], 2));
    for (i = 0; i < DIM; i++)
      t1[i] = t1[i] / t1mag;

    /* Calculate second tangent vector */
    t2[0] = n[1] * t1[2] - n[2] * t1[1];
    t2[1] = n[2] * t1[0] - n[0] * t1[2];
    t2[2] = n[0] * t1[1] - n[1] * t1[0];
    t2mag = sqrt(pow(t2[0], 2) + pow(t2[1], 2) + pow(t2[2], 2));
    for (i = 0; i < DIM; i++)
      t2[i] = t2[i] / t2mag;

    /* Calculate velocity components */
    for (i = 0; i < DIM; i++)
      veloU[i] = t1[i] * vel1 + t2[i] * vel2;
  } else {
    GOMA_EH(GOMA_ERROR, "Not a supported velocity-function model");
  }

  if (mp->VeloLFunctionModel == CONSTANT) {
    veloL[0] = mp->veloL[0];
    veloL[1] = mp->veloL[1];
    veloL[2] = mp->veloL[2];
  } else if (mp->VeloLFunctionModel == LINEAR_TIME) {
    veloL[0] = mp->u_veloL_function_constants[0] + mp->u_veloL_function_constants[3] * time;
    veloL[1] = mp->u_veloL_function_constants[1] + mp->u_veloL_function_constants[4] * time;
    veloL[2] = 0.;
  } else if (mp->VeloLFunctionModel == USER) {
    GOMA_EH(GOMA_ERROR,
            "USER is no longer a supported lower velocity-function model (see SLIDER_POLY_TIME).");
  } else if (mp->VeloLFunctionModel == SLIDER_POLY_TIME) {
    // Define variables and initialize
    int i;
    dbl np, time_scale, tn;

    // Read in parameters
    np = mp->len_u_veloL_function_constants;
    time_scale = mp->u_veloL_function_constants[0];
    tn = time_scale * time;

    // Assemble Note: This model assumes X-direction ONLY!!!!!!
    veloL[0] = mp->u_veloL_function_constants[1];
    for (i = 2; i < np; i++) {
      veloL[0] += mp->u_veloL_function_constants[i] * pow(tn, i - 1);
    }

    veloL[1] = 0.;
    veloL[2] = 0.;
  } else if (mp->VeloLFunctionModel == ROLL) {
    if (mp->HeightLFunctionModel == ROLL) {
      R = mp->u_heightL_function_constants[0];
      /*  origin and direction of rotation axis	*/
      origin[0] = mp->u_heightL_function_constants[1];
      origin[1] = mp->u_heightL_function_constants[2];
      origin[2] = mp->u_heightL_function_constants[3];
      dir_angle[0] = mp->u_heightL_function_constants[4];
      dir_angle[1] = mp->u_heightL_function_constants[5];
      dir_angle[2] = mp->u_heightL_function_constants[6];
      H_dot = mp->u_heightL_function_constants[7];
      origin[2] += H_dot * time;

      /*  find intersection of axis with normal plane - i.e., locate point on
              axis that intersects plane normal to axis that contains local point. */

      cos_denom = (SQUARE(dir_angle[0]) + SQUARE(dir_angle[1]) + SQUARE(dir_angle[2]));
      t = (dir_angle[0] * (fv->x[0] - origin[0]) + dir_angle[1] * (fv->x[1] - origin[1]) +
           dir_angle[2] * (-origin[2])) /
          cos_denom;
      axis_pt[0] = origin[0] + dir_angle[0] * t;
      axis_pt[1] = origin[1] + dir_angle[1] * t;
      axis_pt[2] = origin[2] + dir_angle[2] * t;

      /*  compute radial direction	*/

      dist = sqrt(SQUARE(fv->x[0] - axis_pt[0]) + SQUARE(fv->x[1] - axis_pt[1]));
      if (dist > R) {
        veloL[0] = 0.0;
        veloL[1] = 0.0;
        veloL[2] = 0.0;
      } else {
        hgt = axis_pt[2] + sqrt(SQUARE(R) - SQUARE(dist));
        t = (dir_angle[0] * (fv->x[0] - origin[0]) + dir_angle[1] * (fv->x[1] - origin[1]) +
             dir_angle[2] * (hgt - origin[2])) /
            cos_denom;
        axis_pt[0] = origin[0] + dir_angle[0] * t;
        axis_pt[1] = origin[1] + dir_angle[1] * t;
        axis_pt[2] = origin[2] + dir_angle[2] * t;

        /*  compute radius and radial direction	*/

        rad_dir[0] = (fv->x[0] - axis_pt[0]) / R;
        rad_dir[1] = (fv->x[1] - axis_pt[1]) / R;
        rad_dir[2] = (hgt - axis_pt[2]) / R;

        /* compute velocity direction as perpendicular to both axis and radial
                direction.  Positive direction is determined by right hand rule */

        v_dir[0] = dir_angle[1] * rad_dir[2] - dir_angle[2] * rad_dir[1];
        v_dir[1] = dir_angle[2] * rad_dir[0] - dir_angle[0] * rad_dir[2];
        v_dir[2] = dir_angle[0] * rad_dir[1] - dir_angle[1] * rad_dir[0];

        omega = mp->u_veloL_function_constants[0];
        veloL[0] = omega * R * v_dir[0];
        veloL[1] = omega * R * v_dir[1];
        veloL[2] = omega * R * v_dir[2];
      }
    } else {
      GOMA_WH(GOMA_ERROR, "VelocityL and HeightL ROLL functions don't match.");
    }
  } else if (mp->VeloLFunctionModel == TANGENTIAL_ROTATE) {
    dbl n[DIM], v[DIM], t1[DIM], t2[DIM];
    dbl vel1, vel2, t1mag, t2mag;
    int i, j;
    /* dbl H; */
    dbl H_U, H_L, dH_U_dtime, dH_L_dtime, dH_U_dp, dH_U_ddh;
    dbl dH_U_dX[DIM], dH_L_dX[DIM];
    dbl thetax, thetay;
    dbl n2[DIM] = {0.0};
    dbl R[DIM][DIM];

    /* Import material parameters */
    for (i = 0; i < DIM; i++)
      v[i] = mp->u_veloL_function_constants[i];
    for (i = 0; i < DIM; i++)
      n[i] = fv->snormal[i];
    vel1 = mp->u_veloL_function_constants[3];
    vel2 = mp->u_veloL_function_constants[4];

    /* Rotate normal vector according to height function model slope */
    /* NOTE:  This functionality is not quite complete yet.  I'm
     * only confidient that it works when the shell is oriented with
     * the normal pointing in the z or -z direction, i.e. the shell
     * is the x-y plane.  Will have to see if it works in other
     * situations.  Regardless, if dH_L_dX is zero, this function will work
     * in any situation.  If you need dH_L_dX on a curved surface
     * we really need to think about what dH_L_dX really means. --SAR */

    /* H = height_function_model(&H_U, &dH_U_dtime, &H_L, &dH_L_dtime, dH_U_dX, dH_L_dX, &dH_U_dp,
     * &dH_U_ddh, time, delta_t); */

    height_function_model(&H_U, &dH_U_dtime, &H_L, &dH_L_dtime, dH_U_dX, dH_L_dX, &dH_U_dp,
                          &dH_U_ddh, time, delta_t);
    thetax = atan(dH_L_dX[1]);
    thetay = -atan(dH_L_dX[0]);
    R[0][0] = cos(thetay);
    R[0][1] = 0.0;
    R[0][2] = sin(thetay);
    R[1][0] = sin(thetax) * sin(thetay);
    R[1][1] = cos(thetax);
    R[1][2] = -sin(thetax) * cos(thetay);
    R[2][0] = -cos(thetax) * sin(thetay);
    R[2][1] = sin(thetax);
    R[2][2] = cos(thetax) * cos(thetay);
    for (i = 0; i < DIM; i++) {
      for (j = 0; j < DIM; j++) {
        n2[i] += R[i][j] * n[j];
      }
    }
    for (i = 0; i < DIM; i++)
      n[i] = n2[i];

    /* Calculate first tangent vector */
    t1[0] = n[1] * v[2] - n[2] * v[1];
    t1[1] = n[2] * v[0] - n[0] * v[2];
    t1[2] = n[0] * v[1] - n[1] * v[0];
    t1mag = sqrt(pow(t1[0], 2) + pow(t1[1], 2) + pow(t1[2], 2));
    for (i = 0; i < DIM; i++)
      t1[i] = t1[i] / t1mag;

    /* Calculate second tangent vector */
    t2[0] = n[1] * t1[2] - n[2] * t1[1];
    t2[1] = n[2] * t1[0] - n[0] * t1[2];
    t2[2] = n[0] * t1[1] - n[1] * t1[0];
    t2mag = sqrt(pow(t2[0], 2) + pow(t2[1], 2) + pow(t2[2], 2));
    for (i = 0; i < DIM; i++)
      t2[i] = t2[i] / t2mag;

    /* Calculate velocity components */
    for (i = 0; i < DIM; i++)
      veloL[i] = t1[i] * vel1 + t2[i] * vel2;
  } else {
    GOMA_EH(GOMA_ERROR, "Not a supported velocity-function model");
  }

  /* Note that this is the relative velocity of the upper and lower surfaces */
  speed = sqrt(pow((veloU[0] - veloL[0]), 2) + pow((veloU[1] - veloL[1]), 2) +
               pow((veloU[2] - veloL[2]), 2));

  return (speed);
}

/*****************************************************************************/
double
film_evaporation_model(double C,             /* Suspension concentration */
                       double *dEvapRate_dC, /* Sensitivity of evap rate w.r.t. concentration */
                       double H,             /* Film thickness */
                       double *dEvapRate_dH) /* Sensitivity of evap rate w.r.t. film thickness */
/******************************************************************************
 *
 *  A function which computes the rate of evaporation of solvent inside thin film
 *  This model is used for the lubrication capability
 *
 *  Kris Tjiptowidjojo (January 25 2010)
 *
 *
 ******************************************************************************/

{
  double EvapRate = 0.0;
  double nexp;
  double E_0;
  double CMax;

  double H_star = 0.0;
  double F;
  double width;
  double alpha;
  double step_func = 1.0;
  double dstep_func_dH = 0.0;

  if (mp->FilmEvapModel == CONSTANT) {
    EvapRate = mp->FilmEvap;
    *dEvapRate_dC = 0.;
  }

  else if (mp->FilmEvapModel == CONC_POWER) {
    E_0 = mp->u_FilmEvap_function_constants[0];
    nexp = mp->u_FilmEvap_function_constants[1];
    CMax = mp->u_FilmEvap_function_constants[2];

    /*   Implement a cutoff concentration that caps the evaporation
         Here I have chosen a cutoff value 1% below maximum packing. However,
         this may be too low and I will have to adjust it.
    */

    if (C > 0.0 && C < 0.95 * CMax) {
      EvapRate = E_0 * pow(1.0 - C / CMax, nexp);

      *dEvapRate_dC = -nexp * (E_0 / CMax) * pow(1.0 - C / CMax, nexp - 1.);
    }

    else if (C <= 0.) {
      EvapRate = 0.0;
      *dEvapRate_dC = 0.0;
    }

    else if (C >= 0.95 * CMax) {

      if (C < 0.999 * CMax) {
        EvapRate = E_0 * pow(1.0 - C / CMax, nexp);
      } else {
        EvapRate = 0.0;
      }

      /*       *dEvapRate_dC = -nexp * (E_0/CMax) * pow(1.0 - 0.95, nexp - 1.); */
      *dEvapRate_dC = 0.0;
    }

  } else {
    GOMA_EH(GOMA_ERROR, "Not a supported film evaporation model");
  }

  /* Pre-multiply the evaporation function with heaviside function
     so that the evaporation dies off at H = H_star */

  if (mp->DisjPressModel == TWO_TERM) {
    H_star = mp->u_DisjPress_function_constants[3];
  } else if ((mp->DisjPressModel == TWO_TERM_EXT_CA) || (mp->DisjPressModel == ONE_TERM)) {
    H_star = mp->u_DisjPress_function_constants[2];
  }

  width = 0.2 * H_star;
  alpha = 0.5 * width;

  F = H - (H_star + alpha);

  if (alpha != 0.) {

    if (F < (-1. * alpha)) {
      step_func = 0.0;
      dstep_func_dH = 0.0;
    } else if (F > alpha) {
      step_func = 1.0;
      dstep_func_dH = 0.0;
    } else {
      step_func = 0.5 * (1. + F / alpha + sin(M_PIE * F / alpha) / M_PIE);
      dstep_func_dH = 0.5 * (1. / alpha + cos(M_PIE * F / alpha) / alpha);
    }
  }

  *dEvapRate_dH = EvapRate * dstep_func_dH;
  *dEvapRate_dC *= step_func;
  EvapRate *= step_func;

  return (EvapRate);
}

/*****************************************************************************/
double disjoining_pressure_model(
    double H,                             /* Film thickness or interfacial separation */
    double grad_H[DIM],                   /* Film slope */
    int *n_dof,                           /* Numbers of degrees of freedom array */
    int *dof_map,                         /* Map of DOFs */
    double grad_DisjPress[DIM],           /* Disjoining pressure gradient */
    double dgrad_DisjPress_dH1[DIM][MDE], /* Sensitivity of disjoining pressure w.r.t. film thicknes
                                             - grad_phi_j terms */
    double dgrad_DisjPress_dH2[DIM][MDE], /* Sensitivity of disjoining pressure w.r.t. film thicknes
                                             - phi_j terms */
    double dgrad_DisjPress_dH[DIM]
                             [MDE]) /* Sensitivity of disjoining pressure w.r.t. film thickness */
/******************************************************************************
 *
 *  A function which computes disjoining pressure inside thin film
 *  This model is used for the lubrication capability
 *
 *  Kris Tjiptowidjojo (January 26 2010)
 *
 *
 ******************************************************************************/

{
  int i, j;
  double DisjPress = 0.0;
  double grad_II_H[DIM];
  double angle;
  double grad_angle[DIM];
  double B = 0;
  double dB_dangle = 0;
  double f = 0;
  double df_dH = 0, d2f_dH2 = 0;
  double nexp;
  double mexp;
  double H_star;
  double factor;

  double phi_j, grad_phi_j[DIM], grad_II_phi_j[DIM], d_grad_II_phi_j_dmesh[DIM][DIM][MDE];

  Inn(grad_H, grad_II_H);

  memset(grad_angle, 0.0, sizeof(double) * DIM);

  memset(grad_DisjPress, 0.0, sizeof(double) * DIM);
  memset(dgrad_DisjPress_dH1, 0.0, sizeof(double) * DIM * MDE);
  memset(dgrad_DisjPress_dH2, 0.0, sizeof(double) * DIM * MDE);
  memset(dgrad_DisjPress_dH, 0.0, sizeof(double) * DIM * MDE);

  /************** CALCULATE DISJOINING PRESSURE AND ITS SENSITIVITIES **************/

  if (mp->DisjPressModel == CONSTANT) {
    B = 0.0;
    f = 0.0;
    DisjPress = mp->DisjPress;
    dB_dangle = 0.0;
    df_dH = 0.0;
    d2f_dH2 = 0.0;
    return (DisjPress);
  }

  else if (mp->DisjPressModel == TWO_TERM) {
    angle = mp->u_DisjPress_function_constants[0];
    nexp = mp->u_DisjPress_function_constants[1];
    mexp = mp->u_DisjPress_function_constants[2];
    H_star = mp->u_DisjPress_function_constants[3];
    factor = mp->u_DisjPress_function_constants[4];

    B = (mp->surface_tension / H_star) * (nexp - 1.) * (mexp - 1.) *
        (1. - cos(angle * M_PIE / 180.)) / (factor * (nexp - 1.) - (mexp - 1.));

    f = pow(H_star / H, nexp) - factor * pow(H_star / H, mexp);

    DisjPress = B * f;

    dB_dangle = (mp->surface_tension / H_star) * (nexp - 1.) * (mexp - 1.) *
                sin(angle * M_PIE / 180.) * (M_PIE / 180.) / (factor * (nexp - 1.) - (mexp - 1.));

    df_dH = -nexp * pow(H_star, nexp) / pow(H, nexp + 1.) +
            factor * mexp * pow(H_star, mexp) / pow(H, mexp + 1.);

    d2f_dH2 = nexp * (nexp + 1.) * pow(H_star, nexp) / pow(H, nexp + 2.) -
              factor * mexp * (mexp + 1.) * pow(H_star, mexp) / pow(H, mexp + 2.);

  }

  else if (mp->DisjPressModel == TWO_TERM_EXT_CA) {
    if (!(efv->ev)) {
      GOMA_EH(GOMA_ERROR,
              "Model TWO_TERM_EXT_CA requires contact angle input from external file !");
    }

    angle = fv->external_field[0];

    Inn(fv->grad_ext_field[0], grad_angle);

    nexp = mp->u_DisjPress_function_constants[0];
    mexp = mp->u_DisjPress_function_constants[1];
    H_star = mp->u_DisjPress_function_constants[2];
    factor = mp->u_DisjPress_function_constants[3];

    B = (mp->surface_tension / H_star) * (nexp - 1.) * (mexp - 1.) *
        (1. - cos(angle * M_PIE / 180.)) / (factor * (nexp - 1.) - (mexp - 1.));

    f = pow(H_star / H, nexp) - factor * pow(H_star / H, mexp);

    DisjPress = B * f;

    dB_dangle = (mp->surface_tension / H_star) * (nexp - 1.) * (mexp - 1.) *
                sin(angle * M_PIE / 180.) * (M_PIE / 180.) / (factor * (nexp - 1.) - (mexp - 1.));

    df_dH = -nexp * pow(H_star, nexp) / pow(H, nexp + 1.) +
            factor * mexp * pow(H_star, mexp) / pow(H, mexp + 1.);

    d2f_dH2 = nexp * (nexp + 1.) * pow(H_star, nexp) / pow(H, nexp + 2.) -
              factor * mexp * (mexp + 1.) * pow(H_star, mexp) / pow(H, mexp + 2.);

  }

  else if (mp->DisjPressModel == ONE_TERM) {
    B = mp->u_DisjPress_function_constants[0];
    nexp = mp->u_DisjPress_function_constants[1];
    H_star = mp->u_DisjPress_function_constants[2];

    f = pow(H_star / H, nexp);

    DisjPress = B * f;

    dB_dangle = 0.0;

    df_dH = -nexp * pow(H_star, nexp) / pow(H, nexp + 1.);

    d2f_dH2 = nexp * (nexp + 1.) * pow(H_star, nexp) / pow(H, nexp + 2.);
  }

  else {
    GOMA_EH(GOMA_ERROR, "Not a supported disjoining pressure model");
  }

  /************** CALCULATE DISJOINING PRESSURE GRADIENT AND ITS SENSITIVITIES **************/

  for (i = 0; i < DIM; i++) {
    grad_DisjPress[i] = dB_dangle * f * grad_angle[i] + B * df_dH * grad_II_H[i];
  }

  for (i = 0; i < DIM; i++) {
    for (j = 0; j < ei[pg->imtrx]->dof[SHELL_FILMH]; j++) {
      dgrad_DisjPress_dH1[i][j] = B * df_dH;
    }
  }

  for (i = 0; i < DIM; i++) {
    for (j = 0; j < ei[pg->imtrx]->dof[SHELL_FILMH]; j++) {
      dgrad_DisjPress_dH2[i][j] = dB_dangle * df_dH * grad_angle[i] + B * d2f_dH2 * grad_II_H[i];
    }
  }

  for (i = 0; i < DIM; i++) {
    grad_DisjPress[i] = dB_dangle * f * grad_angle[i] + B * df_dH * grad_II_H[i];

    for (j = 0; j < ei[pg->imtrx]->dof[SHELL_FILMH]; j++) {
      ShellBF(SHELL_FILMH, j, &phi_j, grad_phi_j, grad_II_phi_j, d_grad_II_phi_j_dmesh,
              n_dof[MESH_DISPLACEMENT1], dof_map);

      dgrad_DisjPress_dH[i][j] = dB_dangle * df_dH * phi_j * grad_angle[i] +
                                 B * d2f_dH2 * phi_j * grad_II_H[i] + B * df_dH * grad_II_phi_j[i];
    }
  }

  return (DisjPress);
}

/******************************************************************************/

double diffusion_coefficient_model(double mu,              /* Viscosity of the liquid medium */
                                   double *dDiffCoeff_dmu) /* diffusion coefficient's sensitivity
                                                               w.r.t. viscosity */

/******************************************************************************
 *
 *  A function which computes particles diffusion of coefficient inside thin film
 *  This model is used for the lubrication capability
 *
 *  Kris Tjiptowidjojo (March 9 2010)
 *
 *
 ******************************************************************************/

{
  double DiffCoeff = 0.0;
  double Boltz;
  double Temp;
  double R_part;

  if (mp->DiffCoeffModel == CONSTANT) {
    DiffCoeff = mp->DiffCoeff;
    *dDiffCoeff_dmu = 0.;
  }

  else if (mp->DiffCoeffModel == STOKES_EINSTEIN) {
    Boltz = mp->u_DiffCoeff_function_constants[0];
    Temp = mp->u_DiffCoeff_function_constants[1];
    R_part = mp->u_DiffCoeff_function_constants[2];

    DiffCoeff = (Boltz * Temp) / (6. * M_PIE * R_part * mu);
    *dDiffCoeff_dmu = -DiffCoeff / mu;

  }

  else {
    GOMA_EH(GOMA_ERROR, "Not a supported diffusion coefficient model");
  }

  return (DiffCoeff);
}
/*****************************************************************************/
double porous_shell_closed_porosity_model() {
  /******************************************************************************
   *
   *  This function computes the porosity of a structured porous shell based on
   *  either a constant value or reading data from a file.  This model is
   *  used with the porous shell capability in assemble_porous_shell.
   *
   *  Scott A. Roberts (sarober@sandia.gov) - March 2010
   *
   ******************************************************************************/
  dbl phi = 0.0;

  if (mp->PorousShellClosedPorosityModel == CONSTANT) {
    phi = mp->PorousShellClosedPorosity;

  } else if (mp->PorousShellClosedPorosityModel == EXTERNAL_FIELD) {
    phi = mp->u_PorousShellClosedPorosity_function_constants[0] *
          (fv->external_field[mp->por_shell_closed_porosity_ext_field_index] + 0.01);

  } else if (mp->PorousShellClosedPorosityModel == MULTI_MODE) {
    phi = mp->u_PorousShellClosedPorosity_function_constants[0];

  } else {
    GOMA_EH(GOMA_ERROR, "Not a supported porosity model");
  }
  return (phi);
}
/* END of porous_shell_porosity_model */

/*****************************************************************************/
double porous_shell_porosity_model(int k) {
  /******************************************************************************
   *
   *  This function computes the porosity of a structured porous shell based on
   *  either a constant value or reading data from a file.  This model is
   *  used with the porous shell capability in assemble_porous_shell.
   *
   *  Scott A. Roberts (sarober@sandia.gov) - March 2010
   *
   ******************************************************************************/
  dbl phi = 0.0;

  if (mp->PorousShellPorosityModel[k] == CONSTANT) {
    phi = mp->PorousShellPorosity[k];

  } else if (mp->PorousShellPorosityModel[k] == EXTERNAL_FIELD) {
    phi = mp->u_PorousShellPorosity[k][0] *
          (fv->external_field[mp->por_shell_porosity_ext_field_index[k]] + 0.01);
  } else {
    GOMA_EH(GOMA_ERROR, "Not a supported porosity model");
  }
  return (phi);
}
/* END of porous_shell_porosity_model */

/*****************************************************************************/
double porous_shell_closed_radius_model() {
  /******************************************************************************
   *
   *  This function computes the radius of a structured porous shell based on
   *  either a constant value or an external field.  Used with the function
   *  assemble_porous_shell.
   *
   *  Scott A. Roberts (sarober@sandia.gov) - March 2010
   *
   ******************************************************************************/
  dbl r = 0.0;

  if (mp->PorousShellClosedRadiusModel == CONSTANT) {
    r = mp->PorousShellClosedRadius;

  } else if (mp->PorousShellClosedRadiusModel == EXTERNAL_FIELD) {
    r = mp->u_PorousShellClosedRadius_function_constants[0] *
        fv->external_field[mp->por_shell_closed_radius_ext_field_index];

  } else if (mp->PorousShellClosedRadiusModel == MULTI_MODE) {
    r = 1.0;

  } else {
    GOMA_EH(GOMA_ERROR, "Not a supported radius model");
  }

  return (r);
}
/* END of porous_shell_radius_model */

/*****************************************************************************/
double porous_shell_closed_height_model() {
  /******************************************************************************
   *
   *  This function computes the height of a structured porous shell based on
   *  either a constant value or an external field.  Used with the function
   *  assemble_porous_shell.
   *
   *  Scott A. Roberts (sarober@sandia.gov) - March 2010
   *
   ******************************************************************************/
  dbl H = 0.0;
  if (mp->PorousShellClosedHeightModel == CONSTANT) {
    H = mp->PorousShellClosedHeight;
  } else if (mp->PorousShellClosedHeightModel == EXTERNAL_FIELD) {
    H = mp->u_PorousShellClosedHeight_function_constants[0] *
        fv->external_field[mp->por_shell_closed_height_ext_field_index];
  } else {
    GOMA_EH(GOMA_ERROR, "Not a supported height model");
  }
  return (H);
}
/* END of porous_shell_height_model */

/*****************************************************************************/

double porous_shell_height_model(int ipore) {
  /******************************************************************************
   *
   *  This function computes the height of a structured porous shell based on
   *  either a constant value or an external field.  Used with the function
   *  assemble_porous_shell.
   *
   *  Scott A. Roberts (sarober@sandia.gov) - March 2010
   *
   ******************************************************************************/
  dbl H = 0.0;
  if (mp->PorousShellHeightModel[ipore] == CONSTANT) {
    H = mp->PorousShellHeight[ipore];
  } else if (mp->PorousShellHeightModel[ipore] == EXTERNAL_FIELD) {
    H = mp->u_PorousShellHeight[ipore][0] *
        fv->external_field[mp->por_shell_height_ext_field_index[ipore]];
  } else {
    GOMA_EH(GOMA_ERROR, "Not a supported height model");
  }
  return (H);
}
/* END of porous_shell_height_model */

/*****************************************************************************/

double porous_shell_cross_perm_model(int ipore) {
  /******************************************************************************
   *
   *  This function computes the cross permeability  of an open porous shell based on
   *  either a constant value or an external field.  Used with the function
   *  assemble_porous_shell_open and assemble_porous_shell_open_2.
   *
   *  Kristianto Tjiptowidjojo (tjiptowi@unm.edu) - April 2017
   *
   *  Updated in December 4 2019 to be used with saturation formulation of porous shell
   *
   ******************************************************************************/
  dbl kappa = 0.0;

  if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN]) /* Pressure formulation */
  {
    if (mp->PorousShellCrossKappaModel == CONSTANT) {
      kappa = mp->PorousShellCrossKappa;
    } else if (mp->PorousShellCrossKappaModel == EXTERNAL_FIELD) {
      GOMA_EH(mp->Xperm_external_field_index, "Cross Permeability external field not found!");
      kappa = mp->PorousShellCrossKappa = mp->u_PorousShellCrossKappa_function_constants[0] *
                                          fv->external_field[mp->Xperm_external_field_index];
      if (pd->TimeIntegration == TRANSIENT) {
        mp_old->PorousShellCrossKappa = mp->u_PorousShellCrossKappa_function_constants[0] *
                                        fv->external_field[mp->Xperm_external_field_index];
      }
    } else {
      GOMA_EH(GOMA_ERROR, "Unrecognized Porous Shell Cross Permeability  model");
    }
  } else if (pd->e[pg->imtrx][R_SHELL_SAT_1]) /* Saturation formulation */
  {
    if (mp->PorousShellCrossPermeabilityModel[ipore] == CONSTANT) {
      kappa = mp->PorousShellCrossPermeability[ipore];
    } else if (mp->PorousShellCrossPermeabilityModel[ipore] == EXTERNAL_FIELD) {
      GOMA_EH(mp->por_shell_cross_permeability_ext_field_index[ipore],
              "Cross Permeability external field not found!");
      kappa = mp->PorousShellCrossPermeability[ipore] =
          mp->u_PorousShellCrossPermeability[ipore][0] *
          fv->external_field[mp->por_shell_cross_permeability_ext_field_index[ipore]];
    } else {
      GOMA_EH(GOMA_ERROR, "Unrecognized Porous Shell Cross Permeability  model");
    }
  }
  return (kappa);
}
/* END of porous_shell_cross_perm_model */

/*****************************************************************************/

double porous_shell_rel_perm_model(int ipore, double saturation) {
  /******************************************************************************
   *
   *  This function computes the relative permeability of an open porous shell.
   *  Used with the function assemble_porous_shell_saturation.
   *
   *  Kristianto Tjiptowidjojo (tjiptowi@unm.edu) - October 2018
   *
   ******************************************************************************/
  double k_rel = 0.0;
  double s_eff, d_s_eff, a1, factor, factor2;
  double expon2, sat_min, sat_max, viscosity, lambda;
  int i_rel_perm_ev;
  double scale;

  switch (mp->PorousShellRelPermModel[ipore]) {

  case CONSTANT:

    k_rel = mp->PorousShellRelPerm[ipore];
    mp->d_PorousShellRelPerm[ipore][SHELL_SAT_1] = 0.0;
    mp->d_PorousShellRelPerm[ipore][SHELL_SAT_2] = 0.0;
    mp->d_PorousShellRelPerm[ipore][SHELL_SAT_3] = 0.0;

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
      mp->d_PorousShellRelPerm[ipore][SHELL_SAT_1] = 0.0;
      mp->d_PorousShellRelPerm[ipore][SHELL_SAT_2] = 0.0;
      mp->d_PorousShellRelPerm[ipore][SHELL_SAT_3] = 0.0;
    }

    /*
     *  Clip the relative permeability at one -> it can never be
     *  greater than one.  Actually, note that if somehow s_eff is
     *  very close to 1.0 and fails this test, then you are dividing
     *  by zero as factor=1.0 below.    Now and then GOMA aborts due
     *  to this.
     */
    else if (s_eff >= 0.99999) {
      k_rel = mp->PorousShellRelPerm[ipore] = 1.0 / viscosity;
      mp->d_PorousShellRelPerm[ipore][SHELL_SAT_1] = 0.0;
      mp->d_PorousShellRelPerm[ipore][SHELL_SAT_2] = 0.0;
      mp->d_PorousShellRelPerm[ipore][SHELL_SAT_3] = 0.0;
    }

    /*
     *  Otherwise, apply Van Genuchten formula
     */
    else {
      expon2 = 1.0 / lambda;
      factor = pow(s_eff, expon2);
      factor2 = pow(1.0 - factor, lambda);
      a1 = 1.0 - factor2;
      k_rel = mp->PorousShellRelPerm[ipore] = sqrt(s_eff) * a1 * a1 / viscosity;

      /*
       *   Calculate Jacobian entries
       */
      if (a1 == 0.0) {
        mp->d_PorousShellRelPerm[ipore][SHELL_SAT_1] = 0.0;
        mp->d_PorousShellRelPerm[ipore][SHELL_SAT_2] = 0.0;
        mp->d_PorousShellRelPerm[ipore][SHELL_SAT_3] = 0.0;
      }
      d_s_eff = 1.0 / (sat_max - sat_min);
      switch (ipore) {
      case 0:
        mp->d_PorousShellRelPerm[ipore][SHELL_SAT_1] =
            d_s_eff * k_rel * (0.5 / s_eff + 2. * factor2 * factor / ((1.0 - factor) * s_eff * a1));
        break;

      case 1:
        mp->d_PorousShellRelPerm[ipore][SHELL_SAT_2] =
            d_s_eff * k_rel * (0.5 / s_eff + 2. * factor2 * factor / ((1.0 - factor) * s_eff * a1));
        break;

      case 2:
        mp->d_PorousShellRelPerm[ipore][SHELL_SAT_3] =
            d_s_eff * k_rel * (0.5 / s_eff + 2. * factor2 * factor / ((1.0 - factor) * s_eff * a1));
        break;
      }
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
    mp->d_PorousShellRelPerm[ipore][SHELL_SAT_1] = 0.0;
    mp->d_PorousShellRelPerm[ipore][SHELL_SAT_2] = 0.0;
    mp->d_PorousShellRelPerm[ipore][SHELL_SAT_3] = 0.0;

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
      mp->d_PorousShellRelPerm[ipore][SHELL_SAT_1] = 0.0;
      mp->d_PorousShellRelPerm[ipore][SHELL_SAT_2] = 0.0;
      mp->d_PorousShellRelPerm[ipore][SHELL_SAT_3] = 0.0;
    }

    /*
     *  Clip the relative permeability at one -> it can never be
     *  greater than one.  Actually, note that if somehow s_eff is
     *  very close to 1.0 and fails this test, then you are dividing
     *  by zero as factor=1.0 below.    Now and then GOMA aborts due
     *  to this.
     */
    else if (s_eff >= 0.99999) {
      k_rel = mp->PorousShellRelPerm[ipore] = 1.0 / viscosity;
      mp->d_PorousShellRelPerm[ipore][SHELL_SAT_1] = 0.0;
      mp->d_PorousShellRelPerm[ipore][SHELL_SAT_2] = 0.0;
      mp->d_PorousShellRelPerm[ipore][SHELL_SAT_3] = 0.0;
    }

    /*
     *  Otherwise, apply Van Genuchten formula
     */
    else {
      expon2 = 1.0 / lambda;
      factor = pow(s_eff, expon2);
      factor2 = pow(1.0 - factor, lambda);
      a1 = 1.0 - factor2;
      k_rel = mp->PorousShellRelPerm[ipore] = sqrt(s_eff) * a1 * a1 / viscosity;

      /*
       *   Calculate Jacobian entries
       */
      if (a1 == 0.0) {
        mp->d_PorousShellRelPerm[ipore][SHELL_SAT_1] = 0.0;
        mp->d_PorousShellRelPerm[ipore][SHELL_SAT_2] = 0.0;
        mp->d_PorousShellRelPerm[ipore][SHELL_SAT_3] = 0.0;
      }
      d_s_eff = 1.0 / (sat_max - sat_min);
      switch (ipore) {
      case 0:
        mp->d_PorousShellRelPerm[ipore][SHELL_SAT_1] =
            d_s_eff * k_rel * (0.5 / s_eff + 2. * factor2 * factor / ((1.0 - factor) * s_eff * a1));
        break;

      case 1:
        mp->d_PorousShellRelPerm[ipore][SHELL_SAT_2] =
            d_s_eff * k_rel * (0.5 / s_eff + 2. * factor2 * factor / ((1.0 - factor) * s_eff * a1));
        break;

      case 2:
        mp->d_PorousShellRelPerm[ipore][SHELL_SAT_3] =
            d_s_eff * k_rel * (0.5 / s_eff + 2. * factor2 * factor / ((1.0 - factor) * s_eff * a1));
        break;
      }
    }

    break;

  default:
    GOMA_EH(GOMA_ERROR, "Unrecognized Porous Shell Relative Permeability model");
    break;
  }

  return (k_rel);
}
/* END of porous_shell_rel_perm_model */

/*****************************************************************************/

void porous_shell_permeability_model(int ipore) {
  /******************************************************************************
   *
   *  This function computes the permeability  of an open porous shell
   *  Used with the function assemble_porous_shell_saturation.
   *
   *  Kristianto Tjiptowidjojo (tjiptowi@unm.edu) - November 2018
   *
   ******************************************************************************/
  int iext_field = -1;
  int aa, bb;
  double K11 = 0.0, K22 = 0.0, K33 = 0.0, a[DIM][DIM];

  memset(mp->PorousShellPermTensor[ipore], 0, sizeof(double) * DIM * DIM);

  if (mp->PorousShellPermeabilityModel[ipore] == CONSTANT) {
    /* Do nothing */
  } else if (mp->PorousShellPermeabilityModel[ipore] == EXTERNAL_FIELD) {
    iext_field = mp->por_shell_permeability_ext_field_index[ipore];
    GOMA_EH(iext_field, "Porous shell permeability external field not found!");
    mp->PorousShellPermeability[ipore] =
        mp->u_PorousShellPermeability[ipore][0] * fv->external_field[iext_field];
  } else if (mp->PorousShellPermeabilityModel[ipore] == ORTHOTROPIC) {
    /* preload some things */
    K11 = mp->u_PorousShellPermeability[ipore][0];
    K22 = mp->u_PorousShellPermeability[ipore][1];
    K33 = mp->u_PorousShellPermeability[ipore][2];

    /*  N.B.   there are 3 base vectors at the rest state that define the
     * orientation of the material relative to the orthotropic directions
     * If we have sheet laying on the x-y plane, then this matrix of vectors
     * is the Idemfactor.  This will most likely be the case.
     */

    a[0][0] = mp->u_PorousShellPermeability[ipore][3];
    a[0][1] = mp->u_PorousShellPermeability[ipore][4];
    a[0][2] = mp->u_PorousShellPermeability[ipore][5];

    a[1][0] = mp->u_PorousShellPermeability[ipore][6];
    a[1][1] = mp->u_PorousShellPermeability[ipore][7];
    a[1][2] = mp->u_PorousShellPermeability[ipore][8];

    a[2][0] = mp->u_PorousShellPermeability[ipore][9];
    a[2][1] = mp->u_PorousShellPermeability[ipore][10];
    a[2][2] = mp->u_PorousShellPermeability[ipore][11];

    for (aa = 0; aa < DIM; aa++) {
      for (bb = 0; bb < DIM; bb++) {
        mp->PorousShellPermTensor[ipore][aa][bb] =
            a[0][aa] * a[0][bb] * K11 + a[1][aa] * a[1][bb] * K22 + a[2][aa] * a[2][bb] * K33;
      }
    }
  } else {
    GOMA_EH(GOMA_ERROR, "Unrecognized Porous Shell Permeability  model");
  }
}
/* END of porous_shell_permeability_model */

/*****************************************************************************/

/******************************************************************************/
void dynamic_contact_angle_model(double *cos_caU,     // cos(theta) for upper
                                 double *cos_caL,     // cos(theta) for lower
                                 double V,            // wetting line velocity
                                 double *d_cos_caU_V, // sensitivity for upper
                                 double *d_cos_caL_V  // sensitivity for lower
) {
  /*****************************************************************************
   * This function calculates the cosine of the contact angle in the lubrication
   * model and associated shell equations for dynamic contact angles.  Calculates
   * both the lower and upper dynamic contact angles.  Requires the contact
   * line velocity.
   *****************************************************************************/

  /* Parameters */
  dbl theta0, kBT, lambda, vl;
  dbl C1, C2;
  dbl sigma, mu;
  dbl cosU, cosL;
  dbl cosU_V, cosL_V;
  dbl K0;

  /* Load parameters */
  sigma = mp->surface_tension;
  mu = mp->viscosity;

  /* Upper contact angle */
  switch (mp->DcaUFunctionModel) {

  case CONSTANT:
    cosU = cos(mp->dcaU * M_PIE / 180.0);
    cosU_V = 0.0;
    break;

  case DYNAMIC_CA:
    mp->dcaU = mp->u_dcaU_function_constants[0];
    theta0 = mp->u_dcaU_function_constants[0] * M_PIE / 180.0;
    kBT = mp->u_dcaU_function_constants[1];
    lambda = mp->u_dcaU_function_constants[2];
    vl = mp->u_dcaU_function_constants[3];
    K0 = kBT / (mu * vl) * exp(-sigma * pow(lambda, 2) / kBT * (1 + cos(theta0)));
    C1 = 2 * kBT / (sigma * pow(lambda, 2));
    C2 = 1 / (2 * lambda * K0);
    cosU = cos(theta0) - C1 * asinh(C2 * V);
    cosU_V = -C1 * C2 / sqrt(1 + pow(C2 * V, 2));
    break;

  case DYNAMIC_LINEAR_CA:
    mp->dcaU = mp->u_dcaU_function_constants[0];
    theta0 = mp->u_dcaU_function_constants[0] * M_PIE / 180.0;
    kBT = mp->u_dcaU_function_constants[1];
    lambda = mp->u_dcaU_function_constants[2];
    vl = mp->u_dcaU_function_constants[3];
    K0 = kBT / (mu * vl) * exp(-sigma * pow(lambda, 2) / kBT * (1 + cos(theta0)));
    cosU = cos(theta0) - kBT / (K0 * sigma * pow(lambda, 3)) * V;
    cosU_V = -kBT / (K0 * sigma * pow(lambda, 3));
    break;

  default:
    GOMA_EH(GOMA_ERROR, "Wrong upper contact angle model");
    cosU = 0.0;
    cosU_V = 0.0;
    break;
  }
  *cos_caU = cosU;
  *d_cos_caU_V = cosU_V;

  /* Lower contact angle */
  switch (mp->DcaLFunctionModel) {

  case CONSTANT:
    cosL = cos(mp->dcaL * M_PIE / 180.0);
    cosL_V = 0.0;
    break;

  case DYNAMIC_CA:
    mp->dcaL = mp->u_dcaL_function_constants[0];
    theta0 = mp->u_dcaL_function_constants[0] * M_PIE / 180.0;
    kBT = mp->u_dcaL_function_constants[1];
    lambda = mp->u_dcaL_function_constants[2];
    vl = mp->u_dcaL_function_constants[3];
    K0 = kBT / (mu * vl) * exp(-sigma * pow(lambda, 2) / kBT * (1 + cos(theta0)));
    C1 = 2 * kBT / (sigma * pow(lambda, 2));
    C2 = 1 / (2 * lambda * K0);
    cosL = cos(theta0) - C1 * asinh(C2 * V);
    cosL_V = -C1 * C2 / sqrt(1 + pow(C2 * V, 2));
    break;

  case DYNAMIC_LINEAR_CA:
    mp->dcaL = mp->u_dcaL_function_constants[0];
    theta0 = mp->u_dcaL_function_constants[0] * M_PIE / 180.0;
    kBT = mp->u_dcaL_function_constants[1];
    lambda = mp->u_dcaL_function_constants[2];
    vl = mp->u_dcaL_function_constants[3];
    K0 = kBT / (mu * vl) * exp(-sigma * pow(lambda, 2) / kBT * (1 + cos(theta0)));
    cosL = cos(theta0) - kBT / (K0 * sigma * pow(lambda, 3)) * V;
    cosL_V = -kBT / (K0 * sigma * pow(lambda, 3));
    break;

  default:
    GOMA_EH(GOMA_ERROR, "Wrong lower contact angle model");
    cosL = 0.0;
    cosL_V = 0.0;
    break;
  }
  *cos_caL = cosL;
  *d_cos_caL_V = cosL_V;

  return;
}
/*** END OF dynamic_contact_angle_model ***/

/******************************************************************************/
void porous_shell_open_source_model(
    double j_1_2[MDE],                 // Flux between porous layers 1 and 2
    double j_2_3[MDE],                 // Flux between porous layers 2 and 3
    double dj_1_2[MAX_POR_SHELL][MDE], // Sensitivity of the flux between porous layers 1 and 2
    double dj_2_3[MAX_POR_SHELL][MDE]  // Sensitivity of the flux between porous layers 2 and 3
    )
/*****************************************************************************
 * This function calculates inter-layer fluxes amongst porous shell layers.
 * As of now, each layer is assumed to be stacked on top one another.
 * i.e. Layer 3 on top of layer 2 on top of layer 1.
 *
 *
 * Kristianto Tjiptowidjojo   (tjiptowi@unm.edu)  October 2019
 *
 *****************************************************************************/
{
  int ipore, var, j;

  int porous_shell_var[MAX_POR_SHELL];
  porous_shell_var[0] = SHELL_SAT_1;
  porous_shell_var[1] = SHELL_SAT_2;
  porous_shell_var[2] = SHELL_SAT_3;

  dbl H[MAX_POR_SHELL] = {0.0};     // Pore height (vertical)
  dbl kappa[MAX_POR_SHELL] = {0.0}; // Cross permeability
  dbl mu = mp->viscosity;           // Viscosity

  dbl sat_nodes[MAX_POR_SHELL][MDE] = {{0.0}};
  dbl cap_pres_nodes[MAX_POR_SHELL][MDE] = {{0.0}};
  dbl d_cap_pres_nodes_dS[MAX_POR_SHELL][MDE] = {{0.0}};

  dbl Hside_1[MDE] = {0.0};
  dbl dHside_1_dS[MDE] = {0.0};
  dbl Hside_1_square[MDE] = {0.0};
  dbl dHside_1_square_dS[MDE] = {0.0};
  dbl sat_max_1 = 0.99;
  dbl width_1 = 0.24;
  dbl alpha_1 = 0.5 * width_1;
  dbl sat_center_1 = sat_max_1 - alpha_1;
  dbl sat_normalized_1[MDE] = {0.0};

  dbl Hside_2[MDE] = {0.0};
  dbl dHside_2_dS[MDE] = {0.0};
  dbl sat_min_2 = 0.5;
  dbl width_2 = 0.05;
  dbl alpha_2 = 0.5 * width_2;
  dbl sat_center_2 = sat_min_2 - alpha_2;
  dbl sat_normalized_2[MDE] = {0.0};

  dbl Hside_3[MDE] = {0.0};
  dbl dHside_3_dS[MDE] = {0.0};
  dbl Hside_3_square[MDE] = {0.0};
  dbl dHside_3_square_dS[MDE] = {0.0};
  dbl sat_max_3 = 0.99;
  dbl width_3 = 0.24;
  dbl alpha_3 = 0.5 * width_3;
  dbl sat_center_3 = sat_max_3 - alpha_3;
  dbl sat_normalized_3[MDE] = {0.0};

  /* Extract all of the necessary information */
  for (ipore = 0; ipore < pd->Num_Porous_Shell_Eqn; ipore++) {
    var = porous_shell_var[ipore];
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        switch (ipore) {
        case 0:
          sat_nodes[ipore][j] = *esp->sh_sat_1[j];
          sat_normalized_1[j] = sat_nodes[ipore][j] - sat_center_1;
          break;
        case 1:
          sat_nodes[ipore][j] = *esp->sh_sat_2[j];
          sat_normalized_2[j] = sat_nodes[ipore][j] - sat_center_2;
          break;
        case 2:
          sat_nodes[ipore][j] = *esp->sh_sat_3[j];
          sat_normalized_3[j] = sat_nodes[ipore][j] - sat_center_3;
          break;
        }
        cap_pres_nodes[ipore][j] = load_cap_pres(ipore, j, -1, sat_nodes[ipore][j]);
        d_cap_pres_nodes_dS[ipore][j] = mp->d_cap_pres[var];
      }
      H[ipore] = porous_shell_height_model(ipore);
      kappa[ipore] = porous_shell_cross_perm_model(ipore);
    }
  }

  /* Apply heaviside function to deactivate flux at S_1 > 0.99*/

  for (j = 0; j < ei[pg->imtrx]->dof[SHELL_SAT_1]; j++) {
    if (sat_nodes[0][j] >= sat_max_1) {
      Hside_1[j] = 0.0;
      dHside_1_dS[j] = 0.0;
      Hside_1_square[j] = 0.0;
      dHside_1_square_dS[j] = 0.0;
    } else if (sat_nodes[0][j] <= (sat_max_1 - width_1)) {
      Hside_1[j] = 1.0;
      dHside_1_dS[j] = 0.0;
      Hside_1_square[j] = 1.0;
      dHside_1_square_dS[j] = 0.0;
    } else {
      Hside_1[j] = 1.0 - 0.5 * (1. + sat_normalized_1[j] / alpha_1 -
                                sin(M_PIE * sat_normalized_1[j] / alpha_1) / M_PIE);
      dHside_1_dS[j] =
          -0.5 * (1.0 / alpha_1 + cos(M_PIE * sat_normalized_1[j] / alpha_1) / alpha_1);
      Hside_1_square[j] = Hside_1[j] * Hside_1[j];
      dHside_1_square_dS[j] = 2.0 * Hside_1[j] * dHside_1_dS[j];
    }
  }

  /* Apply heaviside function to activate flux at S_2 > S_min*/
  for (j = 0; j < ei[pg->imtrx]->dof[SHELL_SAT_2]; j++) {
    if (sat_nodes[1][j] >= sat_min_2) {
      Hside_2[j] = 1.0;
      dHside_2_dS[j] = 0.0;
    } else if (sat_nodes[1][j] <= (sat_min_2 - width_2)) {
      Hside_2[j] = 0.0;
      dHside_2_dS[j] = 0.0;
    } else {
      Hside_2[j] = 0.5 * (1. + sat_normalized_2[j] / alpha_2 +
                          sin(M_PIE * sat_normalized_2[j] / alpha_2) / M_PIE);
      dHside_2_dS[j] = 0.5 * (1.0 / alpha_2 + cos(M_PIE * sat_normalized_2[j] / alpha_2) / alpha_2);
    }
  }

  /* Apply heaviside function to deactivate flux at S_3 > 0.99*/
  if (pd->e[pg->imtrx][R_SHELL_SAT_3]) {
    for (j = 0; j < ei[pg->imtrx]->dof[SHELL_SAT_3]; j++) {
      if (sat_nodes[2][j] >= sat_max_3) {
        Hside_3[j] = 0.0;
        dHside_3_dS[j] = 0.0;
        Hside_3_square[j] = 0.0;
        dHside_3_square_dS[j] = 0.0;
      } else if (sat_nodes[2][j] <= (sat_max_3 - width_3)) {
        Hside_3[j] = 1.0;
        dHside_3_dS[j] = 0.0;
        Hside_3_square[j] = 1.0;
        dHside_3_square_dS[j] = 0.0;
      } else {
        Hside_3[j] = 1.0 - 0.5 * (1. + sat_normalized_3[j] / alpha_3 -
                                  sin(M_PIE * sat_normalized_3[j] / alpha_3) / M_PIE);
        dHside_3_dS[j] =
            -0.5 * (1.0 / alpha_3 + cos(M_PIE * sat_normalized_3[j] / alpha_3) / alpha_3);
        Hside_3_square[j] = Hside_3[j] * Hside_3[j];
        dHside_3_square_dS[j] = 2.0 * Hside_3[j] * dHside_3_dS[j];
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
  if (dj_1_2 != NULL) {
    for (j = 0; j < ei[pg->imtrx]->dof[SHELL_SAT_1]; j++) {
      dj_1_2[0][j] = Hside_1_square[j] * Hside_2[j] * (kappa[0] / mu) * d_cap_pres_nodes_dS[0][j] /
                         (2.0 * H[0]) +
                     dHside_1_square_dS[j] * Hside_2[j] * (kappa[0] / mu) *
                         (cap_pres_nodes[0][j] - cap_pres_nodes[1][j]) / (2.0 * H[0]);
      dj_1_2[0][j] += Hside_1_square[j] * Hside_2[j] * (kappa[1] / mu) * d_cap_pres_nodes_dS[0][j] /
                          (2.0 * H[1]) +
                      dHside_1_square_dS[j] * Hside_2[j] * (kappa[1] / mu) *
                          (cap_pres_nodes[0][j] - cap_pres_nodes[1][j]) / (2.0 * H[1]);

      dj_1_2[1][j] = Hside_1_square[j] * Hside_2[j] * (kappa[0] / mu) *
                         (-d_cap_pres_nodes_dS[1][j]) / (2.0 * H[0]) +
                     Hside_1_square[j] * dHside_2_dS[j] * (kappa[0] / mu) *
                         (cap_pres_nodes[0][j] - cap_pres_nodes[1][j]) / (2.0 * H[0]);
      dj_1_2[1][j] += Hside_1_square[j] * Hside_2[j] * (kappa[1] / mu) *
                          (-d_cap_pres_nodes_dS[1][j]) / (2.0 * H[1]) +
                      Hside_1_square[j] * dHside_2_dS[j] * (kappa[1] / mu) *
                          (cap_pres_nodes[0][j] - cap_pres_nodes[1][j]) / (2.0 * H[1]);
    }
  }

  if (pd->e[pg->imtrx][R_SHELL_SAT_3]) {
    if (dj_2_3 != NULL) {
      for (j = 0; j < ei[pg->imtrx]->dof[SHELL_SAT_3]; j++) {

        dj_2_3[1][j] = Hside_2[j] * Hside_3_square[j] * (kappa[1] / mu) *
                           (-d_cap_pres_nodes_dS[1][j]) / (2.0 * H[1]) +
                       dHside_2_dS[j] * Hside_3_square[j] * (kappa[1] / mu) *
                           (cap_pres_nodes[2][j] - cap_pres_nodes[1][j]) / (2.0 * H[1]);
        dj_2_3[1][j] += Hside_2[j] * Hside_3_square[j] * (kappa[2] / mu) *
                            (-d_cap_pres_nodes_dS[1][j]) / (2.0 * H[2]) +
                        dHside_2_dS[j] * Hside_3_square[j] * (kappa[2] / mu) *
                            (cap_pres_nodes[2][j] - cap_pres_nodes[1][j]) / (2.0 * H[2]);

        dj_2_3[2][j] = Hside_2[j] * Hside_3_square[j] * (kappa[1] / mu) *
                           d_cap_pres_nodes_dS[2][j] / (2.0 * H[1]) +
                       Hside_2[j] * dHside_3_square_dS[j] * (kappa[1] / mu) *
                           (cap_pres_nodes[2][j] - cap_pres_nodes[1][j]) / (2.0 * H[1]);
        dj_2_3[2][j] += Hside_2[j] * Hside_3_square[j] * (kappa[2] / mu) *
                            d_cap_pres_nodes_dS[2][j] / (2.0 * H[2]) +
                        Hside_2[j] * dHside_3_square_dS[j] * (kappa[2] / mu) *
                            (cap_pres_nodes[2][j] - cap_pres_nodes[1][j]) / (2.0 * H[2]);
      }
    }
  }
  return;
}

int lubrication_fluid_source(double *flux,                           /* Fluid flux */
                             double d_flux[MAX_VARIABLE_TYPES][MDE], /* Fluid flux sensitivities */
                             int *n_dof /* Array containing numbers of DOF */
) {
  /******************************************************************************
   *
   *  This function computes fluid flux into or out from lubrication shell.
   *  Right now it only supports for flux between continuum and lubrication layers
   *  Used with the functions assemble_lubrication and assemble_film.
   *
   *  Kristianto Tjiptowidjojo (tjiptowi@unm.edu) - February 2021
   *
   ******************************************************************************/
  double lub_press = 0.0;
  double kappa = 0.0;
  double mu = 0.0;
  double L = 0.0;
  int j = 0;
  double phi_j;

  /* Determine the appropriate lubrication pressure DOF */
  if (pd->v[pg->imtrx][LUBP]) {
    lub_press = fv->lubp;
  } else if (pd->v[pg->imtrx][SHELL_FILMP]) {
    lub_press = fv->sh_fp;
  } else {
    GOMA_EH(GOMA_ERROR, "Cannot find appropriate lubrication pressure");
  }

  /* Evaluate lubrication fluid source */
  if (mp->LubSourceModel == CONSTANT) {
    *flux = mp->lubsource;
  } else if (mp->LubSourceModel == CONTINUUM_FLUID) {
    kappa = mp->u_lubsource_function_constants[0];
    mu = mp->u_lubsource_function_constants[1];
    L = mp->u_lubsource_function_constants[2];

    *flux = (kappa / mu / L) * (fv->P - lub_press);
    if (d_flux != NULL) {
      for (j = 0; j < n_dof[PRESSURE];
           j++) // Use n_dof since PRESSURE only exists on continuum block
      {
        phi_j = bf[PRESSURE]->phi[j];
        d_flux[PRESSURE][j] = (kappa / mu / L) * phi_j;
      }

      if (pd->v[pg->imtrx][LUBP]) {
        for (j = 0; j < ei[pg->imtrx]->dof[LUBP]; j++) {
          phi_j = bf[LUBP]->phi[j];
          d_flux[LUBP][j] = (kappa / mu / L) * (-phi_j);
        }
      }

      if (pd->v[pg->imtrx][SHELL_FILMP]) {
        for (j = 0; j < ei[pg->imtrx]->dof[SHELL_FILMP]; j++) {
          phi_j = bf[SHELL_FILMP]->phi[j];
          d_flux[SHELL_FILMP][j] = (kappa / mu / L) * (-phi_j);
        }
      }
    }
  } else {
    GOMA_EH(GOMA_ERROR, "Unsupported Lubrication Fluid Source Model");
  }

  return (1);
}
/* END of lubrication_fluid_source */
/*****************************************************************************/
/* END of file mm_std_models_shell.c */
