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
 *$Id: user_pre.c,v 5.1 2007-09-18 18:53:49 prschun Exp $
 */

#include <stdio.h>

/* GOMA include files */

#include "mm_as.h"
#include "mm_eh.h"
#include "std.h"
#include "user_pre.h"

#define GOMA_USER_PRE_C

/*********** R O U T I N E S   I N   T H I S   F I L E ************************
 *
 *       NAME            TYPE            CALLED_BY
 *    user_init_object                    object_distance
 *    user_mat_init                       rf_util.c/
 *    ------------             ---------               --------------
 *************************/
double
user_surf_object(int *int_params MAYBE_UNUSED, dbl *param MAYBE_UNUSED, dbl *r MAYBE_UNUSED) {
  double distance = 0;
  static int warning = 0;

  /*
  int num_params;

  num_params = int_params[0];
  */

  /**********************************************************/

  /* Comment out our remove this line if using this routine */
  if (warning == 0) {
    DPRINTF(stderr, "\n\n#############\n"
                    "# WARNING!! #  No user_defined surface object model implemented"
                    "\n#############\n");
    warning = 1;
  }

  return distance;
} /* End of routine user_init_object */

double user_mat_init(const int var,
                     const int node,
                     const double init_value,
                     const double p[],
                     const double xpt[],
                     const int mn MAYBE_UNUSED,
                     const double var_vals[] MAYBE_UNUSED) {

  double value = 0;
  /* Set this to a nonzero value if using this routine */
  static int warning = -1;

  if (warning == 0) {
    DPRINTF(stderr, "\n\n#############\n"
                    "# WARNING!! #  No user_defined material initialization model implemented"
                    "\n#############\n");
    warning = 1;
  } else if (var == TEMPERATURE && 1) {
    /* 1D Heat Equation with slab centerline at top of print */
    double xpt0[DIM], e_time, distz, T_below, T_init;
    double alpha, speed, ht, sum, xn, exp_arg;
    int n_terms = 4, nt, dir;
    for (dir = 0; dir < DIM; dir++) {
      xpt0[dir] = p[dir];
    }
    alpha = p[DIM];
    speed = p[DIM + 1];
    ht = p[DIM + 2];
    T_below = p[DIM + 3];
    T_init = init_value;
    e_time = 0.; /* negative elapsed time since deposition  */
    for (dir = 0; dir < DIM; dir++) {
      e_time += efv->ext_fld_ndl_val[dir][node] * (xpt[dir] - xpt0[dir]);
    }
    e_time /= speed;

    sum = 0.;
    for (nt = 0; nt < n_terms; nt++) {
      xn = 0.5 + ((double)nt);
      exp_arg = e_time * alpha * SQUARE(xn * M_PIE / ht);
      exp_arg = fmin(exp_arg, 0.0);        /* avoid big exp arguments for t < 0 */
      distz = xpt0[2] + 0.5 * ht - xpt[2]; /* distance from top */
      sum += exp(exp_arg) * cos(M_PIE / ht * distz * xn) * pow(-1., nt) / xn;
    }
    sum *= 2. / M_PIE;
    value = fmin(T_below - (T_below - T_init) * sum, T_init);
    value = fmax(value, T_below); /* Bound to physically realistic values */
  } else if (var == TEMPERATURE && 0) {
    /* 2D Heat Equation of whole slab with T_amb on 3 sides */
    double xpt0[DIM], e_time, distv, distz, T_below, T_init, T_amb;
    double alpha, speed, ht, width, sum, xn, xm, exp_arg, sum1;
    int n_terms = 4, nt, mt, dir;
    for (dir = 0; dir < DIM; dir++) {
      xpt0[dir] = p[dir];
    }
    alpha = p[DIM];
    speed = p[DIM + 1];
    ht = p[DIM + 2];
    width = p[DIM + 3];
    T_amb = p[DIM + 4];
    T_below = p[DIM + 5];
    T_init = init_value;

    e_time = 0.; /* negative elapsed time since deposition  */
    for (dir = 0; dir < DIM; dir++) {
      e_time += efv->ext_fld_ndl_val[dir][node] * (xpt[dir] - xpt0[dir]);
    }
    e_time /= speed;
    /*  "y-coord" value */
    distv = -efv->ext_fld_ndl_val[1][node] * (xpt[0] - xpt0[0]) +
            +efv->ext_fld_ndl_val[0][node] * (xpt[1] - xpt0[1]) + 0.5 * width;
    distz = xpt[2] - xpt0[2] + 0.5 * ht; /*  distance from bottom  */
    /* homogeneous boundary solution */
    sum = 0.;
    for (nt = 0; nt < n_terms; nt++) {
      xn = 0.5 + ((double)nt);
      for (mt = 0; mt < n_terms; mt++) {
        xm = 0.5 + ((double)mt);
        exp_arg = e_time * 2. * alpha * SQUARE(M_PIE) * (SQUARE(xn / ht) + SQUARE(xm / width));
        exp_arg = fmin(exp_arg, 0.0);
        sum += exp(exp_arg) * sin(2. * M_PIE * xn * distz / ht) *
               sin(2. * M_PIE * xm * distv / width) / (xn * xm);
      }
    }
    sum *= 4. / SQUARE(M_PIE);
    /* Add heterogeneous BC on bottom, i.e. SS solution */
    sum1 = 0.;
    for (nt = 0; nt < n_terms; nt++) {
      xn = 0.5 + ((double)nt);
      sum1 += sin(2. * M_PIE * xn * distv / width) * sinh(2. * M_PIE * xn * (1. - distz / ht)) /
              sinh(2. * M_PIE * xn * width / ht);
    }
    sum1 *= 2. / M_PIE * (T_amb - T_below);
    value = fmin(T_amb - (T_amb - T_init) * (sum + sum1), T_init);
    value = fmax(value, T_amb);
  } else if (var >= MESH_DISPLACEMENT1 && var <= MESH_DISPLACEMENT3) {
    value = init_value + p[0] * xpt[0] + p[1] * xpt[1] + p[2] * xpt[2];

  } else {
    GOMA_EH(GOMA_ERROR, "Not a supported usermat initialization condition ");
  }
  return value;
} /* End of routine user_mat_init */

int user_initialize(const int var,
                    double *x,
                    const double init_value,
                    const double p[],
                    const double xpt[],
                    const double var_vals[] MAYBE_UNUSED) {

  double value = 0;
  int i, var_somewhere, idv, mn;
  /* Set this to a nonzero value if using this routine */
  static int warning = -1;

  if (warning == 0) {
    DPRINTF(stderr, "\n\n#############\n"
                    "# WARNING!! #  No user_defined material initialization model implemented"
                    "\n#############\n");
    warning = 1;
  }

  if (var > -1) {
    if (upd->vp[pg->imtrx][var] > -1) {
      for (i = 0; i < DPI_ptr->num_owned_nodes; i++) {
        var_somewhere = FALSE;
        for (mn = 0; mn < upd->Num_Mat; mn++) {
          idv = Index_Solution(i, var, 0, 0, mn, pg->imtrx);
          if (idv != -1) {
            var_somewhere = TRUE;
            break;
          }
        }
        if (var_somewhere) {
          if (var == TEMPERATURE) {
            double dist, alpha, speed, ht, T_below, T_init, sum, xn, exp_arg;
            int n_terms = 4, nt, dir;
            alpha = p[0];
            speed = p[1];
            ht = p[2];
            T_below = p[3];
            T_init = init_value;
            sum = 0.;
            for (nt = 0; nt < n_terms; nt++) {
              xn = 0.5 + ((double)nt);
              dist = 0.;
              for (dir = 0; dir < DIM; dir++) {
                dist += SQUARE(xpt[dir]);
              }
              dist = sqrt(dist);
              exp_arg = dist * alpha * SQUARE(xn * M_PIE / ht) / speed;
              sum += exp(exp_arg) * cos(M_PIE / ht * dist * xn) * 2. / M_PIE * pow(-1., nt) / xn;
            }
            value = fmin(T_below - (T_below - T_init) * sum, T_init);
            value = fmax(value, T_below);
            if (value < 0) {
              fprintf(stderr, "Trouble, negative temperature! %g %g %g %g\n", value, sum, exp_arg,
                      dist);
            }
            x[idv] = value;
          } else {
            GOMA_EH(GOMA_ERROR, "Not a supported user initialization condition ");
          }
        }
      }
    }
  }

  return 1;
} /* End of routine user_initialize */

/*****************************************************************************/
/* End of file user_pre.c*/
/*****************************************************************************/
