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

#include "table.h"
#include "mm_as.h"
#include "mm_eh.h"
#include "mm_fill_terms.h"
#include "rf_bc.h"
#include "rf_fem_const.h"

void apply_table_mp(double *func, struct Data_Table *table) {
  int i;
  double interp_val, var1[1], slope, temp;

  if (pd->gv[TEMPERATURE]) {
    temp = fv->T;
  } else {
    temp = upd->Process_Temperature;
  }

  for (i = 0; i < table->columns - 1; i++) {
    if (strcmp(table->t_name[i], "TEMPERATURE") == 0) {
      var1[i] = temp;
    } else if (strcmp(table->t_name[i], "MASS_FRACTION") == 0) {
      var1[i] = fv->c[table->species_eq];
    } else if (strcmp(table->t_name[i], "WAVELENGTH") == 0) {
      const double c0 = 1.0 / (sqrt(upd->Free_Space_Permittivity * upd->Free_Space_Permeability));
      dbl freq = upd->EM_Frequency;
      dbl lambda0 = c0 / freq;
      var1[i] = lambda0;
    } else if (strcmp(table->t_name[i], "FAUX_PLASTIC") == 0) {
      GOMA_EH(GOMA_ERROR, "Oops, I shouldn't be using this call for FAUX_PLASTIC.");
    } else {
      GOMA_EH(GOMA_ERROR, "Material Table Model Error-Unknown Function Column ");
    }
  }

  interp_val = interpolate_table(table, var1, &slope, NULL);
  *func = interp_val;
}

double interpolate_table(struct Data_Table *table, double x[], double *sloper, double dfunc_dx[])
/*
 *      A general routine that uses data supplied in a Data_Table
        structure to compute the ordinate of the table that corresponds
        to the abscissa supplied in x.

        Author: Thomas A. Baer, Org 9111
        Date  : July 16, 1998
        Revised: raroach October 12, 1999

    Parameters:
        table = pointer to Data_Table structure
        x     = array of abscissa(s) where ordinate is to be evaluated
        slope = d(ordinate)/dx

    returns:
        value of ordinate at x and the slope there too.
 */
{
  int i, N, iinter, istartx, istarty;
  double func = 0, y1, y2, y3, y4, tt, uu;
  double *t, *t2, *t3, *f;
  double cee, phi[3], xleft, xright;
  int ngrid1, ngrid2, ngrid3;

  N = table->tablelength - 1;
  t = table->t;
  t2 = table->t2;
  t3 = table->t3;
  f = table->f;
  uu = 0.0;

  switch (table->interp_method) {
  case LINEAR: /* This is knucklehead linear interpolation scheme */

    /*  maybe someday you would prefer a more logical fem-based interpolation */

    if (0) {
      for (i = 0; i < N; i++) {
        cee = (x[0] - t[i]) / (t[i + 1] - t[i]);
        if ((cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0) || (cee > 1.0 && i == N - 1)) {
          phi[0] = -cee + 1.;
          phi[1] = cee;
          table->slope[0] = f[i] * phi[0] + f[i + 1] * phi[1];
          table->slope[1] = f[N + 1 + i] * phi[0] + f[N + 2 + i] * phi[1];
          break;
        }
      }
      table->slope[2] = 0.0;
    } else {
      if (x[0] < t[0]) {
        table->slope[0] = (f[1] - f[0]) / (t[1] - t[0]);
        func = f[0] + (table->slope[0]) * (x[0] - t[0]);
      }

      for (i = 0; x[0] >= t[i] && i < N; i++) {
        if (x[0] >= t[i] && x[0] < t[i + 1]) {
          table->slope[0] = (f[i + 1] - f[i]) / (t[i + 1] - t[i]);
          func = f[i] + (table->slope[0]) * (x[0] - t[i]);
        }
      }
      if (x[0] >= t[N]) {
        table->slope[0] = (f[N] - f[N - 1]) / (t[N] - t[N - 1]);
        func = f[N] + (table->slope[0]) * (x[0] - t[N]);
      }
      table->slope[1] = 0.0;
      *sloper = table->slope[0];
    }
    break;

  case QUADRATIC: /* quadratic lagrangian interpolation scheme */

    if (table->columns == 3) {
      for (i = 0; i < N; i += 2) {
        cee = (x[0] - t[i]) / (t[i + 2] - t[i]);
        if ((cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0) || (cee > 1.0 && i == N - 2)) {
          phi[0] = 2. * cee * cee - 3. * cee + 1.;
          phi[1] = -4. * cee * cee + 4. * cee;
          phi[2] = 2. * cee * cee - cee;
          table->slope[0] = f[i] * phi[0] + f[i + 1] * phi[1] + f[i + 2] * phi[2];
          table->slope[1] = f[N + 1 + i] * phi[0] + f[N + 2 + i] * phi[1] + f[N + 3 + i] * phi[2];
          break;
        }
      }
      table->slope[2] = 0.0;
    } else {
      for (i = 0; i < N; i += 2) {
        cee = (x[0] - t[i]) / (t[i + 2] - t[i]);
        if ((cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0) || (cee > 1.0 && i == N - 2)) {
          phi[0] = 2. * cee * cee - 3. * cee + 1.;
          phi[1] = -4. * cee * cee + 4. * cee;
          phi[2] = 2. * cee * cee - cee;
          func = f[i] * phi[0] + f[i + 1] * phi[1] + f[i + 2] * phi[2];
          phi[0] = 4. * cee - 3.;
          phi[1] = -8. * cee + 4.;
          phi[2] = 4. * cee - 1.;
          table->slope[0] =
              (f[i] * phi[0] + f[i + 1] * phi[1] + f[i + 2] * phi[2]) / (t[i + 2] - t[i]);
          break;
        }
      }
      table->slope[1] = 0.0;
      *sloper = table->slope[0];
    }
    break;

  case QUAD_GP: /* quadratic lagrangian interpolation scheme */

    if (table->columns == 3) {
      for (i = 0; i < N; i += 3) {
        xleft =
            (5. + sqrt(15.)) / 6. * t[i] - 2. / 3. * t[i + 1] + (5. - sqrt(15.)) / 6. * t[i + 2];
        xright =
            (5. - sqrt(15.)) / 6. * t[i] - 2. / 3. * t[i + 1] + (5. + sqrt(15.)) / 6. * t[i + 2];
        cee = (x[0] - xleft) / (xright - xleft);
        if ((cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0) || (cee > 1.0 && i == N - 3)) {
          phi[0] = (20. * cee * cee - 2. * (sqrt(15.) + 10.) * cee + sqrt(15.) + 5.) / 6.;
          phi[1] = (-10. * cee * cee + 10. * cee - 1.) * 2. / 3.;
          phi[2] = (20. * cee * cee + 2. * (sqrt(15.) - 10.) * cee - sqrt(15.) + 5.) / 6.;
          table->slope[0] = f[i] * phi[0] + f[i + 1] * phi[1] + f[i + 2] * phi[2];
          table->slope[1] = f[N + 1 + i] * phi[0] + f[N + 2 + i] * phi[1] + f[N + 3 + i] * phi[2];
          break;
        }
      }
      table->slope[2] = 0.0;
    } else {
      for (i = 0; i < N; i += 3) {
        xleft =
            (5. + sqrt(15.)) / 6. * t[i] - 2. / 3. * t[i + 1] + (5. - sqrt(15.)) / 6. * t[i + 2];
        xright =
            (5. - sqrt(15.)) / 6. * t[i] - 2. / 3. * t[i + 1] + (5. + sqrt(15.)) / 6. * t[i + 2];
        cee = (x[0] - xleft) / (xright - xleft);
        if ((cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0) || (cee > 1.0 && i == N - 3)) {
          phi[0] = (20. * cee * cee - 2. * (sqrt(15.) + 10.) * cee + sqrt(15.) + 5.) / 6.;
          phi[1] = (-10. * cee * cee + 10. * cee - 1.) * 2. / 3.;
          phi[2] = (20. * cee * cee + 2. * (sqrt(15.) - 10.) * cee - sqrt(15.) + 5.) / 6.;
          func = f[i] * phi[0] + f[i + 1] * phi[1] + f[i + 2] * phi[2];
          phi[0] = (20. * cee - sqrt(15.) + 10.) / 3.;
          phi[1] = (-40. * cee + 20.) / 3.;
          phi[2] = (20. * cee + sqrt(15.) - 10.) / 3.;
          table->slope[0] =
              (f[i] * phi[0] + f[i + 1] * phi[1] + f[i + 2] * phi[2]) / (xright - xleft);
          break;
        }
      }
      table->slope[1] = 0.0;
      *sloper = table->slope[0];
    }
    break;

  case BIQUADRATIC: /* biquadratic lagrangian interpolation scheme */

    ngrid1 = table->tablelength / table->ngrid;
    if (table->columns == 5) {
      table->slope[0] = quad_isomap_invert(x[0], x[1], 0, t, t2, NULL, f, ngrid1, table->ngrid, 1,
                                           2, 2, dfunc_dx);
      table->slope[1] = quad_isomap_invert(x[0], x[1], 0, t, t2, NULL, &f[N + 1], ngrid1,
                                           table->ngrid, 1, 2, 2, dfunc_dx);
      table->slope[2] = quad_isomap_invert(x[0], x[1], 0, t, t2, NULL, &f[2 * N + 2], ngrid1,
                                           table->ngrid, 1, 2, 2, dfunc_dx);
    } else {
      func = quad_isomap_invert(x[0], x[1], 0, t, t2, NULL, f, ngrid1, table->ngrid, 1, 2, 2,
                                dfunc_dx);
    }
    break;

  case BILINEAR: /* BILINEAR Interpolation Scheme */
    /*Find Interval of Different Values of Abscissa #1*/
    iinter = 0;
    for (i = 0; i < N; i++) {
      if (table->t[i] != table->t[i + 1] && iinter == 0) {
        iinter = i + 1;
      }
    }
    if (iinter == 1) {
      fprintf(stderr, " MP Interpolate Error - Need more than 1 point per set");
      GOMA_EH(GOMA_ERROR, "Table interpolation not implemented");
    }

    istartx = iinter;

    for (i = iinter; t[i] <= x[0] && i < N - iinter + 1; i = i + iinter) {
      istartx = i + iinter;
    }

    istarty = istartx;

    for (i = istartx + 1; t2[i] <= x[1] && i < istartx + iinter - 1; i++) {
      istarty = i;
    }

    y1 = f[istarty];
    y2 = f[istarty + 1];
    y3 = f[istarty + 1 - iinter];
    y4 = f[istarty - iinter];

    tt = (x[1] - t2[istarty]) / (t2[istarty + 1] - t2[istarty]);
    uu = (x[0] - t[istarty]) / (t[istarty + 1 - iinter] - t[istarty]);

    func = (1. - tt) * (1. - uu) * y1 + tt * (1. - uu) * y2 + tt * uu * y3 + (1. - tt) * uu * y4;

    table->slope[1] = (f[istarty + 1] - f[istarty]) / (t2[istarty + 1] - t2[istarty]);
    table->slope[0] =
        (f[istarty + 1] - f[istarty + 1 - iinter]) / (t[istarty + 1] - t[istarty + 1 - iinter]);
    *sloper = table->slope[0];
    if (dfunc_dx != NULL) {
      dfunc_dx[0] = table->slope[0];
      dfunc_dx[1] = table->slope[1];
    }
    break;

  case TRILINEAR: /* trilinear lagrangian interpolation scheme */

    ngrid1 = table->ngrid;
    ngrid2 = table->ngrid2 / table->ngrid;
    ngrid3 = table->tablelength / table->ngrid2;
    func =
        quad_isomap_invert(x[0], x[1], x[2], t, t2, t3, f, ngrid1, ngrid2, ngrid3, 1, 3, dfunc_dx);
    break;

  case TRIQUADRATIC: /* triquadratic lagrangian interpolation scheme */

    ngrid1 = table->ngrid;
    ngrid2 = table->ngrid2 / table->ngrid;
    ngrid3 = table->tablelength / table->ngrid2;
    func =
        quad_isomap_invert(x[0], x[1], x[2], t, t2, t3, f, ngrid1, ngrid2, ngrid3, 2, 3, dfunc_dx);
    break;

  default:
    GOMA_EH(GOMA_ERROR, "Table interpolation order not implemented");
  }

  return (func);
}

/***************************************************************************/

double
table_distance_search(struct Data_Table *table, double x[], double *sloper, double dfunc_dx[])
/*
 *      A routine to find the shortest distance from a point
 *      to a table represented curve

        Author: Robert B. Secor
        Date  : December 5, 2014
        Revised:

    Parameters:
        table = pointer to Data_Table structure
        x     = array of point coordinates
        slope = d(ordinate)/dx

    returns:
        value of distance
 */
{
  int i, N, iinter, istartx, istarty, basis, ordinate, i_min = 0, elem, iter;
  double func = 0, y1, y2, y3, y4, tt, uu;
  double *t, *t2, *t3, *f;
  double cee, phi[3], xleft, xright, delta, epstol = 0.0001, update, point;
  double phic[3] = {0, 0, 0};
  int ngrid1, ngrid2, ngrid3;
  double dist, dist_min = BIG_PENALTY;

  N = table->tablelength - 1;
  t = table->t;
  t2 = table->t2;
  t3 = table->t3;
  f = table->f;
  if (table->t_index[0] == MESH_POSITION1) {
    basis = 0;
    ordinate = 1;
  } else if (table->t_index[0] == MESH_POSITION2) {
    basis = 1;
    ordinate = 0;
  } else {
    basis = 0;
    ordinate = 1;
  }
  uu = 0.0;

  switch (table->interp_method) {
  case LINEAR: /* This is knucklehead linear interpolation scheme */

    /*  maybe someday you would prefer a more logical fem-based interpolation */

    if (0) {
      for (i = 0; i < N; i++) {
        cee = (x[0] - t[i]) / (t[i + 1] - t[i]);
        if ((cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0) || (cee > 1.0 && i == N - 1)) {
          phi[0] = -cee + 1.;
          phi[1] = cee;
          table->slope[0] = f[i] * phi[0] + f[i + 1] * phi[1];
          table->slope[1] = f[N + 1 + i] * phi[0] + f[N + 2 + i] * phi[1];
          break;
        }
      }
      table->slope[2] = 0.0;
    } else {
      if (x[0] < t[0]) {
        table->slope[0] = (f[1] - f[0]) / (t[1] - t[0]);
        func = f[0] + (table->slope[0]) * (x[0] - t[0]);
      }

      for (i = 0; x[0] >= t[i] && i < N; i++) {
        if (x[0] >= t[i] && x[0] < t[i + 1]) {
          table->slope[0] = (f[i + 1] - f[i]) / (t[i + 1] - t[i]);
          func = f[i] + (table->slope[0]) * (x[0] - t[i]);
        }
      }
      if (x[0] >= t[N]) {
        table->slope[0] = (f[N] - f[N - 1]) / (t[N] - t[N - 1]);
        func = f[N] + (table->slope[0]) * (x[0] - t[N]);
      }
      table->slope[1] = 0.0;
      *sloper = table->slope[0];
    }
    break;

  case QUADRATIC: /* quadratic lagrangian interpolation scheme */

    if (table->columns == 3) {
      for (i = 0; i < N; i += 2) {
        cee = (x[0] - t[i]) / (t[i + 2] - t[i]);
        if ((cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0) || (cee > 1.0 && i == N - 2)) {
          phi[0] = 2. * cee * cee - 3. * cee + 1.;
          phi[1] = -4. * cee * cee + 4. * cee;
          phi[2] = 2. * cee * cee - cee;
          table->slope[0] = f[i] * phi[0] + f[i + 1] * phi[1] + f[i + 2] * phi[2];
          table->slope[1] = f[N + 1 + i] * phi[0] + f[N + 2 + i] * phi[1] + f[N + 3 + i] * phi[2];
          break;
        }
      }
      table->slope[2] = 0.0;
    } else {
      /* search through points for minimum distance*/
      dist_min = dist = sqrt(SQUARE(x[basis] - t[0]) + SQUARE(x[ordinate] - f[0]));
      for (i = 1; i < N; i++) {
        dist = sqrt(SQUARE(x[basis] - t[i]) + SQUARE(x[ordinate] - f[i]));
        if (dist < dist_min) {
          i_min = i;
          dist_min = dist;
        }
      }
      elem = i_min / 2 - 1;
      cee = 0.;
      for (iter = 0; iter < 10; iter++) {
        phi[0] = 2. * cee * cee - 3. * cee + 1.;
        phi[1] = -4. * cee * cee + 4. * cee;
        phi[2] = 2. * cee * cee - cee;
        func = f[2 * elem] * phi[0] + f[2 * elem + 1] * phi[1] + f[2 * elem + 2] * phi[2];
        phic[0] = 4. * cee - 3.;
        phic[1] = -8. * cee + 4.;
        phic[2] = 4. * cee - 1.;
        delta = t[2 * elem + 2] - t[2 * elem];
        *sloper = (f[2 * elem] * phic[0] + f[2 * elem + 1] * phic[1] + f[2 * elem + 2] * phic[2]);
        point = t[2 * elem] + cee * delta;
        update = ((point - x[basis]) * delta + (func - x[ordinate]) * (*sloper)) /
                 (SQUARE(delta) + SQUARE(*sloper));
        cee -= update;
        if (fabs(update) < epstol)
          break;
        if (cee < 0.0 && elem > 0) {
          cee += 1.0;
          elem -= 1;
        }
        if (cee > 1.0 && elem < ((N - 1) / 2 - 1)) {
          cee -= 1.0;
          elem += 1;
        }
      }
      dist_min = sqrt(SQUARE(x[basis] - point) + SQUARE(x[ordinate] - func));
      table->slope[0] = dfunc_dx[basis] = (point - x[basis]) / dist_min;
      dfunc_dx[ordinate] = (func - x[ordinate]) / dist_min;
      table->slope[1] = point;
      table->slope[2] = func;
      *sloper = table->slope[0];
    }
    break;

  case QUAD_GP: /* quadratic lagrangian interpolation scheme */

    if (table->columns == 3) {
      for (i = 0; i < N; i += 3) {
        xleft =
            (5. + sqrt(15.)) / 6. * t[i] - 2. / 3. * t[i + 1] + (5. - sqrt(15.)) / 6. * t[i + 2];
        xright =
            (5. - sqrt(15.)) / 6. * t[i] - 2. / 3. * t[i + 1] + (5. + sqrt(15.)) / 6. * t[i + 2];
        cee = (x[0] - xleft) / (xright - xleft);
        if ((cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0) || (cee > 1.0 && i == N - 3)) {
          phi[0] = (20. * cee * cee - 2. * (sqrt(15.) + 10.) * cee + sqrt(15.) + 5.) / 6.;
          phi[1] = (-10. * cee * cee + 10. * cee - 1.) * 2. / 3.;
          phi[2] = (20. * cee * cee + 2. * (sqrt(15.) - 10.) * cee - sqrt(15.) + 5.) / 6.;
          table->slope[0] = f[i] * phi[0] + f[i + 1] * phi[1] + f[i + 2] * phi[2];
          table->slope[1] = f[N + 1 + i] * phi[0] + f[N + 2 + i] * phi[1] + f[N + 3 + i] * phi[2];
          break;
        }
      }
      table->slope[2] = 0.0;
    } else {
      for (i = 0; i < N; i += 3) {
        xleft =
            (5. + sqrt(15.)) / 6. * t[i] - 2. / 3. * t[i + 1] + (5. - sqrt(15.)) / 6. * t[i + 2];
        xright =
            (5. - sqrt(15.)) / 6. * t[i] - 2. / 3. * t[i + 1] + (5. + sqrt(15.)) / 6. * t[i + 2];
        cee = (x[0] - xleft) / (xright - xleft);
        if ((cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0) || (cee > 1.0 && i == N - 3)) {
          phi[0] = (20. * cee * cee - 2. * (sqrt(15.) + 10.) * cee + sqrt(15.) + 5.) / 6.;
          phi[1] = (-10. * cee * cee + 10. * cee - 1.) * 2. / 3.;
          phi[2] = (20. * cee * cee + 2. * (sqrt(15.) - 10.) * cee - sqrt(15.) + 5.) / 6.;
          func = f[i] * phi[0] + f[i + 1] * phi[1] + f[i + 2] * phi[2];
          phi[0] = (20. * cee - sqrt(15.) + 10.) / 3.;
          phi[1] = (-40. * cee + 20.) / 3.;
          phi[2] = (20. * cee + sqrt(15.) - 10.) / 3.;
          table->slope[0] =
              (f[i] * phi[0] + f[i + 1] * phi[1] + f[i + 2] * phi[2]) / (xright - xleft);
          break;
        }
      }
      table->slope[1] = 0.0;
      *sloper = table->slope[0];
    }
    break;

  case BIQUADRATIC: /* biquadratic lagrangian interpolation scheme */

    ngrid1 = table->tablelength / table->ngrid;
    if (table->columns == 5) {
      table->slope[0] = quad_isomap_invert(x[0], x[1], 0, t, t2, NULL, f, ngrid1, table->ngrid, 1,
                                           2, 2, dfunc_dx);
      table->slope[1] = quad_isomap_invert(x[0], x[1], 0, t, t2, NULL, &f[N + 1], ngrid1,
                                           table->ngrid, 1, 2, 2, dfunc_dx);
      table->slope[2] = quad_isomap_invert(x[0], x[1], 0, t, t2, NULL, &f[2 * N + 2], ngrid1,
                                           table->ngrid, 1, 2, 2, dfunc_dx);
    } else {
      func = quad_isomap_invert(x[0], x[1], 0, t, t2, NULL, f, ngrid1, table->ngrid, 1, 2, 2,
                                dfunc_dx);
    }
    break;

  case BILINEAR: /* BILINEAR Interpolation Scheme */
    /*Find Interval of Different Values of Abscissa #1*/
    iinter = 0;
    for (i = 0; i < N; i++) {
      if (table->t[i] != table->t[i + 1] && iinter == 0) {
        iinter = i + 1;
      }
    }
    if (iinter == 1) {
      fprintf(stderr, " MP Interpolate Error - Need more than 1 point per set");
      GOMA_EH(GOMA_ERROR, "Table interpolation not implemented");
    }

    istartx = iinter;

    for (i = iinter; t[i] <= x[0] && i < N - iinter + 1; i = i + iinter) {
      istartx = i + iinter;
    }

    istarty = istartx;

    for (i = istartx + 1; t2[i] <= x[1] && i < istartx + iinter - 1; i++) {
      istarty = i;
    }

    y1 = f[istarty];
    y2 = f[istarty + 1];
    y3 = f[istarty + 1 - iinter];
    y4 = f[istarty - iinter];

    tt = (x[1] - t2[istarty]) / (t2[istarty + 1] - t2[istarty]);
    uu = (x[0] - t[istarty]) / (t[istarty + 1 - iinter] - t[istarty]);

    func = (1. - tt) * (1. - uu) * y1 + tt * (1. - uu) * y2 + tt * uu * y3 + (1. - tt) * uu * y4;

    table->slope[1] = (f[istarty + 1] - f[istarty]) / (t2[istarty + 1] - t2[istarty]);
    table->slope[0] =
        (f[istarty + 1] - f[istarty + 1 - iinter]) / (t[istarty + 1] - t[istarty + 1 - iinter]);
    *sloper = table->slope[0];
    if (dfunc_dx != NULL) {
      dfunc_dx[0] = table->slope[0];
      dfunc_dx[1] = table->slope[1];
    }
    break;

  case TRILINEAR: /* trilinear lagrangian interpolation scheme */

    ngrid1 = table->ngrid;
    ngrid2 = table->ngrid2 / table->ngrid;
    ngrid3 = table->tablelength / table->ngrid2;
    func =
        quad_isomap_invert(x[0], x[1], x[2], t, t2, t3, f, ngrid1, ngrid2, ngrid3, 1, 3, dfunc_dx);
    break;

  case TRIQUADRATIC: /* triquadratic lagrangian interpolation scheme */

    ngrid1 = table->ngrid;
    ngrid2 = table->ngrid2 / table->ngrid;
    ngrid3 = table->tablelength / table->ngrid2;
    func =
        quad_isomap_invert(x[0], x[1], x[2], t, t2, t3, f, ngrid1, ngrid2, ngrid3, 2, 3, dfunc_dx);
    break;

  default:
    GOMA_EH(GOMA_ERROR, "Table interpolation order not implemented");
  }

  return (dist_min);
}

void apply_table_bc(double *func,
                    double d_func[MAX_VARIABLE_TYPES + MAX_CONC],
                    struct Boundary_Condition *BC_Type,
                    double time_value) {
  /*
        Compute the difference between the current value of the
        field in question and its value obtained by interpolating
        a table of data points.  The abscissa of the table data is
        restricted to one of the three coordinates.  Return the
        the difference and its sensitivity.

        Author : Thomas A. Baer, Org 9111
        Date   : July 16, 1998

    Parameters:
        func = pointer to double that carries back the residual value
        d_func = array of sensitivities of residual wrt to all variables.
        BC_Type = pointer to Boundary Condition structure, i.e. &(BC_Types[bc_input_id]

   */

  int var, basis;
  double slope, interp_val, x_table[2];
  double dfunc_dx[3];

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  basis = BC_Type->table->t_index[0];

  /*  Setup Dummy variable to pass array to interpolation function */
  if (basis != -1)
    x_table[0] = fv->x[basis];
  else
    x_table[0] = time_value;

  if (BC_Type->table->interp_method == BIQUADRATIC || BC_Type->table->interp_method == BILINEAR) {
    x_table[1] = fv->x[BC_Type->table->t_index[1]];
  }

  interp_val = interpolate_table(BC_Type->table, x_table, &slope, dfunc_dx);
  interp_val *= BC_Type->table->yscale;
  slope *= BC_Type->table->yscale;

  var = BC_Type->table->f_index;

  switch (var) {
  case VELOCITY1:
    *func = fv->v[0] - interp_val;
    d_func[var] = 1.0;
    break;
  case VELOCITY2:
    *func = fv->v[1] - interp_val;
    d_func[var] = 1.0;
    break;
  case VELOCITY3:
    *func = fv->v[2] - interp_val;
    d_func[var] = 1.0;
    break;
  case TEMPERATURE:
    *func = fv->T - interp_val;
    d_func[var] = 1.0;
    break;
  case MESH_DISPLACEMENT1:
    *func = fv->d[0] - interp_val;
    d_func[var] = 1.0;
    break;
  case MESH_DISPLACEMENT2:
    *func = fv->d[1] - interp_val;
    d_func[var] = 1.0;
    break;
  case MESH_DISPLACEMENT3:
    *func = fv->d[2] - interp_val;
    d_func[var] = 1.0;
    break;
  case MESH_POSITION1:
    *func = fv->x[0] - interp_val;
    d_func[MESH_DISPLACEMENT1] = 1.0;
    break;
  case MESH_POSITION2:
    *func = fv->x[1] - interp_val;
    d_func[MESH_DISPLACEMENT2] = 1.0;
    break;
  case MESH_POSITION3:
    *func = fv->x[2] - interp_val;
    d_func[MESH_DISPLACEMENT3] = 1.0;
    break;
  case MASS_FRACTION:
    *func = fv->c[BC_Type->species_eq] - interp_val;
    d_func[MAX_VARIABLE_TYPES + BC_Type->species_eq] = 1.0;
    break;
  case PRESSURE:
    *func = fv->P - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS11:
    *func = fv->S[0][0][0] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS12:
    *func = fv->S[0][0][1] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS22:
    *func = fv->S[0][1][1] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS13:
    *func = fv->S[0][0][2] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS23:
    *func = fv->S[0][1][2] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS33:
    *func = fv->S[0][2][2] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS11_1:
    *func = fv->S[1][0][0] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS12_1:
    *func = fv->S[1][0][1] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS22_1:
    *func = fv->S[1][1][1] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS13_1:
    *func = fv->S[1][0][2] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS23_1:
    *func = fv->S[1][1][2] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS33_1:
    *func = fv->S[1][2][2] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS11_2:
    *func = fv->S[2][0][0] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS12_2:
    *func = fv->S[2][0][1] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS22_2:
    *func = fv->S[2][1][1] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS13_2:
    *func = fv->S[2][0][2] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS23_2:
    *func = fv->S[2][1][2] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS33_2:
    *func = fv->S[2][2][2] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS11_3:
    *func = fv->S[3][0][0] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS12_3:
    *func = fv->S[3][0][1] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS22_3:
    *func = fv->S[3][1][1] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS13_3:
    *func = fv->S[3][0][2] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS23_3:
    *func = fv->S[3][1][2] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS33_3:
    *func = fv->S[3][2][2] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS11_4:
    *func = fv->S[4][0][0] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS12_4:
    *func = fv->S[4][0][1] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS22_4:
    *func = fv->S[4][1][1] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS13_4:
    *func = fv->S[4][0][2] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS23_4:
    *func = fv->S[4][1][2] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS33_4:
    *func = fv->S[4][2][2] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS11_5:
    *func = fv->S[5][0][0] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS12_5:
    *func = fv->S[5][0][1] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS22_5:
    *func = fv->S[5][1][1] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS13_5:
    *func = fv->S[5][0][2] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS23_5:
    *func = fv->S[5][1][2] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS33_5:
    *func = fv->S[5][2][2] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS11_6:
    *func = fv->S[6][0][0] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS12_6:
    *func = fv->S[6][0][1] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS22_6:
    *func = fv->S[6][1][1] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS13_6:
    *func = fv->S[6][0][2] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS23_6:
    *func = fv->S[6][1][2] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS33_6:
    *func = fv->S[6][2][2] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS11_7:
    *func = fv->S[7][0][0] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS12_7:
    *func = fv->S[7][0][1] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS22_7:
    *func = fv->S[7][1][1] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS13_7:
    *func = fv->S[7][0][2] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS23_7:
    *func = fv->S[7][1][2] - interp_val;
    d_func[var] = 1.0;
    break;
  case POLYMER_STRESS33_7:
    *func = fv->S[7][2][2] - interp_val;
    d_func[var] = 1.0;
    break;

  default:
    GOMA_EH(GOMA_ERROR, "Variable not yet implemented in TABLE_BC");
    break;
  }

  /* And, at the last, we account for the dependence of the interpolation
   * on the basis coordinate
   */

  if (BC_Type->table->interp_method == BIQUADRATIC || BC_Type->table->interp_method == BILINEAR) {
    d_func[R_MESH1 + BC_Type->table->t_index[0]] -= dfunc_dx[0] * BC_Type->BC_Data_Float[0];
    d_func[R_MESH1 + BC_Type->table->t_index[1]] -= dfunc_dx[1] * BC_Type->BC_Data_Float[0];
  } else {
    if (basis != -1 && pd->e[pg->imtrx][R_MESH1 + basis]) {
      d_func[R_MESH1 + basis] -= slope;
    }
  }
}

double interpolate_table_sat(struct Data_Table *table, double x[DIM])

{
  int i, N, Np1, iinter, istartx, istarty, iad;
  double func = 0, y1, y2, y3, y4, tt, uu;
  double *t, *t2, *f;

  N = table->tablelength - 1;
  Np1 = table->tablelength;
  t = table->t;
  t2 = table->t2;
  f = table->f;
  uu = 0.0;

  switch (table->interp_method) { /* switch */
  case LINEAR:                    /* This is knucklehead linear interpolation scheme */
    /* check if absorb or desorb */
    if (x[1] > 0) {
      iad = 1;
    } else {
      iad = 0;
    }

    if (x[0] < t[0]) {
      table->slope[0] = (f[1 + Np1 * iad] - f[0 + Np1 * iad]) / (t[1] - t[0]);
      func = f[0 + Np1 * iad] + (table->slope[0]) * (x[0] - t[0]);
    }

    for (i = 0; x[0] >= t[i] && i < N; i++) {
      if (x[0] >= t[i] && x[0] < t[i + 1]) {
        table->slope[0] = (f[i + 1 + Np1 * iad] - f[i + Np1 * iad]) / (t[i + 1] - t[i]);
        func = f[i + Np1 * iad] + (table->slope[0]) * (x[0] - t[i]);
      }
    }
    if (x[0] >= t[N]) {
      table->slope[0] = (f[N + Np1 * iad] - f[N - 1 + Np1 * iad]) / (t[N] - t[N - 1]);
      func = f[N + Np1 * iad] + (table->slope[0]) * (x[0] - t[N]);
    }
    table->slope[1] = 0.0;
    break;
  case BILINEAR: /* BILINEAR Interpolation Scheme */
    /* check if absorb or desorb */
    if (x[2] > 0) {
      iad = 1;
    } else {
      iad = 0;
    }

    /*Find Interval of Different Values of Abscissa #1*/
    iinter = 0;
    for (i = 0; i < N; i++) {
      if (table->t[i] != table->t[i + 1] && iinter == 0) {
        iinter = i + 1;
      }
    }
    if (iinter == 1) {
      fprintf(stderr, " MP Interpolate Error - Need more than 1 point per set");
      GOMA_EH(GOMA_ERROR, "Table interpolation not implemented");
    }

    istartx = iinter;

    for (i = 0; t[i] <= x[0] && i < N + 1; i = i + iinter) {
      istartx = i;
    }

    istarty = istartx;

    for (i = istartx + 1; t2[i] <= x[1] && i < istartx + iinter - 1; i++) {
      istarty = i;
    }

    y1 = f[istarty + Np1 * iad];
    y2 = f[istarty + 1 + Np1 * iad];
    y3 = f[istarty + 1 + iinter + Np1 * iad];
    y4 = f[istarty + iinter + Np1 * iad];

    tt = (x[1] - t2[istarty]) / (t2[istarty + 1] - t2[istarty]);
    uu = (x[0] - t[istarty]) / (t[istarty + 1 + iinter] - t[istarty]);

    func = (1. - tt) * (1. - uu) * y1 + tt * (1. - uu) * y2 + tt * uu * y3 + (1. - tt) * uu * y4;

    table->slope[1] =
        (1.0 - uu) * ((f[istarty + 1 + Np1 * iad] - f[istarty + Np1 * iad]) /
                      (t2[istarty + 1] - t2[istarty])) +
        uu * ((f[istarty + 1 + iinter + Np1 * iad] - f[istarty + iinter + Np1 * iad]) /
              (t2[istarty + 1 + iinter] - t2[istarty + iinter]));
    /* slope0 not needed */
    table->slope[0] = (f[istarty + 1 + Np1 * iad] - f[istarty + 1 + iinter + Np1 * iad]) /
                      (t[istarty + 1] - t[istarty + 1 + iinter]);
    break;
  default:
    GOMA_EH(GOMA_ERROR, "Table interpolation order for Saturation not implemented");
  }

  return (func);
}