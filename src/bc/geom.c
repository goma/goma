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

#include "bc/geom.h"
#include "mm_as.h"
#include "mm_eh.h"
#include "rf_fem.h"
#include "user_bc.h"

void moving_plane(int ielem_dim, double *func, double d_func[], dbl *aa, double time) {
  fplane(ielem_dim, func, d_func, aa);
}

void fplane(int ielem_dim,
            double *func,
            double d_func[], /* dimensioned [MAX_VARIABLE_TYPES+MAX_CONC] */
            dbl *aa)         /*  function parameters from data card  */
{
  /**************************** EXECUTION BEGINS *******************************/
  int i;

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  *func = (double)aa[3];
  for (i = 0; i < ielem_dim; i++) {
    d_func[MESH_DISPLACEMENT1 + i] = aa[i];
    *func += (double)aa[i] * fv->x[i];
  }

} /* END of routine fplane                                                   */

void f_fillet(const int ielem_dim,
              double *func,
              double d_func[],     /* dimensioned [MAX_VARIABLE_TYPES+MAX_CONC] */
              const double *p,     /*  function parameters from data card  */
              const int num_const) /* number of passed parameters   */
{
  /**************************** EXECUTION BEGINS *******************************/
  double pt[DIM], side_th[2], rad, center[DIM], alpha, theta_mid, theta_avg, tmp;
  double theta, siderad = -1., circ[DIM] = {0., 0., 0.}, dsign, beta = 0, theta_side[2];
  double cham_ang = -1., cham_rotn = 0., gamma = 0.;
  int iside = 0, chamfer = 0, i, dim = 2;

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  if (num_const < 5) {
    GOMA_EH(GOMA_ERROR, "Need at least 5 parameters for 2D fillet geometry bc!\n");
  }

  pt[0] = p[0];
  pt[1] = p[1];
  side_th[0] = p[2];
  side_th[1] = p[3];
  rad = p[4];
  if (num_const >= 7) {
    siderad = p[5];
    iside = ((int)p[6]);
  }
  if (num_const >= 8) {
    chamfer = ((int)p[7]);
  }
  if (num_const >= 9) {
    cham_ang = p[8];
  }
  if (ielem_dim > dim)
    GOMA_WH(-1, "FILLET_BC: Only z-invariant geometry available for now.\n");

  /**  find center of die face  **/

  alpha = 0.5 * (side_th[1] - side_th[0]);
  theta_avg = 0.5 * (side_th[1] + side_th[0]);
  if (iside) {
    dsign = ((double)(3 - 2 * iside));
    beta = side_th[2 - iside] - dsign * asin((rad - siderad * cos(2. * alpha)) / (rad + siderad));
    for (i = 0; i < dim; i++) {
      circ[i] = pt[i] + dsign * siderad * sin(side_th[iside - 1] - 0.5 * M_PIE * i);
      center[i] = circ[i] + (rad + siderad) * cos(beta - 0.5 * M_PIE * i);
    }
  } else {
    for (i = 0; i < dim; i++) {
      center[i] = pt[i] + (rad / sin(alpha)) * cos(theta_avg - 0.5 * M_PIE * i);
    }
  }

  /**   compute angle of point on curve from arc center **/

  theta = atan2(fv->x[1] - center[1], fv->x[0] - center[0]);
  theta = theta > side_th[1] - 1.5 * M_PIE ? theta : theta + 2 * M_PIE;
  theta_mid = atan2(center[1] - pt[1], center[0] - pt[0]);
  theta_side[0] = side_th[0] - 0.5 * M_PIE;
  theta_side[1] = side_th[1];
  if (iside == 1)
    theta_side[0] = beta - M_PIE;
  if (iside == 2)
    theta_side[1] = beta + 0.5 * M_PIE;
  if (cham_ang > 0) {
    cham_rotn = cham_ang - (theta_mid + 0.5 * M_PIE);
    gamma = atan(-2. * cos(alpha) * sin(cham_rotn) / cos(alpha + cham_rotn));
  }

  /**  use different f depending on theta  **/

  if ((theta_side[0] - gamma) <= theta && theta <= theta_avg) {
    if (iside == 1) {
      *func = SQUARE(fv->x[0] - circ[0]) + SQUARE(fv->x[1] - circ[1]) - SQUARE(siderad);
      d_func[MESH_DISPLACEMENT1] = 2. * (fv->x[0] - circ[0]);
      d_func[MESH_DISPLACEMENT2] = 2. * (fv->x[1] - circ[1]);
    } else {
      *func = (fv->x[1] - pt[1]) * cos(side_th[0]) - (fv->x[0] - pt[0]) * sin(side_th[0]);
      d_func[MESH_DISPLACEMENT1] = -sin(side_th[0]);
      d_func[MESH_DISPLACEMENT2] = cos(side_th[0]);
    }

  } else if (theta_avg <= theta && (theta - 0.5 * M_PIE) <= (theta_side[1] + 0 * gamma)) {
    if (iside == 2) {
      *func = SQUARE(fv->x[0] - circ[0]) + SQUARE(fv->x[1] - circ[1]) - SQUARE(siderad);
      d_func[MESH_DISPLACEMENT1] = 2. * (fv->x[0] - circ[0]);
      d_func[MESH_DISPLACEMENT2] = 2. * (fv->x[1] - circ[1]);
    } else {
      *func = (fv->x[1] - pt[1]) * cos(side_th[1]) - (fv->x[0] - pt[0]) * sin(side_th[1]);
      d_func[MESH_DISPLACEMENT1] = -sin(side_th[1]);
      d_func[MESH_DISPLACEMENT2] = cos(side_th[1]);
    }

  } else {
    if (chamfer) {
      tmp = theta_mid + cham_rotn;
      *func = (fv->x[1] - center[1]) * sin(tmp) + (fv->x[0] - center[0]) * cos(tmp) +
              rad * sin(alpha - cham_rotn);
      d_func[MESH_DISPLACEMENT1] = cos(tmp);
      d_func[MESH_DISPLACEMENT2] = sin(tmp);
    } else {
      *func = SQUARE(fv->x[0] - center[0]) + SQUARE(fv->x[1] - center[1]) - SQUARE(rad);
      d_func[MESH_DISPLACEMENT1] = 2. * (fv->x[0] - center[0]);
      d_func[MESH_DISPLACEMENT2] = 2. * (fv->x[1] - center[1]);
    }
  }

  if (ielem_dim == 3)
    d_func[MESH_DISPLACEMENT3] = 0.0;

} /* END of routine f_fillet                                                   */

void f_double_rad(const int ielem_dim,
                  double *func,
                  double d_func[],     /* dimensioned [MAX_VARIABLE_TYPES+MAX_CONC] */
                  const double *p,     /*  function parameters from data card  */
                  const int num_const) /* number of passed parameters   */
{
  /**************************** EXECUTION BEGINS *******************************/
  double xpt1, ypt1, theta1, rad1, xcen1, ycen1, alpha1;
  double xpt2, ypt2, theta2, rad2, xcen2, ycen2, alpha2;
  double theta1m, theta2m, th1, th2, th2t, curv_mid, rad_curv;
  double beta = 0, xcirc, ycirc, dist1, dist2, dist_mid;
  int is_curved = 0;

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  if (num_const < 8)
    GOMA_EH(GOMA_ERROR, "Need at least 8 parameters for Double Rad lip geometry bc!\n");

  xpt1 = p[0];
  ypt1 = p[1];
  theta1 = p[2];
  rad1 = p[3];
  xpt2 = p[4];
  ypt2 = p[5];
  theta2 = p[6];
  rad2 = p[7];
  if (num_const >= 8) {
    curv_mid = p[8];
  } else {
    curv_mid = 0.0;
  }

  is_curved = DOUBLE_NONZERO(curv_mid);

  /*  slope of middle line                */

  theta1m = atan2(ypt2 - ypt1, xpt2 - xpt1);
  theta1m = theta1m > theta1 ? theta1m : theta1m + 2 * M_PIE;
  theta2m = atan2(ypt1 - ypt2, xpt1 - xpt2);
  alpha1 = 0.5 * (theta1m - theta1);
  alpha2 = 0.5 * (theta2 - theta2m);

  xcen1 = xpt1 + (rad1 / sin(alpha1)) * cos(theta1 + alpha1);
  ycen1 = ypt1 + (rad1 / sin(alpha1)) * sin(theta1 + alpha1);
  xcen2 = xpt2 + (rad2 / sin(alpha2)) * cos(theta2m + alpha2);
  ycen2 = ypt2 + (rad2 / sin(alpha2)) * sin(theta2m + alpha2);

  if (is_curved) {
    rad_curv = 1. / curv_mid;
    dist_mid = sqrt(SQUARE(xpt1 - xpt2) + SQUARE(ypt1 - ypt2));
    beta = asin(0.5 * dist_mid * curv_mid);
    xcirc = 0.5 * (xpt1 + xpt2) + rad_curv * cos(beta) * sin(theta1m);
    ycirc = 0.5 * (ypt1 + ypt2) - rad_curv * cos(beta) * cos(theta1m);
    /**   Shift fillet centers based on curvature  **/
    /*  Using approximate distance from 90 degree corner for simplicity */

    dist1 = rad1 + sqrt((rad_curv - 0.5 * dist_mid) * (rad_curv + 0.5 * dist_mid - 2 * rad1)) -
            sqrt((rad_curv - 0.5 * dist_mid) * (rad_curv + 0.5 * dist_mid));
    dist2 = rad2 + sqrt((rad_curv - 0.5 * dist_mid) * (rad_curv + 0.5 * dist_mid - 2 * rad2)) -
            sqrt((rad_curv - 0.5 * dist_mid) * (rad_curv + 0.5 * dist_mid));
#if 0
     dist1 = 0;  dist2 = 0;
#endif
    xcen1 -= dist1 * cos(theta1);
    ycen1 -= dist1 * sin(theta1);
    xcen2 -= dist2 * cos(theta2);
    ycen2 -= dist2 * sin(theta2);
#if 0
fprintf(stderr,"arc distances %g %g \n",dist1,dist2);
fprintf(stderr,"rads %g %g %g %g\n",rad1, rad2, rad_curv,beta);
fprintf(stderr,"circle %g %g %g %g\n",xcirc,ycirc,xcen2, ycen2);
#endif
  }

  /**   compute angle of point on curve from arc center **/

  th1 = atan2(fv->x[1] - ycen1, fv->x[0] - xcen1);
  th2 = atan2(fv->x[1] - ycen2, fv->x[0] - xcen2);
  th2t = th2 > 0.0 ? th2 : th2 + 2 * M_PIE;

  /**  use different f depending on theta  **/

  if ((theta1 - 0.5 * M_PIE) <= th1 && th1 <= (theta1 + alpha1)) {
    *func = (fv->x[1] - ypt1) * cos(theta1) - (fv->x[0] - xpt1) * sin(theta1);
    d_func[MESH_DISPLACEMENT1] = -sin(theta1);
    d_func[MESH_DISPLACEMENT2] = cos(theta1);
    /*fprintf(stderr,"DR case 1 %g %g %g %g %g %g\n",*func,fv->x[0],fv->x[1],th1,th2,th2t);*/
  } else if ((theta2 - alpha2) <= th2t && (th2t - 0.5 * M_PIE) <= theta2) {
    *func = (fv->x[1] - ypt2) * cos(theta2) - (fv->x[0] - xpt2) * sin(theta2);
    d_func[MESH_DISPLACEMENT1] = -sin(theta2);
    d_func[MESH_DISPLACEMENT2] = cos(theta2);
    /*fprintf(stderr,"DR case 2 %g %g %g %g %g %g\n",*func,fv->x[0],fv->x[1],th1,th2,th2t);*/
  } else if (theta2m <= (th1 + 0.5 * M_PIE - beta) && th1 <= (theta1 - 0.5 * M_PIE)) {
    *func = SQUARE(fv->x[0] - xcen1) + SQUARE(fv->x[1] - ycen1) - SQUARE(rad1);
    d_func[MESH_DISPLACEMENT1] = 2. * (fv->x[0] - xcen1);
    d_func[MESH_DISPLACEMENT2] = 2. * (fv->x[1] - ycen1);
    /*fprintf(stderr,"DR case 3 %g %g %g %g %g %g\n",*func,fv->x[0],fv->x[1],th1,th2,th2t);*/
  } else if ((theta2m - 0.5 * M_PIE - beta) <= th2 &&
             ((th1 + 0.5 * M_PIE - beta) <= theta2m ||
              (th1 + 0.5 * M_PIE - beta) <= (theta1m - M_PIE))) {
    if (is_curved) {
      *func = SQUARE(fv->x[0] - xcirc) + SQUARE(fv->x[1] - ycirc) - SQUARE(rad_curv);
      d_func[MESH_DISPLACEMENT1] = 2. * (fv->x[0] - xcirc);
      d_func[MESH_DISPLACEMENT2] = 2. * (fv->x[1] - ycirc);
    } else {
      *func = (fv->x[1] - ypt1) * cos(theta1m) - (fv->x[0] - xpt1) * sin(theta1m);
      d_func[MESH_DISPLACEMENT1] = -sin(theta1m);
      d_func[MESH_DISPLACEMENT2] = cos(theta1m);
    }
    /*fprintf(stderr,"DR case 4 %g %g %g %g %g %g %g
     * %g\n",*func,fv->x[0],fv->x[1],th1,th2,th2t,theta1m,theta2m);*/
  } else if ((th2t - 0.5 * M_PIE) >= theta2 && (theta2m - 0.5 * M_PIE - beta) >= th2) {
    *func = SQUARE(fv->x[0] - xcen2) + SQUARE(fv->x[1] - ycen2) - SQUARE(rad2);
    d_func[MESH_DISPLACEMENT1] = 2. * (fv->x[0] - xcen2);
    d_func[MESH_DISPLACEMENT2] = 2. * (fv->x[1] - ycen2);
    /*fprintf(stderr,"DR case 5 %g %g %g %g %g %g\n",*func,fv->x[0],fv->x[1],th1,th2,th2t);*/
  } else {
    fprintf(stderr, "Double Rad case not found... %g %g %g %g %g\n", fv->x[0], fv->x[1], th1, th2,
            th2t);
  }

  if (ielem_dim == 3)
    d_func[MESH_DISPLACEMENT3] = 0.0;

} /* END of routine f_double_rad                                             */

void rotate_line(dbl origin[DIM], dbl point[DIM], dbl angle, dbl xrot[DIM]) {
  dbl ca = cos(angle);
  dbl sa = sin(angle);
  xrot[0] = origin[0] + ca * (point[0] - origin[0]) - sa * (point[1] - origin[1]);
  xrot[1] = origin[1] + sa * (point[0] - origin[0]) + ca * (point[1] - origin[1]);
}

void fillet_center(dbl origin[DIM],
                   dbl p1[DIM],
                   dbl p2[DIM],
                   dbl rad,
                   dbl xcen[DIM],
                   dbl start1[DIM],
                   dbl start2[DIM]) {
  dbl line1[DIM], line2[DIM];

  line1[0] = origin[0] - p1[0];
  line1[1] = origin[1] - p1[1];

  line2[0] = origin[0] - p2[0];
  line2[1] = origin[1] - p2[1];

  dbl angle1 = atan2(line1[1], line1[0]);
  dbl angle2 = atan2(line2[1], line2[0]);

  start1[0] = origin[0] - rad * cos(angle1);
  start1[1] = origin[1] - rad * sin(angle1);

  start2[0] = origin[0] - rad * cos(angle2);
  start2[1] = origin[1] - rad * sin(angle2);

  dbl x1[DIM] = {origin[0] - rad * cos(angle1), origin[1] - rad * sin(angle1)};
  dbl x2[DIM] = {origin[0] - rad * cos(angle2), origin[1] - rad * sin(angle2)};

  dbl mid[DIM] = {(x1[0] + x2[0]) / 2.0, (x1[1] + x2[1]) / 2.0};

  dbl midline[DIM] = {origin[0] - mid[0], origin[1] - mid[1]};

  dbl angle_mid = atan2(midline[1], midline[0]);

  dbl midlength = sqrt(SQUARE(midline[0]) + SQUARE(midline[1]));

  xcen[0] = mid[0] - midlength * cos(angle_mid);
  xcen[1] = mid[1] - midlength * sin(angle_mid);
}

bool in_between_vectors_2D(const dbl v1[DIM], const dbl v2[DIM], const dbl v3[DIM]) {
  // z component of cross product v1 x v3 * v1 x v2
  dbl cross1 = (v1[1] * v3[0] - v1[0] * v3[1]) * (v1[1] * v2[0] - v1[0] * v2[1]);
  // z component of cross product v2 x v3 * v2 x v1
  dbl cross2 = (v2[1] * v3[0] - v2[0] * v3[1]) * (v2[1] * v1[0] - v2[0] * v1[1]);
  // both have same sign then the point is between the two vectors
  return cross1 * cross2 > 0.0;
}

bool near_circle(const dbl xcen[DIM], const dbl rad, const dbl point[DIM], const dbl tol) {
  // check if point is within the circle radius plus tolerance
  dbl dist = sqrt(SQUARE(point[0] - xcen[0]) + SQUARE(point[1] - xcen[1]));
  return (dist <= rad + tol);
}

/*****************************************************************************/
/* This function is used to create a double fillet geometry boundary condition.
 * It is used in the case of a die with two fillets on the edges.
 * This is a modified version of the f_double_rad function with more bookkeeping.
 * f_double_rad is being kept for backward compatibility.
 *
 * Original modifications done by Chance.
 */
void f_double_fillet(const int ielem_dim,
                     double *func,
                     double d_func[],     /* dimensioned [MAX_VARIABLE_TYPES+MAX_CONC] */
                     const double *p,     /*  function parameters from data card  */
                     const int num_const) /* number of passed parameters   */
{
  /**************************** EXECUTION BEGINS *******************************/
  double xpt1, ypt1, theta1, rad1, xcen1, ycen1, alpha1;
  double xpt2, ypt2, theta2, rad2, xcen2, ycen2, alpha2;
  double xi1, yi1, xf1, yf1, xi2, yi2, xf2, yf2, xint, yint;
  double s1, f1, s2, f2, th1ub;
  double theta1m, theta2m, th1, th2, th2t, curv_mid, rad_curv;
  double beta = 0, xcirc, ycirc, dist1, dist2, dist_mid;
  double atol = 1.0e-9;
  int is_curved = 0;

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  if (num_const < 8)
    GOMA_EH(GOMA_ERROR, "Need at least 8 parameters for Double Rad lip geometry bc!\n");
  xpt1 = p[0];
  ypt1 = p[1];
  theta1 = p[2];
  rad1 = p[3];
  xpt2 = p[4];
  ypt2 = p[5];
  theta2 = p[6];
  rad2 = p[7];
  if (num_const >= 8) {
    curv_mid = p[8];
  } else {
    curv_mid = 0.0;
  }

  is_curved = DOUBLE_NONZERO(curv_mid);

  /*  slope of middle line                */

  theta1m = atan2(ypt2 - ypt1, xpt2 - xpt1);
  theta1m = theta1m > theta1 ? theta1m : theta1m + 2 * M_PIE;
  theta2m = atan2(ypt1 - ypt2, xpt1 - xpt2);
  alpha1 = 0.5 * (theta1m - theta1);
  alpha2 = 0.5 * (theta2 - theta2m);

  xcen1 = xpt1 + (rad1 / sin(alpha1)) * cos(theta1 + alpha1);
  ycen1 = ypt1 + (rad1 / sin(alpha1)) * sin(theta1 + alpha1);
  xcen2 = xpt2 + (rad2 / sin(alpha2)) * cos(theta2m + alpha2);
  ycen2 = ypt2 + (rad2 / sin(alpha2)) * sin(theta2m + alpha2);

  if (is_curved) {
    rad_curv = 1. / curv_mid;
    dist_mid = sqrt(SQUARE(xpt1 - xpt2) + SQUARE(ypt1 - ypt2));
    beta = asin(0.5 * dist_mid * curv_mid);
    xcirc = 0.5 * (xpt1 + xpt2) + rad_curv * cos(beta) * sin(theta1m);
    ycirc = 0.5 * (ypt1 + ypt2) - rad_curv * cos(beta) * cos(theta1m);
    /**   Shift fillet centers based on curvature  **/
    /*  Using approximate distance from 90 degree corner for simplicity */

    dist1 = rad1 + sqrt((rad_curv - 0.5 * dist_mid) * (rad_curv + 0.5 * dist_mid - 2 * rad1)) -
            sqrt((rad_curv - 0.5 * dist_mid) * (rad_curv + 0.5 * dist_mid));
    dist2 = rad2 + sqrt((rad_curv - 0.5 * dist_mid) * (rad_curv + 0.5 * dist_mid - 2 * rad2)) -
            sqrt((rad_curv - 0.5 * dist_mid) * (rad_curv + 0.5 * dist_mid));
#if 0
     dist1 = 0;  dist2 = 0;
#endif
    xcen1 -= dist1 * cos(theta1);
    ycen1 -= dist1 * sin(theta1);
    xcen2 -= dist2 * cos(theta2);
    ycen2 -= dist2 * sin(theta2);
#if 0
fprintf(stderr,"arc distances %g %g \n",dist1,dist2);
fprintf(stderr,"rads %g %g %g %g\n",rad1, rad2, rad_curv,beta);
fprintf(stderr,"circle %g %g %g %g\n",xcirc,ycirc,xcen2, ycen2);
#endif
  }

  xf1 = xpt1 + rad1 / tan(alpha1) * cos(theta1);
  yf1 = ypt1 + rad1 / tan(alpha1) * sin(theta1);
  xi1 = xpt1 - rad1 / tan(alpha1) * cos(theta1m - M_PIE);
  yi1 = ypt1 - rad1 / tan(alpha1) * sin(theta1m - M_PIE);
  xi2 = xcen2 + rad2 * cos(theta2 + 0.5 * M_PIE);
  yi2 = ycen2 + rad2 * sin(theta2 + 0.5 * M_PIE);
  xf2 = xpt2 + rad2 / tan(alpha2) * cos(theta2m);
  yf2 = ypt2 + rad2 / tan(alpha2) * sin(theta2m);

  s1 = atan2(yi1 - ycen1, xi1 - xcen1) > 0 ? atan2(yi1 - ycen1, xi1 - xcen1)
                                           : atan2(yi1 - ycen1, xi1 - xcen1) + 2.0 * M_PIE;
  f1 = atan2(yf1 - ycen1, xf1 - xcen1) > 0 ? atan2(yf1 - ycen1, xf1 - xcen1)
                                           : atan2(yf1 - ycen1, xf1 - xcen1) + 2.0 * M_PIE;
  s2 = atan2(yi2 - ycen2, xi2 - xcen2) > 0 ? atan2(yi2 - ycen2, xi2 - xcen2)
                                           : atan2(yi2 - ycen2, xi2 - xcen2) + 2.0 * M_PIE;
  f2 = atan2(yf2 - ycen2, xf2 - xcen2) > 0 ? atan2(yf2 - ycen2, xf2 - xcen2)
                                           : atan2(yf2 - ycen2, xf2 - xcen2) + 2.0 * M_PIE;

  /**   compute angle of point on curve from arc center **/

  th1 = atan2(fv->x[1] - ycen1, fv->x[0] - xcen1);
  th2 = atan2(fv->x[1] - ycen2, fv->x[0] - xcen2);
  th2t = th2 > 0.0 ? th2 : th2 + 2 * M_PIE;

  if (sqrt((theta1 - 0.5 * M_PIE) * (theta1 - 0.5 * M_PIE)) < atol ||
      sqrt((theta1 + 0.5 * M_PIE) * (theta1 + 0.5 * M_PIE)) < atol) {
    if (sqrt((theta2 - 0.5 * M_PIE) * (theta2 - 0.5 * M_PIE)) < atol ||
        sqrt((theta2 + 0.5 * M_PIE) * (theta2 + 0.5 * M_PIE)) < atol) {
      th1ub = 0.5 * M_PIE;
    } else {
      xint = xpt1;
      yint = (xint - xpt2) * tan(theta2) + ypt2;
      th1ub = atan2(ycen1 - yint, xcen1 - xint);
    }
  } else if (sqrt((theta2 - 0.5 * M_PIE) * (theta2 - 0.5 * M_PIE)) < atol ||
             sqrt((theta2 + 0.5 * M_PIE) * (theta2 + 0.5 * M_PIE)) < atol) {
    xint = xpt2;
    yint = (xint - xpt1) * tan(theta1) + ypt1;
    th1ub = atan2(ycen1 - yint, xcen1 - xint);
  } else if (sqrt((theta2 - theta1) * (theta2 - theta1)) < atol) {
    th1ub = theta1;
  } else {
    th1ub = theta1;
  }

  /**  use different f depending on theta  **/

  if ((theta1 - 0.5 * M_PIE) <= th1 && (th1 <= th1ub)) {
    *func = (fv->x[1] - ypt1) * cos(theta1) - (fv->x[0] - xpt1) * sin(theta1);
    d_func[MESH_DISPLACEMENT1] = -sin(theta1);
    d_func[MESH_DISPLACEMENT2] = cos(theta1);
    /*fprintf(stderr,"DR case 1 %g %g %g %g %g %g\n",*func,fv->x[0],fv->x[1],th1,th2,th2t);*/
  } else if ((theta2 - alpha2) <= th2t && (th2t - 0.5 * M_PIE) <= theta2) {
    *func = (fv->x[1] - ypt2) * cos(theta2) - (fv->x[0] - xpt2) * sin(theta2);
    d_func[MESH_DISPLACEMENT1] = -sin(theta2);
    d_func[MESH_DISPLACEMENT2] = cos(theta2);
    /*fprintf(stderr,"DR case 2 %g %g %g %g %g %g\n",*func,fv->x[0],fv->x[1],th1,th2,th2t);*/
  } else if (theta2m <= (th1 + 0.5 * M_PIE - beta) && th1 <= (theta1 - 0.5 * M_PIE)) {
    *func = SQUARE(fv->x[0] - xcen1) + SQUARE(fv->x[1] - ycen1) - SQUARE(rad1);
    d_func[MESH_DISPLACEMENT1] = 2. * (fv->x[0] - xcen1);
    d_func[MESH_DISPLACEMENT2] = 2. * (fv->x[1] - ycen1);
    /*fprintf(stderr,"DR case 3 %g %g %g %g %g %g\n",*func,fv->x[0],fv->x[1],th1,th2,th2t);*/
  } else if ((theta2m - 0.5 * M_PIE - beta) <= th2 &&
             ((th1 + 0.5 * M_PIE - beta) <= theta2m ||
              (th1 + 0.5 * M_PIE - beta) <= (theta1m - M_PIE))) {
    if (is_curved) {
      *func = SQUARE(fv->x[0] - xcirc) + SQUARE(fv->x[1] - ycirc) - SQUARE(rad_curv);
      d_func[MESH_DISPLACEMENT1] = 2. * (fv->x[0] - xcirc);
      d_func[MESH_DISPLACEMENT2] = 2. * (fv->x[1] - ycirc);
    } else {
      *func = (fv->x[1] - ypt1) * cos(theta1m) - (fv->x[0] - xpt1) * sin(theta1m);
      d_func[MESH_DISPLACEMENT1] = -sin(theta1m);
      d_func[MESH_DISPLACEMENT2] = cos(theta1m);
    }
    /*fprintf(stderr,"DR case 4 %g %g %g %g %g %g %g
     * %g\n",*func,fv->x[0],fv->x[1],th1,th2,th2t,theta1m,theta2m);*/
  } else if ((th2t - 0.5 * M_PIE) >= theta2 && ((theta2m - 0.5 * M_PIE - beta) >= th2)) {
    *func = SQUARE(fv->x[0] - xcen2) + SQUARE(fv->x[1] - ycen2) - SQUARE(rad2);
    d_func[MESH_DISPLACEMENT1] = 2. * (fv->x[0] - xcen2);
    d_func[MESH_DISPLACEMENT2] = 2. * (fv->x[1] - ycen2);
    /*fprintf(stderr,"DR case 5 %g %g %g %g %g %g\n",*func,fv->x[0],fv->x[1],th1,th2,th2t);*/
  } else {
    if (s2 <= th2 && th2 <= f2) {
      *func = SQUARE(fv->x[0] - xcen2) + SQUARE(fv->x[1] - ycen2) - SQUARE(rad2);
      d_func[MESH_DISPLACEMENT1] = 2. * (fv->x[0] - xcen2);
      d_func[MESH_DISPLACEMENT2] = 2. * (fv->x[1] - ycen2);
      /*fprintf(stderr,"DR case 5 %g %g %g %g %g %g\n",*func,fv->x[0],fv->x[1],th1,th2,th2t);*/
    } else if (s1 <= th1 && th1 <= f1) {
      *func = (fv->x[1] - ypt1) * cos(theta1m) - (fv->x[0] - xpt1) * sin(theta1m);
      d_func[MESH_DISPLACEMENT1] = -sin(theta1m);
      d_func[MESH_DISPLACEMENT2] = cos(theta1m);
    } else if (f1 <= th1 && th1 <= th1ub) {
      *func = (fv->x[1] - ypt1) * cos(theta1) - (fv->x[0] - xpt1) * sin(theta1);
      d_func[MESH_DISPLACEMENT1] = -sin(theta1);
      d_func[MESH_DISPLACEMENT2] = cos(theta1);
    } else {
      fprintf(stderr, "Double Rad case not found... %g %g %g %g %g\n", fv->x[0], fv->x[1], th1, th2,
              th2t);
    }
  }

  if (ielem_dim == 3)
    d_func[MESH_DISPLACEMENT3] = 0.0;

} /* END of routine f_double_fillet */

/*****************************************************************************/
/* This function is used to create a double fillet geometry boundary condition.
 * It is used in the case of a die with two fillets on the edges.
 *
 * This is a more geometry-based approach instead of angle based like f_double_rad and
 * f_double_fillet
 */
void f_double_fillet_geom_based(const int ielem_dim,
                                double *func,
                                double d_func[],     /* dimensioned [MAX_VARIABLE_TYPES+MAX_CONC] */
                                const double *p,     /*  function parameters from data card  */
                                const int num_const) /* number of passed parameters   */
{
  /**************************** EXECUTION BEGINS *******************************/

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  if (num_const < 8)
    GOMA_EH(GOMA_ERROR, "Need at least 8 parameters for Double Rad lip geometry bc!\n");
  dbl pt1_x = p[0];
  dbl pt1_y = p[1];
  dbl theta_1 = p[2];
  dbl radius_1 = p[3];
  dbl pt2_x = p[4];
  dbl pt2_y = p[5];
  dbl theta_2 = p[6];
  dbl radius_2 = p[7];
  dbl curv_mid = 0.0;
  if (num_const >= 8) {
    curv_mid = p[8];
  }

  bool is_curved = DOUBLE_NONZERO(curv_mid);
  if (is_curved) {
    GOMA_EH(GOMA_ERROR, "DOUBLE_FILLET_GEOM_BASED does not support curved centers yet.\n");
  }

  dbl xdir = pt1_x - pt2_x;
  dbl ydir = pt1_y - pt2_y;
  dbl sangle = atan2(ydir, xdir);

  // fillet 1
  dbl fillet_cen[DIM];
  dbl start_main[DIM];
  dbl start_edge[DIM];
  dbl pt1[DIM] = {pt1_x, pt1_y};
  dbl pt2[DIM] = {pt2_x, pt2_y};
  dbl rot1[DIM];
  rotate_line(pt1, pt2, theta_1, rot1);
  fillet_center(pt1, pt2, rot1, radius_1, fillet_cen, start_main, start_edge);
  dbl xcen1 = fillet_cen[0];
  dbl ycen1 = fillet_cen[1];

  // fillet 2
  dbl fillet_cen2[DIM];
  dbl rot2[DIM];
  dbl end_main[DIM];
  dbl end_edge[DIM];
  rotate_line(pt2, pt1, theta_2, rot2);
  fillet_center(pt2, pt1, rot2, radius_2, fillet_cen2, end_main, end_edge);
  dbl xcen2 = fillet_cen2[0];
  dbl ycen2 = fillet_cen2[1];

#if 0
  static int once = 1;
  if (once) {

    printf("Double Fillet BC: pt1 (%g, %g), pt2 (%g, %g)\n", pt1_x, pt1_y, pt2_x, pt2_y);
    printf("Double Fillet BC: xcen1 (%g, %g), xcen2 (%g, %g)\n", xcen1, ycen1, xcen2, ycen2);
    printf("Double Fillet BC: radius_1 %g, radius_2 %g\n", radius_1, radius_2);
    printf("Double Fillet BC: start_main (%g, %g), start_edge (%g, %g)\n", start_main[0],
           start_main[1], start_edge[0], start_edge[1]);
    printf("Double Fillet BC: end_main (%g, %g), end_edge (%g, %g)\n", end_main[0], end_main[1],
           end_edge[0], end_edge[1]);
    once = 0;
  }
#endif

  // 5 regions
  dbl fdist[5];
  dbl fdx[5];
  dbl fdy[5];

  // fillet 1
  fdist[0] = radius_1 - sqrt(SQUARE(fv->x[0] - xcen1) + SQUARE(fv->x[1] - ycen1));
  fdx[0] = -(fv->x[0] - xcen1) / sqrt(SQUARE(fv->x[0] - xcen1) + SQUARE(fv->x[1] - ycen1));
  fdy[0] = -(fv->x[1] - ycen1) / sqrt(SQUARE(fv->x[0] - xcen1) + SQUARE(fv->x[1] - ycen1));

  // fillet 2
  fdist[1] = radius_2 - sqrt(SQUARE(fv->x[0] - xcen2) + SQUARE(fv->x[1] - ycen2));
  fdx[1] = -(fv->x[0] - xcen2) / sqrt(SQUARE(fv->x[0] - xcen2) + SQUARE(fv->x[1] - ycen2));
  fdy[1] = -(fv->x[1] - ycen2) / sqrt(SQUARE(fv->x[0] - xcen2) + SQUARE(fv->x[1] - ycen2));

  // main slope
  // check if point is on the main slope
  fdist[2] = (fv->x[1] - pt1_y) * cos(sangle) + (fv->x[0] - pt1_x) * sin(sangle);
  fdx[2] = sin(sangle);
  fdy[2] = cos(sangle);

  // edge 1
  dbl slope_edge = atan2(start_edge[1] - pt1_y, start_edge[0] - pt1_x);
  fdist[3] = ((fv->x[1] - start_edge[1]) * cos(slope_edge)) +
             ((fv->x[0] - start_edge[0]) * sin(slope_edge));
  fdx[3] = sin(slope_edge);
  fdy[3] = cos(slope_edge);

  // edge 2
  dbl slope_edge2 = atan2(end_edge[1] - pt2_y, end_edge[0] - pt2_x);
  fdist[4] =
      ((fv->x[1] - end_edge[1]) * cos(slope_edge2)) + ((fv->x[0] - end_edge[0]) * sin(slope_edge2));
  fdx[4] = sin(slope_edge2);
  fdy[4] = cos(slope_edge2);

  dbl fil1[DIM] = {fv->x[0] - xcen1, fv->x[1] - ycen1};
  dbl edge1[DIM] = {start_edge[0] - xcen1, start_edge[1] - ycen1};
  dbl main1[DIM] = {start_main[0] - xcen1, start_main[1] - ycen1};

  dbl fil2[DIM] = {fv->x[0] - xcen2, fv->x[1] - ycen2};
  dbl edge2[DIM] = {end_edge[0] - xcen2, end_edge[1] - ycen2};
  dbl main2[DIM] = {end_main[0] - xcen2, end_main[1] - ycen2};

  int min_idx = 0;
  if (in_between_vectors_2D(edge1, main1, fil1) &&
      near_circle(fillet_cen, radius_1, fv->x, 1.5 * radius_1)) {
    min_idx = 0;
  } else if (in_between_vectors_2D(edge2, main2, fil2) &&
             near_circle(fillet_cen2, radius_2, fv->x, 1.5 * radius_2)) {
    min_idx = 1;
  } else {
    for (int i = 2; i < 5; i++) {
      if (fabs(fdist[i]) < fabs(fdist[min_idx])) {
        min_idx = i;
      }
    }
  }

  *func = fdist[min_idx];

  d_func[MESH_DISPLACEMENT1] = fdx[min_idx];
  d_func[MESH_DISPLACEMENT2] = fdy[min_idx];

  if (ielem_dim == 3)
    d_func[MESH_DISPLACEMENT3] = 0.0;

} /* END of routine f_double_fillet_geom_based */

void f_roll_fluid(int ielem_dim,
                  double *func,
                  double d_func[],     /* dimensioned [MAX_VARIABLE_TYPES+MAX_CONC] */
                  const double *p,     /*  function parameters from data card  */
                  const int num_const, /* number of passed parameters   */
                  double *xsurf)       /* number of passed parameters   */
{
  /**************************** EXECUTION BEGINS *******************************/
  double roll_rad;     /* roll radius */
  double origin[3];    /* roll axis origin (x,y,z) */
  double dir_angle[3]; /* axis direction angles */
  double coord[3];     /* current coordinates */
  double axis_pt[3], rad_dir[3], d_dist[3], dist, R, factor, t;
  double omega, v_dir[3], v_roll[3];
  double velo_avg = 0.0, pgrad = 0.;
  double v_solid = 0., res, jac, delta, flow, eps = 1.0e-8, viscinv;
  double jacinv, thick;
  int Pflag = TRUE;
  double pg_factor = 1.0, tang_sgn = 1.0, v_mag = 0.;
  ;

  int j, var;
#if 0
  int jvar,k;
  double dthick_dV, dthick_dP;
#endif

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  if (num_const < 7)
    GOMA_EH(GOMA_ERROR, "Need at least 7 parameters for Roll geometry bc!\n");

  roll_rad = p[0];
  origin[0] = p[1];
  origin[1] = p[2];
  origin[2] = p[3];
  dir_angle[0] = p[4];
  dir_angle[1] = p[5];
  dir_angle[2] = p[6];

  /* calculate distance from interface surface to solid surface for repulsion calculations */

  coord[0] = fv->x[0];
  coord[1] = fv->x[1];
  if (ielem_dim == 3) {
    coord[2] = fv->x[2];
  } else {
    coord[2] = 0.0;
  }

  /*  find intersection of axis with normal plane - i.e., locate point on
   *          axis that intersects plane normal to axis that contains local point. */

  factor = SQUARE(dir_angle[0]) + SQUARE(dir_angle[1]) + SQUARE(dir_angle[2]);
  t = (dir_angle[0] * (coord[0] - origin[0]) + dir_angle[1] * (coord[1] - origin[1]) +
       dir_angle[2] * (coord[2] - origin[2])) /
      factor;
  axis_pt[0] = origin[0] + dir_angle[0] * t;
  axis_pt[1] = origin[1] + dir_angle[1] * t;
  axis_pt[2] = origin[2] + dir_angle[2] * t;

  /*  compute radius and radial direction */

  R = sqrt(SQUARE(coord[0] - axis_pt[0]) + SQUARE(coord[1] - axis_pt[1]) +
           SQUARE(coord[2] - axis_pt[2]));
  rad_dir[0] = (coord[0] - axis_pt[0]) / R;
  rad_dir[1] = (coord[1] - axis_pt[1]) / R;
  rad_dir[2] = (coord[2] - axis_pt[2]) / R;
  dist = R - roll_rad;
  d_dist[0] = rad_dir[0] * (1. - SQUARE(dir_angle[0]) / factor) +
              rad_dir[1] * (-dir_angle[1] * dir_angle[0] / factor) +
              rad_dir[2] * (-dir_angle[2] * dir_angle[0] / factor);
  d_dist[1] = rad_dir[1] * (1. - SQUARE(dir_angle[1]) / factor) +
              rad_dir[0] * (-dir_angle[0] * dir_angle[1] / factor) +
              rad_dir[2] * (-dir_angle[2] * dir_angle[1] / factor);
  d_dist[2] = rad_dir[2] * (1. - SQUARE(dir_angle[2]) / factor) +
              rad_dir[0] * (-dir_angle[0] * dir_angle[2] / factor) +
              rad_dir[1] * (-dir_angle[1] * dir_angle[2] / factor);

  if (num_const < 10)
    GOMA_WH(-1, "ROLL_FLUID: Less than 10 parameters - reverting to roll surface!\n");

  omega = p[7];
  /* compute velocity direction as perpendicular to both axis and radial
   *         direction.  Positive direction is determined by right hand rule */

  v_dir[0] = dir_angle[1] * rad_dir[2] - dir_angle[2] * rad_dir[1];
  v_dir[1] = dir_angle[2] * rad_dir[0] - dir_angle[0] * rad_dir[2];
  v_dir[2] = dir_angle[0] * rad_dir[1] - dir_angle[1] * rad_dir[0];

  v_roll[0] = omega * roll_rad * v_dir[0];
  v_roll[1] = omega * roll_rad * v_dir[1];
  v_roll[2] = omega * roll_rad * v_dir[2];

  if (TimeIntegration == TRANSIENT && pd->gv[R_MESH1]) {
    /* Add the mesh motion to the substrate velocity */
    v_roll[0] += fv_dot->x[0];
    v_roll[1] += fv_dot->x[1];
    v_roll[2] += fv_dot->x[2];
  }

  /* quantities specific to FLUID bcs   */

  if (num_const > 8 && p[9] >= 0.0) {
    dist = 0.;
    for (var = 0; var < pd->Num_Dim; var++) {
      /* Uses undeformed node position */
      dist += SQUARE(fv->x0[var] - xsurf[var]);
    }
    dist /= SQUARE(p[10]);
    /*if(dist < 10)fprintf(stderr,"roll_fl %g %g %g\n",fv->x0[0],xsurf[0],dist);
     */

    Pflag = (int)p[11];
    velo_avg = 0.0;
    pgrad = 0.;
    v_mag = 0.;
    for (j = 0; j < pd->Num_Dim; j++) {
      velo_avg += fv->stangent[0][j] * (v_roll[j] + fv->v[j]);
      v_solid += fv->stangent[0][j] * v_roll[j];
      v_mag += SQUARE(v_roll[j]);
      if (Pflag) {
        pgrad += fv->stangent[0][j] * fv->grad_P[j];
      }
    }
    v_mag = sqrt(v_mag);
    tang_sgn = v_solid / v_mag;
    tang_sgn = (double)SGN(v_solid / v_mag);
    velo_avg *= 0.5;
    /* sometimes the tangent/normals flip causing havoc....*/
    if (v_solid < 0) {
      GOMA_WH(-1, "fvelo_slip: normals and tangents have flipped! - try CONTACT_LINE model\n");
      velo_avg *= tang_sgn;
      v_solid *= tang_sgn;
      pgrad *= tang_sgn;
    }

    pg_factor = 1.0;
    if (dist < 10.0) {
      pg_factor = 1.0 - exp(-dist);
    }
    pgrad *= pg_factor;

    flow = MAX(0., p[9] * v_solid);
    viscinv = 1. / p[8];
    thick = flow / velo_avg;
    j = 0;
    do {
      res = -CUBE(thick) * viscinv * pgrad / 12. + thick * velo_avg - flow;
      jac = -0.25 * SQUARE(thick) * viscinv * pgrad + velo_avg;
      jacinv = 1.0 / jac;
      delta = -res * jacinv;
      thick += delta;
      j++;
    } while (fabs(delta) > eps && j < 20);
#if 0
      dthick_dV = -0.5*jacinv;     /*  1/h*derivative  */
      dthick_dP = CUBE(thick)*viscinv/12.*jacinv;
#endif
#if 0
fprintf(stderr,"slip %d %g %g %g %g\n",Pflag,fv->x[0],thick,flow/v_solid,velo_avg);
fprintf(stderr,"more %g %g %g %g\n",res,jac, dthick_dV,dthick_dP);
#endif
    thick = 0.;
    *func = dist - thick;
    d_func[MESH_DISPLACEMENT1] = d_dist[0];
    d_func[MESH_DISPLACEMENT2] = d_dist[1];
    d_func[MESH_DISPLACEMENT3] = d_dist[2];
#if 0
    for (jvar=0; jvar<pd->Num_Dim; jvar++)
      {
        var = VELOCITY1 + jvar;
        for (k=0; k<pd->Num_Dim; k++)
          {
           d_func[var] += -thick*dthick_dV*fv->stangent[0][k];
          }
       }
#endif
#if 0
/* Mesh motion Jacobian entries   */
        for (jvar=0; jvar<ei->ielem_dim; jvar++)
          {
            var = MESH_DISPLACEMENT1 + jvar;
            if (pd->v[pg->imtrx][var])
              {
                    for (k = 0; k < pd->Num_Dim; k++)
                      {
                        d_func[var] += -thick*dthick_dV*fv->v[k]
                                *fv->stangent[0][k];
                        if(Pflag)
                          {
                          d_func[var] += -dthick_dP*pg_factor*fv->grad_P[k]*fv->stangent[0][k];
                          }
                      }
              }
          }

#endif
#if 0
   var = PRESSURE;
    if (pd->v[pg->imtrx][var])
      {
        if(Pflag )
          {
               for (k = 0; k < pd->Num_Dim; k++)
                  {
                    d_func[var] += -dthick_dP*pg_factor*fv->stangent[0][k];
                  }
          }
      }
#endif
  } else {
    *func = dist;
    d_func[MESH_DISPLACEMENT1] = d_dist[0];
    d_func[MESH_DISPLACEMENT2] = d_dist[1];
    d_func[MESH_DISPLACEMENT3] = d_dist[2];
  }

} /* END of routine f_roll_fluid                                                   */

#ifdef FEATURE_ROLLON_PLEASE
#include "feature_rollon.h"
#endif

/*****************************************************************************/
/*       Functions for user defined solid boundary description               */
/*****************************************************************************/

void fspline(int ielem_dim,
             double *func,
             double d_func[], /* dimensioned [MAX_VARIABLE_TYPES + MAX_CONC] */
             double p[],      /* parameters to parameterize temperature eqn model*/
             double time)     /* time  at which bc's are evaluated     */
{
  if (af->Assemble_LSA_Mass_Matrix)
    return;

  d_func[MESH_DISPLACEMENT1] = dfncd1(fv->x[0], fv->x[1], fv->x[2], p, time);

  d_func[MESH_DISPLACEMENT2] = dfncd2(fv->x[0], fv->x[1], fv->x[2], p, time);

  if (ielem_dim == 3)
    d_func[MESH_DISPLACEMENT3] = dfncd3(fv->x[0], fv->x[1], fv->x[2], p, time);

  *func = fnc(fv->x[0], fv->x[1], fv->x[2], p, time);

} /* END of routine fspline                                                  */

void fspline_rs(int ielem_dim,
                double *func,
                double d_func[], /* dimensioned [MAX_VARIABLE_TYPES + MAX_CONC] */
                double p[],      /* parameters to parameterize temperature eqn model*/
                double time)     /* time  at which bc's are evaluated     */
{
  if (af->Assemble_LSA_Mass_Matrix)
    return;

  d_func[SOLID_DISPLACEMENT1] = dfncd1(fv->x[0], fv->x[1], fv->x[2], p, time);

  d_func[SOLID_DISPLACEMENT2] = dfncd2(fv->x[0], fv->x[1], fv->x[2], p, time);

  if (ielem_dim == 3)
    d_func[SOLID_DISPLACEMENT3] = dfncd3(fv->x[0], fv->x[1], fv->x[2], p, time);

  *func = fnc(fv->x[0], fv->x[1], fv->x[2], p, time);

} /* END of routine fspline_rs */