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
#ifndef lint
static char *cvs_contbord_id =
  "$Id: loca_bord.c,v 5.1 2007-09-18 18:53:42 prschun Exp $";
#endif
*/
/*
 -----------------------------------------------------------------------------
   LOCA 1.0: Library of Continuation Algorithms
   Copyright (C) 2001, Sandia National Laboratories
 -----------------------------------------------------------------------------
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "loca_const.h"
#include "loca_util_const.h"

/******************************************************************************
 *
 *           GLOBALS DEFINED IN THIS FILE
 *
 *       Function                type                       Called By
 *       --------------        ---------                ----------------------
 *
 *
 ******************************************************************************/

static double scalar_perturbation(double param, double perturb);
extern double scaled_dot_prod(double *, double *, double *, int);

int AGS_option = 0;

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int arc_length_bordering_alg(
    double *x, double *delta_x, struct con_struct *con, double reltol, double abstol)

/*
 * Bordering algorithm for arc-length continuation is performed.
 * Notation is from page 181 of JNS thesis. The vector "w" is
 * already calculated (eq 7-21b), and in array delta_x.
 */

{
  /* shorthand for long con sub-structure */
  struct general_info_struct *cgi = &(con->general_info);
  struct private_info_struct *cpi = &(con->private_info);

  double *dx2, y;
  int i;
  double arc_length_equation = 0.0, param_update, scaled_resid;
  double al_eq_soln_contrib, al_eq_param_contrib;

  /*************************** BEGIN EXECUTION *******************************/

  /*
   * Next, "u" is calculated and stored as dx2. (eq. 7-21a) This is exactly
   * the same linear solve as a continuation predictor step.
   */

  dx2 = alloc_vec();

  calc_rhs_continuation(ARC_CONT_SOL2, x, dx2, NULL, NULL, NULL, cgi->param, cgi->perturb, NULL,
                        cgi->numUnks, cgi->numOwnedUnks);
  i = linear_solver_conwrap(dx2, CHECK_JACOBIAN, NULL);

  /*
   * Calculate the arc-length equation ("g" in bordering algorithm notation)
   */

  al_eq_soln_contrib = scaled_dot_prod(x, cpi->x_tang, cpi->scale_vec, cgi->numOwnedUnks) -
                       scaled_dot_prod(cpi->x_old, cpi->x_tang, cpi->scale_vec, cgi->numOwnedUnks);

  al_eq_param_contrib = cpi->dp_ds * (cgi->param - cpi->param_old);

  arc_length_equation = al_eq_soln_contrib + al_eq_param_contrib - cpi->arc_step;

  /* Calculate "y", the -update to the continuation parameter (JNS eq. 7-22) */

  y = (arc_length_equation -
       scaled_dot_prod(cpi->x_tang, delta_x, cpi->scale_vec, cgi->numOwnedUnks)) /
      (cpi->dp_ds - scaled_dot_prod(cpi->x_tang, dx2, cpi->scale_vec, cgi->numOwnedUnks));

  /*
   * y is subtracted because the minus sign in Newton's method has been
   * switched to the update term -- the same reason delta_x is subtracted
   * in update_soln
   */

  cgi->param -= y;
  assign_parameter_conwrap(cgi->param);

  /* Modify delta_x, which will be returned to the nonlinear solver (eq.7-23)*/

  for (i = 0; i < cgi->numUnks; i++)
    delta_x[i] -= y * dx2[i];

  free_vec(&dx2);

  /*
   * Check whether or not the arc-length equation and continuation param
   * are converged
   */

  param_update = fabs(y) / (reltol * fabs(cgi->param) + abstol);
  scaled_resid = fabs(arc_length_equation) / (reltol * fabs(cpi->arc_step) + abstol);

  if (cgi->printproc > 7) {
    printf("\n\t\tFraction of arc-step due to soln, param = %g %g\n",
           al_eq_soln_contrib / cpi->arc_step, al_eq_param_contrib / cpi->arc_step);
  }
  if (cgi->printproc > 4) {
    printf("\n\tArc Length Continuation: Convergence Criteria\n");
    printf("\tVariable     Scaled Update (<1)  Unscaled Update  New Value\n");
    printf("\t***********************************************************\n");
    printf("\tparameter    %e        %e    %g\n", param_update, -y, cgi->param);
    printf("\tarc length   %e        %e    %g\n", scaled_resid, arc_length_equation,
           cpi->arc_step + arc_length_equation);
    printf("\t***********************************************************\n");
  }

  if ((param_update < 1.0) && (scaled_resid < 1.0))
    return (TRUE);
  else
    return (FALSE);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int turning_point_alg(
    double *x, double *delta_x, struct con_struct *con, double reltol, double abstol)

/*
 * Algorithm for locking on to a turning point.
 * Theory currently comes from a TEX document of Louis Romero. (AGS 1/98)
 * Lines labeled  SCALED  are additions for new scaling (section 4 of same
 * write-up). The SCALED lines can be commented out to recover unscaled version.
 */
#define SCALE_TP
{
  /* shorthand for long con sub-structure */
  struct general_info_struct *cgi = &(con->general_info);
  struct private_info_struct *cpi = &(con->private_info);

  double *a, *b, *c, *d, *x_tmp, dt_p;
#ifdef SCALE_TP
  double a_big, b_big, c_big, d_big;
#endif
  int i;
  double param_update, r_update, vecnorm, gnum_unks, tmp, RayQ;
  double *y = cpi->x_tang;      /* null vector */
  double *phi = cpi->scale_vec; /* vector for normalizing the null vector */
  static int first = TRUE;

  /* This flag can be 0 or 1 and effects scaling */
  /* Seems to work better as 0 for TP calcs      */
  AGS_option = 0;

  /*************************** BEGIN EXECUTION *******************************/

  /* Allocate arrays for turning point bordering algorithm */

  a = alloc_vec();
  b = alloc_vec();
  c = alloc_vec();
  d = alloc_vec();
  x_tmp = alloc_vec();

  /* construct "a" vector from delta_x */

  for (i = 0; i < cgi->numUnks; i++)
    a[i] = -delta_x[i];

  /*
   * Next, "b" is calculated. This is exactly
   * the same linear solve as a continuation predictor step except
   * the bif_param is perturbed, not the continuation param.
   */

  calc_rhs_continuation(TP_CONT_SOL2, x, b, NULL, NULL, NULL, con->turning_point_info.bif_param,
                        cgi->perturb, NULL, cgi->numUnks, cgi->numOwnedUnks);

  i = linear_solver_conwrap(b, CHECK_JACOBIAN, NULL);

  /*
   * First time through this routine, initialize null vector y and
   * norm vector phi with scaled b.  First Newton iter of subsequent
   * continuation steps, set guess for y and phi equal to converged
   * null vector of previous converged step.
   * If using a previous null vector, it has already been read in.
   */

  if (first) {
    /* Commented out by EDW 1/28/2002:
        vecnorm = ltransnorm(b, b);

        if (vecnorm <= 0.0) {
          if (cgi->printproc > 4)
            printf("\n\t### Error tp alg: vecnorm = %g\n\n", vecnorm);
          exit(-1);
        }

        for (i=0; i < cgi->numUnks; i++) phi[i] = b[i] / sqrt(vecnorm);
    */

    if (cgi->nv_restart == FALSE) {

      calc_scale_vec_conwrap(x, phi, cgi->numUnks);

      for (i = 0; i < cgi->numUnks; i++)
        y[i] = phi[i];

      vecnorm = ltransnorm(b, phi);

      for (i = 0; i < cgi->numUnks; i++)
        y[i] = b[i] / vecnorm;

    }

    /* Load in old null vector if specified */
    else {

      if (cgi->printproc > 7)
        printf("Loading in previous null vector ...\n");

      for (i = 0; i < cgi->numUnks; i++)
        y[i] = con->turning_point_info.nv[i];

      calc_scale_vec_conwrap(x, phi, cgi->numUnks);

      vecnorm = ltransnorm(y, phi);

      for (i = 0; i < cgi->numUnks; i++)
        y[i] /= vecnorm;
    }

    first = FALSE;
  }

  else {

    /* This section is optional, changing scale vector between Newton iters */

    calc_scale_vec_conwrap(x, phi, cgi->numUnks);

    vecnorm = ltransnorm(y, phi);
    if (cgi->printproc > 7)
      printf("\t\tRescaling r_vec by %g to make it length 1\n", vecnorm);

    for (i = 0; i < cgi->numUnks; i++)
      y[i] /= vecnorm;
  }

  /* Commented out by EDW 1/28/2002:
    else if (cpi->first_iter) {

      vecnorm = ltransnorm(y, y);

      for (i=0; i < cgi->numUnks; i++) phi[i] = y[i] / sqrt(vecnorm);
      for (i=0; i < cgi->numUnks; i++) y[i] = phi[i];
    }
  */

  /* Rescale a and b vectors as in Louie's write-up, section 4 */

#ifdef SCALE_TP
  b_big = 1.0 / ltransnorm(b, phi);    /*SCALED*/
  a_big = -ltransnorm(a, phi) * b_big; /*SCALED*/
  for (i = 0; i < cgi->numUnks; i++)
    a[i] += a_big * b[i];              /*SCALED*/
  for (i = 0; i < cgi->numUnks; i++)
    b[i] *= b_big;                     /*SCALED*/
#endif

  /* Next, "c" is calculated as a function of a and y. */

  calc_rhs_continuation(TP_CONT_SOL3, x, c, a, phi, x_tmp, con->turning_point_info.bif_param,
                        cgi->perturb, y, cgi->numUnks, cgi->numOwnedUnks);

  /* Get null vector residual now. */

  RayQ = null_vector_resid(0.0, 0.0, y, NULL, FALSE);

  i = linear_solver_conwrap(c, SAME_BUT_UNSCALED_JACOBIAN, x_tmp);

  /* Next, "d" is calculated as a function of b and y. */

  calc_rhs_continuation(TP_CONT_SOL4, x, d, b, phi, x_tmp, con->turning_point_info.bif_param,
                        cgi->perturb, y, cgi->numUnks, cgi->numOwnedUnks);

  i = linear_solver_conwrap(d, SAME_BUT_UNSCALED_JACOBIAN, x_tmp);

  /*
   * Calculate the updates to bif_param (stored in dt_p),
   * y (stored in c), and x (stored in -delta_x).
   */

  /* Rescale c and d vectors as in Louie's write-up, section 4 */

#ifdef SCALE_TP
  d_big = 1.0 / ltransnorm(d, phi);    /*SCALED*/
  c_big = -ltransnorm(c, phi) * d_big; /*SCALED*/
  for (i = 0; i < cgi->numUnks; i++)
    c[i] += c_big * d[i];              /*SCALED*/
  for (i = 0; i < cgi->numUnks; i++)
    d[i] *= d_big;                     /*SCALED*/
#endif

  dt_p = ((1.0 - AGS_option * ltransnorm(y, phi)) - ltransnorm(c, phi)) / ltransnorm(d, phi);

  for (i = 0; i < cgi->numUnks; i++)
    c[i] += dt_p * d[i] + (AGS_option - 1) * y[i];

    /* dt_p change meaning from Beta to alpha here */

#ifdef SCALE_TP
  dt_p = dt_p * d_big + c_big; /*SCALED*/
#endif

  for (i = 0; i < cgi->numUnks; i++)
    delta_x[i] = -a[i];

  for (i = 0; i < cgi->numUnks; i++)
    delta_x[i] -= dt_p * b[i];

    /*dt_p change meaning from alpha to dt_p here */

#ifdef SCALE_TP
  dt_p = a_big + dt_p * b_big; /*SCALED*/
#endif

  /* If delta_r is not nearly zero, truncate x and r steps for safety */

  /*
    if (fabs(ltransnorm(c, phi)) > 0.00001)  {
      vecnorm = 0.00001/fabs(ltransnorm(c, phi));
      if (cgi->printproc > 7)
        printf("Roundoff causing nonzero delta_r, rescaling updates by %g\n",
                vecnorm);
      for (i=0; i < cgi->numUnks; i++) {
        delta_x[i] *= vecnorm;
        c[i] *= vecnorm;
      }
    }
  */

  /*
   * Update the bif_param and the null vector y
   */

  con->turning_point_info.bif_param += dt_p;
  assign_bif_parameter_conwrap(con->turning_point_info.bif_param);

  for (i = 0; i < cgi->numUnks; i++)
    y[i] += c[i];

  /*
   * Check whether or not the arc-length equation and continuation param
   * are converged
   */

  param_update = fabs(dt_p) / (reltol * fabs(con->turning_point_info.bif_param) + abstol);

  /* get update norm of y, to check for convergence */

  r_update = 0.0;
  for (i = 0; i < cgi->numUnks; i++) {
    tmp = c[i] / (fabs(y[i]) * reltol + abstol);
    r_update += tmp * tmp;
  }
  gnum_unks = gsum_double_conwrap(1.0 * cgi->numOwnedUnks);
  r_update = sqrt(gsum_double_conwrap(r_update) / gnum_unks);

  if (RayQ == -1.0) {
    if (cgi->printproc > 1)
      printf("\n\tRayleigh Quotient Error: zero denom\n");
  } else {
    if (cgi->printproc > 4)
      printf("\n\tRayleigh Quotient before updates: %g\n", RayQ);
  }

  if (cgi->printproc > 4) {
    printf("\n\tTurning Point Continuation: Convergence Criteria\n");
    printf("\tVariable     Scaled Update (<1)  Unscaled Update  New Value\n");
    printf("\t***********************************************************\n");
    printf("\tparameter    %e        %e    %g\n", param_update, dt_p,
           con->turning_point_info.bif_param);
    printf("\tNull vector  %e\n", r_update);
    printf("\t***********************************************************\n");
  }

  free_vec(&a);
  free_vec(&b);
  free_vec(&c);
  free_vec(&d);
  free_vec(&x_tmp);

  /* return convergence status of the parameter and null vector */

  if ((param_update < 1.0) && (r_update < 10.0))
    return (TRUE);
  else
    return (FALSE);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
static int
pitchfork_alg(double *x, double *delta_x, struct con_struct *con, double reltol, double abstol)

/*
 * Algorithm for locking on to a Pitchfork bifurcation
 * This assumes that cpi->x_tang contains an asymmetric vector
 * the first time through, and an initial guess for the
 * null vector every time (with both meanings the first time).
 */

{
  /* shorthand for long con sub-structure */
  struct general_info_struct *cgi = &(con->general_info);
  struct private_info_struct *cpi = &(con->private_info);

  double *a, *b, *c, *d, *e, *f, *x_tmp, dt_p;
  int i;
  double param_update, r_update, vecnorm, gnum_unks, tmp, RayQ;
  double *phi = cpi->x_tang;
  double eps_update, d_eps, ipa, ipb, ipc, ipx, ltd, lte, ltf;
  double t1, t2, t3, t4;
  double *psi = con->pitchfork_info.psi;
  static int first = TRUE;
  static double eps;

  /* This flag can be 0 or 1 and effects scaling */
  /* Seems to work better as 1 for PF calcs      */
  AGS_option = 1;

  /*************************** BEGIN EXECUTION *******************************/

  /* first time through this routine, set psi asymmetric vector to con-tang
     and set epsilon "asymmetry variable" to 0.0 */

  if (first) {
    if (cgi->printproc > 7)
      printf("\n\tIn pitchfork alg, AGS_option=%d\n", AGS_option);
    eps = 0.0;
    tmp = sqrt(ip(psi, psi));
    if (tmp == 0.0) {
      if (cgi->printproc > 1) {
        printf("ERROR in pitchfork alg: Asymmetric vector must be supplied\n");
      }
      exit(-1);
    }
    for (i = 0; i < cgi->numUnks; i++)
      psi[i] = psi[i] / tmp;

    /* psi is used as phi vector as well */

    vec_copy(psi, phi);

    first = FALSE;

    t1 = ip(x, psi);
    t2 = ip(x, x);
    t3 = ip(psi, psi);
    t4 = t1 / sqrt(t2 * t3);
    if (cgi->printproc > 7) {
      printf("\tPitchfork Alg: On input <x,psi> = %g; Should be 0.0.\n", t4);
    }
  }

  /* calculate variable averages, to be used in ltransnorm calls below */

  calc_scale_vec_conwrap(x, cpi->scale_vec, cgi->numOwnedUnks);

  /* Make sure the null vector phi is length 1 to begin with */

  vecnorm = ltransnorm(phi, cpi->scale_vec);
  if (cgi->printproc > 7)
    printf("\n\t\tRescaling phi by %g to make it length 1\n", vecnorm);

  for (i = 0; i < cgi->numUnks; i++)
    phi[i] /= vecnorm;

  /* Allocate arrays for pitchfork bordering algorithm */

  a = alloc_vec();
  b = alloc_vec();
  c = alloc_vec();
  d = alloc_vec();
  e = alloc_vec();
  f = alloc_vec();
  x_tmp = alloc_vec();

  /* Begin real stuff */
  /* construct "a" vector from delta_x */

  for (i = 0; i < cgi->numUnks; i++)
    a[i] = -delta_x[i];

  /*
   * Next, "b" is calculated. This is exactly
   * the same linear solve as a continuation predictor step except
   * the bif_param is perturbed, not the continuation param.
   */

  calc_rhs_continuation(TP_CONT_SOL2, x, b, NULL, NULL, NULL, con->pitchfork_info.bif_param,
                        cgi->perturb, NULL, cgi->numUnks, cgi->numOwnedUnks);

  i = linear_solver_conwrap(b, CHECK_JACOBIAN, NULL);

  /* Next, "c" is calculated using just psi as rhs */

  for (i = 0; i < cgi->numUnks; i++)
    c[i] = -psi[i];

  i = linear_solver_conwrap(c, OLD_JACOBIAN, NULL);

  /* Next, "d" is calculated as a function of a and phi. */

  calc_rhs_continuation(TP_CONT_SOL3, x, d, a, cpi->scale_vec, x_tmp, con->pitchfork_info.bif_param,
                        cgi->perturb, phi, cgi->numUnks, cgi->numOwnedUnks);

  /* Get null vector residual now. */

  RayQ = null_vector_resid(0.0, 0.0, phi, NULL, FALSE);

  i = linear_solver_conwrap(d, SAME_BUT_UNSCALED_JACOBIAN, x_tmp);

  if (AGS_option == 1)
    for (i = 0; i < cgi->numUnks; i++)
      d[i] += phi[i];

  /* Next, "e" is calculated as a function of b and phi. */

  calc_rhs_continuation(TP_CONT_SOL4, x, e, b, cpi->scale_vec, x_tmp, con->pitchfork_info.bif_param,
                        cgi->perturb, phi, cgi->numUnks, cgi->numOwnedUnks);

  i = linear_solver_conwrap(e, SAME_BUT_UNSCALED_JACOBIAN, x_tmp);

  /* Next, "f" is calculated as a function of c and phi. */

  calc_rhs_continuation(TP_CONT_SOL3, x, f, c, cpi->scale_vec, x_tmp, con->pitchfork_info.bif_param,
                        cgi->perturb, phi, cgi->numUnks, cgi->numOwnedUnks);

  i = linear_solver_conwrap(f, SAME_BUT_UNSCALED_JACOBIAN, x_tmp);

  if (AGS_option == 1)
    for (i = 0; i < cgi->numUnks; i++)
      f[i] += phi[i];

  /*
   * Calculate the updates to bif_param (stored in dt_p),
   * phi (stored in d), and x (stored in -delta_x).
   */

  ipx = ip(x, psi);
  ipa = ip(a, psi);
  ipb = ip(b, psi);
  ipc = ip(c, psi);
  ltd = ltransnorm(d, cpi->scale_vec);
  lte = ltransnorm(e, cpi->scale_vec);
  ltf = ltransnorm(f, cpi->scale_vec);

  d_eps = -eps + ((ipx + ipa) * lte + ipb * (1.0 - ltd)) / (ipb * ltf - ipc * lte);

  dt_p = (1.0 - ltd - ltf * (d_eps + eps)) / lte;

  /* Negative of update vector here */
  for (i = 0; i < cgi->numUnks; i++)
    delta_x[i] = -a[i] - c[i] * (eps + d_eps) - b[i] * dt_p;

  /* use c space for delta phi */

  for (i = 0; i < cgi->numUnks; i++)
    c[i] = -phi[i] + d[i] + e[i] * dt_p + f[i] * (eps + d_eps);

  /*
   * Update the bif_param, eps,  and the null vector phi
   */

  /*
  if (cgi->printproc > 7) printf("writing x, n, psi, delta_x, delta_n\n");
    solution_output_conwrap(x,       1, 0, 5, 0);
    solution_output_conwrap(phi,   1, 0, 6, 0);
    solution_output_conwrap(psi,     1, 0, 7, 0);
    solution_output_conwrap(delta_x, 1, 0, 8, 0);
    solution_output_conwrap(c,       1, 0, 9, 0);
  */

  eps += d_eps;

  con->pitchfork_info.bif_param += dt_p;
  assign_bif_parameter_conwrap(con->pitchfork_info.bif_param);

  for (i = 0; i < cgi->numUnks; i++)
    phi[i] += c[i];

  /*
   * Check whether or not the continuation param, null vector, and eps
   * are converged
   */

  param_update = fabs(dt_p) / (reltol * fabs(con->pitchfork_info.bif_param) + abstol);
  eps_update = fabs(d_eps) / (reltol * fabs(eps) + abstol);

  /* get update norm of phi, to check for convergence */

  r_update = 0.0;
  for (i = 0; i < cgi->numUnks; i++) {
    tmp = c[i] / (fabs(phi[i]) * reltol + abstol);
    r_update += tmp * tmp;
  }
  gnum_unks = gsum_double_conwrap(1.0 * cgi->numOwnedUnks);
  r_update = sqrt(gsum_double_conwrap(r_update) / gnum_unks);

  if (RayQ == -1.0) {
    if (cgi->printproc > 1)
      printf("\n\tRayleigh Quotient Error: zero denom\n");
  } else {
    if (cgi->printproc > 4)
      printf("\n\tRayleigh Quotient before updates: %g\n", RayQ);
  }

  if (cgi->printproc > 4) {
    printf("\n\tPitchfork Continuation: Convergence Criteria\n");
    printf("\tVariable     Scaled Update (<1)  Unscaled Update  New Value\n");
    printf("\t***********************************************************\n");
    printf("\tparameter    %e        %e    %g\n", param_update, dt_p,
           con->pitchfork_info.bif_param);
    printf("\teps          %e        %e    %g\n", eps_update, d_eps, eps);
    printf("\tNull vector  %e\n", r_update);
    printf("\t***********************************************************\n");
  }

  free_vec(&a);
  free_vec(&b);
  free_vec(&c);
  free_vec(&d);
  free_vec(&e);
  free_vec(&f);
  free_vec(&x_tmp);

  /* return convergence status of the parameter, eps and null vector */

  if (param_update < 1.0 && eps_update < 1.0 && r_update < 10.0)
    return (TRUE);
  else
    return (FALSE);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
static int
hopf_alg(double *x, double *delta_x, struct con_struct *con, double reltol, double abstol)

/*
 * Algorithm for locking on to and tracking a Hopf point.
 */

{
  /* shorthand for long con sub-structure */
  struct general_info_struct *cgi = &(con->general_info);
  struct private_info_struct *cpi = &(con->private_info);
  struct hopf_info_struct *chi = &(con->hopf_info);

  double *a, *b, *c, *d, *e, *f, *g, *h, *x_tmp, *M_1, *M_2, dt_p;
  int i;
  double x_update, param_update, gnum_unks, tmp, RayQ, tmp2, alpha, beta;
  double omega_update, y_update, z_update, y_tmp, z_tmp;
  static int first = TRUE, first2 = TRUE;
  double test, count;

  double ltc, ltd, lte, ltf, ltg, lth, delta_omega, delta_y, delta_z;

  /*************************** BEGIN EXECUTION *******************************/

  /* Allocate arrays for hopf bordering algorithm */
  a = alloc_vec();
  b = alloc_vec();
  c = alloc_vec();
  d = alloc_vec();
  e = alloc_vec();
  f = alloc_vec();
  g = alloc_vec();
  h = alloc_vec();
  x_tmp = alloc_vec();

  /*
   * First time through this routine, check (1) that eigenvectors are
   * non-zero and (2) the functional dependence of the Mass Matrix.
   */

  if (first) {
    if (cgi->printproc > 4)
      printf("\n\tEntering Hopf Algorithm\n");
    tmp = sqrt(dp(chi->y_vec, chi->y_vec));
    tmp2 = sqrt(dp(chi->z_vec, chi->z_vec));
    if ((tmp == 0.0) || (tmp2 == 0.0)) {
      if (cgi->printproc > 1) {
        printf("ERROR in Hopf alg: complex vector pair must be supplied\n");
      }
      exit(-1);
    }

    /*
     * Test to see if the Mass Matrix is a function of the
     * parameter or of the solution vector. Set mass_param and mass_x
     * accordingly - 1=true, 0=false
     */
    if (chi->mass_flag == 1) {
      cpi->mass_param = 1;
      cpi->mass_x = 1;
    } else if (chi->mass_flag == 2) {
      cpi->mass_param = 0;
      cpi->mass_x = 1;
    } else if (chi->mass_flag == 3) {
      cpi->mass_param = 1;
      cpi->mass_x = 0;
    } else if (chi->mass_flag == 4) {
      cpi->mass_param = 0;
      cpi->mass_x = 0;
    } else if (chi->mass_flag == 5) {

      /* Heuristic determination */
      cpi->mass_x = 0;
      cpi->mass_param = 0;

      if (cgi->printproc > 7)
        printf("\n\tHopf Continuation: Determining "
               "Mass Matrix Dependencies\n");
      mass_matrix_fill_conwrap(x, x_tmp);
      mass_matvec_mult_conwrap(chi->y_vec, a);
      matvec_mult_conwrap(chi->z_vec, c);

      for (i = 0; i < cgi->numUnks; i++)
        x_tmp[i] = x[i] + scalar_perturbation(x[i], cgi->perturb);
      mass_matrix_fill_conwrap(x_tmp, f);
      mass_matvec_mult_conwrap(chi->y_vec, b);
      matvec_mult_conwrap(chi->z_vec, d);
      count = 0.0;
      for (i = 0; i < cgi->numOwnedUnks; i++) {
        if (d[i] - c[i] != 0.0) {
          count += 1.0;
          e[i] = chi->omega * (b[i] - a[i]) / (d[i] - c[i]);
        } else
          e[i] = 0.0;
      }
      test = sqrt(dp(e, e)) / count;
      /* printf("heur:M/x=%e\n", test); */
      if (test > 0.01)
        cpi->mass_x = 1;

      dt_p = scalar_perturbation(chi->bif_param, cgi->perturb);
      assign_bif_parameter_conwrap(chi->bif_param + dt_p);
      mass_matrix_fill_conwrap(x, x_tmp);
      mass_matvec_mult_conwrap(chi->y_vec, b);
      matvec_mult_conwrap(chi->z_vec, d);
      assign_bif_parameter_conwrap(chi->bif_param);
      count = 0.0;
      for (i = 0; i < cgi->numOwnedUnks; i++) {
        if (d[i] - c[i] != 0.0) {
          count += 1.0;
          e[i] = chi->omega * (b[i] - a[i]) / (d[i] - c[i]);
        } else
          e[i] = 0.0;
      }
      count = gsum_double_conwrap(count);
      test = sqrt(dp(e, e)) / count;
      /* printf("heur:M/param=%e\n", test); */
      if (test > 0.01)
        cpi->mass_param = 1;

      if (cgi->printproc > 7) {
        if (cpi->mass_x == 1)
          printf("\n\tMass Matrix is a Function of Solution Vector!\n");
        else
          printf("\n\tMass Matrix is Independent of Solution Vector!\n");
        if (cpi->mass_param == 1)
          printf("\tMass Matrix is a Function of Bifurcation Parameter!\n");
        else
          printf("\tMass Matrix is Independent of Bifurcation Parameter!\n");
      }
    } else {
      if (cgi->printproc > 1)
        printf("\n\tERROR: Mass Matrix Derivatives Flag not set in\n");
      if (cgi->printproc > 1)
        printf("\t solve_continuation() \n");
      exit(-1);
    }

    first = FALSE;
  }

  if (cgi->printproc > 7) {
    if (cpi->mass_x == 1)
      printf("\n\tdM/dx is included in Komplex solves\n");
    else
      printf("\n\tdM/dx is NOT included in Komplex solves\n");
    if (cpi->mass_param == 1)
      printf("\tdM/d(param) is included in Komplex solves\n");
    else
      printf("\tdM/d(param) is NOT included in Komplex solves\n");
  }

  /* If Mass Matrix is not constant, allocate extra work arrays */
  if (cpi->mass_param || cpi->mass_x) {
    M_1 = alloc_vec();
    M_2 = alloc_vec();
  }

  if (first2) {
    /* calculate variable averages, to be used in ltransnorm calls below */
    /* calc_scale_vec_conwrap(x, cpi->scale_vec, cgi->numUnks); */
    for (i = 0; i < cgi->numOwnedUnks; i++)
      cpi->scale_vec[i] = chi->y_vec[i];

    /* Make sure the eigenvectors y_vec and z_vec are orthogonal */
    alpha = ltransnorm(chi->y_vec, cpi->scale_vec);
    beta = ltransnorm(chi->z_vec, cpi->scale_vec);
    for (i = 0; i < cgi->numOwnedUnks; i++) {
      y_tmp = chi->y_vec[i];
      z_tmp = chi->z_vec[i];
      chi->y_vec[i] = (alpha * y_tmp + beta * z_tmp) / (alpha * alpha + beta * beta);
      chi->z_vec[i] = -(beta * y_tmp - alpha * z_tmp) / (alpha * alpha + beta * beta);
    }
    first2 = FALSE;
  }

  /* construct "a" vector from delta_x */

  if (cgi->printproc > 4)
    printf("\n\tHopf Continuation: Constructed *a* vector!\n");
  for (i = 0; i < cgi->numUnks; i++)
    a[i] = -delta_x[i];

  /*
   * Next, "b" is calculated. This is exactly
   * the same linear solve as a continuation predictor step except
   * the bif_param is perturbed, not the continuation param.
   */

  if (cgi->printproc > 4)
    printf("\n\tHopf Continuation: Calculating *b* vector!\n\n");

  calc_rhs_continuation(TP_CONT_SOL2, x, b, NULL, NULL, NULL, chi->bif_param, cgi->perturb, NULL,
                        cgi->numUnks, cgi->numOwnedUnks);

  i = linear_solver_conwrap(b, CHECK_JACOBIAN, NULL);

  /* Fill the Mass Matrix for RHS and Komplex Solves*/
  mass_matrix_fill_conwrap(x, x_tmp);

  /*
   * Next, "c" and "d" are calculated as a function of y_vec, z_vec,
   * and the Mass Matrix M. c and d hold the rhs vector on input to
   * solver and have solution on exit.
   */
  if (cgi->printproc > 4)
    printf("\n\tHopf Continuation: Calculating *c* and *d* vectors"
           " (Komplex solve)!\n");

  calc_rhs_continuation(HP_CONT_SOL3, x, c, chi->z_vec, chi->y_vec, x_tmp, chi->bif_param,
                        cgi->perturb, d, cgi->numUnks, cgi->numOwnedUnks);

  i = komplex_linear_solver_conwrap(c, d, NEW_JACOBIAN, &chi->omega, x_tmp);

  /*
   * Next, "e" and "f" are calculated as a function of a, y_vec, and z_vec.
   * e and f hold the rhs vector on input to  solver and have solution
   * on exit.
   */
  if (cgi->printproc > 4)
    printf("\n\tHopf Continuation: Calculating *e* and *f* vectors"
           " (Komplex solve)!\n");
  calc_rhs_continuation(TP_CONT_SOL3, x, e, a, cpi->scale_vec, x_tmp, chi->bif_param, cgi->perturb,
                        chi->y_vec, cgi->numUnks, cgi->numOwnedUnks);
  calc_rhs_continuation(TP_CONT_SOL3, x, f, a, cpi->scale_vec, x_tmp, chi->bif_param, cgi->perturb,
                        chi->z_vec, cgi->numUnks, cgi->numOwnedUnks);

  /* Get null vector residual now. */

  RayQ = null_vector_resid(0.0, chi->omega, chi->y_vec, chi->z_vec, TRUE);

  if (cpi->mass_x) {
    calc_rhs_continuation(HP_CONT_DMDX, x, M_1, a, cpi->scale_vec, x_tmp, chi->bif_param,
                          cgi->perturb, chi->z_vec, cgi->numUnks, cgi->numOwnedUnks);
    calc_rhs_continuation(HP_CONT_DMDX, x, M_2, a, cpi->scale_vec, x_tmp, chi->bif_param,
                          cgi->perturb, chi->y_vec, cgi->numUnks, cgi->numOwnedUnks);
    for (i = 0; i < cgi->numUnks; i++) {
      e[i] += M_1[i] * chi->omega;
      f[i] += -M_2[i] * chi->omega;
    }
  }
  i = komplex_linear_solver_conwrap(e, f, OLD_JACOBIAN, &chi->omega, x_tmp);

  /*
   * Next, "g" and "h" are calculated as a function of a, y_vec, and z_vec.
   * g and h hold the rhs vector on input to  solver and have solution
   * on exit.
   */
  if (cgi->printproc > 4)
    printf("\n\tHopf Continuation: Calculating *g* and *h* vectors"
           " (Komplex solve)!\n");
  calc_rhs_continuation(TP_CONT_SOL4, x, g, b, cpi->scale_vec, x_tmp, chi->bif_param, cgi->perturb,
                        chi->y_vec, cgi->numUnks, cgi->numOwnedUnks);
  calc_rhs_continuation(TP_CONT_SOL4, x, h, b, cpi->scale_vec, x_tmp, chi->bif_param, cgi->perturb,
                        chi->z_vec, cgi->numUnks, cgi->numOwnedUnks);
  if (cpi->mass_x) {
    calc_rhs_continuation(HP_CONT_DMDX, x, M_1, b, cpi->scale_vec, x_tmp, chi->bif_param,
                          cgi->perturb, chi->z_vec, cgi->numUnks, cgi->numOwnedUnks);
    calc_rhs_continuation(HP_CONT_DMDX, x, M_2, b, cpi->scale_vec, x_tmp, chi->bif_param,
                          cgi->perturb, chi->y_vec, cgi->numUnks, cgi->numOwnedUnks);
    for (i = 0; i < cgi->numUnks; i++) {
      g[i] += M_1[i] * chi->omega;
      h[i] += -M_2[i] * chi->omega;
    }
  }
  if (cpi->mass_param) {
    vec_copy(x, M_1);
    vec_copy(x, M_2);
    calc_rhs_continuation(HP_CONT_DMDPARAM, x, M_1, a, cpi->scale_vec, x_tmp, chi->bif_param,
                          cgi->perturb, chi->z_vec, cgi->numUnks, cgi->numOwnedUnks);
    calc_rhs_continuation(HP_CONT_DMDPARAM, x, M_2, a, cpi->scale_vec, x_tmp, chi->bif_param,
                          cgi->perturb, chi->y_vec, cgi->numUnks, cgi->numOwnedUnks);
    for (i = 0; i < cgi->numUnks; i++) {
      g[i] += M_1[i] * chi->omega;
      h[i] += -M_2[i] * chi->omega;
    }
  }
  i = komplex_linear_solver_conwrap(g, h, OLD_JACOBIAN_DESTROY, &chi->omega, x_tmp);

  if (cgi->printproc > 4)
    printf("\n\tHopf Continuation: Finished Komplex solves, "
           "Updating solution \n\n");

  /*
   * Calculate the updates and scaled updates to bif_param,
   * omega, y_vec, z_vec, and x.
   */

  /* Calculate Update parameters */
  ltc = ltransnorm(c, cpi->scale_vec);
  ltd = ltransnorm(d, cpi->scale_vec);
  lte = ltransnorm(e, cpi->scale_vec);
  ltf = ltransnorm(f, cpi->scale_vec);
  ltg = ltransnorm(g, cpi->scale_vec);
  lth = ltransnorm(h, cpi->scale_vec);

  /* Update turning point param */
  dt_p = (ltc * ltf - lte * ltd + ltd) / (ltd * ltg - ltc * lth);
  chi->bif_param += dt_p;
  assign_bif_parameter_conwrap(chi->bif_param);
  param_update = fabs(dt_p) / (reltol * fabs(chi->bif_param) + abstol);

  /* Update imaginary eigenvalue omega */
  delta_omega = (lth * dt_p + ltf) / ltd;
  chi->omega += delta_omega;
  omega_update = fabs(delta_omega) / (reltol * fabs(chi->omega) + abstol);

  /* Update delta_x - state variables: Negative of update vector here */
  x_update = 0.0;
  for (i = 0; i < cgi->numOwnedUnks; i++)
    delta_x[i] = -(a[i] + dt_p * b[i]);
  for (i = 0; i < cgi->numOwnedUnks; i++)
    x_update += fabs(delta_x[i]) / (fabs(x[i]) * reltol + abstol) * fabs(delta_x[i]) /
                (fabs(x[i]) * reltol + abstol);

  /* Update eigenvectors for real and imaginary parts */
  y_update = 0.0;
  z_update = 0.0;
  for (i = 0; i < cgi->numOwnedUnks; i++) {

    e[i] = -chi->y_vec[i] + e[i] + g[i] * dt_p - c[i] * delta_omega;
    f[i] = -chi->z_vec[i] + f[i] + h[i] * dt_p - d[i] * delta_omega;
    delta_y = e[i];
    delta_z = f[i];
    y_update += fabs(delta_y) / (fabs(chi->y_vec[i]) * reltol + abstol) * fabs(delta_y) /
                (fabs(chi->y_vec[i]) * reltol + abstol);
    z_update += fabs(delta_z) / (fabs(chi->z_vec[i]) * reltol + abstol) * fabs(delta_z) /
                (fabs(chi->z_vec[i]) * reltol + abstol);
    chi->y_vec[i] += delta_y;
    chi->z_vec[i] += delta_z;
  }
  gnum_unks = gsum_double_conwrap(1.0 * cgi->numOwnedUnks);
  y_update = sqrt(gsum_double_conwrap(y_update) / gnum_unks);
  z_update = sqrt(gsum_double_conwrap(z_update) / gnum_unks);
  x_update = sqrt(gsum_double_conwrap(x_update) / gnum_unks);

  /*
   * Check whether or not continuation param, omega and eigenvectors
   * are converged
   */

  if (RayQ == -1.0) {
    if (cgi->printproc > 1)
      printf("\n\tRayleigh Quotient Error: zero denom\n");
  } else {
    if (cgi->printproc > 4)
      printf("\n\tRayleigh Quotient before updates: %g\n", RayQ);
  }

  if (cgi->printproc > 4) {
    printf("\n\tHopf Continuation: Convergence Criteria\n");
    printf("\tVariable     Scaled Update (<1)  Unscaled Update  New Value\n");
    printf("\t***********************************************************\n");
    printf("\tX vector     %e\n", x_update);
    printf("\tparameter    %e        %e    %g\n", param_update, dt_p, chi->bif_param);
    printf("\tomega        %e        %e    %g\n", omega_update, delta_omega, chi->omega);
    printf("\tY vector     %e\n", y_update);
    printf("\tZ vector     %e\n", z_update);
    printf("\t***********************************************************\n");
  }

  free_vec(&a);
  free_vec(&b);
  free_vec(&c);
  free_vec(&d);
  free_vec(&e);
  free_vec(&f);
  free_vec(&g);
  free_vec(&h);
  free_vec(&x_tmp);
  if (cpi->mass_param || cpi->mass_x) {
    free_vec(&M_1);
    free_vec(&M_2);
  }

  /* return convergence status of the parameters and vectors */
  if ((x_update < 1.0) && (param_update < 1.0) && (omega_update < 1.0) && (y_update < 100.0) &&
      (z_update < 100.0)) {
    first2 = TRUE;
    return (TRUE);
  } else
    return (FALSE);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
static int phase_transition_alg(
    double *x, double *delta_x, struct con_struct *con, double reltol, double abstol)

/*
 * Algorithm for locking on to a phase transition.
 * This involves finding two solutions that exist at the
 * same parameters and at the same energy.
 */

{
  /* shorthand for long con sub-structure */
  struct general_info_struct *cgi = &(con->general_info);

  double *a, *b, *c, *d, *x2, g, dg_dtp, dt_p;
  double *x_tmp, *x2_tmp, dgac, dgbd, eps;
  double tmp = 0, gnum_unks, x2_update, param_update;
  int i;

  /*************************** BEGIN EXECUTION *******************************/

  x2 = con->phase_transition_info.x2;
  a = alloc_vec();
  b = alloc_vec();
  c = alloc_vec();
  d = alloc_vec();
  x_tmp = alloc_vec();
  x2_tmp = alloc_vec();

  for (i = 0; i < cgi->numUnks; i++) {
    a[i] = -delta_x[i];
  }

  /* Construct b vector using tangent calculation */

  calc_rhs_continuation(TP_CONT_SOL2, x, b, NULL, NULL, NULL, con->phase_transition_info.bif_param,
                        cgi->perturb, NULL, cgi->numUnks, cgi->numOwnedUnks);

  i = linear_solver_conwrap(b, OLD_JACOBIAN, NULL);

  /* Now do Newton iteration on second vector for c */

  matrix_residual_fill_conwrap(x2, c, RHS_MATRIX);

  i = linear_solver_conwrap(c, NEW_JACOBIAN, NULL);

  for (i = 0; i < cgi->numUnks; i++)
    c[i] *= -1.0;

  /* Construct d vector using tangent calculation */

  calc_rhs_continuation(TP_CONT_SOL2, x, d, NULL, NULL, NULL, con->phase_transition_info.bif_param,
                        cgi->perturb, NULL, cgi->numUnks, cgi->numOwnedUnks);

  i = linear_solver_conwrap(d, CHECK_JACOBIAN, NULL);

  /* Now start bordering alg for equal-energy constraint
   * g is energy difference to be driven to zero,
   * dg_dtp is the derivative of the energy w.r.t. the parameter
   * dgac is the directional derivative in of g the a:c direction
   * dgbd is the directional derivative of g in the b:d direction
   */

  g = free_energy_diff_conwrap(x, x2);

  dt_p = scalar_perturbation(con->phase_transition_info.bif_param, cgi->perturb);

  assign_bif_parameter_conwrap(con->phase_transition_info.bif_param + dt_p);

  dg_dtp = (free_energy_diff_conwrap(x, x2) - g) / dt_p;

  assign_bif_parameter_conwrap(con->phase_transition_info.bif_param);

  eps = cgi->perturb *
        sqrt(dp(x, x) / (dp(a, a) + cgi->perturb) + dp(x2, x2) / (dp(c, c) + cgi->perturb));

  for (i = 0; i < cgi->numUnks; i++) {
    x_tmp[i] = x[i] + eps * a[i];
    x2_tmp[i] = x2[i] + eps * c[i];
  }
  dgac = (free_energy_diff_conwrap(x_tmp, x2_tmp) - g) / eps;

  eps = cgi->perturb *
        sqrt(dp(x, x) / (dp(b, b) + cgi->perturb) + dp(x2, x2) / (dp(d, d) + cgi->perturb));

  for (i = 0; i < cgi->numUnks; i++) {
    x_tmp[i] = x[i] + eps * b[i];
    x2_tmp[i] = x2[i] + eps * d[i];
  }
  dgbd = (free_energy_diff_conwrap(x_tmp, x2_tmp) - g) / eps;

  dt_p = -(g + dgac) / (dgbd + dg_dtp);

  /* update continuation parameter */

  con->phase_transition_info.bif_param += dt_p;
  assign_bif_parameter_conwrap(con->phase_transition_info.bif_param);

  /* update x2, checking for convergence */

  for (i = 0; i < cgi->numUnks; i++)
    delta_x[i] = c[i] + dt_p * d[i];

  /* get update norm to check for convergence */

  x2_update = 0.0;
  for (i = 0; i < cgi->numUnks; i++) {
    tmp += delta_x[i] / (fabs(x2[i]) * reltol + abstol);
    x2_update += tmp * tmp;
  }
  gnum_unks = gsum_double_conwrap(1.0 * cgi->numOwnedUnks);
  x2_update = sqrt(gsum_double_conwrap(x2_update / gnum_unks));

  param_update = fabs(dt_p) / (reltol * fabs(con->phase_transition_info.bif_param) + abstol);

  /* update x2 and modify delta_x for update in usual Newton routine */

  for (i = 0; i < cgi->numUnks; i++)
    x2[i] += delta_x[i];
  for (i = 0; i < cgi->numUnks; i++)
    delta_x[i] = -a[i] - dt_p * b[i];

  free_vec(&a);
  free_vec(&b);
  free_vec(&c);
  free_vec(&d);
  free_vec(&x_tmp);
  free_vec(&x2_tmp);

  if (cgi->printproc > 4) {
    printf("\tPhase transition algorithm scaled parameter update  : %g\n", param_update);
    printf("\tPhase transition algorithm scaled solution #2 update: %g\n", x2_update);
  }

  /* x2_update and param_update must be < 1 for convergence */

  if (x2_update < 1.0 && param_update < 1.0)
    return (TRUE);
  else
    return (FALSE);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static double scalar_perturbation(double param, double perturb)
/* Routine for picking a perturbation to a scalar for simple forward
 * difference derivatives.
 */
{
  double delta;

  delta = perturb * (perturb + fabs(param));
  if (delta < 1.0e-8) {
    delta = 1.0e-8;
  }
  return delta;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void calc_rhs_continuation(int rhs_type,
                           double *x,
                           double *resid_vector,
                           double *ab_vec,
                           double *scale_vec,
                           double *x_tmp,
                           double param,
                           double perturb,
                           double *r_vec,
                           int numUnks,
                           int numOwnedUnks)

/* routine to pre-calculate the non-standard right-hand-sides for solves */
/* of continuation runs. This routine returns the matrix fill time.      */
/* The matrix at the current solution is also calculated                 */

{
  double abdp, dc_p, dc_p1, *resid_delta;
  int i;

  /* Allocate and initialize resid_delta */

  vec_init(resid_vector);
  resid_delta = alloc_vec();

  /*
   * For the first 2 options, resid_delta is the right hand side at the
   * perturbed value of the continuation parameter
   */

  if (rhs_type == CONT_TANGENT) {

    dc_p = scalar_perturbation(param, perturb);

    assign_parameter_conwrap(param + dc_p);

    matrix_residual_fill_conwrap(x, resid_delta, RHS_ONLY);

    assign_parameter_conwrap(param);

    matrix_residual_fill_conwrap(x, resid_vector, RHS_MATRIX);

    for (i = 0; i < numUnks; i++)
      resid_vector[i] = (resid_delta[i] - resid_vector[i]) / dc_p;
  }

  else if (rhs_type == ARC_CONT_SOL2) {

    dc_p = scalar_perturbation(param, perturb);

    assign_parameter_conwrap(param + dc_p);

    matrix_residual_fill_conwrap(x, resid_delta, RHS_ONLY);

    assign_parameter_conwrap(param);

    matrix_residual_fill_conwrap(x, resid_vector, RHS_MATRIX_SAVE);

    for (i = 0; i < numUnks; i++)
      resid_vector[i] = (resid_delta[i] - resid_vector[i]) / dc_p;
  }

  /*
   * For TP_CONT_SOL2, resid_delta is the right hand side at the
   * perturbed value of the turning point parameter
   */

  else if (rhs_type == TP_CONT_SOL2) {

    dc_p = scalar_perturbation(param, perturb);

    assign_bif_parameter_conwrap(param + dc_p);

    matrix_residual_fill_conwrap(x, resid_delta, RHS_ONLY);

    assign_bif_parameter_conwrap(param);

    matrix_residual_fill_conwrap(x, resid_vector, RHS_MATRIX_SAVE);

    for (i = 0; i < numUnks; i++)
      resid_vector[i] = -(resid_delta[i] - resid_vector[i]) / dc_p;
  }

  /*
   * For TP_CONT_SOL3, resid_delta is the matrix at perturbed value of
   * solution in the direction of ab_vec multiplied by the vector r_vec
   */

  else if (rhs_type == TP_CONT_SOL3) {

    abdp = dp(ab_vec, ab_vec);

    if (fabs(abdp) < 1.0e-99) {
      printf(" ab_vec too small: %g, try a larger perturbation!\n", abdp);
      exit(-1);
    } else {
      dc_p1 = perturb * (perturb + sqrt(dp(x, x) / dp(ab_vec, ab_vec)));
    }

    for (i = 0; i < numUnks; i++)
      x_tmp[i] = x[i] + dc_p1 * ab_vec[i];

    matrix_residual_fill_conwrap(x_tmp, resid_delta, RHS_MATRIX);

    matvec_mult_conwrap(r_vec, resid_delta);

    matrix_residual_fill_conwrap(x, resid_vector, RECOVER_MATRIX);

    matvec_mult_conwrap(r_vec, resid_vector);

    /*
    printf("In SOL3: n.j.n = %g  %g\n", dp(r_vec, resid_vector),
    dp(r_vec, resid_vector,numOwnedUnks)/dp(r_vec, r_vec));
    */

    for (i = 0; i < numUnks; i++)
      resid_vector[i] = -AGS_option * resid_vector[i] - (resid_delta[i] - resid_vector[i]) / dc_p1;
  }

  /*
   * For TP_CONT_SOL4, resid_delta is the matrix at perturbed value of
   * solution in the direction of ab_vec multiplied by the vector r_vec
   * plus the matrix perturbed in the direction of the bif_parameter
   * multiplied by the vector r_vec.
   */

  else if (rhs_type == TP_CONT_SOL4) {

    dc_p = scalar_perturbation(param, perturb);

    abdp = dp(ab_vec, ab_vec);

    if (fabs(abdp) < 1.0e-99) {
      printf(" ab_vec too small: %g, try a larger perturbation!\n", abdp);
      exit(-1);
    } else {
      dc_p1 = perturb * (perturb + sqrt(dp(x, x) / dp(ab_vec, ab_vec)));
    }

    for (i = 0; i < numUnks; i++)
      x_tmp[i] = x[i] + dc_p1 * ab_vec[i];

    matrix_residual_fill_conwrap(x_tmp, resid_delta, MATRIX_ONLY);

    matvec_mult_conwrap(r_vec, x_tmp);

    assign_bif_parameter_conwrap(param + dc_p);

    matrix_residual_fill_conwrap(x, resid_delta, MATRIX_ONLY);

    assign_bif_parameter_conwrap(param);

    matvec_mult_conwrap(r_vec, resid_delta);

    /*
     * Two numerical differences with different perturbations are summed
     * together
     */

    for (i = 0; i < numUnks; i++)
      resid_delta[i] = resid_delta[i] / dc_p + x_tmp[i] / dc_p1;

    /* also sum together the perturbations */

    dc_p = 1.0 / (1.0 / dc_p1 + 1.0 / dc_p);

    matrix_residual_fill_conwrap(x, resid_vector, RECOVER_MATRIX);

    matvec_mult_conwrap(r_vec, resid_vector);

    for (i = 0; i < numUnks; i++)
      resid_vector[i] = -resid_delta[i] + resid_vector[i] / dc_p;

  }

  else if (rhs_type == HP_CONT_SOL3) {
    mass_matvec_mult_conwrap(ab_vec, resid_vector);
    mass_matvec_mult_conwrap(scale_vec, resid_delta);
    for (i = 0; i < numUnks; i++)
      r_vec[i] = -resid_delta[i];

  }

  else if (rhs_type == HP_CONT_DMDX) {

    abdp = dp(ab_vec, ab_vec);

    if (fabs(abdp) < 1.0e-99) {
      printf(" ab_vec too small: %g, try a larger perturbation !\n", abdp);
      exit(-1);
    } else {
      dc_p1 = perturb * (perturb + sqrt(dp(x, x) / dp(ab_vec, ab_vec)));
    }

    for (i = 0; i < numUnks; i++)
      x_tmp[i] = x[i] + dc_p1 * ab_vec[i];
    mass_matrix_fill_conwrap(x_tmp, resid_delta);
    mass_matvec_mult_conwrap(r_vec, resid_delta);
    mass_matrix_fill_conwrap(x, resid_vector);
    mass_matvec_mult_conwrap(r_vec, resid_vector);

    for (i = 0; i < numUnks; i++)
      resid_vector[i] = -(resid_delta[i] - resid_vector[i]) / dc_p1;
  }

  else if (rhs_type == HP_CONT_DMDPARAM) {

    dc_p = scalar_perturbation(param, perturb);
    assign_bif_parameter_conwrap(param + dc_p);
    mass_matrix_fill_conwrap(x, resid_vector);
    mass_matvec_mult_conwrap(r_vec, resid_delta);
    assign_bif_parameter_conwrap(param);
    mass_matrix_fill_conwrap(x, resid_vector);
    mass_matvec_mult_conwrap(r_vec, resid_vector);

    for (i = 0; i < numUnks; i++)
      resid_vector[i] = -(resid_delta[i] - resid_vector[i]) / dc_p;
  }

  free_vec(&resid_delta);

  return;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int continuation_hook(double *x, double *delta_x, void *con_void, double reltol, double abstol) {
  /* return TRUE for non-continuation probelms, as well as zero and
   * first order continuation runs
   */

  int converged = TRUE;
  struct con_struct *con = (struct con_struct *)con_void;

  /* shorthand for long con sub-structure */
  struct general_info_struct *cgi = &(con->general_info);
  struct private_info_struct *cpi = &(con->private_info);

  /*
   * If arc_length continuation is being done, then branch here
   * to do the rest of the bordering algorithm. This involves a
   * second matrix fill and linear solve. (See JNS thesis
   * eqns 7.18-7.19 for equations, 7.21-7.23 for method.
   */

  if (cgi->method == ARC_LENGTH_CONTINUATION && cpi->step_num > 0)
    converged = arc_length_bordering_alg(x, delta_x, con, reltol, abstol);

  /*
   * If turning point calculations are being done, then branch here
   * to do the rest of the algorithm. This involves 6 more matrix fills
   * and 3 more linear solves in the current implementation.
   */

  else if (cgi->method == TURNING_POINT_CONTINUATION)
    converged = turning_point_alg(x, delta_x, con, reltol, abstol);

  /*
   * If pitchfork bifurcation calculations are being done, then branch
   * here to do the rest of the algorithm.
   */

  else if (cgi->method == PITCHFORK_CONTINUATION)
    converged = pitchfork_alg(x, delta_x, con, reltol, abstol);

  /*
   * If Hopf bifurcation calculations are being done, then branch
   * here to do the rest of the algorithm.
   */

  else if (cgi->method == HOPF_CONTINUATION)
    converged = hopf_alg(x, delta_x, con, reltol, abstol);

  /*
   * If phase_transition tracking is being done, then branch here
   * to do the rest of the algorithm. This involves
   * 3 more linear solves in the current implementation.
   */

  else if (cgi->method == PHASE_TRANSITION_CONTINUATION)
    converged = phase_transition_alg(x, delta_x, con, reltol, abstol);

  /* No boredering alg needed for zero and first order continuations */

  else if (cgi->method == ZERO_ORDER_CONTINUATION || cgi->method == FIRST_ORDER_CONTINUATION ||
           cgi->method == LOCA_LSA_ONLY ||
           (cgi->method == ARC_LENGTH_CONTINUATION && cpi->step_num == 0)) {
    converged = TRUE;
  }

  /* perform error check */

  else {
    if (cgi->printproc > 1)
      printf("\nERROR continuation_hook: Unknown method %d\n", cgi->method);
    exit(-1);
  }

  /* Turn off flag that is true for first Newton iteration of each solve */

  cpi->first_iter = FALSE;

  /* Return flag indicating convergence of the Newton iter */

  return converged;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
