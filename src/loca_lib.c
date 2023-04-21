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
 -----------------------------------------------------------------------------
   LOCA 1.0: Library of Continuation Algorithms
   Copyright (C) 2001, Sandia National Laboratories

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 -----------------------------------------------------------------------------
*/
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "el_quality.h"
#include "exo_struct.h"
#include "loca_const.h"
#include "loca_util_const.h"
#include "sl_util_structs.h"
#include "std.h"

/*****************************************************************************/
/********************* Static Functions Used *********************************/
/*****************************************************************************/

struct arc_scale_struct {
  double dx0;         /* Initial scale factor for solution             */
  double dx_fac;      /* scale factor for solution                     */
  double dx_fac_max;  /* maximum allowable value for arc->dx_fac       */
  double umag2;       /* square of scaled solution vector magnitude    */
  double dx_fac_old;  /* scale factor from previous step               */
  double dp_ds_limit; /* dp_ds value at which to default to maximum    */
                      /* value of arc->dx_fac = arc->dx_fac_max        */
  double dp_ds_goal;  /* dp_ds value to reset to by rescaling          */
  double dp_ds_old;   /* dp_ds value from previous step                */
  double ds_fac;      /* arc length scale factor                       */
};

static double solution_scale(struct con_struct *con, struct arc_scale_struct *arc);
double scaled_dot_prod(double *x, double *y, double *sc, int n);
static void print_cont_step1(int order, double step, double step_old, struct con_struct *con);
static void print_cont_step2(int order, double step, struct con_struct *con);
static void print_cont_step_fail(int order, double step, struct con_struct *con);
static double simple_step_control(int num_newt_conv, int max_Newton_steps, double step_ctrl);
static void print_line(char *charstr, int ntimes);
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int con_lib(struct con_struct *con, Exo_DB *exo, struct GomaLinearSolverData *ams)

/*****************************************************************************
 *
 *  Input Variables:
 *
 ******************************************************************************/

{

  /* Local Variables */

  int n;     /* Loop index                          */
  int order; /* Continuation order flag:
                       0 - zero-order continuation
                       1 - first-order continuation
                       2 - arc-length continuation
                     This flag is always 0 on the first
                     solve, and 0 for turning point or any
                     other special continuation */
  int i;
  int num_newt_conv = 0; /* Number of newton iterations to reach
                                 convergence for last nonlinear solve
                                 -- used to pick next step size
                                 ALSO error flag, when < 0            */
  const char *yo = "con_lib";
  int sn_old = 0, sn_new = 0; /* Sign of cpi->dp_ds, to check for a
                                 turning point                        */

  int tan_flag; /* Set to zero when tang_factor is smaller
                   than step limit specified: step will be
                   halved and repeated even though converged */

  /* These quantities are used
     for arc length step control:              */
  double delta_s, end_passed, max_step, temp_step;
  double *x_tang_old, *x2_old, *y_vec_old, *z_vec_old, bif_param_old = 0.0;
  double ds_ratio = 1.0, first_arc_step = 0.0, arc_step_adj = 0.0;
  double tang_factor = 1.0;
  double step;                       /* step size of continuation parameter                  */
  double step_old;                   /* old step size of continuation parameter              */
  double omega_old = 0.0;            /* old Hopf tracking frequency                          */
  struct arc_scale_struct arc_scale; /* allocate scaling params in struct   */
  struct arc_scale_struct *arc = &(arc_scale); /* pointer to arc_scale struct */

  /* shorthand to some commonly used structures */
  struct general_info_struct *cgi = &(con->general_info);
  struct stepping_info_struct *csi = &(con->stepping_info);
  struct arclength_info_struct *cai = &(con->arclength_info);
  struct private_info_struct *cpi = &(con->private_info);
  double *x = cgi->x;
  memset((void **)arc, 0, sizeof(struct arc_scale_struct));

  /******************************* First Executable Statment *****************/
  if (cgi->printproc > 1)
    printf("\tLOCA v1.0, Copyright 2001 Sandia Corporation\n\n");

  /* Send vector length information to utilities file, so routines such
   * as vector allocs and dot products won't need length information    */

  initialize_util_routines(cgi->numOwnedUnks, cgi->numUnks);

  /*
   * Initialize arrays to store predicted and old solutions. Save the current
   * solution, x, into cpi->x_old
   */

  cpi->x_old = alloc_vec();
  cpi->scale_vec = alloc_vec();
  cpi->param_old = cgi->param;
  step = 0.0;
  step_old = 0.0;
  vec_copy(x, cpi->x_old);

  switch (cgi->method) {
  case ZERO_ORDER_CONTINUATION:
  case FIRST_ORDER_CONTINUATION:
    cpi->x_tang = alloc_vec();
    break;
  case ARC_LENGTH_CONTINUATION:
    cpi->arc_step = 0.0;
    cpi->x_tang = alloc_vec();
    x_tang_old = alloc_vec();
    break;
  case TURNING_POINT_CONTINUATION:
  case PITCHFORK_CONTINUATION:
    cpi->x_tang = alloc_vec();
    x_tang_old = alloc_vec();
    break;
  case HOPF_CONTINUATION:
    omega_old = con->hopf_info.omega;
    y_vec_old = alloc_vec();
    z_vec_old = alloc_vec();
    break;
  case PHASE_TRANSITION_CONTINUATION:
    x2_old = alloc_vec();
    break;
  }

  /*
   * Initialize variables used in arc length step control
   */

  if (cgi->method == ARC_LENGTH_CONTINUATION) {
    arc->dx_fac = 1.0;
    if ((cai->dp_ds2_goal) < 1.0e-6) {
      arc->dx_fac = 100.0;
    }
    arc->dx0 = arc->dx_fac;
    arc->dx_fac_max = 1e+8;
    arc->dp_ds_goal = sqrt(cai->dp_ds2_goal);
    arc->dp_ds_limit = arc->dx_fac_max * arc->dp_ds_goal /
                       sqrt(1.0 + cai->dp_ds2_goal * (arc->dx_fac_max * arc->dx_fac_max - 1.0));
    if (cai->dp_ds_max < arc->dp_ds_goal)
      cai->dp_ds_max = arc->dp_ds_goal;
  }

  /* Adjust the BCs/Properties/whatever that the con param really represents */

  assign_parameter_conwrap(cgi->param);
  if (cgi->method == TURNING_POINT_CONTINUATION)
    assign_bif_parameter_conwrap(con->turning_point_info.bif_param);
  else if (cgi->method == PITCHFORK_CONTINUATION)
    assign_bif_parameter_conwrap(con->pitchfork_info.bif_param);
  else if (cgi->method == HOPF_CONTINUATION)
    assign_bif_parameter_conwrap(con->hopf_info.bif_param);
  else if (cgi->method == PHASE_TRANSITION_CONTINUATION)
    assign_bif_parameter_conwrap(con->phase_transition_info.bif_param);

  /* In tp_continuation, perturb initial guess off of potential singularity*/

  if (cgi->method == TURNING_POINT_CONTINUATION || cgi->method == PITCHFORK_CONTINUATION) {
    if (cgi->printproc > 4)
      printf("\tcon_lib: Adding random"
             " perturbation for continuation\n");

    perturb_solution_conwrap(x, cpi->x_old, cpi->scale_vec, cgi->numOwnedUnks);
  }

  /*
   * Print out general time integration information
   */

  if (cgi->printproc > 1) {
    printf("\n");
    print_line("~", 80);
    print_line("~", 80);
    printf("%s: Start Continuation\n", yo);
    printf("\tInitial step size = %e \n", csi->first_step);
    printf("\tMax number of continuation steps = %d\n", csi->max_steps);
    printf("\tMax parameter value = %g\n", csi->max_param);
    print_line("~", 80);
    print_line("~", 80);
    printf("\n");
  }

  /***************************** CONTINUATION LOOP **************************/

  /*
   * Initialize the time step counter to 0. Set order flag to zero-order
   * continuation through first solution.
   */

  cpi->step_num = order = 0;
  cpi->nstep = csi->base_step;
  csi->last_step = FALSE;

  /*
   * Loop through a number of continuation steps - note the loop index may not
   * represent the actual step counter, due to failed steps.  The convention
   * here is that solution 0 is not a step, so there will be
   * csi->max_steps+1 solutions, numbered 0 through
   * con->stepping_info.max_steps.
   */

  for (n = 0; n <= csi->max_steps; n++) {

    /*
     * Print out an initial statement about the step.
     */

    if (cgi->printproc > 1)
      print_cont_step1(order, step, step_old, con);

    /*
     *  Set flag for new solve to detect first Newton iter
     */

    cpi->first_iter = TRUE;

    /*
     * Solve the system of equations at the current step.
     * Note - x is considered to be updated, on return from this
     * solution.
     */

    num_newt_conv = nonlinear_solver_conwrap(x, (void *)con, cpi->step_num, cgi->param, step);

    /*
     * If tan_factor changes too much, tan_flag tells to halve step & reset.
     */

    tan_flag = TRUE;

    /*
     * Check for convergence
     */

    if (num_newt_conv < 0) {

      /*
       * Convergence Failure!
       *
       * If initial guess did not converge, abort entire continuation run
       */

      if (cpi->step_num == 0) {
        n = csi->max_steps; /* Force IO and exit  */
        cpi->step_num = -1; /* Set failure flag */

        if (cgi->printproc > 1) {
          printf("\n\t %s: INITIAL GUESS DID NOT CONVERGE:", yo);
          printf("\n\t\t\t     ABORTING CONTINUATION RUN\n");
        }
      }

      /*
       * If this convergence failure wasn't the first or last step, cut step
       * size in half and calculate a new initial guess.  New guess is
       * the old solution plus the tangent to the previous prediction
       * times the halved step size.
       */

      else {

        if (n < csi->max_steps) {
          step *= 0.5;

          /* Check for step size too small, abort if below min_delta_p */
          if (fabs(step) < csi->min_delta_p) {
            n = csi->max_steps;
            cpi->step_num = -1;
            if (cgi->printproc > 1) {
              printf("\n\t %s: CONTINUATION STEP SIZE TOO SMALL:", yo);
              printf("\n\t\t\t     ABORTING CONTINUATION RUN\n");
            }
          }

          cgi->param = cpi->param_old + step;
          assign_parameter_conwrap(cgi->param);

          switch (order) {
          case 0:
            vec_copy(cpi->x_old, x);
            switch (cgi->method) {
            case TURNING_POINT_CONTINUATION:
              vec_copy(x_tang_old, cpi->x_tang);
              assign_bif_parameter_conwrap(bif_param_old);
              con->turning_point_info.bif_param = bif_param_old;
              break;
            case PITCHFORK_CONTINUATION:
              vec_copy(x_tang_old, cpi->x_tang);
              assign_bif_parameter_conwrap(bif_param_old);
              con->pitchfork_info.bif_param = bif_param_old;
              break;
            case PHASE_TRANSITION_CONTINUATION:
              vec_copy(x2_old, con->phase_transition_info.x2);
              assign_bif_parameter_conwrap(bif_param_old);
              con->phase_transition_info.bif_param = bif_param_old;
              break;
            case HOPF_CONTINUATION:
              vec_copy(y_vec_old, con->hopf_info.y_vec);
              vec_copy(z_vec_old, con->hopf_info.z_vec);
              assign_bif_parameter_conwrap(bif_param_old);
              con->hopf_info.bif_param = bif_param_old;
              con->hopf_info.omega = omega_old;
              break;
            }
            break;
          case 1:
            for (i = 0; i < cgi->numUnks; i++)
              x[i] = cpi->x_old[i] - step * cpi->x_tang[i];
            break;
          case 2:
            cpi->arc_step *= 0.5;
            for (i = 0; i < cgi->numUnks; i++)
              x[i] = cpi->x_old[i] + cpi->arc_step * cpi->x_tang[i];
            break;
          }

          /*
           * Also, reset last_step flag to continue to final parameter value.
           */

          csi->last_step = FALSE;
        }

        /*
         * If it was the last step, however, reset to previous solution.
         */

        else {
          cgi->param = cpi->param_old;
          assign_parameter_conwrap(cgi->param);
          step = 0.0;
          vec_copy(x, cpi->x_old);
          csi->last_step = TRUE;
        }

        /*
         * Print out failure message
         */

        if (cgi->printproc > 1)
          print_cont_step_fail(order, step, con);
      }
    }

    else {

      /*
       * Solver did Converge!!
       */

      /*
       * Check to see if parameter value passed end value (either direction)
       */

      end_passed = (cgi->param - csi->max_param) * (cpi->param_old - csi->max_param);
      if ((cpi->step_num != 0 && end_passed <= 0.0) || n == csi->max_steps)
        csi->last_step = TRUE;

      /*
       * For arc-length continuation using the tangent factor,
       * calculate cpi->x_tang and new tang_factor in advance,
       * then compare to value for previous step.  If this
       * difference exceeds macpi->x_tang_step, treat this as a failed
       * step and halve arc_step as above.
       */

      if (!csi->last_step && cgi->method == ARC_LENGTH_CONTINUATION) {
        if (cgi->printproc > 4) {
          printf("\n\tDoing Pseudo Arc-length continuation --");
          printf("\n\tCalculating tangent vector by one linear "
                 "solve\n");
        }

        if (cpi->step_num == 0)
          step = csi->first_step;

        calc_rhs_continuation(CONT_TANGENT, x, cpi->x_tang, NULL, NULL, NULL, cgi->param,
                              cgi->perturb, NULL, cgi->numUnks, cgi->numOwnedUnks);

        i = linear_solver_conwrap(cpi->x_tang, NEW_JACOBIAN, NULL);

        if (cpi->step_num > 0)
          tang_factor =
              fabs(scaled_dot_prod(x_tang_old, cpi->x_tang, cpi->scale_vec, cgi->numOwnedUnks)) /
              sqrt(scaled_dot_prod(x_tang_old, x_tang_old, cpi->scale_vec, cgi->numOwnedUnks) *
                   scaled_dot_prod(cpi->x_tang, cpi->x_tang, cpi->scale_vec, cgi->numOwnedUnks));
        tang_factor = exp(cai->tang_exp * log(tang_factor));
        if (cgi->printproc > 7)
          printf(" Tangent factor is %9g\n", tang_factor);

        /* Repeat step if tang_factor is too small */

        if (cpi->step_num > 1 && tang_factor < cai->tang_step_limit) {

          if (cgi->printproc > 7)
            printf(" Step limit exceeded: Retrying ...");
          tan_flag = FALSE;
          step *= 0.5;
          cgi->param = cpi->param_old + step;
          assign_parameter_conwrap(cgi->param);
          cpi->arc_step *= 0.5;
          for (i = 0; i < cgi->numUnks; i++)
            x[i] = cpi->x_old[i] + cpi->arc_step * x_tang_old[i];
          vec_copy(x_tang_old, cpi->x_tang);
        }
      }

      /*
       * If tan_flag has not been set to zero, proceed with continuation
       */

      if (tan_flag == TRUE) {

        /*
         * Print out final results of a successful time step
         */

        if (order == 2)
          step = cgi->param - cpi->param_old;

        if (cgi->printproc > 4)
          print_cont_step2(order, step, con);

        /*
         * If first continuation step, set to value from input file.  If
         * controlled steps, use # of Newton iters to pick next step size
         * For arc-length continuation, impose maximum step size as
         * approximated by arc_step and dp_ds.
         * Note:  without time step control, step size can never increase.
         */

        step_old = step;
        if (cpi->step_num == 0)
          step = csi->first_step;
        else {

          /* normal step control */
          if (!csi->last_step && csi->step_ctrl > 0.0) {
            if (order == 2)
              cpi->arc_step *=
                  simple_step_control(num_newt_conv, csi->max_newton_its, csi->step_ctrl);
            else {
              step *= simple_step_control(num_newt_conv, csi->max_newton_its, csi->step_ctrl);
              if (fabs(step) > csi->max_delta_p)
                step = ((step > 0) ? csi->max_delta_p : -csi->max_delta_p);
            }
          }
          /* for constant step runs where the step has been cut, let it
           * increase again with step control of 0.5
           */
          else if (order < 2 && fabs(step) < fabs(csi->first_step)) {
            step *= simple_step_control(num_newt_conv, csi->max_newton_its, 0.5);
          } else if (order == 2 && fabs(cpi->arc_step) < arc_step_adj) {
            cpi->arc_step *= simple_step_control(num_newt_conv, csi->max_newton_its, 0.5);
            if (cpi->arc_step > arc_step_adj)
              cpi->arc_step = arc_step_adj;
          }
        }

        delta_s = ((cgi->method == ARC_LENGTH_CONTINUATION) ? cpi->arc_step : step);

        /*
         * Output information at the end of every successful time step
         * Depending on the solution method, there can be up to three
         * solution vectors and parameter values to write out at a solution
         */

        switch (cgi->method) {
        case ZERO_ORDER_CONTINUATION:
        case FIRST_ORDER_CONTINUATION:
        case ARC_LENGTH_CONTINUATION:
          solution_output_conwrap(1, x, cgi->param, NULL, 0.0, NULL, 0.0, cpi->step_num,
                                  num_newt_conv, con);
          break;
        case TURNING_POINT_CONTINUATION:
          solution_output_conwrap(2, x, cgi->param, cpi->x_tang, con->turning_point_info.bif_param,
                                  NULL, 0.0, cpi->step_num, num_newt_conv, con);
          break;
        case PITCHFORK_CONTINUATION:
          solution_output_conwrap(2, x, cgi->param, cpi->x_tang, con->pitchfork_info.bif_param,
                                  NULL, 0.0, cpi->step_num, num_newt_conv, con);
          break;
        case HOPF_CONTINUATION:
          solution_output_conwrap(3, x, cgi->param, con->hopf_info.y_vec, con->hopf_info.bif_param,
                                  con->hopf_info.z_vec, con->hopf_info.omega, cpi->step_num,
                                  num_newt_conv, con);
          break;
        case PHASE_TRANSITION_CONTINUATION:
          solution_output_conwrap(2, x, cgi->param, con->phase_transition_info.x2,
                                  con->phase_transition_info.bif_param, NULL, 0.0, cpi->step_num,
                                  num_newt_conv, con);
          break;
        }

        /* Check element quality */
        int good_mesh = element_quality(exo, x, ams->proc_config);
        if (!good_mesh) {
          goto free_and_clear;
        }

        /*
         * Check current parameter value against the maximum.
         */

        if (csi->last_step) {
          n = csi->max_steps; /* Force IO and exit  */
          if (cgi->printproc > 1) {
            printf("Simulation completed continuation in %d steps\n", cpi->nstep);
            printf("Final Parameter Value: %e\n\n", cgi->param);
            if (end_passed <= 0.) {
              printf("\n\n\t I will continue no more!\n\t No more continuation for you!\n");
            }
          }
        }

        if (n < csi->max_steps) {

          /*
           * Finally, its time to do some continuation, since the previous step
           * converged and it wasn't the last step.
           */

          if (cgi->printproc > 4) {
            printf("\n");
            print_line("~", 80);
            printf("\nCalculating initial guess for next continuation "
                   "step\n");
          }

          /*
           * Set the continuation order 0, 1 (Euler-Newton), or 2 (Arc-length)
           * Always do only zero order continuation of a turning point
           * or any other special continuation.
           */

          switch (cgi->method) {
          case FIRST_ORDER_CONTINUATION:
            order = 1;
            break;
          case ARC_LENGTH_CONTINUATION:
            order = 2;
            break;
          default:
            order = 0;
          }

          /*
           * Possibly adjust the step value for this step so it hits maximum
           * value exactly (zero or first order continuation).
           * For arc length continuation, this is done after solution scaling.
           */

          if (order != 2 && cpi->step_num != 0) {
            temp_step = delta_s;
            end_passed =
                ((cgi->param + temp_step - csi->max_param) * (cgi->param - csi->max_param));

            /*
             * If end_passed <= 0, next step would take param past end value
             */

            if (end_passed <= 0) {
              step = csi->max_param - cgi->param;
              csi->last_step = TRUE;
              if (cgi->printproc > 7)
                fprintf(stderr, "\n\t ******** LAST PATH STEP!\n");
            }
          }

          /*
           * Calculate the tangent to the solution, cpi->x_tang, for the
           * current step. This is trivial for 0-order continuation, requires
           * 1 linear solve for 1st order continuation and for arc-length
           * continuation. Use tangent to predict new solution in x.
           */

          switch (order) {

          case 0:
            if (cgi->printproc > 4) {
              printf("\n   Doing Zeroth-order continuation --");
              printf("\n   previous solution used as initial guess\n");
            }

            /*
             * NO definition of cpi->x_tang needed for zero order.
             * Don't set it to zero because that will mess up the
             * turning point and phase transition tracking algorithms,
             * which use that space for other purposes.
             */

            /*
             * Save the old solution, before overwriting with the new solution
             * for use in restarting after failed steps
             */

            vec_copy(x, cpi->x_old);
            switch (cgi->method) {
            case TURNING_POINT_CONTINUATION:
              vec_copy(cpi->x_tang, x_tang_old);
              bif_param_old = con->turning_point_info.bif_param;
              break;
            case PITCHFORK_CONTINUATION:
              vec_copy(cpi->x_tang, x_tang_old);
              bif_param_old = con->pitchfork_info.bif_param;
              break;
            case PHASE_TRANSITION_CONTINUATION:
              vec_copy(con->phase_transition_info.x2, x2_old);
              bif_param_old = con->phase_transition_info.bif_param;
              break;
            case HOPF_CONTINUATION:
              vec_copy(con->hopf_info.y_vec, y_vec_old);
              vec_copy(con->hopf_info.z_vec, z_vec_old);
              bif_param_old = con->hopf_info.bif_param;
              omega_old = con->hopf_info.omega;
              break;
            }

            /* perturb guess off of singularity in tp_continuation */

            if (cgi->method == TURNING_POINT_CONTINUATION ||
                cgi->method == PITCHFORK_CONTINUATION) {
              if (cgi->printproc > 4)
                printf("\tcon_lib: Adding random"
                       " perturbation for continuation\n");
              perturb_solution_conwrap(x, cpi->x_old, cpi->scale_vec, cgi->numOwnedUnks);
            }

            break;

          case 1:
            if (cgi->printproc > 4) {
              printf("\n   Doing First-order continuation --");
              printf("\n   calculating tangent vector by one linear "
                     "solve\n");
            }

            /*
             * Choose perturbation for numerical derivative of Residuals w.r.t
             * continuation parameter, and solve for the tangent vector as in
             * eq. 7.13 in John's thesis.  The continuation parameter and
             * perturbation, cgi->param and delta_param, are passed to the
             * linear solver in the spots for the time and CJ, which aren't'
             * needed for steady problems.
             */

            calc_rhs_continuation(CONT_TANGENT, x, cpi->x_tang, NULL, NULL, NULL, cgi->param,
                                  cgi->perturb, NULL, cgi->numUnks, cgi->numOwnedUnks);

            linear_solver_conwrap(cpi->x_tang, NEW_JACOBIAN, NULL);

            /*
             * Save the old solution, before overwriting with the new solution
             */

            vec_copy(x, cpi->x_old);

            /*
             * Multiply the tangent vector, cpi->x_tang initially, by the step
             * length, and add to x, to obtain an initial guess at
             * the next parameter value.
             */

            for (i = 0; i < cgi->numUnks; i++) {
              x[i] -= cpi->x_tang[i] * step;
            }

            break;

          case 2: /* Arclength ccontinuation predictor and step control */

            /* cpi->x_tang vector found above. */

            /* Function "solution_scale" rescales solution as necessary. */

            ds_ratio = solution_scale(con, arc);
            if (cgi->printproc > 7) {
              printf(" Solution scale factor is %9g\n", arc->dx_fac);
            }

            /* Adjust arc_step for current scale factor. */
            cpi->arc_step /= ds_ratio;
            arc_step_adj = fabs(first_arc_step * arc->ds_fac);

            /* Adjust arc_step for current tangent factor. */
            cpi->arc_step *= tang_factor;

            /* Also reduce arc_step if too large. */
            max_step = fabs(csi->max_delta_p / cpi->dp_ds);
            if (cpi->arc_step > max_step)
              cpi->arc_step = max_step;

            /*
             * Readjust the step value for this step so it hits maximum
             * or end value (approximately) for arc length continuation.
             */

            if (cpi->step_num != 0) {
              temp_step = delta_s * cpi->dp_ds;
              if (step < 0)
                temp_step *= -1.0;
              end_passed =
                  ((cgi->param + temp_step - csi->max_param) * (cgi->param - csi->max_param));

              /* If end_passed < 0, next step would take param past end value */

              if (end_passed < 0) {
                temp_step = (csi->max_param - cgi->param);
                if (step < 0)
                  temp_step *= -1.0;
                cpi->arc_step = fabs(temp_step / cpi->dp_ds);
                csi->last_step = TRUE;
                if (cgi->printproc > 7)
                  fprintf(stderr, "\n\t ******** LAST PATH STEP!\n");
              }
            }

            /*
             * If this is the first step, pick cpi->arc_step so that this step
             * will progress the parameter by approximately csi->step
             */

            if (cpi->step_num == 0) {

              if (step < 0) {
                cpi->dp_ds *= -1.0;
                sn_old = -1;
              } else
                sn_old = 1;

              cpi->arc_step = step / cpi->dp_ds;
              first_arc_step = cpi->arc_step;
            } else {

              /*
               * Pick sign of cpi->dp_ds according to eq. 7.14b in JNS thesis
               * NOTE: -1.0 factor multiplying solution terms is because
               * cpi->x_tang is currently the negative of the tangent vector.
               *
               * and check if a turning point was passed --
               */

              if (-1.0 * (scaled_dot_prod(cpi->x_tang, x, cpi->scale_vec, cgi->numOwnedUnks) -
                          scaled_dot_prod(cpi->x_tang, cpi->x_old, cpi->scale_vec,
                                          cgi->numOwnedUnks)) +
                      cgi->param - cpi->param_old <
                  0.0) {

                cpi->dp_ds *= -1.0;
                sn_new = -1;
              } else
                sn_new = 1;

              if ((cgi->printproc > 1) && sn_old != sn_new)
                printf("\n\n\tA turning point was passed !!!!!!!\n");

              sn_old = sn_new;
            }

            /*
             * Save the old solution, before overwriting with the new solution
             */
            vec_copy(x, cpi->x_old);

            /*
             * Calculate prediction for next step from Eqs. 7.15&7.16 in JNS
             * thesis (leaving cpi->x_tang = u_dot).
             */

            for (i = 0; i < cgi->numUnks; i++) {
              cpi->x_tang[i] *= -cpi->dp_ds;
              x[i] += cpi->x_tang[i] * cpi->arc_step;
            }
            step = cpi->dp_ds * cpi->arc_step;
            vec_copy(cpi->x_tang, x_tang_old);

            break;
          }

          /*
           * Increment the continuation parameter.  Update the
           * BCs/Properties/whatever that the continuation parameter really
           * represents.
           */

          cpi->param_old = cgi->param;
          cgi->param += step;
          assign_parameter_conwrap(cgi->param);

          /*
           * Increment the step counter. Print final message.
           */

          cpi->step_num++;
          cpi->nstep++;

        } /* END of:  if (n < csi->max_steps) */

      } /* END of:  if (tan_flag == TRUE) */

    } /* END of else section for converged solves */

  } /* END of loop over continuation step attempts --- for (n = 0; ... --- */

  /*********************CLEAN-UP AREA*****************************************/
free_and_clear:
  /*
   * Free auxillary vectors no matter what happened
   */

  free_vec(&cpi->x_old);
  free_vec(&cpi->scale_vec);

  switch (cgi->method) {
  case ZERO_ORDER_CONTINUATION:
  case FIRST_ORDER_CONTINUATION:
    free_vec(&cpi->x_tang);
    break;
  case ARC_LENGTH_CONTINUATION:
  case TURNING_POINT_CONTINUATION:
  case PITCHFORK_CONTINUATION:
    if (!cgi->nv_save)
      free_vec(&cpi->x_tang);
    free_vec(&x_tang_old);
    break;
  case PHASE_TRANSITION_CONTINUATION:
    free_vec(&x2_old);
    break;
  case HOPF_CONTINUATION:
    free_vec(&y_vec_old);
    free_vec(&z_vec_old);
    break;
  }

  /*
   * Send back the overall result of the time step
   */

  return cpi->step_num;

} /**************** END of solve_continuation () *****************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double solution_scale(struct con_struct *con, struct arc_scale_struct *arc) {

  double ds_ratio = 1.0;
  double umag;
  int i;

  /* shorthand to some commonly used variables */
  struct general_info_struct *cgi = &(con->general_info);
  struct private_info_struct *cpi = &(con->private_info);
  struct arclength_info_struct *cai = &(con->arclength_info);

  /* Calculate average of each variable for scaling of arc-length eq */

  calc_scale_vec_conwrap(cgi->x, cpi->scale_vec, cgi->numUnks);

  /* Get dot product before solution scaling. */

  umag = arc->dx0 * arc->dx0 *
         scaled_dot_prod(cpi->x_tang, cpi->x_tang, cpi->scale_vec, cgi->numOwnedUnks);

  /* Adjust cpi->scale_vec by current scale factor arc->dx_fac. */

  for (i = 0; i < cgi->numUnks; i++)
    cpi->scale_vec[i] *= arc->dx_fac;

  arc->umag2 = scaled_dot_prod(cpi->x_tang, cpi->x_tang, cpi->scale_vec, cgi->numOwnedUnks);

  /*
   * Calculate deriv of parameter w.r.t arc length, Eq.7.14a in JNS
   * thesis.  Save actual value in arc->dp_ds_act.
   */

  cpi->dp_ds = 1.0 / sqrt(1.0 + arc->umag2);
  arc->dp_ds_old = cpi->dp_ds;

  /*
   * On the first step, set arc->dx_fac to the value which will correspond
   * to the desired dp_ds value.
   */

  if (cpi->step_num == 0 && cai->dp_ds2_goal > 1e-6) {
    arc->dx_fac = cpi->dp_ds * sqrt((1.0 - cai->dp_ds2_goal) / (1.0 - cpi->dp_ds * cpi->dp_ds)) /
                  arc->dp_ds_goal;
    arc->dx0 = arc->dx_fac;
    for (i = 0; i < cgi->numUnks; i++)
      cpi->scale_vec[i] *= arc->dx_fac;
    arc->umag2 = scaled_dot_prod(cpi->x_tang, cpi->x_tang, cpi->scale_vec, cgi->numOwnedUnks);
    cpi->dp_ds = 1.0 / sqrt(1.0 + arc->umag2);
    ds_ratio = cpi->dp_ds / arc->dp_ds_old;
  }

  /*
   * If dp_ds is too large, calculate factor to adjust arc_step
   * such that (dp)^2 / (ds)^2 doesn't exceed dp_ds2_goal.
   */

  if (cpi->step_num > 0 && cai->dp_ds2_goal > 1e-6 && cpi->dp_ds > cai->dp_ds_max) {
    if (cgi->printproc > 7)
      printf("  dp_ds out of limits at %g: Rescaling...\n", cpi->dp_ds);

    /* Calculate scale factor for cpi->scale_vec (arc->dx_fac). */
    /* Limit arc->dx_fac to arc->dx_fac_max to avoid division by zero. */

    arc->dx_fac_old = arc->dx_fac;
    if (cpi->dp_ds > arc->dp_ds_limit)
      arc->dx_fac = arc->dx_fac_max;

    else
      arc->dx_fac *= cpi->dp_ds * sqrt((1.0 - cai->dp_ds2_goal) / (1.0 - cpi->dp_ds * cpi->dp_ds)) /
                     arc->dp_ds_goal;

    /* Multiply cpi->scale_vec through by arc->dx_fac. */

    for (i = 0; i < cgi->numUnks; i++)
      cpi->scale_vec[i] *= arc->dx_fac / arc->dx_fac_old;

    /* Recalculate unknown dot product (arc->umag2) and dp_ds. */

    arc->umag2 = scaled_dot_prod(cpi->x_tang, cpi->x_tang, cpi->scale_vec, cgi->numOwnedUnks);
    cpi->dp_ds = 1.0 / sqrt(1.0 + arc->umag2);
    ds_ratio = cpi->dp_ds / arc->dp_ds_old;
  }

  /* Get arc->ds_fac */

  arc->ds_fac = 1.0 / (cpi->dp_ds * sqrt(1.0 + umag));

  /* Return ds_ratio */

  return ds_ratio;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double scaled_dot_prod(double *x, double *y, double *scale_vec, int n) {
  int i;
  double sum = 0;

  for (i = 0; i < n; i++)
    sum += x[i] * y[i] * scale_vec[i] * scale_vec[i];
  return gsum_double_conwrap(sum);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
static double simple_step_control(int num_newt_conv, int max_Newton_steps, double step_ctrl)
/*
 * Function to calculate the increase in time step for the pseudo time-
 * integration option based on:
 *
 *    num_newt_conv,      the number of Newton iterations the last step
 *                        required to reach convergence, and
 *    max_Newton_steps,   the maximum number of Newton steps allowed.
 *    step_ctrl           aggressiveness of step size routine,
 *                        0.0 for constant step, 2.0 is very big
 *
 * This simplistic function will increase the time step if the last step
 * converged. It is a quadratic function in the ratio of the number of
 * Newton steps taken compared to the maximum number of steps allowed
 * up to a maximum of 1+aggressiveness times the previous step.
 */

{

  double factor;

  factor = (max_Newton_steps - (double)num_newt_conv) / (max_Newton_steps - 1.0);

  return (1.0 + step_ctrl * factor * factor);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
static void print_cont_step1(int order, double step, double step_old, struct con_struct *con)

/*
 * Print out for relevant time step information
 */

{
  /* shorthand to some commonly used variables */
  struct general_info_struct *cgi = &(con->general_info);
  struct private_info_struct *cpi = &(con->private_info);

  const char *string = "";

  if (cgi->method == TURNING_POINT_CONTINUATION)
    string = "Zero-order Turning Point Continuation";
  else if (cgi->method == PHASE_TRANSITION_CONTINUATION)
    string = "Zero-order Phase Transition Continuation";
  else if (cgi->method == PITCHFORK_CONTINUATION)
    string = "Zero-order Pitchfork Continuation";
  else if (cgi->method == HOPF_CONTINUATION)
    string = "Zero-order Hopf Continuation";
  else if (order == 0)
    string = "Zero-order Continuation";
  else if (order == 1)
    string = "First-order Continuation";
  else if (order == 2)
    string = "Pseudo Arc-length Continuation";

  printf("\n");
  print_line("~", 80);
  printf("\nStart of Step: %5d   Continuation Param = %9.5g from "
         "%9.5g\n",
         cpi->nstep, cgi->param, cgi->param - step);
  printf("\tContinuation method = %s\n", string);
  printf("\tdelta_c_p        = %8.5e ", step);
  printf("\n\tdelta_c_p_old    = %8.5e\n", step_old);

  if (order == 2)
    printf("\tdelta_s = %g\n", cpi->arc_step);
  printf("\n");
  print_line("~", 80);

} /*************** END print_cont_step1 **************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void print_cont_step2(int order, double step, struct con_struct *con)

/*
 * Print out for relevant continuation step information
 */

{
  /* shorthand to some commonly used variables */
  struct general_info_struct *cgi = &(con->general_info);
  struct private_info_struct *cpi = &(con->private_info);

  printf("\n");
  print_line("~", 80);
  printf("Continuation Step Number %5d was a success: Param = %10g\n", cpi->nstep, cgi->param);
  if (cgi->method == TURNING_POINT_CONTINUATION)
    printf("\tTurning Point located at:  %10g  %10g\n", cgi->param,
           con->turning_point_info.bif_param);
  else if (cgi->method == PHASE_TRANSITION_CONTINUATION)
    printf("\tPhase Transition located at:  %10g  %10g\n", cgi->param,
           con->phase_transition_info.bif_param);
  else if (cgi->method == PITCHFORK_CONTINUATION)
    printf("\tPitchfork Bifurcation located at:  %10g  %10g\n", cgi->param,
           con->pitchfork_info.bif_param);
  else if (cgi->method == HOPF_CONTINUATION)
    printf("\tHopf Bifurcation located at:  %10g  %10g  %10g\n", cgi->param,
           con->hopf_info.bif_param, con->hopf_info.omega);
  else
    printf("\tOrder   = %d\n", order);

  if (order == 2)
    printf("\t delta_s = %g", cpi->arc_step);
  printf("\tdelta_c_p = %g\n", step);
  printf("\n");
  print_line("~", 80);
  printf("\n");

} /************* END of print_cont_step2 () **********************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void print_cont_step_fail(int order, double step, struct con_struct *con)

/*
 * Print Out descriptive information on why the current step failed
 */

{
  /* shorthand to some commonly used variables */
  struct general_info_struct *cgi = &(con->general_info);
  struct private_info_struct *cpi = &(con->private_info);

  printf("\n");
  print_line("~", 80);
  printf("\tContinuation Step Number %5d experienced a convergence "
         "failure\n",
         cpi->nstep);
  printf("\tin the non-linear or linear solver\n");
  printf("\t\tValue of parameter at failed step    = %g\n", cgi->param);

  if (order < 2) {
    printf("\t\tdelta_c_p of the failed step         = %g\n", 2.0 * step);
    printf("\t\tNext value of delta_c_p              = %g\n", step);
  } else {
    printf("\t\tdelta_s of the failed step         = %g\n", 2.0 * cpi->arc_step);
    printf("\t\tNext value of delta_s              = %g\n", cpi->arc_step);
  }
  printf("\n");
  print_line("~", 80);

} /*************** END of print_cont_step_fail () ****************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
static void print_line(char *charstr, int ntimes)

{
  int i;
  for (i = 0; i < ntimes - 1; i++)
    printf("%c", *charstr);
  printf("\n");
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
