/************************************************************************ *
 * Goma - Multiphysics finite element software                             *
 * Sandia National Laboratories                                            *
 *                                                                         *
 * Copyright (c) 2014 Sandia Corporation.                                  *
 *                                                                         *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,  *
 * the U.S. Government retains certain rights in this software.            *
 *                                                                         *
 * This software is distributed under the GNU General Public License.	  *
 \************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "std.h"

#include "exo_struct.h"
#include "rf_fem_const.h"
#include "rf_vars_const.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "rf_node_const.h"
#include "usr_print.h"
#include "sl_amesos_interface.h"
#include "brk_utils.h"

#include "sl_epetra_interface.h"
#include "sl_epetra_util.h"

#define _RF_SOLVE_SEGREGATED_C
#include "goma.h"
#include "mm_solve_linear_segregated.h"
#include "el_quality.h"

#define ROUND_TO_ONE 0.9999999
/*************************************************************************************
 *  solve_problem_segregated
 *
 *      	 Routine that controls solution of overall FEM problem in segregated fashion.
 *        This routine will control time integration, order of matrices to be solved,
 *        matrices fills and nonlinear solver and results output.
 *
 *        Started as a clone from solve_problem routine
 *
 * Revision History
 * ================
 * 26 November 2014 - Kristianto Tjiptowidjojo - Creation.
 *
 * ******************************************************************************/

void solve_problem_segregated(Exo_DB *exo, /* ptr to the finite element mesh database  */
Dpi *dpi, /* distributed processing information       */
dbl *te_out) /* te_out - return actual end time */
{

  /*
   * Sparse matrix storage vectors
   * (MSR format.  See "SPARSKIT: a basic tool kit for sparse matrix
   * computations" by Youcef Saad)
   */
  int **ija = NULL; /* column pointer array                     */
  double **a = NULL; /* nonzero array                            */
  double **a_old = NULL; /* nonzero array                            */
  int **ija_attic = NULL; /* storage for external dofs                */

  /*
   * Variables
   */

  static double **x = NULL; /* solution vector                   */

  static double **x_old = NULL; /* old solution vector               */
  static double **x_older = NULL; /* older solution vector             */
  static double **x_oldest = NULL; /* oldest solution vector saved      */
  static double **xdot = NULL; /* current time derivative of soln   */
  static double **xdot_old = NULL; /* old time derivative of soln       */
  static double **xdot_older = NULL; /* old time derivative of soln       */
  double **x_pred = NULL;

  double **x_update = NULL; /* update at last iteration          */
  double **resid_vector = NULL; /* residual                          */
  double **resid_vector_sens = NULL; /* residual sensitivity              */
  double **scale = NULL; /* scale vector for modified newton  */

  int *node_to_fill = NULL;

  static struct Aztec_Linear_Solver_System **ams;

  /*
   * Variables for time integration
   */
  double *timeValueRead;
  double timeValueReadTrans = 0.0;
  double time_print, i_print;
  double theta = 0.0, time;
  static double time1 = 0.0; /* Current time that the simulation is trying  to find the solution for */
  double delta_t_new, delta_t, delta_t_old, delta_t_older, delta_t_oldest;
  char tspstring[MAX_FNL];
  int n;
  int nt;
  int last_renorm_nt; /* time step at which last renorm occured   */
  int time_step_reform; /* counter for jacobian reformation stride  */
  int converged = TRUE; /* success or failure of Newton iteration   */
  static int nprint = 0;
  int step_fix = 0;
  int const_delta_t, const_delta_ts, step_print;
  int success_dt;
  int failed_recently_countdown = 0;

  /*
   * Local variables
   */
  int good_mesh;
  int error, err;
  static double **gvec = NULL;
  static double ****gvec_elem = NULL;
  static struct Results_Description **rd;

  double *gv = NULL; /* Global variable values for ExoII database */
  int *tnv; /* total number of nodal variables and kinds */
  int *tev; /* total number of elem variables and kinds  */
  int *tnv_post; /* total number of nodal variables and kinds
   for post processing                       */
  int *tev_post; /* total number of elem variables and kinds
   for post processing                       */

  unsigned int matrix_systems_mask;

  int i;
  int inewton;
  int *numProcUnknowns;

#ifdef RELAX_ON_TRANSIENT_PLEASE
  int relax_bit = TRUE; /* Enables relaxation after a transient convergence failure*/
#else
  int relax_bit = FALSE;
#endif

  static int callnum = 1; /* solve_problem_segregated call counter */

  static const char yo[] = "solve_problem_segregated"; /* So my name is in a string.        */

  timeValueRead = calloc(upd->Total_Num_Matrices, sizeof(double));
  tnv = malloc(sizeof(int) * upd->Total_Num_Matrices);
  tev = malloc(sizeof(int) * upd->Total_Num_Matrices);
  tnv_post = malloc(sizeof(int) * upd->Total_Num_Matrices);
  tev_post = malloc(sizeof(int) * upd->Total_Num_Matrices);
  gvec_elem = malloc(sizeof(double ***) * upd->Total_Num_Matrices);

  if (Num_Proc > 1 && tran->fix_freq > 0) {
    step_fix = 1; /* Always fix on the first timestep to match print frequency */
  }

  /*
   *  Malloc the space for the results description structure and set all
   *  of that space to naught initially.
   *  Do this only once if in library mode.
   */
  if (callnum == 1) {
    rd = malloc(sizeof(struct Results_Description *) * upd->Total_Num_Matrices);
    for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
      rd[pg->imtrx] = alloc_struct_1(struct Results_Description, 1);
    }
  }

  for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
    /*
     * Some preliminaries to help setup EXODUS II database output.
     */

    tnv[pg->imtrx] = cnt_nodal_vars(); /*  tnv_post is calculated in load_nodal_tkn*/
    tev[pg->imtrx] = cnt_elem_vars(); /*  tev_post is calculated in load_elem_tkn*/

    if (tnv[pg->imtrx] < 0) {
      DPRINTF(stderr, "%s:\tbad tnv.\n", yo);
      EH(-1, "\t");
    }

    if (tev[pg->imtrx] < 0) {
      DPRINTF(stderr, "%s:\tMaybe bad tev? See goma design committee ;) \n",
          yo);
      EH(-1, "\t");
    }

    rd[pg->imtrx]->nev = 0; /* number element variables in results */
    rd[pg->imtrx]->ngv = 0; /* number global variables in results  */
    rd[pg->imtrx]->nhv = 0; /* number history variables in results */

    rd[pg->imtrx]->ngv = 5; /* number global variables in results
     * see load_global_var_info for names
     */

    error = load_global_var_info(rd[pg->imtrx], 0, "CONV");
    error = load_global_var_info(rd[pg->imtrx], 1, "NEWT_IT");
    error = load_global_var_info(rd[pg->imtrx], 2, "MAX_IT");
    error = load_global_var_info(rd[pg->imtrx], 3, "CONVRATE");
    error = load_global_var_info(rd[pg->imtrx], 4, "MESH_VOLUME");

    gv = alloc_dbl_1(rd[pg->imtrx]->ngv, 0.0);

    /*
     *  Load output nodal types, kinds and names into the structure
     *  which will be used to define what's in the output file.
     */
    error = load_nodal_tkn(rd[pg->imtrx], &tnv[pg->imtrx],
        &tnv_post[pg->imtrx]);
    if (error != 0) {
      DPRINTF(stderr, "%s:  problem with load_nodal_tkn()\n", yo);
      EH(-1, "\t");
    }

    /*
     *  Load output element var types, kinds and names into the structure
     *  which will be used to define what's in the output file.
     */
    error = load_elem_tkn(rd[pg->imtrx], exo, tev[pg->imtrx],
        &tev_post[pg->imtrx]);
    if (error != 0) {
      DPRINTF(stderr, "%s:  problem with load_elem_tkn()\n", yo);
      EH(-1, "\t");
    }

    /*
     * Write out the names of the nodal variables that we will be sending to
     * the EXODUS II output file later - do only once if in library mode.
     */

    if (callnum == 1) {
      gvec_elem[pg->imtrx] = (double ***) alloc_ptr_1(exo->num_elem_blocks);
      if ((tev[pg->imtrx] + tev_post[pg->imtrx]) > 0) {
        for (i = 0; i < exo->num_elem_blocks; i++) {
          gvec_elem[pg->imtrx][i] = (double **) alloc_ptr_1(
              tev[pg->imtrx] + tev_post[pg->imtrx]);
        }
      }
    }
  }

  if (callnum == 1) {
    wr_result_prelim_exo_segregated(rd, exo, ExoFileOut, gvec_elem);
  }

  /*
   * This gvec workhorse transports output variables as nodal based vectors
   * that are gather from the solution vector. Note: it is NOT a global
   * vector at all and only carries this processor's nodal variables to
   * the exodus database.
   */

  gvec = malloc(sizeof(double *) * upd->Total_Num_Matrices);
  for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
    asdv(&gvec[pg->imtrx], Num_Node);
  }

  /*
   * Allocate space and manipulate for all the nodes that this processor
   * is aware of...
   */

  numProcUnknowns = malloc(upd->Total_Num_Matrices * sizeof(int));
  for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
    numProcUnknowns[pg->imtrx] = NumUnknowns[pg->imtrx]
        + NumExtUnknowns[pg->imtrx];
  }

  resid_vector = malloc(upd->Total_Num_Matrices * sizeof(double *));
  resid_vector_sens = malloc(upd->Total_Num_Matrices * sizeof(double *));
  scale = malloc(upd->Total_Num_Matrices * sizeof(double *));

  for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
    asdv(&(resid_vector[pg->imtrx]), numProcUnknowns[pg->imtrx]);
    asdv(&(resid_vector_sens[pg->imtrx]), numProcUnknowns[pg->imtrx]);
    asdv(&(scale[pg->imtrx]), numProcUnknowns[pg->imtrx]);
  }

  /*
   * Allocate Aztec structures and initialize all elements to zero
   */
  if (callnum == 1) {
    ams = (struct Aztec_Linear_Solver_System **) alloc_ptr_1(
        upd->Total_Num_Matrices);
    for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
      ams[pg->imtrx] = alloc_struct_1(struct Aztec_Linear_Solver_System, 1);
    }
  }

  for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
    AZ_set_proc_config(ams[pg->imtrx]->proc_config, MPI_COMM_WORLD);
  }

  /* Allocate solution arrays on first call only */
  if (callnum == 1) {
    x = malloc(upd->Total_Num_Matrices * sizeof(double *));
    x_old = malloc(upd->Total_Num_Matrices * sizeof(double *));
    x_older = malloc(upd->Total_Num_Matrices * sizeof(double *));
    x_oldest = malloc(upd->Total_Num_Matrices * sizeof(double *));
    xdot = malloc(upd->Total_Num_Matrices * sizeof(double *));
    xdot_old = malloc(upd->Total_Num_Matrices * sizeof(double *));
    xdot_older = malloc(upd->Total_Num_Matrices * sizeof(double *));
    x_pred = malloc(upd->Total_Num_Matrices * sizeof(double *));

    for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
      x[pg->imtrx] = alloc_dbl_1(numProcUnknowns[pg->imtrx], 0.0);
      x_old[pg->imtrx] = alloc_dbl_1(numProcUnknowns[pg->imtrx], 0.0);
      x_older[pg->imtrx] = alloc_dbl_1(numProcUnknowns[pg->imtrx], 0.0);
      x_oldest[pg->imtrx] = alloc_dbl_1(numProcUnknowns[pg->imtrx], 0.0);
      xdot[pg->imtrx] = alloc_dbl_1(numProcUnknowns[pg->imtrx], 0.0);
      xdot_old[pg->imtrx] = alloc_dbl_1(numProcUnknowns[pg->imtrx], 0.0);
      xdot_older[pg->imtrx] = alloc_dbl_1(numProcUnknowns[pg->imtrx], 0.0);
      x_pred[pg->imtrx] = alloc_dbl_1(numProcUnknowns[pg->imtrx], 0.0);
    }
  }

  x_update = malloc(upd->Total_Num_Matrices * sizeof(double *));
  for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
    x_update[pg->imtrx] = alloc_dbl_1(
        numProcUnknowns[pg->imtrx] + numProcUnknowns[pg->imtrx], 0.0);
  }

  /* Allocate sparse matrix */

  ija = malloc(upd->Total_Num_Matrices * sizeof(int *));
  ija_attic = malloc(upd->Total_Num_Matrices * sizeof(int *));
  a = malloc(upd->Total_Num_Matrices * sizeof(double *));
  a_old = malloc(upd->Total_Num_Matrices * sizeof(double *));

  if (strcmp(Matrix_Format, "epetra") == 0) {
    err = check_compatible_solver();
    EH(err,
        "Incompatible matrix solver for epetra, epetra supports amesos and aztecoo solvers.");
    check_parallel_error("Matrix format / Solver incompatibility");
    for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
      ams[pg->imtrx]->RowMatrix = EpetraCreateRowMatrix(
          num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx]);
      EpetraCreateGomaProblemGraph(ams[pg->imtrx], exo, dpi);
    }
  } else if (strcmp(Matrix_Format, "msr") == 0) {

    log_msg("alloc_MSR_sparse_arrays...")
;    for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
      alloc_MSR_sparse_arrays(&(ija[pg->imtrx]), &(a[pg->imtrx]),
          &(a_old[pg->imtrx]), 0, node_to_fill, exo, dpi);

      /*
       * An attic to store external dofs column names is needed when
       * running in parallel.
       */

      alloc_extern_ija_buffer(num_universe_dofs[pg->imtrx],
          num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx],
          ija[pg->imtrx], &(ija_attic[pg->imtrx]));

      /*
       * Any necessary one time initialization of the linear
       * solver package (Aztec).
       */
      ams[pg->imtrx]->bindx = ija[pg->imtrx];
      ams[pg->imtrx]->val = a[pg->imtrx];
      ams[pg->imtrx]->belfry = ija_attic[pg->imtrx];
      ams[pg->imtrx]->val_old = a_old[pg->imtrx];

      /*
       * These point to nowhere since we're using MSR instead of VBR
       * format.
       */

      ams[pg->imtrx]->indx = NULL;
      ams[pg->imtrx]->bpntr = NULL;
      ams[pg->imtrx]->rpntr = NULL;
      ams[pg->imtrx]->cpntr = NULL;

      ams[pg->imtrx]->npn = dpi->num_internal_nodes + dpi->num_boundary_nodes;
      ams[pg->imtrx]->npn_plus = dpi->num_internal_nodes
          + dpi->num_boundary_nodes + dpi->num_external_nodes;

      ams[pg->imtrx]->npu = num_internal_dofs[pg->imtrx]
          + num_boundary_dofs[pg->imtrx];
      ams[pg->imtrx]->npu_plus = num_universe_dofs[pg->imtrx];

      ams[pg->imtrx]->nnz = ija[pg->imtrx][num_internal_dofs[pg->imtrx]
          + num_boundary_dofs[pg->imtrx]] - 1;
      ams[pg->imtrx]->nnz_plus = ija[pg->imtrx][num_universe_dofs[pg->imtrx]];

    }
  } else {
    EH(-1, "Attempted to allocate unknown sparse matrix format");
  }

  /* Read initial values from exodus file */
  for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
    init_vec(x[pg->imtrx], cx, exo, dpi, NULL, 0, &timeValueRead[pg->imtrx]);
  }

  /***************************************************************************
   *            STEADY STATE SOLUTION PROCEDURE
   ***************************************************************************/
  if (TimeIntegration == STEADY) {

    theta = 0.0; /* for steady problems. theta def in rf_fem.h */
    delta_t = 0.0;

    /* Right now, we are solving segregated problems in numerical order fashion
     * Matrix 0, then 1, then so on. In the future we could change that.
     */
    for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {

      find_and_set_Dirichlet(x[pg->imtrx], xdot[pg->imtrx], exo, dpi);

      matrix_systems_mask = 1;

      log_msg("sl_init()...")
;      sl_init(matrix_systems_mask, ams, exo, dpi, cx);

#ifdef PARALLEL
      /*
       * Make sure the solver was properly initialized on all processors.
       */
      check_parallel_error("Solver initialization problems");
#endif /* PARALLEL */

      err = solve_linear_segregated(ams[pg->imtrx], x[pg->imtrx], delta_t,
          theta, x_old[pg->imtrx], x_older[pg->imtrx], xdot[pg->imtrx],
          xdot_old[pg->imtrx], resid_vector[pg->imtrx], x_update[pg->imtrx],
          scale[pg->imtrx], &converged, &nprint, gv, time1, exo, dpi, cx, 0,
          &time_step_reform);
      /*      err = solve_nonlinear_problem(ams[pg->imtrx], x[pg->imtrx], delta_t,
       theta, x_old[pg->imtrx], x_older[pg->imtrx], xdot[pg->imtrx],
       xdot_old[pg->imtrx], resid_vector[pg->imtrx], x_update[pg->imtrx],
       scale[pg->imtrx], &converged, &nprint, tev[pg->imtrx],
       tev_post[pg->imtrx], gv, rd[pg->imtrx], NULL, NULL, gvec[pg->imtrx],
       gvec_elem, time1, exo, dpi, cx, 0, &time_step_reform, is_steady_state,
       NULL, NULL, time1, NULL,
       NULL, NULL, NULL); */

    }

    write_solution_segregated(ExoFileOut, resid_vector, x, x_old, xdot,
        xdot_old, tev, tev_post, gv, rd, gvec, gvec_elem, &nprint, delta_t,
        theta, time1, NULL, exo, dpi);
  }
  /********************************************************************************
   *                            Transient solution process 
   ********************************************************************************/
  else {
    if (Debug_Flag && ProcID == 0) {
      fprintf(stderr, "MaxTimeSteps: %d \tTimeMax: %f\n", MaxTimeSteps,
          TimeMax);
      fprintf(stderr, "solving transient problem\n");
    }

    /*
     *  Transfer information from the Transient_Information structure to local variables
     */
    Delta_t0 = tran->Delta_t0;
    Delta_t_min = tran->Delta_t_min;
    Delta_t_max = tran->Delta_t_max;
    MaxTimeSteps = tran->MaxTimeSteps;
    TimeMax = tran->TimeMax;
    eps = tran->eps;
#ifndef COUPLED_FILL
    exp_subcycle = tran->exp_subcycle;
#endif /* not COUPLED_FILL */   

    // Determine if we are using a constant time step or not
    if (Delta_t0 < 0.0) {
      Delta_t0 = -Delta_t0;
      const_delta_t = 1;
    } else {
      const_delta_t = 0;
    }

    time = time1 = tran->init_time; /* Allow non-zero initial time */
    tran->time_value = tran->time_value_old = time1;
    tran->delta_t_old = Delta_t0;
    if (Delta_t0 > Delta_t_max)
      Delta_t0 = Delta_t_max;
    delta_t = delta_t_old = delta_t_older = Delta_t0;
    tran->delta_t = delta_t; /*Load this up for use in load_fv_mesh_derivs */
    tran->delta_t_avg = delta_t;

    /*
     *  Allocate space for prediction vector to be saved here,
     *  since it is only used locally
     */
    for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
        x_pred[pg->imtrx] = alloc_dbl_1(numProcUnknowns[pg->imtrx], 0.0);


      /*
       * set boundary conditions on the initial conditions to the
       * transient calculation.
       *  NOTE -> At this point, xdot[] is set to zero. Therefore,
       *          there may be serious errors in the specification
       *          of the boundary condition at this point. Some
       *          ODE solvers actually solve an initial problem for
       *          the evaluation of xdot[] at t = 0+. This algorithm
       *          perhaps could be introduced here.
       */
      find_and_set_Dirichlet(x[pg->imtrx], xdot[pg->imtrx], exo, dpi);

      /*
       * Before propagating x_n back into the historical records in
       * x_n-1, x_n-2 and x_n-3, make sure external dofs are the best
       * they can be. That is, ask the processors that are supposed to
       * know...
       */
      exchange_dof(cx, dpi, x[pg->imtrx]);

      /*
       * Now copy the initial solution, x[], into the history solutions
       * x_old[], etc. Note, xdot[] = xdot_old[] = 0 at this point,
       * which is in agreement with the specification of the history
       * solutions.
       */
      dcopy1(numProcUnknowns[pg->imtrx], x[pg->imtrx], x_old[pg->imtrx]);
      dcopy1(numProcUnknowns[pg->imtrx], x_old[pg->imtrx], x_older[pg->imtrx]);
      dcopy1(numProcUnknowns[pg->imtrx], x_older[pg->imtrx], x_oldest[pg->imtrx]);
    }

    /* initialize the counters for when to print out data */
    time_print = time;
    step_print = 1;
    matrix_systems_mask = 1;

    af->Sat_hyst_reevaluate = FALSE;

    /*
     * Now, just pass pointer to ams structure with all Aztec stuff
     * bundled inside. At the other end, extract "ija" and "a" as
     * appropriate, but the other items are there now, too.
     * Do this only once if in library mode.
     */
    if (callnum == 1)
      for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
        sl_init(matrix_systems_mask, ams, exo, dpi, cx);
      }

    /*
     * make sure the Aztec was properly initialized
     */
    check_parallel_error("Aztec Initialization");

    /* Set the number of successful time steps, nt, to zero */
    nt = 0;
    time_step_reform = Time_Jacobian_Reformation_stride;
    const_delta_ts = const_delta_t;
    last_renorm_nt = 0;

    /*******************************************************************
     *  TOP OF THE TIME STEP LOOP -> Loop over time steps whether
     *                               they be successful or not
     *******************************************************************/
    for (n = 0; n < MaxTimeSteps; n++) {
      /*
       * Calculate the absolute time for the current step, time1
       */
      time1 = time + delta_t;

      if (time1 > TimeMax) {
        DPRINTF(stderr, "\t\tLAST TIME STEP!\n");
        time1 = TimeMax;
        delta_t = time1 - time;
        tran->delta_t = delta_t;
        tran->delta_t_avg = 0.25
            * (delta_t + delta_t_old + delta_t_older + delta_t_oldest);
      }
      tran->time_value = time1;

      for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {

        /*
         * What is known at this exact point in the code:
         *
         *
         *  At time = time, x_old[] = the solution
         *                  xdot_old[] = derivative of the solution
         *  At time = time - delta_t_old:
         *                  x_older[] = the solution
         *                  xdot_older[] = derivative of the solution
         *  At time = time -  delta_t_old -  delta_t_older
         *                  x_oldest[] = the solution
         *                  xdot_oldest[] = derivative of the solution
         *  The value of x[] and xdot[] contain ambivalent information
         *  at this point.
         *
         *  We seek the solution at time = time1 = time + delta_t
         *  by first obtaining a predicted solution x_pred[] with
         *  associated xdot[], and then solving a corrected solution,
         *  x[], with associated time derivative, xdot[].
         *
         *  Note, we may be here due to a failed time step. In this
         *  case x[] and xdot[] will be filled with garbage. For a
         *  previously completed time step, x[] and xdot[] will be
         *  equal to x_old[] and xdot_old[].
         */

        /*
         * SMD 1/24/11
         * If external field is time_dep update the current solution,
         * x_old, to the values of the external variables at that time point.
         */
        if (efv->ev) {
          timeValueReadTrans = time;
          int w;
          for (w = 0; w < efv->Num_external_field; w++) {
            if (strcmp(efv->field_type[w], "transient") == 0) {
              err = rd_trans_vectors_from_exoII(x_old[pg->imtrx], efv->file_nm[w],
                  w, n, &timeValueReadTrans, cx, dpi);
              if (err != 0) {
                DPRINTF(stderr, "%s: err from rd_trans_vectors_from_exoII\n", yo);
              }
            }
          }
        }

        /*
         * Get started with forward/Backward Euler predictor-corrector
         * to damp out any bad things
         */
        if ((nt - last_renorm_nt) == 0) {
          theta = 0.0;
          const_delta_t = 1.0;

        } else if ((nt - last_renorm_nt) >= 3) {
          /* Now revert to the scheme input by the user */
          theta = tran->theta;
          const_delta_t = const_delta_ts;
          /*
           * If the previous step failed due to a convergence error
           * or time step truncation error, then revert to a
           * Backwards-Euler method to restart the calculation
           * using a smaller time step.
           * -> standard ODE solver trick (HKM -> Haven't
           *    had time to benchmark this. Will leave it commented
           *    out).
           *
           *  if (!converged || !success_dt) {
           *    theta = 0.0;
           *  }
           */
        }

        /* Reset the node->DBC[] arrays to -1 where set
         * so that the boundary conditions are set correctly
         * at each time step.
         */
        nullify_dirichlet_bcs();

        find_and_set_Dirichlet(x[pg->imtrx], xdot[pg->imtrx], exo, dpi);

        if (ProcID == 0) {
          if (theta == 0.0)
            strcpy(tspstring, "(BE)");
          else if (theta == 0.5)
            strcpy(tspstring, "(CN)");
          else if (theta == 1.0)
            strcpy(tspstring, "(FE)");
          else
            sprintf(tspstring, "(TSP %3.1f)", theta);
          fprintf(stderr, "\n=> Try for soln at t=%g with dt=%g [%d for %d] %s\n",
              time1, delta_t, nt, n, tspstring);
          log_msg("Predicting try at t=%g, dt=%g [%d for %d so far] %s",
          time1, delta_t, nt, n, tspstring)
  ;      }

        /*
         * Predict the solution, x[], and its derivative, xdot[],
         * at the new time, time1, using the old solution, xdot_old[],
         * And its derivatives at the old time, time.
         */

        predict_solution(numProcUnknowns[pg->imtrx], delta_t, delta_t_old,
            delta_t_older, theta, x[pg->imtrx], x_old[pg->imtrx],
            x_older[pg->imtrx], x_oldest[pg->imtrx], xdot[pg->imtrx],
            xdot_old[pg->imtrx], xdot_older[pg->imtrx]);

        /*
         * Now, that we have a predicted solution for the current
         * time, x[], exchange the degrees of freedom to update the
         * ghost node information.
         */
        exchange_dof(cx, dpi, x[pg->imtrx]);
        exchange_dof(cx, dpi, xdot[pg->imtrx]);

        /*
         *  Set dirichlet conditions in some places. Note, I believe
         *  this step can change the solution vector
         */
        find_and_set_Dirichlet(x[pg->imtrx], xdot[pg->imtrx], exo, dpi);

        /*
         *  HKM -> I don't know if this extra exchange operation
         *         is needed or not. It was originally in the
         *         algorithm. It may be needed if find_and_set..()
         *         changes the solution vector. However, it would
         *         seem to me that we could get rid of the duplication
         *         of effort here.
         *         -> I also added an exchange of xdot[], because
         *            if x[] is needed to be exchanged, then xdot[] must
         *            be exchanged as well.
         */

        exchange_dof(cx, dpi, x[pg->imtrx]);
        exchange_dof(cx, dpi, xdot[pg->imtrx]);

        /*
         * Save the predicted solution for the time step
         * norm calculation to be carried out after convergence
         * of the nonlinear implicit problem
         */
        dcopy1(numProcUnknowns[pg->imtrx], x[pg->imtrx], x_pred[pg->imtrx]);

        /*
         *  Solve the nonlinear problem. If we achieve convergence,
         *  set the flag, converged, to true on return. If not
         *  set the flag to false.

         */
        err = solve_nonlinear_problem(ams[pg->imtrx], x[pg->imtrx], delta_t,
               theta, x_old[pg->imtrx], x_older[pg->imtrx], xdot[pg->imtrx],
               xdot_old[pg->imtrx], resid_vector[pg->imtrx], x_update[pg->imtrx],
               scale[pg->imtrx], &converged, &nprint, tev[pg->imtrx],
               tev_post[pg->imtrx], gv, rd[pg->imtrx], NULL, NULL, gvec[pg->imtrx],
               gvec_elem, time1, exo, dpi, cx, 0, &time_step_reform, 0,
               NULL, NULL, time1, NULL,
               NULL, NULL, NULL);
        /*
        err = solve_linear_segregated(ams[pg->imtrx], x[pg->imtrx], delta_t,
            theta, x_old[pg->imtrx], x_older[pg->imtrx], xdot[pg->imtrx],
            xdot_old[pg->imtrx], resid_vector[pg->imtrx], x_update[pg->imtrx],
            scale[pg->imtrx], &converged, &nprint, gv, time1, exo, dpi, cx, n,
            &time_step_reform);
*/
        if (err == -1)
          converged = FALSE;
        inewton = err;
        evpl_glob[0]->update_flag = 0; /*See get_evp_stress_tensor for description */
        af->Sat_hyst_reevaluate = FALSE; /*See load_saturation for description*/

        /*
         * HKM -> I do not know if these operations are needed. I added
         *        an exchange of xdot[] here, because if x[] is exchanged
         *        then xdot needs to be exchanged as well.
         */

        exchange_dof(cx, dpi, x[pg->imtrx]);
        exchange_dof(cx, dpi, xdot[pg->imtrx]);

        if (!converged) break;
      }

      if (converged)
        af->Sat_hyst_reevaluate = TRUE; /*see load_saturation */

      /* Check element quality */
      good_mesh = element_quality(exo, x[pg->imtrx], ams[0]->proc_config);

      /*
       * Check the time step truncation error.
       * If too large, set the success_dt flag to FALSE. We will
       * then not accept the current time step, reduced delta_t,
       * and retry.
       */
      if (converged) {
        /* Assume we want the minimum delta_t_new */
        delta_t_new = 1e20;
        for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
          double mat_dt_new;
          mat_dt_new = time_step_control(delta_t, delta_t_old, const_delta_t,
              x[pg->imtrx], x_pred[pg->imtrx], x_old[pg->imtrx], NULL, NULL, eps,
              &success_dt, tran->use_var_norm);

          if (mat_dt_new < delta_t_new) {
            delta_t_new = mat_dt_new;
          }
        }
        if (const_delta_t) {
          success_dt = TRUE;
          delta_t_new = delta_t;
        } else if (failed_recently_countdown > 0) {
          delta_t_new = delta_t;
          failed_recently_countdown--;
        } else if (delta_t_new > Delta_t_max) {
          delta_t_new = Delta_t_max;
          /*        } else if ( !success_dt && delta_t_new < tran->resolved_delta_t_min ) {*/
        } else if (delta_t_new < tran->resolved_delta_t_min) {
          /*          if ( delta_t > tran->resolved_delta_t_min ) {  */
          /* fool algorithm into using delta_t = tran->resolved_delta_t_min */
          delta_t_new = tran->resolved_delta_t_min;
          success_dt = TRUE;
          DPRINTF(stderr, "\n\tminimum resolved step limit!\n");
        }
      }

      if (converged && success_dt) {
        nt += 1;
        time = time1;

        /* Determine whether to print out the data or not */
        i_print = 0;
        if (tran->print_freq == 0) {
          if ((time > time_print)
              || (fabs(time - time_print) < (1.e-4 * tran->print_delt))) {
            if (tran->print_delt2 < 0.) {
              i_print = 1;
              time_print += tran->print_delt;
            } else {
              if (time < tran->print_delt2_time) {
                i_print = 1;
                time_print += tran->print_delt;
              } else {
                i_print = 1;
                time_print += tran->print_delt2;
              }
            }
          }
        } else {
          if (nt == step_print) {
            i_print = 1;
            step_print += tran->print_freq;
          }
        }

        if (time1 >= (ROUND_TO_ONE * TimeMax))
          i_print = 1;

        error = 0;
        if (i_print) {
          if (Write_Intermediate_Solutions == 0) {
            write_solution_segregated(ExoFileOut, resid_vector, x, x_old, xdot,
                    xdot_old, tev, tev_post, gv, rd, gvec, gvec_elem, &nprint, delta_t,
                    theta, time1, NULL, exo, dpi);
            nprint++;
          }
        } /* if(i_print) */
        evpl_glob[0]->update_flag = 1;

        /* Fix output if current time step matches frequency */
        if (step_fix != 0 && nt == step_fix) {
#ifdef PARALLEL
          /* Barrier because fix needs both files to be finished printing
           and fix always occurs on the same timestep as printing */
          MPI_Barrier(MPI_COMM_WORLD);
#endif
          if (ProcID == 0 && Brk_Flag == 1) {
            fix_output();
          }
          /* Fix step is relative to print step */
          step_fix += tran->fix_freq * tran->print_freq;
        }

        /*
         * Adjust the time step if the new time will be larger than the
         * next printing time.
         */
        if (tran->print_freq == 0 && success_dt) {
          if ((time + 1.2 * delta_t_new >= time_print) && (time_print > time)) {
            delta_t_new = time_print - time;
            DPRINTF(stderr,
                "reset delta_t = %g to maintain printing frequency\n",
                delta_t_new);
            if (delta_t_new <= 0)
              EH(-1, "error with time-step printing control");
          } else if (time >= time_print) {
            if (delta_t_new != tran->print_delt) {
              delta_t_new = tran->print_delt;
              DPRINTF(stderr,
                  "reset delta_t = %g to maintain printing frequency\n",
                  delta_t_new);
              if (delta_t_new <= 0) {
                EH(-1, "error with time-step printing control");
              }
            }
          }
        }

        /*
         *   save xdot to xdot_old for next time step
         */
        for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
          dcopy1(numProcUnknowns[pg->imtrx], xdot_old[pg->imtrx],
              xdot_older[pg->imtrx]);
          if (tran->solid_inertia)
            dcopy1(numProcUnknowns[pg->imtrx], tran->xdbl_dot,
                tran->xdbl_dot_old);
          dcopy1(numProcUnknowns[pg->imtrx], xdot[pg->imtrx],
              xdot_old[pg->imtrx]);
          dcopy1(numProcUnknowns[pg->imtrx], x_older[pg->imtrx],
              x_oldest[pg->imtrx]);
          dcopy1(numProcUnknowns[pg->imtrx], x_old[pg->imtrx],
              x_older[pg->imtrx]);
          dcopy1(numProcUnknowns[pg->imtrx], x[pg->imtrx], x_old[pg->imtrx]);
        }
        delta_t_oldest = delta_t_older;
        delta_t_older = delta_t_old;
        delta_t_old = delta_t;
        tran->delta_t_old = delta_t_old;
        tran->time_value_old = time;
        delta_t = delta_t_new;
        tran->delta_t = delta_t; /*load up for use in load_fv_mesh_derivs*/
        tran->delta_t_avg = 0.25
            * (delta_t + delta_t_old + delta_t_older + delta_t_oldest);

        if (time1 >= (ROUND_TO_ONE * TimeMax)) {
          DPRINTF(stderr, "\t\tout of time!\n");
          if (Anneal_Mesh) {
            /*
             * Transform the node point coordinates according to the
             * displacements and write out all the results using the
             * displaced coordinates. Set the displacement field to
             * zero, too.
             */
            for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {

              err = anneal_mesh(x[pg->imtrx], tev[pg->imtrx], tev_post[pg->imtrx],
                  gv, rd[pg->imtrx], time1, exo, dpi);
              EH(err, "anneal_mesh() bad return.");
            }
          }
          goto free_and_clear;
        }
        if (!good_mesh)
          goto free_and_clear;

      } /*  if(converged && success_dt) */

      else /* not converged or unsuccessful time step */
      {
        /* Set bit TRUE in next line to enable retries for failed first timestep*/
        if (relax_bit && nt == 0 && n < 15) {
          if (inewton == -1) {
            DPRINTF(stderr,
                "\nHmm... trouble on first step \n  Let's try some more relaxation  \n");
            if ((damp_factor1 <= 1. && damp_factor1 >= 0.)
                && (damp_factor2 <= 1. && damp_factor2 >= 0.)
                && (damp_factor3 <= 1. && damp_factor3 >= 0.)) {
              custom_tol1 *= 0.01;
              custom_tol2 *= 0.01;
              custom_tol3 *= 0.01;
              DPRINTF(stderr, "  custom tolerances %g %g %g  \n", custom_tol1,
                  custom_tol2, custom_tol3);
            } else {
              damp_factor1 *= 0.5;
              DPRINTF(stderr, "  damping factor %g  \n", damp_factor1);
            }
          } else {
            DPRINTF(stderr,
                "\nHmm... could not converge on first step\n Let's try some more iterations\n");
            dcopy1(numProcUnknowns[pg->imtrx], x[pg->imtrx], x_old[pg->imtrx]);

            if ((damp_factor1 <= 1. && damp_factor1 >= 0.)
                && (damp_factor2 <= 1. && damp_factor2 >= 0.)
                && (damp_factor3 <= 1. && damp_factor3 >= 0.)) {
              custom_tol1 *= 100.;
              custom_tol2 *= 100.;
              custom_tol3 *= 100.;
              DPRINTF(stderr, "  custom tolerances %g %g %g  \n", custom_tol1,
                  custom_tol2, custom_tol3);
            } else {
              damp_factor1 *= 2.0;
              damp_factor1 = MIN(damp_factor1, 1.0);
              DPRINTF(stderr, "  damping factor %g  \n", damp_factor1);
            }
          }
        } else if (delta_t
            < tran->resolved_delta_t_min / tran->time_step_decelerator) {
          DPRINTF(stderr, "\n\tminimum resolved step limit!\n");
          delta_t_oldest = delta_t_older;
          delta_t_older = delta_t_old;
          delta_t_old = delta_t;
          tran->delta_t_old = delta_t_old;
          tran->time_value_old = time;
          delta_t = tran->resolved_delta_t_min;
          tran->delta_t = delta_t;
          tran->delta_t_avg = 0.25
              * (delta_t + delta_t_old + delta_t_older + delta_t_oldest);
          time1 = time + delta_t;
          tran->time_value = time1;
        } else {
          DPRINTF(stderr, "\n\tlast time step failed, dt *= %g for next try!\n",
              tran->time_step_decelerator);

          delta_t *= tran->time_step_decelerator;
          tran->delta_t = delta_t;
          tran->delta_t_avg = 0.25
              * (delta_t + delta_t_old + delta_t_older + delta_t_oldest);
          time1 = time + delta_t;
          tran->time_value = time1;
          evpl_glob[0]->update_flag = 2;
          af->Sat_hyst_reevaluate = 0;

          /* if specified with "Steps of constant delta_t after failure"
           use a constant delta_t to help the painful recovery
           */
          failed_recently_countdown = tran->const_dt_after_failure;
        }

      }

      if (delta_t <= Delta_t_min) {
        DPRINTF(stderr, "\n\tdelta_t = %e < %e\n\n", delta_t, Delta_t_min);

        DPRINTF(stderr, "time step too small, I'm giving up!\n");
        break;
      }

    } /* end of time step loop */
  } /* end of if steady else transient */
  free_and_clear: return;
}
