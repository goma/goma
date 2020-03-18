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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "std.h"
#include "brk_utils.h"
#include "exo_struct.h"
#include "mm_as_structs.h"
#include "rf_fem_const.h"
#include "rf_node_const.h"
#include "sl_epetra_interface.h"
#include "sl_epetra_util.h"
#include "az_aztec.h"
#include "dp_comm.h"
#include "dp_types.h"
#include "dp_utils.h"
#include "el_geom.h"
#include "mm_as.h"
#include "mm_augc_util.h"
#include "mm_bc.h"
#include "mm_eh.h"
#include "mm_fill_ls.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_util.h"
#include "mm_flux.h"
#include "mm_more_utils.h"
#include "mm_mp.h"
#include "mm_mp_structs.h"
#include "mm_post_def.h"
#include "mm_post_proc.h"
#include "mm_sol_nonlinear.h"
#include "mm_unknown_map.h"
#include "mpi.h"
#include "rf_allo.h"
#include "rf_bc.h"
#include "rf_bc_const.h"
#include "rf_fem.h"
#include "rf_io.h"
#include "rf_io_const.h"
#include "rf_io_structs.h"
#include "rf_mp.h"
#include "rf_solve.h"
#include "rf_solver.h"
#include "rf_util.h"
#include "sl_matrix_util.h"
#include "sl_util.h"
#include "sl_util_structs.h"
#include "wr_exo.h"
#include "wr_soln.h"
#include "dpi.h"

#define GOMA_RF_SOLVE_SEGREGATED_C
#include "rf_solve_segregated.h"

#define ROUND_TO_ONE 0.9999999
static int discard_previous_time_step(int num_unks, double *x, double *x_old,
                                      double *x_older, double *x_oldest,
                                      double *xdot, double *xdot_old,
                                      double *xdot_older);

static double vector_distance_vel(int size, double *vec1, double *vec2); 
static double vector_distance_pres(int size, double *vec1, double *vec2);

double vector_distance_squared(int size, double *vec1, double *vec2,
                               int ignore_pressure, int imtrx) {
  double distance_sq = 0;
#ifdef PARALLEL
  double global_distance = 0;
#endif
  double x;
  int i;
  for (i = 0; i < size; i++) {
    if (ignore_pressure && idv[imtrx][i][0] == PRESSURE)
      continue;
    x = (vec1[i] - vec2[i]);
    distance_sq += x * x;
  }
#ifdef PARALLEL
  MPI_Allreduce(&distance_sq, &global_distance, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  distance_sq = global_distance;
#endif /* PARALLEL */
  return distance_sq;
}

double vector_distance(int size, double *vec1, double *vec2) {
  double distance = 0;
#ifdef PARALLEL
  double global_distance = 0;
#endif
  double x;
  int i;
  for (i = 0; i < size; i++) {
    x = (vec1[i] - vec2[i]);
    distance += x * x;
  }
#ifdef PARALLEL
  MPI_Allreduce(&distance, &global_distance, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  distance = global_distance;
#endif /* PARALLEL */
  return sqrt(distance);
}

static double vector_distance_vel(int size, double *vec1, double *vec2) {
  double distance = 0;
#ifdef PARALLEL
  double global_distance = 0;
#endif
  double x;
  int i;
  for (i = 0; i < size; i++) {
    if (idv[0][i][0] == VELOCITY1 || idv[0][i][0] == VELOCITY2) {
      x = (vec1[i] - vec2[i]);
      distance += x * x;
    }
  }
#ifdef PARALLEL
  MPI_Allreduce(&distance, &global_distance, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  distance = global_distance;
#endif /* PARALLEL */
  return sqrt(distance);
}

static double vector_distance_pres(int size, double *vec1, double *vec2) {
  double distance = 0;
#ifdef PARALLEL
  double global_distance = 0;
#endif
  double x;
  int i;
  for (i = 0; i < size; i++) {
    if (idv[0][i][0] == PRESSURE) {
      x = (vec1[i] - vec2[i]);
      distance += x * x;
    }
  }
#ifdef PARALLEL
  MPI_Allreduce(&distance, &global_distance, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  distance = global_distance;
#endif /* PARALLEL */
  return sqrt(distance);
}
/*************************************************************************************
 *  solve_problem_segregated
 *
 *      	 Routine that controls solution of overall FEM problem in
 * segregated fashion. This routine will control time integration, order of
 * matrices to be solved, matrices fills and nonlinear solver and results
 * output.
 *
 *        Started as a clone from solve_problem routine
 *
 * Revision History
 * ================
 * 26 November 2014 - Kristianto Tjiptowidjojo - Creation.
 *
 * ******************************************************************************/

void solve_problem_segregated(
    Exo_DB *exo, /* ptr to the finite element mesh database  */
    Dpi *dpi,    /* distributed processing information       */
    dbl *te_out) /* te_out - return actual end time */
{

  /*
   * Sparse matrix storage vectors
   * (MSR format.  See "SPARSKIT: a basic tool kit for sparse matrix
   * computations" by Youcef Saad)
   */
  int **ija = NULL;       /* column pointer array                     */
  double **a = NULL;      /* nonzero array                            */
  double **a_old = NULL;  /* nonzero array                            */
  int **ija_attic = NULL; /* storage for external dofs                */

  /*
   * Variables
   */

  static double **x = NULL; /* solution vector                   */

  static double **x_old = NULL;      /* old solution vector               */
  static double **x_older = NULL;    /* older solution vector             */
  static double **x_oldest = NULL;   /* oldest solution vector saved      */
  static double **xdot = NULL;       /* current time derivative of soln   */
  static double **xdot_old = NULL;   /* old time derivative of soln       */
  static double **xdot_older = NULL; /* old time derivative of soln       */
  double **x_previous = NULL;

  double **x_pred = NULL;

  double **x_update = NULL;          /* update at last iteration          */
  double **resid_vector = NULL;      /* residual                          */
  double **resid_vector_sens = NULL; /* residual sensitivity              */
  double **scale = NULL;             /* scale vector for modified newton  */
  double **delta_x = NULL;

  int *node_to_fill = NULL;

  static struct Aztec_Linear_Solver_System **ams;
  /*
   * Variables for time integration
   */
  double *timeValueRead;
  double timeValueReadTrans = 0.0;
  double time_print, i_print;
  double theta = 0.0, time = 0;
  static double time1 = 0.0; /* Current time that the simulation is trying  to
                                find the solution for */
  double delta_t_new = 0, delta_t = 0, delta_t_old = 0, delta_t_older = 0,
         delta_t_oldest = 0;
  char tspstring[MAX_FNL];
  int n;
  int nt;
  int last_renorm_nt;   /* time step at which last renorm occured   */
  int time_step_reform; /* counter for jacobian reformation stride  */
  int converged = TRUE; /* success or failure of Newton iteration   */
  static int nprint = 0;
  int step_fix = 0;
  int const_delta_t, const_delta_ts, step_print;
  int success_dt;
  int failed_recently_countdown = 0;
  int num_total_nodes;
  /*
   * Local variables
   */
  int good_mesh;
  int error, err = 0;
  static double **gvec = NULL;
  static double ****gvec_elem = NULL;
  static struct Results_Description **rd;

  double *gv = NULL; /* Global variable values for ExoII database */
  int *tnv;          /* total number of nodal variables and kinds */
  int *tev;          /* total number of elem variables and kinds  */
  int *tnv_post;     /* total number of nodal variables and kinds
       for post processing                       */
  int *tev_post;     /* total number of elem variables and kinds
       for post processing                       */

  unsigned int matrix_systems_mask = 1;

  int i;
  int inewton;
  int *numProcUnknowns;

#ifdef RELAX_ON_TRANSIENT_PLEASE
  int relax_bit =
      TRUE; /* Enables relaxation after a transient convergence failure*/
#else
  int relax_bit = FALSE;
#endif

  int totalnAC;
  int *matrix_nAC;
  int iAC;                      /* Counter                                  */
  double **x_AC = NULL;         /* Solution vector of extra unknowns          */
  double **x_AC_old = NULL;     /* old solution vector of extra unknowns      */
  double **x_AC_older = NULL;   /* older solution vector of extra unknowns    */
  double **x_AC_oldest = NULL;  /* oldest solution vector of extra unknowns   */
  double **x_AC_dot = NULL;     /* current time derivative of extra unknowns  */
  double **x_AC_dot_old = NULL; /* old time derivative of extra unknowns      */
  double **x_AC_dot_older =
      NULL;                  /* Older time derivative of extra unknowns    */
  double **x_AC_pred = NULL; /* predicted extraunknowns */

  struct AC_Information **matrix_augc;

  int did_renorm; /* Flag indicating if we renormalized.       */
  int Renorm_Now =
      FALSE; /* Flag forcing renormalization regardless of gradient */
  int Fill_Matrix = 0;
  int timestep_subcycle = 0;

  double time2 = 0.0;
  struct Extended_Shape_Fcn_Basics **matrix_xfem = NULL;
  tran->solid_inertia = 0;
  static int callnum = 1; /* solve_problem_segregated call counter */

  static const char yo[] =
      "solve_problem_segregated"; /* So my name is in a string.        */

  timeValueRead = calloc(upd->Total_Num_Matrices, sizeof(double));
  tnv = malloc(sizeof(int) * upd->Total_Num_Matrices);
  tev = malloc(sizeof(int) * upd->Total_Num_Matrices);
  tnv_post = malloc(sizeof(int) * upd->Total_Num_Matrices);
  tev_post = malloc(sizeof(int) * upd->Total_Num_Matrices);
  gvec_elem = malloc(sizeof(double ***) * upd->Total_Num_Matrices);

  if (Num_Proc > 1 && tran->fix_freq > 0) {
    step_fix =
        1; /* Always fix on the first timestep to match print frequency */
  }
  pg->matrices = malloc(sizeof(struct Matrix_Data) * upd->Total_Num_Matrices);
  num_total_nodes = dpi->num_universe_nodes;
  /*
   *            BEGIN EXECUTION
   */

  int is_steady_state = (TimeIntegration == STEADY) ? TRUE : FALSE;

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

    tnv[pg->imtrx] =
        cnt_nodal_vars(); /*  tnv_post is calculated in load_nodal_tkn*/
    tev[pg->imtrx] =
        cnt_elem_vars(); /*  tev_post is calculated in load_elem_tkn*/

    if (tnv[pg->imtrx] < 0) {
      DPRINTF(stderr, "%s:\tbad tnv.\n", yo);
      EH(GOMA_ERROR, "\t");
    }

    if (tev[pg->imtrx] < 0) {
      DPRINTF(stderr, "%s:\tMaybe bad tev? See goma design committee ;) \n",
              yo);
      EH(GOMA_ERROR, "\t");
    }

    rd[pg->imtrx]->nev = 0; /* number element variables in results */
    rd[pg->imtrx]->ngv = 0; /* number global variables in results  */
    rd[pg->imtrx]->nhv = 0; /* number history variables in results */

    rd[pg->imtrx]->ngv = 5 + nAC; /* number global variables in results
                                   * see load_global_var_info for names
                                   */

    error = load_global_var_info(rd[pg->imtrx], 0, "CONV");
    error = load_global_var_info(rd[pg->imtrx], 1, "NEWT_IT");
    error = load_global_var_info(rd[pg->imtrx], 2, "MAX_IT");
    error = load_global_var_info(rd[pg->imtrx], 3, "CONVRATE");
    error = load_global_var_info(rd[pg->imtrx], 4, "MESH_VOLUME");

    if (rd[pg->imtrx]->ngv > MAX_NGV)
      EH(GOMA_ERROR, "Augmenting condition values overflowing MAX_NGV.  Change and "
             "rerun .");

    if (nAC > 0) {
      char name[20];

      for (i = 0; i < nAC; i++) {
        sprintf(name, "AUGC_%d", i + 1);
        error = load_global_var_info(rd[pg->imtrx], 5 + i, name);
      }
    }

    if (gv == NULL) {
      gv = alloc_dbl_1(rd[pg->imtrx]->ngv, 0.0);
    }

    /*
     *  Load output nodal types, kinds and names into the structure
     *  which will be used to define what's in the output file.
     */
    error =
        load_nodal_tkn(rd[pg->imtrx], &tnv[pg->imtrx], &tnv_post[pg->imtrx]);
    if (error != 0) {
      DPRINTF(stderr, "%s:  problem with load_nodal_tkn()\n", yo);
      EH(GOMA_ERROR, "\t");
    }

    /*
     *  Load output element var types, kinds and names into the structure
     *  which will be used to define what's in the output file.
     */
    error =
        load_elem_tkn(rd[pg->imtrx], exo, tev[pg->imtrx], &tev_post[pg->imtrx]);
    if (error != 0) {
      DPRINTF(stderr, "%s:  problem with load_elem_tkn()\n", yo);
      EH(GOMA_ERROR, "\t");
    }

    /*
     * Write out the names of the nodal variables that we will be sending to
     * the EXODUS II output file later - do only once if in library mode.
     */

    if (callnum == 1) {
      gvec_elem[pg->imtrx] = (double ***)alloc_ptr_1(exo->num_elem_blocks);
      if ((tev[pg->imtrx] + tev_post[pg->imtrx]) > 0) {
        for (i = 0; i < exo->num_elem_blocks; i++) {
          gvec_elem[pg->imtrx][i] =
              (double **)alloc_ptr_1(tev[pg->imtrx] + tev_post[pg->imtrx]);
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
    numProcUnknowns[pg->imtrx] =
        NumUnknowns[pg->imtrx] + NumExtUnknowns[pg->imtrx];
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
    ams = (struct Aztec_Linear_Solver_System **)alloc_ptr_1(
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

    timestep_subcycle = 0;
    for (int imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
      double sub_delta_t = 1.0;
      int num_sub_time_steps = 1;
      if (pg->matrix_subcycle_count[imtrx] < 1) {
        EH(GOMA_ERROR, "Subcycle count expected to be > 0");
      } else {
        num_sub_time_steps = pg->matrix_subcycle_count[imtrx];
        sub_delta_t = 1.0 / (pg->matrix_subcycle_count[imtrx]);
      }

      if (num_sub_time_steps > 1 &&
          (upd->SegregatedSubcycles > 1 ||
           (ls != NULL && ls->SubcyclesAfterRenorm > 1))) {
        EH(GOMA_ERROR, "Full Subcycling is not supported with time subcycling of "
               "matrices");
      }

      pg->delta_t_fraction[imtrx] = sub_delta_t;
      if (num_sub_time_steps > 1) {
        timestep_subcycle = 1;
      }
    }

    if (timestep_subcycle) {
      pg->sub_step_solutions =
          malloc(upd->Total_Num_Matrices * sizeof(struct Matrix_Data));
      for (int imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
        if (pg->matrix_subcycle_count[imtrx] > 1) {

          pg->sub_step_solutions[imtrx].x =
              alloc_dbl_1(numProcUnknowns[imtrx], 0.0);
          pg->sub_step_solutions[imtrx].x_old =
              alloc_dbl_1(numProcUnknowns[imtrx], 0.0);
          pg->sub_step_solutions[imtrx].x_older =
              alloc_dbl_1(numProcUnknowns[imtrx], 0.0);
          pg->sub_step_solutions[imtrx].x_oldest =
              alloc_dbl_1(numProcUnknowns[imtrx], 0.0);
          pg->sub_step_solutions[imtrx].xdot =
              alloc_dbl_1(numProcUnknowns[imtrx], 0.0);
          pg->sub_step_solutions[imtrx].xdot_old =
              alloc_dbl_1(numProcUnknowns[imtrx], 0.0);
          pg->sub_step_solutions[imtrx].xdot_older =
              alloc_dbl_1(numProcUnknowns[imtrx], 0.0);
          pg->sub_step_solutions[imtrx].x_update =
              alloc_dbl_1(numProcUnknowns[imtrx], 0.0);
        }
        }
      }

      x = malloc(upd->Total_Num_Matrices * sizeof(double *));
      x_old = malloc(upd->Total_Num_Matrices * sizeof(double *));
      x_older = malloc(upd->Total_Num_Matrices * sizeof(double *));
      x_oldest = malloc(upd->Total_Num_Matrices * sizeof(double *));
      xdot = malloc(upd->Total_Num_Matrices * sizeof(double *));
      xdot_old = malloc(upd->Total_Num_Matrices * sizeof(double *));
      xdot_older = malloc(upd->Total_Num_Matrices * sizeof(double *));
      x_pred = malloc(upd->Total_Num_Matrices * sizeof(double *));
      delta_x = malloc(upd->Total_Num_Matrices * sizeof(double *));
      x_previous = malloc(upd->Total_Num_Matrices * sizeof(double *));
      for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
        x[pg->imtrx] = alloc_dbl_1(numProcUnknowns[pg->imtrx], 0.0);
        x_old[pg->imtrx] = alloc_dbl_1(numProcUnknowns[pg->imtrx], 0.0);
        x_older[pg->imtrx] = alloc_dbl_1(numProcUnknowns[pg->imtrx], 0.0);
        x_oldest[pg->imtrx] = alloc_dbl_1(numProcUnknowns[pg->imtrx], 0.0);
        xdot[pg->imtrx] = alloc_dbl_1(numProcUnknowns[pg->imtrx], 0.0);
        xdot_old[pg->imtrx] = alloc_dbl_1(numProcUnknowns[pg->imtrx], 0.0);
        xdot_older[pg->imtrx] = alloc_dbl_1(numProcUnknowns[pg->imtrx], 0.0);
        x_pred[pg->imtrx] = alloc_dbl_1(numProcUnknowns[pg->imtrx], 0.0);
        delta_x[pg->imtrx] = alloc_dbl_1(numProcUnknowns[pg->imtrx], 0.0);
        x_previous[pg->imtrx] = alloc_dbl_1(numProcUnknowns[pg->imtrx], 0.0);
      }
    }

    x_update = malloc(upd->Total_Num_Matrices * sizeof(double *));
    for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
      x_update[pg->imtrx] = alloc_dbl_1(
          numProcUnknowns[pg->imtrx] + numProcUnknowns[pg->imtrx], 0.0);
    }

    for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
      pg->matrices[pg->imtrx].ams = ams[pg->imtrx];
      pg->matrices[pg->imtrx].x = x[pg->imtrx];
      pg->matrices[pg->imtrx].x_old = x_old[pg->imtrx];
      pg->matrices[pg->imtrx].x_older = x_older[pg->imtrx];
      pg->matrices[pg->imtrx].xdot = xdot[pg->imtrx];
      pg->matrices[pg->imtrx].xdot_old = xdot_old[pg->imtrx];
      pg->matrices[pg->imtrx].x_update = x_update[pg->imtrx];
      pg->matrices[pg->imtrx].scale = scale[pg->imtrx];
      pg->matrices[pg->imtrx].resid_vector = resid_vector[pg->imtrx];
    }

    /* Allocate sparse matrix */

    ija = malloc(upd->Total_Num_Matrices * sizeof(int *));
    ija_attic = malloc(upd->Total_Num_Matrices * sizeof(int *));
    a = malloc(upd->Total_Num_Matrices * sizeof(double *));
    a_old = malloc(upd->Total_Num_Matrices * sizeof(double *));

    if (strcmp(Matrix_Format, "epetra") == 0) {
      err = check_compatible_solver();
      EH(err,
         "Incompatible matrix solver for epetra, epetra supports amesos and "
         "aztecoo solvers.");
      check_parallel_error("Matrix format / Solver incompatibility");
      for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
        ams[pg->imtrx]->RowMatrix = EpetraCreateRowMatrix(
            num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx]);
        EpetraCreateGomaProblemGraph(ams[pg->imtrx], exo, dpi);
      }
    } else if (strcmp(Matrix_Format, "msr") == 0) {

      log_msg("alloc_MSR_sparse_arrays...");
      for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
        alloc_MSR_sparse_arrays(&(ija[pg->imtrx]), &(a[pg->imtrx]),
                                &(a_old[pg->imtrx]), 0, node_to_fill, exo, dpi);

        /*
         * An attic to store external dofs column names is needed when
         * running in parallel.
         */

        alloc_extern_ija_buffer(num_universe_dofs[pg->imtrx],
                                num_internal_dofs[pg->imtrx] +
                                    num_boundary_dofs[pg->imtrx],
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
        ams[pg->imtrx]->npn_plus = dpi->num_internal_nodes +
                                   dpi->num_boundary_nodes +
                                   dpi->num_external_nodes;

        ams[pg->imtrx]->npu =
            num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx];
        ams[pg->imtrx]->npu_plus = num_universe_dofs[pg->imtrx];

        ams[pg->imtrx]->nnz = ija[pg->imtrx][num_internal_dofs[pg->imtrx] +
                                             num_boundary_dofs[pg->imtrx]] -
                              1;
        ams[pg->imtrx]->nnz_plus = ija[pg->imtrx][num_universe_dofs[pg->imtrx]];
      }
    } else {
      EH(GOMA_ERROR, "Attempted to allocate unknown sparse matrix format");
    }

    double *global_x_AC = NULL;

    if (nAC > 0) {
      global_x_AC = calloc(sizeof(double), nAC);
    }

    /* Read initial values from exodus file */
    for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
      init_vec(x[pg->imtrx], cx[pg->imtrx], exo, dpi, global_x_AC, nAC,
               &timeValueRead[pg->imtrx]);

      if (upd->ep[pg->imtrx][POLYMER_STRESS11] >= 0 && Conformation_Flag == 1) {
        initial_guess_stress_to_log_conf(x[pg->imtrx], num_total_nodes);
      }
    }

    if (TimeIntegration != STEADY) {
      if (tran->init_time < 0.0) {
        tran->init_time = timeValueRead[0];
        DPRINTF(stdout, "\n Initial Simulation Time Has been set to %g\n",
                timeValueRead[0]);
      }
    }

    for (iAC = 0; iAC < nAC; iAC++) {
      if (augc[iAC].Type != AC_USERBC && augc[iAC].Type != AC_FLUX) {
        EH(GOMA_ERROR, "Can only use BC and flux AC's in segregated solve");
      }
    }

    totalnAC = nAC;
    matrix_augc =
        malloc(sizeof(struct AC_Information *) * upd->Total_Num_Matrices);
    matrix_nAC = calloc(sizeof(int), upd->Total_Num_Matrices);

    for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
      matrix_nAC[pg->imtrx] = 0;
      matrix_augc[pg->imtrx] = malloc(sizeof(struct AC_Information) * totalnAC);
    }

    int **invACidx = NULL;

    if (nAC > 0) {
      invACidx = malloc(sizeof(int *) * upd->Total_Num_Matrices);
      int i;
      for (i = 0; i < upd->Total_Num_Matrices; i++) {
        invACidx[i] = malloc(sizeof(int) * nAC);
      }
    }

    for (iAC = 0; iAC < nAC; iAC++) {
      if (augc[iAC].Type == AC_USERBC) {
        int ibc = augc[iAC].BCID;
        int eqn = BC_Types[ibc].equation;

        if (!(eqn >= V_FIRST && eqn < V_LAST)) {
          EH(GOMA_ERROR, "AC BC not associated with an equation, not supported in "
                 "segregated solve");
        }

        int found = FALSE;
        int imtrx;
        // find matrix that has that equation
        for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
          if (upd->vp[imtrx][eqn] > -1) {
            found = TRUE;
            matrix_augc[imtrx][matrix_nAC[imtrx]] = augc[iAC];
            invACidx[imtrx][matrix_nAC[imtrx]] = iAC;
            matrix_nAC[imtrx]++;
          }
        }

        if (!found) {
          EH(GOMA_ERROR, "Could not associate BC AC with a matrix");
        }
      } else if (augc[iAC].Type == AC_FLUX) {
        int found = FALSE;
        int imtrx;
        for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
          if (upd->vp[imtrx][R_MOMENTUM1] > -1) {
            WH(-1, "Associating FLUX AC with momentum matrix");
            found = TRUE;
            matrix_augc[imtrx][matrix_nAC[imtrx]] = augc[iAC];
            invACidx[imtrx][matrix_nAC[imtrx]] = iAC;
            matrix_nAC[imtrx]++;
          }
        }

        if (!found) {
          EH(GOMA_ERROR, "Could not associate FLUX AC with momentum matrix");
        }
      } else {
        EH(GOMA_ERROR, "AC type not supported");
      }
    }

    for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
      for (iAC = 0; iAC < matrix_nAC[pg->imtrx]; iAC++) {
        matrix_augc[pg->imtrx][iAC].d_evol_dx =
            alloc_dbl_1(numProcUnknowns[pg->imtrx], 0.0);
        matrix_augc[pg->imtrx][iAC].d_lsvel_dx =
            alloc_dbl_1(numProcUnknowns[pg->imtrx], 0.0);
        matrix_augc[pg->imtrx][iAC].d_lsvol_dx =
            alloc_dbl_1(numProcUnknowns[pg->imtrx], 0.0);
      }
    }

    /* Allocate AC unknown arrays on the first call */

    x_AC = malloc(sizeof(double *) * upd->Total_Num_Matrices);
    x_AC_old = malloc(sizeof(double *) * upd->Total_Num_Matrices);
    x_AC_older = malloc(sizeof(double *) * upd->Total_Num_Matrices);
    x_AC_oldest = malloc(sizeof(double *) * upd->Total_Num_Matrices);
    x_AC_dot = malloc(sizeof(double *) * upd->Total_Num_Matrices);
    x_AC_dot_old = malloc(sizeof(double *) * upd->Total_Num_Matrices);
    x_AC_dot_older = malloc(sizeof(double *) * upd->Total_Num_Matrices);
    x_AC_pred = malloc(sizeof(double *) * upd->Total_Num_Matrices);
    if (totalnAC > 0) {
      for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
        if (matrix_nAC[pg->imtrx] > 0) {
          x_AC[pg->imtrx] = alloc_dbl_1(matrix_nAC[pg->imtrx], 0.0);
          x_AC_old[pg->imtrx] = alloc_dbl_1(matrix_nAC[pg->imtrx], 0.0);
          x_AC_older[pg->imtrx] = alloc_dbl_1(matrix_nAC[pg->imtrx], 0.0);
          x_AC_oldest[pg->imtrx] = alloc_dbl_1(matrix_nAC[pg->imtrx], 0.0);
          x_AC_dot[pg->imtrx] = alloc_dbl_1(matrix_nAC[pg->imtrx], 0.0);
          x_AC_dot_old[pg->imtrx] = alloc_dbl_1(matrix_nAC[pg->imtrx], 0.0);
          x_AC_dot_older[pg->imtrx] = alloc_dbl_1(matrix_nAC[pg->imtrx], 0.0);
          x_AC_pred[pg->imtrx] = alloc_dbl_1(matrix_nAC[pg->imtrx], 0.0);

          // use initialization
          for (iAC = 0; iAC < matrix_nAC[pg->imtrx]; iAC++) {
            x_AC[pg->imtrx][iAC] = global_x_AC[invACidx[pg->imtrx][iAC]];
          }
        } else {
          x_AC[pg->imtrx] = NULL;
          x_AC_old[pg->imtrx] = NULL;
          x_AC_older[pg->imtrx] = NULL;
          x_AC_oldest[pg->imtrx] = NULL;
          x_AC_dot[pg->imtrx] = NULL;
          x_AC_dot_old[pg->imtrx] = NULL;
          x_AC_dot_older[pg->imtrx] = NULL;
          x_AC_pred[pg->imtrx] = NULL;
        }
      }
    }

    dcopy1(totalnAC, global_x_AC, &(gv[5]));
    /***************************************************************************
     *            STEADY STATE SOLUTION PROCEDURE
     ***************************************************************************/
    if (TimeIntegration == STEADY) {

      theta = 0.0; /* for steady problems. theta def in rf_fem.h */
      delta_t = 0.0;

      for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {

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
        nullify_dirichlet_bcs();
        find_and_set_Dirichlet(x[pg->imtrx], xdot[pg->imtrx], exo, dpi);

        /*
         * Before propagating x_n back into the historical records in
         * x_n-1, x_n-2 and x_n-3, make sure external dofs are the best
         * they can be. That is, ask the processors that are supposed to
         * know...
         */
        exchange_dof(cx[pg->imtrx], dpi, x[pg->imtrx], pg->imtrx);

        /*
         * Now copy the initial solution, x[], into the history solutions
         * x_old[], etc. Note, xdot[] = xdot_old[] = 0 at this point,
         * which is in agreement with the specification of the history
         * solutions.
         */
        dcopy1(numProcUnknowns[pg->imtrx], x[pg->imtrx], x_old[pg->imtrx]);
        dcopy1(numProcUnknowns[pg->imtrx], x_old[pg->imtrx],
               x_older[pg->imtrx]);
        dcopy1(numProcUnknowns[pg->imtrx], x_older[pg->imtrx],
               x_oldest[pg->imtrx]);
      }
      /* Outer loop for number of steady state steps */
      matrix_systems_mask = 1;
      for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
        sl_init(matrix_systems_mask, ams, exo, dpi, cx[pg->imtrx]);
      }

      /*
       * make sure the Aztec was properly initialized
       */
      check_parallel_error("Solver initialization problems");

      for (n = 0; n < tran->MaxSteadyStateSteps; n++) {

        for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
          dcopy1(numProcUnknowns[pg->imtrx], x_older[pg->imtrx],
                 x_oldest[pg->imtrx]);
          dcopy1(numProcUnknowns[pg->imtrx], x_old[pg->imtrx],
                 x_older[pg->imtrx]);
          dcopy1(numProcUnknowns[pg->imtrx], x[pg->imtrx], x_old[pg->imtrx]);
        }

        for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {

          nullify_dirichlet_bcs();

          find_and_set_Dirichlet(x[pg->imtrx], xdot[pg->imtrx], exo, dpi);

          if (ProcID == 0) {
            printf("\n===================== SOLVING MATRIX %d Step %d "
                   "===========================\n\n",
                   pg->imtrx, n);
          }

          nAC = matrix_nAC[pg->imtrx];
          augc = matrix_augc[pg->imtrx];

          err = solve_nonlinear_problem(
              ams[pg->imtrx], x[pg->imtrx], delta_t, theta, x_old[pg->imtrx],
              x_older[pg->imtrx], xdot[pg->imtrx], xdot_old[pg->imtrx],
              resid_vector[pg->imtrx], x_update[pg->imtrx], scale[pg->imtrx],
              &converged, &nprint, tev[pg->imtrx], tev_post[pg->imtrx], gv,
              rd[pg->imtrx], NULL, NULL, gvec[pg->imtrx], gvec_elem[pg->imtrx],
              time1, exo, dpi, cx[pg->imtrx], 0, &time_step_reform,
              is_steady_state, x_AC[pg->imtrx], x_AC_dot[pg->imtrx], time1,
              NULL, NULL, NULL, NULL);

          if (err == -1)
            converged = FALSE;

          if (!converged)
            break;

          if (nAC > 0) {
            DPRINTF(stdout, "\n------------------------------\n");
            DPRINTF(stdout, "Augmenting Conditions:    %4d\n", nAC);
            DPRINTF(stdout, "Number of extra unknowns: %4d\n\n", nAC);

            for (iAC = 0; iAC < nAC; iAC++) {
              if (augc[iAC].Type == AC_USERBC) {
                DPRINTF(stdout, "\tBC[%4d] DF[%4d]=% 10.6e\n", augc[iAC].BCID,
                        augc[iAC].DFID, x_AC[pg->imtrx][iAC]);
              }
              /*		  else if(augc[iAC].Type == AC_USERMAT ||
                                augc[iAC].Type == AC_FLUX_MAT )
                                {
                                DPRINTF(stderr, "\tMT[%4d] MP[%4d]=% 10.6e\n",
                                augc[iAC].MTID, augc[iAC].MPID, x_AC[iAC]);
                                }
                                else if(augc[iAC].Type == AC_VOLUME)
                                {
                                evol_local = augc[iAC].evol;
                                #ifdef PARALLEL
                                if( Num_Proc > 1 ) {
                                MPI_Allreduce( &evol_local, &evol_global, 1,
                                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                                evol_local = evol_global;
                                }
                                #endif
                                DPRINTF(stderr, "\tMT[%4d] VC[%4d]=%10.6e
                 Param=%10.6e\n", augc[iAC].MTID, augc[iAC].VOLID, evol_local,
                                x_AC[iAC]);
                                }
                                else if(augc[iAC].Type == AC_LS_VEL)
                                {
                                evol_local = augc[iAC].lsvol;
                                lsvel_local = augc[iAC].lsvel;


                                #ifdef PARALLEL
                                if( Num_Proc > 1 ) {
                                MPI_Allreduce( &evol_local, &evol_global, 1,
                                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                                MPI_Allreduce( &lsvel_local, &lsvel_global, 1,
                                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                                evol_local = evol_global;
                                lsvel_local = lsvel_global;
                                }
                                #endif

                                lsvel_local = lsvel_local / evol_local;

                                DPRINTF(stderr, "\tMT[%4d] LSVEL
                 phase[%4d]=%10.6e Param=%10.6e\n", augc[iAC].MTID,
                 augc[iAC].LSPHASE, lsvel_local, x_AC[iAC]);
                                }
                                else if(augc[iAC].Type == AC_FLUX)
                                {p
                                DPRINTF(stderr, "\tBC[%4d] DF[%4d]=%10.6e\n",
                                augc[iAC].BCID, augc[iAC].DFID, x_AC[iAC]);
                                } */
            }
          }
        }

        if (converged) {

          for (i = 0; matrix_nAC[pg->imtrx] > 0 && i < matrix_nAC[pg->imtrx];
               i++) {
            gv[5 + invACidx[pg->imtrx][i]] = x_AC[pg->imtrx][i];
          }

          double distance = 0;

          for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices;
               pg->imtrx++) {
            find_and_set_Dirichlet(x[pg->imtrx], xdot[pg->imtrx], exo, dpi);
            double d1 = vector_distance(NumUnknowns[pg->imtrx], x[pg->imtrx],
                                        x_old[pg->imtrx]);
            if (pg->imtrx == 0) {
              double vel = vector_distance_vel(NumUnknowns[pg->imtrx],
                                               x[pg->imtrx], x_old[pg->imtrx]);
              double pres = vector_distance_pres(
                  NumUnknowns[pg->imtrx], x[pg->imtrx], x_old[pg->imtrx]);
              if (ProcID == 0) {
                printf("M%d distance = %g, %g, %g\n", i, d1, vel, pres);
              }
            } else {
              if (ProcID == 0) {
                printf("M%d distance = %g \n", i, d1);
              }
            }
            distance +=
                vector_distance_squared(NumUnknowns[pg->imtrx], x[pg->imtrx],
                                        x_old[pg->imtrx], FALSE, pg->imtrx);
          }

          distance = sqrt(distance);

          if (ProcID == 0) {
            printf("\nL_2 x changes %g\n", distance);
          }

          //	write_solution_segregated(ExoFileOut, resid_vector, x, x_old,
          // xdot, 				  xdot_old, tev, tev_post, gv,
          // rd, gvec, gvec_elem, &nprint, delta_t,
          // theta, (double) nprint, NULL, exo, dpi); 	nprint++;

          if (distance < tran->steady_state_tolerance) {

            pg->imtrx = 0;
            for (i = 0; i < nn_post_fluxes; i++) {
              (void)evaluate_flux(
                  exo, dpi, pp_fluxes[i]->ss_id, pp_fluxes[i]->flux_type,
                  pp_fluxes[i]->flux_type_name, pp_fluxes[i]->blk_id,
                  pp_fluxes[i]->species_number, pp_fluxes[i]->flux_filenm,
                  pp_fluxes[i]->profile_flag, x[pg->imtrx], xdot[pg->imtrx],
                  NULL, delta_t_old, time, 1);
            }

            write_solution_segregated(ExoFileOut, resid_vector, x, x_old, xdot,
                                      xdot_old, tev_post, gv, rd, gvec, &nprint,
                                      delta_t, theta, time1, NULL, exo, dpi);

            if (ProcID == 0) {
              printf("\n Steady state reached \n");
            }
            goto free_and_clear;
          }
        } else {
          break;
        }
      }
      if (ProcID == 0) {
        printf("Failed to solve steady state.");
      }
      goto free_and_clear;

    }
    /********************************************************************************
     *                            Transient solution process
     ********************************************************************************/
    else {

      int renorm_subcycle_count = 0;
      if (Debug_Flag && ProcID == 0) {
        fprintf(stdout, "MaxTimeSteps: %d \tTimeMax: %f\n", tran->MaxTimeSteps,
                tran->TimeMax);
        fprintf(stdout, "solving transient problem\n");
      }

      nt = 0;
      /*
       *  Transfer information from the Transient_Information structure to local
       * variables
       */
      double delta_t0 = tran->Delta_t0;
      double delta_t_min = tran->Delta_t_min;
      double delta_t_max = tran->Delta_t_max;
      int MaxTimeSteps = tran->MaxTimeSteps;
      double TimeMax = tran->TimeMax;
      eps = tran->eps;
#ifndef COUPLED_FILL
      exp_subcycle = tran->exp_subcycle;
#endif /* not COUPLED_FILL */

      // Determine if we are using a constant time step or not
      if (delta_t0 < 0.0) {
        delta_t0 = -delta_t0;
        const_delta_t = 1;
      } else {
        const_delta_t = 0;
      }

      time = time1 = tran->init_time; /* Allow non-zero initial time */
      tran->time_value = tran->time_value_old = time1;
      tran->delta_t_old = delta_t0;
      if (delta_t0 > delta_t_max)
        delta_t0 = delta_t_max;
      delta_t = delta_t_old = delta_t_older = delta_t0;
      tran->delta_t = delta_t; /*Load this up for use in load_fv_mesh_derivs */
      tran->delta_t_avg = delta_t;

      /*
       *  Allocate space for prediction vector to be saved here,
       *  since it is only used locally
       */
      for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
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

        if (nt == 0 && pg->imtrx == 0) {
          xfem = NULL;
          if (upd->XFEM) {
            matrix_xfem = calloc(upd->Total_Num_Matrices,
                                 sizeof(struct Extended_Shape_Fcn_Basics *));
            int imtrx;
            for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
              matrix_xfem[imtrx] =
                  alloc_struct_1(struct Extended_Shape_Fcn_Basics, 1);
              matrix_xfem[imtrx]->ielem = -1;
              matrix_xfem[imtrx]->tot_vol =
                  alloc_dbl_1(numProcUnknowns[imtrx], 0.0);
              matrix_xfem[imtrx]->active_vol =
                  alloc_dbl_1(numProcUnknowns[imtrx], 0.0);
            }
            if (ls == NULL) {
              EH(GOMA_ERROR, "Currently, XFEM requires traditional level set (not pf)");
            }
          }
        }

        if (upd->XFEM) {
          xfem = matrix_xfem[pg->imtrx];
        }

        if (upd->ep[pg->imtrx][FILL] > -1 &&
            nt == 0) { /*  Start of LS initialization */

          Fill_Matrix = pg->imtrx;
          ls->MatrixNum = pg->imtrx;
#ifndef COUPLED_FILL
          EH(GOMA_ERROR, "Segregated not setup for COUPLED_FILL undefined");
#endif /* not COUPLED_FILL */

          if (ls != NULL || pfd != NULL) {

            int eqntype = ls->Init_Method;
            /* This is a temporary loc for this allocation */

            switch (ls->Evolution) {
            case LS_EVOLVE_ADVECT_EXPLICIT:
              DPRINTF(stdout,
                      "\n\t Using decoupled / subcycling for FILL equation.\n");
              break;
            case LS_EVOLVE_ADVECT_COUPLED:
              DPRINTF(stdout, "\n\t Using Coupled Level Set evolution!\n");
              break;
            case LS_EVOLVE_SLAVE:
              DPRINTF(stdout, "\n\t USING SLAVE LEVEL SET INTERFACE\n");
              break;
            case LS_EVOLVE_SEMILAGRANGIAN:
              DPRINTF(stdout,
                      "\n\t Using semi-Lagrangian Level Set Evolution\n");
              break;
            default:
              EH(GOMA_ERROR, "Level Set Evolution scheme not found \n");
            }

            if (ls->Length_Scale < 0.0)
              EH(GOMA_ERROR,
                 "\tError: a Level Set Length Scale needs to be specified\n");

            if (ls->Integration_Depth > 0 || ls->SubElemIntegration ||
                ls->AdaptIntegration) {

              if (ls->Integration_Depth > 0) {
                int first_elem;

                first_elem = find_first_elem_with_var(exo, LS);

                if (first_elem != -1) {
                  load_ei(first_elem, exo, 0, pg->imtrx);

                  Subgrid_Tree = create_shape_fcn_tree(ls->Integration_Depth);
                  DPRINTF(stdout, "\n\tSubgrid Integration of level set "
                                  "interface active.\n");
                }
              } else if (ls->SubElemIntegration) {
                DPRINTF(stdout, "\n\tSubelement Integration of level set "
                                "interface active.\n");
              } else if (ls->AdaptIntegration) {
                DPRINTF(stdout, "\n\tAdaptive Integration of level set "
                                "interface active.\n");
                DPRINTF(stdout, "\tAdaptive Integration Interface Order = %d\n",
                        ls->Adaptive_Order);
              }
              Subgrid_Int.ip_total = 0;
              Subgrid_Int.s = NULL;
              Subgrid_Int.wt = NULL;
            }

            switch (eqntype) {
            case PROJECT:

              DPRINTF(stdout, "\n\t Projection level set initialization \n");

              EH(GOMA_ERROR, "Use of \"PROJECT\" is obsolete.");

              break;

            case EXO_READ:

              DPRINTF(stdout, "\t\t Level set read from exodus database \n");

              break;

            case SURFACES:

              DPRINTF(stdout,
                      "\n\t\t Surface object level set initialization : ");

              /* parallel synchronization of initialization surfaces */
              if (Num_Proc > 1) {
                if (!ls->init_surf_list)
                  ls->init_surf_list = create_surf_list();
                assemble_Global_surf_list(ls->init_surf_list);
              }

              surf_based_initialization(x[pg->imtrx], NULL, NULL, exo,
                                        num_total_nodes, ls->init_surf_list, 0.,
                                        0., 0.);

              DPRINTF(stdout, "- done \n");

              break;

            default:
              WH(-1, "Level Set Initialization method not found \n");
            } /* end of switch( eqntype )  */

            exchange_dof(cx[pg->imtrx], dpi, x[pg->imtrx], pg->imtrx);

            if (converged) {
              switch (ls->Renorm_Method) {

              case HUYGENS:
              case HUYGENS_C:
              case HUYGENS_MASS_ITER:
                Renorm_Now =
                    (ls->Force_Initial_Renorm ||
                     (ls->Renorm_Freq != 0 && ls->Renorm_Countdown == 0));

                did_renorm = huygens_renormalization(
                    x[pg->imtrx], num_total_nodes, exo, cx[pg->imtrx], dpi,
                    num_fill_unknowns, numProcUnknowns[pg->imtrx], time1,
                    Renorm_Now);

                break;

              case CORRECT:

                EH(GOMA_ERROR, "Use of \"CORRECT\" is obsolete.");
                break;
              default:
                if (ls->Evolution == LS_EVOLVE_ADVECT_EXPLICIT ||
                    ls->Evolution == LS_EVOLVE_ADVECT_COUPLED)
                  WH(-1, "No level set renormalization is on.\n");
              } /* end of switch(ls->Renorm_Method ) */
            }
            /*
             * More initialization needed. Have to set those field variables
             * that initially are indexed by level set function.  For example,
             * species concentration and temperature.
             */
            if (ls->Num_Var_Init > 0)
              ls_var_initialization(x, exo, dpi, cx);

            /* 	  DPRINTF(stderr, "Done with ls_var_initialization.\n"); */

            /*
             * Now check to see if we need to build a surface represent.
             * on each time step.  Initialize the structures if so.
             */
            {
              int build = FALSE, ibc = 0;

              while (!build && ibc < Num_BC) {
                build = (BC_Types[ibc].BC_Name == LS_INLET_BC);
                build = build || (BC_Types[ibc].BC_Name == LS_ADC_BC);
                ibc++;
              }

              /* Here we create space for an isosurface list that is updated
               * every time step.
               */

              if (build && ls->last_surf_list == NULL) {
                struct LS_Surf *tmp_surf = NULL;
                struct LS_Surf_Iso_Data *tmp_data;

                ls->last_surf_list = create_surf_list();

                tmp_surf = create_surf(LS_SURF_ISOSURFACE);
                tmp_data = (struct LS_Surf_Iso_Data *)tmp_surf->data;
                tmp_data->isovar = FILL;
                tmp_data->isoval = 0.0;

                append_surf(ls->last_surf_list, tmp_surf);
              }

            } /* matches int build */

          } /* end of ls != NULL */
          dcopy1(numProcUnknowns[pg->imtrx], x[pg->imtrx], x_old[pg->imtrx]);
          dcopy1(numProcUnknowns[pg->imtrx], x[pg->imtrx], x_older[pg->imtrx]);
          dcopy1(numProcUnknowns[pg->imtrx], x[pg->imtrx], x_oldest[pg->imtrx]);

          exchange_dof(cx[pg->imtrx], dpi, x[pg->imtrx], pg->imtrx);
          exchange_dof(cx[pg->imtrx], dpi, x_old[pg->imtrx], pg->imtrx);
          exchange_dof(cx[pg->imtrx], dpi, x_oldest[pg->imtrx], pg->imtrx);

          if (ls != NULL && ls->last_surf_list != NULL) {
            /* Find the interface surf at the last full time step (x_old =
             * tmp_x) for use during this explicit step */
            x_static = x[pg->imtrx];
            x_old_static = x_old[pg->imtrx];
            xdot_static = xdot[pg->imtrx];
            xdot_old_static = xdot_old[pg->imtrx];

            create_subsurfs(ls->last_surf_list, x_old[pg->imtrx], exo);
          }
        }
        /*
         * Before propagating x_n back into the historical records in
         * x_n-1, x_n-2 and x_n-3, make sure external dofs are the best
         * they can be. That is, ask the processors that are supposed to
         * know...
         */
        exchange_dof(cx[pg->imtrx], dpi, x[pg->imtrx], pg->imtrx);

        /*
         * Now copy the initial solution, x[], into the history solutions
         * x_old[], etc. Note, xdot[] = xdot_old[] = 0 at this point,
         * which is in agreement with the specification of the history
         * solutions.
         */
        dcopy1(numProcUnknowns[pg->imtrx], x[pg->imtrx], x_old[pg->imtrx]);
        dcopy1(numProcUnknowns[pg->imtrx], x_old[pg->imtrx],
               x_older[pg->imtrx]);
        dcopy1(numProcUnknowns[pg->imtrx], x_older[pg->imtrx],
               x_oldest[pg->imtrx]);
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
          sl_init(matrix_systems_mask, ams, exo, dpi, cx[pg->imtrx]);
        }

      nt = 0;

      /*
       * make sure the Aztec was properly initialized
       */
      check_parallel_error("Aztec Initialization");

      /* Set the number of successful time steps, nt, to zero */
      nt = 0;
      time_step_reform = Time_Jacobian_Reformation_stride;
      const_delta_ts = const_delta_t;
      last_renorm_nt = 0;

      if (Write_Initial_Solution) {
        write_solution_segregated(ExoFileOut, resid_vector, x, x_old, xdot,
                                  xdot_old, tev_post, gv, rd, gvec, &nprint,
                                  delta_t, theta, time1, NULL, exo, dpi);
        nprint++;
      }

      // Save x into older saves
      for (int imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
        exchange_dof(cx[imtrx], dpi, x[imtrx], imtrx);
        dcopy1(numProcUnknowns[imtrx], x[imtrx], x_old[imtrx]);
        dcopy1(numProcUnknowns[imtrx], x[imtrx], x_older[imtrx]);
        dcopy1(numProcUnknowns[imtrx], x[imtrx], x_oldest[imtrx]);
      }

      /*******************************************************************
       *  TOP OF THE TIME STEP LOOP -> Loop over time steps whether
       *                               they be successful or not
       *******************************************************************/
      for (n = 0; n < MaxTimeSteps; n++) {

        tran->step = n;

        for (int subcycle = 0; subcycle < upd->SegregatedSubcycles ||
                               subcycle < renorm_subcycle_count;
             subcycle++) {

          for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices;
               pg->imtrx++) {
            /*
             * Calculate the absolute time for the current step, time1
             */
            time1 = time + delta_t;

            if (time1 > TimeMax) {
              DPRINTF(stdout, "\t\tLAST TIME STEP!\n");
              time1 = TimeMax;
              delta_t = time1 - time;
              tran->delta_t = delta_t;
              tran->delta_t_avg = 0.25 * (delta_t + delta_t_old +
                                          delta_t_older + delta_t_oldest);
            }
            tran->time_value = time1;

            if (upd->XFEM) {
              xfem = matrix_xfem[pg->imtrx];
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

            if (pg->matrix_subcycle_count[pg->imtrx] > 1) {
              double sub_time = time;

              if (n == 0) {
                pg->sub_delta_t[pg->imtrx] =
                    pg->delta_t_fraction[pg->imtrx] * delta_t;
                pg->sub_delta_t_old[pg->imtrx] =
                    pg->delta_t_fraction[pg->imtrx] * delta_t;
                pg->sub_delta_t_older[pg->imtrx] =
                    pg->delta_t_fraction[pg->imtrx] * delta_t;
                dcopy1(numProcUnknowns[pg->imtrx], x[pg->imtrx],
                       pg->sub_step_solutions[pg->imtrx].x);
                dcopy1(numProcUnknowns[pg->imtrx], x_old[pg->imtrx],
                       pg->sub_step_solutions[pg->imtrx].x_old);
                dcopy1(numProcUnknowns[pg->imtrx], x_older[pg->imtrx],
                       pg->sub_step_solutions[pg->imtrx].x_older);
                dcopy1(numProcUnknowns[pg->imtrx], xdot[pg->imtrx],
                       pg->sub_step_solutions[pg->imtrx].xdot);
                dcopy1(numProcUnknowns[pg->imtrx], xdot_old[pg->imtrx],
                       pg->sub_step_solutions[pg->imtrx].xdot_old);
                exchange_dof(cx[pg->imtrx], dpi,
                             pg->sub_step_solutions[pg->imtrx].x, pg->imtrx);

                exchange_dof(cx[pg->imtrx], dpi,
                             pg->sub_step_solutions[pg->imtrx].x_old,
                             pg->imtrx);

                exchange_dof(cx[pg->imtrx], dpi,
                             pg->sub_step_solutions[pg->imtrx].x_older,
                             pg->imtrx);

                exchange_dof(cx[pg->imtrx], dpi,
                             pg->sub_step_solutions[pg->imtrx].xdot, pg->imtrx);

                exchange_dof(cx[pg->imtrx], dpi,
                             pg->sub_step_solutions[pg->imtrx].xdot_old,
                             pg->imtrx);
              }
              for (int sub_time_step = 0;
                   sub_time_step < pg->matrix_subcycle_count[pg->imtrx];
                   sub_time_step++) {

                pg->sub_delta_t[pg->imtrx] =
                    pg->delta_t_fraction[pg->imtrx] * delta_t;

                double time1 = sub_time + pg->sub_delta_t[pg->imtrx];
                tran->time_value = time1;

                /*
                 * SMD 1/24/11
                 * If external field is time_dep update the current solution,
                 * x_old, to the values of the external variables at that time
                 * point.
                 */
                if (efv->ev) {
                  timeValueReadTrans = sub_time;
                  int w;
                  for (w = 0; w < efv->Num_external_field; w++) {
                    if (strcmp(efv->field_type[w], "transient") == 0) {
                      err = rd_trans_vectors_from_exoII(
                          pg->sub_step_solutions[pg->imtrx].x_old,
                          efv->file_nm[w], w, n, &timeValueReadTrans,
                          cx[pg->imtrx], dpi);
                      if (err != 0) {
                        DPRINTF(stderr,
                                "%s: err from rd_trans_vectors_from_exoII\n",
                                yo);
                      }
                    }
                  }
                }

                /* Reset the node->DBC[] arrays to -1 where set
                 * so that the boundary conditions are set correctly
                 * at each time step.
                 */
                nullify_dirichlet_bcs();

                find_and_set_Dirichlet(pg->sub_step_solutions[pg->imtrx].x,
                                       pg->sub_step_solutions[pg->imtrx].xdot,
                                       exo, dpi);

                exchange_dof(cx[pg->imtrx], dpi,
                             pg->sub_step_solutions[pg->imtrx].x, pg->imtrx);


                exchange_dof(cx[pg->imtrx], dpi,
                             pg->sub_step_solutions[pg->imtrx].xdot, pg->imtrx);


                if (ProcID == 0) {
                  if (theta == 0.0)
                    strcpy(tspstring, "(BE)");
                  else if (theta == 0.5)
                    strcpy(tspstring, "(CN)");
                  else if (theta == 1.0)
                    strcpy(tspstring, "(FE)");
                  else
                    sprintf(tspstring, "(TSP %3.1f)", theta);
                  fprintf(
                      stdout,
                      "\n=> Try for soln at t=%g with dt=%g [%d for %d] %s\n",
                      time1, pg->sub_delta_t[pg->imtrx], nt, n, tspstring);
                  log_msg("Predicting try at t=%g, dt=%g [%d for %d so far] %s",
                          time1, pg->sub_delta_t[pg->imtrx], nt, n, tspstring);
                }

                /*
                 * Predict the solution, x[], and its derivative, xdot[],
                 * at the new time, time1, using the old solution, xdot_old[],
                 * And its derivatives at the old time, time.
                 */
                if (upd->SegregatedSolve && pg->imtrx == 0) {
                  EH(GOMA_ERROR, "Segregated pressure not supported with sub time stepping");
                } else {
                  predict_solution(
                      numProcUnknowns[pg->imtrx], pg->sub_delta_t[pg->imtrx],
                      pg->sub_delta_t_old[pg->imtrx],
                      pg->sub_delta_t_older[pg->imtrx], theta,
                      pg->sub_step_solutions[pg->imtrx].x,
                      pg->sub_step_solutions[pg->imtrx].x_old,
                      pg->sub_step_solutions[pg->imtrx].x_older,
                      pg->sub_step_solutions[pg->imtrx].x_oldest,
                      pg->sub_step_solutions[pg->imtrx].xdot,
                      pg->sub_step_solutions[pg->imtrx].xdot_old,
                      pg->sub_step_solutions[pg->imtrx].xdot_older);
                }

                if (ls != NULL && ls->Evolution == LS_EVOLVE_SLAVE &&
                    pg->imtrx == Fill_Matrix) {
                  surf_based_initialization(pg->sub_step_solutions[pg->imtrx].x,
                                            NULL, NULL, exo, num_total_nodes,
                                            ls->init_surf_list, time1, theta,
                                            delta_t);
                }

                /*
                 * Now, that we have a predicted solution for the current
                 * time, x[], exchange the degrees of freedom to update the
                 * ghost node information.
                 */
                exchange_dof(cx[pg->imtrx], dpi,
                             pg->sub_step_solutions[pg->imtrx].x, pg->imtrx);
                exchange_dof(cx[pg->imtrx], dpi,
                             pg->sub_step_solutions[pg->imtrx].xdot, pg->imtrx);

                if (matrix_nAC[pg->imtrx] > 0) {
                  EH(GOMA_ERROR,
                     "Augmenting conditions not supported for sub time cycles");
                }

                /*
                 *  Set dirichlet conditions in some places. Note, I believe
                 *  this step can change the solution vector
                 */
                find_and_set_Dirichlet(pg->sub_step_solutions[pg->imtrx].x,
                                       pg->sub_step_solutions[pg->imtrx].xdot,
                                       exo, dpi);

                exchange_dof(cx[pg->imtrx], dpi,
                             pg->sub_step_solutions[pg->imtrx].x, pg->imtrx);
                exchange_dof(cx[pg->imtrx], dpi,
                             pg->sub_step_solutions[pg->imtrx].xdot, pg->imtrx);
                /*
                 * Save the predicted solution for the time step
                 * norm calculation to be carried out after convergence
                 * of the nonlinear implicit problem
                 */
                dcopy1(numProcUnknowns[pg->imtrx],
                       pg->sub_step_solutions[pg->imtrx].x, x_pred[pg->imtrx]);

                DPRINTF(
                    stdout,
                    "\n--------------------- SOLVING MATRIX %d sub time step "
                    "%d/%d "
                    "--------------------\n\n",
                    pg->imtrx + 1, sub_time_step + 1,
                    pg->matrix_subcycle_count[pg->imtrx]);

                nAC = matrix_nAC[pg->imtrx];
                augc = matrix_augc[pg->imtrx];

                err = solve_nonlinear_problem(
                    ams[pg->imtrx], pg->sub_step_solutions[pg->imtrx].x,
                    pg->sub_delta_t[pg->imtrx], theta,
                    pg->sub_step_solutions[pg->imtrx].x_old,
                    pg->sub_step_solutions[pg->imtrx].x_older,
                    pg->sub_step_solutions[pg->imtrx].xdot,
                    pg->sub_step_solutions[pg->imtrx].xdot_old,
                    resid_vector[pg->imtrx],
                    pg->sub_step_solutions[pg->imtrx].x_update,
                    scale[pg->imtrx], &converged, &nprint, tev[pg->imtrx],
                    tev_post[pg->imtrx], gv, rd[pg->imtrx], NULL, NULL,
                    gvec[pg->imtrx], gvec_elem[pg->imtrx], time1, exo, dpi,
                    cx[pg->imtrx], 0, &time_step_reform, 0, x_AC[pg->imtrx],
                    x_AC_dot[pg->imtrx], time1, NULL, NULL, NULL, NULL);

                if (err == -1) {
                  converged = FALSE;
                }
                if (!converged) {
                  /* Copy previous solution values if failed timestep */
                  dcopy1(numProcUnknowns[pg->imtrx], pg->sub_step_solutions[pg->imtrx].x_old, pg->sub_step_solutions[pg->imtrx].x);
                  break;
                }

                if (pd_glob[0]->v[pg->imtrx][MOMENT0] ||
                    pd_glob[0]->v[pg->imtrx][MOMENT1] ||
                    pd_glob[0]->v[pg->imtrx][MOMENT2] ||
                    pd_glob[0]->v[pg->imtrx][MOMENT3]) {
                  /*     Floor values to 0 */
                  int floored_values = 0;
                  for (int var = MOMENT0; var <= MOMENT3; var++) {
                    for (i = 0; i < num_total_nodes; i++) {
                      if (pd_glob[0]->v[pg->imtrx][var]) {
                        int j = Index_Solution(i, var, 0, 0, -1, pg->imtrx);

                        if (j != -1 && x[pg->imtrx][j] < 0) {
                          pg->sub_step_solutions[pg->imtrx].x[j] = 0.0;
                          floored_values++;
                        }
                      }
                    }
                  }

                  int global_floored = 0;
                  MPI_Allreduce(&floored_values, &global_floored, 1, MPI_INT,
                                MPI_SUM, MPI_COMM_WORLD);

                  P0PRINTF("Floored %d moment values\n", global_floored);
                }



                sub_time += pg->sub_delta_t[pg->imtrx];
                // change delta t's as this time step may have changed the
                // current delta t
                pg->sub_delta_t_older[pg->imtrx] =
                    pg->sub_delta_t_old[pg->imtrx];
                pg->sub_delta_t_old[pg->imtrx] = pg->sub_delta_t[pg->imtrx];
                // update sub solutions
                dcopy1(numProcUnknowns[pg->imtrx],
                       pg->sub_step_solutions[pg->imtrx].xdot_old,
                       pg->sub_step_solutions[pg->imtrx].xdot_older);
                dcopy1(numProcUnknowns[pg->imtrx],
                       pg->sub_step_solutions[pg->imtrx].xdot,
                       pg->sub_step_solutions[pg->imtrx].xdot_old);
                dcopy1(numProcUnknowns[pg->imtrx],
                       pg->sub_step_solutions[pg->imtrx].x_older,
                       pg->sub_step_solutions[pg->imtrx].x_oldest);
                dcopy1(numProcUnknowns[pg->imtrx],
                       pg->sub_step_solutions[pg->imtrx].x_old,
                       pg->sub_step_solutions[pg->imtrx].x_older);
                dcopy1(numProcUnknowns[pg->imtrx],
                       pg->sub_step_solutions[pg->imtrx].x,
                       pg->sub_step_solutions[pg->imtrx].x_old);

                exchange_dof(cx[pg->imtrx], dpi,
                             pg->sub_step_solutions[pg->imtrx].x, pg->imtrx);

                exchange_dof(cx[pg->imtrx], dpi,
                             pg->sub_step_solutions[pg->imtrx].x_old,
                             pg->imtrx);

                exchange_dof(cx[pg->imtrx], dpi,
                             pg->sub_step_solutions[pg->imtrx].x_older,
                             pg->imtrx);

                exchange_dof(cx[pg->imtrx], dpi,
                             pg->sub_step_solutions[pg->imtrx].xdot, pg->imtrx);

                exchange_dof(cx[pg->imtrx], dpi,
                             pg->sub_step_solutions[pg->imtrx].xdot_old,
                             pg->imtrx);
              }
              dcopy1(numProcUnknowns[pg->imtrx],
                     pg->sub_step_solutions[pg->imtrx].x, x[pg->imtrx]);
              // update actual values
              for (i = 0; i < numProcUnknowns[pg->imtrx]; i++) {
                xdot[pg->imtrx][i] =
                    (1.0 + 2.0 * theta) / delta_t *
                        (x[pg->imtrx][i] - x_old[pg->imtrx][i]) -
                    (2.0 * theta) * xdot_old[pg->imtrx][i];
              }


            } else {
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
               * x_old, to the values of the external variables at that time
               * point.
               */
              if (efv->ev) {
                timeValueReadTrans = time;
                int w;
                for (w = 0; w < efv->Num_external_field; w++) {
                  if (strcmp(efv->field_type[w], "transient") == 0) {
                    err = rd_trans_vectors_from_exoII(
                        x_old[pg->imtrx], efv->file_nm[w], w, n,
                        &timeValueReadTrans, cx[pg->imtrx], dpi);
                    if (err != 0) {
                      DPRINTF(stderr,
                              "%s: err from rd_trans_vectors_from_exoII\n", yo);
                    }
                  }
                }
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
                fprintf(stdout,
                        "\n=> Try for soln at t=%g with dt=%g [%d for %d] %s\n",
                        time1, delta_t, nt, n, tspstring);
                log_msg("Predicting try at t=%g, dt=%g [%d for %d so far] %s",
                        time1, delta_t, nt, n, tspstring);
              }

              /*
               * Predict the solution, x[], and its derivative, xdot[],
               * at the new time, time1, using the old solution, xdot_old[],
               * And its derivatives at the old time, time.
               */
              if (subcycle == 0) {
                if (upd->SegregatedSolve && pg->imtrx == 0) {
                  predict_solution_u_star(numProcUnknowns[pg->imtrx], delta_t,
                                          delta_t_old, delta_t_older, theta, x,
                                          x_old, x_older, x_oldest);
                } else {
                  predict_solution(
                      numProcUnknowns[pg->imtrx], delta_t, delta_t_old,
                      delta_t_older, theta, x[pg->imtrx], x_old[pg->imtrx],
                      x_older[pg->imtrx], x_oldest[pg->imtrx], xdot[pg->imtrx],
                      xdot_old[pg->imtrx], xdot_older[pg->imtrx]);
                }
              }

              if (ls != NULL && ls->Evolution == LS_EVOLVE_SLAVE) {
                surf_based_initialization(x[pg->imtrx], NULL, NULL, exo,
                                          num_total_nodes, ls->init_surf_list,
                                          time1, theta, delta_t);
              }

              /*
               * Now, that we have a predicted solution for the current
               * time, x[], exchange the degrees of freedom to update the
               * ghost node information.
               */
              exchange_dof(cx[pg->imtrx], dpi, x[pg->imtrx], pg->imtrx);
              exchange_dof(cx[pg->imtrx], dpi, xdot[pg->imtrx], pg->imtrx);

              if (matrix_nAC[pg->imtrx] > 0 && subcycle == 0) {

                predict_solution(matrix_nAC[pg->imtrx], delta_t, delta_t_old,
                                 delta_t_older, theta, x_AC[pg->imtrx],
                                 x_AC_old[pg->imtrx], x_AC_older[pg->imtrx],
                                 x_AC_oldest[pg->imtrx], x_AC_dot[pg->imtrx],
                                 x_AC_dot_old[pg->imtrx],
                                 x_AC_dot_older[pg->imtrx]);

                for (iAC = 0; iAC < matrix_nAC[pg->imtrx]; iAC++) {
                  update_parameterAC(iAC, x[pg->imtrx], xdot[pg->imtrx],
                                     x_AC[pg->imtrx], cx[pg->imtrx], exo, dpi);
                  augc[iAC].tmp2 = x_AC_dot[pg->imtrx][iAC];
                  augc[iAC].tmp3 = x_AC_old[pg->imtrx][iAC];
                }
              }

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

              exchange_dof(cx[pg->imtrx], dpi, x[pg->imtrx], pg->imtrx);
              exchange_dof(cx[pg->imtrx], dpi, xdot[pg->imtrx], pg->imtrx);

              /*
               * Save the predicted solution for the time step
               * norm calculation to be carried out after convergence
               * of the nonlinear implicit problem
               */
              if (subcycle == 0) {
                dcopy1(numProcUnknowns[pg->imtrx], x[pg->imtrx],
                       x_pred[pg->imtrx]);
              }

              /*
               *  Solve the nonlinear problem. If we achieve convergence,
               *  set the flag, converged, to true on return. If not
               *  set the flag to false.

               */

              if (ProcID == 0 && upd->SegregatedSubcycles == 1 &&
                  renorm_subcycle_count <= 1) {
                printf("\n========================== SOLVING MATRIX %d "
                       "===========================\n\n",
                       pg->imtrx + 1);
              } else if (ProcID == 0 && (upd->SegregatedSubcycles > 1 ||
                                         renorm_subcycle_count > 1)) {
                int num_subcycles =
                    MAX(upd->SegregatedSubcycles, renorm_subcycle_count);
                printf(
                    "\n================== SOLVING MATRIX %d Subcycle %d/%d "
                    "===================\n\n",
                    pg->imtrx + 1, subcycle + 1, num_subcycles);
              }

              nAC = matrix_nAC[pg->imtrx];
              augc = matrix_augc[pg->imtrx];

              err = solve_nonlinear_problem(
                  ams[pg->imtrx], x[pg->imtrx], delta_t, theta,
                  x_old[pg->imtrx], x_older[pg->imtrx], xdot[pg->imtrx],
                  xdot_old[pg->imtrx], resid_vector[pg->imtrx],
                  x_update[pg->imtrx], scale[pg->imtrx], &converged, &nprint,
                  tev[pg->imtrx], tev_post[pg->imtrx], gv, rd[pg->imtrx], NULL,
                  NULL, gvec[pg->imtrx], gvec_elem[pg->imtrx], time1, exo, dpi,
                  cx[pg->imtrx], 0, &time_step_reform, 0, x_AC[pg->imtrx],
                  x_AC_dot[pg->imtrx], time1, NULL, NULL, NULL, NULL);
            } // sub-time loop if else

            if (pd_glob[0]->v[pg->imtrx][MOMENT0] ||
                pd_glob[0]->v[pg->imtrx][MOMENT1] ||
                pd_glob[0]->v[pg->imtrx][MOMENT2] ||
                pd_glob[0]->v[pg->imtrx][MOMENT3]) {
              /*     Floor values to 0 */
              int floored_values = 0;
              for (int var = MOMENT0; var <= MOMENT3; var++) {
                for (i = 0; i < num_total_nodes; i++) {
                  if (pd_glob[0]->v[pg->imtrx][var]) {
                    int j = Index_Solution(i, var, 0, 0, -1, pg->imtrx);

                    if (j != -1 && x[pg->imtrx][j] < 0) {
                      x[pg->imtrx][j] = 0.0;
                      floored_values++;
                    }
                  }
                }
              }

              int global_floored = 0;
              MPI_Allreduce(&floored_values, &global_floored, 1, MPI_INT,
                            MPI_SUM, MPI_COMM_WORLD);

              P0PRINTF("Floored %d moment values\n", global_floored);
            }

            /*
              err = solve_linear_segregated(ams[pg->imtrx], x[pg->imtrx],
              delta_t, theta, x_old[pg->imtrx], x_older[pg->imtrx],
              xdot[pg->imtrx], xdot_old[pg->imtrx], resid_vector[pg->imtrx],
              x_update[pg->imtrx], scale[pg->imtrx], &converged, &nprint, gv,
              time1, exo, dpi, cx, n, &time_step_reform);
            */
            if (err == -1) {
              converged = FALSE;
            }
            if (!converged) {
              /* Copy previous solution values if failed timestep */
              for (int imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
                dcopy1(numProcUnknowns[imtrx], x_old[imtrx], x[imtrx]);
              }
            }
            inewton = err;
            evpl_glob[0]->update_flag =
                0; /*See get_evp_stress_tensor for description */
            af->Sat_hyst_reevaluate =
                FALSE; /*See load_saturation for description*/

            if (converged) {
              for (i = 0;
                   matrix_nAC[pg->imtrx] > 0 && i < matrix_nAC[pg->imtrx];
                   i++) {
                gv[5 + invACidx[pg->imtrx][i]] = x_AC[pg->imtrx][i];
              }

              if (nAC > 0) {
                DPRINTF(stdout, "\n------------------------------\n");
                DPRINTF(stdout, "Augmenting Conditions:    %4d\n", nAC);
                DPRINTF(stdout, "Number of extra unknowns: %4d\n\n", nAC);

                for (iAC = 0; iAC < nAC; iAC++) {
                  if (augc[iAC].Type == AC_USERBC) {
                    DPRINTF(stdout, "\tBC[%4d] DF[%4d]=% 10.6e\n",
                            augc[iAC].BCID, augc[iAC].DFID,
                            x_AC[pg->imtrx][iAC]);
                    /* temporary printing */
#if 0
		    if( (int)augc[iAC].DataFlt[1] == 6)
		      {
			DPRINTF(stderr, "\tBC[%4d] DF[%4d]=% 10.6e\n", augc[iAC].DFID, 0, BC_Types[augc[iAC].DFID].BC_Data_Float[0]);
			DPRINTF(stderr, "\tBC[%4d] DF[%4d]=% 10.6e\n", augc[iAC].DFID, 2, BC_Types[augc[iAC].DFID].BC_Data_Float[2]);
			DPRINTF(stderr, "\tBC[%4d] DF[%4d]=% 10.6e\n", augc[iAC].DFID, 3, BC_Types[augc[iAC].DFID].BC_Data_Float[3]);
			augc[iAC].DataFlt[5] += augc[iAC].DataFlt[6];
			DPRINTF(stderr, "\tAC[%4d] DF[%4d]=% 10.6e\n", iAC, 5, augc[iAC].DataFlt[5]);

		      }
		    if( (int)augc[iAC].DataFlt[1] == 61)
		      {
			augc[iAC].DataFlt[5] += augc[iAC].DataFlt[6];
			DPRINTF(stderr, "\tAC[%4d] DF[%4d]=% 10.6e\n", iAC, 5, augc[iAC].DataFlt[5]);

		      }
#endif
                  }
                  /*	      else if (augc[iAC].Type == AC_USERMAT ||
                     augc[iAC].Type == AC_FLUX_MAT)
                                {
                                DPRINTF(stderr, "\tMT[%4d] MP[%4d]=% 10.6e\n",
                     augc[iAC].MTID, augc[iAC].MPID, x_AC[iAC]);
                                }
                                else if (augc[iAC].Type == AC_VOLUME)
                                {
                                evol_local = augc[iAC].evol;
                                #ifdef PARALLEL
                                if (Num_Proc > 1) {
                                MPI_Allreduce(&evol_local, &evol_global, 1,
                     MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); evol_local =
                     evol_global;
                                }
                                #endif
                                DPRINTF(stderr, "\tMT[%4d] VC[%4d]=%10.6e
                     Param=%10.6e\n", augc[iAC].MTID, augc[iAC].VOLID,
                     evol_local, x_AC[iAC]);
                                }
                                else if (augc[iAC].Type == AC_FLUX)
                                {
                                DPRINTF(stderr, "\tBC[%4d] DF[%4d]=%10.6e\n",
                     augc[iAC].BCID, augc[iAC].DFID, x_AC[iAC]);
                                }
                                else if (augc[iAC].Type == AC_POSITION)
                                {
                                DPRINTF(stderr, "\tNodeSet[%4d]_Pos = %10.6e
                     F_bal = %10.6e VC[%4d] Param=%10.6e\n", augc[iAC].MTID,
                     augc[iAC].evol, augc[iAC].lm_resid, augc[iAC].VOLID,
                     x_AC[iAC]);
                                } */
                }
              }
            }

            /*
             * HKM -> I do not know if these operations are needed. I added
             *        an exchange of xdot[] here, because if x[] is exchanged
             *        then xdot needs to be exchanged as well.
             */

            exchange_dof(cx[pg->imtrx], dpi, x[pg->imtrx], pg->imtrx);
            exchange_dof(cx[pg->imtrx], dpi, xdot[pg->imtrx], pg->imtrx);

            if (!converged)
              goto finish_step;
          }
        } // subcycle loop

        // reset renorm subcycle after we have completed the appropriate number
        // of subcycles
        renorm_subcycle_count = 0;

      finish_step:

        if (converged)
          af->Sat_hyst_reevaluate = TRUE; /*see load_saturation */

        /* Check element quality */
        /*good_mesh = element_quality(exo, x[pg->imtrx], ams[0]->proc_config);
         */
        good_mesh = 1;

        /*
         * Check the time step truncation error.
         * If too large, set the success_dt flag to FALSE. We will
         * then not accept the current time step, reduced delta_t,
         * and retry.
         */
        if (converged) {
          int num_success = 0;
          /* Assume we want the minimum delta_t_new */
          delta_t_new = 1e20;
          for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices;
               pg->imtrx++) {
            double mat_dt_new = 1e20;

            if (pg->time_step_control_disabled[pg->imtrx]) {
              success_dt = 1;
            } else {
              nAC = matrix_nAC[pg->imtrx];
              augc = matrix_augc[pg->imtrx];
              mat_dt_new = time_step_control(
                  delta_t, delta_t_old, const_delta_t, x[pg->imtrx],
                  x_pred[pg->imtrx], x_old[pg->imtrx], x_AC[pg->imtrx],
                  x_AC_pred[pg->imtrx], eps, &success_dt, tran->use_var_norm);
            }

            if (upd->SegregatedSolve && (pg->imtrx == 1 || pg->imtrx == 3)) {
              success_dt = 1;
            }

            num_success += success_dt ? 1 : 0;
            if (upd->SegregatedSolve) {
              if (pg->imtrx == 0) {
                delta_t_new = mat_dt_new;
              }
            } else {
              delta_t_new = MIN(mat_dt_new, delta_t_new);
            }
          }

          if (num_success != upd->Total_Num_Matrices) {
            success_dt = 0;
          }

          if (const_delta_t) {
            success_dt = TRUE;
            delta_t_new = delta_t;
          } else if (failed_recently_countdown > 0) {
            delta_t_new = delta_t;
            failed_recently_countdown--;
          } else if (delta_t_new > delta_t_max) {
            delta_t_new = delta_t_max;
            /*        } else if ( !success_dt && delta_t_new <
             * tran->resolved_delta_t_min ) {*/
          } else if (delta_t_new < tran->resolved_delta_t_min) {
            /*          if ( delta_t > tran->resolved_delta_t_min ) {  */
            /* fool algorithm into using delta_t = tran->resolved_delta_t_min */
            delta_t_new = tran->resolved_delta_t_min;
            success_dt = TRUE;
            DPRINTF(stdout, "\n\tminimum resolved step limit!\n");
          }

          pg->imtrx = Fill_Matrix;
          if (ls != NULL && tran->Courant_Limit != 0.) {
            double Courant_dt;
            Courant_dt =
                tran->Courant_Limit * pg->matrix_subcycle_count[Fill_Matrix] *
                Courant_Time_Step(x[pg->imtrx], x_old[pg->imtrx],
                                  x_older[pg->imtrx], xdot[pg->imtrx],
                                  xdot_old[pg->imtrx], resid_vector[pg->imtrx],
                                  ams[pg->imtrx]->proc_config, exo);
            if (Courant_dt > 0. && Courant_dt < delta_t_new) {
              DPRINTF(stdout, "\nCourant Limit requires dt <= %g\n",
                      Courant_dt);
              delta_t_new = Courant_dt;
            }
          }
        }

        if (converged && success_dt) {
          nt += 1;
          time = time1;

          /* Determine whether to print out the data or not */
          i_print = 0;
          if (tran->print_freq == 0) {
            if ((time > time_print) ||
                (fabs(time - time_print) < (1.e-4 * tran->print_delt))) {
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

          pg->imtrx = 0;
          for (i = 0; i < nn_post_fluxes; i++) {
            (void)evaluate_flux(
                exo, dpi, pp_fluxes[i]->ss_id, pp_fluxes[i]->flux_type,
                pp_fluxes[i]->flux_type_name, pp_fluxes[i]->blk_id,
                pp_fluxes[i]->species_number, pp_fluxes[i]->flux_filenm,
                pp_fluxes[i]->profile_flag, x[pg->imtrx], xdot[pg->imtrx], NULL,
                delta_t_old, time, 1);
          }

          for (i = 0; i < nn_volume; i++) {
            evaluate_volume_integral(
                exo, dpi, pp_volume[i]->volume_type, pp_volume[i]->volume_name,
                pp_volume[i]->blk_id, pp_volume[i]->species_no,
                pp_volume[i]->volume_fname, pp_volume[i]->params,
                pp_volume[i]->num_params, NULL, x[pg->imtrx], xdot[pg->imtrx],
                delta_t, time1, 1);
          }

          if (time1 >= (ROUND_TO_ONE * TimeMax))
            i_print = 1;

          error = 0;
          if (i_print) {
            if (Write_Intermediate_Solutions == 0) {
              write_solution_segregated(
                  ExoFileOut, resid_vector, x, x_old, xdot, xdot_old, tev_post,
                  gv, rd, gvec, &nprint, delta_t, theta, time1, NULL, exo, dpi);
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
            if ((time + 1.2 * delta_t_new >= time_print) &&
                (time_print > time)) {
              delta_t_new = time_print - time;
              DPRINTF(stdout,
                      "reset delta_t = %g to maintain printing frequency\n",
                      delta_t_new);
              if (delta_t_new <= 0)
                EH(GOMA_ERROR, "error with time-step printing control");
            } else if (time >= time_print) {
              if (delta_t_new != tran->print_delt) {
                delta_t_new = tran->print_delt;
                DPRINTF(stdout,
                        "reset delta_t = %g to maintain printing frequency\n",
                        delta_t_new);
                if (delta_t_new <= 0) {
                  EH(GOMA_ERROR, "error with time-step printing control");
                }
              }
            }
          }

          /* Check if steady state has been reached */
          if (tran->march_to_steady_state && n > 0) {
            /* for now check last two matrices */
            int steady_state_reached = TRUE;
            double max_distance = 0;

            for (i = 0; i < upd->Total_Num_Matrices; i++) {
              double distance = vector_distance(NumUnknowns[i], x[i], x_old[i]);
              if (distance > tran->steady_state_tolerance) {
                steady_state_reached = FALSE;
              }
              if (distance > max_distance)
                max_distance = distance;
            }

            if (ProcID == 0) {
              printf("\nMaximum Delta x %g\n", max_distance);
            }

            if (steady_state_reached) {
              if (ProcID == 0) {
                printf("\n Steady state reached \n");
              }
              goto free_and_clear;
            }
          }

          delta_t_oldest = delta_t_older;
          delta_t_older = delta_t_old;
          delta_t_old = delta_t;
          tran->delta_t_old = delta_t_old;
          tran->time_value_old = time;
          delta_t = delta_t_new;
          tran->delta_t = delta_t; /*load up for use in load_fv_mesh_derivs*/
          tran->delta_t_avg =
              0.25 * (delta_t + delta_t_old + delta_t_older + delta_t_oldest);

          if (time1 >= (ROUND_TO_ONE * TimeMax)) {
            DPRINTF(stdout, "\t\tout of time!\n");
            if (Anneal_Mesh) {
              /*
               * Transform the node point coordinates according to the
               * displacements and write out all the results using the
               * displaced coordinates. Set the displacement field to
               * zero, too.
               */
              for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices;
                   pg->imtrx++) {

                err = anneal_mesh(x[pg->imtrx], tev[pg->imtrx],
                                  tev_post[pg->imtrx], gv, rd[pg->imtrx], time1,
                                  exo, dpi);
                EH(err, "anneal_mesh() bad return.");
              }
            }
            goto free_and_clear;
          }
          if (!good_mesh)
            goto free_and_clear;

          if (converged && ls != NULL) {
            pg->imtrx = Fill_Matrix;
            int ibc, ls_adc_event;
            /* Resolve LS_ADC boundaries ( Attach/Dewet/Coalesce ) */
            for (ibc = 0; ibc < Num_BC; ibc++) {
              ls_adc_event = FALSE;
              switch (BC_Types[ibc].BC_Name) {
              case LS_ADC_OLD_BC:
                resolve_ls_adc_old(&(BC_Types[ibc]), exo, x[pg->imtrx], delta_t,
                                   &ls_adc_event, nt);
                break;
              case LS_ADC_BC:
                resolve_ls_adc(ls->last_surf_list->start->subsurf_list,
                               &(BC_Types[ibc]), exo, x[pg->imtrx], delta_t,
                               &ls_adc_event, nt);

                break;
              default:
                break;
              }
#ifdef PARALLEL
              if (ls_adc_event) {
                exchange_dof(cx[pg->imtrx], dpi, x[pg->imtrx], 0);
              }
#endif
              if (ls_adc_event && tran->Restart_Time_Integ_After_Renorm) {
                /* like a restart */
                for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices;
                     pg->imtrx++) {
                  discard_previous_time_step(
                      numProcUnknowns[pg->imtrx], x[pg->imtrx],
                      x_old[pg->imtrx], x_older[pg->imtrx], x_oldest[pg->imtrx],
                      xdot[pg->imtrx], xdot_old[pg->imtrx],
                      xdot_older[pg->imtrx]);
                }

                last_renorm_nt = nt;
                if (delta_t_new > fabs(delta_t0))
                  delta_t_new *= tran->time_step_decelerator;

                pg->imtrx = 0;
              }
            }

            /* Check for renormalization  */

            ls->Renorm_Countdown -= 1;
            switch (ls->Renorm_Method) {

            case HUYGENS:
            case HUYGENS_C:
            case HUYGENS_MASS_ITER:
              Renorm_Now =
                  (ls->Renorm_Freq != 0 && ls->Renorm_Countdown == 0) ||
                  ls_adc_event == TRUE;

              pg->imtrx = Fill_Matrix;
              did_renorm = huygens_renormalization(
                  x[pg->imtrx], num_total_nodes, exo, cx[pg->imtrx], dpi,
                  num_fill_unknowns, numProcUnknowns[pg->imtrx], time2,
                  Renorm_Now);
              if (did_renorm) {
                exchange_dof(cx[pg->imtrx], dpi, x[pg->imtrx], pg->imtrx);

                renorm_subcycle_count = ls->SubcyclesAfterRenorm;
              }
              if (did_renorm && tran->Restart_Time_Integ_After_Renorm) {
                /* like a restart */
                for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices;
                     pg->imtrx++) {
                  discard_previous_time_step(
                      numProcUnknowns[pg->imtrx], x[pg->imtrx],
                      x_old[pg->imtrx], x_older[pg->imtrx], x_oldest[pg->imtrx],
                      xdot[pg->imtrx], xdot_old[pg->imtrx],
                      xdot_older[pg->imtrx]);
                }
                last_renorm_nt = nt;
                if (delta_t_new > fabs(delta_t0))
                  delta_t_new *= tran->time_step_decelerator;
              }
              pg->imtrx = 0;
              break;

            case CORRECT:
              EH(GOMA_ERROR, "Use of \"CORRECT\" is obsolete.");
              break;
            default:
              break;
            }

            if (ls->Sat_Hyst_Renorm_Lockout > 0) {
              af->Sat_hyst_reevaluate = FALSE;
              ls->Sat_Hyst_Renorm_Lockout -= 1;
            }
          }

          /*
           *   save xdot to xdot_old for next time step
           */
          for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices;
               pg->imtrx++) {
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

            if (matrix_nAC[pg->imtrx] > 0) {
              dcopy1(matrix_nAC[pg->imtrx], x_AC_dot_old[pg->imtrx],
                     x_AC_dot_older[pg->imtrx]);
              dcopy1(matrix_nAC[pg->imtrx], x_AC_dot[pg->imtrx],
                     x_AC_dot_old[pg->imtrx]);
              dcopy1(matrix_nAC[pg->imtrx], x_AC_older[pg->imtrx],
                     x_AC_oldest[pg->imtrx]);
              dcopy1(matrix_nAC[pg->imtrx], x_AC_old[pg->imtrx],
                     x_AC_older[pg->imtrx]);
              dcopy1(matrix_nAC[pg->imtrx], x_AC[pg->imtrx],
                     x_AC_old[pg->imtrx]);
            }
          }

        }    /*  if(converged && success_dt) */
        else /* not converged or unsuccessful time step */
        {
          /* Set bit TRUE in next line to enable retries for failed first
           * timestep*/
          if (relax_bit && nt == 0 && n < 15) {
            if (inewton == -1) {
              DPRINTF(stdout,
                      "\nHmm... trouble on first step \n  Let's try some "
                      "more relaxation  \n");
              if ((damp_factor1 <= 1. && damp_factor1 >= 0.) &&
                  (damp_factor2 <= 1. && damp_factor2 >= 0.) &&
                  (damp_factor3 <= 1. && damp_factor3 >= 0.)) {
                custom_tol1 *= 0.01;
                custom_tol2 *= 0.01;
                custom_tol3 *= 0.01;
                DPRINTF(stdout, "  custom tolerances %g %g %g  \n", custom_tol1,
                        custom_tol2, custom_tol3);
              } else {
                damp_factor1 *= 0.5;
                DPRINTF(stdout, "  damping factor %g  \n", damp_factor1);
              }
            } else {
              DPRINTF(stdout,
                      "\nHmm... could not converge on first step\n Let's "
                      "try some more iterations\n");
              dcopy1(numProcUnknowns[pg->imtrx], x[pg->imtrx],
                     x_old[pg->imtrx]);

              if ((damp_factor1 <= 1. && damp_factor1 >= 0.) &&
                  (damp_factor2 <= 1. && damp_factor2 >= 0.) &&
                  (damp_factor3 <= 1. && damp_factor3 >= 0.)) {
                custom_tol1 *= 100.;
                custom_tol2 *= 100.;
                custom_tol3 *= 100.;
                DPRINTF(stdout, "  custom tolerances %g %g %g  \n", custom_tol1,
                        custom_tol2, custom_tol3);
              } else {
                damp_factor1 *= 2.0;
                damp_factor1 = MIN(damp_factor1, 1.0);
                DPRINTF(stdout, "  damping factor %g  \n", damp_factor1);
              }
            }
          } else if (delta_t <
                     tran->resolved_delta_t_min / tran->time_step_decelerator) {
            DPRINTF(stdout, "\n\tminimum resolved step limit!\n");
            delta_t_oldest = delta_t_older;
            delta_t_older = delta_t_old;
            delta_t_old = delta_t;
            tran->delta_t_old = delta_t_old;
            tran->time_value_old = time;
            delta_t = tran->resolved_delta_t_min;
            tran->delta_t = delta_t;
            tran->delta_t_avg =
                0.25 * (delta_t + delta_t_old + delta_t_older + delta_t_oldest);
            time1 = time + delta_t;
            tran->time_value = time1;
          } else {
            DPRINTF(stdout,
                    "\n\tlast time step failed, dt *= %g for next try!\n",
                    tran->time_step_decelerator);

            delta_t *= tran->time_step_decelerator;
            tran->delta_t = delta_t;
            tran->delta_t_avg =
                0.25 * (delta_t + delta_t_old + delta_t_older + delta_t_oldest);
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

        if (delta_t <= delta_t_min) {
          DPRINTF(stdout, "\n\tdelta_t = %e < %e\n\n", delta_t, delta_t_min);

          DPRINTF(stdout, "time step too small, I'm giving up!\n");
          break;
        }
      } /* end of time step loop */
    }   /* end of if steady else transient */
  free_and_clear:

    if (timestep_subcycle) {
      for (int imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
        if (pg->matrix_subcycle_count[imtrx] > 1) {
          free(pg->sub_step_solutions[imtrx].x);
          free(pg->sub_step_solutions[imtrx].x_old);
          free(pg->sub_step_solutions[imtrx].x_older);
          free(pg->sub_step_solutions[imtrx].x_oldest);
          free(pg->sub_step_solutions[imtrx].xdot);
          free(pg->sub_step_solutions[imtrx].xdot_old);
          free(pg->sub_step_solutions[imtrx].xdot_older);
          free(pg->sub_step_solutions[imtrx].x_update);
        }
      }
      free(pg->sub_step_solutions);
    }

    for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
      free(x[pg->imtrx]);
      free(x_old[pg->imtrx]);
      free(x_older[pg->imtrx]);
      free(x_oldest[pg->imtrx]);
      free(xdot[pg->imtrx]);
      free(xdot_old[pg->imtrx]);
      free(xdot_older[pg->imtrx]);
      free(x_pred[pg->imtrx]);
      free(delta_x[pg->imtrx]);
      free(x_previous[pg->imtrx]);
      free(x_update[pg->imtrx]);
    }

    free(x);
    free(x_old);
    free(x_older);
    free(x_oldest);
    free(xdot);
    free(xdot_old);
    free(xdot_older);
    free(x_pred);
    free(delta_x);
    free(x_previous);
    free(x_update);

    free(timeValueRead);
    for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
      free(rd[pg->imtrx]);
      if ((tev[pg->imtrx] + tev_post[pg->imtrx]) > 0) {
        for (i = 0; i < exo->num_elem_blocks; i++) {
          free(gvec_elem[pg->imtrx][i]);
        }
      }
      free(gvec_elem[pg->imtrx]);
    }
    free(rd);
    free(gvec_elem);
    free(tnv);
    free(tev);
    free(tnv_post);
    free(tev_post);

    for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
      free(resid_vector[pg->imtrx]);
      free(resid_vector_sens[pg->imtrx]);
      free(scale[pg->imtrx]);
    }

    free(resid_vector);
    free(resid_vector_sens);
    free(scale);

    for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
      free(gvec[pg->imtrx]);
      free(matrix_augc[pg->imtrx]);
    }
    free(matrix_augc);
    free(gvec);

    free(numProcUnknowns);

    free(ija);
    free(ija_attic);
    free(a);
    free(a_old);

    if (totalnAC > 0) {
      for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
        if (matrix_nAC[pg->imtrx] > 0) {
          free(x_AC[pg->imtrx]);
          free(x_AC_old[pg->imtrx]);
          free(x_AC_older[pg->imtrx]);
          free(x_AC_oldest[pg->imtrx]);
          free(x_AC_dot[pg->imtrx]);
          free(x_AC_dot_old[pg->imtrx]);
          free(x_AC_dot_older[pg->imtrx]);
          free(x_AC_pred[pg->imtrx]);
        }
      }
    }
    free(matrix_nAC);
    free(x_AC);
    free(x_AC_old);
    free(x_AC_older);
    free(x_AC_oldest);
    free(x_AC_dot);
    free(x_AC_dot_old);
    free(x_AC_dot_older);
    free(x_AC_pred);
    free(gv);

    return;
  }

  void predict_solution_u_star(int N, dbl delta_t, dbl delta_t_old,
                               dbl delta_t_older, dbl theta_arg, dbl **x,
                               dbl **x_old, dbl **x_older, dbl **x_oldest) {
    int i;
    dbl c1, c2;

    c1 = delta_t * (1.0 + theta_arg * delta_t / delta_t_old);
    c2 = theta_arg * (delta_t * delta_t) / (delta_t_old);
    for (i = 0; i < N; i++) {
      x[0][i] = x_old[2][i] + c1 * (x_old[0][i] - x_older[2][i]) / delta_t_old -
                c2 * (x_older[0][i] - x_oldest[2][i]) / delta_t_older;
    }
  }

  static int discard_previous_time_step(
      int num_unks, double *x, double *x_old, double *x_older, double *x_oldest,
      double *xdot, double *xdot_old, double *xdot_older) {

    dcopy1(num_unks, x, x_old);
    dcopy1(num_unks, x_old, x_older);
    dcopy1(num_unks, x_older, x_oldest);

    /* also need to kill xdot(s) */
    memset(xdot, 0, sizeof(double) * num_unks);
    memset(xdot_old, 0, sizeof(double) * num_unks);
    memset(xdot_older, 0, sizeof(double) * num_unks);

    return (0);
  }
