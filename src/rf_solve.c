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
 *$Id: rf_solve.c,v 5.21 2010-03-17 22:23:54 hkmoffa Exp $
 */

/*
 * Revision history has goneaway from this place.
 */

#include "rf_solve.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ac_particles.h"
#include "ac_stability_util.h"
#include "az_aztec.h"
#include "brkfix/fix.h"
#include "decomp_interface.h"
#include "dp_comm.h"
#include "dp_types.h"
#include "dp_utils.h"
#include "dpi.h"
#include "el_elm.h"
#include "el_elm_info.h"
#include "el_geom.h"
#include "exo_struct.h"
#include "linalg/sparse_matrix.h"
#include "mm_as.h"
#include "mm_as_structs.h"
#include "mm_augc_util.h"
#include "mm_bc.h"
#include "mm_eh.h"
#include "mm_fill_ls.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_stress.h"
#include "mm_fill_util.h"
#include "mm_flux.h"
#include "mm_more_utils.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "mm_post_def.h"
#include "mm_post_proc.h"
#include "mm_sol_nonlinear.h"
#include "mm_unknown_map.h"
#include "mm_viscosity.h"
#include "mpi.h"
#include "polymer_time_const.h"
#include "rd_exo.h"
#include "rd_mesh.h"
#include "rf_allo.h"
#include "rf_bc.h"
#include "rf_bc_const.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_io.h"
#include "rf_io_const.h"
#include "rf_io_structs.h"
#include "rf_masks.h"
#include "rf_mp.h"
#include "rf_node_const.h"
#include "rf_solve_segregated.h"
#include "rf_solver.h"
#include "rf_util.h"
#include "sl_auxutil.h"
#include "sl_matrix_util.h"
#include "sl_petsc.h"
#include "sl_petsc_complex.h"
#include "sl_util.h"
#include "sl_util_structs.h"
#include "std.h"
#include "usr_print.h"
#include "wr_dpi.h"
#include "wr_exo.h"
#include "wr_soln.h"

#define GOMA_RF_SOLVE_C
#include "el_quality.h"

#ifdef GOMA_ENABLE_OMEGA_H
#include "adapt/omega_h_interface.h"
#endif

/*
 * Global variables defined in this file.
 */

extern struct elem_side_bc_struct ***First_Elem_Side_BC_Array;
extern struct elem_edge_bc_struct ***First_Elem_Edge_BC_Array;

#define ROUND_TO_ONE 0.9999999

/*
 * Declarations of static functions defined in this file.
 */

static void predict_solution_newmark(int,       /* N */
                                     double,    /* delta_t */
                                     double[],  /* x */
                                     double[],  /* x_old */
                                     double[],  /* xdot */
                                     double[]); /* xdot_old */

int discard_previous_time_step(
    int, double *, double *, double *, double *, double *, double *, double *);

static void shift_nodal_values(int, double, double *, int);

extern FSUB_TYPE dsyev_(char *JOBZ,
                        char *UPLO,
                        int *N,
                        double *A,
                        int *LDA,
                        double *W,
                        double *WORK,
                        int *LWORK,
                        int *INFO,
                        int len_jobz,
                        int len_uplo);

extern FSUB_TYPE dsysv_(char *JOBZ,
                        int *N,
                        int *N_RHS,
                        double *A,
                        int *LDA,
                        int *PVT,
                        double *B,
                        int *LDB,
                        double *WORK,
                        int *LWORK,
                        int *INFO,
                        int len_jobz);
// C = A X B
void slow_square_dgemm(
    int transpose_b, int N, double A[DIM][DIM], double B[DIM][DIM], double C[DIM][DIM]) {
  int i, j, k;
  double tmp;
  if (transpose_b) {
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
        tmp = 0;
        for (k = 0; k < N; k++) {
          tmp += A[i][k] * B[j][k];
        }
        C[i][j] = tmp;
      }
    }
  } else {
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
        tmp = 0;
        for (k = 0; k < N; k++) {
          tmp += A[i][k] * B[k][j];
        }
        C[i][j] = tmp;
      }
    }
  }
}

void initial_guess_stress_to_log_conf(double *x, int num_total_nodes) {
  int a, b;

  double s[DIM][DIM];
  double log_s[DIM][DIM];
  int s_idx[2][2];
  int N = VIM;
  int LDA = N;
  int node, v, i, j;

  int INFO;
  int LWORK = 20;
  double WORK[LWORK];
  double A[DIM * DIM];
  dbl gamma_dot[DIM][DIM];

  int mode, mn;
  double lambda = 0;
  double mup;
  int v_s[MAX_MODES][DIM][DIM];

  stress_eqn_pointer(v_s);

  for (mn = 0; mn < upd->Num_Mat; mn++) // mn = Material Number
  {
    for (mode = 0; mode < vn_glob[mn]->modes; mode++) {
      ve[mode] = ve_glob[mn][mode];

      for (node = 0; node < num_total_nodes; node++) {
        memset(WORK, 0, sizeof(double) * LWORK);
        memset(A, 0.0, sizeof(double) * DIM * DIM);

        for (a = 0; a < 2; a++) {
          for (b = 0; b < 2; b++) {
            v = v_s[mode][a][b];
            s_idx[a][b] = Index_Solution(node, v, 0, 0, -2, pg->imtrx);
          }
        }

        mup = viscosity(ve[mode]->gn, gamma_dot, NULL);
        lambda = polymer_time_const(ve[mode]->time_const_st, gamma_dot, NULL);

        // skip node if stress variables not found
        if (s_idx[0][0] == -1 || s_idx[0][1] == -1 || s_idx[1][1] == -1)
          continue;

        // get stress tensor from initial guess
        s[0][0] = x[s_idx[0][0]];
        s[0][1] = x[s_idx[0][1]];
        s[1][0] = x[s_idx[0][1]];
        s[1][1] = x[s_idx[1][1]];

        // Convert stress to c
        for (i = 0; i < DIM; i++) {
          for (j = 0; j < DIM; j++) {
            s[i][j] = (lambda / mup) * s[i][j] + (double)delta(i, j);
          }
        }

        // convert to column major
        for (i = 0; i < VIM; i++) {
          for (j = 0; j < VIM; j++) {
            A[i * VIM + j] = s[j][i];
          }
        }

        double W[DIM];

        // eig solver
        dsyev_("V", "U", &N, A, &LDA, W, WORK, &LWORK, &INFO, 1, 1);

        double U[DIM][DIM];

        // transpose (revert to row major)
        for (i = 0; i < VIM; i++) {
          for (j = 0; j < VIM; j++) {
            U[i][j] = A[j * VIM + i];
          }
        }

        // Take log of diagonal
        double D[DIM][DIM];
        for (i = 0; i < VIM; i++) {
          for (j = 0; j < VIM; j++) {
            if (i == j) {
              D[i][j] = log(W[i]);
            } else {
              D[i][j] = 0.0;
            }
          }
        }

        /* matrix multiplication, the slow way */
        slow_square_dgemm(0, VIM, U, D, log_s);

        // multiply by transpose
        slow_square_dgemm(1, VIM, log_s, U, D);

        for (i = 0; i < VIM; i++) {
          for (j = 0; j < VIM; j++) {
            log_s[i][j] = D[i][j];
          }
        }

        x[s_idx[0][0]] = log_s[0][0];
        x[s_idx[0][1]] = log_s[0][1];
        x[s_idx[1][1]] = log_s[1][1];

      } /* Loop over nodes */
    }   /* Loop over modes */
  }     /* Loop over materials */
}

void solve_problem(Exo_DB *exo, /* ptr to the finite element mesh database  */
                   Dpi *dpi,    /* distributed processing information       */
                   dbl *te_out) /* te_out - return actual end time */

/*
    Routine that controls solution of overall FEM reacting
    flow problem. This routine will control time integration,
    matrix fills and nonlinear solver and results output.
*/

{
  /*
   * Sparse matrix storage vectors
   * (MSR format.  See "SPARSKIT: a basic tool kit for sparse matrix
   * computations" by Youcef Saad)
   */
  int *ija = NULL;              /* column pointer array                     */
  double *a = NULL;             /* nonzero array                            */
  double *a_old = NULL;         /* nonzero array                            */
  static double *x = NULL;      /* solution vector                          */
  static double *x_save = NULL; /* solution vector for reset */

  int iAC;                              /* Counter                                  */
  static double *x_AC = NULL;           /* Solution vector of extra unknowns          */
  static double *x_AC_old = NULL;       /* old solution vector of extra unknowns      */
  static double *x_AC_older = NULL;     /* older solution vector of extra unknowns    */
  static double *x_AC_oldest = NULL;    /* oldest solution vector of extra unknowns   */
  static double *x_AC_dot = NULL;       /* current time derivative of extra unknowns  */
  static double *x_AC_dot_old = NULL;   /* old time derivative of extra unknowns      */
  static double *x_AC_dot_older = NULL; /* Older time derivative of extra unknowns    */
  static double *x_AC_pred = NULL;      /* predicted extraunknowns */

  int *ija_attic = NULL; /* storage for external dofs                  */

  int eb_indx, ev_indx;

  /*
   * Variables
   */
  double *x_pred = NULL;            /* prediction of solution vector     */
  static double *x_old = NULL;      /* old solution vector               */
  static double *x_older = NULL;    /* older solution vector             */
  static double *x_oldest = NULL;   /* oldest solution vector saved      */
  static double *xdot = NULL;       /* current time derivative of soln   */
  static double *xdot_save = NULL;  /* current time derivative of soln for reset  */
  static double *xdot_old = NULL;   /* old time derivative of soln       */
  static double *xdot_older = NULL; /* old time derivative of soln       */

  double *x_sens = NULL;    /* solution sensitivity                     */
  double **x_sens_p = NULL; /* solution sensitivity for parameters      */
  int num_pvector = 0;      /* number of solution sensitivity vectors   */
#ifdef GOMA_ENABLE_OMEGA_H
  int adapt_step = 0;
#endif
  int last_adapt_nt = 0;

  /* sparse variables for fill equation subcycling */

  static struct GomaLinearSolverData *ams[NUM_ALSS] = {NULL};

  /* "sl_util_structs.h" */

  double *x_update = NULL; /* update at last iteration                 */

  double *resid_vector = NULL;      /* residual                                 */
  double *resid_vector_sens = NULL; /* residual sensitivity                  */

  double *scale = NULL; /* scale vector for modified newton         */

  int *node_to_fill = NULL;

  char tspstring[MAX_FNL]; /* literal representation of time step
                            * parameter, theta [0=BE,.5=CN,1=FE] but
                            * any float possible...                    */

  int n;                 /* total number of time steps attempted     */
  int nt;                /* total number of successful time steps    */
  int last_renorm_nt;    /* time step at which last renorm occured   */
  int time_step_reform;  /* counter for jacobian reformation stride  */
  int converged = TRUE;  /* success or failure of Newton iteration   */
  int success_dt = TRUE; /* success or failure of time step          */
  int failed_recently_countdown = 0;
  int i, num_total_nodes;
  int numProcUnknowns;
  int const_delta_t, const_delta_ts, step_print;
  int step_fix = 0, i_fix = 0; /* What step to fix the problem on */
  int good_mesh = TRUE;
  int w; /* counter for looping external variables */
  static int nprint = 0;
  double time_print, i_print;
  double theta = 0.0, time;
  static double time1 =
      0.0; /* Current time that the simulation is trying  to find the solution for */
#ifdef LIBRARY_MODE
  static double delta_t_save = 0.0;
#endif
  double delta_t, delta_t_new = 0.0;
  double delta_t_old, delta_t_older, delta_t_oldest = 0.0;
  double timeValueRead = 0.0; /* time value read from an exodus input file
                               * used to initialize the solution vector */
  double timeValueReadTrans = 0.0;
  //  static double time_value = 0.0;
  double eps;

  double time2 = 0.0;
  /*
   * Other local variables...
   */

  int error, err, is_steady_state, inewton;
  int *gindex = NULL, gsize;
  int *p_gsize;
  static double *gvec = NULL;
  static double ***gvec_elem = NULL;
  FILE *file = NULL;
  static struct Results_Description *rd;
  struct Level_Set_Data *ls_old;

  int tnv;      /* total number of nodal variables and kinds */
  int tev;      /* total number of elem variables and kinds  */
  int tnv_post; /* total number of nodal variables and kinds
                   for post processing                       */
  int tev_post; /* total number of elem variables and kinds
                   for post processing                       */

  double *gv;                   /* Global variable values */
  double *x_pp = NULL;          /* Post-proc variables for export */
  static int *xp_id = NULL;     /* Post-proc variable ID */
  double *base_p_por = NULL;    /* Base values for porosity updates */
  double *base_p_liq = NULL;    /* Base values for porosity updates */
  int update_porosity = FALSE;  /* Flag for external porosity updates */
  int update_etch_area = FALSE; /* Flag for etch area fraction updates */
#ifdef HAVE_FRONT
  int max_unk_elem, one, three; /* variables used as mf_setup arguments      */
#endif
  unsigned int matrix_systems_mask;
  int did_renorm;         /* Flag indicating if we renormalized.       */
  int Renorm_Now = FALSE; /* Flag forcing renormalization regardless of gradient */
  double evol_local = 0.0;
  double lsvel_local = 0.0;
#ifdef PARALLEL
  double evol_global = 0.0;
  double lsvel_global = 0.0;
#endif /* PARALLEL */

  static int callnum = 1; /* solve_problem call counter */
  int last_call = TRUE;   /* Indicates final rf_solve call */
#ifdef LIBRARY_MODE
  int last_step = FALSE; /* Indicates final time step on this call */
#endif
  double damp_factor_org[2] = {damp_factor1, damp_factor2};
#ifdef RELAX_ON_TRANSIENT_PLEASE
  int relax_bit = TRUE; /* Enables relaxation after a transient convergence failure*/
  double toler_org[3] = {custom_tol1, custom_tol2, custom_tol3};
#else
  int relax_bit = FALSE;
#endif
  int no_relax_retry = ceil(2. - log(damp_factor_org[0]) / log(2.));
  int nonconv_roll = 0;
  int use_custom_damp =
      ((damp_factor1 <= 1. && damp_factor1 >= 0.) && (damp_factor2 <= 1. && damp_factor2 >= 0.) &&
       (damp_factor3 <= 1. && damp_factor3 >= 0.));

  static const char yo[] = "solve_problem"; /* So my name is in a string.        */

  /*
   * 		BEGIN EXECUTION
   */

  /* Set step_fix only if parallel run and only if fix freq is enabled*/
  if (Num_Proc > 1 && tran->fix_freq > 0) {
    step_fix = 1; /* Always fix on the first timestep to match print frequency */
  }

  tran->time_value = time1;

  is_steady_state = (TimeIntegration == STEADY) ? TRUE : FALSE;

  p_gsize = &gsize;

#ifdef LIBRARY_MODE
  fprintf(stdout, "  Commencing call #%3d from ANIMAS to solve_problem\n", callnum);
  if (libio->animas_step != -1)
    last_call = FALSE;

  /* Determine if external porosity updates are required */
  if (Num_Var_In_Type[pg->imtrx][POR_LIQ_PRES] && efv->ev_porous_index > -1) {
    update_porosity = TRUE;
    fprintf(stdout, " External porosity field %d will be updated.\n", efv->ev_porous_index);
  }
  if (libio->goma_first == 1)
    fprintf(stdout, "  Goma goes first");
  if (libio->goma_first == 0)
    fprintf(stdout, "  Goma goes second");

#endif

  /* Determine if external area fraction updates are required */
  if ((Num_Var_In_Type[pg->imtrx][MASS_FRACTION]) && (efv->ev_etch_area > -1) &&
      (efv->ev_etch_depth > -1)) {
    update_etch_area = TRUE;
    fprintf(stderr, " External fields etch area %d and etch depth %d will be updated.\n",
            efv->ev_etch_area, efv->ev_etch_depth);
  }

  if (Unlimited_Output && strlen(Soln_OutFile)) {
    file = fopen(Soln_OutFile, "w");
    if (file == NULL) {
      fprintf(stdout, "%s:  opening soln file, %s, for writing\n", yo, Soln_OutFile);
      GOMA_EH(GOMA_ERROR, "Can not open solution file\n");
    }
  }

  /* set problem flags for writing exodus data base and solving problem */

  /*
   * Some preliminaries to help setup EXODUS II database output.
   */

  tnv = cnt_nodal_vars();
  /*  tnv_post is calculated in load_nodal_tkn*/
  tev = cnt_elem_vars(exo);
  /*  tev_post is calculated in load_elem_tkn*/

  if (tnv < 0) {
    DPRINTF(stderr, "%s:\tbad tnv.\n", yo);
    GOMA_EH(GOMA_ERROR, "\t");
  }

  if (tev < 0) {
    DPRINTF(stderr, "%s:\tMaybe bad tev? See goma design committee ;) \n", yo);
    GOMA_EH(GOMA_ERROR, "\t");
  }

  /*
   * When using overlap AC's for fluid/solid overlapping grid problems,
   * create the additional set of AC constraints here.
   */
  if (Do_Overlap && augc[nAC - 1].Type == AC_OVERLAP) {
    err = create_overlap_acs(exo, nAC - 1);
    GOMA_EH(err, "Problem with create_overlap_acs!");
  }
  if (nAC > 0 && augc[0].Type == AC_PERIODIC) {
    err = create_periodic_acs(exo);
    GOMA_EH(err, "Problem with create_periodic_acs!");
  }

  /*
   *  Malloc the space for the results description structure and set all
   *  of that space to naught initially.
   *  Do this only once if in library mode.
   */
  if (callnum == 1) {
    rd = alloc_struct_1(struct Results_Description, 1);
#ifdef LIBRARY_MODE
    xp_id = alloc_int_1(MAX_EXTERNAL_FIELD, 0);
#endif
  }

  rd->nev = 0;       /* number element variables in results */
  rd->nhv = 0;       /* number history variables in results */
  rd->ngv = 6 + nAC; /* number global variables in results
                        see load_global_var_info for names*/

  error = load_global_var_info(rd, 0, "CONV");
  error = load_global_var_info(rd, 1, "NEWT_IT");
  error = load_global_var_info(rd, 2, "MAX_IT");
  error = load_global_var_info(rd, 3, "CONVORDER");
  error = load_global_var_info(rd, 4, "CONVRATE");
  error = load_global_var_info(rd, 5, "MESH_VOLUME");

  if (rd->ngv > MAX_NGV)
    GOMA_EH(GOMA_ERROR, "Augmenting condition values overflowing MAX_NGV.  Change and rerun .");

  if (callnum == 1) {
    Spec_source_inventory = Dmatrix_birth(upd->Num_Mat, upd->Max_Num_Species_Eqn + 1);
    Spec_source_lumped_mass = alloc_dbl_1(exo->num_nodes, 0.0);
  }

  if (nAC > 0) {
    char name[20];

    for (i = 0; i < nAC; i++) {
      sprintf(name, "AUGC_%d", i + 1);
      error = load_global_var_info(rd, 6 + i, name);
    }
  }

  gv = alloc_dbl_1(rd->ngv, 0.0);

  /*
   *  Load output nodal types, kinds and names into the structure
   *  which will be used to define what's in the output file.
   */
  error = load_nodal_tkn(rd, &tnv, &tnv_post);
  if (error != 0) {
    DPRINTF(stderr, "%s:  problem with load_nodal_tkn()\n", yo);
    GOMA_EH(GOMA_ERROR, "\t");
  }

  /*
   * Post-processing vars are stored in the static local array xp_id,
   * and must be used to reset Export_XP_ID on calls after the first.
   * This is a temporary fix!
   */
#ifdef LIBRARY_MODE
  if (callnum == 1) {
    for (i = 0; i < MAX_EXTERNAL_FIELD; i++)
      xp_id[i] = Export_XP_ID[i];
  } else {
    for (i = 0; i < MAX_EXTERNAL_FIELD; i++)
      Export_XP_ID[i] = xp_id[i];
  }
#endif

  /* For retrieving post-processing variables */
  if (Num_Export_XP > 0 && POROUS_SATURATION != -1) {
    asdv(&x_pp, Num_Export_XS * exo->num_nodes);
  }

  /*
   *  Load output element var types, kinds and names into the structure
   *  which will be used to define what's in the output file.
   */
  error = load_elem_tkn(rd, exo, tev, &tev_post);
  if (error != 0) {
    DPRINTF(stderr, "%s:  problem with load_elem_tkn()\n", yo);
    GOMA_EH(GOMA_ERROR, "\t");
  }
#ifdef PARALLEL
  check_parallel_error("Results file error");
#endif /* PARALLEL */

  /*
   * Write out the names of the nodal variables that we will be sending to
   * the EXODUS II output file later - do only once if in library mode.
   */

  if (callnum == 1) {
    gvec_elem = (double ***)alloc_ptr_1(exo->num_elem_blocks);
    if ((tev + tev_post) > 0) {
      for (i = 0; i < exo->num_elem_blocks; i++) {
        gvec_elem[i] = (double **)alloc_ptr_1(tev + tev_post);
      }
    }
    wr_result_prelim_exo(rd, exo, ExoFileOut, gvec_elem);
  }

  /*
   * This gvec workhorse transports output variables as nodal based vectors
   * that are gather from the solution vector. Note: it is NOT a global
   * vector at all and only carries this processor's nodal variables to
   * the exodus database.
   */

  asdv(&gvec, Num_Node);

  /*
   * Allocate space and manipulate for all the nodes that this processor
   * is aware of...
   */

  num_total_nodes = dpi->num_universe_nodes;

  numProcUnknowns = NumUnknowns[pg->imtrx] + NumExtUnknowns[pg->imtrx];

  asdv(&resid_vector, numProcUnknowns);
  asdv(&resid_vector_sens, numProcUnknowns);
  asdv(&scale, numProcUnknowns);

  /*
   * Allocate Aztec structures and initialize all elements to zero
   */
  if (callnum == 1) {
    for (i = 0; i < NUM_ALSS; i++) {
      ams[i] = alloc_struct_1(struct GomaLinearSolverData, 1);
    }
  }
#ifdef MPI
  AZ_set_proc_config(ams[0]->proc_config, MPI_COMM_WORLD);
#else  /* MPI */
  AZ_set_proc_config(ams[0]->proc_config, 0);
#endif /* MPI */

  /* Allocate solution arrays on first call only */
  if (callnum == 1) {
    x_save = alloc_dbl_1(numProcUnknowns, 0.0);
    x = alloc_dbl_1(numProcUnknowns, 0.0);
    x_old = alloc_dbl_1(numProcUnknowns, 0.0);
    x_older = alloc_dbl_1(numProcUnknowns, 0.0);
    x_oldest = alloc_dbl_1(numProcUnknowns, 0.0);
    xdot_save = alloc_dbl_1(numProcUnknowns, 0.0);
    xdot = alloc_dbl_1(numProcUnknowns, 0.0);
    xdot_old = alloc_dbl_1(numProcUnknowns, 0.0);
    xdot_older = alloc_dbl_1(numProcUnknowns, 0.0);
  }
  x_update = alloc_dbl_1(numProcUnknowns + numProcUnknowns, 0.0);

  /* Initialize solid inertia flag */
  set_solid_inertia();

  if (tran->solid_inertia) {
    tran->xdbl_dot = (double *)array_alloc(1, numProcUnknowns, sizeof(double));
    tran->xdbl_dot_old = (double *)array_alloc(1, numProcUnknowns, sizeof(double));

    for (i = 0; i < numProcUnknowns; i++) {
      tran->xdbl_dot[i] = 0.;
      tran->xdbl_dot_old[i] = 0.;
    }

    /*set these in input deck when ready, or if you feel compelled
      to change them.  PRS 10/3/2001 */
    tran->newmark_gamma = 0.9;
    tran->newmark_beta = 0.49;
  }

  node_to_fill = alloc_int_1(num_total_nodes, 0);

  pg->matrices = malloc(sizeof(struct Matrix_Data));
  pg->matrices[pg->imtrx].ams = ams[JAC];
  pg->matrices[pg->imtrx].x = x;
  pg->matrices[pg->imtrx].x_old = x_old;
  pg->matrices[pg->imtrx].x_older = x_older;
  pg->matrices[pg->imtrx].xdot = xdot;
  pg->matrices[pg->imtrx].xdot_old = xdot_old;
  pg->matrices[pg->imtrx].x_update = x_update;
  pg->matrices[pg->imtrx].scale = scale;
  pg->matrices[pg->imtrx].resid_vector = resid_vector;

  /* Allocate sparse matrix */

  ams[JAC]->GomaMatrixData = NULL;
  if ((strcmp(Matrix_Format, "tpetra") == 0) || (strcmp(Matrix_Format, "epetra") == 0)) {
    err = check_compatible_solver();
    GOMA_EH(err, "Incompatible matrix solver for tpetra, tpetra supports stratimikos");
    check_parallel_error("Matrix format / Solver incompatibility");
    GomaSparseMatrix goma_matrix;
    goma_error err = GomaSparseMatrix_CreateFromFormat(&goma_matrix, Matrix_Format);
    GOMA_EH(err, "GomaSparseMatrix_CreateFromFormat");
    int local_nodes = Num_Internal_Nodes + Num_Border_Nodes + Num_External_Nodes;
    err = GomaSparseMatrix_SetProblemGraph(
        goma_matrix, num_internal_dofs[pg->imtrx], num_boundary_dofs[pg->imtrx],
        num_external_dofs[pg->imtrx], local_nodes, Nodes, MaxVarPerNode, Matilda, Inter_Mask, exo,
        dpi, cx[pg->imtrx], pg->imtrx, Debug_Flag, ams[JAC]);
    GOMA_EH(err, "GomaSparseMatrix_SetProblemGraph");
    ams[JAC]->GomaMatrixData = goma_matrix;
#ifdef GOMA_ENABLE_PETSC
#if PETSC_USE_COMPLEX
  } else if (strcmp(Matrix_Format, "petsc_complex") == 0) {
    err = check_compatible_solver();
    GOMA_EH(err, "Incompatible matrix solver for petsc, solver must be petsc");
    check_parallel_error("Matrix format / Solver incompatibility");
    pg->imtrx = 0;
    goma_error err = goma_setup_petsc_matrix_complex(
        ams[JAC], exo, dpi, x, x_old, xdot, xdot_old, num_internal_dofs[pg->imtrx],
        num_boundary_dofs[pg->imtrx], num_external_dofs[pg->imtrx], pg->imtrx);
    GOMA_EH(err, "goma_setup_petsc_matrix_complex");
#else
  } else if (strcmp(Matrix_Format, "petsc") == 0) {
    err = check_compatible_solver();
    GOMA_EH(err, "Incompatible matrix solver for petsc, solver must be petsc");
    check_parallel_error("Matrix format / Solver incompatibility");
    pg->imtrx = 0;
    goma_error err = goma_setup_petsc_matrix(
        ams[JAC], exo, dpi, x, x_old, xdot, xdot_old, num_internal_dofs[pg->imtrx],
        num_boundary_dofs[pg->imtrx], num_external_dofs[pg->imtrx], pg->imtrx);
    GOMA_EH(err, "goma_setup_petsc_matrix");
#endif
#endif
  } else if (strcmp(Matrix_Format, "msr") == 0) {
    log_msg("alloc_MSR_sparse_arrays...");
    alloc_MSR_sparse_arrays(&ija, &a, &a_old, 0, node_to_fill, exo, dpi);
    /*
     * An attic to store external dofs column names is needed when
     * running in parallel.
     */
    alloc_extern_ija_buffer(num_universe_dofs[pg->imtrx],
                            num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx], ija,
                            &ija_attic);
    /*
     * Any necessary one time initialization of the linear
     * solver package (Aztec).
     */
    ams[JAC]->bindx = ija;
    ams[JAC]->val = a;
    ams[JAC]->belfry = ija_attic;
    ams[JAC]->val_old = a_old;

    /*
     * These point to nowhere since we're using MSR instead of VBR
     * format.
     */

    ams[JAC]->indx = NULL;
    ams[JAC]->bpntr = NULL;
    ams[JAC]->rpntr = NULL;
    ams[JAC]->cpntr = NULL;

    ams[JAC]->npn = dpi->num_internal_nodes + dpi->num_boundary_nodes;
    ams[JAC]->npn_plus =
        dpi->num_internal_nodes + dpi->num_boundary_nodes + dpi->num_external_nodes;

    ams[JAC]->npu = num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx];
    ams[JAC]->npu_plus = num_universe_dofs[pg->imtrx];

    ams[JAC]->nnz = ija[num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx]] - 1;
    ams[JAC]->nnz_plus = ija[num_universe_dofs[pg->imtrx]];

  } else if (strcmp(Matrix_Format, "vbr") == 0) {
    log_msg("alloc_VBR_sparse_arrays...");
    alloc_VBR_sparse_arrays(ams[JAC], exo, dpi);
    ija_attic = NULL;
    ams[JAC]->belfry = ija_attic;

    a = ams[JAC]->val;
    if (!save_old_A)
      a_old = ams[JAC]->val_old = NULL;
  } else {
    GOMA_EH(GOMA_ERROR, "Attempted to allocate unknown sparse matrix format: %s", Matrix_Format);
  }
#ifdef GOMA_ENABLE_PETSC
#if !(PETSC_USE_COMPLEX)
  if (upd->petsc_solve_post_proc && rd->TotalNVPostOutput) {
    goma_setup_petsc_post_proc_matrix(exo, dpi, x, x_old, xdot, xdot_old);
  }
#endif
#endif

  /*
   * allocate memory for Volume Constraint Jacobian. ACS 2/99
   * and level set velocity constraint Jacobian. AMG 4/02
   */
  if (nAC > 0) {
    for (iAC = 0; iAC < nAC; iAC++) {
      augc[iAC].d_evol_dx = alloc_dbl_1(numProcUnknowns, 0.0);
      augc[iAC].d_lsvel_dx = alloc_dbl_1(numProcUnknowns, 0.0);
      augc[iAC].d_lsvol_dx = alloc_dbl_1(numProcUnknowns, 0.0);
    }
  }

  /*
   *  if compute parameter sensitivities, allocate space for solution
   *  sensitivity vectors
   */

  for (i = 0; i < nn_post_fluxes_sens; i++) {
    num_pvector = MAX(num_pvector, pp_fluxes_sens[i]->vector_id);
  }
  for (i = 0; i < nn_post_data_sens; i++) {
    num_pvector = MAX(num_pvector, pp_data_sens[i]->vector_id);
  }

  if ((nn_post_fluxes_sens + nn_post_data_sens) > 0) {
    num_pvector++;
    num_pvector = MAX(num_pvector, 2);
    x_sens_p = Dmatrix_birth(num_pvector, numProcUnknowns);
    x_sens = alloc_dbl_1(numProcUnknowns, 0.0);
  } else {
    x_sens_p = NULL;
  }

  /* Allocate AC unknown arrays on the first call */
  if (nAC > 0 && callnum == 1) {
    x_AC = alloc_dbl_1(nAC, 0.0);
    x_AC_old = alloc_dbl_1(nAC, 0.0);
    x_AC_older = alloc_dbl_1(nAC, 0.0);
    x_AC_oldest = alloc_dbl_1(nAC, 0.0);
    x_AC_dot = alloc_dbl_1(nAC, 0.0);
    x_AC_dot_old = alloc_dbl_1(nAC, 0.0);
    x_AC_dot_older = alloc_dbl_1(nAC, 0.0);
    x_AC_pred = alloc_dbl_1(nAC, 0.0);
  }

  /* Set initial guess from an input exodus file or other method on the first call only */
  if (callnum == 1) {
    init_vec(x, cx[0], exo, dpi, x_AC, nAC, &timeValueRead);

    /*
     *  Determine if we should use this time as the initial time in the simulation
     */
#ifndef ALLOW_NEGATIVE_TIMES_PLEASE
    if (TimeIntegration != STEADY) {
      if (tran->init_time < 0.0) {
        tran->init_time = timeValueRead;
        DPRINTF(stdout, "\n Initial Simulation Time Has been set to %g\n", timeValueRead);
      }
    }
#endif
  }

  if (Conformation_Flag == 1) // If mapping is needed for log-conformation tensor
  {
    initial_guess_stress_to_log_conf(x, num_total_nodes);
  }

  if ((callnum == 1) && (pmv_hyst != NULL)) {
    error = init_pmv_hyst(exo);
    GOMA_EH(error, "Error in initiating Porous Media Variables Hysteresis Struct");
  }
  /* Load external fields from import vectors xnv_in & xev_in */
#ifdef LIBRARY_MODE
  /* First check if porosity updates are necessary */
  if (update_porosity) {
    base_p_por = alloc_dbl_1(num_total_nodes, 0.0);
    base_p_liq = alloc_dbl_1(num_total_nodes, 0.0);

    /* Load starting porous liquid pressures from solution vector */
    for (i = 0; i < num_total_nodes; i++) {
      j = Index_Solution(i, POR_LIQ_PRES, 0, 0, -1, pg->imtrx);
      if (j > -1)
        base_p_liq[i] = x[j];
    }
  }

  /* Now load imports, and also base_p_por if applicable */
  error = load_import_fields(base_p_por, exo, callnum);
  GOMA_EH(error, "Problem with load_import_fields!");
#else

  /****************************Anneal from external***********************/
  if (efv->ev_porous_decouple) {
    anneal_mesh_with_external_field(exo);
  }

#endif

  dcopy1(nAC, x_AC, &(gv[5]));

  if (Output_Variable_Stats) {
    err = variable_stats(x, timeValueRead, FALSE);
    GOMA_EH(err, "Problem with variable_stats!");
    if (ProcID == 0)
      fflush(stdout);
  }

  /***************************************************************************
   *            STEADY STATE SOLUTION PROCEDURE
   ***************************************************************************/
  if (TimeIntegration == STEADY) {

    theta = 0.0; /* for steady problems. theta def in rf_fem.h */
    delta_t = 0.0;

    find_and_set_Dirichlet(x, xdot, exo, dpi);

    matrix_systems_mask = 1;

    log_msg("sl_init()...");
    sl_init(matrix_systems_mask, ams, exo, dpi, cx[0]);
    if (nAC > 0 || nn_post_fluxes_sens > 0 || nn_post_data_sens > 0)
      ams[JAC]->options[AZ_keep_info] = 1;

      /*
       * Now, just pass pointer to ams structure with all Aztec stuff
       * bundled inside. At the other end, extract "ija" and "a" as
       * appropriate, but the other items are there now, too.
       */

#ifdef PARALLEL

    /*
     * Make sure the solver was properly initialized on all processors.
     */
    check_parallel_error("Solver initialization problems");
#endif /* PARALLEL */

    if (nEQM > 0) {
      DPRINTF(stdout, "\nINITIAL ELEMENT QUALITY CHECK---\n");
      good_mesh = element_quality(exo, x, ams[0]->proc_config);
    }

    err = solve_nonlinear_problem(ams[JAC], x, delta_t, theta, x_old, x_older, xdot, xdot_old,
                                  resid_vector, x_update, scale, &converged, &nprint, tev, tev_post,
                                  gv, rd, gindex, p_gsize, gvec, gvec_elem, time1, exo, dpi, cx[0],
                                  0, &time_step_reform, is_steady_state, x_AC, x_AC_dot, time1,
                                  resid_vector_sens, x_sens, x_sens_p, NULL);

    if (!converged) {
      DPRINTF(stderr, "\n\tFailure: could not converge on the steady state solution.\n");
      DPRINTF(stderr, "\tSorry, I really tried...\n");
      if (Linear_Stability) {
        DPRINTF(stderr, "\tThere's no point in solving the eigensystem (but I'll do it anyway).\n");
      }
    }

    log_msg("Returning from solve_nonlinear_problem with %d", err);
    GOMA_EH(err, "Problem from solve_nonlinear_problem.");

    /* Check element quality */
    good_mesh = element_quality(exo, x, ams[0]->proc_config);

    if (file != NULL) {
      error = write_ascii_soln(x, resid_vector, numProcUnknowns, x_AC, nAC, 0.0, file);
    }

    GOMA_EH(error, "Error writing ASCII soln file.");

    if (Write_Intermediate_Solutions == 0) {
      nprint = 0;

      write_solution(ExoFileOut, resid_vector, x, x_sens_p, x_old, xdot, xdot_old, tev, tev_post,
                     gv, rd, gvec, gvec_elem, &nprint, delta_t, theta, time1, x_pp, exo, dpi);
    } /* end of if Write Intermediate Solutions */

    /* Print out values of extra unknowns from augmenting conditions */
    if (nAC > 0) {
      DPRINTF(stdout, "\n------------------------------\n");
      DPRINTF(stdout, "Augmenting Conditions:    %4d\n", nAC);
      DPRINTF(stdout, "Number of extra unknowns: %4d\n\n", nAC);

      for (iAC = 0; iAC < nAC; iAC++) {
        evol_local = augc[iAC].evol;
#ifdef PARALLEL
        if (Num_Proc > 1 && (augc[iAC].Type == AC_VOLUME || augc[iAC].Type == AC_POSITION ||
                             augc[iAC].Type == AC_ANGLE || augc[iAC].Type == AC_POSITION_MT)) {
          MPI_Allreduce(&evol_local, &evol_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
          evol_local = evol_global;
        }
#endif
        if (augc[iAC].Type == AC_USERBC) {
          DPRINTF(stdout, "\tBC[%4d] DF[%4d]=% 10.6e\n", augc[iAC].BCID, augc[iAC].DFID, x_AC[iAC]);
        } else if (augc[iAC].Type == AC_USERMAT || augc[iAC].Type == AC_FLUX_MAT) {
          DPRINTF(stdout, "\tMT[%4d] MP[%4d]=% 10.6e\n", augc[iAC].MTID, augc[iAC].MPID, x_AC[iAC]);
        } else if (augc[iAC].Type == AC_VOLUME) {
          DPRINTF(stdout, "\tMT[%4d] VC[%4d]=%10.6e Param=%10.6e\n", augc[iAC].MTID,
                  augc[iAC].VOLID, evol_local, x_AC[iAC]);
        } else if (augc[iAC].Type == AC_POSITION || augc[iAC].Type == AC_POSITION_MT) {
          DPRINTF(stdout, "\tNodeSet[%4d]_Pos = %10.6e F_bal = %10.6e MT[%4d] Param=%10.6e\n",
                  augc[iAC].MTID, evol_local, augc[iAC].lm_resid, augc[iAC].VOLID, x_AC[iAC]);
        } else if (augc[iAC].Type == AC_ANGLE) {
          evol_local = augc[iAC].lm_resid + augc[iAC].CONSTV;
          DPRINTF(stdout, "\tNodeSet[%4d]_Ang = %g F_bal = %6.3e MT[%4d] Param=%6.3e\n",
                  augc[iAC].MTID, evol_local, augc[iAC].lm_resid, augc[iAC].VOLID, x_AC[iAC]);
        } else if (augc[iAC].Type == AC_LS_VEL) {
          evol_local = augc[iAC].lsvol;
          lsvel_local = augc[iAC].lsvel;

#ifdef PARALLEL
          if (Num_Proc > 1) {
            MPI_Allreduce(&evol_local, &evol_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&lsvel_local, &lsvel_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            evol_local = evol_global;
            lsvel_local = lsvel_global;
          }
#endif /* PARALLEL */

          lsvel_local = lsvel_local / evol_local;

          DPRINTF(stdout, "\tMT[%4d] LSVEL phase[%4d]=%10.6e Param=%10.6e\n", augc[iAC].MTID,
                  augc[iAC].LSPHASE, lsvel_local, x_AC[iAC]);
        } else if (augc[iAC].Type == AC_FLUX) {
          DPRINTF(stdout, "\tBC[%4d] DF[%4d]=%10.6e\n", augc[iAC].BCID, augc[iAC].DFID, x_AC[iAC]);
        }
      }
    }

    /* Integrate fluxes, forces  (only if converged)
     */
    if (converged) {
      for (i = 0; i < nn_post_fluxes; i++) {
        (void)evaluate_flux(exo, dpi, pp_fluxes[i]->ss_id, pp_fluxes[i]->flux_type,
                            pp_fluxes[i]->flux_type_name, pp_fluxes[i]->blk_id,
                            pp_fluxes[i]->species_number, pp_fluxes[i]->flux_filenm,
                            pp_fluxes[i]->profile_flag, x, xdot, NULL, delta_t, time1, 1);
      }

      /* Compute flux, force sensitivities
       */
      for (i = 0; i < nn_post_fluxes_sens; i++) {
        (void)evaluate_flux_sens(exo, dpi, pp_fluxes_sens[i]->ss_id, pp_fluxes_sens[i]->flux_type,
                                 pp_fluxes_sens[i]->flux_type_name, pp_fluxes_sens[i]->blk_id,
                                 pp_fluxes_sens[i]->species_number, pp_fluxes_sens[i]->sens_type,
                                 pp_fluxes_sens[i]->sens_id, pp_fluxes_sens[i]->sens_flt,
                                 pp_fluxes_sens[i]->sens_flt2, pp_fluxes_sens[i]->vector_id,
                                 pp_fluxes_sens[i]->flux_filenm, pp_fluxes_sens[i]->profile_flag, x,
                                 xdot, x_sens_p, delta_t, time1, 1);
      }

      /*
       * Compute global volumetric quantities
       */
      for (i = 0; i < nn_volume; i++) {
        evaluate_volume_integral(exo, dpi, pp_volume[i]->volume_type, pp_volume[i]->volume_name,
                                 pp_volume[i]->blk_id, pp_volume[i]->species_no,
                                 pp_volume[i]->volume_fname, pp_volume[i]->params,
                                 pp_volume[i]->num_params, NULL, x, xdot, delta_t, time1, 1);
      }
      if (Output_Variable_Stats) {
        err = variable_stats(x, time1, Output_Variable_Regression);
        GOMA_EH(err, "Problem with variable_stats!");
        if (ProcID == 0)
          fflush(stdout);
      }

      DPRINTF(stderr, "\t\tconverged SS!\n");
    } /* if converged */

    if (Anneal_Mesh) {
      /*
       * Transform the node point coordinates according to the
       * displacements and write out all the results using the
       * displaced coordinates. Set the displacement field to
       * zero, too.
       */

      err = anneal_mesh(x, tev, tev_post, gv, rd, time1, exo, dpi);
      GOMA_EH(err, "anneal_mesh() bad return.");
    }

    if (Linear_Stability) {
      err =
          solve_stability_problem(ams[JAC], x, delta_t, theta, resid_vector, x_old, x_older, xdot,
                                  xdot_old, x_update, &converged, &nprint, tnv, tnv_post, tev,
                                  tev_post, rd, gindex, p_gsize, gvec, gvec_elem, time1, exo, dpi);
      GOMA_EH(err, "Problem from solve_stability_problem.");
      DPRINTF(stderr, "\t\tcompleted LSA!\n");
    }

    if (Particle_Dynamics) {
      /* If this was steady-state we need to ensure the fv_old* values are the same as fv. */
      dcopy1(NumUnknowns[pg->imtrx], x, x_old);
      dcopy1(NumUnknowns[pg->imtrx], x, x_older);
      initialize_particles(exo, x, x_old, xdot, xdot_old, resid_vector);
      for (n = 0; n < Particle_Max_Time_Steps; n++) {
        time = (n + 1) * Particle_Output_Time_Step;
        DPRINTF(stdout, "\nComputing particles for time %g (%2.0f%% done)\n", time,
                (dbl)n / (dbl)Particle_Max_Time_Steps * 100.0);
        err = compute_particles(exo, x, x_old, xdot, xdot_old, resid_vector, time,
                                Particle_Output_Time_Step, n);
        GOMA_EH(err, "Error performing particle calculations.");
      }
    }
  } /* if(steady) */

  /********************************************************************************
   *                            Transient solution process
   ********************************************************************************/
  else {
    if (Debug_Flag && ProcID == 0) {
      fprintf(stdout, "MaxTimeSteps: %d \tTimeMax: %f\n", tran->MaxTimeSteps, tran->TimeMax);
      fprintf(stdout, "solving transient problem\n");
    }

    /*
     *  Transfer information from the Transient_Information structure to local variables
     */
    double delta_t0 = tran->Delta_t0;
    double delta_t_min = tran->Delta_t_min;
    double delta_t_max = tran->Delta_t_max;
    double max_time_steps = tran->MaxTimeSteps;
    double time_max = tran->TimeMax;
    eps = tran->eps;

    // Determine if we are using a constant time step or not
    if (delta_t0 < 0.0) {
      delta_t0 = -delta_t0;
      const_delta_t = 1;
    } else {
      const_delta_t = 0;
    }

    /*
     * If this is a second Goma call in LIBRARY_MODE, then set up step
     * size and maximum steps as follows:
     * solve_steps = 0: solve until end time passed in, use Goma max. steps
     * solve_steps < 0: Take |solve| steps, reset delta_t to initial value.
     * solve_steps > 0: Take |solve| steps, restore delta_t from last call.
     */
#ifdef LIBRARY_MODE
    if (libio->solve_steps > 0) {
      MaxTimeSteps = libio->solve_steps;
      if (callnum > 1)
        Delta_t0 = delta_t_save;
    } else if (libio->solve_steps < 0) {
      MaxTimeSteps = -libio->solve_steps;
      if (callnum > 1)
        Delta_t0 = delta_t_save;
    } else {
      if (callnum > 1 && delta_t_save * libio->decelerator > Delta_t0) {
        Delta_t0 = delta_t_save * libio->decelerator;
      }
    }
#endif

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
    x_pred = alloc_dbl_1(numProcUnknowns, 0.0);
    x_pred_static = x_pred;

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
    find_and_set_Dirichlet(x, xdot, exo, dpi);

    /*
     * Before propagating x_n back into the historical records in
     * x_n-1, x_n-2 and x_n-3, make sure external dofs are the best
     * they can be. That is, ask the processors that are supposed to
     * know...
     */
    exchange_dof(cx[0], dpi, x, 0);

    /*
     * Now copy the initial solution, x[], into the history solutions
     * x_old[], etc. Note, xdot[] = xdot_old[] = 0 at this point,
     * which is in agreement with the specification of the history
     * solutions.
     */
    dcopy1(numProcUnknowns, x, x_old);
    dcopy1(numProcUnknowns, x_old, x_older);
    dcopy1(numProcUnknowns, x_older, x_oldest);

    if (nAC > 0) {
      dcopy1(nAC, x_AC, x_AC_old);
      dcopy1(nAC, x_AC_old, x_AC_older);
      dcopy1(nAC, x_AC_older, x_AC_oldest);
    }

    /* initialize the counters for when to print out data */
    time_print = time;
    step_print = 1;
    matrix_systems_mask = 1;

    /*set some other misc. action flags */
    af->Sat_hyst_reevaluate = FALSE;

    /* Call prefront (or mf_setup) if necessary */

    /*
     * Now, just pass pointer to ams structure with all Aztec stuff
     * bundled inside. At the other end, extract "ija" and "a" as
     * appropriate, but the other items are there now, too.
     * Do this only once if in library mode.
     */
    if (callnum == 1)
      sl_init(matrix_systems_mask, ams, exo, dpi, cx[0]);

    /*
     * make sure the Aztec was properly initialized
     */
    check_parallel_error("Aztec Initialization");

    /* Set the number of successful time steps, nt, to zero */
    nt = 0;
    time_step_reform = Time_Jacobian_Reformation_stride;
    const_delta_ts = const_delta_t;
    last_renorm_nt = 0;

    if (Particle_Dynamics)
      initialize_particles(exo, x, x_old, xdot, xdot_old, resid_vector);

    /*
     * Write out the initial solution to an ascii file
     * and to the exodus output file, if requested to do so by
     * an optional flag in the input file
     *  -> Helpful in debugging what's going on.
     */
    if (Write_Initial_Solution) {
      if (file != NULL) {
        error = write_ascii_soln(x, resid_vector, numProcUnknowns, x_AC, nAC, time, file);
        if (error != 0)
          DPRINTF(stderr, "%s:  error writing ASCII soln file\n", yo);
      }
      (void)write_solution(ExoFileOut, resid_vector, x, x_sens_p, x_old, xdot, xdot_old, tev,
                           tev_post, gv, rd, gvec, gvec_elem, &nprint, delta_t, theta, time, x_pp,
                           exo, dpi);
      nprint++;
    }

    /*
     * In order to write an updated final solution, one extra call
     * is made to solve_problem just to call write_solution.
     * Now, reset MaxTimeSteps to skip the time step loop.
     */
#ifdef LIBRARY_MODE
    if (last_call && libio->goma_first == 1 && libio->print_flag >= 0)
      MaxTimeSteps = 0;
#endif

    /* Initial element quality check (if requested) */
    if (nEQM > 0) {
      DPRINTF(stdout, "\nINITIAL ELEMENT QUALITY CHECK---\n");
      good_mesh = element_quality(exo, x, ams[0]->proc_config);
    }

    /*******************************************************************
     *  TOP OF THE TIME STEP LOOP -> Loop over time steps whether
     *                               they be successful or not
     *******************************************************************/
    for (n = 0; n < max_time_steps; n++) {
      /*
       * Calculate the absolute time for the current step, time1
       */
      time1 = time + delta_t;
#ifdef LIBRARY_MODE
      delta_t_save = delta_t;
#endif
      if (time1 > time_max) {
        DPRINTF(stdout, "\t\tLAST TIME STEP!\n");
        time1 = time_max;
        delta_t = time1 - time;
        tran->delta_t = delta_t;
        tran->delta_t_avg = 0.25 * (delta_t + delta_t_old + delta_t_older + delta_t_oldest);
#ifdef LIBRARY_MODE
        last_step = TRUE;
#endif
      }
      tran->time_value = time1;
#ifdef LIBRARY_MODE
      if (n == (MaxTimeSteps - 1))
        last_step = TRUE;
#endif
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
        for (w = 0; w < efv->Num_external_field; w++) {
          if (strcmp(efv->field_type[w], "transient") == 0) {
            err = rd_trans_vectors_from_exoII(x_old, efv->file_nm[w], w, n, &timeValueReadTrans,
                                              exo, cx[0], dpi);
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
      if ((nt - last_renorm_nt) == 0 || (nt - last_adapt_nt) == 0) {
        theta = 0.0;
        const_delta_t = 1.0;

      } else if ((nt - last_renorm_nt) >= 3 && (nt - last_adapt_nt) >= 2) {
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

      find_and_set_Dirichlet(x, xdot, exo, dpi);

      if (nt == 0) {
        xfem = NULL;
        if (upd->XFEM) {
          xfem = alloc_struct_1(struct Extended_Shape_Fcn_Basics, 1);
          xfem->ielem = -1;
          xfem->tot_vol = alloc_dbl_1(numProcUnknowns, 0.0);
          xfem->active_vol = alloc_dbl_1(numProcUnknowns, 0.0);
          if (ls == NULL) {
            GOMA_EH(GOMA_ERROR, "Currently, XFEM requires traditional level set (not pf)");
          }
        }
      }

      /*
       * Initial Start up (t=0) for the FILL/LEVEL_SET equations
       */

      if (upd->ep[pg->imtrx][FILL] > -1 && nt == 0) { /*  Start of LS initialization */

        if (ls != NULL || pfd != NULL) {

          int eqntype = ls->Init_Method;
          /* This is a temporary loc for this allocation */

          switch (ls->Evolution) {
          case LS_EVOLVE_ADVECT_EXPLICIT:
            DPRINTF(stdout, "\n\t Using decoupled / subcycling for FILL equation.\n");
            break;
          case LS_EVOLVE_ADVECT_COUPLED:
            DPRINTF(stdout, "\n\t Using Coupled Level Set evolution!\n");
            break;
          case LS_EVOLVE_SLAVE:
            DPRINTF(stdout, "\n\t USING SLAVE LEVEL SET INTERFACE\n");
            break;
          case LS_EVOLVE_SEMILAGRANGIAN:
            DPRINTF(stdout, "\n\t Using semi-Lagrangian Level Set Evolution\n");
            break;
          default:
            GOMA_EH(GOMA_ERROR, "Level Set Evolution scheme not found \n");
          }

          if (ls->Length_Scale < 0.0)
            GOMA_EH(GOMA_ERROR, "\tError: a Level Set Length Scale needs to be specified\n");

          if (ls->Integration_Depth > 0 || ls->SubElemIntegration || ls->AdaptIntegration) {

            if (ls->Integration_Depth > 0) {
              int first_elem;

              first_elem = find_first_elem_with_var(exo, LS);

              if (first_elem != -1) {
                load_ei(first_elem, exo, 0, pg->imtrx);

                Subgrid_Tree = create_shape_fcn_tree(ls->Integration_Depth);
                DPRINTF(stdout, "\n\tSubgrid Integration of level set interface active.\n");
              }
            } else if (ls->SubElemIntegration) {
              DPRINTF(stdout, "\n\tSubelement Integration of level set interface active.\n");
            } else if (ls->AdaptIntegration) {
              DPRINTF(stdout, "\n\tAdaptive Integration of level set interface active.\n");
              DPRINTF(stdout, "\tAdaptive Integration Interface Order = %d\n", ls->Adaptive_Order);
            }
            Subgrid_Int.ip_total = 0;
            Subgrid_Int.s = NULL;
            Subgrid_Int.wt = NULL;
          }

          switch (eqntype) {
          case PROJECT:

            DPRINTF(stdout, "\n\t Projection level set initialization \n");
            GOMA_EH(GOMA_ERROR, "Use of \"PROJECT\" is obsolete.");

            break;

          case EXO_READ:

            DPRINTF(stdout, "\t\t Level set read from exodus database \n");

            break;

          case SURFACES:

            DPRINTF(stdout, "\n\t\t Surface object level set initialization : ");

            /* parallel synchronization of initialization surfaces */
            if (Num_Proc > 1) {
              if (!ls->init_surf_list)
                ls->init_surf_list = create_surf_list();
              assemble_Global_surf_list(ls->init_surf_list);
            }

            surf_based_initialization(x, NULL, NULL, exo, num_total_nodes, ls->init_surf_list, 0.,
                                      0., 0.);

#ifdef DEBUG_PARALLEL
            {
              FILE *data_out;
              int II, je;

              if (ProcID == 0)
                data_out = fopen("Proc0.dat", "w");
              if (ProcID == 1)
                data_out = fopen("Proc1.dat", "w");

              fprintf(data_out, "num_total_nodes:  %d \n", num_total_nodes);

              for (II = 0; II < num_total_nodes; II++) {
                je = Index_Solution(II, LS, 0, 0, -1, pg->imtrx);

                if (je != -1) {
                  fprintf(data_out, "%d  %lf \n", II, x[je]);
                }
              }

              fflush(data_out);
              fclose(data_out);
            }
#endif

            DPRINTF(stdout, "- done \n");

            break;

          case SM_OBJECT:
            GOMA_EH(GOMA_ERROR, "CGM not supported, SM_OBJECT level set initialization");
            break;

          default:
            GOMA_WH(-1, "Level Set Initialization method not found \n");
          } /* end of switch( eqntype )  */

          exchange_dof(cx[0], dpi, x, 0);

          if (converged) {
            switch (ls->Renorm_Method) {

            case HUYGENS:
            case HUYGENS_C:
            case HUYGENS_MASS_ITER:
              Renorm_Now =
                  (ls->Force_Initial_Renorm || (ls->Renorm_Freq != 0 && ls->Renorm_Countdown == 0));

              did_renorm =
                  huygens_renormalization(x, num_total_nodes, exo, cx[0], dpi, num_fill_unknowns,
                                          numProcUnknowns, time1, Renorm_Now);

#ifndef PHASE_COUPLED_FILL
              if (did_renorm) {
                /*get_fill_vector(num_total_nodes,x,xf,node_to_fill);*/
              }
#endif /* not COUPLED_FILL */

              break;

            case CORRECT:
              GOMA_EH(GOMA_ERROR, "Use of \"CORRECT\" is obsolete.");
              break;
            default:
              if (ls->Evolution == LS_EVOLVE_ADVECT_EXPLICIT ||
                  ls->Evolution == LS_EVOLVE_ADVECT_COUPLED)
                GOMA_WH(-1, "No level set renormalization is on.\n");
            } /* end of switch(ls->Renorm_Method ) */
          }
          /*
           * More initialization needed. Have to set those field variables that initially are
           * indexed by level set function.  For example, species  concentration and temperature.
           */
          if (ls->Num_Var_Init > 0)
            ls_var_initialization(&x, exo, dpi, cx);

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

        if (converged) // avoid death spiral on initial failure
        {
          dcopy1(numProcUnknowns, x, x_old);
          dcopy1(numProcUnknowns, x, x_older);
          dcopy1(numProcUnknowns, x, x_oldest);
        }

        exchange_dof(cx[0], dpi, x, 0);
        exchange_dof(cx[0], dpi, x_old, 0);
        exchange_dof(cx[0], dpi, x_oldest, 0);
      }

      ls_old = ls;
      if (upd->vp[pg->imtrx][PHASE1] > -1 && nt == 0) { /* Start of Phase Function initialization */

        if (pfd != NULL) {
          struct Level_Set_Data *ls_save = ls;
          if (upd->vp[pg->imtrx][PHASE1] > -1) {
            switch (pfd->ls[0]->Evolution) {
            case LS_EVOLVE_ADVECT_EXPLICIT:
              DPRINTF(stdout,
                      "\n\t Using decoupled / subcycling for FILL equation for R_PHASE0.\n");
              break;
            case LS_EVOLVE_ADVECT_COUPLED:
              DPRINTF(stdout, "\n\t Using Coupled Level Set evolution! for R_PHASE0\n");
              break;
            case LS_EVOLVE_SLAVE:
              DPRINTF(stdout, "\n\t USING SLAVE LEVEL SET INTERFACE for R_PHASE0\n");
              break;
            case LS_EVOLVE_SEMILAGRANGIAN:
              DPRINTF(stdout, "\n\t Using semi-Lagrangian Level Set Evolution for R_PHASE0\n");
              break;
            default:
              GOMA_EH(GOMA_ERROR, "PHASE Function Evolution scheme not found \n");
              break;
            }
          }

          for (i = 0; i < pfd->num_phase_funcs; i++) {
            ls = pfd->ls[i]; /* This is a crucial step.  It allows us to use the Level Set Machinery
                                to initialize the phase functions
                                BUT... Make sure to set it back to what we start with when through.
                                Many places the test ls != NULL is used to thread the code
                             */

            switch (ls->Init_Method) {
            case SURFACES:
              DPRINTF(stdout, "\n\t\t Surface object initialization for phase function: %d", i);

              /* parallel synchronization of initialization surfaces */
              if (Num_Proc > 1) {
                if (!ls->init_surf_list)
                  ls->init_surf_list = create_surf_list();
                assemble_Global_surf_list(ls->init_surf_list);
              }

              surf_based_initialization(x, NULL, NULL, exo, num_total_nodes, ls->init_surf_list, 0.,
                                        0., 0.);
              break;
            case EXO_READ:
              DPRINTF(stdout, "\n\t\t Exodus file read initialization for phase function fields");
              break;
            } /* end of switch(ls->Init_Method ) */
          }   /* end of i<pfd->num_phase_funcs */

          ls = ls_save; /* OK, used the level set routines now be nice and point
                           ls back to where you found it */

          /* Allocate arrays for constraint Jacobian contributions */
          if (pfd->Use_Constraint == TRUE) {
            int i;
            pfd->jac_info->d_pf_lm = alloc_dbl_1(numProcUnknowns, 0.0);
            pfd->jac_info->d_lm_pf = alloc_dbl_1(numProcUnknowns, 0.0);

            if (pfd->ls[0]->Init_Method != EXO_READ) {
              for (i = 0; i < pfd->num_phase_funcs; i++) {
                shift_nodal_values(PHASE1 + i, -pfd->shift[i], x, num_total_nodes);
              }
            }
          }

        } /* end of pfd!=NULL */
        ls = ls_old;

        dcopy1(numProcUnknowns, x, x_old);
        dcopy1(numProcUnknowns, x, x_older);
        dcopy1(numProcUnknowns, x, x_oldest);
        exchange_dof(cx[0], dpi, x, 0);
        exchange_dof(cx[0], dpi, x_old, 0);
        exchange_dof(cx[0], dpi, x_oldest, 0);

      } /* end of phase function initialization */

      if (ls != NULL && ls->last_surf_list != NULL) {
        /* Find the interface surf at the last full time step (x_old = tmp_x)
         * for use during this explicit step */

        create_subsurfs(ls->last_surf_list, x_old, exo);
      }

      if (ProcID == 0) {
        if (theta == 0.0)
          strcpy(tspstring, "(BE)");
        else if (theta == 0.5)
          strcpy(tspstring, "(CN)");
        else if (theta == 1.0)
          strcpy(tspstring, "(FE)");
        else
          sprintf(tspstring, "(TSP %3.1f)", theta);
        fprintf(stdout, "\n=> Try for soln at t=%g with dt=%g [%d for %d] %s\n", time1, delta_t, nt,
                n, tspstring);
        log_msg("Predicting try at t=%g, dt=%g [%d for %d so far] %s", time1, delta_t, nt, n,
                tspstring);
      }

      /*
       * Predict the solution, x[], and its derivative, xdot[],
       * at the new time, time1, using the old solution, xdot_old[],
       * And its derivatives at the old time, time.
       */

      if (!nonconv_roll) {
        predict_solution(numProcUnknowns, delta_t, delta_t_old, delta_t_older, theta, x, x_old,
                         x_older, x_oldest, xdot, xdot_old, xdot_older);

        if (tran->solid_inertia) {
          predict_solution_newmark(num_total_nodes, delta_t, x, x_old, xdot, xdot_old);
          exchange_dof(cx[0], dpi, tran->xdbl_dot, 0);
        }
      } else {
        DPRINTF(stderr, "skipping predict_solution at time: %g %d\n", time1, nonconv_roll);
      }

#ifdef LASER_RAYTRACE
      if (ls != NULL) {
        double(*point0)[DIM] = NULL;
        double(*point1)[DIM] = NULL;
        int *owning_elem = NULL;
        int facet, num_facets;

        num_facets = generate_facet_list(&point0, &point1, &owning_elem, x, exo);
        for (facet = 0; facet < num_facets; ++facet) {
          fprintf(stdout,
                  "FACET %d: owning element = %d, First point = (%g,%g), Second point = (%g,%g)\n",
                  facet, owning_elem[facet], point0[facet][0], point0[facet][1], point1[facet][0],
                  point1[facet][1]);
        }
      }
#endif

      if (ls != NULL && ls->Evolution == LS_EVOLVE_SLAVE) {
        surf_based_initialization(x, NULL, NULL, exo, num_total_nodes, ls->init_surf_list, time1,
                                  theta, delta_t);
      }
      if (pfd != NULL) {
        ls_old = ls; /*the ol' switcheroo*/
        ls = pfd->ls[0];
        if (ls->Evolution == LS_EVOLVE_SLAVE) {
          surf_based_initialization(x, NULL, NULL, exo, num_total_nodes, ls->init_surf_list, time1,
                                    theta, delta_t);
        }
        ls = ls_old;
      }

      /* Now go back and correct all those dofs using XFEM */
      if (xfem != NULL) {
        xfem_predict(num_total_nodes, numProcUnknowns, delta_t, delta_t_old, delta_t_older, theta,
                     x, x_old, x_older, x_oldest, xdot, xdot_old, xdot_older);
      }

      /*
       * Now, that we have a predicted solution for the current
       * time, x[], exchange the degrees of freedom to update the
       * ghost node information.
       */
      exchange_dof(cx[0], dpi, x, 0);
      exchange_dof(cx[0], dpi, xdot, 0);

      if (nAC > 0) {

        if (!nonconv_roll) {
          predict_solution(nAC, delta_t, delta_t_old, delta_t_older, theta, x_AC, x_AC_old,
                           x_AC_older, x_AC_oldest, x_AC_dot, x_AC_dot_old, x_AC_dot_older);
        }

        for (iAC = 0; iAC < nAC; iAC++) {
          update_parameterAC(iAC, x, xdot, x_AC, cx[0], exo, dpi);
          augc[iAC].tmp2 = x_AC_dot[iAC];
          augc[iAC].tmp3 = x_AC_old[iAC];
        }
      }

      /*
       *  Set dirichlet conditions in some places. Note, I believe
       *  this step can change the solution vector
       */
      find_and_set_Dirichlet(x, xdot, exo, dpi);

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

      exchange_dof(cx[0], dpi, x, 0);
      exchange_dof(cx[0], dpi, xdot, 0);
      if (tran->solid_inertia)
        exchange_dof(cx[0], dpi, tran->xdbl_dot, 0);

      /*
       * Save the predicted solution for the time step
       * norm calculation to be carried out after convergence
       * of the nonlinear implicit problem
       */
      dcopy1(numProcUnknowns, x, x_pred);
      if (nAC > 0)
        dcopy1(nAC, x_AC, x_AC_pred);

#ifdef GOMA_ENABLE_OMEGA_H
      if ((tran->ale_adapt || (ls != NULL && ls->adapt)) && tran->theta != 0) {
        GOMA_EH(GOMA_ERROR, "Error theta time step parameter = %g only 0.0 supported", tran->theta);
      }
      if ((tran->ale_adapt || (ls != NULL && ls->adapt)) && pg->imtrx == 0 &&
          (nt == 0 || ((ls != NULL && nt % ls->adapt_freq == 0) ||
                       (tran->ale_adapt && nt % tran->ale_adapt_freq == 0)))) {
        if (last_adapt_nt == nt && adapt_step > 0) {
          adapt_step--;
        }
        last_adapt_nt = nt;
        adapt_mesh_omega_h(ams, exo, dpi, &x, &x_old, &x_older, &xdot, &xdot_old, &x_oldest,
                           &resid_vector, &x_update, &scale, adapt_step);
        adapt_step++;
        num_total_nodes = dpi->num_universe_nodes;
        num_total_nodes = dpi->num_universe_nodes;
        numProcUnknowns = NumUnknowns[pg->imtrx] + NumExtUnknowns[pg->imtrx];
        if (nt == 0) {
          if (ls->Num_Var_Init > 0)
            ls_var_initialization(&x, exo, dpi, cx);
        }
        x_save = realloc(x_save, sizeof(double) * numProcUnknowns);
        xdot_save = realloc(xdot_save, sizeof(double) * numProcUnknowns);
        exchange_dof(cx[0], dpi, x, 0);
        dcopy1(numProcUnknowns, x, x_old);
        dcopy1(numProcUnknowns, x, x_save);
        dcopy1(numProcUnknowns, x_old, x_older);
        dcopy1(numProcUnknowns, x_older, x_oldest);
        dcopy1(numProcUnknowns, xdot, xdot_save);
        realloc_dbl_1(&x_pred, numProcUnknowns, 0);
        realloc_dbl_1(&gvec, Num_Node, 0);
        realloc_dbl_1(&xdot_older, numProcUnknowns, 0);
        x_pred_static = x_pred;
        memset(xdot, 0, sizeof(double) * numProcUnknowns);
        memset(xdot_older, 0, sizeof(double) * numProcUnknowns);
        memset(x_pred, 0, sizeof(double) * numProcUnknowns);
        memset(resid_vector, 0, sizeof(double) * numProcUnknowns);
        memset(scale, 0, sizeof(double) * numProcUnknowns);
        memset(x_update, 0, sizeof(double) * (numProcUnknowns + numProcUnknowns));
        dcopy1(numProcUnknowns, xdot, xdot_old);
        wr_result_prelim_exo(rd, exo, ExoFileOut, gvec_elem);
        nprint = 0;
        //        (void) write_solution(ExoFileOut, resid_vector, x, x_sens_p,
        //                              x_old, xdot, xdot_old, tev, tev_post, gv,
        //                              rd, gvec, gvec_elem,
        //                              &nprint, delta_t, theta, 0, x_pp,
        //                              exo, dpi);
        //        nprint++;
        nullify_dirichlet_bcs();
        find_and_set_Dirichlet(x, xdot, exo, dpi);
      }
#endif

      numProcUnknowns = NumUnknowns[pg->imtrx] + NumExtUnknowns[pg->imtrx];
      /*
       *  Solve the nonlinear problem. If we achieve convergence,
       *  set the flag, converged, to true on return. If not
       *  set the flag to false.
       */
      dcopy1(numProcUnknowns, x, x_save);
      dcopy1(numProcUnknowns, xdot, xdot_save);
      tran->current_theta = theta;
      err = solve_nonlinear_problem(ams[JAC], x, delta_t, theta, x_old, x_older, xdot, xdot_old,
                                    resid_vector, x_update, scale, &converged, &nprint, tev,
                                    tev_post, gv, rd, gindex, p_gsize, gvec, gvec_elem, time1, exo,
                                    dpi, cx[0], n, &time_step_reform, is_steady_state, x_AC,
                                    x_AC_dot, time1, resid_vector_sens, x_sens, x_sens_p, NULL);
      if (err == -1)
        converged = FALSE;
      inewton = err;
      evpl_glob[0]->update_flag = 0;   /*See get_evp_stress_tensor for description */
      af->Sat_hyst_reevaluate = FALSE; /*See load_saturation for description*/
#ifdef RESET_TRANSIENT_RELAXATION_PLEASE
      /* Set TRUE to disable relaxation on timesteps after the first*/
      /* For transient, reset the Newton damping factors after a
       *   successful time step
       */
      if (nt > 0 && converged) {
        damp_factor2 = -1.;
        damp_factor1 = 1.0;
        no_relax_retry = 0;
      }
#endif
      if (converged)
        nonconv_roll = 0;

      /*
       * HKM -> I do not know if these operations are needed. I added
       *        an exchange of xdot[] here, because if x[] is exchanged
       *        then xdot needs to be exchanged as well.
       */

      exchange_dof(cx[0], dpi, x, 0);
      exchange_dof(cx[0], dpi, xdot, 0);

      if (!converged) {
        if (inewton < Max_Newton_Steps) {
          DPRINTF(stderr, "copying x_save for time: %g %d\n", time1, nonconv_roll);
          dcopy1(numProcUnknowns, x_save, x);
          dcopy1(numProcUnknowns, xdot_save, xdot);
        } else if (!relax_bit) {
          dcopy1(numProcUnknowns, x_save, x);
          dcopy1(numProcUnknowns, xdot_save, xdot);
        }
      }

      if (converged)
        af->Sat_hyst_reevaluate = TRUE; /*see load_saturation */

      if ((af->Sat_hyst_reevaluate) && (pmv_hyst != NULL)) {
        /* Determine what curve to follow and if switch is in order */
        err = evaluate_sat_hyst_criterion_nodal(x, xdot, exo);
      }

      /* Check element quality */
      good_mesh = element_quality(exo, x, ams[0]->proc_config);

      /*
       * Check the time step truncation error.
       * If too large, set the success_dt flag to FALSE. We will
       * then not accept the current time step, reduced delta_t,
       * and retry.
       */
      if (converged) {

        if (upd->ep[pg->imtrx][TURB_K] >= 0 || upd->ep[pg->imtrx][TURB_OMEGA] >= 0) {
          /*     Floor values to 0 */
          int floored_values = 0;
          for (int var = TURB_K; var <= TURB_OMEGA; var++) {
            for (int mn = 0; mn < upd->Num_Mat; mn++) {
              if (pd_glob[mn]->v[pg->imtrx][var]) {
                for (i = 0; i < num_total_nodes; i++) {
                  int j = Index_Solution(i, var, 0, 0, mn, pg->imtrx);

                  if (j != -1 && x[j] < 0) {
                    x[j] = 0.0;
                    floored_values++;
                  }
                }
              }
            }
          }

          int global_floored = 0;
          MPI_Allreduce(&floored_values, &global_floored, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

          if (global_floored > 0)
            P0PRINTF("Floored %d values\n", global_floored);
        }

        delta_t_new = time_step_control(delta_t, delta_t_old, const_delta_t, x, x_pred, x_old, x_AC,
                                        x_AC_pred, eps, &success_dt, tran->use_var_norm);
        if (const_delta_t) {
          success_dt = TRUE;
          delta_t_new = delta_t;
        } else if (failed_recently_countdown > 0) {
          delta_t_new = delta_t;
          failed_recently_countdown--;
        } else if (delta_t_new > delta_t_max) {
          delta_t_new = delta_t_max;
        } else if (delta_t_new < tran->resolved_delta_t_min) {
          /* fool algorithm into using delta_t = tran->resolved_delta_t_min */
          delta_t_new = tran->resolved_delta_t_min;
          success_dt = TRUE;
          DPRINTF(stderr, "\n\tminimum resolved step limit! - step control\n");
          /*     if(!success_dt)delta_t /= tran->time_step_decelerator;
               tran->delta_t  = delta_t;
               tran->delta_t_avg = 0.25*(delta_t+delta_t_old+delta_t_older
                                           +delta_t_oldest);  */
        } else {
          /* accept any converged solution with
             delta_t <= tran->resolved_delta_t_min
          */
          /*success_dt = TRUE;
          delta_t_new = delta_t;*/
        }

        if (ls != NULL && tran->Courant_Limit != 0.) {
          double Courant_dt;
          Courant_dt =
              tran->Courant_Limit * Courant_Time_Step(x, x_old, x_older, xdot, xdot_old,
                                                      resid_vector, ams[0]->proc_config, exo);
          if (Courant_dt > 0. && Courant_dt < delta_t_new) {
            DPRINTF(stdout, "\nCourant Limit requires dt <= %g\n", Courant_dt);
            delta_t_new = Courant_dt;
          }
        }
      }

      if (converged && success_dt) {
        if (Filter_Species) {
          err = filter_conc(num_total_nodes, x, filter_species_material_number, c_min, c_max);
        }
        nt += 1;
        time = time1;

        /* Determine whether to print out the data or not */
        i_print = 0;
        i_fix = 0;
        if (tran->print_freq == 0) {
          if ((time > time_print) || (fabs(time - time_print) < (1.e-4 * tran->print_delt))) {
            if (tran->print_delt2 < 0.) {
              i_print = 1;
              i_fix = 1;
              time_print += tran->print_delt;
            } else {
              if (time < tran->print_delt2_time) {
                i_print = 1;
                i_fix = 1;
                time_print += tran->print_delt;
              } else {
                i_print = 1;
                i_fix = 1;
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

        /* Check if steady state has been reached */
        if (tran->march_to_steady_state && n > 0) {
          /* for now check last two matrices */
          int steady_state_reached = TRUE;
          double max_distance = 0;

          double distance = vector_distance(NumUnknowns[0], x, x_old);
          if (distance > tran->steady_state_tolerance) {
            steady_state_reached = FALSE;
          }
          max_distance = distance;

          if (ProcID == 0) {
            printf("\nDelta x %g\n", max_distance);
          }

          if (steady_state_reached) {
            if (ProcID == 0) {
              printf("\n Steady state reached \n");
            }
            goto free_and_clear;
          }
        }

        if (time1 >= (ROUND_TO_ONE * time_max))
          i_print = 1;

        /* Dump out user specified information to separate file.
         */
        err = usr_print(&time, delta_t, x, NULL, -1);
        GOMA_EH(err, "usr_print");

        /* Particle calculations.  time = time at *beginning* of
         * current timestep, n = timestep. */
        if (Particle_Dynamics)
          err = compute_particles(exo, x, x_old, xdot, xdot_old, resid_vector, time, delta_t, n);
        GOMA_EH(err, "Error performing particle calculations.");

        if (update_etch_area && converged) {
          err = advance_etch_area_ext_field(n, num_total_nodes, delta_t, x);
          GOMA_EH(err, "Problem with advance_etch_area_ext_field()!");
        }

        error = 0;
        if (i_print) {
          if (file != NULL) {
            error = write_ascii_soln(x, resid_vector, numProcUnknowns, x_AC, nAC, time, file);
            if (error != 0) {
              fprintf(stderr, "%s:  error writing ASCII soln file\n", yo);
            }
          }
          if (Write_Intermediate_Solutions == 0) {
#ifdef LIBRARY_MODE
            lib_print = FALSE;
            if (libio->goma_first == 0) {
              if (libio->print_flag == 0 && last_call && last_step)
                lib_print = TRUE;
              if (libio->print_flag == 1 && last_step)
                lib_print = TRUE;
              if (libio->print_flag == 2)
                lib_print = TRUE;
              fprintf(stderr, "JAS FIRST...animas_step = %d, lib_print = %d\n", libio->animas_step,
                      lib_print);
            } else {
              if (libio->print_flag == 2 && !last_step)
                lib_print = TRUE;
            }

            if (lib_print)
            /* Write out the full solution */
            {
              write_solution(ExoFileOut, resid_vector, x, x_sens_p, x_old, xdot, xdot_old, tev,
                             tev_post, gv, rd, gvec, gvec_elem, &nprint, delta_t, theta, time, x_pp,
                             exo, dpi);
              nprint++;
            } else if (Num_Export_XP > 0)
            /* Just update the post-processing fields */
            {
              post_process_nodal(x, x_sens_p, x_old, xdot, xdot_old, resid_vector, nprint + 1,
                                 &time_value, delta_t, theta, x_pp, exo, dpi, rd, NULL);
            }
#else
            write_solution(ExoFileOut, resid_vector, x, x_sens_p, x_old, xdot, xdot_old, tev,
                           tev_post, gv, rd, gvec, gvec_elem, &nprint, delta_t, theta, time, x_pp,
                           exo, dpi);
            nprint++;
#endif

            if (ls != NULL && ls->Interface_Output == TRUE) {
              print_point_list(x, exo, ls->output_file, time);
            }
          }
          /* Print out values of extra unknowns from augmenting conditions */
          if (nAC > 0) {
            DPRINTF(stdout, "\n------------------------------\n");
            DPRINTF(stdout, "Augmenting Conditions:    %4d\n", nAC);
            DPRINTF(stdout, "Number of extra unknowns: %4d\n\n", nAC);

            for (iAC = 0; iAC < nAC; iAC++) {
              evol_local = augc[iAC].evol;
#ifdef PARALLEL
              if (Num_Proc > 1 &&
                  (augc[iAC].Type == AC_VOLUME || augc[iAC].Type == AC_POSITION ||
                   augc[iAC].Type == AC_ANGLE || augc[iAC].Type == AC_POSITION_MT)) {
                MPI_Allreduce(&evol_local, &evol_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                evol_local = evol_global;
              }
#endif
              if (augc[iAC].Type == AC_USERBC) {
                DPRINTF(stdout, "\tBC[%4d] DF[%4d]=% 10.6e\n", augc[iAC].BCID, augc[iAC].DFID,
                        x_AC[iAC]);
              } else if (augc[iAC].Type == AC_USERMAT || augc[iAC].Type == AC_FLUX_MAT) {
                DPRINTF(stdout, "\tMT[%4d] MP[%4d]=% 10.6e\n", augc[iAC].MTID, augc[iAC].MPID,
                        x_AC[iAC]);
              } else if (augc[iAC].Type == AC_VOLUME) {
                DPRINTF(stdout, "\tMT[%4d] VC[%4d]=%10.6e Param=%10.6e\n", augc[iAC].MTID,
                        augc[iAC].VOLID, evol_local, x_AC[iAC]);
              } else if (augc[iAC].Type == AC_FLUX) {
                DPRINTF(stdout, "\tBC[%4d] DF[%4d]=%10.6e\n", augc[iAC].BCID, augc[iAC].DFID,
                        x_AC[iAC]);
              } else if (augc[iAC].Type == AC_POSITION || augc[iAC].Type == AC_POSITION_MT) {
                DPRINTF(stdout, "\tNodeSet[%4d]_Pos = %10.6e F_bal = %10.6e MT[%4d] Param=%10.6e\n",
                        augc[iAC].MTID, evol_local, augc[iAC].lm_resid, augc[iAC].VOLID, x_AC[iAC]);
              } else if (augc[iAC].Type == AC_ANGLE) {
                evol_local = augc[iAC].lm_resid + augc[iAC].CONSTV;
                DPRINTF(stdout, "\tNodeSet[%4d]_Ang = %g F_bal = %6.3e MT[%4d] Param=%6.3e\n",
                        augc[iAC].MTID, evol_local, augc[iAC].lm_resid, augc[iAC].VOLID, x_AC[iAC]);
              }
            }
          }

        } /* if(i_print) */
        evpl_glob[0]->update_flag = 1;

        /* Fix output if current time step matches frequency */
#ifdef PARALLEL
        if ((step_fix != 0 && nt == step_fix) || ((i_fix == 1) && (tran->fix_freq > 0))) {
          /* Barrier because fix needs both files to be finished printing
             and fix always occurs on the same timestep as printing */
          MPI_Barrier(MPI_COMM_WORLD);
          fix_output();
          /* Fix step is relative to print step */
          step_fix += tran->fix_freq * tran->print_freq;
        }
#endif
        /*
         * Adjust the time step if the new time will be larger than the
         * next printing time.
         */
        if (tran->print_freq == 0 && success_dt) {
          if ((time + 1.2 * delta_t_new >= time_print) && (time_print > time)) {
            delta_t_new = time_print - time;
            DPRINTF(stdout, "reset delta_t = %g to maintain printing frequency\n", delta_t_new);
            if (delta_t_new <= 0)
              GOMA_EH(GOMA_ERROR, "error with time-step printing control");
          } else if (time >= time_print) {
            if (delta_t_new != tran->print_delt) {
              delta_t_new = tran->print_delt;
              DPRINTF(stdout, "reset delta_t = %g to maintain printing frequency\n", delta_t_new);
              if (delta_t_new <= 0) {
                GOMA_EH(GOMA_ERROR, "error with time-step printing control");
              }
            }
          }
        }

#ifdef LIBRARY_MODE
        /*
         * When porosity is passed in from another code, do a
         * routine update to the external field.
         */
        if (update_porosity && converged) {
          error = advance_porosity_ev(n, num_total_nodes, x, base_p_por, base_p_liq);
          GOMA_EH(error, "Problem with advance_porosity_ev()!");
        }
#endif

        if (converged && ls != NULL) {
          int ibc, ls_adc_event = FALSE;
          /* Resolve LS_ADC boundaries ( Attach/Dewet/Coalesce ) */
          for (ibc = 0; ibc < Num_BC; ibc++) {
            ls_adc_event = FALSE;
            switch (BC_Types[ibc].BC_Name) {
            case LS_ADC_OLD_BC:
              resolve_ls_adc_old(&(BC_Types[ibc]), exo, x, delta_t, &ls_adc_event, nt);
              break;
            case LS_ADC_BC:
              resolve_ls_adc(ls->last_surf_list->start->subsurf_list, &(BC_Types[ibc]), exo, x,
                             delta_t, &ls_adc_event, nt);

              break;
            default:
              break;
            }
#ifdef PARALLEL
            if (ls_adc_event) {
              exchange_dof(cx[0], dpi, x, 0);
            }
#endif
            if (ls_adc_event && tran->Restart_Time_Integ_After_Renorm) {
              /* like a restart */
              discard_previous_time_step(numProcUnknowns, x, x_old, x_older, x_oldest, xdot,
                                         xdot_old, xdot_older);
              last_renorm_nt = nt;
              if (delta_t_new > fabs(delta_t0))
                delta_t_new *= tran->time_step_decelerator;
            }
          }

          /* Check for renormalization  */

          ls->Renorm_Countdown -= 1;

          switch (ls->Renorm_Method) {

          case HUYGENS:
          case HUYGENS_C:
          case HUYGENS_MASS_ITER:
            Renorm_Now =
                (ls->Renorm_Freq != 0 && ls->Renorm_Countdown == 0) || ls_adc_event == TRUE;

            did_renorm =
                huygens_renormalization(x, num_total_nodes, exo, cx[0], dpi, num_fill_unknowns,
                                        numProcUnknowns, time2, Renorm_Now);
            if (did_renorm) {
              exchange_dof(cx[0], dpi, x, 0);
            }
            if (did_renorm && tran->Restart_Time_Integ_After_Renorm) {
              /* like a restart */
              discard_previous_time_step(numProcUnknowns, x, x_old, x_older, x_oldest, xdot,
                                         xdot_old, xdot_older);
              last_renorm_nt = nt;
              if (delta_t_new > fabs(delta_t0))
                delta_t_new *= tran->time_step_decelerator;
            }

            break;

          case CORRECT:
            GOMA_EH(GOMA_ERROR, "Use of \"CORRECT\" is obsolete.");
            break;
          default:
            break;
          }

          if (ls->Sat_Hyst_Renorm_Lockout > 0) {
            af->Sat_hyst_reevaluate = FALSE;
            ls->Sat_Hyst_Renorm_Lockout -= 1;
          }
        }

        if (converged && pfd != NULL) {
          struct Level_Set_Data *ls_save = ls;

          for (i = 0; i < pfd->num_phase_funcs; i++) {

            ls = pfd->ls[i];
            ls->Renorm_Countdown -= 1;
            switch (ls->Renorm_Method) {

            case HUYGENS:
            case HUYGENS_C:
            case HUYGENS_MASS_ITER:
              Renorm_Now = (ls->Renorm_Freq != 0 && ls->Renorm_Countdown == 0);

              did_renorm =
                  huygens_renormalization(x, num_total_nodes, exo, cx[0], dpi, num_fill_unknowns,
                                          numProcUnknowns, time2, Renorm_Now);
              if (did_renorm) {
                exchange_dof(cx[0], dpi, x, 0);
              }
              if (did_renorm && tran->Restart_Time_Integ_After_Renorm) {
                /* like a restart */
                discard_previous_time_step(numProcUnknowns, x, x_old, x_older, x_oldest, xdot,
                                           xdot_old, xdot_older);
                last_renorm_nt = nt;
                if (delta_t_new > fabs(delta_t0))
                  delta_t_new *= tran->time_step_decelerator;
              }
              break;
            default:
              break;
            }
          }

          ls = ls_save;
        }
        /*
         *   save xdot to xdot_old for next time step
         */
        dcopy1(numProcUnknowns, xdot_old, xdot_older);
        if (tran->solid_inertia)
          dcopy1(numProcUnknowns, tran->xdbl_dot, tran->xdbl_dot_old);
        dcopy1(numProcUnknowns, xdot, xdot_old);
        dcopy1(numProcUnknowns, x_older, x_oldest);
        dcopy1(numProcUnknowns, x_old, x_older);
        dcopy1(numProcUnknowns, x, x_old);
        delta_t_oldest = delta_t_older;
        delta_t_older = delta_t_old;
        delta_t_old = delta_t;
        tran->delta_t_old = delta_t_old;
        tran->time_value_old = time;
        delta_t = delta_t_new;
        tran->delta_t = delta_t; /*load up for use in load_fv_mesh_derivs*/
        tran->delta_t_avg = 0.25 * (delta_t + delta_t_old + delta_t_older + delta_t_oldest);

        if (nAC > 0) {
          dcopy1(nAC, x_AC_dot_old, x_AC_dot_older);
          dcopy1(nAC, x_AC_dot, x_AC_dot_old);
          dcopy1(nAC, x_AC_older, x_AC_oldest);
          dcopy1(nAC, x_AC_old, x_AC_older);
          dcopy1(nAC, x_AC, x_AC_old);
        }

        /* Integrate fluxes, forces
         */
        for (i = 0; i < nn_post_fluxes; i++) {
          (void)evaluate_flux(exo, dpi, pp_fluxes[i]->ss_id, pp_fluxes[i]->flux_type,
                              pp_fluxes[i]->flux_type_name, pp_fluxes[i]->blk_id,
                              pp_fluxes[i]->species_number, pp_fluxes[i]->flux_filenm,
                              pp_fluxes[i]->profile_flag, x, xdot, NULL, delta_t_old, time, 1);
        }

        /* Compute flux, force sensitivities
         */
        for (i = 0; i < nn_post_fluxes_sens; i++) {
          (void)evaluate_flux_sens(exo, dpi, pp_fluxes_sens[i]->ss_id, pp_fluxes_sens[i]->flux_type,
                                   pp_fluxes_sens[i]->flux_type_name, pp_fluxes_sens[i]->blk_id,
                                   pp_fluxes_sens[i]->species_number, pp_fluxes_sens[i]->sens_type,
                                   pp_fluxes_sens[i]->sens_id, pp_fluxes_sens[i]->sens_flt,
                                   pp_fluxes_sens[i]->sens_flt2, pp_fluxes_sens[i]->vector_id,
                                   pp_fluxes_sens[i]->flux_filenm, pp_fluxes_sens[i]->profile_flag,
                                   x, xdot, x_sens_p, delta_t_old, time, 1);
        }

        if (Output_Variable_Stats) {
          if (time >= time_max) {
            err = variable_stats(x, time, Output_Variable_Regression);
          } else {
            err = variable_stats(x, time, FALSE);
          }
          GOMA_EH(err, "Problem with variable_stats!");
          if (ProcID == 0)
            fflush(stdout);
        }

#ifdef REACTION_PRODUCT_EFV
        for (i = 0; i < exo->num_nodes; i++) {
          if (efv->ev && nt > 1) {
            int ef;
            if (fabs(Spec_source_lumped_mass[i]) > DBL_SMALL) {
              for (ef = 0; ef < efv->Num_external_field; ef++) {
                efv->ext_fld_ndl_val[ef][i] *= Spec_source_lumped_mass[i];
              }
            }
          }
        }
        memset(Spec_source_lumped_mass, 0.0, sizeof(double) * exo->num_nodes);
#endif
        for (i = 0; i < nn_volume; i++) {
          evaluate_volume_integral(exo, dpi, pp_volume[i]->volume_type, pp_volume[i]->volume_name,
                                   pp_volume[i]->blk_id, pp_volume[i]->species_no,
                                   pp_volume[i]->volume_fname, pp_volume[i]->params,
                                   pp_volume[i]->num_params, NULL, x, xdot, delta_t_old, time, 1);
        }

#ifdef REACTION_PRODUCT_EFV
        for (i = 0; i < exo->num_nodes; i++) {
          if (efv->ev && nt > 1) {
            int ef;
            for (ef = 0; ef < efv->Num_external_field; ef++) {
              if (fabs(Spec_source_lumped_mass[i]) > DBL_SMALL) {
                efv->ext_fld_ndl_val[ef][i] /= Spec_source_lumped_mass[i];
              } else {
                efv->ext_fld_ndl_val[ef][i] = 0.0;
              }
            }
          }
        }
        if (efv->ev && i_print) {
          error = 0;
          if (file != NULL) {
            error = write_ascii_soln(x, resid_vector, numProcUnknowns, x_AC, nAC, time, file);
            if (error != 0) {
              fprintf(stderr, "%s:  error writing ASCII soln file\n", yo);
            }
          }
          if (Write_Intermediate_Solutions == 0) {
            write_solution(ExoFileOut, resid_vector, x, x_sens_p, x_old, xdot, xdot_old, tev,
                           tev_post, gv, rd, gvec, gvec_elem, &nprint, delta_t, theta, time, x_pp,
                           exo, dpi);
            nprint++;
          }
        }
#endif

        if (time1 >= (ROUND_TO_ONE * time_max)) {
          DPRINTF(stdout, "\t\tout of time!\n");
          if (Anneal_Mesh) {
            /*
             * Transform the node point coordinates according to the
             * displacements and write out all the results using the
             * displaced coordinates. Set the displacement field to
             * zero, too.
             */
            err = anneal_mesh(x, tev, tev_post, gv, rd, time1, exo, dpi);
            GOMA_EH(err, "anneal_mesh() bad return.");
          }
          goto free_and_clear;
        }
        if (!good_mesh)
          goto free_and_clear;

      } /*  if(converged && success_dt) */

      else /* not converged or unsuccessful time step */
      {
        if (relax_bit && (nonconv_roll < no_relax_retry)) {
          /*success_dt = TRUE;  */
#ifdef RESET_TRANSIENT_RELAXATION_PLEASE
          nonconv_roll++;
          if (nonconv_roll == 1 && nt != 0) {
            damp_factor1 = damp_factor_org[0];
            damp_factor2 = damp_factor_org[1];
            custom_tol1 = toler_org[0];
            custom_tol2 = toler_org[1];
            custom_tol3 = toler_org[2];
          }
#endif
          if (inewton == -1) {
            DPRINTF(stdout, "\nHmm... trouble on this step \n  Let's try some more relaxation %d\n",
                    no_relax_retry - nonconv_roll);
            if (use_custom_damp) {
              custom_tol1 *= 0.01;
              custom_tol2 *= 0.01;
              custom_tol3 *= 0.01;
              DPRINTF(stdout, "  custom tolerances %g %g %g  \n", custom_tol1, custom_tol2,
                      custom_tol3);
            } else {
              damp_factor1 *= 0.5;
              DPRINTF(stdout, "  damping factor %g  \n", damp_factor1);
            }
          } else if (!converged && (inewton == Max_Newton_Steps)) {
            DPRINTF(stdout,
                    "\nHmm... could not converge on this step\nLet's try some more iterations %d\n",
                    no_relax_retry - nonconv_roll);
            dcopy1(numProcUnknowns, x, x_old);
            if (nAC > 0) {
              dcopy1(nAC, x_AC, x_AC_old);
            }
            if (use_custom_damp) {
              if (nonconv_roll > 1 || nt == 0) {
                custom_tol1 *= 100.;
                custom_tol2 *= 100.;
                custom_tol3 *= 100.;
                DPRINTF(stdout, "  custom tolerances %g %g %g  \n", custom_tol1, custom_tol2,
                        custom_tol3);
              }
            } else {
              if (nonconv_roll > 1 || nt == 0) {
                damp_factor1 *= 2.0;
                damp_factor1 = MIN(damp_factor1, 1.0);
              }
              DPRINTF(stdout, "  damping factor %g  \n", damp_factor1);
            }
          } else {
            DPRINTF(stderr, "\n\tlast time step failed, dt *= %g for next try!\n",
                    tran->time_step_decelerator);

            damp_factor1 = MAX(damp_factor1, 1.0);
            delta_t *= tran->time_step_decelerator;
            tran->delta_t = delta_t;
            tran->delta_t_avg = 0.25 * (delta_t + delta_t_old + delta_t_older + delta_t_oldest);
            time1 = time + delta_t;
            tran->time_value = time1;
            evpl_glob[0]->update_flag = 2;
            af->Sat_hyst_reevaluate = 0;

            /* if specified with "Steps of constant delta_t after failure"
               use a constant delta_t to help the painful recovery
             */
            failed_recently_countdown = tran->const_dt_after_failure;
          }
        } else if (converged &&
                   delta_t < tran->resolved_delta_t_min / tran->time_step_decelerator) {
          DPRINTF(stderr, "\n\tminimum resolved step limit! - not converged\n");
          delta_t_oldest = delta_t_older;
          delta_t_older = delta_t_old;
          delta_t_old = delta_t;
          tran->delta_t_old = delta_t_old;
          tran->time_value_old = time;
          delta_t = tran->resolved_delta_t_min;
          tran->delta_t = delta_t;
          tran->delta_t_avg = 0.25 * (delta_t + delta_t_old + delta_t_older + delta_t_oldest);
          time1 = time + delta_t;
          tran->time_value = time1;
        } else {
          DPRINTF(stdout, "\n\tlast time step failed, dt *= %g for next try!\n",
                  tran->time_step_decelerator);

          delta_t *= tran->time_step_decelerator;
          tran->delta_t = delta_t;
          tran->delta_t_avg = 0.25 * (delta_t + delta_t_old + delta_t_older + delta_t_oldest);
          time1 = time + delta_t;
          tran->time_value = time1;
          evpl_glob[0]->update_flag = 2;
          af->Sat_hyst_reevaluate = 0;

          /* if specified with "Steps of constant delta_t after failure"
             use a constant delta_t to help the painful recovery
           */
          failed_recently_countdown = tran->const_dt_after_failure;
        }

#define DEBUG_SUBELEMENT_CONFIGURATION 0
#if DEBUG_SUBELEMENT_CONFIGURATION
        if (ls != NULL && ls->SubElemIntegration) {
          DPRINTF(stdout, "predicted subelement configuration:\n");
          subelement_mesh_output(x_pred, exo);
          DPRINTF(stdout, "failed solution subelement configuration:\n");
          subelement_mesh_output(x, exo);
        }
#endif
      }

      if (delta_t <= delta_t_min) {
        DPRINTF(stderr, "\n\tdelta_t = %e < %e\n\n", delta_t, delta_t_min);

        DPRINTF(stdout, "time step too small, I'm giving up!\n");
        break;
      }

    } /* end of time step loop */

    /*  I moved this section up so that it is actually executed - RBS
  if (Anneal_Mesh)
    {
     * Transform the node point coordinates according to the
     * displacements and write out all the results using the
     * displaced coordinates. Set the displacement field to
     * zero, too.
      err = anneal_mesh(x, tev, tev_post, gv, rd, time_value, exo, dpi);
      GOMA_EH(err, "anneal_mesh() bad return.");
    }
     */

  } /* end of if steady else transient */

free_and_clear:

/* If exporting variables to another code, save them now! */
#ifdef LIBRARY_MODE
  callnum++;
  *te_out = time;
  if (Num_Export_XS > 0 || Num_Export_XP > 0) {
    err = load_export_vars(num_total_nodes, x, x_pp);
    GOMA_EH(err, "Problem with saving export variables");
  }
#endif

  /* Free a bunch of variables that aren't needed anymore */

  for (i = 0; i < Num_ROT; i++) {
    safer_free((void **)&(ROT_Types[i].elems));
  }

  /*
   * Curiously, these bananas were allocated end to end - one free does it all.
   */

  safer_free((void **)&ROT_Types);
  safer_free((void **)&node_to_fill);

  safer_free((void **)&resid_vector);
  safer_free((void **)&resid_vector_sens);
  safer_free((void **)&scale);
  if (last_call)
    safer_free((void **)&x);
  if (last_call)
    safer_free((void **)&xp_id);
  safer_free((void **)&gv);
  safer_free((void **)&x_pp);
  if (update_porosity) {
    safer_free((void **)&base_p_por);
    safer_free((void **)&base_p_liq);
  }

  if (nAC > 0 && last_call) {
    safer_free((void **)&x_AC);
    safer_free((void **)&x_AC_old);
    safer_free((void **)&x_AC_older);
    safer_free((void **)&x_AC_oldest);
    safer_free((void **)&x_AC_dot);
    safer_free((void **)&x_AC_dot_old);
    safer_free((void **)&x_AC_dot_older);
    safer_free((void **)&x_AC_pred);
  }

  safer_free((void **)&x_pred);

  if (last_call) {
    safer_free((void **)&x_save);
    safer_free((void **)&xdot_save);
    safer_free((void **)&x_old);
    safer_free((void **)&x_older);
    safer_free((void **)&x_oldest);
    safer_free((void **)&xdot);
    safer_free((void **)&xdot_old);
    safer_free((void **)&xdot_older);
  }

  if (tran->solid_inertia) {
    free(tran->xdbl_dot);
    free(tran->xdbl_dot_old);
  }
  safer_free((void **)&x_update);

  Dmatrix_death(Spec_source_inventory, upd->Num_Mat, upd->Max_Num_Species_Eqn + 1);
  safer_free((void **)&Spec_source_lumped_mass);
  if ((nn_post_data_sens + nn_post_fluxes_sens) > 0) {
    safer_free((void **)&x_sens);
    Dmatrix_death(x_sens_p, num_pvector, numProcUnknowns);
  }

  if (last_call) {
    for (i = 0; i < MAX_NUMBER_MATLS; i++) {
      for (n = 0; n < MAX_MODES; n++) {
        safer_free((void **)&(ve_glob[i][n]->gn));
        safer_free((void **)&(ve_glob[i][n]->time_const_st));
        safer_free((void **)&(ve_glob[i][n]));
      }
      safer_free((void **)&(vn_glob[i]));
    }
  }

  safer_free((void **)&ve);

  if (num_universe_dofs[pg->imtrx] !=
      (num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx])) {
    safer_free((void **)&ija_attic);
  }

  if (last_call) {
    GomaSparseMatrix goma_matrix = ams[JAC]->GomaMatrixData;
    GomaSparseMatrix_Destroy(&goma_matrix);
#ifdef GOMA_ENABLE_PETSC
    if (strcmp(Matrix_Format, "petsc") == 0) {
      err = goma_petsc_free_matrix(ams[JAC]);
      GOMA_EH(err, "free petsc matrix");
    }
#endif

    sl_free(matrix_systems_mask, ams);

    for (i = 0; i < NUM_ALSS; i++) {
      safer_free((void **)&(ams[i]));
    }

    safer_free((void **)&gvec);

    i = 0;
    for (eb_indx = 0; eb_indx < exo->num_elem_blocks; eb_indx++) {
      for (ev_indx = 0; ev_indx < rd->nev; ev_indx++) {
        if (exo->elem_var_tab[i++] == 1) {
          safer_free((void **)&(gvec_elem[eb_indx][ev_indx]));
        }
      }
      safer_free((void **)&(gvec_elem[eb_indx]));
    }

    safer_free((void **)&gvec_elem);
    safer_free((void **)&rd);
    for (int i = 0; i < num_total_nodes; i++) {
      free(Local_Offset[0][i]);
      free(Dolphin[0][i]);
    }
    free(Dolphin[0]);
    free(Local_Offset[0]);
    safer_free((void **)&Local_Offset);
    safer_free((void **)&Dolphin);
  }

  if (ls != NULL) {
    if (ls->init_surf_list != NULL)
      free_surf_list(&(ls->init_surf_list));
  }

  free_shape_fcn_tree(Subgrid_Tree);

  if (file != NULL)
    fclose(file);

  return;
} /* END of routine solve_problem()  */
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

void predict_solution(int N,
                      double delta_t,
                      double delta_t_old,
                      double delta_t_older,
                      double theta_arg,
                      double x[],
                      double x_old[],
                      double x_older[],
                      double x_oldest[],
                      double xdot[],
                      double xdot_old[],
                      double xdot_older[])

/*
 *    Function to calculate the predicted solution vector, x_pred_n for the
 *    (n+1)th time step.
 *    This routine can be used by a first order - forward Euler / backward
 *    Euler predictor / corrector method or for a second order
 *    Adams-Bashforth / Trapezoidal Rule predictor / corrector method.
 *           See Nachos
 *    documentation Sand86-1816 and Gresho, Lee, Sani LLNL report UCRL - 83282
 *    for more information.
 *
 *    9/28/94 RRR - It is now modified so it will be a variable order predictor
 *    with the trapezoid rule being 2nd order and everything else being somewhat
 *    lower all the way down to first order.
 *
 *  variables:
 *
 *    on input:
 *
 *       N          - number of unknowns
 *      delta_t     - magnitude of time step at time n     (i.e., = t_n+1 - t_n)
 *      delta_t_old - magnitude of time step at time n - 1 (i.e., = t_n - t_n-1)
 *      delta_t_older - magnitude of time step at time n - 2 (i.e., = t_n-1 - t_n-2)
 *      x_old[]     - solution vector at time n
 *      x_older[]     - solution vector at time n-1
 *      x_oldest[]     - solution vector at time n-2
 *      xdot_old[]      - acceleration vector from the predictor at time n
 *      xdot_older[]  - acceleration vector from the predictor at time n - 1
 *
 *   on output:
 *
 *      x[]         - predicted solution vector at time n + 1
 *
 */
{
  int i;
  double c1, c2, c3;

  if (theta_arg == 0.5) {
    c1 = delta_t * (delta_t + delta_t_old) / delta_t_older / (delta_t_older + delta_t_old);
    c2 = -delta_t * (delta_t + delta_t_old + delta_t_older) / (delta_t_old * delta_t_older);
    c3 = (delta_t + delta_t_old + delta_t_older) * (delta_t + delta_t_old) / delta_t_old /
         (delta_t_older + delta_t_old);

    for (i = 0; i < N; i++) {
      x[i] = c3 * x_old[i] + c2 * x_older[i] + c1 * x_oldest[i];
    }
  } else {
    c1 = delta_t * (1.0 + theta_arg * delta_t / delta_t_old);
    c2 = theta_arg * (delta_t * delta_t) / (delta_t_old);
    for (i = 0; i < N; i++) {
      x[i] = x_old[i] + c1 * xdot_old[i] - c2 * xdot_older[i];
    }
  }

  /*
   * Apply any additional intelligence as to what the predicted solution
   * x[] should be equal to here. For example, if there any stability or
   * phase boundary constraints, they should be installed at this location
   * via a function call. Also, one might want to implement time dependent
   * Dirichlet conditions here.
   */

  /*
   * Update xdot[] to be consistent with the predicted x[] calculated
   * above
   */
  for (i = 0; i < N; i++) {
    xdot[i] =
        (1.0 + 2.0 * theta_arg) / delta_t * (x[i] - x_old[i]) - (2.0 * theta_arg) * xdot_old[i];
  }

} /* END of routine predict_solution  */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/**************************************************************************/

static void predict_solution_newmark(
    int N, double delta_t, double x[], double x_old[], double xdot[], double xdot_old[])

/*
 *    Function to calculate the predicted solution vector, x_pred_n for the
 *    (n+1)th time step.
 *    This routine is set up for the newmark-beta scheme as described by
 *    TJR Hughes' book on FEM.    Note that only MESH displacements are
 *    updated in DYNAMIC_LAGRANGIAN mesh motion materials.  One could
 *    add to this if need be.
 *
 *       N          - number of Nodes
 *      delta_t     - magnitude of time step at time n     (i.e., = t_n+1 - t_n)
 *      x_old[]     - solution vector at time n
 *      xdot_old[]      - acceleration vector from the predictor at time n
 *
 *   on output:
 *
 *      x[]         - predicted solution vector at time n + 1
 *
 */
{
  int i, j, imat, num_mat, mat_index, a;
  double c1, c2, c3;
  UMI_LIST_STRUCT *curr_mat_list;
  NODE_INFO_STRUCT *node_ptr;

  c1 = delta_t;
  c2 = delta_t * delta_t * (1. - 2. * tran->newmark_beta) / 2.;
  c3 = delta_t * (1. - tran->newmark_gamma);
  for (i = 0; i < N; i++) {

    node_ptr = Nodes[i];
    curr_mat_list = &(node_ptr->Mat_List);
    num_mat = curr_mat_list->Length;

    for (imat = 0; imat < num_mat; imat++) {
      mat_index = (curr_mat_list->List)[imat];
      if (pd_glob[mat_index]->MeshMotion == DYNAMIC_LAGRANGIAN) {
        for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {
          j = Index_Solution(i, R_MESH1 + a, 0, 0, -1, pg->imtrx);
          x[j] = x_old[j] + c1 * xdot_old[j] + c2 * tran->xdbl_dot_old[j];
          xdot[j] = xdot_old[j] + c3 * tran->xdbl_dot_old[j];
        }
      }
    }
  }
} /* END of routine predict_solution_newmark  */

int discard_previous_time_step(int num_unks,
                               double *x,
                               double *x_old,
                               double *x_older,
                               double *x_oldest,
                               double *xdot,
                               double *xdot_old,
                               double *xdot_older) {

  dcopy1(num_unks, x, x_old);
  dcopy1(num_unks, x_old, x_older);
  dcopy1(num_unks, x_older, x_oldest);

  /* also need to kill xdot(s) */
  memset(xdot, 0, sizeof(double) * num_unks);
  memset(xdot_old, 0, sizeof(double) * num_unks);
  memset(xdot_older, 0, sizeof(double) * num_unks);

  return (0);
}

/*****************************************************************************/
/*****************************************************************************/
/****************************************************************************/

int anneal_mesh(double x[],
                int tev,
                int tev_post,
                double *glob_vars_val,
                struct Results_Description *rd,
                double time_value,
                Exo_DB *exo,
                Dpi *dpi)

/* anneal_mesh -- output results with displaced coordinates
 *
 *    Function to dump an exoII output file with the mesh annealed:
 *    i.e. the nodal coordinates updated with
 *    (reference state + displacement
 *
 * Created: 1997/08 Randy Schunk
 *
 * Revised: 1997/09/22 15:52 MDT pasacki@sandia.gov
 */
{
  int dim;
  int displacement_somewhere; /* boolean */
  int i;
  int j;
  int m; /* matl index */
  int num_nodes, rd_nnv_save, rd_nev_save;
  int p;
  int ielem;
  int imtrx;
  int e_start, e_end;
  int *moved;
  /*  int num_local_nodes; */
  int gnn;
  int ln;
  int var;
  int dofs[DIM];
  /*  int dof_list[MAX_VARIABLE_TYPES+MAX_CONC][MDE]; */

  double phi[MDE];
  double displacement[DIM];
  double ***gvec_elem;
  double *x_file;
  double d[DIM][MDE];

  double *savex = NULL, *savey = NULL, *savez = NULL; /* temporary placeholders while
                                                       * anneal coordinates are
                                                       * written...                */
  double **new_coord;

  double *nodal_result_vector; /* temporarily hold one nodal variable
                                * prior to writing out into EXODUS II file */
  char afilename[MAX_FNL];

  FILE *anneal_dat;

  int numProcUnknowns = NumUnknowns[pg->imtrx] + NumExtUnknowns[pg->imtrx];

  asdv(&x_file, numProcUnknowns);

  strcpy(afilename, anneal_file);

  dcopy1(numProcUnknowns, x, x_file);

  /*
   * Return immediately if there are no mesh displacements.
   */

  displacement_somewhere = FALSE;

  for (m = 0; m < upd->Num_Mat; m++) {
    for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
      displacement_somewhere |= (pd_glob[m]->e[imtrx][R_MESH1]);
    }
  }

  if (!displacement_somewhere) {
    GOMA_WH(-1, "Attempt to anneal w/ no active displacement eqn anywhere!");
    return (0);
  }

  dim = exo->num_dim;
  num_nodes = exo->num_nodes;

  new_coord = (double **)calloc(dim, sizeof(double *));

  for (p = 0; p < dim; p++) {
    new_coord[p] = (double *)calloc(num_nodes, sizeof(double));
  }

  if (dim > 0) {
    dcopy1(exo->base_mesh->num_nodes, exo->base_mesh->x_coord, new_coord[0]);
  }

  if (dim > 1) {
    dcopy1(exo->base_mesh->num_nodes, exo->base_mesh->y_coord, new_coord[1]);
  }

  if (dim > 2) {
    dcopy1(exo->base_mesh->num_nodes, exo->base_mesh->z_coord, new_coord[2]);
  }

  nodal_result_vector = (double *)calloc(num_nodes, sizeof(double));

  moved = (int *)calloc(num_nodes, sizeof(int));

  memset(moved, 0, sizeof(int) * num_nodes);

  /*
   * Loop through nodes, find displacement, and add it into
   * the coordinate
   */

  e_start = exo->eb_ptr[0];
  e_end = exo->eb_ptr[exo->num_elem_blocks];

  for (ielem = e_start; ielem < e_end; ielem++) {
    load_elem_dofptr(ielem, exo, x, x_file, x, x, 1);

    memset(d, 0, sizeof(double) * DIM * MDE);
    memset(dofs, 0, sizeof(int) * DIM);

    for (ln = 0; ln < ei[pg->imtrx]->num_local_nodes; ln++) {
      double factor = 1.0;
      double xi[3] = {0.0, 0.0, 0.0};

      find_nodal_stu(ln, ei[pg->imtrx]->ielem_type, xi, xi + 1, xi + 2);

      gnn = exo->elem_node_list[exo->elem_node_pntr[ielem] + ln];

      if (exo->ghost_node_to_base[gnn] == -1)
        continue;

      memset(displacement, 0, sizeof(double) * DIM);

      if (moved[gnn] != 1) {
        for (p = 0; p < DIM; p++) {
          var = MESH_DISPLACEMENT1 + p;

          for (i = 0; i < ei[pg->imtrx]->dof[var]; i++) {
            phi[i] = newshape(xi, ei[pg->imtrx]->ielem_type, PSI, ei[pg->imtrx]->dof_list[var][i],
                              ei[pg->imtrx]->ielem_shape, pd->i[pg->imtrx][var], i);
          }

          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              displacement[p] += *esp->d[p][j] * phi[j];
              *esp_old->d[p][j] = 0.0;
            }
            moved[gnn] = 1;
          } else
            displacement[p] = 0.0;
        }
      }

      for (p = 0; p < dim; p++)
        new_coord[p][exo->ghost_node_to_base[gnn]] += factor * displacement[p];
    }
  }

  /*
   * Temporarily save the old coordinates and assign the new
   * displaced coordinates for the annealed file...
   */

  if (dim > 0) {
    savex = exo->base_mesh->x_coord;
    exo->base_mesh->x_coord = &new_coord[0][0];
  }

  if (dim > 1) {
    savey = exo->base_mesh->y_coord;
    exo->base_mesh->y_coord = &new_coord[1][0];
  }

  if (dim > 2) {
    savez = exo->base_mesh->z_coord;
    exo->base_mesh->z_coord = &new_coord[2][0];
  }

  if (Num_Proc > 1)
    multiname(afilename, ProcID, Num_Proc);

  one_base(exo, Num_Proc); /* cause node numbers, etc. to start at 1 */

  wr_mesh_exo(exo, afilename, 0);

  /*
   * Return internal EXODUS II data to 0 based node and element numbers
   * for convenience.
   */
  zero_base(exo);

  /*
   * Grab back the undisplaced coordinates that were originally read.
   */

  if (dim > 0)
    exo->base_mesh->x_coord = savex;

  if (dim > 1)
    exo->base_mesh->y_coord = savey;

  if (dim > 2)
    exo->base_mesh->z_coord = savez;

  if ((tev + tev_post) != 0) {
    gvec_elem = (double ***)array_alloc(2, exo->num_elem_blocks, tev + tev_post, sizeof(double *));
  } else {
    gvec_elem = NULL;
  }

  /*
   * Write the names of the results variables to the annealed file.
   */
  rd_nnv_save = rd->nnv;
  rd_nev_save = rd->nev;
  rd->nnv = rd->TotalNVSolnOutput;
  rd->nev = tev;

  wr_result_prelim_exo(rd, exo, afilename, gvec_elem);

  rd->nnv = rd_nnv_save; /* Return values to normal */
  rd->nev = rd_nev_save;

  if (Num_Proc > 1)
    wr_dpi(dpi, afilename);

  for (i = 0; i < rd->TotalNVSolnOutput; i++) {
    /*
     * if nodal variable is a displacement, set it to zero
     */
    if (rd->nvtype[i] == MESH_DISPLACEMENT1 || rd->nvtype[i] == MESH_DISPLACEMENT2 ||
        rd->nvtype[i] == MESH_DISPLACEMENT3) {
      init_vec_value(nodal_result_vector, 0.0, num_nodes);
    } else {
      extract_nodal_vec(x, rd->nvtype[i], rd->nvkind[i], rd->nvmatID[i], nodal_result_vector, exo,
                        FALSE, time_value);
    }
    wr_nodal_result_exo(exo, afilename, nodal_result_vector, i + 1, 1, time_value);
  }

  /* Now pick up the element variables */
  for (i = 0; i < tev; i++) {
    wr_elem_result_exo(exo, afilename, gvec_elem, i, 1, time_value, rd);
  }

  /* And finally the global variables */

  if (rd->ngv > 0) {
    wr_global_result_exo(exo, afilename, 1, rd->ngv, glob_vars_val);
  }

  if (Num_Proc == 1) {
    anneal_dat = fopen("anneal.dat", "w");

    (void)write_ascii_soln(x_file, NULL, numProcUnknowns, x, 0, 0.0, anneal_dat);

    fclose(anneal_dat);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (!Skip_Fix && Num_Proc > 1 && ProcID == 0) {
    DPRINTF(stdout, "\nFixing exodus file %s\n", anneal_file);
    fix_exo_file(Num_Proc, anneal_file);
  }

  /*
   * Free up any temporarily allocated memory.
   */

  for (p = 0; p < dim; p++) {
    safer_free((void **)&(new_coord[p]));
  }
  safer_free((void **)&new_coord);
  safer_free((void **)&nodal_result_vector);
  safer_free((void **)&x_file);
  safer_free((void **)&moved);
  safer_free((void **)&gvec_elem);

  return (0);
} /* END of routine anneal_mesh */
/*****************************************************************************/

/*****************************************************************************/
/*****************************************************************************/

int anneal_mesh_with_external_field(const Exo_DB *exo)

/* anneal_mesh_with_external_field --
 *         anneal mesh before solving problem with displaced coordinates
 *
 * Created: 2004/01 Randy Schunk
 *
 */
{
  int dim;
  int i;
  int num_nodes;
  int p;

  double displacement[DIM], X_old[DIM], X_new[DIM];
  double **new_coord;

  dim = exo->num_dim;
  num_nodes = exo->num_nodes;

  new_coord = (double **)calloc(dim, sizeof(double *));

  for (p = 0; p < dim; p++) {
    new_coord[p] = (double *)calloc(num_nodes, sizeof(double));
    dcopy1(num_nodes, Coor[p], new_coord[p]);
  }

  for (i = 0; i < num_nodes; i++) {
    /*
     * Default is full anneal: X_new = X_old + displacement
     */

    X_old[0] = exo->x_coord[i];
    X_old[1] = exo->y_coord[i];
    if (dim > 2) {
      X_old[2] = exo->z_coord[i];
    }
    for (p = 0; p < dim; p++) {

      /* Special case: this material should not have displacement */
      /* equations defined for it because we are    */
      /* bringing them in from external field vars  */
      /* We are assuming that the interpolation level of the efv  */
      /* is the same as that of the goma transport equation, which*/
      /* for now is just porosity.  We are annealing the mesh     */
      /* because we want to have the right geometry for evaluating*/
      /* all geometry-dependent building blocks for the Galerkin  */
      /* weighted residual, viz. grad operators, mapping for gauss*/
      /*integration, etc.                                        */

      /*N.B. assumes that the first external field variables in the
       *input file are dmx, dmy, dmz, PERIOD!
       */

      displacement[p] = efv->ext_fld_ndl_val[p][i];
    }

    anneal_map(dim, X_old, displacement, X_new);

    for (p = 0; p < dim; p++) {
      new_coord[p][i] = X_new[p];
    }
  }

  /* Now we need to make sure that all geometry-related calculations are
   * performed with these updated coordinates.
   */

  /*
   * Temporarily save the old coordinates and assign the new
   * displaced coordinates for the annealed file...
   */

  Coor[0] = &new_coord[0][0]; /*Old Deprecated holder that is stil
                               *used
                               */
  Coor[1] = &new_coord[1][0];

  if (dim > 2) {
    Coor[2] = &new_coord[2][0];
  }

  // adjust base mesh
  for (int i = 0; i < num_nodes; i++) {
    int base_node = exo->ghost_node_to_base[i];
    if (base_node != -1) {
      exo->base_mesh->x_coord[base_node] = Coor[0][i];
      if (dim > 1) {
        exo->base_mesh->y_coord[base_node] = Coor[1][i];
      }
      if (dim > 2) {
        exo->base_mesh->z_coord[base_node] = Coor[2][i];
      }
    }
  }

  for (p = 0; p < dim; p++) {
    safer_free((void **)&(new_coord[p]));
  }
  safer_free((void **)&new_coord);

  return (0);
} /* END of routine anneal_mesh_with_external_field */

/*****************************************************************************/
/*****************************************************************************/

void shift_nodal_values(int var, double shift, double *x, int num_total_nodes) {
  int I, ie, matID = -1;

  for (I = 0; I < num_total_nodes; I++) {
    ie = Index_Solution(I, var, 0, 0, matID, pg->imtrx);
    if (ie != -1) {
      x[ie] += shift;
    }
  }
  return;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* load_export_vars -- save requested solution and post-processing vars */
#ifdef LIBRARY_MODE
int load_export_vars(const int nodes, dbl x[], dbl *x_pp)

{
  double var;
  int i, j, k, ivar;
  int err = 0;
  int scount = 0, pcount = 0;
  /*
    int nsoln[MAX_EXTERNAL_FIELD];
    int npost[MAX_EXTERNAL_FIELD];
  */

  /* Ensure that exports were requested */
  if (Num_Export_XS == 0 && Num_Export_XP == 0) {
    fprintf(stderr, "No export variables have been requested!\n");
    return 0;
  }

  /*
  for (i=0; i<MAX_EXTERNAL_FIELD; i++) {
  nsoln[i] = Export_XS_ID[i];
  npost[i] = Export_XP_ID[i]; }
  */

  /* Ensure good arrays were passed in */
  /*
    if (nsoln == NULL || npost == NULL)
      {
        fprintf(stderr, "Invalid variable ID array: nsoln or npost");
        return -1;
      }
  */
  if (libio->xsoln == NULL || libio->xpost == NULL) {
    fprintf(stderr, "Invalid variable storage array: xsoln or xpost");
    return -1;
  }

  /* Retrieve solution variables in the order requested */
  for (i = 0; i < Num_Export_XS; i++) {
    ivar = Export_XS_ID[i];
    if (Num_Var_In_Type[pg->imtrx][ivar] == 0) {
      fprintf(stderr, "Inactive variable specified: ID = %d\n", ivar);
      return -1;
    }
    for (j = 0; j < nodes; j++) {
      k = Index_Solution(j, ivar, 0, 0, -1, pg->imtrx);
      var = ((k >= 0) ? x[k] : 0.0);
      libio->xsoln[scount] = var;
      scount++;
    }
  }

  /* Retrieve post-processing variables in the order requested */
  if (Num_Export_XP > 0 && x_pp != NULL) {
    for (i = 0; i < Num_Export_XP; i++) {
      for (j = 0; j < nodes; j++) {
        libio->xpost[pcount] = x_pp[pcount];
        pcount++;
      }
    }
  }

  return 0;
}
#endif

/*****************************************************************************/

int variable_stats(double *x, const double time, const int coord_linear) {
  int i, var, idv, mn;
  double max[MAX_VARIABLE_TYPES], min[MAX_VARIABLE_TYPES];
  double sqr[MAX_VARIABLE_TYPES], mean[MAX_VARIABLE_TYPES];
  int ncp[MAX_VARIABLE_TYPES]; /* number of ea var contributing to norm */
  int var_somewhere;           /* boolean */
  double crd_sum[DIM], crd_sqr[2 * DIM], crd_cross[DIM][MAX_VARIABLE_TYPES];
  double fcn_sum[MAX_VARIABLE_TYPES], fcn_sqr[MAX_VARIABLE_TYPES];
  int crd_ncp[DIM], fcn_ncp[MAX_VARIABLE_TYPES];
  int pdim = pd->Num_Dim;
#ifdef PARALLEL
  double max_buf[MAX_VARIABLE_TYPES];  /* accumulated over all procs */
  double min_buf[MAX_VARIABLE_TYPES];  /* accumulated over all procs */
  double mean_buf[MAX_VARIABLE_TYPES]; /* accumulated over all procs */
  double sqr_buf[MAX_VARIABLE_TYPES];  /* accumulated over all procs */
  int ncp_buf[MAX_VARIABLE_TYPES];
  double crd_sum_buf[DIM];                       /* accumulated over all procs */
  double crd_sqr_buf[2 * DIM];                   /* accumulated over all procs */
  double crd_cross_buf[DIM][MAX_VARIABLE_TYPES]; /* accumulated over all procs */
  double fcn_sum_buf[MAX_VARIABLE_TYPES];        /* accumulated over all procs */
  double fcn_sqr_buf[MAX_VARIABLE_TYPES];        /* accumulated over all procs */
  int crd_ncp_buf[MAX_VARIABLE_TYPES];
  int fcn_ncp_buf[MAX_VARIABLE_TYPES];
#endif

  memset(ncp, 0, sizeof(int) * MAX_VARIABLE_TYPES);
  init_vec_value(sqr, 0.0, MAX_VARIABLE_TYPES);
  init_vec_value(mean, 0.0, MAX_VARIABLE_TYPES);
  init_vec_value(min, DBL_MAX, MAX_VARIABLE_TYPES);
  init_vec_value(max, -DBL_MAX, MAX_VARIABLE_TYPES);
  if (coord_linear) {
    memset(crd_ncp, 0, sizeof(int) * DIM);
    memset(fcn_ncp, 0, sizeof(int) * MAX_VARIABLE_TYPES);
    memset(crd_sum, 0.0, sizeof(double) * DIM);
    memset(crd_sqr, 0.0, sizeof(double) * 2 * DIM);
    memset(crd_cross, 0.0, sizeof(double) * DIM * MAX_VARIABLE_TYPES);
    memset(fcn_sum, 0.0, sizeof(double) * MAX_VARIABLE_TYPES);
    memset(fcn_sqr, 0.0, sizeof(double) * MAX_VARIABLE_TYPES);
  }

  if (TimeIntegration == STEADY) {
    DPRINTF(stdout, "\nVariable Stats:  parm=%g \n", time);
  } else {
    DPRINTF(stdout, "\nVariable Stats:  time=%g \n", time);
  }
  DPRINTF(stdout, "\t   min        max        mean       stdev     #values\n");
  DPRINTF(stdout, "=============================================================\n");
  for (var = 0; var < MAX_VARIABLE_TYPES; var++) {
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
          if (x[idv] > max[var])
            max[var] = x[idv];
          if (x[idv] < min[var])
            min[var] = x[idv];
          mean[var] += x[idv];
          sqr[var] += SQUARE(x[idv]);
          ncp[var]++;
        }
      }
    }
  }
  if (coord_linear) {
    int dir, dirj, dvar, idx;
    double coord, coordj;
    /*  Coordinate stuff first  */
    for (i = 0; i < DPI_ptr->num_owned_nodes; i++) {
      for (dir = 0; dir < pdim; dir++) {
        var = MESH_DISPLACEMENT1 + dir;
        if (upd->vp[pg->imtrx][var] > -1) {
          idv = Index_Solution(i, var, 0, 0, -1, pg->imtrx);
          coord = Coor[dir][i] + x[idv];
        } else {
          coord = Coor[dir][i];
        }
        crd_sum[dir] += coord;
        crd_sqr[dir] += SQUARE(coord);
        if (dir == (var - MESH_DISPLACEMENT1)) {
          crd_ncp[dir]++;
        }
        for (dirj = dir + 1; dirj < pdim; dirj++) {
          var = MESH_DISPLACEMENT1 + dirj;
          if (upd->vp[pg->imtrx][var] > -1) {
            idv = Index_Solution(i, var, 0, 0, -1, pg->imtrx);
            coordj = Coor[dirj][i] + x[idv];
          } else {
            coordj = Coor[dirj][i];
          }
          crd_sqr[pdim + dir + dirj - 1] += coord * coordj;
        }
      }
    }

    for (i = 0; i < DPI_ptr->num_owned_nodes; i++) {
      for (dir = 0; dir < pdim; dir++) {
        dvar = MESH_DISPLACEMENT1 + dir;
        if (upd->vp[pg->imtrx][dvar] > -1) {
          idx = Index_Solution(i, dvar, 0, 0, -1, pg->imtrx);
          coord = Coor[dir][i] + x[idx];
          for (var = 0; var < MAX_VARIABLE_TYPES; var++) {
            if (upd->vp[pg->imtrx][var] > -1) {
              var_somewhere = FALSE;
              for (mn = 0; mn < upd->Num_Mat; mn++) {
                idv = Index_Solution(i, var, 0, 0, mn, pg->imtrx);
                if (idv != -1) {
                  var_somewhere = TRUE;
                  break;
                }
              }
              if (var_somewhere) {
                crd_cross[dir][var] += coord * x[idv];
              }
            }
          }
        }
      }

      for (var = 0; var < MAX_VARIABLE_TYPES; var++) {
        if (upd->vp[pg->imtrx][var] > -1) {
          var_somewhere = FALSE;
          for (mn = 0; mn < upd->Num_Mat; mn++) {
            idv = Index_Solution(i, var, 0, 0, mn, pg->imtrx);
            if (idv != -1) {
              var_somewhere = TRUE;
              break;
            }
          }
          if (var_somewhere) {
            fcn_sum[var] += x[idv];
            fcn_sqr[var] += SQUARE(x[idv]);
            fcn_ncp[var]++;
          }
        }
      }
    }
  }
/*
 * Find global quantities
 */
#ifdef PARALLEL
  MPI_Reduce((void *)max, (void *)max_buf, MAX_VARIABLE_TYPES, MPI_DOUBLE, MPI_MAX, 0,
             MPI_COMM_WORLD);
  MPI_Reduce((void *)min, (void *)min_buf, MAX_VARIABLE_TYPES, MPI_DOUBLE, MPI_MIN, 0,
             MPI_COMM_WORLD);
  MPI_Reduce((void *)mean, (void *)mean_buf, MAX_VARIABLE_TYPES, MPI_DOUBLE, MPI_SUM, 0,
             MPI_COMM_WORLD);
  MPI_Reduce((void *)sqr, (void *)sqr_buf, MAX_VARIABLE_TYPES, MPI_DOUBLE, MPI_SUM, 0,
             MPI_COMM_WORLD);
  MPI_Reduce((void *)ncp, (void *)ncp_buf, MAX_VARIABLE_TYPES, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (coord_linear) {
    int j;
    MPI_Reduce((void *)crd_sum, (void *)crd_sum_buf, DIM, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce((void *)crd_sqr, (void *)crd_sqr_buf, 2 * DIM, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce((void *)crd_cross, (void *)crd_cross_buf, DIM * MAX_VARIABLE_TYPES, MPI_DOUBLE,
               MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce((void *)fcn_sum, (void *)fcn_sum_buf, MAX_VARIABLE_TYPES, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce((void *)fcn_sqr, (void *)fcn_sqr_buf, MAX_VARIABLE_TYPES, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce((void *)crd_ncp, (void *)crd_ncp_buf, DIM, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce((void *)fcn_ncp, (void *)fcn_ncp_buf, MAX_VARIABLE_TYPES, MPI_INT, MPI_SUM, 0,
               MPI_COMM_WORLD);
    for (i = 0; i < pdim; i++) {
      crd_sum[i] = crd_sum_buf[i];
      crd_sqr[i] = crd_sqr_buf[i];
      crd_sqr[pdim + i] = crd_sqr_buf[pdim + i];
      crd_ncp[i] = crd_ncp_buf[i];
      for (j = 0; j < MAX_VARIABLE_TYPES; j++) {
        crd_cross[i][j] = crd_cross_buf[i][j];
      }
    }
    for (i = 0; i < MAX_VARIABLE_TYPES; i++) {
      fcn_sum[i] = fcn_sum_buf[i];
      fcn_sqr[i] = fcn_sqr_buf[i];
      fcn_ncp[i] = fcn_ncp_buf[i];
    }
  }

  for (i = 0; i < MAX_VARIABLE_TYPES; i++) {
    max[i] = max_buf[i];
    min[i] = min_buf[i];
    mean[i] = mean_buf[i];
    sqr[i] = sqr_buf[i];
    ncp[i] = ncp_buf[i];
  }
#endif
  if (ProcID == 0) {
    for (var = 0; var < MAX_VARIABLE_TYPES; var++) {
      double std_arg = 0.;
      if (upd->vp[pg->imtrx][var] > -1) {
        mean[var] /= ncp[var];
        std_arg = sqr[var] - ncp[var] * SQUARE(mean[var]);
        if (std_arg > 0.) {
          std_arg = sqrt(std_arg / (double)ncp[var]);
        }
        DPRINTF(stdout, "%s \t%8g  %8g  %8g  %8g  %8d\n", Var_Name[var].name2, min[var], max[var],
                mean[var], std_arg, ncp[var]);
      }
    }
  }
  DPRINTF(stdout, "============================================================= \n");
  if (ProcID == 0 && coord_linear) {
    int dir, dirj, j, k;
    double reg_A[4][4], reg_RHS[4][MAX_PROB_VAR];
    int n_rhs, reg_N = 0, reg_order = 4;
    double db_reg_N;

    for (dir = 0; dir < pdim; dir++) {
      reg_N = MAX(reg_N, crd_ncp[dir]);
    }
    db_reg_N = (double)reg_N;
    n_rhs = 0; /* Since all processors collected, try some reduced arrays on Proc 1*/
    for (var = 0; var < MAX_VARIABLE_TYPES; var++) {
      if (fcn_ncp[var] == reg_N) {
        reg_RHS[0][n_rhs] = fcn_sum[var];
        for (dir = 0; dir < pdim; dir++) {
          reg_RHS[dir + 1][n_rhs] = crd_cross[dir][var];
        }
        n_rhs++;
      }
    }

    memset(reg_A, 0.0, sizeof(double) * 16);
    reg_A[0][0] = db_reg_N;
    for (dir = 0; dir < pdim; dir++) {
      reg_A[dir + 1][dir + 1] = crd_sqr[dir];
      reg_A[0][dir + 1] = crd_sum[dir];
      reg_A[dir + 1][0] = crd_sum[dir];
      for (dirj = dir + 1; dirj < pdim; dirj++) {
        reg_A[dir + 1][dirj + 1] = crd_sqr[pdim + dir + dirj - 1];
        reg_A[dirj + 1][dir + 1] = crd_sqr[pdim + dir + dirj - 1];
      }
    }
    // linear solve for linear coordinate regression
    // convert to column major
    if (0) {
      double reg_W[4], A[16], RHS[4 * MAX_PROB_VAR];
      int reg_pvt[4], LDA = 4, info;
      for (i = 0; i < reg_order; i++) {
        for (j = 0; j < reg_order; j++) {
          A[i * reg_order + j] = reg_A[j][i];
        }
      }
      for (i = 0; i < n_rhs; i++) {
        for (j = 0; j < reg_order; j++) {
          RHS[i * reg_order + j] = reg_RHS[j][i];
        }
      }

      dsysv_("U", &reg_order, &n_rhs, A, &LDA, reg_pvt, RHS, &LDA, reg_W, &LDA, &info, 1);
      if (info > 0)
        fprintf(stderr, "Linear matrix singular at row %d\n", info);
      if (info < 0)
        fprintf(stderr, "Illegal value for Linear matrix %d\n", info);
      // transpose (revert to row major)
      for (i = 0; i < reg_order; i++) {
        for (j = 0; j < reg_order; j++) {
          reg_A[i][j] = A[j * reg_order + i];
        }
      }
      for (i = 0; i < reg_order; i++) {
        for (j = 0; j < n_rhs; j++) {
          reg_RHS[i][j] = RHS[j * n_rhs + i];
        }
      }
    } else {
      double pivot, sum;
      /* Gauss Elimination - no pivoting */
      for (i = 0; i < reg_order - 1; i++) {
        for (j = i + 1; j < reg_order; j++) {
          pivot = reg_A[j][i] / reg_A[i][i];
          for (k = i + 1; k < reg_order; k++) {
            reg_A[j][k] -= pivot * reg_A[i][k];
          }
          for (k = 0; k < n_rhs; k++) {
            reg_RHS[j][k] -= pivot * reg_RHS[i][k];
          }
        }
      }
      for (k = 0; k < n_rhs; k++) {
        reg_RHS[reg_order - 1][k] /= reg_A[reg_order - 1][reg_order - 1];
        for (i = reg_order - 2; i >= 0; i--) {
          sum = 0.;
          for (j = i + 1; j < reg_order; j++) {
            sum += reg_A[i][j] * reg_RHS[j][k];
          }
          reg_RHS[i][k] = (reg_RHS[i][k] - sum) / reg_A[i][i];
        }
      }
    }

    if (TimeIntegration == STEADY) {
      DPRINTF(stdout, "\nLinear Regression:  parm=%g \n", time);
    } else {
      DPRINTF(stdout, "\nLinear Regression:  time=%g \n", time);
    }
    DPRINTF(stdout, "var  \tFcn = A + B*x + C*y + D*z  \tcorrelation \trmse\n");
    DPRINTF(stdout, "=============================================================\n");
    /* construct linear matrix solution */

    n_rhs = 0;
    for (var = 0; var < MAX_VARIABLE_TYPES; var++) {
      if (upd->vp[pg->imtrx][var] > -1) {
        double sigma_fcn, sum_lin, sigma_lin, sum_lin_sqr, cross_sum, rmse, corr;
        /* correlation coefficient calculation  */
        if (fcn_ncp[var] == reg_N) {
          sigma_fcn = db_reg_N * fcn_sqr[var] - SQUARE(fcn_sum[var]);
          sum_lin = db_reg_N * reg_RHS[0][n_rhs];
          for (j = 0; j < pdim; j++) {
            sum_lin += reg_RHS[j + 1][n_rhs] * crd_sum[j];
          }
          sum_lin_sqr = reg_RHS[0][n_rhs] * (2. * sum_lin - reg_RHS[0][n_rhs] * db_reg_N);
          for (j = 0; j < pdim; j++) {
            sum_lin_sqr += SQUARE(reg_RHS[j + 1][n_rhs]) * crd_sqr[j];
            for (k = j + 1; k < pdim; k++) {
              sum_lin_sqr +=
                  2. * reg_RHS[j + 1][n_rhs] * reg_RHS[k + 1][n_rhs] * crd_sqr[pdim + j + k - 1];
            }
          }
          sigma_lin = db_reg_N * sum_lin_sqr - SQUARE(sum_lin);
          cross_sum = reg_RHS[0][n_rhs] * fcn_sum[var];
          for (j = 0; j < pdim; j++) {
            cross_sum += reg_RHS[j + 1][n_rhs] * crd_cross[j][var];
          }
          corr = (db_reg_N * cross_sum - fcn_sum[var] * sum_lin) / sqrt(sigma_fcn * sigma_lin);
          rmse = sqrt((fcn_sqr[var] - 2. * cross_sum + sum_lin_sqr) / db_reg_N);
          DPRINTF(stdout, "%s \t%7g  %7g  %7g  %7g  \t%g  \t%g\n", Var_Name[var].name2,
                  reg_RHS[0][n_rhs], reg_RHS[1][n_rhs], reg_RHS[2][n_rhs], reg_RHS[3][n_rhs], corr,
                  rmse);
          n_rhs++;
        }
      }
    }
    DPRINTF(stdout, "=============================================================\n");
  }
  return (1);
}
/*****************************************************************************/

/*****************************************************************************/
/*  end of file rf_solve.c  */
/*****************************************************************************/
