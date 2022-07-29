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

/* NEWTON ITERATIVE SOLUTION DRIVER WITH SENSITIVITIES WRT SPECIFIED
 * FLUX AND/OR DATA PARAMETER AND/OR CONTINUATION PARAMETER
 */

#include "sl_epetra_interface.h"
#include "sl_util_structs.h"

#define GOMA_MM_SOL_NONLINEAR_C
/* Needed to declare POSIX function drand48 */
#define _XOPEN_SOURCE
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ac_stability_util.h"
#include "bc/rotate.h"
#include "bc/rotate_coordinates.h"
#include "dp_comm.h"
#include "dp_types.h"
#include "dp_utils.h"
#include "dpi.h"
#include "exo_struct.h"
#include "loca_const.h"
#include "md_timer.h"
#include "mm_as.h"
#include "mm_as_structs.h"
#include "mm_augc_util.h"
#include "mm_fill.h"
#include "mm_fill_aux.h"
#include "mm_fill_ls.h"
#include "mm_fill_solid.h"
#include "mm_fill_util.h"
#include "mm_input.h"
#include "mm_more_utils.h"
#include "mm_mp.h"
#include "mm_mp_structs.h"
#include "mm_numjac.h"
#include "mm_post_def.h"
#include "mm_post_proc.h"
#include "mm_sol_nonlinear.h"
#include "mm_unknown_map.h"
#include "rf_allo.h"
#include "rf_bc.h"
#include "rf_bc_const.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_io.h"
#include "rf_io_const.h"
#include "rf_io_structs.h"
#include "rf_mp.h"
#include "rf_node_const.h"
#include "rf_solver.h"
#include "rf_solver_const.h"
#include "rf_util.h"
#include "sl_amesos_interface.h"
#include "sl_auxutil.h"
#include "sl_aztecoo_interface.h"
#include "sl_lu.h"
#include "sl_matrix_util.h"
#include "sl_petsc.h"
#include "sl_stratimikos_interface.h"
#include "sl_umf.h"
#include "sl_util.h"
#include "std.h"
#include "wr_exo.h"
#include "wr_side_data.h"

/*
 * EDW: The prototype for function "mf_sol_lineqn" has been moved
 * to mm_sol_nonlinear.h so that LOCA has access to it.
 */

/*
 * This variable is used for more than one linear solver package.
 */

static int first_linear_solver_call = TRUE;

/*
 * Default: do not attempt to use Harwell MA28 linear solver. Kundert's is
 *          more robust and Harwell has a better successor to MA28 that you
 *          can buy with money.
 *
 * #define HARWELL
 */

#ifdef PARALLEL
#ifndef MPI
#define MPI /* otherwise az_aztec.h trounces MPI_Request */
#endif
#endif

#include "az_aztec.h"

static int soln_sens                /* mm_sol_nonlinear.c                        */
    (double,                        /* lambda - parameter                        */
     double[],                      /* x - soln vector                           */
     double[],                      /* xdot - dxdt predicted for new time step   */
     double,                        /* delta_s - step                            */
     Exo_DB *,                      /* exo                                       */
     Dpi *,                         /* dpi                                       */
     Comm_Ex *,                     /* cx                                        */
     double[],                      /* res_p                                     */
     int,                           /* numProcUnknowns                           */
     int[],                         /* ija                                       */
     double[],                      /* a                                         */
     double[],                      /* x_old                                     */
     double[],                      /* x_older                                   */
     double[],                      /* xdot_old                                  */
     double[],                      /* x_update                                  */
     double,                        /* delta_t                                   */
     double,                        /* theta                                     */
     double,                        /* time_value                                */
     int,                           /* num_total_nodes                           */
     double[],                      /* res_m                                     */
     int,                           /* Factor_Flag                               */
     int,                           /* matr_form                                 */
     double[],                      /* resid_vector_sens                         */
     double[],                      /* x_sens                                    */
     double **,                     /* x_sens_p                                  */
     double[],                      /* scale                                     */
     struct GomaLinearSolverData *, /* ams - ptrs to Aztec linear    *
                                     * systems                       */
     int,                           /* first_linear_solver_call                  */
     int,                           /* Norm_below_tolerance                      */
     int,                           /* Rate_above_tolerance                      */
     int,                           /* vector_id                                 */
     int,                           /* sens_type                                 */
     int,                           /* sens_id                                   */
     int,                           /* sens_flt                                  */
     int,                           /* sens_flt2                                  */
     int *,                         /* frontal solver variables */
     int *,                         /*  fs->ncod  */
     double *,                      /*  fs->bc    */
     double *,                      /*  smallpiv  */
     double *,                      /*  singpiv   */
     int *,                         /*  iautopiv  */
     int *,                         /*  iscale    */
     double *,                      /*  scaling_max */
     double *,                      /*  h_elem_avg  */
     double *,
     int,     /* UMF_system_id */
     char[]); /* calling purpose */

/*
 * The one place place these global variables are defined.
 * Declared frequently via std.h, etc. every where else.
 */

int neg_elem_volume = FALSE;
int neg_elem_volume_global = FALSE;

int neg_lub_height = FALSE;
int neg_lub_height_global = FALSE;

int zero_detJ = FALSE;
int zero_detJ_global = FALSE;

/*

   GOMA NON-LINEAR EQUATION SOLVER

   FOR EACH LINEAR SUB-SOLVE,

   STANDARD SOLVE:

             [ J ] [ x ]  =  [ f ]

   IF HAVE AUGMENTING CONDITIONS, SOLVE:

             [ J  B ] [ x ]    [ f ]
    [M][v] = [      ] [   ] =  [   ]
             [ C  D ] [ y ]    [ g ]

   J: JACOBIAN MATRIX, dR/dq
   q: VECTOR OF GOMA UNKNOWNS
   R: RESIDUAL VECTOR = -f
   B: dR/dp
   p: VECTOR OF UNKNOWNS ASSOCIATED WITH AUGMENTING CONDITIONS
   C: dA/dq
   A: JACOBIAN MATRIX OF AUGMENTING CONDITIONS
   D: dA/dp
   f: f = -R
   g: RESIDUAL VECTOR OF AUGMENTING CONDITIONS

   AUGMENTING CONDITIONS BY IAN GATES
   2/98 - 10/98

    IF CONTINUATION, DO:

    J dq/dp = - dR/dp

    IN soln_sens AFTER CONVERGED

    J: JACOBIAN MATRIX, dR/dq
    q: VECTOR OF GOMA UNKNOWNS
    R: RESIDUAL VECTOR
    p: PARAMETER

    SENSITIVITIES AND AUGMENTING CONDITIONS
    IDG 2/98 - 10/98

    FLUX AND DATA PARAMETER SENSITIVITIES
    RBS 4/99 - 6/99

    MODIFIED FOR HUNTING SENSITIVITIES
    IDG 7/99 - 7/99

*/

int solve_nonlinear_problem(struct GomaLinearSolverData *ams,
                            /* ptrs to Aztec linear systems */
                            double x[],        /* soln vector on this proc */
                            double delta_t,    /* time step size */
                            double theta,      /* parameter to vary time
                                                * integration from
                                                *   explicit (theta = 1) to
                                                *   implicit (theta = 0) */
                            double x_old[],    /* soln vector @ previous time */
                            double x_older[],  /* soln vector @ previous, previous time */
                            double xdot[],     /* dxdt predicted for new time */
                            double xdot_old[], /* dxdt for previous time */
                            double resid_vector[],
                            double x_update[],
                            double scale[],        /*Scale factor held for modified newton
                                                    *resolves */
                            int *converged,        /* whether the Newton iteration
                                                    * has converged (out) */
                            int *nprint,           /* counter for time step number */
                            int tev,               /* total number elem variables to
                                                    * output to EXODUS II file */
                            int tev_post,          /* extra element post processing
                                                    * results */
                            double *glob_var_vals, /* global variable values */
                            RESULTS_DESCRIPTION_STRUCT *rd,
                            /* details about post proc vars */
                            int *gindex,
                            int *gsize,
                            double *gvec,
                            double ***gvec_elem,
                            double time_value,
                            Exo_DB *exo,
                            Dpi *dpi,
                            Comm_Ex *cx,
                            int nt,
                            int *time_step_reform,
                            int is_steady_state,
                            double x_AC[], /* updating evp quantities */
                            double x_AC_dot[],
                            double lambda,
                            double *resid_vector_sens,
                            double *x_sens,
                            double **x_sens_p, /*  solution sensitivities */
                            void *con_ptr) {

  static int prev_matrix = 0;

  double *a = ams->val;  /* nonzero values of a CMSR matrix */
  int *ija = ams->bindx; /* column pointer array into matrix "a"*/

  int *ija_save = ams->belfry; /* to hide external row/eqns from Aztec */

  static int dofs_hidden = FALSE; /* boolean indicating when dofs are hidden */

  int i, j, k;

  int numProcUnknowns; /* number of degrees of freedom this processor
                        * sees (internal+boundary+external) but no
                        * more. */

  int GNumUnknowns; /* Global number of unknowns in the */
                    /* system    */
  int inewton;      /* Newton iteration counter */

  int return_value; /* nonzero if things screw up ...  */

  int print_damp_factor; /* Always printed after update norms */

  int print_visc_sens; /* Always printed after update norms */

  double Norm[5][3]; /* Global norms... */
                     /*   [0][0] == residual, L_oo norm */
                     /*   [0][1] == residual, L_1  norm */
                     /*   [0][2] == residual, L_2  norm */
                     /*   [1][0] == correction, L_oo norm */
                     /*   [1][1] == correction, L_1  norm */
                     /*   [1][2] == correction, L_2  norm */
                     /*   [2][0] == AC residual, L_oo norm */
                     /*   [2][1] == AC residual, L_1  norm */
                     /*   [2][2] == AC residual, L_2  norm */
                     /*   [3][0] == AC correction, L_oo norm */
                     /*   [3][1] == AC correction, L_1  norm */
                     /*   [3][2] == AC correction, L_2  norm */
                     /*   [4][0] == solution, L_oo norm */
                     /*   [4][1] == solution, L_1  norm */
                     /*   [4][2] == solution, L_2  norm */

  double Norm_r[2][3];                     /* Relative (to solution vector) norms... */
                                           /*   [0][0] == correction, L_oo norm */
                                           /*   [0][1] == correction, L_1  norm */
                                           /*   [0][2] == correction, L_2  norm */
                                           /*   [1][0] == AC correction, L_oo norm */
                                           /*   [1][1] == AC correction, L_1  norm */
                                           /*   [1][2] == AC correction, L_2  norm */
  double Resid_Norm_stack[3];              /* Place holder for last residual norms   */
  double Soln_Norm_stack[3];               /* Place holder for last update norms   */
  double Conv_order = 0, Soln_order = 0;   /* Order of convergence  */
  double Conv_rate = 0, Soln_rate = 0;     /* Convergence rates, i.e. neg. semilog slope*/
  double AConv_order = 0, ASoln_order = 0; /* AC counterparts  */
  double AConv_rate = 0, ASoln_rate = 0;
  double AC_Resid_Norm_stack[3]; /* Place holder for last residual norms   */
  double AC_Soln_Norm_stack[3];  /* Place holder for last update norms   */
  int Norm_below_tolerance;      /* Boolean for modified newton test*/
  int Rate_above_tolerance;      /* Boolean for modified newton test*/
  int step_reform;               /* counter for Jacobian reformation */

  double Reltol = 1.0e-2, Abstol = 1.0e-6; /* LOCA convergence criteria */
  int continuation_converged = TRUE;
  int num_total_nodes = dpi->num_universe_nodes;
  /* Number of nodes that each processor is
   * responsible for                           */
  dbl h_elem_avg; /* global average element size for PSPG */
  dbl U_norm;     /* global average velocity for PSPG */

  double delta_s = 0.0; /* STEP */
  int *ncod = NULL;
  dbl *bc = NULL; /* Dummy pointers for non-Front cases */

  static char yo[] = "solve_nonlinear_problem"; /* routine identifier */

  /*
   * The following vectors all have lengths equal to the number of local
   * unknowns - i.e., equations for which the current processor is responsible.
   */

  double *delta_x; /* update */

  int error, why;
  int num_unk_r, num_unk_x;

  char dofname_r[80];
  char dofname_nr[80];
  char dofname_x[80];
  char dofname_g[80];
  char dofname_ry[80];
  char dofname_y[80];

  static int UMF_system_id; /* Used to uniquely identify the
                             * "regular" system to solve from the
                             * other UMF systems. */
  static int Factor_Flag;   /* Legend::
                     * --------
                     * 1 -- this is the first LU,
                     * 2 -- >1st LU, similar matrix structure,
                     * 3 -- new rhs, use already factored LU
  
                     UMFPACK USAGE:
  
                     0   : LOAD MATRIX AND ANALYSIS/DECOMPOSITION
                           AND BACK SUBSTITUTION
                     1   : LOAD MATRIX AND DECOMPOSITION USING PAST ANALYSIS
                           AND BACK SUBSTITUTION
                     > 2 : BACK SUBSTITUTION ONLY
  
                     */

  int matr_form = 0; /* 1: MSR FORMAT MATRIX FOR UMFPACK DRIVER */

  int i_post;
  int err = 0;

  dbl total_mesh_volume; /* global variable value */

  int local_order;      /* of the unknowns that this processor owns */
  int local_order_plus; /* and if the external rows are included */
  int local_nnz;        /* number of nonzero matrix entries that
                         * are owned by this processor */
  int local_nnz_plus;   /* and if the external rows are included */

  int global_order;      /* order of the global system */
  int global_order_plus; /* and if the external rows are overincluded */
  int global_nnz;        /* a sum of the number of nonzero matrix
                          * entries owned by each processor */
  int global_nnz_plus;   /* a sum that overincludes the external rows */

  char stringer[80];    /* holding format of num linear solve itns */
  char stringer_AC[80]; /* holding format of num AC linear solve itns */

  dbl a_start; /* mark start of assembly */
  dbl a_end;   /* mark end of assembly */

  char ctod[80]; /* hold current time of day */

  dbl s_start; /* mark start of solve */
  dbl s_end;   /* mark end of solve */

  /*
   * With the advent of the Benner's frontal solver it is more difficult to
   * profile the CPU time used by the assembly and by the solver since they
   * are more closely intertwined by the frontal solver calling element level
   * assembly.
   *
   * As some sop toward ameliorating this situation, count up the total
   * CPU time for both assembly and solving here and then access accumulated
   * totals for just assembly that get suitably initialized and incremented
   * in matrix_fill().
   */

  double asmslv_time;
  double slv_time;

  int iAC, jAC;
  int num_unk_g, num_unk_y;
  dbl *res_p, *res_m;
  dbl *gAC, *hAC, *tAC, *yAC;
  dbl **bAC = NULL, **cAC = NULL, **dAC = NULL, **sAC = NULL, **wAC = NULL;
  dbl *tempAC;

  int iiAC, imaxAC = 0;
  int *ordAC = NULL;
  dbl bigAC, dumAC, sumAC;

  /* frontal solver parameters */
  int mf_resolve;
  dbl smallpiv;
  dbl singpiv;
  int iautopiv;
  int iscale; /* you will have to turn this off for resolves */
  dbl scaling_max;

  dbl ac_start = 0; /* mark start of AC assembly */
  dbl ac_end = 0;   /* mark end of AC assembly */

  dbl sc_start = 0; /* mark start of AC solve */
  dbl sc_end = 0;   /* mark end of AC solve */

  double param_val;
  int sens_vec_ct;
  char sens_caller[40]; /* string containing caller of soln_sens */

  int linear_solver_blk;               /* count calls to AZ_solve() */
  int total_ls_its = 0;                /* linear solver iteration counter */
  int num_linear_solve_blks;           /* one pass for now */
  int matrix_solved;                   /* boolean */
  dbl gamma, beta;                     /* copies of Newmark-beta time integration
                                          parameters */
  int alc_with_acs = FALSE;            /* Flag to handle arc length eqn as an AC   */
  dbl equation_ALC, delta_y_ALC = 0.0; /* Arc length algorithm parameters */
  struct con_struct *con = con_ptr;
  struct Level_Set_Data *ls_old;

  /* Had to add these for newmark-beta updates based on material */
  UMI_LIST_STRUCT *curr_mat_list;
  NODE_INFO_STRUCT *node_ptr;
  int num_mat, imat, mat_index, b;
  bool solve_skipped = false;

  if (pg->imtrx != prev_matrix) {
    first_linear_solver_call = TRUE;
    prev_matrix = pg->imtrx;
  }

  /*
   * Begin executable statements...
   */

  sens_caller[0] = '\0';
  return_value = 0; /* Set the return value to failure
                     * until convergence is achieved */

  /*
   * INITIALIZE
   */

  num_unk_r = num_unk_x = num_unk_g = num_unk_y = 0;

  memset(Norm, 0, 3 * 5 * sizeof(double));
  memset(Norm_r, 0, 3 * 2 * sizeof(double));

  if (Debug_Flag) {
    DPRINTF(stderr, "%s: max Newton itns = %d\n", yo, Max_Newton_Steps);
    DPRINTF(stderr, "%s: convergence tol = %.1e %.1e \n", yo, Epsilon[pg->imtrx][0],
            Epsilon[pg->imtrx][2]);
  }

  /*
   * The number of nonzero matrix entries that this processor owns.
   * External rows are not counted, though ija[] does stretch out to
   * include them too for the parallel case. In the parallel case,
   * watch out to use the bigger ija[] for assembly purposes, but the
   * smaller one for Aztec (after processing via hide_external()).
   */
  numProcUnknowns = NumUnknowns[pg->imtrx] + NumExtUnknowns[pg->imtrx];
  if (strcmp(Matrix_Format, "msr") == 0) {
    local_order = ams->npn;
    local_order_plus = ams->npn_plus;
    local_nnz_plus = ams->nnz_plus;
    local_nnz = ams->nnz;
    canine_chaos(local_order, local_order_plus, local_nnz, local_nnz_plus, &global_order,
                 &global_order_plus, &global_nnz, &global_nnz_plus);

    /*
     * Setup some legacy variables.
     */

    NZeros = local_nnz;
    GNZeros = global_nnz;
    GNumUnknowns = global_order;

    if (Debug_Flag) {
      DPRINTF(stderr, "%s: At matrix %d, NumUnknowns  = %d\n", yo, (pg->imtrx + 1),
              NumUnknowns[pg->imtrx]);
      DPRINTF(stderr, "%s: NZeros       = %d\n", yo, NZeros);
      DPRINTF(stderr, "%s: GNZeros      = %d\n", yo, GNZeros);
      DPRINTF(stderr, "%s: GNumUnknowns = %d\n", yo, GNumUnknowns);
    }

    /*
     *  matrix_stats (a, ija, NumUnknowns[pg->imtrx], &NZeros, &GNZeros, &GNumUnknowns);
     */
  }

  asdv(&delta_x, numProcUnknowns);
  asdv(&res_p, numProcUnknowns);
  asdv(&res_m, numProcUnknowns);

  /*
   * Initialize augmenting condition arrays if needed
   */
  if (nAC > 0) {

    /*
     * Set alc_with_acs flag
     * Exclude arc length equation on first step (consistent with LOCA)
     */
    if (loca_in->Cont_Order == 2) {
      alc_with_acs = TRUE;
      if (con->private_info.step_num == 0)
        nAC--;
    }

    asdv(&gAC, numProcUnknowns);
    asdv(&hAC, numProcUnknowns);
    asdv(&tAC, numProcUnknowns);
    asdv(&yAC, numProcUnknowns);
    asdv(&tempAC, numProcUnknowns);

    bAC = Dmatrix_birth(nAC, numProcUnknowns);
    cAC = Dmatrix_birth(nAC, numProcUnknowns);
    dAC = Dmatrix_birth(nAC, nAC);
    sAC = Dmatrix_birth(nAC, nAC);
    wAC = Dmatrix_birth(nAC, numProcUnknowns);

    ordAC = Ivector_birth(nAC);
  }

  /*
   * Initial conditions at the beginning of the Newton iteration loop...
   * -> Set the convervence flag to false
   */
  *converged = FALSE;
  inewton = 0;
  if (Max_Newton_Steps <= 0)
    *converged = TRUE;

  /*
   * EDW: LOCA, when used, controls UMF_system_id as it calls the solver
   * in other places as well.
   */

  if (con_ptr != NULL)
    UMF_system_id = LOCA_UMF_ID;

  if (Time_Jacobian_Reformation_stride > 1) {
    if (nt == *time_step_reform || nt == 0) {
      /* here we are simply fooling the tolerance checker to
         reform based on iteration number or time step number */
      Norm_below_tolerance = FALSE;
      Rate_above_tolerance = FALSE;
      init_vec_value(scale, 1.0, numProcUnknowns);

      if (nt != 0) {
        *time_step_reform += Time_Jacobian_Reformation_stride;
      }
    } else {
      Norm_below_tolerance = TRUE;
      Rate_above_tolerance = TRUE;
    }
  } else {
    Norm_below_tolerance = FALSE;
    Rate_above_tolerance = FALSE;
  }

  Resid_Norm_stack[2] = Resid_Norm_stack[1] = Resid_Norm_stack[0] = 0.1;
  Soln_Norm_stack[2] = Soln_Norm_stack[1] = Soln_Norm_stack[0] = 0.1;
  AC_Resid_Norm_stack[2] = AC_Resid_Norm_stack[1] = AC_Resid_Norm_stack[0] = 0.1;
  AC_Soln_Norm_stack[2] = AC_Soln_Norm_stack[1] = AC_Soln_Norm_stack[0] = 0.1;
  step_reform = Newt_Jacobian_Reformation_stride;

  /*
   * Make sure every processor is "finished commenting" about the preliminary
   * problem setup - keep the iteration free of clutter caused by buffered
   * stdouts and stderrs from multiple processors...
   */
  fflush(stdout);
  fflush(stderr);
#ifdef PARALLEL
  (void)MPI_Barrier(MPI_COMM_WORLD);
#endif /* PARALLEL */
  fflush(stdout);
  fflush(stderr);

  if (TimeIntegration == STEADY) {

    DPRINTF(stdout, "\n\n");
    if (Solver_Output_Format & 1)
      DPRINTF(stdout, "         ");
    if (Solver_Output_Format & 2)
      DPRINTF(stdout, "    ");
    if ((Solver_Output_Format & 4) && (Solver_Output_Format & 8) && (Solver_Output_Format & 16)) {
      DPRINTF(stdout, "  R e s i d u a l        ");
    } else {
      DPRINTF(stdout, "  Residual       ");
    }
    if ((Solver_Output_Format & 32) && (Solver_Output_Format & 64) &&
        (Solver_Output_Format & 128)) {
      DPRINTF(stdout, " C o r r e c t i o n\n");
    } else {
      DPRINTF(stdout, " Correction\n");
    }
  } else {

    DPRINTF(stdout,
            "\n    N e w t o n  C o n v e r g e n c e  - I m p l i c i t   T i m e   S t e p\n");
  }

  DPRINTF(stdout, "\n");
  if (Solver_Output_Format & 1)
    DPRINTF(stdout, "   ToD   ");
  if (Solver_Output_Format & 2)
    DPRINTF(stdout, "itn ");
  if (Solver_Output_Format & 4)
    DPRINTF(stdout, "  L_oo  ");
  if (Solver_Output_Format & 8)
    DPRINTF(stdout, "  L_1   ");
  if (Solver_Output_Format & 16)
    DPRINTF(stdout, "  L_2   ");
  if (Solver_Output_Format & 32)
    DPRINTF(stdout, "  L_oo  ");
  if (Solver_Output_Format & 64)
    DPRINTF(stdout, "  L_1   ");
  if (Solver_Output_Format & 128)
    DPRINTF(stdout, "  L_2   ");
  if (Solver_Output_Format & 256)
    DPRINTF(stdout, "lis ");
  if (Solver_Output_Format & 512)
    DPRINTF(stdout, " asm/slv (sec)  ");
  if (Solver_Output_Format & 1024)
    DPRINTF(stdout, "ROrd/Rt ");
  if (Solver_Output_Format & 2048)
    DPRINTF(stdout, "UOrd/Rt ");
  DPRINTF(stdout, "\n");
  if (Solver_Output_Format & 1)
    DPRINTF(stdout, "-------- ");
  if (Solver_Output_Format & 2)
    DPRINTF(stdout, "--- ");
  if (Solver_Output_Format & 4)
    DPRINTF(stdout, "------- ");
  if (Solver_Output_Format & 8)
    DPRINTF(stdout, "------- ");
  if (Solver_Output_Format & 16)
    DPRINTF(stdout, "------- ");
  if (Solver_Output_Format & 32)
    DPRINTF(stdout, "------- ");
  if (Solver_Output_Format & 64)
    DPRINTF(stdout, "------- ");
  if (Solver_Output_Format & 128)
    DPRINTF(stdout, "------- ");
  if (Solver_Output_Format & 256)
    DPRINTF(stdout, "--- ");
  if (Solver_Output_Format & 512)
    DPRINTF(stdout, "--------------- ");
  if (Solver_Output_Format & 1024)
    DPRINTF(stdout, "------- ");
  if (Solver_Output_Format & 2048)
    DPRINTF(stdout, "-------");
  DPRINTF(stdout, "\n");

  if (Num_ROT == 0 /*&& inewton == 0*/ && exo->num_dim == 3) {
    setup_rotated_bc_nodes(exo, dpi, BC_Types, Num_BC, x);
  }

  /*********************************************************************************
   *
   *                         Top of the Newton Iteration Loop
   *
   *********************************************************************************/
  while ((!(*converged)) && (inewton < Max_Newton_Steps)) {
    init_vec_value(resid_vector, 0.0, numProcUnknowns);
    init_vec_value(delta_x, 0.0, numProcUnknowns);
    /* Zero matrix values */
    if (strcmp(Matrix_Format, "epetra") == 0) {
      EpetraPutScalarRowMatrix(ams->RowMatrix, 0.0);
    } else {
      init_vec_value(a, 0.0, ams->nnz);
    }
    get_time(ctod);

    /*
     * Mark the CPU time for combined assembly and solve purposes...
     */

    log_msg("%s: Newton iteration %d", yo, inewton);

#ifndef ALE_DCA_INFO_PLEASE
    if (Solver_Output_Format & 1)
      DPRINTF(stdout, "%s ", ctod);
    if (inewton < 10) {
      if (Solver_Output_Format & 2)
        DPRINTF(stdout, "[%d] ", inewton);
    } else if (inewton < 100) {
      if (Solver_Output_Format & 2)
        DPRINTF(stdout, "%d] ", inewton);
    } else {
      if (Solver_Output_Format & 2)
        DPRINTF(stdout, "%d ", inewton);
    }
#endif

    /*
     * If doing continuation in parallel, restore external matrix rows here.
     */
    if (Num_Proc > 1 && strcmp(Matrix_Format, "msr") == 0 && dofs_hidden) {
      show_external(num_universe_dofs[pg->imtrx],
                    (num_universe_dofs[pg->imtrx] - num_external_dofs[pg->imtrx]), ija, ija_save,
                    a);
      dofs_hidden = FALSE;
    }

    /* For Overlap AC algorithm, reset kinematic AC residuals here */
    if (Do_Overlap && nAC > 0) {
      err = assign_overlap_acs(x, exo);
      if (err == -1) {
        /* error can be caused by neg_elem_volume or from insufficient AC's
         * check for neg_elem_volume and if so, return, otherwise must abort
         */
        check_for_neg_elem_volume(augc[nAC - 1].solid_eb + 1, x, resid_vector, x_old, x_older, xdot,
                                  xdot_old, x_update, &delta_t, &theta, &time_value, exo);
        if (neg_elem_volume) {
          return_value = -1;
          goto free_and_clear;
        } else {
          GOMA_EH(GOMA_ERROR, "Insufficient number of AC's set aside for overlap");
        }
      }
      for (iAC = 0; iAC < nAC; iAC++)
        augc[iAC].lm_resid = 0.0;
    }

    if (pfd != NULL) {
      if (pfd->jac_info != NULL) {
        pfd->Constraint_Integral = 0.0;
        init_vec_value(pfd->jac_info->d_pf_lm, 0.0, numProcUnknowns);
        init_vec_value(pfd->jac_info->d_lm_pf, 0.0, numProcUnknowns);
      }
    }

    /* get global element size and velocity norm if needed for PSPG or Cont_GLS */
    if (upd->matrix_index[VELOCITY1] == pg->imtrx &&
        ((PSPG && Num_Var_In_Type[pg->imtrx][PRESSURE]) ||
         (Cont_GLS && Num_Var_In_Type[pg->imtrx][VELOCITY1]))) {
      h_elem_avg = global_h_elem_siz(x, x_old, xdot, resid_vector, exo, dpi);
      U_norm = global_velocity_norm(x, exo, dpi);
    } else {
      h_elem_avg = 0.;
      U_norm = 0.;
    }

    if (Debug_Flag < 0 && Debug_Flag != -4 &&
        (Debug_Flag != -5 || inewton >= Max_Newton_Steps - 1)) {
      /* Former block 0 in mm_fill.c. Here is some initialization */
      num_total_nodes = dpi->num_universe_nodes;

      /*
       * NORMAL, TANGENT and OTHER Vectors required for ROTATION are calculated
       * ahead of time so we don't run into anomolous behavior due to neclaced
       * elements, junction points, . . .
       */

      if (Num_ROT > 0)
        calculate_all_rotation_vectors(exo, x);
      else if (Use_2D_Rotation_Vectors == TRUE)
        calculate_2D_rotation_vectors(exo, x);

      numerical_jacobian(ams, x, resid_vector, delta_t, theta, x_old, x_older, xdot, xdot_old,
                         x_update, num_total_nodes, First_Elem_Side_BC_Array[pg->imtrx], Debug_Flag,
                         time_value, exo, dpi, &h_elem_avg, &U_norm);

      DPRINTF(stdout, "%s: numerical Jacobian done....\n", yo);
      P0PRINTF("\n-done\n\n");
      exit(0);
    } else {

      if (!Norm_below_tolerance || !Rate_above_tolerance) {
        init_vec_value(resid_vector, 0.0, numProcUnknowns);
        init_vec_value(a, 0.0, (NZeros + 1));
        af->Assemble_Residual = TRUE;
        af->Assemble_Jacobian = TRUE;
        af->Assemble_LSA_Jacobian_Matrix = FALSE;
        af->Assemble_LSA_Mass_Matrix = FALSE;
      } else {
        init_vec_value(resid_vector, 0.0, numProcUnknowns);
        af->Assemble_Residual = TRUE;
        af->Assemble_Jacobian = FALSE;
        af->Assemble_LSA_Jacobian_Matrix = FALSE;
        af->Assemble_LSA_Mass_Matrix = FALSE;
      }
      a_start = ut();
      a_end = a_start;

      /* Former block 0 in mm_fill.c. Here is some initialization */
      num_total_nodes = dpi->num_universe_nodes;

      /* If using volume change metric for element quality, reset here */
      if (nEQM > 0 && eqm->do_vol && af->Assemble_Residual) {
        eqm->vol_sum = 0.0;
        eqm->vol_low = 9999.9;
        eqm->vol_count = 0;
      }

      /*
       * NORMAL, TANGENT and OTHER Vectors required for ROTATION are calculated
       * ahead of time so we don't run into anomolous behavior due to neclaced
       * elements, junction points, . . .
       */
      if (Num_ROT > 0)
        calculate_all_rotation_vectors(exo, x);
      else if (Use_2D_Rotation_Vectors == TRUE)
        calculate_2D_rotation_vectors(exo, x);

      /* Initialize volume constraint, not the most efficient way. */

      numProcUnknowns = NumUnknowns[pg->imtrx] + NumExtUnknowns[pg->imtrx];
      if (nAC > 0) {

        for (iAC = 0; iAC < nAC; iAC++) {
          init_vec_value(augc[iAC].d_evol_dx, 0.0, numProcUnknowns);
          init_vec_value(augc[iAC].d_lsvel_dx, 0.0, numProcUnknowns);
          init_vec_value(augc[iAC].d_lsvol_dx, 0.0, numProcUnknowns);
          augc[iAC].evol = 0.;
          augc[iAC].lsvel = 0.;
          augc[iAC].lsvol = 0.;
        }
      }

      /* Exchange dof before matrix fill so parallel information
         is properly communicated */
      exchange_dof(cx, dpi, x, pg->imtrx);

      err = matrix_fill_full(ams, x, resid_vector, x_old, x_older, xdot, xdot_old, x_update,
                             &delta_t, &theta, First_Elem_Side_BC_Array[pg->imtrx], &time_value,
                             exo, dpi, &num_total_nodes, &h_elem_avg, &U_norm, NULL);

      a_end = ut();
      if (err == -1) {
        return_value = -1;
        goto free_and_clear;
      }

      bool enable_numerical_jacobian = false;
      for (int mn = 0; mn < upd->Num_Mat; mn++) {
        if (((vn_glob[mn]->evssModel == LOG_CONF || vn_glob[mn]->evssModel == LOG_CONF_GRADV) &&
             pd_glob[mn]->v[pg->imtrx][POLYMER_STRESS11]) ||
            ((pd_glob[mn]->v[pg->imtrx][EM_E1_REAL] && pd_glob[mn]->v[pg->imtrx][EM_H1_REAL]))) {
          enable_numerical_jacobian = true;
          break;
        }
      }

      if (enable_numerical_jacobian && af->Assemble_Jacobian == TRUE) {
        err = numerical_jacobian_compute_stress(ams, x, resid_vector, delta_t, theta, x_old,
                                                x_older, xdot, xdot_old, x_update, num_total_nodes,
                                                First_Elem_Side_BC_Array[pg->imtrx], Debug_Flag,
                                                time_value, exo, dpi, &h_elem_avg, &U_norm);
        a_end = ut();
      }

      if (err == -1) {
        return_value = -1;
        goto free_and_clear;
      }

      /* Scale matrix first to get rid of problems with
       * penalty parameter. In front option this is done
       * within the solver
       */
      if (!Norm_below_tolerance || !Rate_above_tolerance) {
        row_sum_scaling_scale(ams, resid_vector, scale);
      } else {
        vector_scaling(NumUnknowns[pg->imtrx], resid_vector, scale);
      }

      if (err == -1) {
        return_value = -1;
        goto free_and_clear;
      }
    }

    /*
     *
     *   DEAL WITH AUGMENTING CONDITION MATRICES
     *
     */

    if (nAC > 0) {
      MF_Args mf_args;

      /* OK, copy all the data into the mf_args structure.  A pointer
         to this guy is sent down to aug cond routines.  Cleans things
         up a bit and makes it easier to call matrix_fill multiple times. */

      mf_args.ams = ams;
      mf_args.x = x;
      mf_args.resid = resid_vector;
      mf_args.x_old = x_old;
      mf_args.x_older = x_older;
      mf_args.xdot = xdot;
      mf_args.xdot_old = xdot_old;
      mf_args.x_update = x_update;
      mf_args.delta_t = &delta_t;
      mf_args.theta_ = &theta;
      mf_args.first_elem_side_bc = First_Elem_Side_BC_Array[pg->imtrx];
      mf_args.time = &time_value;
      mf_args.exo = exo;
      mf_args.dpi = dpi;
      mf_args.num_total_nodes = &num_total_nodes;
      mf_args.h_elem_avg = &h_elem_avg;
      mf_args.U_norm = &U_norm;
      mf_args.estifm = NULL;

      ac_start = ut();

      for (iAC = 0; iAC < nAC; iAC++) {
        memset(cAC[iAC], 0, sizeof(double) * numProcUnknowns);
        memset(bAC[iAC], 0, sizeof(double) * numProcUnknowns);
      }
      memset(gAC, 0, sizeof(double) * numProcUnknowns);

      for (iAC = 0; iAC < nAC; iAC++) {
        switch (augc[iAC].Type) {
        case AC_USERBC:
        case AC_USERMAT:
          user_aug_cond(iAC, nAC, x_AC, bAC, cAC, dAC, gAC, numProcUnknowns, x_sens_p, cx,
                        &mf_args);
          break;

        case AC_FLUX:
        case AC_FLUX_MAT:
        case AC_VOLUME:
        case AC_LS_VEL:
        case AC_POSITION:
        case AC_ANGLE:
          std_aug_cond(iAC, nAC, x_AC, bAC, cAC, dAC, gAC, numProcUnknowns, cx, &mf_args);
          break;

        case AC_OVERLAP:
          /* This case handled below. */
          break;

        case AC_LGRM:
        case AC_PF_CONSTRAINT:
          std_lgr_cond(iAC, nAC, x_AC, bAC, cAC, dAC, gAC, resid_vector, scale, numProcUnknowns, cx,
                       &mf_args);
          break;

        case AC_PERIODIC:
          periodic_bc_cond(iAC, nAC, x_AC, bAC, cAC, dAC, gAC, resid_vector, scale, numProcUnknowns,
                           cx, &mf_args);
          break;

        case AC_ARC_LENGTH:
          alc_aug_cond(iAC, nAC, x_AC, bAC, cAC, dAC, gAC, &equation_ALC, numProcUnknowns, cx,
                       con_ptr, &mf_args);
          break;

        default:
          GOMA_EH(GOMA_ERROR, "Unknown augmenting condition type");
          break;
        }
      }

      /*
       * This function is called just once per iteration to fill in all OVERLAP
       * type AC constraints for fluid/solid problems with overlapping grids
       */
      if (augc[nAC - 1].Type == AC_OVERLAP) {
        /* as of 3/23/05 this routine must use pfd[0] struct */
        ls_old = ls;
        if (pfd->ls[0] != NULL) {
          ls = pfd->ls[0];
          overlap_aug_cond(ija, a, x_AC, gAC, bAC, cAC, dAC, cx, &mf_args);
          ls = ls_old;
        }
      }

      for (iAC = 0; iAC < nAC; iAC++) {
        for (j = 0; j < NumUnknowns[pg->imtrx]; j++) {
          bAC[iAC][j] /= scale[j];
        }
      }

      ac_end = ut();

#if 1
      /* row sum scale AC equations */
      row_sum_scaling_scale_AC(cAC, dAC, gAC, nAC);
#endif

      Norm[2][0] = Loo_norm_1p(gAC, nAC, &num_unk_g, dofname_g);
      Norm[2][1] = L1_norm_1p(gAC, nAC);
      Norm[2][2] = L2_norm_1p(gAC, nAC);
    }

    /*
     * We already know the resid vector now, so let's print out its norms
     * now instead of after a big long matrix elimination...
     */

    print_damp_factor = FALSE;
    print_visc_sens = FALSE;
    Norm[0][0] = Loo_norm(resid_vector, NumUnknowns[pg->imtrx], &num_unk_r, dofname_r);
    Norm[0][1] = L1_norm(resid_vector, NumUnknowns[pg->imtrx]);
    Norm[0][2] = L2_norm(resid_vector, NumUnknowns[pg->imtrx]);

    log_msg("%-38s = %23.16e", "residual norm (L_oo)", Norm[0][0]);
    log_msg("%-38s = %23.16e", "residual norm (L_1)", Norm[0][1]);
    log_msg("%-38s = %23.16e", "residual norm (L_2)", Norm[0][2]);
#ifdef ALE_DCA_INFO_PLEASE
    if (Solver_Output_Format & 1)
      DPRINTF(stdout, "%s ", ctod);
    if (inewton < 10) {
      if (Solver_Output_Format & 2)
        DPRINTF(stdout, "[%d] ", inewton);
    } else if (inewton < 100) {
      if (Solver_Output_Format & 2)
        DPRINTF(stdout, "%d] ", inewton);
    } else {
      if (Solver_Output_Format & 2)
        DPRINTF(stdout, "%d ", inewton);
    }
#endif

    if (Solver_Output_Format & 4)
      DPRINTF(stdout, "%7.1e ", Norm[0][0]);
    if (Solver_Output_Format & 8)
      DPRINTF(stdout, "%7.1e ", Norm[0][1]);
    if (Solver_Output_Format & 16)
      DPRINTF(stdout, "%7.1e ", Norm[0][2]);

    fflush(stdout);

    if ((inewton > 0) && (Norm[0][2] < Epsilon[pg->imtrx][0]) &&
        (Norm[0][0] < Epsilon[pg->imtrx][0]) && (Norm[2][2] < Epsilon[pg->imtrx][0]) &&
        (Norm[2][0] < Epsilon[pg->imtrx][0])) {
      *converged = (Epsilon[pg->imtrx][2] > 1.0) ? TRUE : FALSE;
    }
    if (!isfinite(Norm[0][0]) || !isfinite(Norm[0][1]) || !isfinite(Norm[0][2])) {
      DPRINTF(stderr, "\nNaN found in residual norms\n");
      return_value = -1;
      goto free_and_clear;
    }

#ifdef DEBUG_NORM
    if (fabs(resid_vector[num_unk_r]) == fabs(Norm[0][0])) {
      printf("P_%d: dofname[%d] = ", ProcID, num_unk_r);
      if (dofname)
        printf("%s", dofname[num_unk_r]);
      else
        printf(" (unk)");
      printf(", resid = %g\n", resid_vector[num_unk_r]);
    }
#endif /* DEBUG_NORM */

    /*
     * Before starting the solve of the linear system, fix up the matrix
     * to exclude external rows. However, we'll save all the column name
     * values so we can reincorporate them back in before the next assembly
     * step.
     */
    if (Num_Proc > 1 && strcmp(Matrix_Format, "msr") == 0) {
      hide_external(num_universe_dofs[pg->imtrx], NumUnknowns[pg->imtrx], ija, ija_save, a);
      dofs_hidden = TRUE;
    }

#ifdef DEBUG_JACOBIAN
    if (inewton < 2) {
      if (strcmp(Matrix_Format, "msr") == 0) {
        print_msr_matrix(num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx], ija, a, x);
        print_array(ija, ija[ija[0] - 1], "ija", type_int, ProcID);
      } else {
        print_vbr_matrix(ams, exo, dpi, NumUnknowns[pg->imtrx]);
      }
    }
#endif /* DEBUG_JACOBIAN */
       /*
        *       DUMP THE MATRIX if requested to do so
        */
#ifdef MATRIX_DUMP
    if (strcmp(Matrix_Format, "msr") == 0) {
      matrix_dump_msr(ams, exo, dpi, x);
    } else if (strcmp(Matrix_Format, "vbr") == 0) {
      matrix_dump_vbr(ams, exo, dpi, x);
    } else if (strcmp(Matrix_Format, "front") == 0) {
      if (ams->Number_Jac_Dump != 0) {
        WH(-1, "Matrix dump requested, but not supported for frontal solver");
      }
    }
#endif /* MATRIX_DUMP */
    /*************************************************************************
     *             SOLVE THE LINEAR SYSTEM
     *************************************************************************/
    s_start = ut();
    s_end = s_start;

    if (*converged) {
      /*  If we're going to skip the Solve step, compute and save some useful
          information first...*/

      Resid_Norm_stack[0] = Resid_Norm_stack[1];
      Resid_Norm_stack[1] = Resid_Norm_stack[2];
      Resid_Norm_stack[2] = Norm[0][2] / BIG_PENALTY;
      if (nAC) {
        AC_Resid_Norm_stack[0] = AC_Resid_Norm_stack[1];
        AC_Resid_Norm_stack[1] = AC_Resid_Norm_stack[2];
        AC_Resid_Norm_stack[2] = Norm[2][2] / BIG_PENALTY;
      }

      if (inewton && (Resid_Norm_stack[2] > 0) && (Resid_Norm_stack[2] != 1) &&
          (Resid_Norm_stack[1] > 0) && (Resid_Norm_stack[1] != 1) &&
          (inewton <= 1 || ((Resid_Norm_stack[0] > 0) && (Resid_Norm_stack[0] != 1)))) {
#if 1
        if (inewton <= 1) {
          Conv_order = log10(Resid_Norm_stack[2]) / log10(Resid_Norm_stack[1]);
          Conv_rate = log10(Resid_Norm_stack[1]) - log10(Resid_Norm_stack[2]);
        } else {
          Conv_order = log10(Resid_Norm_stack[2] / Resid_Norm_stack[1]) /
                       log10(Resid_Norm_stack[1] / Resid_Norm_stack[0]);
          Conv_rate = -0.5 * log10(Resid_Norm_stack[0]) + 2 * log10(Resid_Norm_stack[1]) -
                      1.5 * log10(Resid_Norm_stack[2]);
        }
#else
        Conv_order = log10(Resid_Norm_stack[2]) / log10(Resid_Norm_stack[1]);
        Conv_rate = log10(Resid_Norm_stack[1]) - log10(Resid_Norm_stack[2]);
#endif
      }

      if (nAC && inewton && (AC_Resid_Norm_stack[2] > 0) && (AC_Resid_Norm_stack[2] != 1) &&
          (AC_Resid_Norm_stack[1] > 0) && (AC_Resid_Norm_stack[1] != 1) &&
          (inewton <= 1 || ((AC_Resid_Norm_stack[0] > 0) && (AC_Resid_Norm_stack[0] != 1)))) {
#if 1
        if (inewton <= 1) {
          AConv_order = log10(AC_Resid_Norm_stack[2]) / log10(AC_Resid_Norm_stack[1]);
          AConv_rate = log10(AC_Resid_Norm_stack[1]) - log10(AC_Resid_Norm_stack[2]);
        } else if (DOUBLE_NONZERO(1.0 - AC_Resid_Norm_stack[1] / AC_Resid_Norm_stack[0])) {
          AConv_order = log10(AC_Resid_Norm_stack[2] / AC_Resid_Norm_stack[1]) /
                        log10(AC_Resid_Norm_stack[1] / AC_Resid_Norm_stack[0]);
          AConv_rate = -0.5 * log10(AC_Resid_Norm_stack[0]) + 2 * log10(AC_Resid_Norm_stack[1]) -
                       1.5 * log10(AC_Resid_Norm_stack[2]);
        } else {
          AConv_order = log10(AC_Resid_Norm_stack[2]) / log10(AC_Resid_Norm_stack[1]);
          AConv_rate = log10(AC_Resid_Norm_stack[1]) - log10(AC_Resid_Norm_stack[2]);
        }
#else
        AConv_order = log10(AC_Resid_Norm_stack[2]) / log10(AC_Resid_Norm_stack[1]);
        AConv_rate = log10(AC_Resid_Norm_stack[1]) - log10(AC_Resid_Norm_stack[2]);
#endif
      }
      if (glob_var_vals != NULL) {
        glob_var_vals[0] = (double)*converged;
        glob_var_vals[1] = (double)inewton;
        glob_var_vals[2] = (double)Max_Newton_Steps;
        glob_var_vals[3] = Conv_order;
        glob_var_vals[4] = Conv_rate;
        total_mesh_volume = 0;

        for (i = 0; nAC > 0 && i < nAC; i++) {
          total_mesh_volume += augc[i].evol;
          glob_var_vals[6 + i] = x_AC[i];
        }

        glob_var_vals[5] = (double)total_mesh_volume;
      }
      solve_skipped = true;
      goto skip_solve;
    }

    switch (Linear_Solver) {
    case UMFPACK2:
    case UMFPACK2F:
      if (strcmp(Matrix_Format, "msr"))
        GOMA_EH(GOMA_ERROR, "ERROR: umfpack solver needs msr matrix format");

      if (!Norm_below_tolerance || !Rate_above_tolerance)
        Factor_Flag = 1;
      else
        Factor_Flag = 3;

      if (first_linear_solver_call) {
        Factor_Flag = 0;
        UMF_system_id = -1;
      }

      /* Force refactorization if UMFPACK2F */
      if (Linear_Solver == UMFPACK2F) {
        Factor_Flag = 0;
      }
      matr_form = 1;

      UMF_system_id = SL_UMF(UMF_system_id, &first_linear_solver_call, &Factor_Flag, &matr_form,
                             &NumUnknowns[pg->imtrx], &NZeros, &ija[0], &ija[0], &a[0],
                             &resid_vector[0], &delta_x[0]);

      first_linear_solver_call = FALSE;

      if (!Norm_below_tolerance || !Rate_above_tolerance)
        Factor_Flag = 1;
      else
        Factor_Flag = 3;

      LOCA_UMF_ID = UMF_system_id;
      strcpy(stringer, " 1 ");
      break;

    case SPARSE13a:
      if (strcmp(Matrix_Format, "msr")) {
        GOMA_EH(GOMA_ERROR, "ERROR: lu solver needs msr matrix format");
      }
      dcopy1(NumUnknowns[pg->imtrx], resid_vector, delta_x);
      if (!Norm_below_tolerance || !Rate_above_tolerance) {
        lu(NumUnknowns[pg->imtrx], NumExtUnknowns[pg->imtrx], NZeros, a, ija, delta_x,
           (first_linear_solver_call ? 1 : 2));
        first_linear_solver_call = FALSE;
      } else {
        lu(NumUnknowns[pg->imtrx], NumExtUnknowns[pg->imtrx], NZeros, a, ija, delta_x, 3);
        first_linear_solver_call = FALSE;
      }
      /*
       * Note that sl_lu has static variables to keep track of
       * first call or not.
       */
      strcpy(stringer, " 1 ");
      break;

    case AZTEC:
      /*
       * Initialization is now performed up in
       * solve_problem() in rf_solve.c, not here in the
       * midst of the Newton iteration.
       */
      if (first_linear_solver_call) {
        ams->options[AZ_pre_calc] = AZ_calc;
      } else {
        if (strcmp(Matrix_Factorization_Reuse, "calc") == 0) {
          /*
           * Gonna start from scratch even though I've cooked a
           * preconditioner in the kitchen all day? Well, then
           * you won't need the leftover pieces from all my
           * hard preparation last time around.
           */

          AZ_free_memory(ams->data_org[AZ_name]);

          ams->options[AZ_pre_calc] = AZ_calc;

        } else if (strcmp(Matrix_Factorization_Reuse, "recalc") == 0) {
          ams->options[AZ_pre_calc] = AZ_recalc;
        } else if (strcmp(Matrix_Factorization_Reuse, "reuse") == 0) {
          ams->options[AZ_pre_calc] = AZ_reuse;
        } else {
          GOMA_EH(GOMA_ERROR, "Unknown factorization reuse specification.");
        }
      }

      linear_solver_blk = 0;     /* count calls to AZ_solve() */
      num_linear_solve_blks = 1; /* upper limit to AZ_solve() calls */
      matrix_solved = FALSE;
      while ((!matrix_solved) && (linear_solver_blk < num_linear_solve_blks)) {
        /*
         * Someday the user may want to do fancy heuristics based
         * on all kinds of cost functions, artificial intelligence
         * neural networks, etc.
         *
         * For the linear system "Ax=b", we have
         *    A -- indx, bindx(ija), rpntr, cpntr, bpntr, val(a)
         *    x -- delta_x, newton correction vector
         *    b -- resid_vector, newton residual equation vector
         */
        if (!Norm_below_tolerance || !Rate_above_tolerance) {
          /* Save old A before Aztec rescales it */
          if (save_old_A)
            dcopy1(NZeros, ams->val, ams->val_old);
        } else {
          /*Recover last A*/
          if (save_old_A)
            dcopy1(NZeros, ams->val_old, ams->val);
        }

        AZ_solve(delta_x, resid_vector, ams->options, ams->params, ams->indx, ams->bindx,
                 ams->rpntr, ams->cpntr, ams->bpntr, ams->val, ams->data_org, ams->status,
                 ams->proc_config);

        first_linear_solver_call = FALSE;

        if (Debug_Flag > 0) {
          dump_aztec_status(ams->status);
        }

        why = (int)ams->status[AZ_why];
        aztec_stringer(why, ams->status[AZ_its], &stringer[0]);

        matrix_solved = (ams->status[AZ_why] == AZ_normal);
        linear_solver_blk++;
      }

      /* Necessary anymore?
       * 	  if (options[AZ_pre_calc] == AZ_calc) {
       * 	  FREE the memory used in storing preconditioner info
       *   - unless using the RE_USE option
       * 	      AZ_clear(data_org[AZ_name]);
       * 	    AZ_free_memory(data_org[AZ_name]);
       * 	  }
       */
      break;

    case AMESOS:

      if (strcmp(Matrix_Format, "msr") == 0) {
        amesos_solve_msr(Amesos_Package, ams, delta_x, resid_vector, 1, pg->imtrx);
      } else if (strcmp(Matrix_Format, "epetra") == 0) {
        amesos_solve_epetra(Amesos_Package, ams, delta_x, resid_vector, pg->imtrx);
      } else {
        GOMA_EH(
            GOMA_ERROR,
            " Sorry, only MSR and Epetra matrix formats are currently supported with the Amesos "
            "solver suite\n");
      }
      strcpy(stringer, " 1 ");
      break;

    case AZTECOO:
      if (strcmp(Matrix_Format, "epetra") == 0) {
        aztecoo_solve_epetra(ams, delta_x, resid_vector);
        why = (int)ams->status[AZ_why];
        aztec_stringer(why, ams->status[AZ_its], &stringer[0]);
        matrix_solved = (ams->status[AZ_why] == AZ_normal);
      } else {
        GOMA_EH(GOMA_ERROR,
                "Sorry, only Epetra matrix formats are currently supported with the AztecOO solver "
                "suite\n");
      }
      break;

#ifdef GOMA_ENABLE_PETSC
    case PETSC_SOLVER:
      if (strcmp(Matrix_Format, "petsc") == 0) {
        int its;
        petsc_solve(ams, delta_x, resid_vector, &its);
        exchange_dof(cx, dpi, delta_x, pg->imtrx);
        matrix_solved = 1;
        char itsstring[10];
        itsstring[9] = '\0';
        snprintf(itsstring, 9, "%d", its);
        strcpy(stringer, itsstring);
      } else {
        GOMA_EH(GOMA_ERROR,
                "Sorry, only petsc matrix formats are currently supported with the petsc solver\n");
      }
      break;
#endif
    case STRATIMIKOS:
      if (strcmp(Matrix_Format, "epetra") == 0) {
        int iterations;
        int err =
            stratimikos_solve(ams, delta_x, resid_vector, &iterations, Stratimikos_File, pg->imtrx);
        if (err) {
          GOMA_EH(err, "Error in stratimikos solve");
          check_parallel_error("Error in solve - stratimikos");
        }
        aztec_stringer(AZ_normal, iterations, &stringer[0]);
      } else {
        GOMA_EH(GOMA_ERROR,
                "Sorry, only Epetra matrix formats are currently supported with the Stratimikos "
                "interface\n");
      }
      break;

    default:
      GOMA_EH(GOMA_ERROR, "That linear solver package is not implemented.");
      break;
    }
    s_end = ut();
    /**************************************************************************
     *        END OF LINEAR SYSTEM SOLVE SECTION
     **************************************************************************/
    /*
     *
     *          SENSITIVITIES FOR AC's
     *
     */

    if (nAC > 0) {
      sc_start = ut();

      /*
       * LOOP OVER NUMBER OF AUGMENTING CONDITIONS
       */
      for (iAC = 0; iAC < nAC; iAC++) {

        switch (Linear_Solver) {
        case UMFPACK2:
        case UMFPACK2F:
          matr_form = 1;

          /* MMH: I don't believe that a refactorization is EVER
           * intended here.  However, that was the option selected
           * when Linear_Solver is UMFPACK2F.  I'm guessing that
           * although this isn't, in fact, what you want to happen, it
           * is what was previously done before I consolidated
           * UMFPACK2 and UMFPACK2F. */
          if (Linear_Solver == UMFPACK2F)
            Factor_Flag = 0;
          else
            Factor_Flag = 3;

          if (first_linear_solver_call)
            GOMA_EH(GOMA_ERROR, "Solving for AC's BEFORE a regular solve");

#ifdef DEBUG_SL_UMF
          printf("%s: entering SL_UMF for augmenting conditions solve\n", yo);
          fflush(stdout);
#endif /* DEBUG_SL_UMF */
          UMF_system_id = SL_UMF(UMF_system_id, &first_linear_solver_call, &Factor_Flag, &matr_form,
                                 &NumUnknowns[pg->imtrx], &NZeros, &ija[0], &ija[0], &a[0],
                                 &bAC[iAC][0], &wAC[iAC][0]);
#ifdef DEBUG_SL_UMF
          printf("%s: returning from SL_UMF\n", yo);
          fflush(stdout);
#endif /* DEBUG_SL_UMF */

          /*  */
          strcpy(stringer_AC, " 1 ");
          LOCA_UMF_ID = UMF_system_id;
          break;

        case SPARSE13a:
          dcopy1(NumUnknowns[pg->imtrx], &bAC[iAC][0], &tAC[0]);
          lu(NumUnknowns[pg->imtrx], NumExtUnknowns[pg->imtrx], NZeros, a, ija, &tAC[0], 3);
          dcopy1(NumUnknowns[pg->imtrx], &tAC[0], &wAC[iAC][0]);
          strcpy(stringer_AC, " 1 ");
          break;

        case AMESOS:
          if (strcmp(Matrix_Format, "msr") == 0) {
            amesos_solve_msr(Amesos_Package, ams, &wAC[iAC][0], &bAC[iAC][0], 0, pg->imtrx);
          } else if (strcmp(Matrix_Format, "epetra") == 0) {
            amesos_solve_epetra(Amesos_Package, ams, &wAC[iAC][0], &bAC[iAC][0], pg->imtrx);
          } else {
            GOMA_EH(GOMA_ERROR,
                    " Sorry, only MSR and Epetra matrix formats are currently supported with the "
                    "Amesos solver suite\n");
          }
          strcpy(stringer_AC, " 1 ");
          break;

        case AZTEC:
          /*
           * Initialization is now performed up in
           * solve_problem() in rf_solve.c, not here in the
           * midst of the Newton iteration.
           */
          ams->options[AZ_pre_calc] = (ams->options[AZ_keep_info] ? AZ_reuse : AZ_calc);

          linear_solver_blk = 0;     /* count calls to AZ_solve() */
          num_linear_solve_blks = 1; /* upper limit to AZ_solve() calls */
          matrix_solved = FALSE;
          while ((!matrix_solved) && (linear_solver_blk < num_linear_solve_blks)) {
            /*
             * Someday the user may want to do fancy heuristics based
             * on all kinds of cost functions, artificial intelligence
             * neural networks, etc.
             *
             * For the linear system "Ax=b", we have
             *    A -- indx, bindx(ija), rpntr, cpntr, bpntr, val(a)
             *    x -- delta_x, newton correction vector
             *    b -- resid_vector, newton residual equation vector
             */
            AZ_solve(&wAC[iAC][0], &bAC[iAC][0], ams->options, ams->params, ams->indx, ams->bindx,
                     ams->rpntr, ams->cpntr, ams->bpntr, ams->val, ams->data_org, ams->status,
                     ams->proc_config);

            if (Debug_Flag > 0) {
              dump_aztec_status(ams->status);
            }

            why = (int)ams->status[AZ_why];
            aztec_stringer(why, ams->status[AZ_its], &stringer_AC[0]);

            matrix_solved = (ams->status[AZ_why] == AZ_normal);
            linear_solver_blk++;
          }

          /* Necessary anymore?
           * 	  if (options[AZ_pre_calc] == AZ_calc) {
           * 	  FREE the memory used in storing preconditioner info
           *   - unless using the RE_USE option
           * 	      AZ_clear(data_org[AZ_name]);
           * 	    AZ_free_memory(data_org[AZ_name]);
           * 	  }
           */

          break;
        case AZTECOO:
          if (strcmp(Matrix_Format, "epetra") == 0) {
            aztecoo_solve_epetra(ams, &wAC[iAC][0], &bAC[iAC][0]);
            why = (int)ams->status[AZ_why];
            aztec_stringer(why, ams->status[AZ_its], &stringer_AC[0]);
            matrix_solved = (ams->status[AZ_why] == AZ_normal);
          } else {
            GOMA_EH(GOMA_ERROR,
                    "Sorry, only Epetra matrix formats are currently supported with the AztecOO "
                    "solver suite\n");
          }
          break;
        case STRATIMIKOS:
          if (strcmp(Matrix_Format, "epetra") == 0) {
            int iterations;
            int err = stratimikos_solve(ams, &wAC[iAC][0], &bAC[iAC][0], &iterations,
                                        Stratimikos_File, pg->imtrx);
            GOMA_EH(err, "Error in stratimikos solve");
            if (iterations == -1) {
              strcpy(stringer, "err");
            } else {
              aztec_stringer(AZ_normal, iterations, &stringer[0]);
            }
          } else {
            GOMA_EH(
                GOMA_ERROR,
                "Sorry, only Epetra matrix formats are currently supported with the Stratimikos "
                "interface\n");
          }
          break;
        default:
          GOMA_EH(GOMA_ERROR, "That linear solver package is not implemented.");
          break;
        }
      } /* END AC LOOP */

      /**************************************************************************
       *  END OF AUGMENTED CONDITIONS SECTION
       **************************************************************************/
      /*
       *
       *               SCHUR COMPLEMENT OF J IN [M]
       *
       */

      for (iAC = 0; iAC < nAC; iAC++) {
        for (jAC = 0; jAC < nAC; jAC++) {
          sAC[iAC][jAC] = dot_product(NumUnknowns[pg->imtrx], &cAC[iAC][0], &wAC[jAC][0]);
        }
      }

#ifdef PARALLEL
      if (Num_Proc > 1) {
        for (iAC = 0; iAC < nAC; iAC++) {
          for (jAC = 0; jAC < nAC; jAC++)
            tempAC[jAC] = sAC[iAC][jAC];
          MPI_Allreduce(tempAC, sAC[iAC], nAC, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        }
      }
#endif /* PARALLEL */

      for (iAC = 0; iAC < nAC; iAC++) {
        for (jAC = 0; jAC < nAC; jAC++) {
          sAC[iAC][jAC] = dAC[iAC][jAC] - sAC[iAC][jAC];
        }
      }

      /*
       *
       *    GET UPDATE ON AC's
       *
       */
      for (iAC = 0; iAC < nAC; iAC++) {
        tAC[iAC] = dot_product(NumUnknowns[pg->imtrx], &cAC[iAC][0], &delta_x[0]);
      }
#ifdef PARALLEL
      if (Num_Proc > 1) {
        for (jAC = 0; jAC < nAC; jAC++)
          tempAC[jAC] = tAC[jAC];
        MPI_Allreduce(tempAC, tAC, nAC, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      }
#endif /* PARALLEL */

      v2sum(nAC, &yAC[0], 1.0, &gAC[0], -1.0, &tAC[0]);

      for (iAC = 0; iAC < nAC; iAC++) {
        ordAC[iAC] = iAC;
      }
      for (i = 0; i < nAC; i++) {
        bigAC = 0.0;
        for (j = 0; j < nAC; j++) {
          if (fabs(sAC[i][j]) > bigAC) {
            bigAC = fabs(sAC[i][j]);
          }
        }
        if (bigAC == 0.0) {
          puts("\n E: AC: sAC Singular matrix.");
          return_value = -1;
          goto free_and_clear;
        }
        hAC[i] = 1.0 / bigAC;
      }
      for (j = 0; j < nAC; j++) {
        for (i = 0; i < j; i++) {
          sumAC = sAC[i][j];
          for (k = 0; k < i; k++) {
            sumAC -= sAC[i][k] * sAC[k][j];
          }
          sAC[i][j] = sumAC;
        }
        bigAC = 0.0;
        for (i = j; i < nAC; i++) {
          sumAC = sAC[i][j];
          for (k = 0; k < j; k++) {
            sumAC -= sAC[i][k] * sAC[k][j];
          }
          sAC[i][j] = sumAC;
          dumAC = hAC[i] * fabs(sumAC);
          if (dumAC >= bigAC) {
            bigAC = dumAC;
            imaxAC = i;
          }
        }
        if (j != imaxAC) {
          for (k = 0; k < nAC; k++) {
            dumAC = sAC[imaxAC][k];
            sAC[imaxAC][k] = sAC[j][k];
            sAC[j][k] = dumAC;
          }
          hAC[imaxAC] = hAC[j];
        }
        ordAC[j] = imaxAC;
        if (sAC[j][j] == 0.0) {
          sAC[j][j] = 1.0e-20;
        }

        dumAC = 1.0 / sAC[j][j];
        for (i = j + 1; i < nAC; i++) {
          sAC[i][j] *= dumAC;
        }
      }
      iiAC = -1;
      for (i = 0; i < nAC; i++) {
        k = ordAC[i];
        sumAC = yAC[k];
        yAC[k] = yAC[i];
        if (iiAC != -1) {
          for (j = iiAC; j < i; j++) {
            sumAC -= sAC[i][j] * yAC[j];
          }
        } else if (sumAC != 0.0) {
          iiAC = i;
        }
        yAC[i] = sumAC;
      }
      for (i = (nAC - 1); i >= 0; i--) {
        sumAC = yAC[i];
        if (i < (nAC - 1)) {
          for (j = i + 1; j < nAC; j++) {
            sumAC -= sAC[i][j] * yAC[j];
          }
        }
        yAC[i] = sumAC / sAC[i][i];
      }

      /*

      SOLVE FOR REMAINING UPDATES

      */

      for (i = 0; i < NumUnknowns[pg->imtrx]; i++) {
        tAC[i] = 0.0;
        for (iAC = 0; iAC < nAC; iAC++) {
          tAC[i] += wAC[iAC][i] * yAC[iAC];
        }
      }
      /**************************************************************************
       *          END OF SHUR COMPLIMENT SECTION
       **************************************************************************/
      /*
       *
       *   UPDATE GOMA UNKNOWNS DUE TO CHANGES IN AUGMENTED CONDITIONS
       *
       */
      for (i = 0; i < NumUnknowns[pg->imtrx]; i++) {
        delta_x[i] -= tAC[i];
      }
      sc_end = ut();

      /* Set arc length parameter update */
      if (alc_with_acs)
        delta_y_ALC = yAC[nAC - 1];

    } /* END OF if (nAC > 0) */
    /**************************************************************************
     *         END OF AUGMENTED SYSTEM SOLVE SECTION
     **************************************************************************/

    Norm[1][0] = Loo_norm(delta_x, NumUnknowns[pg->imtrx], &num_unk_x, dofname_x);
    Norm[1][1] = L1_norm(delta_x, NumUnknowns[pg->imtrx]);
    Norm[1][2] = L2_norm(delta_x, NumUnknowns[pg->imtrx]);
    Norm_r[0][0] = Loo_norm_r(delta_x, x, NumUnknowns[pg->imtrx], &num_unk_x, dofname_nr);
    Norm_r[0][1] = L1_norm_r(delta_x, x, NumUnknowns[pg->imtrx]);
    Norm_r[0][2] = L2_norm_r(delta_x, x, NumUnknowns[pg->imtrx]);

    if (nAC > 0) {
      Norm[3][0] = Loo_norm_1p(yAC, nAC, &num_unk_y, dofname_y);
      Norm[3][1] = L1_norm_1p(yAC, nAC);
      Norm[3][2] = L2_norm_1p(yAC, nAC);
      Norm_r[1][0] = Loo_norm_r_1p(yAC, x_AC, nAC, &num_unk_y, dofname_ry);
      Norm_r[1][1] = L1_norm_r_1p(yAC, x_AC, nAC);
      Norm_r[1][2] = L2_norm_r_1p(yAC, x_AC, nAC);
      AC_Resid_Norm_stack[0] = AC_Resid_Norm_stack[1];
      AC_Resid_Norm_stack[1] = AC_Resid_Norm_stack[2];
      AC_Resid_Norm_stack[2] = Norm[2][2] / BIG_PENALTY;
      AC_Soln_Norm_stack[0] = AC_Soln_Norm_stack[1];
      AC_Soln_Norm_stack[1] = AC_Soln_Norm_stack[2];
      AC_Soln_Norm_stack[2] = Norm_r[1][2] / BIG_PENALTY;
    }

    /* Save some norm info for modified Newton stuff */

    Resid_Norm_stack[0] = Resid_Norm_stack[1];
    Resid_Norm_stack[1] = Resid_Norm_stack[2];
    Resid_Norm_stack[2] = Norm[0][2] / BIG_PENALTY;
    Soln_Norm_stack[0] = Soln_Norm_stack[1];
    Soln_Norm_stack[1] = Soln_Norm_stack[2];
    Soln_Norm_stack[2] = Norm_r[0][2] / BIG_PENALTY;
    /* Compute rate_of_convergence.  Really what we want
                   here is dln(norm)/dnorm   FIGURE IT OUT! */
    /* Asymptotically we should have Norm_(i+1) ~ lambda * Norm_i ^ alpha
       or log[Norm_(i+1)] = log(lambda) + alpha*log(Norm_i),
       where alpha defines the convergence order.
       For <= 2 iterations, assume lambda=1 ==> alpha = log[Norm_(i+1)]/log[Norm_i]
       Otherwise, eliminate lambda using 3 iteration results to get:
            alpha = log[Norm_(i+1)/Norm_i]/log[Norm_i/Norm_(i-1)].
    */

    if (inewton && (Resid_Norm_stack[2] > 0) && (Resid_Norm_stack[2] != 1) &&
        (Resid_Norm_stack[1] > 0) && (Resid_Norm_stack[1] != 1) &&
        (inewton <= 1 || ((Resid_Norm_stack[0] > 0) && (Resid_Norm_stack[0] != 1)))) {
#if 1
      if (inewton <= 1) {
        Conv_order = log10(Resid_Norm_stack[2]) / log10(Resid_Norm_stack[1]);
        Soln_order = log10(Soln_Norm_stack[2]) / log10(Soln_Norm_stack[1]);
        Conv_rate = log10(Resid_Norm_stack[1]) - log10(Resid_Norm_stack[2]);
        Soln_rate = log10(Soln_Norm_stack[1]) - log10(Soln_Norm_stack[2]);
      } else {
        Conv_order = log10(Resid_Norm_stack[2] / Resid_Norm_stack[1]) /
                     log10(Resid_Norm_stack[1] / Resid_Norm_stack[0]);
        Soln_order = log10(Soln_Norm_stack[2] / Soln_Norm_stack[1]) /
                     log10(Soln_Norm_stack[1] / Soln_Norm_stack[0]);
        Conv_rate = -0.5 * log10(Resid_Norm_stack[0]) + 2 * log10(Resid_Norm_stack[1]) -
                    1.5 * log10(Resid_Norm_stack[2]);
        Soln_rate = -0.5 * log10(Soln_Norm_stack[0]) + 2 * log10(Soln_Norm_stack[1]) -
                    1.5 * log10(Soln_Norm_stack[2]);
      }
#else
      Conv_order = log10(Resid_Norm_stack[2]) / log10(Resid_Norm_stack[1]);
      Soln_order = log10(Soln_Norm_stack[2]) / log10(Soln_Norm_stack[1]);
      Conv_rate = log10(Resid_Norm_stack[1]) - log10(Resid_Norm_stack[2]);
      Soln_rate = log10(Soln_Norm_stack[1]) - log10(Soln_Norm_stack[2]);
#endif
    }
    if (nAC && inewton && (AC_Resid_Norm_stack[2] > 0) && (AC_Resid_Norm_stack[2] != 1) &&
        (AC_Resid_Norm_stack[1] > 0) && (AC_Resid_Norm_stack[0] != 1) &&
        (inewton <= 1 || ((AC_Resid_Norm_stack[0] > 0) && (AC_Resid_Norm_stack[0] != 1)))) {
#if 1
      if (inewton <= 1 && DOUBLE_NONZERO(log10(AC_Resid_Norm_stack[1]))) {
        AConv_order = log10(AC_Resid_Norm_stack[2]) / log10(AC_Resid_Norm_stack[1]);
        ASoln_order = log10(AC_Soln_Norm_stack[2]) / log10(AC_Soln_Norm_stack[1]);
        AConv_rate = log10(AC_Resid_Norm_stack[1]) - log10(AC_Resid_Norm_stack[2]);
        ASoln_rate = log10(AC_Soln_Norm_stack[1]) - log10(AC_Soln_Norm_stack[2]);
      } else if (AC_Resid_Norm_stack[0] > 0 && AC_Resid_Norm_stack[1] > 0 &&
                 DOUBLE_NONZERO(log10(AC_Resid_Norm_stack[1] / AC_Resid_Norm_stack[0]))) {
        AConv_order = log10(AC_Resid_Norm_stack[2] / AC_Resid_Norm_stack[1]) /
                      log10(AC_Resid_Norm_stack[1] / AC_Resid_Norm_stack[0]);
        ASoln_order = log10(AC_Soln_Norm_stack[2] / AC_Soln_Norm_stack[1]) /
                      log10(AC_Soln_Norm_stack[1] / AC_Soln_Norm_stack[0]);
        AConv_rate = -0.5 * log10(AC_Resid_Norm_stack[0]) + 2 * log10(AC_Resid_Norm_stack[1]) -
                     1.5 * log10(AC_Resid_Norm_stack[2]);
        ASoln_rate = -0.5 * log10(AC_Soln_Norm_stack[0]) + 2 * log10(AC_Soln_Norm_stack[1]) -
                     1.5 * log10(AC_Soln_Norm_stack[2]);
      }
#else
      AConv_order = log10(AC_Resid_Norm_stack[2]) / log10(AC_Resid_Norm_stack[1]);
      ASoln_order = log10(AC_Soln_Norm_stack[2]) / log10(AC_Soln_Norm_stack[1]);
      AConv_rate = log10(AC_Resid_Norm_stack[1]) - log10(AC_Resid_Norm_stack[2]);
      ASoln_rate = log10(AC_Soln_Norm_stack[1]) - log10(AC_Soln_Norm_stack[2]);
#endif
    }

    /* fail if we didn't get a finite solution */
    if (!isfinite(Norm[1][0]) || !isfinite(Norm[1][1]) || !isfinite(Norm[1][2]) ||
        !isfinite(Norm[0][0]) || !isfinite(Norm[0][1]) || !isfinite(Norm[0][2])) {
      return_value = -1;
      goto free_and_clear;
    }

    if (modified_newton && modified_newt_norm_tol != 0. && convergence_rate_tolerance != 0.) {

      /* Compute convergence level and compare with tolerance */
      Norm_below_tolerance = (Resid_Norm_stack[2] < modified_newt_norm_tol);

      Rate_above_tolerance = (Conv_rate > convergence_rate_tolerance);
    }

    if (Newt_Jacobian_Reformation_stride > 1) {
      if (inewton == step_reform - 1) {
        /* here we are simply fooling the tolerance checker to
           reform based on iteration number or time step number */
        Norm_below_tolerance = FALSE;
        Rate_above_tolerance = FALSE;
        step_reform += Newt_Jacobian_Reformation_stride;
      } else {
        Norm_below_tolerance = TRUE;
        Rate_above_tolerance = TRUE;
      }
    } else if (Time_Jacobian_Reformation_stride != 1) {
      /* do nothing different*/
    }

    log_msg("%-38s = %23.16e", "correction norm (L_oo)", Norm[1][0]);
    log_msg("%-38s = %23.16e", "correction norm (L_1)", Norm[1][1]);
    log_msg("%-38s = %23.16e", "correction norm (L_2)", Norm[1][2]);

    if (Write_Intermediate_Solutions || (Iout == 1)) {

#define debug_subelement_decomposition 1
#if debug_subelement_decomposition
      if (ls != NULL && ls->SubElemIntegration)
        subelement_mesh_output(x, exo);
#endif

      /*
       * HKM -> should add a few sync() statements here
       */
      sync_processors();
      print_sync_start(FALSE);

      print_sync_end(FALSE);
    }

    fflush(NULL);

    /*
     * If the solution is to be written at intermediate stages of the
     * Newton iteration, the new default is to write out the initial
     * guess immediately BEFORE the first Newton iteration. This is
     * Polly's idea and, of course, it's a good one. This way, the
     * Residual you see corresponds to the x that was input, including
     * the initial guess.
     */

    if (Write_Intermediate_Solutions && Unlimited_Output) {
      error = wr_soln_vec(x, resid_vector, NumUnknowns[pg->imtrx], inewton);
      GOMA_EH(error, "problem from wr_soln_vec");
    }

    if ((damp_factor1 <= 1. && damp_factor1 >= 0.) && (damp_factor2 <= 1. && damp_factor2 >= 0.) &&
        (damp_factor3 <= 1. && damp_factor3 >= 0.)) {
      print_damp_factor = TRUE;
      if (!Visc_Sens_Copy)
        Include_Visc_Sens = FALSE;
      if (Norm[0][0] > custom_tol3) {
        damp_factor = damp_factor3;
        /*
                        DPRINTF(stderr, "\n Invoking Damping factor %f\n", damp_factor3);
        */
      } else if (Norm[0][0] > custom_tol2) {
        damp_factor = damp_factor2;
        /*
                        DPRINTF(stderr, "\n Invoking Damping factor %f\n", damp_factor2);
        */
      } else if (Norm[0][0] > custom_tol1) {
        damp_factor = damp_factor1;
        /*
                        DPRINTF(stderr, "\n Invoking Damping factor %f\n", damp_factor1);
        */
      } else {
        damp_factor = 1.00;
        /*
                        DPRINTF(stderr, " \n Invoking Damping factor %f", damp_factor);
        */
        if (!Visc_Sens_Copy) {
          Include_Visc_Sens = TRUE;
          print_visc_sens = TRUE;
          /*
                                   DPRINTF(stderr, " Invoking Viscosity Sensitivities");
          */
        }
        /*
                        DPRINTF(stderr, "\n");
        */
      }
    } else {
      /*default damping factor case */
      if (damp_factor2 == -1.)
        damp_factor = damp_factor1;
      if (!Visc_Sens_Copy) {
        if (Visc_Sens_Factor * inewton < Max_Newton_Steps) {
          Include_Visc_Sens = FALSE;
        } else {
          Include_Visc_Sens = TRUE;
          print_visc_sens = TRUE;
          /*
                                   DPRINTF(stderr, " Invoking Viscosity Sensitivities\n");
          */
        }
      }
    }

    if (damp_factor <= 1.0e-06) {
      damp_factor = damp_factor1;
    }
    dcopy1(numProcUnknowns, delta_x, x_update);

    /*
     * For LOCA continuation algorithms, invoke any
     * needed bordering algorithms here.
     * If arc length algorithm has been done as an AC,
     * just get the convergence status.
     */

    if (con_ptr != NULL) {
      if (alc_with_acs)
        continuation_converged =
            arc_length_status(con_ptr, equation_ALC, delta_y_ALC, Reltol, Abstol);
      else
        continuation_converged = continuation_hook(x, delta_x, con_ptr, Reltol, Abstol);
    }

    /*******************************************************************
     *
     *   UPDATE GOMA UNKNOWNS
     *
     *******************************************************************/
    for (i = 0; i < NumUnknowns[pg->imtrx]; i++) {
      x[i] -= damp_factor * var_damp[idv[pg->imtrx][i][0]] * delta_x[i];
    }
    exchange_dof(cx, dpi, x, pg->imtrx);
    if (pd->TimeIntegration != STEADY) {
      for (i = 0; i < NumUnknowns[pg->imtrx]; i++) {
        xdot[i] -=
            damp_factor * var_damp[idv[pg->imtrx][i][0]] * delta_x[i] * (1.0 + 2 * theta) / delta_t;
      }
      exchange_dof(cx, dpi, xdot, pg->imtrx);

      /* Now go back and correct all those dofs in solid regions undergoing newmark-beta
       * transient scheme */
      if (tran->solid_inertia) {
        beta = tran->newmark_beta;
        gamma = tran->newmark_gamma;
        for (i = 0; i < num_total_nodes; i++) {
          node_ptr = Nodes[i];
          curr_mat_list = &(node_ptr->Mat_List);
          num_mat = curr_mat_list->Length;

          for (imat = 0; imat < num_mat; imat++) {
            mat_index = (curr_mat_list->List)[imat];
            if (pd_glob[mat_index]->MeshMotion == DYNAMIC_LAGRANGIAN) {
              for (b = 0; b < ei[pg->imtrx]->ielem_dim; b++) {
                j = Index_Solution(i, R_MESH1 + b, 0, 0, -1, pg->imtrx);
                xdot[j] = (x[j] - x_old[j]) * (gamma / beta) / delta_t +
                          tran->xdbl_dot_old[j] *
                              ((1 - gamma) - gamma * (1. - 2. * beta) / 2. / beta) * delta_t +
                          xdot_old[j] * (1. - gamma / beta);

                tran->xdbl_dot[j] = (x[j] - x_old[j]) / beta / delta_t / delta_t -
                                    xdot_old[j] / beta / delta_t -
                                    tran->xdbl_dot_old[j] * (1. - 2. * beta) / 2. / beta;
              }
            }
          }
        }
        exchange_dof(cx, dpi, xdot, pg->imtrx);
        exchange_dof(cx, dpi, tran->xdbl_dot, pg->imtrx);
      }

      /* Now go back and correct all those dofs that use XFEM */
      if (xfem != NULL) {
        if (upd->Total_Num_Matrices > 0) {
          if (pg->imtrx == ls->MatrixNum) {
            xfem_correct(num_total_nodes, x, xdot, x_old, xdot_old, delta_x, theta, delta_t);
          } else {
            xfem_correct(num_total_nodes, x, xdot, x_old, xdot_old, NULL, theta, delta_t);
          }
        } else {
          xfem_correct(num_total_nodes, x, xdot, x_old, xdot_old, delta_x, theta, delta_t);
        }
        exchange_dof(cx, dpi, x, pg->imtrx);
        exchange_dof(cx, dpi, xdot, pg->imtrx);
      }
    }

    /* implicit embedded level set surfaces*/
    if (ls != NULL && ls->Evolution == LS_EVOLVE_SLAVE) {
      surf_based_initialization(x, delta_x, xdot, exo, num_total_nodes, ls->init_surf_list,
                                time_value, theta, delta_t);
      exchange_dof(cx, dpi, x, pg->imtrx);
      exchange_dof(cx, dpi, xdot, pg->imtrx);
    }
    if (pfd != NULL) {
      ls_old = ls;
      ls = pfd->ls[0];

      if (ls->Evolution == LS_EVOLVE_SLAVE) {
        surf_based_initialization(x, delta_x, xdot, exo, num_total_nodes, ls->init_surf_list,
                                  time_value, theta, delta_t);
        exchange_dof(cx, dpi, x, pg->imtrx);
        exchange_dof(cx, dpi, xdot, pg->imtrx);
      }
      ls = ls_old;
    }

    /*
     *
     *   UPDATE THE AUGMENTED UNKNOWNS
     *
     */
    if (nAC > 0) {
      if (Debug_Flag > 0) {
        DPRINTF(stdout, "\n------------------------------\n");
        DPRINTF(stdout, "Augmenting Conditions:    %4d\n", nAC);
        DPRINTF(stdout, "Number of extra unknowns: %4d\n\n", nAC);
      }
      for (iAC = 0; iAC < nAC; iAC++) {
        x_AC[iAC] -= damp_factor * yAC[iAC];
        if (pd->TimeIntegration != STEADY) {
          x_AC_dot[iAC] -= damp_factor * yAC[iAC] * (1.0 + 2 * theta) / delta_t;
          augc[iAC].tmp2 = x_AC_dot[iAC];
        }
        update_parameterAC(iAC, x, xdot, x_AC, cx, exo, dpi);
        augc[iAC].tmp1 = x_AC[iAC];

        /*
         * PRINT OUT VALUES OF EXTRA UNKNOWNS
         * FROM AUGMENTING CONDITIONS
         */
        if (Debug_Flag > 0) {
          if (augc[iAC].Type == AC_USERBC) {
            DPRINTF(stdout, "\tBC[%4d] DF[%4d]=% 10.6e Update=% 10.6e\n", augc[iAC].BCID,
                    augc[iAC].DFID, x_AC[iAC], damp_factor * yAC[iAC]);
          } else if (augc[iAC].Type == AC_USERMAT || augc[iAC].Type == AC_FLUX_MAT) {
            DPRINTF(stdout, "\tMT[%4d] MP[%4d]=% 10.6e Update=% 10.6e\n", augc[iAC].MTID,
                    augc[iAC].MPID, x_AC[iAC], damp_factor * yAC[iAC]);
          } else if (augc[iAC].Type == AC_VOLUME) {
            DPRINTF(stdout, "\tMT[%4d] VC[%4d]=%10.6e Param=%10.6e\n", augc[iAC].MTID,
                    augc[iAC].VOLID, augc[iAC].evol, x_AC[iAC]);
          } else if (augc[iAC].Type == AC_FLUX) {
            DPRINTF(stdout, "\tBC[%4d] DF[%4d]=% 10.6e Update=% 10.6e\n", augc[iAC].BCID,
                    augc[iAC].DFID, x_AC[iAC], damp_factor * yAC[iAC]);
          } else if (augc[iAC].Type == AC_LGRM) {
            DPRINTF(stdout, "\tAC[%d], Lagrange Multiplier=%10.6e Update=%10.6e\n", iAC, x_AC[iAC],
                    damp_factor * yAC[iAC]);
          } else if (augc[iAC].Type == AC_ARC_LENGTH) {
            DPRINTF(stdout, "\tAC[%d], Arc Length Parameter=%10.6e Update=%10.6e\n", iAC, x_AC[iAC],
                    damp_factor * yAC[iAC]);
          } else if (augc[iAC].Type == AC_OVERLAP) {
            DPRINTF(stdout, "\tAC[%d], Elem %d Side %d  Dim %d:  LM=%10.6e  Update=%10.6e\n", iAC,
                    augc[iAC].lm_elem, augc[iAC].lm_side, augc[iAC].lm_dim, x_AC[iAC],
                    damp_factor * yAC[iAC]);
          } else if (augc[iAC].Type == AC_PERIODIC) {
            DPRINTF(stdout, "\tAC[%d], Elem %d Side %d  Var %s:  LM=%10.6e  Update=%10.6e\n", iAC,
                    augc[iAC].lm_elem, augc[iAC].lm_side, Var_Name[augc[iAC].VAR].name1, x_AC[iAC],
                    damp_factor * yAC[iAC]);
          } else if (augc[iAC].Type == AC_POSITION) {
            DPRINTF(stdout, "\tMT[%4d] XY[%4d]=%10.6e Param=%10.6e\n", augc[iAC].MTID,
                    augc[iAC].VOLID, augc[iAC].evol, x_AC[iAC]);
          } else if (augc[iAC].Type == AC_ANGLE) {
            DPRINTF(stdout, "\tMT[%4d] ANG[%4d]=%10.6e Param=%10.6e\n", augc[iAC].VOLID,
                    augc[iAC].MTID, augc[iAC].evol, x_AC[iAC]);
          }
        }
      }
    }
    /********************************************************************
     *  CHECK IF CONVERGED
     *   (EDW: Modified to check bordering algorithm convergence
     *         as required by LOCA)
     ********************************************************************/

    *converged =
        ((Norm[0][2] < Epsilon[pg->imtrx][0]) && (Norm[2][2] < Epsilon[pg->imtrx][0]) &&
         (Norm[0][0] < Epsilon[pg->imtrx][0]) && (Norm[0][2] < Epsilon[pg->imtrx][0]) &&
         (Norm_r[0][2] < Epsilon[pg->imtrx][2]) && (Norm_r[1][2] < Epsilon[pg->imtrx][2]) &&
         (Norm_r[0][0] < Epsilon[pg->imtrx][2]) && (Norm_r[1][0] < Epsilon[pg->imtrx][2]) &&
         (continuation_converged));

    /*******************************************************************
     *
     * copy out global unknowns if array as been allocated for them
     *
     *
     *******************************************************************/

    if (glob_var_vals != NULL) {
      glob_var_vals[0] = (double)*converged;
      glob_var_vals[1] = (double)inewton;
      glob_var_vals[2] = (double)Max_Newton_Steps;
      glob_var_vals[3] = Conv_order;
      glob_var_vals[4] = Conv_rate;
      total_mesh_volume = 0;

      for (i = 0; nAC > 0 && i < nAC; i++) {
        total_mesh_volume += augc[i].evol;
        glob_var_vals[6 + i] = x_AC[i];
      }

      glob_var_vals[5] = (double)total_mesh_volume;
    }

    /********************************************************************
     *
     *    OPTIONALLY,
     *    WRITE OUT THE INTERMEDIATE SOLUTION AT EACH NEWTON ITERATION
     *
     ********************************************************************/

  skip_solve:

    if (Write_Intermediate_Solutions) {
      if (TimeIntegration == STEADY) {
        time_value = (double)*nprint + 1.;
      }

      /* First nodal vars */
      for (i = 0; i < rd->TotalNVSolnOutput; i++) {
        extract_nodal_vec(x, rd->nvtype[i], rd->nvkind[i], rd->nvmatID[i], gvec, exo, FALSE,
                          time_value);
        wr_nodal_result_exo(exo, ExoFileOut, gvec, i + 1, *nprint + 1, time_value);
      }
      /* Add additional user-specified post processing variables */
      if (rd->TotalNVPostOutput > 0) {
        post_process_nodal(x, x_sens_p, x_old, xdot, xdot_old, resid_vector, *nprint + 1,
                           &time_value, delta_t, theta, NULL, exo, dpi, rd, ExoFileOut, pg->imtrx);
      }
      /* Write out time derivatives if requested */
      if (TIME_DERIVATIVES != -1 && (TimeIntegration != STEADY)) {
        for (i = 0; i < rd->TotalNVSolnOutput; i++) {
          i_post = rd->TotalNVSolnOutput + rd->TotalNVPostOutput + i;
          extract_nodal_vec(xdot, rd->nvtype[i_post], rd->nvkind[i_post], rd->nvmatID[i_post], gvec,
                            exo, TRUE, time_value);
          wr_nodal_result_exo(exo, ExoFileOut, gvec, i_post + 1, *nprint + 1, time_value);
        }
      }

      /* Now element vars */
      for (i = 0; i < tev; i++) {
        extract_elem_vec(x, i, rd->evtype[i], gvec_elem, exo, 0);
        wr_elem_result_exo(exo, ExoFileOut, gvec_elem, i, *nprint + 1, time_value, rd);
      }
      /* Add additional user-specified post processing variables */
      if (tev_post > 0) {
        post_process_elem(x, x_old, xdot, xdot_old, resid_vector, tev, tev_post, gvec_elem,
                          *nprint + 1, &time_value, delta_t, exo, dpi, rd);

        /* Write out time derivatives if requested */
        if (TIME_DERIVATIVES != -1 && (TimeIntegration != STEADY)) {
          for (i = 0; i < tev; i++) {
            i_post = tev_post + i;
            extract_elem_vec(xdot, i_post, rd->evtype[i_post], gvec_elem, exo, 0);
            wr_elem_result_exo(exo, ExoFileOut, gvec_elem, i_post, *nprint + 1, time_value, rd);
          }
        }
      }

      /* and global variables.  The case glob_var_vals == NULL is caught in the function */
      wr_global_result_exo(exo, ExoFileOut, *nprint + 1, 6 + nAC, glob_var_vals);

      *nprint = *nprint + 1;
    }

    /*
     * Back to normal. Readjust the ija[] pointers to mimic the bloated
     * system for this processor that pretends like it really cares about
     * filling in matrix entries for external degrees of freedom.
     *
     * If we are doing continuation then after convergence one must
     * defer this event until after performing the sensitivity calculation
     */
    if (Num_Proc > 1 && strcmp(Matrix_Format, "msr") == 0) {
      if ((!*converged) &&
          ((Continuation != ALC_NONE) || (nn_post_fluxes_sens > 0) || (nn_post_data_sens > 0))) {
        show_external(num_universe_dofs[pg->imtrx],
                      (num_universe_dofs[pg->imtrx] - num_external_dofs[pg->imtrx]), ija, ija_save,
                      a);
        dofs_hidden = FALSE;
      }
    }

    if (solve_skipped) {
      if (Solver_Output_Format & 32)
        DPRINTF(stdout, "        ");
      if (Solver_Output_Format & 64)
        DPRINTF(stdout, "        ");
      if (Solver_Output_Format & 128)
        DPRINTF(stdout, "        ");
      if (Solver_Output_Format & 256)
        DPRINTF(stdout, "skp ");
    } else {
      if (Epsilon[pg->imtrx][2] > 1) {
        if (Solver_Output_Format & 32)
          DPRINTF(stdout, "%7.1e ", Norm[1][0]);
        if (Solver_Output_Format & 64)
          DPRINTF(stdout, "%7.1e ", Norm[1][1]);
        if (Solver_Output_Format & 128)
          DPRINTF(stdout, "%7.1e ", Norm[1][2]);
        if (Solver_Output_Format & 256)
          DPRINTF(stdout, "%s ", stringer);
      } else {
        if (Solver_Output_Format & 32)
          DPRINTF(stdout, "%7.1e ", Norm_r[0][0]);
        if (Solver_Output_Format & 64)
          DPRINTF(stdout, "%7.1e ", Norm_r[0][1]);
        if (Solver_Output_Format & 128)
          DPRINTF(stdout, "%7.1e ", Norm_r[0][2]);
        if (Solver_Output_Format & 256)
          DPRINTF(stdout, "%s ", stringer);
      }
    }

    asmslv_time = (a_end - a_start);
    slv_time = (s_end - s_start);
    if (Solver_Output_Format & 512)
      DPRINTF(stdout, "%7.1e/%7.1e ", asmslv_time, slv_time);
    if (Solver_Output_Format & 1024)
      DPRINTF(stdout, "%.1f/%.1f ", Conv_order, Conv_rate);
    if (Solver_Output_Format & 2048)
      DPRINTF(stdout, "%.1f/%.1f ", Soln_order, Soln_rate);
    DPRINTF(stdout, "\n");

    if (Write_Intermediate_Solutions || (Iout == 1)) {
      if (dofname_r[0] != '\0') {
        fprintf(stdout, "L_oo cause: R->(%s)", dofname_r);
      }
      if (dofname_x[0] != '\0') {
        fprintf(stdout, "DelX->(%s)\n", dofname_x);
      }
    }

    /* print damping factor and/or viscosity sens message here */
    if (print_damp_factor)
      DPRINTF(stdout, " Invoking damping factor %f\n", damp_factor);
    if (print_visc_sens)
      DPRINTF(stdout, " Invoking Viscosity Sensitivities\n");

    if (nAC > 0) {
      if (Solver_Output_Format & 1)
        DPRINTF(stdout, "         ");
      if (Solver_Output_Format & 2)
        DPRINTF(stdout, " AC ");
      if (Solver_Output_Format & 4)
        DPRINTF(stdout, "%7.1e ", Norm[2][0]);
      if (Solver_Output_Format & 8)
        DPRINTF(stdout, "%7.1e ", Norm[2][1]);
      if (Solver_Output_Format & 16)
        DPRINTF(stdout, "%7.1e ", Norm[2][2]);
      if (Epsilon[pg->imtrx][2] > 1) {
        if (Solver_Output_Format & 32)
          DPRINTF(stdout, "%7.1e ", Norm[3][0]);
        if (Solver_Output_Format & 64)
          DPRINTF(stdout, "%7.1e ", Norm[3][1]);
        if (Solver_Output_Format & 128)
          DPRINTF(stdout, "%7.1e ", Norm[3][2]);
      } else {
        if (Solver_Output_Format & 32)
          DPRINTF(stdout, "%7.1e ", Norm_r[1][0]);
        if (Solver_Output_Format & 64)
          DPRINTF(stdout, "%7.1e ", Norm_r[1][1]);
        if (Solver_Output_Format & 128)
          DPRINTF(stdout, "%7.1e ", Norm_r[1][2]);
      }
      if (Solver_Output_Format & 256)
        DPRINTF(stdout, "%s ", stringer_AC);
      asmslv_time = (ac_end - ac_start);
      slv_time = (sc_end - sc_start);
      if (Solver_Output_Format & 512)
        DPRINTF(stdout, "%7.1e/%7.1e ", asmslv_time, slv_time);
      if (Solver_Output_Format & 1024)
        DPRINTF(stdout, "%.1f/%.1f ", AConv_order, AConv_rate);
      if (Solver_Output_Format & 2048)
        DPRINTF(stdout, "%.1f/%.1f ", ASoln_order, ASoln_rate);
      DPRINTF(stdout, "\n");
    }

    inewton++;
    af->Sat_hyst_reevaluate = FALSE; /*only want this true
                                       for first iteration*/
  }                                  /* End of loop over newton iterations */

  /**********************************************************************/
  /**********************************************************************
   *
   * This is the end of Newton iteration loop. Either we have converged or
   * we have exceeded the permissable number of Newton iterations.
   * Here we will write out some of the convergence (or lack thereof) info
   * to the exodusII file so that intelligent decisions can be made regarding
   * automated continuation, from an external source like DAKOTA.
   * What we will write in the form of global variables in the exodus II file is
   *    If converged (0 or 1)
   *    if inewton =  Max_Newton_Steps (0 or 1)
   *    Convergence rate at last steps
   *
   *
   **********************************************************************/

  /**  return number of newton iterations  **/
  return_value = inewton;

  if (!*converged) {
    if (Debug_Flag) {
      DPRINTF(stdout, "\n%s:  Newton iteration FAILED.\n", yo);
    }
  } else {
    if (Debug_Flag) {
      DPRINTF(stdout, "\n%s:  Newton iteration CONVERGED.\n", yo);
    }
  }

  /**
   *     Solution vector norms for detecting turning points
   */
  Norm[4][0] = Loo_norm(x, NumUnknowns[pg->imtrx], &num_unk_x, dofname_x);
  Norm[4][1] = L1_norm(x, NumUnknowns[pg->imtrx]) / ((double)NumUnknowns[pg->imtrx]);
  Norm[4][2] = L2_norm(x, NumUnknowns[pg->imtrx]) / sqrt((double)NumUnknowns[pg->imtrx]);

  if (Solver_Output_Format & 1)
    DPRINTF(stdout, "-------- ");
  if (Solver_Output_Format & 2)
    DPRINTF(stdout, "--- ");
  if (Solver_Output_Format & 4)
    DPRINTF(stdout, "------- ");
  if (Solver_Output_Format & 8)
    DPRINTF(stdout, "------- ");
  if (Solver_Output_Format & 16)
    DPRINTF(stdout, "------- ");
  if (Solver_Output_Format & 32)
    DPRINTF(stdout, "------- ");
  if (Solver_Output_Format & 64)
    DPRINTF(stdout, "------- ");
  if (Solver_Output_Format & 128)
    DPRINTF(stdout, "------- ");
  if (Solver_Output_Format & 256)
    DPRINTF(stdout, "--- ");
  if (Solver_Output_Format & 512)
    DPRINTF(stdout, "--------------- ");
  if (Solver_Output_Format & 1024)
    DPRINTF(stdout, "------- ");
  if (Solver_Output_Format & 2048)
    DPRINTF(stdout, "-------");
  DPRINTF(stdout, "\n");
  DPRINTF(stdout, "scaled solution norms  %10.3e %10.3e %10.3e \n", Norm[4][0], Norm[4][1],
          Norm[4][2]);

  /*
   * COMPUTE SOLUTION SENSITIVITY TO PARAMETERS AS REQUESTED
   */

  sens_vec_ct = -1;

  /*
   *
   *   FLUX SENSITIVITIES
   *
   */
  mf_resolve = 0;
  ncod = NULL;
  bc = NULL;
  smallpiv = 0.;
  singpiv = 0.;
  iautopiv = 0;
  iscale = 0;
  scaling_max = 0.;

  for (i = 0; i < nn_post_fluxes_sens; i++) {
    if (pp_fluxes_sens[i]->vector_id > sens_vec_ct) {
      /*
       *
       *	GET SENSITIVITY PARAMETER VALUE
       *
       */
      retrieve_parameterS(&param_val, x, xdot, pp_fluxes_sens[i]->sens_type,
                          pp_fluxes_sens[i]->sens_id, pp_fluxes_sens[i]->sens_flt,
                          pp_fluxes_sens[i]->sens_flt2, cx, exo, dpi);

      strcat(sens_caller, "Flux Sensitivity");

      err = soln_sens(
          param_val, x, xdot, delta_t, exo, dpi, cx, res_p, numProcUnknowns, ija, a, x_old, x_older,
          xdot_old, x_update, delta_t, theta, time_value, num_total_nodes, res_m, 3, matr_form,
          resid_vector_sens, x_sens, x_sens_p, scale, ams, first_linear_solver_call,
          Norm_below_tolerance, Rate_above_tolerance, pp_fluxes_sens[i]->vector_id,
          pp_fluxes_sens[i]->sens_type, pp_fluxes_sens[i]->sens_id, pp_fluxes_sens[i]->sens_flt,
          pp_fluxes_sens[i]->sens_flt2, &mf_resolve, ncod, bc, &smallpiv, &singpiv, &iautopiv,
          &iscale, &scaling_max, &h_elem_avg, &U_norm, UMF_system_id, sens_caller);
      sens_vec_ct++;
    }
  }

  /*
   *
   *    DATA SENSITIVITIES
   *
   */

  for (i = 0; i < nn_post_data_sens; i++) {
    if (pp_data_sens[i]->vector_id > sens_vec_ct) {
      /*
       *
       *	GET SENSITIVITY PARAMETER VALUE
       *
       */
      retrieve_parameterS(&param_val, x, xdot, pp_data_sens[i]->sens_type, pp_data_sens[i]->sens_id,
                          pp_data_sens[i]->sens_flt, pp_data_sens[i]->sens_flt2, cx, exo, dpi);

      strcat(sens_caller, "Data Sensitivity");

      err = soln_sens(param_val, x, xdot, delta_s, exo, dpi, cx, res_p, numProcUnknowns, ija, a,
                      x_old, x_older, xdot_old, x_update, delta_t, theta, time_value,
                      num_total_nodes, res_m, 3, matr_form, resid_vector_sens, x_sens, x_sens_p,
                      scale, ams, first_linear_solver_call, Norm_below_tolerance,
                      Rate_above_tolerance, pp_data_sens[i]->vector_id, pp_data_sens[i]->sens_type,
                      pp_data_sens[i]->sens_id, pp_data_sens[i]->sens_flt,
                      pp_data_sens[i]->sens_flt2, &mf_resolve, /* frontal solver variables */
                      ncod, bc, &smallpiv, &singpiv, &iautopiv, &iscale, &scaling_max, &h_elem_avg,
                      &U_norm, UMF_system_id, sens_caller);

      sens_vec_ct++;
    }
  }
  /*
   *        OUTPUT DATA TO FILES
   */

  /*  only if converged */

  if (*converged) {

    for (i = 0; i < nn_post_data; i++) {
      err = ns_data_print(pp_data[i], x, exo, time_value, delta_t);
    }
    /*
     *       OUTPUT DATA SENSITIVITY TO FILES
     */

    for (i = 0; i < nn_post_data_sens; i++) {
      err = ns_data_sens_print(pp_data_sens[i], x, x_sens_p, time_value);
    }

  } /* if converged */

  /**

   COMPUTE SOLUTION SENSITIVITY TO THE CONTINUATION PARAMETER

   **/

  /*

  CHECK IF CONTINUATION PARAMETER IS SAME AS SOLUTION
  SENSITIVITY AND THUS SAVE A CALL
  SAME FOR HUNTING SENSITIVITIES

  */

  /*  only if converged */

  if (Continuation > 0 && *converged) {

    switch (Continuation) {
    case ALC_FIRST:

      if (cont->sensvec_id != -1) {
        for (i = 0; i < numProcUnknowns; i++) {
          x_sens[i] = x_sens_p[cont->sensvec_id][i];
        }
      } else {

        strcat(sens_caller, "First Order Continuation");

        err = soln_sens(lambda,  /* parameter */
                        x,       /* soln vector */
                        xdot,    /* dxdt */
                        delta_s, /* step */
                        exo, dpi, cx, res_p, numProcUnknowns, ija, a, x_old, x_older, xdot_old,
                        x_update, delta_t, theta, time_value, num_total_nodes, res_m, 3, matr_form,
                        resid_vector_sens, x_sens, x_sens_p, scale, ams, first_linear_solver_call,
                        Norm_below_tolerance, Rate_above_tolerance, -1, 1, 0, 0, 0,
                        &mf_resolve, /* frontal solver variables */
                        ncod, bc, &smallpiv, &singpiv, &iautopiv, &iscale, &scaling_max,
                        &h_elem_avg, &U_norm, UMF_system_id, sens_caller);
      }
      break;
    case HUN_FIRST:

      strcat(sens_caller, "First Order Hunting Continuation ");

      err = soln_sens(lambda,  /* parameter */
                      x,       /* soln vector */
                      xdot,    /* dxdt */
                      delta_s, /* step */
                      exo, dpi, cx, res_p, numProcUnknowns, ija, a, x_old, x_older, xdot_old,
                      x_update, delta_t, theta, time_value, num_total_nodes, res_m, 3, matr_form,
                      resid_vector_sens, x_sens, x_sens_p, scale, ams, first_linear_solver_call,
                      Norm_below_tolerance, Rate_above_tolerance, -1, 1, 0, 0, 0,
                      &mf_resolve, /* frontal solver variables */
                      ncod, bc, &smallpiv, &singpiv, &iautopiv, &iscale, &scaling_max, &h_elem_avg,
                      &U_norm, UMF_system_id, sens_caller);
      break;
    }

    if (Linear_Solver == AZTEC)
      total_ls_its += err;

  } /* if (Continuation > 0) */

  if (Linear_Solver == AZTEC)
    ams->status[AZ_its] = total_ls_its;
  LOCA_UMF_ID = UMF_system_id;

  /*
   * If using LOCA, there may be another resolve after exiting the
   * nonlinear solver, so defer restoring external matrix rows
   * until the next linear solver call.
   */

free_and_clear:
  if (Num_Proc > 1 && strcmp(Matrix_Format, "msr") == 0) {
    if (Continuation != LOCA && dofs_hidden) {
      show_external(num_universe_dofs[pg->imtrx],
                    (num_universe_dofs[pg->imtrx] - num_external_dofs[pg->imtrx]), ija, ija_save,
                    a);
      dofs_hidden = FALSE;
    }
  }

  safe_free((void *)delta_x);
  safe_free((void *)res_p);
  safe_free((void *)res_m);

  if (nAC > 0) {
    safe_free((void *)gAC);
    safe_free((void *)hAC);
    safe_free((void *)tAC);
    safe_free((void *)yAC);
    safe_free((void *)tempAC);

    Dmatrix_death(bAC, nAC, numProcUnknowns);
    Dmatrix_death(cAC, nAC, numProcUnknowns);
    Dmatrix_death(dAC, nAC, nAC);
    Dmatrix_death(sAC, nAC, nAC);
    Dmatrix_death(wAC, nAC, numProcUnknowns);

    Ivector_death(ordAC, nAC);

    /* Restore arc length equation AC at end of first step. */
    if (alc_with_acs && con->private_info.step_num == 0)
      nAC++;
  }

  return (return_value);
} /* end of solve_nonlinear_problem */
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/* L2_norm -- compute L2 norm of a vector scattered across multiple procs */

double L2_norm(double *vector, int nloc) {
  int i;
  double local_value;
#ifdef PARALLEL
  double global_value;
#endif /* PARALLEL */

  local_value = 0.;
  for (i = 0; i < nloc; i++, vector++) {
    local_value += (*vector) * (*vector);
  }

#ifdef PARALLEL
  MPI_Allreduce(&local_value, &global_value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  local_value = global_value;
#endif /* PARALLEL */

  return (sqrt(local_value));
}
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

double L2_norm_diff(double *v1, double *v2, const int nloc) {
  int i;
  double local_value = 0.0;
#ifdef PARALLEL
  double global_value;
#endif /* PARALLEL */

  for (i = 0; i < nloc; i++, v1++, v2++) {
    local_value += (*v1 - *v2) * (*v1 - *v2);
  }
#ifdef PARALLEL
  MPI_Allreduce(&local_value, &global_value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  local_value = global_value;
#endif /* PARALLEL */

  return (sqrt(local_value));
}
/* L2_norm -- compute L2 norm of a vector scattered across multiple procs */
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/* L2_norm -- relative */

double L2_norm_r(double *vector, double *vecscal, int nloc) {
  int i;
  double local_value;
#ifdef PARALLEL
  double global_value;
#endif /* PARALLEL */

  local_value = 0.;
  for (i = 0; i < nloc; i++, vector++, vecscal++) {
    local_value += (*vector) * (*vector) / (1 + (*vecscal) * (*vecscal));
  }

#ifdef PARALLEL
  MPI_Allreduce(&local_value, &global_value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  local_value = global_value;
#endif /* PARALLEL */

  return (sqrt(local_value));
}
/* L2_norm -- relative */

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/* L1_norm -- compute L1 norm of a vector scattered across multiple procs */

double L1_norm(double *vector, int nloc) {
  int i;
  double local_value;
#ifdef PARALLEL
  double global_value;
#endif /* PARALLEL */

  local_value = 0.;
  for (i = 0; i < nloc; i++, vector++) {
    local_value += fabs(*vector);
  }
#ifdef PARALLEL
  MPI_Allreduce(&local_value, &global_value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  local_value = global_value;
#endif /* PARALLEL */

  return (local_value);
}
/* L1_norm -- compute L1 norm of a vector scattered across multiple procs */

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
double L1_norm_r(double *vector, double *vecscal, int nloc) {
  int i;
  double local_value;
#ifdef PARALLEL
  double global_value;
#endif /* PARALLEL */

  local_value = 0.;
  for (i = 0; i < nloc; i++, vector++, vecscal++) {
    local_value += fabs(*vector) / (1 + fabs(*vecscal));
  }
#ifdef PARALLEL
  MPI_Allreduce(&local_value, &global_value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  local_value = global_value;
#endif /* PARALLEL */

  return (local_value);
}
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

double Loo_norm(double *vector, int nloc, int *num_unk, char *dofname_x)

/*******************************************************************
 *  Loo_norm -- compute an Loo norm of a vector scattered
 *              across multiple procs
 *
 * Input
 * ----------
 *  vector[]   ->  vector whose norm is to be calculated
 *  nloc       ->  local length of vector[] and vecscal[]
 *  dofname_r  ->  previously allocated space -> at least 80 chars
 *
 * Output
 * ----------
 *  return  -> value of the Loo norm
 *  num_unk -> position in the solution vector of largest component
 *            ( if not on this processor, return -1)
 *  dofname_x->  description of the current unknown
 *            ( if not on this processor, return null string )
 *******************************************************************/
{
  int i;
  double local_value = -1.0;
#ifdef PARALLEL
  struct {
    double d;
    int i;
  } in, out;
#endif /* PARALLEL */
  for (i = 0; i < nloc; i++, vector++) {
    if (fabs(*vector) > local_value) {
      local_value = fabs(*vector);
      *num_unk = i;
    }
  }
#ifdef PARALLEL
  in.d = local_value;
  in.i = ProcID;
  MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
  local_value = out.d;
  if (out.i == ProcID) {
    dofname40(*num_unk, dofname_x);
  } else {
    *num_unk = -1;
    if (dofname_x != NULL)
      dofname_x[0] = '\0';
  }
#else  /* PARALLEL */
  dofname40(*num_unk, dofname_x);
#endif /* PARALLEL */
  return (local_value);
}
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

double Loo_norm_r(double *vector, double *vecscal, int nloc, int *num_unk, char *dofname_r)

/*******************************************************************
 *  Loo_norm_r -- compute a scaled Loo norm of a vector scattered
 *                across multiple procs
 *
 * Input
 * ----------
 *  vector[]   ->  vector whose norm is to be calculated
 *  vecscal[]  ->  scaled factor for each component of vector[]
 *  nloc       ->  local length of vector[] and vecscal[]
 *  dofname_r  ->  previously allocated space -> at least 80 chars
 *
 * Output
 * ----------
 *  return  -> value of the scaled Loo norm
 *  num_unk -> position in the solution vector of largest component
 *            ( if not on this processor, return -1)
 *  dofname_x->  description of the current unknown
 *            ( if not on this processor, return null string )
 *******************************************************************/
{
  int i;
  double local_value = -1.0E200, quotient;
#ifdef PARALLEL
  struct {
    double d;
    int i;
  } in, out;
#endif /* PARALLEL */
  for (i = 0; i < nloc; i++, vector++, vecscal++) {
    quotient = fabs(*vector) / (1.0 + fabs(*vecscal));
    if (quotient > local_value) {
      local_value = quotient;
      *num_unk = i;
    }
  }
#ifdef PARALLEL
  in.d = local_value;
  in.i = ProcID;
  MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
  local_value = out.d;
  if (out.i == ProcID) {
    dofname40(*num_unk, dofname_r);
  } else {
    *num_unk = -1;
    if (dofname_r != NULL)
      dofname_r[0] = '\0';
  }
#else  /* PARALLEL */
  dofname40(*num_unk, dofname_r);
#endif /* PARALLEL */
  return (local_value);
}
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/* Loo_norm -- relative */

/* L2_norm_1p  -- compute the L2 norm of a vector on a single processor */
double L2_norm_1p(double *vector, int nloc) {
  int i;
  double local_value;

  local_value = 0.;
  for (i = 0; i < nloc; i++, vector++) {
    local_value += (*vector) * (*vector);
  }

  return (sqrt(local_value));
}

/* L2_norm_diff_1p -- compute the change in L2 norm of a vector
                        on a single processor */
double L2_norm_diff_1p(double *v1, double *v2, const int nloc) {
  int i;
  double local_value = 0.0;

  for (i = 0; i < nloc; i++, v1++, v2++) {
    local_value += (*v1 - *v2) * (*v1 - *v2);
  }

  return (sqrt(local_value));
}
/* end L2_norm_diff_1p */

/* L2_norm_r_1p -- compute the L2 norm of a scaled vector on single proc */

double L2_norm_r_1p(double *vector, double *vecscal, int nloc) {
  int i;
  double local_value;

  local_value = 0.;
  for (i = 0; i < nloc; i++, vector++, vecscal++) {
    local_value += (*vector) * (*vector) / (1 + (*vecscal) * (*vecscal));
  }

  return (sqrt(local_value));
}
/* end L2_norm_1p -- relative */

/* L1_norm_1p -- compute L1 norm of a vector on a single proc */

double L1_norm_1p(double *vector, int nloc) {
  int i;
  double local_value;

  local_value = 0.;
  for (i = 0; i < nloc; i++, vector++) {
    local_value += fabs(*vector);
  }

  return (local_value);
}
/* end L1_norm_1p */

/* L1_norm_r_1p -- compute L1 norm of a scaled vector on a single proc */
double L1_norm_r_1p(double *vector, double *vecscal, int nloc) {
  int i;
  double local_value;

  local_value = 0.;
  for (i = 0; i < nloc; i++, vector++, vecscal++) {
    local_value += fabs(*vector) / (1 + fabs(*vecscal));
  }

  return (local_value);
}
/* L1_norm -- relative */

/* Loo_norm_1p -- compute Loo norm of a vector on a single proc */
double Loo_norm_1p(double *vector, int nloc, int *num_unk, char *dofname_x) {
  int i;
  double local_value = -1e30;

  for (i = 0; i < nloc; i++, vector++) {
    if (fabs(*vector) > local_value) {
      local_value = fabs(*vector);
      *num_unk = i;
    }
  }

  return (local_value);
}
/* end Loo_norm_r_1p */

/* Loo_norm_r_1p -- compute Loo norm of a scaled vector on a single proc */
double Loo_norm_r_1p(double *vector, double *vecscal, int nloc, int *num_unk, char *dofname_r) {
  int i;
  double local_value;
  double quotient;

  local_value = -1.e30;

  for (i = 0; i < nloc; i++, vector++, vecscal++) {
    quotient = fabs(*vector) / (1 + fabs(*vecscal));
    if (quotient > local_value) {
      local_value = quotient;
      *num_unk = i;
    }
  }

  return (local_value);
}
/* end Loo_norm_r_1p */
/*
 * print_array() - dump out a list of array components to a procname-based file
 *
 * Example usage:
 *
 *	print_array(x, n, "x", "%g", ProcID);
 *
 *	print_array(ija, nnz, "ija", "%d", ProcID);
 *
 *
 * Created: 1997/10/28 17:26 MST pasacki@sandia.gov
 * Revised:
 */

void print_array(
    const void *array, const int length, const char *name, const Edt datatype, const int procid) {
  int i;
  char filename[MAX_FNL];
  FILE *f;

  sprintf(filename, "%s_%dof%d", name, procid + 1, Num_Proc);
  f = fopen(filename, "a");

#if 0
  fprintf(f, "\n---------------------------------------------------------\n");
#endif

  switch (datatype) {
  case type_int:
    for (i = 0; i < length; i++) {
#if 0
	  fprintf(f, "%s[%d] = %d\n", name, i, ((int *)array)[i]);
#endif
      fprintf(f, "%d %d\n", i, ((const int *)array)[i]);
    }
    break;

  case type_double:
    for (i = 0; i < length; i++) {
#if 0
	  fprintf(f, "%s[%d] = %g\n", name, i, ((double *)array)[i]);
#endif
      fprintf(f, "%d %23.16e\n", i, ((const double *)array)[i]);
    }
    break;

  default:
    GOMA_EH(GOMA_ERROR, "Cannot handle that datatype.");
    break;
  }

  fclose(f);
  return;
}

/*******************************************************************
 *     soln_sens.c
 *  computes sensitivity of the solution from the Jacobian matrix
 *	using finite differences
 *******************************************************************/

/*

    SOLVE:

    J dq/dp = - dR/dp

    J: JACOBIAN MATRIX, dR/dq
    q: VECTOR OF GOMA UNKNOWNS
    R: RESIDUAL VECTOR
    p: PARAMETER

*/

static int soln_sens(double lambda,  /*  parameter */
                     double x[],     /* soln vector */
                     double xdot[],  /*  dxdt predicted for new time */
                     double delta_s, /* step */
                     Exo_DB *exo,
                     Dpi *dpi,
                     Comm_Ex *cx,
                     double res_p[],
                     int numProcUnknowns,
                     int ija[],
                     double a[],
                     double x_old[],
                     double x_older[],
                     double xdot_old[],
                     double x_update[],
                     double delta_t,
                     double theta,
                     double time_value,
                     int num_total_nodes,
                     double res_m[],
                     int Factor_Flag,
                     int matr_form,
                     double resid_vector_sens[],
                     double x_sens[],
                     double **x_sens_p,
                     double scale[],
                     struct GomaLinearSolverData *ams, /* ptrs to
                                                        * Aztec
                                                        * linear
                                                        * systems
                                                        */
                     int first_linear_solver_call,
                     int Norm_below_tolerance,
                     int Rate_above_tolerance,
                     int vector_id,
                     int sens_type,
                     int sens_id,
                     int sens_flt,
                     int sens_flt2,
                     int *mf_resolve, /* frontal solver variables */
                     int *fsncod,
                     double *fsbc,
                     double *smallpiv,
                     double *singpiv,
                     int *iautopiv,
                     int *iscale,
                     double *scaling_max,
                     double *ptr_h_elem_avg,
                     double *ptr_U_norm,
                     int UMF_system_id,
                     char *sens_caller)

{
  double dlambda, lambda_tmp, hunt_val;

  int i, iHC, iAC;
  dbl a_start; /* mark start of assembly */
  dbl a_end;   /* mark end of assembly */
  int err;

  int linear_solver_blk;     /* count calls to AZ_solve() */
  int num_linear_solve_blks; /* one pass for now */
  int matrix_solved;         /* boolean */
  char stringer[80];         /* holding format of num linear solve itns */

  double h_elem_avg = *ptr_h_elem_avg;
  double U_norm = *ptr_U_norm;

  double time_local = 0.0;
  double time_global = 0.0;

  double fd_factor = FD_FACTOR; /*  finite difference step */

  /*
  static int first_soln_sens_linear_solver_call = 1;
  */
  /*  static int UMF_system_id = -1; */ /* We'll give soln_sens it's own UMF
                                         * system space, since it does its own
                                         * matrix_fill's. */

  exchange_dof(cx, dpi, x, pg->imtrx);

  /*
   *
   * GET SENSITIVITY OF RESIDUAL TO NATURAL PARAMETER
   *
   */
  err = 0;
  a_start = ut();

  dlambda = fd_factor * lambda;
  dlambda = (fabs(dlambda) < fd_factor ? fd_factor : dlambda);

  /*
   * GET RESIDUAL SENSITIVITY
   */

  lambda_tmp = lambda + dlambda;

  if (vector_id == -1) {
    switch (Continuation) {
    case ALC_FIRST:
      update_parameterC(0, lambda_tmp, x, xdot, NULL, delta_s, cx, exo, dpi);
      break;
    case HUN_FIRST:
      for (iHC = 0; iHC < nHC; iHC++) {
        if (abs(hunt[iHC].ramp) == 2) {
          hunt_val = hunt[iHC].BegParameterValue *
                     pow(hunt[iHC].EndParameterValue / hunt[iHC].BegParameterValue, lambda_tmp);
        } else {
          hunt_val = hunt[iHC].BegParameterValue +
                     lambda_tmp * (hunt[iHC].EndParameterValue - hunt[iHC].BegParameterValue);
        }
        update_parameterHC(iHC, hunt_val, x, xdot, NULL, 1.0, cx, exo, dpi);
      }
      break;
    }
  } else {
    update_parameterS(lambda_tmp, x, xdot, sens_type, sens_id, sens_flt, sens_flt2, cx, exo, dpi);
  }

  init_vec_value(res_p, 0.0, numProcUnknowns);
  af->Assemble_Residual = TRUE;
  af->Assemble_Jacobian = FALSE;
  af->Assemble_LSA_Jacobian_Matrix = FALSE;
  af->Assemble_LSA_Mass_Matrix = FALSE;
  if (nAC > 0) {
    for (iAC = 0; iAC < nAC; iAC++) {
      augc[iAC].evol = 0.0;
    }
  }

  err = matrix_fill_full(ams, x, res_p, x_old, x_older, xdot, xdot_old, x_update, &delta_t, &theta,
                         First_Elem_Side_BC_Array[pg->imtrx], &time_value, exo, dpi,
                         &num_total_nodes, &h_elem_avg, &U_norm, NULL);
  if (err == -1)
    return (err);

  lambda_tmp = lambda - dlambda;
  if (vector_id == -1) {
    switch (Continuation) {
    case ALC_FIRST:
      update_parameterC(0, lambda_tmp, x, xdot, NULL, delta_s, cx, exo, dpi);
      break;
    case HUN_FIRST:
      for (iHC = 0; iHC < nHC; iHC++) {
        if (abs(hunt[iHC].ramp) == 2) {
          hunt_val = hunt[iHC].BegParameterValue *
                     pow(hunt[iHC].EndParameterValue / hunt[iHC].BegParameterValue, lambda_tmp);
        } else {
          hunt_val = hunt[iHC].BegParameterValue +
                     lambda_tmp * (hunt[iHC].EndParameterValue - hunt[iHC].BegParameterValue);
        }
        update_parameterHC(iHC, hunt_val, x, xdot, NULL, 1.0, cx, exo, dpi);
      }
      break;
    }
  } else {
    update_parameterS(lambda_tmp, x, xdot, sens_type, sens_id, sens_flt, sens_flt2, cx, exo, dpi);
  }

  init_vec_value(res_m, 0.0, numProcUnknowns);
  af->Assemble_Residual = TRUE;
  af->Assemble_Jacobian = FALSE;
  af->Assemble_LSA_Jacobian_Matrix = FALSE;
  af->Assemble_LSA_Mass_Matrix = FALSE;
  if (nAC > 0) {
    for (iAC = 0; iAC < nAC; iAC++) {
      augc[iAC].evol = 0.0;
    }
  }

  err = matrix_fill_full(ams, x, res_m, x_old, x_older, xdot, xdot_old, x_update, &delta_t, &theta,
                         First_Elem_Side_BC_Array[pg->imtrx], &time_value, exo, dpi,
                         &num_total_nodes, &h_elem_avg, &U_norm, NULL);

  if (err == -1)
    return (err);

  /*
   * RESET PARAMETERS
   */

  if (vector_id == -1) {
    switch (Continuation) {
    case ALC_FIRST:
      update_parameterC(0, lambda, x, xdot, NULL, delta_s, cx, exo, dpi);
      break;
    case HUN_FIRST:
      for (iHC = 0; iHC < nHC; iHC++) {
        if (abs(hunt[iHC].ramp) == 2) {
          hunt_val = hunt[iHC].BegParameterValue *
                     pow(hunt[iHC].EndParameterValue / hunt[iHC].BegParameterValue, lambda);
        }

        else {
          hunt_val = hunt[iHC].BegParameterValue +
                     lambda * (hunt[iHC].EndParameterValue - hunt[iHC].BegParameterValue);
        }
        update_parameterHC(iHC, hunt_val, x, xdot, NULL, 1.0, cx, exo, dpi);
      }
      break;
    }
  } else {
    update_parameterS(lambda, x, xdot, sens_type, sens_id, sens_flt, sens_flt2, cx, exo, dpi);
  }

  lambda_tmp = 0.5 / dlambda;
  v2sum(numProcUnknowns, &resid_vector_sens[0], lambda_tmp, &res_p[0], -lambda_tmp, &res_m[0]);

  if (vector_id != -1) {
    for (i = 0; i < numProcUnknowns; i++) {
      x_sens_p[vector_id + 1][i] = resid_vector_sens[i];
    }
  }
  /*
   *
   * GET SENSITIVITY OF SOLUTION TO NATURAL PARAMETER
   *
   */
#ifdef MATRIX_DUMP_SENS
  if (strcmp(Matrix_Format, "msr") == 0) {
    matrix_dump_msr(ams, exo, dpi, x);
  } else {
    matrix_dump_vbr(ams, exo, dpi, x);
  }
#endif

  vector_scaling(NumUnknowns[pg->imtrx], resid_vector_sens, scale);

  switch (Linear_Solver) {
  case UMFPACK2:
  case UMFPACK2F:
    /*  */
    matr_form = 1;

    /* Factor_Flag is always 3 when entering soln_sens, as far as I
     * can tell... */

    /* MMH: I don't believe that a refactorization is EVER intended
     * here.  However, that was the option selected when
     * Linear_Solver is UMFPACK2F.  I'm guessing that although this
     * isn't, in fact, what you want to happen, it is what was
     * previously done before I consolidated UMFPACK2 and
     * UMFPACK2F. */
    if (Linear_Solver == UMFPACK2F)
      Factor_Flag = 0;

    if (first_linear_solver_call)
      GOMA_EH(GOMA_ERROR, "Solving for AC's BEFORE a regular solve");

      /* MMH: I believe that this system will always be the same
       * structure/system as the one used for the regular solve.  This
       * was the behavior before I consolidated the UMFPACK and
       * UMFPACK2F. */

#ifdef DEBUG_SL_UMF
    printf("%s: entering SL_UMF for soln_sens solve\n", y0);
    fflush(stdout);
#endif
    UMF_system_id =
        SL_UMF(UMF_system_id,
               /*
               &first_soln_sens_linear_solver_call,
               */
               &first_linear_solver_call, &Factor_Flag, &matr_form, &NumUnknowns[pg->imtrx],
               &NZeros, &ija[0], &ija[0], &a[0], &resid_vector_sens[0], &x_sens[0]);
#ifdef DEBUG_SL_UMF
    printf("%s: returning from SL_UMF\n", y0);
    fflush(stdout);
#endif
    /*
    first_soln_sens_linear_solver_call = 0;
    */

    strcpy(stringer, " 1 ");
    break;

  case SPARSE13a:
    dcopy1(NumUnknowns[pg->imtrx], resid_vector_sens, x_sens);
    lu(NumUnknowns[pg->imtrx], NumExtUnknowns[pg->imtrx], NZeros, a, ija, x_sens, 3);
    first_linear_solver_call = FALSE;
    /*
     * Note that sl_lu has static variables to keep track of
     * first call or not.
     */
    strcpy(stringer, " 1 ");
    break;

  case AMESOS:
    if (strcmp(Matrix_Format, "msr") == 0) {
      amesos_solve_msr(Amesos_Package, ams, x_sens, resid_vector_sens, 0, pg->imtrx);
    } else if (strcmp(Matrix_Format, "epetra") == 0) {
      amesos_solve_epetra(Amesos_Package, ams, x_sens, resid_vector_sens, pg->imtrx);
    } else {
      GOMA_EH(GOMA_ERROR,
              " Sorry, only MSR and Epetra matrix formats are currently supported with the Amesos "
              "solver suite\n");
    }
    strcpy(stringer, " 1 ");
    break;
    break;

  case AZTEC:
    /*
     * Initialization is now performed up in
     * solve_problem() in rf_solve.c, not here in the
     * midst of the Newton iteration.
     */
    ams->options[AZ_pre_calc] = AZ_calc;

    linear_solver_blk = 0;     /* count calls to AZ_solve() */
    num_linear_solve_blks = 1; /* upper limit to AZ_solve() calls */
    matrix_solved = FALSE;
    while ((!matrix_solved) && (linear_solver_blk < num_linear_solve_blks)) {
      /*
       * Someday the user may want to do fancy heuristics based
       * on all kinds of cost functions, artificial intelligence
       * neural networks, etc.
       *
       * For the linear system "Ax=b", we have
       *    A -- indx, bindx(ija), rpntr, cpntr, bpntr, val(a)
       *    x -- x_sens, newton correction vector
       *    b -- resid_vector_sens, newton residual equation vector
       */
      AZ_solve(x_sens, resid_vector_sens, ams->options, ams->params, ams->indx, ams->bindx,
               ams->rpntr, ams->cpntr, ams->bpntr, ams->val, ams->data_org, ams->status,
               ams->proc_config);

      if (Debug_Flag > 0) {
        dump_aztec_status(ams->status);
      }

      aztec_stringer((int)ams->status[AZ_why], ams->status[AZ_its], &stringer[0]);

      matrix_solved = (ams->status[AZ_why] == AZ_normal);
      linear_solver_blk++;
    }

    break;
  case AZTECOO:
    if (strcmp(Matrix_Format, "epetra") == 0) {
      aztecoo_solve_epetra(ams, x_sens, resid_vector_sens);
      aztec_stringer((int)ams->status[AZ_why], ams->status[AZ_its], &stringer[0]);
      matrix_solved = (ams->status[AZ_why] == AZ_normal);
    } else {
      GOMA_EH(GOMA_ERROR,
              "Sorry, only Epetra matrix formats are currently supported with the AztecOO solver "
              "suite\n");
    }
    break;

  case STRATIMIKOS:
    if (strcmp(Matrix_Format, "epetra") == 0) {
      int iterations;
      int err = stratimikos_solve(ams, x_sens, resid_vector_sens, &iterations, Stratimikos_File,
                                  pg->imtrx);
      GOMA_EH(err, "Error in stratimikos solve");
      if (iterations == -1) {
        strcpy(stringer, "err");
      } else {
        aztec_stringer(AZ_normal, iterations, &stringer[0]);
      }
    } else {
      GOMA_EH(GOMA_ERROR,
              "Sorry, only Epetra matrix formats are currently supported with the Stratimikos "
              "interface\n");
    }
    break;
  case MA28:
    /*
     * sl_ma28 keeps internal static variables to determine whether
     * it is the first call or not.
     */
#ifdef HARWELL
    err = cmsr_ma28(NumUnknowns[pg->imtrx], NZeros, a, ija, x_sens, resid_vector_sens);
#endif
#ifndef HARWELL
    GOMA_EH(GOMA_ERROR, "That linear solver package is not implemented.");
#endif
    strcpy(stringer, " 1 ");
    break;
  default:
    GOMA_EH(GOMA_ERROR, "That linear solver package is not implemented.");
    break;
  }

  /*
   * GET RIGHT SIGN FOR dx/dlambda;
   * ABOVE WE SOLVED:  J dx/d* = dR/d*
   */

  vchange_sign(numProcUnknowns, &x_sens[0]);

  exchange_dof(cx, dpi, x_sens, pg->imtrx);

  a_end = ut();

  if (vector_id != -1) {
    for (i = 0; i < numProcUnknowns; i++) {
      x_sens_p[vector_id][i] = x_sens[i];
    }
  }

  /*
   * DONE
   */

  time_local = a_end - a_start;
#ifdef PARALLEL
  MPI_Allreduce(&time_local, &time_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  time_local = time_global;
#endif
  fflush(stdout);
  sens_caller[strlen(sens_caller)] = '\0';
  DPRINTF(stdout, "\n\n%s resolve time:  %7.1e\n", sens_caller, time_local);
  sens_caller[0] = '\0';

  return (err);
} /*   end of routine soln_sens()   */
/* end of file mm_sol_nonlinear.c */
