/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2014 Sandia Corporation.                                  *
*                                                                         *
* Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,  *
* the U.S. Government retains certain rights in this software.            *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
\************************************************************************/
 
/*
 *$Id: ac_loca_interface.c,v 5.8 2008-12-19 22:54:25 rbsecor Exp $
 */

#ifdef USE_RCSID
static char rcsid[] =
"$Id: ac_loca_interface.c,v 5.8 2008-12-19 22:54:25 rbsecor Exp $";
#endif

/* DOCUMENTATION FOR CONTINUATION LIBRARY:    10/26/1999
 *    Andrew Salinger  Org. 9221, MS-1111, 845-3523
 *
 * The continuation library currently consists of 9 files. 
 * (1) [ac_loca_interface.c]  The Goma version of the interface
 *     to the continuation library. This file (and only this file)
 *     must be extensively modified for every new code that
 *     uses the library. It consists of solve_continuation, the
 *     top level routine called by the application code, and ~20
 *     wrapper routines, which must be filled for your specific
 *     code. The most noteworthy of these include
 *     nonlinear_solver_conwrap, linear_solver_conwrap,
 *     matrix_residual_fill_conwrap, and assign_parameter_conwrap.
 *     Each wrapper has its own comments.
 * (2) [loca_lib.c]  This is the stepping algorithm for zeroth-order,
 *     first-order, and arc-length continuation. The turning point
 *     tracking algorithm is only zeroth-order. This routine includes
 *     such things as step-size control and predictor steps. Ideally,
 *     this file does not need to be modified when linking to a new
 *     application code.
 * (3) [loca_bord.c] This file contains the bordering algorithms,
 *     currently for arc-length continuation and turning point tracking.
 *     These routines are accessed from within the Newton iteration,
 *     and make frequent calls to the wrapper routines (matrix fill
 *     and solve). Ideally, this file does not need to be modified when
 *     linking to a new application code.
 * (4) [loca_util.c] This file contains utility functions used by
 *     LOCA algorithms such as vector operations (allocation, copying,
 *     dot products) and eigenvector sorting.
 * (5) [loca_const.h] This header file includes definitions of the
 *     continuation structures, define statements for flags used
 *     by the continuation library, and prototypes for the wrapper
 *     routines. Ideally, this file does not need to be modified when
 *     linking to a new application code.
 * (6) [loca_util_const.h] This header file contains prototypes for the
 *     functions in loca_util.c.
 * (7) [loca_eigenvalue.c] This file is the LOCA interface to the
 *     ARPACK/PARPACK eigensolver package.
 * (8) [loca_eigen_cayley.F] This file contains FORTRAN subroutines
 *     which are called while running ARPACK.
 * (9) [loca_eigen_c2f.F] This file contains routines which allow
 *     FORTRAN subroutines to be called from C functions.
 *
 * How to interface to this library:
 * (0) Have a steady-state code that uses Newton's method
 * (1) Call solve_continuation from your code, in the place
 *     where you normally call the steady-state or transient
 *     drivers.
 * (2) In your current nonlinear solver, add "(void *) con_ptr"
 *     to the argument list. Pass NULL in that place for all
 *     calls to the nonlinear solver besides the one from
 *     nonlinear_solver_conwrap below. In your Newton loop,
 *     after the matrix equation has been solved for the
 *     update vector, delta_x, but before the solution vector, x,
 *     has been updated by delta_x, put in the following call:
 *      if (con_ptr != NULL) continuation_converged =
 *         continuation_hook(x, delta_x, con_ptr, Reltol, Abstol);
 *     (Reltol=1.0e-3, Abstol=1.0e-8 are good defaults.)
 *     When checking convergence of your Newton's method, also
 *     check that (continuation_converged == TRUE).
 * (3) Change the contents of all routines in this file. In
 *     solve_continuation, many fields of the con and set_con structures
 *     must be set. Follow the template and the comments.
 *     Also, the passdown structure can be defined and set with 
 *     any information that needs to be passed from the 
 *     solve_continuation routine to the wrapper routines below
 *     that isn't needed by the continuation library. All
 *     of the wrapper routines (all routines in this file besides
 *     solve_continuation, which all have names ending in _conwrap)
 *     must be filled in with the corresponding call to your
 *     application code. Each has its own comments.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>

/* Put include statements for your code here. */

#include "ac_conti.h"
#include "ac_hunt.h"
#include "ac_particles.h"
#include "ac_stability.h"
#include "ac_stability_util.h"
#include "ac_update_parameter.h"
#include "bc_colloc.h"
#include "bc_contact.h"
#include "bc_curve.h"
#include "bc_dirich.h"
#include "bc_integ.h"
#include "bc_rotate.h"
#include "bc_special.h"
#include "bc_surfacedomain.h"
#include "dp_comm.h"
#include "dp_map_comm_vec.h"
#include "dp_types.h"
#include "dp_utils.h"
#include "dp_vif.h"
#include "dpi.h"
#include "el_elm.h"
#include "el_elm_info.h"
#include "el_geom.h"
#include "el_quality.h"
#include "exo_conn.h"
#include "exo_struct.h"
#include "loca_const.h"
#include "loca_util_const.h"
#include "md_timer.h"
#include "mm_as.h"
#include "mm_as_alloc.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_augc_util.h"
#include "mm_bc.h"
#include "mm_chemkin.h"
#include "mm_dil_viscosity.h"
#include "mm_eh.h"
#include "mm_elem_block_structs.h"
#include "mm_fill.h"
#include "mm_fill_aux.h"
#include "mm_fill_fill.h"
#include "mm_fill_jac.h"
#include "mm_fill_ls.h"
#include "mm_fill_porous.h"
#include "mm_fill_potential.h"
#include "mm_fill_pthings.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_rs.h"
#include "mm_fill_shell.h"
#include "mm_fill_solid.h"
#include "mm_fill_species.h"
#include "mm_fill_stress.h"
#include "mm_fill_terms.h"
#include "mm_fill_util.h"
#include "mm_flux.h"
#include "mm_input.h"
#include "mm_interface.h"
#include "mm_more_utils.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "mm_ns_bc.h"
#include "mm_numjac.h"
#include "mm_post_def.h"
#include "mm_post_proc.h"
#include "mm_prob_def.h"
#include "mm_qtensor_model.h"
#include "mm_shell_bc.h"
#include "mm_shell_util.h"
#include "mm_sol_nonlinear.h"
#include "mm_species.h"
#include "mm_std_models.h"
#include "mm_unknown_map.h"
#include "mm_viscosity.h"
#include "rd_dpi.h"
#include "rd_exo.h"
#include "rd_mesh.h"
#include "rf_allo.h"
#include "rf_bc.h"
#include "rf_bc_const.h"
#include "rf_element_storage_const.h"
#include "rf_element_storage_struct.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_fill_const.h"
#include "rf_io.h"
#include "rf_io_const.h"
#include "rf_io_structs.h"
#include "rf_masks.h"
#include "rf_mp.h"
#include "rf_node_const.h"
#include "rf_pre_proc.h"
#include "rf_shape.h"
#include "rf_solve.h"
#include "rf_solver.h"
#include "rf_solver_const.h"
#include "rf_util.h"
#include "rf_vars_const.h"
#include "sl_amesos_interface.h"
#include "sl_aux.h"
#include "sl_auxutil.h"
#include "sl_eggroll.h"
#include "sl_lu.h"
#include "sl_matrix_util.h"
#include "sl_umf.h"
#include "sl_util.h"
#include "std.h"
#include "user_ac.h"
#include "user_bc.h"
#include "user_mp.h"
#include "user_mp_gen.h"
#include "user_post.h"
#include "user_pre.h"
#include "wr_dpi.h"
#include "wr_exo.h"
#include "wr_side_data.h"
#include "wr_soln.h"

#ifdef KOMPLEX
#include "azk_komplex.h"
#endif

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#ifdef HAVE_FRONT
extern int mf_setup
(int *,                   /* nelem_glob */
       int *,                   /* neqn_glob */
       int *,                   /* mxdofel */
       int *,                   /* nfullsum */
       int *,                   /* symflag */
       int *,                   /* nell_order */
       int *,                   /* el_proc_assign */
       int *,                   /* level */
       int *,                   /* nopdof */
       int *,                   /* loc_dof */
       int *,                   /* constraint */
       const char *,            /* cname */
       int *);                 /* allocated */
#endif

/* Define passdown structure: this structure is global to this file, and
 * provides a way to pass variables from the top solve_continuation routine
 * down to the various wrapper routines below, without passing through the
 * continuation library. For instance, the Jacobian matrix is never seen
 * in the continuation routines, but is needed by several of the wrapper
 * routines in this file. The passdown structure is a way make the
 * matrix global to this file.
 */

struct passdown_struct {
  int num_mat_fills;
  int num_res_fills;
  int num_linear_its;
  int num_eigen_its;
  int nvd;
  int sv_flag;
  int *sv_index;
  double *sv_count;

 /* Include arguments to GOMA solver routine */

  struct Action_Flags *af;	/* GOMA action flags */
  struct elem_side_bc_struct **First_Elem_Side_BC_Array;
  struct Aztec_Linear_Solver_System *ams; /* ptrs to Aztec linear systems */
  double *x;     /* soln vector on this proc */
  double delta_t; /* time step size */
  double theta;   /* parameter to vary time integration from 
                   *   explicit (theta = 1) to implicit (theta = 0) */
  double *x_old; /* soln vector @ previous time */
  double *x_older; /* soln vector @ previous, previous time */
  double *x_oldest; /* oldest soln vector saved */
  double *xdot;  /* dxdt predicted for new time */
  double *xdot_old; /* dxdt for previous time */
  double *resid_vector;
  double *x_update; 
  double *scale;  /*Scale factor held for modified newton resolves */
  int *converged; /* whether the Newton iteration has converged (out) */
  int *nprint;    /* counter for time step number */
  int    tev;     /* total number elem variables to
                                             * output to EXODUS II file */
  int tev_post;   /* extra element post processing results */
  RESULTS_DESCRIPTION_STRUCT *rd; /* details about post proc vars */
  AZ_MATRIX *amat;
  AZ_MATRIX *mmat;
  int  *gindex;
  int  *gsize;
  double *gvec;
  double ***gvec_elem;
  double time_value;
  Exo_DB *exo;
  Dpi *dpi;
  Comm_Ex *cx;
  int nt;
  int *path_step_reform;
  int is_steady_state;
  double *x_AC;   /* updating evp quantities */
  double *x_AC_dot;   
  double lambda;
  double *resid_vector_sens;
  double *x_sens;
  double *x_sens_temp;
  double **x_sens_p;  /*  solution sensitivities */
  int     *proc_config;
  int     *options;
  double  *status;
  double  *mass_matrix;  /* mass matrix */
  double  *shifted_matrix;  /* shifted matrix for ARPACK */
  int method;
  int do_3D_of_2D;
  int LSA_flag;
  int last_step;
  FILE *file;	   /* Soln_OutFile */
} passdown;

  static void print_final(double conparam, int step_num, int num_newt_its,
               		  int num_modnewt_its, int total_linear_its);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int do_loca (Comm_Ex *cx,  /* array of communications structures */
             Exo_DB  *exo, /* ptr to the finite element mesh database */
             Dpi     *dpi) /* distributed processing information */

/* Interface routine to the continuation library.
 *
 * Input:
 *    con_par_ptr   Pointer to the initial value of the continuation
 *                  parameter.
 *    proc_config, options, status:  Aztec arrays.
 *
 * Output:
 *    nstep         Pointer to integer for number of steps taken, negative
 *                  if the continuation run failed.
 *
 * Return Value:
 *    nstep is also the return value -- go figure.
 */
{

  /* Define continuation problem arrays */
  int    *ija=NULL;             /* column pointer array			*/
  int    *ija_attic=NULL;       /* column pointer array storage		*/
  int    *node_to_fill=NULL;
  double *a=NULL;               /* nonzero array			*/
  double *a_old=NULL;           /* nonzero array			*/
  double *x=NULL;               /* solution vector			*/
  double *x_old=NULL;           /* old solution vector			*/
  double *x_older=NULL;         /* older solution vector		*/
  double *x_oldest=NULL;        /* oldest solution vector saved		*/
  double *xdot=NULL;            /* current path derivative of soln	*/
  double *xdot_old=NULL;        /* previous solution derivative		*/
  double *x_update=NULL;        /* Newton update to solution		*/
  double *x_AC=NULL;            /* AC extra unknown array		*/
  double *x_AC_old=NULL;        /* AC old extra unknown array		*/
  double *x_AC_dot=NULL;       
  double *x_sens=NULL;          /* solution sensitivity			*/
  double *x_sens_temp=NULL;     /* temp. solution sensitivity		*/
  double **x_sens_p=NULL;       /* Solution sensitivity to parameters	*/
  double *resid_vector=NULL;    /* Newton residual of solution		*/
  double *resid_vector_sens=NULL; /* Newton residual sensitivity	*/
  double *scale=NULL;           /* Scale vector for solution		*/
  double *gvec=NULL;
  double ***gvec_elem=NULL;
  double timeValueRead = 0.0;

  /* Define continuation problem structures */
  RESULTS_DESCRIPTION_STRUCT *rd=NULL;
#ifdef COUPLED_FILL
  struct Aztec_Linear_Solver_System *ams[NUM_ALSS]={NULL}; 
#else /* COUPLED_FILL */
  struct Aztec_Linear_Solver_System *ams[NUM_ALSS]={NULL, NULL}; 
#endif /* COUPLED_FILL */
                 /* sl_util_structs.h */

  /* Define continuation problem files */
  FILE          *cl_aux=NULL;  /* *file=NULL; */

  /* Define local variables */
  int iAC;
  int eb_indx, ev_indx;
  int num_pvector=0;
  int n; /* ni; */
  int nstep;
  int nprint;
  int i, err, error;
/*  int inewton;  */
  int converged;
  int *step_reform = NULL;
  int numProcUnknowns;
  int *gindex = NULL;
  int gsize;
  int *p_gsize;
  int tev, tev_post;
  int tnv, tnv_post;
#ifdef HAVE_FRONT
  int max_unk_elem, one, three; /* used only for HAVE_FRONT */
                                /* but must be declared anyway */
#endif
  unsigned int matrix_systems_mask;
  double *con_par_ptr;
  double lambda, delta_s;
  /* double err_dbl;  */
 /* double evol_local=0.0; */
#ifdef PARALLEL
#endif
  
  /* Define LOCA continuation structures */
  struct con_struct con;

  char *yo = "do_loca";

  /********************** First Executable Statment **********************/

  /* As of 10/26/2001, LOCA works in parallel for all continuation algs. */
#ifdef DEBUG
  DPRINTF(stderr, "%s() begins...\n", yo);
#endif

  /* Perform initialization (as in ac_conti.c) */
  p_gsize = &gsize;
  nstep = 0;
  nprint = 0;

  /*
   * Some preliminaries to help setup EXODUS II database output.
   * tnv_post is calculated in load_nodal_tkn
   * tev_post is calculated in load_elem_tkn
   */
#ifdef DEBUG
  DPRINTF(stderr, "cnt_nodal_vars() begins...\n");
#endif
  tnv = cnt_nodal_vars();
  tev = cnt_elem_vars();
  
#ifdef DEBUG
  DPRINTF(stderr, "Found %d total primitive nodal variables to output.\n", tnv);
  DPRINTF(stderr, "Found %d total primitive elem variables to output.\n", tev);
#endif
  
  if (tnv < 0)
    {
      DPRINTF(stderr, "%s:\tbad tnv.\n", yo);
      EH(-1, "\t");
    }
  
  rd = (struct Results_Description *) 
    smalloc(sizeof(struct Results_Description));

  if (rd == NULL) 
    EH(-1, "Could not grab Results Description.");

  (void) memset((void *) rd, 0, sizeof(struct Results_Description));
  
  rd->nev = 0;                  /* number element variables in results */
  rd->ngv = 0;                  /* number global variables in results */
  rd->nhv = 0;                  /* number history variables in results */

  rd->ngv = 5;                  /* number global variables in results 
                                   see load_global_var_info for names*/
  error = load_global_var_info(rd, 0, "CONV");
  error = load_global_var_info(rd, 1, "NEWT_IT");
  error = load_global_var_info(rd, 2, "MAX_IT");
  error = load_global_var_info(rd, 3, "CONVRATE");
  error = load_global_var_info(rd, 4, "MESH_VOLUME");

  /* load nodal types, kinds, names */
  error = load_nodal_tkn(rd, &tnv, &tnv_post); 
  
  if (error)
    {
      DPRINTF(stderr, "%s:  problem with load_nodal_tkn()\n", yo);
      EH(-1,"\t");
    }

  /* load elem types, names */
  error = load_elem_tkn(rd, exo, tev, &tev_post); 
  if (error)
    {
      DPRINTF(stderr, "%s:  problem with load_elem_tkn()\n", yo);
      EH(-1,"\t");
    }
#ifdef PARALLEL
  check_parallel_error("Results file error");
#endif

  /* 
   * Write out the names of the nodal variables that we will be sending to
   * the EXODUS II output file later.
   */
#ifdef DEBUG
  DPRINTF(stderr, "wr_result_prelim() starts...\n", tnv);
#endif

  gvec_elem = (double ***) smalloc ( (exo->num_elem_blocks)*sizeof(double **));
  for (i = 0; i < exo->num_elem_blocks; i++)
    gvec_elem[i] = (double **) smalloc ( (tev + tev_post)*sizeof(double *));

  wr_result_prelim_exo(rd, exo, ExoFileOut, gvec_elem);

#ifdef DEBUG
  fprintf(stderr, "P_%d: wr_result_prelim_exo() ends...\n", ProcID, tnv);
#endif

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

  numProcUnknowns = NumUnknowns[pg->imtrx] + NumExtUnknowns[pg->imtrx];

  /* allocate memory for Volume Constraint Jacobian */
  if ( nAC > 0)
    for(iAC=0;iAC<nAC;iAC++)
      augc[iAC].d_evol_dx = (double*) malloc(numProcUnknowns*sizeof(double));
  
  asdv(&resid_vector, numProcUnknowns);
  asdv(&resid_vector_sens, numProcUnknowns);
  asdv(&scale, numProcUnknowns);

  /* Allocate Aztec system structure(s) */
  for (i = 0; i < NUM_ALSS; i++) 
    {
      ams[i] = (struct Aztec_Linear_Solver_System *)
        array_alloc(1, 1, sizeof(struct Aztec_Linear_Solver_System )); 
    }

  /* Set Aztec proc_config array (for many different cases) */
#ifdef MPI
  AZ_set_proc_config( ams[0]->proc_config, MPI_COMM_WORLD );
#ifndef COUPLED_FILL
  if( Explicit_Fill ) AZ_set_proc_config( ams[1]->proc_config, MPI_COMM_WORLD );
#endif /* not COUPLED_FILL */
#else /* MPI */
  AZ_set_proc_config( ams[0]->proc_config, 0 );
#ifndef COUPLED_FILL
  if( Explicit_Fill ) AZ_set_proc_config( ams[1]->proc_config, 0 );
#endif /* not COUPLED_FILL */
#endif /* MPI */

  /* allocate space for and initialize solution arrays */
  asdv(&x,        numProcUnknowns);
  asdv(&x_old,    numProcUnknowns);
  asdv(&x_older,  numProcUnknowns);
  asdv(&x_oldest, numProcUnknowns);
  asdv(&xdot,     numProcUnknowns);
  asdv(&xdot_old, numProcUnknowns);
  asdv(&x_update, numProcUnknowns);
  asdv(&x_sens,   numProcUnknowns);
  asdv(&x_sens_temp,   numProcUnknowns);

  /* Initialize solid inertia flag */
  set_solid_inertia();
  
  /* FRIENDLY COMMAND LINE EQUIV */

  if( ProcID == 0 )
   {
      cl_aux = fopen("goma-cl.txt", "w+");

      fprintf(cl_aux, "goma -a -i input ");
      fprintf(cl_aux, "-cb %10.6e ", cont->BegParameterValue);
      fprintf(cl_aux, "-ce %10.6e ", cont->EndParameterValue);
      fprintf(cl_aux, "-cd %10.6e ", cont->Delta_s0);
      fprintf(cl_aux, "-cn %d ", cont->MaxPathSteps);
      fprintf(cl_aux, "-cm %d ", Continuation);
      fprintf(cl_aux, "-ct %d ", cont->upType);

      switch (cont->upType)
        {
        case 1:                 /* BC TYPE */
        case 3:                 /* AC TYPE */
          fprintf(cl_aux, "-c_bc %d ", cont->upBCID);
          fprintf(cl_aux, "-c_df %d ", cont->upDFID);
          break;
        case 2:                 /* MAT TYPE */
          fprintf(cl_aux, "-c_mn %d ", cont->upMTID+1);
          fprintf(cl_aux, "-c_mp %d ", cont->upMPID);
          break;
        case 4:                 /* USER MAT TYPE */
          fprintf(cl_aux, "-c_mn %d ", cont->upMTID+1);
          fprintf(cl_aux, "-c_mp %d ", cont->upMPID);
          fprintf(cl_aux, "-c_md %d ", cont->upMDID);
          break;
        case 5:                 /* USER-DEFINED FUNCTION TYPE */
        case 6:                 /* ANGULAR PARAMETER TYPE */
          /* NOTE:  These are not available via the command line! */
          break;
        default:
          if (loca_in->Cont_Alg != LOCA_LSA_ONLY)
            {
              fprintf(stderr, "%s: Bad cont->upType, %d\n", yo, cont->upType);
              EH(-1,"Bad cont->upType");
            }
          break;                        /* duh */
        }

      fprintf(cl_aux, "\n");

      fclose(cl_aux);
   }
#ifdef PARALLEL
  check_parallel_error("Continuation setup error");
#endif

  /* Call prefront (or mf_setup) if necessary */
  if (Linear_Solver == FRONT)
    {

#ifdef PARALLEL
  if (Num_Proc > 1) EH(-1, "Whoa.  No front allowed with nproc>1");  
  check_parallel_error("Front solver not allowed with nprocs>1");
#endif
          
#ifdef HAVE_FRONT  
      /* Also got to define these because it wants pointers to these numbers */
      max_unk_elem = (MAX_PROB_VAR + MAX_CONC)*MDE;

      one = 1;
      three = 3;

      /* NOTE: We need a overall flag in the vn_glob struct that tells whether FULL_DG
         is on anywhere in domain.  This assumes only one material.  See sl_front_setup for test.
         that test needs to be in the input parser.  */
      if(vn_glob[0]->dg_J_model == FULL_DG) 
        max_unk_elem = (MAX_PROB_VAR + MAX_CONC)*MDE + 4*vn_glob[0]->modes*4*MDE;

       err = mf_setup(&exo->num_elems, 
                     &NumUnknowns[pg->imtrx], 
                     &max_unk_elem, 
                     &three,
                     &one,
                     exo->elem_order_map,
                     fss->el_proc_assign,
                     fss->level,
                     fss->nopdof,
                     fss->ncn,
                     fss->constraint,
                     front_scratch_directory,
                     &fss->ntra); 
      EH(err,"problems in frontal setup ");

#else
      EH(-1,"Don't have frontal solver compiled and linked in");
#endif
    }


  /*
   *  if computing parameter sensitivities, allocate space for solution
   *  sensitivity vectors
   */
        for(i=0;i<nn_post_fluxes_sens;i++)     
          {
            num_pvector=MAX(num_pvector,pp_fluxes_sens[i]->vector_id);
          }
        for(i=0;i<nn_post_data_sens;i++)        
          {
            num_pvector=MAX(num_pvector,pp_data_sens[i]->vector_id);
          }

  if((nn_post_fluxes_sens + nn_post_data_sens) > 0)
    {
      num_pvector++;
      num_pvector = MAX(num_pvector,2);
    /*x_sens_p = Dmatrix_birth(num_pvector,numProcUnknowns);*/
      x_sens_p = (double **) array_alloc(2, num_pvector, numProcUnknowns, sizeof(double));
    }
  else
    x_sens_p = NULL;

  /* Allocate augmenting condition extra unknown arrays */
  if (nAC > 0)
    {
      asdv(&x_AC, nAC);
      asdv(&x_AC_old, nAC);
      asdv(&x_AC_dot, nAC);
    }

  /* Initialize the starting parameter value, lambda */
  if (loca_in->Cont_Alg != LOCA_LSA_ONLY)
    {
      lambda = cont->BegParameterValue;
      delta_s = cont->Delta_s0;
      update_parameterC(0, lambda, x, xdot, x_AC, delta_s, cx, exo, dpi);
    }
  else
    {
      lambda = delta_s = 0.0;
    }
  con_par_ptr = &lambda;

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

  if (strcmp( Matrix_Format, "epetra") == 0) {
    EH(-1, "Error epetra Matrix format not currently supported with loca interface");
  }
  /* Allocate sparse matrix (MSR format) */
  else if( strcmp( Matrix_Format, "msr" ) == 0)
    {
      log_msg("alloc_MSR_sparse_arrays...");
      alloc_MSR_sparse_arrays(&ija, 
                              &a, 
                              &a_old, 
                              0, 
                              node_to_fill, 
                              exo, 
                              dpi);
      /*
       * An attic to store external dofs column names is needed when
       * running in parallel.
       */
      alloc_extern_ija_buffer(num_universe_dofs[pg->imtrx], 
                              num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx], 
                              ija, &ija_attic);
      /*
       * Any necessary one time initialization of the linear
       * solver package (Aztec).
       */
      ams[JAC]->bindx   = ija;
      ams[JAC]->val     = a;
      ams[JAC]->belfry  = ija_attic;
      ams[JAC]->val_old = a_old;
          
      /*
       * These point to nowhere since we're using MSR instead of VBR
       * format.
       */
      ams[JAC]->indx  = NULL;
      ams[JAC]->bpntr = NULL;
      ams[JAC]->rpntr = NULL;
      ams[JAC]->cpntr = NULL;
      ams[JAC]->npn      = dpi->num_internal_nodes + dpi->num_boundary_nodes;
      ams[JAC]->npn_plus = dpi->num_internal_nodes
                         + dpi->num_boundary_nodes + dpi->num_external_nodes;

      ams[JAC]->npu      = num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx];
      ams[JAC]->npu_plus = num_universe_dofs[pg->imtrx];

      ams[JAC]->nnz = ija[num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx] ] - 1;
      ams[JAC]->nnz_plus = ija[num_universe_dofs[pg->imtrx] ];
    }

  /* Allocate sparse matrix (VBR format) */
  else if(  strcmp( Matrix_Format, "vbr" ) == 0)
    {
      log_msg("alloc_VBR_sparse_arrays...");
      alloc_VBR_sparse_arrays (ams[JAC],
                               exo,
                               dpi);
      ija_attic = NULL;
      ams[JAC]->belfry  = ija_attic;

      a = ams[JAC]->val;
      if( !save_old_A ) a_old = ams[JAC]->val_old = NULL;
    }

  /* Allocate sparse matrix (FRONT/ESTIFM format) */
  else if ( strcmp( Matrix_Format, "front") == 0 )
    {
      /* Don't allocate any sparse matrix space when using front */
      ams[JAC]->bindx   = NULL;
      ams[JAC]->val     = NULL;
      ams[JAC]->belfry  = NULL;
      ams[JAC]->val_old = NULL;
      ams[JAC]->indx  = NULL;
      ams[JAC]->bpntr = NULL;
      ams[JAC]->rpntr = NULL;
      ams[JAC]->cpntr = NULL;

    }
  else
    EH(-1,"Attempted to allocate unknown sparse matrix format");

  /* Load initial solution guess */
  init_vec(x, cx, exo, dpi, x_AC, nAC, &timeValueRead);

  /*  if read ACs, update data floats */
  if (nAC > 0 && augc[0].iread == 1)
    {
      for(iAC=0 ; iAC<nAC ; iAC++)
        {
          update_parameterAC(iAC, x, xdot, x_AC, cx, exo, dpi);
        }
    }

  /* Initialize solution sensitivity vectors */
  vzero(numProcUnknowns, &x_sens[0]);
  vzero(numProcUnknowns, &x_sens_temp[0]);

  /* Initialize previous solution vectors */
  dcopy1(numProcUnknowns,x,x_old);
  dcopy1(numProcUnknowns,x_old,x_older);
  dcopy1(numProcUnknowns,x_older,x_oldest);
  if(nAC > 0) dcopy1(nAC,x_AC, x_AC_old);

  /* Initialize linear solver */
  matrix_systems_mask = 1;
  log_msg("sl_init()...");
  sl_init(matrix_systems_mask, ams, exo, dpi, cx);

  /* Make sure the solver was properly initialized on all processors */
#ifdef PARALLEL
  check_parallel_error("Solver initialization problems");
#endif

  ams[JAC]->options[AZ_keep_info] = 1;

  /* set boundary conditions on the initial conditions */
  nullify_dirichlet_bcs();
  find_and_set_Dirichlet(x, xdot, exo, dpi);
  exchange_dof(cx, dpi, x, pg->imtrx);

  /* 
   * Set passdown structure -- variables needed in the argument
   * lists to wrapped routines but not needed in the continuation
   * library.
   */

  NZeros = ams[JAC]->nnz;
  passdown.num_mat_fills  = 0;
  passdown.num_res_fills  = 0;
  passdown.num_linear_its = 0;
  passdown.num_eigen_its = 0;
  passdown.sv_flag = TRUE;
  passdown.nvd = Num_Var_Info_Records;
  passdown.sv_index = NULL;
  passdown.sv_count = NULL;

  /* Initialize mass matrix and shifted matrix as needed */
  if(loca_in->Cont_Alg == HP_CONTINUATION || Linear_Stability)
    {
      if(Linear_Solver == FRONT)
        EH(-1, "Cannot have mass matrix with frontal solver!");
      passdown.mass_matrix = (double *) array_alloc(1, NZeros+5, sizeof(double));
      init_vec_value(passdown.mass_matrix, 0.0, NZeros+5);

  /* Create AZ_MATRIX version for VBR matrix format */
      if (strcmp(Matrix_Format, "vbr") == 0)
        {
          passdown.mmat = AZ_matrix_create(NumUnknowns[pg->imtrx]);
          AZ_set_VBR(passdown.mmat,
                     ams[JAC]->rpntr,
                     ams[JAC]->cpntr,
                     ams[JAC]->bpntr,
                     ams[JAC]->indx,
                     ams[JAC]->bindx,
                     passdown.mass_matrix,
                     ams[JAC]->data_org,
                     0, NULL, AZ_LOCAL);
        }
      else passdown.mmat = NULL;
    }
  else passdown.mass_matrix = NULL;

  /* Arguments to GOMA solver routine */

  passdown.ams		  = ams[JAC];
  passdown.x		  = x;
  passdown.delta_t	  = 0.0;
  passdown.theta	  = 0.0;
  passdown.x_old	  = x_old;
  passdown.x_older	  = x_older;
  passdown.x_oldest	  = x_oldest;
  passdown.xdot		  = xdot;
  passdown.xdot_old	  = xdot_old;
  passdown.resid_vector	  = resid_vector;
  passdown.x_update	  = x_update;
  passdown.scale	  = scale;
  passdown.converged	  = &converged;
  passdown.nprint	  = &nprint;
  passdown.tev		  = tev;
  passdown.tev_post	  = tev_post;
  passdown.rd		  = rd;
  passdown.gindex	  = gindex;
  passdown.gsize	  = p_gsize;
  passdown.gvec	  	  = gvec;
  passdown.gvec_elem	  = gvec_elem;
  passdown.time_value	  = lambda;
  passdown.exo		  = exo;
  passdown.dpi		  = dpi;
  passdown.cx		  = cx;
  passdown.nt		  = 0;
  passdown.path_step_reform = step_reform;
  passdown.is_steady_state  = TRUE;
  passdown.x_AC	  	  = x_AC;
  passdown.x_AC_dot	  = x_AC_dot;
  passdown.lambda	  = lambda;
  passdown.resid_vector_sens = resid_vector_sens;
  passdown.x_sens	  = x_sens;
  passdown.x_sens_temp	  = x_sens_temp;
  passdown.x_sens_p	  = x_sens_p;

  passdown.proc_config    = ams[JAC]->proc_config;
  passdown.options        = ams[JAC]->options;
  passdown.status         = ams[JAC]->status;
  passdown.LSA_flag       = FALSE;
  passdown.last_step      = FALSE;


  /* This is required for VBR matrix-vector multiply calls: */
  if (strcmp(Matrix_Format, "vbr") == 0)
    {
      passdown.amat = AZ_matrix_create(NumUnknowns[pg->imtrx]);
      AZ_set_VBR(passdown.amat,
                 passdown.ams->rpntr,
                 passdown.ams->cpntr,
                 passdown.ams->bpntr,
                 passdown.ams->indx,
                 passdown.ams->bindx,
                 passdown.ams->val,
                 passdown.ams->data_org,
                 0, NULL, AZ_LOCAL);
    }
  else passdown.amat = NULL;

  /* Open ASCII output file Soln_OutFile */
  if(strlen(Soln_OutFile))
    {
      passdown.file = fopen(Soln_OutFile, "w");
      if (passdown.file == NULL) {
        DPRINTF(stdout, "%s:  opening soln file for writing\n", yo);
        EH(-1, "\t");
      }
    }
#ifdef PARALLEL
  check_parallel_error("Soln output file error");
#endif

  /* Determine number of LSA passes to perform */
  /* Open eigenvector ExodusII output files for writing */
  passdown.do_3D_of_2D = FALSE;
  if (Linear_Stability)
    {

  /* Check for ARPACK and/or PARPACK */
#ifndef HAVE_PARPACK
      if (Num_Proc > 1)
        {
          WH(-1, "PARPACK was not compiled in -- cannot do eigensolves!");
          Linear_Stability = LSA_NONE;
        }
#endif
#ifndef HAVE_ARPACK
      WH(-1, "ARPACK was not compiled in -- cannot do eigensolves!");
      Linear_Stability = LSA_NONE;
#endif
    
      if (Linear_Stability == LSA_3D_OF_2D ||
          Linear_Stability == LSA_3D_OF_2D_SAVE)
        {
   /* Temporary error message - EDW */
          passdown.do_3D_of_2D = TRUE;
          n = LSA_number_wave_numbers;
        }
      else n = 1;
      err = create_eigen_outfiles(passdown.exo,
                                  passdown.dpi,
                                  passdown.rd);
      EH(err, "Unable to open eigenvector output files!");
    }

  /********************************************************* 
   *     End of passdown structure Initialization          * 
   *********************************************************/
 
  /********************************************************* 
   *    Begin loading of 'con' structure of input info     * 
   *********************************************************/

  /*
   * Set variables in the con.general_info structure:
   * Initial guess for solution vector,
   * number of unknowns on this processor, 
   * number of unknowns owned by this processor,
   * whether or not this processor prints, and 
   * what continuation method is requested.
   * (con.general_info.method set below)
   */
  
  con.general_info.param        = lambda;
  con.general_info.x            = passdown.x;
  con.general_info.numUnks      = NumUnknowns[pg->imtrx] + NumExtUnknowns[pg->imtrx];
  con.general_info.numOwnedUnks = NumUnknowns[pg->imtrx];
  con.general_info.perturb      = loca_in->perturb;
  if (ProcID == 0)
    {
      if (loca_in->Cont_Alg == LOCA_LSA_ONLY)
        {
          con.general_info.printproc = 5;
        }
      else
        {
          con.general_info.printproc = loca_in->debug;
        }
    }
  else
    {
      con.general_info.printproc = 0;
    }

  /*
   * Set variables in the con.stepping_info structure:
   * These variable control the parameter stepping algorithm
   * (see following arclength_info structure as well)
   */

  con.stepping_info.first_step     = cont->Delta_s0;
  con.stepping_info.base_step      = 1;
  con.stepping_info.max_steps      = cont->MaxPathSteps - 1;
  con.stepping_info.last_step      = FALSE;
  con.stepping_info.max_param      = cont->EndParameterValue;
  con.stepping_info.min_delta_p    = cont->Delta_s_min;
  con.stepping_info.max_delta_p    = cont->Delta_s_max;
  con.stepping_info.step_ctrl      = loca_in->StepAggr;
  con.stepping_info.max_newton_its = Max_Newton_Steps;

  /*
   * The next several con sub-structures are method (strategy) dependent
   * Only one of these sub-structures should be filled for each run
   */

  /*
   * First, initialize arrays in these substructures to NULL.
   * Only those which are needed will be allocated.
   */
  con.turning_point_info.nv = NULL;
  con.pitchfork_info.psi = NULL;
  con.hopf_info.y_vec = NULL;
  con.hopf_info.z_vec = NULL;

  con.general_info.nv_restart = FALSE;
  con.general_info.nv_save = FALSE;
  switch (loca_in->Cont_Alg) {
    case CONTINUATION:
      if (loca_in->Cont_Order == 0)
        con.general_info.method = ZERO_ORDER_CONTINUATION;

      else if (loca_in->Cont_Order == 1)
        con.general_info.method = FIRST_ORDER_CONTINUATION;

      else if (loca_in->Cont_Order == 2)
        con.general_info.method = ARC_LENGTH_CONTINUATION;

  /* Set variables in the con.arclength_info structure: */
	con.arclength_info.dp_ds2_goal	    = loca_in->DpDs2;
	con.arclength_info.dp_ds_max	    = loca_in->DpDsHi;
	con.arclength_info.tang_exp	    = loca_in->Texp;
	con.arclength_info.tang_step_limit = loca_in->MaxTS;
      break;

    case TP_CONTINUATION: 
      con.general_info.method = TURNING_POINT_CONTINUATION;
      if (loca_in->NVSave) con.general_info.nv_save = TRUE;

  /* Set variables in the con.turning_point_info structure: */
      con.turning_point_info.bif_param = loca_in->TPGuess;
      if (loca_in->NVRestart)
        {
          con.turning_point_info.nv = (double *) array_alloc (1,
                                      con.general_info.numUnks, sizeof(double));
          con.general_info.nv_restart = TRUE;
          if (Num_Proc > 1)
            multiname(loca_in->NV_exoII_infile, ProcID, Num_Proc);
          DPRINTF(stdout, "Reading previous null vector ...\n");
          err = rd_vectors_from_exoII(con.turning_point_info.nv,
				      loca_in->NV_exoII_infile, 0, 0, INT_MAX, &timeValueRead);
          if (err != 0)
            {
              DPRINTF(stderr, "do_loca:  err from rd_vectors_from_exoII\n");
              con.general_info.nv_restart = FALSE;
            }
        }
      break;

    case PF_CONTINUATION: 
      con.general_info.method = PITCHFORK_CONTINUATION;

  /* Set variables in the con.pitchfork_info structure: */
      con.pitchfork_info.bif_param = loca_in->TPGuess;
      con.pitchfork_info.psi = (double *)
	array_alloc (1, con.general_info.numUnks, sizeof(double));
      if (Num_Proc > 1)
        multiname(loca_in->NV_exoII_infile, ProcID, Num_Proc);
      DPRINTF(stdout, "Reading previous null vector ...\n");
      err = rd_vectors_from_exoII(con.pitchfork_info.psi,
                                  loca_in->NV_exoII_infile, 0, 0, INT_MAX, &timeValueRead);
      if (err != 0)
        {
          DPRINTF(stderr, "do_loca:  err from rd_vectors_from_exoII\n");
          exit(-1);
        }
      break;

    case HP_CONTINUATION: 
      con.general_info.method = HOPF_CONTINUATION;

  /* Make sure Komplex library is compiled for Hopf tracking */
#ifndef KOMPLEX
      EH(-1, "Hopf Tracking Algorithm Requires Komplex Library!\n"
           "Recompile with KOMPLEX flag set.\n");
#endif

  /* Set variables in the con.hopf_info structure: */
      con.hopf_info.bif_param = loca_in->TPGuess;
      con.hopf_info.omega     = loca_in->omega;
   /* con.hopf_info.mass_flag = loca_in->mass_flag; */
      con.hopf_info.mass_flag = TRUE;

  /* Get Exodus file names for both parts of initial null vector */
      if (Num_Proc > 1)
        {
          multiname(loca_in->NV_exoII_infile, ProcID, Num_Proc);
          multiname(loca_in->NV_imag_infile,  ProcID, Num_Proc);
        }

  /* Allocate arrays for real and imaginary parts of null vector */
      con.hopf_info.y_vec = (double *) array_alloc (1,
                             con.general_info.numUnks, sizeof(double));
      con.hopf_info.z_vec = (double *) array_alloc (1,
                             con.general_info.numUnks, sizeof(double));

  /* Load y_vec and z_vec from these files */
      DPRINTF(stdout, "Reading previous null vector (real part) ...\n");
      err = rd_vectors_from_exoII(con.hopf_info.y_vec, 
				  loca_in->NV_exoII_infile, 0, 0, INT_MAX, &timeValueRead);
      if (err != 0) EH(-1, "do_loca: error reading real part of null vector");
      DPRINTF(stdout, "Reading previous null vector (imaginary part) ...\n");
      err = rd_vectors_from_exoII(con.hopf_info.z_vec, 
				  loca_in->NV_imag_infile, 0, 0, INT_MAX, &timeValueRead);
      if (err != 0) EH(-1, "do_loca: error reading imag. part of null vector");

  /* If using MSR matrix format, instantiate amat (struct AZ_MATRIX).
     VBR case was already handled above. */
      if (strcmp(Matrix_Format, "msr") == 0)
        {
          passdown.amat = AZ_matrix_create(NumUnknowns[pg->imtrx]);
          AZ_set_MSR(passdown.amat, passdown.ams->bindx, passdown.ams->val,
                     passdown.ams->data_org, 0, NULL, AZ_LOCAL);
        }
      break;

/*
    case SQP_OPTIMIZATION:
#ifdef SQP_OPTIMIZER
      EH(-1, "sqp optimization not yet available in Goma!");
      solve_sqp_optimization(con_par_ptr, passdown.mesh, nstep);
#else
      EH(-1, "sqp optimization requested but not compiled in!");
#endif
      break;
*/

  /* This is for a single steady state followed by LSA */
    case LOCA_LSA_ONLY:
      con.general_info.method = LOCA_LSA_ONLY;
      con.general_info.param = 0.0;
      con.stepping_info.max_steps = 0;
      con.stepping_info.last_step = TRUE;
      break;

    default:
      printf("ERROR %s: Unknown Continuation method: %d\n",yo,
			loca_in->Cont_Alg); exit(-1);
  }
  passdown.method = con.general_info.method;

  /* Set variables in the con.eigen_info structure for eigenvalue calcs */

    con.eigen_info.Num_Eigenvalues = 0;
    if (Linear_Stability)
      {
        con.eigen_info.Num_Eigenvalues = eigen->Eigen_NEV_WANT;
        con.eigen_info.Num_Eigenvectors = eigen->Eigen_Record_Modes;
        if (con.eigen_info.Num_Eigenvalues > 0)
          {
            con.eigen_info.Shift_Point[0]  = eigen->Eigen_Cayley_Sigma;
            if (eigen->Eigen_Algorithm == LSA_CAYLEY)
              {
                con.eigen_info.Shift_Point[1] = eigen->Eigen_Cayley_Mu;
              }
            else
              {
                con.eigen_info.Shift_Point[1] = eigen->Eigen_Cayley_Sigma;
              }
            con.eigen_info.Shift_Point[2]  = 1.0;   /* This is a default */
            con.eigen_info.Shift_Point[3]  = eigen->Eigen_SI_Tol_Param;
            con.eigen_info.Arnoldi         = eigen->Eigen_Krylov_Subspace;
            con.eigen_info.Residual_Tol[0] = eigen->Eigen_Relative_Tol;
            con.eigen_info.Residual_Tol[1] = eigen->Eigen_Linear_Tol;
            con.eigen_info.Max_Iter   = eigen->Eigen_Maximum_Outer_Iterations;
            con.eigen_info.Every_n_Steps   = eigen->Eigen_Solve_Freq;
            con.eigen_info.sort            = TRUE;
          }
        else EH(-1, "Number of eigenvalues must be specified!");
      }

  /* Check starting mesh element quality if requested */
    if (nEQM > 0)
      {
        DPRINTF(stdout, "\nINITIAL ELEMENT QUALITY CHECK---\n");
        element_quality(exo, x, ams[0]->proc_config);
      }

  /* First, handle single pass for eigensolver (if LOCA_LSA_ONLY) */

  if (con.general_info.method == LOCA_LSA_ONLY)
    {
      initialize_util_routines(NumUnknowns[pg->imtrx], NumUnknowns[pg->imtrx] + NumExtUnknowns[pg->imtrx]);
      con.private_info.step_num = 0;
      err = nonlinear_solver_conwrap(x, (void *)&con, 0, 0.0, 0.0);
      solution_output_conwrap(1, x, 0.0, NULL, 0.0, NULL, 0.0, 0, err, &con);
  printf("LOCA LSA ONLY HAS FINISHED: \n");
    }

  /* Otherwise, now call continuation library and return */

  else nstep = con_lib(&con);

  *con_par_ptr = con.general_info.param;
  
  if (con.general_info.printproc)
    print_final(con.general_info.param, nstep, passdown.num_mat_fills,
		passdown.num_res_fills, passdown.num_linear_its);

  if (strlen(Soln_OutFile)) fclose(passdown.file);

  /* Write null vector to file if requested */
  if (loca_in->NVSave && (passdown.method == TURNING_POINT_CONTINUATION
                          || passdown.method == PITCHFORK_CONTINUATION) )
    {

  /* Open a new ExodusII file */
      if (Num_Proc > 1)
        multiname(loca_in->NV_exoII_outfile, ProcID, Num_Proc);
      one_base(exo);
      wr_mesh_exo(exo, loca_in->NV_exoII_outfile, 0);
      wr_result_prelim_exo(rd, exo, loca_in->NV_exoII_outfile, NULL);
      if (Num_Proc > 1) wr_dpi(dpi, loca_in->NV_exoII_outfile, 0);
      *passdown.nprint = 0;

  /* Write the null vector to this file */
        write_solution(loca_in->NV_exoII_outfile,
                       passdown.resid_vector,
                       con.private_info.x_tang,
                       passdown.x_sens_p,
                       passdown.x_old,
                       passdown.xdot,
                       passdown.xdot_old,
                       passdown.tev,
                       passdown.tev_post,
                       NULL,
                       passdown.rd,
                       passdown.gvec,
                       passdown.gvec_elem,
                       passdown.nprint,
                       0.0,
                       passdown.theta,
                       0.0,
                       NULL,
                       passdown.exo,
                       passdown.dpi);
      zero_base(exo);
    }

  /* Write null vector to files if requested (Hopf tracking case) */
#ifdef KOMPLEX
  if (loca_in->NVSave && passdown.method == HOPF_CONTINUATION)
    {

  /* Open two new ExodusII files */
      if (Num_Proc > 1)
        {
          multiname(loca_in->NV_exoII_outfile, ProcID, Num_Proc);
          multiname(loca_in->NV_imag_outfile, ProcID, Num_Proc);
        }
      one_base(exo);

  /* First write real part (y_vec) to NV_exoII_outfile */
      wr_mesh_exo(exo, loca_in->NV_exoII_outfile, 0);
      wr_result_prelim_exo(rd, exo, loca_in->NV_exoII_outfile, NULL);
      if (Num_Proc > 1) wr_dpi(dpi, loca_in->NV_exoII_outfile, 0);
      *passdown.nprint = 0;

      write_solution(loca_in->NV_exoII_outfile,
                     passdown.resid_vector,
                     con.hopf_info.y_vec,
                     passdown.x_sens_p, 
                     passdown.x_old,
                     passdown.xdot,
                     passdown.xdot_old,
                     passdown.tev,
                     passdown.tev_post,
		     NULL,
                     passdown.rd,
                     passdown.gindex,
                     passdown.gsize,
                     passdown.gvec,
                     passdown.gvec_elem,
                     passdown.nprint, 
                     0.0,
                     passdown.theta,
                     0.0,
                     NULL,
                     passdown.exo,
                     passdown.dpi);

  /* Then write imaginary part (z_vec) to NV_imag_outfile */
      wr_mesh_exo(exo, loca_in->NV_imag_outfile, 0);
      wr_result_prelim_exo(rd, exo, loca_in->NV_imag_outfile, NULL);
      if (Num_Proc > 1) wr_dpi(dpi, loca_in->NV_imag_outfile, 0);
      *passdown.nprint = 0;

      write_solution(loca_in->NV_imag_outfile,
                     passdown.resid_vector,
                     con.hopf_info.z_vec,
                     passdown.x_sens_p, 
                     passdown.x_old,
                     passdown.xdot,
                     passdown.xdot_old,
                     passdown.tev,
                     passdown.tev_post,
		     NULL,
                     passdown.rd,
                     passdown.gindex,
                     passdown.gsize,
                     passdown.gvec,
                     passdown.gvec_elem,
                     passdown.nprint, 
                     0.0,
                     passdown.theta,
                     0.0,
                     NULL,
                     passdown.exo,
                     passdown.dpi);

      zero_base(exo);
    }
#endif

  /* Deallocate local arrays */

  safer_free( (void **) &ROT_Types);
  safer_free( (void **) &node_to_fill);
  safer_free( (void **) &resid_vector);
  safer_free( (void **) &resid_vector_sens);
  safer_free( (void **) &scale);
  safer_free( (void **) &x);
  safer_free( (void **) &x_old);
  safer_free( (void **) &x_older);
  safer_free( (void **) &x_oldest);
  safer_free( (void **) &xdot);
  safer_free( (void **) &xdot_old);
  safer_free( (void **) &x_update);
  safer_free( (void **) &x_sens);
  safer_free( (void **) &x_sens_temp);
  if (nAC > 0)
    {
      safer_free( (void **) &x_AC);
      safer_free( (void **) &x_AC_old);
      safer_free( (void **) &x_AC_dot);
    }

  if((nn_post_data_sens+nn_post_fluxes_sens) > 0)
		safer_free( (void**) &x_sens_p);
/*		Dmatrix_death(x_sens_p,num_pvector,numProcUnknowns);*/

  sl_free(matrix_systems_mask, ams);
  for (i = 0; i < NUM_ALSS; i++)
    {
      safer_free((void **) &(ams[i]));
    }
  safer_free( (void **) &rd);
  safer_free( (void **) &gvec);
  safer_free( (void **) &cpcc);
  if (tpcc != NULL) safer_free( (void **) &tpcc);

  i = 0;
  for ( eb_indx = 0; eb_indx < exo->num_elem_blocks; eb_indx++ )
    {
      if (exo->elem_var_tab != NULL)
        {
          for ( ev_indx = 0; ev_indx < passdown.rd->nev; ev_indx++ )
            {
	      if (exo->elem_var_tab[i++] == 1)
                {
	          safer_free((void **) &(gvec_elem [eb_indx][ev_indx]));
                }
            }
        }
      safer_free((void **) &(gvec_elem [eb_indx]));
    }

  for(i = 0; i < MAX_NUMBER_MATLS; i++) {
    for(n = 0; n < MAX_MODES; n++) {
      safer_free((void **) &(ve_glob[i][n]->gn));
      safer_free((void **) &(ve_glob[i][n]));
    }
    safer_free((void **) &(vn_glob[i]));
  }

  safer_free( (void **) &gvec_elem);
  safer_free( (void **) &Local_Offset);
  safer_free( (void **) &Dolphin);

  safer_free((void **) &passdown.sv_index);
  safer_free((void **) &passdown.sv_count);
  safer_free((void **) &passdown.mass_matrix);
  if (passdown.amat != NULL) AZ_matrix_destroy(&passdown.amat);
  if (passdown.mmat != NULL) AZ_matrix_destroy(&passdown.mmat);

  safer_free((void **) &con.turning_point_info.nv);
  safer_free((void **) &con.pitchfork_info.psi);
  safer_free((void **) &con.hopf_info.y_vec);
  safer_free((void **) &con.hopf_info.z_vec);

  free(pg->matrices);

  return nstep;
} /**************** END of do_loca() *****************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int nonlinear_solver_conwrap(double *x, void *con_ptr, int step_num,
			     double lambda, double delta_s)
/* Put the call to your nonlinear solver here.
 * Input:
 *    x         solution vector
 *    con_ptr   pointer to continuation structure, cast to (void *)
 *              must be passed to nonlinear solver and then passed
 *              to bordering algorithms.
 *    step_num  Continuation step number
 *
 * Output:
 *
 * Return Value:
 *    num_newt_its  Number of Newton iterations needed for
 *                  convergence, used to pick next step size.
 *                  Negative value means nonlinear solver didn't converge.
 *
 */
{
  struct Aztec_Linear_Solver_System *ams = &(passdown.ams[JAC]);
  int  err;
  int converged;
  int nits=0; /* num_modnewt_its=0;  */
  int i, iAC;
  int iCC = 0, iTC = 0, iUC = 0, nCC = 0;
  double theta=0.0;
  double evol_local=0.0;
#ifdef PARALLEL
  double evol_global=0.0;
#endif
  char *yo = "do_loca";

/* set up boundary conditions */
  nullify_dirichlet_bcs();
  find_and_set_Dirichlet (x, passdown.xdot, passdown.exo, passdown.dpi);
  exchange_dof(passdown.cx, passdown.dpi, x, pg->imtrx);

/* show continuation type */
  if (ProcID == 0) {
    fprintf(stdout, "\n\t----------------------------------");
    switch (passdown.method) {
      case ZERO_ORDER_CONTINUATION:
	DPRINTF (stdout, "\n\tZero order continuation:");
	break;
      case FIRST_ORDER_CONTINUATION:
	DPRINTF (stdout, "\n\tFirst order continuation:");
	break;
      case ARC_LENGTH_CONTINUATION:
	DPRINTF (stdout, "\n\tArc length continuation:");
	break;
      case TURNING_POINT_CONTINUATION:
	DPRINTF (stdout, "\n\tTurning point continuation:");
	break;
      case PITCHFORK_CONTINUATION:
	DPRINTF (stdout, "\n\tPitchfork continuation:");
	break;
      case HOPF_CONTINUATION:
	DPRINTF (stdout, "\n\tHopf continuation:");
	break;
      case LOCA_LSA_ONLY:
        DPRINTF (stdout, "\n\tLinear stability analysis:\n");
	break;
      default:
	DPRINTF (stderr, "%s: Bad Continuation, %d\n", yo, passdown.method);
	exit(-1);
	break;
    }

/* Print step information */
    if (passdown.method != LOCA_LSA_ONLY)
      {
        nCC = cpcc[0].nCC;
        DPRINTF (stdout, "\n\tStep number: %4d of %4d (max)",
			      step_num+1, cont->MaxPathSteps);
            theta = (lambda - cont->BegParameterValue)
                  / (cont->EndParameterValue - cont->BegParameterValue);
        if (nCC > 1 || nUC > 0)
          {
            DPRINTF (stdout, "\n\tAttempting solution at: theta = %g", theta);
          }
        else
          {
            DPRINTF (stdout, "\n\tAttempting solution at:");
          }
        if (cont->upType == 5 && nUC > 0)
          {
            for (iUC=0; iUC<nUC; iUC++)
              {
                switch(cpuc[iUC].Type) {
                  case 1:		/* BC */
    	            DPRINTF (stdout, "\n\tBCID=%3d DFID=%5d",
                             cpuc[iUC].BCID, cpuc[iUC].DFID);
    	            break;
                  case 2:		/* MT */
    	            DPRINTF (stdout, "\n\tMTID=%3d MPID=%5d",
                             cpuc[iUC].MTID, cpuc[iUC].MPID);
  	            break;
                  case 3:		/* AC */
     	            DPRINTF (stdout, "\n\tACID=%3d DFID=%5d",
                             cpuc[iUC].BCID, cpuc[iUC].DFID);
     	            break;
                  case 4:		/* UM */
      	            DPRINTF (stdout, "\n\tMTID=%3d MPID=%5d FLOAT=%3d",
                             cpuc[iUC].MTID, cpuc[iUC].MPID, cpuc[iUC].MDID);
	            break;
                  case 5:               /* UF */
	            break;
                  default:
	            DPRINTF (stderr, "%s: Bad condition type, %d\n",
                             yo, cpcc[iUC].Type);
	            exit(-1);
	            break;
                }
                DPRINTF(stdout, " Parameter= %10.6e delta_s= %10.6e",
                        cpuc[iUC].value, cpuc[iUC].value-cpuc[iUC].old_value);
              }
          }
        else
          {
            for (iCC=0; iCC<nCC; iCC++)
              {
                switch(cpcc[iCC].Type) {
                  case 1:		/* BC */
    	            DPRINTF (stdout, "\n\tBCID=%3d DFID=%5d",
                             cpcc[iCC].BCID, cpcc[iCC].DFID);
    	            break;
                  case 2:		/* MT */
    	            DPRINTF (stdout, "\n\tMTID=%3d MPID=%5d",
                             cpcc[iCC].MTID, cpcc[iCC].MPID);
  	            break;
                  case 3:		/* AC */
     	            DPRINTF (stdout, "\n\tACID=%3d DFID=%5d",
                             cpcc[iCC].BCID, cpcc[iCC].DFID);
     	            break;
                  case 4:		/* UM */
      	            DPRINTF (stdout, "\n\tMTID=%3d MPID=%5d FLOAT=%3d",
                             cpcc[iCC].MTID, cpcc[iCC].MPID, cpcc[iCC].MDID);
	            break;
                  case 5:               /* UF */
	            break;
                  case 6:               /* AN */
                    if (iCC == 0) DPRINTF (stdout, "\n\tAngular");
	            break;
                  default:
	            DPRINTF (stderr, "%s: Bad condition type, %d\n",
                             yo, cpcc[iCC].Type);
	            exit(-1);
	            break;
                 }
               DPRINTF(stdout, " Parameter= %10.6e delta_s= %10.6e",
 /*                    cpcc[iCC].value, cpcc[iCC].ratio * delta_s); */
                       cpcc[iCC].value, cpcc[iCC].value - cpcc[iCC].old_value);
              }
          }
      }                 /* End of "passdown.method != LOCA_LSA_ONLY" block */
  }			/* End of print block */
        
    passdown.theta = tran->theta;
        
#ifdef DEBUG
    DPRINTF(stderr, "%s: starting solve_nonlinear_problem\n", yo);
#endif

    err = solve_nonlinear_problem(ams,
				  x, 
				  passdown.delta_t, 
				  passdown.theta,
				  passdown.x_old,
				  passdown.x_older, 
				  passdown.xdot,
				  passdown.xdot_old,
				  passdown.resid_vector, 
				  passdown.x_update,
				  passdown.scale, 
				 &converged,
				  passdown.nprint,
				  passdown.tev, 
				  passdown.tev_post,
				  NULL,
				  passdown.rd,
				  passdown.gindex,
				  passdown.gsize,
				  passdown.gvec, 
				  passdown.gvec_elem, 
				  lambda,
				  passdown.exo, 
				  passdown.dpi, 
				  passdown.cx, 
				  0, 
				  passdown.path_step_reform,
				  passdown.is_steady_state,
				  passdown.x_AC, 
 				  passdown.x_AC_dot, 
				  lambda, 
				  passdown.resid_vector_sens, 
				  passdown.x_sens_temp,
				  passdown.x_sens_p,
                                  con_ptr);

#ifdef DEBUG
    fprintf(stderr, "%s: returned from solve_nonlinear_problem\n", yo);
#endif

/* Bail out if a deformation gradient occurs */
    if (err == -1)
      {
        converged = FALSE;
        return err;
      }
    nits = err;

/* Write converged solution */
    if (converged)
      {
        if ( Write_Intermediate_Solutions == 0 && Unlimited_Output ) {
#ifdef DEBUG
         DPRINTF(stderr, "%s: write_solution call WIS\n", yo);
#endif
            write_solution(ExoFileOut,
                           passdown.resid_vector,
                           x,
                           passdown.x_sens_p,
                           passdown.x_old,
                           passdown.xdot,
                           passdown.xdot_old,
                           passdown.tev,
                           passdown.tev_post,
                           NULL,
                           passdown.rd,
                           passdown.gvec,
                           passdown.gvec_elem,
                           passdown.nprint,
                           delta_s,
                           passdown.theta,
                           passdown.lambda,
                           NULL,
                           passdown.exo,
                           passdown.dpi);
#ifdef DEBUG
         DPRINTF(stderr, "%s: write_solution end call WIS\n", yo);
#endif
       }

    DPRINTF(stdout,
            "\n\tStep accepted, theta (proportion complete) = %10.6e\n",
            theta);

     /* Save continuation parameter values */
     if (passdown.method != LOCA_LSA_ONLY)
       {
         for (iCC = 0; iCC < nCC; iCC++) cpcc[iCC].old_value = cpcc[iCC].value;
         for (iTC = 0; iTC < nTC; iTC++) tpcc[iTC].old_value = tpcc[iTC].value;
       }

        /* PRINT OUT VALUES OF EXTRA UNKNOWNS FROM AUGMENTING CONDITIONS */

        if (nAC > 0)
          {
            DPRINTF(stdout, "\n------------------------------\n");
            DPRINTF(stdout, "Augmenting Conditions:    %4d\n", nAC);
            DPRINTF(stdout, "Number of extra unknowns: %4d\n\n", nAC);

            for (iAC = 0; iAC < nAC; iAC++)
             {
              if (augc[iAC].Type == AC_USERBC)
               {
                DPRINTF(stdout, "\tBC[%4d] DF[%4d] = %10.6e\n",
                        augc[iAC].BCID, augc[iAC].DFID, passdown.x_AC[iAC]);
               }
              else if (augc[iAC].Type == AC_USERMAT ||
                       augc[iAC].Type == AC_FLUX_MAT )
               {
                DPRINTF(stdout, "\n New MT[%4d] MP[%4d] = %10.6e\n",
                        augc[iAC].MTID, augc[iAC].MPID, passdown.x_AC[iAC]);
               }
              else if(augc[iAC].Type == AC_VOLUME)
               {
                evol_local = augc[iAC].evol;
#ifdef PARALLEL
                if( Num_Proc > 1 ) {
                     MPI_Allreduce( &evol_local, &evol_global, 1,
                                    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                }
                evol_local = evol_global;
#endif
                DPRINTF(stdout, "\tMT[%4d] VC[%4d]=%10.6e Param=%10.6e\n",
                        augc[iAC].MTID, augc[iAC].VOLID, evol_local,
                        passdown.x_AC[iAC]);
               }
	      else if(augc[iAC].Type == AC_POSITION)
               {
                evol_local = augc[iAC].evol;
#ifdef PARALLEL
                if( Num_Proc > 1 ) {
                     MPI_Allreduce( &evol_local, &evol_global, 1,
                                    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                }
                evol_local = evol_global;
#endif
                DPRINTF(stdout, "\tMT[%4d] XY[%4d]=%10.6e Param=%10.6e\n",
                        augc[iAC].MTID, augc[iAC].VOLID, evol_local,
                        passdown.x_AC[iAC]);
               }
              else if(augc[iAC].Type == AC_FLUX)
               {
                DPRINTF(stdout, "\tBC[%4d] DF[%4d]=%10.6e\n",
                        augc[iAC].BCID, augc[iAC].DFID, passdown.x_AC[iAC]);
               }
             }

          }

      /* Check element quality */
      element_quality(passdown.exo, x, ams->proc_config);

          /* INTEGRATE FLUXES, FORCES */

          for (i = 0; i < nn_post_fluxes; i++)
            evaluate_flux (passdown.exo,
			   passdown.dpi,
			   pp_fluxes[i]->ss_id,
			   pp_fluxes[i]->flux_type ,
			   pp_fluxes[i]->flux_type_name ,
			   pp_fluxes[i]->blk_id ,
			   pp_fluxes[i]->species_number,
			   pp_fluxes[i]->flux_filenm,
			   pp_fluxes[i]->profile_flag,
			   x,
			   passdown.xdot,
			   NULL,
			   delta_s,
			   lambda,
			   1);

          /* COMPUTE FLUX, FORCE SENSITIVITIES */

          for (i = 0; i < nn_post_fluxes_sens; i++)
            evaluate_flux_sens (passdown.exo,
				passdown.dpi,
				pp_fluxes_sens[i]->ss_id,
				pp_fluxes_sens[i]->flux_type ,
				pp_fluxes_sens[i]->flux_type_name ,
				pp_fluxes_sens[i]->blk_id ,
				pp_fluxes_sens[i]->species_number,
				pp_fluxes_sens[i]->sens_type,
				pp_fluxes_sens[i]->sens_id,
				pp_fluxes_sens[i]->sens_flt,
				pp_fluxes_sens[i]->sens_flt2,
				pp_fluxes_sens[i]->vector_id,
				pp_fluxes_sens[i]->flux_filenm,
				pp_fluxes_sens[i]->profile_flag,
				x,
				passdown.xdot,
				passdown.x_sens_p,
				delta_s,
				lambda,
				1);

      }   /*  end of if converged block  */

          
  passdown.num_res_fills += nits+1;
  passdown.num_mat_fills += nits;
  if (Linear_Solver == AZTEC) passdown.num_linear_its += ams->status[AZ_its];

  if (!converged) return (-nits);
  else return (nits);

}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int linear_solver_conwrap(double *x, int jac_flag, double *tmp)
/* Put the call to your linear solver here. There are three options
 * about the reuse of the preconditioner. It is always safe to 
 * just solve the matrix from scratch.
 * Input:
 *    x          Right hand side
 *    jac_flag   Flag indicating the status of the Jacobian so that
 *               preconditioners can be used: 
 *               NEW_JACOBIAN:   recalculate preconditioner
 *               OLD_JACOBIAN:   reuse preconditioner
 *               SAME_BUT_UNSCALED_JACOBIAN: Must rescale the matrix and
 *                               can reuse preconditioner. This happens
 *                               when the matrix has been recalculated
 *                               at the same conditions as before.
 *               CHECK_JACOBIAN: Same matrix, but rebuild preconditioner anyway.
 *		  NOTE: For now, the preconditioner is always recalculated.
 *    tmp        Work space array same length as x, only used for
 *               the SAME_BUT_UNSCALED_JACOBIAN option.
 *    rescale    Flag indicating if scale vector needs to be
 *		 recalculated (AZTEC only).
 *
 * Output:
 *    x          Solution vector
 *
 * Return Value:
 *    Negative value means linear solver didn't converge.
 */
{
  struct Aztec_Linear_Solver_System *ams = &(passdown.ams[JAC]);
  static int first_linear_solver_call=FALSE;
  double *a   = ams->val;       /* nonzero values of a CMSR matrix */
  int    *ija = ams->bindx;     /* column pointer array into matrix "a" */
  static int    Factor_Flag;    /* For UMFPACK */
  int           matr_form;      /* 1: MSR FORMAT MATRIX FOR UMFPACK DRIVER */
  int           error = 0;
  int           why = 0;
  char          stringer[80];   /* holding format of num linear solve itns */
  dbl           s_start;        /* mark start of solve */
  dbl           s_end;          /* mark end of solve */

  int   linear_solver_blk;      /* count calls to AZ_solve() */
  int   linear_solver_itns;     /* count cumulative linearsolver iterations */
  int   num_linear_solve_blks;  /* one pass for now */
  int   matrix_solved;          /* boolean */

/* Additional values for frontal solver */
#ifdef HAVE_FRONT
  int mf_resolve;
  dbl smallpiv;
  dbl singpiv;
  int iautopiv;
  int iscale;   /* you will have to turn this off for resolves */
  dbl scaling_max;
  dbl h_elem_avg;                        /* global average element size for PSPG */
  dbl U_norm    ;                        /* global average velocity for PSPG */
#endif

  int numUnks      = NumUnknowns[pg->imtrx] + NumExtUnknowns[pg->imtrx];
  double *xr=NULL;


/* Create temporary vector to pass RHS to solver - EDW 6/1/2000 */
  xr = (double *) array_alloc(1, numUnks, sizeof(double));
  dcopy1(numUnks, x, xr);

/* Rescale RHS for scaled Jacobian */

  if (Linear_Solver != FRONT)
    {
      row_sum_scaling_scale(ams, xr, passdown.scale);
    }

/* Call chosen linear solver */

      s_start = ut();
      switch (Linear_Solver)
        {
        case UMFPACK2:
        case UMFPACK2F:
          if (strcmp(Matrix_Format, "msr"))
            EH(-1,"ERROR: umfpack solver needs msr matrix format");

          Factor_Flag = 1;
          if (Linear_Solver == UMFPACK2F) Factor_Flag = 0;
          /*  */
          matr_form = 1;
          LOCA_UMF_ID = SL_UMF(LOCA_UMF_ID,
                 &first_linear_solver_call, &Factor_Flag, &matr_form, 
                 &NumUnknowns[pg->imtrx], &NZeros, &ija[0], &ija[0], &a[0],
                 &xr[0], &x[0]);
          /*  */
          first_linear_solver_call = FALSE;
          strcpy(stringer, " 1 ");
          break;

        case SPARSE13a:
          if (strcmp(Matrix_Format, "msr"))
            EH(-1,"ERROR: lu solver needs msr matrix format");

	  dcopy1(NumUnknowns[pg->imtrx], xr, x);
          lu(NumUnknowns[pg->imtrx], NumExtUnknowns[pg->imtrx], NZeros, a, ija, x, 2);
          first_linear_solver_call = FALSE;
      /* 
       * Note that sl_lu has static variables to keep track of
       * first call or not.
       */

          strcpy(stringer, " 1 ");
          break;
          
        case AZTEC:

	  /* Set option of preconditioner reuse */
	    
/*    if ( first_linear_solver_call )
        {
          ams->options[AZ_pre_calc] = AZ_calc;
        }
      else
        { */

          if ( strcmp(Matrix_Factorization_Reuse, "calc") == 0 )
            {
              /*
               * Gonna start from scratch even though I've cooked a
               * preconditioner in the kitchen all day? Well, then
               * you won't need the leftover pieces from all my
               * hard preparation last time around.
               */
              
              AZ_free_memory(ams->data_org[AZ_name]); 
              
              ams->options[AZ_pre_calc] = AZ_calc;
              
            }
          else if ( strcmp(Matrix_Factorization_Reuse, "recalc") == 0 )
            {
              ams->options[AZ_pre_calc] = AZ_recalc;
            }
          else if ( strcmp(Matrix_Factorization_Reuse, "reuse") == 0 )
            {
              ams->options[AZ_pre_calc] = AZ_reuse;
            }
          else
            {
              EH(-1, "Unknown factorization reuse specification.");
            }

      /*}*/

          vzero(numUnks, &x[0]);
          linear_solver_blk     = 0; /* count calls to AZ_solve() */
          num_linear_solve_blks = 1; /* upper limit to AZ_solve() calls */
          linear_solver_itns    = 0; /* cumulative number of iterations */
          matrix_solved         = FALSE; 
          while ( ( ! matrix_solved                            ) && 
                  ( linear_solver_blk < num_linear_solve_blks  ) )
            {
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

	      /* Solve the matrix */
              AZ_solve(x, xr, ams->options, ams->params, 
                       ams->indx, ams->bindx, ams->rpntr, ams->cpntr, 
                       ams->bpntr, ams->val, ams->data_org, ams->status, 
                       ams->proc_config);

              first_linear_solver_call = FALSE;

              if ( Debug_Flag > 0 )
                {
                  dump_aztec_status(ams->status);
                }
              
              why = (int)ams->status[AZ_why];
              error = ( why == AZ_normal ? 0 : -1);
              aztec_stringer(why, ams->status[AZ_its], &stringer[0]);
              
              matrix_solved = ( ams->status[AZ_why] == AZ_normal) ;
              linear_solver_blk++;
              linear_solver_itns += ams->status[AZ_its];
              
            } /* End of while loop */

            /* FREE the memory used in storing preconditioner info
             *   - unless using the RE_USE option */

          if (ams->options[AZ_pre_calc] == AZ_calc) {
            AZ_free_memory(ams->data_org[AZ_name]); 
	  }

	  passdown.num_linear_its += linear_solver_itns;

          break;
          
        case AMESOS:

             if( strcmp( Matrix_Format,"msr" ) == 0 ) {
                 amesos_solve_msr( Amesos_Package, ams, x, xr, 1 , pg->imtrx);
             } else if ( strcmp( Matrix_Format,"epetra" ) == 0 ) {
                 amesos_solve_epetra(Amesos_Package, ams, x, xr, pg->imtrx);
             } else {
                 EH(-1," Sorry, only MSR and Epetra matrix formats are currently supported with the Amesos solver suite\n");
             }
        strcpy(stringer, " 1 ");
        break;

        case MA28:
          /*
           * sl_ma28 keeps interntal static variables to determine whether
           * it is the first call or not.
           */
#ifdef HARWELL    
          error = cmsr_ma28 (NumUnknowns[pg->imtrx], NZeros, a, ija, x, xr);
#endif
#ifndef HARWELL
          EH(-1, "That linear solver package is not implemented.");
#endif
          strcpy(stringer, " 1 ");
          break;

        case FRONT:

          if (Num_Proc > 1) EH(-1, "Whoa.  No front allowed with nproc>1");
#ifdef HAVE_FRONT  

/* Initialize frontal solver arguments */

          mf_resolve = 1;
          smallpiv = 1.e-5;
          singpiv = 1.e-14;
          iautopiv = 1;
          iscale = 1;   /* This routine handles resolves only! */
          scaling_max = 1.0;

          /* get global element size and velocity norm if needed for PSPG or Cont_GLS */
          if((PSPG && Num_Var_In_Type[pg->imtrx][PRESSURE]) || (Cont_GLS && Num_Var_In_Type[pg->imtrx][VELOCITY1]))
            {
              h_elem_avg = global_h_elem_siz(x,
					     passdown.x_old,
					     passdown.xdot,
					     passdown.resid_vector,
					     passdown.exo,
					     passdown.dpi);
              U_norm     = global_velocity_norm(x,
						passdown.exo,
						passdown.dpi);
            }
          else
            {
              h_elem_avg = 0.;
              U_norm     = 0.;
            }

            error = mf_solve_lineqn(&mf_resolve, /* re_solve                 */
                                    xr, /* rhs                               */
                                    1, /* nrhs                               */
                                    fss->ncod, /* nsetbc                      */
                                    fss->bc, /* bcvalue                       */
                                    &smallpiv, /* smallpiv                   */
                                    &singpiv, /* singpiv                     */
                                    &iautopiv, /* iautopiv                   */
                                    &iscale, /* iscale                       */
                                    matrix_fill, /* element matrix fill fnc  */
                                    x, /* lhs                                */
        /* This list of args */     &scaling_max, /* scaling max             */
        /* below matrix_fill */     passdown.scale,
        /* pointer is the arg*/     passdown.ams,
        /* list for matrix_fill */  x,
        /* If you change that*/     passdown.resid_vector,
        /* arglist, you must */     passdown.x_old,
        /* change the frontal */    passdown.x_older,
        /* solver.            */    passdown.xdot,
                                    passdown.xdot_old,
                                    passdown.x_update,
                                    &(passdown.delta_t),
                                    &(passdown.theta),
                                    First_Elem_Side_BC_Array,
                                    &(passdown.time_value),
                                    passdown.exo,
                                    passdown.dpi,
                                    &(passdown.dpi->num_universe_nodes),
                                    &h_elem_avg,
                                    &U_norm);
#endif
          break;

        default:
          EH(-1, "That linear solver package is not implemented.");
          break;
        }

/* Report solve time and status */
    s_end = ut();
    if (ProcID == 0)
      {
        if (Linear_Solver == AZTEC)
          {
            if (why == AZ_normal)
              {
                printf("\tResolve time = %7.1e   lits = %s\n",
                  (s_end - s_start), stringer);
              }
            else
              {
                printf("WARNING:  Aztec status was %s !\n", stringer);
              }
          }
        else
          {
            printf(" Resolve_time:%7.1e ", (s_end - s_start) );
          }
      }
      
/* Backup old solutions if not previously done */
  if (passdown.method == FIRST_ORDER_CONTINUATION)
    {
      dcopy1(NumUnknowns[pg->imtrx], passdown.x_older, passdown.x_oldest);
      dcopy1(NumUnknowns[pg->imtrx], passdown.x_old, passdown.x_older);
      dcopy1(NumUnknowns[pg->imtrx], passdown.x, passdown.x_old);
      dcopy1(NumUnknowns[pg->imtrx], passdown.x_sens_temp, passdown.x_sens);
    }

  safe_free( (void *) xr);
  return error;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int komplex_linear_solver_conwrap(double *c, double *d, 
                                  int jac_flag, double *omega, double *tmp)
/* Put the call to your komplex linear solver here. There are three options
 * about the reuse of the preconditioner. It is always safe to
 * just solve the matrix from scratch.
 * Input:
 *    c          Right hand side of real part
 *    d          Right hand side of imaginary part
 *    jac_flag   Flag indicating the status of the Jacobian so that
 *               preconditioners can be used:
 *               NEW_JACOBIAN:   recalculate preconditioner
 *               OLD_JACOBIAN:   reuse preconditioner
 *               OLD_JACOBIAN_DESTROY:   reuse preconditioner,
 *                                       then destroy the preconditioner
 *
 * Output:
 *    c, d       Solution vectors
 *
 * Return Value:
 *    Negative value means linear solver didn't converge.
 */
{
#ifdef KOMPLEX

  /* Declare Variables */
  struct Aztec_Linear_Solver_System *ams = &(passdown.ams[JAC]);
  int numUnks      = NumUnknowns[pg->imtrx] + NumExtUnknowns[pg->imtrx];
  int i;
  int tmp_pre_calc;
  int why;
  dbl           lits;           /* number of linear solver iterations taken */
  char          stringer[80];   /* holding format of num linear solve itns */
  dbl           s_start;        /* mark start of solve */
  dbl           s_end;          /* mark end of solve */
  double *x=NULL, *rhs=NULL;    /* Komplex initial guess and R.H.S. */
  double *a=NULL, *b=NULL;      /* Solution vectors */
  double *wM=NULL;              /* omega*Mass Matrix */
  AZ_MATRIX  *Kmat=NULL;        /* Komplex matrix to be solved. */
  AZ_PRECOND *Prec=NULL;        /* Preconditioner for Komplex */
  AZ_MATRIX *J = passdown.amat; /* Jacobian Matrix */

  /* Allocate and initialize variables */
  wM = (double *) array_alloc(1, ams->nnz+5, sizeof(double));
  a  = (double *) array_alloc(1, numUnks, sizeof(double));
  b  = (double *) array_alloc(1, numUnks, sizeof(double));
  init_vec_value(a, 0.0, numUnks);
  init_vec_value(b, 0.0, numUnks);
/*init_vec_value(tmp, 0.0, numUnks);*/

  /* Reuse preconditioning option  (commented out for now): */
  tmp_pre_calc = passdown.options[AZ_pre_calc];
  /*
  if ((jac_flag == OLD_JACOBIAN)||(jac_flag == OLD_JACOBIAN_DESTROY)) {
    passdown.options[AZ_pre_calc] = AZ_reuse;
  } else if (jac_flag == NEW_JACOBIAN){
    passdown.options[AZ_pre_calc] = AZ_calc;
  } else {
    printf("ERROR: Komplex linear solve conwrap: unknown precond flag=%d\n"
           ,jac_flag);
    exit(-1);
  }
  */
  passdown.options[AZ_pre_calc] = AZ_calc;

  /* Apply scaling to the entire system */
  row_sum_scaling_scale(ams, c, passdown.scale);
  matrix_scaling(ams, passdown.mass_matrix, 1.0, passdown.scale);
  vector_scaling(NumUnknowns[pg->imtrx], d, passdown.scale);

  /*
   * Construct Elements of the Komplex Matrix:
   *
   *    J   wM   a     c
   *                 =
   *   -wM   J   b     d
   */
  s_start = ut();
  for (i=0 ; i<ams->nnz; i++)
    wM[i] = passdown.mass_matrix[i] * (*omega) * (-1.0);

  /* Create the Komplex matrix */
  AZK_create_linsys_ri2k(a, b, c, d, ams->options, ams->params,
                         ams->proc_config, J, wM, &x, &rhs, &Kmat);

  /* Don't overwrite the scaling factors from the "a" vector calculation! */
  Kmat->data_org[AZ_name] = 2;

  /* Create Komplex matrix preconditioner */
  AZK_create_precon(ams->options, ams->params,
                    ams->proc_config, x, rhs, Kmat, &Prec);

  /* Solve the linear system */
  AZ_iterate(x, rhs, ams->options, ams->params, ams->status,
             ams->proc_config, Kmat, Prec, NULL);
/*
  if (passdown.status[AZ_why] != AZ_normal)
    {
      DPRINTF(stderr, "######Warning: linear solver did not converge\n");
    }
*/
  /* Report solve time and status */
  s_end = ut();
  why  = (int)ams->status[AZ_why];
  lits = ams->status[AZ_its];
  aztec_stringer(why, lits, &stringer[0]);

  if (ProcID == 0)
    {
      if (why == AZ_normal)
        {
          printf("\tKomplex solve time = %7.1e   lits = %s\n",
            (s_end - s_start), stringer);
        }
      else
        {
          printf("WARNING:  Aztec Komplex status was %s !\n", stringer);
        }
    }


  /* Break solution into real (c) and imaginary (d) parts. */
  AZK_extract_solution_k2ri(ams->options, ams->params,
                            ams->proc_config, Kmat, Prec, x, c, d);

  /* Clean up Memory */
  AZ_free_memory (Kmat->data_org[AZ_name]);
  AZK_destroy_precon (ams->options, ams->params, ams->proc_config, Kmat, &Prec);  AZK_destroy_linsys (ams->options, ams->params,
                      ams->proc_config, &x, &rhs, &Kmat);
  safe_free((void *) wM);
  safe_free((void *) a);
  safe_free((void *) b);
  passdown.options[AZ_pre_calc] = tmp_pre_calc;

#endif
  return 0;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void matrix_residual_fill_conwrap(double *x, double *rhs, int matflag)
/* Put the call to your matrix/residual fill routine here.
 * Input:
 *    x         Solution vector
 *    matflag   Flag indicating residual (RHS_ONLY), matrix (MATRIX_ONLY),
 *              or both (RHS_MATRIX) are requested.
 *		Also, saves unscaled Jacobian (SAVE_UNSCALED_MATRIX)
 *		or recovers it (RECOVER_UNSCALED_MATRIX).
 *
 * Output:
 *    rhs       Right hand side
 *
 * Return Value:
 *
 * Modified for GOMA matrix_fill function - EDW 6/1/2000 
 */
{

/* Get first and last elements on this processor */
  struct Aztec_Linear_Solver_System *ams = &(passdown.ams[JAC]);
  int save_flag = FALSE;
  double h_elem_avg, U_norm;

/* For eigensolver, Jacobian matrix is assembled along with mass matrix */
  if (passdown.LSA_flag) return;
  af->Assemble_LSA_Jacobian_Matrix = FALSE;
  af->Assemble_LSA_Mass_Matrix = FALSE;

/* If this call is for the first resolve, set save_flag */
  if (matflag == RHS_MATRIX_SAVE)
    {
      if (passdown.method != ARC_LENGTH_CONTINUATION) save_flag = TRUE;
      matflag = RHS_MATRIX;
    }

/* If using frontal solver:  just use old matrix */
  if (Linear_Solver == FRONT)
    {
       matflag = RHS_ONLY;
       save_flag = FALSE;
    }

/* Get global element size and velocity norm if needed for PSPG or Cont_GLS */
  if((PSPG && Num_Var_In_Type[pg->imtrx][PRESSURE]) || (Cont_GLS && Num_Var_In_Type[pg->imtrx][VELOCITY1]))
    {
      h_elem_avg = global_h_elem_siz(x,
				     passdown.x_old,
				     passdown.xdot,
				     passdown.resid_vector,
				     passdown.exo,
				     passdown.dpi);
      U_norm     = global_velocity_norm(x,
					passdown.exo,
					passdown.dpi);
    }
  else
    {
      h_elem_avg = 0.;
      U_norm     = 0.;
    }

  switch (matflag) {
    case RHS_ONLY:
      init_vec_value(passdown.resid_vector, 0.0, NumUnknowns[pg->imtrx]);
      af->Assemble_Residual = TRUE;
      af->Assemble_Jacobian = FALSE;
      passdown.num_res_fills++;
      break;
    case MATRIX_ONLY:
      init_vec_value(ams->val, 0.0, NZeros);
      af->Assemble_Residual = FALSE;
      af->Assemble_Jacobian = TRUE;
      passdown.num_mat_fills++;
      break;
    case RHS_MATRIX:
      init_vec_value(passdown.resid_vector, 0.0, NumUnknowns[pg->imtrx]);
      init_vec_value(ams->val, 0.0, NZeros);
      af->Assemble_Residual = TRUE;
      af->Assemble_Jacobian = TRUE;
      passdown.num_res_fills++;
      passdown.num_mat_fills++;
      break;
  }
  af->Assemble_LSA_Mass_Matrix = FALSE;

/* Recover old Jacobian only */

  if(matflag == RECOVER_MATRIX) dcopy1(NZeros, ams->val_old, ams->val);

/* Remaining cases: perform requested fill */
  else {

    exchange_dof(passdown.cx, passdown.dpi, x, pg->imtrx);

    (void) matrix_fill_full(ams,
			    x,
			    rhs,
			    passdown.x_old,
			    passdown.x_older,
			    passdown.xdot,
			    passdown.xdot_old,
			    passdown.x_update,
			    &(passdown.delta_t),
			    &(passdown.theta),
			    First_Elem_Side_BC_Array[pg->imtrx],
			    &(passdown.time_value),
			    passdown.exo,
			    passdown.dpi,
			    &(passdown.dpi->num_universe_nodes),
			    &(h_elem_avg),
			    &(U_norm),
                            NULL);
  }

/* Save unscaled matrix before the first resolve */
  if (save_flag) dcopy1(NZeros, ams->val, ams->val_old);

}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void mass_matrix_fill_conwrap(double *x, double *rhs)
/* Put the call to your matrix/residual fill routine here.
 * Input:
 *    x         Solution vector
 *    rhs       Right hand side
 *
 * Output:
 *    Creates mass matrix M
 *
 * Return Value:
 */
{
  struct Aztec_Linear_Solver_System *ams = &(passdown.ams[JAC]);
  int i, j, mn, nnodes;
  int passes;
  int *ija = ams->bindx;
  int **e_save=NULL;
  double ***etm_save=NULL;
  double *scale=NULL, *zero=NULL;
  double *jacobian_matrix=NULL;
  double *tmp_matrix=NULL;
  double theta_save;
  double h_elem_avg, U_norm;

/* Initialize arrays */
  zero = (double *) array_alloc(1, NumUnknowns[pg->imtrx], sizeof(double));
  init_vec_value(&zero[0], 0.0, NumUnknowns[pg->imtrx]);

  scale = (double *) array_alloc(1, NumUnknowns[pg->imtrx], sizeof(double));
  init_vec_value(scale, 1.0, NumUnknowns[pg->imtrx]);

/* Initialize tmp_matrix array for 3D of 2D stability if needed */
  if (passdown.do_3D_of_2D)
    {
      passes = 2;
      tmp_matrix = (double *) array_alloc(1, NZeros+1, sizeof(double));
    }
  else passes = 1;

/* Resolve difference between x and x_old to get correct mass matrix terms */
  exchange_dof(passdown.cx, passdown.dpi, x, pg->imtrx);
  dcopy1(NumUnknowns[pg->imtrx], passdown.x, passdown.x_old);
  dcopy1(NumUnknowns[pg->imtrx], passdown.x, passdown.x_older);

  nnodes = passdown.dpi->num_universe_nodes;

/* Save original e and etm arrays */
  theta_save = passdown.theta;
  e_save = (int **)array_alloc(2, upd->Num_Mat, MAX_EQNS, sizeof(int));
  etm_save = (dbl ***)array_alloc(3, upd->Num_Mat, MAX_EQNS,
                                  MAX_TERM_TYPES, sizeof(dbl));
  for(mn = 0; mn < upd->Num_Mat; mn++)
    for(i = 0; i < MAX_EQNS; i++)
      {
        e_save[mn][i] = pd_glob[mn]->e[pg->imtrx][i];
        for(j = 0; j < MAX_TERM_TYPES; j++)
          etm_save[mn][i][j] = pd_glob[mn]->etm[pg->imtrx][i][j];
      }

/* Get global element size and velocity norm if needed for PSPG or Cont_GLS */
  if((PSPG && Num_Var_In_Type[pg->imtrx][PRESSURE]) || (Cont_GLS && Num_Var_In_Type[pg->imtrx][VELOCITY1]))
    {
      h_elem_avg = global_h_elem_siz(x,
                                     passdown.x_old,
                                     passdown.xdot,
                                     passdown.resid_vector,
                                     passdown.exo,
                                     passdown.dpi);
      U_norm     = global_velocity_norm(x,
                                        passdown.exo,
                                        passdown.dpi);
    }
  else
    {
      h_elem_avg = 0.;
      U_norm     = 0.;
    }

  /* Get jacobian matrix - use space allocated for ams->val */
  jacobian_matrix = ams->val;
  
  /* Fill the LSA Jacobian (in 2 passes of using 3D of 2D LSA) */
  for (i=0; i<passes; i++)
    {
  
  /* On first 3D of 2D pass, route matrix fill to tmp_matrix */
      LSA_3D_of_2D_pass = ( (passes == 2) ? i+1 : 0);
      if (LSA_3D_of_2D_pass == 1) ams->val = tmp_matrix;
      else ams->val = jacobian_matrix;
      if (passdown.do_3D_of_2D)
        {
          DPRINTF (stdout, "Assembling LSA Jacobian pass %d ...\n", i+1);
        }
      else if (passdown.LSA_flag)
        {
          DPRINTF (stdout, "Assembling LSA Jacobian ...\n");
        }

  for (j=0; j<NZeros+1; j++) ams->val[j] = 0.0;
  af->Assemble_Residual = TRUE;
  af->Assemble_Jacobian = TRUE;
  af->Assemble_LSA_Jacobian_Matrix = TRUE;
  af->Assemble_LSA_Mass_Matrix = FALSE;

  exchange_dof(passdown.cx, passdown.dpi, x, pg->imtrx);

  (void) matrix_fill_full(ams,
                          x,
                          rhs,
                          passdown.x_old,
                          passdown.x_older,
                          passdown.xdot,
                          passdown.xdot_old,
                          passdown.x_update,
                          &(passdown.delta_t),
                          &(passdown.theta),
                          First_Elem_Side_BC_Array[pg->imtrx],
                          &(passdown.time_value),
                          passdown.exo,
                          passdown.dpi,
                          &(nnodes),
                          &(h_elem_avg),
                          &(U_norm),
                          NULL);
 

    }  /* End of Jacobian pass loop */

  /* Assemble the two passes if doing 3D of 2D */
  if (passdown.do_3D_of_2D)
    {
      for (j=0; j<NZeros+1; j++) jacobian_matrix[j] += tmp_matrix[j];
    }

  /* This call will fill in scale[] with the proper row scales. */
  if (passdown.LSA_flag)
    {
      row_sum_scaling_scale(ams, passdown.resid_vector, scale);
    }

  /* Get mass matrix */

  /* Set up e and etm arrays for mass matrix */
  /* theta = 0 makes the method implicit */
  passdown.theta = 0.0;
  TimeIntegration = TRANSIENT;
  for(mn = 0; mn < upd->Num_Mat; mn++) 
    {
      pd_glob[mn]->TimeIntegration = TRANSIENT;
      for(i = 0; i < MAX_EQNS; i++)
        if(pd_glob[mn]->e[pg->imtrx][i])
          {
            pd_glob[mn]->e[pg->imtrx][i] = T_MASS;
            for (j=0; j<MAX_TERM_TYPES; j++) pd_glob[mn]->etm[pg->imtrx][i][j] = 0.0;
            pd_glob[mn]->etm[pg->imtrx][i][LOG2_MASS] = 1.0;
          }
    }

  /* Fill the LSA Mass matrix (in 2 passes of using 3D of 2D LSA) */
  for (i=0; i<passes; i++)
    {
  
  /* On first 3D of 2D pass, route matrix fill to tmp_matrix */
      LSA_3D_of_2D_pass = ( (passes == 2) ? i+1 : 0);
      if (LSA_3D_of_2D_pass == 1) ams->val = tmp_matrix;
      else ams->val = passdown.mass_matrix;
      if (passdown.do_3D_of_2D)
        {
          DPRINTF (stdout, "Assembling LSA Mass matrix pass %d ...\n", i+1);
        }
      else if (passdown.LSA_flag)
        {
          DPRINTF (stdout, "Assembling LSA Mass matrix ...\n");
        }

  for(j=0; j<NZeros+1; j++) ams->val[j] = 0.0;
  af->Assemble_Residual = TRUE;
  af->Assemble_Jacobian = TRUE;
  af->Assemble_LSA_Jacobian_Matrix = FALSE;
  af->Assemble_LSA_Mass_Matrix = TRUE;
  passdown.delta_t = 1.0;
  tran->delta_t = 1.0;      /*for Newmark-Beta terms in Lagrangian Solid*/
  
  (void) matrix_fill_full(ams,
                          x,
                          rhs,
                          passdown.x_old,
                          passdown.x_older,
                          passdown.xdot,
                          passdown.xdot_old,
                          passdown.x_update,
                          &(passdown.delta_t),
                          &(passdown.theta),
                          First_Elem_Side_BC_Array[pg->imtrx],
                          &(passdown.time_value),
                          passdown.exo,
                          passdown.dpi,
                          &(nnodes),
                          &(h_elem_avg),
                          &(U_norm),
                          NULL);


    }  /* End of mass matrix pass loop */

  /* Assemble the two passes if doing 3D of 2D */
  if (passdown.do_3D_of_2D)
    {
      for (j=0; j<NZeros+1; j++) passdown.mass_matrix[j] += tmp_matrix[j];
    }
  ams->val = jacobian_matrix;

  /* Scale the mass matrix the same as the Jacobian is scaled.
     Also, change the signs */
  if (passdown.LSA_flag)
    {
      matrix_scaling(ams, passdown.mass_matrix, 1.0, scale);
    }
  for (j=0; j<NZeros+1; j++) passdown.mass_matrix[j] *= -1.0;
  dcopy1(NumUnknowns[pg->imtrx], scale, passdown.scale);

  /* Restore original e and etm values */
  TimeIntegration = STEADY;
  passdown.theta = theta_save;
  for(mn = 0; mn < upd->Num_Mat; mn++)
    for(i = 0; i < MAX_EQNS; i++)
      {
        pd_glob[mn]->TimeIntegration = STEADY;
        pd_glob[mn]->e[pg->imtrx][i] = e_save[mn][i];
        for(j = 0; j < MAX_TERM_TYPES; j++)
          pd_glob[mn]->etm[pg->imtrx][i][j] = etm_save[mn][i][j];
      }

  /* Print matrices -- Careful, these can be big disk space hogs!! */
  if(eigen->Eigen_Matrix_Output)
    {
      if (!passdown.do_3D_of_2D) LSA_3D_of_2D_wave_number = -1.0;
      output_stability_matrices(passdown.mass_matrix, jacobian_matrix, ija,
                                nnodes, NumUnknowns[pg->imtrx], NZeros);
    }

  /* Clean up */
  safe_free ( (void *) scale);
  safe_free ( (void *) zero);
  safe_free ( (void *) tmp_matrix);
  safe_free ( (void *) e_save);
  safe_free ( (void *) etm_save);

  return;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void matvec_mult_conwrap(double *x, double *y)
/* Put the call to your matrix-vector multiply here.
 * Input:
 *    x         Vector of length number of unknowns
 *
 * Output:
 *    y         Matrix times x.
 *
 * Return Value:
 */

/* NOTE: Aztec provides a matrix-vector multiply routine for VBR-format matrices.
 *       For now, users without AZTEC will have to use the MSR format when doing
 *       continuation with LOCA - EDW 6/1/2000.
 */

{
  struct Aztec_Linear_Solver_System *ams = &(passdown.ams[JAC]);
  register int j, k, irow, bindx_row;
  int          N, nzeros;

  /* exchange boundary info */
  exchange_dof(passdown.cx, passdown.dpi, x, pg->imtrx);

/* First, handle MSR matrices */ 
  if( strcmp( Matrix_Format, "msr" ) == 0)

/* Note: This routine borrowed from "az_matvec_mult.c"! - EDW */
    {

      N = ams->data_org[AZ_N_internal] + ams->data_org[AZ_N_border];

      for (irow = 0; irow < N; irow++) {

    /* compute diagonal contribution */

        *y = ams->val[irow] * x[irow];

    /* nonzero off diagonal contribution */

        bindx_row = ams->bindx[irow];
        nzeros    = ams->bindx[irow+1] - bindx_row;

        for (j = 0; j < nzeros; j++) {
          k   = bindx_row + j;
          *y += ams->val[k] * x[ams->bindx[k]];
        }
        y++;
      }
    }

  else if (strcmp(Matrix_Format, "vbr") == 0)
    {

/* For VBR matrices, use Aztec matvec mult routine */
      AZ_VBR_matvec_mult(x, y, passdown.amat, passdown.proc_config);
    }

/* Error message if not MSR or VBR */
  else
    {
      EH(-1, "Matrix format must be MSR or VBR!");
    }

}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void mass_matvec_mult_conwrap(double *x, double *y)
/* Put the call to your mass matrix-vector multiply here.
 * Input:
 *    x         Vector of length number of unknowns
 *
 * Output:
 *    y         Mass matrix times x.
 *
 * Return Value:
 */

/*
 * NOTE: Aztec provides a matrix-vector multiply routine for VBR-format matrices.
 *       For now, users without AZTEC will have to use the MSR format when
 *       doing continuation with LOCA - EDW.
 */
{
  struct Aztec_Linear_Solver_System *ams = &(passdown.ams[JAC]);
  double *m = passdown.mass_matrix;
  register int j, k, irow, bindx_row;
  int          N,  nzeros;

  /* exchange boundary info */
  exchange_dof(passdown.cx, passdown.dpi, x, pg->imtrx);

/* First, handle MSR matrices */ 
  if( strcmp( Matrix_Format, "msr" ) == 0)

/* Note: This routine borrowed from "az_matvec_mult.c"! - EDW */
    {

      N = ams->data_org[AZ_N_internal] + ams->data_org[AZ_N_border];

      for (irow = 0; irow < N; irow++) {

    /* compute diagonal contribution */

        *y = m[irow] * x[irow];

    /* nonzero off diagonal contribution */

        bindx_row = ams->bindx[irow];
        nzeros    = ams->bindx[irow+1] - bindx_row;

        for (j = 0; j < nzeros; j++) {
          k   = bindx_row + j;
          *y += m[k] * x[ams->bindx[k]];
        }
        y++;
      }
    }

  else if (strcmp(Matrix_Format, "vbr") == 0)
    {

/* For VBR matrices, use Aztec matvec mult routine */
      AZ_VBR_matvec_mult(x, y, passdown.mmat, passdown.proc_config);
    }

/* Error message if not MSR or VBR */
  else
    {
      EH(-1, "Matrix format must be MSR or VBR!");
    }

  return;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void create_shifted_matrix_conwrap(void)
/* Allocates a sets sparsity pointers for shifted matrix.
 * Only used by eigensolver
 */
{
#ifdef HAVE_ARPACK
  passdown.shifted_matrix = (double *)
                             array_alloc(1, NZeros+5, sizeof(double));
  init_vec_value(passdown.shifted_matrix, 0.0, NZeros+5);
#endif
  return;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void shifted_matrix_fill_conwrap(double sigma)
/* Routine to created shifted matrix  J-sigma M
 * Only used by eigensolver
 */
{
#ifdef HAVE_ARPACK
  struct Aztec_Linear_Solver_System *ams = &(passdown.ams[JAC]);
  int i;

  for(i=0; i<NZeros+1; i++)
    passdown.shifted_matrix[i] = ams->val[i] - sigma * passdown.mass_matrix[i];
#endif

  return;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void shifted_linear_solver_conwrap(double *x, double *y, 
                                   int jac_flag, double tol)
{
#ifdef HAVE_ARPACK
  struct Aztec_Linear_Solver_System *ams = &(passdown.ams[JAC]);
  static int first_linear_solver_call=TRUE;
  static int stab_umf_id = -1;
  double *a;                    /* nonzero values of a CMSR matrix */
  int    *ija = ams->bindx;     /* column pointer array into matrix "a" */
  static int    Factor_Flag;    /* For UMFPACK */
  int           matr_form;      /* 1: MSR FORMAT MATRIX FOR UMFPACK DRIVER */
  dbl           lits;           /* number of linear solver iterations taken */
  char          stringer[80];   /* holding format of num linear solve itns */
  int   linear_solver_blk;      /* count calls to AZ_solve() */
  int   linear_solver_itns;     /* count cumulative linearsolver iterations */
  int   num_linear_solve_blks;  /* one pass for now */
  int   matrix_solved = 0;      /* boolean */
  int tmp_scale, tmp_conv;
  dbl tmp_tol;
  dbl *tmp_a;

/* Allocate a separate system for stability if using UMFPACK */
  if (first_linear_solver_call)
    {
      if (Linear_Solver != UMFPACK2 && Linear_Solver != UMFPACK2F)
        first_linear_solver_call = FALSE;
    }

/* Save some previous settings, point ams->val to the shifted matrix */
  tmp_a = ams->val;
  tmp_scale = ams->options[AZ_scaling];
  tmp_conv  = ams->options[AZ_conv];
  tmp_tol   = ams->params[AZ_tol];

/* Reset these for the eigensolver */
  ams->val = passdown.shifted_matrix;
  a = ams->val;
  ams->options[AZ_scaling] = AZ_none;
  ams->options[AZ_conv] = AZ_noscaled;
  ams->params[AZ_tol] = tol;

/* Proceed with the chosen linear solver */
  switch (Linear_Solver)
    {
    case UMFPACK2:
    case UMFPACK2F:
      if (strcmp(Matrix_Format, "msr"))
	EH(-1,"ERROR: umfpack solver needs msr matrix format");

      Factor_Flag = ( (jac_flag == NEW_JACOBIAN) ? 0 : 3);
      if (Linear_Solver == UMFPACK2F) Factor_Flag = 0;
      /*  */
      matr_form = 1;
      stab_umf_id = SL_UMF(stab_umf_id,
			   &first_linear_solver_call, &Factor_Flag, &matr_form, 
             &NumUnknowns[pg->imtrx], &NZeros, &ija[0], &ija[0], &a[0],
			   &x[0], &y[0]);
      /*  */
      first_linear_solver_call = FALSE;
      strcpy(stringer, " 1 ");
      break;

    case SPARSE13a:
      if (strcmp(Matrix_Format, "msr"))
	EH(-1,"ERROR: lu solver needs msr matrix format");

        dcopy1(NumUnknowns[pg->imtrx], x, y);
        lu(NumUnknowns[pg->imtrx], NumExtUnknowns[pg->imtrx], NZeros, a, ija, y, 2);
      first_linear_solver_call = FALSE;
      /* 
       * Note that sl_lu has static variables to keep track of
       * first call or not.
       */

      strcpy(stringer, " 1 ");
      break;
    case AMESOS:
      if( strcmp( Matrix_Format,"msr" ) == 0 ) {
	amesos_solve_msr( Amesos_Package, ams, y, x, 1, pg->imtrx );
      } else {
	EH(-1," Sorry, only MSR  matrix format supported for loca eigenvalue");
      }
      first_linear_solver_call = FALSE;
      strcpy(stringer, " 1 ");
      break;
	

    case AZTEC:

      /* Set option of preconditioner reuse */
            
      if ( jac_flag == OLD_JACOBIAN )
        {
          ams->options[AZ_pre_calc] = AZ_reuse;
        }
      else
        {
          if ( strcmp(Matrix_Factorization_Reuse, "calc") == 0 )
            {
	      /*
	       * Gonna start from scratch even though I've cooked a
	       * preconditioner in the kitchen all day? Well, then
	       * you won't need the leftover pieces from all my
	       * hard preparation last time around.
	       */
 
              AZ_free_memory(ams->data_org[AZ_name]); 
              
              ams->options[AZ_pre_calc] = AZ_calc;
              
            }
          else if ( strcmp(Matrix_Factorization_Reuse, "recalc") == 0 )
            {
              ams->options[AZ_pre_calc] = AZ_recalc;
            }
          else if ( strcmp(Matrix_Factorization_Reuse, "reuse") == 0 )
            {
              ams->options[AZ_pre_calc] = AZ_reuse;
            }
          else
            {
              EH(-1, "Unknown factorization reuse specification.");
            }
        }

      linear_solver_blk     = 0; /* count calls to AZ_solve() */
      num_linear_solve_blks = 1; /* upper limit to AZ_solve() calls */
      linear_solver_itns    = 0; /* cumulative number of iterations */
      matrix_solved         = FALSE; 
      while ( ( ! matrix_solved                            ) && 
              ( linear_solver_blk < num_linear_solve_blks  ) )
        {
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

	  /* Solve the matrix */
          AZ_solve(y, x, ams->options, ams->params, 
                   ams->indx, ams->bindx, ams->rpntr, ams->cpntr, 
                   ams->bpntr, ams->val, ams->data_org, ams->status, 
                   ams->proc_config);

          first_linear_solver_call = FALSE;

          if ( Debug_Flag > 0 )
            {
              dump_aztec_status(ams->status);
            }
              
          strcpy(stringer, "   ");
          switch ( (int)(ams->status[AZ_why]) )
            {
	    case AZ_normal:
	      lits = ams->status[AZ_its];
	      if ( lits < 1000 )
		{
		  sprintf(stringer, "%3d", (int)lits);
		}
	      else if ( lits < 10000 )
		{
		  sprintf(stringer, "%2dh", (int)(lits/1e2));
		}
	      else if ( lits < 100000 )
		{
		  sprintf(stringer, "%2dk", (int)(lits/1e3));
		}
	      else
		{
		  sprintf(stringer, "%3.0e", lits);
		}
	      break;
	    case AZ_param:
	      strcpy(stringer, "bad");
	      break;
	    case AZ_breakdown:
	      strcpy(stringer, "brk");
	      break;
	    case AZ_loss:
	      strcpy(stringer, "los");
	      break;
	    case AZ_maxits:
	      strcpy(stringer, "max");
	      break;
	    case AZ_ill_cond:
	      strcpy(stringer, "ill");
	      break;
	    default:
	      strcpy(stringer, "???");
	      break;
            }
              
          matrix_solved = ( ams->status[AZ_why] == AZ_normal) ;
          linear_solver_blk++;
          linear_solver_itns += ams->status[AZ_its];
              
        } /* End of while loop */

      passdown.num_eigen_its += linear_solver_itns;
      break;

    case MA28:
      /*
       * sl_ma28 keeps interntal static variables to determine whether
       * it is the first call or not.
       */
#ifdef HARWELL    
        err = cmsr_ma28 (NumUnknowns[pg->imtrx], NZeros, a, ija, y, x);
#endif
#ifndef HARWELL
      EH(-1, "That linear solver package is not implemented.");
#endif
      strcpy(stringer, " 1 ");
      break;

    case FRONT:
      /* Frontal solver cannot be used for eigensolves! */
      EH(-1, "Frontal solver cannot be used for eigensolves!");
      break;

    default:
      EH(-1, "That linear solver package is not implemented for eigensolves");
      break;
    }

/* Restore original settings for the next step */
  ams->val = tmp_a;
  ams->options[AZ_scaling] = tmp_scale;
  ams->options[AZ_conv] = tmp_conv;
  ams->params[AZ_tol] = tmp_tol;
#endif

/* Done */
  return;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void destroy_shifted_matrix_conwrap(void)
{
/* Just deallocate the matrix */
#ifdef HAVE_ARPACK
  safe_free ( (void *) passdown.shifted_matrix);
#endif
  return;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void assign_parameter_conwrap(double param)
/* Put the call to a routine to assign the continuation parameter here.
 * Input:
 *    param     New value of continuation parameter.
 *
 * Output:
 *
 * Return Value:
 */
{
  int iCC;
  double delta, theta;

  if (passdown.method == LOCA_LSA_ONLY) return;

  delta = param - cont->BegParameterValue;

  /*
   * Angular continuation parameters are handled differently
   * Here, param is the continuation angle (in degrees)
   */
  if (cont->upType == 6)  /* AN */
    {
      theta = param * M_PIE / 180.0;
      for (iCC=1; iCC<nCC; iCC++)
        {
          switch (cpcc[iCC].fn_flag) {
            case 0:
              cpcc[iCC].value = cpcc[iCC].coeff_0
                              + cpcc[iCC].coeff_1 * sin(theta);
              break;
            case 1:
              cpcc[iCC].value = cpcc[iCC].coeff_0
                              + cpcc[iCC].coeff_1 * cos(theta);
              break;
            case 2:
              cpcc[iCC].value = cpcc[iCC].coeff_0
                              + cpcc[iCC].coeff_1 * tan(theta);
              break;
            case 3:
              cpcc[iCC].value = cpcc[iCC].coeff_0
                              + cpcc[iCC].coeff_1 * sin(theta)
                              + cpcc[iCC].coeff_2 * cos(theta);
              break;
            }

          update_parameterC(iCC, cpcc[iCC].value, passdown.x,
                            passdown.xdot, passdown.x_AC, cont->Delta_s0,
		            passdown.cx, passdown.exo, passdown.dpi);
        }
    }
  else
    {
      for (iCC=0; iCC<nCC; iCC++)
        {
  /* This is for the variable-exponent special case */
          if (cpcc[iCC].fn_flag == 3)
            {
              cpcc[iCC].value = cpcc[iCC].coeff_0 + cpcc[iCC].coeff_1
                              * pow(param, cpcc[iCC].coeff_2);
            }
  /* This is for all other continuation types (simple linear relation) */
          else
            {
              cpcc[iCC].value = cpcc[iCC].Beg_CC_Value + cpcc[iCC].ratio * delta;
            }

          update_parameterC(iCC, cpcc[iCC].value, passdown.x,
                            passdown.xdot, passdown.x_AC, cont->Delta_s0,
		            passdown.cx, passdown.exo, passdown.dpi);
        }
    }
  return;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void assign_bif_parameter_conwrap(double bif_param)
/* Put the call to a routine to assign the bifurcation parameter here.
 * Input:
 *    bif_param     New value of continuation parameter.
 *
 * Output:
 *
 * Return Value:
 */
{
  int iTC;
  double delta, theta;

  delta = bif_param - tpcc[0].Beg_CC_Value;

  /*
   * Angular continuation parameters are handled differently
   * Here, bif_param is the continuation angle (in degrees)
   */
  if (loca_in->TPupType == 6)  /* AN */
    {
      theta = bif_param * M_PIE / 180.0;
      for (iTC=1; iTC<nTC; iTC++)
        {
          switch (tpcc[iTC].fn_flag) {
            case 0:
              tpcc[iTC].value = tpcc[iTC].coeff_0
                              + tpcc[iTC].coeff_1 * sin(theta);
              break;
            case 1:
              tpcc[iTC].value = tpcc[iTC].coeff_0
                              + tpcc[iTC].coeff_1 * cos(theta);
              break;
            case 2:
              tpcc[iTC].value = tpcc[iTC].coeff_0
                              + tpcc[iTC].coeff_1 * tan(theta);
              break;
            case 3:
              tpcc[iTC].value = tpcc[iTC].coeff_0
                              + tpcc[iTC].coeff_1 * sin(theta)
                              + tpcc[iTC].coeff_2 * cos(theta);
              break;
          }
                                                                                
          update_parameterTP(iTC, tpcc[iTC].value, passdown.x,
                             passdown.xdot, passdown.x_AC, cont->Delta_s0,
                             passdown.cx, passdown.exo, passdown.dpi);
        }
    }
  /* This is for all other continuation types (simple linear relation) */
  else
    {
      for (iTC=0; iTC<nTC; iTC++)
        {
  /* This is for the variable-exponent special case */
          if (tpcc[iTC].fn_flag == 3)
            {
              tpcc[iTC].value = tpcc[iTC].coeff_0 + tpcc[iTC].coeff_1
                              * pow(bif_param, tpcc[iTC].coeff_2);
            }
  /* This is for all other continuation types (simple linear relation) */
          else
            {
              tpcc[iTC].value = tpcc[iTC].Beg_CC_Value + tpcc[iTC].ratio * delta;
            }

          update_parameterTP(iTC, tpcc[iTC].value, passdown.x,
                             passdown.xdot, passdown.x_AC, cont->Delta_s0,
		             passdown.cx, passdown.exo, passdown.dpi);
        }
    }
  return;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void calc_scale_vec_conwrap(double *x, double *scale_vec, int numUnks)
/* Put the call to a routine to calculate a scaling vector here.
 * Input:
 *    x          New value of continuation parameter.
 *    numUnks    Number of unknowns on this proc, the length of x
 *               and scale_vec.
 *
 * Output:
 *    scale_vec  Vector of length number of unknowns used to scale
 *               variables so that one type of unknown (e.g. pressure)
 *               doesn't dominate over others. Used to balance the
 *               variables and the arc-length variable in arc-length
 *               continuation, and for scaling the null vector in
 *               turning point tracking. Using reciprocal of the average
 *               value of that variable type is a good choice. Vector
 *               of all ones should suffice for most problems.
 *
 * Return Value:
 */
{
  static int sv_init = TRUE;	/* Initialize arrays on first call */
  int p=9, ip1=8;
  int i, iunk, idof, index, ivd;
  double *sv_sum=NULL;		/* Running sum of each var type    */
  VARIABLE_DESCRIPTION_STRUCT *vdi;  /* Ptr to current vd struct   */

  /* Calculate variable averages only if passdown.sv_flag is set */
  if (passdown.sv_flag)
    {
  /* Allocate and zero summing array */

  /*
   * Initialize local arrays on first call.
   * sv_index[j] contains index# of vd structure for unknown j.
   * sv_count[k] tallies number of unknowns having the vd structure
   *		with list_index = k.
   * sv_sum[k]   sums unknown values of type k.
   * NOTE: This may not work if different processors have
   *	 different lists of vd's!
   */

      if (sv_init)
        {

  /* Determine if discontinuous (P1) pressure interpolation is being used */
          for (i=0; i<upd->Num_Mat; i++)
	    {

  /* Adjust nvd to allow for separate treatment of the three pressure vars */
              if (pd_glob[i]->i[pg->imtrx][p] == ip1) passdown.nvd += 2;
	    }

  /* Allocate arrays needed to calculate scale vector */
          passdown.sv_index = (int *) array_alloc(1, NumUnknowns[pg->imtrx], sizeof(int));
          passdown.sv_count = (double *) array_alloc(1, passdown.nvd, sizeof(double));
          init_vec_value(passdown.sv_count, 0.0, passdown.nvd);

  /* Loop over only this processor's owned unknowns */
          for (iunk=0; iunk<NumUnknowns[pg->imtrx]; iunk++)
	    {

  /* Get vd structure for this unknown */
	      vdi = Index_Solution_Inv(iunk, NULL, &ivd, NULL, &idof, pg->imtrx);
	      index = vdi->List_Index;

  /* P1: The second and third pressure vars are stored at the back of the array */
  /* KT 02/19/2016: This logic causes memory error when using P1
                    Commenting these lines do not appear to affect the performance
                    of arc length continuation */
//	      if (vdi->Ndof > 1 && idof > 0)
//                {
//	          index_adj = idof - 1;
//	          if (vdi->MatID > 0) index_adj += 2 * vdi->MatID; /* MatID starts @ 0 */
//	          index = Num_Var_Info_Records + index_adj;
//	        }

  /* Assign sv_index and increment sv_count */
	      passdown.sv_index[iunk] = index;
	      passdown.sv_count[index] += 1.0;
	    }

  /* Broadcast sv_count to all processors */
     /*for (i=0; i<passdown.nvd; i++)
	  passdown.sv_count[i] = gsum_double_conwrap(passdown.sv_count[i]);*/
        }

  /* Reset sv_init */
      sv_init = FALSE;

  /* Allocate and zero summing array */
      sv_sum = (double *) array_alloc(1, passdown.nvd, sizeof(double));
      init_vec_value(sv_sum, 0.0, passdown.nvd);

  /* Calculate the scale vector for the current x */
      for (iunk=0; iunk<NumUnknowns[pg->imtrx]; iunk++)
        {

  /* Add each value to the appropriate index of sv_sum */
          index = passdown.sv_index[iunk];
          sv_sum[index] += fabs(x[iunk]);
        }

  /* Get global averages for each unknown type */
      for (index=0; index<passdown.nvd; index++)
        {
    /*sv_sum[index] =
        gsum_double_conwrap(sv_sum[index]) / passdown.sv_count[index];*/
          sv_sum[index] /= passdown.sv_count[index];

  /*
   * Check for variables which are zero (e.g. mesh displacements
   * on the first step after a remesh/remap).
   * Set these averages to unity to avoid division by zero.
   */
          if (sv_sum[index] < 1.0e-10) sv_sum[index] = 1.0;
        }

  /* Assign average reciprocals to scale_vec */
      for (iunk=0; iunk<NumUnknowns[pg->imtrx]; iunk++)
        {
          index = passdown.sv_index[iunk];
          scale_vec[iunk] = 1.0 / sv_sum[index];
        }

      safe_free( (void *) sv_sum);
    } /* End of passdown.sv_flag case */


  /* Otherwise, just set scale_vec to all ones */
  else
    {
      for (i=0; i < NumUnknowns[pg->imtrx]; i++) scale_vec[i] = 1.0;
    }

  return;

}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void get_scale_vec_conwrap(int iteration, int flag, double *x, double *rhs)
/* Needed for rSQP */
/* Will fill in later - EDW */
{
return;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
double gsum_double_conwrap(double sum)
/* Put the call to a routine to calculate a global sum.
 * Just return sum for single processor jobs.
 * Input:
 *    sum     Value of double on this processor to be summed on all procs.
 *
 * Output:
 *
 * Return Value:
 *    The global sum is returned on all processors.
 */
{
  struct Aztec_Linear_Solver_System *ams = &(passdown.ams[JAC]);

  if (Num_Proc > 1) return AZ_gsum_double(sum, ams->proc_config);
  else return sum;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int gmax_int_conwrap(int max)
/* Put the call to a routine to calculate a global sum.
 * Just return sum for single processor jobs.
 * Input:
 *    max     Value of integer on this processor to be maxed on all procs.
 *
 * Output:
 *
 * Return Value:
 *    The global max is returned on all processors.
 *
 * Only used by Eigensolver
 */
{
  if (Num_Proc > 1) return AZ_gmax_int(max, passdown.proc_config);
  else return max;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void random_vector_conwrap(double *x, int numOwnedUnks)
/* Put a routine to calculate a random vector.
 * Input:
 *    numOwnedUnks  Length of owned nodes part of x.
 *
 * Output:
 *    x             Random vector.
 *
 * Used by eigensolver only
 */
{
#ifdef HAVE_ARPACK
  struct Aztec_Linear_Solver_System *ams = &(passdown.ams[JAC]);

  AZ_random_vector(x, ams->data_org, passdown.proc_config);
#endif
return;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void perturb_solution_conwrap(double *x, double *x_old,
		              double *scale_vec, int numOwnedUnks)
/* Put a routine to perturb the solution vector a little bit here.
 * This is to move a converged solution at a singularity off
 * of the singularity before doing continuation. This ain't pretty
 * but has helped convergence on some turning point tracking problems.
 * Input:
 *    x_old         Current solution vector.
 *    scale_vec     Work space for a vector to scale x.
 *    numOwnedUnks  Length of owned nodes part of x, x_old, scale_vec
 *
 * Output:
 *    x             Solution vector perturbed a bit.
 *    
 * Return Value:
 */
{
  int i;
  /* int ivar;
  VARIABLE_DESCRIPTION_STRUCT *vdi;
  */

  fill_dvec_rand (x, numOwnedUnks);

  for (i=0; i<numOwnedUnks; i++) {

/* Determine is this is a mesh displacement variable */
 
/* vdi = Index_Solution_Inv(i, NULL, NULL, NULL, NULL, pg->imtrx);
   ivar = vdi->Variable_Type; */
/* If so, don't perturb it! */
/* EDW: This function is being disabled for now. */
/*  if (ivar >= 5 && ivar <= 7) x[i] = 0.0;
    else x[i] = 2.0 * x[i] - 1.0; */
    x[i] = 0.0;
  }

  calc_scale_vec_conwrap(x_old, scale_vec, numOwnedUnks);

  for (i=0; i<numOwnedUnks; i++) x[i] = x_old[i] + 1.0e-5 * x[i] / scale_vec[i];

  exchange_dof(passdown.cx, passdown.dpi, x, pg->imtrx);

}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void solution_output_conwrap(int num_soln_flag,
                             double *x, double param1, double *x2,
                             double param2, double *x3, double param3,
			     int nt, int nits, struct con_struct *con)
/* Put the call to your solution output (both file and screen) routines here.
 * Input:
 *    x            Solution vector.
 *    param        Parameter value (to output in time stamp location)
 *    soln_type_flag  Flag for which solution vector this is. Usually 0 but
 *                    for Phase Transitions the second solution has this flag=1
 *    step_num+1   Time index to output to (step_num is 0 based).
 *    num_its      Number of Newton iterations used for for convergence
 *
 * Output:
 *
 * Return Value:
 */
{
  int error, i_print;
  static int step_print=0;
  static int n_print=0;
  
#ifdef HAVE_ARPACK
  int i, n;
  int freq = eigen->Eigen_Solve_Freq;
#endif

  char *yo = "do_loca";

/* determine whether to print out the data or not */

  if (num_soln_flag == 0) return;

  passdown.last_step = con->stepping_info.last_step;
  i_print = FALSE;
  if (nt == step_print)
    {
      i_print = TRUE;
      step_print += cont->print_freq;
    }
          
  if (i_print)
    {
      error = write_ascii_soln (x,
				passdown.resid_vector,
				NumUnknowns[pg->imtrx],
				passdown.x_AC,
				nAC,
				param1,
				passdown.file);
      if ( error )
        DPRINTF(stdout, "%s:  error writing ASCII soln file\n", yo);
      if (Write_Intermediate_Solutions == 0 && Unlimited_Output)
        {
            write_solution(ExoFileOut,
                           passdown.resid_vector,
                           x,
                           passdown.x_sens_p,
                           passdown.x_old,
                           passdown.xdot,
                           passdown.xdot_old,
                           passdown.tev,
                           passdown.tev_post,
                           NULL,
                           passdown.rd,
                           passdown.gvec,
                           passdown.gvec_elem,
                           &n_print,
                           param1,
                           passdown.theta,
                           param1,
                           NULL,
                           passdown.exo,
                           passdown.dpi);
          n_print++;
        }
    }
      
  /*
   * Backup old solutions
   * can use previous solutions for prediction one day
   * For first-order algorithm, defer this until after resolve.
   */
  if (passdown.method != FIRST_ORDER_CONTINUATION)
    {
      dcopy1(NumUnknowns[pg->imtrx], passdown.x_older, passdown.x_oldest);
      dcopy1(NumUnknowns[pg->imtrx], passdown.x_old, passdown.x_older);
      dcopy1(NumUnknowns[pg->imtrx], passdown.x, passdown.x_old);
      dcopy1(NumUnknowns[pg->imtrx], passdown.x_sens_temp, passdown.x_sens);
    }

  /* If Eigenvalues requested using ARPACK, calculate them here */
  /* This is also the control point for 3D of 2D stability */
#ifdef HAVE_ARPACK

  /* Decide if eigenvalues will be calculated on this step */
  if (con->eigen_info.Num_Eigenvalues < 1) return;
  i_print = FALSE;
  if (freq > 0 && nt%freq == 0) i_print = TRUE;
  if (freq == -1 && passdown.last_step) i_print = TRUE;

  if (i_print)
    {

/* Set LSA_flag and number of ARPACK calls (one unless doing 3D of 2D) */
      passdown.LSA_flag = TRUE;
      n = (passdown.do_3D_of_2D ? LSA_number_wave_numbers : 1);

  /* Loop over LSA wave numbers, or just zero */
      for (i=0; i<n; i++)
        {

  /* Set current wave number if applicable */
	  
          if (n == 1)
            {
	      if (Linear_Stability == LSA_3D_OF_2D) 
		{
		  EH(-1, "With LOCA, you need to have more than one 3D wave number specified");
		}

              LSA_3D_of_2D_wave_number = 0.0;
            }
          else
            {
              LSA_current_wave_number = i;
              LSA_3D_of_2D_wave_number = LSA_wave_numbers[i];
              DPRINTF (stdout, "\n\t3D of 2D LSA -- Wavenumber %d of %d:   %g\n",
                       LSA_current_wave_number + 1, LSA_number_wave_numbers,
                       LSA_3D_of_2D_wave_number);
            }

  /* Call eigensolver */
          calc_eigenvalues_loca(con);
        }

  /* Now reset LSA_flag */
      passdown.LSA_flag = FALSE;
    }
#endif

}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void eigenvector_output_conwrap(int j, int num_soln_flag, double *xr, double evr,
                                double *xi, double evi, int step_num)
/* Call to write out eigenvectors
 * Input:
 *    j    Eigenvalue number
 *    num_soln_flag  =1 for real vector, real eigenvalue
                     =2 for complex (imaginary part has info)
 *    xr   Real part of eigenvector
 *    evr  Real part of eigenvalue
 *    xi   Imaginary part of eigenvector (NULL if num_soln_flag==1)
 *    evi  Imaginary part of eigenvalue
 *    step_num  integer step number for use in output
 *
 * Output:
 *
 * Return Value:
 */
{
  static int old_step = -1;
  static int nprint = -1;
  int i_print;
/*  int n = LSA_current_wave_number; */
  int freq = eigen->Eigen_Write_Freq;
  char efile[MAX_FNL];

  if (j >= eigen->Eigen_Record_Modes) return;

  /* Decide whether to output eigenvectors on this call */
  i_print = FALSE;
  if (freq > 0 && step_num%freq == 0) i_print = TRUE;
  if (freq == -1 && passdown.last_step) i_print = TRUE;
  if (!i_print) return;

  /* Determine if this is a new step number, increment step counter if so */
  if (step_num != old_step) nprint++;

  /* Get the eigenvector file name */
  strcpy(efile, eigen->Eigen_Output_File);
  get_eigen_outfile_name(efile, j, LSA_current_wave_number);

  /* Write the real vector using the real eigenvalue part as the time stamp */
    write_solution(efile,
                   passdown.resid_vector,
                   xr,
                   passdown.x_sens_p,
                   passdown.x_old,
                   passdown.xdot,
                   passdown.xdot_old,
                   passdown.tev,
                   passdown.tev_post,
                   NULL,
                   passdown.rd,
                   passdown.gvec,
                   passdown.gvec_elem,
                   &nprint,
                   0.0,
                   passdown.theta,
                   evr,
                   NULL,
                   passdown.exo,
                   passdown.dpi);

  /* Write the imag vector using the imag eigenvalue part as the time stamp */
  if (num_soln_flag == 2)
    {
      if (j == (eigen->Eigen_Record_Modes-1) )
        {
          DPRINTF(stderr, "Not enough modes requested - cannot write imaginary part!");
        }
      else

  /* Write imaginary part to next eigenvector file */
        {
          strcpy(efile, eigen->Eigen_Output_File);
          get_eigen_outfile_name(efile, j+1, LSA_current_wave_number);
            write_solution(efile,
                           passdown.resid_vector,
                           xi,
                           passdown.x_sens_p,
                           passdown.x_old,
                           passdown.xdot,
                           passdown.xdot_old,
                           passdown.tev,
                           passdown.tev_post,
                           NULL,
                           passdown.rd,
                           passdown.gvec,
                           passdown.gvec_elem,
                           &nprint,
                           0.0,
                           passdown.theta,
                           evi,
                           NULL,
                           passdown.exo,
                           passdown.dpi);
        }
    }

  if (Debug_Flag > 1)
    {
      DPRINTF(stdout, "Mode %d eigenvector recorded\n", j);
    }

  /* Save current step_num for next call */
  old_step = step_num;
  return;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
double free_energy_diff_conwrap(double *x, double *x2)
/* Call to return the free energy difference betwen two solutions
 * Input:
 *    x    One solution vector
 *    x2   Second solution vector
 *
 * Output:
 *
 * Return Value:
 *    The difference in the free energy beween the two solutions
 */
{
#ifdef TRAMONTO
  return calc_free_energy(NULL,x,1.0,1.0)-calc_free_energy(NULL,x2,1.0,1.0);
#else
  return 0.0;
#endif
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
static void print_final(double param, int step_num, int mat_fills,
		        int res_fills, int linear_its)

/*
 * Print out the final results and counters
 */

{
  printf("\n"); /*print_line("~", 80);*/
  printf("CONTINUATION ROUTINE HAS FINISHED: \n");
  printf("\tEnding Parameter value     = %g\n", param);
  printf("\tNumber of steps            = %d\n", step_num+1);
  printf("\tNumber of Matrix fills     = %d\n", mat_fills);
  printf("\tNumber of Residual fills   = %d\n", res_fills);
  if (Linear_Solver == AZTEC)
       {
        printf("\tNumber of linear solve its = %d\n", linear_its);
        fprintf(stdout,"\tNumber of linear solve its = %d\n", linear_its);
        }
  printf("\n"); 
  DPRINTF(stdout,"\n\n\t I will continue no more!\n\t No more continuation for you!\n");

}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
