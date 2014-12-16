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
 *$Id: rf_solve.c,v 5.21 2010-03-17 22:23:54 hkmoffa Exp $
 */

/*
 * Revision history has goneaway from this place.
 */

#ifdef USE_RCSID
static char rcsid[] = "$Id: rf_solve.c,v 5.21 2010-03-17 22:23:54 hkmoffa Exp $";
#endif

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

#define _RF_SOLVE_C
#include "goma.h"
#include "el_quality.h"

#ifdef HAVE_FRONT
extern int mf_setup
PROTO((int *,			/* nelem_glob */
       int *,			/* neqn_glob */
       int *,			/* mxdofel */
       int *,			/* nfullsum */
       int *,			/* symflag */
       int *,			/* nell_order */
       int *,			/* el_proc_assign */
       int *,			/* level */
       int *,			/* nopdof */
       int *,			/* loc_dof */
       int *,			/* constraint */
       const char *,		/* cname */
       int *));			/* allocated */
#endif

/*
 * Global variables defined in this file.
 */

struct elem_side_bc_struct **First_Elem_Side_BC_Array;
struct elem_edge_bc_struct **First_Elem_Edge_BC_Array;

#define ROUND_TO_ONE 0.9999999


int w;

/*
 * Declarations of static functions defined in this file.
 */

static void predict_solution
PROTO((int ,			/* N */
       double ,			/* delta_t */
       double ,			/* delta_t_old */
       double ,			/* delta_t_older */
       double ,			/* theta */
       double  [],		/* x */
       double  [],		/* x_old */
       double  [],		/* x_older */
       double  [],		/* x_oldest */
       double  [],		/* xdot */
       double  [],	        /* xdot_old */
       double  []));		/* xdot_older */

static void predict_solution_newmark
PROTO((int ,			/* N */
       double ,			/* delta_t */
       double  [],		/* x */
       double  [],		/* x_old */
       double  [],		/* xdot */
       double  []));	        /* xdot_old */

static int discard_previous_time_step
PROTO (( int, 
		 double *,
		 double *,
		 double *,
		 double *,
		 double *,
		 double *, 
		 double *));


static void shift_nodal_values
PROTO(( int,
		double,
		double *,
		int ));

void
solve_problem(Exo_DB *exo,	 /* ptr to the finite element mesh database  */
	      Dpi *dpi,		 /* distributed processing information       */
              dbl *te_out)	 /* te_out - return actual end time */

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
  int    *ija = NULL;            /* column pointer array                     */
  double *a = NULL;              /* nonzero array                            */
  double *a_old = NULL;          /* nonzero array                            */
  static double *x = NULL;       /* solution vector                          */

  int     iAC;                   /* Counter                                  */
  static double *x_AC = NULL;           /* Solution vector of extra unknowns          */
  static double *x_AC_old = NULL;       /* old solution vector of extra unknowns      */
  static double *x_AC_older = NULL;     /* older solution vector of extra unknowns    */
  static double *x_AC_oldest = NULL;    /* oldest solution vector of extra unknowns   */
  static double *x_AC_dot = NULL;       /* current time derivative of extra unknowns  */
  static double *x_AC_dot_old = NULL;   /* old time derivative of extra unknowns      */
  static double *x_AC_dot_older = NULL; /* Older time derivative of extra unknowns    */
  static double *x_AC_pred = NULL;      /* predicted extraunknowns */
  
  int    *ija_attic = NULL;      /* storage for external dofs                  */

  int eb_indx, ev_indx;

  /* 
   * Variables
   */
  double *x_pred = NULL;                /* prediction of solution vector     */
  static double *x_old = NULL;          /* old solution vector               */
  static double *x_older = NULL;        /* older solution vector             */
  static double *x_oldest = NULL;       /* oldest solution vector saved      */
  static double *xdot = NULL;           /* current time derivative of soln   */
  static double *xdot_old = NULL;       /* old time derivative of soln       */
  static double *xdot_older = NULL;     /* old time derivative of soln       */

  double *x_sens = NULL;	 /* solution sensitivity                     */
  double **x_sens_p = NULL;	 /* solution sensitivity for parameters      */
  int num_pvector=0;             /* number of solution sensitivity vectors   */

  /* sparse variables for fill equation subcycling */

#ifdef COUPLED_FILL
  static struct Aztec_Linear_Solver_System *ams[NUM_ALSS]={NULL}; 
#else
  int    *ijaf = NULL;           /* column pointer array for fill equation   */
  double *afill = NULL;	         /* nonzero array  for fill equation         */
  double *afill_old = NULL;      /* nonzero array  for fill equation         */
  double *xf = NULL;             /* solution vector for fill equation        */
  double *rf = NULL;             /* residual vector for fill equation        */
  int    *ijaf_attic = NULL;     /* for hiding external nodes from Aztec     */
  static struct Aztec_Linear_Solver_System *ams[NUM_ALSS] = {NULL, NULL};

  /* 
   * Variables for time integration 
   */
     
  double *xf_old = NULL;         /* old solution vector                      */
  double *xfdot = NULL;          /* current time derivative of soln          */
  double *xfdot_old = NULL;      /* old time derivative of soln              */

  double *xf_save = NULL;
  
  double delta_t_exp = 0.0;
  int    n_exp = 0;

#endif /* COUPLED_FILL */
  
  /* "sl_util_structs.h" */

  double *x_update = NULL;       /* update at last iteration                 */

  double *resid_vector = NULL;   /* residual                                 */
  double *resid_vector_sens = NULL; /* residual sensitivity                  */

  double *scale = NULL; 	 /* scale vector for modified newton         */

  int 	 *node_to_fill = NULL;	

  char   tspstring[MAX_FNL];	 /* literal representation of time step 
				  * parameter, theta [0=BE,.5=CN,1=FE] but
				  * any float possible...                    */

  int	 n;                      /* total number of time steps attempted     */
  int	 nt;                     /* total number of successful time steps    */
  int    last_renorm_nt;         /* time step at which last renorm occured   */
  int	 time_step_reform;       /* counter for jacobian reformation stride  */
  int	 converged = TRUE;       /* success or failure of Newton iteration   */
  int	 success_dt = TRUE;      /* success or failure of time step          */
  int    failed_recently_countdown = 0;
  int    i, num_total_nodes;
  int    numProcUnknowns;
  int    const_delta_t, const_delta_ts, step_print;
  int    step_fix = 0;           /* What step to fix the problem on */
  int    good_mesh = TRUE;
  int    w;                      /* counter for looping external variables */
  static int nprint = 0;
  double time_print, i_print;
  double theta = 0.0, time;
  static double time1 = 0.0;     /* Current time that the simulation is trying  to find the solution for */
#ifdef LIBRARY_MODE
  static double delta_t_save = 0.0;
#endif
  double delta_t, delta_t_new = 0.0;
  double delta_t_old, delta_t_older, delta_t_oldest = 0.0;
  double timeValueRead = 0.0;    /* time value read from an exodus input file 
				  * used to initialize the solution vector */
  double timeValueReadTrans = 0.0;
  //  static double time_value = 0.0;
  double eps;

  double time2 = 0.0;
#ifdef RESET_TRANSIENT_RELAXATION_PLEASE
  double damp_factor_org[2]={damp_factor1,damp_factor2};
#endif
  /*
   * Other local variables...
   */
  
  int         error, err, is_steady_state, inewton;
  int 		*gindex = NULL, gsize;
  int		*p_gsize;
  static double	*gvec=NULL;
  static double	***gvec_elem=NULL;
  FILE *file=NULL; 
  static struct Results_Description  *rd;
  struct Level_Set_Data *ls_old;
  
  int		tnv;		 /* total number of nodal variables and kinds */
  int		tev;		 /* total number of elem variables and kinds  */
  int		tnv_post;	 /* total number of nodal variables and kinds  
				    for post processing                       */
  int		tev_post;	 /* total number of elem variables and kinds 
				    for post processing                       */

  double *gv;                    /* Global variable values */
  double *x_pp=NULL;             /* Post-proc variables for export */
  static int *xp_id=NULL;	 /* Post-proc variable ID */
  double *base_p_por=NULL;	 /* Base values for porosity updates */
  double *base_p_liq=NULL;	 /* Base values for porosity updates */
  int update_porosity=FALSE;	 /* Flag for external porosity updates */
#ifdef HAVE_FRONT  
  int max_unk_elem, one, three;  /* variables used as mf_setup arguments      */
#endif
  unsigned int  matrix_systems_mask;
  int did_renorm;                /* Flag indicating if we renormalized.       */
  int Renorm_Now = FALSE;        /* Flag forcing renormalization regardless of gradient */
  double evol_local=0.0;
  double lsvel_local=0.0;
#ifdef PARALLEL
  double evol_global=0.0;
  double lsvel_global=0.0;
#endif /* PARALLEL */

  static int callnum = 1;	 /* solve_problem call counter */
  int last_call = TRUE;		 /* Indicates final rf_solve call */
#ifdef LIBRARY_MODE
  int last_step = FALSE;	 /* Indicates final time step on this call */
#endif
#ifdef RELAX_ON_TRANSIENT_PLEASE
  int relax_bit = TRUE;	 /* Enables relaxation after a transient convergence failure*/
#else
  int relax_bit = FALSE;	
#endif

  static const char yo[]="solve_problem"; /* So my name is in a string.        */

  /*
   * 		BEGIN EXECUTION
   */

#ifdef DEBUG
  fprintf(stderr, "P_%d solve_problem() begins...\n",ProcID);
#endif /* DEBUG */

  /* Set step_fix only if parallel run and only if fix freq is enabled*/
  if (Num_Proc > 1 && tran->fix_freq > 0) {
    step_fix = 1; /* Always fix on the first timestep to match print frequency */
  }

  tran->time_value = time1;

  is_steady_state = ( TimeIntegration == STEADY ) ? TRUE : FALSE;

  p_gsize = &gsize;

#ifdef LIBRARY_MODE
  fprintf(stderr, "  Commencing call #%3d from ANIMAS to solve_problem\n",
          callnum);
  if (libio->animas_step != -1) last_call = FALSE;

  /* Determine if external porosity updates are required */
  if (Num_Var_In_Type[POR_LIQ_PRES] && efv->ev_porous_index > -1)
  {
    update_porosity = TRUE;
    fprintf(stderr, " External porosity field %d will be updated.\n",
              efv->ev_porous_index);
  }
  if (libio->goma_first == 1) fprintf(stderr, "  Goma goes first");
  if (libio->goma_first == 0) fprintf(stderr, "  Goma goes second");

#endif   
 
  if (Unlimited_Output && strlen(Soln_OutFile))  {
    file = fopen(Soln_OutFile, "w");
    if (file == NULL) {
      fprintf(stderr, "%s:  opening soln file, %s, for writing\n", 
	      yo, Soln_OutFile);
      EH(-1, "Can not open solution file\n");
    }
  }
  
  /* set problem flags for writing exodus data base and solving problem */

  /*
   * Some preliminaries to help setup EXODUS II database output.
   */

#ifdef DEBUG
  fprintf(stderr, "P_%d cnt_nodal_vars() begins...\n",ProcID);
#endif /* DEBUG */

  tnv = cnt_nodal_vars();
  /*  tnv_post is calculated in load_nodal_tkn*/
  tev = cnt_elem_vars();
  /*  tev_post is calculated in load_elem_tkn*/
  
#ifdef DEBUG
  fprintf(stderr, "Found %d total primitive nodal variables to output.\n", tnv);
  fprintf(stderr, "Found %d total primitive elem variables to output.\n", tev);
#endif /* DEBUG */
  
  if (tnv < 0)
  {
    DPRINTF(stderr, "%s:\tbad tnv.\n", yo);
    EH(-1, "\t");
  }
  
  if ( tev < 0 )
  {
    DPRINTF(stderr, "%s:\tMaybe bad tev? See goma design committee ;) \n", yo);
    EH(-1, "\t");
  }

  /*
   * When using overlap AC's for fluid/solid overlapping grid problems,
   * create the additional set of AC constraints here.
   */
  if (Do_Overlap && augc[nAC-1].Type == AC_OVERLAP)
    {
      err = create_overlap_acs(exo, nAC-1);
      EH(err, "Problem with create_overlap_acs!");
    }
  if (nAC > 0 && augc[0].Type == AC_PERIODIC)
    {
      err = create_periodic_acs(exo);
      EH(err, "Problem with create_periodic_acs!");
    }

  /*
   *  Malloc the space for the results description structure and set all
   *  of that space to naught initially.
   *  Do this only once if in library mode.
   */
  if (callnum == 1)
    {
      rd = alloc_struct_1(struct Results_Description, 1);
#ifdef LIBRARY_MODE
      xp_id = alloc_int_1(MAX_EXTERNAL_FIELD, 0);
#endif
    }
  
  rd->nev = 0;			/* number element variables in results */
  rd->ngv = 0;			/* number global variables in results  */
  rd->nhv = 0;			/* number history variables in results */


    rd->ngv = 5 + nAC;			/* number global variables in results 
					   see load_global_var_info for names*/
    error = load_global_var_info(rd, 0, "CONV");
    error = load_global_var_info(rd, 1, "NEWT_IT");
    error = load_global_var_info(rd, 2, "MAX_IT");
    error = load_global_var_info(rd, 3, "CONVRATE");
    error = load_global_var_info(rd, 4, "MESH_VOLUME");

    if ( rd->ngv > MAX_NGV ) 
      EH(-1, "Augmenting condition values overflowing MAX_NGV.  Change and rerun .");

  if ( nAC > 0   )
  {
    char name[7];

    for( i = 0 ; i < nAC ; i++ )
      {
	sprintf(name, "AUGC_%d",i+1);
	error = load_global_var_info(rd, 5 + i, name);
      }
  }

  gv = alloc_dbl_1( rd->ngv, 0.0 );
  
  /*
   *  Load output nodal types, kinds and names into the structure
   *  which will be used to define what's in the output file.
   */
  error = load_nodal_tkn(rd, &tnv, &tnv_post);  
  if (error !=0) {
    DPRINTF(stderr, "%s:  problem with load_nodal_tkn()\n", yo);
    EH(-1,"\t");
  }

  /*
   * Post-processing vars are stored in the static local array xp_id,
   * and must be used to reset Export_XP_ID on calls after the first.
   * This is a temporary fix!
   */ 
#ifdef LIBRARY_MODE
  if (callnum == 1)
    {
      for (i=0; i<MAX_EXTERNAL_FIELD; i++) xp_id[i] = Export_XP_ID[i];
    }
  else
    {
      for (i=0; i<MAX_EXTERNAL_FIELD; i++) Export_XP_ID[i] = xp_id[i];
    }
#endif

  /* For retrieving post-processing variables */
  if (Num_Export_XP > 0 && POROUS_SATURATION != -1)
    {
      asdv(&x_pp, Num_Export_XS * exo->num_nodes);
    }

  /*
   *  Load output element var types, kinds and names into the structure
   *  which will be used to define what's in the output file.
   */
  error = load_elem_tkn(rd, exo, tev, &tev_post);  
  if (error !=0) {
    DPRINTF(stderr, "%s:  problem with load_elem_tkn()\n", yo);
    EH(-1,"\t");
  }
#ifdef PARALLEL
  check_parallel_error("Results file error");
#endif /* PARALLEL */

  /* 
   * Write out the names of the nodal variables that we will be sending to
   * the EXODUS II output file later - do only once if in library mode.
   */

#ifdef DEBUG
  fprintf(stderr, "P_%d %d wr_result_prelim() starts...\n",ProcID, tnv);
#endif /* DEBUG */

  if (callnum == 1)
    {
      gvec_elem = (double ***) alloc_ptr_1(exo->num_elem_blocks);
      if ((tev + tev_post) > 0)
        {
          for (i = 0; i < exo->num_elem_blocks; i++)
            {
              gvec_elem[i] = (double **) alloc_ptr_1(tev + tev_post);
            }
        }
      wr_result_prelim_exo(rd, exo, ExoFileOut, gvec_elem);
    }


#ifdef DEBUG
  fprintf(stderr, "P_%d: %d wr_result_prelim_exo() ends...\n", ProcID, tnv);
#endif /* DEBUG */

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

  numProcUnknowns = NumUnknowns + NumExtUnknowns;

#ifdef DEBUG
  fprintf(stderr, "P_%d: numProcUnknowns = %d (%d+%d)\n",ProcID, numProcUnknowns, 
	  NumUnknowns, NumExtUnknowns);
#endif /* DEBUG */
  
  asdv(&resid_vector, numProcUnknowns);
  asdv(&resid_vector_sens, numProcUnknowns);
  asdv(&scale, numProcUnknowns);

#ifdef DEBUG
  fprintf(stderr, "P_%d: begin Solver allocation\n",ProcID);
  fprintf(stderr, "P_%d: NUM_ALSS=%d\n",ProcID, NUM_ALSS);
#endif /* DEBUG */

  /*
   * Allocate Aztec structures and initialize all elements to zero
   */
  if (callnum == 1)
    {
      for (i = 0; i < NUM_ALSS; i++)
        {
           ams[i] = alloc_struct_1(struct Aztec_Linear_Solver_System, 1);
        }				
    }				
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

#ifdef DEBUG
  fprintf(stderr, "P_%d: Solver allocation complete\n",ProcID);
#endif /* DEBUG */

  /* Allocate solution arrays on first call only */
  if (callnum == 1)
    {
      x          = alloc_dbl_1(numProcUnknowns, 0.0);
      x_old      = alloc_dbl_1(numProcUnknowns, 0.0);
      x_older    = alloc_dbl_1(numProcUnknowns, 0.0);
      x_oldest   = alloc_dbl_1(numProcUnknowns, 0.0);
      xdot       = alloc_dbl_1(numProcUnknowns, 0.0);
      xdot_old   = alloc_dbl_1(numProcUnknowns, 0.0);
      xdot_older = alloc_dbl_1(numProcUnknowns, 0.0);
    }
  x_update = alloc_dbl_1(numProcUnknowns + numProcUnknowns, 0.0);

  /* Initialize solid inertia flag */
  set_solid_inertia();

  if(tran->solid_inertia)
    {
      tran->xdbl_dot = (double *) 
	array_alloc(1, numProcUnknowns, sizeof(double));
      tran->xdbl_dot_old = (double *) 
	array_alloc(1, numProcUnknowns, sizeof(double));

      for (i = 0; i < numProcUnknowns; i++)
	{
	  tran->xdbl_dot[i] = 0.;
	  tran->xdbl_dot_old[i] = 0.;
	}

      /*set these in input deck when ready, or if you feel compelled
	to change them.  PRS 10/3/2001 */
      tran->newmark_gamma = 0.9;
      tran->newmark_beta = 0.49;
    }
  
  node_to_fill = alloc_int_1(num_total_nodes, 0);

  /* Allocate sparse matrix */

  if (strcmp(Matrix_Format, "epetra") == 0) {
    err = check_compatible_solver();
    EH(err, "Incompatible matrix solver for epetra, epetra supports amesos and aztecoo solvers.");
    check_parallel_error("Matrix format / Solver incompatibility");
    ams[JAC]->RowMatrix = EpetraCreateRowMatrix(num_internal_dofs + num_boundary_dofs);
    EpetraCreateGomaProblemGraph(ams[JAC], exo, dpi);
  } else if (strcmp(Matrix_Format, "msr") == 0) {
    log_msg("alloc_MSR_sparse_arrays...");
    alloc_MSR_sparse_arrays(&ija, &a, &a_old, 0, node_to_fill, exo, dpi);
    /*
     * An attic to store external dofs column names is needed when
     * running in parallel.
     */
    alloc_extern_ija_buffer(num_universe_dofs,
                            num_internal_dofs + num_boundary_dofs,
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

    ams[JAC]->npn = dpi->num_internal_nodes + dpi->num_boundary_nodes;
    ams[JAC]->npn_plus = dpi->num_internal_nodes + dpi->num_boundary_nodes
        + dpi->num_external_nodes;

    ams[JAC]->npu = num_internal_dofs + num_boundary_dofs;
    ams[JAC]->npu_plus = num_universe_dofs;

    ams[JAC]->nnz = ija[num_internal_dofs + num_boundary_dofs] - 1;
    ams[JAC]->nnz_plus = ija[num_universe_dofs];

    ams[JAC]->RowMatrix = NULL;

  } else if (strcmp(Matrix_Format, "vbr") == 0) {
    log_msg("alloc_VBR_sparse_arrays...");
    alloc_VBR_sparse_arrays(ams[JAC], exo, dpi);
    ija_attic = NULL;
    ams[JAC]->belfry = ija_attic;

    a = ams[JAC]->val;
    if (!save_old_A)
      a_old = ams[JAC]->val_old = NULL;
  } else if (strcmp(Matrix_Format, "front") == 0) {
    /* Don't allocate any sparse matrix space when using front */
    ams[JAC]->bindx   = NULL;
    ams[JAC]->val     = NULL;
    ams[JAC]->belfry  = NULL;
    ams[JAC]->val_old = NULL;
    ams[JAC]->indx  = NULL;
    ams[JAC]->bpntr = NULL;
    ams[JAC]->rpntr = NULL;
    ams[JAC]->cpntr = NULL;
  } else {
    EH(-1, "Attempted to allocate unknown sparse matrix format");
  }
	  
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

  for (i = 0;i < nn_post_fluxes_sens; i++) {
    num_pvector=MAX(num_pvector,pp_fluxes_sens[i]->vector_id);
  }
  for (i = 0; i < nn_post_data_sens; i++)        {
    num_pvector = MAX(num_pvector, pp_data_sens[i]->vector_id);
  }
  
  if ((nn_post_fluxes_sens + nn_post_data_sens) > 0) {
    num_pvector++;
    num_pvector = MAX(num_pvector,2);
    x_sens_p = Dmatrix_birth(num_pvector,numProcUnknowns);
    x_sens = alloc_dbl_1(numProcUnknowns, 0.0);
  }
  else {
    x_sens_p = NULL;
  }

  /* Allocate AC unknown arrays on the first call */
  if (nAC > 0 && callnum == 1)
    {
      x_AC           = alloc_dbl_1(nAC, 0.0);
      x_AC_old       = alloc_dbl_1(nAC, 0.0);
      x_AC_older     = alloc_dbl_1(nAC, 0.0);
      x_AC_oldest    = alloc_dbl_1(nAC, 0.0);
      x_AC_dot       = alloc_dbl_1(nAC, 0.0);
      x_AC_dot_old   = alloc_dbl_1(nAC, 0.0);
      x_AC_dot_older = alloc_dbl_1(nAC, 0.0);
      x_AC_pred      = alloc_dbl_1(nAC, 0.0);
    }

  /* Set initial guess from an input exodus file or other method on the first call only */
  if (callnum == 1) 
    {
      init_vec(x, cx, exo, dpi, x_AC, nAC, &timeValueRead);

      /*
       *  Determine if we should use this time as the initial time in the simulation
       */
      if (TimeIntegration != STEADY)
	{
	  if (tran->init_time < 0.0)  
	    {
	      tran->init_time = timeValueRead;
	      DPRINTF(stdout, "\n Initial Simulation Time Has been set to %g\n", timeValueRead);
	    }
	}
    }

  /* Load external fields from import vectors xnv_in & xev_in */
#ifdef LIBRARY_MODE
  /* First check if porosity updates are necessary */
  if (update_porosity)
    {
      base_p_por = alloc_dbl_1(num_total_nodes, 0.0);
      base_p_liq = alloc_dbl_1(num_total_nodes, 0.0);

  /* Load starting porous liquid pressures from solution vector */
      for (i=0; i<num_total_nodes; i++)
        {
          j = Index_Solution(i, POR_LIQ_PRES, 0, 0, -1);
          if (j > -1) base_p_liq[i] = x[j];
        }
    }

  /* Now load imports, and also base_p_por if applicable */
  error = load_import_fields(base_p_por, exo, callnum);
  EH(error, "Problem with load_import_fields!");
#else
  /****************************Anneal from external***********************/
  if (efv->ev_porous_decouple)
    {
      anneal_mesh_with_external_field(exo);
    }

#endif

  dcopy1( nAC, x_AC, &(gv[5]) );


  /***************************************************************************
   *            STEADY STATE SOLUTION PROCEDURE
   ***************************************************************************/
  if (TimeIntegration == STEADY) {
#ifdef DEBUG
    fprintf(stderr, "P_%d %d beginning STEADY analysis...\n", ProcID, tnv);
#endif /* DEBUG */
      
    theta = 0.0;        /* for steady problems. theta def in rf_fem.h */
    delta_t = 0.0;
      
    find_and_set_Dirichlet(x, xdot, exo, dpi);
      
    matrix_systems_mask = 1;
      
    log_msg("sl_init()...");
    sl_init(matrix_systems_mask, ams, exo, dpi, cx);
    if( nAC > 0  || 
	nn_post_fluxes_sens > 0 ||
	nn_post_data_sens > 0 ) ams[JAC]->options[AZ_keep_info] = 1;
      
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

#ifdef DEBUG
    fprintf(stderr, "Proc_%d %s: starting solve_nonlinear_problem\n",ProcID, yo);
#endif /* DEBUG */

    /* Call prefront (or mf_setup) if necessary */
    if (Linear_Solver == FRONT) {
#ifdef HAVE_FRONT  

      /* Also got to define these because it wants pointers to these numbers */
      max_unk_elem = (MAX_PROB_VAR + MAX_CONC) * MDE;
      one = 1;
      three = 3;

      /* NOTE: We need a overall flag in the vn_glob struct that tells
       * whether FULL_D is on anywhere in domain. 
       * This assumes only one material.  See sl_front_setup for test.
       * that test needs to be in the input parser. 
       */
      if (vn_glob[0]->dg_J_model == FULL_DG) {
	max_unk_elem = (MAX_PROB_VAR + MAX_CONC)*MDE + 4*vn_glob[0]->modes*4*MDE;
	three = 0; 
      }
	    
#ifdef PARALLEL
      if (Num_Proc > 1) EH(-1, "Whoa.  No front allowed with nproc>1");
      check_parallel_error("Front solver not allowed with nprocs>1");
#endif /* PARALLEL */
	  
      err = mf_setup(&exo->num_elems, &NumUnknowns, 
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

#else /* HAVE_FRONT */
      EH(-1,"Don't have frontal solver compiled and linked in");
#endif /* HAVE_FRONT */
    }
      
    if (nEQM > 0)
      {
        DPRINTF(stderr, "\nINITIAL ELEMENT QUALITY CHECK---\n");
        good_mesh = element_quality(exo, x, ams[0]->proc_config);
      }

    err = solve_nonlinear_problem(ams[JAC], x, delta_t, theta,
				  x_old, x_older, xdot, xdot_old,
				  resid_vector, x_update, scale,  
				  &converged, &nprint, tev, tev_post, gv,
				  rd, gindex, p_gsize, gvec, gvec_elem, 
				  time1, exo, dpi, cx, 0, 
				  &time_step_reform, is_steady_state,
 				  x_AC, x_AC_dot, time1, resid_vector_sens,
				  x_sens, x_sens_p, NULL);

#ifdef DEBUG
    fprintf(stderr, "%s: returned from solve_nonlinear_problem\n", yo);
#endif /* DEBUG */
      
    if (!converged) {
      DPRINTF(stderr, 
	      "\n\tFailure: could not converge on the steady state solution.\n");
      DPRINTF(stderr, "\tSorry, I really tried...\n");
      if (Linear_Stability) {
	DPRINTF(stderr, 
		"\tThere's no point in solving the eigensystem (but I'll do it anyway).\n");
      }
    }

    log_msg("Returning from solve_nonlinear_problem with %d", err);
    EH(err, "Problem from solve_nonlinear_problem.");
  
      /* Check element quality */
      good_mesh = element_quality(exo, x, ams[0]->proc_config);

    if (file != NULL) {
      error = write_ascii_soln(x, resid_vector, numProcUnknowns,
			       x_AC, nAC, 0.0, file);
    }
      
    EH(error, "Error writing ASCII soln file.");
      
    if (Write_Intermediate_Solutions == 0) {
      nprint = 0;

      write_solution(ExoFileOut, resid_vector, x, x_sens_p, 
		     x_old, xdot, xdot_old, tev, tev_post, gv, 
		     rd, gindex, p_gsize, gvec, gvec_elem, 
		     &nprint, delta_t, theta, time1, x_pp, 
		     exo, dpi);
    } /* end of if Write Intermediate Solutions */

    /* Print out values of extra unknowns from augmenting conditions */
    if (nAC > 0)
    {
      DPRINTF(stderr, "\n------------------------------\n");
      DPRINTF(stderr, "Augmenting Conditions:    %4d\n", nAC);
      DPRINTF(stderr, "Number of extra unknowns: %4d\n\n", nAC);

      for(iAC = 0; iAC < nAC; iAC++)
      {
	if(augc[iAC].Type == AC_USERBC)
	{
	  DPRINTF(stderr, "\tBC[%4d] DF[%4d]=% 10.6e\n", 
		  augc[iAC].BCID, augc[iAC].DFID, x_AC[iAC]);
	}
      else if(augc[iAC].Type == AC_USERMAT ||
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
#endif /* PARALLEL */
	  DPRINTF(stderr, "\tMT[%4d] VC[%4d]=%10.6e Param=%10.6e\n", 
		  augc[iAC].MTID, augc[iAC].VOLID, evol_local, 
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
#endif /* PARALLEL */

	  lsvel_local = lsvel_local / evol_local;

	  DPRINTF(stderr, "\tMT[%4d] LSVEL phase[%4d]=%10.6e Param=%10.6e\n", 
		  augc[iAC].MTID, augc[iAC].LSPHASE, lsvel_local, 
		  x_AC[iAC]);
	}
	else if(augc[iAC].Type == AC_FLUX)
	{
	  DPRINTF(stderr, "\tBC[%4d] DF[%4d]=%10.6e\n", 
		  augc[iAC].BCID, augc[iAC].DFID, x_AC[iAC]);
	}
      }
    }

    /* Integrate fluxes, forces  (only if converged)
     */
    if( converged )	{
      for (i = 0; i < nn_post_fluxes; i++) {
	(void) evaluate_flux(exo, dpi, 
			     pp_fluxes[i]->ss_id, 
			     pp_fluxes[i]->flux_type,
			     pp_fluxes[i]->flux_type_name,
			     pp_fluxes[i]->blk_id,
			     pp_fluxes[i]->species_number, 
			     pp_fluxes[i]->flux_filenm,
			     pp_fluxes[i]->profile_flag,
			     x, xdot, NULL, delta_t, time1, 1); 
      }
      
      /* Compute flux, force sensitivities
       */
      for (i = 0; i < nn_post_fluxes_sens; i++) {
	(void) evaluate_flux_sens(exo, dpi, 
				  pp_fluxes_sens[i]->ss_id,
				  pp_fluxes_sens[i]->flux_type,
				  pp_fluxes_sens[i]->flux_type_name,
				  pp_fluxes_sens[i]->blk_id,
				  pp_fluxes_sens[i]->species_number,
				  pp_fluxes_sens[i]->sens_type,
				  pp_fluxes_sens[i]->sens_id,
				  pp_fluxes_sens[i]->sens_flt,
				  pp_fluxes_sens[i]->sens_flt2,
				  pp_fluxes_sens[i]->vector_id,
				  pp_fluxes_sens[i]->flux_filenm,
				  pp_fluxes_sens[i]->profile_flag,
				  x, xdot, x_sens_p, delta_t,
				  time1, 1);
      }
      
      /*
       * Compute global volumetric quantities
       */
#if 1
	/* This section is a kludge to determine minimum and maximum
 		variable values for outputing surfaces of constant concentration
 	*/
 	if( nn_volume )
 	  {
 	  int search_minmax = FALSE;
  	  for (i = 0; i < nn_volume; i++) {
  		if(pp_volume[i]->volume_type == I_SURF_SPECIES &&
				pp_volume[i]->species_no == -1 ) 
 		search_minmax=TRUE;
 		}
 	  if( search_minmax && Num_Proc > 1)
 		  WH(-1,"Species surfaces not recommended in parallel\n");
 	  if( search_minmax && Num_Proc == 1)
  	  {
  		  int inode,offset,idof;
  		  VARIABLE_DESCRIPTION_STRUCT *vd;
  		  double specmax[MAX_CONC],specmin[MAX_CONC],specavg, fraction;
 		  int eq_off = 0;
 
 		  if(pd->e[R_ENERGY] )  eq_off++;
 		  if(pd->e[R_MOMENTUM1] )  eq_off += pd->Num_Dim;
  		  inode = 0;
  		  for(i=0 ; i < (num_internal_dofs+num_boundary_dofs); i++) {
  			  vd = Index_Solution_Inv(i, &inode, NULL, &offset, &idof);
  			  if( vd->Variable_Type == MASS_FRACTION )      {
  				if( i >= pd->Num_Species_Eqn + eq_off ) {
  		                  if (x[i] > specmax[offset-eq_off]) 
 						specmax[offset-eq_off] = x[i];
  				  if (x[i] < specmin[offset-eq_off]) 
 						specmin[offset-eq_off] = x[i];
  		                  } else  {
  	                          specmin[offset-eq_off] = x[i];
  	                          specmax[offset-eq_off] = x[i];
  		                  }
  		          }
  		  }
  		  DPRINTF(stderr,"Species concentration minmax: %d \n",ProcID);
  	  	  for (i = 0; i < pd->Num_Species_Eqn; i++) {
  			DPRINTF(stderr,"%d %g %g \n",i,specmin[i],specmax[i]);
  			}
  	  for (i = 0; i < nn_volume; i++) {
  		if(pp_volume[i]->volume_type == I_SURF_SPECIES) 
  			{
          DPRINTF(stderr,"params %d %d %d\n",i,pp_volume[i]->num_params,pd->Num_Species);
            		if(pp_volume[i]->num_params < pd->Num_Species+1)
                 		{
         DPRINTF(stderr,"params %d %d\n",pp_volume[i]->num_params,pd->Num_Species);
         WH(-1,"SURF_SPECIES parameters should number Num_Species+1");
                		}
 			if( pp_volume[i]->species_no == -1 )
 				{
				fraction = pp_volume[i]->params[pd->Num_Species];
				pp_volume[i]->params[pd->Num_Species] = 
					pp_volume[i]->params[pd->Num_Species-1];
				for(idof=0 ; idof<pd->Num_Species_Eqn; idof++)
					{
					specavg = specmin[idof] + 
					  fraction*(specmax[idof]-specmin[idof]);
					pp_volume[i]->params[pd->Num_Species] +=
					  pp_volume[i]->params[idof]*specavg;
					}
DPRINTF(stderr,"new surface value = %g \n",pp_volume[i]->params[pd->Num_Species]);
				}
 			}
 		}	/* nn_volume loop	*/
 	  }	/*search  */
	}	/*  nn_volume	*/
#endif
	for (i = 0; i < nn_volume; i++ ) {
	  evaluate_volume_integral(exo, dpi,
				   pp_volume[i]->volume_type,
				   pp_volume[i]->volume_name,
				   pp_volume[i]->blk_id,
				   pp_volume[i]->species_no,
				   pp_volume[i]->volume_fname,
				   pp_volume[i]->params,
				   NULL,  x, xdot, delta_t,
				   time1, 1);
	}
    }	/* if converged */
    

    if (Anneal_Mesh) {
      /*
       * Transform the node point coordinates according to the
       * displacements and write out all the results using the
       * displaced coordinates. Set the displacement field to
       * zero, too.
       */

#ifdef DEBUG
      fprintf(stderr, "%s: anneal_mesh()...\n", yo);
#endif /* DEBUG */
      err = anneal_mesh(x, tev, tev_post, gv,  rd, time1, exo, dpi);
#ifdef DEBUG
      fprintf(stderr, "%s: anneal_mesh()-done\n", yo);
#endif /* DEBUG */
      EH(err, "anneal_mesh() bad return.");
    }

    if (Linear_Stability) {
      err = solve_stability_problem(ams[JAC], x, delta_t, theta,
				    resid_vector, x_old, x_older,
				    xdot, xdot_old, x_update,
				    &converged, &nprint, tnv,
				    tnv_post, rd, gindex, p_gsize,
				    gvec, time1, exo, dpi);
      EH(err, "Problem from solve_stability_problem.");
    }
    
    if(Particle_Dynamics)
      {
	/* If this was steady-state we need to ensure the fv_old* values are the same as fv. */
	dcopy1(NumUnknowns, x, x_old);
	dcopy1(NumUnknowns, x, x_older);
	initialize_particles(exo, x, x_old, xdot, xdot_old, resid_vector);
	for(n = 0; n < Particle_Max_Time_Steps; n++)
	  {
	    time = (n+1) * Particle_Output_Time_Step;
	    DPRINTF(stderr, "\nComputing particles for time %g (%2.0f%% done)\n", time, (dbl)n/(dbl)Particle_Max_Time_Steps*100.0);
	    err = compute_particles(exo, x, x_old, xdot, xdot_old, resid_vector, time, Particle_Output_Time_Step, n);
	    EH(err, "Error performing particle calculations.");
	  }
      }
  }  /* if(steady) */

  /********************************************************************************
   *                            Transient solution process 
   ********************************************************************************/
  else
    {
      if (Debug_Flag && ProcID == 0) {
	fprintf(stderr,"MaxTimeSteps: %d \tTimeMax: %f\n",MaxTimeSteps,TimeMax);
	fprintf(stderr,"solving transient problem\n");
      }
    
    /*
     *  Transfer information from the Transient_Information structure to local variables
     */
    Delta_t0     = tran->Delta_t0;
    Delta_t_min  = tran->Delta_t_min;
    Delta_t_max  = tran->Delta_t_max;
    MaxTimeSteps = tran->MaxTimeSteps;
    TimeMax      = tran->TimeMax;
    eps          = tran->eps;
#ifndef COUPLED_FILL
    exp_subcycle = tran->exp_subcycle;
#endif /* not COUPLED_FILL */   

    // Determine if we are using a constant time step or not
    if (Delta_t0 < 0.0 ) {
      Delta_t0	    = -Delta_t0;
      const_delta_t = 1;
    }  else {
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
    if (libio->solve_steps > 0)
      {
        MaxTimeSteps = libio->solve_steps;
        if (callnum > 1) Delta_t0 = delta_t_save;
      }
    else if (libio->solve_steps < 0)
      {
        MaxTimeSteps = -libio->solve_steps;
        if (callnum > 1) Delta_t0 = delta_t_save;
      }
    else
      {
	if (callnum > 1 && delta_t_save*libio->decelerator > Delta_t0 ) 
	  {
	   Delta_t0 = delta_t_save*libio->decelerator;
	  }
      }
#endif

    time   = time1 = tran->init_time; /* Allow non-zero initial time */
    tran->time_value = tran->time_value_old = time1;
    tran->delta_t_old = Delta_t0;
    if (Delta_t0 > Delta_t_max) Delta_t0 = Delta_t_max;
    delta_t = delta_t_old = delta_t_older = Delta_t0;
    tran->delta_t = delta_t;    /*Load this up for use in load_fv_mesh_derivs */
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
    exchange_dof(cx, dpi, x);

    /*
     * Now copy the initial solution, x[], into the history solutions
     * x_old[], etc. Note, xdot[] = xdot_old[] = 0 at this point,
     * which is in agreement with the specification of the history
     * solutions. 
     */
    dcopy1(numProcUnknowns, x,       x_old);
    dcopy1(numProcUnknowns, x_old,   x_older);
    dcopy1(numProcUnknowns, x_older, x_oldest);

    if (nAC > 0) {
      dcopy1(nAC, x_AC,       x_AC_old);
      dcopy1(nAC, x_AC_old,   x_AC_older);
      dcopy1(nAC, x_AC_older, x_AC_oldest);
    }

    /* initialize the counters for when to print out data */
    time_print		= time;
    step_print		= 1;
    matrix_systems_mask	= 1;

    /*set some other misc. action flags */
    af->Sat_hyst_reevaluate = FALSE;

#ifndef COUPLED_FILL
    /*
     * Initialize the level set solution vectors and
     * associated matrix problem 
     */
    if (Explicit_Fill) {
      countmap_vardofs(FILL, num_total_nodes, node_to_fill);

      asdv(&xf,        num_fill_unknowns);
      asdv(&xf_old,    num_fill_unknowns);

      asdv(&xfdot,     num_fill_unknowns);
      asdv(&xfdot_old, num_fill_unknowns);

      alloc_MSR_sparse_arrays(&ijaf, &afill, &afill_old, 1, node_to_fill, 
			      exo, dpi);

      alloc_extern_ija_buffer(num_fill_unknowns, owned_fill_unknowns,
			      ijaf, &ijaf_attic);
      asdv(&rf,        num_fill_unknowns);
      asdv(&xf_save,   num_fill_unknowns);

      ams[FIL]->bindx = ijaf;
      ams[FIL]->val   = afill;

      ams[FIL]->belfry = ijaf_attic;

      ams[FIL]->val_old = NULL;

      ams[FIL]->indx  = NULL;
      ams[FIL]->bpntr = NULL;
      ams[FIL]->rpntr = NULL;
      ams[FIL]->cpntr = NULL;

      matrix_systems_mask += 2;
    }
#endif /* not COUPLED_FILL */
      
    /* Call prefront (or mf_setup) if necessary */
    if (Linear_Solver == FRONT) {
#ifdef HAVE_FRONT 	  

      /* Also have to define these because it wants pointers to these numbers */
          
      max_unk_elem = (MAX_PROB_VAR + MAX_CONC)*MDE;
      one	   = 1;
      three	   = 3;
      /* NOTE: We need a overall flag in the vn_glob struct that tells whether FULL_DG
	 is on anywhere in domain.  This assumes only one material.  See sl_front_setup for test.
	 That test needs to be in the input parser.  */
      if(vn_glob[0]->dg_J_model == FULL_DG) 
	  max_unk_elem = (MAX_PROB_VAR + MAX_CONC)*MDE + 4*vn_glob[0]->modes*4*MDE;
	
      err = mf_setup(&exo->num_elems, 
		     &NumUnknowns, 
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
#else /* HAVE_FRONT */
      EH(-1,"Don't have frontal solver compiled and linked in");
#endif /* HAVE_FRONT */
    }
      
#ifdef COUPLED_FILL
    /*
     * Prior to the primary time stepping loop, do any one time 
     * initialization required for the Aztec linear solver, for
     * the full Jacobian linear system.
     */
#else /* COUPLED_FILL */
    /*
     * Prior to the primary time stepping loop, do any one time 
     * initialization required for the Aztec linear solver, both
     * for the full Jacobian linear system *and* the explicit fill
     * equation linear system.
     */
#endif /* COUPLED_FILL */

    /*
     * Now, just pass pointer to ams structure with all Aztec stuff
     * bundled inside. At the other end, extract "ija" and "a" as
     * appropriate, but the other items are there now, too.
     * Do this only once if in library mode.
     */
    if (callnum == 1) sl_init(matrix_systems_mask, ams, exo, dpi, cx);	
      
    /*
     * make sure the Aztec was properly initialized
     */
    check_parallel_error("Aztec Initialization");
      
    /* Set the number of successful time steps, nt, to zero */
    nt		     = 0;   
    time_step_reform = Time_Jacobian_Reformation_stride;
    const_delta_ts   = const_delta_t;
    last_renorm_nt   = 0;

    if(Particle_Dynamics)
      initialize_particles(exo, x, x_old, xdot, xdot_old, resid_vector);

    /*
     * Write out the initial solution to an ascii file
     * and to the exodus output file, if requested to do so by
     * an optional flag in the input file
     *  -> Helpful in debugging what's going on.
     */
    if (Write_Initial_Solution) {
      if (file != NULL)  { 
	error = write_ascii_soln(x, resid_vector, numProcUnknowns,
				 x_AC, nAC, time, file);
	if (error != 0)
	    DPRINTF(stderr,"%s:  error writing ASCII soln file\n", yo);
      }
      (void) write_solution(ExoFileOut, resid_vector, x, x_sens_p,
			    x_old, xdot, xdot_old, tev, tev_post, gv,
			    rd, gindex, p_gsize, gvec, gvec_elem,
			    &nprint, delta_t, theta, time, x_pp,
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
    if (nEQM > 0)
      {
        DPRINTF(stderr, "\nINITIAL ELEMENT QUALITY CHECK---\n");
        good_mesh = element_quality(exo, x, ams[0]->proc_config);
      }

    /*******************************************************************
     *  TOP OF THE TIME STEP LOOP -> Loop over time steps whether
     *                               they be successful or not
     *******************************************************************/
    for (n = 0; n < MaxTimeSteps; n++)
      {
      /*
       * Calculate the absolute time for the current step, time1
       */
      time1 = time + delta_t;
#ifdef LIBRARY_MODE
      delta_t_save = delta_t;
#endif
      if (time1 > TimeMax) { 
	DPRINTF(stderr, "\t\tLAST TIME STEP!\n"); 
	time1 = TimeMax;
	delta_t = time1 - time;
	tran->delta_t = delta_t;
	tran->delta_t_avg = 0.25*(delta_t+delta_t_old+delta_t_older
					+delta_t_oldest);
#ifdef LIBRARY_MODE
        last_step = TRUE;
#endif
      }
      tran->time_value = time1;
#ifdef LIBRARY_MODE
      if (n == (MaxTimeSteps-1) ) last_step = TRUE;
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
	for (w=0; w<efv->Num_external_field; w++) {
	  if (strcmp(efv->field_type[w], "transient") == 0)
	    {
	      err = rd_trans_vectors_from_exoII(x_old, efv->file_nm[w], w, n, &timeValueReadTrans, cx, dpi);
	      if (err != 0) {
		DPRINTF(stderr,
			"%s: err from rd_trans_vectors_from_exoII\n",yo);
	      }
	    }
	}
      }


      /* 
       * Get started with forward/Backward Euler predictor-corrector 
       * to damp out any bad things
       */
      if ( ( nt - last_renorm_nt) == 0) 
	{ 
	  theta	      = 0.0; 
	  const_delta_t = 1.0;
	  
	} 
      else if ( (nt - last_renorm_nt) >= 3) 
	{
	/* Now revert to the scheme input by the user */
	theta	      = tran->theta; 
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

      if ( nt == 0 )
        {
          xfem = NULL;
          if ( upd->XFEM )
            {
              xfem = alloc_struct_1(struct Extended_Shape_Fcn_Basics, 1);
              xfem->ielem = -1;
              xfem->tot_vol = alloc_dbl_1(numProcUnknowns, 0.0);
              xfem->active_vol =  alloc_dbl_1(numProcUnknowns, 0.0);
              if (ls == NULL)
                {
                  EH(-1,"Currently, XFEM requires traditional level set (not pf)");
                }
            }
        }
        
      /*
       * Initial Start up (t=0) for the FILL/LEVEL_SET equations
       */
    
      if (upd->ep[FILL] > -1  && nt == 0) 
	{ /*  Start of LS initialization */

#ifndef COUPLED_FILL
	/*
	 * Extract the fill unknowns into their own packed vector, xf[],
	 * at t = t_0
	 */
        if (Explicit_Fill) {
	  get_fill_vector(num_total_nodes, x_old, xf, node_to_fill);
	  dcopy1(num_fill_unknowns, xf, xf_old);
	}
#endif /* not COUPLED_FILL */

	if (ls != NULL || pfd != NULL) 
	  {
	
	  int eqntype	   = ls->Init_Method;
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
              EH(-1,"Level Set Evolution scheme not found \n");
          }


	  if ( ls->Length_Scale < 0.0 )
	      EH(-1, "\tError: a Level Set Length Scale needs to be specified\n");
            
	  if( ls->Integration_Depth > 0 || ls->SubElemIntegration || ls->AdaptIntegration )
	    {

              if ( ls->Integration_Depth > 0 )
                {
		  int first_elem;

		  first_elem = find_first_elem_with_var( exo, LS );
		  
		  if( first_elem != -1 )
		    {
                    load_ei(first_elem , exo, 0);

		      Subgrid_Tree = create_shape_fcn_tree ( ls->Integration_Depth );
		      DPRINTF(stdout,"\n\tSubgrid Integration of level set interface active.\n");
		    }
                }
              else if ( ls->SubElemIntegration )
                {
                  DPRINTF(stdout,"\n\tSubelement Integration of level set interface active.\n");
                }
              else if ( ls->AdaptIntegration )
                {
                  DPRINTF(stdout,"\n\tAdaptive Integration of level set interface active.\n");
                  DPRINTF(stdout,"\tAdaptive Integration Interface Order = %d\n",ls->Adaptive_Order);
                }
              Subgrid_Int.ip_total = 0;
              Subgrid_Int.s = NULL;
              Subgrid_Int.wt = NULL;
	    }

	  switch (eqntype) {
	  case  PROJECT :

	      DPRINTF(stderr,"\n\t Projection level set initialization \n");
#ifdef COUPLED_FILL
	      EH(-1,"Use of \"PROJECT\" is obsolete.");
#else /* COUPLED_FILL */
	      init_vec_value (xf, 0.0, num_fill_unknowns);
	      err = integrate_explicit_eqn(ams[FIL], rf, xf, xf_old, xfdot, 
					   xfdot_old, x, x_old, x_oldest,
					   step_size, theta, &time2, 
					   PROJECT, node_to_fill, exo, dpi, cx);
#endif /* COUPLED_FILL */

	      break;
		       
	  case EXO_READ :
		      
	      DPRINTF(stderr, "\t\t Level set read from exodus database \n");

	      break;

	  case SURFACES :

	      DPRINTF(stderr, "\n\t\t Surface object level set initialization : ");

	      /* parallel synchronization of initialization surfaces */
              if ( Num_Proc > 1 )
                {
                  if (!ls->init_surf_list) ls->init_surf_list = create_surf_list();
                  assemble_Global_surf_list( ls->init_surf_list );
                }
              
	      surf_based_initialization( x, NULL, NULL, exo, num_total_nodes,
					 ls->init_surf_list, 0., 0., 0. );

#ifdef DEBUG_PARALLEL
	      {
		FILE *data_out;
		int II, je;

		if(ProcID == 0 ) data_out = fopen("Proc0.dat", "w");
		if(ProcID == 1 ) data_out = fopen("Proc1.dat", "w");

		fprintf(data_out, "num_total_nodes:  %d \n", num_total_nodes);

		for( II=0; II < num_total_nodes; II++)
		  {
		    je =  Index_Solution(II, LS , 0, 0 , -1);

		    if( je != -1)
		      {
			fprintf( data_out, "%d  %lf \n", II, x[je] );
		      }
		  }

		fflush(data_out);
		fclose(data_out);
	      }
#endif
	      
	      DPRINTF(stderr, "- done \n");

#ifndef COUPLED_FILL
	      /*
	       * fill up xf and copy xf to xf_old.
	       * This is initialization so xf == xf_old is the best
	       * we can do.
	       */
              if (Explicit_Fill) {
	        get_fill_vector(num_total_nodes, x, xf, node_to_fill);
	        dcopy1(num_fill_unknowns, xf, xf_old);
	      }
#endif /* not COUPLED_FILL */

	      break;

	  case SM_OBJECT:
	    cgm_based_initialization(x, num_total_nodes);
#ifndef COUPLED_FILL
              if (Explicit_Fill)
		{
		  get_fill_vector(num_total_nodes, x, xf, node_to_fill);
		  dcopy1(num_fill_unknowns, xf, xf_old);
		}
#endif
	    break;

	  default:
	      WH(-1,"Level Set Initialization method not found \n");
	  } /* end of switch( eqntype )  */

	exchange_dof(cx, dpi, x);

	if (converged) 
	{
	  switch (ls->Renorm_Method) {

	  case HUYGENS :
	  case HUYGENS_C :
            Renorm_Now =  ( ls->Force_Initial_Renorm || (ls->Renorm_Freq != 0 && ls->Renorm_Countdown == 0) );

	    did_renorm = huygens_renormalization(x, num_total_nodes, exo, cx, dpi,  
						 num_fill_unknowns, numProcUnknowns,  
						 time1, Renorm_Now );
	    
#ifndef COUPLED_FILL
	    if ( did_renorm )
	      {
		get_fill_vector(num_total_nodes,x,xf,node_to_fill);
	      }
#endif /* not COUPLED_FILL */
#ifndef PHASE_COUPLED_FILL
	    if ( did_renorm )
	      {
		/*get_fill_vector(num_total_nodes,x,xf,node_to_fill);*/
	      }
#endif /* not COUPLED_FILL */

	      break;

	  case CORRECT :
#ifdef COUPLED_FILL
	    EH(-1,"Use of \"CORRECT\" is obsolete.");
#else /* COUPLED_FILL */
	    {
	      double step_size =  ls->Length_Scale/10.0;
	      int    num_steps = 15;
	      
	      put_fill_vector(num_total_nodes, x_old, xf, node_to_fill);
	      exchange_dof(cx, dpi, x_old);
	      put_fill_vector(num_total_nodes, x_oldest, xf, node_to_fill);
	      exchange_dof(cx, dpi, x_oldest);
	      correct_level_set(ams[FIL], xf, rf, x, x_old, x_oldest, node_to_fill,
				num_total_nodes, num_fill_unknowns, step_size,
				theta, num_steps, CORRECT, exo, dpi, cx);
	    }
#endif /* COUPLED_FILL */
	    break;	
	  default:
	      if ( ls->Evolution == LS_EVOLVE_ADVECT_EXPLICIT ||
                   ls->Evolution == LS_EVOLVE_ADVECT_COUPLED )
                 WH(-1,"No level set renormalization is on.\n");
	  } /* end of switch(ls->Renorm_Method ) */
	}
	  /*
	   * More initialization needed. Have to set those field variables that initially are indexed by
	   * level set function.  For example, species  concentration and temperature.
	   */
	  if( ls->Num_Var_Init > 0 ) 
	    ls_var_initialization ( x, exo, dpi, cx );

	  /* 	  DPRINTF(stderr, "Done with ls_var_initialization.\n"); */

	  /*
	   * Now check to see if we need to build a surface represent.
	   * on each time step.  Initialize the structures if so.
	   */
	  {
	    int build = FALSE, ibc = 0;

	    while ( !build && ibc < Num_BC)
	      {
		build = ( BC_Types[ibc].BC_Name == LS_INLET_BC );
		build = build || ( BC_Types[ibc].BC_Name == LS_ADC_BC );
		ibc++;
	      }
		  
#ifndef COUPLED_FILL
		build = build || ( ls->Evolution == LS_EVOLVE_SEMILAGRANGIAN );
#endif
		  
	/* Here we create space for an isosurface list that is updated
	 * every time step.  
	 */
	  
	    if( build && ls->last_surf_list == NULL ) 
	      {
		struct LS_Surf *tmp_surf = NULL;
		struct LS_Surf_Iso_Data *tmp_data;

		ls->last_surf_list = create_surf_list();

		tmp_surf = create_surf( LS_SURF_ISOSURFACE );
		tmp_data = ( struct LS_Surf_Iso_Data *) tmp_surf->data;
		tmp_data->isovar = FILL;
		tmp_data->isoval = 0.0;

		append_surf( ls->last_surf_list, tmp_surf );
	      }
            
	  } /* matches int build */

	}  /* end of ls != NULL */
      


#ifndef COUPLED_FILL
        if (Explicit_Fill) {
          put_fill_vector(num_total_nodes, x, xf, node_to_fill);
	}
#endif /* not COUPLED_FILL */
	dcopy1(numProcUnknowns, x, x_old);
	dcopy1(numProcUnknowns, x, x_older);
	dcopy1(numProcUnknowns, x, x_oldest);

	exchange_dof(cx, dpi, x);
	exchange_dof(cx, dpi, x_old);
	exchange_dof(cx, dpi, x_oldest);

	}


      ls_old = ls;
      if(upd->vp[PHASE1] > -1 && nt == 0 )
	{ /* Start of Phase Function initialization */
		  
		  
		  
		  if (pfd != NULL)
		  {
			struct Level_Set_Data *ls_save =ls;			
			if (upd->vp[PHASE1] > -1)
			  {
			    switch (pfd->ls[0]->Evolution) {
			    case LS_EVOLVE_ADVECT_EXPLICIT:
			      DPRINTF(stdout, "\n\t Using decoupled / subcycling for FILL equation for R_PHASE0.\n");
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
			      EH(-1,"PHASE Function Evolution scheme not found \n");
			      break;
			    }
			  }

			for (i = 0; i < pfd->num_phase_funcs;i++)
			  {	  
			    ls = pfd->ls[i];   /* This is a crucial step.  It allows us to use the Level Set Machinery 
						  to initialize the phase functions
						  BUT... Make sure to set it back to what we start with when through. 
						  Many places the test ls != NULL is used to thread the code 
					       */
				  
			    switch (ls->Init_Method)
			      {
			      case SURFACES:
				DPRINTF(stderr, "\n\t\t Surface object initialization for phase function: %d", i);

				/* parallel synchronization of initialization surfaces */
				if ( Num_Proc > 1 )
				  {
				    if (!ls->init_surf_list) ls->init_surf_list = create_surf_list();
				    assemble_Global_surf_list( ls->init_surf_list );
				  }
						  
				surf_based_initialization(x, NULL, NULL, exo, num_total_nodes,
							  ls->init_surf_list, 0., 0., 0. );
				break;
			      case EXO_READ:
				DPRINTF(stderr, "\n\t\t Exodus file read initialization for phase function fields");	  
				break;
			      } /* end of switch(ls->Init_Method ) */
			  } /* end of i<pfd->num_phase_funcs */

			ls = ls_save;  /* OK, used the level set routines now be nice and point 
					  ls back to where you found it */


			/* Allocate arrays for constraint Jacobian contributions */
			if (pfd->Use_Constraint == TRUE)
			  {
			    int i;
			    pfd->jac_info->d_pf_lm = alloc_dbl_1(numProcUnknowns, 0.0);
			    pfd->jac_info->d_lm_pf = alloc_dbl_1(numProcUnknowns, 0.0);
			    
			    if (pfd->ls[0]->Init_Method != EXO_READ)
			      {
				for (i = 0; i < pfd->num_phase_funcs; i++)
				  {
				    shift_nodal_values(PHASE1+i, -pfd->shift[i], x, num_total_nodes);
				  }
			      }
			  }



		  } /* end of pfd!=NULL */
	  ls = ls_old;


	  dcopy1(numProcUnknowns, x, x_old);
	  dcopy1(numProcUnknowns, x, x_older);
	  dcopy1(numProcUnknowns, x, x_oldest);
	  exchange_dof(cx, dpi, x);
	  exchange_dof(cx, dpi, x_old);
	  exchange_dof(cx, dpi, x_oldest);

	} /* end of phase function initialization */


	if( ls != NULL && ls->last_surf_list != NULL )
	  {
	    /* Find the interface surf at the last full time step (x_old = tmp_x)
	     * for use during this explicit step */

	    create_subsurfs(  ls->last_surf_list, x_old, exo );
	  }
      
#ifndef COUPLED_FILL
      /* 
       * Subcycling of the fill vector
       *  HKM -> Note, took out the dependence of the converged flag
       *         and switched initializations to x_old[]. I don't know
       *         what was done before for failed time step cases. I suspect
       *         this point was not covered, as the solution from the
       *         fill subcycling will obviously depend on delta_t which
       *         will have changed if a time step truncation error failure
       *         occurred.
       */
      if (Explicit_Fill && nt) {
	double *tmp_x, *tmp_x_old, *tmp_xdot;   /* temporary vector pointers */


	/*
	 * Gather the fill equation unknowns at t = t_n into a vector, xf[]
	 * and save them.  Note that we start from scratch each subcycling sequence
	 * So xf_old should be equal to xf and xfdot and xfdot_old should be zero 
	 * IMHO, it is unnecessary work to try to extract info from previous successul
	 * subcycling sequence, given that step sizes might have changed and the like (TAB)
	 */
	get_fill_vector(num_total_nodes, x_old, xf, node_to_fill);
	dcopy1(num_fill_unknowns, xf, xf_save);
	dcopy1(num_fill_unknowns, xf, xf_old);
	memset( xfdot, 0, sizeof(double)*num_fill_unknowns );
	memset( xfdot_old, 0, sizeof(double)*num_fill_unknowns );
	

	/*
	 * Calculate the delta_t that will be used in the subcycle
	 */
	delta_t_exp =  delta_t/ (double) exp_subcycle ;



	DPRINTF(stderr, 
		"\n\tsubcycling on Fill equation:  dt_exp = %e\n",
		delta_t_exp);
	log_msg("fill subcycle step %d, time = %g", n_exp, time2);

	/*
	 * Section to subcyle the fill equations
	 * -> we do this before advancing the other equations from n to n+1
	 *    We do this even before having a prediction of the other equations!
	 *    Therefore, it should be noted that x and xdot are not filled with
	 *    anything of significance at this point in the code.
	 *    -> a better algorithm may be to move the predictor ahead of this
	 *       subcycle section, in order to make use of predicted solution
	 *       in the subcycle.
	 */

	/* OK.  right before we start lets allocate some temporary vectors for holding
	   the fill equation as it evolves during the subcycling steps.  */

	tmp_x =  alloc_dbl_1(numProcUnknowns, 0.0);
	tmp_x_old =  alloc_dbl_1(numProcUnknowns, 0.0);
	tmp_xdot =  alloc_dbl_1(numProcUnknowns, 0.0);

	/* But this isn't confusing enough so let's copy x_old into tmp_x, 
	   xdot_old into tmp_xdot and x_older into tmp_x_old.  Note that at this point
	   in the code x_old, x_older contain "established" solutions which
	   we need later on in the integrate_explicit_eqn routine 
	*/

	dcopy1( numProcUnknowns,  x_old, tmp_x ) ;
	dcopy1( numProcUnknowns,  x_older, tmp_x_old ) ;
	dcopy1( numProcUnknowns,  xdot_old, tmp_xdot ) ;




	for (n_exp = 0; n_exp < exp_subcycle; n_exp++) {
	  /*
	   * always use BE for level set integration 
	   */
	  double _theta = (ls != NULL) ? 0.0 : theta;
	  /*
	   * time2 is the current time in the subcylce. The last time
	   * of the subcycle will be equal to time1, calculated above.
	   */
	  time2 = time + delta_t_exp * (n_exp + 1);

	  if ( ls == NULL || ls->Evolution == LS_EVOLVE_ADVECT_EXPLICIT )
	    {
	      DPRINTF(stderr,"\n\tsubcycling time step: %d  time = %e\n", 
		      n_exp, time2);


	  /*
	   *  Integrate the FILL equations for one time step
	   *  HKM -> At this point, for the case where the previous time
	   *         step was successful, x[] is equal to x_old[], except
	   *         for the FILL unknowns, which are being updated due
	   *         to the subcycling.
	   */
	      err = integrate_explicit_eqn(ams[FIL], rf, xf, xf_old, xfdot,
					   xfdot_old, tmp_x, tmp_x_old, tmp_xdot,
					   delta_t_exp, _theta, &time2, ADVECT,
					   node_to_fill,  exo, dpi, cx);

	      log_msg("fill subcycle step %d, time = %g", n_exp, time2);
	    }
	  else if ( ls->Evolution == LS_EVOLVE_SEMILAGRANGIAN )
	    {
	      DPRINTF(stderr,"\n\tsemi_lagrange time step: %d  time = %e\n", 
		      n_exp, time2);

	      semi_lagrange_step( num_total_nodes, numProcUnknowns, num_fill_unknowns, x, xf, xf_old,
				  xfdot, xfdot_old, node_to_fill, delta_t_exp, theta, exo, dpi, cx) ;

	      log_msg("semi_lagrange step %d, time = %g", n_exp, time2);
	    }


	      
	  /*
	   * Update the main solution vector with the new FILL values
	   * and their derivatives
	   */
	  put_fill_vector(num_total_nodes, tmp_x, xf, node_to_fill);
	  exchange_dof(cx, dpi, tmp_x);
	  put_fill_vector(num_total_nodes, tmp_x_old, xf, node_to_fill);
	  exchange_dof(cx, dpi, tmp_x_old);
	  put_fill_vector(num_total_nodes, tmp_xdot, xfdot, node_to_fill);
	  exchange_dof(cx, dpi, tmp_xdot);

	  dcopy1(num_fill_unknowns, xf, xf_old);

	  dcopy1(num_fill_unknowns, xfdot, xfdot_old);
	}

	/* Now ditch the temporary vectors */

	safer_free( ( void **) &tmp_x );
	safer_free( ( void **) &tmp_x_old );
	safer_free( ( void **) &tmp_xdot );

	
      } /* End of the FILL variable subcycling in time */
	  

#endif /* not COUPLED_FILL */

      if (ProcID == 0) {
	if (theta == 0.0)
	    strcpy(tspstring, "(BE)");
	else if (theta == 0.5)
	    strcpy(tspstring, "(CN)");
	else if (theta == 1.0)
	    strcpy(tspstring, "(FE)");
	else
	    sprintf(tspstring, "(TSP %3.1f)", theta);
	fprintf(stderr,
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

      predict_solution(numProcUnknowns, delta_t, delta_t_old,
		       delta_t_older,  theta, x, x_old, x_older, 
		       x_oldest, xdot, xdot_old, xdot_older);
      
      if(tran->solid_inertia)
	{
	  predict_solution_newmark(num_total_nodes, delta_t, x, x_old, xdot, xdot_old);
	  exchange_dof(cx, dpi, tran->xdbl_dot);
	}

#ifdef LASER_RAYTRACE
      if (ls != NULL)
	{
	  double (*point0)[DIM] = NULL;
          double (*point1)[DIM] = NULL;
          int * owning_elem = NULL;
          int facet, num_facets;
          
          num_facets = generate_facet_list(&point0, &point1, &owning_elem, x, exo);
          for ( facet = 0; facet < num_facets; ++facet )
            {
              fprintf(stderr, "FACET %d: owning element = %d, First point = (%g,%g), Second point = (%g,%g)\n", 
                      facet,owning_elem[facet],point0[facet][0],point0[facet][1],point1[facet][0],point1[facet][1]);
            }
	}
#endif
      
      if (ls != NULL && ls->Evolution == LS_EVOLVE_SLAVE)
	{
	  surf_based_initialization(x, NULL, NULL, exo, num_total_nodes,
				       ls->init_surf_list, time1, theta, delta_t);
	}
      if (pfd != NULL)
        {
	  ls_old = ls;  /*the ol' switcheroo*/
          ls = pfd->ls[0];
	  if (ls->Evolution == LS_EVOLVE_SLAVE)
	    {
	      surf_based_initialization(x, NULL, NULL, exo, num_total_nodes,
				       ls->init_surf_list, time1, theta, delta_t);
	    }
          ls = ls_old;
        }

      /* Now go back and correct all those dofs using XFEM */
      if(xfem != NULL)
        {
          xfem_predict( num_total_nodes, numProcUnknowns, delta_t, delta_t_old,
		        delta_t_older,  theta, x, x_old, x_older,
		        x_oldest, xdot, xdot_old, xdot_older);
        }
        
      /* 
       * Now, that we have a predicted solution for the current
       * time, x[], exchange the degrees of freedom to update the
       * ghost node information.
       */
      exchange_dof(cx, dpi, x);
      exchange_dof(cx, dpi, xdot);
        
#ifdef DEBUG
      if (nt == 0) {
	print_array(x, numProcUnknowns, "x_A", type_double, ProcID);
      }
#endif /* DEBUG */

      if (nAC > 0) {

	predict_solution(nAC, delta_t, delta_t_old, delta_t_older,
			 theta, x_AC, x_AC_old, x_AC_older, x_AC_oldest,
			 x_AC_dot, x_AC_dot_old, x_AC_dot_older);

	for(iAC = 0; iAC < nAC; iAC++)
	  {
	    update_parameterAC(iAC, x, xdot, x_AC, cx, exo, dpi);
 	    augc[iAC].tmp2 = x_AC_dot[iAC];
 	    augc[iAC].tmp3 = x_AC_old[iAC];
	  }
      }

      /*
       *  Set dirichlet conditions in some places. Note, I believe
       *  this step can change the solution vector
       */
      find_and_set_Dirichlet(x, xdot, exo, dpi); 


#ifndef COUPLED_FILL
      /* 
       * Now, that we have a predicted solution for the current
       * time, x[], exchange the degrees of freedom to update the
       * ghost node information.
       */
      if (Explicit_Fill) {
	/* The predicted fill values from the previous call are not need 
	   since ( in theory ) we have updated the fill field during 
	   the explicit loop.  Here we overwrite the predicted fill values;
	*/

	put_fill_vector(num_total_nodes, x, xf, node_to_fill);
      }
#endif /* not COUPLED_FILL */

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

      exchange_dof(cx, dpi, x);
      exchange_dof(cx, dpi, xdot);
      if(tran->solid_inertia)  exchange_dof(cx, dpi, tran->xdbl_dot);

      /*
       * Save the predicted solution for the time step
       * norm calculation to be carried out after convergence
       * of the nonlinear implicit problem
       */
      dcopy1(numProcUnknowns, x, x_pred);
      if ( nAC > 0 ) dcopy1( nAC, x_AC, x_AC_pred );

#ifdef RESET_TRANSIENT_RELAXATION_PLEASE  
  /* Set TRUE to disable relaxation on timesteps after the first*/
      /* For transient, reset the Newton damping factors after a
       *   successful time step
       */
      if (nt > 0)
      {
              if (converged && (nt - last_renorm_nt) > 0 )
              { damp_factor2 = -1.;  damp_factor1 = 1.0;}
              else
              {
               damp_factor2 = damp_factor_org[1];
               damp_factor1 = damp_factor_org[0];
              }
      }
#endif
      /*
       *  Solve the nonlinear problem. If we achieve convergence,
       *  set the flag, converged, to true on return. If not
       *  set the flag to false.
       */
      err = solve_nonlinear_problem(ams[JAC], x, delta_t, theta, x_old, 
				    x_older, xdot, xdot_old, resid_vector,  
				    x_update, scale, &converged, &nprint,
				    tev, tev_post, gv, rd, gindex, p_gsize, 
				    gvec, gvec_elem, time1, exo, dpi, cx, 
				    n, &time_step_reform, is_steady_state,
 				    x_AC, x_AC_dot, time1, resid_vector_sens,
				    x_sens, x_sens_p, NULL);
      if (err == -1) converged = FALSE;
      inewton = err;
      evpl_glob[0]->update_flag = 0; /*See get_evp_stress_tensor for description */
      af->Sat_hyst_reevaluate = FALSE; /*See load_saturation for description*/

#ifdef DEBUG
      if (nt == 0) {
	print_array(x, numProcUnknowns, "x_B", type_double, ProcID);
      }
#endif /* DEBUG */

      /*
       * HKM -> I do not know if these operations are needed. I added
       *        an exchange of xdot[] here, because if x[] is exchanged
       *        then xdot needs to be exchanged as well.
       */

      exchange_dof(cx, dpi, x);
      exchange_dof(cx, dpi, xdot);

      if (converged) af->Sat_hyst_reevaluate = TRUE;  /*see load_saturation */

      /* Check element quality */
      good_mesh = element_quality(exo, x, ams[0]->proc_config);

      /*
       * Check the time step truncation error.
       * If too large, set the success_dt flag to FALSE. We will
       * then not accept the current time step, reduced delta_t,
       * and retry.
       */
      if (converged) {
	
	delta_t_new = time_step_control(delta_t, delta_t_old, const_delta_t,
					x, x_pred, x_old, x_AC, x_AC_pred,
					eps, &success_dt, tran->use_var_norm);
	if (const_delta_t) {
	  success_dt  = TRUE;
	  delta_t_new = delta_t;
	} else if ( failed_recently_countdown > 0 ) {
          delta_t_new = delta_t;
          failed_recently_countdown--;
	} else if (delta_t_new > Delta_t_max) {
	  delta_t_new = Delta_t_max;
        } else if ( !success_dt && delta_t_new < tran->resolved_delta_t_min ) {
          if ( delta_t > tran->resolved_delta_t_min ) {
            /* fool algorithm into using delta_t = tran->resolved_delta_t_min */
            delta_t = tran->resolved_delta_t_min;
            delta_t /= tran->time_step_decelerator;
	    tran->delta_t  = delta_t;
	    tran->delta_t_avg = 0.25*(delta_t+delta_t_old+delta_t_older
					+delta_t_oldest);
          } else {
            /* accept any converged solution with
               delta_t <= tran->resolved_delta_t_min */
            success_dt = TRUE;
            delta_t_new = delta_t;
          }
        }
        
        if ( ls != NULL && tran->Courant_Limit != 0. ) {
          double Courant_dt;
          Courant_dt = tran->Courant_Limit *
                       Courant_Time_Step( x, x_old, x_older, xdot, xdot_old,
                                        resid_vector, ams[0]->proc_config, exo );
          if ( Courant_dt > 0. && Courant_dt < delta_t_new ) {
            DPRINTF(stderr,"\nCourant Limit requires dt <= %g\n",Courant_dt);
            delta_t_new = Courant_dt;
          }
        }
        
      }

      if (converged && success_dt) {
	if (Filter_Species) {
	  err = filter_conc(num_total_nodes, x, filter_species_material_number, 
			    c_min, c_max ); 
	}
	nt  += 1;
	time = time1;
  
	/* Determine whether to print out the data or not */
	i_print = 0;
	if (tran->print_freq == 0) {
	  if ((time > time_print) || 
	      (fabs(time - time_print) < (1.e-4 * tran->print_delt)))
	  {
	    if (tran->print_delt2 < 0.) {
	      i_print	  = 1;
	      time_print += tran->print_delt;
	    } else {
	      if (time < tran->print_delt2_time) {
		i_print	    = 1;
		time_print += tran->print_delt;
	      } else {
		i_print	    = 1;
		time_print += tran->print_delt2;
	      }
	    }
	  }
	} else {
	  if (nt == step_print) {
	    i_print	= 1;
	    step_print += tran->print_freq;
	  }
	}

	if (time1 >= (ROUND_TO_ONE * TimeMax)) i_print = 1;

	/* Dump out user specified information to separate file.
	 */
	err = usr_print(&time, delta_t, x, NULL, -1); 
	EH(err, "usr_print");

	/* Particle calculations.  time = time at *beginning* of
	 * current timestep, n = timestep. */ 
	if(Particle_Dynamics)
	  err = compute_particles(exo, x, x_old, xdot, xdot_old, resid_vector, time, delta_t, n);
	EH(err, "Error performing particle calculations.");

	error = 0;
	if (i_print) {
	  if (file != NULL) {
	    error = write_ascii_soln(x, resid_vector, numProcUnknowns,
				     x_AC, nAC, time, file);
	    if (error != 0) {
	      fprintf(stderr, 
		      "%s:  error writing ASCII soln file\n", yo);
	    }
	  }
	  if (Write_Intermediate_Solutions == 0) {
#ifdef LIBRARY_MODE
            lib_print = FALSE;
            if (libio->goma_first == 0)
              {
                if (libio->print_flag == 0 && last_call && last_step)
                lib_print = TRUE;
                if (libio->print_flag == 1 && last_step) lib_print = TRUE;
                if (libio->print_flag == 2) lib_print = TRUE;
              fprintf(stderr, "JAS FIRST...animas_step = %d, lib_print = %d\n",
                      libio->animas_step, lib_print);
              }
            else
              {
                if (libio->print_flag == 2 && !last_step) lib_print = TRUE;
              }

            if ( lib_print )
              /* Write out the full solution */
              {
                write_solution(ExoFileOut, resid_vector, x, x_sens_p,
                               x_old, xdot, xdot_old, tev, tev_post, gv,
                               rd, gindex, p_gsize, gvec, gvec_elem,
                               &nprint, delta_t, theta, time, x_pp, exo, dpi);
                nprint++;
              }
            else if (Num_Export_XP > 0)
              /* Just update the post-processing fields */
              {
                post_process_nodal(x, x_sens_p, x_old, xdot, xdot_old,
                                   resid_vector, nprint+1, &time_value,
                                   delta_t, theta, x_pp, exo, dpi, rd, NULL);
              }
#else
	    write_solution(ExoFileOut, resid_vector, x, x_sens_p,
			   x_old, xdot, xdot_old, tev, tev_post, gv,
			   rd, gindex, p_gsize, gvec, gvec_elem,
			   &nprint, delta_t, theta, time, x_pp, exo, dpi);
	    nprint++;
#endif

          if (ls !=NULL && ls->Interface_Output == TRUE)
	      {
              print_point_list(x, exo, ls->output_file, time);
	      }

	  }
	  /* Print out values of extra unknowns from augmenting conditions */
	  if (nAC > 0) {
	    DPRINTF(stderr, "\n------------------------------\n");
	    DPRINTF(stderr, "Augmenting Conditions:    %4d\n", nAC);
	    DPRINTF(stderr, "Number of extra unknowns: %4d\n\n", nAC);
		      
	    for (iAC = 0; iAC < nAC; iAC++) {
	      if (augc[iAC].Type == AC_USERBC) 
		{
		  DPRINTF(stderr, "\tBC[%4d] DF[%4d]=% 10.6e\n", augc[iAC].BCID, augc[iAC].DFID, x_AC[iAC]);
		}
	      else if (augc[iAC].Type == AC_USERMAT || augc[iAC].Type == AC_FLUX_MAT)
		{
		  DPRINTF(stderr, "\tMT[%4d] MP[%4d]=% 10.6e\n", augc[iAC].MTID, augc[iAC].MPID, x_AC[iAC]);
		}
	      else if (augc[iAC].Type == AC_VOLUME)
		{
		  evol_local = augc[iAC].evol;
#ifdef PARALLEL
		  if (Num_Proc > 1) {
		    MPI_Allreduce(&evol_local, &evol_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);		  
		    evol_local = evol_global;
		  }
#endif	  
		  DPRINTF(stderr, "\tMT[%4d] VC[%4d]=%10.6e Param=%10.6e\n", augc[iAC].MTID, augc[iAC].VOLID, evol_local, x_AC[iAC]);
		} 
	      else if (augc[iAC].Type == AC_FLUX) 
		{
		  DPRINTF(stderr, "\tBC[%4d] DF[%4d]=%10.6e\n", augc[iAC].BCID, augc[iAC].DFID, x_AC[iAC]);
		}
	      else if (augc[iAC].Type == AC_POSITION)
		{
		  DPRINTF(stderr, "\tNodeSet[%4d]_Pos = %10.6e F_bal = %10.6e VC[%4d] Param=%10.6e\n", 
			  augc[iAC].MTID, augc[iAC].evol, augc[iAC].lm_resid, augc[iAC].VOLID, x_AC[iAC]);
		}
	    }
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
          step_fix += tran->fix_freq*tran->print_freq;
        }
	/* 
	 * Adjust the time step if the new time will be larger than the
	 * next printing time.
	 */
	if (tran->print_freq == 0 && success_dt) {
	  if ((time + 1.2 * delta_t_new >= time_print)
	      && (time_print > time)) {
	    delta_t_new = time_print - time;
	    DPRINTF(stderr, 
		    "reset delta_t = %g to maintain printing frequency\n"
		    , delta_t_new);
	    if (delta_t_new <= 0) 
		EH(-1, "error with time-step printing control");
	  } else if(time >= time_print) {
	    if (delta_t_new != tran->print_delt) {
	      delta_t_new = tran->print_delt;
	      DPRINTF(stderr, 
		      "reset delta_t = %g to maintain printing frequency\n"
		      , delta_t_new);
	      if (delta_t_new <= 0) {
		EH(-1, "error with time-step printing control");
	      }
	    }
	  }
	}

#ifdef LIBRARY_MODE
        /*
         * When porosity is passed in from another code, do a
         * routine update to the external field.
         */
        if (update_porosity && converged)
          {
            error = advance_porosity_ev(n, num_total_nodes,
                                        x, base_p_por, base_p_liq);
            EH(error, "Problem with advance_porosity_ev()!");
          }
#endif


	if (converged && ls != NULL) 
	{
	  int ibc, ls_adc_event;
		/* Resolve LS_ADC boundaries ( Attach/Dewet/Coalesce ) */
		for( ibc=0;ibc<Num_BC;ibc++)
		{
		  ls_adc_event = FALSE;
			switch ( BC_Types[ibc].BC_Name )
			{
			case LS_ADC_OLD_BC:
				resolve_ls_adc_old ( &(BC_Types[ibc]), exo, x, delta_t, &ls_adc_event, nt );
				break;
			case LS_ADC_BC:
				resolve_ls_adc ( ls->last_surf_list->start->subsurf_list, &(BC_Types[ibc]), exo, x, delta_t, &ls_adc_event, nt );

				break;
			default:
				break;
			}
#ifdef PARALLEL
		    if (ls_adc_event) 
		      {
			exchange_dof( cx, dpi, x );
		      }
#endif
                  if ( ls_adc_event && tran->Restart_Time_Integ_After_Renorm )
                         {
                      /* like a restart */
                      discard_previous_time_step(numProcUnknowns, x, x_old,x_older,
                                      x_oldest, xdot,xdot_old,xdot_older);
                      last_renorm_nt = nt;
                      if ( delta_t_new > fabs(Delta_t0) )
                              delta_t_new *= tran->time_step_decelerator;
                        }

		}

		/* Check for renormalization  */


		ls->Renorm_Countdown -= 1;
		
		switch (ls->Renorm_Method) 
		{
			
			case HUYGENS:
			case HUYGENS_C:
				
				Renorm_Now = ( ls->Renorm_Freq != 0 && 
					ls->Renorm_Countdown == 0 ) 
					|| ls_adc_event == TRUE;
				
				did_renorm = huygens_renormalization(x, num_total_nodes,
					 exo, cx, dpi, num_fill_unknowns, numProcUnknowns,
					  time2, Renorm_Now);
				if ( did_renorm )
					{ exchange_dof(cx, dpi, x);	}
				if ( did_renorm && tran->Restart_Time_Integ_After_Renorm )
				{
					/* like a restart */
					discard_previous_time_step(numProcUnknowns, x, x_old,x_older,x_oldest, xdot,xdot_old,xdot_older);
					last_renorm_nt = nt;
					if ( delta_t_new > fabs(Delta_t0) ) delta_t_new *= tran->time_step_decelerator;
				}

					break;
				
			case CORRECT:
				EH(-1,"Use of \"CORRECT\" is obsolete.");
				break;
			default:
				break;
		}
		
		if( ls->Sat_Hyst_Renorm_Lockout > 0 ) {
		  af->Sat_hyst_reevaluate = FALSE;
		  ls->Sat_Hyst_Renorm_Lockout -= 1;
		}
		
	}

	if( converged && pfd != NULL )
	  {
	    struct Level_Set_Data *ls_save = ls;
		
	    for( i=0; i<pfd->num_phase_funcs ; i++)
	      {
			
		ls = pfd->ls[i];
		ls->Renorm_Countdown -= 1;
		switch (ls->Renorm_Method) 
		  {
				
		  case HUYGENS:
		  case HUYGENS_C:
					
		    Renorm_Now = ( ls->Renorm_Freq != 0 && 
				   ls->Renorm_Countdown == 0 );
				
		    did_renorm = huygens_renormalization(x, num_total_nodes, exo,
							 cx, dpi, num_fill_unknowns, numProcUnknowns,
							 time2, Renorm_Now);
		    if ( did_renorm )
		      { exchange_dof(cx, dpi, x); }
		    if ( did_renorm && tran->Restart_Time_Integ_After_Renorm )
		      {
			/* like a restart */
			discard_previous_time_step(numProcUnknowns, x, x_old,x_older,
						   x_oldest, xdot,xdot_old,xdot_older);
			last_renorm_nt = nt;
			if ( delta_t_new > fabs(Delta_t0) ) 
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
	if(tran->solid_inertia)
	  dcopy1(numProcUnknowns, tran->xdbl_dot,    tran->xdbl_dot_old);
	dcopy1(numProcUnknowns, xdot,     xdot_old);
	dcopy1(numProcUnknowns, x_older,  x_oldest);
	dcopy1(numProcUnknowns, x_old,    x_older);
	dcopy1(numProcUnknowns, x,        x_old);
	delta_t_oldest = delta_t_older;
	delta_t_older  = delta_t_old;
	delta_t_old    = delta_t;
        tran->delta_t_old = delta_t_old;
        tran->time_value_old = time;
	delta_t        = delta_t_new;
	tran->delta_t  = delta_t;  /*load up for use in load_fv_mesh_derivs*/
	tran->delta_t_avg = 0.25*(delta_t+delta_t_old+delta_t_older
					+delta_t_oldest);

	if (nAC > 0) {
	  dcopy1(nAC, x_AC_dot_old, x_AC_dot_older);
	  dcopy1(nAC, x_AC_dot,     x_AC_dot_old);
	  dcopy1(nAC, x_AC_older,   x_AC_oldest);
	  dcopy1(nAC, x_AC_old,     x_AC_older);
	  dcopy1(nAC, x_AC,         x_AC_old);
	}

	/* Integrate fluxes, forces  
	 */
	for (i = 0; i < nn_post_fluxes; i++) {
	  (void) evaluate_flux(exo, dpi, 
			       pp_fluxes[i]->ss_id, 
			       pp_fluxes[i]->flux_type ,
			       pp_fluxes[i]->flux_type_name ,
			       pp_fluxes[i]->blk_id , 
			       pp_fluxes[i]->species_number, 
			       pp_fluxes[i]->flux_filenm,
			       pp_fluxes[i]->profile_flag,
			       x, xdot, NULL, delta_t, time, 1); 
	}

	/* Compute flux, force sensitivities
	 */
	for (i = 0; i < nn_post_fluxes_sens; i++) {
	  (void) evaluate_flux_sens(exo, dpi,
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
				    x,xdot,x_sens_p,delta_t,time,1);
	}
#if 1
	/* This section is a kludge to determine minimum and maximum
		variable values for outputing surfaces of constant concentration
	*/
	if( nn_volume )
	  {
	  int search_minmax = FALSE;
 	  for (i = 0; i < nn_volume; i++) {
 		if(pp_volume[i]->volume_type == I_SURF_SPECIES &&
				pp_volume[i]->species_no == -1 ) 
		search_minmax=TRUE;
		}
	  if( search_minmax && Num_Proc > 1)
		  WH(-1,"Species surfaces not recommended in parallel\n");
	  if( search_minmax && Num_Proc == 1)
 	  {
 		  int inode,offset,idof;
 		  VARIABLE_DESCRIPTION_STRUCT *vd;
 		  double specmax[MAX_CONC],specmin[MAX_CONC],specavg, fraction;
		  int eq_off = 0;

		  if(pd->e[R_ENERGY] )  eq_off++;
		  if(pd->e[R_MOMENTUM1] )  eq_off += pd->Num_Dim;
 		  inode = 0;
 		  for(i=0 ; i < (num_internal_dofs+num_boundary_dofs); i++) {
 			  vd = Index_Solution_Inv(i, &inode, NULL, &offset, &idof);
 			  if( vd->Variable_Type == MASS_FRACTION )      {
 				if( i >= pd->Num_Species_Eqn + eq_off ) {
 		                  if (x[i] > specmax[offset-eq_off]) 
						specmax[offset-eq_off] = x[i];
 				  if (x[i] < specmin[offset-eq_off]) 
						specmin[offset-eq_off] = x[i];
 		                  } else  {
 	                          specmin[offset-eq_off] = x[i];
 	                          specmax[offset-eq_off] = x[i];
 		                  }
 		          }
 		  }
 		  DPRINTF(stderr,"Species concentration minmax: %d \n",ProcID);
 	  	  for (i = 0; i < pd->Num_Species_Eqn; i++) {
 			DPRINTF(stderr,"%d %g %g \n",i,specmin[i],specmax[i]);
 			}
 	  for (i = 0; i < nn_volume; i++) {
 		if(pp_volume[i]->volume_type == I_SURF_SPECIES) 
 			{
         DPRINTF(stderr,"params %d %d %d\n",i,pp_volume[i]->num_params,pd->Num_Species);
           		if(pp_volume[i]->num_params < pd->Num_Species+1)
                		{
         DPRINTF(stderr,"params %d %d\n",pp_volume[i]->num_params,pd->Num_Species);
         WH(-1,"SURF_SPECIES parameters should number Num_Species+1");
                		}
 			if( pp_volume[i]->species_no == -1 )
 				{
				fraction = pp_volume[i]->params[pd->Num_Species];
				pp_volume[i]->params[pd->Num_Species] = 
					pp_volume[i]->params[pd->Num_Species-1];
				for(idof=0 ; idof<pd->Num_Species_Eqn; idof++)
					{
					specavg = specmin[idof] + 
					  fraction*(specmax[idof]-specmin[idof]);
					pp_volume[i]->params[pd->Num_Species] +=
					  pp_volume[i]->params[idof]*specavg;
					}
				DPRINTF(stderr,"new surface value = %g \n",pp_volume[i]->params[pd->Num_Species]);
				}
 			}
 		}	/* nn_volume loop	*/
 	  }	/*search  */
	}	/*  nn_volume	*/
#endif
 	  for (i = 0; i < nn_volume; i++) {
 	    evaluate_volume_integral(exo, dpi, pp_volume[i]->volume_type,
 				     pp_volume[i]->volume_name,
 				     pp_volume[i]->blk_id,
 				     pp_volume[i]->species_no,
 				     pp_volume[i]->volume_fname,
 				     pp_volume[i]->params,
 				     NULL, x, xdot, delta_t, time, 1);
 	  }


	if (time1 >= (ROUND_TO_ONE * TimeMax))  {
	  DPRINTF(stderr,"\t\tout of time!\n");
     	  if (Anneal_Mesh)
	    {
	      /*
	       * Transform the node point coordinates according to the
	       * displacements and write out all the results using the
	       * displaced coordinates. Set the displacement field to
	       * zero, too.
	       */
	      err = anneal_mesh(x, tev, tev_post, gv, rd, time1, exo, dpi);
	      EH(err, "anneal_mesh() bad return.");
	    }
	  goto free_and_clear;
	}
        if (!good_mesh) goto free_and_clear;

      } /*  if(converged && success_dt) */
	  
      else /* not converged or unsuccessful time step */
      {
/* Set bit TRUE in next line to enable retries for failed first timestep*/
        if(relax_bit && nt == 0 && n < 15) {
        DPRINTF(stderr,"\nHmm... could not converge on first step\n Let's try some more iterations\n");
             if(inewton == -1)        {
                       damp_factor1 *= 0.5;
                  }   else  {
                       damp_factor1 *= 2.0;
                     damp_factor1 = MIN(damp_factor1,1.0);
                     dcopy1(numProcUnknowns, x,x_old);
                       if (nAC > 0) {
                                 dcopy1(nAC, x_AC,       x_AC_old);
                                 }
                  }
        DPRINTF(stderr,"\tdamping factor %g \n",damp_factor1);
           } else  {
	DPRINTF(stderr,"\n\tlast time step failed, dt *= %g for next try!\n",
		tran->time_step_decelerator);
	      
	delta_t *= tran->time_step_decelerator;
	tran->delta_t  = delta_t;
	tran->delta_t_avg = 0.25*(delta_t+delta_t_old+delta_t_older
					+delta_t_oldest);
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
        if ( ls != NULL && ls->SubElemIntegration ) {
	  DPRINTF(stderr,"predicted subelement configuration:\n");
	  subelement_mesh_output(x_pred, exo);
	  DPRINTF(stderr,"failed solution subelement configuration:\n");
	  subelement_mesh_output(x, exo);
	}
#endif

#ifndef COUPLED_FILL
	/*
	 *  HKM NOTE: I am not sure why the attempt below to fix up
	 *        x[] is being made. For a failed time step, x[]
	 *        will contain garbage at this point. To retry the step
	 *        we will discard the time step by restarting with
	 *        x_old[] at the top of the time step loop. Therefore,
	 *        I advocate eliminating the next block of code.
	 */
	if (Explicit_Fill && nt) {
	  /* restore original values for failed time step */
	  put_fill_vector(num_total_nodes, x, xf_save, node_to_fill);
	  exchange_dof(cx, dpi, x);
	  put_fill_vector(num_total_nodes, x_old, xf_save, node_to_fill);
	  exchange_dof(cx, dpi, x_old);
	}
#endif /* not COUPLED_FILL */
      }
      
      if (delta_t <= Delta_t_min) {
	DPRINTF(stderr,"\n\tdelta_t = %e < %e\n\n",delta_t, Delta_t_min);
	
	DPRINTF(stderr,"time step too small, I'm giving up!\n");
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
	EH(err, "anneal_mesh() bad return.");
      }
       */
      
  } /* end of if steady else transient */
  
 free_and_clear:

/* If exporting variables to another code, save them now! */
#ifdef LIBRARY_MODE
  callnum++;
  *te_out = time;
  if (Num_Export_XS > 0 || Num_Export_XP > 0)
    {
      err = load_export_vars(num_total_nodes, x, x_pp);
      EH(err, "Problem with saving export variables");
    }
#endif

  /* Free a bunch of variables that aren't needed anymore */

  for (i = 0; i < Num_ROT; i++)
    {
      safer_free((void **) &(ROT_Types[i].elems));
    }
  
  /*
   * Curiously, these bananas were allocated end to end - one free does it all.
   */

  safer_free((void **) &ROT_Types);
  safer_free((void **) &node_to_fill);
  
  safer_free( (void **) &resid_vector); 
  safer_free( (void **) &resid_vector_sens);
  safer_free( (void **) &scale); 
  if (last_call) safer_free( (void **) &x); 
  if (last_call) safer_free( (void **) &xp_id); 
  safer_free( (void **) &gv); 
  safer_free( (void **) &x_pp); 
  if (update_porosity) {
    safer_free( (void **) &base_p_por);
    safer_free( (void **) &base_p_liq);
  }

  if (nAC > 0 && last_call) {
    safer_free( (void **) &x_AC);
    safer_free( (void **) &x_AC_old);
    safer_free( (void **) &x_AC_older);
    safer_free( (void **) &x_AC_oldest);
    safer_free( (void **) &x_AC_dot);
    safer_free( (void **) &x_AC_dot_old);
    safer_free( (void **) &x_AC_dot_older);
    safer_free( (void **) &x_AC_pred);
  }

  safer_free((void **) &x_pred); 

  if (last_call)
    {
      safer_free((void **) &x_old); 
      safer_free((void **) &x_older); 
      safer_free((void **) &x_oldest); 
      safer_free((void **) &xdot); 
      safer_free((void **) &xdot_old);
      safer_free((void **) &xdot_older);
    }

  if(tran->solid_inertia)
    {
      free(tran->xdbl_dot);
      free(tran->xdbl_dot_old);
    }
  safer_free((void **) &x_update); 

  if (Linear_Solver == FRONT) {
    free(fss->bc);
    free(fss->ncod);
    free(fss->ncn);
    free(fss->ncn_base);
    free(fss->mdf);
    free(fss->nopp);
    free(fss->el_proc_assign);
    free(fss->constraint);
    free(fss->level);
    free(fss->perm);
    free(fss->iperm);
    free(fss->mask);
  }

  if ((nn_post_data_sens + nn_post_fluxes_sens) > 0) {
    safer_free((void **) &x_sens);
    Dmatrix_death(x_sens_p,num_pvector,numProcUnknowns);
  }

  if (last_call)
    {
      for(i = 0; i < MAX_NUMBER_MATLS; i++)
        {
        for(n = 0; n < MAX_MODES; n++)
          {
            safer_free((void **) &(ve_glob[i][n]->gn)); 
            safer_free((void **) &(ve_glob[i][n])); 
          }
        safer_free((void **) &(vn_glob[i])); 
      }
  }

  safer_free((void **) &ve); 

  if(num_universe_dofs != (num_internal_dofs + num_boundary_dofs) ) {
    safer_free((void **) &ija_attic);
  }

#ifndef COUPLED_FILL
  if (Explicit_Fill) {
    safer_free((void **) &rf); 
    safer_free((void **) &xf); 
    safer_free((void **) &xf_save); 
    safer_free((void **) &xf_old);
    safer_free((void **) &xfdot); 
    safer_free((void **) &xfdot_old);
    if (dpi->num_universe_nodes != 
	(dpi->num_internal_nodes + dpi->num_boundary_nodes) )
	safer_free((void **) &ijaf_attic);
  }
#endif /* not COUPLED_FILL */

  if (last_call) {
    if (strcmp(Matrix_Format, "epetra") == 0)
    {
      EpetraDeleteRowMatrix(ams[JAC]->RowMatrix);
      if (ams[JAC]->GlobalIDs != NULL)
      {
        free(ams[JAC]->GlobalIDs);
      }
    }

    sl_free(matrix_systems_mask, ams);

    for(i = 0; i < NUM_ALSS; i++)
      {
        safer_free((void **) &(ams[i]));
      }

    safer_free((void **) &gvec);

    i = 0;
    for (eb_indx = 0; eb_indx < exo->num_elem_blocks; eb_indx++)
      {
        for (ev_indx = 0; ev_indx < rd->nev; ev_indx++)
          {
            if (exo->elem_var_tab[i++] == 1)
              {
                safer_free((void **) &(gvec_elem[eb_indx][ev_indx]));
              }
          }
        safer_free((void **) &(gvec_elem[eb_indx]));
      }

    safer_free((void **) &gvec_elem);
    safer_free((void **) &rd);
    safer_free((void **) &Local_Offset);
    safer_free((void **) &Dolphin);
  }

if( ls != NULL )
{
   if ( ls->init_surf_list != NULL ) free_surf_list( &(ls->init_surf_list) );
}

free_shape_fcn_tree( Subgrid_Tree );

  if (file != NULL) fclose(file);
 
#ifdef DEBUG
  fprintf(stderr, "%s: leaving solve_problem()\n", yo);
#endif /* DEBUG */

  return;
} /* END of routine solve_problem()  */
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

static void

predict_solution(int N, double delta_t, double delta_t_old,
		 double delta_t_older, double theta_arg, double x[],
		 double x_old[], double x_older[], double x_oldest[],
		 double xdot[], double xdot_old[], double xdot_older[])

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
    c1 = delta_t * (delta_t + delta_t_old) / delta_t_older /
      (delta_t_older + delta_t_old);
    c2 = -delta_t * (delta_t + delta_t_old + delta_t_older) / 
      (delta_t_old * delta_t_older);
    c3 = (delta_t + delta_t_old + delta_t_older) * (delta_t + delta_t_old) /
      delta_t_old / (delta_t_older + delta_t_old);
      
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
    xdot[i] = (1.0 + 2.0 * theta_arg) / delta_t * (x[i] - x_old[i]) -
 	      (2.0 * theta_arg) * xdot_old[i];
  }

} /* END of routine predict_solution  */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/**************************************************************************/

static void
predict_solution_newmark(int N, double delta_t, 
			 double x[],
			 double x_old[],  
			 double xdot[], 
			 double xdot_old[])

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
  c2 = delta_t*delta_t*(1.-2.*tran->newmark_beta)/2.;
  c3 = delta_t*(1.-tran->newmark_gamma);
  for (i = 0; i < N; i++) {

    node_ptr = Nodes[i];
    curr_mat_list = &(node_ptr->Mat_List);
    num_mat = curr_mat_list->Length;

    for (imat = 0; imat < num_mat; imat++) 
      {
	mat_index= (curr_mat_list->List)[imat];
	if(pd_glob[mat_index]->MeshMotion == DYNAMIC_LAGRANGIAN) 
	  {
	    for (a = 0; a < ei->ielem_dim; a++)
	      {
		j = Index_Solution(i, R_MESH1 + a, 0, 0 , -1);
		x[j] = x_old[j] + c1 * xdot_old[j]
		  + c2 * tran->xdbl_dot_old[j];
		xdot[j] = xdot_old[j] + c3*tran->xdbl_dot_old[j];
	      }
	  }
      }
  }
} /* END of routine predict_solution_newmark  */



static int discard_previous_time_step(int num_unks, 
				      double *x,
				      double *x_old,
				      double *x_older,
				      double *x_oldest,
				      double *xdot,
				      double *xdot_old, 
				      double *xdot_older)
{
  
  dcopy1(num_unks, x,       x_old);
  dcopy1(num_unks, x_old,   x_older);
  dcopy1(num_unks, x_older, x_oldest);
	
  /* also need to kill xdot(s) */
  memset(xdot,0, sizeof(double)*num_unks);
  memset(xdot_old,0, sizeof(double)*num_unks);
  memset(xdot_older,0, sizeof(double)*num_unks);
	
  return (0);
}	

/*****************************************************************************/
/*****************************************************************************/
/****************************************************************************/

int
anneal_mesh(double x[], int tev, int tev_post, double *glob_vars_val,
	    struct Results_Description  *rd,
	    double time_value,  Exo_DB *exo, Dpi *dpi)

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
  int displacement_somewhere;	/* boolean */
  int i;
  int j;
  int m;			/* matl index */
  int num_nodes, rd_nnv_save, rd_nev_save;
  int p;
  int ielem;
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
  double	***gvec_elem;
  double *x_file;
  double d[DIM][MDE];

  double *savex=NULL, *savey=NULL, *savez=NULL; /* temporary placeholders while
						 * anneal coordinates are 
						 * written...                */
  double **new_coord;

  double *nodal_result_vector;	/* temporarily hold one nodal variable 
				 * prior to writing out into EXODUS II file */
  char afilename[MAX_FNL];

  FILE *anneal_dat;


  int numProcUnknowns = NumUnknowns + NumExtUnknowns;

  asdv(&x_file, numProcUnknowns);

  strcpy(afilename, anneal_file);

  dcopy1( numProcUnknowns, x, x_file );

  /*
   * Return immediately if there are no mesh displacements.
   */

  displacement_somewhere = FALSE;

  for(m = 0; m < upd->Num_Mat; m++)
      displacement_somewhere |= ( pd_glob[m]->e[R_MESH1] );

  if ( !displacement_somewhere ) 
  {
    WH(-1, "Attempt to anneal w/ no active displacement eqn anywhere!");
    return(0);
  }

  dim       = exo->num_dim;
  num_nodes = exo->num_nodes;
	      
  new_coord = (double **) calloc(dim, sizeof(double *));

  for(p = 0; p < dim; p++)
  {
    new_coord[p] = (double *) calloc(num_nodes, sizeof(double));
    dcopy1( num_nodes, Coor[p], new_coord[p] );
  }

  nodal_result_vector = (double *) calloc(num_nodes, sizeof(double));

  moved = (int *) calloc( num_nodes, sizeof(int) );

  memset(moved, 0, sizeof(int)*num_nodes);

  /*
   * Loop through nodes, find displacement, and add it into 
   * the coordinate
   */

  e_start = exo->eb_ptr[0];
  e_end   = exo->eb_ptr[exo->num_elem_blocks];

  for(ielem = e_start; ielem < e_end; ielem++)
  {
    load_elem_dofptr(ielem, exo, x, x_file, x, x, x, 1);

    memset(d, 0, sizeof(double  )*DIM*MDE);
    memset(dofs, 0, sizeof(int)*DIM);
    
    for(ln = 0; ln < ei->num_local_nodes; ln++)
    {
      double factor=1.0;
      double xi[3] = {0.0, 0.0, 0.0};
	  
      find_nodal_stu(ln, ei->ielem_type, xi, xi+1, xi+2);

      gnn = exo->elem_node_list[ exo->elem_node_pntr[ielem] + ln ] ;

      memset(displacement, 0, sizeof(double)*DIM);
	  
      if( moved[gnn] != 1 )
      {
	for(p = 0; p < DIM; p++)
	{
	  var = MESH_DISPLACEMENT1 + p;
		  
	  for(i = 0; i < ei->dof[var]; i++)
	  {
	    phi[i] = newshape(xi, 
			      ei->ielem_type, 
			      PSI, 
			      ei->dof_list[var][i], 
			      ei->ielem_shape,
			      pd->i[var],
			      i);
	  }

	  if( pd->v[var] )
	  {
	    for(j = 0; j < ei->dof[var]; j++)
	    {
	      displacement[p] += *esp->d[p][j] * phi[j];
	      *esp_old->d[p][j] = 0.0;
	    }
	    moved[gnn] = 1;
	  }
	  else
	      displacement[p] = 0.0;
	}
      }

      for(p = 0; p < dim; p++)
	  new_coord[p][gnn] += factor*displacement[p];
    }
  }
  
	  
  /*
   * Temporarily save the old coordinates and assign the new
   * displaced coordinates for the annealed file...
   */
  
  if( dim > 0 )
  {
    savex        = exo->x_coord;
    exo->x_coord = &new_coord[0][0];
  }

  if( dim > 1 )
  {
    savey        = exo->y_coord;
    exo->y_coord = &new_coord[1][0];
  }

  if( dim > 2 )
  {
    savez = exo->z_coord;
    exo->z_coord = &new_coord[2][0];
  }

  if( Num_Proc > 1 )
      multiname(afilename, ProcID, Num_Proc);

  one_base(exo);		/* cause node numbers, etc. to start at 1 */

  wr_mesh_exo(exo, afilename, 0);

  /*
   * Grab back the undisplaced coordinates that were originally read.
   */

  if( dim > 0 )
      exo->x_coord = savex;
  
  if( dim > 1 )
      exo->y_coord = savey;
  
  if( dim > 2 )
      exo->z_coord = savez;

  if( (tev + tev_post ) != 0 )
  {
    gvec_elem = (double ***)
	array_alloc(2, exo->num_elem_blocks, tev + tev_post, sizeof(double *));
  }
  else
  {
    gvec_elem = NULL;
  }

  /*
   * Write the names of the results variables to the annealed file.
   */
  rd_nnv_save = rd->nnv;
  rd_nev_save = rd->nev;
  rd->nnv = rd->TotalNVSolnOutput;
  rd->nev = tev;
  
  wr_result_prelim_exo(rd, 
		       exo, 
		       afilename,
		       gvec_elem );

  rd->nnv = rd_nnv_save; /* Return values to normal */
  rd->nev = rd_nev_save;

  if( Num_Proc > 1 )
      wr_dpi(dpi, afilename, 0);
	      
  for (i = 0; i < rd->TotalNVSolnOutput; i++) {
    /* 
     * if nodal variable is a displacement, set it to zero 
     */
    if(rd->nvtype[i] == MESH_DISPLACEMENT1 || 
       rd->nvtype[i] == MESH_DISPLACEMENT2 || 
       rd->nvtype[i] == MESH_DISPLACEMENT3) {
      init_vec_value(nodal_result_vector, 0.0, num_nodes);
    } else {
      extract_nodal_vec(x, rd->nvtype[i], rd->nvkind[i], 
			rd->nvmatID[i],
			nodal_result_vector, exo, FALSE, time_value);
    }
    wr_nodal_result_exo(exo, afilename, nodal_result_vector, i+1, 
			1, time_value);
  }

  /* Now pick up the element variables */
  for (i = 0; i < tev; i++) {
    wr_elem_result_exo(exo, anneal_file, gvec_elem, i, 1,
		       time_value, rd);
  }

  /* And finally the global variables */

  if( rd->ngv > 0  )
    {
      wr_global_result_exo( exo, anneal_file, 1, rd->ngv, glob_vars_val );
    }

  anneal_dat = fopen("anneal.dat", "w");

  (void) write_ascii_soln(x_file, NULL, numProcUnknowns,
		          x, 0, 0.0, anneal_dat);

  fclose(anneal_dat);

  /*
   * Return internal EXODUS II data to 0 based node and element numbers
   * for convenience.
   */

  zero_base(exo);
	  
  /*
   * Free up any temporarily allocated memory.
   */

  for(p = 0; p < dim; p++) {
    safer_free((void **) &(new_coord[p]));
  }
  safer_free((void **) &new_coord);
  safer_free((void **) &nodal_result_vector);
  safer_free((void **) &x_file);
  safer_free((void **) &moved);
  safer_free((void **) &gvec_elem);

  return(0);
} /* END of routine anneal_mesh */
/*****************************************************************************/

/*****************************************************************************/
/*****************************************************************************/

int
anneal_mesh_with_external_field(const Exo_DB *exo)

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

  dim       = exo->num_dim;
  num_nodes = exo->num_nodes;
	      
  new_coord = (double **) calloc(dim, sizeof(double *));

  for(p = 0; p < dim; p++)
  {
    new_coord[p] = (double *) calloc(num_nodes, sizeof(double));
    dcopy1( num_nodes, Coor[p], new_coord[p] );
  }

  for (i=0; i<num_nodes; i++)
    {
      /*
       * Default is full anneal: X_new = X_old + displacement
       */ 

      X_old[0]        = exo->x_coord[i];
      X_old[1]        = exo->y_coord[i];
      if( dim > 2 )
	{
	  X_old[2]        = exo->z_coord[i];
	}
      for (p=0; p<dim; p++)
	{

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
	      
      for ( p=0; p<dim; p++)
	{
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
  
    Coor[0]      = &new_coord[0][0];   /*Old Deprecated holder that is stil
					*used
					*/
    Coor[1]      = &new_coord[1][0];

  if( dim > 2 )
  {
    Coor[2]      = &new_coord[2][0];
  }

  for(p = 0; p < dim; p++) {
    safer_free((void **) &(new_coord[p]));
  }
  safer_free((void **) &new_coord);

  return(0);
} /* END of routine anneal_mesh_with_external_field */


void
shift_nodal_values ( int var,
					 double shift,
					 double *x, 
					 int num_total_nodes )
{
	int I, ie, matID=-1;
	
	for( I=0; I<num_total_nodes; I++ )
	{
		ie = Index_Solution(I, var, 0, 0, matID );
		if( ie != -1 )
		{
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
int
load_export_vars(const int nodes,
                 dbl x[],
                 dbl *x_pp)
                                                                                
{
  double var;
  int i, j, k, ivar;
  int err=0;
  int scount=0, pcount=0;
/*
  int nsoln[MAX_EXTERNAL_FIELD];
  int npost[MAX_EXTERNAL_FIELD];
*/
                                                                                
                                                                                
  /* Ensure that exports were requested */
  if (Num_Export_XS == 0 && Num_Export_XP == 0)
    {
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
  if (libio->xsoln == NULL || libio->xpost == NULL)
    {
      fprintf(stderr, "Invalid variable storage array: xsoln or xpost");
      return -1;
    }
                                                                                
  /* Retrieve solution variables in the order requested */
  for (i=0; i<Num_Export_XS; i++)
    {
      ivar = Export_XS_ID[i];
      if (Num_Var_In_Type[ivar] == 0)
        {
          fprintf(stderr, "Inactive variable specified: ID = %d\n", ivar);
          return -1;
        }
      for (j=0; j<nodes; j++)
        {
          k = Index_Solution(j, ivar, 0, 0, -1);
          var = ( (k>=0) ? x[k] : 0.0);
          libio->xsoln[scount] = var;
          scount++;
        }
    }
                                                                                
  /* Retrieve post-processing variables in the order requested */
  if (Num_Export_XP > 0 && x_pp != NULL)
    {
      for (i=0; i<Num_Export_XP; i++)
        {
          for (j=0; j<nodes; j++)
            {
              libio->xpost[pcount] = x_pp[pcount];
              pcount++;
            }
        }
    }
                                                                                
  return 0;
}
#endif

/*****************************************************************************/
/*  END of file rf_solve.c  */
/*****************************************************************************/
