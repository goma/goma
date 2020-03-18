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

/* directs and controls solution process for hunting with zero and
 * first order  hunting
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define GOMA_AC_HUNT_C
#include "brk_utils.h"
#include "ac_hunt.h"
#include "ac_update_parameter.h"
#include "az_aztec.h"
#include "dp_comm.h"
#include "dp_types.h"
#include "dp_utils.h"
#include "dpi.h"
#include "el_geom.h"
#include "exo_struct.h"
#include "mm_as.h"
#include "mm_as_structs.h"
#include "mm_augc_util.h"
#include "mm_bc.h"
#include "mm_eh.h"
#include "mm_fill_util.h"
#include "mm_flux.h"
#include "mm_more_utils.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "mm_post_def.h"
#include "mm_post_proc.h"
#include "mm_sol_nonlinear.h"
#include "mpi.h"
#include "rf_allo.h"
#include "rf_bc.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_io.h"
#include "rf_io_structs.h"
#include "rf_mp.h"
#include "rf_node_const.h"
#include "rf_solve.h"
#include "rf_solver.h"
#include "rf_solver_const.h"
#include "rf_util.h"
#include "sl_auxutil.h"
#include "sl_util_structs.h"
#include "std.h"
#include "wr_exo.h"
#include "wr_soln.h"

#ifdef HAVE_FRONT
extern int mf_setup
(int *,			/* nelem_glob */
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
       int *);			/* allocated */
#endif

int w;

#include "sl_util.h"		/* defines sl_init() */
#include "el_quality.h"

/*

   HUNTING WITH ZEROTH AND FIRST ORDER CONTINUATION
   WITH AUGMENTING CONDITIONS

   BASED ON solve_problem() IN rf_solve.c

   BY IAN GATES & ROBERT SECOR

   2/98 - 11/99

*/

void
hunt_problem(Comm_Ex *cx,	/* array of communications structures */
	     Exo_DB *exo,	/* ptr to the finite element mesh database */
	     Dpi *dpi)	        /* distributed processing information */
{
  int    *ija=NULL;           /* column pointer array                         */
  double *a=NULL;             /* nonzero array                                */
  double *a_old=NULL;         /* nonzero array                                */
  double *x=NULL;             /* solution vector                              */

  int     iAC;                /* COUNTER                                      */
  double *x_AC = NULL;        /* SOLUTION VECTOR OF EXTRA UNKNOWNS            */
  double *x_AC_old=NULL;      /* old SOLUTION VECTOR OF EXTRA UNKNOWNS        */
  double *x_AC_dot = NULL;

  int     iHC;                /* COUNTER                                      */

  int    *ija_attic=NULL;     /* storage for external dofs                    */

  int eb_indx, ev_indx;

  /*
   * variables for path traversal
   */

  double *x_old=NULL;         /* old solution vector                          */
  double *x_older=NULL;       /* older solution vector                        */
  double *x_oldest=NULL;      /* oldest solution vector saved                 */
  double *xdot=NULL;          /* current path derivative of soln              */
  double *xdot_old=NULL;
  double *x_update=NULL;

  double *x_sens=NULL;        /* solution sensitivity */
  double **x_sens_p=NULL;     /* solution sensitivity for parameters */
  int num_pvector=0;          /*  number of solution sensitivity vectors */
#ifdef COUPLED_FILL
  struct Aztec_Linear_Solver_System *ams[NUM_ALSS]={NULL};
#else /* COUPLED_FILL */
  struct Aztec_Linear_Solver_System *ams[NUM_ALSS]={NULL, NULL};
#endif /* COUPLED_FILL */
                              /* sl_util_structs.h */

  double *resid_vector=NULL;  /* residual */
  double *resid_vector_sens=NULL;    /* residual sensitivity */
  double *scale=NULL;      /* scale vector for modified newton */

  int 	 *node_to_fill = NULL;

  int		n;            /* total number of path steps attempted */
  int		ni;           /* total number of nonlinear solves */
  int		nt;           /* total number of successful path steps */
  int		path_step_reform; /* counter for jacobian reformation stride */
  int		converged;    /* success or failure of Newton iteration */
  int		success_ds;   /* success or failure of path step */

  int           i;

  int           nprint=0, num_total_nodes;

  int           numProcUnknowns;
  int           *const_delta_s=NULL;
  int           step_print;
  double        i_print;
  int           step_fix = 0;      /* What step to fix the problem on */
  int           good_mesh = TRUE;
  double	*path=NULL, *path1=NULL;
  double	*delta_s=NULL, *delta_s_new=NULL, *delta_s_old=NULL;
  double        *delta_s_older=NULL, *delta_s_oldest=NULL;
  double        *hDelta_s0=NULL, *hDelta_s_min=NULL, *hDelta_s_max=NULL;
  double        delta_t;
  double	theta=0.0;
  double        eps;
  double        *lambda=NULL, *lambdaEnd=NULL;
  double	hunt_par, dhunt_par, hunt_par_old;	/* hunting continuation parameter */
  double        dhunt_par_max=1.0, dhunt_par_min=0., dhunt_par_0=0.1;
  double        dhunt_par_new=0.1, dhunt_par_old;
  int           log_ID=-1;
  double        timeValueRead = 0.0;

  /*
   * ALC management variables
   */

  int           alqALC;
  int           *aldALC=NULL;

  /*
   * Other local variables
   */

  int	        error, err, is_steady_state, inewton;
  int 		*gindex = NULL, gsize;
  int		*p_gsize=NULL;
  double	*gvec=NULL;
  double        ***gvec_elem;
  FILE          *file=NULL;
  double 	toler_org[3];
  double        exo_time=0;

  struct Results_Description  *rd=NULL;

  int		tnv;		/* total number of nodal variables and kinds */
  int		tev;		/* total number of elem variables and kinds */
  int		tnv_post;	/* total number of nodal variables and kinds
					   for post processing */
  int		tev_post;	/* total number of elem variables and kinds
					   for post processing */
  double        *gv;

#ifdef HAVE_FRONT
  int max_unk_elem, one, three; /* variables used as mf_setup arguments*/
#endif

  unsigned int
  matrix_systems_mask;

  double evol_local=0.0;
#ifdef PARALLEL
  double evol_global=0.0;
#endif

  /* Set step_fix only if parallel run and only if fix freq is enabled*/
  if (Num_Proc > 1 && cont->fix_freq > 0) {
    step_fix = 1; /* Always fix on the first timestep to match print frequency */
  }


  static char yo[]="hunt_problem";

  /*
   * 		BEGIN EXECUTION
   */

#ifdef DEBUG
  fprintf(stderr, "hunt_problem() begins...\n");
#endif

  toler_org[0] = custom_tol1;
  toler_org[1] = custom_tol2;
  toler_org[2] = custom_tol3;

  is_steady_state = TRUE;

  p_gsize = &gsize;

  /*
   * set aside space for gather global vectors to print to exoII file
   * note: this is temporary
   *
   * For 2D prototype problem:  allocate space for T, dx, dy arrays
   */

  if( strlen( Soln_OutFile)  )
    {
#ifdef DEBUG
      printf("Trying to open \"%s\" for writing.\n", Soln_OutFile);
#endif
      file = fopen(Soln_OutFile, "w");
      if (file == NULL)  {
        DPRINTF(stdout, "%s:  opening soln file for writing\n", yo);
        EH(-1, "\t");
      }
    }
#ifdef PARALLEL
  check_parallel_error("Soln output file error");
#endif

  /*
   * Some preliminaries to help setup EXODUS II database output.
   */

#ifdef DEBUG
  fprintf(stderr, "cnt_nodal_vars() begins...\n");
#endif

  tnv = cnt_nodal_vars();
  /*  tnv_post is calculated in load_nodal_tkn*/
  tev = cnt_elem_vars();
  /*  tev_post is calculated in load_elem_tkn*/

#ifdef DEBUG
  fprintf(stderr, "Found %d total primitive nodal variables to output.\n", tnv);
  fprintf(stderr, "Found %d total primitive elem variables to output.\n", tev);
#endif

  if ( tnv < 0 )
    {
      DPRINTF(stderr, "%s:\tbad tnv.\n", yo);
      EH(-1, "\t");
    }

  if ( tev < 0 )
    {
      DPRINTF(stderr, "%s:\tMaybe bad tev? See goma design committee ;) \n", yo);
/*       exit(-1); */
    }

  rd = (struct Results_Description *)
    smalloc(sizeof(struct Results_Description));

  if (rd == NULL)
    { EH(-1, "Could not grab Results Description."); }
  (void) memset((void *) rd, 0, sizeof(struct Results_Description));

  rd->nev = 0;			/* number element variables in results */
  rd->ngv = 5 + nAC;	        /* number global variables in results
				   see load_global_var_info for names*/
  rd->nhv = 0;			/* number history variables in results */

  if ( is_steady_state == TRUE ) {
    error = load_global_var_info(rd, 0, "CONV");
    error = load_global_var_info(rd, 1, "NEWT_IT");
    error = load_global_var_info(rd, 2, "MAX_IT");
    error = load_global_var_info(rd, 3, "CONVRATE");
    error = load_global_var_info(rd, 4, "MESH_VOLUME");

  }

    if ( nAC > 0   )
    {
      char name[20];

      for( i = 0 ; i < nAC ; i++ )
	{
	  sprintf(name, "AUGC_%d",i+1);
	  error = load_global_var_info(rd, 5 + i, name);
	}
    }

  gv = alloc_dbl_1( rd->ngv, 0.0 );

  /* load nodal types, kinds, names */
  error = load_nodal_tkn( rd,
			  &tnv,
			  &tnv_post); /* load nodal types, kinds, names */

  if (error !=0)
    {
      DPRINTF(stderr, "%s:  problem with load_nodal_tkn()\n", yo);
      EH(-1,"\t");
    }

  /* load elem types, names */
  error = load_elem_tkn( rd,
			 exo,
			 tev,
			 &tev_post); /* load elem types, names */

  if ( error !=0 )
    {
      DPRINTF(stderr, "%s:  problem with load_elem_tkn()\n", yo);
      EH(-1,"\t");
    }

  /*
   * Write out the names of the nodal variables that we will be sending to
   * the EXODUS II output file later.
   */

#ifdef DEBUG
  fprintf(stderr, "wr_result_prelim() starts...\n", tnv);
#endif

  gvec_elem = (double ***) smalloc ( (exo->num_elem_blocks)*sizeof(double **));
  for (i = 0; i < exo->num_elem_blocks; i++) {
    gvec_elem[i] = (double **) smalloc ( (tev + tev_post)*sizeof(double *));
  }

  wr_result_prelim_exo( rd,
                        exo,
                        ExoFileOut,
                        gvec_elem );

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

  num_total_nodes = dpi->num_universe_nodes;

  numProcUnknowns = NumUnknowns[pg->imtrx] + NumExtUnknowns[pg->imtrx];

  /* allocate memory for Volume Constraint Jacobian. ACS 2/99 */

  if ( nAC > 0)
    {
      for(iAC=0;iAC<nAC;iAC++) {
	augc[iAC].d_evol_dx = (double*) malloc(numProcUnknowns*sizeof(double));
      } }

  asdv(&resid_vector, numProcUnknowns);
  asdv(&resid_vector_sens, numProcUnknowns);
  asdv(&scale, numProcUnknowns);

  for (i=0;i<NUM_ALSS;i++)
    {
      ams[i] = (struct Aztec_Linear_Solver_System *)
	array_alloc(1, 1, sizeof(struct Aztec_Linear_Solver_System ));
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

  /*
   * allocate space for and initialize solution arrays
   */

  asdv(&x,        numProcUnknowns);
  asdv(&x_old,    numProcUnknowns);
  asdv(&x_older,  numProcUnknowns);
  asdv(&x_oldest, numProcUnknowns);
  asdv(&xdot,     numProcUnknowns);
  asdv(&xdot_old, numProcUnknowns);
  asdv(&x_update, numProcUnknowns);

  asdv(&x_sens, numProcUnknowns);

  /*
   * Initialize solid inertia flag
   */
  set_solid_inertia();

  /*
   * ALLOCATE ALL THOSE WORK VECTORS FOR HUNTING
   */

  asdv(&lambda,         nHC);
  asdv(&lambdaEnd,      nHC);
  asdv(&path,           nHC);
  asdv(&path1,          nHC);
  asdv(&hDelta_s0,      nHC);
  asdv(&hDelta_s_min,   nHC);
  asdv(&hDelta_s_max,   nHC);
  asdv(&delta_s,        nHC);
  asdv(&delta_s_new,    nHC);
  asdv(&delta_s_old,    nHC);
  asdv(&delta_s_older,  nHC);
  asdv(&delta_s_oldest, nHC);

  aldALC        = Ivector_birth(nHC);
  const_delta_s = Ivector_birth(nHC);

  if (nAC > 0)
  {
    asdv(&x_AC, nAC);
    asdv(&x_AC_old, nAC);
    asdv(&x_AC_dot, nAC);
  }


  /*

   HUNTING BY ZERO AND FIRST ORDER CONTINUATION

  */

  alqALC = 1;

  delta_t = 0.0;
  tran->delta_t = 0.0;      /*for Newmark-Beta terms in Lagrangian Solid*/

  nprint = 0;

  MaxPathSteps      = cont->MaxPathSteps;
  eps               = cont->eps;

  for (iHC=0;iHC<nHC;iHC++) {

    const_delta_s[iHC] = 0;

    lambda[iHC]       = hunt[iHC].BegParameterValue;
    lambdaEnd[iHC]    = hunt[iHC].EndParameterValue;
    if(hunt[iHC].ramp == 2 )
         { if(log_ID == -1) log_ID = iHC; }

    if ((lambdaEnd[iHC]-lambda[iHC]) > 0.0)
         { aldALC[iHC] = +1; }
    else
         { aldALC[iHC] = -1; }

    if (hunt[iHC].ramp == 1) {
      hunt[iHC].Delta_s0 = fabs(lambdaEnd[iHC]-lambda[iHC])/((double)(MaxPathSteps-1));
      const_delta_s[iHC] = 1;
    }

    hDelta_s0[iHC]     = hunt[iHC].Delta_s0;
    hDelta_s_min[iHC]  = hunt[iHC].Delta_s_min;
    hDelta_s_max[iHC]  = hunt[iHC].Delta_s_max;

    path[iHC] = path1[iHC] = lambda[iHC];

    if (Debug_Flag && ProcID == 0) {
      fprintf(stdout,"MaxPathSteps: %d \tlambdaEnd: %f\n", MaxPathSteps, lambdaEnd[iHC]);
      fprintf(stdout,"continuation in progress\n");
    }

    if (hDelta_s0[iHC] > hDelta_s_max[iHC])
    {
      hDelta_s0[iHC] = hDelta_s_max[iHC];
    }

    delta_s[iHC] = delta_s_old[iHC] = delta_s_older[iHC] = hDelta_s0[iHC];

    /*
     * ADJUST NATURAL PARAMETER
     */

    update_parameterHC(iHC, path1[iHC], x, xdot, x_AC, delta_s[iHC], cx, exo, dpi);
  }

  /*  define continuation parameter */

  iHC = MAX(0,log_ID);
  dhunt_par = 0.0;
  if(hunt[iHC].EndParameterValue == hunt[iHC].BegParameterValue)
 	{	hunt_par = 1.0;	}
  else
 	{
          if(hunt[iHC].ramp == 2 )
              {
	         hunt_par = (log10(path1[iHC])-log10(hunt[iHC].BegParameterValue))
	             /(log10(hunt[iHC].EndParameterValue) - log10(hunt[iHC].BegParameterValue));
	         dhunt_par_min = log10(1.0+aldALC[iHC]*hunt[iHC].Delta_s_min/hunt[iHC].BegParameterValue)
	             /(log10(hunt[iHC].EndParameterValue) - log10(hunt[iHC].BegParameterValue));
	         dhunt_par_max = log10(1.0+aldALC[iHC]*hunt[iHC].Delta_s_max/hunt[iHC].BegParameterValue)
	             /(log10(hunt[iHC].EndParameterValue) - log10(hunt[iHC].BegParameterValue));
	         dhunt_par_0 = log10(1.0+aldALC[iHC]*hunt[iHC].Delta_s0/hunt[iHC].BegParameterValue)
	             /(log10(hunt[iHC].EndParameterValue) - log10(hunt[iHC].BegParameterValue));
              }  else   {
	         hunt_par = (path1[iHC]-hunt[iHC].BegParameterValue)
	             /(hunt[iHC].EndParameterValue - hunt[iHC].BegParameterValue);
	         dhunt_par_min = aldALC[iHC]*hunt[iHC].Delta_s_min
	             /(hunt[iHC].EndParameterValue - hunt[iHC].BegParameterValue);
	         dhunt_par_max = aldALC[iHC]*hunt[iHC].Delta_s_max
	             /(hunt[iHC].EndParameterValue - hunt[iHC].BegParameterValue);
	         dhunt_par_0 = aldALC[iHC]*hunt[iHC].Delta_s0
	             /(hunt[iHC].EndParameterValue - hunt[iHC].BegParameterValue);
              }
           hunt_par=fabs(hunt_par);
 	}
  hunt_par_old = hunt_par;
  dhunt_par = dhunt_par_old = dhunt_par_0;
    if (dhunt_par_0 > dhunt_par_max)
    { dhunt_par_0 = dhunt_par_max; }

  /* Call prefront (or mf_setup) if necessary */
  if (Linear_Solver == FRONT)
  {
    if (Num_Proc > 1) EH(-1, "Whoa.  No front allowed with nproc>1");

#ifdef HAVE_FRONT
    /* Also got to define these because it wants pointers to these numbers */

    max_unk_elem = (MAX_PROB_VAR + MAX_CONC)*MDE;
    one = 1;
    three = 3;

    /* NOTE: We need a overall flag in the vn_glob struct that tells whether FULL_DG
       is on anywhere in domain.  This assumes only one material.  See sl_front_setup for test.
       that test needs to be in the input parser.  */

    if(vn_glob[0]->dg_J_model == FULL_DG)
    {
      max_unk_elem = (MAX_PROB_VAR + MAX_CONC)*MDE + 4*vn_glob[0]->modes*4*MDE;
    }

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
         *  if compute parameter sensitivities, allocate space for solution
         *  sensitivity vectors
         */

        for(i=0;i<nn_post_fluxes_sens;i++)      {
          num_pvector=MAX(num_pvector,pp_fluxes_sens[i]->vector_id);}
        for(i=0;i<nn_post_data_sens;i++)        {
          num_pvector=MAX(num_pvector,pp_data_sens[i]->vector_id);}

  if((nn_post_fluxes_sens + nn_post_data_sens) > 0)
  {
    num_pvector++;
    num_pvector = MAX(num_pvector,2);
        x_sens_p = Dmatrix_birth(num_pvector,numProcUnknowns);
  }
  else
  {
    x_sens_p = NULL;
  }

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

  if( strcmp( Matrix_Format, "msr" ) == 0)
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
    ams[JAC]->npn_plus = dpi->num_internal_nodes + dpi->num_boundary_nodes + dpi->num_external_nodes;

    ams[JAC]->npu      = num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx];
    ams[JAC]->npu_plus = num_universe_dofs[pg->imtrx];

    ams[JAC]->nnz = ija[num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx]] - 1;
    ams[JAC]->nnz_plus = ija[num_universe_dofs[pg->imtrx]];

  }
  else if(  strcmp( Matrix_Format, "vbr" ) == 0)
  {
    log_msg("alloc_VBR_sparse_arrays...");
    alloc_VBR_sparse_arrays ( ams[JAC],
			      exo,
			      dpi);
    ija_attic = NULL;
    ams[JAC]->belfry  = ija_attic;

    a = ams[JAC]->val;
    if( !save_old_A ) a_old = ams[JAC]->val_old;
  }
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
  {
    EH(-1,"Attempted to allocate unknown sparse matrix format");
  }

  init_vec(x, cx, exo, dpi, x_AC, nAC, &timeValueRead);

/*  if read ACs, update data floats */
  if (nAC > 0)
  {
    if(augc[0].iread == 1)
      {
	for(iAC=0 ; iAC<nAC ; iAC++)
	  { update_parameterAC(iAC, x, xdot, x_AC, cx, exo, dpi); }
      }
  }


  /*
       * set boundary conditions on the initial conditions
       */

  find_and_set_Dirichlet(x, xdot, exo, dpi);

  exchange_dof(cx, dpi, x, pg->imtrx);

  dcopy1(numProcUnknowns,x,x_old);
  dcopy1(numProcUnknowns,x_old,x_older);
  dcopy1(numProcUnknowns,x_older,x_oldest);

  if( nAC > 0)
  {
    dcopy1(nAC,x_AC, x_AC_old);}

  /*
       * initialize the counters for when to print out data
       */

  step_print = 1;

  matrix_systems_mask = 1;

  log_msg("sl_init()...");
  sl_init(matrix_systems_mask, ams, exo, dpi, cx);

#ifdef PARALLEL
  /*
  * Make sure the solver was properly initialized on all processors.
  */
  check_parallel_error("Solver initialization problems");
#endif

      ams[JAC]->options[AZ_keep_info] = 1;

    DPRINTF(stdout, "\nINITIAL ELEMENT QUALITY CHECK---\n");
    good_mesh = element_quality(exo, x, ams[0]->proc_config);

  /*
       * set the number of successful path steps to zero
       */

  nt = 0;

  /*
       * LOOP THROUGH PARAMETER UNTIL MAX NUMBER
       * OF STEPS SURPASSED
       */

  if (nAC > 0) {
    dcopy1( nAC, x_AC, &(gv[5]) );
  }

  for (n=0;n<MaxPathSteps;n++) {

    alqALC = 1;

    for (iHC=0;iHC<nHC;iHC++) {

      switch (aldALC[iHC]) {
      case -1: /* REDUCING PARAMETER DIRECTION */
	  if (path1[iHC] <= lambdaEnd[iHC]) {
	    alqALC = -1;
	    path1[iHC] = lambdaEnd[iHC];
	    delta_s[iHC] = path[iHC]-path1[iHC];
	  }
	  break;
      case +1: /* RISING PARAMETER DIRECTION */
	  if (path1[iHC] >= lambdaEnd[iHC]) {
	    alqALC = -1;
	    path1[iHC] = lambdaEnd[iHC];
	    delta_s[iHC] = path1[iHC]-path[iHC];
	  }
	  break;
      }
    }   /*  end of iHC loop */

      /*
       * ADJUST NATURAL PARAMETER
       */

    for (iHC=0;iHC<nHC;iHC++) {
      update_parameterHC(iHC, path1[iHC], x, xdot, x_AC, delta_s[iHC], cx, exo, dpi);
    }   /*  end of iHC loop */

        iHC = MAX(0,log_ID);
  	if(hunt[iHC].EndParameterValue == hunt[iHC].BegParameterValue)
 		{	hunt_par = 1.0;	}
	else
 		{
                  if(hunt[iHC].ramp == 2 )
                       {
	                 hunt_par = (log10(path1[iHC])-log10(hunt[iHC].BegParameterValue))
	                     /(log10(hunt[iHC].EndParameterValue) - log10(hunt[iHC].BegParameterValue));
                         if(hunt_par > 0)	{
	                        dhunt_par = (log10(path1[iHC])-log10(path[iHC]))
	      /(log10(hunt[iHC].EndParameterValue) - log10(hunt[iHC].BegParameterValue));
                                }
                       }  else   {
		         hunt_par = (path1[iHC]-hunt[iHC].BegParameterValue)
		             /(hunt[iHC].EndParameterValue - hunt[iHC].BegParameterValue)  ;
                         if(hunt_par > 0)	{
		                 dhunt_par = (path1[iHC]-path[iHC])
		                       /(hunt[iHC].EndParameterValue - hunt[iHC].BegParameterValue)  ;
                                }
                       }
/*                  hunt_par=fabs(hunt_par);*/
 		}
	  for (iHC=0;iHC<nHC;iHC++) {
              if(hunt[iHC].ramp == 2)
                {
                 delta_s[iHC] = -hunt[iHC].BegParameterValue *
                     pow(hunt[iHC].EndParameterValue/hunt[iHC].BegParameterValue,hunt_par-dhunt_par);
                 delta_s[iHC] += hunt[iHC].BegParameterValue *
                     pow(hunt[iHC].EndParameterValue/hunt[iHC].BegParameterValue,hunt_par);
                 }
              else  {
                 delta_s[iHC] = dhunt_par*
           (hunt[iHC].EndParameterValue - hunt[iHC].BegParameterValue);
                 }
              }

    /*
     * IF STEP CHANGED, REDO FIRST ORDER PREDICTION
     */
 if (hunt_par >= 1.0) { alqALC = -1;  }

    if(alqALC == -1)
    {
      DPRINTF(stdout,"\n\t ******** LAST PATH STEP!\n");
      dcopy1(numProcUnknowns,x_old,x);

      switch (Continuation) {
      case HUN_ZEROTH:
          break;
      case  HUN_FIRST:
          v1add(numProcUnknowns, &x[0], dhunt_par, &x_sens[0]);
	  break;
      }
    }

    /*
     * reset Dirichlet condition Mask, node->DBC to -1 where it
     * is set in order for Dirichlet conditions to be
     * set appropriately for each path step
     */

    nullify_dirichlet_bcs();

    find_and_set_Dirichlet (x, xdot, exo, dpi);

    exchange_dof(cx, dpi, x, pg->imtrx);

    if(ProcID ==0) {
      DPRINTF(stdout, "\n\t----------------------------------");
      switch (Continuation) {
      case HUN_ZEROTH:
          DPRINTF(stdout, "\n\tZero Order Hunting:");
	  break;
      case  HUN_FIRST:
          DPRINTF(stdout, "\n\tFirst Order Hunting:");
	  break; }
      DPRINTF(stdout, "\n\tStep number: %4d of %4d (max)", n+1, MaxPathSteps);
      DPRINTF(stdout, "\n\tAttempting solution at: theta = %g ;  step = %g",hunt_par,dhunt_par);
      for (iHC=0;iHC<nHC;iHC++) {
	switch (hunt[iHC].Type) {
	case 1: /* BC */
            DPRINTF(stdout, "\n\tBCID=%3d DFID=%5d", hunt[iHC].BCID, hunt[iHC].DFID);
	    break;
	case 2: /* MT */
            DPRINTF(stdout, "\n\tMTID=%3d MPID=%5d", hunt[iHC].MTID+1, hunt[iHC].MPID);
	    break;
 	case 3: /* AC */
            DPRINTF(stdout, "\n\tACID=%3d DFID=%5d", hunt[iHC].BCID, hunt[iHC].DFID);
 	    break;
	}
        DPRINTF(stdout, " Parameter= % 10.6e delta_s= %10.6e", path1[iHC], delta_s[iHC]);
      }
    }

    ni = 0;
    do {

#ifdef DEBUG
      fprintf(stderr, "%s: starting solve_nonlinear_problem\n", yo);
#endif
      err = solve_nonlinear_problem(ams[JAC],
				    x,
				    delta_t,
				    theta,
				    x_old,
				    x_older,
				    xdot,
				    xdot_old,
				    resid_vector,
				    x_update,
				    scale,
				    &converged,
				    &nprint,
				    tev,
				    tev_post,
				    gv,
				    rd,
				    gindex,
				    p_gsize,
				    gvec,
				    gvec_elem,
 				    path1[0],
				    exo,
				    dpi,
				    cx,
				    0,
				    &path_step_reform,
				    is_steady_state,
				    x_AC,
 				    x_AC_dot,
				    hunt_par,
				    resid_vector_sens,
				    x_sens,
				    x_sens_p,
                                    NULL);

#ifdef DEBUG
      fprintf(stderr, "%s: returned from solve_nonlinear_problem\n", yo);
#endif

      if (err == -1) converged = 0;
      inewton = err;
      if (converged)
      {
	EH(error, "error writing ASCII soln file."); /* srs need to check */

	if (Write_Intermediate_Solutions == 0) {
#ifdef DEBUG
	  fprintf(stderr, "%s: write_solution call WIS\n", yo);
#endif
          exo_time = aldALC[0]*path1[0];
	  write_solution(ExoFileOut, resid_vector, x, x_sens_p, x_old,
			 xdot, xdot_old, tev, tev_post, gv,  rd,
			 gvec, gvec_elem, &nprint, delta_s[0],
                         theta, exo_time, NULL, exo, dpi);
#ifdef DEBUG
	  fprintf(stderr, "%s: write_solution end call WIS\n", yo);
#endif
	}

	/*
	 * PRINT OUT VALUES OF EXTRA UNKNOWNS
	 * FROM AUGMENTING CONDITIONS
	 */

	if (nAC > 0)
          {

            DPRINTF(stdout, "\n------------------------------\n");
            DPRINTF(stdout, "Augmenting Conditions:    %4d\n", nAC);
            DPRINTF(stdout, "Number of extra unknowns: %4d\n\n", nAC);

            for (iAC = 0; iAC < nAC; iAC++)
             {
              if (augc[iAC].Type == AC_USERBC)
               {
                DPRINTF(stdout, "\tAC[%4d] DF[%4d] = %10.6e\n",
                        augc[iAC].BCID, augc[iAC].DFID, x_AC[iAC]);
               }
              else if (augc[iAC].Type == AC_USERMAT  ||
                       augc[iAC].Type == AC_FLUX_MAT )
               {
                DPRINTF(stdout, "\n MT[%4d] MP[%4d] = %10.6e\n",
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
                DPRINTF(stdout, "\tMT[%4d] VC[%4d]=%10.6e Param=%10.6e\n",
                        augc[iAC].MTID, augc[iAC].VOLID, evol_local,
                        x_AC[iAC]);
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
                        x_AC[iAC]);
               }
               else if(augc[iAC].Type == AC_FLUX)
               {
                DPRINTF(stdout, "\tBC[%4d] DF[%4d]=%10.6e\n",
                        augc[iAC].BCID, augc[iAC].DFID, x_AC[iAC]);
               }
             }
	  }

      /* Check element quality */
      good_mesh = element_quality(exo, x, ams[0]->proc_config);

	/*

	  INTEGRATE FLUXES, FORCES

	*/

	for (i = 0; i < nn_post_fluxes; i++)
	{
	  evaluate_flux ( exo, dpi,
			  pp_fluxes[i]->ss_id,
			  pp_fluxes[i]->flux_type ,
			  pp_fluxes[i]->flux_type_name ,
			  pp_fluxes[i]->blk_id ,
			  pp_fluxes[i]->species_number,
			  pp_fluxes[i]->flux_filenm,
			  pp_fluxes[i]->profile_flag,
			  x,xdot,NULL,delta_s[0],path1[0],1);
	}


	/*
	  COMPUTE FLUX, FORCE SENSITIVITIES
	*/


	for (i = 0; i < nn_post_fluxes_sens; i++)
	{
	  evaluate_flux_sens ( exo, dpi,
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
			       x,xdot,x_sens_p,delta_s[0],path1[0],1);
	}
 	/*
      	 * Compute global volumetric quantities
      	 */
     	 for (i = 0; i < nn_volume; i++ ) {
       		evaluate_volume_integral(exo, dpi,
                                pp_volume[i]->volume_type,
                                pp_volume[i]->volume_name,
                                pp_volume[i]->blk_id,
                                pp_volume[i]->species_no,
                                pp_volume[i]->volume_fname,
                                pp_volume[i]->params,
                                pp_volume[i]->num_params,
                                NULL,  x, xdot, delta_s[0],
                                path1[0], 1);
     		}

      } /* end of if converged block */

      /*
       * INCREMENT COUNTER
       */

      ni++;

      /*
       *
       * DID IT CONVERGE ?
       * IF NOT, REDUCE STEP SIZE AND TRY AGAIN
       *
       */

      if (!converged) {

	if (ni > Max_Newton_Steps/2) {
          DPRINTF(stdout,"\n ************************************\n");
          DPRINTF(stdout," W: Did not converge in Newton steps.\n");
          DPRINTF(stdout,"    Find better initial guess.       \n");
          DPRINTF(stdout," ************************************\n");
          goto free_and_clear;
 	  /*exit(0);  */
	}

        /*
         * ADJUST STEP SIZE - unless failed on first step
         */

        if ( nt != 0 )
        {
	DPRINTF(stderr, "\n\tFailed to converge:\n");

	dhunt_par *= 0.5;
        hunt_par = hunt_par_old + dhunt_par;
	for (iHC=0;iHC<nHC;iHC++) {
            if(hunt[iHC].ramp == 2)
              {
               path1[iHC] = hunt[iHC].BegParameterValue *
                   pow(hunt[iHC].EndParameterValue/hunt[iHC].BegParameterValue,hunt_par);
              }
            else if(hunt[iHC].ramp == 1)
              {
               delta_s[iHC] *= 0.5;
 	       switch (aldALC[iHC]) {
	           case -1:
	               path1[iHC] = path[iHC] - delta_s[iHC];
	               break;
	           case +1:
	               path1[iHC] = path[iHC] + delta_s[iHC];
	               break;
	           }
              }
            else
              {
               path1[iHC] = hunt[iHC].BegParameterValue +
                      hunt_par*(hunt[iHC].EndParameterValue - hunt[iHC].BegParameterValue);
              }
            }

	alqALC = 1;

        DPRINTF(stdout, "Decreasing step-length to %10.6e.\n", dhunt_par);

	if (dhunt_par < dhunt_par_min) {
          DPRINTF(stdout,"\n X: C step-length reduced below minimum.");
          DPRINTF(stdout,"\n theta: %g ;  theta_min: %g",dhunt_par,dhunt_par_min);
          DPRINTF(stdout,"\n    Program terminated.\n");
	  /* This needs to have a return value of 0, indicating
	     * success, for the continuation script to not treat this
	     * as a failed command. */
        goto free_and_clear;
          }
#ifdef PARALLEL
              check_parallel_error("\t");
#endif

	  /*
	   * ADJUST NATURAL PARAMETER
	   */

	for (iHC=0;iHC<nHC;iHC++) {
	  update_parameterHC(iHC, path1[iHC], x, xdot, x_AC, delta_s[iHC], cx, exo, dpi);
	}  /* end of iHC loop  */

        iHC = MAX(0,log_ID);
  	if(hunt[iHC].EndParameterValue == hunt[iHC].BegParameterValue)
 		{	hunt_par = 1.0;	}
	else
 		{
                 if(hunt[iHC].ramp == 2)
                    {
	             hunt_par = (log10(path1[iHC])-log10(hunt[iHC].BegParameterValue))
	                /(log10(hunt[iHC].EndParameterValue) - log10(hunt[iHC].BegParameterValue));
                    }  else   {
	  	     hunt_par = (path1[iHC]-hunt[iHC].BegParameterValue)
	     	        /(hunt[iHC].EndParameterValue - hunt[iHC].BegParameterValue)  ;
                    }
          /*      hunt_par=fabs(hunt_par);  */
	        dhunt_par = hunt_par-hunt_par_old;
 		}

	/*
	 * GET ZERO OR FIRST ORDER PREDICTION
	 */


	switch (Continuation) {
	case HUN_ZEROTH:
	    vcopy(numProcUnknowns, &x[0], 1.0, &x_old[0]);
	    break;
	case  HUN_FIRST:
	    v2sum(numProcUnknowns, &x[0], 1.0, &x_old[0], dhunt_par, &x_sens[0]);
            break;
	}

	/* MMH: Needed to put this in, o/w it may find that the
         * solution and residual HAPPEN to satisfy the convergence
         * criterion for the next newton solve...
         */
        find_and_set_Dirichlet(x, xdot, exo, dpi);

        exchange_dof(cx, dpi, x, pg->imtrx);

	if (nAC > 0)
          {
	    dcopy1(nAC, x_AC_old, x_AC);
	    for(iAC=0 ; iAC<nAC ; iAC++)
	      { update_parameterAC(iAC, x, xdot, x_AC, cx, exo, dpi); }
	  }

                iHC = MAX(0,log_ID);
  		if(hunt[iHC].EndParameterValue == hunt[iHC].BegParameterValue)
 			{	hunt_par = 1.0;	}
  		else
 			{
                         if(hunt[iHC].ramp == 2)
                             {
	                      hunt_par = (log10(path1[iHC])-log10(hunt[iHC].BegParameterValue))
	                        /(log10(hunt[iHC].EndParameterValue) - log10(hunt[iHC].BegParameterValue));
                             }  else   {
	  		      hunt_par = (path1[iHC]-hunt[iHC].BegParameterValue)
	      			/(hunt[iHC].EndParameterValue - hunt[iHC].BegParameterValue)  ;
                              hunt_par=fabs(hunt_par);
 			     }
                        }
 	}
 	else if (inewton == -1)
 	{
        DPRINTF(stdout,"\nHmm... trouble on first step \n  Let's try some more relaxation  \n");
 	      if((damp_factor1 <= 1. && damp_factor1 >= 0.) &&
 	         (damp_factor2 <= 1. && damp_factor2 >= 0.) &&
        		 (damp_factor3 <= 1. && damp_factor3 >= 0.))
 		{
 		custom_tol1 *= 0.01;
 		custom_tol2 *= 0.01;
 		custom_tol3 *= 0.01;
        DPRINTF(stdout,"  custom tolerances %g %g %g  \n",custom_tol1,custom_tol2,custom_tol3);
 		}
 		else
 		{
 		damp_factor1 *= 0.5;
        DPRINTF(stdout,"  damping factor %g  \n",damp_factor1);
 		}

 	    vcopy(numProcUnknowns, &x[0], 1.0, &x_old[0]);

 	/* MMH: Needed to put this in, o/w it may find that the
          * solution and residual HAPPEN to satisfy the convergence
          * criterion for the next newton solve...
          */
         find_and_set_Dirichlet(x, xdot, exo, dpi);

         exchange_dof(cx, dpi, x, pg->imtrx);


 	if (nAC > 0)
          {
 	    dcopy1(nAC, x_AC_old, x_AC);
 	    for(iAC=0 ; iAC<nAC ; iAC++)
 	      { update_parameterAC(iAC, x, xdot, x_AC, cx, exo, dpi); }
 	  }

 	}
 	else
 	{
        DPRINTF(stdout,"\nHmm... could not converge on first step\n Let's try some more iterations\n");
 	      if((damp_factor1 <= 1. && damp_factor1 >= 0.) &&
 	         (damp_factor2 <= 1. && damp_factor2 >= 0.) &&
        		 (damp_factor3 <= 1. && damp_factor3 >= 0.))
 		{
                if(hunt[0].BCID == -1)
                   {
	            if (Write_Intermediate_Solutions == 0) {
	                 write_solution(ExoFileOut, resid_vector, x, x_sens_p, x_old,
			        xdot, xdot_old, tev, tev_post, gv,  rd,
			        gvec, gvec_elem, &nprint, delta_s[0],
 			        theta, custom_tol1, NULL, exo, dpi);
	                 nprint++;
 	                 }
 	         }
 		custom_tol1 *= 100.;
 		custom_tol2 *= 100.;
 		custom_tol3 *= 100.;
        DPRINTF(stdout,"  custom tolerances %g %g %g  \n",custom_tol1,custom_tol2,custom_tol3);
 		}
 		else
 		{
                if(hunt[0].BCID == -1)
                   {
	            if (Write_Intermediate_Solutions == 0) {
        DPRINTF(stdout,"  writing solution %g  \n",damp_factor1);
	                 write_solution(ExoFileOut, resid_vector, x, x_sens_p, x_old,
			        xdot, xdot_old, tev, tev_post, gv,  rd,
			        gvec, gvec_elem, &nprint, delta_s[0],
 			        theta, damp_factor1, NULL, exo, dpi);
	                 nprint++;
 	                 }
 	         }
 		damp_factor1 *= 2.0;
		damp_factor1 = MIN(damp_factor1,1.0);
        DPRINTF(stdout,"  damping factor %g  \n",damp_factor1);
 		}
 	  }


      }  /* end of !converged */

    } while (converged == 0);

    /*
     * CONVERGED
     */
    nt++;
    custom_tol1 = toler_org[0];
    custom_tol2 = toler_org[1];
    custom_tol3 = toler_org[2];
    damp_factor1 = 1.0;
    DPRINTF(stdout,
	    "\n\tStep accepted, theta (proportion complete) = %10.6e\n",
	    hunt_par);
    for (iHC=0;iHC<nHC;iHC++) {
      switch (hunt[iHC].Type) {
      case 1:		/* BC */
          DPRINTF(stdout, "\tStep accepted, BCID=%3d DFID=%5d",
		  hunt[iHC].BCID, hunt[iHC].DFID);
	  break;
      case 2:		/* MT */
          DPRINTF(stdout, "\tStep accepted, MTID=%3d MPID=%5d",
		  hunt[iHC].MTID+1, hunt[iHC].MPID);
	  break;
      case 3:		/* AC */
          DPRINTF(stdout, "\tStep accepted, ACID=%3d DFID=%5d",
 		  hunt[iHC].BCID, hunt[iHC].DFID);
 	  break;
      }
      DPRINTF(stdout, " Parameter= % 10.6e\n", path1[iHC]);
    }

    /*
     * check path step error, if too large do not enlarge path step
     */

    iHC = MAX(0,log_ID);
    if ((ni == 1) && (n != 0) && (!const_delta_s[iHC]))
      {
       dhunt_par_new = path_step_control(num_total_nodes,
					     dhunt_par, dhunt_par_old,
					     x,
					     eps,
					     &success_ds,
					     cont->use_var_norm, inewton);
       if (dhunt_par_new > dhunt_par_max) {dhunt_par_new = dhunt_par_max;}
      } else {
	success_ds = 1;
	dhunt_par_new = dhunt_par;
      }

    /*
     * determine whether to print out the data or not
     */

    i_print = 0;
    if (nt == step_print) {
      i_print = 1;
      step_print += cont->print_freq; }

    if (alqALC == -1)
    { i_print = 1; }

    if (i_print) {
      error = write_ascii_soln(x, resid_vector, numProcUnknowns,
 			       x_AC, nAC, path1[0], file);
      if (error) {
	DPRINTF(stderr, "%s:  error writing ASCII soln file\n", yo);
      }
      if ( Write_Intermediate_Solutions == 0 ) {
        exo_time = aldALC[0]*path1[0];
	write_solution(ExoFileOut, resid_vector, x, x_sens_p,
		       x_old, xdot, xdot_old, tev, tev_post,  gv,
                       rd, gvec, gvec_elem, &nprint,
                       delta_s[0], theta, exo_time, NULL, exo, dpi);
	nprint++;
      }
    }

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
        step_fix += cont->fix_freq*cont->print_freq;
      }

    /*
     * backup old solutions
     * can use previous solutions for prediction one day
     */

    dcopy1(numProcUnknowns,x_older,x_oldest);
    dcopy1(numProcUnknowns,x_old,x_older);
    dcopy1(numProcUnknowns,x,x_old);

    dcopy1(nHC,delta_s_older,delta_s_oldest);
    dcopy1(nHC,delta_s_old  ,delta_s_older );
    dcopy1(nHC,delta_s      ,delta_s_old   );
    dcopy1(nHC,delta_s_new  ,delta_s       );
    dhunt_par_old = dhunt_par;
    dhunt_par = dhunt_par_new;
    hunt_par_old = hunt_par;
    if ( nAC > 0) {
      dcopy1(nAC, x_AC, x_AC_old);
    }

    /*
     * INCREMENT/DECREMENT PARAMETER
     */

     hunt_par += dhunt_par;
     if(hunt_par > 1.0)
          {
          dhunt_par = 1.0 - (hunt_par - dhunt_par);
          hunt_par = 1.0;
          }
     for (iHC=0;iHC<nHC;iHC++) {
           if(hunt[iHC].ramp == 2)
             {
              path1[iHC] = hunt[iHC].BegParameterValue*
                   pow(hunt[iHC].EndParameterValue/hunt[iHC].BegParameterValue,hunt_par);
              path[iHC] = hunt[iHC].BegParameterValue*
              pow(hunt[iHC].EndParameterValue/hunt[iHC].BegParameterValue,hunt_par-dhunt_par);
             }
           else
             {
           path1[iHC] = hunt[iHC].BegParameterValue + hunt_par*
                (hunt[iHC].EndParameterValue - hunt[iHC].BegParameterValue);
           path[iHC] = hunt[iHC].BegParameterValue + (hunt_par-dhunt_par)*
                (hunt[iHC].EndParameterValue - hunt[iHC].BegParameterValue);
             }
          }
      /*
       * ADJUST NATURAL PARAMETER
       */

    for (iHC=0;iHC<nHC;iHC++) {
      update_parameterHC(iHC, path1[iHC], x, xdot, x_AC, delta_s[iHC], cx, exo, dpi);
    }  /*  end of iHC loop */

    /*
     * GET FIRST ORDER PREDICTION
     */

    switch (Continuation) {
    case HUN_ZEROTH:
	break;
    case  HUN_FIRST:
	v1add(numProcUnknowns, &x[0], dhunt_par, &x_sens[0]);
        break; }

        if (!good_mesh) goto free_and_clear;

    /*
     *
     * CHECK END CONTINUATION
     *
     */

    if (alqALC == -1)
    { alqALC = 0; }
    else
    { alqALC = 1; }

    if (alqALC == 0) {
      DPRINTF(stdout,"\n\n\t I will continue no more!\n\t No more continuation for you!\n");
      goto free_and_clear;
    }

  } /* n */

      if(n == MaxPathSteps &&
	 aldALC[0] * (lambdaEnd[0] - path[0]) > 0)
	{
	  DPRINTF(stderr,"\n\tFailed to reach end of hunt in maximum number of successful steps (%d).\n\tSorry.\n",
		  MaxPathSteps);
          goto free_and_clear;
 	  /*exit(0);  */
	}
#ifdef PARALLEL
      check_parallel_error("Hunting error");
#endif

  /*
   * DONE CONTINUATION
   */

 free_and_clear:

  /*
   * Transform the node point coordinates according to the
   * displacements and write out all the results using the
   * displaced coordinates. Set the displacement field to
   * zero, too.
   */

  if (Anneal_Mesh) {
#ifdef DEBUG
    fprintf(stderr, "%s: anneal_mesh()...\n", yo);
#endif
    err = anneal_mesh(x, tev, tev_post, NULL, rd, path1[0], exo, dpi);
#ifdef DEBUG
    DPRINTF(stderr, "%s: anneal_mesh()-done\n", yo);
#endif
    EH(err, "anneal_mesh() bad return.");
  }

  /*
   * Free a bunch of variables that aren't needed anymore
   */
  safer_free((void **) &ROT_Types);
  safer_free((void **) &node_to_fill);

  safer_free( (void **) &resid_vector);
  safer_free( (void **) &resid_vector_sens);
  safer_free( (void **) &scale);
  safer_free( (void **) &x);

  if (nAC > 0) {
    safer_free( (void **) &x_AC);
    safer_free( (void **) &x_AC_old);
    safer_free( (void **) &x_AC_dot);
  }

  safer_free( (void **) &x_old);
  safer_free( (void **) &x_older);
  safer_free( (void **) &x_oldest);
  safer_free( (void **) &xdot);
  safer_free( (void **) &xdot_old);
  safer_free( (void **) &x_update);

  safer_free( (void **) &x_sens);

  if((nn_post_data_sens+nn_post_fluxes_sens) > 0)
          Dmatrix_death(x_sens_p,num_pvector,numProcUnknowns);

  for(i = 0; i < MAX_NUMBER_MATLS; i++) {
    for(n = 0; n < MAX_MODES; n++) {
      safer_free((void **) &(ve_glob[i][n]->gn));
      safer_free((void **) &(ve_glob[i][n]));
    }
    safer_free((void **) &(vn_glob[i]));
  }

  sl_free(matrix_systems_mask, ams);

  for (i=0;i<NUM_ALSS;i++) {
    safer_free( (void**) &(ams[i]));
  }

  safer_free( (void **) &gvec);

  safer_free( (void **) &lambda);
  safer_free( (void **) &lambdaEnd);
  safer_free( (void **) &path);
  safer_free( (void **) &path1);
  safer_free( (void **) &hDelta_s0);
  safer_free( (void **) &hDelta_s_min);
  safer_free( (void **) &hDelta_s_max);
  safer_free( (void **) &delta_s);
  safer_free( (void **) &delta_s_new);
  safer_free( (void **) &delta_s_old);
  safer_free( (void **) &delta_s_older);
  safer_free( (void **) &delta_s_oldest);

  Ivector_death(&aldALC[0], nHC);
  Ivector_death(&const_delta_s[0], nHC);

  i = 0;
  for ( eb_indx = 0; eb_indx < exo->num_elem_blocks; eb_indx++ ) {
    for ( ev_indx = 0; ev_indx < rd->nev; ev_indx++ ) {
      if ( exo->elem_var_tab[i++] == 1 ) {
        safer_free ((void **) &(gvec_elem [eb_indx][ev_indx]) );
      }
    }
    safer_free ((void **) &(gvec_elem [eb_indx]));
  }

  safer_free( (void **) &gvec_elem);

  safer_free( (void **) &rd);
  safer_free( (void **) &Local_Offset);
  safer_free( (void **) &Dolphin);

  if( strlen( Soln_OutFile)  )
    {
       fclose(file);
    }

  free(gv);
  free(pg->matrices);

  return;

} /* END of routine hunt_problem  */
/*****************************************************************************/
/*****************************************************************************/
/*  END of file ac_hunt.c  */
/*****************************************************************************/
