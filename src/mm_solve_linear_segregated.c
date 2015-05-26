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

#ifdef USE_RCSID
static char rcsid[] =
  "$Id: mm_sol_nonlinear.c,v 5.19 2010-07-21 21:03:04 sarober Exp $";
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include "std.h"
#include "rf_io_const.h"
#include "rf_io_structs.h"
#include "rf_io.h"
#include "rf_solver.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_mp.h"
#include "rf_allo.h"
#include "rf_bc_const.h"
#include "mm_more_utils.h"
  
#include "rf_solver_const.h"
#include "sl_matrix_util.h"

#include "mm_qp_storage.h"

#include "sl_util_structs.h"

#include "sl_amesos_interface.h"

#include "sl_aztecoo_interface.h"

#define _MM_SOLVE_LINEAR_SEGREGATED_C
#include "mm_solve_linear_segregated.h"
#include "goma.h"

/*
 * EDW: The prototype for function "mf_sol_lineqn" has been moved
 * to mm_sol_nonlinear.h so that LOCA has access to it.
 */


/*
 * This variable is used for more than one linear solver package.
 */

static int first_linear_solver_call=TRUE;


/*
 * Default: do not attempt to use Harwell MA28 linear solver. Kundert's is
 *          more robust and Harwell has a better successor to MA28 that you
 *          can buy with money.
 *
 * #define HARWELL
 */

#ifdef PARALLEL
#ifndef MPI
#define MPI			/* otherwise az_aztec.h trounces MPI_Request */
#endif
#endif

#include "az_aztec.h"

#include "sl_util.h"

#include "el_geom.h"

#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_post_def.h"

#include "mm_eh.h"

int
solve_linear_segregated(struct Aztec_Linear_Solver_System *ams,
			/* ptrs to Aztec linear systems */
			double x[],     /* soln vector on this proc */
			double delta_t, /* time step size */
			double theta,   /* parameter to vary time 
					 * integration from 
					 *   explicit (theta = 1) to 
					 *   implicit (theta = 0) */
			double x_old[], /* soln vector @ previous time */
			double x_older[], /* soln vector @ previous, previous time */
			double xdot[],  /* dxdt predicted for new time */
			double xdot_old[], /* dxdt for previous time */
			double resid_vector[],
			double x_update[], 
			double scale[],  /*Scale factor held for modified newton
					  *resolves */
			int *converged, /* whether the Newton iteration
					 * has converged (out) */
			int *nprint,    /* counter for time step number */
			double *glob_var_vals,   /* global variable values */
			/* details about post proc vars */
			double time_value,
			Exo_DB *exo,
			Dpi *dpi,
			Comm_Ex *cx,
			int nt,
			int *time_step_reform)
{
  static int prev_matrix = 0;

  double *a   = ams->val;	/* nonzero values of a CMSR matrix */
  int    *ija = ams->bindx;	/* column pointer array into matrix "a"*/

  int    *ija_save = ams->belfry; /* to hide external row/eqns from Aztec */

  static int    dofs_hidden=FALSE; /* boolean indicating when dofs are hidden */

  int 	        i, j, k;

  int           numProcUnknowns; /* number of degrees of freedom this processor
				  * sees (internal+boundary+external) but no
				  * more. */

  int           GNumUnknowns;		/* Global number of unknowns in the */
					/* system    */
  int 		inewton;		/* Newton iteration counter */

  int		return_value;	        /* nonzero if things screw up ...  */

  int		print_damp_factor;      /* Always printed after update norms */

  int		print_visc_sens;        /* Always printed after update norms */

  double	Norm[5][3];		/* Global norms... */
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

  double        Norm_r[2][3];           /* Relative (to solution vector) norms... */
                                        /*   [0][0] == correction, L_oo norm */
                                        /*   [0][1] == correction, L_1  norm */
                                        /*   [0][2] == correction, L_2  norm */
                                        /*   [1][0] == AC correction, L_oo norm */
                                        /*   [1][1] == AC correction, L_1  norm */
                                        /*   [1][2] == AC correction, L_2  norm */
  double        Norm_new;                /* Place holder for latest norm */
  double        Norm_old;                /* Place holder for last norm   */
  int           Norm_below_tolerance;    /* Boolean for modified newton test*/
  int           Rate_above_tolerance;    /* Boolean for modified newton test*/
  int continuation_converged = TRUE;
  int num_total_nodes = dpi->num_universe_nodes;
  /* Number of nodes that each processor is
   * responsible for                           */
  dbl h_elem_avg;                        /* global average element size for PSPG */
  dbl U_norm    ;                        /* global average velocity for PSPG */



  static char	yo[]="solve_nonlinear_problem";	/* routine identifier */

  /* 
   * The following vectors all have lengths equal to the number of local
   * unknowns - i.e., equations for which the current processor is responsible.
   */

  double 	*delta_x;		/* update */

  int num_unk_r, num_unk_x; 

  char dofname_r[80];
  char dofname_nr[80];
  char dofname_x[80];

  int       err = 0;

  dbl       total_mesh_volume;   /* global variable value */

  int local_order;		/* of the unknowns that this processor owns */
  int local_order_plus;		/* and if the external rows are included */
  int local_nnz;		/* number of nonzero matrix entries that
				 * are owned by this processor */
  int local_nnz_plus;		/* and if the external rows are included */

  int global_order;		/* order of the global system */
  int global_order_plus;	/* and if the external rows are overincluded */
  int global_nnz;		/* a sum of the number of nonzero matrix 
				 * entries owned by each processor */
  int global_nnz_plus;		/* a sum that overincludes the external rows */

  char		stringer[80];	/* holding format of num linear solve itns */
  dbl		a_start;	/* mark start of assembly */
  dbl		a_end;		/* mark end of assembly */

  char		ctod[80];	/* hold current time of day */

  dbl		s_start;	/* mark start of solve */
  dbl		s_end;		/* mark end of solve */

  double asmslv_start;
  double asmslv_end;
  int           num_unk_g, num_unk_y;
  dbl           *res_p, *res_m;


#ifdef DEBUG_MMH
  int inode, i_Var_Desc, i_offset, idof_eqn, idof_var;
  dbl abs_row_sum, row_sum;
  VARIABLE_DESCRIPTION_STRUCT *vd_eqn, *vd_var;
#endif /* DEBUG_MMH */

  if (pg->imtrx != prev_matrix) {
    first_linear_solver_call = TRUE;
    prev_matrix = pg->imtrx;
  }

  /*
   * Begin executable statements...
   */

  return_value = 0;		/* Set the return value to failure
				 * until convergence is achieved */

  /*
   * INITIALIZE
   */

  num_unk_r = num_unk_x = num_unk_g = num_unk_y = 0; 

  memset(Norm, 0, 3*5 * sizeof(double));
  memset(Norm_r, 0, 3*2 * sizeof(double));

  if (Debug_Flag) {
    DPRINTF(stderr, "%s: max Newton itns = %d\n", 
	    yo, Max_Newton_Steps);
    DPRINTF(stderr, "%s: convergence tol = %.1e %.1e \n",
	    yo, Epsilon[0], Epsilon[2]);
  }

  /*
   * The number of nonzero matrix entries that this processor owns.
   * External rows are not counted, though ija[] does stretch out to 
   * include them too for the parallel case. In the parallel case,
   * watch out to use the bigger ija[] for assembly purposes, but the
   * smaller one for Aztec (after processing via hide_external()).
   */
  numProcUnknowns = NumUnknowns[pg->imtrx] + NumExtUnknowns[pg->imtrx];
  if (strcmp(Matrix_Format, "msr") == 0) 
    {
      local_order      = ams->npn;
      local_order_plus = ams->npn_plus;
      local_nnz_plus   = ams->nnz_plus;
      local_nnz        = ams->nnz;
      canine_chaos(local_order, local_order_plus, local_nnz, local_nnz_plus,
		   &global_order, &global_order_plus, 
		   &global_nnz, &global_nnz_plus);

      /*
       * Setup some legacy variables.
       */

      NZeros       = local_nnz;
      GNZeros      = global_nnz;
      GNumUnknowns = global_order;

      if (Debug_Flag) {
	DPRINTF(stderr, "%s: At matrix %d, NumUnknowns  = %d\n", yo, (pg->imtrx + 1), NumUnknowns[pg->imtrx]);
	DPRINTF(stderr, "%s: NZeros       = %d\n", yo, NZeros);
	DPRINTF(stderr, "%s: GNZeros      = %d\n", yo, GNZeros);
	DPRINTF(stderr, "%s: GNumUnknowns = %d\n", yo, GNumUnknowns);
      }      

      /*
       *  matrix_stats (a, ija, NumUnknowns[pg->imtrx], &NZeros, &GNZeros, &GNumUnknowns);
       */
#ifdef DEBUG
      fprintf(stderr, "P_%d: lo=%d, lo+=%d, lnnz=%d, lnnz+=%d\n", ProcID,
	      local_order, local_order_plus, local_nnz, local_nnz_plus);
      fprintf(stderr, "P_%d: go=%d, go+=%d, gnnz=%d, gnnz+=%d\n", ProcID,
	      global_order, global_order_plus, global_nnz, global_nnz_plus);
#endif /* DEBUG */
    }

  asdv(&delta_x, numProcUnknowns);
  asdv(&res_p, numProcUnknowns);
  asdv(&res_m, numProcUnknowns);

  /*
   * Initial conditions at the beginning of the Newton iteration loop...
   * -> Set the convervence flag to false
   */
  *converged = FALSE;
  inewton    = 0;
  if(Max_Newton_Steps <= 0) *converged = TRUE;

  if (Time_Jacobian_Reformation_stride > 1) {
    Norm_below_tolerance = TRUE;
    Rate_above_tolerance = TRUE;
  } else {
    Norm_below_tolerance = FALSE;
    Rate_above_tolerance = FALSE;
  }

  Norm_old             = 1.e+30;
  Norm_new             = 1.e+30;
  
  /*
   * Make sure every processor is "finished commenting" about the preliminary
   * problem setup - keep the iteration free of clutter caused by buffered
   * stdouts and stderrs from multiple processors...
   */
  fflush(stdout);
  fflush(stderr);
#ifdef PARALLEL  
  (void) MPI_Barrier(MPI_COMM_WORLD);
#endif /* PARALLEL */
  fflush(stdout);
  fflush(stderr);

  if ( TimeIntegration == STEADY )
    {

      DPRINTF(stderr, "\n\n");
      DPRINTF(stderr, 
	      "               R e s i d u a l         C o r r e c t i o n\n");
    }
  else
    {

      DPRINTF(stderr, 
	      "\n    N e w t o n  C o n v e r g e n c e  - I m p l i c i t   T i m e   S t e p\n");
    }

  DPRINTF(stderr, 
	  "\n  ToD    itn   L_oo    L_1     L_2     L_oo    L_1     L_2   lis  asm/slv (sec)\n");
  DPRINTF(stderr, 
	  "-------- --- ------- ------- ------- ------- ------- ------- --- ---------------\n");

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

  asmslv_start = ut(); asmslv_end = asmslv_start;

  log_msg("%s: Newton iteration %d", yo, inewton);

  if ( inewton < 10 )
    {
      DPRINTF(stderr,"%s [%d] ", ctod, inewton );
    }
  else if ( inewton < 100 )
    {
      DPRINTF(stderr,"%s %d] ", ctod, inewton );
    }
  else
    {
      DPRINTF(stderr,"%s %d ", ctod, inewton );	  
    }
     

  h_elem_avg = 0.;
  U_norm     = 0.;
  

  if (!Norm_below_tolerance || !Rate_above_tolerance)
    {
      init_vec_value (resid_vector, 0.0, numProcUnknowns);
      init_vec_value (a, 0.0, ( NZeros + 1 ) );
      af->Assemble_Residual = TRUE;
      af->Assemble_Jacobian = TRUE;
      af->Assemble_LSA_Jacobian_Matrix = FALSE;
      af->Assemble_LSA_Mass_Matrix = FALSE;
    }
  else
    {
      init_vec_value(resid_vector, 0.0, numProcUnknowns);
      af->Assemble_Residual = TRUE;
      af->Assemble_Jacobian = FALSE;
      af->Assemble_LSA_Jacobian_Matrix = FALSE;
      af->Assemble_LSA_Mass_Matrix = FALSE;
    }
  a_start = ut(); a_end = a_start;

  /* Former block 0 in mm_fill.c. Here is some initialization */
  num_total_nodes = dpi->num_universe_nodes;

  /* If using volume change metric for element quality, reset here */
  if (nEQM>0 && eqm->do_vol && af->Assemble_Residual)
    {
      eqm->vol_sum = 0.0;
      eqm->vol_low = 9999.9;
      eqm->vol_count = 0;
    }

  /*
   * NORMAL, TANGENT and OTHER Vectors required for ROTATION are calculated
   * ahead of time so we don't run into anomolous behavior due to neclaced 
   * elements, junction points, . . .
   */
  if (Num_ROT > 0) calculate_all_rotation_vectors(exo, x);
  else  if ( Use_2D_Rotation_Vectors == TRUE ) calculate_2D_rotation_vectors(exo,x);

  /* Initialize volume constraint, not the most efficient way. */

  numProcUnknowns = NumUnknowns[pg->imtrx] + NumExtUnknowns[pg->imtrx];

  /* Exchange dof before matrix fill so parallel information
     is properly communicated */
  exchange_dof(cx,dpi, x);

  err = matrix_fill_full(ams, x, resid_vector, 
			 x_old, x_older, xdot, xdot_old, x_update,
			 &delta_t, &theta, 
			 First_Elem_Side_BC_Array, 
			 &time_value, exo, dpi,
			 &num_total_nodes,
			 &h_elem_avg, &U_norm, NULL);
	      
  a_end = ut();
  if (err == -1) {
    return_value = -1;
    goto free_and_clear;
  }

  /* Scale matrix first to get rid of problems with 
   * penalty parameter. In front option this is done
   * within the solver
   */
  row_sum_scaling_scale(ams, resid_vector, scale);

  if (err == -1) {
    return_value = -1;
    goto free_and_clear;
  }

  /*
   * We already know the resid vector now, so let's print out its norms
   * now instead of after a big long matrix elimination...
   */
  print_damp_factor = FALSE;
  print_visc_sens = FALSE;
  Norm[0][0] = Loo_norm(resid_vector, NumUnknowns[pg->imtrx], &num_unk_r, dofname_r);
  Norm[0][1] = L1_norm (resid_vector, NumUnknowns[pg->imtrx]);
  Norm[0][2] = L2_norm (resid_vector, NumUnknowns[pg->imtrx]);

  log_msg("%-38s = %23.16e", "residual norm (L_oo)", Norm[0][0]);
  log_msg("%-38s = %23.16e", "residual norm (L_1)", Norm[0][1]);
  log_msg("%-38s = %23.16e", "residual norm (L_2)", Norm[0][2]);

  DPRINTF(stderr, "%7.1e %7.1e %7.1e ",  Norm[0][0], Norm[0][1], Norm[0][2]);
		  
  /*
   * Before starting the solve of the linear system, fix up the matrix
   * to exclude external rows. However, we'll save all the column name 
   * values so we can reincorporate them back in before the next assembly 
   * step.
   */
  if (Num_Proc > 1 && strcmp( Matrix_Format, "msr" ) == 0 ) {
    hide_external(num_universe_dofs[pg->imtrx], NumUnknowns[pg->imtrx], ija, ija_save, a);
    dofs_hidden = TRUE;
  }

  /*************************************************************************
   *             SOLVE THE LINEAR SYSTEM
   *************************************************************************/
  s_start = ut(); s_end = s_start;
  linear_solve(ams, x, resid_vector, delta_x, stringer);
  s_end = ut();
  /**************************************************************************
   *        END OF LINEAR SYSTEM SOLVE SECTION
   **************************************************************************/
      
  Norm[1][0] = Loo_norm(delta_x, NumUnknowns[pg->imtrx], &num_unk_x, dofname_x);
  Norm[1][1] = L1_norm (delta_x, NumUnknowns[pg->imtrx]);
  Norm[1][2] = L2_norm (delta_x, NumUnknowns[pg->imtrx]);
  Norm_r[0][0] = Loo_norm_r(delta_x, x, NumUnknowns[pg->imtrx], &num_unk_x,dofname_nr);
  Norm_r[0][1] = L1_norm_r (delta_x, x, NumUnknowns[pg->imtrx]);
  Norm_r[0][2] = L2_norm_r (delta_x, x, NumUnknowns[pg->imtrx]);
      
  /* Save some norm info for modified Newton stuff */
      
  Norm_old = Norm_new;
  Norm_new = Norm[0][1];

  /* fail if we didn't get a finite solution */
  if (!finite(Norm[1][0]) || !finite(Norm[1][1]) || !finite(Norm[1][2]) ||
      !finite(Norm[0][0]) || !finite(Norm[0][1]) || !finite(Norm[0][2])) {
    return_value = -1;
    goto free_and_clear;
  }

  log_msg("%-38s = %23.16e", "correction norm (L_oo)", Norm[1][0]);
  log_msg("%-38s = %23.16e", "correction norm (L_1)", Norm[1][1]);
  log_msg("%-38s = %23.16e", "correction norm (L_2)", Norm[1][2]);

  asmslv_end = ut();

  fflush(NULL);

  if((damp_factor1 <= 1. && damp_factor1 >= 0.) &&
     (damp_factor2 <= 1. && damp_factor2 >= 0.) &&
     (damp_factor3 <= 1. && damp_factor3 >= 0.))
    {
      print_damp_factor = TRUE;
      if( !Visc_Sens_Copy ) Include_Visc_Sens = FALSE;
      if(Norm[0][0] > custom_tol3)
	{
	  damp_factor = damp_factor3;
	  /*
	    DPRINTF(stderr, "\n Invoking Damping factor %f\n", damp_factor3);
	  */
	}
      else if(Norm[0][0] > custom_tol2)
	{
	  damp_factor = damp_factor2;
	  /*
	    DPRINTF(stderr, "\n Invoking Damping factor %f\n", damp_factor2);
	  */
	}
      else if(Norm[0][0] > custom_tol1)
	{
	  damp_factor = damp_factor1;
	  /*
	    DPRINTF(stderr, "\n Invoking Damping factor %f\n", damp_factor1);
	  */
	}
      else
	{
	  damp_factor = 1.00;
	  /*
	    DPRINTF(stderr, " \n Invoking Damping factor %f", damp_factor);
	  */
	  if( !Visc_Sens_Copy ) 
	    {
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
    }
  else
    {
      /*default damping factor case */
      if(damp_factor2 == -1.) damp_factor = damp_factor1; 
      if( !Visc_Sens_Copy )
	{
	  if(Visc_Sens_Factor*inewton < Max_Newton_Steps)
	    { Include_Visc_Sens = FALSE; }
	  else
	    {
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

  /*******************************************************************
   *
   *   UPDATE GOMA UNKNOWNS
   *
   *******************************************************************/
  for (i = 0; i < NumUnknowns[pg->imtrx]; i++) {
    x[i] -= damp_factor * var_damp[idv[pg->imtrx][i][0]] * delta_x[i];
  }
  exchange_dof(cx, dpi, x);
  if (pd->TimeIntegration != STEADY) {
    for (i = 0; i < NumUnknowns[pg->imtrx]; i++) {
      xdot[i] -= damp_factor * var_damp[idv[pg->imtrx][i][0]] * delta_x[i] * (1.0 + 2 * theta) / delta_t;
    }
    exchange_dof(cx, dpi, xdot);
  }
  /*******************************************************************
   *
   * copy out global unknowns if array as been allocated for them 
   *
   *
   *******************************************************************/

  if ( glob_var_vals != NULL )
    {
      glob_var_vals[0] = (double) *converged;
      glob_var_vals[1] = (double) inewton;
      glob_var_vals[2] = (double) Max_Newton_Steps;
      if (Norm_new > 0.0 ) {
	glob_var_vals[3] = (double) (log10(Norm_new)/log10(Norm_old));
      } else {
	glob_var_vals[3] = 0.;
      }
      total_mesh_volume = 0;

      glob_var_vals[4] = (double) total_mesh_volume;
    }

  /********************************************************************
   *  CHECK IF CONVERGED
   *   (EDW: Modified to check bordering algorithm convergence
   *         as required by LOCA)
   ********************************************************************/
      
  *converged = ((Norm[0][2] < Epsilon[0] && Norm[2][2] < Epsilon[0]) &&
		((Norm_r[0][2] + Norm_r[1][2]) < Epsilon[2]) &&
		(continuation_converged));

  /*
   * Back to normal. Readjust the ija[] pointers to mimic the bloated
   * system for this processor that pretends like it really cares about
   * filling in matrix entries for external degrees of freedom.
   *
   * If we are doing continuation then after convergence one must
   * defer this event until after performing the sensitivity calculation
   */
  if (Num_Proc > 1 && strcmp( Matrix_Format, "msr" ) == 0) 
    {
      show_external(num_universe_dofs[pg->imtrx], 
		    (num_universe_dofs[pg->imtrx] - num_external_dofs[pg->imtrx]), 
		    ija, ija_save, a);
      dofs_hidden = FALSE;
    }
	

  if(Epsilon[2] > 1)
    {
      if ( !(*converged) || (  Linear_Solver == FRONT ) || inewton == 0 ) {
	DPRINTF(stderr, "%7.1e %7.1e %7.1e %s ",
		Norm[1][0], Norm[1][1], Norm[1][2], stringer);
      }
      else {
#ifdef SKIP_LAST_SOLVE
	DPRINTF(stderr, "%23c %s ",' ', " ns");			
#else
	DPRINTF(stderr, "%7.1e %7.1e %7.1e %s ",
		Norm[1][0], Norm[1][1], Norm[1][2], stringer);
#endif
      }
    }
  else
    {
      DPRINTF(stderr, "%7.1e %7.1e %7.1e %s ",
	      Norm_r[0][0], Norm_r[0][1], Norm_r[0][2], stringer);
    }
	
	
  DPRINTF(stderr, "%7.1e/%7.1e\n", (a_end-a_start), (s_end-s_start));
  
  /* print damping factor and/or viscosity sens message here */
  if (print_damp_factor) DPRINTF(stderr, " Invoking damping factor %f\n", damp_factor);
  if (print_visc_sens) DPRINTF(stderr, " Invoking Viscosity Sensitivities\n");

  
  inewton++;
  
  return_value = inewton;

  /**
   *     Solution vector norms for detecting turning points
   */
  Norm[4][0] = Loo_norm(x, NumUnknowns[pg->imtrx], &num_unk_x, dofname_x);
  Norm[4][1] = L1_norm (x, NumUnknowns[pg->imtrx])/((double)NumUnknowns[pg->imtrx]);
  Norm[4][2] = L2_norm (x, NumUnknowns[pg->imtrx])/sqrt((double)NumUnknowns[pg->imtrx]);

  DPRINTF(stderr, "scaled solution norms  %13.6e %13.6e %13.6e \n", 
	  Norm[4][0], Norm[4][1], Norm[4][2]);


  if ( Num_Proc > 1 && strcmp( Matrix_Format, "msr" ) == 0 )
    {
      show_external(num_universe_dofs[pg->imtrx],
		    (num_universe_dofs[pg->imtrx] - num_external_dofs[pg->imtrx] ),
		    ija, ija_save, a);
      dofs_hidden = FALSE;
    }

 free_and_clear:

  safe_free( (void *) delta_x);
  safe_free( (void *) res_p);
  safe_free( (void *) res_m);

  return(return_value);

} /* end of solve_nonlinear_problem */


int linear_solve(struct Aztec_Linear_Solver_System *ams,
		 double x[],     /* soln vector on this proc */
		 double resid_vector[],
		 double delta_x[],
		 char stringer[])
{
  int why;
  double *a   = ams->val;	/* nonzero values of a CMSR matrix */
  int    *ija = ams->bindx;	/* column pointer array into matrix "a"*/

  int matrix_solved;

  static int UMF_system_id;
  static int		Factor_Flag;

  switch (Linear_Solver)
    {
    case UMFPACK2:
    case UMFPACK2F:
      if (strcmp(Matrix_Format, "msr"))
	EH(-1,"ERROR: umfpack solver needs msr matrix format");

      Factor_Flag = 1;
      if (first_linear_solver_call)
	{
	  Factor_Flag = 0;
	  UMF_system_id = -1;
	}

      /* Force refactorization if UMFPACK2F */
      if (Linear_Solver == UMFPACK2F) {
	Factor_Flag = 0;
      }
      int matr_form = 1;

      UMF_system_id = SL_UMF(UMF_system_id,
			     &first_linear_solver_call,
			     &Factor_Flag, &matr_form,
			     &NumUnknowns[pg->imtrx], &NZeros, &ija[0],
			     &ija[0], &a[0], &resid_vector[0],
			     &delta_x[0]);

      //      first_linear_solver_call = FALSE;

      LOCA_UMF_ID = UMF_system_id;
      strcpy(stringer, " 1 ");
      break;

    case SPARSE13a:
      if (strcmp(Matrix_Format, "msr")) {
	EH(-1,"ERROR: lu solver needs msr matrix format");
      }
      dcopy1(NumUnknowns[pg->imtrx], resid_vector, delta_x);
      lu (NumUnknowns[pg->imtrx], NumExtUnknowns[pg->imtrx], NZeros, 
	  a, ija, delta_x, (first_linear_solver_call?1:2));
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

      if ( save_old_A ) dcopy1(NZeros,ams->val, ams->val_old);

      AZ_solve(delta_x, resid_vector, ams->options, ams->params, 
	       ams->indx, ams->bindx, ams->rpntr, ams->cpntr, 
	       ams->bpntr, ams->val, ams->data_org, ams->status, 
	       ams->proc_config);

      first_linear_solver_call = FALSE;

      if ( Debug_Flag > 0 ) {
	dump_aztec_status(ams->status);
      }


      why = (int)ams->status[AZ_why];
      aztec_stringer(why, ams->status[AZ_its], &stringer[0]);

      matrix_solved = ( ams->status[AZ_why] == AZ_normal) ;
      break;	
          	  
    case AMESOS:

      if( strcmp( Matrix_Format,"msr" ) == 0 ) {
	amesos_solve_msr( Amesos_Package, ams, delta_x, resid_vector, 1 , pg->imtrx);
      } else if ( strcmp( Matrix_Format,"epetra" ) == 0 ) {
	amesos_solve_epetra(Amesos_Package, ams, delta_x, resid_vector, pg->imtrx);
      } else {
	EH(-1," Sorry, only MSR and Epetra matrix formats are currently supported with the Amesos solver suite\n");
      }
      strcpy(stringer, " 1 ");
      break;

    case AZTECOO:
      if ( strcmp( Matrix_Format,"epetra" ) == 0 ) {
	aztecoo_solve_epetra(ams, delta_x, resid_vector);
	why = (int) ams->status[AZ_why];
	aztec_stringer(why, ams->status[AZ_its], &stringer[0]);
	matrix_solved = (ams->status[AZ_why] == AZ_normal);
      } else {
	EH(-1, "Sorry, only Epetra matrix formats are currently supported with the AztecOO solver suite\n");
      }
      break;

    default:
      EH(-1, "That linear solver package is not implemented.");
      break;
    }

  /**************************************************************************
   *        END OF LINEAR SYSTEM SOLVE SECTION
   **************************************************************************/
  return 0;
}
