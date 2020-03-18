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

/* Constructs and solves linear systems associated with finding 
 * solution stability. */

/*
 *$Id: ac_stability.c,v 5.1 2007-09-18 18:53:40 prschun Exp $
 */

/*
 * Revision history has goneaway from this place.  Major revisions
 * listed below. 
 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define GOMA_AC_STABILITY_C
#include "ac_stability.h"

#include "sl_eggroll.h"
#include "dpi.h"
#include "exo_struct.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "mm_fill.h"
#include "mm_unknown_map.h"
#include "rd_mesh.h"
#include "rf_allo.h"
#include "rf_bc.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_io.h"
#include "rf_mp.h"
#include "rf_solver.h"
#include "rf_util.h"
#include "sl_eggroll_def.h"
#include "sl_matrix_util.h"
#include "sl_util_structs.h"
#include "std.h"

struct Results_Description;
struct elem_side_bc_struct;

struct Results_Description;
struct elem_side_bc_struct;

/* Make LSA_MATRIX_OUTPUT_TOLERANCE negative if you want to include
 * the zeros in the sparse structure.  Otherwise, they are filtered
 * out. */
#define LSA_MATRIX_OUTPUT_TOLERANCE 1e-15
#define LSA_COMPARE 0

int Linear_Stability = LSA_NONE;
int LSA_3D_of_2D_pass = 0;
dbl LSA_3D_of_2D_wave_number = -1.0;
int LSA_number_wave_numbers = 0;
dbl *LSA_wave_numbers = NULL;
int LSA_current_wave_number = 0;

/* Calculate matrices needed for stability analysis
 * and solve for spectrum is required.
 *
 * Modified by Ian Gates Sep 26, 1997.
 *
 * Modified by MMH 1/15/00.  Hacked out store_stability_problem.
 *
 * Modified by MMH 11/28/00.  Changed the call list so we could pass
 * x_dot, x_old, etc., to matrix_fill.  Some of the call list probably
 * isn't making any difference but I'd rather err on the cautious
 * side...
 *
 * Modified by MMH 12/18/00.  3D stability of 2D flow begins!  Broke
 * the (single) function into multiple parts because it was getting
 * too fat...
 */
int
solve_stability_problem(struct Aztec_Linear_Solver_System *ams,                 
			double x[],			/* Value of the solution vector */
			double delta_t,			/* Time step size */
			double theta,			/* Variable time integration parameter
							   explicit (theta = 1) to 
							   implicit (theta = 0) */
			double resid_vector[],
			double x_old[],			/* Value of the old solution vector */
			double x_older[],		/* */
			double xdot[],			/* Value of xdot predicted for new 
							   solution */
			double xdot_old[],		/* */
			double x_update[],		/* */
			int *converged,			/* Whether the Newton has converged */
			int *nprint,			/* Counter for time step number */
			int tnv,			/* Number of nodal results */
			int tnv_post,			/* Number of post processing results */
			struct Results_Description  *rd,
			int *gindex,
			int *p_gsize,
			double *gvec,
			double time_value,
			Exo_DB *exo,			/* Ptr to finite element mesh db */
			Dpi *dpi)			/* Ptr to distributed processing info */
{
  int res;

  fprintf(stderr, "\nWARNING: Linear Stability Analysis is not (and cannot be) compatible\n");
  fprintf(stderr, "  with all user-supplied boundary conditions (i.e., BC cards with\n");
  fprintf(stderr, "  the term USER in them).  If you are using some, and you want LSA,\n");
  fprintf(stderr, "  then you must modify your boundary conditions accordingly.\n\n");

  /* MMH:
   * Here's the trick!  x_old was not the same as x when it came time
   * to do transient calculations.  So things were off by a small
   * amount (a few nodes with < 1% difference, say).  But this small
   * difference caused large errors in the eigenvalues (10%).  I fixed
   * this 11/28/00. */
  dcopy1(NumUnknowns[pg->imtrx],x,x_old);
  dcopy1(NumUnknowns[pg->imtrx],x,x_older);

  if(!(Linear_Stability == LSA_3D_OF_2D || Linear_Stability == LSA_3D_OF_2D_SAVE))
    res = solve_full_stability_problem(ams, x, delta_t, theta, resid_vector, 
				       x_old, x_older, xdot, xdot_old,
				       x_update, converged, nprint, tnv,
				       tnv_post, rd, gindex, p_gsize, gvec,
				       time_value, exo, dpi);
  else
    res = solve_3D_of_2D_stability_problem(ams, x, delta_t, theta,
					   resid_vector, x_old, x_older,
					   xdot, xdot_old, x_update,
					   converged, nprint, tnv, tnv_post,
					   rd, gindex, p_gsize, gvec,
					   time_value, exo, dpi);
  return(res);
}/* END of routine solve_stability_problem */

/* This routine will perform all of the eigenvalue handling for
 * "regular" systems.  That means it is NOT for 3D stability of a 2D
 * flow. */
int
solve_full_stability_problem(struct Aztec_Linear_Solver_System *ams,                 
			     double x[], /* Value of the solution vector */
			     double delta_t, /* Time step size */
			     double theta, /* Variable time integration parameter
					    * explicit (theta = 1) to implicit
					    * (theta = 0) */
			     double resid_vector[],
			     double x_old[], /* Value of the old solution vector */
			     double x_older[], /* */
			     double xdot[],  /* Value of xdot predicted for new 
					      * solution */
			     double xdot_old[], /* */
			     double x_update[], /* */
			     int *converged, /* Whether the Newton has converged */
			     int *nprint, /* Counter for time step number */
			     int tnv, /* Number of nodal results */
			     int tnv_post, /* Number of post processing results */
			     struct Results_Description  *rd,
			     int *gindex,
			     int *p_gsize,
			     double *gvec,
			     double time_value,
			     Exo_DB *exo, /* Ptr to finite element mesh db */
			     Dpi *dpi) /* Ptr to distributed processing info */
{
  int istuff[PARAMETER_ARRAY_LENGTH];
  double dstuff[PARAMETER_ARRAY_LENGTH];

  int i, j, mn;
  int num_total_nodes;		/* Number of nodes that each processor is
				   responsible for */
  double *scale, *zero;
  double double_temp;

  double *jacobian_matrix;
  double *mass_matrix, *mass_matrix_tmpA, *mass_matrix_tmpB;
  int *ija = ams->bindx;	/* This structure is the same for ALL matrices... */

  /* Initialize... */
  zero = calloc(NumUnknowns[pg->imtrx], sizeof(double));
  init_vec_value(&zero[0], 0.0, NumUnknowns[pg->imtrx]);

  scale = calloc(NumUnknowns[pg->imtrx], sizeof(double));
  init_vec_value(scale, 0.0, NumUnknowns[pg->imtrx]);

  mass_matrix = mass_matrix_tmpA = mass_matrix_tmpB = NULL;
  jacobian_matrix = NULL;

  /* Used to let the output routines know that we're not doing normal
   * mode analysis and just create one set of output files... */
  LSA_3D_of_2D_wave_number = -1.0;
  LSA_current_wave_number = 0;

  /* Allocate mass matrix */
  mass_matrix = calloc((NZeros+5), sizeof(double));
  if(LSA_COMPARE)
    {
      mass_matrix_tmpA = calloc((NZeros+5), sizeof(double));
      mass_matrix_tmpB = calloc((NZeros+5), sizeof(double));
    }

  num_total_nodes = dpi->num_universe_nodes;

  /* Get jacobian matrix */
  printf("Assembling J...\n");
  jacobian_matrix = ams->val;
  init_vec_value(&jacobian_matrix[0], 0.0, (NZeros+1));
  
  af->Assemble_Residual = TRUE;
  af->Assemble_Jacobian = TRUE;
  af->Assemble_LSA_Jacobian_Matrix = TRUE;
  af->Assemble_LSA_Mass_Matrix = FALSE;
  
  matrix_fill_full(ams, &x[0], &resid_vector[0], 
		   x_old, x_older, xdot, xdot_old, x_update,
		   &delta_t, &theta, 
		   First_Elem_Side_BC_Array[pg->imtrx], 
		   &time_value, exo, dpi, 
                   &num_total_nodes, &zero[0], &zero[0], NULL);

  /* This call will fill in scale[] with the proper row scales.  I
   * need to get this in terms of the Jacobian, and the apply it to
   * the "mass" matrix, later... */
  row_sum_scaling_scale(ams, resid_vector, scale);

  /* Now compute mass matrix. */

  /* MMH: 
   * theta = 0 makes the method implicit, so that we get the
   * var_{}^{n+1} coefficient. */
  theta = 0.0;

  /* The other way: */
  if(LSA_COMPARE) /* For a comparison, this is the Subtraction method... */
    {
      init_vec_value(&mass_matrix_tmpA[0], 0.0, (NZeros+1));
      init_vec_value(&mass_matrix_tmpB[0], 0.0, (NZeros+1));

      TimeIntegration = TRANSIENT;
      for(mn = 0; mn < upd->Num_Mat; mn++) 
	{
	  pd_glob[mn]->TimeIntegration = TRANSIENT;
	  for(i = 0; i < MAX_EQNS; i++)
	    if(pd_glob[mn]->e[pg->imtrx][i])
	      {
		pd_glob[mn]->e[pg->imtrx][i] |= T_MASS;
		pd_glob[mn]->etm[pg->imtrx][i][LOG2_MASS] = 1.0;
	      }
	}

      printf("Assembling B_tmpA...\n");
      af->Assemble_Residual = TRUE;
      af->Assemble_Jacobian = TRUE;
      af->Assemble_LSA_Jacobian_Matrix = FALSE;
      af->Assemble_LSA_Mass_Matrix = FALSE;
      ams->val = mass_matrix_tmpA;
      delta_t = 1.0;		/* To get the phi_i*phi_j terms
				 * (transient soln.) */
      tran->delta_t = 1.0;      /*for Newmark-Beta terms in Lagrangian Solid*/

      matrix_fill_full(ams, x, resid_vector, 
		       x_old, x_older, xdot, xdot_old, x_update,
		       &delta_t, &theta, 
		       First_Elem_Side_BC_Array[pg->imtrx], 
		       &time_value, exo, dpi, 
		       &num_total_nodes, 
                       &zero[0], &zero[0], NULL);

      printf("Assembling B_tmpB...\n");
      af->Assemble_Residual = TRUE;
      af->Assemble_Jacobian = TRUE;
      af->Assemble_LSA_Jacobian_Matrix = FALSE;
      af->Assemble_LSA_Mass_Matrix = FALSE;
      ams->val = mass_matrix_tmpB;
      delta_t = 1.0e+72;	/* To simulate steady-state soln. */
      tran->delta_t = 1.0e+72;  /*for Newmark-Beta terms in Lagrangian Solid*/
      matrix_fill_full(ams, x, resid_vector, 
		       x_old, x_older, xdot, xdot_old, x_update,
		       &delta_t, &theta, 
		       First_Elem_Side_BC_Array[pg->imtrx], 
		       &time_value, exo, dpi, 
		       &num_total_nodes, 
                       &zero[0], &zero[0], NULL);
      
      for(i = 0; i < NZeros+1; i++)
	{
	  double_temp = mass_matrix_tmpB[i]-mass_matrix_tmpA[i];
	  /*
	    if(double_temp != 0.0)
	    if(fabs(double_temp/mass_matrix_tmpA[i]) < 1e-2)
	    printf("Index %d is near-annihilatory.  value = %g, difference = %g, ratio = %g\n",
	    i, mass_matrix_tmpA[i], double_temp, fabs(double_temp/mass_matrix_tmpA[i]));
	  */
	  mass_matrix_tmpA[i] = double_temp;
	}
    }

  /* Initialize */
  for(mn = 0; mn < upd->Num_Mat; mn++) 
    {
      pd_glob[mn]->TimeIntegration = TRANSIENT;
      for(i = 0; i < MAX_EQNS; i++)
	/* MMH: Note that there is apparently no support for
	 * multiple species equations. */
	if(pd_glob[mn]->e[pg->imtrx][i])
	  {
	    /*
	      printf("Equation %d is on.\n",i);
	    */
	    pd_glob[mn]->e[pg->imtrx][i] = T_MASS;
	    for (j=0;j<MAX_TERM_TYPES;j++)
	      pd_glob[mn]->etm[pg->imtrx][i][j] = 0.0;
	    pd_glob[mn]->etm[pg->imtrx][i][LOG2_MASS] = 1.0;
	  }
    }

  /* Hopefully, with all of these settings, this should compute the
   * mass matrix in mass_matrix and ija */
  init_vec_value(&mass_matrix[0], 0.0, (NZeros+1));

  printf("Assembling B...\n"); 
  af->Assemble_Residual = TRUE;
  af->Assemble_Jacobian = TRUE;
  af->Assemble_LSA_Jacobian_Matrix = FALSE;
  af->Assemble_LSA_Mass_Matrix = TRUE;
  ams->val = mass_matrix;
  delta_t   = 1.0;
  tran->delta_t = 1.0;      /*for Newmark-Beta terms in Lagrangian Solid*/
  matrix_fill_full(ams, x, resid_vector, 
		   x_old, x_older, xdot, xdot_old, x_update,
		   &delta_t, &theta, First_Elem_Side_BC_Array[pg->imtrx], 
		   &time_value, exo, dpi, &num_total_nodes, 
                   &zero[0], &zero[0], NULL);

  /* MMH: Curious... Why NZeros + 1 ??? */
      /*
  for(i = 0; i < NZeros+1; i++)
    {
      mass_matrix[i] = -mass_matrix[i];
      mass_matrix[i] = mass_matrix_tmpA[i];
    }
      */

  /* Now scale the mass matrix with the jacobian scaling we already
   * computed.
   */
  matrix_scaling(ams, mass_matrix, -1.0, scale);
  if (LSA_COMPARE)
    {
      matrix_scaling(ams, mass_matrix_tmpA, -1.0, scale);

      compare_mass_matrices(jacobian_matrix, mass_matrix,
			    mass_matrix_tmpA, ija, NumUnknowns[pg->imtrx]);
    }

  /* Print jacobian and mass matrices */
  if(Linear_Stability == LSA_SAVE || eigen->Eigen_Matrix_Output)
    output_stability_matrices(mass_matrix, jacobian_matrix, ija,
			      num_total_nodes, NumUnknowns[pg->imtrx], NZeros);
  
  if(Linear_Stability != LSA_SAVE)
    {
      /* Call eggroll initiator */
      for(i = 0; i < PARAMETER_ARRAY_LENGTH; i++)
	{
	  istuff[i] = -1;
	  dstuff[i] = -1.0;
	}
      eggroll_init(NumUnknowns[pg->imtrx], NZeros, &istuff[0], &dstuff[0]);
      
      /* Call eggroll solver */
      eggrollwrap(istuff, dstuff, ija, jacobian_matrix, mass_matrix,
		  x, ExoFileOut, ProblemType, delta_t,  theta, x_old, 
		  xdot, xdot_old, resid_vector, converged, nprint, tnv, 
		  tnv_post, rd, gindex, p_gsize, gvec, 
		  time_value, exo, Num_Proc, dpi);
    }
      
  /* Free space, and return ams->val to what it was when we came
   * in... */
  free(mass_matrix);
  if(LSA_COMPARE)
    {
      free(mass_matrix_tmpA);
      free(mass_matrix_tmpB);
    }
  free(scale);
  free(zero);
  ams->val = jacobian_matrix;
  return(0);
}

/* This routine will perform all of the eigenvalue handling for
 * calculations of 3D stability of a 2D flow.
 */
int
solve_3D_of_2D_stability_problem(struct Aztec_Linear_Solver_System *ams,                 
				 double x[], /* Value of the solution vector */
				 double delta_t, /* Time step size */
				 double theta, /* Variable time integration parameter
						* explicit (theta = 1) to implicit
						* (theta = 0) */
				 double resid_vector[],
				 double x_old[], /* Value of the old solution vector */
				 double x_older[], /* */
				 double xdot[],  /* Value of xdot predicted for new 
						  * solution */
				 double xdot_old[], /* */
				 double x_update[], /* */
				 int *converged, /* Whether the Newton has converged */
				 int *nprint, /* Counter for time step number */
				 int tnv, /* Number of nodal results */
				 int tnv_post, /* Number of post processing results */
				 struct Results_Description  *rd,
				 int *gindex,
				 int *p_gsize,
				 double *gvec,
				 double time_value,
				 Exo_DB *exo, /* Ptr to finite element mesh db */
				 Dpi *dpi) /* Ptr to distributed processing info */
{
  int istuff[PARAMETER_ARRAY_LENGTH];
  double dstuff[PARAMETER_ARRAY_LENGTH];

  int i, j, mn, wn;
  int num_total_nodes;		/* Number of nodes that each processor is
				   responsible for */
  double *scale, *zero;

  double *jacobian_matrix, *jacobian_matrix_1, *jacobian_matrix_2;
  double *mass_matrix, *mass_matrix_1, *mass_matrix_2, *mass_matrix_1_tmpA,
    *mass_matrix_1_tmpB, *mass_matrix_2_tmpA, *mass_matrix_2_tmpB;
  int *ija = ams->bindx;	/* This structure is the same for ALL matrices... */

  int **e_save;
  dbl ***etm_save;

  /* Initialize... */
  zero = calloc(NumUnknowns[pg->imtrx], sizeof(double));
  init_vec_value(&zero[0], 0.0, NumUnknowns[pg->imtrx]);

  scale = calloc(NumUnknowns[pg->imtrx], sizeof(double));
  init_vec_value(scale, 0.0, NumUnknowns[pg->imtrx]);

  mass_matrix = mass_matrix_1 = mass_matrix_2 = mass_matrix_1_tmpA =
    mass_matrix_1_tmpB = mass_matrix_2_tmpA = mass_matrix_2_tmpB =
    NULL;
  jacobian_matrix = jacobian_matrix_1 = jacobian_matrix_2 = NULL;

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

  /* Allocate mass matrices */
  mass_matrix = calloc((NZeros+5), sizeof(double));
  mass_matrix_1 = calloc((NZeros+5), sizeof(double));
  mass_matrix_2 = calloc((NZeros+5), sizeof(double));
  if(LSA_COMPARE)
    {
      mass_matrix_1_tmpA = calloc((NZeros+5), sizeof(double));
      mass_matrix_1_tmpB = calloc((NZeros+5), sizeof(double));
      mass_matrix_2_tmpA = calloc((NZeros+5), sizeof(double));
      mass_matrix_2_tmpB = calloc((NZeros+5), sizeof(double));
    }

  /* Allocate jacobian matrices */
  jacobian_matrix = calloc((NZeros+5), sizeof(double));
  /* jacobian_matrix_1 is going to use the ams->val space already
   * allocated... */
  jacobian_matrix_1 = ams->val;
  jacobian_matrix_2 = calloc((NZeros+5), sizeof(double));

  num_total_nodes = dpi->num_universe_nodes;

  for(wn = 0; wn < LSA_number_wave_numbers; wn++)
    {
      LSA_current_wave_number = wn;
      LSA_3D_of_2D_wave_number = LSA_wave_numbers[wn];
      printf("Solving for wave number %d out of %d.  WAVE NUMBER = %g\n", wn
	     + 1, LSA_number_wave_numbers, LSA_3D_of_2D_wave_number); 
      
      /* Get the original pd_glob[mn]->e[pg->imtrx][i] and pd_glob[mn]->etm[pg->imtrx][i][*] 
       * values back. */
      TimeIntegration = STEADY;
      theta = 0.0;
      for(mn = 0; mn < upd->Num_Mat; mn++)
	for(i = 0; i < MAX_EQNS; i++)
	  {
	    pd_glob[mn]->TimeIntegration = STEADY;
	    pd_glob[mn]->e[pg->imtrx][i] = e_save[mn][i];
	    for(j = 0; j < MAX_TERM_TYPES; j++)
	      pd_glob[mn]->etm[pg->imtrx][i][j] = etm_save[mn][i][j];
	  }
      
      /* Jacobian matrix, pass 1 */
      printf("Assembling J (pass 1)...\n");
      LSA_3D_of_2D_pass = 1;
      ams->val = jacobian_matrix_1;
      init_vec_value(&jacobian_matrix_1[0], 0.0, (NZeros+1));
      af->Assemble_Residual = TRUE;
      af->Assemble_Jacobian = TRUE;
      af->Assemble_LSA_Jacobian_Matrix = TRUE;
      af->Assemble_LSA_Mass_Matrix = FALSE;
      
      matrix_fill_full(ams, x, resid_vector, 
		       x_old, x_older, xdot, xdot_old, x_update,
		       &delta_t, &theta, First_Elem_Side_BC_Array[pg->imtrx], 
		       &time_value, exo, dpi, 
                       &num_total_nodes, &zero[0], &zero[0], NULL);

      /* Jacobian matrix, pass 2 */
      printf("Assembling J (pass 2)...\n");
      LSA_3D_of_2D_pass = 2;
      ams->val = jacobian_matrix_2;
      init_vec_value(&jacobian_matrix_2[0], 0.0, (NZeros+1));
      af->Assemble_Residual = TRUE;
      af->Assemble_Jacobian = TRUE;
      af->Assemble_LSA_Jacobian_Matrix = TRUE;
      af->Assemble_LSA_Mass_Matrix = FALSE;
      
      matrix_fill_full(ams, x, resid_vector, 
		       x_old, x_older, xdot, xdot_old, x_update,
		       &delta_t, &theta, 
		       First_Elem_Side_BC_Array[pg->imtrx], 
		       &time_value, exo, dpi, 
                       &num_total_nodes, &zero[0], &zero[0], NULL);

      /* Add the two passes together to get the "real" Jacobian
       * matrix. */
      for(i = 0; i < NZeros+1; i++)
	jacobian_matrix[i] = jacobian_matrix_1[i] +
	  jacobian_matrix_2[i];

      /* This call will fill in scale[] with the proper row scales.  I
       * need to get this in terms of the Jacobian, and the apply it to
       * the "mass" matrix, later... */
      ams->val = jacobian_matrix;
      row_sum_scaling_scale(ams, resid_vector, scale);

      /* Now compute the mass matrix.  This is more complicated than
       * it seems if we want to perform a comparison with the
       * Subtractionism method.  There are, in fact, three different
       * matrices with which we could compare along the way: mass
       * matrix from pass 1, mass matrix from pass 2, and final
       * (Additionism?) mass matrix.  I will construct comparisons
       * for all three...
       *
       * There is an order restriction.  We cannot recover the
       * pd_glob[mn]->e[imtrx][i] (and pd_glob[mn]->etm[imtrx][i][*]) values AFTER
       * we compute the mass matrix with my method, because we set
       * everything to 0 except the T_MASS-ish terms.  So we should
       * compute all of the Subtractionism matrices,
       * mass_matrix_[1,2]_tmp[A,b], first.  Then compute the two
       * mass_matrix_[1,2] passes.  Since we save the pd_glob[mn]
       * structures to reset things over multiple wave number
       * computations, we could refer to the saved values to remove
       * this restriction, but there's less code this way... */

      /* theta = 0 makes the method implicit, so that we get the
       * var_{}^{n+1} coefficient.  We want this for ALL mass matrix
       * computations, but not for the jacobian computations (they
       * should have the mass etm = 0, so it wouldn't matter...). */
      theta = 0.0;

      TimeIntegration = TRANSIENT;
      for(mn = 0; mn < upd->Num_Mat; mn++) 
	{
	  pd_glob[mn]->TimeIntegration = TRANSIENT;
	  for(i = 0; i < MAX_EQNS; i++)
	    if(pd_glob[mn]->e[pg->imtrx][i])
	      {
		pd_glob[mn]->e[pg->imtrx][i] |= T_MASS;
		pd_glob[mn]->etm[pg->imtrx][i][LOG2_MASS] = 1.0;
	      }
	}

      /* Subtractionism method for mass matrix, passes 1 and 2. */
      if(LSA_COMPARE)
	{
	  init_vec_value(&mass_matrix_1_tmpA[0], 0.0, (NZeros+1));
	  init_vec_value(&mass_matrix_1_tmpB[0], 0.0, (NZeros+1));
	  init_vec_value(&mass_matrix_2_tmpA[0], 0.0, (NZeros+1));
	  init_vec_value(&mass_matrix_2_tmpB[0], 0.0, (NZeros+1));

	  /* They all have the same action flag settings... */
	  af->Assemble_Residual = TRUE;
	  af->Assemble_Jacobian = TRUE;
	  af->Assemble_LSA_Jacobian_Matrix = FALSE;
	  af->Assemble_LSA_Mass_Matrix = FALSE;

	  printf("Assembling B_tmpA (pass 1)...\n");
	  LSA_3D_of_2D_pass = 1;
	  ams->val = mass_matrix_1_tmpA;
	  delta_t = 1.0;	/* To get the phi_i*phi_j terms
				 * (transient soln.) */
	  tran->delta_t = 1.0;  /*for Newmark-Beta terms in Lagrangian Solid*/
	  matrix_fill_full(ams, x, resid_vector, x_old,
			   x_older, xdot, xdot_old,
			   x_update, &delta_t, &theta,
			   First_Elem_Side_BC_Array[pg->imtrx], &time_value, exo,
			   dpi, &num_total_nodes, &zero[0], &zero[0],
                           NULL);
	  
	  printf("Assembling B_tmpB (pass 1)...\n");
	  LSA_3D_of_2D_pass = 1;
	  ams->val = mass_matrix_1_tmpB;
	  delta_t = 1.0e+72;	/* To simulate steady-state soln. */
	  tran->delta_t = 1.0e+72;   /*for Newmark-Beta terms in Lagrangian Solid*/
	  matrix_fill_full(ams, x, resid_vector, x_old,
			   x_older, xdot, xdot_old, x_update,
			   &delta_t, &theta, First_Elem_Side_BC_Array[pg->imtrx],
			   &time_value, exo, dpi, &num_total_nodes,
                           &zero[0], &zero[0], NULL);
	  
	  printf("Assembling B_tmpA (pass 2)...\n");
	  LSA_3D_of_2D_pass = 2;
	  ams->val = mass_matrix_2_tmpA;
	  delta_t = 1.0;	/* To get the phi_i*phi_j terms
				 * (transient soln.) */
	  tran->delta_t = 1.0;      /*for Newmark-Beta terms in Lagrangian Solid*/
	  matrix_fill_full(ams, x, resid_vector, x_old,
			   x_older, xdot, xdot_old,
			   x_update, &delta_t, &theta,
			   First_Elem_Side_BC_Array[pg->imtrx], &time_value, exo,
			   dpi, &num_total_nodes, &zero[0], &zero[0],
                           NULL);
	  
	  printf("Assembling B_tmpB (pass 2)...\n");
	  LSA_3D_of_2D_pass = 2;
	  ams->val = mass_matrix_2_tmpB;
	  delta_t = 1.0e+72;	/* To simulate steady-state soln. */
	  tran->delta_t = 1.0e+72;  /*for Newmark-Beta terms in Lagrangian Solid*/
	  matrix_fill_full(ams, &x[0], &resid_vector[0], x_old,
			   x_older, xdot, xdot_old, x_update,
			   &delta_t, &theta, First_Elem_Side_BC_Array[pg->imtrx],
			   &time_value, exo, dpi, &num_total_nodes,
                           &zero[0], &zero[0], NULL);

	  /* Now construct the actual Subtractionism versions of each
	   * pass.  These are stored in mass_matrix_[1,2]_tmpA. */
	  for(i = 0; i < NZeros+1; i++)
	    {
	      mass_matrix_1_tmpA[i] = mass_matrix_1_tmpB[i] - mass_matrix_1_tmpA[i];
	      mass_matrix_2_tmpA[i] = mass_matrix_2_tmpB[i] - mass_matrix_2_tmpA[i];
	    }
	}
      
      /* Now we want to compute both passes of our mass matrix using
       * the action flag.  We already set some of the parameters
       * correcetly above... */

      /* Initialize */
      for(mn = 0; mn < upd->Num_Mat; mn++) 
	{
	  for(i = 0; i < MAX_EQNS; i++)
	    if(pd_glob[mn]->e[pg->imtrx][i])
	      {
		pd_glob[mn]->e[pg->imtrx][i] = T_MASS;
		for (j = 0; j < MAX_TERM_TYPES; j++)
		  pd_glob[mn]->etm[pg->imtrx][i][j] = 0.0;
		pd_glob[mn]->etm[pg->imtrx][i][LOG2_MASS] = 1.0;
	      }
	}

      init_vec_value(&mass_matrix[0], 0.0, (NZeros+1));
      init_vec_value(&mass_matrix_1[0], 0.0, (NZeros+1));
      init_vec_value(&mass_matrix_2[0], 0.0, (NZeros+1));
      
      /* These action flags and delta_t are the same for both
       * passes. */ 
      af->Assemble_Residual = TRUE;
      af->Assemble_Jacobian = TRUE;
      af->Assemble_LSA_Jacobian_Matrix = FALSE;
      af->Assemble_LSA_Mass_Matrix = TRUE;
      delta_t = 1.0;
      tran->delta_t = 1.0;      /*for Newmark-Beta terms in Lagrangian Solid*/

      printf("Assembling B (pass 1)...\n"); 
      LSA_3D_of_2D_pass = 1;
      ams->val = mass_matrix_1;
      matrix_fill_full(ams, x, resid_vector, x_old, x_older,
		       xdot, xdot_old, x_update, &delta_t, &theta,
		       First_Elem_Side_BC_Array[pg->imtrx], &time_value, exo,
		       dpi, &num_total_nodes, &zero[0], &zero[0],
                       NULL);

      printf("Assembling B (pass 2)...\n"); 
      LSA_3D_of_2D_pass = 2;
      ams->val = mass_matrix_2;
      matrix_fill_full(ams, x, resid_vector, x_old, x_older,
		       xdot, xdot_old, x_update, &delta_t, &theta,
		       First_Elem_Side_BC_Array[pg->imtrx], &time_value, exo,
		       dpi, &num_total_nodes, &zero[0], &zero[0],
                       NULL);

      /* Get my mass matrix by adding the two passes together (don't
       * forget the -! */
      for(i = 0; i < NZeros+1; i++)
	mass_matrix[i] = -(mass_matrix_1[i] + mass_matrix_2[i]);

      /* Scale all of the preliminary pieces for comparison, if we're
       * doing one.  Otherwise, just scale the above mass_matrix... */
      if(LSA_COMPARE)
	{
	  for(i = 0; i < NumUnknowns[pg->imtrx]; i++)
	    {
	      mass_matrix_1_tmpA[i] /= scale[i];
	      mass_matrix_2_tmpA[i] /= scale[i];
	      mass_matrix_1[i] /= scale[i];
	      mass_matrix_2[i] /= scale[i];
	      mass_matrix[i] /= scale[i];
	      for(j = ija[i]; j < ija[i+1]; j++)
		{
		  mass_matrix_1_tmpA[j] /= scale[i];
		  mass_matrix_2_tmpA[j] /= scale[i];
		  mass_matrix_1[j] /= scale[i];
		  mass_matrix_2[j] /= scale[i];
		  mass_matrix[j] /= scale[i];
		}
	    }
	  /* Compare 1st pass results. */
	  printf("Comparing pass 1 of mass matrix with Subtractionism method...");
	  compare_mass_matrices(jacobian_matrix, mass_matrix_1,
				mass_matrix_1_tmpA, ija, NumUnknowns[pg->imtrx]);
	  
	  /* Comapre 2nd pass results. */
	  printf("Comparing pass 2 of mass matrix with Subtractionism method...");
	  compare_mass_matrices(jacobian_matrix, mass_matrix_2,
				mass_matrix_2_tmpA, ija, NumUnknowns[pg->imtrx]);

	  /* Get the sum of passes for Subtractionism method.  This
	   * was already done above for my regular method. */
	  for(i = 0; i < NZeros+1; i++)
	    mass_matrix_1_tmpA[i] = mass_matrix_1_tmpA[i] + mass_matrix_1_tmpB[i];

	  /* Compare sum of passes results. */
	  printf("Comparing sum of passes of mass matrix with Subtractionism method...");
	  compare_mass_matrices(jacobian_matrix, mass_matrix,
				mass_matrix_1_tmpA, ija, NumUnknowns[pg->imtrx]);
	}
      else
	for(i = 0; i < NumUnknowns[pg->imtrx]; i++)
	  {
	    mass_matrix[i] /= scale[i];
	    for(j = ija[i]; j < ija[i+1]; j++)
	      mass_matrix[j] /= scale[i];
	  }
      
      /* Print Jacobian and mass matrices */
      if(Linear_Stability == LSA_3D_OF_2D_SAVE || eigen->Eigen_Matrix_Output)
	output_stability_matrices(mass_matrix, jacobian_matrix, ija,
				  num_total_nodes, NumUnknowns[pg->imtrx],
				  NZeros);
      
      if(Linear_Stability != LSA_3D_OF_2D_SAVE)
	{
	  /* Call eggroll initiator */
	  for(i = 0; i < PARAMETER_ARRAY_LENGTH; i++)
	    {
	      istuff[i] = -1;
	      dstuff[i] = -1.0;
	    }
	  printf("Entering eggrollinit...\n"); fflush(stdout);
	  eggroll_init(NumUnknowns[pg->imtrx], NZeros, &istuff[0], &dstuff[0]);
	  
	  /* Call eggroll solver */
	  printf("Entering eggrollwrap...\n"); fflush(stdout);
	  eggrollwrap(istuff, dstuff, ija,
		      jacobian_matrix, mass_matrix, x,
		      ExoFileOut, ProblemType, delta_t, theta,
		      x_old, xdot, xdot_old, resid_vector,
		      converged, nprint, tnv, tnv_post, rd, gindex,
		      p_gsize, gvec, time_value, exo, Num_Proc, dpi); 
	}
    } /* for(wn = 0; wn < LSA_num_wave_numbers; wn++) */
      
  /* Free space, and return ams->val to what it was when we came
   * in... */
  free(jacobian_matrix_2);
  free(jacobian_matrix);
  if(LSA_COMPARE)
    {
      free(mass_matrix_2_tmpB);
      free(mass_matrix_2_tmpA);
      free(mass_matrix_1_tmpB);
      free(mass_matrix_1_tmpA);
    }
  free(mass_matrix_2);
  free(mass_matrix_1);
  free(mass_matrix);
  free(e_save);
  free(etm_save);
  free(scale);
  free(zero);
  ams->val = jacobian_matrix_1;
  return(0);
}

void output_stability_matrices(double *mass_matrix,
			       double *jacobian_matrix,
			       int *ija,
			       int num_total_nodes,
			       int NumUnknowns,
			       int NZeros)
{
  static int fake = 0;

  int i, j, k, num_B_nonzeros, num_J_nonzeros;
  FILE *mass_file, *jacobian_file, *vars_file;
  char LSA_mass_output_filename[256],
    LSA_jacobian_output_filename[256],
    LSA_vars_output_filename[256]; 

  /* Get the filenames.  These will be different if we're doing a
   * normal mode analysis for 3D stability of a 2D flow. */
  if(LSA_3D_of_2D_wave_number == -1.0)
    {
      sprintf(LSA_mass_output_filename, "LSA_mass_coo.out");
      sprintf(LSA_jacobian_output_filename, "LSA_jac_coo.out");
      sprintf(LSA_vars_output_filename, "LSA_vars.out");
    }
  else
    {
      sprintf(LSA_mass_output_filename, "LSA_mass_coo-%g.out",
	      (double)fake);
      sprintf(LSA_jacobian_output_filename, "LSA_jac_coo-%g.out",
	      (double)fake);
      sprintf(LSA_vars_output_filename, "LSA_vars-%g.out",
	      (double)fake++);
    }

/* Multiplex output names in parallel */
  if (Num_Proc > 1)
    {
      multiname(LSA_mass_output_filename, ProcID, Num_Proc);
      multiname(LSA_jacobian_output_filename, ProcID, Num_Proc);
      multiname(LSA_vars_output_filename, ProcID, Num_Proc);
    }

  printf("Writing matrix files (%s, %s, %s)...",
	 LSA_mass_output_filename,
	 LSA_jacobian_output_filename,
	 LSA_vars_output_filename);
  fflush(stdout);

  /* First, count out how many actual nonzeros there are.  This will
   * reduce filesizes substantially (especially for the mass
   * matrix.
   */
  num_J_nonzeros = 0;
  num_B_nonzeros = 0;
  for (i=0; i<NumUnknowns; i++)
    {
      if(fabs(jacobian_matrix[i]) > LSA_MATRIX_OUTPUT_TOLERANCE) num_J_nonzeros++;
      if(fabs(mass_matrix[i]) > LSA_MATRIX_OUTPUT_TOLERANCE) num_B_nonzeros++;
      for (j=ija[i]; j<ija[i+1]; j++)
	{
	  if(fabs(jacobian_matrix[j]) > LSA_MATRIX_OUTPUT_TOLERANCE) num_J_nonzeros++;
	  if(fabs(mass_matrix[j]) > LSA_MATRIX_OUTPUT_TOLERANCE) num_B_nonzeros++;
	}
    }

  mass_file = fopen(LSA_mass_output_filename, "w");
  if(!mass_file)
    EH(GOMA_ERROR, "Could not open mass matrix output file.");
  jacobian_file = fopen(LSA_jacobian_output_filename, "w");
  if(!jacobian_file)
    EH(GOMA_ERROR, "Could not open jacobian matrix output file.");
  
  /* Print out header line containing #rows, #nonzeros.
   */
  fprintf(jacobian_file, "%12d\n%12d\n", num_J_nonzeros, NumUnknowns);
  fprintf(mass_file, "%12d\n%12d\n", num_B_nonzeros, NumUnknowns);

  k = 0;
  for (i=0; i<NumUnknowns; i++)
    {
      if(fabs(mass_matrix[i]) > LSA_MATRIX_OUTPUT_TOLERANCE)
	fprintf(mass_file, " %5d %5d % 17.14e\n", i, i, mass_matrix[i]);
      if(fabs(jacobian_matrix[i]) > LSA_MATRIX_OUTPUT_TOLERANCE)
	fprintf(jacobian_file, " %5d %5d % 17.14e\n", i, i, jacobian_matrix[i]);
      k++;
      for (j=ija[i]; j<ija[i+1]; j++)
	{
	  if(fabs(mass_matrix[j]) > LSA_MATRIX_OUTPUT_TOLERANCE)
	    fprintf(mass_file, " %5d %5d % 17.14e\n", i, ija[j], mass_matrix[j]);
	  if(fabs(jacobian_matrix[j]) > LSA_MATRIX_OUTPUT_TOLERANCE)
	    fprintf(jacobian_file, " %5d %5d % 17.14e\n", i, ija[j], jacobian_matrix[j]);
	  k++;
	}
    }
  fclose(mass_file);
  fclose(jacobian_file);

  /* MMH Find out which components of the solution vector correspond
   * to which unknowns.  I know there's a better way to do this, but
   * I just don't know what it is.  Oh well...
   */
  vars_file = fopen(LSA_vars_output_filename, "w");
  if(!vars_file)
    EH(GOMA_ERROR, "Could not open variable listing output file.");
  for(i = 0; i < num_total_nodes; i++)
    {
      if( (j = Index_Solution(i, VELOCITY1, 0, 0, -1, pg->imtrx)) > -1 )
	fprintf(vars_file, "u %d\n", j);
      if( (j = Index_Solution(i, VELOCITY2, 0, 0, -1, pg->imtrx)) > -1 )
	fprintf(vars_file, "v %d\n", j);
      if( (j = Index_Solution(i, VELOCITY3, 0, 0, -1, pg->imtrx)) > -1 )
	fprintf(vars_file, "w %d\n", j);
      if( (j = Index_Solution(i, MESH_DISPLACEMENT1, 0, 0, -1, pg->imtrx)) > -1 )
	fprintf(vars_file, "x %d\n", j);
      if( (j = Index_Solution(i, MESH_DISPLACEMENT2, 0, 0, -1, pg->imtrx)) > -1 )
	fprintf(vars_file, "y %d\n", j);
      if( (j = Index_Solution(i, PRESSURE, 0, 0, -1, pg->imtrx)) > -1 )
	fprintf(vars_file, "p %d\n", j);
      if( (j = Index_Solution(i, PRESSURE, 0, 1, -1, pg->imtrx)) > -1 )
	fprintf(vars_file, "p %d\n", j);
      if( (j = Index_Solution(i, PRESSURE, 0, 2, -1, pg->imtrx)) > -1 )
	fprintf(vars_file, "p %d\n", j);
    }
  fclose(vars_file);

  puts("done.");
  printf("LSA (in) nnz = %d\n   (out) nnz = %d\n (out) J_nnz = %d\n (out) B_nnz = %d\n\n",
	 NZeros, k, num_J_nonzeros, num_B_nonzeros);
  fflush(stdout);
}

void compare_mass_matrices(double *jacobian_matrix,
			   double *mass_matrix,
			   double *mass_matrix_cmp,
			   int *ija,
			   int NumUnknowns)
{
  int i, j, k;
  int inode, i_Var_Desc, i_offset, idof;
  double *row_max_vector;

  row_max_vector = calloc(NumUnknowns, sizeof(double));

  printf("Checking for differences between two methods.\n");
  if(Iout != 1)
    printf("Iout (Output level) not equal to 1, cannot construct resname/dofname.\n");

  for(i = 0; i < NumUnknowns; i++)
    {
      row_max_vector[i] = fabs(jacobian_matrix[i]);
      if(fabs(mass_matrix[i]) > row_max_vector[i])
	row_max_vector[i] = fabs(mass_matrix[i]);
      if(fabs(mass_matrix_cmp[i]) > row_max_vector[i])
	row_max_vector[i] = fabs(mass_matrix_cmp[i]);
      for(j = ija[i]; j < ija[i+1]; j++)
	{
	  if(fabs(jacobian_matrix[j]) > row_max_vector[i])
	    row_max_vector[i] = fabs(jacobian_matrix[j]);
	  if(fabs(mass_matrix[j]) > row_max_vector[i])
	    row_max_vector[i] = fabs(mass_matrix[j]);
	  if(fabs(mass_matrix_cmp[j]) > row_max_vector[i])
	    row_max_vector[i] = fabs(mass_matrix_cmp[j]);
	}
      /* printf("row_max_vector[%d] = %g\n", i, row_max_vector[i]); */
    }

  /* Check for differences between the two methods. */
  k = 0;
  for(i = 0; i < NumUnknowns; i++)
    {
      if(!Index_Solution_Inv(i, &inode, &i_Var_Desc, &i_offset, &idof, pg->imtrx))
	printf("Uh-oh! Error on Index_Solution_Inv for i = %d\n", i);
      for(j = ija[i]; j < ija[i+1]; j++)
	{
	  if(fabs(mass_matrix[j] - mass_matrix_cmp[j]) > 1e-6 * row_max_vector[i])
	    {
	      printf("\nRow %d (0-based), column %d (0-based)\n", i, j);
	      printf("residual equation = %s\n", resname[pg->imtrx][i]);
	      printf("unknown column = %s\n", dofname[pg->imtrx][ija[NumUnknowns + 1 + k]]);
	      printf("mass_matrix = %g, mass_matrix_cmp = %g, row_max = %g\n",
		     mass_matrix[j], mass_matrix_cmp[j], row_max_vector[i]);
	    }
	  k++;
	}
    }

  free(row_max_vector);
}

void
eggroll_init(const int nj,	/* Number of equations */
 	     const int nnz_j,	/* Number of non-zeroes */
	     int *istuff,	/* Info for eigenvalue extraction */
	     dbl *dstuff)	/* Info for eigenvalue extraction */
/* This function has been moved here from the sl_ files */
{
  int i;

  istuff[0] = eigen->Eigen_Krylov_Subspace; /* dimension of krylov space */
  istuff[1] = nj;		/* number of equations */
  istuff[2] = nnz_j;		/* number of nonzeros */
  istuff[3] = eigen->Eigen_Filter;	/* number of filter steps */
  istuff[4] = eigen->Eigen_Recycle;	/* number of recylce steps */
  istuff[5] = 0;		/* cold start value (?) */
  istuff[6] = eigen->Eigen_NEV_WANT;	/* number of eigenpairs to find */
  istuff[7] = 4;		/* number of shifts given */
  istuff[8] = eigen->Eigen_Maximum_Iterations;	/* maximum number of iterations
*/
  istuff[9] = eigen->Eigen_Record_Modes; /* number of modes to write in Exodus files */

  dstuff[0] = eigen->Eigen_Tolerance;	/* solve to this tolerance.  It's not really doing this -mmh */
  dstuff[1] = 1.0;		/* f vector ? */
  dstuff[2] = 1.0;		/* m sign ? */
  dstuff[3] = eigen->Eigen_IV_Wt;	/* amount of uniform randomness to add to initial vector (0 = none) */
  for (i = 0; i < istuff[7]; i++)
    dstuff[10+i] = eigen->Eigen_Shifts[i]; /* initial shift values. */
}
