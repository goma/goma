
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
 *$Id: mm_numjac.c,v 5.5 2009-04-24 23:42:33 hkmoffa Exp $
 */


/* Standard include files */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

/* GOMA include files */

#include "mm_numjac.h"

#include "std.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_solver.h"
#include "el_elm.h"
#include "mm_eh.h"
#include "rf_vars_const.h"
#include "mm_mp_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "rf_masks.h"
#include "rf_bc_const.h"
#include "rf_solver_const.h"
#include "sl_util.h"
#include "mm_qp_storage.h"
#include "dpi.h"
#include "el_elm_info.h"
#include "exo_struct.h"
#include "mm_fill.h"
#include "mm_fill_aux.h"
#include "mm_fill_ls.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_stress.h"
#include "mm_mp.h"
#include "mm_mp_structs.h"
#include "mpi.h"
#include "rd_mesh.h"
#include "rf_allo.h"
#include "rf_bc.h"
#include "rf_node_const.h"
#include "sl_matrix_util.h"
#include "sl_util_structs.h"
#include "mm_std_models.h"

#define GOMA_MM_NUMJAC_C

static void piksr2		/* mm_numjac.c                               */
(int ,			/* n                                         */
       int [],			/* arr                                       */
       int [],			/* brr                                       */
       dbl []);		/* crr                                       */

#ifdef FORWARD_DIFF_NUMJAC
static void compute_numerical_jacobian_errors
(dbl,			/* analytic value */
       dbl,			/* numerical value */
       dbl [],			/* absolute error */
       dbl []);		/* relative error */
#endif

typedef struct {
  int num_colors;
  int *column_color;
  int *colptr;
  int *rowptr;
  int nnz;
} Coloring;

static Coloring* find_coloring(struct Aztec_Linear_Solver_System *ams,
			int num_unknowns,
			int num_total_nodes,
			Exo_DB *exo,
			Dpi *dpi);

static void free_coloring(Coloring *coloring); 

typedef struct {
  int a_val, b_val;
  dbl c_val;
} data_t;

typedef struct intLinkedList {
  int val;
  struct intLinkedList *next;
} IntLinkedList;

/* Prepend an item to the list, replacing the current head
 * value with the new value */
static void int_linked_list_prepend(IntLinkedList * head, int val)
{
  // This is kinda hacky to not use **list, any ptr to old head items
  // will now have updated val
  IntLinkedList * item = malloc(sizeof(IntLinkedList));
  item->next = head->next;
  item->val = head->val;
  head->next = item;
  head->val = val;
}

/* Check if an item is in the list */
static int item_in_int_linked_list(IntLinkedList * list, int val)
{
  if (list->val == val) {
    return TRUE;
  } else if (list->next == NULL) {
    return FALSE;
  } else {
    return item_in_int_linked_list(list->next, val);
  }
}

/* free memory for the list */
static void free_int_linked_list(IntLinkedList *list)
{
  if (list == NULL) return;
  free_int_linked_list(list->next);
  free(list);
}

/* Find the matrix coloring for finite difference

   Coloring is determined in a greedy fashion,
   might be useful to use actual algorithms in the future
   or an external library (ColPack)

   Color each coloumn such that all columns
   with the same color share no rows containing
   nonzeros in that column
*/
static Coloring* find_coloring(struct Aztec_Linear_Solver_System *ams,
			int num_unknowns,
			int num_total_nodes,
			Exo_DB *exo,
			Dpi *dpi)
{
  int i, j;
  /* for storing transpose in column format */
  int *col = calloc(sizeof(int), ams->bindx[num_unknowns]);
  int *row = calloc(sizeof(int), ams->bindx[num_unknowns]);
  int *colptr = calloc(sizeof(int), num_unknowns+1);
  int *rowidx;
  int *has_row;
  int *column_color;
  Coloring *coloring = malloc(sizeof(Coloring));
#ifdef DEBUG_FD_COLORING
  double t1, t2;
  t1 = MPI_Wtime();
#endif


  /* Find Coo and CSC format, probably not great for giant problems,
   * but we are doing numerical jacobians which
   * is probably the biggest limit to problem size here
   * Storage costs increase by ~ 3 * nnz */
  int nnz = 0;
  for (i = 0; i < num_unknowns; i++) {
    row[nnz] = i;
    col[nnz] = i;
    nnz++;
    for (j = ams->bindx[i]; j < ams->bindx[i+1]; j++) {
      row[nnz] = i;
      col[nnz] = ams->bindx[j];
      nnz++;
    }
  }

  rowidx = calloc(sizeof(int), nnz);
  /* Find CSC format */
  for (i = 0; i < nnz; i++) {
    colptr[col[i] + 1] += 1;
  }

  // fix indexes
  for (i = 0; i < num_unknowns; i++) {
    colptr[i+1] += colptr[i];
  }

  int *counts = calloc(sizeof(int), num_unknowns);

  // set rows for each column
  for (i = 0; i < nnz; i++) {
    j = col[i];
    int offset = colptr[j] + counts[j];
    rowidx[offset] = row[i];
    counts[j] += 1;
  }

  free(counts);
  free(row);
  free(col);
#ifdef DEBUG_FD_COLORING
  t2 = MPI_Wtime();
  printf( "%d Elapsed time is %f\n", ProcID, t2 - t1 );
#endif

  column_color = malloc(sizeof(int) * num_unknowns);

  for (j = 0; j < num_unknowns; j++) {
    column_color[j] = -1;
  }

  /* greedy coloring
   *
   * Assign column to a color only if that
   * color is not already associated with rows
   * in the column
   */
  int num_colored = 0;
  int num_colors = 0;
  has_row = calloc(num_unknowns, sizeof(int));
  // until all columns have a color
  while (num_colored < num_unknowns) {
    // Color ever column we can with the current color
    for (j = 0; j < num_unknowns; j++) {
      if (column_color[j] == -1) {
	// Only valid if all rows are not in this color
	int valid = TRUE;
	for (i = colptr[j]; i < colptr[j+1]; i++) {
	  if (has_row[rowidx[i]]) {
	    valid = FALSE;
	    break;
	  }
	}

	// We can add this column to the current color
	if (valid) {
	  for (i = colptr[j]; i < colptr[j+1]; i++) {
	    has_row[rowidx[i]] = TRUE;
	  }
	  column_color[j] = num_colors;
	  num_colored++;
	}
      }
    }
    memset(has_row, 0, sizeof(int) * num_unknowns);
    num_colors++;
  }

  coloring->num_colors = num_colors;
  coloring->column_color = column_color;
#ifdef DEBUG_FD_COLORING
  t2 = MPI_Wtime();
  printf( "%d Elapsed2 time is %f\n", ProcID, t2 - t1 );
#endif

#if defined(PRINT_FD_COLORING) || defined(DEBUG_FD_COLORING)
  int min = num_unknowns, max = 0;
  double average = ((double) num_unknowns) / ((double) num_colors);

  for (i = 0; i < num_colors; i++) {
    int count = 0;
    for (j = 0; j < num_unknowns; j++) {
      if (column_color[j] == i){
	count++;
      }
    }
    if (count < min) {
      min = count;
    } else if (count > max) {
      max = count;
    }
  }
  // print some stats
  printf("Coloring results:\n");
  printf("Columns %d, Colors %d\n", num_unknowns, num_colors);
  printf("Min %d, Max %d, Average %g\n", min, max, average);
#endif
#ifdef DEBUG_FD_COLORING
  t2 = MPI_Wtime();
  printf( "%d Elapsed3 time is %f\n", ProcID, t2 - t1 );
#endif

  coloring->colptr = colptr;
  coloring->rowptr = rowidx;
  coloring->nnz = nnz;
  free(has_row);
  return coloring;
}

static void free_coloring(Coloring *coloring) {
  free(coloring->colptr);
  free(coloring->rowptr);
  free(coloring->column_color);
  free(coloring);
}

void
numerical_jacobian_compute_stress(struct Aztec_Linear_Solver_System *ams,
		   double x[],	/* Solution vector for the current processor */
		   double resid_vector[],   /* Residual vector for the current
					     * processor */
		   double delta_t, /* time step size */
		   double theta, /* parameter to vary time integration from
				    explicit (theta = 1) to
				    implicit (theta = 0) */
		   double x_old[], /* Value of the old solution vector */
		   double x_older[], /* Value of the real old soln vect */

		   double xdot[], /* Value of xdot predicted for new solution */
		   double xdot_old[], /* Value of xdot at previous time */

		   double x_update[],
		   int num_total_nodes,

		   struct elem_side_bc_struct *first_elem_side_BC_array[],
				/* This is an array of pointers to the first
				   surface integral defined for each element.
				   It has a length equal to the total number
				   of elements defined on the current proc */
		   int Debug_Flag, /* flag for calculating numerical jacobian
				      -1 == calc num jac w/o rescaling
				      -2 == calc num jac w/  rescaling */
		   double time_value, /* current value of time */
		   Exo_DB *exo,	    /* ptr to whole fe mesh */
		   Dpi *dpi,        /* any distributed processing info */
		   double *h_elem_avg,
		   double *U_norm)
/******************************************************************************
  This function compares the analytical jacobian entries calculated in
  matrix_fill the numerical ones approximated by central difference method.

  Author:          K. S. Chen (1511) (based on an earlier version by P. R. Schunk).
  Date:            January 19, 1994

  Updated: M. M. Hopkins.  Mucho optimization and other scaling options.
  Updated: D. R. Noble.  Added bracketed approach to checking Jacobian.
  Debug Option = -1 => unscaled rows
                 -2 => scaled rows
                 -3 => rows scaled by diagonal value

  Forward Difference Jacobian Checker: (USE -DFORWARD_DIFF_NUMJAC)
  If the absolute error exceeds RESIDUAL_TOLERANCE, or the scaled
  error exceeds SCALED_RESIDUAL_TOLERANCE (both defined in
  mm_numjac.h), an error message is reported.  Certain assumptions are
  made if the magnitudes of errors are on the order of
  SCALED_RESIDUAL_TOLERANCE_CUTOFF.  See the function
  compute_numerical_jacobian_errors() for details.

  Bracketed Jacobian Checker:
  Checks if the change in the residual caused by perturbing the solution vector
  is inconsistent with the two analytical jacobians computed for the original
  and perturbed solutions. Should catch any jacobian error that is larger than
  the change in the jacobian between the two solution points.  Should only return
  false positives when there is a change in the sign of the second derivative
  between the two solution points.
******************************************************************************/
{
  int i, j, nnonzero;
  int idx, color;
  int zeroCA = -1;
  double *resid_vector_1, *x_1;
  double dx;
  int *irow, *jcolumn, *nelem;
  int v_s[MAX_MODES][DIM][DIM];
  int mode;
  int my_elem_num, my_node_num;
  IntLinkedList *elem_list = NULL;
  int var_i, var_j;
  double *dx_col;
  double x_scale[MAX_VARIABLE_TYPES];
  int count[MAX_VARIABLE_TYPES];
  double *nj;
  char errstring[256];
  double *resid_vector_save;
  int numProcUnknowns = NumUnknowns[pg->imtrx] + NumExtUnknowns[pg->imtrx];
  // assuming we are only computing coloring once, and that the matrix does not change
  Coloring *coloring = NULL;

  if(strcmp(Matrix_Format, "msr"))
    EH(-1, "Cannot compute numerical jacobian values for non-MSR formats.");

/* calculates the total number of non-zero entries in the analytical jacobian, a[] */
  nnonzero = NZeros+1;

  /* allocate arrays to hold jacobian and vector values */
  irow = (int *) array_alloc(1, nnonzero, sizeof(int));
  jcolumn = (int *) array_alloc(1, nnonzero, sizeof(int));
  nelem = (int *) array_alloc(1, nnonzero, sizeof(int));
  resid_vector_1 =  (double *) array_alloc(1, numProcUnknowns, sizeof(double));
  x_1 =  (double *) array_alloc(1, numProcUnknowns, sizeof(double));
  dx_col = calloc(numProcUnknowns, sizeof(double));
  nj = calloc(sizeof(double), ams->bindx[numProcUnknowns]);
  memcpy(nj, ams->val, ams->bindx[numProcUnknowns]*(sizeof(double)));

  // Load stress equation pointers
  stress_eqn_pointer(v_s);

  resid_vector_save = (double *) array_alloc(1, numProcUnknowns, sizeof(double));
  memcpy(resid_vector_save, resid_vector, numProcUnknowns*(sizeof(double)));

  /* Cannot do this with Front */
  if (Linear_Solver == FRONT) EH(-1,"Cannot use frontal solver with numjac. Use umf or lu");

  /* Initialization */

  /* There are a couple of places in checking the Jacobian numerically
   * that you really need to know the scale of the unknowns in the problem.
   * One is to determine the right size fo a finite difference step and
   * the second is in evaluating the scale of the residual.  So first step
   * is to estimate the scale for all variables in the problem */
  memset(x_scale, 0, (MAX_VARIABLE_TYPES)*sizeof(dbl));
  memset(count, 0, (MAX_VARIABLE_TYPES)*sizeof(int));
  for (i = 0; i < numProcUnknowns; i++)
    {
      var_i = idv[pg->imtrx][i][0];
      count[var_i]++;
      x_scale[var_i] += x[i]*x[i];
    }
  dbl scale_tmp = global_h_elem_siz(x, x_old, xdot, resid_vector, exo, dpi);
  for (i = 0; i < MAX_VARIABLE_TYPES; i++)
    {
      if (count[i]) x_scale[i] = sqrt(x_scale[i]/count[i]);
      /* Now check for bad news.  If x[i] is zero everywhere then,
       * use the element size for displacements and for other
       * quantities assume x is order 1.
       */
      if (x_scale[i] == 0.)
        {
          switch (i)
            {
            case MESH_DISPLACEMENT1:
            case MESH_DISPLACEMENT2:
            case MESH_DISPLACEMENT3:
            case SOLID_DISPLACEMENT1:
            case SOLID_DISPLACEMENT2:
            case SOLID_DISPLACEMENT3:
              x_scale[i] = scale_tmp;
              break;

            default:
              x_scale[i] = 1.;
              break;
            }
        }
    }
  /* for level set problems we have an inherent scale */
  if (ls != NULL && ls->Length_Scale != 0.) x_scale[FILL] = ls->Length_Scale;

  /* copy x vector */
  memcpy(x_1, x, numProcUnknowns*(sizeof(double)));

  // Coloring is made the default as the cost of coloring is small
  coloring = find_coloring(ams, numProcUnknowns, num_total_nodes, exo, dpi);

  /*
   *  now calculate analytical and numerical jacobians at perturbed values
   *  check that the perturbed residuals are consistent with range possible
   *  for range of analytical jacobian values
   */
  for (color = 0; color < coloring->num_colors; color++)       /* loop over each color */
    {
#ifdef DEBUG_FD_COLORING
      double t1, t2, t3, t4;
      t1 = MPI_Wtime();
#endif
      memset(resid_vector_1, 0, NumUnknowns[pg->imtrx]*sizeof(dbl));
      /*
       * Perturb many variables
       */
#ifdef DEBUG_FD_COLORING
      int count = 0;
      int elem_count = 0;
#endif
      for (j = 0; j < numProcUnknowns; j++) {
	if (coloring->column_color[j] == color) {
	  dx = x_scale[idv[pg->imtrx][j][0]] * FD_DELTA_UNKNOWN;
	  if(dx < 1.0E-15) dx = 1.0E-7;
	  x_1[j] = x[j] + dx;

	  if (pd_glob[0]->TimeIntegration != STEADY) {
	    xdot[j] += (x_1[j] - x[j])  * (1.0 + 2 * theta) / delta_t;
	  }

	  dx_col[j] = dx;
#ifdef DEBUG_FD_COLORING
	  count++;
#endif
	  my_node_num = idv[pg->imtrx][j][2];

	  for(i = exo->node_elem_pntr[my_node_num];
	      i < exo->node_elem_pntr[my_node_num+1]; i++) {
	    my_elem_num = exo->node_elem_list[i];
	    if (elem_list == NULL) {
	      elem_list = malloc(sizeof(IntLinkedList));
	      elem_list->val = my_elem_num;
	      elem_list->next = NULL;
#ifdef DEBUG_FD_COLORING
	      elem_count++;
#endif
	    } else if (item_in_int_linked_list(elem_list, my_elem_num)) {
	      EH(-1, "Jacobian elem coloring error, trying to assemble already accounted for element");
	    } else {
	      int_linked_list_prepend(elem_list, my_elem_num);
#ifdef DEBUG_FD_COLORING
	      elem_count++;
#endif
	    }
	  }

	  for(i = exo->node_elem_pntr[my_node_num];
	      i < exo->node_elem_pntr[my_node_num+1]; i++) {
	    my_elem_num = exo->node_elem_list[i];
	    int k;
	    for(k = exo->elem_node_pntr[my_elem_num];
		k < exo->elem_node_pntr[my_elem_num+1]; k++) {
	      int node_num = exo->elem_node_list[k];
	      int l;
	      for(l = exo->node_elem_pntr[node_num];
		  l < exo->node_elem_pntr[node_num+1]; l++) {
		int elem_num = exo->node_elem_list[l];
		if(elem_num == -1) {
		  continue;
		}
		if(!item_in_int_linked_list(elem_list, elem_num)) {
		  int_linked_list_prepend(elem_list, elem_num);
#ifdef DEBUG_FD_COLORING
		  elem_count++;
#endif
		}
	      }
	    }
	  }
	}
      }

#ifdef DEBUG_FD_COLORING
      printf("ColorStats %d [elem = %d], [count = %d]\n", color, elem_count, count);
      t2 = MPI_Wtime();
      printf( "%d Color0 time is %f\n", ProcID, t2 - t1 );
#endif

      af->Assemble_Residual = TRUE;
      af->Assemble_LSA_Jacobian_Matrix = FALSE;
      af->Assemble_LSA_Mass_Matrix = FALSE;
      af->Assemble_Jacobian = FALSE;
      neg_elem_volume = FALSE;
      neg_lub_height = FALSE;
      zero_detJ = FALSE;

      IntLinkedList *elptr;
      for (elptr = elem_list; elptr != NULL; elptr = elptr->next) {
	int ielem = elptr->val;
	int ebn;
	/*First we must calculate the material-referenced element
	 *number so as to be compatible with the ElemStorage struct
	 */
	ebn = find_elemblock_index(ielem, exo);
	int mn = Matilda[ebn];
	if (mn < 0) {
	  continue;
	}

	/*needed for saturation hyst. func. */
	PRS_mat_ielem = ielem - exo->eb_ptr[ebn];

	matrix_fill_stress(ams, x_1, resid_vector_1,
			   x_old, x_older,  xdot, xdot_old, x_update,
			   &delta_t, &theta,
			   first_elem_side_BC_array,
			   &time_value, exo, dpi,
			   &ielem, &num_total_nodes,
			   h_elem_avg, U_norm, NULL, zeroCA);

      }

#ifdef DEBUG_FD_COLORING
      t3 = MPI_Wtime();
      printf( "%d Color1 time is %f\n", ProcID, t3 - t2 );
#endif
      /*
       * Free memory allocated above
       */
      global_qp_storage_destroy();

      for (j = 0; j < numProcUnknowns; j++) {
	if (color == coloring->column_color[j]) {
	  for (idx = coloring->colptr[j]; idx < coloring->colptr[j+1]; idx++) {
	    i = coloring->rowptr[idx];
	    var_i = idv[pg->imtrx][i][0];
	    var_j = idv[pg->imtrx][j][0];

	    for (mode=0; mode < vn->modes; mode++)
	      {
		/* Only for stress terms */
		if (idv[pg->imtrx][i][0] >= v_s[mode][0][0] && idv[pg->imtrx][i][0] <= v_s[mode][2][2])
		  {

		    if (Inter_Mask[pg->imtrx][var_i][var_j]) {

		      int ja = (i == j) ? j : in_list(j, ams->bindx[i], ams->bindx[i+1], ams->bindx);
		      if (ja == -1) {
			sprintf(errstring, "Index not found (%d, %d) for interaction (%d, %d)", i, j, idv[pg->imtrx][i][0], idv[pg->imtrx][j][0]);
			EH(ja, errstring);
		      }
		      nj[ja] = (resid_vector_1[i] - resid_vector[i]) / (dx_col[j]);
		    }
		  }
	      } // Loop over modes
	  }
	}
      }

      /*
       * return solution vector to its original state
       */
      for (j = 0; j < numProcUnknowns; j++) {
	if (coloring->column_color[j] == color) {
	  if (pd_glob[0]->TimeIntegration != STEADY) {
	    xdot[j] -= (x_1[j] - x[j])  * (1.0 + 2 * theta) / delta_t;
	  }
	}
      }
      memcpy(x_1, x, numProcUnknowns*(sizeof(double)));
      free_int_linked_list(elem_list);
      elem_list = NULL;
#ifdef DEBUG_FD_COLORING
      t4 = MPI_Wtime();
      printf( "%d Color2 time is %f\n", ProcID, t4 - t3 );
#endif
    }                          /* End of for (j=0; j<numProcUnknowns; j++) */

#ifdef PARALLEL
  neg_elem_volume_global = FALSE;
  MPI_Allreduce(&neg_elem_volume, &neg_elem_volume_global, 1,
		MPI_INT, MPI_LOR, MPI_COMM_WORLD);
  neg_elem_volume = neg_elem_volume_global;

  neg_lub_height_global = FALSE;
  MPI_Allreduce(&neg_lub_height, &neg_lub_height_global, 1,
		MPI_INT, MPI_LOR, MPI_COMM_WORLD);
  neg_lub_height = neg_lub_height_global;

  zero_detJ_global = FALSE;
  MPI_Allreduce(&zero_detJ, &zero_detJ_global, 1,
		MPI_INT, MPI_LOR, MPI_COMM_WORLD);
  zero_detJ = zero_detJ_global;
#endif

  if (neg_elem_volume) {
    DPRINTF(stderr, "neg_elem_volume triggered \n");
    exit(-1);
  }
  if (neg_lub_height) {
    DPRINTF(stderr, "neg_lub_height triggered \n");
    exit(-1);
  }
  if (zero_detJ) {
    DPRINTF(stderr, "zero_detJ triggered \n");
    exit(-1);
  }


  memcpy(ams->val, nj, ams->bindx[numProcUnknowns]*(sizeof(double)));
  free(nj);
  free(dx_col);
  memcpy(resid_vector, resid_vector_save, numProcUnknowns*(sizeof(double)));
  free(resid_vector_save);
  free_coloring(coloring);
  /* free arrays to hold jacobian and vector values */
  safe_free( (void *) irow) ;
  safe_free( (void *) jcolumn) ;
  safe_free( (void *) nelem) ;
  safe_free( (void *) resid_vector_1) ;
  safe_free( (void *) x_1) ;
}


/*
 * Global variables defined here. Declared frequently via rf_bc.h
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void
numerical_jacobian(struct Aztec_Linear_Solver_System *ams,
		   double x[],	/* Solution vector for the current processor */
		   double resid_vector[],   /* Residual vector for the current
					     * processor */
		   double delta_t, /* time step size */
		   double theta, /* parameter to vary time integration from
				    explicit (theta = 1) to
				    implicit (theta = 0) */
		   double x_old[], /* Value of the old solution vector */
		   double x_older[], /* Value of the real old soln vect */

		   double xdot[], /* Value of xdot predicted for new solution */
		   double xdot_old[], /* Value of xdot at previous time */

		   double x_update[],
		   int num_total_nodes,

		   struct elem_side_bc_struct *first_elem_side_BC_array[],
				/* This is an array of pointers to the first
				   surface integral defined for each element.
				   It has a length equal to the total number
				   of elements defined on the current proc */
		   int Debug_Flag, /* flag for calculating numerical jacobian
				      -1 == calc num jac w/o rescaling
				      -2 == calc num jac w/  rescaling */
		   double time_value, /* current value of time */
		   Exo_DB *exo,	    /* ptr to whole fe mesh */
		   Dpi *dpi,        /* any distributed processing info */
		   double *h_elem_avg,
		   double *U_norm)

/******************************************************************************
  This function compares the analytical jacobian entries calculated in
  matrix_fill the numerical ones approximated by central difference method.

  Author:          K. S. Chen (1511) (based on an earlier version by P. R. Schunk).
  Date:            January 19, 1994

  Updated: M. M. Hopkins.  Mucho optimization and other scaling options.
  Updated: D. R. Noble.  Added bracketed approach to checking Jacobian.
  Debug Option = -1 => unscaled rows
                 -2 => scaled rows
                 -3 => rows scaled by diagonal value

  Forward Difference Jacobian Checker: (USE -DFORWARD_DIFF_NUMJAC)
  If the absolute error exceeds RESIDUAL_TOLERANCE, or the scaled
  error exceeds SCALED_RESIDUAL_TOLERANCE (both defined in
  mm_numjac.h), an error message is reported.  Certain assumptions are
  made if the magnitudes of errors are on the order of
  SCALED_RESIDUAL_TOLERANCE_CUTOFF.  See the function
  compute_numerical_jacobian_errors() for details.

  Bracketed Jacobian Checker:
  Checks if the change in the residual caused by perturbing the solution vector
  is inconsistent with the two analytical jacobians computed for the original
  and perturbed solutions. Should catch any jacobian error that is larger than
  the change in the jacobian between the two solution points.  Should only return
  false positives when there is a change in the sign of the second derivative
  between the two solution points.
******************************************************************************/
{
  int i, j, k, l, m, ii, nn, kount, nnonzero, index;
  int zeroCA;
  double *a = ams->val;
  int *ija = ams->bindx;
  double *aj_diag, *aj_off_diag, *scale;
  double *resid_vector_1, *x_1, resid_scale;
  double resid_min, resid_max, resid_error;
  /* double resid_scaled_error; */
  double dx, delta_min, delta_max;
  double delta_aj_percentage, roundoff, confidence, resid_diff;
  int *irow, *jcolumn, *nelem;

  int num_elems, num_dofs;
  int my_elem_num, my_node_num, elem_num, node_num;
  int elem_already_listed;
  int *output_list;
  int *elem_list, *dof_list;
  NODE_INFO_STRUCT *node;
  NODAL_VARS_STRUCT *nvs;
  VARIABLE_DESCRIPTION_STRUCT *vd;
  int I, var_i, var_j, ibc, bc_input_id, eqn;
  struct elem_side_bc_struct *elem_side_bc;
	  double x_scale[MAX_VARIABLE_TYPES];
	  int count[MAX_VARIABLE_TYPES];
	  int Inter_Mask_save[MAX_VARIABLE_TYPES][MAX_VARIABLE_TYPES];
	#ifdef FORWARD_DIFF_NUMJAC
	  double *nj, nj_err, nj_scaled_err;
	#else
	  double *aj, *aj_1, nj;
	#endif

	#ifdef DEBUG_NUMJAC
	  dbl min_scale, max_scale, abs_min, abs_max;
#endif

  DPRINTF(stderr, "\n Starting Numerical Jacobian Checker\n");
  if(strcmp(Matrix_Format, "msr"))
    EH(-1, "Cannot compute numerical jacobian values for non-MSR formats.");

/* calculates the total number of non-zero entries in the analytical jacobian, a[] */
  nnonzero = NZeros+1;
  nn = ija[NumUnknowns[pg->imtrx]]-ija[0]; /* total number of diagonal entries a[] */

  /* allocate arrays to hold jacobian and vector values */
  irow = (int *) array_alloc(1, nnonzero, sizeof(int));
  jcolumn = (int *) array_alloc(1, nnonzero, sizeof(int));
  nelem = (int *) array_alloc(1, nnonzero, sizeof(int));
  aj_diag =  (double *) array_alloc(1, NumUnknowns[pg->imtrx], sizeof(double));
  aj_off_diag =  (double *) array_alloc(1, nnonzero, sizeof(double));
  resid_vector_1 =  (double *) array_alloc(1, NumUnknowns[pg->imtrx], sizeof(double));
  x_1 =  (double *) array_alloc(1, NumUnknowns[pg->imtrx], sizeof(double));
  scale =  (double *) array_alloc(1, NumUnknowns[pg->imtrx], sizeof(double));
  output_list = (int *)array_alloc(1, NumUnknowns[pg->imtrx], sizeof(int));
  dof_list = (int *)array_alloc(1, NumUnknowns[pg->imtrx], sizeof(int));
  elem_list = (int *)array_alloc(1, ELEM_LIST_SIZE, sizeof(int));
#ifdef FORWARD_DIFF_NUMJAC
  nj =  (double *) array_alloc(1, nnonzero, sizeof(double));
#else
  aj =  (double *) array_alloc(1, NumUnknowns[pg->imtrx], sizeof(double));
  aj_1 =  (double *) array_alloc(1, NumUnknowns[pg->imtrx], sizeof(double));
#endif

  if (aj_off_diag == NULL || scale == NULL) EH(-1, "No room for storage for computing numerical jacobian");

  /* Cannot do this with Front */
  if (Linear_Solver == FRONT) EH(-1,"Cannot use frontal solver with numjac. Use umf or lu");

  /* Initialization */
  memset(aj_off_diag, 0, nnonzero*sizeof(dbl));
  memset(aj_diag, 0, NumUnknowns[pg->imtrx]*sizeof(dbl));

  /* save Inter_Mask[pg->imtrx] away, turn on all entries so that we can make sure
   * that Inter_Mask[pg->imtrx] is being turned on for all entries being used
   */
  for(j =0; j < MAX_VARIABLE_TYPES; j++)
    {
      for(i = 0; i < MAX_VARIABLE_TYPES; i++)
        {
          Inter_Mask_save[j][i] = Inter_Mask[pg->imtrx][j][i];
          Inter_Mask[pg->imtrx][j][i] = 1;
        }
    }

  /* There are a couple of places in checking the Jacobian numerically
   * that you really need to know the scale of the unknowns in the problem.
   * One is to determine the right size fo a finite difference step and
   * the second is in evaluating the scale of the residual.  So first step
   * is to estimate the scale for all variables in the problem */
  memset(x_scale, 0, (MAX_VARIABLE_TYPES)*sizeof(dbl));
  memset(count, 0, (MAX_VARIABLE_TYPES)*sizeof(int));
  for (i = 0; i < NumUnknowns[pg->imtrx]; i++)
    {
      var_i = idv[pg->imtrx][i][0];
      count[var_i]++;
      x_scale[var_i] += x[i]*x[i];
    }
  for (i = 0; i < MAX_VARIABLE_TYPES; i++)
    {
      if (count[i]) x_scale[i] = sqrt(x_scale[i]/count[i]);
      /* Now check for bad news.  If x[i] is zero everywhere then,
       * use the element size for displacements and for other
       * quantities assume x is order 1.
       */
      if (x_scale[i] == 0.)
        {
          switch (i)
            {
            case MESH_DISPLACEMENT1:
            case MESH_DISPLACEMENT2:
            case MESH_DISPLACEMENT3:
            case SOLID_DISPLACEMENT1:
            case SOLID_DISPLACEMENT2:
            case SOLID_DISPLACEMENT3:
              x_scale[i] = global_h_elem_siz(x, x_old, xdot, resid_vector, exo, dpi);
              break;

            default:
              x_scale[i] = 1.;
              break;
            }
        }
    }
  /* for level set problems we have an inherent scale */
  if (ls != NULL && ls->Length_Scale != 0.) x_scale[FILL] = ls->Length_Scale;

  /* copy x vector */
  for (i = 0; i < NumUnknowns[pg->imtrx]; i++)
    {
      x_1[i] = x[i];
    }

  /* first calculate the residual vector corresponding to the solution vector read in
     the initial guess file; also calculate the analytical jacobian entries */
  af->Assemble_Residual = TRUE;
  af->Assemble_Jacobian = TRUE;
  af->Assemble_LSA_Jacobian_Matrix = FALSE;
  af->Assemble_LSA_Mass_Matrix = FALSE;

  DPRINTF(stderr, "Computing analytic entries ...");
  (void) matrix_fill_full(ams, x, resid_vector,
			  x_old, x_older, xdot, xdot_old,x_update,
			  &delta_t, &theta,
			  first_elem_side_BC_array,
			  &time_value, exo, dpi,
			  &num_total_nodes,
			  h_elem_avg, U_norm, NULL);

#ifdef DEBUG_NUMJAC
  DPRINTF(stderr, "Before scaling:\n");
  for(i = 0; i < 20; i++)
    DPRINTF(stderr, "resid[% 2d] = %-10.4g\n", i, resid_vector[i]);
#endif

  if(Debug_Flag == -2)
    {
      /* Scale matrix first to get rid of problems with
       * penalty parameter.
       */
      row_sum_scaling_scale(ams, resid_vector, scale);

#ifdef DEBUG_NUMJAC
      abs_min = 1.0e+10;
      abs_max = 0.0;
      for(i = 0; i < NumUnknowns[pg->imtrx]; i++)
	{
	  if(fabs(scale[i])>abs_max) abs_max=fabs(scale[i]);
	  if(fabs(scale[i])<abs_min) abs_min = fabs(scale[i]);
	}
      DPRINTF(stderr, "abs_min = %g, abs_max = %g\n", abs_min, abs_max);
#endif
    }
  if(Debug_Flag == -3)
    {
      /* Scale matrix by diagonal entry.  This is usually the largest
       * in magnitude.  If this is zero, then perform no scaling.
       */
      for(i = 0; i < NumUnknowns[pg->imtrx]; i++)
	scale[i] = (a[i] == 0.0) ? 1.0 : a[i];
      row_scaling(NumUnknowns[pg->imtrx], a, ija, resid_vector, scale);
    }

#ifdef DEBUG_NUMJAC
  DPRINTF(stderr, "After scaling:\n");
  for(i = 0; i < 20; i++)
    DPRINTF(stderr, "resid[% 2d] = %-10.4g\n", i, resid_vector[i]);
  min_scale = 1.0e+20;
  max_scale = -min_scale;
  DPRINTF(stderr, "Scale vector:\n");
  for(i = 0; i < NumUnknowns[pg->imtrx]; i++)
    {
      DPRINTF(stderr, "scale[% 2d] = %-10.4g\n", i, scale[i]);
      if(scale[i] < min_scale) min_scale = scale[i];
      if(scale[i] > max_scale) max_scale = scale[i];
    }
  DPRINTF(stderr, "min_scale = % 9.4g, max_scale = % 9.4g\n",
	  min_scale, max_scale);
#endif

  /* extract diagonal and off-diagonal elements from the coefficient matrix stored
     in sparse-storage format */
  for (i=0; i<NumUnknowns[pg->imtrx]; i++)
    aj_diag[i] = a[i];                    /* diagonal elements */

  kount=0;                              /* off-diagonal elements */
  for (i=0; i<NumUnknowns[pg->imtrx]; i++)
    {
      nelem[i] = ija[i+1] - ija[i];
      for (k=0; k<nelem[i]; k++)
	{
	  irow[kount]=i;                   /* row # in global jacobian matrix */
	  ii = kount + NumUnknowns[pg->imtrx] + 1;
	  jcolumn[kount]=ija[ii];          /* column # in global jacobian matrix */
	  aj_off_diag[kount] = a[ii];
	  kount=kount+1;
	}
    }

  DPRINTF(stderr, "Sorting nonzeros ...");
  piksr2(nn, jcolumn, irow, aj_off_diag);  /* arrange coefficient matrix columnwise,*/
                                            /* in ascending column number order */
  DPRINTF(stderr, "done\n");

  /*
   *  now calculate analytical and numerical jacobians at perturbed values
   *  check that the perturbed residuals are consistent with range possible
   *  for range of analytical jacobian values
   */
  for (j = 0; j < NumUnknowns[pg->imtrx]; j++)       /* loop over each column */
    {
      /*
       * Perturb one variable at a time
       */

     if ( ls != NULL && ls->Ignore_F_deps && idv[pg->imtrx][j][0] == FILL ) continue;

#ifdef FORWARD_DIFF_NUMJAC
      x_1[j] = x[j] + x_scale[idv[pg->imtrx][j][0]] * DELTA_UNKNOWN;
#else
      dx = x_scale[idv[pg->imtrx][j][0]] * FD_DELTA_UNKNOWN;
      x_1[j] = x[j] + dx;
#endif

      num_elems = 0;
      for(i = 0; i < ELEM_LIST_SIZE; i++)
	elem_list[i] = 0;
      for(i = 0; i < NumUnknowns[pg->imtrx]; i++)
	output_list[i] = FALSE;

      af->Assemble_Residual = TRUE;
      af->Assemble_LSA_Jacobian_Matrix = FALSE;
      af->Assemble_LSA_Mass_Matrix = FALSE;
#ifdef FORWARD_DIFF_NUMJAC
      af->Assemble_Jacobian = FALSE;
#else
      af->Assemble_Jacobian = TRUE;
#endif
      neg_elem_volume = FALSE;
      neg_lub_height = FALSE;
      zero_detJ = FALSE;

      my_node_num = idv[pg->imtrx][j][2];

      /* Which elements to fill?  We need every element that contains
       * this node, plus all of the elements connected to them, even if
       * they are only connected by a node (and not a side).
       *
       * First, we put all the elements containing this node into the
       * list.  It is not possible for repeated elements, here.
       */
      for(i = exo->node_elem_pntr[my_node_num];
	  i < exo->node_elem_pntr[my_node_num+1]; i++)
	{
	  my_elem_num = exo->node_elem_list[i];
	  elem_list[num_elems++] = my_elem_num;
	}

      /* Now we go through each element we have, then each node on
       * those elements, then every element containing those nodes.
       * Add those elements, if they have not already been added.
       */
      for(i = exo->node_elem_pntr[my_node_num];
	  i < exo->node_elem_pntr[my_node_num+1]; i++)
	{
	  my_elem_num = exo->node_elem_list[i];
	  for(k = exo->elem_node_pntr[my_elem_num];
	      k < exo->elem_node_pntr[my_elem_num+1]; k++)
	    {
	      node_num = exo->elem_node_list[k];
	      for(l = exo->node_elem_pntr[node_num];
		  l < exo->node_elem_pntr[node_num+1]; l++)
		{
		  elem_num = exo->node_elem_list[l];
		  if(elem_num == -1)
		    continue;
		  elem_already_listed = FALSE;
		  for(m = 0; m < num_elems; m++)
		    if(elem_list[m] == elem_num)
		      {
			elem_already_listed = TRUE;
			break;
		      }
		  if(!elem_already_listed)
		    elem_list[num_elems++] = elem_num;
		}
	    }
	}

      /* For which variables do we report the numerical vs. analytic
       * jacobians?  Only those that are actually contained in an
       * element that contains our unknown's node.  We need to search
       * for all of the unknowns on all of the nodes on all of these
       * elements (whew!).  We specifically SHOULDN'T compare
       * numerical and analytic jacobians for any nodes except these,
       * because all of those nodes are not fully populated (so that
       * the residuals will come out incorrect for comparison
       * purposes).
       */
      for (i = exo->node_elem_pntr[my_node_num];
	   i < exo->node_elem_pntr[my_node_num+1]; i++) {
	my_elem_num = exo->node_elem_list[i];
	load_ei(my_elem_num, exo, 0, pg->imtrx);
	for (k = exo->elem_node_pntr[my_elem_num];
	     k < exo->elem_node_pntr[my_elem_num+1]; k++) {
	  node_num = exo->elem_node_list[k];
	  node = Nodes[node_num];
	  nvs = node->Nodal_Vars_Info[pg->imtrx];
	  for (l = 0; l < nvs->Num_Var_Desc; l++) {
	    vd = nvs->Var_Desc_List[l];
	    for (m = 0; m < vd->Ndof; m++) {
	      index = node->First_Unknown[pg->imtrx] + nvs->Nodal_Offset[l] + m;
	      output_list[index] = TRUE;
	    }
	  }
	}
      }

      /* make compact list of Eqdof's that will be checked; put diagonal term first */
      dof_list[0] = j;
      num_dofs = 1;
      for (i=0; i<NumUnknowns[pg->imtrx]; i++)
        {
          if (i!=j && output_list[i])
            {
              dof_list[num_dofs++] = i;
            }
        }

      /* compute residual and Jacobian at perturbed soln */
      memset(a, 0, nnonzero*sizeof(dbl));
      memset(resid_vector_1, 0, NumUnknowns[pg->imtrx]*sizeof(dbl));

      if (pd_glob[0]->TimeIntegration != STEADY) {
	xdot[j] += (x_1[j] - x[j])  * (1.0 + 2 * theta) / delta_t;
      }

      if ( xfem != NULL )
        clear_xfem_contribution( ams->npu );

      for (i = 0; i < num_elems; i++) {
	zeroCA = -1;
	if (i == 0) zeroCA = 1;
	load_ei(elem_list[i], exo, 0, pg->imtrx);
	matrix_fill(ams, x_1, resid_vector_1,
		    x_old, x_older,  xdot, xdot_old, x_update,
		    &delta_t, &theta,
		    first_elem_side_BC_array,
		    &time_value, exo, dpi,
		    &elem_list[i], &num_total_nodes,
		    h_elem_avg, U_norm, NULL, zeroCA);
	if( neg_elem_volume ) break;
	if( neg_lub_height ) break;
	if( zero_detJ ) break;
      }

      if ( xfem != NULL )
        check_xfem_contribution( ams->npu, ams, resid_vector_1, x_1, exo );

      /*
       * Free memory allocated above
       */
      global_qp_storage_destroy();

#ifdef PARALLEL
      neg_elem_volume_global = FALSE;
      MPI_Allreduce(&neg_elem_volume, &neg_elem_volume_global, 1,
                     MPI_INT, MPI_LOR, MPI_COMM_WORLD);
      neg_elem_volume = neg_elem_volume_global;

      neg_lub_height_global = FALSE;
      MPI_Allreduce(&neg_lub_height, &neg_lub_height_global, 1,
                     MPI_INT, MPI_LOR, MPI_COMM_WORLD);
      neg_lub_height = neg_lub_height_global;

      zero_detJ_global = FALSE;
      MPI_Allreduce(&zero_detJ, &zero_detJ_global, 1,
                     MPI_INT, MPI_LOR, MPI_COMM_WORLD);
      zero_detJ = zero_detJ_global;
#endif

      if (neg_elem_volume) {
	DPRINTF(stderr, "neg_elem_volume triggered \n");
	exit(-1);
      }
      if (neg_lub_height) {
        DPRINTF(stderr, "neg_lub_height triggered \n");
        exit(-1);
      }
      if (zero_detJ) {
        DPRINTF(stderr, "zero_detJ triggered \n");
        exit(-1);
      }

#ifdef DEBUG_NUMJAC
      DPRINTF(stderr, "For j = % 2d, before scaling:\n", j);
      for (ii = 0; ii < num_dofs; ii++)
        {
	  i = dof_list[ii];
	  DPRINTF(stderr, "resid[% 2d] = %-10.4g\n", i, resid_vector_1[i]);
	}
#endif

      if (Debug_Flag == -2 || Debug_Flag == -3)
	{
          /* Scale to get rid of problems with
	   * penalty parameter.
	   */
#ifdef FORWARD_DIFF_NUMJAC
          vector_scaling(NumUnknowns[pg->imtrx], resid_vector_1, scale);
#else
          row_scaling(NumUnknowns[pg->imtrx], a, ija, resid_vector_1, scale);
#endif
	}

#ifdef DEBUG_NUMJAC
      DPRINTF(stderr, "For j = % 2d, after scaling:\n", j);
      for (ii = 0; ii < num_dofs; ii++)
        {
	  i = dof_list[ii];
	  DPRINTF(stderr, "resid[% 2d] = %-10.4g\n", i, resid_vector_1[i]);
	}
#endif

#ifdef FORWARD_DIFF_NUMJAC
      for (ii = 0; ii < num_dofs; ii++)
        {
	  i = dof_list[ii];
	  if(x[j] != 0.0)
	    nj[i] = (resid_vector_1[i] - resid_vector[i])/(x[j] * DELTA_UNKNOWN);
	  else
	    nj[i] = (resid_vector_1[i] - resid_vector[i])/DELTA_UNKNOWN;
#ifdef DEBUG_NUMJAC
	    DPRINTF(stderr,"x[%02d]=%-15g, resid_vector_1[%02d]=%-15g, resid_vector[%02d]=%-15g",
		    i,x[i],i,resid_vector_1[i],i,resid_vector[i]);
	    DPRINTF(stderr, " nj[%02d]=%-15g\n",i,nj[i]);
#endif
	}

      compute_numerical_jacobian_errors(aj_diag[j], nj[j], &nj_err, &nj_scaled_err);

      /* COMPARISON: analytical vs. numerical --- the diagonal element for column j */
      if(nj_err >= RESIDUAL_TOLERANCE ||
	 nj_scaled_err >= SCALED_RESIDUAL_TOLERANCE)
	{
	  DPRINTF(stderr, "Diag%.22s Var%.22s aj=%-10.4g nj=%-10.4g er=%9.4g rer=%9.4g x=%-10.4g r=%-10.4g\n"
		  ,resname[pg->imtrx][j],dofname[pg->imtrx][j], aj_diag[j], nj[j], nj_err,
		  nj_scaled_err, x_1[j], resid_vector_1[j]);
	}


      /* COMPARISON: analytical vs. numerical ---  the off-diagonal elements for column j */
      for (k=0; k<(ija[NumUnknowns[pg->imtrx]]-ija[0]); k++)
	{
	  if(jcolumn[k] == j)       /* match the column numbers */
	    {
	      for (ii = 0; ii < num_dofs; ii++)
                {
                  i = dof_list[ii];
		  if(irow[k] == i)      /* match the row numbers */
		    {
		      compute_numerical_jacobian_errors(aj_off_diag[k], nj[i], &nj_err, &nj_scaled_err);
		      /* MMH
		       * Added in a condition that the residual has to be
		       * "bigger than zero" (or 1e-6).  This could also
		       * be done as /aj_diag[j]...
		       */
		      if(nj_err >= RESIDUAL_TOLERANCE ||
			 nj_scaled_err >= SCALED_RESIDUAL_TOLERANCE)
			DPRINTF(stderr,
				"  Eq%.42s Var%.42s aj=%-10.4g nj=%-10.4g er=%9.4g rer=%9.4g x=%-10.4g r=%-10.4g\n",
				resname[pg->imtrx][irow[k]], dofname[pg->imtrx][jcolumn[k]],
				aj_off_diag[k], nj[i], nj_err,
				nj_scaled_err, x_1[j], resid_vector_1[j]);
		    }
		}
	    }
	}
#else
      /* BRACKETED JACOBIAN CHECKER */
      /* extract diagonal and off-diagonal elements from the coefficient matrix stored
         in sparse-storage format */
      memset(aj, 0, NumUnknowns[pg->imtrx]*sizeof(dbl));
      memset(aj_1, 0, NumUnknowns[pg->imtrx]*sizeof(dbl));
      for (ii = 0; ii < num_dofs; ii++)
        {
	  i = dof_list[ii];
          if (i == j)
            {
              aj[i] = aj_diag[j];
              aj_1[i] = a[i];
            }
          else
            {
              for (k=0; k<(ija[NumUnknowns[pg->imtrx]]-ija[0]); k++)
                {
                  if ((jcolumn[k] == j) && (irow[k] == i))
                    {
                      aj[i] = aj_off_diag[k];
                    }
                }
              for (k = ija[i]; k< ija[i+1]; k++)
                {
                  if (ija[k] == j)
                    {
                      aj_1[i] = a[k];
                    }
                }
            }
        }

      for (ii = 0; ii < num_dofs; ii++)
        {

          i = dof_list[ii];

          /* compute valid range for resid_vector_1[i] */
          if (dx*aj[i] > dx*aj_1[i])
            {
              delta_min = dx*aj_1[i];
              delta_max = dx*aj[i];
            }
          else
            {
              delta_min = dx*aj[i];
              delta_max = dx*aj_1[i];
            }

          /* attempt to calculate the magnitude of the roundoff error in the
           * residual calculation.  this is estimated to be MAX( J x_scale ) */
          resid_scale = fabs( resid_vector[i] );
          /*DRN-MAX BROKEN? resid_scale = MAX( resid_scale, fabs( a[i] * x_scale[idv[pg->imtrx][i][0]] ) );*/
	  if ( fabs( a[i] * x_scale[idv[pg->imtrx][i][0]] ) > resid_scale ) resid_scale = fabs( a[i] * x_scale[idv[pg->imtrx][i][0]] );
          for (k = ija[i]; k< ija[i+1]; k++)
            {
              var_j = idv[pg->imtrx][ija[k]][0];
              /*DRN-MAX BROKEN? resid_scale = MAX( resid_scale, fabs( a[k] * x_scale[var_j]) );*/
	      if ( fabs( a[k] * x_scale[var_j]) > resid_scale ) resid_scale = fabs( a[k] * x_scale[var_j]);
            }

          roundoff = 1.e-11;
          resid_min = resid_vector[i] + delta_min - roundoff*resid_scale;
          resid_max = resid_vector[i] + delta_max + roundoff*resid_scale;
          resid_diff = delta_max - delta_min;
          resid_diff += 2.*roundoff*resid_scale;

          if (resid_vector_1[i]<resid_min || resid_vector_1[i]>resid_max)
            {
              nj = (resid_vector_1[i] - resid_vector[i]) / dx;
              /* this scaled error examines the size of the deviation to
               * the width of the acceptance band.  The bigger the number
               * the more confidence you should be able to have that the
               * error is significant */
              /*DRN-MAX BROKEN? resid_error = MAX( resid_vector_1[i] - resid_max, resid_min - resid_vector_1[i] );*/
	      if ( resid_vector_1[i]<resid_min ) resid_error = resid_min - resid_vector_1[i];
	      else resid_error = resid_vector_1[i] - resid_max;

              /* resid_scaled_error = resid_error / resid_diff; */

              /* The following is a measure of the percentage of the acceptance band that is
               * due to changes in the jacobian from the unperturbed to perturbed
               * states.  The closer this is to 1, the more confidence you have
               * that the error is significant.  However, this will always be zero
               * for a constant sensitivity.  The more linear the relationship over
               * dx, the smaller this will be. A value of 0 means that one would expect
               * the numerical jacobian to agree with the analytical ones (which in this
               * case are the same) to within the roundoff error.
               * If this value is small and resid_scaled_error is small (meaning order 1)
               * than you may be safe ignoring this entry.  If this is close to 1, you
               * can be fairly confident that this is a signficant error for any value
               * of resid_scaled_error.
               */
              delta_aj_percentage = (delta_max - delta_min) / resid_diff;

              /* Here's an attempt to combine these ideas to give a measure of
               * confidence that we are flagging a true error.  This is the error
               * relative to the expected roundoff error.  I am a little hesitant
               * to use this because of the inherent uncertainty in the scale of the
               * residual error.
               */
              confidence = resid_error /
                           (2.*roundoff*resid_scale);

              if (i==j)
                {
                  DPRINTF(stderr,
                    "Diag%32.32s Var%32.32s x=%-10.4g dx=%-10.4g aj=%-10.4g nj=%-10.4g aj_1=%-10.4g d_aj=%-10.4g conf=%-10.4g\n",
	                  resname[pg->imtrx][i], dofname[pg->imtrx][j], x[j], dx,
			              aj[i], nj, aj_1[i],
                    delta_aj_percentage, confidence );
                }
              else
                {
		              DPRINTF(stderr,
			              "Eq%-32.32s Var%-32.32s x=%-10.2g dx=%-10.2g aj=%-10.4g nj=%-10.4g aj_1=%-10.2g d_aj=%-10.2g conf=%-10.2g\n",
			              resname[pg->imtrx][i], dofname[pg->imtrx][j], x[j], dx,
			              aj[i], nj, aj_1[i],
                    delta_aj_percentage, confidence );
                }

              /* print disclaimer for entries that have much lower confidence levels */
              if (x[j] == 0.)
                {
                  switch (idv[pg->imtrx][j][0])
                    {
                    case MESH_DISPLACEMENT1:
                    case MESH_DISPLACEMENT2:
                    case MESH_DISPLACEMENT3:
                    case SOLID_DISPLACEMENT1:
                    case SOLID_DISPLACEMENT2:
                    case SOLID_DISPLACEMENT3:

                    /* DPRINTF(stderr, "  NOTE: Jacobian errors associated with displacements "
                                      "are less reliable with undeformed initial conditions.\n");
*/
                      break;

                    default:
                      break;
                    }
                }

              /* list all BC's applied at this node and highlight ones with
               * desired sensitivity */
              my_node_num = idv[pg->imtrx][i][2];

              for (k = exo->node_elem_pntr[my_node_num];
                   k < exo->node_elem_pntr[my_node_num+1]; k++)
                {
                  my_elem_num = exo->node_elem_list[k];
                  load_ei(my_elem_num, exo, 0, pg->imtrx);

                  if (first_elem_side_BC_array[my_elem_num] != NULL)
                    {

                      elem_side_bc = first_elem_side_BC_array[my_elem_num];
                      /***************************************************************************
                       *  begining of do while construct which loops over the sides of this
                       *  element that have boundary conditions applied on to them.
                       ***************************************************************************/
                      do
                        {
                          for (ibc = 0; (bc_input_id = (int) elem_side_bc->BC_input_id[ibc]) != -1; ibc++)
                            {
                              var_i = idv[pg->imtrx][i][0];
                              var_j = idv[pg->imtrx][j][0];
                              I = idv[pg->imtrx][i][2];

                              if (in_list( I, 0, elem_side_bc->num_nodes_on_side, elem_side_bc->local_node_id ))
                                {
                                  eqn = var_i;
                                  if (BC_Types[bc_input_id].desc->vector &&
				      (eqn == VELOCITY2 || eqn == VELOCITY3) ) eqn = VELOCITY1;
#ifdef DEBUG_NUMJAC
                                  if (BC_Types[bc_input_id].desc->equation == eqn &&
				      BC_Types[bc_input_id].desc->sens[var_j] )
                                    {
                                      DPRINTF(stderr, "  >>> ");
                                    }
                                  else
                                    {
                                      DPRINTF(stderr, "      ");
                                    }

                                  DPRINTF(stderr, "%s on %sID=%d\n",
                                          BC_Types[bc_input_id].desc->name1,
			                  BC_Types[bc_input_id].Set_Type,
			                  BC_Types[bc_input_id].BC_ID);
#endif
                                }
                            }
                        } while ((elem_side_bc = elem_side_bc->next_side_bc) != NULL);
                    } /* END if (First_Elem_Side_BC_Array[my_elem_num] != NULL) */
                }
            }

          /* check Inter_Mask[pg->imtrx] for missing entries */
          var_i = idv[pg->imtrx][i][0];
          var_j = idv[pg->imtrx][j][0];
          if (!Inter_Mask_save[var_i][var_j])
            {
              /* check to make sure no dependence appears in analytical jacobian */
              if ((aj[i] != 0.) || (aj_1[i] != 0.))
                {
                  DPRINTF(stderr,
                    "Potential dependency error: Inter_Mask[pg->imtrx][Eq%.32s][Var%.32s]=0, but aj=%-10.4g aj_1=%-10.4g\n",
	                  resname[pg->imtrx][i], dofname[pg->imtrx][j], aj[i], aj_1[i] );
                }

              /* check to make sure no dependence appears in numerical jacobian */
              resid_min = resid_vector[i] - roundoff*resid_scale;
              resid_max = resid_vector[i] + roundoff*resid_scale;
              if (resid_vector_1[i]<resid_min || resid_vector_1[i]>resid_max)
                {
                  nj = (resid_vector_1[i] - resid_vector[i]) / dx;
                  DPRINTF(stderr,
                    "Potential dependency error: Inter_Mask[pg->imtrx][Eq%.32s][Var%.32s]=0, but nj=%-10.4g\n",
	                  resname[pg->imtrx][i], dofname[pg->imtrx][j], nj );
                }
            }
        }
#endif

      /*
       * return solution vector to its original state
       */
      if (pd_glob[0]->TimeIntegration != STEADY) {
	xdot[j] -= (x_1[j] - x[j])  * (1.0 + 2 * theta) / delta_t;
      }
      x_1[j] = x[j];
    }                          /* End of for (j=0; j<NumUnknowns[pg->imtrx]; j++) */

  /* free arrays to hold jacobian and vector values */
  safe_free( (void *) irow) ;
  safe_free( (void *) jcolumn) ;
  safe_free( (void *) nelem) ;
  safe_free( (void *) aj_diag) ;
  safe_free( (void *) aj_off_diag) ;
  safe_free( (void *) resid_vector_1) ;
  safe_free( (void *) x_1) ;
  safe_free( (void *) scale) ;
  safe_free( (void *) output_list);
  safe_free( (void *) dof_list);
  safe_free( (void *) elem_list);
#ifdef FORWARD_DIFF_NUMJAC
  safe_free( (void *) nj) ;
#else
  safe_free( (void *) aj) ;
  safe_free( (void *) aj_1) ;
#endif

}                             /*   End of function numerical_jacobian  */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int
intcompare(const void *left, const void *right)

    /*
     *
     */
{
  register int left_a_val, right_a_val;

  left_a_val = ((const data_t *)left)->a_val;
  right_a_val = ((const data_t *)right)->a_val;
  if (left_a_val > right_a_val)
    return 1;
  else if (left_a_val < right_a_val)
    return -1;
  else
    return 0;
}


/* function sorts arrays arr[], brr[], crr[] in ascending values of arr[] */
static void
piksr2(const int n, int arr[], int brr[], dbl crr[])
{
  int i;
  data_t *vals;

  vals = (data_t *)array_alloc(1, n, sizeof(data_t));
  for(i = 0; i < n; i++)
    {
      vals[i].a_val = arr[i];
      vals[i].b_val = brr[i];
      vals[i].c_val = crr[i];
    }
  qsort((data_t *)vals, n, sizeof(data_t), intcompare);
  for(i = 0; i < n; i++)
    {
      arr[i] = vals[i].a_val;
      brr[i] = vals[i].b_val;
      crr[i] = vals[i].c_val;
    }
  safe_free((void *)vals);
}       /* End of function piksr2 */


/* Do something intelligent with our raw errors... For example, if the
 * diagonal value is 1, and aj = 0.0, and nj = 1.0e-20, then this
 * probably isn't an actual error!
 */
#ifdef FORWARD_DIFF_NUMJAC
static void
compute_numerical_jacobian_errors(const dbl exact, /* analytic value */
				  const dbl approx, /* numerical value */
				  dbl *nj_err, /* absolute error */
				  dbl *nj_scaled_err) /* scaled error */
{
  dbl fexact, fapprox;

  fexact = fabs(exact);
  fapprox = fabs(approx);
  *nj_err = fabs(exact-approx);

  /* Case 1: exact is zero, approx is numerically zero. */
  if(fexact == 0.0 && fapprox < SCALED_RESIDUAL_TOLERANCE_CUTOFF)
    {
      *nj_scaled_err = 0.0;
      return;
    }
  /* Case 2: approx is zero, exact is numerically zero. */
  if(fapprox == 0.0 && fexact < SCALED_RESIDUAL_TOLERANCE_CUTOFF)
    {
      *nj_scaled_err = 0.0;
      return;
    }
  /* Case 3: exact is zero, approx is nonzero. */
  if(fexact == 0.0)
    {
      *nj_scaled_err = 1.0e+10;
      return;
      /* Case 3a: cannot ignore size of absolute error */
      /*
      if(*nj_err >= MIXED_RESIDUAL_TOLERANCE)
	*nj_scaled_err = 1.0e+10;
	*/
      /* Case 3b: The absolute error is so small that we ignore it. */
      /*
      else
	*nj_scaled_err = 0.0;
      return;
      */
    }
  /* Case 4: exact is nonzero, report scaled error. */
  *nj_scaled_err = *nj_err / fexact;
  /* Case 5: have a scaled error that causes reporting, but the
   * absolute error is so small that we don't worry about it.
   */
  if(*nj_scaled_err >= SCALED_RESIDUAL_TOLERANCE &&
     *nj_err < MIXED_RESIDUAL_TOLERANCE)
    *nj_scaled_err = 0.0;
  /* Case 6: have an absolute error that causes reporting, but the
   * scaled error is so small that we don't worry about it.
   */
  if(*nj_err >= RESIDUAL_TOLERANCE &&
     *nj_scaled_err < MIXED_SCALED_RESIDUAL_TOLERANCE)
    *nj_err = 0.0;
}
#endif

#ifndef COUPLED_FILL
void
numerical_jacobian_fill(int ijaf[],	/* fill Vector of integer pointers into a matrix */
			double afill[], /* Vector of non-zero entries in the
					 * coefficient matrix */
			double xf[],	/* fill Solution vector for the current processor */
			double rf[],    /* Residual vector for the current
					 * processor */
			double delta_t, /* time step size */
		        double theta,   /* parameter to vary time integration
                                         * from explicit (theta = 1) to
				         *      implicit (theta = 0) */
			double x[],     /* Value current big solution vector holding everything*/
			double x_old[], /* Value of the old solution vector */
			double xdot[],  /* Value of xdot predicted for new solution */


			int Debug_Flag, /* flag for calculating numerical jacobian
					   -3 == calc num jac w/  rescaling */
			int node_to_fill[], /* this is a map from the */
			Exo_DB *exo,	    /* ptr to whole fe mesh */
			Dpi *dpi) /* ptr to parallel info */

/******************************************************************************
  This function compares the analytical jacobian entries calculated in matrix_fill
                         the numerical ones approximated by central difference method.

  Author:          K. S. Chen (1511) (based on an earlier version by P. R. Schunk).
  Date:            January 19, 1994

******************************************************************************/

{
  int i, j, k, ii, nn, kount, nnonzero, num_total_nodes;

  extern struct elem_side_bc_struct **First_Elem_Side_BC_Array;

  const double epsilon = 1.e-07;


  double epsilon1=1.e-3;
  dbl *aj_diag, *aj_off_diag, *nj, *resid_vector_old, *x_save, *scale;
  int *irow, *jcolumn, *nelem;


  fprintf(stderr,"\n Starting Numerical Jacobian Checker for Fill equation\n");

  nnonzero	  = fill_zeros+1;
  nn		  = ijaf[num_fill_unknowns]-ijaf[0]; /* total number of diagonal entries a[] */
  num_total_nodes = Num_Internal_Nodes + Num_Border_Nodes;

/* allocate arrays to hold jacobian and vector values */
  irow		   = (int *) array_alloc(1, nnonzero, sizeof(int));
  jcolumn	   = (int *) array_alloc(1, nnonzero, sizeof(int));
  nelem		   = (int *) array_alloc(1, nnonzero, sizeof(int));
  aj_diag	   = (double *) array_alloc(1, num_fill_unknowns, sizeof(double));
  aj_off_diag	   = (double *) array_alloc(1, nnonzero, sizeof(double));
  nj		   = (double *) array_alloc(1, nnonzero, sizeof(double));
  resid_vector_old = (double *) array_alloc(1, num_fill_unknowns, sizeof(double));
  x_save	   = (double *) array_alloc(1, num_fill_unknowns, sizeof(double));
  scale		   = (double *) array_alloc(1, num_fill_unknowns, sizeof(double));

  if (nj == NULL || scale == NULL) EH(-1, "No room for numerical jacobian arrays");

  if (Debug_Flag == -2) epsilon1 = 1.e-6;

  /* initialization */
  memset(aj_off_diag, 0, nnonzero*sizeof(dbl)); /* off-diagonal analytical jacobian elements */
  memset(nj, 0, nnonzero*sizeof(dbl));          /* numerical jacobian elements */
  memset(aj_diag, 0, num_fill_unknowns*sizeof(dbl));  /* diagonal analytical jacobian elements */


  /* first calculate the residual vector corresponding to the solution vector read in
     the initial guess file; also calculate the analytical jacobian entries */

  af->Assemble_Residual		   = TRUE;
  af->Assemble_Jacobian		   = TRUE;
  af->Assemble_LSA_Jacobian_Matrix = FALSE;
  af->Assemble_LSA_Mass_Matrix	   = FALSE;

  (void) fill_matrix(afill, ijaf, rf, xf, x, x_old,
		     xdot, delta_t, theta, ADVECT, node_to_fill,
		     First_Elem_Side_BC_Array, exo, dpi);

  if (Debug_Flag == -2) {
    /* Scale matrix first to get rid of problems with
       penalt parameter */
    row_sum_scale_MSR(num_fill_unknowns, afill, ijaf, rf, scale);
  }

  /* save solution vector and residual vector before numerical jacobian calculations */

  dcopy1( num_fill_unknowns, xf, x_save);
  dcopy1( num_fill_unknowns, rf, resid_vector_old);

  /* extract diagonal and off-diagonal elements from the coefficient matrix stored
     in sparse-storage format */

  dcopy1(num_fill_unknowns, afill, aj_diag);  /* diagonal elements */

  kount=0;                              /* off-diagonal elements */
  for (i=0; i<num_fill_unknowns; i++)
    {
      nelem[i] = ijaf[i+1] - ijaf[i];
      for (k=0; k<nelem[i]; k++)
	{
	  irow[kount]=i;                   /* row # in global jacobian matrix */
	  ii = kount + num_fill_unknowns + 1;
	  jcolumn[kount]=ijaf[ii];          /* column # in global jacobian matrix */
	  aj_off_diag[kount] = afill[ii];
	  kount=kount+1;
	}
    }

  piksr2(nn, jcolumn, irow, aj_off_diag);  /* arrange coefficient matrix columnwise,*/
  /* in ascending column number order */


  /* calculate numerical jacobian entries columnwise and then compare them with the
     analytical jacobian entries */

  for (j=0; j<num_fill_unknowns; j++)       /* loop over each column */
    {

      xf[j] = x_save[j] + epsilon;     /* perturb one variable at a time */
      /* let big vector know of the change for load_fv */

      put_fill_vector(num_total_nodes, x, xf, node_to_fill);

      for (i=0; i<num_fill_unknowns; i++)
	rf[i] = 0.0;         /* zero residual vector before its calculation */

      af->Assemble_Residual = TRUE;
      af->Assemble_Jacobian = FALSE;
      af->Assemble_LSA_Jacobian_Matrix = FALSE;
      af->Assemble_LSA_Mass_Matrix = FALSE;

      (void) fill_matrix(afill, ijaf, rf, xf, x, x_old,
			 xdot, delta_t, theta, ADVECT,
                         node_to_fill, First_Elem_Side_BC_Array, exo, dpi);

      if (Debug_Flag == -2) {
	/* Scale matrix first to get rid of problems with
	   penalt parameter */
	row_scaling(num_fill_unknowns, afill, ijaf, rf, scale);
      }

      for (i=0; i<num_fill_unknowns; i++)  /* cal numerical jacobian vector for column j */
	nj[i] = (rf[i] - resid_vector_old[i])/epsilon;


      /* COMPARISON: analytical vs. numerical --- the diagonal element for column j */

      if(ABS(aj_diag[j] - nj[j]) > epsilon1)
	fprintf(stderr, " aj=%-10.4g nj=%-10.4g resid=%-12.5g at unknown j = %d\n",
		aj_diag[j], nj[j], rf[j], j);


      /* COMPARISON: analytical vs. numerical ---  the off-diagonal elements for column j */

      for (k=0; k<(ijaf[num_fill_unknowns]-ijaf[0]); k++)
	{
	  if(jcolumn[k] == j)       /* match the column numbers */
	    {
	      for (i=0; i<num_fill_unknowns; i++)
		{
		  if(i == irow[k])      /* match the row numbers */
		    {

		      if(ABS(aj_off_diag[k]-nj[i]) > epsilon1 && aj_off_diag[k] != 1.0e+06)
			fprintf(stderr," aj=%-10.4g nj=%-10.4g at unknown j = %d row i = %d\n",
				aj_off_diag[k], nj[i], j, i);
		    }
		}
	    }
	}

      xf[j] = xf[j] - epsilon;  /* return solution vector to its original state */

    }                          /* End of for (j=0; j<num_fill_unknowns; j++) */

  /* free arrays to hold jacobian and vector values */
  safe_free( (void *) irow) ;
  safe_free( (void *) jcolumn) ;
  safe_free( (void *) nelem) ;
  safe_free( (void *) aj_diag) ;
  safe_free( (void *) aj_off_diag) ;
  safe_free( (void *) nj) ;
  safe_free( (void *) resid_vector_old) ;
  safe_free( (void *) x_save) ;
  safe_free( (void *) scale) ;

}                             /*   End of function numerical_jacobian_fill  */
#endif /* not COUPLED_FILL */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double
calc_numerical_delta(double base_value)

    /*************************************************************************
     *
     *  calc_numerical_delta():
     *
     *   This routine calculates a delta given a variable. It doesn't do
     *   anything elegant, such as checking for a "usual value" for the
     *   variable. Therefore, problems may arise.
     *
     *   The return value is guarranteed to be positive definate.
     *
     *************************************************************************/
{
  return (fabs(base_value) * DELTA_UNKNOWN + 1.0E-10);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int originalChoice_Jacobian = -23;
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void
AF_assemble_Residual_Only(void)

    /*************************************************************************
     *
     *  AF_assemble_Residual_Only
     *
     *  This routine sets the action flag to just calculate the residual
     *  only
     *
     *************************************************************************/
{
  originalChoice_Jacobian = af->Assemble_Jacobian;
  if (af->Assemble_Jacobian) {
     af->Assemble_Jacobian = FALSE;
  }
  if (! af->Assemble_Residual) af->Assemble_Residual = TRUE;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void
AF_restore_Jacobian_Flag(void)

    /*************************************************************************
     *
     *  AF_restore_Jacobian_Flag
     *
     *
     *************************************************************************/
{
  if (originalChoice_Jacobian == -23) {
    EH(-1,"restore_Jacobian_Action_Flag ERROR: no original choice");
  }
  af->Assemble_Jacobian = originalChoice_Jacobian;
  originalChoice_Jacobian = -23;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
