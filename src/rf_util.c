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

/* Needed to declare POSIX function drand48 */
#include "rf_solve.h"
#define _XOPEN_SOURCE

#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bc_contact.h"
#include "dp_comm.h"
#include "dp_types.h"
#include "dpi.h"
#include "el_elm.h"
#include "el_elm_info.h"
#include "el_geom.h"
#include "exo_struct.h"
#include "exodusII.h"
#include "md_timer.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_augc_util.h"
#include "mm_eh.h"
#include "mm_fill_ls.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_util.h"
#include "mm_input.h"
#include "mm_more_utils.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "mm_species.h"
#include "mm_std_models.h"
#include "mm_unknown_map.h"
#include "rd_mesh.h"
#include "rd_pixel_image.h"
#include "rf_allo.h"
#include "rf_bc_const.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_io.h"
#include "rf_io_const.h"
#include "rf_mp.h"
#include "rf_node_const.h"
#include "rf_solver.h"
#include "rf_util.h"
#include "rf_vars_const.h"
#include "std.h"
/************ R O U T I N E S   I N   T H I S   F I L E  **********************

       NAME            		TYPE        		CALL BY
----------------------------------------------------------------------
        log2				int
        sort2_int_double		void
        time_step_control		double
        path_step_control               double
        find_max			int
        find_min			int
        find_inter			int
        init_vec			void
        fill_dvec_rand			void
        random_vec	 	static void		init_vec()
        exchange_pointers		void
        get_remap_index			void
        remap_int_array			void
        remap_double_array		void
        print_global_vec		void
        rd_globals_from_exoII    static int

******************************************************************************/

/******************** PROTOTYPES FOR STATIC FUNCTIONS ************************/
static void read_initial_guess /* rf_util.c                                 */
    (double[],                 /* u                                         */
     const int,                /* np                                        */
     double[],                 /* uAC                                       */
     const int);               /* nAC                                       */

static void inject_nodal_vec(double[],        /* sol_vec - full dof vector for this proc   */
                             const int,       /* var_no - VELOCITY1, etc.                  */
                             const int,       /* k - species index                         */
                             const int,       /* idof - dof #                              */
                             const int,       /* matID - material index to scatter to      */
                             const double[]); /* nodal_vec - condensed node based vector   */

static void inject_elem_vec(double[],          /* sol_vec - full dof vector for this proc   */
                            const int,         /* var_no - VELOCITY1, etc.                  */
                            const int,         /* k - species index                         */
                            const int,         /* idof - dof #                              */
                            const int,         /* matID - material index to scatter to      */
                            const double[],    /* nodal_vec - condensed node based vector   */
                            const Exo_DB *exo, /* exodus database */
                            const int num_elems_blk); /* number of elements in the block. */

static void init_structural_shell_coord(double[]); /* u[] - solution vector */

static void init_shell_normal_unknowns(double[],        /* u[] - solution vector */
                                       const Exo_DB *); /* Exodus database */

/*****************************************************************************/
/*****************************************************************************/

#if 0

/* scatter_double_vector() -- akin to BLAS DCOPY,  y[index[i]] <- x[i]
 * 
 * Uses an additional array, index[], for the indeces of y getting the
 * sequential pieces of x[]. User is responsible for allocating x[], y[], 
 * index[] and filling index[] with meaningful values.
 *
 * Typically, y[] will be any size, while x[] and index[] must accomodate
 * indeces from 0 to n-1.
 *
 * Created: 1997/10/15 14:13 MDT pasacki@sandia.gov
 *
 * Revised:
 */

void
scatter_double_vector(const int    n, 
		      const int    index[],
		      double y[], 
		      const double x[])
{
  int i;

  for ( i=0; i<n; i++)
    {
      y[index[i]] = x[i];
    }

  return;
}


/* gather_double_vector() - akin to BLAS DCOPY,  y[index[i]] -> x[i]
 * 
 * Uses an additional array, index[], for the indeces of y getting the
 * sequential pieces of x[]. User is responsible for allocating x[], y[], 
 * index[] and filling index[] with meaningful values.
 *
 * Typically, y[] will be any size, while x[] and index[] must accomodate
 * indeces from 0 to n-1.
 *
 * Created: 1997/10/15 14:18 MDT pasacki@sandia.gov
 *
 * Revised:
 */

void
gather_double_vector(const int n, 
		     const int index[],
		     const double y[], 
		     double x[])
{
  int i;

  for ( i=0; i<n; i++)
    {
      x[i] = y[index[i]];
    }

  return;
}
#endif

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void sort2_int_int(const int n, int ra[], int rb[])

/*************************************************************************
 *       Numerical Recipies in C source code
 *       modified to have first argument an integer array (JS)
 *
 *       Sorts the array ra[1,..,n] in ascending numerical order
 *       using the heapsort algorithm, while making the corresponding
 *       rearrangement of the array rb[1,..,n].
 *       NOTE: Usually, rb is an index array (rb[i] = i), and thus
 *             rb is used to find the original index of the
 *             elements of ra[].
 *
 *       NOTE: In the usual brain dead Numerical Recipies in C
 *             way, ra[] and rb[] are expected to have base address
 *             indecies of 1. Therefore, you have to decrement them by
 *             one in the argument list, when calling this routine
 *             from a normal C routine.
 **************************************************************************/
{
  int l, j, ir, i, rra, rrb;
  /*
   *  No need to sort if one or fewer items.
   */
  if (n <= 1)
    return;
  l = (n >> 1) + 1;
  ir = n;
  for (;;) {
    if (l > 1) {
      rra = ra[--l];
      rrb = rb[l];
    } else {
      rra = ra[ir];
      rrb = rb[ir];
      ra[ir] = ra[1];
      rb[ir] = rb[1];
      if (--ir == 1) {
        ra[1] = rra;
        rb[1] = rrb;
        return;
      }
    }
    i = l;
    j = l << 1;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j + 1])
        ++j;
      if (rra < ra[j]) {
        ra[i] = ra[j];
        rb[i] = rb[j];
        j += (i = j);
      } else {
        j = ir + 1;
      }
    }
    ra[i] = rra;
    rb[i] = rrb;
  }
}

#ifndef COUPLED_FILL
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void put_fill_vector(const int N, double x[], const double fill_vector[], const int node_to_fill[])

/************************************************************************
 *
 * put_fill_vector:
 *
 *     Insert the value of the FILL variable type unknowns corresponding
 *     to the non-specific material case or the first material back into
 *     the main solution vector.
 *
 *  Input
 *     N = Number of nodes (internal plus boundary and maybe external
 *     fill_vector[] = FILL variable type solution vector
 *     node_to_fill[] = index of the first unknown at a node into the
 *                      FILL variable type solution vector.
 *  Output
 *     x[] = Main solution vector
 *************************************************************************/
{
  int i;     /* i count from 0 to total number of nodes   */
  int ie;    /* ie is the position of the variables of
                interest in the x and xpred vectors       */
  int mn;    /* Material number                           */
  int nvdof; /* variables to keep Index_Solution happy    */
  int ki;    /* counter from 0 to the number of dofs      */

  for (i = 0; i < N; i++) {
    nvdof = Dolphin[pg->imtrx][i][R_FILL]; /* Number of FILL dofs at this node. */
    for (ki = 0; ki < nvdof; ki++) {
      ie = Index_Solution(i, R_FILL, 0, ki, -1, pg->imtrx);
      if (ie != -1) {
        if (ie > NumUnknowns[pg->imtrx]) {
          GOMA_EH(ie, "put_fill_vector");
        } else {
          x[ie] = fill_vector[node_to_fill[i] + ki];
        }
      } else if (ie == -1) {
        mn = first_matID_at_node(i);
        ie = Index_Solution(i, R_FILL, 0, ki, mn, pg->imtrx);
        if (ie != -1) {
          x[ie] = fill_vector[node_to_fill[i] + ki];
        } else {
          GOMA_EH(ie, "put_fill_vector");
        }
      }
    }
  }
} /* END of routine put_fill_vector */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void get_fill_vector(const int N,
                     const double x[],
                     double fill_vector[],
                     const int node_to_fill[]) {
  int i;  /* i count from 0 to total number of nodes   */
  int ie; /* ie is the position of the variables of
             interest in the x and xpred vectors       */
  int mn;
  const int ktype = 0;
  int nvdof; /* Number of R_FILL degrees of freedom at a node */
  int ki;    /* counter from 0 to number of dofs          */

  for (i = 0; i < N; i++) {
    nvdof = Dolphin[pg->imtrx][i][R_FILL];
    for (ki = 0; ki < nvdof; ki++) {
      ie = Index_Solution(i, R_FILL, ktype, ki, -1, pg->imtrx);
      if (ie != -1) {
        fill_vector[node_to_fill[i] + ki] = x[ie];
      } else {
        mn = first_matID_at_node(i);
        ie = Index_Solution(i, R_FILL, ktype, ki, mn, pg->imtrx);
        if (ie != -1) {
          fill_vector[node_to_fill[i] + ki] = x[ie];
        }
      }
    }
  }
} /* END of routine get_fill_vector */
#endif /* not COUPLED_FILL */
/*****************************************************************************/

/* countmap_vardofs() -- fill [nodeindex]->dof map for a given var, return count
 *
 * This routine is a trivial generalization of RRR's count_fill_unknowns().
 * It may have utility for certain postprocessing applications as well as its
 * intended use of counting up the number of fill equation unknowns and making
 * a map from the node index into the fill equation dof number.
 *
 * Notes: map must point to an array of size num_nodes.
 *          Note: Dofs for MASS_FRACTION unknowns with multiple
 *                subvariables are counted as one degree of freedom here,
 *                just as in the original Dolphin array.
 *
 * Created: 1997/10/15 10:53 MDT pasacki@sandia.gov
 *
 * Revised:
 */

int countmap_vardofs(const int varType, const int num_nodes, int *map) {
  int node, nun, count = 0;
  NODAL_VARS_STRUCT *nv;
  if (varType < 0 || varType > MAX_VARIABLE_TYPES - 1) {
    GOMA_EH(-1, "Attempt to count a bogus variable.");
  }
  for (node = 0; node < num_nodes; node++) {
    nv = Nodes[node]->Nodal_Vars_Info[pg->imtrx];
    nun = get_nv_ndofs_modMF(nv, varType);
    if (nun > 0) {
      map[node] = count;
      count += nun;
    } else {
      map[node] = -1;
    }
  }
  return (count);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
#if 0
void
count_fill_unknowns ( int N,
                      int node_to_fill[] )

{
/* LOCAL VARIABLES */
  int i;               /* i count from 0 to total number of nodes   */

  num_fill_unknowns = 0;
  for (i = 0; i < N; i++)
    {
      if ( Dolphin[pg_imtrx][i][FILL] )
	{
	  node_to_fill[i] = num_fill_unknowns;
	  num_fill_unknowns += Dolphin[pg->imtrx][i][FILL];
	}
    }
} /* END of routine count_fill_unknowns */
/*****************************************************************************/



void
filter_fill ( const int N,	/* Filter the fill vector and     */
              double fill_vector[] )       /* clip it (limits are 0. and 1.) */

{
  int i;             /* i count from 0 to number of total fill unknowns  */

  for (i = 0; i < N; i++)
    {
      if(fill_vector[i] >= 1.)
	{
	  fill_vector[i] = 1.;
	}
      else if(fill_vector[i] <= 0.)
	{
	  fill_vector[i] = 0.;
	}
    }
} /* END of routine filter_fill */
#endif

/*****************************************************************************/

/* N - number of nodes x - solution vector at n calculated from
 *  implicit corrector filter_species_material_number is the material
 *  number of the first material containing a hyrdodynamic suspension
 *
 *  Censor the concentration field for SUSPENSION mass flux model.
 *  Notion is to keep the concentrations in the range 0< C < Max_packing.
 *
 *  And while you're at it, censor the shear rate field to be positive.
 */

int filter_conc(const int N, /* number of nodes */
                dbl x[],     /* solution vector */
                const int filter_species_material_number,
                /* 1st material number w/ suspension */
                const dbl cmin, /* lower cutoff */
                const dbl cmax) /* upper cutoff */
{
  const double minimum_shear_rate = 0.0;
  int i, ie;

  if (pd->e[pg->imtrx][R_MASS]) {
    for (i = 0; i < N; i++) {
      if (Dolphin[pg->imtrx][i][R_MASS]) {
        ie = Index_Solution(i, R_MASS, filter_species_material_number, 0, -1, pg->imtrx);
        if (ie != -1) {
          x[ie] = x[ie] > cmin ? x[ie] : cmin;
          if (cmax != 0.0) {
            x[ie] = x[ie] < cmax ? x[ie] : cmax;
          }
        }
      }
    }
  }
  if (pd->e[pg->imtrx][R_SHEAR_RATE]) {
    for (i = 0; i < N; i++) {
      if (Dolphin[pg->imtrx][i][R_SHEAR_RATE]) {
        ie = Index_Solution(i, R_SHEAR_RATE, 0, 0, -1, pg->imtrx);
        if (ie != -1) {
          if (x[ie] < minimum_shear_rate) {
            x[ie] = minimum_shear_rate;
          }
        }
      }
    }
  }
  return (1);
} /* END of routine filter_conc */
/***************************************************************************/

double time_step_control(const double delta_t,
                         const double delta_t_old,
                         const int const_delta_t,
                         const double x[],
                         const double x_pred[],
                         const double x_old[],
                         const double x_AC[],
                         const double x_AC_pred[],
                         const double eps,
                         int *success_dt,
                         const int use_var_norm[])

/**********************************************************************
 *
 * time_step_control()
 *
 *
 * Routine calculates the recommended size of the next time step from
 * the norm of the difference between the predicted and corrected
 * solution. If delta < 0 is input deck, the user desire for delta_t
 * to be constant and it is not altered by this routine.
 *
 * Note, mesh unknowns and pressure unknowns are often pseudo-equilib
 * variables. Thus, they are frequently left out of the calculation
 * of the time-step truncation error calculation.
 *
 * Input
 * ------
 *
 * N               - number of nodes
 * delta_t         - time step size for the current time step (
 *                   i.e., n)
 * delta_t_old     - time step size for the previous time step (n -1)
 * x[]             - solution vector at n calculated from implicit
 *                   corrector
 * x_pred[]        - solution vector at n calculated from explicit
 *                   predictor
 * x_AC            - extra unknown vector at n from implicit corrector
 * x_AC_pred       - extra unknown vector at n from explicit prediction
 * alpha           - time step control parameter
 * beta            - time step control parameter
 * eps              - time step control parameter
 * use_var_norm    - this tells us whether or not to use a variable
 *                   for the norm calculations. 0 means don't use ...
 *                   1 means use. The choice is made in the input deck
 *
 * Output
 *---------
 * success_dt      - 1 if the time step was successful and
 *                   0 if the time step was unsuccessful
 *
 * Return
 *--------
 * delta_t         - the recommended time step for the next iteration
 *                   whether or not the current iteration was
 *                   considered a success.
 *
 * HKM -> Species appear to be handled incorrectly. The effect of
 *        smaller concentration species may be completely neglected
 *        under the current system. Thus, we need to introduce an
 *        ATOL RTOL formalism into the routine below, and collect
 *        contributions according to variable index number rather
 *        than equation number. This should make the routine faster
 *        as well.
 ***********************************************************************/
{
  int i;                          /* i count from 0 to total number of nodes   */
  int num_unknowns;               /* num_unknowns is the actual number of
                                     variables used in the norm calculation    */
  int ncp[MAX_VARIABLE_TYPES];    /* number of ea var contributing to norm */
  double max[MAX_VARIABLE_TYPES]; /* L_infinity norm for variable type [] */
  int eqn;                        /* eqn is the equation being searched for    */
  VARIABLE_DESCRIPTION_STRUCT *vd;
  const double alpha = TIME_STEP_ALPHA;
  const double beta = TIME_STEP_BETA;
  double delta_t_new = 0.0;
  double abs_eps = fabs(eps);
  double Err_norm;
  double e_d, e_v, e_T, e_y, e_P, e_S, e_V, e_AC, e_qs;
  double e_shk, e_sht, e_shd, e_shu, e_F, e_ap, e_extv, e_sh_lub, e_int, e_rheo;
  double scaling;
  double ecp[MAX_VARIABLE_TYPES]; /* error in corrector-predictor for ea var */
  double *ecp_AC = NULL;
  int inode, idof, valid;
#ifdef PARALLEL
  double ecp_buf[MAX_VARIABLE_TYPES]; /* accumulated over all procs */
  int ncp_buf[MAX_VARIABLE_TYPES];    /* accumulated over all procs */
  double max_buf[MAX_VARIABLE_TYPES]; /* accumulated over all procs */
#endif
#ifdef DM_COORD_SCALE_PLEASE
  int bit_DM_scale = TRUE;
#else
  int bit_DM_scale = FALSE;
#endif

  static const char yo[] = "time_step_control";

  /*
   * Set counters to zero
   */
  num_unknowns = 0;
  Err_norm = 0.0;
  init_vec_value(ecp, 0.0, MAX_VARIABLE_TYPES);
  init_vec_value(max, 0.0, MAX_VARIABLE_TYPES);
  memset(ncp, 0, sizeof(int) * MAX_VARIABLE_TYPES);

  if (nAC > 0)
    ecp_AC = alloc_dbl_1(nAC, 0);

  for (i = 0; i < MAX_VARIABLE_TYPES; i++) {
    ecp[i] = 0.0;
    ncp[i] = 0;
    max[i] = 0.0;
  }

  /*
   * Collect deviations from predicted and max values
   * and ncp[eqn] => Number of unknowns in each binned variable type
   */
  for (i = 0; i < (num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx]); i++) {

    valid = TRUE;

    inode = 0;
    vd = Index_Solution_Inv(i, &inode, NULL, NULL, &idof, pg->imtrx);
    eqn = vd->Variable_Type;

#ifdef DEBUG_HKM
    if (eqn != idv[pg->imtrx][i][0]) {
      GOMA_EH(-1, "error in  Index_Solution_Inv mapping");
    }
#endif

    /* ignore corrections caused by changing side of interface */
    if (ls != NULL && xfem != NULL && ls->Length_Scale == 0.) {
      int MatID = vd->MatID;
      int interp;
      double F, Fold;
      if (MatID == -1)
        MatID = 0;
      interp = pd_glob[MatID]->i[pg->imtrx][eqn];

      if (is_xfem_interp(interp)) {
        gnn_distance(inode, x, x_old, NULL, &F, &Fold, NULL);
        if (is_extended_dof(inode, idof, vd, F)) {
          valid = FALSE;
        } else {
          if (sign_change(F, Fold)) {
            valid = FALSE;
          }
        }
      }
    }

    if (valid) {
      ecp[eqn] += SQUARE(x[i] - x_pred[i]);
      ncp[eqn]++;
      /* Set bit TRUE in next line to scale displacements with coordinates
       * instead of displacements to avoid large displacement error when
       * there is near zero displacement - i.e., better timestep control */
      if (bit_DM_scale &&
          (eqn == MESH_DISPLACEMENT1 || eqn == MESH_DISPLACEMENT2 || eqn == MESH_DISPLACEMENT3)) {
        if (fabs(Coor[eqn - MESH_DISPLACEMENT1][inode]) > max[eqn])
          max[eqn] = fabs(Coor[eqn - MESH_DISPLACEMENT1][inode]);
      } else {
        if (fabs(x[i]) > max[eqn])
          max[eqn] = fabs(x[i]);
      }
    }
  }

  if (nAC > 0) {
    for (i = 0; i < nAC; i++) {
      ecp_AC[i] = SQUARE(x_AC[i] - x_AC_pred[i]);
    }
  }
/*
 * Find global quantities
 */
#ifdef PARALLEL
  MPI_Allreduce((void *)ecp, (void *)ecp_buf, MAX_VARIABLE_TYPES, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce((void *)max, (void *)max_buf, MAX_VARIABLE_TYPES, MPI_DOUBLE, MPI_MAX,
                MPI_COMM_WORLD);
  MPI_Allreduce((void *)ncp, (void *)ncp_buf, MAX_VARIABLE_TYPES, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  for (i = 0; i < MAX_VARIABLE_TYPES; i++) {
    ecp[i] = ecp_buf[i];
    max[i] = max_buf[i];
    ncp[i] = ncp_buf[i];
  }
#endif

  /*
   * normalize the L_2 norm by the maximum value of an
   * variable type across the domain.
   */
  if (eps < 0.0) {
    for (i = 0; i < MAX_VARIABLE_TYPES; i++) {
      if (max[i] > DBL_SEMI_SMALL) {
#ifdef VAR_UPDATE_UNITY_SCALE
        ecp[i] = ecp[i] / (1.0 + SQUARE(max[i]));
#else
        ecp[i] = ecp[i] / SQUARE(max[i]);
#endif
      }
    }
    if (nAC > 0) {
      for (i = 0; i < nAC; i++) {
        if (fabs(x_AC[i]) > 0.0)
          ecp_AC[i] /= fabs(x_AC[i]);
      }
    }
  }

  /*
   * Construct a single error norm from each variable's contribution
   * depending on user selection specification in the input deck.
   */
  if (use_var_norm[0]) {
    Err_norm += ecp[MESH_DISPLACEMENT1];
    Err_norm += ecp[MESH_DISPLACEMENT2];
    Err_norm += ecp[MESH_DISPLACEMENT3];
    Err_norm += ecp[MAX_STRAIN];
    Err_norm += ecp[CUR_STRAIN];
    num_unknowns += ncp[MESH_DISPLACEMENT1];
    num_unknowns += ncp[MESH_DISPLACEMENT2];
    num_unknowns += ncp[MESH_DISPLACEMENT3];
    num_unknowns += ncp[MAX_STRAIN];
    num_unknowns += ncp[CUR_STRAIN];
  }

  if (use_var_norm[1]) {
    Err_norm += ecp[VELOCITY1];
    Err_norm += ecp[VELOCITY2];
    Err_norm += ecp[VELOCITY3];
    num_unknowns += ncp[VELOCITY1];
    num_unknowns += ncp[VELOCITY2];
    num_unknowns += ncp[VELOCITY3];

    /*
     * Looks like Matt's particle momentum is considered lumped together
     * with solvent phase momentum equations for the purposes of computing
     * a global time step norm...
     */
    Err_norm += ecp[PVELOCITY1];
    Err_norm += ecp[PVELOCITY2];
    Err_norm += ecp[PVELOCITY3];
    num_unknowns += ncp[PVELOCITY1];
    num_unknowns += ncp[PVELOCITY2];
    num_unknowns += ncp[PVELOCITY3];
  }

  if (use_var_norm[2]) {
    Err_norm += ecp[TEMPERATURE];
    num_unknowns += ncp[TEMPERATURE];
  }

  if (use_var_norm[3]) { /* Maybe someday we can discriminate between
                          * individual species at fault for predictor
                          * missing the corrector. */
    Err_norm += ecp[MASS_FRACTION];
    Err_norm += ecp[POR_LIQ_PRES];
    Err_norm += ecp[POR_GAS_PRES];
    Err_norm += ecp[POR_POROSITY];
    Err_norm += ecp[POR_SATURATION];
    Err_norm += ecp[POR_SINK_MASS];

    num_unknowns += ncp[MASS_FRACTION];
    num_unknowns += ncp[POR_LIQ_PRES];
    num_unknowns += ncp[POR_GAS_PRES];
    num_unknowns += ncp[POR_POROSITY];
    num_unknowns += ncp[POR_SATURATION];
    num_unknowns += ncp[POR_SINK_MASS];
  }

  if (use_var_norm[4]) { /* Pressure, even though there is no dP/dt
                          * term in any of the equations for
                          * incompressible flow... */
    Err_norm += ecp[PRESSURE];
    num_unknowns += ncp[PRESSURE];
  }

  if (use_var_norm[5]) { /* Polymer extra stress contribution */
    Err_norm += ecp[POLYMER_STRESS11];
    Err_norm += ecp[POLYMER_STRESS12];
    Err_norm += ecp[POLYMER_STRESS22];
    Err_norm += ecp[POLYMER_STRESS13];
    Err_norm += ecp[POLYMER_STRESS23];
    Err_norm += ecp[POLYMER_STRESS33];

    num_unknowns += ncp[POLYMER_STRESS11];
    num_unknowns += ncp[POLYMER_STRESS12];
    num_unknowns += ncp[POLYMER_STRESS22];
    num_unknowns += ncp[POLYMER_STRESS13];
    num_unknowns += ncp[POLYMER_STRESS23];
    num_unknowns += ncp[POLYMER_STRESS33];

    /*
     * Hey, lots of multi-mode Giesekus stress equations, too!
     */
    for (eqn = POLYMER_STRESS11_1; eqn <= POLYMER_STRESS33_7; eqn++) {
      Err_norm += ecp[eqn];
      num_unknowns += ncp[eqn];
    }
  }

  if (use_var_norm[6]) {
    Err_norm += ecp[VOLTAGE];
    num_unknowns += ncp[VOLTAGE];
  }

  if (use_var_norm[7]) {
    Err_norm += ecp[SOLID_DISPLACEMENT1];
    Err_norm += ecp[SOLID_DISPLACEMENT2];
    Err_norm += ecp[SOLID_DISPLACEMENT3];

    num_unknowns += ncp[SOLID_DISPLACEMENT1];
    num_unknowns += ncp[SOLID_DISPLACEMENT2];
    num_unknowns += ncp[SOLID_DISPLACEMENT3];
  }

  if (nAC > 0) /* Don't we want to include the AC unknowns in time step control ? */
  {
    if (use_var_norm[9]) { /* answer -> not usually since they are often algebraic constraints */
      for (i = 0; i < nAC; i++) {
        Err_norm += ecp_AC[i];
      }
      num_unknowns += nAC;
    }
  }
  /** add in shell element components  */

  Err_norm += ecp[SURF_CHARGE];
  Err_norm += ecp[SHELL_CURVATURE];
  Err_norm += ecp[SHELL_CURVATURE2];
  Err_norm += ecp[SHELL_TENSION];
  Err_norm += ecp[SHELL_X];
  Err_norm += ecp[SHELL_Y];
  Err_norm += ecp[SHELL_USER];
  Err_norm += ecp[ACOUS_PREAL];
  Err_norm += ecp[ACOUS_PIMAG];
  Err_norm += ecp[ACOUS_REYN_STRESS];
  Err_norm += ecp[SHELL_BDYVELO];
  Err_norm += ecp[SHELL_LUBP];
  Err_norm += ecp[SHELL_TEMPERATURE];
  Err_norm += ecp[SHELL_DELTAH];
  Err_norm += ecp[SHELL_FILMP];
  Err_norm += ecp[SHELL_FILMH];
  Err_norm += ecp[SHELL_PARTC];
  Err_norm += ecp[LIGHT_INTP];
  Err_norm += ecp[LIGHT_INTM];
  Err_norm += ecp[LIGHT_INTD];
  Err_norm += ecp[RESTIME];
  /*    Err_norm      += ecp[EXT_VELOCITY];  */
  Err_norm += ecp[TFMP_PRES];
  Err_norm += ecp[TFMP_SAT];
  Err_norm += ecp[SHELL_SURF_DIV_V];
  Err_norm += ecp[SHELL_SURF_CURV];
  Err_norm += ecp[SHELL_NORMAL1];
  Err_norm += ecp[SHELL_NORMAL2];
  Err_norm += ecp[SHELL_NORMAL3];
  Err_norm += ecp[N_DOT_CURL_V];
  Err_norm += ecp[GRAD_S_V_DOT_N1];
  Err_norm += ecp[GRAD_S_V_DOT_N2];
  Err_norm += ecp[GRAD_S_V_DOT_N3];

  num_unknowns += ncp[SURF_CHARGE];
  num_unknowns += ncp[SHELL_CURVATURE];
  num_unknowns += ncp[SHELL_CURVATURE2];
  num_unknowns += ncp[SHELL_TENSION];
  num_unknowns += ncp[SHELL_X];
  num_unknowns += ncp[SHELL_Y];
  num_unknowns += ncp[SHELL_USER];
  num_unknowns += ncp[ACOUS_PREAL];
  num_unknowns += ncp[ACOUS_PIMAG];
  num_unknowns += ncp[ACOUS_REYN_STRESS];
  num_unknowns += ncp[SHELL_BDYVELO];
  num_unknowns += ncp[SHELL_LUBP];
  num_unknowns += ncp[SHELL_TEMPERATURE];
  num_unknowns += ncp[SHELL_DELTAH];
  num_unknowns += ncp[SHELL_FILMP];
  num_unknowns += ncp[SHELL_FILMH];
  num_unknowns += ncp[SHELL_PARTC];
  num_unknowns += ncp[LIGHT_INTP];
  num_unknowns += ncp[LIGHT_INTM];
  num_unknowns += ncp[LIGHT_INTD];
  num_unknowns += ncp[RESTIME];
  /*    num_unknowns += ncp[EXT_VELOCITY];  */
  num_unknowns += ncp[TFMP_PRES];
  num_unknowns += ncp[TFMP_SAT];
  num_unknowns += ncp[SHELL_SURF_DIV_V];
  num_unknowns += ncp[SHELL_SURF_CURV];
  num_unknowns += ncp[SHELL_NORMAL1];
  num_unknowns += ncp[SHELL_NORMAL2];
  num_unknowns += ncp[SHELL_NORMAL3];
  num_unknowns += ncp[N_DOT_CURL_V];
  num_unknowns += ncp[GRAD_S_V_DOT_N1];
  num_unknowns += ncp[GRAD_S_V_DOT_N2];
  num_unknowns += ncp[GRAD_S_V_DOT_N3];

  if (use_var_norm[8]) /* LS equation is set with special card in Level Set section */
  {
    Err_norm += ecp[FILL];
    num_unknowns += ncp[FILL];
    Err_norm += ecp[PHASE1];
    num_unknowns += ncp[PHASE1];
  }

  if (use_var_norm[9]) {
    Err_norm += ecp[LUBP];
    Err_norm += ecp[LUBP_2];
    Err_norm += ecp[SHELL_SAT_CLOSED];
    Err_norm += ecp[SHELL_PRESS_OPEN];
    Err_norm += ecp[SHELL_PRESS_OPEN_2];
    Err_norm += ecp[SHELL_SAT_1];
    Err_norm += ecp[SHELL_SAT_2];
    Err_norm += ecp[SHELL_SAT_3];
    Err_norm += ecp[SHELL_SAT_GASN];
    Err_norm += ecp[SHELL_LUB_CURV];
    Err_norm += ecp[SHELL_LUB_CURV_2];
    num_unknowns += ncp[LUBP];
    num_unknowns += ncp[LUBP_2];
    num_unknowns += ncp[SHELL_SAT_CLOSED];
    num_unknowns += ncp[SHELL_PRESS_OPEN];
    num_unknowns += ncp[SHELL_PRESS_OPEN_2];
    num_unknowns += ncp[SHELL_SAT_1];
    num_unknowns += ncp[SHELL_SAT_2];
    num_unknowns += ncp[SHELL_SAT_3];
    num_unknowns += ncp[SHELL_SAT_GASN];
    num_unknowns += ncp[SHELL_LUB_CURV];
    num_unknowns += ncp[SHELL_LUB_CURV_2];
  }

#if 0 /* ------------------- maybe someday you'll want these, too... -----*/
  if (use_var_norm["index for shear rate equation"] ) {
    Err_norm      += ecp[SHEAR_RATE];
    num_unknowns += ncp[SHEAR_RATE];
  }
  if (use_var_norm["index for vorticity principle shear directions"] )  {
    Err_norm      += ecp[VORT_DIR1];
    Err_norm      += ecp[VORT_DIR2];
    Err_norm      += ecp[VORT_DIR3];
    Err_norm      += ecp[VORT_LAMBDA];
    num_unknowns += ncp[VORT_DIR1];
    num_unknowns += ncp[VORT_DIR2];
    num_unknowns += ncp[VORT_DIR3];
    num_unknowns += ncp[VORT_LAMBDA];
  }

  if (use_var_norm["index for suspension temperature equation"] ) {
    Err_norm      += ecp[BOND_EVOLUTION];
    num_unknowns += ncp[BOND_EVOLUTION];
  }
#endif

  if (num_unknowns == 0) {
    DPRINTF(stderr, "\"Time step error\" norm includes no active variables!\n");
    DPRINTF(stderr, "You specified (d=%d, v=%d, T=%d, y=%d, P=%d, S=%d, V=%d)\n", use_var_norm[0],
            use_var_norm[1], use_var_norm[2], use_var_norm[3], use_var_norm[4], use_var_norm[5],
            use_var_norm[6]);
    GOMA_EH(-1, "Poorly formed time step norm.");
  }

  scaling = 1.0 / (num_unknowns * (2.0 + delta_t_old / delta_t));
  Err_norm *= scaling;
  Err_norm = sqrt(Err_norm);

  /*
   * Bin the individual variable type errors. Then, scale them.
   */
  e_d = (ecp[MESH_DISPLACEMENT1] + ecp[MESH_DISPLACEMENT2] + ecp[MESH_DISPLACEMENT3]) +
        ecp[MAX_STRAIN] + ecp[CUR_STRAIN];
  e_v = ecp[VELOCITY1] + ecp[VELOCITY2] + ecp[VELOCITY3];
  e_T = ecp[TEMPERATURE];
  e_y = (ecp[MASS_FRACTION] + ecp[POR_LIQ_PRES] + ecp[POR_GAS_PRES] + ecp[POR_POROSITY] +
         ecp[POR_SATURATION] + ecp[POR_SINK_MASS]);
  e_P = ecp[PRESSURE];
  e_S = (ecp[POLYMER_STRESS11] + ecp[POLYMER_STRESS12] + ecp[POLYMER_STRESS22] +
         ecp[POLYMER_STRESS13] + ecp[POLYMER_STRESS23] + ecp[POLYMER_STRESS33]);
  for (eqn = POLYMER_STRESS11_1; eqn <= POLYMER_STRESS33_7; eqn++) {
    e_S += ecp[eqn];
  }
  for (i = 0, e_AC = 0.0; i < nAC; i++)
    e_AC += ecp_AC[i];

  e_V = ecp[VOLTAGE];
  e_qs = ecp[SURF_CHARGE];
  e_shk = ecp[SHELL_CURVATURE] + ecp[SHELL_CURVATURE2];
  e_sht = ecp[SHELL_TENSION];
  e_shd = ecp[SHELL_X] + ecp[SHELL_Y];
  e_shu = ecp[SHELL_USER];
  e_sh_lub = ecp[LUBP] + ecp[LUBP_2] + ecp[SHELL_FILMP] + ecp[SHELL_FILMH] + ecp[SHELL_PARTC] +
             ecp[SHELL_SAT_CLOSED] + ecp[SHELL_SAT_GASN] + ecp[SHELL_TEMPERATURE] +
             ecp[SHELL_DELTAH] + ecp[SHELL_LUB_CURV] + ecp[SHELL_LUB_CURV_2] +
             ecp[SHELL_PRESS_OPEN] + ecp[SHELL_PRESS_OPEN_2] + ecp[SHELL_SAT_1] + ecp[SHELL_SAT_2] +
             ecp[SHELL_SAT_3];
  e_F = ecp[FILL] + ecp[PHASE1];
  e_ap = ecp[ACOUS_PREAL] + ecp[ACOUS_PIMAG] + ecp[ACOUS_REYN_STRESS];
  e_extv = ecp[EXT_VELOCITY];
  e_int = ecp[LIGHT_INTP] + ecp[LIGHT_INTM] + ecp[LIGHT_INTD] + ecp[RESTIME];
  e_rheo = ecp[SHELL_SURF_DIV_V] + ecp[SHELL_SURF_CURV] + ecp[SHELL_NORMAL1] + ecp[SHELL_NORMAL2] +
           ecp[SHELL_NORMAL3] + ecp[N_DOT_CURL_V] + ecp[GRAD_S_V_DOT_N1] + ecp[GRAD_S_V_DOT_N2] +
           ecp[GRAD_S_V_DOT_N3];

  e_d = sqrt(e_d * scaling);
  e_v = sqrt(e_v * scaling);
  e_T = sqrt(e_T * scaling);
  e_y = sqrt(e_y * scaling);
  e_P = sqrt(e_P * scaling);
  e_S = sqrt(e_S * scaling);
  e_V = sqrt(e_V * scaling);
  e_qs = sqrt(e_qs * scaling);
  e_AC = sqrt(e_AC * scaling);
  e_shk = sqrt(e_shk * scaling);
  e_sht = sqrt(e_sht * scaling);
  e_shd = sqrt(e_shd * scaling);
  e_shu = sqrt(e_shu * scaling);
  e_F = sqrt(e_F * scaling);
  e_ap = sqrt(e_ap * scaling);
  e_extv = sqrt(e_extv * scaling);
  e_sh_lub = sqrt(e_sh_lub * scaling);
  e_int = sqrt(e_int * scaling);
  e_rheo = sqrt(e_rheo * scaling);

  /*
   * Print out the breakdown of contributions as well as the user specified
   * selection of which contributions actually get used in computing the
   * overall norm.
   */

  /*
   * Determine whether the step is successful or not from the
   * formula below.
   */
  *success_dt = (Err_norm < beta * abs_eps);
  if (const_delta_t) {
    DPRINTF(stdout, "\nCONSTANT DELTA_T          [");
    if (ncp[MESH_DISPLACEMENT1]) {
      DPRINTF(stdout, "%7.1e", e_d);
    }
    if (ncp[VELOCITY1] || ncp[PVELOCITY1]) {
      DPRINTF(stdout, ", %7.1e", e_v);
    }
    if (ncp[TEMPERATURE]) {
      DPRINTF(stdout, ", %7.1e", e_T);
    }
    if (ncp[MASS_FRACTION] || ncp[POR_LIQ_PRES] || ncp[POR_GAS_PRES] || ncp[POR_POROSITY] ||
        ncp[POR_SATURATION] || ncp[POR_SINK_MASS]) {
      DPRINTF(stdout, ", %7.1e", e_y);
    }
    if (ncp[PRESSURE]) {
      DPRINTF(stdout, ", %7.1e", e_P);
    }
    if (ncp[POLYMER_STRESS11]) {
      DPRINTF(stdout, ", %7.1e", e_S);
    }
    if (ncp[VOLTAGE]) {
      DPRINTF(stdout, ", %7.1e", e_V);
    }
    if (ncp[SURF_CHARGE]) {
      DPRINTF(stdout, ", %7.1e", e_qs);
    }
    if (ncp[SHELL_CURVATURE]) {
      DPRINTF(stdout, ", %7.1e", e_shk);
    }
    if (ncp[SHELL_TENSION]) {
      DPRINTF(stdout, ", %7.1e", e_sht);
    }
    if (ncp[SHELL_X]) {
      DPRINTF(stdout, ", %7.1e", e_shd);
    }
    if (ncp[SHELL_USER]) {
      DPRINTF(stdout, ", %7.1e", e_shu);
    }
    if (ncp[FILL] || ncp[PHASE1]) {
      DPRINTF(stdout, ", %7.1e", e_F);
    }
    if (ncp[ACOUS_PREAL] || ncp[ACOUS_PIMAG]) {
      DPRINTF(stdout, ", %7.1e", e_ap);
    }
    if (ncp[EXT_VELOCITY]) {
      DPRINTF(stdout, ", %7.1e", e_extv);
    }
    if (ncp[LIGHT_INTP] || ncp[LIGHT_INTM] || ncp[LIGHT_INTD] || ncp[RESTIME]) {
      DPRINTF(stdout, ", %7.1e", e_int);
    }
    if (ncp[LUBP] || ncp[LUBP_2] || ncp[SHELL_FILMP] || ncp[SHELL_TEMPERATURE] ||
        ncp[SHELL_DELTAH] || ncp[SHELL_LUB_CURV] || ncp[SHELL_LUB_CURV_2] ||
        ncp[SHELL_SAT_CLOSED] || ncp[SHELL_PRESS_OPEN] || ncp[SHELL_PRESS_OPEN_2] ||
        ncp[SHELL_SAT_1] || ncp[SHELL_SAT_2] || ncp[SHELL_SAT_3]) {
      DPRINTF(stdout, ", %7.1e", e_sh_lub);
    }
    if (ncp[SHELL_NORMAL1] || ncp[SHELL_SURF_CURV] || ncp[GRAD_S_V_DOT_N1]) {
      DPRINTF(stdout, ", %7.1e", e_rheo);
    }
    if (nAC > 0) {
      DPRINTF(stdout, ", %7.1e", e_AC);
    }
    DPRINTF(stdout, "]\n");
    log_msg("Constant delta_t");
  } else if (*success_dt) {
    DPRINTF(stdout, "\nOK  %7.1e < %3g %7.1e [", Err_norm, beta, abs_eps);
    if (ncp[MESH_DISPLACEMENT1]) {
      DPRINTF(stdout, "%7.1e", e_d);
    }
    if (ncp[VELOCITY1] || ncp[PVELOCITY1]) {
      DPRINTF(stdout, ", %7.1e", e_v);
    }
    if (ncp[TEMPERATURE]) {
      DPRINTF(stdout, ", %7.1e", e_T);
    }
    if (ncp[MASS_FRACTION] || ncp[POR_LIQ_PRES] || ncp[POR_GAS_PRES] || ncp[POR_POROSITY] ||
        ncp[POR_SATURATION] || ncp[POR_SINK_MASS]) {
      DPRINTF(stdout, ", %7.1e", e_y);
    }
    if (ncp[PRESSURE]) {
      DPRINTF(stdout, ", %7.1e", e_P);
    }
    if (ncp[POLYMER_STRESS11]) {
      DPRINTF(stdout, ", %7.1e", e_S);
    }
    if (ncp[VOLTAGE]) {
      DPRINTF(stdout, ", %7.1e", e_V);
    }
    if (ncp[SURF_CHARGE]) {
      DPRINTF(stdout, ", %7.1e", e_qs);
    }
    if (ncp[SHELL_CURVATURE]) {
      DPRINTF(stdout, ", %7.1e", e_shk);
    }
    if (ncp[SHELL_TENSION]) {
      DPRINTF(stdout, ", %7.1e", e_sht);
    }
    if (ncp[SHELL_X]) {
      DPRINTF(stdout, ", %7.1e", e_shd);
    }
    if (ncp[SHELL_USER]) {
      DPRINTF(stdout, ", %7.1e", e_shu);
    }
    if (ncp[FILL] || ncp[PHASE1]) {
      DPRINTF(stdout, ", %7.1e", e_F);
    }
    if (ncp[ACOUS_PREAL] || ncp[ACOUS_PIMAG]) {
      DPRINTF(stdout, ", %7.1e", e_ap);
    }
    if (ncp[EXT_VELOCITY]) {
      DPRINTF(stdout, ", %7.1e", e_extv);
    }
    if (ncp[LIGHT_INTP] || ncp[LIGHT_INTM] || ncp[LIGHT_INTD] || ncp[RESTIME]) {
      DPRINTF(stdout, ", %7.1e", e_int);
    }
    if (ncp[LUBP] || ncp[LUBP_2] || ncp[SHELL_FILMP] || ncp[SHELL_TEMPERATURE] ||
        ncp[SHELL_DELTAH] || ncp[SHELL_LUB_CURV] || ncp[SHELL_LUB_CURV_2] ||
        ncp[SHELL_SAT_CLOSED] || ncp[SHELL_PRESS_OPEN] || ncp[SHELL_PRESS_OPEN_2] ||
        ncp[SHELL_SAT_1] || ncp[SHELL_SAT_2] || ncp[SHELL_SAT_3]) {
      DPRINTF(stdout, ", %7.1e", e_sh_lub);
    }
    if (ncp[SHELL_NORMAL1] || ncp[SHELL_SURF_CURV] || ncp[GRAD_S_V_DOT_N1]) {
      DPRINTF(stdout, ", %7.1e", e_rheo);
    }
    if (nAC > 0) {
      DPRINTF(stdout, ", %7.1e", e_AC);
    }
    DPRINTF(stdout, "]\n");
    log_msg("Predictor was OK, %g < %g * %g", Err_norm, beta, eps);
  } else {
    DPRINTF(stdout, "\nYUK %7.1e > %3g %7.1e [", Err_norm, beta, abs_eps);
    if (ncp[MESH_DISPLACEMENT1]) {
      DPRINTF(stdout, "%7.1e", e_d);
    }
    if (ncp[VELOCITY1] || ncp[PVELOCITY1]) {
      DPRINTF(stdout, ", %7.1e", e_v);
    }
    if (ncp[TEMPERATURE]) {
      DPRINTF(stdout, ", %7.1e", e_T);
    }
    if (ncp[MASS_FRACTION] || ncp[POR_LIQ_PRES] || ncp[POR_GAS_PRES] || ncp[POR_POROSITY] ||
        ncp[POR_SATURATION] || ncp[POR_SINK_MASS]) {
      DPRINTF(stdout, ", %7.1e", e_y);
    }
    if (ncp[PRESSURE]) {
      DPRINTF(stdout, ", %7.1e", e_P);
    }
    if (ncp[POLYMER_STRESS11]) {
      DPRINTF(stdout, ", %7.1e", e_S);
    }
    if (ncp[VOLTAGE]) {
      DPRINTF(stdout, ", %7.1e", e_V);
    }
    if (ncp[SURF_CHARGE]) {
      DPRINTF(stdout, ", %7.1e", e_qs);
    }
    if (ncp[SHELL_CURVATURE]) {
      DPRINTF(stdout, ", %7.1e", e_shk);
    }
    if (ncp[SHELL_TENSION]) {
      DPRINTF(stdout, ", %7.1e", e_sht);
    }
    if (ncp[SHELL_X]) {
      DPRINTF(stdout, ", %7.1e", e_shd);
    }
    if (ncp[SHELL_USER]) {
      DPRINTF(stdout, ", %7.1e", e_shu);
    }
    if (ncp[FILL] || ncp[PHASE1]) {
      DPRINTF(stdout, ", %7.1e", e_F);
    }
    if (ncp[ACOUS_PREAL] || ncp[ACOUS_PIMAG]) {
      DPRINTF(stdout, ", %7.1e", e_ap);
    }
    if (ncp[EXT_VELOCITY]) {
      DPRINTF(stdout, ", %7.1e", e_extv);
    }
    if (ncp[LIGHT_INTP] || ncp[LIGHT_INTM] || ncp[LIGHT_INTD] || ncp[RESTIME]) {
      DPRINTF(stdout, ", %7.1e", e_int);
    }
    if (ncp[LUBP] || ncp[LUBP_2] || ncp[SHELL_FILMP] || ncp[SHELL_TEMPERATURE] ||
        ncp[SHELL_DELTAH] || ncp[SHELL_LUB_CURV] || ncp[SHELL_LUB_CURV_2] ||
        ncp[SHELL_SAT_CLOSED] || ncp[SHELL_PRESS_OPEN] || ncp[SHELL_PRESS_OPEN_2] ||
        ncp[SHELL_SAT_1] || ncp[SHELL_SAT_2] || ncp[SHELL_SAT_3]) {
      DPRINTF(stdout, ", %7.1e", e_sh_lub);
    }
    if (ncp[SHELL_NORMAL1] || ncp[SHELL_SURF_CURV] || ncp[GRAD_S_V_DOT_N1]) {
      DPRINTF(stdout, ", %7.1e", e_rheo);
    }
    if (nAC > 0) {
      DPRINTF(stdout, ", %7.1e", e_AC);
    }
    DPRINTF(stdout, "]\n");
    log_msg("Predictor was YUK, %g > %g * %g", Err_norm, beta, eps);
  }
  DPRINTF(stdout, "Predictor-corrector norm: [");
  if (ncp[MESH_DISPLACEMENT1]) {
    DPRINTF(stdout, "  %1d d  ", use_var_norm[0]);
  }
  if (ncp[VELOCITY1] || ncp[PVELOCITY1]) {
    DPRINTF(stdout, ",   %1d v  ", use_var_norm[1]);
  }
  if (ncp[TEMPERATURE]) {
    DPRINTF(stdout, ",   %1d T  ", use_var_norm[2]);
  }
  if (ncp[MASS_FRACTION] || ncp[POR_LIQ_PRES] || ncp[POR_GAS_PRES] || ncp[POR_POROSITY] ||
      ncp[POR_SATURATION] || ncp[POR_SINK_MASS]) {
    DPRINTF(stdout, ",   %1d y  ", use_var_norm[3]);
  }
  if (ncp[PRESSURE]) {
    DPRINTF(stdout, ",   %1d P  ", use_var_norm[4]);
  }
  if (ncp[POLYMER_STRESS11]) {
    DPRINTF(stdout, ",   %1d S  ", use_var_norm[5]);
  }
  if (ncp[VOLTAGE]) {
    DPRINTF(stdout, ",   %1d V  ", use_var_norm[6]);
  }
  if (ncp[SURF_CHARGE]) {
    DPRINTF(stdout, ",   %1d Q  ", 1);
  }
  if (ncp[SHELL_CURVATURE]) {
    DPRINTF(stdout, ",   %1d K  ", 1);
  }
  if (ncp[SHELL_TENSION]) {
    DPRINTF(stdout, ",   %1d ST ", 1);
  }
  if (ncp[SHELL_X]) {
    DPRINTF(stdout, ",   %1d Sd ", 1);
  }
  if (ncp[SHELL_USER]) {
    DPRINTF(stdout, ",   %1d Su ", 1);
  }
  if (ncp[FILL] || ncp[PHASE1]) {
    DPRINTF(stdout, ",   %1d F  ", 1);
  }
  if (ncp[ACOUS_PREAL] || ncp[ACOUS_PIMAG]) {
    DPRINTF(stdout, ",   %1d A  ", 1);
  }
  if (ncp[EXT_VELOCITY]) {
    DPRINTF(stdout, ",   %1d Ev  ", 1);
  }
  if (ncp[LIGHT_INTP] || ncp[LIGHT_INTM] || ncp[LIGHT_INTD] || ncp[RESTIME]) {
    DPRINTF(stdout, ", %1d INT", 1);
  }
  if (ncp[LUBP] || ncp[LUBP_2] || ncp[SHELL_FILMP] || ncp[SHELL_TEMPERATURE] || ncp[SHELL_DELTAH] ||
      ncp[SHELL_LUB_CURV] || ncp[SHELL_LUB_CURV_2] || ncp[SHELL_SAT_CLOSED] ||
      ncp[SHELL_PRESS_OPEN] || ncp[SHELL_PRESS_OPEN_2] || ncp[SHELL_SAT_1] || ncp[SHELL_SAT_2] ||
      ncp[SHELL_SAT_3]) {
    DPRINTF(stdout, ", %1d SHELL", 1);
  }
  if (ncp[SHELL_NORMAL1] || ncp[SHELL_SURF_CURV] || ncp[GRAD_S_V_DOT_N1]) {
    DPRINTF(stdout, ", %1d RHEO", 1);
  }
  if (nAC > 0) {
    DPRINTF(stdout, ",   %1d AC ", use_var_norm[9]);
  }
  DPRINTF(stdout, "]\n");

  log_msg("Predictor details:");
  log_msg("                 %1d  e_d = %g", use_var_norm[0], e_d);
  log_msg("                 %1d  e_v = %g", use_var_norm[1], e_v);
  log_msg("                 %1d  e_T = %g", use_var_norm[2], e_T);
  log_msg("                 %1d  e_y = %g", use_var_norm[3], e_y);
  log_msg("                 %1d  e_P = %g", use_var_norm[4], e_P);
  log_msg("                 %1d  e_S = %g", use_var_norm[5], e_S);
  log_msg("                 %1d  e_V = %g", use_var_norm[6], e_V);
  log_msg("                 %1d  e_sh_lub = %g", use_var_norm[9], e_sh_lub);
  if (nAC > 0)
    log_msg("                 %1d  e_AC = %g", use_var_norm[9], e_AC);

  /*
   * Calculate the time step to be used in the next time
   */
  if (const_delta_t) {
    delta_t_new = delta_t;
  } else if (!*success_dt) {
    /*
     * If the current time iteration had too much time step truncation
     * error, decrease the time step by a factor of two and return
     */
    delta_t_new = delta_t / 2.;
    return (delta_t_new);
  } else {
    if (Err_norm <= 0.0) {
      if (abs_eps > 0.0) {
        delta_t_new = TIME_STEP_GROWTH_CAP * delta_t;
      } else {
        delta_t_new = delta_t;
      }
    } else if (Err_norm > abs_eps) {
      delta_t_new = delta_t * pow(abs_eps / Err_norm, 1. / 3.);
    } else if (Err_norm <= abs_eps && Err_norm >= abs_eps / alpha) {
      delta_t_new = delta_t;
    } else if (Err_norm < abs_eps / alpha) {
      delta_t_new = delta_t * pow(abs_eps / (alpha * Err_norm), 1. / 3.);
    }
  }

  /*
   * ensure time step doesn't suddenly become too large
   */
  if (delta_t_new / delta_t > TIME_STEP_GROWTH_CAP) {
    delta_t_new = TIME_STEP_GROWTH_CAP * delta_t;
  }

  if (nAC > 0)
    safer_free((void **)&ecp_AC);
  return (delta_t_new);
} /* END of routine time_step_control */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double path_step_control(int N,
                         double delta_s_old,
                         double delta_s_older,
                         double x[],
                         double eps,
                         int *success_ds,
                         int use_var_norm[],
                         int inewton)

/*
   N - number of nodes
   delta_s_old - time step size last successful time step (n)
   delta_s_older - time step size two time steps previous (n -1)
   x - solution vector at n calculated from implicit corrector
   alpha - path step control parameter
   beta - path step control parameter
   eps - path step control parameter
   success_ds - 1 if the path step was successful and 0 if it failed
   use_var_norm - this tells us whether or not to use a variable for the norm
                calculations. 0 means don't use ... 1 means use. The choice is
                made in the input deck
   inewton - number of newton iterations

   BASED ON time_step_control()
   Robert Secor
*/
{
  double delta_s, beta;
  int iter_desired;

#ifdef DEBUG
  static const char yo[] = "path_step_control";
#endif

  /* EXTERNAL VARIABLES */

  *success_ds = 0;

  /**  stepsize based on number of newton iterations
   *    for now have desired number of iterations equal
   *    to 1/2 of the max Newton steps -- someday probably
   *    hook up to the input file
   */
  iter_desired = Max_Newton_Steps / 2;
  beta = (pow(2.0, (float)(iter_desired - inewton)) - 1.0) / (double)iter_desired;
  if (beta > 1.0)
    beta = 3.0;
  else
    beta = pow(10.0, beta);
  if (beta > 3.0)
    beta = 3.0;
  if (beta < 0.5)
    beta = 0.5;

  delta_s = beta * delta_s_old;

  *success_ds = 1;

  return (delta_s);

} /* END of routine path_step_control */
/*****************************************************************************/

#if 0
int find_max ( int list_length,
	       int list[] )

     /*
       Function which finds the largest integer from a vector of integers
       
       Author:          Scott Hutchinson (1421)
       Date:            8 December 1992
       Revised:         8 December 1992
       
     */
     
     /******************************************************************************/
{
  /* LOCAL VARIABLES */
  register int i, max;
  
  /*************************** execution begins *******************************/
  
  if (list_length > 0) {
    max = list[0];
    for (i = 1; i < list_length; i++)
      if (list[i] > max) max = list[i];
    return (max);
  } else 
    return (INT_MIN);

} /* END of routine find_max */
/******************************************************************************/

int
find_min ( int list_length,
           int list[] )
     
/*
 *     Function which finds the smallest integer from a vector of integers
 *      
 *      Author:          Scott Hutchinson (1421)
 *      Date:            8 December 1992
 *      Revised:         8 December 1992
 *      
 */
     
{
/* LOCAL VARIABLES */
  register int i, min;
  
  /*************************** execution begins ******************************/
  if (list_length > 0) {
    min = list[0];
    for (i = 1; i < list_length; i++)
      if (list[i] < min)  min = list[i];
    return (min);
  } else
      return (INT_MAX);

} /* END of routine find_min */

/*****************************************************************************/
/*****************************************************************************/

int
find_inter ( int *ipos,
             int inter_ptr[],
             int set1[],
             int set2[],
             int length1,
             int length2  )
/*
 *
 *      Function which finds the intersection of two integer lists
 *      
 *      Author:          Scott Hutchinson (1421)
 *      Date:            8 December 1992
 *      Revised:         8 December 1992
 *      
 */
   
{
/* LOCAL VARIABLES */
  register int i, j, itemp, intersect = 0;
  
  /*************************** execution begins ******************************/

  itemp = *ipos;

  for (i = 0; i < length1; i++)
    for (j = 0; j < length2; j++)
      if (set1[i] == set2[j])
	inter_ptr[itemp++] = i;

  if (itemp == *ipos)
    intersect = -1;

  *ipos = itemp;

  return(intersect);

} /*  END of routine find_inter */
#endif
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void init_vec(
    double u[], Comm_Ex *cx, Exo_DB *exo, Dpi *dpi, double uAC[], int nAC, double *timeValueRead)

/************************************************************************
 *
 * init_vec():
 *
 *	Initialize the solution vector to a specific initial
 *	values determined by the flag "guess_flag".
 *      Also, initialize constant external variable vectors,
 *       say for a CD-type analysis
 *      Also, initialize and convert pixel images and direct mapping onto
 *       the base mesh.   (PRS - 7/7/2011)
 *
 *  NOTE:
 *   Initialize dofs owned by this processor first. Then, do an
 *   exchange_dof() operation to communicate values to ghost nodes
 ************************************************************************/
{
  int ebi; /* element block index */
  int a, i, iAC, j, k, n, w, err, ndof;
  int e_start, e_end, ielem, ielem_type, num_nodes, index, nunks;
  int mn, ipos, interpType, retn, var;
  static const char yo[] = "init_vec";
  double *dum_var;
  double *frac_vec;
  NODAL_VARS_STRUCT *nv;
  int *block_order;
  int block_temp;
  int ebj;
  struct Data_Table *ext_table = NULL;
  double timeValueReadExt = 0.0;

  /***********************************************************************/
  /*                   PHASE 1 -> global initializations                 */
  /***********************************************************************/

  // First thing we do is aug vector with existing state of the system
  // This state may then be overwritten below
  if (nAC > 0) {
    for (i = 0; i < nAC; i++) {
      load_extra_unknownsAC(i, uAC, cx, exo, dpi);
    }
  }

  switch (Guess_Flag) {
  case 0:
    /* Initialize all of the solution vector to zero */
    for (i = 0; i < NumUnknowns[pg->imtrx]; i++) {
      u[i] = 0.0;
    }

    /* Some Special Cases */

    /* Structural shell variables */
    if (upd->vp[pg->imtrx][SHELL_X] && upd->vp[pg->imtrx][SHELL_Y]) {
      init_structural_shell_coord(u);
    }

    /* Shell normal vector unknowns requiring initialization */
    if (upd->vp[pg->imtrx][SHELL_NORMAL1] > -1 && upd->vp[pg->imtrx][SHELL_NORMAL2] > -1) {
      init_shell_normal_unknowns(u, exo);
    }
    break;

  case 1:
    /*
     *  Initialize the solution vector to random
     *  numbers between 0 and 1
     */
    fill_dvec_rand(u, NumUnknowns[pg->imtrx]);
    break;

  case 2:
    /* Initialize all of the solution vector to  the value of 1 */
    for (i = 0; i < NumUnknowns[pg->imtrx]; i++) {
      u[i] = 1.0;
    }
    break;

  case 4:
    /* Initialize a solution vector from an ascii file */
    read_initial_guess(u, NumUnknowns[pg->imtrx], uAC, nAC);
    break;

  case 5:
    /*
     * Initialize the solution vector from a stored solution
     * in the exodus input file
     */
    DPRINTF(stdout, "\nInitial guess read from \"%s\" ...(last soln in file)\n", ExoFile);
    err = rd_vectors_from_exoII(u, ExoFile, 0, 0, INT_MAX, timeValueRead, exo);
    if (err != 0) {
      DPRINTF(stderr, "%s: err fr rd_vectors_from_exoII()\n", yo);
      exit(-1);
    }
    DPRINTF(stdout, "\t\t Values read time plane at time = %g\n", *timeValueRead);

    /*
     * Initialize augmenting conditions for stored solution in exodus input fil
     * As of Mar 18, 2002, brkfix doesn't support global variables, hence, this
     * capability does not exist for parallel operations.
     */
    if (nAC > 0) {
      int ngv;

      /* 	  DPRINTF(stdout, */
      /* 		  "%s:  reading augmenting conditions initial guess from \"%s\" ...\n",yo,
       * ExoFile); */

      ngv = rd_globals_from_exoII(uAC, ExoFile, 6, nAC);
      if (ngv < 0)
        DPRINTF(stderr, "%s: error trying to read augmenting conditions\n", yo);

      /* Update parameters prior to solution with these values */

      for (iAC = 0; iAC + 6 < ngv; iAC++) {
        DPRINTF(stdout, "AUGC[%d] initial guess :%6.3g found in exoII database - reading\n", iAC,
                uAC[iAC]);
        update_parameterAC(iAC, NULL, NULL, uAC, cx, exo, dpi);
      }
    }

    break;

  case 6:
    /*
     * Initialize the solution vector from a stored solution
     * in a different, specified exodus file
     */
    if (ExoTimePlane < INT_MAX) {
      DPRINTF(stdout, "\nInitial guess read from \"%s\" ...(soln plane = %d)\n", ExoAuxFile,
              ExoTimePlane);
    } else {
      DPRINTF(stdout, "\nInitial guess read from \"%s\" ...(last soln in file)\n", ExoAuxFile);
    }
    err = rd_vectors_from_exoII(u, ExoAuxFile, 0, 0, ExoTimePlane, timeValueRead, exo);
    if (err != 0) {
      DPRINTF(stderr, "%s:  err fr rd_vectors_from_exoII()\n", yo);
    }
    DPRINTF(stdout, "\t\t Values read time plane at time = %g\n", *timeValueRead);
    /*
     * check for variables which have bounds applied to overide the values from exoII files
     */
    if (Num_Var_Bound > 0) {
      int idv;
      var = Var_init[0].var;
      if (pd->v[pg->imtrx][var]) {
        for (i = 0; i < DPI_ptr->num_owned_nodes; i++) {
          idv = Index_Solution(i, var, 0, 0, -1, pg->imtrx);
          if (idv != -1) {
            if (u[idv] < Var_init[0].init_val_min)
              u[idv] = Var_init[0].init_val_min;
            if (u[idv] > Var_init[0].init_val_max)
              u[idv] = Var_init[0].init_val_max;
          }
        }
      }
      exchange_dof(cx, dpi, u, pg->imtrx);
    }

    /*
     * Initialize augmenting conditions for storred solution in exodus input file
     * As of Mar 18, 2002, brkfix doesn't support global variables, hence, this
     * capability does not exist for parallel operations.
     */

    if (nAC > 0) {
      int ngv;

      DPRINTF(stdout, "%s:  reading augmenting conditions initial guess from \"%s\" ...\n", yo,
              ExoAuxFile);

      ngv = rd_globals_from_exoII(uAC, ExoAuxFile, 6, nAC);
      if (ngv < 0)
        DPRINTF(stderr, "%s: error trying to read augmenting conditions\n", yo);

      /* Update parameters prior to solution with these values
       * Update only those parameters values that were actually found */

      for (iAC = 0; iAC + 6 < ngv; iAC++) {
        DPRINTF(stdout, "AUGC[%d] initial guess :%6.3g found in exoII database - reading\n", iAC,
                uAC[iAC]);
        update_parameterAC(iAC, NULL, NULL, uAC, cx, exo, dpi);
      }
    }

    break;

  default:
    DPRINTF(stderr, "%s:  unknown Guess_Flag\n", yo);
    exit(-1);
  }
  /***********************************************************************/
  /*       PHASE 1a -> initializations of AC variables from input file   */
  /***********************************************************************/
  if (nAC > 0) {
    for (iAC = 0; iAC < nAC; iAC++) {
      if (augc[iAC].iread == 2) {
        uAC[iAC] = augc[iAC].tmp3;
        DPRINTF(stdout, "AUGC[%d] initialized from input file :%6.3g \n", iAC, uAC[iAC]);
        update_parameterAC(iAC, NULL, NULL, uAC, cx, exo, dpi);
      }
    }
  }

  /***********************************************************************/
  /*       PHASE 2 -> global initializations of specific variable types  */
  /***********************************************************************/
  /*
   * check for variables which have specified initialization to
   * overide the values from input files, etc.
   */
  if (Num_Var_Init > 0) {
    retn = 0;
    dum_var = alloc_dbl_1(DPI_ptr->num_universe_nodes, DBL_NOINIT);
    for (i = Num_Var_Bound; i < Num_Var_Init; i++) {
      switch (Var_init[i].var) {
      case MASS_FRACTION:
        DPRINTF(stdout, "\tSetting %s number %d (variable [%d]) to %g\n",
                Var_Name[Var_init[i].var].name1, Var_init[i].ktype, Var_init[i].var,
                Var_init[i].init_val);

        break;
      case SPECIES_MASS_FRACTION:
      case SPECIES_MOLE_FRACTION:
      case SPECIES_VOL_FRACTION:
        retn = 1;
        continue;
      default:
        DPRINTF(stdout, "\tSetting %s (variable [%d]) to %g\n", Var_Name[Var_init[i].var].name1,
                Var_init[i].var, Var_init[i].init_val);
        break;
      }
      if (1 || Var_init[i].len_u_pars == -1) { /* Disabled for now, still some bugs */
        init_vec_value(dum_var, Var_init[i].init_val, DPI_ptr->num_universe_nodes);
        inject_nodal_vec(u, Var_init[i].var, Var_init[i].ktype, 0, -2, dum_var);
      } else {
        double xpt[DIM] = {0, 0, 0}, var_val[MAX_VARIABLE_TYPES];
        int dir, ierr;
        for (dir = 0; dir < pd->Num_Dim; dir++) {
          xpt[dir] = Coor[dir][i];
        }
        init_vec_value(var_val, 0.0, MAX_VARIABLE_TYPES);
        ierr = user_initialize(Var_init[i].var, u, Var_init[i].init_val, Var_init[i].u_pars, xpt,
                               var_val);
        if (ierr == -1) {
          GOMA_EH(ierr, "Problem with user_initialize...");
        }
      }
    }
    /*
     *    Section to calculate a consistent set of species fractions
     *    before injection into the solution vector.
     *    -> This loop handles the cases that were skipped in the previous
     *       loop.
     *    -> Note we use the problem description for the first material to get
     *       the maximum number of species. This should be changed to a value
     *       specific to the particular material type. However, since we are
     *       at the goma input file deck lvl, the first material will have to
     *       do.
     */
    if (retn == 1) {
      frac_vec = alloc_dbl_1(upd->Max_Num_Species, 0.0);
      retn =
          check_consistent_fraction_vector(Var_init, Num_Var_Init, upd->Max_Num_Species, frac_vec);

      /*
       *   Decide whether we need to do any conversions between the types of
       *   fraction vectors before insertion.
       */

      if (pd_glob[0]->Species_Var_Type != retn) {
        (void)convert_species_var(pd_glob[0]->Species_Var_Type, mp_glob[0], retn, frac_vec, 0.);
      }

      /*
       * Inject the value of the mass fraction vector, frac_vec, into all
       * parts of the solution vector, for all materials. Only need to
       * inject the maximum number of species equations.
       */
      for (i = 0; i < upd->Max_Num_Species_Eqn; i++) {
        init_vec_value(dum_var, frac_vec[i], DPI_ptr->num_owned_nodes);
        inject_nodal_vec(u, MASS_FRACTION, i, 0, -2, dum_var);
      }
      safer_free((void **)&frac_vec);
    }
    safer_free((void **)&dum_var);
  }

  /***********************************************************************/
  /*  PHASE 3 -> initializations of specific variable types              */
  /*             for a specific  material                                */
  /***********************************************************************/

  /*   arrange for element blocks to be processed in a specific order
   *   according to material number since this is important in order
   *   to maintain consistency between serial and parallel execution
   */

  block_order = alloc_int_1(exo->num_elem_blocks, 0);
  for (n = 0; n < exo->num_elem_blocks; n++)
    block_order[n] = n;

  if (Num_Proc > 1) {
    for (n = 0; n < exo->num_elem_blocks; n++) {
      for (k = n + 1; k < exo->num_elem_blocks; k++) {
        if (Matilda[block_order[k]] < Matilda[block_order[n]]) {
          block_temp = block_order[n];
          block_order[n] = block_order[k];
          block_order[k] = block_temp;
        }
      }
    }
  }

  /*
   * Loop over element blocks that exist for this processor, determine
   * which material corresponds to it.
   */

  for (ebj = 0; ebj < exo->num_elem_blocks; ebj++) {
    ebi = block_order[ebj];
    mn = Matilda[ebi];
    if (mn < 0) {
      continue;
    }
    pd = pd_glob[mn];
    mp = mp_glob[mn];
    if (Num_Var_Init_Mat[mn] > 0) {
      e_start = exo->eb_ptr[ebi];
      e_end = exo->eb_ptr[ebi + 1];

      for (j = 0; j < Num_Var_Init_Mat[mn]; j++) {
        if (Var_init_mat[mn][j].len_u_pars < 0) {
          DPRINTF(stdout, "\tSetting MAT %d %s number %d (variable [%d]) to %g\n", mn,
                  Var_Name[Var_init_mat[mn][j].var].name1, Var_init_mat[mn][j].ktype,
                  Var_init_mat[mn][j].var, Var_init_mat[mn][j].init_val);
        }
      }
      /*
       *  Loop over each element in the element block
       */
      for (ielem = e_start; ielem < e_end; ielem++) {
        ielem_type = Elem_Type(exo, ielem);
        num_nodes = elem_info(NNODES, ielem_type);
        index = Proc_Connect_Ptr[ielem];
        /*
         *  Loop over each local node in the element
         */
        for (n = 0; n < num_nodes; n++) {
          i = Proc_Elem_Connect[index++];
          nv = Nodes[i]->Nodal_Vars_Info[pg->imtrx];
          for (j = 0; j < Num_Var_Init_Mat[mn]; j++) {
            int slaved = Var_init_mat[mn][j].slave_block;

            var = Var_init_mat[mn][j].var;
            nunks = get_nv_ndofs_modMF(nv, var);

            if (slaved == 1)
              slaved = Nodes[i]->Mat_List.Length > 1;

            /*
             * We only want to initialize a variable at a node
             * if we have an unknown at that node, and there is
             * a valid interpolation for that variable in the
             * current element block
             */
            if (nunks > 0 && pd->i[pg->imtrx][var] && !slaved &&
                Var_init_mat[mn][j].len_u_pars < 0) {
              /*
               * Check against ktype here to make sure we have
               * an associated unknown
               */
              if (var == MASS_FRACTION && Var_init_mat[mn][j].ktype >= pd->Num_Species_Eqn) {
                continue;
              }
              if (nunks == 1) {
                ipos = Index_Solution(i, var, Var_init_mat[mn][j].ktype, 0, mn, pg->imtrx);
                u[ipos] = Var_init_mat[mn][j].init_val;
              } else {
                /*
                 * Ok, there is more than one degree of freedom for this
                 * variable type at this node. Why? Let's break down
                 * the reason
                 *  HKM -> we can't take this block out, until we
                 *         get rid of the debugging section below.
                 */
                interpType = pd->i[pg->imtrx][var];
                if (interpType == I_P0 || interpType == I_P1) {
                  /*
                   *  For P0 and P1 interpolation, only first dof is set
                   */
                  ipos = Index_Solution(i, var, Var_init_mat[mn][j].ktype, 0, mn, pg->imtrx);
                  u[ipos] = Var_init_mat[mn][j].init_val;
                } else if (interpType == I_PQ1 || interpType == I_PQ2) {
                  /*
                   * For linear and quadratic discontinuous interpolations
                   * all degrees of freedom in the element, which are all
                   * located at the centroid node, are set to the same
                   * value. Thus, find out how many degrees of freedom there
                   * are (also in the Variable_Description structure) and
                   * then set all of them to the initialization value
                   * set in the input deck.
                   */
                  ndof = 4;
                  if (interpType == I_PQ2)
                    ndof = 9;
                  for (k = 0; k < ndof; k++) {
                    ipos = Index_Solution(i, var, Var_init_mat[mn][j].ktype, k, mn, pg->imtrx);
                    u[ipos] = Var_init_mat[mn][j].init_val;
                  }
                } else {
                  /*
                   * For the case of discontinuous variables, we need to know
                   * which unknown to apply the condition to.
                   *
                   *   This currently is determined by a even odd
                   *   scheme wrt the element block id.
                   */
                  ipos = Index_Solution(i, var, Var_init_mat[mn][j].ktype, 0, mn, pg->imtrx);
                  u[ipos] = Var_init_mat[mn][j].init_val;
                }
              }
            }
          } /* end for j<Num_Var_Init_Mat[mn] */
        }   /* end for n<num_nodes */
      }     /* end for e_start=ielem<e_end */
    }       /* end if (Num_Var_Init_Mat[mn] > 0 */
  }
  /*  safer_free((void **) &block_order);  */
  /*
   *  Exchange the degrees of freedom with neighboring processors
   */
  exchange_dof(cx, dpi, u, pg->imtrx);

  /*
   * check for external variables which have to be read in
   * from other exodus files, or passed in from a host code.
   * or fields generated from pixel files
   */

  if (efv->ev) {
    k = 0;
    for (w = 0; w < efv->Num_external_field; w++) {
      if (efv->i[w] == I_TABLE) {
        ext_Tables[k] = setup_table_external(efv->file_nm[w], ext_table, efv->name[w]);
        k++;
      }

      /*
       * EDW NOTE: Using IMPORT or IMPORT_EV in the file_nm field indicates
       *           that the external var is imported into Goma from a
       *           calling program instead of being read in from a file.
       *           This task is handled elsewhere for these cases.
       */
      else if (strcmp(efv->file_nm[w], "IMPORT") != 0 &&
               strcmp(efv->file_nm[w], "IMPORT_EV") != 0) {
        if (!efv->ipix[w]) {
          DPRINTF(stdout, "%s:  reading fixed field \"%s\" from \"%s\" ...\n", yo, efv->name[w],
                  efv->file_nm[w]);
          err = rd_vectors_from_exoII(u, efv->file_nm[w], 1, w, INT_MAX, &timeValueReadExt, exo);
          if (err != 0) {
            DPRINTF(stderr, "%s:  err fr rd_vectors_from_exoII() while reading external fields\n",
                    yo);
          }
        }

        else if (efv->ipix[w] == 1) /*Original pixel to mesh tool */
        {
          DPRINTF(stderr, "\nMapping pixel image to mesh with original algorithm...");
          err = rd_image_to_mesh(w, exo);
        } else if (efv->ipix[w] == 2) /*'Fast' pixel to mesh tool */
        {
          DPRINTF(stderr, "\nMapping pixel image to mesh with 'fast' algorithm...");
          err = rd_image_to_mesh2(w, exo);
        } else {
          GOMA_EH(-1, "something wrong with efv->ipix");
        }
      }
#ifndef LIBRARY_MODE
      else if (strcmp(efv->file_nm[w], "IMPORT") == 0 ||
               strcmp(efv->file_nm[w], "IMPORT_EV") == 0) {
        GOMA_EH(-1, "External fields can only be imported in LIBRARY_MODE!");
      }
#endif
    }
  }

  /*
   *  Check to see if we are performing an TALE analysis, in which case we
   * need to save the displacement fields as they are read in,
   * as the reference state for KIN_DISPLACEMENT BCs
   *
   * Loop over element blocks that exist for this processor, determine
   * if any material is being treated with TALE
   */
  efv->TALE = FALSE;
  for (ebi = 0; ebi < exo->num_elem_blocks; ebi++) {
    mn = Matilda[ebi];
    if (mn < 0) {
      continue;
    }
    if (pd_glob[mn]->MeshMotion == TOTAL_ALE) {
      efv->TALE = TRUE;
    }
  }
  if (efv->TALE) {
    /*
     * If we are doing a TALE analysis,
     *   Allocate and save the initial displacement fields using
     *   fields in the external field variables struct.
     */
    for (a = 0; a < exo->num_dim; a++) {
      efv->init_displacement_ndl_val[a] = alloc_dbl_1(exo->num_nodes, 0.0);
      extract_nodal_vec(u, MESH_DISPLACEMENT1 + a, 0, -2, efv->init_displacement_ndl_val[a], exo,
                        FALSE, 0.);
    }
    for (a = 0; a < exo->num_dim; a++) {
      efv->init_displacement_ndl_val[a + DIM] = alloc_dbl_1(exo->num_nodes, 0.0);
      extract_nodal_vec(u, SOLID_DISPLACEMENT1 + a, 0, -2, efv->init_displacement_ndl_val[DIM + a],
                        exo, FALSE, 0.);
    }
  }
  exchange_dof(cx, dpi, u, pg->imtrx);

  /* User initialization part
   * Loop over element blocks that exist for this processor, determine
   * which material corresponds to it.
   */
  for (ebj = 0; ebj < exo->num_elem_blocks; ebj++) {
    ebi = block_order[ebj];
    mn = Matilda[ebi];
    if (mn < 0) {
      continue;
    }
    pd = pd_glob[mn];
    mp = mp_glob[mn];
    if (Num_Var_Init_Mat[mn] > 0) {
      e_start = exo->eb_ptr[ebi];
      e_end = exo->eb_ptr[ebi + 1];

      for (j = 0; j < Num_Var_Init_Mat[mn]; j++) {
        if (Var_init_mat[mn][j].len_u_pars > 0) {
          DPRINTF(stdout, "\tSetting MAT %d %s number %d (variable [%d]) to USER %g\n", mn,
                  Var_Name[Var_init_mat[mn][j].var].name1, Var_init_mat[mn][j].ktype,
                  Var_init_mat[mn][j].var, Var_init_mat[mn][j].init_val);
          DPRINTF(stdout, "\tUser parameters");
          for (i = 0; i < Var_init_mat[mn][j].len_u_pars; i++) {
            DPRINTF(stdout, "\t %g", Var_init_mat[mn][j].u_pars[i]);
          }
          DPRINTF(stdout, " \n");
        }
      }
      /*
       *  Loop over each element in the element block
       */
      for (ielem = e_start; ielem < e_end; ielem++) {
        ielem_type = Elem_Type(exo, ielem);
        num_nodes = elem_info(NNODES, ielem_type);
        index = Proc_Connect_Ptr[ielem];
        /*
         *  Loop over each local node in the element
         */
        for (n = 0; n < num_nodes; n++) {
          i = Proc_Elem_Connect[index++];
          nv = Nodes[i]->Nodal_Vars_Info[pg->imtrx];
          for (j = 0; j < Num_Var_Init_Mat[mn]; j++) {
            int slaved = Var_init_mat[mn][j].slave_block;

            var = Var_init_mat[mn][j].var;
            nunks = get_nv_ndofs_modMF(nv, var);

            if (slaved == 1)
              slaved = Nodes[i]->Mat_List.Length > 1;

            /*
             * We only want to initialize a variable at a node
             * if we have an unknown at that node, and there is
             * a valid interpolation for that variable in the
             * current element block
             */
            if (nunks > 0 && pd->i[pg->imtrx][var] && !slaved &&
                Var_init_mat[mn][j].len_u_pars >= 0) {
              double xpt[DIM] = {0, 0, 0}, var_val[MAX_VARIABLE_TYPES];
              int dir;
              for (dir = 0; dir < pd->Num_Dim; dir++) {
                xpt[dir] = Coor[dir][i];
              }
              init_vec_value(var_val, 0.0, MAX_VARIABLE_TYPES);
              for (a = 0; a < MAX_VARIABLE_TYPES; a++) {
                if (pd->v[pg->imtrx][a]) {
                  k = Index_Solution(i, a, 0, 0, mn, pg->imtrx);
                  if (k != -1) {
                    var_val[a] = u[k];
                  }
                }
              }
              /*
               * Check against ktype here to make sure we have
               * an associated unknown
               */
              if (var == MASS_FRACTION && Var_init_mat[mn][j].ktype >= pd->Num_Species_Eqn) {
                continue;
              }
              if (nunks == 1) {
                ipos = Index_Solution(i, var, Var_init_mat[mn][j].ktype, 0, mn, pg->imtrx);
                u[ipos] = user_mat_init(var, i, Var_init_mat[mn][j].init_val,
                                        Var_init_mat[mn][j].u_pars, xpt, mn, var_val);
              } else {
                /*
                 * Ok, there is more than one degree of freedom for this
                 * variable type at this node. Why? Let's break down
                 * the reason
                 *  HKM -> we can't take this block out, until we
                 *         get rid of the debugging section below.
                 */
                interpType = pd->i[pg->imtrx][var];
                if (interpType == I_P0 || interpType == I_P1) {
                  /*
                   *  For P0 and P1 interpolation, only first dof is set
                   */
                  ipos = Index_Solution(i, var, Var_init_mat[mn][j].ktype, 0, mn, pg->imtrx);
                  u[ipos] = user_mat_init(var, i, Var_init_mat[mn][j].init_val,
                                          Var_init_mat[mn][j].u_pars, xpt, mn, var_val);
                } else if (interpType == I_PQ1 || interpType == I_PQ2) {
                  /*
                   * For linear and quadratic discontinuous interpolations
                   * all degrees of freedom in the element, which are all
                   * located at the centroid node, are set to the same
                   * value. Thus, find out how many degrees of freedom there
                   * are (also in the Variable_Description structure) and
                   * then set all of them to the initialization value
                   * set in the input deck.
                   */
                  ndof = 4;
                  if (interpType == I_PQ2)
                    ndof = 9;
                  for (k = 0; k < ndof; k++) {
                    ipos = Index_Solution(i, var, Var_init_mat[mn][j].ktype, k, mn, pg->imtrx);
                    u[ipos] = user_mat_init(var, i, Var_init_mat[mn][j].init_val,
                                            Var_init_mat[mn][j].u_pars, xpt, mn, var_val);
                  }
                } else {
                  /*
                   * For the case of discontinuous variables, we need to know
                   * which unknown to apply the condition to.
                   *
                   *   This currently is determined by a even odd
                   *   scheme wrt the element block id.
                   */
                  ipos = Index_Solution(i, var, Var_init_mat[mn][j].ktype, 0, mn, pg->imtrx);
                  u[ipos] = user_mat_init(var, i, Var_init_mat[mn][j].init_val,
                                          Var_init_mat[mn][j].u_pars, xpt, mn, var_val);
                }
              }
            }
          } /* end for j<Num_Var_Init_Mat[mn] */
        }   /* end for n<num_nodes */
      }     /* end for e_start=ielem<e_end */
    }       /* end if (Num_Var_Init_Mat[mn] > 0 */
  }         /* end for element blocks */
  safer_free((void **)&block_order);

  /*
   *  Exchange the degrees of freedom with neighboring processors
   */
  exchange_dof(cx, dpi, u, pg->imtrx);

  return;
} /* END of routine init_vec **************************************************/
/******************************************************************************/
/******************************************************************************/

void init_structural_shell_coord(double u[])
/******************************************************************************
 *
 * init_structural_shell_coord(): Initialize the shell coordinates variables to
 *                                the mesh coordinates for initial guess. If you
 *                                don't do this than the displacement boundary
 *                                conditions applied to connect shell coords to
 *                                the moving mesh displacements will blow up
 *                                upon startup.
 *
 * Input
 * =====
 * u = Array of initial values for the indepenedent variables.
 *
 * Return
 * ======
 * void
 *
 ******************************************************************************/
{
  int node, comp, var, dof;

  for (comp = 0; comp < 2; comp++) /*Note that shells ONLY work in 2D so you should
                                    *never get here for 3D */
  {
    var = SHELL_X + comp;
    if (upd->vp[pg->imtrx][var]) {
      for (node = 0; node < Num_Node; node++) {
        if (Dolphin[pg->imtrx][node][var] > 0) {
          dof = Index_Solution(node, var, 0, 0, -1, pg->imtrx);
          u[dof] = Coor[comp][node];
        }
      } /* for: node(s) */
    }   /* if: var in upd-> */
  }     /* for: comp(onents) */

  return;

} /* End of init_structural_shell_coord() */

void init_shell_normal_unknowns(double x[], const Exo_DB *exo)
/******************************************************************************
 *
 * init_shell_normal_unknowns_3D(): Provides an initialization of the
 *                               shell normal vector unknowns. This is
 *                               necessary because if not provided, the
 *                               the accompanying curvature equation may
 *                               degenerate to a false trivial solution.
 *
 *                               The unknowns are initialized by looping
 *                               over all of elements containing shell
 *                               normal, setup shop at every node,
 *                               compute the normal vector, then copy it to
 *                               solution vectors
 *
/ * Input
* =====
* u = Array of initial values for the indepenedent variables.
* exo = Exodus database
*
* Return
* ======
* void
*
******************************************************************************/
{
  int e_start = 0, e_end = 0, ielem = 0;
  int ielem_type, ielem_dim, iconnect_ptr;
  int ilnode, ignode, num_local_nodes;
  int nxi, nyi, nzi;
  int i, err;
  dbl s, t, u, xi[DIM];

  for (i = 0; i < exo->num_elem_blocks; i++) {
    e_start = exo->eb_ptr[i];
    e_end = exo->eb_ptr[i + 1];
    ielem = exo->eb_ptr[i];
    ielem_type = Elem_Type(exo, ielem);
    ielem_dim = elem_info(NDIM, ielem_type);
    if (ielem_dim == pd->Num_Dim)
      continue;

    /* Loop over all elements */
    for (ielem = e_start; ielem < e_end; ielem++) {
      ielem_type = Elem_Type(exo, ielem);
      load_ei(ielem, exo, 0, pg->imtrx);
      err = load_elem_dofptr(ielem, exo, x, x, x, x, 0);
      GOMA_EH(err, "Can't load elem_dofptr in shell normals initialization");

      ielem_dim = elem_info(NDIM, ielem_type);
      num_local_nodes = elem_info(NNODES, ielem_type);
      iconnect_ptr = Proc_Connect_Ptr[ielem];

      switch (ei[pg->imtrx]->ielem_shape) {
      case LINE_SEGMENT:
      case SHELL:
      case TRISHELL:

        /* Loop over nodes within the element */
        for (ilnode = 0; ilnode < num_local_nodes; ilnode++) {
          /* Find s, t, u, coordinates of each node */
          find_nodal_stu(ilnode, ielem_type, &s, &t, &u);
          xi[0] = s;
          xi[1] = t;
          xi[2] = u;

          setup_shop_at_point(ielem, xi, exo);

          shell_determinant_and_normal(ielem, iconnect_ptr, num_local_nodes, ielem_dim, 1);

          /* Get global node number */
          ignode = Proc_Elem_Connect[iconnect_ptr + ilnode];

          /* Get node number and indices into solution vector for normal dofs */
          nxi = Index_Solution(ignode, SHELL_NORMAL1, 0, 0, -2, pg->imtrx);
          nyi = Index_Solution(ignode, SHELL_NORMAL2, 0, 0, -2, pg->imtrx);

          x[nxi] = fv->snormal[0];
          x[nyi] = fv->snormal[1];

          if (pd->Num_Dim == 3) {
            nzi = Index_Solution(ignode, SHELL_NORMAL3, 0, 0, -2, pg->imtrx);
            x[nzi] = fv->snormal[2];
          }
        }
        break;
      }
    }
  }
  return;
}

static void read_initial_guess(double u[], const int np, double uAC[], const int nAC)

/**************************************************************************
 *
 * read_initial_guess():
 *
 * Read initial guess from ASCII text neutral file (eg, "contin.dat")
 *
 * in:  np		total number of unknowns
 *     nAC          Number of Augmented conditions
 *
 * out:	u       Initial guess to solution vector.
 *         uAC      Initial guess for the value of the augmented
 *                  conditions
 **************************************************************************/
{
  int i, nchar;
  char input[MAX_COMMAND_LINE_LENGTH];
  FILE *file;
  char yo[] = "read_initial_guess";

  if (Debug_Flag > 0) {
    DPRINTF(stderr, "Initial guess read from  \"%s\" ...\n", Init_GuessFile);
  }

  /*
   *  Open the guess file for reading
   */
  file = fopen(Init_GuessFile, "r");
  if (file == NULL) {
    DPRINTF(stderr, "%s:  error opening file \"%s\" for reading\n", yo, Init_GuessFile);
    exit(-1);
  }

  /*
   *  Read the solution one line at a time.
   *  Each solution value is defined as the first number on each line.
   *  The rest of the line is ignored.
   */
  for (i = 0; i < np; i++) {
    nchar = read_line(file, input, FALSE);
    if (nchar <= 0) {
      fprintf(stderr, "%s: line %d of the initial guess file %s had an error, nchar = %d\n", yo, i,
              Init_GuessFile, nchar);
      fprintf(stderr, "%s:\t line = \"%s\"", yo, input);
      GOMA_EH(-1, yo);
    }
    if (!interpret_double(input, u + i)) {
      fprintf(stderr, "%s: line %d of the initial guess file %s had an error, %d\n", yo, i,
              Init_GuessFile, nchar);
      fprintf(stderr, "%s:\t line = \"%s\"", yo, input);
      GOMA_EH(-1, yo);
    }
  }

  /*
   *  Read the augmented conditions initial guess one line at a time
   *  Each augmented condition solution value is defined as the first number
   *  on each line. The rest of the line is ignored.
   */
  if (nAC > 0) {
    if (augc[0].iread == 1) {
      for (i = 0; i < nAC; i++) {
        nchar = read_line(file, input, FALSE);
        if (nchar <= 0) {
          fprintf(stderr, "%s: Aborting(nchar) read of Aug var %d from file %s , %s\n", yo, i,
                  Init_GuessFile, input);
          break;
        }
        if (!interpret_double(input, uAC + i)) {
          fprintf(stderr, "%s: Aborting(dble) read of Aug var %d from file %s , %s\n", yo, i,
                  Init_GuessFile, input);
          break;
        }
      }
    }
  }
  if (Debug_Flag > 0) {
    DPRINTF(stderr, "%s:  Successfully read the initial guess from \"%s\" ...\n", yo,
            Init_GuessFile);
  }
  fclose(file);
} /* END of routine read_initial_guess()  */
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int write_ascii_soln(double *u,     /* Solution vector */
                     double *resid, /* Residual vector (optional */
                     int np,        /* Number of unknowns on the processor */
                     double *uAC,   /* Value of augmented variables */
                     int nAC,       /* Number of augmented conditions */
                     double time,   /* Current time */
                     FILE *file)    /* FILE pointer, pointing to an open
                                     * file */
/*************************************************************************
 *
 *  write_ascii_soln()
 *
 *       Write solution to an ascii output file
 *
 *************************************************************************/
{
  int i;
  if (file == NULL)
    return 0;
  if (TimeIntegration != STEADY)
    fprintf(file, "time = %f\n\n", time);
  if (Continuation != ALC_NONE)
    fprintf(file, "path = %f\n\n", time);
  if (dofname != NULL) {
    if (resid == NULL) {
      for (i = 0; i < np; i++) {
        fprintf(file, "%23.16e %s\n", u[i], dofname[pg->imtrx][i]);
      }
    } else {
      for (i = 0; i < np; i++) {
        fprintf(file, "%23.16e %s %23.16e\n", u[i], dofname[pg->imtrx][i], resid[i]);
      }
    }
  } else {
    if (resid == NULL) {
      for (i = 0; i < np; i++)
        fprintf(file, "%23.16e\n", u[i]);
    } else {
      for (i = 0; i < np; i++)
        fprintf(file, "%23.16e %23.16e\n", u[i], resid[i]);
    }
  }
  /**  print out augmenting variables **/
  if (dofname != NULL) {
    for (i = 0; i < nAC; i++) {
      fprintf(file, "%23.16e Aug_Cond=%d\n", uAC[i], i);
    }
  } else {
    for (i = 0; i < nAC; i++) {
      fprintf(file, "%23.16e %10d\n", uAC[i], i);
    }
  }

  /*
   *  Add a statement to the output file indicating that the current
   *  solution output is done
   */
  fprintf(file, " oy !\n");
  return (0);
}
/******************************************************************************/

/* wr_soln_vec -- write out solution vector to a specified file
 *
 * Notes:	This generalizes the write_initial_guess routine by
 *		permitting one to write the whole vector to a file at
 *		some intermediate Newton iteration.
 *
 *		Write both the solution and the dof-map and residuals.
 *
 * Created:	Tue Apr  5 11:42:23 MDT 1994 pasacki@sandia.gov
 *
 * Modified:	Mon Mar 27 08:56 MST 1995 pasacki@sandia.gov
 */

int wr_soln_vec(double u[],    /* solution vector */
                double r[],    /* residual vector */
                const int np,  /* number of elements in solution vector */
                const int itn) /* Newton iteration we are at  */
{
  int i, status;
  static char fmt[] = "%23.16e %.22s %23.16e\n";
  char fname[MAX_FNL];
  FILE *file;

  status = 0;

  sprintf(fname, "tmp_%d.d", itn);

#ifdef PARALLEL
  multiname(fname, ProcID, Num_Proc);
#endif

  file = fopen(fname, "w");

  if (file == NULL) {
    GOMA_EH(-1, "Problem opening intermediate results file.");
  }

  for (i = 0; i < np; i++) {
    fprintf(file, fmt, u[i], dofname[pg->imtrx][i], r[i]);
  }

  fclose(file);

  return (status);
} /* END of routine wr_soln_vec  */
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

int rd_vectors_from_exoII(double u[],
                          const char *file_nm,
                          const int action_flag,
                          const int variable_no,
                          const int desired_time_step,
                          double *timeValueRead,
                          const Exo_DB *exo)

/*******************************************************************
 *
 * rd_vectors_from_exoII:
 *
 * Read initial guess from exoIIv2 file
 *
 * in:  file_nm     = String name of the exodus file
 *      action_flag = 0 read initial guess for problem
 *                  = 1 read extern auxillary fixed variables
 *                  = 2 read element variables into elem_storage_struct
 *	    variable_no =  Used only when action_flag = 1
 *                     Specifies the number of the external
 *                     variable to be read
 *                     (basically the card number in order)
 *      desired_time_step = time step number to read from
 *
 * out:	u	        initial guess to solution vector.
 *******************************************************************/
{
  int i, error, vdex, num_dim, num_nodes, mn, icount;
  int num_elem, num_elem_blk, num_node_sets, num_side_sets, time_step;
  float version;               /* version number of EXODUS II */
  int exoid;                   /* ID of the open EXODUS II file */
  char title[MAX_LINE_LENGTH]; /* title of the EXODUS II database */
  float ret_float;             /* any returned float */
  char ret_char[3];            /* any returned character */
  int num_vars;                /* number of var_type variables */
  char **var_names = NULL;     /* array containing num_vars variable names */
  int num_elem_vars = 0;
  char **elem_var_names = NULL; /* array containing element variable names */
  int w;                        /* counter for species concentration */
  int var;
  MATRL_PROP_STRUCT *matrl = 0;
  double ftimeValue;
#ifdef DEBUG
  static const char yo[] = "rd_vectors_from_exoII";
#endif

  CPU_word_size = sizeof(double);
  IO_word_size = 0;

  exoid = ex_open(file_nm, EX_READ, &CPU_word_size, &IO_word_size, &version);
  GOMA_EH(exoid, "ex_open");

  error = ex_get_init(exoid, title, &num_dim, &num_nodes, &num_elem, &num_elem_blk, &num_node_sets,
                      &num_side_sets);
  GOMA_EH(error, "ex_get_init for efv or init guess");

  // add some checks to make sure we are probably reading the same mesh
  GOMA_ASSERT_ALWAYS(num_dim == exo->base_mesh->num_dim);
  GOMA_ASSERT_ALWAYS(num_nodes == exo->base_mesh->num_nodes);
  GOMA_ASSERT_ALWAYS(num_elem == exo->base_mesh->num_elems);
  GOMA_ASSERT_ALWAYS(num_elem_blk == exo->base_mesh->num_elem_blocks);

  /*
   * Obtain the number of time steps in the exodus file, time_step,
   * We will read only from the last time step
   */
  error = ex_inquire(exoid, EX_INQ_TIME, &time_step, &ret_float, ret_char);
  GOMA_EH(error, "ex_inquire");

  if (time_step == 0) {
    // early exit
    GOMA_WH(GOMA_ERROR, "Warning no time steps found in %s", file_nm);
    ex_close(exoid);
    return 0;
  }

  /* Figure out what time step to select. Will select the last time
   * step unless the input variable desired_time_step is set lower.
   * The lower limit in exodus is 1.
   */
  if (desired_time_step < time_step) {
    time_step = MAX(1, desired_time_step);
  }

  // Return the value of the time
  error = ex_get_time(exoid, time_step, &ftimeValue);
  if (error == -1) {
    ftimeValue = 0.0;
  }
  if (timeValueRead) {
    *timeValueRead = ftimeValue;
  }

  /* Based on problem type and available info in database, extract
   * appropriate fields
   */

  /*
   * Get the number of nodal variables in the file, and allocate
   * space for storage of their names.
   */
  error = ex_get_variable_param(exoid, EX_NODAL, &num_vars);
  GOMA_EH(error, "ex_get_variable_param nodal");
  error = ex_get_variable_param(exoid, EX_ELEM_BLOCK, &num_elem_vars);
  GOMA_EH(error, "ex_get_var_param elem");

  /* First extract all nodal variable names in exoII database */
  if (num_vars > 0) {
    var_names = alloc_VecFixedStrings(num_vars, (MAX_STR_LENGTH + 1));
    error = ex_get_variable_names(exoid, EX_NODAL, num_vars, var_names);
    GOMA_EH(error, "ex_get_variable_names nodal");
    for (i = 0; i < num_vars; i++)
      strip(var_names[i]);
  } else {
    fprintf(stderr, "Warning: no nodal variables stored in exoII input file.\n");
  }

  /* First extract all element variable names in exoII database */
  if (num_elem_vars > 0) {
    elem_var_names = alloc_VecFixedStrings(num_elem_vars, (MAX_STR_LENGTH + 1));
    error = ex_get_variable_names(exoid, EX_ELEM_BLOCK, num_elem_vars, elem_var_names);
    GOMA_EH(error, "ex_get_variable_names element");
    for (i = 0; i < num_elem_vars; i++)
      strip(elem_var_names[i]);
  }

  /* If action_flag is 0,
   * Cross check problem type (variables requested)
   * and variables available
   */

  if (action_flag == 0) {
    for (var = V_FIRST; var < V_LAST; var++) {
      icount = 0;
      if (Num_Var_In_Type[pg->imtrx][var]) {
        if (var == MASS_FRACTION) {
          for (mn = -1; mn < upd->Num_Mat; mn++) {
            if (mn == -1) {
              for (i = upd->Num_Mat - 1; i >= 0; i--) {
                if (mp_glob[i]->Num_Species == upd->Max_Num_Species) {
                  matrl = mp_glob[i];
                }
              }
            } else {
              matrl = mp_glob[mn];
            }
            for (w = 0; w < matrl->Num_Species_Eqn; w++) {
              error = rd_exoII_nv(u, var, mn, matrl, var_names, num_nodes, num_vars, exoid,
                                  time_step, w);
              if (!error)
                icount++;
            }
          }
        } else {
          for (mn = -1; mn < upd->Num_Mat; mn++) {
            if (mn == -1) {
              for (i = upd->Num_Mat - 1; i >= 0; i--) {
                if (pd_glob[i]->i[pg->imtrx][var]) {
                  matrl = mp_glob[i];
                }
              }
            } else {
              matrl = mp_glob[mn];
            }
            if (mn != -1 && (pd_glob[mn]->i[pg->imtrx][var] == I_P0)) {
              int eb_index = in_list(mn, 0, exo->num_elem_blocks, Matilda);
              if (eb_index != -1 && exo->base_mesh->eb_num_elems[eb_index] > 0) {
                error = rd_exoII_ev(u, var, mn, matrl, elem_var_names,
                                    exo->base_mesh->eb_num_elems[eb_index], num_elem_vars, exoid,
                                    time_step, 0, exo);
              }
            } else if (mn != -1 && (pd_glob[mn]->i[pg->imtrx][var] == I_P1)) {
              int eb_index = in_list(mn, 0, exo->num_elem_blocks, Matilda);
              if (eb_index != -1 && exo->base_mesh->eb_num_elems[eb_index] > 0) {
                int dof = getdofs(type2shape(exo->eb_elem_itype[eb_index]), I_P1);
                for (int i = 0; i < dof; i++) {
                  error = rd_exoII_ev(u, var, mn, matrl, elem_var_names,
                                      exo->base_mesh->eb_num_elems[eb_index], num_elem_vars, exoid,
                                      time_step, i, exo);
                }
              }
            } else {
              error = rd_exoII_nv(u, var, mn, matrl, var_names, num_nodes, num_vars, exoid,
                                  time_step, 0);
            }
            if (!error)
              icount++;
          }
        }
      }
    }
  }

  if (action_flag == 1) {
    if (efv->ev) {
      /*
       * Allocate memory for external field variable arrays
       */
      efv->ext_fld_ndl_val[variable_no] = alloc_dbl_1(exo->num_nodes, 0.0);
      printf("rd_vectors_from_exoII: Allocated field %d for %s at %p\n", variable_no,
             efv->name[variable_no], (void *)efv->ext_fld_ndl_val[variable_no]);
      vdex = -1;
#ifdef REACTION_PRODUCT_EFV
      if (TimeIntegration != STEADY) {
        efv->ext_fld_ndl_val_old[variable_no] = alloc_dbl_1(exo->num_nodes, 0.);
        efv->ext_fld_ndl_val_older[variable_no] = alloc_dbl_1(exo->num_nodes, 0.);
      }
#endif
      for (i = 0; i < num_vars; i++) {
        if (strcmp(var_names[i], efv->name[variable_no]) == 0) {
          vdex = i + 1;
        }
      }
      if (vdex == -1) {
        DPRINTF(stdout, "\n Cannot find external fields in exoII database, setting to null");
      } else {
        dbl *tmp_vector = alloc_dbl_1(num_nodes, 0.0);
        error = ex_get_var(exoid, time_step, EX_NODAL, vdex, 1, num_nodes, tmp_vector);
        for (int k = 0; k < exo->num_nodes; k++) {
          int base_index = exo->ghost_node_to_base[k];
          if (base_index != -1) {
            efv->ext_fld_ndl_val[variable_no][k] = tmp_vector[base_index];
          }
        }

        exchange_node(cx[0], DPI_ptr, efv->ext_fld_ndl_val[variable_no]);
        free(tmp_vector);

        GOMA_EH(error, "ex_get_var nodal");
      }
    }
  }

  safer_free((void **)&var_names);
  safer_free((void **)&elem_var_names);
  error = ex_close(exoid);
  GOMA_EH(error, "ex_close");
  return 0;
} /* end rd_vectors_from_exoII*/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int rd_trans_vectors_from_exoII(double u[],
                                const char *file_nm,
                                const int variable_no,
                                const int desired_time_step,
                                double *timeValueIn,
                                Exo_DB *exo,
                                Comm_Ex *cx,
                                Dpi *dpi)

/*******************************************************************
 *
 * rd_vectors_from_exoII:
 *
 * Read time dependent extenal fields from exoIIv2 file
 *
 * in:  file_nm     = String name of the exodus file
 *	    variable_no =  Used only when action_flag = 1
 *                     Specifies the number of the external
 *                     variable to be read
 *                     (basically the card number in order)
 *      desired_time_step = time step number to read from
 *      timeValueIn = current time value of the solution process
 *                     time value to be matched in the external field
 *
 * out:	u	        initial guess to solution vector.
 *******************************************************************/
{
  int i, k, error, vdex, num_dim, num_nodes;
  int num_elem, num_elem_blk, num_node_sets, num_side_sets, time_step;
  float version;               /* version number of EXODUS II */
  int exoid;                   /* ID of the open EXODUS II file */
  char title[MAX_LINE_LENGTH]; /* title of the EXODUS II database */
  float ret_float;             /* any returned float */
  char ret_char[3];            /* any returned character */
  int num_vars;                /* number of var_type variables */
  char **var_names = NULL;     /* array containing num_vars variable names */
  int num_elem_vars = 0;
  double ftimeValue, time_higher, time_lower;
  double *val_low, *val_high, slope, yint;
  int time_step_read, time_step_higher, time_step_lower, time_step_max;
#ifdef DEBUG
  static const char yo[] = "rd_trans_vectors_from_exoII";
#endif

  CPU_word_size = sizeof(double);
  IO_word_size = 0;

  exoid = ex_open(file_nm, EX_READ, &CPU_word_size, &IO_word_size, &version);
  GOMA_EH(exoid, "ex_open");

  error = ex_get_init(exoid, title, &num_dim, &num_nodes, &num_elem, &num_elem_blk, &num_node_sets,
                      &num_side_sets);
  GOMA_EH(error, "ex_get_init for efv or init guess");

  /*
   * Obtain the number of time steps in the exodus file, time_step,
   */
  error = ex_inquire(exoid, EX_INQ_TIME, &time_step, &ret_float, ret_char);
  GOMA_EH(error, "ex_inquire");

  /* Figure out what time step to select. Will select the closest time
   * step to the current time in the solution procedure.
   * The lower limit in exodus is 1.
   */
  time_step_max = time_step;
  time_step_read = 1;
  time_step_lower = time_step_read;
  time_lower = 0.0;

  error = ex_get_time(exoid, time_step_read, &ftimeValue);
  if (error == -1) {
    ftimeValue = 0.0;
  }

  /* Cycle through the external time values, stop when the time value desired
   * has been passed.
   */
  while ((*timeValueIn > ftimeValue) && (time_step_read < time_step_max)) {
    time_step_lower = time_step_read;
    time_lower = ftimeValue;
    time_step_read++;
    ex_get_time(exoid, time_step_read, &ftimeValue);
  }

  /*if ((*timeValueIn <= ftimeValue) || (time_step_read >= time_step_max))
    {
      time_higher = ftimeValue;
    }
   */
  /* Calculate which external time plane is closest to the desired time */
  time_step_higher = time_step_read;
  time_higher = ftimeValue;

  // Return the value of the time
  /*error = ex_get_time(exoid, time_step, &ftimeValue);
  if (error == -1) {
    ftimeValue = 0.0;
  }
  if (timeValueIn) {
    *timeValueIn = ftimeValue;
  } */

  // Open file to write results for debug info
  /*ofp = fopen("Trans_EFV", "a");
  fprintf(ofp,"Time value passed to function %e\n", *timeValueIn);
  fprintf(ofp, "time_step_higher %d\n", time_step_higher);
  fprintf(ofp, "time_step_lower %d\n", time_step_lower);
  fprintf(ofp, "time_higher %e\n", time_higher);
  fprintf(ofp, "time_lower %e\n", time_lower);
  //fclose(ofp);*/

  /*
   * Get the number of nodal variables in the file, and allocate
   * space for storage of their names.
   */
  error = ex_get_variable_param(exoid, EX_NODAL, &num_vars);
  GOMA_EH(error, "ex_get_variable_param nodal");
  error = ex_get_variable_param(exoid, EX_ELEM_BLOCK, &num_elem_vars);
  GOMA_EH(error, "ex_get_variable_param elem");

  /* First extract all nodal variable names in exoII database */
  if (num_vars > 0) {
    var_names = alloc_VecFixedStrings(num_vars, (MAX_STR_LENGTH + 1));
    error = ex_get_variable_names(exoid, EX_NODAL, num_vars, var_names);
    GOMA_EH(error, "ex_get_variable_names nodal");
    for (i = 0; i < num_vars; i++)
      strip(var_names[i]);
  } else {
    fprintf(stderr, "Warning: no nodal variables stored in exoII input file.\n");
  }

  if (efv->ev) {
    /*
     * Already allocated this in rd_vectors_from_exoII for external field variable arrays
     */
    // efv->ext_fld_ndl_val[variable_no] = alloc_dbl_1(num_nodes, 0.0);
    val_low = alloc_dbl_1(num_nodes, 0.0);
    val_high = alloc_dbl_1(num_nodes, 0.0);

    if (desired_time_step == 0) {
      fprintf(stderr, "rd_trans_vectors_from_exoII: Into existing field %d for %s at %p\n",
              variable_no, efv->name[variable_no], (void *)efv->ext_fld_ndl_val[variable_no]);
    }
    vdex = -1;
    for (i = 0; i < num_vars; i++) {
      if (strcmp(var_names[i], efv->name[variable_no]) == 0) {
        vdex = i + 1;
      }
    }
    if (vdex == -1) {
      DPRINTF(stdout, "\n Cannot find external fields in exoII database, setting to null");
    } else {
      error = ex_get_var(exoid, time_step_lower, EX_NODAL, vdex, 1, num_nodes, val_low);
      error = ex_get_var(exoid, time_step_higher, EX_NODAL, vdex, 1, num_nodes, val_high);
      GOMA_EH(error, "ex_get_var nodal");

      for (k = 0; k < exo->num_nodes; k++) {
        int base_index = exo->ghost_node_to_base[k];
        if (base_index != -1) {
          slope = (val_high[base_index] - val_low[base_index]) / (time_higher - time_lower);
          yint = val_low[base_index] - slope * time_lower;
          efv->ext_fld_ndl_val[variable_no][k] = slope * (*timeValueIn) + yint;
        }
      }
    }
  }

  safer_free((void **)&var_names);
  error = ex_close(exoid);
  GOMA_EH(error, "ex_close");
  // fclose(ofp);

  /*
   *  Exchange the degrees of freedom with neighboring processors
   */
  exchange_node(cx, dpi, efv->ext_fld_ndl_val[variable_no]);

  return 0;
} /* end rd_trans_vectors_from_exoII*/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int rd_exoII_nv(double *u,
                int varType,
                int mn,
                MATRL_PROP_STRUCT *matrl,
                char **var_names,
                int num_nodes,
                int num_vars,
                int exoII_id,
                int time_step,
                int spec)

/*************************************************************************
 *
 * rd_exoII_nv():
 *
 *    Reads an exodus nodal variable matched by the name, Var_exoII,
 * into the goma solution vector, u.
 * Because Exodus has no way to specify nodal variables by material
 * type, this routine will initialize all variables of the variable type,
 * Var_exoII->Index at the node irrespective of what material they
 * are in.
 *************************************************************************/
{
  int vdex = -1, i, error, status = 0;
  char exo_var_name[256], exo_var_desc[256];
  double *variable = NULL;
  assign_var_name(varType, spec, matrl, exo_var_name, exo_var_desc, mn);
  for (i = 0; i < num_vars; i++) {
    if (strcmp(var_names[i], exo_var_name) == 0) {
      vdex = i + 1;
    }
  }
  if (vdex != -1) {
    variable = alloc_dbl_1(num_nodes, 0.0);
    status = vdex;
    DPRINTF(stdout, "Nodal variable %s found in exoII database - reading.\n", exo_var_name);
    error = ex_get_var(exoII_id, time_step, EX_NODAL, vdex, 1, num_nodes, variable);
    GOMA_EH(error, "ex_get_var nodal");
    inject_nodal_vec(u, varType, spec, 0, mn, variable);
    safer_free((void **)&variable);
  }
  return status;
} /* END of routine rd_exoII_util */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int rd_exoII_ev(double *u,
                int varType,
                int mn,
                MATRL_PROP_STRUCT *matrl,
                char **elem_var_names,
                int num_elems_block,
                int num_elem_vars,
                int exoII_id,
                int time_step,
                int spec,
                const Exo_DB *exo)

/*************************************************************************
 *
 * rd_exoII_ev():
 *
 *    Reads an exodus element variable matched by the name, Var_exoII,
 * into the goma solution vector, u.
 * Because Exodus has no way to specify element variables by material
 * type, this routine will initialize all variables of the variable type,
 * Var_exoII->Index at the node irrespective of what material they
 * are in.
 *************************************************************************/
{
  int vdex = -1, i, error, status = 0;
  char exo_var_name[256], exo_var_desc[256];
  double *variable = NULL;
  assign_var_name(varType, spec, matrl, exo_var_name, exo_var_desc, mn);
  for (i = 0; i < num_elem_vars; i++) {
    if (strcmp(elem_var_names[i], exo_var_name) == 0) {
      vdex = i + 1;
    }
  }
  if (vdex != -1) {
    variable =
        alloc_dbl_1(num_elems_block, 0.0); // This should be at for number of elements in a block.
    status = vdex;
    DPRINTF(stdout, "Element variable %s for material %d found in exoII database - reading.\n",
            exo_var_name, mn + 1);
    error = ex_get_var(exoII_id, time_step, EX_ELEM_BLOCK, vdex, mn + 1, num_elems_block, variable);
    GOMA_EH(error, "ex_get_var element");
    inject_elem_vec(u, varType, 0, spec, mn, variable, exo, num_elems_block);
    safer_free((void **)&variable);
  }
  return status;
} /* END of routine rd_exoII_util */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void fill_dvec_rand(double u[], /* vector to fill with doubles in [0,1] */
                    int len)    /* length of vector */
/*
 *        Set the vector u to a random vector.
 *
 *        Author:         John N. Shadid (1421)
 *        Date:           1/21/93
 *        Revised:	 Tue Mar 22 12:22:38 MST 1994 pasacki@sandia.gov
 *
 *
 * Notes: Revised to use the nice drand48() stdlib function and be useful
 *	  for any vector of doubles - no more hardwired implicit assumptions.
 *
 *	  Revised: 1997/08/23 14:47 MDT pasacki@sandia.gov
 */

{
  int i;
  long seedval;

  /*
   * Initialize the seed for drand48()...
   */

  seedval = (long)ut();

  srand48(seedval);

  for (i = 0; i < len; i++) {
    u[i] = drand48();
  }
  return;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

int rd_globals_from_exoII(double u[], const char *file_nm, const int start, const int n)

/*******************************************************************
 *
 * rd_globals_from_exoII:
 *
 * Read initial guess for global variables from exoIIv2 file
 *  in: file_nm            filename  of exoII database file
 *      start              index of starting global variable to be copied to u
 *      n                  read no more than n global variables into u
 *  out:
 *      u                  initial guess for global variables.  Can be subset of
 *                         all global variables present ( start ! = 0 )
 *
 *******************************************************************/
{
  int i, error, num_dim, num_nodes;
  int num_elem, num_elem_blk, num_node_sets, num_side_sets, time_step;
  float version;               /* version number of EXODUS II */
  int exoid;                   /* ID of the open EXODUS II file */
  char title[MAX_LINE_LENGTH]; /* title of the EXODUS II database */
  float ret_float;             /* any returned float */
  char ret_char[3];            /* any returned character */
  int num_global_vars = -1;    /* number of global variables present */
  double *global_vars;         /* global variable values */
#ifdef DEBUG
  static const char yo[] = "rd_globals_from_exoII";
#endif

  CPU_word_size = sizeof(double);
  IO_word_size = 0;

  exoid = ex_open(file_nm, EX_READ, &CPU_word_size, &IO_word_size, &version);
  GOMA_EH(exoid, "ex_open");

  error = ex_get_init(exoid, title, &num_dim, &num_nodes, &num_elem, &num_elem_blk, &num_node_sets,
                      &num_side_sets);
  GOMA_EH(error, "ex_get_init for efv or init guess");

  error = ex_get_variable_param(exoid, EX_GLOBAL, &num_global_vars);

  if (num_global_vars >
      start) /* This catches both no global vars and no augmenting values to be read */
  {

    global_vars = alloc_dbl_1(num_global_vars, 0.0);

    /*
     * Obtain the number of time steps in the exodus file, time_step,
     * We will read only from the last time step
     */
    error = ex_inquire(exoid, EX_INQ_TIME, &time_step, &ret_float, ret_char);
    GOMA_EH(error, "ex_inquire");

    error = ex_get_var(exoid, time_step, EX_GLOBAL, 1, 1, num_global_vars, global_vars);
    GOMA_EH(error, "ex_get_var global");

    /*
     *  Read only the global vars that are there.  No more, no less
     */

    if (start + n < num_global_vars)
      num_global_vars = start + n;

    for (i = 0; (start + i) < num_global_vars; i++)
      u[i] = global_vars[start + i];

    safer_free((void **)&global_vars);
  }
  error = ex_close(exoid);
  GOMA_EH(error, "ex_close");
  return num_global_vars;
}

/************************************************************************/
/************************************************************************/
/************************************************************************/

static void inject_nodal_vec(double sol_vec[],
                             const int varType,
                             const int k,
                             const int idof,
                             const int matID,
                             const double nodal_vec[])

/**********************************************************************
 *
 * inject_nodal_vec:
 *
 * Routine to scatter a variable vector "k" of variable type
 * "var_no" to the solution vector "sol_vec".
 *
 *
 *     Now, for distributed processing, this routine only scatters a
 *     nodal variable into the full dof vector for the nodes that this
 *     processor owns. Forget the "global" part.
 *
 * It first tries to use the generic material unknown. Then, it
 * tries to use the unknown corresponding to the first material
 * index at that node.
 *
 * Author: 	P. R. Schunk (1511, SNL)
 * Date:		10/12/94
 * Revised: 1997/08/27 10:31 MDT pasacki@sandia.gov
 *
 * Parameter List:
 *
 * nodal_vec[] == 	vector containing the solution variable type
 *                    vector index by the local node number on the
 *         		current processor.
 *
 * var_no      ==  	integer variable type which defines
 *		        what variable is to be scattered into the
 *                    solution vector. Additionally, only the
 *                    generic material ID part of the solution
 *                    vector will be overwritten.
 *		        (see rf_fem_const.h - Variable Names)
 *
 *	k        == 	kth subvariable of type "var_no", k = 0 is first
 *		       variable of this type. Subvariables are only
 *                    used for MASS_FRACTION variable types.
 *  idof       ==     Degree of freedom id for this variable. Usually,
 *                    this is equal to zero. However, there are
 *                    some cases of multiple degrees of freedom
 *                    for a variable type at a node.
 * matID       ==     Material ID for the variable to scatter to.
 *                    -2: special value means to scatter to all
 *                        variables of type VariableType no matter
 *                        what material they refer to.
 *                    -1: This number refers to the
 *                        non-specific-to-a-material variable at
 *                        the node.
 *                    >=0: Variable which is specific to the
 *                         material matID.
 *
 * nodal_vec   ==	Processor solution vector
 *******************************************************************/
{
  int index, i, j, vindex;
  NODAL_VARS_STRUCT *nv;
  VARIABLE_DESCRIPTION_STRUCT *vd;
  for (i = 0; i < DPI_ptr->num_universe_nodes; i++) {
    nv = Nodes[i]->Nodal_Vars_Info[pg->imtrx];
    int base_index = EXO_ptr->ghost_node_to_base[i];
    if (base_index != -1) {
      if (matID == -2) {
        for (j = 0; j < (int)nv->Num_Var_Desc_Per_Type[varType]; j++) {
          vindex = nv->Var_Type_Index[varType][j];
          vd = nv->Var_Desc_List[vindex];
#ifdef DEBUG_HKM
          if (idof < 0 || idof > vd->Ndof) {
            fprintf(stderr, "init_vec ERROR: bad idof\n");
            GOMA_EH(-1, "init_vec bad idof");
          }
#endif
          if (k == (int)vd->Subvar_Index) {
            index = (Nodes[i]->First_Unknown[pg->imtrx] + nv->Nodal_Offset[vindex] + idof);
            sol_vec[index] = nodal_vec[base_index];
          }
        }
      } else {
        int ndof = 0;
        if (pd->i[pg->imtrx][varType] == I_PQ1) {
          ndof = 4;
          GOMA_WH(-1, "Using 4 dof injecting variable");
        } else if (pd->i[pg->imtrx][varType] == I_PQ2) {
          ndof = 9;
          GOMA_WH(-1, "Using 9 dof injecting variable");
        }

        if (ndof > 0) {
          int local_dof = 0;
          for (local_dof = 0; local_dof < ndof; local_dof++) {
            index = Index_Solution(i, varType, k, local_dof, matID, pg->imtrx);
            if (index != -1) {
              sol_vec[index] = nodal_vec[base_index];
            }
          }
        } else {
          index = Index_Solution(i, varType, k, idof, matID, pg->imtrx);
          if (index != -1) {
            sol_vec[index] = nodal_vec[base_index];
          }
        }
      }
    }
  }
  return;
}

/************************************************************************/
/************************************************************************/
/************************************************************************/

static void inject_elem_vec(double sol_vec[],
                            const int varType,
                            const int k,
                            const int idof,
                            const int matID,
                            const double elem_vec[],
                            const Exo_DB *exo,
                            const int num_elems_blk)

/**********************************************************************
 *
 * inject_elem_vec:
 *
 * Routine to scatter a variable vector "k" of variable type
 * "var_no" to the solution vector "sol_vec".
 *
 *
 *     Now, for distributed processing, this routine only scatters a
 *     element variable into the full dof vector for the element blocks that this
 *     processor owns. Forget the "global" part.
 *
 * It first tries to use the generic material unknown. Then, it
 * tries to use the unknown corresponding to the first material
 * index at that node.
 *
 * Author: 	P. R. Schunk (1511, SNL)
 * Date:		10/12/94
 * Revised: 1997/08/27 10:31 MDT pasacki@sandia.gov
 *
 * Parameter List:
 *
 * elem_vec[] == 	vector containing the solution variable type
 *                    vector index by the local node number on the
 *         		current processor.
 *
 * var_no      ==  	integer variable type which defines
 *		        what variable is to be scattered into the
 *                    solution vector. Additionally, only the
 *                    generic material ID part of the solution
 *                    vector will be overwritten.
 *		        (see rf_fem_const.h - Variable Names)
 *
 *	k        == 	kth subvariable of type "var_no", k = 0 is first
 *		       variable of this type. Subvariables are only
 *                    used for MASS_FRACTION variable types.
 *  idof       ==     Degree of freedom id for this variable. Usually,
 *                    this is equal to zero. However, there are
 *                    some cases of multiple degrees of freedom
 *                    for a variable type at a node.
 * matID       ==     Material ID for the variable to scatter to.
 *                    -2: special value means to scatter to all
 *                        variables of type VariableType no matter
 *                        what material they refer to.
 *                    -1: This number refers to the
 *                        non-specific-to-a-material variable at
 *                        the node.
 *                    >=0: Variable which is specific to the
 *                         material matID.
 *
 * elem_vec   ==	Processor solution vector
 *******************************************************************/
{
  int e_start, e_end, ielem, ielem_type, num_local_nodes;
  int iconnect_ptr, i, I, index;
  int found_quantity;
  int eb_index = in_list(matID, 0, exo->num_elem_blocks, Matilda);
  GOMA_EH(eb_index, "Trying to read unknown material element block index inject_elem_vec");
  e_start = exo->eb_ptr[eb_index];
  e_end = exo->eb_ptr[eb_index + 1];
  for (ielem = e_start; ielem < e_end; ielem++) {

    ielem_type = Elem_Type(exo, ielem); /* func defd in el_geom.h */
    num_local_nodes = elem_info(NNODES, ielem_type);
    iconnect_ptr = Proc_Connect_Ptr[ielem]; /* find ptr to beginning */
    /* of this element's */
    /* connectivity list */

    /* We're looking at a nodal quantity that should be defined at only 1
       node of this element (aka, the pressure value off the hanging center
       node). If we find it at more than 1 node, we have a serious problem
       so we're leaving. If we don't find any values, we set the element
       value to 0. RRl */
    found_quantity = FALSE;
    /* Only do this for elements with a haning center node, otherwise the
       extraction of the quantity can be ambiguous for non-regular grid models */
    if (ielem_type == BIQUAD_QUAD || ielem_type == TRIQUAD_HEX || ielem_type == C_BILINEAR_QUAD ||
        ielem_type == C_TRILINEAR_HEX) {
      for (i = 0; i < num_local_nodes; i++) {
        I = Proc_Elem_Connect[iconnect_ptr + i];
        /* NOTE: here, the element variables (such as PRESSURE) are being
           extracted from the solution vector coming off of the hanging
           interior nodes, or a given specified node for such a quantity.
           There should never be more than one of this quantity defined
           per element, or we have a problem treating it as an element
           variable. Hence the found_quantity check.                       */
        index = Index_Solution(I, varType, k, idof, matID, pg->imtrx);
        int base_index = exo->eb_ghost_elem_to_base[eb_index][ielem - e_start];
        if (index != -1 && base_index != -1) {
          /* This should be the one node that has our value - set the element
             value to this */
          sol_vec[index] = elem_vec[base_index];
          if (found_quantity == TRUE) {
            fprintf(stderr,
                    "Warning: Too many nodes returning quantities for element variable %s (%s) - "
                    "may not be accurate\n",
                    Exo_Var_Names[varType].name2, Exo_Var_Names[varType].name1);
            exit(-1);
          }
          found_quantity = TRUE;
        }
      }
    } else {
      int i = 0;
      I = Proc_Elem_Connect[iconnect_ptr + i];
      /* NOTE: here, the element variables (such as PRESSURE) are being
         extracted from the solution vector coming off of the hanging
         interior nodes, or a given specified node for such a quantity.
         There should never be more than one of this quantity defined
         per element, or we have a problem treating it as an element
         variable. Hence the found_quantity check.                       */
      index = Index_Solution(I, varType, k, idof, matID, pg->imtrx);
      int base_index = exo->eb_ghost_elem_to_base[eb_index][ielem - e_start];
      if (index != -1 && base_index != -1) {
        /* This should be the one node that has our value - set the element
           value to this */
        sol_vec[index] = elem_vec[ielem - e_start];
        if (found_quantity == TRUE) {
          fprintf(stderr,
                  "Warning: Too many nodes returning quantities for element variable %s (%s) - may "
                  "not be accurate\n",
                  Exo_Var_Names[varType].name2, Exo_Var_Names[varType].name1);
          exit(-1);
        }
        found_quantity = TRUE;
      }
    }
  }

  return;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/
/* build_node_index_var() - given the global node index, node_index_global,
 *                          then find the nodes having dof corresponding
 *                          a var named variable and populate the vectors
 *                          map[count] with the global node number, and
 *                          var_node_list[count] with the local node number
 *
 * Notes: the number of values destined to reside in map and hence
 *        the size of map must already be set, perhaps using countvar_dof()
 *        before allocating the map array
 *
 * Created: 1999/12/22 srsubia@sandia.gov
 *
 * Revised:
 */

int build_node_index_var(const int varType,
                         const int num_nodes,
                         const int *node_index_global,
                         int *map,
                         int *var_node_list) {
  int node, count;
  int imtrx;
  NODAL_VARS_STRUCT *nv;
  if (varType < 0 || varType > MAX_VARIABLE_TYPES - 1) {
    GOMA_EH(-1, "Attempt to count a bogus variable.");
  }
  count = 0;
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    for (node = 0; node < num_nodes; node++) {
      nv = Nodes[node]->Nodal_Vars_Info[imtrx];
      if (nv->Num_Var_Desc_Per_Type[varType] > 0) {
        map[count] = node_index_global[node];
        var_node_list[count] = node;
        count++;
      }
    }
  }
  return count;
}
/***************************************************************************/

/* count_vardofs() -- given the variable type, varType, this function
 *                    returns the number of degrees of freedom in the
 *                    solution vector corresponding to that variable type
 *                    at the nodes from zero to num_nodes.
 *                    Note: Dofs for MASS_FRACTION unknowns
 *                          with multiple subvariables are counted as
 *                          one degree of freedom here, just as in the
 *                          original Dolphin array.
 *
 * Created: 1999/12/22 srsubia@sandia.gov
 *
 * Revised: HKM
 */

int count_vardofs(const int varType, const int num_nodes) {
  int node, nun, count = 0;
  int imtrx;
  NODAL_VARS_STRUCT *nv;
  if (varType < 0 || varType > MAX_VARIABLE_TYPES - 1) {
    GOMA_EH(-1, "Attempt to count a bogus variable.");
  }
  count = 0;
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    for (node = 0; node < num_nodes; node++) {
      nv = Nodes[node]->Nodal_Vars_Info[imtrx];
      nun = get_nv_ndofs_modMF(nv, varType);
      count += nun;
    }
  }
  return (count);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void fprint_strn(FILE *fd, const char *string, const int num)

/************************************************************************
 *
 * print_strn():
 *
 *  Prints a string a repetitive number of times to a named stream.
 ************************************************************************/
{
  int i;
  for (i = 0; i < num; i++)
    fprintf(fd, "%s", string);
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void fprint_line(FILE *fd, const char *string, const int num)

/************************************************************************
 *
 * fprint_line():
 *
 *  Prints a string a repetitive number of times to a named stream.
 *  Then, prints the the end of line character.
 ************************************************************************/
{
  int i;
  for (i = 0; i < num; i++)
    fprintf(fd, "%s", string);
  fprintf(fd, "\n");
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void init_vec_value(double *vector, const double value, const int length)

/************************************************************************
 *
 * init_vec_value():
 *
 *  Initializes a vector of doubles of length, length, to a constant
 *  value, value.
 ************************************************************************/
{
  int i;

  if (vector == NULL) {
    if (Debug_Flag > 1) {
      GOMA_WH(-1, "Warning: attempted to initialize NULL vector.");
    }
    return;
  }

  if (value == 0.0) {
    (void)memset((void *)vector, 0, (sizeof(double) * length));
  } else {
    for (i = 0; i < length; i++)
      *vector++ = value;
  }
  return;
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void dcopy1(const int length, const double *src_vector, double *dst_vector)

/************************************************************************
 *
 * dcopy1():
 *
 *  Copies the contents of one double vector to another. A stride of 1
 *  is assumed.  Note, this
 *  version doesn't do any unrolling and doesn't attempt to use
 *  memcopy. Both of these options may lead to speed increases.
 ************************************************************************/
{
  int i;
  for (i = 0; i < length; i++) {
    *dst_vector++ = *src_vector++;
  }
  return;
}

int find_first_elem_with_var(Exo_DB *e, int var) {
  int eb, mn, found = FALSE, first_elem = -1;

  for (eb = 0; (!found) && eb < e->num_elem_blocks; eb++) {
    mn = Matilda[eb];

    found = pd_glob[mn]->v[pg->imtrx][var];
    if (found)
      first_elem = e->eb_ptr[eb];
  }

  return (first_elem);
}

#ifdef LIBRARY_MODE
int load_import_fields(dbl *base_p_por, const Exo_DB *exo, int callnum) {
  int i, w;
  int nv_start, ev_start, nv_count = 0, ev_count = 0;
  int ne = exo->num_elems;
  int nn = exo->num_nodes;
  int doing_ev = FALSE;
  dbl *ev_tmp = NULL;

  /* Make sure there are fields to load */
  if (Num_Import_NV == 0 && Num_Import_EV == 0)
    return 0;
  if (libio->xnv_in == NULL || libio->xev_in == NULL) {
    fprintf(stderr, "Invalid import vector(s)!\n");
    return -1;
  }
  if (Num_Import_EV > 0)
    ev_tmp = alloc_dbl_1(nn, 0.0);

  /* Loop over all defined external fields */
  for (w = 0; w < efv->Num_external_field; w++) {

    /* Import the nodal vars */
    if (!strcmp(efv->file_nm[w], "IMPORT")) {
      fprintf(stderr, " Importing external #%d nvar field %s...\n", w, efv->name[w]);
      efv->ext_fld_ndl_val[w] = alloc_dbl_1(nn, 0.0);
      nv_start = nn * nv_count;
      if (doing_ev)
        GOMA_WH(-1, "Problem with order of external field cards!\n");
      for (i = 0; i < nn; i++) {
        efv->ext_fld_ndl_val[w][i] = libio->xnv_in[nv_start + i];

        if (base_p_por != NULL && w == efv->ev_porous_index) {
          base_p_por[i] = efv->ext_fld_ndl_val[w][i];
        }
      }
      nv_count++;
    }
  }

  /****************************Anneal from external***********************/
  if (efv->ev_porous_decouple) {
    anneal_mesh_with_external_field(exo);
  }

  for (w = 0; w < efv->Num_external_field; w++) {
    /*
     * Element vars loaded second because they have to projected to nodes.  For this
     *   we may need the external nodal displacements to already have been imported
     */
    if (!strcmp(efv->file_nm[w], "IMPORT_EV")) {
      fprintf(stderr, " Importing external #%d evar field %s...\n", w, efv->name[w]);
      efv->ext_fld_ndl_val[w] = alloc_dbl_1(nn, 0.0);
      doing_ev = TRUE;
      vzero(nn, &ev_tmp[0]);
      interp_ev_to_nodes(exo, ev_tmp, ev_count);
      for (i = 0; i < nn; i++) {
        efv->ext_fld_ndl_val[w][i] = ev_tmp[i];

        /* Save base porosity for updates if applicable */
        if (base_p_por != NULL && w == efv->ev_porous_index) {
          base_p_por[i] = efv->ext_fld_ndl_val[w][i];
        }
      }
      ev_count++;
    }

    /* All other cases: var has already been read in from Exodus file! */
  }

  fprintf(stderr, " Imported %d of %d nodal vars and %d of %d elem vars.\n", nv_count,
          Num_Import_NV, ev_count, Num_Import_EV);
  safer_free((void **)&ev_tmp);
  return 0;
}

void interp_ev_to_nodes(const Exo_DB *exo, dbl *ev_tmp, int iev) {
  int ne1, ne2, nel;
  int i, j, n, el;
  int nn = exo->num_nodes;
  int ne = exo->num_elems;
  dbl sum, avg;

  dbl hsquared[DIM] = {0., 0., 0.}, hh[DIM][DIM], dhh[DIM][MDE];

  dbl h_siz, h_siz_sum = 0.0;

  /* Set index into xev_in */
  n = iev * ne;

  /* Loop over all nodes */
  for (i = 0; i < nn; i++) {
    ne1 = exo->node_elem_pntr[i];
    ne2 = exo->node_elem_pntr[i + 1];
    nel = ne2 - ne1;
    sum = avg = 0.0;

    h_siz_sum = 0.0;

    /* Sum var value over elements containing this node */
    for (j = ne1; j < ne2; j++) {
      el = exo->node_elem_list[j];

      ei[pg->imtrx]->ielem = el;
      ei[pg->imtrx]->ielem_type = Elem_Type(exo, el);
      ei[pg->imtrx]->num_local_nodes = elem_info(NNODES, ei[pg->imtrx]->ielem_type);

      h_elem_siz(hsquared, hh, dhh, FALSE);

      h_siz = hsquared[0] * hsquared[1];
      if (pd->Num_Dim > 2)
        h_siz *= hsquared[2];

      h_siz = sqrt(h_siz);

      h_siz_sum += h_siz;

      sum += libio->xev_in[n + el] * h_siz;
    }
    avg = sum / h_siz_sum;

    /* Write nodal value */
    ev_tmp[i] = avg;
  }

  /* Done */
  return;
}

int advance_porosity_ev(const int time_step, const int nn, dbl *x, dbl *base_p_por, dbl *base_p_liq)

/*
 * Function updates the external field corresponding to imported nodal
 * porosity at each time step (after the first) according to the
 * following constitutive relation:
 *
 *         n+1        n             0             n+1          n
 *    (por)    = (por)  + Cr * (por)   * [ (p_liq)    - (p_liq)  ]
 * where:
 *        Cr = porous ("rock") compressibility of solid
 *       por = porosity
 *     p_liq = porous liquid pressure
 * and superscripts denote the time steps:
 *       n+1 = current
 *         n = first in current Goma call
 *         0 = start of problem (not current Goma call)
 *
 * This update is done here on a node-by-node basis, rather than
 * doing at assembly time.
 */
{
  static int previous_step = -1;
  double cr, por0, delta;
  int i, n, mn, i_por_ev;
  int eq = POR_LIQ_PRES;

  /* Bail out if update not needed */
  if (time_step == previous_step)
    return 0;
  previous_step = time_step;
  if (time_step == 0)
    return 0;

  /* check for valid porosity external field index */
  i_por_ev = efv->ev_porous_index;
  GOMA_EH(i_por_ev, "Porosity external field not found!");

  /* Loop over nodes in problem */
  for (n = 0; n < nn; n++) {

    /*
     * Get block number and DOF index for porous liquid pressure at node
     * EDW note: This will only work for nodes along an element block
     * boundary when the porous material has the smaller material index!
     */
    mn = first_matID_at_node(n);
    i = Index_Solution(n, eq, 0, 0, -1, pg->imtrx);

    /* Proceed if eq is active at node */
    if (i > -1) {
      cr = mp_glob[mn]->porous_compressibility;
      por0 = mp_glob[mn]->initial_porosity;
      delta = cr * por0 * (x[i] - base_p_liq[n]);
      efv->ext_fld_ndl_val[i_por_ev][n] = base_p_por[n] + delta;
    }
  }

  return 0;
}

int advance_porosity_goma_first(
    const int time_step, const int nn, dbl *x, dbl *base_p_por, dbl *base_p_liq)

/*
 * Function updates the external field corresponding to imported nodal
 * porosity at each time step (after the first) according to the
 * following constitutive relation:
 *
 *         n+1        n             0             n+1          n
 *    (por)    = (por)  + Cr * (por)   * [ (p_liq)    - (p_liq)  ]
 * where:
 *        Cr = porous ("rock") compressibility of solid
 *       por = porosity
 *     p_liq = porous liquid pressure
 * and superscripts denote the time steps:
 *       n+1 = current
 *         n = first in current Goma call
 *         0 = start of problem (not current Goma call)
 *
 * This update is done here on a node-by-node basis, rather than
 * doing at assembly time.
 */
{
  static int previous_step = -1;
  double cr, por0, delta;
  int i, n, mn, i_por_ev;
  int eq = POR_LIQ_PRES;

  /* Bail out if update not needed */
  if (time_step == previous_step)
    return 0;
  previous_step = time_step;
  if (time_step == 0)
    return 0;

  /* check for valid porosity external field index */
  i_por_ev = efv->ev_porous_index;
  GOMA_EH(i_por_ev, "Porosity external field not found!");

  /* Loop over nodes in problem */
  for (n = 0; n < nn; n++) {

    /*
     * Get block number and DOF index for porous liquid pressure at node
     * EDW note: This will only work for nodes along an element block
     * boundary when the porous material has the smaller material index!
     */
    mn = first_matID_at_node(n);
    i = Index_Solution(n, eq, 0, 0, -1, pg->imtrx);

    /* Proceed if eq is active at node */
    if (i > -1) {
      cr = mp_glob[mn]->porous_compressibility;
      por0 = mp_glob[mn]->initial_porosity;
      delta = cr * por0 * (x[i] - base_p_liq[n]);
      efv->ext_fld_ndl_val[i_por_ev][n] = base_p_por[n] + delta;
    }
  }

  return 0;
}

int advance_porosity_jas_leads(const int time_step,
                               double current_gtime,
                               dbl starting_gtime,
                               dbl animas_time_step,
                               const int nn,
                               dbl *p_por_final,
                               dbl *p_por_dot) {
  int i, n, mn, i_por_ev;
  int eq = POR_LIQ_PRES;

  /* check for valid porosity external field index */
  i_por_ev = efv->ev_porous_index;
  GOMA_EH(i_por_ev, "Porosity external field not found!");

  /* Loop over nodes in problem */
  for (n = 0; n < nn; n++) {
    efv->ext_fld_ndl_val[i_por_ev][n] =
        p_por_final[n] - p_por_dot[n] * (animas_time_step + current_gtime - starting_gtime);
  }

  return 0;
}

#endif /* LIBRARY_MODE */

int init_pmv_hyst(const Exo_DB *exo)

/**********************************************************************
 * Function to initiate porous media variables hysteresis structure
 * containing history of curve switching and other nodal quantities
 * This update is done here on a node-by-node basis, rather than
 * doing at assembly time.
 ***********************************************************************/

{
  int i_num_switch_ext_field, i_curve_type_ext_field;

  int inode, mn, ipore;

  /* Loop over nodes in problem */
  for (inode = 0; inode < exo->num_nodes; inode++) {

    /*
     * Get block number and DOF index for porous liquid pressure at node
     * EDW note: This will only work for nodes along an element block
     * boundary when the porous material has the smaller material index!
     */
    mn = first_matID_at_node(inode);

    for (ipore = 0; ipore < pd_glob[mn]->Num_Porous_Shell_Eqn; ipore++) {

      /* Initialize switch trigger to OFF */
      pmv_hyst->curve_switch[ipore][inode] = 0;

      /* Initialize number of curve switches from external field */
      i_num_switch_ext_field =
          mp_glob[mn]->por_shell_cap_pres_hyst_num_switch_ext_field_index[ipore];
      if (i_num_switch_ext_field > -1) {
        pmv_hyst->num_switch[ipore][inode] =
            (int)efv->ext_fld_ndl_val[i_num_switch_ext_field][inode];
      }

      /* Initialize curve types from external field */
      i_curve_type_ext_field =
          mp_glob[mn]->por_shell_cap_pres_hyst_curve_type_ext_field_index[ipore];
      if (i_curve_type_ext_field > -1) {
        pmv_hyst->curve_type[ipore][inode] =
            (int)efv->ext_fld_ndl_val[i_curve_type_ext_field][inode];
        pmv_hyst->curve_type_old[ipore][inode] =
            (int)efv->ext_fld_ndl_val[i_curve_type_ext_field][inode];
      }

      /*
       * Initialize imbibition minimum and drainage maximum saturation values from main curve
       * parameters specified in the material file
       */
      pmv_hyst->sat_min_imbibe[ipore][inode] = mp_glob[mn]->u_PorousShellCapPres[ipore][0];
      pmv_hyst->sat_max_drain[ipore][inode] = mp_glob[mn]->u_PorousShellCapPres[ipore][5];

      /* Finally initialize the saturation and capillary pressure values at switching points*/
      if (pmv_hyst->curve_type[ipore][inode] == 1) /* If starting at draining curve */
      {
        pmv_hyst->sat_switch[ipore][inode] =
            mp_glob[mn]->u_PorousShellCapPres[0][4]; /* Switch at the minimum saturation */
        pmv_hyst->cap_pres_switch[ipore][inode] =
            1.0e-12; /* Set to arbitrary small capillary pressure */
      } else if (pmv_hyst->curve_type[ipore][inode] == 0) /* If starting at imbibition curve */
      {
        pmv_hyst->sat_switch[ipore][inode] =
            mp_glob[mn]->u_PorousShellCapPres[0][1]; /* Switch at the maximum saturation */
        pmv_hyst->cap_pres_switch[ipore][inode] =
            1.0e12; /* Set to arbitrary large capillary pressure */
      }
    }
  }

  return 0;
}

/**********************************************************************
 * Function updates the external field corresponding to imported nodal
 * etch area  at each time step (after the first) according to the
 * following relation:
 *
 *              n+1              n                     n+1    n
 *    (etch_area)    = (etch_area)  - etch_rate   * [ t    - t  ]
 * where:
 *        etch rate = Rate of etching in CGS unit
 *        t = time
 *       n+1 = current
 *         n = first in current Goma call
 *
 * This update is done here on a node-by-node basis, rather than
 * doing at assembly time.
 ***********************************************************************/
int advance_etch_area_ext_field(const int time_step, const int nn, const dbl delta_t, dbl *x) {
  int i_etch_area = -1, i_etch_depth = -1;
  int n, i_H2O, i_KOH;
  double rho_H2O, rho_KOH;
  double etch_area_old, etch_depth_old, etch_area, etch_depth;
  double etch_rate;

  i_etch_area = efv->ev_etch_area;
  i_etch_depth = efv->ev_etch_depth;

  /* Loop over nodes in problem */
  for (n = 0; n < nn; n++) {
    /*
     * Load up concentration field - g/cm^3
     */
    i_H2O = Index_Solution(n, MASS_FRACTION, 0, 0, -1, pg->imtrx);
    i_KOH = Index_Solution(n, MASS_FRACTION, 1, 0, -1, pg->imtrx);
    rho_H2O = x[i_H2O];
    rho_KOH = x[i_KOH];

    /*
     * Load up etch area and depth from external field
     */
    etch_area_old = efv->ext_fld_ndl_val[i_etch_area][n];
    etch_depth_old = efv->ext_fld_ndl_val[i_etch_depth][n];

    /*
     * Calculate etch rate
     */
    etch_rate = calc_KOH_Si_etch_rate_100(rho_H2O, rho_KOH, NULL);

    /*
     * Update etch area and depth
     */
    etch_area = etch_area_old - etch_rate * delta_t / (350.0 * 1.0e-7) / sqrt(2.0);
    etch_depth = etch_depth_old + etch_rate * delta_t;

    if (etch_area < 0.0) {
      etch_area = 0.0;
      etch_depth = etch_depth_old;
    }
    efv->ext_fld_ndl_val[i_etch_area][n] = etch_area;
    efv->ext_fld_ndl_val[i_etch_depth][n] = etch_depth;

  } /* End of loop over nodes */
  return 0;
} /* End of advance_etch_area_ext_field */

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int evaluate_sat_hyst_criterion_nodal(const dbl *x_static,
                                      const dbl *x_dot_static,
                                      const Exo_DB *exo)
/*
 * Function to determine when to switch curves on hysteretic models
 * of porous shell formulation
 *
 * This update is done here on a node-by-node basis, rather than
 * doing at assembly time.
 */
{
  int inode, mn, ipore, ie;
  int num_switch = 0;
  int curve_type = -1;
  int num_switch_max = 4; /* Hard set the maximum number of curve switch for now */
  double sat_node, sat_dot_node, cap_pres_node;
  int i_num_switch_ext_field, i_curve_type_ext_field;

  /* Loop over nodes in problem */
  for (inode = 0; inode < exo->num_nodes; inode++) {

    /*
     * Get block number and DOF index for porous liquid pressure at node
     * EDW note: This will only work for nodes along an element block
     * boundary when the porous material has the smaller material index!
     */
    mn = first_matID_at_node(inode);

    for (ipore = 0; ipore < pd_glob[mn]->Num_Porous_Shell_Eqn; ipore++) {

      num_switch = pmv_hyst->num_switch[ipore][inode];

      if (ipore == 0) {
        ie = Index_Solution(inode, R_SHELL_SAT_1, 0, 0, -1, pg->imtrx);
        sat_node = x_static[ie];
        sat_dot_node = x_dot_static[ie];
      } else if (ipore == 1) {
        ie = Index_Solution(inode, R_SHELL_SAT_2, 0, 0, -1, pg->imtrx);
        sat_node = x_static[ie];
        sat_dot_node = x_dot_static[ie];
      } else if (ipore == 2) {
        ie = Index_Solution(inode, R_SHELL_SAT_3, 0, 0, -1, pg->imtrx);
        sat_node = x_static[ie];
        sat_dot_node = x_dot_static[ie];
      }

      curve_type = pmv_hyst->curve_type[ipore][inode];

      /* Call the load_cap_pres with curve switch OFF */
      pmv_hyst->curve_switch[ipore][inode] = 0;
      cap_pres_node = load_cap_pres(ipore, -1, inode, sat_node);

      /*If the accumulation is positive, above a certain
       *tolerance level, and was positive previously, then we will
       *remain on the same wetting curve with same previous curve
       *parameters.   Likewise for negative accumulation with the
       *drying curve.    We will switch if the accumulation rate
       *changes sign and the magnitude is above a certan tolerance
       *level.
       */

      if ((sat_dot_node > 0.0) && (curve_type == 0)) {
        /* We were on imbibition curve, and will remain so */
        pmv_hyst->curve_switch[ipore][inode] = 0;
      }

      else if ((sat_dot_node < 0.0) && (curve_type == 1)) {
        /* We were on a draining curve, and will remain so */
        pmv_hyst->curve_switch[ipore][inode] = 0;
      }

      else if ((sat_dot_node > 0.0) && (curve_type == 1) && (num_switch < num_switch_max)) {
        /*
         * We were on a draining curve but now may
         * potentially switch to a imbibition curve
         */

        if (fabs(sat_dot_node) > mp_glob[mn]->u_PorousShellCapPres[0][9])
        /* If greater than tolerance, switch to imbibition */
        {
          pmv_hyst->curve_switch[ipore][inode] = 1;
          pmv_hyst->num_switch[ipore][inode] += 1;

          pmv_hyst->curve_type_old[ipore][inode] = 1;
          pmv_hyst->curve_type[ipore][inode] = 0;

          pmv_hyst->sat_switch[ipore][inode] = sat_node;
          pmv_hyst->cap_pres_switch[ipore][inode] = cap_pres_node;
        } else
        /* Else, don't switch */
        {
          pmv_hyst->curve_switch[ipore][inode] = 0;
        }
      }

      else if ((sat_dot_node < 0.0) && (curve_type == 0) && (num_switch < num_switch_max)) {
        /*
         * We were on a imbibition curve but now may
         * potentially switch to a draining curve
         */

        if (fabs(sat_dot_node) > mp_glob[mn]->u_PorousShellCapPres[0][9])
        /* If greater than tolerance, switch to draining */
        {
          pmv_hyst->curve_switch[ipore][inode] = 1;
          pmv_hyst->num_switch[ipore][inode] += 1;

          pmv_hyst->curve_type_old[ipore][inode] = 0;
          pmv_hyst->curve_type[ipore][inode] = 1;

          pmv_hyst->sat_switch[ipore][inode] = sat_node;
          pmv_hyst->cap_pres_switch[ipore][inode] = cap_pres_node;
        } else
        /* Else, don't switch */
        {
          pmv_hyst->curve_switch[ipore][inode] = 0;
        }
      }

      /* Finally, update the external field variables */

      /* Update number of curve switches */
      i_num_switch_ext_field =
          mp_glob[mn]->por_shell_cap_pres_hyst_num_switch_ext_field_index[ipore];
      if (i_num_switch_ext_field > -1) {
        efv->ext_fld_ndl_val[i_num_switch_ext_field][inode] =
            (double)pmv_hyst->num_switch[ipore][inode];
      }

      /* Update curve type */
      i_curve_type_ext_field =
          mp_glob[mn]->por_shell_cap_pres_hyst_curve_type_ext_field_index[ipore];
      if (i_curve_type_ext_field > -1) {
        efv->ext_fld_ndl_val[i_curve_type_ext_field][inode] =
            (double)pmv_hyst->curve_type[ipore][inode];
      }
    }
  }
  return (1);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

/* END of file rf_util.c  */
