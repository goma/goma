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

/* mm_fill_ls -- auxiliary routines devoted to initialization and tracking of
 *               level set interface representations
 */

/*
 *$Id: mm_fill_ls.c,v 5.21 2009-11-13 23:20:08 prschun Exp $
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "std.h" /* This needs to be here. */
#include "ac_stability.h"
#include "ac_stability_util.h"
#include "bc_contact.h"
#include "dp_comm.h"
#include "dp_types.h"
#include "dp_utils.h"
#include "dpi.h"
#include "el_elm_info.h"
#include "exo_conn.h"
#include "exo_struct.h"
#include "exodusII.h"
#include "mm_as_alloc.h"
#include "mm_fill_aux.h"
#include "mm_fill_fill.h"
#include "mm_fill_ls.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_shell.h"
#include "mm_fill_terms.h"
#include "mm_fill_util.h"
#include "mm_flux.h"
#include "mm_post_def.h"
#include "mm_qtensor_model.h"
#include "mm_unknown_map.h"
#include "mpi.h"
#include "rd_mesh.h"
#include "rf_node_const.h"
#include "rf_shape.h"
#include "sl_auxutil.h"
#include "sl_util_structs.h"
#include "user_pre.h"

#ifdef PARALLEL
#ifndef MPI
#define MPI /* otherwise az_aztec.h trounces MPI_Request */
#endif
#endif

#include "az_aztec.h"

/* goma include files (of course!) */

#include "el_elm.h"
#include "el_geom.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "rf_allo.h"
#include "rf_bc.h"
#include "rf_bc_const.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_io.h"
#include "rf_io_const.h"
#include "rf_mp.h"
#include "rf_solver.h"
#include "rf_vars_const.h"
#include "mm_eh.h"
#include "mm_mp.h"
#include "mm_mp_structs.h"
#include "sl_util.h"

#define GOMA_MM_FILL_LS_C
#include "sl_epetra_util.h"

/*
static int interface_on_side
( double isoval,
        int ielem_type,
        int id_side ));
*/
static void retrieve_node_coordinates(int, double *, double *, double **);

static struct LS_Surf_List *create_surfs_from_ns(int, double *, Exo_DB *);

static struct LS_Surf_List *create_surfs_from_ss(int, double *, Exo_DB *);

static struct LS_Surf_List *create_surfs_from_iso(int, double, double *,
                                                  Exo_DB *);

static void initialize_sign(int, double *, Exo_DB *);

static double initial_level_set(double, double, double);

static double gradient_norm_err(dbl *, Exo_DB *, Dpi *, dbl);

static int Hrenorm_simplemass(Exo_DB *exo, Comm_Ex *cx, Dpi *dpi, double x[],
                              struct LS_Surf_List *list, int num_total_nodes,
                              int num_ls_unkns, int num_total_unkns,
                              double time);

static int Hrenorm_smolianksi_only(Exo_DB *exo, Comm_Ex *cx, Dpi *dpi,
                                   double x[], struct LS_Surf_List *list,
                                   int num_total_nodes, int num_ls_unkns,
                                   int num_total_unkns, double time);

static int Hrenorm_constrain(Exo_DB *, Comm_Ex *, Dpi *, double[],
                             struct LS_Surf_List *, int, int, int, double);

static double find_LS_mass(const Exo_DB *, const Dpi *, const double *,
                           double *, double[], int);

static void find_quad_facets(struct LS_Surf_List *, int, double, int *,
                             Exo_DB *);

static int stash_node_displacements(double **, int, double *, Exo_DB *);

static int unique_surf(struct LS_Surf_List *, struct LS_Surf *);

static int point_on_ca_boundary(int, Exo_DB *);

static void find_intersections(struct LS_Surf_List *, int, double, Exo_DB *);

static int knockout_point(double[DIM], double[DIM], double);

static struct LS_Surf_List *create_surfs_from_arc(double[DIM], double,
                                                  double[DIM], double);

static double find_arc_region_sign(struct LS_Surf_Arc_Data *, double *);

static void compute_link_level_set(double *, double *, double, int, double,
                                   double *, double *, Integ_Elem *);

static struct LS_Surf *next_surf_or_subsurf(struct LS_Surf_List *, int *);

static struct LS_Surf *create_next_surf_or_subsurf(struct LS_Surf_List *, int,
                                                   int);

static double distance(double *, double *);

#ifdef PARALLEL
static void ddd_add_surf(DDD, struct LS_Surf *);
#endif

static int interface_in_grid(SGRID *);

static void divide_shape_fcn_tree(NTREE *, int);

static void compute_shape_fcn_values(NTREE *);

static void gather_subtree_coords(NTREE *, double *, double (*)[DIM]);

static void load_subtree_coords(int, NTREE *, double (*)[DIM]);

/*
static void print_shape_fcn_tree
( NTREE *)) ;
*/

static void find_grid_LS_value(SGRID *);

static int current_elem_xfem_state(int[], int *, double[], const Exo_DB *);

static double scalar_value_at_local_node(int, int, int, int, int, double *,
                                         Exo_DB *);

/*
static double element_average
( int,
                double [],
                const Exo_DB* ,
                int );
*/

static int vertex_on_element_boundary(double[DIM], double *, int *);

static void copy_distance_function(double *, double **);

static double determine_adc_probability(struct Boundary_Condition *, int,
                                        Exo_DB *, double *, int, int, int *,
                                        double, double *);

static void apply_adc_to_ss(Exo_DB *, double *, int, double);

static void apply_adc_to_elem(Exo_DB *, double *, int, int, int, int *, double,
                              double);

static double determine_nearest_distance(Exo_DB *, double *, int, int *);

static int check_alignment(double);

static struct LS_Surf *closest_other_surf(struct LS_Surf_List *, double *,
                                          Exo_DB *, struct LS_Surf *);

static double find_adc_node(int, double *, Exo_DB *, struct LS_Surf *, double *,
                            int *);

static void apply_adc_function(double *, Exo_DB *, double *, double, double);

static int is_LS_spurious(Exo_DB *, double *, int, double, double *);
static void purge_spurious_LS(double *, Exo_DB *, int);

#define EXPLICIT FALSE
#define MAX_STEP 500
#define SUBELEM_SIG_CROSS_TOL 1.e-6

#ifndef COUPLED_FILL
#define GRADIENT TRUE
#define STREAMWISE FALSE
void semi_lagrange_step(const int num_total_nodes, int num_total_unknowns,
                        int num_fill_unknowns, double x[], double F[],
                        double F_old[], double Fdot[], double Fdot_old[],
                        int node_to_fill[], double delta_t, double theta,
                        Exo_DB *exo, Dpi *dpi, Comm_Ex *cx)

{

  int inode, i;
  int eqn = LS;
  int dim = pd->Num_Dim;

  int global_fill_unknowns;

  int *ie_to_fill, *ext_dof;
  double *F_, *R, *b, *dC;

  double M0 = 0.0, global_LS_flux, M;

  double R_lamda, lamda, d_lamda, R_norm = 1.0, delta_norm = 1.0;

  int max_its;

#ifdef PARALLEL
  int local_fill_unknowns = num_fill_unknowns;

  double local_bTb, local_bTR, local_RTR;
  double global_bTb = 0.0, global_bTR = 0.0, global_RTR = 0.0;

  exchange_dof(cx, dpi, x, pg->imtrx);

#endif

  if (dim > 2)
    EH(-1, "SEMI_LAGRANGE time stepping not implemented in 3D.");
  if (ls->Isosurface_Subsurf_Type != LS_SURF_FACET)
    EH(-1, "Need FACET Isosurface type .");

  global_fill_unknowns = num_fill_unknowns;

  ext_dof = alloc_int_1(num_fill_unknowns, FALSE);

  ie_to_fill = alloc_int_1(num_fill_unknowns, 0);

  dalloc(num_fill_unknowns, F_);
  dalloc(num_fill_unknowns, R);
  dalloc(num_fill_unknowns, b);
  dalloc(num_total_unknowns, dC);

  for (inode = 0; inode < num_total_nodes; inode++) {

    double r_node[DIM];
    double v_node[DIM];
    double r_last[DIM];

    struct LS_Surf *closest;

    double distance;

    if (num_varType_at_node(inode, eqn) == 1) {
      r_node[0] = Coor[0][inode];
      r_node[1] = Coor[1][inode];

      if ((num_varType_at_node(inode, VELOCITY1) == 1) &&
          (num_varType_at_node(inode, VELOCITY2) == 1)) {
        v_node[0] = x[Index_Solution(inode, VELOCITY1, 0, 0, -1, pg->imtrx)];
        v_node[1] = x[Index_Solution(inode, VELOCITY2, 0, 0, -1, pg->imtrx)];
      } else {
        EH(-1, "Need equal interpolation order LS and VELOCITY1, VELOCITY2 for "
               "SEMI_LAGRANGE.\n");
      }

      r_last[0] = r_node[0] - v_node[0] * delta_t;
      r_last[1] = r_node[1] - v_node[1] * delta_t;

      closest = closest_surf(ls->last_surf_list, x, exo, r_last);

      distance = closest->closest_point->distance;

      F_[node_to_fill[inode]] = distance;

      ie_to_fill[node_to_fill[inode]] =
          Index_Solution(inode, ls->var, 0, 0, -2, pg->imtrx);

#ifdef PARALLEL
      /* Here we set up the array ext_dof which tells me whether this dof belong
       * to me or no */
      if (Num_Proc > 1) {
        if (inode > (dpi->num_internal_nodes + dpi->num_boundary_nodes - 1)) {
          ext_dof[node_to_fill[inode]] = TRUE;
          local_fill_unknowns -= 1;
        }
      }
#endif
    }
  }
#ifdef PARALLEL
  if (Num_Proc > 1) {

    global_fill_unknowns = 0;

    MPI_Allreduce(&local_fill_unknowns, &global_fill_unknowns, 1, MPI_INT,
                  MPI_SUM, MPI_COMM_WORLD);
  }

#endif

  M0 = find_LS_mass(exo, dpi, NULL, dC, x, num_total_unknowns);

  global_LS_flux =
      find_LS_global_flux(exo, dpi, NULL, NULL, x, num_total_unknowns);

  M0 -= global_LS_flux * delta_t;

  R_lamda = lamda = 0.0;

  max_its = 10;

  DPRINTF(stderr, "\n\t   L_2 in R   L_2 in dx   |dM| \n");
  DPRINTF(stderr, "\t  ---------  ---------   --------\n");

  while (delta_norm >= 1.e-5 && max_its > 0) {
    int k;
    double bTb, bTR, RTR;

    for (k = 0, bTb = bTR = RTR = 0.0; k < num_fill_unknowns; k++) {
      b[k] = dC[ie_to_fill[k]];
      R[k] = 2.0 * (F[k] - F_[k]) + lamda * b[k];

      if (ext_dof[k] == FALSE) {
        bTb += b[k] * b[k];
        bTR += b[k] * R[k];
        RTR += R[k] * R[k];
      }
    }

#ifdef PARALLEL

    if (Num_Proc > 1) {

      local_bTb = bTb;
      local_bTR = bTR;
      local_RTR = RTR;

      MPI_Allreduce(&local_bTb, &global_bTb, 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD);
      bTb = global_bTb;

      MPI_Allreduce(&local_bTR, &global_bTR, 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD);
      bTR = global_bTR;

      MPI_Allreduce(&local_RTR, &global_RTR, 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD);

      RTR = global_RTR;
    }

#endif

    R_norm = sqrt(fabs(RTR)) / global_fill_unknowns + fabs(R_lamda);

    DPRINTF(stderr, "\t%10.2e ", R_norm);

    d_lamda = (2.0 * R_lamda - bTR) / bTb;

    for (k = 0, delta_norm = 0.0; k < num_fill_unknowns; k++) {
      F[k] += -0.5 * R[k] - 0.5 * b[k] * d_lamda;

      x[ie_to_fill[k]] = F[k];

      if (ext_dof[k] == FALSE)
        delta_norm += pow(-0.5 * R[k] - 0.5 * b[k] * d_lamda, 2.0);
    }

    lamda += d_lamda;

#ifdef PARALLEL
    if (Num_Proc > 1) {
      double global_delta_norm = 0.0, local_delta_norm;

      local_delta_norm = delta_norm;

      MPI_Allreduce(&local_delta_norm, &global_delta_norm, 1, MPI_DOUBLE,
                    MPI_SUM, MPI_COMM_WORLD);
      delta_norm = global_delta_norm;
    }

    exchange_dof(cx, dpi, x, pg->imtrx);
#endif

    delta_norm =
        sqrt(delta_norm) / global_fill_unknowns + sqrt(d_lamda * d_lamda);

    M = find_LS_mass(exo, dpi, NULL, dC, x, num_total_unknowns);

    R_lamda = M - M0;

    max_its--;

    DPRINTF(stderr, "%10.2e   %10.2e \n", delta_norm, fabs(M - M0));
  }

  if (max_its == 0) {
    for (i = 0; i < num_fill_unknowns; i++)
      F[i] = F_[i];
  }

  for (i = 0; i < num_fill_unknowns; i++) {
    Fdot[i] = (1.0 + 2.0 * theta) / delta_t * (F[i] - F_old[i]) -
              2.0 * theta * Fdot_old[i];
  }
}

#endif /* COUPLED_FILL */

/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/

int apply_ls_inlet_bc(double afill[], int ijaf[], double x[], double rf[],
                      int node_to_fill[],
                      struct elem_side_bc_struct *elem_side_bc, Exo_DB *exo) {

  int i;
  int ie;
  int I;

  double r[DIM];
  double distance;

  struct LS_Surf *closest;

  for (i = 0; i < elem_side_bc->num_nodes_on_side; i++) {
    I = elem_side_bc->local_node_id[i];

    retrieve_node_coordinates(I, x, r, NULL);

    closest = closest_surf(ls->last_surf_list, x, exo, r);

    distance = closest->closest_point->distance;

    if (num_varType_at_node(I, ls->var) == 1) {
      ie = Index_Solution(I, ls->var, 0, 0, -1, pg->imtrx);

      if (closest->closest_point->confidence == 0) {
        if (x[ie] * distance < 0.0)
          distance *= -1.0;
      }
      rf[node_to_fill[I]] += BIG_PENALTY * (x[ie] - distance);

      afill[node_to_fill[I]] += BIG_PENALTY;
    }
  }

  return (0);
}

int apply_strong_fill_ca_bc(
    double afill[],       /* Jacobian matrix for fill equation  */
    int ijaf[],           /* pointer to nonzeros in Jacobian matrix   */
    double x[],           /* Solution vector for the current processor    */
    double rf[],          /* Residual vector for the current processor    */
    const double delta_t, /* current time step size                       */
    const double theta_,  /* parameter (0 to 1) to vary time integration
                            ( implicit - 0 to explicit - 1)  */
    int node_to_fill[], const int ielem, /* element number */
    int ielem_type,                      /* element type  */
    const int num_local_nodes, const int ielem_dim, const int iconnect_ptr,
    struct elem_side_bc_struct *elem_side_bc,   /* Pointer to an element side
                                                   boundary condition structure */
    const int num_total_nodes, const double ca, /* contact angle */
    const Exo_DB *exo)

/******************************************************************************
  Function which applies contact angle for fill equation (level set)
  This is applied as a stong integrated condition ON THE FILL EQN
  specified through STRONG_FILL_CA_BC

  Author:         D. R. Noble
  Date:           17 February 2000
  Revised:

******************************************************************************/
{
  /* LOCAL VARIABLES */
  int a, ip, i, j, I, J, dim; /* counters */
  int idof, jdof, ie, je, ja; /* more counters */
  int nodes_per_side;
  int local_elem_node_id[MAX_NODES_PER_SIDE];

  int eqn;
  int err; /* status variable for functions */
  int ip_total;
  int nvdofi, nvdofj, ki, kj;
  dbl grad_F[DIM]; /* gradient of Fill. */

  double phi_i;
  dbl rhs;
  double s, t, u; /* Gaussian quadrature point locations  */
  double xi[DIM]; /* Local element coordinates of Gauss point. */
  double wt; /* Quadrature weights units - ergs/(sec*cm*K) = g*cm/(sec^3*K) */
  double F;
  double alpha;
  double delta_func;
  /***************************************************************************/

  /* Find out the number of surface quadrature points
     -this is assumed independent of the surface */
  ip_total = elem_info(NQUAD_SURF, ielem_type);

  eqn = FILL;
  dim = pd->Num_Dim;

  /* initialize grad_F */
  for (a = 0; a < DIM; a++) {
    grad_F[a] = 0;
  }

  /* Surface integration over element */
  for (ip = 0; ip < ip_total; ip++) {
    /* find the quadrature point locations for current ip */
    find_surf_st(ip, ielem_type, elem_side_bc->id_side, pd->Num_Dim, xi, &s, &t,
                 &u);

    /* find the quadrature weight for current ip */
    wt = Gq_surf_weight(ip, ielem_type);

    /* ****************************************/
    err = load_basis_functions(xi, bfd);
    EH(err, "problem from load_basis_functions");

    err = beer_belly();
    EH(err, "beer_belly");

    /* precalculate variables at  current integration pt.*/
    err = load_fv();
    EH(err, "load_fv");

    err = load_bf_grad();
    EH(err, "load_bf_grad");

    err = load_fv_grads();
    EH(err, "load_fv_grads");

    /* calculate the determinant of the surface jacobian and the normal to
     * the surface all at one time */

    err = get_side_info(ielem_type, elem_side_bc->id_side, &nodes_per_side,
                        local_elem_node_id);
    EH(err, "get_side_info");

    surface_determinant_and_normal(
        ielem, ei[pg->imtrx]->iconnect_ptr, num_local_nodes,
        ei[pg->imtrx]->ielem_dim - 1, elem_side_bc->id_side, nodes_per_side,
        local_elem_node_id);

    do_LSA_mods(LSA_SURFACE);

    /* calculate fill gradient */
    for (a = 0; a < dim; a++) {
      if (!EXPLICIT) {
        grad_F[a] = fv->grad_F[a];
      } else {
        grad_F[a] = fv_old->grad_F[a];
      }
    }
    /* calculate continuous delta function that is only
     * non-zero near the zero level set.
     */
    alpha = 0.5 * ls->Length_Scale;
    F = fv->F;
    /* delta_func = exp( -fabs(F)/alpha )/2.0/alpha; */
    if (fabs(F) < alpha) {
      delta_func = 0.5 * (1. + cos(M_PIE * F / alpha)) / alpha;
    } else {
      delta_func = 0.;
    }
    /*
     * Put local contributions into global right-hand side
     * if it is not a right-hand side variable-it won't get added in
     * (contribution is zero)
     */
    if (af->Assemble_Residual) {
      rhs = 0.;
      for (i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
        I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];

        /* check for multiple dofs */
        nvdofi = Dolphin[pg->imtrx][I][eqn];

        for (ki = 0; ki < nvdofi; ki++) {
          rhs = cos(ca * M_PIE / 180.);

          for (a = 0; a < dim; a++) {
            rhs += grad_F[a] * fv->snormal[a];
          }

          idof = ei[pg->imtrx]->ln_to_first_dof[eqn][i] + ki;
          /* also convert from node number to dof number */
          phi_i = bf[eqn]->phi[idof];

          rf[node_to_fill[I] + ki] +=
              BIG_PENALTY * wt * fv->sdet * phi_i * delta_func * rhs;
        }
      }
    }

    if (af->Assemble_Jacobian) {
      for (i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
        I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
        nvdofi = Dolphin[pg->imtrx][I][eqn];

        for (ki = 0; ki < nvdofi; ki++) {
          ie = node_to_fill[I] + ki;
          idof = ei[pg->imtrx]->ln_to_first_dof[eqn][i] + ki;
          phi_i = bf[eqn]->phi[idof];

          /* derivatives of fill equation wrt to fill variable */
          for (j = 0; j < ei[pg->imtrx]->num_local_nodes; j++) {
            J = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + j];
            nvdofj = Dolphin[pg->imtrx][J][eqn];
            for (kj = 0; kj < nvdofj; kj++) {
              jdof = ei[pg->imtrx]->ln_to_first_dof[eqn][j] + kj;

              je = node_to_fill[J] + kj;
              ja = (ie == je) ? ie : in_list(je, ijaf[ie], ijaf[ie + 1], ijaf);
              EH(ja, "Could not find vbl in sparse matrix.");

              for (a = 0; a < dim; a++) {
                afill[ja] += BIG_PENALTY * wt * fv->sdet * phi_i * delta_func *
                             bf[eqn]->grad_phi[jdof][a] * fv->snormal[a];
              }
            }
          }
        }
      }
    }
  }

  return 0;
}
/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/

/*****************************************************************************/
/******************************************************************************/
/******************************************************************************/

/******************************************************************************
 * huygens_renormalization() : Renormalize the level set distance function.
 *
 * Input
 * =====
 * x                 = The vector of unknowns.
 * num_total_nodes   = Total number of nodes.
 * exo               = Exodus file pointer.
 * cx                =
 * dpi               =
 * num_ls_unknowns   = Total number of FILL unknowns.
 * num_total_unkns   = Total number of unknowns.
 * time              = The, uh, time.
 *
 * Output
 * ======
 * x                 = Modified vector of unknowns (FILL modified)
 *
 * Return
 * ======
 * 0                 = No renormalization done (x not modified).
 * 1                 = Renormalized (x modified).
 *
 ******************************************************************************/
int


huygens_renormalization ( double *x,
		const int num_total_nodes,
		Exo_DB *exo,
		Comm_Ex *cx,
		Dpi    *dpi,
		const int num_ls_unkns,
		const int num_total_unkns,
		const double time,
		const int Renorm_Now )
{
  struct LS_Surf_List *list = NULL;
  struct LS_Surf *isosurf = NULL;
  double ls_err;
  double tolerance = ls->Renorm_Tolerance;
  struct LS_Surf_Iso_Data *s;
  int status = 1;
  double renorm_width = ls->Length_Scale * ls->Control_Width;

  if (ls->Length_Scale == 0.) {
    renorm_width = ls->Control_Width *
                   global_h_elem_siz(x, x_old_static, xdot_static, x, exo, dpi);
  }

  ls_err = gradient_norm_err(x, exo, dpi, renorm_width);

  if (Renorm_Now || ls_err > tolerance) {
    /* Let's make a note of why we're renormalizing. */
    if (ls_err > tolerance) {
      DPRINTF(stdout, "\n\t Gradient norm error exceeds tolerance: %g > %g",
              ls_err, tolerance);
    }
    if (ls->Renorm_Countdown == 0) {
      DPRINTF(
          stdout,
          "\n\t Maximum number of steps without renormalization reached: %d",
          ls->Renorm_Freq);
    }
    DPRINTF(stdout, "\n\t Huygens renormalization : ");

    /* this call cleanses the LS field of "droplets" that surround exactly one
     * node */

    purge_spurious_LS(x, exo, num_total_nodes);

    list = create_surf_list();
    isosurf = create_surf(LS_SURF_ISOSURFACE);
    s = (struct LS_Surf_Iso_Data *)isosurf->data;
    s->isovar = ls->var;
    if (ls->Initial_LS_Displacement != 0.) {
      s->isoval = ls->Initial_LS_Displacement;
      ls->Initial_LS_Displacement = 0.;
    } else {
      s->isoval = 0.;
    }

    append_surf(list, isosurf);

    if (ls->Renorm_Method == HUYGENS) {
      surf_based_initialization(x, NULL, NULL, exo, num_total_nodes, list, time,
                                0., 0.);
    } else if (ls->Renorm_Method == HUYGENS_C) {
      Hrenorm_constrain(exo, cx, dpi, x, list, num_total_nodes, num_ls_unkns,
                        num_total_unkns, time);
    } else if (ls->Renorm_Method == HUYGENS_MASS_ITER) {
      Hrenorm_simplemass(exo, cx, dpi, x, list, num_total_nodes, num_ls_unkns,
                         num_total_unkns, time);
    } else if (ls->Renorm_Method == SMOLIANSKI_ONLY) {
      Hrenorm_smolianksi_only(exo, cx, dpi, x, list, num_total_nodes,
                              num_ls_unkns, num_total_unkns, time);
    } else {
      EH(-1, "You shouldn't actually be here. \n");
    }

    free_surf_list(&list);

    ls->Renorm_Countdown = ls->Renorm_Freq;
    status = 1;

    /*
     * the following flag's function is turn off the potential for saturation
     * hysteresis curve switching for 4 time steps after a renormalization step.
     * This hopefully avoids spurious curve switching related to renormalization
     * induced weirdness
     */

    ls->Sat_Hyst_Renorm_Lockout = 4;

    DPRINTF(stdout, "    done. \n");

  } else if (ls->Renorm_Freq == 0) {
    status = 0;
    DPRINTF(stdout, "\n\t Renormalization is disabled. \n");
  } else {
    status = 0;
    DPRINTF(stdout, "\n\t Renormalization unnecessary ( %g < %g ). \n", ls_err,
            tolerance);
  }

  return (status);
}
#ifndef COUPLED_FILL
/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/

void correct_level_set(struct Aztec_Linear_Solver_System *ams, double xf[],
                       double rf[], double x[], double x_old[],
                       double x_oldest[], int node_to_fill[],
                       int num_total_nodes, int num_fill_unknowns,
                       double init_step_size, double theta, int num_steps,
                       int eqntype, Exo_DB *exo, Dpi *dpi, Comm_Ex *cx) {
  int step = 0, I;
  double *xfdot, *xfdot_old, *xf_old = NULL;

  double time = 0.0;

  double h = ls->Length_Scale;

  double ave_level_grad_err;

  double width = ls->Control_Width;

  double tolerance = ls->Renorm_Tolerance;

  static double step_size = 0.0;

  short int diverge = FALSE;

  xf_old = alloc_dbl_1(num_fill_unknowns, 0.0);
  xfdot = alloc_dbl_1(num_fill_unknowns, 0.0);
  xfdot_old = alloc_dbl_1(num_fill_unknowns, 0.0);
  get_fill_vector(num_total_nodes, x, xf, node_to_fill);
  memset(rf, 0, sizeof(double) * num_fill_unknowns);

  if (step_size == 0.0)
    step_size = init_step_size;

  ave_level_grad_err = gradient_norm_err(x, exo, dpi, width * h);

  if (fabs(ave_level_grad_err) > tolerance) {

    DPRINTF(stderr, "\n\t Correcting \n");
    DPRINTF(stderr, "\n\t\t  step       L1           its  \n");
    DPRINTF(stderr, "\t\t--------  ---------    --------\n");

    while ((fabs(ave_level_grad_err) > tolerance) && step++ < num_steps) {
      int its;

      dcopy1(num_fill_unknowns, xf, xf_old);

      time += step_size;

      its = integrate_explicit_eqn(ams, rf, xf, xf_old, xfdot, xfdot_old, x,
                                   x_old, x_oldest, step_size, theta, &time,
                                   eqntype, node_to_fill, exo, dpi, cx);

      if (its > 0) {
        int i, I, K;
        double last_norm = ave_level_grad_err;

        put_fill_vector(num_total_nodes, x, xf, node_to_fill);
        exchange_dof(cx, dpi, x, pg->imtrx);
        put_fill_vector(num_total_nodes, x_old, xf, node_to_fill);
        exchange_dof(cx, dpi, x_old, pg->imtrx);

        ave_level_grad_err = gradient_norm_err(x, exo, dpi, width * h);

        diverge = ave_level_grad_err >= last_norm;

        DPRINTF(stderr, "\t\t   [%d]     %10.3e      %d       %10.4e\n", step,
                ave_level_grad_err, its, step_size);
      } else if (step_size > init_step_size / 1000.0) {
        time = 0.0;
        step_size /= 10.0;
        step = 0;
        get_fill_vector(num_total_nodes, x_oldest, xf, node_to_fill);
        put_fill_vector(num_total_nodes, x, xf, node_to_fill);
        exchange_dof(cx, dpi, x, pg->imtrx);
        put_fill_vector(num_total_nodes, x_old, xf, node_to_fill);
        exchange_dof(cx, dpi, x_old, pg->imtrx);
        memset(xf_old, 0, sizeof(double) * num_fill_unknowns);

        DPRINTF(stderr, "\t\t  Resetting correcting step size to %10.3e\n",
                step_size);
      } else {
        step = num_steps; /* force the while loop to exit */
        diverge = TRUE;
      }
    }

    if (step >= num_steps && diverge) {
      get_fill_vector(num_total_nodes, x_oldest, xf, node_to_fill);
      put_fill_vector(num_total_nodes, x, xf, node_to_fill);
      exchange_dof(cx, dpi, x, pg->imtrx);
      put_fill_vector(num_total_nodes, x_old, xf, node_to_fill);
      exchange_dof(cx, dpi, x_old, pg->imtrx);

      DPRINTF(stderr, "\n\t\t  Correcting step fails. No correction applied "
                      "this time step \n");

      if (step_size / 2.0 > init_step_size / 1000.0) {
        step_size /= 2.0;
      } else {
        step_size = init_step_size / 1000.0;
      }
    } else {
      if (step_size * 2.0 < init_step_size) {
        step_size *= 2.0;
      } else {
        step_size = init_step_size;
      }
    }
  } else {
    DPRINTF(stderr, "\n\t Correction of level set unnecessay: %f \n",
            ave_level_grad_err);
  }

  safer_free((void **)&xf_old);
  safer_free((void **)&xfdot);
  safer_free((void **)&xfdot_old);
}
#endif /* not COUPLED_FILL */

#if 1

/* assemble_level_correct -- Assemble Jacobian and rhs to level-set distance
 * function correction evolution equation.
 *
 *
 *
 */

#ifndef COUPLED_FILL

int assemble_level_correct(
    double afill[], /* Jacobian matrix for fill equation  */
    int ijaf[],     /* pointer to nonzeros in Jacobian matrix   */
    double rf[],    /* rhs vector   */
    double dt,      /* current time step size */
    double tt,      /* parameter to vary time integration from
                     * explicit (tt = 1) to implicit (tt = 0) */
    int node_to_fill[]) {
  /*
   * Some local variables for convenience...
   */

  int eqn;
  int dim;
  int a, b;

  int i, j;
  int status;
  int I, J;
  int idof, jdof;

  int ie, je, ja;
  int ki, kj, nvdofi, nvdofj;

  dbl F; /* Fill. */
  dbl F_0;
  dbl F_old; /* Fill at last time step. */
  dbl F_dot; /* Fill derivative wrt time. */

  dbl grad_F[DIM];      /* gradient of Fill. */
  dbl grad_F_mag = 0.0; /* magnitude of fill gradient vector */

  dbl grad_F0[DIM];
  dbl grad_F0_mag = 0.0;

  dbl vel_mag = 0.0; /* velocity magnitude */

  dbl w[DIM]; /* correction velocities */
  dbl d_w_dF[DIM][MDE];
  dbl S = 0.0; /* sign of distance function (-1 or +1 ) */
  dbl P1 = 0.0;

  dbl xx[DIM];    /* position field. */
  dbl x_old[DIM]; /* old position field. */
  dbl x_dot[DIM]; /* current position field derivative wrt time. */

  dbl num_diff;

  /*
   * Galerkin weighting functions for i-th energy residuals
   * and some of their derivatives...
   */

  dbl phi_i;

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;
  dbl h3; /* Volume element (scale factors). */
  dbl det_J;
  dbl wt;

  double rhs;

#if 0
  extern dbl vec_dot
  (const int n1,
	 dbl *v1,
	 dbl *v2);
#endif

  dbl alpha = 0.5 * ls->Length_Scale;
  dbl upwind = 0.5;

  status = 0;

  /*
   * Unpack variables from structures for local convenience...
   */

  dim = pd->Num_Dim;

  eqn = R_FILL;

  /*
   * Bail out fast if there's nothing to do...
   */

  if (!pd->e[pg->imtrx][eqn]) {
    return (status);
  }

  wt = fv->wt; /* Gauss point weight. */

  h3 = fv->h3; /* Differential volume element. */

  det_J = bf[eqn]->detJ; /* Really, ought to be mesh eqn. */

  F = fv->F;
  F_old = fv_old->F;
  F_0 =
      fv_dot
          ->F; /* tabaer notes that fv_dot is being misused here to store the
                * level set field from which the all-important sign is found */

  /*  F_dot_old = sqrt( fv->x[0]*fv->x[0] + fv->x[1]*fv->x[1] ) - 1.0;  */

  /*   F_dot    =  (1 + 2. * tt) * (F - F_old)/dt - 2. * tt  * F_dot_old; */

  F_dot = (F - F_old) / dt;

  /*   S =  F_0 / sqrt( F_0 * F_0 + alpha*alpha ); */

  for (a = 0; a < dim; a++) {
    if (!EXPLICIT) {
      grad_F[a] = fv->grad_F[a];
    } else {
      grad_F[a] = fv_old->grad_F[a];
    }

    grad_F_mag += grad_F[a] * grad_F[a];

    vel_mag += fv->v[a] * fv->v[a];

    grad_F0[a] = fv_dot->grad_F[a];

    grad_F0_mag += grad_F0[a] * grad_F0[a];
  }

  grad_F_mag = sqrt(grad_F_mag);
  grad_F0_mag = sqrt(grad_F0_mag);

  vel_mag = sqrt(vel_mag);

  S = F_0 / sqrt(F_0 * F_0 + (grad_F0_mag * alpha * grad_F0_mag * alpha));

  /*
   * Gradient correction
   */

  if (GRADIENT) {
    memset(w, 0, sizeof(double) * DIM);
    memset(d_w_dF, 0, sizeof(double) * DIM * MDE);

    for (a = 0; (grad_F_mag > 1.e-12) && a < dim; a++) {
      w[a] = S * grad_F[a] / grad_F_mag;

      for (j = 0; j < ei[pg->imtrx]->dof[eqn]; j++) {
        P1 = vec_dot(dim, bf[eqn]->grad_phi[j], grad_F);

        d_w_dF[a][j] =
            bf[eqn]->grad_phi[j][a] - P1 * grad_F[a] / grad_F_mag / grad_F_mag;

        d_w_dF[a][j] *= S / grad_F_mag;
      }
    }
  }

  /*
   * Streamwise correction
   */
  if (STREAMWISE) {
    for (a = 0; a < dim; a++) {
      w[a] = S * fv->v[a] / vel_mag;

      for (j = 0; j < ei[pg->imtrx]->dof[eqn]; j++) {
        d_w_dF[a][j] = 0.0;
      }
    }
  }

  /*
   * Put local contributions into global right-hand side
   * if it is not a right-hand side variable-it won't get added in (contribution
   * is zero)
   */
  if (af->Assemble_Residual) {
    for (i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
      I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
      /* check for multiple dofs */
      nvdofi = Dolphin[pg->imtrx][I][eqn];

      for (ki = 0; ki < nvdofi; ki++) {

        /* check to make sure that unknowns are defined at this node,
           otherwise don't add anything to this node */
        idof = ei[pg->imtrx]->ln_to_first_dof[eqn][i] + ki;

        /* also convert from node number to dof number */
        phi_i = bf[eqn]->phi[idof];

        for (a = 0; a < dim; a++) {
          double p = 0;
          p += upwind * alpha * w[a] * bf[eqn]->grad_phi[idof][a];

          phi_i += p;
        }

        rhs = F_dot * phi_i;
        /* try implementing an explicit/implicit method */
        /*      rhs = (F-F_old)/dt * phi_i ; */

        rhs -= S * phi_i;

        for (a = 0; a < dim; a++) {
          rhs += phi_i * w[a] * grad_F[a];
        }

        rf[node_to_fill[I] + ki] += rhs * wt * det_J * h3;
      }
    }
  }

  if (af->Assemble_Jacobian) {
    for (i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
      I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];

      nvdofi = Dolphin[pg->imtrx][I][eqn];
      for (ki = 0; ki < nvdofi; ki++) {

        ie = node_to_fill[I] + ki;
        idof = ei[pg->imtrx]->ln_to_first_dof[eqn][i] + ki;

        phi_i = bf[eqn]->phi[idof];

        for (a = 0; a < dim; a++) {
          double p = 0;
          p += upwind * alpha * w[a] * bf[eqn]->grad_phi[idof][a];

          phi_i += p;
        }

        rhs = F_dot - S + vec_dot(dim, w, grad_F);

        /* derivatives of fill equation wrt to fill variable */
        for (j = 0; j < ei[pg->imtrx]->num_local_nodes; j++) {
          J = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + j];
          nvdofj = Dolphin[pg->imtrx][J][eqn];
          for (kj = 0; kj < nvdofj; kj++)

          {

            jdof = ei[pg->imtrx]->ln_to_first_dof[eqn][j] + kj;
            phi_j = bf[eqn]->phi[jdof];

            je = node_to_fill[J] + kj;
            ja = (ie == je) ? ie : in_list(je, ijaf[ie], ijaf[ie + 1], ijaf);
            EH(ja, "Could not find vbl in sparse matrix.");

            afill[ja] += wt * h3 * det_J * phi_i * phi_j * (1 + 2. * tt) / dt;

            for (a = 0; !EXPLICIT && a < dim; a++) {
              afill[ja] += wt * h3 * det_J * phi_i *
                           (w[a] * bf[eqn]->grad_phi[jdof][a] +
                            d_w_dF[a][jdof] * grad_F[a]);
              afill[ja] += wt * h3 * det_J * rhs * upwind * alpha *
                           (d_w_dF[a][jdof] * bf[eqn]->grad_phi[idof][a]);
            }
          }
        }
      }
    }
  }

  return (status);
} /* end of assemble_level_correct */

#endif /*COUPLED_FILL*/

int assemble_level_project(
    double afill[], /* Jacobian matrix for fill equation  */
    int ijaf[],     /* pointer to nonzeros in Jacobian matrix   */
    double rf[],    /* rhs vector   */
    double dt,      /* current time step size */
    double tt,      /* parameter to vary time integration from
                     * explicit (tt = 1) to implicit (tt = 0) */
    int node_to_fill[]) {
  /*
   * Some local variables for convenience...
   */

  int eqn;

  int i, j;
  int status;
  int I, J;
  int idof, jdof;

  int ie, je, ja;
  int ki, kj, nvdofi, nvdofj;

  dbl F; /* Fill. */
  dbl F_calc;

  /*
   * Galerkin weighting functions for i-th energy residuals
   * and some of their derivatives...
   */

  dbl phi_i;

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;
  dbl h3; /* Volume element (scale factors). */
  dbl det_J;
  dbl wt;

  /* static char yo[] = "assemble_fill"; */

  status = 0;

  /*
   * Unpack variables from structures for local convenience...
   */

  eqn = R_FILL;

  /*
   * Bail out fast if there's nothing to do...
   */

  if (!pd->e[pg->imtrx][eqn]) {
    return (status);
  }

  wt = fv->wt; /* Gauss point weight. */

  h3 = fv->h3; /* Differential volume element. */

  det_J = bf[eqn]->detJ; /* Really, ought to be mesh eqn. */

  F = fv->F;

  F_calc = initial_level_set(fv->x[0], fv->x[1], fv->x[2]);

  /*
   * Put local contributions into global right-hand side
   * if it is not a right-hand side variable-it won't get added in (contribution
   * is zero)
   */
  if (af->Assemble_Residual) {
    for (i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
      I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
      /* check for multiple dofs */
      nvdofi = Dolphin[pg->imtrx][I][eqn];

      for (ki = 0; ki < nvdofi; ki++) {

        /* check to make sure that unknowns are defined at this node,
           otherwise don't add anything to this node */
        idof = ei[pg->imtrx]->ln_to_first_dof[eqn][i] + ki;

        /* also convert from node number to dof number */
        phi_i = bf[eqn]->phi[idof];

        rf[node_to_fill[I] + ki] += wt * phi_i * det_J * h3 * (F - F_calc);
      }
    }
  }

  if (af->Assemble_Jacobian) {
    for (i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
      I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];

      nvdofi = Dolphin[pg->imtrx][I][eqn];
      for (ki = 0; ki < nvdofi; ki++) {
        ie = node_to_fill[I] + ki;
        idof = ei[pg->imtrx]->ln_to_first_dof[eqn][i] + ki;
        phi_i = bf[eqn]->phi[idof];

        /* derivatives of fill equation wrt to fill variable */
        for (j = 0; j < ei[pg->imtrx]->num_local_nodes; j++) {
          J = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + j];
          nvdofj = Dolphin[pg->imtrx][J][eqn];
          for (kj = 0; kj < nvdofj; kj++)

          {

            jdof = ei[pg->imtrx]->ln_to_first_dof[eqn][j] + kj;
            phi_j = bf[eqn]->phi[jdof];

            je = node_to_fill[J] + kj;
            ja = (ie == je) ? ie : in_list(je, ijaf[ie], ijaf[ie + 1], ijaf);
            EH(ja, "Could not find vbl in sparse matrix.");

            afill[ja] += wt * h3 * det_J * phi_i * phi_j;
          }
        }
      }
    }
  }

  return (status);
} /* end of assemble_level_project */

#endif

static double initial_level_set(double x, double y, double z) {

  /* insert distance function here */

  return (0);
}

static double gradient_norm_err(double *x, Exo_DB *exo, Dpi *dpi, double range)

{

  int eb, blk_id;
  double params = 0.5 * range;
  double grad_err = 0.;
  double area = 0.;
  double error;
  double area_scale = pow(range, pd->Num_Dim);

  for (eb = 0; eb < dpi->num_elem_blocks_global; eb++) {
    int mn;
    blk_id = dpi->eb_id_global[eb];
    mn = map_mat_index(blk_id);

    if (pd_glob[mn]->e[pg->imtrx][ls->var])
      grad_err += evaluate_volume_integral(exo, dpi, I_MAG_GRAD_FILL_ERROR,
                                           NULL, blk_id, 0, NULL, &params, 1,
                                           NULL, x, x, 0.0, 0.0, 0);

    if (pd_glob[mn]->e[pg->imtrx][ls->var])
      area +=
          evaluate_volume_integral(exo, dpi, I_LS_ARC_LENGTH, NULL, blk_id, 0,
                                   NULL, &params, 1, NULL, x, x, 0.0, 0.0, 0);
  }

  if (fabs(area) > 0.005 * area_scale)
    error = sqrt(grad_err / area);
  else
    error = 0.0;

  return (error);
}

/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/

void surf_based_initialization(double *x, double *delta_x, double *xdot,
                               Exo_DB *exo, int num_total_nodes,
                               struct LS_Surf_List *list, double time,
                               double theta, double delta_t) {
  int a, I, ie;
  int DeformingMesh = upd->ep[pg->imtrx][R_MESH1];
  double r[DIM];
  double **Disp = NULL;

  struct LS_Surf *closest;
  /* NOTE: if delta_x is not NULL it is set to be the change in x with the
   * same (somewhat strange) sign convention used in mm_sol_nonlinear.c.
   * Similarly if x_dot is not NULL update it.
   * For all cases x is updated to the distance from the surfaces in list.
   */

  /* Handle special case of list of LS_SURF_NS where we need to initialize sign.
   * Note that the current implementation doesn't make sense to combine this
   * with other surfaces.
   */
  if (list->start->type == LS_SURF_NS) {
    struct LS_Surf_NS_Data *s = (struct LS_Surf_NS_Data *)list->start->data;

    if (list->start->next)
      EH(-1, "Current implementation is limited to a single NS surface");

    initialize_sign(s->PosEB_id, x, exo);
  }

  /* if we need to construct subsurfs, do it */

  create_subsurfs(list, x, exo);

  if (DeformingMesh != -1) {
    Disp = (double **)smalloc(DIM * sizeof(double *));

    for (a = 0; a < pd->Num_Dim; a++)
      Disp[a] = (double *)smalloc(num_total_nodes * sizeof(double));

    stash_node_displacements(Disp, num_total_nodes, x, exo);
  }

#define PRINT_LS 0
#if PRINT_LS
  print_surf_list(list, time);
  if (list->start->type == LS_SURF_ISOSURFACE && ls->SubElemIntegration)
    subelement_mesh_output(x, exo);
#endif

  for (I = 0; I < num_total_nodes; I++) {
    retrieve_node_coordinates(I, x, r, Disp);

    ie = Index_Solution(I, ls->var, 0, 0, -2, pg->imtrx);

    if (ie != -1) {
      closest = closest_surf(list, x, exo, r);

      /* Crude support for periodic level set function during redistancing
       */
      if (ls->Periodic_Planes) {
        double rperiodic[3], delta;
        double closest_dist = fabs(closest->closest_point->distance);
        double dist;
        double rbest[3];
        rbest[0] = r[0];
        rbest[1] = r[1];
        rbest[2] = r[2];
        /* X periodicity */
        if (ls->Periodic_Plane_Loc[0] != ls->Periodic_Plane_Loc[1]) {
          delta = ls->Periodic_Plane_Loc[0] - ls->Periodic_Plane_Loc[1];
          rperiodic[0] = r[0] + delta;
          rperiodic[1] = r[1];
          rperiodic[2] = r[2];
          closest = closest_surf(list, x, exo, rperiodic);
          dist = fabs(closest->closest_point->distance);
          if (dist < closest_dist) {
            closest_dist = dist;
            rbest[0] = rperiodic[0];
            rbest[1] = rperiodic[1];
            rbest[2] = rperiodic[2];
          }
          rperiodic[0] = r[0] - delta;
          closest = closest_surf(list, x, exo, rperiodic);
          dist = fabs(closest->closest_point->distance);
          if (dist < closest_dist) {
            closest_dist = dist;
            rbest[0] = rperiodic[0];
            rbest[1] = rperiodic[1];
            rbest[2] = rperiodic[2];
          }
        }
        /* Y periodicity */
        if (ls->Periodic_Plane_Loc[2] != ls->Periodic_Plane_Loc[3]) {
          delta = ls->Periodic_Plane_Loc[2] - ls->Periodic_Plane_Loc[3];
          rperiodic[0] = r[0];
          rperiodic[1] = r[1] + delta;
          rperiodic[2] = r[2];
          closest = closest_surf(list, x, exo, rperiodic);
          dist = fabs(closest->closest_point->distance);
          if (dist < closest_dist) {
            closest_dist = dist;
            rbest[0] = rperiodic[0];
            rbest[1] = rperiodic[1];
            rbest[2] = rperiodic[2];
          }
          rperiodic[1] = r[1] - delta;
          closest = closest_surf(list, x, exo, rperiodic);
          dist = fabs(closest->closest_point->distance);
          if (dist < closest_dist) {
            closest_dist = dist;
            rbest[0] = rperiodic[0];
            rbest[1] = rperiodic[1];
            rbest[2] = rperiodic[2];
          }
        }
        /* Z periodicity */
        if (pd->Num_Dim == 3 &&
            ls->Periodic_Plane_Loc[4] != ls->Periodic_Plane_Loc[5]) {
          delta = ls->Periodic_Plane_Loc[4] - ls->Periodic_Plane_Loc[5];
          rperiodic[0] = r[0];
          rperiodic[1] = r[1];
          rperiodic[2] = r[2] + delta;
          closest = closest_surf(list, x, exo, rperiodic);
          dist = fabs(closest->closest_point->distance);
          if (dist < closest_dist) {
            closest_dist = dist;
            rbest[0] = rperiodic[0];
            rbest[1] = rperiodic[1];
            rbest[2] = rperiodic[2];
          }
          rperiodic[2] = r[2] - delta;
          closest = closest_surf(list, x, exo, rperiodic);
          dist = fabs(closest->closest_point->distance);
          if (dist < closest_dist) {
            closest_dist = dist;
            rbest[0] = rperiodic[0];
            rbest[1] = rperiodic[1];
            rbest[2] = rperiodic[2];
          }
        }
        closest = closest_surf(list, x, exo, rbest);
      }

      /* if closest point is on a boundary on which a contact angle is
       * specified, we want to "extend" the distance function into
       * the boundary to avoid a cusp in the distance function
       */
      if (ls->Contact_Inflection) {
        struct LS_Surf_Closest_Point *cp = closest->closest_point;

        if (cp->inflection) {
          double rmirror[3];

          rmirror[0] = cp->x[0] + (cp->x[0] - r[0]);
          rmirror[1] = cp->x[1] + (cp->x[1] - r[1]);
          rmirror[2] = 0.;

          if (pd->Num_Dim == 3) {
            rmirror[2] = cp->x[2] + (cp->x[2] - r[2]);
          }
          closest = closest_surf(list, x, exo, rmirror);
          closest->closest_point->distance *= -1.;
        }
      }

      /* make sure to preserve sign for certain types of surfaces */
      if (closest->type == LS_SURF_ISOSURFACE) {
        struct LS_Surf_Iso_Data *s = (struct LS_Surf_Iso_Data *)closest->data;

        if ((x[ie] - s->isoval) * closest->closest_point->distance < 0.)
          closest->closest_point->distance *= -1;
      } else if (closest->type == LS_SURF_NS) {
        /* sign set above based on element block, just need to preserve it */
        if (x[ie] * closest->closest_point->distance < 0.)
          closest->closest_point->distance *= -1;
      } else if (closest->type == LS_SURF_POINT) {
        /* this probably shouldn't happen */
        /* assume that we want distance to keep sign */
        if (x[ie] * closest->closest_point->distance < 0.)
          closest->closest_point->distance *= -1;
      } else if (closest->type == LS_SURF_ARC) {
        int sign;

        sign = find_arc_region_sign(closest->data, r);

        closest->closest_point->distance *= sign;
      }

      if (ls != NULL && ls->Huygens_Freeze_Nodes && fabs(time) > 0) {

        int node_is_frozen = 0;
        for (int ielem = exo->node_elem_pntr[I];
             ielem < exo->node_elem_pntr[I + 1]; ielem++) {
          int elem = exo->node_elem_list[ielem];
          if (elem_on_isosurface(elem, x, exo, ls->var, 0)) {
            node_is_frozen = 1;
          }
        }
        if (!node_is_frozen) {
          if (delta_x != NULL)
            delta_x[ie] = x[ie] - closest->closest_point->distance;
          if (xdot != NULL)
            xdot[ie] -= delta_x[ie] * (1.0 + 2 * theta) / delta_t;

          x[ie] = closest->closest_point->distance;
        }

      } else {

        if (delta_x != NULL)
          delta_x[ie] = x[ie] - closest->closest_point->distance;
        if (xdot != NULL)
          xdot[ie] -= delta_x[ie] * (1.0 + 2 * theta) / delta_t;

        x[ie] = closest->closest_point->distance;
      }
    }
  }

  if (DeformingMesh != -1) {

    for (a = 0; a < pd->Num_Dim; a++)
      safe_free((void *)Disp[a]);

    safe_free((void *)Disp);
  }

  return;
}

/***************************************************************************************/
/***************************************************************************************/

void free_surf_list(struct LS_Surf_List **list_p) {
  struct LS_Surf_List *list;
  list = *list_p;

  if (list != NULL) {

    while (list->start != NULL) {
      list->current = list->start->next;

      if (list->start->subsurf_list)
        free_surf_list(&(list->start->subsurf_list));

      safer_free((void **)&(list->start->data));
      safer_free((void **)&(list->start->closest_point));

      safer_free((void **)&(list->start));

      list->start = list->current;
    }

    safer_free((void **)list_p);
  }
}

void free_subsurfs(struct LS_Surf_List *list) {
  struct LS_Surf *surf;

  if (list != NULL) {
    surf = list->start;

    while (surf) {
      if (surf->subsurf_list)
        free_surf_list(&(surf->subsurf_list));
      surf = surf->next;
    }
  }
  return;
}

struct LS_Surf *closest_surf(struct LS_Surf_List *list, double *x, Exo_DB *exo,
                             double r[DIM]) {
  struct LS_Surf *surf, *closest;
  double distance, closest_distance;
  double confidence, closest_confidence;
  double abs_closest_distance, abs_distance;
  double tol = 1.e-5;

  surf = list->start;
  find_surf_closest_point(surf, x, exo, r);
  closest = surf;
  closest_distance = closest->closest_point->distance;
  closest_confidence = closest->closest_point->confidence;

  surf = surf->next;
  while (surf != NULL) {

    find_surf_closest_point(surf, x, exo, r);
    distance = surf->closest_point->distance;
    confidence = surf->closest_point->confidence;

    /* Hopefully not too confusing here.
     * To avoid sign errors we need to permit the possibility of making
     * a very small magnitude error.  We need to catch pick the right
     * surface when two surface are nominally the same distance away
     * but one is more confident than the other about the sign of the
     * distance function.  The tolerance for error in the magnitude is
     * set in the variable tol.  If it is too small, we may pick a surface
     * that has incorrect sign information resulting in a large error in
     * the distance function.  If it is too large, we may make too large
     * or errors in the magnitude of the distance function by picking a
     * surface that is further away just cuz it is more confident about the
     * sign.
     * If we had a narrow band approach this just wouldn't matter since
     * the magnitude errors will only happen relatively far from the interface.
     */

    abs_closest_distance = fabs(closest_distance);
    abs_distance = fabs(distance);
    if (((confidence == closest_confidence) &&
         (abs_distance < abs_closest_distance)) ||
        ((confidence < closest_confidence) &&
         (abs_distance < (1. - tol) * abs_closest_distance)) ||
        ((confidence > closest_confidence) &&
         (abs_distance < (1. + tol) * abs_closest_distance))) {
      closest = surf;
      closest_distance = distance;
      closest_confidence = confidence;
    }

    surf = surf->next;
  }
  return (closest);
}

void find_surf_closest_point(struct LS_Surf *surf, double *x, Exo_DB *exo,
                             double *r) {

  struct LS_Surf_Closest_Point *cp;
  /* the surface has a closest_point structure that has the info about the
   * point on surf that is closest to r
   */

  cp = surf->closest_point;

  /* The closest point confidence is used to avoid sign errors on vertex-based
   * objects.  The default is unity which means that we are completely
   * confident of the sign of this answer.  For objects that have some
   * degree of uncertainty about the sign, the surface will be selected
   * that has the lowest distance magnitude AND the largest confidence
   */

  cp->confidence = 1.; /* default is fully confident of sign */

  /* The element number elem along with the parameter coords xi tell how to
   * interpolate quantities at the nearest point.  Some surface types, like
   * planes, spheres, and circles are mesh objects so you can't directly
   * find these quantities.  This is the default which is denoted by elem=-1.
   * For cases where this mapping is available it should be loaded to facilitate
   * evaluating quantities on the surface at the nearest point
   */

  cp->elem = -1;
  cp->elem_side = -1;

  /* The inflection parameter tells whether the closest point lies on a boundary
   * which has a CA condition applied (NOTE: the ls->Contact_Inflection must
   * also be turn on).  Currently, this is possible for facets or point based
   * surfaces.
   */

  cp->inflection = FALSE;

  switch (surf->type) {
  case LS_SURF_POINT: {
    struct LS_Surf_Point_Data *s = (struct LS_Surf_Point_Data *)surf->data;
    double *p = s->x;
    double distance;

    distance = (r[0] - p[0]) * (r[0] - p[0]) + (r[1] - p[1]) * (r[1] - p[1]);

    if (pd->Num_Dim == 3) {
      distance += (r[2] - p[2]) * (r[2] - p[2]);
    }

    /* note that sign will have to be determined elsewhere */
    cp->confidence = 0.; /* a point gives no indication of sign at all */
    cp->distance = sqrt(distance);
    cp->inflection = s->inflection;

    /* get data necessary for evaluating quantities at this point */
    cp->elem = s->elem;
    cp->xi[0] = s->xi[0];
    cp->xi[1] = s->xi[1];
    cp->xi[2] = s->xi[2];
    cp->x[0] = s->x[0];
    cp->x[1] = s->x[1];
    cp->x[2] = s->x[2];
  }

  break;

  case LS_SURF_PLANE: {
    struct LS_Surf_Plane_Data *s = (struct LS_Surf_Plane_Data *)surf->data;
    double *n = s->n;
    double d = s->d;
    double dot_prod = 0.0;

    dot_prod += r[0] * n[0] + r[1] * n[1];

    if (pd->Num_Dim == 3)
      dot_prod += r[2] * n[2];

    cp->distance = dot_prod - d;

    /* NOTE: it makes no sense to try to interpolate quantities here */

  }

  break;
  case LS_SURF_CIRCLE:
  case LS_SURF_SPHERE: {
    struct LS_Surf_Sphere_Data *s = (struct LS_Surf_Sphere_Data *)surf->data;
    int sign = ((s->r < 0.) ? -1 : 1);
    double *c = s->center;
    double radius = sign * s->r;
    double distance;

    if (surf->type == LS_SURF_SPHERE) {
      distance = (r[2] - c[2]) * (r[2] - c[2]);
    } else {
      distance = 0.0;
    }

    distance += (r[0] - c[0]) * (r[0] - c[0]) + (r[1] - c[1]) * (r[1] - c[1]);

    cp->distance = sign * (sqrt(distance) - radius);

    /* NOTE: it makes no sense to try to interpolate quantities here */

  }

  break;

  case LS_SURF_FACET: {
    int a;
    double ray[2], d[2], ray_normal[2], fraction;
    double ray_mag_squared, d_dot_ray, d_dot_n;
    double distance;

    struct LS_Surf_Facet_Data *s = (struct LS_Surf_Facet_Data *)surf->data;

    if (s->num_points == 2) {
      struct LS_Surf_Point_Data *s1 =
          (struct LS_Surf_Point_Data *)surf->subsurf_list->start->data;
      struct LS_Surf_Point_Data *s2 =
          (struct LS_Surf_Point_Data *)surf->subsurf_list->start->next->data;
      double *p1 = s1->x;
      double *p2 = s2->x;

      ray[0] = p2[0] - p1[0];
      ray[1] = p2[1] - p1[1];
      ray_mag_squared = ray[0] * ray[0] + ray[1] * ray[1];

      /* if the points p1 and p2 are oriented correctly
         (i.e. proceeding counterclockwise around body),
         ray_normal points toward positive distance
         (in my view this is an outward facing normal)
       */

      ray_normal[0] = ray[1];
      ray_normal[1] = -ray[0];
      normalize_really_simple_vector(ray_normal, 2);

      d[0] = r[0] - p1[0];
      d[1] = r[1] - p1[1];

      d_dot_ray = dot_product(2, d, ray);

      if (d_dot_ray <= 0.) /* we are closest to point p1 */
      {
        fraction = 0.;
        cp->inflection = s1->inflection;
      } else {
        if (d_dot_ray >= ray_mag_squared) /* we are closest to point p2 */
        {
          fraction = 1.;
          d[0] -= ray[0];
          d[1] -= ray[1];
          cp->inflection = s2->inflection;
        } else /* we are closest to some point on the segment */
        {
          fraction = d_dot_ray / ray_mag_squared;
          d[0] -= fraction * ray[0];
          d[1] -= fraction * ray[1];
          if (s1->inflection && s2->inflection)
            cp->inflection = TRUE;
        }
      }

      distance = sqrt(d[0] * d[0] + d[1] * d[1]);
      d_dot_n = dot_product(2, d, ray_normal);
      if (d_dot_n < 0)
        distance *= -1.;

      /* put together the mapping for this facet
       */
      cp->elem = s->elem;
      cp->elem_side = s->elem_side;
      cp->x[2] = 0.;
      cp->xi[2] = 0.;
      for (a = 0; a < pd->Num_Dim; a++) {
        cp->x[a] = (1. - fraction) * s1->x[a] + fraction * s2->x[a];
        cp->xi[a] = (1. - fraction) * s1->xi[a] + fraction * s2->xi[a];
      }

      /* confidence based on angle between ray_normal and d */
      if ((distance == 0.) || ((fraction > 0.) && (fraction < 1.))) {
        cp->confidence = 1.;
      } else {
        cp->confidence = d_dot_n / distance;
      }
      cp->distance = distance;
    } else {
      EH(-1, "Facet based surfaces not yet implemented in 3-D");
    }

  } break;

  case LS_SURF_SS:
  case LS_SURF_NS:
  case LS_SURF_ISOSURFACE:
  case LS_SURF_ARC: {
    int a;
    struct LS_Surf *subsurf;
    /* recursive call to find closest subsurface */
    subsurf = closest_surf(surf->subsurf_list, x, exo, r);

    /* inherit closest_point data from closest subsurface */
    cp->elem = subsurf->closest_point->elem;
    cp->elem_side = subsurf->closest_point->elem_side;
    for (a = 0; a < 3; a++) {
      cp->x[a] = subsurf->closest_point->x[a];
      cp->xi[a] = subsurf->closest_point->xi[a];
    }
    cp->distance = subsurf->closest_point->distance;
    cp->confidence = subsurf->closest_point->confidence;
    cp->inflection = subsurf->closest_point->inflection;
  }

  break;

  case LS_SURF_USER: {
    struct LS_Surf_User_Data *s = (struct LS_Surf_User_Data *)surf->data;

    cp->distance = user_surf_object(s->Int_Data, s->Real_Data, r);

  } break;

  default:
    break;
  }
}

static struct LS_Surf_List *create_surfs_from_ns(int ns_id, double *x,
                                                 Exo_DB *exo)
/* construct list of surface objects from node set */
{
  /* first make list of surfaces from specified node sets */
  int i, ns, I, inflection;
  int ns_num_nodes;
  int *ns_node_list = NULL;
  double p[3], xi[3];
  struct LS_Surf_List *list;
  struct LS_Surf *surf;

  list = create_surf_list();

  ns = 0;
  while (exo->ns_id[ns] != ns_id && ns++ < exo->num_node_sets)
    ;

  if (ns >= exo->num_node_sets) {
    if (Num_Proc == 1)
      EH(-1, "Error: cannot find initial level set nodeset for Huygens "
             "initialization \n");
  } else {
    ns_num_nodes = exo->ns_num_nodes[ns];

    ns_node_list = exo->ns_node_list + exo->ns_node_index[ns];

    for (i = 0; i < ns_num_nodes; i++) {
      I = ns_node_list[i];

      memset(p, 0, sizeof(double) * 3);
      memset(xi, 0, sizeof(double) * 3);

      retrieve_node_coordinates(I, x, p, NULL);

      inflection = point_on_ca_boundary(I, exo);

      /* Not supplying elem or xi here since we don't really have elem
       * accessible. */
      surf = create_surf_point(p, -1, xi, inflection);
      append_surf(list, surf);
    }
  }

  return (list);
}

static struct LS_Surf_List *create_surfs_from_ss(int ss_id, double *x,
                                                 Exo_DB *exo)
/* construct list of surface objects from side set */
{
  /* first make list of surfaces from specified node sets */
  int i, n, iss;
  int num_nodes_on_side, ielem, side, ielem_type;
  int nodes_per_side, ielem_dim;
  int iconnect_ptr;
  int local_elem_node_id[MAX_NODES_PER_SIDE];
  int I, inflection, ln, lnn[2];
  double p[3], xi[3];
  struct LS_Surf_List *list;
  struct LS_Surf *surf, *vertex[2];

  /* find the side set index for this side set */
  if ((iss = in_list(ss_id, 0, exo->num_side_sets, exo->ss_id)) == -1) {
    fprintf(stderr,
            "Error in create_subsurfs_from_ss.  Cannot find side set id %d",
            ss_id);
    exit(-1);
  }

  list = create_surf_list();

  /* loop over all elements in this SS */
  for (i = 0; i < exo->ss_num_sides[iss]; i++) {

    ielem = exo->ss_elem_list[exo->ss_elem_index[iss] + i];

    side = exo->ss_side_list[exo->ss_elem_index[iss] + i];


    /* Get the number of nodes on the side of the current element
          NOTE - Currently, this must be the same for all elements
          in the side set, a fairly big limitation        */

    num_nodes_on_side =
        (exo->ss_node_side_index[iss][i + 1] - exo->ss_node_side_index[iss][i]);

    ielem_type = Elem_Type(exo, ielem);
    get_side_info(ielem_type, side, &nodes_per_side, local_elem_node_id);
    iconnect_ptr = exo->elem_ptr[ielem];
    ielem_dim = elem_info(NDIM, ielem_type);

    for (n = 0; n < num_nodes_on_side - 1; n++) {

      if (ielem_dim == 2) {

        if (num_nodes_on_side == 2) {
          lnn[0] = local_elem_node_id[0];
          lnn[1] = local_elem_node_id[1];
        } else if (num_nodes_on_side == 3) {
          if (n == 0) {
            lnn[0] = local_elem_node_id[0];
            lnn[1] = local_elem_node_id[2];
          } else if (n == 1) {
            lnn[0] = local_elem_node_id[2];
            lnn[1] = local_elem_node_id[1];
          }
        }

        for (ln = 0; ln < 2; ln++) {
          I = exo->node_list[iconnect_ptr + lnn[ln]];

          retrieve_node_coordinates(I, x, p, NULL);

          inflection = point_on_ca_boundary(I, exo);
          find_nodal_stu(lnn[ln], ielem_type, xi, xi + 1, xi + 2);

          vertex[ln] = create_surf_point(p, ielem, xi, inflection);
        }

        surf = create_surf_facet_line(vertex[0], vertex[1], ielem, side);
        append_surf(list, surf);

      } else {
        EH(-1, "create_subsurfs_from_ss currently only implemented for 2-D");
      }
    }
  } /* for ( i=0; i<exo->ss_num_sides[iss]; i++) */

  return (list);
}

int stash_node_displacements(double **d, int num_total_nodes, double *x,
                             Exo_DB *exo) {
  int *moved;
  int e_start, e_end, ielem = 0;
  int i, var, p, ln, gnn;
  int dim = pd->Num_Dim;

  double phi[MDE];

  moved = (int *)smalloc(num_total_nodes * sizeof(int));
  memset(moved, 0, num_total_nodes * sizeof(int));

  e_start = exo->eb_ptr[0];
  e_end = exo->eb_ptr[exo->num_elem_blocks];

  for (ielem = e_start; ielem < e_end; ielem++) {
    load_elem_dofptr(ielem, exo, x, x, x, x, 1);

    for (ln = 0; ln < ei[pg->imtrx]->num_local_nodes; ln++) {
      double xi[3] = {0.0, 0.0, 0.0};

      gnn = exo->elem_node_list[exo->elem_node_pntr[ielem] + ln];

      find_nodal_stu(ln, ei[pg->imtrx]->ielem_type, xi, xi + 1, xi + 2);

      if (moved[gnn] != 1) {
        for (p = 0; p < dim; p++) {
          d[p][gnn] = 0.0;
          var = MESH_DISPLACEMENT1 + p;

          if (pd->v[pg->imtrx][var]) {

            for (i = 0; i < ei[pg->imtrx]->dof[var]; i++) {
              phi[i] = newshape(xi, ei[pg->imtrx]->ielem_type, PSI,
                                ei[pg->imtrx]->dof_list[var][i],
                                ei[pg->imtrx]->ielem_shape,
                                pd->i[pg->imtrx][var], i);

              d[p][gnn] += *esp->d[p][i] * phi[i];
            }
            moved[gnn] = 1;
          }
        }
      }
    }
  }
  safe_free((void *)moved);
  return (0);
}

static void retrieve_node_coordinates(int I, double *x, double *p, double **d) {
  int a, ie, var;

  p[0] = Coor[0][I];
  p[1] = Coor[1][I];
  p[2] = (pd->Num_Dim == 3) ? Coor[2][I] : 0.0;

  if (d == NULL) {
    for (a = 0; a < pd->Num_Dim; a++) {
      var = MESH_DISPLACEMENT1 + a;

      if (pd->v[pg->imtrx][var]) {
        ie = Index_Solution(I, var, 0, 0, -1, pg->imtrx);

        if (ie != -1)
          p[a] += x[ie];
      }
    }
  } else {
    p[0] += d[0][I];
    p[1] += d[1][I];
    p[2] += (pd->Num_Dim == 3) ? d[2][I] : 0.0;
  }

  return;
}

static struct LS_Surf_List *create_surfs_from_arc(double c[DIM], double r,
                                                  double n[DIM], double d) {

  struct LS_Surf_List *list;

  double theta = 0;
  int nint = 1200;
  double del_theta = 360.0 / nint;
  double test_point[DIM] = {0., 0., 0.};

  list = create_surf_list();

  while (theta < 360.0) {
    struct LS_Surf *surf;
    struct LS_Surf_Point_Data *s;

    /*  stuff that gets done */
    test_point[0] = c[0] + r * cos(M_PIE * theta / 180.0);
    test_point[1] = c[1] + r * sin(M_PIE * theta / 180.0);

    if (knockout_point(test_point, n, d) == FALSE) {
      surf = create_surf(LS_SURF_POINT);
      s = (struct LS_Surf_Point_Data *)surf->data;

      s->x[0] = test_point[0];
      s->x[1] = test_point[1];
      s->x[2] = 0.0;

      if (unique_surf(list, surf)) {
        append_surf(list, surf);
      } else {
        safe_free(surf);
      }
    }
    theta += del_theta;
  }
  return (list);
}

static int knockout_point(double test[DIM], double normal[DIM], double d) {
  double R;

  R = test[0] * normal[0] + test[1] * normal[1];

  return (R > d);
}

static double find_arc_region_sign(struct LS_Surf_Arc_Data *s, double *r) {

  double *c = s->center;
  double *n = s->n;
  double radius = s->r;
  double d = s->d;
  double R;

  double sign = s->sign;

  if ((r[0] * n[0] + r[1] * n[1]) < d) {
    R = pow((r[0] - c[0]), 2) + pow((r[1] - c[1]), 2);

    if (R > radius * radius) {
      sign *= -1.0;
    }
  }
  return (sign);
}

static struct LS_Surf_List *create_surfs_from_iso(int isovar, double isoval,
                                                  double x[], Exo_DB *exo) {
  int ebi, e, e_start, e_end;

  struct LS_Surf_List *list;
  struct Shape_Fcn_Tree *ntree = NULL;
  int isosurf_exists = FALSE;

  list = create_surf_list();

  for (ebi = 0; ebi < exo->num_elem_blocks; ebi++) {

    ei[pg->imtrx]->elem_blk_id = exo->eb_id[ebi];

    pd = pd_glob[Matilda[ebi]];
    mp = mp_glob[Matilda[ebi]];

    e_start = exo->eb_ptr[ebi];
    e_end = exo->eb_ptr[ebi + 1];

    if (pd->e[pg->imtrx][isovar]) {
      for (e = e_start; e < e_end; e++) {

        if (elem_on_isosurface(e, x, exo, isovar, isoval)) {

          load_elem_dofptr(e, exo, x_static, x_old_static, xdot_static,
                           xdot_old_static, 0);

          bf_mp_init(pd);

          switch (ls->Isosurface_Subsurf_Type) {

          case LS_SURF_POINT:

            if (ls->Search_Option == GRID_SEARCH) {
              SGRID *sgrid;

              if (ntree == NULL) {
                ntree = create_shape_fcn_tree(ls->Grid_Search_Depth);
                /* 			      print_shape_fcn_tree ( ntree ); */
              }

              sgrid = create_search_grid(ntree);

              divide_search_grid(sgrid, ls->Grid_Search_Depth);

              /* 			print_search_grid ( sgrid ); */
              find_grid_intersections(sgrid, list);

              free_search_grid(&sgrid);
              sgrid = NULL;
            } else {
              find_intersections(list, isovar, isoval, exo);
            }
            break;

          case LS_SURF_FACET:
            if (ls->SubElemIntegration && isovar == LS) {
              get_subelement_facets(list, isoval);
            } else {
              find_facets(list, isovar, isoval, exo);
            }
            break;

          default:
            EH(-1, "Unknown subsurface type for isosurface reconstruction");
            break;
          }
        }
      }
    }
  }

  if (ntree != NULL) {
    free_shape_fcn_tree(ntree);
    ntree = NULL;
  }

  isosurf_exists = (list->size == 0) ? FALSE : TRUE;

  if (Num_Proc > 1) {
#ifdef PARALLEL
    int global_bool = FALSE;
    MPI_Allreduce(&isosurf_exists, &global_bool, 1, MPI_INT, MPI_LOR,
                  MPI_COMM_WORLD);
    isosurf_exists = global_bool;
#endif
  }

  if (isosurf_exists == FALSE)
    smooth_stop_with_msg("\n\tUnable to locate points on zero level set "
                         "isosurface.  Stopping Goma.\n");

  return (list);
}

void create_subsurfs(struct LS_Surf_List *list, double *x, Exo_DB *exo) {
  struct LS_Surf *surf = list->start;

  while (surf != NULL) {

    switch (surf->type) {
    case LS_SURF_SS: {
      struct LS_Surf_SS_Data *s = (struct LS_Surf_SS_Data *)surf->data;

      /* delete old list first */
      free_surf_list(&(surf->subsurf_list));

      surf->subsurf_list = create_surfs_from_ss(s->ss_id, x, exo);
    } break;

    case LS_SURF_NS: {
      struct LS_Surf_NS_Data *s = (struct LS_Surf_NS_Data *)surf->data;

      /* delete old list first */
      free_surf_list(&(surf->subsurf_list));

      surf->subsurf_list = create_surfs_from_ns(s->ns_id, x, exo);
    } break;

    case LS_SURF_ISOSURFACE: {
      struct LS_Surf_Iso_Data *s = (struct LS_Surf_Iso_Data *)surf->data;

      /* delete old list first */
      free_surf_list(&(surf->subsurf_list));

      surf->subsurf_list = create_surfs_from_iso(s->isovar, s->isoval, x, exo);
    }

    break;
    case LS_SURF_ARC: {
      struct LS_Surf_Arc_Data *s = (struct LS_Surf_Arc_Data *)surf->data;
      free_surf_list(&(surf->subsurf_list));

      surf->subsurf_list = create_surfs_from_arc(s->center, s->r, s->n, s->d);

      print_surf_list(surf->subsurf_list, 0.0);

    } break;
    default:
      break;
    }

    if (Num_Proc > 1) {
      if (surf->subsurf_list != NULL)
        assemble_Global_surf_list(surf->subsurf_list);
    }

    surf = surf->next;
  }

  return;
}
static void initialize_sign(int PosEB_id, double *x, Exo_DB *e) {
  int ie, i, I;
  int eb;
  int num_nodes_per_blk;

  /*  if( e->num_elem_blocks != 2 && Num_Proc == 1 )
      {
        EH(-1,"Huygens Initialization only functions for two element blocks
     \n");
      }
  */
  for (eb = 0; eb < e->num_elem_blocks; eb++) {
    num_nodes_per_blk = e->eb_num_nodes_per_elem[eb] * e->eb_num_elems[eb];

    for (i = 0; i < num_nodes_per_blk; i++) {
      I = e->eb_conn[eb][i];
      ie = Index_Solution(I, ls->var, 0, 0, -2, pg->imtrx);
      if (ie != -1)
        x[ie] = e->eb_id[eb] == PosEB_id ? 1.0 : -1.0;
    }
  }
}

#ifdef PARALLEL
static void ddd_add_surf(DDD pkg, struct LS_Surf *surf) {
  ddd_add_member(pkg, &(surf->type), 1, MPI_INT);

  switch (surf->type) {
  case LS_SURF_POINT: {
    struct LS_Surf_Point_Data *s = (struct LS_Surf_Point_Data *)surf->data;

    ddd_add_member(pkg, s->x, 3, MPI_DOUBLE);
    ddd_add_member(pkg, &(s->elem), 1, MPI_INT);
    ddd_add_member(pkg, s->xi, 3, MPI_DOUBLE);
    ddd_add_member(pkg, &(s->inflection), 1, MPI_INT);
  } break;

  case LS_SURF_PLANE: {
    struct LS_Surf_Plane_Data *s = (struct LS_Surf_Plane_Data *)surf->data;

    ddd_add_member(pkg, s->n, 3, MPI_DOUBLE);
    ddd_add_member(pkg, &(s->d), 1, MPI_DOUBLE);
  } break;

  case LS_SURF_CIRCLE:
  case LS_SURF_SPHERE: {
    struct LS_Surf_Sphere_Data *s = (struct LS_Surf_Sphere_Data *)surf->data;

    ddd_add_member(pkg, s->center, 3, MPI_DOUBLE);
    ddd_add_member(pkg, &(s->r), 1, MPI_DOUBLE);
  } break;
  case LS_SURF_ARC: {
    struct LS_Surf_Arc_Data *s = (struct LS_Surf_Arc_Data *)surf->data;
    ddd_add_member(pkg, s->center, 3, MPI_DOUBLE);
    ddd_add_member(pkg, &(s->r), 1, MPI_DOUBLE);
    ddd_add_member(pkg, s->n, 3, MPI_DOUBLE);
    ddd_add_member(pkg, &(s->d), 1, MPI_DOUBLE);
    ddd_add_member(pkg, &(s->sign), 1, MPI_DOUBLE);
  } break;
  case LS_SURF_FACET: {
    struct LS_Surf_Facet_Data *s = (struct LS_Surf_Facet_Data *)surf->data;

    ddd_add_member(pkg, &(s->num_points), 1, MPI_INT);
  } break;

  case LS_SURF_SS: {
    struct LS_Surf_SS_Data *s = (struct LS_Surf_SS_Data *)surf->data;

    ddd_add_member(pkg, &(s->ss_id), 1, MPI_INT);
  } break;

  case LS_SURF_NS: {
    struct LS_Surf_NS_Data *s = (struct LS_Surf_NS_Data *)surf->data;

    ddd_add_member(pkg, &(s->ns_id), 1, MPI_INT);
    ddd_add_member(pkg, &(s->PosEB_id), 1, MPI_INT);
  } break;

  case LS_SURF_ISOSURFACE: {
    struct LS_Surf_Iso_Data *s = (struct LS_Surf_Iso_Data *)surf->data;

    ddd_add_member(pkg, &(s->isovar), 1, MPI_INT);
    ddd_add_member(pkg, &(s->isoval), 1, MPI_DOUBLE);
  } break;

  case LS_SURF_USER: {
    struct LS_Surf_User_Data *s = (struct LS_Surf_User_Data *)surf->data;

    ddd_add_member(pkg, s->Int_Data, 5, MPI_INT);
    ddd_add_member(pkg, s->Real_Data, 10, MPI_DOUBLE);
  } break;

  default: {
    EH(-1, "ddd_add_surf called with unknown surface type");
  } break;
  }
}

#endif

void assemble_Global_surf_list(struct LS_Surf_List *list) {
#ifdef PARALLEL

  DDD pkg = NULL;

  int *list_size = NULL;

  struct LS_Surf *surf = NULL;

  int i, n;
  int my_list_size;
  int *levels, *types, level;

#endif

  if (list == NULL)
    EH(-1, "Lists must exist before calling assemble_Global_surf_list");

#ifdef PARALLEL
  list_size = (int *)smalloc(Num_Proc * sizeof(int));

  /* recurse list to find how many surfaces I have */
  /* To start, set list->current to NULL so that next_surf_or_subsurf
   * knows to return first entry in the list.
   */
  list->current = NULL;
  my_list_size = 0;
  do {
    level = 0;
    surf = next_surf_or_subsurf(list, &level);
    if (surf)
      my_list_size++;
  } while (surf);

  /* find out how much everybody is going to tell us
   */
  MPI_Allgather(&my_list_size, 1, MPI_INT, list_size, 1, MPI_INT,
                MPI_COMM_WORLD);

  /* handle one processor at a time */
  for (i = 0; i < Num_Proc; i++) {
    if (list_size[i] > 0) {
      /* get (or create) initial data concerning type and number of subsurfs */
      levels = (int *)smalloc(list_size[i] * sizeof(int));
      types = (int *)smalloc(list_size[i] * sizeof(int));

      pkg = ddd_alloc();

      if (ProcID == i) {
        /* Traverse my original list of surfs getting the level and type. */

        /* again reset the list */
        list->current = NULL;

        for (n = 0; n < list_size[i]; n++) {
          level = 0;
          surf = next_surf_or_subsurf(list, &level);

          if (surf == NULL) {
            EH(-1, "That really shouldn't happen.");
          }

          levels[n] = level;
          types[n] = surf->type;
        }
      }

      /* now broadcast this info */
      MPI_Bcast(levels, list_size[i], MPI_INT, i, MPI_COMM_WORLD);
      MPI_Bcast(types, list_size[i], MPI_INT, i, MPI_COMM_WORLD);

      /* now both sender and receivers march through shoving addresses
       * via ddd_add_member
       */
      if (ProcID == i) {
        /* again reset the list */
        list->current = NULL;
      }

      for (n = 0; n < list_size[i]; n++) {
        if (ProcID == i) {
          level = 0;
          surf = next_surf_or_subsurf(list, &level);
        } else {
          surf = create_next_surf_or_subsurf(list, levels[n], types[n]);
        }

        ddd_add_surf(pkg, surf);
      }

      safer_free((void **)&levels);
      safer_free((void **)&types);

      ddd_set_commit(pkg);
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(MPI_BOTTOM, 1, pkg->new_type, i, MPI_COMM_WORLD);
      ddd_free(pkg);
    }
  }

  if (list->size == 0 && Num_Proc == 1) {
    EH(-1, "No points found on zero level set");
    /*       DPRINTF(stderr,"No points found on zero level set\n"); */
  }

  safer_free((void **)&list_size);
#endif
}

/* (f1 < 0. ) ? (f2 >= 0.) : (f2 < 0.) */
int sign_change(double f1, double f2) {
  if (f1 < 0.) {
    if (f2 < 0.)
      return FALSE;
    return TRUE;
  } else {
    if (f2 < 0.)
      return TRUE;
    return FALSE;
  }
}

/*
 * Following examines and element ielem as to whether the isoval
 * contour of the isovar variable is present in the element.
 */

int elem_on_isosurface(int elem, double x[], const Exo_DB *exo, int isovar,
                       double isoval) {
  int i, I;
  int iconn_ptr = exo->elem_ptr[elem];
  int dofs;
  double f[MDE];
  int ielem_type, ielem_shape, interpolation;
  int ebn, mn;

  ebn = find_elemblock_index(elem, exo);
  mn = Matilda[ebn];
  if (!(pd_glob[mn]->v[pg->imtrx][isovar]))
    return (FALSE);

  ielem_type = Elem_Type(exo, elem);
  ielem_shape = type2shape(ielem_type);
  interpolation = pd_glob[mn]->i[pg->imtrx][isovar];
  dofs = getdofs(ielem_shape, interpolation);

  I = exo->node_list[iconn_ptr + 0];

  f[0] = x[Index_Solution(I, isovar, 0, 0, -2, pg->imtrx)] - isoval;

  /* if there is a zero crossing, somebody must have a different
   * sign than node 0
   */
  for (i = 1; i < dofs; i++) {
    I = exo->node_list[iconn_ptr + i];

    f[i] = x[Index_Solution(I, isovar, 0, 0, -2, pg->imtrx)] - isoval;

    if (sign_change(f[i], f[0])) {
      return (TRUE);
    }
  }

/* DRN - this is somewhat experimental */
#if 0 
  /* look for crossings that don't intersect a side but don't change the sign
     of any nodes
   */
  switch ( interpolation ) {
    case I_Q1:
      return FALSE;
    case I_Q2:
      switch ( ielem_shape ) {
	case QUADRILATERAL:
	  {
	    int iside;
	    int nodes_per_side;
	    int lnn[3];
	    for ( iside=0; iside<4; iside++ )
              {
                get_side_info( ielem_type, iside+1, &nodes_per_side, lnn);
		if ( fabs(f[lnn[2]]) < 0.25*fabs(f[lnn[0]]+f[lnn[1]]) - 0.5*sqrt(f[lnn[0]]*f[lnn[1]]) )
		  {
#if 0
		    DPRINTF(stderr,"Yikes! Zero crossing in elem_on_isosurface for ielem=%d with all nodes on one side!!!\n",elem);
#endif
		    return TRUE;
		  }
	      }
	    return FALSE;
	  }
      }
  }
#endif

  return (FALSE);
}

int elem_near_isosurface(int elem, double x[], const Exo_DB *exo, int isovar,
                         double isoval) {
  int i, I;
  int iconn_ptr = exo->elem_ptr[elem];
  int dofs;
  double f[MDE];
  int ielem_type, ielem_shape, interpolation;
  int ebn, mn;

  ebn = find_elemblock_index(elem, exo);
  mn = Matilda[ebn];
  if (!(pd_glob[mn]->v[pg->imtrx][isovar]))
    return (FALSE);

  ielem_type = Elem_Type(exo, elem);
  ielem_shape = type2shape(ielem_type);
  interpolation = pd_glob[mn]->i[pg->imtrx][isovar];
  dofs = getdofs(ielem_shape, interpolation);

  I = exo->node_list[iconn_ptr + 0];

  f[0] = x[Index_Solution(I, isovar, 0, 0, -2, pg->imtrx)] - isoval;

  dbl hsquared[DIM];
  dbl hh[DIM][DIM];
  dbl dhh_dxnode[DIM][MDE];

  h_elem_siz(hsquared, hh, dhh_dxnode, pd_glob[mn]->gv[R_MESH1]);

  double h_elem = 0;
  for (int a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {
    h_elem += hsquared[a];
  }
  /* This is the size of the element */
  h_elem = sqrt(h_elem / ((double)ei[pg->imtrx]->ielem_dim));
  /* if there is a zero crossing, somebody must have a different
   * sign than node 0
   */
  for (i = 1; i < dofs; i++) {
    I = exo->node_list[iconn_ptr + i];

    f[i] = x[Index_Solution(I, isovar, 0, 0, -2, pg->imtrx)] - isoval;

    if (fabs(f[i]) < (h_elem * 0.5)) {
      return (TRUE);
    }
  }

  return (FALSE);
}

/*
 * Following examines an element, elem, and computes the
 * average value of the given var.
 */
/*
double
element_average ( int elem,
                  double x[],
                  const Exo_DB* exo,
                  int var )
{
  int i, I;
  int iconn_ptr = exo->elem_ptr[elem];
  int dofs;
  int ielem_type, ielem_shape;
  int ebn, mn;
  double sum = 0.;

  ebn = find_elemblock_index(elem, exo);
  mn = Matilda[ebn];
  if (!(pd_glob[mn]->v[pg->imtrx][var])) return (0.);

  ielem_type = Elem_Type(exo, elem);
  ielem_shape  = type2shape(ielem_type);
  dofs = getdofs(ielem_shape, pd_glob[mn]->i[pg->imtrx][var]);

  I = exo->node_list[ iconn_ptr + 0 ];

  for ( i=0 ; i< dofs; i++)
    {
      I = exo->node_list[ iconn_ptr + i ];

      sum += x[Index_Solution(I, var, 0, 0, -2, pg->imtrx)];
    }

  return (sum/dofs);
}
*/
/*
 * Following examines current element as to whether the isoval
 * contour of the isovar variable is present in the element.
 */

int current_elem_on_isosurface(int isovar, double isoval) {
  double *esp;
  double f[MDE];
  int i;

  if (!(pd->gv[isovar]))
    return (FALSE);

  esp = x_static +
        ei[pg->imtrx]->ieqn_ledof[ei[pg->imtrx]->lvdof_to_ledof[isovar][0]];
  f[0] = *esp - isoval;

  /* if there is a zero crossing, somebody must have a different
   * sign than node 0
   */
  for (i = 1; i < ei[pg->imtrx]->dof[isovar]; i++) {
    esp = x_static +
          ei[pg->imtrx]->ieqn_ledof[ei[pg->imtrx]->lvdof_to_ledof[isovar][i]];
    f[i] = *esp - isoval;

    if (sign_change(f[i], f[0])) {
      return (TRUE);
    }
  }

/* DRN - this is somewhat experimental */
#if 1
  /* look for crossings that intersect a side but don't change the sign
     of any nodes
   */
  switch (pd->i[pg->imtrx][isovar]) {
  case I_Q1:
    return FALSE;
  case I_Q2:
    switch (ei[pg->imtrx]->ielem_shape) {
    case QUADRILATERAL: {
      int iside;
      int nodes_per_side;
      int lnn[3];
      for (iside = 0; iside < 4; iside++) {
        get_side_info(ei[pg->imtrx]->ielem_type, iside + 1, &nodes_per_side,
                      lnn);
        if (fabs(f[lnn[2]]) < 0.25 * fabs(f[lnn[0]] + f[lnn[1]]) -
                                  0.5 * sqrt(f[lnn[0]] * f[lnn[1]])) {
#if 0
		    DPRINTF(stderr,"Yikes! Zero crossing in current_elem_on_isosurface for ielem=%d with all nodes on one side!!!\n",ei[pg->imtrx]->ielem);
#endif
          return TRUE;
        }
      }
      return FALSE;
    }
    }
  }
#endif

  return (FALSE);
}
/*
#define BASE_ELEM_SIG_CROSS_TOL 1.e-12
int
significant_element_crossing ( int elem,
                               double x[],
                               const Exo_DB* exo )
{
  int i, I;
  int iconn_ptr = exo->elem_ptr[elem];
  int dofs;
  int ielem_type, ielem_shape;
  int ebn, mn;
  double value;
  double max_neg = 0.;
  double max_pos = 0.;

  ebn = find_elemblock_index(elem, exo);
  mn = Matilda[ebn];
  if (!(pd_glob[mn]->v[pg->imtrx][ls->var])) return (FALSE);

  ielem_type = Elem_Type(exo, elem);
  ielem_shape  = type2shape(ielem_type);
  dofs = getdofs(ielem_shape, pd_glob[mn]->i[pg->imtrx][ls->var]);

  for ( i=0 ; i< dofs; i++)
    {
      I = exo->node_list[ iconn_ptr + i ];

      value = x[Index_Solution(I, ls->var, 0, 0, -2, pg->imtrx)];

      if ( -value > max_neg ) max_neg = -value;
      if ( value > max_pos ) max_pos = value;
    }

  if ( max_neg * max_pos == 0. ) return (FALSE);

  if ( max_neg > max_pos )
    {
      if ( max_pos/max_neg > BASE_ELEM_SIG_CROSS_TOL ) return (TRUE);
    }
  else
    {
      if ( max_neg/max_pos > BASE_ELEM_SIG_CROSS_TOL ) return (TRUE);
    }

  return (FALSE);
}
*/
/*
static int
significant_current_element_crossing ()
{
  int i;
  double value;
  double max_neg = 0.;
  double max_pos = 0.;

  if (!(pd->v[pg->imtrx][ls->var])) return (FALSE);

  for ( i=0 ; i< ei[pg->imtrx]->dof[ls->var]; i++)
    {
          switch (ls->var )
          {
                  case FILL:
                          value =  *(esp->F[i]);
                          break;
                  case PHASE1:
                  case PHASE2:
                  case PHASE3:
                  case PHASE4:
                  case PHASE5:
                          value = *(esp->pF[ls->var-PHASE1][i]);
                          break;
          }

      if ( -value > max_neg ) max_neg = -value;
      if ( value > max_pos ) max_pos = value;
    }

  if ( max_neg * max_pos == 0. ) return (FALSE);

  if ( max_neg > max_pos )
    {
      if ( max_pos/max_neg > BASE_ELEM_SIG_CROSS_TOL ) return (TRUE);
    }
  else
    {
      if ( max_neg/max_pos > BASE_ELEM_SIG_CROSS_TOL ) return (TRUE);
    }

  return (FALSE);
}
*/

static int Hrenorm_simplemass(Exo_DB *exo, Comm_Ex *cx, Dpi *dpi, double x[],
                              struct LS_Surf_List *list, int num_total_nodes,
                              int num_ls_unkns, int num_total_unkns,
                              double time) {
  int I, ie, max_its, i;
  double *dC;
  double *F, *F_, *R, *b;
  double M;
  int global_ls_unkns = num_ls_unkns;
  static double M0;
  static int M0_set = 0;
#ifdef PARALLEL

  double interface_size = 0;
  double c = 0;
  int eb, blk_id;
  double Mold;

  MPI_Allreduce(&num_ls_unkns, &global_ls_unkns, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);

#endif

  dC = F = F_ = R = b = NULL;

  dalloc(num_total_unkns, dC);

  dalloc(num_ls_unkns, F);
  dalloc(num_ls_unkns, F_);
  dalloc(num_ls_unkns, R);
  dalloc(num_ls_unkns, b);

  if ((!M0_set) && (ls->Mass_Value >= 0.0)) {
    if (ls->Mass_Value > 0) {
      M0 = ls->Mass_Value;
    } else {
      M0 = find_LS_mass(exo, dpi, NULL, dC, x, num_total_unkns);
      M0_set = 1;
    }
  } else if (!M0_set) {
    M0 = find_LS_mass(exo, dpi, NULL, dC, x, num_total_unkns);
  }

  surf_based_initialization(x, NULL, NULL, exo, num_total_nodes, list, time, 0.,
                            0.);
#ifdef PARALLEL
  exchange_dof(cx, dpi, x, pg->imtrx);
#endif

  max_its = 20;
  Mold = M0;
  M = find_LS_mass(exo, dpi, NULL, dC, x, num_total_unkns);
  DPRINTF(stdout, "\n\t\t Mass old %g, Mass new %g: \t", M0, M);

  for (i = 0; i < max_its; i++) {
    if (fabs(M - Mold) < 1e-7) {
      break;
    }

    interface_size = 0;
    for (eb = 0; eb < dpi->num_elem_blocks_global; eb++) {
      blk_id = dpi->eb_id_global[eb];

      interface_size +=
          evaluate_volume_integral(exo, dpi, I_LS_ARC_LENGTH, NULL, blk_id, 0,
                                   NULL, NULL, 1, NULL, x, x, 0.0, 0.0, 0);
    }

    c = (M - Mold) / interface_size;
    for (I = 0; I < num_total_nodes; I++) {
      if ((ie = Index_Solution(I, ls->var, 0, 0, -2, pg->imtrx)) != -1) {
        x[ie] += c;
      }
    }
#ifdef PARALLEL
    exchange_dof(cx, dpi, x, pg->imtrx);
#endif

    Mold = M0;
    M = find_LS_mass(exo, dpi, NULL, dC, x, num_total_unkns);

    DPRINTF(stdout, "\n\t\t iter %d, Additive value: %f, new mass %g \n\t\t", i,
            c, M);
  }

#ifdef PARALLEL
  exchange_dof(cx, dpi, x, pg->imtrx);
#endif

  safe_free(dC);
  safe_free(F);
  safe_free(F_);
  safe_free(R);
  safe_free(b);

  return (TRUE);
}

static int Hrenorm_smolianksi_only(Exo_DB *exo, Comm_Ex *cx, Dpi *dpi,
                                   double x[], struct LS_Surf_List *list,
                                   int num_total_nodes, int num_ls_unkns,
                                   int num_total_unkns, double time) {
  int I, ie, max_its, i;
  double *dC;
  double *F, *F_, *R, *b;
  double M;
  int global_ls_unkns = num_ls_unkns;
  static double M0;
  static int M0_set = 0;
#ifdef PARALLEL

  double interface_size = 0;
  double c = 0;
  int eb, blk_id;
  double Mold;

  MPI_Allreduce(&num_ls_unkns, &global_ls_unkns, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);

#endif

  dC = F = F_ = R = b = NULL;

  dalloc(num_total_unkns, dC);

  dalloc(num_ls_unkns, F);
  dalloc(num_ls_unkns, F_);
  dalloc(num_ls_unkns, R);
  dalloc(num_ls_unkns, b);

  if ((!M0_set) && (ls->Mass_Value >= 0.0)) {
    if (ls->Mass_Value > 0) {
      M0 = ls->Mass_Value;
    } else {
      M0 = find_LS_mass(exo, dpi, NULL, dC, x, num_total_unkns);
      M0_set = 1;
    }
  } else if (!M0_set) {
    M0 = find_LS_mass(exo, dpi, NULL, dC, x, num_total_unkns);
  }
#ifdef PARALLEL
  exchange_dof(cx, dpi, x, pg->imtrx);
#endif

  max_its = 20;
  Mold = M0;
  M = find_LS_mass(exo, dpi, NULL, dC, x, num_total_unkns);
  DPRINTF(stdout, "\n\t\t Mass old %g, Mass new %g: \t", M0, M);

  for (i = 0; i < max_its; i++) {
    if (fabs(M - Mold) < 1e-7) {
      break;
    }

    interface_size = 0;
    for (eb = 0; eb < dpi->num_elem_blocks_global; eb++) {
      blk_id = dpi->eb_id_global[eb];

      interface_size +=
          evaluate_volume_integral(exo, dpi, I_LS_ARC_LENGTH, NULL, blk_id, 0,
                                   NULL, NULL, 1, NULL, x, x, 0.0, 0.0, 0);
    }

    c = (M - Mold) / interface_size;
    for (I = 0; I < num_total_nodes; I++) {
      if ((ie = Index_Solution(I, ls->var, 0, 0, -2, pg->imtrx)) != -1) {
        x[ie] += c;
      }
    }
#ifdef PARALLEL
    exchange_dof(cx, dpi, x, pg->imtrx);
#endif

    Mold = M0;
    M = find_LS_mass(exo, dpi, NULL, dC, x, num_total_unkns);

    DPRINTF(stdout, "\n\t\t iter %d, Additive value: %f, new mass %g \n\t\t", i,
            c, M);
  }

#ifdef PARALLEL
  exchange_dof(cx, dpi, x, pg->imtrx);
#endif

  safe_free(dC);
  safe_free(F);
  safe_free(F_);
  safe_free(R);
  safe_free(b);

  return (TRUE);
}

static int Hrenorm_constrain(Exo_DB *exo, Comm_Ex *cx, Dpi *dpi, double x[],
                             struct LS_Surf_List *list, int num_total_nodes,
                             int num_ls_unkns, int num_total_unkns,
                             double time) {
  int I, ie, k, max_its;
  double M0, lamda, R_lamda, norm = 1.0;
  double *dC;
  double *F, *F_, *R, *b;
  int *ie_map;
  double bTb, bTR, RTR, M;
  double d_lamda;
  int global_ls_unkns = num_ls_unkns;
#ifdef PARALLEL
  int *ext_dof = NULL;
  double local_bTb = 0.;
  double local_bTR = 0.;
  double local_RTR = 0.;

  double global_bTb = 0.;
  double global_bTR = 0.;
  double global_RTR = 0.;

  MPI_Allreduce(&num_ls_unkns, &global_ls_unkns, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);

#endif

  dC = F = F_ = R = b = NULL;

  dalloc(num_total_unkns, dC);

  dalloc(num_ls_unkns, F);
  dalloc(num_ls_unkns, F_);
  dalloc(num_ls_unkns, R);
  dalloc(num_ls_unkns, b);

  DPRINTF(stdout, "\n\t\t Mass constraint iteration: \t");

  ie_map = (int *)smalloc(num_ls_unkns * sizeof(int));

#ifdef PARALLEL
  if (Num_Proc > 1)
    ext_dof = (int *)smalloc(num_ls_unkns * sizeof(int));
#endif

  for (I = 0, k = 0; I < num_total_nodes; I++) {
    if ((ie = Index_Solution(I, ls->var, 0, 0, -2, pg->imtrx)) != -1) {
      F[k] = x[ie];
      ie_map[k] = ie;

#ifdef PARALLEL
      /* Here we set up the array ext_dof which */
      if (Num_Proc > 1) {
        if (I > (dpi->num_internal_nodes + dpi->num_boundary_nodes - 1))
          ext_dof[k] = TRUE;
        else
          ext_dof[k] = FALSE;
      }
#endif

      k++;
    }
  }

  if (k != num_ls_unkns) {
    EH(-1, "Error in Hrenorm_constrain. Level set unknowns unaccounted for.");
  }

  M0 = find_LS_mass(exo, dpi, NULL, dC, x, num_total_unkns);

  if (ls->Mass_Value != 0.0)
    M0 = ls->Mass_Value;

  surf_based_initialization(x, NULL, NULL, exo, num_total_nodes, list, time, 0.,
                            0.);

  for (k = 0; k < num_ls_unkns; k++) {
    F_[k] = x[ie_map[k]];
    x[ie_map[k]] = F[k];
  }

  R_lamda = lamda = 0.0;

  max_its = 10;

  while (norm >= 1.e-6 && max_its > 0) {

    for (k = 0, bTb = bTR = RTR = 0.0; k < num_ls_unkns; k++) {
      int external = FALSE;

      b[k] = dC[ie_map[k]];
      R[k] = 2.0 * (F[k] - F_[k]) + lamda * b[k];

#ifdef PARALLEL
      if (Num_Proc > 1)
        external = ext_dof[k];
        /* when accumulating the dot products it is important not to add the
         * external dofs */
#endif

      if (!external) {
        bTb += b[k] * b[k];
        bTR += b[k] * R[k];
        RTR += R[k] * R[k];
      }
    }

#ifdef PARALLEL
    if (Num_Proc > 1) {
      local_bTb = bTb;
      local_bTR = bTR;
      local_RTR = RTR;

      MPI_Allreduce(&local_bTb, &global_bTb, 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD);
      bTb = global_bTb;

      MPI_Allreduce(&local_bTR, &global_bTR, 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD);
      bTR = global_bTR;

      MPI_Allreduce(&local_RTR, &global_RTR, 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD);

      RTR = global_RTR;
    }

#endif
    norm = sqrt(fabs(RTR)) / global_ls_unkns + fabs(R_lamda);

    d_lamda = (2.0 * R_lamda - bTR) / bTb;

    for (k = 0; k < num_ls_unkns; k++) {
      F[k] += -0.5 * R[k] - 0.5 * b[k] * d_lamda;

      x[ie_map[k]] = F[k];
    }
#ifdef PARALLEL
    exchange_dof(cx, dpi, x, pg->imtrx);
#endif

    lamda += d_lamda;

    M = find_LS_mass(exo, dpi, NULL, dC, x, num_total_unkns);

    R_lamda = M - M0;

    max_its--;
    DPRINTF(stdout, ". %e/%e ", M, M0);
  }

  DPRINTF(stdout, "\n\t\t Multiplier value: %f \n\t\t", lamda);

  safe_free(dC);
  safe_free(F);
  safe_free(F_);
  safe_free(R);
  safe_free(b);

  return (TRUE);
}

static double find_LS_mass(const Exo_DB *exo, const Dpi *dpi,
                           const double *params, double *dC, double x[],
                           int num_total_unkns) {
  double M = 0.0;
  int eb, blk_id;
  int num_params = params == NULL ? 0 : 1;

  if (dC != NULL) {
    memset(dC, 0, num_total_unkns * sizeof(double));
  }

  for (eb = 0, M = 0.0; eb < dpi->num_elem_blocks_global; eb++) {
    int mn;
    blk_id = dpi->eb_id_global[eb];
    mn = map_mat_index(blk_id);

    if (pd_glob[mn]->e[pg->imtrx][ls->var])
      M += evaluate_volume_integral(exo, dpi, ls->Mass_Sign, NULL, blk_id, 0,
                                    NULL, params, num_params, dC, x, x, 0.0,
                                    tran->time_value, 0);
  }

  return (M);
}

double find_LS_global_flux(const Exo_DB *exo, const Dpi *dpi,
                           const double *params, double *dC, double x[],
                           int num_total_unkns) {

  double M = 0.0;
  int eb, blk_id;

  if (dC != NULL)
    memset(dC, 0, num_total_unkns * sizeof(double));

  for (eb = 0, M = 0.0; eb < dpi->num_elem_blocks_global; eb++) {
    int mn;
    blk_id = dpi->eb_id_global[eb];
    mn = map_mat_index(blk_id);

    if (pd_glob[mn]->e[pg->imtrx][ls->var])

      M += evaluate_global_flux(
          exo, dpi, (ls->Mass_Sign == I_NEG_FILL ? NEG_LS_FLUX : POS_LS_FLUX),
          blk_id, 0, NULL, NULL, x, 0.0, 0);
  }

  return (M);
}

/******************************************************************************
 * find_LS_vel: computes the average of a velocity component in either the
 * positive or negative phases (w.r.t. level set FILL).  The 'chosen_vel'
 * parameter should be one of I_{NEG,POS}_V{X,Y,Z}.  I copied the find_LS_mass
 * function for this....
 *
 * Input
 * =====
 * exo                = ExodusII pointer.
 * dpi                = Comuncations pointer.
 * params             = If non-NULL, params[0] is the interfacial width to use
 *                      in lieu of ls->Length_Scale. NB:
 *compute_volume_integrand() uses H(F) so this really doesn't matter. chosen_vel
 *= Specifies which velocity component and phase to average. x = The solution
 *vector. num_total_unknowns = Number of total unknowns.
 *
 * Returns
 * =======
 * VV                 = Average velocity component over a given phase.
 *
 ******************************************************************************/
double find_LS_vel(const Exo_DB *exo, const Dpi *dpi, const double *params,
                   const int chosen_vel, double x[], int num_total_unkns) {
  double Vel, Vol;
  int eb, blk_id, mn, phase;
  int num_params = params == NULL ? 0 : 1;

  if (chosen_vel == I_NEG_VX || chosen_vel == I_NEG_VY ||
      chosen_vel == I_NEG_VZ)
    phase = I_NEG_FILL;
  else
    phase = I_POS_FILL;

  Vel = 0.0;
  Vol = 0.0;
  for (eb = 0; eb < dpi->num_elem_blocks_global; eb++) {
    blk_id = dpi->eb_id_global[eb];
    mn = map_mat_index(blk_id);

    if (pd_glob[mn]->e[pg->imtrx][R_LEVEL_SET]) {
      Vel +=
          evaluate_volume_integral(exo, dpi, chosen_vel, NULL, blk_id, 0, NULL,
                                   params, num_params, NULL, x, x, 0.0, 0.0, 0);
      Vol +=
          evaluate_volume_integral(exo, dpi, phase, NULL, blk_id, 0, NULL,
                                   params, num_params, NULL, x, x, 0.0, 0.0, 0);
    }
  }

  return (Vel / Vol);
}

void print_point_list(double *x, Exo_DB *exo, char *filename, double time_value)

{

  struct LS_Surf_List *list;
  struct LS_Surf *surf;

  int default_Isosurface_Subsurf_Type = ls->Isosurface_Subsurf_Type;

  FILE *ofp = NULL;

  if (ProcID == 0 && filename != NULL) /* Only Proc 0 writes */
  {

    if ((ofp = fopen(filename, "w")) == NULL) {
      EH(-1, "Error opening level set output file.\n");
    }
  }

  list = create_surf_list();

  append_surf_isosurf(list, LS, 0.0);

  ls->Isosurface_Subsurf_Type =
      LS_SURF_POINT; /* We are only interested in points at the moment */

  create_subsurfs(list, x, exo);

  surf = list->start->subsurf_list->start;

  list->current = surf;

  DPRINTF(ofp, "Time: %11.5e          X\t          Y\t          Z\n",
          time_value);

  while (surf != NULL) {
    struct LS_Surf_Point_Data *s = (struct LS_Surf_Point_Data *)surf->data;

    DPRINTF(ofp, "                 %11.5e\t%11.5e\t%11.5e\n", s->x[0], s->x[1],
            s->x[2]);

    list->current = surf->next;
    surf = list->current;
  }

  ls->Isosurface_Subsurf_Type = default_Isosurface_Subsurf_Type;

  free_surf_list(&list);

  if (ProcID == 0)
    fclose(ofp);
}

int generate_facet_list(double (**point0)[DIM], double (**point1)[DIM],
                        int **owning_elem, double *x, Exo_DB *exo)

{
  struct LS_Surf_List *list;
  struct LS_Surf_List *facet_list;
  struct LS_Surf *surf, *surf_point;
  int num_facets, count;
  struct LS_Surf_Point_Data *pt;

  int default_Isosurface_Subsurf_Type = ls->Isosurface_Subsurf_Type;

  list = create_surf_list();

  append_surf_isosurf(list, LS, 0.0);

  ls->Isosurface_Subsurf_Type =
      LS_SURF_FACET; /* We are interested in a facet description of interface */

  create_subsurfs(list, x, exo);

  facet_list = list->start->subsurf_list;

  surf = facet_list->start;

  num_facets = 0;

  while (surf != NULL) {
    if (surf->type != LS_SURF_FACET) {
      EH(-1, "Error traversing facet list\n");
    }
    num_facets++;
    surf = surf->next;
  }

  if (*owning_elem != NULL) {
    safe_free((void *)*point0);
    safe_free((void *)*point1);
    safe_free((void *)*owning_elem);
  }

  *point0 = (double(*)[DIM])smalloc(DIM * num_facets * sizeof(double));
  *point1 = (double(*)[DIM])smalloc(DIM * num_facets * sizeof(double));
  *owning_elem = (int *)smalloc(num_facets * sizeof(int));

  count = 0;
  surf = facet_list->start;
  while (surf != NULL) {
    if (surf->type != LS_SURF_FACET) {
      EH(-1, "Error traversing facet list\n");
    }

    struct LS_Surf_Facet_Data *f = (struct LS_Surf_Facet_Data *)surf->data;

    if (f->num_points != 2) {
      EH(-1, "Only linear facets in 2-D are supported.\n");
    }

    (*owning_elem)[count] = f->elem;

    surf_point = surf->subsurf_list->start;

    if (surf_point->type != LS_SURF_POINT) {
      EH(-1, "Error traversing facet point list\n");
    }

    pt = (struct LS_Surf_Point_Data *)surf_point->data;
    (*point0)[count][0] = pt->x[0];
    (*point0)[count][1] = pt->x[1];

    pt = (struct LS_Surf_Point_Data *)surf_point->next->data;
    (*point1)[count][0] = pt->x[0];
    (*point1)[count][1] = pt->x[1];

    count++;
    surf = surf->next;
  }

  ls->Isosurface_Subsurf_Type = default_Isosurface_Subsurf_Type;

  free_surf_list(&list);

  return num_facets;
}

int
print_ls_interface( double *x,
		    Exo_DB *exo,
		    Dpi    *dpi,
		    const double time,
		    char *filenm,
		    int print_all_times )
{
  char output_filenm[MAX_FNL];
  struct LS_Surf_List *list = NULL;
  struct LS_Surf *isosurf = NULL;
  struct LS_Surf_Iso_Data *s;
  FILE *outfile = NULL;
  int status = 0;

  strncpy(output_filenm, filenm, MAX_FNL-1);
  multiname(output_filenm, ProcID, Num_Proc);

  if (print_all_times) {
    outfile = fopen(output_filenm, "a");
  } else {
    outfile = fopen(output_filenm, "w");
  }

  if (outfile == NULL) {
    EH(-1, "Output file for level set interface could not be opened");
  }

  list = create_surf_list();
  isosurf = create_surf( LS_SURF_ISOSURFACE );
  s = (struct LS_Surf_Iso_Data *) isosurf->data;
  s->isovar = ls->var;
  if ( ls->Initial_LS_Displacement != 0. )
    {
      s->isoval = ls->Initial_LS_Displacement;
      ls->Initial_LS_Displacement = 0.;
    }
  else
    {
      s->isoval = 0.;
    }


  append_surf( list, isosurf );

  create_subsurfs( list, x, exo );

  struct LS_Surf *surf;
  surf = list->start->subsurf_list->start;
  while (surf != NULL) {

    switch (surf->type) {
    case LS_SURF_POINT :
      {
	struct LS_Surf_Point_Data *s = (struct LS_Surf_Point_Data *) surf->data;
	double *p = s->x;
	if (print_all_times) {
	  fprintf(outfile, "%g\t%g\t%g\t%d\n", time, p[0], p[1], 0);
	} else {
	  fprintf(outfile, "%g\t%g\t%d\n", p[0], p[1], 0);
	}
      }
      break;
    case LS_SURF_FACET :
      {
	struct LS_Surf_Facet_Data *s = (struct LS_Surf_Facet_Data *) surf->data;

	if (s->num_points == 2)
	  {
	    struct LS_Surf_Point_Data *s1 = (struct LS_Surf_Point_Data *)
	      surf->subsurf_list->start->data;
	    struct LS_Surf_Point_Data *s2 = (struct LS_Surf_Point_Data *)
	      surf->subsurf_list->start->next->data;
	    double *p1 = s1->x;
	    double *p2 = s2->x;
	    if (print_all_times) {
	      fprintf(outfile, "%g\t%g\t%g\t%d\n", time, p1[0], p1[1], 0);
	      fprintf(outfile, "%g\t%g\t%g\t%d\n", time, p2[0], p2[1], 1);
	    } else {
	      fprintf(outfile, "%g\t%g\t%d\n", p1[0], p1[1], 0);
	      fprintf(outfile, "%g\t%g\t%d\n", p2[0], p2[1], 1);
	    }
	  }
	else
	  {
	    EH(-1,"Facet based surfaces not yet implemented in 3-D");
	  }
      }
      break;
    default:
      EH(-1, "Cannot print level set interfaces that are not Points or Facets");
      break;
    }

    surf = surf->next;
  }

  free_surf_list ( &list );
  fclose(outfile);
  return(status);
}

void print_surf_list(struct LS_Surf_List *list, double time) {
  int i = 0;
  FILE *f;
  static FILE *g = NULL;
  struct LS_Surf *surf = list->start;
  int j, level;
  char filename1[80], filename2[80], err_msg[MAX_CHAR_ERR_MSG];

#ifdef PARALLEL
  if (Num_Proc > 1) {
    sprintf(filename1, "level_set%d_of_%d.dat", ProcID + 1, Num_Proc);
    sprintf(filename2, "level_set_all%d_of_%d.dat", ProcID + 1, Num_Proc);
  } else {
    sprintf(filename1, "level_set.dat");
    sprintf(filename2, "level_set_all.dat");
  }
#else
  sprintf(filename1, "level_set.dat");
  sprintf(filename2, "level_set_all.dat");
#endif

  if ((f = fopen(filename1, "w")) == NULL) {
    sprintf(err_msg, "Error opening %s\n", filename1);
    EH(-1, err_msg);
  }

  if (g == NULL) {
    if ((g = fopen(filename2, "w")) == NULL) {
      sprintf(err_msg, "Error opening %s\n", filename2);
      EH(-1, err_msg);
    }
  }

  /* recurse entire list including subsurfs */
  /* To start, set list->current to NULL so that next_surf_or_subsurf
   * knows to return first entry in the list.
   */
  list->current = NULL;

  level = 0;
  surf = next_surf_or_subsurf(list, &level);

  while (surf) {
    /* add tabs to indicate surface hierarchy */
    for (j = 0; j < level; j++) {
      fprintf(f, "\t");
      fprintf(g, "\t");
    }

    switch (surf->type) {
    case LS_SURF_POINT: {

#ifdef LS_SURFPOINT_FILE
      struct LS_Surf_Point_Data *s = (struct LS_Surf_Point_Data *)surf->data;
      fprintf(f, "POINT %d\t %f\t %f\t %f\t %d\n", i, s->x[0], s->x[1], s->x[2],
              s->inflection);
      fprintf(g, "POINT %f\t %d\t %f\t %f\t %f\t %d\n", time, i, s->x[0],
              s->x[1], s->x[2], s->inflection);
      if ((pd->Num_Dim != 3) && (fabs(s->x[2]) > 1.e-8)) {
        printf("huh?\n");
      }
#endif
    }

    break;

    case LS_SURF_PLANE: {
      struct LS_Surf_Plane_Data *s = (struct LS_Surf_Plane_Data *)surf->data;

      fprintf(f, "PLANE %d\t %f\t %f\t %f\t %f\n", i, s->n[0], s->n[1], s->n[2],
              s->d);
      fprintf(g, "PLANE %f\t %d\t %f\t %f\t %f\t %f\n", time, i, s->n[0],
              s->n[1], s->n[2], s->d);

    }

    break;
    case LS_SURF_CIRCLE:
    case LS_SURF_SPHERE: {
      struct LS_Surf_Sphere_Data *s = (struct LS_Surf_Sphere_Data *)surf->data;

      fprintf(f, "CICLE/SPHERE %d\t %f\t %f\t %f\t %f\n", i, s->center[0],
              s->center[1], s->center[2], s->r);
      fprintf(g, "CICLE/SPHERE %f\t %d\t %f\t %f\t %f\t %f\n", time, i,
              s->center[0], s->center[1], s->center[2], s->r);
    }

    break;

    case LS_SURF_FACET: {
      struct LS_Surf_Facet_Data *s = (struct LS_Surf_Facet_Data *)surf->data;

      fprintf(f, "FACET %d\t %d\n", i, s->num_points);
      fprintf(g, "FACET %f\t %d\t %d\t\n", time, i, s->num_points);
    }

    break;

    case LS_SURF_SS: {
      struct LS_Surf_SS_Data *s = (struct LS_Surf_SS_Data *)surf->data;

      fprintf(f, "SS %d\t %d\n", i, s->ss_id);
      fprintf(g, "SS %f\t %d\t %d\n", time, i, s->ss_id);
    } break;

    case LS_SURF_NS: {
      struct LS_Surf_NS_Data *s = (struct LS_Surf_NS_Data *)surf->data;

      fprintf(f, "NS %d\t %d\t %d\n", i, s->ns_id, s->PosEB_id);
      fprintf(g, "NS %f\t %d\t %d\t %d\n", time, i, s->ns_id, s->PosEB_id);
    } break;

    case LS_SURF_ISOSURFACE: {
      struct LS_Surf_Iso_Data *s = (struct LS_Surf_Iso_Data *)surf->data;

      fprintf(f, "Isosurface %d\t %d\t %f\n", i, s->isovar, s->isoval);
      fprintf(g, "Isosurface %f\t %d\t %d\t %f\n", time, i, s->isovar,
              s->isoval);
    }

    break;

    default:
      break;
    }

    i++;

    level = 0;
    surf = next_surf_or_subsurf(list, &level);
  }
  fclose(f);
}

static void find_intersections(struct LS_Surf_List *list, int isovar,
                               double isoval, Exo_DB *exo) {
  double xi[3] = {0., 0., 0.};
  double yi[3] = {0., 0., 0.};
  double x[3] = {0., 0., 0.};
  int i, j, I, J, link;
  int inflection;
  struct LS_Surf *surf;

  switch (ei[pg->imtrx]->ielem_shape) {

  case QUADRILATERAL:
  case SHELL:

    switch (pd->i[pg->imtrx][isovar]) {

    case I_Q1: /* bilinear quadrilateral */
    {
      double links[6][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 2}, {1, 3}};

      for (link = 0; link < 6; link++) {

        i = links[link][0];
        j = links[link][1];
        I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
        J = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + j];

        find_nodal_stu(i, ei[pg->imtrx]->ielem_type, xi, xi + 1, xi + 2);
        find_nodal_stu(j, ei[pg->imtrx]->ielem_type, yi, yi + 1, yi + 2);

        if (find_link_intersection(xi, yi, isovar, isoval, NULL) == TRUE) {
          /* check if crossing is on an edge with a ca condition */
          inflection = FALSE; /* innocent till proven guilty */
          if (ls->Contact_Inflection && point_on_ca_boundary(I, exo) &&
              point_on_ca_boundary(J, exo))
            inflection = TRUE;
          map_local_coordinates(xi, x);
          surf = create_surf_point(x, ei[pg->imtrx]->ielem, xi, inflection);
          if (unique_surf(list, surf)) {
            append_surf(list, surf);
          } else {
            safe_free(surf);
          }
        }
      }
    } break;

    case I_Q2: /* biquadratic quadrilateral */
    {
      double links[20][2] = {{0, 4}, {4, 1}, {1, 5}, {5, 2}, {2, 6},
                             {6, 3}, {3, 7}, {7, 0}, {4, 8}, {8, 6},
                             {7, 8}, {8, 5}, {4, 7}, {0, 8}, {4, 5},
                             {1, 8}, {2, 8}, {5, 6}, {6, 7}, {3, 8}};

      for (link = 0; link < 20; link++) {

        i = links[link][0];
        j = links[link][1];
        I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
        J = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + j];

        find_nodal_stu(i, ei[pg->imtrx]->ielem_type, xi, xi + 1, xi + 2);
        find_nodal_stu(j, ei[pg->imtrx]->ielem_type, yi, yi + 1, yi + 2);

        if (find_link_intersection(xi, yi, isovar, isoval, NULL) == TRUE) {
          /* check if crossing is on an edge with a ca condition */
          inflection = FALSE; /* innocent till proven guilty */
          if (ls->Contact_Inflection && point_on_ca_boundary(I, exo) &&
              point_on_ca_boundary(J, exo))
            inflection = TRUE;
          map_local_coordinates(xi, x);
          surf = create_surf_point(x, ei[pg->imtrx]->ielem, xi, inflection);
          if (unique_surf(list, surf)) {
            append_surf(list, surf);
          } else {
			free(surf->data);
			free(surf->closest_point);
			free(surf);
          }
        }
      }
    } break;

    default:
      EH(-1, "Huygens renormalization not implemented for interpolation");
      break;
    }
    break;

  case HEXAHEDRON:

    switch (pd->i[pg->imtrx][isovar]) {

    case I_Q1: /* trilinear hex */
    {
      double links[28][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {1, 5},
                             {2, 6}, {3, 7}, {4, 5}, {5, 6}, {6, 7}, {7, 4},
                             {0, 2}, {1, 3}, {0, 5}, {1, 4}, {1, 6}, {2, 5},
                             {3, 6}, {2, 7}, {0, 7}, {3, 4}, {4, 6}, {5, 7},
                             {0, 6}, {1, 7}, {2, 4}, {3, 5}};

      for (link = 0; link < 28; link++) {

        i = links[link][0];
        j = links[link][1];
        I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
        J = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + j];

        find_nodal_stu(i, ei[pg->imtrx]->ielem_type, xi, xi + 1, xi + 2);
        find_nodal_stu(j, ei[pg->imtrx]->ielem_type, yi, yi + 1, yi + 2);

        if (find_link_intersection(xi, yi, isovar, isoval, NULL) == TRUE) {
          /* check if crossing is on an edge with a ca condition */
          inflection = FALSE; /* innocent till proven guilty */
          if (ls->Contact_Inflection && point_on_ca_boundary(I, exo) &&
              point_on_ca_boundary(J, exo))
            inflection = TRUE;
          map_local_coordinates(xi, x);
          surf = create_surf_point(x, ei[pg->imtrx]->ielem, xi, inflection);
          if (unique_surf(list, surf)) {
            append_surf(list, surf);
          } else {
            safe_free(surf);
          }
        }
      }
    } break;

    case I_Q2: /* triquadratic hex */
    {
      double links[86][2] = {
          {0, 11},  {11, 21}, {21, 8},  {8, 0},   {12, 23}, {23, 20}, {20, 25},
          {25, 12}, {11, 23}, {21, 20}, {8, 25},  {0, 12},  {11, 3},  {3, 10},
          {10, 21}, {23, 15}, {15, 26}, {26, 20}, {3, 15},  {10, 26}, {21, 9},
          {9, 1},   {1, 8},   {20, 24}, {24, 13}, {13, 25}, {9, 24},  {1, 13},
          {10, 3},  {2, 9},   {26, 14}, {14, 24}, {2, 19},  {4, 19},  {19, 22},
          {22, 16}, {16, 4},  {12, 4},  {23, 19}, {20, 22}, {25, 16}, {19, 7},
          {7, 18},  {18, 22}, {15, 7},  {26, 18}, {22, 17}, {17, 5},  {5, 16},
          {24, 17}, {13, 5},  {18, 6},  {6, 17},  {19, 6},  {0, 20},  {12, 21},
          {23, 8},  {11, 25}, {11, 26}, {3, 20},  {10, 23}, {21, 15}, {8, 24},
          {21, 13}, {9, 25},  {1, 20},  {21, 14}, {10, 24}, {3, 20},  {9, 26},
          {12, 22}, {23, 16}, {20, 4},  {25, 19}, {23, 18}, {15, 22}, {26, 19},
          {20, 7},  {25, 17}, {20, 5},  {24, 16}, {13, 22}, {20, 6},  {26, 17},
          {14, 22}, {24, 18}}; /* This set of links is not all the links that
                                  exist on a 27node brick */

      for (link = 0; link < 86; link++) {

        i = links[link][0];
        j = links[link][1];
        I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
        J = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + j];

        find_nodal_stu(i, ei[pg->imtrx]->ielem_type, xi, xi + 1, xi + 2);
        find_nodal_stu(j, ei[pg->imtrx]->ielem_type, yi, yi + 1, yi + 2);

        if (find_link_intersection(xi, yi, isovar, isoval, NULL) == TRUE) {
          /* check if crossing is on an edge with a ca condition */
          inflection = FALSE; /* innocent till proven guilty */
          if (ls->Contact_Inflection && point_on_ca_boundary(I, exo) &&
              point_on_ca_boundary(J, exo))
            inflection = TRUE;
          map_local_coordinates(xi, x);
          surf = create_surf_point(x, ei[pg->imtrx]->ielem, xi, inflection);
          if (unique_surf(list, surf)) {
            append_surf(list, surf);
          } else {
            safe_free(surf);
          }
        }
      }
    } break;

    default:
      EH(-1, "Huygens renormalization not implemented for interpolation");
      break;
    }
    break;

  case TRIANGLE:
  case TRISHELL:

    switch (pd->i[pg->imtrx][isovar]) {

    case I_Q1: /* bilinear triangular shell */
    case I_Q2: {
      double links[3][2] = {{0, 1}, {1, 2}, {2, 0}};

      for (link = 0; link < 3; link++) {

        i = links[link][0];
        j = links[link][1];
        I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
        J = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + j];

        find_nodal_stu(i, ei[pg->imtrx]->ielem_type, xi, xi + 1, xi + 2);
        find_nodal_stu(j, ei[pg->imtrx]->ielem_type, yi, yi + 1, yi + 2);

        if (find_link_intersection(xi, yi, isovar, isoval, NULL) == TRUE) {
          /* check if crossing is on an edge with a ca condition */
          inflection = FALSE; /* innocent till proven guilty */
          if (ls->Contact_Inflection && point_on_ca_boundary(I, exo) &&
              point_on_ca_boundary(J, exo))
            inflection = TRUE;
          map_local_coordinates(xi, x);
          surf = create_surf_point(x, ei[pg->imtrx]->ielem, xi, inflection);
          if (unique_surf(list, surf)) {
            append_surf(list, surf);
          } else {
            safe_free(surf);
          }
        }
      }
    } break;

    default:
      EH(-1, "Huygens renormalization not implemented for this interpolation "
             "on TRIs");
      break;
    }
    break;

  case TETRAHEDRON:

    switch (pd->i[pg->imtrx][isovar]) {

    case I_Q1: /* trilinear hex */
    {
      int links[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1,3}, {2,3}};

      for (link = 0; link < 6; link++) {

        i = links[link][0];
        j = links[link][1];
        I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
        J = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + j];

        find_nodal_stu(i, ei[pg->imtrx]->ielem_type, xi, xi + 1, xi + 2);
        find_nodal_stu(j, ei[pg->imtrx]->ielem_type, yi, yi + 1, yi + 2);

        if (find_link_intersection(xi, yi, isovar, isoval, NULL) == TRUE) {
          /* check if crossing is on an edge with a ca condition */
          inflection = FALSE; /* innocent till proven guilty */
          if (ls->Contact_Inflection && point_on_ca_boundary(I, exo) &&
              point_on_ca_boundary(J, exo))
            inflection = TRUE;
          map_local_coordinates(xi, x);
          surf = create_surf_point(x, ei[pg->imtrx]->ielem, xi, inflection);
          if (unique_surf(list, surf)) {
            append_surf(list, surf);
          } else {
            safe_free(surf);
          }
        }
      }
    } break;
    default:
      EH(-1, "Huygens renormalization not implemented for this interpolation "
             "on TRIs");
      break;
    }
    break;

  default:
    EH(-1, "Huygens renormalization not implemented for this element shape");
    break;
  }
}

void find_facets(struct LS_Surf_List *list, int isovar, double isoval,
                 Exo_DB *exo) {

  switch (ei[pg->imtrx]->ielem_shape) {

  case QUADRILATERAL:
  case SHELL:

    switch (pd->i[pg->imtrx][isovar]) {
    case I_Q1: /* bilinear quadrilateral */
    {

      int quad[4] = {0, 1, 2, 3};

      find_quad_facets(list, isovar, isoval, quad, exo);

    } break;

    case I_Q2: /* biquadratic quadrilateral */
    {
      int quads[4][4] = {
          {0, 4, 8, 7}, {4, 1, 5, 8}, {8, 5, 2, 6}, {7, 8, 6, 3}};

      find_quad_facets(list, isovar, isoval, quads[0], exo);
      find_quad_facets(list, isovar, isoval, quads[1], exo);
      find_quad_facets(list, isovar, isoval, quads[2], exo);
      find_quad_facets(list, isovar, isoval, quads[3], exo);

    } break;

    default:
      printf("isovar=%d, pd->i[pg->imtrx][isovar]=%d\n", isovar,
             pd->i[pg->imtrx][isovar]);
      EH(-1, "Facet based contouring not implemented for quads with this "
             "interpolation type");
      break;
    }
    break;

  case HEXAHEDRON:
    switch (pd->i[pg->imtrx][isovar]) {

    default:
      EH(-1, "Facet based contouring not implemented for hexes with this "
             "interpolation type");
      break;
    }
    break;

  default:
    EH(-1, "Facet based contouring not implemented for this element shape");
    break;
  }
}

static void find_quad_facets(struct LS_Surf_List *list, int isovar,
                             double isoval, int *quad, Exo_DB *exo) {
  double xi[3], yi[3], x[3], *esp;
  int i, j, I, J;
  int edge, vert_count, edge_sign[4], inflection;
  struct LS_Surf *vertex[4], *surf;
  double val_center, val_I, val_J;

  vert_count = 0;
  val_center = 0.;
  for (edge = 0; edge < 4; edge++) {
    i = quad[edge];
    j = quad[(edge + 1) % 4];
    I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
    J = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + j];

    esp = x_static +
          ei[pg->imtrx]->ieqn_ledof[ei[pg->imtrx]->lvdof_to_ledof[isovar][i]];
    val_I = *esp - isoval;

    esp = x_static +
          ei[pg->imtrx]->ieqn_ledof[ei[pg->imtrx]->lvdof_to_ledof[isovar][j]];
    val_J = *esp - isoval;

    /* keep track of average value for resolving ambiguous configuration */
    val_center += val_I;

    /* an edge is considered to have a crossing if the function
     * changes sign along the edge OR if the first node on the
     * side is zero
     */
    if ((val_I * val_J < 0.) || ((val_I == 0.) && (val_J != 0.))) {
      if (val_J > val_I) {
        edge_sign[vert_count] = 1;
      } else if (val_J < val_I) {
        edge_sign[vert_count] = -1;
      } else {
        edge_sign[vert_count] = 0;
      }

      /* locate crossing */
      find_nodal_stu(i, ei[pg->imtrx]->ielem_type, xi, xi + 1, xi + 2);
      find_nodal_stu(j, ei[pg->imtrx]->ielem_type, yi, yi + 1, yi + 2);
      if (!find_link_intersection(xi, yi, isovar, isoval, NULL)) {
        EH(-1, "Shouldn't be here.");
      }

      /* check if crossing is on an edge with a ca condition */
      inflection = FALSE; /* innocent till proven guilty */
      if (ls->Contact_Inflection && point_on_ca_boundary(I, exo) &&
          point_on_ca_boundary(J, exo))
        inflection = TRUE;

      map_local_coordinates(xi, x);
      vertex[vert_count] =
          create_surf_point(x, ei[pg->imtrx]->ielem, xi, inflection);

      vert_count++;
    }
  }

  switch (vert_count) {

  case 0:
  case 1: {
    return;
  }
    /* no "break;" to avoid the "statement not reached" warning */

  case 2: {
    if (edge_sign[0] > edge_sign[1]) {
      surf = create_surf_facet_line(vertex[0], vertex[1], ei[pg->imtrx]->ielem,
                                    -1);
    } else {
      surf = create_surf_facet_line(vertex[1], vertex[0], ei[pg->imtrx]->ielem,
                                    -1);
    }
    append_surf(list, surf);
  } break;

  case 4: {
    /* Slightly ambiguous situation of all 4 edges having crossings.
     * To handle this, make the center have the same sign it currently
     * has.
     */

    if (edge_sign[0] * val_center < 0.) {
      if (edge_sign[0] > edge_sign[1]) {
        surf = create_surf_facet_line(vertex[0], vertex[1],
                                      ei[pg->imtrx]->ielem, -1);
      } else {
        surf = create_surf_facet_line(vertex[1], vertex[0],
                                      ei[pg->imtrx]->ielem, -1);
      }
      append_surf(list, surf);

      if (edge_sign[2] > edge_sign[3]) {
        surf = create_surf_facet_line(vertex[2], vertex[3],
                                      ei[pg->imtrx]->ielem, -1);
      } else {
        surf = create_surf_facet_line(vertex[3], vertex[2],
                                      ei[pg->imtrx]->ielem, -1);
      }
      append_surf(list, surf);
    } else {
      if (edge_sign[1] > edge_sign[2]) {
        surf = create_surf_facet_line(vertex[1], vertex[2],
                                      ei[pg->imtrx]->ielem, -1);
      } else {
        surf = create_surf_facet_line(vertex[2], vertex[1],
                                      ei[pg->imtrx]->ielem, -1);
      }
      append_surf(list, surf);

      if (edge_sign[3] > edge_sign[0]) {
        surf = create_surf_facet_line(vertex[3], vertex[0],
                                      ei[pg->imtrx]->ielem, -1);
      } else {
        surf = create_surf_facet_line(vertex[0], vertex[3],
                                      ei[pg->imtrx]->ielem, -1);
      }
      append_surf(list, surf);
    }
  } break;

  default: {
    printf("number of edges = %d\n", vert_count);
    EH(-1, "Silly me, I thought this couldn't happen");
  } break;
  }
}

static int point_on_ca_boundary(int I, Exo_DB *exo)
/* check to see if point is on SS with CA condition.
 * we do this so because points on such boundaries should
 * be treated as inflection points in the zero level set
 */
{
  static int num_ss_on_ca_boundary = -1;
  static int iss_on_ca_boundary[20];
  int inflection;
  int ibc, ica, i, J;
  int iss = 0;

  /* first time here, make list of side sets with CA boundary conditions */
  if (num_ss_on_ca_boundary == -1) {
    num_ss_on_ca_boundary = 0;
    for (ibc = 0; ibc < Num_BC; ibc++) {
      if (BC_Types[ibc].BC_Name == FILL_CA_BC) {
        if ((iss = in_list(BC_Types[ibc].BC_ID, 0, exo->num_side_sets,
                           exo->ss_id)) == -1) {
          EH(-1, "Cannot locate SS index in point_on_ca_boundary");
        }
        if (num_ss_on_ca_boundary >= 20)
          EH(-1, "Need to increase array for ss_on_ca_boundary[] in "
                 "point_on_ca_boundary");
        iss_on_ca_boundary[num_ss_on_ca_boundary] = iss;
        num_ss_on_ca_boundary++;
      }
    }
  }

  inflection = FALSE; /* innocent till proven guilty */

  for (ica = 0; ica < num_ss_on_ca_boundary; ica++) {
    iss = iss_on_ca_boundary[ica];

    for (i = 0; i < Proc_SS_Node_Count[iss]; i++) {
      J = exo->ss_node_list[iss][i];
      if (J == I)
        inflection = TRUE;
    }
  }

  return inflection;
}

int find_link_intersection(double *xi, double *yi, int isovar, double isoval,
                           Integ_Elem *subelement) {
  double rlo, rhi, dr, dr_old;
  double F1, F2;
  int a, step = 0;
  double r = 0., F = 0., dF = 0.;
  double temp;
  double tol = 1.e-16;

  compute_link_level_set(xi, yi, 0., isovar, isoval, &F1, &dF, subelement);
  compute_link_level_set(xi, yi, 1., isovar, isoval, &F2, &dF, subelement);
  if (fabs(F1) < tol * fabs(F2)) {
    return (TRUE);
  }; /* that was easy */
  if (fabs(F2) < tol * fabs(F1)) {
    for (a = 0; a < pd->Num_Dim; a++) {
      xi[a] = yi[a];
    }
    return (TRUE);
  };

  if (F1 * F2 >= 0.) {
    /* no zero crossing on this link */
    return (FALSE);
  }

  rlo = 0.;
  rhi = 1.;
  dr_old = 1.;
  dr = dr_old;

  while ((fabs(dr) > tol) && step < MAX_STEP) {
    /* if 1st step or out of range or too slow, use bisection */
    if ((step == 0) || /* first step */
        ((((r - rhi) * dF - F) * ((r - rlo) * dF - F)) >=
         0.) ||                             /* out of range */
        (fabs(2. * F) > fabs(dr_old * dF))) /* too slow */
    {
      if (step == 0) {
        r = -F1 / (F2 - F1);
      } else {
        dr_old = dr;
        dr = 0.5 * (rhi - rlo);
        temp = r;
        r = rlo + dr;
        if (r - temp == 0.)
          break;
        if (r - rlo == 0.)
          break;
      }
    } else /* use Newton step */
    {
      dr_old = dr;
      dr = -F / dF;
      temp = r;
      r += dr;
      if (r - temp == 0.)
        break;
    }

    /* calc F and dF */
    compute_link_level_set(xi, yi, r, isovar, isoval, &F, &dF, subelement);

    /* keep the root bracketed between rlo and rhi */
    if (F * F1 > 0.) {
      if (r > rlo) {
        rlo = r;
      }
    } else {
      if (r < rhi) {
        rhi = r;
      }
    }

    step++;
  }

  if (step >= MAX_STEP) {
    fprintf(
        stdout,
        "The maximum iteration count was exceeded in find_link_intersection!");
    return (FALSE);
  }

  for (a = 0; a < pd->Num_Dim; a++) {
    xi[a] = (1. - r) * xi[a] + r * yi[a];
  }

  return (TRUE);
}

static void compute_link_level_set(double *xi, double *yi, double r, int isovar,
                                   double isoval, double *F, double *dF,
                                   Integ_Elem *subelement) {
  int i, a;
  double ri[3] = {0.0, 0.0, 0.0}, dri[3] = {0.0, 0.0, 0.0};
  double *esp;

  for (a = 0; a < pd->Num_Dim; a++) {
    ri[a] = (1. - r) * xi[a] + r * yi[a];
    dri[a] = yi[a] - xi[a];
  }

  if (subelement == NULL) {
    load_basis_functions(ri, bfd);

    for (i = 0, *F = -isoval, *dF = 0.0; i < ei[pg->imtrx]->dof[isovar]; i++) {

      esp = x_static +
            ei[pg->imtrx]->ieqn_ledof[ei[pg->imtrx]->lvdof_to_ledof[isovar][i]];

      *F += *esp * bf[isovar]->phi[i];

      for (a = 0; a < pd->Num_Dim; a++) {
        *dF += *esp * bf[isovar]->dphidxi[i][a] * dri[a];
      }
    }
  } else {
    /* need to account for mapping from subelement to parent element */
    double pi[3] = {0.0, 0.0, 0.0};
    int a, b, i;
    double J[DIM][DIM];
    int dim = ei[pg->imtrx]->ielem_dim;

    subelement_J(subelement, yi, J, TRUE);
    map_subelement_stu(pi, subelement, ri);

    load_basis_functions(pi, bfd);

    for (i = 0, *F = -isoval, *dF = 0.0; i < ei[pg->imtrx]->dof[isovar]; i++) {

      esp = x_static +
            ei[pg->imtrx]->ieqn_ledof[ei[pg->imtrx]->lvdof_to_ledof[isovar][i]];

      *F += *esp * bf[isovar]->phi[i];

      for (a = 0; a < dim; a++) {
        for (b = 0; b < dim; b++) {
          *dF += *esp * bf[isovar]->dphidxi[i][b] * J[a][b] * dri[a];
        }
      }
    }
  }
}

struct LS_Surf *create_surf(int type)
/* generic surface creator */
{
  struct LS_Surf *surf;

  surf = (struct LS_Surf *)smalloc(sizeof(struct LS_Surf));

  surf->type = type;
  surf->closest_point = (struct LS_Surf_Closest_Point *)smalloc(
      sizeof(struct LS_Surf_Closest_Point));
  surf->subsurf_list = NULL;
  surf->next = NULL;

  /* now for the not so generic part */
  switch (surf->type) {
  case LS_SURF_POINT: {
    surf->data = (void *)smalloc(sizeof(struct LS_Surf_Point_Data));
  } break;

  case LS_SURF_PLANE: {
    surf->data = (void *)smalloc(sizeof(struct LS_Surf_Plane_Data));
  } break;

  case LS_SURF_CIRCLE:
  case LS_SURF_SPHERE: {
    surf->data = (void *)smalloc(sizeof(struct LS_Surf_Sphere_Data));
  } break;

  case LS_SURF_FACET: {
    surf->data = (void *)smalloc(sizeof(struct LS_Surf_Facet_Data));
  } break;

  case LS_SURF_SS: {
    surf->data = (void *)smalloc(sizeof(struct LS_Surf_SS_Data));
  } break;

  case LS_SURF_NS: {
    surf->data = (void *)smalloc(sizeof(struct LS_Surf_NS_Data));
  } break;

  case LS_SURF_ISOSURFACE: {
    surf->data = (void *)smalloc(sizeof(struct LS_Surf_Iso_Data));
  } break;

  case LS_SURF_USER: {
    surf->data = (void *)smalloc(sizeof(struct LS_Surf_User_Data));
  } break;
  case LS_SURF_ARC: {
    surf->data = (void *)smalloc(sizeof(struct LS_Surf_Arc_Data));
  } break;
  default: {
    EH(-1, "create_surf called with unknown surface type");
  } break;
  }

  return surf;
}

void append_surf(struct LS_Surf_List *list, struct LS_Surf *surf)
/* append surface to list */
{
  if (list->size == 0) {
    list->start = surf;
  } else {
    list->end->next = surf;
  }
  list->end = surf;
  list->end->next = NULL;
  list->size++;
}

void append_surf_isosurf(struct LS_Surf_List *list, int isovar, double isoval) {
  struct LS_Surf *tmp_surf = NULL;
  struct LS_Surf_Iso_Data *tmp_data;

  tmp_surf = create_surf(LS_SURF_ISOSURFACE);
  tmp_data = (struct LS_Surf_Iso_Data *)tmp_surf->data;
  tmp_data->isovar = isovar;
  tmp_data->isoval = isoval;

  append_surf(list, tmp_surf);
}

static double distance(double *x, double *y) {
  double dist;

  dist = pow((x[0] - y[0]), 2.0) + pow((x[1] - y[1]), 2.0);

  if (pd->Num_Dim == 3)
    dist += pow((x[2] - y[2]), 2.0);

  return sqrt(dist);
}

int unique_surf(struct LS_Surf_List *list, struct LS_Surf *surf)
/* check if surface is unique or is a duplicate */
{
  int unique = TRUE; /* innocent till proven guilty here */

  switch (surf->type) {
  case LS_SURF_POINT: {
    struct LS_Surf_Point_Data *s = (struct LS_Surf_Point_Data *)surf->data;

    /* step through list looking for duplicates */
    list->current = list->start;
    while (list->current) {
      if (list->current->type == LS_SURF_POINT) {
        struct LS_Surf_Point_Data *sl =
            (struct LS_Surf_Point_Data *)list->current->data;
        if (distance(s->x, sl->x) < ls->Length_Scale * 0.001) {
          unique = FALSE;
          return unique;
        }
      }
      list->current = list->current->next;
    }
  } break;

  default: {
    EH(-1, "unique_surf not implemented for this surface type");
  } break;
  }

  return unique;
}

struct LS_Surf_List *create_surf_list(void)
/* surface list creator */
{
  struct LS_Surf_List *list;

  list = (struct LS_Surf_List *)smalloc(sizeof(struct LS_Surf_List));

  list->size = 0;
  list->start = list->end = list->current = NULL;

  return list;
}

struct LS_Surf *create_surf_point(double *x, int elem, double *xi,
                                  int inflection) {
  struct LS_Surf *surf;
  struct LS_Surf_Point_Data *s;

  surf = create_surf(LS_SURF_POINT);
  s = (struct LS_Surf_Point_Data *)surf->data;

  s->x[0] = x[0];
  s->x[1] = x[1];
  s->x[2] = x[2];
  s->xi[0] = xi[0];
  s->xi[1] = xi[1];
  s->xi[2] = xi[2];

  s->elem = elem;
  s->inflection = inflection;

  return surf;
}

struct LS_Surf *create_surf_facet_line(struct LS_Surf *v1, struct LS_Surf *v2,
                                       int elem, int elem_side) {
  struct LS_Surf *surf;
  struct LS_Surf_Facet_Data *s;

  surf = create_surf(LS_SURF_FACET);
  s = (struct LS_Surf_Facet_Data *)surf->data;

  s->num_points = 2;
  s->elem = elem;
  s->elem_side = elem_side;

  surf->subsurf_list = create_surf_list();

  append_surf(surf->subsurf_list, v1);
  append_surf(surf->subsurf_list, v2);

  return surf;
}

static struct LS_Surf *next_surf_or_subsurf(struct LS_Surf_List *list,
                                            int *level)
/* Recursively descend list
 * On the first call, list->current should be set to NULL.
 * This will cause the first item in the list to be returned
 * and set up for future calls.
 * NOTE: *level is returned as the level at which the surface
 * was found (how many levels of subsurfaces were descended).
 * This is incremented as each level is descended so it must
 * be initialized to zero for each call.
 */
{
  struct LS_Surf *surf;

  if (list->current == NULL) /* this is the first call for this list */
  {
    list->current = list->start;
    surf = list->current;
  } else {
    if (list->current->subsurf_list) /* working on subsurf list */
    {
      (*level)++;
      surf = next_surf_or_subsurf(list->current->subsurf_list, level);
      if (surf != NULL)
        return surf;
      /* we have reached the end of the subsurf list, returns to current list */
      (*level)--;
    }

    list->current = list->current->next;
    surf = list->current;
  }

  /* reset subsurface list so that next call will descend into it correctly */
  if ((surf != NULL) && (surf->subsurf_list != NULL)) {
    surf->subsurf_list->current = NULL;
  }

  return surf;
}

static struct LS_Surf *create_next_surf_or_subsurf(struct LS_Surf_List *list,
                                                   int level, int type)
/* create new list entry
 * level indicates how many levels of subsurfaces need to be descended before
 * adding the new surface
 */
{
  struct LS_Surf *surf;

  if (level > 0) {
    if (list->end->subsurf_list == NULL) {
      if (level > 1)
        EH(-1, "Didn't think this could happen");
      list->end->subsurf_list = create_surf_list();
    }
    surf =
        create_next_surf_or_subsurf(list->end->subsurf_list, level - 1, type);
  } else {
    surf = create_surf(type);
    append_surf(list, surf);
  }
  return surf;
}

/*
 * ls_var_initialization :  Initializes field variables that are indexed to
 * level set These would include species in two-phase transport or temperature
 * in casting.
 *
 *   Input
 *         x  = vector of nodal unknowns
 *         x_old = vector of old nodal unknowns
 */

void ls_var_initialization(double **u, Exo_DB *exo, Dpi *dpi, Comm_Ex **cx) {
  int n, k, ebi, ebj, mn;
  int block_temp, ielem, e_start, e_end;
  int *block_order;

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

  for (ebj = 0; ebj < exo->num_elem_blocks; ebj++) {
    ebi = block_order[ebj];
    mn = Matilda[ebi];
    pd = pd_glob[mn];
    mp = mp_glob[mn];

    // PRS note: hardwired here to R_PHASE1
    if (pd->e[pg->imtrx][R_LEVEL_SET]) {
      int ielem_type, num_nodes, index;
      int i, j, var, nunks, ie;
      NODAL_VARS_STRUCT *nv;

      e_start = exo->eb_ptr[ebi];
      e_end = exo->eb_ptr[ebi + 1];

      for (ielem = e_start; ielem < e_end; ielem++) {
        ielem_type = Elem_Type(exo, ielem);
        num_nodes = elem_info(NNODES, ielem_type);
        index = Proc_Connect_Ptr[ielem];
        ls->Elem_Sign = 0;

        for (n = 0; n < num_nodes; n++) /* local nodes */
        {
          scalar_value_at_local_node(ielem, ielem_type, n, LS, 0, u[pg->imtrx],
                                     exo);
          i = Proc_Elem_Connect[index++];

          for (int imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {

            nv = Nodes[i]->Nodal_Vars_Info[imtrx];

            for (j = Num_Var_Init; j < ls->Num_Var_Init + Num_Var_Init; j++) {

              var = Var_init[j].var;

              nunks = get_nv_ndofs_modMF(nv, var);

              if (nunks > 0 && pd->i[imtrx][var]) {
                if (nunks == 1) {
                  ie = Index_Solution(i, var, Var_init[j].ktype, 0, mn, imtrx);

                  level_set_property(Var_init[j].init_val_minus,
                                     Var_init[j].init_val_plus,
                                     ls->Length_Scale, &(u[imtrx][ie]), NULL);

                } else {
                  EH(-1, " Cannot initialized multiple degrees of freedom at a "
                         "node \n");
                }
              }
            }
          }
        }
      }
    }
  }

  safer_free((void **)&block_order);

  for (int imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    exchange_dof(cx[imtrx], dpi, u[imtrx], imtrx);
  }

  return;
}

static double scalar_value_at_local_node(int ielem, int ielem_type, int lnode,
                                         int var, int ktype, double *u,
                                         Exo_DB *exo) {

  double xi[DIM], val = 0;
  double scr1, scr2;

  if (u == x_static) /* be the least disruptive possible */
  {
    load_elem_dofptr(ielem, exo, x_static, x_old_static, xdot_static,
                     xdot_old_static, 0);
  } else {
    load_elem_dofptr(ielem, exo, u, u, u, u, 0);
  }

  bf_mp_init(pd);

  find_nodal_stu(lnode, ielem_type, xi, xi + 1, xi + 2);

  load_basis_functions(xi, bfd);

  switch (var) {
  case FILL:
    if (pd->v[pg->imtrx][var])
      scalar_fv_fill(esp->F, esp_dot->F, esp_old->F, bf[var]->phi,
                     ei[pg->imtrx]->dof[var], &fv->F, &scr1, &scr2);
    val = fv->F;
    break;
  case TEMPERATURE:
    if (pd->v[pg->imtrx][var])
      scalar_fv_fill(esp->T, esp_dot->T, esp_old->T, bf[var]->phi,
                     ei[pg->imtrx]->dof[var], &fv->T, &scr1, &scr2);
    val = fv->T;
    break;
  default:
    val = 0;
    break;
  }
  return (val);
}

/*
static double find_LS_value_on_side
( int ielem_type ,
        int id_side,
        double s,
        double t,
        double xi[DIM]))
*/
void iso_contour_on_side(double isoval, int dim, int ielem_type, int id_side,
                         int *ip_total, double (**s)[DIM], double **wt) {
  switch (dim) {
  case 2: {
    int local_elem_node_id[MAX_NODES_PER_SIDE];
    int nodes_per_side;
    int i;
    int have_crossing = FALSE;

    get_side_info(ielem_type, id_side, &nodes_per_side, local_elem_node_id);

    for (i = 1; i < nodes_per_side; i++) {
      int lvdof = ei[pg->imtrx]->ln_to_dof[ls->var][local_elem_node_id[i]];

      if (lvdof != -1) {

        if (sign_change(*esp->F[local_elem_node_id[0]] - isoval,
                        *esp->F[local_elem_node_id[i]] - isoval))
          have_crossing = TRUE;
      }
    }

    if (have_crossing) {
      double xi_p[DIM];

      *ip_total = 1;
      *s = (double(*)[DIM])smalloc(DIM * sizeof(double));
      *wt = (double *)smalloc(sizeof(double));
      *wt[0] = 1.0;

      if (sign_change(*esp->F[local_elem_node_id[0]] - isoval,
                      *esp->F[local_elem_node_id[1]] - isoval)) {
        load_surf_st(ielem_type, id_side, pd->Num_Dim, *s[0], -1., 0.);
        load_surf_st(ielem_type, id_side, pd->Num_Dim, xi_p, 1., 0.);

        find_link_intersection(*s[0], xi_p, ls->var, isoval, NULL);
      } else /* unusual case of having two crossings on this side, hmmm, what to
                do? */
      {
        if (nodes_per_side != 3)
          EH(-1, "Unexpected element type!");

        if (sign_change(*esp->F[local_elem_node_id[0]] - isoval,
                        *esp->F[local_elem_node_id[2]] - isoval)) {
          load_surf_st(ielem_type, id_side, pd->Num_Dim, *s[0], -1., 0.);
          load_surf_st(ielem_type, id_side, pd->Num_Dim, xi_p, 0., 0.);

          find_link_intersection(*s[0], xi_p, ls->var, isoval, NULL);
        } else {
          EH(-1, "This really shouldn't happen!");
        }
      }
    } else {
      *ip_total = 0;
      *s = NULL;
      *wt = NULL;
    }
  } break;
  case 3:
    EH(-1, "Third dimension remains problematic.\n");
  }
  return;
}
/*
static int
interface_on_side( double isoval,
                   int ielem_type,
                   int id_side )
{
  int local_elem_node_id[MAX_NODES_PER_SIDE];
  int nodes_per_side;
  int ldof_eqn;
  int i;

  int all_negative = TRUE;
  int all_positive = TRUE;



  get_side_info( ielem_type, id_side, &nodes_per_side, local_elem_node_id );

  for( i=0; ( all_positive || all_negative) &&  i<nodes_per_side ; i++)
    {
      int id;

      id = local_elem_node_id[i];

      if(  (ldof_eqn = ei[pg->imtrx]->ln_to_first_dof[ls->var][id]) != -1 )
        {
          all_negative = ( *(esp->F[ldof_eqn]) - isoval <= 0.0 ) &&
all_negative; all_positive = ( *(esp->F[ldof_eqn]) - isoval > 0.0 ) &&
all_positive;
        }
    }

  return ( !( all_negative || all_positive ) );
}
*/
/*
static double
find_LS_value_on_side( int ielem_type,
                       int id_side,
                       double s,
                       double t,
                       double xi[DIM] )
{
  double value;
  int i;

  load_surf_st( ielem_type, id_side, pd->Num_Dim, xi, s, t );

  load_basis_functions( xi, bfd );

  for(i=0, value=0.0; i<ei[pg->imtrx]->dof[ls->var]; i++ )
    {
      value += *(esp->F[i])*bf[ls->var]->phi[i];
    }
  return(value);
}
*/
/******************************************************************************
 * level_set_interface(): Calculates useful level set functions and their
 *                        derivatives.  Some terms are probably not worth
 *        computing (such as d_H_dgradF[]) but it's easier to be systematic.
 *        The routine tries to be as efficient as possible -- no derivatives
 *        are calculated if we're not near the interface or if do_deriv is
 *FALSE.
 *
 * Input
 * -----
 *   F	      = The value of the level set (distance) function.
 *   grad_F   = Gradient of F.
 *   width    = Width of the interface region; width = 2*alpha.
 *   do_deriv = If TRUE, calculate the derivatives.
 *
 * Output
 * ------
 *   near	     = TRUE if the current posistion (i.e. F) is in the
 *                     interfacial region
 *   H		     = Heaviside function; smooth form of: H=0 for F<0,
 *                     H=1 for F>0
 *   d_H_dF	     = Derivative of H w.r.t. F. Normally this would be
 *                     delta but not here
 *   d_H_dgradF      = Derivative of H w.r.t. the components of grad(F).
 *   delta	     = Delta function: smooth form of: delta(F) = 1 for
 *                     F=0, =0 for F != 0; N.B. This delta has a correction
 *                     for cases where F is not a pure distance function.
 *   d_delta_dF	     = Derivative of delta w.r.t. F
 *   d_delta_dgradF  = Derivative of delta w.r.t. the components of grad(F).
 *   normal	     = Normal vector: typically normal = grad(F); here we use
 *                     normal = grad(F) / |grad(F)| to be safe.
 *   d_normal_dF     = Derivative of normal w.r.t. F.
 *   d_noraml_dgradF = Derivative of normal w.r.t. the components of grad(F).
 *
 * Return
 * ------
 *   0 = Success.
 *  -1 = Failure.
 *
 * Example Usage of Derivative Terms
 * ---------------------------------
 *
 * d_normal_dgradF[a][b] = d normal[a] / d gradF[b]
 *
 * So, to compute Jacobian terms like d normal[a]/d F_j:
 *
 * d normal[a] / d F_j = d_normal_dF[a]        * bf[FILL]->phi[j]
 *                     + d_normal_dgradF[a][0] * bf[FILL]->grad_phi[j][0]
 *                     + d_normal_dgradF[a][1] * bf[FILL]->grad_phi[j][1]
 *                     + d_normal_dgradF[a][2] * bf[FILL]->grad_phi[j][2]
 *
 ******************************************************************************/
int level_set_interface(const double F, const double grad_F[DIM],
                        const double width, const int do_deriv, int *near,
                        double *H, double *d_H_dF, double d_H_dgradF[DIM],
                        double *delta, double *d_delta_dF,
                        double d_delta_dgradF[DIM], double normal[DIM],
                        double d_normal_dF[DIM],
                        double d_normal_dgradF[DIM][DIM]) {
  int a, b;
  int status = 0;
  double mag, maginv;
  double alpha = 0.5 * width;

  /**************************************************
   * Initialize the results.
   **************************************************/
  /* if do_deriv == TRUE, we assume these pointers aren't NULL. */
  if (do_deriv) {
    *d_H_dF = 0.;
    *d_delta_dF = 0.;
    memset(d_H_dgradF, 0, sizeof(double) * DIM);
    memset(d_delta_dgradF, 0, sizeof(double) * DIM);
    memset(d_normal_dF, 0, sizeof(double) * DIM);
    memset(d_normal_dgradF, 0, sizeof(double) * DIM * DIM);
  }

  /* If we're not in the mushy zone, set the essentials and bail out. */
  if (fabs(F) >= alpha) {
    *near = FALSE;
    *H = (F < 0) ? 0. : 1.;
    *delta = 0.;
    memset(normal, 0, sizeof(double) * DIM);
    return (status);
  }

  /**************************************************
   * Calculate the functions for real
   **************************************************/
  mag = 0.0;
  for (a = 0; a < VIM; a++) {
    normal[a] = grad_F[a];
    mag += grad_F[a] * grad_F[a];
  }
  mag = sqrt(mag);
  maginv = (mag == 0.0) ? 1.0 : 1.0 / mag;

  for (a = 0; a < VIM; a++) {
    normal[a] *= maginv;
  }

  *near = TRUE;
  *H = 0.5 * (1. + F / alpha + sin(M_PIE * F / alpha) / M_PIE);
  *delta = 0.5 * (1. + cos(M_PIE * F / alpha)) / alpha * mag;

  /* Bail out if we don't need to calculate the derivatives. */
  if (!do_deriv)
    return (status);

  /**************************************************
   * Calculate the derivatives of the functions
   **************************************************/

  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      d_normal_dgradF[a][b] = -pow(maginv, 3.0) * grad_F[b] * grad_F[a];
    }
    d_normal_dgradF[a][a] += maginv;
  }

  *d_H_dF = 0.5 * (1. + cos(M_PIE * F / alpha)) / alpha;
  *d_delta_dF = -0.5 * sin(M_PIE * F / alpha) * M_PIE / (alpha * alpha) * mag;
  for (a = 0; a < VIM; a++) {
    d_delta_dgradF[a] = *d_H_dF * maginv * grad_F[a];
  }

  return (status);
}

/******************************************************************************
 *
 * level_set_property() : Calculate a general, scalar material property using
 *                        the level set function.
 *
 * Input
 * -----
 *   p0	   = Material property for FILL < 0 ("minus" side)
 *   p1	   = Material property for FILL > 0 ("plus" side)
 *
 * Output
 * ------
 *   pp		  = Material property at the current coordinate.
 *   d_pp_dF[MDE] = Derivative of the material property w.r.t. the FILL
 *                  variable. N.B. If d_pp_dF == NULL, this derivative
 *                  is not calculated.
 *
 * Returns
 * -------
 *   0 = Success.
 *  -1 = Failure.
 *
 ******************************************************************************/
int level_set_property(const double p0, const double p1, const double width,
                       double *pp, double d_pp_dF[MDE]) {
  int var, j;
  int do_deriv;

  /* See if we need to bother with derivatives. */
  do_deriv = d_pp_dF != NULL;

  /* Fetch the level set interfacial functions. */
  load_lsi(width);

  /* Calculate the material property. */
  if (ls->Elem_Sign == -1)
    *pp = p0;
  else if (ls->Elem_Sign == 1)
    *pp = p1;
  else
    *pp = p0 + (p1 - p0) * lsi->H;

  if (ls->Elem_Sign != 0 && do_deriv) {
    var = ls->var;
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      d_pp_dF[j] = 0.;
    }
  }

  /* Bail out if we don't need derivatives or if we're not in the mushy zone. */
  if (!do_deriv || !lsi->near || ls->Elem_Sign != 0)
    return (0);

#ifdef COUPLED_FILL
  load_lsi_derivs();

  /* Calculate the deriviatives of the material property w.r.t. FILL. */
  var = ls->var;
  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
    /* Calculate the Jacobian terms. */
    d_pp_dF[j] = (p1 - p0) * lsi->d_H_dF[j];
  }

#endif /* COUPLED_FILL */

#if 0
  /*
   * NB: H(F) doesn't not (currently) depend on grad_F and hence pp
   *     does not depend on MESH_DISPLACEMENTs.  I'll leave this here 
   *     just in case....  Hence the 'if 0'
   */

  if ( ! pd->v[pg->imtrx][MESH_DISPLACEMENT1] ) return(0);

  for ( b=0; b < VIM; b++ )
    {
      var = MESH_DISPLACEMENT1 + b;
      for ( j=0; j < ei[pg->imtrx]->dof[var]; j++)
	{
	  d_pp_dmesh[b][j] = (p1 - p0) * lsi->d_H_dmesh[b][j];
	}
    }
#endif /* 0 */

  return (0);
}

int level_set_property_offset(const double p0, const double p1,
                              const double width, double *pp,
                              double d_pp_dF[MDE]) {
  int var, j;
  int do_deriv;

  /* See if we need to bother with derivatives. */
  do_deriv = d_pp_dF != NULL;

  /* Fetch the level set interfacial functions. */
  load_lsi_offset(width);

  /* Calculate the material property. */
  if (ls->Elem_Sign == -1)
    *pp = p0;
  else if (ls->Elem_Sign == 1)
    *pp = p1;
  else
    *pp = p0 + (p1 - p0) * lsi->H;

  if (ls->Elem_Sign != 0 && do_deriv) {
    var = ls->var;
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      d_pp_dF[j] = 0.;
    }
  }

  /* Bail out if we don't need derivatives or if we're not in the mushy zone. */
  if (!do_deriv || !lsi->near || ls->Elem_Sign != 0)
    return (0);

#ifdef COUPLED_FILL
  load_lsi_derivs();

  /* Calculate the deriviatives of the material property w.r.t. FILL. */
  var = ls->var;
  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
    /* Calculate the Jacobian terms. */
    d_pp_dF[j] = (p1 - p0) * lsi->d_H_dF[j];
  }

#endif /* COUPLED_FILL */

#if 0
  /*
   * NB: H(F) doesn't not (currently) depend on grad_F and hence pp
   *     does not depend on MESH_DISPLACEMENTs.  I'll leave this here
   *     just in case....  Hence the 'if 0'
   */

  if ( ! pd->v[pg->imtrx][MESH_DISPLACEMENT1] ) return(0);

  for ( b=0; b < VIM; b++ )
    {
      var = MESH_DISPLACEMENT1 + b;
      for ( j=0; j < ei[pg->imtrx]->dof[var]; j++)
        {
          d_pp_dmesh[b][j] = (p1 - p0) * lsi->d_H_dmesh[b][j];
        }
    }
#endif /* 0 */

  return (0);
}

/******************************************************************************
 *
 * ls_transport_property() : Almost identical in function to level_set_property.
 *                           Computes "diffuse" level set Hv function of
 material
 *                           property.  Differs from level_set_property in that
 it
 *                           returns sensitivity wrt F field instead of F nodal
 unknowns.
 *                           It is used primarily for diffusivities, thermal
 conductivities
 *
 *
 * Input
 * -----
 *   p0	   = Material property for FILL < 0 ("minus" side)
 *   p1	   = Material property for FILL > 0 ("plus" side)
 *
 * Output
 * ------
 *   pp		  = Material property at the current coordinate.
 *   d_pp_dF      = Derivative of the material property w.r.t. the FILL
 *                  variable. N.B. If d_pp_dF == NULL, this derivative
 *                  is not calculated.
 *
 * Returns
 * -------mm_fill_ls.c:   else if ( ls->Elem_Sign == 1 ) *pp = p1;

 *   0 = Success.
 *  -1 = Failure.
 *
 ******************************************************************************/
int ls_transport_property(const double p0, const double p1, const double width,
                          double *pp, double *d_pp_dF) {

  int do_deriv;

  /* See if we need to bother with derivatives. */
  do_deriv = d_pp_dF != NULL;

  /* Fetch the level set interfacial functions. */
  load_lsi(width);

  /* Calculate the material property. */
  if (ls->Elem_Sign == -1)
    *pp = p0;
  else if (ls->Elem_Sign == 1)
    *pp = p1;
  else
    *pp = p0 + (p1 - p0) * lsi->H;

  if (do_deriv)
    *d_pp_dF = 0.;

  /* Bail out if we don't need derivatives or if we're not in the mushy zone. */
  if (!do_deriv || !lsi->near || ls->Elem_Sign != 0)
    return (0);

#ifdef COUPLED_FILL
  load_lsi_derivs();

  /* Calculate the deriviatives of the material property w.r.t. FILL. */
  if (ls->Elem_Sign == 0)
    *d_pp_dF = (p1 - p0) * lsi->dH;
  else
    *d_pp_dF = 0.;

#endif /* COUPLED_FILL */

  return (0);
}

void determine_ls_elem_overlap_state(void) {
  /* For level set problems we need to know the characteristics of all of this
   * element */
  /* ls->elem_overlap_state = 0 -> no crossing in element
     ls->elem_overlap_state = 1 -> crossing in element
   */

  if (ls == NULL)
    return;

  ls->elem_overlap_state = FALSE;

  if (fabs(ls->Length_Scale) > 1.0E-20) /* diffuse interface */
  {
    if (current_elem_overlaps_interface(ls->Length_Scale))
      ls->elem_overlap_state = TRUE;
  } else /* sharp interface */
  {
    if (current_elem_on_isosurface(ls->var, 0.))
      ls->elem_overlap_state = TRUE;
  }
}

void load_xfem_for_elem(double x[], const Exo_DB *exo) {
  int eqn;

  /* only evaluate element quantities if new element */
  if (ei[pg->imtrx]->ielem != xfem->ielem || Debug_Flag < 0) {
    xfem->ielem = ei[pg->imtrx]->ielem;
    xfem->elem_state = current_elem_xfem_state(xfem->node_var_state,
                                               &xfem->elem_var_state, x, exo);
    /* fill in bogus values to avoid pathological cases */
    xfem->xi[0] = -2;
    xfem->xi[1] = -2;
    xfem->xi[2] = -2;

    /* check if gradients will be needed for XG-type interpolations */
    xfem->have_XG = FALSE;
    for (eqn = V_FIRST; eqn < V_LAST && !xfem->have_XG; eqn++) {
      if (pd->i[pg->imtrx][eqn] == I_P1_XG ||
          pd->i[pg->imtrx][eqn] == I_Q1_XG ||
          pd->i[pg->imtrx][eqn] == I_Q2_XG ||
          pd->i[pg->imtrx][eqn] == I_Q1_HG ||
          pd->i[pg->imtrx][eqn] == I_Q1_HVG ||
          pd->i[pg->imtrx][eqn] == I_Q2_HG || pd->i[pg->imtrx][eqn] == I_Q2_HVG)
        xfem->have_XG = TRUE;
    }

    /* check if quantities will be needed for discontinuous-type interpolations
     */
    xfem->have_disc = FALSE;
    for (eqn = V_FIRST; eqn < V_LAST && !xfem->have_disc; eqn++) {
      if (pd->i[pg->imtrx][eqn] == I_Q1_HV ||
          pd->i[pg->imtrx][eqn] == I_Q1_HG ||
          pd->i[pg->imtrx][eqn] == I_Q1_HVG ||
          pd->i[pg->imtrx][eqn] == I_Q2_HV ||
          pd->i[pg->imtrx][eqn] == I_Q2_HG || pd->i[pg->imtrx][eqn] == I_Q2_HVG)
        xfem->have_disc = TRUE;
    }
#if 0
      /* possibly dangerous, but try to speed things up by setting
         active_interp_ledof for XFEM vars 
       */
      for (eqn = V_FIRST; eqn < V_LAST; eqn++)
        {
          if ( pd->e[pg->imtrx][eqn] )
	    {
	      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) 
	        {
	          ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];
		  
                  xfem_dof_state( i, pd->i[pg->imtrx][eqn], ei[pg->imtrx]->ielem_shape,
                                  &xfem_active, &extended_dof, &base_interp, &base_dof );
		  if ( extended_dof && !xfem_active ) ei[pg->imtrx]->active_interp_ledof[ledof] = FALSE;
		}
	    }
	}
#endif
  }
}

void load_xfem_for_stu(const double xi[]) {
  int i, F_elem_type = -1;
  int dof_ls;

  /* in XFEM, basis functions depend on distance function */
  /* here, we evaluate these functions */

  /* make sure xfem has been set up for element */
  if (ei[pg->imtrx]->ielem != xfem->ielem) {
    EH(-1, "Must call load_xfem_for_elem before calling load_xfem_for_stu.\n");
  }

  if (xfem->elem_state == 1) {
    /* only evaluate distance if at new point */
    if (xi[0] != xfem->xi[0] || xi[1] != xfem->xi[1] || xi[2] != xfem->xi[2]) {
      double alpha = 0.5 * ls->Length_Scale;

      xfem->xi[0] = xi[0];
      xfem->xi[1] = xi[1];
      xfem->xi[2] = xi[2];

      dof_ls = ei[pg->imtrx]->dof[ls->var];

      if (pd->i[pg->imtrx][ls->var] == I_Q1) {
        F_elem_type = pd->Num_Dim == 2 ? BILINEAR_QUAD : TRILINEAR_HEX;
      } else if (pd->i[pg->imtrx][FILL] == I_Q2) {
        F_elem_type = pd->Num_Dim == 2 ? BIQUAD_QUAD : TRIQUAD_HEX;
      } else if (pd->i[pg->imtrx][FILL] == I_NOTHING) {
        xfem->near = FALSE;
        xfem->F = 2. * ls->Length_Scale;
        xfem->H = 1.;
        xfem->delta = 0.;
        xfem->dF_xi[0] = 0.;
        xfem->dF_xi[1] = 0.;
        xfem->dF_xi[2] = 0.;
        xfem->bf_plus = 0.;
        xfem->F_plus = 0.;
        xfem->grad_bf_plus[0] = 0.;
        xfem->grad_bf_plus[1] = 0.;
        xfem->grad_bf_plus[2] = 0.;
        xfem->grad_F_plus[0] = 0.;
        xfem->grad_F_plus[1] = 0.;
        xfem->grad_F_plus[2] = 0.;
        return;
      } else
        EH(-1, "Unexpected interpolation for FILL");

      xfem->F = 0.;
      for (i = 0; i < dof_ls; i++) {
        xfem->F += shape(xi[0], xi[1], xi[2], F_elem_type, PSI, i) * *esp->F[i];
      }

      if (fabs(xfem->F) < alpha) {
        xfem->near = TRUE;
        xfem->H =
            0.5 * (1. + xfem->F / alpha + sin(M_PIE * xfem->F / alpha) / M_PIE);
        xfem->delta = 0.5 * (1. + cos(M_PIE * xfem->F / alpha)) / alpha;
      } else {
        xfem->near = FALSE;
        xfem->H = (xfem->F < 0.0) ? 0.0 : 1.0;
        xfem->delta = 0.;
        /* DRN: not sure if we want to do this or not */
        /*if( fabs(xfem->F) < DBL_SMALL) xfem->H = 0.5;*/
        if (ls->on_sharp_surf) {
          xfem->H = (ls->Elem_Sign < 0) ? 0.0 : 1.0;
          xfem->delta = 0.;
        }
      }

      /* gradients and other derived quantities */
      xfem->dF_xi[0] = 0.;
      xfem->dF_xi[1] = 0.;
      xfem->dF_xi[2] = 0.;
      xfem->bf_plus = 0.;
      xfem->F_plus = 0.;
      xfem->grad_bf_plus[0] = 0.;
      xfem->grad_bf_plus[1] = 0.;
      xfem->grad_bf_plus[2] = 0.;
      xfem->grad_F_plus[0] = 0.;
      xfem->grad_F_plus[1] = 0.;
      xfem->grad_F_plus[2] = 0.;

      /* drn - this contains somewhat convoluted logic to desparately try to
         keep the costs low but keep from having too many if tests

         basically, we need the following
         xfem->near:      need xfem->dF_xi[]
         xfem->have_XG:   need xfem->dF_xi[] & xfem->F_plus &
         xfem->grad_F_plus[] xfem->have_disc: need xfem->F_plus &
         xfem->grad_F_plus[] & xfem->bf_plus & xfem->grad_bf_plus[]

         the most important thing is to avoid extra calls to shape.  the logic
         below sometimes calculates unnecessary quantities, but should avoid
         unnecessary calls to shape.
       */
      if (xfem->have_XG || xfem->have_disc || xfem->near) {
        double phi, dphi = 0.0;

        if (xfem->have_XG || xfem->have_disc) {
          for (i = 0; i < dof_ls; i++) {
            if (*esp->F[i] >= 0.) {
              phi = shape(xi[0], xi[1], xi[2], F_elem_type, PSI, i);
              xfem->bf_plus += phi;
              xfem->F_plus += phi * *esp->F[i];
            }
          }
        }

        /* 1-D, 2-D, or 3-D (I guess that's everybody!) */
        for (i = 0; i < dof_ls; i++) {
          if (xfem->near || xfem->have_XG || *esp->F[i] >= 0.) {
            dphi = shape(xi[0], xi[1], xi[2], F_elem_type, DPSI_S, i);
            xfem->dF_xi[0] += dphi * *esp->F[i];
          }
          if ((xfem->have_XG || xfem->have_disc) && *esp->F[i] >= 0.) {
            xfem->grad_bf_plus[0] += dphi;
            xfem->grad_F_plus[0] += dphi * *esp->F[i];
          }
        }

        /* 2-D or 3-D */
        if (ei[pg->imtrx]->ielem_dim > 1) {
          for (i = 0; i < dof_ls; i++) {
            if (xfem->near || xfem->have_XG || *esp->F[i] >= 0.) {
              dphi = shape(xi[0], xi[1], xi[2], F_elem_type, DPSI_T, i);
              xfem->dF_xi[1] += dphi * *esp->F[i];
            }
            if ((xfem->have_XG || xfem->have_disc) && *esp->F[i] >= 0.) {
              xfem->grad_bf_plus[1] += dphi;
              xfem->grad_F_plus[1] += dphi * *esp->F[i];
            }
          }
        }

        /* 3-D only */
        if (ei[pg->imtrx]->ielem_dim > 2) {
          for (i = 0; i < dof_ls; i++) {
            if (xfem->near || xfem->have_XG || *esp->F[i] >= 0.) {
              dphi = shape(xi[0], xi[1], xi[2], F_elem_type, DPSI_U, i);
              xfem->dF_xi[2] += dphi * *esp->F[i];
            }
            if ((xfem->have_XG || xfem->have_disc) && *esp->F[i] >= 0.) {
              xfem->grad_bf_plus[2] += dphi;
              xfem->grad_F_plus[2] += dphi * *esp->F[i];
            }
          }
        }
      }
    }
  }
}

void xfem_correct(int num_total_nodes, double x[], double xdot[],
                  double x_old[], double xdot_old[], double delta_x[],
                  double theta_arg, double delta_t) {
  NODE_INFO_STRUCT *node;
  NODAL_VARS_STRUCT *nv;
  VARIABLE_DESCRIPTION_STRUCT *vd;
  int I, ie, idof, lvdesc, var_type;
  int interp;
  int ioffset;

  for (I = 0; I < num_total_nodes; I++) {
    node = Nodes[I];
    nv = node->Nodal_Vars_Info[pg->imtrx];
    ioffset = 0;
    for (lvdesc = 0; lvdesc < nv->Num_Var_Desc; lvdesc++) {
      vd = nv->Var_Desc_List[lvdesc];
      if (MatID == -1)
        MatID = 0;
      var_type = vd->Variable_Type;

      interp = pd_glob[MatID]->i[pg->imtrx][var_type];

      int fill_matrix = pd_glob[MatID]->mi[R_FILL];
      if (upd->Total_Num_Matrices == 1 && fill_matrix < 0) {
        fill_matrix = pg->imtrx;
      } else if (fill_matrix < 0) {
        EH(-1, "Could not find fill matrix");
      }
      if (is_xfem_interp(interp)) {
        double F, F_old, F_prev;

        if (delta_x == NULL) {

          gnn_distance(I, pg->matrices[fill_matrix].x,
                       pg->matrices[fill_matrix].x_old, NULL, &F, &F_old, NULL);

        } else {

          gnn_distance(I, pg->matrices[fill_matrix].x,
                       pg->matrices[fill_matrix].x_old, delta_x, &F, &F_old,
                       &F_prev);

          if (sign_change(F, F_prev)) /* corrector changed sign of node */
          {
            int ie_xfem;
            switch (interp) {
            case I_P0_GP:
            case I_Q1_GP:
            case I_Q2_GP:
            case I_P0_GN:
            case I_Q1_GN:
            case I_Q2_GN:
/* Debugging */
#if 0
                          DPRINTF(stderr, "Debugging: node changed sign in corrector:\n");
			  DPRINTF(stderr, "F = %g, F_prev = %g at node, I = %d, v=%d\n", F, F_prev, I, var_type);
                          DPRINTF(stderr, "x[%d] = %g\n\n", ie, x[ie]);
#endif
              break;
            case I_P0_G:
            case I_Q1_G:
            case I_Q2_G: {
              double temp;
              ie = node->First_Unknown[pg->imtrx] + ioffset;
              ie_xfem = node->First_Unknown[pg->imtrx] + ioffset + 1;

              temp = x[ie];
              x[ie] = x[ie_xfem];
              x[ie_xfem] = temp;
/* Debugging */
#if 0
                          DPRINTF(stderr, "Debugging: node changed sign in corrector:\n");
			  DPRINTF(stderr, "F = %g, F_prev = %g at node, I = %d, v=%d\n", F, F_prev, I, var_type);
                          DPRINTF(stderr, "x[%d] = %g, x[%d+1] = %g\n\n", ie, x[ie], ie, x[ie+1]);
#endif
            } break;
            case I_P0_XV:
            case I_Q1_XV:
            case I_Q2_XV: {
              int sign_prev = -1;
              if (F_prev >= 0.)
                sign_prev = 1;

              ie = node->First_Unknown[pg->imtrx] + ioffset;
              ie_xfem = node->First_Unknown[pg->imtrx] + ioffset + 1;

              x[ie] = x[ie] - sign_prev * x[ie_xfem];
/* Debugging */
#if 0
                          DPRINTF(stderr, "Debugging: node changed sign in corrector:\n");
			  DPRINTF(stderr, "F = %g, F_prev = %g at node, I = %d, v=%d\n", F, F_prev, I, var_type);
                          DPRINTF(stderr, "x[%d] = %g, x[%d+1] = %g\n\n", ie, x[ie], ie, x[ie+1]);
#endif
            } break;
            default:
              WH(-1, "x correction not yet implemented for this type of XFEM");
            }
          }
        }

        if (sign_change(F,
                        F_old)) /* node changed sign from previous time step */
        {
          int ie_xfem;
          switch (interp) {
          case I_P0_GP:
          case I_Q1_GP:
          case I_Q2_GP:
          case I_P0_GN:
          case I_Q1_GN:
          case I_Q2_GN:
/* Debugging */
#if 0
                          DPRINTF(stderr, "Debugging: node has changed sign since previous time step:\n");
			  DPRINTF(stderr, "F = %g, F_old = %g at node, I = %d, v=%d\n", F, F_old, I, var_type);
			  DPRINTF(stderr, "x[%d] = %g\n", ie, x[ie]); 
			  DPRINTF(stderr, "x_old[%d] = %g\n", ie, x_old[ie]); 
                          DPRINTF(stderr, "xdot[%d] = %g\n\n", ie, xdot[ie]);
#endif
            break;
          case I_P0_G:
          case I_Q1_G:
          case I_Q2_G:
            ie = node->First_Unknown[pg->imtrx] + ioffset;
            ie_xfem = node->First_Unknown[pg->imtrx] + ioffset + 1;

            xdot[ie] =
                (1.0 + 2.0 * theta_arg) / delta_t * (x[ie] - x_old[ie_xfem]) -
                (2.0 * theta_arg) * xdot_old[ie_xfem];
            xdot[ie_xfem] =
                (1.0 + 2.0 * theta_arg) / delta_t * (x[ie_xfem] - x_old[ie]) -
                (2.0 * theta_arg) * xdot_old[ie];
            break;
          case I_P0_XV:
          case I_Q1_XV:
          case I_Q2_XV: {
            int sign_old = -1;
            double x_old_i, xdot_old_i;
            if (F_old >= 0.)
              sign_old = 1;

            ie = node->First_Unknown[pg->imtrx] + ioffset;
            ie_xfem = node->First_Unknown[pg->imtrx] + ioffset + 1;

            x_old_i = x_old[ie] - sign_old * x_old[ie_xfem];
            xdot_old_i = xdot_old[ie] - sign_old * xdot_old[ie_xfem];

            xdot[ie] = (1.0 + 2.0 * theta_arg) / delta_t * (x[ie] - x_old_i) -
                       (2.0 * theta_arg) * xdot_old_i;

            xdot[ie_xfem] = (1.0 + 2.0 * theta_arg) / delta_t *
                                (x[ie_xfem] - x_old[ie_xfem]) -
                            (2.0 * theta_arg) * xdot_old[ie_xfem];
/* Debugging */
#if 0
                          DPRINTF(stderr, "Debugging: node has changed sign since previous time step:\n");
			  DPRINTF(stderr, "F = %g, F_old = %g at node, I = %d, v=%d\n", F, F_old, I, var_type);
			  DPRINTF(stderr, "x[%d] = %g, x[%d+1] = %g\n", ie, x[ie], ie, x[ie+1]); 
			  DPRINTF(stderr, "x_old[%d] = %g, x_old[%d+1] = %g\n", ie, x_old[ie], ie, x_old[ie+1]); 
                          DPRINTF(stderr, "xdot[%d] = %g, xdot[%d+1] = %g\n\n", ie, xdot[ie], ie, xdot[ie+1]);
#endif
          } break;
          case I_P1_XG:
          case I_Q1_XG:
          case I_Q2_XG: {

            double delta_g = fabs(F_old) - fabs(F);

            ie = node->First_Unknown[pg->imtrx] + ioffset;
            ie_xfem = node->First_Unknown[pg->imtrx] + ioffset + 1;

            xdot[ie] = (1.0 + 2.0 * theta_arg) / delta_t *
                           (x[ie] - x_old[ie] - delta_g) -
                       (2.0 * theta_arg) * xdot_old[ie];

            xdot[ie_xfem] = (1.0 + 2.0 * theta_arg) / delta_t *
                                (x[ie_xfem] - x_old[ie_xfem]) -
                            (2.0 * theta_arg) * xdot_old[ie_xfem];
          } break;
          default:
            WH(-1, "xdot correction not yet implemented for this type of XFEM");
          }
        } else {
          for (idof = 0; idof < vd->Ndof; idof++) {
            ie = node->First_Unknown[pg->imtrx] + ioffset + idof;
            xdot[ie] = (1.0 + 2.0 * theta_arg) / delta_t * (x[ie] - x_old[ie]) -
                       (2.0 * theta_arg) * xdot_old[ie];
          }
        }
      }

      ioffset += vd->Ndof;
    }
  }
}

void xfem_predict(int num_total_nodes, int numProcUnknowns, double delta_t,
                  double delta_t_old, double delta_t_older, double theta_arg,
                  double x[], double x_old[], double x_older[],
                  double x_oldest[], double xdot[], double xdot_old[],
                  double xdot_older[]) {
  NODE_INFO_STRUCT *node;
  NODAL_VARS_STRUCT *nv;
  VARIABLE_DESCRIPTION_STRUCT *vd;
  int I, ie, lvdesc, var_type;
  int interp;
  double c1, c2, c3 = 0.0;
  int ioffset;

  if (theta_arg == 0.5) {
    c1 = delta_t * (delta_t + delta_t_old) / delta_t_older /
         (delta_t_older + delta_t_old);
    c2 = -delta_t * (delta_t + delta_t_old + delta_t_older) /
         (delta_t_old * delta_t_older);
    c3 = (delta_t + delta_t_old + delta_t_older) * (delta_t + delta_t_old) /
         delta_t_old / (delta_t_older + delta_t_old);
  } else {
    c1 = delta_t * (1.0 + theta_arg * delta_t / delta_t_old);
    c2 = theta_arg * (delta_t * delta_t) / (delta_t_old);
  }

  for (I = 0; I < num_total_nodes; I++) {
    node = Nodes[I];
    nv = node->Nodal_Vars_Info[pg->imtrx];
    ioffset = 0;
    for (lvdesc = 0; lvdesc < nv->Num_Var_Desc; lvdesc++) {
      vd = nv->Var_Desc_List[lvdesc];
      if (MatID == -1)
        MatID = 0;
      var_type = vd->Variable_Type;

      interp = pd_glob[MatID]->i[pg->imtrx][var_type];

      if (is_xfem_interp(interp)) {
        double F, F_old, F_older, F_oldest;

        gnn_distance(I, x, x_old, NULL, &F, &F_old, NULL);
        gnn_distance(I, x_older, x_oldest, NULL, &F_older, &F_oldest, NULL);

        if (sign_change(F, F_old) || sign_change(F, F_older) ||
            (theta_arg == 0.5 && sign_change(F, F_oldest))) {
          int ie_xfem;
          switch (interp) {
          case I_P0_GP:
          case I_Q1_GP:
          case I_Q2_GP:
          case I_P0_GN:
          case I_Q1_GN:
          case I_Q2_GN:
            break;
          case I_P0_G:
          case I_Q1_G:
          case I_Q2_G: {
            ie = node->First_Unknown[pg->imtrx] + ioffset;
            ie_xfem = node->First_Unknown[pg->imtrx] + ioffset + 1;

            if (theta_arg == 0.5) {
              double x_old_local, x_old_other;
              double x_older_local, x_older_other;
              double x_oldest_local, x_oldest_other;

              if (sign_change(F, F_old)) {
                x_old_local = x_old[ie_xfem];
                x_old_other = x_old[ie];
              } else {
                x_old_local = x_old[ie];
                x_old_other = x_old[ie_xfem];
              }
              if (sign_change(F, F_older)) {
                x_older_local = x_older[ie_xfem];
                x_older_other = x_older[ie];
              } else {
                x_older_local = x_older[ie];
                x_older_other = x_older[ie_xfem];
              }
              if (sign_change(F, F_oldest)) {
                x_oldest_local = x_oldest[ie_xfem];
                x_oldest_other = x_oldest[ie];
              } else {
                x_oldest_local = x_oldest[ie];
                x_oldest_other = x_oldest[ie_xfem];
              }

              x[ie] =
                  c3 * x_old_local + c2 * x_older_local + c1 * x_oldest_local;
              x[ie_xfem] =
                  c3 * x_old_other + c2 * x_older_other + c1 * x_oldest_other;
            } else {
              double x_old_local, x_old_other;
              double xdot_old_local, xdot_old_other;
              double xdot_older_local, xdot_older_other;

              if (sign_change(F, F_old)) {
                x_old_local = x_old[ie_xfem];
                x_old_other = x_old[ie];
                xdot_old_local = xdot_old[ie_xfem];
                xdot_old_other = xdot_old[ie];
              } else {
                x_old_local = x_old[ie];
                x_old_other = x_old[ie_xfem];
                xdot_old_local = xdot_old[ie];
                xdot_old_other = xdot_old[ie_xfem];
              }
              if (sign_change(F, F_older)) {
                xdot_older_local = xdot_older[ie_xfem];
                xdot_older_other = xdot_older[ie];
              } else {
                xdot_older_local = xdot_older[ie];
                xdot_older_other = xdot_older[ie_xfem];
              }

              x[ie] = x_old_local + c1 * xdot_old_local - c2 * xdot_older_local;
              x[ie_xfem] =
                  x_old_other + c1 * xdot_old_other - c2 * xdot_older_other;
            }
          } break;
          case I_P0_XV:
          case I_Q1_XV:
          case I_Q2_XV: {
            int sign_old = 0, sign_older = 0, sign_oldest = 0;

            if (sign_change(F, F_old)) {
              if (F_old < 0.)
                sign_old = -1;
              else
                sign_old = 1;
            }
            if (sign_change(F, F_older)) {
              if (F_older < 0.)
                sign_older = -1;
              else
                sign_older = 1;
            }
            if (sign_change(F, F_oldest)) {
              if (F_oldest < 0.)
                sign_oldest = -1;
              else
                sign_oldest = 1;
            }

            ie = node->First_Unknown[pg->imtrx] + ioffset;
            ie_xfem = node->First_Unknown[pg->imtrx] + ioffset + 1;

            if (theta_arg == 0.5) {
              x[ie] = c3 * (x_old[ie] - sign_old * x_old[ie_xfem]) +
                      c2 * (x_older[ie] - sign_older * x_older[ie_xfem]) +
                      c1 * (x_oldest[ie] - sign_oldest * x_oldest[ie_xfem]);

              x[ie_xfem] = c3 * x_old[ie_xfem] + c2 * x_older[ie_xfem] +
                           c1 * x_oldest[ie_xfem];
            } else {
              x[ie] = (x_old[ie] - sign_old * x_old[ie_xfem]) +
                      c1 * (xdot_old[ie] - sign_old * xdot_old[ie_xfem]) -
                      c2 * (xdot_older[ie] - sign_older * xdot_older[ie_xfem]);
              x[ie_xfem] = (x_old[ie_xfem] + c1 * xdot_old[ie_xfem] -
                            c2 * xdot_older[ie_xfem]);
            }
/* Debugging */
#if 0
                        DPRINTF(stderr, "Debugging: node changed sign in predictor:\n");
			DPRINTF(stderr, "F = %g, F_old = %g at node, I = %d, v=%d\n", F, F_old, I, var_type);
			DPRINTF(stderr, "x_old[%d] = %g, x_old[%d+1] = %g\n", ie, x_old[ie], ie, x_old[ie+1]); 
			DPRINTF(stderr, "xdot_old[%d] = %g, xdot_old[%d+1] = %g\n", ie, xdot_old[ie], ie, xdot_old[ie+1]); 
                        DPRINTF(stderr, "x[%d] = %g, x[%d+1] = %g\n\n", ie, x[ie], ie, x[ie+1]);
#endif
          } break;
          default:
            WH(-1, "predictor not yet implemented for this type of XFEM");
          }
        }
      }
      ioffset += vd->Ndof;
    }
  }

  /* fix xdot */
  xfem_correct(num_total_nodes, x, xdot, x_old, xdot_old, NULL, theta_arg,
               delta_t);
}

/******************************************************************************
 * xfem_var_diff: compute difference in var
 *                vdiff = the difference between the value on the other side
 *                        of the interface and the value on this side
 *                phidiff = d(vdiff)/vj
 *                gradvdiff = grad(vdiff)
 *
 ******************************************************************************/
void xfem_var_diff(int var, double *vdiff, double phidiff[MDE],
                   double gradvdiff[DIM]) {
  int interp = pd->i[pg->imtrx][var];
  BASIS_FUNCTIONS_STRUCT *bfv = bf[var];
  int xfem_active, extended_dof, base_interp, base_dof;
  int a, i;
  double *esp;

  /* initialize */
  *vdiff = 0.;
  gradvdiff[0] = 0.;
  gradvdiff[1] = 0.;
  gradvdiff[2] = 0.;
  for (i = 0; i < ei[pg->imtrx]->dof[var]; i++)
    phidiff[i] = 0.;

  switch (interp) {
  case I_P0_XV:
  case I_Q1_XV:
  case I_Q2_XV: {
    int sign = -ls->Elem_Sign;
    for (i = 0; i < ei[pg->imtrx]->dof[var]; i++) {
      xfem_dof_state(i, interp, ei[pg->imtrx]->ielem_shape, &xfem_active,
                     &extended_dof, &base_interp, &base_dof);
      if (extended_dof) {
        esp = x_static +
              ei[pg->imtrx]->ieqn_ledof[ei[pg->imtrx]->lvdof_to_ledof[var][i]];

        *vdiff += bfv->phi[2 * base_dof] * *esp * sign;
        phidiff[i] = bfv->phi[2 * base_dof] * sign;

        for (a = 0; a < VIM; a++) {
          gradvdiff[a] += bfv->grad_phi[2 * base_dof][a] * *esp * sign;
        }
      }
    }
  } break;
  case I_NOTHING:
    break;
  default:
    EH(-1, "xfem_var_diff not yet implemented for this type of XFEM");
  }
}

/******************************************************************************
 * zero_lsi: Initialize a level set interface data structure
 *
 * Author: Pat Notz 10/29/01
 ******************************************************************************/
void zero_lsi(void) {
  /* Sanity checking. */
  if (lsi == NULL) {
    EH(-1, "lsi Level_Set_Interface structure is NULL.");
  }

  if (ls == NULL) {
    EH(-1, "ls Level_Set_Data structure is NULL.");
  }

  lsi->near = FALSE;
  lsi->alpha = 0.0;

  lsi->H = 0.0;
  lsi->delta = 0.0;

  memset(lsi->normal, 0, sizeof(double) * DIM);

  /* This is useful for calculating the above (and other) quantities. */
  lsi->gfmag = 0.0;
  lsi->gfmaginv = 0.0;
}

void zero_lsi_derivs(void) {
  /* Sanity checking. */
  if (lsi == NULL) {
    EH(-1, "lsi Level_Set_Interface structure is NULL.");
  }

  if (ls == NULL) {
    EH(-1, "ls Level_Set_Data structure is NULL.");
  }

  lsi->dH = 0.;
  memset(lsi->d_H_dF, 0, sizeof(double) * MDE);
  memset(lsi->d_H_dmesh, 0, sizeof(double) * DIM * MDE);

  memset(lsi->d_delta_dF, 0, sizeof(double) * MDE);
  memset(lsi->d_delta_dmesh, 0, sizeof(double) * DIM * MDE);

  memset(lsi->d_normal_dF, 0, sizeof(double) * DIM * MDE);
  memset(lsi->d_normal_dmesh, 0, sizeof(double) * DIM * DIM * MDE);

  /* This is useful for calculating the above (and other) quantities. */
  memset(lsi->d_gfmag_dF, 0, sizeof(double) * MDE);
  memset(lsi->d_gfmag_dmesh, 0, sizeof(double) * DIM * MDE);

  memset(lsi->d_gfmaginv_dF, 0, sizeof(double) * MDE);
  memset(lsi->d_gfmaginv_dmesh, 0, sizeof(double) * DIM * MDE);
}

/******************************************************************************
 * load_lsi: Load the level set interface functions into the global
 *           lsi (Level_Set_Interface struct) based on the current state
 *           of the global fv (Field_Variables) data structure (viz., fv->F).
 *
 * Input
 * =====
 * width = The width of the interfacial region.
 *
 * Returns 0 on success.
 *
 * Author: Pat Notz 10/29/01
 ******************************************************************************/
int load_lsi(const double width) {
  double F = 0, alpha, *grad_F = NULL;
  int a, b;
  int i, j, k;

  /* Zero things out. */
  zero_lsi();

  /* Check if we're in the mushy zone. */
  lsi->alpha = 0.5 * width;
  alpha = lsi->alpha;

  copy_distance_function( &F, &grad_F);

  lsi->near = ls->on_sharp_surf || fabs(F) < alpha;

  /* Calculate the interfacial functions we want to know even if not in mushy
   * zone. */

  lsi->gfmag = 0.0;
  for (a = 0; a < VIM; a++) {
    lsi->normal[a] = grad_F[a];
    lsi->gfmag += grad_F[a] * grad_F[a];
  }
  lsi->gfmag = sqrt(lsi->gfmag);
  lsi->gfmaginv = (lsi->gfmag == 0.0) ? 1.0 : 1.0 / lsi->gfmag;

  for (a = 0; a < VIM; a++) {
    lsi->normal[a] *= lsi->gfmaginv;
  }

  /* If we're not in the mushy zone: */
  if (ls->on_sharp_surf) {
    /*lsi->H = ( F < 0.0) ? 0.0 : 1.0 ;*/
    lsi->H = (ls->Elem_Sign < 0) ? 0.0 : 1.0;
    lsi->delta = 1.;
  } else if (!lsi->near) {
    lsi->H = (F < 0.0) ? 0.0 : 1.0;
    lsi->delta = 0.;
  } else {
    lsi->H = 0.5 * (1. + F / alpha + sin(M_PIE * F / alpha) / M_PIE);
    lsi->delta = 0.5 * (1. + cos(M_PIE * F / alpha)) * lsi->gfmag / alpha;
  }

  /**** Shield the operations below since they are very expensive relative to
     the previous operations in the load_lsi routine. Add your variables as
     needed  ********/

  if (pd->gv[LUBP] || pd->gv[LUBP_2] ||
      pd->gv[SHELL_SAT_CLOSED] ||
      pd->gv[SHELL_PRESS_OPEN] ||
      pd->gv[SHELL_PRESS_OPEN_2] ||
      pd->gv[SHELL_SAT_GASN]) {

    /* Evaluate heaviside using FEM basis functions */
    double Hni, d_Hni_dF, Fi;
      double Hni_old, Fi_old;
    int eqn = R_FILL;
    lsi->Hn = 0.0;
      lsi->Hn_old = 0.0;
    memset(lsi->gradHn, 0.0, sizeof(double) * DIM);
      memset(lsi->gradHn_old, 0.0, sizeof(double)*DIM);
    memset(lsi->d_Hn_dF, 0.0, sizeof(double) * MDE);
    memset(lsi->d_gradHn_dF, 0.0, sizeof(double) * DIM * MDE);
    memset(lsi->d_Hn_dmesh, 0.0, sizeof(double) * DIM * MDE);
    memset(lsi->d_gradHn_dmesh, 0.0, sizeof(double) * DIM * DIM * MDE);
    if(pd->gv[LUBP] || pd->gv[SHELL_SAT_CLOSED] || pd->gv[SHELL_PRESS_OPEN ] || pd->gv[SHELL_SAT_GASN] ) 
      {
	for ( i = 0; i < ei[pg->imtrx]->dof[eqn]; i++ ) {
	  Fi = *esp->F[i];
	  if ( fabs(Fi) > lsi->alpha ) {
	    Hni = ( Fi < 0.0 ) ? 0.0 : 1.0;
	    d_Hni_dF = 0.0;
	  } else {
	    Hni      = 0.5 * (1.0 + Fi/lsi->alpha + sin(M_PIE*Fi/lsi->alpha)/M_PIE);
	    d_Hni_dF = 0.5 * (1/lsi->alpha + cos(M_PIE*Fi/lsi->alpha)/lsi->alpha);
	  }
	  lsi->Hn         += Hni      * bf[eqn]->phi[i];
	  lsi->d_Hn_dF[i] += d_Hni_dF * bf[eqn]->phi[i];
	  if (pd->gv[MESH_DISPLACEMENT1]) {
	    for ( b = 0; b < DIM; b++ ) {
	      for ( k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++ ) {
		lsi->d_Hn_dmesh[b][k] += Hni * bf[eqn]->phi[i] * bf[MESH_DISPLACEMENT1]->phi[k];
	      }
	    }
	  }
	  for ( j = 0; j < VIM; j++ ) {
	    lsi->gradHn[j]         += Hni      * bf[eqn]->grad_phi[i][j];
	    lsi->d_gradHn_dF[j][i] += d_Hni_dF * bf[eqn]->grad_phi[i][j];
	    if (pd->gv[MESH_DISPLACEMENT1]) {
	      for ( b = 0; b < DIM; b++ ) {
		for ( k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++ ) {
		  lsi->d_gradHn_dmesh[j][b][k] += Hni * bf[eqn]->d_grad_phi_dmesh[i][j][b][k];
		}
	      }
	    }
	  }

	  Fi_old = *esp_old->F[i];
	  if ( fabs(Fi_old) > lsi->alpha ) {
	    Hni_old = ( Fi_old < 0.0 ) ? 0.0 : 1.0;
	  } else {
	    Hni_old  = 0.5 * (1.0 + Fi_old/lsi->alpha + sin(M_PIE*Fi_old/lsi->alpha)/M_PIE);
	  }
	  lsi->Hn_old += Hni_old * bf[eqn]->phi[i];
	  for ( j = 0; j < VIM; j++ ) {
	    lsi->gradHn_old[j] += Hni_old * bf[eqn]->grad_phi[i][j];
	  }
	}
      }
    else if (pd->gv[LUBP_2] || pd->gv[SHELL_PRESS_OPEN_2])
      {
	eqn = R_PHASE1;
	for ( i = 0; i < ei[pg->imtrx]->dof[eqn]; i++ ) {
	  Fi = *esp->pF[0][i];
	  if ( fabs(Fi) > lsi->alpha ) {
	    Hni = ( Fi < 0.0 ) ? 0.0 : 1.0;
	    d_Hni_dF = 0.0;
	  } else {
	    Hni      = 0.5 * (1.0 + Fi/lsi->alpha + sin(M_PIE*Fi/lsi->alpha)/M_PIE);
	    d_Hni_dF = 0.5 * (1/lsi->alpha + cos(M_PIE*Fi/lsi->alpha)/lsi->alpha);
	  }
	  lsi->Hn         += Hni      * bf[eqn]->phi[i];
	  lsi->d_Hn_dF[i] += d_Hni_dF * bf[eqn]->phi[i];
	  if (pd->gv[MESH_DISPLACEMENT1]) {
	    for ( b = 0; b < DIM; b++ ) {
	      for ( k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++ ) {
		lsi->d_Hn_dmesh[b][k] += Hni * bf[eqn]->phi[i] * bf[MESH_DISPLACEMENT1]->phi[k];
	      }
	    }
	  }
	  for ( j = 0; j < VIM; j++ ) {
	    lsi->gradHn[j]         += Hni      * bf[eqn]->grad_phi[i][j];
	    lsi->d_gradHn_dF[j][i] += d_Hni_dF * bf[eqn]->grad_phi[i][j];
	    if (pd->gv[MESH_DISPLACEMENT1]) {
	      for ( b = 0; b < DIM; b++ ) {
		for ( k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++ ) {
		  lsi->d_gradHn_dmesh[j][b][k] += Hni * bf[eqn]->d_grad_phi_dmesh[i][j][b][k];
		}
	      }
	    }
	  }

	  Fi_old = *esp_old->pF[0][i];
	  if ( fabs(Fi_old) > lsi->alpha ) {
	    Hni_old = ( Fi_old < 0.0 ) ? 0.0 : 1.0;
	  } else {
	    Hni_old  = 0.5 * (1.0 + Fi_old/lsi->alpha + sin(M_PIE*Fi_old/lsi->alpha)/M_PIE);
	  }
	  lsi->Hn_old += Hni_old * bf[eqn]->phi[i];
	  for ( j = 0; j < VIM; j++ ) {
	    lsi->gradHn_old[j] += Hni_old * bf[eqn]->grad_phi[i][j];
	  }
	}
      } 
 
  } /* end of if pd->v[LUBP] || ... etc */

  /************ End of shielding **************************/

  if (fabs(alpha) < 1e-32) {
    alpha = 1e-32;
  }
  lsi->delta_max = lsi->gfmag / alpha;

  return(0);

}

int
load_lsi_old(const double width, struct Level_Set_Interface *lsi_old)
{
  double F_old = 0, alpha, *grad_F_old = NULL;
  int a;

  if (ls->var != FILL) {
    EH(-1, "Unknown level set variable");
  }

  lsi_old->near  = FALSE;
  lsi_old->alpha = 0.0;

  lsi_old->H = 0.0;
  lsi_old->delta = 0.0;

  memset(lsi_old->normal, 0, sizeof(double)*DIM);

  /* This is useful for calculating the above (and other) quantities. */
  lsi_old->gfmag = 0.0;
  lsi_old->gfmaginv = 0.0;

  /* Check if we're in the mushy zone. */
  lsi_old->alpha = 0.5 * width;
  alpha      = lsi_old->alpha;

  F_old = fv_old->F;
  grad_F_old = fv_old->grad_F;

  lsi_old->near  = ls->on_sharp_surf || fabs(F_old) < alpha;

  /* Calculate the interfacial functions we want to know even if not in mushy zone. */

  lsi_old->gfmag = 0.0;
  for ( a=0; a < VIM; a++ )
    {
      lsi_old->normal[a] = grad_F_old[a];
      lsi_old->gfmag    += grad_F_old[a] * grad_F_old[a];
    }
  lsi_old->gfmag = sqrt( lsi_old->gfmag );
  lsi_old->gfmaginv     = ( lsi_old->gfmag == 0.0 ) ? 1.0 : 1.0 / lsi_old->gfmag;

  for ( a=0; a < VIM; a++)
    {
      lsi_old->normal[a] *= lsi_old->gfmaginv;
    }

  /* If we're not in the mushy zone: */
  if ( ls->on_sharp_surf )
    {
      lsi_old->H = ( ls->Elem_Sign < 0 ) ? 0.0 : 1.0 ;
      lsi_old->delta = 1.;
    }
  else if ( ! lsi_old->near )
    {
      lsi_old->H = ( F_old < 0.0) ? 0.0 : 1.0 ;
      lsi_old->delta = 0.;
    }
  else
    {
      lsi_old->H     = 0.5 * (1. + F_old / alpha + sin(M_PIE * F_old / alpha) / M_PIE);
      lsi_old->delta = 0.5 * (1. + cos(M_PIE * F_old / alpha)) * lsi_old->gfmag / alpha;
    }


/**** Shield the operations below since they are very expensive relative to the previous
      operations in the load_lsi routine. Add your variables as needed  ********/

  if (pd->v[pg->imtrx][LUBP]  || pd->v[pg->imtrx][LUBP_2] || pd->v[pg->imtrx][SHELL_SAT_CLOSED] || pd->v[pg->imtrx][SHELL_PRESS_OPEN ] ||
      pd->v[pg->imtrx][SHELL_PRESS_OPEN_2] || pd->v[pg->imtrx][SHELL_SAT_GASN] )
    {
      EH(-1, "No support for LUBP/SHELL_SAT/SHELL_PRESS");
    } /* end of if pd->v[pg->imtrx][LUBP] || ... etc */

/************ End of shielding **************************/

  lsi_old->delta_max = lsi_old->gfmag/alpha;

  return (0);
}

int load_lsi_offset(const double width) {
  double F = 0, alpha, *grad_F = NULL;
  int a, b;
  int i, j, k;

  /* Zero things out. */
  zero_lsi();

  /* Check if we're in the mushy zone. */
  lsi->alpha = 0.5 * width;
  alpha = lsi->alpha;

  copy_distance_function(&F, &grad_F);

  // add offset
  F += 2 * alpha;

  lsi->near = ls->on_sharp_surf || fabs(F) < alpha;

  /* Calculate the interfacial functions we want to know even if not in mushy
   * zone. */

  lsi->gfmag = 0.0;
  for (a = 0; a < VIM; a++) {
    lsi->normal[a] = grad_F[a];
    lsi->gfmag += grad_F[a] * grad_F[a];
  }
  lsi->gfmag = sqrt(lsi->gfmag);
  lsi->gfmaginv = (lsi->gfmag == 0.0) ? 1.0 : 1.0 / lsi->gfmag;

  for (a = 0; a < VIM; a++) {
    lsi->normal[a] *= lsi->gfmaginv;
  }

  /* If we're not in the mushy zone: */
  if (ls->on_sharp_surf) {
    /*lsi->H = ( F < 0.0) ? 0.0 : 1.0 ;*/
    lsi->H = (ls->Elem_Sign < 0) ? 0.0 : 1.0;
    lsi->delta = 1.;
  } else if (!lsi->near) {
    lsi->H = (F < 0.0) ? 0.0 : 1.0;
    lsi->delta = 0.;
  } else {
    lsi->H = 0.5 * (1. + F / alpha + sin(M_PIE * F / alpha) / M_PIE);
    lsi->delta = 0.5 * (1. + cos(M_PIE * F / alpha)) * lsi->gfmag / alpha;
  }

  /**** Shield the operations below since they are very expensive relative to
     the previous operations in the load_lsi routine. Add your variables as
     needed  ********/

  if (pd->v[pg->imtrx][LUBP] || pd->v[pg->imtrx][LUBP_2] ||
      pd->v[pg->imtrx][SHELL_SAT_CLOSED] ||
      pd->v[pg->imtrx][SHELL_PRESS_OPEN] ||
      pd->v[pg->imtrx][SHELL_PRESS_OPEN_2] ||
      pd->v[pg->imtrx][SHELL_SAT_GASN]) {

    /* Evaluate heaviside using FEM basis functions */
    double Hni, d_Hni_dF, Fi;
    int eqn = R_FILL;
    lsi->Hn = 0.0;
    memset(lsi->gradHn, 0.0, sizeof(double) * DIM);
    memset(lsi->d_Hn_dF, 0.0, sizeof(double) * MDE);
    memset(lsi->d_gradHn_dF, 0.0, sizeof(double) * DIM * MDE);
    memset(lsi->d_Hn_dmesh, 0.0, sizeof(double) * DIM * MDE);
    memset(lsi->d_gradHn_dmesh, 0.0, sizeof(double) * DIM * DIM * MDE);

    if (pd->v[pg->imtrx][LUBP] || pd->v[pg->imtrx][SHELL_SAT_CLOSED] ||
        pd->v[pg->imtrx][SHELL_PRESS_OPEN] ||
        pd->v[pg->imtrx][SHELL_SAT_GASN]) {
      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        Fi = *esp->F[i];
        if (fabs(Fi) > lsi->alpha) {
          Hni = (Fi < 0.0) ? 0.0 : 1.0;
          d_Hni_dF = 0.0;
        } else {
          Hni = 0.5 *
                (1.0 + Fi / lsi->alpha + sin(M_PIE * Fi / lsi->alpha) / M_PIE);
          d_Hni_dF = 0.5 * (1 / lsi->alpha +
                            cos(M_PIE * Fi / lsi->alpha) / lsi->alpha);
        }
        lsi->Hn += Hni * bf[eqn]->phi[i];
        lsi->d_Hn_dF[i] += d_Hni_dF * bf[eqn]->phi[i];
        if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
          for (b = 0; b < DIM; b++) {
            for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
              lsi->d_Hn_dmesh[b][k] +=
                  Hni * bf[eqn]->phi[i] * bf[MESH_DISPLACEMENT1]->phi[k];
            }
          }
        }
        for (j = 0; j < VIM; j++) {
          lsi->gradHn[j] += Hni * bf[eqn]->grad_phi[i][j];
          lsi->d_gradHn_dF[j][i] += d_Hni_dF * bf[eqn]->grad_phi[i][j];
          if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
            for (b = 0; b < DIM; b++) {
              for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
                lsi->d_gradHn_dmesh[j][b][k] +=
                    Hni * bf[eqn]->d_grad_phi_dmesh[i][j][b][k];
              }
            }
          }
        }
      }
    } else if (pd->v[pg->imtrx][LUBP_2] ||
               pd->v[pg->imtrx][SHELL_PRESS_OPEN_2]) {
      eqn = R_PHASE1;
      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        Fi = *esp->pF[0][i];
        if (fabs(Fi) > lsi->alpha) {
          Hni = (Fi < 0.0) ? 0.0 : 1.0;
          d_Hni_dF = 0.0;
        } else {
          Hni = 0.5 *
                (1.0 + Fi / lsi->alpha + sin(M_PIE * Fi / lsi->alpha) / M_PIE);
          d_Hni_dF = 0.5 * (1 / lsi->alpha +
                            cos(M_PIE * Fi / lsi->alpha) / lsi->alpha);
        }
        lsi->Hn += Hni * bf[eqn]->phi[i];
        lsi->d_Hn_dF[i] += d_Hni_dF * bf[eqn]->phi[i];
        if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
          for (b = 0; b < DIM; b++) {
            for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
              lsi->d_Hn_dmesh[b][k] +=
                  Hni * bf[eqn]->phi[i] * bf[MESH_DISPLACEMENT1]->phi[k];
            }
          }
        }
        for (j = 0; j < VIM; j++) {
          lsi->gradHn[j] += Hni * bf[eqn]->grad_phi[i][j];
          lsi->d_gradHn_dF[j][i] += d_Hni_dF * bf[eqn]->grad_phi[i][j];
          if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
            for (b = 0; b < DIM; b++) {
              for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
                lsi->d_gradHn_dmesh[j][b][k] +=
                    Hni * bf[eqn]->d_grad_phi_dmesh[i][j][b][k];
              }
            }
          }
        }
      }
    }

  } /* end of if pd->v[pg->imtrx][LUBP] || ... etc */

  /************ End of shielding **************************/

  lsi->delta_max = lsi->gfmag / alpha;

  return (0);
}

/******************************************************************************
 * load_lsi_shell_second: Special case of load_lsi forced to be written to
 *circumvent the shell-shell-friend situation encountered in multilayer shell
 *stacks.  The friend code does not work for this pathological case, so we need
 *to do nonlocal operations and cannot do a setup-shop-at-point approach.
 *  *
 * Input
 * =====
 * width = The width of the interfacial region.
 *
 * Returns 0 on success.
 *
 * Author: Randy Schunk, dubiously, on 8/21/2012
 ******************************************************************************/
int load_lsi_shell_second(const double width) {
  double F = 0, alpha, *grad_F = NULL;
  int a, b;
  int i, j, k;

  /* Zero things out. */
  zero_lsi();

  /* Check if we're in the mushy zone. */
  lsi->alpha = 0.5 * width;
  alpha = lsi->alpha;

  copy_distance_function(&F, &grad_F);

  lsi->near = ls->on_sharp_surf || fabs(F) < alpha;

  /* Calculate the interfacial functions we want to know even if not in mushy
   * zone. */

  lsi->gfmag = 0.0;
  for (a = 0; a < VIM; a++) {
    lsi->normal[a] = grad_F[a];
    lsi->gfmag += grad_F[a] * grad_F[a];
  }
  lsi->gfmag = sqrt(lsi->gfmag);
  lsi->gfmaginv = (lsi->gfmag == 0.0) ? 1.0 : 1.0 / lsi->gfmag;

  for (a = 0; a < VIM; a++) {
    lsi->normal[a] *= lsi->gfmaginv;
  }

  /* If we're not in the mushy zone: */
  if (ls->on_sharp_surf) {
    /*lsi->H = ( F < 0.0) ? 0.0 : 1.0 ;*/
    lsi->H = (ls->Elem_Sign < 0) ? 0.0 : 1.0;
    lsi->delta = 1.;
  } else if (!lsi->near) {
    lsi->H = (F < 0.0) ? 0.0 : 1.0;
    lsi->delta = 0.;
  } else {
    lsi->H = 0.5 * (1. + F / alpha + sin(M_PIE * F / alpha) / M_PIE);
    lsi->delta = 0.5 * (1. + cos(M_PIE * F / alpha)) * lsi->gfmag / alpha;
  }

  /**** Shield the operations below since they are very expensive relative to
     the previous operations in the load_lsi routine. Add your variables as
     needed.   The balance force
        formulation here.   compute at nodes first, and then spread to gauss
     points ********/

  double Hni, d_Hni_dF, Fi;
  int eqn;
  lsi->Hn = 0.0;
  memset(lsi->gradHn, 0.0, sizeof(double) * DIM);
  memset(lsi->d_Hn_dF, 0.0, sizeof(double) * MDE);
  memset(lsi->d_gradHn_dF, 0.0, sizeof(double) * DIM * MDE);
  memset(lsi->d_Hn_dmesh, 0.0, sizeof(double) * DIM * MDE);
  memset(lsi->d_gradHn_dmesh, 0.0, sizeof(double) * DIM * DIM * MDE);
  if (upd->vp[pg->imtrx][LUBP_2] >= 0 ||
      upd->vp[pg->imtrx][SHELL_PRESS_OPEN_2] >= 0) {
    eqn = R_PHASE1;
    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      Fi = *esp->pF[0][i];
      if (fabs(Fi) > lsi->alpha) {
        Hni = (Fi < 0.0) ? 0.0 : 1.0;
        d_Hni_dF = 0.0;
      } else {
        Hni = 0.5 *
              (1.0 + Fi / lsi->alpha + sin(M_PIE * Fi / lsi->alpha) / M_PIE);
        d_Hni_dF =
            0.5 * (1 / lsi->alpha + cos(M_PIE * Fi / lsi->alpha) / lsi->alpha);
      }
      lsi->Hn += Hni * bf[eqn]->phi[i];
      lsi->d_Hn_dF[i] += d_Hni_dF * bf[eqn]->phi[i];
      if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
        for (b = 0; b < DIM; b++) {
          for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
            lsi->d_Hn_dmesh[b][k] +=
                Hni * bf[eqn]->phi[i] * bf[MESH_DISPLACEMENT1]->phi[k];
          }
        }
      }
      for (j = 0; j < VIM; j++) {
        lsi->gradHn[j] += Hni * bf[eqn]->grad_phi[i][j];
        lsi->d_gradHn_dF[j][i] += d_Hni_dF * bf[eqn]->grad_phi[i][j];
        if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
          for (b = 0; b < DIM; b++) {
            for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
              lsi->d_gradHn_dmesh[j][b][k] +=
                  Hni * bf[eqn]->d_grad_phi_dmesh[i][j][b][k];
            }
          }
        }
      }
    }
  } /* end of if upd->vp[pg->imtrx][LUBP] || ... etc */
  else {
    EH(-1, " you shouldn't be in this routine. Go check it out or contact PRS "
           "8/21/2012");
  }

  /************ End of shielding **************************/

  lsi->delta_max = lsi->gfmag / alpha;

  return (0);
}

static void copy_distance_function(double *F, double **grad_F) {
  int offset = 0;

  switch (ls->var) {
  case FILL:
    *F = fv->F;
    *grad_F = fv->grad_F;
    break;
  case PHASE1:
  case PHASE2:
  case PHASE3:
  case PHASE4:
  case PHASE5:
    offset = ls->var - PHASE1;
    *F = fv->pF[offset];
    *grad_F = fv->grad_pF[offset];
    break;
  default:
    EH(-1, " Unknown distance function variable type.\n");
    break;
  }
}

/******************************************************************************
 * load_lsi_derivs Load up the derivatives of the level set interface
 *                  functions into the global lsi (Level_Set_Interface struct)
 *                  based on the current state of the global fv
 *                  (Field_Variables) data structure.  This function should
 *                  only be called after load_lsi().
 *
 * Returns 0 on success.
 *
 * Author: Pat Notz 10/29/01
 ******************************************************************************/
int load_lsi_derivs(void) {
  double F = 0, phi_j, grad_phi_j[DIM], *grad_F = NULL;
  double alpha = lsi->alpha;
  int a, b, j, var;

  /* Zero things out. */
  zero_lsi_derivs();

  copy_distance_function(&F, &grad_F);

  /* Initialize grad_phi_j */
  for (a = 0; a < DIM; a++) {
    grad_phi_j[a] = 0;
  }

  /*
   * If we're here, the pd->v[pg->imtrx]ar[ls->var] is true...
   *
   * Always compute the distance function variable derivs, even for uncoupled
   * fill problems... see Hrenorm_constrain
   *
   * Derivatives w.r.t. distance function variable
   */
  var = ls->var;
  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
    /* Fetch the basis functions. */
    phi_j = bf[var]->phi[j];
    for (a = 0; a < VIM; a++) {
      grad_phi_j[a] = bf[var]->grad_phi[j][a];
    }

    /* Derivative of gfmag. */
    lsi->d_gfmag_dF[j] = 0.0;
    for (a = 0; a < VIM; a++) {
      lsi->d_gfmag_dF[j] += grad_F[a] * grad_phi_j[a] * lsi->gfmaginv;
    }

    /* Derivative of gfmaginv. */
    lsi->d_gfmaginv_dF[j] = (lsi->gfmag == 0.0)
                                ? 0.0
                                : -pow(lsi->gfmaginv, 2.0) * lsi->d_gfmag_dF[j];

    /* Derivative of the normal vector. */
    for (a = 0; a < VIM; a++) {
      lsi->d_normal_dF[a][j] =
          grad_phi_j[a] * lsi->gfmaginv + lsi->d_gfmaginv_dF[j] * grad_F[a];
    }

  } /* for: j */

  /*
   * Derivatives w.r.t. MESH_DISPLACEMENTs
   */
  if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {

    for (b = 0; b < VIM; b++) {
      var = MESH_DISPLACEMENT1 + b;
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        /* Fetch the basis functions. */
        phi_j = bf[var]->phi[j];
        for (a = 0; a < VIM; a++) {
          grad_phi_j[a] = bf[var]->grad_phi[j][a];
        }

        /* gfmag */
        for (a = 0; a < VIM; a++) {
          lsi->d_gfmag_dmesh[b][j] =
              grad_F[a] * fv->d_grad_F_dmesh[a][b][j] * lsi->gfmaginv;
        }

        /* gfmaginv */
        lsi->d_gfmaginv_dmesh[b][j] =
            (lsi->gfmag == 0.0)
                ? 0.0
                : -pow(lsi->gfmaginv, 2.0) * lsi->d_gfmag_dmesh[b][j];

        /* normal */
        for (a = 0; a < VIM; a++) {
          lsi->d_normal_dmesh[a][b][j] =
              fv->d_grad_F_dmesh[a][b][j] * lsi->gfmaginv +
              lsi->d_gfmaginv_dmesh[b][j] * grad_F[a];
        }

      } /* for: j */

    } /* for: b */
  }

  /* DRN: this is required to get path dependence terms right */
  var = FILL;
  if (ls->on_sharp_surf) {
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      /* Fetch the basis functions. */
      phi_j = bf[var]->phi[j];

      /* Derivative of the delta function. */
      lsi->d_delta_dF[j] = lsi->d_gfmag_dF[j] * lsi->gfmaginv;
      lsi->d_H_dF[j] = phi_j * lsi->gfmaginv;
    }
  }

  /* If we're not in the mushy zone, all remaining derivs should be zero. */
  if (ls->on_sharp_surf || !lsi->near)
    return (0);

  lsi->dH = 0.5 * (1.0 / alpha) * (1. + cos(M_PIE * F / alpha));

  /*
   * Derivatives w.r.t. FILL for non-zero alpha
   */
  var = ls->var;
  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
    /* Fetch the basis functions. */
    phi_j = bf[var]->phi[j];

    /* Derivative of the H function. */
    lsi->d_H_dF[j] = phi_j * lsi->dH;

    /* Derivative of the delta function. */
    lsi->d_delta_dF[j] =
        -0.5 * (M_PIE * phi_j / alpha * sin(M_PIE * F / alpha)) * lsi->gfmag /
            alpha +
        0.5 * (1. + cos(M_PIE * F / alpha)) * lsi->d_gfmag_dF[j] / alpha;

  } /* for: j */

  /*
   * Derivatives w.r.t. MESH_DISPLACEMENTs for non-zero alpha
   */
  if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
    for (b = 0; b < VIM; b++) {
      var = MESH_DISPLACEMENT1 + b;
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

        /* In the current implementation, H does not depend on grad_F
           and, hence, it doens't depend on the mesh. */

        /* delta */
        lsi->d_delta_dmesh[b][j] = 0.5 * (1. + cos(M_PIE * F / alpha)) / alpha *
                                   lsi->d_gfmag_dmesh[b][j];

      } /* for: j */

    } /* for: b */
  }

  return (0);
}

/*
 *
 *  load_lsi_adjmatr
 *
 *    Function that computes value of level set function based on values
 *    from an adjacent material under the assumption that the current material
 *    does not contain level set function unknowns
 *
 *    Not that only the value of the level set function and values derived from
 * it (H, delta, d_H_dF
 */

int load_lsi_adjmatr(const double width) {
  double F, alpha, phi_j;
  int j, ln;
  int lvdesc, lvdof, var_type, num_dofs;

  lvdesc = 0;
  while ((var_type = ei[pg->imtrx]->Lvdesc_to_Var_Type[lvdesc]) != FILL)
    lvdesc++;

  num_dofs = ei[pg->imtrx]->Lvdesc_Numdof[lvdesc];

  fv->F = scalar_fv_fill_adjmatrl(esp->F, lvdesc, num_dofs, var_type);

  /* Zero things out. */
  zero_lsi();
  zero_lsi_derivs();

  /* Check if we're in the mushy zone. */

  lsi->alpha = 0.5 * width;
  alpha = lsi->alpha;
  F = fv->F;
  lsi->near = ls->on_sharp_surf || fabs(F) < alpha;

  /*
   * Note that we can't compute grad_F in the adjmatr (Yet)  so all quantities
   * related cannot be determined
   */

  /* If we're not in the mushy zone: */
  if (ls->on_sharp_surf) {
    /*lsi->H = ( F < 0.0) ? 0.0 : 1.0 ;*/
    lsi->H = (ls->Elem_Sign < 0) ? 0.0 : 1.0;
    lsi->delta = 1.;
  } else if (!lsi->near) {
    lsi->H = (F < 0.0) ? 0.0 : 1.0;
    lsi->delta = 0.;
  } else {
    lsi->H = 0.5 * (1. + F / alpha + sin(M_PIE * F / alpha) / M_PIE);
    lsi->delta =
        0.5 * (1. + cos(M_PIE * F / alpha)) /
        alpha; /* this is an approximation since we can't compute gfmag */
    lsi->dH = 0.5 * (1.0 / alpha) * (1. + cos(M_PIE * F / alpha));
    for (j = 0; j < num_dofs; j++) {

      lvdof = ei[pg->imtrx]->Lvdesc_to_lvdof[lvdesc][j];

      ln = ei[pg->imtrx]->dof_list[var_type][lvdof];

      phi_j = bf[var_type]->phi[ln];

      /* Derivative of the H function. */
      lsi->d_H_dF[j] = phi_j * lsi->dH;

      /* Derivative of the delta function. */
      lsi->d_delta_dF[j] =
          -0.5 * (M_PIE * phi_j / alpha * sin(M_PIE * F / alpha)) / alpha;
    }
  }

  return (0);
}

double ls_modulate_property(double p1, double p2, double width, double pm_minus,
                            double pm_plus, double dpdF[MDE], double *factor) {
  double p_plus, p_minus, p;

  p_minus = p1 * pm_plus + p2 * pm_minus;
  p_plus = p1 * pm_minus + p2 * pm_plus;

  level_set_property(p_minus, p_plus, width, &p, dpdF);

  if (ls->Elem_Sign == -1)
    *factor = pm_plus;
  else if (ls->Elem_Sign == 1)
    *factor = pm_minus;
  else
    *factor = pm_plus * (1.0 - lsi->H) + pm_minus * lsi->H;

  return (p);
}

double ls_modulate_property_offset(double p1, double p2, double width,
                                   double pm_minus, double pm_plus,
                                   double dpdF[MDE], double *factor) {
  double p_plus, p_minus, p;

  p_minus = p1 * pm_plus + p2 * pm_minus;
  p_plus = p1 * pm_minus + p2 * pm_plus;

  level_set_property_offset(p_minus, p_plus, width, &p, dpdF);

  if (ls->Elem_Sign == -1)
    *factor = pm_plus;
  else if (ls->Elem_Sign == 1)
    *factor = pm_minus;
  else
    *factor = pm_plus * (1.0 - lsi->H) + pm_minus * lsi->H;

  return (p);
}

static int current_elem_xfem_state(
    int node_var_state[], int *elem_var_state,
    double x[], /* Solution vector for the current processor    */
    const Exo_DB *exo) {
  /* check if there is a dependency on extended unknowns for this element */
  /* we also need to classify the nodes of this element:
     node_var_state[i] == 0 -> this node does not have active enriched dofs in
     any element node_var_state[i] == 1 -> this node has active enriched dofs in
     this element node_var_state[i] == 2 -> this node does not have active
     enriched dofs in this element but does in some other element
   */
  /* we also need to know the status of element vars:
     elem_var_state == 0 -> elem_vars in this element are not active
     elem_var_state == 1 -> elem_vars in this element are active
   */
  /* return = 0 -> no extending fns active in this element
     return = 1 -> at least one of the nodes of this element have node_var_state
     == 1 or elem_vars_state == 1 return = 2 -> at least one of the nodes of
     this element have node_var_state == 2
   */
  int i, j, e, I;
  int elem_state = 0;

  /* turn everything off by default */

  *elem_var_state = 0;
  for (i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
    node_var_state[i] = 0;
  }

  if (ls->Length_Scale != 0.) /* diffuse interface */
  {
    if (ls->elem_overlap_state) {
      /* element vars are active */
      *elem_var_state = 1;
      /* nodal vars are active */
      for (i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
        node_var_state[i] = 1;
      }
    } else {
      /* element vars are *NOT* active */
      *elem_var_state = 0;
      /* nodal vars still might be active */
      for (i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
        /* check all neighboring elements to see if any span the interface */
        I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
        for (j = exo->node_elem_pntr[I];
             j < exo->node_elem_pntr[I + 1] && !node_var_state[i]; j++) {
          e = exo->node_elem_list[j];
          if (elem_overlaps_interface(e, x, exo, ls->Length_Scale)) {
            node_var_state[i] = 2;
          }
        }
      }
    }
  } else /* sharp interface */
  {

    if (ls->elem_overlap_state) {
      /* element vars are active */
      *elem_var_state = 1;
      /* nodal vars are active */
      for (i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
        node_var_state[i] = 1;
      }
    } else {
      /* element vars are *NOT* active */
      *elem_var_state = 0;
      /* nodal vars still might be active */

      for (i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
        /* check all neighboring elements */
        I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
        for (j = exo->node_elem_pntr[I];
             j < exo->node_elem_pntr[I + 1] && !node_var_state[i]; j++) {
          e = exo->node_elem_list[j];
          if (elem_on_isosurface(e, x, exo, ls->var, 0.)) {
            node_var_state[i] = 2;
          }
        }
      }
    }
  }

  for (i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
    if (*elem_var_state == 1 || node_var_state[i] == 1)
      elem_state = 1;
    if (elem_state == 0 && node_var_state[i] == 2)
      elem_state = 2;
  }

  return (elem_state);
}

int is_xfem_interp(const int interp) {
  return (interp == I_P0_G || interp == I_P1_G || interp == I_Q1_G ||
          interp == I_Q2_G || interp == I_P0_GP || interp == I_P1_GP ||
          interp == I_Q1_GP || interp == I_Q2_GP || interp == I_P0_GN ||
          interp == I_P1_GN || interp == I_Q1_GN || interp == I_Q2_GN ||
          interp == I_P0_XV || interp == I_P1_XV || interp == I_Q1_XV ||
          interp == I_Q2_XV || interp == I_P1_XG || interp == I_Q1_XG ||
          interp == I_Q2_XG || interp == I_Q1_HG || interp == I_Q1_HV ||
          interp == I_Q1_HVG || interp == I_Q2_HG || interp == I_Q2_HV ||
          interp == I_Q2_HVG);
}

void xfem_dof_state(
    const int ledof, const int interpolation, const int eshape,
    int *xfem_active, /* flag indicating xfem affects this dof's basis functions
                       */
    int *extended_dof, /* flag indicating if this an extended dof */
    int *base_interp,  /* base interpolation, ie, I_Q1_XG -> I_Q1 */
    int *base_dof)     /* what dof of base_interp does this dof map to */
{
  switch (interpolation) {
  case I_P0_G:
  case I_P0_GP:
  case I_P0_GN:
  case I_P0_XV:
    *base_interp = I_P0;
    break;
  case I_P1_G:
  case I_P1_GP:
  case I_P1_GN:
  case I_P1_XV:
  case I_P1_XG:
    *base_interp = I_P1;
    break;
  case I_Q1_G:
  case I_Q1_GP:
  case I_Q1_GN:
  case I_Q1_XV:
  case I_Q1_XG:
  case I_Q1_HV:
  case I_Q1_HG:
  case I_Q1_HVG:
    *base_interp = I_Q1;
    break;
  case I_Q2_G:
  case I_Q2_GP:
  case I_Q2_GN:
  case I_Q2_XV:
  case I_Q2_XG:
  case I_Q2_HV:
  case I_Q2_HG:
  case I_Q2_HVG:
    *base_interp = I_Q2;
    break;
  default:
    /* assume non extended function and exit now */
    *base_interp = interpolation;
    *extended_dof = FALSE;
    *base_dof = ledof;
    *xfem_active = FALSE;
    return;
    break;
  }

  switch (interpolation) {
  case I_P0_G:
  case I_P0_XV:
  case I_P1_G:
  case I_P1_XV:
  case I_P1_XG:
  case I_Q1_G:
  case I_Q2_G:
  case I_Q1_XV:
  case I_Q1_XG:
  case I_Q2_XV:
  case I_Q2_XG:
    *extended_dof = (ledof % 2 == 1);
    if (*extended_dof)
      *base_dof = (ledof - 1) / 2;
    else
      *base_dof = ledof / 2;
    break;
  case I_Q1_HV:
  case I_Q1_HG:
  case I_Q1_HVG:
  case I_Q2_HV:
  case I_Q2_HG:
  case I_Q2_HVG:
    *extended_dof = (ledof >= getdofs(eshape, *base_interp));
    if (*extended_dof)
      *base_dof =
          ledof -
          getdofs(eshape, *base_interp); /* watch out used differently here */
    else
      *base_dof = ledof;
    break;
  case I_Q1_GP:
  case I_Q2_GP:
    *base_dof = ledof;
    *extended_dof = lnn_distance(*base_dof) < 0.;
    break;
  case I_P0_GP:
  case I_P1_GP:
    *base_dof = ledof;
    *extended_dof = xfem->elem_var_state == 1 || lnn_distance(*base_dof) < 0.;
    break;
  case I_Q1_GN:
  case I_Q2_GN:
    *base_dof = ledof;
    *extended_dof = lnn_distance(*base_dof) >= 0.;
    break;
  case I_P0_GN:
  case I_P1_GN:
    *base_dof = ledof;
    *extended_dof = xfem->elem_var_state == 1 || lnn_distance(*base_dof) >= 0.;
    break;
  default:
    EH(-1, "Unrecognized extended shape function.");
    break;
  }

  switch (interpolation) {
  case I_P0_G:
  case I_P0_XV:
  case I_P1_G:
  case I_P1_XV:
  case I_P1_XG:
  case I_P0_GP:
  case I_P0_GN:
  case I_Q1_HV:
  case I_Q1_HG:
  case I_Q1_HVG:
  case I_Q2_HV:
  case I_Q2_HG:
  case I_Q2_HVG:
    *xfem_active = (xfem->elem_var_state == 1);
    break;
  case I_Q1_G:
  case I_Q2_G:
  case I_Q1_GP:
  case I_Q1_GN:
  case I_Q2_GP:
  case I_Q2_GN:
  case I_Q1_XV:
  case I_Q1_XG:
  case I_Q2_XV:
  case I_Q2_XG:
    *xfem_active = (xfem->node_var_state[*base_dof] == 1);
    break;
  default:
    EH(-1, "Unrecognized extended shape function.");
    break;
  }

  return;
}

int is_extended_dof(const int I, const int idof,
                    VARIABLE_DESCRIPTION_STRUCT *vd, const double F) {
  /* DRN: this function isn't done yet */
  int var = vd->Variable_Type;
  int MatID = vd->MatID;
  int interp;
  int extended_dof = FALSE;

  if (MatID == -1)
    MatID = 0;
  interp = pd_glob[MatID]->i[pg->imtrx][var];

  switch (interp) {
  case I_P0_G:
  case I_P0_XV:
  case I_Q1_G:
  case I_Q2_G:
  case I_Q1_XV:
  case I_Q1_XG:
  case I_Q2_XV:
  case I_Q2_XG:
    extended_dof = (idof == 1);
    break;
  case I_P1_G:
  case I_P1_XV:
  case I_P1_XG:
    extended_dof = FALSE; /*FIX ME*/
    break;
  case I_Q1_HV:
  case I_Q1_HG:
  case I_Q1_HVG:
  case I_Q2_HV:
  case I_Q2_HG:
  case I_Q2_HVG:
    extended_dof = FALSE; /*FIX ME*/
                          /*
                           *extended_dof = (ledof >= getdofs( eshape, *base_interp ));
                           */
    break;
  case I_P0_GP:
  case I_P1_GP:
  case I_Q1_GP:
  case I_Q2_GP:
    extended_dof = (F < 0.);
    break;
  case I_P0_GN:
  case I_P1_GN:
  case I_Q1_GN:
  case I_Q2_GN:
    extended_dof = (F >= 0.);
    break;
  default:
    EH(-1, "Unrecognized extended shape function.");
    break;
  }

  return (extended_dof);
} /* END of function is_extended_dof  */

int dof_incomplete(int node, int elem_type, int interpolation, int eshape) {
  int iside, nodes_per_side;
  int lnn[MDE];
  double f[MDE];
  switch (interpolation) {
  case I_Q1:
    return TRUE;
  case I_Q2:
    switch (eshape) {
    case QUADRILATERAL:
      if (node < 4)
        return TRUE;
      if (node == 8)
        return FALSE;
      for (iside = 0; iside < 4; iside++) {
        get_side_info(elem_type, iside + 1, &nodes_per_side, lnn);
        if (lnn[2] == node) {
          f[0] = lnn_distance(lnn[0]);
          f[1] = lnn_distance(lnn[1]);
          f[2] = lnn_distance(lnn[2]);
          return (!sign_change(f[0], f[1]) && !sign_change(f[1], f[2]));
        }
      }
      EH(-1, "Unexpected Error.");
      break;
    default:
      EH(-1, "Not implemented yet.");
      break;
    }
    break;
  default:
    EH(-1, "Invalid shape function.");
    break;
  }
  return (-1);
} /* END of function dof_incomplete */

/* SEARCH GRID FUNCTIONS */
/*
 * The search grid method for finding points on the zero ls contour
 * is quite simple.  Once an element is determined to have a zero contour
 * section in it, it is subdivided into four (2D) or eight (3D) equal
 * sized smaller subelements ( called grids for lack of a better name ).
 * In s,u, and t space the element is split along planes passing through
 * the element center and parallel to the principle axis.  So, for example,
 * a quad element would generate four grids whose vertices are at the points:
 *  Grid 0:  (-1,-1), (0,-1), (0,0), (-1,0)
 *  Grid 1:  ( 0,-1), (1,-1), (1,0), ( 0,0)
 *  Grid 2:  ( 0,0 ), (1,0),  (1,1), ( 0,1)
 *  Grid 3:  (-1,0 0, (0,0),  (0,1), ( -1,1)
 *
 *  On each grid, the LS function is computed at the vertices ( call them nodes
 * ) and from this those grids through which the interface passes can be
 *  determined.  The other grids are ignored.  Those grids that have an
 * interface are divided once again by precisely the same method.  In fact, they
 * are divided by a recursive call to the function which divided the first
 * element This process continues until a preset level of division is reached.
 *
 *  The Search_Grid_Structure accomodates this division.  Each SGRID structure
 * represents a single grid.  It has members for the s,t,u coordinates of the
 * nodes and value for the LS function at the nodes.  It also has an array of
 * pointers that are populated by the grids that appear when the parent grid is
 * divided. Each of these subgrids might have subgrids of their own and so on.
 * Thus, each SGRID structure might be the starting node of an entire tree of
 * smaller subgrids.
 *
 */

SGRID *create_search_grid(NTREE *ntree) {
  SGRID *new_grid;

  new_grid = (SGRID *)smalloc(sizeof(SGRID));

  new_grid->ei = ei[pg->imtrx];
  /* new_grid->dim = pd->Num_Dim; PRS */
  new_grid->dim = ei[pg->imtrx]->ielem_dim;
  new_grid->level = 0;
  new_grid->num_verts = (new_grid->dim == 2) ? 4 : 8;
  new_grid->tree = ntree;

  switch (new_grid->dim) {
    int j;

  case 2:
    /*       new_grid->xi = ( double (*) [DIM] ) smalloc( sizeof(double)*4*DIM
     * ); */

    /*       new_grid->LS_value = smalloc( sizeof(double)*4 ); */

    for (j = 0; j < 4; j++) {
      find_nodal_stu(j, ei[pg->imtrx]->ielem_type, new_grid->tree->xi[j],
                     new_grid->tree->xi[j] + 1, new_grid->tree->xi[j] + 2);
      new_grid->LS_value[j] = *(esp->F[j]);
    }

    break;

  case 3:
    /*       new_grid->xi = ( double (*) [DIM] ) smalloc( sizeof(double)*8*DIM
     * ); */
    /*       new_grid->LS_value = smalloc( sizeof(double)*8 ); */

    for (j = 0; j < 8; j++) {
      find_nodal_stu(j, ei[pg->imtrx]->ielem_type, new_grid->tree->xi[j],
                     new_grid->tree->xi[j] + 1, new_grid->tree->xi[j] + 2);
      new_grid->LS_value[j] = *(esp->F[j]);
    }

    break;
  default:
    break;
  }
  return (new_grid);
}

/***************************************************************
 * divide_search_grid******************************************
 ***************************************************************
 *   input:   parent  -    pointer to search_grid structure
 *            max_level  - integer setting maximum degree of division requested
 *
 *   Nominally, allocs the subgrid array in the parent search_grid
 *   then populates each sub search grid with coordinates and LS values
 *   computed with respect to the parent.
 *
 *   It checks each sub search grid from the passage of the interface through
 *   it.  If no interface present the sub grid is deallocated.  If the interface
 *   is present, the sub search grid is further divided according to a recursive
 *   call to this same function.  When this function returns, the parent
 *structure pointer stands at the top of a branching tree structure extending
 *downward max_level - parent->level branches. Pretty cool huh ?
 *
 *  returns void.
 *********************************************************************/

void divide_search_grid(SGRID *parent, int max_level) {
  if (parent->level <
      max_level) /* If my level is the max level, return to end the recursion */
  {
    int index;
    int num_subgrids;
    int dim = parent->dim;

    num_subgrids = (dim == 2) ? 4 : 8;

    parent->num_subgrids = num_subgrids;

    /* allocated subgrid structures  */

    parent->subgrids = (SGRID **)smalloc(num_subgrids * sizeof(SGRID *));

    for (index = 0; index < num_subgrids; index++)

    {
      parent->subgrids[index] = (SGRID *)smalloc(sizeof(SGRID));

      /* Inherit so of the structure member values from the parent */

      parent->subgrids[index]->ei = parent->ei;
      parent->subgrids[index]->dim = parent->dim;
      parent->subgrids[index]->level = parent->level + 1;
      parent->subgrids[index]->num_verts = (parent->dim == 2) ? 4 : 8;
      parent->subgrids[index]->num_subgrids = 0;
      parent->subgrids[index]->subgrids = NULL;
      parent->subgrids[index]->tree = parent->tree->subtrees[index];

      find_grid_LS_value(parent->subgrids[index]);

      /*
       * Next statements test whether zero LS contour passes through
       * subgrids[index] search grid.
       */

      if (interface_in_grid(parent->subgrids[index]) == TRUE) {
        /*
         * Divide subgrid further if interface is in grid.  Note that if
         * if this subgrid is at max_level, there will be an immediate return
         */

        divide_search_grid(parent->subgrids[index], max_level);
      } else {
        /*
         * This grid doesn't contain the LS interface and therefore is no longer
         * of interest.  So free it and cap its branch with a NULL pointer.
         */

        free_search_grid(&(parent->subgrids[index]));
      }
    }
  }
  return;
}

static void find_grid_LS_value(SGRID *grid) {
  int i, j;
  int num_fcns, num_verts;
  double(*phi)[MDE];

  num_verts = grid->num_verts;
  num_fcns = grid->tree->num_fcns;
  phi = grid->tree->phi;

  for (i = 0; i < num_verts; i++) {
    grid->LS_value[i] = 0.0;

    for (j = 0; j < num_fcns; j++) {
      grid->LS_value[i] += *(esp->F[j]) * phi[i][j];
    }
  }

  return;
}

#if 1
void map_local_coordinates(double *xi, double *x) {
  int a, j;
  int dim = ei[pg->imtrx]->ielem_dim;
  int ShapeVar = pd->ShapeVar;
  int DeformingMesh = pd->e[pg->imtrx][R_MESH1];
  int mdof = ei[pd->mi[ShapeVar]]->dof[ShapeVar];
  int ln, I;
  int iconnect = Proc_Connect_Ptr[ei[pg->imtrx]->ielem];
  double phi_j;

  if (ei[pg->imtrx]->ielem_shape == SHELL ||
      ei[pg->imtrx]->ielem_shape == TRISHELL) {
    dim = pd->Num_Dim;
  }

  for (a = 0; a < dim; a++) {
    x[a] = 0.0;

    if (!DeformingMesh) {
      for (j = 0; j < mdof; j++) {
        ln = ei[pd->mi[ShapeVar]]->dof_list[ShapeVar][j];

        I = Proc_Elem_Connect[iconnect + ln];

        phi_j = newshape(xi, ei[pg->imtrx]->ielem_type, PSI, ln,
                         ei[pg->imtrx]->ielem_shape,
                         pd->i[pd->mi[ShapeVar]][ShapeVar], j);

        x[a] += Coor[a][I] * phi_j;
      }
    } else {
      for (j = 0; j < mdof; j++) {
        ln = ei[pd->mi[ShapeVar]]->dof_list[ShapeVar][j];

        I = Proc_Elem_Connect[iconnect + ln];

        phi_j = newshape(xi, ei[pg->imtrx]->ielem_type, PSI, ln,
                         ei[pg->imtrx]->ielem_shape,
                         pd->i[pd->mi[ShapeVar]][ShapeVar], j);

        x[a] += (Coor[a][I] + *esp->d[a][j]) * phi_j;
      }
    }
  }

  return;
}
#endif

static int interface_in_grid(SGRID *grid) {
  int i;
  int num_verts = grid->num_verts;

  double val0 = grid->LS_value[0];

  if (val0 == 0.0)
    return (TRUE);

  for (i = 1; i < num_verts; i++) {
    if (grid->LS_value[i] * val0 <= 0.0)
      return (TRUE);
  }

  return FALSE;
}

void free_search_grid(SGRID **ptr_grid) {
  int i;
  SGRID *grid = *ptr_grid;

  if (grid == NULL)
    return;

  if (grid->num_subgrids == 0) {
    safe_free(grid);
  } else {
    for (i = 0; i < grid->num_subgrids; i++) {
      free_search_grid(&(grid->subgrids[i]));
    }

    safe_free(grid->subgrids);

    grid->num_subgrids = 0;

    free_search_grid(&grid);
    *ptr_grid = NULL;
  }
  return;
}

void print_search_grid(SGRID *grid) {
  int i;
  int num_verts;
  double x[DIM];

  int l;

  if (grid == NULL)
    return;

  l = grid->level;
  num_verts = (grid->dim == 2) ? 4 : 8;

  if (grid->num_subgrids == 0) {
    for (i = 0; i < num_verts; i++) {
      l = grid->level;

      while (l-- > 0)
        printf("\t");

      map_local_coordinates(grid->tree->xi[i], x);

      printf("%f\t%f\t%f\n", x[0], x[1], x[2]);
    }
  } else {
    for (i = 0; i < grid->num_subgrids; i++) {
      l = grid->level;
      while (l-- > 0)
        printf("\t");
      printf("Level %d, Subgrid %d \n", grid->level, i);

      print_search_grid(grid->subgrids[i]);
    }
  }
  return;
}

void find_grid_intersections(SGRID *grid, struct LS_Surf_List *list) {
  int index;
  int i, j, l, link;
  double xi[3] = {0., 0., 0.}, yi[3] = {0., 0., 0.};
  double x[3];
  struct LS_Surf *surf;

  if (grid == NULL)
    return;

  if (grid->num_subgrids ==
      0) /* Search only the lowest level grids for intersections */
  {
    switch (grid->dim) {
    case 2: {
      static int links[6][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 2}, {1, 3}};

      for (link = 0; link < 6; link++) {
        i = links[link][0];
        j = links[link][1];

        for (l = 0; l < grid->dim; l++) {
          xi[l] = grid->tree->xi[i][l];
          yi[l] = grid->tree->xi[j][l];
        }

        if (find_link_intersection(xi, yi, ls->var, 0.0, NULL) == TRUE) {

          map_local_coordinates(xi, x);
          surf = create_surf_point(x, ei[pg->imtrx]->ielem, xi, FALSE);

          if (unique_surf(list, surf)) {
            append_surf(list, surf);
          } else {
            safe_free(surf);
          }
        }
      }
    } break;
    case 3: {
      static int links[4][2] = {
          {0, 6},
          {1, 7},
          {2, 4},
          {3, 5}}; /* These are the diagonals.  Saves us having to check
                      for dupes. Expensive in 3D */

      for (link = 0; link < 4; link++) {
        i = links[link][0];
        j = links[link][1];

        for (l = 0; l < grid->dim; l++) {
          xi[l] = grid->tree->xi[i][l];
          yi[l] = grid->tree->xi[j][l];
        }

        if (find_link_intersection(xi, yi, ls->var, 0.0, NULL) == TRUE) {
          map_local_coordinates(xi, x);
          surf = create_surf_point(x, ei[pg->imtrx]->ielem, xi, FALSE);

          append_surf(list, surf);
        }
      }
    } break;
    default:
      break;
    }
  } else {
    for (index = 0; index < grid->num_subgrids; index++) {
      find_grid_intersections(grid->subgrids[index], list);
    }
  }

  return;
}

/*
 *   Shape function tree section
 *
 *   I should really discuss this particular structure in more detail.  It
 * became apparent to me when setting up the Search_Grid recursive structures
 * for both interface reconstruction and adaptive subgrid integration, that I
 * was having to compute locations of grid vertices and values of shape
 * functions at these points  and locations of integration points over and over.
 * All of these quantities can be computed up front on the (s,u,t) master
 * element Hence the notion of the shape function tree structure was aborned.
 *
 *   The NTREE structure is similar to the SGRID structure in that it is a
 * self-similar, recursively-connected tree structure.  That is, the parent
 * structure points to children structures of the same type, but the parent
 *   structure might be the child of a larger parent itself.  A given NTREE
 * structure stores its dimensionality, its division depth, the number of
 * vertices it has, the (s,u,t) coordinates of these vertices, the values of a
 * single shape function interpolant at these vertices, the location and weights
 * of numerical integration points on the the structure.  It knows how many
 * children it has and has pointers to all of them.  It, however, doesn't know
 * who its parent is.  Sad really.
 *
 *   The SGRID structure can be be sparesly populated, that is, some branches
 * terminate prematurely before they reach the maximum division depth set
 * (max_level).  This is because some grid's don't overlap the interface or the
 * "interface region." This is not so with the NTREE structure.  All branches
 * are fully populated to the maximum division depth.  This is because the NTREE
 * structure is to be created up front and used as a resource when constructing
 * the SGRID structures.  One can't know *a priori* which grids any SGRID
 * reconstruction will need, so we have to compute ALL of them.  This can mean a
 * lot of information.  For example, an 3D NTREE structure to a division depth
 * of 4 will have 1 + 8 + 64 + 512 + 4096 = 4681 distinct grids within it.  All
 * of them with 8 vertices along with associated shape function values not to
 * mention the integration point information.  A lot of computations and
 * something you don't want to reconstruct for every element.
 */

NTREE *create_shape_fcn_tree(int max_level) {
  NTREE *tree;

  tree = (NTREE *)smalloc(sizeof(NTREE));

  /* tree->dim = pd->Num_Dim;*/
  tree->dim = ei[pg->imtrx]->ielem_dim;
  tree->level = 0;
  tree->num_verts = (tree->dim == 2) ? 4 : 8;
  tree->bf = bf[ls->var];
  tree->num_fcns = ei[pg->imtrx]->dof[ls->var];
  tree->num_subtrees = 0;
  tree->subtrees = NULL;

  switch (tree->dim) {
    int j;
  case 2:
    tree->xi = (double(*)[DIM])smalloc(sizeof(double) * 4 * DIM);
    tree->phi = (double(*)[MDE])smalloc(sizeof(double) * 4 * MDE);

    for (j = 0; j < 4; j++) {
      find_nodal_stu(j, ei[pg->imtrx]->ielem_type, tree->xi[j], tree->xi[j] + 1,
                     tree->xi[j] + 2);
    }

    tree->s = (double(*)[DIM])smalloc(sizeof(double) * 9 * DIM);
    tree->wt = (double *)smalloc(sizeof(double) * 9);

    break;

  case 3:

    tree->xi = (double(*)[DIM])smalloc(sizeof(double) * 8 * DIM);
    tree->phi = (double(*)[MDE])smalloc(sizeof(double) * 8 * MDE);

    for (j = 0; j < 8; j++) {
      find_nodal_stu(j, ei[pg->imtrx]->ielem_type, tree->xi[j], tree->xi[j] + 1,
                     tree->xi[j] + 2);
    }

    tree->s = (double(*)[DIM])smalloc(sizeof(double) * 8 * DIM);
    tree->wt = (double *)smalloc(sizeof(double) * 8);

    break;
  default:
    break;
  }

  compute_shape_fcn_values(tree);

  divide_shape_fcn_tree(tree, max_level);

  return (tree);
}

static void compute_shape_fcn_values(NTREE *tree)

{
  int i, j, ledof, jdof;
  double ri[DIM];

  for (i = 0; i < tree->num_verts; i++) {
    ri[0] = tree->xi[i][0];
    ri[1] = tree->xi[i][1];
    ri[2] = tree->xi[i][2];

    /* load_basis_functions */
    jdof = 0;
    for (j = 0; j < ei[pg->imtrx]->dof[ls->var]; j++) {
      ledof = ei[pg->imtrx]->lvdof_to_ledof[ls->var][j];
      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
        tree->phi[i][j] = newshape(ri, ei[pg->imtrx]->ielem_type, PSI,
                                   ei[pg->imtrx]->dof_list[ls->var][j],
                                   ei[pg->imtrx]->ielem_shape,
                                   pd->i[pg->imtrx][ls->var], jdof);
        jdof++;
      } else {
        tree->phi[i][j] = 0.0;
      }
    }
  }
  return;
}

#define MAX_SUBTREE_NODES 27

static void divide_shape_fcn_tree(NTREE *parent, int max_level) {
  if (parent->level < max_level) {
    int index;
    int num_subtrees;
    double xi_m[DIM] = {0.0, 0.0, 0.0};
    double subtree_xi[MAX_SUBTREE_NODES][DIM];
    int num_gpts = 0;

    switch (parent->dim) {
    case 3:
      xi_m[2] = (parent->xi[0][2] + parent->xi[4][2]) / 2.0;
	  /* fall through */
    case 2:
      xi_m[0] = (parent->xi[0][0] + parent->xi[1][0]) / 2.0;
      xi_m[1] = (parent->xi[1][1] + parent->xi[2][1]) / 2.0;
    }

    num_subtrees = (parent->dim == 2) ? 4 : 8;

    parent->num_subtrees = num_subtrees;

    parent->subtrees = (NTREE **)smalloc(num_subtrees * (int) sizeof(NTREE *));

    gather_subtree_coords(parent, xi_m, subtree_xi);

    for (index = 0; index < num_subtrees; index++) {
      int num_verts = (parent->dim == 2) ? 4 : 8;

      parent->subtrees[index] = (NTREE *)smalloc(sizeof(NTREE));

      parent->subtrees[index]->dim = parent->dim;
      parent->subtrees[index]->level = parent->level + 1;
      parent->subtrees[index]->num_verts = num_verts;
      parent->subtrees[index]->bf = parent->bf;
      parent->subtrees[index]->num_fcns = parent->num_fcns;
      parent->subtrees[index]->num_subtrees = 0;
      parent->subtrees[index]->subtrees = NULL;

      parent->subtrees[index]->xi =
          (double(*)[DIM])smalloc(num_verts * (int) sizeof(double) * DIM);
      parent->subtrees[index]->phi =
          (double(*)[MDE])smalloc(num_verts * (int) sizeof(double) * MDE);

      load_subtree_coords(index, parent->subtrees[index], subtree_xi);

      if (TRUE) {
        if (parent->dim == 3)
          num_gpts = (parent->subtrees[index]->level < 8) ? 8 : 1;
        if (parent->dim == 2)
          num_gpts = (parent->subtrees[index]->level < 8) ? 4 : 1;

        parent->subtrees[index]->s =
            (double(*)[DIM])smalloc(num_gpts * DIM * sizeof(double));
        parent->subtrees[index]->wt =
            (double *)smalloc(num_gpts * sizeof(double));

        find_tree_integration_pts(parent->subtrees[index], num_gpts);
      }

      compute_shape_fcn_values(parent->subtrees[index]);

      divide_shape_fcn_tree(parent->subtrees[index], max_level);
    }
  }
}

static void gather_subtree_coords(NTREE *tree, double *xi_m,
                                  double (*sub_xi)[DIM]) {

  int i, a;
  double s = 0, t = 0, u = 0;

  if (tree == NULL)
    return;

  for (a = 0; a < DIM; a++) {
    for (i = 0; i < tree->num_verts; i++) {
      sub_xi[i][a] = tree->xi[i][a];
    }
  }

  switch (tree->dim) {
  case 2:

    for (i = 4; i < 9; i++) {
      switch (i) {
      case 4:
        s = xi_m[0];
        t = tree->xi[0][1];
        break;

      case 5:
        s = tree->xi[1][0];
        t = xi_m[1];
        break;

      case 6:
        s = xi_m[0];
        t = tree->xi[2][1];
        break;

      case 7:
        s = tree->xi[3][0];
        t = xi_m[1];
        break;

      case 8:
        s = xi_m[0];
        t = xi_m[1];
        break;
      default:
        break;
      }

      sub_xi[i][0] = s;
      sub_xi[i][1] = t;
      sub_xi[i][2] = 0.0;
    }

    break;

  case 3:
    for (i = 8; i < 27; i++) {
      switch (i) {
      case 8:
        s = xi_m[0];
        t = tree->xi[0][1];
        u = tree->xi[0][2];
        break;
      case 9:
        s = tree->xi[1][0];
        t = xi_m[1];
        u = tree->xi[1][2];
        break;

      case 10:
        s = xi_m[0];
        t = tree->xi[2][1];
        u = tree->xi[2][2];
        break;

      case 11:
        s = tree->xi[3][0];
        t = xi_m[1];
        u = tree->xi[3][2];
        break;

      case 12:
        s = tree->xi[0][0];
        t = tree->xi[0][1];
        u = xi_m[2];
        break;

      case 13:
        s = tree->xi[1][0];
        t = tree->xi[1][1];
        u = xi_m[2];
        break;

      case 14:
        s = tree->xi[2][0];
        t = tree->xi[2][1];
        u = xi_m[2];
        break;

      case 15:
        s = tree->xi[3][0];
        t = tree->xi[3][1];
        u = xi_m[2];
        break;

      case 16:
        s = xi_m[0];
        t = tree->xi[4][1];
        u = tree->xi[4][2];
        break;
      case 17:
        s = tree->xi[5][0];
        t = xi_m[1];
        u = tree->xi[5][2];
        break;

      case 18:
        s = xi_m[0];
        t = tree->xi[6][1];
        u = tree->xi[6][2];
        break;

      case 19:
        s = tree->xi[7][0];
        t = xi_m[1];
        u = tree->xi[7][2];
        break;

      case 20:
        s = xi_m[0];
        t = xi_m[1];
        u = xi_m[2];
        break;

      case 21:
        s = xi_m[0];
        t = xi_m[1];
        u = tree->xi[0][2];
        break;

      case 22:
        s = xi_m[0];
        t = xi_m[1];
        u = tree->xi[4][2];
        break;

      case 23:
        s = tree->xi[0][0];
        t = xi_m[1];
        u = xi_m[2];
        break;

      case 24:
        s = tree->xi[1][0];
        t = xi_m[1];
        u = xi_m[2];
        break;

      case 25:
        s = xi_m[0];
        t = tree->xi[0][1];
        u = xi_m[2];
        break;

      case 26:
        s = xi_m[0];
        t = tree->xi[2][1];
        u = xi_m[2];
        break;

      default:
        break;
      }

      sub_xi[i][0] = s;
      sub_xi[i][1] = t;
      sub_xi[i][2] = u;
    }
  default:
    break;
  }
  return;
}

static void load_subtree_coords(int index, NTREE *tree, double (*sub_xi)[DIM]) {
  int node;
  switch (tree->dim) {
  case 2: {
    static int map2D[4][4] = {
        {0, 4, 8, 7}, {4, 1, 5, 8}, {8, 5, 2, 6}, {7, 8, 6, 3}};

    for (node = 0; node < 4; node++) {
      tree->xi[node][0] = sub_xi[map2D[index][node]][0];
      tree->xi[node][1] = sub_xi[map2D[index][node]][1];
      tree->xi[node][2] = sub_xi[map2D[index][node]][2];
    }

  } break;

  case 3: {
    static int map3D[8][8] = {
        {0, 8, 21, 11, 12, 25, 20, 23},  {8, 1, 9, 21, 25, 13, 24, 20},
        {21, 9, 2, 10, 20, 24, 14, 26},  {11, 21, 10, 3, 23, 20, 26, 15},
        {12, 25, 20, 23, 4, 16, 22, 19}, {25, 13, 24, 20, 16, 5, 17, 22},
        {20, 24, 14, 26, 22, 17, 6, 18}, {23, 20, 26, 15, 19, 22, 18, 7}};

    for (node = 0; node < 8; node++) {
      tree->xi[node][0] = sub_xi[map3D[index][node]][0];
      tree->xi[node][1] = sub_xi[map3D[index][node]][1];
      tree->xi[node][2] = sub_xi[map3D[index][node]][2];
    }

  }

  break;
  default:
    break;
  }
  return;
}

void free_shape_fcn_tree(NTREE *tree) {
  int i;

  if (tree == NULL)
    return;

  if (tree->num_subtrees == 0) {
    safe_free(tree->xi);
    safe_free(tree->phi);
    safe_free(tree->s);
    safe_free(tree->wt);
    safe_free(tree);
  } else {
    for (i = 0; i < tree->num_subtrees; i++) {
      free_shape_fcn_tree(tree->subtrees[i]);
    }

    safe_free(tree->subtrees);

    tree->num_subtrees = 0;

    free_shape_fcn_tree(tree);
  }
  return;
}

/*
static void
print_shape_fcn_tree( NTREE *tree)
{
  int i;
  int num_verts;
  double x[DIM] = {0.,0.,0.};

  int l;

  if ( tree == NULL ) return;

  l = tree->level;
  num_verts = ( tree->dim == 2 ) ? 4 : 8;


  if( tree->num_subtrees == 0 )
    {
      for( i = 0 ; i < num_verts; i++ )
        {
          l = tree->level;

          while ( l-- > 0 ) printf("\t");

          map_local_coordinates( tree->xi[i], x );

          printf( "%lf\t%lf\t%lf\n", x[0], x[1], x[2]);
        }
    }
  else
    {
      for( i = 0; i<tree->num_subtrees; i++ )
        {
          l = tree->level;
          while ( l-- > 0 ) printf("\t");
          printf("Level %d, Subtree %d \n",tree->level, i );


          print_shape_fcn_tree( tree->subtrees[i] );
        }
    }
  return;
}
*/

int build_integration_grid(SGRID *parent, int max_level, double width) {
  int total_gpts = 0;

  if (parent->level <
      max_level) /* If my level is the max level, return to end the recursion */
  {
    int index;
    int num_subgrids;
    int dim = parent->dim;

    num_subgrids = (dim == 2) ? 4 : 8;

    parent->num_subgrids = num_subgrids;

    /* allocated subgrid structures  */

    parent->subgrids = (SGRID **)array_alloc(1, num_subgrids, sizeof(SGRID *));

    for (index = 0; index < num_subgrids; index++)

    {

      parent->subgrids[index] = (SGRID *)array_alloc(1, 1, sizeof(SGRID));
      parent->subgrids[index]->ei = parent->ei;
      parent->subgrids[index]->dim = parent->dim;
      parent->subgrids[index]->level = parent->level + 1;
      parent->subgrids[index]->num_verts = (parent->dim == 2) ? 4 : 8;
      parent->subgrids[index]->num_subgrids = 0;
      parent->subgrids[index]->subgrids = NULL;
      parent->subgrids[index]->tree = parent->tree->subtrees[index];

      find_grid_LS_value(parent->subgrids[index]);

      /*
       * Next statement test whetherelement overlaps interface region
       */

      if (grid_overlaps_interface(parent->subgrids[index], width)) {
        /*
         * if so we divide it again.
         */

        total_gpts +=
            build_integration_grid(parent->subgrids[index], max_level, width);
      } else {
        /* Do nothing.  We want to keep the grid to do integration on later
         * but don't need to divide it further.
         */

        total_gpts += parent->subgrids[index]->tree->num_gpts;
      }
    }
  } else {
    total_gpts = parent->tree->num_gpts;
  }
  return (total_gpts);
}

int find_tree_integration_pts(NTREE *tree, int num_gpts) {
  int i, elem_type;
  double l1 = 0, l2 = 0, l3 = 0, s_m[DIM];
  double s, t, u;

  compute_tree_size(tree, &l1, &l2, &l3, s_m);

  switch (tree->dim) {
  case 2:
    switch (num_gpts) {
    case 9:
    case 4:

      elem_type = num_gpts == 4 ? BILINEAR_QUAD : BIQUAD_QUAD;

      for (i = 0; i < num_gpts; i++) {
        find_stu(i, elem_type, &s, &t, &u);

        tree->s[i][0] = s_m[0] + l1 * s / 2.0;
        tree->s[i][1] = s_m[1] + l2 * t / 2.0;
        tree->s[i][2] = 0.0;

        tree->wt[i] = Gq_weight(i, elem_type);
        tree->wt[i] *= fabs(l1 * l2) / 4.0;

        tree->num_gpts = num_gpts;
      }
      break;

    case 1:
      tree->s[0][0] = s_m[0];
      tree->s[0][1] = s_m[1];
      tree->s[0][2] = 0.0;
      tree->wt[0] = fabs(l1 * l2);
      tree->num_gpts = 1;
      break;

    default:
      break;
    }
    break;
  case 3:
    switch (num_gpts) {
    case 27:
    case 8:
      elem_type = num_gpts == 8 ? TRILINEAR_HEX : TRIQUAD_HEX;

      for (i = 0; i < num_gpts; i++) {
        find_stu(i, elem_type, &s, &t, &u);

        tree->s[i][0] = s_m[0] + l1 * s / 2.0;
        tree->s[i][1] = s_m[1] + l2 * t / 2.0;
        tree->s[i][2] = s_m[2] + l3 * u / 2.0;

        tree->wt[i] = Gq_weight(i, elem_type);
        tree->wt[i] *= fabs(l1 * l2 * l3) / 8.0;

        tree->num_gpts = num_gpts;
      }
      break;

    case 1:
      tree->s[0][0] = s_m[0];
      tree->s[0][1] = s_m[1];
      tree->s[0][2] = s_m[2];
      tree->wt[0] = fabs(l1 * l2 * l3);
      tree->num_gpts = 1;
      break;
    default:
      break;
    }
    break;
  default:
    break;
  }
  return (0);
}

void compute_tree_size(NTREE *tree, double *l1, double *l2, double *l3,
                       double *s_m) {
  switch (tree->dim) {
  case 2:
    *l1 = tree->xi[1][0] - tree->xi[0][0];
    *l2 = tree->xi[2][1] - tree->xi[1][1];
    *l3 = 0.0;

    s_m[0] = (tree->xi[1][0] + tree->xi[0][0]) / 2.0;
    s_m[1] = (tree->xi[2][1] + tree->xi[1][1]) / 2.0;
    s_m[2] = 0.0;
    break;
  case 3:
    *l1 = tree->xi[1][0] - tree->xi[0][0];
    *l2 = tree->xi[2][1] - tree->xi[1][1];
    *l3 = tree->xi[4][2] - tree->xi[0][2];

    s_m[0] = (tree->xi[1][0] + tree->xi[0][0]) / 2.0;
    s_m[1] = (tree->xi[2][1] + tree->xi[1][1]) / 2.0;
    s_m[2] = (tree->xi[4][2] + tree->xi[0][2]) / 2.0;
    break;
  default:
    break;
  }
}

int grid_overlaps_interface(SGRID *grid, double width) {
  double interface_low = -width / 2.0;
  double interface_high = width / 2.0;

  int all_low = TRUE;
  int all_high = TRUE;

  int i;

  for (i = 0; i < grid->num_verts; i++) {
    all_low = (grid->LS_value[i] <= interface_low) && all_low;

    all_high = (grid->LS_value[i] >= interface_high) && all_high;
  }

  return (!(all_low || all_high));
}

int gather_integration_pts(SGRID *grid, double (*s)[DIM], double *wt,
                           int num_gpts) {
  int index;

  if (grid->num_subgrids != 0) {
    for (index = 0; index < grid->num_subgrids; index++) {
      num_gpts = gather_integration_pts(grid->subgrids[index], s, wt, num_gpts);
    }
  } else {
    NTREE *tree = grid->tree;

    for (index = 0; index < tree->num_gpts; index++) {
      s[num_gpts][0] = tree->s[index][0];
      s[num_gpts][1] = tree->s[index][1];
      s[num_gpts][2] = tree->s[index][2];

      wt[num_gpts] = tree->wt[index];

      num_gpts++;
    }
  }
  return (num_gpts);
}

int print_subgrid_integration_pts(double (*s)[DIM], double *wt, int num_gpts) {
  int i;
  double x[DIM];

  FILE *ifp;

#define DEBUG_SUB_INTEGRATION 0
#if DEBUG_SUB_INTEGRATION
  double sum = 0.;
  double test = 0.;
  int nx = 2;
  int ny = 2;
  double answer;
  int ip_total = elem_info(NQUAD, ei[pg->imtrx]->ielem_type);
  double xi[3];
  double orig_sum = 0.;
  double orig_test = 0.;
  double orig_wt;
#endif

  ifp = fopen("subgrid_int_pts", "a");

  for (i = 0; i < num_gpts; i++) {
    map_local_coordinates(s[i], x);

    if (pd->Num_Dim == 2)
      fprintf(ifp, "%lf\t%lf\n", x[0], x[1]);
    else if (pd->Num_Dim == 3)
      fprintf(ifp, "%lf\t%lf\t%lf\n", x[0], x[1], x[2]);

#if DEBUG_SUB_INTEGRATION
    test += (pow(x[0], nx) * pow(x[1], ny)) * wt[i];
    sum += wt[i];
#endif
  }

#if DEBUG_SUB_INTEGRATION
  for (i = 0; i < ip_total; i++) {
    find_stu(i, ei[pg->imtrx]->ielem_type, xi, xi + 1,
             xi + 2); /* find quadrature point */
    orig_wt = Gq_weight(
        i, ei[pg->imtrx]->ielem_type); /* find quadrature weights for */
    map_local_coordinates(xi, x);

    orig_test += (pow(x[0], nx) * pow(x[1], ny)) * orig_wt;
    orig_sum += orig_wt;
  }
  if ((sum - orig_sum) / orig_sum > 1.e-8)
    printf("subgrid total wt = %g, orig total wt = %g, err = %g\n", sum,
           orig_sum, (sum - orig_sum) / orig_sum);
  if ((test - orig_test) / orig_test > 1.e-8)
    printf("subgrid test = %g, orig test = %g, err = %g\n", test, orig_test,
           (test - orig_test) / orig_test);
#endif

  fclose(ifp);
  return 0;
}

int print_subgrid_surface_integration_pts(double (*s)[DIM], double *wt,
                                          int num_gpts) {
  int i;
  double x[DIM];

  FILE *ifp;

  ifp = fopen("subgrid_surface_int_pts", "a");

  for (i = 0; i < num_gpts; i++) {
    map_local_coordinates(s[i], x);

    if (pd->Num_Dim == 2)
      fprintf(ifp, "%lf\t%lf\t%lf\t%lf\n", x[0], x[1], x[2], wt[i]);
    else if (pd->Num_Dim == 3)
      fprintf(ifp, "%lf\t%lf\t%lf\t%lf\n", x[0], x[1], x[2], wt[i]);
  }

  fclose(ifp);
  return 0;
}

int current_elem_overlaps_interface(double width) {

  double interface_low = -width / 2.0;
  double interface_high = width / 2.0;
  double value = 0.0;

  int all_low = TRUE;
  int all_high = TRUE;
  int all_positive = TRUE;
  int all_negative = TRUE;

  int element_status = FALSE;

  int i;
  double absolute_minimum = DBL_MAX;
  double absolute_maximum = 0.;

  if (!(pd->gv[ls->var]))
    return (FALSE);

  if (width == 0.) {
    int on_isosurface = current_elem_on_isosurface(ls->var, 0.);
    return on_isosurface;
  }

  for (i = 0; i < ei[pd->mi[ls->var]]->dof[ls->var]; i++) {
    switch (ls->var) {
    case FILL:
      value = *(esp->F[i]);
      break;
    case PHASE1:
    case PHASE2:
    case PHASE3:
    case PHASE4:
    case PHASE5:
      value = *(esp->pF[ls->var - PHASE1][i]);
      break;
    }

    all_low = (value < interface_low) && all_low;

    all_high = (value >= interface_high) && all_high;

    all_positive = (value >= 0.0) && all_positive;

    all_negative = (value < 0.0) && all_negative;

    if (fabs(value) < absolute_minimum) {
      absolute_minimum = fabs(value);
    }
    if (fabs(value) > absolute_maximum) {
      absolute_maximum = fabs(value);
    }
  }

  if (all_low || all_high) {
    /* This is an unequivocal case.  The element does not overlap */

    element_status = FALSE;
  } else if (!all_positive && !all_negative) {
    /* another unequivocal case.  The zero contour is in the element */

    element_status = TRUE;
  } else {
    /* This is the sticky case.  The high or low contour is in the element
     * but is possible that it is just entering or leaving the element.  In
     * which case there is potential that no integration points will fall within
     * the interface band.  This case must be detected and the element must be
     * labelled false.
     */
#if 0
      double distance_scale, band_width;
      
      distance_scale = 2.0*(absolute_maximum - absolute_minimum);

      band_width = distance_scale*0.5*( 1.0/(pow( 2.0, (double) ls->Integration_Depth ) ) );

      if ( absolute_minimum + band_width >= width/2.0 )
	{
	  element_status = FALSE;

/* 	  fprintf(stderr, "Status check for element : %d  is FALSE \n", ei[pg->imtrx]->ielem +1 ); */
	}
      else
	{
	  element_status = TRUE;
	}
#endif
    element_status = TRUE;
  }

  return (element_status);
}

int elem_overlaps_interface(int elem, double x[], const Exo_DB *exo,
                            double width) {
  int i, I;
  int iconn_ptr = exo->elem_ptr[elem];
  int dofs;
  double value;
  int ielem_type, ielem_shape;
  int ebn, mn;
  double interface_low = -width / 2.0;
  double interface_high = width / 2.0;
  int all_low = TRUE;
  int all_high = TRUE;
  int all_positive = TRUE;
  int all_negative = TRUE;

  int element_status = FALSE;

#if 0
  int  min_node;
#endif

  double absolute_minimum = DBL_MAX;
  double absolute_maximum = 0.;

  ebn = find_elemblock_index(elem, exo);
  mn = Matilda[ebn];
  if (!(pd_glob[mn]->v[pg->imtrx][ls->var]))
    return (FALSE);

  if (width == 0.) {
    int on_isosurface = elem_on_isosurface(elem, x, exo, ls->var, 0.);
    return on_isosurface;
  }

  ielem_type = Elem_Type(exo, elem);
  ielem_shape = type2shape(ielem_type);
  dofs = getdofs(ielem_shape, pd_glob[mn]->i[pg->imtrx][ls->var]);

  for (i = 0; i < dofs; i++) {
    I = exo->node_list[iconn_ptr + i];

    value = x[Index_Solution(I, ls->var, 0, 0, -2, pg->imtrx)];

    all_low = (value < interface_low) && all_low;

    all_high = (value >= interface_high) && all_high;

    all_positive = (value >= 0.0) && all_positive;

    all_negative = (value < 0.0) && all_negative;

    if (fabs(value) < absolute_minimum) {
      absolute_minimum = fabs(value);
#if 0
	  min_node = i;
#endif
    }
    if (fabs(value) > absolute_maximum) {
      absolute_maximum = fabs(value);
    }
  }

  if (all_low || all_high) {
    /* This is an unequivocal case.  The element does not overlap */

    element_status = FALSE;
  } else if (!all_positive && !all_negative) {
    /* another unequivocal case.  The zero contour is in the element */

    element_status = TRUE;
  } else {
#if 0 /* Change other #if's containing min_node if using */
      /* This is the sticky case.  The high or low contour is in the element 
       * but is possible that it is just entering or leaving the element.  In which case
       * there is potential that no integration points will fall within the 
       * interface band.  This case must be detected and the element must be labelled
       * false.
       */

      double xi[DIM], det, distance_scale, band_width;
      int err;
            
      
      find_nodal_stu( min_node, ielem_type, xi, xi+1, xi+2 );

      err = load_basis_functions(xi, bfd);
      EH( err, "problem from load_basis_functions");
      
      err = beer_belly();
      EH( err, "beer_belly");

      det = bf[ls->var]->detJ;

      distance_scale = 2.0*pow( det , 1.0/ ( (double) pd->Num_Dim ) );

      band_width = distance_scale*0.5*( 1.0/(pow( 2.0, (double) ls->Integration_Depth ) ) );
      

      if ( absolute_minimum + band_width >= width/2.0 )
	{
	  element_status = FALSE;

/* 	  fprintf(stderr, "Status check for element : %d  is FALSE \n", ei[pg->imtrx]->ielem +1 ); */

	}
      else if ( absolute_minimum + band_width < width/2.0 )
	{
	  element_status = TRUE;

	}
#endif
#if 0
      double distance_scale, band_width;
      
      distance_scale = 2.0*(absolute_maximum - absolute_minimum);

      band_width = distance_scale*0.5*( 1.0/(pow( 2.0, (double) ls->Integration_Depth ) ) );
      
      if ( absolute_minimum + band_width >= width/2.0 )
	{
	  element_status = FALSE;
	}
      else
	{
	  element_status = TRUE;
	}
#endif
    element_status = TRUE;
  }

  return (element_status);
}

int get_subgrid_integration_pts(NTREE *start_tree, SGRID **start_grid,
                                double (**s)[DIM], double **weight,
                                double width) {

  int num_gpts;

  free_search_grid(start_grid);

  *start_grid = create_search_grid(start_tree);

  num_gpts = build_integration_grid(*start_grid, ls->Integration_Depth, width);

  if (*s != NULL) {
    safe_free((void *)*s);
    safe_free((void *)*weight);
  }

  *s = (double(*)[DIM])smalloc(DIM * num_gpts * sizeof(double));

  *weight = (double *)smalloc(num_gpts * sizeof(double));

  gather_integration_pts(*start_grid, *s, *weight, 0);

  return (num_gpts);
}

int gather_surface_subgrid_integration_pts(SGRID *grid, int id_side,
                                           double surface_centroid[DIM],
                                           double (*s)[DIM], double *wt,
                                           int num_gpts) {
  int index;

  if (grid->num_subgrids != 0) {
    for (index = 0; index < grid->num_subgrids; index++) {
      num_gpts = gather_surface_subgrid_integration_pts(
          grid->subgrids[index], id_side, surface_centroid, s, wt, num_gpts);
    }
  } else {
    int vert_index, side_index;
    int all_done = FALSE;
    int elem_type, num_surf_pts;
    int ip;
    NTREE *tree = grid->tree;
    double l1 = 0, l2 = 0, l3 = 0, s_m[DIM];
    double xi[DIM] = {0., 0., 0.};
    double ss, tt, uu;

    elem_type = tree->dim == 2 ? BILINEAR_QUAD : TRILINEAR_HEX;

    compute_tree_size(tree, &l1, &l2, &l3, s_m);

    for (vert_index = 0; !all_done && vert_index < tree->num_verts;
         vert_index++) {
      if (vertex_on_element_boundary(tree->xi[vert_index], surface_centroid,
                                     &side_index)) {
        s_m[side_index] = surface_centroid[side_index];

        num_surf_pts = elem_info(NQUAD_SURF, elem_type);

        for (ip = 0; ip < num_surf_pts; ip++) {
          find_surf_st(ip, elem_type, id_side, tree->dim, xi, &ss, &tt, &uu);

          xi[side_index] = 0.0;

          s[num_gpts][0] = s_m[0] + l1 * xi[0] / 2.0;
          s[num_gpts][1] = s_m[1] + l2 * xi[1] / 2.0;
          s[num_gpts][2] = s_m[2] + l3 * xi[2] / 2.0;

          if (tree->dim == 2) {

            wt[num_gpts] =
                Gq_surf_weight(ip, elem_type) * (l1 + l2) / (2.0) / 2.0;
          } else if (tree->dim == 3) {
            wt[num_gpts] = Gq_surf_weight(ip, elem_type) *
                           pow(((l1 + l2 + l3) / 3.0), 2.0) / 4.0;
          }

          num_gpts++;
        }

        all_done = TRUE;
      }
    }
  }

#if 0
			if( vertex_on_element_boundary( tree->xi[vert_index], &side_index ) &&
				fabs( tree->xi[vert_index][side_index] -  surface_centroid[side_index] ) < 1.e-12  )
			{
				
				s[num_gpts][0] = tree->xi[vert_index][0];
				s[num_gpts][1] = tree->xi[vert_index][1];
				s[num_gpts][2] = tree->xi[vert_index][2];
				
				compute_tree_size( tree, &l1, &l2, &l3, s_m );
				
				wt[num_gpts] = pow( 2.0, -(tree->dim-1) )* ( l1 + l2 + l3 ) /( (double) tree->dim );
				
				num_gpts++;
			}
#endif

  return (num_gpts);
}

static int vertex_on_element_boundary(double xi[DIM], double *surface_centroid,
                                      int *id) {
  int i;

  for (i = 0; i < DIM; i++) {
    if (fabs(1.0 - fabs(xi[i])) < 1.e-12) {
      if (fabs(xi[i] - surface_centroid[i]) < 1.e-12) {
        *id = i;
        return (TRUE);
      }
    }
  }

  return (FALSE);
}

void assemble_boundary_extension_velocity(double x[], Exo_DB *exo, Dpi *dpi) {

  /*
   * Some local variables for convenience...
   */

  int eqn;
  int dim;

  int i, j;

  int var, peqn, pvar;

  double *val_ptr;

  double xi[DIM] = {0., 0., 0.};
  int exterior_faces[8];
  int num_exterior_faces;
  int id_side;
  int node_to_set[MDE];

  int local_elem_node_id[MAX_NODES_PER_SIDE];
  int nodes_per_side;
  int have_nodes_to_set = FALSE;
  double sign;
  int debug_here = FALSE;
  int closest_node = -1;
  double F[MDE];
  double closest_F;
  int n, a;
  double F_value;

  /*
   * Unpack variables from structures for local convenience...
   */

  dim = pd->Num_Dim;

  eqn = R_EXT_VELOCITY;
  peqn = upd->ep[pg->imtrx][eqn];

  /* this element needs to be set if it does not span the interface
     and there are nodes that lie on an exterior face that has
     inward-facing characteristics
   */
  if (ls->elem_overlap_state)
    return;

  /* check to see if this element has external faces */
  if ((num_exterior_faces = get_exterior_faces(ei[pg->imtrx]->ielem,
                                               exterior_faces, exo, dpi)) == 0)
    return;

  /* check each exterior face for inward facing characteristics */
  for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++)
    node_to_set[i] = FALSE;
  for (j = 0; j < num_exterior_faces; j++) {
    id_side = exterior_faces[j] + 1;

    get_side_info(ei[pg->imtrx]->ielem_type, id_side, &nodes_per_side,
                  local_elem_node_id);

    /* loop over nodes on side */
    for (n = 0; n < nodes_per_side; n++) {
      /* check if characteristic points inward here */
      i = local_elem_node_id[n];
      find_nodal_stu(i, ei[pg->imtrx]->ielem_type, xi, xi + 1, xi + 2);

      setup_shop_at_point(ei[pg->imtrx]->ielem, xi, exo);
      load_lsi(ls->Length_Scale);
      if (pfd == NULL) {
        F_value = fv->F;
      } else {
        F_value = fv->pF[ls->var - PHASE1];
      }
      F[i] = F_value;

      surface_determinant_and_normal(
          ei[pg->imtrx]->ielem, ei[pg->imtrx]->iconnect_ptr,
          ei[pg->imtrx]->num_local_nodes, ei[pg->imtrx]->ielem_dim - 1, id_side,
          nodes_per_side, local_elem_node_id);

      sign = 0.;
      for (a = 0; a < dim; a++) {
        if (F_value < 0.)
          sign -= fv->snormal[a] * lsi->normal[a];
        else
          sign += fv->snormal[a] * lsi->normal[a];
      }

      if (sign < 0.) {
        node_to_set[i] = TRUE;
        have_nodes_to_set = TRUE;
      }
    }
  }

  if (have_nodes_to_set) {

    if (debug_here)
      printf("Current elem #%d has external nodes to constrain.\n",
             ei[pg->imtrx]->ielem + 1);

    /* find closest node on exterior face */
    closest_F = 1.e30;
    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      if (node_to_set[i]) {
        if (fabs(F[i]) < closest_F) {
          closest_node = i;
          closest_F = fabs(F[i]);
        }
      }
    }

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      if (node_to_set[i] && i != closest_node) {
        if (debug_here) {
          int I = ei[pg->imtrx]->gnn_list[eqn][i];
          printf(
              " will setup node #%d with boundary upwind in current elem #%d\n",
              I + 1, ei[pg->imtrx]->ielem + 1);
          printf(" (closest node in this element is #%d)\n",
                 ei[pg->imtrx]->gnn_list[eqn][closest_node] + 1);
        }

        val_ptr =
            x +
            ei[pg->imtrx]->ieqn_ledof[ei[pg->imtrx]->lvdof_to_ledof[eqn][i]];

        lec->R[peqn][i] += BIG_PENALTY * *val_ptr;

        val_ptr =
            x +
            ei[pg->imtrx]
                ->ieqn_ledof[ei[pg->imtrx]->lvdof_to_ledof[eqn][closest_node]];

        lec->R[peqn][i] -= BIG_PENALTY * *val_ptr;

        if (af->Assemble_Jacobian) {
          /*
           * J_ext_v_ext_v
           */
          var = eqn;
          pvar = upd->vp[pg->imtrx][var];

          lec->J[peqn][pvar][i][i] += BIG_PENALTY;
          lec->J[peqn][pvar][i][closest_node] -= BIG_PENALTY;
        }
      }
    }
  }

  return;
}

int assemble_extension_velocity(dbl hsquared[DIM], dbl hh[DIM][DIM],
                                dbl dh_dxnode[DIM][MDE]) {

  int status;

  /*
   * Some local variables for convenience...
   */

  int eqn;
  int dim;
  int a;

  int i, j;

  dbl *grad_F; /* gradient of Fill. */

  dbl S = 0.0; /* sign of distance function (-1 or +1 ) */

  dbl hh_siz;
  dbl h_elem;

  dbl resid;

  /*
   * Galerkin weighting functions for i-th energy residuals
   * and some of their derivatives...
   */

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl *grad_phi_i, *grad_phi_j;
  dbl h3; /* Volume element (scale factors). */
  dbl det_J;
  dbl wt;
  dbl wt_func, d_wt_func;
  dbl tau;

  dbl advection, advection_a, advection_b;
  dbl gradF_gradphi[MDE];        /* gradF.gradphi */
  dbl n_dot_gradphi[MDE];        /* n.gradphi */
  dbl gradext_v_gradphi[MDE];    /* gradext_v.gradphi */
  dbl gradphi_gradphi[MDE][MDE]; /* gradphi.gradphi */

  dbl d_det_J_dmesh_pj, dh3dmesh_pj;

  dbl *grad_ext_v;

  int var, peqn, pvar, p;

  status = 0;

  /*
   * Unpack variables from structures for local convenience...
   */

  dim = pd->Num_Dim;

  eqn = R_EXT_VELOCITY;

  /*
   * Bail out fast if there's nothing to do...
   */

  if (!pd->e[pg->imtrx][eqn]) {
    return (status);
  }

  wt = fv->wt; /* Gauss point weight. */

  h3 = fv->h3; /* Differential volume element. */

  det_J = bf[eqn]->detJ;

  load_lsi(0.);
  load_lsi_derivs();

  if (pfd == NULL) {
    grad_F = fv->grad_F;
  } else {
    grad_F = fv->grad_pF[ls->var - PHASE1];
  }

  hh_siz = 0.;
  for (a = 0; a < dim; a++) {
    hh_siz += hsquared[a];
  }
  /* This is the average value of h**2 in the element */

  hh_siz = hh_siz / ((double)dim);

  /* This is the size of the element */

  h_elem = sqrt(hh_siz);

  /*grad_F = fv->grad_F;  */
  grad_ext_v = fv->grad_ext_v;
#define GRADF_GRADEXTV 1
/* DRN: although the option NORMAL_GRADEXTV seems more logical,
   the normal is discontinuous so it can show pathological
   behavior that GRADF_GRADEXTV may not suffer from
 */
#define NORMAL_GRADEXTV 0
#if GRADF_GRADEXTV
  resid = 0.;
  for (a = 0; a < VIM; a++)
    resid += grad_F[a] * grad_ext_v[a];
#endif
  for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
    gradF_gradphi[i] = 0.;
    grad_phi_i = bf[eqn]->grad_phi[i];
    for (a = 0; a < VIM; a++) {
      gradF_gradphi[i] += grad_F[a] * grad_phi_i[a];
    }
  }
#if NORMAL_GRADEXTV
  resid = 0.;
  for (a = 0; a < VIM; a++)
    resid += lsi->normal[a] * grad_ext_v[a];
#endif
  for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
    n_dot_gradphi[i] = 0.;
    grad_phi_i = bf[eqn]->grad_phi[i];
    for (a = 0; a < VIM; a++) {
      n_dot_gradphi[i] += lsi->normal[a] * grad_phi_i[a];
    }
  }

  tau = 0.5 * h_elem;

  for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
    grad_phi_i = bf[eqn]->grad_phi[i];

    gradext_v_gradphi[i] = 0.;
    for (a = 0; a < VIM; a++) {
      gradext_v_gradphi[i] += grad_ext_v[a] * grad_phi_i[a];
    }
  }

  var = EXT_VELOCITY;
  for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
    grad_phi_i = bf[eqn]->grad_phi[i];

    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      grad_phi_j = bf[var]->grad_phi[j];

      gradphi_gradphi[i][j] = 0.;
      for (a = 0; a < VIM; a++) {
        gradphi_gradphi[i][j] += grad_phi_j[a] * grad_phi_i[a];
      }
    }
  }

  /* sign of LS */
  S = 2. * lsi->H - 1.;
  if (af->Assemble_Residual) {
    peqn = upd->ep[pg->imtrx][eqn];
    /*
     * In the element, there will be contributions to this many equations
     * based on the number of degrees of freedom...
     */

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

#if GRADF_GRADEXTV
      wt_func =
          S * bf[eqn]->phi[i] + tau * gradF_gradphi[i]; /* Petrov-Galerkin */
#endif
#if NORMAL_GRADEXTV
      wt_func =
          S * bf[eqn]->phi[i] + tau * n_dot_gradphi[i]; /* Petrov-Galerkin */
#endif
      grad_phi_i = bf[eqn]->grad_phi[i];

      advection = 0.;
      if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
        advection = resid * wt_func;
        advection *= det_J * wt * h3;
        advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
      }

      lec->R[peqn][i] += advection;
    }
  }

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
#if GRADF_GRADEXTV
      wt_func =
          S * bf[eqn]->phi[i] + tau * gradF_gradphi[i]; /* Petrov-Galerkin */
#endif
#if NORMAL_GRADEXTV
      wt_func =
          S * bf[eqn]->phi[i] + tau * n_dot_gradphi[i]; /* Petrov-Galerkin */
#endif
      grad_phi_i = bf[eqn]->grad_phi[i];

      /*
       * J_ext_v_ext_v
       */
      var = EXT_VELOCITY;
      if (pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          grad_phi_j = bf[var]->grad_phi[j];

          advection = 0.;

          if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
#if GRADF_GRADEXTV
            for (a = 0; a < dim; a++) {
              advection += grad_F[a] * grad_phi_j[a] * wt_func;
            }
#endif
#if NORMAL_GRADEXTV
            for (a = 0; a < dim; a++) {
              advection += lsi->normal[a] * grad_phi_j[a] * wt_func;
            }
#endif
            advection *= det_J * wt * h3;
            advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
          }

          lec->J[peqn][pvar][i][j] += advection;
        }
      }

#ifdef COUPLED_FILL /* Only add Jacobian entries for coupled fill problems. */
      /*
       * J_ext_v_F
       */

      var = ls->var;

      if (pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          grad_phi_j = bf[var]->grad_phi[j];

          d_wt_func = 0.;
#if GRADF_GRADEXTV
          for (a = 0; a < dim; a++) {
            d_wt_func += tau * grad_phi_j[a] * grad_phi_i[a];
          }
#endif
#if NORMAL_GRADEXTV
          for (a = 0; a < dim; a++) {
            d_wt_func += tau * lsi->d_normal_dF[a][j] * grad_phi_i[a];
          }
#endif

          advection = 0.;

          if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
#if GRADF_GRADEXTV
            for (a = 0; a < dim; a++) {
              advection += grad_phi_j[a] * grad_ext_v[a] * wt_func;
            }
#endif
#if NORMAL_GRADEXTV
            for (a = 0; a < dim; a++) {
              advection += lsi->d_normal_dF[a][j] * grad_ext_v[a] * wt_func;
            }
#endif
            advection += resid * d_wt_func;
            advection *= det_J * wt * h3;
            advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
          }

          lec->J[peqn][pvar][i][j] += advection;
        }
      }
#endif /* not COUPLED_FILL */

      /*
       * J_ext_v_dmesh THIS IS A MESS
       */
      for (p = 0; p < dim; p++) {
        var = MESH_DISPLACEMENT1 + p;
        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_det_J_dmesh_pj = bf[eqn]->d_det_J_dm[p][j];

            dh3dmesh_pj = fv->dh3dq[p] * bf[var]->phi[j];

            advection = 0.;

            if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
              /*
               * three parts:
               *    advection_a =
               *    	Int ( d(Vv)/dmesh h3 |Jv| )
               *
               *    advection_b =
               *  (i)	Int ( Vv h3 d(|Jv|)/dmesh )
               *  (ii)      Int ( Vv dh3/dmesh |Jv|   )
               */

              advection_a = S * resid;

              advection_a *= (d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj);

              advection_b = 0.;
              for (a = 0; a < dim; a++) {
#ifdef COUPLED_FILL

                advection_b += fv->d_grad_F_dmesh[a][p][j] * grad_ext_v[a];
#endif /* COUPLED_FILL */
                advection_b += grad_F[a] * fv->d_grad_ext_v_dmesh[a][p][j];
              }
              advection_b *= det_J * h3;

              advection = advection_a + advection_b;

              advection *=
                  wt_func * wt * pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
            }

            lec->J[peqn][pvar][i][j] += advection;
          }
        }
      }
    }
  }

  status = 1;

  return (status);
}

void get_subelement_descriptions(double x[], Exo_DB *exo,
                                 Integ_Elem_Desc_List *list) {
  /* gather subelement information into a single list */
  int ebi, ielem, e_start, e_end;
  Integ_Elem *e;
  double *x_current = x_static;

  list->size = 0;
  list->start = NULL;

  for (ebi = 0; ebi < exo->num_elem_blocks; ebi++) {

    pd = pd_glob[Matilda[ebi]];
    mp = mp_glob[Matilda[ebi]];

    e_start = exo->eb_ptr[ebi];
    e_end = exo->eb_ptr[ebi + 1];

    if (pd->v[pg->imtrx][ls->var]) {
      for (ielem = e_start; ielem < e_end; ielem++) {

        if (elem_on_isosurface(ielem, x, exo, ls->var, 0.)) {

          load_elem_dofptr(ielem, exo, x, x_old_static, xdot_static,
                           xdot_old_static, 0);

          e = create_integ_elements(0.);

          gather_subelement_descriptions(list, e);

          free_integ_elements(e);
        }
      }
    }
  }

  /* restore static pointer */
  x_static = x_current;

  return;
}

void ls_surface_extents(double x[], Exo_DB *exo, double min[3], double max[3]) {
  int iside;
  Integ_Elem_Desc *s;
  Integ_Elem_Desc_List list;
  int nodes_per_side;
  int lnn[MAX_NODES_PER_SIDE];
  int n;
  int first = TRUE;

  /* first step is to gather subelement information into a single list */
  get_subelement_descriptions(x, exo, &list);

  if (list.size == 0 && Num_Proc == 1) {
    DPRINTF(stderr, "No elements found on the interface.\n");
    return;
  }

  /* now that the information has been gathered, we can examine list */
  s = list.start;
  while (s != NULL) {
    for (iside = 0; iside < s->num_sides; iside++) {
      if (s->bc_sides[iside]) /* side is on ls interface */
      {
        if (s->num_nodes == 6) /* only quadratic triangles supported for now */
        {
          get_side_info(QUAD_TRI, iside + 1, &nodes_per_side, lnn);

          if (first) {
            min[0] = s->x[lnn[0]][0];
            max[0] = s->x[lnn[0]][0];
            min[1] = s->x[lnn[0]][1];
            max[1] = s->x[lnn[0]][1];
            min[2] = 0.;
            max[2] = 0.;
            first = FALSE;
          }

          for (n = 0; n < 3; n++) {
            if (s->x[lnn[n]][0] < min[0])
              min[0] = s->x[lnn[n]][0];
            if (s->x[lnn[n]][0] > max[0])
              max[0] = s->x[lnn[n]][0];
            if (s->x[lnn[n]][1] < min[1])
              min[1] = s->x[lnn[n]][1];
            if (s->x[lnn[n]][1] > max[1])
              max[1] = s->x[lnn[n]][1];
          }
        }
      }
    }
    s = s->next;
  }

  /* be sure to leave the place as clean as you found it */
  free_subelement_descriptions(&list.start);

  return;
}

void free_subelement_descriptions(Integ_Elem_Desc **s_ptr) {
  if (s_ptr != NULL && *s_ptr != NULL) {
    Integ_Elem_Desc *s = *s_ptr;
    free_subelement_descriptions(&(s->next));
    safe_free(s->x);
    safe_free(s);
  }
}

void gather_subelement_descriptions(Integ_Elem_Desc_List *list, Integ_Elem *e) {
  int i;
  if (e->num_subelements > 0) {
    for (i = 0; i < e->num_subelements; i++) {
      gather_subelement_descriptions(list, e->subelements[i]);
    }
  } else {
    Integ_Elem_Desc *s;

    s = (Integ_Elem_Desc *)smalloc(sizeof(Integ_Elem_Desc));
    s->num_nodes = e->num_local_nodes;
    s->x = (double(*)[DIM])smalloc(sizeof(double) * s->num_nodes * DIM);
    for (i = 0; i < s->num_nodes; i++) {
      map_local_coordinates(e->xi[i], s->x[i]);
    }
    s->num_sides = e->num_sides;
    s->bc_sides = (int *)smalloc(sizeof(int) * s->num_sides);
    for (i = 0; i < s->num_sides; i++) {
      s->bc_sides[i] = e->bc_sides[i];
    }
    s->next = list->start;
    list->start = s;
    list->size++;
  }
}

void subelement_mesh_output(double x[], Exo_DB *exo) {
  Integ_Elem_Desc *s;
  Integ_Elem_Desc_List list;
  int i, j, ivconn, isconn;
  int *vconn, *sconn;
  int comp_ws = 8, io_ws = 8;
  int exoid, nnodes, nvelems, nselems, nodes_per_elem, nodes_per_side;
  double *coord_x, *coord_y, *coord_z;
  int iside;
  char filename[256];
  char *description = "Subelement decomposition";
  static int count = 0;

  sprintf(filename, "subelements.%.5d.exoII", count);

  /* first step is to gather subelement information into a single list */
  get_subelement_descriptions(x, exo, &list);

  if (list.size == 0 && Num_Proc == 1) {
    DPRINTF(stderr, "No elements found on the interface.\n");
    return;
  }

  DPRINTF(stderr, "Subelement decomposition contains %d elements.\n",
          list.size);

  /* count sublements and their nodes */
  s = list.start;
  nnodes = 0;
  nvelems = 0;
  nselems = 0;
  nodes_per_elem = 0;
  nodes_per_side = 0;

  /* only support output of linear or quadratic triangles for now */
  while (s != NULL && s->num_nodes != 3 && s->num_nodes != 6)
    s = s->next;

  if (s != NULL) {
    nodes_per_elem = s->num_nodes;
    if (nodes_per_elem == 3)
      nodes_per_side = 2;
    else if (nodes_per_elem == 6)
      nodes_per_side = 3;
    else
      EH(-1, "Subelement type not supported.");
  }

  while (s != NULL) {
    if (s->num_nodes == nodes_per_elem) {
      nnodes += s->num_nodes;
      nvelems += 1;
    }
    for (iside = 0; iside < s->num_sides; iside++) {
      if (s->bc_sides[iside])
        nselems += 1;
    }
    s = s->next;
  }

  if (nvelems == 0) {
    DPRINTF(
        stderr,
        "No valid triangles in sublement list, no output will be written.\n");
    free_subelement_descriptions(&list.start);
    return;
  }

  /* no attempt to form true connectivity */
  DPRINTF(stderr, "Creating sublement file %s for %d nodes and %d elements.\n",
          filename, nnodes, nvelems);

  exoid = ex_create(filename, EX_CLOBBER, &comp_ws, &io_ws);
  ex_put_init(exoid, description, 2, nnodes, nvelems + nselems, 2, 2, 0);

  /* gather coordinates into contiguous arrays and form connectivity array */
  coord_x = (double *)smalloc(nnodes * sizeof(double));
  coord_y = (double *)smalloc(nnodes * sizeof(double));
  coord_z = (double *)smalloc(nnodes * sizeof(double));
  vconn = (int *)smalloc(nvelems * nodes_per_elem * sizeof(int));
  sconn = (int *)smalloc(nselems * nodes_per_side * sizeof(int));
  s = list.start;
  i = 0;
  ivconn = 0;
  isconn = 0;
  while (s != NULL) {
    if (s->num_nodes == nodes_per_elem) {
      for (j = 0; j < s->num_nodes; j++) {
        coord_x[i] = s->x[j][0];
        coord_y[i] = s->x[j][1];
        coord_z[i] = s->x[j][2];
        if (j < nodes_per_elem)
          vconn[ivconn++] = i + 1;
        i++;
      }
      for (iside = 0; iside < s->num_sides; iside++) {
        if (s->bc_sides[iside]) {
          int istart = i - s->num_nodes;
          if (iside == 0) {
            sconn[isconn++] = istart + 0 + 1;
            sconn[isconn++] = istart + 1 + 1;
            if (nodes_per_side == 3)
              sconn[isconn++] = istart + 3 + 1;
          } else if (iside == 1) {
            sconn[isconn++] = istart + 1 + 1;
            sconn[isconn++] = istart + 2 + 1;
            if (nodes_per_side == 3)
              sconn[isconn++] = istart + 4 + 1;
          } else if (iside == 2) {
            sconn[isconn++] = istart + 2 + 1;
            sconn[isconn++] = istart + 0 + 1;
            if (nodes_per_side == 3)
              sconn[isconn++] = istart + 5 + 1;
          } else
            EH(-1, "Subelement type not supported.");
        }
      }
    }
    s = s->next;
  }
  ex_put_coord(exoid, coord_x, coord_y, coord_z);

  /*
  ex_put_node_num_map( exoid, nmap );
  ex_put_elem_num_map( exoid, emap );
  */
  ex_put_block(exoid, EX_ELEM_BLOCK, 1, "TRI", nvelems, nodes_per_elem, 0, 0,
               0);
  ex_put_conn(exoid, EX_ELEM_BLOCK, 1, vconn, 0, 0);

  ex_put_block ( exoid, EX_ELEM_BLOCK, 1, "TRI", nvelems, nodes_per_elem, 0, 0, 0 );
  ex_put_conn ( exoid, EX_ELEM_BLOCK, 1, vconn, 0, 0 );
  ex_put_conn(exoid, EX_ELEM_BLOCK, 2, sconn, 0, 0);
  ex_put_block ( exoid, EX_ELEM_BLOCK, 2, "SHEL", nselems, nodes_per_side, 0, 0, 0 );
  ex_put_conn ( exoid, EX_ELEM_BLOCK, 2, sconn, 0, 0 );

  /* dummy nodeset on volume */
  ex_put_set_param( exoid, EX_NODE_SET, 1, ivconn, 0 );
  ex_put_set( exoid, EX_NODE_SET, 1, vconn, NULL );

  /* dummy nodeset on surface */
  ex_put_set_param( exoid, EX_NODE_SET, 2, isconn, 0 );
  ex_put_set( exoid, EX_NODE_SET, 2, sconn, NULL );

  /* no data for now */

  ex_close(exoid);

  free_subelement_descriptions(&list.start);
  safe_free(coord_x);
  safe_free(coord_y);
  safe_free(coord_z);
  safe_free(vconn);
  safe_free(sconn);

  count++;

  return;
}

int get_facet_integration_pts(double (**s)[DIM], double **weight, Exo_DB *exo) {
  struct LS_Surf_List *list;
  struct LS_Surf *surf;
  struct LS_Surf_Point_Data *vert[2];
  struct LS_Surf_Facet_Data *facet_data;
  double t;
  double length, sdet, wt;
  int ip, ip_total, i, num_gpts;

  if (*s != NULL) {
    safe_free((void *)*s);
    safe_free((void *)*weight);
  }

  /* get list of facets in this element */
  list = create_surf_list();

  find_facets(list, ls->var, 0., exo);

  if (Do_Overlap && ls->CrossMeshQuadPoints > 0)
    ip_total = ls->CrossMeshQuadPoints;
  else
    ip_total = 3;

  num_gpts = ip_total * list->size;

  surf = list->start;

  *s = (double(*)[DIM])smalloc(DIM * num_gpts * sizeof(double));
  *weight = (double *)smalloc(num_gpts * sizeof(double));

  /* process each facet */
  i = 0;
  while (surf) {
    facet_data = (struct LS_Surf_Facet_Data *)surf->data;

    if (facet_data->num_points != 2)
      EH(-1, "Only 2 point facets currently supported");

    vert[0] = (struct LS_Surf_Point_Data *)surf->subsurf_list->start->data;
    vert[1] =
        (struct LS_Surf_Point_Data *)surf->subsurf_list->start->next->data;

    length = sqrt(pow((vert[0]->x[0] - vert[1]->x[0]), 2.0) +
                  pow((vert[0]->x[1] - vert[1]->x[1]), 2.0));

    /* loop over gauss points on facet */
    for (ip = 0; ip < ip_total; ip++) {
      find_segment_s_wt(ip, ip_total, &t, &wt);
      (*s)[i][0] =
          vert[0]->xi[0] + 0.5 * (t + 1.) * (vert[1]->xi[0] - vert[0]->xi[0]);
      (*s)[i][1] =
          vert[0]->xi[1] + 0.5 * (t + 1.) * (vert[1]->xi[1] - vert[0]->xi[1]);
      (*s)[i][2] = 0.;

      map_local_coordinates((*s)[i], fv->x);
      load_coordinate_scales(pd->CoordinateSystem, fv);

      /* Line integral term for detJ */
      /* note this isn't quite right for spherical coords */
      sdet = 0.5 * length * fv->h3;

      (*weight)[i] = wt * sdet;
      i++;
    }

    surf = surf->next;
  } /* while (surf) */

  free_surf_list(&list);

  return num_gpts;
}

static int fail_courant_condition(void) {
  int i;
  double max_F;
  double min_F;
  double dF, max_dF = 0.;
  double **F = esp->F;
  double **F_old = esp_old->F;

  if (TimeIntegration == STEADY)
    return FALSE;

  max_F = *F_old[0];
  min_F = *F_old[0];
  for (i = 0; i < ei[pg->imtrx]->dof[FILL]; i++) {
    if (*F_old[i] > max_F)
      max_F = *F_old[i];
    if (*F_old[i] < min_F)
      min_F = *F_old[i];
    dF = fabs(*F[i] - *F_old[i]);
    if (dF > max_dF)
      max_dF = dF;
  }
  if (max_dF > 10. * (max_F - min_F)) {
    fprintf(stderr, "Courant condition failed in elem=%d\n",
            ei[pg->imtrx]->ielem + 1);
    return TRUE;
  }

  return FALSE;
}

#if 0
double
Courant_Time_Step( double x[], double x_old[], double x_older[], 
                   double xdot[], double xdot_old[], 
                   double resid_vector[],
                   Exo_DB *exo )
{
  int ebi, ielem, e_start, e_end;
  int count = 0;
  double dt, min_dt;
  double v_mag2;
  double hsquared[DIM];
  double hhv[DIM][DIM];
  double dhv_dxnode[DIM][MDE];
  double h_elem;
  int dim, wim;
  int a, i;
  
  for ( ebi=0; ebi<exo->num_elem_blocks; ebi++)
    {

      pd  = pd_glob[Matilda[ebi]];
      mp  = mp_glob[Matilda[ebi]];

      e_start = exo->eb_ptr[ebi];
      e_end   = exo->eb_ptr[ebi+1];
      
      dim = pd->Num_Dim;
      wim = dim;

      if (pd->CoordinateSystem == SWIRLING ||
          pd->CoordinateSystem == PROJECTED_CARTESIAN ||
          pd->CoordinateSystem == CARTESIAN_2pt5D)
        wim = wim+1;

      if (ls->var != NULL)
	{
	  if ( pd->v[pg->imtrx][ls->var] )
	    {
	      for( ielem = e_start ; ielem < e_end ; ielem++)
		{

		  load_elem_dofptr(ielem, exo, x, x_old,
				   xdot, xdot_old,
				   resid_vector, 0);
                               
		  h_elem_siz(hsquared, hhv, dhv_dxnode, pd->e[pg->imtrx][R_MESH1]);
              
		  h_elem = 0.;
		  for ( a=0; a<dim; a++) h_elem += hsquared[a];
		  /* This is the size of the element */
		  h_elem = sqrt(h_elem/ ((double )dim));
              
		  if ( pd->v[pg->imtrx][EXT_VELOCITY] )
		    {
		      for ( i=0; i< ei[pg->imtrx]->dof[EXT_VELOCITY]; i++ )
			{
			  if ( *esp->ext_v[i] != 0. )
			    {
			      dt = fabs(h_elem / *esp->ext_v[i]);
			      if ( count++ == 0 || dt < min_dt ) min_dt = dt;
			    }
			}
		    }
		  if ( pd->v[pg->imtrx][VELOCITY1] && tran->Fill_Equation != FILL_EQN_EXT_V )
		    {
		      for ( i=0; i< ei[pg->imtrx]->dof[VELOCITY1]; i++ )
			{
			  if ( is_xfem_interp( pd->i[pg->imtrx][VELOCITY1] ) )
			    {
			      int xfem_active, extended_dof, base_interp, base_dof;
			      xfem_dof_state( i, pd->i[pg->imtrx][VELOCITY1], ei[pg->imtrx]->ielem_shape,
					      &xfem_active, &extended_dof, &base_interp, &base_dof );
			      if ( extended_dof ) continue;
			    }
			  v_mag2 = 0.;
			  for ( a=0; a<wim; a++ )
			    {
			      v_mag2 += *esp->v[a][i] * *esp->v[a][i];
			    }
			  if ( v_mag2 != 0. )
			    {
			      dt = fabs(h_elem / sqrt( v_mag2 ));
			      if ( count++ == 0 || dt < min_dt ) min_dt = dt;
			    }
			}
		    }
		}
	    }
	}
    }
  return (min_dt);
}
#else
double Courant_Time_Step(double x[], double x_old[], double x_older[],
                         double xdot[], double xdot_old[],
                         double resid_vector[], int *proc_config, Exo_DB *exo) {
  int ebi, ielem, e_start, e_end;
  int count = 0;
  double dt, min_dt = 0.;
  double hsquared[DIM];
  double hhv[DIM][DIM];
  double dhv_dxnode[DIM][MDE];
  double h_elem;
  double *xi, xi_array[3];
  double sum, sumv;
  int dim, wim;
  int ip_total, ip;
  double wt, vnorm;
  int a;
  int overlaps_interface;
  int got_interface = FALSE;
  SGRID *element_search_grid = NULL;
  int pass, num_passes;

  for (ebi = 0; ebi < exo->num_elem_blocks; ebi++) {

    pd = pd_glob[Matilda[ebi]];
    mp = mp_glob[Matilda[ebi]];

    e_start = exo->eb_ptr[ebi];
    e_end = exo->eb_ptr[ebi + 1];

    dim = pd->Num_Dim;
    wim = dim;

    if (pd->CoordinateSystem == SWIRLING ||
          pd->CoordinateSystem == PROJECTED_CARTESIAN ||
          pd->CoordinateSystem == CARTESIAN_2pt5D)
      wim = wim + 1;

    if (pd->v[pg->imtrx][ls->var]) {
      for (ielem = e_start; ielem < e_end; ielem++) {
        overlaps_interface = FALSE;
        if (ls->Length_Scale != 0.) /* diffuse interface */
        {
          if (elem_overlaps_interface(ielem, x, exo, ls->Length_Scale))
            overlaps_interface = TRUE;
        } else /* sharp interface */
        {
          if (elem_on_isosurface(ielem, x, exo, ls->var, 0.))
            overlaps_interface = TRUE;
        }
        if (overlaps_interface)
          got_interface = TRUE;
        if (!overlaps_interface)
          continue;

        load_elem_dofptr(ielem, exo, x, x_old, xdot, xdot_old, 0);

        h_elem_siz(hsquared, hhv, dhv_dxnode, pd->e[pg->imtrx][R_MESH1]);

        h_elem = 0.;
        for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++)
          h_elem += hsquared[a];
        /* This is the size of the element */
        h_elem = sqrt(h_elem / ((double)ei[pg->imtrx]->ielem_dim));

        if (ls->Integration_Depth > 0) {
          ip_total = get_subgrid_integration_pts(
              Subgrid_Tree, &element_search_grid, &Subgrid_Int.s,
              &Subgrid_Int.wt, ls->Length_Scale);
        } else if (ls->SubElemIntegration) {
          ip_total = get_subelement_integration_pts(
              &Subgrid_Int.s, &Subgrid_Int.wt, &Subgrid_Int.ip_sign, 0., -1, 0);
        } else {
          ip_total = elem_info(NQUAD, ei[pg->imtrx]->ielem_type);
        }

        if (is_xfem_interp(pd->i[pd->mi[VELOCITY1]][VELOCITY1]))
          num_passes = 2;
        else
          num_passes = 1;

        for (pass = 0; pass < num_passes; pass++) {
          if (num_passes == 2)
            ls->Elem_Sign = -1 + 2 * pass;
          else
            ls->Elem_Sign = 0;

          sum = 0.;
          sumv = 0.;
          for (ip = 0; ip < ip_total; ip++) {
            if (ls->Integration_Depth > 0) {
              xi = Subgrid_Int.s[ip];
              wt = Subgrid_Int.wt[ip];
            } else if (ls->SubElemIntegration) {
              xi = Subgrid_Int.s[ip];
              wt = Subgrid_Int.wt[ip];
              /* on sharp interface */
              ls->on_sharp_surf = TRUE;
            } else {
              xi = &(xi_array[0]);
              find_stu(ip, ei[pg->imtrx]->ielem_type, xi, xi + 1,
                       xi + 2); /* find quadrature point */
              wt = Gq_weight(
                  ip,
                  ei[pg->imtrx]->ielem_type); /* find quadrature weights for */
            }

            setup_shop_at_point(ei[pg->imtrx]->ielem, xi, exo);
            load_lsi(ls->Length_Scale);

            fv->wt = wt * lsi->delta;

            /* prediction of normal velocity */
            vnorm = 0.;
            if (pd->v[pd->mi[EXT_VELOCITY]][EXT_VELOCITY]) {
              vnorm += (2.5 * fv->ext_v - 1.5 * fv_old->ext_v);
            }
            if (pd->v[pd->mi[VELOCITY1]][VELOCITY1] &&
                tran->Fill_Equation != FILL_EQN_EXT_V) {
              for (a = 0; a < VIM; a++) {
                vnorm += (2.5 * fv->v[a] - 1.5 * fv_old->v[a]) * lsi->normal[a];
              }
            }
            sumv += fv->wt * vnorm;
            sum += fv->wt;
          }

          if ((sumv != 0.) && (sum != 0.)) {
            dt = fabs(h_elem * sum / sumv);
            if (count++ == 0 || dt < min_dt)
              min_dt = dt;
          }

        } /* for ( pass = 0; pass < num_passes; pass++ ) */

        if (ls->Integration_Depth > 0) {
          safe_free((void *)Subgrid_Int.s);
          Subgrid_Int.s = NULL;
          safe_free((void *)Subgrid_Int.wt);
          Subgrid_Int.wt = NULL;
          free_search_grid(&element_search_grid);
        }
      }
    }
  }
  /*fprintf(stderr,"Min dt = %g\n",min_dt);*/

  /* Search for global min_dt */
  if (Num_Proc > 1) {
    /* If interface not on this processor, don't allow zero min_dt! */
    if (!got_interface)
      min_dt = 100.0 * tran->Delta_t_max;
    min_dt = AZ_gmin_double(min_dt, proc_config);
  }

  /* restore */
  ls->on_sharp_surf = FALSE;
  ls->Elem_Sign = 0;
  return (min_dt);
}
#endif

int get_subelement_integration_pts(double (**s)[DIM], double **weight,
                                   int **ip_sign, double isoval, int gpt_type,
                                   int sign) {

  Integ_Elem *e;
  int num_gpts;

  if (pd->v[pg->imtrx][LS]) {
    neg_elem_volume |= fail_courant_condition();
  } else {
    neg_elem_volume = 0;
  }

  if (neg_elem_volume)
    return (0);

  e = create_integ_elements(isoval);

  num_gpts = num_subelement_integration_pts(e, gpt_type, sign);

  if (*s != NULL) {
    safe_free((void *)*s);
    *s = NULL;
  }
  if (*weight != NULL) {
    safe_free((void *)*weight);
    *weight = NULL;
  }
  if (ip_sign != NULL && *ip_sign != NULL) {
    safe_free((void *)*ip_sign);
    *ip_sign = NULL;
  }

  *s = (double(*)[DIM])smalloc(DIM * num_gpts * sizeof(double));

  *weight = (double *)smalloc(num_gpts * sizeof(double));

  if (ip_sign != NULL) {
    *ip_sign = (int *)smalloc(num_gpts * sizeof(int));
    gather_subelement_integration_pts(e, *s, *weight, *ip_sign, gpt_type, sign,
                                      0);
  } else {
    gather_subelement_integration_pts(e, *s, *weight, NULL, gpt_type, sign, 0);
  }

  free_integ_elements(e);

  return (num_gpts);
}

void get_subelement_facets(struct LS_Surf_List *list, double isoval) {

  Integ_Elem *e;

  e = create_integ_elements(isoval);

  gather_subelement_facets(list, e);

  free_integ_elements(e);
}

void gather_subelement_facets(struct LS_Surf_List *list, Integ_Elem *e) {
  if (e->num_subelements > 0) {
    int i;
    for (i = 0; i < e->num_subelements; i++) {
      gather_subelement_facets(list, e->subelements[i]);
    }
  } else {
    int nodes_per_side;
    int lnn[2];
    struct LS_Surf *vertex[2], *surf;
    int iside;
    double x[DIM] = {0., 0., 0.};

    for (iside = 0; iside < e->num_sides; iside++) {
      if (e->bc_sides[iside]) {
        get_side_info(e->ielem_type, iside + 1, &nodes_per_side, lnn);

        if (nodes_per_side == 2) {
          map_local_coordinates(e->xi[lnn[0]], x);
          vertex[0] =
              create_surf_point(x, ei[pg->imtrx]->ielem, e->xi[lnn[0]], FALSE);

          map_local_coordinates(e->xi[lnn[1]], x);
          vertex[1] =
              create_surf_point(x, ei[pg->imtrx]->ielem, e->xi[lnn[1]], FALSE);

          surf = create_surf_facet_line(vertex[0], vertex[1],
                                        ei[pg->imtrx]->ielem, -1);
          append_surf(list, surf);
        } else {
          /* DRN - higher order facets would be a much better solution */
          int num_facets = 20;
          double yi[DIM], yi0[DIM], yi1[DIM], xi[DIM];
          int i;
          double s0, s1;

          find_nodal_stu(lnn[0], e->ielem_type, yi0, yi0 + 1, yi0 + 2);
          find_nodal_stu(lnn[1], e->ielem_type, yi1, yi1 + 1, yi1 + 2);

          for (i = 0; i < num_facets; i++) {
            s0 = ((double)(2 * i - num_facets)) / num_facets;
            s1 = ((double)(2 * (i + 1) - num_facets)) / num_facets;

            yi[0] = yi0[0] + 0.5 * (s0 + 1.) * (yi1[0] - yi0[0]);
            yi[1] = yi0[1] + 0.5 * (s0 + 1.) * (yi1[1] - yi0[1]);
            yi[2] = 0.;

            map_subelement_stu(xi, e, yi);
            map_local_coordinates(xi, x);
            vertex[0] = create_surf_point(x, ei[pg->imtrx]->ielem, xi, FALSE);
            if (fabs(x[2]) > 1.e-6) {
              printf("bad1\n");
            }

            yi[0] = yi0[0] + 0.5 * (s1 + 1.) * (yi1[0] - yi0[0]);
            yi[1] = yi0[1] + 0.5 * (s1 + 1.) * (yi1[1] - yi0[1]);
            yi[2] = 0.;

            map_subelement_stu(xi, e, yi);
            map_local_coordinates(xi, x);
            vertex[1] = create_surf_point(x, ei[pg->imtrx]->ielem, xi, FALSE);
            if (fabs(x[2]) > 1.e-6) {
              printf("bad2\n");
            }

            surf = create_surf_facet_line(vertex[0], vertex[1],
                                          ei[pg->imtrx]->ielem, -1);
            append_surf(list, surf);
          }
        }
      }
    }
  }
}

Integ_Elem *create_integ_elements(double isoval) {
  Integ_Elem *e;
  double(*xi)[DIM];
  int i;
  int side_ids[4] = {0, 1, 2, 3};

  e = (Integ_Elem *)smalloc(sizeof(Integ_Elem));

  xi = (double(*)[DIM])smalloc(sizeof(double) * ei[pg->imtrx]->num_local_nodes *
                               DIM);

  for (i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
    find_nodal_stu(i, ei[pg->imtrx]->ielem_type, &xi[i][0], &xi[i][1],
                   &xi[i][2]);
  }

  build_integ_element(e, isoval, ei[pg->imtrx]->ielem_type, xi, side_ids,
                      FALSE);

  safe_free(xi);

  return (e);
}

void free_integ_elements(Integ_Elem *e) {
  int i;
  for (i = 0; i < e->num_subelements; i++) {
    free_integ_elements(e->subelements[i]);
  }
  if (e->num_subelements > 0)
    safe_free(e->subelements);
  if (e->bc_sides != NULL)
    safe_free(e->bc_sides);
  if (e->side_ids != NULL)
    safe_free(e->side_ids);

  safe_free(e->f);
  safe_free(e->xi);
  safe_free(e);
}
/*
int
significant_subelement_crossing( Integ_Elem * e, double tol )
{
  double max_neg = 0.;
  double max_pos = 0.;
  int i;

  for ( i = 0; i < e->num_local_nodes; i++ ) {
    if ( -e->f[i] > max_neg ) max_neg = -e->f[i];
    if ( e->f[i] > max_pos ) max_pos = e->f[i];
  }

  if ( max_neg > max_pos )  {
    if ( max_pos/max_neg > tol ) return(TRUE);
  } else {
    if ( max_neg/max_pos > tol ) return(TRUE);
  }
  return(FALSE);
}
*/
/*
static double
subelement_error(Integ_Elem * e, double isoval)
{
  int i;
  double phi;
  double xi[3], yi[3];
  double f = 0., fbase;
  double fmax = e->f[0];
  double fmin = e->f[0];
  double error;
  double fscale = 0.;
  xi[0] = 0.; xi[1] = 0.; xi[2] = 0.;
  yi[0] = 1./3.; yi[1] = 1./3.; yi[2] = 0.;

  for ( i = 0; i < e->num_local_nodes; i++ )
    {
      phi = shape(yi[0], yi[1], yi[2], e->ielem_type, PSI, i);
      xi[0] += phi * e->xi[i][0];
      xi[1] += phi * e->xi[i][1];
      xi[2] += phi * e->xi[i][2];
      if ( e->f[i] > fmax ) fmax = e->f[i];
      if ( e->f[i] < fmin ) fmin = e->f[i];
      fscale += e->f[i] * e->f[i];
      f += phi * e->f[i];
    }
  fbase = fv_at_point( xi, ls->var ) - isoval;
  error = fabs((fbase - f)/(fmax - fmin));
  error = fabs(fbase - f)/sqrt(fscale/e->num_local_nodes);
  return error;
}
*/
static double subelement_quality(Integ_Elem *e) {
  double yi[3];
  double detJ, mindetJ = 1.e33, sumdetJ = 0.;
  double quality;
  int ip_total, ip;

  ip_total = elem_info(NQUAD, e->ielem_type);
  for (ip = 0; ip < ip_total; ip++) {
    find_stu(ip, e->ielem_type, yi, yi + 1, yi + 2);
    detJ = subelement_detJ(e, yi, TRUE);
    if (ip == 0 || detJ < mindetJ)
      mindetJ = detJ;
    sumdetJ += fabs(detJ);
  }

  if (sumdetJ < 1.e-16)
    quality = 1.;
  else
    quality = mindetJ * e->num_local_nodes / sumdetJ;
  return quality;
}

static int subelement_side_crossing(Integ_Elem *e, int iside, double tol)
/* checks status of side of subelement
   return = 0: no crossing at all of this side of subelement
   return = 1: single significant crossing of this side of subelement
   return = 2: two crossings on this side of the subelement and extent of
   crossing is significant return = -1: at least one crossing of subelement but
   not considered significant
 */
{
  int nodes_per_side;
  int lnn[MAX_NODES_PER_SIDE];
  double max_neg = 0.;
  double max_pos = 0.;
  int num_crossings = 0;
  int i;

  get_side_info(e->ielem_type, iside + 1, &nodes_per_side, lnn);

  if (nodes_per_side == 2) {
    if (sign_change(e->f[lnn[0]], e->f[lnn[1]]))
      num_crossings++;
  } else if (nodes_per_side == 3) {
    if (sign_change(e->f[lnn[0]], e->f[lnn[2]]))
      num_crossings++;
    if (sign_change(e->f[lnn[2]], e->f[lnn[1]]))
      num_crossings++;
    if (num_crossings == 0 &&
        fabs(e->f[lnn[2]]) < 0.25 * fabs(e->f[lnn[0]] + e->f[lnn[1]]) -
                                 0.5 * sqrt(e->f[lnn[0]] * e->f[lnn[1]]))
      num_crossings = 2;
  } else {
    EH(-1, "Subelement type not supported.");
  }

  /* because of all the bad things that can happen when a side has two
     crossings, we return early in this case
   */
  if (num_crossings == 2)
    return (2);

  /* There is one more check we need to make.
   * Namely, even if this side has no crossings, does it
   * exactly run right along the interface?
   * If so, the extrema along this side will all have
   * abs(F) == 0.
   */

  /* find nodal extrema */
  for (i = 0; i < nodes_per_side; i++) {
    if (-e->f[lnn[i]] > max_neg)
      max_neg = -e->f[lnn[i]];
    if (e->f[lnn[i]] > max_pos)
      max_pos = e->f[lnn[i]];
  }

  /* calculate extreme value on sides considering higher order interpolant */
  if (nodes_per_side == 3) {
    double a, b, f;
    a = e->f[lnn[0]] - 2. * e->f[lnn[2]] + e->f[lnn[1]];
    b = e->f[lnn[1]] - e->f[lnn[0]];
    if (fabs(a) > 0.5 * fabs(b)) {
      f = e->f[lnn[2]] - b * b / (8. * a);
      if (-f > max_neg)
        max_neg = -f;
      if (f > max_pos)
        max_pos = f;
    }
  }

#if 0
  fprintf(stderr,"  e->sign = %d, side = %d, max_neg = %g, max_pos = %g\n",e->sign,iside,max_neg,max_pos);
#endif

  /* test for case of running along surface */
  if (max_neg == 0. && max_pos == 0.)
    return (-2);

  /* check for case of no or insignificant crossing */
  if ((e->sign == -1 && max_pos == 0.) || (e->sign == 1 && max_neg == 0.))
    return (-1);

  return (1);
}

void build_integ_element(Integ_Elem *e, double isoval, int ielem_type,
                         double (*xi)[DIM], int *side_ids, int is_conformal) {
  int i, shape;

  int tri_type;
  double sum_f;

  if (xfem != NULL && xfem->have_XG && ei[pg->imtrx]->ielem_type == BIQUAD_QUAD)
    tri_type = QUAD6_TRI;
  else
    tri_type = QUAD_TRI;

  e->ielem_type = ielem_type;
  e->num_local_nodes = elem_info(NNODES, ielem_type);

  e->xi = (double(*)[DIM])smalloc(sizeof(double) * e->num_local_nodes * DIM);
  e->f = (double *)smalloc(sizeof(double) * e->num_local_nodes);
  for (i = 0; i < e->num_local_nodes; i++) {
    e->xi[i][0] = xi[i][0];
    e->xi[i][1] = xi[i][1];
    e->xi[i][2] = xi[i][2];
    e->f[i] = fv_at_point(xi[i], ls->var) - isoval;
  }

  /* compute "sign" of subelement
     while this doesn't make sense for subelems that span interface
     it is needed for discerning if subelem is responsible for
     constructing interfacial elements
   */
  sum_f = 0.;
  for (i = 0; i < e->num_local_nodes; i++)
    sum_f += e->f[i];
  if (sum_f < 0.)
    e->sign = -1;
  else
    e->sign = 1;

  shape = type2shape(ielem_type);
  e->num_sides = shape2sides(shape);

  e->side_ids = (int *)smalloc(sizeof(int) * e->num_sides);
  for (i = 0; i < e->num_sides; i++) {
    e->side_ids[i] = side_ids[i];
  }

  /* For this next part, there are basically, two classes of elements
   *   1. elements that have been specifically created to conform to interface
   *   2. elements that do not necessarily conform to interface
   * These are distinguished based on the flag is_conformal.  For non-conformal
   * elements, we need to look for sub-elements.
   */

  e->num_subelements = 0;
  e->subelements = NULL;

  e->bc_sides = NULL;
  e->is_conformal = is_conformal;

  if (!e->is_conformal) {
    switch (ielem_type) {
    case BIQUAD_QUAD:
    case BILINEAR_QUAD: {
      int job = -1; /* job == -1: no subelements
                       job == 0: create 4 triangles
                     */
      double *f = e->f;
      int side_crossing[4];
      double nodes[6][DIM];
      int side_ids[3];

	    memset(side_crossing, 0, sizeof(int)*4);

      /* determine what we are going to do with this element (set job) */

      if (ls->Integration_Depth > 0) { /* diffuse interface */

        int overlaps_interface = FALSE;
        double alpha = 0.5 * ls->Length_Scale;
        int all_low = TRUE;
        int all_high = TRUE;

        /* see if subelement intersects diffuse interface */
        for (i = 0; i < e->num_local_nodes; i++) {
          all_low = (f[i] < -alpha) && all_low;
          all_high = (f[i] >= alpha) && all_high;
        }
        overlaps_interface = !(all_low || all_high);

        if (overlaps_interface)
          job = 0;

      } else { /* sharp interface */
#if 0
              int have_crossing = FALSE;
	      for ( iside = 0; iside < e->num_sides; iside++ )
	        {
		  side_crossing[iside] = subelement_side_crossing( e, iside, SUBELEM_SIG_CROSS_TOL );
		  if ( side_crossing[iside] > 0 ) have_crossing = TRUE;
		}

              if ( have_crossing ) job = 0;
#endif
#if 1
        /* always cut into triangles */
        job = 0;
#endif
      }

      if (job == 0) { /* create 4 triangles */

        e->num_subelements = 4;
        e->subelements =
            (Integ_Elem **)smalloc(sizeof(Integ_Elem) * e->num_subelements);

        /* triangle 1 */
        nodes[0][0] = e->xi[0][0];
        nodes[0][1] = e->xi[0][1];
        nodes[0][2] = 0.;
        nodes[1][0] = e->xi[1][0];
        nodes[1][1] = e->xi[1][1];
        nodes[1][2] = 0.;
        if (ielem_type == BIQUAD_QUAD) {
          nodes[2][0] = e->xi[8][0];
          nodes[2][1] = e->xi[8][1];
          nodes[2][2] = 0.;
          nodes[3][0] = e->xi[4][0];
          nodes[3][1] = e->xi[4][1];
          nodes[3][2] = 0.;
        } else {
          nodes[2][0] =
              0.25 * (e->xi[0][0] + e->xi[1][0] + e->xi[2][0] + e->xi[3][0]);
          nodes[2][1] =
              0.25 * (e->xi[0][1] + e->xi[1][1] + e->xi[2][1] + e->xi[3][1]);
          nodes[2][2] = 0.;
          nodes[3][0] = 0.5 * (nodes[0][0] + nodes[1][0]);
          nodes[3][1] = 0.5 * (nodes[0][1] + nodes[1][1]);
          nodes[3][2] = 0.;
        }
        nodes[4][0] = 0.5 * (nodes[2][0] + nodes[1][0]);
        nodes[4][1] = 0.5 * (nodes[2][1] + nodes[1][1]);
        nodes[4][2] = 0.;
        nodes[5][0] = 0.5 * (nodes[2][0] + nodes[0][0]);
        nodes[5][1] = 0.5 * (nodes[2][1] + nodes[0][1]);
        nodes[5][2] = 0.;
        side_ids[0] = e->side_ids[0];
        side_ids[1] = -1; /* not on any parent side */
        side_ids[2] = -1; /* not on any parent side */
        e->subelements[0] = (Integ_Elem *)smalloc(sizeof(Integ_Elem));
        build_integ_element(e->subelements[0], isoval, tri_type, nodes,
                            side_ids, FALSE);

        /* triangle 2 */
        nodes[0][0] = e->xi[1][0];
        nodes[0][1] = e->xi[1][1];
        nodes[0][2] = 0.;
        nodes[1][0] = e->xi[2][0];
        nodes[1][1] = e->xi[2][1];
        nodes[1][2] = 0.;
        if (ielem_type == BIQUAD_QUAD) {
          nodes[2][0] = e->xi[8][0];
          nodes[2][1] = e->xi[8][1];
          nodes[2][2] = 0.;
          nodes[3][0] = e->xi[5][0];
          nodes[3][1] = e->xi[5][1];
          nodes[3][2] = 0.;
        } else {
          nodes[2][0] =
              0.25 * (e->xi[0][0] + e->xi[1][0] + e->xi[2][0] + e->xi[3][0]);
          nodes[2][1] =
              0.25 * (e->xi[0][1] + e->xi[1][1] + e->xi[2][1] + e->xi[3][1]);
          nodes[2][2] = 0.;
          nodes[3][0] = 0.5 * (nodes[0][0] + nodes[1][0]);
          nodes[3][1] = 0.5 * (nodes[0][1] + nodes[1][1]);
          nodes[3][2] = 0.;
        }
        nodes[4][0] = 0.5 * (nodes[2][0] + nodes[1][0]);
        nodes[4][1] = 0.5 * (nodes[2][1] + nodes[1][1]);
        nodes[4][2] = 0.;
        nodes[5][0] = 0.5 * (nodes[2][0] + nodes[0][0]);
        nodes[5][1] = 0.5 * (nodes[2][1] + nodes[0][1]);
        nodes[5][2] = 0.;
        side_ids[0] = e->side_ids[1];
        side_ids[1] = -1; /* not on any parent side */
        side_ids[2] = -1; /* not on any parent side */
        e->subelements[1] = (Integ_Elem *)smalloc(sizeof(Integ_Elem));
        build_integ_element(e->subelements[1], isoval, tri_type, nodes,
                            side_ids, FALSE);

        /* triangle 3 */
        nodes[0][0] = e->xi[2][0];
        nodes[0][1] = e->xi[2][1];
        nodes[0][2] = 0.;
        nodes[1][0] = e->xi[3][0];
        nodes[1][1] = e->xi[3][1];
        nodes[1][2] = 0.;
        if (ielem_type == BIQUAD_QUAD) {
          nodes[2][0] = e->xi[8][0];
          nodes[2][1] = e->xi[8][1];
          nodes[2][2] = 0.;
          nodes[3][0] = e->xi[6][0];
          nodes[3][1] = e->xi[6][1];
          nodes[3][2] = 0.;
        } else {
          nodes[2][0] =
              0.25 * (e->xi[0][0] + e->xi[1][0] + e->xi[2][0] + e->xi[3][0]);
          nodes[2][1] =
              0.25 * (e->xi[0][1] + e->xi[1][1] + e->xi[2][1] + e->xi[3][1]);
          nodes[2][2] = 0.;
          nodes[3][0] = 0.5 * (nodes[0][0] + nodes[1][0]);
          nodes[3][1] = 0.5 * (nodes[0][1] + nodes[1][1]);
          nodes[3][2] = 0.;
        }
        nodes[4][0] = 0.5 * (nodes[2][0] + nodes[1][0]);
        nodes[4][1] = 0.5 * (nodes[2][1] + nodes[1][1]);
        nodes[4][2] = 0.;
        nodes[5][0] = 0.5 * (nodes[2][0] + nodes[0][0]);
        nodes[5][1] = 0.5 * (nodes[2][1] + nodes[0][1]);
        nodes[5][2] = 0.;
        side_ids[0] = e->side_ids[2];
        side_ids[1] = -1; /* not on any parent side */
        side_ids[2] = -1; /* not on any parent side */
        e->subelements[2] = (Integ_Elem *)smalloc(sizeof(Integ_Elem));
        build_integ_element(e->subelements[2], isoval, tri_type, nodes,
                            side_ids, FALSE);

        /* triangle 4 */
        nodes[0][0] = e->xi[3][0];
        nodes[0][1] = e->xi[3][1];
        nodes[0][2] = 0.;
        nodes[1][0] = e->xi[0][0];
        nodes[1][1] = e->xi[0][1];
        nodes[1][2] = 0.;
        if (ielem_type == BIQUAD_QUAD) {
          nodes[2][0] = e->xi[8][0];
          nodes[2][1] = e->xi[8][1];
          nodes[2][2] = 0.;
          nodes[3][0] = e->xi[7][0];
          nodes[3][1] = e->xi[7][1];
          nodes[3][2] = 0.;
        } else {
          nodes[2][0] =
              0.25 * (e->xi[0][0] + e->xi[1][0] + e->xi[2][0] + e->xi[3][0]);
          nodes[2][1] =
              0.25 * (e->xi[0][1] + e->xi[1][1] + e->xi[2][1] + e->xi[3][1]);
          nodes[2][2] = 0.;
          nodes[3][0] = 0.5 * (nodes[0][0] + nodes[1][0]);
          nodes[3][1] = 0.5 * (nodes[0][1] + nodes[1][1]);
          nodes[3][2] = 0.;
        }
        nodes[4][0] = 0.5 * (nodes[2][0] + nodes[1][0]);
        nodes[4][1] = 0.5 * (nodes[2][1] + nodes[1][1]);
        nodes[4][2] = 0.;
        nodes[5][0] = 0.5 * (nodes[2][0] + nodes[0][0]);
        nodes[5][1] = 0.5 * (nodes[2][1] + nodes[0][1]);
        nodes[5][2] = 0.;
        side_ids[0] = e->side_ids[3];
        side_ids[1] = -1; /* not on any parent side */
        side_ids[2] = -1; /* not on any parent side */
        e->subelements[3] = (Integ_Elem *)smalloc(sizeof(Integ_Elem));
        build_integ_element(e->subelements[3], isoval, tri_type, nodes,
                            side_ids, FALSE);

      } else { /* no subelements */

        if (ls->Integration_Depth > 0) { /* diffuse interface */
          /* assign sign of element */
          e->sign = 0;

          /* make sure bc_sides isn't used since this interface is diffuse */
          e->bc_sides = NULL;
        } else { /* sharp interface */
          /* see if any sides of this element are on interface */
          e->bc_sides = (int *)smalloc(sizeof(int) * e->num_sides);

          e->bc_sides[0] = (side_crossing[0] == -2);
          e->bc_sides[1] = (side_crossing[1] == -2);
          e->bc_sides[2] = (side_crossing[2] == -2);
          e->bc_sides[3] = (side_crossing[3] == -2);
        }
      }
    } break;

    case QUAD_TRI:
    case QUAD6_TRI: {
      int job = -1; /* job == -1: no subelements
                       job == 0:  create 4 conformal subelements
                       job == 1:  create 4 non-conformal subelements by
                       subdividing each side
                     */
      double *f = e->f;
      int side_crossing[3], iside;
      double nodes[6][DIM];
      int index[6];
      double xi0[DIM], xi1[DIM], xi2[DIM];
      int side_ids[3];
      double f_temp;
      int ii;

      /* initialize index[] */
      for (ii = 0; ii < 6; ii++) {
        index[ii] = 0;
      }

      /* determine what we are going to do with this element (set job) */

      if (ls->Integration_Depth > 0) {
        int overlaps_interface = FALSE;
        double alpha = 0.5 * ls->Length_Scale;
        double detJ, length_scale;
        double scale = pow(2.0, (double)ls->Integration_Depth);
        double yi[DIM],
            center[DIM] = {0.33333333333333333333, 0.33333333333333333333, 0.};
        int all_low = TRUE;
        int all_high = TRUE;

        /* see if subelement intersects diffuse interface */
        for (i = 0; i < e->num_local_nodes; i++) {
          all_low = (f[i] < -alpha) && all_low;
          all_high = (f[i] >= alpha) && all_high;
        }
        overlaps_interface = !(all_low || all_high);

        if (overlaps_interface) {
          /* determine size of subelement compared to interface thickness */
          map_subelement_stu(center, e, yi);
          detJ = subelement_detJ(e, yi, FALSE);
          length_scale = pow(detJ, 1.0 / ((double)pd->Num_Dim));
          if (length_scale > ls->Length_Scale / scale)
            job = 1;
        }

      } else { /* sharp interface */
        int have_crossing = FALSE;
        int two_crossings = FALSE;
#if 0
              fprintf(stderr,"Parent ielem+1 = %d\n",ei[pg->imtrx]->ielem+1);
              for ( i=0; i<e->num_local_nodes; i++ ) {
                fprintf(stderr,"  Node #%d: xi=(%g,%g), F=%g\n",
                        i, xi[i][0], xi[i][1], f[i]);
              }
#endif
        for (iside = 0; iside < e->num_sides; iside++) {
          side_crossing[iside] =
              subelement_side_crossing(e, iside, SUBELEM_SIG_CROSS_TOL);
          if (side_crossing[iside] > 0)
            have_crossing = TRUE;
          if (side_crossing[iside] == 2)
            two_crossings = TRUE;
        }
#if 0
              for ( i=0; i<e->num_sides; i++ ) {
                fprintf(stderr,"  side_crossing[%d] = %d\n", i, side_crossing[i]);
              }
#endif
        if (have_crossing) {
          if (two_crossings) {
            job = 1;
          } else {
            /* otherwise, attempt to create 4 conforming triangle sub-elements
             */
            double temp[DIM];
            int success;

            /* reorder so that first side does not have crossing */
            if (!sign_change(f[0], f[1])) {
              index[0] = 0;
              index[1] = 1;
              index[2] = 2;
              index[3] = 3;
              index[4] = 4;
              index[5] = 5;
            } else if (!sign_change(f[1], f[2])) {
              index[0] = 1;
              index[1] = 2;
              index[2] = 0;
              index[3] = 4;
              index[4] = 5;
              index[5] = 3;
            } else {
              index[0] = 2;
              index[1] = 0;
              index[2] = 1;
              index[3] = 5;
              index[4] = 3;
              index[5] = 4;
            }

            xi0[0] = xi[index[2]][0];
            xi0[1] = xi[index[2]][1];
            xi0[2] = xi[index[2]][2];
            temp[0] = xi[index[1]][0];
            temp[1] = xi[index[1]][1];
            temp[2] = xi[index[1]][2];
            success = find_link_intersection(xi0, temp, ls->var, isoval, NULL);

            xi2[0] = xi[index[2]][0];
            xi2[1] = xi[index[2]][1];
            xi2[2] = xi[index[2]][2];
            temp[0] = xi[index[0]][0];
            temp[1] = xi[index[0]][1];
            temp[2] = xi[index[0]][2];
            success &= find_link_intersection(xi2, temp, ls->var, isoval, NULL);
#if 0
                  /* DRN: I really liked this commented out algorithm, but it isn't quite as 
		     accurate as the current algorithm further below.  The issue is that, although
		     the quadratic triangle has parabolic edges, the underlying biquadratic function
		     has still higher order terms.  They are neglected in this algorithm and it
		     can reduce the accuracy for very fine meshes.
		   */
		  /* the center point isn't quite so easy */
                  /* we would like the center point to be located about halfway between the
                     other two crossings
                     So we will find a root along guaranteed direction and then equalize side
                   */
                  xi1[0] = xi[index[2]][0]; xi1[1] = xi[index[2]][1]; xi1[2] = xi[index[2]][2];
                  temp[0] = xi[index[3]][0]; temp[1] = xi[index[3]][1]; temp[2] = xi[index[3]][2];
                  success &= find_link_intersection( xi1, temp, ls->var, isoval, NULL );

                  if ( success )
                    {
                      /* now that we have xi0 xi1 and xi2, move xi1 to lie approx halfway between
                         xi0 and xi2
                       */
                      double dist0, dist2, mid;
                      temp[0] = xi1[0]; temp[1] = xi1[1]; temp[2] = 0.;
                      
                      dist0 = sqrt((xi1[0] - xi0[0]) * (xi1[0] - xi0[0]) +
                                   (xi1[1] - xi0[1]) * (xi1[1] - xi0[1]));
                      dist2 = sqrt((xi1[0] - xi2[0]) * (xi1[0] - xi2[0]) +
                                   (xi1[1] - xi2[1]) * (xi1[1] - xi2[1]));
                                   
                      if ( dist0 != dist2 )
                        {
                          double err = -1. + 2.*dist2/(dist0+dist2);
                          mid = 0.5*(-1. + sqrt(1.+4.*err*err))/err;
                          
                          xi1[0] = 0.5*mid*(mid-1.)  * xi0[0] +
                                   (1.-mid)*(1.+mid) * temp[0] +
                                   0.5*mid*(mid+1.)  * xi2[0];
                          xi1[1] = 0.5*mid*(mid-1.)  * xi0[1] +
                                   (1.-mid)*(1.+mid) * temp[1] +
                                   0.5*mid*(mid+1.)  * xi2[1];
#if 0
                          DPRINTF(stderr,"before dist0=%g,dist2=%g\n",dist0,dist2);         
                          dist0 = sqrt((xi1[0] - xi0[0]) * (xi1[0] - xi0[0]) +
                                       (xi1[1] - xi0[1]) * (xi1[1] - xi0[1]));
                          dist2 = sqrt((xi1[0] - xi2[0]) * (xi1[0] - xi2[0]) +
                                       (xi1[1] - xi2[1]) * (xi1[1] - xi2[1]));
                          DPRINTF(stderr,"after mid=%g,dist0=%g,dist2=%g\n",mid,dist0,dist2);
#endif
                        }
                    }
#endif
#if 1
            /* the center point isn't quite so easy */
            /* we would like the center point to be located about halfway
               between the other two crossings So we will find a root between
               the midpt and node 2 or node 3 depending on which one is on the
               opposite side of the interface from the midpoint.
             */
            temp[0] = 0.5 * (xi0[0] + xi2[0]);
            temp[1] = 0.5 * (xi0[1] + xi2[1]);
            temp[2] = 0.;
            f_temp = fv_at_point(temp, ls->var) - isoval;
            if (sign_change(f_temp, f[index[2]])) {
              xi1[0] = xi[index[2]][0];
              xi1[1] = xi[index[2]][1];
              xi1[2] = xi[index[2]][2];
              success &=
                  find_link_intersection(xi1, temp, ls->var, isoval, NULL);
            } else if (sign_change(f_temp, f[index[3]])) {
              xi1[0] = xi[index[3]][0];
              xi1[1] = xi[index[3]][1];
              xi1[2] = xi[index[3]][2];
              success &=
                  find_link_intersection(xi1, temp, ls->var, isoval, NULL);
            } else {
              success = FALSE;
            }
#endif

            if (success) {
              job = 0;
            } else {
              job = 1; /* divide into 4 nonconformal subelements */
                       /* (test print lines commented out)
                                           fprintf(stderr,"Failed to decompose QUAD_TRI
                          into subelements.\n");          for ( i=0;
                          i<e->num_local_nodes; i++ ) {          fprintf(stderr,"Element #%d: Node
                          #%d:          xi=(%g,%g), F=%g\n",          ei[pg->imtrx]->ielem, i,
                          xi[index[i]][0],          xi[index[i]][1], f[index[i]]);
                                           }
                                           for ( i=0; i<e->num_sides; i++ ) {
                                             fprintf(stderr,"side_crossing[%d] = %d\n",
                          i, side_crossing[i]);
                                           }
                       */
            }
          }
        }
      }
      if (job == 0) { /* divide into 4 conformal elements */

        int err = FALSE;
        double quality_tol = 0.3;

        e->num_subelements = 4;
        e->subelements =
            (Integ_Elem **)smalloc(sizeof(Integ_Elem) * e->num_subelements);

        /* triangle 1 */
        nodes[0][0] = xi0[0];
        nodes[0][1] = xi0[1];
        nodes[0][2] = 0.;
        nodes[1][0] = xi[index[2]][0];
        nodes[1][1] = xi[index[2]][1];
        nodes[1][2] = 0.;
        nodes[2][0] = xi2[0];
        nodes[2][1] = xi2[1];
        nodes[2][2] = 0.;
        nodes[3][0] = 0.5 * (xi0[0] + xi[index[2]][0]);
        nodes[3][1] = 0.5 * (xi0[1] + xi[index[2]][1]);
        nodes[3][2] = 0.;
        nodes[4][0] = 0.5 * (xi2[0] + xi[index[2]][0]);
        nodes[4][1] = 0.5 * (xi2[1] + xi[index[2]][1]);
        nodes[4][2] = 0.;
        nodes[5][0] = xi1[0];
        nodes[5][1] = xi1[1];
        nodes[5][2] = 0.;

        side_ids[0] = e->side_ids[index[1]];
        side_ids[1] = e->side_ids[index[2]];
        side_ids[2] = -2; /* on isosurface */

        e->subelements[0] = (Integ_Elem *)smalloc(sizeof(Integ_Elem));
        build_integ_element(e->subelements[0], isoval, tri_type, nodes,
                            side_ids, TRUE);
        if (subelement_quality(e->subelements[0]) < quality_tol)
          err = TRUE;

        /* triangle 2 */
        nodes[0][0] = xi0[0];
        nodes[0][1] = xi0[1];
        nodes[0][2] = 0.;
        nodes[1][0] = xi2[0];
        nodes[1][1] = xi2[1];
        nodes[1][2] = 0.;
        nodes[2][0] = xi[index[3]][0];
        nodes[2][1] = xi[index[3]][1];
        nodes[2][2] = 0.;
        nodes[3][0] = xi1[0];
        nodes[3][1] = xi1[1];
        nodes[3][2] = 0.;
        nodes[4][0] = 0.5 * (xi2[0] + xi[index[3]][0]);
        nodes[4][1] = 0.5 * (xi2[1] + xi[index[3]][1]);
        nodes[4][2] = 0.;
        nodes[5][0] = 0.5 * (xi0[0] + xi[index[3]][0]);
        nodes[5][1] = 0.5 * (xi0[1] + xi[index[3]][1]);
        nodes[5][2] = 0.;

        side_ids[0] = -2; /* on isosurface */
        side_ids[1] = -1; /* not on any parent side */
        side_ids[2] = -1; /* not on any parent side */

        e->subelements[1] = (Integ_Elem *)smalloc(sizeof(Integ_Elem));
        build_integ_element(e->subelements[1], isoval, tri_type, nodes,
                            side_ids, TRUE);
        if (subelement_quality(e->subelements[1]) < quality_tol)
          err = TRUE;

        /* triangle 3 */
        nodes[0][0] = xi[index[1]][0];
        nodes[0][1] = xi[index[1]][1];
        nodes[0][2] = 0.;
        nodes[1][0] = xi0[0];
        nodes[1][1] = xi0[1];
        nodes[1][2] = 0.;
        nodes[2][0] = xi[index[3]][0];
        nodes[2][1] = xi[index[3]][1];
        nodes[2][2] = 0.;
        nodes[3][0] = 0.5 * (xi0[0] + xi[index[1]][0]);
        nodes[3][1] = 0.5 * (xi0[1] + xi[index[1]][1]);
        nodes[3][2] = 0.;
        nodes[4][0] = 0.5 * (xi0[0] + xi[index[3]][0]);
        nodes[4][1] = 0.5 * (xi0[1] + xi[index[3]][1]);
        nodes[4][2] = 0.;
        nodes[5][0] = 0.5 * (xi[index[1]][0] + xi[index[3]][0]);
        nodes[5][1] = 0.5 * (xi[index[1]][1] + xi[index[3]][1]);
        nodes[5][2] = 0.;

        side_ids[0] = e->side_ids[index[1]];
        side_ids[1] = -1; /* not on any parent side */
        side_ids[2] = e->side_ids[index[0]];

        e->subelements[2] = (Integ_Elem *)smalloc(sizeof(Integ_Elem));
        build_integ_element(e->subelements[2], isoval, tri_type, nodes,
                            side_ids, TRUE);
        if (subelement_quality(e->subelements[2]) < quality_tol)
          err = TRUE;

        /* triangle 4 */
        nodes[0][0] = xi[index[3]][0];
        nodes[0][1] = xi[index[3]][1];
        nodes[0][2] = 0.;
        nodes[1][0] = xi2[0];
        nodes[1][1] = xi2[1];
        nodes[1][2] = 0.;
        nodes[2][0] = xi[index[0]][0];
        nodes[2][1] = xi[index[0]][1];
        nodes[2][2] = 0.;
        nodes[3][0] = 0.5 * (xi2[0] + xi[index[3]][0]);
        nodes[3][1] = 0.5 * (xi2[1] + xi[index[3]][1]);
        nodes[3][2] = 0.;
        nodes[4][0] = 0.5 * (xi2[0] + xi[index[0]][0]);
        nodes[4][1] = 0.5 * (xi2[1] + xi[index[0]][1]);
        nodes[4][2] = 0.;
        nodes[5][0] = 0.5 * (xi[index[0]][0] + xi[index[3]][0]);
        nodes[5][1] = 0.5 * (xi[index[0]][1] + xi[index[3]][1]);
        nodes[5][2] = 0.;

        side_ids[0] = -1; /* not on any parent side */
        side_ids[1] = e->side_ids[index[2]];
        side_ids[2] = e->side_ids[index[0]];

        e->subelements[3] = (Integ_Elem *)smalloc(sizeof(Integ_Elem));
        build_integ_element(e->subelements[3], isoval, tri_type, nodes,
                            side_ids, TRUE);
        if (subelement_quality(e->subelements[3]) < quality_tol)
          err = TRUE;

        if (err) {
          job = 1;
          err = FALSE;
          /* clear memory for unsuccessful subelements */
          for (i = 0; i < 3; i++) {
            free_integ_elements(e->subelements[i]);
          }
          safe_free(e->subelements);
        }
      }
      if (job == 0) {
      } else if (job == 1) { /* divide into 4 nonconformal elements */

        int i0, i1, i2;

        e->num_subelements = 4;
        e->subelements =
            (Integ_Elem **)smalloc(sizeof(Integ_Elem) * e->num_subelements);

        /* triangle 1*/
        i0 = 0;
        i1 = 3;
        i2 = 5;
        nodes[0][0] = xi[i0][0];
        nodes[0][1] = xi[i0][1];
        nodes[0][2] = 0.;
        nodes[1][0] = xi[i1][0];
        nodes[1][1] = xi[i1][1];
        nodes[1][2] = 0.;
        nodes[2][0] = xi[i2][0];
        nodes[2][1] = xi[i2][1];
        nodes[2][2] = 0.;
        nodes[3][0] = 0.5 * (nodes[0][0] + nodes[1][0]);
        nodes[3][1] = 0.5 * (nodes[0][1] + nodes[1][1]);
        nodes[3][2] = 0.;
        nodes[4][0] = 0.5 * (nodes[1][0] + nodes[2][0]);
        nodes[4][1] = 0.5 * (nodes[1][1] + nodes[2][1]);
        nodes[4][2] = 0.;
        nodes[5][0] = 0.5 * (nodes[2][0] + nodes[0][0]);
        nodes[5][1] = 0.5 * (nodes[2][1] + nodes[0][1]);
        nodes[5][2] = 0.;

        side_ids[0] = e->side_ids[0];
        side_ids[1] = -1; /* not on any parent side */
        side_ids[2] = e->side_ids[2];

        e->subelements[0] = (Integ_Elem *)smalloc(sizeof(Integ_Elem));
        build_integ_element(e->subelements[0], isoval, tri_type, nodes,
                            side_ids, FALSE);

        /* triangle 2*/
        i0 = 1;
        i1 = 4;
        i2 = 3;
        nodes[0][0] = xi[i0][0];
        nodes[0][1] = xi[i0][1];
        nodes[0][2] = 0.;
        nodes[1][0] = xi[i1][0];
        nodes[1][1] = xi[i1][1];
        nodes[1][2] = 0.;
        nodes[2][0] = xi[i2][0];
        nodes[2][1] = xi[i2][1];
        nodes[2][2] = 0.;
        nodes[3][0] = 0.5 * (nodes[0][0] + nodes[1][0]);
        nodes[3][1] = 0.5 * (nodes[0][1] + nodes[1][1]);
        nodes[3][2] = 0.;
        nodes[4][0] = 0.5 * (nodes[1][0] + nodes[2][0]);
        nodes[4][1] = 0.5 * (nodes[1][1] + nodes[2][1]);
        nodes[4][2] = 0.;
        nodes[5][0] = 0.5 * (nodes[2][0] + nodes[0][0]);
        nodes[5][1] = 0.5 * (nodes[2][1] + nodes[0][1]);
        nodes[5][2] = 0.;

        side_ids[0] = e->side_ids[1];
        side_ids[1] = -1; /* not on any parent side */
        side_ids[2] = e->side_ids[0];

        e->subelements[1] = (Integ_Elem *)smalloc(sizeof(Integ_Elem));
        build_integ_element(e->subelements[1], isoval, tri_type, nodes,
                            side_ids, FALSE);

        /* triangle 3*/
        i0 = 2;
        i1 = 5;
        i2 = 4;
        nodes[0][0] = xi[i0][0];
        nodes[0][1] = xi[i0][1];
        nodes[0][2] = 0.;
        nodes[1][0] = xi[i1][0];
        nodes[1][1] = xi[i1][1];
        nodes[1][2] = 0.;
        nodes[2][0] = xi[i2][0];
        nodes[2][1] = xi[i2][1];
        nodes[2][2] = 0.;
        nodes[3][0] = 0.5 * (nodes[0][0] + nodes[1][0]);
        nodes[3][1] = 0.5 * (nodes[0][1] + nodes[1][1]);
        nodes[3][2] = 0.;
        nodes[4][0] = 0.5 * (nodes[1][0] + nodes[2][0]);
        nodes[4][1] = 0.5 * (nodes[1][1] + nodes[2][1]);
        nodes[4][2] = 0.;
        nodes[5][0] = 0.5 * (nodes[2][0] + nodes[0][0]);
        nodes[5][1] = 0.5 * (nodes[2][1] + nodes[0][1]);
        nodes[5][2] = 0.;

        side_ids[0] = e->side_ids[2];
        side_ids[1] = -1; /* not on any parent side */
        side_ids[2] = e->side_ids[1];

        e->subelements[2] = (Integ_Elem *)smalloc(sizeof(Integ_Elem));
        build_integ_element(e->subelements[2], isoval, tri_type, nodes,
                            side_ids, FALSE);

        /* triangle 4*/
        i0 = 3;
        i1 = 4;
        i2 = 5;
        nodes[0][0] = xi[i0][0];
        nodes[0][1] = xi[i0][1];
        nodes[0][2] = 0.;
        nodes[1][0] = xi[i1][0];
        nodes[1][1] = xi[i1][1];
        nodes[1][2] = 0.;
        nodes[2][0] = xi[i2][0];
        nodes[2][1] = xi[i2][1];
        nodes[2][2] = 0.;
        nodes[3][0] = 0.5 * (nodes[0][0] + nodes[1][0]);
        nodes[3][1] = 0.5 * (nodes[0][1] + nodes[1][1]);
        nodes[3][2] = 0.;
        nodes[4][0] = 0.5 * (nodes[1][0] + nodes[2][0]);
        nodes[4][1] = 0.5 * (nodes[1][1] + nodes[2][1]);
        nodes[4][2] = 0.;
        nodes[5][0] = 0.5 * (nodes[2][0] + nodes[0][0]);
        nodes[5][1] = 0.5 * (nodes[2][1] + nodes[0][1]);
        nodes[5][2] = 0.;

        side_ids[0] = -1; /* not on any parent side */
        side_ids[1] = -1; /* not on any parent side */
        side_ids[2] = -1; /* not on any parent side */

        e->subelements[3] = (Integ_Elem *)smalloc(sizeof(Integ_Elem));
        build_integ_element(e->subelements[3], isoval, tri_type, nodes,
                            side_ids, FALSE);

      } else { /* no subelements */
        if (ls->Integration_Depth > 0) {
          /* assign sign of element */
          e->sign = 0;

          /* make sure bc_sides isn't used since this interface is diffuse */
          e->bc_sides = NULL;
        } else {
          /* see if any sides of this element are on interface */
          e->bc_sides = (int *)smalloc(sizeof(int) * e->num_sides);
          e->bc_sides[0] = (side_crossing[0] == -2);
          e->bc_sides[1] = (side_crossing[1] == -2);
          e->bc_sides[2] = (side_crossing[2] == -2);
        }
      }
    } break;
    default:
      EH(-1, "Unsupported element type.");
      break;
    }
  } else /* if ( !is_conformal ) */
  {
    /* these are elems that were designed to conform to interface */
    /* we need to assign element sign and mark bc_sides */
    switch (ielem_type) {
    case QUAD_TRI:
    case QUAD6_TRI: {
      /* assign sign of element */
      e->sign = 1;
      if (side_ids[0] != -2 && side_ids[2] != -2) {
        if (e->f[0] < 0.)
          e->sign = -1;
      } else if (side_ids[0] != -2 && side_ids[1] != -2) {
        if (e->f[1] < 0.)
          e->sign = -1;
      } else if (side_ids[1] != -2 && side_ids[2] != -2) {
        if (e->f[2] < 0.)
          e->sign = -1;
      } else {
        EH(-1, "This really shouldn't happen!");
      }

      /* see if any sides of this element are on interface */
      e->bc_sides = (int *)smalloc(sizeof(int) * e->num_sides);
      e->bc_sides[0] = (e->sign == -1 && side_ids[0] == -2);
      e->bc_sides[1] = (e->sign == -1 && side_ids[1] == -2);
      e->bc_sides[2] = (e->sign == -1 && side_ids[2] == -2);
    } break;
    default:
      EH(-1, "Unsupported element type.");
      break;
    }
  }
}

void map_subelement_stu(double *xi, Integ_Elem *e, double *yi) {
  int i;
  double phi;

  xi[0] = 0.;
  xi[1] = 0.;
  xi[2] = 0.;

  for (i = 0; i < e->num_local_nodes; i++) {
    phi = shape(yi[0], yi[1], yi[2], e->ielem_type, PSI, i);
    xi[0] += phi * e->xi[i][0];
    xi[1] += phi * e->xi[i][1];
    xi[2] += phi * e->xi[i][2];
  }
}

void subelement_J(Integ_Elem *e, double *yi, double J[DIM][DIM],
                  int subelement_map) {
  int a, b, i;
  double dphidyi[MDE][DIM];
  int dim = ei[pg->imtrx]->ielem_dim;

  switch (dim) {
  case 1:
    for (i = 0; i < e->num_local_nodes; i++) {
      dphidyi[i][0] = shape(yi[0], yi[1], yi[2], e->ielem_type, DPSI_S, i);
    }
    break;
  case 2:
    for (i = 0; i < e->num_local_nodes; i++) {
      dphidyi[i][0] = shape(yi[0], yi[1], yi[2], e->ielem_type, DPSI_S, i);
      dphidyi[i][1] = shape(yi[0], yi[1], yi[2], e->ielem_type, DPSI_T, i);
    }
    break;
  case 3:
    for (i = 0; i < e->num_local_nodes; i++) {
      dphidyi[i][0] = shape(yi[0], yi[1], yi[2], e->ielem_type, DPSI_S, i);
      dphidyi[i][1] = shape(yi[0], yi[1], yi[2], e->ielem_type, DPSI_T, i);
      dphidyi[i][2] = shape(yi[0], yi[1], yi[2], e->ielem_type, DPSI_U, i);
    }
    break;
  }

  if (subelement_map) /* Jacobian mapping back to base element coords */
  {
    for (a = 0; a < dim; a++) {
      for (b = 0; b < dim; b++) {
        J[a][b] = 0.;
        for (i = 0; i < e->num_local_nodes; i++) {
          J[a][b] += e->xi[i][b] * dphidyi[i][a];
        }
      }
    }
  } else /* Jacobian mapping back to physical coords */
  {
    double x[MDE][DIM];
    for (i = 0; i < e->num_local_nodes; i++) {
      map_local_coordinates(e->xi[i], x[i]);
    }
    for (a = 0; a < dim; a++) {
      for (b = 0; b < dim; b++) {
        J[a][b] = 0.;
        for (i = 0; i < e->num_local_nodes; i++) {
          J[a][b] += x[i][b] * dphidyi[i][a];
        }
      }
    }
  }
}

double subelement_detJ(Integ_Elem *e, double *yi, int subelement_map) {
  double J[DIM][DIM], detJ = 0.0;
  int dim = ei[pg->imtrx]->ielem_dim;
  int elem_shape;

  subelement_J(e, yi, J, subelement_map);

  switch (dim) {
  case 1:
    detJ = J[0][0];
    break;
  case 2:
    detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    break;
  case 3:
    detJ = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) -
           J[0][1] * (J[1][0] * J[2][2] - J[2][0] * J[1][2]) +
           J[0][2] * (J[1][0] * J[2][1] - J[2][0] * J[1][1]);

    break;
  }

  if (subelement_map) {
    elem_shape = type2shape(e->ielem_type);
    switch (elem_shape) {
    case TRIANGLE:
      detJ *= 1.;
      break;
    case QUADRILATERAL:
      detJ *= 2.;
      break;
    }

    elem_shape = type2shape(ei[pg->imtrx]->ielem_type);
    switch (elem_shape) {
    case TRIANGLE:
      detJ /= 1.;
      break;
    case QUADRILATERAL:
      detJ /= 2.;
      break;
    }
  }

  return detJ;
}

double subelement_surfdet(Integ_Elem *e, double *yi, int side,
                          int subelement_map) {
  double J[DIM][DIM];
  double xi[DIM];
  int dim = ei[pg->imtrx]->ielem_dim;
  int ielem_shape;
  int nodes_per_side;
  int local_elem_node_id[MAX_NODES_PER_SIDE];
  int a, b;
  int DeformingMesh = pd->e[pg->imtrx][R_MESH1];
  struct Basis_Functions *map_bf = bf[pd->ShapeVar];

  /* the requirements here can be a bit different than in for subelement_detJ */
  /* There we computed |dX/dY| where Y is the subelent coord and X is the
     element coord.  For the volume integrals this get multiplied by |dx/dX|
     where x is the physical space.
   */
  /* If we are integrating the parent element along one of its side, we
     operate is the same mode here. (subelement_map == TRUE)
   */
  /* But if we are integrating along the isosurface (subelement_map == FALSE),
     we need the complete |dx/ds| where x is again the physical space
     and s is along the side of our subelement (in y coordinates).  This
     will be used directly by routines assembling embedded surface integrals.
   */

  subelement_J(e, yi, J, subelement_map);

  for (a = 0; a < dim; a++) {
    for (b = 0; b < dim; b++) {
      map_bf->J[a][b] = J[a][b];
    }
  }

  /* Quite possibly a bad idea, but we are going to try to get
     surface_determinant_and_normal to do the rest of the work.
     How is this accomplished?  Well it involves lying to the routine
     about the element properties and then resetting them upon return.
   */
  /*
     First step is to get coordinate scales loaded for coord
     systems like CYLINDRICAL.
     We do this if we are on isosurface so computing complete
     determinant (subelement_map == FALSE).
   */
  if (!subelement_map) {
    map_subelement_stu(xi, e, yi);
    map_local_coordinates(xi, fv->x);
    load_coordinate_scales(pd->CoordinateSystem, fv);
  }

  /* ick, there is an assumption in surface_determinant_and_normal that we
     violate here.  The assumption is that only nodes on the side of the
     element affect the determinant and normal on that side.
     Here we are along an "internal side" that is not a true side of the
     element and therefore these quantities can depend on all displacements
     in the element.

     Also to get moving mesh to work, the dJ's must be fixed just as the J
     was modified above.
   */
  if (DeformingMesh)
    WH(-1, "Deforming mesh not handled correctly by subelement_surfdet");

  /* NOTE: One thing that I found out the hard way is that calls to
     load_basis_functions will be ineffectual while in this state of
     lying about the element shape.  So avoid trying to do anything
     else in the lines between the lie below and the line where the
     correct shape is restored.
   */
  ielem_shape =
      ei[pg->imtrx]
          ->ielem_shape; /* store the truth so we can restore it later */
  ei[pg->imtrx]->ielem_shape = type2shape(e->ielem_type); /* lying to ei */

  get_side_info(e->ielem_type, side, &nodes_per_side, local_elem_node_id);
  surface_determinant_and_normal(
      -1, ei[pg->imtrx]->iconnect_ptr, e->num_local_nodes,
      ei[pg->imtrx]->ielem_dim - 1, side, nodes_per_side, local_elem_node_id);

  /* the lies must end */
  ei[pg->imtrx]->ielem_shape = ielem_shape; /* fixing ei */

  return (fv->sdet);
}

int num_subelement_integration_pts(Integ_Elem *e, int gpt_type, int sign) {
  int num_gpts = 0;
  if (e->num_subelements > 0) {
    int i;
    for (i = 0; i < e->num_subelements; i++) {
      num_gpts +=
          num_subelement_integration_pts(e->subelements[i], gpt_type, sign);
    }
  } else {
    if (sign == 0 || sign == e->sign) {
      if (gpt_type >= 0) /* gauss pts along side of parent element */
      {
        int iside, ip_total;
        ip_total = elem_info(NQUAD_SURF, e->ielem_type);
        for (iside = 0; iside < e->num_sides; iside++) {
          if (e->side_ids[iside] == gpt_type)
            num_gpts += ip_total;
        }
      } else if (gpt_type == -1) /* gauss pts along isosurface */
      {
        int iside, ip_total;
        if (Do_Overlap && ls->CrossMeshQuadPoints > 0)
          ip_total = ls->CrossMeshQuadPoints;
        else
          ip_total = elem_info(NQUAD_SURF, e->ielem_type);
        for (iside = 0; iside < e->num_sides; iside++) {
          if (e->bc_sides[iside])
            num_gpts += ip_total;
        }
      } else if (gpt_type == -2) /* volume gauss pts */
      {
        num_gpts = elem_info(NQUAD, e->ielem_type);
      } else {
        EH(-1, "Incorrect type of sublelement gauss points requested.");
      }
    }
  }

  return num_gpts;
}

int gather_subelement_integration_pts(Integ_Elem *e, double (*s)[DIM],
                                      double *wt, int *ip_sign, int gpt_type,
                                      int sign, int index) {
  if (e->num_subelements > 0) {
    int i;
    for (i = 0; i < e->num_subelements; i++) {
      index = gather_subelement_integration_pts(e->subelements[i], s, wt,
                                                ip_sign, gpt_type, sign, index);
    }
  } else {
    if (sign == 0 || sign == e->sign) {
      int ip_total;
      int ip;
      double yi[DIM];
      double sdet, weight;
      double detJ;

      if (gpt_type >= 0) /* gauss pts along side of parent element */
      {
        double s_i, t_i, u_i;
        int iside;
        ip_total = elem_info(NQUAD_SURF, e->ielem_type);
        for (iside = 0; iside < e->num_sides; iside++) {
          if (e->side_ids[iside] == gpt_type) {
            for (ip = 0; ip < ip_total; ip++) {
              find_surf_st(ip, e->ielem_type, iside + 1, pd->Num_Dim, yi, &s_i,
                           &t_i, &u_i);
              weight = Gq_surf_weight(ip, e->ielem_type);

              map_subelement_stu(s[index], e, yi);
              sdet = subelement_surfdet(e, yi, iside + 1, TRUE);
              wt[index] = weight * sdet;
              if (ip_sign != NULL)
                ip_sign[index] = e->sign;
              index++;
            }
          }
        }
      } else if (gpt_type == -1) /* gauss pts along isosurface */
      {
        double s_i, t_i, u_i;
        int iside;
        if (Do_Overlap && ls->CrossMeshQuadPoints > 0)
          ip_total = ls->CrossMeshQuadPoints;
        else
          ip_total = elem_info(NQUAD_SURF, e->ielem_type);
        for (iside = 0; iside < e->num_sides; iside++) {
          if (e->bc_sides[iside]) {
            for (ip = 0; ip < ip_total; ip++) {
              if (Do_Overlap && ls->CrossMeshQuadPoints > 0) {
                find_subsurf_st_wt(ip, ls->CrossMeshQuadPoints, e->ielem_type,
                                   iside + 1, pd->Num_Dim, yi, &weight);
              } else {
                find_surf_st(ip, e->ielem_type, iside + 1, pd->Num_Dim, yi,
                             &s_i, &t_i, &u_i);
                weight = Gq_surf_weight(ip, e->ielem_type);
              }
              map_subelement_stu(s[index], e, yi);
              sdet = subelement_surfdet(e, yi, iside + 1, FALSE);
              wt[index] = weight * sdet;
              if (ip_sign != NULL)
                ip_sign[index] = 0;
              index++;
            }
          }
        }
      } else if (gpt_type == -2) /* volume gauss pts */
      {
        ip_total = elem_info(NQUAD, e->ielem_type);
        for (ip = 0; ip < ip_total; ip++) {
          find_stu(ip, e->ielem_type, yi, yi + 1, yi + 2);
          map_subelement_stu(s[index], e, yi);
          detJ = subelement_detJ(e, yi, TRUE);
          if (detJ < 0.) {
            if (fabs(detJ) > 0.01)
              neg_elem_volume = TRUE;
#define DEBUG_DRN 0
#if DEBUG_DRN
            printf("for elem=%d, subelement detJ = %g\n", ei[pg->imtrx]->ielem,
                   detJ);
#endif
            detJ = 0.;
          }
          wt[index] = Gq_weight(ip, e->ielem_type) * detJ;
          if (ip_sign != NULL)
            ip_sign[index] = e->sign;
          index++;
        }
      } else {
        EH(-1, "Incorrect type of sublelement gauss points requested.");
      }
    }
  }

  return index;
}

void clear_xfem_contribution(int N) {
  int ie;

  /*	memset( &(xfem->active_vol[0]), 0, sizeof(double)*N);
          memset( &(xfem->tot_vol[0]), 0, sizeof(double)*N);

  */

  for (ie = 0; ie < N; ie++) {
    xfem->active_vol[ie] = 0.;
    xfem->tot_vol[ie] = 0.;
  }
}

void compute_xfem_contribution(int N) {
  int eqn, i, ie, ledof;
  int xfem_active, extended_dof = 0, base_interp, base_dof;
  double dV;

  for (eqn = V_FIRST; eqn < V_LAST; eqn++) {
    if (pd->e[pg->imtrx][eqn] && eqn != R_FILL && eqn != R_EXT_VELOCITY) {
      dV = bf[eqn]->detJ * fv->wt * fv->h3;
      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        xfem_dof_state(i, pd->i[pg->imtrx][eqn], ei[pg->imtrx]->ielem_shape,
                       &xfem_active, &extended_dof, &base_interp, &base_dof);
        /*fprintf(stderr,"compute_xfem_contribution: eqn=%d, i=%d, base_dof=%d,
         * extended_dof=%d\n",eqn,i,base_dof,extended_dof);*/
        if (extended_dof) {
          ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

          if (ei[pg->imtrx]->owned_ledof[ledof]) {
            if (eqn == MASS_FRACTION) {
              int ke;
              for (ke = 0; ke < upd->Max_Num_Species_Eqn; ke++) {
                ie = ei[pg->imtrx]->ieqn_ledof[ledof] + ke;

                if (ie < 0 || ie >= N) {
                  DPRINTF(stderr, "compute_xfem_contribution ie = %d\n", ie);
                  EH(-1, "compute_xfem_contrib, ie out of bounds\n");
                }

                xfem->active_vol[ie] += bf[eqn]->phi[i] * dV;
                xfem->tot_vol[ie] += dV;
              }
            } else {
              ie = ei[pg->imtrx]->ieqn_ledof[ledof];

              if (ie < 0 || ie >= N) {
                DPRINTF(stderr, "compute_xfem_contribution ie = %d\n", ie);
                EH(-1, "compute_xfem_contrib, ie out of bounds\n");
              }

              xfem->active_vol[ie] += bf[eqn]->phi[i] * dV;
              xfem->tot_vol[ie] += dV;
            }
          }
        }
      }
    }
  }
}

void check_xfem_contribution(int N, struct Aztec_Linear_Solver_System *ams,
                             double resid[], double x[], Exo_DB *exo) {
  int irow;
  double eps_standard = 1.e-4;
  double eps_diffusive = 1.e-10;
  double eps;
  int j, k;
  int j_last, ija_row;
  int eqn;
  int *ija = ams->bindx;
  double *a = ams->val;

  if (strcmp(Matrix_Format, "msr") == 0) {
    for (irow = 0; irow < N; irow++) {
      eqn = idv[pg->imtrx][irow][0];
      if (eqn == R_MASS || eqn == R_ENERGY)
        eps = eps_diffusive;
      else
        eps = eps_standard;
      if (fabs(xfem->active_vol[irow]) < eps * xfem->tot_vol[irow]) {
        j_last = ija[irow + 1] - ija[irow];
        ija_row = ija[irow];

#if 0
            /* set value to predicted one */
            a[irow] = 1.;
            resid[irow] = x[irow]-x_pred_static[irow];
#endif
#if 0
            /* set value to zero */
            a[irow] = 1.;
            resid[irow] = x[irow];
#endif
#if 0
            /* set residual to zero */
            a[irow] = 1.;
            resid[irow] = 0.;
#endif
#if 1
        /* set time derivative to zero */
        a[irow] = 1.;
        resid[irow] = x[irow] - x_old_static[irow];
#endif

        for (j = 0; j < j_last; j++) {
          k = ija_row + j;
          a[k] = 0.;
        }

        if (FALSE && xfem->active_vol[irow] != 0.) /* debugging */
        {
          DPRINTF(stderr,
                  "kill partial equation, row=%d, n = %d, active/tot=%g\n",
                  irow, idv[pg->imtrx][irow][2] + 1,
                  fabs(xfem->active_vol[irow]) / xfem->tot_vol[irow]);
        }
      }
    }
  } else if (strcmp(Matrix_Format, "epetra") == 0) {
    for (irow = 0; irow < N; irow++) {
      eqn = idv[pg->imtrx][irow][0];
      if (eqn == R_MASS || eqn == R_ENERGY) {
        eps = eps_diffusive;
      } else {
        eps = eps_standard;
      }
      if (fabs(xfem->active_vol[irow]) < eps * xfem->tot_vol[irow]) {

        EpetraSetDiagonalOnly(ams, ams->GlobalIDs[irow]);
        resid[irow] = x[irow] - x_old_static[irow];

        if (FALSE && xfem->active_vol[irow] != 0.) /* debugging */
        {
          DPRINTF(stderr,
                  "kill partial equation, row=%d, n = %d, active/tot=%g\n",
                  irow, idv[pg->imtrx][irow][2] + 1,
                  fabs(xfem->active_vol[irow]) / xfem->tot_vol[irow]);
        }
      }
    }
  } else {
    EH(-1, "Unsupported matrix format in check_xfem_contribution");
  }

  return;
}

void resolve_ls_adc_old(struct Boundary_Condition *LS_ADC, Exo_DB *exo,
                        double *x, double delta_t, int *adc_event,
                        int time_step) {
  int ssid = LS_ADC->BC_ID, ss_elem_index;
  int nodes_per_side, local_elem_node_id[MAX_NODES_PER_SIDE], ielem_type;
  int i, iss, ielem, side;
  double switch_value = 0.0;

  double P;

  if ((iss = in_list(ssid, 0, exo->num_side_sets, exo->ss_id)) == -1) {
    if (Num_Proc == 1) {
      fprintf(stderr, "Error in resolve_ls_adc.  Cannot find side set id %d",
              ssid);
      exit(-1);
    } else
      /* presumably ssid doesn't exist on this proc.  Move along. Nothing to see
       * here */
      return;
  }

  ss_elem_index = exo->ss_elem_index[iss];

  /* loop over all elements in this SS */
  for (i = 0; *adc_event == FALSE && i < exo->ss_num_sides[iss]; i++) {
    ielem = exo->ss_elem_list[ss_elem_index + i];

    side = exo->ss_side_list[ss_elem_index + i];

    ielem_type = Elem_Type(exo, ielem);
    get_side_info(ielem_type, side, &nodes_per_side, local_elem_node_id);

    if (elem_on_isosurface(ielem, x, exo, LS, 0.0)) {
      int random = rand(), rand_max = RAND_MAX;

      P = determine_adc_probability(LS_ADC, ielem, exo, x, side, nodes_per_side,
                                    local_elem_node_id, delta_t, &switch_value);

      if (P != 0.0 && (random <= rand_max * P)) {
        /*apply_adc_to_elem ( exo, x, ielem, side, nodes_per_side,
         * local_elem_node_id, 0,  switch_value );*/
        *adc_event = TRUE;
      }
    }
  }
  apply_adc_to_ss(exo, x, iss, switch_value);
  return;
}

static double determine_adc_probability(struct Boundary_Condition *ls_adc,
                                        int ielem, Exo_DB *exo, double *x,
                                        int side, int nodes_per_side,
                                        int *local_elem_node_id, double delta_t,
                                        double *switch_value) {
  double P = 1.0;
  double avg_cos, area, surf_area;
  double xi[DIM];
  double capture_angle = ls_adc->BC_Data_Float[0];
  double capture_distance = ls_adc->BC_Data_Float[1];
  double capture_rate = ls_adc->BC_Data_Float[2];
  double average_snormal[DIM];
  double value0;
  int a, i, ip, ip_total;

  load_elem_dofptr(ielem, exo, x, x, x, x, 0);

  value0 = *esp->F[local_elem_node_id[0]];

  for (i = 1; P != 0.0 && i < nodes_per_side; i++) {
    int ln;

    ln = local_elem_node_id[i];
    if (fabs(*esp->F[ln]) < 1.e-12)
      P = 0.0;
    if (*esp->F[ln] * value0 < 0.0)
      P = 0.0;
  }

  if (P == 0.0)
    return (P);

  ip_total = elem_info(NQUAD_SURF, ei[pg->imtrx]->ielem_type);

  average_snormal[0] = 0.0;
  average_snormal[1] = 0.0;
  average_snormal[2] = 0.0;

  for (ip = 0, surf_area = 0.0; ip < ip_total; ip++) {
    double s, t, u;
    find_surf_st(ip, ei[pg->imtrx]->ielem_type, side, ei[pg->imtrx]->ielem_dim,
                 xi, &s, &t, &u);
    fv->wt = Gq_surf_weight(ip, ei[pg->imtrx]->ielem_type);

    load_basis_functions(xi, bfd);

    beer_belly();

    load_fv();

    load_bf_grad();

    load_fv_grads();

    surface_determinant_and_normal(
        ielem, ei[pg->imtrx]->iconnect_ptr, ei[pg->imtrx]->num_local_nodes,
        ei[pg->imtrx]->ielem_dim - 1, side, nodes_per_side, local_elem_node_id);

    for (a = 0; a < pd->Num_Dim; a++)
      average_snormal[a] += fv->snormal[a] / ip_total;

    surf_area += fv->wt * fv->h3 * fv->sdet;
  }

  ip_total = elem_info(NQUAD, ei[pg->imtrx]->ielem_type);

  for (ip = 0, avg_cos = 0.0, area = 0.0; ip < ip_total; ip++) {
    double sum;

    find_stu(ip, ei[pg->imtrx]->ielem_type, &xi[0], &xi[1], &xi[2]);

    fv->wt = Gq_weight(ip, ei[pg->imtrx]->ielem_type);

    load_basis_functions(xi, bfd);

    beer_belly();

    load_fv();

    load_bf_grad();

    load_fv_grads();

    load_lsi(ls->Length_Scale);

    for (a = 0, sum = 0.0; a < pd->Num_Dim; a++) {
      sum += average_snormal[a] * lsi->normal[a];
    }

    avg_cos += fv->wt * sum * fv->h3 * bf[pd->ShapeVar]->detJ;
    area += fv->wt * fv->h3 * bf[pd->ShapeVar]->detJ;
  }

  if (area < 1.e-12)
    EH(-1, "Error zero area element detected in determine_adc_probability.\n");

  avg_cos /= area;

  if (fabs(avg_cos) < fabs(cos(M_PI * capture_angle / 180)))
    P = 0.0;

  if (P != 0.0) {
    value0 =
        determine_nearest_distance(exo, x, nodes_per_side, local_elem_node_id);

    *switch_value = fabs(value0);

    if (fabs(value0) < capture_distance) {
      P = capture_rate * surf_area * delta_t;
      P = P >= 1.0 ? 1.0 : P;
    } else {
      P = capture_rate * sqrt(area) * delta_t *
          exp(1.0 - pow(value0 / capture_distance, 2.0));
    }
    fprintf(stderr, "ADC event possible at elem %d with probability %g\n",
            ielem, P);
  }

  return (P);
}

static void apply_adc_to_ss(Exo_DB *exo, double *x, int iss,
                            double switch_value)

{
  int i, ielem, side, ss_elem_index;
  int nodes_per_side, local_elem_node_id[MAX_NODES_PER_SIDE], ielem_type;
  double start_sign = 123.0;
  int *apply_to_side;

  ss_elem_index = exo->ss_elem_index[iss];
  apply_to_side = alloc_int_1(exo->ss_num_sides[iss], FALSE);

  for (i = 0; i < exo->ss_num_sides[iss]; i++) {
    ielem = exo->ss_elem_list[ss_elem_index + i];

    if (elem_on_isosurface(ielem, x, exo, LS, 0.0)) {
      apply_to_side[i] = TRUE;
    }
  }

  for (i = 0; i < exo->ss_num_sides[iss]; i++) {
    ielem = exo->ss_elem_list[ss_elem_index + i];

    side = exo->ss_side_list[ss_elem_index + i];

    ielem_type = Elem_Type(exo, ielem);
    get_side_info(ielem_type, side, &nodes_per_side, local_elem_node_id);

    if (apply_to_side[i] == TRUE) {
      load_elem_dofptr(ielem, exo, x, x, x, x, 0);

      if (start_sign == 123.0)
        start_sign = sign_of(*esp->F[local_elem_node_id[0]]);

      apply_adc_to_elem(exo, x, ielem, side, nodes_per_side, local_elem_node_id,
                        start_sign, switch_value);
    }
  }
  safe_free((void *)apply_to_side);
  return;
}

static void apply_adc_to_elem(Exo_DB *exo, double *x, int ielem, int side,
                              int nodes_per_side, int *local_elem_node_id,
                              double start_sign, double switch_value)

{
  int i, ln, I;
  double r[DIM];
  struct LS_Surf_List *list;

  if (start_sign == 0.0)
    start_sign = sign_of(*esp->F[local_elem_node_id[0]]);

  list = create_surf_list();

  for (i = 0; i < nodes_per_side; i++) {
    struct LS_Surf *s;

    ln = local_elem_node_id[i];
    I = ei[pg->imtrx]->gnn_list[LS][ln];

    r[0] = Coor[0][I];
    r[1] = Coor[1][I];
    if (pd->Num_Dim == 3)
      r[2] = Coor[2][I];

    s = create_surf_point(r, ei[pg->imtrx]->ielem, r, FALSE);

    if (unique_surf(list, s)) {
      append_surf(list, s);
    } else {
      safe_free(s);
    }
  }

  switch_value *= (-start_sign) * 1.e-2;

  for (i = 0; i < nodes_per_side; i++) {
    ln = local_elem_node_id[i];

    if ((*esp->F[ln]) * start_sign > 0.0) {
      *esp->F[ln] = switch_value;
    }
  }

  for (i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
    struct LS_Surf *s;
    if (in_list(i, 0, nodes_per_side, local_elem_node_id) == -1) {
      I = ei[pg->imtrx]->gnn_list[LS][i];
      r[0] = Coor[0][I];
      r[1] = Coor[1][I];
      if (pd->Num_Dim == 3)
        r[2] = Coor[2][I];

      s = closest_surf(list, x, exo, r);
      if (ei[pg->imtrx]->ln_to_first_dof[LS][i] != -1) {
        *esp->F[i] = (-start_sign) * fabs(s->closest_point->distance);
      }
    }
  }
  DPRINTF(stderr, "ADC event applied to element %d.\n", ielem);

  free_surf_list(&list);

  return;
}

static double determine_nearest_distance(Exo_DB *exo, double *x,
                                         int nodes_per_side,
                                         int *local_elem_node_id) {
  int i, ln, I;
  struct LS_Surf_List *list;
  double nearest_distance;
  double r[DIM] = {0., 0., 0.};

  list = create_surf_list();

  find_intersections(list, LS, 0.0, exo);

  for (i = 0, nearest_distance = 0.0; i < nodes_per_side; i++) {
    struct LS_Surf *s;
    double dist;

    ln = local_elem_node_id[i];
    I = ei[pg->imtrx]->gnn_list[LS][ln];

    r[0] = Coor[0][I];
    r[1] = Coor[1][I];
    if (pd->Num_Dim == 3)
      r[2] = Coor[2][I];

    s = closest_surf(list, x, exo, r);

    dist = fabs(s->closest_point->distance);

    nearest_distance += dist / nodes_per_side;
  }

  free_surf_list(&list);

  return (nearest_distance);
}

struct LS_Surf *resolve_ls_adc(struct LS_Surf_List *Surf_list,
                               struct Boundary_Condition *ls_adc, Exo_DB *exo,
                               double *x, double delta_t, int *adc_event,
                               int time_step) {
  struct LS_Surf *II = Surf_list->start;

  double capture_angle = ls_adc->BC_Data_Float[0];
  double capture_distance = ls_adc->BC_Data_Float[1];
  double capture_rate = ls_adc->BC_Data_Float[2];
  double P_sum = 0.0;

  *adc_event = FALSE;

  /*printf("Testing for contact the following locations :\n");*/

  while ((*adc_event == FALSE) && II != NULL) {
    struct LS_Surf_Point_Data *dd = (struct LS_Surf_Point_Data *)II->data;

    setup_shop_at_point(dd->elem, dd->xi, exo);

    load_lsi(ls->Length_Scale);

    if (check_alignment(capture_angle)) {
      double phase_func = fabs(fv->pF[0]);
      double contact_area;
      struct LS_Surf *nearest_other;

      double P = 0.0;
      int random = rand(), rand_max = RAND_MAX;

      /*printf("\t\t(%lf, %lf)", dd->x[0], dd->x[1] );*/

      nearest_other = closest_other_surf(Surf_list, x, exo, II);

      contact_area = nearest_other->closest_point->distance;

      if (pd->Num_Dim == 3)
        contact_area *= contact_area * M_PI / 4.0;

      if (phase_func < capture_distance) {
        P = capture_rate * contact_area * delta_t;
      } else {
        P = capture_rate * contact_area * delta_t *
            exp(1.0 - pow(phase_func / capture_distance, 3.0));
      }

      P = P >= 1.0 ? 1.0 : P;

      P_sum += P;

      if (P != 0.0 && (random <= rand_max * P)) {
        int I_adc;
        double distance, p[DIM], val0;

        val0 = find_adc_node(ls_adc->BC_ID, x, exo, II, &distance, &I_adc);
        printf("ADC event detected at point (%lf, %lf) with probability %lf "
               "and contact area %lf\n",
               dd->x[0], dd->x[1], P, contact_area);
        printf("ADC contact event to be resolved at node %d, distance %lf with "
               "value %lf\n",
               I_adc, distance, val0);

        if (I_adc != -1) {
          retrieve_node_coordinates(I_adc, x, p, NULL);

          apply_adc_function(x, exo, p, distance, 1.05 * val0);

          *adc_event = TRUE;
        }
      }
    }

    II = II->next;
  }

  if (P_sum != 0.0) {
    DPRINTF(stderr,
            "\tCumulative contact probability for this time step \t ( %lf ) \n",
            P_sum);
  }
  return (II);
}

static int check_alignment(double capture_angle) {
  int dim = pd->Num_Dim;

  double nw[DIM] = {0., 0., 0.};

  double cosine = 0.;

  int is_aligned = FALSE;

  /* So the big assumption here is that the gradient of the phase function field
   * AT THE INTERFACE CURVE can be used as a good approximation to the wall
   * normal This assumption will improve the nearer the level set interface is
   * to the wall itself o'course
   */

  nw[0] = -fv->grad_pF[0][0];
  nw[1] = -fv->grad_pF[0][1];
  if (dim == 3)
    nw[2] = -fv->grad_pF[0][2];

  normalize_really_simple_vector(nw, dim);

  cosine = nw[0] * lsi->normal[0];
  cosine += nw[1] * lsi->normal[1];
  if (dim == 3)
    cosine += nw[2] * lsi->normal[2];

  if (fabs(cosine) >= fabs(cos(M_PI * capture_angle / 180)))
    is_aligned = TRUE;

  return (is_aligned);
}

static struct LS_Surf *closest_other_surf(struct LS_Surf_List *Surf_list,
                                          double *x, Exo_DB *exo,
                                          struct LS_Surf *this_surf) {

  double min_distance = 1.e+30;
  struct LS_Surf *that_surf;
  struct LS_Surf *closest_other = NULL;
  struct LS_Surf_Point_Data *dd = (struct LS_Surf_Point_Data *)this_surf->data;

  that_surf = Surf_list->start;

  while (that_surf != NULL) {
    if (that_surf != this_surf) {
      find_surf_closest_point(that_surf, x, exo, dd->x);

      if (that_surf->closest_point->distance < min_distance) {
        min_distance = that_surf->closest_point->distance;
        closest_other = that_surf;
      }
    }

    that_surf = that_surf->next;
  }

  return (closest_other);
}

static double find_adc_node(int ns_id, double *x, Exo_DB *exo,
                            struct LS_Surf *II, double *dist, int *I_adc) {
  int i;

  int ns, ns_num_nodes;
  int *ns_node_list;
  double p[DIM];

  int ie_adc = -1;

  ns = 0;

  *I_adc = -1;
  *dist = 1.e10;

  while (exo->ns_id[ns] != ns_id && ns++ < exo->num_node_sets)
    ;

  if (ns >= exo->num_node_sets) {
    if (Num_Proc == 1)
      EH(-1, "In find_adc_node:  Can't locate nodeset id \n");
    else
      /* Parallel problem:  likely that this sideset isn't on this proc */
      return (-1.e33);
  }

  ns_num_nodes = exo->ns_num_nodes[ns];

  ns_node_list = exo->ns_node_list + exo->ns_node_index[ns];

  for (i = 0; i < ns_num_nodes; i++) {
    int I, ie;
    I = ns_node_list[i];

    ie = Index_Solution(I, LS, 0, 0, -2, pg->imtrx);

    if (ie != -1) {

      retrieve_node_coordinates(I, x, p, NULL);

      find_surf_closest_point(II, x, exo, p);

      if (II->closest_point->distance < *dist) {
        *dist = II->closest_point->distance;
        *I_adc = I;
        ie_adc = ie;
      }
    }
  }

  if (ie_adc == -1)
    EH(-1, " Error in find_adc_node:  Can't locate level set unknown\n");

  return (x[ie_adc]);
}

static void apply_adc_function(double *x, Exo_DB *exo, double *p0, double d,
                               double val0) {

  int I, ie;
  double p[DIM], func;

  for (I = 0; I < exo->num_nodes; I++) {
    ie = Index_Solution(I, LS, 0, 0, -2, pg->imtrx);

    if (ie != -1) {
      double ll;

      retrieve_node_coordinates(I, x, p, NULL);

      ll = (p[0] - p0[0]) * (p[0] - p0[0]);
      ll += (p[1] - p0[1]) * (p[1] - p0[1]);

      if (pd->Num_Dim == 3)
        ll += (p[2] - p0[2]) * (p[2] - p0[2]);

      ll = sqrt(ll);

      func = val0 * exp(-pow((ll / d), 2.0));

      if (fabs(func) > 1.e-15) {
        x[ie] = x[ie] - func;
      }
    }
  }
  return;
}

static int is_LS_spurious(Exo_DB *exo, double *x, int this_node, double sign,
                          double *avg_value) {
  int j, ie;
  int node;
  int node_ptr = exo->node_node_pntr[this_node];
  int num_nodes = exo->node_node_pntr[this_node + 1] - node_ptr;

  int is_spurious = TRUE;

  double value = 0.0;

  *avg_value = 0.0;

  for (j = 0; is_spurious && j < num_nodes; j++) {
    node = exo->node_node_list[node_ptr + j];

    ie = Index_Solution(node, LS, 0, 0, -2, pg->imtrx);

    if (ie != -1 && node != this_node) {
      value = x[ie];

      *avg_value += value / (num_nodes - 1);

      if (sign * value >= 0.0) {
        is_spurious = FALSE;
        *avg_value = 1.e30;
      }
    }
  }
  return (is_spurious);
}
static void purge_spurious_LS(double *x, Exo_DB *exo, int num_total_nodes) {
  int I, ie;
  double sign, avg_value;

  for (I = 0; I < num_total_nodes; I++) {

    ie = Index_Solution(I, LS, 0, 0, -2, pg->imtrx);

    if (ie != -1) {
      sign = x[ie] >= 0 ? 1.0 : -1.0;

      if (is_LS_spurious(exo, x, I, sign, &avg_value)) {
        x[ie] = avg_value;
      }
    }
  }
  return;
}
