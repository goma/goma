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
 *$Id: mm_fill.c,v 5.36 2010-07-21 16:39:26 hkmoffa Exp $
 */

/* Standard include files */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* GOMA include files */

#include "ac_stability.h"
#include "ac_stability_util.h"
#include "bc/rotate.h"
#include "bc/rotate_coordinates.h"
#include "bc_colloc.h"
#include "bc_contact.h"
#include "bc_curve.h"
#include "bc_dirich.h"
#include "bc_integ.h"
#include "bc_special.h"
#include "dpi.h"
#include "el_elm.h"
#include "el_elm_info.h"
#include "el_geom.h"
#include "exo_struct.h"
#include "load_field_variables.h"
#include "md_timer.h"
#include "mm_as.h"
#include "mm_as_alloc.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "mm_fill.h"
#include "mm_fill_aux.h"
#include "mm_fill_common.h"
#include "mm_fill_em.h"
#include "mm_fill_fill.h"
#include "mm_fill_ls.h"
#include "mm_fill_population.h"
#include "mm_fill_porous.h"
#include "mm_fill_potential.h"
#include "mm_fill_pthings.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_rs.h"
#include "mm_fill_shell.h"
#include "mm_shell_util.h"

#include "dpi.h"
#include "exo_struct.h"

#include "mm_fill_momentum.h"
#include "mm_fill_species.h"
#include "mm_fill_stress.h"
#include "mm_fill_terms.h"
#include "mm_fill_turbulent.h"
#include "mm_fill_util.h"
#include "mm_flux.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "mm_qp_storage.h"
#include "mm_qtensor_model.h"
#include "mm_shell_util.h"
#include "mm_unknown_map.h"
#include "mpi.h"
#include "rd_mesh.h"
#include "rf_allo.h"
#include "rf_bc.h"
#include "rf_bc_const.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_io.h"
#include "rf_masks.h"
#include "rf_mp.h"
#include "rf_node_const.h"
#include "rf_solver.h"
#include "rf_solver_const.h"
#include "sl_epetra_util.h"
#include "sl_util.h"
#include "sl_util_structs.h"
#include "std.h"
#ifdef GOMA_ENABLE_PETSC
#include "sl_petsc.h"
#endif
#include "density.h"
#include "mm_fill_continuity.h"
#include "mm_fill_energy.h"
#include "mm_fill_split.h"
#include "sl_petsc_complex.h"
#include "stdbool.h"
#include "wr_side_data.h"

#define GOMA_MM_FILL_C

/*
 * Global variables defined here. Declared frequently via rf_bc.h
 */

/*
 * Demarcate the CPU consumption of matrix_fill().
 * Play some tricks to accumulate a reasonable depiction of the
 * total assembly time for all the elements and not just one element.
 */

double mm_fill_total;
extern int PRS_mat_ielem;

static void load_lec(Exo_DB *, /* Exodus database pointer */
                     int,      /* element number we are working on */
                     struct GomaLinearSolverData *,
                     double[],  /* Solution vector */
                     double[],  /* Residual vector */
                     double *); /* element stiffness Matrix for frontal solver*/
static void zero_lec(void);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/*************************************************************************
 *
 * matrix_fill_full:
 *
 *
 *  Construct the Jacobian and Residual. This is a wrapper around the
 *  element level fill routine, matrix_fill().
 *
 * Return:
 *  -1 : A negative element volume was encountered somewhere in the
 *       domain. The Matrix is not filled till completion
 *   0 : Successful completion.
 *************************************************************************/

int matrix_fill_full(struct GomaLinearSolverData *ams,
                     double x[],
                     double resid_vector[],
                     double x_old[],
                     double x_older[],
                     double xdot[],
                     double xdot_old[],
                     double x_update[],
                     double *ptr_delta_t,
                     double *ptr_theta,
                     struct elem_side_bc_struct *first_elem_side_BC_array[],
                     double *ptr_time_value,
                     Exo_DB *exo,
                     Dpi *dpi,
                     int *ptr_num_total_nodes,
                     dbl *ptr_h_elem_avg,
                     dbl *ptr_U_norm,
                     dbl *estifm) {
  int ielem = 0, e_start = 0, e_end = 0, ebn = 0;
  char yo[] = "matrix_fill_full";
  int err, err_global;

#define debug_subelement_decomposition 0
#if debug_subelement_decomposition
  if (ls != NULL && ls->SubElemIntegration)
    subelement_mesh_output(x, exo);
  if (ls != NULL && ls->SubElemIntegration) {
    double surf_min[3];
    double surf_max[3];
    ls_surface_extents(x, exo, surf_min, surf_max);
    fprintf(stderr, "Level set surface extents: %g <= x <= %g, %g <= y <= %g\n", surf_min[0],
            surf_max[0], surf_min[1], surf_max[1]);
  }
#endif

  if (xfem != NULL)
    clear_xfem_contribution(ams->npu);

#define DEBUG_LS_INTEGRATION 0
#if DEBUG_LS_INTEGRATION
  if (ls != NULL) {
    ls->Neg_Vol = 0.;
    ls->Pos_Vol = 0.;
    ls->Surface_Area = 0.;
  }
#endif

  /*
   * One could do other loops here. For example, loop over global
   * nodes, precalculating likely quantities.
   */

  /*
   * Loop over all of the elements one a time. Obtain their
   * element contributions to the global matrix
   *
   * -> This would be a good spot to loop over each element block instead
   *    Then allocate memory here for common data elements
   */
  neg_elem_volume = FALSE;
  neg_lub_height = FALSE;
  zero_detJ = FALSE;

  e_start = exo->eb_ptr[0];
  e_end = exo->eb_ptr[exo->num_elem_blocks];

  for (ielem = e_start, ebn = 0; ielem < e_end && !neg_elem_volume && !neg_lub_height && !zero_detJ;
       ielem++) {

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

    err = matrix_fill(ams, x, resid_vector, x_old, x_older, xdot, xdot_old, x_update, ptr_delta_t,
                      ptr_theta, first_elem_side_BC_array, ptr_time_value, exo, dpi, &ielem,
                      ptr_num_total_nodes, ptr_h_elem_avg, ptr_U_norm, estifm, 0);

    if (err)
      break;

    if (neg_elem_volume) {
      log_msg("Negative elem det J in element (%d)", ielem + 1);
      if (ls != NULL && ls->SubElemIntegration)
        subelement_mesh_output(x, exo);
    }

    if (neg_lub_height) {
      log_msg("Negative lubrication height in element (%d)", ielem + 1);
    }

    if (zero_detJ) {
      log_msg("Zero determinant of Jacobian of transformation (%d)", ielem + 1);
    }
  }

  /*
   * Free memory allocated above
   */
  global_qp_storage_destroy();

  /*
   * Now coordinate the processors so that they all know about a negative or zero
   * volume in an element and negative lubrication height
   */
#ifdef PARALLEL
  MPI_Allreduce(&neg_elem_volume, &neg_elem_volume_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  neg_elem_volume = neg_elem_volume_global;

  MPI_Allreduce(&neg_lub_height, &neg_lub_height_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  neg_lub_height = neg_lub_height_global;

  MPI_Allreduce(&zero_detJ, &zero_detJ_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  zero_detJ = zero_detJ_global;

  MPI_Allreduce(&err, &err_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  err = err_global;

#endif

  if (err)
    return -1;

  /*
   * Do common utility tasks, such as printing out the matrix when in
   * debug mode
   */

  if (xfem != NULL)
    check_xfem_contribution(ams->npu, ams, resid_vector, x, exo);

#if DEBUG_LS_INTEGRATION
  if (ls != NULL) {
    DPRINTF(stderr,
            "Debugging integration, Time,Neg_Vol,Pos_Vol,Area = %20.12e %20.12e %20.12e %20.12e\n",
            *ptr_time_value, ls->Neg_Vol, ls->Pos_Vol, ls->Surface_Area);
  }
#endif

  if (neg_elem_volume)
    return -1;
  if (neg_lub_height)
    return -1;
  if (zero_detJ)
    return -1;

  return 0;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int matrix_fill(struct GomaLinearSolverData *ams,
                double x[],            /* Solution vector */
                double resid_vector[], /* Residual vector */
                double x_old[],        /* Solution vector , previous last time step */
                double x_older[],      /* Solution vector , previous prev time step */
                double xdot[],         /* xdot of current solution                  */
                double xdot_old[],     /* xdot_old of current solution              */
                double x_update[],     /* last update vector */
                double *ptr_delta_t,   /* current time step size                    */
                double *ptr_theta,     /* parameter to vary time integration from
                                        * explicit (theta = 1) to implicit
                                        * (theta = 0) */
                struct elem_side_bc_struct *first_elem_side_BC_array[],
                /* This is an array of pointers to the first
                   surface integral defined for each element.
                   It has a length equal to the total number
                   of elements defined on current processor */
                double *ptr_time_value,
                Exo_DB *exo,              /* ptr to EXODUS II finite element mesh db */
                Dpi *dpi,                 /* ptr to distributed processing info */
                int *ptr_ielem,           /* element number */
                int *ptr_num_total_nodes, /* Number of nodes that each processor is
                                             responsible for                              */
                dbl *ptr_h_elem_avg,      /* global average element size for PSPG */
                dbl *ptr_U_norm,          /* global average velocity for PSPG */
                dbl *estifm,              /* element stiffness Matrix for frontal solver*/
                int zeroCA)

/******************************************************************************
Function which fills the FEM stiffness matrices and the right-hand side.

Author:          Harry Moffat, Scott Hutchinson (1421)
Date:            24 November 1992
Revised:         Wed Mar  2 14:05:04 MST 1994 pasacki@sandia.gov
Revised:         9/15/94 RRR
Revised:         Summer 1998, SY Tam (UNM)
******************************************************************************/

{
  double delta_t, theta, time_value, h_elem_avg; /*see arg list */
  int ielem, num_total_nodes;

#if 0
  int status = 0;		/* variable describing status of the matrix */
				/* calculation this is set to -1 if further */
				/* iteration is not recommended*/
#endif

  static int CA_id[MAX_CA];     /*  array of CA conditions */
  static int CA_fselem[MAX_CA]; /*  array of CA free surface elements  */
  static int CA_sselem[MAX_CA]; /*  array of CA solid surface elements */
  static int CA_proc[MAX_CA];   /*  Processor which has each CA */

  int mn;                     /* material block counter */
  int err;                    /* temp variable to hold diagnostic flags.      */
  int ip;                     /* ip is the local quadrature point index       */
  int ip_total;               /* ip_total is the total number of volume
                                 quadrature points in the element             */
  int j;                      /* local index loop counter 	             */
  int i;                      /* Index for the local node number - row        */
  int b, eqn = -1, e, lm_dof; /*convenient params for neatness */
  int I;                      /* Index for global node number - row        */
  int ielem_type = 0;         /* Element type of the current element          */
  int ielem_dim;              /* Element physical dimension                   */
  int num_local_nodes;        /* Number of local basis functions in the
                                 current element                              */
  int iconnect_ptr;           /* Pointer to the beginning of the connectivity
                                 list for the current element                 */
  int ibc;                    /* Index for the boundary condition 	      */
  int var = -1;               /* variable name (TEMPERATURE, etc) */

  double s, t, u; /* Gaussian-quadrature point locations          */

  int call_int, call_col, call_contact, call_shell_grad, call_sharp_int;
  /* int call_special; */
  /* flags for calling boundary
     condition routines */
  int call_rotate;
  int rotate_mesh, rotate_momentum;
  int assemble_rs, make_trivial = 0; /* Flags for Eulerian solid mechanics
                                    and level-set */
  int bct, mode;

  /* ___________________________________________________________________________*/
  /* ___________________________________________________________________________*/

  double xi[DIM]; /* Local element coordinates of Gauss point. */

  double wt = 0.0; /* Quadrature weights
                    units - ergs/(sec*cm*K) = g*cm/(sec^3*K)     */

  double ls_F[MDE];                  /* local copy for adaptive weights  */
  double ad_wtpos[10], ad_wtneg[10]; /*adaptive integration weights  */

#if 0
  dbl h[DIM];                 /* element size variable for PSPG */
  dbl hh[DIM][DIM];
  dbl dh_dxnode[DIM][MDE];

  dbl hsquared[DIM];          /* square of the element size variable for SUPG */
  dbl hhv[DIM][DIM];
  dbl dhv_dxnode[DIM][MDE];
  dbl v_avg[DIM];              /* average element velocity for SUPG */
  dbl dv_dnode[DIM][MDE];      /* nodal derivatives of average element velocity for SUPG */
#endif

  struct Petrov_Galerkin_Data pg_data;

  struct Porous_Media_Terms pm_terms; /*Needed up here for Hysteresis switching criterion*/

  struct elem_side_bc_struct *elem_side_bc;
  /* Pointer to an element side boundary condition
     structure				      */
  struct elem_edge_bc_struct *elem_edge_bc;
  /* Pointer to an element edge boundary condition
     structure				      */
  int bc_input_id; /* Input ID of the surf_int boundary condition  */
  int id_side;     /* side counter for discontinuous Galerkin method */

  /* List to keep track of nodes at which residual and
     Jacobian corrections have been made           */
  /* for mesh equations */
  /* List to keep track of nodes at which residual and
     Jacobian corrections have been made           */
  /* for mesh equations */
  int local_node_list_fs[MDE]; /* list to keep track of nodes at which solid contributions
                                  have been transfered to liquid (fluid-solid boundaries)*/

  int discontinuous_mass;   /* flag that tells you if you are doing Discontinuous Galerkin
                               for the species equations */
  int discontinuous_stress; /* flag that tells you if you are doing Discontinuous Galerkin
                               for the species equations */
  int ielem_type_mass = -1; /* flag to give discontinuous interpolation type */

  int pspg_local = 0;
  int PorousShellOn = 0; /* flag for porous shell equations */
  bool owner = TRUE;

  NODE_INFO_STRUCT *node;
  SGRID *element_search_grid = NULL;

  struct Level_Set_Data *ls_old;

  static double mm_fill_start, mm_fill_end; /* Count CPU time this call. */

  static char yo[] = "matrix_fill"; /* My name to take blame... */

  int *pde;

  /*
   * BEGINNING OF EXECUTABLE STATEMENTS
   */

  mm_fill_start = ut();

  /*
   * Unpack some pointers that were added to help portability of front
   */
  delta_t = *ptr_delta_t;
  theta = *ptr_theta;
  time_value = *ptr_time_value;
  ielem = *ptr_ielem;
  num_total_nodes = *ptr_num_total_nodes;
  h_elem_avg = pg_data.h_elem_avg = *ptr_h_elem_avg;
  pg_data.U_norm = *ptr_U_norm;

  if (Debug_Flag > 1) {
    P0PRINTF("%s: begins\n", yo);
  }

  /******************************************************************************/
  /*                                BLOCK 1                                     */
  /*          LOOP OVER THE ELEMENTS DEFINED ON THE CURRENT PROCESSOR           */
  /*          INITIALIZATION THAT IS DEPENDENT ON THE ELEMENT NUMBER            */
  /******************************************************************************/

  /*
   * For each variable there are generally different degrees of
   * freedom that they and their equations contribute to.
   *
   * For this element, load up arrays that tell, for each variable,
   *
   *	(i) how many degrees of freedom they contribute towords
   *	(ii) the local node number associated with each degree of
   *	     freedom
   *	(iii) pointers in the "esp" structure that tell where
   *	      things are located in the global scheme...
   *		(a) nodal unknowns in this element...
   *		(b) Residual equations receiving contributions in
   *		    this element.
   *		(c) where the Jacobian entries go in the global
   *		    "a" matrix in its MSR format...
   */

  /*   Initialize CA element arrays. These arrays are */
  /*   are used to track the connectivity of elements */
  /*   around contact lines, so that a variable wall  */
  /*   normal can be used when the wall and free surface */
  /*   are in separate elements.                      */
  /*   PRS NOTE:  This is actually a very inefficient */
  /*   initialization, as these arrays could be set up*/
  /*   in one pre-processing step. However, that step */
  /*   requires an element sweep without assembly.    */
  /*   Could and should be done.                      */

  /*
   * Hmmm...the elem_order_map[] array is not always guaranteed to exist.
   * Need to make provision for this, say in rd_exo.c. (and don't forget
   * to free it up at the end, too!)
   *
   * There's a boolean exo->elem_order_map_exists that should be set and
   * checked accordingly to avoid problems in case exo->elem_order_map is NULL
   *
   */

  if (exo->ns_node_len > 0 &&
      ((zeroCA == 1) || ((Linear_Solver != FRONT && ielem == exo->eb_ptr[0]) ||
                         (Linear_Solver == FRONT && ielem == exo->elem_order_map[0] - 1)))) {
    int nsp, nspk, count = -1;
    memset(CA_fselem, -1, sizeof(int) * MAX_CA);
    memset(CA_sselem, -1, sizeof(int) * MAX_CA);
    memset(CA_id, -1, sizeof(int) * MAX_CA);
    memset(CA_proc, -1, sizeof(int) * MAX_CA);
    for (j = 0; j < Num_BC; j++) {
      switch (BC_Types[j].BC_Name) {
      case CA_BC:
      case CA_MOMENTUM_BC:
      case VELO_THETA_HOFFMAN_BC:
      case VELO_THETA_TPL_BC:
      case VELO_THETA_COX_BC:
      case VELO_THETA_SHIK_BC: {
        nsp = match_nsid(BC_Types[j].BC_ID);
        if (nsp != -1) {
          int n_nodes = exo->ns_num_nodes[nsp];
          if (n_nodes > 0) {
            nspk = exo->ns_node_list[exo->ns_node_index[nsp]];
            if (Nodes[nspk]->Proc == ProcID) {
              count++;
              CA_proc[count] = ProcID;
            }
          }
        }
      } break;
      }
    }

    /*
     * Initialize the accumulated CPU time for assembly...
     */
    mm_fill_total = 0;
  }

  /*
   *  load_elem_dofptr:
   *       This routine loads a lot of useful pointers concerning the
   *       current element refering to materials number:
   *             mp = pointer to the material property structure
   *             pd = problem description structure
   *             cr = cr structure
   *             elc
   *             vn = vn structure
   *             ve = ve structure
   *       It also fills in these element based structures:
   *             ei = Element Indeces structures
   *
   * First, make all the pointers for different variables point to
   * a bona fide zero value so that references to undefined variables
   * give a zero by default...
   */
  err = load_elem_dofptr(ielem, exo, x, x_old, xdot, xdot_old, 0);
  GOMA_EH(err, "load_elem_dofptr");

  err = bf_mp_init(pd);
  mn = ei[pg->imtrx]->mn;
  pde = (int *)pd->e[pg->imtrx];

  for (mode = 0; mode < vn->modes; mode++) {
    ve[mode] = ve_glob[mn][mode];
  }

  /* discontinuous Galerkin information */
  /* element type, assumed to be globally consistent for now */

  /* need to take this section out and up above.  Inefficient now as it assumes
     every element is different */

  discontinuous_mass = 0;

  if (pd->i[pg->imtrx][MASS_FRACTION] == I_P1) {
    if (pd->Num_Dim == 2)
      ielem_type_mass = P1_QUAD;
    if (pd->Num_Dim == 3)
      ielem_type_mass = P1_HEX;
    discontinuous_mass = 1;
  } else if (pd->i[pg->imtrx][MASS_FRACTION] == I_P0) {
    if (pd->Num_Dim == 2)
      ielem_type_mass = P0_QUAD;
    if (pd->Num_Dim == 3)
      ielem_type_mass = P0_HEX;
    discontinuous_mass = 1;
  } else if (pd->i[pg->imtrx][MASS_FRACTION] == I_PQ1) {
    if (pd->Num_Dim == 2)
      ielem_type_mass = BILINEAR_QUAD;
    if (pd->Num_Dim == 3)
      GOMA_EH(GOMA_ERROR, "Sorry PQ1 interpolation has not been implemented in 3D yet.");
    discontinuous_mass = 1;
  } else if (pd->i[pg->imtrx][MASS_FRACTION] == I_PQ2) {
    if (pd->Num_Dim == 2)
      ielem_type_mass = BIQUAD_QUAD;
    if (pd->Num_Dim == 3)
      GOMA_EH(GOMA_ERROR, "Sorry PQ2 interpolation has not been implemented in 3D yet.");
    discontinuous_mass = 1;
  } else {
    ielem_type_mass = ielem_type;
  }

  discontinuous_stress = 0;

  if (pd->i[pg->imtrx][POLYMER_STRESS11] == I_P1) {
    discontinuous_stress = 1;
  } else if (pd->i[pg->imtrx][POLYMER_STRESS11] == I_P0) {
    discontinuous_stress = 1;
  } else if (pd->i[pg->imtrx][POLYMER_STRESS11] == I_PQ1) {
    if (pd->Num_Dim == 3)
      GOMA_EH(GOMA_ERROR, "Sorry PQ1 interpolation has not been implemented in 3D yet.");
    discontinuous_stress = 1;
  } else if (pd->i[pg->imtrx][POLYMER_STRESS11] == I_PQ2) {
    if (pd->Num_Dim == 3)
      GOMA_EH(GOMA_ERROR, "Sorry PQ2 interpolation has not been implemented in 3D yet.");
    discontinuous_stress = 1;
  }

  ielem_type = ei[pg->imtrx]->ielem_type; /* element type */

  num_local_nodes = ei[pg->imtrx]->num_local_nodes; /* number of local  basis functions */

  ielem_dim = ei[pg->imtrx]->ielem_dim; /* physical dimension  of this element */

  iconnect_ptr =
      ei[pg->imtrx]
          ->iconnect_ptr; /* find pointer to beginning  of this element's connectivity list */

  /* subgrid or subelement integration setup */
  if (pd->gv[FILL] && ls != NULL && ls->Integration_Depth > 0 && ls->elem_overlap_state) {
    Subgrid_Int.active = TRUE;

#ifdef SUBELEMENT_FOR_SUBGRID
    /* DRN: you can also do subgrid integration with the following.
       it recursive divides element to create small subelements and then
       creates subgrid integration points */
    Subgrid_Int.ip_total = get_subelement_integration_pts(&Subgrid_Int.s, &Subgrid_Int.wt,
                                                          &Subgrid_Int.ip_sign, 0., -2, 0);
    /*DPRINTF(stderr,"DEBUG ielem=%d, ip_total=%d\n",ei[pg->imtrx]->ielem,Subgrid_Int.ip_total);*/
#else
    Subgrid_Int.ip_total = get_subgrid_integration_pts(
        Subgrid_Tree, &element_search_grid, &Subgrid_Int.s, &Subgrid_Int.wt, ls->Length_Scale);
    /*DPRINTF(stderr,"DEBUG ielem=%d, ip_total=%d\n",ei[pg->imtrx]->ielem,Subgrid_Int.ip_total);*/
#endif
#if 0
      print_subgrid_integration_pts ( Subgrid_Int.s, Subgrid_Int.wt, Subgrid_Int.ip_total );
#endif
  } else if (pd->gv[FILL] && ls != NULL && ls->SubElemIntegration && ls->elem_overlap_state) {
    Subgrid_Int.active = TRUE;
    Subgrid_Int.ip_total = get_subelement_integration_pts(&Subgrid_Int.s, &Subgrid_Int.wt,
                                                          &Subgrid_Int.ip_sign, 0., -2, 0);
    if (neg_elem_volume) {
#if 0
	  /* Squawk then turn off the switch. */
          GOMA_WH( -1, "Negative subelement volume detected.");
	  neg_elem_volume = FALSE;
#endif
#if 1
      /* Squawk then fail timestep. */
      fprintf(stderr, "Negative subelement volume in element (%d)!\n", ielem + 1);
      return -1;
#endif
    }
#if 0
      print_subgrid_integration_pts ( Subgrid_Int.s, Subgrid_Int.wt, Subgrid_Int.ip_total );
#endif
  } else if (pd->gv[FILL] && ls != NULL && ls->AdaptIntegration && ls->elem_overlap_state) {
    Subgrid_Int.active = TRUE;
  } else {
    Subgrid_Int.active = FALSE;
  }

  /*
   * Clean out local element contribution (lec) accumulator...
   */
  zero_lec();

  /* XFEM setup */
  /* DRN-now using check_xfem_contribution to do this */
#if 0
  if ( xfem!= NULL ) kill_extended_eqns( x, exo );
#endif
  if (ls != NULL)
    ls->Elem_Sign = 0;

  /******************************************************************************/
  /*                              BLOCK 1.5                                     */
  /*             INITIAL CALCULATIONS AT CENTROID OF THE ELEMENT                */
  /******************************************************************************/

  /*
   * set element size variable to zero if not used -- just to be safe
   */

  memset(pg_data.h, 0, sizeof(double) * DIM);
  memset(pg_data.hh, 0, sizeof(double) * DIM * DIM);
  memset(pg_data.dh_dxnode, 0, sizeof(double) * MDE * DIM);
  memset(pg_data.hsquared, 0, sizeof(double) * DIM);
  memset(pg_data.hhv, 0, sizeof(double) * DIM * DIM);
  memset(pg_data.dhv_dxnode, 0, sizeof(double) * MDE * DIM);
  memset(pg_data.v_avg, 0, sizeof(double) * DIM);
  memset(pg_data.dv_dnode, 0, sizeof(double) * MDE * DIM);
  pg_data.mu_avg = 0.;
  pg_data.rho_avg = 0.;

  /* get element level constants for upwinding and
     stabilized schemes, if necessary */
  /*
   * If PSPG is turned on, then calculate the centroid viscosity
   * for use in the PSPG formulas. Note, we actually call
   * load_basis_functions here. Is this big penalty necessary or
   * can be piggyback on top of one gauss point?
   */
  if (PSPG) {
    if (PSPG == 1) {
      pspg_local = 0;
    }
    /* This is the flag for the standard local PSPG */
    else if (PSPG == 2) {
      pspg_local = 1;
    }
  }

  if ((PSPG || (mp->Mwt_funcModel == SUPG)) && pde[R_PRESSURE] && pde[R_MOMENTUM1]) {
    xi[0] = 0.0;
    xi[1] = 0.0;
    xi[2] = 0.0;
    (void)load_basis_functions(xi, bfd);
    setup_shop_at_point(ielem, xi, exo);
    pg_data.mu_avg = element_viscosity();
    pg_data.rho_avg = density(NULL, time_value);

    //    double mu = 0;
    //    double rho = 0;
    //    for (ip = 0; ip < num_local_nodes; ip ++)
    //      {
    //        find_nodal_stu(ip, ielem_type, &xi[0], &xi[1], &xi[2]);
    //        setup_shop_at_point(ei[pg->imtrx]->ielem, xi, exo);
    //        double gamma[DIM][DIM];
    //        for ( int a=0; a<VIM; a++)
    //          {
    //            for ( int b=0; b<VIM; b++)
    //              {
    //                gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
    //              }
    //          }
    //        mu += viscosity(gn, gamma, NULL);
    //        rho += density(NULL, time_value);
    //      }
    //    pg_data.mu_avg = mu / ((double) num_local_nodes);
    //    pg_data.rho_avg = rho / ((double) num_local_nodes);

    if (pspg_local) {
      h_elem_siz(pg_data.hsquared, pg_data.hhv, pg_data.dhv_dxnode, pde[R_MESH1]);
      element_velocity(pg_data.v_avg, pg_data.dv_dnode, exo);
    }
  }

  if (Cont_GLS && pde[R_PRESSURE] && pde[R_MOMENTUM1]) {
    xi[0] = 0.0;
    xi[1] = 0.0;
    xi[2] = 0.0;
    (void)load_basis_functions(xi, bfd);
    pg_data.mu_avg = element_viscosity();
    pg_data.rho_avg = density(NULL, time_value);

    if (Cont_GLS == 2) {
      h_elem_siz(pg_data.hsquared, pg_data.hhv, pg_data.dhv_dxnode, pde[R_MESH1]);
      element_velocity(pg_data.v_avg, pg_data.dv_dnode, exo);
    }
  }

  if (pd->gv[FILL] || pd->gv[PHASE1]) /* UMR fix for non-FILL problems */
  {
    if (ls != NULL) {
      h_elem_siz(pg_data.hsquared, pg_data.hhv, pg_data.dhv_dxnode, pde[R_MESH1]);
    }
  }

  /*
   * The current SUPG model requires a value for the element's size and the
   * average velocity. Evaluate this here.
   */
  if ((mp->Momentwt_funcModel == SUPG &&
       (pde[R_MOMENT0] || pde[R_MOMENT1] || pde[R_MOMENT2] || pde[R_MOMENT3]) &&
       pd->gv[R_MOMENTUM1]) ||
      (vn->wt_funcModel == SUPG && pde[R_STRESS11] && pd->gv[R_MOMENTUM1]) ||
      (mp->Spwt_funcModel == SUPG && pde[R_MASS] && (pd->gv[R_MOMENTUM1] || pd->gv[R_MESH1])) ||
      (mp->Mwt_funcModel == SUPG && pde[R_MOMENTUM1]) ||
      (mp->Ewt_funcModel == SUPG && pde[R_ENERGY] && (pd->gv[R_MOMENTUM1] || pd->gv[R_MESH1])) ||
      (mp->Ewt_funcModel == SUPG && pde[R_SHELL_ENERGY] && (pde[R_LUBP])) ||
      ((pde[R_FILL]) && ((tran->Fill_Weight_Fcn == FILL_WEIGHT_SUPG) ||
                         (tran->Fill_Weight_Fcn == FILL_WEIGHT_EXPLICIT)))) {
    h_elem_siz(pg_data.hsquared, pg_data.hhv, pg_data.dhv_dxnode, pd->gv[R_MESH1]);
    element_velocity(pg_data.v_avg, pg_data.dv_dnode, exo);
  }

  if (cr->MassFluxModel == HYDRODYNAMIC) {
    /* For shock capturing diffusivity in Phillips model */
    h_elem_siz(pg_data.hsquared, pg_data.hhv, pg_data.dhv_dxnode, pde[R_MESH1]);
  }
  if (mp->DiffusivityModel[0] == HYDRO && mp->QTensorDiffType[0] == CONSTANT) {
    h_elem_siz(pg_data.hsquared, pg_data.hhv, pg_data.dhv_dxnode, pde[R_MESH1]);
    assemble_qtensor(pg_data.hsquared);
  }

  if (cr->MassFluxModel == HYDRODYNAMIC_QTENSOR_OLD) {
    h_elem_siz(pg_data.hsquared, pg_data.hhv, pg_data.dhv_dxnode, pde[R_MESH1]);
    assemble_qtensor(pg_data.hsquared);
  }

  /* Ryan's new qtensor model */
  else if (mp->DiffusivityModel[0] == HYDRO && mp->QtensorExtensionPModel == CONSTANT) {
    h_elem_siz(pg_data.hsquared, pg_data.hhv, pg_data.dhv_dxnode, pde[R_MESH1]);
    assemble_new_qtensor(pg_data.hsquared);
  }

  /*
   * Here we will do a quick element-level manipulation for the use of Lagrange
   * multipliers for fluid-structural surfaces (viz. overset grid). This section may change if
   * the assemble_volume_lagrange_multiplier guts change.
   * Quickly trivialize residuals and jacobian for elements that don't contain
   * an isosurface. For elements not containing a solid/fluid boundary, we will
   * set set the Jacobian and Residual to be such that no updates occur. Actually, this
   * test should be done at a much higher level, and not within a Gauss Integration loop, but
   * we will do this here for now.   PRS: If we go forward with this we need to put this
   * in a subroutine
   */

  ls_old = ls;
  if (!pde[R_MESH1]) {
    if (pfd != NULL)
      ls = pfd->ls[0]; /* this is a major hack */
    if (pde[R_LAGR_MULT1] && !ls->elem_overlap_state) {
      for (b = 0; b < pd->Num_Dim; b++) {
        for (i = 0; i < ei[pg->imtrx]->dof[R_LAGR_MULT1 + b]; i++) {
          if (af->Assemble_Residual) {
            lm_dof = ei[pg->imtrx]->gun_list[R_LAGR_MULT1 + b][i];
            eqn = upd->ep[pg->imtrx][R_LAGR_MULT1 + b];
            var = upd->vp[pg->imtrx][LAGR_MULT1 + b];
            lec->R[LEC_R_INDEX(eqn, i)] = x[lm_dof];
          }

          if (af->Assemble_Jacobian) {
            zero_lec_row(lec->J, eqn, i);
            lec->J[LEC_J_INDEX(eqn, var, i, i)] = 1.0;
          }
        }
      }
    }
  } else {
    if (pde[R_LAGR_MULT1]) {
      for (b = 0; b < pd->Num_Dim; b++) {
        for (i = 0; i < ei[pg->imtrx]->dof[R_LAGR_MULT1 + b]; i++) {
          if (af->Assemble_Residual) {
            eqn = upd->ep[pg->imtrx][R_LAGR_MULT1 + b];
            var = upd->vp[pg->imtrx][LAGR_MULT1 + b];
            lec->R[LEC_R_INDEX(eqn, i)] = 0.;
          }

          if (af->Assemble_Jacobian) {
            zero_lec_row(lec->J, eqn, i);
            lec->J[LEC_J_INDEX(eqn, var, i, i)] = 1.e-15;
          }
        }
      }
    }
  }
  ls = ls_old; /*Make things right again */

  /*
   * Here we will do another quick element-level manipulation for Eulerian Solid
   * mechanics problems using the level-set method and the assemble_real_solid
   * TALE routine.   In this situation we have a fixed mesh (hence no mesh motion
   * equations), and real-solid motion through the mesh. When not in the solid region
   * as delineated by the level-set field we will trivialize the real-solid equations
   * by setting the displacements to zero.  On elements that connect to elements through
   * which the level-set interface cuts, we will not assemble AT ALL, so as to achieve
   * the natural stress condition on the intefacial elements.
   */
  /* assemble fake equation for elements off interface */
  assemble_rs = TRUE;

  if (pde[R_SOLID1]) {
    if (pd->etm[pg->imtrx][R_SOLID1][(LOG2_MASS)])
      GOMA_EH(GOMA_ERROR, "Cannot do real inertia for TALE yet. Remove this line if trying to "
                          "perform EULERIAN solid mechanics");
    eqn = R_SOLID1;
    if (pd->TimeIntegration != STEADY && pd->etm[pg->imtrx][R_SOLID1][(LOG2_MASS)] &&
        !ls->elem_overlap_state && *esp->F[0] >= 0.0) {
      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        make_trivial = TRUE;
        if ((pd->i[pg->imtrx][eqn] == I_Q1) || (pd->i[pg->imtrx][eqn] == I_Q2)) {
          /* check all neighboring elements to see if any
             span the interface */
          I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
          for (j = exo->node_elem_pntr[I]; j < exo->node_elem_pntr[I + 1]; j++) {
            e = exo->node_elem_list[j];
            if (elem_overlaps_interface(e, x, exo, ls->Length_Scale)) {
              make_trivial = FALSE;
              assemble_rs = FALSE;
            }
          }
        }
      }
      if (make_trivial) {
        assemble_rs = FALSE;
        for (b = 0; b < pd->Num_Dim; b++) {
          for (i = 0; i < ei[pg->imtrx]->dof[R_SOLID1 + b]; i++) {
            if (af->Assemble_Residual) {
              eqn = upd->ep[pg->imtrx][R_SOLID1 + b];
              var = upd->vp[pg->imtrx][R_SOLID1 + b];
              lec->R[LEC_R_INDEX(eqn, i)] = 0.;
            }

            if (af->Assemble_Jacobian) {
              zero_lec_row(lec->J, eqn, i);
              lec->J[LEC_J_INDEX(eqn, var, i, i)] = 1.e-15;
            }
          }
        }
      }
    }
  }

  /******************************************************************************/
  /*                              BLOCK 1.6                                     */
  /*                Precalculations above the Gauss point loop                  */
  /******************************************************************************/

  if (pde[R_POR_LIQ_PRES]) {
    if (mp->Porous_Mass_Lump) {
      load_nodal_porous_properties(theta, delta_t);
    }
  }

  if ((pde[R_SHELL_SAT_OPEN]) || (pde[R_SHELL_SAT_OPEN_2])) {
    if (pde[R_SHELL_SAT_OPEN]) {
      load_nodal_shell_porous_properties(theta, delta_t, R_SHELL_SAT_OPEN);
    } else if (pde[R_SHELL_SAT_OPEN_2]) {
      load_nodal_shell_porous_properties(theta, delta_t, R_SHELL_SAT_OPEN_2);
    }
  }

  if (pde[R_PRESSURE]) {
    if (PSPP == 1) {
      err = assemble_projection_stabilization(exo, time_value);
      GOMA_EH(err, "assemble_projection_stabilization");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_projection_stabilization");
      if (err)
        return -1;
#endif
    }

    if (PSPP == 2) {
      err = assemble_PPPS_generalized(exo);
      GOMA_EH(err, "assemble_PPPS_generalized");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_PPPS_generalized");
      if (err)
        return -1;
#endif
    }

    if (PSPP == 3) {
      err = assemble_projection_time_stabilization(exo, time_value, theta, delta_t);
      GOMA_EH(err, "assemble_projection_stabilization");
#ifdef CHECK_FINITE
      CHECKFINITE("assemble_projection_stabilization");
#endif
    }
  }

  /* Initialize the flag for porous shell */
  if ((pde[R_SHELL_SAT_1]) || (pde[R_SHELL_SAT_2]) || (pde[R_SHELL_SAT_3])) {
    PorousShellOn = 1;
  }

  /* Do precalcations for porous shell hysteresis problem */
  if ((mp->PorousMediaType != CONTINUOUS) &&
      ((mp->PorousShellCapPresModel[0] == VAN_GENUCHTEN_HYST) ||
       (mp->PorousShellCapPresModel[0] == VAN_GENUCHTEN_HYST_EXT)) &&
      (PorousShellOn)) {
    //     if (af->Sat_hyst_reevaluate)
    //       {
    //        ebn = find_elemblock_index(ielem, exo);
    //        ielem_local = ielem - exo->eb_ptr[ebn];
    /* Determine what curve to follow and if switch is in order */
    //        err = evaluate_sat_hyst_criterion_nodal(ielem_local);
    //       }
  }

  /******************************************************************************/
  /*                              BLOCK 2.0a                                    */
  /*                   START OF VOLUME INTEGRATION LOOP                         */
  /*                LOOP OVER THE "REGULAR: QUADRATURE POINTS                   */
  /*          INITIALIZATION THAT IS DEPENDENT ON THE QUADRATURE POINT VALUE    */
  /******************************************************************************/

  ip_total = elem_info(NQUAD, ielem_type);

  /* Loop over all the Volume Quadrature integration points */

  if (pde[R_FILL]) /* No need to do this loop if there is no LS/FILL variable OK */
  {
    struct LS_Mass_Lumped_Penalty mass_lumped_penalty;
    if (ls != NULL && ls->Toure_Penalty) {
      assemble_ls_mass_lumped_penalty(&mass_lumped_penalty, ip_total, ielem_type, &pg_data);
    }

    for (ip = 0; ip < ip_total; ip++) {
      MMH_ip = ip;

      find_stu(ip, ielem_type, &s, &t, &u); /* find quadrature point */
      wt = Gq_weight(ip, ielem_type);       /* find quadrature weights for */
      if (ls != NULL)
        ls->Elem_Sign = 0;

      /* locations for current ip */

      /*
       * Local element coordinates...
       */

      xi[0] = s;
      xi[1] = t;
      xi[2] = u;

      fv->wt = wt;

      /*
       *	    PASS wt to assembly routines!!!!! ...
       */

      /*
       * Load up basis function information for each basis function
       * type needed by the current material. This is done in terms
       * of local element coordinates.
       */

      err = load_basis_functions(xi, bfd);
      GOMA_EH(err, "problem from load_basis_functions");

      /*
       * This has elemental Jacobian transformation and some
       * basic mesh derivatives...
       *
       * The physical space derivatives from beer_belly are still
       * "raw", in the sense that they do not yet include the
       * scale factors necessary to make them into *gradients*.
       * That is done in load_fv.
       */

      err = beer_belly();
      GOMA_EH(err, "beer_belly");
      if (neg_elem_volume)
        return -1;
      if (zero_detJ)
        return -1;

      /*
       * Load up field variable values at this Gauss point, but not
       * any gradients of field variables. Do load up the scale factors
       * and derivatives inherent to the coordinate system, such as
       * grad_(e_a) tensor... which is nontrivial in cylindrical
       * coordinates.
       */

      err = load_fv();
      GOMA_EH(err, "load_fv");

      /*
       * Here, load in the final part of the necessary basis function
       * information derivatives in the physical space coordinates.
       * Now we can form the true physical space gradient operators
       * because we have available the scale factors.
       */

      err = load_bf_grad();
      GOMA_EH(err, "load_bf_grad");

      err = load_fv_vector();
      GOMA_EH(err, "load_fv_vector");
      /*
       * Finally, load the mesh derivatives of the gradients of the
       * basis functions with respect to physical space coordinates.
       *
       * Only call this high computational intensity tensor workout if
       * we really need this information...
       */

      if (pd->gv[R_MESH1]) {
        err = load_bf_mesh_derivs();
        GOMA_EH(err, "load_bf_mesh_derivs");
      }

      /*
       * Load up physical space gradients of field variables at this
       * Gauss point.
       */
      err = load_fv_grads();
      GOMA_EH(err, "load_fv_grads");

      if (pd->gv[R_MESH1]) {
        err = load_fv_mesh_derivs(1);
        GOMA_EH(err, "load_fv_mesh_derivs");
      }

      /*
       * Presumably, we have all the pieces we need.
       */

      do_LSA_mods(LSA_VOLUME);

      if (ls != NULL && ls->Evolution == LS_EVOLVE_SLAVE) {
        err = assemble_fill_fake(theta, delta_t);
        GOMA_EH(err, "assemble_fill_fake");
#ifdef CHECK_FINITE
        err = CHECKFINITE("assemble_fill_fake");
        if (err)
          return -1;
#endif /* CHECK_FINITE */
      } else if (tran->Fill_Equation == FILL_EQN_EXT_V && pde[R_EXT_VELOCITY]) {
        ls_old = ls;
        if (pfd != NULL)
          ls = pfd->ls[0];
        err = assemble_fill_ext_v(theta, delta_t, pg_data.hsquared, pg_data.hh, pg_data.dh_dxnode);
        ls = ls_old; /*Make things right again */
        GOMA_EH(err, "assemble_fill_ext_v");
#ifdef CHECK_FINITE
        err = CHECKFINITE("assemble_fill_ext_v");
        if (err)
          return -1;
#endif /* CHECK_FINITE */
      } else if (tran->Fill_Equation == FILL_EQN_ADVECT) {
        err = assemble_fill(theta, delta_t, &pg_data, R_FILL, xi, exo, time_value,
                            &mass_lumped_penalty);
        GOMA_EH(err, "assemble_fill");
#ifdef CHECK_FINITE
        err = CHECKFINITE("assemble_fill");
        if (err)
          return -1;
#endif /* CHECK_FINITE */
        if (pfd != NULL && pde[R_EXT_VELOCITY]) {
          ls_old = ls;
          ls = pfd->ls[0];
          err = assemble_fill(theta, delta_t, &pg_data, R_PHASE1, xi, exo, time_value,
                              &mass_lumped_penalty);
          ls = ls_old; /*Make things right again */
          GOMA_EH(err, "assemble_fill");
#ifdef CHECK_FINITE
          err = CHECKFINITE("assemble_fill");
          if (err)
            return -1;
#endif /* CHECK_FINITE */
        }
      }
    }
    /******************************************************************************/
  }
  /* END  for (ip = 0; ip < ip_total; ip++)                               */

  /******************************************************************************/
  /*                              BLOCK 2.0b                                    */
  /*                   START OF VOLUME INTEGRATION LOOP                         */
  /*                FOR EQNS THAT MIGHT DEPEND ON LEVEL SET FUNCTIONS           */
  /******************************************************************************/

  if (ls == NULL || !ls->elem_overlap_state) {
    /* case 1: normal gauss integration */
    ip_total = elem_info(NQUAD, ielem_type);
  } else if (ls->Integration_Depth == 0 && !ls->SubElemIntegration && !ls->AdaptIntegration) {
    if (xfem != NULL) {
      /* case 2: XFEM with normal gauss integration */
      ip_total = 2 * elem_info(NQUAD, ielem_type);
    } else {
      /* case 3: no XFEM with normal gauss integration */
      ip_total = elem_info(NQUAD, ielem_type);
    }
  } else if (Subgrid_Int.active && ls->Integration_Depth > 0) {
    if (xfem != NULL) {
      /* case 4: XFEM with subgrid integration */
      ip_total = 2 * Subgrid_Int.ip_total;
    } else {
      /* case 5: NO XFEM with subgrid integration */
      ip_total = Subgrid_Int.ip_total;
    }
  } else if (Subgrid_Int.active && ls->SubElemIntegration) {
    /* case 6: subelement integration */
    ip_total = Subgrid_Int.ip_total;
  } else if (Subgrid_Int.active && ls->AdaptIntegration) {
    /* case 7: Adaptive weight integration - same as case 2 and 3 */
    ip_total = elem_info(NQUAD, ielem_type);

    if (xfem != NULL) {
      for (i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
        ls_F[i] = *esp->F[i];
      }
      i = adaptive_weight(ad_wtpos, ip_total, ielem_dim, ls_F, ls->Length_Scale, 2, ielem_type);
      GOMA_WH(i, "problem with adaptive weight routine");
      for (i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
        ls_F[i] = -ls_F[i];
      }
      i = adaptive_weight(ad_wtneg, ip_total, ielem_dim, ls_F, ls->Length_Scale, 2, ielem_type);
      GOMA_WH(i, "problem with adaptive weight routine");
      ip_total = 2 * elem_info(NQUAD, ielem_type);
    }
  } else {
    GOMA_EH(GOMA_ERROR, "Unrecognized integration scheme!");
  }

  /* Loop over all the Volume Quadrature integration points */

  for (ip = 0; ip < ip_total; ip++) {
    MMH_ip = ip;

    if (ls == NULL || !ls->elem_overlap_state) {
      /* case 1: normal gauss integration */
      find_stu(ip, ielem_type, &s, &t, &u); /* find quadrature point */
      wt = Gq_weight(ip, ielem_type);       /* find quadrature weights for */
    } else if (ls->Integration_Depth == 0 && !ls->SubElemIntegration && !ls->AdaptIntegration) {
      if (xfem != NULL) {
        /* case 2: XFEM with normal gauss integration */
        find_stu(ip / 2, ielem_type, &s, &t, &u); /* find quadrature point */
        wt = Gq_weight(ip / 2, ielem_type);       /* find quadrature weights for */
        if (ip % 2 == 0)
          ls->Elem_Sign = -1;
        else
          ls->Elem_Sign = 1;
      } else {
        /* case 3: no XFEM with normal gauss integration */
        find_stu(ip, ielem_type, &s, &t, &u); /* find quadrature point */
        wt = Gq_weight(ip, ielem_type);       /* find quadrature weights for */
        ls->Elem_Sign = 0;
      }
    } else if (Subgrid_Int.active && ls->Integration_Depth > 0) {
      if (xfem != NULL) {
        /* case 4: XFEM with subgrid integration */
        s = Subgrid_Int.s[ip / 2][0];
        t = Subgrid_Int.s[ip / 2][1];
        u = Subgrid_Int.s[ip / 2][2];
        wt = Subgrid_Int.wt[ip / 2];
        if (ip % 2 == 0)
          ls->Elem_Sign = -1;
        else
          ls->Elem_Sign = 1;
      } else {
        /* case 5: NO XFEM with subgrid integration */
        s = Subgrid_Int.s[ip][0];
        t = Subgrid_Int.s[ip][1];
        u = Subgrid_Int.s[ip][2];
        wt = Subgrid_Int.wt[ip];
        ls->Elem_Sign = 0;
      }
    } else if (Subgrid_Int.active && ls->SubElemIntegration) {
      /* case 6: subelement integration */
      s = Subgrid_Int.s[ip][0];
      t = Subgrid_Int.s[ip][1];
      u = Subgrid_Int.s[ip][2];
      wt = Subgrid_Int.wt[ip];
      ls->Elem_Sign = Subgrid_Int.ip_sign[ip];
    } else if (Subgrid_Int.active && ls->AdaptIntegration) {
      /* case 7: Adaptive weight integration - similiar to case 2 and 3 */
      if (xfem != NULL) {
        find_stu(ip / 2, ielem_type, &s, &t, &u); /* find quadrature point */
        if (ip % 2 == 0) {
          wt = ad_wtneg[ip_total / 2 - 1 - ip / 2];
          ls->Elem_Sign = -1;
        } else {
          wt = ad_wtpos[ip_total / 2 - 1 - ip / 2];
          ls->Elem_Sign = 1;
        }
      } else {
        find_stu(ip, ielem_type, &s, &t, &u); /* find quadrature point */
        wt = Gq_weight(ip, ielem_type);       /* find quadrature weights for */
        ls->Elem_Sign = 0;
      }
    } else {
      GOMA_EH(GOMA_ERROR, "Unrecognized integration scheme!");
    }

    /* locations for current ip */

    /*
     * Local element coordinates...
     */

    xi[0] = s;
    xi[1] = t;
    xi[2] = u;

    fv->wt = wt;

    /*
     *	    PASS wt to assembly routines!!!!! ...
     */

    /*
     * Load up basis function information for each basis function
     * type needed by the current material. This is done in terms
     * of local element coordinates.
     */

    err = load_basis_functions(xi, bfd);

    GOMA_EH(err, "problem from load_basis_functions");

    /*
     * This has elemental Jacobian transformation and some
     * basic mesh derivatives...
     *
     * The physical space derivatives from beer_belly are still
     * "raw", in the sense that they do not yet include the
     * scale factors necessary to make them into *gradients*.
     * That is done in load_fv.
     */

    err = beer_belly();
    GOMA_EH(err, "beer_belly");
    if (neg_elem_volume)
      return -1;
    if (zero_detJ)
      return -1;

    /*
     * Load up field variable values at this Gauss point, but not
     * any gradients of field variables. Do load up the scale factors
     * and derivatives inherent to the coordinate system, such as
     * grad_(e_a) tensor... which is nontrivial in cylindrical
     * coordinates.
     */
    err = load_fv();
    GOMA_EH(err, "load_fv");

    /*
     * Here, load in the final part of the necessary basis function
     * information derivatives in the physical space coordinates.
     * Now we can form the true physical space gradient operators
     * because we have available the scale factors.
     */

    err = load_bf_grad();
    GOMA_EH(err, "load_bf_grad");

    err = load_fv_vector();
    GOMA_EH(err, "load_fv_vector");

    /*
     * Finally, load the mesh derivatives of the gradients of the
     * basis functions with respect to physical space coordinates.
     *
     * Only call this high computational intensity tensor workout if
     * we really need this information...
     */

    if (pd->gv[R_MESH1]) {
      err = load_bf_mesh_derivs();
      GOMA_EH(err, "load_bf_mesh_derivs");
    }

    /*
     * Load up physical space gradients of field variables at this
     * Gauss point.
     */
    err = load_fv_grads();
    GOMA_EH(err, "load_fv_grads");

    if (pd->gv[R_MESH1]) {
      err = load_fv_mesh_derivs(1);
      GOMA_EH(err, "load_fv_mesh_derivs");
    }

    /* special section for XFEM with diffuse integration */
    /* this is cases 2 and 4 from above */
    if (ls != NULL && xfem != NULL && ls->elem_overlap_state && !ls->SubElemIntegration &&
        !ls->AdaptIntegration) {
      load_lsi(ls->Length_Scale);
      if (ls->Elem_Sign == -1) {
        wt *= 1. - lsi->H;
        fv->wt *= 1. - lsi->H;
        if (fv->wt == 0.)
          continue;
      } else {
        wt *= lsi->H;
        fv->wt *= lsi->H;
        if (fv->wt == 0.)
          continue;
      }
    }
    /* debugging integration */
#if DEBUG_LS_INTEGRATION
    if (ls != NULL) {
      load_lsi(ls->Length_Scale);
      if (ls->Elem_Sign == -1) {
        ls->Neg_Vol += fv->wt * fv->h3 * bf[pd->ShapeVar]->detJ;
      } else if (ls->Elem_Sign == 1) {
        ls->Pos_Vol += fv->wt * fv->h3 * bf[pd->ShapeVar]->detJ;
      } else {
        ls->Neg_Vol += fv->wt * (1. - lsi->H) * fv->h3 * bf[pd->ShapeVar]->detJ;
        ls->Pos_Vol += fv->wt * lsi->H * fv->h3 * bf[pd->ShapeVar]->detJ;
      }
    }
#endif
    if (xfem != NULL) {
      compute_xfem_contribution(ams->npu);
    }

    /* QUESTION
     * Since we are already at a gauss point and have calculated
     * all the field variables and their derivatives, shouldn't
     * we calculate all the material properties here, rather than
     * in the subroutines?
     * Then, each mp would be calculated only ONCE per gauss point
     * and would be available for all the assembly routines
     *
     * ANSWER
     * Yes.  Write the source!
     */
    computeCommonMaterialProps_gp(time_value);

    /*
     * Presumably, we have all the pieces we need.
     */
    do_LSA_mods(LSA_VOLUME);

    if (vn->evssModel == EVSS_F || vn->evssModel == EVSS_GRADV) {
      err = assemble_stress_fortin(theta, delta_t, &pg_data);
      err = segregate_stress_update(x_update);
      GOMA_EH(err, "assemble_stress_fortin");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_stress_fortin");
      if (err)
        return -1;
#endif
    } else if (vn->evssModel == EVSS_G) {
      err = assemble_stress(theta, delta_t, pg_data.hsquared, pg_data.hhv, pg_data.dhv_dxnode,
                            pg_data.v_avg, pg_data.dv_dnode);
      GOMA_EH(err, "assemble_stress");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_stress");
      if (err)
        return -1;
#endif
    } else if (vn->evssModel == EVSS_L) {
      err = assemble_stress_level_set(theta, delta_t, pg_data.hsquared, pg_data.hhv,
                                      pg_data.dhv_dxnode, pg_data.v_avg, pg_data.dv_dnode);
      GOMA_EH(err, "assemble_stress_level_set");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_stress_level_set");
      if (err)
        return -1;
#endif
    } else if (vn->evssModel == LOG_CONF || vn->evssModel == LOG_CONF_GRADV) {
      err = assemble_stress_log_conf(theta, delta_t, &pg_data);

      GOMA_EH(err, "assemble_stress_log_conf");
      if (err)
        return -1;
      err = segregate_stress_update(x_update);
      GOMA_EH(err, "assemble_stress_update");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_stress_log_conf");
      if (err)
        return -1;
#endif
    } else if (vn->evssModel == SQRT_CONF) {
      err = assemble_stress_sqrt_conf(theta, delta_t, &pg_data);

      GOMA_EH(err, "assemble_stress_sqrt_conf");
      if (err)
        return -1;
      err = segregate_stress_update(x_update);
      GOMA_EH(err, "assemble_stress_update");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_stress_sqrt_conf");
      if (err)
        return -1;
#endif
    } else if (vn->evssModel == CONF) {
      err = assemble_stress_conf(theta, delta_t, &pg_data);

      GOMA_EH(err, "assemble_stress_conf");
      if (err)
        return -1;
      err = segregate_stress_update(x_update);
      GOMA_EH(err, "assemble_stress_update");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_stress_conf");
      if (err)
        return -1;
#endif
    } else if (vn->evssModel == LOG_CONF_TRANSIENT || vn->evssModel == LOG_CONF_TRANSIENT_GRADV) {
      err = assemble_stress_log_conf_transient(theta, delta_t, &pg_data);

      GOMA_EH(err, "assemble_stress_log_conf");
      if (err)
        return -1;
      err = segregate_stress_update(x_update);
      GOMA_EH(err, "assemble_stress_update");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_stress_log_conf");
      if (err)
        return -1;
#endif
    }

    if (pde[R_SHEAR_RATE]) {
      err = assemble_invariant(theta, delta_t);

      GOMA_EH(err, "assemble_invariant");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_invariant");
      if (err)
        return -1;
#endif
    }

    if (pde[R_ENORM]) {
      err = assemble_Enorm();
      GOMA_EH(err, "assemble_Enorm");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_Enorm");
      if (err)
        return -1;
#endif
    }

    if (pde[R_GRADIENT11]) {
      if (gn->ConstitutiveEquation == BINGHAM_MIXED) {
        err = assemble_rate_of_strain(theta, delta_t);
      } else {
        err = assemble_gradient(theta, delta_t);
      }
      GOMA_EH(err, "assemble_gradient");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_gradient");
      if (err)
        return -1;
#endif
    }

    if (pde[R_MESH1] && !pde[R_SHELL_CURVATURE] && !pde[R_SHELL_TENSION]) {
      err = assemble_mesh(time_value, theta, delta_t, ielem, ip, ip_total);
      GOMA_EH(err, "assemble_mesh");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_mesh");
      if (err)
        return -1;
#endif
      if (neg_elem_volume)
        return -1;
    }

    /*
     *  Load up porous media properties if necessary.
     *  These properties must be load here since they
     *  are needed to assemble the mesh Jacobian.
     *  Note that as of 4/06 we need to call this after mesh equation
     *  assembly as the anisotropic permeabilities depend on fv->deform_grad
     */
    if ((mp->PorousMediaType != CONTINUOUS) &&
        !(PorousShellOn)) /* Exclude porous shell from this operation */
    {
      err = load_porous_properties();
      GOMA_EH(err, "load_porous_properties");
      if ((mp->PorousMediaType == POROUS_UNSATURATED ||
           mp->PorousMediaType == POROUS_SHELL_UNSATURATED ||
           mp->PorousMediaType == POROUS_TWO_PHASE) &&
          (mp->SaturationModel == TANH_HYST)) {
        if (af->Sat_hyst_reevaluate) {
          /* gotta evaluate the various terms of the equation first */
          err = get_porous_part_sat_terms(&pm_terms, theta, delta_t);
          GOMA_EH(err, "problem in getting the partially-saturated porous  terms");

          /* Determine what curve to follow and if switch is in order */
          err = evaluate_sat_hyst_criterion(ip, PRS_mat_ielem, &pm_terms, theta, delta_t);

          /* Now that you have re-evaluated the switching parameter,
           * load up the saturation again in this case.
           */
          err = load_porous_properties();
          GOMA_EH(err, "load_porous_properties");
        }
      }
    }

    if (pde[R_MASS]) {
      if (pd->MassFluxModel == FICKIAN_SHELL) {
        err = assemble_shell_species(time_value, theta, delta_t, xi, &pg_data, exo);
      } else {
        err = assemble_mass_transport(time_value, theta, delta_t, &pg_data);
        GOMA_EH(err, "assemble_mass_transport");
#ifdef CHECK_FINITE
        err = CHECKFINITE("assemble_mass_transport");
        if (err)
          return -1;
#endif
        if (neg_elem_volume)
          return -1;
      }
    }

    if (pde[R_POR_LIQ_PRES] || pde[R_POR_SATURATION]) {
      err = assemble_porous_transport(time_value, theta, delta_t);
      GOMA_EH(err, "assemble_porous");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_porous");
      if (err)
        return -1;
#endif
      if (neg_elem_volume)
        return -1;
    }

    if (pde[R_SOLID1] && assemble_rs) {
      err = assemble_real_solid(time_value, theta, delta_t);
      GOMA_EH(err, "assemble_mesh");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_mesh");
      if (err)
        return -1;
#endif
    }

    if (pde[R_ENERGY]) {
      err = assemble_energy(time_value, theta, delta_t, &pg_data);
      GOMA_EH(err, "assemble_energy");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_energy");
      if (err)
        return -1;
#endif
    }

    if (pde[R_POTENTIAL]) {
      err = assemble_potential(time_value, theta, delta_t);
      GOMA_EH(err, "assemble_potential");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_potential");
      if (err)
        return -1;
#endif
    }

    if (pde[R_ACOUS_PREAL]) {
      err = assemble_acoustic(time_value, theta, delta_t, &pg_data, R_ACOUS_PREAL, ACOUS_PREAL);
      GOMA_EH(err, "assemble_acoustic");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_acoustic");
      if (err)
        return -1;
#endif
    }

    if (pde[R_ACOUS_PIMAG]) {
      err = assemble_acoustic(time_value, theta, delta_t, &pg_data, R_ACOUS_PIMAG, ACOUS_PIMAG);
      GOMA_EH(err, "assemble_acoustic");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_acoustic");
      if (err)
        return -1;
#endif
    }

    if (pde[R_ACOUS_REYN_STRESS]) {
      err = assemble_acoustic_reynolds_stress(time_value, theta, delta_t, &pg_data);
      GOMA_EH(err, "assemble_acoustic_reynolds_stress");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_acoustic_reynolds_stress");
      if (err)
        return -1;
#endif
    }

    if (pde[R_LIGHT_INTP]) {
      err = assemble_poynting(time_value, theta, delta_t, &pg_data, R_LIGHT_INTP, LIGHT_INTP);
      GOMA_EH(err, "assemble_poynting");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_poynting");
      if (err)
        return -1;
#endif
    }

    if (pde[R_LIGHT_INTM]) {
      err = assemble_poynting(time_value, theta, delta_t, &pg_data, R_LIGHT_INTM, LIGHT_INTM);
      GOMA_EH(err, "assemble_poynting");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_poynting");
      if (err)
        return -1;
#endif
    }

    if (pde[R_LIGHT_INTD]) {
      err = assemble_poynting(time_value, theta, delta_t, &pg_data, R_LIGHT_INTD, LIGHT_INTD);
      GOMA_EH(err, "assemble_poynting");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_poynting");
      if (err)
        return -1;
#endif
    }

    if (pde[R_RESTIME]) {
      err = assemble_poynting(time_value, theta, delta_t, &pg_data, R_RESTIME, RESTIME);
      GOMA_EH(err, "assemble_poynting");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_poynting");
      if (err)
        return -1;
#endif
    }
    if (((pde[R_EM_E1_REAL] && !pde[R_EM_H1_REAL]) || (pde[R_EM_E2_REAL] && !pde[R_EM_H2_REAL]) ||
         (pde[R_EM_E3_REAL] && !pde[R_EM_H3_REAL])) &&
        bf[EM_E1_REAL]->interpolation != I_N1) {
      err = assemble_ewave_curlcurl(time_value, theta, delta_t, R_EM_E1_REAL, EM_E1_REAL);
      GOMA_EH(err, "assemble_ewave");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_ewave");
      if (err)
        return -1;
#endif
    } else if (pde[R_EM_E1_REAL] && bf[EM_E1_REAL]->interpolation == I_N1) {
      err = assemble_ewave_nedelec(time_value);
      GOMA_EH(err, "assemble_ewave_nedelec");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_emwave");
      if (err)
        return -1;
#endif
    } else if (pde[R_EM_E1_REAL]) {
      err = assemble_emwave(time_value, theta, delta_t, &pg_data, R_EM_E1_REAL, EM_E1_REAL,
                            EM_E1_IMAG);
      GOMA_EH(err, "assemble_emwave");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_emwave");
      if (err)
        return -1;
#endif
    }
#if 0
    if (pde[R_EM_CONT_REAL]) {
      err = assemble_em_continuity();
      GOMA_EH(err, "assemble_em_continuity");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_em_continuity");
      if (err)
        return -1;
#endif
      if (neg_elem_volume)
        return -1;
    }

    if (pde[R_EM_E2_REAL] && !pde[R_EM_H2_REAL]) {
      //        err = assemble_ewave_tensor_bf(time_value, theta, delta_t,
      //                                R_EM_E2_REAL, EM_E2_REAL);
      //        GOMA_EH( err, "assemble_ewave");
      // #ifdef CHECK_FINITE
      //        err = CHECKFINITE("assemble_ewave");
      //        if (err) return -1;
      // #endif
    } else if (pde[R_EM_E2_REAL]) {
      err = assemble_emwave(time_value, theta, delta_t, &pg_data, R_EM_E2_REAL, EM_E2_REAL,
                            EM_E2_IMAG);
      GOMA_EH(err, "assemble_emwave");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_emwave");
      if (err)
        return -1;
#endif
    }

    if (pde[R_EM_E3_REAL] && !pde[R_EM_H3_REAL]) {
      //        err = assemble_ewave_tensor_bf(time_value, theta, delta_t,
      //                                R_EM_E3_REAL, EM_E3_REAL);
      //        GOMA_EH( err, "assemble_ewave");
      // #ifdef CHECK_FINITE
      //        err = CHECKFINITE("assemble_ewave");
      //        if (err) return -1;
      // #endif
    } else if (pde[R_EM_E3_REAL]) {
      err = assemble_emwave(time_value, theta, delta_t, &pg_data, R_EM_E3_REAL, EM_E3_REAL,
                            EM_E3_IMAG);
      GOMA_EH(err, "assemble_emwave");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_emwave");
      if (err)
        return -1;
#endif
    }

    if (pde[R_EM_E1_IMAG] && !pde[R_EM_H1_IMAG]) {
      //        err = assemble_ewave_tensor_bf(time_value, theta, delta_t,
      //                                R_EM_E1_IMAG, EM_E1_IMAG);
      //        GOMA_EH( err, "assemble_ewave");
      // #ifdef CHECK_FINITE
      //        err = CHECKFINITE("assemble_ewave");
      //        if (err) return -1;
      // #endif
    } else if (pde[R_EM_E1_IMAG]) {
      err = assemble_emwave(time_value, theta, delta_t, &pg_data, R_EM_E1_IMAG, EM_E1_IMAG,
                            EM_E1_REAL);
      GOMA_EH(err, "assemble_emwave");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_emwave");
      if (err)
        return -1;
#endif
    }

    if (pde[R_EM_E2_IMAG] && !pde[R_EM_H2_IMAG]) {
      //        err = assemble_ewave_tensor_bf(time_value, theta, delta_t,
      //                                R_EM_E2_IMAG, EM_E2_IMAG);
      //        GOMA_EH( err, "assemble_ewave");
      // #ifdef CHECK_FINITE
      //        err = CHECKFINITE("assemble_ewave");
      //        if (err) return -1;
      // #endif
    } else if (pde[R_EM_E2_IMAG]) {
      err = assemble_emwave(time_value, theta, delta_t, &pg_data, R_EM_E2_IMAG, EM_E2_IMAG,
                            EM_E2_REAL);
      GOMA_EH(err, "assemble_emwave");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_emwave");
      if (err)
        return -1;
#endif
    }

    if (pde[R_EM_E3_IMAG] && !pde[R_EM_H3_IMAG]) {
      //        err = assemble_ewave_tensor_bf(time_value, theta, delta_t,
      //                                R_EM_E3_IMAG, EM_E3_IMAG);
      //        GOMA_EH( err, "assemble_ewave");
      // #ifdef CHECK_FINITE
      //        err = CHECKFINITE("assemble_ewave");
      //        if (err) return -1;
      // #endif
    } else if (pde[R_EM_E3_IMAG]) {
      err = assemble_emwave(time_value, theta, delta_t, &pg_data, R_EM_E3_IMAG, EM_E3_IMAG,
                            EM_E3_REAL);
      GOMA_EH(err, "assemble_emwave");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_emwave");
      if (err)
        return -1;
#endif
    }

    if (pde[R_EM_H1_REAL]) {
      err = assemble_emwave(time_value, theta, delta_t, &pg_data, R_EM_H1_REAL, EM_H1_REAL,
                            EM_H1_IMAG);
      GOMA_EH(err, "assemble_emwave");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_emwave");
      if (err)
        return -1;
#endif
    }

    if (pde[R_EM_H2_REAL]) {
      err = assemble_emwave(time_value, theta, delta_t, &pg_data, R_EM_H2_REAL, EM_H2_REAL,
                            EM_H2_IMAG);
      GOMA_EH(err, "assemble_emwave");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_emwave");
      if (err)
        return -1;
#endif
    }

    if (pde[R_EM_H3_REAL]) {
      err = assemble_emwave(time_value, theta, delta_t, &pg_data, R_EM_H3_REAL, EM_H3_REAL,
                            EM_H3_IMAG);
      GOMA_EH(err, "assemble_emwave");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_emwave");
      if (err)
        return -1;
#endif
    }

    if (pde[R_EM_H1_IMAG]) {
      err = assemble_emwave(time_value, theta, delta_t, &pg_data, R_EM_H1_IMAG, EM_H1_IMAG,
                            EM_H1_REAL);
      GOMA_EH(err, "assemble_emwave");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_emwave");
      if (err)
        return -1;
#endif
    }

    if (pde[R_EM_H2_IMAG]) {
      err = assemble_emwave(time_value, theta, delta_t, &pg_data, R_EM_H2_IMAG, EM_H2_IMAG,
                            EM_H2_REAL);
      GOMA_EH(err, "assemble_emwave");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_emwave");
      if (err)
        return -1;
#endif
    }

    if (pde[R_EM_H3_IMAG]) {
      err = assemble_emwave(time_value, theta, delta_t, &pg_data, R_EM_H3_IMAG, EM_H3_IMAG,
                            EM_H3_REAL);
      GOMA_EH(err, "assemble_emwave");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_emwave");
      if (err)
        return -1;
#endif
    }
#endif
    if (pde[R_POR_SINK_MASS]) {
      err = assemble_pore_sink_mass(time_value, theta, delta_t);
      GOMA_EH(err, "assemble_pore_sink_mass");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_pore_sink_mass");
      if (err)
        return -1;
#endif
    }

    if (pde[R_EFIELD1]) {
      err = assemble_electric_field();
      GOMA_EH(err, "assemble_electric_field");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_electric_field");
      if (err)
        return -1;
#endif
    }

    if (pde[R_SURF_CHARGE]) {
      err = assemble_surface_charge(time_value, theta, delta_t, wt, xi, exo, R_SURF_CHARGE);
      GOMA_EH(err, "assemble_surface_charge");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_surface_charge");
      if (err)
        return -1;
#endif
    }

    if (pde[R_SHELL_USER]) {
      err = assemble_surface_charge(time_value, theta, delta_t, wt, xi, exo, R_SHELL_USER);
      GOMA_EH(err, "assemble_surface_charge");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_surface_charge");
      if (err)
        return -1;
#endif
    }

    if (pde[R_SHELL_BDYVELO]) {
      err = assemble_surface_charge(time_value, theta, delta_t, wt, xi, exo, R_SHELL_BDYVELO);
      GOMA_EH(err, "assemble_surface_charge");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_surface_charge");
      if (err)
        return -1;
#endif
    }

    if (pde[R_SHELL_LUBP]) {
      GOMA_EH(GOMA_ERROR, "SHELL_LUBP routine not available yet");
    }

    if (pde[R_LUBP]) {
      err = assemble_lubrication(R_LUBP, time_value, theta, delta_t, xi, exo);
      GOMA_EH(err, "assemble_lubrication");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_lubrication");
      if (err)
        return -1;
#endif
      if (neg_lub_height)
        return -1;
    }

    if (pde[R_LUBP_2]) {
      err = assemble_lubrication(R_LUBP_2, time_value, theta, delta_t, xi, exo);
      GOMA_EH(err, "assemble_lubrication");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_lubrication");
      if (err)
        return -1;
#endif
      if (neg_lub_height)
        return -1;
    }

    if (pde[R_MAX_STRAIN]) {
      err = assemble_max_strain();
      GOMA_EH(err, "assemble_max_strain");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_max_strain");
      if (err)
        return -1;
#endif
    }

    if (pde[R_CUR_STRAIN]) {
      err = assemble_cur_strain();
      GOMA_EH(err, "assemble_cur_strain");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_cur_strain");
      if (err)
        return -1;
#endif
    }

    if (pde[R_SHELL_LUB_CURV]) {
      err = assemble_lubrication_curvature(time_value, theta, delta_t, &pg_data, xi, exo);
      GOMA_EH(err, "assemble_lubrication_curvature");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_lubrication_curvature");
      if (err)
        return -1;
#endif
    }

    if (pde[R_SHELL_LUB_CURV_2]) {
      ls_old = ls;
      ls = pfd->ls[0];
      err = assemble_lubrication_curvature_2(time_value, theta, delta_t, &pg_data, xi, exo);
      GOMA_EH(err, "assemble_lubrication_curvature_2");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_lubrication_curvature_2");
      if (err)
        return -1;
#endif
      ls = ls_old;
    }

    if (pde[R_SHELL_ENERGY]) {
      err = assemble_shell_energy(time_value, theta, delta_t, xi, &pg_data, exo);
      GOMA_EH(err, "assemble_shell_energy");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_shell_energy");
      if (err)
        return -1;
#endif
    }

    if (pde[R_SHELL_DELTAH]) {
      err = assemble_shell_deltah(time_value, theta, delta_t, xi, exo);
      GOMA_EH(err, "assemble_shell_deltah");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_shell_deltah");
      if (err)
        return -1;
#endif
    }

    /* Both SHELL_FILMP and SHELL_FILMH have to be activated to solve film profile equation */

    if (pde[R_SHELL_FILMP] && pde[R_SHELL_FILMH]) {
      if (ei[pg->imtrx]->ielem_dim == 1) {
        err = assemble_film_1D(time_value, theta, delta_t, xi, exo);
        GOMA_EH(err, "assemble_film_1D");
#ifdef CHECK_FINITE
        err = CHECKFINITE("assemble_film_1D");
        if (err)
          return -1;
#endif
      } else {
        err = assemble_film(time_value, theta, delta_t, xi, exo);
        GOMA_EH(err, "assemble_film");
#ifdef CHECK_FINITE
        err = CHECKFINITE("assemble_film");
        if (err)
          return -1;
#endif
      }
    } else if ((!(pde[R_SHELL_FILMP])) && pde[R_SHELL_FILMH]) {
      GOMA_EH(-1, "Both SHELL_FILMP and SHELL_FILMH must be activated !");
    } else if (pde[R_SHELL_FILMP] && (!(pde[R_SHELL_FILMH]))) {
      GOMA_EH(-1, "Both SHELL_FILMP and SHELL_FILMH must be activated !");
    }

    /* Film particles equation has to be activated along with either SHELL_FILMP and SHELL_FILMH or
     * LUBP */

    if (pde[R_SHELL_FILMP] && pde[R_SHELL_FILMH] && pde[R_SHELL_PARTC]) {
      err = assemble_film_particles(time_value, theta, delta_t, xi, &pg_data, exo);
      GOMA_EH(err, "assemble_film_particles");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_film_particles");
      if (err)
        return -1;
#endif
    }

    else if (pde[R_LUBP] && pde[R_SHELL_PARTC]) {
      err = assemble_film_particles(time_value, theta, delta_t, xi, &pg_data, exo);
      GOMA_EH(err, "assemble_film_particles");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_film_particles");
      if (err)
        return -1;
#endif
    }

    else if ((!(pde[R_SHELL_FILMP])) && pde[R_SHELL_FILMH] && pde[R_SHELL_PARTC]) {
      GOMA_EH(-1, " SHELL_PARTC requires SHELL_FILMP and SHELL_FILMH !");
    }

    else if (pde[R_SHELL_FILMP] && !(pde[R_SHELL_FILMH]) && pde[R_SHELL_PARTC]) {
      GOMA_EH(-1, " SHELL_PARTC requires SHELL_FILMP and SHELL_FILMH !");
    }

    if (pde[R_SHELL_SAT_CLOSED]) {
      err = assemble_porous_shell_closed(theta, delta_t, xi, exo);
      GOMA_EH(err, "assemble_porous_shell_closed");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_porous_shell_closed");
      if (err)
        return -1;
#endif
    }

    if (pde[R_SHELL_SAT_GASN]) {
      if (!pde[R_SHELL_SAT_CLOSED])
        GOMA_EH(-1, "SHELL_SAT_GASN required SHELL_SAT_CLOSED!");
      err = assemble_porous_shell_gasn(theta, delta_t, xi, exo);
      GOMA_EH(err, "assemble_porous_shell_gasn");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_porous_shell_gasn");
      if (err)
        return -1;
#endif
    }

    if (pde[R_SHELL_SAT_OPEN]) {
      err = assemble_porous_shell_open(theta, delta_t, xi, exo);
      GOMA_EH(err, "assemble_porous_shell_open");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_porous_shell_open");
      if (err)
        return -1;
#endif
    }

    if (pde[R_SHELL_SAT_OPEN_2]) {
      err = assemble_porous_shell_open_2(theta, delta_t, xi, exo);
      GOMA_EH(err, "assemble_porous_shell_open_2");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_porous_shell_open_2");
      if (err)
        return -1;
#endif
    }

    if ((pde[R_SHELL_SAT_1]) || (pde[R_SHELL_SAT_2]) || (pde[R_SHELL_SAT_3])) {
      err = assemble_porous_shell_saturation(theta, delta_t, xi, exo);
      GOMA_EH(err, "assemble_porous_shell_saturation");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_porous_shell_saturation");
      if (err)
        return -1;
#endif
    }

    if (pde[R_SHELL_ANGLE1]) {
      err = assemble_shell_angle(time_value, theta, delta_t, xi, exo);
      GOMA_EH(err, "assemble_shell_angle");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_shell_angle");
      if (err)
        return -1;
#endif
    }

    if (pde[R_N_DOT_CURL_V]) {
      err = assemble_shell_surface_rheo_pieces(time_value, theta, delta_t, xi, exo);
      GOMA_EH(err, "assemble_shell_surface_rheo_pieces");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_shell_surface_rheo_pieces");
      if (err)
        return -1;
#endif
    }

    /* Shell structure with both sh_K and sh_tens: */
    if (pde[R_SHELL_CURVATURE] && pde[R_SHELL_TENSION]) {
      err = assemble_shell_structure(time_value, theta, delta_t, wt, xi, exo);
      GOMA_EH(err, "assemble_shell_structure");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_shell_structure");
      if (err)
        return -1;
#endif
      if (pde[R_MESH1]) {
        err = assemble_shell_coordinates(time_value, theta, delta_t, wt, xi, exo);
        GOMA_EH(err, "assemble_shell_coordinates");
#ifdef CHECK_FINITE
        err = CHECKFINITE("assemble_shell_coordinates");
        if (err)
          return -1;
#endif
      }
    }

    /* Shell structure with only sh_tens, not sh_K */
    else if (!pde[R_SHELL_CURVATURE] && pde[R_SHELL_TENSION]) {
      err = assemble_shell_tension(time_value, theta, delta_t, wt, xi, exo);
      GOMA_EH(err, "assemble_shell_tension");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_shell_tension");
      if (err)
        return -1;
#endif
      if (pde[R_MESH1]) {
        err = assemble_shell_coordinates(time_value, theta, delta_t, wt, xi, exo);
        GOMA_EH(err, "assemble_shell_coordinates");
#ifdef CHECK_FINITE
        err = CHECKFINITE("assemble_shell_coordinates");
        if (err)
          return -1;
#endif
      }
    }

    /* Shell structure with only sh_K, not sh_tens is verboten! */
    else if (pde[R_MESH1] && pde[R_SHELL_CURVATURE] && !pde[R_SHELL_TENSION] &&
             !pde[R_SHELL_CURVATURE2]) {
      GOMA_EH(-1, "Must have both SHELL_TENSION AND SHELL_CURVATURE eqs");
    }

    if (pde[R_SHELL_DIFF_FLUX]) {
      err = assemble_shell_diffusion(time_value, theta, delta_t, wt, xi, exo);
      GOMA_EH(err, "assemble_shell_diffusion");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_shell_diffusion");
      if (err)
        return -1;
#endif
    }

    /* Web structure with both sh_K and sh_tens: */
    if (pde[R_SHELL_CURVATURE] && pde[R_SHELL_TENSION] &&
        (mp->FSIModel == FSI_SHELL_ONLY || mp->FSIModel == FSI_SHELL_ONLY_MESH)) {
      err = assemble_shell_web_structure(time_value, theta, delta_t, wt, xi, exo);
      GOMA_EH(err, "assemble_shell_web_structure");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_shell_web_structure");
      if (err)
        return -1;
#endif
      if (pde[R_MESH1]) {
        err = assemble_shell_web_coordinates(time_value, theta, delta_t, wt, xi, exo);
        GOMA_EH(err, "assemble_shell_web_coordinates");
#ifdef CHECK_FINITE
        err = CHECKFINITE("assemble_shell_web_coordinates");
        if (err)
          return -1;
#endif
      }
    }

    if ((pde[R_SHELL_DIFF_CURVATURE] || pde[R_SHELL_NORMAL1]) && !(pde[R_SHELL_NORMAL3]) &&
        !(pde[R_SHELL_CURVATURE])) {
      if (!pde[R_SHELL_NORMAL1] || !pde[R_SHELL_NORMAL2]) {
        GOMA_EH(GOMA_ERROR,
                "Both SHELL_NORMAL1 and SHELL_NORMAL2 required with SHELL_DIFF_CURVATURE eqn!");
      }
      err = assemble_shell_geometry(time_value, theta, delta_t, wt, xi, exo);
      GOMA_EH(err, "assemble_shell_geometry");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_shell_geometry");
      if (err)
        return -1;
#endif
    }

    if ((pde[R_SHELL_NORMAL1] && pde[R_SHELL_NORMAL2] && pde[R_SHELL_NORMAL3]) ||
        (pde[R_MESH1] && pde[R_SHELL_NORMAL1] && pde[R_SHELL_NORMAL2])) {

      err = assemble_shell_normal(xi, exo);
      GOMA_EH(err, "assemble_shell_normal");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_shell_normal");
      if (err)
        return -1;
#endif
    }

    if (pde[R_SHELL_CURVATURE] && pde[R_SHELL_CURVATURE2]) {
      err = assemble_shell_curvature(xi, exo);
      GOMA_EH(err, "assemble_shell_curvature");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_shell_curvature");
      if (err)
        return -1;
#endif
    }

    if ((pde[R_MESH1] && pde[R_SHELL_NORMAL1] && pde[R_SHELL_NORMAL2] && pde[R_SHELL_NORMAL3])) {
      err = assemble_shell_mesh(time_value, theta, delta_t, xi, exo);
      GOMA_EH(err, "assemble_shell_mesh");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_shell_mesh");
      if (err)
        return -1;
#endif
    }

    if (pde[R_TFMP_MASS] && pde[R_TFMP_BOUND]) {
      err = assemble_shell_tfmp(time_value, theta, delta_t, xi, &pg_data, exo);
      GOMA_EH(err, "assemble_shell_tfmp");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_shell_tfmp");
      if (err) {
        return -1;
      }
#endif
    }
    if (!pde[R_TFMP_MASS] && pde[R_TFMP_BOUND]) {
      err = assemble_shell_lubrication(time_value, theta, delta_t, xi, exo);
      GOMA_EH(err, "assemble_shell_lubrication");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_shell_lubrication");
      if (err) {
        return -1;
      }
#endif

      if (neg_lub_height) {
        GOMA_WH(
            -1,
            "returning from matrix fill because neg_lub_height after assemble_shell_lubrication");
        return -1;
      }
    }

    if (pde[R_MOMENTUM1]) {
      if (upd->SegregatedSolve) {
        err = assemble_momentum_segregated(time_value, theta, delta_t, &pg_data);
        GOMA_EH(err, "assemble_momentum");
#ifdef CHECK_FINITE
        CHECKFINITE("assemble_momentum");
#endif
      } else {
        err = assemble_momentum(time_value, theta, delta_t, h_elem_avg, &pg_data, xi, exo);
        GOMA_EH(err, "assemble_momentum");
#ifdef CHECK_FINITE
        CHECKFINITE("assemble_momentum");
#endif
      }
    }

    if (pde[R_PMOMENTUM1]) {
      err = assemble_pmomentum(time_value, theta, delta_t);
      GOMA_EH(err, "assemble_pmomentum");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_pmomentum");
      if (err)
        return -1;
#endif
    }

    if (pde[R_EDDY_NU]) {
      err = assemble_spalart_allmaras(time_value, theta, delta_t, &pg_data);
      GOMA_EH(err, "assemble_spalart_allmaras");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_spalart_allmaras");
      if (err)
        return -1;
#endif
    }

    if (pde[R_MOMENT0] || pde[R_MOMENT1] || pde[R_MOMENT2] || pde[R_MOMENT3]) {
      err = assemble_moments(time_value, theta, delta_t, &pg_data);
      GOMA_EH(err, "assemble_moments");
#ifdef CHECK_FINITE
      CHECKFINITE("assemble_moments");
#endif
    }
    if (pde[R_DENSITY_EQN]) {
      err = assemble_density();
      GOMA_EH(err, "assemble_density");
#ifdef CHECK_FINITE
      CHECKFINITE("assemble_density");
#endif
    }

    if (pde[R_FILL]) {
      if (tran->Fill_Equation == FILL_EQN_EIKONAL) {
        err = assemble_fill_gradf(theta, delta_t, pg_data.hsquared, pg_data.hh, pg_data.dh_dxnode);
        GOMA_EH(err, "assemble_fill_gradf");
#ifdef CHECK_FINITE
        err = CHECKFINITE("assemble_fill_gradf");
        if (err)
          return -1;
#endif /* CHECK_FINITE */
      }
    }

    if (pde[R_EXT_VELOCITY]) {
      ls_old = ls;
      if (pfd != NULL)
        ls = pfd->ls[0];
      assemble_extension_velocity(pg_data.hsquared, pg_data.hh, pg_data.dh_dxnode);
      ls = ls_old; /*Make things right again */
    }

    if (pde[R_CURVATURE]) {
      if (pde[R_NORMAL1])
        err = assemble_div_normals();
      else
        err = assemble_curvature();

      GOMA_EH(err, "assemble curvature projection");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble curvature projection");
      if (err)
        return -1;
#endif
    }

    if (pde[R_NORMAL1]) {
      err = assemble_normals();
      GOMA_EH(err, "assemble_normals");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_normals");
      if (err)
        return -1;
#endif
    }

    if (pde[R_PRESSURE]) {
      if (upd->SegregatedSolve) {
        err = assemble_continuity_segregated(time_value, theta, delta_t, &pg_data);
        GOMA_EH(err, "assemble_continuity");
#ifdef CHECK_FINITE
        CHECKFINITE("assemble_continuity");
#endif
        if (neg_elem_volume)
          return -1;
      } else {
        err = assemble_continuity(time_value, theta, delta_t, &pg_data);
        GOMA_EH(err, "assemble_continuity");
#ifdef CHECK_FINITE
        CHECKFINITE("assemble_continuity");
#endif
        if (neg_elem_volume)
          return -1;
      }
    }

    if (pde[R_PSTAR]) {
      if (upd->SegregatedSolve) {
        err = assemble_pstar(time_value, theta, delta_t, &pg_data);
        GOMA_EH(err, "assemble_pstar");
#ifdef CHECK_FINITE
        CHECKFINITE("assemble_pstar");
#endif
        if (neg_elem_volume)
          return -1;
      } else {
        GOMA_EH(-1, "Only implemented for Segregated Solve");
      }
    }

    if (pde[R_USTAR]) {
      if (upd->SegregatedSolve) {
        err = assemble_ustar(time_value, theta, delta_t, &pg_data);
        GOMA_EH(err, "assemble_ustar");
#ifdef CHECK_FINITE
        CHECKFINITE("assemble_ustar");
#endif
        if (neg_elem_volume)
          return -1;
      } else {
        GOMA_EH(-1, "Only implemented for Segregated Solve");
      }
    }

    if (pde[R_VORT_DIR1]) /* Then R_VORT_DIR2 and R_VORT_DIR3 should be on*/
    {
      err = assemble_vorticity_direction();
      GOMA_EH(err, "assemble_vorticity_direction");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_vorticity_direction");
      if (err)
        return -1;
#endif
    }

    if (pde[R_BOND_EVOLUTION]) {
      err = assemble_bond_evolution(time_value, theta, delta_t);
      GOMA_EH(err, "assemble_bond_evolution");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_bond_evolution");
      if (err)
        return -1;
#endif
    }

#ifdef PHASE_COUPLED_FILL
    if (pde[R_PHASE1] && !pde[R_EXT_VELOCITY]) {
      ls_old = ls;
      if (pfd != NULL)
        ls = pfd->ls[0];

      err = assemble_phase_function(time_value, theta, delta_t, xi, exo);
      GOMA_EH(err, "assemble_phase_functions");
      ls = ls_old;

#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_phase_functions");
      if (err)
        return -1;
#endif /* CHECK_FINITE */
      if (pfd->Use_Constraint == TRUE) {
        err = assemble_pf_constraint(delta_t, &(pfd->Constraint_Integral), augc[0].lm_value,
                                     pfd->jac_info->d_pf_lm, pfd->jac_info->d_lm_pf);
        GOMA_EH(err, " assemble_pf_constraint \n");
      }
    }
#else /* PHASE_COUPLED_FILL */
    if (pde[R_PHASE1]) {
      ls_old = ls;
      ls = pfd->ls[0];
      err = assemble_fill_fake(theta, delta_t);
      GOMA_EH(err, "assemble_fill_fake");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_fill_fake");
      if (err)
        return -1;
#endif /* CHECK_FINITE*/
      ls = ls_old;
    }
#endif /*PHASE_COUPLED_FILL */

    if (pd->VolumeIntegral > -1) {
      if (Num_Proc > 1 && dpi->elem_owner[ielem] != ProcID)
        owner = FALSE;
      err = assemble_volume(owner);
      GOMA_EH(err, "assemble_volume");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_volume");
      if (err)
        return -1;
#endif
    }

    if (pd->LSVelocityIntegral > -1) {
      if (Num_Proc > 1 && dpi->elem_owner[ielem] != ProcID)
        owner = FALSE;
      err = assemble_LSvelocity(owner, ielem);
      GOMA_EH(err, "assemble_LSvelocity");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_LSvelocity");
      if (err)
        return -1;
#endif
    }

    if (pde[LAGR_MULT1]) {
      /* Hmm, think you want to move this to the surface site on apply imbedded
       * BC first. Maybe keep this here as a blank template for volume integrals
       */
      /*err = assemble_volume_lagrange_multiplier( time_value, theta, delta_t ); */
      GOMA_EH(err, "assemble_volume_lagrange_multiplier");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_volume_lagrange_multiplier");
      if (err)
        return -1;
#endif
    }

    /******************************************************************************/
  }
  /* END  for (ip = 0; ip < ip_total; ip++)                               */

  if (pd->gv[R_LEVEL_SET] && ls != NULL)
    apply_embedded_colloc_bc(ielem, x, delta_t, theta, time_value, exo, dpi);

  if (pde[R_EXT_VELOCITY]) {
    ls_old = ls;
    if (pfd != NULL)
      ls = pfd->ls[0]; /* this is a major hack */
    assemble_boundary_extension_velocity(x, exo, dpi);
    ls = ls_old; /*Make things right again */
  }

  /*First apply BCS to embedded primary level-set surface. */
  if (ls != NULL) {
    if (pd->gv[ls->var] && ls->elem_overlap_state) {
      if (ls->Length_Scale == 0. || Do_Overlap) {
        apply_embedded_bc(ielem, x, delta_t, theta, time_value, &pg_data, -1, NULL, NULL, NULL,
                          exo);
#ifdef CHECK_FINITE
        err = CHECKFINITE("apply_embedded_bc");
        if (err)
          return -1;
#endif
      } else {
        apply_distributed_sources(ielem, ls->Length_Scale, x, exo, delta_t, theta, time_value,
                                  &pg_data, -1, NULL, NULL, NULL);
#ifdef CHECK_FINITE
        err = CHECKFINITE("apply_distributed_sources");
        if (err)
          return -1;
#endif
      }
    }
  }

  /*Now apply embedded bcs to pfd level-set surfaces like overset grids */
  if (pfd != NULL) {
    ls_old = ls;
    ls = pfd->ls[0];
    if (pfd->ls[0]->Evolution == LS_EVOLVE_SLAVE)

    {
      if (pd->gv[ls->var] && ls->elem_overlap_state) {
        if (ls->Length_Scale == 0. || Do_Overlap) {
          apply_embedded_bc(ielem, x, delta_t, theta, time_value, &pg_data, -1, NULL, NULL, NULL,
                            exo);
#ifdef CHECK_FINITE
          err = CHECKFINITE("apply_embedded_bc");
#endif
          if (err)
            return -1;
        } else {
          apply_distributed_sources(ielem, ls->Length_Scale, x, exo, delta_t, theta, time_value,
                                    &pg_data, -1, NULL, NULL, NULL);
#ifdef CHECK_FINITE
          err = CHECKFINITE("apply_distributed_sources");
          if (err)
            return -1;
#endif
        }
      }
    } else if (pd->gv[R_PHASE1] && !pd->gv[R_FILL]) {
      apply_distributed_sources(ielem, ls->Length_Scale, x, exo, delta_t, theta, time_value,
                                &pg_data, -1, NULL, NULL, NULL);
    } else if (pd->gv[R_PHASE1] && pd->gv[R_EXT_VELOCITY]) {
      if (pd->gv[ls->var] && ls->elem_overlap_state) {
        apply_embedded_bc(ielem, x, delta_t, theta, time_value, &pg_data, -1, NULL, NULL, NULL,
                          exo);
#ifdef CHECK_FINITE
        err = CHECKFINITE("apply_embedded_bc");
        if (err)
          return -1;
#endif
      }
    }
    ls = ls_old;
  }

  /**************************************************************************/
  /*                          BLOCK 2'                                      */
  /*                   START OF SURFACE INTEGRATION LOOP                    */
  /*         LOOP OVER THE NUMBER OF SURFACE QUADRATURE POINTS              */
  /*             For Discontinuous Galerkin implementations                 */
  /**************************************************************************/

  /* Loop over all the surface Quadrature integration points */

  if (discontinuous_mass) {
    int neighbor;
    int index, face;

    for (face = 0; face < ei[pg->imtrx]->num_sides; face++) {

      index = exo->elem_elem_pntr[ielem] + face;
      neighbor = exo->elem_elem_list[index];

      if (Num_Proc > 1) {
        GOMA_WH(-1, "DG doesn't seem to work in parallel");
      }

      id_side = face + 1;

      err = assemble_surface_species(exo, x, delta_t, theta, ielem_type, ielem_type_mass, id_side,
                                     neighbor, ielem, num_local_nodes);
      GOMA_EH(err, "assemble_surface_species");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_surface_species");
      if (err)
        return -1;
#endif
    }
  }

  if (discontinuous_stress) {
    int neighbor;
    int index, face;

    var = 4 * lec->max_dof * (MAX_PROB_VAR + MAX_CONC) * lec->max_dof;
    memset(lec->J_stress_neighbor, 0, sizeof(double) * var);

    if (Num_Proc > 1) {
      GOMA_EH(GOMA_ERROR, "Discontinuous Galerkin for stress not implemented in parallel");
      return -1;
    }

    for (face = 0; face < ei[pg->imtrx]->num_sides; face++) {

      index = exo->elem_elem_pntr[ielem] + face;
      neighbor = exo->elem_elem_list[index];

      id_side = face + 1;

      /*      if(neighbor != -1) */
      {
        err = assemble_surface_stress(exo, x, ams, x_update, delta_t, theta, ielem_type,
                                      ielem_type_mass, id_side, neighbor, ielem, num_local_nodes);
        GOMA_EH(err, "assemble_surface_stress");
#ifdef CHECK_FINITE
        err = CHECKFINITE("assemble_surface_stress");
        if (err)
          return -1;
#endif
      }
    }
  }

  /**************************************************************************
   *                          BLOCK 3 - Weak SURFACE Boundary Conditions
   *                   START OF SURFACE INTEGRATION LOOP
   * apply BC's that are added into the BOUNDARY term of the integrated-by-parts
   * weighted residual equations.  Do this before rotating mesh or momentum equations
   **************************************************************************/
  /******************************************************************************/
  if (first_elem_side_BC_array[ielem] != NULL) {
    elem_side_bc = first_elem_side_BC_array[ielem];
    /***************************************************************************
     *  begining of do while construct which loops over the sides of this
     *  element that have boundary conditions applied on to them.
     ***************************************************************************/
    do {

      call_int = 0;
      call_shell_grad = 0;
      call_sharp_int = FALSE;

      int call_nedelec = 0;

      for (ibc = 0; (bc_input_id = (int)elem_side_bc->BC_input_id[ibc]) != -1; ibc++) {
        if (BC_Types[bc_input_id].desc->method == WEAK_INT_SURF) {
          call_int = 1;
        }
        if (BC_Types[bc_input_id].desc->method == WEAK_INT_NEDELEC) {
          call_nedelec = 1;
        }
        if ((BC_Types[bc_input_id].desc->method == WEAK_SHELL_GRAD ||
             BC_Types[bc_input_id].desc->method == STRONG_SHELL_GRAD) &&
            (num_elem_friends[ielem] > 0) && (ielem_dim > 1)) {
          call_shell_grad = TRUE;
        }
        if (BC_Types[bc_input_id].desc->method == WEAK_SHARP_INT)
          call_sharp_int = TRUE;
      }

      if (call_int) {
        err = apply_integrated_bc(x, resid_vector, delta_t, theta, &pg_data, ielem, ielem_type,
                                  num_local_nodes, ielem_dim, iconnect_ptr, elem_side_bc,
                                  num_total_nodes, WEAK_INT_SURF, time_value, element_search_grid,
                                  exo);
        GOMA_EH(err, " apply_integrated_bc");
#ifdef CHECK_FINITE
        err = CHECKFINITE("apply_integrated_bc");
        if (err)
          return -1;
#endif
        if (neg_elem_volume)
          return -1;
      }

      if (call_shell_grad) {
        err = apply_shell_grad_bc(x, resid_vector, delta_t, theta, pg_data.h_elem_avg, pg_data.h,
                                  pg_data.mu_avg, pg_data.U_norm, ielem, ielem_type,
                                  num_local_nodes, ielem_dim, iconnect_ptr, elem_side_bc,
                                  num_total_nodes, WEAK_SHELL_GRAD, time_value, exo);
        GOMA_EH(err, " apply_shell_grad_bc");
#ifdef CHECK_FINITE
        err = CHECKFINITE("apply_shell_grad_bc");
        if (err)
          return -1;
#endif
        if (neg_elem_volume)
          return -1;
      }
      if (call_sharp_int) {
        err = apply_sharp_integrated_bc(x, resid_vector, time_value, delta_t, theta,
                                        pg_data.hsquared, ielem, ielem_type, num_local_nodes,
                                        ielem_dim, iconnect_ptr, elem_side_bc, WEAK_SHARP_INT, exo);
      }
      GOMA_EH(err, " apply_sharp_integrated_bc");
#ifdef CHECK_FINITE
      err = CHECKFINITE("apply_sharp_integrated_bc");
      if (err)
        return -1;
#endif
      if (neg_elem_volume)
        return -1;
      if (call_nedelec) {
        err = apply_nedelec_bc(x, resid_vector, delta_t, theta, &pg_data, ielem, ielem_type,
                               num_local_nodes, ielem_dim, iconnect_ptr, elem_side_bc,
                               num_total_nodes, WEAK_INT_NEDELEC, time_value, element_search_grid,
                               exo);
        GOMA_EH(err, " apply_integrated_bc");
#ifdef CHECK_FINITE
        err = CHECKFINITE("apply_integrated_bc");
        if (err)
          return -1;
#endif
        if (neg_elem_volume)
          return -1;
      }

    } while ((elem_side_bc = elem_side_bc->next_side_bc) != NULL);
  } /* END if (First_Elem_Side_BC_Array[ielem] != NULL) */

  /**************************************************************************
   *                          BLOCK 4 - Weak EDGE Boundary Conditions
   *                   START OF CURVE INTEGRATION LOOP
   * apply BC's that are added into the integrated-by-parts part of the BOUNDARY
   * term of the integrated-by-parts weighted residual equations.   This only applies
   * for conditions like capillary that require end-of-the-interface fluxes to
   * be added back in (i.e. if transport equations is integrated by parts twice)
   * Do this before rotating mesh or momentum equations
   **************************************************************************/
  if (First_Elem_Edge_BC_Array[pg->imtrx][ielem] != NULL) {
    elem_edge_bc = First_Elem_Edge_BC_Array[pg->imtrx][ielem];

    /******************************************************************************/
    do { /* begining of do while construct
          * which loops over the sides of this element that have boundary
          * conditions */
      /****************************************************************************/

      /*
       *  Set flags for subroutines to call for each boundary condition
       *  on this side
       */
      call_int = 0;
      for (ibc = 0; (bc_input_id = (int)elem_edge_bc->BC_input_id[ibc]) != -1; ibc++) {
        bct = BC_Types[bc_input_id].desc->method;
        if (bct == WEAK_INT_EDGE)
          call_int = 1;
      }
      /*
       * BOUNDARY CONDITIONS OF TYPE 1 - Integrated along edge -
       *        only do the weak integrated conditions here
       */
      if (call_int) {
        err = apply_integrated_curve_bc(x, resid_vector, delta_t, theta, ielem, ielem_type,
                                        num_local_nodes, ielem_dim, iconnect_ptr, elem_edge_bc,
                                        num_total_nodes, WEAK_INT_EDGE, exo);
        GOMA_EH(err, " apply_integrated_curve_bc");
#ifdef CHECK_FINITE
        err = CHECKFINITE("apply_integrated_curve_bc");
        if (err)
          return -1;
#endif
        if (neg_elem_volume)
          return -1;
      }
      /****************************************************************************/
    } while ((elem_edge_bc = elem_edge_bc->next_edge_bc) != NULL);
    /* END of do  while () construct				      */
    /******************************************************************************/
  } /* END if (First_Elem_edge_BC_Array[ielem] != NULL) 		      */

  /**************************************************************************
   *                          BLOCK 5 - Weak POINT Boundary Conditions
   *                   START OF POINT INTEGRATION LOOP  ?!?
   *  are we ever going to need a condition like this?  Probably not, because
   *  integrals over a point in 3D are zero... unless it's a singular point
   *  for now we'll assume this can't happen
   *
   * 	Found one - In a 1D shell, the point represents the element edge, for
   *  which a boundary term needs to be added to the element integral.
   **************************************************************************/

  /******************************************************************************/
  /**************************************************************************
   *                          BLOCK 6 - Weak Shifty Boundary Conditions
   *   This section is for boundary conditions which equate the BOUNDARY terms
   *   of two different residual equations, e.g. FLUID_SOLID and SOLID_FLUID
   *   type BC's.
   *   Do this before rotating mesh or momentum equations
   **************************************************************************/
  /******************************************************************************/
  if (first_elem_side_BC_array[ielem] != NULL) {
    elem_side_bc = first_elem_side_BC_array[ielem];

    /******************************************************************************/
    do { /* begining of do while construct
          * which loops over the sides of this element that have boundary
          * conditions */
      /******************************************************************************/
      /* Determine which types of bc's apply on this side */
      /*call_special = 0; */
      for (ibc = 0; (bc_input_id = (int)elem_side_bc->BC_input_id[ibc]) != -1; ibc++) {
        bct = BC_Types[bc_input_id].desc->method;
        /*if (bct == WEAK_SHIFT) call_special = 1; */
      } /* end of loop over boundary condition number */

      /* MMH notes: This next call was commented out when I updated
       * the BC's for LSA.  If someone wants to be able to do the
       * following call, with LSA, then they need to be careful.
       *
       * If uncommenting also uncomment other call_special references
       */
      /* if (call_special) shift_residuals_and_jacobians(ija, a, x, resid_vector, */
      /* 					      	first_elem_side_BC_array) */
      /****************************************************************************/
    } while ((elem_side_bc = elem_side_bc->next_side_bc) != NULL);
    /* END of do  while () construct				      */
    /******************************************************************************/
  } /* END if (First_Elem_Side_BC_Array[ielem] != NULL) 		      */
  /******************************************************************************/

  /**************************************************************************
   *                          BLOCK 7
   * rotate equations here before dirichlet or strong conditions
   * and do any shifting of residual and jacobian contributions
   **************************************************************************/
  /*
   * Make sure there are BC's at this node and that we have vector equations
   */
  /*if (  first_elem_side_BC_array[ielem] != NULL && */
  if ((pde[R_MOMENTUM1] || pde[R_MESH1]) && Num_ROT == 0) {
    if (ielem_dim == 2) {
      call_rotate = 0;
      /* determine if rotation is needed */
      for (i = 0; (call_rotate == 0) && (i < num_local_nodes); i++) {
        /*
         * To address a particular residual equation, map the local
         * elemental node number i into a processor node index I.
         * Then, make sure this node is owned by this processor
         */
        I = Proc_Elem_Connect[iconnect_ptr + i];
        if (I < (dpi->num_internal_nodes + dpi->num_boundary_nodes)) {
          if (in_list(I, 0, num_mom_rotate[pg->imtrx], mom_rotate_node[pg->imtrx]) != -1)
            call_rotate = 1;
          if (in_list(I, 0, num_mesh_rotate[pg->imtrx], mesh_rotate_node[pg->imtrx]) != -1)
            call_rotate = 1;
        }
      }
      /* MMH: I am assuming that apply_rotated_bc will work properly
       * in LSA b/c it doesn't "add" anything, only rearranges it,
       * and we supposedly have the correct contributions by now.
       */
      if (call_rotate) {
        if (Use_2D_Rotation_Vectors == TRUE)
          rotate_eqns_at_node_2D(iconnect_ptr, ielem_dim, num_local_nodes, ams);
        else if (first_elem_side_BC_array[ielem] != NULL) {
          err = apply_rotated_bc(resid_vector, first_elem_side_BC_array, ielem, ielem_type,
                                 num_local_nodes, ielem_dim, iconnect_ptr, num_total_nodes, exo);
          GOMA_EH(err, " apply_rotated_bc");
#ifdef CHECK_FINITE
          err = CHECKFINITE("apply_rotated_bc");
          if (err)
            return -1;
#endif
        }
      }
    }
  }

  if ((pde[R_MOMENTUM1] || pde[R_MESH1]) && Num_ROT > 0) {
    int id_mesh, id_mom; /* local temporary things */

    /* Make sure there are BC's at this node and that we have vector equations*/

    /* determine if rotation is needed */
    for (i = 0; i < num_local_nodes; i++) {
      id_mesh = ei[pg->imtrx]->ln_to_dof[MESH_DISPLACEMENT1][i];
      id_mom = ei[pg->imtrx]->ln_to_dof[VELOCITY1][i];
      /*
       *  To address a particular residual equation, map the local
       *   elemental node number i into a global index I
       *  You only have to rotate equations at nodes that you own?!?
       */
      I = Proc_Elem_Connect[iconnect_ptr + i];
      if (I < (dpi->num_internal_nodes + dpi->num_boundary_nodes)) {
        /*
         * Check for rotation of mesh equations or momentum equations
         * work on MESH first
         */
        rotate_mesh = 0;
        if (ROT_list[I] != NULL && (ROT_list[I][VECT_EQ_MESH] != -1) && pde[R_MESH1] &&
            id_mesh > -1)
          rotate_mesh = 1;

        rotate_momentum = 0;

        if (ROT_list[I] != NULL && (ROT_list[I][VECT_EQ_MOM] != -1) && pde[R_MOMENTUM1] &&
            id_mom > -1)
          rotate_momentum = 1;

        /* Now rotate the appropriate equations and sensitivities*/
        if (rotate_mesh) {
          /* Call mesh rotation routine */
          /* MMH: I am assuming that rotate_mesh_eqn will work
           * properly in LSA b/c it doesn't "add" anything, only
           * rearranges it, and we supposedly have the correct
           * contributions by now.  I am worried, though, with this
           * call b/c there is direct access to a[] in here, instead
           * of through lec->J's...
           */
          rotate_mesh_eqn(id_mesh, I, iconnect_ptr, ielem_dim, ams);
        }

        if (rotate_momentum) {
          /* Call momentum rotation routine */
          /* MMH: See comment above for rotate_mesh_eqn. */
          rotate_momentum_eqn(id_mom, I, iconnect_ptr, ielem_dim, ams);
        }
      }
    } /* end of loop over nodes */
  }   /* end of if Num_ROT > 0 */

  if (pde[R_MOMENTUM1] && goma_automatic_rotations.automatic_rotations) {
    /* determine if rotation is needed */
    for (i = 0; i < num_local_nodes; i++) {
      int id_mom = ei[pg->imtrx]->ln_to_dof[VELOCITY1][i];
      /*
       *  To address a particular residual equation, map the local
       *   elemental node number i into a global index I
       *  You only have to rotate equations at nodes that you own?!?
       */
      I = Proc_Elem_Connect[iconnect_ptr + i];
      if (I < (dpi->num_internal_nodes + dpi->num_boundary_nodes)) {

        int rotate_momentum = 0;

        if (goma_automatic_rotations.rotation_nodes[I]
                .eqn_is_rotated[vector_equation_from_equation(R_MOMENTUM1)] &&
            pde[R_MOMENTUM1]) {
          rotate_momentum = 1;
        }

        if (rotate_momentum) {
          rotate_momentum_auto(id_mom, I, ielem_dim, ams);
        }
      }
    } /* end of loop over nodes */
  }   /* end of if Num_ROT > 0 */
  if (pde[R_MESH1] && goma_automatic_rotations.automatic_rotations) {
    /* determine if rotation is needed */
    for (i = 0; i < num_local_nodes; i++) {
      int id_mesh = ei[pg->imtrx]->ln_to_dof[R_MESH1][i];
      /*
       *  To address a particular residual equation, map the local
       *   elemental node number i into a global index I
       *  You only have to rotate equations at nodes that you own?!?
       */
      I = Proc_Elem_Connect[iconnect_ptr + i];
      if (I < (dpi->num_internal_nodes + dpi->num_boundary_nodes)) {

        int rotate_mesh = 0;

        if (goma_automatic_rotations.rotation_nodes[I]
                .eqn_is_rotated[vector_equation_from_equation(R_MESH1)] &&
            pde[R_MESH1]) {
          rotate_mesh = 1;
        }

        if (rotate_mesh) {
          rotate_mesh_auto(id_mesh, I, ielem_dim, ams);
        }
      }
    } /* end of loop over nodes */
  }   /* end of if Num_ROT > 0 */

  /******************************************************************************/
  /*                              BLOCK 9                                       */
  /* Start of STRONG Boundary Conditions (other than Dirichlet). There are two  */
  /* types depending on whether they are integrated or pointwise collocation    */
  /******************************************************************************/
  /*
   *   initialize node list for fluid/solid boundaries
   */
  for (i = 0; i < num_local_nodes; i++) {
    local_node_list_fs[i] = -1;
  }

  /******************************************************************************/
  if (first_elem_side_BC_array[ielem] != NULL) {
    /****************************************************************************/
    /* We may want to figure out which equations are rotated and which boundary
     * conditions are applied in a strong sense (replacing equations as opposed
     * to the weak sense in which boundary conditions are added to the existing
     * equations), and do the rotation and zeroing (or subtracting local
     * contributions) now, before we get into the heavy duty boundary
     * condition stuff
     */
    elem_side_bc = first_elem_side_BC_array[ielem];

    /*****************************************************************************/
    do { /* begining of do while construct */
      /* which loops over the sides of this element that have boundary
         conditions */

      /*
       *  Set flags for subroutines to call for each boundary condition
       *  on this side
       */
      call_int = 0;
      call_col = 0;
      call_contact = 0.;
      int call_nedelec = 0.;
      for (ibc = 0; (bc_input_id = (int)elem_side_bc->BC_input_id[ibc]) != -1; ibc++) {
        bct = BC_Types[bc_input_id].desc->method;
        if (bct == STRONG_INT_SURF)
          call_int = 1;
        if (bct == STRONG_INT_NEDELEC)
          call_nedelec = 1;
        if (bct == COLLOCATE_SURF)
          call_col = 1;
        if (bct == CONTACT_SURF)
          call_contact = 1;
      }
      /*
       * Major change here 6/10/98 to accomodate frontal solver.  Here the
       * FLUID_SOLID/SOLID_FLUID BCs actually use local element contribution
       * from the lec structure.  So, before the apply_integrated_bc was called
       * first and usually a companion NO_SLIP condition obliterated the fluid
       * momentum contributions to lec and so when the fluid-solid condition
       * followed, there were no stresses left to balance.  HEre we simply
       * switched the order of the calls below (i.e. call_col befor call_int)
       * to avoid this.  Cross your fingers!
       */

      /*
       * BOUNDARY CONDITIONS OF TYPE 1 - Pointwise Collocation at the
       * nodal points
       */
      if (call_col) {
        err = apply_point_colloc_bc(resid_vector, delta_t, theta, ielem, ip_total, ielem_type,
                                    num_local_nodes, ielem_dim, iconnect_ptr, elem_side_bc,
                                    num_total_nodes, local_node_list_fs, time_value, exo);
        GOMA_EH(err, " apply_point_colloc_bc");
#ifdef CHECK_FINITE
        err = CHECKFINITE("apply_point_colloc_bc");
        if (err)
          return -1;
#endif
        if (neg_elem_volume)
          return -1;
      }

      /*
       * BOUNDARY CONDITIONS OF TYPE 2 - Integrated over the surface
       * - only do the strong integrated conditions here
       */
      if (call_int) {
        err = apply_integrated_bc(x, resid_vector, delta_t, theta, &pg_data, ielem, ielem_type,
                                  num_local_nodes, ielem_dim, iconnect_ptr, elem_side_bc,
                                  num_total_nodes, STRONG_INT_SURF, time_value, element_search_grid,
                                  exo);
        GOMA_EH(err, " apply_integrated_bc");
#ifdef CHECK_FINITE
        err = CHECKFINITE("apply_integrated_bc");
        if (err)
          return -1;
#endif
        /*printf("Element: %d, ID_side: %d \n", ei[pg->imtrx]->ielem, elem_side_bc->id_side );
          for(i3=0; i3< (int)  elem_side_bc->num_nodes_on_side; i3++)
          {
          id3 = (int) elem_side_bc->local_elem_node_id[i3]; I3 =  I = Proc_Elem_Connect[iconnect_ptr
          + id3]; printf("\tI: %d, R[0][%d]: %7.4g, R[1][%d]: %7.4g\n", I3, id3, lec->R[0][id3],
          id3,  lec->R[1][id3]);
          }*/

        if (neg_elem_volume)
          return -1;
      }

      /*
       * BOUNDARY CONDITIONS which don't really fit in with the above types
       */
      if (call_int || call_col) {
        /* add call to Gibb's inequality condition is evaluated
         * within apply_special_bc, and potentially contact lines
         * are suddenly fixed here
         */
        err =
            apply_special_bc(ams, x, resid_vector, x_old, x_older, xdot, xdot_old, delta_t, theta,
                             first_elem_side_BC_array, ielem, ip_total, ielem_type, num_local_nodes,
                             ielem_dim, iconnect_ptr, elem_side_bc, num_total_nodes, SPECIAL, CA_id,
                             CA_fselem, CA_sselem, exo, time_value);
        GOMA_EH(err, " apply_special_bc");
#ifdef CHECK_FINITE
        err = CHECKFINITE("apply_special_bc");
        if (err)
          return -1;
#endif
        if (neg_elem_volume)
          return -1;
      }

      if (call_contact) {
        struct Level_Set_Data *ls_save = ls;

        /* Special contact conditions are applied here. These are boundary
           conditions that require a contact search or a mesh search, such
           as that which might happen with overlapping grids.   */

        /* Cludge her for BCs that are tied to R_PHASE0 */
        if (pfd != NULL) {
          ls = pfd->ls[0];
        } else {
          GOMA_EH(GOMA_ERROR,
                  "YOU cannot apply CONTACT_SURF BCs in mm_names.h with FILL field. R_PHASE only");
        }

        err = apply_contact_bc(x, resid_vector, delta_t, theta, pg_data.h_elem_avg, pg_data.h,
                               pg_data.mu_avg, pg_data.U_norm, first_elem_side_BC_array, ielem,
                               ielem_type, num_local_nodes, ielem_dim, iconnect_ptr, elem_side_bc,
                               num_total_nodes, CONTACT_SURF, -1, NULL, NULL, NULL, NULL,
                               time_value, exo);
        GOMA_EH(err, " apply_contact_bc");
#ifdef CHECK_FINITE
        err = CHECKFINITE("apply_contact_bc");
        if (err)
          return -1;
#endif
        if (neg_elem_volume)
          return -1;

        ls = ls_save;
      }
      if (call_nedelec) {
        err = apply_nedelec_bc(x, resid_vector, delta_t, theta, &pg_data, ielem, ielem_type,
                               num_local_nodes, ielem_dim, iconnect_ptr, elem_side_bc,
                               num_total_nodes, STRONG_INT_NEDELEC, time_value, element_search_grid,
                               exo);
        GOMA_EH(err, " apply_integrated_bc");
#ifdef CHECK_FINITE
        err = CHECKFINITE("apply_integrated_bc");
        if (err)
          return -1;
#endif
        /*printf("Element: %d, ID_side: %d \n", ei[pg->imtrx]->ielem, elem_side_bc->id_side );
          for(i3=0; i3< (int)  elem_side_bc->num_nodes_on_side; i3++)
          {
          id3 = (int) elem_side_bc->local_elem_node_id[i3]; I3 =  I = Proc_Elem_Connect[iconnect_ptr
          + id3]; printf("\tI: %d, R[0][%d]: %7.4g, R[1][%d]: %7.4g\n", I3, id3, lec->R[0][id3],
          id3,  lec->R[1][id3]);
          }*/

        if (neg_elem_volume)
          return -1;
      }
      /****************************************************************************/
    } while ((elem_side_bc = elem_side_bc->next_side_bc) != NULL);
    /* END of do  while () construct				      */
    /******************************************************************************/
  } /* END if (First_Elem_Side_BC_Array[ielem] != NULL) 		      */
  /******************************************************************************/

  /******************************************************************************/
  if (First_Elem_Edge_BC_Array[pg->imtrx][ielem] != NULL) {
    /******************************************************************************/
    elem_edge_bc = First_Elem_Edge_BC_Array[pg->imtrx][ielem];

    /****************************************************************************/
    do { /* begining of do while construct */
      /* which loops over the sides of this element that have boundary
         conditions */
      /****************************************************************************/

      /*
       *  Set flags for subroutines to call for each boundary
       *  condition on this side
       */
      call_int = 0;
      call_col = 0;
      for (ibc = 0; (bc_input_id = (int)elem_edge_bc->BC_input_id[ibc]) != -1; ibc++) {
        bct = BC_Types[bc_input_id].desc->method;
        if (bct == STRONG_INT_EDGE)
          call_int = 1;
        if (bct == COLLOCATE_EDGE)
          call_col = 1;
      }

      /*
       * BOUNDARY CONDITIONS OF TYPE 1 - Integrated over the surface
       * - only do the strong integrated conditions here
       */
      if (call_int) {
        err = apply_integrated_curve_bc(x, resid_vector, delta_t, theta, ielem, ielem_type,
                                        num_local_nodes, ielem_dim, iconnect_ptr, elem_edge_bc,
                                        num_total_nodes, STRONG_INT_EDGE, exo);
        GOMA_EH(err, " apply_integrated_curve_bc");
#ifdef CHECK_FINITE
        err = CHECKFINITE("apply_integrated_curve_bc");
        if (err)
          return -1;
#endif
        if (neg_elem_volume)
          return -1;
      }

      /*
       * BOUNDARY CONDITIONS OF TYPE 2 - Pointwise Collocation at the
       * nodal points
       */
      if (call_col) {
        err = apply_point_colloc_edge_bc(x, x_old, x_older, xdot, xdot_old, resid_vector, delta_t,
                                         theta, ielem, ip_total, ielem_type, num_local_nodes,
                                         ielem_dim, iconnect_ptr, elem_edge_bc, num_total_nodes,
                                         local_node_list_fs, time_value);
        GOMA_EH(err, " apply_point_colloc_bc");
#ifdef CHECK_FINITE
        err = CHECKFINITE("apply_point_colloc_bc");
        if (err)
          return -1;
#endif
        if (neg_elem_volume)
          return -1;
      }

      /****************************************************************************/
    } while ((elem_edge_bc = elem_edge_bc->next_edge_bc) != NULL);
    /* END of do  while () construct				      */
    /******************************************************************************/
  } /* END if (First_Elem_Edge_BC_Array[ielem] != NULL) 		      */

  /********************************************************************************
   *                          BLOCK 8
   * Now zero out the rows that correspond to Dirichlet Conditions
   *  i.e. conditions where the value of the solution vector is known and
   *       so the residual is set to zero and jacobian has a one on the diagonal
   *       (sometimes one is not on the diagonal, as in the case
   *        of distinguishing conditions)
   ********************************************************************************/

  err = put_dirichlet_in_matrix(x, num_total_nodes);
  GOMA_EH(err, " put_dirichlet_in_matrix");
#ifdef CHECK_FINITE
  err = CHECKFINITE("put_dirichlet_in_matrix");
  if (err)
    return -1;
#endif

  /* PRS test code for shell endpoint conditions.  4/21/2004 */
  if (pde[R_MESH1] && pde[R_SHELL_TENSION]) {
    for (i = 0; i < num_local_nodes; i++) {
      I = Proc_Elem_Connect[iconnect_ptr + i];
      node = Nodes[I];
      /*Now apply condition */
      if (node->DBSH_SLOPE_X == 1) {
        eqn = R_MESH1;
        lec->R[LEC_R_INDEX(upd->ep[pg->imtrx][eqn], i)] = 1.0 * BIG_PENALTY;
      }
      if (node->DBSH_SLOPE_Y == 1) {
        eqn = R_MESH2;
        lec->R[LEC_R_INDEX(upd->ep[pg->imtrx][eqn], i)] = 0.0 * BIG_PENALTY;
      }
    }
  }

  if ((Linear_Stability == LSA_3D_OF_2D || Linear_Stability == LSA_3D_OF_2D_SAVE) &&
      !(af->Assemble_LSA_Jacobian_Matrix) && !(af->Assemble_LSA_Mass_Matrix))
    /* Cheap, but effective!  Zero out the equation rows for
     * R_MOMENTUM3, and the columns for sensitivities to VELOCITY3,
     * then stick in a Dirichlet condition.  The end effect should be
     * that our original 2D steady-state problem has an appended
     * system equivalent to the identity problem, "I x = 0". */
    for (i = 0; i < lec->max_dof; i++) {
      eqn = upd->ep[pg->imtrx][R_MOMENTUM3];
      var = upd->vp[pg->imtrx][VELOCITY3];
      zero_lec_row(lec->J, eqn, i);
      zero_lec_column(lec->J, var, i);
      lec->J[LEC_J_INDEX(eqn, var, i, i)] = 1.0;
      lec->R[LEC_R_INDEX(eqn, i)] = 0.0;
    }
  /* If we're doing normal mode LSA, and we're on wavenumber 0, then
   * the jacobian contributions to the w-momentum equation are
   * (rightly) all zero.  This causes a zero row, which we don't
   * really want.  We want to emulate the Dirichlet = 0 "boundary
   * condition" for all internal nodes, so we put a 1 on the diagonal
   * for the same reason(s) as above.  We don't want this to be
   * subtracted out, so we only apply it on one pass. */
  if ((Linear_Stability == LSA_3D_OF_2D || Linear_Stability == LSA_3D_OF_2D_SAVE) &&
      af->Assemble_LSA_Jacobian_Matrix && LSA_3D_of_2D_wave_number == 0 && LSA_3D_of_2D_pass == 1)
    for (i = 0; i < lec->max_dof; i++)
      lec->J[LEC_J_INDEX(R_MOMENTUM3, VELOCITY3, i, i)] = 1.0;

  if (ls != NULL && ls->Ignore_F_deps) {
    int peqn, pvar;
    for (eqn = V_FIRST; eqn < V_LAST; eqn++) {
      peqn = upd->ep[pg->imtrx][eqn];
      if (peqn != -1 && eqn != R_FILL) {
        var = FILL;
        pvar = upd->vp[pg->imtrx][var];
        for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] = 0.;
            /*
              if ( lec->J[LEC_J_INDEX(peqn,pvar,i,j)] != 0. )
              DPRINTF(stderr,"lec->J[eqn=%d][var=FILL][%d][%d] =
              %g\n",eqn,i,j,lec->J[LEC_J_INDEX(peqn,pvar,i,j)]);
            */
          }
        }
      }
    }
  }
  if (pfd != NULL && pfd->ls[0]->Evolution == LS_EVOLVE_SLAVE) {
    int peqn, pvar;
    for (eqn = V_FIRST; eqn < V_LAST; eqn++) {
      peqn = upd->ep[pg->imtrx][eqn];
      if (peqn != -1 && eqn != R_PHASE1) {
        var = PHASE1;
        pvar = upd->vp[pg->imtrx][var];
        for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] = 0.;
          }
        }
      }
    }
  }

  /*
   * Load local element stiffness matrix (lec) into global matrix, depending
   * on solver used - front vs. everything else, and matrix storage format
   * MSR vs VBR.
   */
#if 0
  if (first_elem_side_BC_array[ielem] != NULL) 
    {
      elem_side_bc = first_elem_side_BC_array[ielem];
	
      do {
	int id,jd;
	for (ibc = 0; 
	     (bc_input_id = (int) elem_side_bc->BC_input_id[ibc]) != -1; ibc++) 
	  {
	    if(BC_Types[bc_input_id].desc->method == STRONG_INT_SURF ) {
				
	      printf("\nElement: %d, ID_side: %d\n", ei[pg->imtrx]->ielem,  elem_side_bc->id_side);
				
	      for( i=0; i < (int) elem_side_bc->num_nodes_on_side; i++)
		{	
		  id = (int) elem_side_bc->local_elem_node_id[i];
		  I = Proc_Elem_Connect[iconnect_ptr + id];
					
		  for(eqn=0; eqn < 2; eqn++)
		    {
		      if( id != ei[pg->imtrx]->ln_to_first_dof[eqn][id] ) { printf("Error found elem: %d, I:%d, eqn %d\n",  ei[pg->imtrx]->ielem, I, eqn); }
		      printf("%d \t", I );
		      for( var=0; var< 2; var++) 
			{
			  for( j=0; j < (int) elem_side_bc->num_nodes_on_side; j++)
			    {	
			      jd = (int) elem_side_bc->local_elem_node_id[j];		
			      printf("J[%d][%d][%d][%d]: %9.3g ", eqn, var, id, jd, lec->J[eqn][var][id][jd] );
			    }
			}
                      printf("\tR[%d][%d]: %7.4g\n", eqn, id, lec->R[LEC_R_INDEX(eqn,id)] );
		    }
					
		}
	    }
	  }
      }
      while ( (elem_side_bc =  elem_side_bc->next_side_bc) != NULL );
    }
#endif

  load_lec(exo, ielem, ams, x, resid_vector, estifm);

  /*  if( pfd != NULL && pfd->Use_Constraint == TRUE )
      {
      load_pf_constraint(  pf_constraint, d_pf_lm, d_lm_pf );
      }
  */

  /******************************************************************************/
  /*                              BLOCK 13                                      */
  /*                            CLEANUP BLOCK                                   */
  /******************************************************************************/

  /*
   *  Good spot to print out the element stiffness matrix and residual or to
   *  free any memory that was allocated on the element level
   */
  mm_fill_end = ut();
  mm_fill_total += mm_fill_end - mm_fill_start;

  /*
   *  Free up Subgrid integration arrays
   */

  if (Subgrid_Int.active) {
    safe_free((void *)Subgrid_Int.s);
    Subgrid_Int.s = NULL;
    safe_free((void *)Subgrid_Int.wt);
    Subgrid_Int.wt = NULL;
    free_search_grid(&element_search_grid);
  }

  if (Debug_Flag > 1) {
    P0PRINTF("%s: ends\n", yo);
    MMH_ip = -1;
  }

  if (Proc_NS_List_Length > 0 &&
      ((Linear_Solver != FRONT && ielem == exo->eb_ptr[exo->num_elem_blocks] - 1) ||
       (Linear_Solver == FRONT && ielem == exo->elem_order_map[exo->num_elem_blocks] - 1))) {
    if (zeroCA == 0) {
      int count = 0, Num_CAs_done = 0;
      for (j = 0; j < MAX_CA; j++) {
        if (CA_id[j] == -2)
          Num_CAs_done++;
        if (CA_proc[j] == ProcID)
          count++;
      }

      if (count != Num_CAs_done) {
        GOMA_WH(-1, "\nNot all contact angle conditions were applied!\n");
        for (j = 0; j < count; j++) {
          fprintf(stderr, "CA:%d ID:%d fselem:%d sselem:%d Proc:%d\n", j, CA_id[j], CA_fselem[j],
                  CA_sselem[j], CA_proc[j]);
        }
        fprintf(stderr, "Count=%d  Done=%d\n", count, Num_CAs_done);
      }
    }
  }

  return 0;
} /*   END OF matrix_fill                                                     */
/******************************************************************************/

/* matrix_fill_stress is called from numerical_jacobian_compute_stress function
 * and is mainly used for the log-conformation formulation for viscoelastic stress.
 * It is a copy of the matrix_fill function, but only calls assemble_stress.
 */
int matrix_fill_stress(struct GomaLinearSolverData *ams,
                       double x[],            /* Solution vector */
                       double resid_vector[], /* Residual vector */
                       double x_old[],        /* Solution vector , previous last time step */
                       double x_older[],      /* Solution vector , previous prev time step */
                       double xdot[],         /* xdot of current solution                  */
                       double xdot_old[],     /* xdot_old of current solution              */
                       double x_update[],     /* last update vector */
                       double *ptr_delta_t,   /* current time step size                    */
                       double *ptr_theta,     /* parameter to vary time integration from
                                               * explicit (theta = 1) to implicit
                                               * (theta = 0) */
                       struct elem_side_bc_struct *first_elem_side_BC_array[],
                       /* This is an array of pointers to the first
                          surface integral defined for each element.
                          It has a length equal to the total number
                          of elements defined on current processor */
                       double *ptr_time_value,
                       Exo_DB *exo,              /* ptr to EXODUS II finite element mesh db */
                       Dpi *dpi,                 /* ptr to distributed processing info */
                       int *ptr_ielem,           /* element number */
                       int *ptr_num_total_nodes, /* Number of nodes that each processor is
                                                    responsible for                              */
                       dbl *ptr_h_elem_avg,      /* global average element size for PSPG */
                       dbl *ptr_U_norm,          /* global average velocity for PSPG */
                       dbl *estifm,              /* element stiffness Matrix for frontal solver*/
                       int zeroCA)               /* boolean to zero zeroCA */

{
  double delta_t, theta, time_value; /*see arg list */
  int ielem, num_total_nodes;

#if 0
  int status = 0;		/* variable describing status of the matrix */
				/* calculation this is set to -1 if further */
				/* iteration is not recommended*/
#endif

  static int CA_id[MAX_CA];     /*  array of CA conditions */
  static int CA_fselem[MAX_CA]; /*  array of CA free surface elements  */
  static int CA_sselem[MAX_CA]; /*  array of CA solid surface elements */
  static int CA_proc[MAX_CA];   /*  Processor which has each CA */

  int mn;                     /* material block counter */
  int err;                    /* temp variable to hold diagnostic flags.      */
  int ip;                     /* ip is the local quadrature point index       */
  int ip_total;               /* ip_total is the total number of volume
                                 quadrature points in the element             */
  int j;                      /* local index loop counter 	             */
  int i;                      /* Index for the local node number - row        */
  int b, eqn = -1, e, lm_dof; /*convenient params for neatness */
  int I;                      /* Index for global node number - row        */
  int ielem_type = 0;         /* Element type of the current element          */
  int ielem_dim;              /* Element physical dimension                   */
  int num_local_nodes;        /* Number of local basis functions in the
                                 current element                              */
  int iconnect_ptr;           /* Pointer to the beginning of the connectivity
                                 list for the current element                 */
  int ibc;                    /* Index for the boundary condition 	      */
  int var = -1;               /* variable name (TEMPERATURE, etc) */

  double s, t, u; /* Gaussian-quadrature point locations          */

  int call_int, call_col, call_contact, call_shell_grad, call_sharp_int;
  /* int call_special; */
  /* flags for calling boundary
     condition routines */
  int call_rotate;
  int rotate_mesh, rotate_momentum;
  int make_trivial = 0; /* Flags for Eulerian solid mechanics
                                    and level-set */
  int bct, mode;

  /* ___________________________________________________________________________*/
  /* ___________________________________________________________________________*/

  double xi[DIM]; /* Local element coordinates of Gauss point. */

  double wt = 0.0; /* Quadrature weights
                    units - ergs/(sec*cm*K) = g*cm/(sec^3*K)     */

  struct Petrov_Galerkin_Data pg_data;

  struct elem_side_bc_struct *elem_side_bc;
  /* Pointer to an element side boundary condition
     structure				      */
  struct elem_edge_bc_struct *elem_edge_bc;
  /* Pointer to an element edge boundary condition
     structure				      */
  int bc_input_id; /* Input ID of the surf_int boundary condition  */
  int id_side;     /* side counter for discontinuous Galerkin method */

  /* List to keep track of nodes at which residual and
     Jacobian corrections have been made           */
  /* for mesh equations */
  /* List to keep track of nodes at which residual and
     Jacobian corrections have been made           */
  /* for mesh equations */
  int local_node_list_fs[MDE]; /* list to keep track of nodes at which solid contributions
                                  have been transfered to liquid (fluid-solid boundaries)*/

  int discontinuous_stress; /* flag that tells you if you are doing Discontinuous Galerkin
                               for the species equations */
  int ielem_type_mass = -1; /* flag to give discontinuous interpolation type */

  int pspg_local = 0;

  NODE_INFO_STRUCT *node;
  SGRID *element_search_grid = NULL;

  struct Level_Set_Data *ls_old;

  static double mm_fill_start, mm_fill_end; /* Count CPU time this call. */

  static char yo[] = "matrix_fill"; /* My name to take blame... */

  int *pde;

  /*
   * BEGINNING OF EXECUTABLE STATEMENTS
   */

  mm_fill_start = ut();

  /*
   * Unpack some pointers that were added to help portability of front
   */
  delta_t = *ptr_delta_t;
  theta = *ptr_theta;
  time_value = *ptr_time_value;
  ielem = *ptr_ielem;
  num_total_nodes = *ptr_num_total_nodes;
  pg_data.U_norm = *ptr_U_norm;

  if (Debug_Flag > 1) {
    P0PRINTF("%s: begins\n", yo);
  }

  /******************************************************************************/
  /*                                BLOCK 1                                     */
  /*          LOOP OVER THE ELEMENTS DEFINED ON THE CURRENT PROCESSOR           */
  /*          INITIALIZATION THAT IS DEPENDENT ON THE ELEMENT NUMBER            */
  /******************************************************************************/

  if (Proc_NS_List_Length > 0 &&
      ((zeroCA == 1) || ((Linear_Solver != FRONT && ielem == exo->eb_ptr[0]) ||
                         (Linear_Solver == FRONT && ielem == exo->elem_order_map[0] - 1)))) {
    int nsp, nspk, count = -1;
    memset(CA_fselem, -1, sizeof(int) * MAX_CA);
    memset(CA_sselem, -1, sizeof(int) * MAX_CA);
    memset(CA_id, -1, sizeof(int) * MAX_CA);
    memset(CA_proc, -1, sizeof(int) * MAX_CA);
    for (j = 0; j < Num_BC; j++) {
      switch (BC_Types[j].BC_Name) {
      case CA_BC:
      case CA_MOMENTUM_BC:
      case VELO_THETA_HOFFMAN_BC:
      case VELO_THETA_TPL_BC:
      case VELO_THETA_COX_BC:
      case VELO_THETA_SHIK_BC:
        nsp = match_nsid(BC_Types[j].BC_ID);
        nspk = Proc_NS_List[Proc_NS_Pointers[nsp]];
        if (nsp != -1 && Nodes[nspk]->Proc == ProcID) {
          count++;
          CA_proc[count] = ProcID;
        }
        break;
      }
    }

    /*
     * Initialize the accumulated CPU time for assembly...
     */
    mm_fill_total = 0;
  }

  err = load_elem_dofptr(ielem, exo, x, x_old, xdot, xdot_old, 0);
  GOMA_EH(err, "load_elem_dofptr");

  err = bf_mp_init(pd);
  mn = ei[pg->imtrx]->mn;
  pde = (int *)pd->e;

  for (mode = 0; mode < vn->modes; mode++) {
    ve[mode] = ve_glob[mn][mode];
  }

  /* discontinuous Galerkin information */
  /* element type, assumed to be globally consistent for now */

  /* need to take this section out and up above.  Inefficient now as it assumes
     every element is different */

  discontinuous_stress = 0;

  if (pd->i[pg->imtrx][POLYMER_STRESS11] == I_P1) {
    discontinuous_stress = 1;
  } else if (pd->i[pg->imtrx][POLYMER_STRESS11] == I_P0) {
    discontinuous_stress = 1;
  } else if (pd->i[pg->imtrx][POLYMER_STRESS11] == I_PQ1) {
    if (pd->Num_Dim == 3)
      GOMA_EH(GOMA_ERROR, "Sorry PQ1 interpolation has not been implemented in 3D yet.");
    discontinuous_stress = 1;
  } else if (pd->i[pg->imtrx][POLYMER_STRESS11] == I_PQ2) {
    if (pd->Num_Dim == 3)
      GOMA_EH(GOMA_ERROR, "Sorry PQ2 interpolation has not been implemented in 3D yet.");
    discontinuous_stress = 1;
  }

  ielem_type = ei[pg->imtrx]->ielem_type; /* element type */

  num_local_nodes = ei[pg->imtrx]->num_local_nodes; /* number of local  basis functions */

  ielem_dim = ei[pg->imtrx]->ielem_dim; /* physical dimension  of this element */

  iconnect_ptr =
      ei[pg->imtrx]
          ->iconnect_ptr; /* find pointer to beginning  of this element's connectivity list */

  Subgrid_Int.active = FALSE;

  /*
   * Clean out local element contribution (lec) accumulator...
   */
  zero_lec();

  /* XFEM setup */
  /* DRN-now using check_xfem_contribution to do this */

  /******************************************************************************/
  /*                              BLOCK 1.5                                     */
  /*             INITIAL CALCULATIONS AT CENTROID OF THE ELEMENT                */
  /******************************************************************************/

  /*
   * set element size variable to zero if not used -- just to be safe
   */

  memset(pg_data.h, 0, sizeof(double) * DIM);
  memset(pg_data.hh, 0, sizeof(double) * DIM * DIM);
  memset(pg_data.dh_dxnode, 0, sizeof(double) * MDE * DIM);
  memset(pg_data.hsquared, 0, sizeof(double) * DIM);
  memset(pg_data.hhv, 0, sizeof(double) * DIM * DIM);
  memset(pg_data.dhv_dxnode, 0, sizeof(double) * MDE * DIM);
  memset(pg_data.v_avg, 0, sizeof(double) * DIM);
  memset(pg_data.dv_dnode, 0, sizeof(double) * MDE * DIM);
  pg_data.mu_avg = 0.;
  pg_data.rho_avg = 0.;

  if (PSPG) {
    if (PSPG == 1) {
      pspg_local = 0;
    }
    /* This is the flag for the standard local PSPG */
    else if (PSPG == 2) {
      pspg_local = 1;
    }
  }

  if (PSPG && pde[R_PRESSURE] && pde[R_MOMENTUM1]) {
    xi[0] = 0.0;
    xi[1] = 0.0;
    xi[2] = 0.0;
    (void)load_basis_functions(xi, bfd);
    pg_data.mu_avg = element_viscosity();
    pg_data.rho_avg = density(NULL, time_value);

    if (pspg_local) {
      h_elem_siz(pg_data.hsquared, pg_data.hhv, pg_data.dhv_dxnode, pde[R_MESH1]);
      element_velocity(pg_data.v_avg, pg_data.dv_dnode, exo);
    }
  }

  if (Cont_GLS && pde[R_PRESSURE] && pde[R_MOMENTUM1]) {
    xi[0] = 0.0;
    xi[1] = 0.0;
    xi[2] = 0.0;
    (void)load_basis_functions(xi, bfd);
    pg_data.mu_avg = element_viscosity();
    pg_data.rho_avg = density(NULL, time_value);

    if (Cont_GLS == 2) {
      h_elem_siz(pg_data.hsquared, pg_data.hhv, pg_data.dhv_dxnode, pde[R_MESH1]);
      element_velocity(pg_data.v_avg, pg_data.dv_dnode, exo);
    }
  }

  /*
   * The current SUPG model requires a value for the element's size and the
   * average velocity. Evaluate this here.
   */
  if ((vn->wt_funcModel == SUPG && pde[R_STRESS11] && pde[R_MOMENTUM1]) ||
      (mp->Spwt_funcModel == SUPG && pde[R_MASS] && (pde[R_MOMENTUM1] || pde[R_MESH1])) ||
      (mp->Mwt_funcModel == SUPG && pde[R_MOMENTUM1]) ||
      (mp->Ewt_funcModel == SUPG && pde[R_ENERGY] && (pde[R_MOMENTUM1] || pde[R_MESH1])) ||
      (mp->Ewt_funcModel == SUPG && pde[R_SHELL_ENERGY] && (pde[R_LUBP]))) {
    h_elem_siz(pg_data.hsquared, pg_data.hhv, pg_data.dhv_dxnode, pde[R_MESH1]);
    element_velocity(pg_data.v_avg, pg_data.dv_dnode, exo);
  }

  if (cr->MassFluxModel == HYDRODYNAMIC) {
    /* For shock capturing diffusivity in Phillips model */
    h_elem_siz(pg_data.hsquared, pg_data.hhv, pg_data.dhv_dxnode, pde[R_MESH1]);
  }
  if (mp->DiffusivityModel[0] == HYDRO && mp->QTensorDiffType[0] == CONSTANT) {
    h_elem_siz(pg_data.hsquared, pg_data.hhv, pg_data.dhv_dxnode, pde[R_MESH1]);
    assemble_qtensor(pg_data.hsquared);
  }

  if (cr->MassFluxModel == HYDRODYNAMIC_QTENSOR_OLD) {
    h_elem_siz(pg_data.hsquared, pg_data.hhv, pg_data.dhv_dxnode, pde[R_MESH1]);
    assemble_qtensor(pg_data.hsquared);
  }

  /* Ryan's new qtensor model */
  else if (mp->DiffusivityModel[0] == HYDRO && mp->QtensorExtensionPModel == CONSTANT) {
    h_elem_siz(pg_data.hsquared, pg_data.hhv, pg_data.dhv_dxnode, pde[R_MESH1]);
    assemble_new_qtensor(pg_data.hsquared);
  }

  ls_old = ls;
  if (!pde[R_MESH1]) {
    if (pfd != NULL)
      ls = pfd->ls[0]; /* this is a major hack */
    if (pde[R_LAGR_MULT1] && !ls->elem_overlap_state) {
      for (b = 0; b < pd->Num_Dim; b++) {
        for (i = 0; i < ei[pg->imtrx]->dof[R_LAGR_MULT1 + b]; i++) {
          if (af->Assemble_Residual) {
            lm_dof = ei[pg->imtrx]->gun_list[R_LAGR_MULT1 + b][i];
            eqn = upd->ep[pg->imtrx][R_LAGR_MULT1 + b];
            var = upd->vp[pg->imtrx][LAGR_MULT1 + b];
            lec->R[LEC_R_INDEX(eqn, i)] = x[lm_dof];
          }

          if (af->Assemble_Jacobian) {
            zero_lec_row(lec->J, eqn, i);
            lec->J[LEC_J_INDEX(eqn, var, i, i)] = 1.0;
          }
        }
      }
    }
  } else {
    if (pde[R_LAGR_MULT1]) {
      for (b = 0; b < pd->Num_Dim; b++) {
        for (i = 0; i < ei[pg->imtrx]->dof[R_LAGR_MULT1 + b]; i++) {
          if (af->Assemble_Residual) {
            eqn = upd->ep[pg->imtrx][R_LAGR_MULT1 + b];
            var = upd->vp[pg->imtrx][LAGR_MULT1 + b];
            lec->R[LEC_R_INDEX(eqn, i)] = 0.;
          }

          if (af->Assemble_Jacobian) {
            zero_lec_row(lec->J, eqn, i);
            lec->J[LEC_J_INDEX(eqn, var, i, i)] = 1.e-15;
          }
        }
      }
    }
  }
  ls = ls_old; /*Make things right again */

  /*
   * Here we will do another quick element-level manipulation for Eulerian Solid
   * mechanics problems using the level-set method and the assemble_real_solid
   * TALE routine.   In this situation we have a fixed mesh (hence no mesh motion
   * equations), and real-solid motion through the mesh. When not in the solid region
   * as delineated by the level-set field we will trivialize the real-solid equations
   * by setting the displacements to zero.  On elements that connect to elements through
   * which the level-set interface cuts, we will not assemble AT ALL, so as to achieve
   * the natural stress condition on the intefacial elements.
   */
  /* assemble fake equation for elements off interface */

  if (pde[R_SOLID1]) {
    if (pd->etm[pg->imtrx][R_SOLID1][(LOG2_MASS)])
      GOMA_EH(GOMA_ERROR, "Cannot do real inertia for TALE yet. Remove this line if trying to "
                          "perform EULERIAN solid mechanics");
    eqn = R_SOLID1;
    if (pd->TimeIntegration != STEADY && pd->etm[pg->imtrx][R_SOLID1][(LOG2_MASS)] &&
        !ls->elem_overlap_state && *esp->F[0] >= 0.0) {
      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        make_trivial = TRUE;
        if ((pd->i[pg->imtrx][eqn] == I_Q1) || (pd->i[pg->imtrx][eqn] == I_Q2)) {
          /* check all neighboring elements to see if any
             span the interface */
          I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
          for (j = exo->node_elem_pntr[I]; j < exo->node_elem_pntr[I + 1]; j++) {
            e = exo->node_elem_list[j];
            if (elem_overlaps_interface(e, x, exo, ls->Length_Scale)) {
              make_trivial = FALSE;
            }
          }
        }
      }
      if (make_trivial) {
        for (b = 0; b < pd->Num_Dim; b++) {
          for (i = 0; i < ei[pg->imtrx]->dof[R_SOLID1 + b]; i++) {
            if (af->Assemble_Residual) {
              eqn = upd->ep[pg->imtrx][R_SOLID1 + b];
              var = upd->vp[pg->imtrx][R_SOLID1 + b];
              lec->R[LEC_R_INDEX(eqn, i)] = 0.;
            }

            if (af->Assemble_Jacobian) {
              zero_lec_row(lec->J, eqn, i);
              lec->J[LEC_J_INDEX(eqn, var, i, i)] = 1.e-15;
            }
          }
        }
      }
    }
  }

  /******************************************************************************/
  /*                              BLOCK 1.6                                     */
  /*                Precalculations above the Gauss point loop                  */
  /******************************************************************************/

  if (pde[R_POR_LIQ_PRES]) {
    if (mp->Porous_Mass_Lump) {
      load_nodal_porous_properties(theta, delta_t);
    }
  }

  if ((pde[R_SHELL_SAT_OPEN]) || (pde[R_SHELL_SAT_OPEN_2])) {
    if (pde[R_SHELL_SAT_OPEN]) {
      load_nodal_shell_porous_properties(theta, delta_t, R_SHELL_SAT_OPEN);
    } else if (pde[R_SHELL_SAT_OPEN_2]) {
      load_nodal_shell_porous_properties(theta, delta_t, R_SHELL_SAT_OPEN_2);
    }
  }

  ip_total = elem_info(NQUAD, ielem_type);

  /* Loop over all the Volume Quadrature integration points */

  for (ip = 0; ip < ip_total; ip++) {
    MMH_ip = ip;

    /* case 1: normal gauss integration */
    find_stu(ip, ielem_type, &s, &t, &u); /* find quadrature point */
    wt = Gq_weight(ip, ielem_type);       /* find quadrature weights for */

    /* locations for current ip */

    /*
     * Local element coordinates...
     */

    xi[0] = s;
    xi[1] = t;
    xi[2] = u;

    fv->wt = wt;

    /*
     *	    PASS wt to assembly routines!!!!! ...
     */

    /*
     * Load up basis function information for each basis function
     * type needed by the current material. This is done in terms
     * of local element coordinates.
     */

    err = load_basis_functions(xi, bfd);
    GOMA_EH(err, "problem from load_basis_functions");

    err = beer_belly();
    GOMA_EH(err, "beer_belly");
    if (neg_elem_volume)
      return -1;
    if (zero_detJ)
      return -1;

    err = load_bf_grad();
    GOMA_EH(err, "load_bf_grad");

    err = load_fv();
    GOMA_EH(err, "load_fv");

    if (pde[R_MESH1] || pd->v[pg->imtrx][R_MESH1]) {
      err = load_bf_mesh_derivs();
      GOMA_EH(err, "load_bf_mesh_derivs");
    }

    /*
     * Load up physical space gradients of field variables at this
     * Gauss point.
     */
    err = load_fv_grads();
    GOMA_EH(err, "load_fv_grads");

    if (pde[R_MESH1] || pd->v[pg->imtrx][R_MESH1]) {
      err = load_fv_mesh_derivs(1);
      GOMA_EH(err, "load_fv_mesh_derivs");
    }

    computeCommonMaterialProps_gp(time_value);

    /*
     * Presumably, we have all the pieces we need.
     */
    do_LSA_mods(LSA_VOLUME);

    if (vn->evssModel == LOG_CONF || vn->evssModel == LOG_CONF_GRADV) {
      err = assemble_stress_log_conf(theta, delta_t, &pg_data);
      if (err)
        return -1;
      err = segregate_stress_update(x_update);
      GOMA_EH(err, "assemble_stress_log_conf");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_stress_log_conf");
      if (err)
        return -1;
#endif
    } else if (vn->evssModel == SQRT_CONF) {
      err = assemble_stress_sqrt_conf(theta, delta_t, &pg_data);

      GOMA_EH(err, "assemble_stress_sqrt_conf");
      if (err)
        return -1;
      err = segregate_stress_update(x_update);
      GOMA_EH(err, "assemble_stress_update");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_stress_log_conf");
      if (err)
        return -1;
#endif
    } else if (vn->evssModel == CONF) {
      err = assemble_stress_conf(theta, delta_t, &pg_data);

      GOMA_EH(err, "assemble_stress_sqrt_conf");
      if (err)
        return -1;
      err = segregate_stress_update(x_update);
      GOMA_EH(err, "assemble_stress_update");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_stress_log_conf");
      if (err)
        return -1;
#endif
    } else if (vn->evssModel == LOG_CONF_TRANSIENT || vn->evssModel == LOG_CONF_TRANSIENT_GRADV) {
      err = assemble_stress_log_conf_transient(theta, delta_t, &pg_data);

      GOMA_EH(err, "assemble_stress_log_conf");
      if (err)
        return -1;
      err = segregate_stress_update(x_update);
      GOMA_EH(err, "assemble_stress_update");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_stress_log_conf");
      if (err)
        return -1;
#endif
    }
    if (pde[R_EM_E1_REAL]) {
      err = assemble_emwave(time_value, theta, delta_t, &pg_data, R_EM_E1_REAL, EM_E1_REAL,
                            EM_E1_IMAG);
      GOMA_EH(err, "assemble_emwave");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_emwave");
      if (err)
        return -1;
#endif
    }

    if (pde[R_EM_E2_REAL]) {
      err = assemble_emwave(time_value, theta, delta_t, &pg_data, R_EM_E2_REAL, EM_E2_REAL,
                            EM_E2_IMAG);
      GOMA_EH(err, "assemble_emwave");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_emwave");
      if (err)
        return -1;
#endif
    }

    if (pde[R_EM_E3_REAL]) {
      err = assemble_emwave(time_value, theta, delta_t, &pg_data, R_EM_E3_REAL, EM_E3_REAL,
                            EM_E3_IMAG);
      GOMA_EH(err, "assemble_emwave");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_emwave");
      if (err)
        return -1;
#endif
    }

    if (pde[R_EM_E1_IMAG]) {
      err = assemble_emwave(time_value, theta, delta_t, &pg_data, R_EM_E1_IMAG, EM_E1_IMAG,
                            EM_E1_REAL);
      GOMA_EH(err, "assemble_emwave");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_emwave");
      if (err)
        return -1;
#endif
    }

    if (pde[R_EM_E2_IMAG]) {
      err = assemble_emwave(time_value, theta, delta_t, &pg_data, R_EM_E2_IMAG, EM_E2_IMAG,
                            EM_E2_REAL);
      GOMA_EH(err, "assemble_emwave");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_emwave");
      if (err)
        return -1;
#endif
    }

    if (pde[R_EM_E3_IMAG]) {
      err = assemble_emwave(time_value, theta, delta_t, &pg_data, R_EM_E3_IMAG, EM_E3_IMAG,
                            EM_E3_REAL);
      GOMA_EH(err, "assemble_emwave");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_emwave");
      if (err)
        return -1;
#endif
    }

    if (pde[R_EM_H1_REAL]) {
      err = assemble_emwave(time_value, theta, delta_t, &pg_data, R_EM_H1_REAL, EM_H1_REAL,
                            EM_H1_IMAG);
      GOMA_EH(err, "assemble_emwave");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_emwave");
      if (err)
        return -1;
#endif
    }

    if (pde[R_EM_H2_REAL]) {
      err = assemble_emwave(time_value, theta, delta_t, &pg_data, R_EM_H2_REAL, EM_H2_REAL,
                            EM_H2_IMAG);
      GOMA_EH(err, "assemble_emwave");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_emwave");
      if (err)
        return -1;
#endif
    }

    if (pde[R_EM_H3_REAL]) {
      err = assemble_emwave(time_value, theta, delta_t, &pg_data, R_EM_H3_REAL, EM_H3_REAL,
                            EM_H3_IMAG);
      GOMA_EH(err, "assemble_emwave");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_emwave");
      if (err)
        return -1;
#endif
    }

    if (pde[R_EM_H1_IMAG]) {
      err = assemble_emwave(time_value, theta, delta_t, &pg_data, R_EM_H1_IMAG, EM_H1_IMAG,
                            EM_H1_REAL);
      GOMA_EH(err, "assemble_emwave");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_emwave");
      if (err)
        return -1;
#endif
    }

    if (pde[R_EM_H2_IMAG]) {
      err = assemble_emwave(time_value, theta, delta_t, &pg_data, R_EM_H2_IMAG, EM_H2_IMAG,
                            EM_H2_REAL);
      GOMA_EH(err, "assemble_emwave");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_emwave");
      if (err)
        return -1;
#endif
    }

    if (pde[R_EM_H3_IMAG]) {
      err = assemble_emwave(time_value, theta, delta_t, &pg_data, R_EM_H3_IMAG, EM_H3_IMAG,
                            EM_H3_REAL);
      GOMA_EH(err, "assemble_emwave");
#ifdef CHECK_FINITE
      err = CHECKFINITE("assemble_emwave");
      if (err)
        return -1;
#endif
    }

    /******************************************************************************/
  }
  /* END  for (ip = 0; ip < ip_total; ip++)                               */

  /**************************************************************************/
  /*                          BLOCK 2'                                      */
  /*                   START OF SURFACE INTEGRATION LOOP                    */
  /*         LOOP OVER THE NUMBER OF SURFACE QUADRATURE POINTS              */
  /*             For Discontinuous Galerkin implementations                 */
  /**************************************************************************/

  if (discontinuous_stress) {
    int neighbor;
    int index, face;

    var = 4 * lec->max_dof * (MAX_PROB_VAR + MAX_CONC) * lec->max_dof;
    memset(lec->J_stress_neighbor, 0, sizeof(double) * var);

    for (face = 0; face < ei[pg->imtrx]->num_sides; face++) {

      index = exo->elem_elem_pntr[ielem] + face;
      neighbor = exo->elem_elem_list[index];

      id_side = face + 1;

      /*      if(neighbor != -1) */
      {
        err = assemble_surface_stress(exo, x, ams, x_update, delta_t, theta, ielem_type,
                                      ielem_type_mass, id_side, neighbor, ielem, num_local_nodes);
        GOMA_EH(err, "assemble_surface_stress");
#ifdef CHECK_FINITE
        err = CHECKFINITE("assemble_surface_stress");
        if (err)
          return -1;
#endif
      }
    }
  }

  /**************************************************************************
   *                          BLOCK 3 - Weak SURFACE Boundary Conditions
   *                   START OF SURFACE INTEGRATION LOOP
   * apply BC's that are added into the BOUNDARY term of the integrated-by-parts
   * weighted residual equations.  Do this before rotating mesh or momentum equations
   **************************************************************************/
  /******************************************************************************/
  if (first_elem_side_BC_array[ielem] != NULL) {
    elem_side_bc = first_elem_side_BC_array[ielem];
    /***************************************************************************
     *  begining of do while construct which loops over the sides of this
     *  element that have boundary conditions applied on to them.
     ***************************************************************************/
    do {

      call_int = 0;
      call_shell_grad = 0;
      call_sharp_int = FALSE;

      for (ibc = 0; (bc_input_id = (int)elem_side_bc->BC_input_id[ibc]) != -1; ibc++) {
        if (BC_Types[bc_input_id].desc->method == WEAK_INT_SURF) {
          call_int = 1;
        }
        if ((BC_Types[bc_input_id].desc->method == WEAK_SHELL_GRAD ||
             BC_Types[bc_input_id].desc->method == STRONG_SHELL_GRAD) &&
            (num_elem_friends[ielem] > 0) && (ielem_dim > 1)) {
          call_shell_grad = TRUE;
        }
        if (BC_Types[bc_input_id].desc->method == WEAK_SHARP_INT)
          call_sharp_int = TRUE;
      }

      if (call_int) {
        err = apply_integrated_bc(x, resid_vector, delta_t, theta, &pg_data, ielem, ielem_type,
                                  num_local_nodes, ielem_dim, iconnect_ptr, elem_side_bc,
                                  num_total_nodes, WEAK_INT_SURF, time_value, element_search_grid,
                                  exo);
        GOMA_EH(err, " apply_integrated_bc");
#ifdef CHECK_FINITE
        err = CHECKFINITE("apply_integrated_bc");
        if (err)
          return -1;
#endif
        if (neg_elem_volume)
          return -1;
      }

      if (call_shell_grad) {
        err = apply_shell_grad_bc(x, resid_vector, delta_t, theta, pg_data.h_elem_avg, pg_data.h,
                                  pg_data.mu_avg, pg_data.U_norm, ielem, ielem_type,
                                  num_local_nodes, ielem_dim, iconnect_ptr, elem_side_bc,
                                  num_total_nodes, WEAK_SHELL_GRAD, time_value, exo);
        GOMA_EH(err, " apply_shell_grad_bc");
#ifdef CHECK_FINITE
        err = CHECKFINITE("apply_shell_grad_bc");
        if (err)
          return -1;
#endif
        if (neg_elem_volume)
          return -1;
      }
      if (call_sharp_int) {
        err = apply_sharp_integrated_bc(x, resid_vector, time_value, delta_t, theta,
                                        pg_data.hsquared, ielem, ielem_type, num_local_nodes,
                                        ielem_dim, iconnect_ptr, elem_side_bc, WEAK_SHARP_INT, exo);
      }
      GOMA_EH(err, " apply_sharp_integrated_bc");
#ifdef CHECK_FINITE
      err = CHECKFINITE("apply_sharp_integrated_bc");
      if (err)
        return -1;
#endif
      if (neg_elem_volume)
        return -1;

    } while ((elem_side_bc = elem_side_bc->next_side_bc) != NULL);
  } /* END if (First_Elem_Side_BC_Array[ielem] != NULL) */

  /**************************************************************************
   *                          BLOCK 4 - Weak EDGE Boundary Conditions
   *                   START OF CURVE INTEGRATION LOOP
   * apply BC's that are added into the integrated-by-parts part of the BOUNDARY
   * term of the integrated-by-parts weighted residual equations.   This only applies
   * for conditions like capillary that require end-of-the-interface fluxes to
   * be added back in (i.e. if transport equations is integrated by parts twice)
   * Do this before rotating mesh or momentum equations
   **************************************************************************/
  if (First_Elem_Edge_BC_Array[pg->imtrx][ielem] != NULL) {
    elem_edge_bc = First_Elem_Edge_BC_Array[pg->imtrx][ielem];

    /******************************************************************************/
    do { /* begining of do while construct
          * which loops over the sides of this element that have boundary
          * conditions */
      /****************************************************************************/

      /*
       *  Set flags for subroutines to call for each boundary condition
       *  on this side
       */
      call_int = 0;
      for (ibc = 0; (bc_input_id = (int)elem_edge_bc->BC_input_id[ibc]) != -1; ibc++) {
        bct = BC_Types[bc_input_id].desc->method;
        if (bct == WEAK_INT_EDGE)
          call_int = 1;
      }
      /*
       * BOUNDARY CONDITIONS OF TYPE 1 - Integrated along edge -
       *        only do the weak integrated conditions here
       */
      if (call_int) {
        err = apply_integrated_curve_bc(x, resid_vector, delta_t, theta, ielem, ielem_type,
                                        num_local_nodes, ielem_dim, iconnect_ptr, elem_edge_bc,
                                        num_total_nodes, WEAK_INT_EDGE, exo);
        GOMA_EH(err, " apply_integrated_curve_bc");
#ifdef CHECK_FINITE
        err = CHECKFINITE("apply_integrated_curve_bc");
        if (err)
          return -1;
#endif
        if (neg_elem_volume)
          return -1;
      }
      /****************************************************************************/
    } while ((elem_edge_bc = elem_edge_bc->next_edge_bc) != NULL);
    /* END of do  while () construct				      */
    /******************************************************************************/
  } /* END if (First_Elem_edge_BC_Array[ielem] != NULL) 		      */

  /**************************************************************************
   *                          BLOCK 5 - Weak POINT Boundary Conditions
   *                   START OF POINT INTEGRATION LOOP  ?!?
   *  are we ever going to need a condition like this?  Probably not, because
   *  integrals over a point in 3D are zero... unless it's a singular point
   *  for now we'll assume this can't happen
   **************************************************************************/

  /******************************************************************************/
  /**************************************************************************
   *                          BLOCK 6 - Weak Shifty Boundary Conditions
   *   This section is for boundary conditions which equate the BOUNDARY terms
   *   of two different residual equations, e.g. FLUID_SOLID and SOLID_FLUID
   *   type BC's.
   *   Do this before rotating mesh or momentum equations
   **************************************************************************/
  /******************************************************************************/
  if (first_elem_side_BC_array[ielem] != NULL) {
    elem_side_bc = first_elem_side_BC_array[ielem];

    /******************************************************************************/
    do { /* begining of do while construct
          * which loops over the sides of this element that have boundary
          * conditions */
      /******************************************************************************/
      /* Determine which types of bc's apply on this side */
      /*call_special = 0; */
      for (ibc = 0; (bc_input_id = (int)elem_side_bc->BC_input_id[ibc]) != -1; ibc++) {
        bct = BC_Types[bc_input_id].desc->method;
        /*if (bct == WEAK_SHIFT) call_special = 1; */
      } /* end of loop over boundary condition number */

    } while ((elem_side_bc = elem_side_bc->next_side_bc) != NULL);
    /* END of do  while () construct				      */
    /******************************************************************************/
  } /* END if (First_Elem_Side_BC_Array[ielem] != NULL) 		      */
  /******************************************************************************/

  /**************************************************************************
   *                          BLOCK 7
   * rotate equations here before dirichlet or strong conditions
   * and do any shifting of residual and jacobian contributions
   **************************************************************************/
  /*
   * Make sure there are BC's at this node and that we have vector equations
   */
  /*if (  first_elem_side_BC_array[ielem] != NULL && */
  if ((pde[R_MOMENTUM1] || pde[R_MESH1]) && Num_ROT == 0) {
    if (ielem_dim == 2) {
      call_rotate = 0;
      /* determine if rotation is needed */
      for (i = 0; (call_rotate == 0) && (i < num_local_nodes); i++) {
        /*
         * To address a particular residual equation, map the local
         * elemental node number i into a processor node index I.
         * Then, make sure this node is owned by this processor
         */
        I = Proc_Elem_Connect[iconnect_ptr + i];
        if (I < (dpi->num_internal_nodes + dpi->num_boundary_nodes)) {
          if (in_list(I, 0, num_mom_rotate[pg->imtrx], mom_rotate_node[pg->imtrx]) != -1)
            call_rotate = 1;
          if (in_list(I, 0, num_mesh_rotate[pg->imtrx], mesh_rotate_node[pg->imtrx]) != -1)
            call_rotate = 1;
        }
      }
      /* MMH: I am assuming that apply_rotated_bc will work properly
       * in LSA b/c it doesn't "add" anything, only rearranges it,
       * and we supposedly have the correct contributions by now.
       */
      if (call_rotate) {
        if (Use_2D_Rotation_Vectors == TRUE)
          rotate_eqns_at_node_2D(iconnect_ptr, ielem_dim, num_local_nodes, ams);
        else if (first_elem_side_BC_array[ielem] != NULL) {
          err = apply_rotated_bc(resid_vector, first_elem_side_BC_array, ielem, ielem_type,
                                 num_local_nodes, ielem_dim, iconnect_ptr, num_total_nodes, exo);
          GOMA_EH(err, " apply_rotated_bc");
#ifdef CHECK_FINITE
          err = CHECKFINITE("apply_rotated_bc");
          if (err)
            return -1;
#endif
        }
      }
    }
  }

  if ((pde[R_MOMENTUM1] || pde[R_MESH1]) && Num_ROT > 0) {
    int id_mesh, id_mom; /* local temporary things */

    /* Make sure there are BC's at this node and that we have vector equations*/

    /* determine if rotation is needed */
    for (i = 0; i < num_local_nodes; i++) {
      id_mesh = ei[pg->imtrx]->ln_to_dof[MESH_DISPLACEMENT1][i];
      id_mom = ei[pg->imtrx]->ln_to_dof[VELOCITY1][i];
      /*
       *  To address a particular residual equation, map the local
       *   elemental node number i into a global index I
       *  You only have to rotate equations at nodes that you own?!?
       */
      I = Proc_Elem_Connect[iconnect_ptr + i];
      if (I < (dpi->num_internal_nodes + dpi->num_boundary_nodes)) {
        /*
         * Check for rotation of mesh equations or momentum equations
         * work on MESH first
         */
        rotate_mesh = 0;
        if (ROT_list[I] != NULL && (ROT_list[I][VECT_EQ_MESH] != -1) && pde[R_MESH1] &&
            id_mesh > -1)
          rotate_mesh = 1;

        rotate_momentum = 0;

        if (ROT_list[I] != NULL && (ROT_list[I][VECT_EQ_MOM] != -1) && pde[R_MOMENTUM1] &&
            id_mom > -1)
          rotate_momentum = 1;

        /* Now rotate the appropriate equations and sensitivities*/
        if (rotate_mesh) {
          /* Call mesh rotation routine */
          /* MMH: I am assuming that rotate_mesh_eqn will work
           * properly in LSA b/c it doesn't "add" anything, only
           * rearranges it, and we supposedly have the correct
           * contributions by now.  I am worried, though, with this
           * call b/c there is direct access to a[] in here, instead
           * of through lec->J's...
           */
          rotate_mesh_eqn(id_mesh, I, iconnect_ptr, ielem_dim, ams);
        }

        if (rotate_momentum) {
          /* Call momentum rotation routine */
          /* MMH: See comment above for rotate_mesh_eqn. */
          rotate_momentum_eqn(id_mom, I, iconnect_ptr, ielem_dim, ams);
        }
      }
    } /* end of loop over nodes */
  }   /* end of if Num_ROT > 0 */

  if (pde[R_MOMENTUM1] && goma_automatic_rotations.automatic_rotations) {
    int id_mom; /* local temporary things */
    /* determine if rotation is needed */
    for (i = 0; i < num_local_nodes; i++) {
      id_mom = ei[pg->imtrx]->ln_to_dof[VELOCITY1][i];
      /*
       *  To address a particular residual equation, map the local
       *   elemental node number i into a global index I
       *  You only have to rotate equations at nodes that you own?!?
       */
      I = Proc_Elem_Connect[iconnect_ptr + i];
      if (I < (dpi->num_internal_nodes + dpi->num_boundary_nodes)) {

        rotate_momentum = 0;

        if (goma_automatic_rotations.rotation_nodes[I].is_rotated && pde[R_MOMENTUM1]) {
          rotate_momentum = 1;
        }

        if (rotate_momentum) {
          rotate_momentum_auto(id_mom, I, ielem_dim, ams);
        }
      }
    } /* end of loop over nodes */
  }   /* end of if Num_ROT > 0 */

  /******************************************************************************/
  /*                              BLOCK 9                                       */
  /* Start of STRONG Boundary Conditions (other than Dirichlet). There are two  */
  /* types depending on whether they are integrated or pointwise collocation    */
  /******************************************************************************/
  /*
   *   initialize node list for fluid/solid boundaries
   */
  for (i = 0; i < num_local_nodes; i++) {
    local_node_list_fs[i] = -1;
  }

  /******************************************************************************/
  if (first_elem_side_BC_array[ielem] != NULL) {
    /****************************************************************************/
    /* We may want to figure out which equations are rotated and which boundary
     * conditions are applied in a strong sense (replacing equations as opposed
     * to the weak sense in which boundary conditions are added to the existing
     * equations), and do the rotation and zeroing (or subtracting local
     * contributions) now, before we get into the heavy duty boundary
     * condition stuff
     */
    elem_side_bc = first_elem_side_BC_array[ielem];

    /*****************************************************************************/
    do { /* begining of do while construct */
      /* which loops over the sides of this element that have boundary
         conditions */

      /*
       *  Set flags for subroutines to call for each boundary condition
       *  on this side
       */
      call_int = 0;
      call_col = 0;
      call_contact = 0.;
      for (ibc = 0; (bc_input_id = (int)elem_side_bc->BC_input_id[ibc]) != -1; ibc++) {
        bct = BC_Types[bc_input_id].desc->method;
        if (bct == STRONG_INT_SURF)
          call_int = 1;
        if (bct == COLLOCATE_SURF)
          call_col = 1;
        if (bct == CONTACT_SURF)
          call_contact = 1;
      }
      if (call_col) {
        err = apply_point_colloc_bc(resid_vector, delta_t, theta, ielem, ip_total, ielem_type,
                                    num_local_nodes, ielem_dim, iconnect_ptr, elem_side_bc,
                                    num_total_nodes, local_node_list_fs, time_value, exo);
        GOMA_EH(err, " apply_point_colloc_bc");
#ifdef CHECK_FINITE
        err = CHECKFINITE("apply_point_colloc_bc");
        if (err)
          return -1;
#endif
        if (neg_elem_volume)
          return -1;
      }

      /*
       * BOUNDARY CONDITIONS OF TYPE 2 - Integrated over the surface
       * - only do the strong integrated conditions here
       */
      if (call_int) {
        err = apply_integrated_bc(x, resid_vector, delta_t, theta, &pg_data, ielem, ielem_type,
                                  num_local_nodes, ielem_dim, iconnect_ptr, elem_side_bc,
                                  num_total_nodes, STRONG_INT_SURF, time_value, element_search_grid,
                                  exo);
        GOMA_EH(err, " apply_integrated_bc");
#ifdef CHECK_FINITE
        err = CHECKFINITE("apply_integrated_bc");
        if (err)
          return -1;
#endif
        /*printf("Element: %d, ID_side: %d \n", ei[pg->imtrx]->ielem, elem_side_bc->id_side );
          for(i3=0; i3< (int)  elem_side_bc->num_nodes_on_side; i3++)
          {
          id3 = (int) elem_side_bc->local_elem_node_id[i3]; I3 =  I = Proc_Elem_Connect[iconnect_ptr
          + id3]; printf("\tI: %d, R[0][%d]: %7.4g, R[1][%d]: %7.4g\n", I3, id3, lec->R[0][id3],
          id3,  lec->R[1][id3]);
          }*/

        if (neg_elem_volume)
          return -1;
      }

      /*
       * BOUNDARY CONDITIONS which don't really fit in with the above types
       */
      if (0 && (call_int || call_col)) {
        /* add call to Gibb's inequality condition is evaluated
         * within apply_special_bc, and potentially contact lines
         * are suddenly fixed here
         */
        err =
            apply_special_bc(ams, x, resid_vector, x_old, x_older, xdot, xdot_old, delta_t, theta,
                             first_elem_side_BC_array, ielem, ip_total, ielem_type, num_local_nodes,
                             ielem_dim, iconnect_ptr, elem_side_bc, num_total_nodes, SPECIAL, CA_id,
                             CA_fselem, CA_sselem, exo, time_value);
        GOMA_EH(err, " apply_special_bc");
#ifdef CHECK_FINITE
        err = CHECKFINITE("apply_special_bc");
        if (err)
          return -1;
#endif
        if (neg_elem_volume)
          return -1;
      }

      if (call_contact) {
        struct Level_Set_Data *ls_save = ls;

        /* Special contact conditions are applied here. These are boundary
           conditions that require a contact search or a mesh search, such
           as that which might happen with overlapping grids.   */

        /* Cludge her for BCs that are tied to R_PHASE0 */
        if (pfd != NULL) {
          ls = pfd->ls[0];
        } else {
          GOMA_EH(GOMA_ERROR,
                  "YOU cannot apply CONTACT_SURF BCs in mm_names.h with FILL field. R_PHASE only");
        }

        err = apply_contact_bc(x, resid_vector, delta_t, theta, pg_data.h_elem_avg, pg_data.h,
                               pg_data.mu_avg, pg_data.U_norm, first_elem_side_BC_array, ielem,
                               ielem_type, num_local_nodes, ielem_dim, iconnect_ptr, elem_side_bc,
                               num_total_nodes, CONTACT_SURF, -1, NULL, NULL, NULL, NULL,
                               time_value, exo);
        GOMA_EH(err, " apply_contact_bc");
#ifdef CHECK_FINITE
        err = CHECKFINITE("apply_contact_bc");
        if (err)
          return -1;
#endif
        if (neg_elem_volume)
          return -1;

        ls = ls_save;
      }
      /****************************************************************************/
    } while ((elem_side_bc = elem_side_bc->next_side_bc) != NULL);
    /* END of do  while () construct				      */
    /******************************************************************************/
  } /* END if (First_Elem_Side_BC_Array[ielem] != NULL) 		      */
  /******************************************************************************/

  /******************************************************************************/
  if (First_Elem_Edge_BC_Array[pg->imtrx][ielem] != NULL) {
    /******************************************************************************/
    elem_edge_bc = First_Elem_Edge_BC_Array[pg->imtrx][ielem];

    /****************************************************************************/
    do { /* begining of do while construct */
      /* which loops over the sides of this element that have boundary
         conditions */
      /****************************************************************************/

      /*
       *  Set flags for subroutines to call for each boundary
       *  condition on this side
       */
      call_int = 0;
      call_col = 0;
      for (ibc = 0; (bc_input_id = (int)elem_edge_bc->BC_input_id[ibc]) != -1; ibc++) {
        bct = BC_Types[bc_input_id].desc->method;
        if (bct == STRONG_INT_EDGE)
          call_int = 1;
        if (bct == COLLOCATE_EDGE)
          call_col = 1;
      }

      /*
       * BOUNDARY CONDITIONS OF TYPE 1 - Integrated over the surface
       * - only do the strong integrated conditions here
       */
      if (call_int) {
        err = apply_integrated_curve_bc(x, resid_vector, delta_t, theta, ielem, ielem_type,
                                        num_local_nodes, ielem_dim, iconnect_ptr, elem_edge_bc,
                                        num_total_nodes, STRONG_INT_EDGE, exo);
        GOMA_EH(err, " apply_integrated_curve_bc");
#ifdef CHECK_FINITE
        err = CHECKFINITE("apply_integrated_curve_bc");
        if (err)
          return -1;
#endif
        if (neg_elem_volume)
          return -1;
      }

      /*
       * BOUNDARY CONDITIONS OF TYPE 2 - Pointwise Collocation at the
       * nodal points
       */
      if (call_col) {
        err = apply_point_colloc_edge_bc(x, x_old, x_older, xdot, xdot_old, resid_vector, delta_t,
                                         theta, ielem, ip_total, ielem_type, num_local_nodes,
                                         ielem_dim, iconnect_ptr, elem_edge_bc, num_total_nodes,
                                         local_node_list_fs, time_value);
        GOMA_EH(err, " apply_point_colloc_bc");
#ifdef CHECK_FINITE
        err = CHECKFINITE("apply_point_colloc_bc");
        if (err)
          return -1;
#endif
        if (neg_elem_volume)
          return -1;
      }

      /****************************************************************************/
    } while ((elem_edge_bc = elem_edge_bc->next_edge_bc) != NULL);
    /* END of do  while () construct				      */
    /******************************************************************************/
  } /* END if (First_Elem_Edge_BC_Array[ielem] != NULL) 		      */

  /********************************************************************************
   *                          BLOCK 8
   * Now zero out the rows that correspond to Dirichlet Conditions
   *  i.e. conditions where the value of the solution vector is known and
   *       so the residual is set to zero and jacobian has a one on the diagonal
   *       (sometimes one is not on the diagonal, as in the case
   *        of distinguishing conditions)
   ********************************************************************************/

  err = put_dirichlet_in_matrix(x, num_total_nodes);
  GOMA_EH(err, " put_dirichlet_in_matrix");
#ifdef CHECK_FINITE
  err = CHECKFINITE("put_dirichlet_in_matrix");
  if (err)
    return -1;
#endif

  /* PRS test code for shell endpoint conditions.  4/21/2004 */
  if (pde[R_MESH1] && pde[R_SHELL_TENSION]) {
    for (i = 0; i < num_local_nodes; i++) {
      I = Proc_Elem_Connect[iconnect_ptr + i];
      node = Nodes[I];
      /*Now apply condition */
      if (node->DBSH_SLOPE_X == 1) {
        eqn = R_MESH1;
        lec->R[LEC_R_INDEX(upd->ep[pg->imtrx][eqn], i)] = 1.0 * BIG_PENALTY;
      }
      if (node->DBSH_SLOPE_Y == 1) {
        eqn = R_MESH2;
        lec->R[LEC_R_INDEX(upd->ep[pg->imtrx][eqn], i)] = 0.0 * BIG_PENALTY;
      }
    }
  }

  if ((Linear_Stability == LSA_3D_OF_2D || Linear_Stability == LSA_3D_OF_2D_SAVE) &&
      !(af->Assemble_LSA_Jacobian_Matrix) && !(af->Assemble_LSA_Mass_Matrix))
    /* Cheap, but effective!  Zero out the equation rows for
     * R_MOMENTUM3, and the columns for sensitivities to VELOCITY3,
     * then stick in a Dirichlet condition.  The end effect should be
     * that our original 2D steady-state problem has an appended
     * system equivalent to the identity problem, "I x = 0". */
    for (i = 0; i < lec->max_dof; i++) {
      eqn = upd->ep[pg->imtrx][R_MOMENTUM3];
      var = upd->vp[pg->imtrx][VELOCITY3];
      zero_lec_row(lec->J, eqn, i);
      zero_lec_column(lec->J, var, i);
      lec->J[LEC_J_INDEX(eqn, var, i, i)] = 1.0;
      lec->R[LEC_R_INDEX(eqn, i)] = 0.0;
    }
  /* If we're doing normal mode LSA, and we're on wavenumber 0, then
   * the jacobian contributions to the w-momentum equation are
   * (rightly) all zero.  This causes a zero row, which we don't
   * really want.  We want to emulate the Dirichlet = 0 "boundary
   * condition" for all internal nodes, so we put a 1 on the diagonal
   * for the same reason(s) as above.  We don't want this to be
   * subtracted out, so we only apply it on one pass. */
  if ((Linear_Stability == LSA_3D_OF_2D || Linear_Stability == LSA_3D_OF_2D_SAVE) &&
      af->Assemble_LSA_Jacobian_Matrix && LSA_3D_of_2D_wave_number == 0 && LSA_3D_of_2D_pass == 1)
    for (i = 0; i < lec->max_dof; i++)
      lec->J[LEC_J_INDEX(R_MOMENTUM3, VELOCITY3, i, i)] = 1.0;

  /*
   * Load local element stiffness matrix (lec) into global matrix, depending
   * on solver used - front vs. everything else, and matrix storage format
   * MSR vs VBR.
   */

  load_lec(exo, ielem, ams, x, resid_vector, estifm);

  /*  if( pfd != NULL && pfd->Use_Constraint == TRUE )
      {
      load_pf_constraint(  pf_constraint, d_pf_lm, d_lm_pf );
      }
  */

  /******************************************************************************/
  /*                              BLOCK 13                                      */
  /*                            CLEANUP BLOCK                                   */
  /******************************************************************************/

  /*
   *  Good spot to print out the element stiffness matrix and residual or to
   *  free any memory that was allocated on the element level
   */
  mm_fill_end = ut();
  mm_fill_total += mm_fill_end - mm_fill_start;

  /*
   *  Free up Subgrid integration arrays
   */

  if (Subgrid_Int.active) {
    safe_free((void *)Subgrid_Int.s);
    Subgrid_Int.s = NULL;
    safe_free((void *)Subgrid_Int.wt);
    Subgrid_Int.wt = NULL;
    free_search_grid(&element_search_grid);
  }

  if (Debug_Flag > 1) {
    P0PRINTF("%s: ends\n", yo);
    MMH_ip = -1;
  }
  if ((Linear_Solver != FRONT && ielem == exo->eb_ptr[exo->num_elem_blocks] - 1) ||
      (Linear_Solver == FRONT && ielem == exo->elem_order_map[exo->num_elem_blocks] - 1)) {
    if (zeroCA == 0) {
      int count = 0, Num_CAs_done = 0;
      for (j = 0; j < MAX_CA; j++) {
        if (CA_id[j] == -2)
          Num_CAs_done++;
        if (CA_proc[j] == ProcID)
          count++;
      }

      if (count != Num_CAs_done) {
        GOMA_WH(-1, "\nNot all contact angle conditions were applied!\n");
        for (j = 0; j < count; j++) {
          fprintf(stderr, "CA:%d ID:%d fselem:%d sselem:%d Proc:%d\n", j, CA_id[j], CA_fselem[j],
                  CA_sselem[j], CA_proc[j]);
        }
        fprintf(stderr, "Count=%d  Done=%d\n", count, Num_CAs_done);
      }
    }
  }

  return 0;
} /*   END OF matrix_fill_stress                                                     */

static void load_lec(Exo_DB *exo, /* ptr to EXODUS II finite element mesh db */
                     int ielem,   /* Element number we are working on */
                     struct GomaLinearSolverData *ams,
                     double x[],            /* Solution vector */
                     double resid_vector[], /* Residual vector */
                     dbl *estifm)           /* element stiffness Matrix for frontal solver*/

/**************************************************************************
 *
 * load_lec was designed to take the local stiffness matrix, lec, and
 * load it  into the global matrix in three different ways depending
 * on user specification.  These ways are 1) Standard GOMA MSR format
 * 2) Frontal solver format where the  solver itself compiles the global
 * matrix 3)  VBR format .
 *
 * Written by: RRR 6/99
 * Revised:
 *************************************************************************/
{
  int e, v, i, j, pe, pv;
  int dofs;
  int I, J, K;
  int gnn, ie, ke, kv, nvdof;
  int je, ja, ledof;
  int je_new;
  struct Element_Indices *ei_ptr;
#ifdef DEBUG_LEC
  char lec_name[256], ler_name[256];
  FILE *llll, *rrrr;
  static int lec_it = -1;
  int Print_Zeroes = TRUE;
  lec_it++;

  sprintf(lec_name, "lec_dump_%d_%d.txt", DPI_ptr->elem_index_global[ei[pg->imtrx]->ielem], ProcID);
  sprintf(ler_name, "ler_dump_%d_%d.txt", DPI_ptr->elem_index_global[ei[pg->imtrx]->ielem], ProcID);
  llll = fopen(lec_name, "a");
  rrrr = fopen(ler_name, "a");
  fprintf(rrrr, "------------------------------------------------------\n");
  fprintf(rrrr, "local element = %d, Proc = %d\n", ei[pg->imtrx]->ielem, ProcID);
  fprintf(rrrr, "lec_it = %d\n", lec_it);
  fprintf(rrrr, "global element = %d\n", DPI_ptr->elem_index_global[ei[pg->imtrx]->ielem]);
  fprintf(rrrr, "\nGlobal_NN Proc_NN  Equation    idof    Proc_SolnNum     ResidValue\n");
#endif

  if (strcmp(Matrix_Format, "epetra") == 0) {
    EpetraLoadLec(ielem, ams, resid_vector);
  }
#ifdef GOMA_ENABLE_PETSC
#if PETSC_USE_COMPLEX
  else if (strcmp(Matrix_Format, "petsc_complex") == 0) {
    petsc_load_lec_complex(ielem, ams, resid_vector);
  }
#else
  else if (strcmp(Matrix_Format, "petsc") == 0) {
    petsc_load_lec(ielem, ams, resid_vector);
  }
#endif
#endif
  else {
    /* load up matrix in MSR format */
    if (strcmp(Matrix_Format, "msr") == 0) {
      double *a = ams->val;
      int *ija = ams->bindx;

      for (e = V_FIRST; e < V_LAST; e++) {
        pe = upd->ep[pg->imtrx][e];
        if (pe != -1) {
          if (e == R_MASS) {
            for (ke = 0; ke < upd->Max_Num_Species_Eqn; ke++) {
              pe = MAX_PROB_VAR + ke;
              dofs = ei[pg->imtrx]->dof[e];
              for (i = 0; i < dofs; i++) {
                /*
                 * Check to see whether this dof is on a
                 * node owned by this processor. If it isn't
                 * there is no need to fill in the processor
                 * residual and Jacobian row for this variable.
                 */
                ledof = ei[pg->imtrx]->lvdof_to_ledof[e][i];
                if (ei[pg->imtrx]->owned_ledof[ledof]) {
                  gnn = ei[pg->imtrx]->gnn_list[e][i];
                  nvdof = ei[pg->imtrx]->Baby_Dolphin[e][i];
                  je_new = ei[pg->imtrx]->ieqn_ledof[ledof] + ke;
                  ie = Index_Solution(gnn, e, ke, nvdof, ei[pg->imtrx]->matID_ledof[ledof],
                                      pg->imtrx);

                  resid_vector[ie] += lec->R[LEC_R_INDEX(MAX_PROB_VAR + ke, i)];
#ifdef DEBUG_LEC
                  {
                    if (fabs(lec->R[LEC_R_INDEX(MAX_PROB_VAR + ke, i)]) > DBL_SMALL ||
                        Print_Zeroes) {
                      fprintf(rrrr, "%7d %7d MF%-3d %9d -  %12d - %10.3g\n",
                              DPI_ptr->node_index_global[gnn], gnn, ke, i, ie,
                              lec->R[LEC_R_INDEX(MAX_PROB_VAR + ke, i)]);
                    }
                  }
#endif

                  if (af->Assemble_Jacobian) {
                    for (v = V_FIRST; v < V_LAST; v++) {
                      pv = upd->vp[pg->imtrx][v];
                      if (pv != -1 && (Inter_Mask[pg->imtrx][e][v])) {
                        ei_ptr = ei[pg->imtrx];
                        if (ei[pg->imtrx]->owningElementForColVar[v] != ielem) {
                          if (ei[pg->imtrx]->owningElementForColVar[v] != -1) {
                            ei_ptr = ei[pg->imtrx]->owningElement_ei_ptr[v];
                            if (ei_ptr == 0) {
                              printf("ei_ptr == 0\n");
                              exit(-1);
                            }
                          }
                        }
                        if (v == MASS_FRACTION) {
                          for (kv = 0; kv < upd->Max_Num_Species_Eqn; kv++) {
                            pv = MAX_PROB_VAR + kv;
                            for (j = 0; j < ei_ptr->dof[v]; j++) {
                              ledof = ei_ptr->lvdof_to_ledof[v][j];
                              je_new = ei_ptr->ieqn_ledof[ledof] + kv;
                              je = Index_Solution(ei_ptr->gnn_list[v][j], v, kv,
                                                  ei_ptr->Baby_Dolphin[v][j],
                                                  ei_ptr->matID_ledof[ledof], pg->imtrx);
                              if (je != je_new) {
                                /*
                                 * HKM -> another special case. Until we delineate
                                 *        the subspecies as individual local
                                 *        dofs, je_new will be wrong here. je
                                 *        is correct.
                                 */
                                if (Nodes[ei_ptr->gnn_list[v][j]]->Mat_List.Length < 2) {
                                }
                              }
                              GOMA_EH(je, "Bad var index.");
                              ja = (ie == je) ? ie : in_list(je, ija[ie], ija[ie + 1], ija);
                              GOMA_EH(ja, "Could not find vbl in sparse matrix.");
                              a[ja] += lec->J[LEC_J_INDEX(pe, pv, i, j)];

#ifdef DEBUG_LEC
                              {
                                if (fabs(lec->J[LEC_J_INDEX(pe, pv, i, j)]) > DBL_SMALL ||
                                    Print_Zeroes) {
                                  fprintf(llll, "%9d %9d %9d %9d -  %12d - %10.3g\n", pe, pv, i, j,
                                          ja, lec->J[LEC_J_INDEX(pe, pv, i, j)]);
                                }
                              }
#endif
                            }
                          }
                        } else {
                          pv = upd->vp[pg->imtrx][v];
                          kv = 0;
                          for (j = 0; j < ei_ptr->dof[v]; j++) {
                            ledof = ei_ptr->lvdof_to_ledof[v][j];
                            je_new = ei_ptr->ieqn_ledof[ledof];
                            je = Index_Solution(ei_ptr->gnn_list[v][j], v, kv,
                                                ei_ptr->Baby_Dolphin[v][j],
                                                ei_ptr->matID_ledof[ledof], pg->imtrx);
                            if (je != je_new) {
                              fprintf(stderr, "Oh fiddlesticks: je = %d, je_new = %d\n", je,
                                      je_new);
                              GOMA_EH(GOMA_ERROR, "LEC Indexing error");
                            }
                            GOMA_EH(je, "Bad var index.");
                            ja = (ie == je) ? ie : in_list(je, ija[ie], ija[ie + 1], ija);
                            GOMA_EH(ja, "Could not find vbl in sparse matrix.");
                            a[ja] += lec->J[LEC_J_INDEX(pe, pv, i, j)];
#ifdef DEBUG_LEC
                            {
                              if (fabs(lec->J[LEC_J_INDEX(pe, pv, i, j)]) > DBL_SMALL ||
                                  Print_Zeroes) {
                                fprintf(llll, "%9d %9d %9d %9d -  %12d - %10.3g\n", pe, pv, i, j,
                                        ja, lec->J[LEC_J_INDEX(pe, pv, i, j)]);
                              }
                            }
#endif
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          } else {
            pe = upd->ep[pg->imtrx][e];
            dofs = ei[pg->imtrx]->dof[e];
            for (i = 0; i < dofs; i++) {
              /*
               * Check to see whether this dof is on a
               * node owned by this processor. If it isn't
               * there is no need to fill in the processor
               * residual and Jacobian row for this variable.
               */
              ledof = ei[pg->imtrx]->lvdof_to_ledof[e][i];
              if (ei[pg->imtrx]->owned_ledof[ledof]) {
                ie = ei[pg->imtrx]->gun_list[e][i];
                resid_vector[ie] += lec->R[LEC_R_INDEX(pe, i)];
#ifdef DEBUG_LEC
                {
                  if (fabs(lec->R[LEC_R_INDEX(pe, i)]) > DBL_SMALL || Print_Zeroes) {

                    fprintf(rrrr, "%9d %9d -  %12d - %10.3g\n", pe, i, ie,
                            lec->R[LEC_R_INDEX(pe, i)]);
                  }
                }
#endif

                if (af->Assemble_Jacobian) {
                  for (v = V_FIRST; v < V_LAST; v++) {
                    pv = upd->vp[pg->imtrx][v];
                    if (pv != -1 && (Inter_Mask[pg->imtrx][e][v])) {
                      if (v == MASS_FRACTION) {

                        ei_ptr = ei[pg->imtrx];
                        if (ei[pg->imtrx]->owningElementForColVar[v] != ielem) {
                          if (ei[pg->imtrx]->owningElementForColVar[v] != -1) {
                            ei_ptr = ei[pg->imtrx]->owningElement_ei_ptr[v];
                            if (ei_ptr == 0) {
                              GOMA_EH(GOMA_ERROR, "ei_ptr == 0\n");
                              exit(-1);
                            }
                          }
                        }
                        for (kv = 0; kv < upd->Max_Num_Species_Eqn; kv++) {
                          pv = MAX_PROB_VAR + kv;
                          for (j = 0; j < ei_ptr->dof[v]; j++) {
                            ledof = ei_ptr->lvdof_to_ledof[v][j];
                            je_new = ei_ptr->ieqn_ledof[ledof] + kv;
                            je = Index_Solution(ei_ptr->gnn_list[v][j], v, kv,
                                                ei_ptr->Baby_Dolphin[v][j],
                                                ei_ptr->matID_ledof[ledof], pg->imtrx);
                            if (je != je_new) {
                              /*
                               * HKM -> another special case. Until we delineate
                               *        the subspecies as individual local
                               *        dofs, je_new will be wrong here. je
                               *        is correct.
                               */
                              if (Nodes[ei[pg->imtrx]->gnn_list[v][j]]->Mat_List.Length < 2) {
                              }
                            }
                            GOMA_EH(je, "Bad var index.");
                            ja = (ie == je) ? ie : in_list(je, ija[ie], ija[ie + 1], ija);
                            GOMA_EH(ja, "Could not find vbl in sparse matrix.");
                            a[ja] += lec->J[LEC_J_INDEX(pe, pv, i, j)];

#ifdef DEBUG_LEC
                            {
                              if (fabs(lec->J[LEC_J_INDEX(pe, pv, i, j)]) > DBL_SMALL ||
                                  Print_Zeroes) {
                                fprintf(llll, "%9d %9d %9d %9d -  %12d - %10.3g\n", pe, pv, i, j,
                                        ja, lec->J[LEC_J_INDEX(pe, pv, i, j)]);
                              }
                            }
#endif
                          }
                        }
                      } else {
                        kv = 0;
                        // NEW IMPLEMENTATION
                        ei_ptr = ei[pg->imtrx];
                        if (ei[pg->imtrx]->owningElementForColVar[v] != ielem) {
                          if (ei[pg->imtrx]->owningElementForColVar[v] != -1) {
                            ei_ptr = ei[pg->imtrx]->owningElement_ei_ptr[v];
                            if (ei_ptr == 0) {
                              GOMA_EH(GOMA_ERROR, "ei slave pointer is null");
                              exit(-1);
                            }
                          }
                        }
                        for (j = 0; j < ei_ptr->dof[v]; j++) {
                          ledof = ei_ptr->lvdof_to_ledof[v][j];
                          je_new = ei_ptr->ieqn_ledof[ledof];
                          je = Index_Solution(ei_ptr->gnn_list[v][j], v, kv,
                                              ei_ptr->Baby_Dolphin[v][j],
                                              ei_ptr->matID_ledof[ledof], pg->imtrx);
                          if (je != je_new) {
                            fprintf(stderr, "Oh fiddlesticks: je = %d, je_new = %d\n", je, je_new);
                            GOMA_EH(GOMA_ERROR, "LEC Indexing error");
                          }
                          GOMA_EH(je, "Bad var index.");
                          ja = (ie == je) ? ie : in_list(je, ija[ie], ija[ie + 1], ija);
                          GOMA_EH(ja, "Could not find vbl in sparse matrix.");
                          a[ja] += lec->J[LEC_J_INDEX(pe, pv, i, j)];
#ifdef DEBUG_LEC
                          {
                            if (fabs(lec->J[LEC_J_INDEX(pe, pv, i, j)]) > DBL_SMALL ||
                                Print_Zeroes) {
                              fprintf(llll, "%9d %9d %9d %9d -  %12d - %10.3g\n", pe, pv, i, j, ja,
                                      lec->J[LEC_J_INDEX(pe, pv, i, j)]);
                            }
                          }
#endif
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    } /* Matrix_Format == "msr" */
    /* load up matrix in VBR format */
    else if (strcmp(Matrix_Format, "vbr") == 0) {
      double *a = ams->val;
      int *bindx = ams->bindx;
      int *bpntr = ams->bpntr;
      int *rpntr = ams->rpntr;
      int *indx = ams->indx;

      int row_dofs;
      int idofs, jdofs;
      int blk_row, blk_col;
      int lec_row, lec_col;
      double *a_ptr;

      for (i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
        /* I is the global row block index, also global node number */
        I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
        /* if I is an external node to this processor, skip all this stuff */
        if (I < ams->npn) {
          row_dofs = rpntr[I + 1] - rpntr[I];
          for (e = V_FIRST; e < V_LAST; e++) {
            pe = upd->ep[pg->imtrx][e];
            if (pe != -1) {
              if (e == R_MASS) {
                for (ke = 0; ke < upd->Max_Num_Species_Eqn; ke++) {
                  pe = MAX_PROB_VAR + ke;
                  blk_row = Local_Offset[pg->imtrx][I][e] + ke * Dolphin[pg->imtrx][I][e];
                  lec_row = ei[pg->imtrx]->ln_to_first_dof[e][i];
                  for (idofs = 0; idofs < Dolphin[pg->imtrx][I][e]; idofs++) {
                    *(resid_vector + rpntr[I] + blk_row + idofs) +=
                        lec->R[LEC_R_INDEX(pe, lec_row + idofs)];

                    if (af->Assemble_Jacobian) {
                      for (v = V_FIRST; v < V_LAST; v++) {
                        pv = upd->vp[pg->imtrx][v];
                        if ((pv != -1) && (Inter_Mask[pg->imtrx][e][v])) {
                          // NEW
                          ei_ptr = ei[pg->imtrx];
                          if (ei[pg->imtrx]->owningElementForColVar[v] != ielem) {
                            if (ei[pg->imtrx]->owningElementForColVar[v] != -1) {
                              ei_ptr = ei[pg->imtrx]->owningElement_ei_ptr[v];
                              if (ei_ptr == 0) {
                                GOMA_EH(GOMA_ERROR, "ei_ptr == 0\n");
                                exit(-1);
                              }
                            }
                          }
                          for (j = 0; j < ei_ptr->num_local_nodes; j++) {
                            J = Proc_Elem_Connect[ei_ptr->iconnect_ptr +
                                                  j]; /* J is the global column block index */
                            K = in_list(J, bpntr[I], bpntr[I + 1],
                                        bindx); /* K is the block index */
                            GOMA_EH(K, " Can't locate column index in bindx ");
                            a_ptr =
                                a + indx[K]; /* a_ptr points to first entry in val of Kth block */

                            if (v == MASS_FRACTION) {
                              for (kv = 0; kv < upd->Max_Num_Species_Eqn; kv++) {
                                pv = MAX_PROB_VAR + kv;
                                blk_col =
                                    Local_Offset[pg->imtrx][J][v] + kv * Dolphin[pg->imtrx][J][v];
                                lec_col = ei_ptr->ln_to_first_dof[v][j];
                                for (jdofs = 0; jdofs < Dolphin[pg->imtrx][J][v]; jdofs++) {
                                  *(a_ptr + row_dofs * (blk_col + jdofs) + blk_row + idofs) +=
                                      lec->J[LEC_J_INDEX(pe, pv, lec_row + idofs, lec_col + jdofs)];
                                }
                              }
                            } else {
                              blk_col = Local_Offset[pg->imtrx][J][v];
                              lec_col = ei_ptr->ln_to_first_dof[v][j];
                              for (jdofs = 0; jdofs < Dolphin[pg->imtrx][J][v]; jdofs++) {
                                *(a_ptr + row_dofs * (blk_col + jdofs) + blk_row + idofs) +=
                                    lec->J[LEC_J_INDEX(pe, pv, lec_row + idofs, lec_col + jdofs)];
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              } else {
                blk_row = Local_Offset[pg->imtrx][I][e];
                lec_row = ei[pg->imtrx]->ln_to_first_dof[e][i];
                for (idofs = 0; idofs < Dolphin[pg->imtrx][I][e]; idofs++) {
                  *(resid_vector + rpntr[I] + blk_row + idofs) +=
                      lec->R[LEC_R_INDEX(pe, lec_row + idofs)];
                  if (af->Assemble_Jacobian) {
                    for (v = V_FIRST; v < V_LAST; v++) {
                      pv = upd->vp[pg->imtrx][v];
                      if ((pv != -1) && (Inter_Mask[pg->imtrx][e][v])) {
                        // NEW
                        ei_ptr = ei[pg->imtrx];
                        if (ei[pg->imtrx]->owningElementForColVar[v] != ielem) {
                          if (ei[pg->imtrx]->owningElementForColVar[v] != -1) {
                            ei_ptr = ei[pg->imtrx]->owningElement_ei_ptr[v];
                            if (ei_ptr == 0) {
                              printf("ei_ptr == 0\n");
                              exit(-1);
                            }
                          }
                        }
                        for (j = 0; j < ei_ptr->num_local_nodes; j++) {
                          J = Proc_Elem_Connect[ei_ptr->iconnect_ptr +
                                                j]; /* J is the global column block index */
                          K = in_list(J, bpntr[I], bpntr[I + 1], bindx); /* K is the block index */
                          GOMA_EH(K, " Can't locate column index in bindx ");
                          a_ptr = a + indx[K]; /* a_ptr points to first entry in val of Kth block */

                          if (v == MASS_FRACTION) {
                            for (kv = 0; kv < upd->Max_Num_Species_Eqn; kv++) {
                              pv = MAX_PROB_VAR + kv;
                              blk_col =
                                  Local_Offset[pg->imtrx][J][v] + kv * Dolphin[pg->imtrx][J][v];
                              lec_col = ei_ptr->ln_to_first_dof[v][j];
                              for (jdofs = 0; jdofs < Dolphin[pg->imtrx][J][v]; jdofs++) {
                                *(a_ptr + row_dofs * (blk_col + jdofs) + blk_row + idofs) +=
                                    lec->J[LEC_J_INDEX(pe, pv, lec_row + idofs, lec_col + jdofs)];
                              }
                            }
                          } else {
                            /*  All variables of the same type at the same node are
                             *  contiguous within the solution vector.
                             *  Find the local offset first
                             */
                            blk_col = Local_Offset[pg->imtrx][J][v];
                            /*  For the current node, j, find the index for the
                             *  first dof of type v in the LESM.
                             */
                            lec_col = ei_ptr->ln_to_first_dof[v][j];
                            /*  Dolphin contains the number of dofs of type v at the
                             *  global node J
                             */
                            for (jdofs = 0; jdofs < Dolphin[pg->imtrx][J][v]; jdofs++) {
                              *(a_ptr + row_dofs * (blk_col + jdofs) + blk_row + idofs) +=
                                  lec->J[LEC_J_INDEX(pe, pv, lec_row + idofs, lec_col + jdofs)];
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
#ifdef DEBUG_LEC
  fclose(llll);
  fclose(rrrr);
#endif
}
/****************************************************************************/

static void zero_lec(void)

/**************************************************************************
 *
 * zero_lec()
 *
 *  This routine zeroes the local element stiffness vector and Jacobian.
 **************************************************************************/
{
  memset(lec->R, 0, MAX_LOCAL_VAR_DESC * lec->max_dof * sizeof(dbl));
  if (af->Assemble_Jacobian) {
    memset(lec->J, 0,
           MAX_LOCAL_VAR_DESC * MAX_LOCAL_VAR_DESC * lec->max_dof * lec->max_dof * sizeof(dbl));
    memset(lec->J_stress_neighbor, 0,
           4 * lec->max_dof * MAX_LOCAL_VAR_DESC * lec->max_dof * sizeof(dbl));
  }
}
/****************************************************************************/

int checkfinite(const char *file, const int line, const char *message) {
  /* look through lec->R and lec-J for NAN's or anything else inappropriate */
  int eqn, var, peqn, pvar;
  int i;
  int j = 0;
  int all_finite = TRUE;
  struct Element_Indices *ei_ptr;
  int ielem = ei[pg->imtrx]->ielem;
  /*   fprintf(stdout,"Hi. We are in checkfinite.");  */
  /*   if ( !isfinite(x/y)) */
  /*     { */
  /*       fprintf(stderr,"It's working for 0/0.\n");  */
  /*     } */
  /*   if ( !isfinite(1/y)) */
  /*     { */
  /*       fprintf(stderr,"It's working for 1/0.\n");  */
  /*     } */
  for (eqn = V_FIRST; eqn < V_LAST; eqn++) {
    peqn = upd->ep[pg->imtrx][eqn];
    if (peqn != -1) {
      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        j = finite(lec->R[LEC_R_INDEX(peqn, i)]);
        if (!j) {
          fprintf(stderr, "lec->R[%s][edof=%d] = %g\n", EQ_Name[eqn].name1, i,
                  lec->R[LEC_R_INDEX(peqn, i)]);
          all_finite = FALSE;
          goto non_finite_entry;
        }
        for (var = V_FIRST; var < V_LAST; var++) {
          pvar = upd->vp[pg->imtrx][var];
          if (pvar != -1) {
            ei_ptr = ei[pg->imtrx];
            if (ei[pg->imtrx]->owningElementForColVar[var] != ielem) {
              if (ei[pg->imtrx]->owningElementForColVar[var] != -1) {
                ei_ptr = ei[pg->imtrx]->owningElement_ei_ptr[var];
              }
            }
            for (j = 0; j < ei_ptr->dof[var]; j++) {
              if (!isfinite(lec->J[LEC_J_INDEX(peqn, pvar, i, j)])) {
                fprintf(stderr, "lec->J[%s][%s][edof=%d][vdof=%d] = %g\n", EQ_Name[eqn].name1,
                        Var_Name[var].name1, i, j, lec->J[LEC_J_INDEX(peqn, pvar, i, j)]);
                all_finite = FALSE;
                goto non_finite_entry;
              }
            }
          }
        }
      }
    }
  }

non_finite_entry:
  if (all_finite != TRUE) {
#ifndef PARALLEL
    fprintf(stderr, "Non-finite entry in lec found at  %s:%d: %s\n", file, line, message);
#else
    fprintf(stderr, "P_%d: Non-finite entry in lec found at %s:%d: %s\n", ProcID, file, line,
            message);
#endif
    return -1;
  }

  return 0;
}

int load_pf_constraint(double pf_constraint, double d_pf_lm[][MDE], double d_lm_pf[][MDE]) {
  int i, j, ledof;
  int ie;
  int var;

  pfd->Constraint_Integral += pf_constraint;

  for (i = 0; i < pfd->num_phase_funcs; i++) {
    var = PHASE1 + i;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        ledof = ei[pg->imtrx]->lvdof_to_ledof[var][j];
        if (ei[pg->imtrx]->owned_ledof[ledof]) {
          ie = ei[pg->imtrx]->gun_list[var][j];
          pfd->jac_info->d_pf_lm[ie] += d_pf_lm[i][j];
          pfd->jac_info->d_lm_pf[ie] += d_lm_pf[i][j];
        }
      }
    }
  }
  return (TRUE);
}
/****************************************************************************/
