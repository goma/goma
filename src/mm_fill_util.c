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
 *$Id: mm_fill_util.c,v 5.24 2010-02-26 21:40:14 prschun Exp $
 */

#ifdef USE_RCSID
static char rcsid[] =
    "$Id: mm_fill_util.c,v 5.24 2010-02-26 21:40:14 prschun Exp $";
#endif

/*
 * Added load_elem_dofptr stuff to setup elemental level indeces pointing
 * to each dof for ea variable, also pre loaded pointers
 */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "el_elm.h"
#include "el_geom.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_mp_const.h"
#include "rf_allo.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_io.h"
#include "rf_masks.h"
#include "rf_mp.h"
#include "rf_solver.h"
#include "rf_vars_const.h"
#include "std.h"
#include "mm_eh.h"
#include "mm_mp.h"
#include "mm_mp_structs.h"
#include "bc_contact.h"
#include "dp_utils.h"
#include "dpi.h"
#include "el_elm_info.h"
#include "exo_struct.h"
#include "mm_fill_ls.h"
#include "mm_fill_util.h"
#include "mm_post_def.h"
#include "rd_mesh.h"
#include "rf_node_const.h"
#include "rf_shape.h"
#include "sl_util_structs.h"

#define GOMA_MM_FILL_UTIL_C

/*********** R O U T I N E S   I N   T H I S   F I L E *************************
 *
 *       NAME			TYPE			CALLED_BY
 *    ------------             ---------               --------------
 *
 *  beer_belly()
 *  calc_surf_det ()             double		mm_fill.c:  matrix_fill
 *  calc_surf_normal ()          void		mm_fill.c:  matrix_fill
 *  calc_surf_tangent ()         void		mm_fill.c:  matrix_fill
 *  fill_surf_shape ()		void		mm_fill.c:  matrix_fill
 *  alloc_MSR_sparse_arrays ()   void		rf_solve.c: solve_problem
 *  alloc_VBR_sparse_arrays ()   void		rf_solve.c: solve_problem
 *  find_problem_graph_fill ()   static int	alloc_MSR_sparse_arrays
 *  find_MSR_problem_graph ()    static int	alloc_MSR_sparse_arrays
 *  find_VBR_problem graph ()    static int      alloc_VBR_sparse_arrays
 *  fill_variable_vector()       static int      find_MSR_problem_graph
 *  set_diag_to_one ()           void		mm_fill.c:  matrix_fill
 *  set_diag_to_large ()         void		mm_fill.c:  matrix_fill
 *
 *******************************************************************************/

/*
 * Prototypes declarations of static functions in defined in this file.
 */

static int
find_MSR_problem_graph(int *[], /* ija - column pointer array                */
                       int,     /* itotal_nodes - number of nodes this
                                 * processor is resposible for               */
                       Exo_DB *); /* exo - ptr to FE EXODUS II database */

static int
find_problem_graph_fill(int *[], /* ija - column pointer array                */
                        int,     /* itotal_nodes - number of nodes this
                                  * processor is resposible for               */
                        int,     /* max_neigh_elem - maximum number of elements
                                  * containing a given node                   */
                        int[],   /* node_to_fill                              */
                        Exo_DB *); /* exo - ptr to FE EXODUS II database */

static int find_VBR_problem_graph(int *[], int *[], int *[], int *[], int *[],
                                  int, int, Exo_DB *);

static int var_if_interp_type_enabled(PROBLEM_DESCRIPTION_STRUCT *pd_ptr,
                               int interp_type); 
static dbl yzbeta1(dbl scale, int dim, dbl Y, dbl Z, dbl d_Z[MDE], dbl beta,
                       dbl u, dbl d_u[MDE], dbl grad_u[DIM],
                       dbl d_grad_u[MDE][DIM], dbl h_elem, int interp_eqn,
                       dbl deriv[MDE]); 

static dbl yzbeta2(dbl scale, dbl Y, dbl Z, dbl d_Z[MDE], dbl deriv[MDE], dbl h_elem, int interp_eqn);

static const dbl DIFFUSION_EPSILON = 1e-8;

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* beer_belly() -- find J, inverse, and mesh derivatives for all basis functions
 *
 * Created:	Tue Mar 15 12:29:51 MST 1994 pasacki@sandia.gov
 *
 *
 *
 * Revised:	Wed Mar 23 07:28:04 MST 1994 pasacki@sandia.gov
 *
 * 		-- Added provision for no moving mesh cases, with
 *		   choice for the shape function based on first active
 *		   variable encountered.
 *
 *		   Or should it be based on the maximum number of
 *		   nodes per element in the mesh?
 *		   No, because the mesh
 *		   may have 9 nodes/quad, but no variables interpolated
 *		   with anything greater than 4 nodes/quad...
 *
 *		   I'll just do it based on the first active eqn/vbl, but
 *		   the results are your responsibility if you choose
 *		   a nonLagrangian interpolation polynomial for your element.
 *
 * Revised:	9/19/94 by RRR for transient old mesh stuff
 *
 *   Notes for shells: For shells, dim is different than pdim. Also, the
 *   number of mesh unknowns is equal to the number of shell mesh equations,
 *   not the number of parent-element mesh equations.
 *   It would appear that MapBF is set wrong for shell equations. However,
 *   it doesn't make a difference.
 */
int beer_belly(void) {
  int status = 0, i, j, k, n, t, dim, pdim, mdof, index, node, si;
  int DeformingMesh, ShapeVar;
  struct Basis_Functions *MapBf;
  size_t v_length;
  dbl f, g, sum;
  int imtrx = upd->matrix_index[pd->ShapeVar];

  static int is_initialized = FALSE;
  static int elem_blk_id_save = -123;

  dim = ei[imtrx]->ielem_dim;
  pdim = pd->Num_Dim;
  int elem_type = ei[imtrx]->ielem_type;
  int elem_shape = type2shape(elem_type);

  ShapeVar = pd->ShapeVar;

  /* If this is a shell element, it may be a deforming mesh
   * even if there are no mesh equations on the shell block.
   * The ei[imtrx]->deforming_mesh flag is TRUE for shell elements when
   * there are mesh equations active on either the shell block or
   * any neighboring bulk block.
   */

  if (pd->gv[MESH_DISPLACEMENT1]) {
    DeformingMesh = ei[upd->matrix_index[MESH_DISPLACEMENT1]]->deforming_mesh;
  } else {
    DeformingMesh = ei[imtrx]->deforming_mesh;
  }

  if ((si = in_list(pd->IntegrationMap, 0, Num_Interpolations,
                    Unique_Interpolations)) == -1) {
    EH(-1, "Seems to be a problem finding the IntegrationMap interpolation.");
  }
  MapBf = bfd[si];

  mdof = ei[imtrx]->dof[ShapeVar];

  /*
   * For every type "t" of unique basis function used in this problem,
   * initialize appropriate arrays...
   */

  if (ei[imtrx]->elem_blk_id != elem_blk_id_save) {
    is_initialized = FALSE;
  }

  for (t = 0; t < Num_Basis_Functions; t++) {
    bfd[t]->detJ = 0.0;
    v_length = DIM * DIM * sizeof(double);
    memset(&(bfd[t]->J[0][0]), 0, v_length);
    memset(&(bfd[t]->B[0][0]), 0, v_length);

    if (DeformingMesh && !is_initialized) {
      v_length = DIM * MDE * sizeof(double);
      memset(&(bfd[t]->d_det_J_dm[0][0]), 0, v_length);

      v_length = DIM * DIM * DIM * MDE * sizeof(double);
      memset(&(bfd[t]->dJ[0][0][0][0]), 0, v_length);
      memset(&(bfd[t]->dB[0][0][0][0]), 0, v_length);
    }
  }

  if (!is_initialized) {
    is_initialized = TRUE;
    elem_blk_id_save = ei[imtrx]->elem_blk_id;
  }

  /*
   * For convenience, while we are here, interpolate to find physical space
   * location using the mesh basis function.
   *
   * Generally, the other basis functions will give various different estimates
   * for position, depending on whether the shapes are mapped sub/iso/super
   * parametrically...
   */
  for (i = 0; i < VIM; i++) {
    fv->x[i] = 0.;
    fv_old->x[i] = 0.;
  }

  /*
   * NOTE: pdim is the number of coordinates, which may differ from
   * the element dimension (dim), as for shell elements!
   */
  for (i = 0; i < pdim; i++) {
    if (DeformingMesh) {
      for (k = 0; k < ei[upd->matrix_index[R_MESH1]]->dof[R_MESH1]; k++) {
        node = ei[upd->matrix_index[R_MESH1]]->dof_list[R_MESH1][k];

        index = Proc_Elem_Connect[Proc_Connect_Ptr[ei[upd->matrix_index[R_MESH1]]->ielem] +
                                  node];

        fv->x[i] += (Coor[i][index] + *esp->d[i][k]) * bf[R_MESH1]->phi[k];

        fv_old->x[i] +=
            (Coor[i][index] + *esp_old->d[i][k]) * bf[R_MESH1]->phi[k];
      }
    } else {
      for (k = 0; k < mdof; k++) {
        node = ei[imtrx]->dof_list[ShapeVar][k];

        index = Proc_Elem_Connect[Proc_Connect_Ptr[ei[imtrx]->ielem] + node];

        fv->x[i] += Coor[i][index] * MapBf->phi[k];

        fv_old->x[i] = fv->x[i];
      }
    }
  }

  /*
   * Elemental Jacobian is now affected by mesh displacement of nodes
   * from their initial nodal point coordinates...
   */

  for (i = 0; i < dim; i++) {
    for (j = 0; j < pdim; j++) {
      if (DeformingMesh) {
        for (k = 0; k < ei[upd->matrix_index[R_MESH1]]->dof[R_MESH1]; k++) {
          node = ei[upd->matrix_index[R_MESH1]]->dof_list[R_MESH1][k];

          index =
              Proc_Elem_Connect[Proc_Connect_Ptr[ei[upd->matrix_index[R_MESH1]]->ielem] +
                                node];

          MapBf->J[i][j] +=
              (Coor[j][index] + *esp->d[j][k]) * bf[R_MESH1]->dphidxi[k][i];
        }
      } else {
        for (k = 0; k < mdof; k++) {
          node = ei[imtrx]->dof_list[ShapeVar][k];
          index = Proc_Elem_Connect[Proc_Connect_Ptr[ei[imtrx]->ielem] + node];
          MapBf->J[i][j] += Coor[j][index] * bf[ShapeVar]->dphidxi[k][i];
        }
      }
    }
  }

  /*
   * Now copy mapping Jacobians from the MapBf basis function to the other Basis
   * functions
   *
   */

  for (t = 0; t < Num_Basis_Functions; t++) {
    if (t != si) {
      for (i = 0; i < dim; i++) {
        for (j = 0; j < pdim; j++) {
          bfd[t]->J[i][j] = MapBf->J[i][j];
        }
      }
    }
  }

  /*
   * Here we do the Jacobian, Inverse Jacobian and Jacobian determinant
   * only once for the basis functions identified as being the mapping
   * basis function.  At the end, we copy the results over to the
   * other order interpolants.
   */

  /* But first we pull a fast one.  If this is a 3D shell element, we up dim to dim+1 so
   *  that the 3D case is executed, but only after we populate the 3rd column of J with 
   * arbitrary nonzero constants so as to keep J full rank
   
   * AMC - This might work for the gradients once they are in the plane of the element,
   * but the gradients with respect to global coordinate system do indeed depend on
   * the values assigned to the right most column of the 
   * Jacobian of the mapping (MapBf->J).
   * 
   * There are problems that depend on the rightmost column being the values
   * computed above, so use this block with caution.
   */
  if(elem_shape == SHELL
     || elem_shape == TRISHELL
     || (mp->ehl_integration_kind == SIK_S))
    {
      dim++;
      for (t = 0; t < Num_Basis_Functions; t++)
	{ 
	  for (j = 0; j < pdim; j++)
	    {
	      bfd[t]->J[pd->Num_Dim-1][j] = MapBf->J[pd->Num_Dim-1][j] = (j+1)*1.0;
	    }
	}

    /*Real Quick check on Jacobian to make sure this arbitrary assignment
     *didn't screw things up. Note that the detJ in the shell case can be
     *negative, but it is important to point out that we are not using it for
     *for integration, but only as a crutch for inversion of J */
      if (pd->Num_Dim == 3) {
      MapBf->detJ = MapBf->J[0][0] * ( MapBf->J[1][1] * MapBf->J[2][2]
				       -MapBf->J[1][2] * MapBf->J[2][1])
	- MapBf->J[0][1] * ( MapBf->J[1][0] * MapBf->J[2][2]
			     -MapBf->J[2][0] * MapBf->J[1][2])
	+ MapBf->J[0][2] * ( MapBf->J[1][0] * MapBf->J[2][1]
			     -MapBf->J[2][0] * MapBf->J[1][1]);
      }
      if (pd->Num_Dim == 2) {
        MapBf->detJ = MapBf->J[0][0] * MapBf->J[1][1]
	- MapBf->J[0][1] * MapBf->J[1][0];
      }

    if (fabs(MapBf->detJ) < 1.e-10) {
      zero_detJ = TRUE;
#ifdef PARALLEL
      fprintf(stderr, "\nP_%d: Uh-oh, detJ =  %e\n", ProcID, fabs(MapBf->detJ));
#else
      fprintf(stderr, "\n Uh-oh, detJ =  %e\n", fabs(MapBf->detJ));
#endif
      return (2);
    }
  }

  /* Compute inverse of Jacobian for only the MapBf right now */

  /*
   * Wiggly mesh derivatives..
   */
  switch (dim) {
  case 1:
    /*
     * NOTE: This is the case for 1D shell elements, which must be
     * handled differently from other element types! Oddly,for a 2D shell in
     * a 3D problem, we skip to case 3.  Long story, but it is cleaner that way.
     */

    sum = 0.0;
    for (j = 0; j < pdim; j++) {
      sum += MapBf->J[0][j] * MapBf->J[0][j];
      if (fabs(MapBf->J[0][j]) > DBL_SMALL) {
        // MapBf->B[0][j] = 1.0 / MapBf->J[0][j];
        MapBf->B[j][0] = 1.0 / MapBf->J[0][j];
      }
    }
    MapBf->detJ = sqrt(sum);

    if (DeformingMesh) {
      for (j = 0; j < pdim; j++) {
        for (k = 0; k < ei[upd->matrix_index[R_MESH1]]->dof[R_MESH1]; k++) {
          MapBf->dJ[0][j][j][k] = bf[R_MESH1]->dphidxi[k][0];
          MapBf->d_det_J_dm[j][k] =
              bf[R_MESH1]->dphidxi[k][0] * MapBf->J[0][j] / MapBf->detJ;
          if (fabs(MapBf->J[0][j]) > DBL_SMALL) {
            MapBf->dB[j][0][j][k] =
                -bf[R_MESH1]->dphidxi[k][0] / (MapBf->J[0][j] * MapBf->J[0][j]);
            // MapBf->dB[0][j] [j][k] = -bf[R_MESH1]->dphidxi[k][0]
            //	/ (MapBf->J[0][j] * MapBf->J[0][j]);
          }
        }
      }
    }
    break;

    case 2:
      dim = ei[pg->imtrx]->ielem_dim;
      MapBf->detJ    =  MapBf->J[0][0] * MapBf->J[1][1]
	- MapBf->J[0][1] * MapBf->J[1][0];

    MapBf->B[0][0] = MapBf->J[1][1] / MapBf->detJ;
    MapBf->B[0][1] = -MapBf->J[0][1] / MapBf->detJ;
    MapBf->B[1][0] = -MapBf->J[1][0] / MapBf->detJ;
    MapBf->B[1][1] = MapBf->J[0][0] / MapBf->detJ;

    /*
     * Derivatives of elemental Jacobian matrix with respect
     * to each degree of freedom for mesh displacement components.
     */

    if (DeformingMesh) {
      for (i = 0; i < dim; i++) {
        for (j = 0; j < pdim; j++) {
          for (n = 0; n < ei[upd->matrix_index[R_MESH1]]->dof[R_MESH1]; n++) {
            MapBf->dJ[i][j][j][n] = bf[R_MESH1]->dphidxi[n][i];
          }
        }
      }

      /*
       * Derivatives of the determinant of the elemental Jacobian
       * matrix with respect to each degree of freedom for
       * mesh displacement.
       */

      for (k = 0; k < pdim; k++) {
        for (n = 0; n < ei[upd->matrix_index[R_MESH1]]->dof[R_MESH1]; n++) {
          MapBf->d_det_J_dm[k][n] = MapBf->dJ[0][0][k][n] * MapBf->J[1][1] +
                                    MapBf->J[0][0] * MapBf->dJ[1][1][k][n] -
                                    MapBf->dJ[0][1][k][n] * MapBf->J[1][0] -
                                    MapBf->J[0][1] * MapBf->dJ[1][0][k][n];
        }
      }

      /*
       * Derivatives of entries of the inverse of the elemental
       * Jacobian matrix with respect to each degree of freedom for
       * each component of mesh displacement.
       */

      f = 1. / (MapBf->detJ);

      for (k = 0; k < pdim; k++) {
        for (n = 0; n < ei[upd->matrix_index[R_MESH1]]->dof[R_MESH1]; n++) {
          g = -(f * f) * MapBf->d_det_J_dm[k][n];
          MapBf->dB[0][0][k][n] =
              MapBf->dJ[1][1][k][n] * f + MapBf->J[1][1] * g;

          MapBf->dB[0][1][k][n] =
              -MapBf->dJ[0][1][k][n] * f - MapBf->J[0][1] * g;

          MapBf->dB[1][0][k][n] =
              -MapBf->dJ[1][0][k][n] * f - MapBf->J[1][0] * g;

          MapBf->dB[1][1][k][n] =
              MapBf->dJ[0][0][k][n] * f + MapBf->J[0][0] * g;
        }
      }
    }

    break;

  case 3:

    /* Now that we are here, reset dim for the shell case */
    dim = ei[imtrx]->ielem_dim;

    MapBf->detJ = MapBf->J[0][0] * (MapBf->J[1][1] * MapBf->J[2][2] -
                                    MapBf->J[1][2] * MapBf->J[2][1]) -
                  MapBf->J[0][1] * (MapBf->J[1][0] * MapBf->J[2][2] -
                                    MapBf->J[2][0] * MapBf->J[1][2]) +
                  MapBf->J[0][2] * (MapBf->J[1][0] * MapBf->J[2][1] -
                                    MapBf->J[2][0] * MapBf->J[1][1]);

    MapBf->B[0][0] =
        (MapBf->J[1][1] * MapBf->J[2][2] - MapBf->J[2][1] * MapBf->J[1][2]) /
        (MapBf->detJ);

    MapBf->B[0][1] =
        -(MapBf->J[0][1] * MapBf->J[2][2] - MapBf->J[2][1] * MapBf->J[0][2]) /
        (MapBf->detJ);

    MapBf->B[0][2] =
        (MapBf->J[0][1] * MapBf->J[1][2] - MapBf->J[1][1] * MapBf->J[0][2]) /
        (MapBf->detJ);

    MapBf->B[1][0] =
        -(MapBf->J[1][0] * MapBf->J[2][2] - MapBf->J[2][0] * MapBf->J[1][2]) /
        (MapBf->detJ);

    MapBf->B[1][1] =
        (MapBf->J[0][0] * MapBf->J[2][2] - MapBf->J[2][0] * MapBf->J[0][2]) /
        (MapBf->detJ);

    MapBf->B[1][2] =
        -(MapBf->J[0][0] * MapBf->J[1][2] - MapBf->J[1][0] * MapBf->J[0][2]) /
        (MapBf->detJ);

    MapBf->B[2][0] =
        (MapBf->J[1][0] * MapBf->J[2][1] - MapBf->J[1][1] * MapBf->J[2][0]) /
        (MapBf->detJ);

    MapBf->B[2][1] =
        -(MapBf->J[0][0] * MapBf->J[2][1] - MapBf->J[2][0] * MapBf->J[0][1]) /
        (MapBf->detJ);

    MapBf->B[2][2] =
        (MapBf->J[0][0] * MapBf->J[1][1] - MapBf->J[1][0] * MapBf->J[0][1]) /
        (MapBf->detJ);

    if (DeformingMesh) {
      /*
       * Derivatives of elemental Jacobian matrix with respect
       * to each degree of freedom for mesh displacement components.
       */

      for (i = 0; i < dim; i++) {
        for (j = 0; j < pdim; j++) {
          for (n = 0; n < ei[upd->matrix_index[R_MESH1]]->dof[R_MESH1]; n++) {
            MapBf->dJ[i][j][j][n] = bf[R_MESH1]->dphidxi[n][i];
          }
        }
      }

      /*
       * Derivatives of the determinant of the elemental Jacobian
       * matrix with respect to each degree of freedom for
       * mesh displacement.
       */
      for (k = 0; k < pdim; k++) {
        for (n = 0; n < mdof; n++) {
          MapBf->d_det_J_dm[k][n] =
              MapBf->dJ[0][0][k][n] * (MapBf->J[1][1] * MapBf->J[2][2] -
                                       MapBf->J[1][2] * MapBf->J[2][1])

              + MapBf->J[0][0] * (MapBf->dJ[1][1][k][n] * MapBf->J[2][2] +
                                  MapBf->J[1][1] * MapBf->dJ[2][2][k][n] -
                                  MapBf->dJ[1][2][k][n] * MapBf->J[2][1] -
                                  MapBf->J[1][2] * MapBf->dJ[2][1][k][n])

              - MapBf->dJ[0][1][k][n] * (MapBf->J[1][0] * MapBf->J[2][2] -
                                         MapBf->J[2][0] * MapBf->J[1][2])

              - MapBf->J[0][1] * (MapBf->dJ[1][0][k][n] * MapBf->J[2][2] +
                                  MapBf->J[1][0] * MapBf->dJ[2][2][k][n] -
                                  MapBf->dJ[2][0][k][n] * MapBf->J[1][2] -
                                  MapBf->J[2][0] * MapBf->dJ[1][2][k][n])

              + MapBf->dJ[0][2][k][n] * (MapBf->J[1][0] * MapBf->J[2][1] -
                                         MapBf->J[2][0] * MapBf->J[1][1])

              + MapBf->J[0][2] * (MapBf->dJ[1][0][k][n] * MapBf->J[2][1] +
                                  MapBf->J[1][0] * MapBf->dJ[2][1][k][n] -
                                  MapBf->dJ[2][0][k][n] * MapBf->J[1][1] -
                                  MapBf->J[2][0] * MapBf->dJ[1][1][k][n]);
        }
      }

      /*
       * Derivatives of entries of the inverse of the
       * elemental Jacobian matrix with respect to
       * mesh displacement components at different dofs
       * of the element
       */
      f = 1. / (MapBf->detJ);

      for (k = 0; k < pdim; k++) {
        for (n = 0; n < ei[upd->matrix_index[R_MESH1]]->dof[R_MESH1]; n++) {
          g = -(f * f) * MapBf->d_det_J_dm[k][n];

          MapBf->dB[0][0][k][n] = (MapBf->dJ[1][1][k][n] * MapBf->J[2][2] +
                                   MapBf->J[1][1] * MapBf->dJ[2][2][k][n] -
                                   MapBf->dJ[2][1][k][n] * MapBf->J[1][2] -
                                   MapBf->J[2][1] * MapBf->dJ[1][2][k][n]) *
                                      f +
                                  (MapBf->J[1][1] * MapBf->J[2][2] -
                                   MapBf->J[2][1] * MapBf->J[1][2]) *
                                      g;

          MapBf->dB[0][1][k][n] = -(MapBf->dJ[0][1][k][n] * MapBf->J[2][2] +
                                    MapBf->J[0][1] * MapBf->dJ[2][2][k][n] -
                                    MapBf->dJ[2][1][k][n] * MapBf->J[0][2] -
                                    MapBf->J[2][1] * MapBf->dJ[0][2][k][n]) *
                                      f -
                                  (MapBf->J[0][1] * MapBf->J[2][2] -
                                   MapBf->J[2][1] * MapBf->J[0][2]) *
                                      g;

          MapBf->dB[0][2][k][n] = (MapBf->dJ[0][1][k][n] * MapBf->J[1][2] +
                                   MapBf->J[0][1] * MapBf->dJ[1][2][k][n] -
                                   MapBf->dJ[1][1][k][n] * MapBf->J[0][2] -
                                   MapBf->J[1][1] * MapBf->dJ[0][2][k][n]) *
                                      f +
                                  (MapBf->J[0][1] * MapBf->J[1][2] -
                                   MapBf->J[1][1] * MapBf->J[0][2]) *
                                      g;

          MapBf->dB[1][0][k][n] = -(MapBf->dJ[1][0][k][n] * MapBf->J[2][2] +
                                    MapBf->J[1][0] * MapBf->dJ[2][2][k][n] -
                                    MapBf->dJ[2][0][k][n] * MapBf->J[1][2] -
                                    MapBf->J[2][0] * MapBf->dJ[1][2][k][n]) *
                                      f -
                                  (MapBf->J[1][0] * MapBf->J[2][2] -
                                   MapBf->J[2][0] * MapBf->J[1][2]) *
                                      g;

          MapBf->dB[1][1][k][n] = (MapBf->dJ[0][0][k][n] * MapBf->J[2][2] +
                                   MapBf->J[0][0] * MapBf->dJ[2][2][k][n] -
                                   MapBf->dJ[2][0][k][n] * MapBf->J[0][2] -
                                   MapBf->J[2][0] * MapBf->dJ[0][2][k][n]) *
                                      f +
                                  (MapBf->J[0][0] * MapBf->J[2][2] -
                                   MapBf->J[2][0] * MapBf->J[0][2]) *
                                      g;

          MapBf->dB[1][2][k][n] = -(MapBf->dJ[0][0][k][n] * MapBf->J[1][2] +
                                    MapBf->J[0][0] * MapBf->dJ[1][2][k][n] -
                                    MapBf->dJ[1][0][k][n] * MapBf->J[0][2] -
                                    MapBf->J[1][0] * MapBf->dJ[0][2][k][n]) *
                                      f -
                                  (MapBf->J[0][0] * MapBf->J[1][2] -
                                   MapBf->J[1][0] * MapBf->J[0][2]) *
                                      g;

          MapBf->dB[2][0][k][n] = (MapBf->dJ[1][0][k][n] * MapBf->J[2][1] +
                                   MapBf->J[1][0] * MapBf->dJ[2][1][k][n] -
                                   MapBf->dJ[1][1][k][n] * MapBf->J[2][0] -
                                   MapBf->J[1][1] * MapBf->dJ[2][0][k][n]) *
                                      f +
                                  (MapBf->J[1][0] * MapBf->J[2][1] -
                                   MapBf->J[1][1] * MapBf->J[2][0]) *
                                      g;

          MapBf->dB[2][1][k][n] = -(MapBf->dJ[0][0][k][n] * MapBf->J[2][1] +
                                    MapBf->J[0][0] * MapBf->dJ[2][1][k][n] -
                                    MapBf->dJ[2][0][k][n] * MapBf->J[0][1] -
                                    MapBf->J[2][0] * MapBf->dJ[0][1][k][n]) *
                                      f -
                                  (MapBf->J[0][0] * MapBf->J[2][1] -
                                   MapBf->J[2][0] * MapBf->J[0][1]) *
                                      g;

          MapBf->dB[2][2][k][n] = (MapBf->dJ[0][0][k][n] * MapBf->J[1][1] +
                                   MapBf->J[0][0] * MapBf->dJ[1][1][k][n] -
                                   MapBf->dJ[1][0][k][n] * MapBf->J[0][1] -
                                   MapBf->J[1][0] * MapBf->dJ[0][1][k][n]) *
                                      f +
                                  (MapBf->J[0][0] * MapBf->J[1][1] -
                                   MapBf->J[1][0] * MapBf->J[0][1]) *
                                      g;
        }
      }
    }
    break;

  default:
    EH(-1, "Bad dim.");
    break;
  }

  /* Now copy all this stuff from the MapBf to the other basis functions */

  for (t = 0; t < Num_Basis_Functions; t++) {
    if (t != si) {
      bfd[t]->detJ = MapBf->detJ;
      for (i = 0; i < pdim; i++) {
        for (j = 0; j < pdim; j++) {
          bfd[t]->B[i][j] = MapBf->B[i][j];

          if (DeformingMesh) {
            for (k = 0; k < pdim; k++) {
              for (n = 0; n < ei[upd->matrix_index[R_MESH1]]->dof[R_MESH1]; n++) {
                bfd[t]->d_det_J_dm[k][n] = MapBf->d_det_J_dm[k][n];

                bfd[t]->dJ[i][j][k][n] = MapBf->dJ[i][j][k][n];
                bfd[t]->dB[i][j][k][n] = MapBf->dB[i][j][k][n];
              }
            }
          }
        }
      }
    }
  }

  return (status);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void calc_surf_tangent(
    const int ielem,                /* current element number               */
    const int iconnect_ptr,         /* Ptr to beginning of connectivity
                                       list for current element      */
    const int nodes_per_elem,       /* num nodes in the element    */
    const int ielem_surf_dim,       /* physical dimension of element
                                     * surface (0, 1, 2)           */
    const int num_nodes_on_side,    /* num nodes side of elem   */
    const int local_elem_node_id[]) /* local element node
                                     * numbers on elem side */

/*
 * Function which calculates the components of unit surface tangent
 * vector(s) in the current element. NB!
 * This routine rests on having
 * already computed snormal and dsnormal_dx from routine
 * calc_surf_normal() routine.
 *
 *       Author:          P. R. Schunk (1511)
 *
 * Nomenclature:
 *
 *  fv->snormal[icoord]
 *  fv->stangent[itan][icoord]
 *     where
 *       itan = the id of the tangent (there are two for 3D probs)
 *       icoord = coordinate dimension of the tangent or normal vector
 *
 *        In 2D normal x tangent = k, these vectors satisfy the rh rule.
 *
 *  fv->dsnormal_dx[icoord][icoord_dx][idof_dx]
 *  fv->dstangent_dx[itan][icoord][icoord_dx][idof_dx]
 *     where
 *       icoord_dx = coordinate dimension of the displacement unknown
 *       idof_dx  = degree of freedom of the displacement unknown
 *                  within the element
 */
{
  int dalpha_dx[MAX_PDIM];
  int i, id, j, k, l, inode, ldof, ShapeVar;
  double alpha;

  ShapeVar = pd->ShapeVar;

  memset(fv->stangent, 0, 2 * MAX_PDIM * sizeof(double));
  memset(fv->dstangent_dx, 0, 2 * MAX_PDIM * MAX_PDIM * MDE * sizeof(double));

  switch (ielem_surf_dim) {

  case 0:
    fv->stangent[0][0] = 1.;
    fv->stangent[1][0] = 0.;
    for (j = 0; j < nodes_per_elem; j++)
      fv->dstangent_dx[0][0][0][j] = 0.;
    break;

  case 1:
    /*
     *  In 2D n x t = k, these vectors satisfy the rh rule.
     */
    fv->stangent[0][0] = -fv->snormal[1];
    fv->stangent[0][1] =  fv->snormal[0];
    fv->stangent[1][2] =  1.0;
    for (j=0 ; j < nodes_per_elem; j++) {
      fv->dstangent_dx[0][0][0][j]=0.;
      fv->dstangent_dx[0][0][1][j]=0.;
      fv->dstangent_dx[0][1][0][j]=0.;
      fv->dstangent_dx[0][1][1][j]=0.;
    }
    for (i = 0; i < num_nodes_on_side; i++) {
      id = (int)local_elem_node_id[i];
      inode = Proc_Elem_Connect[iconnect_ptr + id];
      ldof = ei[pd->mi[ShapeVar]]->ln_to_dof[ShapeVar][id];

      if (num_varType_at_node(inode, MESH_DISPLACEMENT1)) {
        fv->dstangent_dx[0][0][1][ldof] = -fv->dsnormal_dx[1][1][ldof];
        fv->dstangent_dx[0][0][0][ldof] = -fv->dsnormal_dx[1][0][ldof];
        fv->dstangent_dx[0][1][1][ldof] = fv->dsnormal_dx[0][1][ldof];
        fv->dstangent_dx[0][1][0][ldof] = fv->dsnormal_dx[0][0][ldof];
      }
    }
    break;

  case 2:

    /* First tangent defined as
     *    t1 = [(n X i) + (n X j)] / norm[n X (i+j)]
     */
    alpha = sqrt(2. * fv->snormal[2] * fv->snormal[2] +
                 SQUARE(fv->snormal[0] - fv->snormal[1]));
    fv->stangent[0][0] = -fv->snormal[2] / alpha;
    fv->stangent[0][1] = fv->snormal[2] / alpha;
    fv->stangent[0][2] = (fv->snormal[0] - fv->snormal[1]) / alpha;
    if (SQUARE(alpha) <= 1.e-7) {
      fprintf(stderr,
              "Warning - norm of tangent on surface of element %d Jacobian is "
              "dangerously small %e\n",
              ielem, alpha);
    }

    /* First tangent defined by seed vector gives t1 = tseed - n(tseed DOT n)
     * normalized */
    /* not yet implemented */

    /* if in an element adjacent to one edge, make 1st tangent be perpendicular
     * to that edge */
    /* not yet implemented */

    /* Second tangent defined as t2 = n X t1 automatically normalized because t1
     * and n are unit vectors */

    fv->stangent[1][0] = (fv->snormal[1] * fv->stangent[0][2] -
                          fv->snormal[2] * fv->stangent[0][1]);
    fv->stangent[1][1] = (fv->snormal[2] * fv->stangent[0][0] -
                          fv->snormal[0] * fv->stangent[0][2]);
    fv->stangent[1][2] = (fv->snormal[0] * fv->stangent[0][1] -
                          fv->snormal[1] * fv->stangent[0][0]);

    for (j = 0; j < nodes_per_elem; j++) {
      for (k = 0; k < ielem_surf_dim + 1; k++) {
        for (l = 0; l < ielem_surf_dim + 1; l++) {
          fv->dstangent_dx[0][k][l][j] = 0.;
          fv->dstangent_dx[1][k][l][j] = 0.;
        }
      }
    }
    for (i = 0; i < num_nodes_on_side; i++) {
      id = (int)local_elem_node_id[i];
      inode = Proc_Elem_Connect[iconnect_ptr + id];
      ldof = ei[pd->mi[ShapeVar]]->ln_to_dof[ShapeVar][id];
      if (num_varType_at_node(inode, MESH_DISPLACEMENT1)) {
        for (j = 0; j < ielem_surf_dim + 1; j++) {
          dalpha_dx[j] = (2. * fv->snormal[2] * fv->dsnormal_dx[2][j][ldof] +
                          (fv->snormal[0] - fv->snormal[1]) *
                              (fv->dsnormal_dx[0][j][ldof] -
                               fv->dsnormal_dx[1][j][ldof])) /
                         alpha;

          fv->dstangent_dx[0][0][j][ldof] =
              (-alpha * fv->dsnormal_dx[2][j][ldof] +
               fv->snormal[2] * dalpha_dx[j]) /
              SQUARE(alpha);

          fv->dstangent_dx[0][1][j][ldof] =
              (alpha * fv->dsnormal_dx[2][j][ldof] -
               fv->snormal[2] * dalpha_dx[j]) /
              SQUARE(alpha);

          fv->dstangent_dx[0][2][j][ldof] =
              (alpha *
                   (fv->dsnormal_dx[0][j][ldof] - fv->dsnormal_dx[1][j][ldof]) -
               (fv->snormal[0] - fv->snormal[1]) * dalpha_dx[j]) /
              SQUARE(alpha);

          /* derivatives of nXt1 */
          fv->dstangent_dx[1][0][j][ldof] =
              (fv->dsnormal_dx[1][j][ldof] * fv->stangent[0][2] -
               fv->dsnormal_dx[2][j][ldof] * fv->stangent[0][1] +
               fv->snormal[1] * fv->dstangent_dx[0][2][j][ldof] -
               fv->snormal[2] * fv->dstangent_dx[0][1][j][ldof]);

          fv->dstangent_dx[1][1][j][ldof] =
              (fv->dsnormal_dx[2][j][ldof] * fv->stangent[0][0] -
               fv->dsnormal_dx[0][j][ldof] * fv->stangent[0][2] +
               fv->snormal[2] * fv->dstangent_dx[0][0][j][ldof] -
               fv->snormal[0] * fv->dstangent_dx[0][2][j][ldof]);

          fv->dstangent_dx[1][2][j][ldof] =
              (fv->dsnormal_dx[0][j][ldof] * fv->stangent[0][1] -
               fv->dsnormal_dx[1][j][ldof] * fv->stangent[0][0] +
               fv->snormal[0] * fv->dstangent_dx[0][1][j][ldof] -
               fv->snormal[1] * fv->dstangent_dx[0][0][j][ldof]);
        }
      }
    }
    break;
  }
} /* calc_surf_tangent */
/********************************************************************** */
/*********************************************************************** */
/*********************************************************************** */

void calc_tangent_from_seed(struct Rotation_Vectors *tangent,
                            struct Rotation_Vectors *normal,
                            const double seed[DIM], const int dim)
/*
 * Function which calculates the components of unit surface tangent vector
 * in the current element. NB! This routine rests on having already computed
 * snormal and dsnormal_dx from routine calc_surf_normal routine.
 *       Author:          P. R. Schunk (1511)
 *
 * Initialize some arrays for fixed grid defaults. 950306pas
 *
 * This form of the routine calculates tangents from a seed vector RAC 08/12/96
 * This form of the routine calculates tangents from a seed vector RAC 09/04/96
 *
 */

{
  int p, q, b, j;

  tangent->ok = 1;

  /* First tanget defined by seed vector gives
   *  t1 = tseed - n(tseed DOT n) normalized
   */
  for (p = 0; p < dim; p++) {
    tangent->vector[p] = seed[p];
    for (q = 0; q < dim; q++) {
      tangent->vector[p] -= normal->vector[p] * normal->vector[q] * seed[q];
    }
  }

  if (af->Assemble_Jacobian && pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
    tangent->d_vector_n = normal->d_vector_n;
    for (j = 0; j < tangent->d_vector_n; j++) {
      tangent->d_vector_J[j] = normal->d_vector_J[j];

      for (b = 0; b < dim; b++) {
        for (p = 0; p < dim; p++) {
          for (q = 0; q < dim; q++) {
            tangent->d_vector_dx[p][b][j] -=
                normal->d_vector_dx[p][b][j] * normal->vector[q] * seed[q] +
                normal->vector[p] * normal->d_vector_dx[q][b][j] * seed[q];
          }
        }
      }
    }
  }

  /* now normalize it */
  simple_normalize_vector(tangent, dim);

  return;
} /* calc_tangent_from_seed */
/*********************************************************************** */
/*********************************************************************** */
/*ARGSUSED*/
void calc_tangent_along_basis(struct Rotation_Vectors *tangent,
                              struct Rotation_Vectors *normal, const int dim,
                              const int id_side, const int num_nodes_on_side,
                              const int local_elem_node_id[MAX_NODES_PER_SIDE])
/*
 * Function which calculates the components of unit surface tangent vector
 * in the current element. NB! This routine rests on having already computed
 * snormal and dsnormal_dx from routine calc_surf_normal routine.
 *       Author:          P. R. Schunk (1511)
 *
 * Initialize some arrays for fixed grid defaults. 950306pas
 *
 * This form of the routine calculates tangents from a seed vector RAC 08/12/96
 * This form of the routine calculates tangents from a seed vector RAC 09/04/96
 *
 */

{
  /* TAB certifies that this function conforms to the exo/patran side numbering
   * convention 11/10/98. */

  int i_basis = -1;
  int p;
  int i, id, inode, ldof;

  const int ShapeVar = pd->ShapeVar;
  /***************************** execution begins
   * ******************************/

  /* default to first basis direction for calculating unseeded tangent */
  switch (id_side) {
  case 1:
    i_basis = 0;
    break;
  case 2:
    i_basis = 1;
    break;
  case 3:
    i_basis = 0;
    break;
  case 4:
    i_basis = 1;
    break;
  case 5:
    i_basis = 0;
    break;
  case 6:
    i_basis = 0;
    break;
  }

  tangent->ok = 1;
  /* find direction along the parametric basis  */
  tangent->d_vector_n = num_nodes_on_side;
  for (i = 0; i < tangent->d_vector_n; i++) {
    id = (int)local_elem_node_id[i];
    inode = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + id];
    tangent->d_vector_J[i] = inode;
    ldof = ei[pd->mi[ShapeVar]]->ln_to_dof[ShapeVar][id];
    if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
      if (num_varType_at_node(inode, MESH_DISPLACEMENT1)) {
        for (p = 0; p < dim; p++) {
          tangent->vector[p] += bf[ShapeVar]->dphidxi[ldof][i_basis] *
                                (Coor[p][inode] + *esp->d[p][ldof]);
        }
      }
    } else {
      for (p = 0; p < dim; p++) {
        tangent->vector[p] +=
            bf[ShapeVar]->dphidxi[ldof][i_basis] * Coor[p][inode];
      }
    }
  }

  if (af->Assemble_Jacobian && pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
    for (i = 0; i < tangent->d_vector_n; i++) {
      id = (int)local_elem_node_id[i];
      inode = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + id];
      ldof = ei[pg->imtrx]->ln_to_dof[ShapeVar][id];
      if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
        if (num_varType_at_node(inode, MESH_DISPLACEMENT1)) {
          for (p = 0; p < dim; p++) {
            tangent->d_vector_dx[p][p][i] +=
                bf[ShapeVar]->dphidxi[ldof][i_basis];
          }
        }
      }
    }
  }

  /* now normalize it */
  simple_normalize_vector(tangent, dim);

  return;
} /* calc_tangent_from_basis */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void cross_vectors(struct Rotation_Vectors *tangent2,
                   struct Rotation_Vectors *tangent1,
                   struct Rotation_Vectors *normal, const int dim)
/*
 * Function which calculates the cross product of two vectors
 *
 *  tangent2 = normal X tangent1
 *
 * RAC 09/04/96
 */

{
  int p, b, j;
  double sum;

  /***************************** execution begins
   * ******************************/
  tangent2->ok = 1;

  /* Second Tangent T2 = N X T1 */
  tangent2->vector[0] = -(normal->vector[1] * tangent1->vector[2] -
                          normal->vector[2] * tangent1->vector[1]);
  tangent2->vector[1] = -(normal->vector[2] * tangent1->vector[0] -
                          normal->vector[0] * tangent1->vector[2]);
  if (dim == 3) {
    tangent2->vector[2] = -(normal->vector[0] * tangent1->vector[1] -
                            normal->vector[1] * tangent1->vector[0]);
  }

  if (af->Assemble_Jacobian && pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
    tangent2->d_vector_n = normal->d_vector_n;
    for (j = 0; j < tangent2->d_vector_n; j++) {
      tangent2->d_vector_J[j] = normal->d_vector_J[j];

      for (b = 0; b < dim; b++) {

        /* Second Tangent T2 = N X T1 */
        tangent2->d_vector_dx[0][b][j] =
            -(normal->d_vector_dx[1][b][j] * tangent1->vector[2] -
              normal->d_vector_dx[2][b][j] * tangent1->vector[1]) -
            (normal->vector[1] * tangent1->d_vector_dx[2][b][j] -
             normal->vector[2] * tangent1->d_vector_dx[1][b][j]);
        tangent2->d_vector_dx[1][b][j] =
            -(normal->d_vector_dx[2][b][j] * tangent1->vector[0] -
              normal->d_vector_dx[0][b][j] * tangent1->vector[2]) -
            (normal->vector[2] * tangent1->d_vector_dx[0][b][j] -
             normal->vector[0] * tangent1->d_vector_dx[2][b][j]);
        if (dim == 3) {
          tangent2->d_vector_dx[2][b][j] =
              -(normal->d_vector_dx[0][b][j] * tangent1->vector[1] -
                normal->d_vector_dx[1][b][j] * tangent1->vector[0]) -
              (normal->vector[0] * tangent1->d_vector_dx[1][b][j] -
               normal->vector[1] * tangent1->d_vector_dx[0][b][j]);
        }
      }
    }
  }

  for (p = 0, sum = 0.; p < dim; p++) {
    sum += tangent1->vector[p] * tangent2->vector[p];
  }
  if (fabs(sum) > 1e-8) {
    printf(
        "2nd Tangent not perpendicular to 1st tangent (%g) vector %g,%g,%g\n",
        sum, tangent1->vector[0], tangent1->vector[1], tangent1->vector[2]);
  }
  for (p = 0, sum = 0.; p < dim; p++) {
    sum += normal->vector[p] * tangent2->vector[p];
  }
  if (fabs(sum) > 1e-8) {
    printf("2nd Tangent not perpendicular to normal (%g) vector %g,%g,%g\n",
           sum, tangent1->vector[0], tangent1->vector[1], tangent1->vector[2]);
  }
  return;
} /* cross_vectors */
/********************************************************************** */
/********************************************************************** */

void calc_unseeded_edge_tangents(
    struct Rotation_Vectors *tangent,
    const int iconnect_ptr,        /* Ptr to beginning of
                                    * connectivity list for
                                    * current element */
    const int dim,                 /* physical dimension of the
                                    * surface of the element
                                    * i.e., (0, 1, 2) */
    const int id_side,             /* ID of the side of the
                                    * element*/
    const int id_edge,             /* ID of the edge of the
                                    * element*/
    const int num_nodes_on_edge,   /* number of nodes
                                    * on the edge of
                                    * the element */
    const int edge_elem_node_id[], /* vector of the
                                    * local element
                                    * node numbers
                                    * on the edge of
                                    * the element */
    const int param_dir)           /* parametric direction of
                                    * the edge in local
                                    * coordinates */

/*
 * Function which calculates the components of unit surface tangent vector
 * in the current element. NB! This routine rests on having already computed
 * snormal and dsnormal_dx from routine calc_surf_normal routine.
 *       Author:          P. R. Schunk (1511)
 *
 * Initialize some arrays for fixed grid defaults. 950306pas
 *
 * This form of the routine calculates tangents along an edge of an element
 * from a seed vector RAC 08/13/96
 *
 */

{
  int j, inode, ldof;
  const int ShapeVar = pd->ShapeVar;
  const int DeformingMesh = pd->e[pg->imtrx][R_MESH1];
  int p;
  int j_id;
  double sign;

  tangent->ok = 1;

  /* determine sign on parametric curve to make Normal, line tangent, and
   * binormal (outward pointing normal to egde) a right-handed basis */
  sign = 0.;
  switch (id_edge) {
  case (1):
    switch (id_side) {
    case (4):
      sign = 1.;
      break;
    case (5):
      sign = -1.;
      break;
    default:
      EH(-1, "Side not connected to edge");
      break;
    }
    break;

  case (2):
    switch (id_side) {
    case (5):
      sign = 1.;
      break;
    case (2):
      sign = -1.;
      break;
    default:
      EH(-1, "Side not connected to edge");
      break;
    }
    break;

  case (3):
    switch (id_side) {
    case (5):
      sign = 1.;
      break;
    case (1):
      sign = -1.;
      break;
    default:
      EH(-1, "Side not connected to edge");
      break;
    }
    break;

  case (4):
    switch (id_side) {
    case (3):
      sign = 1.;
      break;
    case (5):
      sign = -1.;
      break;
    default:
      EH(-1, "Side not connected to edge");
      break;
    }
    break;

  case (5):
    switch (id_side) {
    case (6):
      sign = 1.;
      break;
    case (4):
      sign = -1.;
      break;
    default:
      EH(-1, "Side not connected to edge");
      break;
    }
    break;

  case (6):
    switch (id_side) {
    case (2):
      sign = 1.;
      break;
    case (6):
      sign = -1.;
      break;
    default:
      EH(-1, "Side not connected to edge");
      break;
    }
    break;

  case (7):
    switch (id_side) {
    case (1):
      sign = 1.;
      break;
    case (6):
      sign = -1.;
      break;
    default:
      EH(-1, "Side not connected to edge");
      break;
    }
    break;

  case (8):
    switch (id_side) {
    case (6):
      sign = 1.;
      break;
    case (3):
      sign = -1.;
      break;
    default:
      EH(-1, "Side not connected to edge");
      break;
    }
    break;

  case (9):
    switch (id_side) {
    case (1):
      sign = 1.;
      break;
    case (4):
      sign = -1.;
      break;
    default:
      EH(-1, "Side not connected to edge");
      break;
    }
    break;

  case (10):
    switch (id_side) {
    case (4):
      sign = 1.;
      break;
    case (3):
      sign = -1.;
      break;
    default:
      EH(-1, "Side not connected to edge");
      break;
    }
    break;

  case (11):
    switch (id_side) {
    case (2):
      sign = 1.;
      break;
    case (1):
      sign = -1.;
      break;
    default:
      EH(-1, "Side not connected to edge");
      break;
    }
    break;

  case (12):
    switch (id_side) {
    case (3):
      sign = 1.;
      break;
    case (2):
      sign = -1.;
      break;
    default:
      EH(-1, "Side not connected to edge");
      break;
    }
    break;

  default:
    EH(-1, "Edge not found");
    break;
  }

  /* -------------------------------------------------
   *   an EDGE
   * -------------------------------------------------*/

  /* calculate unit vector in parametric direction of local coord param_dir */
  for (j = 0; j < num_nodes_on_edge; j++) {
    j_id = (int)edge_elem_node_id[j];
    inode = Proc_Elem_Connect[iconnect_ptr + j_id];
    ldof = ei[pd->mi[ShapeVar]]->ln_to_dof[ShapeVar][j_id];
    if (DeformingMesh) {
      if (num_varType_at_node(inode, MESH_DISPLACEMENT1)) {
        for (p = 0; p < dim; p++) {
          tangent->vector[p] += sign * bf[ShapeVar]->dphidxi[ldof][param_dir] *
                                (Coor[p][inode] + *esp->d[p][ldof]);
        }
      }
    } else {
      for (p = 0; p < dim; p++) {
        tangent->vector[p] +=
            sign * bf[ShapeVar]->dphidxi[ldof][param_dir] * Coor[p][inode];
      }
    }
  }

  if (af->Assemble_Jacobian && pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
    tangent->d_vector_n = num_nodes_on_edge;
    for (j = 0; j < num_nodes_on_edge; j++) {
      j_id = (int)edge_elem_node_id[j];
      inode = Proc_Elem_Connect[iconnect_ptr + j_id];
      tangent->d_vector_J[j] = inode;
      ldof = ei[pg->imtrx]->ln_to_dof[ShapeVar][j_id];
      if (num_varType_at_node(inode, MESH_DISPLACEMENT1)) {
        for (p = 0; p < dim; p++) {
          tangent->d_vector_dx[p][p][j] +=
              sign * bf[ShapeVar]->dphidxi[ldof][param_dir];
        }
      }
    }
  }

  /* now normalize it */
  simple_normalize_vector(tangent, dim);

} /* calc_unseeded_edge_tangents */
/*****************************************************************************/
/********************************************************************** */

void calc_unseeded_edge_tangents_TET(
    struct Rotation_Vectors *tangent,
    const int iconnect_ptr,        /* Ptr to beginning of
                                    * connectivity list for
                                    * current element */
    const int dim,                 /* physical dimension of the
                                    * surface of the element
                                    * i.e., (0, 1, 2) */
    const int id_side,             /* ID of the side of the
                                    * element*/
    const int id_edge,             /* ID of the edge of the
                                    * element*/
    const int num_nodes_on_edge,   /* number of nodes
                                    * on the edge of
                                    * the element */
    const int edge_elem_node_id[], /* vector of the
                                    * local element
                                    * node numbers
                                    * on the edge of
                                    * the element */
    const int param_dir)           /* parametric direction of
                                    * the edge in local
                                    * coordinates */

/*
 * Function which calculates the edge tangent for tet elements using information
 * from calls to find_id_edge_TET and find_id_side
 *       Author:          P. R. Schunk (1516)
 *
 * This form of the routine calculates tangents along an edge of an element
 * from a seed vector RAC 08/13/96 PRS 02/12/2012
 *
 */

{
  int j, inode, ldof;
  const int ShapeVar = pd->ShapeVar;
  const int DeformingMesh = pd->e[pg->imtrx][R_MESH1];
  int p;
  int j_id;
  double sign;

  tangent->ok = 1;

  /* determine sign on parametric curve to make Normal, line tangent, and
   * binormal (outward pointing normal to egde) a right-handed basis */
  sign = 0.;
  switch (id_edge) {
  case (1):
    switch (id_side) {
    case (1):
      sign = -1.;
      break;
    case (4):
      sign = 1.;
      break;
    default:
      EH(-1, "Side not connected to edge");
      break;
    }
    break;

  case (2):
    switch (id_side) {
    case (2):
      sign = 1.;
      break;
    case (4):
      sign = -1.;
      break;
    default:
      EH(-1, "Side not connected to edge");
      break;
    }
    break;

  case (3):
    switch (id_side) {
    case (3):
      sign = 1.;
      break;
    case (4):
      sign = -1.;
      break;
    default:
      EH(-1, "Side not connected to edge");
      break;
    }
    break;

  case (4):
    switch (id_side) {
    case (1):
      sign = 1.;
      break;
    case (3):
      sign = -1.;
      break;
    default:
      EH(-1, "Side not connected to edge");
      break;
    }
    break;

  case (5):
    switch (id_side) {
    case (1):
      sign = 1.;
      break;
    case (2):
      sign = -1.;
      break;
    default:
      EH(-1, "Side not connected to edge");
      break;
    }
    break;

  case (6):
    switch (id_side) {
    case (2):
      sign = 1.;
      break;
    case (3):
      sign = -1.;
      break;
    default:
      EH(-1, "Side not connected to edge");
      break;
    }
    break;

  default:
    EH(-1, "Edge not found for TETS");
    break;
  }

  /* -------------------------------------------------
   *   an EDGE
   * -------------------------------------------------*/

  /* calculate unit vector in parametric direction of local coord param_dir */
  int param_dir_tet = 0;
  for (j = 0; j < num_nodes_on_edge; j++) {
    j_id = (int)edge_elem_node_id[j];
    inode = Proc_Elem_Connect[iconnect_ptr + j_id];
    ldof = ei[pd->mi[ShapeVar]]->ln_to_dof[ShapeVar][j_id];
    if (id_edge == 1 || id_edge == 3 || id_edge == 4) {
      if (DeformingMesh) {
        if (num_varType_at_node(inode, MESH_DISPLACEMENT1)) {
          for (p = 0; p < dim; p++) {
            tangent->vector[p] += sign *
                                  bf[ShapeVar]->dphidxi[ldof][param_dir] *
                                  (Coor[p][inode] + *esp->d[p][ldof]);
          }
        }
      } else {
        for (p = 0; p < dim; p++) {
          tangent->vector[p] +=
              sign * bf[ShapeVar]->dphidxi[ldof][param_dir] * Coor[p][inode];
        }
      }
    } else if (id_edge == 2) {
      if (j_id == 1) {
        param_dir_tet = 0;
      }
      if (j_id == 2) {
        param_dir_tet = 1;
        sign *= -1;
      }
      if (DeformingMesh) {
        if (num_varType_at_node(inode, MESH_DISPLACEMENT1)) {

          for (p = 0; p < dim; p++) {
            tangent->vector[p] += sign *
                                  bf[ShapeVar]->dphidxi[ldof][param_dir_tet] *
                                  (Coor[p][inode] + *esp->d[p][ldof]);
          }
        }
      } else {
        for (p = 0; p < dim; p++) {
          tangent->vector[p] += sign *
                                bf[ShapeVar]->dphidxi[ldof][param_dir_tet] *
                                Coor[p][inode];
        }
      }
    } else if (id_edge == 6) {
      if (j_id == 2) {
        param_dir_tet = 1;
      }
      if (j_id == 3) {
        param_dir_tet = 2;
        sign *= -1;
      }
      if (DeformingMesh) {
        if (num_varType_at_node(inode, MESH_DISPLACEMENT1)) {

          for (p = 0; p < dim; p++) {
            tangent->vector[p] += sign *
                                  bf[ShapeVar]->dphidxi[ldof][param_dir_tet] *
                                  (Coor[p][inode] + *esp->d[p][ldof]);
          }
        }
      } else {
        for (p = 0; p < dim; p++) {
          tangent->vector[p] += sign *
                                bf[ShapeVar]->dphidxi[ldof][param_dir_tet] *
                                Coor[p][inode];
        }
      }
    } else if (id_edge == 5) {
      if (j_id == 1) {
        param_dir_tet = 0;
      }
      if (j_id == 3) {
        param_dir_tet = 2;
        sign *= -1;
      }
      if (DeformingMesh) {
        if (num_varType_at_node(inode, MESH_DISPLACEMENT1)) {

          for (p = 0; p < dim; p++) {
            tangent->vector[p] += sign *
                                  bf[ShapeVar]->dphidxi[ldof][param_dir_tet] *
                                  (Coor[p][inode] + *esp->d[p][ldof]);
          }
        }
      } else {
        for (p = 0; p < dim; p++) {
          tangent->vector[p] += sign *
                                bf[ShapeVar]->dphidxi[ldof][param_dir_tet] *
                                Coor[p][inode];
        }
      }
    } else {
      EH(-1, "Cannot find id_edge in unseeded_TET");
    }
  }

  if (af->Assemble_Jacobian && pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
    tangent->d_vector_n = num_nodes_on_edge;
    for (j = 0; j < num_nodes_on_edge; j++) {
      j_id = (int)edge_elem_node_id[j];
      inode = Proc_Elem_Connect[iconnect_ptr + j_id];
      tangent->d_vector_J[j] = inode;
      ldof = ei[pd->mi[ShapeVar]]->ln_to_dof[ShapeVar][j_id];
      if (num_varType_at_node(inode, MESH_DISPLACEMENT1)) {
        if (id_edge == 1 || id_edge == 3 || id_edge == 4) {
          for (p = 0; p < dim; p++) {
            tangent->d_vector_dx[p][p][j] +=
                sign * bf[ShapeVar]->dphidxi[ldof][param_dir];
          }
        } else if (id_edge == 2) {
          if (j_id == 1) {
            param_dir_tet = 0;
          }
          if (j_id == 2) {
            param_dir_tet = 1;
            sign *= -1;
          }
          for (p = 0; p < dim; p++) {
            tangent->d_vector_dx[p][p][j] +=
                sign * bf[ShapeVar]->dphidxi[ldof][param_dir_tet];
          }
        } else if (id_edge == 6) {
          if (j_id == 2) {
            param_dir_tet = 1;
          }
          if (j_id == 3) {
            param_dir_tet = 2;
            sign *= -1;
          }
          for (p = 0; p < dim; p++) {
            tangent->d_vector_dx[p][p][j] +=
                sign * bf[ShapeVar]->dphidxi[ldof][param_dir_tet];
          }
        } else if (id_edge == 5) {
          if (j_id == 1) {
            param_dir_tet = 0;
          }
          if (j_id == 3) {
            param_dir_tet = 2;
            sign *= -1;
          }
          for (p = 0; p < dim; p++) {
            tangent->d_vector_dx[p][p][j] +=
                sign * bf[ShapeVar]->dphidxi[ldof][param_dir_tet];
          }
        }
      }
    }
  }

  /* now normalize it */
  simple_normalize_vector(tangent, dim);

} /* calc_unseeded_edge_tangents_TET */
/*****************************************************************************/
/*****************************************************************************/

void simple_normalize_vector(struct Rotation_Vectors *vector, const int dim) {
  /* Local variables */
  int p, q, j;
  double sum, sqsum;
  double d_sum[MAX_PDIM][MDE];

  for (p = 0, sum = 0.; p < dim; p++)
    sum += vector->vector[p] * vector->vector[p];
  sqsum = sqrt(sum);

  if (af->Assemble_Jacobian && pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
    for (j = 0; j < vector->d_vector_n; j++) {
      for (p = 0; p < dim; p++) {
        d_sum[p][j] = 0.;
        for (q = 0; q < dim; q++) {
          d_sum[p][j] += 2. * vector->vector[q] * vector->d_vector_dx[q][p][j];
        }
      }
    }

    for (j = 0; j < vector->d_vector_n; j++) {
      for (q = 0; q < dim; q++) {
        for (p = 0; p < dim; p++) {
          vector->d_vector_dx[q][p][j] /= sqsum;
          vector->d_vector_dx[q][p][j] -=
              0.5 * vector->vector[q] * d_sum[p][j] / sum / sqsum;
        }
      }
    }
  }

  for (p = 0; p < dim; p++)
    vector->vector[p] /= sqsum;

  /* check the sum */
  for (p = 0, sum = 0.; p < dim; p++)
    sum += vector->vector[p] * vector->vector[p];
  if (fabs(sum - 1.) > 1e-10)
    EH(-1, "Bad normalization");
  return;
}
/********************************************************************************/
/********************************************************************************/
/********************************************************************************/

int load_bf_grad(void)

/********************************************************************************
 *
 *  load_bf_grad():
 *
 * This supersedes the functionality of map_shape(), loading
 * up the physical space gradients of each different basis
 * function.
 *
 * Also, load up the physical space gradient of basis functions
 * multiplied by the different coordinate unit vectors. This is
 * easy for Cartesian coordinates, but not for general curvilinear
 * coordinates.
 *
 * Default to looping over all 3 indeces (0,1,2) even for two dimensional
 * problems (by using VIM). This should help us avoid bad memory references
 * for axisymmetric problems where hoop stress terms are required.
 *
 * Added general scale factor capability so that grad_phi is converted from
 * raw form to actual form. Thus,
 *
 * Incoming: grad_phi= J^-1 grad_phi_local
 *
 * Suppose we have Cartesian coordinate system:
 * then: q1=x, q2=y, q3=z and scale factors h1=h2=h3=1
 *
 * Axisymmetric: q1=z, q2=r, q3=theta and scale factors are 1,1,(1/r)
 * Note that the last scale factor will have nontrivial contributions to
 * mesh derivatives.
 * Also, note that the unit vectors in axisymmetric case have nontrivial
 * contributions, from d(e_r)/dtheta
 *
 *
 *           d(phi)
 *           ------
 *            d q1
 *
 *           d(phi)
 *           ------   =
 *            d q2
 *
 *           d(phi)
 *           ------
 *            d q3
 * Returns:
 *		0 -- if things went ok
 *	       -1 -- if something went wrong
 *
 *
 * Revised: Tue Feb 14 09:07 MST 1995 pasacki@sandia.gov
 *
 * Revised: 9/2/99 MMH Added curly things
 *
 * Note: Axisymmetric convention (z,r,theta) are coordinates.
 ********************************************************************************/
{
  int i, k, a, b, p, dofs = 0, v, vi, status, siz;
  struct Basis_Functions *bfv;

#ifdef DO_NOT_UNROLL
  int WIM;

  if ((pd->CoordinateSystem == CARTESIAN) ||
      (pd->CoordinateSystem == CYLINDRICAL)) {
    WIM = pd->Num_Dim;
  } else {
    WIM = VIM;
  }
#endif

  status = 0;

  /* zero array for initialization */
  /*  v_length = DIM*DIM*DIM*MDE;
      init_vec_value(zero_array, 0., v_length); */

  for (b = 0; b < Num_Basis_Functions; b++) {
    /*
     * This is a kludge until I can fix up the bfd's with the same
     * local node dof_list that the ei's currently enjoy.
     *
     * Find a variable "v" that has the current interpolation
     * scheme in the material
     */
    int imtrx;
    for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
      for (vi = 0, v = -1; v == -1 && vi < MAX_VARIABLE_TYPES; vi++) {

        if (pd->i[imtrx][vi] == bfd[b]->interpolation) {
          v = vi;
        }
      }

      /*
       *  Don't calculate basis function if interpolation doesn't
       *  correspond to a variable
       *  OR if element shape doesn't match the current element block.
       */
      if (v != -1 && bfd[b]->element_shape == ei[pg->imtrx]->ielem_shape) {
        bfv = bf[v];
        dofs = ei[imtrx]->dof[v];

        /* initialize variables */
        siz = pd->Num_Dim * MDE * sizeof(double);
        /* memset(&(bfv->d_phi[0][0]),0,siz); */

        /*
         * First load up components of the *raw* derivative vector "d_phi"
         */
        switch (pd->Num_Dim) {
        case 1:
          for (i = 0; i < dofs; i++) {
            bfv->d_phi[i][0] = (bfv->B[0][0] * bfv->dphidxi[i][0] +
                                bfv->B[0][1] * bfv->dphidxi[i][1]);
            bfv->d_phi[i][1] = 0.0;
            bfv->d_phi[i][2] = 0.0;
          }
          break;
        case 2:
          for (i = 0; i < dofs; i++) {
            bfv->d_phi[i][0] = (bfv->B[0][0] * bfv->dphidxi[i][0] +
                                bfv->B[0][1] * bfv->dphidxi[i][1]);
            bfv->d_phi[i][1] = (bfv->B[1][0] * bfv->dphidxi[i][0] +
                                bfv->B[1][1] * bfv->dphidxi[i][1]);
            bfv->d_phi[i][2] = 0.0;
          }
          break;
        case 3:
          for (i = 0; i < dofs; i++) {
            bfv->d_phi[i][0] = (bfv->B[0][0] * bfv->dphidxi[i][0] +
                                bfv->B[0][1] * bfv->dphidxi[i][1] +
                                bfv->B[0][2] * bfv->dphidxi[i][2]);
            bfv->d_phi[i][1] = (bfv->B[1][0] * bfv->dphidxi[i][0] +
                                bfv->B[1][1] * bfv->dphidxi[i][1] +
                                bfv->B[1][2] * bfv->dphidxi[i][2]);
            bfv->d_phi[i][2] = (bfv->B[2][0] * bfv->dphidxi[i][0] +
                                bfv->B[2][1] * bfv->dphidxi[i][1] +
                                bfv->B[2][2] * bfv->dphidxi[i][2]);
          }
          break;
        default:
          EH(-1, "Unexpected Dimension");
          break;
        }

        /*
         * Now, patch up the physical space gradient of this prototype
         * scalar function so scale factors are included.
         */

        /*	memset(&(bfv->grad_phi[0][0]),0,size1);  */

        for (i = 0; i < dofs; i++) {
#ifdef DO_NOT_UNROLL
          for (p = 0; p < WIM; p++) {
            bfv->grad_phi[i][p] = (bfv->d_phi[i][p]) / (fv->h[p]);
          }
#else
          bfv->grad_phi[i][0] = (bfv->d_phi[i][0]) / (fv->h[0]);
          bfv->grad_phi[i][1] = (bfv->d_phi[i][1]) / (fv->h[1]);
          if (VIM == 3)
            bfv->grad_phi[i][2] = (bfv->d_phi[i][2]) / (fv->h[2]);
#endif
        }

        /*
         * Second load up components of the tensor grad(phi_i e_a) where
         * e_a are the unit vectors in this orthogonal coordinate system.
         *
         * In Cartesian coordinates:
         *
         * 		grad(phi_i e_a)_j,k = d(phi_i)/dx_j * delta(a, k)
         *
         * Thus, these can be formed easily from the components of the
         * gradient of a scalar that were formed above.
         *
         * Now, use VIM to loop over 3 indeces in vectors and tensors
         * so that all terms for axisymmetric problems, including the notorious
         * hoop stress term, are properly accounted for.
         */

        /*
         * Initialize...
         */

        /* memset(&(bfv->grad_phi_e[0][0][0][0]),0,size2); */

        for (i = 0; i < dofs; i++) {
#ifdef DO_NOT_UNROLL
          for (p = 0; p < VIM; p++) {
            for (a = 0; a < VIM; a++) {
              for (q = 0; q < VIM; q++) {
                if (q == a)
                  bfv->grad_phi_e[i][a][p][a] = bfv->grad_phi[i][p];
                else
                  bfv->grad_phi_e[i][a][p][q] = 0.0;
              }
            }
          }
#else

          bfv->grad_phi_e[i][0][0][0] = bfv->grad_phi[i][0];
          bfv->grad_phi_e[i][0][0][1] = 0.0;
          bfv->grad_phi_e[i][1][1][0] = 0.0;
          bfv->grad_phi_e[i][1][1][1] = bfv->grad_phi[i][1];
          bfv->grad_phi_e[i][0][1][0] = bfv->grad_phi[i][1];
          bfv->grad_phi_e[i][0][1][1] = 0.0;
          bfv->grad_phi_e[i][1][0][0] = 0.0;
          bfv->grad_phi_e[i][1][0][1] = bfv->grad_phi[i][0];

          if (VIM == 3) {
            bfv->grad_phi_e[i][2][0][0] = 0.0;
            bfv->grad_phi_e[i][2][0][1] = 0.0;
            bfv->grad_phi_e[i][2][0][2] = bfv->grad_phi[i][0];

            bfv->grad_phi_e[i][2][1][0] = 0.0;
            bfv->grad_phi_e[i][2][1][1] = 0.0;
            bfv->grad_phi_e[i][2][1][2] = bfv->grad_phi[i][1];

            bfv->grad_phi_e[i][2][2][0] = 0.0;
            bfv->grad_phi_e[i][2][2][1] = 0.0;
            bfv->grad_phi_e[i][2][2][2] = bfv->grad_phi[i][2];

            bfv->grad_phi_e[i][1][0][2] = 0.0;
            bfv->grad_phi_e[i][1][1][2] = 0.0;
            bfv->grad_phi_e[i][1][2][0] = 0.0;
            bfv->grad_phi_e[i][1][2][1] = bfv->grad_phi[i][2];
            bfv->grad_phi_e[i][1][2][2] = 0.0;

            bfv->grad_phi_e[i][0][0][2] = 0.0;
            bfv->grad_phi_e[i][0][1][2] = 0.0;
            bfv->grad_phi_e[i][0][2][0] = bfv->grad_phi[i][2];
            bfv->grad_phi_e[i][0][2][1] = 0.0;
            bfv->grad_phi_e[i][0][2][2] = 0.0;
          }
#endif

          /* } */

          if (pd->CoordinateSystem != CARTESIAN) {

            /*   for ( i=0; i<dofs; i++)
                 { */

#ifdef DO_NOT_UNROLL
            for (a = 0; a < WIM; a++) {
              for (p = 0; p < VIM; p++) {
                for (q = 0; q < VIM; q++) {
                  bfv->grad_phi_e[i][a][p][q] +=
                      bfv->phi[i] * fv->grad_e[a][p][q];
                }
              }
            }
#else
            /*
              bfv->grad_phi_e[i][#] [0][0] += bfv->phi[i] * fv->grad_e[#][0][0];
              bfv->grad_phi_e[i][#] [1][1] += bfv->phi[i] * fv->grad_e[#][1][1];
              bfv->grad_phi_e[i][#] [0][1] += bfv->phi[i] * fv->grad_e[#][0][1];
              bfv->grad_phi_e[i][#] [1][0] += bfv->phi[i] * fv->grad_e[#][1][0];

              if(VIM == 3 )
              {
              bfv->grad_phi_e[i][#] [2][2] += bfv->phi[i] * fv->grad_e[#][2][2];
              bfv->grad_phi_e[i][#] [2][1] += bfv->phi[i] * fv->grad_e[#][2][1];
              bfv->grad_phi_e[i][#] [2][0] += bfv->phi[i] * fv->grad_e[#][2][0];
              bfv->grad_phi_e[i][#] [1][2] += bfv->phi[i] * fv->grad_e[#][1][2];
              bfv->grad_phi_e[i][#] [0][2] += bfv->phi[i] * fv->grad_e[#][0][2];
              }
            */

            bfv->grad_phi_e[i][0][0][0] += bfv->phi[i] * fv->grad_e[0][0][0];
            bfv->grad_phi_e[i][0][1][1] += bfv->phi[i] * fv->grad_e[0][1][1];
            bfv->grad_phi_e[i][0][0][1] += bfv->phi[i] * fv->grad_e[0][0][1];
            bfv->grad_phi_e[i][0][1][0] += bfv->phi[i] * fv->grad_e[0][1][0];

            bfv->grad_phi_e[i][1][0][0] += bfv->phi[i] * fv->grad_e[1][0][0];
            bfv->grad_phi_e[i][1][1][1] += bfv->phi[i] * fv->grad_e[1][1][1];
            bfv->grad_phi_e[i][1][0][1] += bfv->phi[i] * fv->grad_e[1][0][1];
            bfv->grad_phi_e[i][1][1][0] += bfv->phi[i] * fv->grad_e[1][1][0];

            if (VIM == 3) {
              bfv->grad_phi_e[i][0][2][2] += bfv->phi[i] * fv->grad_e[0][2][2];
              bfv->grad_phi_e[i][0][2][1] += bfv->phi[i] * fv->grad_e[0][2][1];
              bfv->grad_phi_e[i][0][2][0] += bfv->phi[i] * fv->grad_e[0][2][0];
              bfv->grad_phi_e[i][0][1][2] += bfv->phi[i] * fv->grad_e[0][1][2];
              bfv->grad_phi_e[i][0][0][2] += bfv->phi[i] * fv->grad_e[0][0][2];

              bfv->grad_phi_e[i][1][2][2] += bfv->phi[i] * fv->grad_e[1][2][2];
              bfv->grad_phi_e[i][1][2][1] += bfv->phi[i] * fv->grad_e[1][2][1];
              bfv->grad_phi_e[i][1][2][0] += bfv->phi[i] * fv->grad_e[1][2][0];
              bfv->grad_phi_e[i][1][1][2] += bfv->phi[i] * fv->grad_e[1][1][2];
              bfv->grad_phi_e[i][1][0][2] += bfv->phi[i] * fv->grad_e[1][0][2];

              bfv->grad_phi_e[i][2][0][0] += bfv->phi[i] * fv->grad_e[2][0][0];
              bfv->grad_phi_e[i][2][1][1] += bfv->phi[i] * fv->grad_e[2][1][1];
              bfv->grad_phi_e[i][2][0][1] += bfv->phi[i] * fv->grad_e[2][0][1];
              bfv->grad_phi_e[i][2][1][0] += bfv->phi[i] * fv->grad_e[2][1][0];

              bfv->grad_phi_e[i][2][2][2] += bfv->phi[i] * fv->grad_e[2][2][2];
              bfv->grad_phi_e[i][2][2][1] += bfv->phi[i] * fv->grad_e[2][2][1];
              bfv->grad_phi_e[i][2][2][0] += bfv->phi[i] * fv->grad_e[2][2][0];
              bfv->grad_phi_e[i][2][1][2] += bfv->phi[i] * fv->grad_e[2][1][2];
              bfv->grad_phi_e[i][2][0][2] += bfv->phi[i] * fv->grad_e[2][0][2];
            }

#endif
          }
        }

        /* MMH: Construct the curl of phi_i * e_a.
         *
         *   curl_phi_e[i][a][p] = e_p . curl(phi_i * e_a)
         *                           __
         * or:                   = ( \/ x (phi_i e_a) )_p
         *
         *                           /\
         *                           ||
         *                           ||
         *                            How Geekly!
         */

        /* This is currently only used as a post-processing variable
         * by me.  Since it takes so much of load_bf_grad()'s time,
         * we're only going to do real work if the vorticty vector
         * post processing variable is turned on.  This isn't even the
         * right place to do that (since it is only post-processing),
         * but some day someone is going to want a vorticity-based
         * variable.  I had one, but it is not presently being used.
         */

        if (CURL_V != -1) {
          siz = DIM * DIM * MDE * sizeof(double);
          memset(&(bfv->curl_phi_e[0][0][0]), 0, siz);

          for (i = 0; i < dofs; i++) {
            /* MMH: Always compute all three components of the vorticity
             * vector.  If we are really in 2D, then only the third
             * component of vorticity is useful.
             */
            for (a = 0; a < DIM; a++) {
              for (p = 0; p < DIM; p++) {
                for (k = 0; k < DIM; k++) { /* VIM */
                  bfv->curl_phi_e[i][a][p] +=
                      permute(p, k, a) * bfv->grad_phi[i][k];
                }
                bfv->curl_phi_e[i][a][p] += bfv->phi[i] * fv->curl_e[a][p];
              }
            }
          }
        }

      } /* end of if v */
    }   /* end of basis function loop. */
  }

  return (status);
}
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/* load_bf_mesh_derivs() -- load mesh derivatives of bf gradients in realspace
 *
 * needs:
 *		bf[var]->... filled with regular gradients
 *		             filled with B, and dB
 *
 *
 * input: (assume the following parts of the bf arrays are loaded)
 * ------
 *		ei			-- Element Indeces information
 *		pd			-- Problem Description information
 *
 *	bf[v]->
 *
 *		grad_phi[i][j]		--	e_j . grad(phi_i)
 *
 *		grad_phi_e[i][j][k][l]	--	e_k e_l : grad(phi_i e_j )
 *
 *		dB[i][j] [k][l]		--	d ( B_i,j ) / d d_k,l
 *
 * where:
 *      B   == inverse Jacobian matrix
 *	e_i == unit vector in this orthogonal coordinate system
 *
 * output:
 * --------
 *
 *	bf[v]->
 *              d_d_phi_dmesh[i][p] [b][j]   --  d ( d(phi_i)/dp )
 *						 -------------------
 *						 d ( d_b,j )
 *
 *
 *		d_grad_phi_dmesh[i][p] [b][j] -- d ( e_p . grad(phi_i) )
 *						 -----------------------
 *						 d ( d_b,j )
 *
 *		d_grad_phi_e_dmesh[i][j][k][l][m][n] --
 *
 *						d ( e_k e_l : grad(phi_i e_j) )
 *						-------------------------------
 *						d ( d_m,n )
 *
 *
 *	dphixdm[i][j][k][n] - deriv of (deriv of phi[i] wrt to global coord[j])
 *                            wrt mesh disp. k at node n
 *
 * Note!!!!
 *		d( )/dmesh_b,j is just a derivative while
 *		grad( )	       says something about the coordinate system
 *
 *	Implication: make sure when you use d ( ) / d ( d_theta, j)
 *		     that it's really what you need and dont'
 *
 *
 * return values:
 *        0                 - everything went OK
 *       -1                 - something went wrong
 *
 *
 * Created:	Mon Feb 20 15:47 MST 1995 pasacki@sandia.gov
 */
int load_bf_mesh_derivs(void) {
  int a, b, p, q, i, j, bix, v, siz;
  int dim;       /* number of spatial dimensions */
  int wim;       /* looping variable useful for elliptical polar
                    coordinates, equal to the number velocity unknowns */
  int dimNonSym; /* Number of dimensions of the curvilinear coordinate
                  * system which has nonzero gradients. Coordinates such
                  * as the theta component in the cylindrical coordinates
                  * which have assumed symmetry don't count. Therefore
                  * dimNonSym = 2 for cylindrical and swirling */
  int vdofs;     /* degrees of freedom for variable that is */
  /* interpolated using bfl */
  int mdofs; /* degrees of freedom for mesh displacement */
  /* unknowns interpolated using bfm */
  int mn = ei[pg->imtrx]->mn;
  dbl f, g[DIM] = {0};      /* Temporary variables for convenience. */
  dbl g2[DIM] = {0};
  dbl phi_m[MDE], phi_l[MDE]; /* load shapefunctions into local variables */

  BASIS_FUNCTIONS_STRUCT *bfl; /* Basis function of current interest */
  BASIS_FUNCTIONS_STRUCT *bfm; /* Basis function for mesh displacement */
  BASIS_FUNCTIONS_STRUCT *bf_ptr;
  int imtrx;

  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    /*
     * If we aren't solving for displacements, then do a quick
     * return
     */
    if (!ei[imtrx]->deforming_mesh) {
      continue;
    }

    dim = pd->Num_Dim;
    dimNonSym = dim;
    if (pd->CoordinateSystem == CARTESIAN ||
        pd->CoordinateSystem == CYLINDRICAL ||
        pd->CoordinateSystem == PROJECTED_CARTESIAN) {
      wim = dim;
    } else if (pd->CoordinateSystem == SWIRLING ||
             pd->CoordinateSystem == CARTESIAN_2pt5D) {
      wim = 3;
    } else {
      /* MMH: What makes it here??? */
      wim = VIM;
    }

    /*
     * Preload the number of degrees of freedom in the mesh
     * for the current element, mdofs, and the pointer to the
     * basis function structure, bfm
     */
    mdofs = ei[imtrx]->dof[R_MESH1];
    bfm = bf[MESH_DISPLACEMENT1];
    /*
     * preload mesh shape functions
     */
    for (j = 0; j < mdofs; j++) {
      phi_m[j] = bfm->phi[j];
    }

    for (p = 0; p < VIM; p++) {
      g[p] = 1.0 / fv->h[p];
      g2[p] = g[p] * g[p];
    }

    /*
     * Loop over the number of different unique
     * basis functions defined in this material
     */
    for (bix = 0; bix < Num_Basis_Functions; bix++) {
      bf_ptr = bfd[bix];

      /*
       * Find a variable "v" that has the same interpolation as
       * the current unique basis function...then we can use
       * convenient bf[v] as an alias for bfd[bix]...
       */
      v = bf_ptr->Var_Type_MatID[mn];

      /*
       *  Don't calculate basis function if interpolation doesn't
       *  correspond to a variable
       *  OR if element shape doesn't match the current element block.
       */
      if (v != -1 && bf_ptr->element_shape == ei[imtrx]->ielem_shape) {
        /*
         * Store the number of dofs in the interpolation, vdofs.
         * Store bfl
         */
        vdofs = ei[imtrx]->dof[v];
        bfl = bf[v];

        /* preload shape functions */
        for (i = 0; i < vdofs; i++) {
          phi_l[i] = bfl->phi[i];
        }

        /*
         * Evaluate the sensitivity of derivatives of the fundamental scalar
         * quantity (the basis function, phi_i) in the coordinate directions
         * with respect to mesh displacement unknowns
         *
         * i.e.,
         *		d ( d_phi[i][p] )
         *		--------------------
         *		d ( d_b,j )
         */

        /* initialize variables */

        /* siz = sizeof(double)*DIM*DIM*MDE*MDE;
           memset(&(bfl->d_d_phi_dmesh[0][0][0][0]), 0, siz); */

        switch (dim) {
        case 2:
#ifdef DO_NOT_UNROLL
          for (i = 0; i < vdofs; i++) {
            for (p = 0; p < dim; p++) {
              memset(&(bfl->d_d_phi_dmesh[i][p][0][0]), 0,
                     mdofs * sizeof(double));
              memset(&(bfl->d_d_phi_dmesh[i][p][1][0]), 0,
                     mdofs * sizeof(double));
              for (j = 0; j < mdofs; j++) {
                bfl->d_d_phi_dmesh[i][p][0][j] +=
                    (bfl->dB[p][0][0][j] * bfl->dphidxi[i][0] +
                     bfl->dB[p][1][0][j] * bfl->dphidxi[i][1]);
                bfl->d_d_phi_dmesh[i][p][1][j] +=
                    (bfl->dB[p][0][1][j] * bfl->dphidxi[i][0] +
                     bfl->dB[p][1][1][j] * bfl->dphidxi[i][1]);
              }
            }
          }
#else
          for (i = 0; i < vdofs; i++) {
            for (j = 0; j < mdofs; j++) {

              bfl->d_d_phi_dmesh[i][0][0][j] =
                  (bfl->dB[0][0][0][j] * bfl->dphidxi[i][0] +
                   bfl->dB[0][1][0][j] * bfl->dphidxi[i][1]);

              bfl->d_d_phi_dmesh[i][1][0][j] =
                  (bfl->dB[1][0][0][j] * bfl->dphidxi[i][0] +
                   bfl->dB[1][1][0][j] * bfl->dphidxi[i][1]);

              bfl->d_d_phi_dmesh[i][0][1][j] =
                  (bfl->dB[0][0][1][j] * bfl->dphidxi[i][0] +
                   bfl->dB[0][1][1][j] * bfl->dphidxi[i][1]);
              bfl->d_d_phi_dmesh[i][1][1][j] =
                  (bfl->dB[1][0][1][j] * bfl->dphidxi[i][0] +
                   bfl->dB[1][1][1][j] * bfl->dphidxi[i][1]);
            }
          }
#endif
          break;
        case 3:

#ifdef DO_NOT_UNROLL
          for (i = 0; i < vdofs; i++) {
            for (p = 0; p < dim; p++) {
              memset(&(bfl->d_d_phi_dmesh[i][p][0][0]), 0,
                     mdofs * sizeof(double));
              memset(&(bfl->d_d_phi_dmesh[i][p][1][0]), 0,
                     mdofs * sizeof(double));
              memset(&(bfl->d_d_phi_dmesh[i][p][2][0]), 0,
                     mdofs * sizeof(double));

              for (j = 0; j < mdofs; j++) {
                bfl->d_d_phi_dmesh[i][p][0][j] +=
                    (bfl->dB[p][0][0][j] * bfl->dphidxi[i][0] +
                     bfl->dB[p][1][0][j] * bfl->dphidxi[i][1] +
                     bfl->dB[p][2][0][j] * bfl->dphidxi[i][2]);
                bfl->d_d_phi_dmesh[i][p][1][j] +=
                    (bfl->dB[p][0][1][j] * bfl->dphidxi[i][0] +
                     bfl->dB[p][1][1][j] * bfl->dphidxi[i][1] +
                     bfl->dB[p][2][1][j] * bfl->dphidxi[i][2]);
                bfl->d_d_phi_dmesh[i][p][2][j] +=
                    (bfl->dB[p][0][2][j] * bfl->dphidxi[i][0] +
                     bfl->dB[p][1][2][j] * bfl->dphidxi[i][1] +
                     bfl->dB[p][2][2][j] * bfl->dphidxi[i][2]);
              }
            }
          }
#else
          for (i = 0; i < vdofs; i++) {
            for (j = 0; j < mdofs; j++) {

              bfl->d_d_phi_dmesh[i][0][0][j] =
                  (bfl->dB[0][0][0][j] * bfl->dphidxi[i][0] +
                   bfl->dB[0][1][0][j] * bfl->dphidxi[i][1] +
                   bfl->dB[0][2][0][j] * bfl->dphidxi[i][2]);
              bfl->d_d_phi_dmesh[i][1][0][j] =
                  (bfl->dB[1][0][0][j] * bfl->dphidxi[i][0] +
                   bfl->dB[1][1][0][j] * bfl->dphidxi[i][1] +
                   bfl->dB[1][2][0][j] * bfl->dphidxi[i][2]);
              bfl->d_d_phi_dmesh[i][2][0][j] =
                  (bfl->dB[2][0][0][j] * bfl->dphidxi[i][0] +
                   bfl->dB[2][1][0][j] * bfl->dphidxi[i][1] +
                   bfl->dB[2][2][0][j] * bfl->dphidxi[i][2]);

              bfl->d_d_phi_dmesh[i][0][1][j] =
                  (bfl->dB[0][0][1][j] * bfl->dphidxi[i][0] +
                   bfl->dB[0][1][1][j] * bfl->dphidxi[i][1] +
                   bfl->dB[0][2][1][j] * bfl->dphidxi[i][2]);
              bfl->d_d_phi_dmesh[i][1][1][j] =
                  (bfl->dB[1][0][1][j] * bfl->dphidxi[i][0] +
                   bfl->dB[1][1][1][j] * bfl->dphidxi[i][1] +
                   bfl->dB[1][2][1][j] * bfl->dphidxi[i][2]);
              bfl->d_d_phi_dmesh[i][2][1][j] =
                  (bfl->dB[2][0][1][j] * bfl->dphidxi[i][0] +
                   bfl->dB[2][1][1][j] * bfl->dphidxi[i][1] +
                   bfl->dB[2][2][1][j] * bfl->dphidxi[i][2]);

              bfl->d_d_phi_dmesh[i][0][2][j] =
                  (bfl->dB[0][0][2][j] * bfl->dphidxi[i][0] +
                   bfl->dB[0][1][2][j] * bfl->dphidxi[i][1] +
                   bfl->dB[0][2][2][j] * bfl->dphidxi[i][2]);
              bfl->d_d_phi_dmesh[i][1][2][j] =
                  (bfl->dB[1][0][2][j] * bfl->dphidxi[i][0] +
                   bfl->dB[1][1][2][j] * bfl->dphidxi[i][1] +
                   bfl->dB[1][2][2][j] * bfl->dphidxi[i][2]);
              bfl->d_d_phi_dmesh[i][2][2][j] =
                  (bfl->dB[2][0][2][j] * bfl->dphidxi[i][0] +
                   bfl->dB[2][1][2][j] * bfl->dphidxi[i][1] +
                   bfl->dB[2][2][2][j] * bfl->dphidxi[i][2]);
            }
          }

#endif
          break;
        default:
          EH(-1, "Unknown dimension.");
          break;
        }

        /*
         * Evaluate the sensitivity of gradients of the fundamental scalar
         * quantity (the basisfunction, phi_i) in this coordinate system
         * with respect to mesh displacement unknowns
         *
         * i.e.,
         *		d ( grad_phi[i][p] )
         *		--------------------
         *		d ( d_b,j )
         *
         * Initialize the "ipbj" component, then sum over the two dot products
         * the two terms that contribute.
         */

        /* initialize variables */

        /* siz = sizeof(double)*DIM*DIM*MDE*MDE;
           memset(&(bfl->d_grad_phi_dmesh[0][0][0][0]), 0, siz); */
#ifdef DO_NOT_UNROLL
        for (i = 0; i < vdofs; i++) {
          for (p = 0; p < dimNonSym; p++) {
            for (b = 0; b < dimNonSym; b++) {
              memset(&(bfl->d_grad_phi_dmesh[i][p][b][0]), 0,
                     mdofs * sizeof(double));
              f = -g2[p] * (fv->hq[p][b]);
              for (j = 0; j < mdofs; j++) {
                for (q = 0; q < dimNonSym; q++) {
                  bfl->d_grad_phi_dmesh[i][p][b][j] +=
                      (f * phi_m[j] * bfl->B[p][q] +
                       g[p] * bfl->dB[p][q][b][j]) *
                      (bfl->dphidxi[i][q]);
                }
              }
            }
          }
        }
#else
        for (i = 0; i < vdofs; i++) {

          /*  memset( &(bfl->d_grad_phi_dmesh[i][0] [0][0]), 0,
             mdofs*sizeof(double) ); memset( &(bfl->d_grad_phi_dmesh[i][1]
             [1][0]), 0, mdofs*sizeof(double) ); memset(
             &(bfl->d_grad_phi_dmesh[i][1] [0][0]), 0, mdofs*sizeof(double) );
              memset( &(bfl->d_grad_phi_dmesh[i][0] [1][0]), 0,
             mdofs*sizeof(double) );

              if( dim == 3 )
              {
              memset( &(bfl->d_grad_phi_dmesh[i][2] [2][0]), 0,
             mdofs*sizeof(double) ); memset( &(bfl->d_grad_phi_dmesh[i][2]
             [1][0]), 0, mdofs*sizeof(double) ); memset(
             &(bfl->d_grad_phi_dmesh[i][2] [0][0]), 0, mdofs*sizeof(double) );
              memset( &(bfl->d_grad_phi_dmesh[i][1] [2][0]), 0,
             mdofs*sizeof(double) ); memset( &(bfl->d_grad_phi_dmesh[i][0]
             [2][0]), 0, mdofs*sizeof(double) );
              }
          */
          for (j = 0; j < mdofs; j++) {

            /* (p,b) = (0,0) */

            f = -g2[0] * (fv->hq[0][0]);

            bfl->d_grad_phi_dmesh[i][0][0][j] =
                (f * phi_m[j] * bfl->B[0][0] + g[0] * bfl->dB[0][0][0][j]) *
                (bfl->dphidxi[i][0]);

            bfl->d_grad_phi_dmesh[i][0][0][j] +=
                (f * phi_m[j] * bfl->B[0][1] + g[0] * bfl->dB[0][1][0][j]) *
                (bfl->dphidxi[i][1]);

            /* (p,b) = (0,1) */
            f = -g2[0] * (fv->hq[0][1]);

            bfl->d_grad_phi_dmesh[i][0][1][j] =
                (f * phi_m[j] * bfl->B[0][0] + g[0] * bfl->dB[0][0][1][j]) *
                (bfl->dphidxi[i][0]);

            bfl->d_grad_phi_dmesh[i][0][1][j] +=
                (f * phi_m[j] * bfl->B[0][1] + g[0] * bfl->dB[0][1][1][j]) *
                (bfl->dphidxi[i][1]);

            /* (p,b) = (1,0) */
            f = -g2[1] * (fv->hq[1][0]);

            bfl->d_grad_phi_dmesh[i][1][0][j] =
                (f * phi_m[j] * bfl->B[1][0] + g[1] * bfl->dB[1][0][0][j]) *
                (bfl->dphidxi[i][0]);

            bfl->d_grad_phi_dmesh[i][1][0][j] +=
                (f * phi_m[j] * bfl->B[1][1] + g[1] * bfl->dB[1][1][0][j]) *
                (bfl->dphidxi[i][1]);

            /* (p,b) = (1,1) */

            f = -g2[1] * (fv->hq[1][1]);
            bfl->d_grad_phi_dmesh[i][1][1][j] =
                (f * phi_m[j] * bfl->B[1][0] + g[1] * bfl->dB[1][0][1][j]) *
                (bfl->dphidxi[i][0]);

            bfl->d_grad_phi_dmesh[i][1][1][j] +=
                (f * phi_m[j] * bfl->B[1][1] + g[1] * bfl->dB[1][1][1][j]) *
                (bfl->dphidxi[i][1]);

            if (dimNonSym == 3) {
              /* (p,b) = (0,0) */
              bfl->d_grad_phi_dmesh[i][0][0][j] +=
                  (f * phi_m[j] * bfl->B[0][2] + g[0] * bfl->dB[0][2][0][j]) *
                  (bfl->dphidxi[i][2]);
              /* (p,b) = (0,1) */
              bfl->d_grad_phi_dmesh[i][0][1][j] +=
                  (f * phi_m[j] * bfl->B[0][2] + g[0] * bfl->dB[0][2][1][j]) *
                  (bfl->dphidxi[i][2]);
              /* (p,b) = (1,0) */
              bfl->d_grad_phi_dmesh[i][1][0][j] +=
                  (f * phi_m[j] * bfl->B[1][2] + g[1] * bfl->dB[1][2][0][j]) *
                  (bfl->dphidxi[i][2]);
              /* (p,b) = (1,1) */
              bfl->d_grad_phi_dmesh[i][1][1][j] +=
                  (f * phi_m[j] * bfl->B[1][2] + g[1] * bfl->dB[1][2][1][j]) *
                  (bfl->dphidxi[i][2]);

              /* (p,b) = (2,2) */
              f = -g2[2] * (fv->hq[2][2]);
              bfl->d_grad_phi_dmesh[i][2][2][j] =
                  (f * phi_m[j] * bfl->B[2][0] + g[2] * bfl->dB[2][0][2][j]) *
                  (bfl->dphidxi[i][0]);

              bfl->d_grad_phi_dmesh[i][2][2][j] +=
                  (f * phi_m[j] * bfl->B[2][1] + g[2] * bfl->dB[2][1][2][j]) *
                  (bfl->dphidxi[i][1]);

              bfl->d_grad_phi_dmesh[i][2][2][j] +=
                  (f * phi_m[j] * bfl->B[2][2] + g[2] * bfl->dB[2][2][2][j]) *
                  (bfl->dphidxi[i][2]);

              /* (p,b) = (2,1) */
              f = -g2[2] * (fv->hq[2][1]);
              bfl->d_grad_phi_dmesh[i][2][1][j] =
                  (f * phi_m[j] * bfl->B[2][0] + g[2] * bfl->dB[2][0][1][j]) *
                  (bfl->dphidxi[i][0]);

              bfl->d_grad_phi_dmesh[i][2][1][j] +=
                  (f * phi_m[j] * bfl->B[2][1] + g[2] * bfl->dB[2][1][1][j]) *
                  (bfl->dphidxi[i][1]);

              bfl->d_grad_phi_dmesh[i][2][1][j] +=
                  (f * phi_m[j] * bfl->B[2][2] + g[2] * bfl->dB[2][2][1][j]) *
                  (bfl->dphidxi[i][2]);
              /* (p,b) = (2,0) */
              f = -g2[2] * (fv->hq[2][0]);
              bfl->d_grad_phi_dmesh[i][2][0][j] =
                  (f * phi_m[j] * bfl->B[2][0] + g[2] * bfl->dB[2][0][0][j]) *
                  (bfl->dphidxi[i][0]);

              bfl->d_grad_phi_dmesh[i][2][0][j] +=
                  (f * phi_m[j] * bfl->B[2][1] + g[2] * bfl->dB[2][1][0][j]) *
                  (bfl->dphidxi[i][1]);

              bfl->d_grad_phi_dmesh[i][2][0][j] +=
                  (f * phi_m[j] * bfl->B[2][2] + g[2] * bfl->dB[2][2][0][j]) *
                  (bfl->dphidxi[i][2]);

              /* (p,b) = (1,2) */
              f = -g2[1] * (fv->hq[1][2]);
              bfl->d_grad_phi_dmesh[i][1][2][j] =
                  (f * phi_m[j] * bfl->B[1][0] + g[1] * bfl->dB[1][0][2][j]) *
                  (bfl->dphidxi[i][0]);

              bfl->d_grad_phi_dmesh[i][1][2][j] +=
                  (f * phi_m[j] * bfl->B[1][1] + g[1] * bfl->dB[1][1][2][j]) *
                  (bfl->dphidxi[i][1]);

              bfl->d_grad_phi_dmesh[i][1][2][j] +=
                  (f * phi_m[j] * bfl->B[1][2] + g[1] * bfl->dB[1][2][2][j]) *
                  (bfl->dphidxi[i][2]);

              /* (p,b) = (0,2) */
              f = -g2[0] * (fv->hq[0][2]);
              bfl->d_grad_phi_dmesh[i][0][2][j] =
                  (f * phi_m[j] * bfl->B[0][0] + g[0] * bfl->dB[0][0][2][j]) *
                  (bfl->dphidxi[i][0]);

              bfl->d_grad_phi_dmesh[i][0][2][j] +=
                  (f * phi_m[j] * bfl->B[0][1] + g[0] * bfl->dB[0][1][2][j]) *
                  (bfl->dphidxi[i][1]);

              bfl->d_grad_phi_dmesh[i][0][2][j] +=
                  (f * phi_m[j] * bfl->B[0][2] + g[0] * bfl->dB[0][2][2][j]) *
                  (bfl->dphidxi[i][2]);
            }
          }
        }
#endif

        /*
         * Evaluate the sensitivity of gradients of the fundamental vector
         * quantity (the basisfunction, phi_i times unit vectors e_a) in
         * this coordinate system with respect to mesh displacement unknowns
         *
         * i.e.,
         *		d ( grad_phi_e[i][a] [p][q])
         *		--------------------
         *		d ( d_b,j )
         */

        /* initialize variables */

        siz = sizeof(double) * DIM * DIM * DIM * DIM * MDE * MDE;
        memset(&(bfl->d_grad_phi_e_dmesh[0][0][0][0][0][0]), 0, siz);

        /* for Cartesian coordinates we have a nice vanilla derivative
         *   -> Note, that not all entries are filled in below. We rely on
         * initial zeroing to zero the entries for cartesian coordinates (or the
         * memset above)
         */
#ifdef DO_NO_UNROLL
        for (i = 0; i < vdofs; i++) {
          for (a = 0; a < wim; a++) {
            for (p = 0; p < dim; p++) {
              for (b = 0; b < dim; b++) {
                for (j = 0; j < mdofs; j++) {
                  bfl->d_grad_phi_e_dmesh[i][a][p][a][b][j] =
                      bfl->d_grad_phi_dmesh[i][p][b][j];
                }
              }
            }
          }
        }
#else
        for (i = 0; i < vdofs; i++) {
          for (j = 0; j < mdofs; j++) {
            /*	  bfl->d_grad_phi_e_dmesh[i][0] [p][0] [b][j] =
               bfl->d_grad_phi_dmesh[i][p] [b][j];
                  bfl->d_grad_phi_e_dmesh[i][1] [p][1] [b][j] =
               bfl->d_grad_phi_dmesh[i][p] [b][j];*/

            /* (p,b) = (0,0) */
            bfl->d_grad_phi_e_dmesh[i][0][0][0][0][j] =
                bfl->d_grad_phi_dmesh[i][0][0][j];
            bfl->d_grad_phi_e_dmesh[i][1][0][1][0][j] =
                bfl->d_grad_phi_dmesh[i][0][0][j];

            /* (p,b) = (1,1)*/
            bfl->d_grad_phi_e_dmesh[i][0][1][0][1][j] =
                bfl->d_grad_phi_dmesh[i][1][1][j];
            bfl->d_grad_phi_e_dmesh[i][1][1][1][1][j] =
                bfl->d_grad_phi_dmesh[i][1][1][j];

            /* (p,b) = (0,1)*/
            bfl->d_grad_phi_e_dmesh[i][0][0][0][1][j] =
                bfl->d_grad_phi_dmesh[i][0][1][j];
            bfl->d_grad_phi_e_dmesh[i][1][0][1][1][j] =
                bfl->d_grad_phi_dmesh[i][0][1][j];

            /* (p,b) = (1,0)*/
            bfl->d_grad_phi_e_dmesh[i][0][1][0][0][j] =
                bfl->d_grad_phi_dmesh[i][1][0][j];
            bfl->d_grad_phi_e_dmesh[i][1][1][1][0][j] =
                bfl->d_grad_phi_dmesh[i][1][0][j];

            if (wim == 3) {
              /*  bfl->d_grad_phi_e_dmesh[i][2] [p][2] [b][j] =
               * bfl->d_grad_phi_dmesh[i][p] [b][j];*/
              /* (p,b) = (0,0) */
              bfl->d_grad_phi_e_dmesh[i][2][0][2][0][j] =
                  bfl->d_grad_phi_dmesh[i][0][0][j];
              /* (p,b) = (1,1)*/
              bfl->d_grad_phi_e_dmesh[i][2][1][2][1][j] =
                  bfl->d_grad_phi_dmesh[i][1][1][j];
              /* (p,b) = (0,1)*/
              bfl->d_grad_phi_e_dmesh[i][2][0][2][1][j] =
                  bfl->d_grad_phi_dmesh[i][0][1][j];
              /* (p,b) = (1,0)*/
              bfl->d_grad_phi_e_dmesh[i][2][1][2][0][j] =
                  bfl->d_grad_phi_dmesh[i][1][0][j];

              if (dimNonSym == 3) {
                /* (p,b) = (2,2)*/
                bfl->d_grad_phi_e_dmesh[i][0][2][0][2][j] =
                    bfl->d_grad_phi_dmesh[i][2][2][j];
                bfl->d_grad_phi_e_dmesh[i][1][2][1][2][j] =
                    bfl->d_grad_phi_dmesh[i][2][2][j];
                bfl->d_grad_phi_e_dmesh[i][2][2][2][2][j] =
                    bfl->d_grad_phi_dmesh[i][2][2][j];
                /* (p,b) = (2,0)*/
                bfl->d_grad_phi_e_dmesh[i][0][2][0][0][j] =
                    bfl->d_grad_phi_dmesh[i][2][0][j];
                bfl->d_grad_phi_e_dmesh[i][1][2][1][0][j] =
                    bfl->d_grad_phi_dmesh[i][2][0][j];
                bfl->d_grad_phi_e_dmesh[i][2][2][2][0][j] =
                    bfl->d_grad_phi_dmesh[i][2][0][j];
                /* (p,b) = (2,1)*/
                bfl->d_grad_phi_e_dmesh[i][0][2][0][1][j] =
                    bfl->d_grad_phi_dmesh[i][2][1][j];
                bfl->d_grad_phi_e_dmesh[i][1][2][1][1][j] =
                    bfl->d_grad_phi_dmesh[i][2][1][j];
                bfl->d_grad_phi_e_dmesh[i][2][2][2][1][j] =
                    bfl->d_grad_phi_dmesh[i][2][1][j];
                /* (p,b) = (1,2)*/
                bfl->d_grad_phi_e_dmesh[i][0][1][0][2][j] =
                    bfl->d_grad_phi_dmesh[i][1][2][j];
                bfl->d_grad_phi_e_dmesh[i][1][1][1][2][j] =
                    bfl->d_grad_phi_dmesh[i][1][2][j];
                bfl->d_grad_phi_e_dmesh[i][2][1][2][2][j] =
                    bfl->d_grad_phi_dmesh[i][1][2][j];
                /* (p,b) = (0,2)*/
                bfl->d_grad_phi_e_dmesh[i][0][0][0][2][j] =
                    bfl->d_grad_phi_dmesh[i][0][2][j];
                bfl->d_grad_phi_e_dmesh[i][1][0][1][2][j] =
                    bfl->d_grad_phi_dmesh[i][0][2][j];
                bfl->d_grad_phi_e_dmesh[i][2][0][2][2][j] =
                    bfl->d_grad_phi_dmesh[i][0][2][j];
              }
            }
          }
        }

#endif

        /* add more involved pieces as necessary */
        if (pd->CoordinateSystem != CARTESIAN) {

          /*  for ( i=0; i<vdofs; i++) {
              for ( a=0; a<dim; a++) {
              for ( p=0; p<dim; p++) {
              for ( b=0; b<dim; b++)  {
              for ( j=0; j<mdofs; j++) {
              bfl->d_grad_phi_e_dmesh[i][a] [p][a] [b][j]
              = bfl->d_grad_phi_dmesh[i][p] [b][j];
              }
              }
              }
              }
              } */

          for (i = 0; i < vdofs; i++) {
            for (a = 0; a < wim; a++) {
              for (p = 0; p < VIM; p++) {
                for (q = 0; q < VIM; q++) {
                  for (b = 0; b < dim; b++) {
                    for (j = 0; j < mdofs; j++) {
                      if (q != a) {
                        if (bfl->d_grad_phi_e_dmesh[i][a][p][q][b][j] != 0.0) {
                          printf("we shouldn't be here\n");
                          exit(-1);
                        }
                        bfl->d_grad_phi_e_dmesh[i][a][p][q][b][j] = 0.0;
                      }
                      if (dim < VIM && (p == VIM || q == VIM)) {
                        if (bfl->d_grad_phi_e_dmesh[i][a][p][q][b][j] != 0.0) {
                          printf("we shouldn't be here\n");
                          exit(-1);
                        }
                        bfl->d_grad_phi_e_dmesh[i][a][p][q][b][j] = 0.0;
                      }

                      bfl->d_grad_phi_e_dmesh[i][a][p][q][b][j] +=
                          phi_l[i] * fv->d_grad_e_dq[a][p][q][b] * phi_m[j];
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
  return 0;

} /* END of load_bf_mesh_derivs */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int load_basis_functions(const double xi[], /*  [DIM]               */
                         struct Basis_Functions **bfa) /* ptr to basis function
                                                        * * array of interest */

/************************************************************************
 *
 * load_basis_functions():
 *
 *    Calculates the values of all the basis functions active in the
 * current element. It also calculates the derivatives of the basis
 * functions wrt local element coordinates.
 *
 ************************************************************************/
{
  int b, i, v, jdof, ledof;
  int mn = ei[pg->imtrx]->mn;

  int imtrx;

  BASIS_FUNCTIONS_STRUCT *bf_ptr;
  /*
   * Load basis functions and derivatives in the unit elements for each
   * kind of unique basis function that we have...
   */
  for (b = 0; b < Num_Basis_Functions; b++) {
    for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
      bf_ptr = bfa[b];

      /*
       * Look up the variable "v" that has the the right interpolation
       * in the current material as the current basis function
       */

      v = bf_ptr->Var_Type_MatID[mn];

      if (!pd->v[imtrx][v])
        continue;

      /*
       * don't calculate basis function if interpolation doesn't
       *  correspond to any variable in the current material
       *  OR if element shape doesn't match the current element block.
       */
      if (v != -1 && bf_ptr->element_shape == ei[imtrx]->ielem_shape) {
        /*
         * Now, case the dimensionality and look up basis functions
         * and their derivatives at the quadrature point
         * using elemental coordinates
         */
        jdof = 0;
        switch (pd->Num_Dim) {
        case 1:
          for (i = 0; i < ei[imtrx]->dof[v]; i++) {
            ledof = ei[imtrx]->lvdof_to_ledof[v][i];
            if (ei[imtrx]->active_interp_ledof[ledof]) {
              bf_ptr->phi[i] = newshape(
                  xi, ei[imtrx]->ielem_type, PSI, ei[imtrx]->dof_list[v][i],
                  bf_ptr->element_shape, bf_ptr->interpolation, jdof);
              bf_ptr->dphidxi[i][0] = newshape(
                  xi, ei[imtrx]->ielem_type, DPSI_S, ei[imtrx]->dof_list[v][i],
                  bf_ptr->element_shape, bf_ptr->interpolation, jdof);
              jdof++;
            } else {
              bf_ptr->phi[i] = 0.0;
              bf_ptr->dphidxi[i][0] = 0.0;
            }
          }
          break;

        case 2:
          for (i = 0; i < ei[imtrx]->dof[v]; i++) {
            ledof = ei[imtrx]->lvdof_to_ledof[v][i];
            if (ei[imtrx]->active_interp_ledof[ledof]) {
              bf_ptr->phi[i] = newshape(
                  xi, ei[imtrx]->ielem_type, PSI, ei[imtrx]->dof_list[v][i],
                  bf_ptr->element_shape, bf_ptr->interpolation, jdof);
              bf_ptr->dphidxi[i][0] = newshape(
                  xi, ei[imtrx]->ielem_type, DPSI_S, ei[imtrx]->dof_list[v][i],
                  bf_ptr->element_shape, bf_ptr->interpolation, jdof);
              bf_ptr->dphidxi[i][1] = newshape(
                  xi, ei[imtrx]->ielem_type, DPSI_T, ei[imtrx]->dof_list[v][i],
                  bf_ptr->element_shape, bf_ptr->interpolation, jdof);
              jdof++;
            } else {
              bf_ptr->phi[i] = 0.0;
              bf_ptr->dphidxi[i][0] = 0.0;
              bf_ptr->dphidxi[i][1] = 0.0;
            }
          }
          break;

        case 3:
          for (i = 0; i < ei[imtrx]->dof[v]; i++) {
            ledof = ei[imtrx]->lvdof_to_ledof[v][i];
            if (ei[imtrx]->active_interp_ledof[ledof]) {
              bf_ptr->phi[i] = newshape(
                  xi, ei[imtrx]->ielem_type, PSI, ei[imtrx]->dof_list[v][i],
                  bf_ptr->element_shape, bf_ptr->interpolation, jdof);
              bf_ptr->dphidxi[i][0] = newshape(
                  xi, ei[imtrx]->ielem_type, DPSI_S, ei[imtrx]->dof_list[v][i],
                  bf_ptr->element_shape, bf_ptr->interpolation, jdof);
              bf_ptr->dphidxi[i][1] = newshape(
                  xi, ei[imtrx]->ielem_type, DPSI_T, ei[imtrx]->dof_list[v][i],
                  bf_ptr->element_shape, bf_ptr->interpolation, jdof);
              bf_ptr->dphidxi[i][2] = newshape(
                  xi, ei[imtrx]->ielem_type, DPSI_U, ei[imtrx]->dof_list[v][i],
                  bf_ptr->element_shape, bf_ptr->interpolation, jdof);
              jdof++;
            } else {
              bf_ptr->phi[i] = 0.0;
              bf_ptr->dphidxi[i][0] = 0.0;
              bf_ptr->dphidxi[i][1] = 0.0;
              bf_ptr->dphidxi[i][2] = 0.0;
            }
          }
          break;
        }
      }
    }
  }
  return (0);
} /* END of routine load_basis_functions */
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/* asdv() -- allocate, set a double vector
 *
 * Description: a pointer to a pointer to double is input. Space for
 *              a vector of doubles is allocated as a one dimensional array
 *              and each element of the array is set to zero.
 *
 * Return: void
 *
 * Created: 1997/01/29 09:07 MST pasacki@sandia.gov
 */

void asdv(double **v,  /* vector to be allocated */
          const int n) /* number of elements in vector */
{
  int i;
  double *p;
  const double zero = (double)0;

  *v = (double *)smalloc(n * sizeof(double));

  for (i = 0, p = *v; i < n; i++, p++) {
    *p = zero;
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void alloc_extern_ija_buffer(const int n,    /* Order full system    (in) */
                             const int m,    /* Order trim system    (in) */
                             const int *ija, /* Orig column pointers (in) */
                             int **ptr_ija_attic) /* where to hide (out) */

/* alloc_extern_ija_buffer() -- get space for external rows hiding from aztec
 *
 * Created: 1997/11/03 12:47 MST pasacki@sandia.gov
 *
 * Revised:
 */
{
  int size;
  if (n == m)
    return;

  /*
   * How big must it be? Big enough to save all the diagonal entries and
   * the column names for all of the external rows past m-1.
   */
  size = n + 1 + (ija[n] - ija[m]);
  *ptr_ija_attic = alloc_int_1(size, 0);
  return;
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void alloc_MSR_sparse_arrays(
    int **ptr_ija,      /* column pointer array */
    double **ptr_a,     /* values of nonzero matrix entries */
    double **ptr_a_old, /* backup array of matrix nonzeros */
    const int Fill,     /* flag to allocate spaces for either
                         * the totally coupled problem or just
                         * explicit nodal fill eqn */
    int node_to_fill[], /* node index gives fill dof */
    Exo_DB *exo,        /* ptr to the whole mesh */
    Dpi *dpi)           /* distributed processing info */

/******************************************************************************

Function which allocates the sparse matrix storage vectors as well as other
vectors associated with the problem solution.


Author:          Scott Hutchinson (1421)
Date:            11 January 1993
Revised:         11/19/96 RRR
Revised:         1997/10/15 08:11 MDT pasacki@sandia.gov
Revised:	   1997/10/28 16:48 MST pasacki@sandia.gov

*******************************************************************************/
{
  int itotal_nodes;
  int max_neigh_elem = -1;

#ifdef PARALLEL
  int local_value;
#endif

  int num_total_fill_unknowns = 0;
  int num_total_unknowns;

  /*
   * I prefer to think of ija and a as specific types of variables and
   * changed the definition in the argument list to more accurately reflect
   * what is really going on. The old notation is unnecessarily obtuse.
   */

  int *ija;
  double *a;
  double *a_old;

  int *ija_temp; /* temporary vector used in determining ija and
                    its size                                     */
  int nnz;       /* number of nonzero entries in 'a'             */

  itotal_nodes = dpi->num_universe_nodes; /* include internal, boundary and
                                           * external nodes that are in the
                                           * universe of this processor */

  if (Num_Dim == 1) {
    max_neigh_elem = 2;
  } else if (Num_Dim == 2) {
    max_neigh_elem = MAX_SUR_ELEM_2D; /* 20 deg angles give 18 elems... */
  } else if (Num_Dim == 3) {
    max_neigh_elem = MAX_SUR_ELEM_3D; /* 20 deg solid angle gives 60 elems */
  }

  if (Fill) {
    nnz = find_problem_graph_fill(&ija_temp, itotal_nodes, max_neigh_elem,
                                  node_to_fill, exo);
  } else {
    nnz = find_MSR_problem_graph(&ija_temp, itotal_nodes, exo);
  }

#ifndef DEBUG_HKM
  if (Debug_Flag)
#endif
  {
    printf("Processor %d nnz = %d\n", ProcID, nnz);
  }

  /* Allocate and fill the true ija. */

  /*
   * It looks like ija is allocated twice.  The second allocation
   * means the first is lost -- a huge memory leak.  I'll comment the
   * first allocation.
   */
  /* ija = (int *) array_alloc(1, nnz, sizeof(int)); */

  ija = alloc_int_1(nnz, INT_NOINIT);
  memcpy(ija, ija_temp, nnz * sizeof(int));

#ifdef DEBUG_IJA
  /*
   * Dump out ija[]...
   */
  {
    FILE *fff;
    fff = fopen("ija_dump2.txt", "w");
    for (i = 0; i < nnz; i++) {
      fprintf(fff, "P_%d: ija[%d] = %d\n", ProcID, i, ija[i]);
    }
    fclose(fff);
  }
#endif

  /* free temporary arrays */
  safer_free((void **)&ija_temp);

  /*
   * Allocate the nonzero storage vector and initialize it to zero
   */
  a = alloc_dbl_1(nnz, 0.0);

  /*   num_fill_unknowns  = dpi->num_universe_nodes; */

  num_total_unknowns = num_universe_dofs[pg->imtrx];

  /* set boolean for save old A matrix for data sensitivities,
     flux sensitivities, modified newton schemes, augmenting conditions,
     and continuation */

  save_old_A = FALSE;
  if (modified_newton || nAC > 0 || nn_post_data_sens > 0 ||
      nn_post_fluxes_sens > 0 || Continuation > 0) {
    save_old_A = TRUE;
  }

  if (save_old_A) {
    a_old = alloc_dbl_1(nnz, 0.0);
  } else {
    a_old = NULL;
  }

#ifdef PARALLEL
  /*
   *  Calculate the total number of fill unknowns and regular unknowns in
   *  the global problem by summing up contributions from each
   *  processor.
   */

  num_total_fill_unknowns =
      gsum_Int(internal_fill_unknowns + boundary_fill_unknowns);
  num_total_unknowns =
      gsum_Int(num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx]);
  dpi->num_dofs_global = num_total_unknowns;

  /*
   * The number of nonzero matrix entries are only those owned by this
   * processor, so let's not count any belonging to external nodes or dofs.
   */

  if (Fill) {
    local_value = (ija[internal_fill_unknowns + boundary_fill_unknowns] -
                   external_fill_unknowns - 1);
  } else {
    local_value =
        (ija[num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx]] -
         num_external_dofs[pg->imtrx] - 1);
  }
  nnz = gsum_Int(local_value);
#else
  dpi->num_dofs_global = num_total_unknowns;
  num_total_fill_unknowns = num_fill_unknowns;
#endif

  if (Fill) {
    DPRINTF(stdout, "\n%-30s= %d\n", "Number of fill unknowns",
            num_total_fill_unknowns);
    DPRINTF(stdout, "\n%-30s= %d\n", "Number of filmatrix nonzeroes", nnz);
  } else {
    DPRINTF(stdout, "\n%-30s= %d\n", "Number of unknowns", num_total_unknowns);
    DPRINTF(stdout, "\n%-30s= %d\n", "Number of matrix nonzeroes", nnz);
  }

  *ptr_ija = ija;
  *ptr_a = a;
  *ptr_a_old = a_old;

  return;
}
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

void alloc_VBR_sparse_arrays(struct Aztec_Linear_Solver_System *ams,
                             Exo_DB *exo, Dpi *dpi)
/*
 *
 */
{
  int nnz, nnz_tmp;
  int row_nodes, col_nodes;

  row_nodes = dpi->num_internal_nodes + dpi->num_boundary_nodes;
  col_nodes = dpi->num_universe_nodes;

  nnz = find_VBR_problem_graph(&(ams->indx), &(ams->bindx), &(ams->bpntr),
                               &(ams->rpntr), &(ams->cpntr), row_nodes,
                               col_nodes, exo);

  if (Debug_Flag) {
    printf("Processor %d nnz = %d\n", ProcID, nnz);
  }

  ams->val = alloc_dbl_1(nnz + 1, 0.0);
  /* Why the "+ 1" ?  Later on just one more entry is need
   * for compatibility with the way MSR worked.  Actually,
   * it is related to the idiotic ija_extern kludge */

  ams->nnz = nnz;
  ams->nnz_plus = nnz;

  ams->npn = row_nodes;
  ams->npn_plus = col_nodes;

  ams->npu = ams->rpntr[row_nodes];
  ams->npu_plus = ams->cpntr[col_nodes];

  /* set boolean for save old A matrix for data sensitivities,
     flux sensitivities, modified newton schemes, augmenting conditions,
     and continuation */

  save_old_A = FALSE;
  if (modified_newton || nAC > 0 || nn_post_data_sens > 0 ||
      nn_post_fluxes_sens > 0 || Continuation > 0) {
    save_old_A = TRUE;
  }

  if (save_old_A) {
    asdv(&(ams->val_old), nnz);
  }

  /*
   *   calculate the total number of degrees of freedom in the problem
   *   by summing the contributions from each processor
   */
  dpi->num_dofs_global =
      gsum_Int(num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx]);
  nnz_tmp = gsum_Int(nnz);

  DPRINTF(stdout, "\n%-30s= %d\n", "Number of unknowns", dpi->num_dofs_global);
  DPRINTF(stdout, "\n%-30s= %d\n", "Number of matrix nonzeroes", nnz_tmp);
  return;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int
find_problem_graph_fill(int *ija[],         /* column pointer array */
                        int itotal_nodes,   /* number of nodes this processor
                                               is resposible for */
                        int max_neigh_elem, /* maximum number of elements
                                               containing a given node */
                        int node_to_fill[], Exo_DB *exo)

/*
 *  Function which determines the graph of the sparse problem using the MSR
 *  format.  It also returns the exact number of nonzeroes in the problem.
 *  This routine is written for an unknown map of 1 dof per node.
 *
 *  ---------------------------------------------------------------------------
 *
 */
{
  int j, itmp, inode, ielem, iunknown;
  int inter_unknown, inter_node_unknowns, inter_node, inum_unknowns;
  int nz_ptr = 0;
  int irow_index = 0;
  int icol_index;
  int nz_temp = 0;
  /*
   * Surrounding element variables...
   */
  int *nsur_elem; /* array telling how many elements surround */
  /* each and every node in the mesh...*/
  int **isur_elem; /* array telling the names of the elements */
  /* that surround each and every element... */

  /*
   * Allocate arrays and initialize...
   * start allocating ija to a small estimate - reallocate later
   */
  nz_temp = num_fill_unknowns * Max_NP_Elem;
  *ija = alloc_int_1(nz_temp, -1);
#ifdef DEBUG
  printf("Initial ija allocation to %d\n", nz_temp);
#endif
  nsur_elem = alloc_int_1(itotal_nodes, 0);
  isur_elem = alloc_int_2(itotal_nodes, max_neigh_elem, -1);

  /*
   * Find surrounding elements for each node, that is, the elements which
   * contain each node.  Keep track of the names of elements which
   * contain each node in isur_elem[][] and the number of elements that
   * contain each node in nsur_elem[]...
   */

  for (ielem = 0; ielem < Num_Internal_Elems; ielem++) {
    for (j = 0; j < elem_info(NNODES, Elem_Type(exo, ielem)); j++) {
      inode = Proc_Elem_Connect[Proc_Connect_Ptr[ielem] + j];
      if (inode < itotal_nodes) {
        /*
         * internal or border node -> fill
         */

        isur_elem[inode][nsur_elem[inode]++] = ielem;

        if (nsur_elem[inode] >= max_neigh_elem) {
          fprintf(stderr, "nsur_elem exceed maximum of %d\n", max_neigh_elem);
          exit(-1);
        }
      }
    }
  }

  /* set the number of unknowns per node to one */

  nz_ptr = num_fill_unknowns + 1;

  for (inode = 0; inode < itotal_nodes; inode++) {
    /*
     * Loop over the unknowns associated with this node,
     * which is one for this case.
     */

    inum_unknowns = num_varType_at_node(inode, FILL);

    for (iunknown = 0; iunknown < inum_unknowns; iunknown++) {
      (*ija)[irow_index] = nz_ptr;

      /*
       * Loop over the nodes which are determined to have an interaction
       * with the current row node
       */
      for (j = exo->node_node_pntr[inode]; j < exo->node_node_pntr[inode + 1];
           j++) {
        inter_node = exo->node_node_list[j];

        /*
         * Loop over the unknowns associated with the
         * interacting node and see if there should be
         * an interaction.
         *
         */

        inter_node_unknowns = num_varType_at_node(inter_node, FILL);

        for (inter_unknown = 0; inter_unknown < inter_node_unknowns;
             inter_unknown++) {
          icol_index = node_to_fill[inter_node] + inter_unknown;

          if (icol_index > -1) {
            if (in_list(icol_index, (*ija)[irow_index], nz_ptr, *ija) < 0 &&
                icol_index != irow_index) {
              if (nz_ptr + 1 == nz_temp) {
                /* reallocate ija larger!! */
                nz_temp += num_fill_unknowns * Max_NP_Elem;
#ifdef DEBUG
                printf("Reallocate ija to %d\n", nz_temp);
#endif

                *ija = (int *)realloc((void *)*ija, nz_temp * sizeof(int));
                if (*ija == NULL)
                  EH(-1, "No space for ija_temp");
                for (itmp = nz_ptr; itmp < nz_temp; (*ija)[itmp++] = -1)
                  ;
              }
              (*ija)[nz_ptr++] = icol_index;
            }
          }
        }
      }
      irow_index++;
    }
  }

  (*ija)[irow_index] = nz_ptr;

  safer_free((void **)&nsur_elem);
  safer_free((void **)&isur_elem);
  return (nz_ptr);
} /* END of find_problem_graph_fill  */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

static int
find_MSR_problem_graph(int *ija[],       /* column pointer array            */
                       int itotal_nodes, /* number nodes this processor     */
                       Exo_DB *exo)      /* ptr to FE db                    */

{
  int j, itmp, inode, i1, i2, eb1;
  int iunknown, inter_unknown, inter_node, row_num_unknowns, col_num_unknowns;
  int nz_ptr, irow_index = 0;
  int icol_index, rowVarType, colVarType;
  int *ija_ptr;
  int add_var = 0, nz_temp;
  NODE_INFO_STRUCT *nodeCol;
  NODAL_VARS_STRUCT *nv, *nvCol;
  double avgVarPerNode, avgNodeToNodeConn;

  int *inode_varType, *inode_matID;
  int *inter_node_varType, *inter_node_matID;
#ifdef DEBUG_GRAPH
  char *yo = "find_MSR_problem_graph: ";
#endif
#ifdef DEBUG_GRAPH
  FILE *ifp = NULL;
  ifp = fopen("msr_format.txt", "w");
#endif /* DEBUG_GRAPH */

  /*
   * Calculate the average number of degrees of freedom per node
   */
  avgVarPerNode = (double)num_universe_dofs[pg->imtrx] / (double)itotal_nodes;
  /*
   * Estimate the average connectivity of a node
   */
  avgNodeToNodeConn =
      (double)exo->node_node_pntr[itotal_nodes] / (double)itotal_nodes;
  /*
   *  Start allocating ija to a slightly larger value than
   *  needed by a good estimate.
   *  -> Check for null pointer because this is a large
   *     allocation
   */
  nz_temp =
      num_universe_dofs[pg->imtrx] * avgNodeToNodeConn * avgVarPerNode * 1.1;

  *ija = ija_ptr = alloc_int_1(nz_temp, -1);
  if (*ija == NULL) {
    EH(-1, "No space for ija_temp");
  }

  inode_varType = alloc_int_1(MaxVarPerNode, INT_NOINIT);
  inode_matID = alloc_int_1(MaxVarPerNode, INT_NOINIT);
  inter_node_varType = alloc_int_1(MaxVarPerNode, INT_NOINIT);
  inter_node_matID = alloc_int_1(MaxVarPerNode, INT_NOINIT);

  /*
   * Initialize the pointer into the jacobian, nz_ptr,
   * corresponding to the
   * start of off-diagonal entries for a particular row.
   * Row zero will have its first entry at position
   * num_universe_dofs + 1, because the diagonal is storred first
   * and then the total number of entries is storred at
   * ija[num_universe_dofs].
   */
  nz_ptr = num_universe_dofs[pg->imtrx] + 1;
#ifdef DEBUG_GRAPH
  fprintf(ifp, "P_%d: nz_ptr = %d\n", ProcID, nz_ptr);
#endif

  /*
   * loop over all of the nodes on this processor
   */
  for (inode = 0; inode < itotal_nodes; inode++) {
    nv = Nodes[inode]->Nodal_Vars_Info[pg->imtrx];
    /*
     * Fill the vector list which points to the unknowns defined at this
     * node...
     */
    row_num_unknowns = fill_variable_vector(inode, inode_varType, inode_matID);
    /*
     * Do a check against the number of unknowns at this
     * node stored in the global array
     */
    if (row_num_unknowns != nv->Num_Unknowns) {
      EH(-1, "Inconsistency counting unknowns.");
    }

    /*
     * Loop over the unknowns defined at this row node
     */
    for (iunknown = 0; iunknown < row_num_unknowns; iunknown++) {
      /*
       * Retrieve the var type of the current unknown
       */
      rowVarType = inode_varType[iunknown];

      /*
       * Store the beginning location for off-diagonal jacobian
       * entries for the current unknown
       */
      ija_ptr[irow_index] = nz_ptr;
#ifdef DEBUG_GRAPH
      fprintf(ifp, "node = %d, unk = %d, irow_index = %d, nz_ptr = %d\n", inode,
              iunknown, irow_index, nz_ptr);
#endif
      /*
       * Loop over the nodes which are determined to have an interaction
       * with the current row node
       */
      for (j = exo->node_node_pntr[inode]; j < exo->node_node_pntr[inode + 1];
           j++) {
        inter_node = exo->node_node_list[j];
        nodeCol = Nodes[inter_node];
        nvCol = nodeCol->Nodal_Vars_Info[pg->imtrx];

        /*
         * fill the vector list which points to the unknowns
         * defined at this interaction node
         */
        col_num_unknowns = fill_variable_vector(inter_node, inter_node_varType,
                                                inter_node_matID);
        if (col_num_unknowns != nvCol->Num_Unknowns) {
          EH(-1, "Inconsistency counting unknowns.");
        }

        /*
         * Loop over the unknowns associated with the
         * interacting node and see if there should be an interaction.
         */

        for (inter_unknown = 0; inter_unknown < col_num_unknowns;
             inter_unknown++) {

          /*
           * HKM ->
           *  The entire process below is designed to find out
           *  what unknown we are processing. This coding is very
           *  convoluted and fraught with pitfalls. It should
           *  be rewritten using the structured approach in MPSalsa
           *  to variable indentification.
           */
          colVarType = inter_node_varType[inter_unknown];

          /*
           * Query the Interaction mask to determine if a jacobian entry
           * should be created
           */
          add_var = Inter_Mask[pg->imtrx][rowVarType][colVarType];

          /* The following code should be activated when solving DG viscoelastic
           * problems with full Jacobian treatment of upwind element stress
           * terms
           */
          if (exo->centroid_list[inode] != -1 && inode != inter_node &&
              exo->centroid_list[inter_node] != -1) {
            eb1 = exo->elem_eb[exo->centroid_list[inode]];

            if (vn_glob[Matilda[eb1]]->dg_J_model == FULL_DG) {
              i1 = pd_glob[Matilda[eb1]]->i[pg->imtrx][rowVarType];
              i2 = pd_glob[Matilda[eb1]]->i[pg->imtrx][colVarType];

              if ((rowVarType == colVarType) &&
                  (i1 == I_P0 || i1 == I_P1 || i1 == I_PQ1 || i1 == I_PQ2) &&
                  (i2 == I_P0 || i2 == I_P1 || i2 == I_PQ1 || i2 == I_PQ2) &&
                  (rowVarType != PRESSURE) &&
                  (rowVarType > VELOCITY_GRADIENT33 ||
                   rowVarType < VELOCITY_GRADIENT11)) {
                add_var = Inter_Mask[pg->imtrx][rowVarType][colVarType];
              } else {
                add_var = 0;
              }
            }
          }

          if (Debug_Flag < 0)
            add_var = TRUE; /* add all vars for checking jacobian */

          if (add_var) {
            /*
             * Determine the equation number of the current unknown
             */
            icol_index = nodeCol->First_Unknown[pg->imtrx] + inter_unknown;

            if (icol_index != irow_index) {
              /*
               * Check to See if we need to increase storage
               */
              if (nz_ptr + 1 == nz_temp) {

                /* reallocate ija larger */
                nz_temp += (int)(2.0 * (itotal_nodes - inode) *
                                 ((double)nz_ptr / (double)inode));
#ifdef DEBUG_GRAPH
                printf("%s Reallocate ija size to %d from %d\n", yo, nz_temp,
                       nz_ptr + 1);
#endif
                *ija = ija_ptr =
                    (int *)realloc((void *)*ija, nz_temp * sizeof(int));
                if (*ija == NULL)
                  EH(-1, "No space for ija_temp");
                for (itmp = nz_ptr; itmp < nz_temp; ija_ptr[itmp++] = -1)
                  ;
              }
              ija_ptr[nz_ptr++] = icol_index;
            }
          }
        }
      }
      irow_index++;
    }
  }

  ija_ptr[irow_index] = nz_ptr;

#ifdef DEBUG_GRAPH
  printf("%s Final size of ija is %d\n", yo, nz_ptr);
#endif

  safer_free((void **)&inode_varType);
  safer_free((void **)&inode_matID);
  safer_free((void **)&inter_node_varType);
  safer_free((void **)&inter_node_matID);

#ifdef DEBUG_GRAPH
  fclose(ifp);
#endif
  return (nz_ptr);
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

static int find_VBR_problem_graph(int *indx[], int *bindx[], int *bpntr[],
                                  int *rpntr[], int *cpntr[], int row_nodes,
                                  int col_nodes, Exo_DB *exo)

/*
 *
 */
{
  int i, j, k, nnz_blocs = 0; /* number of non-zero blocks */
  int nnz = 0, row_ptr = 0, col_ptr = 0;
  NODAL_VARS_STRUCT *nvRow, *nvCol;
#ifdef DEBUG
  char filename[80];
  FILE *of;
#endif

  *bpntr = alloc_int_1(row_nodes + 1, INT_NOINIT);
  *rpntr = alloc_int_1(row_nodes + 1, INT_NOINIT);
  *cpntr = alloc_int_1(col_nodes + 1, INT_NOINIT);

  for (i = 0; i < row_nodes + 1; i++) {
    (*bpntr)[i] = exo->node_node_pntr[i];
  }

  nnz_blocs = (*bpntr)[row_nodes];

  *bindx = alloc_int_1(nnz_blocs, INT_NOINIT);
  *indx = alloc_int_1(nnz_blocs + 1, INT_NOINIT);

  for (i = 0; i < nnz_blocs; i++) {
    (*bindx)[i] = exo->node_node_list[i];
  }

  for (i = 0; i < row_nodes; i++) {
    nvRow = Nodes[i]->Nodal_Vars_Info[pg->imtrx];
    for (k = (*bpntr)[i]; k < (*bpntr)[i + 1]; k++) {
      j = (*bindx)[k];
      nvCol = Nodes[j]->Nodal_Vars_Info[pg->imtrx];
      (*indx)[k] = nnz;
      nnz += nvRow->Num_Unknowns * nvCol->Num_Unknowns;
    }
  }

  (*indx)[nnz_blocs] = nnz;

  for (i = 0; i < row_nodes; i++) {
    nvRow = Nodes[i]->Nodal_Vars_Info[pg->imtrx];
    (*rpntr)[i] = row_ptr;
    row_ptr += nvRow->Num_Unknowns;
  }
  (*rpntr)[row_nodes] = row_ptr;

  for (j = 0; j < col_nodes; j++) {
    nvCol = Nodes[j]->Nodal_Vars_Info[pg->imtrx];
    (*cpntr)[j] = col_ptr;
    col_ptr += nvCol->Num_Unknowns;
  }
  (*cpntr)[col_nodes] = col_ptr;

#ifdef DEBUG

  sprintf(filename, "vbr%d_of_%d", ProcID + 1, Num_Proc);
  of = fopen(filename, "w");

  fprintf(of, "P_%d nnz=%d\n\n", ProcID, nnz);
  fprintf(of, "P_%d nnz_blocs=%d\n\n", ProcID, nnz_blocs);
  fprintf(of, "P_%d row_nodes=%d\n\n", ProcID, row_nodes);
  fprintf(of, "P_%d col_nodes=%d\n\n", ProcID, col_nodes);

  for (i = 0; i < nnz_blocs; i++) {
    fprintf(of, "P_%d bindx[%d]=%d\n", ProcID, i, exo->node_node_list[i]);
  }
  fprintf(of, "\n");

  for (i = 0; i < nnz_blocs + 1; i++) {
    fprintf(of, "P_%d indx[%d]=%d\n", ProcID, i, (*indx)[i]);
  }

  fprintf(of, "\n");
  for (i = 0; i < row_nodes + 1; i++) {
    fprintf(of, "P_%d bpntr[%d]=%d\n", ProcID, i, (*bpntr)[i]);
  }

  fprintf(of, "\n");
  for (i = 0; i < row_nodes + 1; i++) {
    fprintf(of, "P_%d rpntr[%d]=%d\n", ProcID, i, (*rpntr)[i]);
  }

  fprintf(of, "\n");
  for (j = 0; j < col_nodes + 1; j++) {
    fprintf(of, "P_%d cpntr[%d]=%d\n", ProcID, j, (*cpntr)[j]);
  }
#endif

  return (nnz);
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int fill_variable_vector(int inode, int ivec_varType[], int ivec_matID[])

/************************************************************************
 *
 * fill_variable_vector():
 *
 * Function which fills vectors of variable types and matid for the
 * variables defined at a given node.
 *
 ************************************************************************/
{
  int l, m, index = 0;
  NODE_INFO_STRUCT *node;
  NODAL_VARS_STRUCT *nvs;
  VARIABLE_DESCRIPTION_STRUCT *vd;
  node = Nodes[inode];
  nvs = node->Nodal_Vars_Info[pg->imtrx];
  for (l = 0; l < nvs->Num_Var_Desc; l++) {
    vd = nvs->Var_Desc_List[l];
    for (m = 0; m < vd->Ndof; m++) {
      ivec_varType[index] = vd->Variable_Type;
      ivec_matID[index] = vd->MatID;
      index++;
    }
  }
  if (index > MaxVarPerNode) {
    EH(-1, "fill_variable_vector: index is greater than MaxVarPerNode");
  }
  return index;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void zero_lec_row(
    double local_J[MAX_PROB_VAR + MAX_CONC][MAX_PROB_VAR + MAX_CONC][MDE][MDE],
    int eqn_type, /* Eqn Type of row to be zeroed     */
    int ldof)     /* Local dof of that equation type  */

/*
  Function which zero's the row of the local element stiffness matrix
  corresponding to equation of type eqn_type and the local dof number ldof

*/

{

  /* Local variables */

  int i_var; /* Row index into the global stiffness matrix   */

  /***************************** EXECUTION BEGINS
   * *******************************/

  for (i_var = 0; i_var < MAX_PROB_EQN + MAX_CONC; i_var++) {
    memset(local_J[eqn_type][i_var][ldof], 0, sizeof(double) * MDE);
  }

} /* END of routine zero_lec_row                                        */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/* Strange function, eh?  This routine will zero out the column of the
 * local element stiffness matrix corresponding to variables of type
 * var_type, and local degree of freedom ldof.
 *
 * Author: Matt Hopkins, 12/7/00
 */
void zero_lec_column(
    double local_J[MAX_PROB_VAR + MAX_CONC][MAX_PROB_VAR + MAX_CONC][MDE][MDE],
    int var_type, /* Variable type of column to be zeroed */
    int ldof)     /* Local dof of that variable type */
{
  int eqn, dof;

  for (eqn = 0; eqn < MAX_PROB_EQN + MAX_CONC; eqn++)
    for (dof = 0; dof < MDE; dof++)
      local_J[eqn][var_type][dof][ldof] = 0.0;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 * Function the returns the index in the VBR a array of the block
 * associated with the Ith block row and Jth block column */

int find_VBR_index(const int I, /* Block row index */
                   const int J, /* Block column index */
                   struct Aztec_Linear_Solver_System *ams) {
  int K = in_list(J, ams->bpntr[I], ams->bpntr[I + 1], ams->bindx);

  return (ams->indx[K]);

} /* END of find_VBR_index */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 * newshape() -- this is a wrapper routine for the Shadid's shape() routine.
 *		  We check for special basis functions (now, just pressure)
 *		  and otherwise hand off the real work to the old routine.
 *
 * Created:	Mon Mar 14 16:08:23 MST 1994 pasacki@sandia.gov
 *
 */

double newshape(const double xi[],    /* local coordinates    */
                const int Ielem_type, /* element type */
                const int Iquant, /* desired quantity (phi, phi_s, etc.     */
                const int Inode, /* current element node                      */
                const int eshape,        /* element shape        */
                const int interpolation, /* interpolation */
                const int ledof) /* Typically, this is the local node num     *
                                  * but for pressure basis functions, this    *
                                  * can be a dof at the centroid node         */
/***************************************************************************
 *
 *  newshape():
 *
 *	Routine to calculate the value of the shape function psi, and
 *	it's partial derivatives, dpsi_s, dpsi_t and dpsi_u if applicable
 *	at the point (s,t,u) on the master element [-1,1]x[-1,1]x[-1,1].
 *
 *  Input
 * --------
 *   xi[DIM] -> This is the vector of local element coordinates. The
 *              basis functions and derivatives of the basis functions
 *              are evaluated at the position in the element indicated
 *              by this vector.
 *
 *	Created:	Tue Mar 15 07:38:13 MST 1994 pasacki@sandia.gov
 *
 *	Revised:	Tue Mar 29 08:44:51 MST 1994 pasacki@sandia.gov
 *
 * Revised:        Tue Aug 25 09:55:21 NST 2009 prschun@sandia.gov
 *************************************************************************/
{
  double value = 0.0;
  const double s = xi[0];
  const double t = xi[1];
  const double u = xi[2];

  /*
   * This wrapper routine catches the case of I_P1 basis functions for
   * pressure and hands off all other cases to Shadid's routines.
   */
  switch (eshape) {

  case LINE_SEGMENT:
    if (interpolation == I_Q1) {
      switch (Iquant) {
      case PSI:
        value = ((Inode == 1) ? 0.5 * (1.0 + s) : 0.5 * (1.0 - s));
        break;
      case DPSI_S:
        value = ((Inode == 1) ? 0.5 : -0.5);
        break;
      default:
        /*EH(-1, "Not a valid Iquant choice!");*/ /*Clever way out of this, for
                                                     Iquant=2 */
        break;
      }
    } else if (interpolation == I_Q2) {
      switch (Iquant) {
      case PSI: {
        switch (Inode) {
        case 0:
          value = 0.5 * s * (s - 1.0);
          break;
        case 1:
          value = 0.5 * s * (s + 1.0);
          break;
        case 2:
          value = 1.0 - s * s;
          break;
        default:
          break;
        }
        break;
      }
      case DPSI_S: {
        switch (Inode) {
        case 0:
          value = s - 0.5;
          break;
        case 1:
          value = s + 0.5;
          break;
        case 2:
          value = -2.0 * s;
          break;
        default:
          break;
        }
        break;
      }
      default:
        break;
      }
    } else {
      EH(-1, "Only LINEAR_BAR or QUAD_BAR 1D elements supported now!");
    }
    break;

  case SHELL:
    if (is_xfem_interp(interpolation)) {
      value = extended_shape(xi, Ielem_type, Iquant, Inode, eshape,
                             interpolation, ledof);
    } else if (interpolation == I_Q1 || interpolation == I_Q1_D) {
      /* because we may want to solve with Q1 basis functions on a BIQUAD_QUAD
       * mesh */
      if (Iquant < 3) {
        value = shape(s, t, u, BILINEAR_QUAD, Iquant, Inode);
      } else {
        value = 0;
      }
    } else if ((interpolation == I_Q2) || (interpolation == I_Q2_D) ||
               (interpolation == I_Q2_D_LSA)) {
      /*
       * because we may want to solve with Q2 with discontinuous
       * boundaries BIQUAD_QUAD mesh.
       * Note here that Inode is coming from dof_list[v][ldof] which is
       * {0 1 1 2 2 3 4 5 5 6 7 8} for a biquadratic surface element
       */
      if (Iquant < 3) {
        value = shape(s, t, u, BIQUAD_QUAD, Iquant, Inode);
      } else {
        value = 0;
      }

      /* now simply zero out the second entry of each double number using
         baby_dolphin, but to it up above */

    } else {
      EH(-1, "Only Q1 and Q2 elements allowed for shells");
    }
    break;

  case TRIANGLE:
    if (interpolation == I_Q1) {
      value = shape(s, t, u, LINEAR_TRI, Iquant, Inode);
    } else if (interpolation == I_Q2) {
      value = shape(s, t, u, QUAD_TRI, Iquant, Inode);
    } else {
      EH(-1, "Don't recognize this basis type for triangles");
    }

    break;

  case TRISHELL:
    if (interpolation == I_Q1) {
      if (Iquant < 3) {
        value = shape(s, t, u, BILINEAR_TRISHELL, Iquant, Inode);
      } else {
        value = 0;
      }
    } else {
      EH(-1, "Don't recognize this basis type for linear triangular shells");
    }

    break;

  case QUADRILATERAL:
    if (is_xfem_interp(interpolation)) {
      value = extended_shape(xi, Ielem_type, Iquant, Inode, eshape,
                             interpolation, ledof);
    } else if (interpolation == I_P0) {
      /*
       * Shape function for peicewise constant pressures
       */
      switch (Iquant) {
      case PSI:
        value = 1.;
        break;


      case DPSI_S:
        value = 0.;
        break;

      case DPSI_T:
        value = 0.;
        break;

      default:
        EH(-1, "Bad pressure quantity specified.");
        break;
      }
    } else if (interpolation == I_P1) {
      /*
       * Shape function for linear discontinuous pressures
       */
      switch (Iquant) {
      case PSI:
        switch (ledof) {
        case 0:
          value = 1.;
          break;

        case 1:
          value = s;
          break;

        case 2:
          value = t;
          break;
        }
        break;

      case DPSI_S:
        switch (ledof) {
        case 0:
          value = 0.;
          break;

        case 1:
          value = 1.;
          break;

        case 2:
          value = 0.;
          break;
        }
        break;

      case DPSI_T:
        switch (ledof) {
        case 0:
          value = 0.;
          break;

        case 1:
          value = 0.;
          break;

        case 2:
          value = 1.;
          break;
        }
        break;

      default:
        EH(-1, "Bad pressure quantity specified.");
        break;
      }
    } else if (interpolation == I_Q1 || interpolation == I_Q1_D) {
      /* because we may want to solve with Q1 basis functions on a BIQUAD_QUAD
       * mesh */
      value = shape(s, t, u, BILINEAR_QUAD, Iquant, Inode);
    } else if ((interpolation == I_Q2) || (interpolation == I_Q2_D) ||
               (interpolation == I_Q2_D_LSA)) {
      /*
       * because we may want to solve with Q2 with discontinuous
       * boundaries BIQUAD_QUAD mesh.
       * Note here that Inode is coming from dof_list[v][ldof] which is
       * {0 1 1 2 2 3 4 5 5 6 7 8} for a biquadratic surface element
       */
      value = shape(s, t, u, BIQUAD_QUAD, Iquant, Inode);

      /* now simply zero out the second entry of each double number using
         baby_dolphin, but to it up above */

    } else if (interpolation == I_SP) {
      /* Subparametric interpolation which is linear on interior and quadratic
       * on all external surfaces */

      if (Inode >= 8)
        EH(-1, "This element has at least 8 nodes!");

      switch (Iquant) {  /* select quantity */
      case PSI:          /* shape function */
        switch (Inode) { /* select specific shape function */
        case 0:
          if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 4]]->EDGE !=
                  1 &&
              Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 7]]->EDGE !=
                  1) {
            value = 0.25 * (1. - s) * (1. - t);
          } else {
            value = 0.25 *
                    ((1. - s) * s * (1. - t) * t + (1 - s * s) * (1 - t * t));
            if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 4]]
                    ->EDGE != 1) {
              value += 0.25 * (1. - s * s) * (t - 1.) * t;
            }
            if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 7]]
                    ->EDGE != 1) {
              value += 0.25 * (1. - t * t) * (s - 1.) * s;
            }
          }
          break;
        case 1:
          if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 4]]->EDGE !=
                  1 &&
              Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 5]]->EDGE !=
                  1) {
            value = 0.25 * (1. + s) * (1. - t);
          } else {
            value =
                0.25 * ((1. + s) * s * (t - 1) * t + (1 - s * s) * (1 - t * t));
            if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 4]]
                    ->EDGE != 1) {
              value += 0.25 * (1. - s * s) * (t - 1.) * t;
            }
            if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 5]]
                    ->EDGE != 1) {
              value += 0.25 * (1. - t * t) * (s + 1.) * s;
            }
          }
          break;
        case 2:
          if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 5]]->EDGE !=
                  1 &&
              Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 6]]->EDGE !=
                  1) {
            value = 0.25 * (1. + s) * (1. + t);
          } else {
            value =
                0.25 * ((1. + s) * s * (t + 1) * t + (1 - s * s) * (1 - t * t));
            if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 6]]
                    ->EDGE != 1) {
              value += 0.25 * (1. - s * s) * (t + 1.) * t;
            }
            if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 5]]
                    ->EDGE != 1) {
              value += 0.25 * (1. - t * t) * (s + 1.) * s;
            }
          }
          break;
        case 3:
          if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 6]]->EDGE !=
                  1 &&
              Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 7]]->EDGE !=
                  1) {
            value = 0.25 * (1. - s) * (1. + t);
          } else {
            value =
                0.25 * ((s - 1) * s * (t + 1) * t + (1 - s * s) * (1 - t * t));
            if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 6]]
                    ->EDGE != 1) {
              value += 0.25 * (1. - s * s) * (t + 1.) * t;
            }
            if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 7]]
                    ->EDGE != 1) {
              value += 0.25 * (1. - t * t) * (s - 1.) * s;
            }
          }
          break;
          /* can only get to these cases on external boundaries */
        case 4:
          if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 4]]->EDGE !=
              1)
            EH(-1, "Subparametric node not on edge ");
          value = 0.5 * (1. - s * s) * (t - 1.) * t;
          break;
        case 5:
          if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 5]]->EDGE !=
              1)
            EH(-1, "Subparametric node not on edge ");
          value = 0.5 * (1. - t * t) * (s + 1.) * s;
          break;
        case 6:
          if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 6]]->EDGE !=
              1)
            EH(-1, "Subparametric node not on edge ");
          value = 0.5 * (1. - s * s) * (t + 1.) * t;
          break;
        case 7:
          if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 7]]->EDGE !=
              1)
            EH(-1, "Subparametric node not on edge ");
          value = 0.5 * (1. - t * t) * (s - 1.) * s;
          break;
        }
        break;

      case DPSI_S:       /* partial of shape fn w.r.t. s */
        switch (Inode) { /* select specific shape function */
        case 0:
          if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 4]]->EDGE !=
                  1 &&
              Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 7]]->EDGE !=
                  1) {
            value = -0.25 * (1. - t);
          } else {
            value = 0.25 * ((1. - 2 * s) * (1. - t) * t - 2 * s * (1 - t * t));
            if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 4]]
                    ->EDGE != 1) {
              value -= 0.5 * s * (t - 1.) * t;
            }
            if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 7]]
                    ->EDGE != 1) {
              value += 0.25 * (1. - t * t) * (2. * s - 1.);
            }
          }
          break;
        case 1:
          if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 4]]->EDGE !=
                  1 &&
              Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 5]]->EDGE !=
                  1) {
            value = 0.25 * (1. - t);
          } else {
            value = 0.25 * ((1. + 2 * s) * (t - 1) * t - 2 * s * (1 - t * t));
            if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 4]]
                    ->EDGE != 1) {
              value -= 0.5 * s * (t - 1.) * t;
            }
            if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 5]]
                    ->EDGE != 1) {
              value += 0.25 * (1. - t * t) * (2 * s + 1.);
            }
          }
          break;
        case 2:
          if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 5]]->EDGE !=
                  1 &&
              Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 6]]->EDGE !=
                  1) {
            value = 0.25 * (1. + t);
          } else {
            value = 0.25 * ((1. + 2 * s) * (t + 1) * t - 2 * s * (1 - t * t));
            if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 6]]
                    ->EDGE != 1) {
              value -= 0.5 * s * (t + 1.) * t;
            }
            if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 5]]
                    ->EDGE != 1) {
              value += 0.25 * (1. - t * t) * (2 * s + 1.);
            }
          }
          break;
        case 3:
          if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 6]]->EDGE !=
                  1 &&
              Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 7]]->EDGE !=
                  1) {
            value = -0.25 * (1. + t);
          } else {
            value = 0.25 * ((2 * s - 1) * (t + 1) * t - 2 * s * (1 - t * t));
            if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 6]]
                    ->EDGE != 1) {
              value -= 0.5 * s * (t + 1.) * t;
            }
            if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 7]]
                    ->EDGE != 1) {
              value += 0.25 * (1. - t * t) * (2 * s - 1.);
            }
          }
          break;
          /* can only get to these cases on external boundaries */
        case 4:
          if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 4]]->EDGE !=
              1)
            EH(-1, "Subparametric node not on edge ");
          value = -s * (t - 1.) * t;
          break;
        case 5:
          if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 5]]->EDGE !=
              1)
            EH(-1, "Subparametric node not on edge ");
          value = 0.5 * (1. - t * t) * (2 * s + 1.);
          break;
        case 6:
          if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 6]]->EDGE !=
              1)
            EH(-1, "Subparametric node not on edge ");
          value = -s * (t + 1.) * t;
          break;
        case 7:
          if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 7]]->EDGE !=
              1)
            EH(-1, "Subparametric node not on edge ");
          value = 0.5 * (1. - t * t) * (2 * s - 1.);
          break;
        }
        break;

      case DPSI_T:       /* partial of shape fn w.r.t. t */
        switch (Inode) { /* select specific shape function */
        case 0:
          if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 4]]->EDGE !=
                  1 &&
              Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 7]]->EDGE !=
                  1) {
            value = -0.25 * (1. - s);
          } else {
            value = 0.25 * ((1. - s) * s * (1. - 2 * t) - (1 - s * s) * 2 * t);
            if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 4]]
                    ->EDGE != 1) {
              value += 0.25 * (1. - s * s) * (2 * t - 1.);
            }
            if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 7]]
                    ->EDGE != 1) {
              value -= 0.5 * t * (s - 1.) * s;
            }
          }
          break;
        case 1:
          if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 4]]->EDGE !=
                  1 &&
              Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 5]]->EDGE !=
                  1) {
            value = -0.25 * (1. + s);
          } else {
            value = 0.25 * ((1. + s) * s * (2 * t - 1) - (1 - s * s) * 2 * t);
            if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 4]]
                    ->EDGE != 1) {
              value += 0.25 * (1. - s * s) * (2 * t - 1.);
            }
            if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 5]]
                    ->EDGE != 1) {
              value -= 0.5 * t * (s + 1.) * s;
            }
          }
          break;
        case 2:
          if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 5]]->EDGE !=
                  1 &&
              Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 6]]->EDGE !=
                  1) {
            value = 0.25 * (1. + s);
          } else {
            value = 0.25 * ((1. + s) * s * (2 * t + 1) - (1 - s * s) * 2 * t);
            if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 6]]
                    ->EDGE != 1) {
              value += 0.25 * (1. - s * s) * (2 * t + 1.);
            }
            if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 5]]
                    ->EDGE != 1) {
              value -= 0.5 * t * (s + 1.) * s;
            }
          }
          break;
        case 3:
          if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 6]]->EDGE !=
                  1 &&
              Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 7]]->EDGE !=
                  1) {
            value = 0.25 * (1. - s);
          } else {
            value = 0.25 * ((s - 1) * s * (2 * t + 1) - (1 - s * s) * 2 * t);
            if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 6]]
                    ->EDGE != 1) {
              value += 0.25 * (1. - s * s) * (2 * t + 1.);
            }
            if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 7]]
                    ->EDGE != 1) {
              value -= 0.5 * t * (s - 1.) * s;
            }
          }
          break;
          /* can only get to these cases on external boundaries */
        case 4:
          if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 4]]->EDGE !=
              1)
            EH(-1, "Subparametric node not on edge ");
          value = 0.5 * (1. - s * s) * (2 * t - 1.);
          break;
        case 5:
          if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 5]]->EDGE !=
              1)
            EH(-1, "Subparametric node not on edge ");
          value = -t * (s + 1.) * s;
          break;
        case 6:
          if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 6]]->EDGE !=
              1)
            EH(-1, "Subparametric node not on edge ");
          value = 0.5 * (1. - s * s) * (2 * t + 1.);
          break;
        case 7:
          if (Nodes[Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 7]]->EDGE !=
              1)
            EH(-1, "Subparametric node not on edge ");
          value = -t * (s - 1.) * s;
          break;
        }
        break;

      default:
        fprintf(stderr, "Bad BIQUAD_QUAD case: %d!\n", Iquant);
        EH(-1, "Bad selection of phi,dphids, etc.");
        break;
      }
    } else if (interpolation == I_PQ1) {
      value = shape(s, t, u, BILINEAR_QUAD, Iquant, ledof);
    } else if (interpolation == I_PQ2) {
      value = shape(s, t, u, BIQUAD_QUAD, Iquant, ledof);
    } else {
      value = shape(s, t, u, Ielem_type, Iquant, Inode);
    }
    break;
  case TETRAHEDRON:
    if (interpolation == I_Q1 || interpolation == I_Q1_D) {
      value = shape(s, t, u, LINEAR_TET, Iquant, Inode);
    } else if (interpolation == I_P0) {
      // Shape function for peicewise constant interpolation
      switch (Iquant) {
      case PSI:
        value = 1.;
        break;

      case DPSI_S:
        value = 0.;
        break;

      case DPSI_T:
        value = 0.;
        break;

      case DPSI_U:
        value = 0.;
        break;

      default:
        EH(-1, "Bad pressure quantity specified.");
        break;
      }
    } else {
      EH(-1, "Don't recognize this basis type for Linear tetrahedron");
    }

    break;

  case HEXAHEDRON:
    if (is_xfem_interp(interpolation)) {
      value = extended_shape(xi, Ielem_type, Iquant, Inode, eshape,
                             interpolation, ledof);
    } else if (interpolation == I_P0) {
      /*
       * Shape function for peicewise constant interpolation
       */
      switch (Iquant) {
      case PSI:
        value = 1.;
        break;

      case DPSI_S:
        value = 0.;
        break;

      case DPSI_T:
        value = 0.;
        break;

      case DPSI_U:
        value = 0.;
        break;

      default:
        EH(-1, "Bad pressure quantity specified.");
        break;
      }
    } else if (interpolation == I_P1) {
      /*
       * Shape function for linear discontinuous pressures
       */
      switch (Iquant) {
      case PSI:
        switch (ledof) {
        case 0:
          value = 1.;
          break;

        case 1:
          value = s;
          break;

        case 2:
          value = t;
          break;

        case 3:
          value = u;
          break;
        }
        break;

      case DPSI_S:
        switch (ledof) {
        case 0:
        case 2:
        case 3:
          value = 0.;
          break;

        case 1:
          value = 1.;
          break;
        }
        break;

      case DPSI_T:
        switch (ledof) {
        case 0:
        case 1:
        case 3:
          value = 0.;
          break;

        case 2:
          value = 1.;
          break;
        }
        break;

      case DPSI_U:
        switch (ledof) {
        case 0:
        case 1:
        case 2:
          value = 0.;
          break;

        case 3:
          value = 1.;
          break;
        }
        break;

      default:
        EH(-1, "Bad pressure quantity specified.");
        break;
      }
    }

    else if (interpolation == I_Q1 || interpolation == I_Q1_D) {
      /* because we may want to solve with Q1 basis functions on a TRIQUAD_HEX
       * mesh */
      value = shape(s, t, u, TRILINEAR_HEX, Iquant, Inode);
    } else {
      value = shape(s, t, u, Ielem_type, Iquant, Inode);
    }
    break;

  default:
    WH(-1, "Called shape() without determining the proper element type");
    value = shape(s, t, u, Ielem_type, Iquant, Inode);
    if (ei[pg->imtrx]->ielem > 7)
      printf("      Shape() call returned %g!\n", value);
    break;
  }

  return (value);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* shape function for XFEM */
double
extended_shape(const double xi[],    /* local coordinates    */
               const int Ielem_type, /* element type */
               const int Iquant, /* desired quantity (phi, phi_s, etc.        */
               const int Inode,  /* current element node                      */
               const int eshape, /* element shape                             */
               const int interpolation, /* interpolation */
               const int ledof) /* Typically, this is the local node num     *
                                 * but for pressure basis functions, this    *
                                 * can be a dof at the centroid node         */
{
  int xfem_active = FALSE; /* innocent till proven guilty, how american! */
  int base_interp, base_dof, extended_dof;
  double F_i;
  double value = 0.0;

  load_xfem_for_stu(xi);

  xfem_dof_state(ledof, interpolation, eshape, &xfem_active, &extended_dof,
                 &base_interp, &base_dof);

  if (xfem_active) {
    double H = xfem->H;
    double delta = xfem->delta;

    if (ls->Elem_Sign == -1) {
      H = 0.;
      delta = 0.;
    } else if (ls->Elem_Sign == 1) {
      H = 1.;
      delta = 0.;
    }

    switch (interpolation) {
    case I_P0_G:
    case I_P1_G:
    case I_Q1_G:
    case I_Q2_G:
      /*F_i = dof_distance( Ielem_type, base_interp, Inode );*/
      F_i = lnn_distance(Inode);
      value = newshape(xi, Ielem_type, Iquant, Inode, eshape, base_interp,
                       base_dof);
      if ((!extended_dof && F_i < 0.) || (extended_dof && F_i >= 0.)) {
        value *= 1. - H;
        if (delta != 0. &&
            (Iquant == DPSI_S || Iquant == DPSI_T || Iquant == DPSI_U)) {
          double psi = newshape(xi, Ielem_type, PSI, Inode, eshape, base_interp,
                                base_dof);
          if (Iquant == DPSI_S)
            value += psi * delta * xfem->dF_xi[0];
          if (Iquant == DPSI_T)
            value += psi * delta * xfem->dF_xi[1];
          if (Iquant == DPSI_U)
            value += psi * delta * xfem->dF_xi[2];
        }
      } else {
        value *= H;
        if (delta != 0. &&
            (Iquant == DPSI_S || Iquant == DPSI_T || Iquant == DPSI_U)) {
          double psi = newshape(xi, Ielem_type, PSI, Inode, eshape, base_interp,
                                base_dof);
          if (Iquant == DPSI_S)
            value -= psi * delta * xfem->dF_xi[0];
          if (Iquant == DPSI_T)
            value -= psi * delta * xfem->dF_xi[1];
          if (Iquant == DPSI_U)
            value -= psi * delta * xfem->dF_xi[2];
        }
      }
      break;
    case I_P0_GP:
    case I_P1_GP:
    case I_Q1_GP:
    case I_Q2_GP:
      value = newshape(xi, Ielem_type, Iquant, Inode, eshape, base_interp,
                       base_dof);
      value *= H;
      if (delta != 0. &&
          (Iquant == DPSI_S || Iquant == DPSI_T || Iquant == DPSI_U)) {
        double psi =
            newshape(xi, Ielem_type, PSI, Inode, eshape, base_interp, base_dof);
        if (Iquant == DPSI_S)
          value -= psi * delta * xfem->dF_xi[0];
        if (Iquant == DPSI_T)
          value -= psi * delta * xfem->dF_xi[1];
        if (Iquant == DPSI_U)
          value -= psi * delta * xfem->dF_xi[2];
      }
      break;
    case I_P0_GN:
    case I_P1_GN:
    case I_Q1_GN:
    case I_Q2_GN:
      value = newshape(xi, Ielem_type, Iquant, Inode, eshape, base_interp,
                       base_dof);
      value *= 1. - H;
      if (delta != 0. &&
          (Iquant == DPSI_S || Iquant == DPSI_T || Iquant == DPSI_U)) {
        double psi =
            newshape(xi, Ielem_type, PSI, Inode, eshape, base_interp, base_dof);
        if (Iquant == DPSI_S)
          value -= psi * delta * xfem->dF_xi[0];
        if (Iquant == DPSI_T)
          value -= psi * delta * xfem->dF_xi[1];
        if (Iquant == DPSI_U)
          value -= psi * delta * xfem->dF_xi[2];
      }
      break;

#if 0
	  /* this is the chessa enrichment */
        case I_P1_XG:
        case I_Q1_XG:
        case I_Q2_XG:
          {
            value = newshape( xi, Ielem_type, Iquant, Inode, eshape, base_interp, base_dof );

            if ( extended_dof )
              {
                double S = 2. * H - 1.;

                /*F_i = dof_distance( Ielem_type, base_interp, Inode );*/
	        F_i = lnn_distance( Inode );
                value *=  S * xfem->F - fabs( F_i );

                if (Iquant == DPSI_S || Iquant == DPSI_T || Iquant == DPSI_U)
                  {
                    double psi = newshape( xi, Ielem_type, PSI, Inode, eshape, base_interp, base_dof );
                    if ( Iquant == DPSI_S ) value += psi * (2.*xfem->F*delta + S) * xfem->dF_xi[0];
                    if ( Iquant == DPSI_T ) value += psi * (2.*xfem->F*delta + S) * xfem->dF_xi[1];
                    if ( Iquant == DPSI_U ) value += psi * (2.*xfem->F*delta + S) * xfem->dF_xi[2];
                  }
              }
          }
          break;
#endif
#if 1
      /* this is the moes enrichment */
    case I_P1_XG:
    case I_Q1_XG:
    case I_Q2_XG: {
      value = newshape(xi, Ielem_type, Iquant, Inode, eshape, base_interp,
                       base_dof);

      if (extended_dof) {
        value *= 2. * (xfem->F_plus - H * xfem->F);

        if (Iquant == DPSI_S || Iquant == DPSI_T || Iquant == DPSI_U) {
          double psi = newshape(xi, Ielem_type, PSI, Inode, eshape, base_interp,
                                base_dof);
          if (Iquant == DPSI_S)
            value +=
                psi * 2. *
                (xfem->grad_F_plus[0] - (xfem->F * delta + H) * xfem->dF_xi[0]);
          if (Iquant == DPSI_T)
            value +=
                psi * 2. *
                (xfem->grad_F_plus[1] - (xfem->F * delta + H) * xfem->dF_xi[1]);
          if (Iquant == DPSI_U)
            value +=
                psi * 2. *
                (xfem->grad_F_plus[2] - (xfem->F * delta + H) * xfem->dF_xi[2]);
        }
      }
    } break;
#endif

    case I_Q1_HV:
    case I_Q2_HV: {
      if (extended_dof) {
        if (Iquant == PSI)
          value = H - xfem->bf_plus;
        else if (Iquant == DPSI_S)
          value = delta * xfem->dF_xi[0] - xfem->grad_bf_plus[0];
        else if (Iquant == DPSI_T)
          value = delta * xfem->dF_xi[1] - xfem->grad_bf_plus[1];
        else if (Iquant == DPSI_U)
          value = delta * xfem->dF_xi[2] - xfem->grad_bf_plus[2];
      } else {
        value = newshape(xi, Ielem_type, Iquant, Inode, eshape, base_interp,
                         base_dof);
      }
    } break;

    case I_Q1_HG:
    case I_Q2_HG: {
      if (extended_dof) {
        if (Iquant == PSI)
          value = xfem->F * H - xfem->F_plus;
        else if (Iquant == DPSI_S)
          value = (H + xfem->F * delta) * xfem->dF_xi[0] - xfem->grad_F_plus[0];
        else if (Iquant == DPSI_T)
          value = (H + xfem->F * delta) * xfem->dF_xi[1] - xfem->grad_F_plus[1];
        else if (Iquant == DPSI_U)
          value = (H + xfem->F * delta) * xfem->dF_xi[2] - xfem->grad_F_plus[2];
      } else {
        value = newshape(xi, Ielem_type, Iquant, Inode, eshape, base_interp,
                         base_dof);
      }
    } break;

    case I_Q1_HVG:
    case I_Q2_HVG: {
      if (extended_dof) {
        if (base_dof == 0) {
          if (Iquant == PSI)
            value = H - xfem->bf_plus;
          else if (Iquant == DPSI_S)
            value = delta * xfem->dF_xi[0] - xfem->grad_bf_plus[0];
          else if (Iquant == DPSI_T)
            value = delta * xfem->dF_xi[1] - xfem->grad_bf_plus[1];
          else if (Iquant == DPSI_U)
            value = delta * xfem->dF_xi[2] - xfem->grad_bf_plus[2];
        } else {
          if (Iquant == PSI)
            value = xfem->F * H - xfem->F_plus;
          else if (Iquant == DPSI_S)
            value =
                (H + xfem->F * delta) * xfem->dF_xi[0] - xfem->grad_F_plus[0];
          else if (Iquant == DPSI_T)
            value =
                (H + xfem->F * delta) * xfem->dF_xi[1] - xfem->grad_F_plus[1];
          else if (Iquant == DPSI_U)
            value =
                (H + xfem->F * delta) * xfem->dF_xi[2] - xfem->grad_F_plus[2];
        }
      } else {
        value = newshape(xi, Ielem_type, Iquant, Inode, eshape, base_interp,
                         base_dof);
      }
    } break;

    case I_P0_XV:
    case I_P1_XV:
    case I_Q1_XV:
    case I_Q2_XV: {
      value = newshape(xi, Ielem_type, Iquant, Inode, eshape, base_interp,
                       base_dof);

      if (extended_dof) {
        double H_i = 0.0;

        /*F_i = dof_distance( Ielem_type, base_interp, Inode );*/
        F_i = lnn_distance(Inode);

        if (F_i >= 0.)
          H_i = 1.;
        value *= H - H_i;
        if (delta != 0. &&
            (Iquant == DPSI_S || Iquant == DPSI_T || Iquant == DPSI_U)) {
          double psi = newshape(xi, Ielem_type, PSI, Inode, eshape, base_interp,
                                base_dof);
          if (Iquant == DPSI_S)
            value += psi * delta * xfem->dF_xi[0];
          if (Iquant == DPSI_T)
            value += psi * delta * xfem->dF_xi[1];
          if (Iquant == DPSI_U)
            value += psi * delta * xfem->dF_xi[2];
        }
      }
    } break;
    }
  } else {
    if (extended_dof) {
      value = 0.;
#if 0
	  if ( interpolation == I_P0_GP ||
	       interpolation == I_P1_GP ||
	       interpolation == I_Q1_GP ||
	       interpolation == I_Q2_GP )
	    {
	      if ( lnn_distance( Inode ) >= 0. )
	        value = newshape( xi, Ielem_type, Iquant, Inode, eshape, base_interp, base_dof );
            }
	  if ( interpolation == I_P0_GN ||
	       interpolation == I_P1_GN ||
	       interpolation == I_Q1_GN ||
	       interpolation == I_Q2_GN )
	    {
	      if ( lnn_distance( Inode ) < 0. )
	        value = newshape( xi, Ielem_type, Iquant, Inode, eshape, base_interp, base_dof );
            }
#endif
    } else {
      value = newshape(xi, Ielem_type, Iquant, Inode, eshape, base_interp,
                       base_dof);
    }
  }
  return value;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int
calc_shearrate(dbl *gammadot,	/* strain rate invariant */
	       dbl gamma_dot[DIM][DIM], /* strain rate tensor */
	       dbl d_gd_dv[DIM][MDE],
	       dbl d_gd_dmesh[DIM][MDE])
{
  int mdofs = 0;
  int p, q, a, b;
  int vdofs, i, j, v;

  dbl grad_phi_e_gam[MDE][DIM] [DIM][DIM]; /* transpose of grad(phi_i ea) tensor
					      + grad(phi_i ea) tensor */
  dbl d_gamma_dot_dmesh [DIM][DIM] [DIM][MDE]; /* d/dmesh(grad_v)T */

  int status = 1;


  /* Zero out sensitivities */

  if(d_gd_dv != NULL) memset(d_gd_dv, 0, sizeof(double)*DIM*MDE);
  if(d_gd_dmesh != NULL) memset(d_gd_dmesh, 0, sizeof(double)*DIM*MDE);


  *gammadot = 0.;
  /* get gamma_dot invariant for viscosity calculations */
  for ( a=0; a<VIM; a++)
    {
      for ( b=0; b<VIM; b++)
	{
	  *gammadot +=  gamma_dot[a][b] * gamma_dot[b][a];
	}
    }
  
  *gammadot  =  sqrt(0.5*fabs(*gammadot)); 
  
  /* get stuff for Jacobian entries */
  v = VELOCITY1;
  vdofs = ei[pg->imtrx]->dof[v];
  
  if ( d_gd_dmesh != NULL || d_gd_dv != NULL)
  {
  if ( pd->v[pg->imtrx][R_MESH1] )
    {
      mdofs = ei[pg->imtrx]->dof[R_MESH1];
    }
  
  for ( p=0; p<VIM; p++)
    {
      for ( q=0; q<VIM; q++)
	{
	  for ( a=0; a<VIM; a++)
	    {
	      for ( i=0; i<vdofs; i++)
		{
		  grad_phi_e_gam[i][a] [p][q] =
		    bf[v]->grad_phi_e[i][a] [p][q]
		    + bf[v]->grad_phi_e[i][a] [q][p]  ;
		}
	    }
	}
    }
  }
  
  /*
   * d( gamma_dot )/dmesh
   */
  
  if ( pd->v[pg->imtrx][R_MESH1] && d_gd_dmesh != NULL)
    {
      
      for ( p=0; p<VIM; p++)
	{
	  for ( q=0; q<VIM; q++)
	    {
	      for ( b=0; b<VIM; b++)
		{
		  for ( j=0; j<mdofs; j++)
		    {
		      
		      d_gamma_dot_dmesh[p][q] [b][j] =
			fv->d_grad_v_dmesh[p][q] [b][j] +
			fv->d_grad_v_dmesh[q][p] [b][j] ;
		    }
		}
	    }
	}
      
      /*
       * d( gammadot )/dmesh
       */
      
      if(*gammadot != 0.)
	{
	  for ( b=0; b<VIM; b++)
	    {
	      for ( j=0; j<mdofs; j++)
		{
		  d_gd_dmesh [b][j] = 0.;
		  for ( p=0; p<VIM; p++)
		    {
		      for ( q=0; q<VIM; q++)
			{
			  d_gd_dmesh [b][j] +=
			    0.5 * d_gamma_dot_dmesh[p][q] [b][j] 
			    * gamma_dot[q][p]
			    / *gammadot; 
		      
			}
		    }
		}
	    }
	}
    }
  
  /*
   * d( gammadot )/dv
   */
  
  if(*gammadot != 0. && d_gd_dv != NULL)
    {
      for ( a=0; a<VIM; a++)
	{
	  for ( i=0; i<vdofs; i++)
	    {
	      d_gd_dv[a][i] = 0.;
	      for ( p=0; p<VIM; p++)
		{
		  for ( q=0; q<VIM; q++)
		    {
		      d_gd_dv[a][i] +=
			0.5 * grad_phi_e_gam[i][a] [p][q] 
			* gamma_dot[q][p]
			/ *gammadot; 
		    }
		}
	    }
	}
    }
  return(status);
}

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

static int var_if_interp_type_enabled(PROBLEM_DESCRIPTION_STRUCT *pd_ptr,
                               int interp_type) {
  int imtrx;
  int var;
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    if ((var = in_list(interp_type, 0, MAX_VARIABLE_TYPES, pd_ptr->i[imtrx])) !=
        -1) {
      return var;
    }
  }
  return -1;
}

void determine_ShapeVar(PROBLEM_DESCRIPTION_STRUCT *pd_ptr)

/*********************************************************************
 *
 * determine_ShapeVar():
 *
 *  Part 6 --
 *
 *   Now figure out the mapping order.
 *
 *   - IntegrationMap variable is determined here. This is the
 *                    interpolation type that will be used in the
 *                    formulation of the isoparametric mapping.
 *   - ShapeVar       variable also. This is the variable which
 *                    has the interpolant specified by IntegrationMap.
 ***********************************************************************/
{
  if (pd_ptr->IntegrationMap == ISOPARAMETRIC) {

    int imtrx;
    int mesh_matrix = -1;

    for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
      if (pd_ptr->v[imtrx][R_MESH1]) {
        mesh_matrix = imtrx;
        break;
      }
    }

    /*
     *  For shell variables e[] may be zero, but there may be a v[].
     */
    if (mesh_matrix >= 0) {
      /*
       * If deforming mesh always make displacement the mapping
       * interpolation
       */
      pd_ptr->IntegrationMap = pd_ptr->i[mesh_matrix][R_MESH1];
      pd_ptr->ShapeVar = R_MESH1;
    } else {
      if ((pd_ptr->ShapeVar = var_if_interp_type_enabled(pd_ptr, I_Q2)) != -1) {
        pd_ptr->IntegrationMap = I_Q2;
      } else if ((pd_ptr->ShapeVar =
                      var_if_interp_type_enabled(pd_ptr, I_Q1)) != -1) {
        pd_ptr->IntegrationMap = I_Q1;
      } else if ((pd_ptr->ShapeVar =
                      var_if_interp_type_enabled(pd_ptr, I_SP)) != -1) {
        pd_ptr->IntegrationMap = I_SP;
      } else {
        pd_ptr->ShapeVar = pd_ptr->m[pg->imtrx][0];
        pd_ptr->IntegrationMap = pd_ptr->i[pg->imtrx][pd_ptr->ShapeVar];
      }
    }
  } else {
    /*
     * This is for the case where the Q1, Q2 or SP has been specified
     * as the mapping order
     */
    pd_ptr->ShapeVar = in_list(pd_ptr->IntegrationMap, 0, MAX_VARIABLE_TYPES,
                               pd_ptr->i[pg->imtrx]);
    if (pd_ptr->ShapeVar == -1) {
      EH(-1,
         "Error: Specified Element Mapping has no corresponding variable .");
    }
    if (pd_ptr->e[pg->imtrx][R_MESH1] &&
        pd_ptr->IntegrationMap != pd_ptr->i[0][R_MESH1]) {
      WH(-1, "Warning: Having the Element Mapping differ from the "
             "displacement interpolation may produce unexpected results.");
    }
  }
}
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

void determine_ProjectionVar(PROBLEM_DESCRIPTION_STRUCT *pd_ptr)

/*********************************************************************
 *
 * determine_ProjectionVar():
 *
 *   Now figure out the which variable's basis functions will be used
 *   for the projection operator for this element block.
 *
 *   The projection interpolation function should be an interpolation
 *   which spans all of the local nodes in an element.
 *
 *  Order: Q2, I_Q2_D, I_SP, I_Q1, I_Q1_D
 *
 * Note from MMH to self: Since I_Q2_LSA and I_Q2_D_LSA are only
 * used when another velocity has I_Q2 or I_Q2_D, I don't have to
 * add checks for I_Q2_LSA and I_Q2_D_LSA here...
 *
 ***********************************************************************/
{
  int var;
  pd_ptr->ProjectionVar = -1;
  if (((var = var_if_interp_type_enabled(pd_ptr, I_Q2)) != -1) ||
      ((var = var_if_interp_type_enabled(pd_ptr, I_Q2_D)) != -1) ||
      ((var = var_if_interp_type_enabled(pd_ptr, I_Q2_G)) != -1) ||
      ((var = var_if_interp_type_enabled(pd_ptr, I_Q2_GP)) != -1) ||
      ((var = var_if_interp_type_enabled(pd_ptr, I_Q2_GN)) != -1) ||
      ((var = var_if_interp_type_enabled(pd_ptr, I_Q2_XV)) != -1) ||
      ((var = var_if_interp_type_enabled(pd_ptr, I_SP)) != -1) ||
      ((var = var_if_interp_type_enabled(pd_ptr, I_Q1)) != -1) ||
      ((var = var_if_interp_type_enabled(pd_ptr, I_Q1_D)) != -1) ||
      ((var = var_if_interp_type_enabled(pd_ptr, I_Q1_G)) != -1) ||
      ((var = var_if_interp_type_enabled(pd_ptr, I_Q1_GP)) != -1) ||
      ((var = var_if_interp_type_enabled(pd_ptr, I_Q1_GN)) != -1) ||
      ((var = var_if_interp_type_enabled(pd_ptr, I_Q1_XV)) != -1)) {
    pd_ptr->ProjectionVar = var;
  } else {
    P0PRINTF("Warning: No suitable basis function was found for a Projection "
             "Operation\n");
  }
#ifdef DEBUG_HKM
  P0PRINTF("Projection Variable set to %d with interp type = %d\n",
           pd_ptr->ProjectionVar, pd_ptr->i[pg->imtrx][pd_ptr->ProjectionVar]);
#endif
}

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

void set_solid_inertia(void)
/* Initializes tran->solid_inertia flag */
{
  int mn;

  tran->solid_inertia = FALSE;
  for (mn = 0; mn < upd->Num_Mat; mn++) {
    if (pd->TimeIntegration != STEADY &&
        pd_glob[mn]->etm[pg->imtrx][R_MESH1][(LOG2_MASS)] &&
        pd_glob[mn]->MeshMotion == DYNAMIC_LAGRANGIAN) {
      tran->solid_inertia = TRUE;
    } else if (pd->TimeIntegration != STEADY &&
               pd_glob[mn]->etm[pg->imtrx][R_SOLID1][(LOG2_MASS)] &&
               pd_glob[mn]->MeshMotion == TOTAL_ALE) {
      tran->solid_inertia = TRUE;
    }
  }

  return;
}
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double calc_tensor_invariant(dbl T[DIM][DIM],       // Original tensor
                             dbl d_TI_dT[DIM][DIM], // Sensitivities
                             int Ii) // Which invariant to calculate
{

  // Define variables
  int i, j, a, b;
  dbl TI, M, tr_T;
  dbl Tdev[DIM][DIM], d_Tdev_dT[DIM][DIM][DIM][DIM];

  // Initialize
  TI = 0.0;
  memset(d_TI_dT, 0, sizeof(double) * DIM * DIM);
  memset(Tdev, 0, sizeof(double) * DIM * DIM);
  memset(d_Tdev_dT, 0, sizeof(double) * DIM * DIM * DIM * DIM);

  // Error checking
  if (DIM != 3)
    EH(-1, "calc_tensor_invariants() has not been instrumented for 1 or 2 "
           "dimensions");

  // Select invariant and calculate
  if (Ii == 1) { // First invariant

    for (i = 0; i < DIM; i++) {
      TI += T[i][i];
      d_TI_dT[i][i] = 1.0;
    }

  } else if (Ii == 2) { // Second invariant

    TI += T[0][0] * T[1][1] + T[1][1] * T[2][2] + T[0][0] * T[2][2];
    TI -= T[0][1] * T[1][0] + T[1][2] * T[2][1] + T[0][2] * T[2][0];

    d_TI_dT[0][0] = T[1][1] + T[2][2];
    d_TI_dT[0][1] = -T[1][0];
    d_TI_dT[0][2] = -T[2][0];
    d_TI_dT[1][0] = -T[0][1];
    d_TI_dT[1][1] = T[0][0] + T[2][2];
    d_TI_dT[1][2] = -T[2][1];
    d_TI_dT[2][0] = -T[0][2];
    d_TI_dT[2][1] = -T[1][2];
    d_TI_dT[2][2] = T[0][0] + T[1][1];

  } else if (Ii == 3) { // Third invariant

    TI += T[0][0] * (T[1][1] * T[2][2] - T[1][2] * T[2][1]);
    TI -= T[0][1] * (T[1][0] * T[2][2] - T[1][2] * T[2][0]);
    TI += T[0][2] * (T[1][0] * T[2][1] - T[1][1] * T[2][0]);

    d_TI_dT[0][0] = T[1][1] * T[2][2] - T[1][2] * T[2][1];
    d_TI_dT[0][1] = T[1][2] * T[2][0] - T[1][0] * T[2][2];
    d_TI_dT[0][2] = T[1][0] * T[2][1] - T[1][1] * T[2][0];
    d_TI_dT[1][0] = T[0][2] * T[2][1] - T[0][1] * T[2][2];
    d_TI_dT[1][1] = T[0][0] * T[2][2] - T[0][2] * T[2][0];
    d_TI_dT[1][2] = T[0][1] * T[2][0] - T[0][0] * T[2][1];
    d_TI_dT[2][0] = T[0][1] * T[1][2] - T[0][2] * T[1][1];
    d_TI_dT[2][1] = T[0][2] * T[1][0] - T[0][0] * T[1][2];
    d_TI_dT[2][2] = T[0][0] * T[1][1] - T[0][1] * T[1][0];

  } else if (Ii == 4) { // von Mises invariant, sqrt(3*II(T_dev))

    // Calculate deviatoric tensor component
    tr_T = 0.0;
    for (i = 0; i < DIM; i++)
      tr_T += T[i][i];
    for (i = 0; i < DIM; i++) {
      for (j = 0; j < DIM; j++) {
        Tdev[i][j] += T[i][j];
        d_Tdev_dT[i][j][i][j] += 1.0;
        if (i == j) {
          Tdev[i][j] -= tr_T / 3.0;
          for (a = 0; a < DIM; a++) {
            d_Tdev_dT[i][j][a][a] -= 1.0 / 3.0;
          }
        }
      }
    }

    // Calculate second invariant
    for (i = 0; i < DIM; i++) {
      for (j = 0; j < DIM; j++) {
        TI += 0.5 * Tdev[i][j] * Tdev[i][j];
        for (a = 0; a < DIM; a++) {
          for (b = 0; b < DIM; b++) {
            d_TI_dT[a][b] += 1.0 * Tdev[i][j] * d_Tdev_dT[i][j][a][b];
          }
        }
      }
    }

    // Convert to von Mises invariant
    if (TI > 1e-8) {
      M = sqrt(3.0) / (2.0 * sqrt(TI));
    } else {
      M = 0.0;
    }
    TI = sqrt(3 * TI);
    for (i = 0; i < DIM; i++) {
      for (j = 0; j < DIM; j++) {
        d_TI_dT[i][j] *= M;
      }
    }

  } else {
    EH(-1, "You requested an invariant that I don't know how to calculate.");
  }

  // Return value of invariant
  return TI;

} // End of calc_tensor_invariants()

void supg_tau_shakib(SUPG_terms *supg_terms, int dim, double dt,
                     int interp_eqn) {
  double G[DIM][DIM];

  for (int i = 0; i < DIM; i++) {
    for (int j = 0; j < DIM; j++) {
      G[i][j] = 0;
      for (int k = 0; k < DIM; k++) {
        G[i][j] += bf[interp_eqn]->B[k][i] * bf[interp_eqn]->B[k][j];
      }
    }
  }

  double v_d_gv = 0;
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      v_d_gv += fv->v[i] * G[i][j] * fv->v[j];
    }
  }

  double d_v_d_gv[DIM][MDE];
  for (int a = 0; a < dim; a++) {
    for (int k = 0; k < ei[pg->imtrx]->dof[VELOCITY1]; k++) {
      d_v_d_gv[a][k] = 0.0;
      for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
          d_v_d_gv[a][k] +=
              delta(a, i) * bf[VELOCITY1 + a]->phi[k] * G[i][j] * fv->v[j] +
              delta(a, j) * fv->v[i] * G[i][j] * bf[VELOCITY1 + a]->phi[k];
        }
      }
    }
  }

  double beta = 2 / sqrt(15);
  supg_terms->supg_tau = beta / (sqrt(4 / (dt * dt) + v_d_gv));

  for (int a = 0; a < dim; a++) {
    for (int k = 0; k < ei[pg->imtrx]->dof[VELOCITY1]; k++) {
      supg_terms->d_supg_tau_dv[a][k] = -0.5 * d_v_d_gv[a][k] *
                                        (1 / (4 / (dt * dt) + v_d_gv)) *
                                        supg_terms->supg_tau;
    }
  }

  for (int a = 0; a < dim; a++) {
    if (pd->e[pg->imtrx][MESH_DISPLACEMENT1 + a]) {
      EH(-1,
         "Mesh displacement derivatives not implemented for shakib supg_tau");
    }
  }
}

void supg_tau_gauss_point(SUPG_terms *supg_terms, int dim, dbl diffusivity,
                          PG_DATA *pg_data) {
  double vnorm = 0;

  for (int i = 0; i < VIM; i++) {
    vnorm += fv->v[i] * fv->v[i];
  }
  vnorm = sqrt(vnorm);

  double hk = 0;
  for (int a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {
    hk += pg_data->hsquared[a];
  }
  /* This is the size of the element */
  hk = sqrt(hk / ((double)ei[pg->imtrx]->ielem_dim));

  double D = diffusivity;

  double hk_dX[DIM][MDE];
  for (int a = 0; a < dim; a++) {
    for (int j = 0; j < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1 + a]; j++) {
      double tmp = 0;
      for (int b = 0; b < dim; b++) {
        tmp += (2 * pg_data->hhv[b][a] * pg_data->dhv_dxnode[b][j]) /
               (2 * sqrt(pg_data->hsquared[b]));
      }
      hk_dX[a][j] = tmp / dim;
    }
  }

  double Pek = 0.5 * vnorm * hk / (D + DBL_EPSILON);

  double eta = Pek;
  double eta_dX[DIM][MDE];
  double eta_dV[DIM][MDE];
  if (Pek > 1) {
    eta = 1;
    for (int i = 0; i < DIM; i++) {
      for (int j = 0; j < MDE; j++) {
        eta_dX[i][j] = 0;
        eta_dV[i][j] = 0;
      }
    }
  } else {
    for (int i = 0; i < DIM; i++) {
      for (int j = 0; j < MDE; j++) {
        if (pd->e[pg->imtrx][VELOCITY1 + i]) {
          eta_dV[i][j] = 0.5 * 0.5 * hk * fv->v[i] * bf[VELOCITY1 + i]->phi[j] /
                         (vnorm * D);
        }

        if (pd->e[pg->imtrx][MESH_DISPLACEMENT1 + i]) {
          eta_dX[i][j] = 0.5 * vnorm * hk_dX[i][j] / D;
        }
      }
    }
  }

  if (vnorm > 0) {
    supg_terms->supg_tau = 0.5 * hk * eta / vnorm;

    for (int a = 0; a < VIM; a++) {
      int var = VELOCITY1 + a;
      for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        supg_terms->d_supg_tau_dv[a][j] = 0.5 * hk * eta * fv->v[a] *
                                              bf[var]->phi[j] /
                                              (-vnorm * vnorm * vnorm) +
                                          0.5 * hk * eta_dV[a][j] / vnorm;
      }

      var = MESH_DISPLACEMENT1 + a;
      for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        supg_terms->d_supg_tau_dX[a][j] =
            0.5 * hk_dX[a][j] * eta / vnorm + 0.5 * hk * eta_dX[a][j] / vnorm;
      }
    }

  } else {
    supg_terms->supg_tau = 0;
    for (int i = 0; i < DIM; i++) {
      for (int j = 0; j < MDE; j++) {
        supg_terms->d_supg_tau_dv[i][j] = 0.0;
      }
      for (int j = 0; j < MDE; j++) {
        supg_terms->d_supg_tau_dX[i][j] = 0.0;
      }
    }
  }
}

void supg_tau(SUPG_terms *supg_terms, int dim, dbl diffusivity,
              PG_DATA *pg_data, double dt, int shakib, int interp_eqn) {
  if (shakib) {
    supg_tau_shakib(supg_terms, dim, dt, interp_eqn);
  } else {
    supg_tau_gauss_point(supg_terms, dim, diffusivity, pg_data);
  }
}

static dbl yzbeta1(dbl scale, int dim, dbl Y, dbl Z, dbl d_Z[MDE], dbl beta,
                       dbl u, dbl d_u[MDE], dbl grad_u[DIM],
                       dbl d_grad_u[MDE][DIM], dbl h_elem, int interp_eqn,
                       dbl deriv[MDE]) {

//  static const dbl EPSILON = 1e-10;
  dbl Y_inv = 1.0 / Y;
//  dbl resid_scale = Y_inv * Z + EPSILON;
  dbl inner = 0;
  for (int i = 0; i < dim; i++) {
    inner += Y_inv * Y_inv * grad_u[i] * grad_u[i];
  }
  for (int i = 0; i < MDE; i++) {
    deriv[i] = 0;
  }
  return 1 * fabs(Z) * (1.0/(sqrt(inner+1e-12))) * h_elem * 0.5;

}

static dbl yzbeta2(dbl scale, dbl Y, dbl Z, dbl d_Z[MDE], dbl deriv[MDE], dbl h_elem, int interp_eqn)
{
//  static const dbl EPSILON = 1e-10;
  for (int k = 0; k < ei[pg->imtrx]->dof[interp_eqn]; k++) {
    deriv[k] = 0;
  }
  dbl yzbeta = 1.0 * fabs(Z) * h_elem * h_elem * 0.025;
  return yzbeta;
}

dbl yzbeta(dbl scale, int dim, dbl Y, dbl Z, dbl d_Z[MDE], dbl beta,
                       dbl u, dbl d_u[MDE], dbl grad_u[DIM],
                       dbl d_grad_u[MDE][DIM], dbl h_elem, int interp_eqn,
                       dbl deriv[MDE]) {
  dbl Y_inv = 1.0 / Y;

  dbl gradunit[DIM];
  dbl grad_u_norm = 0;

  for (int i = 0; i < dim; i++) {
    grad_u_norm  += grad_u[i] * grad_u[i];
  }
  grad_u_norm = sqrt(grad_u_norm) + DBL_EPSILON;
  dbl inv_grad_u_norm = 1 / grad_u_norm;

  for (int j = 0; j < ei[pd->mi[interp_eqn]]->dof[interp_eqn]; j++) {
    dbl inner = 0;
    for (int i = 0; i < dim; i++) {
      inner += Y_inv * grad_u[i] * Y_inv * grad_u[i];
    }
    inner += DIFFUSION_EPSILON;

    dbl d_inner = 0;
    for (int i = 0; i < dim; i++) {
      d_inner += 2 * Y_inv * d_grad_u[j][i] * Y_inv * grad_u[i];
    }


//    dbl scalar_part = Y_inv * u + DIFFUSION_EPSILON;

//    dbl p = 1-beta;
//    dbl q = (beta/2.0 - 1);


    dbl d_grad_u_norm = 0;
    for (int i = 0; i < dim; i++) {
      d_grad_u_norm += bf[interp_eqn]->grad_phi[j][i] * grad_u[i];
    }
    d_grad_u_norm *= inv_grad_u_norm;

    dbl d_inv_grad_u_norm = - inv_grad_u_norm * inv_grad_u_norm * d_grad_u_norm;

    for (int i = 0; i < dim; i++) {
      gradunit[i] = grad_u[i] * inv_grad_u_norm;
    }

    dbl h_dc = 0;
    for (int i = 0; i < ei[pd->mi[interp_eqn]]->dof[interp_eqn]; i++) {
      for (int j =0; j < dim; j++) {
        h_dc += fabs(gradunit[j] * d_grad_u[i][j]);
      }
    }

    //h_dc = 2 / h_dc;

//    dbl d_h_dc = 0;


    //deriv[j] = (Z+DIFFUSION_EPSILON) / (fabs(Z+DIFFUSION_EPSILON)) * d_Z[j] * pow(inner,q);
    //deriv[j] += fabs(Z+DIFFUSION_EPSILON) * d_inner * q * pow(inner, q-1);

    deriv[j] = d_inv_grad_u_norm;
  }
  return inv_grad_u_norm;
}


dbl yzbeta_model(int model, dbl scale, dbl beta, int dim, dbl Y,
                             dbl Z, dbl d_Z[MDE], dbl u, dbl d_u[MDE],
                             dbl grad_u[DIM], dbl d_grad_u[MDE][DIM],
                             dbl h_elem, int interp_eqn, dbl deriv[MDE]) {
  dbl dc = 0;
  switch (model) {
  case YZBETA_ONE:
    dc = yzbeta1(scale, dim, Y, Z, d_Z, 1.0, u, d_u, grad_u, d_grad_u,
                      h_elem, interp_eqn, deriv);
    break;
  case YZBETA_TWO:
    dc = yzbeta2(scale, Y, Z, d_Z, deriv, h_elem, interp_eqn);
    break;
  case YZBETA_MIXED: {
    dbl deriv1[MDE] = {0};
    dbl deriv2[MDE] = {0};
    dbl dc1, dc2;
    dc1 = yzbeta1(scale, dim, Y, Z, d_Z, 1.0, u, d_u, grad_u, d_grad_u,
                    h_elem, interp_eqn, deriv);
    dc2 = yzbeta2(scale, Y, Z, d_Z, deriv, h_elem, interp_eqn);
    for (int j = 0; j < ei[pd->mi[interp_eqn]]->dof[interp_eqn]; j++) {
      deriv[j] = 0.5 * (deriv1[j] + deriv2[j]);
    }
    dc = 0.5 * (dc1 + dc2);
  } break;
  case YZBETA_CUSTOM:
    dc = yzbeta(scale, dim, Y, Z, d_Z, beta, u, d_u, grad_u, d_grad_u,
                      h_elem, interp_eqn, deriv);
    break;
  default:
    EH(-1, "Unknown YZBETA Model");
    break;
  }

  return dc;
}
void get_supg_tau(SUPG_terms *supg_terms,
                  int dim,
                  dbl diffusivity,
                  PG_DATA *pg_data)
{
  double vnorm = 0;

  for (int i = 0; i < VIM; i++) {
    vnorm += fv->v[i]*fv->v[i];
  }
  vnorm = sqrt(vnorm);

  double hk = 0;
  for (int i = 0; i < dim; i++) {
    hk += sqrt(pg_data->hsquared[i]);
  }

  hk /= (double) dim;

  double D = diffusivity;

  double hk_dX[DIM][MDE];
  for (int a = 0; a < dim; a++)
    {
      for (int j = 0; j < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1+a]; j++)
        {
          double tmp = 0;
          for (int b = 0; b < dim; b++)
            {
              tmp += (2*pg_data->hhv[b][a] * pg_data->dhv_dxnode[b][j])/(2*sqrt(pg_data->hsquared[b]));
            }
          hk_dX[a][j] = tmp/dim;
        }
    }

  double Pek = 0.5 * vnorm * hk / D;

  double eta = Pek;
  double eta_dX[DIM][MDE];
  double eta_dV[DIM][MDE];
  if (Pek > 1) {
    eta = 1;
    for (int i = 0; i < DIM; i++)
    {
      for (int j = 0; j < MDE; j++)
      {
        eta_dX[i][j] = 0;
        eta_dV[i][j] = 0;
      }
    }
  }
  else
  {
    for (int i = 0; i < DIM; i++)
    {
      for (int j = 0; j < MDE; j++)
      {
        if (pd->e[VELOCITY1+i])
        {
          eta_dV[i][j] = 0.5 * 0.5 * hk * fv->v[i]*bf[VELOCITY1+i]->phi[j] / (vnorm*D);

        }

        if (pd->e[MESH_DISPLACEMENT1+i])
        {
          eta_dX[i][j] = 0.5 * vnorm * hk_dX[i][j] / D;

        }
      }
    }
  }

  if (vnorm > 0) {
    supg_terms->supg_tau = 0.5 * hk * eta / vnorm;

    for (int a = 0; a < VIM; a++)
      {
        int var = VELOCITY1 + a;
        for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++)
          {
            supg_terms->d_supg_tau_dv[a][j] = 0.5*hk*eta*fv->v[a]*bf[var]->phi[j] /
                (- vnorm*vnorm*vnorm) + 0.5 * hk * eta_dV[a][j] / vnorm;
          }

        var = MESH_DISPLACEMENT1 + a;
        for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++)
          {
            supg_terms->d_supg_tau_dX[a][j] = 0.5 * hk_dX[a][j] * eta / vnorm + 0.5 * hk * eta_dX[a][j] / vnorm;
          }
      }


  } else {
    supg_terms->supg_tau = 0;
    for (int i = 0; i < DIM; i++)
      {
        for (int j = 0; j < MDE; j++)
          {
            supg_terms->d_supg_tau_dv[i][j] = 0.0;
          }
        for (int j = 0; j < MDE; j++)
          {
            supg_terms->d_supg_tau_dX[i][j] = 0.0;
          }
      }
  }

}
