
/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) Goma Developers, National Technology & Engineering        *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/

#include "ad/util.h"
#include "ad/structs.h"
#include "el_elm_info.h"
#include "el_geom.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_eh.h"
#include "mm_fill_stress.h"
#include "mm_mp.h"
#include "rf_fem.h"

// GLOBAL
AD_Field_Variables *ad_fv = NULL;

int ad_calc_shearrate(ADType &gammadot,             /* strain rate invariant */
                      ADType gamma_dot[DIM][DIM]) { /* strain rate tensor */
  gammadot = 0.;
  int vdim = VIM;
  if (pd->gv[FILM_HEIGHT])
    vdim = 3;
  /* get gamma_dot invariant for viscosity calculations */
  for (int a = 0; a < vdim; a++) {
    for (int b = 0; b < vdim; b++) {
      gammadot += gamma_dot[a][b] * gamma_dot[b][a];
    }
  }
  // if (pd->gv[FILM_HEIGHT]) {
  //   gammadot += -(ad_fv->grad_v[0][0] + ad_fv->grad_v[1][1])  * -(ad_fv->grad_v[0][0] +
  //   ad_fv->grad_v[1][1]);
  // }

  gammadot = sqrt(0.5 * fabs(gammadot) + 1e-14);
  return 0;
}
int ad_beer_belly(void) {
  if (ad_fv == NULL) {
    ad_fv = new AD_Field_Variables();
  }
  int status = 0, i, j, k, dim, pdim, mdof, index, node, si;
  int DeformingMesh, ShapeVar;
  struct Basis_Functions *MapBf;
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

  if ((si = in_list(pd->IntegrationMap, 0, Num_Interpolations, Unique_Interpolations)) == -1) {
    GOMA_EH(GOMA_ERROR, "Seems to be a problem finding the IntegrationMap interpolation.");
  }
  MapBf = bfd[si];

  mdof = ei[imtrx]->dof[ShapeVar];

  if (MapBf->interpolation == I_N1) {
    mdof = MapBf->shape_dof;
  }

  /*
   * For every type "t" of unique basis function used in this problem,
   * initialize appropriate arrays...
   */

  if (ei[imtrx]->elem_blk_id != elem_blk_id_save) {
    is_initialized = FALSE;
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
    ad_fv->x[i] = 0.;
  }

  /*
   * NOTE: pdim is the number of coordinates, which may differ from
   * the element dimension (dim), as for shell elements!
   */
  for (i = 0; i < pdim; i++) {
    if (DeformingMesh) {
      for (k = 0; k < ei[upd->matrix_index[R_MESH1]]->dof[R_MESH1]; k++) {
        node = ei[upd->matrix_index[R_MESH1]]->dof_list[R_MESH1][k];

        index = Proc_Elem_Connect[Proc_Connect_Ptr[ei[upd->matrix_index[R_MESH1]]->ielem] + node];

        ad_fv->x[i] += (Coor[i][index] + ADType(ad_fv->total_ad_variables,
                                                ad_fv->offset[R_MESH1 + i] + k, *esp->d[i][k])) *
                       bf[R_MESH1]->phi[k];
      }
    } else {
      for (k = 0; k < mdof; k++) {
        node = MapBf->interpolation == I_N1 ? k : ei[imtrx]->dof_list[ShapeVar][k];

        index = Proc_Elem_Connect[Proc_Connect_Ptr[ei[imtrx]->ielem] + node];

        ad_fv->x[i] += Coor[i][index] * MapBf->phi[k];
      }
    }
  }

  /*
   * Elemental Jacobian is now affected by mesh displacement of nodes
   * from their initial nodal point coordinates...
   */

  for (i = 0; i < dim; i++) {
    for (j = 0; j < pdim; j++) {
      ad_fv->J[i][j] = 0.0;
      if (DeformingMesh) {
        for (k = 0; k < ei[upd->matrix_index[R_MESH1]]->dof[R_MESH1]; k++) {
          node = ei[upd->matrix_index[R_MESH1]]->dof_list[R_MESH1][k];

          index = Proc_Elem_Connect[Proc_Connect_Ptr[ei[upd->matrix_index[R_MESH1]]->ielem] + node];

          ad_fv->J[i][j] +=
              (Coor[j][index] +
               ADType(ad_fv->total_ad_variables, ad_fv->offset[R_MESH1 + j] + k, *esp->d[j][k])) *
              bf[R_MESH1]->dphidxi[k][i];
        }
      } else {
        for (k = 0; k < mdof; k++) {
          node = MapBf->interpolation == I_N1 ? k : ei[imtrx]->dof_list[ShapeVar][k];
          index = Proc_Elem_Connect[Proc_Connect_Ptr[ei[imtrx]->ielem] + node];
          ad_fv->J[i][j] += Coor[j][index] * bf[ShapeVar]->dphidxi[k][i];
        }
      }
    }
  }

  if (elem_shape == SHELL || elem_shape == TRISHELL || (mp->ehl_integration_kind == SIK_S)) {
    dim++;
    for (j = 0; j < pdim; j++) {
      ad_fv->J[pd->Num_Dim - 1][j] = MapBf->J[pd->Num_Dim - 1][j] = (j + 1) * 1.0;
    }

    /*Real Quick check on Jacobian to make sure this arbitrary assignment
     *didn't screw things up. Note that the detJ in the shell case can be
     *negative, but it is important to point out that we are not using it for
     *for integration, but only as a crutch for inversion of J */
    if (pd->Num_Dim == 3) {
      ad_fv->detJ =
          ad_fv->J[0][0] * (ad_fv->J[1][1] * ad_fv->J[2][2] - ad_fv->J[1][2] * ad_fv->J[2][1]) -
          ad_fv->J[0][1] * (ad_fv->J[1][0] * ad_fv->J[2][2] - ad_fv->J[2][0] * ad_fv->J[1][2]) +
          ad_fv->J[0][2] * (ad_fv->J[1][0] * ad_fv->J[2][1] - ad_fv->J[2][0] * ad_fv->J[1][1]);
    }
    if (pd->Num_Dim == 2) {
      ad_fv->detJ = ad_fv->J[0][0] * ad_fv->J[1][1] - ad_fv->J[0][1] * ad_fv->J[1][0];
    }

    if (fabs(ad_fv->detJ) < 1.e-10) {
      zero_detJ = TRUE;
#ifdef PARALLEL
      fprintf(stderr, "\nP_%d: Uh-oh, detJ =  %e\n", ProcID, fabs(ad_fv->detJ.val()));
#else
      fprintf(stderr, "\n Uh-oh, detJ =  %e\n", fabs(ad_fv->detJ));
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
    GOMA_EH(GOMA_ERROR, "dim = 1 not implemented in ad_beer_belly");
    break;

  case 2:
    dim = ei[pg->imtrx]->ielem_dim;
    ad_fv->detJ = ad_fv->J[0][0] * ad_fv->J[1][1] - ad_fv->J[0][1] * ad_fv->J[1][0];

    ad_fv->B[0][0] = ad_fv->J[1][1] / ad_fv->detJ;
    ad_fv->B[0][1] = -ad_fv->J[0][1] / ad_fv->detJ;
    ad_fv->B[1][0] = -ad_fv->J[1][0] / ad_fv->detJ;
    ad_fv->B[1][1] = ad_fv->J[0][0] / ad_fv->detJ;

    break;

  case 3:

    /* Now that we are here, reset dim for the shell case */
    dim = ei[imtrx]->ielem_dim;

    ad_fv->detJ =
        ad_fv->J[0][0] * (ad_fv->J[1][1] * ad_fv->J[2][2] - ad_fv->J[1][2] * ad_fv->J[2][1]) -
        ad_fv->J[0][1] * (ad_fv->J[1][0] * ad_fv->J[2][2] - ad_fv->J[2][0] * ad_fv->J[1][2]) +
        ad_fv->J[0][2] * (ad_fv->J[1][0] * ad_fv->J[2][1] - ad_fv->J[2][0] * ad_fv->J[1][1]);

    ad_fv->B[0][0] =
        (ad_fv->J[1][1] * ad_fv->J[2][2] - ad_fv->J[2][1] * ad_fv->J[1][2]) / (ad_fv->detJ);

    ad_fv->B[0][1] =
        -(ad_fv->J[0][1] * ad_fv->J[2][2] - ad_fv->J[2][1] * ad_fv->J[0][2]) / (ad_fv->detJ);

    ad_fv->B[0][2] =
        (ad_fv->J[0][1] * ad_fv->J[1][2] - ad_fv->J[1][1] * ad_fv->J[0][2]) / (ad_fv->detJ);

    ad_fv->B[1][0] =
        -(ad_fv->J[1][0] * ad_fv->J[2][2] - ad_fv->J[2][0] * ad_fv->J[1][2]) / (ad_fv->detJ);

    ad_fv->B[1][1] =
        (ad_fv->J[0][0] * ad_fv->J[2][2] - ad_fv->J[2][0] * ad_fv->J[0][2]) / (ad_fv->detJ);

    ad_fv->B[1][2] =
        -(ad_fv->J[0][0] * ad_fv->J[1][2] - ad_fv->J[1][0] * ad_fv->J[0][2]) / (ad_fv->detJ);

    ad_fv->B[2][0] =
        (ad_fv->J[1][0] * ad_fv->J[2][1] - ad_fv->J[1][1] * ad_fv->J[2][0]) / (ad_fv->detJ);

    ad_fv->B[2][1] =
        -(ad_fv->J[0][0] * ad_fv->J[2][1] - ad_fv->J[2][0] * ad_fv->J[0][1]) / (ad_fv->detJ);

    ad_fv->B[2][2] =
        (ad_fv->J[0][0] * ad_fv->J[1][1] - ad_fv->J[1][0] * ad_fv->J[0][1]) / (ad_fv->detJ);

    break;

  default:
    GOMA_EH(GOMA_ERROR, "Bad dim.");
    break;
  }
  return (status);
}

int ad_load_bf_grad(void) {
  int i, a, p, dofs = 0, status;
  struct Basis_Functions *bfv;

#ifdef DO_NOT_UNROLL
  int WIM;

  if ((pd->CoordinateSystem == CARTESIAN) || (pd->CoordinateSystem == CYLINDRICAL)) {
    WIM = pd->Num_Dim;
  } else {
    WIM = VIM;
  }
#endif

  status = 0;

  /* zero array for initialization */
  /*  v_length = DIM*DIM*DIM*MDE;
      init_vec_value(zero_array, 0., v_length); */
  if (ad_fv->basis.empty()) {
    ad_fv->basis.resize(V_LAST);
  }

  for (int v = V_FIRST; v < V_LAST; v++) {
    if (pd->gv[v]) {

      bfv = bf[v];
      dofs = ei[upd->matrix_index[v]]->dof[v];
      if (bfv->interpolation == I_N1) {
        dofs = bfv->shape_dof;
      }

      /* initialize variables */
      /* memset(&(bfv->d_phi[0][0]),0,siz); */

      /*
       * First load up components of the *raw* derivative vector "d_phi"
       */
      switch (pd->Num_Dim) {
      case 1:
        for (i = 0; i < dofs; i++) {
          ad_fv->basis[v].d_phi[i][0] =
              (ad_fv->B[0][0] * bfv->dphidxi[i][0] + ad_fv->B[0][1] * bfv->dphidxi[i][1]);
          ad_fv->basis[v].d_phi[i][1] = 0.0;
          ad_fv->basis[v].d_phi[i][2] = 0.0;
        }
        break;
      case 2:
        for (i = 0; i < dofs; i++) {
          ad_fv->basis[v].d_phi[i][0] =
              (ad_fv->B[0][0] * bfv->dphidxi[i][0] + ad_fv->B[0][1] * bfv->dphidxi[i][1]);
          ad_fv->basis[v].d_phi[i][1] =
              (ad_fv->B[1][0] * bfv->dphidxi[i][0] + ad_fv->B[1][1] * bfv->dphidxi[i][1]);
          ad_fv->basis[v].d_phi[i][2] = 0.0;
        }
        break;
      case 3:
        for (i = 0; i < dofs; i++) {
          ad_fv->basis[v].d_phi[i][0] =
              (ad_fv->B[0][0] * bfv->dphidxi[i][0] + ad_fv->B[0][1] * bfv->dphidxi[i][1] +
               ad_fv->B[0][2] * bfv->dphidxi[i][2]);
          ad_fv->basis[v].d_phi[i][1] =
              (ad_fv->B[1][0] * bfv->dphidxi[i][0] + ad_fv->B[1][1] * bfv->dphidxi[i][1] +
               ad_fv->B[1][2] * bfv->dphidxi[i][2]);
          ad_fv->basis[v].d_phi[i][2] =
              (ad_fv->B[2][0] * bfv->dphidxi[i][0] + ad_fv->B[2][1] * bfv->dphidxi[i][1] +
               ad_fv->B[2][2] * bfv->dphidxi[i][2]);
        }
        break;
      default:
        GOMA_EH(GOMA_ERROR, "Unexpected Dimension");
        break;
      }

      /*
       * Now, patch up the physical space gradient of this prototype
       * scalar function so scale factors are included.
       */

      /*	memset(&(bfv->grad_phi[0][0]),0,size1);  */

      for (i = 0; i < dofs; i++) {
        for (p = 0; p < WIM; p++) {
          ad_fv->basis[v].grad_phi[i][p] = (ad_fv->basis[v].d_phi[i][p]) / (fv->h[p]);
        }
      }

      for (i = 0; i < dofs; i++) {
        for (p = 0; p < VIM; p++) {
          for (a = 0; a < VIM; a++) {
            for (int q = 0; q < VIM; q++) {
              if (q == a)
                ad_fv->basis[v].grad_phi_e[i][a][p][a] = ad_fv->basis[v].grad_phi[i][p];
              else
                ad_fv->basis[v].grad_phi_e[i][a][p][q] = 0.0;
            }
          }
        }

        /* } */

        if (pd->CoordinateSystem != CARTESIAN) {
          GOMA_EH(GOMA_ERROR, "Only Cartesian coordinate system is supported, ad_load_bf_grad");
        }
      }
    } /* end of if v */
  } /* end of basis function loop. */

  return (status);
}

static inline ADType set_ad_or_dbl(dbl val, int eqn, int dof) {
  ADType tmp;
  if (ad_fv->total_ad_variables > 0 && af->Assemble_Jacobian == TRUE) {
    if (pd->gv[eqn]) {
      tmp = ADType(ad_fv->total_ad_variables, ad_fv->offset[eqn] + dof, val);
    } else {
      tmp = val;
    }
  } else {
    tmp = val;
  }
  return tmp;
}

extern "C" void fill_ad_field_variables() {
  if (ad_fv == NULL) {
    ad_fv = new AD_Field_Variables();
  }
  ad_fv->ielem = ei[pg->imtrx]->ielem;
  int num_ad_variables = 0;
  for (int i = V_FIRST; i < V_LAST; i++) {
    ad_fv->offset[i] = 0;
    // if (af->Assemble_Jacobian == TRUE) {
    if (pd->gv[i]) {
      ad_fv->offset[i] = num_ad_variables;
      num_ad_variables += ei[upd->matrix_index[i]]->dof[i];
    }
    // }
  }

  ad_fv->total_ad_variables = num_ad_variables;

  ad_beer_belly();
  ad_load_bf_grad();
  for (int p = 0; p < WIM; p++) {
    ad_fv->x_dot[p] = 0;
  }
  if (pd->gv[R_MESH1]) {
    for (int p = 0; p < WIM; p++) {
      for (int i = 0; i < ei[upd->matrix_index[R_MESH1 + p]]->dof[R_MESH1 + p]; i++) {
        ad_fv->d[p] += set_ad_or_dbl(*esp->d[p][i], R_MESH1 + p, i) * bf[R_MESH1 + p]->phi[i];
        if (pd->TimeIntegration != STEADY) {
          ADType udot = set_ad_or_dbl(*esp_dot->d[p][i], R_MESH1 + p, i);
          if (af->Assemble_Jacobian == TRUE) {
            udot.fastAccessDx(ad_fv->offset[R_MESH1 + p] + i) =
                (1. + 2. * tran->current_theta) / tran->delta_t;
          }
          ad_fv->x_dot[p] += udot * bf[R_MESH1 + p]->phi[i];
        } else {
          ad_fv->x_dot[p] = 0;
        }
      }
    }
  }

  if (pd->gv[VELOCITY1]) {
    for (int p = 0; p < WIM; p++) {
      ad_fv->v[p] = 0;
      ad_fv->v_dot[p] = 0;
      for (int i = 0; i < ei[upd->matrix_index[VELOCITY1 + p]]->dof[VELOCITY1 + p]; i++) {
        ad_fv->v[p] += set_ad_or_dbl(*esp->v[p][i], VELOCITY1 + p, i) * bf[VELOCITY1 + p]->phi[i];
        if (pd->TimeIntegration != STEADY) {
          ADType udot = set_ad_or_dbl(*esp_dot->v[p][i], VELOCITY1 + p, i);
          if (af->Assemble_Jacobian == TRUE) {
            udot.fastAccessDx(ad_fv->offset[VELOCITY1 + p] + i) =
                (1. + 2. * tran->current_theta) / tran->delta_t;
          }
          ad_fv->v_dot[p] += udot * bf[VELOCITY1 + p]->phi[i];
        } else {
          ad_fv->v_dot[p] = 0;
        }
      }
    }
    for (int p = 0; p < VIM; p++) {
      for (int q = 0; q < VIM; q++) {
        ad_fv->grad_v[p][q] = 0;

        for (int r = 0; r < WIM; r++) {
          for (int i = 0; i < ei[upd->matrix_index[VELOCITY1 + r]]->dof[VELOCITY1 + r]; i++) {
            ad_fv->grad_v[p][q] += set_ad_or_dbl(*esp->v[r][i], VELOCITY1 + r, i) *
                                   ad_fv->basis[VELOCITY1 + r].grad_phi_e[i][r][p][q];
          }
        }
      }
    }
  }

  if (pd->gv[SHEAR_RATE]) {
    ad_fv->SH = 0;
    for (int i = 0; i < ei[upd->matrix_index[SHEAR_RATE]]->dof[SHEAR_RATE]; i++) {
      ad_fv->SH += ADType(num_ad_variables, ad_fv->offset[SHEAR_RATE] + i, *esp->SH[i]) *
                   bf[SHEAR_RATE]->phi[i];
    }

    for (int q = 0; q < pd->Num_Dim; q++) {
      ad_fv->grad_eddy_nu[q] = 0;

      for (int i = 0; i < ei[upd->matrix_index[SHEAR_RATE]]->dof[SHEAR_RATE]; i++) {
        ad_fv->grad_eddy_nu[q] +=
            ADType(num_ad_variables, ad_fv->offset[SHEAR_RATE] + i, *esp->SH[i]) *
            ad_fv->basis[SHEAR_RATE].grad_phi[i][q];
      }
    }
  }

  if (pd->gv[SHELL_SAT_1]) {
    ad_fv->sh_sat_1 = 0;
    for (int i = 0; i < ei[upd->matrix_index[SHELL_SAT_1]]->dof[SHELL_SAT_1]; i++) {
      ad_fv->sh_sat_1 +=
          ADType(num_ad_variables, ad_fv->offset[SHELL_SAT_1] + i, *esp->sh_sat_1[i]) *
          bf[SHELL_SAT_1]->phi[i];
    }

    for (int q = 0; q < pd->Num_Dim; q++) {
      ad_fv->grad_sh_sat_1[q] = 0;

      for (int i = 0; i < ei[upd->matrix_index[SHELL_SAT_1]]->dof[SHELL_SAT_1]; i++) {
        ad_fv->grad_sh_sat_1[q] +=
            ADType(num_ad_variables, ad_fv->offset[SHELL_SAT_1] + i, *esp->sh_sat_1[i]) *
            ad_fv->basis[SHELL_SAT_1].grad_phi[i][q];
      }
    }
  }

  if (pd->gv[SHELL_SAT_2]) {
    ad_fv->sh_sat_2 = 0;
    for (int i = 0; i < ei[upd->matrix_index[SHELL_SAT_2]]->dof[SHELL_SAT_2]; i++) {
      ad_fv->sh_sat_2 +=
          ADType(num_ad_variables, ad_fv->offset[SHELL_SAT_2] + i, *esp->sh_sat_2[i]) *
          bf[SHELL_SAT_2]->phi[i];
    }

    for (int q = 0; q < pd->Num_Dim; q++) {
      ad_fv->grad_sh_sat_2[q] = 0;

      for (int i = 0; i < ei[upd->matrix_index[SHELL_SAT_2]]->dof[SHELL_SAT_2]; i++) {
        ad_fv->grad_sh_sat_2[q] +=
            ADType(num_ad_variables, ad_fv->offset[SHELL_SAT_2] + i, *esp->sh_sat_2[i]) *
            ad_fv->basis[SHELL_SAT_2].grad_phi[i][q];
      }
    }
  }

  if (pd->gv[SHELL_SAT_3]) {
    ad_fv->sh_sat_3 = 0;
    for (int i = 0; i < ei[upd->matrix_index[SHELL_SAT_3]]->dof[SHELL_SAT_3]; i++) {
      ad_fv->sh_sat_3 +=
          ADType(num_ad_variables, ad_fv->offset[SHELL_SAT_3] + i, *esp->sh_sat_3[i]) *
          bf[SHELL_SAT_3]->phi[i];
    }

    for (int q = 0; q < pd->Num_Dim; q++) {
      ad_fv->grad_sh_sat_3[q] = 0;

      for (int i = 0; i < ei[upd->matrix_index[SHELL_SAT_3]]->dof[SHELL_SAT_3]; i++) {
        ad_fv->grad_sh_sat_3[q] +=
            ADType(num_ad_variables, ad_fv->offset[SHELL_SAT_3] + i, *esp->sh_sat_3[i]) *
            ad_fv->basis[SHELL_SAT_3].grad_phi[i][q];
      }
    }
  }

  if (pd->gv[EDDY_NU]) {
    ad_fv->eddy_nu = 0;
    ad_fv->eddy_nu_dot = 0;
    for (int i = 0; i < ei[upd->matrix_index[EDDY_NU]]->dof[EDDY_NU]; i++) {
      ad_fv->eddy_nu += ADType(num_ad_variables, ad_fv->offset[EDDY_NU] + i, *esp->eddy_nu[i]) *
                        bf[EDDY_NU]->phi[i];

      if (pd->TimeIntegration != STEADY) {
        ADType ednudot = ADType(num_ad_variables, ad_fv->offset[EDDY_NU] + i, *esp_dot->eddy_nu[i]);
        ednudot.fastAccessDx(ad_fv->offset[EDDY_NU] + i) =
            (1. + 2. * tran->current_theta) / tran->delta_t;
        ad_fv->eddy_nu_dot += ednudot * bf[EDDY_NU]->phi[i];
      } else {
        ad_fv->eddy_nu_dot = 0;
      }
    }

    for (int q = 0; q < pd->Num_Dim; q++) {
      ad_fv->grad_eddy_nu[q] = 0;

      for (int i = 0; i < ei[upd->matrix_index[EDDY_NU]]->dof[EDDY_NU]; i++) {
        ad_fv->grad_eddy_nu[q] +=
            ADType(num_ad_variables, ad_fv->offset[EDDY_NU] + i, *esp->eddy_nu[i]) *
            ad_fv->basis[EDDY_NU].grad_phi[i][q];
      }
    }
  }

  if (pd->gv[TURB_K]) {
    ad_fv->turb_k = 0;
    ad_fv->turb_k_dot = 0;
    for (int i = 0; i < ei[upd->matrix_index[TURB_K]]->dof[TURB_K]; i++) {
      ad_fv->turb_k +=
          ADType(num_ad_variables, ad_fv->offset[TURB_K] + i, *esp->turb_k[i]) * bf[TURB_K]->phi[i];

      if (pd->TimeIntegration != STEADY) {
        ADType ednudot = ADType(num_ad_variables, ad_fv->offset[TURB_K] + i, *esp_dot->turb_k[i]);
        ednudot.fastAccessDx(ad_fv->offset[TURB_K] + i) =
            (1. + 2. * tran->current_theta) / tran->delta_t;
        ad_fv->turb_k_dot += ednudot * bf[TURB_K]->phi[i];
      } else {
        ad_fv->turb_k_dot = 0;
      }
    }

    for (int q = 0; q < pd->Num_Dim; q++) {
      ad_fv->grad_turb_k[q] = 0;

      for (int i = 0; i < ei[upd->matrix_index[TURB_K]]->dof[TURB_K]; i++) {
        ad_fv->grad_turb_k[q] +=
            ADType(num_ad_variables, ad_fv->offset[TURB_K] + i, *esp->turb_k[i]) *
            ad_fv->basis[TURB_K].grad_phi[i][q];
      }
    }
  }

  if (pd->gv[TURB_OMEGA]) {
    ad_fv->turb_omega = 0;
    ad_fv->turb_omega_dot = 0;
    for (int i = 0; i < ei[upd->matrix_index[TURB_OMEGA]]->dof[TURB_OMEGA]; i++) {
      ad_fv->turb_omega +=
          ADType(num_ad_variables, ad_fv->offset[TURB_OMEGA] + i, *esp->turb_omega[i]) *
          bf[TURB_OMEGA]->phi[i];

      if (pd->TimeIntegration != STEADY) {
        ADType ednudot =
            ADType(num_ad_variables, ad_fv->offset[TURB_OMEGA] + i, *esp_dot->turb_omega[i]);
        ednudot.fastAccessDx(ad_fv->offset[TURB_OMEGA] + i) =
            (1. + 2. * tran->current_theta) / tran->delta_t;
        ad_fv->turb_omega_dot += ednudot * bf[TURB_OMEGA]->phi[i];
      } else {
        ad_fv->turb_omega_dot = 0;
      }
    }

    for (int q = 0; q < pd->Num_Dim; q++) {
      ad_fv->grad_turb_omega[q] = 0;

      for (int i = 0; i < ei[upd->matrix_index[TURB_OMEGA]]->dof[TURB_OMEGA]; i++) {
        ad_fv->grad_turb_omega[q] +=
            ADType(num_ad_variables, ad_fv->offset[TURB_OMEGA] + i, *esp->turb_omega[i]) *
            ad_fv->basis[TURB_OMEGA].grad_phi[i][q];
      }
    }
  }
  if (pd->gv[FILM_HEIGHT]) {
    ad_fv->film_height = 0;
    ad_fv->film_height_dot = 0;
    for (int i = 0; i < ei[upd->matrix_index[FILM_HEIGHT]]->dof[FILM_HEIGHT]; i++) {
      ad_fv->film_height +=
          ADType(num_ad_variables, ad_fv->offset[FILM_HEIGHT] + i, *esp->film_height[i]) *
          bf[FILM_HEIGHT]->phi[i];

      if (pd->TimeIntegration != STEADY) {
        ADType ednudot =
            ADType(num_ad_variables, ad_fv->offset[FILM_HEIGHT] + i, *esp_dot->film_height[i]);
        ednudot.fastAccessDx(ad_fv->offset[FILM_HEIGHT] + i) =
            (1. + 2. * tran->current_theta) / tran->delta_t;
        ad_fv->film_height_dot += ednudot * bf[FILM_HEIGHT]->phi[i];
      } else {
        ad_fv->film_height_dot = 0;
      }
    }

    for (int q = 0; q < pd->Num_Dim; q++) {
      ad_fv->grad_film_height[q] = 0;

      for (int i = 0; i < ei[upd->matrix_index[FILM_HEIGHT]]->dof[FILM_HEIGHT]; i++) {
        ad_fv->grad_film_height[q] +=
            ADType(num_ad_variables, ad_fv->offset[FILM_HEIGHT] + i, *esp->film_height[i]) *
            ad_fv->basis[FILM_HEIGHT].grad_phi[i][q];
      }
    }
  }

  if (pd->gv[PRESSURE]) {
    ad_fv->P = 0;
    for (int i = 0; i < ei[upd->matrix_index[PRESSURE]]->dof[PRESSURE]; i++) {
      ad_fv->P += set_ad_or_dbl(*esp->P[i], PRESSURE, i) * bf[PRESSURE]->phi[i];
    }
    for (int q = 0; q < pd->Num_Dim; q++) {
      ad_fv->grad_P[q] = 0;

      for (int i = 0; i < ei[upd->matrix_index[PRESSURE]]->dof[PRESSURE]; i++) {
        ad_fv->grad_P[q] +=
            set_ad_or_dbl(*esp->P[i], PRESSURE, i) * ad_fv->basis[PRESSURE].grad_phi[i][q];
      }
    }
  }

  if (pd->gv[TEMPERATURE]) {
    ad_fv->T = 0;
    for (int i = 0; i < ei[upd->matrix_index[TEMPERATURE]]->dof[TEMPERATURE]; i++) {
      ad_fv->T += set_ad_or_dbl(*esp->T[i], TEMPERATURE, i) * bf[TEMPERATURE]->phi[i];
    }
  }

  if (pd->gv[POLYMER_STRESS11]) {
    int v_s[MAX_MODES][DIM][DIM];
    stress_eqn_pointer(v_s);
    int sdim = VIM;
    if (pd->gv[FILM_HEIGHT])
      sdim = 3;
    for (int mode = 0; mode < vn->modes; mode++) {
      for (int p = 0; p < sdim; p++) {
        for (int q = 0; q < sdim; q++) {
          ad_fv->S[mode][p][q] = 0;
          ad_fv->S_dot[mode][p][q] = 0;
        }
      }
    }
    for (int mode = 0; mode < vn->modes; mode++) {
      for (int p = 0; p < sdim; p++) {
        for (int q = 0; q < sdim; q++) {
          if (p <= q) {
            int v = v_s[mode][p][q];
            if (pd->gv[v]) {
              int dofs = ei[upd->matrix_index[v]]->dof[v];
              for (int i = 0; i < dofs; i++) {
                ad_fv->S[mode][p][q] +=
                    ADType(num_ad_variables, ad_fv->offset[v] + i, *esp->S[mode][p][q][i]) *
                    bf[v]->phi[i];
                if (pd->TimeIntegration != STEADY) {
                  ADType sdot =
                      ADType(num_ad_variables, ad_fv->offset[v] + i, *esp_dot->S[mode][p][q][i]);
                  sdot.fastAccessDx(ad_fv->offset[v] + i) =
                      (1. + 2. * tran->current_theta) / tran->delta_t;
                  ad_fv->S_dot[mode][p][q] += sdot * bf[v]->phi[i];
                } else {
                  ad_fv->S_dot[mode][p][q] = 0;
                }
              }
            }
            /* form the entire symmetric stress matrix for the momentum equation */
            ad_fv->S[mode][q][p] = ad_fv->S[mode][p][q];
            ad_fv->S_dot[mode][q][p] = ad_fv->S_dot[mode][p][q];
          }
          for (int r = 0; r < sdim; r++) {
            ad_fv->grad_S[mode][r][p][q] = 0.;
            int v = v_s[mode][p][q];
            if (pd->gv[v]) {
              int dofs = ei[upd->matrix_index[v]]->dof[v];

              for (int i = 0; i < dofs; i++) {
                if (p <= q) {
                  ad_fv->grad_S[mode][r][p][q] +=
                      ADType(num_ad_variables, ad_fv->offset[v] + i, *esp->S[mode][p][q][i]) *
                      ad_fv->basis[v].grad_phi[i][r];
                } else {
                  ad_fv->grad_S[mode][r][p][q] +=
                      ADType(num_ad_variables, ad_fv->offset[v] + i, *esp->S[mode][q][p][i]) *
                      ad_fv->basis[v].grad_phi[i][r];
                }
              }
            }
          }
        }
      }
      for (int r = 0; r < pd->Num_Dim; r++) {
        ad_fv->div_S[mode][r] = 0.0;

        for (int q = 0; q < pd->Num_Dim; q++) {
          ad_fv->div_S[mode][r] += ad_fv->grad_S[mode][q][q][r];
        }
      }
    }
  }
  for (int p = 0; pd->gv[VELOCITY_GRADIENT11] && p < 3; p++) {
    int v_g[DIM][DIM];
    v_g[0][0] = VELOCITY_GRADIENT11;
    v_g[0][1] = VELOCITY_GRADIENT12;
    v_g[1][0] = VELOCITY_GRADIENT21;
    v_g[1][1] = VELOCITY_GRADIENT22;
    v_g[0][2] = VELOCITY_GRADIENT13;
    v_g[1][2] = VELOCITY_GRADIENT23;
    v_g[2][0] = VELOCITY_GRADIENT31;
    v_g[2][1] = VELOCITY_GRADIENT32;
    v_g[2][2] = VELOCITY_GRADIENT33;
    for (int q = 0; q < 3; q++) {
      int v = v_g[p][q];
      ad_fv->G[p][q] = 0;
      if (pd->gv[v]) {
        int dofs = ei[upd->matrix_index[v]]->dof[v];
        for (int i = 0; i < dofs; i++) {
          ad_fv->G[p][q] +=
              ADType(num_ad_variables, ad_fv->offset[v] + i, *esp->G[p][q][i]) * bf[v]->phi[i];
        }
      }
    }
    for (int p = 0; p < 3; p++) {
      for (int q = 0; q < 3; q++) {
        int v = v_g[p][q];
        for (int r = 0; r < VIM; r++) {
          ad_fv->grad_G[r][p][q] = 0.0;
          if (pd->gv[v]) {
            int dofs = ei[upd->matrix_index[v]]->dof[v];
            for (int i = 0; i < dofs; i++) {
              ad_fv->grad_G[r][p][q] +=
                  ADType(num_ad_variables, ad_fv->offset[v] + i, *esp->G[p][q][i]) *
                  bf[v]->grad_phi[i][r];
            }
          }
        }
      }
    }

    /*
     * div(G) - this is a vector!
     */
    for (int r = 0; r < pd->Num_Dim; r++) {
      ad_fv->div_G[r] = 0.0;
      for (int q = 0; q < pd->Num_Dim; q++) {
        ad_fv->div_G[r] += ad_fv->grad_G[q][q][r];
      }
    }
  }

  // if (ei[pg->imtrx]->ielem == 418) {
  //   printf("ad_fv->P = %.15f\n", ad_fv->P.val());
  // }

#if 0
  // check field variables
  for (int p = 0; p < VIM; p++) {
    if (fabs(ad_fv->v[p].val() - fv->v[p]) > 1e-14) {
      printf("diff in fv->v[%d] %.12f != %.12f\n", p, ad_fv->v[p].val(), fv->v[p]);
    }
    for (int q = 0; q < VIM; q++) {
      if (fabs(ad_fv->grad_v[p][q].val() - fv->grad_v[p][q]) > 1e-12) {
        printf("diff in fv->grad_v[%d][%d] %.12f != %.12f\n", p, q, ad_fv->grad_v[p][q].val(),
               fv->grad_v[p][q]);
      }
    }
  }
  if (fabs(ad_fv->eddy_nu.val() - fv->eddy_nu) > 1e-14) {
    printf("diff in fv->eddy_nu %.12f != %.12f\n", ad_fv->eddy_nu.val(), fv->eddy_nu);
  }
  if (fabs(ad_fv->eddy_nu_dot.val() - fv_dot->eddy_nu) > 1e-14) {
    printf("diff in fv->eddy_nu_dot %.12f != %.12f\n", ad_fv->eddy_nu_dot.val(), fv_dot->eddy_nu);
  }
  for (int p = 0; p < pd->Num_Dim; p++) {
    if (fabs(ad_fv->grad_eddy_nu[p].val() - fv->grad_eddy_nu[p]) > 1e-14) {
      printf("diff in fv->grad_eddy_nu[%d] %.12f != %.12f\n", p, ad_fv->grad_eddy_nu[p].val(),
             fv->grad_eddy_nu[p]);
    }
  }
#endif
}