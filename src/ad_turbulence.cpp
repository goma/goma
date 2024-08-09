
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

#include <cstddef>
#ifdef GOMA_ENABLE_SACADO
#include "Sacado.hpp"

#include "ad_turbulence.h"
extern "C" {
#include "mm_fill_stabilization.h"
/* GOMA include files */
#include "mm_eh.h"
#include "mm_fill_stress.h"
#define GOMA_AD_TURBULENCE_CPP
#include "density.h"
#include "el_elm.h"
#include "el_elm_info.h"
#include "el_geom.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_mp.h"
#include "mm_mp_structs.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "std.h"
}

extern ADType ad_viscosity(struct Generalized_Newtonian *gn_local, ADType gamma_dot[DIM][DIM]);
AD_Field_Variables *ad_fv = NULL;

ADType ad_sa_viscosity(struct Generalized_Newtonian *gn_local);

/*  _______________________________________________________________________  */

static int calc_vort_mag(ADType &vort_mag, ADType omega[DIM][DIM]) {

  vort_mag = 0.;
  /* get gamma_dot invariant for viscosity calculations */
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      vort_mag += omega[a][b] * omega[a][b];
    }
  }

  vort_mag = sqrt(0.5 * (vort_mag) + 1e-14);
  return 0;
}

ADType ad_dcdd(int dim, int eqn, const ADType grad_U[DIM], ADType gs_inner_dot[DIM]) {
  ADType tmp = 0.0;
  ADType s[DIM] = {0.0};
  ADType r[DIM] = {0.0};
  for (int w = 0; w < dim; w++) {
    tmp += (ad_fv->v[w] - ad_fv->x_dot[w]) * (ad_fv->v[w] - ad_fv->x_dot[w]);
  }
  tmp = 1.0 / (sqrt(tmp + 1e-32));
  for (int w = 0; w < dim; w++) {
    s[w] = (ad_fv->v[w] - ad_fv->x_dot[w]) * tmp;
  }
  ADType mags = 0;
  for (int w = 0; w < dim; w++) {
    mags += (grad_U[w] * grad_U[w]);
  }
  mags = 1.0 / (sqrt(mags + 1e-32));
  for (int w = 0; w < dim; w++) {
    r[w] = grad_U[w] * mags;
  }

  ADType he = 0.0;
  for (int q = 0; q < ei[pg->imtrx]->dof[eqn]; q++) {
    ADType tmp = 0;
    for (int w = 0; w < dim; w++) {
      tmp += ad_fv->basis[eqn].grad_phi[q][w] * ad_fv->basis[eqn].grad_phi[q][w];
    }
    he += 1.0 / sqrt(std::max(tmp, 1e-20));
  }

  tmp = 0;
  for (int q = 0; q < ei[pg->imtrx]->dof[eqn]; q++) {
    for (int w = 0; w < dim; w++) {
      tmp += fabs(r[w] * ad_fv->basis[eqn].grad_phi[q][w]);
    }
  }
  ADType hrgn = 1.0 / (tmp + 1e-14);

  ADType magv = 0.0;
  for (int q = 0; q < VIM; q++) {
    magv += ad_fv->v[q] * ad_fv->v[q];
  }
  magv = sqrt(magv + 1e-32);

  ADType tau_dcdd = 0.5 * he * (1.0 / (mags + 1e-16)) * hrgn * hrgn;
  // ADType tau_dcdd = he * (1.0 / mags) * hrgn * hrgn / lambda;
  // printf("%g %g ", supg_tau.val(), tau_dcdd.val());
  // tau_dcdd = 1 / sqrt(1.0 / (supg_tau * supg_tau + 1e-32) +
  //                     1.0 / (tau_dcdd * tau_dcdd + 1e-32));
  // printf("%g \n ", tau_dcdd.val());
  ADType ss[DIM][DIM] = {{0.0}};
  ADType rr[DIM][DIM] = {{0.0}};
  ADType rdots = 0.0;
  for (int w = 0; w < dim; w++) {
    for (int z = 0; z < dim; z++) {
      ss[w][z] = s[w] * s[z];
      rr[w][z] = r[w] * r[z];
    }
    rdots += r[w] * s[w];
  }

  ADType inner_tensor[DIM][DIM] = {{0.0}};
  for (int w = 0; w < dim; w++) {
    for (int z = 0; z < dim; z++) {
      inner_tensor[w][z] = rr[w][z] - rdots * rdots * ss[w][z];
    }
  }

  for (int w = 0; w < dim; w++) {
    ADType tmp = 0.;
    for (int z = 0; z < dim; z++) {
      tmp += grad_U[w] * inner_tensor[w][z];
    }
    gs_inner_dot[w] = tmp;
    // gs_inner_dot[w] = grad_s[w][ii][jj];
  }
  return tau_dcdd;
}

template <typename scalar>
static int ad_calc_sa_S(scalar &S,                /* strain rate invariant */
                        scalar omega[DIM][DIM]) { /* strain rate tensor */
  S = 0.;
  /* get gamma_dot invariant for viscosity calculations */
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      S += omega[a][b] * omega[a][b];
    }
  }

  S = sqrt(0.5 * (S) + 1e-16);

  return 0;
}

int ad_calc_shearrate(ADType &gammadot,             /* strain rate invariant */
                      ADType gamma_dot[DIM][DIM]) { /* strain rate tensor */
  gammadot = 0.;
  /* get gamma_dot invariant for viscosity calculations */
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      gammadot += gamma_dot[a][b] * gamma_dot[b][a];
    }
  }

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
    GOMA_EH(GOMA_ERROR, "This is not implemented for shell elements. ad_beer_belly");
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
  }   /* end of basis function loop. */

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

  if (pd->gv[POLYMER_STRESS11]) {
    int v_s[MAX_MODES][DIM][DIM];
    stress_eqn_pointer(v_s);
    for (int mode = 0; mode < vn->modes; mode++) {
      for (int p = 0; p < VIM; p++) {
        for (int q = 0; q < VIM; q++) {
          ad_fv->S[mode][p][q] = 0;
          ad_fv->S_dot[mode][p][q] = 0;
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
          for (int r = 0; r < VIM; r++) {
            ad_fv->grad_S[mode][r][p][q] = 0.;
            int v = v_s[mode][p][q];
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
      for (int r = 0; r < pd->Num_Dim; r++) {
        ad_fv->div_S[mode][r] = 0.0;

        for (int q = 0; q < pd->Num_Dim; q++) {
          ad_fv->div_S[mode][r] += ad_fv->grad_S[mode][q][q][r];
        }
      }
    }
  }
  for (int p = 0; pd->gv[VELOCITY_GRADIENT11] && p < VIM; p++) {
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
    for (int q = 0; q < VIM; q++) {
      int v = v_g[p][q];
      if (pd->gv[v]) {
        ad_fv->G[p][q] = 0;
        int dofs = ei[upd->matrix_index[v]]->dof[v];
        for (int i = 0; i < dofs; i++) {
          ad_fv->G[p][q] +=
              ADType(num_ad_variables, ad_fv->offset[v] + i, *esp->G[p][q][i]) * bf[v]->phi[i];
        }
      }
    }
    for (int p = 0; p < VIM; p++) {
      for (int q = 0; q < VIM; q++) {
        int v = v_g[p][q];
        for (int r = 0; r < VIM; r++) {
          ad_fv->grad_G[r][p][q] = 0.0;
          int dofs = ei[upd->matrix_index[v]]->dof[v];
          for (int i = 0; i < dofs; i++) {
            ad_fv->grad_G[r][p][q] +=
                ADType(num_ad_variables, ad_fv->offset[v] + i, *esp->G[p][q][i]) *
                bf[v]->grad_phi[i][r];
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
#if 1
void ad_supg_tau_shakib(ADType &supg_tau, int dim, dbl dt, ADType diffusivity, int interp_eqn) {
  dbl G[DIM][DIM];

  get_metric_tensor(bf[interp_eqn]->B, dim, ei[pg->imtrx]->ielem_type, G);

  ADType v_d_gv = 0;
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      v_d_gv += fabs(ad_fv->v[i] * G[i][j] * ad_fv->v[j]);
    }
  }

  ADType diff_g_g = 0;
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      diff_g_g += G[i][j] * G[i][j];
    }
  }
  diff_g_g *= 9 * diffusivity * diffusivity;

  if (dt > 0) {
    supg_tau = 1.0 / (sqrt(4 / (dt * dt) + v_d_gv + diff_g_g));
  } else {
    supg_tau = 1.0 / (sqrt(v_d_gv + diff_g_g) + 1e-14);
  }
}
#else
void ad_supg_tau_shakib(ADType &supg_tau, int dim, dbl dt, ADType diffusivity, int interp_eqn) {
  ADType h_e = 0;
  for (int i = 0; i < ei[upd->matrix_index[interp_eqn]]->dof[interp_eqn]; i++) {
    ADType tmp = 0;
    for (int j = 0; j < pd->Num_Dim; j++) {
      tmp += ad_fv->basis[interp_eqn].grad_phi[i][j] * ad_fv->basis[interp_eqn].grad_phi[i][j];
    }
    h_e += std::sqrt(tmp);
  }
  h_e = 1 / h_e;

  ADType h_ugn = 0;
  ADType vmag = 0;
  for (int j = 0; j < pd->Num_Dim; j++) {
    vmag += (ad_fv->v[j] - ad_fv->x_dot[j]) * (ad_fv->v[j] - ad_fv->x_dot[j]);
  }
  vmag = std::sqrt(vmag + 1e-32);
  for (int i = 0; i < ei[upd->matrix_index[interp_eqn]]->dof[interp_eqn]; i++) {
    ADType tmp = 0;
    for (int j = 0; j < pd->Num_Dim; j++) {
      tmp += std::abs((ad_fv->v[j] - ad_fv->x_dot[j]) * ad_fv->basis[interp_eqn].grad_phi[i][j]);
    }
    h_ugn += tmp;
  }
  if (vmag > 1e-14) {
    h_ugn = vmag / (h_ugn + 1e-16);
  } else {
    h_ugn = 0;
  }

  ADType u_e = 0;
  for (int i = 0; i < ei[upd->matrix_index[VELOCITY1]]->dof[VELOCITY1]; i++) {
    ADType tmp = 0;
    for (int j = 0; j < WIM; j++) {
      ADType u_i =
          ADType(ad_fv->total_ad_variables, ad_fv->offset[VELOCITY1 + j] + i, *esp->v[j][i]);
      ADType xdot_i = 0;
      if (pd->gv[R_MESH1 + j]) {
        xdot_i =
            ADType(ad_fv->total_ad_variables, ad_fv->offset[R_MESH1 + j] + i, *esp_dot->d[j][i]);
        xdot_i.fastAccessDx(ad_fv->offset[R_MESH1 + j] + i) =
            (1. + 2. * tran->current_theta) / tran->delta_t;
      }
      tmp += SQUARE(u_i - xdot_i);
    }
    u_e += sqrt(tmp + 1e-32);
  }
  u_e /= ei[upd->matrix_index[VELOCITY1]]->dof[VELOCITY1];

  if (dt > 0) {
    supg_tau = 1 / sqrt(4.0 / (dt * dt) + h_ugn * h_ugn / (u_e * u_e));
  } else {
    supg_tau = h_ugn / u_e;
  }
}
#endif

void ad_get_metric_tensor(ADType B[DIM][DIM], int dim, int element_type, ADType G[DIM][DIM]) {
  dbl adjustment[DIM][DIM] = {{0}};
  const dbl invroot3 = 0.577350269189626;
  const dbl tetscale = 0.629960524947437; // 0.5 * cubroot(2)

  switch (element_type) {
  case LINEAR_TRI:
    adjustment[0][0] = (invroot3) * 2;
    adjustment[0][1] = (invroot3) * -1;
    adjustment[1][0] = (invroot3) * -1;
    adjustment[1][1] = (invroot3) * 2;
    break;
  case LINEAR_TET:
    adjustment[0][0] = tetscale * 2;
    adjustment[0][1] = tetscale * 1;
    adjustment[0][2] = tetscale * 1;
    adjustment[1][0] = tetscale * 1;
    adjustment[1][1] = tetscale * 2;
    adjustment[1][2] = tetscale * 1;
    adjustment[2][0] = tetscale * 1;
    adjustment[2][1] = tetscale * 1;
    adjustment[2][2] = tetscale * 2;
    break;
  default:
    adjustment[0][0] = 1.0;
    adjustment[1][1] = 1.0;
    adjustment[2][2] = 1.0;
    break;
  }

  // G = B * adjustment * B^T where B = J^-1

  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      G[i][j] = 0;
      for (int k = 0; k < dim; k++) {
        for (int m = 0; m < dim; m++) {
          G[i][j] += B[i][k] * adjustment[k][m] * B[j][m];
        }
      }
    }
  }
}
#if 1
void ad_only_tau_momentum_shakib(ADType &tau, int dim, dbl dt, int pspg_scale) {
  ADType G[DIM][DIM];
  dbl inv_rho = 1.0;
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;
  ADType gamma[DIM][DIM];
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      gamma[i][j] = ad_fv->grad_v[i][j] + ad_fv->grad_v[j][i];
    }
  }

  if (pspg_scale) {
    dbl rho = density(d_rho, dt);
    if (rho > 0.0) {
      inv_rho = 1.0 / rho;
    }
  }

  ad_get_metric_tensor(ad_fv->B, dim, ei[pg->imtrx]->ielem_type, G);

  ADType v_d_gv = 0;
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      v_d_gv += fabs((ad_fv->v[i] - ad_fv->x_dot[i]) * G[i][j] * (ad_fv->v[j] - ad_fv->x_dot[j]));
    }
  }

  ADType mu = ad_viscosity(gn, gamma);
  for (int mode = 0; mode < vn->modes; mode++) {
    ADType mup = ad_viscosity(ve[mode]->gn, gamma);
    mu += mup;
  }

  ADType coeff = (12.0 * mu * mu);

  ADType diff_g_g = 0;
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      diff_g_g += coeff * G[i][j] * G[i][j];
    }
  }

  if (pd->TimeIntegration != STEADY) {
    tau = inv_rho / (sqrt(4 / (dt * dt) + v_d_gv + diff_g_g));
  } else {
    tau = inv_rho / (sqrt(v_d_gv + diff_g_g) + 1e-14);
  }
}
#else
void ad_only_tau_momentum_shakib(ADType &tau, int dim, dbl dt, int pspg_scale) {
  ADType G[DIM][DIM];
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;
  ADType gamma[DIM][DIM];
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      gamma[i][j] = ad_fv->grad_v[i][j] + ad_fv->grad_v[j][i];
    }
  }

  dbl rho = 1.;
  if (pspg_scale) {
    rho = density(d_rho, dt);
    if (rho < 1e-14) {
      rho = 1.;
    }
  }

  int interp_eqn = VELOCITY1;

  ADType mu = ad_viscosity(gn, gamma);
  for (int mode = 0; mode < vn->modes; mode++) {
    ADType mup = ad_viscosity(ve[mode]->gn, gamma);
    mu += mup;
  }

  ADType sugn1 = 0;
  for (int j = 0; j < ei[pg->imtrx]->dof[interp_eqn]; j++) {
    for (int i = 0; i < VIM; i++) {
      sugn1 += std::abs(ad_fv->v[i] * ad_fv->basis[interp_eqn].grad_phi[j][i]);
    }
  }
  sugn1 = 1.0 / std::max(sugn1, 1e-20);

  ADType sugn2 = 0;
  sugn2 = dt / 2;

  ADType r[DIM] = {0.};
  ADType v_norm;
  for (int i = 0; i < dim; i++) {
    v_norm = ad_fv->v[i] * ad_fv->v[i];
  }
  v_norm = sqrt(std::max(v_norm, 1e-20));

  ADType norm_grad_v_norm = 0;
  ADType grad_v_norm[DIM] = {0.};
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < VIM; j++) {
      grad_v_norm[i] += ad_fv->v[i] * ad_fv->grad_v[i][j];
    }
    grad_v_norm[i] /= v_norm;
    norm_grad_v_norm += grad_v_norm[i] * grad_v_norm[i];
  }

  norm_grad_v_norm = sqrt(std::max(norm_grad_v_norm, 1e-20));
  for (int i = 0; i < dim; i++) {
    r[i] = grad_v_norm[i] / norm_grad_v_norm;
  }

  ADType h_rgn = 0;
  for (int j = 0; j < ei[pg->imtrx]->dof[interp_eqn]; j++) {
    for (int i = 0; i < VIM; i++) {
      h_rgn += std::abs(r[i] * ad_fv->basis[interp_eqn].grad_phi[j][i]);
    }
  }
  h_rgn = 2.0 / std::max(h_rgn, 1e-20);

  ADType sugn3 = std::max(h_rgn * h_rgn * rho / (4 * mu), 1e-20);

  if (pd->TimeIntegration != STEADY) {
    tau = (1 / rho) * 1.0 / sqrt(1 / (sugn1 * sugn1) + 1 / (sugn2 * sugn2) + 1 / (sugn3 * sugn3));
  } else {
    tau = (1 / rho) * 1.0 / sqrt(1 / (sugn1 * sugn1) + 1 / (sugn3 * sugn3));
  }
}
#endif

extern "C" void
ad_tau_momentum_shakib(momentum_tau_terms *tau_terms, int dim, dbl dt, int pspg_scale) {
  dbl G[DIM][DIM];
  dbl inv_rho = 1.0;
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;

  if (pspg_scale) {
    dbl rho = density(d_rho, dt);
    inv_rho = 1.0 / rho;
  }

  int interp_eqn = VELOCITY1;
  get_metric_tensor(bf[interp_eqn]->B, dim, ei[pg->imtrx]->ielem_type, G);

  ADType v_d_gv = 0;
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      v_d_gv += fabs(ad_fv->v[i] * G[i][j] * ad_fv->v[j]);
    }
  }
  ADType gamma[DIM][DIM];
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      gamma[i][j] = ad_fv->grad_v[i][j] + ad_fv->grad_v[j][i];
    }
  }

  ADType mu = ad_viscosity(gn, gamma);

  ADType coeff = (12.0 * mu * mu);

  ADType diff_g_g = 0;
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      diff_g_g += coeff * G[i][j] * G[i][j];
    }
  }

  dbl d_v_d_gv[DIM][MDE];
  for (int a = 0; a < dim; a++) {
    for (int k = 0; k < ei[pg->imtrx]->dof[VELOCITY1]; k++) {
      d_v_d_gv[a][k] = 0.0;
      for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
          d_v_d_gv[a][k] += delta(a, i) * bf[VELOCITY1 + a]->phi[k] * G[i][j] * fv->v[j] +
                            delta(a, j) * fv->v[i] * G[i][j] * bf[VELOCITY1 + a]->phi[k];
        }
      }
    }
  }
  ADType tau;

  if (pd->TimeIntegration != STEADY) {
    tau = inv_rho / (sqrt(4 / (dt * dt) + v_d_gv + diff_g_g));
    tau_terms->tau = tau.val();
  } else {
    tau = inv_rho / (sqrt(v_d_gv + diff_g_g) + 1e-14);
    tau_terms->tau = tau.val();
  }

  for (int a = 0; a < dim; a++) {
    for (int k = 0; k < ei[pg->imtrx]->dof[VELOCITY1]; k++) {
      tau_terms->d_tau_dv[a][k] = tau.dx(ad_fv->offset[VELOCITY1 + a] + k);
    }
  }

#if 0
  if (pd->e[pg->imtrx][MESH_DISPLACEMENT1]) {
    dbl dG[DIM][DIM][DIM][MDE];
    get_metric_tensor_deriv(bf[MESH_DISPLACEMENT1]->B, bf[MESH_DISPLACEMENT1]->dB, dim,
                            MESH_DISPLACEMENT1, ei[pg->imtrx]->ielem_type, dG);
    for (int a = 0; a < dim; a++) {
      for (int k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1 + a]; k++) {
        dbl v_d_gv_dx = 0;
        for (int i = 0; i < dim; i++) {
          for (int j = 0; j < dim; j++) {
            v_d_gv_dx += fv->v[i] * dG[i][j][a][k] * fv->v[j];
          }
        }

        dbl diff_g_g_dx = 0;
        for (int i = 0; i < dim; i++) {
          for (int j = 0; j < dim; j++) {
            diff_g_g_dx += 2 * coeff * dG[i][j][a][k] * G[i][j];
          }
        }
        tau_terms->d_tau_dX[a][k] = inv_rho * -0.5 *
                                    (v_d_gv_dx + diff_g_g_dx + d_mu->X[a][k] * d_diff_g_g_dmu) *
                                    supg_tau_cubed;
      }
    }
  }
  if (pd->e[pg->imtrx][TEMPERATURE]) {
    for (int k = 0; k < ei[pg->imtrx]->dof[TEMPERATURE]; k++) {
      tau_terms->d_tau_dT[k] = inv_rho * -0.5 * (d_mu->T[k] * d_diff_g_g_dmu) * supg_tau_cubed;
    }
  }
  if (pd->e[pg->imtrx][PRESSURE]) {
    for (int k = 0; k < ei[pg->imtrx]->dof[PRESSURE]; k++) {
      tau_terms->d_tau_dP[k] = inv_rho * -0.5 * (d_mu->P[k] * d_diff_g_g_dmu) * supg_tau_cubed;
    }
  }
  if (pd->e[pg->imtrx][FILL]) {
    for (int k = 0; k < ei[pg->imtrx]->dof[FILL]; k++) {
      tau_terms->d_tau_dF[k] = inv_rho * -0.5 * (d_mu->F[k] * d_diff_g_g_dmu) * supg_tau_cubed;
    }
  }
  if (pd->e[pg->imtrx][BOND_EVOLUTION]) {
    for (int k = 0; k < ei[pg->imtrx]->dof[BOND_EVOLUTION]; k++) {
      tau_terms->d_tau_dnn[k] = inv_rho * -0.5 * (d_mu->nn[k] * d_diff_g_g_dmu) * supg_tau_cubed;
    }
  }
#endif
  if (pd->e[pg->imtrx][EDDY_NU]) {
    for (int k = 0; k < ei[pg->imtrx]->dof[EDDY_NU]; k++) {
      tau_terms->d_tau_dEDDY_NU[k] = tau.dx(ad_fv->offset[EDDY_NU] + k);
    }
  }
#if 0
  if (pd->e[pg->imtrx][TURB_K]) {
    for (int k = 0; k < ei[pg->imtrx]->dof[TURB_K]; k++) {
      tau_terms->d_tau_dturb_k[k] =
          inv_rho * -0.5 * (d_mu->turb_k[k] * d_diff_g_g_dmu) * supg_tau_cubed;
    }
  }
  if (pd->e[pg->imtrx][TURB_OMEGA]) {
    for (int k = 0; k < ei[pg->imtrx]->dof[TURB_OMEGA]; k++) {
      tau_terms->d_tau_dturb_omega[k] =
          inv_rho * -0.5 * (d_mu->turb_omega[k] * d_diff_g_g_dmu) * supg_tau_cubed;
    }
  }
  if (pd->e[pg->imtrx][MASS_FRACTION]) {
    for (int w = 0; w < pd->Num_Species_Eqn; w++) {
      for (int k = 0; k < ei[pg->imtrx]->dof[MASS_FRACTION]; k++) {
        tau_terms->d_tau_dC[w][k] =
            inv_rho * -0.5 * (d_mu->C[w][k] * d_diff_g_g_dmu) * supg_tau_cubed;
      }
    }
  }
#endif
}

extern "C" void ad_sa_wall_func(double func[DIM],
                                double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE]) {
  // kind of hacky, near wall velocity is velocity at central node
  ADType unw[DIM];
  ADType eddy_nw;
  unw[0] = ADType(ad_fv->total_ad_variables, ad_fv->offset[VELOCITY1] + 8, *esp->v[0][8]);
  unw[1] = ADType(ad_fv->total_ad_variables, ad_fv->offset[VELOCITY2] + 8, *esp->v[1][8]);
  eddy_nw = ADType(ad_fv->total_ad_variables, ad_fv->offset[EDDY_NU] + 8, *esp->eddy_nu[8]);

  ADType mu = 0;
  dbl scale = 1.0;
  DENSITY_DEPENDENCE_STRUCT d_rho;
  if (gn->ConstitutiveEquation == TURBULENT_SA_DYNAMIC) {
    scale = density(&d_rho, tran->time_value);
  }
  int negative_mu_e = FALSE;
  if (fv_old->eddy_nu < 0) {
    negative_mu_e = TRUE;
  }

  double mu_newt = mp->viscosity;
  ADType fv1 = 1.0;
  if (negative_mu_e) {
    mu = 0;
  } else {

    ADType mu_e = eddy_nw;
    ADType cv1 = 7.1;
    ADType chi = mu_e / mu_newt;
    fv1 = pow(chi, 3) / (pow(chi, 3) + pow(cv1, 3));

    mu = scale * (mu_e * fv1);
    if (mu > 1e3 * mu_newt) {
      mu = 1e3 * mu_newt;
    }
  }

  ADType ut = 0;
  for (int i = 0; i < 2; i++) {
    ut += unw[i] * fv->stangent[0][i];
  }

  ADType normgv = 0;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      normgv += ad_fv->grad_v[i][j] * ad_fv->grad_v[i][j];
    }
  }
  normgv = std::sqrt(normgv);

  ADType mu_t = scale * ut * ut / (normgv + 1e-12);

  ADType mu_tt = std::max(mu_t - mu, 0);

  ADType omega[DIM][DIM];
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      omega[a][b] = fv_old->grad_v[a][b] + fv_old->grad_v[b][a];
    }
  }

  ADType S = 0.;
  /* get gamma_dot invariant for viscosity calculations */
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      S += omega[a][b] * omega[a][b];
    }
  }

  S = sqrt(0.5 * (S) + 1e-16);

  ADType eddy_nu = std::min(std::min(5, 1e-1 * scale * S * mp->viscosity), fv_old->eddy_nu * 1.5);

  func[0] = fv->eddy_nu - eddy_nu.val();

  // printf("eddynu = %g %g %g\n", eddy_nu.val(), fv->eddy_nu, func[0]);

  // for (int p = 0; p < WIM; p++) {
  //   for (int i = 0; i < ei[pg->imtrx]->dof[VELOCITY1 + p]; i++) {
  //     d_func[0][VELOCITY1 + p][i] = eddy_nu.dx(ad_fv->v_offset[p] + i);
  //   }
  // }

  // for (int i = 0; i < ei[pg->imtrx]->dof[EDDY_NU]; i++) {
  //   d_func[0][EDDY_NU][i] = eddy_nu.dx(ad_fv->eddy_nu_offset + i);
  // }
}

extern "C" void ad_omega_wall_func(double func[DIM],
                                   double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE]) {

  ADType d = std::max(fv->wall_distance, 0.0004);
  dbl nu = 1.5e-5;
  dbl beta1 = 0.075;
  ADType r = ad_fv->turb_omega - 6 * nu / (beta1 * d * d);
  func[0] = r.val();
  for (int j = 0; j < ei[pg->imtrx]->dof[TURB_OMEGA]; j++) {
    d_func[0][TURB_OMEGA][j] = r.dx(ad_fv->offset[TURB_OMEGA] + j);
  }
}

ADType ad_sa_viscosity(struct Generalized_Newtonian *gn_local) {
  ADType mu = 0;
  dbl scale = 1.0;
  DENSITY_DEPENDENCE_STRUCT d_rho;
  if (gn_local->ConstitutiveEquation == TURBULENT_SA_DYNAMIC) {
    scale = density(&d_rho, tran->time_value);
  }
  int negative_mu_e = FALSE;
  if (fv_old->eddy_nu < 0) {
    negative_mu_e = TRUE;
  }

  double mu_newt = mp->viscosity;
  if (negative_mu_e) {
    mu = mu_newt;
  } else {

    ADType mu_e = ad_fv->eddy_nu;
    ADType cv1 = 7.1;
    ADType chi = mu_e / mu_newt;
    ADType fv1 = pow(chi, 3) / (pow(chi, 3) + pow(cv1, 3));

    mu = scale * (mu_newt + (mu_e * fv1));
    if (mu > 1e3 * mu_newt) {
      mu = 1e3 * mu_newt;
    }
  }

  return mu;
}

extern "C" dbl ad_sa_viscosity(struct Generalized_Newtonian *gn_local,
                               VISCOSITY_DEPENDENCE_STRUCT *d_mu) {
  ADType mu = 0;
  dbl scale = 1.0;
  DENSITY_DEPENDENCE_STRUCT d_rho;
  if (gn_local->ConstitutiveEquation == TURBULENT_SA_DYNAMIC) {
    scale = density(&d_rho, tran->time_value);
  }
  int negative_mu_e = FALSE;
  if (fv_old->eddy_nu < 0) {
    negative_mu_e = TRUE;
  }

  double mu_newt = mp->viscosity;
  if (negative_mu_e) {
    mu = mu_newt;
  } else {

    ADType mu_e = ad_fv->eddy_nu;
    ADType cv1 = 7.1;
    ADType chi = mu_e / mu_newt;
    ADType fv1 = pow(chi, 3) / (pow(chi, 3) + pow(cv1, 3));

    mu = scale * (mu_newt + (mu_e * fv1));
    if (mu > 1e3 * mu_newt) {
      mu = 1e3 * mu_newt;
    }

    if (d_mu != NULL) {
      for (int j = 0; j < ei[pg->imtrx]->dof[EDDY_NU]; j++) {
        d_mu->eddy_nu[j] = mu.dx(ad_fv->offset[EDDY_NU] + j);
      }
    }
  }

  return mu.val();
}

/* assemble_spalart_allmaras -- assemble terms (Residual & Jacobian) for conservation
 *                              of eddy viscosity for Spalart Allmaras turbulent flow model
 *
 *  Kessels, P. C. J. "Finite element discretization of the Spalart-Allmaras
 *  turbulence model." (2016).
 *
 *  Spalart, Philippe, and Steven Allmaras. "A one-equation turbulence model for
 *  aerodynamic flows." 30th aerospace sciences meeting and exhibit. 1992.
 *
 * in:
 *      time value
 *      theta (time stepping parameter, 0 for BE, 0.5 for CN)
 *      time step size
 *      Streamline Upwind Petrov Galerkin (PG) data structure
 *
 * out:
 *      lec -- gets loaded up with local contributions to resid, Jacobian
 *
 * Created:     August 2022 kristianto.tjiptowidjojo@averydennison.com
 * Modified:    June 2023 Weston Ortiz
 *
 */

extern "C" int ad_assemble_spalart_allmaras(dbl time_value, /* current time */
                                            dbl tt, /* parameter to vary time integration from
                                                       explicit (tt = 1) to implicit (tt = 0)    */
                                            dbl dt, /* current time step size                    */
                                            const PG_DATA *pg_data) {

  //! WIM is the length of the velocity vector
  int i, j, a, b;
  int eqn, var, peqn, pvar;
  int *pdv = pd->v[pg->imtrx];

  int status = 0;

  eqn = EDDY_NU;
  double d_area = fv->wt * bf[eqn]->detJ * fv->h3;

  /* Get Eddy viscosity at Gauss point */
  ADType mu_e = ad_fv->eddy_nu;

  int negative_sa = false;
  // Previous workaround, see comment for negative_Se below
  //
  //   int transient_run = FALSE;
  //   if (pd->TimeIntegration != STEADY) {
  //     transient_run = true;
  //   }
  //   // Use old values for equation switching for transient runs
  //   // Seems to work reasonably well.
  //   if (transient_run && (fv_old->eddy_nu < 0)) {
  //     negative_sa = true;
  //   } else if (!transient_run && (mu_e < 0)) {
  //     // Kris thinks it might work with switching equations in steady state
  //     negative_sa = true;
  //   }
  if (mu_e < 0) {
    negative_sa = true;
  }

  /* Get fluid viscosity */
  double mu_newt = mp->viscosity;

  /* Rate of rotation tensor  */
  ADType omega[DIM][DIM];
  double omega_old[DIM][DIM];
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      omega[a][b] = (ad_fv->grad_v[a][b] - ad_fv->grad_v[b][a]);
      omega_old[a][b] = (fv_old->grad_v[a][b] - fv_old->grad_v[b][a]);
    }
  }

  /* Vorticity */
  ADType S = 0.0;
  double S_old = 0;
  ad_calc_sa_S(S, omega);
  ad_calc_sa_S(S_old, omega_old);

  double d = fv->wall_distance;
  /* Get distance from nearest wall */
  if (d < 1.0e-6)
    d = 1.0e-6;

  /* Model coefficients (constants) */
  double cb1 = 0.1355;
  double cb2 = 0.622;
  double cv1 = 7.1;
  double cv2 = 0.7;
  double cv3 = 0.9;
  double sigma = (2.0 / 3.0);
  double cw2 = 0.3;
  double cw3 = 2.0;
  double cn1 = 16;
  double kappa = 0.41;
  double cw1 = (cb1 / kappa / kappa) + (1.0 + cb2) / sigma;

  /* More model coefficients (depends on mu_e) */
  ADType chi = mu_e / mu_newt;
  ADType fv1 = pow(chi, 3) / (pow(chi, 3) + pow(cv1, 3));
  ADType fv2 = 1.0 - (chi) / (1.0 + chi * fv1);
  ADType fn = 1.0;
  if (negative_sa) {
    fn = (cn1 + pow(chi, 3.0)) / (cn1 - pow(chi, 3));
  }
  ADType Sbar = (mu_e * fv2) / (kappa * kappa * d * d);
  int negative_Se = false;
  // I tried to use the old values for equation switching for transient runs but
  // I end up getting floating point errors. because of Se going very far below
  // zero I'm trying to use current values instead and hope that Newton's method
  // will converge with the switching equations
  // previously:
  // . double Sbar_old = (fv_old->eddy_nu * fv2) / (kappa * kappa * d * d);
  //   if (transient_run && (Sbar_old < -cv2 * S_old)) {
  //     negative_Se = true;
  //   } else if (!transient_run && (Sbar < -cv2 * S)) {
  // .   negative_Se = true;
  //   }
  if (Sbar < -cv2 * S) {
    negative_Se = true;
  }
  ADType S_e = S + Sbar;
  if (negative_Se) {
    S_e = S + S * (cv2 * cv2 * S + cv3 * Sbar) / ((cv3 - 2 * cv2) * S - Sbar);
  }
  double r_max = 10.0;
  ADType r = 0.0;
  if (fabs(S_e) > 1.0e-6) {
    r = mu_e / (kappa * kappa * d * d * S_e);
  } else {
    r = r_max;
  }
  if (r >= r_max) {
    r = r_max;
  }
  // Arbitrary limit to avoid floating point errors should only hit this when
  // S_e is very small and either mu_e or S_e are negative.  Which means we are
  // already trying to alleviate the issue.
  if (r < -100) {
    r = -100;
  }
  ADType g = r + cw2 * (pow(r, 6) - r);
  ADType fw_inside = (1.0 + pow(cw3, 6)) / (pow(g, 6) + pow(cw3, 6));
  ADType fw = g * pow(fw_inside, (1.0 / 6.0));

  dbl supg = 1.;
  ADType supg_tau = 0;
  if (mp->Mwt_funcModel == GALERKIN) {
    supg = 0.;
  } else if (mp->Mwt_funcModel == SUPG || mp->Mwt_funcModel == SUPG_GP ||
             mp->Mwt_funcModel == SUPG_SHAKIB) {
    supg = mp->Mwt_func;
    ad_supg_tau_shakib(supg_tau, pd->Num_Dim, dt, mu_newt, EDDY_NU);
  }

  /*
   * Residuals_________________________________________________________________
   */

  std::vector<ADType> resid(ei[pg->imtrx]->dof[eqn]);
  for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
    resid[i] = 0;
  }

  if (af->Assemble_Residual) {
    /*
     * Assemble residual for eddy viscosity
     */
    eqn = EDDY_NU;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ADType wt_func = bf[eqn]->phi[i];

      if (supg > 0) {
        if (supg != 0.0) {
          for (int p = 0; p < VIM; p++) {
            wt_func += supg * supg_tau * ad_fv->v[p] * bf[eqn]->grad_phi[i][p];
          }
        }
      }

      /* Assemble mass term */
      ADType mass = 0.0;
      if (pd->TimeIntegration != STEADY) {
        if (pd->e[pg->imtrx][eqn] & T_MASS) {
          mass += ad_fv->eddy_nu_dot * wt_func * d_area;
          mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
        }
      }

      /* Assemble advection term */
      ADType adv = 0;
      for (int p = 0; p < VIM; p++) {
        adv += ad_fv->v[p] * ad_fv->grad_eddy_nu[p];
      }
      adv *= wt_func * d_area;
      adv *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

      double neg_c = 1.0;
      if (negative_sa) {
        neg_c = -1.0;
      }

      /* Assemble source terms */
      ADType src_1 = cb1 * S_e * mu_e;
      ADType src_2 = neg_c * cw1 * fw * (mu_e * mu_e) / (d * d);
      ADType src = -src_1 + src_2;
      src *= wt_func * d_area;
      src *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

      /* Assemble diffusion terms */
      ADType diff_1 = 0.0;
      ADType diff_2 = 0.0;
      for (int p = 0; p < VIM; p++) {
        diff_1 += bf[eqn]->grad_phi[i][p] * (mu_newt + mu_e * fn) * ad_fv->grad_eddy_nu[p];
        diff_2 += wt_func * cb2 * ad_fv->grad_eddy_nu[p] * ad_fv->grad_eddy_nu[p];
      }
      ADType diff = (1.0 / sigma) * (diff_1 - diff_2);
      diff *= d_area;
      diff *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];

      resid[i] += mass + adv + src + diff;
      lec->R[LEC_R_INDEX(peqn, i)] += mass.val() + adv.val() + src.val() + diff.val();
    } /* end of for (i=0,ei[pg->imtrx]->dofs...) */
  }   /* end of if assemble residual */

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    eqn = EDDY_NU;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

      /* Sensitivity w.r.t. eddy viscosity */
      var = EDDY_NU;
      if (pdv[var]) {
        pvar = upd->vp[pg->imtrx][var];

        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += resid[i].dx(ad_fv->offset[EDDY_NU] + j);
        } /* End of loop over j */
      }   /* End of if the variable is active */

      /* Sensitivity w.r.t. velocity */
      for (b = 0; b < VIM; b++) {
        var = VELOCITY1 + b;
        if (pdv[var]) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += resid[i].dx(ad_fv->offset[var] + j);

          } /* End of loop over j */
        }   /* End of if the variale is active */
      }     /* End of loop over velocity components */

    } /* End of loop over i */
  }   /* End of if assemble Jacobian */
  return (status);
}

ADType turb_omega_wall_bc(void) {
  double rho = density(NULL, tran->time_value);
  double beta1 = 0.075;
  double nu = mp->viscosity / rho;
  double Dy = 0.0;

  if (upd->turbulent_info->use_internal_wall_distance) {
    fv->wall_distance = 0.;
    if (pd->gv[pd->ShapeVar]) {
      int dofs = ei[upd->matrix_index[pd->ShapeVar]]->dof[pd->ShapeVar];
      for (int i = 0; i < dofs; i++) {
        dbl d =
            upd->turbulent_info
                ->wall_distances[ei[upd->matrix_index[pd->ShapeVar]]->gnn_list[pd->ShapeVar][i]];
        if (d > Dy) {
          Dy = d;
        }
      }
    }
  } else {
    GOMA_EH(GOMA_ERROR, "Unimplemented wall distance for turb_omega_wall_bc\n");
  }
  double omega_wall = 10.0 * (6 * nu) / (beta1 * Dy * Dy);
  return omega_wall;
}

static void
calc_blending_functions(dbl rho, DENSITY_DEPENDENCE_STRUCT *d_rho, ADType &F1, ADType &F2) {
  dbl d = fv->wall_distance;
  ADType k = ad_fv->turb_k;
  ADType omega = ad_fv->turb_omega;
  dbl nu = mp->viscosity / rho;
  dbl beta_star = 0.09;
  dbl sigma_omega2 = 0.856;
  ADType ddkdomega = 0;
  for (int i = 0; i < pd->Num_Dim; i++) {
    ddkdomega += ad_fv->grad_turb_k[i] * ad_fv->grad_turb_omega[i];
  }
  ADType CD_komega = std::max(2 * rho * sigma_omega2 * ddkdomega / (omega + 1e-20), 1e-20);
  ADType C500 = 500 * nu / (d * d * omega + 1e-20);

  ADType arg11 = 0;
  if (k > 0) {
    ADType arg11 = sqrt(k) / (beta_star * omega * d + 1e-20);
  }
  ADType arg12 = 4 * rho * sigma_omega2 * k / (CD_komega * d * d + 1e-20);
  ADType arg1 = std::min(std::max(arg11, C500), arg12);
  // ADType arg1 = std::min(std::max(sqrt(std::max(k+1e-20, 0)) / (beta_star * omega * d + 1e-20),
  // C500),
  //                        4 * rho * sigma_omega2 * k / (CD_komega * d * d + 1e-20));
  ADType arg2 = std::max(2 * arg11, C500);

  F1 = tanh(arg1 * arg1 * arg1 * arg1);
  F2 = tanh(arg2 * arg2);
#if 0
  if (d_F1 != NULL) {
    dbl d_dot_k_omega_dk[MDE];
    dbl d_dot_k_omega_domega[MDE];
    for (int i = 0; i < pd->Num_Dim; i++) {
      for (int j = 0; j < ei[pg->imtrx]->dof[TURB_K]; j++) {
        d_dot_k_omega_dk[j] = bf[TURB_K]->grad_phi[j][i] * fv->grad_turb_omega[i];
      }
      for (int j = 0; j < ei[pg->imtrx]->dof[TURB_OMEGA]; j++) {
        d_dot_k_omega_domega[j] = bf[TURB_OMEGA]->grad_phi[j][i] * fv->grad_turb_k[i];
      }
    }

    dbl d_CD_komega_dk[MDE] = {0.};
    dbl d_CD_komega_domega[MDE] = {0.};
    if (CD_komega > 1e-10) {
      for (int j = 0; j < ei[pg->imtrx]->dof[TURB_K]; j++) {
        d_CD_komega_dk[j] = 2 * rho * sigma_omega2 * d_dot_k_omega_dk[j] / omega;
      }
      for (int j = 0; j < ei[pg->imtrx]->dof[TURB_K]; j++) {
        d_CD_komega_domega[j] =
            2 * rho * sigma_omega2 * d_dot_k_omega_domega[j] / omega -
            2 * rho * sigma_omega2 * dot_k_omega * bf[TURB_OMEGA]->phi[j] / (omega * omega);
      }
    }

    
  }
#endif
}

/* assemble_turb_k -- assemble terms (Residual & Jacobian) for conservation
 *
 * SST-2003m turbulence model
 *
 * in:
 *      time value
 *      theta (time stepping parameter, 0 for BE, 0.5 for CN)
 *      time step size
 *      Streamline Upwind Petrov Galerkin (PG) data structure
 *
 * out:
 *      lec -- gets loaded up with local contributions to resid, Jacobian
 *
 * Created:    July 2023 Weston Ortiz
 *
 */
extern "C" int ad_assemble_turb_k(dbl time_value, /* current time */
                                  dbl tt,         /* parameter to vary time integration from
                                                     explicit (tt = 1) to implicit (tt = 0)    */
                                  dbl dt,         /* current time step size                    */
                                  const PG_DATA *pg_data) {

  //! WIM is the length of the velocity vector
  int i;
  int eqn, peqn;
  int *pdv = pd->v[pg->imtrx];

  int status = 0;

  eqn = TURB_K;

  dbl d_area = fv->wt * bf[eqn]->detJ * fv->h3;

  dbl mu = mp->viscosity;
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;
  dbl rho = density(d_rho, time_value);
  ADType F1 = 0;
  ADType F2 = 0;
  calc_blending_functions(rho, d_rho, F1, F2);

  ADType SI;
  ADType gamma_dot[DIM][DIM];
  for (int i = 0; i < DIM; i++) {
    for (int j = 0; j < DIM; j++) {
      gamma_dot[i][j] = (ad_fv->grad_v[i][j] + ad_fv->grad_v[j][i]);
    }
  }
  ad_calc_shearrate(SI, gamma_dot);
  dbl a1 = 0.31;

  /* Rate of rotation tensor  */
  ADType omega[DIM][DIM];
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      omega[a][b] = (ad_fv->grad_v[a][b] - ad_fv->grad_v[b][a]);
    }
  }

  /* Vorticity */
  ADType Omega = 0.0;
  calc_vort_mag(Omega, omega);

  dbl x_dot[DIM] = {0.};
  if (pd->gv[R_MESH1]) {
    for (int i = 0; i < DIM; i++) {
      x_dot[i] = fv_dot->x[i];
    }
  }

  dbl supg = 1.;
  ADType supg_tau = 0;
  if (mp->SAwt_funcModel == GALERKIN) {
    supg = 0.;
  } else if (mp->SAwt_funcModel == SUPG || mp->SAwt_funcModel == SUPG_GP ||
             mp->SAwt_funcModel == SUPG_SHAKIB) {
    supg = mp->SAwt_func;
    ad_supg_tau_shakib(supg_tau, pd->Num_Dim, dt, mu, TURB_K);
  }

  dbl beta_star = 0.09;
  dbl sigma_k1 = 0.85;
  dbl sigma_k2 = 1.0;

  // blended values
  ADType sigma_k = F1 * sigma_k1 + (1 - F1) * sigma_k2;

  ADType mu_t = rho * a1 * fv_old->turb_k / (std::max(a1 * ad_fv->turb_omega, Omega * F2) + 1e-16);
  ADType P = mu_t * SI * SI;
  ADType Plim = std::min(P, 20 * beta_star * rho * ad_fv->turb_omega * fv_old->turb_k);

  std::vector<ADType> resid(ei[pg->imtrx]->dof[eqn]);
  for (int i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
    resid[i] = 0;
  }
  /*
   * Residuals_________________________________________________________________
   */
  if (af->Assemble_Residual) {
    /*
     * Assemble residual for eddy viscosity
     */
    eqn = TURB_K;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ADType wt_func = bf[eqn]->phi[i];

      if (supg > 0) {
        if (supg != 0.0) {
          for (int p = 0; p < VIM; p++) {
            wt_func += supg * supg_tau * ad_fv->v[p] * bf[eqn]->grad_phi[i][p];
          }
        }
      }

      /* Assemble mass term */
      ADType mass = 0.0;
      if (pd->TimeIntegration != STEADY) {
        if (pd->e[pg->imtrx][eqn] & T_MASS) {
          mass += rho * ad_fv->turb_k_dot * wt_func * d_area;
          mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
        }
      }

      /* Assemble advection term */
      ADType adv = 0;
      for (int p = 0; p < VIM; p++) {
        adv += rho * (ad_fv->v[p] - x_dot[p]) * ad_fv->grad_turb_k[p];
      }
      adv *= wt_func * d_area;
      adv *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

      /* Assemble source terms */
      ADType src = Plim - beta_star * rho * ad_fv->turb_omega * fv_old->turb_k;
      src *= -wt_func * d_area;
      src *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

      /* Assemble diffusion terms */
      ADType diff = 0.0;
      for (int p = 0; p < VIM; p++) {
        diff += bf[eqn]->grad_phi[i][p] * (mu + mu_t * sigma_k) * ad_fv->grad_turb_k[p];
      }
      diff *= d_area;
      diff *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];

      resid[i] += mass + adv + src + diff;
      lec->R[LEC_R_INDEX(peqn, i)] += mass.val() + adv.val() + src.val() + diff.val();
    } /* end of for (i=0,ei[pg->imtrx]->dofs...) */
  }   /* end of if assemble residual */

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    eqn = TURB_K;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      for (int var = V_FIRST; var < V_LAST; var++) {

        /* Sensitivity w.r.t. velocity */
        if (pdv[var]) {
          int pvar = upd->vp[pg->imtrx][var];

          for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            // J = &(lec->J[LEC_J_INDEX(peqn, pvar, ii, 0)]);
            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += resid[i].dx(ad_fv->offset[var] + j);

          } /* End of loop over j */
        }   /* End of if the variale is active */
      }
    }
  } /* End of if assemble Jacobian */
  return (status);
}

ADType ad_only_turb_k_omega_sst_viscosity(void) {
  ADType mu = 0;
  double mu_newt = mp->viscosity;
  dbl rho;
  rho = density(NULL, tran->time_value);
  ADType F1 = 0;
  ADType F2 = 0;
  calc_blending_functions(rho, NULL, F1, F2);

  ADType omega[DIM][DIM];
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      omega[a][b] = (ad_fv->grad_v[a][b] - ad_fv->grad_v[b][a]);
    }
  }
  /* Vorticity */
  ADType Omega = 0.0;
  calc_vort_mag(Omega, omega);
  dbl a1 = 0.31;
  ADType mu_t = rho * a1 * ad_fv->turb_k / (std::max(a1 * ad_fv->turb_omega, Omega * F2) + 1e-16);

  mu = mu_newt + mu_t;
  return mu;
}

extern "C" dbl ad_turb_k_omega_sst_viscosity(VISCOSITY_DEPENDENCE_STRUCT *d_mu) {
  ADType mu = 0;
  double mu_newt = mp->viscosity;
  dbl rho;
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;
  rho = density(d_rho, tran->time_value);
  ADType F1 = 0;
  ADType F2 = 0;
  calc_blending_functions(rho, d_rho, F1, F2);

  ADType omega[DIM][DIM];
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      omega[a][b] = (ad_fv->grad_v[a][b] - ad_fv->grad_v[b][a]);
    }
  }
  /* Vorticity */
  ADType Omega = 0.0;
  calc_vort_mag(Omega, omega);
  dbl a1 = 0.31;
  ADType mu_t = rho * a1 * ad_fv->turb_k / (std::max(a1 * ad_fv->turb_omega, Omega * F2) + 1e-16);

  mu = mu_newt + mu_t;
  if (d_mu != NULL) {
    for (int j = 0; j < ei[pg->imtrx]->dof[TURB_OMEGA]; j++) {
      d_mu->turb_omega[j] = mu.dx(ad_fv->offset[TURB_OMEGA] + j);
    }
    for (int j = 0; j < ei[pg->imtrx]->dof[TURB_K]; j++) {
      d_mu->turb_k[j] = mu.dx(ad_fv->offset[TURB_K] + j);
    }
    for (int b = 0; b < VIM; b++) {
      int var = VELOCITY1 + b;
      for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_mu->v[b][j] = mu.dx(ad_fv->offset[var] + j);
      }
    }
  }
  return mu.val();
}

/* assemble_turb_omega -- assemble terms (Residual & Jacobian) for conservation
 *
 * k-omega SST
 *
 * in:
 *      time value
 *      theta (time stepping parameter, 0 for BE, 0.5 for CN)
 *      time step size
 *      Streamline Upwind Petrov Galerkin (PG) data structure
 *
 * out:
 *      lec -- gets loaded up with local contributions to resid, Jacobian
 *
 * Created:    July 2023 Weston Ortiz
 *
 */
extern "C" int ad_assemble_turb_omega(dbl time_value, /* current time */
                                      dbl tt,         /* parameter to vary time integration from
                                                         explicit (tt = 1) to implicit (tt = 0)    */
                                      dbl dt, /* current time step size                    */
                                      const PG_DATA *pg_data) {

  //! WIM is the length of the velocity vector
  int i;
  int eqn, peqn;
  int *pdv = pd->v[pg->imtrx];

  int status = 0;

  eqn = TURB_OMEGA;

  dbl d_area = fv->wt * bf[eqn]->detJ * fv->h3;

  dbl mu = mp->viscosity;
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;
  dbl rho = density(d_rho, time_value);
  ADType F1 = 0;
  ADType F2 = 0;
  calc_blending_functions(rho, d_rho, F1, F2);

  ADType SI;
  ADType gamma_dot[DIM][DIM];
  for (int i = 0; i < DIM; i++) {
    for (int j = 0; j < DIM; j++) {
      gamma_dot[i][j] = (ad_fv->grad_v[i][j] + ad_fv->grad_v[j][i]);
    }
  }

  ADType omega[DIM][DIM];
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      omega[a][b] = (ad_fv->grad_v[a][b] - ad_fv->grad_v[b][a]);
    }
  }
  /* Vorticity */
  ADType Omega = 0.0;
  calc_vort_mag(Omega, omega);

  ad_calc_shearrate(SI, gamma_dot);
  dbl a1 = 0.31;

  dbl x_dot[DIM] = {0.};
  if (pd->gv[R_MESH1]) {
    for (int i = 0; i < DIM; i++) {
      x_dot[i] = fv_dot->x[i];
    }
  }

  dbl supg = 1.;
  ADType supg_tau;
  if (mp->SAwt_funcModel == GALERKIN) {
    supg = 0.;
  } else if (mp->SAwt_funcModel == SUPG || mp->SAwt_funcModel == SUPG_GP ||
             mp->SAwt_funcModel == SUPG_SHAKIB) {
    supg = mp->SAwt_func;
    ad_supg_tau_shakib(supg_tau, pd->Num_Dim, dt, mu, TURB_OMEGA);
  }

  dbl beta_star = 0.09;
  dbl sigma_k1 = 0.85;
  dbl sigma_k2 = 1.0;
  dbl sigma_omega1 = 0.5;
  dbl sigma_omega2 = 0.856;
  dbl beta1 = 0.075;
  dbl beta2 = 0.0828;
  dbl gamma1 = beta1 / beta_star;
  dbl gamma2 = beta2 / beta_star;

  // blended values
  ADType gamma = F1 * gamma1 + (1 - F1) * gamma2;
  ADType sigma_k = F1 * sigma_k1 + (1 - F1) * sigma_k2;
  ADType sigma_omega = F1 * sigma_omega1 + (1 - F1) * sigma_omega2;
  ADType beta = F1 * beta1 + (1 - F1) * beta2;

  ADType mu_t = rho * a1 * ad_fv->turb_k / (std::max(a1 * fv_old->turb_omega, Omega * F2) + 1e-16);
  ADType P = mu_t * SI * SI;
  ADType Plim = std::min(P, 10 * beta_star * rho * fv_old->turb_omega * ad_fv->turb_k);

  std::vector<ADType> resid(ei[pg->imtrx]->dof[eqn]);
  for (int i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
    resid[i] = 0;
  }
  /*
   * Residuals_________________________________________________________________
   */
  if (af->Assemble_Residual) {
    /*
     * Assemble residual for eddy viscosity
     */
    eqn = TURB_OMEGA;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ADType wt_func = bf[eqn]->phi[i];

      if (supg > 0) {
        if (supg != 0.0) {
          for (int p = 0; p < VIM; p++) {
            wt_func += supg * supg_tau * ad_fv->v[p] * bf[eqn]->grad_phi[i][p];
          }
        }
      }

      /* Assemble mass term */
      ADType mass = 0.0;
      if (pd->TimeIntegration != STEADY) {
        if (pd->e[pg->imtrx][eqn] & T_MASS) {
          mass += rho * ad_fv->turb_omega_dot * wt_func * d_area;
          mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
        }
      }

      /* Assemble advection term */
      ADType adv = 0;
      for (int p = 0; p < VIM; p++) {
        adv += rho * (ad_fv->v[p] - x_dot[p]) * ad_fv->grad_turb_omega[p];
      }
      adv *= wt_func * d_area;
      adv *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

      /* Assemble source terms */
      ADType src1 = (gamma * rho / (mu_t + 1e-16)) * Plim -
                    beta * rho * fv_old->turb_omega * fv_old->turb_omega;
      ADType src2 = 0;
      for (int p = 0; p < pd->Num_Dim; p++) {
        src2 += ad_fv->grad_turb_k[p] * fv_old->grad_turb_omega[p];
      }
      src2 *= 2 * (1 - F1) * rho * sigma_omega2 / (fv_old->turb_omega + 1e-16);
      ADType src = src1 + src2;
      src *= -wt_func * d_area;
      src *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

      /* Assemble diffusion terms */
      ADType diff = 0.0;
      for (int p = 0; p < VIM; p++) {
        diff += bf[eqn]->grad_phi[i][p] * (mu + mu_t * sigma_omega) * ad_fv->grad_turb_omega[p];
      }
      diff *= d_area;
      diff *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];

      resid[i] += mass + adv + src + diff;
      lec->R[LEC_R_INDEX(peqn, i)] += mass.val() + adv.val() + src.val() + diff.val();
    } /* end of for (i=0,ei[pg->imtrx]->dofs...) */
  }   /* end of if assemble residual */

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    eqn = TURB_OMEGA;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      for (int var = V_FIRST; var < V_LAST; var++) {

        /* Sensitivity w.r.t. velocity */
        if (pdv[var]) {
          int pvar = upd->vp[pg->imtrx][var];

          for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            // J = &(lec->J[LEC_J_INDEX(peqn, pvar, ii, 0)]);
            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += resid[i].dx(ad_fv->offset[var] + j);

          } /* End of loop over j */
        }   /* End of if the variale is active */
      }
    }
  } /* End of if assemble Jacobian */
  return (status);
}

ADType clipping_func(ADType x, ADType x_s, ADType x_e) {
  ADType psi = 0;
  ADType theta = M_PI_2 * (2 * x - (x_s + x_e)) / (x_s - x_e);

  if (x < x_s) {
    psi = 0;
  } else if (x > x_e) {
    psi = 1;
  } else {
    psi = 0.5 * (1 + std::sin(theta));
  }

  return psi;
}

ADType ad_only_turb_k_omega_viscosity(void) {
  ADType mu = 0;
  double mu_newt = mp->viscosity;
  dbl k_inf = upd->turbulent_info->k_inf;
  ADType psi = clipping_func(ad_fv->turb_k, 0, 10 * k_inf);
  ADType psi_neg = clipping_func(-ad_fv->turb_k, 0, 10 * k_inf);
  dbl rho = density(NULL, tran->time_value);

  ADType mu_t = rho * psi * ad_fv->turb_k / (std::exp(ad_fv->turb_omega));
  mu = mu_newt + mu_t;
  return mu;
}

extern "C" int ad_assemble_turb_k_modified(dbl time_value, /* current time */
                                           dbl tt, /* parameter to vary time integration from
                                                      explicit (tt = 1) to implicit (tt = 0)    */
                                           dbl dt, /* current time step size                    */
                                           const PG_DATA *pg_data) {

  //! WIM is the length of the velocity vector
  int i;
  int eqn, peqn;
  int *pdv = pd->v[pg->imtrx];

  int status = 0;

  eqn = TURB_K;

  dbl d_area = fv->wt * bf[eqn]->detJ * fv->h3;

  dbl omega_inf = std::exp(upd->turbulent_info->omega_inf);
  dbl k_inf = upd->turbulent_info->k_inf;

  // logarithmic formulation of k-omega model
  ADType omega = std::exp(ad_fv->turb_omega);
  ADType psi = clipping_func(ad_fv->turb_k, 0, 10 * k_inf);
  ADType psi_neg = clipping_func(-fv_old->turb_k, 0, 10 * k_inf);
  ADType k = psi * ad_fv->turb_k;
  // ADType k = psi * fv_old->turb_k;

  dbl sigma_k = 0.5;
  dbl beta_star = 9.0 / 100.0;

  ADType SI;
  ADType gamma_dot[DIM][DIM];
  for (int i = 0; i < DIM; i++) {
    for (int j = 0; j < DIM; j++) {
      gamma_dot[i][j] = (ad_fv->grad_v[i][j] + ad_fv->grad_v[j][i]);
    }
  }

  ADType Omega_tens[DIM][DIM];
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      Omega_tens[a][b] = (ad_fv->grad_v[a][b] - ad_fv->grad_v[b][a]);
    }
  }
  /* Vorticity */
  ADType Omega = 0.0;
  calc_vort_mag(Omega, Omega_tens);

  ad_calc_shearrate(SI, gamma_dot);

  dbl rho = density(NULL, time_value);

  ADType mu_t = rho * k / (omega + 1e-16);
  // ADType mu_t_kdiff = rho * k / (omega + 1e-16) - rho * psi_neg * ad_fv->turb_k / omega_inf;
  ADType mu_t_kdiff = rho * k / (omega + 1e-16) - rho * psi_neg * fv_old->turb_k / omega_inf;
  dbl mu = mp->viscosity;

  ADType P = mu_t * Omega * Omega;

  P = std::min(P, 10 * beta_star * rho * omega * std::max(k, 0));

  dbl supg = 1.;
  ADType supg_tau;
  if (mp->SAwt_funcModel == GALERKIN) {
    supg = 0.;
  } else if (mp->SAwt_funcModel == SUPG || mp->SAwt_funcModel == SUPG_GP ||
             mp->SAwt_funcModel == SUPG_SHAKIB) {
    supg = mp->SAwt_func;
    ad_supg_tau_shakib(supg_tau, pd->Num_Dim, dt, mu_t, TURB_OMEGA);
  }

  std::vector<ADType> resid(ei[pg->imtrx]->dof[TURB_K]);
  for (int i = 0; i < ei[pg->imtrx]->dof[TURB_K]; i++) {
    resid[i] = 0;
  }
  /*
   * Residuals_________________________________________________________________
   */
  if (af->Assemble_Residual) {
    eqn = TURB_K;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ADType wt_func = bf[eqn]->phi[i];

      if (supg > 0) {
        if (supg != 0.0) {
          for (int p = 0; p < VIM; p++) {
            wt_func += supg * supg_tau * ad_fv->v[p] * bf[eqn]->grad_phi[i][p];
          }
        }
      }

      /* Assemble mass term */
      ADType mass = 0.0;
      if (pd->TimeIntegration != STEADY) {
        if (pd->e[pg->imtrx][eqn] & T_MASS) {
          mass += rho * ad_fv->turb_k_dot * wt_func * d_area;
          mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
        }
      }

      /* Assemble advection term */
      ADType adv = 0;
      for (int p = 0; p < VIM; p++) {
        adv += rho * (ad_fv->v[p] - ad_fv->x_dot[p]) * ad_fv->grad_turb_k[p];
      }
      adv *= wt_func * d_area;
      adv *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

      /* Assemble source terms */
      ADType src1 = P - beta_star * rho * omega * std::max(k, 0);
      ADType src2 = 0;

      // ADType src2 = -beta_star * rho * omega_inf * psi_neg * ad_fv->turb_k;
      // ADType src2 = -beta_star * rho * omega_inf * psi_neg * ad_fv->turb_k;
      ADType src = (src1 + src2);
      // src = src1;
      src *= -wt_func * d_area;
      src *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

      /* Assemble diffusion terms */
      ADType diff = 0.0;
      for (int p = 0; p < VIM; p++) {
        diff += bf[eqn]->grad_phi[i][p] * (mu + mu_t_kdiff * sigma_k) * ad_fv->grad_turb_k[p];
      }
      diff *= d_area;
      diff *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];

      resid[i] -= mass + adv + src + diff;
      lec->R[LEC_R_INDEX(peqn, i)] -= mass.val() + adv.val() + src.val() + diff.val();
    } /* end of for (i=0,ei[pg->imtrx]->dofs...) */
  }   /* end of if assemble residual */

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    eqn = TURB_K;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      for (int var = V_FIRST; var < V_LAST; var++) {

        /* Sensitivity w.r.t. velocity */
        if (pdv[var]) {
          int pvar = upd->vp[pg->imtrx][var];

          for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            // J = &(lec->J[LEC_J_INDEX(peqn, pvar, ii, 0)]);
            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += resid[i].dx(ad_fv->offset[var] + j);

          } /* End of loop over j */
        }   /* End of if the variale is active */
      }
    }
  } /* End of if assemble Jacobian */
  return (status);
}
extern "C" int ad_assemble_turb_omega_modified(dbl time_value, /* current time */
                                               dbl tt, /* parameter to vary time integration from
                                                          explicit (tt = 1) to implicit (tt = 0) */
                                               dbl dt, /* current time step size */
                                               const PG_DATA *pg_data) {

  //! WIM is the length of the velocity vector
  int i;
  int eqn, peqn;
  int *pdv = pd->v[pg->imtrx];

  int status = 0;

  eqn = TURB_OMEGA;

  dbl d_area = fv->wt * bf[eqn]->detJ * fv->h3;

  dbl omega_inf = std::exp(upd->turbulent_info->omega_inf);
  dbl k_inf = upd->turbulent_info->k_inf;

  // logarithmic formulation of k-omega model
  ADType omega = std::exp(ad_fv->turb_omega);
  ADType psi = clipping_func(ad_fv->turb_k, 0, 10 * k_inf);
  ADType psi_neg = clipping_func(-ad_fv->turb_k, 0, 10 * k_inf);
  ADType k = psi * ad_fv->turb_k;

  dbl sigma_omega = 0.5;
  dbl beta_star = 9.0 / 100.0;
  dbl beta = 3.0 / 40.;
  dbl gamma = 5.0 / 9.0;

  ADType SI;
  ADType gamma_dot[DIM][DIM];
  for (int i = 0; i < DIM; i++) {
    for (int j = 0; j < DIM; j++) {
      gamma_dot[i][j] = (ad_fv->grad_v[i][j] + ad_fv->grad_v[j][i]);
    }
  }

  ADType Omega_tens[DIM][DIM];
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      Omega_tens[a][b] = (ad_fv->grad_v[a][b] - ad_fv->grad_v[b][a]);
    }
  }
  /* Vorticity */
  ADType Omega = 0.0;
  calc_vort_mag(Omega, Omega_tens);

  ad_calc_shearrate(SI, gamma_dot);

  dbl rho = density(NULL, time_value);

  ADType mu_t = rho * k / (omega + 1e-16);
  ADType mu_t_kdiff = rho * k / (omega + 1e-16) - rho * psi_neg * ad_fv->turb_k / omega_inf;
  dbl mu = mp->viscosity;

  // ADType P = mu_t * SI * SI;
  ADType P = mu_t * Omega * Omega;
  P = std::min(P, 10 * beta_star * rho * omega * k);

  dbl supg = 1.;
  ADType supg_tau;
  if (mp->SAwt_funcModel == GALERKIN) {
    supg = 0.;
  } else if (mp->SAwt_funcModel == SUPG || mp->SAwt_funcModel == SUPG_GP ||
             mp->SAwt_funcModel == SUPG_SHAKIB) {
    supg = mp->SAwt_func;
    ad_supg_tau_shakib(supg_tau, pd->Num_Dim, dt, mu, TURB_OMEGA);
  }

  std::vector<ADType> resid(ei[pg->imtrx]->dof[TURB_OMEGA]);
  for (int i = 0; i < ei[pg->imtrx]->dof[TURB_OMEGA]; i++) {
    resid[i] = 0;
  }
  /*
   * Residuals_________________________________________________________________
   */
  if (af->Assemble_Residual) {
    /*
     * Assemble residual for eddy viscosity
     */
    eqn = TURB_OMEGA;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ADType wt_func = bf[eqn]->phi[i];

      if (supg > 0) {
        if (supg != 0.0) {
          for (int p = 0; p < VIM; p++) {
            wt_func += supg * supg_tau * ad_fv->v[p] * bf[eqn]->grad_phi[i][p];
          }
        }
      }

      /* Assemble mass term */
      ADType mass = 0.0;
      if (pd->TimeIntegration != STEADY) {
        if (pd->e[pg->imtrx][eqn] & T_MASS) {
          mass += rho * ad_fv->turb_omega_dot * wt_func * d_area;
          mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
        }
      }

      /* Assemble advection term */
      ADType adv = 0;
      for (int p = 0; p < VIM; p++) {
        adv += rho * (ad_fv->v[p] - ad_fv->x_dot[p]) * ad_fv->grad_turb_omega[p];
      }
      adv *= wt_func * d_area;
      adv *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

      /* Assemble source terms */
      ADType src1 = (gamma / (k + 1e-16)) * P - beta * rho * omega;

      ADType src2 = 0;
      for (int p = 0; p < pd->Num_Dim; p++) {
        src2 += (mu + sigma_omega * mu_t) * ad_fv->grad_turb_omega[p] * ad_fv->grad_turb_omega[p];
      }
      ADType src = src1 + src2;
      src *= -wt_func * d_area;
      src *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

      /* Assemble diffusion terms */
      ADType diff = 0.0;
      for (int p = 0; p < VIM; p++) {
        diff += bf[eqn]->grad_phi[i][p] * (mu + mu_t * sigma_omega) * ad_fv->grad_turb_omega[p];
      }
      diff *= d_area;
      diff *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];

      resid[i] -= mass + adv + src + diff;
      lec->R[LEC_R_INDEX(peqn, i)] -= mass.val() + adv.val() + src.val() + diff.val();
    } /* end of for (i=0,ei[pg->imtrx]->dofs...) */

  } /* end of if assemble residual */

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    eqn = TURB_OMEGA;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      for (int var = V_FIRST; var < V_LAST; var++) {

        /* Sensitivity w.r.t. velocity */
        if (pdv[var]) {
          int pvar = upd->vp[pg->imtrx][var];

          for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            // J = &(lec->J[LEC_J_INDEX(peqn, pvar, ii, 0)]);
            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += resid[i].dx(ad_fv->offset[var] + j);

          } /* End of loop over j */
        }   /* End of if the variale is active */
      }
    }
  } /* End of if assemble Jacobian */
  return (status);
}
extern "C" int
ad_assemble_turb_k_omega_modified(dbl time_value, /* current time */
                                  dbl tt,         /* parameter to vary time integration from
                                                     explicit (tt = 1) to implicit (tt = 0)    */
                                  dbl dt,         /* current time step size                    */
                                  const PG_DATA *pg_data) {

  //! WIM is the length of the velocity vector
  int i;
  int eqn, peqn;
  int *pdv = pd->v[pg->imtrx];

  int status = 0;

  eqn = TURB_OMEGA;

  dbl d_area = fv->wt * bf[eqn]->detJ * fv->h3;

  dbl omega_inf = upd->turbulent_info->omega_inf;
  dbl k_inf = upd->turbulent_info->k_inf;

  // logarithmic formulation of k-omega model
  ADType omega = std::exp(ad_fv->turb_omega);
  ADType psi = clipping_func(ad_fv->turb_k, 0, 10 * k_inf);
  ADType psi_neg = clipping_func(-ad_fv->turb_k, 0, 10 * k_inf);
  ADType k = psi * ad_fv->turb_k;

  dbl sigma_k = 0.5;
  dbl sigma_omega = 0.5;
  dbl beta_star = 9.0 / 100.0;
  dbl beta = 3.0 / 40.;
  dbl gamma = 5.0 / 9.0;

  ADType SI;
  ADType gamma_dot[DIM][DIM];
  for (int i = 0; i < DIM; i++) {
    for (int j = 0; j < DIM; j++) {
      gamma_dot[i][j] = (ad_fv->grad_v[i][j] + ad_fv->grad_v[j][i]);
    }
  }

  ADType Omega_tens[DIM][DIM];
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      Omega_tens[a][b] = (ad_fv->grad_v[a][b] - ad_fv->grad_v[b][a]);
    }
  }
  /* Vorticity */
  ADType Omega = 0.0;
  calc_vort_mag(Omega, Omega_tens);

  ad_calc_shearrate(SI, gamma_dot);

  dbl rho = density(NULL, time_value);

  ADType mu_t = rho * k / (omega);
  // ADType mu_t_kdiff = rho * k / (omega);// - rho * psi_neg * ad_fv->turb_k / omega_inf;
  ADType mu_t_kdiff = rho * k / (omega)-rho * psi_neg * ad_fv->turb_k / omega_inf;
  dbl mu = mp->viscosity;

  ADType P = mu_t * Omega * Omega;

  dbl supg = 1.;
  // if (mp->SAwt_funcModel == GALERKIN) {
  // supg = 0.;
  // } else if (mp->SAwt_funcModel == SUPG || mp->SAwt_funcModel == SUPG_GP ||
  //  mp->SAwt_funcModel == SUPG_SHAKIB) {
  // supg = mp->SAwt_func;
  // }

  std::vector<std::vector<ADType>> resid(2);
  resid[0].resize(ei[pg->imtrx]->dof[TURB_OMEGA]);
  resid[1].resize(ei[pg->imtrx]->dof[TURB_K]);
  for (int i = 0; i < ei[pg->imtrx]->dof[TURB_OMEGA]; i++) {
    resid[0][i] = 0;
  }
  for (int i = 0; i < ei[pg->imtrx]->dof[TURB_K]; i++) {
    resid[1][i] = 0;
  }
  /*
   * Residuals_________________________________________________________________
   */
  if (af->Assemble_Residual) {
    ADType gs_inner_dot[DIM];
    ADType supg_tau_w, supg_tau_k;
    ad_supg_tau_shakib(supg_tau_w, pd->Num_Dim, dt, mu + sigma_omega * mu_t, TURB_OMEGA);
    ad_supg_tau_shakib(supg_tau_k, pd->Num_Dim, dt, mu + sigma_k * mu_t, TURB_K);
    ADType vshock_k = ad_dcdd(pd->Num_Dim, TURB_K, ad_fv->grad_turb_k, gs_inner_dot);
    vshock_k = std::min(supg_tau_k, vshock_k);
    ADType vshock_w = ad_dcdd(pd->Num_Dim, TURB_OMEGA, ad_fv->grad_turb_omega, gs_inner_dot);
    vshock_k = std::min(supg_tau_w, vshock_w);
    /*
     * Assemble residual for eddy viscosity
     */
    eqn = TURB_OMEGA;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ADType wt_func = bf[eqn]->phi[i];

      if (supg > 0) {
        if (supg != 0.0) {
          for (int p = 0; p < VIM; p++) {
            wt_func += supg * supg_tau_w * ad_fv->v[p] * bf[eqn]->grad_phi[i][p];
          }
        }
      }

      /* Assemble mass term */
      ADType mass = 0.0;
      if (pd->TimeIntegration != STEADY) {
        if (pd->e[pg->imtrx][eqn] & T_MASS) {
          mass += rho * ad_fv->turb_omega_dot * wt_func * d_area;
          mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
        }
      }

      /* Assemble advection term */
      ADType adv = 0;
      for (int p = 0; p < VIM; p++) {
        adv += rho * (ad_fv->v[p] - ad_fv->x_dot[p]) * ad_fv->grad_turb_omega[p];
      }
      adv *= wt_func * d_area;
      adv *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

      /* Assemble source terms */
      ADType src1 = (gamma / (k + 1e-16)) * std::min(P, 10 * beta_star * rho * omega * k) -
                    beta * rho * omega;

      ADType src2 = 0;
      for (int p = 0; p < pd->Num_Dim; p++) {
        src2 += (mu + sigma_omega * mu_t) * ad_fv->grad_turb_omega[p] * ad_fv->grad_turb_omega[p];
      }
      ADType src = src1 + src2;
      src *= -wt_func * d_area;
      src *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

      /* Assemble diffusion terms */
      ADType diff = 0.0;
      for (int p = 0; p < VIM; p++) {
        diff += bf[eqn]->grad_phi[i][p] * (mu + mu_t * sigma_omega + vshock_w) *
                ad_fv->grad_turb_omega[p];
      }
      diff *= d_area;
      diff *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];

      resid[0][i] -= mass + adv + src + diff;
      lec->R[LEC_R_INDEX(peqn, i)] -= mass.val() + adv.val() + src.val() + diff.val();
    } /* end of for (i=0,ei[pg->imtrx]->dofs...) */

    eqn = TURB_K;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ADType wt_func = bf[eqn]->phi[i];

      if (supg > 0) {
        if (supg != 0.0) {
          for (int p = 0; p < VIM; p++) {
            wt_func += supg * supg_tau_k * ad_fv->v[p] * bf[eqn]->grad_phi[i][p];
          }
        }
      }

      /* Assemble mass term */
      ADType mass = 0.0;
      if (pd->TimeIntegration != STEADY) {
        if (pd->e[pg->imtrx][eqn] & T_MASS) {
          mass += rho * ad_fv->turb_k_dot * wt_func * d_area;
          mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
        }
      }

      /* Assemble advection term */
      ADType adv = 0;
      for (int p = 0; p < VIM; p++) {
        adv += rho * (ad_fv->v[p] - ad_fv->x_dot[p]) * ad_fv->grad_turb_k[p];
      }
      adv *= wt_func * d_area;
      adv *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

      /* Assemble source terms */
      ADType src1 = std::min(P, 10 * beta_star * rho * omega * k) - beta_star * rho * omega * k;
      // ADType src2 = 0;

      ADType src2 = -beta_star * rho * omega_inf * psi_neg * ad_fv->turb_k;
      ADType src = src1 + src2;
      src *= -wt_func * d_area;
      src *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

      /* Assemble diffusion terms */
      ADType diff = 0.0;
      for (int p = 0; p < VIM; p++) {
        diff += bf[eqn]->grad_phi[i][p] * (mu + mu_t_kdiff * sigma_k + vshock_k) *
                ad_fv->grad_turb_k[p];
      }
      diff *= d_area;
      diff *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];

      resid[1][i] -= mass + adv + src + diff;
      lec->R[LEC_R_INDEX(peqn, i)] -= mass.val() + adv.val() + src.val() + diff.val();
    } /* end of for (i=0,ei[pg->imtrx]->dofs...) */
  }   /* end of if assemble residual */

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    eqn = TURB_OMEGA;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      for (int var = V_FIRST; var < V_LAST; var++) {

        /* Sensitivity w.r.t. velocity */
        if (pdv[var]) {
          int pvar = upd->vp[pg->imtrx][var];

          for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            // J = &(lec->J[LEC_J_INDEX(peqn, pvar, ii, 0)]);
            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += resid[0][i].dx(ad_fv->offset[var] + j);

          } /* End of loop over j */
        }   /* End of if the variale is active */
      }
    }
    eqn = TURB_K;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      for (int var = V_FIRST; var < V_LAST; var++) {

        /* Sensitivity w.r.t. velocity */
        if (pdv[var]) {
          int pvar = upd->vp[pg->imtrx][var];

          for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            // J = &(lec->J[LEC_J_INDEX(peqn, pvar, ii, 0)]);
            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += resid[1][i].dx(ad_fv->offset[var] + j);

          } /* End of loop over j */
        }   /* End of if the variale is active */
      }
    }
  } /* End of if assemble Jacobian */
  return (status);
}

void compute_sst_blending(ADType &F1, ADType &F2) {
  dbl mu = mp->viscosity;
  dbl rho = density(NULL, tran->time_value);

  ADType grad_k_dot_grad_omega = 0.;
  for (int i = 0; i < DIM; i++) {
    grad_k_dot_grad_omega += ad_fv->grad_turb_k[i] * ad_fv->grad_turb_omega[i];
  }
  ADType omega = std::max(1e-20, ad_fv->turb_omega);
  ADType k = std::max(1e-20, ad_fv->turb_k);
  dbl sigma_omega2 = 0.856;
  dbl beta_star = 0.09;
  dbl d = fv->wall_distance;

  ADType CD_kw = std::max(2 * rho * sigma_omega2 * (1 / omega) * grad_k_dot_grad_omega, 1e-10);

  ADType arg1 =
      std::min(std::max(std::sqrt(k) / (beta_star * omega * d), 500 * (mu / rho) / (d * d * omega)),
               4 * rho * sigma_omega2 * k / (CD_kw * d * d));

  ADType arg2 =
      std::max(2 * std::sqrt(k) / (beta_star * omega * d), 500 * (mu / rho) / (d * d * omega));

  F1 = std::tanh(arg1 * arg1 * arg1 * arg1);

  F2 = std::tanh(arg2 * arg2);
}

ADType sst_viscosity(const ADType &Omega, const ADType &F2) {
  ADType omega = std::max(1e-20, ad_fv->turb_omega);
  ADType k = std::max(1e-20, ad_fv->turb_k);
  dbl rho = density(NULL, tran->time_value);
  dbl a1 = 0.31;
  ADType mu_turb = mp->viscosity + rho * a1 * k / (std::max(a1 * omega, Omega * F2));
  return mu_turb;
}

ADType ad_yzbeta(int dim, int eqn, ADType Y, ADType Z, const ADType grad_U[DIM]) {

  ADType invY = 1.0 / Y;
  ADType yzbeta = 0.0;

  ADType js = 0.0;
  for (int i = 0; i < dim; i++) {
    js += grad_U[i] * grad_U[i];
  }
  js = sqrt(std::max(js, 1e-20));

  ADType hdc = 0.0;
  for (int i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
    for (int a = 0; a < dim; a++) {
      hdc += std::abs(grad_U[a] / js) * ad_fv->basis[eqn].grad_phi[i][a];
    }
  }
  hdc = 2.0 / std::max(hdc, 1e-20);

  ADType inner = 0.0;
  for (int i = 0; i < dim; i++) {
    inner += invY * invY * grad_U[i] * grad_U[i];
  }
  inner = 1.0 / sqrt(std::max(inner, 1e-20));

  ADType yzbeta1 = std::abs(invY * Z) * inner;
  ADType yzbeta2 = std::abs(invY * Z) * std::pow(hdc / 2.0, 2.0);
  return 0.5 * (yzbeta1 + yzbeta2);
}

#if 1
extern "C" int ad_assemble_k_omega_sst_modified(dbl time_value, /* current time */
                                                dbl tt, /* parameter to vary time integration from
                                                           explicit (tt = 1) to implicit (tt = 0) */
                                                dbl dt, /* current time step size */
                                                const PG_DATA *pg_data) {

  //! WIM is the length of the velocity vector
  int i;
  int eqn, peqn;
  int *pdv = pd->v[pg->imtrx];

  int status = 0;

  eqn = TURB_K;

  dbl supg = 1.0;

  ADType G[DIM][DIM];

  ad_get_metric_tensor(ad_fv->B, pd->Num_Dim, ei[pg->imtrx]->ielem_type, G);

  ADType v_d_gv = 0;
  for (int i = 0; i < pd->Num_Dim; i++) {
    for (int j = 0; j < pd->Num_Dim; j++) {
      v_d_gv += fabs((ad_fv->v[i] - ad_fv->x_dot[i]) * G[i][j] * (ad_fv->v[j] - ad_fv->x_dot[j]));
    }
  }

  dbl rho = density(NULL, time_value);

  ADType W[DIM][DIM];
  for (int i = 0; i < DIM; i++) {
    for (int j = 0; j < DIM; j++) {
      W[i][j] = 0.5 * (ad_fv->grad_v[i][j] - ad_fv->grad_v[j][i]);
    }
  }

  ADType Omega = 0.0;
  for (int i = 0; i < DIM; i++) {
    for (int j = 0; j < DIM; j++) {
      Omega += W[i][j] * W[i][j];
    }
  }
  Omega = sqrt(std::max(Omega, 1e-20));

  dbl mu = mp->viscosity;

  ADType grad_k_dot_grad_omega = 0.;
  for (int i = 0; i < pd->Num_Dim; i++) {
    grad_k_dot_grad_omega += ad_fv->grad_turb_k[i] * ad_fv->grad_turb_omega[i];
  }

  ADType omega = std::max(1e-20, ad_fv->turb_omega);
  ADType k = std::max(1e-20, ad_fv->turb_k);
  dbl sigma_omega1 = 0.5;
  dbl sigma_omega2 = 0.856;
  dbl sigma_k1 = 0.85;
  dbl sigma_k2 = 1.0;
  dbl beta1 = 0.075;
  dbl beta2 = 0.0828;
  dbl beta_star = 0.09;
  // dbl kappa = 0.41;

  dbl gamma1 = 5.0 / 9.0; // beta1 / beta_star - sigma_omega1 * kappa * kappa / sqrt(beta_star);
  dbl gamma2 = 0.44;      // beta2 / beta_star - sigma_omega2 * kappa * kappa / sqrt(beta_star);

  ADType F1, F2;
  compute_sst_blending(F1, F2);

  ADType mu_turb = sst_viscosity(Omega, F2);

  ADType sigma_k = sigma_k1 * F1 + sigma_k2 * (1 - F1);
  ADType sigma_omega = sigma_omega1 * F1 + sigma_omega2 * (1 - F1);
  ADType beta = beta1 * F1 + beta2 * (1 - F1);

  ADType gamma = gamma1 * F1 + gamma2 * (1 - F1);

  ADType P = mu_turb * Omega * Omega;
  ADType Plim = std::min(P, 10 * beta_star * rho * omega * k);

  ADType d_area = fv->wt * ad_fv->detJ * fv->h3;
  std::vector<std::vector<ADType>> resid(2);
  resid[0].resize(ei[pg->imtrx]->dof[TURB_K]);
  resid[1].resize(ei[pg->imtrx]->dof[TURB_OMEGA]);
  for (int i = 0; i < ei[pg->imtrx]->dof[TURB_K]; i++) {
    resid[0][i] = 0;
  }
  for (int i = 0; i < ei[pg->imtrx]->dof[TURB_OMEGA]; i++) {
    resid[1][i] = 0;
  }
  /*
   * Residuals_________________________________________________________________
   */
  if (af->Assemble_Residual) {
    /*
     * Assemble residual for eddy viscosity
     */
    eqn = TURB_K;
    peqn = upd->ep[pg->imtrx][eqn];

    ADType Z_k = rho * ad_fv->turb_k_dot;
    for (int i = 0; i < pd->Num_Dim; i++) {
      Z_k += rho * ad_fv->v[i] * ad_fv->grad_turb_k[i];
    }
    ADType coeff = 12 * (mu + mu_turb * sigma_k) * (mu + mu_turb * sigma_k);

    ADType diff_g_g = 0;
    for (int i = 0; i < pd->Num_Dim; i++) {
      for (int j = 0; j < pd->Num_Dim; j++) {
        diff_g_g += coeff * G[i][j] * G[i][j];
      }
    }

    ADType sugn1 = 0;
    for (int i = 0; i < pd->Num_Dim; i++) {
      for (int j = 0; j < ei[pg->imtrx]->dof[eqn]; j++) {
        sugn1 += std::abs(ad_fv->v[i] * bf[eqn]->grad_phi[j][i]);
      }
    }
    sugn1 = 1.0 / std::max(sugn1, 1e-20);

    ADType sugn2 = dt / 2.0;

    ADType r[DIM] = {0};
    ADType norm_grad_k = 0;
    for (int i = 0; i < pd->Num_Dim; i++) {
      norm_grad_k = ad_fv->grad_turb_k[i] * ad_fv->grad_turb_k[i];
    }
    norm_grad_k = sqrt(norm_grad_k);
    for (int i = 0; i < pd->Num_Dim; i++) {
      r[i] = std::abs(ad_fv->grad_turb_k[i]) / norm_grad_k;
    }

    ADType h_rgn = 0;
    for (int i = 0; i < pd->Num_Dim; i++) {
      for (int j = 0; j < ei[pg->imtrx]->dof[eqn]; j++) {
        h_rgn += std::abs(r[i] * bf[eqn]->grad_phi[j][i]);
      }
    }
    h_rgn = 2.0 / std::max(h_rgn, 1e-20);

    ADType sugn3 = h_rgn * h_rgn * rho / (4 * (mu + mu_turb * sigma_k));
    sugn3 = 0;

    ADType tau_supg_k = 1.0 / (sqrt(4 / (dt * dt) + v_d_gv + diff_g_g));
    // ADType tau_supg_k = 1.0 / sqrt(1 / (sugn1 * sugn1) + 1 / (sugn2 * sugn2) + 1 / (sugn3 *
    // sugn3));

    // ADType vshock_k =
    // ad_yzbeta(pd->Num_Dim, TURB_K, upd->turbulent_info->k_inf, Z_k, ad_fv->grad_turb_k);

    ADType gs_inner_dot[DIM];
    ADType grad_turb_k[DIM] = {fv_old->grad_turb_k[0], fv_old->grad_turb_k[1],
                               fv_old->grad_turb_k[2]};
    ADType vshock_k = 0 * ad_dcdd(pd->Num_Dim, TURB_K, grad_turb_k, gs_inner_dot);

    vshock_k = std::min(tau_supg_k, vshock_k);
    // if (fv->wall_distance < 1) {
    //   // larger diffusion near wall
    //   vshock_k = std::max(vshock_k, 1e-2*std::exp(-fv->wall_distance * 20));
    // }

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ADType wt_func = bf[eqn]->phi[i];

      if (supg > 0) {
        if (supg != 0.0) {
          for (int p = 0; p < VIM; p++) {
            wt_func += supg * tau_supg_k * ad_fv->v[p] * bf[eqn]->grad_phi[i][p];
          }
        }
      }

      /* Assemble mass term */
      ADType mass = 0.0;
      if (pd->TimeIntegration != STEADY) {
        if (pd->e[pg->imtrx][eqn] & T_MASS) {
          mass += rho * ad_fv->turb_k_dot;
          mass *= wt_func * d_area;
          mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
        }
      }

      /* Assemble advection term */
      ADType adv = 0;
      for (int p = 0; p < VIM; p++) {
        adv += rho * (ad_fv->v[p] - ad_fv->x_dot[p]) * ad_fv->grad_turb_k[p];
      }
      adv *= wt_func * d_area;
      adv *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

      /* Assemble source terms */
      // Production - Destruction
      ADType src = Plim - beta_star * rho * omega * fv_old->turb_k;
      src *= -wt_func * d_area;
      src *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

      /* Assemble diffusion terms */
      ADType diff = 0.0;
      for (int p = 0; p < VIM; p++) {
        diff +=
            bf[eqn]->grad_phi[i][p] * (mu + mu_turb * sigma_k + vshock_k) * ad_fv->grad_turb_k[p];
        // + bf[eqn]->grad_phi[i][p] * vshock_k * gs_inner_dot[p];
      }
      diff *= d_area;
      diff *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];

      resid[0][i] -= mass + adv + src + diff;
      lec->R[LEC_R_INDEX(peqn, i)] -= mass.val() + adv.val() + src.val() + diff.val();
    } /* end of for (i=0,ei[pg->imtrx]->dofs...) */

    eqn = TURB_OMEGA;
    if (0 && fv->wall_distance < 0.005) {
      dbl d = std::max(0.0004, fv->wall_distance);
      dbl beta1 = 0.075;
      dbl nu = 1.57e-5;
      dbl ow = 10 * 6 * nu / (beta1 * d * d);
      peqn = upd->ep[pg->imtrx][eqn];
      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        ADType wt_func = bf[eqn]->phi[i];

        /* Assemble mass term */
        ADType mass = 0.0;
        if (pd->TimeIntegration != STEADY) {
          if (pd->e[pg->imtrx][eqn] & T_MASS) {
            mass += 1e1 * ad_fv->turb_omega * wt_func * d_area;
            mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
          }
        }

        /* Assemble advection term */
        ADType src = 0;
        for (int p = 0; p < VIM; p++) {
          src -= 1e1 * ow;
        }
        src *= wt_func * d_area;
        src *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

        /* Assemble diffusion terms */
        ADType diff = 0.0;
        // for (int p = 0; p < VIM; p++) {
        //   diff +=
        //       bf[eqn]->grad_phi[i][p] * (mu + mu_turb * sigma_omega) * ad_fv->grad_turb_omega[p];
        //   // + bf[eqn]->grad_phi[i][p] * vshock_w * gs_inner_dot[p];
        // }
        diff *= d_area;
        diff *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];

        resid[1][i] -= mass + src + diff;
        lec->R[LEC_R_INDEX(peqn, i)] -= mass.val() + src.val() + diff.val();
      } /* end of for (i=0,ei[pg->imtrx]->dofs...) */
    }
    {
      peqn = upd->ep[pg->imtrx][eqn];
      ADType Z_w = rho * ad_fv->turb_omega_dot;
      for (int i = 0; i < pd->Num_Dim; i++) {
        Z_w += rho * ad_fv->v[i] * ad_fv->grad_turb_omega[i];
      }

      coeff = 12 * (mu + mu_turb * sigma_omega) * (mu + mu_turb * sigma_omega);

      diff_g_g = 0;
      for (int i = 0; i < pd->Num_Dim; i++) {
        for (int j = 0; j < pd->Num_Dim; j++) {
          diff_g_g += coeff * G[i][j] * G[i][j];
        }
      }

      ADType norm_grad_omega = 0;
      for (int i = 0; i < pd->Num_Dim; i++) {
        norm_grad_omega = ad_fv->grad_turb_omega[i] * ad_fv->grad_turb_omega[i];
      }
      norm_grad_omega = sqrt(norm_grad_omega);
      for (int i = 0; i < pd->Num_Dim; i++) {
        r[i] = std::abs(ad_fv->grad_turb_omega[i]) / norm_grad_omega;
      }

      h_rgn = 0;
      for (int i = 0; i < pd->Num_Dim; i++) {
        for (int j = 0; j < ei[pg->imtrx]->dof[eqn]; j++) {
          h_rgn += std::abs(r[i] * bf[eqn]->grad_phi[j][i]);
        }
      }
      h_rgn = 2.0 / std::max(h_rgn, 1e-20);

      sugn3 = h_rgn * h_rgn / (4 * (mu + mu_turb * sigma_omega));
      // ADType tau_supg_w =
      // 1.0 / sqrt(1 / (sugn1 * sugn1) + 1 / (sugn2 * sugn2) + 1 / (sugn3 * sugn3));
      ADType tau_supg_w = 1.0 / (sqrt(4 / (dt * dt) + v_d_gv + diff_g_g));

      // ADType vshock_w = ad_yzbeta(pd->Num_Dim, TURB_OMEGA, upd->turbulent_info->omega_inf, Z_w,
      // ad_fv->grad_turb_omega);

      ADType grad_turb_w[DIM] = {fv_old->grad_turb_omega[0], fv_old->grad_turb_omega[1],
                                 fv_old->grad_turb_omega[2]};
      ADType vshock_w = 1 * ad_dcdd(pd->Num_Dim, TURB_OMEGA, grad_turb_w, gs_inner_dot);

      vshock_w = std::min(tau_supg_w, vshock_w);

      // if (fv->wall_distance < 1) {
      //   // larger diffusion near wall
      //   vshock_w = std::max(vshock_w, 1e-2*std::exp(-fv->wall_distance * 20));
      // }

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        ADType wt_func = bf[eqn]->phi[i];

        if (supg > 0) {
          if (supg != 0.0) {
            for (int p = 0; p < VIM; p++) {
              wt_func += supg * tau_supg_w * ad_fv->v[p] * bf[eqn]->grad_phi[i][p];
            }
          }
        }

        /* Assemble mass term */
        ADType mass = 0.0;
        if (pd->TimeIntegration != STEADY) {
          if (pd->e[pg->imtrx][eqn] & T_MASS) {
            mass += rho * ad_fv->turb_omega_dot * wt_func * d_area;
            mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
          }
        }

        /* Assemble advection term */
        ADType adv = 0;
        for (int p = 0; p < VIM; p++) {
          adv += rho * (ad_fv->v[p] - ad_fv->x_dot[p]) * ad_fv->grad_turb_omega[p];
        }
        adv *= wt_func * d_area;
        adv *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

        /* Assemble source terms */
        ADType src = (gamma * rho / mu_turb) * Plim - beta * rho * omega * omega +
                     2 * (1 - F1) * (rho * sigma_omega2 / omega) * grad_k_dot_grad_omega;
        // penalty for near wall
        // if (d < 0.4) {
        //   if (ad_fv->turb_omega > 2.422) {
        //   ADType red = 100*(ad_fv->turb_omega - 2.422);
        //   if (std::abs(red) > 10*std::abs(src)) {
        //     red = 10*std::abs(src) * red/std::abs(red);
        //   }
        //   src = -red;
        //   }
        // }
        src *= -wt_func * d_area;
        src *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

        /* Assemble diffusion terms */
        ADType diff = 0.0;
        for (int p = 0; p < VIM; p++) {
          diff += bf[eqn]->grad_phi[i][p] * (mu + mu_turb * sigma_omega + vshock_w) *
                  ad_fv->grad_turb_omega[p];
          // + bf[eqn]->grad_phi[i][p] * vshock_w * gs_inner_dot[p];
        }
        diff *= d_area;
        diff *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];

        resid[1][i] -= mass + adv + src + diff;
        lec->R[LEC_R_INDEX(peqn, i)] -= mass.val() + adv.val() + src.val() + diff.val();
      } /* end of for (i=0,ei[pg->imtrx]->dofs...) */
    }
  } /* end of if assemble residual */

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    eqn = TURB_OMEGA;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      for (int var = V_FIRST; var < V_LAST; var++) {

        /* Sensitivity w.r.t. velocity */
        if (pdv[var]) {
          int pvar = upd->vp[pg->imtrx][var];

          for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            // J = &(lec->J[LEC_J_INDEX(peqn, pvar, ii, 0)]);
            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += resid[1][i].dx(ad_fv->offset[var] + j);

          } /* End of loop over j */
        }   /* End of if the variale is active */
      }
    }
    eqn = TURB_K;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      for (int var = V_FIRST; var < V_LAST; var++) {

        /* Sensitivity w.r.t. velocity */
        if (pdv[var]) {
          int pvar = upd->vp[pg->imtrx][var];

          for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            // J = &(lec->J[LEC_J_INDEX(peqn, pvar, ii, 0)]);
            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += resid[0][i].dx(ad_fv->offset[var] + j);

          } /* End of loop over j */
        }   /* End of if the variale is active */
      }
    }
  } /* End of if assemble Jacobian */
  return (status);
}
int ad_assemble_invariant(double tt, /* parameter to vary time integration from
                                      * explicit (tt = 1) to implicit (tt = 0)    */
                          double dt) /*  time step size                          */
{
  int dim;
  int p;

  int eqn;
  int peqn;
  int i;
  int status;

  dbl h3 = fv->h3; /* Volume element (scale factors). */
  /*
   * Galerkin weighting functions for i-th and a-th momentum residuals
   * and some of their derivatives...
   */

  ADType wt_func;

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl wt;

  status = 0;

  /*
   * Unpack variables from structures for local convenience...
   */

  dim = pd->Num_Dim;

  /*
   * Bail out fast if there's nothing to do...
   */

  if (!pd->e[pg->imtrx][eqn = R_SHEAR_RATE]) {
    return (status);
  }

  peqn = upd->ep[pg->imtrx][eqn];

  wt = fv->wt; /* Numerical integration weight */

  ADType det_J = ad_fv->detJ; /* Really, ought to be mesh eqn. */

  ADType omega[DIM][DIM];
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      omega[a][b] = ad_fv->grad_v[a][b] + ad_fv->grad_v[b][a];
    }
  }

  ADType S = 0.;
  /* get gamma_dot invariant for viscosity calculations */
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      S += omega[a][b] * omega[a][b];
    }
  }

  if (S > 1e-20) {
    S = sqrt(0.5 * S);
  }

  /*
   * Residuals_________________________________________________________________
   */

  std::vector<ADType> resid(ei[pg->imtrx]->dof[eqn]);
  for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
    resid[i] = 0;
  }
  if (af->Assemble_Residual) {
    /*
     * Assemble the second_invariant equation
     */

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

      wt_func = bf[eqn]->phi[i];

      ADType advection = 0.;

      if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
        advection = -S;
        advection *= wt_func * det_J * wt * h3;
        advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
      }

      /*
       * Diffusion Term..
       */

      /* OK this really isn't a diffusion term.  Its really a
       * filtering term.  But it looks like a diffusion operator.No?
       */

      ADType diffusion = 0.;

      if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
        for (p = 0; p < dim; p++) {
          diffusion += ad_fv->basis[eqn].grad_phi[i][p] * ad_fv->grad_SH[p];
        }

        diffusion *= det_J * wt * h3;
        diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
      }

      /*
       * Source term...
       */

      ADType source = 0;

      if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
        source += ad_fv->SH;
        source *= wt_func * det_J * h3 * wt;
        source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
      }

      resid[i] += advection + source + diffusion;
      lec->R[LEC_R_INDEX(peqn, i)] += advection.val() + source.val() + diffusion.val();
    }
  }

  /*
   * Jacobian terms_________________________________________________________________
   */

  if (af->Assemble_Jacobian) {

    eqn = SHEAR_RATE;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      for (int var = V_FIRST; var < V_LAST; var++) {

        /* Sensitivity w.r.t. velocity */
        if (pd->v[pg->imtrx][var]) {
          int pvar = upd->vp[pg->imtrx][var];

          for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            // J = &(lec->J[LEC_J_INDEX(peqn, pvar, ii, 0)]);
            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += resid[i].dx(ad_fv->offset[var] + j);

          } /* End of loop over j */
        }   /* End of if the variale is active */
      }
    }
  } /* end of if(af,, */

  return (status);

} /* END of assemble_invariant */
#endif
#endif