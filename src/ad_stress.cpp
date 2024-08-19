#include <Sacado.hpp>

#include "ad_momentum.h"
#include "ad_turbulence.h"

extern "C" {

/* GOMA include files */
#include "ad_turbulence.h"
#include "density.h"
#include "el_elm.h"
#include "exo_struct.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "mm_fill_aux.h"
#include "mm_fill_ls.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_rs.h"
#include "mm_fill_solid.h"
#include "mm_fill_species.h"
#include "mm_fill_stabilization.h"
#include "mm_fill_stress.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "mm_qtensor_model.h"
#include "mm_viscosity.h"
#include "polymer_time_const.h"
#include "rf_allo.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_solver.h"
#include "std.h"
#include "user_mp.h"
}

ADType ad_vec_dot(const int n1, ADType *v1, ADType *v2) {
  int i;
  ADType rc = 0.0;

  for (i = 0; i < n1; i++) {
    rc += *v1 * *v2;
    v1++;
    v2++;
  }
  return (rc);
}
int ad_tensor_dot(ADType t1[DIM][DIM],
                  ADType t2[DIM][DIM],
                  ADType t1_dot_t2[DIM][DIM],
                  const int dim) {
  int i, j, k;
  int status;
  ADType v1[DIM];
  ADType v2[DIM];

  for (k = 0; k < dim; k++) {
    for (i = 0; i < dim; i++) {
      v1[i] = t1[k][i];
    }
    for (j = 0; j < dim; j++) {
      for (i = 0; i < dim; i++) {
        v2[i] = t2[i][j];
      }
      t1_dot_t2[k][j] = ad_vec_dot(dim, v1, v2);
    }
  }

  status = 1;
  return (status);
}
void ad_compute_a_dot_b(ADType b[DIM][DIM], ADType G[DIM][DIM], ADType a_dot_b[DIM][DIM]) {

  if (VIM == 2) {

    ADType a12 =
        ((b[0][1] * G[0][0] - b[0][0] * G[0][1]) + (b[1][1] * G[1][0] - b[1][0] * G[1][1])) /
        (b[0][0] + b[1][1] + 1e-16);

    ADType a[DIM][DIM] = {{0., a12, 0.}, {-a12, 0., 0.}, {0., 0., 0.}};

    ad_tensor_dot(a, b, a_dot_b, VIM);

  } else { // VIM = 3
    ADType D =
        -b[0][1] * (b[0][1] * (b[0][0] + b[1][1]) + b[0][2] * b[1][2]) -
        b[0][2] * (b[0][1] * b[1][2] + b[0][2] * (b[0][0] + b[2][2])) +
        (b[1][1] + b[2][2]) * (-b[1][2] * b[1][2] + (b[0][0] + b[1][1]) * (b[0][0] + b[2][2])) +
        1e-16;
    ADType invD = 1.0 / D;

    ADType a12 = invD * (-pow(b[0][1], 2) + (b[0][0] + b[2][2]) * (b[1][1] + b[2][2])) *
                     (G[0][0] * b[0][1] - G[0][1] * b[0][0] + G[1][0] * b[1][1] -
                      G[1][1] * b[0][1] + G[2][0] * b[1][2] - G[2][1] * b[0][2]) -
                 (b[0][1] * b[0][2] + b[1][2] * (b[1][1] + b[2][2])) *
                     (G[0][0] * b[0][2] - G[0][2] * b[0][0] + G[1][0] * b[1][2] -
                      G[1][2] * b[0][1] + G[2][0] * b[2][2] - G[2][2] * b[0][2]) +
                 (b[0][1] * b[1][2] + b[0][2] * (b[0][0] + b[2][2])) *
                     (G[0][1] * b[0][2] - G[0][2] * b[0][1] + G[1][1] * b[1][2] -
                      G[1][2] * b[1][1] + G[2][1] * b[2][2] - G[2][2] * b[1][2]);

    ADType a13 =
        invD *
        ((-pow(b[0][2], 2) + (b[0][0] + b[1][1]) * (b[1][1] + b[2][2])) *
             (G[0][0] * b[0][2] - G[0][2] * b[0][0] + G[1][0] * b[1][2] - G[1][2] * b[0][1] +
              G[2][0] * b[2][2] - G[2][2] * b[0][2]) +
         (-b[0][1] * b[0][2] - b[1][2] * (b[1][1] + b[2][2])) *
             (G[0][0] * b[0][1] - G[0][1] * b[0][0] + G[1][0] * b[1][1] - G[1][1] * b[0][1] +
              G[2][0] * b[1][2] - G[2][1] * b[0][2]) -
         (b[0][1] * (b[0][0] + b[1][1]) + b[0][2] * b[1][2]) *
             (G[0][1] * b[0][2] - G[0][2] * b[0][1] + G[1][1] * b[1][2] - G[1][2] * b[1][1] +
              G[2][1] * b[2][2] - G[2][2] * b[1][2])) /
        (-b[0][1] * (b[0][1] * (b[0][0] + b[1][1]) + b[0][2] * b[1][2]) -
         b[0][2] * (b[0][1] * b[1][2] + b[0][2] * (b[0][0] + b[2][2])) +
         (b[1][1] + b[2][2]) * (-pow(b[1][2], 2) + (b[0][0] + b[1][1]) * (b[0][0] + b[2][2])));

    ADType a23 = invD * (-pow(b[1][2], 2) + (b[0][0] + b[1][1]) * (b[0][0] + b[2][2])) *
                     (G[0][1] * b[0][2] - G[0][2] * b[0][1] + G[1][1] * b[1][2] -
                      G[1][2] * b[1][1] + G[2][1] * b[2][2] - G[2][2] * b[1][2]) +
                 (b[0][1] * b[1][2] + b[0][2] * (b[0][0] + b[2][2])) *
                     (G[0][0] * b[0][1] - G[0][1] * b[0][0] + G[1][0] * b[1][1] -
                      G[1][1] * b[0][1] + G[2][0] * b[1][2] - G[2][1] * b[0][2]) -
                 (b[0][1] * (b[0][0] + b[1][1]) + b[0][2] * b[1][2]) *
                     (G[0][0] * b[0][2] - G[0][2] * b[0][0] + G[1][0] * b[1][2] -
                      G[1][2] * b[0][1] + G[2][0] * b[2][2] - G[2][2] * b[0][2]);

    ADType a[DIM][DIM] = {
        {0.0, a12, a13},
        {-a12, 0.0, a23},
        {-a13, -a23, 0.0},
    };

    ad_tensor_dot(a, b, a_dot_b, VIM);
  }
}

int ad_sqrt_conf_source(int mode, ADType b[DIM][DIM], ADType source_term[DIM][DIM]) {
  ADType binv[DIM][DIM];
  ADType d_binv_db[DIM][DIM][DIM][DIM];
  if (VIM == 2) {
    ADType det = b[0][0] * b[1][1] - b[0][1] * b[0][1] + 1e-16;
    binv[0][0] = b[1][1] / det;
    binv[0][1] = -b[0][1] / det;
    binv[1][0] = -b[0][1] / det;
    binv[1][1] = b[0][0] / det;

  } else if (VIM == 3) {
    ADType det = b[0][0] * (b[1][1] * b[2][2] - b[1][2] * b[2][1]) -
                 b[0][1] * (b[1][0] * b[2][2] - b[2][0] * b[1][2]) +
                 b[0][2] * (b[1][0] * b[2][1] - b[2][0] * b[1][1]) + 1e-16;

    binv[0][0] = (b[1][1] * b[2][2] - b[2][1] * b[1][2]) / (det);

    binv[0][1] = -(b[0][1] * b[2][2] - b[2][1] * b[0][2]) / (det);

    binv[0][2] = (b[0][1] * b[1][2] - b[1][1] * b[0][2]) / (det);

    binv[1][0] = -(b[1][0] * b[2][2] - b[2][0] * b[1][2]) / (det);

    binv[1][1] = (b[0][0] * b[2][2] - b[2][0] * b[0][2]) / (det);

    binv[1][2] = -(b[0][0] * b[1][2] - b[1][0] * b[0][2]) / (det);

    binv[2][0] = (b[1][0] * b[2][1] - b[1][1] * b[2][0]) / (det);

    binv[2][1] = -(b[0][0] * b[2][1] - b[2][0] * b[0][1]) / (det);

    binv[2][2] = (b[0][0] * b[1][1] - b[1][0] * b[0][1]) / (det);

  } else {
    GOMA_EH(GOMA_ERROR, "Unknown VIM = %d for SQRT conformation tensor", VIM);
  }

  switch (vn->ConstitutiveEquation) {
  case OLDROYDB: {
    for (int ii = 0; ii < VIM; ii++) {
      for (int jj = 0; jj < VIM; jj++) {
        source_term[ii][jj] = -0.5 * (binv[ii][jj] - b[ii][jj]);
      }
    }
  } break;
  case PTT: {

    ADType trace = 0;
    for (int i = 0; i < VIM; i++) {
      for (int j = 0; j < VIM; j++) {
        trace += b[i][j] * b[i][j];
      }
    }

    ADType Z = 1.0;
    ADType dZ_dtrace = 0;

    // PTT exponent
    eps = ve[mode]->eps;

    if (vn->ptt_type == PTT_LINEAR) {
      Z = 1 + eps * (trace - (double)VIM);
    } else if (vn->ptt_type == PTT_EXPONENTIAL) {
      const double exp_max = 700;
      ADType inner = eps * (trace - (double)VIM);
      if ((inner > exp_max) || (inner < -exp_max)) {
        GOMA_WH_MANY(GOMA_ERROR, "Exponential overflow in PTT_EXPONENTIAL");
        return GOMA_ERROR;
      }
      Z = exp(eps * (trace - (double)VIM));
    } else {
      GOMA_EH(GOMA_ERROR, "Unrecognized PTT Form %d", vn->ptt_type);
    }

    for (int ii = 0; ii < VIM; ii++) {
      for (int jj = 0; jj < VIM; jj++) {
        source_term[ii][jj] = -0.5 * Z * (binv[ii][jj] - b[ii][jj]);
      }
    }
  } break;
  default:
    GOMA_EH(GOMA_ERROR, "Unknown Constitutive equation form for SQRT_CONF");
    break;
  }

  return GOMA_SUCCESS;
}
void ad_load_modal_pointers(
    int ve_mode, /* mode number */
    dbl tt,
    dbl dt,
    ADType s[DIM][DIM],           /* stress tensor for mode ve_mode */
    ADType s_dot[DIM][DIM],       /* stress tensor time derivative for mode ve_mode */
    ADType grad_s[DIM][DIM][DIM]) /* grad of stress tensor for mode ve_mode */

{
  int a, b, p; /* indeces for dimensions */
  /* load up things we need in the assembly routine for each mode in turn*/

  /* put stress in a nice working array */

  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      s[a][b] = ad_fv->S[ve_mode][a][b];
      if (pd->TimeIntegration != STEADY) {
        s_dot[a][b] = ad_fv->S_dot[ve_mode][a][b];
      } else {
        s_dot[a][b] = 0.;
      }
    }
  }

  for (p = 0; p < VIM; p++) {
    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        grad_s[p][a][b] = ad_fv->grad_S[ve_mode][p][a][b];
      }
    }
  }
}

int ad_assemble_stress_sqrt_conf(dbl tt, /* parameter to vary time integration from
                                          * explicit (tt = 1) to implicit (tt = 0) */
                                 dbl dt, /* current time step size */
                                 PG_DATA *pg_data) {
  int dim, q, w;

  int eqn;
  int peqn, pvar;
  int evss_gradv = 0;

  int i, j, status, mode;
  ADType v[DIM];     /* Velocity field. */
  ADType x_dot[DIM]; /* current position field derivative wrt time. */
  dbl h3;            /* Volume element (scale factors). */

  ADType grad_v[DIM][DIM];
  ADType gamma[DIM][DIM]; /* Shear-rate tensor based on velocity */
  dbl dgamma[DIM][DIM];   /* Shear-rate tensor based on velocity */
  ADType det_J;           /* determinant of element Jacobian */

  int err;
  dbl alpha = 0;  /* This is the Geisekus mobility parameter */
  dbl lambda = 0; /* polymer relaxation constant */
  double xi;
  double d_xi_dF[MDE];
  dbl eps = 0; /* This is the PTT elongation parameter */
  double d_eps_dF[MDE];
  /*
   *
   * Note how carefully we avoid refering to d(phi[i])/dx[j] and refer instead
   * to the j-th component of grad_phi[j][i] so that this vector can be loaded
   * up with components that may be different in non Cartesian coordinate
   * systems.
   *
   * We will, however, insist on *orthogonal* coordinate systems, even if we
   * might permit them to be curvilinear.
   *
   * Assume all components of velocity are interpolated with the same kind
   * of basis function.
   */

  /*
   * Petrov-Galerkin weighting functions for i-th and ab-th stress residuals
   * and some of their derivatives...
   */

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl wt;

  /* Variables for stress */

  int R_s[MAX_MODES][DIM][DIM];
  int v_s[MAX_MODES][DIM][DIM];

  ADType b[DIM][DIM];     /* stress tensor */
  ADType b_dot[DIM][DIM]; /* stress tensor from last time step */
  ADType grad_b[DIM][DIM][DIM];

  ADType g[DIM][DIM];  /* velocity gradient tensor */
  ADType gt[DIM][DIM]; /* transpose of velocity gradient tensor */

  /* dot product tensors */

  ADType b_dot_g[DIM][DIM];

  const bool saramitoEnabled =
      (vn->ConstitutiveEquation == SARAMITO_OLDROYDB || vn->ConstitutiveEquation == SARAMITO_PTT ||
       vn->ConstitutiveEquation == SARAMITO_GIESEKUS);

  if (saramitoEnabled) {
    GOMA_EH(GOMA_ERROR, "Saramito not available for SQRT_CONF");
  }

  /*  shift function */
  dbl at = 0.0;
  dbl wlf_denom;

  /* advective terms are precalculated */
  ADType v_dot_del_b[DIM][DIM];
  ADType x_dot_del_b[DIM][DIM];

  /* SUPG variables */
  dbl supg = 0;

  if (vn->evssModel == EVSS_GRADV) {
    evss_gradv = 1;
  }

  status = 0;

  eqn = R_STRESS11;

  /*
   * Bail out fast if there's nothing to do...
   */

  if (!pd->e[pg->imtrx][eqn]) {
    return (status);
  }

  /*
   * Unpack variables from structures for local convenience...
   */

  dim = pd->Num_Dim;

  wt = fv->wt;

  det_J = ad_fv->detJ; /* Really, ought to be mesh eqn. */

  h3 = fv->h3; /* Differential volume element (scales). */

  /* load eqn and variable number in tensor form */
  (void)stress_eqn_pointer(v_s);
  (void)stress_eqn_pointer(R_s);

  /*
   * Field variables...
   */
  for (int a = 0; a < WIM; a++) {
    v[a] = ad_fv->v[a];

    /* note, these are zero for steady calculations */
    x_dot[a] = 0.0;
    if (pd->TimeIntegration != STEADY && pd->gv[MESH_DISPLACEMENT1 + a]) {
      x_dot[a] = ad_fv->x_dot[a];
    }
  }

  /*
   * In Cartesian coordinates, this velocity gradient tensor will
   * have components that are...
   *
   * 			grad_v[a][b] = d v_b
   *				       -----
   *				       d x_a
   */

  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      grad_v[a][b] = ad_fv->grad_v[a][b];
    }
  }

  /* load up shearrate tensor based on velocity */
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      gamma[a][b] = grad_v[a][b] + grad_v[b][a];
      dgamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
    }
  }

  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      if (evss_gradv) {
        g[a][b] = ad_fv->grad_v[a][b];
        gt[a][b] = ad_fv->grad_v[b][a];
      } else {
        g[a][b] = ad_fv->G[a][b];
        gt[b][a] = g[a][b];
      }
    }
  }

  if (vn->wt_funcModel == GALERKIN) {
    supg = 0.;
  } else if (vn->wt_funcModel == SUPG) {
    supg = vn->wt_func;
  }

  ADType supg_tau;
  if (supg != 0.) {
    ad_supg_tau_shakib(supg_tau, dim, dt, 1e-8, eqn);
  }
  /* end Petrov-Galerkin addition */
  dbl yzbeta_factor = 0.0;
  dbl beta[2] = {1.0, 2.0};
  if (vn->shockcaptureModel == SC_YZBETA) {
    yzbeta_factor = vn->shockcapture;
  } else if (vn->shockcaptureModel == SC_DCDD) {
  } else if (vn->shockcaptureModel != SC_NONE) {
    GOMA_EH(GOMA_ERROR, "Unknown shock capture model, only YZBETA supported for SQRT_CONF");
  }

  /*  shift factor  */
  if (pd->gv[TEMPERATURE]) {
    if (vn->shiftModel == CONSTANT) {
      at = vn->shift[0];
    } else if (vn->shiftModel == MODIFIED_WLF) {
      wlf_denom = vn->shift[1] + fv->T - mp->reference[TEMPERATURE];
      if (wlf_denom != 0.) {
        at = exp(vn->shift[0] * (mp->reference[TEMPERATURE] - fv->T) / wlf_denom);
      } else {
        at = 1.;
      }
    }
  } else {
    at = 1.;
  }

  /* Begin loop over modes */
  for (mode = 0; mode < vn->modes; mode++) {

    ad_load_modal_pointers(mode, tt, dt, b, b_dot, grad_b);

    /* precalculate advective terms of form (v dot del tensor)*/

    /*
     * Stress tensor...(Note "anti-BSL" sign convention on deviatoric stress)
     */
    for (int ii = 0; ii < VIM; ii++) {
      for (int jj = 0; jj < VIM; jj++) {
        v_dot_del_b[ii][jj] = 0.;
        x_dot_del_b[ii][jj] = 0.;
        for (q = 0; q < WIM; q++) {
          v_dot_del_b[ii][jj] += v[q] * grad_b[q][ii][jj];
          x_dot_del_b[ii][jj] += x_dot[q] * grad_b[q][ii][jj];
        }
      }
    }

    if (saramitoEnabled == TRUE) {
      GOMA_EH(GOMA_ERROR, "Saramito not enabled sqrt");
    }

    double d_alpha_dF[MDE];
    /* get Geisekus mobility parameter */
    if (ve[mode]->alphaModel == CONSTANT) {
      alpha = ve[mode]->alpha;
    } else if (ls != NULL && ve[mode]->alphaModel == VE_LEVEL_SET) {
      double pos_alpha = ve[mode]->pos_ls.alpha;
      double neg_alpha = ve[mode]->alpha;
      double width = ls->Length_Scale;
      err = level_set_property(neg_alpha, pos_alpha, width, &alpha, d_alpha_dF);
      GOMA_EH(err, "level_set_property() failed for mobility parameter.");
    } else {
      GOMA_EH(GOMA_ERROR, "Unknown mobility parameter model");
    }

    /* get time constant */
    lambda = polymer_time_const(ve[mode]->time_const_st, dgamma, NULL);

    xi = 0;
    if (ve[mode]->xiModel == CONSTANT) {
      xi = ve[mode]->xi;
    } else if (ls != NULL && ve[mode]->xiModel == VE_LEVEL_SET) {
      double pos_xi = ve[mode]->pos_ls.xi;
      double neg_xi = ve[mode]->xi;
      double width = ls->Length_Scale;
      err = level_set_property(neg_xi, pos_xi, width, &xi, d_xi_dF);
      GOMA_EH(err, "level_set_property() failed for ptt xi parameter.");
    } else {
      GOMA_EH(GOMA_ERROR, "Unknown PTT Xi parameter model");
    }

    if (DOUBLE_NONZERO(xi)) {
      GOMA_EH(GOMA_ERROR, "PTT Xi parameter currently required to be 0 for SQRT_CONF");
    }

    if (ve[mode]->epsModel == CONSTANT) {
      eps = ve[mode]->eps;
    } else if (ls != NULL && ve[mode]->epsModel == VE_LEVEL_SET) {
      double pos_eps = ve[mode]->pos_ls.eps;
      double neg_eps = ve[mode]->eps;
      double width = ls->Length_Scale;
      err = level_set_property(neg_eps, pos_eps, width, &eps, d_eps_dF);
      GOMA_EH(err, "level_set_property() failed for ptt epsilon parameter.");
    } else {
      GOMA_EH(GOMA_ERROR, "Unknown PTT Epsilon parameter model");
    }

    (void)ad_tensor_dot(b, g, b_dot_g, VIM);

    ADType a_dot_b[DIM][DIM];

    ad_compute_a_dot_b(b, g, a_dot_b);

    ADType source_term[DIM][DIM];
    goma_error err = ad_sqrt_conf_source(mode, b, source_term);
    if (err) {
      return 1;
    }
    std::vector<std::vector<std::vector<ADType>>> resid(VIM);
    for (int a = 0; a < VIM; a++) {
      resid[a].resize(VIM);
      for (int b = 0; b < VIM; b++) {
        resid[a][b].resize(ei[pg->imtrx]->dof[v_s[mode][a][b]]);
        for (i = 0; i < ei[pg->imtrx]->dof[v_s[mode][a][b]]; i++) {
          resid[a][b][i] = 0;
        }
      }
    }
    /*
     * Residuals_________________________________________________________________
     */

    if (af->Assemble_Residual) {
      /*
       * Assemble each component "ab" of the polymer stress equation...
       */
      for (int ii = 0; ii < VIM; ii++) {
        for (int jj = 0; jj < VIM; jj++) {

          if (ii <= jj) /* since the stress tensor is symmetric, only assemble the upper half */
          {
            eqn = R_s[mode][ii][jj];

            /*
             * In the element, there will be contributions to this many equations
             * based on the number of degrees of freedom...
             */

            for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
              ADType wt_func = bf[eqn]->phi[i];
              /* add Petrov-Galerkin terms as necessary */
              if (supg != 0.) {
                for (w = 0; w < dim; w++) {
                  // wt_func += supg * supg_tau * ad_fv->v[w];// * ad_fv->basis[eqn].grad_phi[i][w];
                  wt_func += supg * supg_tau * ad_fv->v[w] * ad_fv->basis[eqn].grad_phi[i][w];
                }
              }

              ADType mass = 0.;

              if (pd->TimeIntegration != STEADY) {
                if (pd->e[pg->imtrx][eqn] & T_MASS) {
                  mass = b_dot[ii][jj];
                  mass *= wt_func * at * lambda * det_J * wt;
                  mass *= h3;
                  mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                }
              }

              ADType advection = 0.;
              if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                if (DOUBLE_NONZERO(lambda)) {

                  advection += v_dot_del_b[ii][jj] - x_dot_del_b[ii][jj];
                  advection -= b_dot_g[ii][jj];
                  advection -= a_dot_b[ii][jj];
                  advection *= wt_func * at * lambda * det_J * wt * h3;
                  advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                }
              }

              ADType diffusion = 0.;
              if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                if (vn->shockcaptureModel == SC_YZBETA) {
                  ADType Z = b_dot[ii][jj];
                  Z += 1e-32 + v_dot_del_b[ii][jj] - x_dot_del_b[ii][jj];
                  Z -= b_dot_g[ii][jj];
                  Z -= a_dot_b[ii][jj];
                  Z *= at * lambda;
                  Z += source_term[ii][jj];

                  ADType Y_inv = 1.0;
                  ADType hdc = 0;
                  ADType js = 0;
                  for (int k = 0; k < ei[pg->imtrx]->dof[eqn]; k++) {
                    for (int p = 0; p < VIM; p++) {
                      js += fabs(grad_b[p][ii][jj] * ad_fv->basis[eqn].grad_phi[k][p]);
                    }
                  }
                  hdc = 1 / (js + 1e-16);

                  ADType inner = 1e-32;
                  for (int p = 0; p < VIM; p++) {
                    inner += Y_inv * grad_b[p][ii][jj] * grad_b[p][ii][jj];
                  }

                  ADType kdc = 0;
                  for (int ib = 0; ib < 1; ib++) {
                    ADType bt = beta[ib];
                    if (bt > 1.0) {
                      kdc += fabs(Y_inv * Z) * pow(hdc, bt) * pow(fabs(b[ii][jj]), 1 - bt);
                    } else {
                      kdc += fabs(Y_inv * Z) * hdc * pow(inner, bt / 2 - 1);
                    }
                    // kdc += pow(inner, bt / 2 - 1);
                  }
                  kdc *= 0.5;
                  for (int r = 0; r < VIM; r++) {
                    diffusion += kdc * grad_b[r][ii][jj] * ad_fv->basis[eqn].grad_phi[i][r];
                  }

                  diffusion *= yzbeta_factor * det_J * wt * h3;
                  diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                } else if (vn->shockcaptureModel == SC_DCDD) {
                  dbl dcdd_factor = vn->shockcapture;
                  ADType tmp = 0.0;
                  ADType s[DIM] = {0.0};
                  ADType r[DIM] = {0.0};
                  for (int w = 0; w < dim; w++) {
                    tmp += (v[w] - x_dot[w]) * (v[w] - x_dot[w]);
                  }
                  tmp = 1.0 / (sqrt(tmp + 1e-32));
                  for (int w = 0; w < dim; w++) {
                    s[w] = (v[w] - x_dot[w]) * tmp;
                  }
                  ADType mags = 0;
                  for (int w = 0; w < dim; w++) {
                    mags += (grad_b[w][ii][jj] * grad_b[w][ii][jj]);
                  }
                  mags = 1.0 / (sqrt(mags + 1e-32));
                  for (int w = 0; w < dim; w++) {
                    r[w] = grad_b[w][ii][jj] * mags;
                  }

                  ADType he = 0.0;
                  for (int q = 0; q < ei[pg->imtrx]->dof[eqn]; q++) {
                    ADType tmp = 0;
                    for (int w = 0; w < dim; w++) {
                      tmp += ad_fv->basis[eqn].grad_phi[q][w] * ad_fv->basis[eqn].grad_phi[q][w];
                    }
                    he += 1.0 / sqrt(tmp);
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
                    magv += v[q] * v[q];
                  }
                  magv = sqrt(magv + 1e-32);

                  ADType tau_dcdd = 0.5 * he * (1.0 / (mags + 1e-16)) * hrgn * hrgn;
                  // ADType tau_dcdd = he * (1.0 / mags) * hrgn * hrgn / lambda;
                  // printf("%g %g ", supg_tau.val(), tau_dcdd.val());
                  // tau_dcdd = 1 / sqrt(1.0 / (supg_tau * supg_tau + 1e-32) +
                  //                     1.0 / (tau_dcdd * tau_dcdd + 1e-32));
                  tau_dcdd = std::min(supg_tau, tau_dcdd);
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

                  ADType gs_inner_dot[DIM] = {0.0};
                  for (int w = 0; w < dim; w++) {
                    ADType tmp = 0.;
                    for (int z = 0; z < dim; z++) {
                      tmp += grad_b[w][ii][jj] * inner_tensor[w][z];
                    }
                    gs_inner_dot[w] = tmp;
                    // gs_inner_dot[w] = grad_s[w][ii][jj];
                  }

                  for (int w = 0; w < dim; w++) {
                    // diffusion += tau_dcdd * grad_b[w][ii][jj] * ad_fv->basis[eqn].grad_phi[i][w];
                    diffusion += tau_dcdd * gs_inner_dot[w] * ad_fv->basis[eqn].grad_phi[i][w];
                  }
                  diffusion *= dcdd_factor * det_J * wt * h3;
                }
              }

              /*
               * Source term...
               */

              ADType source = 0.;
              if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                // consider whether saramitoCoeff should multiply here
                source = source_term[ii][jj];
                source *= wt_func * det_J * h3 * wt;

                source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
              }

              /*
               * Add contributions to this residual (globally into Resid, and
               * locally into an accumulator)
               */

              lec->R[LEC_R_INDEX(upd->ep[pg->imtrx][eqn], i)] +=
                  mass.val() + advection.val() + diffusion.val() + source.val();
              resid[ii][jj][i] += mass + advection + diffusion + source;
            }
          }
        }
      }
    }

    /*
     * Jacobian terms...
     */

    if (af->Assemble_Jacobian) {
      for (int ii = 0; ii < VIM; ii++) {
        for (int jj = 0; jj < VIM; jj++) {
          if (ii <= jj) /* since the stress tensor is symmetric, only assemble the upper half */
          {
            eqn = R_s[mode][ii][jj];
            peqn = upd->ep[pg->imtrx][eqn];

            for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
              for (int var = V_FIRST; var < V_LAST; var++) {
                if (pd->v[pg->imtrx][var]) {
                  pvar = upd->vp[pg->imtrx][var];
                  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                    lec->J[LEC_J_INDEX(peqn, pvar, i, j)] +=
                        resid[ii][jj][i].dx(ad_fv->offset[var] + j);
                  }
                }
              }
            }
          }
        }
      }
    }

  } /* End loop over modes */

  return (status);
}
