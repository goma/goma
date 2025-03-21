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

/* This routine contains all the viscosity models along with the routine
 * that calls them. If you add a new viscosity model, it would go in here.
 */

/* Standard include files */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* GOMA include files */

#include "density.h"
#include "el_elm.h"
#include "mm_as.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "mm_fill_ls.h"
#include "mm_fill_terms.h"
#include "mm_fill_util.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "polymer_time_const.h"
#include "rf_allo.h"
#include "rf_bc_const.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_solver.h"
#include "std.h"
#include "user_mp.h"
#include "user_mp_gen.h"

#define GOMA_MM_VISCOSITY_C

double polymer_time_const(struct PolymerTimeConstants *lambda_st,
                          dbl gamma_dot[DIM][DIM],
                          POLYMER_TIME_CONST_DEPENDENCE_STRUCT *d_lambda) {
  int err;
  dbl lambda = 0.;

  /* Zero out sensitivities */
  zeroStructures(d_lambda, 1);
  if (lambda_st->ConstitutiveEquation == CONSTANT) {
    return lambda_st->lambda0;
  }

  switch (lambda_st->ConstitutiveEquation) {
  case CONSTANT:
    return lambda_st->lambda0;
  case POWER_LAW:
    return power_law_time_const(lambda_st, gamma_dot, d_lambda);
  case CARREAU:
    return carreau_polymer_time_const(lambda_st, gamma_dot, d_lambda);
  case VE_LEVEL_SET:
    if (ls != NULL) {
      double pos_mup = lambda_st->pos_ls_lambda;
      double neg_mup = lambda_st->lambda0;
      double width = ls->Length_Scale;
      if (d_lambda != NULL)
        err = level_set_property(neg_mup, pos_mup, width, &lambda, d_lambda->F);
      else
        err = level_set_property(neg_mup, pos_mup, width, &lambda, NULL);
      GOMA_EH(err, "level_set_property() failed for polymer viscosity.");
      return lambda;
    } else {
      GOMA_EH(GOMA_ERROR, "Expected level set to be enabled for VE_LEVEL_SET");
    }
    break;
  default:
    GOMA_EH(GOMA_ERROR, "Invalid polymer time constant model");
    break;
  }
  return lambda;
}

double power_law_time_const(struct PolymerTimeConstants *lambda_st,
                            dbl gamma_dot[DIM][DIM], /* strain rate tensor */
                            POLYMER_TIME_CONST_DEPENDENCE_STRUCT *d_lambda) {

  int a, b;
  int mdofs = 0, vdofs;

  int i, j;

  dbl gammadot; /* strain rate invariant */

  dbl d_gd_dv[DIM][MDE];    /* derivative of strain rate invariant
                               wrt velocity */
  dbl d_gd_dmesh[DIM][MDE]; /* derivative of strain rate invariant
                               wrt mesh */

  dbl val;
  dbl lambda0;
  dbl nexp;
  dbl offset;
  dbl lambda = 0.;

  vdofs = ei[pg->imtrx]->dof[VELOCITY1];

  if (pd->e[pg->imtrx][R_MESH1]) {
    mdofs = ei[pg->imtrx]->dof[R_MESH1];
  }

  calc_shearrate(&gammadot, gamma_dot, d_gd_dv, d_gd_dmesh);

  /* calculate power law viscosity
     mu = mu0 * (offset + gammadot)**(nexp-1))                 */
  lambda0 = lambda_st->lambda0;
  nexp = lambda_st->nexp;
  offset = 0.00001;

  val = pow(gammadot + offset, nexp - 1.);
  lambda = lambda0 * val;
  /*
   * d( mu )/dmesh
   */

  val = pow(gammadot + offset, nexp - 2.);

  if (d_lambda != NULL) {
    d_lambda->gd = lambda0 * (nexp - 1.0) * val;
  }

  if (d_lambda != NULL && pd->e[pg->imtrx][R_MESH1]) {
    for (b = 0; b < VIM; b++) {
      for (j = 0; j < mdofs; j++) {
        if (DOUBLE_NONZERO(gammadot) && Include_Visc_Sens) {

          d_lambda->X[b][j] = d_lambda->gd * d_gd_dmesh[b][j];

        } else {
          /* printf("\ngammadot is zero in viscosity function");*/
          d_lambda->X[b][j] = 0.0;
        }
      }
    }
  }

  /*
   * d( mu )/dv
   */
  if (d_lambda != NULL && pd->e[pg->imtrx][R_MOMENTUM1]) {
    for (a = 0; a < VIM; a++) {
      for (i = 0; i < vdofs; i++) {
        if (DOUBLE_NONZERO(gammadot) && Include_Visc_Sens) {
          d_lambda->v[a][i] = d_lambda->gd * d_gd_dv[a][i];
        } else {
          d_lambda->v[a][i] = 0.0;
        }
      }
    }
  }

  return (lambda);
}

double carreau_polymer_time_const(struct PolymerTimeConstants *lambda_st,
                                  dbl gamma_dot[DIM][DIM], /* strain rate tensor */
                                  POLYMER_TIME_CONST_DEPENDENCE_STRUCT *d_lambda) {

  int a, b;
  int mdofs = 0, vdofs;

  int i, j;

  dbl gammadot; /* strain rate invariant */

  dbl d_gd_dv[DIM][MDE];    /* derivative of strain rate invariant
                               wrt velocity */
  dbl d_gd_dmesh[DIM][MDE]; /* derivative of strain rate invariant
                               wrt mesh */

  dbl val, val1, val2;
  dbl lambda = 0.;
  dbl lambda0;
  dbl lambdainf;
  dbl nexp;
  dbl aexp;
  dbl carreau_lambda;

  vdofs = ei[pg->imtrx]->dof[VELOCITY1];

  if (pd->e[pg->imtrx][R_MESH1]) {
    mdofs = ei[pg->imtrx]->dof[R_MESH1];
  }

  calc_shearrate(&gammadot, gamma_dot, d_gd_dv, d_gd_dmesh);

  lambda0 = lambda_st->lambda0;
  nexp = lambda_st->nexp;
  lambdainf = lambda_st->lambdainf;
  aexp = lambda_st->aexp;
  carreau_lambda = lambda_st->carreau_lambda;

  if (DOUBLE_NONZERO(gammadot)) {
    val2 = pow(carreau_lambda * gammadot, aexp);
  } else {
    val2 = 0.;
  }
  val = pow(1. + val2, (nexp - 1.) / aexp);
  lambda = lambdainf + (lambda0 - lambdainf) * val;

  /* gammadot = 0.0; */
  /* this effectively turns off the viscosity Jac terms */

  if (DOUBLE_NONZERO(gammadot)) {
    val = pow(carreau_lambda * gammadot, aexp - 1.);
  } else {
    val = 0.;
  }
  val1 = pow(1. + val2, (nexp - 1. - aexp) / aexp);

  if (d_lambda != NULL)
    d_lambda->gd = (lambda0 - lambdainf) * (nexp - 1.) * carreau_lambda * val * val1;

  /*
   * d( lambda )/dmesh
   */
  if (d_lambda != NULL && pd->e[pg->imtrx][R_MESH1]) {
    for (b = 0; b < VIM; b++) {
      for (j = 0; j < mdofs; j++) {
        if (DOUBLE_NONZERO(gammadot) && Include_Visc_Sens) {
          d_lambda->X[b][j] = d_lambda->gd * d_gd_dmesh[b][j];
        } else {
          /* printf("\ngammadot is zero in viscosity function");*/
          d_lambda->X[b][j] = 0.0;
        }
      }
    }
  }

  /*
   * d( lambda )/dv
   */
  if (d_lambda != NULL && pd->e[pg->imtrx][R_MOMENTUM1]) {
    for (a = 0; a < VIM; a++) {
      for (i = 0; i < vdofs; i++) {
        if (DOUBLE_NONZERO(gammadot) && Include_Visc_Sens) {
          d_lambda->v[a][i] = d_lambda->gd * d_gd_dv[a][i];
        } else {
          d_lambda->v[a][i] = 0.0;
        }
      }
    }
  }
  return (lambda);
}
