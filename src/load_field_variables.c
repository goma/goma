/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2023 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/

#include "load_field_variables.h"
#include "bc_colloc.h"
#include "el_geom.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_shell.h"
#include "mm_fill_stress.h"
#include "mm_mp.h"
#include "mm_post_def.h"
#include "rf_fem.h"
#include "std.h"

/***************************************************************************/
/****************************************************************************/
/****************************************************************************/

void scalar_fv_fill(double **base,
                    double **base_dot,
                    double **base_old,
                    double *phiv,
                    int dofs,
                    double *val,
                    double *val_dot,
                    double *val_old)
/*
 *   scalar_fv_fill:   Kernal operation for calculating the value of a
 *                     scalar field variable at the quadrature point from
 *                     the basis functions.
 *                     Note for the optimal results, this function should
 *                     be inlined during the compile.
 */
{
  int i;
  *val = *val_dot = *val_old = 0.0;
  if (pd->TimeIntegration == STEADY) {
    for (i = 0; i < dofs; i++) {
      *val += *(base[i]) * phiv[i];
    }
  } else {
    for (i = 0; i < dofs; i++) {
      *val += *(base[i]) * phiv[i];
      *val_dot += *(base_dot[i]) * phiv[i];
      *val_old += *(base_old[i]) * phiv[i];
    }
  }
}

/***************************************************************************/
/****************************************************************************/
/****************************************************************************/

void grad_scalar_fv_fill(double **base, double (*grad_phiv)[DIM], int dofs, double *grad_val)
/*
 */
{
  int i;

  grad_val[0] = grad_val[1] = grad_val[2] = 0.0;
  for (i = 0; i < dofs; i++) {
    grad_val[0] += *(base[i]) * grad_phiv[i][0];
    grad_val[1] += *(base[i]) * grad_phiv[i][1];
    if (VIM == 3)
      grad_val[2] += *(base[i]) * grad_phiv[i][2];
  }
}

void grad_vector_fv_fill(double ***base,
                         double (*grad_phiv)[DIM][DIM][DIM],
                         int dofs,
                         double (*grad_val)[DIM]) {
  int r, i;

  double base_off;
  double(*grad_phiv_off)[DIM];

  memset(grad_val, 0, DIM * DIM * sizeof(double));

  for (r = 0; r < VIM; r++) {
    for (i = 0; i < dofs; i++) {
      base_off = *base[r][i];
      grad_phiv_off = grad_phiv[i][r];

      grad_val[0][0] += base_off * grad_phiv_off[0][0];
      grad_val[1][1] += base_off * grad_phiv_off[1][1];
      grad_val[0][1] += base_off * grad_phiv_off[0][1];
      grad_val[1][0] += base_off * grad_phiv_off[1][0];

      if (VIM == 3) {
        grad_val[2][2] += base_off * grad_phiv_off[2][2];
        grad_val[2][1] += base_off * grad_phiv_off[2][1];
        grad_val[2][0] += base_off * grad_phiv_off[2][0];
        grad_val[1][2] += base_off * grad_phiv_off[1][2];
        grad_val[0][2] += base_off * grad_phiv_off[0][2];
      }

      /*	for ( p=0; p<VIM; p++)
              {
                      for ( q=0; q<VIM; q++)
                      {
                              (grad_val[p][q] += *base[r][i] * grad_phiv[i][r] [p][q];)
                              grad_val[p][q] += base_off * grad_phiv_off[p][q];
                      }
              }*/
    }
  }
}

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

int load_fv(void)

/*******************************************************************************
 * load_fv() -- load up values of all relevant field variables at the
 *              current gauss pt
 *
 * input: (assume the appropriate parts of esp, bf, ei, pd are filled in)
 * ------
 *
 *
 * output: ( this routine fills in the following parts of the fv structure)
 * -------
 *		fv->
 *			T -- temperature	(scalar)
 *			v -- velocity		(vector)
 *			d -- mesh displacement	(vector)
 *			c -- concentration	(multiple scalars)
 *		      por -- porous media       (multiple scalars)
 *			P -- pressure		(scalar)
 *			S -- polymer stress	(tensor)
 *			G -- velocity gradient	(tensor)
 *                     pv -- particle velocity  (vector)
 *                     pG -- particle velocity gradient (tensor)
 *              mp->StateVector[] is filled in as well, for
 *                    pertinent entries that make up the specification of
 *                    the state of the material.
 *
 * NOTE: To accommodate shell elements, this function has been modified
 *       so that fv variables are not zeroed out when they are active
 *       on an element block other than the current one.
 *       The check done for variable v is then:
 *          if ( pd->v[pg->imtrx][v] || upd->vp[pg->imtrx][v] == -1 )
 *       In many cases below, this conditional zeroing is done in
 *       a separate small loop before the main one.
 *
 *
 * Return values:
 *		0 -- if things went OK
 *	       -1 -- if a problem occurred
 *
 * Created:	Fri Mar 18 06:44:41 MST 1994 pasacki@sandia.gov
 *
 * Modified:
 ***************************************************************************/
{
  int v;    /* variable type indicator */
  int i;    /* index */
  int p, q; /* dimension indeces */
  int dim;
  int dofs; /* degrees of freedom for a var in the elem */
  int w;    /* concentration species and porous media vars counter */
  int node, index;
  int N;
  int mode;
  int status = 0;
  int v_s[MAX_MODES][DIM][DIM], v_g[DIM][DIM];
  double rho, *stateVector = mp->StateVector;
  BASIS_FUNCTIONS_STRUCT *bfv;
  int transient_run = pd->TimeIntegration != STEADY;
  int *pdgv = pd->gv;

  status = 0;

  /* load eqn and variable number in tensor form */
  if (pdgv[POLYMER_STRESS11]) {
    status = stress_eqn_pointer(v_s);
    GOMA_EH(status, "stress_eqn_pointer(v_s)");
  }
  if (pdgv[VELOCITY_GRADIENT11]) {
    v_g[0][0] = VELOCITY_GRADIENT11;
    v_g[0][1] = VELOCITY_GRADIENT12;
    v_g[1][0] = VELOCITY_GRADIENT21;
    v_g[1][1] = VELOCITY_GRADIENT22;
    v_g[0][2] = VELOCITY_GRADIENT13;
    v_g[1][2] = VELOCITY_GRADIENT23;
    v_g[2][0] = VELOCITY_GRADIENT31;
    v_g[2][1] = VELOCITY_GRADIENT32;
    v_g[2][2] = VELOCITY_GRADIENT33;
  }

  /*
   * Since it is possible to have a 1D element in a 2D problem,
   * (i.e. shell elements), use the problem dimension instead.
   */
  dim = pd->Num_Dim;

  /*
   * Temperature...
   *    HKM -> Introduced a function to handle the core scalar fill
   *           operation; it should be inlined for maximum speed.
   *           Added the concept that the reference temperature should
   *           be used, when no temperature variable exists in the solution
   *           vector.
   */

  if (pdgv[TEMPERATURE]) {
    v = TEMPERATURE;
    scalar_fv_fill(esp->T, esp_dot->T, esp_old->T, bf[v]->phi, ei[upd->matrix_index[v]]->dof[v],
                   &(fv->T), &(fv_dot->T), &(fv_old->T));
    stateVector[TEMPERATURE] = fv->T;
  } /*else if (upd->vp[pg->imtrx][v] == -1) {
      fv->T     = mp->reference[TEMPERATURE];
      fv_old->T = 0.;
      fv_dot->T = 0.;
      } */

  /*
   * Fill...
   */

  if (pdgv[FILL]) {
    v = FILL;
    scalar_fv_fill(esp->F, esp_dot->F, esp_old->F, bf[v]->phi, ei[upd->matrix_index[v]]->dof[v],
                   &(fv->F), &(fv_dot->F), &(fv_old->F));
    stateVector[v] = fv->F;
    fv_dot_old->F = 0;
    for (int i = 0; i < ei[upd->matrix_index[v]]->dof[v]; i++) {
      if (upd->Total_Num_Matrices > 1) {
        //        fv_dot_old->F += *(pg->matrices[upd->matrix_index[v]].xdot_old -
        //        pg->matrices[upd->matrix_index[v]].xdot +
        //                                  esp_dot->F[i]) * bf[v]->phi[i];
      } else {
        fv_dot_old->F += *(xdot_old_static - xdot_static + esp_dot->F[i]) * bf[v]->phi[i];
      }
    }
  } /*else if (upd->vp[pg->imtrx][v] == -1) {
      fv->F = fv_old->F = fv_dot->F = 0.;
      } */

  /*
   * Suspension Temperature...
   */

  if (pdgv[BOND_EVOLUTION]) {
    v = BOND_EVOLUTION;
    scalar_fv_fill(esp->nn, esp_dot->nn, esp_old->nn, bf[v]->phi, ei[upd->matrix_index[v]]->dof[v],
                   &(fv->nn), &(fv_dot->nn), &(fv_old->nn));
    stateVector[v] = fv->nn;
  } /*else if (upd->vp[pg->imtrx][v] == -1)  {
      fv->nn = fv_old->nn = fv_dot->nn = 0.;
      } */

  /*
   * Voltage...
   */

  if (pdgv[VOLTAGE]) {
    v = VOLTAGE;
    scalar_fv_fill(esp->V, esp_dot->V, esp_old->V, bf[v]->phi, ei[upd->matrix_index[v]]->dof[v],
                   &(fv->V), &(fv_dot->V), &(fv_old->V));
    stateVector[VOLTAGE] = fv->V;
  } /*else if ( upd->vp[pg->imtrx][v] == -1 ) {
      fv->V = fv_old->V = fv_dot->V = 0.;
      }*/

  /*
   * Surface charge density...
   */

  if (pdgv[SURF_CHARGE]) {
    v = SURF_CHARGE;
    scalar_fv_fill(esp->qs, esp_dot->qs, esp_old->qs, bf[v]->phi, ei[upd->matrix_index[v]]->dof[v],
                   &(fv->qs), &(fv_dot->qs), &(fv_old->qs));
    stateVector[SURF_CHARGE] = fv->qs;
  } /*else if ( upd->vp[pg->imtrx][v] == -1 ) {
      fv->qs = fv_old->qs = fv_dot->qs = 0.;
      } */

  /*
   * Structural shell curvature
   */

  if (pdgv[SHELL_CURVATURE]) {
    v = SHELL_CURVATURE;
    scalar_fv_fill(esp->sh_K, esp_dot->sh_K, esp_old->sh_K, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->sh_K), &(fv_dot->sh_K), &(fv_old->sh_K));
    stateVector[SHELL_CURVATURE] = fv->sh_K;
  } /*else if ( upd->vp[pg->imtrx][v] == -1 ) {
      fv->sh_K = fv_old->sh_K = fv_dot->sh_K = 0.;
      } */

  if (pdgv[SHELL_CURVATURE2]) {
    v = SHELL_CURVATURE2;
    scalar_fv_fill(esp->sh_K2, esp_dot->sh_K2, esp_old->sh_K2, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->sh_K2), &(fv_dot->sh_K2),
                   &(fv_old->sh_K2));
    stateVector[SHELL_CURVATURE2] = fv->sh_K2;
  } /*else if ( upd->vp[pg->imtrx][v] == -1 ) {
      fv->sh_K2 = fv_old->sh_K2 = fv_dot->sh_K2 = 0.;
      } */

  /*
   * Structural shell tension
   */

  if (pdgv[SHELL_TENSION]) {
    v = SHELL_TENSION;
    scalar_fv_fill(esp->sh_tens, esp_dot->sh_tens, esp_old->sh_tens, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->sh_tens), &(fv_dot->sh_tens),
                   &(fv_old->sh_tens));
    stateVector[SHELL_TENSION] = fv->sh_tens;
  } /*else if ( updgvp[v] == -1 ) {
      fv->sh_tens = fv_old->sh_tens = fv_dot->sh_tens = 0.;
      }*/
  /*
   * Structural shell x coordinate
   */

  if (pdgv[SHELL_X]) {
    v = SHELL_X;
    scalar_fv_fill(esp->sh_x, esp_dot->sh_x, esp_old->sh_x, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->sh_x), &(fv_dot->sh_x), &(fv_old->sh_x));
    stateVector[SHELL_X] = fv->sh_x;
  } /*else if ( upd->vp[pg->imtrx][v] == -1 ) {
      fv->sh_x = fv_old->sh_x = fv_dot->sh_x    = 0.;
      } */

  /*
   * Structural shell y coordinate
   */

  if (pdgv[SHELL_Y]) {
    v = SHELL_Y;
    scalar_fv_fill(esp->sh_y, esp_dot->sh_y, esp_old->sh_y, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->sh_y), &(fv_dot->sh_y), &(fv_old->sh_y));
    stateVector[SHELL_Y] = fv->sh_y;
  } /* else if ( upd->vp[pg->imtrx][v] == -1 ) {
       fv->sh_y = fv_old->sh_y = fv_dot->sh_y    = 0.;
       } */

  /*
   * Shell user
   */

  if (pdgv[SHELL_USER]) {
    v = SHELL_USER;
    scalar_fv_fill(esp->sh_u, esp_dot->sh_u, esp_old->sh_u, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->sh_u), &(fv_dot->sh_u), &(fv_old->sh_u));
    stateVector[SHELL_USER] = fv->sh_u;
  } /* else if ( upd->vp[pg->imtrx][v] == -1 ) {
       fv->sh_u = fv_old->sh_u = fv_dot->sh_u    = 0.;
       } */

  /*
   * Shear rate from second invariant of rate of strain tensor
   */

  if (pdgv[SHEAR_RATE]) {
    v = SHEAR_RATE;
    fv->SH = 0.0;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (i = 0; i < dofs; i++) {
      fv->SH += *esp->SH[i] * bf[v]->phi[i];
    }
  }
  /*  else if ( upd->vp[pg->imtrx][v] == -1)
      {
      fv->SH = 0.0;
      } */

  /* Square of the norm of the potential field, |E|^2 */

  if (pdgv[ENORM]) {
    v = ENORM;
    fv->Enorm = 0.0;
    fv_old->Enorm = 0.0;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (i = 0; i < dofs; i++) {
      fv->Enorm += *esp->Enorm[i] * bf[v]->phi[i];
      fv_old->Enorm += *esp_old->Enorm[i] * bf[v]->phi[i];
    }
  } /*  else if ( upd->vp[pg->imtrx][v] == -1) {
        fv->Enorm = 0.0;
        fv_old->Enorm = 0.0;
        } */

  /*
   * Curvature of level set function
   */

  if (pdgv[CURVATURE]) {
    v = CURVATURE;
    scalar_fv_fill(esp->H, esp_dot->H, esp_old->H, bf[v]->phi, ei[upd->matrix_index[v]]->dof[v],
                   &(fv->H), &(fv_dot->H), &(fv_old->H));
  }
  /*  else if ( upd->vp[pg->imtrx][v] == -1)
      {
      fv->H = fv_dot->H = fv_old->H = 0.0;
      }
  */

  /*
   *  Normal to level set function
   */

  if (pdgv[NORMAL1]) {
    v = NORMAL1;
    scalar_fv_fill(esp->n[0], esp_dot->n[0], esp_old->n[0], bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->n[0]), &(fv_dot->n[0]), &(fv_old->n[0]));
  }
  /*  else if ( upd->vp[pg->imtrx][v] == -1)
      {
      fv->n[0] = fv_dot->n[0] = fv_old->n[0] = 0.0;
      } */

  if (pdgv[NORMAL2]) {
    v = NORMAL2;
    scalar_fv_fill(esp->n[1], esp_dot->n[1], esp_old->n[1], bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->n[1]), &(fv_dot->n[1]), &(fv_old->n[1]));
  }
  /* else if ( upd->vp[pg->imtrx][v] == -1)
     {
     fv->n[1] = fv_dot->n[1] = fv_old->n[1] = 0.0;
     }*/

  if (pdgv[NORMAL3]) {
    v = NORMAL3;
    scalar_fv_fill(esp->n[2], esp_dot->n[2], esp_old->n[2], bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->n[2]), &(fv_dot->n[2]), &(fv_old->n[2]));
  }
  /*  else if ( upd->vp[pg->imtrx][v] == -1)
      {
      fv->n[2] = fv_dot->n[2] = fv_old->n[2] = 0.0;
      } */

  /*
   *  shell element orientation angles
   */

  if (pdgv[SHELL_ANGLE1]) {
    v = SHELL_ANGLE1;
    fv->sh_ang[0] = 0.;
    for (i = 0; i < ei[upd->matrix_index[v]]->dof[v]; i++) {
      if ((*esp->sh_ang[0][i] - *esp->sh_ang[0][0]) > M_PIE)
        fv->sh_ang[0] += bf[v]->phi[i] * (*esp->sh_ang[0][i] - 2. * M_PIE);
      else if ((*esp->sh_ang[0][i] - *esp->sh_ang[0][0]) < -M_PIE)
        fv->sh_ang[0] += bf[v]->phi[i] * (*esp->sh_ang[0][i] + 2. * M_PIE);
      else
        fv->sh_ang[0] += bf[v]->phi[i] * *esp->sh_ang[0][i];
    }
    /* surely no one will be using these */
    fv_dot->sh_ang[0] = fv_old->sh_ang[0] = 0.0;
  } /*
      else if ( upd->vp[pg->imtrx][v] == -1)
      {
      fv->sh_ang[0] = fv_dot->sh_ang[0] = fv_old->sh_ang[0] = 0.0;
      }*/

  if (pdgv[SHELL_ANGLE2]) {
    v = SHELL_ANGLE2;
    fv->sh_ang[1] = 0.;
    for (i = 0; i < ei[upd->matrix_index[v]]->dof[v]; i++) {
      if ((*esp->sh_ang[1][i] - *esp->sh_ang[1][0]) > M_PIE)
        fv->sh_ang[1] += bf[v]->phi[i] * (*esp->sh_ang[1][i] - 2. * M_PIE);
      else if ((*esp->sh_ang[1][i] - *esp->sh_ang[1][0]) < -M_PIE)
        fv->sh_ang[1] += bf[v]->phi[i] * (*esp->sh_ang[1][i] + 2. * M_PIE);
      else
        fv->sh_ang[1] += bf[v]->phi[i] * *esp->sh_ang[1][i];
    }
    /* surely no one will be using these */
    fv_dot->sh_ang[1] = fv_old->sh_ang[1] = 0.0;
  }
  /*  else if ( upd->vp[pg->imtrx][v] == -1)
      {
      fv->sh_ang[1] = fv_dot->sh_ang[1] = fv_old->sh_ang[1] = 0.0;
      }*/

  /*
   * Surface Rheo shell piece
   */

  if (pdgv[SHELL_SURF_DIV_V]) {
    v = SHELL_SURF_DIV_V;
    scalar_fv_fill(esp->div_s_v, esp_dot->div_s_v, esp_old->div_s_v, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->div_s_v), &(fv_dot->div_s_v),
                   &(fv_old->div_s_v));
    stateVector[SHELL_SURF_DIV_V] = fv->div_s_v;
  } /*else if ( upd->vp[pg->imtrx][v] == -1 ) {
      fv->div_s_v = fv_old->div_s_v = fv_dot->div_s_v    = 0.;
      }*/

  if (pdgv[SHELL_SURF_CURV]) {
    v = SHELL_SURF_CURV;
    scalar_fv_fill(esp->curv, esp_dot->curv, esp_old->curv, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->curv), &(fv_dot->curv), &(fv_old->curv));
    stateVector[SHELL_SURF_CURV] = fv->curv;
  } /*else if ( upd->vp[pg->imtrx][v] == -1 ) {
      fv->curv = fv_old->curv = fv_dot->curv    = 0.;
      }*/

  if (pdgv[N_DOT_CURL_V]) {
    v = N_DOT_CURL_V;
    scalar_fv_fill(esp->n_dot_curl_s_v, esp_dot->n_dot_curl_s_v, esp_old->n_dot_curl_s_v,
                   bf[v]->phi, ei[upd->matrix_index[v]]->dof[v], &(fv->n_dot_curl_s_v),
                   &(fv_dot->n_dot_curl_s_v), &(fv_old->n_dot_curl_s_v));
    stateVector[N_DOT_CURL_V] = fv->n_dot_curl_s_v;
  } /*else if ( upd->vp[pg->imtrx][v] == -1 ) {
      fv->n_dot_curl_s_v = fv_old->n_dot_curl_s_v = fv_dot->n_dot_curl_s_v    = 0.;
      }*/

  if (pdgv[GRAD_S_V_DOT_N1]) {
    v = GRAD_S_V_DOT_N1;
    scalar_fv_fill(esp->grad_v_dot_n[0], esp_dot->grad_v_dot_n[0], esp_old->grad_v_dot_n[0],
                   bf[v]->phi, ei[upd->matrix_index[v]]->dof[v], &(fv->grad_v_dot_n[0]),
                   &(fv_dot->grad_v_dot_n[0]), &(fv_old->grad_v_dot_n[0]));
    stateVector[GRAD_S_V_DOT_N1] = fv->grad_v_dot_n[0];
  } /* else if ( upd->vp[pg->imtrx][v] == -1 ) {
       fv->grad_v_dot_n[0] = fv_old->grad_v_dot_n[0] = fv_dot->grad_v_dot_n[0]    = 0.;
       }*/

  if (pdgv[GRAD_S_V_DOT_N2]) {
    v = GRAD_S_V_DOT_N2;
    scalar_fv_fill(esp->grad_v_dot_n[1], esp_dot->grad_v_dot_n[1], esp_old->grad_v_dot_n[1],
                   bf[v]->phi, ei[upd->matrix_index[v]]->dof[v], &(fv->grad_v_dot_n[1]),
                   &(fv_dot->grad_v_dot_n[1]), &(fv_old->grad_v_dot_n[1]));
    stateVector[GRAD_S_V_DOT_N2] = fv->grad_v_dot_n[1];
  } /* else if ( upd->vp[pg->imtrx][v] == -1 ) {
       fv->grad_v_dot_n[1] = fv_old->grad_v_dot_n[1] = fv_dot->grad_v_dot_n[1]    = 0.;


       } */

  if (pdgv[GRAD_S_V_DOT_N3]) {
    v = GRAD_S_V_DOT_N3;
    scalar_fv_fill(esp->grad_v_dot_n[2], esp_dot->grad_v_dot_n[2], esp_old->grad_v_dot_n[2],
                   bf[v]->phi, ei[upd->matrix_index[v]]->dof[v], &(fv->grad_v_dot_n[2]),
                   &(fv_dot->grad_v_dot_n[2]), &(fv_old->grad_v_dot_n[2]));
    stateVector[GRAD_S_V_DOT_N3] = fv->grad_v_dot_n[2];
  } /*else if ( upd->vp[pg->imtrx][v] == -1 ) {
      fv->grad_v_dot_n[2] = fv_old->grad_v_dot_n[2] = fv_dot->grad_v_dot_n[2]    = 0.;
      }*/

  if (pdgv[SHELL_DIFF_FLUX]) {
    v = SHELL_DIFF_FLUX;
    scalar_fv_fill(esp->sh_J, esp_dot->sh_J, esp_old->sh_J, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->sh_J), &(fv_dot->sh_J), &(fv_old->sh_J));
    stateVector[SHELL_DIFF_FLUX] = fv->sh_J;
  } /*else if ( upd->vp[pg->imtrx][v] == -1 ) {
      fv->sh_J = fv_old->sh_J = fv_dot->sh_J = 0.;
      }*/

  if (pdgv[SHELL_DIFF_CURVATURE]) {
    v = SHELL_DIFF_CURVATURE;
    scalar_fv_fill(esp->sh_Kd, esp_dot->sh_Kd, esp_old->sh_Kd, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->sh_Kd), &(fv_dot->sh_Kd),
                   &(fv_old->sh_Kd));
    stateVector[SHELL_DIFF_CURVATURE] = fv->sh_Kd;
  } /* else if ( upd->vp[pg->imtrx][v] == -1 ) {
       fv->sh_Kd = fv_old->sh_Kd = fv_dot->sh_Kd = 0.;
       }*/

  if (pdgv[SHELL_NORMAL1]) {
    v = SHELL_NORMAL1;
    scalar_fv_fill(esp->n[0], esp_dot->n[0], esp_old->n[0], bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->n[0]), &(fv_dot->n[0]), &(fv_old->n[0]));
    stateVector[SHELL_NORMAL1] = fv->n[0];
  } /*else if ( upd->vp[pg->imtrx][v] == -1 ) {
      fv->n[0] = fv_old->n[0] = fv_dot->n[0] = 0.;
      }*/

  if (pdgv[SHELL_NORMAL2]) {
    v = SHELL_NORMAL2;
    scalar_fv_fill(esp->n[1], esp_dot->n[1], esp_old->n[1], bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->n[1]), &(fv_dot->n[1]), &(fv_old->n[1]));
    stateVector[SHELL_NORMAL2] = fv->n[1];
  } /*else if ( upd->vp[pg->imtrx][v] == -1 ) {
      fv->n[1] = fv_old->n[1] = fv_dot->n[1] = 0.;
      }*/

  if (pdgv[SHELL_NORMAL3]) {
    v = SHELL_NORMAL3;
    scalar_fv_fill(esp->n[2], esp_dot->n[2], esp_old->n[2], bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->n[2]), &(fv_dot->n[2]), &(fv_old->n[2]));
    stateVector[SHELL_NORMAL3] = fv->n[2];
  } /*else if ( upd->vp[pg->imtrx][v] == -1 ) {
      fv->n[2] = fv_old->n[2] = fv_dot->n[2] = 0.;
      }*/

  if (pdgv[SHELL_NORMAL1] && pdgv[SHELL_NORMAL2] && pdgv[SHELL_NORMAL3]) {
    memset(fv->d_n_dxi, 0.0, sizeof(double) * DIM * DIM);
    for (i = 0; i < ei[upd->matrix_index[v]]->dof[SHELL_NORMAL1]; i++) {
      fv->d_n_dxi[0][0] += *esp->n[0][i] * bf[SHELL_NORMAL1]->dphidxi[i][0];
      fv->d_n_dxi[1][0] += *esp->n[1][i] * bf[SHELL_NORMAL2]->dphidxi[i][0];
      fv->d_n_dxi[2][0] += *esp->n[2][i] * bf[SHELL_NORMAL3]->dphidxi[i][0];

      fv->d_n_dxi[0][1] += *esp->n[0][i] * bf[SHELL_NORMAL1]->dphidxi[i][1];
      fv->d_n_dxi[1][1] += *esp->n[1][i] * bf[SHELL_NORMAL2]->dphidxi[i][1];
      fv->d_n_dxi[2][1] += *esp->n[2][i] * bf[SHELL_NORMAL3]->dphidxi[i][1];
    }
  }

  /*
   *	Acoustic Pressure
   */

  if (pdgv[ACOUS_PREAL]) {
    v = ACOUS_PREAL;
    scalar_fv_fill(esp->apr, esp_dot->apr, esp_old->apr, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->apr), &(fv_dot->apr), &(fv_old->apr));
    stateVector[ACOUS_PREAL] = fv->apr;
  } /*else if (upd->vp[pg->imtrx][v] == -1) {
      fv->apr = fv_old->apr = fv_dot->apr = 0.;
      }*/

  if (pdgv[ACOUS_PIMAG]) {
    v = ACOUS_PIMAG;
    scalar_fv_fill(esp->api, esp_dot->api, esp_old->api, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->api), &(fv_dot->api), &(fv_old->api));
    stateVector[ACOUS_PIMAG] = fv->api;
  } /*else if (upd->vp[pg->imtrx][v] == -1) {
      fv->api = fv_old->api = fv_dot->api = 0.;
      }*/

  if (pdgv[ACOUS_REYN_STRESS]) {
    v = ACOUS_REYN_STRESS;
    scalar_fv_fill(esp->ars, esp_dot->ars, esp_old->ars, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->ars), &(fv_dot->ars), &(fv_old->ars));
    stateVector[ACOUS_REYN_STRESS] = fv->ars;
  } /*else if (upd->vp[pg->imtrx][v] == -1) {
      fv->ars = fv_old->ars = fv_dot->ars = 0.;
      }*/

  if (pdgv[SHELL_BDYVELO]) {
    v = SHELL_BDYVELO;
    scalar_fv_fill(esp->sh_bv, esp_dot->sh_bv, esp_old->sh_bv, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->sh_bv), &(fv_dot->sh_bv),
                   &(fv_old->sh_bv));
    stateVector[SHELL_BDYVELO] = fv->sh_bv;
  }

  if (pdgv[SHELL_LUBP]) {
    v = SHELL_LUBP;
    scalar_fv_fill(esp->sh_p, esp_dot->sh_p, esp_old->sh_p, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->sh_p), &(fv_dot->sh_p), &(fv_old->sh_p));
    stateVector[SHELL_LUBP] = fv->sh_p;
  }

  if (pdgv[LUBP]) {
    v = LUBP;
    scalar_fv_fill(esp->lubp, esp_dot->lubp, esp_old->lubp, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->lubp), &(fv_dot->lubp), &(fv_old->lubp));
    stateVector[LUBP] = fv->lubp;
  }

  if (pdgv[LUBP_2]) {
    v = LUBP_2;
    scalar_fv_fill(esp->lubp_2, esp_dot->lubp_2, esp_old->lubp_2, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->lubp_2), &(fv_dot->lubp_2),
                   &(fv_old->lubp_2));
    stateVector[LUBP_2] = fv->lubp_2;
  }

  if (pdgv[SHELL_FILMP]) {
    v = SHELL_FILMP;
    scalar_fv_fill(esp->sh_fp, esp_dot->sh_fp, esp_old->sh_fp, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->sh_fp), &(fv_dot->sh_fp),
                   &(fv_old->sh_fp));
    stateVector[SHELL_FILMP] = fv->sh_fp;
  }

  if (pdgv[SHELL_FILMH]) {
    v = SHELL_FILMH;
    scalar_fv_fill(esp->sh_fh, esp_dot->sh_fh, esp_old->sh_fh, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->sh_fh), &(fv_dot->sh_fh),
                   &(fv_old->sh_fh));
    stateVector[SHELL_FILMH] = fv->sh_fh;
  }

  if (pdgv[SHELL_PARTC]) {
    v = SHELL_PARTC;
    scalar_fv_fill(esp->sh_pc, esp_dot->sh_pc, esp_old->sh_pc, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->sh_pc), &(fv_dot->sh_pc),
                   &(fv_old->sh_pc));
    stateVector[SHELL_PARTC] = fv->sh_pc;
  }

  if (pdgv[SHELL_SAT_CLOSED]) {
    v = SHELL_SAT_CLOSED;
    scalar_fv_fill(esp->sh_sat_closed, esp_dot->sh_sat_closed, esp_old->sh_sat_closed, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->sh_sat_closed), &(fv_dot->sh_sat_closed),
                   &(fv_old->sh_sat_closed));
    stateVector[SHELL_SAT_CLOSED] = fv->sh_sat_closed;
  }
  if (pdgv[SHELL_PRESS_OPEN]) {
    v = SHELL_PRESS_OPEN;
    scalar_fv_fill(esp->sh_p_open, esp_dot->sh_p_open, esp_old->sh_p_open, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->sh_p_open), &(fv_dot->sh_p_open),
                   &(fv_old->sh_p_open));
    stateVector[SHELL_PRESS_OPEN] = fv->sh_p_open;
  }
  if (pdgv[SHELL_PRESS_OPEN_2]) {
    v = SHELL_PRESS_OPEN_2;
    scalar_fv_fill(esp->sh_p_open_2, esp_dot->sh_p_open_2, esp_old->sh_p_open_2, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->sh_p_open_2), &(fv_dot->sh_p_open_2),
                   &(fv_old->sh_p_open_2));
    stateVector[SHELL_PRESS_OPEN_2] = fv->sh_p_open_2;
  }
  if (pdgv[SHELL_SAT_1]) {
    v = SHELL_SAT_1;
    scalar_fv_fill(esp->sh_sat_1, esp_dot->sh_sat_1, esp_old->sh_sat_1, bf[v]->phi,
                   ei[pg->imtrx]->dof[v], &(fv->sh_sat_1), &(fv_dot->sh_sat_1),
                   &(fv_old->sh_sat_1));
    stateVector[SHELL_SAT_1] = fv->sh_sat_1;
  }
  if (pdgv[SHELL_SAT_2]) {
    v = SHELL_SAT_2;
    scalar_fv_fill(esp->sh_sat_2, esp_dot->sh_sat_2, esp_old->sh_sat_2, bf[v]->phi,
                   ei[pg->imtrx]->dof[v], &(fv->sh_sat_2), &(fv_dot->sh_sat_2),
                   &(fv_old->sh_sat_2));
    stateVector[SHELL_SAT_2] = fv->sh_sat_2;
  }
  if (pdgv[SHELL_SAT_3]) {
    v = SHELL_SAT_3;
    scalar_fv_fill(esp->sh_sat_3, esp_dot->sh_sat_3, esp_old->sh_sat_3, bf[v]->phi,
                   ei[pg->imtrx]->dof[v], &(fv->sh_sat_3), &(fv_dot->sh_sat_3),
                   &(fv_old->sh_sat_3));
    stateVector[SHELL_SAT_3] = fv->sh_sat_3;
  }
  if (pdgv[SHELL_TEMPERATURE]) {
    v = SHELL_TEMPERATURE;
    scalar_fv_fill(esp->sh_t, esp_dot->sh_t, esp_old->sh_t, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->sh_t), &(fv_dot->sh_t), &(fv_old->sh_t));
    stateVector[SHELL_TEMPERATURE] = fv->sh_t;
  }
  if (pdgv[SHELL_DELTAH]) {
    v = SHELL_DELTAH;
    scalar_fv_fill(esp->sh_dh, esp_dot->sh_dh, esp_old->sh_dh, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->sh_dh), &(fv_dot->sh_dh),
                   &(fv_old->sh_dh));
    stateVector[SHELL_DELTAH] = fv->sh_dh;
  }
  if (pdgv[SHELL_LUB_CURV]) {
    v = SHELL_LUB_CURV;
    scalar_fv_fill(esp->sh_l_curv, esp_dot->sh_l_curv, esp_old->sh_l_curv, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->sh_l_curv), &(fv_dot->sh_l_curv),
                   &(fv_old->sh_l_curv));
    stateVector[SHELL_LUB_CURV] = fv->sh_l_curv;
  }
  if (pdgv[SHELL_LUB_CURV_2]) {
    v = SHELL_LUB_CURV_2;
    scalar_fv_fill(esp->sh_l_curv_2, esp_dot->sh_l_curv_2, esp_old->sh_l_curv_2, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->sh_l_curv_2), &(fv_dot->sh_l_curv_2),
                   &(fv_old->sh_l_curv_2));
    stateVector[SHELL_LUB_CURV] = fv->sh_l_curv;
  }
  if (pdgv[SHELL_SAT_GASN]) {
    v = SHELL_SAT_GASN;
    scalar_fv_fill(esp->sh_sat_gasn, esp_dot->sh_sat_gasn, esp_old->sh_sat_gasn, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->sh_sat_gasn), &(fv_dot->sh_sat_gasn),
                   &(fv_old->sh_sat_gasn));
    stateVector[SHELL_SAT_GASN] = fv->sh_sat_gasn;
  }
  if (pdgv[SHELL_SHEAR_TOP]) {
    v = SHELL_SHEAR_TOP;
    scalar_fv_fill(esp->sh_shear_top, esp_dot->sh_shear_top, esp_old->sh_shear_top, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->sh_shear_top), &(fv_dot->sh_shear_top),
                   &(fv_old->sh_shear_top));
    stateVector[SHELL_SHEAR_TOP] = fv->sh_shear_top;
  }
  if (pdgv[SHELL_SHEAR_BOT]) {
    v = SHELL_SHEAR_BOT;
    scalar_fv_fill(esp->sh_shear_bot, esp_dot->sh_shear_bot, esp_old->sh_shear_bot, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->sh_shear_bot), &(fv_dot->sh_shear_bot),
                   &(fv_old->sh_shear_bot));
    stateVector[SHELL_SHEAR_BOT] = fv->sh_shear_bot;
  }
  if (pdgv[SHELL_CROSS_SHEAR]) {
    v = SHELL_CROSS_SHEAR;
    scalar_fv_fill(esp->sh_cross_shear, esp_dot->sh_cross_shear, esp_old->sh_cross_shear,
                   bf[v]->phi, ei[upd->matrix_index[v]]->dof[v], &(fv->sh_cross_shear),
                   &(fv_dot->sh_cross_shear), &(fv_old->sh_cross_shear));
    stateVector[SHELL_CROSS_SHEAR] = fv->sh_cross_shear;
  }
  if (pdgv[MAX_STRAIN]) {
    v = MAX_STRAIN;
    scalar_fv_fill(esp->max_strain, esp_dot->max_strain, esp_old->max_strain, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->max_strain), &(fv_dot->max_strain),
                   &(fv_old->max_strain));
    stateVector[MAX_STRAIN] = fv->max_strain;
  }
  if (pdgv[CUR_STRAIN]) {
    v = CUR_STRAIN;
    scalar_fv_fill(esp->cur_strain, esp_dot->cur_strain, esp_old->cur_strain, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->cur_strain), &(fv_dot->cur_strain),
                   &(fv_old->cur_strain));
    stateVector[CUR_STRAIN] = fv->cur_strain;
  }

  if (pdgv[EDDY_NU]) {
    v = EDDY_NU;
    scalar_fv_fill(esp->eddy_nu, esp_dot->eddy_nu, esp_old->eddy_nu, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->eddy_nu), &(fv_dot->eddy_nu),
                   &(fv_old->eddy_nu));
    stateVector[EDDY_NU] = fv->eddy_nu;
  }

  if (pdgv[LIGHT_INTP]) {
    v = LIGHT_INTP;
    scalar_fv_fill(esp->poynt[0], esp_dot->poynt[0], esp_old->poynt[0], bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->poynt[0]), &(fv_dot->poynt[0]),
                   &(fv_old->poynt[0]));
    stateVector[LIGHT_INTP] = fv->poynt[0];
  }
  if (pdgv[LIGHT_INTM]) {
    v = LIGHT_INTM;
    scalar_fv_fill(esp->poynt[1], esp_dot->poynt[1], esp_old->poynt[1], bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->poynt[1]), &(fv_dot->poynt[1]),
                   &(fv_old->poynt[1]));
    stateVector[LIGHT_INTM] = fv->poynt[1];
  }
  if (pdgv[LIGHT_INTD]) {
    v = LIGHT_INTD;
    scalar_fv_fill(esp->poynt[2], esp_dot->poynt[2], esp_old->poynt[2], bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->poynt[2]), &(fv_dot->poynt[2]),
                   &(fv_old->poynt[2]));
    stateVector[LIGHT_INTD] = fv->poynt[2];
  }
  if (pdgv[RESTIME]) {
    v = RESTIME;
    scalar_fv_fill(esp->restime, esp_dot->restime, esp_old->restime, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->restime), &(fv_dot->restime),
                   &(fv_old->restime));
    stateVector[RESTIME] = fv->restime;
  }
  /*
   *	Porous sink mass
   */

  if (pdgv[POR_SINK_MASS]) {
    v = POR_SINK_MASS;
    scalar_fv_fill(esp->sink_mass, esp_dot->sink_mass, esp_old->sink_mass, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->sink_mass), &(fv_dot->sink_mass),
                   &(fv_old->sink_mass));
    stateVector[POR_SINK_MASS] = fv->sink_mass;
  } /*else if (upd->vp[pg->imtrx][v] == -1) {
      fv->sink_mass = fv_old->sink_mass = fv_dot->sink_mass = 0.;
      }*/

  /*
   * Pressure...
   */

  if (pdgv[PRESSURE]) {
    v = PRESSURE;
    scalar_fv_fill(esp->P, esp_dot->P, esp_old->P, bf[v]->phi, ei[upd->matrix_index[v]]->dof[v],
                   &(fv->P), &(fv_dot->P), &(fv_old->P));
    stateVector[v] = fv->P + upd->Pressure_Datum;
  }

  if (pdgv[PSTAR]) {
    v = PSTAR;
    scalar_fv_fill(esp->P_star, esp_dot->P_star, esp_old->P_star, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->P_star), &(fv_dot->P_star),
                   &(fv_old->P_star));
    stateVector[v] = fv->P_star + upd->Pressure_Datum;
  }

  if (pdgv[EM_CONT_REAL]) {
    v = EM_CONT_REAL;
    scalar_fv_fill(esp->epr, esp_dot->epr, esp_old->epr, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->epr), &(fv_dot->epr), &(fv_old->epr));
    stateVector[v] = fv->epr;
  }
  if (pdgv[EM_CONT_IMAG]) {
    v = EM_CONT_IMAG;
    scalar_fv_fill(esp->epi, esp_dot->epi, esp_old->epi, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->epi), &(fv_dot->epi), &(fv_old->epi));
    stateVector[v] = fv->epi;
  }
  /*
   *  Set the state vector pressure datum to include an additional
   *  uniform pressure datum.
   */

  /*
   * Mesh displacement (vector)...
   * and positions (vector)
   */

  /*
   * Insure default value of 3rd coordinate is zero for low
   * dimensional problems...DIM=3, while dim can be less...
   */

  for (p = 0; p < dim; p++) {
    fv->x0[p] = 0.0;
    fv->x[p] = 0.0;
    fv_old->x[p] = 0.0;
    fv_dot->x[p] = 0.0;
    fv_dot_old->x[p] = 0.0;

    fv->d[p] = 0.0;
    fv_old->d[p] = 0.0;
    fv_dot->d[p] = 0.0;
    fv_dot_old->d[p] = 0.0;

    if (tran->solid_inertia) {
      fv_dot_dot->x[p] = 0;
      fv_dot_dot->d[p] = 0;
    }

    /*
     * If this is a shell element, mesh displacements may not be
     * defined on this element block even if the mesh is deforming.
     * In this case, displacements are defined on a neighbor block
     * and ei[upd->matrix_index[v]]->deforming_mesh will be TRUE, so that the true
     * displaced coordinates can be loaded here.
     */
    if (ei[pg->imtrx]->deforming_mesh) {
      /*
       * ShapeVar will always be mesh displacement where it is defined.
       * Otherwise (e.g. for shell elements), it will be set correctly.
       */
      bfv = bf[pd->ShapeVar];
      v = pd->ShapeVar;
      dofs = ei[upd->matrix_index[v]]->dof[v];

      for (i = 0; i < dofs; i++) {
        node = ei[upd->matrix_index[v]]->dof_list[R_MESH1][i];
        index = Proc_Elem_Connect[Proc_Connect_Ptr[ei[upd->matrix_index[v]]->ielem] + node];
        fv->d[p] += *esp->d[p][i] * bfv->phi[i];
        fv->x[p] += (Coor[p][index] + *esp->d[p][i]) * bfv->phi[i];
        fv->x0[p] += Coor[p][index] * bfv->phi[i];

        /*
         * HKM -> We calculate the old xdot's by doing a little pointer
         *        arithmetic. Basically we use esp_dot as an offset to
         *        the base and then correct for the base address.
         */
        if (transient_run) {
          fv_old->d[p] += *esp_old->d[p][i] * bfv->phi[i];
          fv_old->x[p] += (Coor[p][index] + *esp_old->d[p][i]) * bf[v]->phi[i];
          fv_dot->d[p] += *esp_dot->d[p][i] * bfv->phi[i];
          fv_dot->x[p] += *esp_dot->d[p][i] * bfv->phi[i];
          if (tran->solid_inertia) {
            fv_dot_dot->d[p] += *esp_dbl_dot->d[p][i] * bfv->phi[i];
            fv_dot_dot->x[p] += *esp_dbl_dot->d[p][i] * bfv->phi[i];
          }

          if (upd->Total_Num_Matrices > 1) {
            fv_dot_old->x[p] += *(pg->matrices[upd->matrix_index[v]].xdot_old -
                                  pg->matrices[upd->matrix_index[v]].xdot + esp_dot->d[p][i]) *
                                bfv->phi[i];
            fv_dot_old->d[p] += *(pg->matrices[upd->matrix_index[v]].xdot_old -
                                  pg->matrices[upd->matrix_index[v]].xdot + esp_dot->d[p][i]) *
                                bfv->phi[i];
          } else {
            fv_dot_old->x[p] += *(xdot_old_static - xdot_static + esp_dot->d[p][i]) * bfv->phi[i];
            fv_dot_old->d[p] += *(xdot_old_static - xdot_static + esp_dot->d[p][i]) * bfv->phi[i];
          }
          if (tran->solid_inertia) {
            fv_dot_dot_old->d[p] +=
                *(x_dbl_dot_old_static - x_dbl_dot_static + esp_dbl_dot->d[p][i]) * bfv->phi[i];
          }
        } else {
          fv_old->d[p] = fv->d[p];
          fv_old->x[p] = fv->x[p]; /* Fixed grid stays fixed thru time. */
        }
      }
    }

    else /* Zero these only if not using mesh displacements: */
    {
      v = pd->ShapeVar;
      dofs = ei[upd->matrix_index[v]]->dof[v];

      if (bf[v]->interpolation == I_N1) {
        dofs = bf[v]->shape_dof;
      }

      fv->d[p] = 0.0;
      fv_old->d[p] = 0.0;
      fv_dot->d[p] = 0.0;
      fv_dot_dot->d[p] = 0.0;

      fv->x[p] = 0.;
      for (i = 0; i < dofs; i++) {
        if (bf[v]->interpolation == I_N1) {
          node = i;
        } else {
          node = ei[upd->matrix_index[v]]->dof_list[v][i];
        }
        index = Proc_Elem_Connect[Proc_Connect_Ptr[ei[upd->matrix_index[v]]->ielem] + node];
        fv->x[p] += (Coor[p][index]) * bf[v]->phi[i];
      }
      fv_old->x[p] = fv->x[p]; /* Fixed grid stays fixed thru time. */
      fv->x0[p] = fv->x[p];
      fv_dot->x[p] = 0.0;
      fv_dot_dot->x[p] = 0.0;
    }
  }

  /*
   * SOLID displacement (vector)...
   * and positions (vector)
   */
  for (p = 0; pdgv[SOLID_DISPLACEMENT1] && p < dim; p++) {
    v = SOLID_DISPLACEMENT1 + p;
    fv->d_rs[p] = 0.;
    fv_old->d_rs[p] = 0.;
    fv_dot->d_rs[p] = 0.;
    fv_dot_dot->d_rs[p] = 0.;
    fv_dot_old->d_rs[p] = 0.;
    fv_dot_dot_old->d_rs[p] = 0.;

    if (pdgv[v]) {
      dofs = ei[upd->matrix_index[v]]->dof[v];

      for (i = 0; i < dofs; i++) {
        node = ei[upd->matrix_index[v]]->dof_list[R_SOLID1][i];
        index = Proc_Elem_Connect[Proc_Connect_Ptr[ei[upd->matrix_index[v]]->ielem] + node];
        fv->d_rs[p] += *esp->d_rs[p][i] * bf[v]->phi[i];

        if (pd->TimeIntegration != STEADY) {
          fv_old->d_rs[p] += *esp_old->d_rs[p][i] * bf[v]->phi[i];
          fv_dot->d_rs[p] += *esp_dot->d_rs[p][i] * bf[v]->phi[i];
          if (tran->solid_inertia) {
            fv_dot_dot->d_rs[p] += *esp_dbl_dot->d_rs[p][i] * bf[v]->phi[i];
          }

          if (upd->Total_Num_Matrices > 1) {
            fv_dot_old->d_rs[p] +=
                *(pg->matrices[upd->matrix_index[v]].xdot_old -
                  pg->matrices[upd->matrix_index[v]].xdot + esp_dot->d_rs[p][i]) *
                bf[v]->phi[i];
          } else {
            fv_dot_old->d_rs[p] +=
                *(xdot_old_static - xdot_static + esp_dot->d_rs[p][i]) * bf[v]->phi[i];
          }

          if (tran->solid_inertia) {
            fv_dot_dot_old->d_rs[p] +=
                *(x_dbl_dot_old_static - x_dbl_dot_static + esp_dbl_dot->d_rs[p][i]) *
                bf[v]->phi[i];
          }
        } else {
          fv_old->d_rs[p] = fv->d_rs[p];
        }
      }
    }
  }

  status = load_coordinate_scales(pd->CoordinateSystem, fv);
  GOMA_EH(status, "load_coordinate_scales(fv)");

  /*
   * Velocity (vector)...
   */

  /*
   * Default: all velocities are zero...
   */
  for (p = 0; p < WIM; p++) {
    v = VELOCITY1 + p;
    if (pdgv[v]) {
      dofs = ei[upd->matrix_index[v]]->dof[v];
      fv->v[p] = 0.;
      fv_old->v[p] = 0.;
      fv_dot->v[p] = 0.;
      for (i = 0; i < dofs; i++) {
        fv->v[p] += *esp->v[p][i] * bf[v]->phi[i];
        if (pd->TimeIntegration != STEADY) {
          fv_old->v[p] += *esp_old->v[p][i] * bf[v]->phi[i];
          fv_dot->v[p] += *esp_dot->v[p][i] * bf[v]->phi[i];
        }
      }
    }
    stateVector[VELOCITY1 + p] = fv->v[p];
  }

  for (p = 0; p < WIM; p++) {
    v = USTAR + p;
    if (pdgv[v]) {
      dofs = ei[upd->matrix_index[v]]->dof[v];
      fv->v_star[p] = 0.;
      fv_old->v_star[p] = 0.;
      fv_dot->v_star[p] = 0.;
      for (i = 0; i < dofs; i++) {
        fv->v_star[p] += *esp->v_star[p][i] * bf[v]->phi[i];
        if (pd->TimeIntegration != STEADY) {
          fv_old->v_star[p] += *esp_old->v_star[p][i] * bf[v]->phi[i];
          fv_dot->v_star[p] += *esp_dot->v_star[p][i] * bf[v]->phi[i];
        }
      }
    }
    stateVector[VELOCITY1 + p] = fv->v_star[p];
  }

  /*
   * Particle velocity (vector)...
   */

  /*
   * Default: all velocities are zero...
   */
  /*  for ( p=WIM; p<DIM; p++)
      {
      v = PVELOCITY1 + p;
      if ( pd->v[pg->imtrx][v] || (upd->vp[pg->imtrx][v] == -1) )
      {
      fv->pv[p]     = 0.;
      fv_old->pv[p] = 0.;
      fv_dot->pv[p] = 0.;
      }
      } */

  for (p = 0; pdgv[PVELOCITY1] && p < WIM; p++) {
    v = PVELOCITY1 + p;
    if (pdgv[v]) {
      fv->pv[p] = 0.0;
      fv_old->pv[p] = 0.0;
      fv_dot->pv[p] = 0.0;

      dofs = ei[upd->matrix_index[v]]->dof[v];
      for (i = 0; i < dofs; i++) {
        fv->pv[p] += *esp->pv[p][i] * bf[v]->phi[i];

        if (pd->TimeIntegration != STEADY) {
          fv_old->pv[p] += *esp_old->pv[p][i] * bf[v]->phi[i];
          fv_dot->pv[p] += *esp_dot->pv[p][i] * bf[v]->phi[i];
        }
      }
    }
  }

  /* Extension velocity */

  v = EXT_VELOCITY;
  if (pdgv[EXT_VELOCITY]) {
    scalar_fv_fill(esp->ext_v, esp_dot->ext_v, esp_old->ext_v, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->ext_v), &(fv_dot->ext_v),
                   &(fv_old->ext_v));
    stateVector[EXT_VELOCITY] = fv->ext_v;
  }

  /* Electric Field */

  /*
   * Default: all velocities are zero...
   */
  for (p = 0; pdgv[EFIELD1] && p < WIM; p++) {
    v = EFIELD1 + p;

    if (pdgv[v]) {
      fv->E_field[p] = 0.;

      dofs = ei[upd->matrix_index[v]]->dof[v];
      for (i = 0; i < dofs; i++) {
        fv->E_field[p] += *esp->E_field[p][i] * bf[v]->phi[i];
      }
    }
  }

  /* Phase functions */

  for (p = 0; pdgv[PHASE1] && p < pfd->num_phase_funcs; p++) {
    v = PHASE1 + p;
    if (pdgv[v]) {
      scalar_fv_fill(esp->pF[p], esp_dot->pF[p], esp_old->pF[p], bf[v]->phi,
                     ei[upd->matrix_index[v]]->dof[v], &(fv->pF[p]), &(fv_dot->pF[p]),
                     &(fv_old->pF[p]));
    }
  }

  /*
   * Polymer Stress (tensor)...
   */

  for (mode = 0; pdgv[POLYMER_STRESS11] && mode < vn->modes; mode++) {
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        v = v_s[mode][p][q];
        if (pdgv[v] || (upd->vp[pg->imtrx][v] == -1)) {
          /* good default behavior */
          fv->S[mode][p][q] = 0.;
          fv_old->S[mode][p][q] = 0.;
          fv_dot->S[mode][p][q] = 0.;
        }
      }
    }
  }

  for (mode = 0; pdgv[POLYMER_STRESS11] && mode < vn->modes; mode++) {
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        if (p <= q) {
          v = v_s[mode][p][q];
          if (pdgv[v]) {
            dofs = ei[upd->matrix_index[v]]->dof[v];
            for (i = 0; i < dofs; i++) {
              fv->S[mode][p][q] += *esp->S[mode][p][q][i] * bf[v]->phi[i];
              if (pd->TimeIntegration != STEADY) {
                fv_old->S[mode][p][q] += *esp_old->S[mode][p][q][i] * bf[v]->phi[i];
                fv_dot->S[mode][p][q] += *esp_dot->S[mode][p][q][i] * bf[v]->phi[i];
              }
            }
          }
          /* form the entire symmetric stress matrix for the momentum equation */
          fv->S[mode][q][p] = fv->S[mode][p][q];
          fv_old->S[mode][q][p] = fv_old->S[mode][p][q];
          fv_dot->S[mode][q][p] = fv_dot->S[mode][p][q];
        }
      }
    }
  }

  /*
   * Velocity Gradient (tensor)...
   */

  for (p = 0; pdgv[VELOCITY_GRADIENT11] && p < VIM; p++) {
    if (gn->ConstitutiveEquation == BINGHAM_MIXED) {
      for (q = 0; q < VIM; q++) {
        if (q >= p) {
          v = v_g[p][q];
          if (pdgv[v]) {
            fv->G[p][q] = fv_old->G[p][q] = fv_dot->G[p][q] = 0.0;
            if (q > p) {
              fv->G[q][p] = fv_old->G[q][p] = fv_dot->G[q][p] = 0.0;
            }
            dofs = ei[upd->matrix_index[v]]->dof[v];
            for (i = 0; i < dofs; i++) {
              fv->G[p][q] += *esp->G[p][q][i] * bf[v]->phi[i];
              if (q > p) {
                fv->G[q][p] += *esp->G[p][q][i] * bf[v]->phi[i];
              }
              if (pd->TimeIntegration != STEADY) {
                fv_old->G[p][q] += *esp_old->G[p][q][i] * bf[v]->phi[i];
                if (q > p) {
                  fv_dot->G[q][p] += *esp_dot->G[p][q][i] * bf[v]->phi[i];
                }
              }
            }
          }
        }
      }
    } else {
      for (q = 0; q < VIM; q++) {
        v = v_g[p][q];
        if (pdgv[v]) {
          fv->G[p][q] = fv_old->G[p][q] = fv_dot->G[p][q] = 0.0;
          dofs = ei[upd->matrix_index[v]]->dof[v];
          for (i = 0; i < dofs; i++) {
            fv->G[p][q] += *esp->G[p][q][i] * bf[v]->phi[i];
            if (pd->TimeIntegration != STEADY) {
              fv_old->G[p][q] += *esp_old->G[p][q][i] * bf[v]->phi[i];
              fv_dot->G[p][q] += *esp_dot->G[p][q][i] * bf[v]->phi[i];
            }
          }
        }
      }
    }
  }

  /*
   * Species Unknown Variable
   *      Modifications for non-dilute systems:
   *      -> We need to calculate the mole fraction of the last species
   *         in the mechanism, even if there isn't an equation
   *         for it. This is done via the sum MF = 1
   *         constraint, here. Generalization to other equations
   *         of state should be added here in the future.
   */

  if (pdgv[MASS_FRACTION]) {
    v = MASS_FRACTION;
    if (pd->Num_Species_Eqn != pd->Num_Species) {
      N = pd->Num_Species - 1;

      fv->c[N] = fv_old->c[N] = fv_dot->c[N] = 0.;
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        scalar_fv_fill(esp->c[w], esp_dot->c[w], esp_old->c[w], bf[v]->phi,
                       ei[upd->matrix_index[v]]->dof[v], &(fv->c[w]), &(fv_dot->c[w]),
                       &(fv_old->c[w]));
        fv->c[N] -= fv->c[w];
        fv_old->c[N] -= fv_old->c[w];
        fv_dot->c[N] -= fv_dot->c[w];
        stateVector[SPECIES_UNK_0 + w] = fv->c[w];
      }
      switch (mp->Species_Var_Type) {
      case SPECIES_UNDEFINED_FORM:
      case SPECIES_MOLE_FRACTION:
      case SPECIES_MASS_FRACTION:
      case SPECIES_VOL_FRACTION:
        fv->c[N] += 1.0;
        stateVector[SPECIES_UNK_0 + N] = fv->c[N];
        break;
      case SPECIES_CONCENTRATION:
        /*
         *  HKM:  These are currently problematic cases that
         *  argue for including the last species as a dependent
         *  variable in the solution vector. The reason is that
         *  a small nonlinear iteration must be performed here
         *  for calculation of the total species concentration
         *  and the species concentration of the last species
         *  in the mechanism. An equation of state has to
         *  be assumed, and then inverted.
         *
         *   C = Sum_i=1,N[C_i],
         *
         *  where C is determined from the equation of state
         *  as a function of T, P, [other_state_variables], and
         *  C_1 to C_N.
         *
         *  If the problem is
         *  dilute, then the whole process is circumvented:
         *     C = C_N.
         *  or at least C isn't a function of the last species
         *  in the mechanism:
         *     C = C(T, P, C_1, ..., C_N-1)
         *
         *  The dilute case is carried out below. However,
         *  note that c[N] = 0 is assumed in the function
         *  call. If c[N] is used, then this routine will
         *  give the wrong answer, until a nonlinear iteration
         *  loop is installed below.
         *
         * -> The same is true for the SPECIES_DENSITY
         *    loop.
         */
        switch (mp->DensityModel) {
        case DENSITY_CONSTANT_LAST_CONC:
          fv->c[N] = mp->u_density[0];
          break;
        default:
          rho = calc_concentration(mp, FALSE, NULL);
          fv->c[N] += rho;
        }
        stateVector[SPECIES_UNK_0 + N] = fv->c[N];
        break;
      case SPECIES_DENSITY:
        /* note, this won't work for time dependent density models, but I got tired ... RRR*/
        rho = calc_density(mp, FALSE, NULL, 0.);
        fv->c[N] += rho;
        stateVector[SPECIES_UNK_0 + N] = fv->c[N];
        break;
      case SPECIES_CAP_PRESSURE:
        break;
      default:
        break;
      }
      mp->StateVector_speciesVT = upd->Species_Var_Type;
    } else {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        scalar_fv_fill(esp->c[w], esp_dot->c[w], esp_old->c[w], bf[v]->phi,
                       ei[upd->matrix_index[v]]->dof[v], &(fv->c[w]), &(fv_dot->c[w]),
                       &(fv_old->c[w]));
        fv_dot_old->c[w] = 0;
        for (int i = 0; i < ei[upd->matrix_index[v]]->dof[v]; i++) {
          fv_dot_old->c[w] += *(pg->matrices[upd->matrix_index[v]].xdot_old -
                                pg->matrices[upd->matrix_index[v]].xdot + esp_dot->c[w][i]) *
                              bf[v]->phi[i];
        }
        stateVector[SPECIES_UNK_0 + w] = fv->c[w];
      }
    }
  } /* HKM->
     *  (we may want to insert reference concentrations here
     *   if available)
     */

  /*else if (upd->vp[pg->imtrx][v] == -1) {

  for (w = 0; w < pd->Num_Species; w++) {
  fv->c[w]     = 0.;
  fv_old->c[w] = 0.;
  fv_dot->c[w] = 0.;
  stateVector[SPECIES_UNK_0+w] = fv->c[w];
  }
  } */

  /*
   * Porous media Variables
   */

  /*  if (pd->v[pg->imtrx][v] || (upd->vp[pg->imtrx][v] == -1) )
      {
      fv->p_liq=0.;
      fv_old->p_liq = 0.;
      fv_dot->p_liq = 0.;
      fv_dot_old->p_liq = 0.0;
      } */

  if (pdgv[POR_LIQ_PRES]) {
    v = POR_LIQ_PRES;
    fv->p_liq = 0.;
    fv_old->p_liq = 0.;
    fv_dot->p_liq = 0.;
    fv_dot_old->p_liq = 0.0;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (i = 0; i < dofs; i++) {
      fv->p_liq += *esp->p_liq[i] * bf[v]->phi[i];
      if (pd->TimeIntegration != STEADY) {
        fv_old->p_liq += *esp_old->p_liq[i] * bf[v]->phi[i];
        fv_dot->p_liq += *esp_dot->p_liq[i] * bf[v]->phi[i];

        if (upd->Total_Num_Matrices > 1) {
          fv_dot_old->p_liq += *(pg->matrices[upd->matrix_index[v]].xdot_old -
                                 pg->matrices[upd->matrix_index[v]].xdot + esp_dot->p_liq[i]) *
                               bf[v]->phi[i];
        } else {
          fv_dot_old->p_liq += *(xdot_old_static - xdot_static + esp_dot->p_liq[i]) * bf[v]->phi[i];
        }
      }
    }
  }

  /*  if (pd->v[pg->imtrx][v] || (upd->vp[pg->imtrx][v] == -1) )
      {
      fv->p_gas=0.;
      fv_old->p_gas = 0.;
      fv_dot->p_gas = 0.;
      fv_dot_old->p_gas = 0.0;
      }
  */
  if (pdgv[POR_GAS_PRES]) {
    v = POR_GAS_PRES;
    fv->p_gas = 0.;
    fv_old->p_gas = 0.;
    fv_dot->p_gas = 0.;
    fv_dot_old->p_gas = 0.0;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (i = 0; i < dofs; i++) {
      fv->p_gas += *esp->p_gas[i] * bf[v]->phi[i];
      if (pd->TimeIntegration != STEADY) {
        fv_old->p_gas += *esp_old->p_gas[i] * bf[v]->phi[i];
        fv_dot->p_gas += *esp_dot->p_gas[i] * bf[v]->phi[i];
        if (upd->Total_Num_Matrices > 1) {
          fv_dot_old->p_gas += *(pg->matrices[upd->matrix_index[v]].xdot_old -
                                 pg->matrices[upd->matrix_index[v]].xdot + esp_dot->p_gas[i]) *
                               bf[v]->phi[i];
        } else {
          fv_dot_old->p_gas += *(xdot_old_static - xdot_static + esp_dot->p_gas[i]) * bf[v]->phi[i];
        }
      }
    }
  }

  /*  if (pd->v[pg->imtrx][v] || (upd->vp[pg->imtrx][v] == -1) )
      {
      fv->porosity=0.;
      fv_old->porosity = 0.;
      fv_dot->porosity = 0.;
      fv_dot_old->porosity = 0.0;
      }
  */
  if (pdgv[POR_POROSITY]) {
    v = POR_POROSITY;
    fv->porosity = 0.;
    fv_old->porosity = 0.;
    fv_dot->porosity = 0.;
    fv_dot_old->porosity = 0.0;

    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (i = 0; i < dofs; i++) {
      fv->porosity += *esp->porosity[i] * bf[v]->phi[i];
      if (pd->TimeIntegration != STEADY) {
        fv_old->porosity += *esp_old->porosity[i] * bf[v]->phi[i];
        fv_dot->porosity += *esp_dot->porosity[i] * bf[v]->phi[i];
        if (upd->Total_Num_Matrices > 1) {
          fv_dot_old->porosity +=
              *(pg->matrices[upd->matrix_index[v]].xdot_old -
                pg->matrices[upd->matrix_index[v]].xdot + esp_dot->porosity[i]) *
              bf[v]->phi[i];
        } else {
          fv_dot_old->porosity +=
              *(xdot_old_static - xdot_static + esp_dot->porosity[i]) * bf[v]->phi[i];
        }
      }
    }
  }

  if (pdgv[POR_TEMP]) {
    v = POR_TEMP;
    fv->T = 0.;
    fv_old->T = 0.;
    fv_dot->T = 0.;
    fv_dot_old->T = 0.0;

    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (i = 0; i < dofs; i++) {
      fv->T += *esp->T[i] * bf[v]->phi[i];
      if (pd->TimeIntegration != STEADY) {
        fv_old->T += *esp_old->T[i] * bf[v]->phi[i];
        fv_dot->T += *esp_dot->T[i] * bf[v]->phi[i];

        if (upd->Total_Num_Matrices > 1) {
          fv_dot_old->T += *(pg->matrices[upd->matrix_index[v]].xdot_old -
                             pg->matrices[upd->matrix_index[v]].xdot + esp_dot->T[i]) *
                           bf[v]->phi[i];
        } else {
          fv_dot_old->T += *(xdot_old_static - xdot_static + esp_dot->T[i]) * bf[v]->phi[i];
        }
      }
    }
  }

  /*
   * Vorticity principle flow direction
   */
  for (p = 0; pdgv[VORT_DIR1] && p < DIM; p++) {
    v = VORT_DIR1 + p;
    /*if (pd->v[pg->imtrx][v] || upd->vp[pg->imtrx][v] == -1) fv->vd[p] = 0.0; */

    if (pdgv[v]) {
      dofs = ei[upd->matrix_index[v]]->dof[v];
      fv->vd[p] = 0.0;
      for (i = 0; i < dofs; i++)
        fv->vd[p] += *esp->vd[p][i] * bf[v]->phi[i];
    }
  }

  /*
   * Lagrange Multiplier Field
   */
  for (p = 0; pdgv[LAGR_MULT1] && p < DIM; p++) {
    v = LAGR_MULT1 + p;
    /*      if (pd->v[pg->imtrx][v] || upd->vp[pg->imtrx][v] == -1)
            {
            fv->lm[p] = 0.0;
            fv_old->lm[p] = 0.0;
            }
    */
    if (pdgv[v]) {
      fv->lm[p] = 0.0;
      fv_old->lm[p] = 0.0;

      dofs = ei[upd->matrix_index[v]]->dof[v];
      for (i = 0; i < dofs; i++) {
        fv->lm[p] += *esp->lm[p][i] * bf[v]->phi[i];
        if (pd->TimeIntegration != STEADY) {
          fv_old->lm[p] += *esp_old->lm[p][i] * bf[v]->phi[i];
        }
      }
    }
  }

  /*
   * Eigenvalue associated with vd.
   */

  /*  if (pd->v[pg->imtrx][v] || upd->vp[pg->imtrx][v] == -1) fv->vlambda = 0.0; */

  if (pdgv[VORT_LAMBDA]) {
    v = VORT_LAMBDA;
    fv->vlambda = 0.0;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (i = 0; i < dofs; i++)
      fv->vlambda += *esp->vlambda[i] * bf[v]->phi[i];
  }

  for (p = 0; pdgv[MOMENT0] && p < MAX_MOMENTS; p++) {
    v = MOMENT0 + p;
    if (pdgv[v]) {
      scalar_fv_fill(esp->moment[p], esp_dot->moment[p], esp_old->moment[p], bf[v]->phi,
                     ei[upd->matrix_index[v]]->dof[v], &(fv->moment[p]), &(fv_dot->moment[p]),
                     &(fv_old->moment[p]));

      fv_dot_old->moment[p] = 0;
      for (int i = 0; i < ei[upd->matrix_index[v]]->dof[v]; i++) {
        if (upd->Total_Num_Matrices > 1) {
          fv_dot_old->moment[p] +=
              *(pg->matrices[upd->matrix_index[v]].xdot_old -
                pg->matrices[upd->matrix_index[v]].xdot + esp_dot->moment[p][i]) *
              bf[v]->phi[i];
        } else {
          fv_dot_old->moment[p] +=
              *(xdot_old_static - xdot_static + esp_dot->moment[p][i]) * bf[v]->phi[i];
        }
      }
    }
  }

  if (pdgv[DENSITY_EQN]) {
    v = DENSITY_EQN;
    scalar_fv_fill(esp->rho, esp_dot->rho, esp_old->rho, bf[v]->phi,
                   ei[upd->matrix_index[v]]->dof[v], &(fv->rho), &(fv_dot->rho), &(fv_old->rho));
    stateVector[v] = fv->rho;
  }

  if (pdgv[TFMP_PRES]) {

    v = TFMP_PRES;
    scalar_fv_fill(esp->tfmp_pres, esp_dot->tfmp_pres, esp_old->tfmp_pres, bf[v]->phi,
                   ei[pg->imtrx]->dof[v], &(fv->tfmp_pres), &(fv_dot->tfmp_pres),
                   &(fv_old->tfmp_pres));
    stateVector[v] = fv->tfmp_pres;
  }
  if (pdgv[TFMP_SAT]) {
    v = TFMP_SAT;
    scalar_fv_fill(esp->tfmp_sat, esp_dot->tfmp_sat, esp_old->tfmp_sat, bf[v]->phi,
                   ei[pg->imtrx]->dof[v], &(fv->tfmp_sat), &(fv_dot->tfmp_sat),
                   &(fv_old->tfmp_sat));
    stateVector[v] = fv->tfmp_sat;
  }

  /*
   * EM Wave Vectors...
   */
  if (pdgv[EM_E1_REAL] || pdgv[EM_E2_REAL] || pdgv[EM_E3_REAL]) {
    v = EM_E1_REAL;
    if (bf[v]->interpolation == I_N1) {
      for (p = 0; p < DIM; p++) {
        fv->em_er[p] = 0.0;
        fv_old->em_er[p] = 0.0;
        fv_dot->em_er[p] = 0.0;

        if (pdgv[v]) {
          dofs = ei[pg->imtrx]->dof[v];
          for (i = 0; i < dofs; i++) {
            fv->em_er[p] += *esp->em_er[0][i] * bf[v]->phi_e[i][p];

            if (pd->TimeIntegration != STEADY) {
              fv_old->em_er[p] += *esp_old->em_er[0][i] * bf[v]->phi_e[i][p];
              fv_dot->em_er[p] += *esp_dot->em_er[0][i] * bf[v]->phi_e[i][p];
            }
          }
        }
      }

    } else {
      for (p = 0; p < DIM; p++) {
        v = EM_E1_REAL + p;
        fv->em_er[p] = 0.0;
        fv_old->em_er[p] = 0.0;
        fv_dot->em_er[p] = 0.0;

        if (pdgv[v]) {
          dofs = ei[pg->imtrx]->dof[v];
          for (i = 0; i < dofs; i++) {
            fv->em_er[p] += *esp->em_er[p][i] * bf[v]->phi[i];

            if (pd->TimeIntegration != STEADY) {
              fv_old->em_er[p] += *esp_old->em_er[p][i] * bf[v]->phi[i];
              fv_dot->em_er[p] += *esp_dot->em_er[p][i] * bf[v]->phi[i];
            }
          }
        }
      }
    }
  }

  if (pdgv[EM_E1_IMAG] || pdgv[EM_E2_IMAG] || pdgv[EM_E3_IMAG]) {
    v = EM_E1_IMAG;
    if (bf[v]->interpolation == I_N1) {
      for (p = 0; p < DIM; p++) {
        fv->em_ei[p] = 0.0;
        fv_old->em_ei[p] = 0.0;
        fv_dot->em_ei[p] = 0.0;

        if (pdgv[v]) {
          dofs = ei[pg->imtrx]->dof[v];
          for (i = 0; i < dofs; i++) {
            fv->em_ei[p] += *esp->em_ei[0][i] * bf[v]->phi_e[i][p];

            if (pd->TimeIntegration != STEADY) {
              fv_old->em_ei[p] += *esp_old->em_ei[0][i] * bf[v]->phi_e[i][p];
              fv_dot->em_ei[p] += *esp_dot->em_ei[0][i] * bf[v]->phi_e[i][p];
            }
          }
        }
      }

    } else {
      for (p = 0; p < DIM; p++) {
        v = EM_E1_IMAG + p;
        fv->em_ei[p] = 0.0;
        fv_old->em_ei[p] = 0.0;
        fv_dot->em_ei[p] = 0.0;

        if (pdgv[v]) {
          dofs = ei[pg->imtrx]->dof[v];
          for (i = 0; i < dofs; i++) {
            fv->em_ei[p] += *esp->em_ei[p][i] * bf[v]->phi[i];

            if (pd->TimeIntegration != STEADY) {
              fv_old->em_ei[p] += *esp_old->em_ei[p][i] * bf[v]->phi[i];
              fv_dot->em_ei[p] += *esp_dot->em_ei[p][i] * bf[v]->phi[i];
            }
          }
        }
      }
    }
  }

  if (pdgv[EM_H1_REAL] || pdgv[EM_H2_REAL] || pdgv[EM_H3_REAL]) {
    for (p = 0; p < DIM; p++) {
      v = EM_H1_REAL + p;
      fv->em_hr[p] = 0.0;
      fv_old->em_hr[p] = 0.0;
      fv_dot->em_hr[p] = 0.0;
      if (pdgv[v]) {
        dofs = ei[pg->imtrx]->dof[v];
        for (i = 0; i < dofs; i++) {
          fv->em_hr[p] += *esp->em_hr[p][i] * bf[v]->phi[i];

          if (pd->TimeIntegration != STEADY) {
            fv_old->em_hr[p] += *esp_old->em_hr[p][i] * bf[v]->phi[i];
            fv_dot->em_hr[p] += *esp_dot->em_hr[p][i] * bf[v]->phi[i];
          }
        }
      }
    }
  }

  if (pdgv[EM_H1_IMAG] || pdgv[EM_H2_IMAG] || pdgv[EM_H3_IMAG]) {
    for (p = 0; p < DIM; p++) {
      v = EM_H1_IMAG + p;
      fv->em_hi[p] = 0.0;
      fv_old->em_hi[p] = 0.0;
      fv_dot->em_hi[p] = 0.0;
      if (pdgv[v]) {
        dofs = ei[pg->imtrx]->dof[v];
        for (i = 0; i < dofs; i++) {
          fv->em_hi[p] += *esp->em_hi[p][i] * bf[v]->phi[i];

          if (pd->TimeIntegration != STEADY) {
            fv_old->em_hi[p] += *esp_old->em_hi[p][i] * bf[v]->phi[i];
            fv_dot->em_hi[p] += *esp_dot->em_hi[p][i] * bf[v]->phi[i];
          }
        }
      }
    }
  }

  // internal wall distance calculations
  if (upd->turbulent_info->use_internal_wall_distance) {
    fv->wall_distance = 0.;
    if (pdgv[pd->ShapeVar]) {
      dofs = ei[pg->imtrx]->dof[pd->ShapeVar];
      for (i = 0; i < dofs; i++) {
        fv->wall_distance +=
            upd->turbulent_info->wall_distances[ei[pg->imtrx]->gnn_list[pd->ShapeVar][i]] *
            bf[pd->ShapeVar]->phi[i];
      }
    }
  }

  /*
   * External...
   */
  if (efv->ev) {
    int table_id = 0;
    v = EXTERNAL;
    for (w = 0; w < efv->Num_external_field; w++) {
      dofs = ei[pg->imtrx]->dof_ext[w];
      fv->external_field[w] = 0.;
      fv_old->external_field[w] = 0.;
      fv_dot->external_field[w] = 0.;

      if (efv->i[w] != I_TABLE) {
        for (i = 0; i < dofs; i++) {
          fv->external_field[w] += *evp->external_field[w][i] * bfex[w]->phi[i];

          if (pd->TimeIntegration != STEADY) {
            fv_old->external_field[w] += *evp->external_field[w][i] * bfex[w]->phi[i];
            fv_dot->external_field[w] += 0.;
          }
        }
      } else {
        double slope;
        fv->external_field[w] = interpolate_table(ext_Tables[table_id], fv->x, &slope, NULL);
        table_id++;
      }

      /*
       * If the variable name is velocity, and the momentum equations are not active,
       * load the external_fields into the velocity fv for use in Advection-diffusion analysis
       */
      if (strcmp(efv->name[w], "VX") == 0 && !pd->v[pg->imtrx][VELOCITY1]) {
        fv->v[0] = fv->external_field[w];
        fv_old->v[0] = fv_old->external_field[w];
        fv_dot->v[0] = fv_dot->external_field[w];
      }
      if (strcmp(efv->name[w], "VY") == 0 && !pd->v[pg->imtrx][VELOCITY2]) {
        fv->v[1] = fv->external_field[w];
        fv_old->v[1] = fv_old->external_field[w];
        fv_dot->v[1] = fv_dot->external_field[w];
      }
      if (strcmp(efv->name[w], "VZ") == 0 && !pd->v[pg->imtrx][VELOCITY3]) {
        fv->v[2] = fv->external_field[w];
        fv_old->v[2] = fv_old->external_field[w];
        fv_dot->v[2] = fv_dot->external_field[w];
      }
      /*
       * If the variable name is mesh displacement, and the mesh
       * equations are not active, load the external_fields into
       * the mesh fv.  PRS NOTE: I don't think we will need these
       * for the decoupled JAS/GOMA, as the displacments are
       * swallowed into the mesh
       */
      if (strcmp(efv->name[w], "DMX") == 0 && !pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
        fv->d[0] = fv->external_field[w];
        fv_old->d[0] = fv_old->external_field[w];
        fv_dot->d[0] = fv_dot->external_field[w];
      }
      if (strcmp(efv->name[w], "DMY") == 0 && !pd->v[pg->imtrx][MESH_DISPLACEMENT2]) {
        fv->d[1] = fv->external_field[w];
        fv_old->d[1] = fv_old->external_field[w];
        fv_dot->d[1] = fv_dot->external_field[w];
      }
      if (strcmp(efv->name[w], "DMZ") == 0 && !pd->v[pg->imtrx][MESH_DISPLACEMENT3]) {
        fv->d[2] = fv->external_field[w];
        fv_old->d[2] = fv_old->external_field[w];
        fv_dot->d[2] = fv_dot->external_field[w];
      }
      /*
       * If the variable name is porosity, and the porosity equation is not active,
       * load the external_fields into the porosity fv
       */
      if (strcmp(efv->name[w], "P_POR") == 0 && !pd->v[pg->imtrx][POR_POROSITY]) {
        fv->porosity = fv->external_field[w];
        fv_old->porosity = fv_old->external_field[w];
        fv_dot->porosity = fv_dot->external_field[w];
      }
      if (strcmp(efv->name[w], "F") == 0 && !pd->v[pg->imtrx][FILL]) {
        fv->F = fv->external_field[w];
        fv_old->F = fv_old->external_field[w];
        fv_dot->F = fv_dot->external_field[w];
      }
      if (strcmp(efv->name[w], "F1") == 0 && !pd->v[pg->imtrx][PHASE1]) {
        fv->pF[0] = fv->external_field[w];
        fv_old->pF[0] = fv_old->external_field[w];
        fv_dot->pF[0] = fv_dot->external_field[w];
      }
    }
    if (!upd->turbulent_info->use_internal_wall_distance) {
      int i_d_wall = mp->dist_wall_ext_field_index;
      double d = fv->external_field[i_d_wall];
      fv->wall_distance = d;
    }
  } else {
    for (w = 0; w < efv->Num_external_field; w++) {
      fv->external_field[w] = 0.;
      fv_old->external_field[w] = 0.;
    }
  }

  /* initial displacements for TALE anneals. all for KINEMATIC DISPLACEMENT BC */
  if (efv->TALE) {
    for (w = 0; w < dim; w++) {
      v = MESH_DISPLACEMENT1 + w;
      dofs = ei[pg->imtrx]->dof[v];
      fv->initial_displacements[w] = 0.;
      for (i = 0; i < dofs; i++) {
        fv->initial_displacements[w] += *evp->initial_displacements[w][i] * bf[v]->phi[i];
      }
      v = SOLID_DISPLACEMENT1 + w;
      dofs = ei[pg->imtrx]->dof[v];
      fv->initial_displacements[w + DIM] = 0.;
      for (i = 0; i < dofs; i++) {
        fv->initial_displacements[w + DIM] +=
            *evp->initial_displacements[w + DIM][i] * bf[v]->phi[i];
      }
    }
  } else {
    for (w = 0; w < dim; w++) {
      fv->initial_displacements[w] = 0.;
      fv->initial_displacements[w + DIM] = 0.;
    }
  }

  return (status);
}

int load_fv_vector(void)

/*******************************************************************************
 * load_fv() -- load up values of all relevant field variables at the
 *              current gauss pt
 *
 * input: (assume the appropriate parts of esp, bf, ei, pd are filled in)
 * ------
 *
 *
 * output: ( this routine fills in the following parts of the fv structure)
 * -------
 *		fv->
 *			T -- temperature	(scalar)
 *			v -- velocity		(vector)
 *			d -- mesh displacement	(vector)
 *			c -- concentration	(multiple scalars)
 *		      por -- porous media       (multiple scalars)
 *			P -- pressure		(scalar)
 *			S -- polymer stress	(tensor)
 *			G -- velocity gradient	(tensor)
 *                     pv -- particle velocity  (vector)
 *                     pG -- particle velocity gradient (tensor)
 *              mp->StateVector[] is filled in as well, for
 *                    pertinent entries that make up the specification of
 *                    the state of the material.
 *
 * NOTE: To accommodate shell elements, this function has been modified
 *       so that fv variables are not zeroed out when they are active
 *       on an element block other than the current one.
 *       The check done for variable v is then:
 *          if ( pd->v[pg->imtrx][v] || upd->vp[pg->imtrx][v] == -1 )
 *       In many cases below, this conditional zeroing is done in
 *       a separate small loop before the main one.
 *
 *
 * Return values:
 *		0 -- if things went OK
 *	       -1 -- if a problem occurred
 *
 * Created:	Fri Mar 18 06:44:41 MST 1994 pasacki@sandia.gov
 *
 * Modified:
 ***************************************************************************/
{
  int v; /* variable type indicator */
  int i; /* index */
  int p;
  int dofs; /* degrees of freedom for a var in the elem */
  int status = 0;
  int *pdgv = pd->gv;

  status = 0;

  /*
   * EM Wave Vectors...
   */
  if (pdgv[EM_E1_REAL] || pdgv[EM_E2_REAL] || pdgv[EM_E3_REAL]) {
    v = EM_E1_REAL;
    if (bf[v]->interpolation == I_N1) {
      for (p = 0; p < DIM; p++) {
        fv->em_er[p] = 0.0;
        fv_old->em_er[p] = 0.0;
        fv_dot->em_er[p] = 0.0;

        if (pdgv[v]) {
          dofs = ei[pg->imtrx]->dof[v];
          for (i = 0; i < dofs; i++) {
            fv->em_er[p] += *esp->em_er[0][i] * bf[v]->phi_e[i][p];

            if (pd->TimeIntegration != STEADY) {
              fv_old->em_er[p] += *esp_old->em_er[0][i] * bf[v]->phi_e[i][p];
              fv_dot->em_er[p] += *esp_dot->em_er[0][i] * bf[v]->phi_e[i][p];
            }
          }
        }
      }
    }
  }

  if (pdgv[EM_E1_IMAG] || pdgv[EM_E2_IMAG] || pdgv[EM_E3_IMAG]) {
    v = EM_E1_IMAG;
    if (bf[v]->interpolation == I_N1) {
      for (p = 0; p < DIM; p++) {
        fv->em_ei[p] = 0.0;
        fv_old->em_ei[p] = 0.0;
        fv_dot->em_ei[p] = 0.0;

        if (pdgv[v]) {
          dofs = ei[pg->imtrx]->dof[v];
          for (i = 0; i < dofs; i++) {
            fv->em_ei[p] += *esp->em_ei[0][i] * bf[v]->phi_e[i][p];

            if (pd->TimeIntegration != STEADY) {
              fv_old->em_ei[p] += *esp_old->em_ei[0][i] * bf[v]->phi_e[i][p];
              fv_dot->em_ei[p] += *esp_dot->em_ei[0][i] * bf[v]->phi_e[i][p];
            }
          }
        }
      }
    }
  }

  return (status);
}

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
int load_fv_grads(void)
/*******************************************************************************
 *
 * load_fv_grads() -- load relevant field variable gradients at this gauss pt
 *
 * input: (assume the appropriate parts of esp, bf, ei, pd are filled in)
 * ------
 *
 *
 * output: ( this routine fills in the following parts of the fv structure)
 * -------
 *		fv->
 *			grad_T -- temperature gradient
 *			grad_P -- pressure gradient
 *			grad_c -- species unknown gradient
 *			grad_F -- fill gradient
 *			grad_V -- voltage potential gradient
 *			div_v  -- divergence of velocity
 *			grad_v -- velocity gradient tensor
 *                      curl_v -- curl of velocity, a.k.a. vorticity
 *			div_d  -- divergence of displacement (dilatation)
 *			grad_d -- gradient of displacement ("strain")
 *                      grad_X -- gradient of the Eulerian solid reference state.
 *			grad_S -- gradient of the polymer stress tensor
 *			grad_G -- gradient of the velocity gradient tensor
 *                      grad_pv -- gradient of particle velocity
 *                      grad_p_liq, grad_p_gas, grad_porosity
 *                             -- gradient of porous media variables
 *                      grad_n -- gradient of level set normal vector
 *
 * NOTE: To accommodate shell elements, this function has been modified
 *       so that fv variables are not zeroed out when they are active
 *       on an element block other than the current one.
 *       The check done for variable v is then:
 *          if ( pd->v[pg->imtrx][v] || upd->vp[pg->imtrx][v] == -1 )
 *       In many cases below, this conditional zeroing is done in
 *       a separate small loop before the main one.
 *
 * Return values:
 *		0 -- if things went OK
 *	       -1 -- if a problem occurred
 *
 * Created:	Fri Mar 18 07:36:14 MST 1994 pasacki@sandia.gov
 *
 * Modified:	Tue Feb 21 11:08 MST 1995 pasacki@sandia.gov
 ***************************************************************************/
{
  int v;          /* variable type indicator */
  int i, a;       /* index */
  int p, q, r, s; /* dimension index */
  int dofs;       /* degrees of freedom */
  int w;          /* concentration species */
  int dim = pd->Num_Dim;
  int mode; /* modal counter */
  int status = 0;
  int transient_run = (pd->TimeIntegration != STEADY);
  BASIS_FUNCTIONS_STRUCT *bfn;

  /* Use a static flag so unused grads are zero on first call, but are not zero subsequently
   *  This is for efficieny
   */
  static int zero_unused_grads = FALSE;

  /*
   * grad(T)
   */
  if (pd->gv[TEMPERATURE]) {
    v = TEMPERATURE;
    for (p = 0; p < VIM; p++)
      fv->grad_T[p] = 0.0;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      for (i = 0; i < dofs; i++) {
        fv->grad_T[p] += *esp->T[i] * bf[v]->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][TEMPERATURE] == -1) {
    for (p = 0; p < VIM; p++)
      fv->grad_T[p] = 0.0;
  }

  /*
   * grad(P)
   */

  if (pd->gv[PRESSURE]) {
    v = PRESSURE;
    dofs = ei[upd->matrix_index[v]]->dof[v];
#ifdef DO_NO_UNROLL
    for (p = 0; p < VIM; p++) {
      fv->grad_P[p] = 0.0;

      for (i = 0; i < dofs; i++) {
        fv->grad_P[p] += *esp->P[i] * bf[v]->grad_phi[i][p];
      }
    }
#else
    grad_scalar_fv_fill(esp->P, bf[v]->grad_phi, dofs, fv->grad_P);
    // if(transient_run && segregated)
    if (transient_run) {
      grad_scalar_fv_fill(esp_old->P, bf[v]->grad_phi, dofs, fv_old->grad_P);
    }
#endif

  } else if (zero_unused_grads && upd->vp[pg->imtrx][PRESSURE] == -1) {
    for (p = 0; p < VIM; p++) {
      fv->grad_P[p] = 0.0;
      // if(transient_run && segregated)
      if (transient_run) {
        fv_old->grad_P[p] = 0.0;
      }
    }
  }

  if (pd->gv[PSTAR]) {
    v = PSTAR;
    dofs = ei[upd->matrix_index[v]]->dof[v];
#ifdef DO_NO_UNROLL
    for (p = 0; p < VIM; p++) {
      fv->grad_P[p] = 0.0;

      for (i = 0; i < dofs; i++) {
        fv->grad_P[p] += *esp->P[i] * bf[v]->grad_phi[i][p];
      }
    }
#else
    grad_scalar_fv_fill(esp->P_star, bf[v]->grad_phi, dofs, fv->grad_P_star);
    // if(transient_run && segregated)
    if (transient_run) {
      grad_scalar_fv_fill(esp_old->P_star, bf[v]->grad_phi, dofs, fv_old->grad_P_star);
    }
#endif

  } else if (zero_unused_grads && upd->vp[pg->imtrx][PSTAR] == -1) {
    for (p = 0; p < VIM; p++) {
      fv->grad_P[p] = 0.0;
      // if(transient_run && segregated)
      if (transient_run) {
        fv_old->grad_P[p] = 0.0;
      }
    }
  }

  /*
   * grad(nn)
   */

  if (pd->gv[BOND_EVOLUTION]) {
    v = BOND_EVOLUTION;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_nn[p] = 0.;
      for (i = 0; i < dofs; i++) {
        fv->grad_nn[p] += *esp->nn[i] * bf[v]->grad_phi[i][p];
      }
    }
  }

  /*
   * grad(F)
   */

  if (pd->gv[FILL]) {
    v = FILL;
    dofs = ei[upd->matrix_index[v]]->dof[v];
#ifdef DO_NO_UNROLL
    for (p = 0; p < VIM; p++) {
      fv->grad_F[p] = 0.0;
      fv_old->grad_F[p] = 0.0;
      fv_dot->grad_F[p] = 0.0;

      for (i = 0; i < dofs; i++) {
        fv->grad_F[p] += *esp->F[i] * bf[v]->grad_phi[i][p];
        if (transient_run) {
          /* keep this only for VOF/Taylor-Galerkin stuff */
          fv_old->grad_F[p] += *esp_old->F[i] * bf[v]->grad_phi[i][p];

          fv_dot->grad_F[p] += *esp_dot->F[i] * bf[v]->grad_phi[i][p];
        }
      }
    }
#else
    grad_scalar_fv_fill(esp->F, bf[v]->grad_phi, dofs, fv->grad_F);

    if (transient_run) {
      grad_scalar_fv_fill(esp_old->F, bf[v]->grad_phi, dofs, fv_old->grad_F);
      grad_scalar_fv_fill(esp_dot->F, bf[v]->grad_phi, dofs, fv_dot->grad_F);
    }
#endif

  } else if (zero_unused_grads && upd->vp[pg->imtrx][FILL] == -1) {
    for (p = 0; p < VIM; p++) {
      fv->grad_F[p] = 0.0;
      fv_old->grad_F[p] = 0.0;
      fv_dot->grad_F[p] = 0.0;
    }
  }

  /*
   * grad(H)
   */
  if (pd->gv[CURVATURE]) {
    v = CURVATURE;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_H[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_H[p] += *esp->H[i] * bf[v]->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][CURVATURE] == -1) {
    for (p = 0; p < VIM; p++) {
      fv->grad_H[p] = 0.0;
    }
  }

  /*
   * grad(V)
   */
  if (pd->gv[VOLTAGE]) {
    v = VOLTAGE;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_V[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_V[p] += *esp->V[i] * bf[v]->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][VOLTAGE] == -1) {
    for (p = 0; p < VIM; p++)
      fv->grad_V[p] = 0.0;
  }

  /*
   * grad(Enorm)
   */
  if (pd->gv[ENORM]) {
    v = ENORM;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_Enorm[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_Enorm[p] += *esp->Enorm[i] * bf[v]->grad_phi[i][p];
      }
    }
  }

  /*
   * grad(qs)
   */
  if (pd->gv[SURF_CHARGE]) {
    v = SURF_CHARGE;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_qs[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        // HKM -> I changed dphidxi to grad_phi below. There didn't seem
        //        to be a case where the raw element derivative (dphidxi), which
        //        doesn't even take into account of the size of the element, should
        //        be used at this level.
        //   fv->grad_qs[p] += *esp->qs[i] * bf[v]->dphidxi[i][p];
        fv->grad_qs[p] += *esp->qs[i] * bf[v]->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][SURF_CHARGE] == -1) {
    for (p = 0; p < VIM - 1; p++) {
      fv->grad_qs[p] = 0.0;
    }
  }

  /*
   * grad(SH)
   */

  if (pd->gv[SHEAR_RATE]) {
    v = SHEAR_RATE;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_SH[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_SH[p] += *esp->SH[i] * bf[v]->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][SHEAR_RATE] == -1) {
    for (p = 0; p < VIM; p++)
      fv->grad_SH[p] = 0.0;
  }

  /*
   * grad(d)
   */

  if (pd->gv[MESH_DISPLACEMENT1] &&
      !ShellElementParentElementCoverageForVariable[MESH_DISPLACEMENT1]) {
    v = MESH_DISPLACEMENT1;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        fv->grad_d[p][q] = 0.0;
        fv_old->grad_d[p][q] = 0.0;
        for (r = 0; r < VIM; r++) {
          for (i = 0; i < dofs; i++) {
            fv->grad_d[p][q] += *esp->d[r][i] * bf[v]->grad_phi_e[i][r][p][q];
            fv_old->grad_d[p][q] += *esp_old->d[r][i] * bf[v]->grad_phi_e[i][r][p][q];
          }
        }
      }
    }
  } else {
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        fv->grad_d[p][q] = 0.0;
        fv_old->grad_d[p][q] = 0.0;
        fv->grad_d_dot[p][q] = 0.0;
      }
    }
  }

  /*
   * grad(d_dot)
   */
  if (pd->gv[MESH_DISPLACEMENT1] &&
      !ShellElementParentElementCoverageForVariable[MESH_DISPLACEMENT1]) {
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        v = MESH_DISPLACEMENT1 + p;
        if (pd->gv[v]) {
          dofs = ei[upd->matrix_index[v]]->dof[v];
          fv->grad_d_dot[p][q] = 0.;
          for (r = 0; r < dim; r++) {
            for (i = 0; i < dofs; i++) {
              fv->grad_d_dot[p][q] += *esp_dot->d[r][i] * bf[v]->grad_phi_e[i][r][p][q];
            }
          }
        }
      }
    }
  }

  /*
   * grad(d_rs)
   */
  if (pd->gv[SOLID_DISPLACEMENT1]) {
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        v = SOLID_DISPLACEMENT1 + p;

        if (pd->gv[v]) {
          dofs = ei[upd->matrix_index[v]]->dof[v];
          fv->grad_d_rs[p][q] = 0.0;
          for (r = 0; r < dim; r++) {
            for (i = 0; i < dofs; i++) {
              fv->grad_d_rs[p][q] += *esp->d_rs[r][i] * bf[v]->grad_phi_e[i][r][p][q];
            }
          }
        }
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][SOLID_DISPLACEMENT1] == -1) {
    memset(fv->grad_d_rs, 0, sizeof(double) * VIM * VIM);
  }

  /*
   * div(d)
   */

  /*
   * div(d_dot)
   */

  if (pd->gv[MESH_DISPLACEMENT1] &&
      !ShellElementParentElementCoverageForVariable[MESH_DISPLACEMENT1]) {
    v = MESH_DISPLACEMENT1;
    fv->div_d_dot = 0.0;
    fv->div_d = 0.0;

    for (p = 0; p < VIM; p++) {
      fv->div_d_dot += fv->grad_d_dot[p][p];
      fv->div_d += fv->grad_d[p][p];
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][MESH_DISPLACEMENT1] == -1) {
    fv->div_d_dot = 0.0;
    fv->div_d = 0.0;
  }

  /*
   * div(d_rs)
   */

  if (pd->gv[SOLID_DISPLACEMENT1]) {
    for (p = 0; p < VIM; p++) {
      fv->div_d_rs += fv->grad_d_rs[p][p];
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][SOLID_DISPLACEMENT1] == -1) {
    fv->div_d_rs = 0.;
  }

  /////+(*
  ////+(* grad(v)
  ////+(*/
  if (pd->gv[VELOCITY1]) {
    v = VELOCITY1;
    dofs = ei[upd->matrix_index[v]]->dof[v];

    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        fv->grad_v[p][q] = 0.0;
        fv_old->grad_v[p][q] = 0.0;
        for (r = 0; r < WIM; r++) {
          for (i = 0; i < dofs; i++) {
            fv->grad_v[p][q] += (*esp->v[r][i]) * bf[v]->grad_phi_e[i][r][p][q];
            if (pd->TimeIntegration != STEADY) {
              fv_old->grad_v[p][q] += (*esp_old->v[r][i]) * bf[v]->grad_phi_e[i][r][p][q];
            }
          }
        }
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][VELOCITY1] == -1) {
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        fv->grad_v[p][q] = 0.0;

        // if(segregated)
        fv_old->grad_v[p][q] = 0.0;
      }
    }
  }

  if (pd->gv[USTAR]) {
    v = USTAR;
    dofs = ei[upd->matrix_index[v]]->dof[v];

    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        fv->grad_v_star[p][q] = 0.0;
        for (r = 0; r < WIM; r++) {
          for (i = 0; i < dofs; i++) {
            fv->grad_v_star[p][q] += (*esp->v_star[r][i]) * bf[v]->grad_phi_e[i][r][p][q];
          }
        }
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][USTAR] == -1) {
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        fv->grad_v_star[p][q] = 0.0;
      }
    }
  }

  /*
   * grad(pv), particle velocity gradients.
   */
  if (pd->gv[PVELOCITY1]) {
    v = PVELOCITY1;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        for (r = 0; r < WIM; r++) {
          for (i = 0; i < dofs; i++) {
            fv->grad_pv[p][q] += *esp->pv[r][i] * bf[v]->grad_phi_e[i][r][p][q];
          }
        }
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][PVELOCITY1] == -1) {
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        fv->grad_pv[p][q] = 0.0;
      }
    }
  }

  /*
   * grad(ext_v), extension velocity gradients.
   */

  if (pd->gv[EXT_VELOCITY]) {
    v = EXT_VELOCITY;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_ext_v[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_ext_v[p] += *esp->ext_v[i] * bf[v]->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][EXT_VELOCITY] == -1) {
    for (p = 0; p < VIM; p++) {
      fv->grad_ext_v[p] = 0.0;
    }
  }

  /*
   *  grad( n ) .  Level set or shell normal vector gradient
   */
  if (pd->gv[NORMAL1] || pd->gv[SHELL_NORMAL1]) {
    if (pd->gv[NORMAL1])
      v = NORMAL1;
    if (pd->gv[SHELL_NORMAL1])
      v = SHELL_NORMAL1;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        fv->grad_n[p][q] = 0.0;
        fv_old->grad_n[p][q] = 0.0;
        for (r = 0; r < dim; r++) {
          for (i = 0; i < dofs; i++) {
            fv->grad_n[p][q] += *(esp->n[r][i]) * bf[v]->grad_phi_e[i][r][p][q];
            fv_old->grad_n[p][q] += *(esp_old->n[r][i]) * bf[v]->grad_phi_e[i][r][p][q];
          }
        }
      }
    }
  } else if ((zero_unused_grads && upd->vp[pg->imtrx][NORMAL1] == -1) ||
             upd->vp[pg->imtrx][SHELL_NORMAL1] == -1) {
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        fv->grad_n[p][q] = 0.0;
        fv_old->grad_n[p][q] = 0.0;
      }
    }
  }

  /*
   *  div(n)
   */
  if (pd->gv[NORMAL1] || pd->gv[SHELL_NORMAL1]) {
    fv->div_n = 0.0;
    fv_old->div_n = 0.0;
    for (p = 0; p < VIM; p++) {
      fv->div_n += fv->grad_n[p][p];
      fv_old->div_n += fv_old->grad_n[p][p];
    }
  }

  /*
   * div_s(n)
   */
  if (pd->gv[NORMAL1] || pd->gv[SHELL_NORMAL1]) {
    fv->div_s_n = 0.0;
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        fv->div_s_n += (delta(p, q) * fv->grad_n[p][q] - fv->n[p] * fv->n[q] * fv->grad_n[q][p]);
      }
    }
  }

  // Calculation of the surface curvature dyadic
  //             = - del_s (n)  = - (I - n n ) grad n
  if (pd->gv[SHELL_NORMAL1]) {
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        fv->surfCurvatureDyadic[p][q] = -fv->grad_n[p][q];
        for (a = 0; a < VIM; a++) {
          fv->surfCurvatureDyadic[p][q] += fv->n[p] * fv->n[a] * fv->grad_n[a][q];
        }
      }
    }
  }

  /*
   * grad(E_field), Electric field gradients.
   */
  if (pd->gv[EFIELD1]) {
    v = EFIELD1;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    v = EFIELD1;

    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        fv->grad_E_field[p][q] = 0.0;
        for (r = 0; r < WIM; r++) {
          for (i = 0; i < dofs; i++) {
            fv->grad_E_field[p][q] += *esp->E_field[r][i] * bf[v]->grad_phi_e[i][r][p][q];
          }
        }
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][EFIELD1] != -1) {
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        fv->grad_E_field[p][q] = 0.0;
      }
    }
  }

  /* Phase Functions */

  if (pd->gv[PHASE1]) {
    for (r = 0; r < pfd->num_phase_funcs; r++) {
      v = PHASE1 + r;
      dofs = ei[upd->matrix_index[v]]->dof[v];

      for (p = 0; p < VIM; p++) {
        fv->grad_pF[r][p] = 0.0;

        for (i = 0; i < dofs; i++) {
          fv->grad_pF[r][p] += *esp->pF[r][i] * bf[v]->grad_phi[i][p];
        }
      }
    }
  }

  /*
   * curl(v)
   */
  if (CURL_V != -1 && !InShellElementWithParentElementCoverage) {
    v = VELOCITY1;
    bfn = bf[v];

    if (pd->gv[v] || upd->vp[pg->imtrx][v] != -1) {
      for (p = 0; p < DIM; p++)
        fv->curl_v[p] = 0.0;
    }

    if (pd->gv[v]) {
      /* Always compute all three components of the vorticity
       * vector.  If we are really in 2D, then the output routine will
       * only output the 3rd component. DIM is equal to 3.
       */
      dofs = ei[upd->matrix_index[v]]->dof[VELOCITY1];
      for (p = 0; p < DIM; p++) {
        for (i = 0; i < dofs; i++) {
          for (a = 0; a < VIM; a++) /* VIM */
          {
            fv->curl_v[p] += *esp->v[a][i] * bfn->curl_phi_e[i][a][p];
          }
        }
      }
    }
  }

  if (pd->gv[EM_E1_REAL]) {
    v = EM_E1_REAL;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    bfn = bf[v];

    if (bfn->interpolation == I_N1) {
      for (p = 0; p < DIM; p++) {
        fv->curl_em_er[p] = 0.0;
        for (i = 0; i < dofs; i++) {
          fv->curl_em_er[p] += *esp->em_er[0][i] * bfn->curl_phi[i][p];
        }
      }
    } else {
      for (p = 0; p < DIM; p++) {
        fv->curl_em_er[p] = 0.0;
        for (i = 0; i < dofs; i++) {
          for (a = 0; a < dim; a++) {
            fv->curl_em_er[p] += *esp->em_er[a][i] * bfn->curl_phi_e[i][a][p];
          }
        }
      }
    }
  }

  if (pd->gv[EM_E1_IMAG]) {
    v = EM_E1_IMAG;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    bfn = bf[v];
    if (bfn->interpolation == I_N1) {
      for (p = 0; p < DIM; p++) {
        fv->curl_em_ei[p] = 0.0;
        for (i = 0; i < dofs; i++) {
          fv->curl_em_ei[p] += *esp->em_ei[0][i] * bfn->curl_phi[i][p];
        }
      }
    } else {
      for (p = 0; p < DIM; p++) {
        fv->curl_em_ei[p] = 0.0;
        for (i = 0; i < dofs; i++) {
          for (a = 0; a < dim; a++) {
            fv->curl_em_ei[p] += *esp->em_ei[a][i] * bfn->curl_phi_e[i][a][p];
          }
        }
      }
    }
  }

  /*
   * div(v)
   */

  if (pd->gv[VELOCITY1]) {
    fv->div_v = 0.0;
    if (VIM == 3 && pd->CoordinateSystem != CARTESIAN_2pt5D)
      fv_old->div_v = 0.0;
    for (a = 0; a < VIM; a++) {
      fv->div_v += fv->grad_v[a][a];
      // if(segregated)
      fv_old->div_v += fv_old->grad_v[a][a];
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][VELOCITY1] == -1) {
    fv->div_v = 0.0;
    fv_old->div_v = 0.0;
  }

  if (pd->gv[USTAR]) {
    fv->div_v_star = 0.0;
    for (a = 0; a < VIM; a++) {
      fv->div_v_star += fv->grad_v_star[a][a];
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][USTAR] == -1) {
    fv->div_v_star = 0.0;
  }

  /*
   * div(pv)
   */

  if (pd->gv[PVELOCITY1]) {
    v = PVELOCITY1;
    fv->div_pv = 0.0;
    for (p = 0; p < VIM; p++) {
      fv->div_pv += fv->grad_pv[p][p];
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][PVELOCITY1] == -1) {
    fv->div_pv = 0.0;
  }

  /*
   * Species unknown ...
   */

  if (pd->gv[MASS_FRACTION]) {
    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      v = MASS_FRACTION;
      dofs = ei[upd->matrix_index[v]]->dof[v];
      for (p = 0; p < VIM; p++) {
        fv->grad_c[w][p] = 0.;
        fv_old->grad_c[w][p] = 0.;

        for (i = 0; i < dofs; i++) {
          fv->grad_c[w][p] += *esp->c[w][i] * bf[v]->grad_phi[i][p];
          if (pd->TimeIntegration != STEADY) {
            /* keep this only for VOF/Taylor-Galerkin stuff */
            fv_old->grad_c[w][p] += *esp_old->c[w][i] * bf[v]->grad_phi[i][p];
          }
        }
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][MASS_FRACTION] == -1) {
    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      for (p = 0; p < VIM; p++) {
        fv->grad_c[w][p] = 0.;
        fv_old->grad_c[w][p] = 0.;
      }
    }
  }

  /*
   * Porous media unknowns ...
   */

  if (pd->gv[POR_LIQ_PRES]) {
    v = POR_LIQ_PRES;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_p_liq[p] = 0.0;
      fv_old->grad_p_liq[p] = 0.0;

      for (i = 0; i < dofs; i++) {
        fv->grad_p_liq[p] += *esp->p_liq[i] * bf[v]->grad_phi[i][p];
        if (pd->TimeIntegration != STEADY) {
          /* keep this only for VOF/Taylor-Galerkin stuff */
          fv_old->grad_p_liq[p] += *esp_old->p_liq[i] * bf[v]->grad_phi[i][p];
        }
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][POR_LIQ_PRES] == -1) {
    for (p = 0; p < VIM; p++) {
      fv->grad_p_liq[p] = 0.0;
      fv_old->grad_p_liq[p] = 0.0;
    }
  }

  if (pd->gv[POR_GAS_PRES]) {
    v = POR_GAS_PRES;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_p_gas[p] = 0.0;
      fv_old->grad_p_gas[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_p_gas[p] += *esp->p_gas[i] * bf[v]->grad_phi[i][p];
        if (pd->TimeIntegration != STEADY) {
          /* keep this only for VOF/Taylor-Galerkin stuff */
          fv_old->grad_p_gas[p] += *esp_old->p_gas[i] * bf[v]->grad_phi[i][p];
        }
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][POR_GAS_PRES] == -1) {
    for (p = 0; p < VIM; p++) {
      fv->grad_p_gas[p] = 0.0;
      fv_old->grad_p_gas[p] = 0.0;
    }
  }

  if (pd->gv[POR_POROSITY]) {
    v = POR_POROSITY;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_porosity[p] = 0.0;
      fv_old->grad_porosity[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_porosity[p] += *esp->porosity[i] * bf[v]->grad_phi[i][p];
        if (pd->TimeIntegration != STEADY) {
          /* keep this only for VOF/Taylor-Galerkin stuff */
          fv_old->grad_porosity[p] += *esp_old->porosity[i] * bf[v]->grad_phi[i][p];
        }
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][POR_POROSITY] == -1) {
    for (p = 0; p < VIM; p++) {
      fv->grad_porosity[p] = 0.0;
      fv_old->grad_porosity[p] = 0.0;
    }
  }

  if (pd->gv[POR_TEMP]) {
    v = POR_TEMP;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_T[p] = 0.;
      fv_old->grad_T[p] = 0.;
      for (i = 0; i < dofs; i++) {
        fv->grad_T[p] += *esp->T[i] * bf[v]->grad_phi[i][p];
        if (pd->TimeIntegration != STEADY) {
          /* keep this only for VOF/Taylor-Galerkin stuff */
          fv_old->grad_T[p] += *esp_old->T[i] * bf[v]->grad_phi[i][p];
        }
      }
    }
  }

  /*
   * grad(S), mode 0
   */

  if (pd->gv[POLYMER_STRESS11]) {
    v = POLYMER_STRESS11;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (mode = 0; mode < vn->modes; mode++) {
      for (p = 0; p < VIM; p++) {
        for (q = 0; q < VIM; q++) {
          for (r = 0; r < VIM; r++) {
            fv->grad_S[mode][r][p][q] = 0.;

            for (i = 0; i < dofs; i++) {
              if (p <= q) {
                fv->grad_S[mode][r][p][q] += *esp->S[mode][p][q][i] * bf[v]->grad_phi[i][r];
              } else {
                fv->grad_S[mode][r][p][q] += *esp->S[mode][q][p][i] * bf[v]->grad_phi[i][r];
              }
            }
          }
        }
      }
    }

    /*
     * div(S) - this is a vector!
     */
    for (mode = 0; mode < vn->modes; mode++) {
      for (r = 0; r < dim; r++) {
        fv->div_S[mode][r] = 0.0;

        for (q = 0; q < dim; q++) {
          fv->div_S[mode][r] += fv->grad_S[mode][q][q][r];
        }
      }
    }

    if (pd->CoordinateSystem != CARTESIAN) {
      for (mode = 0; mode < vn->modes; mode++) {
        for (s = 0; s < VIM; s++) {
          for (r = 0; r < VIM; r++) {
            for (p = 0; p < VIM; p++) {
              fv->div_S[mode][s] += fv->S[mode][p][s] * fv->grad_e[p][r][s];
            }
          }
        }
        for (s = 0; s < VIM; s++) {
          for (r = 0; r < VIM; r++) {
            for (q = 0; q < VIM; q++) {
              fv->div_S[mode][s] += fv->S[mode][r][q] * fv->grad_e[q][r][s];
            }
          }
        }
      }
    }

  } else if (zero_unused_grads && upd->vp[pg->imtrx][POLYMER_STRESS11] == -1) {
    for (mode = 0; mode < vn->modes; mode++) {
      for (p = 0; p < VIM; p++) {
        fv->div_S[mode][p] = 0.0;
        for (q = 0; q < VIM; q++) {
          for (r = 0; r < VIM; r++) {
            fv->grad_S[mode][r][p][q] = 0.;
          }
        }
      }
    }
  }

  /*
   * grad(G)
   */

  if (pd->gv[VELOCITY_GRADIENT11]) {
    v = VELOCITY_GRADIENT11;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        for (r = 0; r < VIM; r++) {
          fv->grad_G[r][p][q] = 0.0;
          fv->grad_Gt[r][p][q] = 0.0;

          for (i = 0; i < dofs; i++) {
            fv->grad_G[r][p][q] += *esp->G[p][q][i] * bf[v]->grad_phi[i][r];
            fv->grad_Gt[r][p][q] += *esp->G[q][p][i] * bf[v]->grad_phi[i][r];
          }
        }
      }
    }

    /*
     * div(G) - this is a vector!
     */
    for (r = 0; r < dim; r++) {
      fv->div_G[r] = 0.0;
      for (q = 0; q < dim; q++) {
        fv->div_G[r] += fv->grad_G[q][q][r];
      }
    }

    if (pd->CoordinateSystem != CARTESIAN) {
      for (s = 0; s < VIM; s++) {
        for (r = 0; r < VIM; r++) {
          for (p = 0; p < VIM; p++) {
            fv->div_G[s] += fv->G[p][s] * fv->grad_e[p][r][s];
          }
        }
      }

      for (s = 0; s < VIM; s++) {
        for (r = 0; r < VIM; r++) {
          for (q = 0; q < VIM; q++) {
            fv->div_G[s] += fv->G[r][q] * fv->grad_e[q][r][s];
          }
        }
      }
    }

    /*
     * div(Gt) - this is a vector and the divergence of the
     * transpose of G.
     */
    for (r = 0; r < dim; r++) {
      fv->div_Gt[r] = 0.0;
      for (q = 0; q < dim; q++) {
        fv->div_Gt[r] += fv->grad_G[q][r][q];
      }
    }

    if (pd->CoordinateSystem != CARTESIAN) {
      for (s = 0; s < VIM; s++) {
        for (r = 0; r < VIM; r++) {
          for (p = 0; p < VIM; p++) {
            fv->div_Gt[s] += fv->G[s][p] * fv->grad_e[p][r][s];
          }
        }
      }

      for (s = 0; s < VIM; s++) {
        for (r = 0; r < VIM; r++) {
          for (q = 0; q < VIM; q++) {
            fv->div_Gt[s] += fv->G[q][r] * fv->grad_e[q][r][s];
          }
        }
      }
    }

  } else if (zero_unused_grads && upd->vp[pg->imtrx][VELOCITY_GRADIENT11] == -1) {
    for (p = 0; p < VIM; p++) {
      fv->div_G[p] = 0.0;
      fv->div_Gt[p] = 0.0;
      for (q = 0; q < VIM; q++) {
        for (r = 0; r < VIM; r++) {
          fv->grad_G[r][p][q] = 0.;
        }
      }
    }
  }

  /*
   * grad(vd)
   */
  if (pd->gv[VORT_DIR1]) {
    v = VORT_DIR1;
    ;
    dofs = ei[upd->matrix_index[v]]->dof[v];

    for (p = 0; p < DIM; p++) {
      for (q = 0; q < DIM; q++) {
        fv->grad_vd[p][q] = 0.0;
        for (r = 0; r < DIM; r++) {
          for (i = 0; i < dofs; i++) {
            fv->grad_vd[p][q] += *esp->vd[r][i] * bf[v]->grad_phi_e[i][r][p][q];
          }
        }
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][VORT_DIR1] == -1) {
    for (p = 0; p < DIM; p++) {
      for (q = 0; q < DIM; q++) {
        fv->grad_vd[p][q] = 0.0;
      }
    }
  }
  /*
   * div(vd)
   */

  if (pd->gv[VORT_DIR1]) {
    fv->div_vd = 0.0;
    for (p = 0; p < VIM; p++) {
      fv->div_vd += fv->grad_vd[p][p];
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][VORT_DIR1] == -1)
    fv->div_vd = 0.0;

  /*
   * grad_sh_K
   * Gradient of curvature in structural shell
   */

  if (pd->v[pg->imtrx][SHELL_CURVATURE]) {
    v = SHELL_CURVATURE;
    dofs = ei[pg->imtrx]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_sh_K[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_sh_K[p] += *esp->sh_K[i] * bf[v]->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][SHELL_CURVATURE] == -1) {
    for (p = 0; p < VIM; p++)
      fv->grad_sh_K[p] = 0.0;
  }

  if (pd->v[pg->imtrx][SHELL_CURVATURE2]) {
    v = SHELL_CURVATURE2;
    dofs = ei[pg->imtrx]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_sh_K2[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_sh_K2[p] += *esp->sh_K2[i] * bf[v]->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][SHELL_CURVATURE2] == -1) {
    for (p = 0; p < VIM; p++)
      fv->grad_sh_K2[p] = 0.0;
  }

  /*
   * grad(n_dot_curl_s_v)
   *   This is the normal gradient of a scalar field defined on a shell.
   */
  if (pd->gv[N_DOT_CURL_V]) {
    v = N_DOT_CURL_V;
    bfn = bf[v];
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_n_dot_curl_s_v[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_n_dot_curl_s_v[p] += *esp->n_dot_curl_s_v[i] * bfn->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][N_DOT_CURL_V] == -1) {
    for (p = 0; p < VIM; p++) {
      fv->grad_n_dot_curl_s_v[p] = 0.0;
    }
  }

  /*
   * grad(div_s_v)
   *   This is the normal gradient of a scalar field defined on a shell.
   *        The scalar field is the surface divergence of the velocity
   */
  if (pd->gv[SHELL_SURF_DIV_V]) {
    v = SHELL_SURF_DIV_V;
    bfn = bf[v];
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_div_s_v[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_div_s_v[p] += *esp->div_s_v[i] * bfn->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][SHELL_SURF_DIV_V] == -1) {
    for (p = 0; p < VIM; p++) {
      fv->grad_div_s_v[p] = 0.0;
    }
  }

  /*
   * grad(curv)
   *   This is the normal gradient of a scalar field defined on a shell.
   *        The scalar field is the mean curvature of the surface
   */
  if (pd->gv[SHELL_SURF_CURV]) {
    v = SHELL_SURF_CURV;
    bfn = bf[v];
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_curv[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_curv[p] += *esp->curv[i] * bfn->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][SHELL_SURF_CURV] == -1) {
    for (p = 0; p < VIM; p++) {
      fv->grad_curv[p] = 0.0;
    }
  }

  /*
   * grad(grad_s_v_dot_n)
   *   This is the normal gradient of a scalar field defined on a shell.
   *        The scalar field is the vector component of the grad_s_v_dot_n
   */
  if (pd->gv[GRAD_S_V_DOT_N1]) {
    for (r = 0; r < dim; r++) {
      v = GRAD_S_V_DOT_N1 + r;
      bfn = bf[v];
      dofs = ei[upd->matrix_index[v]]->dof[v];
      for (p = 0; p < VIM; p++) {
        fv->serialgrad_grad_s_v_dot_n[r][p] = 0.0;
        for (i = 0; i < dofs; i++) {
          fv->serialgrad_grad_s_v_dot_n[r][p] += *(esp->grad_v_dot_n[r][i]) * bfn->grad_phi[i][p];
        }
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][GRAD_S_V_DOT_N1] == -1) {
    for (r = 0; r < dim; r++) {
      for (p = 0; p < VIM; p++) {
        fv->serialgrad_grad_s_v_dot_n[r][p] = 0.0;
      }
    }
  }

  /*
   * grad(sh_J)
   */
  if (pd->gv[SHELL_DIFF_FLUX]) {
    v = SHELL_DIFF_FLUX;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_sh_J[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        // HKM -> I changed dphidxi to grad_phi below. There didn't seem
        //        to be a case where the raw element derivative (dphidxi) which
        //        doesn't even take into account of the size of the element should
        //        be used at this level.
        //      fv->grad_sh_J[p] += *esp->sh_J[i] * bf[v]->dphidxi[i][p];
        fv->grad_sh_J[p] += *esp->sh_J[i] * bf[v]->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][SHELL_DIFF_FLUX] == -1) {
    for (p = 0; p < VIM; p++)
      fv->grad_sh_J[p] = 0.0;
  }

  /* grad(EDDY_NU)
   *
   */
  if (pd->gv[EDDY_NU]) {
    v = EDDY_NU;
    bfn = bf[v];
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_eddy_nu[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_eddy_nu[p] += *esp->eddy_nu[i] * bfn->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][EDDY_NU] == -1) {
    for (p = 0; p < VIM; p++) {
      fv->grad_eddy_nu[p] = 0.0;
    }
  }

  /*
   * grad(APR)
   */

  if (pd->gv[ACOUS_PREAL]) {
    v = ACOUS_PREAL;
    dofs = ei[upd->matrix_index[v]]->dof[v];

    for (p = 0; p < VIM; p++) {
      fv->grad_apr[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_apr[p] += *esp->apr[i] * bf[v]->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][ACOUS_PREAL] == -1) {
    for (p = 0; p < VIM; p++)
      fv->grad_apr[p] = 0.0;
  }

  if (pd->gv[ACOUS_PIMAG]) {
    v = ACOUS_PIMAG;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_api[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_api[p] += *esp->api[i] * bf[v]->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][ACOUS_PIMAG] == -1) {
    for (p = 0; p < VIM; p++)
      fv->grad_api[p] = 0.0;
  }

  if (pd->gv[ACOUS_REYN_STRESS]) {
    v = ACOUS_REYN_STRESS;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_ars[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_ars[p] += *esp->ars[i] * bf[v]->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][ACOUS_REYN_STRESS] == -1) {
    for (p = 0; p < VIM; p++)
      fv->grad_ars[p] = 0.0;
  }

  if (pd->gv[SHELL_BDYVELO]) {
    v = SHELL_BDYVELO;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_sh_bv[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_sh_bv[p] += *esp->sh_bv[i] * bf[v]->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][SHELL_BDYVELO] == -1) {
    for (p = 0; p < VIM; p++)
      fv->grad_sh_bv[p] = 0.0;
  }

  if (pd->gv[SHELL_LUBP]) {
    v = SHELL_LUBP;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_sh_p[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_sh_p[p] += *esp->sh_p[i] * bf[v]->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][SHELL_LUBP] == -1) {
    for (p = 0; p < VIM; p++)
      fv->grad_sh_p[p] = 0.0;
  }

  if (pd->gv[LUBP]) {
    v = LUBP;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_lubp[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_lubp[p] += *esp->lubp[i] * bf[v]->grad_phi[i][p];
        fv_old->grad_lubp[p] += *esp_old->lubp[i] * bf[v]->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][LUBP] == -1) {
    for (p = 0; p < VIM; p++)
      fv->grad_lubp[p] = 0.0;
  }

  if (pd->gv[LUBP_2]) {
    v = LUBP_2;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_lubp_2[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_lubp_2[p] += *esp->lubp_2[i] * bf[v]->grad_phi[i][p];
        fv_old->grad_lubp_2[p] += *esp_old->lubp_2[i] * bf[v]->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][LUBP_2] == -1) {
    for (p = 0; p < VIM; p++)
      fv->grad_lubp_2[p] = 0.0;
  }

  if (pd->gv[SHELL_LUB_CURV]) {
    v = SHELL_LUB_CURV;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_sh_l_curv[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_sh_l_curv[p] += *esp->sh_l_curv[i] * bf[v]->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][SHELL_LUB_CURV] == -1) {
    for (p = 0; p < VIM; p++)
      fv->grad_sh_l_curv[p] = 0.0;
  }

  if (pd->gv[SHELL_LUB_CURV_2]) {
    v = SHELL_LUB_CURV_2;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_sh_l_curv_2[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_sh_l_curv_2[p] += *esp->sh_l_curv_2[i] * bf[v]->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][SHELL_LUB_CURV_2] == -1) {
    for (p = 0; p < VIM; p++)
      fv->grad_sh_l_curv_2[p] = 0.0;
  }

  if (pd->gv[SHELL_PRESS_OPEN]) {
    v = SHELL_PRESS_OPEN;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_sh_p_open[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_sh_p_open[p] += *esp->sh_p_open[i] * bf[v]->grad_phi[i][p];
        fv_old->grad_sh_p_open[p] += *esp_old->sh_p_open[i] * bf[v]->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][SHELL_PRESS_OPEN] == -1) {
    for (p = 0; p < VIM; p++)
      fv->grad_sh_p_open[p] = 0.0;
  }

  if (pd->gv[SHELL_PRESS_OPEN_2]) {
    v = SHELL_PRESS_OPEN_2;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_sh_p_open_2[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_sh_p_open_2[p] += *esp->sh_p_open_2[i] * bf[v]->grad_phi[i][p];
        fv_old->grad_sh_p_open_2[p] += *esp_old->sh_p_open_2[i] * bf[v]->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][SHELL_PRESS_OPEN_2] == -1) {
    for (p = 0; p < VIM; p++)
      fv->grad_sh_p_open_2[p] = 0.0;
  }

  if (pd->gv[SHELL_SAT_1]) {
    v = SHELL_SAT_1;
    dofs = ei[pg->imtrx]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_sh_sat_1[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_sh_sat_1[p] += *esp->sh_sat_1[i] * bf[v]->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][SHELL_SAT_1] == -1) {
    for (p = 0; p < VIM; p++)
      fv->grad_sh_sat_1[p] = 0.0;
  }
  if (pd->gv[SHELL_FILMP]) {
    v = SHELL_FILMP;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_sh_fp[p] = 0.0;
      fv_old->grad_sh_fp[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_sh_fp[p] += *esp->sh_fp[i] * bf[v]->grad_phi[i][p];
        fv_old->grad_sh_fp[p] += *esp_old->sh_fp[i] * bf[v]->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][SHELL_FILMP] == -1) {
    for (p = 0; p < VIM; p++)
      fv->grad_sh_fp[p] = fv_old->grad_sh_fp[p] = 0.0;
  }

  if (pd->gv[SHELL_FILMH]) {
    v = SHELL_FILMH;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_sh_fh[p] = 0.0;
      fv_old->grad_sh_fh[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_sh_fh[p] += *esp->sh_fh[i] * bf[v]->grad_phi[i][p];
        fv_old->grad_sh_fh[p] += *esp_old->sh_fh[i] * bf[v]->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][SHELL_FILMH] == -1) {
    for (p = 0; p < VIM; p++)
      fv->grad_sh_fh[p] = fv_old->grad_sh_fh[p] = 0.0;
  }

  if (pd->gv[SHELL_PARTC]) {
    v = SHELL_PARTC;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_sh_pc[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_sh_pc[p] += *esp->sh_pc[i] * bf[v]->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][SHELL_PARTC] == -1) {
    for (p = 0; p < VIM; p++)
      fv->grad_sh_pc[p] = 0.0;
  }

  if (pd->gv[SHELL_TEMPERATURE]) {
    v = SHELL_TEMPERATURE;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_sh_t[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_sh_t[p] += *esp->sh_t[i] * bf[v]->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][SHELL_TEMPERATURE] == -1) {
    for (p = 0; p < VIM; p++)
      fv->grad_sh_t[p] = 0.0;
  }

  if (pd->gv[LIGHT_INTP]) {
    v = LIGHT_INTP;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_poynt[0][p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_poynt[0][p] += *esp->poynt[0][i] * bf[v]->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][LIGHT_INTP] == -1) {
    for (p = 0; p < VIM; p++)
      fv->grad_poynt[0][p] = 0.0;
  }

  if (pd->gv[LIGHT_INTM]) {
    v = LIGHT_INTM;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_poynt[1][p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_poynt[1][p] += *esp->poynt[1][i] * bf[v]->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][LIGHT_INTM] == -1) {
    for (p = 0; p < VIM; p++)
      fv->grad_poynt[1][p] = 0.0;
  }

  if (pd->gv[LIGHT_INTD]) {
    v = LIGHT_INTD;
    dofs = ei[upd->matrix_index[v]]->dof[v];
    for (p = 0; p < VIM; p++) {
      fv->grad_poynt[2][p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_poynt[2][p] += *esp->poynt[2][i] * bf[v]->grad_phi[i][p];
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][LIGHT_INTD] == -1) {
    for (p = 0; p < VIM; p++)
      fv->grad_poynt[2][p] = 0.0;
  }

  if (pd->gv[MOMENT0]) {
    for (r = 0; r < MAX_MOMENTS; r++) {
      v = MOMENT0 + r;
      dofs = ei[pg->imtrx]->dof[v];

      for (p = 0; p < VIM; p++) {
        fv->grad_moment[r][p] = 0.0;
        fv_old->grad_moment[r][p] = 0.0;

        for (i = 0; i < dofs; i++) {
          fv->grad_moment[r][p] += *esp->moment[r][i] * bf[v]->grad_phi[i][p];
          fv_old->grad_moment[r][p] += *esp_old->moment[r][i] * bf[v]->grad_phi[i][p];
        }
      }
    }
  }

  if (pd->gv[DENSITY_EQN]) {
    v = DENSITY_EQN;
    dofs = ei[pg->imtrx]->dof[v];
#ifdef DO_NO_UNROLL
    for (p = 0; p < VIM; p++) {
      fv->grad_rho[p] = 0.0;

      for (i = 0; i < dofs; i++) {
        fv->grad_rho[p] += *esp->rho[i] * bf[v]->grad_phi[i][p];
      }
    }
#else
    grad_scalar_fv_fill(esp->rho, bf[v]->grad_phi, dofs, fv->grad_rho);
#endif
  } else if (zero_unused_grads && upd->vp[pg->imtrx][DENSITY_EQN] == -1) {
    for (p = 0; p < VIM; p++)
      fv->grad_rho[p] = 0.0;
  }

  if (pd->gv[TFMP_PRES]) {
    v = TFMP_PRES;
    dofs = ei[pg->imtrx]->dof[v];
#ifdef DO_NO_UNROLL
    for (p = 0; p < VIM; p++) {
      fv->grad_tfmp_pres[p] = 0.0;

      for (i = 0; i < dofs; i++) {
        fv->grad_tfmp_pres[p] += *esp->tfmp_pres[i] * bf[v]->grad_phi[i][p];
        fv_old->grad_tfmp_pres[p] += *esp_old->tfmp_pres[i] * bf[v]->grad_phi[i][p];
      }
    }
#else
    grad_scalar_fv_fill(esp->tfmp_pres, bf[v]->grad_phi, dofs, fv->grad_tfmp_pres);
#endif
  } else if (zero_unused_grads && upd->vp[pg->imtrx][TFMP_PRES] == -1) {
    for (p = 0; p < VIM; p++)
      fv->grad_tfmp_pres[p] = 0.0;
  }

  if (pd->gv[TFMP_SAT]) {
    v = TFMP_SAT;
    dofs = ei[pg->imtrx]->dof[v];
#ifdef DO_NO_UNROLL
    for (p = 0; p < VIM; p++) {
      fv->grad_tfmp_sat[p] = 0.0;

      for (i = 0; i < dofs; i++) {
        fv->grad_tfmp_sat[p] += *esp->tfmp_sat[i] * bf[v]->grad_phi[i][p];
        fv_old->grad_tfmp_sat[p] += *esp_old->tfmp_sat[i] * bf[v]->grad_phi[i][p];
      }
    }
#else
    grad_scalar_fv_fill(esp->tfmp_sat, bf[v]->grad_phi, dofs, fv->grad_tfmp_sat);
#endif
  } else if (zero_unused_grads && upd->vp[pg->imtrx][TFMP_SAT] == -1) {
    for (p = 0; p < VIM; p++)
      fv->grad_tfmp_sat[p] = 0.0;
  }

  if (pd->gv[RESTIME]) {
    v = RESTIME;
    dofs = ei[pg->imtrx]->dof[v];
#ifdef DO_NO_UNROLL
    for (p = 0; p < VIM; p++) {
      fv->grad_restime[p] = 0.0;
      fv_old->grad_restime[p] = 0.0;
      for (i = 0; i < dofs; i++) {
        fv->grad_restime[p] += *esp->restime[i] * bf[v]->grad_phi[i][p];
        fv_old->grad_restime[p] += *esp_old->restime[i] * bf[v]->grad_phi[i][p];
      }
    }
#else
    grad_scalar_fv_fill(esp->restime, bf[v]->grad_phi, dofs, fv->grad_restime);
#endif
  } else if (zero_unused_grads && upd->vp[pg->imtrx][RESTIME] == -1) {
    for (p = 0; p < VIM; p++)
      fv->grad_restime[p] = 0.0;
  }

  /*
   * EM Wave Vector Gradients
   */
  if (pd->gv[EM_E1_REAL] || pd->gv[EM_E2_REAL] || pd->gv[EM_E3_REAL]) {
    v = EM_E1_REAL;
    if (bf[v]->interpolation != I_N1) {
      dofs = ei[upd->matrix_index[v]]->dof[v];

      // grad_vector_fv_fill(esp->em_er, bf[v]->grad_phi_e, dofs, fv->grad_em_er);
      for (p = 0; p < dim; p++) {
        for (q = 0; q < dim; q++) {
          fv->grad_em_er[p][q] = 0.0;
          fv_old->grad_em_er[p][q] = 0.0;
          for (r = 0; r < dim; r++) {
            for (i = 0; i < dofs; i++) {
              fv->grad_em_er[p][q] += (*esp->em_er[r][i]) * bf[v]->grad_phi_e[i][r][p][q];
              if (pd->TimeIntegration != STEADY) {
                fv_old->grad_em_er[p][q] += (*esp_old->em_er[r][i]) * bf[v]->grad_phi_e[i][r][p][q];
              }
            }
          }
        }
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][EM_E1_REAL] == -1 &&
             upd->vp[pg->imtrx][EM_E2_REAL] == -1 && upd->vp[pg->imtrx][EM_E3_REAL] == -1) {
    for (p = 0; p < DIM; p++) {
      for (q = 0; q < DIM; q++) {
        fv->grad_em_er[p][q] = 0.0;
      }
    }
  }

  if (pd->gv[EM_E1_IMAG] || pd->gv[EM_E2_IMAG] || pd->gv[EM_E3_IMAG]) {
    v = EM_E1_IMAG;
    if (bf[v]->interpolation != I_N1) {
      dofs = ei[upd->matrix_index[v]]->dof[v];

      // grad_vector_fv_fill(esp->em_ei, bf[v]->grad_phi_e, dofs, fv->grad_em_ei);
      for (p = 0; p < dim; p++) {
        for (q = 0; q < dim; q++) {
          fv->grad_em_ei[p][q] = 0.0;
          fv_old->grad_em_ei[p][q] = 0.0;
          for (r = 0; r < dim; r++) {
            for (i = 0; i < dofs; i++) {
              fv->grad_em_ei[p][q] += (*esp->em_ei[r][i]) * bf[v]->grad_phi_e[i][r][p][q];
              if (pd->TimeIntegration != STEADY) {
                fv_old->grad_em_ei[p][q] += (*esp_old->em_ei[r][i]) * bf[v]->grad_phi_e[i][r][p][q];
              }
            }
          }
        }
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][EM_E1_IMAG] == -1 &&
             upd->vp[pg->imtrx][EM_E2_IMAG] == -1 && upd->vp[pg->imtrx][EM_E3_IMAG] == -1) {
    for (p = 0; p < DIM; p++) {
      for (q = 0; q < DIM; q++) {
        fv->grad_em_ei[p][q] = 0.0;
      }
    }
  }
  if (pd->gv[EM_H1_REAL] || pd->gv[EM_H2_REAL] || pd->gv[EM_H3_REAL]) {
    v = EM_H1_REAL;
    dofs = ei[upd->matrix_index[v]]->dof[v];

    // grad_vector_fv_fill(esp->em_hr, bf[v]->grad_phi_e, dofs, fv->grad_em_hr);
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        fv->grad_em_hr[p][q] = 0.0;
        fv_old->grad_em_hr[p][q] = 0.0;
        for (r = 0; r < WIM; r++) {
          for (i = 0; i < dofs; i++) {
            fv->grad_em_hr[p][q] += (*esp->em_hr[r][i]) * bf[v]->grad_phi_e[i][r][p][q];
            if (pd->TimeIntegration != STEADY) {
              fv_old->grad_em_hr[p][q] += (*esp_old->em_hr[r][i]) * bf[v]->grad_phi_e[i][r][p][q];
            }
          }
        }
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][EM_H1_REAL] == -1 &&
             upd->vp[pg->imtrx][EM_H2_REAL] == -1 && upd->vp[pg->imtrx][EM_H3_REAL] == -1) {
    for (p = 0; p < DIM; p++) {
      for (q = 0; q < DIM; q++) {
        fv->grad_em_hr[p][q] = 0.0;
      }
    }
  }
  if (pd->gv[EM_H1_IMAG] || pd->gv[EM_H2_IMAG] || pd->gv[EM_H3_IMAG]) {
    v = EM_H1_IMAG;
    dofs = ei[upd->matrix_index[v]]->dof[v];

    // grad_vector_fv_fill(esp->em_hi, bf[v]->grad_phi_e, dofs, fv->grad_em_hi);
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        fv->grad_em_hi[p][q] = 0.0;
        fv_old->grad_em_hi[p][q] = 0.0;
        for (r = 0; r < WIM; r++) {
          for (i = 0; i < dofs; i++) {
            fv->grad_em_hi[p][q] += (*esp->em_hi[r][i]) * bf[v]->grad_phi_e[i][r][p][q];
            if (pd->TimeIntegration != STEADY) {
              fv_old->grad_em_hi[p][q] += (*esp_old->em_hi[r][i]) * bf[v]->grad_phi_e[i][r][p][q];
            }
          }
        }
      }
    }
  } else if (zero_unused_grads && upd->vp[pg->imtrx][EM_H1_IMAG] == -1 &&
             upd->vp[pg->imtrx][EM_H2_IMAG] == -1 && upd->vp[pg->imtrx][EM_H3_IMAG] == -1) {
    for (p = 0; p < DIM; p++) {
      for (q = 0; q < DIM; q++) {
        fv->grad_em_hi[p][q] = 0.0;
      }
    }
  }
  /*
   * External
   */
  if (efv->ev) {
    v = EXTERNAL;
    for (w = 0; w < efv->Num_external_field; w++) {
      dofs = ei[pg->imtrx]->dof_ext[w];
      if (efv->i[w] != I_TABLE) {
        if (strcmp(efv->name[w], "F1") == 0 && !pd->gv[PHASE1]) {
          /* load up the gradient of the the phase function variables
           * Currently, this is the only field whose gradient is computed when
           * is applied as externally
           */
          grad_scalar_fv_fill(evp->external_field[w], bfex[v]->grad_phi, dofs, fv->grad_pF[0]);
        } else {
          for (p = 0; p < VIM; p++) {
            fv->grad_ext_field[w][p] = 0.0;
            for (i = 0; i < dofs; i++) {
              fv->grad_ext_field[w][p] += *evp->external_field[w][i] * bfex[w]->grad_phi[i][p];
            }
          }
        }
      }
    }
  }

  zero_unused_grads = FALSE;
  return (status);
}

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

int load_fv_mesh_derivs(int okToZero)

/*******************************************************************************
 *
 * load_fv_mesh_derivs() -- ld mesh derivatives of field var grads at gauss pt
 *
 * input: (assume the appropriate parts of esp, bf, ei, pd are filled in)
 * ------
 *
 *
 * output: ( this routine fills in the following parts of the fv structure)
 * -------
 *	fv->
 *  		d_grad_T_dmesh[p] [b][j] 	d( e_p . grad(T) )
 * 						------------------
 * 						    d ( d_b,j )
 *
 *		d_grad_P_dmesh[p] [b][j] 	d( e_p . grad(P) )
 * 						------------------
 * 						    d ( d_b,j )
 *
 *		d_grad_c_dmesh[p][w] [b][j] 	d( e_p . grad(c_w) )
 * 						------------------
 * 						    d ( d_b,j )
 *
 *		d_grad_p_liq_dmesh[p] [b][j] 	d( e_p . grad(p_liq) )
 * 						------------------
 * 						    d ( d_b,j )
 *
 *		d_div_v_dmesh [b][j] 		d( grad.v )
 * 						-----------
 * 						d ( d_b,j )
 *
 *		d_grad_v_dmesh[p][q] [b][j] 	d( e_p e_q : grad(v) )
 * 						------------------
 * 						    d ( d_b,j )
 *
 *		d_grad_X_dmesh[p][q] [b][j] 	d( e_p e_q : grad(X) )
 * 						------------------
 * 						    d ( d_b,j )
 *
 *		d_grad_pv_dmesh[p][q] [b][j] 	d( e_p e_q : grad(pv) )
 * 						------------------
 * 						    d ( d_b,j )
 *
 *		d_div_d_dmesh [b][j] 		d( grad.d )
 * 						-----------
 * 						d ( d_b,j )
 *
 *		d_grad_d_dmesh[p][q] [b][j] 	d( e_p e_q : grad(d) )
 * 						------------------
 * 						    d ( d_b,j )
 *
 *		d_grad_S_dmesh[p][q][r] [b][j] 	d( e_p e_q e_r : grad(S) )
 * 						------------------
 * 						    d ( d_b,j )
 *
 *		d_grad_G_dmesh[p][q][r] [b][j] 	d( e_p e_q e_r : grad(G) )
 * 						------------------
 * 						    d ( d_b,j )
 *
 *		d_grad_Gt_dmesh[p][q][r] [b][j] d( e_p e_q e_r : grad(G)^t )
 * 						------------------
 * 						    d ( d_b,j )
 *
 *		d_grad_vd_dmesh[p][q] [b][j] 	d( e_p e_q : grad(vd) )
 * 						------------------
 * 						    d ( d_b,j )
 *
 * NOTE: To accommodate shell elements, this function has been modified
 *       so that fv variables are not zeroed out when they are active
 *       on an element block other than the current one.
 *       The check done for variable v is then:
 *          if ( pd->v[pg->imtrx][v] || upd->vp[pg->imtrx][v] == -1 )
 *       In many cases below, this conditional zeroing is done in
 *       a separate small loop before the main one.
 *
 *
 * Return values:
 *		0 -- if things went OK
 *	       -1 -- if a problem occurred
 *
 *
 * Notes:
 * ------
 *		[1] There is an extra Langrangian contribution for the mesh
 *		    displacment unknowns, in addition to the parts that the
 *		    have to do with spatial gradient operators over a
 *		    deforming mesh that all variables have.
 *
 *
 * Created:	Fri Mar 18 10:18:03 MST 1994 pasacki@sandia.gov
 *
 * Modified:	Tue Feb 21 11:16 MST 1995 pasacki@sandia.gov
 ****************************************************************************/
{
  int v;    /* variable type indicator */
  int i, j; /* index */
  int b;
  int p, q, r, s; /* dimension index */
  int vdofs = 0;  /* degrees of freedom for var  in the elem */
  int mdofs;      /* degrees of freedom for mesh in the elem */
  int w;          /* concentration species */
  int status;
  int dimNonSym; /* # of dimensions that don't have symmetry */
  int dim;       /* # dimensions in the physical mesh */
  unsigned int siz;
  int mode; /* modal counter */
  int transient_run, discontinuous_porous_media, ve_solid;

  dbl T_i, v_ri, P_i, d_ri, F_i, d_dot_ri;
  dbl em_er_ri, em_ei_ri, em_hr_ri, em_hi_ri;

  struct Basis_Functions *bfm, *bfv; /* For mesh variables. */

  static int is_initialized = FALSE;
  int VIMis3;

  status = 0;

  VIMis3 = (VIM == 3) ? TRUE : FALSE;

  /*
   * If this is not a deforming mesh problem, then
   * leave now.
   */
  int imtrx;
  int mesh_on = FALSE;
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    if (ei[imtrx]->deforming_mesh) {
      mesh_on = TRUE;
    }
  }
  if (!mesh_on) {
    return status;
  }

  dim = pd->Num_Dim;
  dimNonSym = pd->Num_Dim;

  mdofs = ei[upd->matrix_index[R_MESH1]]->dof[R_MESH1];

  bfm = bf[R_MESH1];

  transient_run = pd->TimeIntegration != STEADY;
  discontinuous_porous_media = mp->PorousMediaType != CONTINUOUS;
  ve_solid = cr->MeshFluxModel == KELVIN_VOIGT;

  okToZero = FALSE;

  /*
   * d(grad(T))/dmesh
   */

  if (pd->gv[TEMPERATURE]) {
    v = TEMPERATURE;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];
#ifdef DO_NOT_UNROLL
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_T_dmesh[0][0][0]), 0, siz);
    for (i = 0; i < vdofs; i++) {
      T_i = *esp->T[i];
      for (p = 0; p < dimNonSym; p++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_T_dmesh[p][b][j] += T_i * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
#else
    for (j = 0; j < mdofs; j++) {
      T_i = *esp->T[0];

      fv->d_grad_T_dmesh[0][0][j] = T_i * bfv->d_grad_phi_dmesh[0][0][0][j];
      fv->d_grad_T_dmesh[1][1][j] = T_i * bfv->d_grad_phi_dmesh[0][1][1][j];
      fv->d_grad_T_dmesh[1][0][j] = T_i * bfv->d_grad_phi_dmesh[0][1][0][j];
      fv->d_grad_T_dmesh[0][1][j] = T_i * bfv->d_grad_phi_dmesh[0][0][1][j];

      if (dimNonSym == 3) {
        fv->d_grad_T_dmesh[2][2][j] = T_i * bfv->d_grad_phi_dmesh[0][2][2][j];
        fv->d_grad_T_dmesh[2][0][j] = T_i * bfv->d_grad_phi_dmesh[0][2][0][j];
        fv->d_grad_T_dmesh[2][1][j] = T_i * bfv->d_grad_phi_dmesh[0][2][1][j];
        fv->d_grad_T_dmesh[0][2][j] = T_i * bfv->d_grad_phi_dmesh[0][0][2][j];
        fv->d_grad_T_dmesh[1][2][j] = T_i * bfv->d_grad_phi_dmesh[0][1][2][j];
      }

      for (i = 1; i < vdofs; i++) {
        T_i = *esp->T[i];

        fv->d_grad_T_dmesh[0][0][j] += T_i * bfv->d_grad_phi_dmesh[i][0][0][j];
        fv->d_grad_T_dmesh[1][1][j] += T_i * bfv->d_grad_phi_dmesh[i][1][1][j];
        fv->d_grad_T_dmesh[1][0][j] += T_i * bfv->d_grad_phi_dmesh[i][1][0][j];
        fv->d_grad_T_dmesh[0][1][j] += T_i * bfv->d_grad_phi_dmesh[i][0][1][j];

        if (dimNonSym == 3) {
          fv->d_grad_T_dmesh[2][2][j] += T_i * bfv->d_grad_phi_dmesh[i][2][2][j];
          fv->d_grad_T_dmesh[2][0][j] += T_i * bfv->d_grad_phi_dmesh[i][2][0][j];
          fv->d_grad_T_dmesh[2][1][j] += T_i * bfv->d_grad_phi_dmesh[i][2][1][j];
          fv->d_grad_T_dmesh[0][2][j] += T_i * bfv->d_grad_phi_dmesh[i][0][2][j];
          fv->d_grad_T_dmesh[1][2][j] += T_i * bfv->d_grad_phi_dmesh[i][1][2][j];
        }
      }
    }
#endif
  } else if (upd->vp[pg->imtrx][TEMPERATURE] != -1 && !is_initialized && okToZero) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_T_dmesh[0][0][0]), 0, siz);
  }

  int mom;
  for (mom = 0; mom < MAX_MOMENTS; mom++) {
    v = MOMENT0 + mom;
    if (pd->gv[v]) {
      bfv = bf[v];
      vdofs = ei[pg->imtrx]->dof[v];

      siz = sizeof(double) * DIM * DIM * MDE;
      memset(&(fv->d_grad_moment_dmesh[mom][0][0][0]), 0, siz);
      for (i = 0; i < vdofs; i++) {
        for (p = 0; p < dimNonSym; p++) {
          for (b = 0; b < dim; b++) {
            for (j = 0; j < mdofs; j++) {
              fv->d_grad_moment_dmesh[mom][p][b][j] +=
                  (*esp->moment[mom][i]) * bfv->d_grad_phi_dmesh[i][p][b][j];
            }
          }
        }
      }

    } else if (upd->vp[pg->imtrx][v] != -1 && !is_initialized && okToZero) {
      siz = sizeof(double) * DIM * DIM * MDE;
      memset(&(fv->d_grad_moment_dmesh[mom][0][0][0]), 0, siz);
    }
  }

  /*
   * d(grad(V))/dmesh
   */
  if (pd->gv[VOLTAGE]) {
    v = VOLTAGE;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];
#ifdef DO_NOT_UNROLL
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_V_dmesh[0][0][0]), 0, siz);
    for (p = 0; p < dimNonSym; p++) {
      for (b = 0; b < dim; b++) {
        for (j = 0; j < mdofs; j++) {
          for (i = 0; i < vdofs; i++) {
            fv->d_grad_V_dmesh[p][b][j] += *esp->V[i] * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
#else
    for (j = 0; j < mdofs; j++) {
      fv->d_grad_V_dmesh[0][0][j] = *esp->V[0] * bfv->d_grad_phi_dmesh[0][0][0][j];
      fv->d_grad_V_dmesh[1][1][j] = *esp->V[0] * bfv->d_grad_phi_dmesh[0][1][1][j];
      fv->d_grad_V_dmesh[1][0][j] = *esp->V[0] * bfv->d_grad_phi_dmesh[0][1][0][j];
      fv->d_grad_V_dmesh[0][1][j] = *esp->V[0] * bfv->d_grad_phi_dmesh[0][0][1][j];

      if (dimNonSym == 3) {
        fv->d_grad_V_dmesh[2][2][j] = *esp->V[0] * bfv->d_grad_phi_dmesh[0][2][2][j];
        fv->d_grad_V_dmesh[2][0][j] = *esp->V[0] * bfv->d_grad_phi_dmesh[0][2][0][j];
        fv->d_grad_V_dmesh[2][1][j] = *esp->V[0] * bfv->d_grad_phi_dmesh[0][2][1][j];
        fv->d_grad_V_dmesh[0][2][j] = *esp->V[0] * bfv->d_grad_phi_dmesh[0][0][2][j];
        fv->d_grad_V_dmesh[1][2][j] = *esp->V[0] * bfv->d_grad_phi_dmesh[0][1][2][j];
      }

      for (i = 1; i < vdofs; i++) {
        fv->d_grad_V_dmesh[0][0][j] += *esp->V[i] * bfv->d_grad_phi_dmesh[i][0][0][j];
        fv->d_grad_V_dmesh[1][1][j] += *esp->V[i] * bfv->d_grad_phi_dmesh[i][1][1][j];
        fv->d_grad_V_dmesh[1][0][j] += *esp->V[i] * bfv->d_grad_phi_dmesh[i][1][0][j];
        fv->d_grad_V_dmesh[0][1][j] += *esp->V[i] * bfv->d_grad_phi_dmesh[i][0][1][j];

        if (dimNonSym == 3) {
          fv->d_grad_V_dmesh[2][2][j] += *esp->V[i] * bfv->d_grad_phi_dmesh[i][2][2][j];
          fv->d_grad_V_dmesh[2][0][j] += *esp->V[i] * bfv->d_grad_phi_dmesh[i][2][0][j];
          fv->d_grad_V_dmesh[2][1][j] += *esp->V[i] * bfv->d_grad_phi_dmesh[i][2][1][j];
          fv->d_grad_V_dmesh[0][2][j] += *esp->V[i] * bfv->d_grad_phi_dmesh[i][0][2][j];
          fv->d_grad_V_dmesh[1][2][j] += *esp->V[i] * bfv->d_grad_phi_dmesh[i][1][2][j];
        }
      }
    }
#endif
  } else if (upd->vp[pg->imtrx][VOLTAGE] != -1 && !is_initialized && okToZero) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_V_dmesh[0][0][0]), 0, siz);
  }

  /*
   * d(grad(sh_K))/dmesh
   */

  if (pd->v[pg->imtrx][SHELL_CURVATURE]) {
    v = SHELL_CURVATURE;
    bfv = bf[v];
    vdofs = ei[pg->imtrx]->dof[v];
#ifdef DO_NOT_UNROLL
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_sh_K_dmesh[0][0][0]), 0, siz);
    for (i = 0; i < vdofs; i++) {
      T_i = *esp->sh_K[i];
      for (p = 0; p < dimNonSym; p++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_sh_K_dmesh[p][b][j] += T_i * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
#else
    for (j = 0; j < mdofs; j++) {
      T_i = *esp->sh_K[0];

      fv->d_grad_sh_K_dmesh[0][0][j] = T_i * bfv->d_grad_phi_dmesh[0][0][0][j];
      fv->d_grad_sh_K_dmesh[1][1][j] = T_i * bfv->d_grad_phi_dmesh[0][1][1][j];
      fv->d_grad_sh_K_dmesh[1][0][j] = T_i * bfv->d_grad_phi_dmesh[0][1][0][j];
      fv->d_grad_sh_K_dmesh[0][1][j] = T_i * bfv->d_grad_phi_dmesh[0][0][1][j];

      if (dimNonSym == 3) {
        fv->d_grad_sh_K_dmesh[2][2][j] = T_i * bfv->d_grad_phi_dmesh[0][2][2][j];
        fv->d_grad_sh_K_dmesh[2][0][j] = T_i * bfv->d_grad_phi_dmesh[0][2][0][j];
        fv->d_grad_sh_K_dmesh[2][1][j] = T_i * bfv->d_grad_phi_dmesh[0][2][1][j];
        fv->d_grad_sh_K_dmesh[0][2][j] = T_i * bfv->d_grad_phi_dmesh[0][0][2][j];
        fv->d_grad_sh_K_dmesh[1][2][j] = T_i * bfv->d_grad_phi_dmesh[0][1][2][j];
      }

      for (i = 1; i < vdofs; i++) {
        T_i = *esp->sh_K[i];

        fv->d_grad_sh_K_dmesh[0][0][j] += T_i * bfv->d_grad_phi_dmesh[i][0][0][j];
        fv->d_grad_sh_K_dmesh[1][1][j] += T_i * bfv->d_grad_phi_dmesh[i][1][1][j];
        fv->d_grad_sh_K_dmesh[1][0][j] += T_i * bfv->d_grad_phi_dmesh[i][1][0][j];
        fv->d_grad_sh_K_dmesh[0][1][j] += T_i * bfv->d_grad_phi_dmesh[i][0][1][j];

        if (dimNonSym == 3) {
          fv->d_grad_sh_K_dmesh[2][2][j] += T_i * bfv->d_grad_phi_dmesh[i][2][2][j];
          fv->d_grad_sh_K_dmesh[2][0][j] += T_i * bfv->d_grad_phi_dmesh[i][2][0][j];
          fv->d_grad_sh_K_dmesh[2][1][j] += T_i * bfv->d_grad_phi_dmesh[i][2][1][j];
          fv->d_grad_sh_K_dmesh[0][2][j] += T_i * bfv->d_grad_phi_dmesh[i][0][2][j];
          fv->d_grad_sh_K_dmesh[1][2][j] += T_i * bfv->d_grad_phi_dmesh[i][1][2][j];
        }
      }
    }
#endif
  } else if (upd->vp[pg->imtrx][SHELL_CURVATURE] != -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_sh_K_dmesh[0][0][0]), 0, siz);
  }

  if (pd->v[pg->imtrx][SHELL_CURVATURE2]) {
    v = SHELL_CURVATURE2;
    bfv = bf[v];
    vdofs = ei[pg->imtrx]->dof[v];
#ifdef DO_NOT_UNROLL
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_sh_K2_dmesh[0][0][0]), 0, siz);
    for (i = 0; i < vdofs; i++) {
      T_i = *esp->sh_K2[i];
      for (p = 0; p < dimNonSym; p++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_sh_K2_dmesh[p][b][j] += T_i * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
#else
    for (j = 0; j < mdofs; j++) {
      T_i = *esp->sh_K2[0];

      fv->d_grad_sh_K2_dmesh[0][0][j] = T_i * bfv->d_grad_phi_dmesh[0][0][0][j];
      fv->d_grad_sh_K2_dmesh[1][1][j] = T_i * bfv->d_grad_phi_dmesh[0][1][1][j];
      fv->d_grad_sh_K2_dmesh[1][0][j] = T_i * bfv->d_grad_phi_dmesh[0][1][0][j];
      fv->d_grad_sh_K2_dmesh[0][1][j] = T_i * bfv->d_grad_phi_dmesh[0][0][1][j];

      if (dimNonSym == 3) {
        fv->d_grad_sh_K2_dmesh[2][2][j] = T_i * bfv->d_grad_phi_dmesh[0][2][2][j];
        fv->d_grad_sh_K2_dmesh[2][0][j] = T_i * bfv->d_grad_phi_dmesh[0][2][0][j];
        fv->d_grad_sh_K2_dmesh[2][1][j] = T_i * bfv->d_grad_phi_dmesh[0][2][1][j];
        fv->d_grad_sh_K2_dmesh[0][2][j] = T_i * bfv->d_grad_phi_dmesh[0][0][2][j];
        fv->d_grad_sh_K2_dmesh[1][2][j] = T_i * bfv->d_grad_phi_dmesh[0][1][2][j];
      }

      for (i = 1; i < vdofs; i++) {
        T_i = *esp->sh_K2[i];

        fv->d_grad_sh_K2_dmesh[0][0][j] += T_i * bfv->d_grad_phi_dmesh[i][0][0][j];
        fv->d_grad_sh_K2_dmesh[1][1][j] += T_i * bfv->d_grad_phi_dmesh[i][1][1][j];
        fv->d_grad_sh_K2_dmesh[1][0][j] += T_i * bfv->d_grad_phi_dmesh[i][1][0][j];
        fv->d_grad_sh_K2_dmesh[0][1][j] += T_i * bfv->d_grad_phi_dmesh[i][0][1][j];

        if (dimNonSym == 3) {
          fv->d_grad_sh_K2_dmesh[2][2][j] += T_i * bfv->d_grad_phi_dmesh[i][2][2][j];
          fv->d_grad_sh_K2_dmesh[2][0][j] += T_i * bfv->d_grad_phi_dmesh[i][2][0][j];
          fv->d_grad_sh_K2_dmesh[2][1][j] += T_i * bfv->d_grad_phi_dmesh[i][2][1][j];
          fv->d_grad_sh_K2_dmesh[0][2][j] += T_i * bfv->d_grad_phi_dmesh[i][0][2][j];
          fv->d_grad_sh_K2_dmesh[1][2][j] += T_i * bfv->d_grad_phi_dmesh[i][1][2][j];
        }
      }
    }
#endif
  } else if (upd->vp[pg->imtrx][SHELL_CURVATURE2] != -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_sh_K2_dmesh[0][0][0]), 0, siz);
  }

  /*
   * d(grad(apr))/dmesh
   */
  if (pd->gv[ACOUS_PREAL]) {
    v = ACOUS_PREAL;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];
#ifdef DO_NOT_UNROLL
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_apr_dmesh[0][0][0]), 0, siz);
    for (i = 0; i < vdofs; i++) {
      T_i = *esp->apr[i];
      for (p = 0; p < dimNonSym; p++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_apr_dmesh[p][b][j] += T_i * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
#else
    for (j = 0; j < mdofs; j++) {
      T_i = *esp->apr[0];

      fv->d_grad_apr_dmesh[0][0][j] = T_i * bfv->d_grad_phi_dmesh[0][0][0][j];
      fv->d_grad_apr_dmesh[1][1][j] = T_i * bfv->d_grad_phi_dmesh[0][1][1][j];
      fv->d_grad_apr_dmesh[1][0][j] = T_i * bfv->d_grad_phi_dmesh[0][1][0][j];
      fv->d_grad_apr_dmesh[0][1][j] = T_i * bfv->d_grad_phi_dmesh[0][0][1][j];

      if (dimNonSym == 3) {
        fv->d_grad_apr_dmesh[2][2][j] = T_i * bfv->d_grad_phi_dmesh[0][2][2][j];
        fv->d_grad_apr_dmesh[2][0][j] = T_i * bfv->d_grad_phi_dmesh[0][2][0][j];
        fv->d_grad_apr_dmesh[2][1][j] = T_i * bfv->d_grad_phi_dmesh[0][2][1][j];
        fv->d_grad_apr_dmesh[0][2][j] = T_i * bfv->d_grad_phi_dmesh[0][0][2][j];
        fv->d_grad_apr_dmesh[1][2][j] = T_i * bfv->d_grad_phi_dmesh[0][1][2][j];
      }

      for (i = 1; i < vdofs; i++) {
        T_i = *esp->apr[i];

        fv->d_grad_apr_dmesh[0][0][j] += T_i * bfv->d_grad_phi_dmesh[i][0][0][j];
        fv->d_grad_apr_dmesh[1][1][j] += T_i * bfv->d_grad_phi_dmesh[i][1][1][j];
        fv->d_grad_apr_dmesh[1][0][j] += T_i * bfv->d_grad_phi_dmesh[i][1][0][j];
        fv->d_grad_apr_dmesh[0][1][j] += T_i * bfv->d_grad_phi_dmesh[i][0][1][j];

        if (dimNonSym == 3) {
          fv->d_grad_apr_dmesh[2][2][j] += T_i * bfv->d_grad_phi_dmesh[i][2][2][j];
          fv->d_grad_apr_dmesh[2][0][j] += T_i * bfv->d_grad_phi_dmesh[i][2][0][j];
          fv->d_grad_apr_dmesh[2][1][j] += T_i * bfv->d_grad_phi_dmesh[i][2][1][j];
          fv->d_grad_apr_dmesh[0][2][j] += T_i * bfv->d_grad_phi_dmesh[i][0][2][j];
          fv->d_grad_apr_dmesh[1][2][j] += T_i * bfv->d_grad_phi_dmesh[i][1][2][j];
        }
      }
    }
#endif
  } else if (upd->vp[pg->imtrx][ACOUS_PREAL] != -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_apr_dmesh[0][0][0]), 0, siz);
  }

  if (pd->gv[ACOUS_PIMAG]) {
    v = ACOUS_PIMAG;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];
#ifdef DO_NOT_UNROLL
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_api_dmesh[0][0][0]), 0, siz);
    for (i = 0; i < vdofs; i++) {
      T_i = *esp->api[i];
      for (p = 0; p < dimNonSym; p++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_api_dmesh[p][b][j] += T_i * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
#else
    for (j = 0; j < mdofs; j++) {
      T_i = *esp->api[0];

      fv->d_grad_api_dmesh[0][0][j] = T_i * bfv->d_grad_phi_dmesh[0][0][0][j];
      fv->d_grad_api_dmesh[1][1][j] = T_i * bfv->d_grad_phi_dmesh[0][1][1][j];
      fv->d_grad_api_dmesh[1][0][j] = T_i * bfv->d_grad_phi_dmesh[0][1][0][j];
      fv->d_grad_api_dmesh[0][1][j] = T_i * bfv->d_grad_phi_dmesh[0][0][1][j];

      if (dimNonSym == 3) {
        fv->d_grad_api_dmesh[2][2][j] = T_i * bfv->d_grad_phi_dmesh[0][2][2][j];
        fv->d_grad_api_dmesh[2][0][j] = T_i * bfv->d_grad_phi_dmesh[0][2][0][j];
        fv->d_grad_api_dmesh[2][1][j] = T_i * bfv->d_grad_phi_dmesh[0][2][1][j];
        fv->d_grad_api_dmesh[0][2][j] = T_i * bfv->d_grad_phi_dmesh[0][0][2][j];
        fv->d_grad_api_dmesh[1][2][j] = T_i * bfv->d_grad_phi_dmesh[0][1][2][j];
      }

      for (i = 1; i < vdofs; i++) {
        T_i = *esp->api[i];

        fv->d_grad_api_dmesh[0][0][j] += T_i * bfv->d_grad_phi_dmesh[i][0][0][j];
        fv->d_grad_api_dmesh[1][1][j] += T_i * bfv->d_grad_phi_dmesh[i][1][1][j];
        fv->d_grad_api_dmesh[1][0][j] += T_i * bfv->d_grad_phi_dmesh[i][1][0][j];
        fv->d_grad_api_dmesh[0][1][j] += T_i * bfv->d_grad_phi_dmesh[i][0][1][j];

        if (dimNonSym == 3) {
          fv->d_grad_api_dmesh[2][2][j] += T_i * bfv->d_grad_phi_dmesh[i][2][2][j];
          fv->d_grad_api_dmesh[2][0][j] += T_i * bfv->d_grad_phi_dmesh[i][2][0][j];
          fv->d_grad_api_dmesh[2][1][j] += T_i * bfv->d_grad_phi_dmesh[i][2][1][j];
          fv->d_grad_api_dmesh[0][2][j] += T_i * bfv->d_grad_phi_dmesh[i][0][2][j];
          fv->d_grad_api_dmesh[1][2][j] += T_i * bfv->d_grad_phi_dmesh[i][1][2][j];
        }
      }
    }
#endif
  } else if (upd->vp[pg->imtrx][ACOUS_PIMAG] != -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_api_dmesh[0][0][0]), 0, siz);
  }

  if (pd->gv[ACOUS_REYN_STRESS]) {
    v = ACOUS_REYN_STRESS;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];
#ifdef DO_NOT_UNROLL
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_ars_dmesh[0][0][0]), 0, siz);
    for (i = 0; i < vdofs; i++) {
      T_i = *esp->ars[i];
      for (p = 0; p < dimNonSym; p++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_ars_dmesh[p][b][j] += T_i * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
#else
    for (j = 0; j < mdofs; j++) {
      T_i = *esp->ars[0];

      fv->d_grad_ars_dmesh[0][0][j] = T_i * bfv->d_grad_phi_dmesh[0][0][0][j];
      fv->d_grad_ars_dmesh[1][1][j] = T_i * bfv->d_grad_phi_dmesh[0][1][1][j];
      fv->d_grad_ars_dmesh[1][0][j] = T_i * bfv->d_grad_phi_dmesh[0][1][0][j];
      fv->d_grad_ars_dmesh[0][1][j] = T_i * bfv->d_grad_phi_dmesh[0][0][1][j];

      if (dimNonSym == 3) {
        fv->d_grad_ars_dmesh[2][2][j] = T_i * bfv->d_grad_phi_dmesh[0][2][2][j];
        fv->d_grad_ars_dmesh[2][0][j] = T_i * bfv->d_grad_phi_dmesh[0][2][0][j];
        fv->d_grad_ars_dmesh[2][1][j] = T_i * bfv->d_grad_phi_dmesh[0][2][1][j];
        fv->d_grad_ars_dmesh[0][2][j] = T_i * bfv->d_grad_phi_dmesh[0][0][2][j];
        fv->d_grad_ars_dmesh[1][2][j] = T_i * bfv->d_grad_phi_dmesh[0][1][2][j];
      }

      for (i = 1; i < vdofs; i++) {
        T_i = *esp->ars[i];

        fv->d_grad_ars_dmesh[0][0][j] += T_i * bfv->d_grad_phi_dmesh[i][0][0][j];
        fv->d_grad_ars_dmesh[1][1][j] += T_i * bfv->d_grad_phi_dmesh[i][1][1][j];
        fv->d_grad_ars_dmesh[1][0][j] += T_i * bfv->d_grad_phi_dmesh[i][1][0][j];
        fv->d_grad_ars_dmesh[0][1][j] += T_i * bfv->d_grad_phi_dmesh[i][0][1][j];

        if (dimNonSym == 3) {
          fv->d_grad_ars_dmesh[2][2][j] += T_i * bfv->d_grad_phi_dmesh[i][2][2][j];
          fv->d_grad_ars_dmesh[2][0][j] += T_i * bfv->d_grad_phi_dmesh[i][2][0][j];
          fv->d_grad_ars_dmesh[2][1][j] += T_i * bfv->d_grad_phi_dmesh[i][2][1][j];
          fv->d_grad_ars_dmesh[0][2][j] += T_i * bfv->d_grad_phi_dmesh[i][0][2][j];
          fv->d_grad_ars_dmesh[1][2][j] += T_i * bfv->d_grad_phi_dmesh[i][1][2][j];
        }
      }
    }
#endif
  } else if (upd->vp[pg->imtrx][ACOUS_REYN_STRESS] != -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_ars_dmesh[0][0][0]), 0, siz);
  }

  if (pd->gv[SHELL_BDYVELO]) {
    v = SHELL_BDYVELO;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];
#ifdef DO_NOT_UNROLL
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_sh_bv_dmesh[0][0][0]), 0, siz);
    for (i = 0; i < vdofs; i++) {
      T_i = *esp->sh_bv[i];
      for (p = 0; p < dimNonSym; p++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_sh_bv_dmesh[p][b][j] += T_i * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
#else
    for (j = 0; j < mdofs; j++) {
      T_i = *esp->sh_bv[0];

      fv->d_grad_sh_bv_dmesh[0][0][j] = T_i * bfv->d_grad_phi_dmesh[0][0][0][j];
      fv->d_grad_sh_bv_dmesh[1][1][j] = T_i * bfv->d_grad_phi_dmesh[0][1][1][j];
      fv->d_grad_sh_bv_dmesh[1][0][j] = T_i * bfv->d_grad_phi_dmesh[0][1][0][j];
      fv->d_grad_sh_bv_dmesh[0][1][j] = T_i * bfv->d_grad_phi_dmesh[0][0][1][j];

      if (dimNonSym == 3) {
        fv->d_grad_sh_bv_dmesh[2][2][j] = T_i * bfv->d_grad_phi_dmesh[0][2][2][j];
        fv->d_grad_sh_bv_dmesh[2][0][j] = T_i * bfv->d_grad_phi_dmesh[0][2][0][j];
        fv->d_grad_sh_bv_dmesh[2][1][j] = T_i * bfv->d_grad_phi_dmesh[0][2][1][j];
        fv->d_grad_sh_bv_dmesh[0][2][j] = T_i * bfv->d_grad_phi_dmesh[0][0][2][j];
        fv->d_grad_sh_bv_dmesh[1][2][j] = T_i * bfv->d_grad_phi_dmesh[0][1][2][j];
      }

      for (i = 1; i < vdofs; i++) {
        T_i = *esp->sh_bv[i];

        fv->d_grad_sh_bv_dmesh[0][0][j] += T_i * bfv->d_grad_phi_dmesh[i][0][0][j];
        fv->d_grad_sh_bv_dmesh[1][1][j] += T_i * bfv->d_grad_phi_dmesh[i][1][1][j];
        fv->d_grad_sh_bv_dmesh[1][0][j] += T_i * bfv->d_grad_phi_dmesh[i][1][0][j];
        fv->d_grad_sh_bv_dmesh[0][1][j] += T_i * bfv->d_grad_phi_dmesh[i][0][1][j];

        if (dimNonSym == 3) {
          fv->d_grad_sh_bv_dmesh[2][2][j] += T_i * bfv->d_grad_phi_dmesh[i][2][2][j];
          fv->d_grad_sh_bv_dmesh[2][0][j] += T_i * bfv->d_grad_phi_dmesh[i][2][0][j];
          fv->d_grad_sh_bv_dmesh[2][1][j] += T_i * bfv->d_grad_phi_dmesh[i][2][1][j];
          fv->d_grad_sh_bv_dmesh[0][2][j] += T_i * bfv->d_grad_phi_dmesh[i][0][2][j];
          fv->d_grad_sh_bv_dmesh[1][2][j] += T_i * bfv->d_grad_phi_dmesh[i][1][2][j];
        }
      }
    }
#endif
  } else if (upd->vp[pg->imtrx][SHELL_BDYVELO] != -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_sh_bv_dmesh[0][0][0]), 0, siz);
  }

  if (pd->gv[LUBP]) {
    v = LUBP;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];
#ifdef DO_NOT_UNROLL
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_lubp_dmesh[0][0][0]), 0, siz);
    for (i = 0; i < vdofs; i++) {
      T_i = *esp->lubp[i];
      for (p = 0; p < dimNonSym; p++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_lubp_dmesh[p][b][j] += T_i * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
#else
    for (j = 0; j < mdofs; j++) {
      T_i = *esp->lubp[0];

      fv->d_grad_lubp_dmesh[0][0][j] = T_i * bfv->d_grad_phi_dmesh[0][0][0][j];
      fv->d_grad_lubp_dmesh[1][1][j] = T_i * bfv->d_grad_phi_dmesh[0][1][1][j];
      fv->d_grad_lubp_dmesh[1][0][j] = T_i * bfv->d_grad_phi_dmesh[0][1][0][j];
      fv->d_grad_lubp_dmesh[0][1][j] = T_i * bfv->d_grad_phi_dmesh[0][0][1][j];

      if (dimNonSym == 3) {
        fv->d_grad_lubp_dmesh[2][2][j] = T_i * bfv->d_grad_phi_dmesh[0][2][2][j];
        fv->d_grad_lubp_dmesh[2][0][j] = T_i * bfv->d_grad_phi_dmesh[0][2][0][j];
        fv->d_grad_lubp_dmesh[2][1][j] = T_i * bfv->d_grad_phi_dmesh[0][2][1][j];
        fv->d_grad_lubp_dmesh[0][2][j] = T_i * bfv->d_grad_phi_dmesh[0][0][2][j];
        fv->d_grad_lubp_dmesh[1][2][j] = T_i * bfv->d_grad_phi_dmesh[0][1][2][j];
      }

      for (i = 1; i < vdofs; i++) {
        T_i = *esp->lubp[i];

        fv->d_grad_lubp_dmesh[0][0][j] += T_i * bfv->d_grad_phi_dmesh[i][0][0][j];
        fv->d_grad_lubp_dmesh[1][1][j] += T_i * bfv->d_grad_phi_dmesh[i][1][1][j];
        fv->d_grad_lubp_dmesh[1][0][j] += T_i * bfv->d_grad_phi_dmesh[i][1][0][j];
        fv->d_grad_lubp_dmesh[0][1][j] += T_i * bfv->d_grad_phi_dmesh[i][0][1][j];

        if (dimNonSym == 3) {
          fv->d_grad_lubp_dmesh[2][2][j] += T_i * bfv->d_grad_phi_dmesh[i][2][2][j];
          fv->d_grad_lubp_dmesh[2][0][j] += T_i * bfv->d_grad_phi_dmesh[i][2][0][j];
          fv->d_grad_lubp_dmesh[2][1][j] += T_i * bfv->d_grad_phi_dmesh[i][2][1][j];
          fv->d_grad_lubp_dmesh[0][2][j] += T_i * bfv->d_grad_phi_dmesh[i][0][2][j];
          fv->d_grad_lubp_dmesh[1][2][j] += T_i * bfv->d_grad_phi_dmesh[i][1][2][j];
        }
      }
    }
#endif
  } else if (upd->vp[pg->imtrx][LUBP] != -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_lubp_dmesh[0][0][0]), 0, siz);
  }

  if (pd->gv[LUBP_2]) {
    v = LUBP_2;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];
#ifdef DO_NOT_UNROLL
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_lubp_2_dmesh[0][0][0]), 0, siz);
    for (i = 0; i < vdofs; i++) {
      T_i = *esp->lubp_2[i];
      for (p = 0; p < dimNonSym; p++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_lubp_2_dmesh[p][b][j] += T_i * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
#else
    for (j = 0; j < mdofs; j++) {
      T_i = *esp->lubp_2[0];

      fv->d_grad_lubp_2_dmesh[0][0][j] = T_i * bfv->d_grad_phi_dmesh[0][0][0][j];
      fv->d_grad_lubp_2_dmesh[1][1][j] = T_i * bfv->d_grad_phi_dmesh[0][1][1][j];
      fv->d_grad_lubp_2_dmesh[1][0][j] = T_i * bfv->d_grad_phi_dmesh[0][1][0][j];
      fv->d_grad_lubp_2_dmesh[0][1][j] = T_i * bfv->d_grad_phi_dmesh[0][0][1][j];

      if (dimNonSym == 3) {
        fv->d_grad_lubp_2_dmesh[2][2][j] = T_i * bfv->d_grad_phi_dmesh[0][2][2][j];
        fv->d_grad_lubp_2_dmesh[2][0][j] = T_i * bfv->d_grad_phi_dmesh[0][2][0][j];
        fv->d_grad_lubp_2_dmesh[2][1][j] = T_i * bfv->d_grad_phi_dmesh[0][2][1][j];
        fv->d_grad_lubp_2_dmesh[0][2][j] = T_i * bfv->d_grad_phi_dmesh[0][0][2][j];
        fv->d_grad_lubp_2_dmesh[1][2][j] = T_i * bfv->d_grad_phi_dmesh[0][1][2][j];
      }

      for (i = 1; i < vdofs; i++) {
        T_i = *esp->lubp_2[i];

        fv->d_grad_lubp_2_dmesh[0][0][j] += T_i * bfv->d_grad_phi_dmesh[i][0][0][j];
        fv->d_grad_lubp_2_dmesh[1][1][j] += T_i * bfv->d_grad_phi_dmesh[i][1][1][j];
        fv->d_grad_lubp_2_dmesh[1][0][j] += T_i * bfv->d_grad_phi_dmesh[i][1][0][j];
        fv->d_grad_lubp_2_dmesh[0][1][j] += T_i * bfv->d_grad_phi_dmesh[i][0][1][j];

        if (dimNonSym == 3) {
          fv->d_grad_lubp_2_dmesh[2][2][j] += T_i * bfv->d_grad_phi_dmesh[i][2][2][j];
          fv->d_grad_lubp_2_dmesh[2][0][j] += T_i * bfv->d_grad_phi_dmesh[i][2][0][j];
          fv->d_grad_lubp_2_dmesh[2][1][j] += T_i * bfv->d_grad_phi_dmesh[i][2][1][j];
          fv->d_grad_lubp_2_dmesh[0][2][j] += T_i * bfv->d_grad_phi_dmesh[i][0][2][j];
          fv->d_grad_lubp_2_dmesh[1][2][j] += T_i * bfv->d_grad_phi_dmesh[i][1][2][j];
        }
      }
    }
#endif
  } else if (upd->vp[pg->imtrx][LUBP_2] != -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_lubp_2_dmesh[0][0][0]), 0, siz);
  }

  if (pd->gv[SHELL_PRESS_OPEN]) {
    v = SHELL_PRESS_OPEN;
    bfv = bf[v];
#ifdef DO_NOT_UNROLL
    vdofs = ei[upd->matrix_index[v]]->dof[v];
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_sh_p_open_dmesh[0][0][0]), 0, siz);
    for (i = 0; i < vdofs; i++) {
      T_i = *esp->sh_p_open[i];
      for (p = 0; p < dimNonSym; p++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_sh_p_open_dmesh[p][b][j] += T_i * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
#else
    for (j = 0; j < mdofs; j++) {
      T_i = *esp->sh_p_open[0];

      fv->d_grad_sh_p_open_dmesh[0][0][j] = T_i * bfv->d_grad_phi_dmesh[0][0][0][j];
      fv->d_grad_sh_p_open_dmesh[1][1][j] = T_i * bfv->d_grad_phi_dmesh[0][1][1][j];
      fv->d_grad_sh_p_open_dmesh[1][0][j] = T_i * bfv->d_grad_phi_dmesh[0][1][0][j];
      fv->d_grad_sh_p_open_dmesh[0][1][j] = T_i * bfv->d_grad_phi_dmesh[0][0][1][j];

      if (dimNonSym == 3) {
        fv->d_grad_sh_p_open_dmesh[2][2][j] = T_i * bfv->d_grad_phi_dmesh[0][2][2][j];
        fv->d_grad_sh_p_open_dmesh[2][0][j] = T_i * bfv->d_grad_phi_dmesh[0][2][0][j];
        fv->d_grad_sh_p_open_dmesh[2][1][j] = T_i * bfv->d_grad_phi_dmesh[0][2][1][j];
        fv->d_grad_sh_p_open_dmesh[0][2][j] = T_i * bfv->d_grad_phi_dmesh[0][0][2][j];
        fv->d_grad_sh_p_open_dmesh[1][2][j] = T_i * bfv->d_grad_phi_dmesh[0][1][2][j];
      }

      for (i = 1; i < vdofs; i++) {
        T_i = *esp->sh_p_open[i];

        fv->d_grad_sh_p_open_dmesh[0][0][j] += T_i * bfv->d_grad_phi_dmesh[i][0][0][j];
        fv->d_grad_sh_p_open_dmesh[1][1][j] += T_i * bfv->d_grad_phi_dmesh[i][1][1][j];
        fv->d_grad_sh_p_open_dmesh[1][0][j] += T_i * bfv->d_grad_phi_dmesh[i][1][0][j];
        fv->d_grad_sh_p_open_dmesh[0][1][j] += T_i * bfv->d_grad_phi_dmesh[i][0][1][j];

        if (dimNonSym == 3) {
          fv->d_grad_sh_p_open_dmesh[2][2][j] += T_i * bfv->d_grad_phi_dmesh[i][2][2][j];
          fv->d_grad_sh_p_open_dmesh[2][0][j] += T_i * bfv->d_grad_phi_dmesh[i][2][0][j];
          fv->d_grad_sh_p_open_dmesh[2][1][j] += T_i * bfv->d_grad_phi_dmesh[i][2][1][j];
          fv->d_grad_sh_p_open_dmesh[0][2][j] += T_i * bfv->d_grad_phi_dmesh[i][0][2][j];
          fv->d_grad_sh_p_open_dmesh[1][2][j] += T_i * bfv->d_grad_phi_dmesh[i][1][2][j];
        }
      }
    }
#endif
  } else if (upd->vp[pg->imtrx][SHELL_PRESS_OPEN] != -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_sh_p_open_dmesh[0][0][0]), 0, siz);
  }

  if (pd->gv[SHELL_PRESS_OPEN_2]) {
    v = SHELL_PRESS_OPEN_2;
    bfv = bf[v];
#ifdef DO_NOT_UNROLL
    vdofs = ei[upd->matrix_index[v]]->dof[v];
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_sh_p_open_2_dmesh[0][0][0]), 0, siz);
    for (i = 0; i < vdofs; i++) {
      T_i = *esp->sh_p_open_2[i];
      for (p = 0; p < dimNonSym; p++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_sh_p_open_2_dmesh[p][b][j] += T_i * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
#else
    for (j = 0; j < mdofs; j++) {
      T_i = *esp->sh_p_open_2[0];

      fv->d_grad_sh_p_open_2_dmesh[0][0][j] = T_i * bfv->d_grad_phi_dmesh[0][0][0][j];
      fv->d_grad_sh_p_open_2_dmesh[1][1][j] = T_i * bfv->d_grad_phi_dmesh[0][1][1][j];
      fv->d_grad_sh_p_open_2_dmesh[1][0][j] = T_i * bfv->d_grad_phi_dmesh[0][1][0][j];
      fv->d_grad_sh_p_open_2_dmesh[0][1][j] = T_i * bfv->d_grad_phi_dmesh[0][0][1][j];

      if (dimNonSym == 3) {
        fv->d_grad_sh_p_open_2_dmesh[2][2][j] = T_i * bfv->d_grad_phi_dmesh[0][2][2][j];
        fv->d_grad_sh_p_open_2_dmesh[2][0][j] = T_i * bfv->d_grad_phi_dmesh[0][2][0][j];
        fv->d_grad_sh_p_open_2_dmesh[2][1][j] = T_i * bfv->d_grad_phi_dmesh[0][2][1][j];
        fv->d_grad_sh_p_open_2_dmesh[0][2][j] = T_i * bfv->d_grad_phi_dmesh[0][0][2][j];
        fv->d_grad_sh_p_open_2_dmesh[1][2][j] = T_i * bfv->d_grad_phi_dmesh[0][1][2][j];
      }

      for (i = 1; i < vdofs; i++) {
        T_i = *esp->sh_p_open_2[i];

        fv->d_grad_sh_p_open_2_dmesh[0][0][j] += T_i * bfv->d_grad_phi_dmesh[i][0][0][j];
        fv->d_grad_sh_p_open_2_dmesh[1][1][j] += T_i * bfv->d_grad_phi_dmesh[i][1][1][j];
        fv->d_grad_sh_p_open_2_dmesh[1][0][j] += T_i * bfv->d_grad_phi_dmesh[i][1][0][j];
        fv->d_grad_sh_p_open_2_dmesh[0][1][j] += T_i * bfv->d_grad_phi_dmesh[i][0][1][j];

        if (dimNonSym == 3) {
          fv->d_grad_sh_p_open_2_dmesh[2][2][j] += T_i * bfv->d_grad_phi_dmesh[i][2][2][j];
          fv->d_grad_sh_p_open_2_dmesh[2][0][j] += T_i * bfv->d_grad_phi_dmesh[i][2][0][j];
          fv->d_grad_sh_p_open_2_dmesh[2][1][j] += T_i * bfv->d_grad_phi_dmesh[i][2][1][j];
          fv->d_grad_sh_p_open_2_dmesh[0][2][j] += T_i * bfv->d_grad_phi_dmesh[i][0][2][j];
          fv->d_grad_sh_p_open_2_dmesh[1][2][j] += T_i * bfv->d_grad_phi_dmesh[i][1][2][j];
        }
      }
    }
#endif
  } else if (upd->vp[pg->imtrx][SHELL_PRESS_OPEN_2] != -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_sh_p_open_2_dmesh[0][0][0]), 0, siz);
  }

  if (pd->gv[SHELL_FILMP]) {
    v = SHELL_FILMP;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];
#ifdef DO_NOT_UNROLL
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_sh_fp_dmesh[0][0][0]), 0, siz);
    for (i = 0; i < vdofs; i++) {
      T_i = *esp->sh_fp[i];
      for (p = 0; p < dimNonSym; p++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_sh_fp_dmesh[p][b][j] += T_i * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
#else
    for (j = 0; j < mdofs; j++) {
      T_i = *esp->sh_fp[0];

      fv->d_grad_sh_fp_dmesh[0][0][j] = T_i * bfv->d_grad_phi_dmesh[0][0][0][j];
      fv->d_grad_sh_fp_dmesh[1][1][j] = T_i * bfv->d_grad_phi_dmesh[0][1][1][j];
      fv->d_grad_sh_fp_dmesh[1][0][j] = T_i * bfv->d_grad_phi_dmesh[0][1][0][j];
      fv->d_grad_sh_fp_dmesh[0][1][j] = T_i * bfv->d_grad_phi_dmesh[0][0][1][j];

      if (dimNonSym == 3) {
        fv->d_grad_sh_fp_dmesh[2][2][j] = T_i * bfv->d_grad_phi_dmesh[0][2][2][j];
        fv->d_grad_sh_fp_dmesh[2][0][j] = T_i * bfv->d_grad_phi_dmesh[0][2][0][j];
        fv->d_grad_sh_fp_dmesh[2][1][j] = T_i * bfv->d_grad_phi_dmesh[0][2][1][j];
        fv->d_grad_sh_fp_dmesh[0][2][j] = T_i * bfv->d_grad_phi_dmesh[0][0][2][j];
        fv->d_grad_sh_fp_dmesh[1][2][j] = T_i * bfv->d_grad_phi_dmesh[0][1][2][j];
      }

      for (i = 1; i < vdofs; i++) {
        T_i = *esp->sh_fp[i];

        fv->d_grad_sh_fp_dmesh[0][0][j] += T_i * bfv->d_grad_phi_dmesh[i][0][0][j];
        fv->d_grad_sh_fp_dmesh[1][1][j] += T_i * bfv->d_grad_phi_dmesh[i][1][1][j];
        fv->d_grad_sh_fp_dmesh[1][0][j] += T_i * bfv->d_grad_phi_dmesh[i][1][0][j];
        fv->d_grad_sh_fp_dmesh[0][1][j] += T_i * bfv->d_grad_phi_dmesh[i][0][1][j];

        if (dimNonSym == 3) {
          fv->d_grad_sh_fp_dmesh[2][2][j] += T_i * bfv->d_grad_phi_dmesh[i][2][2][j];
          fv->d_grad_sh_fp_dmesh[2][0][j] += T_i * bfv->d_grad_phi_dmesh[i][2][0][j];
          fv->d_grad_sh_fp_dmesh[2][1][j] += T_i * bfv->d_grad_phi_dmesh[i][2][1][j];
          fv->d_grad_sh_fp_dmesh[0][2][j] += T_i * bfv->d_grad_phi_dmesh[i][0][2][j];
          fv->d_grad_sh_fp_dmesh[1][2][j] += T_i * bfv->d_grad_phi_dmesh[i][1][2][j];
        }
      }
    }
#endif
  } else if (upd->vp[pg->imtrx][SHELL_FILMP] != -1) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_sh_fp_dmesh[0][0][0]), 0, siz);
  }

  if (pd->gv[SHELL_FILMH]) {
    v = SHELL_FILMH;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];
#ifdef DO_NOT_UNROLL
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_sh_fh_dmesh[0][0][0]), 0, siz);
    for (i = 0; i < vdofs; i++) {
      T_i = *esp->sh_fh[i];
      for (p = 0; p < dimNonSym; p++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_sh_fh_dmesh[p][b][j] += T_i * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
#else
    for (j = 0; j < mdofs; j++) {
      T_i = *esp->sh_fh[0];

      fv->d_grad_sh_fh_dmesh[0][0][j] = T_i * bfv->d_grad_phi_dmesh[0][0][0][j];
      fv->d_grad_sh_fh_dmesh[1][1][j] = T_i * bfv->d_grad_phi_dmesh[0][1][1][j];
      fv->d_grad_sh_fh_dmesh[1][0][j] = T_i * bfv->d_grad_phi_dmesh[0][1][0][j];
      fv->d_grad_sh_fh_dmesh[0][1][j] = T_i * bfv->d_grad_phi_dmesh[0][0][1][j];

      if (dimNonSym == 3) {
        fv->d_grad_sh_fh_dmesh[2][2][j] = T_i * bfv->d_grad_phi_dmesh[0][2][2][j];
        fv->d_grad_sh_fh_dmesh[2][0][j] = T_i * bfv->d_grad_phi_dmesh[0][2][0][j];
        fv->d_grad_sh_fh_dmesh[2][1][j] = T_i * bfv->d_grad_phi_dmesh[0][2][1][j];
        fv->d_grad_sh_fh_dmesh[0][2][j] = T_i * bfv->d_grad_phi_dmesh[0][0][2][j];
        fv->d_grad_sh_fh_dmesh[1][2][j] = T_i * bfv->d_grad_phi_dmesh[0][1][2][j];
      }

      for (i = 1; i < vdofs; i++) {
        T_i = *esp->sh_fh[i];

        fv->d_grad_sh_fh_dmesh[0][0][j] += T_i * bfv->d_grad_phi_dmesh[i][0][0][j];
        fv->d_grad_sh_fh_dmesh[1][1][j] += T_i * bfv->d_grad_phi_dmesh[i][1][1][j];
        fv->d_grad_sh_fh_dmesh[1][0][j] += T_i * bfv->d_grad_phi_dmesh[i][1][0][j];
        fv->d_grad_sh_fh_dmesh[0][1][j] += T_i * bfv->d_grad_phi_dmesh[i][0][1][j];

        if (dimNonSym == 3) {
          fv->d_grad_sh_fh_dmesh[2][2][j] += T_i * bfv->d_grad_phi_dmesh[i][2][2][j];
          fv->d_grad_sh_fh_dmesh[2][0][j] += T_i * bfv->d_grad_phi_dmesh[i][2][0][j];
          fv->d_grad_sh_fh_dmesh[2][1][j] += T_i * bfv->d_grad_phi_dmesh[i][2][1][j];
          fv->d_grad_sh_fh_dmesh[0][2][j] += T_i * bfv->d_grad_phi_dmesh[i][0][2][j];
          fv->d_grad_sh_fh_dmesh[1][2][j] += T_i * bfv->d_grad_phi_dmesh[i][1][2][j];
        }
      }
    }
#endif
  } else if (upd->vp[pg->imtrx][SHELL_FILMH] != -1) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_sh_fh_dmesh[0][0][0]), 0, siz);
  }

  if (pd->gv[SHELL_PARTC]) {
    v = SHELL_PARTC;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];
#ifdef DO_NOT_UNROLL
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_sh_pc_dmesh[0][0][0]), 0, siz);
    for (i = 0; i < vdofs; i++) {
      T_i = *esp->sh_pc[i];
      for (p = 0; p < dimNonSym; p++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_sh_pc_dmesh[p][b][j] += T_i * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
#else
    for (j = 0; j < mdofs; j++) {
      T_i = *esp->sh_pc[0];

      fv->d_grad_sh_pc_dmesh[0][0][j] = T_i * bfv->d_grad_phi_dmesh[0][0][0][j];
      fv->d_grad_sh_pc_dmesh[1][1][j] = T_i * bfv->d_grad_phi_dmesh[0][1][1][j];
      fv->d_grad_sh_pc_dmesh[1][0][j] = T_i * bfv->d_grad_phi_dmesh[0][1][0][j];
      fv->d_grad_sh_pc_dmesh[0][1][j] = T_i * bfv->d_grad_phi_dmesh[0][0][1][j];

      if (dimNonSym == 3) {
        fv->d_grad_sh_pc_dmesh[2][2][j] = T_i * bfv->d_grad_phi_dmesh[0][2][2][j];
        fv->d_grad_sh_pc_dmesh[2][0][j] = T_i * bfv->d_grad_phi_dmesh[0][2][0][j];
        fv->d_grad_sh_pc_dmesh[2][1][j] = T_i * bfv->d_grad_phi_dmesh[0][2][1][j];
        fv->d_grad_sh_pc_dmesh[0][2][j] = T_i * bfv->d_grad_phi_dmesh[0][0][2][j];
        fv->d_grad_sh_pc_dmesh[1][2][j] = T_i * bfv->d_grad_phi_dmesh[0][1][2][j];
      }

      for (i = 1; i < vdofs; i++) {
        T_i = *esp->sh_pc[i];

        fv->d_grad_sh_pc_dmesh[0][0][j] += T_i * bfv->d_grad_phi_dmesh[i][0][0][j];
        fv->d_grad_sh_pc_dmesh[1][1][j] += T_i * bfv->d_grad_phi_dmesh[i][1][1][j];
        fv->d_grad_sh_pc_dmesh[1][0][j] += T_i * bfv->d_grad_phi_dmesh[i][1][0][j];
        fv->d_grad_sh_pc_dmesh[0][1][j] += T_i * bfv->d_grad_phi_dmesh[i][0][1][j];

        if (dimNonSym == 3) {
          fv->d_grad_sh_pc_dmesh[2][2][j] += T_i * bfv->d_grad_phi_dmesh[i][2][2][j];
          fv->d_grad_sh_pc_dmesh[2][0][j] += T_i * bfv->d_grad_phi_dmesh[i][2][0][j];
          fv->d_grad_sh_pc_dmesh[2][1][j] += T_i * bfv->d_grad_phi_dmesh[i][2][1][j];
          fv->d_grad_sh_pc_dmesh[0][2][j] += T_i * bfv->d_grad_phi_dmesh[i][0][2][j];
          fv->d_grad_sh_pc_dmesh[1][2][j] += T_i * bfv->d_grad_phi_dmesh[i][1][2][j];
        }
      }
    }
#endif
  } else if (upd->vp[pg->imtrx][SHELL_PARTC] != -1) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_sh_pc_dmesh[0][0][0]), 0, siz);
  }

  if (pd->gv[SHELL_TEMPERATURE]) {
    v = SHELL_TEMPERATURE;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];
#ifdef DO_NOT_UNROLL
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_sh_t_dmesh[0][0][0]), 0, siz);
    for (i = 0; i < vdofs; i++) {
      T_i = *esp->sh_t[i];
      for (p = 0; p < dimNonSym; p++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_sh_t_dmesh[p][b][j] += T_i * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
#else
    for (j = 0; j < mdofs; j++) {
      T_i = *esp->sh_t[0];

      fv->d_grad_sh_t_dmesh[0][0][j] = T_i * bfv->d_grad_phi_dmesh[0][0][0][j];
      fv->d_grad_sh_t_dmesh[1][1][j] = T_i * bfv->d_grad_phi_dmesh[0][1][1][j];
      fv->d_grad_sh_t_dmesh[1][0][j] = T_i * bfv->d_grad_phi_dmesh[0][1][0][j];
      fv->d_grad_sh_t_dmesh[0][1][j] = T_i * bfv->d_grad_phi_dmesh[0][0][1][j];

      if (dimNonSym == 3) {
        fv->d_grad_sh_t_dmesh[2][2][j] = T_i * bfv->d_grad_phi_dmesh[0][2][2][j];
        fv->d_grad_sh_t_dmesh[2][0][j] = T_i * bfv->d_grad_phi_dmesh[0][2][0][j];
        fv->d_grad_sh_t_dmesh[2][1][j] = T_i * bfv->d_grad_phi_dmesh[0][2][1][j];
        fv->d_grad_sh_t_dmesh[0][2][j] = T_i * bfv->d_grad_phi_dmesh[0][0][2][j];
        fv->d_grad_sh_t_dmesh[1][2][j] = T_i * bfv->d_grad_phi_dmesh[0][1][2][j];
      }

      for (i = 1; i < vdofs; i++) {
        T_i = *esp->sh_t[i];

        fv->d_grad_sh_t_dmesh[0][0][j] += T_i * bfv->d_grad_phi_dmesh[i][0][0][j];
        fv->d_grad_sh_t_dmesh[1][1][j] += T_i * bfv->d_grad_phi_dmesh[i][1][1][j];
        fv->d_grad_sh_t_dmesh[1][0][j] += T_i * bfv->d_grad_phi_dmesh[i][1][0][j];
        fv->d_grad_sh_t_dmesh[0][1][j] += T_i * bfv->d_grad_phi_dmesh[i][0][1][j];

        if (dimNonSym == 3) {
          fv->d_grad_sh_t_dmesh[2][2][j] += T_i * bfv->d_grad_phi_dmesh[i][2][2][j];
          fv->d_grad_sh_t_dmesh[2][0][j] += T_i * bfv->d_grad_phi_dmesh[i][2][0][j];
          fv->d_grad_sh_t_dmesh[2][1][j] += T_i * bfv->d_grad_phi_dmesh[i][2][1][j];
          fv->d_grad_sh_t_dmesh[0][2][j] += T_i * bfv->d_grad_phi_dmesh[i][0][2][j];
          fv->d_grad_sh_t_dmesh[1][2][j] += T_i * bfv->d_grad_phi_dmesh[i][1][2][j];
        }
      }
    }
#endif
  } else if (upd->vp[pg->imtrx][SHELL_TEMPERATURE] != -1) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_sh_t_dmesh[0][0][0]), 0, siz);
  }

  if (pd->gv[EDDY_NU]) {
    v = EDDY_NU;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];
#ifdef DO_NOT_UNROLL
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_eddy_nu_dmesh[0][0][0]), 0, siz);
    for (i = 0; i < vdofs; i++) {
      T_i = *esp->eddy_nu[i];
      for (p = 0; p < dimNonSym; p++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_eddy_nu_dmesh[p][b][j] += T_i * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
#else
    for (j = 0; j < mdofs; j++) {
      T_i = *esp->eddy_nu[0];

      fv->d_grad_eddy_nu_dmesh[0][0][j] = T_i * bfv->d_grad_phi_dmesh[0][0][0][j];
      fv->d_grad_eddy_nu_dmesh[1][1][j] = T_i * bfv->d_grad_phi_dmesh[0][1][1][j];
      fv->d_grad_eddy_nu_dmesh[1][0][j] = T_i * bfv->d_grad_phi_dmesh[0][1][0][j];
      fv->d_grad_eddy_nu_dmesh[0][1][j] = T_i * bfv->d_grad_phi_dmesh[0][0][1][j];

      if (dimNonSym == 3) {
        fv->d_grad_eddy_nu_dmesh[2][2][j] = T_i * bfv->d_grad_phi_dmesh[0][2][2][j];
        fv->d_grad_eddy_nu_dmesh[2][0][j] = T_i * bfv->d_grad_phi_dmesh[0][2][0][j];
        fv->d_grad_eddy_nu_dmesh[2][1][j] = T_i * bfv->d_grad_phi_dmesh[0][2][1][j];
        fv->d_grad_eddy_nu_dmesh[0][2][j] = T_i * bfv->d_grad_phi_dmesh[0][0][2][j];
        fv->d_grad_eddy_nu_dmesh[1][2][j] = T_i * bfv->d_grad_phi_dmesh[0][1][2][j];
      }

      for (i = 1; i < vdofs; i++) {
        T_i = *esp->eddy_nu[i];

        fv->d_grad_eddy_nu_dmesh[0][0][j] += T_i * bfv->d_grad_phi_dmesh[i][0][0][j];
        fv->d_grad_eddy_nu_dmesh[1][1][j] += T_i * bfv->d_grad_phi_dmesh[i][1][1][j];
        fv->d_grad_eddy_nu_dmesh[1][0][j] += T_i * bfv->d_grad_phi_dmesh[i][1][0][j];
        fv->d_grad_eddy_nu_dmesh[0][1][j] += T_i * bfv->d_grad_phi_dmesh[i][0][1][j];

        if (dimNonSym == 3) {
          fv->d_grad_eddy_nu_dmesh[2][2][j] += T_i * bfv->d_grad_phi_dmesh[i][2][2][j];
          fv->d_grad_eddy_nu_dmesh[2][0][j] += T_i * bfv->d_grad_phi_dmesh[i][2][0][j];
          fv->d_grad_eddy_nu_dmesh[2][1][j] += T_i * bfv->d_grad_phi_dmesh[i][2][1][j];
          fv->d_grad_eddy_nu_dmesh[0][2][j] += T_i * bfv->d_grad_phi_dmesh[i][0][2][j];
          fv->d_grad_eddy_nu_dmesh[1][2][j] += T_i * bfv->d_grad_phi_dmesh[i][1][2][j];
        }
      }
    }
#endif
  } else if (upd->vp[pg->imtrx][EDDY_NU] != -1) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_eddy_nu_dmesh[0][0][0]), 0, siz);
  }

  if (pd->gv[LIGHT_INTP]) {
    v = LIGHT_INTP;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_poynt_dmesh[0][0][0][0]), 0, siz);
    for (i = 0; i < vdofs; i++) {
      T_i = *esp->poynt[0][i];
      for (p = 0; p < dimNonSym; p++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_poynt_dmesh[0][p][b][j] += T_i * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
  } else if (upd->vp[pg->imtrx][LIGHT_INTP] != -1) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_poynt_dmesh[0][0][0][0]), 0, siz);
  }

  if (pd->gv[LIGHT_INTM]) {
    v = LIGHT_INTM;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_poynt_dmesh[1][0][0][0]), 0, siz);
    for (i = 0; i < vdofs; i++) {
      T_i = *esp->poynt[1][i];
      for (p = 0; p < dimNonSym; p++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_poynt_dmesh[1][p][b][j] += T_i * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
  } else if (upd->vp[pg->imtrx][LIGHT_INTM] != -1) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_poynt_dmesh[1][0][0][0]), 0, siz);
  }

  if (pd->gv[LIGHT_INTD]) {
    v = LIGHT_INTD;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_poynt_dmesh[2][0][0][0]), 0, siz);
    for (i = 0; i < vdofs; i++) {
      T_i = *esp->poynt[2][i];
      for (p = 0; p < dimNonSym; p++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_poynt_dmesh[2][p][b][j] += T_i * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
  } else if (upd->vp[pg->imtrx][LIGHT_INTD] != -1) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_poynt_dmesh[2][0][0][0]), 0, siz);
  }

  if (pd->gv[TFMP_PRES]) {
    v = TFMP_PRES;
    bfv = bf[v];
    vdofs = ei[pg->imtrx]->dof[v];
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_tfmp_pres_dmesh[0][0][0]), 0, siz);
    for (i = 0; i < vdofs; i++) {
      T_i = *esp->tfmp_pres[i];
      for (p = 0; p < dimNonSym; p++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_tfmp_pres_dmesh[p][b][j] += T_i * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
  } else if (upd->vp[pg->imtrx][TFMP_PRES] != -1) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_tfmp_pres_dmesh[0][0][0]), 0, siz);
  }

  if (pd->gv[TFMP_SAT]) {
    v = TFMP_SAT;
    bfv = bf[v];
    vdofs = ei[pg->imtrx]->dof[v];
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_tfmp_sat_dmesh[0][0][0]), 0, siz);
    for (i = 0; i < vdofs; i++) {
      T_i = *esp->tfmp_sat[i];
      for (p = 0; p < dimNonSym; p++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_tfmp_sat_dmesh[p][b][j] += T_i * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
  } else if (upd->vp[pg->imtrx][TFMP_PRES] != -1) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_tfmp_pres_dmesh[0][0][0]), 0, siz);
  }

  if (pd->gv[RESTIME]) {
    v = RESTIME;
    bfv = bf[v];
    vdofs = ei[pg->imtrx]->dof[v];
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_restime_dmesh[0][0][0]), 0, siz);
    for (i = 0; i < vdofs; i++) {
      T_i = *esp->restime[i];
      for (p = 0; p < dimNonSym; p++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_restime_dmesh[p][b][j] += T_i * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
  } else if (upd->vp[pg->imtrx][RESTIME] != -1) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_restime_dmesh[0][0][0]), 0, siz);
  }

  if (pd->gv[SHELL_LUB_CURV]) {
    v = SHELL_LUB_CURV;
    bfv = bf[v];
#ifdef DO_NOT_UNROLL
    vdofs = ei[upd->matrix_index[v]]->dof[v];
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_sh_l_curv_dmesh[0][0][0]), 0, siz);
    for (i = 0; i < vdofs; i++) {
      T_i = *esp->sh_l_curv[i];
      for (p = 0; p < dimNonSym; p++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_sh_l_curv_dmesh[p][b][j] += T_i * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
#else
    for (j = 0; j < mdofs; j++) {
      T_i = *esp->sh_l_curv[0];

      fv->d_grad_sh_l_curv_dmesh[0][0][j] = T_i * bfv->d_grad_phi_dmesh[0][0][0][j];
      fv->d_grad_sh_l_curv_dmesh[1][1][j] = T_i * bfv->d_grad_phi_dmesh[0][1][1][j];
      fv->d_grad_sh_l_curv_dmesh[1][0][j] = T_i * bfv->d_grad_phi_dmesh[0][1][0][j];
      fv->d_grad_sh_l_curv_dmesh[0][1][j] = T_i * bfv->d_grad_phi_dmesh[0][0][1][j];

      if (dimNonSym == 3) {
        fv->d_grad_sh_l_curv_dmesh[2][2][j] = T_i * bfv->d_grad_phi_dmesh[0][2][2][j];
        fv->d_grad_sh_l_curv_dmesh[2][0][j] = T_i * bfv->d_grad_phi_dmesh[0][2][0][j];
        fv->d_grad_sh_l_curv_dmesh[2][1][j] = T_i * bfv->d_grad_phi_dmesh[0][2][1][j];
        fv->d_grad_sh_l_curv_dmesh[0][2][j] = T_i * bfv->d_grad_phi_dmesh[0][0][2][j];
        fv->d_grad_sh_l_curv_dmesh[1][2][j] = T_i * bfv->d_grad_phi_dmesh[0][1][2][j];
      }

      for (i = 1; i < vdofs; i++) {
        T_i = *esp->sh_l_curv[i];

        fv->d_grad_sh_l_curv_dmesh[0][0][j] += T_i * bfv->d_grad_phi_dmesh[i][0][0][j];
        fv->d_grad_sh_l_curv_dmesh[1][1][j] += T_i * bfv->d_grad_phi_dmesh[i][1][1][j];
        fv->d_grad_sh_l_curv_dmesh[1][0][j] += T_i * bfv->d_grad_phi_dmesh[i][1][0][j];
        fv->d_grad_sh_l_curv_dmesh[0][1][j] += T_i * bfv->d_grad_phi_dmesh[i][0][1][j];

        if (dimNonSym == 3) {
          fv->d_grad_sh_l_curv_dmesh[2][2][j] += T_i * bfv->d_grad_phi_dmesh[i][2][2][j];
          fv->d_grad_sh_l_curv_dmesh[2][0][j] += T_i * bfv->d_grad_phi_dmesh[i][2][0][j];
          fv->d_grad_sh_l_curv_dmesh[2][1][j] += T_i * bfv->d_grad_phi_dmesh[i][2][1][j];
          fv->d_grad_sh_l_curv_dmesh[0][2][j] += T_i * bfv->d_grad_phi_dmesh[i][0][2][j];
          fv->d_grad_sh_l_curv_dmesh[1][2][j] += T_i * bfv->d_grad_phi_dmesh[i][1][2][j];
        }
      }
    }
#endif
  } else if (upd->vp[pg->imtrx][SHELL_LUB_CURV] != -1) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_sh_l_curv_dmesh[0][0][0]), 0, siz);
  }

  if (pd->gv[SHELL_LUB_CURV_2]) {
    v = SHELL_LUB_CURV_2;
    bfv = bf[v];
#ifdef DO_NOT_UNROLL
    vdofs = ei[upd->matrix_index[v]]->dof[v];
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_sh_l_curv_2_dmesh[0][0][0]), 0, siz);
    for (i = 0; i < vdofs; i++) {
      T_i = *esp->sh_l_curv_2[i];
      for (p = 0; p < dimNonSym; p++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_sh_l_curv_2_dmesh[p][b][j] += T_i * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
#else
    for (j = 0; j < mdofs; j++) {
      T_i = *esp->sh_l_curv_2[0];

      fv->d_grad_sh_l_curv_2_dmesh[0][0][j] = T_i * bfv->d_grad_phi_dmesh[0][0][0][j];
      fv->d_grad_sh_l_curv_2_dmesh[1][1][j] = T_i * bfv->d_grad_phi_dmesh[0][1][1][j];
      fv->d_grad_sh_l_curv_2_dmesh[1][0][j] = T_i * bfv->d_grad_phi_dmesh[0][1][0][j];
      fv->d_grad_sh_l_curv_2_dmesh[0][1][j] = T_i * bfv->d_grad_phi_dmesh[0][0][1][j];

      if (dimNonSym == 3) {
        fv->d_grad_sh_l_curv_2_dmesh[2][2][j] = T_i * bfv->d_grad_phi_dmesh[0][2][2][j];
        fv->d_grad_sh_l_curv_2_dmesh[2][0][j] = T_i * bfv->d_grad_phi_dmesh[0][2][0][j];
        fv->d_grad_sh_l_curv_2_dmesh[2][1][j] = T_i * bfv->d_grad_phi_dmesh[0][2][1][j];
        fv->d_grad_sh_l_curv_2_dmesh[0][2][j] = T_i * bfv->d_grad_phi_dmesh[0][0][2][j];
        fv->d_grad_sh_l_curv_2_dmesh[1][2][j] = T_i * bfv->d_grad_phi_dmesh[0][1][2][j];
      }

      for (i = 1; i < vdofs; i++) {
        T_i = *esp->sh_l_curv[i];

        fv->d_grad_sh_l_curv_2_dmesh[0][0][j] += T_i * bfv->d_grad_phi_dmesh[i][0][0][j];
        fv->d_grad_sh_l_curv_2_dmesh[1][1][j] += T_i * bfv->d_grad_phi_dmesh[i][1][1][j];
        fv->d_grad_sh_l_curv_2_dmesh[1][0][j] += T_i * bfv->d_grad_phi_dmesh[i][1][0][j];
        fv->d_grad_sh_l_curv_2_dmesh[0][1][j] += T_i * bfv->d_grad_phi_dmesh[i][0][1][j];

        if (dimNonSym == 3) {
          fv->d_grad_sh_l_curv_2_dmesh[2][2][j] += T_i * bfv->d_grad_phi_dmesh[i][2][2][j];
          fv->d_grad_sh_l_curv_2_dmesh[2][0][j] += T_i * bfv->d_grad_phi_dmesh[i][2][0][j];
          fv->d_grad_sh_l_curv_2_dmesh[2][1][j] += T_i * bfv->d_grad_phi_dmesh[i][2][1][j];
          fv->d_grad_sh_l_curv_2_dmesh[0][2][j] += T_i * bfv->d_grad_phi_dmesh[i][0][2][j];
          fv->d_grad_sh_l_curv_2_dmesh[1][2][j] += T_i * bfv->d_grad_phi_dmesh[i][1][2][j];
        }
      }
    }
#endif
  } else if (upd->vp[pg->imtrx][SHELL_LUB_CURV_2] != -1) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_sh_l_curv_2_dmesh[0][0][0]), 0, siz);
  }

  /*
   * d(grad(qs))/dmesh
   */

  if (pd->gv[SURF_CHARGE]) {
    v = SURF_CHARGE;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];
#ifdef DO_NOT_UNROLL
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_qs_dmesh[0][0][0]), 0, siz);
    for (p = 0; p < dimNonSym; p++) {
      for (b = 0; b < dim; b++) {
        for (j = 0; j < mdofs; j++) {
          for (i = 0; i < vdofs; i++) {
            fv->d_grad_qs_dmesh[p][b][j] += *esp->qs[i] * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
#else
    for (j = 0; j < mdofs; j++) {
      fv->d_grad_qs_dmesh[0][0][j] = *esp->qs[0] * bfv->d_grad_phi_dmesh[0][0][0][j];
      fv->d_grad_qs_dmesh[1][1][j] = *esp->qs[0] * bfv->d_grad_phi_dmesh[0][1][1][j];
      fv->d_grad_qs_dmesh[1][0][j] = *esp->qs[0] * bfv->d_grad_phi_dmesh[0][1][0][j];
      fv->d_grad_qs_dmesh[0][1][j] = *esp->qs[0] * bfv->d_grad_phi_dmesh[0][0][1][j];

      if (dimNonSym == 3) {
        fv->d_grad_qs_dmesh[2][2][j] = *esp->qs[0] * bfv->d_grad_phi_dmesh[0][2][2][j];
        fv->d_grad_qs_dmesh[2][0][j] = *esp->qs[0] * bfv->d_grad_phi_dmesh[0][2][0][j];
        fv->d_grad_qs_dmesh[2][1][j] = *esp->qs[0] * bfv->d_grad_phi_dmesh[0][2][1][j];
        fv->d_grad_qs_dmesh[0][2][j] = *esp->qs[0] * bfv->d_grad_phi_dmesh[0][0][2][j];
        fv->d_grad_qs_dmesh[1][2][j] = *esp->qs[0] * bfv->d_grad_phi_dmesh[0][1][2][j];
      }

      for (i = 1; i < vdofs; i++) {
        fv->d_grad_qs_dmesh[0][0][j] += *esp->qs[i] * bfv->d_grad_phi_dmesh[i][0][0][j];
        fv->d_grad_qs_dmesh[1][1][j] += *esp->qs[i] * bfv->d_grad_phi_dmesh[i][1][1][j];
        fv->d_grad_qs_dmesh[1][0][j] += *esp->qs[i] * bfv->d_grad_phi_dmesh[i][1][0][j];
        fv->d_grad_qs_dmesh[0][1][j] += *esp->qs[i] * bfv->d_grad_phi_dmesh[i][0][1][j];

        if (dimNonSym == 3) {
          fv->d_grad_qs_dmesh[2][2][j] += *esp->qs[i] * bfv->d_grad_phi_dmesh[i][2][2][j];
          fv->d_grad_qs_dmesh[2][0][j] += *esp->qs[i] * bfv->d_grad_phi_dmesh[i][2][0][j];
          fv->d_grad_qs_dmesh[2][1][j] += *esp->qs[i] * bfv->d_grad_phi_dmesh[i][2][1][j];
          fv->d_grad_qs_dmesh[0][2][j] += *esp->qs[i] * bfv->d_grad_phi_dmesh[i][0][2][j];
          fv->d_grad_qs_dmesh[1][2][j] += *esp->qs[i] * bfv->d_grad_phi_dmesh[i][1][2][j];
        }
      }
    }
#endif
  } else if (upd->vp[pg->imtrx][SURF_CHARGE] != -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_qs_dmesh[0][0][0]), 0, siz);
  }

  /*
   * d(grad(sh_J))/dmesh
   */

  /* This is needed for bulk assembly of shell diffusion KBC */
  if (pd->gv[SHELL_DIFF_FLUX] && ei[upd->matrix_index[SHELL_DIFF_FLUX]]->dof[SHELL_DIFF_FLUX] > 0) {
    v = SHELL_DIFF_FLUX;
    bfv = bf[v];
#ifdef DO_NOT_UNROLL
    vdofs = ei[upd->matrix_index[v]]->dof[v];
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_sh_J_dmesh[0][0][0]), 0, siz);
    for (p = 0; p < dimNonSym; p++) {
      for (b = 0; b < dim; b++) {
        for (j = 0; j < mdofs; j++) {
          for (i = 0; i < vdofs; i++) {
            fv->d_grad_sh_J_dmesh[p][b][j] += *esp->sh_J[i] * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
#else
    for (j = 0; j < mdofs; j++) {
      fv->d_grad_sh_J_dmesh[0][0][j] = *esp->sh_J[0] * bfv->d_grad_phi_dmesh[0][0][0][j];
      fv->d_grad_sh_J_dmesh[1][1][j] = *esp->sh_J[0] * bfv->d_grad_phi_dmesh[0][1][1][j];
      fv->d_grad_sh_J_dmesh[1][0][j] = *esp->sh_J[0] * bfv->d_grad_phi_dmesh[0][1][0][j];
      fv->d_grad_sh_J_dmesh[0][1][j] = *esp->sh_J[0] * bfv->d_grad_phi_dmesh[0][0][1][j];

      if (dimNonSym == 3) {
        fv->d_grad_sh_J_dmesh[2][2][j] = *esp->sh_J[0] * bfv->d_grad_phi_dmesh[0][2][2][j];
        fv->d_grad_sh_J_dmesh[2][0][j] = *esp->sh_J[0] * bfv->d_grad_phi_dmesh[0][2][0][j];
        fv->d_grad_sh_J_dmesh[2][1][j] = *esp->sh_J[0] * bfv->d_grad_phi_dmesh[0][2][1][j];
        fv->d_grad_sh_J_dmesh[0][2][j] = *esp->sh_J[0] * bfv->d_grad_phi_dmesh[0][0][2][j];
        fv->d_grad_sh_J_dmesh[1][2][j] = *esp->sh_J[0] * bfv->d_grad_phi_dmesh[0][1][2][j];
      }

      for (i = 1; i < vdofs; i++) {
        fv->d_grad_sh_J_dmesh[0][0][j] += *esp->sh_J[i] * bfv->d_grad_phi_dmesh[i][0][0][j];
        fv->d_grad_sh_J_dmesh[1][1][j] += *esp->sh_J[i] * bfv->d_grad_phi_dmesh[i][1][1][j];
        fv->d_grad_sh_J_dmesh[1][0][j] += *esp->sh_J[i] * bfv->d_grad_phi_dmesh[i][1][0][j];
        fv->d_grad_sh_J_dmesh[0][1][j] += *esp->sh_J[i] * bfv->d_grad_phi_dmesh[i][0][1][j];

        if (dimNonSym == 3) {
          fv->d_grad_sh_J_dmesh[2][2][j] += *esp->sh_J[i] * bfv->d_grad_phi_dmesh[i][2][2][j];
          fv->d_grad_sh_J_dmesh[2][0][j] += *esp->sh_J[i] * bfv->d_grad_phi_dmesh[i][2][0][j];
          fv->d_grad_sh_J_dmesh[2][1][j] += *esp->sh_J[i] * bfv->d_grad_phi_dmesh[i][2][1][j];
          fv->d_grad_sh_J_dmesh[0][2][j] += *esp->sh_J[i] * bfv->d_grad_phi_dmesh[i][0][2][j];
          fv->d_grad_sh_J_dmesh[1][2][j] += *esp->sh_J[i] * bfv->d_grad_phi_dmesh[i][1][2][j];
        }
      }
    }
#endif
  } else if (upd->vp[pg->imtrx][SHELL_DIFF_FLUX] != -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_sh_J_dmesh[0][0][0]), 0, siz);
  }

  /*
   * d(grad(SH))/dmesh
   */
  if (pd->gv[SHEAR_RATE]) {
    v = SHEAR_RATE;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];
#ifdef DO_NOT_UNROLL
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(fv->d_grad_SH_dmesh, 0, siz);

    for (p = 0; p < dimNonSym; p++) {
      for (b = 0; b < dim; b++) {
        for (j = 0; j < mdofs; j++) {
          for (i = 0; i < vdofs; i++) {
            fv->d_grad_SH_dmesh[p][b][j] += *esp->SH[i] * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
#else
    for (j = 0; j < mdofs; j++) {
      fv->d_grad_SH_dmesh[0][0][j] = *esp->SH[0] * bfv->d_grad_phi_dmesh[0][0][0][j];
      fv->d_grad_SH_dmesh[1][1][j] = *esp->SH[0] * bfv->d_grad_phi_dmesh[0][1][1][j];
      fv->d_grad_SH_dmesh[1][0][j] = *esp->SH[0] * bfv->d_grad_phi_dmesh[0][1][0][j];
      fv->d_grad_SH_dmesh[0][1][j] = *esp->SH[0] * bfv->d_grad_phi_dmesh[0][0][1][j];

      if (dimNonSym == 3) {
        fv->d_grad_SH_dmesh[2][2][j] = *esp->SH[0] * bfv->d_grad_phi_dmesh[0][2][2][j];
        fv->d_grad_SH_dmesh[2][0][j] = *esp->SH[0] * bfv->d_grad_phi_dmesh[0][2][0][j];
        fv->d_grad_SH_dmesh[2][1][j] = *esp->SH[0] * bfv->d_grad_phi_dmesh[0][2][1][j];
        fv->d_grad_SH_dmesh[0][2][j] = *esp->SH[0] * bfv->d_grad_phi_dmesh[0][0][2][j];
        fv->d_grad_SH_dmesh[1][2][j] = *esp->SH[0] * bfv->d_grad_phi_dmesh[0][1][2][j];
      }

      for (i = 1; i < vdofs; i++) {
        fv->d_grad_SH_dmesh[0][0][j] += *esp->SH[i] * bfv->d_grad_phi_dmesh[i][0][0][j];
        fv->d_grad_SH_dmesh[1][1][j] += *esp->SH[i] * bfv->d_grad_phi_dmesh[i][1][1][j];
        fv->d_grad_SH_dmesh[1][0][j] += *esp->SH[i] * bfv->d_grad_phi_dmesh[i][1][0][j];
        fv->d_grad_SH_dmesh[0][1][j] += *esp->SH[i] * bfv->d_grad_phi_dmesh[i][0][1][j];

        if (dimNonSym == 3) {
          fv->d_grad_api_dmesh[2][2][j] += *esp->SH[i] * bfv->d_grad_phi_dmesh[i][2][2][j];
          fv->d_grad_api_dmesh[2][0][j] += *esp->SH[i] * bfv->d_grad_phi_dmesh[i][2][0][j];
          fv->d_grad_api_dmesh[2][1][j] += *esp->SH[i] * bfv->d_grad_phi_dmesh[i][2][1][j];
          fv->d_grad_api_dmesh[0][2][j] += *esp->SH[i] * bfv->d_grad_phi_dmesh[i][0][2][j];
          fv->d_grad_api_dmesh[1][2][j] += *esp->SH[i] * bfv->d_grad_phi_dmesh[i][1][2][j];
        }
      }
    }
#endif
  } else if (upd->vp[pg->imtrx][SHEAR_RATE] != -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(fv->d_grad_SH_dmesh, 0, siz);
  }

  /*
   * d(grad(F))/dmesh
   */
  if (pd->gv[FILL]) {
    v = FILL;
    vdofs = ei[upd->matrix_index[v]]->dof[v];
#ifdef DO_NO_UNROLL
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(fv->d_grad_F_dmesh, 0, siz);
    for (p = 0; p < dimNonSym; p++) {
      for (b = 0; b < dim; b++) {
        for (j = 0; j < mdofs; j++) {
          for (i = 0; i < vdofs; i++) {
            fv->d_grad_F_dmesh[p][b][j] += *esp->F[i] * bf[v]->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
#else
    for (j = 0; j < mdofs; j++) {
      F_i = *esp->F[0];

      fv->d_grad_F_dmesh[0][0][j] = F_i * bf[v]->d_grad_phi_dmesh[0][0][0][j];
      fv->d_grad_F_dmesh[1][1][j] = F_i * bf[v]->d_grad_phi_dmesh[0][1][1][j];
      fv->d_grad_F_dmesh[1][0][j] = F_i * bf[v]->d_grad_phi_dmesh[0][1][0][j];
      fv->d_grad_F_dmesh[0][1][j] = F_i * bf[v]->d_grad_phi_dmesh[0][0][1][j];

      if (dimNonSym == 3) {
        fv->d_grad_F_dmesh[2][2][j] = F_i * bf[v]->d_grad_phi_dmesh[0][2][2][j];
        fv->d_grad_F_dmesh[2][0][j] = F_i * bf[v]->d_grad_phi_dmesh[0][2][0][j];
        fv->d_grad_F_dmesh[2][1][j] = F_i * bf[v]->d_grad_phi_dmesh[0][2][1][j];
        fv->d_grad_F_dmesh[0][2][j] = F_i * bf[v]->d_grad_phi_dmesh[0][0][2][j];
        fv->d_grad_F_dmesh[1][2][j] = F_i * bf[v]->d_grad_phi_dmesh[0][1][2][j];
      }

      for (i = 1; i < vdofs; i++) {
        F_i = *esp->F[i];

        fv->d_grad_F_dmesh[0][0][j] += F_i * bf[v]->d_grad_phi_dmesh[i][0][0][j];
        fv->d_grad_F_dmesh[1][1][j] += F_i * bf[v]->d_grad_phi_dmesh[i][1][1][j];
        fv->d_grad_F_dmesh[1][0][j] += F_i * bf[v]->d_grad_phi_dmesh[i][1][0][j];
        fv->d_grad_F_dmesh[0][1][j] += F_i * bf[v]->d_grad_phi_dmesh[i][0][1][j];

        if (dimNonSym == 3) {
          fv->d_grad_F_dmesh[2][2][j] += F_i * bf[v]->d_grad_phi_dmesh[i][2][2][j];
          fv->d_grad_F_dmesh[2][0][j] += F_i * bf[v]->d_grad_phi_dmesh[i][2][0][j];
          fv->d_grad_F_dmesh[2][1][j] += F_i * bf[v]->d_grad_phi_dmesh[i][2][1][j];
          fv->d_grad_F_dmesh[0][2][j] += F_i * bf[v]->d_grad_phi_dmesh[i][0][2][j];
          fv->d_grad_F_dmesh[1][2][j] += F_i * bf[v]->d_grad_phi_dmesh[i][1][2][j];
        }
      }
    }
#endif

  } else if (upd->vp[pg->imtrx][FILL] != -1) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(fv->d_grad_F_dmesh, 0, siz);
  }

  /*
   * d(grad(P))/dmesh
   */

  if (pd->gv[PRESSURE]) {
    v = PRESSURE;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];
#ifdef DO_NO_UNROLL
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(fv->d_grad_P_dmesh, 0, siz);
    for (i = 0; i < vdofs; i++) {
      P_i = *esp->P[i];
      for (p = 0; p < dimNonSym; p++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_P_dmesh[p][b][j] += P_i * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
#else
    for (j = 0; j < mdofs; j++) {
      P_i = *esp->P[0];
      fv->d_grad_P_dmesh[0][0][j] = P_i * bfv->d_grad_phi_dmesh[0][0][0][j];
      fv->d_grad_P_dmesh[1][1][j] = P_i * bfv->d_grad_phi_dmesh[0][1][1][j];
      fv->d_grad_P_dmesh[1][0][j] = P_i * bfv->d_grad_phi_dmesh[0][1][0][j];
      fv->d_grad_P_dmesh[0][1][j] = P_i * bfv->d_grad_phi_dmesh[0][0][1][j];

      if (dimNonSym == 3) {
        fv->d_grad_P_dmesh[2][2][j] = P_i * bfv->d_grad_phi_dmesh[0][2][2][j];
        fv->d_grad_P_dmesh[2][0][j] = P_i * bfv->d_grad_phi_dmesh[0][2][0][j];
        fv->d_grad_P_dmesh[2][1][j] = P_i * bfv->d_grad_phi_dmesh[0][2][1][j];
        fv->d_grad_P_dmesh[1][2][j] = P_i * bfv->d_grad_phi_dmesh[0][1][2][j];
        fv->d_grad_P_dmesh[0][2][j] = P_i * bfv->d_grad_phi_dmesh[0][0][2][j];
      }

      for (i = 1; i < vdofs; i++) {
        P_i = *esp->P[i];
        fv->d_grad_P_dmesh[0][0][j] += P_i * bfv->d_grad_phi_dmesh[i][0][0][j];
        fv->d_grad_P_dmesh[1][1][j] += P_i * bfv->d_grad_phi_dmesh[i][1][1][j];
        fv->d_grad_P_dmesh[1][0][j] += P_i * bfv->d_grad_phi_dmesh[i][1][0][j];
        fv->d_grad_P_dmesh[0][1][j] += P_i * bfv->d_grad_phi_dmesh[i][0][1][j];

        if (dimNonSym == 3) {
          fv->d_grad_P_dmesh[2][2][j] += P_i * bfv->d_grad_phi_dmesh[i][2][2][j];
          fv->d_grad_P_dmesh[2][0][j] += P_i * bfv->d_grad_phi_dmesh[i][2][0][j];
          fv->d_grad_P_dmesh[2][1][j] += P_i * bfv->d_grad_phi_dmesh[i][2][1][j];
          fv->d_grad_P_dmesh[1][2][j] += P_i * bfv->d_grad_phi_dmesh[i][1][2][j];
          fv->d_grad_P_dmesh[0][2][j] += P_i * bfv->d_grad_phi_dmesh[i][0][2][j];
        }
      }
    }
#endif
  } else if (upd->vp[pg->imtrx][PRESSURE] != -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(fv->d_grad_P_dmesh, 0, siz);
  }

  /*
   * d(grad(nn))/dmesh
   */

  if (pd->gv[BOND_EVOLUTION]) {
    v = BOND_EVOLUTION;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];
#ifdef DO_NOT_UNROLL
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(fv->d_grad_nn_dmesh, 0, siz);
    for (i = 0; i < vdofs; i++) {
      for (p = 0; p < dim; p++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_nn_dmesh[p][b][j] += *esp->nn[i] * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
#else
    for (j = 0; j < mdofs; j++) {
      fv->d_grad_nn_dmesh[0][0][j] = *esp->nn[0] * bfv->d_grad_phi_dmesh[0][0][0][j];
      fv->d_grad_nn_dmesh[1][1][j] = *esp->nn[0] * bfv->d_grad_phi_dmesh[0][1][1][j];
      fv->d_grad_nn_dmesh[1][0][j] = *esp->nn[0] * bfv->d_grad_phi_dmesh[0][1][0][j];
      fv->d_grad_nn_dmesh[0][1][j] = *esp->nn[0] * bfv->d_grad_phi_dmesh[0][0][1][j];

      if (dimNonSym == 3) {
        fv->d_grad_nn_dmesh[2][2][j] = *esp->nn[0] * bfv->d_grad_phi_dmesh[0][2][2][j];
        fv->d_grad_nn_dmesh[2][0][j] = *esp->nn[0] * bfv->d_grad_phi_dmesh[0][2][0][j];
        fv->d_grad_nn_dmesh[2][1][j] = *esp->nn[0] * bfv->d_grad_phi_dmesh[0][2][1][j];
        fv->d_grad_nn_dmesh[0][2][j] = *esp->nn[0] * bfv->d_grad_phi_dmesh[0][0][2][j];
        fv->d_grad_nn_dmesh[1][2][j] = *esp->nn[0] * bfv->d_grad_phi_dmesh[0][1][2][j];
      }

      for (i = 1; i < vdofs; i++) {
        fv->d_grad_nn_dmesh[0][0][j] += *esp->nn[i] * bfv->d_grad_phi_dmesh[i][0][0][j];
        fv->d_grad_nn_dmesh[1][1][j] += *esp->nn[i] * bfv->d_grad_phi_dmesh[i][1][1][j];
        fv->d_grad_nn_dmesh[1][0][j] += *esp->nn[i] * bfv->d_grad_phi_dmesh[i][1][0][j];
        fv->d_grad_nn_dmesh[0][1][j] += *esp->nn[i] * bfv->d_grad_phi_dmesh[i][0][1][j];

        if (dimNonSym == 3) {
          fv->d_grad_nn_dmesh[2][2][j] += *esp->nn[i] * bfv->d_grad_phi_dmesh[i][2][2][j];
          fv->d_grad_nn_dmesh[2][0][j] += *esp->nn[i] * bfv->d_grad_phi_dmesh[i][2][0][j];
          fv->d_grad_nn_dmesh[2][1][j] += *esp->nn[i] * bfv->d_grad_phi_dmesh[i][2][1][j];
          fv->d_grad_nn_dmesh[0][2][j] += *esp->nn[i] * bfv->d_grad_phi_dmesh[i][0][2][j];
          fv->d_grad_nn_dmesh[1][2][j] += *esp->nn[i] * bfv->d_grad_phi_dmesh[i][1][2][j];
        }
      }
    }
#endif
  }

  /*
   * d(grad(c_w))/dmesh
   */
  if (pd->gv[MASS_FRACTION]) {
    v = MASS_FRACTION;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];
#ifdef DO_NOT_UNROLL
    siz = sizeof(double) * DIM * DIM * MAX_CONC * MDE;
    memset(fv->d_grad_c_dmesh, 0, siz);

    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      for (p = 0; p < dim; p++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            for (i = 0; i < vdofs; i++) {
              fv->d_grad_c_dmesh[p][w][b][j] += *esp->c[w][i] * bfv->d_grad_phi_dmesh[i][p][b][j];
            }
          }
        }
      }
    }
#else
    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      for (j = 0; j < mdofs; j++) {
        fv->d_grad_c_dmesh[0][w][0][j] = *esp->c[w][0] * bfv->d_grad_phi_dmesh[0][0][0][j];
        fv->d_grad_c_dmesh[1][w][1][j] = *esp->c[w][0] * bfv->d_grad_phi_dmesh[0][1][1][j];
        fv->d_grad_c_dmesh[1][w][0][j] = *esp->c[w][0] * bfv->d_grad_phi_dmesh[0][1][0][j];
        fv->d_grad_c_dmesh[0][w][1][j] = *esp->c[w][0] * bfv->d_grad_phi_dmesh[0][0][1][j];

        if (dimNonSym == 3) {
          fv->d_grad_c_dmesh[2][w][2][j] = *esp->c[w][0] * bfv->d_grad_phi_dmesh[0][2][2][j];
          fv->d_grad_c_dmesh[2][w][0][j] = *esp->c[w][0] * bfv->d_grad_phi_dmesh[0][2][0][j];
          fv->d_grad_c_dmesh[2][w][1][j] = *esp->c[w][0] * bfv->d_grad_phi_dmesh[0][2][1][j];
          fv->d_grad_c_dmesh[0][w][2][j] = *esp->c[w][0] * bfv->d_grad_phi_dmesh[0][0][2][j];
          fv->d_grad_c_dmesh[1][w][2][j] = *esp->c[w][0] * bfv->d_grad_phi_dmesh[0][1][2][j];
        }

        for (i = 1; i < vdofs; i++) {
          fv->d_grad_c_dmesh[0][w][0][j] += *esp->c[w][i] * bfv->d_grad_phi_dmesh[i][0][0][j];
          fv->d_grad_c_dmesh[1][w][1][j] += *esp->c[w][i] * bfv->d_grad_phi_dmesh[i][1][1][j];
          fv->d_grad_c_dmesh[1][w][0][j] += *esp->c[w][i] * bfv->d_grad_phi_dmesh[i][1][0][j];
          fv->d_grad_c_dmesh[0][w][1][j] += *esp->c[w][i] * bfv->d_grad_phi_dmesh[i][0][1][j];

          if (dimNonSym == 3) {
            fv->d_grad_c_dmesh[2][w][2][j] += *esp->c[w][i] * bfv->d_grad_phi_dmesh[i][2][2][j];
            fv->d_grad_c_dmesh[2][w][0][j] += *esp->c[w][i] * bfv->d_grad_phi_dmesh[i][2][0][j];
            fv->d_grad_c_dmesh[2][w][1][j] += *esp->c[w][i] * bfv->d_grad_phi_dmesh[i][2][1][j];
            fv->d_grad_c_dmesh[0][w][2][j] += *esp->c[w][i] * bfv->d_grad_phi_dmesh[i][0][2][j];
            fv->d_grad_c_dmesh[1][w][2][j] += *esp->c[w][i] * bfv->d_grad_phi_dmesh[i][1][2][j];
          }
        }
      }
    }
#endif
  } else if (upd->vp[pg->imtrx][MASS_FRACTION] != -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * MAX_CONC * MDE;
    memset(fv->d_grad_c_dmesh, 0, siz);
  }

  /*
   * d(grad(porous_media_variables))/dmesh
   */
  if (pd->gv[POR_LIQ_PRES]) {
    v = POR_LIQ_PRES;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];
#ifdef DO_NO_UNROLL
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(fv->d_grad_p_liq_dmesh, 0, siz);

    for (p = 0; p < dim; p++) {
      for (b = 0; b < dim; b++) {
        for (j = 0; j < mdofs; j++) {
          for (i = 0; i < vdofs; i++) {
            fv->d_grad_p_liq_dmesh[p][b][j] += *esp->p_liq[i] * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
#else
    for (j = 0; j < mdofs; j++) {
      P_i = *esp->p_liq[0];
      fv->d_grad_p_liq_dmesh[0][0][j] = P_i * bfv->d_grad_phi_dmesh[0][0][0][j];
      fv->d_grad_p_liq_dmesh[1][1][j] = P_i * bfv->d_grad_phi_dmesh[0][1][1][j];
      fv->d_grad_p_liq_dmesh[0][1][j] = P_i * bfv->d_grad_phi_dmesh[0][0][1][j];
      fv->d_grad_p_liq_dmesh[1][0][j] = P_i * bfv->d_grad_phi_dmesh[0][1][0][j];

      if (dim == 3) {
        fv->d_grad_p_liq_dmesh[2][2][j] = P_i * bfv->d_grad_phi_dmesh[0][2][2][j];
        fv->d_grad_p_liq_dmesh[2][0][j] = P_i * bfv->d_grad_phi_dmesh[0][2][0][j];
        fv->d_grad_p_liq_dmesh[2][1][j] = P_i * bfv->d_grad_phi_dmesh[0][2][1][j];
        fv->d_grad_p_liq_dmesh[0][2][j] = P_i * bfv->d_grad_phi_dmesh[0][0][2][j];
        fv->d_grad_p_liq_dmesh[1][2][j] = P_i * bfv->d_grad_phi_dmesh[0][1][2][j];
      }

      for (i = 1; i < vdofs; i++) {
        P_i = *esp->p_liq[i];
        /* fv->d_grad_p_liq_dmesh[p] [b][j] +=  P_i * bfv->d_grad_phi_dmesh[i][p] [b][j];*/
        fv->d_grad_p_liq_dmesh[0][0][j] += P_i * bfv->d_grad_phi_dmesh[i][0][0][j];
        fv->d_grad_p_liq_dmesh[1][1][j] += P_i * bfv->d_grad_phi_dmesh[i][1][1][j];
        fv->d_grad_p_liq_dmesh[0][1][j] += P_i * bfv->d_grad_phi_dmesh[i][0][1][j];
        fv->d_grad_p_liq_dmesh[1][0][j] += P_i * bfv->d_grad_phi_dmesh[i][1][0][j];

        if (dim == 3) {
          fv->d_grad_p_liq_dmesh[2][2][j] += P_i * bfv->d_grad_phi_dmesh[i][2][2][j];
          fv->d_grad_p_liq_dmesh[2][0][j] += P_i * bfv->d_grad_phi_dmesh[i][2][0][j];
          fv->d_grad_p_liq_dmesh[2][1][j] += P_i * bfv->d_grad_phi_dmesh[i][2][1][j];
          fv->d_grad_p_liq_dmesh[0][2][j] += P_i * bfv->d_grad_phi_dmesh[i][0][2][j];
          fv->d_grad_p_liq_dmesh[1][2][j] += P_i * bfv->d_grad_phi_dmesh[i][1][2][j];
        }
      }
    }
#endif

  } else if (upd->vp[pg->imtrx][POR_LIQ_PRES] != -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(fv->d_grad_p_liq_dmesh, 0, siz);
  }

  if (pd->gv[POR_GAS_PRES]) {
    v = POR_GAS_PRES;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];
#ifdef DO_NOT_UNROLL
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(fv->d_grad_p_gas_dmesh, 0, siz);

    for (p = 0; p < dim; p++) {
      for (b = 0; b < dim; b++) {
        for (j = 0; j < mdofs; j++) {
          for (i = 0; i < vdofs; i++) {
            fv->d_grad_p_gas_dmesh[p][b][j] += *esp->p_gas[i] * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
#else
    for (j = 0; j < mdofs; j++) {
      fv->d_grad_p_gas_dmesh[0][0][j] = *esp->p_gas[0] * bfv->d_grad_phi_dmesh[0][0][0][j];
      fv->d_grad_p_gas_dmesh[1][1][j] = *esp->p_gas[0] * bfv->d_grad_phi_dmesh[0][1][1][j];
      fv->d_grad_p_gas_dmesh[1][0][j] = *esp->p_gas[0] * bfv->d_grad_phi_dmesh[0][1][0][j];
      fv->d_grad_p_gas_dmesh[0][1][j] = *esp->p_gas[0] * bfv->d_grad_phi_dmesh[0][0][1][j];

      if (dimNonSym == 3) {
        fv->d_grad_p_gas_dmesh[2][2][j] = *esp->p_gas[0] * bfv->d_grad_phi_dmesh[0][2][2][j];
        fv->d_grad_p_gas_dmesh[2][0][j] = *esp->p_gas[0] * bfv->d_grad_phi_dmesh[0][2][0][j];
        fv->d_grad_p_gas_dmesh[2][1][j] = *esp->p_gas[0] * bfv->d_grad_phi_dmesh[0][2][1][j];
        fv->d_grad_p_gas_dmesh[0][2][j] = *esp->p_gas[0] * bfv->d_grad_phi_dmesh[0][0][2][j];
        fv->d_grad_p_gas_dmesh[1][2][j] = *esp->p_gas[0] * bfv->d_grad_phi_dmesh[0][1][2][j];
      }

      for (i = 1; i < vdofs; i++) {
        fv->d_grad_p_gas_dmesh[0][0][j] += *esp->p_gas[i] * bfv->d_grad_phi_dmesh[i][0][0][j];
        fv->d_grad_p_gas_dmesh[1][1][j] += *esp->p_gas[i] * bfv->d_grad_phi_dmesh[i][1][1][j];
        fv->d_grad_p_gas_dmesh[1][0][j] += *esp->p_gas[i] * bfv->d_grad_phi_dmesh[i][1][0][j];
        fv->d_grad_p_gas_dmesh[0][1][j] += *esp->p_gas[i] * bfv->d_grad_phi_dmesh[i][0][1][j];

        if (dimNonSym == 3) {
          fv->d_grad_p_gas_dmesh[2][2][j] += *esp->p_gas[i] * bfv->d_grad_phi_dmesh[i][2][2][j];
          fv->d_grad_p_gas_dmesh[2][0][j] += *esp->p_gas[i] * bfv->d_grad_phi_dmesh[i][2][0][j];
          fv->d_grad_p_gas_dmesh[2][1][j] += *esp->p_gas[i] * bfv->d_grad_phi_dmesh[i][2][1][j];
          fv->d_grad_p_gas_dmesh[0][2][j] += *esp->p_gas[i] * bfv->d_grad_phi_dmesh[i][0][2][j];
          fv->d_grad_p_gas_dmesh[1][2][j] += *esp->p_gas[i] * bfv->d_grad_phi_dmesh[i][1][2][j];
        }
      }
    }
#endif
  } else if (upd->vp[pg->imtrx][POR_GAS_PRES] != -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(fv->d_grad_p_gas_dmesh, 0, siz);
  }

  if (pd->gv[POR_POROSITY]) {
    v = POR_POROSITY;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];
#ifdef DO_NOT_UNROLL
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(fv->d_grad_porosity_dmesh, 0, siz);

    for (p = 0; p < dim; p++) {
      for (b = 0; b < dim; b++) {
        for (j = 0; j < mdofs; j++) {
          for (i = 0; i < vdofs; i++) {
            fv->d_grad_porosity_dmesh[p][b][j] +=
                *esp->porosity[i] * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
#else
    for (j = 0; j < mdofs; j++) {
      fv->d_grad_porosity_dmesh[0][0][j] = *esp->porosity[0] * bfv->d_grad_phi_dmesh[0][0][0][j];
      fv->d_grad_porosity_dmesh[1][1][j] = *esp->porosity[0] * bfv->d_grad_phi_dmesh[0][1][1][j];
      fv->d_grad_porosity_dmesh[1][0][j] = *esp->porosity[0] * bfv->d_grad_phi_dmesh[0][1][0][j];
      fv->d_grad_porosity_dmesh[0][1][j] = *esp->porosity[0] * bfv->d_grad_phi_dmesh[0][0][1][j];

      if (dimNonSym == 3) {
        fv->d_grad_porosity_dmesh[2][2][j] = *esp->porosity[0] * bfv->d_grad_phi_dmesh[0][2][2][j];
        fv->d_grad_porosity_dmesh[2][0][j] = *esp->porosity[0] * bfv->d_grad_phi_dmesh[0][2][0][j];
        fv->d_grad_porosity_dmesh[2][1][j] = *esp->porosity[0] * bfv->d_grad_phi_dmesh[0][2][1][j];
        fv->d_grad_porosity_dmesh[0][2][j] = *esp->porosity[0] * bfv->d_grad_phi_dmesh[0][0][2][j];
        fv->d_grad_porosity_dmesh[1][2][j] = *esp->porosity[0] * bfv->d_grad_phi_dmesh[0][1][2][j];
      }

      for (i = 1; i < vdofs; i++) {
        fv->d_grad_porosity_dmesh[0][0][j] += *esp->porosity[i] * bfv->d_grad_phi_dmesh[i][0][0][j];
        fv->d_grad_porosity_dmesh[1][1][j] += *esp->porosity[i] * bfv->d_grad_phi_dmesh[i][1][1][j];
        fv->d_grad_porosity_dmesh[1][0][j] += *esp->porosity[i] * bfv->d_grad_phi_dmesh[i][1][0][j];
        fv->d_grad_porosity_dmesh[0][1][j] += *esp->porosity[i] * bfv->d_grad_phi_dmesh[i][0][1][j];

        if (dimNonSym == 3) {
          fv->d_grad_porosity_dmesh[2][2][j] +=
              *esp->porosity[i] * bfv->d_grad_phi_dmesh[i][2][2][j];
          fv->d_grad_porosity_dmesh[2][0][j] +=
              *esp->porosity[i] * bfv->d_grad_phi_dmesh[i][2][0][j];
          fv->d_grad_porosity_dmesh[2][1][j] +=
              *esp->porosity[i] * bfv->d_grad_phi_dmesh[i][2][1][j];
          fv->d_grad_porosity_dmesh[0][2][j] +=
              *esp->porosity[i] * bfv->d_grad_phi_dmesh[i][0][2][j];
          fv->d_grad_porosity_dmesh[1][2][j] +=
              *esp->porosity[i] * bfv->d_grad_phi_dmesh[i][1][2][j];
        }
      }
    }
#endif
  } else if (upd->vp[pg->imtrx][POR_POROSITY] != -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(fv->d_grad_porosity_dmesh, 0, siz);
  }

  /*
   * d( grad(v))/dmesh
   */
  if (pd->gv[VELOCITY1]) {
    v = VELOCITY1;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];

#ifdef DO_NO_UNROLL
    siz = sizeof(double) * DIM * VIM * DIM * MDE;
    memset(fv->d_grad_v_dmesh, 0, siz);

    for (r = 0; r < WIM; r++) {
      for (i = 0; i < vdofs; i++) {
        v_ri = *esp->v[r][i];
        for (p = 0; p < VIM; p++) {
          for (q = 0; q < VIM; q++) {
            for (b = 0; b < dim; b++) {
              for (j = 0; j < mdofs; j++) {
                fv->d_grad_v_dmesh[p][q][b][j] += v_ri * bfv->d_grad_phi_e_dmesh[i][r][p][q][b][j];
              }
            }
          }
        }
      }
    }
#else

    v_ri = *esp->v[0][0];

    for (b = 0; b < dim; b++) {
      for (j = 0; j < mdofs; j++) {
        fv->d_grad_v_dmesh[0][0][b][j] = 0.0;
        fv->d_grad_v_dmesh[1][1][b][j] = 0.0;
        fv->d_grad_v_dmesh[1][0][b][j] = 0.0;
        fv->d_grad_v_dmesh[0][1][b][j] = 0.0;
        if (VIMis3) {
          fv->d_grad_v_dmesh[2][2][b][j] = 0.0;
          fv->d_grad_v_dmesh[2][0][b][j] = 0.0;
          fv->d_grad_v_dmesh[2][1][b][j] = 0.0;
          fv->d_grad_v_dmesh[0][2][b][j] = 0.0;
          fv->d_grad_v_dmesh[1][2][b][j] = 0.0;
        }
      }
    }

    for (r = 0; r < WIM; r++) {
      for (i = 0; i < vdofs; i++) {
        v_ri = *esp->v[r][i];
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_v_dmesh[0][0][b][j] += v_ri * bfv->d_grad_phi_e_dmesh[i][r][0][0][b][j];
            fv->d_grad_v_dmesh[1][1][b][j] += v_ri * bfv->d_grad_phi_e_dmesh[i][r][1][1][b][j];
            fv->d_grad_v_dmesh[1][0][b][j] += v_ri * bfv->d_grad_phi_e_dmesh[i][r][1][0][b][j];
            fv->d_grad_v_dmesh[0][1][b][j] += v_ri * bfv->d_grad_phi_e_dmesh[i][r][0][1][b][j];

            if (VIMis3) {
              fv->d_grad_v_dmesh[2][2][b][j] += v_ri * bfv->d_grad_phi_e_dmesh[i][r][2][2][b][j];
              fv->d_grad_v_dmesh[2][0][b][j] += v_ri * bfv->d_grad_phi_e_dmesh[i][r][2][0][b][j];
              fv->d_grad_v_dmesh[2][1][b][j] += v_ri * bfv->d_grad_phi_e_dmesh[i][r][2][1][b][j];
              fv->d_grad_v_dmesh[0][2][b][j] += v_ri * bfv->d_grad_phi_e_dmesh[i][r][0][2][b][j];
              fv->d_grad_v_dmesh[1][2][b][j] += v_ri * bfv->d_grad_phi_e_dmesh[i][r][1][2][b][j];
            }
          }
        }
      }
    }
#endif
    /*
     * d( div(v) )/dmesh
     */
    /*
     * There is currently no need for div(pv), so this isn't cloned.
     */
    siz = sizeof(double) * DIM * MDE;
    memset(fv->d_div_v_dmesh, 0, siz);

    for (b = 0; b < dim; b++) {
      for (j = 0; j < mdofs; j++) {
        for (p = 0; p < VIM; p++) {
          fv->d_div_v_dmesh[b][j] += fv->d_grad_v_dmesh[p][p][b][j];
        }
      }
    }
  } else if (upd->vp[pg->imtrx][VELOCITY1] != -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * DIM * MDE;
    memset(fv->d_grad_v_dmesh, 0, siz);
    siz = sizeof(double) * DIM * MDE;
    memset(fv->d_div_v_dmesh, 0, siz);
  }

  /*
   * d( grad(pv))/dmesh
   */
  if (pd->gv[PVELOCITY1]) {
    v = PVELOCITY1;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];
    siz = sizeof(double) * DIM * DIM * DIM * MDE;
    memset(fv->d_grad_pv_dmesh, 0, siz);

    for (r = 0; r < WIM; r++) {
      for (i = 0; i < vdofs; i++) {
        v_ri = *esp->pv[r][i];
        for (p = 0; p < VIM; p++) {
          for (q = 0; q < VIM; q++) {
            for (b = 0; b < dim; b++) {
              for (j = 0; j < mdofs; j++) {
                fv->d_grad_pv_dmesh[p][q][b][j] += v_ri * bfv->d_grad_phi_e_dmesh[i][r][p][q][b][j];
              }
            }
          }
        }
      }
    }
  } else if (upd->vp[pg->imtrx][PVELOCITY1] != -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * DIM * MDE;
    memset(fv->d_grad_pv_dmesh, 0, siz);
  }

  /*
   * d( grad(ext_v))/dmesh
   */
  if (pd->gv[EXT_VELOCITY]) {
    v = EXT_VELOCITY;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(fv->d_grad_ext_v_dmesh, 0, siz);

    for (i = 0; i < vdofs; i++) {
      v_ri = *esp->ext_v[i];
      for (p = 0; p < VIM; p++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_ext_v_dmesh[p][b][j] += v_ri * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
  } else if (upd->vp[pg->imtrx][EXT_VELOCITY] != -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(fv->d_grad_ext_v_dmesh, 0, siz);
  }

  /*
   * d(grad(n))/dmesh  -- ls normal vector wrt mesh
   *
   * This is the gradient of the surface normal with respect to the mesh displacements.
   * This is a tensor quantity.
   *     d_grad_n_dmesh[p][q] [b][j]  -
   *
   *           p is the first coordinate of the tensor (nominally the derivative)
   *           q is the second cordinate of the tensor (nominally the direction of n)
   *           b is the mesh displacement coordinate.
   *           j is the mesh displacement degree of freedom in the current element
   *             (note for SHELL_NORMAL_1 this element refers to the shell element dof)
   *
   */
  if (pd->gv[NORMAL1] || pd->gv[SHELL_NORMAL1]) {
    if (pd->gv[NORMAL1])
      v = NORMAL1;
    if (pd->gv[SHELL_NORMAL1])
      v = SHELL_NORMAL1;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];
    siz = sizeof(double) * DIM * DIM * DIM * MDE;
    memset(fv->d_grad_n_dmesh, 0, siz);
    siz = sizeof(double) * DIM * MDE;
    memset(fv->d_div_n_dmesh, 0, siz);

    if (v == NORMAL1) {
      // HKM -> looks good according to d_grad_v_dmesh[][][][]
      for (p = 0; p < VIM; p++) {
        for (q = 0; q < VIM; q++) {
          for (b = 0; b < dim; b++) {
            for (j = 0; j < mdofs; j++) {
              for (r = 0; r < dim; r++) {
                for (i = 0; i < vdofs; i++) {
                  v_ri = *esp->n[r][i];
                  fv->d_grad_n_dmesh[p][q][b][j] +=
                      v_ri * bfv->d_grad_phi_e_dmesh[i][r][p][q][b][j];
                }
              }
            }
          }
        }
      }
    }

    else if (v == SHELL_NORMAL1) {
      for (r = 0; r < dim; r++) {
        for (i = 0; i < vdofs; i++) {
          v_ri = *esp->n[r][i];
          for (p = 0; p < VIM; p++) {
            for (b = 0; b < dim; b++) {
              for (j = 0; j < mdofs; j++) {
                fv->d_grad_n_dmesh[p][r][b][j] += v_ri * bfv->d_grad_phi_dmesh[i][p][b][j];
              }
            }
          }
        }
      }
    }

    for (b = 0; b < dim; b++) {
      for (j = 0; j < mdofs; j++) {
        for (p = 0; p < VIM; p++) {
          fv->d_div_n_dmesh[b][j] += fv->d_grad_n_dmesh[p][p][b][j];
        }
      }
    }
  } else if (((upd->vp[pg->imtrx][NORMAL1] != -1) && (upd->vp[pg->imtrx][SHELL_NORMAL1] != -1)) &&
             (okToZero)) {
    siz = sizeof(double) * DIM * VIM * DIM * MDE;
    memset(fv->d_grad_n_dmesh, 0, siz);
    siz = sizeof(double) * DIM * MDE;
    memset(fv->d_div_n_dmesh, 0, siz);
  }

  /*
   * d( grad(E_field))/dmesh
   */

  if (pd->gv[EFIELD1]) {
    v = EFIELD1;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];

    siz = sizeof(double) * DIM * DIM * DIM * MDE;
    memset(fv->d_grad_E_field_dmesh, 0, siz);

    for (r = 0; r < WIM; r++) {
      for (i = 0; i < vdofs; i++) {
        v_ri = *esp->E_field[r][i];
        for (p = 0; p < VIM; p++) {
          for (q = 0; q < VIM; q++) {
            for (b = 0; b < dim; b++) {
              for (j = 0; j < mdofs; j++) {
                fv->d_grad_E_field_dmesh[p][q][b][j] +=
                    v_ri * bfv->d_grad_phi_e_dmesh[i][r][p][q][b][j];
              }
            }
          }
        }
      }
    }
  } else if (upd->vp[pg->imtrx][EFIELD1] != -1 && okToZero) {
    siz = sizeof(double) * DIM * VIM * DIM * MDE;
    memset(fv->d_grad_E_field_dmesh, 0, siz);
  }

  /*
   * d( grad(em_er))/dmesh
   */
  if (pd->v[pg->imtrx][EM_E1_REAL] || pd->v[pg->imtrx][EM_E2_REAL] ||
      pd->v[pg->imtrx][EM_E3_REAL]) {
    v = EM_E1_REAL;
    bfv = bf[v];
    vdofs = ei[pg->imtrx]->dof[v];

#ifdef DO_NO_UNROLL
    siz = sizeof(double) * DIM * VIM * DIM * MDE;
    memset(fv->d_grad_em_er_dmesh, 0, siz);

    for (r = 0; r < WIM; r++) {
      for (i = 0; i < vdofs; i++) {
        em_er_ri = *esp->em_er[r][i];
        for (p = 0; p < VIM; p++) {
          for (q = 0; q < VIM; q++) {
            for (b = 0; b < dim; b++) {
              for (j = 0; j < mdofs; j++) {
                fv->d_grad_em_er_dmesh[p][q][b][j] +=
                    em_er_ri * bfv->d_grad_phi_e_dmesh[i][r][p][q][b][j];
              }
            }
          }
        }
      }
    }
#else

    em_er_ri = *esp->em_er[0][0];

    for (b = 0; b < dim; b++) {
      for (j = 0; j < mdofs; j++) {
        fv->d_grad_em_er_dmesh[0][0][b][j] = 0.0;
        fv->d_grad_em_er_dmesh[1][1][b][j] = 0.0;
        fv->d_grad_em_er_dmesh[1][0][b][j] = 0.0;
        fv->d_grad_em_er_dmesh[0][1][b][j] = 0.0;
        if (VIMis3) {
          fv->d_grad_em_er_dmesh[2][2][b][j] = 0.0;
          fv->d_grad_em_er_dmesh[2][0][b][j] = 0.0;
          fv->d_grad_em_er_dmesh[2][1][b][j] = 0.0;
          fv->d_grad_em_er_dmesh[0][2][b][j] = 0.0;
          fv->d_grad_em_er_dmesh[1][2][b][j] = 0.0;
        }
      }
    }

    for (r = 0; r < WIM; r++) {
      for (i = 0; i < vdofs; i++) {
        em_er_ri = *esp->em_er[r][i];
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_em_er_dmesh[0][0][b][j] +=
                em_er_ri * bfv->d_grad_phi_e_dmesh[i][r][0][0][b][j];
            fv->d_grad_em_er_dmesh[1][1][b][j] +=
                em_er_ri * bfv->d_grad_phi_e_dmesh[i][r][1][1][b][j];
            fv->d_grad_em_er_dmesh[1][0][b][j] +=
                em_er_ri * bfv->d_grad_phi_e_dmesh[i][r][1][0][b][j];
            fv->d_grad_em_er_dmesh[0][1][b][j] +=
                em_er_ri * bfv->d_grad_phi_e_dmesh[i][r][0][1][b][j];

            if (VIMis3) {
              fv->d_grad_em_er_dmesh[2][2][b][j] +=
                  em_er_ri * bfv->d_grad_phi_e_dmesh[i][r][2][2][b][j];
              fv->d_grad_em_er_dmesh[2][0][b][j] +=
                  em_er_ri * bfv->d_grad_phi_e_dmesh[i][r][2][0][b][j];
              fv->d_grad_em_er_dmesh[2][1][b][j] +=
                  em_er_ri * bfv->d_grad_phi_e_dmesh[i][r][2][1][b][j];
              fv->d_grad_em_er_dmesh[0][2][b][j] +=
                  em_er_ri * bfv->d_grad_phi_e_dmesh[i][r][0][2][b][j];
              fv->d_grad_em_er_dmesh[1][2][b][j] +=
                  em_er_ri * bfv->d_grad_phi_e_dmesh[i][r][1][2][b][j];
            }
          }
        }
      }
    }
#endif

  } else if (upd->vp[pg->imtrx][EM_E1_REAL] == -1 && upd->vp[pg->imtrx][EM_E2_REAL] == -1 &&
             upd->vp[pg->imtrx][EM_E3_REAL] == -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * DIM * MDE;
    memset(fv->d_grad_em_er_dmesh, 0, siz);
  }

  /*
   * d( grad(em_ei))/dmesh
   */
  if (pd->v[pg->imtrx][EM_E1_IMAG] || pd->v[pg->imtrx][EM_E2_IMAG] ||
      pd->v[pg->imtrx][EM_E3_IMAG]) {
    v = EM_E1_IMAG;
    bfv = bf[v];
    vdofs = ei[pg->imtrx]->dof[v];

#ifdef DO_NO_UNROLL
    siz = sizeof(double) * DIM * VIM * DIM * MDE;
    memset(fv->d_grad_em_ei_dmesh, 0, siz);

    for (r = 0; r < WIM; r++) {
      for (i = 0; i < vdofs; i++) {
        em_ei_ri = *esp->em_ei[r][i];
        for (p = 0; p < VIM; p++) {
          for (q = 0; q < VIM; q++) {
            for (b = 0; b < dim; b++) {
              for (j = 0; j < mdofs; j++) {
                fv->d_grad_em_ei_dmesh[p][q][b][j] +=
                    em_ei_ri * bfv->d_grad_phi_e_dmesh[i][r][p][q][b][j];
              }
            }
          }
        }
      }
    }
#else

    em_ei_ri = *esp->em_ei[0][0];

    for (b = 0; b < dim; b++) {
      for (j = 0; j < mdofs; j++) {
        fv->d_grad_em_ei_dmesh[0][0][b][j] = 0.0;
        fv->d_grad_em_ei_dmesh[1][1][b][j] = 0.0;
        fv->d_grad_em_ei_dmesh[1][0][b][j] = 0.0;
        fv->d_grad_em_ei_dmesh[0][1][b][j] = 0.0;
        if (VIMis3) {
          fv->d_grad_em_ei_dmesh[2][2][b][j] = 0.0;
          fv->d_grad_em_ei_dmesh[2][0][b][j] = 0.0;
          fv->d_grad_em_ei_dmesh[2][1][b][j] = 0.0;
          fv->d_grad_em_ei_dmesh[0][2][b][j] = 0.0;
          fv->d_grad_em_ei_dmesh[1][2][b][j] = 0.0;
        }
      }
    }

    for (r = 0; r < WIM; r++) {
      for (i = 0; i < vdofs; i++) {
        em_ei_ri = *esp->em_ei[r][i];
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_em_ei_dmesh[0][0][b][j] +=
                em_ei_ri * bfv->d_grad_phi_e_dmesh[i][r][0][0][b][j];
            fv->d_grad_em_ei_dmesh[1][1][b][j] +=
                em_ei_ri * bfv->d_grad_phi_e_dmesh[i][r][1][1][b][j];
            fv->d_grad_em_ei_dmesh[1][0][b][j] +=
                em_ei_ri * bfv->d_grad_phi_e_dmesh[i][r][1][0][b][j];
            fv->d_grad_em_ei_dmesh[0][1][b][j] +=
                em_ei_ri * bfv->d_grad_phi_e_dmesh[i][r][0][1][b][j];

            if (VIMis3) {
              fv->d_grad_em_ei_dmesh[2][2][b][j] +=
                  em_ei_ri * bfv->d_grad_phi_e_dmesh[i][r][2][2][b][j];
              fv->d_grad_em_ei_dmesh[2][0][b][j] +=
                  em_ei_ri * bfv->d_grad_phi_e_dmesh[i][r][2][0][b][j];
              fv->d_grad_em_ei_dmesh[2][1][b][j] +=
                  em_ei_ri * bfv->d_grad_phi_e_dmesh[i][r][2][1][b][j];
              fv->d_grad_em_ei_dmesh[0][2][b][j] +=
                  em_ei_ri * bfv->d_grad_phi_e_dmesh[i][r][0][2][b][j];
              fv->d_grad_em_ei_dmesh[1][2][b][j] +=
                  em_ei_ri * bfv->d_grad_phi_e_dmesh[i][r][1][2][b][j];
            }
          }
        }
      }
    }
#endif

  } else if (upd->vp[pg->imtrx][EM_E1_IMAG] == -1 && upd->vp[pg->imtrx][EM_E2_IMAG] == -1 &&
             upd->vp[pg->imtrx][EM_E3_IMAG] == -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * DIM * MDE;
    memset(fv->d_grad_em_ei_dmesh, 0, siz);
  }

  /*
   * d( grad(em_hr))/dmesh
   */
  if (pd->v[pg->imtrx][EM_H1_REAL] || pd->v[pg->imtrx][EM_H2_REAL] ||
      pd->v[pg->imtrx][EM_H3_REAL]) {
    v = EM_H1_REAL;
    bfv = bf[v];
    vdofs = ei[pg->imtrx]->dof[v];

#ifdef DO_NO_UNROLL
    siz = sizeof(double) * DIM * VIM * DIM * MDE;
    memset(fv->d_grad_em_hr_dmesh, 0, siz);

    for (r = 0; r < WIM; r++) {
      for (i = 0; i < vdofs; i++) {
        em_hr_ri = *esp->em_hr[r][i];
        for (p = 0; p < VIM; p++) {
          for (q = 0; q < VIM; q++) {
            for (b = 0; b < dim; b++) {
              for (j = 0; j < mdofs; j++) {
                fv->d_grad_em_hr_dmesh[p][q][b][j] +=
                    em_hr_ri * bfv->d_grad_phi_e_dmesh[i][r][p][q][b][j];
              }
            }
          }
        }
      }
    }
#else

    em_hr_ri = *esp->em_hr[0][0];

    for (b = 0; b < dim; b++) {
      for (j = 0; j < mdofs; j++) {
        fv->d_grad_em_hr_dmesh[0][0][b][j] = 0.0;
        fv->d_grad_em_hr_dmesh[1][1][b][j] = 0.0;
        fv->d_grad_em_hr_dmesh[1][0][b][j] = 0.0;
        fv->d_grad_em_hr_dmesh[0][1][b][j] = 0.0;
        if (VIMis3) {
          fv->d_grad_em_hr_dmesh[2][2][b][j] = 0.0;
          fv->d_grad_em_hr_dmesh[2][0][b][j] = 0.0;
          fv->d_grad_em_hr_dmesh[2][1][b][j] = 0.0;
          fv->d_grad_em_hr_dmesh[0][2][b][j] = 0.0;
          fv->d_grad_em_hr_dmesh[1][2][b][j] = 0.0;
        }
      }
    }

    for (r = 0; r < WIM; r++) {
      for (i = 0; i < vdofs; i++) {
        em_hr_ri = *esp->em_hr[r][i];
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_em_hr_dmesh[0][0][b][j] +=
                em_hr_ri * bfv->d_grad_phi_e_dmesh[i][r][0][0][b][j];
            fv->d_grad_em_hr_dmesh[1][1][b][j] +=
                em_hr_ri * bfv->d_grad_phi_e_dmesh[i][r][1][1][b][j];
            fv->d_grad_em_hr_dmesh[1][0][b][j] +=
                em_hr_ri * bfv->d_grad_phi_e_dmesh[i][r][1][0][b][j];
            fv->d_grad_em_hr_dmesh[0][1][b][j] +=
                em_hr_ri * bfv->d_grad_phi_e_dmesh[i][r][0][1][b][j];

            if (VIMis3) {
              fv->d_grad_em_hr_dmesh[2][2][b][j] +=
                  em_hr_ri * bfv->d_grad_phi_e_dmesh[i][r][2][2][b][j];
              fv->d_grad_em_hr_dmesh[2][0][b][j] +=
                  em_hr_ri * bfv->d_grad_phi_e_dmesh[i][r][2][0][b][j];
              fv->d_grad_em_hr_dmesh[2][1][b][j] +=
                  em_hr_ri * bfv->d_grad_phi_e_dmesh[i][r][2][1][b][j];
              fv->d_grad_em_hr_dmesh[0][2][b][j] +=
                  em_hr_ri * bfv->d_grad_phi_e_dmesh[i][r][0][2][b][j];
              fv->d_grad_em_hr_dmesh[1][2][b][j] +=
                  em_hr_ri * bfv->d_grad_phi_e_dmesh[i][r][1][2][b][j];
            }
          }
        }
      }
    }
#endif

  } else if (upd->vp[pg->imtrx][EM_H1_REAL] == -1 && upd->vp[pg->imtrx][EM_H2_REAL] == -1 &&
             upd->vp[pg->imtrx][EM_H3_REAL] == -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * DIM * MDE;
    memset(fv->d_grad_em_hr_dmesh, 0, siz);
  }

  /*
   * d( grad(em_hi))/dmesh
   */
  if (pd->v[pg->imtrx][EM_H1_IMAG] || pd->v[pg->imtrx][EM_H2_IMAG] ||
      pd->v[pg->imtrx][EM_H3_IMAG]) {
    v = EM_H1_IMAG;
    bfv = bf[v];
    vdofs = ei[pg->imtrx]->dof[v];

#ifdef DO_NO_UNROLL
    siz = sizeof(double) * DIM * VIM * DIM * MDE;
    memset(fv->d_grad_em_hi_dmesh, 0, siz);

    for (r = 0; r < WIM; r++) {
      for (i = 0; i < vdofs; i++) {
        em_hi_ri = *esp->em_hi[r][i];
        for (p = 0; p < VIM; p++) {
          for (q = 0; q < VIM; q++) {
            for (b = 0; b < dim; b++) {
              for (j = 0; j < mdofs; j++) {
                fv->d_grad_em_hi_dmesh[p][q][b][j] +=
                    em_hi_ri * bfv->d_grad_phi_e_dmesh[i][r][p][q][b][j];
              }
            }
          }
        }
      }
    }
#else

    em_hi_ri = *esp->em_hi[0][0];

    for (b = 0; b < dim; b++) {
      for (j = 0; j < mdofs; j++) {
        fv->d_grad_em_hi_dmesh[0][0][b][j] = 0.0;
        fv->d_grad_em_hi_dmesh[1][1][b][j] = 0.0;
        fv->d_grad_em_hi_dmesh[1][0][b][j] = 0.0;
        fv->d_grad_em_hi_dmesh[0][1][b][j] = 0.0;
        if (VIMis3) {
          fv->d_grad_em_hi_dmesh[2][2][b][j] = 0.0;
          fv->d_grad_em_hi_dmesh[2][0][b][j] = 0.0;
          fv->d_grad_em_hi_dmesh[2][1][b][j] = 0.0;
          fv->d_grad_em_hi_dmesh[0][2][b][j] = 0.0;
          fv->d_grad_em_hi_dmesh[1][2][b][j] = 0.0;
        }
      }
    }

    for (r = 0; r < WIM; r++) {
      for (i = 0; i < vdofs; i++) {
        em_hi_ri = *esp->em_hi[r][i];
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_grad_em_hi_dmesh[0][0][b][j] +=
                em_hi_ri * bfv->d_grad_phi_e_dmesh[i][r][0][0][b][j];
            fv->d_grad_em_hi_dmesh[1][1][b][j] +=
                em_hi_ri * bfv->d_grad_phi_e_dmesh[i][r][1][1][b][j];
            fv->d_grad_em_hi_dmesh[1][0][b][j] +=
                em_hi_ri * bfv->d_grad_phi_e_dmesh[i][r][1][0][b][j];
            fv->d_grad_em_hi_dmesh[0][1][b][j] +=
                em_hi_ri * bfv->d_grad_phi_e_dmesh[i][r][0][1][b][j];

            if (VIMis3) {
              fv->d_grad_em_hi_dmesh[2][2][b][j] +=
                  em_hi_ri * bfv->d_grad_phi_e_dmesh[i][r][2][2][b][j];
              fv->d_grad_em_hi_dmesh[2][0][b][j] +=
                  em_hi_ri * bfv->d_grad_phi_e_dmesh[i][r][2][0][b][j];
              fv->d_grad_em_hi_dmesh[2][1][b][j] +=
                  em_hi_ri * bfv->d_grad_phi_e_dmesh[i][r][2][1][b][j];
              fv->d_grad_em_hi_dmesh[0][2][b][j] +=
                  em_hi_ri * bfv->d_grad_phi_e_dmesh[i][r][0][2][b][j];
              fv->d_grad_em_hi_dmesh[1][2][b][j] +=
                  em_hi_ri * bfv->d_grad_phi_e_dmesh[i][r][1][2][b][j];
            }
          }
        }
      }
    }
#endif

  } else if (upd->vp[pg->imtrx][EM_H1_IMAG] == -1 && upd->vp[pg->imtrx][EM_H2_IMAG] == -1 &&
             upd->vp[pg->imtrx][EM_H3_IMAG] == -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * DIM * MDE;
    memset(fv->d_grad_em_hi_dmesh, 0, siz);
  }
  /*
   * d( grad(d) )/dmesh
   * d( grad(d_dot)/dmesh
   */

  /* Yes, I know it should be each "+p", */
  /* but this way we can collect three */
  /* components of vectors and tensors for */
  /* less than 3 independent spatial dimensions*/

  if (pd->gv[MESH_DISPLACEMENT1] == 1) {
    v = MESH_DISPLACEMENT1;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];

#ifdef DO_NO_UNROLL
    siz = sizeof(double) * DIM * MDE;
    memset(fv->d_div_d_dmesh, 0, siz);
    siz = sizeof(double) * DIM * DIM * DIM * MDE;
    memset(fv->d_grad_d_dmesh, 0, siz);
    memset(fv->d_grad_d_dot_dmesh, 0, siz);
    for (r = 0; r < dim; r++) {
      for (i = 0; i < vdofs; i++) {
        d_ri = *esp->d[r][i];
        for (p = 0; p < VIM; p++) {
          for (q = 0; q < VIM; q++) {
            for (b = 0; b < dim; b++) {
              for (j = 0; j < mdofs; j++) {
                /*  fv->d_grad_d_dmesh[p][q] [b][j] = 0.; */
                fv->d_grad_d_dmesh[p][q][b][j] += d_ri * bfv->d_grad_phi_e_dmesh[i][r][p][q][b][j];

                fv->d_grad_d_dot_dmesh[p][q][b][j] +=
                    *esp_dot->d[r][i] * bfv->d_grad_phi_e_dmesh[i][r][p][q][b][j];
              }
            }
          }
        }
      }
    }
#else

    d_ri = *esp->d[0][0];
    d_dot_ri = *esp_dot->d[0][0];

    for (b = 0; b < dim; b++) {
      for (j = 0; j < mdofs; j++) {
        /* The first element of the (r,i) matrix is split out so that it can be initialized by
         * direct assignment This saves having to keep doing memsets on d_grad_d_mesh */

        /*	  fv->d_grad_d_dmesh[p][q] [b][j] += d_ri * bfv->d_grad_phi_e_dmesh[i][r] [p][q]
         * [b][j]; */

        fv->d_grad_d_dmesh[0][0][b][j] = d_ri * bfv->d_grad_phi_e_dmesh[0][0][0][0][b][j];
        fv->d_grad_d_dmesh[1][1][b][j] = d_ri * bfv->d_grad_phi_e_dmesh[0][0][1][1][b][j];
        fv->d_grad_d_dmesh[1][0][b][j] = d_ri * bfv->d_grad_phi_e_dmesh[0][0][1][0][b][j];
        fv->d_grad_d_dmesh[0][1][b][j] = d_ri * bfv->d_grad_phi_e_dmesh[0][0][0][1][b][j];

        /*  fv->d_grad_d_dot_dmesh[p][q] [b][j] += d_dot_ri * bfv->d_grad_phi_e_dmesh[i][r] [p][q]
         * [b][j];*/

        fv->d_grad_d_dot_dmesh[0][0][b][j] = d_dot_ri * bfv->d_grad_phi_e_dmesh[0][0][0][0][b][j];
        fv->d_grad_d_dot_dmesh[1][1][b][j] = d_dot_ri * bfv->d_grad_phi_e_dmesh[0][0][1][1][b][j];
        fv->d_grad_d_dot_dmesh[1][0][b][j] = d_dot_ri * bfv->d_grad_phi_e_dmesh[0][0][1][0][b][j];
        fv->d_grad_d_dot_dmesh[0][1][b][j] = d_dot_ri * bfv->d_grad_phi_e_dmesh[0][0][0][1][b][j];

        if (VIMis3) {
          fv->d_grad_d_dmesh[2][2][b][j] = d_ri * bfv->d_grad_phi_e_dmesh[0][0][2][2][b][j];
          fv->d_grad_d_dmesh[2][0][b][j] = d_ri * bfv->d_grad_phi_e_dmesh[0][0][2][0][b][j];
          fv->d_grad_d_dmesh[2][1][b][j] = d_ri * bfv->d_grad_phi_e_dmesh[0][0][2][1][b][j];
          fv->d_grad_d_dmesh[0][2][b][j] = d_ri * bfv->d_grad_phi_e_dmesh[0][0][0][2][b][j];
          fv->d_grad_d_dmesh[1][2][b][j] = d_ri * bfv->d_grad_phi_e_dmesh[0][0][1][2][b][j];

          fv->d_grad_d_dot_dmesh[2][2][b][j] = d_dot_ri * bfv->d_grad_phi_e_dmesh[0][0][2][2][b][j];
          fv->d_grad_d_dot_dmesh[2][0][b][j] = d_dot_ri * bfv->d_grad_phi_e_dmesh[0][0][2][0][b][j];
          fv->d_grad_d_dot_dmesh[2][1][b][j] = d_dot_ri * bfv->d_grad_phi_e_dmesh[0][0][2][1][b][j];
          fv->d_grad_d_dot_dmesh[0][2][b][j] = d_dot_ri * bfv->d_grad_phi_e_dmesh[0][0][0][2][b][j];
          fv->d_grad_d_dot_dmesh[1][2][b][j] = d_dot_ri * bfv->d_grad_phi_e_dmesh[0][0][1][2][b][j];
        }
      }
    }

    for (r = 0; r < dim; r++) {
      for (i = 0; i < vdofs; i++) {
        d_ri = *esp->d[r][i];
        d_dot_ri = *esp_dot->d[r][i];

        for (b = 0; (r + i) && (b < dim);
             b++) /* the (r+i) test is needed to make sure only the (r=0,i=0) element is excluded */
        {
          for (j = 0; j < mdofs; j++) {

            /*	  fv->d_grad_d_dmesh[p][q] [b][j] += d_ri * bfv->d_grad_phi_e_dmesh[i][r] [p][q]
             * [b][j]; */

            fv->d_grad_d_dmesh[0][0][b][j] += d_ri * bfv->d_grad_phi_e_dmesh[i][r][0][0][b][j];
            fv->d_grad_d_dmesh[1][1][b][j] += d_ri * bfv->d_grad_phi_e_dmesh[i][r][1][1][b][j];
            fv->d_grad_d_dmesh[1][0][b][j] += d_ri * bfv->d_grad_phi_e_dmesh[i][r][1][0][b][j];
            fv->d_grad_d_dmesh[0][1][b][j] += d_ri * bfv->d_grad_phi_e_dmesh[i][r][0][1][b][j];

            /*  fv->d_grad_d_dot_dmesh[p][q] [b][j] += d_dot_ri * bfv->d_grad_phi_e_dmesh[i][r]
             * [p][q] [b][j];*/

            fv->d_grad_d_dot_dmesh[0][0][b][j] +=
                d_dot_ri * bfv->d_grad_phi_e_dmesh[i][r][0][0][b][j];
            fv->d_grad_d_dot_dmesh[1][1][b][j] +=
                d_dot_ri * bfv->d_grad_phi_e_dmesh[i][r][1][1][b][j];
            fv->d_grad_d_dot_dmesh[1][0][b][j] +=
                d_dot_ri * bfv->d_grad_phi_e_dmesh[i][r][1][0][b][j];
            fv->d_grad_d_dot_dmesh[0][1][b][j] +=
                d_dot_ri * bfv->d_grad_phi_e_dmesh[i][r][0][1][b][j];

            if (VIMis3) {
              fv->d_grad_d_dmesh[2][2][b][j] += d_ri * bfv->d_grad_phi_e_dmesh[i][r][2][2][b][j];
              fv->d_grad_d_dmesh[2][0][b][j] += d_ri * bfv->d_grad_phi_e_dmesh[i][r][2][0][b][j];
              fv->d_grad_d_dmesh[2][1][b][j] += d_ri * bfv->d_grad_phi_e_dmesh[i][r][2][1][b][j];
              fv->d_grad_d_dmesh[0][2][b][j] += d_ri * bfv->d_grad_phi_e_dmesh[i][r][0][2][b][j];
              fv->d_grad_d_dmesh[1][2][b][j] += d_ri * bfv->d_grad_phi_e_dmesh[i][r][1][2][b][j];

              fv->d_grad_d_dot_dmesh[2][2][b][j] +=
                  d_dot_ri * bfv->d_grad_phi_e_dmesh[i][r][2][2][b][j];
              fv->d_grad_d_dot_dmesh[2][0][b][j] +=
                  d_dot_ri * bfv->d_grad_phi_e_dmesh[i][r][2][0][b][j];
              fv->d_grad_d_dot_dmesh[2][1][b][j] +=
                  d_dot_ri * bfv->d_grad_phi_e_dmesh[i][r][2][1][b][j];
              fv->d_grad_d_dot_dmesh[0][2][b][j] +=
                  d_dot_ri * bfv->d_grad_phi_e_dmesh[i][r][0][2][b][j];
              fv->d_grad_d_dot_dmesh[1][2][b][j] +=
                  d_dot_ri * bfv->d_grad_phi_e_dmesh[i][r][1][2][b][j];
            }
          }
        }
      }
    }
#endif

#ifdef DO_NO_UNROLL
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            /*
             * Special added Lagrangian piece for d grad_d / d mesh...
             */
            fv->d_grad_d_dmesh[p][q][b][j] += bfv->grad_phi_e[j][b][p][q];

            /*N.B. PRS: need to actually furbish this for Newmark Beta schemes */
            if (transient_run && (discontinuous_porous_media || ve_solid)) {
              fv->d_grad_d_dot_dmesh[p][q][b][j] +=
                  (1. + 2. * tran->theta) * bfv->grad_phi_e[j][b][p][q] / tran->delta_t;
            }
          }
        }
      }
    }
#else

    for (b = 0; b < dim; b++) {
      for (j = 0; j < mdofs; j++) {
        /*
         * Special added Lagrangian piece for d grad_d / d mesh...
         */
        /*  fv->d_grad_d_dmesh[p][q] [b][j] += bfv->grad_phi_e[j][b] [p][q];*/

        fv->d_grad_d_dmesh[0][0][b][j] += bfv->grad_phi_e[j][b][0][0];
        fv->d_grad_d_dmesh[1][1][b][j] += bfv->grad_phi_e[j][b][1][1];
        fv->d_grad_d_dmesh[1][0][b][j] += bfv->grad_phi_e[j][b][1][0];
        fv->d_grad_d_dmesh[0][1][b][j] += bfv->grad_phi_e[j][b][0][1];

        if (VIMis3) {
          fv->d_grad_d_dmesh[2][2][b][j] += bfv->grad_phi_e[j][b][2][2];
          fv->d_grad_d_dmesh[2][0][b][j] += bfv->grad_phi_e[j][b][2][0];
          fv->d_grad_d_dmesh[2][1][b][j] += bfv->grad_phi_e[j][b][2][1];
          fv->d_grad_d_dmesh[0][2][b][j] += bfv->grad_phi_e[j][b][0][2];
          fv->d_grad_d_dmesh[1][2][b][j] += bfv->grad_phi_e[j][b][1][2];
        }

        /*N.B. PRS: need to actually furbish this for Newmark Beta schemes */
        if (transient_run && (discontinuous_porous_media || ve_solid)) {
          /* fv->d_grad_d_dot_dmesh[p][q] [b][j] += (1.+2.*tran->theta)*bfv->grad_phi_e[j][b]
           * [p][q]/tran->delta_t;*/

          fv->d_grad_d_dot_dmesh[0][0][b][j] +=
              (1. + 2. * tran->theta) * bfv->grad_phi_e[j][b][0][0] / tran->delta_t;
          fv->d_grad_d_dot_dmesh[1][1][b][j] +=
              (1. + 2. * tran->theta) * bfv->grad_phi_e[j][b][1][1] / tran->delta_t;
          fv->d_grad_d_dot_dmesh[1][0][b][j] +=
              (1. + 2. * tran->theta) * bfv->grad_phi_e[j][b][1][0] / tran->delta_t;
          fv->d_grad_d_dot_dmesh[0][1][b][j] +=
              (1. + 2. * tran->theta) * bfv->grad_phi_e[j][b][0][1] / tran->delta_t;

          if (VIMis3) {

            fv->d_grad_d_dot_dmesh[2][2][b][j] +=
                (1. + 2. * tran->theta) * bfv->grad_phi_e[j][b][2][2] / tran->delta_t;
            fv->d_grad_d_dot_dmesh[2][0][b][j] +=
                (1. + 2. * tran->theta) * bfv->grad_phi_e[j][b][2][0] / tran->delta_t;
            fv->d_grad_d_dot_dmesh[2][1][b][j] +=
                (1. + 2. * tran->theta) * bfv->grad_phi_e[j][b][2][1] / tran->delta_t;
            fv->d_grad_d_dot_dmesh[0][2][b][j] +=
                (1. + 2. * tran->theta) * bfv->grad_phi_e[j][b][0][2] / tran->delta_t;
            fv->d_grad_d_dot_dmesh[1][2][b][j] +=
                (1. + 2. * tran->theta) * bfv->grad_phi_e[j][b][1][2] / tran->delta_t;
          }
        }
      }
    }

#endif

    /*
     * d( div(d) )/dmesh [ just trace(grad(d)) ]
     */
#ifdef DO_NO_UNROLL
    for (b = 0; b < dim; b++) {
      for (j = 0; j < mdofs; j++) {
        for (p = 0; p < VIM; p++) {
          fv->d_div_d_dmesh[b][j] += fv->d_grad_d_dmesh[p][p][b][j];
        }
      }
    }
#else

    fv->d_div_d_dmesh[0][0] = fv->d_grad_d_dmesh[0][0][0][0];
    fv->d_div_d_dmesh[1][0] = fv->d_grad_d_dmesh[1][1][1][0];
    fv->d_div_d_dmesh[0][0] = fv->d_grad_d_dmesh[1][1][0][0];
    fv->d_div_d_dmesh[1][0] = fv->d_grad_d_dmesh[0][0][1][0];

    if (VIMis3) {
      fv->d_div_d_dmesh[2][0] = fv->d_grad_d_dmesh[2][2][2][0];
      fv->d_div_d_dmesh[2][0] = fv->d_grad_d_dmesh[0][0][2][0];
      fv->d_div_d_dmesh[2][0] = fv->d_grad_d_dmesh[1][1][2][0];
      fv->d_div_d_dmesh[0][0] = fv->d_grad_d_dmesh[2][2][0][0];
      fv->d_div_d_dmesh[1][0] = fv->d_grad_d_dmesh[2][2][1][0];
    }

    for (j = 1; j < mdofs; j++) {

      /* fv->d_div_d_dmesh[b][j] += fv->d_grad_d_dmesh[p][p] [b][j]; */

      fv->d_div_d_dmesh[0][j] += fv->d_grad_d_dmesh[0][0][0][j];
      fv->d_div_d_dmesh[1][j] += fv->d_grad_d_dmesh[1][1][1][j];
      fv->d_div_d_dmesh[0][j] += fv->d_grad_d_dmesh[1][1][0][j];
      fv->d_div_d_dmesh[1][j] += fv->d_grad_d_dmesh[0][0][1][j];

      if (VIMis3) {
        fv->d_div_d_dmesh[2][j] += fv->d_grad_d_dmesh[2][2][2][j];
        fv->d_div_d_dmesh[2][j] += fv->d_grad_d_dmesh[0][0][2][j];
        fv->d_div_d_dmesh[2][j] += fv->d_grad_d_dmesh[1][1][2][j];
        fv->d_div_d_dmesh[0][j] += fv->d_grad_d_dmesh[2][2][0][j];
        fv->d_div_d_dmesh[1][j] += fv->d_grad_d_dmesh[2][2][1][j];
      }
    }

#endif

  } else if (upd->vp[pg->imtrx][MESH_DISPLACEMENT1] != -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * DIM * MDE;
    memset(fv->d_grad_d_dmesh, 0, siz);
    memset(fv->d_grad_d_dot_dmesh, 0, siz);
    siz = sizeof(double) * DIM * MDE;
    memset(fv->d_div_d_dmesh, 0, siz);
  }

  /*
   * d( grad(d_rs) )/dmesh and d( grad(d_rs) ) /dd_rs
   */

  /* Yes, I know it should be each "+p", */
  /* but this way we can collect three */
  /* components of vectors and tensors for */
  /* less than 3 independent spatial dimensions*/

  if (pd->gv[SOLID_DISPLACEMENT1]) {
    v = SOLID_DISPLACEMENT1;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];
    siz = sizeof(double) * DIM * DIM * DIM * MDE;
    memset(fv->d_grad_d_rs_dmesh, 0, siz);
    siz = sizeof(double) * DIM * MDE;
    memset(fv->d_div_d_rs_dmesh, 0, siz);

    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            /*  fv->d_grad_d_rs_dmesh[p][q] [b][j] = 0.; */
            for (r = 0; r < dim; r++) {
              for (i = 0; i < vdofs; i++) {
                fv->d_grad_d_rs_dmesh[p][q][b][j] +=
                    *esp->d_rs[r][i] * bfv->d_grad_phi_e_dmesh[i][r][p][q][b][j];
              }
            }
          }
        }
      }
    }

    /*
     * d( div(d_rs) )/dmesh [ just trace(grad(d_rs)) ]
     */
    for (b = 0; b < VIM; b++) {
      for (j = 0; j < mdofs; j++) {
        for (p = 0; p < VIM; p++) {
          fv->d_div_d_rs_dmesh[b][j] += fv->d_grad_d_rs_dmesh[p][p][b][j];
        }
      }
    }
  } else if (upd->vp[pg->imtrx][SOLID_DISPLACEMENT1] != -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * DIM * MDE;
    memset(fv->d_grad_d_rs_dmesh, 0, siz);
    siz = sizeof(double) * DIM * MDE;
    memset(fv->d_div_d_rs_dmesh, 0, siz);
  }

  /*
   * d(grad(S))/dmesh
   */

  if (pd->gv[POLYMER_STRESS11]) {
    v = POLYMER_STRESS11;
    vdofs = ei[upd->matrix_index[v]]->dof[v];
    bfv = bf[v];
    siz = sizeof(double) * MAX_MODES * DIM * DIM * DIM * DIM * MDE;
    memset(fv->d_grad_S_dmesh, 0, siz);
    siz = sizeof(double) * MAX_MODES * DIM * DIM * MDE;
    memset(fv->d_div_S_dmesh, 0, siz);

    for (mode = 0; mode < vn->modes; mode++) {
      for (p = 0; p < VIM; p++) {
        for (q = 0; q < VIM; q++) {
          for (r = 0; r < VIM; r++) {
            for (b = 0; b < VIM; b++) {
              for (j = 0; j < mdofs; j++) {
                for (i = 0; i < vdofs; i++) {
                  if (p <= q) {
                    fv->d_grad_S_dmesh[mode][r][p][q][b][j] +=
                        *esp->S[mode][p][q][i] * bfv->d_grad_phi_dmesh[i][r][b][j];
                  } else {
                    fv->d_grad_S_dmesh[mode][r][p][q][b][j] +=
                        *esp->S[mode][q][p][i] * bfv->d_grad_phi_dmesh[i][r][b][j];
                  }
                }
              }
            }
          }
        }
      }

      for (r = 0; r < VIM; r++) {
        for (q = 0; q < VIM; q++) {
          for (b = 0; b < VIM; b++) {
            for (j = 0; j < mdofs; j++) {
              fv->d_div_S_dmesh[mode][r][b][j] += fv->d_grad_S_dmesh[mode][q][q][r][b][j];
            }
          }
        }
      }

      if (pd->CoordinateSystem != CARTESIAN) {
        for (s = 0; s < VIM; s++) {
          for (r = 0; r < VIM; r++) {
            for (p = 0; p < VIM; p++) {
              for (q = 0; q < VIM; q++) {
                for (b = 0; b < VIM; b++) {
                  for (j = 0; j < mdofs; j++) {
                    fv->d_div_S_dmesh[mode][s][b][j] +=
                        fv->S[mode][p][q] *
                        (fv->d_grad_e_dq[p][r][q][b] * bfm->phi[j] * (double)delta(s, q) +
                         fv->d_grad_e_dq[q][p][s][b] * bfm->phi[j] * (double)delta(r, p));
                  }
                }
              }
            }
          }
        }
      }
    } /* end of modal loop */
  } else if (upd->vp[pg->imtrx][POLYMER_STRESS11] != -1 && okToZero) {
    siz = sizeof(double) * MAX_MODES * DIM * DIM * DIM * DIM * MDE;
    memset(fv->d_grad_S_dmesh, 0, siz);
    siz = sizeof(double) * MAX_MODES * DIM * DIM * MDE;
    memset(fv->d_div_S_dmesh, 0, siz);
  }

  /*
   * d(grad(G))/dmesh
   */

  if (pd->gv[VELOCITY_GRADIENT11]) {
    v = VELOCITY_GRADIENT11;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];

    siz = sizeof(double) * DIM * DIM * DIM * DIM * MDE;
    memset(fv->d_grad_G_dmesh, 0, siz);
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(fv->d_div_G_dmesh, 0, siz);
    memset(fv->d_div_Gt_dmesh, 0, siz);

    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        for (r = 0; r < VIM; r++) {
          for (b = 0; b < VIM; b++) {
            for (j = 0; j < mdofs; j++) {
              for (i = 0; i < vdofs; i++) {
                fv->d_grad_G_dmesh[r][p][q][b][j] +=
                    *esp->G[p][q][i] * bfv->d_grad_phi_dmesh[i][r][b][j];
              }
            }
          }
        }
      }
    }

    for (r = 0; r < VIM; r++) {
      for (q = 0; q < VIM; q++) {
        for (b = 0; b < VIM; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_div_G_dmesh[r][b][j] += fv->d_grad_G_dmesh[q][q][r][b][j];
          }
        }
      }
    }

    if (pd->CoordinateSystem != CARTESIAN) {
      for (s = 0; s < VIM; s++) {
        for (r = 0; r < VIM; r++) {
          for (p = 0; p < VIM; p++) {
            for (q = 0; q < VIM; q++) {
              for (b = 0; b < dim; b++) {
                for (j = 0; j < mdofs; j++) {
                  fv->d_div_G_dmesh[s][b][j] +=
                      fv->G[p][q] *
                      (fv->d_grad_e_dq[p][r][q][b] * bfm->phi[j] * (double)delta(s, q) +
                       fv->d_grad_e_dq[q][p][s][b] * bfm->phi[j] * (double)delta(r, p));
                }
              }
            }
          }
        }
      }
    }

    for (r = 0; r < VIM; r++) {
      for (q = 0; q < VIM; q++) {
        for (b = 0; b < VIM; b++) {
          for (j = 0; j < mdofs; j++) {
            fv->d_div_Gt_dmesh[r][b][j] += fv->d_grad_G_dmesh[q][r][q][b][j];
          }
        }
      }
    }

    if (pd->CoordinateSystem != CARTESIAN) {
      for (s = 0; s < VIM; s++) {
        for (r = 0; r < VIM; r++) {
          for (p = 0; p < VIM; p++) {
            for (q = 0; q < VIM; q++) {
              for (b = 0; b < VIM; b++) {
                for (j = 0; j < mdofs; j++) {
                  fv->d_div_G_dmesh[s][b][j] +=
                      fv->G[q][p] *
                      (fv->d_grad_e_dq[p][r][q][b] * bfm->phi[j] * (double)delta(s, q) +
                       fv->d_grad_e_dq[q][p][s][b] * bfm->phi[j] * (double)delta(r, p));
                }
              }
            }
          }
        }
      }
    }

  } else if (upd->vp[pg->imtrx][VELOCITY_GRADIENT11] != -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * DIM * DIM * MDE;
    memset(fv->d_grad_G_dmesh, 0, siz);
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(fv->d_div_G_dmesh, 0, siz);
    memset(fv->d_div_Gt_dmesh, 0, siz);
  }

  /*
   * d(grad(vd))/dmesh
   */
  if (pd->gv[VORT_DIR1]) {
    v = VORT_DIR1;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];

    siz = sizeof(double) * DIM * DIM * DIM * MDE;
    memset(fv->d_grad_vd_dmesh, 0, siz);

    for (r = 0; r < DIM; r++) {
      for (i = 0; i < vdofs; i++) {
        v_ri = *esp->vd[r][i];
        for (p = 0; p < DIM; p++) {
          for (q = 0; q < DIM; q++) {
            for (b = 0; b < VIM; b++) {
              for (j = 0; j < mdofs; j++) {
                fv->d_grad_vd_dmesh[p][q][b][j] += v_ri * bfv->d_grad_phi_e_dmesh[i][r][p][q][b][j];
              }
            }
          }
        }
      }
    }

    /*
     * d( div(vd) )/dmesh
     */

    for (b = 0; b < VIM; b++) {
      for (j = 0; j < mdofs; j++) {
        for (p = 0; p < DIM; p++) {
          fv->d_div_vd_dmesh[b][j] += fv->d_grad_vd_dmesh[p][p][b][j];
        }
      }
    }
  } else if (upd->vp[pg->imtrx][VORT_DIR1] != -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * DIM * MDE;
    memset(fv->d_grad_vd_dmesh, 0, siz);
  }

  /*
   *   d(n_dot_curl_s_v)/dmesh
   *        This is carried out on a shell
   */
  if (pd->gv[N_DOT_CURL_V]) {
    v = N_DOT_CURL_V;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(fv->d_grad_n_dot_curl_s_v_dmesh, 0, siz);
    for (p = 0; p < VIM; p++) {
      for (b = 0; b < dim; b++) {
        for (j = 0; j < mdofs; j++) {
          for (i = 0; i < vdofs; i++) {
            fv->d_grad_n_dot_curl_s_v_dmesh[p][b][j] +=
                (*esp->n_dot_curl_s_v[i] * bfv->d_grad_phi_dmesh[i][p][b][j]);
          }
        }
      }
    }
  } else if (upd->vp[pg->imtrx][N_DOT_CURL_V] != -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(fv->d_grad_n_dot_curl_s_v_dmesh, 0, siz);
  }

  /*
   *   d(grad_div_s_v_dmesh[b][jvar][jShell]);
   *        This is carried out on a shell
   */
  if (pd->gv[SHELL_SURF_DIV_V]) {
    v = SHELL_SURF_DIV_V;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(fv->d_grad_div_s_v_dmesh, 0, siz);
    for (p = 0; p < VIM; p++) {
      for (b = 0; b < dim; b++) {
        for (j = 0; j < mdofs; j++) {
          for (i = 0; i < vdofs; i++) {
            fv->d_grad_div_s_v_dmesh[p][b][j] +=
                (*esp->div_s_v[i] * bfv->d_grad_phi_dmesh[i][p][b][j]);
          }
        }
      }
    }
  } else if (upd->vp[pg->imtrx][SHELL_SURF_DIV_V] != -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(fv->d_grad_div_s_v_dmesh, 0, siz);
  }

  /*
   *   d_grad_curv_dmesh[b][jvar][jShell])
   *        This is carried out on a shell
   */
  if (pd->gv[SHELL_SURF_CURV]) {
    v = SHELL_SURF_CURV;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(fv->d_grad_curv_dmesh, 0, siz);
    for (p = 0; p < VIM; p++) {
      for (b = 0; b < dim; b++) {
        for (j = 0; j < mdofs; j++) {
          for (i = 0; i < vdofs; i++) {
            fv->d_grad_curv_dmesh[p][b][j] += (*esp->curv[i] * bfv->d_grad_phi_dmesh[i][p][b][j]);
          }
        }
      }
    }
  } else if (upd->vp[pg->imtrx][SHELL_SURF_CURV] != -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(fv->d_grad_curv_dmesh, 0, siz);
  }

  /*
   *   d_serialgrad_div_s_v_dmesh[b][jvar][jShell]);
   *        This is carried out on a shell
   */
  if (pd->gv[GRAD_S_V_DOT_N1]) {
    for (r = 0; r < dim; r++) {
      v = GRAD_S_V_DOT_N1 + r;
      bfv = bf[v];
      vdofs = ei[upd->matrix_index[v]]->dof[v];
      siz = sizeof(double) * DIM * DIM * DIM * MDE;
      memset(fv->d_serialgrad_grad_s_v_dot_n_dmesh, 0, siz);
      for (p = 0; p < VIM; p++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < mdofs; j++) {
            for (i = 0; i < vdofs; i++) {
              fv->d_serialgrad_grad_s_v_dot_n_dmesh[r][p][b][j] +=
                  (*esp->grad_v_dot_n[r][i] * bfv->d_grad_phi_dmesh[i][p][b][j]);
            }
          }
        }
      }
    }
  } else if (upd->vp[pg->imtrx][GRAD_S_V_DOT_N1] != -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * DIM * MDE;
    memset(fv->d_serialgrad_grad_s_v_dot_n_dmesh, 0, siz);
  }

  /*
   * d(grad(sh_J))/dmesh
   */
  if (pd->gv[SHELL_DIFF_FLUX]) {
    v = SHELL_DIFF_FLUX;
    bfv = bf[v];
    vdofs = ei[upd->matrix_index[v]]->dof[v];

    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_sh_J_dmesh[0][0][0]), 0, siz);

    for (p = 0; p < VIM; p++) {
      for (b = 0; b < dim; b++) {
        for (j = 0; j < mdofs; j++) {
          for (i = 0; i < vdofs; i++) {
            fv->d_grad_sh_J_dmesh[p][b][j] += *esp->sh_J[i] * bfv->d_grad_phi_dmesh[i][p][b][j];
          }
        }
      }
    }
  } else if (upd->vp[pg->imtrx][SHELL_DIFF_FLUX] == -1 && okToZero) {
    siz = sizeof(double) * DIM * DIM * MDE;
    memset(&(fv->d_grad_sh_J_dmesh[0][0][0]), 0, siz);
  }

  /*
   * Now load up the sensitivity of the volume element to the
   * mesh variables...
   *
   * d ( h3 ) / d ( mesh_bj )
   */
  if (pd->gv[MESH_DISPLACEMENT1] == 1) {
    for (b = 0; b < VIM; b++) {
      for (j = 0; j < mdofs; j++) {
        fv->dh3dmesh[b][j] = fv->dh3dq[b] * bfm->phi[j];
      }
    }
  }

  return (status);
}