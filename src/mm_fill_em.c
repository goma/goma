/************************************************************************ *
 * Goma - Multiphysics finite element software                             *
 * Sandia National Laboratories                                            *
 *                                                                         *
 * Copyright (c) 2019 GOMA                                                 *
 *                                                                         *
 * Authors: Robert Secor and Andrew Cochrane                               *
 *                                                                         *
 * This software is distributed under the GNU General Public License.      *
\************************************************************************/

/* Routines included in this file
 *
 * assemble_emwave -- assemble terms for EM frequency domain wave equations
 *
 */

#include "mm_qtensor_model.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

/* GOMA include files */
#define GOMA_MM_FILL_EM_C
#include <complex.h>
#ifdef I
#undef I
#endif

#include "el_elm.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "mm_fill_em.h"
#include "mm_fill_ls.h"
#include "mm_fill_terms.h"
#include "mm_mp.h"
#include "mm_mp_structs.h"
#include "rf_bc_const.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "std.h"

/*  _______________________________________________________________________  */

/* assemble_emwave -- assemble terms (Residual &| Jacobian) for EM harmonic
 *                    wave equations
 *
 * in:
 *     ei -- pointer to Element Indecesstructure
 *     pd -- pointer to Problem Descriptionstructure
 *     af -- pointer to Action Flagstructure
 *     bf -- pointer to Basis Functionstructure
 *     fv -- pointer to Field Variablestructure
 *       fv_old -- pointer to old Diet Field Variablestructure
 *       fv_dot -- pointer to dot Diet Field Variablestructure
 *     cr -- pointer to Constitutive Relationstructure
 *     md -- pointer to Mesh Derivativestructure
 *     me -- pointer to Material Entitystructure
 *
 * out:
 *     a   -- gets loaded up with proper contribution
 *     lec -- gets loaded up with local contributions to resid, Jacobian
 *     r   -- residual RHS vector
 *
 * Created: Thu October 26, 2018 - Robert Secor
 *
 * Note: currently we do a "double load" into the addresses in the global
 *       "a" matrix and resid vector held in "esp", as well as into the
 *       local accumulators in "lec".
 *
 */

int assemble_emwave(double time, /* present time value */
                    double tt,   /* parameter to vary time integration from
                                  * explicit (tt = 1) to implicit (tt = 0) */
                    double dt,   /* current time step size */
                    const PG_DATA *pg_data,
                    const int em_eqn, /* emwave eqn id and var id	*/
                    const int em_var,
                    const int em_conjvar) {
  int eqn, var, peqn, pvar, dim, p, q, b, w, i, j, status;
  int dir = 0; /* identity of conjugate variable  */

  dbl EMF = 0, EMF_conj = 0; /* acoustic pressure	*/
  dbl omega, emf_coeff = 0, conj_coeff = 0;
  dbl emf_coeff_dn = 0, conj_coeff_dn = 0;
  dbl emf_coeff_dk = 0, conj_coeff_dk = 0;
  // dbl mag_permeability=12.57e-07;  // H/m
  dbl mag_permeability = 12.57e-07; // H/m
  int cross_field_var;
  dbl cross_field[DIM];

  struct emwave_stabilization em_stab;
  em_stab.em_eqn = em_eqn;
  em_stab.em_var = em_var;
  em_stab.type = EM_STAB_DPHI_DIV; // enum supports phi_div, dphi_div,
                                   // divphi_div, phi_divsquared and
                                   // dphi_divsquared

  dbl n; /* Refractive index. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_n_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_n = &d_n_struct;

  dbl k; /* Acoustic wave number. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  dbl advection; /* For terms and their derivatives */

  dbl diffusion;
  dbl diff_a, diff_b, diff_c, diff_d;

  /*
   * Galerkin weighting functions for i-th energy residuals
   * and some of their derivatives...
   */

  dbl phi_i;
  dbl grad_phi_i[DIM];

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;

  dbl h3;          /* Volume element (scale factors). */
  dbl dh3dmesh_bj; /* Sensitivity to (b,j) mesh dof. */

  dbl det_J;

  dbl d_det_J_dmeshbj;        /* for specified (b,j) mesh dof */
  dbl dgrad_phi_i_dmesh[DIM]; /* ditto.  */
  dbl wt;

  /* initialize grad_phi_i */
  for (i = 0; i < DIM; i++) {
    grad_phi_i[i] = 0;
  }

  /*   static char yo[] = "assemble_acoustic";*/
  status = 0;
  /*
   * Unpack variables from structures for local convenience...
   */
  dim = pd->Num_Dim;
  eqn = em_eqn;
  /*
   * Bail out fast if there's nothing to do...
   */
  if (!pd->e[pg->imtrx][eqn]) {
    return (status);
  }

  wt = fv->wt;           /* Gauss point weight. */
  h3 = fv->h3;           /* Differential volume element. */
  det_J = bf[eqn]->detJ; /* Really, ought to be mesh eqn. */

  /*
   * Material property constants, etc. Any variations at this Gauss point
   * are calculated with user-defined subroutine options.
   *
   * Properties to consider here are R and k. R and k we will allow to vary
   * with temperature, spatial coordinates, and species concentration.
   */

  omega = upd->EM_Frequency;
  n = refractive_index(d_n, time);

  k = extinction_index(d_k, time);

  // Compute complex impedance
  complex cpx_refractive_index, cpx_rel_permittivity,
      cpx_permittivity; //, impedance;

  cpx_refractive_index = n + _Complex_I * k; // k > 0 is extinction
  cpx_rel_permittivity = SQUARE(cpx_refractive_index);
  cpx_permittivity = cpx_rel_permittivity * mp->permittivity;

  // impedance = csqrt(mag_permeability/rel_permittivity);

  switch (em_var) {
  case EM_E1_REAL:
    EMF = fv->em_er[0];
    EMF_conj = fv->em_ei[0];
    dir = 0;
    break;
  case EM_E1_IMAG:
    EMF = fv->em_ei[0];
    EMF_conj = fv->em_er[0];
    dir = 0;
    break;
  case EM_E2_REAL:
    EMF = fv->em_er[1];
    EMF_conj = fv->em_ei[1];
    dir = 1;
    break;
  case EM_E2_IMAG:
    EMF = fv->em_ei[1];
    EMF_conj = fv->em_er[1];
    dir = 1;
    break;
  case EM_E3_REAL:
    EMF = fv->em_er[2];
    EMF_conj = fv->em_ei[2];
    dir = 2;
    break;
  case EM_E3_IMAG:
    EMF = fv->em_ei[2];
    EMF_conj = fv->em_er[2];
    dir = 2;
    break;
  case EM_H1_REAL:
    EMF = fv->em_hr[0];
    EMF_conj = fv->em_hi[0];
    dir = 0;
    break;
  case EM_H1_IMAG:
    EMF = fv->em_hi[0];
    EMF_conj = fv->em_hr[0];
    dir = 0;
    break;
  case EM_H2_REAL:
    EMF = fv->em_hr[1];
    EMF_conj = fv->em_hi[1];
    dir = 1;
    break;
  case EM_H2_IMAG:
    EMF = fv->em_hi[1];
    EMF_conj = fv->em_hr[1];
    dir = 1;
    break;
  case EM_H3_REAL:
    EMF = fv->em_hr[2];
    EMF_conj = fv->em_hi[2];
    dir = 2;
    break;
  case EM_H3_IMAG:
    EMF = fv->em_hi[2];
    EMF_conj = fv->em_hr[2];
    dir = 2;
    break;
  default:
    GOMA_EH(GOMA_ERROR, "Invalid EM variable name.\n");
    return -1;
  }
  switch (em_var) {
  case EM_E1_REAL:
  case EM_E2_REAL:
  case EM_E3_REAL:
    emf_coeff = -omega * cimag(cpx_permittivity);
    conj_coeff = -omega * creal(cpx_permittivity);
    // emf_coeff = 0;
    // conj_coeff = 0;
    emf_coeff_dn = -omega * cimag(cpx_permittivity) / n;
    emf_coeff_dk = -omega * cimag(cpx_permittivity) / k;
    conj_coeff_dn = -omega * creal(cpx_permittivity) / n;
    conj_coeff_dk = -omega * creal(cpx_permittivity) / k;
    cross_field_var = EM_H1_REAL;
    for (p = 0; p < VIM; p++) {
      cross_field[p] = fv->em_hr[p];
    }
    em_stab.stabilization_field_var = EM_E1_REAL;
    calc_emwave_stabilization_term(&em_stab, 1.0);
    break;
  case EM_E1_IMAG:
  case EM_E2_IMAG:
  case EM_E3_IMAG:
    emf_coeff = -omega * cimag(cpx_permittivity);
    conj_coeff = omega * creal(cpx_permittivity);
    // emf_coeff = 0;
    // conj_coeff = 0;
    emf_coeff_dn = -omega * cimag(cpx_permittivity) / n;
    emf_coeff_dk = -omega * cimag(cpx_permittivity) / k;
    conj_coeff_dn = omega * creal(cpx_permittivity) / n;
    conj_coeff_dk = omega * creal(cpx_permittivity) / k;
    cross_field_var = EM_H1_IMAG;
    for (p = 0; p < VIM; p++) {
      cross_field[p] = fv->em_hi[p];
    }
    em_stab.stabilization_field_var = EM_E1_IMAG;
    calc_emwave_stabilization_term(&em_stab, 1.0);
    break;
  case EM_H1_REAL:
  case EM_H2_REAL:
  case EM_H3_REAL:
    emf_coeff = 0;
    conj_coeff = -omega * mag_permeability;
    // conj_coeff = 0;
    emf_coeff_dn = 0;
    emf_coeff_dk = 0;
    conj_coeff_dn = 0;
    conj_coeff_dk = 0;
    cross_field_var = EM_E1_REAL;
    for (p = 0; p < VIM; p++) {
      cross_field[p] = fv->em_er[p];
    }
    em_stab.stabilization_field_var = EM_H1_REAL;
    calc_emwave_stabilization_term(&em_stab, 1.0);
    break;
  case EM_H1_IMAG:
  case EM_H2_IMAG:
  case EM_H3_IMAG:
    emf_coeff = 0;
    conj_coeff = omega * mag_permeability;
    // conj_coeff = 0;
    emf_coeff_dn = 0;
    emf_coeff_dk = 0;
    conj_coeff_dn = 0;
    conj_coeff_dk = 0;
    cross_field_var = EM_E1_IMAG;
    for (p = 0; p < VIM; p++) {
      cross_field[p] = fv->em_ei[p];
    }
    em_stab.stabilization_field_var = EM_H1_IMAG;
    calc_emwave_stabilization_term(&em_stab, 1.0);
    break;
  default:
    GOMA_EH(GOMA_ERROR, "assemble_emwave must be called with a usable em_var\n");
    return -1;
  }
  /*
   * Residuals___________________________________________________________
   */

  if (af->Assemble_Residual) {
    eqn = em_eqn;
    peqn = upd->ep[pg->imtrx][eqn];
    var = em_var;
    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

      /* this is an optimization for xfem */
      if (xfem != NULL) {
        int xfem_active, extended_dof, base_interp, base_dof;
        xfem_dof_state(i, pd->i[pg->imtrx][eqn], ei[pg->imtrx]->ielem_shape, &xfem_active,
                       &extended_dof, &base_interp, &base_dof);
        if (extended_dof && !xfem_active)
          continue;
      }
      phi_i = bf[eqn]->phi[i];

      advection = 0.;
      if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
        advection += emf_coeff * EMF;
        advection += conj_coeff * EMF_conj;
        advection *= phi_i * h3;
        advection *= det_J * wt;
        advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
      }

      diffusion = 0.;
      if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
        for (p = 0; p < VIM; p++) {
          grad_phi_i[p] = bf[eqn]->grad_phi[i][p];
        }

        for (p = 0; p < VIM; p++) {
          for (q = 0; q < VIM; q++) {
            diffusion -= permute(p, q, dir) * grad_phi_i[p] * cross_field[q];
          }
        }

        diffusion += em_stab.residual_term[i];

        diffusion *= det_J * wt;
        diffusion *= h3;
        diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
      }

      lec->R[LEC_R_INDEX(peqn, i)] += advection + diffusion;
    }
  }

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    eqn = em_eqn;
    peqn = upd->ep[pg->imtrx][eqn];
    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

      // this is an optimization for xfem
      if (xfem != NULL) {
        int xfem_active, extended_dof, base_interp, base_dof;
        xfem_dof_state(i, pd->i[pg->imtrx][eqn], ei[pg->imtrx]->ielem_shape, &xfem_active,
                       &extended_dof, &base_interp, &base_dof);
        if (extended_dof && !xfem_active)
          continue;
      }

      phi_i = bf[eqn]->phi[i];

      /*
       * Set up some preliminaries that are needed for the (a,i)
       * equation for bunches of (b,j) column variables...
       */

      for (p = 0; p < VIM; p++) {
        grad_phi_i[p] = bf[eqn]->grad_phi[i][p];
      }

      /*
       * EMF
       */
      var = em_var;
      if (pd->gv[var]) {
        pvar = upd->vp[pg->imtrx][var];
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];

          advection = 0.;
          if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
            advection += phi_i * emf_coeff * phi_j * det_J * wt;
            advection *= h3;
            advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
          }

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += advection;
        }
      }
      /*
       *  EMF_conj
       */
      var = em_conjvar;
      if (pd->gv[var]) {
        pvar = upd->vp[pg->imtrx][var];
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];

          advection = 0.;
          if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
            advection += phi_i * conj_coeff * phi_j * det_J * wt;
            advection *= h3;
            advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
          }

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += advection;
        }
      }
      /*
       *  cross_field
       */
      for (b = 0; b < dim; b++) {
        var = cross_field_var + b;
        if (pd->gv[var]) {
          pvar = upd->vp[pg->imtrx][var];
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            // for a cross product, this isn't quite right
            // but it will work as along as all scalar
            // components of the vector field have the
            // same basis functions
            phi_j = bf[var]->phi[j];

            diffusion = 0.;
            if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
              for (p = 0; p < VIM; p++) {
                grad_phi_i[p] = bf[eqn]->grad_phi[i][p];
              }

              for (p = 0; p < VIM; p++) {
                for (q = 0; q < VIM; q++) {
                  diffusion -= permute(p, q, dir) * grad_phi_i[p] * delta(q, b) * phi_j;
                }
              }

              // diffusion += stabilization_coefficient*phi_i
              //     *bf[var]->grad_phi[j][b];

              diffusion *= det_J * wt;
              diffusion *= h3;
              diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];

              lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += diffusion;
            }
          }
        }
      }
      /*
       *  stabilization field
       */
      for (b = 0; b < dim; b++) {
        var = em_stab.stabilization_field_var + b;
        if (pd->gv[var]) {
          pvar = upd->vp[pg->imtrx][var];
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            diffusion = 0.;
            if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
              diffusion += em_stab.jacobian_term[i][b][j];

              diffusion *= det_J * wt;
              diffusion *= h3;
              diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];

              lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += diffusion;
            }
          }
        }
      }
      /*
       * J_e_T
       */
      var = TEMPERATURE;
      if (pd->gv[var]) {
        pvar = upd->vp[pg->imtrx][var];
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];

          advection = 0.;
          if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
            advection += phi_i *
                         (EMF * (emf_coeff_dn * d_n->T[j] + emf_coeff_dk * d_k->T[j]) +
                          EMF_conj * (conj_coeff_dn * d_n->T[j] + conj_coeff_dk * d_k->T[j])) *
                         det_J * wt;
            advection *= h3;
            advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
          }

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += advection;
        }
      }

      /*
       * J_e_d
       */
      for (b = 0; b < dim; b++) {
        var = MESH_DISPLACEMENT1 + b;
        if (pd->gv[var]) {
          pvar = upd->vp[pg->imtrx][var];
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            dh3dmesh_bj = fv->dh3dq[b] * bf[var]->phi[j];

            d_det_J_dmeshbj = bf[eqn]->d_det_J_dm[b][j];

            advection = 0.;
            if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
              advection +=
                  phi_i *
                  (EMF * (emf_coeff_dn * d_n->X[b][j] + emf_coeff_dk * d_k->X[b][j]) +
                   EMF_conj * (conj_coeff_dn * d_n->X[b][j] + conj_coeff_dk * d_k->X[b][j])) *
                  det_J * h3 * wt;
              advection += phi_i * (emf_coeff * EMF + conj_coeff * EMF_conj) *
                           (d_det_J_dmeshbj * h3 + det_J * dh3dmesh_bj) * wt;
              advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
            }

            /*
             * multiple parts:
             * 	diff_a = Int(...d(grad_phi_i)/dmesh.q h3 |Jv|)
             *	diff_b = Int(...grad_phi_i.d(q)/dmesh h3 |Jv|)
             *	diff_c = Int(...grad_phi_i.q h3 d(|Jv|)/dmesh)
             *	diff_d = Int(...grad_phi_i.q dh3/dmesh |Jv|  )
             */
            diffusion = 0.;
            if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
              diff_a = 0.;
              for (p = 0; p < dim; p++) {
                dgrad_phi_i_dmesh[p] = bf[eqn]->d_grad_phi_dmesh[i][p][b][j];

                diff_a += dgrad_phi_i_dmesh[p];
              }
              diff_a *= det_J * h3 * wt;

              diff_b = 0.;
              for (p = 0; p < VIM; p++) {
                diff_b += grad_phi_i[p];
              }
              diff_b *= det_J * h3 * wt;

              diff_c = 0.;
              for (p = 0; p < dim; p++) {
                diff_c += grad_phi_i[p];
              }
              diff_c *= d_det_J_dmeshbj * h3 * wt;

              diff_d = 0.;
              for (p = 0; p < dim; p++) {
                diff_d += grad_phi_i[p];
              }
              diff_d *= det_J * dh3dmesh_bj * wt;

              diffusion = diff_a + diff_b + diff_c + diff_d;

              diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
            }

            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += advection + diffusion;
          }
        }
      }

      /*
       * J_e_c
       */
      var = MASS_FRACTION;
      if (pd->e[pg->imtrx][eqn] && pd->gv[var]) {
        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            advection = 0.;
            if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
              advection +=
                  phi_i *
                  (EMF * (emf_coeff_dn * d_n->C[w][j] + emf_coeff_dk * d_k->C[w][j]) +
                   EMF_conj * (conj_coeff_dn * d_n->C[w][j] + conj_coeff_dk * d_k->C[w][j])) *
                  det_J * wt;
              advection *= h3;
              advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
            }

            lec->J[LEC_J_INDEX(peqn, MAX_PROB_VAR + w, i, j)] += advection;
          }
        }
      }
    }
  }

  return (status);
} /* end of assemble_em */

/* assemble_ewave -- assemble terms (Residual &| Jacobian) for EM harmonic
 *                    wave equations
 *                   Substitue H for curl(E) so that goma solves
 *                   1 complex vector equation and 1 complex vector primitive
 *
 * in:
 *     ei -- pointer to Element Indecesstructure
 *     pd -- pointer to Problem Descriptionstructure
 *     af -- pointer to Action Flagstructure
 *     bf -- pointer to Basis Functionstructure
 *     fv -- pointer to Field Variablestructure
 *       fv_old -- pointer to old Diet Field Variablestructure
 *       fv_dot -- pointer to dot Diet Field Variablestructure
 *     cr -- pointer to Constitutive Relationstructure
 *     md -- pointer to Mesh Derivativestructure
 *     me -- pointer to Material Entitystructure
 *
 * out:
 *     a   -- gets loaded up with proper contribution
 *     lec -- gets loaded up with local contributions to resid, Jacobian
 *     r   -- residual RHS vector
 *
 * Created: Wednesday April 22, 2020 - Andrew Cochrane
 *
 *
 */

int assemble_ewave(double time,      // present time
                   double tt,        // time integration method parameter
                   double dt,        // current time step size
                   const int em_eqn, // eqn id
                   const int em_var  //  variable id - should match me_eqn
) {
  int eqn, peqn, i;
  dbl mag_permeability = mp->magnetic_permeability;
  double omega, re_coeff, im_coeff;
  dbl n; /* Refractive index. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_n_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_n = &d_n_struct;

  dbl k; /* Acoustic wave number. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  int field_index;
  double curl_field[DIM];

  dbl advection; /* For terms and their derivatives */

  dbl diffusion;

  dbl phi_i;
  dbl grad_phi_i[DIM];

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl h3; /* Volume element (scale factors). */
  dbl det_J;
  dbl wt;

  struct emwave_stabilization em_stab;
  em_stab.em_eqn = em_eqn;
  em_stab.em_var = em_var;
  em_stab.type = EM_STAB_NONE; // enum supports phi_div, dphi_div,
                               // divphi_div, phi_divsquared and
                               // dphi_divsquared

  /* initialize grad_phi_i */
  for (i = 0; i < DIM; i++) {
    grad_phi_i[i] = 0;
  }
  // initialize curl field
  for (int p = 0; p < DIM; p++) {
    curl_field[p] = 0.0;
  }

  eqn = em_eqn;
  /*
   * Bail out fast if there's nothing to do...
   */
  if (!pd->e[pg->imtrx][eqn]) {
    return (1);
  }

  wt = fv->wt;           /* Gauss point weight. */
  h3 = fv->h3;           /* Differential volume element. */
  det_J = bf[eqn]->detJ; /* Really, ought to be mesh eqn. */

  /*
   * Material property constants, etc. Any variations at this Gauss point
   * are calculated with user-defined subroutine options.
   *
   * Properties to consider here are R and k. R and k we will allow to vary
   * with temperature, spatial coordinates, and species concentration.
   */

  omega = upd->EM_Frequency;
  n = refractive_index(d_n, time);

  k = extinction_index(d_k, time);

  // Compute complex material properties
  complex double cpx_refractive_index, cpx_rel_permittivity,
      cpx_permittivity; //, impedance;
  double r_elperm, i_elperm;

  cpx_refractive_index = n + _Complex_I * k; // k > 0 is extinction
  cpx_rel_permittivity = SQUARE(cpx_refractive_index);
  cpx_permittivity = cpx_rel_permittivity * mp->permittivity;

  // assumed to be constant in an element block
  r_elperm = creal(cpx_permittivity);
  i_elperm = cimag(cpx_permittivity);

  switch (em_eqn) {
  case R_EM_E1_REAL:
  case R_EM_E2_REAL:
  case R_EM_E3_REAL:
    field_index = em_eqn - R_EM_E1_REAL;
    re_coeff = -omega * omega * mag_permeability * r_elperm;
    im_coeff = omega * omega * mag_permeability * i_elperm;
    for (int p = 0; p < DIM; p++) {
      for (int q = 0; q < DIM; q++) {
        for (int r = 0; r < DIM; r++) {
          curl_field[p] += permute(p, q, r) * fv->grad_em_er[q][r];
        }
      }
    }
    em_stab.stabilization_field_var = EM_E1_REAL;
    calc_emwave_stabilization_term(&em_stab, 1.0);
    break;
  case R_EM_E1_IMAG:
  case R_EM_E2_IMAG:
  case R_EM_E3_IMAG:
    field_index = em_eqn - R_EM_E1_IMAG;
    re_coeff = -omega * omega * mag_permeability * i_elperm;
    im_coeff = -omega * omega * mag_permeability * r_elperm;
    for (int p = 0; p < DIM; p++) {
      for (int q = 0; q < DIM; q++) {
        for (int r = 0; r < DIM; r++) {

          curl_field[p] += permute(p, q, r) * fv->grad_em_ei[q][r];
        }
      }
    }
    em_stab.stabilization_field_var = EM_E1_IMAG;
    calc_emwave_stabilization_term(&em_stab, 1.0);
    break;
  default:
    field_index = 0;
    re_coeff = 0;
    im_coeff = 0;
    GOMA_EH(GOMA_ERROR, "assemble_ewave must be called with only electric field for em_eqn!");
    break;
  }
  if (af->Assemble_Residual) {
    eqn = em_eqn;
    peqn = upd->ep[pg->imtrx][eqn];
    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      phi_i = bf[eqn]->phi[i];

      advection = 0.;
      if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
        advection += re_coeff * fv->em_er[field_index];
        advection += im_coeff + fv->em_ei[field_index];
        advection *= phi_i * h3;
        advection *= det_J * wt;
        advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
      }
      diffusion = 0.;
      if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
        for (int p = 0; p < VIM; p++) {
          grad_phi_i[p] = bf[eqn]->grad_phi[i][p];
        }

        for (int p = 0; p < VIM; p++) {
          for (int q = 0; q < VIM; q++) {
            diffusion -=
                delta(field_index, p) * permute(field_index, p, q) * grad_phi_i[p] * curl_field[q];
          }
        }

        diffusion += em_stab.residual_term[i];

        diffusion *= det_J * wt;
        diffusion *= h3;
        diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
      }

      lec->R[LEC_R_INDEX(peqn, i)] += advection + diffusion;
    }
  }
  return (0);
} // end of assemble_ewave

/* assemble_ewave_curlcurl -- assemble terms (Residual &| Jacobian) for EM harmonic
 *                    wave equations
 *                   Substitue H for curl(E) so that goma solves
 *                   1 complex vector equation and 1 complex vector primitive
 *
 * Note: The product rule for cross products is used for the integration by parts
 *       div(A cross B) = curl(A) dot B - A dot curl(B)
 *       div[phi cross curl(E)] = curl(phi) dot curl(E) - phi dot curl(curl(E))
 *
 *
 * in:
 *     ei -- pointer to Element Indecesstructure
 *     pd -- pointer to Problem Descriptionstructure
 *     af -- pointer to Action Flagstructure
 *     bf -- pointer to Basis Functionstructure
 *     fv -- pointer to Field Variablestructure
 *       fv_old -- pointer to old Diet Field Variablestructure
 *       fv_dot -- pointer to dot Diet Field Variablestructure
 *     cr -- pointer to Constitutive Relationstructure
 *     md -- pointer to Mesh Derivativestructure
 *     me -- pointer to Material Entitystructure
 *
 * out:
 *     a   -- gets loaded up with proper contribution
 *     lec -- gets loaded up with local contributions to resid, Jacobian
 *     r   -- residual RHS vector
 *
 * Created: Wednesday April 22, 2020 - Andrew Cochrane
 * Re-written: June 16, 2020 - Weston Ortiz
 * Modified for testing with method of manufactured solutions: Jun 30, 2020 - Andrew Cochrane
 *
 */

int assemble_ewave_curlcurl(double time,      // present time
                            double tt,        // time integration method parameter
                            double dt,        // current time step size
                            const int em_eqn, // eqn id
                            const int em_var  //  variable id - should match me_eqn
) {
  int eqn;
  dbl mag_permeability = mp->magnetic_permeability;
  double omega, re_coeff, im_coeff;
  /*
  struct emwave_stabilization em_stab;
  em_stab.em_eqn = em_eqn;
  em_stab.em_var = em_var;
  em_stab.type = EM_STAB_DPHI_DIV; // enum supports phi_div, dphi_div,
                           // divphi_div, phi_divsquared and
                           // dphi_divsquared
  */
  eqn = em_eqn;
  /*
   * Bail out fast if there's nothing to do...
   * But we might have the wrong eqn
   */
  for (int a = 0; a < DIM; a++) {
    if (pd->e[pg->imtrx][R_EM_E1_REAL + a]) {
      eqn = R_EM_E1_REAL + a;
      break;
    }
  }
  if (!pd->e[pg->imtrx][eqn]) {
    return (1);
  }

  omega = upd->EM_Frequency;
  dbl n; /* Refractive index. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_n_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_n = &d_n_struct;

  dbl k; /* Extinction coefficient */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  n = refractive_index(d_n, time);
  k = extinction_index(d_k, time);

  // Compute complex material properties
  complex double cpx_refractive_index, cpx_rel_permittivity,
      cpx_permittivity; //, impedance;
  double r_elperm, i_elperm;

  cpx_refractive_index = n + _Complex_I * k; // k > 0 is extinction
  cpx_rel_permittivity = SQUARE(cpx_refractive_index);
  cpx_permittivity = cpx_rel_permittivity * mp->permittivity;

  // assumed to be constant in an element block
  r_elperm = creal(cpx_permittivity);
  i_elperm = cimag(cpx_permittivity);
  re_coeff = omega * omega * mag_permeability * r_elperm;
  im_coeff = omega * omega * mag_permeability * i_elperm;

  int reqn, ieqn;
  double radvection_etm, rdiffusion_etm, rsource_etm;
  double iadvection_etm, idiffusion_etm, isource_etm;
  // double rmass_etm;
  // double imass_etm;
  double stab_scale = 4.0;
  if (af->Assemble_Residual) {
    for (int i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      for (int a = 0; a < DIM; a++) {
        reqn = R_EM_E1_REAL + a;
        ieqn = R_EM_E1_IMAG + a;
        // rmass_etm = pd->etm[pg->imtrx][reqn][(LOG2_MASS)];
        radvection_etm = pd->etm[pg->imtrx][reqn][(LOG2_ADVECTION)];
        rdiffusion_etm = pd->etm[pg->imtrx][reqn][(LOG2_DIFFUSION)];
        rsource_etm = pd->etm[pg->imtrx][reqn][(LOG2_SOURCE)];
        // imass_etm = pd->etm[pg->imtrx][ieqn][(LOG2_MASS)];
        iadvection_etm = pd->etm[pg->imtrx][ieqn][(LOG2_ADVECTION)];
        idiffusion_etm = pd->etm[pg->imtrx][ieqn][(LOG2_DIFFUSION)];
        isource_etm = pd->etm[pg->imtrx][ieqn][(LOG2_SOURCE)];

        if (pd->e[pg->imtrx][R_EM_E1_REAL + a]) {
          int eqn_real = EM_E1_REAL + a;
          int eqn_imag = EM_E1_IMAG + a;
          int peqn_real = upd->ep[pg->imtrx][eqn_real];
          int peqn_imag = upd->ep[pg->imtrx][eqn_imag];

          double diffusion_real = 0.0;
          double diffusion_imag = 0.0;
          for (int q = 0; q < DIM; q++) {
            diffusion_real += bf[eqn_real]->curl_phi_e[i][a][q] * fv->curl_em_er[q];
            diffusion_imag += bf[eqn_imag]->curl_phi_e[i][a][q] * fv->curl_em_ei[q];
            // diffusion_real += bf->[eqn_real]->ref_phi_e[i][a][q] * fv->em_er[q];
            // diffusion_imag += bf->[eqn_imag]->ref_phi_e[i][a][q] * fv->em_ei[q];

            // R = curl(E)
            // diffusion_real += delta(a,q)*bf[eqn_real]->phi[i]*fv->curl_em_er[a];
            // diffusion_imag += delta(a,q)*bf[eqn_imag]->phi[i]*fv->curl_em_ei[a];
          }
          // R = E
          // diffusion_real = bf[eqn_real]->phi[i]*fv->em_er[a];
          // diffusion_imag = bf[eqn_imag]->phi[i]*fv->em_ei[a];
          double stab_real = 0.0;
          double stab_imag = 0.0;
          double lagr_term_real = 0.0;
          double lagr_term_imag = 0.0;

          if (pd->gv[R_EM_CONT_REAL] && pd->gv[R_EM_CONT_IMAG]) {
            lagr_term_real =
                bf[eqn_real]->grad_phi[i][a] * (r_elperm * fv->epr - i_elperm * fv->epi);
            lagr_term_imag =
                bf[eqn_real]->grad_phi[i][a] * (r_elperm * fv->epi + i_elperm * fv->epr);
            lagr_term_real = bf[eqn_real]->grad_phi[i][a] * (fv->epr);
            lagr_term_imag = bf[eqn_real]->grad_phi[i][a] * (fv->epi);
          } else {
            stab_real = stab_scale * bf[eqn_real]->grad_phi[i][a] *
                        (fv->grad_em_er[0][0] + fv->grad_em_er[1][1] + fv->grad_em_er[2][2]);
            stab_imag = stab_scale * bf[eqn_real]->grad_phi[i][a] *
                        (fv->grad_em_ei[0][0] + fv->grad_em_ei[1][1] + fv->grad_em_ei[2][2]);
          }
          // double advection_real = 0.0;
          // double advection_imag = 0.0;
          double advection_real =
              -bf[eqn_real]->phi[i] * (re_coeff * fv->em_er[a] - im_coeff * fv->em_ei[a]);
          double advection_imag =
              -bf[eqn_imag]->phi[i] * (re_coeff * fv->em_ei[a] + im_coeff * fv->em_er[a]);
          double src_real = 0.0;
          double src_imag = 0.0;
          switch (a) {
          case 0:
            // src_term = bf[eqn_real]->phi[i] * (-fv->x[1]*fv->x[1]*r_elperm*r_elperm - 3);
            //- stab_scale * bf[eqn_real]->grad_phi[i][a] * ( -fv->x[0]);

            // src term for R = E
            // src_term = -bf[eqn_real]->phi[i] * (fv->x[1]*fv->x[1]);

            // src term for R = curl(E)

            break;
          case 1:
            // src_term = bf[eqn_real]->phi[i] * (fv->x[0]*fv->x[1]*r_elperm*r_elperm);
            //- stab_scale * bf[eqn_real]->grad_phi[i][a] * ( -fv->x[0]);

            // src term for R = E
            // src_term = bf[eqn_real]->phi[i] * (fv->x[0]*fv->x[1]);

            // src term for R = curl(E)

            break;
          case 2:
            // src_term = - stab_scale * bf[eqn_real]->grad_phi[i][a] * ( -fv->x[0]);

            // src term for R = curl(E)
            // src_term = bf[eqn_real]->phi[i]*3.0*fv->x[1];

            break;
          }
          lec->R[LEC_R_INDEX(peqn_real, i)] +=
              (advection_real * radvection_etm + diffusion_real * rdiffusion_etm + stab_real +
               src_real * rsource_etm + lagr_term_real) *
              bf[eqn_real]->detJ * fv->wt * fv->h3;

          lec->R[LEC_R_INDEX(peqn_imag, i)] +=
              (advection_imag * iadvection_etm + diffusion_imag * idiffusion_etm + stab_imag +
               src_imag * isource_etm + lagr_term_imag) *
              bf[eqn_imag]->detJ * fv->wt * fv->h3;
        }
      }
    }
  }
  if (af->Assemble_Jacobian) {
    for (int i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      for (int a = 0; a < DIM; a++) {
        reqn = R_EM_E1_REAL + a;
        ieqn = R_EM_E1_IMAG + a;
        // rmass_etm = pd->etm[pg->imtrx][reqn][(LOG2_MASS)];
        radvection_etm = pd->etm[pg->imtrx][reqn][(LOG2_ADVECTION)];
        rdiffusion_etm = pd->etm[pg->imtrx][reqn][(LOG2_DIFFUSION)];
        rsource_etm = pd->etm[pg->imtrx][reqn][(LOG2_SOURCE)];
        // imass_etm = pd->etm[pg->imtrx][ieqn][(LOG2_MASS)];
        iadvection_etm = pd->etm[pg->imtrx][ieqn][(LOG2_ADVECTION)];
        idiffusion_etm = pd->etm[pg->imtrx][ieqn][(LOG2_DIFFUSION)];
        isource_etm = pd->etm[pg->imtrx][ieqn][(LOG2_SOURCE)];
        if (pd->e[pg->imtrx][R_EM_E1_REAL + a]) {
          int eqn_real = EM_E1_REAL + a;
          int eqn_imag = EM_E1_IMAG + a;
          int peqn_real = upd->ep[pg->imtrx][eqn_real];
          int peqn_imag = upd->ep[pg->imtrx][eqn_imag];

          // Sensitivity to real parts of electric field
          for (int b = 0; b < DIM; b++) {
            if (pd->e[pg->imtrx][R_EM_E1_REAL + b]) {
              int var = EM_E1_REAL + b;
              int pvar_real = upd->vp[pg->imtrx][var];
              for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                double diffusion_real = 0;
                double diffusion_imag = 0;
                for (int q = 0; q < DIM; q++) {
                  diffusion_real +=
                      bf[eqn_real]->curl_phi_e[i][a][q] * bf[eqn_real]->curl_phi_e[j][b][q];
                }
                // R = curl(E)
                // diffusion_real += bf[eqn_real]->phi[i]*bf[eqn_real]->curl_phi_e[j][b][a];

                // Jacobian contribution for R = E
                // diffusion_real += delta(a,b)*bf[eqn_real]->phi[i]*bf[eqn_real]->phi[j];

                double stab_real = 0;
                double stab_imag = 0;
                if (!(pd->gv[R_EM_CONT_REAL] && pd->gv[R_EM_CONT_IMAG])) {
                  stab_real = stab_scale * bf[eqn_real]->grad_phi[i][a] *
                              (bf[var]->grad_phi_e[j][b][0][0] + bf[var]->grad_phi_e[j][b][1][1] +
                               bf[var]->grad_phi_e[j][b][2][2]);
                }

                // double advection_real = 0;
                // double advection_imag = 0;
                double advection_real =
                    -bf[eqn_real]->phi[i] * (delta(a, b) * re_coeff * bf[var]->phi[j]);
                double advection_imag =
                    -bf[eqn_imag]->phi[i] * (delta(a, b) * im_coeff * bf[var]->phi[j]);

                lec->J[LEC_J_INDEX(peqn_real, pvar_real, i, j)] +=
                    (diffusion_real * rdiffusion_etm + advection_real * radvection_etm +
                     stab_real) *
                    bf[eqn_real]->detJ * fv->wt * fv->h3;

                lec->J[LEC_J_INDEX(peqn_imag, pvar_real, i, j)] +=
                    (diffusion_imag * idiffusion_etm + advection_imag * iadvection_etm +
                     stab_imag) *
                    bf[eqn_imag]->detJ * fv->wt * fv->h3;
              }
            }
          }

          // Sensitivity to imaginary parts of electric field
          for (int b = 0; b < DIM; b++) {
            if (pd->e[pg->imtrx][R_EM_E1_REAL + b]) {
              int var = EM_E1_IMAG + b;
              int pvar_imag = upd->vp[pg->imtrx][var];
              for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                double diffusion_real = 0;
                double diffusion_imag = 0;
                for (int q = 0; q < DIM; q++) {
                  diffusion_imag +=
                      bf[eqn_imag]->curl_phi_e[i][a][q] * bf[var]->curl_phi_e[j][b][q];
                }

                // R = curl(E)
                // diffusion_imag += bf[eqn_imag]->phi[i]*bf[eqn_imag]->curl_phi_e[j][b][a];

                // Jacobian Contribution for R = E
                // diffusion_imag += delta(a,b)*bf[eqn_imag]->phi[i]*bf[eqn_imag]->phi[j];

                // Jacobian Contribution for R = curl(E)
                // diffusion_real += bf[eqn_real]->phi[i]*bf[eqn_real]->curl_phi_e[j][b][a];

                double stab_real = 0;
                double stab_imag = 0;
                if (!(pd->gv[R_EM_CONT_REAL] && pd->gv[R_EM_CONT_IMAG])) {
                  stab_imag = stab_scale * bf[eqn_real]->grad_phi[i][a] *
                              (bf[var]->grad_phi_e[j][b][0][0] + bf[var]->grad_phi_e[j][b][1][1] +
                               bf[var]->grad_phi_e[j][b][2][2]);
                }

                // double advection_real = 0;
                // double advection_imag = 0;
                double advection_real =
                    bf[eqn_real]->phi[i] * (delta(a, b) * im_coeff * bf[var]->phi[j]);
                double advection_imag =
                    -bf[eqn_imag]->phi[i] * (delta(a, b) * re_coeff * bf[var]->phi[j]);
                lec->J[LEC_J_INDEX(peqn_real, pvar_imag, i, j)] +=
                    (diffusion_real * rdiffusion_etm + advection_real * radvection_etm +
                     stab_real) *
                    bf[eqn_real]->detJ * fv->wt * fv->h3;

                lec->J[LEC_J_INDEX(peqn_imag, pvar_imag, i, j)] +=
                    (diffusion_imag * idiffusion_etm + advection_imag * iadvection_etm +
                     stab_imag) *
                    bf[eqn_imag]->detJ * fv->wt * fv->h3;
              }
            }
          }

          {
            int var = EM_CONT_REAL;
            if (pd->e[pg->imtrx][var]) {
              int pvar = upd->vp[pg->imtrx][var];
              for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                // dbl lagr_term_real = bf[eqn]->grad_phi[i][a] * r_elperm * bf[var]->phi[j];
                // dbl lagr_term_imag = bf[eqn]->grad_phi[i][a] * i_elperm * bf[var]->phi[j];
                dbl lagr_term_real = bf[eqn]->grad_phi[i][a] * bf[var]->phi[j];
                dbl lagr_term_imag = 0;
                lec->J[LEC_J_INDEX(peqn_real, pvar, i, j)] +=
                    (lagr_term_real)*bf[eqn_real]->detJ * fv->wt * fv->h3;
                lec->J[LEC_J_INDEX(peqn_imag, pvar, i, j)] +=
                    (lagr_term_imag)*bf[eqn_real]->detJ * fv->wt * fv->h3;
              }
            }
          }

          {
            int var = EM_CONT_IMAG;
            if (pd->e[pg->imtrx][var]) {
              int pvar = upd->vp[pg->imtrx][var];
              for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                // dbl lagr_term_real = bf[eqn]->grad_phi[i][a] * i_elperm * bf[var]->phi[j];
                // dbl lagr_term_imag = bf[eqn]->grad_phi[i][a] * r_elperm * bf[var]->phi[j];
                dbl lagr_term_real = 0;
                dbl lagr_term_imag = bf[eqn]->grad_phi[i][a] * bf[var]->phi[j];
                lec->J[LEC_J_INDEX(peqn_real, pvar, i, j)] +=
                    (lagr_term_real)*bf[eqn_real]->detJ * fv->wt * fv->h3;
                lec->J[LEC_J_INDEX(peqn_imag, pvar, i, j)] +=
                    (lagr_term_imag)*bf[eqn_real]->detJ * fv->wt * fv->h3;
              }
            }
          }
        }
      }
    }
  }
  return (0);
} // end of assemble_ewave_curlcurl

/* assemble_ewave_laplacian -- assemble terms (Residual &| Jacobian) for EM harmonic
 *                   wave equations
 *                   Substitue H for curl(E) so that goma solves
 *                   Assume div(E) = 0 and use identity curl(curl(A) = grad(div(A)) - Laplacian(A)
 *                   1 complex vector equation and 1 complex vector primitive
 *
 * in:
 *     ei -- pointer to Element Indecesstructure
 *     pd -- pointer to Problem Descriptionstructure
 *     af -- pointer to Action Flagstructure
 *     bf -- pointer to Basis Functionstructure
 *     fv -- pointer to Field Variablestructure
 *       fv_old -- pointer to old Diet Field Variablestructure
 *       fv_dot -- pointer to dot Diet Field Variablestructure
 *     cr -- pointer to Constitutive Relationstructure
 *     md -- pointer to Mesh Derivativestructure
 *     me -- pointer to Material Entitystructure
 *
 * out:
 *     a   -- gets loaded up with proper contribution
 *     lec -- gets loaded up with local contributions to resid, Jacobian
 *     r   -- residual RHS vector
 *
 * Created: Wednesday April 22, 2020 - Andrew Cochrane
 *
 *
 */

int assemble_ewave_laplacian(double time,      // present time
                             double tt,        // time integration method parameter
                             double dt,        // current time step size
                             const int em_eqn, // eqn id
                             const int em_var  //  variable id - should match me_eqn
) {
  int eqn;
  dbl mag_permeability = mp->magnetic_permeability;
  double omega, re_coeff, im_coeff;
  /*
  struct emwave_stabilization em_stab;
  em_stab.em_eqn = em_eqn;
  em_stab.em_var = em_var;
  em_stab.type = EM_STAB_NONE; // enum supports phi_div, dphi_div,
                           // divphi_div, phi_divsquared and
                           // dphi_divsquared
  */
  eqn = em_eqn;
  /*
   * Bail out fast if there's nothing to do...
   * But we might have the wrong eqn
   */
  for (int a = 0; a < DIM; a++) {
    if (pd->e[pg->imtrx][R_EM_E1_REAL + a]) {
      eqn = R_EM_E1_REAL + a;
      break;
    }
  }
  if (!pd->e[pg->imtrx][eqn]) {
    return (1);
  }

  omega = upd->EM_Frequency;
  dbl n; /* Refractive index. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_n_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_n = &d_n_struct;

  dbl k; /* Acoustic wave number. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  n = refractive_index(d_n, time);
  k = extinction_index(d_k, time);

  // Compute complex material properties
  complex double cpx_refractive_index, cpx_rel_permittivity,
      cpx_permittivity; //, impedance;
  double r_elperm, i_elperm;

  cpx_refractive_index = n + _Complex_I * k; // k > 0 is extinction
  cpx_rel_permittivity = SQUARE(cpx_refractive_index);
  cpx_permittivity = cpx_rel_permittivity * mp->permittivity;

  // assumed to be constant in an element block
  r_elperm = creal(cpx_permittivity);
  i_elperm = cimag(cpx_permittivity);
  re_coeff = omega * omega * mag_permeability * r_elperm;
  im_coeff = omega * omega * mag_permeability * i_elperm;

  double stab_scale = 0.0;
  if (af->Assemble_Residual) {
    for (int i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      for (int a = 0; a < DIM; a++) {
        if (pd->e[pg->imtrx][R_EM_E1_REAL + a]) {
          int eqn_real = EM_E1_REAL + a;
          int eqn_imag = EM_E1_IMAG + a;
          int peqn_real = upd->ep[pg->imtrx][eqn_real];
          int peqn_imag = upd->ep[pg->imtrx][eqn_imag];

          double diffusion_real = 0.0;
          double diffusion_imag = 0.0;

          for (int q = 0; q < DIM; q++) {
            // R = - (gradphi)_a dot (gradE)_a
            diffusion_real += -bf[eqn_real]->grad_phi[i][q] * fv->grad_em_er[q][a];
            diffusion_imag += -bf[eqn_imag]->grad_phi[i][q] * fv->grad_em_ei[q][a];
          }
          // diffusion_real = -fv->grad_em_er[a][0];
          //  R = E
          // diffusion_real = bf[eqn_real]->phi[i]*fv->em_er[a];
          // diffusion_imag = bf[eqn_imag]->phi[i]*fv->em_ei[a];
          double stab_real = 0.0;
          double stab_imag = 0.0;
          double lagr_term_real = 0.0;
          double lagr_term_imag = 0.0;

          if (pd->gv[R_EM_CONT_REAL] && pd->gv[R_EM_CONT_IMAG]) {
            // lagr_term_real = bf[eqn_real]->grad_phi[i][a] * (r_elperm * fv->epr + i_elperm *
            // fv->epi); lagr_term_imag = bf[eqn_real]->grad_phi[i][a] * (r_elperm * fv->epi +
            // i_elperm * fv->epr); lagr_term_real = bf[eqn_real]->grad_phi[i][a] * (fv->epr);
            // lagr_term_imag = bf[eqn_real]->grad_phi[i][a] * (fv->epi);
          } else {
            stab_real = stab_scale * bf[eqn_real]->grad_phi[i][a] *
                        (fv->grad_em_er[0][0] + fv->grad_em_er[1][1] + fv->grad_em_er[2][2]);
            stab_imag = stab_scale * bf[eqn_real]->grad_phi[i][a] *
                        (fv->grad_em_ei[0][0] + fv->grad_em_ei[1][1] + fv->grad_em_ei[2][2]);
          }

          double advection_real = 0.0;
          double advection_imag = 0.0;
          advection_real =
              bf[eqn_real]->phi[i] * (re_coeff * fv->em_er[a] - im_coeff * fv->em_ei[a]);
          advection_imag =
              bf[eqn_imag]->phi[i] * (re_coeff * fv->em_ei[a] + im_coeff * fv->em_er[a]);

          double src_term = 0.0;
          switch (a) {
          case 0:
            // src_term = bf[eqn_real]->phi[i] * (-fv->x[1]*fv->x[1]*r_elperm*r_elperm - 3);
            //- stab_scale * bf[eqn_real]->grad_phi[i][a] * ( -fv->x[0]);

            // src term for R = E
            // src_term = -bf[eqn_real]->phi[i] * (fv->x[1]*fv->x[1]);

            // src term for R = curl(E)

            break;
          case 1:
            // src_term = bf[eqn_real]->phi[i] * (fv->x[0]*fv->x[1]*r_elperm*r_elperm);
            //- stab_scale * bf[eqn_real]->grad_phi[i][a] * ( -fv->x[0]);

            // src term for R = E
            // src_term = bf[eqn_real]->phi[i] * (fv->x[0]*fv->x[1]);

            // src term for R = curl(E)

            break;
          case 2:
            // src_term = - stab_scale * bf[eqn_real]->grad_phi[i][a] * ( -fv->x[0]);

            // src term for R = curl(E)
            // src_term = bf[eqn_real]->phi[i]*3.0*fv->x[1];

            break;
          }
          lec->R[LEC_R_INDEX(peqn_real, i)] +=
              (advection_real + diffusion_real + stab_real + src_term + lagr_term_real) *
              bf[eqn_real]->detJ * fv->wt * fv->h3;
          lec->R[LEC_R_INDEX(peqn_imag, i)] +=
              (advection_imag + diffusion_imag + stab_imag + lagr_term_imag) * bf[eqn_imag]->detJ *
              fv->wt * fv->h3;
        }
      }
    }
  }
  if (af->Assemble_Jacobian) {
    for (int i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      for (int a = 0; a < DIM; a++) {
        if (pd->e[pg->imtrx][R_EM_E1_REAL + a]) {
          int eqn_real = EM_E1_REAL + a;
          int eqn_imag = EM_E1_IMAG + a;
          int peqn_real = upd->ep[pg->imtrx][eqn_real];
          int peqn_imag = upd->ep[pg->imtrx][eqn_imag];

          // Sensitivity to real parts of electric field
          for (int b = 0; b < DIM; b++) {
            if (pd->e[pg->imtrx][R_EM_E1_REAL + b]) {
              int var = EM_E1_REAL + b;
              int pvar_real = upd->vp[pg->imtrx][var];
              for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                double diffusion_real = 0;
                double diffusion_imag = 0;
                for (int q = 0; q < DIM; q++) {
                  // diffusion_real += bf[eqn_real]->curl_phi_e[i][a][q] *
                  // bf[eqn_real]->curl_phi_e[j][b][q];
                  diffusion_real +=
                      -delta(a, b) * bf[eqn_real]->grad_phi[i][q] * bf[eqn_real]->grad_phi[j][q];
                }

                // diffusion_real += -bf[eqn_real]->grad_phi[j][b];

                // R = curl(E)
                // diffusion_real += bf[eqn_real]->phi[i]*bf[eqn_real]->curl_phi_e[j][b][a];

                // Jacobian contribution for R = E
                // diffusion_real += delta(a,b)*bf[eqn_real]->phi[i]*bf[eqn_real]->phi[j];

                double stab_real = 0;
                double stab_imag = 0;
                if (!(pd->gv[R_EM_CONT_REAL] && pd->gv[R_EM_CONT_IMAG])) {
                  stab_real = stab_scale * bf[eqn_real]->grad_phi[i][a] *
                              (bf[var]->grad_phi_e[j][b][0][0] + bf[var]->grad_phi_e[j][b][1][1] +
                               bf[var]->grad_phi_e[j][b][2][2]);
                }

                double advection_real = 0;
                double advection_imag = 0;

                advection_real = bf[eqn_real]->phi[i] * (delta(a, b) * re_coeff * bf[var]->phi[j]);
                advection_imag = bf[eqn_imag]->phi[i] * (delta(a, b) * im_coeff * bf[var]->phi[j]);

                lec->J[LEC_J_INDEX(peqn_real, pvar_real, i, j)] +=
                    (diffusion_real + advection_real + stab_real) * bf[eqn_real]->detJ * fv->wt *
                    fv->h3;
                lec->J[LEC_J_INDEX(peqn_imag, pvar_real, i, j)] +=
                    (diffusion_imag + advection_imag + stab_imag) * bf[eqn_imag]->detJ * fv->wt *
                    fv->h3;
              }
            }
          }

          // Sensitivity to imaginary parts of electric field
          for (int b = 0; b < DIM; b++) {
            if (pd->e[pg->imtrx][R_EM_E1_REAL + b]) {
              int var = EM_E1_IMAG + b;
              int pvar_imag = upd->vp[pg->imtrx][var];
              for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                double diffusion_real = 0;
                double diffusion_imag = 0;
                for (int q = 0; q < DIM; q++) {
                  diffusion_imag +=
                      -delta(a, b) * bf[eqn_imag]->grad_phi[i][q] * bf[eqn_imag]->grad_phi[j][q];
                }

                // R = curl(E)
                // diffusion_imag += bf[eqn_imag]->phi[i]*bf[eqn_imag]->curl_phi_e[j][b][a];

                // Jacobian Contribution for R = E
                // diffusion_imag += delta(a,b)*bf[eqn_imag]->phi[i]*bf[eqn_imag]->phi[j];

                // Jacobian Contribution for R = curl(E)
                // diffusion_real += bf[eqn_real]->phi[i]*bf[eqn_real]->curl_phi_e[j][b][a];

                double stab_real = 0;
                double stab_imag = 0;
                if (!(pd->gv[R_EM_CONT_REAL] && pd->gv[R_EM_CONT_IMAG])) {
                  stab_imag = stab_scale * bf[eqn_real]->grad_phi[i][a] *
                              (bf[var]->grad_phi_e[j][b][0][0] + bf[var]->grad_phi_e[j][b][1][1] +
                               bf[var]->grad_phi_e[j][b][2][2]);
                }

                double advection_real = 0;
                double advection_imag = 0;

                advection_real = bf[eqn_real]->phi[i] * (delta(a, b) * -im_coeff * bf[var]->phi[j]);
                advection_imag = bf[eqn_imag]->phi[i] * (delta(a, b) * re_coeff * bf[var]->phi[j]);
                lec->J[LEC_J_INDEX(peqn_real, pvar_imag, i, j)] +=
                    (diffusion_real + advection_real + stab_real) * bf[eqn_real]->detJ * fv->wt *
                    fv->h3;
                lec->J[LEC_J_INDEX(peqn_imag, pvar_imag, i, j)] +=
                    (diffusion_imag + advection_imag + stab_imag) * bf[eqn_imag]->detJ * fv->wt *
                    fv->h3;
              }
            }
          }

          {
            int var = EM_CONT_REAL;
            if (pd->e[pg->imtrx][var]) {
              int pvar = upd->vp[pg->imtrx][var];
              for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                // dbl lagr_term_real = bf[eqn]->grad_phi[i][a] * r_elperm * bf[var]->phi[j];
                // dbl lagr_term_imag = bf[eqn]->grad_phi[i][a] * i_elperm * bf[var]->phi[j];
                dbl lagr_term_real = bf[eqn]->grad_phi[i][a] * bf[var]->phi[j];
                dbl lagr_term_imag = 0;
                lec->J[LEC_J_INDEX(peqn_real, pvar, i, j)] +=
                    (lagr_term_real)*bf[eqn_real]->detJ * fv->wt * fv->h3;
                lec->J[LEC_J_INDEX(peqn_imag, pvar, i, j)] +=
                    (lagr_term_imag)*bf[eqn_real]->detJ * fv->wt * fv->h3;
              }
            }
          }

          {
            int var = EM_CONT_IMAG;
            if (pd->e[pg->imtrx][var]) {
              int pvar = upd->vp[pg->imtrx][var];
              for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                // dbl lagr_term_real = bf[eqn]->grad_phi[i][a] * i_elperm * bf[var]->phi[j];
                // dbl lagr_term_imag = bf[eqn]->grad_phi[i][a] * r_elperm * bf[var]->phi[j];
                dbl lagr_term_real = 0;
                dbl lagr_term_imag = bf[eqn]->grad_phi[i][a] * bf[var]->phi[j];
                lec->J[LEC_J_INDEX(peqn_real, pvar, i, j)] +=
                    (lagr_term_real)*bf[eqn_real]->detJ * fv->wt * fv->h3;
                lec->J[LEC_J_INDEX(peqn_imag, pvar, i, j)] +=
                    (lagr_term_imag)*bf[eqn_real]->detJ * fv->wt * fv->h3;
              }
            }
          }
        }
      }
    }
  }
  return (0);
} // end of assemble_ewave_laplacian

int assemble_em_continuity() {
  if (!(pd->e[pg->imtrx][R_EM_CONT_REAL] && pd->e[pg->imtrx][R_EM_CONT_IMAG])) {
    return -1;
  }

  dbl n; /* Refractive index. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_n_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_n = &d_n_struct;

  dbl k; /* Acoustic wave number. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  n = refractive_index(d_n, 0);
  k = extinction_index(d_k, 0);

  // Compute complex material properties
  complex double cpx_refractive_index, cpx_rel_permittivity,
      cpx_permittivity; //, impedance;
  double r_elperm, i_elperm;

  cpx_refractive_index = n + _Complex_I * k; // k > 0 is extinction
  cpx_rel_permittivity = SQUARE(cpx_refractive_index);
  cpx_permittivity = cpx_rel_permittivity * mp->permittivity;

  // assumed to be constant in an element block
  r_elperm = creal(cpx_permittivity);
  i_elperm = cimag(cpx_permittivity);

  dbl div_E_real = 0;
  dbl div_E_imag = 0;

  dbl dim = pd->Num_Dim;
  double factor = 1e3;
  for (int i = 0; i < dim; i++) {
    div_E_real += fv->grad_em_er[i][i];
    div_E_imag += fv->grad_em_ei[i][i];
  }
  if (af->Assemble_Residual) {
    for (int i = 0; i < ei[pg->imtrx]->dof[R_EM_CONT_REAL]; i++) {
      int eqn_real = R_EM_CONT_REAL;
      int eqn_imag = R_EM_CONT_IMAG;
      int peqn_real = upd->ep[pg->imtrx][eqn_real];
      int peqn_imag = upd->ep[pg->imtrx][eqn_imag];

      double advection_real =
          -factor * bf[eqn_real]->phi[i] * (r_elperm * div_E_real + i_elperm * div_E_imag);
      double advection_imag =
          -factor * bf[eqn_imag]->phi[i] * (r_elperm * div_E_imag + i_elperm * div_E_real);

      lec->R[LEC_R_INDEX(peqn_real, i)] += (advection_real)*bf[eqn_real]->detJ * fv->wt * fv->h3;
      lec->R[LEC_R_INDEX(peqn_imag, i)] += (advection_imag)*bf[eqn_imag]->detJ * fv->wt * fv->h3;
    }
  }
  if (af->Assemble_Jacobian) {
    for (int i = 0; i < ei[pg->imtrx]->dof[R_EM_CONT_REAL]; i++) {
      int eqn_real = EM_CONT_REAL;
      int eqn_imag = EM_CONT_IMAG;
      int peqn_real = upd->ep[pg->imtrx][eqn_real];
      int peqn_imag = upd->ep[pg->imtrx][eqn_imag];

      for (int b = 0; b < dim; b++) {
        int var = EM_E1_REAL + b;
        int pvar_real = upd->ep[pg->imtrx][var];
        for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          double advection_real = 0;
          double advection_imag = 0;
          for (int q = 0; q < dim; q++) {
            advection_real +=
                -factor * bf[eqn_real]->phi[i] * r_elperm * bf[var]->grad_phi_e[j][b][q][q];
            advection_imag +=
                -factor * bf[eqn_imag]->phi[i] * i_elperm * bf[var]->grad_phi_e[j][b][q][q];
          }

          lec->J[LEC_J_INDEX(peqn_real, pvar_real, i, j)] +=
              (advection_real)*bf[eqn_real]->detJ * fv->wt * fv->h3;
          lec->J[LEC_J_INDEX(peqn_imag, pvar_real, i, j)] +=
              (advection_imag)*bf[eqn_imag]->detJ * fv->wt * fv->h3;
        }
      }

      for (int b = 0; b < dim; b++) {
        int var = EM_E1_IMAG + b;
        int pvar_imag = upd->ep[pg->imtrx][var];
        for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          double advection_real = 0;
          double advection_imag = 0;
          for (int q = 0; q < dim; q++) {
            advection_real +=
                -factor * bf[eqn_real]->phi[i] * i_elperm * bf[var]->grad_phi_e[j][b][q][q];
            advection_imag +=
                -factor * bf[eqn_imag]->phi[i] * r_elperm * bf[var]->grad_phi_e[j][b][q][q];
          }

          lec->J[LEC_J_INDEX(peqn_real, pvar_imag, i, j)] +=
              (advection_real)*bf[eqn_real]->detJ * fv->wt * fv->h3;
          lec->J[LEC_J_INDEX(peqn_imag, pvar_imag, i, j)] +=
              (advection_imag)*bf[eqn_imag]->detJ * fv->wt * fv->h3;
        }
      }
    }
  }
  return (0);
} // end of assemble_em_continuity

void calc_emwave_stabilization_term(struct emwave_stabilization *em_stab,
                                    double stabilization_coefficient) {
  double complex grad_stabilization_field[DIM][DIM] = {{0.0}};
  dbl mag_permeability = 12.57e-07; // H/m

  /*
   *
   */
  switch (em_stab->type) {
  case EM_STAB_NONE:
    for (int i = 0; i < ei[pg->imtrx]->dof[em_stab->em_eqn]; i++) {
      em_stab->residual_term[i] = 0.0;
      for (int b = 0; b < DIM; b++) {
        for (int j = 0; j < ei[pg->imtrx]->dof[em_stab->em_var]; j++) {
          em_stab->jacobian_term[i][b][j] = 0.0;
        }
      }
    }
    return;
    break;
  case EM_STAB_PHI_DIV:
  case EM_STAB_DPHI_DIV:
  case EM_STAB_DIVPHI_DIV:
    switch (em_stab->stabilization_field_var) {
    case EM_E1_REAL:
      for (int p = 0; p < VIM; p++) {
        for (int q = 0; q < VIM; q++) {
          grad_stabilization_field[p][q] = fv->grad_em_er[p][q];
        }
      }
      break;
    case EM_E1_IMAG:
      for (int p = 0; p < VIM; p++) {
        for (int q = 0; q < VIM; q++) {
          grad_stabilization_field[p][q] = fv->grad_em_ei[p][q];
        }
      }
      break;
    case EM_H1_REAL:
      for (int p = 0; p < VIM; p++) {
        for (int q = 0; q < VIM; q++) {
          grad_stabilization_field[p][q] = fv->grad_em_hr[p][q];
        }
      }
      break;
    case EM_H1_IMAG:
      for (int p = 0; p < VIM; p++) {
        for (int q = 0; q < VIM; q++) {
          grad_stabilization_field[p][q] = fv->grad_em_hi[p][q];
        }
      }
      break;
    default:
      GOMA_EH(GOMA_ERROR, "Cannot have unset stabilization_field_var");
      break;
    }
    break;

  case EM_STAB_PHI_DIVSQUARED:
  case EM_STAB_DPHI_DIVSQUARED:
    switch (em_stab->stabilization_field_var) {
    case EM_H1_REAL:
    case EM_H1_IMAG:
      for (int p = 0; p < VIM; p++) {
        for (int q = 0; q < VIM; q++) {
          grad_stabilization_field[p][q] = fv->grad_em_er[p][q] + _Complex_I * fv->grad_em_ei[p][q];
        }
      }
      break;
    case EM_E1_REAL:
    case EM_E1_IMAG:
      for (int p = 0; p < VIM; p++) {
        for (int q = 0; q < VIM; q++) {
          grad_stabilization_field[p][q] = fv->grad_em_hr[p][q] + _Complex_I * fv->grad_em_hi[p][q];
        }
      }
      break;
    default:
      GOMA_EH(GOMA_ERROR, "Must use the first EQN_VAR for a given 3-valued field, e.g. EM_E1_REAL");
    }
    break;

  default:
    // WH(-1,"Cannot use calc_emwave_stabilization without defining type in the struct");
    // return -1;
    break;
  }
  int cartesian_index;
  double complex div_stabilization_field;
  double complex div_stabilization_field_squared;
  switch (em_stab->type) {
  case EM_STAB_PHI_DIV:

    div_stabilization_field = 0.0;

    for (int p = 0; p < VIM; p++) {
      div_stabilization_field += grad_stabilization_field[p][p];
    }

    break;

  case EM_STAB_DPHI_DIV:

    div_stabilization_field = 0.0;

    // need the index that corresponds to the x, y or z
    // of the current residual
    cartesian_index = (em_stab->em_eqn - R_EM_E1_REAL) % 3;

    for (int p = 0; p < VIM; p++) {
      div_stabilization_field += grad_stabilization_field[p][p];
    }

    for (int i = 0; i < ei[pg->imtrx]->dof[em_stab->em_eqn]; i++) {
      em_stab->residual_term[i] = bf[em_stab->em_eqn]->grad_phi[i][cartesian_index] *
                                  stabilization_coefficient * creal(div_stabilization_field);

      for (int b = 0; b < DIM; b++) {
        for (int j = 0; j < ei[pg->imtrx]->dof[em_stab->em_var]; j++) {
          em_stab->jacobian_term[i][b][j] =
              bf[em_stab->em_eqn]->grad_phi[i][cartesian_index] * stabilization_coefficient *
              bf[em_stab->stabilization_field_var + b]->grad_phi[j][b];
        }
      }
    }
    break;

  case EM_STAB_DIVPHI_DIV:
    div_stabilization_field = 0.0;

    double div_phi[MDE] = {0.0};

    for (int p = 0; p < VIM; p++) {
      div_stabilization_field += grad_stabilization_field[p][p];
    }
    for (int i = 0; i < ei[pg->imtrx]->dof[em_stab->em_eqn]; i++) {
      for (int p = 0; p < DIM; p++) {
        div_phi[i] += bf[em_stab->em_eqn]->grad_phi[i][p];
      }
    }

    for (int i = 0; i < ei[pg->imtrx]->dof[em_stab->em_eqn]; i++) {
      em_stab->residual_term[i] =
          div_phi[i] * stabilization_coefficient * creal(div_stabilization_field);

      for (int b = 0; b < DIM; b++) {
        for (int j = 0; j < ei[pg->imtrx]->dof[em_stab->em_var]; j++) {
          em_stab->jacobian_term[i][b][j] =
              div_phi[i] * stabilization_coefficient *
              bf[em_stab->stabilization_field_var + b]->grad_phi[j][b];
        }
      }
    }
    break;

  case EM_STAB_PHI_DIVSQUARED:
    div_stabilization_field = 0.0;

    for (int p = 0; p < VIM; p++) {
      div_stabilization_field += grad_stabilization_field[p][p];
    }
    div_stabilization_field_squared = div_stabilization_field * div_stabilization_field;
    switch (em_stab->stabilization_field_var) {
    case EM_E1_REAL:
      for (int i = 0; i < ei[pg->imtrx]->dof[em_stab->em_eqn]; i++) {
        em_stab->residual_term[i] = bf[em_stab->em_eqn]->phi[i] * stabilization_coefficient *
                                    creal(div_stabilization_field_squared);
        //*mp->permittivity;
        for (int b = 0; b < DIM; b++) {
          for (int j = 0; j < ei[pg->imtrx]->dof[em_stab->em_var]; j++) {
            em_stab->jacobian_term[i][b][j] =
                bf[em_stab->em_eqn]->phi[i] * stabilization_coefficient * 2.0 *
                div_stabilization_field * bf[em_stab->stabilization_field_var + b]->grad_phi[j][b];
          }
        }
      }
      break;
    case EM_H1_REAL:
      for (int i = 0; i < ei[pg->imtrx]->dof[em_stab->em_eqn]; i++) {
        em_stab->residual_term[i] = bf[em_stab->em_eqn]->phi[i] * stabilization_coefficient *
                                    creal(div_stabilization_field_squared);
        //*mag_permeability;
        for (int b = 0; b < DIM; b++) {
          for (int j = 0; j < ei[pg->imtrx]->dof[em_stab->em_var]; j++) {
            em_stab->jacobian_term[i][b][j] =
                bf[em_stab->em_eqn]->phi[i] * stabilization_coefficient * 2.0 *
                div_stabilization_field * bf[em_stab->stabilization_field_var + b]->grad_phi[j][b];
          }
        }
      }
      break;
    case EM_E1_IMAG:
      for (int i = 0; i < ei[pg->imtrx]->dof[em_stab->em_eqn]; i++) {
        em_stab->residual_term[i] = bf[em_stab->em_eqn]->phi[i] * stabilization_coefficient *
                                    cimag(div_stabilization_field_squared);
        //*mp->permittivity;
        for (int b = 0; b < DIM; b++) {
          for (int j = 0; j < ei[pg->imtrx]->dof[em_stab->em_var]; j++) {
            em_stab->jacobian_term[i][b][j] =
                bf[em_stab->em_eqn]->phi[i] * stabilization_coefficient * 2.0 *
                div_stabilization_field * bf[em_stab->stabilization_field_var + b]->grad_phi[j][b];
          }
        }
      }
      break;
    case EM_H1_IMAG:
      for (int i = 0; i < ei[pg->imtrx]->dof[em_stab->em_eqn]; i++) {
        em_stab->residual_term[i] = bf[em_stab->em_eqn]->phi[i] * stabilization_coefficient *
                                    cimag(div_stabilization_field_squared);
        //*mag_permeability;
        for (int b = 0; b < DIM; b++) {
          for (int j = 0; j < ei[pg->imtrx]->dof[em_stab->em_var]; j++) {
            em_stab->jacobian_term[i][b][j] =
                bf[em_stab->em_eqn]->phi[i] * stabilization_coefficient * 2.0 *
                div_stabilization_field * bf[em_stab->stabilization_field_var + b]->grad_phi[j][b];
          }
        }
      }
      break;
    default:
      GOMA_EH(GOMA_ERROR, "must set the stabilization field var with type==phi_divsquared");
      return;
      break;
    }
    break;

  case EM_STAB_DPHI_DIVSQUARED:
    div_stabilization_field = 0.0;

    // need the index that corresponds to the x, y or z
    // of the current residual
    cartesian_index = (em_stab->em_eqn - R_EM_E1_REAL) % 3;

    for (int p = 0; p < VIM; p++) {
      div_stabilization_field += grad_stabilization_field[p][p];
    }
    div_stabilization_field_squared = div_stabilization_field * div_stabilization_field;
    switch (em_stab->stabilization_field_var) {
    case EM_E1_REAL:
      for (int i = 0; i < ei[pg->imtrx]->dof[em_stab->em_eqn]; i++) {
        em_stab->residual_term[i] = bf[em_stab->em_eqn]->grad_phi[i][cartesian_index] *
                                    stabilization_coefficient *
                                    creal(div_stabilization_field_squared) / mag_permeability;
        //*mp->permittivity;
        for (int b = 0; b < DIM; b++) {
          for (int j = 0; j < ei[pg->imtrx]->dof[em_stab->em_var]; j++) {
            em_stab->jacobian_term[i][b][j] =
                bf[em_stab->em_eqn]->phi[i] * stabilization_coefficient * 2.0 *
                div_stabilization_field * bf[em_stab->stabilization_field_var + b]->grad_phi[j][b];
          }
        }
      }
      break;
    case EM_H1_REAL:
      for (int i = 0; i < ei[pg->imtrx]->dof[em_stab->em_eqn]; i++) {
        em_stab->residual_term[i] = bf[em_stab->em_eqn]->grad_phi[i][cartesian_index] *
                                    stabilization_coefficient *
                                    creal(div_stabilization_field_squared);
        //*mag_permeability;
        for (int b = 0; b < DIM; b++) {
          for (int j = 0; j < ei[pg->imtrx]->dof[em_stab->em_var]; j++) {
            em_stab->jacobian_term[i][b][j] =
                bf[em_stab->em_eqn]->phi[i] * stabilization_coefficient * 2.0 *
                div_stabilization_field * bf[em_stab->stabilization_field_var + b]->grad_phi[j][b];
          }
        }
      }
      break;
    case EM_E1_IMAG:
      for (int i = 0; i < ei[pg->imtrx]->dof[em_stab->em_eqn]; i++) {
        em_stab->residual_term[i] = bf[em_stab->em_eqn]->grad_phi[i][cartesian_index] *
                                    stabilization_coefficient *
                                    cimag(div_stabilization_field_squared) / mag_permeability;
        //*mp->permittivity;
        for (int b = 0; b < DIM; b++) {
          for (int j = 0; j < ei[pg->imtrx]->dof[em_stab->em_var]; j++) {
            em_stab->jacobian_term[i][b][j] =
                bf[em_stab->em_eqn]->phi[i] * stabilization_coefficient * 2.0 *
                div_stabilization_field * bf[em_stab->stabilization_field_var + b]->grad_phi[j][b];
          }
        }
      }
      break;
    case EM_H1_IMAG:
      for (int i = 0; i < ei[pg->imtrx]->dof[em_stab->em_eqn]; i++) {
        em_stab->residual_term[i] = bf[em_stab->em_eqn]->grad_phi[i][cartesian_index] *
                                    stabilization_coefficient *
                                    cimag(div_stabilization_field_squared);
        //*mag_permeability;
        for (int b = 0; b < DIM; b++) {
          for (int j = 0; j < ei[pg->imtrx]->dof[em_stab->em_var]; j++) {
            em_stab->jacobian_term[i][b][j] =
                bf[em_stab->em_eqn]->phi[i] * stabilization_coefficient * 2.0 *
                div_stabilization_field * bf[em_stab->stabilization_field_var + b]->grad_phi[j][b];
          }
        }
      }
      break;
    default:
      GOMA_EH(GOMA_ERROR, "must set the stabilization field var with type==phi_divsquared");
      return;
      break;
    }
    break;
  case EM_STAB_NONE:
  default:
    for (int i = 0; i < ei[pg->imtrx]->dof[em_stab->em_eqn]; i++) {
      em_stab->residual_term[i] = 0.0;
      for (int b = 0; b < DIM; b++) {
        for (int j = 0; j < ei[pg->imtrx]->dof[em_stab->em_var]; j++) {
          em_stab->jacobian_term[i][b][j] = 0.0;
        }
      }
    }

    break;
  }

  return;
}

/* Cross the first two complex[DIM] vectors and return their
 * cross product.
 * TODO: create a return struct that contains sensitivities to input
 */
void complex_cross_vectors(const complex *v0, /* v0 */
                           const complex *v1, /* v1 */
                           complex *v2)       /* v2 = v0 x v1 */
{
  int i, j, k;

  memset(v2, 0, DIM * sizeof(complex));
  for (i = 0; i < DIM; i++)
    for (j = 0; j < DIM; j++)
      for (k = 0; k < DIM; k++)
        v2[k] += permute(i, j, k) * v0[i] * v1[j];
} // end of complex_cross_vectors

// returns true if permittivity model is matrix form
bool relative_permittivity_model(complex double *permittivity_out,
                                 complex double *permittivity_matrix) {
  complex double permittivity = 0;
  switch (mp->PermittivityModel) {
  case CONSTANT:
    permittivity = mp->permittivity;
    break;
  case COMPLEX_CONSTANT:
    permittivity = mp->permittivity - _Complex_I * mp->permittivity_imag;
    break;
  case REFRACTIVE_INDEX: {
    dbl n = refractive_index(NULL, 0.0);
    dbl k = extinction_index(NULL, 0.0);
    permittivity = cpow((n - _Complex_I * k), 2);
  } break;
  case RADIAL_PML: {
    dbl amp = mp->u_permittivity[0];
    dbl power = mp->u_permittivity[1];
    dbl cond_max = mp->u_permittivity[2];
    dbl pml_inner_radius = mp->u_permittivity[3];
    dbl pml_outer_radius = mp->u_permittivity[4];

    const double c0 = 1.0 / (sqrt(upd->Free_Space_Permittivity * upd->Free_Space_Permeability));
    const double eps0 = upd->Free_Space_Permittivity;
    const double mu0 = upd->Free_Space_Permeability;
    const double nu0 = sqrt(mu0 / eps0);

    dbl x = fv->x[0];
    dbl y = fv->x[1];
    dbl z = fv->x[2];
    dbl freq = upd->EM_Frequency;
    dbl lambda0 = c0 / freq;
    dbl k0 = 2 * M_PI / lambda0;
    dbl d = pml_outer_radius - pml_inner_radius;

    dbl L_mag = sqrt(x * x + y * y + z * z);
    dbl L = MAX(L_mag - pml_inner_radius, 0.0);
    if (L <= 1e-8) {

      permittivity_matrix[0] = 1.0;
      permittivity_matrix[1] = 1.0;
      permittivity_matrix[2] = 1.0;
      return true;
    } else {
      dbl L_c[DIM] = {d, d, d};
      dbl pml_x = x * L / d;
      dbl pml_y = y * L / d;
      dbl pml_z = z * L / d;

      dbl amax[DIM];
      amax[0] = L_c[0] < 1e-8 ? 1 : 1 + amp * pow(fabs(pml_x / L_c[0]), power);
      amax[1] = L_c[1] < 1e-8 ? 1 : 1 + amp * pow(fabs(pml_y / L_c[1]), power);
      amax[2] = L_c[2] < 1e-8 ? 1 : 1 + amp * pow(fabs(pml_z / L_c[2]), power);

      dbl cond[DIM];
      cond[0] = L_c[0] < 1e-8 ? 0
                              : cond_max * sin(fabs(M_PI * pml_x / (2 * L_c[0]))) *
                                    sin(fabs(M_PI * pml_x / (2 * L_c[0])));
      cond[1] = L_c[1] < 1e-8 ? 0
                              : cond_max * sin(fabs(M_PI * pml_y / (2 * L_c[1]))) *
                                    sin(fabs(M_PI * pml_y / (2 * L_c[1])));
      cond[2] = L_c[2] < 1e-8 ? 0
                              : cond_max * sin(fabs(M_PI * pml_z / (2 * L_c[2]))) *
                                    sin(fabs(M_PI * pml_z / (2 * L_c[2])));

      complex double s[DIM] = {amax[0] * (1 + _Complex_I * nu0 * cond[0] / k0),
                               amax[1] * (1 + _Complex_I * nu0 * cond[1] / k0),
                               amax[2] * (1 + _Complex_I * nu0 * cond[2] / k0)};

      permittivity_matrix[0] = s[1] * s[2] / s[0];
      permittivity_matrix[1] = s[0] * s[2] / s[1];
      permittivity_matrix[2] = s[0] * s[1] / s[2];
      return true;
    }
  } break;
  default:
    GOMA_EH(GOMA_ERROR, "Unknown permitivity model");
  }

  *permittivity_out = permittivity;
  return false;
}

// returns true if permittivity model is matrix form
bool relative_permeability_model(complex double *permeability_out,
                                 complex double *permeability_matrix) {
  complex double permeability = 0;
  switch (mp->PermeabilityModel) {
  case CONSTANT:
    permeability = mp->permeability;
    break;
  case COMPLEX_CONSTANT:
    permeability = mp->permeability + _Complex_I * mp->permeability_imag;
    break;

  case RADIAL_PML: {
    dbl amp = mp->u_permeability[0];
    dbl power = mp->u_permeability[1];
    dbl cond_max = mp->u_permeability[2];
    dbl pml_inner_radius = mp->u_permeability[3];
    dbl pml_outer_radius = mp->u_permeability[4];

    const double c0 = 1.0 / (sqrt(upd->Free_Space_Permittivity * upd->Free_Space_Permeability));
    const double eps0 = upd->Free_Space_Permittivity;
    const double mu0 = upd->Free_Space_Permeability;
    const double nu0 = sqrt(mu0 / eps0);

    dbl x = fv->x[0];
    dbl y = fv->x[1];
    dbl z = fv->x[2];
    dbl freq = upd->EM_Frequency;
    dbl lambda0 = c0 / freq;
    dbl k0 = 2 * M_PI / lambda0;

    dbl d = pml_outer_radius - pml_inner_radius;

    dbl L_mag = sqrt(x * x + y * y + z * z);
    dbl L = MAX(L_mag - pml_inner_radius, 0.0);
    if (L <= 1e-8) {
      permeability_matrix[0] = 1.0;
      permeability_matrix[1] = 1.0;
      permeability_matrix[2] = 1.0;
      return true;
    } else {
      dbl L_c[DIM] = {d, d, d};
      dbl pml_x = x * L / d;
      dbl pml_y = y * L / d;
      dbl pml_z = z * L / d;

      dbl amax[DIM];
      amax[0] = L_c[0] < 1e-8 ? 1 : 1 + amp * pow(fabs(pml_x / L_c[0]), power);
      amax[1] = L_c[1] < 1e-8 ? 1 : 1 + amp * pow(fabs(pml_y / L_c[1]), power);
      amax[2] = L_c[2] < 1e-8 ? 1 : 1 + amp * pow(fabs(pml_z / L_c[2]), power);

      dbl cond[DIM];
      cond[0] = L_c[0] < 1e-8 ? 0
                              : cond_max * sin(fabs(M_PI * pml_x / (2 * L_c[0]))) *
                                    sin(fabs(M_PI * pml_x / (2 * L_c[0])));
      cond[1] = L_c[1] < 1e-8 ? 0
                              : cond_max * sin(fabs(M_PI * pml_y / (2 * L_c[1]))) *
                                    sin(fabs(M_PI * pml_y / (2 * L_c[1])));
      cond[2] = L_c[2] < 1e-8 ? 0
                              : cond_max * sin(fabs(M_PI * pml_z / (2 * L_c[2]))) *
                                    sin(fabs(M_PI * pml_z / (2 * L_c[2])));

      complex double s[DIM] = {amax[0] * (1 + _Complex_I * nu0 * cond[0] / k0),
                               amax[1] * (1 + _Complex_I * nu0 * cond[1] / k0),
                               amax[2] * (1 + _Complex_I * nu0 * cond[2] / k0)};

      permeability_matrix[0] = s[1] * s[2] / s[0];
      permeability_matrix[1] = s[0] * s[2] / s[1];
      permeability_matrix[2] = s[0] * s[1] / s[2];
      return true;
    }
  } break;
  default:
    GOMA_EH(GOMA_ERROR, "Unknown permeability model");
  }

  *permeability_out = permeability;
  return false;
}

int incident_wave(
    dbl x, dbl y, dbl z, dbl omega, complex double wave[DIM], complex double curl_wave[DIM]) {

  switch (mp->IncidentWaveModel) {
  case CONSTANT: {
    GOMA_WH(GOMA_ERROR,
            "Incident wave being called but was not set or is set to CONSTANT, will not be used");
    for (int i = 0; i < DIM; i++) {
      wave[i] = 0;
      curl_wave[i] = 0;
    }
  } break;
  case EM_INC_PLANE_Z_WAVE: {
    for (int i = 0; i < DIM; i++) {
      wave[i] = 0;
      curl_wave[i] = 0;
    }

    dbl E0 = mp->u_incident_wave[0];

    wave[0] = E0 * cexp(-_Complex_I * omega * z);
    curl_wave[1] = -E0 * _Complex_I * omega * cexp(-_Complex_I * omega * z);
  } break;
  default:
    GOMA_EH(GOMA_ERROR, "Unknown Incident wave model %d", mp->IncidentWaveModel);
  }

  return 0;
}

/* assemble_ewave_nedelec
 *
 *  curl curl E - E = source
 */
int assemble_ewave_nedelec(dbl time) {
  int eqn_real = EM_E1_REAL;
  int eqn_imag = EM_E1_IMAG;

  const double c0 = 1.0 / (sqrt(upd->Free_Space_Permittivity * upd->Free_Space_Permeability));
  /*
   * Bail out fast if there's nothing to do...
   * But we might have the wrong eqn
   */
  if (!pd->e[pg->imtrx][eqn_real] || !pd->e[pg->imtrx][eqn_imag]) {
    return (-1);
  }

  dbl x = fv->x[0];
  dbl y = fv->x[1];
  dbl z = fv->x[2];
  dbl freq = upd->EM_Frequency;
  dbl lambda0 = c0 / freq;
  dbl k0 = 2 * M_PI / lambda0;

  dbl sigma = mp->electrical_conductivity;
  complex double wave[3];
  complex double curl_wave[3];
  incident_wave(x, y, z, k0, wave, curl_wave);

  complex double permeability_matrix[DIM]; // diagonal matrix if exists
  complex double permittivity_matrix[DIM]; // diagonal matrix if exists
  complex double permittivity;
  complex double permeability;
  bool permeability_is_matrix = relative_permeability_model(&permeability, permeability_matrix);
  bool permittivity_is_matrix = relative_permittivity_model(&permittivity, permittivity_matrix);

  int reqn = R_EM_E1_REAL;
  int peqn_real = upd->ep[pg->imtrx][reqn];
  int ieqn = R_EM_E1_IMAG;
  int peqn_imag = upd->ep[pg->imtrx][ieqn];
  if (af->Assemble_Residual) {
    for (int i = 0; i < ei[pg->imtrx]->dof[eqn_real]; i++) {
      complex double diffusion = 0.0;

      for (int q = 0; q < DIM; q++) {
        if (permeability_is_matrix) {
          diffusion += bf[ieqn]->curl_phi[i][q] * (1.0 / permeability_matrix[q]) *
                       (fv->curl_em_er[q] + fv->curl_em_ei[q] * _Complex_I);
        } else {
          diffusion += bf[ieqn]->curl_phi[i][q] * (1.0 / permeability) *
                       (fv->curl_em_er[q] + fv->curl_em_ei[q] * _Complex_I);
        }
      }

      complex double advection = 0;

      for (int q = 0; q < DIM; q++) {
        if (permittivity_is_matrix) {
          advection -= k0 * k0 * bf[ieqn]->phi_e[i][q] *
                       (permittivity_matrix[q] * (fv->em_er[q] + fv->em_ei[q] * _Complex_I));
        } else {
          advection -= k0 * k0 * bf[ieqn]->phi_e[i][q] *
                       (((permittivity - _Complex_I * sigma / k0) *
                         (fv->em_er[q] + fv->em_ei[q] * _Complex_I)));
        }
      }
      dbl source = 0;

      lec->R[LEC_R_INDEX(peqn_real, i)] +=
          creal(diffusion + advection + source) * bf[eqn_real]->detJ * fv->wt * fv->h3;
      lec->R[LEC_R_INDEX(peqn_imag, i)] +=
          cimag(diffusion + advection + source) * bf[eqn_imag]->detJ * fv->wt * fv->h3;
    }
  }

  if (af->Assemble_Jacobian) {
    for (int i = 0; i < ei[pg->imtrx]->dof[eqn_real]; i++) {
      int var = EM_E1_REAL;
      int pvar_real = upd->vp[pg->imtrx][var];
      for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

        complex double diffusion = 0.0;

        for (int q = 0; q < DIM; q++) {
          if (permeability_is_matrix) {
            diffusion += bf[ieqn]->curl_phi[i][q] * (1.0 / permeability_matrix[q]) *
                         (bf[var]->curl_phi[j][q]);
          } else {
            diffusion +=
                bf[ieqn]->curl_phi[i][q] * (1.0 / permeability) * (bf[var]->curl_phi[j][q]);
          }
        }

        complex double advection = 0;

        for (int q = 0; q < DIM; q++) {
          if (permittivity_is_matrix) {
            advection -=
                k0 * k0 * bf[ieqn]->phi_e[i][q] * (permittivity_matrix[q]) * (bf[var]->phi_e[j][q]);
          } else {
            advection -= k0 * k0 * bf[ieqn]->phi_e[i][q] *
                         (permittivity - _Complex_I * sigma / k0) * (bf[var]->phi_e[j][q]);
          }
        }

        lec->J[LEC_J_INDEX(peqn_real, pvar_real, i, j)] +=
            creal(diffusion + advection) * bf[eqn_real]->detJ * fv->wt * fv->h3;
        lec->J[LEC_J_INDEX(peqn_imag, pvar_real, i, j)] +=
            cimag(diffusion + advection) * bf[eqn_real]->detJ * fv->wt * fv->h3;
      }

      var = EM_E1_IMAG;
      int pvar_imag = upd->vp[pg->imtrx][var];
      for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

        complex double diffusion = 0.0;

        for (int q = 0; q < DIM; q++) {
          if (permeability_is_matrix) {
            diffusion += bf[ieqn]->curl_phi[i][q] * (1.0 / permeability_matrix[q]) *
                         (bf[var]->curl_phi[j][q] * _Complex_I);
          } else {
            diffusion += bf[ieqn]->curl_phi[i][q] * (1.0 / permeability) *
                         (bf[var]->curl_phi[j][q] * _Complex_I);
          }
        }

        complex double advection = 0;

        for (int q = 0; q < DIM; q++) {
          if (permittivity_is_matrix) {
            advection -= k0 * k0 * bf[ieqn]->phi_e[i][q] * (permittivity_matrix[q]) *
                         (bf[var]->phi_e[j][q] * _Complex_I);
          } else {
            advection -= k0 * k0 * bf[ieqn]->phi_e[i][q] *
                         (permittivity - _Complex_I * sigma / k0) *
                         (bf[var]->phi_e[j][q] * _Complex_I);
          }
        }

        lec->J[LEC_J_INDEX(peqn_real, pvar_imag, i, j)] +=
            creal(diffusion + advection) * bf[eqn_real]->detJ * fv->wt * fv->h3;
        lec->J[LEC_J_INDEX(peqn_imag, pvar_imag, i, j)] +=
            cimag(diffusion + advection) * bf[eqn_real]->detJ * fv->wt * fv->h3;
      }
    }
  }
  return (0);
} // end of assemble_ewave_curlcurl

int em_mms_force(dbl x, dbl y, dbl z, complex double force[DIM]) {

  //  force[0] =
  //      x * y * (1 - y*y) * (1 - z*z) + 2 * x * y * (1 - z*z);
  //  force[1] =    y*y * (1 - x*x) * (1 - z*z) + (1 - y*y) * (2 - x*x - z*z);
  //  force[2] =    y * z * (1 - x*x) * (1 - y*y) + 2 * y * z * (1 - x*x);
  // force[0] = 0;
  // force[1] = 0;
  // force[2] = -y;
  // force[0] = -1;
  // force[1] = 0;
  // force[2] = 0;
  // force[0] = y*y;
  // force[1] = x*y;
  // force[2] = 0;
  // force[0] = -cos(z);
  // force[1] = -cos(x);
  // force[2] = -cos(y);
  // force[0] = sin(y);
  // force[1] = sin(z);
  // force[2] = sin(x);
  x *= 0.005;
  y *= 0.005;
  z *= 0.005;
  complex double I = _Complex_I;
  force[0] = 1.000025 * sin(y);
  force[1] = 1.000025 * sin(z);
  force[2] = 1.000025 * sin(x);
  force[0] += I * 1.000025 * sin(y);
  force[1] += I * 1.000025 * sin(z);
  force[2] += I * 1.000025 * sin(x);
  // force[0] = (-2*y*y - 2*z*z + (1-y*y) * (1-z*z) + 4);
  // force[1] = (-2*x*x - 2*z*z + (1-x*x) * (1-z*z) + 4);
  // force[2] = (-2*x*x - 2*y*y + (1-x*x) * (1-y*y) + 4);
  // force[0] = -0.75 * cexp(-I * (x * sin(M_PI / 9.0) + y / 2.0));
  // force[1] = 0.75 * cexp(-I * (x * sin(M_PI / 9.0) + z / 2.0)) +
  //           cexp(-I * (x * sin(M_PI / 9.0) + z / 2.0)) * sin(M_PI / 9) * sin(M_PI / 9) -
  //           0.5 * cexp(-I * (x * sin(M_PI / 9.0) + y / 2.0)) * sin(M_PI / 9);
  // force[2] = -0.75 * cexp(-I * (x * sin(M_PI / 9.0) + y / 2.0)) +
  //           cexp(-I * (x * sin(M_PI / 9.0) + y / 2.0)) * sin(M_PI / 9);

  // force[0] = -0.75 * cexp(-I * (x + y / 2.0));
  // force[1] = x*x*cexp(-I*x*z/2)/4 + z*z*cexp(-I*z*z/2)/4 - cexp(-I*x*z/2) - cexp(-I*(x+y/2))/2;
  // force[2] = cexp(-I*(x+y/2))/4;

  // force[0] = 0;
  // force[1] = -sin(x)*cos(y);
  // force[2] = -I*sin(z)*cos(x);

  // dbl force[DIM] = {0, 1, 1};
  // dbl force[DIM] = {cos(x)*sin(y), sin(y)*cos(z), sin(x)*cos(z)};
  return 0;
}

int em_mms_exact(dbl x, dbl y, dbl z, complex double exact[DIM]) {
  // exact[0] = (1.0 - y*y)*(1.0 - z*z);
  // exact[1] = (1.0 - x*x)*(1.0 - z*z);
  // exact[2] = (1.0 - x*x)*(1.0 - y*y);
  // exact[0] = y*y;
  // exact[1] = x*y;
  // exact[2] = 0;
  x *= 0.005;
  y *= 0.005;
  z *= 0.005;
  complex double I = _Complex_I;
  exact[0] = sin(y);
  exact[1] = sin(z);
  exact[2] = sin(x);
  exact[0] += I * sin(y);
  exact[1] += I * sin(z);
  exact[2] += I * sin(x);
  // exact[0] = cexp(-I * (x * sin(M_PI / 9.0) + y / 2.0));
  // exact[1] = cexp(-I * (x * sin(M_PI / 9.0) + z / 2.0));
  // exact[2] = cexp(-I * (x * sin(M_PI / 9.0) + y / 2.0));
  // exact[0] = cos(x)*sin(y) + sin(x)*cos(z)*I;
  // exact[1] = 0;
  // exact[2] = 0;
  return 0;
}
