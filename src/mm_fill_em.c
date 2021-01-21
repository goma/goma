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

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

/* GOMA include files */
#define GOMA_MM_FILL_EM_C
#include "std.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_masks.h"
#include "rf_io_const.h"
#include "rf_io_structs.h"
#include "rf_io.h"
#include "rf_mp.h"
#include "el_elm.h"
#include "el_geom.h"
#include "rf_bc_const.h"
#include "rf_solver_const.h"
#include "rf_solver.h"
#include "rf_fill_const.h"
#include "rf_vars_const.h"
#include "mm_mp_const.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_eh.h"
#include "mm_std_models.h"
#include "mm_std_models_shell.h"


#include "mm_mp.h"
#include "mm_mp_structs.h"
#include "mm_fill_terms.h"


#include "goma.h"
#include "mm_species.h"
#include "rf_allo.h"
#include "mm_fill_em.h"
#include "mm_fill_util.h"


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

int
assemble_emwave(double time,	/* present time value */
                  double tt,	/* parameter to vary time integration from
                                 * explicit (tt = 1) to implicit (tt = 0) */
                  double dt,	/* current time step size */
                  const PG_DATA *pg_data,
                  const int em_eqn,	/* emwave eqn id and var id	*/
                  const int em_var,
                  const int em_conjvar ) {
  int eqn, var, peqn, pvar, dim, p, q, b, w, i, j, status;
  int dir=0;			/* identity of conjugate variable  */

  dbl EMF = 0, EMF_conj = 0;		/* acoustic pressure	*/
  dbl omega, emf_coeff=0, conj_coeff=0;
  dbl emf_coeff_dn=0, conj_coeff_dn=0;
  dbl emf_coeff_dk=0, conj_coeff_dk=0;
  //dbl mag_permeability=12.57e-07;  // H/m
  dbl mag_permeability=12.57e-07;  // H/m
  int cross_field_var;
  dbl cross_field[DIM];

  struct emwave_stabilization em_stab;
  em_stab.em_eqn = em_eqn;
  em_stab.em_var = em_var;
  em_stab.type = EM_STAB_DPHI_DIV; // enum supports phi_div, dphi_div,
                           // divphi_div, phi_divsquared and
                           // dphi_divsquared


  dbl n;				/* Refractive index. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_n_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_n = &d_n_struct;

  dbl k;				/* Acoustic wave number. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  dbl advection;			/* For terms and their derivatives */

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

  dbl h3;			/* Volume element (scale factors). */
  dbl dh3dmesh_bj;		/* Sensitivity to (b,j) mesh dof. */

  dbl det_J;

  dbl d_det_J_dmeshbj;			/* for specified (b,j) mesh dof */
  dbl dgrad_phi_i_dmesh[DIM];		/* ditto.  */
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
  dim   = pd->Num_Dim;
  eqn   = em_eqn;
  /*
   * Bail out fast if there's nothing to do...
   */
  if ( ! pd->e[eqn] )
    {
      return(status);
    }

  wt = fv->wt;				/* Gauss point weight. */
  h3 = fv->h3;			/* Differential volume element. */
  det_J = bf[eqn]->detJ;		/* Really, ought to be mesh eqn. */

  /*
   * Material property constants, etc. Any variations at this Gauss point
   * are calculated with user-defined subroutine options.
   *
   * Properties to consider here are R and k. R and k we will allow to vary
   * with temperature, spatial coordinates, and species concentration.
   */

  omega = upd->Acoustic_Frequency;
  n = refractive_index( d_n, time );

  k = extinction_index( d_k, time );

  // Compute complex impedance
  complex cpx_refractive_index, cpx_rel_permittivity,
      cpx_permittivity;//, impedance;

  cpx_refractive_index = n + _Complex_I*k; // k > 0 is extinction
  cpx_rel_permittivity = SQUARE(cpx_refractive_index);
  cpx_permittivity = cpx_rel_permittivity*mp->permittivity;

  //impedance = csqrt(mag_permeability/rel_permittivity);

  switch(em_var) {
    case EM_E1_REAL:
         EMF = fv->em_er[0];  EMF_conj = fv->em_ei[0]; dir=0;
         break;
    case EM_E1_IMAG:
         EMF = fv->em_ei[0];  EMF_conj = fv->em_er[0]; dir=0;
         break;
    case EM_E2_REAL:
         EMF = fv->em_er[1];  EMF_conj = fv->em_ei[1]; dir=1;
         break;
    case EM_E2_IMAG:
         EMF = fv->em_ei[1];  EMF_conj = fv->em_er[1]; dir=1;
         break;
    case EM_E3_REAL:
         EMF = fv->em_er[2];  EMF_conj = fv->em_ei[2]; dir=2;
         break;
    case EM_E3_IMAG:
         EMF = fv->em_ei[2];  EMF_conj = fv->em_er[2]; dir=2;
         break;
    case EM_H1_REAL:
         EMF = fv->em_hr[0];  EMF_conj = fv->em_hi[0]; dir=0;
         break;
    case EM_H1_IMAG:
         EMF = fv->em_hi[0];  EMF_conj = fv->em_hr[0]; dir=0;
         break;
    case EM_H2_REAL:
         EMF = fv->em_hr[1];  EMF_conj = fv->em_hi[1]; dir=1;
         break;
    case EM_H2_IMAG:
         EMF = fv->em_hi[1];  EMF_conj = fv->em_hr[1]; dir=1;
         break;
    case EM_H3_REAL:
         EMF = fv->em_hr[2];  EMF_conj = fv->em_hi[2]; dir=2;
         break;
    case EM_H3_IMAG:
         EMF = fv->em_hi[2];  EMF_conj = fv->em_hr[2]; dir=2;
         break;
    default:
         EH(-1,"Invalid EM variable name.\n");
         return -1;
  }
  switch(em_var) {
  case EM_E1_REAL:
  case EM_E2_REAL:
  case EM_E3_REAL:
         emf_coeff = -omega*cimag(cpx_permittivity);
         conj_coeff = -omega*creal(cpx_permittivity);
         //emf_coeff = 0;
         //conj_coeff = 0;
         emf_coeff_dn = -omega*cimag(cpx_permittivity)/n;
         emf_coeff_dk = -omega*cimag(cpx_permittivity)/k;
         conj_coeff_dn = -omega*creal(cpx_permittivity)/n;
         conj_coeff_dk = -omega*creal(cpx_permittivity)/k;
         cross_field_var = EM_H1_REAL;
         for ( p=0; p<VIM; p++) {
           cross_field[p] = fv->em_hr[p];
         }
         em_stab.stabilization_field_var = EM_E1_REAL;
         calc_emwave_stabilization_term(&em_stab, 1.0);
         break;
    case EM_E1_IMAG:
    case EM_E2_IMAG:
    case EM_E3_IMAG:
         emf_coeff = -omega*cimag(cpx_permittivity);
         conj_coeff = omega*creal(cpx_permittivity);
         //emf_coeff = 0;
         //conj_coeff = 0;
         emf_coeff_dn = -omega*cimag(cpx_permittivity)/n;
         emf_coeff_dk = -omega*cimag(cpx_permittivity)/k;
         conj_coeff_dn = omega*creal(cpx_permittivity)/n;
         conj_coeff_dk = omega*creal(cpx_permittivity)/k;
         cross_field_var = EM_H1_IMAG;
         for ( p=0; p<VIM; p++) {
           cross_field[p] = fv->em_hi[p];
         }
         em_stab.stabilization_field_var = EM_E1_IMAG;
         calc_emwave_stabilization_term(&em_stab, 1.0);
         break;
    case EM_H1_REAL:
    case EM_H2_REAL:
    case EM_H3_REAL:
         emf_coeff = 0;
         conj_coeff = -omega*mag_permeability;
         //conj_coeff = 0;
         emf_coeff_dn = 0; emf_coeff_dk = 0; conj_coeff_dn = 0; conj_coeff_dk = 0;
         cross_field_var = EM_E1_REAL;
         for ( p=0; p<VIM; p++) {
           cross_field[p] = fv->em_er[p];
         }
         em_stab.stabilization_field_var = EM_H1_REAL;
         calc_emwave_stabilization_term(&em_stab, 1.0);
         break;
    case EM_H1_IMAG:
    case EM_H2_IMAG:
    case EM_H3_IMAG:
         emf_coeff = 0;
         conj_coeff = omega*mag_permeability;
         //conj_coeff = 0;
         emf_coeff_dn = 0; emf_coeff_dk = 0; conj_coeff_dn = 0; conj_coeff_dk = 0;
         cross_field_var = EM_E1_IMAG;
         for ( p=0; p<VIM; p++) {
           cross_field[p] = fv->em_ei[p];
         }
         em_stab.stabilization_field_var = EM_H1_IMAG;
         calc_emwave_stabilization_term(&em_stab, 1.0);
         break;
    default:
      EH(-1, "assemble_emwave must be called with a usable em_var\n");
      return -1;
    }
  /*
   * Residuals___________________________________________________________
   */

  if ( af->Assemble_Residual )
    {
      eqn = em_eqn;
      peqn = upd->ep[eqn];
      var = em_var;
      for ( i=0; i<ei->dof[eqn]; i++)
        {

          /* this is an optimization for xfem */
          if ( xfem != NULL )
            {
              int xfem_active, extended_dof, base_interp, base_dof;
              xfem_dof_state( i, pd->i[eqn], ei->ielem_shape,
                              &xfem_active, &extended_dof, &base_interp, &base_dof );
              if ( extended_dof && !xfem_active ) continue;
            }
          phi_i = bf[eqn]->phi[i];

          advection = 0.;
          if ( pd->e[eqn] & T_ADVECTION )
            {
              advection += emf_coeff*EMF;
              advection += conj_coeff*EMF_conj;
              advection *= phi_i*h3;
              advection *= det_J*wt;
              advection *= pd->etm[eqn][(LOG2_ADVECTION)];
            }

          diffusion = 0.;
          if ( pd->e[eqn] & T_DIFFUSION )
            {
              for ( p=0; p<VIM; p++)
                {
                  grad_phi_i[p] = bf[eqn]->grad_phi[i] [p];
                }

              for ( p=0; p<VIM; p++) {
                for ( q=0; q<VIM; q++) {
                  diffusion -= permute(p,q,dir)*grad_phi_i[p]*cross_field[q];
                }
              }

              diffusion += em_stab.residual_term[i];

              diffusion *= det_J * wt;
              diffusion *= h3;
              diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
            }

          lec->R[LEC_R_INDEX(peqn,i)] += advection + diffusion;


        }
    }


  /*
   * Jacobian terms...
   */

  if ( af->Assemble_Jacobian ) {
    eqn   = em_eqn;
    peqn = upd->ep[eqn];
    for ( i=0; i<ei->dof[eqn]; i++) {

      // this is an optimization for xfem
      if ( xfem != NULL ) {
        int xfem_active, extended_dof, base_interp, base_dof;
        xfem_dof_state( i, pd->i[eqn], ei->ielem_shape,
                        &xfem_active, &extended_dof, &base_interp, &base_dof );
        if ( extended_dof && !xfem_active ) continue;
      }

      phi_i = bf[eqn]->phi[i];

      /*
       * Set up some preliminaries that are needed for the (a,i)
       * equation for bunches of (b,j) column variables...
       */

      for ( p=0; p<VIM; p++) {
        grad_phi_i[p] = bf[eqn]->grad_phi[i][p];
      }

      /*
       * EMF
       */
      var = em_var;
      if ( pd->v[var] ) {
        pvar = upd->vp[var];
        for ( j=0; j<ei->dof[var]; j++) {
          phi_j = bf[var]->phi[j];

          advection = 0.;
          if ( pd->e[eqn] & T_ADVECTION ) {
            advection += phi_i * emf_coeff*phi_j * det_J*wt;
            advection *= h3;
            advection *= pd->etm[eqn][(LOG2_ADVECTION)];
          }

          lec->J[LEC_J_INDEX(peqn,pvar,i,j)] += advection;
        }
      }
      /*
       *  EMF_conj
       */
      var = em_conjvar;
      if ( pd->v[var] ) {
        pvar = upd->vp[var];
        for ( j=0; j<ei->dof[var]; j++) {
          phi_j = bf[var]->phi[j];

          advection = 0.;
          if ( pd->e[eqn] & T_ADVECTION ) {
            advection += phi_i * conj_coeff*phi_j * det_J*wt;
            advection *= h3;
            advection *= pd->etm[eqn][(LOG2_ADVECTION)];
          }

          lec->J[LEC_J_INDEX(peqn,pvar,i,j)] += advection;
        }
      }
      /*
       *  cross_field
       */
      for ( b=0; b<dim; b++) {
        var = cross_field_var + b;
        if ( pd->v[var] ) {
          pvar = upd->vp[var];
          for ( j=0; j<ei->dof[var]; j++) {
            // for a cross product, this isn't quite right
            // but it will work as along as all scalar
            // components of the vector field have the
            // same basis functions
            phi_j = bf[var]->phi[j];

            diffusion = 0.;
            if ( pd->e[eqn] & T_DIFFUSION ) {
              for ( p=0; p<VIM; p++) {
                grad_phi_i[p] = bf[eqn]->grad_phi[i] [p];
              }

              for ( p=0; p<VIM; p++) {
                for ( q=0; q<VIM; q++) {
                  diffusion -= permute(p,q,dir)*grad_phi_i[p]*delta(q,b)*phi_j;
                }
              }

              //diffusion += stabilization_coefficient*phi_i
              //    *bf[var]->grad_phi[j][b];

              diffusion *= det_J * wt;
              diffusion *= h3;
              diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];

              lec->J[LEC_J_INDEX(peqn,pvar,i,j)] += diffusion;
            }
          }
        }
      }
      /*
       *  stabilization field
       */
      for ( b=0; b<dim; b++) {
        var = em_stab.stabilization_field_var + b;
        if ( pd->v[var] ) {
          pvar = upd->vp[var];
          for ( j=0; j<ei->dof[var]; j++) {
            diffusion = 0.;
            if ( pd->e[eqn] & T_DIFFUSION ) {
              diffusion += em_stab.jacobian_term[i][b][j];

              diffusion *= det_J * wt;
              diffusion *= h3;
              diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];

              lec->J[LEC_J_INDEX(peqn,pvar,i,j)] += diffusion;
            }
          }
        }
      }
          /*
           * J_e_T
           */
          var = TEMPERATURE;
          if ( pd->v[var] )
            {
              pvar = upd->vp[var];
              for ( j=0; j<ei->dof[var]; j++)
                {
                  phi_j = bf[var]->phi[j];

                  advection = 0.;
                  if ( pd->e[eqn] & T_ADVECTION )
                    {
                      advection += phi_i * (EMF*(emf_coeff_dn*d_n->T[j]+emf_coeff_dk*d_k->T[j])
                               +EMF_conj*(conj_coeff_dn*d_n->T[j]+conj_coeff_dk*d_k->T[j])) * det_J*wt;
                      advection *= h3;
                      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
                    }

                  lec->J[LEC_J_INDEX(peqn,pvar,i,j)] += advection;
                }
            }

          /*
           * J_e_d
           */
          for ( b=0; b<dim; b++)
            {
              var = MESH_DISPLACEMENT1+b;
              if ( pd->v[var] )
                {
                  pvar = upd->vp[var];
                  for ( j=0; j<ei->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];

                      dh3dmesh_bj = fv->dh3dq[b] * bf[var]->phi[j];

                      d_det_J_dmeshbj = bf[eqn]->d_det_J_dm[b][j];

                      advection = 0.;
                  if ( pd->e[eqn] & T_ADVECTION )
                    {
                      advection += phi_i * (EMF*(emf_coeff_dn*d_n->X[b][j]+emf_coeff_dk*d_k->X[b][j])
                          +EMF_conj*(conj_coeff_dn*d_n->X[b][j]+conj_coeff_dk*d_k->X[b][j])) * det_J*h3*wt;
                      advection += phi_i*(emf_coeff*EMF+conj_coeff*EMF_conj) *
                                        (d_det_J_dmeshbj*h3 + det_J*dh3dmesh_bj)*wt;
                      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
                    }

                          /*
                           * multiple parts:
                           * 	diff_a = Int(...d(grad_phi_i)/dmesh.q h3 |Jv|)
                           *	diff_b = Int(...grad_phi_i.d(q)/dmesh h3 |Jv|)
                           *	diff_c = Int(...grad_phi_i.q h3 d(|Jv|)/dmesh)
                           *	diff_d = Int(...grad_phi_i.q dh3/dmesh |Jv|  )
                           */
                      diffusion = 0.;
                      if ( pd->e[eqn] & T_DIFFUSION )
                        {
                          diff_a = 0.;
                          for ( p=0; p<dim; p++)
                            {
                              dgrad_phi_i_dmesh[p]
                                = bf[eqn]->d_grad_phi_dmesh[i][p] [b][j];

                              diff_a += dgrad_phi_i_dmesh[p] ;
                            }
                          diff_a *= det_J * h3 * wt;

                          diff_b = 0.;
                          for ( p=0; p<VIM; p++)
                            {
                              diff_b +=  grad_phi_i[p];
                            }
                          diff_b *= det_J * h3 * wt;

                          diff_c = 0.;
                          for ( p=0; p<dim; p++)
                            {
                              diff_c += grad_phi_i[p] ;
                            }
                          diff_c *= d_det_J_dmeshbj * h3 * wt;

                          diff_d = 0.;
                          for ( p=0; p<dim; p++)
                            {
                              diff_d += grad_phi_i[p] ;
                            }
                          diff_d *= det_J * dh3dmesh_bj * wt;

                          diffusion = diff_a + diff_b + diff_c + diff_d;

                          diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
                        }


                      lec->J[LEC_J_INDEX(peqn,pvar,i,j)] += advection + diffusion;
                    }
                }
            }

          /*
           * J_e_c
           */
          var = MASS_FRACTION;
          if ( pd->e[eqn] && pd->v[var] )
            {
              for ( w=0; w<pd->Num_Species_Eqn; w++)
                {
                  for ( j=0; j<ei->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];

                      advection = 0.;
                  if ( pd->e[eqn] & T_ADVECTION )
                    {
                      advection += phi_i * (EMF*(emf_coeff_dn*d_n->C[w][j]+emf_coeff_dk*d_k->C[w][j])
                               +EMF_conj*(conj_coeff_dn*d_n->C[w][j]+conj_coeff_dk*d_k->C[w][j])) * det_J*wt;
                      advection *= h3;
                      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
                    }

                      lec->J[LEC_J_INDEX(peqn, MAX_PROB_VAR + w, i, j)] += advection;
                    }
                }
            }

        }
    }

  return(status);
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

int assemble_ewave(double time, // present time
                   double tt,   // time integration method parameter
                   double dt,   // current time step size
                   const int em_eqn, // eqn id
                   const int em_var //  variable id - should match me_eqn
                   ) {
  int eqn, peqn, i;
  dbl mag_permeability=mp->magnetic_permeability;
  double omega, re_coeff, im_coeff;
  dbl n;				/* Refractive index. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_n_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_n = &d_n_struct;

  dbl k;				/* Acoustic wave number. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  int field_index;
  double curl_field[DIM];

  dbl advection;			/* For terms and their derivatives */

  dbl diffusion;

  dbl phi_i;
  dbl grad_phi_i[DIM];

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl h3;			/* Volume element (scale factors). */
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
  for (int p=0; p<DIM; p++) {
    curl_field [p] = 0.0;
  }

  eqn   = em_eqn;
  /*
   * Bail out fast if there's nothing to do...
   */
  if ( ! pd->e[eqn] )
    {
      return(1);
    }

  wt = fv->wt;				/* Gauss point weight. */
  h3 = fv->h3;			/* Differential volume element. */
  det_J = bf[eqn]->detJ;		/* Really, ought to be mesh eqn. */

  /*
   * Material property constants, etc. Any variations at this Gauss point
   * are calculated with user-defined subroutine options.
   *
   * Properties to consider here are R and k. R and k we will allow to vary
   * with temperature, spatial coordinates, and species concentration.
   */

  omega = upd->Acoustic_Frequency;
  n = refractive_index( d_n, time );

  k = extinction_index( d_k, time );

  // Compute complex material properties
  complex double cpx_refractive_index, cpx_rel_permittivity,
      cpx_permittivity;//, impedance;
  double r_elperm, i_elperm;

  cpx_refractive_index = n + _Complex_I*k; // k > 0 is extinction
  cpx_rel_permittivity = SQUARE(cpx_refractive_index);
  cpx_permittivity = cpx_rel_permittivity*mp->permittivity;

  // assumed to be constant in an element block
  r_elperm = creal(cpx_permittivity);
  i_elperm = cimag(cpx_permittivity);



  switch(em_eqn) {
  case R_EM_E1_REAL:
  case R_EM_E2_REAL:
  case R_EM_E3_REAL:
    field_index = em_eqn - R_EM_E1_REAL;
    re_coeff = -omega*omega*mag_permeability*r_elperm;
    im_coeff =  omega*omega*mag_permeability*i_elperm;
    for (int p=0; p<DIM; p++) {
      for (int q=0; q<DIM; q++) {
        for (int r=0; r<DIM; r++) {
          curl_field[p] += permute(p,q,r)*fv->grad_em_er[q][r];
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
    re_coeff = -omega*omega*mag_permeability*i_elperm;
    im_coeff = -omega*omega*mag_permeability*r_elperm;
    for (int p=0; p<DIM; p++) {
      for (int q=0; q<DIM; q++) {
        for (int r=0; r<DIM; r++) {

          curl_field[p] += permute(p,q,r)*fv->grad_em_ei[q][r];
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
    EH(-1, "assemble_ewave must be called with only electric field for em_eqn!");
    break;
  }
  if ( af->Assemble_Residual ) {
    eqn = em_eqn;
    peqn = upd->ep[eqn];
    for ( i=0; i<ei->dof[eqn]; i++) {
      phi_i = bf[eqn]->phi[i];

      advection = 0.;
      if ( pd->e[eqn] & T_ADVECTION )
        {
          advection += re_coeff*fv->em_er[field_index];
          advection += im_coeff+fv->em_ei[field_index];
          advection *= phi_i*h3;
          advection *= det_J*wt;
          advection *= pd->etm[eqn][(LOG2_ADVECTION)];
        }
      diffusion = 0.;
      if ( pd->e[eqn] & T_DIFFUSION )
        {
          for (int p=0; p<VIM; p++)
            {
              grad_phi_i[p] = bf[eqn]->grad_phi[i] [p];
            }

          for (int p=0; p<VIM; p++) {
            for (int  q=0; q<VIM; q++) {
              diffusion -= delta(field_index,p)*permute(field_index,p,q)*grad_phi_i[p]*curl_field[q];
            }
          }

          diffusion += em_stab.residual_term[i];

          diffusion *= det_J * wt;
          diffusion *= h3;
          diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
        }

      lec->R[LEC_R_INDEX(peqn,i)] += advection + diffusion;
    }
  }
  return(0);
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

int assemble_ewave_curlcurl(double time, // present time
                   double tt,   // time integration method parameter
                   double dt,   // current time step size
                   const int em_eqn, // eqn id
                   const int em_var //  variable id - should match me_eqn
                   ) {
  int eqn;
  dbl mag_permeability=mp->magnetic_permeability;
  double omega, re_coeff, im_coeff;
  /*
  struct emwave_stabilization em_stab;
  em_stab.em_eqn = em_eqn;
  em_stab.em_var = em_var;
  em_stab.type = EM_STAB_DPHI_DIV; // enum supports phi_div, dphi_div,
                           // divphi_div, phi_divsquared and
                           // dphi_divsquared
  */
  eqn   = em_eqn;
  /*
   * Bail out fast if there's nothing to do...
   * But we might have the wrong eqn
   */
  for (int a=0; a<DIM; a++){
    if (pd->e[R_EM_E1_REAL+a]) {
      eqn = R_EM_E1_REAL+a;
      break;
    }
  }
  if ( ! pd->e[eqn] )
    {
      return(1);
    }

  omega = upd->Acoustic_Frequency;
  dbl n;				/* Refractive index. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_n_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_n = &d_n_struct;

  dbl k;				/* Extinction coefficient */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  n = refractive_index( d_n, time );
  k = extinction_index( d_k, time );

  // Compute complex material properties
  complex double cpx_refractive_index, cpx_rel_permittivity,
      cpx_permittivity;//, impedance;
  double r_elperm, i_elperm;

  cpx_refractive_index = n + _Complex_I*k; // k > 0 is extinction
  cpx_rel_permittivity = SQUARE(cpx_refractive_index);
  cpx_permittivity = cpx_rel_permittivity*mp->permittivity;

  // assumed to be constant in an element block
  r_elperm = creal(cpx_permittivity);
  i_elperm = cimag(cpx_permittivity);
  re_coeff =  omega*omega*mag_permeability*r_elperm;
  im_coeff =  omega*omega*mag_permeability*i_elperm;

  int reqn, ieqn;
  double radvection_etm, rdiffusion_etm, rsource_etm;
  double iadvection_etm, idiffusion_etm, isource_etm;
  //double rmass_etm;
  //double imass_etm;
  double stab_scale = 4.0;
  if ( af->Assemble_Residual ) {
    for (int i = 0; i < ei->dof[eqn]; i++) {
      for (int a = 0; a < DIM; a++) {
        reqn = R_EM_E1_REAL + a;
        ieqn = R_EM_E1_IMAG + a;
        //rmass_etm = pd->etm[reqn][(LOG2_MASS)];
        radvection_etm = pd->etm[reqn][(LOG2_ADVECTION)];
        rdiffusion_etm = pd->etm[reqn][(LOG2_DIFFUSION)];
        rsource_etm = pd->etm[reqn][(LOG2_SOURCE)];
        //imass_etm = pd->etm[ieqn][(LOG2_MASS)];
        iadvection_etm = pd->etm[ieqn][(LOG2_ADVECTION)];
        idiffusion_etm = pd->etm[ieqn][(LOG2_DIFFUSION)];
        isource_etm = pd->etm[ieqn][(LOG2_SOURCE)];

        if ( pd->e[R_EM_E1_REAL+a] ) {
          int eqn_real = EM_E1_REAL+a;
          int eqn_imag = EM_E1_IMAG+a;
          int peqn_real = upd->ep[eqn_real];
          int peqn_imag = upd->ep[eqn_imag];

          double diffusion_real = 0.0;
          double diffusion_imag = 0.0;
          for (int q = 0; q < DIM; q++) {
            diffusion_real += bf[eqn_real]->curl_phi_e[i][a][q] * fv->curl_em_er[q];
            diffusion_imag += bf[eqn_imag]->curl_phi_e[i][a][q] * fv->curl_em_ei[q];
            //diffusion_real += bf->[eqn_real]->phi_e[i][a][q] * fv->em_er[q];
            //diffusion_imag += bf->[eqn_imag]->phi_e[i][a][q] * fv->em_ei[q];

            // R = curl(E)
            //diffusion_real += delta(a,q)*bf[eqn_real]->phi[i]*fv->curl_em_er[a];
            //diffusion_imag += delta(a,q)*bf[eqn_imag]->phi[i]*fv->curl_em_ei[a];
          }
          // R = E
          //diffusion_real = bf[eqn_real]->phi[i]*fv->em_er[a];
          //diffusion_imag = bf[eqn_imag]->phi[i]*fv->em_ei[a];
          double stab_real = 0.0;
          double stab_imag = 0.0;
          double lagr_term_real = 0.0;
          double lagr_term_imag = 0.0;

          if (pd->v[R_EM_CONT_REAL] && pd->v[R_EM_CONT_IMAG]) {
            lagr_term_real = bf[eqn_real]->grad_phi[i][a] * (r_elperm * fv->epr - i_elperm * fv->epi);
            lagr_term_imag = bf[eqn_real]->grad_phi[i][a] * (r_elperm * fv->epi + i_elperm * fv->epr);
            lagr_term_real = bf[eqn_real]->grad_phi[i][a] * (fv->epr);
            lagr_term_imag = bf[eqn_real]->grad_phi[i][a] * (fv->epi);
          } else {
            stab_real = stab_scale * bf[eqn_real]->grad_phi[i][a] * (fv->grad_em_er[0][0] + fv->grad_em_er[1][1] + fv->grad_em_er[2][2]);
            stab_imag = stab_scale * bf[eqn_real]->grad_phi[i][a] * (fv->grad_em_ei[0][0] + fv->grad_em_ei[1][1] + fv->grad_em_ei[2][2]);
          }
          //double advection_real = 0.0;
          //double advection_imag = 0.0;
          double advection_real = -bf[eqn_real]->phi[i] * (re_coeff * fv->em_er[a] - im_coeff * fv->em_ei[a]);
          double advection_imag = -bf[eqn_imag]->phi[i] * (re_coeff * fv->em_ei[a] + im_coeff * fv->em_er[a]);
          double src_real = 0.0;
          double src_imag = 0.0;
          switch(a) {
          case 0:
            //src_term = bf[eqn_real]->phi[i] * (-fv->x[1]*fv->x[1]*r_elperm*r_elperm - 3);
                //- stab_scale * bf[eqn_real]->grad_phi[i][a] * ( -fv->x[0]);

            // src term for R = E
            //src_term = -bf[eqn_real]->phi[i] * (fv->x[1]*fv->x[1]);

            // src term for R = curl(E)

            break;
          case 1:
            //src_term = bf[eqn_real]->phi[i] * (fv->x[0]*fv->x[1]*r_elperm*r_elperm);
                //- stab_scale * bf[eqn_real]->grad_phi[i][a] * ( -fv->x[0]);

            // src term for R = E
            //src_term = bf[eqn_real]->phi[i] * (fv->x[0]*fv->x[1]);

            // src term for R = curl(E)

            break;
          case 2:
            //src_term = - stab_scale * bf[eqn_real]->grad_phi[i][a] * ( -fv->x[0]);

            // src term for R = curl(E)
            //src_term = bf[eqn_real]->phi[i]*3.0*fv->x[1];

            break;
          }
          lec->R[LEC_R_INDEX(peqn_real,i)] += (advection_real*radvection_etm
                                   + diffusion_real*rdiffusion_etm
                                   + stab_real
                                   + src_real*rsource_etm
                                   + lagr_term_real)
              * bf[eqn_real]->detJ * fv->wt * fv->h3;

          lec->R[LEC_R_INDEX(peqn_imag,i)] += (advection_imag*iadvection_etm
                                   + diffusion_imag*idiffusion_etm
                                   + stab_imag
                                   + src_imag*isource_etm
                                   + lagr_term_imag)
              * bf[eqn_imag]->detJ * fv->wt * fv->h3;
        }
      }
    }
  }
  if ( af->Assemble_Jacobian ) {
    for (int i = 0; i < ei->dof[eqn]; i++) {
      for (int a = 0; a < DIM; a++) {
        reqn = R_EM_E1_REAL + a;
        ieqn = R_EM_E1_IMAG + a;
        //rmass_etm = pd->etm[reqn][(LOG2_MASS)];
        radvection_etm = pd->etm[reqn][(LOG2_ADVECTION)];
        rdiffusion_etm = pd->etm[reqn][(LOG2_DIFFUSION)];
        rsource_etm = pd->etm[reqn][(LOG2_SOURCE)];
        //imass_etm = pd->etm[ieqn][(LOG2_MASS)];
        iadvection_etm = pd->etm[ieqn][(LOG2_ADVECTION)];
        idiffusion_etm = pd->etm[ieqn][(LOG2_DIFFUSION)];
        isource_etm = pd->etm[ieqn][(LOG2_SOURCE)];
        if ( pd->e[R_EM_E1_REAL+a] ) {
          int eqn_real = EM_E1_REAL+a;
          int eqn_imag = EM_E1_IMAG+a;
          int peqn_real = upd->ep[eqn_real];
          int peqn_imag = upd->ep[eqn_imag];

          // Sensitivity to real parts of electric field
          for (int b = 0; b < DIM; b++) {
            if ( pd->e[R_EM_E1_REAL+b] ) {
              int var = EM_E1_REAL + b;
              int pvar_real = upd->vp[var];
              for (int j = 0; j < ei->dof[var]; j++) {
                double diffusion_real = 0;
                double diffusion_imag = 0;
                for (int q = 0; q < DIM; q++) {
                  diffusion_real += bf[eqn_real]->curl_phi_e[i][a][q] * bf[eqn_real]->curl_phi_e[j][b][q];
                }
                // R = curl(E)
                //diffusion_real += bf[eqn_real]->phi[i]*bf[eqn_real]->curl_phi_e[j][b][a];

                // Jacobian contribution for R = E
                //diffusion_real += delta(a,b)*bf[eqn_real]->phi[i]*bf[eqn_real]->phi[j];

                double stab_real = 0;
                double stab_imag = 0;
                if (!(pd->v[R_EM_CONT_REAL] && pd->v[R_EM_CONT_IMAG])) {
                  stab_real = stab_scale * bf[eqn_real]->grad_phi[i][a] * (bf[var]->grad_phi_e[j][b][0][0] +
                      bf[var]->grad_phi_e[j][b][1][1]  + bf[var]->grad_phi_e[j][b][2][2]);
                }

                //double advection_real = 0;
                //double advection_imag = 0;
                double advection_real = -bf[eqn_real]->phi[i] * (delta(a,b) * re_coeff * bf[var]->phi[j]);
                double advection_imag = -bf[eqn_imag]->phi[i] * (delta(a,b) * im_coeff * bf[var]->phi[j]);

                lec->J[LEC_J_INDEX(peqn_real,pvar_real,i,j)] += (diffusion_real*rdiffusion_etm
                                                       + advection_real*radvection_etm
                                                       + stab_real)
                    * bf[eqn_real]->detJ * fv->wt * fv->h3;

                lec->J[LEC_J_INDEX(peqn_imag,pvar_real,i,j)] += (diffusion_imag*idiffusion_etm
                                                       + advection_imag*iadvection_etm
                                                       + stab_imag)
                    * bf[eqn_imag]->detJ * fv->wt * fv->h3;
              }
            }
          }

          // Sensitivity to imaginary parts of electric field
          for (int b = 0; b < DIM; b++) {
            if ( pd->e[R_EM_E1_REAL+b] ) {
              int var = EM_E1_IMAG + b;
              int pvar_imag = upd->vp[var];
              for (int j = 0; j < ei->dof[var]; j++) {
                double diffusion_real = 0;
                double diffusion_imag = 0;
                for (int q = 0; q < DIM; q++) {
                  diffusion_imag += bf[eqn_imag]->curl_phi_e[i][a][q] * bf[var]->curl_phi_e[j][b][q];

                }

                // R = curl(E)
                //diffusion_imag += bf[eqn_imag]->phi[i]*bf[eqn_imag]->curl_phi_e[j][b][a];

                // Jacobian Contribution for R = E
                //diffusion_imag += delta(a,b)*bf[eqn_imag]->phi[i]*bf[eqn_imag]->phi[j];

                // Jacobian Contribution for R = curl(E)
                //diffusion_real += bf[eqn_real]->phi[i]*bf[eqn_real]->curl_phi_e[j][b][a];

                double stab_real = 0;
                double stab_imag = 0;
                if (!(pd->v[R_EM_CONT_REAL] && pd->v[R_EM_CONT_IMAG])) {
                  stab_imag = stab_scale * bf[eqn_real]->grad_phi[i][a] * (bf[var]->grad_phi_e[j][b][0][0] +
                      bf[var]->grad_phi_e[j][b][1][1]  + bf[var]->grad_phi_e[j][b][2][2]);
                }

                //double advection_real = 0;
                //double advection_imag = 0;
                double advection_real = bf[eqn_real]->phi[i] * (delta(a,b) * im_coeff * bf[var]->phi[j]);
                double advection_imag = -bf[eqn_imag]->phi[i] * (delta(a,b) * re_coeff * bf[var]->phi[j]);
                lec->J[LEC_J_INDEX(peqn_real,pvar_imag,i,j)] += (diffusion_real*rdiffusion_etm
                                                       + advection_real*radvection_etm
                                                       + stab_real)
                    * bf[eqn_real]->detJ * fv->wt * fv->h3;

                lec->J[LEC_J_INDEX(peqn_imag,pvar_imag,i,j)] += (diffusion_imag*idiffusion_etm
                                                       + advection_imag*iadvection_etm
                                                       + stab_imag)
                    * bf[eqn_imag]->detJ * fv->wt * fv->h3;
              }
            }
          }

          {
            int var = EM_CONT_REAL;
            if (pd->e[var]) {
              int pvar = upd->vp[var];
              for (int j = 0; j < ei->dof[var]; j++) {
                //dbl lagr_term_real = bf[eqn]->grad_phi[i][a] * r_elperm * bf[var]->phi[j];
                //dbl lagr_term_imag = bf[eqn]->grad_phi[i][a] * i_elperm * bf[var]->phi[j];
                dbl lagr_term_real = bf[eqn]->grad_phi[i][a] * bf[var]->phi[j];
                dbl lagr_term_imag = 0;
                lec->J[LEC_J_INDEX(peqn_real,pvar,i,j)] += (lagr_term_real) * bf[eqn_real]->detJ * fv->wt * fv->h3;
                lec->J[LEC_J_INDEX(peqn_imag,pvar,i,j)] += (lagr_term_imag) * bf[eqn_real]->detJ * fv->wt * fv->h3;
              }
            }
          }

          {
            int var = EM_CONT_IMAG;
            if (pd->e[var]) {
              int pvar = upd->vp[var];
              for (int j = 0; j < ei->dof[var]; j++) {
                //dbl lagr_term_real = bf[eqn]->grad_phi[i][a] * i_elperm * bf[var]->phi[j];
                //dbl lagr_term_imag = bf[eqn]->grad_phi[i][a] * r_elperm * bf[var]->phi[j];
                dbl lagr_term_real = 0;
                dbl lagr_term_imag = bf[eqn]->grad_phi[i][a] * bf[var]->phi[j];
                lec->J[LEC_J_INDEX(peqn_real,pvar,i,j)] += (lagr_term_real) * bf[eqn_real]->detJ * fv->wt * fv->h3;
                lec->J[LEC_J_INDEX(peqn_imag,pvar,i,j)] += (lagr_term_imag) * bf[eqn_real]->detJ * fv->wt * fv->h3;
              }
            }
          }
        }
      }
    }
  }
  return(0);
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

int assemble_ewave_laplacian(double time, // present time
                   double tt,   // time integration method parameter
                   double dt,   // current time step size
                   const int em_eqn, // eqn id
                   const int em_var //  variable id - should match me_eqn
                   ) {
  int eqn;
  dbl mag_permeability=mp->magnetic_permeability;
  double omega, re_coeff, im_coeff;
  /*
  struct emwave_stabilization em_stab;
  em_stab.em_eqn = em_eqn;
  em_stab.em_var = em_var;
  em_stab.type = EM_STAB_NONE; // enum supports phi_div, dphi_div,
                           // divphi_div, phi_divsquared and
                           // dphi_divsquared
  */
  eqn   = em_eqn;
  /*
   * Bail out fast if there's nothing to do...
   * But we might have the wrong eqn
   */
  for (int a=0; a<DIM; a++){
    if (pd->e[R_EM_E1_REAL+a]) {
      eqn = R_EM_E1_REAL+a;
      break;
    }
  }
  if ( ! pd->e[eqn] )
    {
      return(1);
    }

  omega = upd->Acoustic_Frequency;
  dbl n;				/* Refractive index. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_n_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_n = &d_n_struct;

  dbl k;				/* Acoustic wave number. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  n = refractive_index( d_n, time );
  k = extinction_index( d_k, time );

  // Compute complex material properties
  complex double cpx_refractive_index, cpx_rel_permittivity,
      cpx_permittivity;//, impedance;
  double r_elperm, i_elperm;

  cpx_refractive_index = n + _Complex_I*k; // k > 0 is extinction
  cpx_rel_permittivity = SQUARE(cpx_refractive_index);
  cpx_permittivity = cpx_rel_permittivity*mp->permittivity;

  // assumed to be constant in an element block
  r_elperm = creal(cpx_permittivity);
  i_elperm = cimag(cpx_permittivity);
  re_coeff =  omega*omega*mag_permeability*r_elperm;
  im_coeff =  omega*omega*mag_permeability*i_elperm;

  double stab_scale = 0.0;
  if ( af->Assemble_Residual ) {
    for (int i = 0; i < ei->dof[eqn]; i++) {
      for (int a = 0; a < DIM; a++) {
        if ( pd->e[R_EM_E1_REAL+a] ) {
          int eqn_real = EM_E1_REAL+a;
          int eqn_imag = EM_E1_IMAG+a;
          int peqn_real = upd->ep[eqn_real];
          int peqn_imag = upd->ep[eqn_imag];

          double diffusion_real = 0.0;
          double diffusion_imag = 0.0;

          for (int q = 0; q < DIM; q++) {
            // R = - (gradphi)_a dot (gradE)_a
            diffusion_real += -bf[eqn_real]->grad_phi[i][q]*fv->grad_em_er[q][a];
            diffusion_imag += -bf[eqn_imag]->grad_phi[i][q]*fv->grad_em_ei[q][a];
          }
          //diffusion_real = -fv->grad_em_er[a][0];
          // R = E
          //diffusion_real = bf[eqn_real]->phi[i]*fv->em_er[a];
          //diffusion_imag = bf[eqn_imag]->phi[i]*fv->em_ei[a];
          double stab_real = 0.0;
          double stab_imag = 0.0;
          double lagr_term_real = 0.0;
          double lagr_term_imag = 0.0;

          if (pd->v[R_EM_CONT_REAL] && pd->v[R_EM_CONT_IMAG]) {
            //lagr_term_real = bf[eqn_real]->grad_phi[i][a] * (r_elperm * fv->epr + i_elperm * fv->epi);
            //lagr_term_imag = bf[eqn_real]->grad_phi[i][a] * (r_elperm * fv->epi + i_elperm * fv->epr);
            //lagr_term_real = bf[eqn_real]->grad_phi[i][a] * (fv->epr);
            //lagr_term_imag = bf[eqn_real]->grad_phi[i][a] * (fv->epi);
          } else {
            stab_real = stab_scale * bf[eqn_real]->grad_phi[i][a] * (fv->grad_em_er[0][0] + fv->grad_em_er[1][1] + fv->grad_em_er[2][2]);
            stab_imag = stab_scale * bf[eqn_real]->grad_phi[i][a] * (fv->grad_em_ei[0][0] + fv->grad_em_ei[1][1] + fv->grad_em_ei[2][2]);
          }

          double advection_real = 0.0;
          double advection_imag = 0.0;
          advection_real = bf[eqn_real]->phi[i] * (re_coeff * fv->em_er[a] - im_coeff * fv->em_ei[a]);
          advection_imag = bf[eqn_imag]->phi[i] * (re_coeff * fv->em_ei[a] + im_coeff * fv->em_er[a]);

          double src_term = 0.0;
          switch(a) {
          case 0:
            //src_term = bf[eqn_real]->phi[i] * (-fv->x[1]*fv->x[1]*r_elperm*r_elperm - 3);
                //- stab_scale * bf[eqn_real]->grad_phi[i][a] * ( -fv->x[0]);

            // src term for R = E
            //src_term = -bf[eqn_real]->phi[i] * (fv->x[1]*fv->x[1]);

            // src term for R = curl(E)

            break;
          case 1:
            //src_term = bf[eqn_real]->phi[i] * (fv->x[0]*fv->x[1]*r_elperm*r_elperm);
                //- stab_scale * bf[eqn_real]->grad_phi[i][a] * ( -fv->x[0]);

            // src term for R = E
            //src_term = bf[eqn_real]->phi[i] * (fv->x[0]*fv->x[1]);

            // src term for R = curl(E)

            break;
          case 2:
            //src_term = - stab_scale * bf[eqn_real]->grad_phi[i][a] * ( -fv->x[0]);

            // src term for R = curl(E)
            //src_term = bf[eqn_real]->phi[i]*3.0*fv->x[1];

            break;
          }
          lec->R[LEC_R_INDEX(peqn_real,i)] += (advection_real + diffusion_real + stab_real + src_term + lagr_term_real) * bf[eqn_real]->detJ * fv->wt * fv->h3;
          lec->R[LEC_R_INDEX(peqn_imag,i)] += (advection_imag + diffusion_imag + stab_imag + lagr_term_imag) * bf[eqn_imag]->detJ * fv->wt * fv->h3;
        }
      }
    }
  }
  if ( af->Assemble_Jacobian ) {
    for (int i = 0; i < ei->dof[eqn]; i++) {
      for (int a = 0; a < DIM; a++) {
        if ( pd->e[R_EM_E1_REAL+a] ) {
          int eqn_real = EM_E1_REAL+a;
          int eqn_imag = EM_E1_IMAG+a;
          int peqn_real = upd->ep[eqn_real];
          int peqn_imag = upd->ep[eqn_imag];

          // Sensitivity to real parts of electric field
          for (int b = 0; b < DIM; b++) {
            if ( pd->e[R_EM_E1_REAL+b] ) {
              int var = EM_E1_REAL + b;
              int pvar_real = upd->vp[var];
              for (int j = 0; j < ei->dof[var]; j++) {
                double diffusion_real = 0;
                double diffusion_imag = 0;
                for (int q = 0; q < DIM; q++) {
                  //diffusion_real += bf[eqn_real]->curl_phi_e[i][a][q] * bf[eqn_real]->curl_phi_e[j][b][q];
                  diffusion_real += -delta(a,b)*bf[eqn_real]->grad_phi[i][q]*bf[eqn_real]->grad_phi[j][q];
                }

                //diffusion_real += -bf[eqn_real]->grad_phi[j][b];

                // R = curl(E)
                //diffusion_real += bf[eqn_real]->phi[i]*bf[eqn_real]->curl_phi_e[j][b][a];

                // Jacobian contribution for R = E
                //diffusion_real += delta(a,b)*bf[eqn_real]->phi[i]*bf[eqn_real]->phi[j];

                double stab_real = 0;
                double stab_imag = 0;
                if (!(pd->v[R_EM_CONT_REAL] && pd->v[R_EM_CONT_IMAG])) {
                  stab_real = stab_scale * bf[eqn_real]->grad_phi[i][a] * (bf[var]->grad_phi_e[j][b][0][0] +
                      bf[var]->grad_phi_e[j][b][1][1]  + bf[var]->grad_phi_e[j][b][2][2]);
                }

                double advection_real = 0;
                double advection_imag = 0;

                advection_real = bf[eqn_real]->phi[i] * (delta(a,b) * re_coeff * bf[var]->phi[j]);
                advection_imag = bf[eqn_imag]->phi[i] * (delta(a,b) * im_coeff * bf[var]->phi[j]);

                lec->J[LEC_J_INDEX(peqn_real,pvar_real,i,j)] += (diffusion_real + advection_real + stab_real) * bf[eqn_real]->detJ * fv->wt * fv->h3;
                lec->J[LEC_J_INDEX(peqn_imag,pvar_real,i,j)] += (diffusion_imag + advection_imag + stab_imag) * bf[eqn_imag]->detJ * fv->wt * fv->h3;
              }
            }
          }

          // Sensitivity to imaginary parts of electric field
          for (int b = 0; b < DIM; b++) {
            if ( pd->e[R_EM_E1_REAL+b] ) {
              int var = EM_E1_IMAG + b;
              int pvar_imag = upd->vp[var];
              for (int j = 0; j < ei->dof[var]; j++) {
                double diffusion_real = 0;
                double diffusion_imag = 0;
                for (int q = 0; q < DIM; q++) {
                  diffusion_imag += -delta(a,b)*bf[eqn_imag]->grad_phi[i][q]*bf[eqn_imag]->grad_phi[j][q];

                }

                // R = curl(E)
                //diffusion_imag += bf[eqn_imag]->phi[i]*bf[eqn_imag]->curl_phi_e[j][b][a];

                // Jacobian Contribution for R = E
                //diffusion_imag += delta(a,b)*bf[eqn_imag]->phi[i]*bf[eqn_imag]->phi[j];

                // Jacobian Contribution for R = curl(E)
                //diffusion_real += bf[eqn_real]->phi[i]*bf[eqn_real]->curl_phi_e[j][b][a];

                double stab_real = 0;
                double stab_imag = 0;
                if (!(pd->v[R_EM_CONT_REAL] && pd->v[R_EM_CONT_IMAG])) {
                  stab_imag = stab_scale * bf[eqn_real]->grad_phi[i][a] * (bf[var]->grad_phi_e[j][b][0][0] +
                      bf[var]->grad_phi_e[j][b][1][1]  + bf[var]->grad_phi_e[j][b][2][2]);
                }

                double advection_real = 0;
                double advection_imag = 0;

                advection_real = bf[eqn_real]->phi[i] * (delta(a,b) *  -im_coeff * bf[var]->phi[j]);
                advection_imag = bf[eqn_imag]->phi[i] * (delta(a,b) * re_coeff * bf[var]->phi[j]);
                lec->J[LEC_J_INDEX(peqn_real,pvar_imag,i,j)] += (diffusion_real + advection_real + stab_real) * bf[eqn_real]->detJ * fv->wt * fv->h3;
                lec->J[LEC_J_INDEX(peqn_imag,pvar_imag,i,j)] += (diffusion_imag + advection_imag + stab_imag) * bf[eqn_imag]->detJ * fv->wt * fv->h3;
              }
            }
          }

          {
            int var = EM_CONT_REAL;
            if (pd->e[var]) {
              int pvar = upd->vp[var];
              for (int j = 0; j < ei->dof[var]; j++) {
                //dbl lagr_term_real = bf[eqn]->grad_phi[i][a] * r_elperm * bf[var]->phi[j];
                //dbl lagr_term_imag = bf[eqn]->grad_phi[i][a] * i_elperm * bf[var]->phi[j];
                dbl lagr_term_real = bf[eqn]->grad_phi[i][a] * bf[var]->phi[j];
                dbl lagr_term_imag = 0;
                lec->J[LEC_J_INDEX(peqn_real,pvar,i,j)] += (lagr_term_real) * bf[eqn_real]->detJ * fv->wt * fv->h3;
                lec->J[LEC_J_INDEX(peqn_imag,pvar,i,j)] += (lagr_term_imag) * bf[eqn_real]->detJ * fv->wt * fv->h3;
              }
            }
          }

          {
            int var = EM_CONT_IMAG;
            if (pd->e[var]) {
              int pvar = upd->vp[var];
              for (int j = 0; j < ei->dof[var]; j++) {
                //dbl lagr_term_real = bf[eqn]->grad_phi[i][a] * i_elperm * bf[var]->phi[j];
                //dbl lagr_term_imag = bf[eqn]->grad_phi[i][a] * r_elperm * bf[var]->phi[j];
                dbl lagr_term_real = 0;
                dbl lagr_term_imag = bf[eqn]->grad_phi[i][a] * bf[var]->phi[j];
                lec->J[LEC_J_INDEX(peqn_real,pvar,i,j)] += (lagr_term_real) * bf[eqn_real]->detJ * fv->wt * fv->h3;
                lec->J[LEC_J_INDEX(peqn_imag,pvar,i,j)] += (lagr_term_imag) * bf[eqn_real]->detJ * fv->wt * fv->h3;
              }
            }
          }
        }
      }
    }
  }
  return(0);
} // end of assemble_ewave_laplacian


int assemble_em_continuity() {
  if (! (pd->e[R_EM_CONT_REAL] && pd->e[R_EM_CONT_IMAG])) {
    return -1;
  }

  dbl n;				/* Refractive index. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_n_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_n = &d_n_struct;

  dbl k;				/* Acoustic wave number. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  n = refractive_index( d_n, 0 );
  k = extinction_index( d_k, 0 );

  // Compute complex material properties
  complex double cpx_refractive_index, cpx_rel_permittivity,
      cpx_permittivity;//, impedance;
  double r_elperm, i_elperm;

  cpx_refractive_index = n + _Complex_I*k; // k > 0 is extinction
  cpx_rel_permittivity = SQUARE(cpx_refractive_index);
  cpx_permittivity = cpx_rel_permittivity*mp->permittivity;

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
  if ( af->Assemble_Residual ) {
    for (int i = 0; i < ei->dof[R_EM_CONT_REAL]; i++) {
      int eqn_real = R_EM_CONT_REAL;
      int eqn_imag = R_EM_CONT_IMAG;
      int peqn_real = upd->ep[eqn_real];
      int peqn_imag = upd->ep[eqn_imag];

      double advection_real = - factor*bf[eqn_real]->phi[i] * (r_elperm * div_E_real + i_elperm * div_E_imag);
      double advection_imag = - factor*bf[eqn_imag]->phi[i] * (r_elperm * div_E_imag + i_elperm * div_E_real);

      lec->R[LEC_R_INDEX(peqn_real,i)] += (advection_real) * bf[eqn_real]->detJ * fv->wt * fv->h3;
      lec->R[LEC_R_INDEX(peqn_imag,i)] += (advection_imag) * bf[eqn_imag]->detJ * fv->wt * fv->h3;
    }
  }
  if ( af->Assemble_Jacobian ) {
    for (int i = 0; i < ei->dof[R_EM_CONT_REAL]; i++) {
      int eqn_real = EM_CONT_REAL;
      int eqn_imag = EM_CONT_IMAG;
      int peqn_real = upd->ep[eqn_real];
      int peqn_imag = upd->ep[eqn_imag];

      for (int b = 0; b < dim; b++) {
        int var = EM_E1_REAL + b;
        int pvar_real = upd->ep[var];
        for (int j = 0; j < ei->dof[var]; j++) {
          double advection_real = 0;
          double advection_imag = 0;
          for (int q = 0; q < dim; q++) {
            advection_real += - factor*bf[eqn_real]->phi[i] * r_elperm * bf[var]->grad_phi_e[j][b][q][q];
            advection_imag += - factor*bf[eqn_imag]->phi[i] * i_elperm * bf[var]->grad_phi_e[j][b][q][q];
          }

          lec->J[LEC_J_INDEX(peqn_real,pvar_real,i,j)] += (advection_real) * bf[eqn_real]->detJ * fv->wt * fv->h3;
          lec->J[LEC_J_INDEX(peqn_imag,pvar_real,i,j)] += (advection_imag) * bf[eqn_imag]->detJ * fv->wt * fv->h3;
        }
      }

      for (int b = 0; b < dim; b++) {
        int var = EM_E1_IMAG + b;
        int pvar_imag = upd->ep[var];
        for (int j = 0; j < ei->dof[var]; j++) {
          double advection_real = 0;
          double advection_imag = 0;
          for (int q = 0; q < dim; q++) {
            advection_real += - factor*bf[eqn_real]->phi[i] * i_elperm * bf[var]->grad_phi_e[j][b][q][q];
            advection_imag += - factor*bf[eqn_imag]->phi[i] * r_elperm * bf[var]->grad_phi_e[j][b][q][q];
          }

          lec->J[LEC_J_INDEX(peqn_real,pvar_imag,i,j)] += (advection_real) * bf[eqn_real]->detJ * fv->wt * fv->h3;
          lec->J[LEC_J_INDEX(peqn_imag,pvar_imag,i,j)] += (advection_imag) * bf[eqn_imag]->detJ * fv->wt * fv->h3;
        }
      }
    }
  }
  return(0);
} // end of assemble_em_continuity

////////////////////////////////////////////////////////////////////////////////
/// Boundary Conditions ////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int apply_em_farfield_direct_vec(double func[DIM],
                double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                double xi[DIM],        /* Local stu coordinates */
                const int bc_name,
                double *bc_data) {
  /***********************************************************************
   * TODO AMC: rewrite this description
   * apply_em_farfield_direct():
   *
   *  Function which evaluates the expression specifying the
   *  a plane-wave directly incident (parallel to normal) on
   *  the boundary.
   *
   *  n_bound CROSS E = -n_bound CROSS E
   *                  * ( eta_2 * kappa_2)
   *                  / ( omega * mu_2 )
   *                  * ( 1 - 2*I)
   *
   *   func = -
   *
   *
   *  The boundary condition EM_DIRECT_BC employs this
   *  function.
   *
   *
   * Input:
   *
   *  em_eqn   = which equation is this applied to
   *  em_var   = which variable is this value sensitive to
   *
   * Output:
   *
   *  func[0] = value of the function mentioned above
   *  d_func[0][varType][lvardof] =
   *              Derivate of func[0] wrt
   *              the variable type, varType, and the local variable
   *              degree of freedom, lvardof, corresponding to that
   *              variable type.
   *
   *   Author: Andrew Cochrane (10/9/2019)
   *
   ********************************************************************/

  int var;
  dbl mag_permeability=12.57e-07; // H/m
  double n1, n2;				/* Refractive index */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_n1_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_n1 = &d_n1_struct;

  double k1, k2;				/* Extinction Coefficient */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k1_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k1 = &d_k1_struct;

  double normal[DIM]; // surface normal vector

  //double complex E1[DIM], H1[DIM]; // complex fields inside domain

  // Need material properties for both sides of interface

  // use mp for inside (subscript 1) ..

  double omega = upd->Acoustic_Frequency;
  n1 = refractive_index( d_n1, 0.0 );

  k1 = extinction_index( d_k1, 0.0 );
  // Compute complex impedance
  complex cpx_refractive_index1, cpx_rel_permittivity1,
      cpx_permittivity1, impedance1, kappa1;

  cpx_refractive_index1 = n1 + _Complex_I*k1;
  cpx_rel_permittivity1 = SQUARE(cpx_refractive_index1);
  cpx_permittivity1 = cpx_rel_permittivity1*mp->permittivity;

  impedance1 = csqrt(mag_permeability/cpx_permittivity1);
  kappa1 = omega*impedance1*cpx_permittivity1;

  // use BC input for outside (subscript 2)
  n2 = bc_data[0];
  k2 = bc_data[1];
  //n2 = 1.000293; // air (wikipedia 2019)
  //k2 = 0;
  // Compute complex impedance
  complex cpx_refractive_index2, cpx_rel_permittivity2,
      cpx_permittivity2, impedance2;

  cpx_refractive_index2 = n2 + _Complex_I*k2;
  cpx_rel_permittivity2 = SQUARE(cpx_refractive_index2);
  cpx_permittivity2 = cpx_rel_permittivity2*mp->permittivity;

  impedance2 = csqrt(mag_permeability/cpx_permittivity2);

  // need Surface Normal vector
  for (int p=0; p<DIM; p++) {
    normal[p] = fv->snormal[p];
  }

  complex Gamma, tau, incidentE[DIM];
  complex incidentH[DIM] = {0.0};
  Gamma = (impedance2 - impedance1)/(impedance2 + impedance1);
  tau = (2.0*impedance2)/(impedance2 + impedance1);

  double complex reduction_factor;

  switch (bc_name) {
    case EM_ER_FARFIELD_DIRECT_BC:
    case EM_EI_FARFIELD_DIRECT_BC:
      reduction_factor = -_Complex_I*tau/kappa1/(1+Gamma);
      break;
    case EM_HR_FARFIELD_DIRECT_BC:
    case EM_HI_FARFIELD_DIRECT_BC:
      reduction_factor = -_Complex_I*tau/kappa1/(1-Gamma);
      break;
    default:
      reduction_factor = 0.0;
      break;
  }


  incidentE[0] = bc_data[2] + _Complex_I*bc_data[5];
  incidentE[1] = bc_data[3] + _Complex_I*bc_data[6];
  incidentE[2] = bc_data[4] + _Complex_I*bc_data[7];

  for (int p=0; p<DIM; p++) {
    for (int q=0; q<DIM; q++) {
      for (int r=0; r<DIM; r++) {
        incidentH[p] += permute(p,q,r)*normal[q]*incidentE[r]/impedance2;
      }
    }
  }

  // construct curls and sensitivities
  // assuming normal is purely real.
  double Re_curl_E[DIM] = {0.0};
  double Im_curl_E[DIM] = {0.0};
  double Re_curl_H[DIM] = {0.0};
  double Im_curl_H[DIM] = {0.0};
  double d_dERb_Re_curl_E[DIM][DIM][MDE] = {{{0.0}}};
  double d_dEIb_Im_curl_E[DIM][DIM][MDE] = {{{0.0}}};
  double d_dHRb_Re_curl_H[DIM][DIM][MDE] = {{{0.0}}};
  double d_dHIb_Im_curl_H[DIM][DIM][MDE] = {{{0.0}}};

  for (int p=0; p<pd->Num_Dim; p++) {
    for (int q=0; q<pd->Num_Dim; q++) {
      for (int r=0; r<pd->Num_Dim; r++) {
        Re_curl_E[p] += permute(p,q,r)*fv->grad_em_er[r][q];
        Im_curl_E[p] += permute(p,q,r)*fv->grad_em_ei[r][q];
        Re_curl_H[p] += permute(p,q,r)*fv->grad_em_hr[r][q];
        Im_curl_H[p] += permute(p,q,r)*fv->grad_em_hi[r][q];
        // assuming all variables have same degrees of freedom

        for (int b=0; b<ei->dof[EM_E1_REAL]; b++) {
          d_dERb_Re_curl_E[p][r][b] += permute(p,q,r)
                                     *bf[EM_E1_REAL + r]->grad_phi[b][q];
        }
        for (int b=0; b<ei->dof[EM_E1_IMAG]; b++) {
          d_dEIb_Im_curl_E[p][r][b] += permute(p,q,r)
                                      *bf[EM_E1_IMAG + r]->grad_phi[b][q];
        }
        for (int b=0; b<ei->dof[EM_H1_REAL]; b++) {
          d_dHRb_Re_curl_H[p][r][b] += permute(p,q,r)
                                      *bf[EM_H1_REAL + r]->grad_phi[b][q];
        }
        for (int b=0; b<ei->dof[EM_H1_IMAG]; b++) {
          d_dHIb_Im_curl_H[p][r][b] += permute(p,q,r)
                                      *bf[EM_H1_IMAG + r]->grad_phi[b][q];

        }
      }
    }
  }

  complex cpx_func[DIM] = {0.0};

  //double real, imag;
  switch(bc_name) {
    case EM_ER_FARFIELD_DIRECT_BC:
    case EM_EI_FARFIELD_DIRECT_BC:
      for (int p=0; p<DIM; p++) {
        cpx_func[p] += (Re_curl_E[p]+ _Complex_I*Im_curl_E[p])
                      *reduction_factor;
      }
      for (int p=0; p<pd->Num_Dim; p++) {
        for (int q=0; q<DIM; q++) {
          for (int r=0; r<DIM; r++) {
            cpx_func[p] += permute(p,q,r)*normal[q]*incidentE[r];
          }
        }
      }
      break;
    case EM_HR_FARFIELD_DIRECT_BC:
    case EM_HI_FARFIELD_DIRECT_BC:
      for (int p=0; p<DIM; p++) {
        cpx_func[p] += (Re_curl_H[p] + _Complex_I*Im_curl_H[p])
                      *reduction_factor;
      }
      for (int p=0; p<pd->Num_Dim; p++) {
        for (int q=0; q<DIM; q++) {
          for (int r=0; r<DIM; r++) {
            cpx_func[p] -= permute(p,q,r)*normal[q]*incidentH[r];
          }
        }
      }
      break;
  }

  switch(bc_name) {
    case EM_ER_FARFIELD_DIRECT_BC:

      for (int p=0; p<DIM; p++) {
        func[p] = creal(cpx_func[p]);
        //func[0] = fv->grad_em_er[1][2] - fv->grad_em_er[2][1];
        /*
        for(int q=0; q<DIM; q++) {
          for (int r=0; r<DIM; r++){
            func[p] += permute(p,q,r)*fv->grad_em_er[r][q];
          }
        }*/
      }
      //eqn = R_EM_H*_REAL;
      var = EM_E1_REAL;
    //  real = 1.0;
    //  imag = 0.0;
      break;
    case EM_EI_FARFIELD_DIRECT_BC:
      for (int p=0; p<DIM; p++) {
        func[p] = cimag(cpx_func[p]);
      }
      //eqn = R_EM_H*_IMAG;
      var = EM_E1_REAL;
      //real = 0.0;
      //imag = 1.0;
      break;
    case EM_HR_FARFIELD_DIRECT_BC:

      for (int p=0; p<DIM; p++) {
        func[p] = creal(cpx_func[p]);
      }

/*
      for (int p=0; p<DIM; p++) {
        func[p] = Im_n_x_H[p];
      }
      */
      //eqn = R_EM_E*_REAL;
      var = EM_H1_REAL;
      //real = 1.0;
      //imag = 0.0;
      break;
    case EM_HI_FARFIELD_DIRECT_BC:
      for (int p=0; p<DIM; p++) {
        func[p] = cimag(cpx_func[p]);
      }
      //eqn = R_EM_E*_IMAG;
      var = EM_H1_REAL;
      //real = 0.0;
      //imag = 1.0;
      break;
    default:
      var = 0;
      //real = 0;
      //imag = 0;
      reduction_factor = 0;
      EH(-1, "Must call apply_em_farfield_direct with an applicable BC_NAME");
      return -1;
      break;
  }

  if(af->Assemble_Jacobian) {
/*
    for (int j=0; j< ei->dof[EM_E2_REAL]; j++){
          d_func[0][EM_E2_REAL][j] = -bf[EM_E2_REAL]->grad_phi[j][2];
          d_func[0][EM_E3_REAL][j] =  bf[EM_E3_REAL]->grad_phi[j][1];
                 0        1           -        1                  2
                 0        2           +        2                  1
                 p        r           +        r                  q
    }
  */
    /*
    for (int p=0; p<DIM; p++) {

      //for (int q=0; q<DIM; q++) {
        for (int r=0; r<DIM; r++) {
          for (int j=0; j<ei->dof[EM_E1_REAL + r]; j++) {
            d_func[p][EM_E1_REAL + r][j] += -d_dERb_Re_curl_E[p][r][j];
          }
        }
      //}
    }
//    d_dERb_Re_curl_E[p][r][b]
    */
    for (int p=0; p<pd->Num_Dim; p++) {
      //for (int q=0; q<pd->Num_Dim; q++) {
        for (int g=0; g<pd->Num_Dim; g++) {
          int gvar = var + g;
          for (int j=0; j<ei->dof[gvar]; j++) {
            switch (bc_name){
              case EM_ER_FARFIELD_DIRECT_BC:
                //d_func[p][gvar][j] += d_dERb_Re_curl_E[p][g][j]*creal(reduction_factor);
                //d_func[p][gvar+3][j] -= d_dEIb_Im_curl_E[p][g][j]*cimag(reduction_factor);
                d_func[p][gvar][j] += -d_dERb_Re_curl_E[p][g][j]*creal(reduction_factor);
                d_func[p][gvar+3][j] += d_dERb_Re_curl_E[p][g][j]*cimag(reduction_factor);
                break;

              case EM_EI_FARFIELD_DIRECT_BC:
                d_func[p][gvar][j] += -d_dERb_Re_curl_E[p][g][j]*cimag(reduction_factor);
                d_func[p][gvar+3][j] += -d_dEIb_Im_curl_E[p][g][j]*creal(reduction_factor);
                break;

              case EM_HR_FARFIELD_DIRECT_BC:
                d_func[p][gvar][j] += -d_dHRb_Re_curl_H[p][g][j]*creal(reduction_factor);
                d_func[p][gvar+3][j] += d_dHIb_Im_curl_H[p][g][j]*cimag(reduction_factor);
                //d_func[p][gvar][j] += d_dHRb_Re_n_x_H[p][g][j];
                //d_func[p][gvar+3][j] += d_dERb_Re_n_x_E[p][g][j];
                break;

              case EM_HI_FARFIELD_DIRECT_BC:
                d_func[p][gvar][j] += -d_dHRb_Re_curl_H[p][g][j]*cimag(reduction_factor);
                d_func[p][gvar+3][j] += -d_dHIb_Im_curl_H[p][g][j]*creal(reduction_factor);
                break;

            }
          }
        //}
      }
    }
  }
  return 0;
} // end of apply_em_direct_vec





/* AMC TODO:  Here's the working cross product w/sensitivity
  for (int p=0; p<pd->Num_Dim; p++) {
    for (int q=0; q<pd->Num_Dim; q++) {
      for (int r=0; r<pd->Num_Dim; r++) {
        //func[p] += fv->em_er[q] + fv->em_ei[q] + fv->em_hr[q] + fv->em_hi[q];
        func[p] += permute(p,q,r)* creal(normal[q]*E1[r]);
      }
    }
  }

  if(af->Assemble_Jacobian) {

    for (int p=0; p<pd->Num_Dim; p++) {
      for (int q=0; q<pd->Num_Dim; q++) {
        for (int g=0; g<3; g++) {
          int gvar = EM_E1_REAL + g;
          for (int j=0; j<ei->dof[gvar]; j++) {


            //d_func[p][gvar][j] = phi_j;
            d_func[p][gvar][j] += permute(q,p,g)*creal(normal[q])*phi_j;

          }
        }
      }
    }
*/



  /*
    switch(bc_name) {
      case EM_ER_FARFIELD_DIRECT_BC:
      case EM_EI_FARFIELD_DIRECT_BC:
      case EM_HR_FARFIELD_DIRECT_BC:
      case EM_HI_FARFIELD_DIRECT_BC:
        for (int g=0; g<pd->Num_Dim; g++){
          int gvar = var + g;
          for (int j=0; j<ei->dof[gvar]; j++) {
            double phi_j = bf[gvar]->phi[j];
            for (int p=0; p<pd->Num_Dim; p++) {
              for (int q=0; q<pd->Num_Dim; q++) {
                //for (int r=0; r<pd->Num_Dim; r++) {
                  d_func[p][gvar][j] += real*creal(permute(p,q,g)
                                                   *reduction_factor
                                                   *normal[q]
                                                   *phi_j)
                                     +  imag*cimag(permute(p,q,g)
                                                   *reduction_factor
                                                   *normal[q]
                                                   *phi_j);
                //}
              }
            }
          }
        }
        break;

    }

  }*/

int apply_em_free_vec(double func[DIM],
                double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                double xi[DIM],        /* Local stu coordinates */
                const int bc_name) {

  double normal[DIM]; // surface normal vector

  // need Surface Normal vector
  for (int p=0; p<DIM; p++) {
    normal[p] = fv->snormal[p];
  }

  // Evaluate n cross E or n cross H
  switch (bc_name) {
    case EM_ER_FREE_BC:
      for (int p=0; p<DIM; p++) {
        for (int q=0; q<DIM; q++) {
          for (int r=0; r<DIM; r++) {
            func[p] += permute(p,q,r)*normal[q]*fv->em_er[r];
          }
        }
      }
      break;
    case EM_EI_FREE_BC:
      for (int p=0; p<DIM; p++) {
        for (int q=0; q<DIM; q++) {
          for (int r=0; r<DIM; r++) {
            func[p] += permute(p,q,r)*normal[q]*fv->em_ei[r];
          }
        }
      }
      break;
    case EM_HR_FREE_BC:
      for (int p=0; p<DIM; p++) {
        for (int q=0; q<DIM; q++) {
          for (int r=0; r<DIM; r++) {
            func[p] += permute(p,q,r)*normal[q]*fv->em_hr[r];
          }
        }
      }
      break;
    case EM_HI_FREE_BC:
      for (int p=0; p<DIM; p++) {
        for (int q=0; q<DIM; q++) {
          for (int r=0; r<DIM; r++) {
            func[p] += permute(p,q,r)*normal[q]*fv->em_hi[r];
          }
        }
      }
      break;
  }
  return 0;
}

int apply_em_sommerfeld_vec(double func[DIM],
                double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                double xi[DIM],        /* Local stu coordinates */
                const int bc_name,
                double *bc_data) {
  /***********************************************************************
   * TODO AMC: rewrite this description
   * apply_em_sommerfeld():
   *
   *  Function for specifying incoming plane wave and outgoing scattered
   *  energy using a radiative boundary condition sometimes called
   *  Sommerfeld radiation condition
   *
   *  n_bound CROSS E = j/kappa * Curl(E - E_i)
   *                  + n_bound CROSS E_i
   *
   *   func =
   *
   *
   *  The boundary condition EM_ employs this
   *  function.
   *
   *
   * Input:
   *
   *  em_eqn   = which equation is this applied to
   *  em_var   = which variable is this value sensitive to
   *
   * Output:
   *
   *  func[0] = value of the function mentioned above
   *  d_func[0][varType][lvardof] =
   *              Derivate of func[0] wrt
   *              the variable type, varType, and the local variable
   *              degree of freedom, lvardof, corresponding to that
   *              variable type.
   *
   *   Author: Andrew Cochrane (10/9/2019)
   *
   ********************************************************************/
  dbl mag_permeability=12.57e-07; // H/m

  double impedance = sqrt(mag_permeability/mp->permittivity);
  double omega = upd->Acoustic_Frequency;
  double n[DIM]; // surface normal vector

  //double complex E1[DIM], H1[DIM]; // complex fields inside domain

  // This BC assumes that the boundary is far from the subject and
  // the material properties are the same on both sides
  // Need the wave number

  double kappa = omega*sqrt(mp->permittivity*mag_permeability);

  // need Surface Normal vector
  for (int p=0; p<DIM; p++) {
    n[p] = fv->snormal[p];
  }

  // polarization P and propagation direction k of incident plane wave
  // from input deck
  complex double P[DIM];
  double k[DIM];

  P[0] = bc_data[0] + _Complex_I*bc_data[3];
  P[1] = bc_data[1] + _Complex_I*bc_data[4];
  P[2] = bc_data[2] + _Complex_I*bc_data[5];

  k[0] = bc_data[6];
  k[1] = bc_data[7];
  k[2] = bc_data[8];

  // normalize k and use wavenumber kappa
  double k_mag = sqrt(k[0]*k[0] + k[1]*k[1] + k[2]*k[2]);

  k[0] = k[0]/k_mag;
  k[1] = k[1]/k_mag;
  k[2] = k[2]/k_mag;

  // Compute E_i
  // E_i = P*exp(i(k DOT x - wt)
  double x[DIM];
  x[0] = fv->x[0];
  x[1] = fv->x[1];
  x[2] = fv->x[2];


  double kappa_dot_x = 0.0;
  for (int q=0; q<DIM; q++){
    kappa_dot_x += k[q]*x[q];
  }
  kappa_dot_x *= kappa;

  complex double E_i[DIM] = {0.0};
  complex double CurlE_i[DIM] = {0.0};
  complex double nCrossE_i[DIM] = {0.0};
  complex double CurlE[DIM] = {0.0};
  complex double H_i[DIM] = {0.0};
  complex double CurlH_i[DIM] = {0.0};
  complex double nCrossH_i[DIM] = {0.0};
  complex double CurlH[DIM] = {0.0};

  for (int p=0; p<DIM; p++) {
    E_i[p] += P[p] * cexp(_Complex_I*kappa_dot_x);
  }

  switch (bc_name) {
    case EM_ER_SOMMERFELD_BC:
    case EM_EI_SOMMERFELD_BC:

      // Compute Curl(E_i)
      // [Curl(E_i)]_e = sum_{f,g} [permute(e,f,g)*[d_dx]_f([E_i]_g)]
      for (int e=0; e<DIM; e++) {
        for (int f=0; f<DIM; f++) {
          for (int g=0; g<DIM; g++) {
            CurlE_i[e] += permute(e,f,g)
                       *_Complex_I*kappa*k[f]*P[g]
                       *cexp(_Complex_I*kappa_dot_x);
          }
        }
      }

      // Compute n CROSS E_i
      for (int e=0; e<DIM; e++) {
        for (int f=0; f<DIM; f++) {
          for (int g=0; g<DIM; g++) {
            nCrossE_i[e] += permute(e,f,g)*n[f]*E_i[g];
          }
        }
      }

      // Compute Curl E
      for (int e=0; e<DIM; e++) {
        for (int f=0; f<DIM; f++) {
          for (int g=0; g<DIM; g++) {
            CurlE[e] += permute(e,f,g)
                        *(fv->grad_em_er[g][f]
                          + _Complex_I*fv->grad_em_ei[g][f]);
          }
        }
      }
      break;

    case EM_HR_SOMMERFELD_BC:
    case EM_HI_SOMMERFELD_BC:
      for (int e=0; e<DIM; e++) {
        for (int f=0; f<DIM; f++) {
          for (int g=0; g<DIM; g++) {
            H_i[e] += permute(e,f,g)*k[f]
                      *E_i[g]/impedance;
          }
        }
      }

      // Compute Curl(H_i)
      // [Curl(H_i)]_e = sum_{f,g} [permute(e,f,g)*[d_dx]_f([H_i]_g)]
      for (int e=0; e<DIM; e++) {
        for (int f=0; f<DIM; f++) {
          for (int g=0; g<DIM; g++) {
            CurlH_i[e] += permute(e,f,g)
                       *_Complex_I*kappa*k[f]*H_i[g];
          }
        }
      }

      // Compute n CROSS H_i
      for (int e=0; e<DIM; e++) {
        for (int f=0; f<DIM; f++) {
          for (int g=0; g<DIM; g++) {
            nCrossH_i[e] += permute(e,f,g)*n[f]*H_i[g];
          }
        }
      }

      // Compute Curl H
      for (int e=0; e<DIM; e++) {
        for (int f=0; f<DIM; f++) {
          for (int g=0; g<DIM; g++) {
            CurlE[e] += permute(e,f,g)
                        *(fv->grad_em_hr[g][f]
                          + _Complex_I*fv->grad_em_hi[g][f]);
          }
        }
      }
      break;
  }

  // Residual Components

  switch (bc_name) {
    case EM_ER_SOMMERFELD_BC:
      for (int p=0; p<DIM; p++) {
        func[p] -= creal(_Complex_I/kappa
                        *(CurlE[p] - CurlE_i[p])
                        + nCrossE_i[p]
                        );
      }
      break;
    case EM_EI_SOMMERFELD_BC:
      for (int p=0; p<DIM; p++) {
        func[p] -= cimag(_Complex_I/kappa
                        *(CurlE[p] - CurlE_i[p])
                        + nCrossE_i[p]
                        );
      }
      break;
    case EM_HR_SOMMERFELD_BC:
      for (int p=0; p<DIM; p++) {
        func[p] -= creal(_Complex_I/kappa
                        *(CurlH[p] - CurlH_i[p])
                        + nCrossH_i[p]
                        );
      }
      break;
    case EM_HI_SOMMERFELD_BC:
      for (int p=0; p<DIM; p++) {
        func[p] -= cimag(_Complex_I/kappa
                        *(CurlH[p] - CurlH_i[p])
                        + nCrossH_i[p]
                        );
      }
      break;
  }
  return 0;
} // end of apply_em_sommerfeld_vec

int apply_ewave_planewave_vec(double func[DIM],
                double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                double xi[DIM],        /* Local stu coordinates */
                const int bc_name,
                double *bc_data) {
  /***********************************************************************
   * TODO AMC: rewrite this description + add Jacobians
   * apply_ewave_planewave_vec():
   *
   *  Function for specifying surface integral terms for plane wave test.
   *
   *  Surface Integral = n_bound CROSS Curl E_i dS
   *
   *   func =
   *
   *
   *  The boundary condition E_PLANEWAVE employs this
   *  function.
   *
   *
   * Input:
   *
   *  em_eqn   = which equation is this applied to
   *  em_var   = which variable is this value sensitive to
   *
   * Output:
   *
   *  func[0] = value of the function mentioned above
   *  d_func[0][varType][lvardof] =
   *              Derivate of func[0] wrt
   *              the variable type, varType, and the local variable
   *              degree of freedom, lvardof, corresponding to that
   *              variable type.
   *
   *   Author: Andrew Cochrane (5/19/2020)
   *
   ********************************************************************/

  //double impedance = sqrt(mp->magnetic_permeability/mp->permittivity);
  double omega = upd->Acoustic_Frequency;
  double n[DIM]; // surface normal vector

  //double complex E1[DIM], H1[DIM]; // complex fields inside domain

  // This BC assumes that the boundary is far from the subject and
  // the material properties are the same on both sides
  // Need the wave number

  double kappa = omega*sqrt(mp->permittivity*mp->magnetic_permeability);

  // need Surface Normal vector
  for (int p=0; p<DIM; p++) {
    n[p] = fv->snormal[p];
  }

  // polarization P, propagation direction k of incident plane wave
  // as well as origin offset from input deck
  complex double P[DIM];
  double k[DIM];
  double x_orig[DIM];

  P[0] = bc_data[0] + _Complex_I*bc_data[3];
  P[1] = bc_data[1] + _Complex_I*bc_data[4];
  P[2] = bc_data[2] + _Complex_I*bc_data[5];

  k[0] = bc_data[6];
  k[1] = bc_data[7];
  k[2] = bc_data[8];

  x_orig[0] = bc_data[9];
  x_orig[1] = bc_data[10];
  x_orig[2] = bc_data[11];

  // normalize k and use wavenumber kappa
  double k_mag = sqrt(k[0]*k[0] + k[1]*k[1] + k[2]*k[2]);

  k[0] = k[0]/k_mag;
  k[1] = k[1]/k_mag;
  k[2] = k[2]/k_mag;

  double x[DIM];
  x[0] = fv->x[0] - x_orig[0];
  x[1] = fv->x[1] - x_orig[1];
  x[2] = fv->x[2] - x_orig[2];


  double kappa_dot_x = 0.0;
  for (int q=0; q<DIM; q++){
    kappa_dot_x += k[q]*x[q];
  }
  kappa_dot_x *= kappa;

  complex double E[DIM] = {0.0};
  complex double CurlE[DIM] = {0.0};
  complex double nCrossCurlE[DIM] = {0.0};

  for (int p=0; p<DIM; p++) {
    E[p] += P[p] * cexp(_Complex_I*kappa_dot_x);
  }

    // Compute Curl(E)
    // [Curl(E)]_e = sum_{f,g} [permute(e,f,g)*[d_dx]_f([E]_g)]
    for (int e=0; e<DIM; e++) {
      for (int f=0; f<DIM; f++) {
        for (int g=0; g<DIM; g++) {
          CurlE[e] += permute(e,f,g)
                     *kappa*k[f]*P[g]
                     *cexp(_Complex_I*kappa_dot_x);
        }
      }
    }

    // Compute n CROSS Curl(E)
    for (int e=0; e<DIM; e++) {
      for (int f=0; f<DIM; f++) {
        for (int g=0; g<DIM; g++) {
          nCrossCurlE[e] += permute(e,f,g)*n[f]*CurlE[g];
        }
      }
    }

  // Residual Components

  switch (bc_name) {
    case E_ER_PLANEWAVE_BC:
      for (int p=0; p<DIM; p++) {
        func[p] += creal(nCrossCurlE[p]);
      }
      break;
    case E_EI_PLANEWAVE_BC:
    for (int p=0; p<DIM; p++) {
      func[p] += cimag(nCrossCurlE[p]);
    }
    break;
  }
  return 0;
} // end of apply_ewave_planewave_vec

int apply_ewave_curlcurl_farfield_vec(double func[DIM],
                double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                double xi[DIM],        /* Local stu coordinates */
                double time,           // present time
                const int bc_name,
                double *bc_data) {
  /***********************************************************************
   * TODO AMC: rewrite this description + add Jacobians
   * apply_ewave_curlcurl_farfield_vec():
   *
   *  Function for specifying surface integral terms for curlcurlE eqn.
   *
   *  Surface Integral = phi dot n_bound CROSS Curl E_1 dS
   *
   *   Derived Expression hasn't accounted for difference in polarization
   *   between E1 and thefarfield incident wave - for now assume they're
   *   the same and in the direction of the farfield incident wave
   *
   *
   *
   *   func = psi dot [ -i kappa_1 (E_1 - 2*E_f) eta_1/eta_2 ]
   *
   *
   *  The boundary condition E_FARFIELD employs this
   *  function.
   *
   *
   * Input:
   *
   *  em_eqn   = which equation is this applied to
   *  em_var   = which variable is this value sensitive to
   *
   * Output:
   *
   *  func[0] = value of the function mentioned above
   *  d_func[0][varType][lvardof] =
   *              Derivate of func[0] wrt
   *              the variable type, varType, and the local variable
   *              degree of freedom, lvardof, corresponding to that
   *              variable type.
   *
   *   Author: Andrew Cochrane (9/21/2020)
   *
   ********************************************************************/

  // need Surface Normal vector
/* maybe not
  double n_surf[DIM]; // surface normal vector
  for (int p=0; p<DIM; p++) {
    n_surf[p] = fv->snormal[p];
  }
*/
  // setup wave and material properties
  double omega = upd->Acoustic_Frequency;
  dbl n_1;				/* Refractive index. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_n_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_n = &d_n_struct;

  dbl k_1;				/* extinction coefficient. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  n_1 = refractive_index( d_n, time );
  k_1 = extinction_index( d_k, time );

  // Compute complex material properties
  complex double cpx_refractive_index_1, cpx_rel_permittivity_1,
      cpx_permittivity_1;//, impedance;

  cpx_refractive_index_1 = n_1 + _Complex_I*k_1; // k > 0 is extinction
  cpx_rel_permittivity_1 = SQUARE(cpx_refractive_index_1);
  cpx_permittivity_1 = cpx_rel_permittivity_1*mp->permittivity; // better set permittivity to vacuum in input deck?

  double complex E_1[DIM]; // complex field inside domain

  for (int p=0; p<DIM; p++) {
    E_1[p] = fv->em_er[p] + _Complex_I*fv->em_ei[p];
  }

  // Need the wave number
  complex double kappa_1 = omega*csqrt(cpx_permittivity_1*mp->magnetic_permeability);
  // Need the impedance
  complex double eta_1 = csqrt(mp->magnetic_permeability/cpx_permittivity_1);

  // polarization P of incident plane wave
  complex double P[DIM];

  // permittivity of the far-field medium
  double n_2, k_2;

  P[0] = bc_data[0] + _Complex_I*bc_data[3];
  P[1] = bc_data[1] + _Complex_I*bc_data[4];
  P[2] = bc_data[2] + _Complex_I*bc_data[5];

  // Need the impedance of the far-field material (superstrate)
  // use refractive index and extinction coefficient
  n_2 = bc_data[6];
  k_2 = bc_data[7];

  // Compute complex material properties
  complex double cpx_rel_permittivity_2, cpx_refractive_index_2,
      cpx_permittivity_2;

  cpx_refractive_index_2 = n_2 + _Complex_I*k_2; // k > 0 is extinction
  cpx_rel_permittivity_2 = SQUARE(cpx_refractive_index_2);
  cpx_permittivity_2 = cpx_rel_permittivity_2*mp->permittivity; // better set permittivity to vacuum in input deck?
  // Need the impedance
  complex double eta_2 = csqrt(mp->magnetic_permeability/cpx_permittivity_2);

  complex double cpx_coeff = -_Complex_I*kappa_1*eta_1/eta_2;

  // Residual Components

  switch (bc_name) {
    case E_ER_FARFIELD_BC:

    for (int p=0; p<pd->Num_Dim; p++) {
        func[p] = creal(cpx_coeff*(E_1[p] - 2*P[p]));
        for (int q=0; q<pd->Num_Dim; q++) {
          int qRvar = EM_E1_REAL + q;
          int qIvar = EM_E1_IMAG + q;
          for (int j=0; j<ei->dof[EM_E1_REAL]; j++) {
            double Rphi_j = bf[EM_E1_REAL+q]->phi[j];
            double Iphi_j = bf[EM_E1_IMAG+q]->phi[j];
            // d(Re[f{z}])/d(Re[z]) = Re(f'(z))
            d_func[p][qRvar][j] = delta(p,q)
                *creal(cpx_coeff*Rphi_j);

            // d(Re[f{z}])/d(Im[z]) = -Im(f'(z))
            d_func[p][qIvar][j] = -delta(p,q)
                *cimag(cpx_coeff*Iphi_j);
          }
        }
      }
      /*/
    for (int p=0; p<DIM; p++) {
        //func[p] += creal((E_1[p]));
      //func[p] = fv->em_er[p];
      func[p] = 1.0;
        for (int q=0; q<pd->Num_Dim; q++) {
          int qRvar = EM_E1_REAL + q;
          int qIvar = EM_E1_IMAG + q;
          for (int j=0; j<ei->dof[EM_E1R_BC]; j++) {
            double Rphi_j = bf[EM_E1_REAL+q]->phi[j];
            double Iphi_j = bf[EM_E1_IMAG+q]->phi[j];
            // d(Re[f{z}])/d(Re[z]) = Re(f'(z))
            //d_func[p][qRvar][j] = delta(p,q)*Rphi_j;
            d_func[p][qRvar][j] = 0.0;
                // *creal(Rphi_j);

            // d(Re[f{z}])/d(Im[z]) = -Im(f'(z))
            d_func[p][qIvar][j] += -delta(p,q)
                *0.0;
                // *cimag(Iphi_j);
          }
        }
      }
    */
      break;
    case E_EI_FARFIELD_BC:

    for (int p=0; p<DIM; p++) {
      func[p] += cimag(cpx_coeff*(E_1[p] - 2*P[p]));
        for (int q=0; q<pd->Num_Dim; q++) {
          int qRvar = EM_E1_REAL + q;
          int qIvar = EM_E1_IMAG + q;
          for (int j=0; j<ei->dof[EM_E1_IMAG]; j++) {
            double Rphi_j = bf[EM_E1_REAL+q]->phi[j];
            double Iphi_j = bf[EM_E1_IMAG+q]->phi[j];

            // d(Im[f{z}])/d(Re[z]) = Im(f'(z))
            d_func[p][qRvar][j] = delta(p,q)
                *cimag(cpx_coeff*Rphi_j);

            // d(Im[f{z}])/d(Im[z]) = Re(f'(z))
            d_func[p][qIvar][j] = delta(p,q)
                *creal(cpx_coeff*Iphi_j);
          }
        }
      }
      /*/
      for (int p=0; p<pd->Num_Dim; p++) {
        func[p] += cimag(cpx_coeff*(E_1[p]));
        //func[p] = fv->em_ei[p];
          for (int q=0; q<pd->Num_Dim; q++) {
            int qRvar = EM_E1_REAL + q;
            int qIvar = EM_E1_IMAG + q;
            for (int j=0; j<ei->dof[EM_E1_IMAG]; j++) {
              double Rphi_j = bf[EM_E1_REAL+q]->phi[j];
              double Iphi_j = bf[EM_E1_IMAG+q]->phi[j];
              // d(Im[f{z}])/d(Re[z]) = Im(f'(z))
              d_func[p][qRvar][j] = delta(p,q)
                  *cimag(cpx_coeff*Rphi_j);

              // d(Im[f{z}])/d(Im[z]) = Re(f'(z))
              d_func[p][qIvar][j] = delta(p,q)
                  *creal(cpx_coeff*Iphi_j);
            }
          }
        }
*/
    break;
  }

  return 0;
} // end of apply_ewave_farfield_vec

int apply_ewave_mms_vec(double func[DIM],
                double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                double xi[DIM],        /* Local stu coordinates */
                const int bc_name,
                double *bc_data) {
  /***********************************************************************
   * TODO AMC: rewrite this description + add Jacobians
   * apply_ewave_planewave_vec():
   *
   *  Function for specifying surface integral terms for plane wave test.
   *
   *  Surface Integral = n_bound CROSS Curl E_mms dS
   *
   *   func =
   *
   *
   *  The boundary condition E_MMS employs this
   *  function.
   *
   *
   * Input:
   *
   *  em_eqn   = which equation is this applied to
   *  em_var   = which variable is this value sensitive to
   *
   * Output:
   *
   *  func[0] = value of the function mentioned above
   *  d_func[0][varType][lvardof] =
   *              Derivate of func[0] wrt
   *              the variable type, varType, and the local variable
   *              degree of freedom, lvardof, corresponding to that
   *              variable type.
   *
   *   Author: Andrew Cochrane (7/15/2020)
   *
   ********************************************************************/

  double omega = upd->Acoustic_Frequency;
  double n[DIM]; // surface normal vector

  //double complex E1[DIM], H1[DIM]; // complex fields inside domain

  // This BC assumes that the boundary is far from the subject and
  // the material properties are the same on both sides
  // Need the wave number

  double kappa = omega*sqrt(mp->permittivity*mp->magnetic_permeability);

  // need Surface Normal vector
  for (int p=0; p<DIM; p++) {
    n[p] = fv->snormal[p];
  }

  // polarization P, propagation direction k of incident plane wave
  // as well as origin offset from input deck
  complex double P[DIM];
  double k[DIM];
  double x_orig[DIM];

  P[0] = bc_data[0] + _Complex_I*bc_data[3];
  P[1] = bc_data[1] + _Complex_I*bc_data[4];
  P[2] = bc_data[2] + _Complex_I*bc_data[5];

  k[0] = bc_data[6];
  k[1] = bc_data[7];
  k[2] = bc_data[8];

  x_orig[0] = bc_data[9];
  x_orig[1] = bc_data[10];
  x_orig[2] = bc_data[11];

  // normalize k and use wavenumber kappa
  double k_mag = sqrt(k[0]*k[0] + k[1]*k[1] + k[2]*k[2]);

  k[0] = k[0]/k_mag;
  k[1] = k[1]/k_mag;
  k[2] = k[2]/k_mag;

  double x[DIM];
  x[0] = fv->x[0] - x_orig[0];
  x[1] = fv->x[1] - x_orig[1];
  x[2] = fv->x[2] - x_orig[2];

  x[0] = fv->x[0];
  x[1] = fv->x[1];
  x[2] = fv->x[2];

  double kappa_dot_x = 0.0;
  for (int q=0; q<DIM; q++){
    kappa_dot_x += k[q]*x[q];
  }
  kappa_dot_x *= kappa;

  complex double E[DIM] = {0.0};
  complex double CurlE[DIM] = {0.0};
  complex double nCrossCurlE[DIM] = {0.0};

  for (int p=0; p<DIM; p++) {
    E[p] += P[p] * cexp(_Complex_I*kappa_dot_x);
  }

    // Compute Curl(E)
    // [Curl(E)]_e = sum_{f,g} [permute(e,f,g)*[d_dx]_f([E]_g)]
    for (int e=0; e<DIM; e++) {
      for (int f=0; f<DIM; f++) {
        for (int g=0; g<DIM; g++) {
          CurlE[e] += permute(e,f,g)
                     *kappa*k[f]*P[g]
                     *cexp(_Complex_I*kappa_dot_x);
        }
      }
    }

    // Compute n CROSS Curl(E)
    for (int e=0; e<DIM; e++) {
      for (int f=0; f<DIM; f++) {
        for (int g=0; g<DIM; g++) {
          nCrossCurlE[e] += permute(e,f,g)*n[f]*CurlE[g];
        }
      }
    }

  // Residual Components

  switch (bc_name) {
    case E_ER_PLANEWAVE_BC:
      for (int p=0; p<DIM; p++) {
        func[p] += creal(nCrossCurlE[p]);
      }
      break;
    case E_EI_PLANEWAVE_BC:
    for (int p=0; p<DIM; p++) {
      func[p] += cimag(nCrossCurlE[p]);
    }
    break;
  }
  return 0;
} // end of apply_ewave_mms_vec

int apply_ewave_2D(double func[DIM],
                double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                double xi[DIM],        /* Local stu coordinates */
                const int bc_name) {
  /***********************************************************************
   * TODO AMC: rewrite this description + add Jacobians
   * apply_ewave_2D():
   *
   *  strongly integrate the dEz_dz = 0
   *
   *  Surface Integral = d/dx_3( x_3 dot E) = 0 dS
   *
   *   func_r[0,1,2] = fv->grad_em_er[0,1,2][2]
   *   func_i[0,1,2] = fv->grad_em_ei[0,1,2][2]
   *
   *  The boundary condition E_PLANEWAVE employs this
   *  function.
   *
   *
   * Input:
   *
   *
   *
   * Output:
   *
   *  func[0,1,2] = value of the function mentioned above
   *  d_func[0][varType][lvardof] =
   *              Derivate of func[0] wrt
   *              the variable type, varType, and the local variable
   *              degree of freedom, lvardof, corresponding to that
   *              variable type.
   *
   *   Author: Andrew Cochrane (6/22/2020)
   *
   ********************************************************************/
  int var;
  double dphi_j_dx3;
  // Residual Components

  switch (bc_name) {
    case E_ER_2D_BC:
      //func[0] += fv->grad_em_er[2][0];
      //func[1] += fv->grad_em_er[2][1];
      func[2] += fv->grad_em_er[2][2];
      break;
    case E_EI_2D_BC:
      //func[0] += fv->grad_em_ei[2][0];
      //func[1] += fv->grad_em_ei[2][1];
      func[2] += fv->grad_em_ei[2][2];
    break;
  }

  // Jacobian Components
  if (af->Assemble_Jacobian) {
    switch (bc_name) {
    case E_ER_2D_BC:
      for (int k = 2; k<DIM; k++) {
        var = EM_E1_REAL + k;
        for (int j=0; j< ei->dof[var]; j++) {
          dphi_j_dx3 = bf[var]->grad_phi[j][2];
          d_func[k][var][j] += dphi_j_dx3;
        }
      }
      break;
    case E_EI_2D_BC:
      for (int k = 2; k<DIM; k++) {
        var = EM_E1_IMAG + k;
        for (int j=0; j< ei->dof[var]; j++) {
          dphi_j_dx3 = bf[var]->grad_phi[j][2];
          d_func[k][var][j] += dphi_j_dx3;
        }
      }
      break;
    }
  }
  return 0;
} // end of apply_ewave_2D


void
calc_emwave_stabilization_term(struct emwave_stabilization *em_stab,
                               double stabilization_coefficient
                               ){
  double complex grad_stabilization_field[DIM][DIM] = {{0.0}};
  dbl mag_permeability=12.57e-07;  // H/m


  /*
   *
   */
  switch(em_stab->type) {
    case EM_STAB_NONE:
      for( int i=0; i<ei->dof[em_stab->em_eqn]; i++){
        em_stab->residual_term[i] = 0.0;
        for (int b=0; b<DIM; b++){
          for (int j=0; j<ei->dof[em_stab->em_var]; j++){
            em_stab->jacobian_term[i][b][j]
                = 0.0;
          }
        }
      }
      return;
      break;
    case EM_STAB_PHI_DIV:
    case EM_STAB_DPHI_DIV:
    case EM_STAB_DIVPHI_DIV:
      switch(em_stab->stabilization_field_var) {
        case EM_E1_REAL:
          for ( int p=0; p<VIM; p++) {
            for ( int q=0; q<VIM; q++) {
              grad_stabilization_field[p][q] = fv->grad_em_er[p][q];
            }
          }
          break;
        case EM_E1_IMAG:
          for ( int p=0; p<VIM; p++) {
            for ( int q=0; q<VIM; q++) {
              grad_stabilization_field[p][q] = fv->grad_em_ei[p][q];
            }
          }
          break;
        case EM_H1_REAL:
          for ( int p=0; p<VIM; p++) {
            for ( int q=0; q<VIM; q++) {
              grad_stabilization_field[p][q] = fv->grad_em_hr[p][q];
            }
          }
          break;
        case EM_H1_IMAG:
          for ( int p=0; p<VIM; p++) {
            for ( int q=0; q<VIM; q++) {
              grad_stabilization_field[p][q] = fv->grad_em_hi[p][q];
            }
          }
          break;
        default:
          EH(-1,"Cannot have unset stabilization_field_var");
          break;
      }
      break;

    case EM_STAB_PHI_DIVSQUARED:
    case EM_STAB_DPHI_DIVSQUARED:
      switch(em_stab->stabilization_field_var){
        case EM_H1_REAL:
        case EM_H1_IMAG:
          for ( int p=0; p<VIM; p++) {
            for ( int q=0; q<VIM; q++) {
              grad_stabilization_field[p][q] = fv->grad_em_er[p][q] + _Complex_I*fv->grad_em_ei[p][q];
            }
          }
          break;
        case EM_E1_REAL:
        case EM_E1_IMAG:
          for ( int p=0; p<VIM; p++) {
            for ( int q=0; q<VIM; q++) {
              grad_stabilization_field[p][q] = fv->grad_em_hr[p][q] + _Complex_I*fv->grad_em_hi[p][q];
            }
          }
          break;
        default:
          EH(-1,"Must use the first EQN_VAR for a given 3-valued field, e.g. EM_E1_REAL");
      }
      break;

    default:
      //WH(-1,"Cannot use calc_emwave_stabilization without defining type in the struct");
      //return -1;
      break;

  }
  int cartesian_index;
  double complex div_stabilization_field;
  double complex div_stabilization_field_squared;
  switch(em_stab->type) {
    case EM_STAB_PHI_DIV:

      div_stabilization_field = 0.0;

      for (int p=0; p<VIM; p++) {
        div_stabilization_field += grad_stabilization_field[p][p];
      }


      break;

    case EM_STAB_DPHI_DIV:

      div_stabilization_field = 0.0;

      // need the index that corresponds to the x, y or z
      // of the current residual
      cartesian_index = (em_stab->em_eqn - R_EM_E1_REAL)%3;

      for (int p=0; p<VIM; p++) {
        div_stabilization_field += grad_stabilization_field[p][p];
      }

      for( int i=0; i<ei->dof[em_stab->em_eqn]; i++){
        em_stab->residual_term[i]
            = bf[em_stab->em_eqn]->grad_phi[i][cartesian_index]
              *stabilization_coefficient
              *creal(div_stabilization_field);

        for (int b=0; b<DIM; b++){
          for (int j=0; j<ei->dof[em_stab->em_var]; j++){
            em_stab->jacobian_term[i][b][j]
                = bf[em_stab->em_eqn]->grad_phi[i][cartesian_index]
                  *stabilization_coefficient
                  *bf[em_stab->stabilization_field_var + b]->grad_phi[j][b];
          }
        }
      }
      break;

     case EM_STAB_DIVPHI_DIV:
      div_stabilization_field = 0.0;

      double div_phi[MDE] = {0.0};

      for (int p=0; p<VIM; p++) {
        div_stabilization_field += grad_stabilization_field[p][p];
      }
      for (int i=0; i<ei->dof[em_stab->em_eqn]; i++) {
        for (int p=0; p<DIM; p++){
          div_phi[i] += bf[em_stab->em_eqn]->grad_phi[i][p];
        }
      }

      for( int i=0; i<ei->dof[em_stab->em_eqn]; i++){
        em_stab->residual_term[i]
            = div_phi[i]
              *stabilization_coefficient
              *creal(div_stabilization_field);

        for (int b=0; b<DIM; b++){
          for (int j=0; j<ei->dof[em_stab->em_var]; j++){
            em_stab->jacobian_term[i][b][j]
                = div_phi[i]
                  *stabilization_coefficient
                  *bf[em_stab->stabilization_field_var + b]->grad_phi[j][b];
          }
        }
      }
     break;

    case EM_STAB_PHI_DIVSQUARED:
      div_stabilization_field = 0.0;

      for (int p=0; p<VIM; p++) {
        div_stabilization_field += grad_stabilization_field[p][p];
      }
      div_stabilization_field_squared = div_stabilization_field
                                        *div_stabilization_field;
      switch(em_stab->stabilization_field_var) {
        case EM_E1_REAL:
          for( int i=0; i<ei->dof[em_stab->em_eqn]; i++){
            em_stab->residual_term[i]
                = bf[em_stab->em_eqn]->phi[i]
                  *stabilization_coefficient
                  *creal(div_stabilization_field_squared);
                  //*mp->permittivity;
            for (int b=0; b<DIM; b++){
              for (int j=0; j<ei->dof[em_stab->em_var]; j++){
                em_stab->jacobian_term[i][b][j]
                    = bf[em_stab->em_eqn]->phi[i]
                      *stabilization_coefficient
                      *2.0
                      *div_stabilization_field
                      *bf[em_stab->stabilization_field_var + b]->grad_phi[j][b];
              }
            }
          }
        break;
        case EM_H1_REAL:
          for( int i=0; i<ei->dof[em_stab->em_eqn]; i++){
            em_stab->residual_term[i]
                = bf[em_stab->em_eqn]->phi[i]
                  *stabilization_coefficient
                  *creal(div_stabilization_field_squared);
                  //*mag_permeability;
            for (int b=0; b<DIM; b++){
              for (int j=0; j<ei->dof[em_stab->em_var]; j++){
                em_stab->jacobian_term[i][b][j]
                    = bf[em_stab->em_eqn]->phi[i]
                      *stabilization_coefficient
                      *2.0
                      *div_stabilization_field
                      *bf[em_stab->stabilization_field_var + b]->grad_phi[j][b];
              }
            }
          }
          break;
        case EM_E1_IMAG:
          for( int i=0; i<ei->dof[em_stab->em_eqn]; i++){
            em_stab->residual_term[i]
                = bf[em_stab->em_eqn]->phi[i]
                  *stabilization_coefficient
                  *cimag(div_stabilization_field_squared);
                  //*mp->permittivity;
            for (int b=0; b<DIM; b++){
              for (int j=0; j<ei->dof[em_stab->em_var]; j++){
                em_stab->jacobian_term[i][b][j]
                    = bf[em_stab->em_eqn]->phi[i]
                      *stabilization_coefficient
                      *2.0
                      *div_stabilization_field
                      *bf[em_stab->stabilization_field_var + b]->grad_phi[j][b];
              }
            }
          }
        break;
        case EM_H1_IMAG:
          for( int i=0; i<ei->dof[em_stab->em_eqn]; i++){
            em_stab->residual_term[i]
                = bf[em_stab->em_eqn]->phi[i]
                  *stabilization_coefficient
                  *cimag(div_stabilization_field_squared);
                  //*mag_permeability;
            for (int b=0; b<DIM; b++){
              for (int j=0; j<ei->dof[em_stab->em_var]; j++){
                em_stab->jacobian_term[i][b][j]
                    = bf[em_stab->em_eqn]->phi[i]
                      *stabilization_coefficient
                      *2.0
                      *div_stabilization_field
                      *bf[em_stab->stabilization_field_var + b]->grad_phi[j][b];
              }
            }
          }
        break;
        default:
          EH(-1,"must set the stabilization field var with type==phi_divsquared");
          return;
          break;
      }
      break;

    case EM_STAB_DPHI_DIVSQUARED:
      div_stabilization_field = 0.0;

      // need the index that corresponds to the x, y or z
      // of the current residual
      cartesian_index = (em_stab->em_eqn - R_EM_E1_REAL)%3;

      for (int p=0; p<VIM; p++) {
        div_stabilization_field += grad_stabilization_field[p][p];
      }
      div_stabilization_field_squared = div_stabilization_field
                                        *div_stabilization_field;
      switch(em_stab->stabilization_field_var) {
        case EM_E1_REAL:
          for( int i=0; i<ei->dof[em_stab->em_eqn]; i++){
            em_stab->residual_term[i]
                = bf[em_stab->em_eqn]->grad_phi[i][cartesian_index]
                  *stabilization_coefficient
                  *creal(div_stabilization_field_squared)
                  /mag_permeability;
              //*mp->permittivity;
            for (int b=0; b<DIM; b++){
              for (int j=0; j<ei->dof[em_stab->em_var]; j++){
                em_stab->jacobian_term[i][b][j]
                    = bf[em_stab->em_eqn]->phi[i]
                      *stabilization_coefficient
                      *2.0
                      *div_stabilization_field
                      *bf[em_stab->stabilization_field_var + b]->grad_phi[j][b];
              }
            }
          }
        break;
        case EM_H1_REAL:
          for( int i=0; i<ei->dof[em_stab->em_eqn]; i++){
            em_stab->residual_term[i]
                = bf[em_stab->em_eqn]->grad_phi[i][cartesian_index]
                  *stabilization_coefficient
                  *creal(div_stabilization_field_squared);
                  //*mag_permeability;
            for (int b=0; b<DIM; b++){
              for (int j=0; j<ei->dof[em_stab->em_var]; j++){
                em_stab->jacobian_term[i][b][j]
                    = bf[em_stab->em_eqn]->phi[i]
                      *stabilization_coefficient
                      *2.0
                      *div_stabilization_field
                      *bf[em_stab->stabilization_field_var + b]->grad_phi[j][b];
              }
            }
          }
        break;
        case EM_E1_IMAG:
          for( int i=0; i<ei->dof[em_stab->em_eqn]; i++){
            em_stab->residual_term[i]
                = bf[em_stab->em_eqn]->grad_phi[i][cartesian_index]
                  *stabilization_coefficient
                  *cimag(div_stabilization_field_squared)
                  /mag_permeability;
                  //*mp->permittivity;
            for (int b=0; b<DIM; b++){
              for (int j=0; j<ei->dof[em_stab->em_var]; j++){
                em_stab->jacobian_term[i][b][j]
                    = bf[em_stab->em_eqn]->phi[i]
                      *stabilization_coefficient
                      *2.0
                      *div_stabilization_field
                      *bf[em_stab->stabilization_field_var + b]->grad_phi[j][b];
              }
            }
          }
        break;
        case EM_H1_IMAG:
          for( int i=0; i<ei->dof[em_stab->em_eqn]; i++){
            em_stab->residual_term[i]
                = bf[em_stab->em_eqn]->grad_phi[i][cartesian_index]
                  *stabilization_coefficient
                  *cimag(div_stabilization_field_squared);
                  //*mag_permeability;
            for (int b=0; b<DIM; b++){
              for (int j=0; j<ei->dof[em_stab->em_var]; j++){
                em_stab->jacobian_term[i][b][j]
                    = bf[em_stab->em_eqn]->phi[i]
                      *stabilization_coefficient
                      *2.0
                      *div_stabilization_field
                      *bf[em_stab->stabilization_field_var + b]->grad_phi[j][b];
              }
            }
          }
        break;
        default:
          EH(-1,"must set the stabilization field var with type==phi_divsquared");
          return;
          break;
      }
      break;
    case EM_STAB_NONE:
    default:
      for( int i=0; i<ei->dof[em_stab->em_eqn]; i++){
        em_stab->residual_term[i] = 0.0;
        for (int b=0; b<DIM; b++){
          for (int j=0; j<ei->dof[em_stab->em_var]; j++){
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
void
complex_cross_vectors(const complex *v0, /* v0 */
                      const complex *v1, /* v1 */
                      complex *v2) /* v2 = v0 x v1 */
{
  int i, j, k;

  memset(v2, 0, DIM * sizeof(complex));
  for(i = 0; i < DIM; i++)
    for(j = 0; j < DIM; j++)
      for(k = 0; k < DIM; k++)
        v2[k] += permute(i,j,k) * v0[i] * v1[j];
} // end of complex_cross_vectors


#undef I


