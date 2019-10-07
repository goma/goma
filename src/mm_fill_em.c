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
                  const int em_conjvar )
{
  int eqn, var, peqn, pvar, dim, p, q, b, w, i, j, status;
  int dir=0;			/* identity of conjugate variable  */

  dbl EMF = 0, EMF_conj = 0;		/* acoustic pressure	*/
  dbl omega, emf_coeff=0, conj_coeff=0;
  dbl emf_coeff_dn=0, conj_coeff_dn=0;
  dbl emf_coeff_dk=0, conj_coeff_dk=0;
  dbl mag_permeability=1.4e-07;
  int cross_field_var;
  dbl cross_field[DIM];

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

  cpx_refractive_index = n + I*k;
  cpx_rel_permittivity = SQUARE(cpx_refractive_index);
  cpx_permittivity = cpx_rel_permittivity*mp->permittivity;

  //impedance = csqrt(mag_permeability/rel_permittivity);

  switch(em_var)
   {
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
  switch(em_var)
   {
    case EM_E1_REAL:
    case EM_E2_REAL:
    case EM_E3_REAL:
         emf_coeff = omega*cimag(cpx_permittivity);
         conj_coeff = omega*creal(cpx_permittivity);
         emf_coeff_dn = omega*cimag(cpx_permittivity)/n;
         emf_coeff_dk = omega*cimag(cpx_permittivity)/k;
         conj_coeff_dn = omega*creal(cpx_permittivity)/n;
         conj_coeff_dk = omega*creal(cpx_permittivity)/k;
         cross_field_var = EM_H1_REAL;
         for ( p=0; p<VIM; p++)   {cross_field[p] = fv->em_hr[p];}
         break;
    case EM_E1_IMAG:
    case EM_E2_IMAG:
    case EM_E3_IMAG:
         emf_coeff = omega*cimag(cpx_permittivity);
         conj_coeff = -omega*creal(cpx_permittivity);
         emf_coeff_dn = omega*cimag(cpx_permittivity)/n;
         emf_coeff_dk = omega*cimag(cpx_permittivity)/k;
         conj_coeff_dn = -omega*creal(cpx_permittivity)/n;
         conj_coeff_dk = omega*creal(cpx_permittivity)/k;
         cross_field_var = EM_H1_IMAG;
         for ( p=0; p<VIM; p++)   {cross_field[p] = fv->em_hi[p];}
         break;
    case EM_H1_REAL:
    case EM_H2_REAL:
    case EM_H3_REAL:
         emf_coeff = 0;
         conj_coeff = -omega*mag_permeability;
         emf_coeff_dn = 0; emf_coeff_dk = 0; conj_coeff_dn = 0; conj_coeff_dk = 0;
         cross_field_var = EM_E1_REAL;
         for ( p=0; p<VIM; p++)   {cross_field[p] = fv->em_er[p];}
         break;
    case EM_H1_IMAG:
    case EM_H2_IMAG:
    case EM_H3_IMAG:
         emf_coeff = 0;
         conj_coeff = omega*mag_permeability;
         emf_coeff_dn = 0; emf_coeff_dk = 0; conj_coeff_dn = 0; conj_coeff_dk = 0;
         cross_field_var = EM_E1_IMAG;
         for ( p=0; p<VIM; p++)   {cross_field[p] = fv->em_ei[p];}
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

              diffusion *= det_J * wt;
              diffusion *= h3;
              diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
            }

          lec->R[peqn][i] += advection +  diffusion;


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

          lec->J[peqn][pvar][i][j] += advection;
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

          lec->J[peqn][pvar][i][j] += advection;
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

              diffusion *= det_J * wt;
              diffusion *= h3;
              diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];

              lec->J[peqn][pvar][i][j] += diffusion;
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

                  lec->J[peqn][pvar][i][j] += advection;
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


                      lec->J[peqn][pvar][i][j] += advection + diffusion;
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

                      lec->J[peqn][MAX_PROB_VAR + w][i][j] += advection;
                    }
                }
            }

        }
    }

  return(status);
} /* end of assemble_em */
