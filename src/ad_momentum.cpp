#define GOMA_AD_MOMENTUM_CPP
#include <Sacado.hpp>

#include "ad_momentum.h"
#include "ad_stress.h"
#include "ad_turbulence.h"

extern "C" {

/* GOMA include files */
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
ADType ad_ls_modulate_property(
    const ADType &p1, const ADType &p2, double width, double pm_minus, double pm_plus) {
  ADType p_plus, p_minus, p;

  p_minus = p1 * pm_plus + p2 * pm_minus;
  p_plus = p1 * pm_minus + p2 * pm_plus;

  /* Fetch the level set interfacial functions. */
  load_lsi(width);

  /* Calculate the material property. */
  if (ls->Elem_Sign == -1)
    p = p_minus;
  else if (ls->Elem_Sign == 1)
    p = p_plus;
  else
    p = p_minus + (p_plus - p_minus) * lsi->H;

  return (p);
}

int ad_ls_modulate_viscosity(
    ADType &mu1, double mu2, double width, double pm_minus, double pm_plus, const int model) {

  if (model == RATIO) {
    GOMA_EH(GOMA_ERROR, "Invalid Viscosity Model ls_modulate");
  }
  mu1 = ad_ls_modulate_property(mu1, mu2, width, pm_minus, pm_plus);

  return (1);
}

ADType ad_bingham_viscosity(struct Generalized_Newtonian *gn_local,
                            ADType gamma_dot[DIM][DIM]) { /* strain rate tensor */

  ADType gammadot; /* strain rate invariant */

  ADType val1;
  ADType visc_cy;
  ADType yield, shear;
  ADType mu = 0.;
  dbl mu0;
  dbl muinf;
  dbl nexp;
  dbl atexp;
  dbl aexp;
  dbl at_shift;
  dbl lambda;
  dbl tau_y = 0.0;
  dbl fexp;
  dbl temp;
#if MELTING_BINGHAM
  dbl tmelt;
#endif

  ad_calc_shearrate(gammadot, gamma_dot);

  mu0 = gn_local->mu0;
  nexp = gn_local->nexp;
  muinf = gn_local->muinf;
  aexp = gn_local->aexp;
  atexp = gn_local->atexp;
  lambda = gn_local->lam;
  if (gn_local->tau_yModel == CONSTANT) {
    tau_y = gn_local->tau_y;
  } else {
    GOMA_EH(GOMA_ERROR, "Invalid Yield Stress Model");
  }
  fexp = gn_local->fexp;

  if (pd->gv[TEMPERATURE]) {
    temp = fv->T;
  } else {
    temp = upd->Process_Temperature;
  }

  if (DOUBLE_NONZERO(temp) && DOUBLE_NONZERO(mp->reference[TEMPERATURE])) {
    /* normal, non-melting version */
    at_shift = exp(atexp * (1. / temp - 1. / mp->reference[TEMPERATURE]));
    if (!isfinite(at_shift)) {
      at_shift = DBL_MAX;
    }
  } else {
    at_shift = 1.;
  }

  if (DOUBLE_NONZERO(at_shift * lambda * gammadot)) {
    shear = std::pow(at_shift * lambda * gammadot, aexp);
    val1 = std::pow(at_shift * lambda * gammadot, aexp - 1.);
  } else {
    shear = 0.;
  }

  if (DOUBLE_NONZERO(gammadot) && DOUBLE_NONZERO(at_shift)) {
    yield = tau_y * (1. - exp(-at_shift * fexp * gammadot)) / (at_shift * gammadot);
  } else {
    yield = tau_y * fexp;
  }

  visc_cy = pow(1. + shear, (nexp - 1.) / aexp);

  mu = at_shift * (muinf + (mu0 - muinf + yield) * visc_cy);

  return (mu);
}

extern "C" dbl ad_viscosity_wrap(struct Generalized_Newtonian *gn_local) {
  ADType gamma[DIM][DIM];
  for (int i = 0; i < DIM; i++) {
    for (int j = 0; j < DIM; j++) {
      gamma[i][j] = ad_fv->grad_v[i][j] + ad_fv->grad_v[j][i];
    }
  }

  auto mu = ad_viscosity(gn, gamma);
  return mu.val();
}
ADType ad_viscosity(struct Generalized_Newtonian *gn_local, ADType gamma_dot[DIM][DIM]) {
  int err;
  ADType mu = 0.;

  /* this section is for all Newtonian models */
  if (gn_local->ConstitutiveEquation == NEWTONIAN) {
    if (mp->ViscosityModel == CONSTANT) {
      /*  mu   = gn_local->mu0; corrected for auto continuation 3/01 */
      if (gn_local->ConstitutiveEquation == CONSTANT) {
        mu = gn_local->mu0;
      } else {
        mu = mp->viscosity;
      }
      mp_old->viscosity = mu.val();
    } else {
      GOMA_EH(GOMA_ERROR, "Unrecognized viscosity model for Newtonian fluid");
    }
  } else if (gn_local->ConstitutiveEquation == CONSTANT) {
    mu = gn_local->mu0;
    mp_old->viscosity = mu.val();
    /*Sensitivities were already set to zero */
  } else if (gn_local->ConstitutiveEquation == TURBULENT_K_OMEGA) {
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
    ADType F1, F2;
    compute_sst_blending(F1, F2);
    mu = sst_viscosity(Omega, F2);
    // mu = ad_only_turb_k_omega_viscosity();
  } else if (gn_local->ConstitutiveEquation == TURBULENT_SA ||
             gn_local->ConstitutiveEquation == TURBULENT_SA_DYNAMIC) {
    mu = ad_sa_viscosity(gn_local);
  } else if (gn_local->ConstitutiveEquation == BINGHAM) {
    mu = ad_bingham_viscosity(gn_local, gamma_dot);
  } else {
    GOMA_EH(GOMA_ERROR, "Unrecognized viscosity model for non-Newtonian fluid");
  }

  if (ls != NULL && gn_local->ConstitutiveEquation != VE_LEVEL_SET &&
      mp->ViscosityModel != LEVEL_SET && mp->ViscosityModel != LS_QUADRATIC && mp->mp2nd != NULL &&
      (mp->mp2nd->ViscosityModel == CONSTANT || mp->mp2nd->ViscosityModel == RATIO)) {
    err = ad_ls_modulate_viscosity(mu, mp->mp2nd->viscosity, ls->Length_Scale,
                                   (double)mp->mp2nd->viscositymask[0],
                                   (double)mp->mp2nd->viscositymask[1], mp->mp2nd->ViscosityModel);
    GOMA_EH(err, "ls_modulate_viscosity");
  }
  return (mu);
}
/* ad_assemble_momentum -- assemble terms (Residual &| Jacobian) for momentum eqns
 *
 * in:
 * 	ei -- pointer to Element Indeces	structure
 *	pd -- pointer to Problem Description	structure
 *	af -- pointer to Action Flag		structure
 *	bf -- pointer to Basis Function		structure
 *	fv -- pointer to Field Variable		structure
 *	cr -- pointer to Constitutive Relation	structure
 *	md -- pointer to Mesh Derivative	structure
 *	me -- pointer to Material Entity	structure
 * 	ija -- vector of pointers into the a matrix
 * 	a  -- global Jacobian matrix
 * 	R  -- global residual vector
 *
 * out:
 *	a   -- gets loaded up with proper contribution
 *	lec -- gets loaded up with local contributions to resid, Jacobian
 * 	r  -- residual RHS vector
 *
 * Created:	Wed Dec  8 14:03:06 MST 1993 pasacki@sandia.gov
 *
 * Revised:	Sun Feb 27 06:53:12 MST 1994 pasacki@sandia.gov
 *
 *
 * Note: currently we do a "double load" into the addresses in the global
 *       "a" matrix and resid vector held in "esp", as well as into the
 *       local accumulators in "lec".
 *
 */

int ad_assemble_momentum(dbl time,       /* current time */
                         dbl tt,         /* parameter to vary time integration from
                                            explicit (tt = 1) to implicit (tt = 0) */
                         dbl dt,         /* current time step size */
                         dbl h_elem_avg, /* average global element size for PSPG*/
                         const PG_DATA *pg_data,
                         double xi[DIM], /* Local stu coordinates */
                         const Exo_DB *exo) {

  int dim;
  int i, p, q, a;

  int ledof, eqn, ii, peqn;

  int *pde = pd->e[pg->imtrx];
  int *pdv = pd->v[pg->imtrx];

  int status;
  struct Basis_Functions *bfm;

  dbl h3; /* Volume element (scale factors). */

  /* field variables */

  dbl rho; /* Density. */

  ADType f[DIM]; /* Body force. */

  ADType d_area;

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
   * Galerkin weighting functions for i-th and a-th momentum residuals
   * and some of their derivatives...
   */
  dbl phi_i;
  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  ADType Pi[DIM][DIM];

  dbl wt;

  int transient_run = pd->TimeIntegration != STEADY;
  int mass_on;
  int advection_on = 0;
  int source_on = 0;
  int diffusion_on = 0;

  dbl mass_etm, advection_etm, diffusion_etm, source_etm;

  double *R;

  /*
   * Petrov-Galerkin weighting functions for i-th residuals
   * and some of their derivatives...
   */

  int *n_dof = NULL;

  status = 0;

  /*
   * Unpack variables from structures for local convenience...
   */

  eqn = R_MOMENTUM1;

  /*
   * Bail out fast if there's nothing to do...
   */
  if (!pd->e[pg->imtrx][eqn]) {
    return (status);
  }

  dim = pd->Num_Dim;

  wt = fv->wt;

  h3 = fv->h3; /* Differential volume element (scales). */

  d_area = ad_fv->detJ * wt * h3;

  dbl supg = 0.;

  if (mp->Mwt_funcModel == GALERKIN) {
    supg = 0.;
  } else if (mp->Mwt_funcModel == SUPG || mp->Mwt_funcModel == SUPG_GP ||
             mp->Mwt_funcModel == SUPG_SHAKIB) {
    supg = mp->Mwt_func;
  }

  /*** Density ***/
  rho = density(NULL, time);

  ADType tau;
  if (supg != 0.) {
    ad_only_tau_momentum_shakib(tau, dim, dt, FALSE);
  }
  /* end Petrov-Galerkin addition */

  /*
   * Material property constants, etc. Any variations for this
   * Gauss point were evaluated in load_material_properties().
   */

  if (pd->e[pg->imtrx][eqn] & T_POROUS_BRINK) {
    GOMA_EH(GOMA_ERROR, "Porous Brinkmann term not enabled for autodiff momentum");
  }

  /*   eqn = R_MOMENTUM1; */
  /*
   * Field variables...
   */

  /*
   * Calculate the momentum stress tensor at the current gauss point
   */
  ad_fluid_stress(Pi);

  ad_momentum_source_term(f, time);

  /*
   * Residuals_________________________________________________________________
   */
  std::vector<std::vector<ADType>> resid(WIM);
  for (a = 0; a < WIM; a++) {
    resid[a].resize(ei[pg->imtrx]->dof[eqn + a]);
    for (i = 0; i < ei[pg->imtrx]->dof[eqn + a]; i++) {
      resid[a][i] = 0;
    }
  }

  if (af->Assemble_Residual) {
    /*
     * Assemble each component "a" of the momentum equation...
     */
    for (a = 0; a < WIM; a++) {
      eqn = R_MOMENTUM1 + a;
      peqn = upd->ep[pg->imtrx][eqn];
      bfm = bf[eqn];

      mass_on = pde[eqn] & T_MASS;
      advection_on = pde[eqn] & T_ADVECTION;
      diffusion_on = pde[eqn] & T_DIFFUSION;
      source_on = pde[eqn] & T_SOURCE;

      mass_etm = pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
      advection_etm = pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
      diffusion_etm = pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
      source_etm = pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

      /*
       * In the element, there will be contributions to this many equations
       * based on the number of degrees of freedom...
       */

      R = &(lec->R[LEC_R_INDEX(peqn, 0)]);

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];
        if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
          /*
           *  Here is where we figure out whether the row is to placed in
           *  the normal spot (e.g., ii = i), or whether a boundary condition
           *  require that the volumetric contribution be stuck in another
           *  ldof pertaining to the same variable type.
           */
          ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

          phi_i = bf[eqn]->phi[i];
          /* only use Petrov Galerkin on advective term - if required */
          ADType wt_func = phi_i;

          /* add Petrov-Galerkin terms as necessary */
          ADType wt_func_supg = phi_i;
          if (supg != 0.) {
            for (p = 0; p < dim; p++) {
              wt_func += supg * tau * ad_fv->v[p] * bfm->grad_phi[i][p];
            }
          }

          /* this is an optimization for xfem */
          if (xfem != NULL) {
            GOMA_EH(GOMA_ERROR, "xfem not configured for AD Momentum");
          }

          ADType mass = 0.;
          if (transient_run) {
            if (mass_on) {
              mass = ad_fv->v_dot[a] * rho;
              mass *= -wt_func * d_area;
              mass *= mass_etm;
            }
          }

          ADType advection = 0.;
          if (advection_on) {
#ifdef DO_NO_UNROLL
            for (p = 0; p < WIM; p++) {
              advection += (v[p] - x_dot[p]) * grad_v[p][a];
            }
#else
            advection += (ad_fv->v[0] - ad_fv->x_dot[0]) * ad_fv->grad_v[0][a];
            advection += (ad_fv->v[1] - ad_fv->x_dot[1]) * ad_fv->grad_v[1][a];
            if (WIM == 3)
              advection += (ad_fv->v[2] - ad_fv->x_dot[2]) * ad_fv->grad_v[2][a];
#endif

            advection *= rho;
            advection *= -wt_func * d_area;
            advection *= advection_etm;
          }

          ADType diffusion = 0.;
          if (diffusion_on) {
            for (p = 0; p < VIM; p++) {
              for (q = 0; q < VIM; q++) {
                diffusion += ad_fv->basis[eqn].grad_phi_e[i][a][p][q] * Pi[q][p];
              }
            }
            diffusion *= -d_area;
            diffusion *= diffusion_etm;
          }

          ADType graddiv = 0;
          if (0) {
            ADType gamma[DIM][DIM];
            for (int i = 0; i < DIM; i++) {
              for (int j = 0; j < DIM; j++) {
                gamma[i][j] = ad_fv->grad_v[i][j] + ad_fv->grad_v[j][i];
              }
            }
            ADType mu = ad_viscosity(gn, gamma);
            ADType div_grad_phi = 0;
            ADType div_v = 0;
            for (int w = 0; w < dim; w++) {
              div_grad_phi += ad_fv->basis[eqn].grad_phi_e[i][a][w][w] *
                              ad_fv->basis[eqn].grad_phi_e[i][a][w][w];
              div_v += ad_fv->grad_v[w][w] * ad_fv->grad_v[w][w];
            }
            graddiv = (2.0 / 3.0) * mu * div_grad_phi * mu * div_v;
            graddiv *= -d_area * h3 * wt;
          }

          /*
           * Source term...
           */
          ADType source = 0.0;
          if (source_on) {
            source += f[a];
            source *= wt_func * d_area;
            source *= source_etm;
          }

          /*
           * Add contributions to this residual (globally into Resid, and
           * locally into an accumulator)
           */

          /*lec->R[LEC_R_INDEX(peqn,ii)] += mass + advection + porous + diffusion + source;*/
          R[ii] += mass.val() + advection.val() + diffusion.val() + source.val() + graddiv.val();
          resid[a][ii] += mass + advection + diffusion + source + graddiv;
          // if (ei[pg->imtrx]->ielem == 418) {
          //   printf("diff = %.15f\n", diffusion.val());
          // }
        } /*end if (active_dofs) */
      }   /* end of for (i=0,ei[pg->imtrx]->dofs...) */
    }
  }

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    for (a = 0; a < WIM; a++) {
      eqn = R_MOMENTUM1 + a;
      peqn = upd->ep[pg->imtrx][eqn];

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];
        if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
          ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

          for (int var = V_FIRST; var < V_LAST; var++) {

            /* Sensitivity w.r.t. velocity */
            if (pdv[var]) {
              int pvar = upd->vp[pg->imtrx][var];

              for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                // J = &(lec->J[LEC_J_INDEX(peqn, pvar, ii, 0)]);
                lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += resid[a][ii].dx(ad_fv->offset[var] + j);

              } /* End of loop over j */
            }   /* End of if the variale is active */
          }

        } /* end of if(active_dofs) */
      }   /* End of loop over i */
    }     /* End of if assemble Jacobian */
  }
  safe_free((void *)n_dof);
  return (status);
}

void ad_ve_polymer_stress(ADType gamma[DIM][DIM], ADType stress[DIM][DIM]) {
#if 1

  dbl dgamma[DIM][DIM];
  for (int i = 0; i < VIM; i++) {
    for (int j = 0; j < VIM; j++) {
      stress[i][j] = 0;
      dgamma[i][j] = gamma[i][j].val();
    }
  }
  switch (vn->evssModel) {
  case SQRT_CONF: {
    for (int mode = 0; mode < vn->modes; mode++) {
      /* get polymer viscosity */
      ADType mup = ad_viscosity(ve[mode]->gn, gamma);
      dbl lambda = polymer_time_const(ve[mode]->time_const_st, dgamma, NULL);

      ADType bdotb[DIM][DIM];
      ADType b[DIM][DIM];
      for (int ii = 0; ii < VIM; ii++) {
        for (int jj = 0; jj < VIM; jj++) {
          if (ii <= jj) {
            b[ii][jj] = ad_fv->S[mode][ii][jj];
            b[jj][ii] = b[ii][jj];
          }
        }
      }

      ad_tensor_dot(b, b, bdotb, VIM);

      for (int ii = 0; ii < VIM; ii++) {
        for (int jj = 0; jj < VIM; jj++) {
          stress[ii][jj] += -(mup / lambda) * (delta(ii, jj) - bdotb[ii][jj]);
        }
      }
    }

  } break;
  default: // Regular stress formulations
  {
    for (int mode = 0; mode < vn->modes; mode++) {
      for (int i = 0; i < VIM; i++) {
        for (int j = 0; j < VIM; j++) {
          stress[i][j] += ad_fv->S[mode][i][j];
        }
      }
    }
  } break;
  }
#endif
}

/*
 * Calculate the total stress tensor for a fluid at a single gauss point
 *  This includes the diagonal pressure contribution
 *
 *  Pi = stress tensor
 *  d_Pi = dependence of the stress tensor on the independent variables
 */
void ad_fluid_stress(ADType Pi[DIM][DIM]) {

  /*
   * Variables for vicosity and derivative
   */
  ADType mu = 0.0;

  /* polymer viscosity and derivatives */

  ADType mup;

  /*  shift function */
  ADType at = 1.0;

  /* solvent viscosity and derivatives */

  ADType mus = 0.0;

  /* numerical "adaptive" viscosity and derivatives */

  ADType mu_num;
  if (pd->gv[TEMPERATURE]) {
    GOMA_EH(GOMA_ERROR, "Temperature not yet implemented ad_fluid_stress");
  }

  /*
   * Field variables...
   */

  dbl evss_f = 0;
  if (pd->gv[POLYMER_STRESS11] && is_evss_f_model(vn->evssModel)) {
    evss_f = 1.0;
  }

  double Heaviside = 1;
  if (ls != NULL && ls->ghost_stress) {
    load_lsi(ls->Length_Scale);
    switch (ls->ghost_stress) {
    case LS_OFF:
      Heaviside = 1;
      break;
    case LS_POSITIVE:
      Heaviside = lsi->H;
      break;
    case LS_NEGATIVE:
      Heaviside = 1 - lsi->H;
      break;
    default:
      GOMA_EH(GOMA_ERROR, "Unknown Level Set Ghost Stress value");
      break;
    }
  }

  ADType gamma_cont[DIM][DIM];
  ADType gamma[DIM][DIM];
  if (evss_f) {
    for (int a = 0; a < VIM; a++) {
      for (int b = 0; b < VIM; b++) {
        gamma_cont[a][b] = ad_fv->G[a][b] + ad_fv->G[b][a];
      }
    }
  } else {
    for (int a = 0; a < VIM; a++) {
      for (int b = 0; b < VIM; b++) {
        gamma_cont[a][b] = 0;
      }
    }
  }

  /* load up shear rate tensor based on velocity */
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      gamma[a][b] = ad_fv->grad_v[a][b] + ad_fv->grad_v[b][a];
    }
  }

  // TODO AD
  mu = ad_viscosity(gn, gamma);
  if (pd->gv[POLYMER_STRESS11]) {
    mus = ad_viscosity(gn, gamma);

    /* This is the adaptive viscosity from Sun et al., 1999.
     * The term multiplies the continuous and discontinuous
     * shear-rate, so it should cancel out and not affect the
     * solution, other than increasing the stability of the
     * algorithm in areas of high shear and stress.
     */

    mu_num = 1;
    // if (DOUBLE_NONZERO(vn->eps)) {
    //   for (int mode = 0; mode < vn->modes; mode++) {
    //     for (a = 0; a < VIM; a++) {
    //       for (b = 0; b < VIM; b++) {
    //         s[a][b] += fv->S[mode][a][b];
    //       }
    //     }
    //   }

    //   mu_num = numerical_viscosity(s, gamma_cont, d_mun_dS, d_mun_dG);

    mu = mu_num * mus;

    for (int mode = 0; mode < vn->modes; mode++) {
      /* get polymer viscosity */
      // TODO AD
      mup = ad_viscosity(ve[mode]->gn, gamma);

      mu += Heaviside * mu_num * at * mup;

    } // for mode
  }   // if POLYMER_STRESS

  /*
   * Calculate the dilational viscosity, if necessary
   */
  if (mp->DilationalViscosityModel != DILVISCM_KAPPAWIPESMU) {
    GOMA_EH(GOMA_ERROR, "Dilational viscosity not enabled AD fluid stress");
  }

  /*
   * Viscoelastic Stress contributions
   */
  ADType polymer_stress[DIM][DIM];
  if (pd->gv[POLYMER_STRESS11]) {
    ad_ve_polymer_stress(gamma, polymer_stress);
  }

  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      Pi[a][b] = -ad_fv->P * (double)delta(a, b) + mu * gamma[a][b];
    }
  }

  if (pd->gv[POLYMER_STRESS11]) {
    for (int a = 0; a < VIM; a++) {
      for (int b = 0; b < VIM; b++) {
        // TODO : derivative may be missing here...
        Pi[a][b] += -evss_f * (mu - mus) * gamma_cont[a][b] + Heaviside * polymer_stress[a][b];
      }
    }
  }
  if (gn->ConstitutiveEquation == BINGHAM_MIXED) {
    GOMA_EH(-1, "BINGHAM_MIXED not enabled for AD fluid stress");
  }
}
/******************************************************************************
 * momentum_source_term(): Computes the body force term for the momentum balance.
 *
 * Input
 * -----
 *   f    == Body force
 *   df->T == Derivative w.r.t. temperature
 *   df->X == Derivative w.r.t. mesh displacements
 *   df->C == Derivative w.r.t. concentration
 *   df->v == Derivative w.r.t. velocity
 *   df->F == Derivative w.r.t. FILL
 *   df->E == Derivative w.r.t. electric field
 *
 ******************************************************************************/
int ad_momentum_source_term(ADType f[DIM], /* Body force. */
                            dbl time) {
  int a;
  int eqn;
  const int dim = pd->Num_Dim;
  int status = 0;

  /* initialize everything to zero */
  for (int a = 0; a < VIM; a++) {
    f[a] = 0;
  }

  /****Momentum Source Model******/
  if (mp->MomentumSourceModel == CONSTANT) {
    int force_dim = dim;
    if (pd->CoordinateSystem == CARTESIAN_2pt5D) {
      force_dim = 3;
    }
    for (int a = 0; a < force_dim; a++) {
      eqn = R_MOMENTUM1 + a;
      if (pd->e[upd->matrix_index[eqn]][eqn] & T_SOURCE) {
        f[a] = mp->momentum_source[a];
      }
    }
  } else if (mp->MomentumSourceModel == LEVEL_SET && pd->gv[FILL]) {
    DENSITY_DEPENDENCE_STRUCT d_rho_struct; /* density dependence */
    DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;
    double rho = density(d_rho, time);

    for (a = 0; a < pd->Num_Dim; a++) {
      f[a] = mp->momentum_source[a] * rho;
    }
  } else if (mp->MomentumSourceModel == VARIABLE_DENSITY) {
    if (mp->DensityModel == SOLVENT_POLYMER) {
      double rho = density(NULL, time);
      for (a = 0; a < dim; a++) {
        eqn = R_MOMENTUM1 + a;
        if (pd->e[upd->matrix_index[eqn]][eqn] & T_SOURCE) {
          f[a] = rho * mp->momentum_source[a];
        }
      }
    } else {
      GOMA_EH(GOMA_ERROR, "Unknown density model for variable density");
    }
  } else {
    GOMA_EH(GOMA_ERROR, "Unknown Navier-Stokes Source Model for AD");
  }

  if (ls != NULL && mp->mp2nd != NULL && mp->MomentumSourceModel != LEVEL_SET &&
      mp->mp2nd->MomentumSourceModel == CONSTANT && (pd->e[pg->imtrx][R_MOMENTUM1] & T_SOURCE)) {
    GOMA_EH(GOMA_ERROR, "Unknown Navier-Stokes Source Model for AD");
  }

  return (status);
}

int ad_calc_pspg(ADType pspg[DIM],
                 dbl time_value, /* current time */
                 dbl tt,         /* parameter to vary time integration from
                                                    explicit (tt = 1) to implicit (tt = 0)    */
                 dbl dt,         /* current time step size                    */
                 const PG_DATA *pg_data) {
  const dbl h_elem_avg = pg_data->h_elem_avg;

  int dim;
  int p, a, b;

  int pspg_global;
  int pspg_local;

  dbl *v = fv->v; /* Velocity field. */
  ADType *grad_P = ad_fv->grad_P;

  /*
   * Variables for vicosity and derivative
   */
  ADType mu;

  /*
   * density and sensitivity terms
   */
  dbl rho;
  dbl rho_t;

  ADType f[DIM]; /* Body force. */

  /*
   * Interpolation functions...
   */

  /*
   * Variables for Pressure Stabilization Petrov-Galerkin...
   */
  int meqn, meqn1;

  ADType mass, advection, diffusion, source, porous;

  ADType momentum[DIM]; /* momentum residual for PSPG */
  ADType x_dot[DIM];
  ADType v_dot[DIM];
  ADType *grad_v[DIM];
  ADType div_s[DIM];
  ADType div_G[DIM];

  /* variables for Brinkman porous flow */
  ADType por = 0, por2 = 0, per = 0, vis = 0, dvis_dT[MDE], sc = 0, speed = 0;

  ADType h_elem = 0;
  ADType Re;
  ADType tau_pspg = 0.;
  ADType tau_pspg1 = 0.;

  ADType hh_siz, vv_speed;

  ADType gamma[DIM][DIM]; /* shrearrate tensor based on velocity */
  dbl dgamma[DIM][DIM];   /* shrearrate tensor based on velocity */

  int mode;

  /* For particle momentum model.
   */
  int species; /* species number for particle phase,  */
  dbl ompvf;   /* 1 - partical volume fraction */

  dim = pd->Num_Dim;

  /* initialize */
  for (a = 0; a < DIM; a++)
    pspg[a] = 0.;

  /* This is the flag for the standard global PSPG */
  if (PSPG == 1) {
    pspg_global = TRUE;
    pspg_local = FALSE;
  }
  /* This is the flag for the standard local PSPG */
  else if (PSPG == 2) {
    pspg_global = FALSE;
    pspg_local = TRUE;
  } else if (PSPG == 3) { // Shakib
    pspg_global = FALSE;
    pspg_local = FALSE;
  } else {
    return 0;
  }

  int v_s[MAX_MODES][DIM][DIM];

  if (pd->gv[POLYMER_STRESS11]) {
    stress_eqn_pointer(v_s);
  }

  /* initialize dependencies */

  if (pd->e[pg->imtrx][R_PMOMENTUM1]) {
    species = (int)mp->u_density[0];
    ompvf = 1.0 - fv->c[species];
  } else {
    species = -1;
    ompvf = 1.0;
  }

  // Global average for pspg_global's element size
  h_elem = h_elem_avg;

  /*** Density ***/
  rho = density(NULL, time_value);

  ADType pspg_tau;

  if (pspg_global) {
    GOMA_EH(GOMA_ERROR, "PSPG Global not enabled for AD");
  } else if (pspg_local) {
    GOMA_EH(GOMA_ERROR, "PSPG Local not enabled for AD");
  } else if (PSPG == 3) { // shakib
    ad_only_tau_momentum_shakib(pspg_tau, pd->Num_Dim, dt, TRUE);
  }

  for (a = 0; a < VIM; a++)
    grad_v[a] = ad_fv->grad_v[a];

  /* load up shearrate tensor based on velocity */
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      gamma[a][b] = grad_v[a][b] + grad_v[b][a];
      dgamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
    }
  }

  /*
   * get viscosity for velocity second derivative/diffusion
   * term in PSPG stuff
   */
  mu = ad_viscosity(gn, gamma);

  /* get variables we will need for momentum residual */

  for (a = 0; a < WIM; a++) {
    if (pd->TimeIntegration != STEADY && pd->v[pg->imtrx][MESH_DISPLACEMENT1 + a]) {
      x_dot[a] = ad_fv->x_dot[a];
    } else {
      x_dot[a] = 0.;
    }

    if (pd->TimeIntegration != STEADY) {
      v_dot[a] = ad_fv->v_dot[a];
    } else {
      v_dot[a] = 0.;
    }
  }

  for (p = 0; p < WIM; p++)
    div_s[p] = 0.;

  if (pd->gv[POLYMER_STRESS11]) {
    if (vn->evssModel == LOG_CONF || vn->evssModel == LOG_CONF_GRADV ||
        vn->evssModel == LOG_CONF_TRANSIENT || vn->evssModel == LOG_CONF_TRANSIENT_GRADV) {
      for (mode = 0; mode < vn->modes; mode++) {
        dbl lambda = 0.0;
        lambda = polymer_time_const(ve[mode]->time_const_st, dgamma, NULL);
        dbl mup = viscosity(ve[mode]->gn, dgamma, NULL);
        int dofs = ei[upd->matrix_index[v_s[mode][0][0]]]->dof[v_s[mode][0][0]];
        dbl grad_S[DIM][DIM][DIM] = {{{0.0}}};
        dbl s[MDE][DIM][DIM];
        dbl exp_s[MDE][DIM][DIM] = {{{0.0}}};
        dbl eig_values[DIM];
        dbl R[DIM][DIM];
        for (int k = 0; k < dofs; k++) {
          if (pg->imtrx == upd->matrix_index[POLYMER_STRESS11] &&
              (vn->evssModel == LOG_CONF_TRANSIENT_GRADV || vn->evssModel == LOG_CONF_TRANSIENT)) {
            for (int i = 0; i < VIM; i++) {
              for (int j = 0; j < VIM; j++) {
                if (j >= i) {
                  s[k][i][j] = *esp_old->S[mode][i][j][k];
                } else {
                  s[k][i][j] = *esp_old->S[mode][j][i][k];
                }
              }
            }
          } else {
            for (int i = 0; i < VIM; i++) {
              for (int j = 0; j < VIM; j++) {
                if (j >= i) {
                  s[k][i][j] = *esp->S[mode][i][j][k];
                } else {
                  s[k][i][j] = *esp->S[mode][j][i][k];
                }
              }
            }
          }
          compute_exp_s(s[k], exp_s[k], eig_values, R);
        }
        for (int p = 0; p < VIM; p++) {
          for (int q = 0; q < VIM; q++) {
            for (int r = 0; r < VIM; r++) {
              grad_S[r][p][q] = 0.;
              for (int i = 0; i < dofs; i++) {
                if (p <= q) {
                  grad_S[r][p][q] += exp_s[i][p][q] * bf[POLYMER_STRESS11]->grad_phi[i][r];
                } else {
                  grad_S[r][p][q] += exp_s[i][q][p] * bf[POLYMER_STRESS11]->grad_phi[i][r];
                }
              }
            }
          }
        }
        dbl div_exp_s[DIM];
        for (int r = 0; r < dim; r++) {
          div_exp_s[r] = 0.0;

          for (int q = 0; q < dim; q++) {
            div_exp_s[r] += grad_S[q][q][r];
          }
        }
        for (p = 0; p < WIM; p++) {
          div_s[p] += (mup / lambda) * div_exp_s[p];
        }
      }
    } else if (vn->evssModel == SQRT_CONF) {
      for (mode = 0; mode < vn->modes; mode++) {
        int dofs = ei[upd->matrix_index[v_s[mode][0][0]]]->dof[v_s[mode][0][0]];
        dbl lambda = polymer_time_const(ve[mode]->time_const_st, dgamma, NULL);
        dbl mup = viscosity(ve[mode]->gn, dgamma, NULL);
        ADType b[MDE][DIM][DIM];
        ADType bdotb[MDE][DIM][DIM];
        for (int i = 0; i < VIM; i++) {
          for (int j = 0; j < VIM; j++) {

            for (int k = 0; k < dofs; k++) {
              if (j >= i) {
                b[k][i][j] = ADType(ad_fv->total_ad_variables, ad_fv->offset[v_s[mode][i][j]] + k,
                                    *esp->S[mode][i][j][k]);
              } else {
                b[k][i][j] = ADType(ad_fv->total_ad_variables, ad_fv->offset[v_s[mode][i][j]] + k,
                                    *esp->S[mode][j][i][k]);
              }
            }
          }
        }
        for (int i = 0; i < dofs; i++) {
          ad_tensor_dot(b[i], b[i], bdotb[i], VIM);
        }

        ADType grad_S[DIM][DIM][DIM];
        for (int p = 0; p < VIM; p++) {
          for (int q = 0; q < VIM; q++) {
            for (int r = 0; r < VIM; r++) {
              grad_S[r][p][q] = 0.;
              for (int i = 0; i < dofs; i++) {

                if (p <= q) {
                  grad_S[r][p][q] += bdotb[i][p][q] * ad_fv->basis[v_s[mode][p][q]].grad_phi[i][r];
                } else {
                  grad_S[r][p][q] += bdotb[i][q][p] * ad_fv->basis[v_s[mode][q][p]].grad_phi[i][r];
                }
              }
            }
          }
        }
        ADType div_bdotb[DIM];
        for (int r = 0; r < dim; r++) {
          div_bdotb[r] = 0.0;

          for (int q = 0; q < dim; q++) {
            div_bdotb[r] += grad_S[q][q][r];
          }
        }
        for (p = 0; p < WIM; p++) {
          div_s[p] += (mup / lambda) * div_bdotb[p];
        }
      }
      // for (mode = 0; mode < vn->modes; mode++) {
      //   dbl lambda = 0.0;
      //   if (ve[mode]->time_constModel == CONSTANT) {
      //     lambda = ve[mode]->time_const;
      //   }
      //   dbl mup = viscosity(ve[mode]->gn, dgamma, NULL);
      //   ADType b[DIM][DIM];
      //   ADType divb[DIM];
      //   for (int i = 0; i < VIM; i++) {
      //     divb[i] = ad_fv->div_S[mode][i];
      //     for (int j = 0; j < VIM; j++) {
      //       b[i][j] = ad_fv->S[mode][i][j];
      //     }
      //   }
      //   ADType divbdotb[DIM];
      //   ADType bddotgradb[DIM];
      //   for (int i = 0; i < VIM; i++) {
      //     divbdotb[i] = 0;
      //     bddotgradb[i] = 0;
      //     for (int j = 0; j < VIM; j++) {
      //       divbdotb[i] += divb[j] * b[j][i];
      //       for (int k = 0; k < VIM; k++) {
      //         bddotgradb[i] += b[i][j] * ad_fv->grad_S[mode][i][j][k];
      //       }
      //     }
      //   }

      //   for (int i = 0; i < VIM; i++) {
      //     div_s[i] += -(mup / lambda) * -1.0 * (divbdotb[i] + bddotgradb[i]);
      //   }
      // }
    } else {
      for (p = 0; p < WIM; p++) {
        for (mode = 0; mode < vn->modes; mode++) {
          div_s[p] += fv->div_S[mode][p];
        }
      }
    }
  }

  if (cr->MassFluxModel == DM_SUSPENSION_BALANCE || cr->MassFluxModel == HYDRODYNAMIC_QTENSOR_OLD ||
      cr->MassFluxModel == HYDRODYNAMIC_QTENSOR)
    GOMA_EH(GOMA_ERROR, "Particle stress not enabled for AD");

  if (pd->gv[VELOCITY_GRADIENT11]) {
    for (p = 0; p < WIM; p++) {
      div_G[p] = ad_fv->div_G[p];
    }
  } else {
    for (p = 0; p < WIM; p++) {
      div_G[p] = 0.;
    }
  }

  if (pd->e[upd->matrix_index[R_MOMENTUM1]][R_MOMENTUM1] & T_POROUS_BRINK) {
    GOMA_EH(GOMA_ERROR, "Porous Brinkman not enabled for AD");
  } else {
    por = 1.;
    por2 = 1.;
    per = 1.;
    vis = mp->viscosity;
    sc = 0.;
  }

  /* for porous media stuff */
  speed = 0.0;
  for (a = 0; a < WIM; a++) {
    speed += v[a] * v[a];
  }
  speed = sqrt(speed);

  /* get momentum source term */
  ad_momentum_source_term(f, time_value);

  if (pd->gv[R_PMOMENTUM1]) {
    rho_t = ompvf * rho;
    meqn1 = R_PMOMENTUM1;
  } else {
    rho_t = rho;
    meqn1 = R_MOMENTUM1;
  }

  for (a = 0; a < WIM; a++) {
    meqn = meqn1 + a;

    mass = 0.;
    if ((pd->e[upd->matrix_index[meqn]][meqn] & T_MASS) && (pd->TimeIntegration != STEADY)) {
      mass = rho_t * v_dot[a] / por;
      mass *= pd->etm[upd->matrix_index[meqn]][meqn][(LOG2_MASS)];
    }

    advection = 0.;
    if (pd->e[upd->matrix_index[meqn]][meqn] & T_ADVECTION) {
      for (p = 0; p < WIM; p++) {
        advection += rho_t * (v[p] - x_dot[p]) * grad_v[p][a] / por2;
      }
      advection *= pd->etm[upd->matrix_index[meqn]][meqn][(LOG2_ADVECTION)];
    }

    diffusion = 0.;
    if (pd->e[upd->matrix_index[meqn]][meqn] & T_DIFFUSION) {
      diffusion = grad_P[a] - div_s[a];
      /*diffusion  -= div_tau_p[a]  */
      diffusion -= mu * div_G[a];
      diffusion *= pd->etm[upd->matrix_index[meqn]][meqn][(LOG2_DIFFUSION)];
    }

    source = 0.;
    if (pd->e[upd->matrix_index[meqn]][meqn] & T_SOURCE) {
      source = -f[a];
      source *= pd->etm[upd->matrix_index[meqn]][meqn][(LOG2_SOURCE)];
    }

    porous = 0.;
    if (pd->e[upd->matrix_index[meqn]][meqn] & T_POROUS_BRINK) {
    }

    momentum[a] = mass + advection + diffusion + source + porous;
    pspg[a] = pspg_tau * momentum[a];
  }

  return 0;
}

int ad_assemble_continuity(dbl time_value, /* current time */
                           dbl tt,         /* parameter to vary time integration from
                                              explicit (tt = 1) to implicit (tt = 0)    */
                           dbl dt,         /* current time step size                    */
                           const PG_DATA *pg_data) {
  int a;

  int eqn;
  int peqn, pvar;

  int i, j;
  int status, err;

  dbl time = 0.0; /*  RSL 6/6/02  */

  ADType det_J;
  dbl h3;
  dbl wt;
  ADType d_area;

  /*
   * Galerkin weighting functions...
   */

  dbl phi_i;

  /*
   * Variables for Pressure Stabilization Petrov-Galerkin...
   */
  int meqn;
  int v_s[MAX_MODES][DIM][DIM];

  int *pdv = pd->v[pg->imtrx];

  ADType pspg[DIM];

  struct Species_Conservation_Terms s_terms;
  int w0 = 0;

  /* For particle momentum model.
   */
  int advection_on = 0;
  int ion_reactions_on = 0, electrode_kinetics_on = 0;
  int lagrangian_mesh_motion = 0, total_ale_on = 0;
  int hydromassflux_on = 0, suspensionsource_on = 0;
  int total_ale_and_velo_off = 0;

  dbl advection_etm;

  status = 0;

  /*
   * Unpack variables from structures for local convenience...
   */

  eqn = R_PRESSURE;
  peqn = upd->ep[pg->imtrx][eqn];

  ADType div_v = 0;

  for (a = 0; a < VIM; a++) {
    div_v += ad_fv->grad_v[a][a];
  }

  /*
   * Bail out fast if there's nothing to do...
   */

  if (!pd->e[pg->imtrx][eqn]) {
    return (status);
  }

  if (pd->gv[POLYMER_STRESS11]) {
    err = stress_eqn_pointer(v_s);
  }

  wt = fv->wt;
  det_J = ad_fv->detJ; /* Really, ought to be mesh eqn. */
  h3 = fv->h3;         /* Differential volume element (scales). */

  d_area = wt * det_J * h3;

  /*
   * Get the deformation gradients and tensors if needed
   */

  lagrangian_mesh_motion = (cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN);
  electrode_kinetics_on = (mp->SpeciesSourceModel[0] == ELECTRODE_KINETICS);
  ion_reactions_on = (mp->SpeciesSourceModel[0] == ION_REACTIONS);
  total_ale_on = (cr->MeshMotion == TOTAL_ALE);
  hydromassflux_on = (cr->MassFluxModel == HYDRODYNAMIC);
  suspensionsource_on = (mp->MomentumSourceModel == SUSPENSION);

  if (lagrangian_mesh_motion && pd->gv[R_MESH1]) {
    err = belly_flop(elc->lame_mu);
    GOMA_EH(err, "error in belly flop");
    if (err == 2)
      return (err);
  }

  if (total_ale_on && !pd->gv[VELOCITY1]) {
    total_ale_and_velo_off = 1;
  }

  if (total_ale_and_velo_off && pd->gv[R_SOLID1]) {
    err = belly_flop_rs(elc_rs->lame_mu);
    GOMA_EH(err, "error in belly flop for real solid");
    if (err == 2)
      return (err);
  }

  if (PSPG) {
    ad_calc_pspg(pspg, time_value, tt, dt, pg_data);
  }

  if (mp->MomentumSourceModel == SUSPENSION_PM || electrode_kinetics_on) /*  RSL 7/25/00  */
  {
    err = get_continuous_species_terms(&s_terms, 0.0, tt, dt, pg_data->hsquared);
    GOMA_EH(err, "problem in getting the species terms");
  }

  if (ion_reactions_on) /*  RSL 3/19/01 and 6/6/02  */
  {
    zero_structure(&s_terms, sizeof(struct Species_Conservation_Terms), 1);
    err = get_continuous_species_terms(&s_terms, time, tt, dt, pg_data->hsquared);
    GOMA_EH(err, "problem in getting the species terms");
  }

  if ((hydromassflux_on) && (mp->DensityModel == SUSPENSION) && (suspensionsource_on)) {
    /*
     * Compute hydrodynamic/sedimentation flux and sensitivities.
     */

    w0 =
        (int)mp->u_density[0]; /* This is the species number that is transported HYDRODYNAMICally */

    hydro_flux(&s_terms, w0, tt, dt, pg_data->hsquared);
  }

  advection_on = pd->e[pg->imtrx][eqn] & T_ADVECTION;

  advection_etm = pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

  dbl ls_disable_pspg = 1;

  std::vector<ADType> resid(ei[pg->imtrx]->dof[eqn]);
  for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
    resid[i] = 0;
  }
  if (af->Assemble_Residual) {
    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

      phi_i = bf[eqn]->phi[i];

      /*
       *  Mass Terms: drhodt terms (usually though problem dependent)
       */

      /*
       *  Advection:
       *    This term refers to the standard del dot v .
       *
       *    int (phi_i div_v d_omega)
       *
       *   Note density is not multiplied into this term normally
       */
      ADType advection = 0.0;
      if (advection_on) {
        if (pd->gv[VELOCITY1]) /* then must be solving fluid mechanics in this material */
        {

          /*
           * Standard incompressibility constraint means we have
           * a solenoidal velocity field
           */

          advection = div_v;

          advection *= phi_i * d_area;
          advection *= advection_etm;
        }
      }

      ADType pressure_stabilization = 0.0;
      if (PSPG) {
        for (a = 0; a < WIM; a++) {
          meqn = R_MOMENTUM1 + a;
          if (pd->gv[meqn]) {
            pressure_stabilization += ad_fv->basis[eqn].grad_phi[i][a] * pspg[a];
          }
        }
        pressure_stabilization *= d_area * ls_disable_pspg;
      }

      /*
       *  Add up the individual contributions and sum them into the local element
       *  contribution for the total continuity equation for the ith local unknown
       */
      lec->R[LEC_R_INDEX(peqn, i)] += advection.val() + pressure_stabilization.val();
      resid[i] = advection + pressure_stabilization;
    }
  }
  if (af->Assemble_Jacobian) {
    eqn = PRESSURE;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

      /* Sensitivity w.r.t. velocity */
      for (int var = V_FIRST; var < V_LAST; var++) {
        if (pdv[var]) {
          pvar = upd->vp[pg->imtrx][var];
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += resid[i].dx(ad_fv->offset[var] + j);
          } /* End of loop over j */
        }   /* End of if the variale is active */
      }

    } /* End of loop over i */
  }   /* End of if assemble Jacobian */
  return 0;
}
