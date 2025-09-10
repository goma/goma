#ifndef GOMA_ENABLE_SACADO
extern "C" {
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "std.h"
}

extern "C" dbl
fluidity_viscosity(int fluidity_species, /* integer associated with conc eqn for bond */
                   VISCOSITY_DEPENDENCE_STRUCT *d_mu) {
  GOMA_EH(GOMA_ERROR, "Fluidity viscosity requires Sacado\n");
  return 0.;
}

extern "C" int fluidity_source(int species_no, struct Species_Conservation_Terms *st, dbl *params) {
  GOMA_EH(GOMA_ERROR, "Fluidity source requires Sacado\n");
  return -1;
}
extern "C" void fluidity_equilibrium_surf(double func[DIM],
                                          double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                                          int fluidity_species,
                                          const double theta,
                                          const double dt) {
  GOMA_EH(GOMA_ERROR, "Fluidity equilibrium surface requires Sacado\n");
}
#else // GOMA_ENABLE_SACADO

#include "Sacado.hpp"
#include "models/fluidity.h"

extern "C" {
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "mm_fill_util.h"
#include "mm_mp.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "std.h"
}

using sfad = Sacado::Fad::SFad<double, 10>;

// Function to calculate equilibrium fluidity phi_eq as a function of stress sigma
sfad calculate_phi_eq(sfad sigma, dbl phi_inf, dbl phi_0, dbl K, dbl n) {
  sfad term = pow(sigma / K, 1.0 / n) * (1.0 / sigma);
  sfad phi_eq = term / ((phi_inf - phi_0) + term);
  return phi_eq;
}

// Smoothed Heaviside function using tanh
// m controls the width of the transition
sfad smoothed_heaviside(sfad x, double m) { return 0.5 * (1.0 + tanh(m * x)); }

/**
 * Calculate the fluidity viscosity from a fluidity species. The species is expected to be
 * a normalized fluidity between 0 and 1
 */
extern "C" dbl
fluidity_viscosity(int fluidity_species, /* integer associated with conc eqn for bond */
                   VISCOSITY_DEPENDENCE_STRUCT *d_mu) {
  /* Local Variables */
  if (!pd->gv[MASS_FRACTION]) {
    GOMA_EH(GOMA_ERROR, "Fluidity viscosity but no Species Concentration Equation\n");
  }

  dbl *params = mp->u_species_source[fluidity_species];

  dbl phi_0 = params[0];
  dbl phi_inf = params[1];

  dbl phi_star = fv->c[fluidity_species];

  bool disable_deriv = false;
  if (phi_star < 0) {
    phi_star = 0;
    disable_deriv = true;
  } else if (phi_star > 1) {
    phi_star = 1;
    disable_deriv = true;
  }

  /* Calculate the fluidity */
  dbl phi = phi_0 + (phi_inf - phi_0) * phi_star;
  if (d_mu != NULL && !disable_deriv) {
    for (int j = 0; j < ei[pg->imtrx]->dof[MASS_FRACTION]; j++)
      d_mu->C[fluidity_species][j] = bf[MASS_FRACTION]->phi[j] * -(phi_inf - phi_0) / (phi * phi);
  }

  return 1.0 / phi;
}

sfad sfad_fluidity_source(int species_no, dbl *params) {
  dbl phi_0 = params[0];
  dbl phi_inf = params[1];
  dbl K = params[2];
  dbl n = params[3];
  dbl tc = params[4];
  dbl sigma_y = params[5];
  dbl m = params[6];
  dbl m_y = params[7];

  sfad phi_star(10, 0, fv->c[species_no]);
  sfad grad_v[DIM][DIM];
  for (int i = 0; i < VIM; i++) {
    for (int j = 0; j < VIM; j++) {
      grad_v[i][j] = fv->grad_v[i][j];
      int idx = 1 + i * VIM + j;
      grad_v[i][j].fastAccessDx(idx) = 1.0;
    }
  }

  sfad phi_star_c = std::min(1, std::max(0, phi_star));

  sfad phi = phi_0 + (phi_inf - phi_0) * phi_star_c;
  sfad mu = 1 / phi;

  sfad tau[DIM][DIM] = {{0.}};
  for (int i = 0; i < pd->Num_Dim; i++) {
    for (int j = 0; j < pd->Num_Dim; j++) {
      tau[i][j] = mu * (grad_v[i][j] + grad_v[j][i]);
    }
  }

  sfad sigma = 0.0;
  for (int i = 0; i < DIM; i++) {
    for (int j = 0; j < DIM; j++) {
      sigma += tau[i][j] * tau[i][j];
    }
  }
  sigma = max(sigma, 1e-32);
  sigma = sqrt(0.5 * sigma);

  // Calculate phi_eq
  sfad term = pow((sigma - sigma_y) / K, 1.0 / n) * (1.0 / sigma);
  sfad phi_eq =
      max(1e-32, smoothed_heaviside(sigma - sigma_y, m_y) * term / ((phi_inf - phi_0) + term));

  // Calculate t_a and s
  // avalanche time
  sfad phi_eq_floor = std::max(phi_eq, 1e-5);
  sfad phi_eq_floor2 = std::max(phi_eq, 1e-10);
  sfad ta = 59.2 * pow(1 - phi_eq_floor2, 1.1) / (pow(phi_eq_floor2, 0.4));
  // sfad ta = ta;
  sfad s = 8.0 / (exp(phi_eq_floor / 0.09) - 1.0) + 1.2;

  sfad H = smoothed_heaviside(phi_star - phi_eq_floor2, m);

  // Breakdown dynamics
  sfad F_breakdown = (s / (ta * phi_eq_floor)) *
                     pow(std::max(0, phi_eq_floor2 - phi_star_c), (s + 1) / s) *
                     pow(max(0, phi_star_c), (s - 1.0) / s);

  // Buildup dynamics
  sfad F_buildup = -(phi_star - phi_eq_floor2) / tc;

  // Calculate F
  sfad F = H * F_buildup + (1.0 - H) * F_breakdown;

  return F;
}

/*
 * Calculate the fluidity source term in a Goma 7 compatible way
 * Wraps sfad_fluidity_source to use AD for Jacobian calculations
 */
extern "C" int fluidity_source(int species_no, struct Species_Conservation_Terms *st, dbl *params) {
  sfad F = sfad_fluidity_source(species_no, params);
  int eqn = MASS_FRACTION;
  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
    st->MassSource[species_no] = F.val();

    /* Jacobian entries for source term */
    int var = MASS_FRACTION;
    if (pd->v[pg->imtrx][var]) {
      for (int j = 0; j < ei[pg->imtrx]->dof[MASS_FRACTION]; j++) {
        st->d_MassSource_dc[species_no][species_no][j] = bf[MASS_FRACTION]->phi[j] * F.dx(0);
      }
    }

    var = VELOCITY1;
    if (pd->v[pg->imtrx][var]) {
      for (int a = 0; a < VIM; a++) {
        for (int p = 0; p < VIM; p++) {
          for (int q = 0; q < VIM; q++) {
            for (int j = 0; j < ei[pg->imtrx]->dof[VELOCITY1 + a]; j++) {
              st->d_MassSource_dv[species_no][a][j] +=
                  F.dx(1 + p * VIM + q) * bf[VELOCITY1 + a]->grad_phi_e[j][a][p][q];
            }
          }
        }
      }
    }
    var = MESH_DISPLACEMENT1;
    if (pd->v[pg->imtrx][var]) {
      for (int a = 0; a < pd->Num_Dim; a++) {
        for (int p = 0; p < pd->Num_Dim; p++) {
          for (int q = 0; q < pd->Num_Dim; q++) {
            for (int j = 0; j < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1 + a]; j++) {
              st->d_MassSource_dmesh[species_no][a][j] +=
                  F.dx(1 + p * VIM + q) * fv->d_grad_v_dmesh[p][q][a][j];
            }
          }
        }
      }
    }
  }
  return 0;
}

/*
 * fluidity_equilibrium_surf
 *
 * This function applies the equilibrium fluidity to a boundary
 * The functional form is dc/dt = F
 * where F is the fluidity source term
 *
 * Jacobian entries are calculated using AD from Sacado
 *
 */
extern "C" void fluidity_equilibrium_surf(double func[DIM],
                                          double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                                          int fluidity_species,
                                          const double theta,
                                          const double dt) {

  int eqn = R_MASS;
  sfad F = sfad_fluidity_source(fluidity_species, mp->u_species_source[fluidity_species]);
  sfad c_dot = fv_dot->c[fluidity_species];
  c_dot.fastAccessDx(0) = (1 + 2 * theta) / dt;
  sfad func_fad = 0;
  if (pd->e[pg->imtrx][eqn] & T_MASS) {
    if (pd->TimeIntegration != STEADY) {
      func_fad = -c_dot;
    }
  }

  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
    func_fad += F;
  }

  func[0] = func_fad.val();
  /* Jacobian entries for source term */
  int var = MASS_FRACTION;
  if (pd->v[pg->imtrx][var]) {
    for (int j = 0; j < ei[pg->imtrx]->dof[MASS_FRACTION]; j++) {
      d_func[0][MAX_VARIABLE_TYPES + fluidity_species][j] =
          bf[MASS_FRACTION]->phi[j] * func_fad.dx(0);
    }
  }

  var = VELOCITY1;
  if (pd->v[pg->imtrx][var]) {
    for (int a = 0; a < VIM; a++) {
      for (int p = 0; p < VIM; p++) {
        for (int q = 0; q < VIM; q++) {
          for (int j = 0; j < ei[pg->imtrx]->dof[VELOCITY1 + a]; j++) {
            d_func[0][VELOCITY1 + a][j] +=
                func_fad.dx(1 + p * VIM + q) * bf[VELOCITY1 + a]->grad_phi_e[j][a][p][q];
          }
        }
      }
    }
  }
}

#endif // GOMA_ENABLE_SACADO