#include "Sacado.hpp"

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
  sfad d_term = -(n - 1.0) * pow(sigma / K, 1.0 / n) / (n * sigma * sigma);
  sfad phi_eq = term / ((phi_inf - phi_0) + term);
  return phi_eq;
}

// Smoothed Heaviside function using tanh
sfad smoothed_heaviside(sfad x, double m) { return 0.5 * (1.0 + tanh(m * x)); }

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

  /* Calculate the fluidity */
  dbl phi = phi_0 + (phi_inf - phi_0) * phi_star;
  if (d_mu != NULL) {
    for (int j = 0; j < ei[pg->imtrx]->dof[MASS_FRACTION]; j++)
      d_mu->C[fluidity_species][j] = bf[MASS_FRACTION]->phi[j] * -(phi_inf - phi_0) / (phi * phi);
  }

  return 1.0 / phi;
}

// Main function to test the implementation
extern "C" int fluidity_source(int species_no, struct Species_Conservation_Terms *st, dbl *params) {

  dbl phi_0 = params[0];
  dbl phi_inf = params[1];
  dbl K = params[2];
  dbl n = params[3];
  dbl tc = params[4];

  sfad phi_star(10, 0, fv->c[species_no]);
  sfad grad_v[DIM][DIM];
  for (int i = 0; i < VIM; i++) {
    for (int j = 0; j < VIM; j++) {
      grad_v[i][j] = fv->grad_v[i][j];
      int idx = 1 + i * VIM + j;
      grad_v[i][j].fastAccessDx(idx) = 1.0;
    }
  }

  sfad phi = phi_0 + (phi_inf - phi_0) * phi_star;
  sfad mu = 1/phi;

  sfad tau[DIM][DIM];
  for (int i = 0; i < DIM; i++) {
    for (int j = 0; j < DIM; j++) {
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
  sfad term = pow(sigma / K, 1.0 / n) * (1.0 / sigma);
  sfad d_term = -(n - 1.0) * pow(sigma / K, 1.0 / n) / (n * sigma * sigma);
  sfad phi_eq = max(1e-32, term / ((phi_inf - phi_0) + term));

  // Calculate t_a and s
  // avalanche time
  sfad ta = 59.2 * pow(1 - phi_eq, 1.1) / (pow(phi_eq, 0.4));
  sfad s = 8.0 / (exp(phi_eq / 0.09) - 1.0) + 1.2;

  sfad H = smoothed_heaviside(phi_star - phi_eq, 500);

  // Breakdown dynamics
  sfad F_breakdown =
      (s / (ta * phi_eq)) * pow(max(0,phi_eq - phi_star), (s + 1) / s) * pow(phi_star, (s - 1.0) / s);

  // Buildup dynamics
  sfad F_buildup = -(phi_star - phi_eq) / tc;

  // Calculate F
  sfad F = H * F_buildup + (1.0 - H) * F_breakdown;
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
  }
  return 0;
}