#ifndef GOMA_AD_TURBULENCE_H
#define GOMA_AD_TURBULENCE_H

#ifdef GOMA_ENABLE_SACADO

#ifdef __cplusplus
#include <Sacado.hpp>
extern "C" {
#include "el_elm.h"
#include "mm_mp_const.h"
#include "std.h"
}
using ADType = Sacado::Fad::DFad<double>;
void ad_supg_tau_shakib(ADType &supg_tau, int dim, dbl dt, ADType diffusivity, int interp_eqn);
struct AD_Basis {
  ADType d_phi[MDE][DIM];                /* d_phi[i][a]    = d(phi_i)/d(q_a) */
  ADType grad_phi[MDE][DIM];             /* grad_phi[i][a] = e_a . grad(phi_i) */
  ADType grad_phi_e[MDE][DIM][DIM][DIM]; /* grad_phi_e[i][a][p][q] */
                                         /* = (e_p e_q): grad(phi_i e_a) */
  ADType curl_phi_e[MDE][DIM][DIM];
};
struct AD_Field_Variables {
  AD_Field_Variables() = default;
  std::vector<AD_Basis> basis;
  ADType detJ;
  ADType J[DIM][DIM];
  ADType B[DIM][DIM];
  ADType v[DIM];
  ADType v_dot[DIM];
  ADType x[DIM];
  ADType d[DIM];
  ADType x_dot[DIM];
  ADType grad_v[DIM][DIM];
  ADType G[DIM][DIM];
  ADType grad_G[DIM][DIM][DIM];
  ADType div_G[DIM];
  ADType S[MAX_MODES][DIM][DIM];
  ADType S_dot[MAX_MODES][DIM][DIM];
  ADType grad_S[MAX_MODES][DIM][DIM][DIM];
  ADType div_S[MAX_MODES][DIM];
  ADType grad_SH[DIM];
  ADType P;
  ADType SH;
  ADType grad_P[DIM];
  ADType eddy_nu;
  ADType eddy_nu_dot;
  ADType grad_eddy_nu[DIM];
  ADType turb_k;
  ADType turb_k_dot;
  ADType grad_turb_k[DIM];
  ADType turb_omega;
  ADType turb_omega_dot;
  ADType grad_turb_omega[DIM];
  int total_ad_variables;
  int ielem;
  int offset[V_LAST];
};

extern AD_Field_Variables *ad_fv;
int ad_calc_shearrate(ADType &gammadot,            /* strain rate invariant */
                      ADType gamma_dot[DIM][DIM]); /* strain rate tensor */

void ad_only_tau_momentum_shakib(ADType &tau, int dim, dbl dt, int pspg_scale);
ADType ad_sa_viscosity(struct Generalized_Newtonian *gn_local);
ADType ad_only_turb_k_omega_viscosity(void);
void compute_sst_blending(ADType &F1, ADType &F2);
ADType sst_viscosity(const ADType &Omega, const ADType &F2);
extern "C" {
#endif

#include "mm_as_structs.h"
#include "mm_fill_stabilization.h"
#include "mm_mp_structs.h"
#include "std.h"

void ad_tau_momentum_shakib(momentum_tau_terms *tau_terms, int dim, dbl dt, int pspg_scale);

int ad_assemble_turb_k(dbl time_value, /* current time */
                       dbl tt,         /* parameter to vary time integration from
                                          explicit (tt = 1) to implicit (tt = 0)    */
                       dbl dt,         /* current time step size                    */
                       const PG_DATA *pg_data);

void ad_sa_wall_func(double func[DIM], double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE]);
dbl ad_turb_k_omega_sst_viscosity(VISCOSITY_DEPENDENCE_STRUCT *d_mu);

int ad_assemble_turb_omega(dbl time_value, /* current time */
                           dbl tt,         /* parameter to vary time integration from
                                              explicit (tt = 1) to implicit (tt = 0)    */
                           dbl dt,         /* current time step size                    */
                           const PG_DATA *pg_data);
dbl ad_sa_viscosity(struct Generalized_Newtonian *gn_local, VISCOSITY_DEPENDENCE_STRUCT *d_mu);
void fill_ad_field_variables();
int ad_assemble_spalart_allmaras(dbl time_value, /* current time */
                                 dbl tt,         /* parameter to vary time integration from
                                                    explicit (tt = 1) to implicit (tt = 0)    */
                                 dbl dt,         /* current time step size                    */
                                 const PG_DATA *pg_data);
int ad_assemble_turb_k_modified(dbl time_value, /* current time */
                                dbl tt,         /* parameter to vary time integration from
                                                   explicit (tt = 1) to implicit (tt = 0)    */
                                dbl dt,         /* current time step size                    */
                                const PG_DATA *pg_data);
int ad_assemble_turb_omega_modified(dbl time_value, /* current time */
                                    dbl tt,         /* parameter to vary time integration from
                                                       explicit (tt = 1) to implicit (tt = 0)    */
                                    dbl dt,         /* current time step size                    */
                                    const PG_DATA *pg_data);
int ad_assemble_turb_k_omega_modified(dbl time_value, /* current time */
                                      dbl tt,         /* parameter to vary time integration from
                                                         explicit (tt = 1) to implicit (tt = 0)    */
                                      dbl dt, /* current time step size                    */
                                      const PG_DATA *pg_data);
int ad_assemble_k_omega_sst_modified(dbl time_value, /* current time */
                                     dbl tt,         /* parameter to vary time integration from
                                                        explicit (tt = 1) to implicit (tt = 0)    */
                                     dbl dt,         /* current time step size                    */
                                     const PG_DATA *pg_data);
int ad_assemble_invariant(double tt,  /* parameter to vary time integration from
                                       * explicit (tt = 1) to implicit (tt = 0)    */
                          double dt); /*  time step size                          */
void ad_omega_wall_func(double func[DIM], double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE]);

#ifdef __cplusplus
}
#endif

#endif

#endif // GOMA_AD_TURBULENCE_H