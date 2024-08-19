#ifndef GOMA_MM_FILL_STABILIZATION_H
#define GOMA_MM_FILL_STABILIZATION_H

#include "el_elm.h"
#include "mm_as_structs.h"
#include "std.h"

enum pspg_type {
  PSPG_NONE = 0,
  PSPG_GLOBAL = 1,
  PSPG_LOCAL = 2,
  PSPG_SHAKIB = 3,
};

typedef struct {
  dbl supg_tau;
  dbl d_supg_tau_dv[DIM][MDE];
  dbl d_supg_tau_dX[DIM][MDE];
} SUPG_terms;

typedef struct {
  dbl tau;
  dbl d_tau_dv[DIM][MDE];
  dbl d_tau_dX[DIM][MDE];
  dbl d_tau_dT[MDE];           /* temperature dependence. */
  dbl d_tau_dC[MAX_CONC][MDE]; /* conc dependence. */
  dbl d_tau_dP[MDE];           /* pressure dependence. */
  dbl d_tau_dF[MDE];           /* FILL dependence. */
  dbl d_tau_dnn[MDE];          /* bond concentration dependence */
  dbl d_tau_dEDDY_NU[MDE];     /* Turbulent viscosity */
  dbl d_tau_dturb_k[MDE];      /* Turbulent viscosity */
  dbl d_tau_dturb_omega[MDE];  /* Turbulent viscosity */
} momentum_tau_terms;

void supg_tau_shakib(SUPG_terms *supg_terms, int dim, dbl dt, dbl diffusivity, int interp_eqn);
void tau_momentum_shakib(momentum_tau_terms *tau_terms, int dim, dbl dt, int pspg_scale);

void get_supg_tau(SUPG_terms *supg_terms, int dim, dbl diffusivity, PG_DATA *pg_data);

void supg_tau_gauss_point(SUPG_terms *supg_terms, int dim, dbl diffusivity, const PG_DATA *pg_data);

void supg_tau(SUPG_terms *supg_terms,
              int dim,
              dbl diffusivity,
              const PG_DATA *pg_data,
              dbl dt,
              int shakib,
              int interp_eqn);

dbl yzbeta(dbl scale,
           int dim,
           dbl Y,
           dbl Z,
           dbl d_Z[MDE],
           dbl beta,
           dbl u,
           dbl d_u[MDE],
           dbl grad_u[DIM],
           dbl d_grad_u[MDE][DIM],
           dbl h_elem,
           int interp_eqn,
           dbl deriv[MDE]);

dbl yzbeta_model(int model,
                 dbl scale,
                 dbl beta,
                 int dim,
                 dbl Y,
                 dbl Z,
                 dbl d_Z[MDE],
                 dbl u,
                 dbl d_u[MDE],
                 dbl grad_u[DIM],
                 dbl d_grad_u[MDE][DIM],
                 dbl h_elem,
                 int interp_eqn,
                 dbl deriv[MDE]);

void get_metric_tensor(dbl B[DIM][DIM], int dim, int element_type, dbl G[DIM][DIM]);

int calc_pspg(dbl pspg[DIM],
              PSPG_DEPENDENCE_STRUCT *d_pspg,
              dbl time_value, /* current time */
              dbl tt,         /* parameter to vary time integration from
                                                 explicit (tt = 1) to implicit (tt = 0)    */
              dbl dt,         /* current time step size                    */
              const PG_DATA *pg_data);

int calc_cont_gls /* mm_fill_terms.c                           */
    (dbl *,
     CONT_GLS_DEPENDENCE_STRUCT *,
     dbl, /* current time                              */
     const PG_DATA *);

bool is_evss_f_model(int model);
#endif // GOMA_MM_FILL_STABILIZATION_H
