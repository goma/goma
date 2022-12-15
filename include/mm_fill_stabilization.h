#ifndef GOMA_MM_FILL_STABILIZATION_H
#define GOMA_MM_FILL_STABILIZATION_H

#include "el_elm.h"
#include "mm_as_structs.h"
#include "std.h"

typedef struct {
  dbl supg_tau;
  dbl d_supg_tau_dv[DIM][MDE];
  dbl d_supg_tau_dX[DIM][MDE];
} SUPG_terms;

typedef struct {
  dbl supg_tau;
  dbl d_supg_tau_dv[DIM][MDE];
  dbl d_supg_tau_dX[DIM][MDE];
  dbl d_supg_tau_dC[MAX_CONC][MDE];
  dbl d_supg_tau_dP[MDE];
  dbl d_supg_tau_dF[MDE];
  dbl d_supg_tau_dT[MDE];
  dbl d_supg_tau_dnn[MDE];
} SUPG_momentum_terms;

void supg_tau_momentum_shakib(SUPG_momentum_terms *supg_terms, int dim, dbl dt);

void supg_tau_shakib(SUPG_terms *supg_terms, int dim, dbl dt, dbl diffusivity, int interp_eqn);

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

#endif // GOMA_MM_FILL_STABILIZATION_H
