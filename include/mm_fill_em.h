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

#ifndef GOMA_MM_FILL_EM_H
#define GOMA_MM_FILL_EM_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MM_FILL_EM_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_MM_FILL_EM_C
#define EXTERN extern
#endif

enum stabilization_type {
  EM_STAB_NONE,
  EM_STAB_PHI_DIV,
  EM_STAB_DPHI_DIV,
  EM_STAB_DIVPHI_DIV,
  EM_STAB_PHI_DIVSQUARED,
  EM_STAB_DPHI_DIVSQUARED
};

struct emwave_stabilization {
  enum stabilization_type type;
  int em_eqn;
  int em_var;
  int stabilization_field_var;
  double residual_term[MDE];
  double jacobian_term[MDE][DIM][MDE];
};

EXTERN int assemble_emwave /* mm_fill_em.c                           */
    (double,               /* time - present time value         */
     double,               /* tt - parameter to vary time integration
                            * from explicit (tt = 1) to
                            * implicit (tt = 0)                   */
     double,               /* dt - current time step size               */
     const PG_DATA *,      /* dvc_dnode                                 */
     const int,            /*  Light intensity eqn id and var id     */
     const int,            /*  Light intensity eqn id and var id     */
     const int);

EXTERN int apply_em_farfield_direct_vec               /* mm_fill_em.c                           */
    (double[DIM],                                     // func
     double[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], // d_func
     double[DIM],                                     // xi
     const int,                                       // bc_name
     double *);

EXTERN int apply_em_sommerfeld_vec                    /* mm_fill_em.c                           */
    (double[DIM],                                     // func
     double[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], // d_func
     double[DIM],                                     // xi
     const int,                                       // bc_name
     double *);

EXTERN int apply_em_free_vec                          /* mm_fill_em.c                           */
    (double[DIM],                                     // func
     double[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], // d_func
     double[DIM],                                     // xi
     const int);                                      // bc_name

EXTERN int apply_ewave_planewave_vec                  /* mm_fill_em.c                           */
    (double[DIM],                                     // func
     double[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], // d_func
     double[DIM],                                     // xi
     const int,                                       // bc_name
     double *);

EXTERN int apply_ewave_curlcurl_farfield_vec          /* mm_fill_em.c                           */
    (double[DIM],                                     // func
     double[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], // d_func
     double[DIM],                                     // xi
     double time,                                     // present time
     const int,                                       // bc_name
     double *);

EXTERN int apply_ewave_2D                             /* mm_fill_em.c                           */
    (double[DIM],                                     // func
     double[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], // d_func
     double[DIM],                                     // xi
     const int);                                      // bc_name

EXTERN void calc_emwave_stabilization_term(struct emwave_stabilization *, double);

int assemble_ewave_curlcurl(double time,       // present time
                            double tt,         // time integration method parameter
                            double dt,         // current time step size
                            const int em_eqn,  // eqn id
                            const int em_var); //  variable id - should match me_eqn

int assemble_ewave_laplacian(double time,       // present time
                             double tt,         // time integration method parameter
                             double dt,         // current time step size
                             const int em_eqn,  // eqn id
                             const int em_var); //  variable id - should match me_eqn

int assemble_em_continuity(void);
#endif /* GOMA_MM_FILL_EM_H */
