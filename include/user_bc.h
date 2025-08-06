/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2022 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/

#ifndef GOMA_USER_BC_H
#define GOMA_USER_BC_H

#include "el_elm.h"
#include "rf_fem_const.h"
#include "std.h"

extern dbl velo_vary_fnc(const int,   /* velo_condition                            */
                         const dbl,   /* x1                                        */
                         const dbl,   /* x2                                        */
                         const dbl,   /* x3                                        */
                         const dbl[], /* p                                         */
                         const dbl);  /* time                                      */

extern dbl dvelo_vary_fnc_d1(const int,   /* velo_condition                            */
                             const dbl,   /* x1                                        */
                             const dbl,   /* x2                                        */
                             const dbl,   /* x3                                        */
                             const dbl[], /* p                                         */
                             const dbl);  /* time                                      */

extern dbl dvelo_vary_fnc_d2(const int,   /* velo_condition                            */
                             const dbl,   /* x1                                        */
                             const dbl,   /* x2                                        */
                             const dbl,   /* x3                                        */
                             const dbl[], /* p                                         */
                             const dbl);  /* time                                      */

extern dbl dvelo_vary_fnc_d3(const int,   /* velo_condition                            */
                             const dbl,   /* x1                                        */
                             const dbl,   /* x2                                        */
                             const dbl,   /* x3                                        */
                             const dbl[], /* p                                         */
                             const dbl);  /* time                                      */

extern dbl fnc(const dbl,   /* x1                                        */
               const dbl,   /* x2                                        */
               const dbl,   /* x3                                        */
               const dbl[], /* p                                         */
               const dbl);  /* time                                      */

extern dbl dfncd1(const dbl,   /* x1                                        */
                  const dbl,   /* x2                                        */
                  const dbl,   /* x3                                        */
                  const dbl[], /* p                                         */
                  const dbl);  /* time                                      */

extern dbl dfncd2(const dbl,   /* x1                                        */
                  const dbl,   /* x2                                        */
                  const dbl,   /* x3                                        */
                  const dbl[], /* p                                         */
                  const dbl);  /* time                                      */

extern dbl dfncd3(const dbl,   /* x1                                        */
                  const dbl,   /* x2                                        */
                  const dbl,   /* x3                                        */
                  const dbl[], /* p                                         */
                  const dbl);  /* time                                      */

extern void quser_surf(double[DIM], /* func                                      */
                       double[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
                       double[],   /* p - parameterize heat transfer model      */
                       const dbl); /* time                                      */

extern void tuser(double *,       /* func                                      */
                  double[],       /* d_func - [MAX_VARIABLE_TYPES + MAX_CONC]  */
                  const double[], /* u_bc - parameterize temperature eqn model */
                  const double);  /* time */

extern void yuser_surf(double *, /* func                                      */
                       double[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func             */
                       const int,      /* species                                   */
                       const double[], /* u_bc - to parameterize species eqn model  */
                       const double);  /* time */

extern void y2_electroneutrality_surf(double *, /* func                                      */
                                      double[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func */
                                      const int, /* species                                   */
                                      const double[], /* u_bc - to parameterize species eqn model */
                                      const double);  /* time */

extern void uuser_surf(double[DIM], /* func                                      */
                       double[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
                       double[],   /* u_bc - parameterize u velocity model      */
                       const dbl); /* time                                      */

extern void vuser_surf(double[DIM], /* func                                      */
                       double[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
                       double[],   /* u_bc - parameterize u velocity model      */
                       const dbl); /* time                                      */

extern void wuser_surf(double[DIM], /* func                                      */
                       double[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
                       double[],   /* u_bc - parameterize u velocity model      */
                       const dbl); /* time  */

extern void uuser_colloc_surf(double *,       /* func                                      */
                              double[],       /* d_func           */
                              const double[], /* u_bc - parameterize u velocity model      */
                              const int,      /* Node ID */
                              const dbl);     /* time  */

extern void vuser_colloc_surf(double *,       /* func                                      */
                              double[],       /* d_func           */
                              const double[], /* u_bc - parameterize v velocity model      */
                              const int,      /* Node ID */
                              const dbl);     /* time  */

extern void wuser_colloc_surf(double *,       /* func                                      */
                              double[],       /* d_func           */
                              const double[], /* u_bc - parameterize w velocity model      */
                              const int,      /* Node ID */
                              const dbl);     /* time  */

extern void dx_user_surf(double *,       /* func                                      */
                         double[],       /* d_func           */
                         const double[], /* u_bc - parameterize u velocity model      */
                         const dbl);     /* time  */
extern void dy_user_surf(double *,       /* func                                      */
                         double[],       /* d_func           */
                         const double[], /* u_bc - parameterize u velocity model      */
                         const dbl);     /* time  */
extern void dz_user_surf(double *,       /* func                                      */
                         double[],       /* d_func           */
                         const double[], /* u_bc - parameterize u velocity model      */
                         const dbl);     /* time  */

extern void p_liq_user_surf(double *,       /* func                                      */
                            double[],       /* d_func           */
                            const double[], /* u_bc - parameterize u velocity model      */
                            const dbl);     /* time  */

extern void
shell_p_open_user_surf(double *func, double d_func[], const double u_bc[], const double time);

extern void fn_dot_T_user(double[DIM], /* func                                      */
                          double[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
                          const double[], /* u_bc                                      */
                          const dbl);     /* time                                      */

extern void flow_n_dot_T_user(double[DIM], /* func                                      */
                              double[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func */
                              const double[], /* u_BC - Parameters from input deck         */
                              const dbl);     /* time                                      */

extern double var_CA_user(double,         /* Ca_local */
                          int,            /* num */
                          const double *, /* a - parameter list from user */
                          double *);      /* d_cos_CA_Ca_local */

extern int user_gibbs_criterion(const double[MAX_PDIM], /* fsnormal - Vector of free surface normal
                                                         * components */
                                const double[MAX_PDIM], /* ssnormal - Vector of solid surface normal
                                                         * components */
                                const int,       /* imodel - Flag which tracks which model    */
                                int *,           /* ipin  - Flag for pinned or not            */
                                const double[]); /* p - User defined parameter list, or model
                                                  * spec. list                                */

extern void force_user_surf /* user_bc.c                                 */
    (double[DIM],           /* func                                      */
     double[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
     double[], /* p - user parameter list                   */
     dbl);     /* time                                      */

extern void volt_user_surf /* user_bc.c                                 */
    (double[DIM],          /* func                                      */
     double[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
     double[],   /* p - user parameter list                   */
     const dbl); /* time                                      */

extern void current_user_surf /* user_bc.c                                 */
    (double[DIM],             /* func                                      */
     double[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
     double[],   /* p - user parameter list                   */
     const dbl); /* time                                      */

extern void mass_flux_user_surf(double[MAX_CONC],
                                double[MAX_CONC][MAX_VARIABLE_TYPES + MAX_CONC],
                                const int,
                                const double[],
                                const double);

#endif /* GOMA_USER_BC_H */
