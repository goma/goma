/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2025 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/

#ifndef GOMA_BC_GEOM_H
#define GOMA_BC_GEOM_H

#include "std.h"

void moving_plane(int ielem_dim, double *func, double d_func[], dbl *aa, double time);

void fplane(int ielem_dim,
            double *func,
            double d_func[], /* dimensioned [MAX_VARIABLE_TYPES+MAX_CONC] */
            dbl *aa);        /*  function parameters from data card  */

void f_fillet(const int ielem_dim,
              double *func,
              double d_func[],      /* dimensioned [MAX_VARIABLE_TYPES+MAX_CONC] */
              const double *p,      /*  function parameters from data card  */
              const int num_const); /* number of passed parameters   */

void f_double_rad(const int ielem_dim,
                  double *func,
                  double d_func[],      /* dimensioned [MAX_VARIABLE_TYPES+MAX_CONC] */
                  const double *p,      /*  function parameters from data card  */
                  const int num_const); /* number of passed parameters   */

void f_double_fillet(const int ielem_dim,
                     double *func,
                     double d_func[],      /* dimensioned [MAX_VARIABLE_TYPES+MAX_CONC] */
                     const double *p,      /*  function parameters from data card  */
                     const int num_const); /* number of passed parameters   */

void f_double_fillet_geom_based(const int ielem_dim,
                                double *func,
                                double d_func[], /* dimensioned [MAX_VARIABLE_TYPES+MAX_CONC] */
                                const double *p, /*  function parameters from data card  */
                                const int num_const); /* number of passed parameters   */

void f_roll_fluid(int ielem_dim,
                  double *func,
                  double d_func[],     /* dimensioned [MAX_VARIABLE_TYPES+MAX_CONC] */
                  const double *p,     /*  function parameters from data card  */
                  const int num_const, /* number of passed parameters   */
                  double *xsurf);      /* number of passed parameters   */

#ifdef FEATURE_ROLLON_PLEASE
EXTERN void f_feature_rollon(const int, /* ielem_dim */
                             double *,  /* func */
                             double[],  /* d_func - dimensioned [MAX_VARIABLE_TYPES+MAX_CONC] */
                             const double *, /* p - function parameters from data card  */
                             const int,      /* number of parameters from bc card  */
                             const int,      /* geometry model id  */
                             const double);  /* time - time at which BC's are evaluated  */
#endif

void fspline(int ielem_dim,
             double *func,
             double d_func[], /* dimensioned [MAX_VARIABLE_TYPES + MAX_CONC] */
             double p[],      /* parameters to parameterize temperature eqn model*/
             double time);    /* time  at which bc's are evaluated     */

void fspline_rs(int ielem_dim,
                double *func,
                double d_func[], /* dimensioned [MAX_VARIABLE_TYPES + MAX_CONC] */
                double p[],      /* parameters to parameterize temperature eqn model*/
                double time);    /* time  at which bc's are evaluated     */

#endif /* GOMA_BC_GEOM_H */