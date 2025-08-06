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

#ifndef GOMA_BC_GENERALIZED_DIRICHLET_H
#define GOMA_BC_GENERALIZED_DIRICHLET_H

int fgeneralized_dirichlet(double *func,
                           double d_func[],        /* MAX_VARIABLE_TYPES + MAX_CONC */
                           const int gd_condition, /* denoting which condition
                                                    * applied */
                           const int bc_input_id,
                           const double tt,  /* parameter to vary time integration
                                              * from explicit (tt = 1) to
                                              * implicit (tt = 0) */
                           const double dt); /* current time step size          */

int evaluate_time_func(const double current_time,
                       double *f_time, /* computed time function */
                       const int bc_input_id);

#endif /* GOMA_BC_GENERALIZED_DIRICHLET_H */