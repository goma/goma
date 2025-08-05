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

#ifndef GOMA_LOAD_VARIABLE_H
#define GOMA_LOAD_VARIABLE_H

int load_variable(double *x_var,        /* variable value */
                  double *d_x_var,      /* sensitivities of variable value */
                  int jvar,             /* variable number */
                  int wspec,            /* species number */
                  double tt,            /* parameter to vary time integration from
                                           explicit (tt = 1) to implicit (tt = 0) */
                  double dt,            /* current time step size */
                  double d_vect_var[]); /* vector sensitivities  */

#endif /* GOMA_LOAD_VARIABLE_H */