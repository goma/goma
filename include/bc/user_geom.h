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

#ifndef GOMA_BC_USER_GEOM_H
#define GOMA_BC_USER_GEOM_H

#include "std.h"

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

#endif /* GOMA_BC_USER_GEOM_H */