/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) Goma Developers, National Technology & Engineering        *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/

#ifndef GOMA_SL_MUMPS_H
#define GOMA_SL_MUMPS_H

#include "mm_eh.h"
#include "sl_util_structs.h"
#include "std.h"

goma_error mumps_solve(struct GomaLinearSolverData *data, dbl *x, dbl *rhs);

#endif // GOMA_SL_MUMPS_H