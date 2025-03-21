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

#ifndef GOMA_POLYMER_TIME_CONST_H
#define GOMA_POLYMER_TIME_CONST_H

#include "el_elm.h"
#include "mm_as_structs.h"
#include "mm_mp_structs.h"
#include "mm_qtensor_model.h"
#include "mm_std_models.h"
#include "std.h"
#ifdef EXTERN
#undef EXTERN
#endif

dbl power_law_time_const(struct PolymerTimeConstants *lambda_st,
                         dbl gamma_dot[DIM][DIM], /* strain rate tensor */
                         POLYMER_TIME_CONST_DEPENDENCE_STRUCT *d_lambda);
dbl carreau_polymer_time_const(struct PolymerTimeConstants *lambda_st,
                               dbl gamma_dot[DIM][DIM], /* strain rate tensor */
                               POLYMER_TIME_CONST_DEPENDENCE_STRUCT *d_lambda);
dbl polymer_time_const(struct PolymerTimeConstants *,
                       dbl[DIM][DIM],
                       POLYMER_TIME_CONST_DEPENDENCE_STRUCT *);

#endif /* GOMA_POLYMER_TIME_CONST_H */
