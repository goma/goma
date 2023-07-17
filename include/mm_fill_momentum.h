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
#ifndef GOMA_MM_FILL_MOMENTUM_H
#define GOMA_MM_FILL_MOMENTUM_H

#include <stdbool.h>

#include "el_elm.h"
#include "exo_struct.h"
#include "mm_as_structs.h"
#include "mm_fill_common.h"
#include "rf_fem_const.h"
#include "std.h"

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MM_FILL_MOMENTUM_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_MM_FILL_MOMENTUM_C
#define EXTERN extern
#endif

EXTERN int assemble_momentum /* mm_fill_terms.c                           */
    (double,                 /* time - present time value                 */
     dbl,                    /* tt - parm to vary time integration from
                              * explicit (tt = 1) to implicit (tt = 0)    */
     dbl,                    /* dt - current time step size               */
     dbl,                    /* h_elem_avg - average global element size  */
     const PG_DATA *,        /* dvc_dnode                                 */
     dbl xi[DIM],            /* Local stu coords */
     const Exo_DB *exo);     /* ExodusII database struct pointer */

EXTERN void fluid_stress(double[DIM][DIM],            /* Pi[DIM][DIM] */
                         STRESS_DEPENDENCE_STRUCT *); /* d_Pi         */

void ve_polymer_stress(double gamma[DIM][DIM],
                       double stress[DIM][DIM],
                       STRESS_DEPENDENCE_STRUCT *d_stress);

int momentum_source_term /* mm_fill_terms.c                           */
    (dbl[DIM],           /* f - Body force.                           */
     MOMENTUM_SOURCE_DEPENDENCE_STRUCT *,
     double);
#endif /* GOMA_MM_FILL_MOMENTUM_H */
