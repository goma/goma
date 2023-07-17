/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2023 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/

#ifndef GOMA_DENSITY_H
#define GOMA_DENSITY_H

#include "el_elm.h"
#include "exo_struct.h"
#include "mm_as_structs.h"
#include "mm_elem_block_structs.h"
#include "mm_fill_common.h"
#include "mm_mp_const.h"
#include "rf_bc_const.h"
#include "rf_fem_const.h"
#include "rf_io_const.h"
#include "rf_vars_const.h"
#include "sl_util_structs.h"
#include "std.h"
#include <stdbool.h>

struct density_dependence {
  double T[MDE];
  double C[MAX_CONC][MDE];
  double F[MDE];
  double pf[MAX_PHASE_FUNC][MDE]; /* phase function */
  double moment[MAX_MOMENTS][MDE];
  double rho[MDE];
};
typedef struct density_dependence DENSITY_DEPENDENCE_STRUCT; /* struct for d_rho */

double density                    /* mm_fill_terms.c                           */
    (DENSITY_DEPENDENCE_STRUCT *, /* density dependence */
     double);
#endif // GOMA_DENSITY_H
