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
#ifndef GOMA_MM_FILL_LS_CAPILLARY_BCS_H
#define GOMA_MM_FILL_LS_CAPILLARY_BCS_H

#include "el_elm.h"
#include "exo_struct.h"
#include "mm_as_structs.h"
#include "mm_fill_common.h"
#include "rf_fem_const.h"
#include "std.h"
#include <stdbool.h>
int assemble_csf_tensor(void);
int assemble_div_n_source(void);
int assemble_div_s_n_source(void);
int assemble_cap_hysing(double dt, double scale);
int assemble_cap_denner_diffusion(double dt, double scale);
int assemble_cap_denner_diffusion_n(double dt, double scale);
#endif // GOMA_MM_FILL_LS_CAPILLARY_BCS_H
