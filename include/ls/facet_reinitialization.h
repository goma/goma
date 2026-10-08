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

#ifndef GOMA_FACET_BASED_REINIT_HPP
#define GOMA_FACET_BASED_REINIT_HPP
#ifdef __cplusplus
extern "C" {
#endif

#include "dp_types.h"
#include "dpi.h"
#include "exo_struct.h"

void facet_based_reinitialization(
    double *x, Exo_DB *exo, Comm_Ex *cx, Dpi *dpi, int num_total_nodes, double time);

#ifdef __cplusplus
}
#endif
#endif // GOMA_FACET_BASED_REINIT_HPP