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

#ifndef GOMA_DP_GHOST_H
#define GOMA_DP_GHOST_H

#ifdef __cplusplus
extern "C" {
#endif

#include "dpi.h"
#include "exo_struct.h"
#include "mm_eh.h"

goma_error generate_ghost_elems(Exo_DB *exo, Dpi *dpi);
goma_error setup_ghost_to_base(Exo_DB *exo, Dpi *dpi);

#ifdef __cplusplus
};
#endif

#endif // GOMA_DP_GHOST_H
