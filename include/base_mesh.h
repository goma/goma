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

#ifndef GOMA_BASE_MESH_H
#define GOMA_BASE_MESH_H

#include "dpi.h"
#include "exo_struct.h"
#include "mm_eh.h"

goma_error setup_base_mesh(Dpi *dpi, Exo_DB *exo, int num_proc);
goma_error free_base_mesh(Exo_DB *exo);

#endif // GOMA_BASE_MESH_H
