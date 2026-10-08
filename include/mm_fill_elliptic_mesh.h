
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

#ifndef GOMA_MM_FILL_ELLIPTIC_MESH_H
#define GOMA_MM_FILL_ELLIPTIC_MESH_H

#include "el_elm.h"
#include "rf_fem_const.h"
#include "std.h"

int assemble_elliptic_mesh(void);

void assemble_essential_elliptic_mesh(dbl func[DIM],
                                      dbl d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                                      int bc_name,
                                      dbl M);
#endif /* GOMA_MM_FILL_ELLIPTIC_MESH_H */
