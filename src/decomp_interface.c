/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2026 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "decomp_interface.h"
#include "mm_as.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "rf_io.h"
#include "std.h"

#define MAX_DECOMP_COMMAND 2048

void decompose_exodus_files(void) {
  int i;

  if (strcmp(ExoAuxFile, "") != 0) {
    if (Debug_Flag) {
      DPRINTF(stdout, "Decomposing exodus file %s\n", ExoAuxFile);
    }
    //    ioss_decompose_mesh(ExoAuxFile);
  }

  if (efv->Num_external_field != 0) {
    for (i = 0; i < efv->Num_external_field; i++) {
      if (Debug_Flag) {
        DPRINTF(stdout, "Decomposing exodus file %s\n", efv->file_nm[i]);
      }
      //     ioss_decompose_mesh(efv->file_nm[i]);
    }
  }
  if (Debug_Flag) {
    DPRINTF(stdout, "Decomposing exodus file %s\n", ExoFile);
  }
  // ioss_decompose_mesh(ExoFile);
}
