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

#ifndef GOMA_ADAPT_RESETUP_PROBLEM_H
#define GOMA_ADAPT_RESETUP_PROBLEM_H

#include "dpi.h"
#include "exo_struct.h"
#include "mm_as_structs.h"

struct GomaLinearSolverData;

int resetup_problem(Exo_DB *exo, Dpi *dpi); /* ptr to the finite element mesh database */
int resetup_matrix(struct GomaLinearSolverData **ams, Exo_DB *exo, Dpi *dpi);

#endif // GOMA_ADAPT_RESETUP_PROBLEM_H

// vim: expandtab sw=2 ts=8
