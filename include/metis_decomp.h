/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2022 Sandia Corporation.                                  *
*                                                                         *
* Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,  *
* the U.S. Government retains certain rights in this software.            *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
\************************************************************************/

#ifndef GOMA_METIS_DECOMP_H
#define GOMA_METIS_DECOMP_H

#include "mm_eh.h"

goma_error goma_metis_decomposition(char *filenames[], int n_files);

#endif // GOMA_METIS_DECOMP_H

