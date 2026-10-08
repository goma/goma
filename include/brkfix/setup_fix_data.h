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

#ifndef GOMA_SETUP_FIX_DATA_H
#define GOMA_SETUP_FIX_DATA_H

#ifdef __cplusplus
extern "C" {
#endif
void setup_fix_data(const char *mono_name, int num_procs, struct fix_data *fd, int *pmax);
void free_fix_data(struct fix_data *fd);
#ifdef __cplusplus
}
#endif
#endif // GOMA_SETUP_FIX_DATA_H
