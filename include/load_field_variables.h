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

#ifndef GOMA_LOAD_FIELD_VARIABLES_H
#define GOMA_LOAD_FIELD_VARIABLES_H

int load_fv /* mm_fill_terms.c                           */
    (void);

int load_fv_all(void);

int load_fv_vector(void);

int load_fv_grads /* mm_fill_terms.c                           */
    (void);

int load_fv_grads_all /* mm_fill_terms.c                           */
    (void);

int load_fv_mesh_derivs /* mm_fill_terms.c                           */
    (int);              /* okToZero - Turns on zeroing in the function
                                  This is usually turned on except when accumulating */

#endif /* GOMA_LOAD_FIELD_VARIABLES_H */