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

#ifndef GOMA_EL_QUALITY_H
#define GOMA_EL_QUALITY_H

#include "el_elm_info.h"
#include "exo_struct.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_EL_QUALITY_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_EL_QUALITY_C
#define EXTERN extern
#endif

EXTERN int element_quality(Exo_DB *, /* Exodus database structure */
                           double *, /* Solution vector */
                           int *);   /* proc_config array */

#endif                               /* GOMA_EL_QUALITY_H */
