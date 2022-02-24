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

/* ppi.h - prototype declarations for ppi.c
 */

#ifndef GOMA_PPI_H
#define GOMA_PPI_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_PPI_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_PPI_C
#define EXTERN extern
#endif

#ifdef GOMA_PPI_C
const char filter[] = "sed -e 's/#.*$//' -e '/^[ 	]*$/d' -e 's/[ 	]*$//'";
#define TEMP_PREFIX "tmp."
#endif

EXTERN void brk_pre_process(char *); /* input filename */

#endif /* GOMA_PPI_H */
