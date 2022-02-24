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

/* rd_in.h - prototype declarations for rd_in.c
 */

#ifndef GOMA_RD_IN_H
#define GOMA_RD_IN_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_RD_IN_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_RD_IN_C
#define EXTERN extern
#endif

EXTERN void rd_input(char *,    /* filename for input */
                     Exo_DB *,  /* EXODUS II finite element database */
                     Bevm ****, /* basic eqnvar multiplicity */
                     int ****,  /* eqnvar dependencies */
                     int ****,  /* local node/dof existences */
                     int **);   /* number of basic eqnvars in ea eb */

#endif /* GOMA_RD_IN_H */
