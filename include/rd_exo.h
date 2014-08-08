/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2014 Sandia Corporation.                                  *
*                                                                         *
* Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,  *
* the U.S. Government retains certain rights in this software.            *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
\************************************************************************/
 
#ifndef _RD_EXO_H
#define _RD_EXO_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _RD_EXO_C
#define EXTERN
#
#endif

#ifndef _RD_EXO_C
#define EXTERN extern
#endif

EXTERN int rd_exo
PROTO((Exo_DB *,
       const char *,
       const int ,
       const int ));

int 
brk_rd_exo(Exo_DB *x,		/* def'd in exo_struct.h */
           char *fn,
           int verbosity,
           int task);

EXTERN int free_exo
PROTO((Exo_DB *));		/* pointer to EXODUS II FE db structure */

EXTERN void zero_base
PROTO((Exo_DB *));

EXTERN void one_base
PROTO((Exo_DB *));

EXTERN int fence_post		/* rd_exo.c                                  */
PROTO((const int ,		/* val    - integer whose category we seek   */
       int *,			/* array  - where to look                    */
       const int ));		/* length - how far to search in array       */

EXTERN void init_exo_struct	/* rd_exo.c */
PROTO((Exo_DB *));		/* ptr to FE database described exo_struct.h */

EXTERN void free_exo_ev		/* rd_exo.c -- free up big elem var results */
PROTO((Exo_DB *));		/* ptr to FE database described exo_struct.h */

EXTERN void free_exo_gv		/* rd_exo.c -- free up big glob var results */
PROTO((Exo_DB *));		/* ptr to FE database described exo_struct.h */

EXTERN void free_exo_nv		/* rd_exo.c -- free up big node var results */
PROTO((Exo_DB *));		/* ptr to FE database described exo_struct.h */

EXTERN void alloc_exo_ev	/* rd_exo.c -- get mem big elemvar results */
PROTO((Exo_DB *,		/* ptr to FE database described exo_struct.h */
       const int ));		/* num_timeplanes */

EXTERN void alloc_exo_gv	/* rd_exo.c -- get mem big globvar results */
PROTO((Exo_DB *,		/* ptr to FE database described exo_struct.h */
       const int ));		/* num_timeplanes */

EXTERN void alloc_exo_nv	/* rd_exo.c -- get mem big nodevar results */
PROTO((Exo_DB *,		/* x - ptr to fe database */
       const int ,		/* num_timeplanes */
       const int ));		/* num_nodal_vars */

#endif /* _RD_EXO_H */
