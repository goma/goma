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

/* rd_exo.h - prototype declarations for rd_exo.c
 */

#ifndef _RD_EXO_H
#define _RD_EXO_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _RD_EXO_C
#define EXTERN /* do nothing */
#endif

#ifndef _RD_EXO_C
#define EXTERN extern
#endif

EXTERN void alloc_init_exo_nv_indeces /* rd_exo.c                            */
PROTO((Exo_DB *));		/* exo                                       */

EXTERN void free_exo_ev		/* rd_exo.c -- free up big elem var results */
PROTO((Exo_DB *));		/* ptr to FE database described exo_struct.h */

EXTERN void free_exo_gv		/* rd_exo.c -- free up big glob var results */
PROTO((Exo_DB *));		/* ptr to FE database described exo_struct.h */

EXTERN void free_exo_nv		/* rd_exo.c -- free up big node var results */
PROTO((Exo_DB *));		/* ptr to FE database described exo_struct.h */

EXTERN void alloc_exo_ev	/* rd_exo.c -- get mem big elemvar results */
PROTO((Exo_DB *,		/* ptr to FE database described exo_struct.h */
       int ));			/* num_timeplanes */

EXTERN void alloc_exo_gv	/* rd_exo.c -- get mem big globvar results */
PROTO((Exo_DB *,		/* ptr to FE database described exo_struct.h */
       int ));			/* num_timeplanes */

EXTERN void alloc_exo_nv	/* rd_exo.c -- get mem big nodevar results */
PROTO((Exo_DB *));		/* ptr to FE database described exo_struct.h */

EXTERN void init_exo_struct	/* rd_exo.c */
PROTO((Exo_DB *));		/* ptr to FE database described exo_struct.h */


EXTERN int rd_exo
PROTO((Exo_DB *,
       char *,
       int ,
       int ));

EXTERN int free_exo
PROTO((Exo_DB *));		/* pointer to EXODUS II FE db structure */

#endif /* _RD_EXO_H */
