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
 
#ifndef _MM_AS_ALLOC_H
#define _MM_AS_ALLOC_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _MM_AS_ALLOC_C
#define EXTERN /* do nothing */
#endif

#ifndef _MM_AS_ALLOC_C
#define EXTERN extern
#endif

EXTERN int pd_alloc
PROTO((void ));

EXTERN int efv_alloc
PROTO((void ));

EXTERN int mp_alloc
PROTO((void ));

EXTERN int cr_alloc
PROTO((void ));

EXTERN int gn_alloc
PROTO((void ));

EXTERN int ve_alloc
PROTO((void ));

EXTERN int elc_alloc
PROTO((void ));

EXTERN int evp_alloc
PROTO((void ));

EXTERN int evp_tensor_alloc
PROTO((Exo_DB *));		/* exo - ptr to std FE db */

EXTERN int elc_rs_alloc
PROTO((void ));

EXTERN int tran_alloc
PROTO((void ));

EXTERN int libio_alloc
PROTO((void ));

EXTERN int eigen_alloc
PROTO((void ));

EXTERN int cont_alloc
PROTO((void ));

EXTERN int loca_alloc
PROTO((void ));

EXTERN int assembly_alloc
PROTO((Exo_DB *));

EXTERN int bf_init
PROTO((Exo_DB *));

EXTERN int bf_mp_init
PROTO((struct Problem_Description *)); /* pd - std ptr to global beast */

EXTERN void bf_reset
PROTO((void));

#endif /* _MM_AS_ALLOC_H */
