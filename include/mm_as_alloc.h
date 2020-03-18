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
 
#ifndef GOMA_MM_AS_ALLOC_H
#define GOMA_MM_AS_ALLOC_H

#include "exo_struct.h"
#include "md_timer.h"

struct Problem_Description;
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MM_AS_ALLOC_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_MM_AS_ALLOC_C
#define EXTERN extern
#endif

EXTERN int pd_alloc
(void );

EXTERN int efv_alloc
(void );

EXTERN int mp_alloc
(void );

EXTERN int cr_alloc
(void );

EXTERN int gn_alloc
(void );

EXTERN int ve_alloc
(void );

EXTERN int elc_alloc
(void );

EXTERN int evp_alloc
(void );

EXTERN int evp_tensor_alloc
(Exo_DB *);		/* exo - ptr to std FE db */

EXTERN int elc_rs_alloc
(void );

EXTERN int tran_alloc
(void );

EXTERN int libio_alloc
(void );

EXTERN int eigen_alloc
(void );

EXTERN int cont_alloc
(void );

EXTERN int loca_alloc
(void );

EXTERN int assembly_alloc
(Exo_DB *);

EXTERN int bf_init
(Exo_DB *);

EXTERN int bf_mp_init
(struct Problem_Description *); /* pd - std ptr to global beast */

EXTERN void bf_reset
(void);

#endif /* GOMA_MM_AS_ALLOC_H */
