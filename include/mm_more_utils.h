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
 
/*
 * mm_more_utils.h -- prototype declarations for mm_more_utils.c
 */

#ifndef GOMA_MM_MORE_UTILS_H
#define GOMA_MM_MORE_UTILS_H

#include "exo_struct.h"
#include "mm_input.h"
#include "std.h"

struct Results_Description;
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MM_MORE_UTILS_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_MM_MORE_UTILS_C
#define EXTERN extern
#endif

EXTERN int cnt_nodal_vars	/* mm_more_utils.c                           */
(void );

EXTERN int cnt_elem_vars	/* mm_more_utils.c                           */
(void );

EXTERN int goal_post_nodal	/* mm_more_utils.c                           */
(const int );		/* var  */

EXTERN int goal_post_elem	/* mm_more_utils.c                           */
(const int );		/* var  */

EXTERN int set_nv_tkud		/* mm_more_utils.c                           */
(struct Results_Description *, /* r - just node results for Exodus II  */
       const int ,		/* i - index                                 */
       const int ,		/* v - variable index                        */
       const int ,		/* k - kind                                  */
       const int ,              /* matIndex - material index                 */
       const char *,		/* name - name (short)                       */
       const char *,		/* unit - units (unused)                     */
       const char *,		/* desc - name (long, unused)                */
       const int);             /* derivative - indicating time derivative   */
       
EXTERN int set_ev_tkud		/* mm_more_utils.c                           */
(struct Results_Description *, /* r - just elem results for Exodus II  */
       const int ,		/* i - index                                 */
       const int ,		/* v - variable index                        */
       const char *,		/* name - name (short)                       */
       const char *,		/* unit - units (unused)                     */
       const char *,		/* desc - name (long, unused)                */
       const int);             /* derivative - indicating time derivative   */

EXTERN int load_global_var_info	/* mm_more_utils.c                           */
(struct Results_Description *, /* r - global results for Exodus II    */
       const int ,		/* i - index                                 */
       const char *);		/* name - name (short)                       */

EXTERN void sum_total_stress	/* mm_more_utils.c                           */
(double [],		/* sol_vec                                   */
      int ,			/* var_no                                    */
      int ,			/* k                                         */
      double [],		/* nodal_vec                                 */
      Exo_DB *);		/* exo                                       */

EXTERN void extract_nodal_vec	/* mm_more_utils.c                           */
(double [],		/* sol_vec                                   */
      int ,			/* var_no                                    */
      int ,			/* k                                         */
      int ,                     /* matIndex                                  */
      double [],		/* nodal_vec                                 */
      Exo_DB *,  		/* exo                                       */
      int,                      /* timeDerivative                            */
       dbl );                  /* current time                              */

EXTERN void extract_nodal_eb_vec /* mm_more_utils.c                          */
(double [],		/* sol_vec                                   */
       int ,			/* var_no                                    */
       int ,			/* ktype                                     */
       int ,                    /* matIndex                                  */
       int ,                    /* eb_index                                  */ 
       double [],		/* nodal_vec                                 */
       Exo_DB *,	        /* exo                                       */
       int,	                /* timeDeriviative                           */
       double );	        /* time                                      */

EXTERN void extract_elem_vec	/* mm_more_utils.c                           */
(const double [],		/* sol_vec                                   */
       const int ,		/* ev_indx                                   */
       const int ,		/* var_no                                    */
       double ***,		/* gvec_elem                                 */
       const Exo_DB * );	/* exo                                       */

EXTERN void anneal_map		/* mm_more_utils.c                           */
(const int ,		/* dim                                       */
       const double [],		/* X_old                                     */
       const double [],		/* displacement                              */
       double []);		/* X_new                                     */

EXTERN int get_new_coord 
(  double *[],
	  double *,
	  const Exo_DB * );

EXTERN void elements_attached_to_NS
( int *,
		int ,
		Exo_DB * ) ;

#endif /* GOMA_MM_MORE_UTILS_H */
