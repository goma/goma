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
 
#ifndef _DP_VIF_H
#define _DP_VIF_H


#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _DP_VIF_C
#define EXTERN /* do nothing */
#endif

#ifndef _DP_VIF_C
#define EXTERN extern
#endif

EXTERN void noahs_raven			/* dp_vif.c */
PROTO((void ));

EXTERN void noahs_ark			/* dp_vif.c */
PROTO((void ));

EXTERN void noahs_dove			/* dp_vif.c */
PROTO((void ));

EXTERN void raven_landing		/* dp_vif.c */
PROTO((void ));

EXTERN void ark_landing			/* dp_vif.c */
PROTO((void ));

#endif /* _DP_VIF_H */
