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
 

#ifndef _RF_PRE_PROC_H
#define _RF_PRE_PROC_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _RF_PRE_PROC_C
#define EXTERN
#
#endif

#ifndef _RF_PRE_PROC_C
#define EXTERN extern
#endif

EXTERN void pre_process		/* rf_pre_proc.c */
PROTO((Exo_DB *exo));

#endif /* _RF_PRE_PROC_H */
