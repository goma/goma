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
 

#ifndef GOMA_RF_PRE_PROC_H
#define GOMA_RF_PRE_PROC_H

#include "ac_particles.h"
#include "exo_struct.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_RF_PRE_PROC_C
#define EXTERN
#
#endif

#ifndef GOMA_RF_PRE_PROC_C
#define EXTERN extern
#endif

EXTERN void pre_process		/* rf_pre_proc.c */
(Exo_DB *exo);

#endif /* GOMA_RF_PRE_PROC_H */
