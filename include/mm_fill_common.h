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
 

#ifndef GOMA_MM_FILL_COMMON_H
#define GOMA_MM_FILL_COMMON_H

#include "mm_std_models.h"
#include "std.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MM_FILL_COMMON_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_MM_FILL_COMMON_C
#define EXTERN extern
#endif


EXTERN void computeCommonMaterialProps_gp
(
       const dbl time  //Current time (sec)
);


#endif 
