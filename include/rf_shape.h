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
 

#ifndef GOMA_RF_SHAPE_H
#define GOMA_RF_SHAPE_H

#include "rf_pre_proc.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_RF_SHAPE_C
#define EXTERN
#
#endif

#ifndef GOMA_RF_SHAPE_C
#define EXTERN extern
#endif

EXTERN double shape		/* rf_shape.c                                */
(const double ,		/* s - quadrature point coordinates          */
       const double ,		/* t                                         */
       const double ,		/* u                                         */
       const int ,		/* Ielem_type - element type                 */
       const int ,		/* Iquant - desired quantity (phi, phi_s, ...*/
       const int );		/* Inode - current element node              */

#endif /* GOMA_RF_SHAPE_H */
