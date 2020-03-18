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
 * $Id: wr_side_data.h,v 5.1 2007-09-18 18:53:49 prschun Exp $
 */


#ifndef GOMA_WR_SIDE_DATA_H
#define GOMA_WR_SIDE_DATA_H

#include "exo_struct.h"
#include "mm_post_def.h"
#include "wr_exo.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_WR_SIDE_DATA_C
#define EXTERN
#
#endif

#ifndef GOMA_WR_SIDE_DATA_C
#define EXTERN extern
#endif

EXTERN int ns_data_print
(pp_Data *,		/* post processing information */
       double [],		/* solution vector */
       const Exo_DB *,		/* handle to EXODUS II info */
       const double ,		/* current time */
       const double );		/* current time step size */

EXTERN int ns_data_sens_print
( const struct Post_Processing_Data_Sens *,
        const double [],	/* solution vector */
        double** ,		/* sensitivity vector */
        const double );	/* current time */

EXTERN int match_nsid
( int );

EXTERN int psid2nn
( int );

EXTERN int nsid2nn
( int );

#endif /* GOMA_WR_SIDE_DATA_H */
