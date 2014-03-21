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


#ifndef _WR_SIDE_DATA_H
#define _WR_SIDE_DATA_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _WR_SIDE_DATA_C
#define EXTERN
#
#endif

#ifndef _WR_SIDE_DATA_C
#define EXTERN extern
#endif

EXTERN int ns_data_print
PROTO((pp_Data *,		/* post processing information */
       double [],		/* solution vector */
       const Exo_DB *,		/* handle to EXODUS II info */
       const double ,		/* current time */
       const double ));		/* current time step size */

EXTERN int ns_data_sens_print
PROTO(( const struct Post_Processing_Data_Sens *,
        const double [],	/* solution vector */
        double** ,		/* sensitivity vector */
        const double ));	/* current time */

EXTERN int match_nsid
PROTO(( int ));

EXTERN int psid2nn
PROTO(( int ));

EXTERN int nsid2nn
PROTO(( int ));

#endif /* _WR_SIDE_DATA_H */
