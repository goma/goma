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
 
#ifndef _BC_DIRICH_H
#define _BC_DIRICH_H


#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _BC_DIRICH_C
#define EXTERN /* do nothing */
#endif

#ifndef _BC_DIRICH_C
#define EXTERN extern
#endif

EXTERN int put_dirichlet_in_matrix
(double [],		/* x - Solution vector                       */
       const int);		/* num_total_nodes */


#endif /* _BC_DIRICH_H */
