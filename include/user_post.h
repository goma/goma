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
 

#ifndef _USER_POST_H
#define _USER_POST_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _USER_POST_C
#define EXTERN
#
#endif

#ifndef _USER_POST_C
#define EXTERN extern
#endif

EXTERN double user_post
PROTO((dbl *));			/* param - ptr to user-defined list          */

EXTERN int usr_ptracking
PROTO((FILE * ,			/*  jfp - filename for output                */
       const int ,		/*  part_id - particle id - starts at 1      */
       const double [],		/*  part_x - current particle coords         */
       const double [],		/*  part_v - current particle velocity       */
       const double [],		/*  part_xold - past coords                  */
       const double [],		/*  part_vold - past velocity                */
       const int ,		/*  heading - flag for writing headings      */
       const double ,		/*  time_value - porticle time               */
       const double ));		/*  time_step - time step                    */

#endif /* _USER_POST_H */
