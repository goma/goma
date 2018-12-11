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
 

/* Header file for most of the structures and #define'ed parameters
 * for the suite of sl_eggroll*.c source files.  These files constitute a
 * generalized eigenvalue problem solver with reverse communication.
 *
 * Originally written by Ian Gates.
 *
 * Modification history:
 *   - July 23, 1997, first version.
 *   - Jan 12, 2000, MMH rearranging and cleaning.
 */

#ifndef _SL_EGGROLL_DEF_H
#define _SL_EGGROLL_DEF_H

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define IDG_EPS    1.0e-14
#define IDG_BIG    1.0e+14
#define IDG_STEP   1.0e-08

#define PARAMETER_ARRAY_LENGTH 20

/* Eigenstuff structure.
 */
typedef struct
{ 
  int n; 
  dbl *e;
  dbl *r;
  dbl *i;
  dbl *x;
  dbl **u;
} EV;

/* Postscript structure
 */
typedef struct 
{
  dbl unused;
  dbl xmin;
  dbl xmax;
  dbl zmin;
  dbl zmax;
  dbl xsc;
  dbl zsc;
  dbl rotn;
  dbl pbox_x1;
  dbl pbox_x2;
  dbl pbox_z1;
  dbl pbox_z2;
  dbl trans_x;
  dbl trans_z;
  dbl bbox_x1;
  dbl bbox_x2;
  dbl bbox_z1;
  dbl bbox_z2;
} PS;

#endif
