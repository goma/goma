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
 
/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 * $RCSfile: mm_chemkin.h,v $
 * $Author: prschun $
 * $Date: 2007-09-18 18:53:42 $
 * $Revision: 5.1 $
 * $Name: not supported by cvs2svn $
 *====================================================================*/

#ifndef _MM_CHEMKIN_H
#define _MM_CHEMKIN_H

#ifdef USE_CHEMKIN
#include "cpc_defs.h"
#include "ck_chemkin_const.h"
#endif

extern int Chemkin_Needed;

#ifdef USE_CHEMKIN
extern int ck_decide_vol_chem    (CPC_VOLDOMAIN_STRUCT *);
extern int chemkin_mat_prop_init (MATRL_PROP_STRUCT *, int,
				 	 PROBLEM_DESCRIPTION_STRUCT *);
extern void chemkin_initialize_mp (void);

/*
 * externals in the mm_placid.c file
 */
extern int placid(int, int, double, double[], double, double , 
                  double [], int, double, double, double [], 
                  CK_NAME[]);

/*
 *   Defined constants for the ifunc parameter to placid
 *     -> See comments in mm_placid.c for their definitions
 */
#define SFLUX_INITIALIZE 1
#define SFLUX_RESIDUAL   2
#define SFLUX_JACOBIAN   3
#define SFLUX_TRANSIENT  4

/*
 *   Defined constants for the bulkFunc parameter to placid
 */
#define BULK_DEPOSITION  1
#define BULK_ETCH        2

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

#else
extern void chemkin_not_linked (char *errString);
#endif


#endif

