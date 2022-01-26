/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2022 Sandia Corporation.                                  *
*                                                                         *
* Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,  *
* the U.S. Government retains certain rights in this software.            *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
\************************************************************************/

/* mm_mp.h -- include file for structure pointing to material properties
 *
 * Created:	Sat Mar 19 15:51:22 MST 1994 pasacki@sandia.gov
 *
 * Revised:
 *
 * $Id: mm_mp.h,v 5.3 2008-01-25 14:58:26 hkmoffa Exp $
 */
#ifndef GOMA_MM_MP_H
#define GOMA_MM_MP_H

#include "mm_mp_structs.h"

extern struct Material_Properties *mp, **mp_glob, *mp_old;

//! Declaration of the global temporary pointer for viscosity related props
/*!
 *  this is used a temporary pointer variable
 */
extern GEN_NEWT_STRUCT *gn;

//! Declaration for a vector of global permament pointers for viscosity related props
/*!
 *   This loops over ??
 *   It has a length of
 */
extern GEN_NEWT_STRUCT **gn_glob;

extern struct Elastic_Constitutive *elc, **elc_glob, *elc_rs, **elc_rs_glob;

extern struct Viscoelastic_Constitutive **ve, ***ve_glob;

extern struct Viscoelastic_Nonmodal *vn, **vn_glob;

extern struct Viscoplastic_Constitutive *evpl, **evpl_glob;

extern struct Variable_Initialization Var_init[MAX_VARIABLE_TYPES + MAX_CONC];

extern struct Variable_Initialization Var_init_mat[MAX_NUMBER_MATLS][MAX_VARIABLE_TYPES + MAX_CONC];

/*
 *  Some Prototypes
 *
 */
/*      mm_matrl.c */

extern void matrl_prop_print(MATRL_PROP_STRUCT *, int);
extern int goma_mat_prop_init(MATRL_PROP_STRUCT *, int, PROBLEM_DESCRIPTION_STRUCT *);
extern double calc_density(MATRL_PROP_STRUCT *, int, PROPERTYJAC_STRUCT *, double time);
extern double calc_concentration(MATRL_PROP_STRUCT *, int, PROPERTYJAC_STRUCT *);
extern void load_properties(MATRL_PROP_STRUCT *, double);

#endif
