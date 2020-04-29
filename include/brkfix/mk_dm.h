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

/* mk_dm.h - prototype declarations for mk_dm.c
 */

#ifndef GOMA_MK_DM_H
#define GOMA_MK_DM_H

#include "brkfix/brkfix_types.h"
#include "brkfix/nodesc.h"
#include "exo_struct.h"
#include "wr_side_data.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MK_DM_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_MK_DM_C
#define EXTERN extern
#endif

EXTERN void make_goma_dofmap /* mk_dm.c */
    (Exo_DB *x, Bevm ***mult, int ***Lucky, int *num_basic_eqnvars, int *node_kind, int *node_dof0,
     Node_Description **pnd, int *nkn);			/* nkn - actual number of node descriptions */

#endif /* GOMA_MK_DM_H */
