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

EXTERN void make_goma_dofmap	/* mk_dm.c */
(Exo_DB *,		/* full mesh description */
       Bevm ***,		/* mult - basic eqnvar multiplicities */
       int ***,			/* evd - eqnvar dependencies */
       int ***,			/* Lucky - local nodal dof existence profs */
       int *,			/* num_basic_eqnvars - for ea elemblock */
       int *,			/* node_kind - list for every global node */
       int *,			/* node_dof0 - name of first dof for ea node */
       Node_Description **,	/* pnd - pointer array to node descriptions */
       int *);			/* nkn - actual number of node descriptions */

#endif /* GOMA_MK_DM_H */
