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

/* rd_in.h - prototype declarations for rd_in.c
 */

#ifndef GOMA_RD_IN_H
#define GOMA_RD_IN_H

#include "brkfix/brkfix_types.h"
#include "exo_struct.h"
#include "wr_side_data.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_RD_IN_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_RD_IN_C
#define EXTERN extern
#endif

EXTERN void rd_input
(char *,			/* filename for input */
       Exo_DB *,		/* EXODUS II finite element database */
       Bevm ****,		/* basic eqnvar multiplicity */
       int ****,		/* eqnvar dependencies */
       int ****,		/* local node/dof existences */
       int **);		/* number of basic eqnvars in ea eb */

#endif /* GOMA_RD_IN_H */



