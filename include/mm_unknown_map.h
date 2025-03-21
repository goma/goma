/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2022 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/

/*
 * mm_unknown_map.h -- prototype declarations for mm_unknown_map.c
 */

/*
 * $Id: mm_unknown_map.h,v 5.1 2007-09-18 18:53:46 prschun Exp $
 */

#ifndef GOMA_MM_UNKNOWN_MAP_H
#define GOMA_MM_UNKNOWN_MAP_H

#include "dp_types.h"
#include "dpi.h"
#include "exo_struct.h"
#include "rf_vars_const.h"

extern void setup_local_nodal_vars(Exo_DB *, /* exo - ptr to FE EXODUS II database        */
                                   Dpi *);   /* dpi - ptr to parallel info                */

extern void setup_external_nodal_vars(Exo_DB *, /* exo - ptr to FE EXODUS II database        */
                                      Dpi *,    /* dpi - ptr to parallel info                */
                                      Comm_Ex **);

extern int find_MaxUnknownNode(void);

extern void print_vars_at_nodes(void);

extern void set_unknown_map(Exo_DB *, /* exo - ptr to FE EXODUS II database        */
                            Dpi *);   /* dpi - ptr to parallel info                */

extern int Index_Solution(const int nodeNum,     /* Global Node Number                        */
                          const int varType,     /* Variable Type                             */
                          const int subvarIndex, /* Subvar index                              */
                          const int iNdof,       /* Local nodal degree of freedom. This is    *
                                                  * equal to zero, except for centroid        *
                                                  * pressures, and hermite cubic interpolated *
                                                  * variables.                                */
                          const int matID,       /* material id -> or -1 if the specification *
                                                  * of the material doesn't matter in this    *
                                                  * instance. It would matter if the soln     *
                                                  * number at this node varied depending on   *
                                                  * the material id                           */
                          const int imtrx        /* Matrix ID */
);

extern int variable_type_nodalInterp(int);

extern VARIABLE_DESCRIPTION_STRUCT *Index_Solution_Inv(const int, int *, int *, int *, int *, int);

extern void dofname40(const int, char *);

#endif /* __MM_UNKNOWN_MAP_H */
