/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2023 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/
#ifndef GOMA_MM_FILL_CONTINUITY_H
#define GOMA_MM_FILL_CONTINUITY_H

#include "el_elm.h"
#include "exo_struct.h"
#include "mm_as_structs.h"
#include "mm_fill_common.h"
#include "mm_fill_energy.h"
#include "rf_fem_const.h"
#include "std.h"
#include <stdbool.h>
int assemble_continuity /* mm_fill_terms.c                           */
    (dbl,               /* time_value */
     dbl,               /* tt - to vary time integration from
                           explicit (tt = 1) to implicit (tt = 0)    */
     dbl,               /* dt - current time step size               */
     const PG_DATA *);
double FoamVolumeSource(double,
                        double,
                        double,
                        double[DIM][MDE],
                        double[MDE],
                        double[DIM][MDE],
                        double[MAX_CONC][MDE],
                        double[MDE]);
double REFVolumeSource(
    double, double, double, double[DIM][MDE], double[MDE], double[DIM][MDE], double[MAX_CONC][MDE]);
int assemble_projection_time_stabilization(Exo_DB *exo, double time, double tt, double dt);
int assemble_projection_stabilization(Exo_DB *, double);
int assemble_PPPS_generalized(Exo_DB *);
#endif // GOMA_MM_FILL_CONTINUITY_H
