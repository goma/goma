/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2025 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/

#ifndef GOMA_TABLE_H
#define GOMA_TABLE_H

#include "el_elm.h"
#include "rf_bc.h"
#include "rf_bc_const.h"
#include "rf_fem_const.h"

void apply_table_mp(double *func, struct Data_Table *table);

double interpolate_table(struct Data_Table *table, double x[], double *sloper, double dfunc_dx[]);

double
table_distance_search(struct Data_Table *table, double x[], double *sloper, double dfunc_dx[]);

void apply_table_bc(double *func,
                    double d_func[MAX_VARIABLE_TYPES + MAX_CONC],
                    struct Boundary_Condition *BC_Type,
                    double time_value);

double interpolate_table_sat(struct Data_Table *table, double x[DIM]);

#endif /* GOMA_TABLE_H */