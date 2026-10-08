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

#ifndef GOMA_WR_SOLN_H
#define GOMA_WR_SOLN_H

#include "dpi.h"
#include "exo_struct.h"
#include "rf_io_structs.h"
#include "std.h"

struct Results_Description;

extern void write_solution(char output_file[],
                           double resid_vector[],
                           double x[],
                           double **x_sens_p,
                           double x_old[],
                           double xdot[],
                           double xdot_old[],
                           int tev,
                           int tev_post,
                           double *gv,
                           struct Results_Description *rd,
                           double *gvec,
                           double ***gvec_elem,
                           int *nprint,
                           dbl delta_t,
                           dbl theta,
                           dbl time_value,
                           dbl *x_pp,
                           Exo_DB *exo,
                           Dpi *dpi); /* ptr to mesh */

extern void write_solution_segregated(char output_file[],
                                      double **resid_vector,
                                      double **x,
                                      double **x_old,
                                      double **xdot,
                                      double **xdot_old,
                                      int *tev_post,
                                      double *gv,
                                      struct Results_Description **rd,
                                      double **gvec,
                                      double ****gvec_elem,
                                      int *nprint,
                                      dbl delta_t,
                                      dbl theta,
                                      dbl time_value,
                                      dbl *x_pp,
                                      Exo_DB *exo,
                                      Dpi *dpi);

#endif
/*******************************************************************************/
/*  END of file wr_soln.h  */
/*******************************************************************************/
