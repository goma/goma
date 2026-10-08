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

#ifndef GOMA_MM_EM_BC_H
#define GOMA_MM_EM_BC_H

#include "el_elm.h"
#include "std.h"

int apply_em_farfield_direct_vec                      /* mm_fill_em.c                           */
    (double[DIM],                                     // func
     double[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], // d_func
     double[DIM],                                     // xi
     const int,                                       // bc_name
     double *);

int apply_em_sommerfeld_vec                           /* mm_fill_em.c                           */
    (double[DIM],                                     // func
     double[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], // d_func
     double[DIM],                                     // xi
     const int,                                       // bc_name
     double *);

int apply_em_free_vec                                 /* mm_fill_em.c                           */
    (double[DIM],                                     // func
     double[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], // d_func
     double[DIM],                                     // xi
     const int);                                      // bc_name

int apply_ewave_planewave_vec                         /* mm_fill_em.c                           */
    (double[DIM],                                     // func
     double[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], // d_func
     double[DIM],                                     // xi
     const int,                                       // bc_name
     double *);

int apply_ewave_curlcurl_farfield_vec                 /* mm_fill_em.c                           */
    (double[DIM],                                     // func
     double[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], // d_func
     double[DIM],                                     // xi
     double time,                                     // present time
     const int,                                       // bc_name
     double *);

int apply_ewave_2D                                    /* mm_fill_em.c                           */
    (double[DIM],                                     // func
     double[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], // d_func
     double[DIM],                                     // xi
     const int);                                      // bc_name

int apply_ewave_nedelec_farfield(double func[DIM],
                                 double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                                 double xi[DIM], /* Local stu coordinates */
                                 double time,    // present time
                                 const int bc_name,
                                 double *bc_data);

void em_absorbing_bc_nedelec(int bc_name,
                             dbl *func,
                             double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE]);

void em_mms_nedelec_bc(int bc_name,
                       dbl *func,
                       double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE]);

#endif // GOMA_MM_EM_BC_H
