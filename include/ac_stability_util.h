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

#ifndef GOMA_AC_STABILITY_UTIL_H
#define GOMA_AC_STABILITY_UTIL_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_AC_STABILITY_UTIL_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_AC_STABILITY_UTIL_C
#define EXTERN extern
#endif

#include "ac_stability.h"
#include "ac_update_parameter.h"
#include "dp_types.h"
#include "dpi.h"
#include "exo_struct.h"
#include "rf_io_structs.h"

EXTERN void do_LSA_mods(int);

EXTERN void modify_basis_and_weight_functions_for_LSA_3D_of_2D(void);

EXTERN void modify_bf_mesh_derivs_for_LSA_3D_of_2D(void);

EXTERN void modify_fv_mesh_derivs_for_LSA_3D_of_2D(void);

EXTERN void modify_normal_vector_for_LSA_3D_of_2D(void);

EXTERN int create_eigen_outfiles(Exo_DB *, Dpi *, RESULTS_DESCRIPTION_STRUCT *, double ***);

EXTERN void get_eigen_outfile_name(char *, int, int);

EXTERN int do_loca(Comm_Ex *, Exo_DB *, Dpi *);

EXTERN int anneal_mesh_LSA /* ac_stability_util.c */
    (double[],             /* x - solution vector */
     Exo_DB *,             /* exo - entire mesh desc. */
     double **,            /* Saved mesh coordinates */
     double **);           /* Saved displacement fields */

EXTERN int unanneal_mesh_LSA /* ac_stabililty_util.c */
    (double[],               /* x - solution vector */
     Exo_DB *,               /* exo - entire mesh desc. */
     double **,              /* Saved mesh coordinates */
     double **);             /* Saved displacement fields */

EXTERN void add_displacement_LSA /* ac_stabililty_util.c */
    (double[],                   /* x - eigenvector */
     Exo_DB *,                   /* exo - entire mesh desc. */
     double **);                 /* Saved displacement fields */

EXTERN void undo_add_displacement_LSA /* ac_stabililty_util.c */
    (double[],                        /* x - eigenvector */
     Exo_DB *,                        /* exo - entire mesh desc. */
     double **);                      /* Saved displacement fields */

#endif /* GOMA_AC_STABILITY_UTIL_H */
