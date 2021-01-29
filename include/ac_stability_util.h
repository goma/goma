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
 

#ifndef _AC_STABILITY_UTIL_H
#define _AC_STABILITY_UTIL_H


#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _AC_STABILITY_UTIL_C
#define EXTERN /* do nothing */
#endif

#ifndef _AC_STABILITY_UTIL_C
#define EXTERN extern
#endif

EXTERN void do_LSA_mods
PROTO((int));

EXTERN void modify_basis_and_weight_functions_for_LSA_3D_of_2D
PROTO(( void ));

EXTERN void modify_bf_mesh_derivs_for_LSA_3D_of_2D
PROTO(( void ));

EXTERN void modify_fv_mesh_derivs_for_LSA_3D_of_2D
PROTO(( void ));

EXTERN void modify_normal_vector_for_LSA_3D_of_2D
PROTO(( void ));

EXTERN int create_eigen_outfiles PROTO((Exo_DB *, Dpi *, RESULTS_DESCRIPTION_STRUCT *, double ***));

EXTERN void get_eigen_outfile_name
PROTO((char *, int, int));

EXTERN int do_loca
PROTO((Comm_Ex *, Exo_DB *, Dpi *));

EXTERN int anneal_mesh_LSA      /* ac_stability_util.c */
PROTO(( double [],              /* x - solution vector */
        Exo_DB *,               /* exo - entire mesh desc. */
        double **,              /* Saved mesh coordinates */
        double **));            /* Saved displacement fields */

EXTERN int unanneal_mesh_LSA    /* ac_stabililty_util.c */
PROTO(( double [],              /* x - solution vector */
        Exo_DB *,               /* exo - entire mesh desc. */
        double **,              /* Saved mesh coordinates */
        double **));            /* Saved displacement fields */

EXTERN void add_displacement_LSA /* ac_stabililty_util.c */
PROTO(( double [],              /* x - eigenvector */
        Exo_DB *,               /* exo - entire mesh desc. */
        double **));            /* Saved displacement fields */

EXTERN void undo_add_displacement_LSA /* ac_stabililty_util.c */
PROTO(( double [],              /* x - eigenvector */
        Exo_DB *,               /* exo - entire mesh desc. */
        double **));            /* Saved displacement fields */

#endif /* _AC_STABILITY_UTIL_H */
