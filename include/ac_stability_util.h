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

#endif /* _AC_STABILITY_UTIL_H */
