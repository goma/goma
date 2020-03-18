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

#include "dp_types.h"
#include "exo_struct.h"
#include "dpi.h"
#include "rf_io_structs.h"
#include "ac_stability.h"

EXTERN void do_LSA_mods
(int);

EXTERN void modify_basis_and_weight_functions_for_LSA_3D_of_2D
( void );

EXTERN void modify_bf_mesh_derivs_for_LSA_3D_of_2D
( void );

EXTERN void modify_fv_mesh_derivs_for_LSA_3D_of_2D
( void );

EXTERN void modify_normal_vector_for_LSA_3D_of_2D
( void );

EXTERN int create_eigen_outfiles
(Exo_DB *, Dpi *, RESULTS_DESCRIPTION_STRUCT *);

EXTERN void get_eigen_outfile_name
(char *, int, int);

EXTERN int do_loca
(Comm_Ex *, Exo_DB *, Dpi *);

#endif /* GOMA_AC_STABILITY_UTIL_H */
