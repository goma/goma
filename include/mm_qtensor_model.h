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
 


#ifndef GOMA_MM_QTENSOR_MODEL_H
#define GOMA_MM_QTENSOR_MODEL_H

#include "el_elm.h"
#include "mm_std_models.h"
#include "std.h"

struct Species_Conservation_Terms;
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MM_QTENSOR_MODEL_C
#define EXTERN
#endif

#ifndef GOMA_MM_QTENSOR_MODEL_C
#define EXTERN extern
#endif

EXTERN int assemble_vorticity_direction(void);

EXTERN int hydro_qtensor_flux	/* mm_qtensor_model.c */
(struct Species_Conservation_Terms *, /* st */
       int );			/* w - species number */

EXTERN int hydro_qtensor_flux_new /* mm_qtensor_model.c */
(struct Species_Conservation_Terms *, /* st */
       int );			/* w - species number */

EXTERN dbl normalize_really_simple_vector /* mm_qtensor_model.c */
(dbl *,			/* vector to normalize */
       int );			/* dimension */

EXTERN dbl really_simple_vector_magnitude /*mm_qtensor_model.c */
( dbl *,
	int );

EXTERN void cross_really_simple_vectors	/* mm_qtensor_model.c */
(const dbl *,		/* v1 */
       const dbl *,		/* v2 */
       dbl *);			/* v3 = v2 x v1 */

EXTERN void compute_principle_directions /* mm_qtensor_model.c */
(dbl *,			/* flow direction vector */
       dbl *,			/* normal direction vector */
       dbl *,			/* vorticity direction vector */
       int );			/* print toggle */

EXTERN void find_eigenvector	/* mm_qtensor_model.c */
(dbl [3][3],			/* 3 x 3 matrix */
       dbl ,			/* eigenvalue */
       dbl *,			/* eigenvector */
       int );			/* print toggle */

EXTERN void find_super_special_eigenvector(dbl [DIM][DIM], /* mm_qtensor_model.c, tensor to diagonalize */
					   dbl *,  /* eigenvector of vorticity */
					   dbl *,   /* eigenvector corresponding to largest eigenvalue */
					   dbl *,   /* eigenvector corresponding to smallest eigenvalue */
					   dbl *,   /* eigenvector corresponding to middle eigenvalue */
					   dbl *,  /* eigenvalue of vorticity direction */
					   int);   /* print toggle  */

EXTERN void find_eigenvalues_eigenvectors(dbl [3][3], dbl *, dbl *, dbl *,
					  dbl *, dbl *, dbl *);  /* mm_qtensor_model.c */

EXTERN void diagonalize_symmetric_tensor /* mm_qtensor_model.c */
(dbl [3][3],		/* tensor to diagonalize. */
       dbl *,			/* flow direction vector. */
       dbl *,			/* norm direction vector. */
       dbl *,			/* vorticity direction vector. */
       dbl *,			/* eigenvalues */
       int );			/* print toggle */

EXTERN void diagonalize_rate_of_deformation_tensor /* mm_qtensor_model.c */
(dbl [3][3],		/* tensor to diagonalize. */
       dbl *,			/* flow direction vector. */
       dbl *,			/* norm direction vector. */
       dbl *,			/* vorticity direction vector. */
       dbl *,			/* eigenvalues */
       int );			/* print toggle */

EXTERN void compute_VQVt_directly /* mm_qtensor_model.c */
(dbl [3][3],		/* rate of deformation tensor */
       dbl [3][3],		/* VQVt */
       int );			/* print toggle */

EXTERN void please_work		/* mm_qtensor_model.c */
(dbl [3][3],		/* rate of deformation tensor */
       dbl *,			/* flow direction vector */
       dbl *,			/* normal direction vector */
       dbl *,			/* vorticity direction vector */
       int );			/* print toggle */

EXTERN void get_characteristic_eq_coeffs
(dbl [3][3],		/* E, rate of deformation tensor */
       dbl *,			/* a0 */
       dbl *,			/* a1 */
       dbl *);			/* a2 */

EXTERN int bias_eigenvector_to(dbl *,    /* eigenvector  */
			       dbl *);    /* reference vector */

EXTERN void assemble_qtensor
(dbl *);

EXTERN void assemble_new_qtensor /* Ryan's qtensor */
(dbl *);


extern int MMH_ip;

/*don't ask about this skeleton. I need the local material-referenced element number deep in the bowels of the assembly routines, and passing exo struct down there is brutal. */
int PRS_mat_ielem; 

#define QTENSOR_SMALL_DBL 1.0e-14

#define MAGIC_VECTOR_0 1.0
#define MAGIC_VECTOR_1 0.0
#define MAGIC_VECTOR_2 0.0
#define BAD_MAGIC_VECTOR_RELATIVE_TOLERANCE 0.01

#define FINITE_DELTA 0.05

extern dbl vort_dir[MDE][DIM];	/* vorticity direction for each gauss point */
extern dbl qtensor[MDE][DIM][DIM]; /* I - 1/2 v^t v for each gauss point*/
extern dbl div_qtensor[MDE][DIM]; /* div(I-1/2 v^tv) for each gauss point */

#endif /* GOMA_MM_QTENSOR_MODEL_H */
