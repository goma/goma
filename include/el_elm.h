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
 

/*
 *$Id: el_elm.h,v 5.3 2009-11-13 23:20:07 prschun Exp $
 */

#ifndef _EL_ELM_H
#define _EL_ELM_H

#include "rf_fem_const.h"

typedef
enum type_elem {BILINEAR_QUAD = 0,
		C_BILINEAR_QUAD, /* 5-node quad.  4 corner nodes + centroid node*/
		S_BIQUAD_QUAD,
		BIQUAD_QUAD,
		BIQUAD_QUAD_LS,	/* uses higher order integration */
		P1_QUAD,	/* added 94/03/15 pas */
		P0_QUAD,	/* added 98/02/09 rrr */

		TRILINEAR_HEX,
		C_TRILINEAR_HEX, /* 9 node hex.  8 corner nodes + centroid node */
		S_TRIQUAD_HEX,
		TRIQUAD_HEX,
		P1_HEX,
		P0_HEX,

		LINEAR_BAR,     /* 1D, linear elements for shell elements*/
		QUAD_BAR,       /* 1D, quadratic elements for shell elements*/

                LINEAR_TRI,
                QUAD_TRI,
                QUAD6_TRI,      /* quadratic triangle with 6th order quadrature */

		TRILINEAR_TET,   /* 3D linear tetrahedral element type */

		BILINEAR_SHELL,  /* 2D, linear elements for 3D shells */
		BIQUAD_SHELL,    /* 2D, quadratic elements for 3D shells */
		BILINEAR_TRISHELL,  /* 2D, linear triangular elements for 3D shells */
                P1_SHELL,	/* 2D, discontinuous linear elements for 3D shells  */
                P0_SHELL	/* 2D, discontinuous constant elements for 3D shells */
}
Type_Elem;

/*
 * Define element shapes...
 * Element_Shape	LINE_SEGMENT	(1D)
 *
 *		  	TRIANGLE	(2D)
 *		  	QUADRILATERAL	(2D)
 *
 *		  	TETRAHEDRON	(3D)
 *		  	PRISM		(3D) (wedge)
 *		  	PYRAMID		(3D) (square bottomed)
 *		  	HEXAHEDRON	(3D)
 *
 * Added: 	Fri Dec 17 12:43:18 MST 1993 pasacki@sandia.gov
 */

#define LINE_SEGMENT	0
#define TRIANGLE	1
#define QUADRILATERAL	2
#define TETRAHEDRON	3
#define PRISM		4
#define HEXAHEDRON	5
#define PYRAMID		6
#define SHELL           7
#define TRISHELL        8



/* define element data "request for information" types */

#define NNODES            1
#define NQUAD             2
#define NDIM              3
#define NQUAD_SURF        4
#define NQUAD_EDGE        5

/* define shape function information types */

#define PSI               0
#define DPSI_S            1
#define DPSI_T            2
#define DPSI_U            3

/* define maximum quantities */

#ifndef MAX_SUR_ELEM_2D
#define MAX_SUR_ELEM_2D  8   /* Maximum number of elements            */
                             /* surrounding (containing )a given node */
#endif

#ifndef AVG_SUR_ELEM_2D
#define AVG_SUR_ELEM_2D  4   /* A more realistic estimate of number of elements */
                             /* surrounding (containing )a given node */
#endif

#ifndef MAX_SUR_ELEM_3D
#define MAX_SUR_ELEM_3D  50  /* A more realistic estimate of number of elements */
#endif
                             /* surrounding (containing )a given node */
#ifndef AVG_SUR_ELEM_3D
#define AVG_SUR_ELEM_3D  8   /* Maximum number of elements            */
                             /* surrounding (containing )a given node */
#endif

#ifndef MAX_SURF_GP
#define MAX_SURF_GP       9   /* Maximum number of surface quad points  */
#endif
                          
#define OUTER_SPACE	(-1)	/* Special value for element id flags */
				/* "non elements" for the element-element */
				/* connectivity. Later, other values will */
				/* indicate "element on another processor" */

#define ANOTHER_PROC	(-2)	/* here it is. */

#define UNASSIGNED_YET	(-3)	/* for the initial fill */

#define MAX_EPN		MAX_SUR_ELEM_3D	/* convenient alias for maximum
					 * number of elements per node */

#ifndef MAX_NODES_PER_SIDE	
#define MAX_NODES_PER_SIDE	(9)
#endif

#ifndef MAX_PDIM
#define MAX_PDIM         3   /* Maximum physical problem dimension    */
#endif

/*
 * This constant is so widely used that a shorter alias is also available...
 */

#ifndef DIM
#define DIM		MAX_PDIM
#endif

char   element_type[ 32 + 1]; /* Hardwired to be equivalent to MAX_STR_LENGTH
                                 in exodusII.h */
int *listel;                  /* Pointer to element order map from exoII  */

extern int ID_Side_Quad[];	/* defined in exo_conn.c */

extern int ID_Side_Hex[];	/* defined in exo_conn.c */

/******************************* PROTOTYPES FOR el_elm_info.c ****************/
#if 0				/* died 1999/10/04 */
extern int
elem_info
   PROTO (( int,                        /* info  */
	    int  ));                    /* ielem_type  */

extern int
imax_quad
   PROTO (( int,                        /*  n  */
	    int [] ));                  /* elem_type[]  */

extern void
find_stu
   PROTO (( int,                        /* iquad  */
	    int,                        /* ielem_type  */
	    double *,                   /* *s  */
            double *,                   /* *t  */
            double *  ));               /* *u  */

extern void
find_surf_st
   PROTO (( int,                        /* iquad  */
	    int,                        /* ielem_type  */ 
	    int,                        /* iside  */
	    int,                        /* dim  */
            double [],                  /* xi[MAX_PDIM]  */
	    double *,                   /* *s  */
            double *  ));               /* *t  */

extern void
find_nodal_stu
   PROTO (( int,                        /* inode   */
	    int,                        /* ielem_type  */ 
	    double *,                   /* *s   */
            double *,                   /* *t   */
            double * ));                /* *u   */

extern double
Gq_weight
   PROTO (( int,                        /* iquad  */
	    int  ));                    /* ielem_type */

extern double
Gq_surf_weight
   PROTO (( int,                        /* iquad  */
	    int  ));                    /* ielem_type */

extern int
in_list
   PROTO (( int,                        /* ivalue  */
	    int,                        /* ibegin  */
            int,                        /* iend  */
            int [] ));                  /* ivector[]  */

extern int
get_type
   PROTO (( char [],                    /* string[]  */
	    int  ));                    /* nodes  */

extern double
shape
   PROTO (( double,                     /*  s   */
            double,                     /*  t   */
            double,                     /*  u   */
	    int,                        /*  Ielem_type  */ 
	    int,                        /*  Iquant  */
	    int  ));                    /*  Inode   */

extern int
type2shape
   PROTO (( int ));			/* element_type */

extern int
shape2sides
   PROTO (( int ));			/* element_shape */
#endif /* #if 0 since 1999/10/04 ... */

#endif

