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
 *$Id: mm_post_proc.h,v 5.1 2007-09-18 18:53:45 prschun Exp $
 */

/* 
 * This include file declares and transfers options
 * to the post_processing routines.
 */

#ifndef GOMA_MM_POST_PROC_H
#define GOMA_MM_POST_PROC_H

#include <stdio.h>

#include "dpi.h"
#include "exo_struct.h"
#include "mm_more_utils.h"
#include "rf_io_structs.h"
#include "std.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MM_POST_PROC_C
#define EXTERN
#endif

#ifndef GOMA_MM_POST_PROC_C
#define EXTERN extern
#endif

/*
 * Prototypes for functions in mm_post_proc.c
 */
EXTERN void post_process_global(double *x,	 /* Solution vector for the current processor */
				Exo_DB *exo,
				Dpi *dpi,
				double time);

EXTERN void post_process_nodal  /* mm_post_proc.c                            */
(double [],               /* x                                         */
       double **,               /* x_sens_p                                  */
       double [],               /* x_old                                     */
       double [],               /* xdot                                      */
       double [],               /* xdot_old                                  */
       double [],               /* resid_vector                              */
       int ,                    /* ts                                        */
       double *,                /* time                                      */
       double ,                 /* delta_t                                   */
       double,                  /* theta                                     */
       double *,                /* x_pp - for retrieving post-process vars   */
       Exo_DB *,                /* exo - ptr to fem db                       */
       Dpi *,                   /* dpi - ptr to distrib inf                  */
       RESULTS_DESCRIPTION_STRUCT *,  /* exodus description of variables   */
       char [],  /* exodus filename   */
       int ); /* matrix offset */

EXTERN void post_process_elem   /* mm_post_proc.c                            */
(double [],               /* x - soln vector                           */
       double [],               /* x_old - soln vector at previous time step */
       double [],               /* xdot - time derivative of soln vector     */
       double [],               /* xdot_old                                  */
       double [],               /* resid_vector - Residual vector            */
       const int ,              /* tev                                       */
       const int ,              /* tev_post                                  */
       double ***,              /* gvec_elem - Triply indexed array          *
                                 * containing element variable values on     *
                                 * return. Convention:                       *
                                 * [elemblock_index][elemvar_index]          *
                                 *      [element_index(inblock)]             */
       const int ,              /* ts                                        */
       const double *,          /* time_ptr                                  */
       const double ,           /* delta_t                                   */
       Exo_DB * const,          /* exo                                       */
       Dpi * const,             /* dpi                                       */
       struct Results_Description *);

EXTERN void rd_post_process_specs /* mm_post_proc.c                          */
(FILE *,                  /* ifp - input file pointer (strm to input)  */
       char *);                /* input - latest buffer of read values      */

EXTERN int load_nodal_tkn       /* mm_post_proc.c                            */
(struct Results_Description *, /* rd                                   */
       int *,                   /* tnv                                       */
       int *);                 /* tnv_post                                  */

EXTERN int load_elem_tkn        /* mm_post_proc.c                            */
(struct Results_Description *, /* rd                                   */
       const Exo_DB *,           /* Exodus II database struct */
      int ,                     /* tev                                       */
      int *);                  /* tev_post                                  */

EXTERN int find_id_edge         /* mm_post_proc.c */
(const int ,              /* ielem */
       const int ,              /* num_nodes_on_edge */
       const int [],            /* local_edge_node_list */
       int [],                  /* id_local_elem_coord */
       int *,                   /* param_dir - direction of parametric
                                 * edge curve */
       const Exo_DB *);        /* exo */

EXTERN int find_id_edge_TET         /* mm_post_proc.c */
(const int ,              /* ielem */
       const int ,              /* num_nodes_on_edge */
       const int [],            /* local_edge_node_list */
       int [],                  /* id_local_elem_coord */
       int *,                   /* param_dir - direction of parametric
                                 * edge curve */
       const Exo_DB *);        /* exo */

EXTERN int count_nodes_on_SS    /* mm_post_proc.c                            */
(const int ,              /* ss_id - SS id of Primary Side Set         */
       const int ,              /* ss_id2 - SS id of 2nd Side Set for edges  */
       const int ,              /* ss_id3 - SS id of 3rd Side Set for vtces  */
       const int ,              /* iconnect_ptr                              */
       const int ,              /* ielem                                     */
       const int ,              /* num_local_nodes                           */
       int [MDE],               /* local_ss_node_list                        */
       int [MDE]);             /* local_elem_node_id                        */

extern int find_id_elem		/* mm_post_proc_util.c */
(const dbl ,		/* x_coordinate */
       const dbl ,		/* y_coordinate */
       const dbl ,		/* z_coordinate */
       dbl [],          /*  solution vector  */
       const Exo_DB *,
       const int,      /*Starting element for contiguous search */
       const int);    /*ending element for contiguous search */

extern int invert_isoparametric_map   /*mm_post_proc_util.c */
( int *,              /* Current element id */
        const dbl [],              /* x_coordinate */
        dbl [],                  /* s isoparametric coordinate (output)*/
        const Exo_DB *,
        dbl [],                 /*  x - solution vector  */
	int *);		/*  velocity basis fcns  */


extern int elem_order_for_nodal_connect(int *, const Exo_DB *);
extern int check_elem_order(const int *, const Exo_DB *);
#endif
