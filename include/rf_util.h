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
 *$Id: rf_util.h,v 5.4 2010-03-09 21:57:55 sarober Exp $
 */

/*
 * Revision history is gone.
 */

#ifndef GOMA_RF_UTIL_H
#define GOMA_RF_UTIL_H

#if 0
extern void scatter_double_vector
(const int ,
       const int [],
       double [],
       const double []);

extern void gather_double_vector
(const int ,		/* n */
       const int [],		/* index */
       const double [],		/* y */
       double []);		/* x */

#endif

extern int ilog2i (unsigned int);

extern int countmap_vardofs	/* rf_util.c */
(const int ,	       	/* variable - VELOCITY1, etc.            (in) */
       const int ,	       	/* num_nodes                             (in) */
       int *);			/* map - map[node_index] = dof_index of variable
				 * at that node                         (out) */

extern void init_vec		/* rf_util.c */
(double [],		/* u - soln vector */
       Comm_Ex *,               /* cx - communications structure */
       Exo_DB *,                /* exo - ptr to finite element mesh database */
       Dpi *,                   /* dpi - distributed processing struct */
       double [],               /* uAC - AC variables */
       int,                     /*  nAC */
       double *timeValueRead); /* timeValueRead time value read from an input file */

extern void fill_dvec_rand
(double [],		/* u   -- some vector of doubles */
       int);			/* len -- length of vector u */

extern int rd_exoII_nv
( double *,
        int ,
	int ,
	MATRL_PROP_STRUCT *,
        char **,
        int ,
        int ,                  /* number of variables in EXODUS II file */
        int ,                  /* ID of the open EXODUS II file */
        int ,                  /* 1-based */
        int );                /* this is zero or the species number */

extern double time_step_control	/* rf_util.c                                 */
(const double ,		/* delta_t_old                               */
       const double ,		/* delta_t_older                             */
       const int ,		/* const_delta_t                             */
       const double [],		/* x                                         */
       const double [],		/* x_pred                                    */
       const double [],		/* x_old                                     */
       const double [],         /* x_AC                                      */
       const double [],         /* x_AC_pred                                 */
       const double ,		/* eps                                       */
       int *,			/* success_dt                                */
       const int []);		/* use_var_norm                              */

extern int dump_stability_matrices /* rf_util.c                              */
(const int [],		/* ija - CMSR column pointers into matrix    */
       const double []);	/* a - Jacobian matrix                       */

extern int filter_conc		/* rf_util.c                                 */
(const int ,		/* N - number of nodes                       */
       dbl [],			/* x - solution vector                       */
       const int ,		/* filter_species_material_number -          */
                                /* 1st material number w/suspension          */
       const dbl ,		/* cmin - lower cutoff                       */
       const dbl );		/* cmax - upper cutoff                       */

#ifndef COUPLED_FILL
extern void get_fill_vector	/* rf_util.c                                 */
(const int ,		/* N                                         */
       const double [],		/* x                                         */
       double [],		/* fill_vector                               */
       const int []);		/* node_to_fill                              */
#endif /* not COUPLED_FILL */
extern void sort2_int_int       /* rf_util.c                                 */
(const int,
       int [],
       int []);

#ifndef COUPLED_FILL
extern void put_fill_vector	/* rf_util.c                                 */
(const int ,		/* N                                         */
       double [],		/* x                                         */
       const double [],		/* fill_vector                               */
       const int []);		/* node_to_fill                              */
#endif /* not COUPLED_FILL */

extern int wr_soln_vec
( double [],		/* u - solution vector                       */
	double [],		/* r - residual vector                       */
	const int ,		/* np - number elements in solution vector   */
	const int );		/* itn - Newton iteration we are at          */

extern int write_ascii_soln
( double *,		/* u */
	double *,              /* Resid */
	int ,			/* np */
	double *,		/* uAC */
	int ,			/* nAC */
	double ,		/* time */
	FILE * );		/* file */

extern double path_step_control	/* rf_util.c                                 */
(int ,			/* N                                         */
       double ,			/* delta_s_old                               */
       double ,			/* delta_s_older                             */
       double [],		/* x                                         */
       double ,			/* eps                                       */
       int *,			/* success_ds                                */
       int [],			/* use_var_norm                              */
       int );			/* inewton                                   */

extern int rd_globals_from_exoII /* rf_util.c */
( double *,                /* field to be read */
        const char *,            /* file containing the field values */
        const int,               /* beginning of global vars */
        const int );            /* number of global vars to read */

extern int build_node_index_var /* rf_util.c                                 */
(const int ,              /* name of variable to be indexed            */
       const int ,              /* upper limit on range of nodes to be       */
                                /* searched in Dolphin array                 */
       const int *,             /* global node index array                   */
       int *,                   /* array of global node numbers              */
       int * );                /* array of local node numbers               */

extern int count_vardofs        /* rf_util.c                                 */
(const int ,              /* name of variable to be counted            */
       const int  );           /* upper limit on range of nodes to be       */
                                /* searched in Dolphin array                 */

extern void fprint_strn         /* rf_util.c                                 */
(FILE *,
       const char *,            /* string                                    */
       const int);             /* num -> repetition number                  */

extern void fprint_line         /* rf_util.c                                 */
(FILE *,
       const char *,            /* string                                    */
       const int);             /* num -> repetition number                  */

extern int rd_vectors_from_exoII /* rf_util.c                                */
(double [],		/* u - solution vector                       */
	const char *,		/* file_nm - name of EXODUS II file	     */
	const int ,		/* action_flag -                             *
                                       0 -- read initial guess for problem   *
                                       1 -- read extern aux fixed variables  */
	const int ,		/* variable_no                               *
                                 * Used only when action_flag = 1	     *
                                 * Specifies the number of the               *
                                 * external variable to be read              *
                                 * (basically the card number in order)	     */
       const int,               /* Time plane to read from                   *
				 * The default value is INT_MAX, which       *
				 * implies the last time plane in the        *
				 * exodus file.                              */
       double *timeValueRead); /* Value of the time in the time plane read  */

extern int rd_trans_vectors_from_exoII /* rf_util.c                          */
(double [],		/* u - solution vector                       */
	const char *,		/* file_nm - name of EXODUS II file	     */
	const int ,		/* variable_no                               *
                                 * Used only when action_flag = 1	     *
                                 * Specifies the number of the               *
                                 * external variable to be read              *
                                 * (basically the card number in order)	     */
       const int,               /* Time plane to read from                   *
				 * The default value is INT_MAX, which       *
				 * implies the last time plane in the        *
				 * exodus file.                              */
       double *timeValueRead,    /* Value of the time in the time plane read  */
       Comm_Ex *,                /* cx - communications structure */
       Dpi *);                  /* dpi - distributed processing struct */


extern void init_vec_value(double *, const double, const int);
extern void dcopy1(const int, const double *, double *);

extern int find_first_elem_with_var
 ( Exo_DB *, 
	 int  ) ;

#ifdef LIBRARY_MODE
extern int load_import_fields
(dbl *base_p_por,			/* Base porosity for updates */
       const Exo_DB *exo,  		/* Ptr to ExodusII database */
       int callnum);                   /* Indicates first solver call */

extern void interp_ev_to_nodes
(const Exo_DB *exo,		/* Ptr to Exodus database */
       dbl *ev_tmp,			/* Destination array */
       int iev);			/* Import elem var ID */
#endif

extern int advance_porosity_ev
(const int time_step,             /* Time step number in current pass */
       const int nn,                    /* Number of nodes on this proc */
       dbl *x,                          /* Solution at current time step */
       dbl *base_p_por,                 /* Porosity at start of Goma call */
       dbl *base_p_liq);               /* Porous liquid pressure
                                           at start of Goma call */


extern void read_porosity_data
(
       const Exo_DB *exo		/* Ptr to Exodus database */
       );

/*
 * Extern statements for variables defined in rf_util.c
 */

extern Strcpy_rtn strcpy_rtn;

#endif

