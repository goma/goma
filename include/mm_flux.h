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
 

#ifndef _MM_FLUX_H
#define _MM_FLUX_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _MM_FLUX_C
#define EXTERN /* do nothing */
#endif

#ifndef _MM_FLUX_C
#define EXTERN extern
#endif


EXTERN double evaluate_flux	/* mm_flux.c                                 */
PROTO((const Exo_DB *,		/* exo - ptr to basic exodus ii mesh info    */
       const Dpi *,             /* distributed processing info */
       const int ,		/* side_set_id - which SSID to evaluate flux */
       const int ,		/* quantity - HEAT_FLUX, FORCE_NORMAL, etc.  */
       const char * ,           /* quantity string  */
       const int ,		/* mat_id - material identification          */
       const int ,		/* species_id - species identification       */
       const char *,		/* filenm - File name pointer                */
       const int ,              /* flux profile print control flag           */
       const double [],		/* x - solution vector                       */
       const double [],		/* xdot - dx/dt vector                       */
             double [],         /* J_AC - augmenting condition sensititives */
       const double ,		/* delta_t - time-step size                  */
       const double ,              /* time_value - current time                 */
       const int ));            /* print flag                                */


EXTERN double evaluate_global_flux
PROTO((const Exo_DB *,
       const Dpi *,
       const int,
       const int,
       const int,
       const double *,
       double *,
       const double [],
       const double,
       const int ));

EXTERN double evaluate_flux_sens        /* mm_flux.c                                 */
PROTO((const Exo_DB *,          /* exo - ptr to basic exodus ii mesh info    */
       const Dpi *,             /* distributed processing info */
       const int ,              /* side_set_id - which SSID to evaluate flux */
       const int ,              /* quantity - HEAT_FLUX, FORCE_NORMAL, etc.  */
       const char * ,           /* quantity string  */
       const int ,              /* mat_id - material identification          */
       const int ,              /* species_id - species identification       */
       const int ,              /*  sensitivity type (1 or 2) = (BC or MT)   */
       const int ,              /* id of sensitivity variable         */
       const int ,              /* id of sensitivity variable float        */
       const int ,              /* id of sensitivity variable float        */
       const int ,              /* vector id of sensitivity variable         */
       const char *,            /* filenm - File name pointer                */
       const int ,              /* flux sens profile print control flag    */
       const double [],         /* x - solution vector                       */
       const double [],         /* xdot -dx/dt vector                    */
       double **,               /* x_sens_p - sensitivity vector     */
       const double ,           /* delta_t - time-step size                  */
       const dbl ,              /* time_value - current time                 */
       const int ));            /* print flag                                */

EXTERN double evaluate_volume_integral 
PROTO((const Exo_DB *,
       const Dpi *,
       const int,
       const char *,
       const int,
       const int,
       const char *,
       const double *,
       const int,
       double [],
       const double [],
       const double [],
       const double,
       const double,  
       const int ));

EXTERN int compute_volume_integrand
PROTO((const int,
       const int,
       const int,
       const double *,
       const int,
       double *,
       double *,
       const int,
       const double,
       const double,
       double [],
       const Exo_DB *));

EXTERN void compute_surface_integrand
PROTO((const int,
       int,
       const int,
       const double *,
       double *,
       double [] ));

EXTERN int adaptive_weight
PROTO((double *, 
       const int,
       const int,
       const double *,
       const double,
       const int,
       const int ));

EXTERN int solve_quadratic
PROTO((const double,
       const double, 
       const double, 
       double * ));

#ifndef NO_CHEBYSHEV_PLEASE
EXTERN int chebyshev_coeff_2DQ
PROTO((const int ,
       const double *,
       double [][2],
       int ,
       double *,
       double *,
       int *,
       const int * ));

EXTERN void heaviside_chev_moments_2DQ
PROTO (( const int,
		 const int, 
		 double *, 
		 const double *,
		 const int,
		 const double * ));
EXTERN void surfdet_chev_coeff_2DQ
PROTO (( const int ,
		 const double *,
		 const double *,
		 double *,
		 const int,
		 const int,
		 const double *  ));

EXTERN void delta_chev_moments_2DQ
PROTO (( const int,
		 double *, 
		 const double *,
		 const double *,
		 const int,
		 const double * ));

#endif

EXTERN int interface_crossing_1DQ
PROTO((const double *,
       double [2]));

EXTERN int interface_crossing_2DQ
PROTO((const double *,
       double [][2],
       int *,
       int *,
       double [][MAX_PDIM] ));

EXTERN void interface_inclination_2DQ
PROTO(( const double *,
	const int ,
	double *,
	const int *,
	double [][2] ));

EXTERN int interface_crossing_3DL
PROTO (( const double *,
	 double [][2],
	 int *,
	 double [][MAX_PDIM] ));



/*  Chebyshev Polynomial data structure  */

#ifdef _MM_FLUX_C
#ifndef _MM_POST_PROC_UTIL_C
#define MAX_CHEV  5

struct Chebyshev_Polynomial
{ int order;
  double root[MAX_CHEV];
  double cosval[MAX_CHEV*MAX_CHEV];
};

typedef struct Chebyshev_Polynomial CHEV_POLY_STRUCT;

extern CHEV_POLY_STRUCT chevpoly[];
struct Chebyshev_Polynomial chevpoly[3] = {
{3, 
{0.866025403784438646763723170753, 0, -0.866025403784438646763723170753},
{1,1,1,0.866025403784438646763723170753, 0, -0.866025403784438646763723170753
,0.5, -1., 0.5}
 },
{4, 
{0.9238795325112867561281831893968, 0.38268343236508977172845998403
, -0.3826834323650897717284599840304, -0.923879532511286756128183189397},
{1,1,1,1, 0.9238795325112867561281831893968, 0.38268343236508977172845998403
, -0.3826834323650897717284599840304, -0.923879532511286756128183189397
, 0.707106781186547524400844362105, -0.707106781186547524400844362105
, -0.707106781186547524400844362105, 0.707106781186547524400844362105
, 0.38268343236508977172845998403, -0.923879532511286756128183189397
, 0.923879532511286756128183189397, -0.38268343236508977172845998403}
},
{5,
{0.951056516295153572116439333379, 0.587785252292473129168705954639, 0
, -0.587785252292473129168705954639, -0.951056516295153572116439333379},
{1,1,1,1,1, 0.951056516295153572116439333379, 0.587785252292473129168705954639
, 0, -0.587785252292473129168705954639, -0.951056516295153572116439333379
, 0.809016994374947424102293417183, -0.309016994374947424102293417183, -1.
, -0.309016994374947424102293417183, 0.809016994374947424102293417183
, 0.587785252292473129168705954639, -0.951056516295153572116439333379, 0
, 0.951056516295153572116439333379, -0.587785252292473129168705954639
, 0.309016994374947424102293417183, -0.809016994374947424102293417183, 1.
, -0.809016994374947424102293417183, 0.309016994374947424102293417183}
}
};
#endif
#endif /* _MM_FLUX_C */


#endif /* _MM_FLUX_H */
