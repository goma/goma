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

#ifndef GOMA_MM_FLUX_H
#define GOMA_MM_FLUX_H

#include "dpi.h"
#include "el_elm.h"
#include "exo_struct.h"
#include "mm_fill_terms.h"
#include "std.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MM_FLUX_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_MM_FLUX_C
#define EXTERN extern
#endif

EXTERN double evaluate_flux(const Exo_DB *exo,      /* ptr to basic exodus ii mesh information */
                            const Dpi *dpi,         /* distributed processing info */
                            const int side_set_id,  /* on which SSID to evaluate flux */
                            const int quantity,     /* to pick HEAT_FLUX, FORCE_NORMAL, etc. */
                            const char *qtity_str,  /* quantity string */
                            const int blk_id,       /* material identification */
                            const int species_id,   /* species identification */
                            const char *filenm,     /* File name pointer */
                            const int profile_flag, /*  flag for printing flux profiles  */
                            dbl *x,                 /* solution vector */
                            dbl *xdot,              /* dx/dt vector */
                            double J_AC[],         /* vector for augmenting condition sensitivities.
                                    May be NULL */
                            const double delta_t,  /* time-step size */
                            const dbl time_value,  /* current time */
                            const int print_flag); /*  flag for printing results,1=print*/

EXTERN double evaluate_global_flux(const Exo_DB *exo,
                                   const Dpi *dpi,
                                   const int quantity,
                                   const int blk_id,
                                   const int species_id,
                                   const double *params,
                                   double *J_AC,
                                   double x[],
                                   const dbl time_value,
                                   const int print_flag);

EXTERN double evaluate_flux_sens(const Exo_DB *exo, /* ptr to basic exodus ii mesh information */
                                 const Dpi *dpi,    /* distributed processing info */
                                 const int side_set_id,  /* on which SSID to evaluate flux */
                                 const int quantity,     /* to pick HEAT_FLUX, FORCE_NORMAL, etc. */
                                 const char *qtity_str,  /* quantity string */
                                 const int mat_id,       /* material identification */
                                 const int species_id,   /* species identification */
                                 const int sens_type,    /*  sensitivity type */
                                 const int sens_id,      /* sensitivity id */
                                 const int sens_flt,     /*  sensitivity float number */
                                 const int sens_flt2,    /*  sensitivity float number for UM */
                                 const int vector_id,    /* sensitivity id */
                                 const char *filenm,     /* File name pointer */
                                 const int profile_flag, /*  flux sens print flag  */
                                 double x[],             /* solution vector */
                                 double xdot[],          /* solution vector */
                                 double **x_sens_p,      /* sensitivity vector */
                                 const double delta_t,   /* time-step size */
                                 const double time_value, /* current time */
                                 const int print_flag);   /*  printing control flag */

EXTERN double
evaluate_volume_integral(const Exo_DB *exo,        /* ptr to basic exodus ii mesh information */
                         const Dpi *dpi,           /* distributed processing info */
                         const int quantity,       /* to pick VOLUME, DISSIPATION, etc. */
                         const char *quantity_str, /* volume integral name */
                         const int blk_id,         /* material identification */
                         const int species_id,     /* species identification */
                         const char *filenm,       /* File name pointer */
                         const double *params,
                         const int num_params,
                         double *J_AC,         /* Pointer to AC sensitivity vector, may be NULL */
                         double x[],           /* solution vector */
                         double xdot[],        /* dx/dt vector */
                         const double delta_t, /* time-step size */
                         const double time_value, /* current time */
                         const int print_flag);   /*  flag for printing results,1=print*/

EXTERN int compute_volume_integrand(const int,
                                    const int,
                                    const int,
                                    const int,
                                    const double *,
                                    const int,
                                    double *,
                                    double *,
                                    const int,
                                    const double,
                                    const double,
                                    double[],
                                    const Exo_DB *);

EXTERN void
compute_surface_integrand(const int, int, const int, const double *, double *, double[]);

EXTERN int
adaptive_weight(double *, const int, const int, const double *, const double, const int, const int);

EXTERN int solve_quadratic(const double, const double, const double, double *);

#ifndef NO_CHEBYSHEV_PLEASE
EXTERN int chebyshev_coeff_2DQ(
    const int, const double *, double[6][2], int, double *, double *, int *, const int *);

EXTERN void heaviside_chev_moments_2DQ(
    const int, const int, double *, const double *, const int, const double *);
EXTERN void surfdet_chev_coeff_2DQ(
    const int, const double *, const double *, double *, const int, const int, const double *);

EXTERN void delta_chev_moments_2DQ(
    const int, double *, const double *, const double *, const int, const double *);

#endif

EXTERN int interface_crossing_1DQ(const double *, double[2]);

EXTERN int interface_crossing_2DQ(const double *, double[6][2], int *, int *, double[12][MAX_PDIM]);

EXTERN void
interface_inclination_2DQ(const double *, const int, double *, const int *, double[6][2]);

EXTERN int interface_crossing_3DL(const double *, double[12][2], int *, double[12][MAX_PDIM]);

/*  Chebyshev Polynomial data structure  */

#ifdef GOMA_MM_FLUX_C
#ifndef GOMA_MM_POST_PROC_UTIL_C
#define MAX_CHEV 5

struct Chebyshev_Polynomial {
  int order;
  double root[MAX_CHEV];
  double cosval[MAX_CHEV * MAX_CHEV];
};

typedef struct Chebyshev_Polynomial CHEV_POLY_STRUCT;

extern CHEV_POLY_STRUCT chevpoly[];
struct Chebyshev_Polynomial chevpoly[3] = {
    {3,
     {0.866025403784438646763723170753, 0, -0.866025403784438646763723170753},
     {1, 1, 1, 0.866025403784438646763723170753, 0, -0.866025403784438646763723170753, 0.5, -1.,
      0.5}},
    {4,
     {0.9238795325112867561281831893968, 0.38268343236508977172845998403,
      -0.3826834323650897717284599840304, -0.923879532511286756128183189397},
     {1, 1, 1, 1, 0.9238795325112867561281831893968, 0.38268343236508977172845998403,
      -0.3826834323650897717284599840304, -0.923879532511286756128183189397,
      0.707106781186547524400844362105, -0.707106781186547524400844362105,
      -0.707106781186547524400844362105, 0.707106781186547524400844362105,
      0.38268343236508977172845998403, -0.923879532511286756128183189397,
      0.923879532511286756128183189397, -0.38268343236508977172845998403}},
    {5,
     {0.951056516295153572116439333379, 0.587785252292473129168705954639, 0,
      -0.587785252292473129168705954639, -0.951056516295153572116439333379},
     {1,
      1,
      1,
      1,
      1,
      0.951056516295153572116439333379,
      0.587785252292473129168705954639,
      0,
      -0.587785252292473129168705954639,
      -0.951056516295153572116439333379,
      0.809016994374947424102293417183,
      -0.309016994374947424102293417183,
      -1.,
      -0.309016994374947424102293417183,
      0.809016994374947424102293417183,
      0.587785252292473129168705954639,
      -0.951056516295153572116439333379,
      0,
      0.951056516295153572116439333379,
      -0.587785252292473129168705954639,
      0.309016994374947424102293417183,
      -0.809016994374947424102293417183,
      1.,
      -0.809016994374947424102293417183,
      0.309016994374947424102293417183}}};
#endif
#endif /* GOMA_MM_FLUX_C */

#endif /* GOMA_MM_FLUX_H */
