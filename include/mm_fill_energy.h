/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2023 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/

#ifndef GOMA_MM_FILL_ENERGY_H
#define GOMA_MM_FILL_ENERGY_H

#include "el_elm.h"
#include "exo_struct.h"
#include "mm_as_structs.h"
#include "mm_fill_common.h"
#include "rf_fem_const.h"
#include "std.h"
#include <stdbool.h>

/* struct for d_k */
struct conductivity_dependence {
  double X[DIM][MDE];      /* mesh dependence. */
  double T[MDE];           /* temperature dependence. */
  double C[MAX_CONC][MDE]; /* conc dependence. */
  double F[MDE];           /* FILL dependence. */
  double moment[MAX_MOMENTS][MDE];
};
typedef struct conductivity_dependence CONDUCTIVITY_DEPENDENCE_STRUCT;

/* struct for d_Cp */
struct heat_capacity_dependence {
  double v[DIM][MDE];      /* velocity dependence. */
  double X[DIM][MDE];      /* mesh dependence. */
  double T[MDE];           /* temperature dependence. */
  double C[MAX_CONC][MDE]; /* conc dependence. */
  double V[MDE];           /* voltage dependence. */
  double F[MDE];           /* FILL dependence. */
};
typedef struct heat_capacity_dependence HEAT_CAPACITY_DEPENDENCE_STRUCT;

/* struct for d_h */
struct heat_source_dependence {
  double v[DIM][MDE];                 /* velocity dependence. */
  double X[DIM][MDE];                 /* mesh dependence. */
  double T[MDE];                      /* temperature dependence. */
  double C[MAX_CONC][MDE];            /* conc dependence. */
  double V[MDE];                      /* voltage dependence. */
  double S[MAX_MODES][DIM][DIM][MDE]; /* stress mode dependence. */
  double F[MDE];                      /* level set field dependence */
  double P[MDE];                      /* acoustic pressure dependence  */
  double APR[MDE];                    /* acoustic pressure dependence  */
  double API[MDE];                    /* acoustic pressure dependence  */
  double INT[MDE];                    /* acoustic pressure dependence  */
  double EM_ER[DIM][MDE];             /* time-harmonic electromagnetic dependence */
  double EM_EI[DIM][MDE];             /* time-harmonic electromagnetic dependence */
  double rst[MDE];                    /* residence time field dependence  */
};
typedef struct heat_source_dependence HEAT_SOURCE_DEPENDENCE_STRUCT;

int assemble_energy /* mm_fill_terms.c                           */
    (double,        /* time - present time value                 */
     double,        /* tt - parameter to vary time integration
                     * from explicit (tt = 1) to
                     * implicit (tt = 0)                   */
     double,        /* dt - current time step size        */
     const PG_DATA *);
double conductivity /* mm_fill_terms.c             */
    (CONDUCTIVITY_DEPENDENCE_STRUCT *, dbl);
double heat_capacity /* mm_fill_terms.c                  */
    (HEAT_CAPACITY_DEPENDENCE_STRUCT *, dbl);
double ls_modulate_thermalconductivity(
    double, double, double, double, double, CONDUCTIVITY_DEPENDENCE_STRUCT *);
double
ls_modulate_heatcapacity(double, double, double, double, double, HEAT_CAPACITY_DEPENDENCE_STRUCT *);
void heat_flux(double[DIM],                   /* q[DIM] */
               HEAT_FLUX_DEPENDENCE_STRUCT *, /* dq     */
               double);
double heat_source /* mm_fill_terms.c                  */
    (HEAT_SOURCE_DEPENDENCE_STRUCT *,
     double, /* time - present time value                 */
     double, /* tt - parameter to vary time integration
              * from explicit (tt = 1) to
              * implicit (tt = 0)                         */
     double);
int ls_modulate_heatsource(
    double *, double, double, double, double, HEAT_SOURCE_DEPENDENCE_STRUCT *);
int assemble_ls_latent_heat_source(
    double, double, double, double, double, int, struct Boundary_Condition *);
double visc_diss_acoustic_source(HEAT_SOURCE_DEPENDENCE_STRUCT *,
                                 dbl *, /* param - General multipliers   */
                                 int);
double em_diss_heat_source(HEAT_SOURCE_DEPENDENCE_STRUCT *,
                           dbl *, /* param - General multipliers   */
                           int);
double em_diss_e_curlcurl_source(HEAT_SOURCE_DEPENDENCE_STRUCT *,
                                 dbl *, /* param - General multipliers   */
                                 int);
#endif // GOMA_MM_FILL_ENERGY_H
