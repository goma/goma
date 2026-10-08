/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2022 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Authors: Robert Secor, Andrew Cochrane, and Weston Ortiz                *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/

#ifndef GOMA_MM_FILL_EM_H
#define GOMA_MM_FILL_EM_H

#include "el_elm.h"
#include "mm_as_structs.h"
#include "mm_fill_terms.h"
#include "rf_fem_const.h"
#ifdef EXTERN
#undef EXTERN
#endif

#include <complex.h>
#undef I

#ifdef GOMA_MM_FILL_EM_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_MM_FILL_EM_C
#define EXTERN extern
#endif

enum em_incident_wave_type {
  // NONE/CONSTANT in input deck
  EM_INC_NONE = 1,
  // Start at 2 to avoid CONSTANT=1
  EM_INC_PLANE_Z_WAVE = 2,
};

enum stabilization_type {
  EM_STAB_NONE,
  EM_STAB_PHI_DIV,
  EM_STAB_DPHI_DIV,
  EM_STAB_DIVPHI_DIV,
  EM_STAB_PHI_DIVSQUARED,
  EM_STAB_DPHI_DIVSQUARED
};

struct emwave_stabilization {
  enum stabilization_type type;
  int em_eqn;
  int em_var;
  int stabilization_field_var;
  double residual_term[MDE];
  double jacobian_term[MDE][DIM][MDE];
};

EXTERN int assemble_emwave /* mm_fill_em.c                           */
    (double,               /* time - present time value         */
     double,               /* tt - parameter to vary time integration
                            * from explicit (tt = 1) to
                            * implicit (tt = 0)                   */
     double,               /* dt - current time step size               */
     const PG_DATA *,      /* dvc_dnode                                 */
     const int,            /*  Light intensity eqn id and var id     */
     const int,            /*  Light intensity eqn id and var id     */
     const int);

EXTERN void calc_emwave_stabilization_term(struct emwave_stabilization *, double);

int assemble_ewave_curlcurl(double time,       // present time
                            double tt,         // time integration method parameter
                            double dt,         // current time step size
                            const int em_eqn,  // eqn id
                            const int em_var); //  variable id - should match me_eqn

int assemble_ewave_laplacian(double time,       // present time
                             double tt,         // time integration method parameter
                             double dt,         // current time step size
                             const int em_eqn,  // eqn id
                             const int em_var); //  variable id - should match me_eqn

int assemble_em_continuity(void);
int assemble_ewave_nedelec(dbl time);
int em_mms_force(dbl x, dbl y, dbl z, complex double force[DIM]);
int em_mms_exact(dbl x, dbl y, dbl z, complex double exact[DIM]);
int incident_wave(
    dbl x, dbl y, dbl z, dbl omega, complex double wave[DIM], complex double curl_wave[DIM]);
bool relative_permittivity_model(complex double *permittivity_out,
                                 complex double *permittivity_matrix);
int apply_ewave_nedelec_farfield(double func[DIM],
                                 double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                                 double xi[DIM], /* Local stu coordinates */
                                 double time,    // present time
                                 const int bc_name,
                                 double *bc_data);
#endif /* GOMA_MM_FILL_EM_H */
