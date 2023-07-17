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

#ifndef GOMA_MM_FILL_TERMS_H
#define GOMA_MM_FILL_TERMS_H

#include <stdbool.h>

#include "el_elm.h"
#include "exo_struct.h"
#include "mm_as_structs.h"
#include "mm_fill_common.h"
#include "mm_fill_energy.h"
#include "rf_fem_const.h"
#include "std.h"

struct Boundary_Condition;
struct Data_Table;
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MM_FILL_TERMS_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_MM_FILL_TERMS_C
#define EXTERN extern
#endif

EXTERN int assemble_mesh /* mm_fill_terms.c                           */
    (double,             /* time                                      */
     double,             /* tt                                        */
     double,             /* dt                                        */
     int,                /* ielem - current element number            */
     int,                /* ip - current integration point            */
     int);               /* ip_total - total gauss integration points */

EXTERN int assemble_momentum_path_dependence(
    double,           /* time */
    double,           /* tt, parameter to vary time integration from
                       * explicit (tt = 1) to implicit (tt = 0)    */
    double,           /* dt, current time step size                */
    const PG_DATA *); /* PG data needed for continuity stabilization */

#if 0 
EXTERN int assemble_continuity_path_dependence
(double ,                 /* time_value                                */
       double ,                 /* tt, parameter to vary time integration from *
                                 * explicit (tt = 1) to implicit (tt = 0)    */
       double ,                 /* dt, current time step size                */
       double ,                 /* h_elem_avg,  average global element size  *
                                 * for PSPG, taken constant wrt to Jacobian  */
       double [],               /* hsquared[DIM], element size info for PSPG */
       double [][DIM],		/* hh[DIM][DIM], not currently used, but     *
                                 * left in just in case they're needed later */
       double [][MDE],		/* dh_dxnode[DIM][MDE]                       */
       double ,		        /* U_norm, global velocity norm for PSPG calcs */
       double ,		        /* mu_avg, element viscosity for PSPG calcs  */
       double ,		        /* rho_avg,element density for PSPG calcs    */
       double [],		/* v_avg[DIM], element velocity for PSPG calcs */
       double [][MDE]);        /* dv_dnode[DIM][MDE],deriv.velocity wrt nodal variables */

#endif

EXTERN int assemble_volume /* mm_fill_terms.c                           */
    (bool);

EXTERN int assemble_curvature(void); /* mm_fill_terms.c                         */

EXTERN int assemble_normals /* mm_fill_terms.c                         */
    (void);

EXTERN int assemble_ls_momentum_source(void);

EXTERN int apply_ls_momentum_source(void);

EXTERN int assemble_div_normals /* mm_fill_terms.c                         */
    (void);

EXTERN int assemble_LSvelocity /* mm_fill_terms.c                           */
    (bool, int);

EXTERN int assemble_acoustic /* mm_fill_terms.c                           */
    (double,                 /* time - present time value         */
     double,                 /* tt - parameter to vary time integration
                              * from explicit (tt = 1) to
                              * implicit (tt = 0)                   */
     double,                 /* dt - current time step size               */
     const PG_DATA *,        /* dvc_dnode                                 */
     const int,              /*  acoustic eqn id and var id		     */
     const int);

EXTERN int assemble_poynting /* mm_fill_terms.c                           */
    (double,                 /* time - present time value         */
     double,                 /* tt - parameter to vary time integration
                              * from explicit (tt = 1) to
                              * implicit (tt = 0)                   */
     double,                 /* dt - current time step size               */
     const PG_DATA *,        /* dvc_dnode                                 */
     const int,              /*  Light intensity eqn id and var id		     */
     const int);

EXTERN void restime_nobc_surf /* mm_fill_terms.c                           */
    (double func[MAX_PDIM], double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE]);

EXTERN int assemble_acoustic_reynolds_stress /* mm_fill_terms.c */
    (double,                                 /* time */
     double,                                 /* tt */
     double,                                 /* dt */
     const PG_DATA *);                       /* dvc_dnode */

EXTERN int assemble_pore_sink_mass /* mm_fill_terms.c                           */
    (double,                       /* time - present time value         */
     double,                       /* tt - parameter to vary time integration
                                    * from explicit (tt = 1) to
                                    * implicit (tt = 0)                   */
     double);                      /* dt - current time step size               */

EXTERN int assemble_acoustic_energy /* mm_fill_terms.c                 */
    (double,                        /* time - present time value         */
     double,                        /* tt - parameter to vary time integration
                                     * from explicit (tt = 1) to
                                     * implicit (tt = 0)                 */
     double,                        /* dt - current time step size               */
     const PG_DATA *);              /* dvc_dnode                                 */

EXTERN double acoustic_impedance             /* mm_fill_terms.c             */
    (CONDUCTIVITY_DEPENDENCE_STRUCT *, dbl); /* time */

EXTERN double wave_number                    /* mm_fill_terms.c             */
    (CONDUCTIVITY_DEPENDENCE_STRUCT *, dbl); /* time */

EXTERN double acoustic_absorption            /* mm_fill_terms.c             */
    (CONDUCTIVITY_DEPENDENCE_STRUCT *, dbl); /* time */

EXTERN double refractive_index               /* mm_fill_terms.c             */
    (CONDUCTIVITY_DEPENDENCE_STRUCT *, dbl); /* time */

EXTERN double light_absorption               /* mm_fill_terms.c             */
    (CONDUCTIVITY_DEPENDENCE_STRUCT *, dbl); /* time */

EXTERN double extinction_index               /* mm_fill_terms.c             */
    (CONDUCTIVITY_DEPENDENCE_STRUCT *, dbl); /* time */

EXTERN int ls_modulate_momentumsource(
    double[DIM], double[DIM], double, double, double, MOMENTUM_SOURCE_DEPENDENCE_STRUCT *);

EXTERN void apply_table_mp(double *func, struct Data_Table *table);

EXTERN dbl solidification_permeability(dbl, dbl[MAX_CONC][MDE]); /* continuous surface tension */

EXTERN int continuous_surface_tension(double,
                                      double[DIM][DIM], /* continuous surface tension            */
                                      double[DIM][DIM][MDE], /* derivative w.r.t. FILL       */
                                      double[DIM][DIM][DIM][MDE]); /* d with respect to mesh */

EXTERN double quad_isomap_invert(const double,   /*  coordinate1  */
                                 const double,   /*  coordinate2  */
                                 const double,   /*  coordinate3  */
                                 const double[], /* grid points of coordinate1 */
                                 const double[], /* grid points of coordinate2 */
                                 const double[], /* grid points of coordinate3 */
                                 const double[], /* function values at grid points */
                                 const int,      /* # of grid points in direction 1 */
                                 const int,      /* # of grid points in direction 2 */
                                 const int,      /* # of grid points in direction 3 */
                                 const int,      /* element order(2=biquadratic, 1=bilinear) */
                                 const int,      /* element dimension */
                                 double[]);      /* gradient array */

extern void load_matrl_statevector(MATRL_PROP_STRUCT *);

EXTERN int apply_distributed_sources(int,
                                     double,
                                     double[],
                                     Exo_DB *,
                                     double,
                                     double,
                                     double,
                                     const PG_DATA *,
                                     int,
                                     double *,
                                     double **,
                                     double **);

EXTERN int assemble_curvature_with_normals_source(void);

EXTERN int assemble_curvature_source(void);

EXTERN void grad_vector_fv_fill(double ***, double (*)[DIM][DIM][DIM], int, double (*)[DIM]);

EXTERN void grad_scalar_fv_fill(double **, double (*)[DIM], int, double *);

EXTERN void
scalar_fv_fill(double **, double **, double **, double *, int, double *, double *, double *);

EXTERN double scalar_fv_fill_adjmatrl(double **, int, int, int);

EXTERN int assemble_q_source(double);

EXTERN int assemble_qlaser_source(const double[], double);

EXTERN int assemble_qvapor_source(const double[]);

EXTERN int assemble_qrad_source(double, double, double, double);

EXTERN int assemble_t_source(double, double);

EXTERN int assemble_cont_t_source(double *);

EXTERN int assemble_ls_yflux_source(
    int, double, double, double, double, double, int, struct Boundary_Condition *);

EXTERN int assemble_ars_source(double, double);

EXTERN int assemble_cont_vel_source(double *, Exo_DB *);

EXTERN int assemble_extv_kinematic(dbl, dbl, dbl, int, struct Boundary_Condition *);

EXTERN int assemble_p_source(double, const int);

EXTERN int assemble_precoil_source(const double[]);

EXTERN int assemble_uvw_source(int, double);

EXTERN int assemble_extension_velocity_path_dependence(void);

EXTERN int assemble_LM_source(double *, int, double *, double **, double **, double[], Exo_DB *);

EXTERN int assemble_interface_extension_velocity_sic(int);

EXTERN int assemble_eik_kinematic(dbl, dbl, dbl, int, struct Boundary_Condition *);

EXTERN int assemble_fill_path_dependence(void);

EXTERN int assemble_energy_path_dependence(double, /* time - present time value                 */
                                           double, /* tt - parameter to vary time integration
                                                    * from explicit (tt = 1) to
                                                    * implicit (tt = 0)                         */
                                           double,
                                           const PG_DATA *); /* dvc_dnode */

EXTERN int assemble_continuity_path_dependence(
    dbl,
    dbl, /* parameter to vary time integration from explicit (tt = 1) to implicit (tt = 0)    */
    dbl, /* current time step size                    */
    const PG_DATA *); /* deriv. velocity wrt nodal variables   */

EXTERN void acoustic_flux(double[DIM],                       /* q[DIM] */
                          ACOUSTIC_FLUX_DEPENDENCE_STRUCT *, /* dq     */
                          double,                            /* time   */
                          const int,                         /* acoustic eqn id and var id	*/
                          const int);

EXTERN int assemble_pf_capillary(double *);

EXTERN int assemble_max_strain(void);

EXTERN int assemble_cur_strain(void);

int assemble_ls_stress_jump(double viscosity_scale, double stress_scale, int heaviside_type);

#endif /* GOMA_MM_FILL_TERMS_H */
