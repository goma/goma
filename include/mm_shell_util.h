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

/*
 *$Id: mm_shell_util.h,v 5.2 2010-05-21 21:40:35 sarober Exp $
 */

#ifndef GOMA_MM_SHELL_UTIL_H
#define GOMA_MM_SHELL_UTIL_H

#include "el_elm.h"
#include "exo_struct.h"
#include "rf_fem_const.h"

struct Generalized_Newtonian;

extern int **elem_friends;
extern int *num_elem_friends;
extern int num_shell_blocks;

extern void init_shell_element_blocks(const Exo_DB *exo); /* ExodusII database struct pointer */

extern int is_shell_element_type(const int elem_type); /* Element type index */

extern int is_shell_element(const int elem,     /* Element index */
                            const Exo_DB *exo); /* ExodusII database struct pointer */

extern int is_shell_block(const int block_id, /* Element block ID */
                          const Exo_DB *exo); /* ExodusII database struct pointer */

extern int solve_2x2 /* Solves A*x = b */
    (double a[DIM][DIM], double b[DIM], double x[DIM]);

extern int solve_3x3 /* Solves A*x = b */
    (double a[DIM][DIM], double b[DIM], double x[DIM]);

extern int find_stu_from_xyz /* Inverts isoparametrix map to get xi */
    (const int,              /* Element of interest */
     const double[DIM],      /* xyz coords in */
     double[DIM],            /* stu coords out */
     const Exo_DB *exo);     /* Ptr to Exodus structure */

extern int bulk_side_id_and_stu /* Determines bulk elem side shared by shell */
    (const int,                 /* Bulk element number */
     const int,                 /* Shell element number */
     const double[DIM],         /* Shell stu coords in */
     double[DIM],               /* Bulk stu coords out */
     const Exo_DB *exo);        /* Ptr to Exodus structure */

extern int load_neighbor_var_data          /* Creates neighbor assembly structures */
    (int el1,                              /* Local element number */
     int el2,                              /* Neighbor element number */
     int *ndofs,                           /* Neighbor variable DOF's */
     int *dof_map,                         /* Shell-bulk local DOF conversion */
     int ndofptr[MAX_VARIABLE_TYPES][MDE], /* neighbor dof pointer arrays	*/
     const int id_side,                    /* id_side for bulk BC's */
     double xi[DIM],                       /* Local element stu coords */
     const Exo_DB *exo);                   /* Ptr to Exodus structure */

extern int find_stu_on_shell /* Converts bulk stu to shell stu */
    (const int bulk_elem,    /* Bulk element number */
     const int id_side,      /* Bulk element side ID (Exo/Patran) */
     const int shell_elem,   /* Shell element number */
     const int bulk_dim,     /* 2D or 3D bulk element */
     const double s,         /* Bulk xi[0] */
     const double t,         /* Bulk xi[1] */
     const double u,         /* Bulk xi[2] */
     double xi2[DIM],        /* Shell stu coords */
     const Exo_DB *exo);

extern int find_stu_on_bulk /* Converts shell stu to bulk stu */
    (const int id_side,     /* Bulk element side ID (Exo/Patran) */
     const int bulk_dim,    /* 2D or 3D bulk element */
     const double s,        /* Shell xi[0] */
     const double t,        /* Shell xi[1] */
     double xi2[DIM]);      /* Bulk stu coords */

extern int shell_normal_div_s            /* Loads shell normal surface divergence */
    (double *div_s_nv,                   /* Surface divergence of normal (scalar) */
     double d_div_s_nv_dnv[DIM][MDE],    /* Self-sensitivities */
     double d_div_s_nv_dmesh[DIM][MDE]); /* Mesh sensitivities */

extern void shell_tangents(double t0[DIM],
                           double t1[DIM],
                           double dt0_dx[DIM][DIM][MDE],
                           double dt1_dx[DIM][DIM][MDE],
                           double dt0_dnormal[DIM][DIM][MDE],
                           double dt1_dnormal[DIM][DIM][MDE]);
extern void shell_tangents_isoparametric(double t0[DIM],
                                         double t1[DIM],
                                         double dt0_dx[DIM][DIM][MDE],
                                         double dt1_dx[DIM][DIM][MDE]);

extern void shell_tangents_seeded(double t0[DIM],
                                  double t1[DIM],
                                  double dt0_dnormal[DIM][DIM][MDE],
                                  double dt1_dnormal[DIM][DIM][MDE]);

extern void shell_stress_tensor(double TT[DIM][DIM],
                                double dTT_dx[DIM][DIM][DIM][MDE],
                                double dTT_dnormal[DIM][DIM][DIM][MDE]);

extern void shell_moment_tensor(double M[DIM][DIM],
                                double dM_dx[DIM][DIM][DIM][MDE],
                                double dM_dnormal[DIM][DIM][DIM][MDE],
                                double dM_curv1[DIM][DIM][MDE],
                                double dM_curv2[DIM][DIM][MDE]);

extern void lubrication_shell_initialize(int *n_dof,        // Degrees of freedom
                                         int *dof_map,      // Map of DOFs
                                         int id_side,       // Side ID
                                         double xi[DIM],    // Local STU coordinates
                                         const Exo_DB *exo, // Exodus database
                                         int use_def        // Use deformed normal anyway
);

extern void Inn(double v[DIM], // Input vector
                double w[DIM]  // Rotated output vector
);

extern void ShellRotate(double v[DIM],                // Input vector
                        double dv_dmx[DIM][DIM][MDE], // Input vector mesh sensitivity
                        double w[DIM],                // Rotated output vector
                        double dw_dmx[DIM][DIM][MDE], // Rotated output vector mesh sensitivity
                        int ndof                      // Number of DOFs for mesh equations
);

extern void
calculate_lub_q_v(const int EQN, double time, double dt, double xi[DIM], const Exo_DB *exo);

extern void calculate_lub_q_v_old(
    const int EQN, double time_old, double dt_old, double xi[DIM], const Exo_DB *exo);

extern void ShellBF(int ev, // Equation or variable to fetch basis functions
                    int ii, // Integer for which DOF to fetch
                    double *phi_i,
                    double grad_phi_i[DIM],
                    double gradII_phi_i[DIM],
                    double d_gradII_phi_i_dx[DIM][DIM][MDE],
                    int ndof,
                    int dof_map[MDE] // Map of DOFs
);

extern void ShellBF_2d_bar(int ev, // Equation or variable to fetch basis functions
                           int ii, // Integer for which DOF to fetch
                           double gradII_phi_i[DIM],
                           double d_gradII_phi_i_dx[DIM][DIM][MDE]);

extern double shell_saturation_pressure_curve(double P, double *dSdP, double *dSdP_P);

extern void calculate_lub_q_v_nonnewtonian(double time, double dt);

extern void calculate_lub_q_v_nonnewtonian_sens(double mu,
                                                double H,
                                                double veloU[DIM],
                                                double veloL[DIM],
                                                double grad_II_P[DIM],
                                                double phi_j,
                                                double grad_phi_j[DIM],
                                                double shear_top_plus,
                                                double shear_top_minus,
                                                double shear_bot_plus,
                                                double shear_bot_minus,
                                                double eps_top,
                                                double eps_bot);

extern int lub_viscosity_integrate(const double strs,
                                   const double H,
                                   double *flow_mag,
                                   double *dq_gradp,
                                   double *dq_dh,
                                   double *srate,
                                   double *pre_P,
                                   double *mu_star,
                                   double *dq_dT);

extern double
lub_viscos_fcn(const struct Generalized_Newtonian *gn, const double gammadot, double *visd);

extern double lub2D_crsrate(const struct Generalized_Newtonian *gn_local,
                            const double rate,
                            const double strs,
                            const double tolerance);

extern int lub2D_viscint_2D(const struct Generalized_Newtonian *gn_local,
                            const double rate0,
                            const double rate1,
                            const double Kconst,
                            double *xint1,
                            double *xint2,
                            double *xint1dK,
                            double *xint2dK,
                            const double tolerance);

extern double lub2D_flowint_2D(const struct Generalized_Newtonian *gn_loc,
                               const double rate0,
                               const double rate1,
                               const double Kconst,
                               const double tolerance);

extern int lub2D_numjac(const double h,
                        const double pgrad,
                        const double rate0,
                        const double rate1,
                        const double Kconst,
                        const double epstol,
                        const double ev[DIM],
                        const double fv[DIM],
                        const double veloU[DIM],
                        const double veloL[DIM],
                        struct Generalized_Newtonian *gn_loc,
                        double xj[3][3]);

extern int lub2D_flow2D(struct Generalized_Newtonian *gn_loc,
                        const double gradp[DIM],
                        const double ev[DIM],
                        const double dev_dpg[DIM][DIM],
                        double qfl[DIM],
                        double dqfl_dp[DIM][DIM],
                        double dqfl_dh[DIM],
                        const double h,
                        double *wrate,
                        const double veloL[DIM],
                        const double veloU[DIM],
                        int iguess,
                        double *wstrs,
                        double relax_rate);

#endif /* GOMA_MM_FILL_SHELL_H */
