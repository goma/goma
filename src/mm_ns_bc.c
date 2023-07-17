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
 *$Id: mm_ns_bc.c,v 5.26 2010-07-21 16:39:27 hkmoffa Exp $
 */

/* Standard include files */

#include "mm_ns_bc.h"

#include <math.h>
#include <mm_fill_stabilization.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* GOMA include files */
#include "ac_stability.h"
#include "az_aztec.h"
#include "bc_colloc.h"
#include "density.h"
#include "el_elm.h"
#include "el_geom.h"
#include "exo_struct.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_elem_block_structs.h"
#include "mm_fill_energy.h"
#include "mm_fill_jac.h"
#include "mm_fill_ls.h"
#include "mm_fill_momentum.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_rs.h"
#include "mm_fill_solid.h"
#include "mm_fill_species.h"
#include "mm_fill_stress.h"
#include "mm_fill_terms.h"
#include "mm_input.h"
#include "mm_interface.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "mm_qp_storage.h"
#include "mm_qtensor_model.h"
#include "mm_shell_util.h"
#include "mm_species.h"
#include "mm_unknown_map.h"
#include "mm_viscosity.h"
#include "mpi.h"
#include "rd_mesh.h"
#include "rf_allo.h"
#include "rf_bc.h"
#include "rf_bc_const.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_mp.h"
#include "rf_solver.h"
#include "rf_vars_const.h"
#include "sl_auxutil.h"
#include "std.h"
#include "user_bc.h"
#include "user_mp.h"
#include "wr_side_data.h"

#define eps(i, j, k) ((i - j) * (j - k) * (k - i) / 2)

#define GOMA_MM_NS_BC_C

/*
 * Global variables defined here. Declared frequently via rf_bc.h
 */

/*********** R O U T I N E S   I N   T H I S   F I L E ***********************
 *
 *						-All routines in this
 *				NOTE:		 file are called only by
 *						 mm_fill.c: matrix_fill
 *						 except possibly for static
 *						 functions.
 *
 *       NAME			TYPE			CALLED BY
 *  -----------------------------------------------------------------
 *  fvelo_normal_bc                      void
 *  fvelo_normal_lub_bc                  void
 *  fmesh_etch_bc                        void
 *  sdc_stefan_flow                      void
 *  sdc_stefan_volume_flow               void
 *  mass_flux_surface                  double
 *  vol_flux_surface                   double
 *  fvelo_normal_edge_bc                 void
 *  fvelo_normal_disc_bc                 void
 *  fvelo_tangent_edge_bc                void
 *  fvelo_tangent_3d                     void
 *  fvelo_tangential_bc                  void
 *  fvelo_tangential_solid_bc            void
 *  fvelo_tangential_ls_bc               void
 *  f_bc                                 void
 *  fvelo_slip_level                     void
 *  kin_bc                               void
 *  check_for_contact                    void
 *  load_surface_tension                 void
 *  elec_surf_stress                     void
 *  fn_dot_T                             void
 *  apply_repulsion                      void
 *  apply_vapor_recoil                   void
 *  flow_n_dot_T_hydro                   void
 *  flow_n_dot_T_var_density             void
 *  hydrostatic_n_dot_T                  void
 *  flow_n_dot_T_nobc                    void
 *  flow_n_dot_T_gradv                   void
 *  PSPG_consistency_bc                  void
 *  fapply_CA                            void
 *  fapply_var_CA                        void
 *  fapply_var_CA_user                   void
 *  evaluate_gibbs_criterion             int
 *  fapply_moving_CA                     void
 *  fapply_ST                            void
 *  apply_ST_scalar                      void
 *  apply_ST_3D                          void
 *  apply_ST_scalar_3D                   void
 *  apply_CA_FILL                        void
 *  ftmelt_bc                            void
 *  continuous_tangent_velocity          void
 *  continuous_normal_velocity           void
 *  discontinuous_velocity               void
 *  fnormal_stress_bc                    void
 *  fvelo_slip_electrokinetic_bc         void
 *  fvelo_electrokinetic_3d              void
 *  qside_directional                    void
 *  fvelo_airfilm                        void
 *
 ******************************************************************************/

/*  EXTERNAL VARIABLES  */

extern double dsigma_dx[DIM][MDE];

static double slip_coefficient(
    const double, const double, const double, const double, double *, const double, const double);

extern FSUB_TYPE dgemv_(char *TRANS,
                        int *M,
                        int *N,
                        double *alpha,
                        double *A,
                        int *LDA,
                        double *X,
                        int *INCX,
                        double *beta,
                        int *Y,
                        int *INCY);

/*
 *  Applies end slope nat'l SHEET_ENDSLOPE boundary condition on TENSION_SHEET boundary condition
 */

void apply_SES(double func[MAX_PDIM],
               double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
               struct elem_side_bc_struct *elem_side_bc,
               double sheet_tension,
               int inode,
               double XR,
               double YR,
               int id,
               int iconnect_ptr,
               int id_side) {
  int b, i, j, I, local_id, ldof, i_basis, ldofm;
  int eqn, var;
  int nf, ie;

  double XN, YN;
  double dY_dX, sign;

  double dX_dxi = 0., dY_dxi = 0.;
  double dX_dS, detJ;
  double dX_dS_dmesh[MDE][DIM];
  double ddetJ_dmesh[MDE][DIM];
  double dY_dX_dmesh[MDE][DIM];

  eqn = MESH_DISPLACEMENT1;

  i_basis = 1 - id_side % 2;

  ldof = ei[pg->imtrx]->ln_to_first_dof[eqn][id];

  XN = Coor[0][inode] + *esp->d[0][ldof];
  YN = Coor[1][inode] + *esp->d[1][ldof];

  memset(dY_dX_dmesh, 0, sizeof(double) * MDE * DIM);

  /*      if( fv->snormal[1] < 0.0 ) sheet_tension = -1.0 * fabs( sheet_tension);
          else sheet_tension = fabs(sheet_tension);        */

  /* this is a total kludge to change the sign of
   * the curvature convention so that a positive pressure
   * will always cause the surface to bulge away from the
   * the domain for boundaries in which the domain is
   * above or below the surface */

  /* Determine if sh_tens is active on neighboring shell block */
  if (num_shell_blocks != 0)
    nf = num_elem_friends[ei[pg->imtrx]->ielem];
  else
    nf = 0;

  /* If so, set up assembly to include variable shell tension */
  if (nf == 1 && upd->vp[pg->imtrx][SHELL_TENSION]) {
    ie = Index_Solution(inode, SHELL_TENSION, 0, 0, -1, pg->imtrx);
    sheet_tension = x_static[ie];
  }

  if (XN > XR) {
    dY_dX = (YN - YR) / (XN - XR);
    sign = -1.0;
    dY_dX_dmesh[ldof][0] = -dY_dX / (XN - XR);
    dY_dX_dmesh[ldof][1] = 1.0 / (XN - XR);
  } else {
    dY_dX = (YR - YN) / (XR - XN);
    sign = 1.0;
    dY_dX_dmesh[ldof][0] = dY_dX / (XR - XN);
    dY_dX_dmesh[ldof][1] = -1.0 / (XR - XN);
  }

  for (i = 0; i < (int)elem_side_bc->num_nodes_on_side; i++) {
    local_id = (int)elem_side_bc->local_elem_node_id[i];
    I = Proc_Elem_Connect[iconnect_ptr + local_id];
    ldof = ei[pg->imtrx]->ln_to_first_dof[eqn][local_id];

    if (ldof >= 0) {
      dX_dxi += (Coor[0][I] + *esp->d[0][ldof]) * bf[eqn]->dphidxi[ldof][i_basis];
      dY_dxi += (Coor[1][I] + *esp->d[1][ldof]) * bf[eqn]->dphidxi[ldof][i_basis];
    }
  }

  if (dX_dxi < 0.0)
    detJ = -sqrt(dX_dxi * dX_dxi + dY_dxi * dY_dxi);
  else
    detJ = sqrt(dX_dxi * dX_dxi + dY_dxi * dY_dxi);

  dX_dS = dX_dxi / detJ;

  memset(dX_dS_dmesh, 0, DIM * MDE * sizeof(double));
  memset(ddetJ_dmesh, 0, DIM * MDE * sizeof(double));

  for (i = 0; i < (int)elem_side_bc->num_nodes_on_side; i++) {
    local_id = (int)elem_side_bc->local_elem_node_id[i];
    ldof = ei[pg->imtrx]->ln_to_first_dof[eqn][local_id];

    if (ldof >= 0) {
      ddetJ_dmesh[ldof][0] = dX_dxi * bf[eqn]->dphidxi[ldof][i_basis] / detJ;
      ddetJ_dmesh[ldof][1] = dY_dxi * bf[eqn]->dphidxi[ldof][i_basis] / detJ;

      dX_dS_dmesh[ldof][0] = bf[eqn]->dphidxi[ldof][i_basis] / detJ;
      dX_dS_dmesh[ldof][0] += -(dX_dxi / detJ / detJ) * ddetJ_dmesh[ldof][0];

      dX_dS_dmesh[ldof][1] = -(dX_dxi / detJ / detJ) * ddetJ_dmesh[ldof][1];
    }
  }

  func[0] = sign * sheet_tension * dY_dX * dX_dS;

  if (af->Assemble_Jacobian) {
    for (j = 0; j < (int)elem_side_bc->num_nodes_on_side; j++) {
      local_id = (int)elem_side_bc->local_elem_node_id[j];

      for (b = 0; b < ei[pg->imtrx]->ielem_dim; b++) {
        var = MESH_DISPLACEMENT1 + b;
        ldof = ei[pg->imtrx]->ln_to_first_dof[var][local_id];

        if (ldof >= 0 && pd->v[pg->imtrx][var]) {

          d_func[0][var][ldof] = dY_dX_dmesh[ldof][b] * dX_dS;
          d_func[0][var][ldof] += dY_dX * dX_dS_dmesh[ldof][b];
          d_func[0][var][ldof] *= sign * sheet_tension;
        }
      }
    }

    ldofm = ei[pg->imtrx]->ln_to_first_dof[R_SHELL_TENSION][id];
    var = SHELL_TENSION;
    if (nf == 1 && upd->vp[pg->imtrx][var]) {
      d_func[0][var][ldofm] = sign * dY_dX * dX_dS;
    }
  }

  return;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void ls_attach_bc(double func[DIM],
                  double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                  const double v_attach) {
  int i, j, a, dim = pd->Num_Dim;
  double length_scale;
  double dot_prod;
  double d_dot_prod_dF[MDE];
  double n_i[DIM] = {0., 0., 0.};

  length_scale = dim == 2 ? 2.0 * fv->sdet : 2.0 * pow(fv->sdet, 0.5);

  load_lsi(2.0 * length_scale);

  if (lsi->near) {
    memset(d_dot_prod_dF, 0, MDE * sizeof(double));

    for (i = 0; i < dim; i++)
      n_i[i] = fv->grad_F[i];
    normalize_really_simple_vector(n_i, dim);

    for (a = 0, dot_prod = 0; a < VIM; a++) {
      dot_prod += n_i[a] * fv->snormal[a];

      for (j = 0; j < ei[pg->imtrx]->dof[LS]; j++) {
        d_dot_prod_dF[j] += lsi->d_normal_dF[a][j] * fv->snormal[a];
      }
    }

    func[0] -= v_attach * lsi->delta * pow(dot_prod, 4.0) * length_scale;

    if (af->Assemble_Jacobian) {
      int var;

      var = LS;

      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_func[0][var][j] -= v_attach * lsi->d_delta_dF[j] * pow(dot_prod, 4.0) * length_scale;
        d_func[0][var][j] -=
            v_attach * lsi->delta * 4.0 * pow(dot_prod, 3.0) * d_dot_prod_dF[j] * length_scale;
      }
    }
  }
  return;
}

/*****************************************************************************/
/****************************************************************************/
void fvelo_normal_bc(double func[DIM],
                     double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                     const double vnormal,            /* normal velocity */
                     const int contact,               /* should we allow ls contact */
                     const double x_dot[MAX_PDIM],    /* Bad name, says Phil!
                                                       * -mesh velocity vector */
                     const dbl tt,                    /* parameter to vary time integration from
                                                         explicit (tt = 1) to implicit (tt = 0) */
                     const dbl dt,                    /* current value of the time step */
                     const int type,                  /*  bc id  */
                     const double length,             /*  interface zone half-width  */
                     const double shift,              /*  interface zone shift for fake bc */
                     const double leak_angle_degrees) /*contact angle to start gas leaking */

/***********************************************************************
 *
 * fvelo_normal_bc():
 *
 *  Function which evaluates the expression specifying the
 *  kinematic boundary condition at a quadrature point on a side
 *  of an element.
 *
 *         func =   - n . (vnormal) + n . (v - xdot)
 *
 *  The boundary conditions, KINEMATIC_PETROV_BC,  KINEMATIC_BC:
 *  VELO_NORMAL_BC, and KINEMATIC_COLLOC_BC
 *  employ this function. For internal boundaries, vnormal is set
 *  to zero,
 *  and this function is evaluated on both sides of the interface to
 *  establish mass continuity across the interface.
 *
 *  Note: this bc is almost exactly the same as fvelo_normal_bc_disc
 *        except that the latter bc contains the density as a
 *        multiplicative constant.
 *
 * Input:
 *
 *  vnormal = specified on the bc card as the first float
 *  v       = mass average velocity at the guass point
 *  xdot    = velocity of the mesh (i.e., the interface if the mesh
 *            is pegged at the inteface)
 *
 * Output:
 *
 *  func[0] = value of the function mentioned above
 *  d_func[0][varType][lvardof] =
 *              Derivate of func[0] wrt
 *              the variable type, varType, and the local variable
 *              degree of freedom, lvardof, corresponding to that
 *              variable type.
 *
 *   Author: P. R. Schunk    (3/24/94)
 *   Revised: RAC  (5/24/95)
 ********************************************************************/
{
  int j, kdir, var, p;
  double phi_j, vcontact;
  double penalty = 1.0, d_penalty_dF = 0.0;
  double dot_prod;
  const double cos_leak = cos((180 - leak_angle_degrees) * (M_PIE / 180));
  const double leak_width = sin((180 - leak_angle_degrees) * (M_PIE / 180)) * sin(10 * M_PIE / 180);

  if (af->Assemble_LSA_Mass_Matrix) {
    for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
      var = MESH_DISPLACEMENT1 + kdir;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          d_func[0][var][j] -= phi_j * fv->snormal[kdir];
        }
      }
    }
    return;
  }

  vcontact = 0.;

  if (contact && (type != VELO_NORMAL_LS_BC && type != VELO_NORMAL_LS_PETROV_BC &&
                  type != VELO_NORMAL_LS_COLLOC_BC)) {
    double length_scale;
    double sign_gas = mp->mp2nd->densitymask[1] - mp->mp2nd->densitymask[0];

#if 0
      length_scale = pd->Num_Dim == 2 ? fv->sdet : pow( fv->sdet , 0.5);
#else
    length_scale = pd->Num_Dim == 2 ? fv->sdet : pow(fv->sdet, 0.5);
    if (upd->CoordinateSystem == CYLINDRICAL || upd->CoordinateSystem == SWIRLING) {
      length_scale = pd->Num_Dim == 2 ? fv->sdet / fv->h3 : pow(fv->sdet / fv->h3, 0.5);
    }
#endif

    load_lsi(ls->Length_Scale);

    /* only apply if in vicinity of boundary and OUTSIDE interface */
    if (sign_gas * lsi->H > 0.0 && fabs(fv->F) < 2.0 * length_scale * ls->Contact_Tolerance) {
      dot_prod = 0.0;
      for (p = 0; p < pd->Num_Dim; p++)
        dot_prod += fv->snormal[p] * lsi->normal[p];

      if (dot_prod > cos_leak) {
        func[0] = 0.0;
        return;
      }
    }
  }

  /* Calculate the residual contribution	*/
  func[0] = -vnormal - vcontact;
  for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
    func[0] += (fv->v[kdir] - x_dot[kdir]) * fv->snormal[kdir];
  }

  /*  modification for "fake" gas outlet  */
  if (contact && (type == VELO_NORMAL_LS_BC || type == VELO_NORMAL_LS_PETROV_BC ||
                  type == VELO_NORMAL_LS_COLLOC_BC)) {
    const double penalty_gas = 0.0, penalty_liq = 1.;
    const double gas_log = log(1. / BIG_PENALTY), liq_log = 0.;
    const int log_scale = 0;
    double factor = 0, Hfcn, d_Hfcn_dF, d_factor_dF = 0, F_prime;
    int visc_sens = mp->mp2nd->viscositymask[1] - mp->mp2nd->viscositymask[0];
    double dvisc_sens = (double)visc_sens;

    Hfcn = 0.;
    d_Hfcn_dF = 0.;
    F_prime = fv->F / length + dvisc_sens * shift;
    if (F_prime < -1.0) {
      Hfcn = 0.;
      d_Hfcn_dF = 0.;
    } else if (F_prime < 1.) {
      Hfcn = 0.5 * (1. + F_prime + sin(M_PIE * F_prime) / M_PIE);
      d_Hfcn_dF = 0.5 * (1. + cos(M_PIE * F_prime)) / length;
    } else {
      Hfcn = 1.;
      d_Hfcn_dF = 0.;
    }
    switch (visc_sens) {
    case 1:
      factor = 1. - Hfcn;
      d_factor_dF = -d_Hfcn_dF;
      break;
    case -1:
      factor = Hfcn;
      d_factor_dF = d_Hfcn_dF;
      break;
    }
    if (func[0] < 0.0 && 0) {
      factor = 1.;
      d_factor_dF = 0.;
    }
    if (log_scale) {
      penalty = exp((liq_log - gas_log) * factor + gas_log);
      d_penalty_dF = penalty * (liq_log - gas_log) * d_factor_dF * (-dvisc_sens);
    } else {
      penalty = (penalty_liq - penalty_gas) * factor + penalty_gas;
      d_penalty_dF = (penalty_liq - penalty_gas) * d_factor_dF * (-dvisc_sens);
    }
    if (1 && fabs(F_prime) < 1.0) {
      dot_prod = 0.0;
      for (p = 0; p < pd->Num_Dim; p++)
        dot_prod += fv->snormal[p] * lsi->normal[p];

      if (fabs(dot_prod - cos_leak) < leak_width) {
        F_prime = (dot_prod - cos_leak) / leak_width;
        Hfcn = 0.5 * (1. + F_prime + sin(M_PIE * F_prime) / M_PIE);
        penalty = Hfcn * penalty_gas + (1. - Hfcn) * penalty;
        d_penalty_dF *= (1. - Hfcn);
      } else if (dot_prod >= (cos_leak + leak_width)) {
        penalty = penalty_gas;
        d_penalty_dF = 0.;
      }
    }
    if (lsi->near && 0)
      fprintf(stderr, "vn_ls %g %g %g %g\n", fv->x[0], penalty, factor, func[0]);
  } /* end modifications for fake gas outlet */

  func[0] *= penalty;

  if (af->Assemble_Jacobian) {

    for (kdir = 0; kdir < pd->Num_Dim; kdir++) {

      for (p = 0; p < pd->Num_Dim; p++) {
        var = MESH_DISPLACEMENT1 + p;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            d_func[0][var][j] +=
                penalty * (fv->v[kdir] - x_dot[kdir]) * fv->dsnormal_dx[kdir][p][j];
            if (TimeIntegration != 0 && p == kdir) {
              d_func[0][var][j] +=
                  penalty * (-(1. + 2. * tt) * phi_j / dt) * fv->snormal[kdir] * delta(p, kdir);
            }
          }
        }
      }

      var = VELOCITY1 + kdir;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          d_func[0][var][j] += penalty * phi_j * fv->snormal[kdir];
        }
      }

      var = PVELOCITY1 + kdir;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          d_func[0][var][j] += penalty * phi_j * fv->snormal[kdir];
        }
      }

      if (contact) {
        var = ls->var;
        if ((type == VELO_NORMAL_LS_BC || type == VELO_NORMAL_LS_PETROV_BC ||
             type == VELO_NORMAL_LS_COLLOC_BC) &&
            pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_func[0][var][j] +=
                (fv->v[kdir] - x_dot[kdir]) * fv->snormal[kdir] * d_penalty_dF * bf[var]->phi[j];
          }
        }
      }

    } /* for: kdir */
  }   /* end of if Assemble_Jacobian */

} /* END of routine fvelo_normal_bc  */
/*****************************************************************************/
/*****************************************************************************/
/****************************************************************************/

void fvelo_normal_lub_bc(double func[DIM],
                         double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                         const int id_side,            /* Side ID */
                         const double x_dot[MAX_PDIM], /* Bad name, says Phil!
                                                        * -mesh velocity vector */
                         const double tt,              /* parameter to vary time integration from
                                                          explicit (tt = 1) to implicit (tt = 0) */
                         const double dt,              /* current value of the time step */
                         double xi[DIM],               /* Local stu coordinates */
                         const Exo_DB *exo,            /* ExodusII database struct pointer */
                         const double param[])         /* Flux parameters  */
/***********************************************************************
 *
 * fvelo_normal_lub_bc():
 *
 *  Function which evaluates the expression specifying the
 *  boundary normal motion  at a quadrature point on a side
 *  of an element due to flux to/from shell lubrication layer.
 *
 *         func  =   - n . (lubflux) + n . (v - xdot)
 *
 *  The boundary conditions VELO_NORMAM_LUB_BC
 *  employ this function. vnormal is typically dictated by current surface
 *  orientation as well as surface chemistry.
 *
 *  Note: this function initially is cloned from fvelo_normal_bc
 *
 * Input:
 *
 *  param[] = specified on the bc card as the first integer
 *
 *
 * Output:
 *
 *  func[0] = value of the function mentioned above
 *  d_func[0][varType][lvardof] =
 *              Derivate of func[0] wrt
 *              the variable type, varType, and the local variable
 *              degree of freedom, lvardof, corresponding to that
 *              variable type.
 *
 *   Author: Kristianto Tjiptowidjojo    (03/11/2021)
 ********************************************************************/

{
  int j, kdir, var, p, q;
  int *n_dof = NULL;
  int dof_map[MDE];
  double phi_j;
  double bound_normal[DIM], d_bd_normal_dx[DIM][DIM][MDE];
  int el1, el2, nf, nbr_type, nbr_dim;
  double lubflux = param[0];

  /* Basic data for local element */
  el1 = ei[pg->imtrx]->ielem;

  nf = num_elem_friends[el1];
  if (nf == 0) {
    GOMA_EH(GOMA_ERROR, "Where is my bulk element neighbor?");
  } else if (nf == 1) {
    el2 = elem_friends[el1][0];

    /* Determine if this element is bulk with shell friend(s) */
    nbr_type = Elem_Type(exo, el2);
    nbr_dim = elem_info(NDIM, nbr_type);

    /* NOTE: this will not work for an edge of a 3D element yet */
    /* return if not on the shell element */
    if (ei[pg->imtrx]->ielem_dim > nbr_dim)
      return;
  } else {
    GOMA_EH(GOMA_ERROR, "Not ready for multiple neighbors yet...");
  }
  /* Save the boundary normal vector */
  for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
    bound_normal[kdir] = fv->snormal[kdir];
    for (j = 0; j < pd->Num_Dim; j++) {
      for (var = 0; var < pd->Num_Dim; var++) {
        d_bd_normal_dx[kdir][j][var] = fv->dsnormal_dx[kdir][j][var];
      }
    }
  }
  /*
   * Proceed if we are on the shell element
   * Prepare geometry
   */
  n_dof = (int *)array_alloc(1, MAX_VARIABLE_TYPES, sizeof(int));
  lubrication_shell_initialize(n_dof, dof_map, id_side, xi, exo, 0);

  /* Calculate the flow rate and its sensitivties */

  calculate_lub_q_v(R_LUBP, tran->time_value, dt, xi, exo);

  /***** CALCULATE RESIDUAL CONTRIBUTION ********************/
  func[0] = lubflux;
  for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
    func[0] += (LubAux->q[kdir] - LubAux->H * x_dot[kdir]) * bound_normal[kdir];
  }

  /***** CALCULATE JACOBIAN CONTRIBUTION ********************/
  if (af->Assemble_Jacobian) {

    for (kdir = 0; kdir < pd->Num_Dim; kdir++) {

      for (p = 0; p < pd->Num_Dim; p++) {
        var = MESH_DISPLACEMENT1 + p;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            d_func[0][var][j] +=
                (LubAux->q[kdir] - LubAux->H * x_dot[kdir]) * d_bd_normal_dx[kdir][p][j];
            d_func[0][var][j] +=
                (LubAux->dq_dx[kdir][p][j] - LubAux->dH_dmesh[p][j] * x_dot[kdir]) *
                bound_normal[kdir];
            if (TimeIntegration != 0 && p == kdir) {
              d_func[0][var][j] +=
                  (-(1. + 2. * tt) * phi_j / dt) * bound_normal[kdir] * delta(p, kdir);
            }
          }
        }
      }
    } /* for: kdir */

    var = LUBP;
    if (n_dof[var] > 0) {
      double grad_phi_j[DIM], grad_II_phi_j[DIM], d_grad_II_phi_j_dmesh[DIM][DIM][MDE];
      for (j = 0; j < n_dof[var]; j++) {
        phi_j = bf[var]->phi[j];
        /* Prepare basis functions */
        ShellBF(var, j, &phi_j, grad_phi_j, grad_II_phi_j, d_grad_II_phi_j_dmesh,
                n_dof[MESH_DISPLACEMENT1], dof_map);
        for (p = 0; p < pd->Num_Dim; p++) {
          d_func[0][var][j] += LubAux->dq_dp2[p][j] * phi_j * grad_II_phi_j[p] * bound_normal[p];
          for (q = 0; q < pd->Num_Dim; q++) {
            d_func[0][var][j] +=
                LubAux->dq_dgradp[p][q][j] * grad_II_phi_j[q] * grad_II_phi_j[p] * bound_normal[p];
          }
        }
      }
    }

    var = SHELL_FILMP;
    if (n_dof[var] > 0) {
      double grad_phi_j[DIM], grad_II_phi_j[DIM], d_grad_II_phi_j_dmesh[DIM][DIM][MDE];
      for (j = 0; j < n_dof[var]; j++) {
        phi_j = bf[var]->phi[j];
        /* Prepare basis functions */
        ShellBF(var, j, &phi_j, grad_phi_j, grad_II_phi_j, d_grad_II_phi_j_dmesh,
                n_dof[MESH_DISPLACEMENT1], dof_map);
        for (p = 0; p < pd->Num_Dim; p++) {
          d_func[0][var][j] += LubAux->dq_dp2[p][j] * phi_j * grad_II_phi_j[p] * bound_normal[p];
          for (q = 0; q < pd->Num_Dim; q++) {
            d_func[0][var][j] +=
                LubAux->dq_dgradp[p][q][j] * grad_II_phi_j[q] * grad_II_phi_j[p] * bound_normal[p];
          }
        }
      }
    }

  } /* end of if Assemble_Jacobian */

  safe_free((void *)n_dof);

} /* END of routine fvelo_normal_lub_bc  */
/*****************************************************************************/
/*****************************************************************************/
/****************************************************************************/

void fmesh_etch_bc(double *func,
                   double d_func[MAX_VARIABLE_TYPES + MAX_CONC],
                   const int etch_plane,      /* Etch plane */
                   const int id,              /* Local node ID */
                   const dbl x_dot[MAX_PDIM], /* Mesh velocity */
                   const dbl tt,              /* parameter to vary time integration from
                                                 explicit (tt = 1) to implicit (tt = 0) */
                   const dbl dt)              /* current value of the time step */
/***********************************************************************
 *
 * fmesh_etch_bc():
 *
 *  Function which evaluates the expression specifying the
 *  boundary normal motion  at a quadrature point on a side
 *  of an element due to etching reaction.
 *
 *         func =  - etch_rate + n . xdot
 *
 *  The boundary conditions MOVING_PLANE_ETCH_BC
 *  employ this function. vnormal is typically dictated by current surface
 *  orientation as well as surface chemistry.
 *
 *  Note: this function initially is cloned from fvelo_normal_bc
 *
 * Input:
 *
 *  etch_plane = specified on the bc card as the first integer
 *
 *
 * Output:
 *
 *  func[0] = value of the function mentioned above
 *  d_func[0][varType][lvardof] =
 *              Derivate of func[0] wrt
 *              the variable type, varType, and the local variable
 *              degree of freedom, lvardof, corresponding to that
 *              variable type.
 *
 *   Author: Kristianto Tjiptowidjojo    (02/27/2017)
 ********************************************************************/
{
  int a, b, w;
  int dim = pd->Num_Dim;
  int var;

  double etch_rate = 0.0;
  double d_etch_rate_d_C[2] = {0.0};
  double etch_rate_sens[MAX_CONC] = {0.0};

  /* Right now we only consider KOH wet etching of silicon at 100 plane */
  if (etch_plane == 100) {
    /* Get etch rate */
    double rho_H2O = fv->c[0];
    double rho_KOH = fv->c[1];
    etch_rate = calc_KOH_Si_etch_rate_100(rho_H2O, rho_KOH, d_etch_rate_d_C);
    etch_rate_sens[0] = d_etch_rate_d_C[0];
    etch_rate_sens[1] = d_etch_rate_d_C[1];
  }

  /* Set the residual */
  for (a = 0; a < dim; a++) {
    *func += fv->snormal[a] * x_dot[a];
  }
  *func -= etch_rate;

  /***** NOW FIND SENSITIVITIES *****/

  /* Mesh sensitivities */
  for (b = 0; b < dim; b++) {
    var = MESH_DISPLACEMENT1 + b;
    for (a = 0; a < dim; a++) {
      d_func[var] += fv->dsnormal_dx[a][b][id] * x_dot[a];
      d_func[var] += fv->snormal[a] * (1.0 + 2.0 * tt) / dt * delta(a, b);
    }
  }

  /* Concentration sensitivities */
  var = MASS_FRACTION;
  if (pd->v[pg->imtrx][var]) {
    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      d_func[MAX_VARIABLE_TYPES + w] -= etch_rate_sens[w];
    }
  }

  return;
} /* END of routine fmesh_etch_bc */
/*****************************************************************************/
/*****************************************************************************/
/****************************************************************************/

void fvelo_tangential_ls_bc(double func[DIM],
                            double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                            const double vtangent,        /* vtangent velocity */
                            const double x_dot[MAX_PDIM], /* Bad name, says Phil!
                                                           * -mesh velocity vector */
                            const dbl tt,                 /* parameter to vary time integration from
                                                             explicit (tt = 1) to implicit (tt = 0) */
                            const dbl dt,                 /* current value of the time step */
                            const double length,          /*  interface zone half-width  */
                            const double shift,           /*  interface zone shift for fake bc */
                            const double leak_angle_degrees) /*contact angle to start gas leaking */

/***********************************************************************
 *
 * fvelo_tangential_ls_bc():
 *
 *  Function which evaluates the expression specifying the
 *  tangential boundary condition at a quadrature point on a side
 *  of an element.
 *
 *         func =   - t . (vnormal) + t . (v - xdot)
 *
 *
 *  Note: this bc is almost exactly the same as fvelo_normal_bc
 *
 * Input:
 *
 *  vtangent= specified on the bc card as the first float
 *  v       = mass average velocity at the guass point
 *  xdot    = velocity of the mesh (i.e., the interface if the mesh
 *            is pegged at the inteface)
 *
 * Output:
 *
 *  func[0] = value of the function mentioned above
 *  d_func[0][varType][lvardof] =
 *              Derivate of func[0] wrt
 *              the variable type, varType, and the local variable
 *              degree of freedom, lvardof, corresponding to that
 *              variable type.
 *
 *   Author: Slah Jendoubi   (6/28/2012)
 *   Revised:
 ********************************************************************/
{
  int j, kdir, var, p;
  double phi_j, vcontact;
  double penalty = 1.0, d_penalty_dF = 0.0;
  const double penalty_gas = 0.0, penalty_liq = 1.;
  const double gas_log = log(1. / BIG_PENALTY), liq_log = 0.;
  const int log_scale = 0;
  double factor = 0, Hfcn, d_Hfcn_dF, d_factor_dF = 0, F_prime;
  int visc_sens = mp->mp2nd->viscositymask[1] - mp->mp2nd->viscositymask[0];
  double dvisc_sens = (double)visc_sens;

  if (af->Assemble_LSA_Mass_Matrix) {
    for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
      var = MESH_DISPLACEMENT1 + kdir;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          d_func[0][var][j] -= phi_j * fv->stangent[0][kdir];
        }
      }
    }
    return;
  }

  vcontact = 0.;

  /* Calculate the residual contribution	*/
  func[0] = -vtangent - vcontact;
  for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
    func[0] += (fv->v[kdir] - x_dot[kdir]) * fv->stangent[0][kdir];
  }

  /*  modification for "fake" gas outlet  */

  Hfcn = 0.;
  d_Hfcn_dF = 0.;
  F_prime = fv->F / length + dvisc_sens * shift;
  if (F_prime < -1.0) {
    Hfcn = 0.;
    d_Hfcn_dF = 0.;
  } else if (F_prime < 1.) {
    Hfcn = 0.5 * (1. + F_prime + sin(M_PIE * F_prime) / M_PIE);
    d_Hfcn_dF = 0.5 * (1. + cos(M_PIE * F_prime)) / length;
  } else {
    Hfcn = 1.;
    d_Hfcn_dF = 0.;
  }
  switch (visc_sens) {
  case 1:
    factor = 1. - Hfcn;
    d_factor_dF = -d_Hfcn_dF;
    break;
  case -1:
    factor = Hfcn;
    d_factor_dF = d_Hfcn_dF;
    break;
  }
  if (log_scale) {
    penalty = exp((liq_log - gas_log) * factor + gas_log);
    d_penalty_dF = penalty * (liq_log - gas_log) * d_factor_dF * (-dvisc_sens);
  } else {
    penalty = (penalty_liq - penalty_gas) * factor + penalty_gas;
    d_penalty_dF = (penalty_liq - penalty_gas) * d_factor_dF * (-dvisc_sens);
  }

  func[0] *= penalty;

  if (af->Assemble_Jacobian) {

    for (kdir = 0; kdir < pd->Num_Dim; kdir++) {

      for (p = 0; p < pd->Num_Dim; p++) {
        var = MESH_DISPLACEMENT1 + p;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            d_func[0][var][j] +=
                penalty * (fv->v[kdir] - x_dot[kdir]) * fv->dstangent_dx[0][kdir][p][j];
            if (TimeIntegration != 0 && p == kdir) {
              d_func[0][var][j] +=
                  penalty * (-(1. + 2. * tt) * phi_j / dt) * fv->stangent[0][kdir] * delta(p, kdir);
            }
          }
        }
      }

      var = VELOCITY1 + kdir;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          d_func[0][var][j] += penalty * phi_j * fv->stangent[0][kdir];
        }
      }

      var = PVELOCITY1 + kdir;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          d_func[0][var][j] += penalty * phi_j * fv->stangent[0][kdir];
        }
      }

      var = ls->var;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_func[0][var][j] +=
              (fv->v[kdir] - x_dot[kdir]) * fv->stangent[0][kdir] * d_penalty_dF * bf[var]->phi[j];
        }
      }

    } /* for: kdir */
  }   /* end of if Assemble_Jacobian */

} /* END of routine fvelo_tangential_ls_bc  */

/*****************************************************************************/

double sdc_stefan_flow(JACOBIAN_VAR_DESC_STRUCT *func_jac,
                       BOUNDARY_CONDITION_STRUCT *bc,
                       int ip,
                       ELEM_SIDE_BC_STRUCT *elem_side_bc,
                       const double x_dot[MAX_PDIM],
                       const dbl time,
                       const dbl tt,
                       const dbl dt,
                       const int intf_id)

/***********************************************************************
 *
 * sdc_stefan_flow():
 *
 *  Function which evaluates the expression specifying the
 *  kinematic boundary condition at a quadrature point on a side
 *  of an element.
 *
 *         func =   n . rho * (v - xdot) + Sum_k_1_to_N[W_k * S_k]
 *
 *  This function is evaluated on just one side of the interface
 *  to establish the link with Stefan flow
 *
 *
 * Input:
 *
 *  v       = mass average velocity at the current surface
 *            quadrature point
 *  xdot    = velocity of the mesh (i.e., the interface if the mesh
 *            is pegged at the inteface)
 *  tt   = parameter to vary time integration from
 *                     explicit (tt = 1) to implicit (tt = 0)
 *  dt      = current value of the time step
 *
 * Output:
 *
 *  func[0] = value of the function mentioned above
 *  d_func[0][varType][lvardof] =
 *              Derivate of func[0] wrt
 *              the variable type, varType, and the local variable
 *              degree of freedom, lvardof, corresponding to that
 *              variable type.
 *
 ********************************************************************/
{
  int pos, pos_mp, k, lvdesc, index_is;
  double func_value, *mw, *Sk, tmp;
  MATRL_PROP_STRUCT *mp2, *mp_b, *mp_a;
  BOUNDARY_CONDITION_STRUCT *bc_rxn;
  int have_T = FALSE, speciesVT, wspec = 0;
  INTERFACE_SOURCE_STRUCT *is, **is_hdl = 0;
  VARIABLE_DESCRIPTION_STRUCT *vd;

  if (af->Assemble_LSA_Mass_Matrix)
    return 0.0;

  mw = mp->molecular_weight;

  /*
   * Find the identity of the second material, mp2, the material on
   * the side of the boundary which this one sided boundary condition
   * is not being applied to.
   */
  mp2 = mp_glob[bc->BC_matrl_index_2];
  if (mp == mp2) {
    mp2 = mp_glob[bc->BC_matrl_index_1];
    if (mp == mp2) {
      GOMA_EH(GOMA_ERROR, "cant find second material");
    }
  }

  /*
   * Obtain the pointer to the boundary condition that contains
   * the reaction rate information for this boundary condition
   */
  bc_rxn = BC_Types + bc->BC_Data_Int[2];

  /*
   * We have two materials, mp and mp2. However, we need to find
   * their proper odering in the interfacial reaction source
   * term. This is obtained from the reaction boundary condition
   * card.
   */
  mp_a = mp_glob[map_mat_index(bc_rxn->BC_Data_Int[1])];
  mp_b = mp_glob[map_mat_index(bc_rxn->BC_Data_Int[2])];

  /*
   * Calculate the residual contribution, possibly branch on
   * the global value of the species unknown type.
   *
   *  First calculate the interfacial mass flux component
   *  We calculate the Jacobian contributions here too.
   *  They are storred in func_jac.
   */
  func_value = mass_flux_surface(func_jac, x_dot, tt, dt);

  /*
   * Process the Source term
   */
  switch (bc->Storage_ID) {
  case VL_EQUIL_PRXN_BC:
    is_hdl =
        (INTERFACE_SOURCE_STRUCT **)side_qp_storage_findalloc(VL_EQUIL_PRXN_BC, ip, elem_side_bc);
    break;
  case IS_EQUIL_PRXN_BC:
    is_hdl =
        (INTERFACE_SOURCE_STRUCT **)side_qp_storage_findalloc(IS_EQUIL_PRXN_BC, ip, elem_side_bc);
    break;
  default:
    GOMA_EH(GOMA_ERROR, "ERROR");
  }
  if (*is_hdl == NULL) {
    *is_hdl = alloc_struct_1(INTERFACE_SOURCE_STRUCT, Num_Interface_Srcs);
  }
  is = *is_hdl;

  /*
   * Find the starting location of the species unknowns for the
   * current material, pos_mp, in the
   * interfacial source term vector
   */
  pos_mp = 0;
  vd = is[intf_id].Var_List[0];
  if (vd->MatID != mp->MatID) {
    if (vd->MatID != mp2->MatID) {
      GOMA_EH(GOMA_ERROR, "unclear materials -> Matid = -1?");
    }
    k = mp2->Num_Species;
    vd = is[intf_id].Var_List[k];
    if (vd->MatID != mp->MatID) {
      GOMA_EH(GOMA_ERROR, "unclear interfacial source term ordering");
    }
    pos_mp = mp2->Num_Species;
  }

  for (wspec = 0; wspec < mp_a->Num_Species; wspec++) {
    if (!is[intf_id].Processed[wspec]) {
      is[intf_id].Processed[wspec] = TRUE;
      /*  bc_rxn = BC_Types + bc->BC_Data_Int[2]+wspec; */
      /*
       * Fill up the Var_Value[] list. Possibly change the type
       * of the species unknown vector at the same time.
       */
      is_masstemp_fillin(is, mp_a, mp_b, SPECIES_CONCENTRATION, time, intf_id);

      /*
       * Now, call the source routine and possibly do the
       * jacobian. The source term is defined as the
       *
       */
      is[intf_id].Do_Jac = af->Assemble_Jacobian;
      switch (bc->BC_Data_Int[1]) {
      case VL_EQUIL_PRXN_BC:
        source_vle_prxn(is, bc_rxn, mp_a, mp_b, have_T, intf_id);
        break;
      case IS_EQUIL_PRXN_BC:
        source_is_equil_prxn(is, bc_rxn, mp_a, mp_b, have_T, intf_id);
        break;
      default:
        GOMA_EH(GOMA_ERROR, "ERROR");
      }

      /*
       * We now need to possibly change the species variable type
       * for the source term and the dependent variable for the
       * Jacobian entries.
       */
      speciesVT = upd->Species_Var_Type;
      if (upd->Species_Var_Type == SPECIES_UNDEFINED_FORM) {
        speciesVT = SPECIES_MASS_FRACTION;
      }
      if (speciesVT != is[intf_id].SpeciesVT) {
        is_change1_speciesVT(is, wspec, 0, mp_a, speciesVT, time, intf_id);
        is_change1_speciesVT(is, wspec, mp_a->Num_Species, mp_b, speciesVT, time, intf_id);
        is_change1_speciesVT(is, mp_a->Num_Species + wspec, 0, mp_a, speciesVT, time, intf_id);
        is_change1_speciesVT(is, mp_a->Num_Species + wspec, mp_a->Num_Species, mp_b, speciesVT,
                             time, intf_id);
        convert_species_var(speciesVT, mp_a, is[intf_id].SpeciesVT, is[intf_id].Var_Value, time);
        convert_species_var(speciesVT, mp_b, is[intf_id].SpeciesVT,
                            is[intf_id].Var_Value + mp_a->Num_Species, time);
        is[intf_id].SpeciesVT = speciesVT;
        /*
         * Possibly Convert the format of what's held constant
         * during the partial derivatives wrt species variable
         * to one in which the sum of the mass fractions are
         * constant is a constraint
         */
        is_change1_lastspecies(is, wspec, 0, mp_a, intf_id);
        is_change1_lastspecies(is, wspec, mp_a->Num_Species, mp_b, intf_id);
        is_change1_lastspecies(is, mp_a->Num_Species + wspec, 0, mp_a, intf_id);
        is_change1_lastspecies(is, mp_a->Num_Species + wspec, mp_a->Num_Species, mp_b, intf_id);
      }
    }

    /*
     * Add in the source terms
     */
    Sk = is[intf_id].SourceTerm + pos_mp;
    func_value += mw[wspec] * Sk[wspec];
  }

  if (af->Assemble_Jacobian) {

    /*
     * Store the actual value of the function first
     */
    func_jac->Func_Value = func_value;

    /*
     * Make sure funcjac is big enough to accept all of the terms below
     * -> one time reallocs are better than multiple reallocs
     */
    jacobianVD_realloc(&func_jac, is[intf_id].Num_Terms + pd->Num_Dim,
                       pd->Num_Dim * ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]);

    /*
     * Next, the dependence of the source term sum on the
     * state variables for both the current material and the
     * material on the other side of the interface.
     */
    for (index_is = 0; index_is < is[intf_id].Num_Terms; index_is++) {
      /*
       * Look up what variable description this dependence is for
       */
      vd = is[intf_id].Var_List[index_is];
      /*
       * Use this vd structure to find the lvdesc index for
       * the corresponding variable description
       */
      lvdesc = ei[pg->imtrx]->VDindex_to_Lvdesc[vd->List_Index];
      if (lvdesc >= 0) {
        for (k = 0, tmp = 0.0; k < mp->Num_Species; k++) {
          pos = pos_mp + k;
          tmp += mw[k] * is[intf_id].JacMatrix[pos][index_is];
        }
        /*
         * Transfer the pertinent information to the
         *  Jacobian_Var_Desc structure
         */
        jacobianVD_addEntry(func_jac, lvdesc, tmp);
      }
    }
  } /* end of if Assemble_Jacobian */
  return func_value;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double sdc_stefan_volume_flow(JACOBIAN_VAR_DESC_STRUCT *func_jac,
                              BOUNDARY_CONDITION_STRUCT *bc,
                              int ip,
                              ELEM_SIDE_BC_STRUCT *elem_side_bc,
                              const double x_dot[MAX_PDIM],
                              const dbl tt,
                              const dbl dt)

/*************************************************************************
 *
 * sdc_stefan_volume_flow():
 *
 *  Function which evaluates the expression specifying the
 *  kinematic boundary condition at a quadrature point on a side
 *  of an element.
 *
 *         func =   n . rho * (v - xdot) + Sum_k_1_to_N[W_k * S_k]
 *
 *  This function is evaluated on just one side of the interface
 *  to establish the link with Stefan flow
 *
 *
 * Input:
 *
 *  v       = mass average velocity at the current surface
 *            quadrature point
 *  xdot    = velocity of the mesh (i.e., the interface if the mesh
 *            is pegged at the inteface)
 *  tt   = parameter to vary time integration from
 *                     explicit (tt = 1) to implicit (tt = 0)
 *  dt      = current value of the time step
 *
 * Output:
 *
 *  func[0] = value of the function mentioned above
 *  d_func[0][varType][lvardof] =
 *              Derivate of func[0] wrt
 *              the variable type, varType, and the local variable
 *              degree of freedom, lvardof, corresponding to that
 *              variable type.
 *
 ************************************************************************/
{
  return 0;
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

double mass_flux_surface(JACOBIAN_VAR_DESC_STRUCT *func_jac,
                         const double x_dot[MAX_PDIM],
                         const dbl tt,
                         const dbl dt)

/************************************************************************
 *
 *  mass_flux_surface():
 *
 *  Function which evaluates the expression specifying the
 *  mass flux through a surface relative to the velocity of that
 *  surface for the current material and side of the current element.
 *
 *         func =   n . rho * (v - xdot)
 *
 *  Input
 * ---------
 *     func_jac -> structure that must be malloced before calling.
 *     tt       -> Parameter in time derivative
 *     dt       -> current time step
 *     xdot[]   -> Current velocity of the interface
 *  Output
 * ---------
 *     return : value of func above
 *     func_jac -> Dependence of func on the dependent variables in the
 *                 problem, if needed.
 *
 ************************************************************************/
{
  int j, kdir, var, p, lvdesc;
  double *phi_ptr, func_value, n_dot_diffv;
  double tmp, tmp0, tmp1, tmp2, tmpt = 0.0;
  PROPERTYJAC_STRUCT *pj;
  VARIABLE_DESCRIPTION_STRUCT *vd;

  /*
   * Calculate the residual contribution, possibly branch on
   * the global value of the species unknown type.
   */
  n_dot_diffv = 0;
  for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
    n_dot_diffv += (fv->v[kdir] - x_dot[kdir]) * fv->snormal[kdir];
  }
  func_value = mp->density * n_dot_diffv;

  if (af->Assemble_Jacobian) {

    /*
     *  Zero the entries in func_jac
     *    -> This is a bottom level function
     */
    func_jac->Num_lvdesc = 0;
    func_jac->Num_lvdof = 0;
    /*
     * Store the actual value of the function first
     */
    func_jac->Func_Value = func_value;

    pj = mp->DensityJac;
    /*
     * Make sure funcjac is big enough to accept all of the terms below
     * -> one time reallocs are better than multiple reallocs
     */
    jacobianVD_realloc(&func_jac, pd->Num_Dim + pj->Num_Terms,
                       pd->Num_Dim * ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]);

    /*
     * First add in the mesh displacement dependencies -
     *  HKM Note: The number of nonzero terms is sparse wrt the
     *            total number of terms (e.g., 3 out of 9
     *            for biquadratic interpolations). Only the nodes
     *            which are actually on the element side have
     *            nonzero entries. Perhaps, we could make use of this
     *            feature to cut down the amount of work below.
     *  HKM Note: for LSA we don't add these in (best guess)
     */
    if (!af->Assemble_LSA_Mass_Matrix) {
      for (p = 0; p < pd->Num_Dim; p++) {
        var = MESH_DISPLACEMENT1 + p;
        phi_ptr = bf[var]->phi;
        if (ei[pg->imtrx]->Num_Lvdesc_Per_Var_Type[var] == 1) {
          lvdesc = ei[pg->imtrx]->Lvdesc_First_Var_Type[var];
          switch (pd->Num_Dim) {
          case 2:
            tmp0 = mp->density * (fv->v[0] - x_dot[0]);
            tmp1 = mp->density * (fv->v[1] - x_dot[1]);
            if (TimeIntegration) {
              tmpt = mp->density * (-(1. + 2. * tt) / dt) * fv->snormal[p];
            }
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              tmp = (fv->dsnormal_dx[0][p][j] * tmp0 + fv->dsnormal_dx[1][p][j] * tmp1);
              if (TimeIntegration) {
                tmp += tmpt * phi_ptr[j];
              }
              if (DOUBLE_NONZERO(tmp)) {
                jacobianLVDOF_addNewEntry(func_jac, var, j, tmp);
              }
              /* else { printf("there are zero terms\n"); } */
            }
            break;
          case 3:
            tmp0 = mp->density * (fv->v[0] - x_dot[0]);
            tmp1 = mp->density * (fv->v[1] - x_dot[1]);
            tmp2 = mp->density * (fv->v[2] - x_dot[2]);
            tmpt = mp->density * (-(1. + 2. * tt) / dt) * fv->snormal[p];
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              tmp = (fv->dsnormal_dx[0][p][j] * tmp0 + fv->dsnormal_dx[1][p][j] * tmp1 +
                     fv->dsnormal_dx[2][p][j] * tmp2);
              if (TimeIntegration) {
                tmp += tmpt * phi_ptr[j];
              }
              if (DOUBLE_NONZERO(tmp)) {
                jacobianLVDOF_addNewEntry(func_jac, var, j, tmp);
              }
            }
            break;
          default:
            GOMA_EH(GOMA_ERROR, "Not Covered");
          }
        }
      }
    }
    /*
     * Next, add in the velocity dependencies
     */
    if (VARIABLE_IN_THE_EB_SOLN_VECTOR(pd, pg->imtrx, VELOCITY1)) {
      for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
        tmp = fv->snormal[kdir] * mp->density;
        lvdesc = ei[pg->imtrx]->Lvdesc_First_Var_Type[VELOCITY1 + kdir];
        jacobianVD_addNewEntry(func_jac, lvdesc, tmp);
      }
    }

    /*
     * Next, the dependence of density on the state variables for
     * the current material.
     */
    for (j = 0; j < pj->Num_Terms; j++) {
      /*
       * Look up what variable description this dependence is for
       */
      vd = pj->Var_List[j];
      /*
       * Use this vd structure to find the lvdesc index for
       * the corresponding variable description
       */
      lvdesc = ei[pg->imtrx]->VDindex_to_Lvdesc[vd->List_Index];
      if (lvdesc >= 0) {
        /*
         * Transfer the pertinent information to the
         *  Jacobian_Var_Desc structure
         */
        jacobianVD_addNewEntry(func_jac, lvdesc, n_dot_diffv * pj->JacVector[j]);
      }
    }
  } /* end of if Assemble_Jacobian */
  return func_value;
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

double vol_flux_surface(JACOBIAN_VAR_DESC_STRUCT *func_jac,
                        const double x_dot[MAX_PDIM],
                        const dbl tt,
                        const dbl dt)

/************************************************************************
 *  vol_flux_surface():
 *
 *  Function which evaluates the expression specifying the
 *  mass averaged volume flux through a surface relative to the
 *  velocity of that surface for the current material and side of
 *  the current element.
 *
 *         func =   n . (v - xdot)
 *
 *  Input
 * ---------
 *     func_jac -> structure that must be malloced before calling.
 *     tt       -> Parameter in time derivative
 *     dt       -> current time step
 *     xdot[]   -> Current velocity of the interface
 *     v        -> average velocity in the material (usually it is
 *                 a mass average)
 *  Output
 * ---------
 *     return : value of func above
 *     func_jac -> Dependence of func on the dependent variables in the
 *                 problem, if needed.
 *
 ************************************************************************/
{
  int j, kdir, var, p, lvdesc;
  double *phi_ptr, func_value, n_dot_diffv;
  double tmp, tmp0, tmp1, tmp2, tmpt = 0.0;

  /*
   * Calculate the residual contribution, possibly branch on
   * the global value of the species unknown type.
   */
  n_dot_diffv = 0;
  for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
    n_dot_diffv += (fv->v[kdir] - x_dot[kdir]) * fv->snormal[kdir];
  }
  func_value = n_dot_diffv;

  if (af->Assemble_Jacobian) {

    /*
     *  Zero the entries in func_jac
     *    -> This is a bottom level function
     */
    func_jac->Num_lvdesc = 0;
    func_jac->Num_lvdof = 0;
    /*
     * Store the actual value of the function first
     */
    func_jac->Func_Value = func_value;

    /*
     * Make sure funcjac is big enough to accept all of the terms below
     * -> one time reallocs are better than multiple reallocs
     */
    jacobianVD_realloc(&func_jac, pd->Num_Dim,
                       pd->Num_Dim * ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]);

    /*
     * First add in the mesh displacement dependencies -
     *  HKM Note: The number of nonzero terms is sparse wrt the
     *            total number of terms (e.g., 3 out of 9
     *            for biquadratic interpolations). Only the nodes
     *            which are actually on the element side have
     *            nonzero entries. Perhaps, we could make use of this
     *            feature to cut down the amount of work below.
     *  HKM Note: for LSA we don't add these in (best guess)
     */
    if (!af->Assemble_LSA_Mass_Matrix) {
      for (p = 0; p < pd->Num_Dim; p++) {
        var = MESH_DISPLACEMENT1 + p;
        if (ei[pg->imtrx]->Num_Lvdesc_Per_Var_Type[var] == 1) {
          phi_ptr = bf[var]->phi;
          lvdesc = ei[pg->imtrx]->Lvdesc_First_Var_Type[var];
          switch (pd->Num_Dim) {
          case 2:
            tmp0 = (fv->v[0] - x_dot[0]);
            tmp1 = (fv->v[1] - x_dot[1]);
            if (TimeIntegration) {
              tmpt = (-(1. + 2. * tt) / dt) * fv->snormal[p];
            }
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              tmp = (fv->dsnormal_dx[0][p][j] * tmp0 + fv->dsnormal_dx[1][p][j] * tmp1);
              if (TimeIntegration) {
                tmp += tmpt * phi_ptr[j];
              }
              if (DOUBLE_NONZERO(tmp)) {
                jacobianLVDOF_addNewEntry(func_jac, var, j, tmp);
              }
              /* else { printf("there are zero terms\n"); } */
            }
            break;
          case 3:
            tmp0 = (fv->v[0] - x_dot[0]);
            tmp1 = (fv->v[1] - x_dot[1]);
            tmp2 = (fv->v[2] - x_dot[2]);
            if (TimeIntegration) {
              tmpt = (-(1. + 2. * tt) / dt) * fv->snormal[p];
            }
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              tmp = (fv->dsnormal_dx[0][p][j] * tmp0 + fv->dsnormal_dx[1][p][j] * tmp1 +
                     fv->dsnormal_dx[2][p][j] * tmp2);
              if (TimeIntegration) {
                tmp += tmpt * phi_ptr[j];
              }
              if (DOUBLE_NONZERO(tmp)) {
                jacobianLVDOF_addNewEntry(func_jac, var, j, tmp);
              }
            }
            break;
          default:
            GOMA_EH(GOMA_ERROR, "Not Covered");
          }
        }
      }
    }
    /*
     * Next, add in the velocity dependencies
     */
    if (ei[pg->imtrx]->Lvdesc_First_Var_Type[VELOCITY1] >= 0) {
      for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
        tmp = fv->snormal[kdir];
        lvdesc = ei[pg->imtrx]->Lvdesc_First_Var_Type[VELOCITY1 + kdir];
        jacobianVD_addNewEntry(func_jac, lvdesc, tmp);
      }
    }

  } /* end of if Assemble_Jacobian */
  return func_value;
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void fvelo_normal_edge_bc(double func[DIM],
                          double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                          const double vnormal,         /* normal velocity */
                          const double x_dot[MAX_PDIM], /* mesh velocity         */
                          const dbl tt,                 /* parameter to vary time integ from
                                                           explct (tt = 1) to implicit (tt = 0)  */
                          const dbl dt)                 /* current value of the time step         */

/****************************************************************************
 *
 *  Function which evaluates  n.(v - xdot) = vnormal and its sensitivities
 *  along a mesh edge.  Here n is the normal to the edge in the plane of
 *  the primary (first) surface used in identifying the edge.
 *  Note that here stangent[0][p] is the binormal to the edge as defined
 *  by the normal of the of the primary surface and the tangent vector
 *  along the edge
 *  Usage:
 *       VELO_NORMAL_EDGE 1 2 vn a
 *     In the case of free surface intercepting a plane, surface 1
 *     should be the plane surface and surface 2 the free surface.
 *     vn is the normal velocity imposed on the fluid on the
 *     line and a is an optional parameter that is not yet used.
 *     If a is not specified it will be set to -1.
 *
 *  Note that, like VELO_NORMAL, this condition replaces the normal
 *  component of the momentum equation
 *
 *
 *            Author: T.A. Baer   (10/31/94) Adapted from fvelo_normal_bc
 *
 *****************************************************************************/
{
  int j, kdir, var, p;
  double phi_j;

  if (af->Assemble_LSA_Mass_Matrix) {
    for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
      var = MESH_DISPLACEMENT1 + kdir;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          d_func[0][var][j] -= phi_j * fv->stangent[0][kdir];
        }
      }
    }
    return;
  }

  if (af->Assemble_Jacobian) {
    for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
      for (p = 0; p < pd->Num_Dim; p++) {
        var = MESH_DISPLACEMENT1 + p;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            d_func[0][var][j] += (fv->v[kdir] - x_dot[kdir]) * fv->dstangent_dx[0][kdir][p][j];
            if (TimeIntegration != 0 && p == kdir) {
              d_func[0][var][j] +=
                  (-(1. + 2. * tt) * phi_j / dt) * fv->stangent[0][kdir] * delta(p, kdir);
            }
          }
        }
      }

      var = VELOCITY1 + kdir;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          d_func[0][var][j] += phi_j * fv->stangent[0][kdir];
        }
      }
    }
  } /* end of if Assemble_Jacobian */

  /* Calculate the residual contribution	*/
  func[0] = -vnormal;
  for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
    func[0] += (fv->v[kdir] - x_dot[kdir]) * fv->stangent[0][kdir];
  }
} /* END of routine fvelo_normal_edge_bc  */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void fvelo_normal_disc_bc(double func[DIM],
                          double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                          const double vnormal,         /* normal velocity */
                          const double x_dot[MAX_PDIM], /* mesh velocity vector */
                          const dbl tt,                 /* parameter to vary time integration from
                                                         * explicit (tt = 1) to implicit (tt = 0) */
                          const dbl dt)                 /* current value of the time step         */

/************************************************************************
 *
 * fvelo_normal_disc_bc():
 *
 *  Function which evaluates the expression specifying the
 *  kinematic boundary condition at a quadrature point on a side
 *  of an element.
 *
 *         func =   - n . (rho *  vnormal) + n . (v - xdot) * rho
 *
 *  The boundary conditions, KINEMATIC_DISC and VELO_NORMAL_DISC_BC,
 *  employ this function. For internal boundaries, vnormal is set
 *  to zero,
 *  and this function is evaluated on both sides of the interface to
 *  establish mass continuity across the interface.
 *
 * Input:
 *
 *  vnormal = specified on the bc card
 *  v       = mass average velocity
 *  xdot    = velocity of the mesh (i.e., the interface if the mesh
 *            is pegged at the inteface)
 *  rho     = density of the fluid.
 *
 * Output:
 *
 *     func[0] = value of the function mentioned above
 *     d_func[0][varType][lvardof] =
 *              Derivate of func[0] wrt
 *              the variable type, varType, and the local variable
 *              degree of freedom, lvardof, corresponding to that
 *              variable type.
 *
 *   Author: P. R. Schunk    (3/24/94)
 *   Revised: RAC  (5/24/95)
 ***********************************************************************/
{
  int j, kdir, var, p;
  double phi_j;

  if (af->Assemble_LSA_Mass_Matrix) {
    for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
      var = MESH_DISPLACEMENT1 + kdir;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          d_func[0][var][j] -= mp->density * phi_j * fv->snormal[kdir];
        }
      }
    }
    return;
  }

  if (af->Assemble_Jacobian) {
    for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
      for (p = 0; p < pd->Num_Dim; p++) {
        var = MESH_DISPLACEMENT1 + p;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            d_func[0][var][j] +=
                mp->density * (fv->v[kdir] - x_dot[kdir]) * fv->dsnormal_dx[kdir][p][j];
            if (TimeIntegration != 0 && p == kdir) {
              d_func[0][var][j] +=
                  mp->density * (-(1. + 2. * tt) * phi_j / dt) * fv->snormal[kdir] * delta(p, kdir);
            }
          }
        }
      }

      var = VELOCITY1 + kdir;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          d_func[0][var][j] += mp->density * phi_j * fv->snormal[kdir];
        }
      }

      var = PVELOCITY1 + kdir;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          d_func[0][var][j] += mp->density * phi_j * fv->snormal[kdir];
        }
      }

      /*
       * HKM -> Note, we will need to add additional terms for the dependence
       *        of the density on the independent variables
       */
    }
  } /* end of if Assemble_Jacobian */

  /* Calculate the residual contribution	*/
  func[0] = -vnormal * mp->density;
  for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
    func[0] += mp->density * (fv->v[kdir] - x_dot[kdir]) * fv->snormal[kdir];
  }
} /* END of routine fvelo_normal_disc_bc  */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void fvelo_tangent_edge_bc(double func[DIM],
                           double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                           dbl *pvv,                     /* velocity of substrate */
                           const double x_dot[MAX_PDIM], /* mesh velocity */
                           const dbl tt,                 /* parameter to vary time integration from
                                                            explicit (tt = 1) to implicit (tt=0)  */
                           const dbl dt)                 /* current value of the time step        */

/******************************************************************************
 *
 *  Function which evaluates  t.(v - xdot - vweb) = 0 and its sensitivities along a mesh edge
 *  Here t is the tangent to the edge in the plane of the primary (first) surface used in
 *identifying the edge and vweb is the velocity of the a moving substrate.
 *
 *  Note that here stangent[1][p] is the tangent to the edge
 *  Usage:
 *       VELO_TANGENT_EDGE 1 2 vx vy vz
 *     In the case of free surface intercepting a plane, surface 1 should be the free surface
 *     and surface 2 the plane surface.  vx, vy, vz are the velocity components of the moving
 *     substrate.
 *
 *  This condition should replace the tangential component of the momentum equation
 *
 *
 *            Author: T.A. Baer   (10/31/94) Adapted from fvelo_normal_edge_bc
 *
 ******************************************************************************/
{
  int j, kdir, var, p;
  double phi_j, vv[DIM];

  if (af->Assemble_LSA_Mass_Matrix) {
    return;
  }
  vv[0] = (double)pvv[0];
  vv[1] = (double)pvv[1];
  vv[2] = (double)pvv[2];

  if (af->Assemble_Jacobian) {
    for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
      for (p = 0; p < pd->Num_Dim; p++) {
        var = MESH_DISPLACEMENT1 + p;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            d_func[0][var][j] += (fv->v[kdir] - vv[kdir]) * fv->dstangent_dx[1][kdir][p][j];
          }
        }
      }

      var = VELOCITY1 + kdir;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          d_func[0][var][j] += phi_j * fv->stangent[1][kdir];
        }
      }
    }
  } /* end of if Assemble_Jacobian */

  /* Calculate the residual contribution	*/
  func[0] = 0.0;
  for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
    func[0] += (fv->v[kdir] - vv[kdir]) * fv->stangent[1][kdir];
  }
} /* END of routine fvelo_tangent_edge_bc  */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void fvelo_tangent_3d(double func[MAX_PDIM],
                      double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                      const double v_s[MAX_PDIM], /* free surface velocity x_dot */
                      const double v_tangent,     /* specified tangential speed */
                      const double tx,            /* x-component of surface tangent */
                      const double ty,            /* y-component of surface tangent */
                      const double tz,            /* z-component of surface tangent */
                      const dbl tt,               /* parameter to vary time integration from *
                                                     explicit (tt = 1) to implicit (tt = 0)  */
                      const dbl dt)               /* current value of the time step   */

/******************************************************************************
 *
 *  Function which applies a tangent velocity on a surface in 3d space
 *   as     (n X t_sp) . (v - vs) = V
 *
 *   where t_sp is a specified tangent vector and n is the normal.  This
 *   works for flat surfaces and surfaces where at least component of curvature
 *   is zero.   Need to furbish for varying tangents in second direction.
 *
 *            Author: P. R. Schunk    (8/19/99)
 *
 ******************************************************************************/
{
  int j, var;
  int a, b, c, p; /* counters                   */
  double t[MAX_PDIM];
  double t_gen[MAX_PDIM], d_t_gen_dx[MAX_PDIM][MAX_PDIM][MDE];
  double phi_j;

  t[0] = tx;
  t[1] = ty;
  t[2] = tz;

  /***************************** EXECUTION BEGINS ******************************/

  t_gen[0] = 0.;
  t_gen[1] = 0.;
  t_gen[2] = 0.;
  memset(d_t_gen_dx, 0, MAX_PDIM * MAX_PDIM * MDE * sizeof(double));
  for (a = 0; a < DIM; a++) {
    for (b = 0; b < DIM; b++) {
      for (c = 0; c < DIM; c++) {
        t_gen[a] += eps(a, b, c) * fv->snormal[b] * t[c];

        for (p = 0; p < pd->Num_Dim; p++) {
          var = MESH_DISPLACEMENT1 + p;
          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_t_gen_dx[a][p][j] += eps(a, b, c) * fv->dsnormal_dx[b][p][j] * t[c];
            }
          }
        }
      }
    }
  }

  if (af->Assemble_LSA_Mass_Matrix) {
    /* Shouldn't this just be pd->Num_Dim? */
    for (a = 0; a < MAX_PDIM; a++)
      for (p = 0; p < pd->Num_Dim; p++) {
        var = MESH_DISPLACEMENT1 + p;
        if (pd->v[pg->imtrx][var])
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            d_func[0][var][j] -= phi_j * t_gen[a] * delta(a, p);
          }
      }
    return;
  }

  if (af->Assemble_Jacobian) {
    for (a = 0; a < MAX_PDIM; a++) {
      for (p = 0; p < pd->Num_Dim; p++) {
        var = MESH_DISPLACEMENT1 + p;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            d_func[0][var][j] += d_t_gen_dx[a][p][j] * (fv->v[a] - v_s[a]);
            if (TimeIntegration != 0 && p == a) {
              d_func[0][var][j] += (-(1. + 2. * tt) * phi_j / dt) * t_gen[a] * delta(a, p);
            }
          }
        }
      }

      var = VELOCITY1 + a;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          d_func[0][var][j] += phi_j * t_gen[a];
        }
      }
    }

  } /* end of if Assemble_Jacobian */

  /* Calculate the residual contribution					     */

  *func = -v_tangent;
  for (a = 0; a < DIM; a++) {
    *func += t_gen[a] * (fv->v[a] - v_s[a]);
  }

} /* END of routine fvelo_tangential_3d  */

/****************************************************************************/

void fvelo_tangential_bc(double func[],
                         double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                         double x[],             /* Solution vector for the current processor */
                         const int dcl_node,     /* globnodenum dynamiccontact line   */
                         const double vtangent,  /* specified input velocity from
                                                  * input deck                     */
                         const double beta,      /* coefficient that scales xdot slip  */
                         const double alpha,     /* coefficient that scales the       *
                                                  * exponent  of position dependent   *
                                                  *  slip velocity                    */
                         const double vtang_dot, /* surface acceleration          */
                         const double xsurf[MAX_PDIM], /* coordinates of surface  *
                                                        * Gauss point, i.e.       *
                                                        * current position        */
                         const double x_dot[MAX_PDIM], /* mesh velocity vector    */
                         const dbl tt,                 /* parameter to vary time integration from *
                                                          explicit (tt = 1) to implicit (tt = 0)  */
                         const dbl dt,                 /* current value of the time step          */
                         const int bc_type,            /* bc identifier          */
                         const int id_side,
                         double xi[DIM],
                         const Exo_DB *exo,
                         const double time_value,
                         const double u_par[],
                         const int n_par)

/******************************************************************************
 *
 *  Function which evaluates the tangential velocity boundary condition
 *        t.(v - vs - beta*xdot)=0
 *
 *            Author: P. R. Schunk    (5/17/94)
 *            Revised: RAC 5/25/95
 *
 ******************************************************************************/

{
  int j, icount, p, kdir, var;
  double phi_j;
  double dist;           /* distance btw current position and dynamic
                            contact lines                                */
  double xdcl[MAX_PDIM]; /* coordinates of dynamic contact lines         */
  dbl ead;               /* "exp( -distance_to_singularity / L_slip )" */
  double vtang_user = 0, vel_user[MAX_PDIM] = {0, 0, 0}, vel_user_dx[MAX_PDIM][MAX_PDIM];

  /***************************** EXECUTION BEGINS *******************************/

  /*This section is for position dependent slip    */
  if (alpha != 0.) {

    if (dcl_node == -1)
      GOMA_EH(GOMA_ERROR, "problem with nodeset specified in VELO_TANGENT BC");

    /* calculate position of dynamic contact line from dcl_node               */
    for (icount = 0; icount < pd->Num_Dim; icount++) {
      /*
       * Do the distance calculation based on undeformed geometry ...
       * this save us Jacobian entries in uncharted sections of the
       * A matrix
       */
      xdcl[icount] = (Coor[icount][dcl_node]);
    }
    dist = 0.;
    /* calculate distance from dcl to current Gauss point */
    for (icount = 0; icount < pd->Num_Dim; icount++) {
      dist += (xsurf[icount] - xdcl[icount]) * (xsurf[icount] - xdcl[icount]);
    }
    dist = sqrt(dist);
  }
  /* turn this off if we don't have position depend slip
     alpha will be zero, so dist and xdcl can be anything as
     long as they are defined                                */

  else {
    dist = 1.;
    xdcl[0] = xdcl[1] = xdcl[2] = 0.;
  }

  /*
   * Protect exponentially decaying slip from underflow far away from the
   * singularity.
   */
  if (alpha * dist > 40) {
    ead = (dbl)0;
  } else {
    ead = exp(-alpha * dist);
  }

  if (af->Assemble_LSA_Mass_Matrix) {
    for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
      var = MESH_DISPLACEMENT1 + kdir;
      if (pd->v[pg->imtrx][var])
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          d_func[0][var][j] -= phi_j * beta * ead * fv->stangent[0][kdir];
        }
    }
    return;
  }

  /* sum the contributions to the global stiffness matrix */
  if (af->Assemble_Jacobian) {
    for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
      for (p = 0; p < pd->Num_Dim; p++) {
        var = MESH_DISPLACEMENT1 + p;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            d_func[0][var][j] +=
                (fv->v[kdir] - x_dot[kdir] * beta * ead) * fv->dstangent_dx[0][kdir][p][j];

            if (TimeIntegration != 0 && p == kdir) {
              d_func[0][var][j] +=
                  (-(1. + 2. * tt) * phi_j / dt * beta * ead) * fv->stangent[0][kdir];
            }
          }
          /*
           * add in Jacobian contributions for dynamic contact line node
           */

          /* find local dof number of DCL
             NOTE: if DCL is not in this element, you'll need to something
             special here - see RRR or PAS */
          /*	  debug_var = 0.;
                  if( (TimeIntegration != 0)	&&
                  (beta != 0.)			&&
                  (alpha != 0.)			&&
                  (debug_var != 0.) )
                  {

                  jdcl = -1;
                  for  ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
                  {
                  if (ei[pg->imtrx]->gnn_list[var][j] == dcl_node) jdcl = j;
                  }
                  GOMA_EH(jdcl, "Dynamic contact line node not found in element");

                  d_func[0][p][jdcl] -= beta * x_dot[kdir]
                  * exp(-dist * alpha) * (xsurf[p]- xdcl[p]) / dist * alpha;
                  } */
        }
      }

      var = VELOCITY1 + kdir;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          d_func[0][var][j] += phi_j * fv->stangent[0][kdir];
        }
      }
      var = PVELOCITY1 + kdir;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          d_func[0][var][j] += phi_j * fv->stangent[0][kdir];
        }
      }
    }
  } /* end of if Assemble_Jacobian */

  /* Calculate the residual contribution	*/
  if (pd->Num_Dim == 2) {
    /*  Velo-tangent bc	*/
    if (bc_type != VELO_TANGENT_USER_BC) {
      *func -= vtangent * (1.0 + vtang_dot * time_value);
    } else {
#ifdef FEATURE_ROLLON_PLEASE
#include "feature_rollon_velo.h"
#else
      vel_user[0] = velo_vary_fnc(UVARY_BC, fv->x[0], fv->x[1], fv->x[2], u_par, time_value);
      vel_user[1] = velo_vary_fnc(VVARY_BC, fv->x[0], fv->x[1], fv->x[2], u_par, time_value);
      if (af->Assemble_Jacobian) {
        if (pd->e[pg->imtrx][R_MESH1]) {
          vel_user_dx[0][0] =
              dvelo_vary_fnc_d1(UVARY_BC, fv->x[0], fv->x[1], fv->x[2], u_par, time_value);
          vel_user_dx[0][1] =
              dvelo_vary_fnc_d2(UVARY_BC, fv->x[0], fv->x[1], fv->x[2], u_par, time_value);
          vel_user_dx[1][0] =
              dvelo_vary_fnc_d1(VVARY_BC, fv->x[0], fv->x[1], fv->x[2], u_par, time_value);
          vel_user_dx[1][1] =
              dvelo_vary_fnc_d2(VVARY_BC, fv->x[0], fv->x[1], fv->x[2], u_par, time_value);
        }
      }
#endif
      vtang_user = fv->stangent[0][0] * vel_user[0] + fv->stangent[0][1] * vel_user[1];
      *func -= vtang_user;
      if (af->Assemble_Jacobian) {
        for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
          for (p = 0; p < pd->Num_Dim; p++) {
            var = MESH_DISPLACEMENT1 + p;
            if (pd->v[pg->imtrx][var]) {
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                phi_j = bf[var]->phi[j];
                d_func[0][var][j] +=
                    (fv->v[kdir] - x_dot[kdir] * beta * ead) * fv->dstangent_dx[0][kdir][p][j];
                d_func[0][var][j] -= (vel_user[kdir] * fv->dstangent_dx[0][kdir][p][j] +
                                      fv->stangent[0][kdir] * vel_user_dx[kdir][p] * phi_j);
              }
            }
          }
        }
      }
    }
    /*  Acoustic Velo-streaming bc	*/
    if (bc_type == VELO_STREAMING_BC) {
      double velo_stream, d_phi_dxi[MDE];
      double dx_dxi[DIM], sign;
      int el1 = ei[pg->imtrx]->ielem;
      int el2, nf, i;
      int *n_dof = NULL;
      int node, index, dof_map[MDE];
      int n_dofptr[MAX_VARIABLE_TYPES][MDE];
      double *n_esp, shell_var[MDE];
      /* See if there is a friend for this element */
      nf = num_elem_friends[el1];
      if (nf == 0)
        GOMA_EH(GOMA_ERROR, "no element friends");
      el2 = elem_friends[el1][0];
      if (nf != 1)
        GOMA_WH(-1, "WARNING: Not set up for more than one element friend!");
      n_dof = (int *)array_alloc(1, MAX_VARIABLE_TYPES, sizeof(int));
      load_neighbor_var_data(el1, el2, n_dof, dof_map, n_dofptr, id_side, xi, exo);
      for (kdir = 0; kdir < n_dof[SHELL_BDYVELO]; kdir++) {
        d_phi_dxi[kdir] = bf[SHELL_BDYVELO]->dphidxi[kdir][0];
        n_esp = x_static + n_dofptr[SHELL_BDYVELO][kdir];
        shell_var[kdir] = *n_esp;
      }
      dx_dxi[0] = dx_dxi[1] = 0.;
      if (ei[pg->imtrx]->deforming_mesh) {
        for (i = 0; i < n_dof[SHELL_BDYVELO]; i++) {
          kdir = dof_map[i];
          node = ei[pg->imtrx]->dof_list[R_MESH1][kdir];
          index = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + node];
          dx_dxi[0] += (Coor[0][index] + *esp->d[0][kdir]) * d_phi_dxi[i];
          dx_dxi[1] += (Coor[1][index] + *esp->d[1][kdir]) * d_phi_dxi[i];
        }
      } else {
        for (i = 0; i < n_dof[SHELL_BDYVELO]; i++) {
          kdir = dof_map[i];
          node = ei[pg->imtrx]->dof_list[R_MESH1][kdir];
          index = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + node];
          dx_dxi[0] += Coor[0][index] * d_phi_dxi[i];
          dx_dxi[1] += Coor[1][index] * d_phi_dxi[i];
        }
      }
      sign = (dx_dxi[0] * fv->stangent[0][0] + dx_dxi[1] * fv->stangent[0][1]) > 0 ? 1 : -1;

      velo_stream = 0;
      for (kdir = 0; kdir < n_dof[SHELL_BDYVELO]; kdir++) {
        velo_stream += shell_var[kdir] * bf[SHELL_BDYVELO]->dphidxi[kdir][0] * fv->h3 / fv->sdet;
      }
      *func += velo_stream * sign * 3. / (8. * upd->Acoustic_Frequency);
      if (af->Assemble_Jacobian) {
        var = SHELL_BDYVELO;
        if (pd->v[pg->imtrx][var]) {
          for (kdir = 0; kdir < n_dof[var]; kdir++) {
            /*  Find the right variable in the bulk context	*/
            for (el1 = 0; el1 < n_dof[var]; el1++) {
              if (dof_map[kdir] == ei[pg->imtrx]->dof_list[var][el1])
                el2 = el1;
            }
            d_func[0][var][el2] +=
                bf[var]->dphidxi[kdir][0] * sign * 3. / (8. * upd->Acoustic_Frequency);
          }
        }
      }
    }
    for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
      *func += (fv->v[kdir] - x_dot[kdir] * beta * ead) * fv->stangent[0][kdir];
    }
  } else if (pd->Num_Dim == 3) {
    GOMA_EH(GOMA_ERROR, "Velo-tangent in 3D is ill-defined");
  }
}
/* END of routine fvelo_tangential_bc  */

/****************************************************************************/

/****************************************************************************/
/*************************************************************************
 *
 *  Function which evaluates the slip velocity boundary condition due to
 *  charged species transport along a charged EDL.  This is also known as
 *  electroosmotic pumping.  In its current form, the slip velocity is only
 *  tied to the potential gradient tangent to the charged wall.  The equation
 *  is also known as the Helmholtz-Smoluchowski equation. In capillary flow,
 *  the flow velocity along the wall is a function of only axial potential
 *  gradient, zeta potential, permittivity, and fluid viscosity.
 *
 *     v(p) = -eta*zeta*grad_V(p)/mu
 *
 *            Author: A.C. Sun    (02/14/03)
 ************************************************************************/

void fvelo_slip_electrokinetic_bc(double func[MAX_PDIM],
                                  double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                                  double permittivity, /* permittivity of the transporting medium */
                                  double zeta_potential) /* zeta_potential at the wall */
{

  /* Local variables */

  int a, j, var, jvar, dim, kdir;
  double phi_j, gradV[MAX_PDIM], vs;
  double eta, zetaV; /* permittivity and zeta potential */

  eta = permittivity;
  zetaV = zeta_potential;

  dim = pd->Num_Dim;
  vs = 0.;
  for (a = 0; a < dim; a++) {
    gradV[a] = fv->grad_V[a];
    vs -= gradV[a] * fv->stangent[0][a];
  }
  vs *= eta * zetaV / mp->viscosity;

  if (af->Assemble_Jacobian) {
    for (jvar = 0; jvar < pd->Num_Dim; jvar++) {
      var = VELOCITY1 + jvar;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          d_func[0][var][j] -= phi_j * fv->stangent[0][jvar];
        }
      }
    }
  }

  /* Calculate the residual contribution	*/
  if (pd->Num_Dim == 2) {
    *func += vs;
    for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
      *func -= (fv->v[kdir]) * fv->stangent[0][kdir];
    }
  } else if (pd->Num_Dim == 3) {
    GOMA_EH(GOMA_ERROR, "Velo-electrokinetic in 3D is ill-defined");
  }

} /* END of routine fvelo_slip_electrokinetic_bc  */
/******************************************************************************/
/****************************************************************************/
/****************************************************************************/

void fvelo_electrokinetic_3d(
    double func[MAX_PDIM],
    double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
    const double v_s[MAX_PDIM],  /* free surface velocity x_dot */
    const double permittivity,   /* permittivity of the transporting medium */
    const double zeta_potential, /* zeta_potential at the wall */
    const double tx,             /* x-component of surface tangent */
    const double ty,             /* y-component of surface tangent */
    const double tz,             /* z-component of surface tangent */
    const dbl tt,                /* parameter to vary time integration from *
                                    explicit (tt = 1) to implicit (tt = 0)  */
    const dbl dt)                /* current value of the time step   */

/******************************************************************************
 *
 *  Function which applies a tangent velocity on a surface in 3d space
 *   as     (n X t_sp) . (v - vs) = V
 *
 *   where t_sp is a specified tangent vector and n is the normal.  This
 *   works for flat surfaces and surfaces where at least component of curvature
 *   is zero.   Need to furbish for varying tangents in second direction.
 *
 *            Author: P. R. Schunk    (8/19/99)
 *                    A.C. Sun (5/8/03)
 ******************************************************************************/
{
  int j, var;
  int a, b, c, p; /* counters                   */
  double t[MAX_PDIM], gradV[MAX_PDIM], v_tangent[MAX_PDIM];
  double t_gen[MAX_PDIM], d_t_gen_dx[MAX_PDIM][MAX_PDIM][MDE];
  double phi_j;
  double eta, zetaV; /* permittivity and zeta potential at the wall*/

  t[0] = tx;
  t[1] = ty;
  t[2] = tz;

  /***************************** EXECUTION BEGINS ******************************/

  t_gen[0] = 0.;
  t_gen[1] = 0.;
  t_gen[2] = 0.;
  v_tangent[0] = 0.;
  v_tangent[1] = 0.;
  v_tangent[2] = 0.;
  memset(d_t_gen_dx, 0, MAX_PDIM * MAX_PDIM * MDE * sizeof(double));

  eta = permittivity;
  zetaV = zeta_potential;

  for (a = 0; a < DIM; a++) {
    for (b = 0; b < DIM; b++) {
      for (c = 0; c < DIM; c++) {
        t_gen[a] += eps(a, b, c) * fv->snormal[b] * t[c];

        for (p = 0; p < pd->Num_Dim; p++) {
          var = MESH_DISPLACEMENT1 + p;
          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_t_gen_dx[a][p][j] += eps(a, b, c) * fv->dsnormal_dx[b][p][j] * t[c];
            }
          }
        }
      }
    }
  }

  if (af->Assemble_LSA_Mass_Matrix) {
    /* Shouldn't this just be pd->Num_Dim? */
    for (a = 0; a < MAX_PDIM; a++)
      for (p = 0; p < pd->Num_Dim; p++) {
        var = MESH_DISPLACEMENT1 + p;
        if (pd->v[pg->imtrx][var])
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            d_func[0][var][j] -= phi_j * t_gen[a] * delta(a, p);
          }
      }
    return;
  }

  if (af->Assemble_Jacobian) {
    for (a = 0; a < MAX_PDIM; a++) {
      for (p = 0; p < pd->Num_Dim; p++) {
        var = MESH_DISPLACEMENT1 + p;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            d_func[0][var][j] += d_t_gen_dx[a][p][j] * (fv->v[a] - v_s[a]);
            if (TimeIntegration != 0 && p == a) {
              d_func[0][var][j] += (-(1. + 2. * tt) * phi_j / dt) * t_gen[a] * delta(a, p);
            }
          }
        }
      }

      var = VELOCITY1 + a;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          d_func[0][var][j] += phi_j * t_gen[a];
        }
      }
    }

  } /* end of if Assemble_Jacobian */

  /* Calculate the residual contribution					     */

  *func = 0.;
  for (a = 0; a < DIM; a++) {
    gradV[a] = fv->grad_V[a];
    v_tangent[a] = gradV[a] * eta * zetaV / mp->viscosity;
    *func += t_gen[a] * (v_tangent[a] + fv->v[a] - v_s[a]);
  }

} /* END of routine fvelo_electrokinetic_3d  */
/****************************************************************************/
/****************************************************************************/
/******************************************************************************
 *
 * fvelo_tangential_solid_bc() - match tangential velocity at fluid/solid bnd
 *
 * Function which evaluates the tangential velocity boundary condition at a
 * fluid - Lagrangian solid boundary.
 *
 *        t.(v - v_sfs.F)=0   for VELO_TANGENT_SOLID   (strong form)
 *        t.(v - v_sfs.F)/beta = t.n.T    for VELO_SLIP_SOLID (weak form)
 *
 * Note, beta can be fixed, for a pre-defined slip region app, or varible
 * to account for position-independent slip results.
 *
 *            Author: P. R. Schunk    (11/4/98))
 *            Revised: RAC 5/25/95
 *
 *****************************************************************************/

void fvelo_tangential_solid_bc(double func[],
                               double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                               double x[],
                               double x_dot[],    /* mesh velocity vector                      */
                               double x_rs_dot[], /* real solid velocity vector */
                               const int type,
                               const int eb_mat_solid,
                               const int eb_mat_fluid, /* element block ID's for solid and fluid */
                               dbl beta,
                               const int dcl_node,           /*   node id for DCL  */
                               double alpha,                 /* extent of slip  */
                               const double xsurf[MAX_PDIM], /* coordinates of surface Gauss  *
                                                              * point, i.e. current position  */
                               const dbl tt, /* parameter to vary time integration method */
                               const dbl dt) /* current value of the time step */
{

  int j, w, j_id;
  int a, p, kdir, jvar, var = -1;
  double phi_j;
  int dim, err;
  double betainv;
  int icount;
  double dist; /* distance btw current position and dynamic
                     contact lines                                */
  double v_solid_mesh[DIM];
  double xdcl[MAX_PDIM]; /* coordinates of dynamic contact lines         */
  double disp;           /* DCL point displacement */

  double vconv[MAX_PDIM];     /*Calculated convection velocity */
  double vconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;

  /***************************** EXECUTION BEGINS *******************************/
  /* Calculate the residual contribution	*/

  if (type == VELO_TANGENT_SOLID_BC) {
    beta = 1.0;
    alpha = 0;
    dist = 1;
  } else if (type != VELO_SLIP_SOLID_BC) {
    GOMA_EH(GOMA_ERROR, "BC must be VELO_TANGENT_SOLID or VELO_SLIP_SOLID");
  }
  if (alpha == 0.) {
    betainv = 1. / beta;
    dist = 1;
  } else if (alpha != 0.) {

    for (icount = 0; icount < pd->Num_Dim; icount++) {
      disp = x[Index_Solution(dcl_node, MESH_DISPLACEMENT1 + icount, 0, 0, -1, pg->imtrx)];
      xdcl[icount] = (Coor[icount][dcl_node] + disp);
    }

    dist = 0.;
    for (icount = 0; icount < pd->Num_Dim; icount++) {
      /*UNDEFORMED* dist +=
       * (xsurf[icount]-Coor[icount][dcl_node])*(xsurf[icount]-Coor[icount][dcl_node]); */
      /*DEFORMED* dist += (fv->x[icount]-xdcl[icount])*(fv->x[icount]-xdcl[icount]); */
      dist += (xsurf[icount] - Coor[icount][dcl_node]) * (xsurf[icount] - Coor[icount][dcl_node]);
    }
    dist = sqrt(dist);
  }

  else {
    dist = 1.;
    xdcl[0] = xdcl[1] = xdcl[2] = 0.;
  }

  /* for exponentially decaying slip, max out betainv when equivalent to
   * STRONG IC for no slip
   */
  if (alpha * dist > log(beta * BIG_PENALTY)) {
    betainv = BIG_PENALTY;
  } else {
    betainv = exp(alpha * dist) / beta;
  }

  dim = pd->Num_Dim;
  if (dim != 2 && type == VELO_TANGENT_SOLID_BC)
    GOMA_EH(GOMA_ERROR, "VELO_TANGENT_SOLID not implemented in 3D yet");

  if (af->Assemble_LSA_Mass_Matrix) {
    if (Current_EB_ptr->Elem_Blk_Id == eb_mat_solid)
      if (mp->PorousMediaType == POROUS_SATURATED || mp->PorousMediaType == CONTINUOUS)
        for (jvar = 0; jvar < dim; jvar++) {
          if (pd->MeshMotion == LAGRANGIAN || pd->MeshMotion == DYNAMIC_LAGRANGIAN)
            var = MESH_DISPLACEMENT1 + jvar;
          else if (pd->MeshMotion == TOTAL_ALE)
            var = SOLID_DISPLACEMENT1 + jvar;
          else
            GOMA_EH(GOMA_ERROR, "Bad pd->MeshMotion");
          if (pd->v[pg->imtrx][var])
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];
              d_func[0][var][j] += phi_j * fv->stangent[0][jvar] * betainv;
            }
        }
    return;
  }

  if (Current_EB_ptr->Elem_Blk_Id == eb_mat_fluid) {
    if (type == VELO_TANGENT_SOLID_BC) {
      for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
        { *func += fv->v[kdir] * fv->stangent[0][kdir] * betainv; }
      }
    } else {
      for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
        func[kdir] -= fv->v[kdir] * betainv;
      }
    }
    /* sum the contributions to the global stiffness matrix */
    if (af->Assemble_Jacobian) {
      if (type == VELO_TANGENT_SOLID_BC) {
        for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
          for (p = 0; p < pd->Num_Dim; p++) {
            var = MESH_DISPLACEMENT1 + p;
            if (pd->v[pg->imtrx][var]) {
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                phi_j = bf[var]->phi[j];
                d_func[0][var][j] += (fv->v[kdir]) * fv->dstangent_dx[0][kdir][p][j] * betainv;
              }
            }
          }

          var = VELOCITY1 + kdir;
          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];
              d_func[0][var][j] += phi_j * fv->stangent[0][kdir] * betainv;
            }
          }
        }
      } else {
        for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
          var = VELOCITY1 + kdir;
          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];
              d_func[kdir][var][j] -= phi_j * betainv;
            }
          }
        }
      }
    } /* end of if Assemble_Jacobian */
  } else if (Current_EB_ptr->Elem_Blk_Id == eb_mat_solid) {
    /*first load up the correct reference velocity and do this in the
     solid phase because you may not have it in the liquid phase, as in
     TALE */
    if (pd->MeshMotion == LAGRANGIAN || pd->MeshMotion == DYNAMIC_LAGRANGIAN) {
      for (a = 0; a < DIM; a++)
        v_solid_mesh[a] = x_dot[a];
      err = belly_flop(elc->lame_mu);
      GOMA_EH(err, "error in belly flop");
      err = get_convection_velocity(vconv, vconv_old, d_vconv, dt, tt);

    } else if (pd->MeshMotion == TOTAL_ALE) {
      for (a = 0; a < DIM; a++)
        v_solid_mesh[a] = x_rs_dot[a];
      err = belly_flop_rs(elc_rs->lame_mu);
      GOMA_EH(err, "error in belly flop");
      err = get_convection_velocity_rs(vconv, vconv_old, d_vconv, dt, tt);

    } else {
      GOMA_EH(GOMA_ERROR,
              "Shouldn't be in this section of velo_tangent_solid with an arbitrary solid");
    }

    if (mp->PorousMediaType == POROUS_SATURATED || mp->PorousMediaType == CONTINUOUS ||
        mp->PorousMediaType == POROUS_UNSATURATED) {
      if (type == VELO_TANGENT_SOLID_BC) {
        for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
          *func += (vconv[kdir] + v_solid_mesh[kdir]) * fv->stangent[0][kdir] * betainv;
        }
      } else {
        for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
          func[kdir] += (vconv[kdir] + v_solid_mesh[kdir]) * betainv;
        }
      }
      if (af->Assemble_Jacobian) {

        /* sum the contributions to the global stiffness matrix */

        if (type == VELO_TANGENT_SOLID_BC) {
          var = TEMPERATURE;
          for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
            if (pd->v[pg->imtrx][var]) {
              phi_j = bf[var]->phi[j_id];
              for (kdir = 0; kdir < dim; kdir++) {
                /*     d( )/dT        */
                d_func[0][var][j_id] += d_vconv->T[kdir][j_id] * fv->stangent[0][kdir] * betainv;
              }
            }
          }

          var = MASS_FRACTION;
          for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
            if (pd->v[pg->imtrx][var]) {
              phi_j = bf[var]->phi[j_id];
              for (kdir = 0; kdir < dim; kdir++) {
                for (w = 0; w < pd->Num_Species_Eqn; w++) {
                  /*     d( )/dC        */
                  d_func[0][MAX_VARIABLE_TYPES + w][j_id] +=
                      d_vconv->C[kdir][w][j_id] * fv->stangent[0][kdir] * betainv;
                }
              }
            }
          }

          /* for real and pseudo-solid displ. there are too different
             pieces.  First the non-time-dependent pieces */

          for (jvar = 0; jvar < dim; jvar++) {
            var = MESH_DISPLACEMENT1 + jvar;
            for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
              if (pd->v[pg->imtrx][var]) {
                phi_j = bf[var]->phi[j_id];
                for (kdir = 0; kdir < dim; kdir++) {
                  /*     d( )/dx        */
                  d_func[0][var][j_id] +=
                      ((vconv[kdir] + v_solid_mesh[kdir]) * fv->dstangent_dx[0][kdir][jvar][j_id] +
                       d_vconv->X[kdir][jvar][j_id] * fv->stangent[0][kdir]) *
                      betainv;
                }
              }
            }
          }
          for (jvar = 0; jvar < dim; jvar++) {
            var = SOLID_DISPLACEMENT1 + jvar;
            for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
              if (pd->v[pg->imtrx][var]) {
                phi_j = bf[var]->phi[j_id];
                for (kdir = 0; kdir < dim; kdir++) {
                  /*     d( )/dx        */
                  d_func[0][var][j_id] +=
                      d_vconv->rs[kdir][jvar][j_id] * fv->stangent[0][kdir] * betainv;
                }
              }
            }
          }
          /* Now the x_dot and x_rs_dot term */

          for (jvar = 0; jvar < dim; jvar++) {
            if (pd->MeshMotion == LAGRANGIAN || pd->MeshMotion == DYNAMIC_LAGRANGIAN) {
              var = MESH_DISPLACEMENT1 + jvar;
            } else if (pd->MeshMotion == TOTAL_ALE) {
              var = SOLID_DISPLACEMENT1 + jvar;
            }
            if (pd->v[pg->imtrx][var]) {
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                phi_j = bf[var]->phi[j];
                if (TimeIntegration != 0) {
                  d_func[0][var][j] +=
                      ((1. + 2. * tt) * phi_j / dt) * fv->stangent[0][jvar] * betainv;
                }
              }
            }
          }

          for (jvar = 0; jvar < dim; jvar++) {
            var = VELOCITY1 + jvar;
            for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
              if (pd->v[pg->imtrx][var]) {
                phi_j = bf[var]->phi[j_id];
                for (kdir = 0; kdir < dim; kdir++) {
                  d_func[0][var][j_id] +=
                      d_vconv->v[kdir][jvar][j_id] * fv->stangent[0][kdir] * betainv;
                }
              }
            }
          }
        } else {
          var = TEMPERATURE;
          for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
            if (pd->v[pg->imtrx][var]) {
              phi_j = bf[var]->phi[j_id];
              for (kdir = 0; kdir < dim; kdir++) {
                /*     d( )/dT        */
                d_func[kdir][var][j_id] += d_vconv->T[kdir][j_id] * betainv;
              }
            }
          }

          var = MASS_FRACTION;
          for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
            if (pd->v[pg->imtrx][var]) {
              phi_j = bf[var]->phi[j_id];
              for (kdir = 0; kdir < dim; kdir++) {
                for (w = 0; w < pd->Num_Species_Eqn; w++) {
                  /*     d( )/dC        */
                  d_func[kdir][MAX_VARIABLE_TYPES + w][j_id] += d_vconv->C[kdir][w][j_id] * betainv;
                }
              }
            }
          }
          /* for real and pseudo-solid displ. there are too different
             pieces.  First the non-time-dependent pieces */

          for (jvar = 0; jvar < dim; jvar++) {
            var = MESH_DISPLACEMENT1 + jvar;
            for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
              if (pd->v[pg->imtrx][var]) {
                phi_j = bf[var]->phi[j_id];
                for (kdir = 0; kdir < dim; kdir++) {
                  /*     d( )/dx        */
                  d_func[kdir][var][j_id] += d_vconv->X[kdir][jvar][j_id] * betainv;
                }
              }
            }
          }
          for (jvar = 0; jvar < dim; jvar++) {
            var = SOLID_DISPLACEMENT1 + jvar;
            for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
              if (pd->v[pg->imtrx][var]) {
                phi_j = bf[var]->phi[j_id];
                for (kdir = 0; kdir < dim; kdir++) {
                  /*     d( )/dx        */
                  d_func[kdir][var][j_id] += d_vconv->rs[kdir][jvar][j_id] * betainv;
                }
              }
            }
          }

          /* Now the x_dot and x_rs_dot term */

          for (jvar = 0; jvar < dim; jvar++) {
            if (pd->MeshMotion == LAGRANGIAN || pd->MeshMotion == DYNAMIC_LAGRANGIAN) {
              var = MESH_DISPLACEMENT1 + jvar;
            } else if (pd->MeshMotion == TOTAL_ALE) {
              var = SOLID_DISPLACEMENT1 + jvar;
            }
            if (pd->v[pg->imtrx][var]) {
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                phi_j = bf[var]->phi[j];
                if (TimeIntegration != 0) {
                  d_func[jvar][var][j] += ((1. + 2. * tt) * phi_j / dt) * betainv;
                }
              }
            }
          }

          for (jvar = 0; jvar < dim; jvar++) {
            var = VELOCITY1 + jvar;
            for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
              if (pd->v[pg->imtrx][var]) {
                phi_j = bf[var]->phi[j_id];
                for (kdir = 0; kdir < dim; kdir++) {
                  d_func[kdir][var][j_id] += d_vconv->v[kdir][jvar][j_id] * betainv;
                }
              }
            }
          }
        }
      } /*  if Assemble_Jacobian  */
    }
    /*      else if (mp->PorousMediaType == POROUS_UNSATURATED)
            {
              GOMA_EH(GOMA_ERROR,"VELO_TANGENT_SOLID not for  POROUS_UNSATURATED MEDIA. USE
       NO_SLIP");
            }  */
    else
      GOMA_EH(GOMA_ERROR, "bad media type in VELO_TANGENT_SOLID");
  }
}
/* END of routine fvelo_tangential_solid_bc  */
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

/******************************************************************************
 *
 * fvelo_normal_solid_bc() - match normal velocity at fluid/solid bnd
 *
 * Function which evaluates the normal velocity boundary condition at a
 * fluid - Lagrangian solid boundary.
 *
 *        n.(v - v_sfs.F)=0   for VELO_NORMAL_SOLID   (strong form)
 *
 *
 *            Author: P. R. Schunk    (11/4/98))
 *            Revised: TAB (11/2001) cloned from velo_tangential_solid_bc
 *
 *****************************************************************************/
void fvelo_normal_solid_bc(double func[],
                           double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                           double x[MAX_PDIM],
                           double x_dot[DIM],    /* mesh velocity vector                      */
                           double x_rs_dot[DIM], /* real solid velocity vector */
                           const int type,       /* BC_Name */
                           const int eb_mat_solid,
                           const int eb_mat_fluid, /* element block ID's for solid and fluid */
                           const dbl tt,           /* parameter to vary time integration method */
                           const dbl dt)           /* current value of the time step */
{

  int j, w, j_id;
  int a, p, kdir, jvar, var = -1;
  double phi_j;
  int dim, err;

  double v_solid_mesh[DIM];

  double vconv[MAX_PDIM];     /*Calculated convection velocity */
  double vconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;
  /* Calculate the residual contribution	*/

  dim = pd->Num_Dim;
  if (dim != 2)
    GOMA_EH(GOMA_ERROR, "VELO_NORMAL_SOLID not implemented in 3D yet");

  if (af->Assemble_LSA_Mass_Matrix) {
    if (Current_EB_ptr->Elem_Blk_Id == eb_mat_solid)
      if (mp->PorousMediaType == POROUS_SATURATED || mp->PorousMediaType == CONTINUOUS)
        for (jvar = 0; jvar < dim; jvar++) {
          if (pd->MeshMotion == LAGRANGIAN || pd->MeshMotion == DYNAMIC_LAGRANGIAN)
            var = MESH_DISPLACEMENT1 + jvar;
          else if (pd->MeshMotion == TOTAL_ALE)
            var = SOLID_DISPLACEMENT1 + jvar;
          else
            GOMA_EH(GOMA_ERROR, "Bad pd->MeshMotion");
          if (pd->v[pg->imtrx][var])
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];
              d_func[0][var][j] += phi_j * fv->snormal[jvar];
            }
        }
    return;
  }

  if (Current_EB_ptr->Elem_Blk_Id == eb_mat_fluid) {
    for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
      *func += fv->v[kdir] * fv->snormal[kdir];
    }
    /* sum the contributions to the global stiffness matrix */
    if (af->Assemble_Jacobian) {
      for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
        for (p = 0; p < pd->Num_Dim; p++) {
          var = MESH_DISPLACEMENT1 + p;
          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];
              d_func[0][var][j] += (fv->v[kdir]) * fv->dsnormal_dx[kdir][p][j];
            }
          }
        }

        var = VELOCITY1 + kdir;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            d_func[0][var][j] += phi_j * fv->snormal[kdir];
          }
        }
      }
    } /* end of if Assemble_Jacobian */
  } else if (Current_EB_ptr->Elem_Blk_Id == eb_mat_solid) {
    /*first load up the correct reference velocity and do this in the
     solid phase because you may not have it in the liquid phase, as in
     TALE */
    if (pd->MeshMotion == LAGRANGIAN || pd->MeshMotion == DYNAMIC_LAGRANGIAN) {
      for (a = 0; a < DIM; a++)
        v_solid_mesh[a] = x_dot[a];
      err = belly_flop(elc->lame_mu);
      GOMA_EH(err, "error in belly flop");
      err = get_convection_velocity(vconv, vconv_old, d_vconv, dt, tt);

    } else if (pd->MeshMotion == TOTAL_ALE) {
      for (a = 0; a < DIM; a++)
        v_solid_mesh[a] = x_rs_dot[a];
      err = belly_flop_rs(elc_rs->lame_mu);
      GOMA_EH(err, "error in belly flop");
      err = get_convection_velocity_rs(vconv, vconv_old, d_vconv, dt, tt);

    } else {
      GOMA_EH(GOMA_ERROR,
              "Shouldn't be in this section of velo_tangent_solid with an arbitrary solid");
    }

    if (mp->PorousMediaType == POROUS_SATURATED || mp->PorousMediaType == CONTINUOUS) {
      for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
        *func += (vconv[kdir] + v_solid_mesh[kdir]) * fv->snormal[kdir];
      }

      if (af->Assemble_Jacobian) {

        /* sum the contributions to the global stiffness matrix */

        var = TEMPERATURE;
        for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
          if (pd->v[pg->imtrx][var]) {
            phi_j = bf[var]->phi[j_id];
            for (kdir = 0; kdir < dim; kdir++) {
              /*     d( )/dT        */
              d_func[0][var][j_id] += d_vconv->T[kdir][j_id] * fv->snormal[kdir];
            }
          }
        }

        var = MASS_FRACTION;
        for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
          if (pd->v[pg->imtrx][var]) {
            phi_j = bf[var]->phi[j_id];
            for (kdir = 0; kdir < dim; kdir++) {
              for (w = 0; w < pd->Num_Species_Eqn; w++) {
                /*     d( )/dC        */
                d_func[0][MAX_VARIABLE_TYPES + w][j_id] +=
                    d_vconv->C[kdir][w][j_id] * fv->snormal[kdir];
              }
            }
          }
        }

        /* for real and pseudo-solid displ. there are too different
           pieces.  First the non-time-dependent pieces */

        for (jvar = 0; jvar < dim; jvar++) {
          var = MESH_DISPLACEMENT1 + jvar;
          for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
            if (pd->v[pg->imtrx][var]) {
              phi_j = bf[var]->phi[j_id];
              for (kdir = 0; kdir < dim; kdir++) {
                /*     d( )/dx        */
                d_func[0][var][j_id] +=
                    (vconv[kdir] + v_solid_mesh[kdir]) * fv->dsnormal_dx[kdir][jvar][j_id] +
                    d_vconv->X[kdir][jvar][j_id] * fv->snormal[kdir];
              }
            }
          }
        }

        /* Now the x_dot and x_rs_dot term */

        for (jvar = 0; jvar < dim; jvar++) {
          if (pd->MeshMotion == LAGRANGIAN || pd->MeshMotion == DYNAMIC_LAGRANGIAN) {
            var = MESH_DISPLACEMENT1 + jvar;
          } else if (pd->MeshMotion == TOTAL_ALE) {
            var = SOLID_DISPLACEMENT1 + jvar;
          }
          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];
              if (TimeIntegration != 0) {
                d_func[0][var][j] += ((1. + 2. * tt) * phi_j / dt) * fv->snormal[jvar];
              }
            }
          }
        }

        for (jvar = 0; jvar < dim; jvar++) {
          var = VELOCITY1 + jvar;
          for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
            if (pd->v[pg->imtrx][var]) {
              phi_j = bf[var]->phi[j_id];
              d_func[0][var][j_id] += d_vconv->v[kdir][jvar][j_id] * fv->snormal[jvar];
            }
          }
        }
      }
    } else if (mp->PorousMediaType == POROUS_UNSATURATED) {
      GOMA_EH(GOMA_ERROR, "VELO_TANGENT_SOLID not for  POROUS_UNSATURATED MEDIA. USE NO_SLIP");
    } else
      GOMA_EH(GOMA_ERROR, "bad media type in VELO_TANGENT_SOLID");
  }
}
/* END of routine fvelo_normal_solid_bc  */

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

void fvelo_slip_bc(double func[MAX_PDIM],
                   double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                   double x[],
                   const int type,      /* whether rotational or not */
                   const int max_float, /* Max float number from input file */
                   double bc_float[MAX_BC_FLOAT_DATA],
                   const int dcl_node,           /*   node id for DCL  */
                   const double xsurf[MAX_PDIM], /* coordinates of surface Gauss  *
                                                  * point, i.e. current position  */
                   const double tt,              /* parameter in time stepping alg           */
                   const double dt)              /* current time step value                  */

/*************************************************************************
 *
 *  Function which evaluates the slip velocity boundary condition as an
 *  integraged boundary condition applied as a traction condition on the
 *  momentum equations.
 *
 *      t . (v -vs) = 0
 *
 *
 *  where v is the actual velocity = x_dot + v
 *        vs is the substrate velocity
 *
 *            Author: P. R. Schunk    (5/17/94)
 *            Revised: 6/1/95 RAC
 ************************************************************************/
{
  double beta = bc_float[0]; /* Navier slip coefficient from input deck */
  /* velocity components of solid surface on
   * which slip condition is applied */
  double vsx = bc_float[1];
  double vsy = bc_float[2];
  double vsz = bc_float[3];
  double alpha = bc_float[4];
  int a, j, var, jvar, p, dim;
  double phi_j, vs[MAX_PDIM];
  double slip_dir[MAX_PDIM], vslip[MAX_PDIM], vrel[MAX_PDIM], vrel_dotn;
  double X_0[3], omega, velo_avg = 0., pgrad = 0., thick = 0., dthick_dV = 0., dthick_dP = 0.;
  int Pflag = TRUE;
  double pg_factor = 1.0, tang_sgn = 1.0;

  int icount;
  double dist;    /* distance btw current position and dynamic CL */
  double betainv; /* inverse of slip coefficient */
  /* double dot_prod; */
  double d_betainv_dvslip_mag, d_betainv_dP;
  double d_betainv_dF[MDE];
  double sign;
  int tang_slip_only;
#define PRESSURE_DEPENDENT_SLIP 0
#if PRESSURE_DEPENDENT_SLIP
  double vslip_mag;
#endif

#define TANGENT_SLIP_ONLY 0
#if TANGENT_SLIP_ONLY
  tang_slip_only = TRUE;
#else
  tang_slip_only = FALSE;
#endif

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  dim = pd->Num_Dim;

  vs[0] = vsx;
  vs[1] = vsy;
  vs[2] = vsz;

  /* dependence of slip coefficient on pressure and magnitude of slip */
  d_betainv_dvslip_mag = 0.;
  d_betainv_dP = 0.;

  if (beta != 0.) {
    betainv = 1. / beta;
  } else {
    betainv = 0.;
  }

  if (ls != NULL) {
    load_lsi(ls->Length_Scale);
    if (af->Assemble_Jacobian)
      load_lsi_derivs();
    /* dot_prod = fabs( dot_product( pd->Num_Dim, lsi->normal, fv->snormal) ); */
  }

  if (type == VELO_SLIP_FILL_BC || type == VELO_SLIP_ROT_FILL_BC) {
    /* This is for level sets. */
    if (ls == NULL) {
      /* This is for volume of fluid. */
      if (fabs(fv->F) > 0.25) {
        betainv = LITTLE_PENALTY;
      } else {
        for (a = 0; a < dim; a++) {
          beta += fv->grad_F[a] * fv->grad_F[a];
        }
        betainv = 1. / beta;
      }
    }
  }

  /*
   *   Redefine the Substrate velocity from scratch
   *   for rotational velocity, if needed.
   *   HKM -> Note I have not put in x_dot into these equations
   *          even though it might be needed here.
   */
  /* NB This is really specified on for 2D
     In 3D it spins with omega pointing in z
     direction about x0,y0  */
  /* Note: positive omega is CLOCKWISE */
  if (type == VELO_SLIP_ROT_BC || type == VELO_SLIP_ROT_FILL_BC || type == VELO_SLIP_ROT_FLUID_BC) {
    double factor = 1.0, current_rad, rad_input = bc_float[5];
    omega = vsx;
    X_0[0] = vsy;
    X_0[1] = vsz;
    if (rad_input > 0) {
      current_rad = sqrt(SQUARE(fv->x[1] - X_0[1]) + SQUARE(fv->x[0] - X_0[0]));
      factor = rad_input / current_rad;
    } else {
      factor = 1.0;
    }
    vs[0] = factor * omega * (fv->x[1] - X_0[1]);
    vs[1] = -factor * omega * (fv->x[0] - X_0[0]);
    vs[2] = 0.;
  } /* if: VELO_SLIP_ROT_BC */

  memset(vrel, 0, sizeof(double) * MAX_PDIM);
  for (p = 0; p < pd->Num_Dim; p++) {
    vrel[p] = (vs[p] - fv->v[p]);
  }

  if (tang_slip_only) {
    for (p = 0; p < pd->Num_Dim; p++) {
      vslip[p] = 0.0;
      for (a = 0; a < pd->Num_Dim; a++) {
        vslip[p] += vrel[a] * ((double)delta(a, p) - fv->snormal[a] * fv->snormal[p]);
      }
    }
    vrel_dotn = 0.;
    for (a = 0; a < pd->Num_Dim; a++) {
      vrel_dotn += vslip[a] * fv->snormal[a];
    }
  } else {
    for (p = 0; p < pd->Num_Dim; p++) {
      vslip[p] = vrel[p];
    }
  }

  for (p = 0; p < pd->Num_Dim; p++) {
    slip_dir[p] = vslip[p];
  }

#if PRESSURE_DEPENDENT_SLIP
  vslip_mag = normalize_really_simple_vector(slip_dir, pd->Num_Dim);

  if (fv->P <= 0.) {
    betainv = 0.;
    d_betainv_dvslip_mag = 0.;
    d_betainv_dP = 0.;
    memset(d_betainv_dF, 0, sizeof(double) * MDE);
  } else {
    /* limit the force to the strong integrated condition limit */
    if (vslip_mag <= beta * fv->P / BIG_PENALTY) {
      betainv = BIG_PENALTY;
      d_betainv_dvslip_mag = 0.;
      d_betainv_dP = 0.;
      memset(d_betainv_dF, 0, sizeof(double) * MDE);
    } else {
      betainv = beta * fv->P / vslip_mag;
      d_betainv_dvslip_mag = -betainv / vslip_mag;
      d_betainv_dP = beta / vslip_mag;
      memset(d_betainv_dF, 0, sizeof(double) * MDE);
    }
  }
#endif /* PRESSURE_DEPENDENT_SLIP */

  /****   compute slip parameter   ***/

  /***************************************************************************
   *                     COMPUTE SLIP PARAMETER
   *
   * This section is for position dependent slip. Calculate position
   * of dynamic contact line from reference node do the distance calculation
   * based on undeformed geometry... This saves us Jacobian entries in
   * uncharted sections of the A matrix
   *
   * Calculate distance from dcl to current Gauss point.  Turn this
   * off if we don't have position depend slip.  alpha will be zero,
   * so dist can be anything as long as it is
   * defined. Protect exponentially decaying slip from underflow far
   * away from the singularity.
   *
   ***************************************************************************/

  if (type == VELO_SLIP_FLUID_BC || type == VELO_SLIP_ROT_FLUID_BC) {
    dist = 0.;
    for (icount = 0; icount < pd->Num_Dim; icount++) {
      /* Uses undeformed node position */
      dist += SQUARE(xsurf[icount] - bc_float[max_float + 1 + icount]);
    }
    dist /= SQUARE(bc_float[6]);
  } else if (alpha != 0.) {
    if ((type == VELO_SLIP_FILL_BC || type == VELO_SLIP_ROT_FILL_BC) && ls != NULL) {
      dist = fabs(fv->F);
      /*dist = fabs(fv->F) / sqrt(1.0 - dot_prod*dot_prod);*/
    } else {
      /* Coord position in bc_float, from BC_Data_Float */
      dist = 0.;
      for (icount = 0; icount < pd->Num_Dim; icount++) {
        /* Uses undeformed node position */
        dist += (xsurf[icount] - bc_float[icount + max_float + 1]) *
                (xsurf[icount] - bc_float[icount + max_float + 1]);
      }
      dist = sqrt(dist);
    }
  } else {
    dist = 1.0;
  }

  /* for exponentially decaying slip, max out betainv when equivalent to
   * STRONG IC for no slip
   */
  if (alpha * dist > log(beta * BIG_PENALTY)) {
    betainv = BIG_PENALTY;
    d_betainv_dvslip_mag = 0.;
    d_betainv_dP = 0.;
    memset(d_betainv_dF, 0, sizeof(double) * MDE);
  } else {
    betainv *= exp(alpha * dist);
    memset(d_betainv_dF, 0, sizeof(double) * MDE);
    if ((type == VELO_SLIP_FILL_BC || type == VELO_SLIP_ROT_FILL_BC) && ls != NULL) {
      /*
       * NB, here: dist = |F| = |FILL|
       * Before  : betainv = 1 / beta
       * Now     : betainv = exp(alpha*|F|) / beta
       *
       * Need to update the sensitivity of betainv w.r.t. FILL.
       *
       * dist = |F| = sign(F) * F = sign(F) * SUM_j ( phi_j * F_j )
       *
       * If F == 0, then exp(alpha*dist) == 1 so we'll just take sign = 0.0.
       */
      sign = fv->F == 0.0 ? 0.0 : (fv->F > 0.0 ? +1.0 : -1.0);
      for (j = 0; j < ei[pg->imtrx]->dof[FILL]; j++) {
        d_betainv_dF[j] = betainv * sign * alpha * bf[FILL]->phi[j];
      }
#if 0
	  DPRINTF(stderr,"VELO_SLIP_FILL x=(%g,%g), fv->F=%g, betainv=%g, alpha*dist=%g\n", 
		  fv->x[0], fv->x[1], fv->F, betainv, alpha*dist);
#endif
    }
  }

  /* quantities specific to FLUID bcs   */

  Pflag = (int)bc_float[7];
  velo_avg = 0.0;
  pgrad = 0.;
  if (type == VELO_SLIP_FLUID_BC || type == VELO_SLIP_ROT_FLUID_BC) {
    double v_solid = 0., res, jac, delta, flow, eps = 1.0e-8, viscinv;
    double jacinv, v_mag = 0.;
    tang_sgn = 0.;
    for (p = 0; p < pd->Num_Dim; p++) {
      tang_sgn += fv->stangent[0][p] * vs[p];
      velo_avg += fv->stangent[0][p] * (vs[p] + fv->v[p]);
      v_solid += fv->stangent[0][p] * vs[p];
      v_mag += SQUARE(vs[p]);
      if (Pflag) {
        pgrad += fv->stangent[0][p] * fv->grad_P[p];
      }
    }
    v_mag = sqrt(v_mag);
    tang_sgn = v_solid / v_mag;
    tang_sgn = (double)SGN(v_solid / v_mag);
    velo_avg *= 0.5;
    /* sometimes the tangent/normals flip causing havoc....*/
    if (v_solid < 0) {
      GOMA_WH(-1, "fvelo_slip: normals and tangents have flipped! - try CONTACT_LINE model\n");
      velo_avg *= tang_sgn;
      v_solid *= tang_sgn;
      pgrad *= tang_sgn;
    }
#if 1
    if (dist < 10.0) {
      pg_factor = 1.0 - exp(-dist);
    } else {
      pg_factor = 1.0;
    }
#else
    if (dist < 1.0) {
      pg_factor = 0.0;
    } else if (dist < 2.0) {
      pg_factor = (dist - 1.0);
    } else {
      pg_factor = 1.0;
    }
#endif
    pgrad *= pg_factor;

    flow = MAX(0., bc_float[4] * v_solid);
    viscinv = 1. / bc_float[0];
    thick = flow / velo_avg;
    j = 0;
    do {
      res = -CUBE(thick) * viscinv * pgrad / 12. + thick * velo_avg - flow;
      jac = -0.25 * SQUARE(thick) * viscinv * pgrad + velo_avg;
      jacinv = 1.0 / jac;
      delta = -res * jacinv;
      thick += delta;
      j++;
    } while (fabs(delta) > eps && j < 20);
    if (thick > DBL_SMALL) {
      betainv = bc_float[0] / thick;
      dthick_dV = -0.5 * jacinv; /*  1/h*derivative  */
      dthick_dP = CUBE(thick) * viscinv / 12. * jacinv;
    } else {
      betainv = LITTLE_PENALTY;
    }
#if 0
fprintf(stderr,"slip %d %g %g %g %g\n",Pflag,fv->x[0],thick,flow/v_solid,velo_avg);
fprintf(stderr,"more %g %g %g %g\n",res,jac,betainv, dthick_dV);
#endif
  }
  /* Calculate the residual contribution. */
  for (p = 0; p < pd->Num_Dim; p++) {
    func[p] += betainv * vslip[p];
  }
  if (Pflag && (type == VELO_SLIP_FLUID_BC || type == VELO_SLIP_ROT_FLUID_BC)) {
    for (p = 0; p < pd->Num_Dim; p++) {
      func[p] += 0.5 * thick * pg_factor * fv->grad_P[p];
    }
  }

  /*
   * The Jacobian calculation is split into two parts.
   * The first part computes the derivatives of func for the
   * whole vslip vector.  The second part, which is done only
   * if tang_slip_only == TRUE, computes the derivatives of the
   * normal component of vslip which are subtracted off.
   */

  if (af->Assemble_Jacobian) {

    for (jvar = 0; jvar < pd->Num_Dim; jvar++) {
      var = VELOCITY1 + jvar;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          /* Main dependence of the velocity on velocity unknowns */
          d_func[jvar][var][j] += (-betainv) * phi_j;
          /* contribution from dependence of betainv on the magnitude of slip */
          for (p = 0; p < pd->Num_Dim; p++) {
            d_func[p][var][j] += (-d_betainv_dvslip_mag) * phi_j * (slip_dir[jvar] * vslip[p]);
          }
        }
      }
    }
    if ((type == VELO_SLIP_FLUID_BC || type == VELO_SLIP_ROT_FLUID_BC)) {
      for (jvar = 0; jvar < pd->Num_Dim; jvar++) {
        var = VELOCITY1 + jvar;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            d_func[jvar][var][j] +=
                -betainv * dthick_dV * tang_sgn * fv->stangent[0][jvar] * vslip[jvar] * phi_j;
            /* don't think we need -- seems OK 4-12-2017 */
            if (Pflag)
              d_func[jvar][var][j] += 0.5 * thick * dthick_dV * tang_sgn * fv->stangent[0][jvar] *
                                      pg_factor * fv->grad_P[jvar] * phi_j;
          }
        }
      }

      /* Mesh motion Jacobian entries   */
      for (jvar = 0; jvar < ei[pg->imtrx]->ielem_dim; jvar++) {
        var = MESH_DISPLACEMENT1 + jvar;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            for (p = 0; p < pd->Num_Dim; p++) {
              d_func[p][var][j] +=
                  -betainv * dthick_dV * vslip[p] * tang_sgn * fv->dstangent_dx[0][p][jvar][j];
              if (Pflag) {
                d_func[p][var][j] += 0.5 * thick * dthick_dV * pg_factor * fv->grad_P[p] *
                                     tang_sgn * fv->dstangent_dx[0][p][jvar][j];
                d_func[p][var][j] += 0.5 * pg_factor * fv->grad_P[p] * dthick_dP * tang_sgn *
                                     fv->dstangent_dx[0][p][jvar][j];
                d_func[p][var][j] += 0.5 * pg_factor * fv->d_grad_P_dmesh[p][jvar][j] * thick;
                d_func[p][var][j] +=
                    vslip[p] * dthick_dP * tang_sgn * fv->dstangent_dx[0][p][jvar][j] * (-betainv);
              }
            }
          }
        }
      }
    }

    for (jvar = 0; jvar < pd->Num_Dim; jvar++) {
      var = PVELOCITY1 + jvar;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          d_func[jvar][var][j] += (-betainv) * phi_j;
          /* contribution from dependence of betainv on the magnitude of slip */
          for (p = 0; p < pd->Num_Dim; p++) {
            d_func[p][var][j] += (-d_betainv_dvslip_mag) * phi_j * (slip_dir[jvar] * vslip[p]);
          }
        }
      }
    }

    var = PRESSURE;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        phi_j = bf[var]->phi[j];
        for (p = 0; p < pd->Num_Dim; p++) {
          d_func[p][var][j] += d_betainv_dP * vslip[p] * phi_j;
        }
      }
      if (Pflag && (type == VELO_SLIP_FLUID_BC || type == VELO_SLIP_ROT_FLUID_BC)) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (p = 0; p < pd->Num_Dim; p++) {
            phi_j = bf[var]->grad_phi[j][p];
            d_func[p][var][j] += 0.5 * thick * phi_j;
            d_func[p][var][j] +=
                0.5 * pg_factor * fv->grad_P[p] * dthick_dP * tang_sgn * fv->stangent[0][p] * phi_j;
            d_func[p][var][j] +=
                vslip[p] * (-betainv) * dthick_dP * tang_sgn * fv->stangent[0][p] * phi_j;
          }
        }
      }
    }

    var = FILL;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        for (p = 0; p < pd->Num_Dim; p++) {
          d_func[p][var][j] += d_betainv_dF[j] * vslip[p];
        }
      }
    }

    /*
     * Special mesh derivatives for level set problems on a moving mesh.
     * --> delta(F) has moving mesh derivatives
     */
    if ((type == VELO_SLIP_FILL_BC || type == VELO_SLIP_ROT_FILL_BC) &&
        pd->v[pg->imtrx][MESH_DISPLACEMENT1] && ls != NULL) {
      for (jvar = 0; jvar < ei[pg->imtrx]->ielem_dim; jvar++) {
        var = MESH_DISPLACEMENT1 + jvar;
        if (pd->v[pg->imtrx][var] && lsi->near) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            for (p = 0; p < pd->Num_Dim; p++) {
              d_func[p][var][j] += -func[p] * lsi->d_delta_dmesh[jvar][j] / lsi->delta;
            }
          }
        }
      }
    }

    /*
     * Extra terms if we are limiting slip to be tangent to surface.
     * specifically, there are mesh dependencies and additional velocity
     * terms.
     */
    if (tang_slip_only) {
      for (jvar = 0; jvar < ei[pg->imtrx]->ielem_dim; jvar++) {
        var = MESH_DISPLACEMENT1 + jvar;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            for (p = 0; p < pd->Num_Dim; p++) {
              for (a = 0; a < pd->Num_Dim; a++) {
                /* Derivatives of the normal vector. */
                d_func[p][var][j] += (-betainv) * vrel[a] *
                                     (fv->snormal[jvar] * fv->dsnormal_dx[a][jvar][j] +
                                      fv->snormal[a] * fv->dsnormal_dx[p][jvar][j]);

                /* Contribution from dependence of betainv on the magnitude of slip. */
                d_func[p][var][j] += (-d_betainv_dvslip_mag) * vrel_dotn * vslip[p] *
                                     (slip_dir[a] * fv->dsnormal_dx[a][jvar][j]);
              }
            }
          }
        }
      }

      for (jvar = 0; jvar < pd->Num_Dim; jvar++) {
        var = VELOCITY1 + jvar;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            for (p = 0; p < pd->Num_Dim; p++) {
              d_func[p][var][j] += (-betainv) * phi_j * (-fv->snormal[jvar] * fv->snormal[p]);
            }
          }
        }
      }

      for (jvar = 0; jvar < pd->Num_Dim; jvar++) {
        var = PVELOCITY1 + jvar;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            for (p = 0; p < pd->Num_Dim; p++) {
              d_func[p][var][j] += (-betainv) * phi_j * (-fv->snormal[jvar] * fv->snormal[p]);
            }
          }
        }
      }
    }
  }

  return;

} /* END of routine fvelo_slip_bc  */

void fvelo_slip_power_bc(double func[MAX_PDIM],
                         double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                         const int type,      /* whether rotational or not */
                         const int max_float, /* Max float number from input file */
                         double bc_float[MAX_BC_FLOAT_DATA],
                         const double tt, /* parameter in time stepping alg           */
                         const double dt) /* current time step value                  */

/*************************************************************************
 *
 *   n . T = beta * ( t . (v - vs) )^m
 *
 *   Where t is a tangent in the wall direction
 *
 ************************************************************************/
{
  double beta = bc_float[0]; /* Navier slip coefficient from input deck */
  /* velocity components of solid surface on
   * which slip condition is applied */
  double vsx = bc_float[1];
  double vsy = bc_float[2];
  double vsz = bc_float[3];
  double expon = bc_float[4];
  int j, var, jvar, p;
  double phi_j, vs[MAX_PDIM];
  double betainv; /* inverse of slip coefficient */
  double vslip[MAX_PDIM];
  double tangent_dot_vslip = 0.0;
  double tangent[MAX_PDIM];

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  vs[0] = vsx;
  vs[1] = vsy;
  vs[2] = vsz;

  int constant_tangent = FALSE;

  if (max_float == 8) {
    tangent[0] = bc_float[5];
    tangent[1] = bc_float[6];
    tangent[2] = bc_float[7];
    constant_tangent = TRUE;
  } else if (pd->Num_Dim == 3 && type == VELO_SLIP_POWER_BC) {
    GOMA_EH(GOMA_ERROR, "Must provide constant tangent for VELO_SLIP_POWER");
  } else {
    tangent[0] = fv->stangent[0][0];
    tangent[1] = fv->stangent[0][1];
    tangent[2] = 0.0;
  }

  if (beta != 0.) {
    betainv = 1. / beta;
  } else {
    betainv = 0.;
  }

  memset(vslip, 0, sizeof(double) * MAX_PDIM);
  for (p = 0; p < pd->Num_Dim; p++) {
    vslip[p] = (fv->v[p] - vs[p]);
  }

  /* Calculate the residual contribution. */
  for (p = 0; p < pd->Num_Dim; p++) {
    tangent_dot_vslip += tangent[p] * vslip[p];
  }

  double tdotv_power = pow(tangent_dot_vslip, expon);

  if (type == VELO_SLIP_POWER_BC) {
    // tangential compoment only
    func[0] = -betainv * tdotv_power;

    if (af->Assemble_Jacobian) {

      for (jvar = 0; jvar < pd->Num_Dim; jvar++) {
        var = VELOCITY1 + jvar;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            /* Main dependence of the velocity on velocity unknowns */
            if (fabs(tangent_dot_vslip) > 1e-15) {
              d_func[0][var][j] +=
                  -betainv * expon * tangent[jvar] * phi_j * tdotv_power / tangent_dot_vslip;
            }
          }
        }
      }

      /* Mesh motion Jacobian entries   */
      if (constant_tangent == FALSE) {
        for (jvar = 0; jvar < ei[pg->imtrx]->ielem_dim; jvar++) {
          var = MESH_DISPLACEMENT1 + jvar;
          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              for (p = 0; p < pd->Num_Dim; p++) {
                if (fabs(tangent_dot_vslip) > 0) {
                  d_func[0][var][j] += -betainv * expon * fv->dstangent_dx[0][p][jvar][j] *
                                       vslip[p] * tdotv_power / tangent_dot_vslip;
                  ;
                }
              }
            }
          }
        }
      }
    }
  } else if (type == VELO_SLIP_POWER_CARD_BC) {
    double vslip_expon[MAX_PDIM];
    for (p = 0; p < VIM; p++) {
      vslip_expon[p] = pow(vslip[p], expon);
      func[p] = -betainv * vslip_expon[p];
    }
    if (af->Assemble_Jacobian) {
      for (jvar = 0; jvar < pd->Num_Dim; jvar++) {
        var = VELOCITY1 + jvar;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            /* Main dependence of the velocity on velocity unknowns */
            if (fabs(vslip[jvar]) > 0) {
              d_func[jvar][var][j] += -betainv * expon * phi_j * vslip_expon[jvar] / vslip[jvar];
            }
          }
        }
      }
    }
  } else {
    GOMA_EH(GOMA_ERROR, "Unknown type for fvelo_slip_power_bc");
  }

  return;
}
/**
 * Exchanges coordinates for the reference node needed for calculation
 * in fvelo_slip_bc()
 *
 * Returns 0 if success
 * Returns -1 if error
 */
int exchange_fvelo_slip_bc_info(int ibc /* Index into BC_Types for VELO_SLIP_BC */) {
  int i;
#ifdef PARALLEL
  int mpi_error;
#endif
  /* if velo slip has */
  int velo_slip_root = 0;
  /* Offset for where to place coordinates in BC_Data_Float */
  int float_offset;
  double node_coord[DIM]; /* temporary buffer for node coordinates */

  float_offset = BC_Types[ibc].max_DFlt + 1;
  /* Skip this if calculation is not needed */
  if (BC_Types[ibc].BC_Data_Int[0] == -1) {
    return 0;
  }

  /* find which processor has the right data */
  if (BC_Types[ibc].BC_Data_Int[0] != -2) {
    velo_slip_root = ProcID;
  }

#ifdef PARALLEL
  mpi_error = MPI_Allreduce(MPI_IN_PLACE, &velo_slip_root, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (mpi_error != MPI_SUCCESS) {
    GOMA_EH(GOMA_ERROR, "Error in MPI Allreduce");
    return -1;
  }
#endif /* #ifdef PARALLEL */

  if (ProcID == velo_slip_root) {
    int node = abs(BC_Types[ibc].BC_Data_Int[0]);
    /* find coordinate position of reference node */
    for (i = 0; i < pd->Num_Dim; i++) {
      node_coord[i] = Coor[i][node];
    }
  }

#ifdef PARALLEL
  /* Communicate the new values in BC_Data_Float */
  mpi_error = MPI_Bcast(&node_coord, pd->Num_Dim, MPI_DOUBLE, velo_slip_root, MPI_COMM_WORLD);

  if (mpi_error != MPI_SUCCESS) {
    GOMA_EH(GOMA_ERROR, "Error in MPI Allreduce");
    return -1;
  }
#endif

  /* set BC_Data_Float values */
  if (1 || BC_Types[ibc].BC_Name != ROLL_FLUID_BC) {
    for (i = 0; i < pd->Num_Dim; i++) {
      BC_Types[ibc].BC_Data_Float[float_offset + i] = node_coord[i];
    }
  } else {
    for (i = 0; i < pd->Num_Dim; i++) {
      BC_Types[ibc].u_BC[float_offset + i] = node_coord[i];
    }

    fprintf(stderr, "exchange %d %g %g\n", ibc, node_coord[0], node_coord[1]);
  }
  return 0;
}

/****************************************************************************/
/***************************************************************************/

void fvelo_airfilm_bc(double func[MAX_PDIM],
                      double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                      double x[],
                      const int type, /* whether rotational or not */
                      double bc_float[MAX_BC_FLOAT_DATA],
                      const int dcl_node,           /*   node id for DCL  */
                      const double xsurf[MAX_PDIM], /* coordinates of surface Gauss  *
                                                     * point, i.e. current position  */
                      const double tt,              /* parameter in time stepping alg           */
                      const double dt)              /* current time step value                  */

/*************************************************************************
 *
 *  Function which evaluates a "slip-like" velocity boundary condition
 *  except that it mimics an entrained air film
 *
 *            Author: R. B. Secor    (10/2/2015)
 *            Revised:
 ************************************************************************/
{
  double gas_mu = bc_float[0]; /* Navier slip coefficient from input deck */
  /* velocity components of solid surface on
   * which slip condition is applied */
  double vsx = bc_float[1];
  double vsy = bc_float[2];
  double vsz = bc_float[3];
  double gas_flow = bc_float[4];
  int j, var, jvar, p, dim;
  double phi_j, vs[MAX_PDIM];
  double vnet[MAX_PDIM], vrel[MAX_PDIM];
  double X_0[3], omega;

  double vnet_mag, gradP_mag; /* inverse of slip coefficient */

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  dim = pd->Num_Dim;

  vs[0] = vsx;
  vs[1] = vsy;
  vs[2] = vsz;

  if (TimeIntegration == TRANSIENT && pd->e[pg->imtrx][R_MESH1]) {
    /* Add the mesh motion to the substrate velocity */
    vs[0] += fv_dot->x[0];
    vs[1] += fv_dot->x[1];
    vs[2] += fv_dot->x[2];
  }

  /*
   *   Redefine the Substrate velocity from scratch
   *   for rotational velocity, if needed.
   *   HKM -> Note I have not put in x_dot into these equations
   *          even though it might be needed here.
   */
  /* NB This is really specified on for 2D
     In 3D it spins with omega pointing in z
     direction about x0,y0  */
  /* Note: positive omega is CLOCKWISE */
  if (type == AIR_FILM_ROT_BC) {
    omega = vsx;
    X_0[0] = vsy;
    X_0[1] = vsz;
    vs[0] = omega * (fv->x[1] - X_0[1]);
    vs[1] = -omega * (fv->x[0] - X_0[0]);
    vs[2] = 0.;
  } /* if: AIR_FILM_ROT_BC */

  memset(vrel, 0, sizeof(double) * MAX_PDIM);
  memset(vnet, 0, sizeof(double) * MAX_PDIM);
  vnet_mag = 0.;
  gradP_mag = 0.;
  for (p = 0; p < dim; p++) {
    vrel[p] = (vs[p] - fv->v[p]);
    vnet[p] = (vs[p] + fv->v[p]);
    vnet_mag += vnet[p];
    gradP_mag += fv->grad_P[p];
  }

  /* Calculate the residual contribution. */
  /*  Drop through special cases for the solution of
   *    the air film thickness  */
  if (gradP_mag < DBL_SMALL) {
    if (fabs(gas_flow) < DBL_SMALL) {
      for (p = 0; p < dim; p++) {
        func[p] += gas_mu * vrel[p];
      }
      if (af->Assemble_Jacobian) {

        for (jvar = 0; jvar < dim; jvar++) {
          var = VELOCITY1 + jvar;
          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];
              d_func[jvar][var][j] += (-gas_mu) * phi_j;
            }
          }
        }
      }
    } else {
      for (p = 0; p < dim; p++) {
        func[p] += gas_mu * vrel[p] * vnet[p] / gas_flow;
      }
      if (af->Assemble_Jacobian) {

        for (jvar = 0; jvar < dim; jvar++) {
          var = VELOCITY1 + jvar;
          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];
              for (p = 0; p < dim; p++) {
                d_func[p][var][j] += (gas_mu / gas_flow) * (vrel[p] - vnet[p]) * phi_j;
              }
            }
          }
        }
      }
    }
  } else if (vnet_mag < DBL_SMALL) {
    if (fabs(gas_flow) < DBL_SMALL) {
      for (p = 0; p < dim; p++) {
        func[p] += gas_mu * vrel[p];
      }
      if (af->Assemble_Jacobian) {

        for (jvar = 0; jvar < dim; jvar++) {
          var = VELOCITY1 + jvar;
          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];
              d_func[jvar][var][j] += (-gas_mu) * phi_j;
            }
          }
        }
      }
    } else {
      for (p = 0; p < dim; p++) {
        func[p] += pow(1.5 * gas_mu * gas_flow * SQUARE(fv->grad_P[p]), 1. / 3.) +
                   vrel[p] * pow(fv->grad_P[p] * SQUARE(gas_mu) / 12. / gas_flow, 1. / 3.);
      }
      if (af->Assemble_Jacobian) {

        for (jvar = 0; jvar < dim; jvar++) {
          var = VELOCITY1 + jvar;
          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];
              for (p = 0; p < dim; p++) {
                d_func[p][var][j] +=
                    pow(fv->grad_P[p] * SQUARE(gas_mu) / 12. / gas_flow, 1. / 3.) * (-phi_j);
              }
            }
          }
        }
        var = PRESSURE;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            for (p = 0; p < dim; p++) {
              d_func[p][var][j] +=
                  (2. / 3. * pow(1.5 * gas_mu * gas_flow / fv->grad_P[p], 1. / 3.) +
                   1. / 3. * vrel[p] * pow(SQUARE(gas_mu) / 12. / gas_flow, 1. / 3.)) *
                  bf[var]->grad_phi[j][p];
            }
          }
        }
        for (jvar = 0; jvar < dim; jvar++) {
          var = MESH_DISPLACEMENT1 + jvar;
          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];
              for (p = 0; p < dim; p++) {
                d_func[p][var][j] +=
                    (2. / 3. * pow(1.5 * gas_mu * gas_flow / fv->grad_P[p], 1. / 3.) +
                     1. / 3. * vrel[p] * pow(SQUARE(gas_mu) / 12. / gas_flow, 1. / 3.)) *
                    fv->d_grad_P_dmesh[p][jvar][j];
              }
            }
          }
        }
      }
    }
  } else if (fabs(gas_flow) < DBL_SMALL) {
    if (vnet_mag < DBL_SMALL) {
      for (p = 0; p < dim; p++) {
        func[p] += gas_mu * vrel[p];
      }
      if (af->Assemble_Jacobian) {

        for (jvar = 0; jvar < dim; jvar++) {
          var = VELOCITY1 + jvar;
          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];
              d_func[jvar][var][j] += (-gas_mu) * phi_j;
            }
          }
        }
      }
    } else {
      for (p = 0; p < dim; p++) {
        func[p] += sqrt(gas_mu * fabs(fv->grad_P[p] / vnet[p]) *
                        (9 * SQUARE(vnet[p]) + SQUARE(vrel[p])) / 6);
      }
      if (af->Assemble_Jacobian) {

        for (jvar = 0; jvar < dim; jvar++) {
          var = VELOCITY1 + jvar;
          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];
              for (p = 0; p < dim; p++) {
                d_func[p][var][j] +=
                    sqrt(gas_mu * fabs(fv->grad_P[p])) * phi_j *
                    (18. * vnet[p] * vrel[p] - SQUARE(vrel[p]) + 9 * SQUARE(vnet[p])) /
                    (6 * SQUARE(vrel[p])) * 0.5 *
                    sqrt(6 * vnet[p] / (9 * SQUARE(vnet[p]) + SQUARE(vrel[p])));
              }
            }
          }
        }
        var = PRESSURE;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            for (p = 0; p < dim; p++) {
              d_func[p][var][j] +=
                  sqrt(gas_mu / vnet[p] * (9 * SQUARE(vnet[p]) + SQUARE(vrel[p])) / 6) * 0.5 /
                  sqrt(fabs(fv->grad_P[p])) * bf[var]->grad_phi[j][p];
            }
          }
        }
        for (jvar = 0; jvar < dim; jvar++) {
          var = MESH_DISPLACEMENT1 + jvar;
          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];
              for (p = 0; p < dim; p++) {
                d_func[p][var][j] +=
                    sqrt(gas_mu / vnet[p] * (9 * SQUARE(vnet[p]) + SQUARE(vrel[p])) / 6) * 0.5 /
                    sqrt(fabs(fv->grad_P[p])) * fv->d_grad_P_dmesh[p][jvar][j];
              }
            }
          }
        }
      }
    }
  } else /*  General Solution */
  {
    for (p = 0; p < dim; p++) {
      func[p] += gas_mu * vrel[p] * vnet[p] / gas_flow;
    }
  }

  return;

} /* END of routine fvelo_airfilm_bc  */
/****************************************************************************/
/****************************************************************************/

#if 1
void fvelo_slip_level(double func[MAX_PDIM],
                      double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                      const int type, /* whether rotational or not */
                      double width,
                      double beta_inside,        /* slip coefficient from input deck */
                      const double vsx,          /* velocity components of solid  */
                      const double vsy,          /* surface on which slip condition   */
                      const double vsz,          /* is applied           */
                      const double beta_outside, /* slip coeff away from of interface*/
                      const double gas_phase_factor,
                      const double contact_fraction,
                      const double tau,
                      const double tt,
                      const double dt) {
  int j, var, jvar, p, q, b;
  double phi_j, vs[MAX_PDIM];
  double beta, betainv;
  double d_beta_dF[MDE];
  double X_0[3], omega;

  dbl Pi[DIM][DIM];
  STRESS_DEPENDENCE_STRUCT d_Pi_struct;
  STRESS_DEPENDENCE_STRUCT *d_Pi = &d_Pi_struct;

  /************************* EXECUTION BEGINS *******************************/

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  if (!af->Assemble_Jacobian)
    d_Pi = NULL;

  load_lsi(width);
  if (af->Assemble_Jacobian) {
    load_lsi_derivs();
    memset(d_beta_dF, 0, MDE * sizeof(double));
  }

  beta = slip_coefficient(beta_inside, beta_outside, fv->F, width, d_beta_dF, gas_phase_factor,
                          contact_fraction);

  betainv = 1.0 / beta;

  vs[0] = vsx;
  vs[1] = vsy;
  vs[2] = vsz;

  if (type == VELO_SLIP_LS_ROT_BC) {
    omega = vsx;
    X_0[0] = vsy;
    X_0[1] = vsz;
    vs[0] = omega * (fv->x[1] - X_0[1]);
    vs[1] = -omega * (fv->x[0] - X_0[0]);
    vs[2] = 0.;
  }

  if (TimeIntegration == TRANSIENT && pd->e[pg->imtrx][R_MESH1]) {
    /* Add the mesh motion to the substrate velocity */
    vs[0] += fv_dot->x[0];
    vs[1] += fv_dot->x[1];
    vs[2] += fv_dot->x[2];
  }

  /* compute stress tensor and its derivatives */
  if (type == VELO_SLIP_LEVEL_SIC_BC) {
    fluid_stress(Pi, d_Pi);
  }

  if (af->Assemble_Jacobian) {
    for (jvar = 0; jvar < pd->Num_Dim; jvar++) {
      var = VELOCITY1 + jvar;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          d_func[jvar][var][j] += (-betainv) * (tau * (1. + 2. * tt) * phi_j / dt + phi_j);
        }
      }
    }

#if 1
    if (type == VELO_SLIP_LEVEL_SIC_BC) {
      for (p = 0; p < pd->Num_Dim; p++) {
        for (q = 0; q < pd->Num_Dim; q++) {
          for (b = 0; b < VIM; b++) {
            var = VELOCITY1 + b;
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_func[p][var][j] -= fv->snormal[q] * d_Pi->v[p][q][b][j];
            }
          }
        }
      }
      var = PRESSURE;
      if (pd->v[pg->imtrx][var]) {
        for (p = 0; p < pd->Num_Dim; p++) {
          for (q = 0; q < pd->Num_Dim; q++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_func[p][var][j] -= fv->snormal[q] * d_Pi->P[p][q][j];
            }
          }
        }
      }
    }
#endif

    var = LS;
    if (pd->v[pg->imtrx][var]) {
      for (p = 0; p < pd->Num_Dim; p++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];

          d_func[p][var][j] += (d_beta_dF[j] / beta / beta) * (fv->v[p] - vs[p]);
        }
      }
    }

#if 1
    if (type == VELO_SLIP_LEVEL_SIC_BC && pd->v[pg->imtrx][var]) {
      for (p = 0; p < pd->Num_Dim; p++) {
        for (q = 0; q < pd->Num_Dim; q++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_func[p][var][j] -= fv->snormal[q] * d_Pi->F[p][q][j];
          }
        }
      }
    }
#endif
  } /* end of of Assemble Jacobian		*/

  /* Calculate the residual contribution	*/
  for (p = 0; p < pd->Num_Dim; p++) {
    func[p] += (-betainv) * (tau * fv_dot->v[p] + fv->v[p] - vs[p]);
    if (type == VELO_SLIP_LEVEL_SIC_BC) {
      for (q = 0; q < pd->Num_Dim; q++) {
        func[p] -= fv->snormal[q] * Pi[p][q];
      }
    }
  }

} /* END of routine fvelo_slip_level  */

static double slip_coefficient(const double beta0,
                               const double beta1,
                               const double F,
                               const double width,
                               double *d_beta_dF,
                               const double gas_phase_factor,
                               const double contact_fraction) {
  double beta;
  int j;
  double sign_gas = mp->mp2nd->densitymask[1] - mp->mp2nd->densitymask[0];
  /* 	const double gas_phase_factor = 8;    multiplier for enhanced gas slip */
  /* 	const double contact_fraction = 1e-6; fraction for contact condition */
  const double beta0_log = log(beta0);
  const double beta1_log = log(beta1);

  /* interpolate based on log scale since the slip coefficient normally
     varies over many orders of magnitude  */

  beta = exp((beta0_log - beta1_log) * lsi->delta / lsi->delta_max + beta1_log);
  if (af->Assemble_Jacobian) {
    for (j = 0; j < ei[pg->imtrx]->dof[LS]; j++) {
      d_beta_dF[j] = beta * (beta0_log - beta1_log) * lsi->d_delta_dF[j] / lsi->delta_max;
    }
  }

  load_lsi(gas_phase_factor * width);
  if (af->Assemble_Jacobian)
    load_lsi_derivs();

  if (sign_gas * F >= contact_fraction * width && lsi->near) /* and add if it is "near" */
  /* This implies complete slip if gas and it is near contact
     such that pressure gradient will drive gas out */
  {
    beta = exp((beta0_log - beta1_log) * lsi->delta / lsi->delta_max + beta1_log);
    if (af->Assemble_Jacobian) {
      for (j = 0; j < ei[pg->imtrx]->dof[LS]; j++) {
        d_beta_dF[j] = beta * (beta0_log - beta1_log) * lsi->d_delta_dF[j] / lsi->delta_max;
      }
    }
  }

  /* restore LS variables  */
  if (gas_phase_factor != 1.) {
    load_lsi(width);
    if (af->Assemble_Jacobian)
      load_lsi_derivs();
  }

  return (beta);
}

#endif
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

#if 0

int
check_for_contact( Exo_DB *e,
		   double *x,
		   int ssid )
{
  int contact = FALSE;
  int k=0,i,I,ie,p;
  double F = fv->F;
  double *xl = fv->x;

  if ( F == 0.0 ) { contact = TRUE; return; }

  while ( e->ss_id[k] != ssid ) k++;

  for( i=0; i< Proc_SS_Node_Count[k]; i++)
    {
      I = e->ss_node_list[k][i];
      ie = Index_Solution(I,R_FILL, 0, 0, -2, pg->imtrx);
      if(  F*x[ie] <= 0.0 )
	{
	  double distance = 0.0;

	  for( p=0; p<pd->Num_Dim; p++)
	    {
	      distance +=  pow( xl[p] - Coor[p][I], 2.0);
	    }

	  if ( pow(distance, 0.5) <= 2.0*ls->Length_Scale )
	    {
	      contact = TRUE;
	    }
	}
    }
return (contact);
}
#endif
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void load_surface_tension(double dsigma_dx[][MDE]) /* dimensions [DIM][MDE] */

{
  int a, b, p, j;
  double dil, fact;
  /* load surface tension */
  if (mp->SurfaceTensionModel == CONSTANT) {
    return;
  } else if (mp->SurfaceTensionModel == USER) {
    usr_surface_tension(mp->u_surface_tension);
  } else if (mp->SurfaceTensionModel == DILATION) {
    belly_flop(elc->lame_mu);
    if (neg_elem_volume)
      return;
    /* calculate the surface dilation  tt..F */
    dil = 1.;
    for (a = 0; a < DIM; a++) {
      for (b = 0; b < DIM; b++) {
        GOMA_EH(GOMA_ERROR, "Need to get correct deformation gradient!!");
        dil += fv->stangent[0][a] * fv->stangent[0][b] * fv->deform_grad[a][b];
      }
    }
    mp->surface_tension = mp->u_surface_tension[0] + mp->u_surface_tension[1] * (1 - 1 / fabs(dil));

    for (p = 0; p < DIM; p++) {
      for (j = 0; j < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; j++) {
        dsigma_dx[p][j] = 0.;
      }
    }

    fact = 1.;
    if (dil < 0)
      fact = -1;
    for (a = 0; a < DIM; a++) {
      for (b = 0; b < DIM; b++) {
        for (p = 0; p < DIM; p++) {
          for (j = 0; j < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; j++) {
            dsigma_dx[p][j] +=
                mp->u_surface_tension[1] * fact / dil / dil *
                (fv->stangent[0][a] * fv->stangent[0][b] * fv->d_deform_grad_dx[a][b][p][j] +
                 fv->dstangent_dx[0][a][p][j] * fv->stangent[0][b] * fv->deform_grad[a][b] +
                 fv->stangent[0][a] * fv->dstangent_dx[0][b][p][j] * fv->deform_grad[a][b]);
          }
        }
      }
    }
  } else if (mp->SurfaceTensionModel == GIBBS_ISOTHERM) {
    double temp;
    int w;
    /* get temperature  */
    if (pd->e[pg->imtrx][TEMPERATURE]) {
      temp = fv->T;
    } else {
      temp = upd->Process_Temperature;
    }

    for (p = 0; p < DIM; p++) {
      for (j = 0; j < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; j++) {
        dsigma_dx[p][j] = 0.;
      }
    }
    mp->surface_tension =
        mp->u_surface_tension[0] * pow((mp->u_surface_tension[1] - temp) /
                                           (mp->u_surface_tension[1] - mp->reference[TEMPERATURE]),
                                       mp->u_surface_tension[2]);
    mp->d_surface_tension[TEMPERATURE] =
        mp->u_surface_tension[0] *
        pow((mp->u_surface_tension[1] - temp) /
                (mp->u_surface_tension[1] - mp->reference[TEMPERATURE]),
            mp->u_surface_tension[2] - 1.) *
        (-mp->u_surface_tension[2] / (mp->u_surface_tension[1] - mp->reference[TEMPERATURE]));

    if (mp->len_u_surface_tension > 3 && pd->Num_Species_Eqn > 0) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        if (fv->c[w] < mp->u_surface_tension[5 + 3 * w]) {
          mp->surface_tension -= mp->u_surface_tension[3 + 3 * w] *
                                 log(1.0 + fv->c[w] / mp->u_surface_tension[4 + 3 * w]);
          mp->d_surface_tension[MAX_VARIABLE_TYPES + w] =
              -mp->u_surface_tension[3 + 3 * w] / (fv->c[w] + mp->u_surface_tension[4 + 3 * w]);
        } else {
          mp->surface_tension -=
              mp->u_surface_tension[3 + 3 * w] *
              log(1.0 + mp->u_surface_tension[5 + 3 * w] / mp->u_surface_tension[4 + 3 * w]);
        }
      }
    }
  } else {
    GOMA_EH(GOMA_ERROR, "Surface tension model not defined");
  }
  return;
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

/******************************************************************************
 *
 * elec_surf_stress():
 *
 * Calculates the contributions of the Maxwell (electric) stress to the
 * traction condition.  This function assumes the fluid is incompressible
 * and that the permittivity is constant and homogenous within a material.
 *
 * Input
 * =====
 * cfunc	= Stores contributions to the residuls
 * d_cfunc	= Stores contributions to the Jacobians
 * id_side	=
 * elem_side_bc	=
 * iconnect_ptr	=
 * perm		= Permittivity
 * etm		= Equation term multiplier -- Note, the boundary equation term
 *              = multiplier for the momentum equation is not applied to this
 *              = BC (see bc_integ.c) so use this instead.
 *
 * Background:
 * ===========
 *   Tm	  = Mechanical stress tensor (e.g., -p I + 2 mu D).
 *   Te	  = Maxwell / Electrical stress tensor.
 *   I	  = Unit tensor.
 *   E	  = Electric field.
 *   K	  = dielectric constant.
 *   eo	  = Permittivity of free space.
 *   e	  = Permittivity (e = K eo).
 *   (i)  = "Inner" fluid.
 *   (o)  = "Outer" fluid.
 *   n	  = normal vector (points into (o)).
 *   t    = tangent vector.
 *   .	  = Contraction (a.b = "a dot b")
 *   s	  = Surface tension.
 *   Grad = Surface gradient: Grad(a) = grad(a) - nn.grad(a)
 *   Div  = Surface divergence: Div(a) = div(a) - nn.div(a)
 *   H	  = Mean curvature: 2H = -Div(n)
 *   V    = Voltage.
 *
 *   Electric stress: Te = eEE - 0.5 e E.E I
 *   Total stress:    T  = Tm + Te
 *
 * The full traction condition then becomes:
 *
 *   n . ( T(o) - T(i) ) = -2H s n - Grad(s)
 *
 *                        ||
 *                        ||
 *                        \/
 *
 *   n . ( Tm(o) - Tm(i) ) = -n . ( Te(o) - Te(i) ) - 2H s n - Grad(s)
 *
 * The left hand side of the latter equation is what appears in the
 * boundary terms of the momentum residuals.  Thus, in addition to the
 * curvature and surface tension gradient contributions, we need to
 * compute the components of the electric field, E = -grad(V).  Since
 * grad(V) is not continuous across the interface (though V is), this
 * function needs to be called from the elements on each side of the
 * interface.
 *
 * This formulation and implementation is valid for materials that are
 * conducting, insulating or semi-insulating.  (Semi-insulators are
 * not yet implemented in Goma.)  Note: n.Te.t = 0 for conductor /
 * dielectric and dielectric / dielectric interfaces., i.e., n.Te.t != 0
 * for semi-insulating materials only (assuming no permanant charge
 * distribution in the dielectric/dielectric case).
 *
 * Author: pknotz (PKN)
 * Date  : Aug-Sep 2001
 *
 ******************************************************************************/
void elec_surf_stress(double cfunc[MDE][DIM],
                      double d_cfunc[MDE][DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                      const int id_side,
                      struct elem_side_bc_struct *elem_side_bc,
                      const int iconnect_ptr,
                      const double perm,
                      const double etm,
                      const int bc_type) {
  int dcomp;     /* Mesh component: 0=X, 1=Y, 2=Z. */
  int mcomp;     /* Momentum equation component: 0=X, 1=Y, 2=Z. */
  int inode;     /* Element-side node index. */
  int ldof;      /* Local DOF index for the residual equation. */
  int j, p;      /* Loop counters. */
  double E[DIM]; /* The electric field. */
  int eqn = -1;  /* Equation index. */
  int id;        /* Local node number. */
  int I;         /* Global node number. */
  int var;       /* Variable index. */
  const int dim = ei[pg->imtrx]->ielem_dim;

  /* Do not check to see if this equation is active. */
  if (bc_type == ELEC_TRACTION_BC) {
    eqn = VELOCITY1;
  } else if (bc_type == ELEC_TRACTION_SOLID_BC) {
    eqn = MESH_DISPLACEMENT1;
  } else
    GOMA_EH(GOMA_ERROR, "invalid ELEC_TRACTION type");

  /* The E field. The sign doesn't matter here, but let's be anal. */
  for (p = 0; p < VIM; p++)
    E[p] = -fv->grad_V[p];

  /* First, compute the contribution to the residual equations. */
  if (af->Assemble_Residual) {
    for (inode = 0; inode < (int)elem_side_bc->num_nodes_on_side; inode++) {
      id = (int)elem_side_bc->local_elem_node_id[inode];
      I = Proc_Elem_Connect[iconnect_ptr + id];
      ldof = ei[pg->imtrx]->ln_to_first_dof[eqn][id];

      /* Dolphin[node][variable] = dof */
      /* if(Dolphin[pg->imtrx][I][VELOCITY1] > 0) */
      if (ldof >= 0) {

        /* Loop over the components of the momentum equation. */
        for (mcomp = 0; mcomp < dim; mcomp++) {

          /* Do the contractions. */
          for (p = 0; p < VIM; p++) {
            cfunc[ldof][mcomp] -= (perm * fv->snormal[p] * E[p] * E[mcomp] -
                                   0.5 * perm * E[p] * E[p] * fv->snormal[mcomp]) *
                                  bf[eqn + mcomp]->phi[ldof] * etm;
          }

        } /* for: momentum components */

      } /* if: lof >= 0 */

    } /* for: node loop */

  } /* if: Residual assembly */

  /* If we're done, bail out. */
  if (!af->Assemble_Jacobian)
    return;

  /* Now, compute the Jacobian contributions. */

  for (inode = 0; inode < elem_side_bc->num_nodes_on_side; inode++) {
    id = elem_side_bc->local_elem_node_id[inode];
    I = Proc_Elem_Connect[iconnect_ptr + id];
    ldof = ei[pg->imtrx]->ln_to_first_dof[eqn][id];

    /* if(ldof >= 0) */
    if (Dolphin[pg->imtrx][I][eqn] > 0) {

      /**********************************************************************
       * Derivatives w.r.t. the voltage.
       **********************************************************************/
      var = VOLTAGE;
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

        /* Loop over the components of the momentum equation. */
        for (mcomp = 0; mcomp < dim; mcomp++) {

          /* Do the contractions. */
          for (p = 0; p < VIM; p++) {
            d_cfunc[ldof][mcomp][var][j] -=
                (perm * fv->snormal[p] * E[p] * (-bf[var]->grad_phi[j][mcomp]) +
                 perm * fv->snormal[p] * (-bf[var]->grad_phi[j][p]) * E[mcomp] -
                 perm * (-bf[var]->grad_phi[j][p]) * E[p] * fv->snormal[mcomp]) *
                bf[eqn + mcomp]->phi[ldof] * etm;

          } /* for: Contraction */

        } /* for: Momentum components */

      } /* for: Voltage DoFs */

      /**********************************************************************
       * Derivatives w.r.t. the mesh coordinates (displacements).
       **********************************************************************/
      for (dcomp = 0; dcomp < dim; dcomp++) {
        var = MESH_DISPLACEMENT1 + dcomp;

        /* If the variable isn't defined, don't do the Jacobians. */
        if (!pd->v[pg->imtrx][var])
          continue;

        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          /* Loop over the components of the momentum equation. */
          for (mcomp = 0; mcomp < dim; mcomp++) {

            /* Do the contractions. */
            for (p = 0; p < VIM; p++) {
              d_cfunc[ldof][mcomp][var][j] -=
                  (perm * fv->dsnormal_dx[p][dcomp][j] * E[p] * E[mcomp] +
                   perm * fv->snormal[p] * (-fv->d_grad_V_dmesh[p][dcomp][j]) * E[mcomp] +
                   perm * fv->snormal[p] * E[p] * (-fv->d_grad_V_dmesh[mcomp][dcomp][j]) -
                   perm * (-fv->d_grad_V_dmesh[p][dcomp][j]) * E[p] * fv->snormal[mcomp] -
                   0.5 * perm * E[p] * E[p] * fv->dsnormal_dx[mcomp][dcomp][j]) *
                  bf[eqn + mcomp]->phi[ldof] * etm;
            }

          } /* for: momentum components */

        } /* for: j = 0,... */

      } /* for: dcomp = 0,... */

    } /* if: lof >= 0 */

  } /* for: inode = 0, ... */

  return;
}

void sheet_tension(double cfunc[MDE][DIM],
                   double d_cfunc[MDE][DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                   const int id_side, /* ID of the side of the element */
                   double sheet_tension,
                   struct elem_side_bc_struct *elem_side_bc,
                   const int iconnect_ptr,
                   double xi[DIM],
                   const Exo_DB *exo)

{
  int b, i, j, k, id, I, ldof, ldofm, i_basis;
  int eqn, var;
  int p, q;
  int *n_dof = NULL;
  int nf, el1, var_sh_tens = FALSE;
  int n_dofptr[MAX_VARIABLE_TYPES][MDE], dof_map[MDE], T_dof[MDE];

  double dY_dxi, dX_dxi, dT_dxi, detJ;
  double dY_dS, dX_dS, dT_dS = 0.0;
  double dX_dS_dmesh[MDE][DIM], dY_dS_dmesh[MDE][DIM];
  double dT_dS_dmesh[MDE][DIM], ddetJ_dmesh[MDE][DIM];
  double sign = 0.0, surf_sign = 0.0;

  double Pi[DIM][DIM];
  STRESS_DEPENDENCE_STRUCT d_Pi;
  double HL; /* hydrodynamic loading */
  double dHL_dX[DIM][MDE];
  double dHL_dv[DIM][MDE];
  double dHL_dP[MDE];
  double T[MDE];

  double dphidxi_shell[MDE], phi_shell[MDE];

  /* Determine if sh_tens is active on neighboring shell block */
  if (num_shell_blocks != 0)
    nf = num_elem_friends[ei[pg->imtrx]->ielem];
  else
    nf = 0;

  /* If so, set up assembly to include variable shell tension */
  if (nf == 1 && upd->vp[pg->imtrx][SHELL_TENSION]) {
    var_sh_tens = TRUE;
    el1 = elem_friends[ei[pg->imtrx]->ielem][0];
    n_dof = (int *)array_alloc(1, MAX_VARIABLE_TYPES, sizeof(int));
    load_neighbor_var_data(ei[pg->imtrx]->ielem, el1, n_dof, dof_map, n_dofptr, id_side, xi, exo);
    sheet_tension = fv->sh_tens; /* Overwrite input with value at this GP */
  }

  eqn = MESH_DISPLACEMENT1;

  i_basis = 1 - id_side % 2;

  dX_dxi = 0.0;
  dY_dxi = 0.0;
  dT_dxi = 0.0;
  memset(dX_dS_dmesh, 0, DIM * MDE * sizeof(double));
  memset(dY_dS_dmesh, 0, DIM * MDE * sizeof(double));
  if (var_sh_tens)
    memset(dT_dS_dmesh, 0, DIM * MDE * sizeof(double));
  memset(ddetJ_dmesh, 0, DIM * MDE * sizeof(double));
  memset(Pi, 0, DIM * DIM * sizeof(double));
  memset(&d_Pi, 0, sizeof(STRESS_DEPENDENCE_STRUCT));

  /*  if( fv->snormal[1] < 0.0 ) sheet_tension = -1.0 * fabs( sheet_tension);
      else sheet_tension = fabs(sheet_tension);    */
  if (fv->snormal[1] < 0.0)
    surf_sign = -1.0;
  else
    surf_sign = 1.0;
  /* this is a total kludge to change the sign of
   * the curvature convention so that a positive pressure
   * will always cause the surface to bulge away from the
   * the domain for boundaries in which the domain is
   * above or below the surface */

  for (i = 0; i < (int)elem_side_bc->num_nodes_on_side; i++) {
    id = (int)elem_side_bc->local_elem_node_id[i];
    I = Proc_Elem_Connect[iconnect_ptr + id];
    ldof = ei[pg->imtrx]->ln_to_first_dof[eqn][id];

    /* Grab nodal sh_tens values if needed */

    if (ldof >= 0) {
      dX_dxi += (Coor[0][I] + *esp->d[0][ldof]) * bf[eqn]->dphidxi[ldof][i_basis];
      dY_dxi += (Coor[1][I] + *esp->d[1][ldof]) * bf[eqn]->dphidxi[ldof][i_basis];
    }
  }

  if (var_sh_tens) {

    for (i = 0; i < (int)elem_side_bc->num_nodes_on_side; i++) {
      id = (int)elem_side_bc->local_elem_node_id[i];
      ldof = ei[pg->imtrx]->ln_to_first_dof[R_SHELL_TENSION][id];
      ldofm = ei[pg->imtrx]->ln_to_first_dof[R_MESH1][id];
      I = Proc_Elem_Connect[iconnect_ptr + id];
      phi_shell[i] = bf[R_MESH1]->phi[ldofm];
      dphidxi_shell[i] = bf[R_MESH1]->dphidxi[ldofm][i_basis];

      /* Grab nodal sh_tens values if needed */
      T_dof[i] = Index_Solution(I, SHELL_TENSION, 0, 0, -1, pg->imtrx);

      if (T_dof[i] != -1) {
        T[i] = x_static[T_dof[i]];
        dT_dxi += T[i] * dphidxi_shell[i];
      }
    }
  }

  sign = dX_dxi < 0.0 ? -1.0 : 1.0; /* This is the sign of dS_dxi */

  detJ = sqrt(dX_dxi * dX_dxi + dY_dxi * dY_dxi);

  if (detJ < 1.e-6) {
    GOMA_EH(GOMA_ERROR, "error in sheet_tension.");
  }

  dY_dS = sign * dY_dxi / detJ;
  dX_dS = sign * dX_dxi / detJ;
  if (var_sh_tens)
    dT_dS = sign * dT_dxi / detJ;

  if (dX_dS <= 0.0)
    printf("Detect zero dX_dS in element %d\n", ei[pg->imtrx]->ielem);

  for (i = 0; i < (int)elem_side_bc->num_nodes_on_side; i++) {
    id = (int)elem_side_bc->local_elem_node_id[i];
    ldof = ei[pg->imtrx]->ln_to_first_dof[eqn][id];

    if (ldof >= 0) {
      ddetJ_dmesh[ldof][0] = dX_dxi * bf[eqn]->dphidxi[ldof][i_basis] / detJ;
      ddetJ_dmesh[ldof][1] = dY_dxi * bf[eqn]->dphidxi[ldof][i_basis] / detJ;

      dX_dS_dmesh[ldof][0] = bf[eqn]->dphidxi[ldof][i_basis] / detJ;
      dX_dS_dmesh[ldof][0] += -(dX_dxi / detJ / detJ) * ddetJ_dmesh[ldof][0];
      dX_dS_dmesh[ldof][0] *= sign;

      dX_dS_dmesh[ldof][1] = -sign * (dX_dxi / detJ / detJ) * ddetJ_dmesh[ldof][1];

      dY_dS_dmesh[ldof][0] = -sign * (dY_dxi / detJ / detJ) * ddetJ_dmesh[ldof][0];

      dY_dS_dmesh[ldof][1] = bf[eqn]->dphidxi[ldof][i_basis] / detJ;
      dY_dS_dmesh[ldof][1] += -(dY_dxi / detJ / detJ) * ddetJ_dmesh[ldof][1];
      dY_dS_dmesh[ldof][1] *= sign;
    }

    if (var_sh_tens) {
      dT_dS_dmesh[ldof][0] = -sign * (dT_dxi / detJ / detJ) * ddetJ_dmesh[ldof][0];
      dT_dS_dmesh[ldof][1] = -sign * (dT_dxi / detJ / detJ) * ddetJ_dmesh[ldof][1];
    }
  }

  fluid_stress(Pi, &d_Pi);

  HL = 0.0;
  memset(dHL_dv, 0, sizeof(double) * DIM * MDE);
  memset(dHL_dX, 0, sizeof(double) * DIM * MDE);
  memset(dHL_dP, 0, sizeof(double) * MDE);

  for (p = 0; p < ei[pg->imtrx]->ielem_dim; p++) {
    for (q = 0; q < ei[pg->imtrx]->ielem_dim; q++) {
      HL += surf_sign * fv->snormal[p] * fv->snormal[q] * Pi[q][p];

      for (j = 0; j < ei[pg->imtrx]->dof[VELOCITY1]; j++) {
        for (b = 0; b < ei[pg->imtrx]->ielem_dim; b++) {

          dHL_dv[b][j] += surf_sign * fv->snormal[p] * fv->snormal[q] * d_Pi.v[q][p][b][j];
        }
      }

      for (j = 0; j < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; j++) {
        for (b = 0; b < ei[pg->imtrx]->ielem_dim; b++) {
          dHL_dX[b][j] += surf_sign * fv->dsnormal_dx[p][b][j] * fv->snormal[q] * Pi[q][p];
          dHL_dX[b][j] += surf_sign * fv->snormal[p] * fv->dsnormal_dx[q][b][j] * Pi[q][p];
          dHL_dX[b][j] += surf_sign * fv->snormal[p] * fv->snormal[q] * d_Pi.X[q][p][b][j];
        }
      }

      for (j = 0; j < ei[pg->imtrx]->dof[PRESSURE]; j++) {
        dHL_dP[j] += surf_sign * fv->snormal[p] * fv->snormal[q] * d_Pi.P[q][p][j];
      }
    }
  }

  if (pd->v[pg->imtrx][eqn]) {
    for (i = 0; i < (int)elem_side_bc->num_nodes_on_side; i++) {
      id = (int)elem_side_bc->local_elem_node_id[i];
      ldof = ei[pg->imtrx]->ln_to_first_dof[eqn][id];

      if (ldof >= 0) {
        cfunc[ldof][0] = -sheet_tension * sign * bf[eqn]->dphidxi[ldof][i_basis] * dY_dS / detJ;
        if (var_sh_tens) {
          cfunc[ldof][0] -= bf[eqn]->phi[ldof] * dT_dS * dY_dS;
        }

        cfunc[ldof][0] -= bf[eqn]->phi[ldof] * (HL)*dX_dS;
      }
    }
  }

  if (af->Assemble_Jacobian) {

    for (i = 0; i < (int)elem_side_bc->num_nodes_on_side; i++) {
      id = (int)elem_side_bc->local_elem_node_id[i];
      ldof = ei[pg->imtrx]->ln_to_first_dof[eqn][id];

      if (ldof >= 0) {

        for (b = 0; b < ei[pg->imtrx]->ielem_dim; b++) {
          var = MESH_DISPLACEMENT1 + b;

          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_cfunc[ldof][0][var][j] =
                  -(bf[eqn]->dphidxi[ldof][i_basis] * dY_dS) * ddetJ_dmesh[j][b] / detJ / detJ;
              d_cfunc[ldof][0][var][j] +=
                  bf[eqn]->dphidxi[ldof][i_basis] * dY_dS_dmesh[j][b] / detJ;
              d_cfunc[ldof][0][var][j] *= -sheet_tension * sign;

              d_cfunc[ldof][0][var][j] -= bf[eqn]->phi[ldof] * (HL)*dX_dS_dmesh[j][b];
              d_cfunc[ldof][0][var][j] -= bf[eqn]->phi[ldof] * (dHL_dX[b][j]) * dX_dS;

              /* if (var_sh_tens)
                 {
                   d_cfunc[ldof][0][var][j] -= (dT_dS_dmesh[j][b]
                   * sign * dY_dS * bf[eqn]->dphidxi[ldof][i_basis] / detJ);
                   d_cfunc[ldof][0][var][j] -= (dY_dS_dmesh[j][b]
                   * sign * dT_dS * bf[eqn]->dphidxi[ldof][i_basis] / detJ);
                                     }*/
              if (var_sh_tens) {
                d_cfunc[ldof][0][var][j] -= bf[eqn]->phi[ldof] * (dT_dS_dmesh[j][b] * dY_dS);
                d_cfunc[ldof][0][var][j] -= bf[eqn]->phi[ldof] * (dY_dS_dmesh[j][b] * dT_dS);
                /* d_cfunc[ldof][0][var][j] += bf[eqn]->phi[ldof]*(dT_dS * dY_dS/detJ/detJ
                 * )*ddetJ_dmesh[j][b];*/
              }
            }
          }
        }

        var = PRESSURE;

        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_cfunc[ldof][0][var][j] = -bf[eqn]->phi[ldof] * (dHL_dP[j]) * dX_dS;
          }
        }

        var = VELOCITY1;
        if (pd->v[pg->imtrx][var]) {

          for (b = 0; b < ei[pg->imtrx]->ielem_dim; b++) {
            var = VELOCITY1 + b;

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_cfunc[ldof][0][var][j] = -bf[eqn]->phi[ldof] * (dHL_dv[b][j]) * dX_dS;
            }
          }
        }

        var = POLYMER_STRESS11;
        if (pd->v[pg->imtrx][var]) {
          GOMA_EH(GOMA_ERROR,
                  "TENSION_SHEET BC is not instrumented for polymeric fluid interactions\n");
        }

        if (var_sh_tens) {
          var = SHELL_TENSION;

          for (k = 0; k < (int)elem_side_bc->num_nodes_on_side; k++) {
            id = (int)elem_side_bc->local_elem_node_id[k];
            I = Proc_Elem_Connect[iconnect_ptr + id];
            j = ei[pg->imtrx]->ln_to_first_dof[R_SHELL_TENSION][id];

            /*  d_cfunc[ldof][0][var][j] = -bf[var]->phi[j] * dY_dS* sign *
               bf[eqn]->dphidxi[ldof][i_basis] / detJ; d_cfunc[ldof][0][var][j]  -=
               bf[eqn]->phi[ldof] *  sign*bf[var]->dphidxi[j][0] * dY_dS /detJ ;   */

            if (T_dof[k] != -1) {
              d_cfunc[ldof][0][var][j] =
                  -phi_shell[k] * sign * bf[eqn]->dphidxi[ldof][i_basis] * dY_dS / detJ;
              d_cfunc[ldof][0][var][j] -=
                  bf[eqn]->phi[ldof] * sign * dphidxi_shell[k] * dY_dS / detJ;
            }
          }
        }
      }
    }
  }
  return;
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void fn_dot_T(double cfunc[MDE][DIM],
              double d_cfunc[MDE][DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
              const int id_side,    /* ID of the side of the element             */
              const double sigma,   /* surface tension                           */
              const double pb[DIM], /* applied pressure                          */
              struct elem_side_bc_struct *elem_side_bc,
              const int iconnect_ptr,
              double dsigma_dx[DIM][MDE])

/*******************************************************************************
 *
 *  Function which calculates the calculates the capillary free surface stress
 *  balance:
 *        2H*sigma*n + pb[0]*n + pb[1]*t1 + pb[2]*t2 (not enabled yet) = n.T
 *  This vector condition is to be added on component wise to the momentum equations.
 *  pb is the applied pressure as in a vacuum or a forcing function
 *
 *******************************************************************************/

{
  /*    TAB certifies that this function conforms to the exo/patran side numbering convention
   * 11/10/98. */
  int j, i, id, var, a, eqn, I, ldof, w, dim;
  int p, q, jvar; /* Degree of freedom counter                 */
                  /***************************** EXECUTION BEGINS ******************************/
  /* Based on current element id_side, choose the correct curvature sign */

  /* EDW: For 3D of 2D LSA; all three components are needed */
  dim = ei[pg->imtrx]->ielem_dim;
  if (Linear_Stability == LSA_3D_OF_2D || Linear_Stability == LSA_3D_OF_2D_SAVE)
    dim = VIM;

  eqn = VELOCITY1;

  if (af->Assemble_Jacobian) {
    for (i = 0; i < (int)elem_side_bc->num_nodes_on_side; i++) {
      id = (int)elem_side_bc->local_elem_node_id[i];
      I = Proc_Elem_Connect[iconnect_ptr + id];
      ldof = ei[pg->imtrx]->ln_to_first_dof[eqn][id];

      /* Calculate the residual contribution from surface gradient of basis function */
      if (Dolphin[pg->imtrx][I][VELOCITY1] > 0) {

        /* Calculate the residual contribution from surface gradient of basis function */
        /*
         *  Evaluate sensitivity to displacements d()/dx
         */
        for (jvar = 0; jvar < ei[pg->imtrx]->ielem_dim; jvar++) {
          var = MESH_DISPLACEMENT1 + jvar;
          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              for (a = 0; a < dim; a++) {
                d_cfunc[ldof][a][var][j] -=
                    (pb[0] * fv->dsnormal_dx[a][jvar][j]) * bf[VELOCITY1 + a]->phi[ldof];
                d_cfunc[ldof][a][var][j] -=
                    (pb[1] * fv->dstangent_dx[0][a][jvar][j]) * bf[VELOCITY1 + a]->phi[ldof];
                if (ei[pg->imtrx]->ielem_dim == DIM)
                  d_cfunc[ldof][a][var][j] -=
                      (pb[2] * fv->dstangent_dx[0][a][jvar][j]) * bf[VELOCITY1 + a]->phi[ldof];
                for (p = 0; p < VIM; p++) {

                  d_cfunc[ldof][a][var][j] -=
                      sigma * mp->surface_tension *
                      bf[VELOCITY1 + a]->d_grad_phi_e_dmesh[ldof][a][p][p][jvar][j];

                  for (q = 0; q < VIM; q++) {
                    d_cfunc[ldof][a][var][j] +=
                        sigma * mp->surface_tension *
                        (bf[VELOCITY1 + a]->d_grad_phi_e_dmesh[ldof][a][q][p][jvar][j] *
                             fv->snormal[p] * fv->snormal[q] +
                         bf[VELOCITY1 + a]->grad_phi_e[ldof][a][q][p] *
                             fv->dsnormal_dx[p][jvar][j] * fv->snormal[q] +
                         bf[VELOCITY1 + a]->grad_phi_e[ldof][a][q][p] * fv->snormal[p] *
                             fv->dsnormal_dx[q][jvar][j]);
                  }
                }
                if (mp->SurfaceTensionModel == DILATION) {
                  for (p = 0; p < VIM; p++) {
                    d_cfunc[ldof][a][var][j] -=
                        sigma * dsigma_dx[jvar][j] * bf[VELOCITY1 + a]->grad_phi_e[ldof][a][p][p];
                    for (q = 0; q < VIM; q++) {
                      d_cfunc[ldof][a][var][j] += sigma * dsigma_dx[jvar][j] *
                                                  bf[VELOCITY1 + a]->grad_phi_e[ldof][a][q][p] *
                                                  fv->snormal[p] * fv->snormal[q];
                    }
                  }
                }
              }
            }
          }
        }
        /*
         *  Evaluate sensitivity to conc and temp (if surface tension is
         *    variable
         */
        var = TEMPERATURE;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            for (a = 0; a < dim; a++) {
              for (p = 0; p < VIM; p++) {
                d_cfunc[ldof][a][var][j] -= sigma * mp->d_surface_tension[TEMPERATURE] *
                                            bf[VELOCITY1 + a]->grad_phi_e[ldof][a][p][p] *
                                            bf[var]->phi[j];
                for (q = 0; q < VIM; q++) {
                  d_cfunc[ldof][a][var][j] += sigma * mp->d_surface_tension[TEMPERATURE] *
                                              bf[VELOCITY1 + a]->grad_phi_e[ldof][a][q][p] *
                                              fv->snormal[p] * fv->snormal[q] * bf[var]->phi[j];
                }
              }
            }
          }
        }
        var = MASS_FRACTION;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            for (w = 0; w < pd->Num_Species_Eqn; w++) {
              for (a = 0; a < dim; a++) {
                for (p = 0; p < VIM; p++) {
                  d_cfunc[ldof][a][var][j] -=
                      sigma * mp->d_surface_tension[MAX_VARIABLE_TYPES + w] *
                      bf[VELOCITY1 + a]->grad_phi_e[ldof][a][p][p] * bf[var]->phi[j];
                  for (q = 0; q < VIM; q++) {
                    d_cfunc[ldof][a][var][j] += sigma *
                                                mp->d_surface_tension[MAX_VARIABLE_TYPES + w] *
                                                bf[VELOCITY1 + a]->grad_phi_e[ldof][a][q][p] *
                                                fv->snormal[p] * fv->snormal[q] * bf[var]->phi[j];
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  eqn = VELOCITY1;
  if (pd->v[pg->imtrx][eqn]) {
    for (i = 0; i < (int)elem_side_bc->num_nodes_on_side; i++) {
      id = (int)elem_side_bc->local_elem_node_id[i];
      I = Proc_Elem_Connect[iconnect_ptr + id];
      ldof = ei[pg->imtrx]->ln_to_first_dof[eqn][id];

      /* Calculate the residual contribution from surface gradient of basis function */
      if (ldof >= 0) {
        /*
         * Equate the boundary stress in the fluid mechanics momentum equations to the
         * external pressure plus capillarity plus the surface tangent at ends of the interface
         *
         * evaluate capillarity as = sigma (I - nn):grad(e_a phi)  this should apply in 2D or 3D
         */

        for (a = 0; a < dim; a++) {
          cfunc[ldof][a] -= pb[0] * fv->snormal[a] * bf[VELOCITY1 + a]->phi[ldof];
          cfunc[ldof][a] -= pb[1] * fv->stangent[0][a] * bf[VELOCITY1 + a]->phi[ldof];
          if (ei[pg->imtrx]->ielem_dim == DIM)
            cfunc[ldof][a] -= pb[2] * fv->stangent[1][a] * bf[VELOCITY1 + a]->phi[ldof];
          for (p = 0; p < VIM; p++) {
            cfunc[ldof][a] -=
                sigma * mp->surface_tension * bf[VELOCITY1 + a]->grad_phi_e[ldof][a][p][p];
            for (q = 0; q < VIM; q++) {
              cfunc[ldof][a] += sigma * mp->surface_tension *
                                bf[VELOCITY1 + a]->grad_phi_e[ldof][a][q][p] * fv->snormal[p] *
                                fv->snormal[q];
            }
          }
        }
      }
    }
  } /* if velocity variable is defined.   */
  else {
    GOMA_EH(GOMA_ERROR, "Bad coord system in capillary or liquid momentum equations not defined");
  }

  return;

} /* END of routine fn_dot_T                                                */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void apply_repulsion(double cfunc[MDE][DIM],
                     double d_cfunc[MDE][DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                     const double pr, /* coefficient for repulsion force to ensure
                                       * no penetration of the solid boundary by
                                       * the free surface                        */
                     const double ap, /* a coefficient in plane equation */
                     const double bp, /* b coefficient in plane equation */
                     const double cp, /* c coefficient in plane equation */
                     const double dp, /* d coefficient in plane equation */
                     struct elem_side_bc_struct *elem_side_bc,
                     const int iconnect_ptr)

/******************************************************************************
 *
 *  Function which calculates the surface repulsion from a
 *  general plane boundary
 *  pr/dist**2 is the force applied to repulse the free surface from the solid boundaries
 *  dist is defined as the distance between the free surface and the solid wall.
 *
 ******************************************************************************/

{

  /* Local variables */

  int j, i, id, var, a, eqn, I, ldof;
  int jvar; /* Degree of freedom counter                 */

  double d_dist[DIM]; /* distance from surface to wall             */
  double dist = 0;    /* squared distance from surface to wall     */
  double repexp = 2.; /* exponent of disance in repulsion term     */
  double denom, factor;
  /***************************** EXECUTION BEGINS ******************************/
  /* if pr is 0. => we don't want free surface/wall repulsion and just
     ensure that dist and dist2 are nonzero so nothing bad happens */
  if (pr == 0)
    return;

  d_dist[0] = 0;
  d_dist[1] = 0;
  d_dist[2] = 0;

  denom = sqrt(ap * ap + bp * bp + cp * cp);
  factor = ap * fv->x[0] + bp * fv->x[1] + cp * fv->x[2] + dp;
  dist = fabs(factor) / denom;
  d_dist[0] = SGN(factor) * ap / denom;
  d_dist[1] = SGN(factor) * bp / denom;
  d_dist[2] = SGN(factor) * cp / denom;

  eqn = VELOCITY1;

  if (af->Assemble_Jacobian) {
    for (i = 0; i < (int)elem_side_bc->num_nodes_on_side; i++) {
      id = (int)elem_side_bc->local_elem_node_id[i];
      I = Proc_Elem_Connect[iconnect_ptr + id];
      ldof = ei[pg->imtrx]->ln_to_dof[eqn][id];
      if (Dolphin[pg->imtrx][I][VELOCITY1] > 0) {

        /*
         *  Evaluate sensitivity to displacements d()/dx
         */
        for (jvar = 0; jvar < VIM; jvar++) {
          var = MESH_DISPLACEMENT1 + jvar;
          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              for (a = 0; a < VIM; a++) {

                d_cfunc[ldof][a][var][j] -= (pr / pow(dist, repexp) * fv->dsnormal_dx[a][jvar][j] -
                                             repexp * pr / pow(dist, repexp + 1) * fv->snormal[a] *
                                                 d_dist[jvar] * bf[var]->phi[j]) *
                                            bf[eqn]->phi[ldof];
              }
            }
          }
        }
      }
    }
  }
  /*fprintf(stderr,"repulse %g %g %g %g
   * %g\n",dist,pr,pr/pow(dist,repexp),fv->snormal[0],fv->snormal[1]);  */

  eqn = VELOCITY1;
  for (i = 0; i < (int)elem_side_bc->num_nodes_on_side; i++) {
    id = (int)elem_side_bc->local_elem_node_id[i];
    I = Proc_Elem_Connect[iconnect_ptr + id];
    ldof = ei[pg->imtrx]->ln_to_dof[eqn][id];
    if (Dolphin[pg->imtrx][I][VELOCITY1] > 0) {

      for (a = 0; a < VIM; a++) {
        cfunc[ldof][a] -= pr / pow(dist, repexp) * fv->snormal[a] * bf[eqn]->phi[ldof];
      }
    }

  } /* end of for (i = 0; i < (int) elem_side_bc->num_nodes_on_side */

} /* END of routine apply_repulsion                                          */
/*****************************************************************************/

void apply_repulsion_roll(double cfunc[MDE][DIM],
                          double d_cfunc[MDE][DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                          double x[],                /* solution vector*/
                          const double roll_rad,     /* roll radius */
                          const double origin[3],    /* roll axis origin (x,y,z) */
                          const double dir_angle[3], /* axis direction angles */
                          const double omega,        /* roll rotation rate  */
                          const double P_rep,        /* repulsion coefficient */
                          const double hscale,       /* repulsion length scale */
                          const double repexp,       /* repulsive force exponent */
                          const double gas_visc,     /* inverse slip coefficient  */
                          const double exp_scale,    /* DCL exculsion zone scale   */
                          const int dcl_node,        /* DCL NS id  */
                          struct elem_side_bc_struct *elem_side_bc,
                          const int iconnect_ptr)

/******************************************************************************
 *
 *  Function which calculates the surface repulsion from a rotating roll
 *  with a slip velocity condition
 *
 ******************************************************************************/

{

  /* Local variables */

  int dim, jvar; /* Degree of freedom counter */

  double d_dist[DIM]; /* distance derivatives  */
  double dist = 1e12; /* squared distance from surface to wall     */
  double factor;
  double coord[3] = {0, 0, 0};
  double axis_pt[3], rad_dir[3], v_dir[3], v_roll[3], t, R;
  double force = 0.0, d_force = 0.0, inv_slip = 0.0, d_inv_slip = 0.0;
  double sheara = 0., d_sheara = 0., shearb = 0., d_shearb = 0.;
  double t_veloc[2] = {0, 0}, dt_veloc_dx[2][MAX_PDIM][MAX_PDIM][MDE];
  double n_veloc = 0, dn_veloc_dx[MAX_PDIM][MAX_PDIM][MDE];
  double dsheara_dx[MAX_PDIM][MAX_PDIM][MDE], dshearb_dx[MAX_PDIM][MAX_PDIM][MDE];
  double dcl_dist, mod_factor = 1., point[3] = {0, 0, 0};
  int nsp, k;

  int j, i, id, var, a, eqn, I, ldof;

  /***************************** EXECUTION BEGINS ******************************/
  /* if pr is 0. => we don't want free surface/wall repulsion and just
     ensure that dist and dist2 are nonzero so nothing bad happens */
  if (P_rep == 0 && gas_visc == 0)
    return;

  eqn = VELOCITY1;
  if (af->Assemble_LSA_Mass_Matrix)
    return;

  /*  initialize variables */
  dim = pd->Num_Dim;

  /* calculate distance from free surface to solid surface for repulsion calculations */

  coord[0] = fv->x[0];
  coord[1] = fv->x[1];
  if (dim == 3) {
    coord[2] = fv->x[2];
  } else {
    coord[2] = 0.0;
  }

  /*  find intersection of axis with normal plane - i.e., locate point on
   *          axis that intersects plane normal to axis that contains local point. */

  factor = SQUARE(dir_angle[0]) + SQUARE(dir_angle[1]) + SQUARE(dir_angle[2]);
  t = (dir_angle[0] * (coord[0] - origin[0]) + dir_angle[1] * (coord[1] - origin[1]) +
       dir_angle[2] * (coord[2] - origin[2])) /
      factor;
  axis_pt[0] = origin[0] + dir_angle[0] * t;
  axis_pt[1] = origin[1] + dir_angle[1] * t;
  axis_pt[2] = origin[2] + dir_angle[2] * t;

  /*  compute radius and radial direction */

  R = sqrt(SQUARE(coord[0] - axis_pt[0]) + SQUARE(coord[1] - axis_pt[1]) +
           SQUARE(coord[2] - axis_pt[2]));
  rad_dir[0] = (coord[0] - axis_pt[0]) / R;
  rad_dir[1] = (coord[1] - axis_pt[1]) / R;
  rad_dir[2] = (coord[2] - axis_pt[2]) / R;
  dist = R - roll_rad;
  d_dist[0] = rad_dir[0] * (1. - SQUARE(dir_angle[0]) / factor) +
              rad_dir[1] * (-dir_angle[1] * dir_angle[0] / factor) +
              rad_dir[2] * (-dir_angle[2] * dir_angle[0] / factor);
  d_dist[1] = rad_dir[1] * (1. - SQUARE(dir_angle[1]) / factor) +
              rad_dir[0] * (-dir_angle[0] * dir_angle[1] / factor) +
              rad_dir[2] * (-dir_angle[2] * dir_angle[1] / factor);
  d_dist[2] = rad_dir[2] * (1. - SQUARE(dir_angle[2]) / factor) +
              rad_dir[0] * (-dir_angle[0] * dir_angle[2] / factor) +
              rad_dir[1] * (-dir_angle[1] * dir_angle[2] / factor);

  /* compute velocity direction as perpendicular to both axis and radial
   *         direction.  Positive direction is determined by right hand rule */

  v_dir[0] = dir_angle[1] * rad_dir[2] - dir_angle[2] * rad_dir[1];
  v_dir[1] = dir_angle[2] * rad_dir[0] - dir_angle[0] * rad_dir[2];
  v_dir[2] = dir_angle[0] * rad_dir[1] - dir_angle[1] * rad_dir[0];

  v_roll[0] = omega * roll_rad * v_dir[0];
  v_roll[1] = omega * roll_rad * v_dir[1];
  v_roll[2] = omega * roll_rad * v_dir[2];

  /* DCL exclusion zone  */
  if (dcl_node != -1) {
    nsp = match_nsid(dcl_node);
    k = Proc_NS_List[Proc_NS_Pointers[nsp]];
    for (j = 0; j < Proc_NS_Count[nsp]; j++) {
      k = Proc_NS_List[Proc_NS_Pointers[nsp] + j];
      i = Index_Solution(k, MESH_DISPLACEMENT1, 0, 0, -1, pg->imtrx);
      GOMA_EH(i, "Could not resolve index_solution.");
      for (a = 0; a < dim; a++) {
        point[a] = Coor[a][k] + x[i + a];
      }
    }
    dcl_dist = sqrt(SQUARE(coord[0] - point[0]) + SQUARE(coord[1] - point[1]) +
                    SQUARE(coord[2] - point[2]));
  } else {
    dcl_dist = 0.;
    GOMA_WH(-1, "No DCL node for CAP_REPULSE_ROLL....\n");
  }
  /*  modifying function for DCL  */
  if (dcl_node != -1) {
    mod_factor = 1. - exp(-dcl_dist / exp_scale);
  } else {
    mod_factor = 1.;
  }

/*  repulsion function  */
#if 0
           force = -P_rep*mod_factor/pow(dist/hscale, repexp); 
           d_force = P_rep*mod_factor*repexp/pow(dist/hscale, repexp+1)/hscale;
#else
  force = -P_rep * mod_factor * (pow(hscale / dist, repexp) - pow(hscale / dist, repexp / 2));
  d_force = P_rep * mod_factor * (repexp / dist) *
            (pow(hscale / dist, repexp) - 0.5 * pow(hscale / dist, repexp / 2));
#endif
  /*  slip velocity function function  */
  inv_slip = -gas_visc * mod_factor / pow(dist / hscale, repexp);
  d_inv_slip = gas_visc * mod_factor * repexp / pow(dist / hscale, repexp + 1) / hscale;
  t_veloc[0] = t_veloc[1] = 0.;
  for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {
    t_veloc[0] += fv->stangent[0][a] * (fv->v[a] - v_roll[a]);
    n_veloc += fv->snormal[a] * (fv->v[a] - v_roll[a]);
    sheara += 0.5 * dist * d_force * d_dist[a] * fv->stangent[0][a];
    d_sheara +=
        0.5 * fv->stangent[0][a] * (d_force * d_dist[a] - d_force * (repexp + 1.) * SQUARE(hscale));

    for (jvar = 0; jvar < ei[pg->imtrx]->ielem_dim; jvar++) {
      var = MESH_DISPLACEMENT1 + jvar;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          dt_veloc_dx[0][a][jvar][j] = fv->dstangent_dx[0][a][jvar][j] * (fv->v[a] - v_roll[a]);
          dn_veloc_dx[a][jvar][j] = fv->dsnormal_dx[a][jvar][j] * (fv->v[a] - v_roll[a]);
          dn_veloc_dx[a][jvar][j] = 0.;
          dsheara_dx[a][jvar][j] =
              0.5 * dist * d_force * d_dist[a] * fv->dstangent_dx[0][a][jvar][j];
        }
      }
    }
  }
  if (dim == 3) {
    for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {
      t_veloc[1] += fv->stangent[1][a] * (fv->v[a] - v_roll[a]);
      shearb += 0.5 * dist * d_force * d_dist[a] * fv->stangent[1][a];
      d_shearb += 0.5 * fv->stangent[1][a] *
                  (d_force * d_dist[a] - d_force * (repexp + 1.) * SQUARE(hscale));

      for (jvar = 0; jvar < ei[pg->imtrx]->ielem_dim; jvar++) {
        var = MESH_DISPLACEMENT1 + jvar;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            dt_veloc_dx[1][a][jvar][j] = fv->dstangent_dx[1][a][jvar][j] * (fv->v[a] - v_roll[a]);
            dshearb_dx[a][jvar][j] =
                0.5 * dist * d_force * d_dist[a] * fv->dstangent_dx[1][a][jvar][j];
          }
        }
      }
    }
  }

  if (af->Assemble_Jacobian) {
    for (i = 0; i < (int)elem_side_bc->num_nodes_on_side; i++) {
      id = (int)elem_side_bc->local_elem_node_id[i];
      I = Proc_Elem_Connect[iconnect_ptr + id];
      ldof = ei[pg->imtrx]->ln_to_dof[eqn][id];
      if (Dolphin[pg->imtrx][I][VELOCITY1] > 0) {

        /*
         *  Evaluate sensitivity to displacements d()/dx
         */
        for (jvar = 0; jvar < ei[pg->imtrx]->ielem_dim; jvar++) {
          var = MESH_DISPLACEMENT1 + jvar;
          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {

                d_cfunc[ldof][a][var][j] += (n_veloc * force * fv->dsnormal_dx[a][jvar][j] +
                                             (n_veloc * d_force + force * dn_veloc_dx[a][jvar][j]) *
                                                 fv->snormal[a] * d_dist[jvar] * bf[var]->phi[j]) *
                                            bf[eqn]->phi[ldof];
                d_cfunc[ldof][a][var][j] +=
                    ((t_veloc[0] * inv_slip + sheara) * fv->dstangent_dx[0][a][jvar][j] +
                     (t_veloc[0] * d_inv_slip + d_sheara) * fv->stangent[0][a] * d_dist[jvar] *
                         bf[var]->phi[j] +
                     fv->stangent[0][a] *
                         (inv_slip * dt_veloc_dx[0][a][jvar][j] + dsheara_dx[a][jvar][j])) *
                    bf[eqn]->phi[ldof];
                if (dim == 3)
                  d_cfunc[ldof][a][var][j] +=
                      ((t_veloc[1] * inv_slip + shearb) * fv->dstangent_dx[1][a][jvar][j] +
                       (t_veloc[1] * d_inv_slip + d_shearb) * fv->stangent[1][a] * d_dist[jvar] *
                           bf[var]->phi[j] +
                       fv->stangent[1][a] *
                           (inv_slip * dt_veloc_dx[1][a][jvar][j] + dshearb_dx[a][jvar][j])) *
                      bf[eqn]->phi[ldof];
              }
            }
          }
        }
        /*
         *  Evaluate sensitivity to velocities
         */
        for (jvar = 0; jvar < ei[pg->imtrx]->ielem_dim; jvar++) {
          var = VELOCITY1 + jvar;
          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {
                d_cfunc[ldof][a][var][j] += fv->snormal[jvar] * bf[var]->phi[j] * force *
                                            fv->snormal[a] * bf[eqn]->phi[ldof];
                d_cfunc[ldof][a][var][j] += fv->stangent[0][jvar] * bf[var]->phi[j] * inv_slip *
                                            fv->stangent[0][a] * bf[eqn]->phi[ldof];

                if (dim == 3)
                  d_cfunc[ldof][a][var][j] += fv->stangent[1][jvar] * bf[var]->phi[j] * inv_slip *
                                              fv->stangent[1][a] * bf[eqn]->phi[ldof];
              }
            }
          }
        }
      }
    }
  }

  eqn = VELOCITY1;
  for (i = 0; i < (int)elem_side_bc->num_nodes_on_side; i++) {
    id = (int)elem_side_bc->local_elem_node_id[i];
    I = Proc_Elem_Connect[iconnect_ptr + id];
    ldof = ei[pg->imtrx]->ln_to_dof[eqn][id];
    if (Dolphin[pg->imtrx][I][VELOCITY1] > 0) {

      for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {
        cfunc[ldof][a] += force * n_veloc * fv->snormal[a] * bf[eqn]->phi[ldof];
        cfunc[ldof][a] +=
            (inv_slip * t_veloc[0] + sheara) * fv->stangent[0][a] * bf[eqn]->phi[ldof];
      }
      if (dim == 3) {
        for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {
          cfunc[ldof][a] +=
              (inv_slip * t_veloc[1] + shearb) * fv->stangent[1][a] * bf[eqn]->phi[ldof];
        }
      }
    }

  } /* end of for (i = 0; i < (int) elem_side_bc->num_nodes_on_side */

} /* END of routine apply_repulsion_roll                                          */
/*****************************************************************************/

void apply_repulsion_user(double cfunc[MDE][DIM],
                          double d_cfunc[MDE][DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                          const double roll_rad,     /* roll radius */
                          const double origin[3],    /* roll axis origin (x,y,z) */
                          const double dir_angle[3], /* axis direction angles */
                          const double hscale,       /* repulsion length scale */
                          const double repexp,       /* repulsive force exponent */
                          const double P_rep,        /* repulsion coefficient */
                          const double betainv,      /* inverse slip coefficient  */
                          const double omega,        /* roll rotation rate  */
                          struct elem_side_bc_struct *elem_side_bc,
                          const int iconnect_ptr)

/******************************************************************************
 *
 *  Function which calculates the surface repulsion from a rotating roll
 *  with a slip velocity condition
 *
 ******************************************************************************/

{

  /* Local variables */

  int dim, jvar; /* Degree of freedom counter */

  double d_dist[DIM]; /* distance derivatives  */
  double dist = 1e12; /* squared distance from surface to wall     */
  double force = 0.0, d_force = 0.0, inv_slip = 0.0, d_inv_slip = 0.0;
  double t_veloc[2], dt_veloc_dx[2][MAX_PDIM][MAX_PDIM][MDE];

  int j, i, id, var, a, eqn, I, ldof;
  double time = 0.;

  /***************************** EXECUTION BEGINS ******************************/
  /* if pr is 0. => we don't want free surface/wall repulsion and just
     ensure that dist and dist2 are nonzero so nothing bad happens */
  if (P_rep == 0 && betainv == 0)
    return;

  eqn = VELOCITY1;
  if (af->Assemble_LSA_Mass_Matrix)
    return;

  /*  initialize variables */
  dim = pd->Num_Dim;

  /* calculate distance from free surface to solid surface for repulsion calculations

        coord[0] = fv->x[0];
        coord[1] = fv->x[1];
        if( dim == 3)
          { coord[2] = fv->x[2];}
        else
          { coord[2] = 0.0;}*/

  dist = fnc(fv->x[0], fv->x[1], fv->x[2], BC_Types[2].u_BC, time);
  d_dist[0] = dfncd1(fv->x[0], fv->x[1], fv->x[2], BC_Types[2].u_BC, time);
  d_dist[1] = dfncd2(fv->x[0], fv->x[1], fv->x[2], BC_Types[2].u_BC, time);

  if (dim == 3)
    d_dist[2] = dfncd3(fv->x[0], fv->x[1], fv->x[2], BC_Types[2].u_BC, time);

  /*  repulsion function  */
  force = -P_rep / pow(dist / hscale, repexp);
  d_force = P_rep * repexp / pow(dist / hscale, repexp + 1) / hscale;
  /*  slip velocity function function  */
  inv_slip = -betainv / pow(dist / hscale, repexp);
  d_inv_slip = betainv * repexp / pow(dist / hscale, repexp + 1) / hscale;
  t_veloc[0] = t_veloc[1] = 0.;
  for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {
    t_veloc[0] += fv->stangent[0][a] * (fv->v[a]);
    for (jvar = 0; jvar < ei[pg->imtrx]->ielem_dim; jvar++) {
      var = MESH_DISPLACEMENT1 + jvar;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          dt_veloc_dx[0][a][jvar][j] = fv->dstangent_dx[0][a][jvar][j] * (fv->v[a]);
        }
      }
    }
  }
  if (dim == 3) {
    for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {
      t_veloc[1] += fv->stangent[1][a] * (fv->v[a]);
      for (jvar = 0; jvar < ei[pg->imtrx]->ielem_dim; jvar++) {
        var = MESH_DISPLACEMENT1 + jvar;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            dt_veloc_dx[1][a][jvar][j] = fv->dstangent_dx[1][a][jvar][j] * (fv->v[a]);
          }
        }
      }
    }
  }

  if (af->Assemble_Jacobian) {
    for (i = 0; i < (int)elem_side_bc->num_nodes_on_side; i++) {
      id = (int)elem_side_bc->local_elem_node_id[i];
      I = Proc_Elem_Connect[iconnect_ptr + id];
      ldof = ei[pg->imtrx]->ln_to_dof[eqn][id];
      if (Dolphin[pg->imtrx][I][VELOCITY1] > 0) {

        /*
         *  Evaluate sensitivity to displacements d()/dx
         */
        for (jvar = 0; jvar < ei[pg->imtrx]->ielem_dim; jvar++) {
          var = MESH_DISPLACEMENT1 + jvar;
          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {

                d_cfunc[ldof][a][var][j] +=
                    (force * fv->dsnormal_dx[a][jvar][j] +
                     d_force * fv->snormal[a] * d_dist[jvar] * bf[var]->phi[j]) *
                    bf[eqn]->phi[ldof];
                d_cfunc[ldof][a][var][j] +=
                    (t_veloc[0] * inv_slip * fv->dstangent_dx[0][a][jvar][j] +
                     t_veloc[0] * d_inv_slip * fv->stangent[0][a] * d_dist[jvar] * bf[var]->phi[j] +
                     inv_slip * fv->stangent[0][a] * dt_veloc_dx[0][a][jvar][j]) *
                    bf[eqn]->phi[ldof];
                if (dim == 3)
                  d_cfunc[ldof][a][var][j] +=
                      (t_veloc[1] * inv_slip * fv->dstangent_dx[1][a][jvar][j] +
                       t_veloc[1] * d_inv_slip * fv->stangent[1][a] * d_dist[jvar] *
                           bf[var]->phi[j] +
                       inv_slip * fv->stangent[1][a] * dt_veloc_dx[1][a][jvar][j]) *
                      bf[eqn]->phi[ldof];
              }
            }
          }
        }
        /*
         *  Evaluate sensitivity to velocities
         */
        for (jvar = 0; jvar < ei[pg->imtrx]->ielem_dim; jvar++) {
          var = VELOCITY1 + jvar;
          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {

                d_cfunc[ldof][a][var][j] += fv->stangent[0][jvar] * bf[var]->phi[j] * inv_slip *
                                            fv->stangent[0][a] * bf[eqn]->phi[ldof];

                if (dim == 3)
                  d_cfunc[ldof][a][var][j] += fv->stangent[1][jvar] * bf[var]->phi[j] * inv_slip *
                                              fv->stangent[1][a] * bf[eqn]->phi[ldof];
              }
            }
          }
        }
      }
    }
  }

  eqn = VELOCITY1;
  for (i = 0; i < (int)elem_side_bc->num_nodes_on_side; i++) {
    id = (int)elem_side_bc->local_elem_node_id[i];
    I = Proc_Elem_Connect[iconnect_ptr + id];
    ldof = ei[pg->imtrx]->ln_to_dof[eqn][id];
    if (Dolphin[pg->imtrx][I][VELOCITY1] > 0) {

      for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {
        cfunc[ldof][a] += force * fv->snormal[a] * bf[eqn]->phi[ldof];
        cfunc[ldof][a] += inv_slip * fv->stangent[0][a] * t_veloc[0] * bf[eqn]->phi[ldof];
      }
      if (dim == 3) {
        for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {
          cfunc[ldof][a] += inv_slip * fv->stangent[1][a] * t_veloc[1] * bf[eqn]->phi[ldof];
        }
      }
    }

  } /* end of for (i = 0; i < (int) elem_side_bc->num_nodes_on_side */

} /* END of routine apply_repulsion_user                                          */
/*****************************************************************************/

void apply_repulsion_table(double cfunc[MDE][DIM],
                           double d_cfunc[MDE][DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                           double x[],             /* solution vector*/
                           const double hscale,    /* repulsion length scale */
                           const double repexp,    /* repulsive force exponent */
                           const double P_rep,     /* repulsion coefficient */
                           const double gas_visc,  /* gas phase viscosity  */
                           const double exp_scale, /* DCL exculsion zone scale  */
                           const double v_wall[3], /* surface velocity */
                           const int dcl_node,     /* DCL NS id  */
                           struct elem_side_bc_struct *elem_side_bc,
                           const int iconnect_ptr)

/******************************************************************************
 *
 *  Function which calculates the surface repulsion from a rotating roll
 *  with a slip velocity condition
 *
 ******************************************************************************/

{

  /* Local variables */

  int dim, jvar; /* Degree of freedom counter */

  double d_dist[DIM], d_tfcn[DIM]; /* distance derivatives  */
  double dist;                     /* squared distance from surface to wall     */
  double coord[3] = {0, 0, 0};
  double force = 0.0, d_force = 0.0, inv_slip = 0.0, d_inv_slip = 0.0;
  double sheara = 0., d_sheara = 0., shearb = 0., d_shearb = 0.;
  double t_veloc[2] = {0, 0}, dt_veloc_dx[2][MAX_PDIM][MAX_PDIM][MDE];
  double n_veloc = 0, dn_veloc_dx[MAX_PDIM][MAX_PDIM][MDE];
  double dsheara_dx[MAX_PDIM][MAX_PDIM][MDE], dshearb_dx[MAX_PDIM][MAX_PDIM][MDE];

  int j, i, id, var, a, eqn, I, ldof;
  double point[3] = {0, 0, 0}, dcl_dist, slope, mod_factor = 1.;
  int nsp, k, bc_table_id = -1;
  BOUNDARY_CONDITION_STRUCT *bc_tab;

  /***************************** EXECUTION BEGINS ******************************/
  /* if pr is 0. => we don't want free surface/wall repulsion and just
     ensure that dist and dist2 are nonzero so nothing bad happens */
  if (P_rep == 0 && gas_visc == 0)
    return;

  eqn = VELOCITY1;
  if (af->Assemble_LSA_Mass_Matrix)
    return;

  /*  initialize variables */
  dim = pd->Num_Dim;

  /* calculate distance from free surface to solid surface for repulsion calculations */

  coord[0] = fv->x[0];
  coord[1] = fv->x[1];
  if (dim == 3) {
    coord[2] = fv->x[2];
  } else {
    coord[2] = 0.0;
  }

  for (a = 0; a < Num_BC; a++) {
    if (BC_Types[a].BC_Name == GD_TABLE_BC) {
      bc_table_id = a;
    }
  }
  if (bc_table_id == -1) {
    GOMA_EH(GOMA_ERROR, "GD_TABLE id not found for CAP_REPULSE_TABLE\n");
  }
  bc_tab = BC_Types + bc_table_id;
  dist = table_distance_search(bc_tab->table, coord, &slope, d_tfcn);
  if (bc_tab->table->t_index[0] == MESH_POSITION1) {
    d_dist[0] = (coord[0] - bc_tab->table->slope[1]) / dist;
    d_dist[1] = (coord[1] - bc_tab->table->slope[2]) / dist;
  } else {
    d_dist[0] = (coord[0] - bc_tab->table->slope[2]) / dist;
    d_dist[1] = (coord[1] - bc_tab->table->slope[1]) / dist;
  }
  d_dist[2] = 0.;
  if (dcl_node != -1) {
    nsp = match_nsid(dcl_node);
    k = Proc_NS_List[Proc_NS_Pointers[nsp]];
    for (j = 0; j < Proc_NS_Count[nsp]; j++) {
      k = Proc_NS_List[Proc_NS_Pointers[nsp] + j];
      i = Index_Solution(k, MESH_DISPLACEMENT1, 0, 0, -1, pg->imtrx);
      GOMA_EH(i, "Could not resolve index_solution.");
      for (a = 0; a < dim; a++) {
        point[a] = Coor[a][k] + x[i + a];
      }
    }
    dcl_dist = sqrt(SQUARE(coord[0] - point[0]) + SQUARE(coord[1] - point[1]) +
                    SQUARE(coord[2] - point[2]));
  } else {
    dcl_dist = 0.;
    GOMA_WH(-1, "No DCL node for CAP_REPULSE_TABLE....\n");
  }
  /*  modifying function for DCL  */
  mod_factor = 1. - exp(-dcl_dist / exp_scale);
  /*  repulsion function  */
  force = -P_rep * mod_factor / pow(dist / hscale, repexp);
  d_force = P_rep * mod_factor * repexp / pow(dist / hscale, repexp + 1) / hscale;
  /*  slip velocity function function  */
  inv_slip = -gas_visc * mod_factor / dist;
  d_inv_slip = gas_visc * mod_factor / SQUARE(dist);

  for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {
    t_veloc[0] += fv->stangent[0][a] * (fv->v[a] - v_wall[a]);
    n_veloc += fv->snormal[a] * (fv->v[a] - v_wall[a]);
    sheara += 0.5 * dist * d_force * d_dist[a] * fv->stangent[0][a];
    d_sheara +=
        0.5 * fv->stangent[0][a] * (d_force * d_dist[a] - d_force * (repexp + 1.) * SQUARE(hscale));
    for (jvar = 0; jvar < ei[pg->imtrx]->ielem_dim; jvar++) {
      var = MESH_DISPLACEMENT1 + jvar;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          dt_veloc_dx[0][a][jvar][j] = fv->dstangent_dx[0][a][jvar][j] * (fv->v[a] - v_wall[a]);
          dn_veloc_dx[a][jvar][j] = fv->dsnormal_dx[a][jvar][j] * (fv->v[a] - v_wall[a]);
          dn_veloc_dx[a][jvar][j] = 0.;
          dsheara_dx[a][jvar][j] =
              0.5 * dist * d_force * d_dist[a] * fv->dstangent_dx[0][a][jvar][j];
        }
      }
    }
  }
  if (dim == 3) {
    for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {
      t_veloc[1] += fv->stangent[1][a] * (fv->v[a] - v_wall[a]);
      shearb += 0.5 * dist * d_force * d_dist[a] * fv->stangent[1][a];
      d_shearb += 0.5 * fv->stangent[1][a] *
                  (d_force * d_dist[a] - d_force * (repexp + 1.) * SQUARE(hscale));
      for (jvar = 0; jvar < ei[pg->imtrx]->ielem_dim; jvar++) {
        var = MESH_DISPLACEMENT1 + jvar;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            dt_veloc_dx[1][a][jvar][j] = fv->dstangent_dx[1][a][jvar][j] * (fv->v[a] - v_wall[a]);
            dshearb_dx[a][jvar][j] =
                0.5 * dist * d_force * d_dist[a] * fv->dstangent_dx[1][a][jvar][j];
          }
        }
      }
    }
  }

  if (af->Assemble_Jacobian) {
    for (i = 0; i < (int)elem_side_bc->num_nodes_on_side; i++) {
      id = (int)elem_side_bc->local_elem_node_id[i];
      I = Proc_Elem_Connect[iconnect_ptr + id];
      ldof = ei[pg->imtrx]->ln_to_dof[eqn][id];
      if (Dolphin[pg->imtrx][I][VELOCITY1] > 0) {

        /*
         *  Evaluate sensitivity to displacements d()/dx
         */
        for (jvar = 0; jvar < ei[pg->imtrx]->ielem_dim; jvar++) {
          var = MESH_DISPLACEMENT1 + jvar;
          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {

                d_cfunc[ldof][a][var][j] += (n_veloc * force * fv->dsnormal_dx[a][jvar][j] +
                                             (n_veloc * d_force + force * dn_veloc_dx[a][jvar][j]) *
                                                 fv->snormal[a] * d_dist[jvar] * bf[var]->phi[j]) *
                                            bf[eqn]->phi[ldof];
                d_cfunc[ldof][a][var][j] +=
                    ((t_veloc[0] * inv_slip + sheara) * fv->dstangent_dx[0][a][jvar][j] +
                     (t_veloc[0] * d_inv_slip + d_sheara) * fv->stangent[0][a] * d_dist[jvar] *
                         bf[var]->phi[j] +
                     fv->stangent[0][a] *
                         (inv_slip * dt_veloc_dx[0][a][jvar][j] + dsheara_dx[a][jvar][j])) *
                    bf[eqn]->phi[ldof];
                if (dim == 3)
                  d_cfunc[ldof][a][var][j] +=
                      ((t_veloc[1] * inv_slip + shearb) * fv->dstangent_dx[1][a][jvar][j] +
                       (t_veloc[1] * d_inv_slip + d_shearb) * fv->stangent[1][a] * d_dist[jvar] *
                           bf[var]->phi[j] +
                       fv->stangent[1][a] *
                           (inv_slip * dt_veloc_dx[1][a][jvar][j] + dshearb_dx[a][jvar][j])) *
                      bf[eqn]->phi[ldof];
              }
            }
          }
        }
        /*
         *  Evaluate sensitivity to velocities
         */
        for (jvar = 0; jvar < ei[pg->imtrx]->ielem_dim; jvar++) {
          var = VELOCITY1 + jvar;
          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {

                d_cfunc[ldof][a][var][j] += fv->snormal[jvar] * bf[var]->phi[j] * force *
                                            fv->snormal[a] * bf[eqn]->phi[ldof];
                d_cfunc[ldof][a][var][j] += fv->stangent[0][jvar] * bf[var]->phi[j] * inv_slip *
                                            fv->stangent[0][a] * bf[eqn]->phi[ldof];

                if (dim == 3)
                  d_cfunc[ldof][a][var][j] += fv->stangent[1][jvar] * bf[var]->phi[j] * inv_slip *
                                              fv->stangent[1][a] * bf[eqn]->phi[ldof];
              }
            }
          }
        }
      }
    }
  }

  eqn = VELOCITY1;
  for (i = 0; i < (int)elem_side_bc->num_nodes_on_side; i++) {
    id = (int)elem_side_bc->local_elem_node_id[i];
    I = Proc_Elem_Connect[iconnect_ptr + id];
    ldof = ei[pg->imtrx]->ln_to_dof[eqn][id];
    if (Dolphin[pg->imtrx][I][VELOCITY1] > 0) {

      for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {
        cfunc[ldof][a] += force * n_veloc * fv->snormal[a] * bf[eqn]->phi[ldof];
        cfunc[ldof][a] +=
            (inv_slip * t_veloc[0] + sheara) * fv->stangent[0][a] * bf[eqn]->phi[ldof];
      }
      if (dim == 3) {
        for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {
          cfunc[ldof][a] +=
              (inv_slip * t_veloc[1] + shearb) * fv->stangent[1][a] * bf[eqn]->phi[ldof];
        }
      }
    }

  } /* end of for (i = 0; i < (int) elem_side_bc->num_nodes_on_side */

} /* END of routine apply_repulsion_table  */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void apply_vapor_recoil(double cfunc[MDE][DIM],
                        double d_cfunc[MDE][DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                        const double T_boil,  /* boiling temperature */
                        const double T_melt,  /* melting temperature */
                        const double T_ref,   /*  reference temperature */
                        const double P_scale, /* pressure scale */
                        const double T_scale, /* Temperature scale */
                        struct elem_side_bc_struct *elem_side_bc,
                        const int iconnect_ptr)

/******************************************************************************
 *
 *  Function which calculates the surface recoil from a
 *  evaportating metal alloy component or water.   This routine is currently
 *  hardwired for iron or water, with the constants.   Need to generalize by
 *  sending all of these coefficients in
 *
 ******************************************************************************/

{

  /* Local variables */
  int j, i, id, var, a, eqn, I, ldof;
  int jvar; /* Degree of freedom counter              */

  double pr, dprtemp, theta;
  double pabl_c0 = 0.0, pabl_c1 = 0.0, pabl_c2 = 0.0, pabl_c3 = 0.0;

  /***************************** EXECUTION BEGINS *******************************/

  /* Mike Kanouff's curve-fit for ablation pressure in Pascals               */
  /*                                              assume iron if T_boil>2000 */
  /*                                                      ice if T_boil<2000 */
  if (T_boil > 2000.0 * T_scale) {
    theta = fv->T - T_boil;
    if (theta > 0.0 && theta <= 170.0 * T_scale) {
      pabl_c0 = 0.0; /* iron coefficients */
      pabl_c1 = 1.8272e-4 * 1.0133e5 * (1.0 / T_scale);
      pabl_c2 = -1.9436e-6 * 1.0133e5 * (1.0 / T_scale) * (1.0 / T_scale);
      pabl_c3 = 1.5732e-8 * 1.0133e5 * (1.0 / T_scale) * (1.0 / T_scale) * (1.0 / T_scale);
    }
    if (theta > 170.0 * T_scale) {
      pabl_c0 = 0.0; /* iron coefficients */
      pabl_c1 = -5.7333e-4 * 1.0133e5 * (1.0 / T_scale);
      pabl_c2 = 4.5500e-6 * 1.0133e5 * (1.0 / T_scale) * (1.0 / T_scale);
      pabl_c3 = 2.3022e-9 * 1.0133e5 * (1.0 / T_scale) * (1.0 / T_scale) * (1.0 / T_scale);
    }
    /*pabl_c0 =  0.0;          */ /* MKS coefficients for iron */
    /*pabl_c1 =  3.723086e+02*(1.0/T_scale);*/
    /*pabl_c2 =  -6.328050e-02*(1.0/T_scale)*(1.0/T_scale);*/
    /*pabl_c3 =  5.559470e-04*(1.0/T_scale)*(1.0/T_scale)*(1.0/T_scale);*/
  } else {
    pabl_c0 = 0.0; /* MKS coefficients for ice  */
    pabl_c1 = 3.294180e+03 * (1.0 / T_scale);
    pabl_c2 = -7.726940e+00 * (1.0 / T_scale) * (1.0 / T_scale);
    pabl_c3 = 5.480973e-01 * (1.0 / T_scale) * (1.0 / T_scale) * (1.0 / T_scale);
  }

  /* Calculate ablation pressure */
  theta = fv->T - T_boil;

  if (theta < 0.) {
    theta = 0.;
    pr = 0.;
    dprtemp = 0;
  } else {
    pr = P_scale *
         (pabl_c0 + pabl_c1 * theta + pabl_c2 * pow(theta, 2.0) + pabl_c3 * pow(theta, 3.0));
    dprtemp = P_scale * (pabl_c1 + 2. * pabl_c2 * theta + 3. * pabl_c3 * pow(theta, 2.0));
  }

  eqn = VELOCITY1;

  if (af->Assemble_Jacobian) {
    for (i = 0; i < (int)elem_side_bc->num_nodes_on_side; i++) {
      id = (int)elem_side_bc->local_elem_node_id[i];
      I = Proc_Elem_Connect[iconnect_ptr + id];
      ldof = ei[pg->imtrx]->ln_to_first_dof[eqn][id];

      /* Calculate the residual contribution from surface gradient of basis function */
      if (Dolphin[pg->imtrx][I][VELOCITY1] > 0) {

        /* Calculate the residual contribution from surface gradient of basis function */
        /*
         *  Evaluate sensitivity to displacements d()/dx
         */
        for (jvar = 0; jvar < ei[pg->imtrx]->ielem_dim; jvar++) {
          var = MESH_DISPLACEMENT1 + jvar;
          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {
                d_cfunc[ldof][a][var][j] -= (pr)*fv->dsnormal_dx[a][jvar][j] * bf[eqn]->phi[ldof];
              }
            }
          }
        }
        /*
         *  Evaluate sensitivity to conc and temp (if surface tension is
         *    variable
         */
        var = TEMPERATURE;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {

              d_cfunc[ldof][a][var][j] -= dprtemp * fv->snormal[a] * bf[eqn]->phi[ldof];
            }
          }
        }
      }
    }
  }

  eqn = VELOCITY1;
  if (pd->v[pg->imtrx][eqn]) {
    for (i = 0; i < (int)elem_side_bc->num_nodes_on_side; i++) {
      id = (int)elem_side_bc->local_elem_node_id[i];
      I = Proc_Elem_Connect[iconnect_ptr + id];
      ldof = ei[pg->imtrx]->ln_to_first_dof[eqn][id];

      /* Calculate the residual contribution */
      if (ldof >= 0) {
        /*
         * Add to the boundary stress in the fluid mechanics momentum equations
         */

        for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {
          cfunc[ldof][a] -= pr * fv->snormal[a] * bf[eqn]->phi[ldof];
        }
      }
    }
  } /* if velocity variable is defined.   */
  else {
    GOMA_EH(GOMA_ERROR, "Bad coord system in capillary or liquid momentum equations not defined");
  }

} /* END of routine apply_vapor_recoil                                          */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void flow_n_dot_T_hydro(double func[DIM],
                        double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                        const double a,
                        const double b,
                        const double c,
                        const double d) /* 4 parms describing pressure variation */

/************************************************************************
 *
 * flow_n_dot_T_hydro()
 *
 *  Function which calculates the contribution to a the normal stress
 *  from a specified pressure along an interface.
 *     The specified pressure can be a linear function of the spatial
 *  coordinates.
 *     Specifically, this function returns :
 *
 *                   func[p] = - Pressure * Surface_Normal_dot_dir[p]
 *
 *  The Jacobian dependence of func[p] wrt to the mesh displacement
 *  unknowns are returned in d_func[p][var_type][j].
 ************************************************************************/
{
  int j, var, p, jvar;
  double press, coeff[DIM];

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  if (af->Assemble_Jacobian) {
    coeff[0] = a;
    coeff[1] = b;
    coeff[2] = c;
    press = a * fv->x[0] + b * fv->x[1] + c * fv->x[2] + d;
    for (jvar = 0; jvar < ei[pg->imtrx]->ielem_dim; jvar++) {
      var = MESH_DISPLACEMENT1 + jvar;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (p = 0; p < pd->Num_Dim; p++) {
            d_func[p][var][j] -= (bf[var]->phi[j] * coeff[jvar] * fv->snormal[p] +
                                  press * fv->dsnormal_dx[p][jvar][j]);
          }
        }
      }
    }
  }

  /*
   * Calculate the pressure at current gauss point via a linear function
   * of the spatial coordinates.
   */
  press = a * fv->x[0] + b * fv->x[1] + c * fv->x[2] + d;
  for (p = 0; p < pd->Num_Dim; p++) {
    func[p] = -press * fv->snormal[p];
  }
} /* END of routine flow_n_dot_T_hydro                                       */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void flow_n_dot_T_var_density(
    double func[DIM],
    double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
    const double a,    /* 1 param describing reference pressure*/
    const double time) /* Time is required for density and momentum_source_term*/

/************************************************************************
 *
 * flow_n_dot_T_var_density()
 *
 *  Function which calculates the contribution to a the normal stress
 *  from a specified pressure along an interface.
 *  The specified pressure is:
 *                   P = P_0 + rho*g*x
 *  where the user specifies P_0, rho is density, and g is taken from
 *  momentum_source_term.
 *     Specifically, this function returns :
 *
 *                   func[p] = - Pressure * Surface_Normal_dot_dir[p]
 *
 *  The Jacobian dependence of func[p] wrt to the mesh displacement, temp,
 *  and concentration unknowns are returned in d_func[p][var_type][j].
 ************************************************************************/
{
  int j, var, p, jvar, k, w;
  double press, rho, f[DIM], f_dot_x;
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;
  MOMENTUM_SOURCE_DEPENDENCE_STRUCT df_struct; /* Body force dependence */
  MOMENTUM_SOURCE_DEPENDENCE_STRUCT *df = &df_struct;

  rho = density(d_rho, time);              // Density
  (void)momentum_source_term(f, df, time); // Momentum source term for gravity vector

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  if (af->Assemble_Jacobian) {
    press = rho * (f[0] * fv->x[0] + f[1] * fv->x[1] + f[2] * fv->x[2]) + a;
    for (jvar = 0; jvar < ei[pg->imtrx]->ielem_dim; jvar++) {
      var = MESH_DISPLACEMENT1 + jvar; // Jacobian wrt mesh displacement
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (p = 0; p < pd->Num_Dim; p++) {
            d_func[p][var][j] -= rho * f[jvar] * bf[var]->phi[j] * fv->snormal[p]; // dx term
            for (k = 0; k < pd->Num_Dim; k++) {
              d_func[p][var][j] -= rho * df->X[k][jvar][j] * fv->x[k] *
                                   fv->snormal[p]; // df term, should be zero if f=gravity
            }
            d_func[p][var][j] -= press * fv->dsnormal_dx[p][jvar][j]; // dn term
          }
        }
      }
    }

    f_dot_x = f[0] * fv->x[0] + f[1] * fv->x[1] + f[2] * fv->x[2];

    var = TEMPERATURE; // Jacobian wrt TEMPERATURE
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      for (p = 0; p < pd->Num_Dim; p++) {
        d_func[p][var][j] -= d_rho->T[j] * f_dot_x * fv->snormal[p];
      }
    }

    var = MASS_FRACTION; // Jacobian wrt concentration
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      for (p = 0; p < pd->Num_Dim; p++) {
        for (w = 0; w < pd->Num_Species; w++) {
          d_func[p][MAX_VARIABLE_TYPES + w][j] -= d_rho->C[w][j] * f_dot_x * fv->snormal[p];
        }
      }
    }
  }

  /*
   * Calculate the pressure at current gauss point
   */
  press = rho * (f[0] * fv->x[0] + f[1] * fv->x[1] + f[2] * fv->x[2]) + a;
  for (p = 0; p < pd->Num_Dim; p++) {
    func[p] = -press * fv->snormal[p];
  }
} /* END of routine flow_n_dot_T_var_density                                       */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#if 0

/* Deprecated March 2002 by TAB */
void
hydrostatic_n_dot_T(double *func,
		    double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE])
    
    /*************************************************************************
     *
     * hydrostatic_n_dot_T()
     *
     *  Function which calculates the contribution to a the normal stress
     *  from the hyrostatic pressure. In other words, if you use this
     *  function, you are taking out the implicit specification of the
     *  hydrostatic pressure at an interface.
     *
     *     Specifically, this function returns :
     *
     *                   func[p] = - Pressure * Surface_Normal_dot_dir[p]
     *
     *  Note, the pressure used above is the actual pressure unknown.
     *
     *  The Jacobian dependence of func[p] wrt to the the pressure and the
     *  mesh displacement unknowns are returned in d_func[p][var_type][j].
     ************************************************************************/
{
  int j, var, p, jvar;

  if (af->Assemble_LSA_Mass_Matrix) return;

  if (af->Assemble_Jacobian) {
    var = PRESSURE;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
	for (p = 0; p < pd->Num_Dim; p++) {
	  d_func[p][var][j] -= bf[var]->phi[j] * fv->snormal[p];
	}
      }
    }
    for (jvar = 0; jvar < ei[pg->imtrx]->ielem_dim; jvar++) {
      var = MESH_DISPLACEMENT1 + jvar;
      if (pd->v[pg->imtrx][var]) {
	for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
	  for (p = 0; p < pd->Num_Dim; p++){
	    d_func[0][var][j] -= fv->P * fv->dsnormal_dx[p][jvar][j];
	  }
	}
      }
    }
  }
  /*  calculate pressure at current gauss point */
  for (p = 0; p < pd->Num_Dim; p++) {
    func[p] = -fv->P * fv->snormal[p];
  }
  GOMA_WH(-1, "Hydrostatic Symmetry condition may cause singular matrix!");
     
} /* END of routine hydrostatic_n_dot_T                                      */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
#endif

void flow_n_dot_T_nobc(double func[DIM],
                       double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                       const double pdatum, /* pressure datum from input card */
                       const int iflag,     /* -1 to use pdatum, otherwise use P */
                       const double time)   /* time value if needed	*/

/****************************************************************************
 *
 *  Function which uses the fully-developed fluid stresses
 *		 in the fluid stress bcs
 *
 *****************************************************************************/
{
  int j, var, p, q;
  int b, c, mode, w;
  double P = fv->P;

  /*
   * Variables for stress tensor and derivative
   */

  dbl Pi[DIM][DIM];
  STRESS_DEPENDENCE_STRUCT d_Pi_struct;
  STRESS_DEPENDENCE_STRUCT *d_Pi = &d_Pi_struct;

  int v_s[MAX_MODES][DIM][DIM];
  int v_g[DIM][DIM];

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  (void)stress_eqn_pointer(v_s);

  v_g[0][0] = VELOCITY_GRADIENT11;
  v_g[0][1] = VELOCITY_GRADIENT12;
  v_g[1][0] = VELOCITY_GRADIENT21;
  v_g[1][1] = VELOCITY_GRADIENT22;
  v_g[0][2] = VELOCITY_GRADIENT13;
  v_g[1][2] = VELOCITY_GRADIENT23;
  v_g[2][0] = VELOCITY_GRADIENT31;
  v_g[2][1] = VELOCITY_GRADIENT32;
  v_g[2][2] = VELOCITY_GRADIENT33;

  if (!af->Assemble_Jacobian)
    d_Pi = NULL; /* save cost if we aren't assembling the Jacobian */

  /* if we are using the datum for pressure, put it into fv->P
     and we'll fix it after calculating the stress tensor
   */
  if (iflag == -1)
    fv->P = pdatum;

  /* hydrostatic part for int == -2  ; func[] is set to n.p_hydrostatic in
          flow_n_dot_T_var_density. "var_density" is a misnomer - density would need
          to be spatially invariate since pressure is relative to p(origin). Even
          density(T) or density(P) would imply T or P is spatially invariate
  */
  if (iflag == -2) {
    flow_n_dot_T_var_density(func, d_func, pdatum, time);
    fv->P = 0;
  }

  /* compute stress tensor and its derivatives */
  fluid_stress(Pi, d_Pi);

  /* now is the time to clean up, so, if using the datum for pressure, fix fv->P
   */
  /* Ahh...no.  This seems like a bad idea.
    if(iflag == -1) fv->P = P;*/

  if (af->Assemble_Jacobian) {
    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var]) {

      for (p = 0; p < pd->Num_Dim; p++) {
        for (q = 0; q < pd->Num_Dim; q++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_func[p][var][j] += fv->snormal[q] * d_Pi->T[p][q][j];
          }
        }
      }
    }

    var = BOND_EVOLUTION;
    if (pd->v[pg->imtrx][var]) {

      for (p = 0; p < pd->Num_Dim; p++) {
        for (q = 0; q < pd->Num_Dim; q++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_func[p][var][j] += fv->snormal[q] * d_Pi->nn[p][q][j];
          }
        }
      }
    }

    var = EDDY_NU;
    if (pd->v[pg->imtrx][var]) {

      for (p = 0; p < pd->Num_Dim; p++) {
        for (q = 0; q < pd->Num_Dim; q++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_func[p][var][j] += fv->snormal[q] * d_Pi->eddy_nu[p][q][j];
          }
        }
      }
    }

    var = RESTIME;
    if (pd->v[pg->imtrx][var]) {

      for (p = 0; p < pd->Num_Dim; p++) {
        for (q = 0; q < pd->Num_Dim; q++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_func[p][var][j] += fv->snormal[q] * d_Pi->degrade[p][q][j];
          }
        }
      }
    }

    var = FILL;
    if (pd->v[pg->imtrx][var]) {
      for (p = 0; p < pd->Num_Dim; p++) {
        for (q = 0; q < pd->Num_Dim; q++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_func[p][var][j] += fv->snormal[q] * d_Pi->F[p][q][j];
          }
        }
      }
    }

    var = MASS_FRACTION;
    if (pd->v[pg->imtrx][var]) {
      for (p = 0; p < pd->Num_Dim; p++) {
        for (q = 0; q < pd->Num_Dim; q++) {
          for (w = 0; w < pd->Num_Species_Eqn; w++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_func[p][MAX_VARIABLE_TYPES + w][j] += fv->snormal[q] * d_Pi->C[p][q][w][j];
            }
          }
        }
      }
    }

    if (iflag != -1) {
      var = PRESSURE;
      if (pd->v[pg->imtrx][var]) {
        for (p = 0; p < pd->Num_Dim; p++) {
          for (q = 0; q < pd->Num_Dim; q++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_func[p][var][j] += fv->snormal[q] * d_Pi->P[p][q][j];
            }
          }
        }
      }
    }

    if (pd->v[pg->imtrx][VELOCITY1]) {
      for (p = 0; p < pd->Num_Dim; p++) {
        for (q = 0; q < pd->Num_Dim; q++) {
          for (b = 0; b < VIM; b++) {
            var = VELOCITY1 + b;
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_func[p][var][j] += fv->snormal[q] * d_Pi->v[p][q][b][j];
            }
          }
        }
      }
    }

    if (pd->v[pg->imtrx][VORT_DIR1]) {
      for (p = 0; p < pd->Num_Dim; p++) {
        for (q = 0; q < pd->Num_Dim; q++) {
          for (b = 0; b < VIM; b++) {
            var = VORT_DIR1 + b;
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_func[p][var][j] += fv->snormal[q] * d_Pi->vd[p][q][b][j];
            }
          }
        }
      }
    }

    if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
      for (p = 0; p < pd->Num_Dim; p++) {
        for (q = 0; q < pd->Num_Dim; q++) {
          for (b = 0; b < VIM; b++) {
            var = MESH_DISPLACEMENT1 + b;
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_func[p][var][j] +=
                  fv->snormal[q] * d_Pi->X[p][q][b][j] + fv->dsnormal_dx[q][b][j] * Pi[p][q];
            }
          }
        }
      }
    }

    if (pd->v[pg->imtrx][POLYMER_STRESS11]) {
      for (p = 0; p < pd->Num_Dim; p++) {
        for (q = 0; q < pd->Num_Dim; q++) {
          for (mode = 0; mode < vn->modes; mode++) {
            for (b = 0; b < VIM; b++) {
              for (c = 0; c < VIM; c++) {
                var = v_s[mode][b][c];
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)

                {
                  d_func[p][var][j] += fv->snormal[q] * d_Pi->S[p][q][mode][b][c][j];
                }
              }
            }
          }
        }
      }
    }

    if (gn->ConstitutiveEquation == BINGHAM_MIXED ||
        (pd->v[pg->imtrx][POLYMER_STRESS11] &&
         (vn->evssModel == EVSS_F || vn->evssModel == SQRT_CONF || vn->evssModel == LOG_CONF ||
          vn->evssModel == EVSS_GRADV || vn->evssModel == CONF ||
          vn->evssModel == LOG_CONF_GRADV))) {
      for (p = 0; p < pd->Num_Dim; p++) {
        for (q = 0; q < pd->Num_Dim; q++) {
          for (b = 0; b < VIM; b++) {
            for (c = 0; c < VIM; c++) {
              var = v_g[b][c];
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                d_func[p][var][j] -= fv->snormal[q] * d_Pi->g[p][q][b][c][j];
              }
            }
          }
        }
      }
    }

  } /*  end of if Assemble_Jacobian  */

  /*  load in stresses */

  for (p = 0; p < pd->Num_Dim; p++) {
    if (iflag != -2) {
      func[p] = 0.;
    }
    for (q = 0; q < pd->Num_Dim; q++) {
      func[p] += fv->snormal[q] * Pi[p][q];
    }
  }
  if (iflag == -1 || iflag == -2)
    fv->P = P;

} /* END of routine flow_n_dot_T_nobc                                       */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void flow_n_dot_T_gradv(double func[DIM],
                        double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                        const double pdatum, /* pressure datum from input card */
                        const int iflag)     /* -1 to use pdatum, otherwise use P  */

/****************************************************************************
 *
 *  Function which uses the fully-developed fluid stresses
 *		 in the fluid stress bcs
 *
 *****************************************************************************/
{
  int i, j, var, p, q;
  int a, b;
  double press;

  /*
   * Variables for vicosity and derivative
   */
  dbl gamma[DIM][DIM];
  dbl mu;
  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct; /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  if (iflag != -1) {
    press = fv->P;
  } else {
    press = pdatum;
  }

  /**
        compute gammadot, viscosity
   **/
  for (i = 0; i < VIM; i++) {
    for (j = 0; j < VIM; j++) {
      gamma[i][j] = fv->grad_v[i][j] + fv->grad_v[j][i];
    }
  }

  mu = viscosity(gn, gamma, d_mu);
  if (af->Assemble_Jacobian) {

    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        for (p = 0; p < pd->Num_Dim; p++) {
          for (q = 0; q < pd->Num_Dim; q++) {
            d_func[p][var][j] += fv->snormal[q] * d_mu->T[j] * fv->grad_v[q][p];
          }
        }
      }
    }

    if (iflag != -1) {
      var = PRESSURE;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (p = 0; p < pd->Num_Dim; p++) {
            d_func[p][var][j] += -bf[var]->phi[j] * fv->snormal[p];
          }
        }
      }
    }

    if (pd->v[pg->imtrx][VELOCITY1]) {
      for (a = 0; a < VIM; a++) {
        var = VELOCITY1 + a;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (p = 0; p < pd->Num_Dim; p++) {
            for (q = 0; q < pd->Num_Dim; q++) {
              d_func[p][var][j] += fv->snormal[q] * (mu * bf[var]->grad_phi_e[j][a][q][p] +
                                                     fv->grad_v[q][p] * d_mu->v[a][j]);
            }
          }
        }
      }
    }
    if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
      for (b = 0; b < VIM; b++) {
        var = MESH_DISPLACEMENT1 + b;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (p = 0; p < pd->Num_Dim; p++) {
            d_func[p][var][j] += -press * fv->dsnormal_dx[p][b][j];

            for (q = 0; q < pd->Num_Dim; q++) {
              d_func[p][var][j] += fv->snormal[q] * (mu * fv->d_grad_v_dmesh[q][p][b][j] +
                                                     fv->grad_v[q][p] * d_mu->X[b][j]) +
                                   fv->dsnormal_dx[q][b][j] * mu * fv->grad_v[q][p];
            }
          }
        }
      }
    }

  } /*  end of if Assemble_Jacobian  */

  /*  load in stresses */

  for (p = 0; p < pd->Num_Dim; p++) {
    func[p] = -press * fv->snormal[p];
    for (q = 0; q < pd->Num_Dim; q++) {
      func[p] += fv->snormal[q] * mu * fv->grad_v[q][p];
    }
  }

} /* END of routine flow_n_dot_T_gradv                                      */
/*****************************************************************************/

void flow_n_dot_T_gradv_t(double func[DIM],
                          double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                          const double pdatum, /* pressure datum from input card */
                          const int iflag)     /* -1 to use pdatum, otherwise use P  */

/****************************************************************************
 *
 *  Function which uses the fully-developed fluid stresses
 *		 in the fluid stress bcs
 *
 *  Imposes n . nabla v = 0 naturally
 *
 *****************************************************************************/
{
  int j, var, p, q;
  int b, c, mode, w;
  double P = fv->P;

  /*
   * Variables for stress tensor and derivative
   */

  dbl Pi[DIM][DIM];
  STRESS_DEPENDENCE_STRUCT d_Pi_struct;
  STRESS_DEPENDENCE_STRUCT *d_Pi = &d_Pi_struct;

  int v_s[MAX_MODES][DIM][DIM];
  int v_g[DIM][DIM];

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  (void)stress_eqn_pointer(v_s);

  v_g[0][0] = VELOCITY_GRADIENT11;
  v_g[0][1] = VELOCITY_GRADIENT12;
  v_g[1][0] = VELOCITY_GRADIENT21;
  v_g[1][1] = VELOCITY_GRADIENT22;
  v_g[0][2] = VELOCITY_GRADIENT13;
  v_g[1][2] = VELOCITY_GRADIENT23;
  v_g[2][0] = VELOCITY_GRADIENT31;
  v_g[2][1] = VELOCITY_GRADIENT32;
  v_g[2][2] = VELOCITY_GRADIENT33;

  if (!af->Assemble_Jacobian)
    d_Pi = NULL; /* save cost if we aren't assembling the Jacobian */

  /* if we are using the datum for pressure, put it into fv->P
     and we'll fix it after calculating the stress tensor
   */
  if (iflag == -1)
    fv->P = pdatum;

  /* hydrostatic part for int == -2  ; func[] is set to n.p_hydrostatic in
          flow_n_dot_T_var_density. "var_density" is a misnomer - density would need
          to be spatially invariate since pressure is relative to p(origin). Even
          density(T) or density(P) would imply T or P is spatially invariate
  */
  // if (iflag == -2) {
  //   flow_n_dot_T_var_density(func, d_func, pdatum, time);
  //   fv->P = 0;
  // }

  bool evss_f = false;

  if (pd->gv[POLYMER_STRESS11]) {
    evss_f = vn->evssModel == LOG_CONF || vn->evssModel == LOG_CONF_GRADV ||
             vn->evssModel == CONF || vn->evssModel == SQRT_CONF || vn->evssModel == EVSS_F ||
             vn->evssModel == EVSS_GRADV;
  }
  /* compute stress tensor and its derivatives */
  fluid_stress(Pi, d_Pi);

  VISCOSITY_DEPENDENCE_STRUCT d_mus_struct; /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mus = &d_mus_struct;
  VISCOSITY_DEPENDENCE_STRUCT d_mup_struct; /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mup = &d_mup_struct;
  zeroStructures(d_mup, 1);
  dbl gamma[DIM][DIM]; // Shear rate tensor based on velocity
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
    }
  }
  dbl mus = viscosity(gn, gamma, d_mus);
  dbl mup = 0;
  if (evss_f) {
    for (mode = 0; mode < vn->modes; mode++) {
      /*  shift factor  */
      dbl at = 0.0;
      dbl d_at_dT[MDE];
      dbl wlf_denom;
      if (pd->e[pg->imtrx][TEMPERATURE]) {
        if (vn->shiftModel == CONSTANT) {
          at = vn->shift[0];
          for (j = 0; j < ei[pg->imtrx]->dof[TEMPERATURE]; j++) {
            d_at_dT[j] = 0.;
          }
        } else if (vn->shiftModel == MODIFIED_WLF) {
          wlf_denom = vn->shift[1] + fv->T - mp->reference[TEMPERATURE];
          if (wlf_denom != 0.) {
            at = exp(vn->shift[0] * (mp->reference[TEMPERATURE] - fv->T) / wlf_denom);
            for (j = 0; j < ei[pg->imtrx]->dof[TEMPERATURE]; j++) {
              d_at_dT[j] = -at * vn->shift[0] * vn->shift[1] / (wlf_denom * wlf_denom) *
                           bf[TEMPERATURE]->phi[j];
            }
          } else {
            at = 1.;
          }
          for (j = 0; j < ei[pg->imtrx]->dof[TEMPERATURE]; j++) {
            d_at_dT[j] = 0.;
          }
        }
      } else {
        at = 1.;
      }
      int a, j, w;
      VISCOSITY_DEPENDENCE_STRUCT d_muptmp_struct; /* viscosity dependence */
      VISCOSITY_DEPENDENCE_STRUCT *d_muptmp = &d_muptmp_struct;
      // Polymer viscosity
      dbl muptmp = viscosity(ve[mode]->gn, gamma, d_muptmp);

      mup += muptmp;
      var = VELOCITY1;
      if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
        for (a = 0; a < WIM; a++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_mup->v[a][j] += at * d_muptmp->v[a][j];
          }
        }
      }

      var = MESH_DISPLACEMENT1;
      if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
        for (a = 0; a < pd->Num_Dim; a++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_mup->X[a][j] += at * d_muptmp->X[a][j];
          }
        }
      }

      var = TEMPERATURE;
      if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_mup->T[j] += (at * d_muptmp->T[j] + mup * d_at_dT[j]);
        }
      }

      var = BOND_EVOLUTION;
      if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_mup->nn[j] += at * d_muptmp->nn[j];
        }
      }

      var = RESTIME;
      if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_mup->degrade[j] += at * d_muptmp->degrade[j];
        }
      }

      var = FILL;
      if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_mup->F[j] += at * d_muptmp->F[j];
        }
      }

      var = PRESSURE;
      if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_mup->P[j] += at * d_muptmp->P[j];
        }
      }

      var = MASS_FRACTION;
      if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_mup->C[w][j] += at * d_muptmp->C[w][j];
          }
        }
      }
    } // for mode
  }   // if evss_f
  /* now is the time to clean up, so, if using the datum for pressure, fix fv->P
   */
  /* Ahh...no.  This seems like a bad idea.
    if(iflag == -1) fv->P = P;*/

  if (af->Assemble_Jacobian) {
    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var]) {

      for (p = 0; p < pd->Num_Dim; p++) {
        for (q = 0; q < pd->Num_Dim; q++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_func[p][var][j] += fv->snormal[q] * d_Pi->T[p][q][j] -
                                 ((d_mus->T[j] + d_mup->T[j]) * fv->grad_v[p][q] -
                                  evss_f * d_mup->T[j] * (fv->grad_v[p][q] - fv->G[p][q]));
          }
        }
      }
    }

    var = BOND_EVOLUTION;
    if (pd->v[pg->imtrx][var]) {

      for (p = 0; p < pd->Num_Dim; p++) {
        for (q = 0; q < pd->Num_Dim; q++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_func[p][var][j] += fv->snormal[q] * d_Pi->nn[p][q][j];
          }
        }
      }
    }

    var = RESTIME;
    if (pd->v[pg->imtrx][var]) {

      for (p = 0; p < pd->Num_Dim; p++) {
        for (q = 0; q < pd->Num_Dim; q++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_func[p][var][j] += fv->snormal[q] * d_Pi->degrade[p][q][j];
          }
        }
      }
    }

    var = FILL;
    if (pd->v[pg->imtrx][var]) {
      for (p = 0; p < pd->Num_Dim; p++) {
        for (q = 0; q < pd->Num_Dim; q++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_func[p][var][j] += fv->snormal[q] * d_Pi->F[p][q][j] -
                                 ((d_mus->F[j] + d_mup->F[j]) * fv->grad_v[p][q] -
                                  evss_f * d_mup->F[j] * (fv->grad_v[p][q] - fv->G[p][q]));
          }
        }
      }
    }

    var = MASS_FRACTION;
    if (pd->v[pg->imtrx][var]) {
      for (p = 0; p < pd->Num_Dim; p++) {
        for (q = 0; q < pd->Num_Dim; q++) {
          for (w = 0; w < pd->Num_Species_Eqn; w++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_func[p][MAX_VARIABLE_TYPES + w][j] +=
                  fv->snormal[q] * d_Pi->C[p][q][w][j] -
                  ((d_mus->C[w][j] + d_mup->C[w][j]) * fv->grad_v[p][q] -
                   evss_f * d_mup->C[w][j] * (fv->grad_v[p][q] - fv->G[p][q]));
            }
          }
        }
      }
    }

    if (iflag != -1) {
      var = PRESSURE;
      if (pd->v[pg->imtrx][var]) {
        for (p = 0; p < pd->Num_Dim; p++) {
          for (q = 0; q < pd->Num_Dim; q++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_func[p][var][j] += fv->snormal[q] * d_Pi->P[p][q][j] -
                                   ((d_mus->P[j] + d_mup->P[j]) * fv->grad_v[p][q] -
                                    evss_f * d_mup->P[j] * (fv->grad_v[p][q] - fv->G[p][q]));
            }
          }
        }
      }
    }

    if (pd->v[pg->imtrx][VELOCITY1]) {
      for (p = 0; p < pd->Num_Dim; p++) {
        for (q = 0; q < pd->Num_Dim; q++) {
          for (b = 0; b < VIM; b++) {
            var = VELOCITY1 + b;
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_func[p][var][j] +=
                  fv->snormal[q] * (d_Pi->v[p][q][b][j] -
                                    ((d_mus->v[b][j] + evss_f * d_mup->v[b][j]) * d_mup->v[b][j] *
                                     fv->grad_v[p][q]) -
                                    ((mus + evss_f * mup) * bf[var]->grad_phi_e[j][b][p][q]));
            }
          }
        }
      }
    }

    if (pd->v[pg->imtrx][VORT_DIR1]) {
      for (p = 0; p < pd->Num_Dim; p++) {
        for (q = 0; q < pd->Num_Dim; q++) {
          for (b = 0; b < VIM; b++) {
            var = VORT_DIR1 + b;
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_func[p][var][j] += fv->snormal[q] * d_Pi->vd[p][q][b][j];
            }
          }
        }
      }
    }

    if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
      for (p = 0; p < pd->Num_Dim; p++) {
        for (q = 0; q < pd->Num_Dim; q++) {
          for (b = 0; b < VIM; b++) {
            var = MESH_DISPLACEMENT1 + b;
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_func[p][var][j] += fv->snormal[q] * d_Pi->X[p][q][b][j] +
                                   fv->dsnormal_dx[q][b][j] * Pi[p][q] -
                                   ((d_mus->X[b][j] + d_mup->X[b][j]) * fv->grad_v[p][q] -
                                    evss_f * d_mup->X[b][j] * (fv->grad_v[p][q] - fv->G[p][q])) -
                                   ((mus + mup) * fv->d_grad_v_dmesh[p][q][b][j]);
            }
          }
        }
      }
    }

    if (pd->v[pg->imtrx][POLYMER_STRESS11]) {
      for (p = 0; p < pd->Num_Dim; p++) {
        for (q = 0; q < pd->Num_Dim; q++) {
          for (mode = 0; mode < vn->modes; mode++) {
            for (b = 0; b < VIM; b++) {
              for (c = 0; c < VIM; c++) {
                var = v_s[mode][b][c];
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)

                {
                  d_func[p][var][j] += fv->snormal[q] * d_Pi->S[p][q][mode][b][c][j];
                }
              }
            }
          }
        }
      }
    }

    if (gn->ConstitutiveEquation == BINGHAM_MIXED ||
        (pd->v[pg->imtrx][POLYMER_STRESS11] &&
         (vn->evssModel == EVSS_F || vn->evssModel == SQRT_CONF || vn->evssModel == LOG_CONF ||
          vn->evssModel == EVSS_GRADV || vn->evssModel == CONF ||
          vn->evssModel == LOG_CONF_GRADV))) {
      for (p = 0; p < pd->Num_Dim; p++) {
        for (q = 0; q < pd->Num_Dim; q++) {
          for (b = 0; b < VIM; b++) {
            for (c = 0; c < VIM; c++) {
              var = v_g[b][c];
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                d_func[p][var][j] -=
                    fv->snormal[q] * (d_Pi->g[p][q][b][c][j] -
                                      evss_f * mup * delta(b, p) * delta(c, q) * bf[var]->phi[j]);
              }
            }
          }
        }
      }
    }

  } /*  end of if Assemble_Jacobian  */

  /*  load in stresses */

  for (p = 0; p < pd->Num_Dim; p++) {
    if (iflag != -2) {
      func[p] = 0.;
    }
    for (q = 0; q < pd->Num_Dim; q++) {
      func[p] +=
          fv->snormal[q] *
          (Pi[p][q] - (mus * fv->grad_v[p][q] + evss_f * mup * (fv->grad_v[p][q] - fv->G[p][q])));
    }
  }
  if (iflag == -1 || iflag == -2)
    fv->P = P;

} /* END of routine flow_n_dot_T_gradv                                      */

/* FLOW_GRADV_SIC Strongly integrated condition */
void flow_n_dot_T_gradv_sic(double func[DIM],
                            double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                            const double pdatum, /* pressure datum from input card */
                            const int iflag)     /* -1 to use pdatum, otherwise use P  */

/****************************************************************************
 *
 *  Function which uses the fully-developed fluid stresses
 *		 in the fluid stress bcs
 *
 *****************************************************************************/
{
  int i, j, var, p, q;
  int a, b;

  /*
   * Variables for vicosity and derivative
   */
  dbl gamma[DIM][DIM];
  dbl mu;
  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct; /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  /**
        compute gammadot, viscosity
   **/
  for (i = 0; i < VIM; i++) {
    for (j = 0; j < VIM; j++) {
      gamma[i][j] = fv->grad_v[i][j] + fv->grad_v[j][i];
    }
  }

  mu = viscosity(gn, gamma, d_mu);
  if (af->Assemble_Jacobian) {

    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        for (p = 0; p < pd->Num_Dim; p++) {
          for (q = 0; q < pd->Num_Dim; q++) {
            d_func[p][var][j] += fv->snormal[q] * d_mu->T[j] * fv->grad_v[p][q];
          }
        }
      }
    }

    if (iflag != -1) {
      var = PRESSURE;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (p = 0; p < pd->Num_Dim; p++) {
            d_func[p][var][j] += 0.;
          }
        }
      }
    }

    if (pd->v[pg->imtrx][VELOCITY1]) {
      for (a = 0; a < VIM; a++) {
        var = VELOCITY1 + a;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (p = 0; p < pd->Num_Dim; p++) {
            for (q = 0; q < pd->Num_Dim; q++) {
              d_func[p][var][j] += fv->snormal[q] * (mu * bf[var]->grad_phi_e[j][a][p][q] +
                                                     fv->grad_v[p][q] * d_mu->v[a][j]);
            }
          }
        }
      }
    }
    if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
      for (b = 0; b < VIM; b++) {
        var = MESH_DISPLACEMENT1 + b;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (p = 0; p < pd->Num_Dim; p++) {
            d_func[p][var][j] += 0.;

            for (q = 0; q < pd->Num_Dim; q++) {
              d_func[p][var][j] += fv->snormal[q] * (mu * fv->d_grad_v_dmesh[p][q][b][j] +
                                                     fv->grad_v[p][q] * d_mu->X[b][j]) +
                                   fv->dsnormal_dx[q][b][j] * mu * fv->grad_v[p][q];
            }
          }
        }
      }
    }

  } /*  end of if Assemble_Jacobian  */

  /*  load in stresses */

  for (p = 0; p < pd->Num_Dim; p++) {
    for (q = 0; q < pd->Num_Dim; q++) {
      func[p] += fv->snormal[q] * mu * fv->grad_v[p][q];
    }
  }

} /* END of routine flow_n_dot_T_gradv_sic                                   */
/*****************************************************************************/

void stress_no_v_dot_gradS(double func[MAX_MODES][6],
                           double d_func[MAX_MODES][6][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                           const double dt,
                           const double tt)

/****************************************************************************
 *
 *  Function which imposes fully-developed fluid stresses
 *  in the velocity direction
 *
 *  It's a clone for assemble_stress_fortin with the v dot grad S term removed
 *
 *  For more details see:
 *
 *    Xueying Xie and Matteo Pasquali J. Non-Newtonian Fluid Mech. 122 (2004) 159 - 176
 *
 *  Kristianto Tjiptowidjojo (tjiptowi@unm.edu)
 *****************************************************************************/
{

  int j, eqn, var, a, b, p, q, mode, w, k;
  int R_s[MAX_MODES][DIM][DIM];
  int v_s[MAX_MODES][DIM][DIM];
  int inv_v_s[DIM][DIM];
  int v_g[DIM][DIM];
  int dim = pd->Num_Dim;

  dbl grad_v[DIM][DIM]; /* Velocity gradient based on velocity - discontinuous across element */
  dbl gamma[DIM][DIM];  /* Shear-rate tensor based on velocity */

  dbl s[DIM][DIM];     /* stress tensor */
  dbl trace = 0.0;     /* trace of the stress tensor */
  dbl s_dot[DIM][DIM]; /* stress tensor from last time step */
  dbl g[DIM][DIM];     /* velocity gradient tensor */
  dbl gt[DIM][DIM];    /* transpose of velocity gradient tensor */

  dbl g_dot[DIM][DIM];  /* velocity gradient tensor time derivative */
  dbl gt_dot[DIM][DIM]; /* transpose of velocity gradient tensor time derivative */

  /* dot product tensors */

  dbl s_dot_s[DIM][DIM];
  dbl s_dot_g[DIM][DIM];
  dbl s_dot_gt[DIM][DIM];
  dbl g_dot_s[DIM][DIM];
  dbl gt_dot_s[DIM][DIM];

  dbl g_dot_g[DIM][DIM];
  dbl gt_dot_g[DIM][DIM];
  dbl gt_dot_gt[DIM][DIM];

  /* polymer viscosity and derivatives */
  dbl mup;
  VISCOSITY_DEPENDENCE_STRUCT d_mup_struct;
  VISCOSITY_DEPENDENCE_STRUCT *d_mup = &d_mup_struct;
  dbl d_mup_dv_pj;

  dbl saramitoCoeff = 1.;
  SARAMITO_DEPENDENCE_STRUCT d_saramito_struct;
  SARAMITO_DEPENDENCE_STRUCT *d_saramito = &d_saramito_struct;
  /* 2nd polymer viscosity -- for modified Jefreys model*/
  dbl mupJeff;

  const bool jeffreysEnabled = (vn->ConstitutiveEquation == MODIFIED_JEFFREYS);
  // todo: will want to parse necessary parameters... for now hard code
  const bool saramitoEnabled =
      (vn->ConstitutiveEquation == SARAMITO_OLDROYDB || vn->ConstitutiveEquation == SARAMITO_PTT ||
       vn->ConstitutiveEquation == SARAMITO_GIESEKUS);

  /*  shift function */
  dbl at = 0.0;
  dbl d_at_dT[MDE];
  dbl wlf_denom;

  /* constitutive equation parameters */
  dbl alpha;      /* This is the Geisekus mobility parameter */
  dbl lambda = 0; /* polymer relaxation constant */
  dbl ucwt, lcwt; /* Upper convected derviative weight, Lower convected derivative weight */
  dbl eps;        /* This is the PTT elongation parameter */
  dbl Z = 1.0;    /* This is the factor appearing in front of the stress tensor in PTT */
  dbl dZ_dtrace = 0.0;

  dbl lambda1 = 0;
  dbl lambda2 = 0;
  dbl elasticMod;

  /* ETMs*/
  dbl mass, advection, source;

  dbl phi_j;

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  /* A bookeeping array from stress[i][j] into its position in func */
  inv_v_s[0][0] = 0; /* S11 */
  inv_v_s[0][1] = 1; /* S12 */
  inv_v_s[1][1] = 2; /* S22 */
  inv_v_s[0][2] = 3; /* S13 */
  inv_v_s[1][2] = 4; /* S23 */
  inv_v_s[2][2] = 5; /* S33 */

  (void)stress_eqn_pointer(v_s);
  (void)stress_eqn_pointer(R_s);

  v_g[0][0] = VELOCITY_GRADIENT11;
  v_g[0][1] = VELOCITY_GRADIENT12;
  v_g[1][0] = VELOCITY_GRADIENT21;
  v_g[1][1] = VELOCITY_GRADIENT22;
  v_g[0][2] = VELOCITY_GRADIENT13;
  v_g[1][2] = VELOCITY_GRADIENT23;
  v_g[2][0] = VELOCITY_GRADIENT31;
  v_g[2][1] = VELOCITY_GRADIENT32;
  v_g[2][2] = VELOCITY_GRADIENT33;

  /*
   * Load up Field variables...
   */

  /****  Velocity gradient, and shear rate ****/
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      grad_v[a][b] = fv->grad_v[a][b];
    }
  }
  /* load up shear rate tensor based on velocity */
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      gamma[a][b] = grad_v[a][b] + grad_v[b][a];
    }
  }
  /* Continuous velocity gradient field and its transpose */
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      g[a][b] = fv->G[a][b];
      gt[b][a] = g[a][b];
    }
  }

  /*  Time-temperature shift factor  */
  if (pd->e[pg->imtrx][TEMPERATURE]) {
    if (vn->shiftModel == CONSTANT) {
      at = vn->shift[0];
      for (j = 0; j < ei[pg->imtrx]->dof[TEMPERATURE]; j++) {
        d_at_dT[j] = 0.;
      }
    } else if (vn->shiftModel == MODIFIED_WLF) {
      wlf_denom = vn->shift[1] + fv->T - mp->reference[TEMPERATURE];
      if (wlf_denom != 0.) {
        at = exp(vn->shift[0] * (mp->reference[TEMPERATURE] - fv->T) / wlf_denom);
        for (j = 0; j < ei[pg->imtrx]->dof[TEMPERATURE]; j++) {
          d_at_dT[j] =
              -at * vn->shift[0] * vn->shift[1] / (wlf_denom * wlf_denom) * bf[TEMPERATURE]->phi[j];
        }
      } else {
        at = 1.;
      }
      for (j = 0; j < ei[pg->imtrx]->dof[TEMPERATURE]; j++) {
        d_at_dT[j] = 0.;
      }
    }
  } else {
    at = 1.;
  }

  // if a modified Jeffreys model is being run, load time derivative of velocity gradient
  if (jeffreysEnabled) {
    (void)tensor_dot(g, g, g_dot_g, VIM);
    (void)tensor_dot(gt, gt, gt_dot_gt, VIM);
    (void)tensor_dot(gt, g, gt_dot_g, VIM);

    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {

        if (pd->TimeIntegration != STEADY) {
          g_dot[a][b] = fv_dot->G[a][b];
        } else {
          g_dot[a][b] = 0.;
        }
        gt_dot[b][a] = g_dot[a][b];
      }
    }
  }

  /* Begin loop over modes */
  for (mode = 0; mode < vn->modes; mode++) {

    /*
     * Load polymeric stress tensor for each mode
     */
    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        s[a][b] = fv->S[mode][a][b];
        if (pd->TimeIntegration != STEADY) {
          s_dot[a][b] = fv_dot->S[mode][a][b];
        } else {
          s_dot[a][b] = 0.;
        }
      }
    }

    /* Calculate trace of stress tensor */
    trace = 0.0;
    for (a = 0; a < VIM; a++) {
      trace += s[a][a];
    }

    /* get polymer viscosity */
    mup = viscosity(ve[mode]->gn, gamma, d_mup);

    if (saramitoEnabled == TRUE) {
      compute_saramito_model_terms(&saramitoCoeff, d_saramito, s, ve[mode]->gn, FALSE);
    } else {
      saramitoCoeff = 1.;
      d_saramito->tau_y = 0;

      for (int i = 0; i < VIM; ++i) {
        for (int j = 0; j < VIM; ++j) {
          d_saramito->s[i][j] = 0;
        }
      }
    }
    /* get Geisekus mobility parameter */
    alpha = ve[mode]->alpha;

    /* get time constant */
    if (ve[mode]->time_constModel == CONSTANT) {
      lambda = ve[mode]->time_const;
    } else if (ve[mode]->time_constModel == CARREAU || ve[mode]->time_constModel == POWER_LAW) {
      lambda = mup / ve[mode]->time_const;
    }

    ucwt = 1.0 - ve[mode]->xi / 2.0;
    lcwt = ve[mode]->xi / 2.0;

    eps = ve[mode]->eps;

    if (jeffreysEnabled) {
      mupJeff = ve[mode]->muJeffreys;
      // if the modified Jeffreys model is used, the parsed value of lambda is the
      // elastic modulus rather than the time consant
      elasticMod = lambda;
      lambda1 = mup / elasticMod; // mup/G
      lambda2 = mupJeff / elasticMod;
      lambda = lambda1 + lambda2;
    }

    Z = 1.0;
    dZ_dtrace = 0;
    if (vn->ConstitutiveEquation == PTT) {
      if (vn->ptt_type == PTT_LINEAR) {
        Z = 1 + eps * lambda * trace / mup;
        dZ_dtrace = eps * lambda / mup;
      } else if (vn->ptt_type == PTT_EXPONENTIAL) {
        Z = exp(eps * lambda * trace / mup);
        dZ_dtrace = Z * eps * lambda / mup;
      } else {
        GOMA_EH(GOMA_ERROR, "Unrecognized PTT Form %d", vn->ptt_type);
      }
    }
    /* get tensor dot products for future use */

    if (alpha != 0.)
      (void)tensor_dot(s, s, s_dot_s, VIM);

    if (ucwt != 0.) {
      (void)tensor_dot(s, g, s_dot_g, VIM);
      (void)tensor_dot(gt, s, gt_dot_s, VIM);
    }

    if (lcwt != 0.) {
      (void)tensor_dot(s, gt, s_dot_gt, VIM);
      (void)tensor_dot(g, s, g_dot_s, VIM);
    }

    /**** Assemble func *****/

    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {

        if (a <= b) /* since the stress tensor is symmetric, only assemble the upper half */
        {
          eqn = R_s[mode][a][b];
          k = inv_v_s[a][b];

          /* Mass term */
          mass = 0.;
          if (pd->TimeIntegration != STEADY) {
            mass = at * lambda * s_dot[a][b];
            mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
          }

          /* Advection term - minus v_dot_grad_S and x_dot_del_S terms */
          advection = 0.;
          if (lambda != 0.) {
            if (ucwt != 0.)
              advection -= ucwt * (gt_dot_s[a][b] + s_dot_g[a][b]);
            if (lcwt != 0.)
              advection += lcwt * (s_dot_gt[a][b] + g_dot_s[a][b]);
            advection *= at * lambda;
            advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
          }

          /* Source term */
          source = 0.;
          source += saramitoCoeff * Z * s[a][b] - at * mup * (g[a][b] + gt[a][b]);

          if (DOUBLE_NONZERO(alpha)) {
            dbl source1 = (s_dot_s[a][b] / mup);

            source1 *= alpha * lambda * saramitoCoeff;
            source += source1;
          }

          if (jeffreysEnabled) {
            source -= mup * lambda2 *
                      (g_dot[a][b] + gt_dot[a][b] -
                       (g_dot_g[a][b] + 2 * gt_dot_g[a][b] + gt_dot_gt[a][b]));
          }
          source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

          func[mode][k] += mass + advection + source;
        }
      }
    }

    /**** Assemble d_func *****/

    if (af->Assemble_Jacobian) {
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          if (a <= b) /* since the stress tensor is symmetric, only assemble the upper half */
          {
            eqn = R_s[mode][a][b];
            k = inv_v_s[a][b];

            /* Sensitivities w.r.t. velocity - J_S_v */
            var = VELOCITY1;
            if (pd->v[pg->imtrx][var]) {
              for (p = 0; p < dim; p++) {
                var = VELOCITY1 + p;

                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  d_mup_dv_pj = d_mup->v[p][j];

                  /* mass term */
                  mass = 0.;

                  /* advection term */
                  advection = 0.;

                  /* source term */
                  source = -at * d_mup_dv_pj * (g[a][b] + gt[a][b]);
                  if (alpha != 0.0) {
                    source -=
                        saramitoCoeff * alpha * lambda * d_mup_dv_pj * s_dot_s[a][b] / (mup * mup);
                  }

                  source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

                  /* Load them up */
                  d_func[mode][k][var][j] += mass + advection + source;
                }
              }
            } /* End of J_S_v */

            /* Sensitivities w.r.t. pressure - J_S_p */
            var = PRESSURE;
            if (pd->v[pg->imtrx][var]) {
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                /* mass term */
                mass = 0.;

                /* advection term */
                advection = 0.;

                /* source term */
                source = -at * d_mup->P[j] * (g[a][b] + gt[a][b]);
                if (alpha != 0.0) {
                  source -=
                      saramitoCoeff * alpha * lambda * d_mup->P[j] * s_dot_s[a][b] / (mup * mup);
                }
                source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

                /* Load them up */
                d_func[mode][k][var][j] += mass + advection + source;
              }
            } /* End of J_S_p */

            /* Sensitivities w.r.t. temperature - J_S_T */
            var = TEMPERATURE;
            if (pd->v[pg->imtrx][var]) {
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

                /* mass term */
                mass = 0.;
                if (pd->TimeIntegration != STEADY) {
                  mass = s_dot[a][b] * d_at_dT[j] * lambda;
                  mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                }

                /* advection term */
                advection = 0.;
                if (lambda != 0.) {
                  if (ucwt != 0.)
                    advection -= ucwt * (gt_dot_s[a][b] + s_dot_g[a][b]);
                  if (lcwt != 0.)
                    advection += lcwt * (s_dot_gt[a][b] + g_dot_s[a][b]);
                  advection *= d_at_dT[j] * lambda;
                  advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                }

                /* source term */
                source = -(g[a][b] + gt[a][b]) * (at * d_mup->T[j] + mup * d_at_dT[j]);
                if (alpha != 0.0) {
                  source -=
                      saramitoCoeff * alpha * lambda * d_mup->T[j] * s_dot_s[a][b] / (mup * mup);
                }
                source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

                /* Load them up */
                d_func[mode][k][var][j] += mass + advection + source;
              }
            } /* End of J_S_T */

            /* Sensitivities w.r.t. mesh displacement - J_S_d */
            for (p = 0; p < dim; p++) {
              var = MESH_DISPLACEMENT1 + p;
              if (pd->v[pg->imtrx][var]) {

                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

                  /* mass term */
                  mass = 0.0;

                  /* advection term */
                  advection = 0.0;

                  /* source term */
                  source = 0.0;
                  source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

                  /* Load them up */
                  d_func[mode][k][var][j] += mass + advection + source;
                }
              }
            } /* End of J_S_d */

            /* Sensitivities w.r.t. mass fraction - J_S_c */
            var = MASS_FRACTION;
            if (pd->v[pg->imtrx][var]) {
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                for (w = 0; w < pd->Num_Species_Eqn; w++) {

                  /* mass term */
                  mass = 0.0;

                  /* advection term */
                  advection = 0.0;

                  /* source term */
                  source = 0.0;
                  source = -at * d_mup->C[w][j] * (g[a][b] + gt[a][b]);
                  if (alpha != 0.) {
                    source -= saramitoCoeff * alpha * lambda * d_mup->C[w][j] * s_dot_s[a][b] /
                              (mup * mup);
                  }
                  source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

                  /* Load them up */
                  d_func[mode][k][MAX_VARIABLE_TYPES + w][j] += mass + advection + source;
                }
              }
            } /* End of J_S_c */
            /* Sensitivities w.r.t. continuous velocity gradient - J_S_G */
            var = VELOCITY_GRADIENT11;
            if (pd->v[pg->imtrx][var]) {
              for (p = 0; p < VIM; p++) {
                for (q = 0; q < VIM; q++) {
                  var = v_g[p][q];

                  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                    phi_j = bf[var]->phi[j];

                    /* mass term */
                    mass = 0.0;

                    /* advection term */
                    advection = 0.0;
                    if (lambda != 0.) {
                      advection -= ucwt * phi_j *
                                   (s[p][b] * (double)delta(a, q) + s[a][p] * (double)delta(b, q));
                      advection += lcwt * phi_j *
                                   (s[a][q] * (double)delta(p, b) + s[q][b] * (double)delta(a, p));
                      advection *= at * lambda;
                    }
                    advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

                    /* source term */
                    source = -at * mup * phi_j *
                             ((double)delta(a, p) * (double)delta(b, q) +
                              (double)delta(b, p) * (double)delta(a, q));
                    source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

                    /* Load them up */
                    d_func[mode][k][var][j] += mass + advection + source;
                  }
                }
              }
            } /* End of J_S_G */

            /* Sensitivities w.r.t. polymeric stress - J_S_S */
            var = POLYMER_STRESS11;
            if (pd->v[pg->imtrx][var]) {
              for (p = 0; p < VIM; p++) {
                for (q = 0; q < VIM; q++) {
                  var = v_s[mode][p][q];

                  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                    phi_j = bf[var]->phi[j];

                    /* mass term */
                    mass = 0.0;
                    if (pd->TimeIntegration != STEADY) {
                      mass =
                          (1. + 2. * tt) * phi_j / dt * (double)delta(a, p) * (double)delta(b, q);
                      mass *= at * lambda;
                    }
                    mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];

                    /* advection term */
                    advection = 0.0;
                    if (lambda != 0.0) {
                      advection -= phi_j * ucwt *
                                   (gt[a][p] * (double)delta(b, q) + g[q][b] * (double)delta(a, p));
                      advection += phi_j * lcwt *
                                   (gt[q][b] * (double)delta(p, a) + g[a][p] * (double)delta(q, b));
                      advection *= at * lambda;
                    }
                    advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

                    /* source term */
                    source = saramitoCoeff * Z * phi_j * (double)delta(a, p) * (double)delta(b, q);
                    if (p == q)
                      source += s[a][b] * dZ_dtrace * phi_j;
                    if (p <= q) {
                      source += phi_j * d_saramito->s[p][q] * Z * s[a][b];
                    }

                    if (alpha != 0.) {
                      source += saramitoCoeff * phi_j * alpha * lambda *
                                (s[q][b] * (double)delta(a, p) + s[a][p] * (double)delta(b, q)) /
                                mup;
                      if (p <= q) {
                        source +=
                            phi_j * d_saramito->s[p][q] * alpha * lambda * s_dot_s[a][b] / mup;
                      }
                    }
                    source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

                    /* Load them up */
                    d_func[mode][k][var][j] += mass + advection + source;
                  }
                }
              }
            } /* End of J_S_S */
          }   /* End of if a == b */
        }     /* End of loop over dimension "b" */
      }       /* End of loop over dimension "a" */
    }         /*  end of if Assemble_Jacobian  */
  }           /* End of loop over modes */
} /* END of routine stress_no_v_dot_gradS                                    */
/*****************************************************************************/

void stress_no_v_dot_gradS_logc(double func[MAX_MODES][6],
                                double d_func[MAX_MODES][6][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                                const double dt,
                                const double tt)

/****************************************************************************
 *
 *  Function which imposes fully-developed fluid stresses
 *  in the velocity direction
 *
 *  It's a clone for assemble_stress_fortin with the v dot grad S term removed
 *
 *  For more details see:
 *
 *    Xueying Xie and Matteo Pasquali J. Non-Newtonian Fluid Mech. 122 (2004) 159 - 176
 *
 *  Kristianto Tjiptowidjojo (tjiptowi@unm.edu)
 *
 *  This routine is adjusted for the log-conformation tensor
 *****************************************************************************/
{

  int i, j, eqn, a, b, mode, w, k;
  int siz;
  int R_s[MAX_MODES][DIM][DIM];
  int v_s[MAX_MODES][DIM][DIM];
  int inv_v_s[DIM][DIM];
  int logc_gradv = 0;

  dbl grad_v[DIM][DIM]; /* Velocity gradient based on velocity - discontinuous across element */
  dbl gamma[DIM][DIM];  /* Shear-rate tensor based on velocity */

  dbl s[DIM][DIM];     /* stress tensor (log-conformation tensor) */
  dbl exp_s[DIM][DIM]; // Exponential of log_conf
  dbl trace = 0.0;     /* trace of the stress tensor */
  dbl s_dot[DIM][DIM]; /* stress tensor from last time step */
  dbl gt[DIM][DIM];    /* transpose of velocity gradient tensor */

  dbl source_term1[DIM][DIM];
  dbl advection_term1[DIM][DIM];

  /* polymer viscosity and derivatives */
  dbl mup;
  VISCOSITY_DEPENDENCE_STRUCT d_mup_struct;
  VISCOSITY_DEPENDENCE_STRUCT *d_mup = &d_mup_struct;

  /*  shift function */
  dbl at = 0.0;
  dbl wlf_denom;

  // Decomposition of velocity vector
  dbl eig_values[DIM];
  dbl R1[DIM][DIM];
  dbl R1_T[DIM][DIM];
  dbl Rt_dot_gradv[DIM][DIM];
  dbl D[DIM][DIM];
  dbl D_dot_D[DIM][DIM];
  dbl M1[DIM][DIM];

  dbl tmp1[DIM][DIM], tmp2[DIM][DIM], tmp3[DIM][DIM];

  /* constitutive equation parameters */
  dbl alpha;      /* This is the Geisekus mobility parameter */
  dbl lambda = 0; /* polymer relaxation constant */
  dbl eps;        /* This is the PTT elongation parameter */
  dbl Z = 1.0;    /* This is the factor appearing in front of the stress tensor in PTT */

  /* ETMs*/
  dbl mass, advection, source;

  dbl d_lambda;

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  /* A bookeeping array from stress[i][j] into its position in func */
  inv_v_s[0][0] = 0; /* S11 */
  inv_v_s[0][1] = 1; /* S12 */
  inv_v_s[1][1] = 2; /* S22 */
  inv_v_s[0][2] = 3; /* S13 */
  inv_v_s[1][2] = 4; /* S23 */
  inv_v_s[2][2] = 5; /* S33 */

  (void)stress_eqn_pointer(v_s);
  (void)stress_eqn_pointer(R_s);

  memset(exp_s, 0, sizeof(double) * DIM * DIM);
  memset(R1, 0, sizeof(double) * DIM * DIM);
  memset(eig_values, 0, sizeof(double) * DIM);

  if (vn->evssModel == LOG_CONF_GRADV) {
    logc_gradv = 1;
  }

  /*
   * Load up Field variables...
   */

  /****  Velocity gradient, and shear rate ****/
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      grad_v[a][b] = fv->grad_v[a][b];
    }
  }
  /* load up shear rate tensor based on velocity */
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      gamma[a][b] = grad_v[a][b] + grad_v[b][a];
    }
  }
  /* Continuous velocity gradient field and its transpose */
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      gt[b][a] = fv->G[a][b];
    }
  }

  /*  Time-temperature shift factor  */
  if (pd->e[pg->imtrx][TEMPERATURE]) {
    if (vn->shiftModel == CONSTANT) {
      at = vn->shift[0];
    } else if (vn->shiftModel == MODIFIED_WLF) {
      wlf_denom = vn->shift[1] + fv->T - mp->reference[TEMPERATURE];
      if (wlf_denom != 0.) {
        at = exp(vn->shift[0] * (mp->reference[TEMPERATURE] - fv->T) / wlf_denom);
      } else {
        at = 1.;
      }
    }
  } else {
    at = 1.;
  }

  /* Begin loop over modes */
  for (mode = 0; mode < vn->modes; mode++) {

    /*
     * Load polymeric stress tensor for each mode
     */
    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        s[a][b] = fv->S[mode][a][b];
        if (pd->TimeIntegration != STEADY) {
          s_dot[a][b] = fv_dot->S[mode][a][b];
        } else {
          s_dot[a][b] = 0.;
        }
      }
    }

    // Polymer viscosity
    mup = viscosity(ve[mode]->gn, gamma, d_mup);

    // Giesekus mobility parameter
    alpha = ve[mode]->alpha;

    // Polymer time constant
    if (ve[mode]->time_constModel == CONSTANT) {
      lambda = ve[mode]->time_const;
    } else if (ve[mode]->time_constModel == CARREAU || ve[mode]->time_constModel == POWER_LAW) {
      lambda = mup / ve[mode]->time_const;
    }

#ifdef ANALEIG_PLEASE
    analytical_exp_s(s, exp_s, eig_values, R1, NULL);
#else
    compute_exp_s(s, exp_s, eig_values, R1);
#endif

    // Decompose velocity gradient

    memset(D, 0, sizeof(double) * DIM * DIM);
    D[0][0] = eig_values[0];
    D[1][1] = eig_values[1];
    if (VIM > 2) {
      D[2][2] = eig_values[2];
    }
    (void)tensor_dot(D, D, D_dot_D, VIM);

    // Decompose velocity gradient

    memset(M1, 0, sizeof(double) * DIM * DIM);
    memset(R1_T, 0, sizeof(double) * DIM * DIM);

    for (i = 0; i < VIM; i++) {
      for (j = 0; j < VIM; j++) {
        R1_T[i][j] = R1[j][i];
      }
    }

    for (i = 0; i < VIM; i++) {
      for (j = 0; j < VIM; j++) {
        Rt_dot_gradv[i][j] = 0.;
        for (w = 0; w < VIM; w++) {
          if (logc_gradv) {
            Rt_dot_gradv[i][j] += R1_T[i][w] * grad_v[j][w];
          } else {
            Rt_dot_gradv[i][j] += R1_T[i][w] * gt[w][j];
          }
        }
      }
    }

    for (i = 0; i < VIM; i++) {
      for (j = 0; j < VIM; j++) {
        M1[i][j] = 0.;
        for (w = 0; w < VIM; w++) {
          M1[i][j] += Rt_dot_gradv[i][w] * R1[w][j];
        }
      }
    }

    // Predetermine advective terms
    trace = eig_values[0] + eig_values[1];
    if (VIM > 2) {
      trace += eig_values[2];
    }

    // PTT exponent
    eps = ve[mode]->eps;

    // PTT
    Z = 1;
    if (vn->ConstitutiveEquation == PTT) {
      if (vn->ptt_type == PTT_LINEAR) {
        Z = 1 + eps * (trace - (double)VIM);
      } else if (vn->ptt_type == PTT_EXPONENTIAL) {
        Z = exp(eps * (trace - (double)VIM));
      } else {
        GOMA_EH(GOMA_ERROR, "Unrecognized PTT Form %d", vn->ptt_type);
      }
    }

    siz = sizeof(double) * DIM * DIM;
    memset(tmp1, 0, siz);
    memset(tmp2, 0, siz);
    memset(tmp3, 0, siz);
    memset(advection_term1, 0, siz);
    memset(source_term1, 0, siz);

    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        if (a != b) {
          d_lambda = eig_values[b] - eig_values[a];
          if (fabs(d_lambda) > 1.e-8) {
            tmp1[a][b] += (log(eig_values[b]) - log(eig_values[a])) / d_lambda;
            tmp1[a][b] *= (eig_values[a] * M1[b][a] + eig_values[b] * M1[a][b]);
          } else {
            tmp1[a][b] += M1[a][b] + M1[b][a];
          }
        }
        if (a == b) {
          source_term1[a][b] += Z * (1.0 - D[a][a]) / lambda;
          if (alpha != 0) {
            source_term1[a][b] += alpha * (2.0 * D[a][a] - 1.0 - D_dot_D[a][a]) / lambda;
          }
          source_term1[a][b] /= eig_values[a];
          source_term1[a][b] += 2.0 * M1[a][a];
        }
      }
    }

    (void)tensor_dot(R1, tmp1, tmp2, VIM);
    (void)tensor_dot(tmp2, R1_T, advection_term1, VIM);
    (void)tensor_dot(R1, source_term1, tmp3, VIM);
    (void)tensor_dot(tmp3, R1_T, source_term1, VIM);

    /**** Assemble func *****/

    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {

        if (a <= b) /* since the stress tensor is symmetric, only assemble the upper half */
        {
          eqn = R_s[mode][a][b];
          k = inv_v_s[a][b];

          /* Mass term */
          mass = 0.;
          if (pd->TimeIntegration != STEADY) {
            mass = at * s_dot[a][b];
            mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
          }

          /* Advection term - minus v_dot_grad_S and x_dot_del_S terms */
          advection = 0.;
          if (lambda != 0.) {
            advection -= advection_term1[a][b];
            advection *= at * pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
          }

          /* Source term */
          source = 0.;
          source -= source_term1[a][b];
          source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

          func[mode][k] += mass + advection + source;
        }
      }
    }
  } /* End of loop over modes */
} /* END of routine stress_no_v_dot_gradS_logc                                 */

void stress_no_v_dot_gradS_sqrt(double func[MAX_MODES][6],
                                double d_func[MAX_MODES][6][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                                const double dt,
                                const double tt) {

  int dim, p, q, w, k;

  int eqn;
  int evss_gradv = 0;

  int i, j, mode;
  dbl v[DIM];      /* Velocity field. */
  dbl x_dot[DIM];  /* current position field derivative wrt time. */
  dbl h3;          /* Volume element (scale factors). */
  dbl dh3dmesh_pj; /* Sensitivity to (p,j) mesh dof. */

  dbl grad_v[DIM][DIM];
  dbl gamma[DIM][DIM]; /* Shear-rate tensor based on velocity */
  dbl det_J;           /* determinant of element Jacobian */

  dbl d_det_J_dmesh_pj; /* for specific (p,j) mesh dof */

  dbl mass; /* For terms and their derivatives */
  dbl mass_a, mass_b;
  dbl advection;
  dbl advection_a, advection_b, advection_c, advection_d;
  dbl diffusion;
  dbl source;
  dbl source1;
  dbl source_a = 0, source_b = 0, source_c = 0;
  int err;
  dbl alpha = 0;  /* This is the Geisekus mobility parameter */
  dbl lambda = 0; /* polymer relaxation constant */
  dbl d_lambda_dF[MDE];
  double xi;
  double d_xi_dF[MDE];
  dbl eps = 0; /* This is the PTT elongation parameter */
  double d_eps_dF[MDE];
  /*
   *
   * Note how carefully we avoid refering to d(phi[i])/dx[j] and refer instead
   * to the j-th component of grad_phi[j][i] so that this vector can be loaded
   * up with components that may be different in non Cartesian coordinate
   * systems.
   *
   * We will, however, insist on *orthogonal* coordinate systems, even if we
   * might permit them to be curvilinear.
   *
   * Assume all components of velocity are interpolated with the same kind
   * of basis function.
   */

  /*
   * Petrov-Galerkin weighting functions for i-th and ab-th stress residuals
   * and some of their derivatives...
   */

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;
  dbl wt;

  /* Variables for stress */

  int R_s[MAX_MODES][DIM][DIM];
  int v_s[MAX_MODES][DIM][DIM];
  int v_g[DIM][DIM];

  dbl b[DIM][DIM];     /* stress tensor */
  dbl b_dot[DIM][DIM]; /* stress tensor from last time step */
  dbl grad_b[DIM][DIM][DIM];
  dbl d_grad_s_dmesh[DIM][DIM][DIM][DIM]
                    [MDE]; /* derivative of grad of stress tensor for mode ve_mode */

  dbl g[DIM][DIM];  /* velocity gradient tensor */
  dbl gt[DIM][DIM]; /* transpose of velocity gradient tensor */

  /* dot product tensors */

  dbl s_dot_s[DIM][DIM];
  dbl b_dot_g[DIM][DIM];

  /* polymer viscosity and derivatives */
  dbl mup;
  VISCOSITY_DEPENDENCE_STRUCT d_mup_struct;
  VISCOSITY_DEPENDENCE_STRUCT *d_mup = &d_mup_struct;

  // todo: will want to parse necessary parameters... for now hard code
  const bool saramitoEnabled =
      (vn->ConstitutiveEquation == SARAMITO_OLDROYDB || vn->ConstitutiveEquation == SARAMITO_PTT ||
       vn->ConstitutiveEquation == SARAMITO_GIESEKUS);

  dbl saramitoCoeff = 1.;

  dbl d_mup_dv_pj;
  dbl d_mup_dmesh_pj;

  /*  shift function */
  dbl at = 0.0;
  dbl d_at_dT[MDE];
  dbl wlf_denom;

  /* advective terms are precalculated */
  dbl v_dot_del_b[DIM][DIM];
  dbl x_dot_del_b[DIM][DIM];

  dbl d_xdotdels_dm;

  dbl d_vdotdels_dm;
  int inv_v_s[DIM][DIM];

  if (vn->evssModel == EVSS_GRADV) {
    evss_gradv = 1;
  }

  eqn = R_STRESS11;
  inv_v_s[0][0] = 0; /* S11 */
  inv_v_s[0][1] = 1; /* S12 */
  inv_v_s[1][1] = 2; /* S22 */
  inv_v_s[0][2] = 3; /* S13 */
  inv_v_s[1][2] = 4; /* S23 */
  inv_v_s[2][2] = 5; /* S33 */

  /*
   * Bail out fast if there's nothing to do...
   */

  if (!pd->e[pg->imtrx][eqn]) {
    return;
  }

  /*
   * Unpack variables from structures for local convenience...
   */

  dim = pd->Num_Dim;

  wt = fv->wt;

  det_J = bf[eqn]->detJ; /* Really, ought to be mesh eqn. */

  h3 = fv->h3; /* Differential volume element (scales). */

  /* load eqn and variable number in tensor form */
  (void)stress_eqn_pointer(v_s);
  (void)stress_eqn_pointer(R_s);

  v_g[0][0] = VELOCITY_GRADIENT11;
  v_g[0][1] = VELOCITY_GRADIENT12;
  v_g[1][0] = VELOCITY_GRADIENT21;
  v_g[1][1] = VELOCITY_GRADIENT22;
  v_g[0][2] = VELOCITY_GRADIENT13;
  v_g[1][2] = VELOCITY_GRADIENT23;
  v_g[2][0] = VELOCITY_GRADIENT31;
  v_g[2][1] = VELOCITY_GRADIENT32;
  v_g[2][2] = VELOCITY_GRADIENT33;

  /*
   * Field variables...
   */
  for (int a = 0; a < WIM; a++) {
    v[a] = fv->v[a];

    /* note, these are zero for steady calculations */
    x_dot[a] = 0.0;
    if (pd->TimeIntegration != STEADY && pd->gv[MESH_DISPLACEMENT1 + a]) {
      x_dot[a] = fv_dot->x[a];
    }
  }

  /*
   * In Cartesian coordinates, this velocity gradient tensor will
   * have components that are...
   *
   * 			grad_v[a][b] = d v_b
   *				       -----
   *				       d x_a
   */

  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      grad_v[a][b] = fv->grad_v[a][b];
    }
  }

  /* load up shearrate tensor based on velocity */
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      gamma[a][b] = grad_v[a][b] + grad_v[b][a];
    }
  }

  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      if (evss_gradv) {
        g[a][b] = fv->grad_v[a][b];
        gt[a][b] = fv->grad_v[b][a];
      } else {
        g[a][b] = fv->G[a][b];
        gt[b][a] = g[a][b];
      }
    }
  }

  /*  shift factor  */
  if (pd->gv[TEMPERATURE]) {
    if (vn->shiftModel == CONSTANT) {
      at = vn->shift[0];
      for (j = 0; j < ei[pg->imtrx]->dof[TEMPERATURE]; j++) {
        d_at_dT[j] = 0.;
      }
    } else if (vn->shiftModel == MODIFIED_WLF) {
      wlf_denom = vn->shift[1] + fv->T - mp->reference[TEMPERATURE];
      if (wlf_denom != 0.) {
        at = exp(vn->shift[0] * (mp->reference[TEMPERATURE] - fv->T) / wlf_denom);
        for (j = 0; j < ei[pg->imtrx]->dof[TEMPERATURE]; j++) {
          d_at_dT[j] =
              -at * vn->shift[0] * vn->shift[1] / (wlf_denom * wlf_denom) * bf[TEMPERATURE]->phi[j];
        }
      } else {
        at = 1.;
      }
      for (j = 0; j < ei[pg->imtrx]->dof[TEMPERATURE]; j++) {
        d_at_dT[j] = 0.;
      }
    }
  } else {
    at = 1.;
  }

  /* Begin loop over modes */
  for (mode = 0; mode < vn->modes; mode++) {

    load_modal_pointers(mode, tt, dt, b, b_dot, grad_b, d_grad_s_dmesh);

    dbl source_term[DIM][DIM];
    dbl d_source_term_db[DIM][DIM][DIM][DIM];
    sqrt_conf_source(mode, b, source_term, d_source_term_db);

    /* precalculate advective terms of form (v dot del tensor)*/

    /*
     * Stress tensor...(Note "anti-BSL" sign convention on deviatoric stress)
     */
    for (int ii = 0; ii < VIM; ii++) {
      for (int jj = 0; jj < VIM; jj++) {
        v_dot_del_b[ii][jj] = 0.;
        x_dot_del_b[ii][jj] = 0.;
        for (q = 0; q < WIM; q++) {
          v_dot_del_b[ii][jj] += v[q] * grad_b[q][ii][jj];
          x_dot_del_b[ii][jj] += x_dot[q] * grad_b[q][ii][jj];
        }
      }
    }

    /* get polymer viscosity */
    mup = viscosity(ve[mode]->gn, gamma, d_mup);

    if (saramitoEnabled == TRUE) {
      GOMA_EH(GOMA_ERROR, "Saramito not enabled sqrt");
    }

    double d_alpha_dF[MDE];
    /* get Geisekus mobility parameter */
    if (ve[mode]->alphaModel == CONSTANT) {
      alpha = ve[mode]->alpha;
    } else if (ls != NULL && ve[mode]->alphaModel == VE_LEVEL_SET) {
      double pos_alpha = ve[mode]->pos_ls.alpha;
      double neg_alpha = ve[mode]->alpha;
      double width = ls->Length_Scale;
      err = level_set_property(neg_alpha, pos_alpha, width, &alpha, d_alpha_dF);
      GOMA_EH(err, "level_set_property() failed for mobility parameter.");
    } else {
      GOMA_EH(GOMA_ERROR, "Unknown mobility parameter model");
    }

    /* get time constant */
    if (ve[mode]->time_constModel == CONSTANT) {
      lambda = ve[mode]->time_const;
    } else if (ve[mode]->time_constModel == CARREAU || ve[mode]->time_constModel == POWER_LAW) {
      lambda = mup / ve[mode]->time_const;
    } else if (ls != NULL && ve[mode]->time_constModel == VE_LEVEL_SET) {
      double pos_lambda = ve[mode]->pos_ls.time_const;
      double neg_lambda = ve[mode]->time_const;
      double width = ls->Length_Scale;
      err = level_set_property(neg_lambda, pos_lambda, width, &lambda, d_lambda_dF);
      GOMA_EH(err, "level_set_property() failed for polymer time constant.");
    }

    xi = 0;
    if (ve[mode]->xiModel == CONSTANT) {
      xi = ve[mode]->xi;
    } else if (ls != NULL && ve[mode]->xiModel == VE_LEVEL_SET) {
      double pos_xi = ve[mode]->pos_ls.xi;
      double neg_xi = ve[mode]->xi;
      double width = ls->Length_Scale;
      err = level_set_property(neg_xi, pos_xi, width, &xi, d_xi_dF);
      GOMA_EH(err, "level_set_property() failed for ptt xi parameter.");
    } else {
      GOMA_EH(GOMA_ERROR, "Unknown PTT Xi parameter model");
    }

    if (ve[mode]->epsModel == CONSTANT) {
      eps = ve[mode]->eps;
    } else if (ls != NULL && ve[mode]->epsModel == VE_LEVEL_SET) {
      double pos_eps = ve[mode]->pos_ls.eps;
      double neg_eps = ve[mode]->eps;
      double width = ls->Length_Scale;
      err = level_set_property(neg_eps, pos_eps, width, &eps, d_eps_dF);
      GOMA_EH(err, "level_set_property() failed for ptt epsilon parameter.");
    } else {
      GOMA_EH(GOMA_ERROR, "Unknown PTT Epsilon parameter model");
    }

    /* get tensor dot products for future use */

    (void)tensor_dot(b, g, b_dot_g, VIM);

    dbl a_dot_b[DIM][DIM];
    dbl d_a_dot_b_db[DIM][DIM][DIM][DIM];
    dbl d_a_dot_b_dG[DIM][DIM][DIM][DIM];

    compute_a_dot_b(b, g, a_dot_b, d_a_dot_b_db, d_a_dot_b_dG);
    /*
     * Residuals_________________________________________________________________
     */

    if (af->Assemble_Residual) {
      /*
       * Assemble each component "ab" of the polymer stress equation...
       */
      for (int ii = 0; ii < VIM; ii++) {
        for (int jj = 0; jj < VIM; jj++) {

          if (ii <= jj) /* since the stress tensor is symmetric, only assemble the upper half */
          {
            eqn = R_s[mode][ii][jj];
            k = inv_v_s[ii][jj];

            /*
             * In the element, there will be contributions to this many equations
             * based on the number of degrees of freedom...
             */

            for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
              mass = 0.;

              if (pd->TimeIntegration != STEADY) {
                if (pd->e[pg->imtrx][eqn] & T_MASS) {
                  mass = b_dot[ii][jj];
                  mass *= at * lambda * det_J * wt;
                  mass *= h3;
                  mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                }
              }

              advection = 0.;
              if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                if (DOUBLE_NONZERO(lambda)) {

                  advection -= b_dot_g[ii][jj];
                  advection -= a_dot_b[ii][jj];
                  advection *= at * lambda * det_J * wt * h3;
                  advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                }
              }

              diffusion = 0.;
              if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                diffusion *= det_J * wt * h3;
                diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
              }

              /*
               * Source term...
               */

              source = 0.;
              if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                // consider whether saramitoCoeff should multiply here
                source += source_term[ii][jj];

                source *= det_J * h3 * wt;

                source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
              }

              /*
               * Add contributions to this residual (globally into Resid, and
               * locally into an accumulator)
               */

              func[mode][k] += mass + advection + diffusion + source;
            }
          }
        }
      }
    }

    /*
     * Jacobian terms...
     */

    if (af->Assemble_Jacobian) {
      dbl R_source, R_advection; /* Places to put the raw residual portions
                                    instead of constantly recalcing them */
      for (int ii = 0; ii < VIM; ii++) {
        for (int jj = 0; jj < VIM; jj++) {
          if (ii <= jj) /* since the stress tensor is symmetric, only assemble the upper half */
          {
            eqn = R_s[mode][ii][jj];
            k = inv_v_s[ii][jj];

            R_advection = v_dot_del_b[ii][jj] - x_dot_del_b[ii][jj];
            R_advection += b_dot_g[ii][jj] + a_dot_b[ii][jj];

            R_source = delta(ii, jj) - b[ii][jj];

            if (DOUBLE_NONZERO(alpha))
              R_source += alpha * lambda * (s_dot_s[ii][jj] / mup);
            R_source *= saramitoCoeff;
            R_source += -at * mup * (g[ii][jj] + gt[ii][jj]);

            for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
              int var;

              /*
               * J_S_T
               */

              var = TEMPERATURE;
              if (pd->v[pg->imtrx][var]) {
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  mass = 0.;

                  if (pd->TimeIntegration != STEADY) {
                    if (pd->e[pg->imtrx][eqn] & T_MASS) {
                      mass = b_dot[ii][jj];
                      mass *= d_at_dT[j] * lambda * det_J * wt;
                      mass *= h3;
                      mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                    }
                  }

                  advection = 0.;
                  if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                    if (DOUBLE_NONZERO(lambda)) {
                      advection -= (a_dot_b[ii][jj] + b_dot_g[ii][jj]);
                      advection *= d_at_dT[j] * lambda * det_J * wt * h3;
                      advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                    }
                  }

                  source = 0.;
                  source1 = 0.;
                  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                    source = -(g[ii][jj] + gt[ii][jj]) * (at * d_mup->T[j] + mup * d_at_dT[j]);

                    if (DOUBLE_NONZERO(alpha)) {
                      source1 -= s_dot_s[ii][jj] / (mup * mup) * d_mup->T[j];
                      source1 *= lambda * alpha * saramitoCoeff;
                      source += source1;
                    }
                    source *= det_J * wt * h3;
                    source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                  }

                  d_func[mode][k][var][j] += mass + advection + source;
                }
              }

              /*
               * J_S_v
               */
              for (p = 0; p < WIM; p++) {
                var = VELOCITY1 + p;
                if (pd->v[pg->imtrx][var]) {
                  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                    phi_j = bf[var]->phi[j];
                    d_mup_dv_pj = d_mup->v[p][j];

                    mass = 0.;

                    if (pd->TimeIntegration != STEADY) {
                      if (pd->e[pg->imtrx][eqn] & T_MASS) {

                        mass *=
                            pd->etm[pg->imtrx][eqn][(LOG2_MASS)] * at * lambda * det_J * wt * h3;
                      }
                    }

                    advection = 0.;
                    if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                      if (DOUBLE_NONZERO(lambda)) {
                        advection *= at * lambda * det_J * wt * h3;
                        advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                      }
                    }

                    diffusion = 0.;
                    if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                      diffusion *= det_J * wt * h3;
                      diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                    }

                    source = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                      source_c = -at * d_mup_dv_pj * (g[ii][jj] + gt[ii][jj]);
                      if (evss_gradv) {
                        if (pd->CoordinateSystem != CYLINDRICAL) {
                          source_c -= at * mup *
                                      (bf[VELOCITY1 + ii]->grad_phi_e[j][p][ii][jj] +
                                       bf[VELOCITY1 + jj]->grad_phi_e[j][p][jj][ii]);
                        } else {
                          source_c -= at * mup *
                                      (bf[VELOCITY1]->grad_phi_e[j][p][ii][jj] +
                                       bf[VELOCITY1]->grad_phi_e[j][p][jj][ii]);
                        }
                      }

                      source_a = 0.;
                      if (DOUBLE_NONZERO(alpha)) {
                        source_a = -s_dot_s[ii][jj] / (mup * mup);
                        source_a *= saramitoCoeff * alpha * lambda * d_mup_dv_pj;
                      }

                      source_b = 0.;
                      source = source_a + source_b + source_c;
                      source *= det_J * wt * h3;
                      source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                    }

                    d_func[mode][k][var][j] += mass + advection + diffusion + source;
                  }
                }
              }

              /*
               * J_S_c
               */
              var = MASS_FRACTION;
              if (pd->v[pg->imtrx][var]) {
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  for (w = 0; w < pd->Num_Species_Eqn; w++) {

                    source = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                      source_a = -at * d_mup->C[w][j] * (g[ii][jj] + gt[ii][jj]);

                      source_b = 0.;
                      if (DOUBLE_NONZERO(alpha)) {
                        source_b -= s_dot_s[ii][jj] / (mup * mup);
                        source_b *= alpha * lambda * saramitoCoeff * d_mup->C[w][j];
                      }
                      source = source_a + source_b;
                      source *= det_J * wt * h3;
                      source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                    }

                    if (w > 1) {
                      GOMA_EH(GOMA_ERROR, "Need more arrays for each species.");
                    }

                    d_func[mode][k][var][j] += source;
                  }
                }
              }

              /*
               * J_S_P
               */
              var = PRESSURE;
              if (pd->v[pg->imtrx][var]) {
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  source = 0.;
                  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                    source_a += -at * d_mup->P[j] * (g[ii][jj] + gt[ii][jj]);

                    source_b = 0.;
                    if (DOUBLE_NONZERO(alpha)) {
                      source_b -= (s_dot_s[ii][jj] / (mup * mup));
                      source_b *= d_mup->P[j] * alpha * lambda * saramitoCoeff;
                    }
                    source = source_a + source_b;
                    source *= det_J * wt * h3;
                    source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                  }

                  d_func[mode][k][var][j] += source;
                }
              }

              /*
               * J_S_d
               */
              for (p = 0; p < dim; p++) {
                var = MESH_DISPLACEMENT1 + p;
                if (pd->v[pg->imtrx][var]) {
                  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                    phi_j = bf[var]->phi[j];
                    d_det_J_dmesh_pj = bf[eqn]->d_det_J_dm[p][j];
                    dh3dmesh_pj = fv->dh3dq[p] * bf[var]->phi[j];
                    d_mup_dmesh_pj = d_mup->X[p][j];

                    mass = 0.;
                    mass_a = 0.;
                    mass_b = 0.;
                    if (pd->TimeIntegration != STEADY) {
                      if (pd->e[pg->imtrx][eqn] & T_MASS) {
                        mass_a = b_dot[ii][jj];
                        mass_a *= (d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj);

                        mass = mass_a + mass_b;
                        mass *= at * lambda * wt * pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                      }
                    }

                    advection = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                      if (DOUBLE_NONZERO(lambda)) {
                        /*
                         * Four parts:
                         *    advection_a =
                         *    	Int ( ea.(v-xdot).d(Vv)/dmesh h3 |Jv| )
                         *
                         *    advection_b =
                         *  (i)	Int ( ea.(v-xdot).Vv h3 d(|Jv|)/dmesh )
                         *  (ii)  Int ( ea.(v-xdot).d(Vv)/dmesh h3 |Jv| )
                         *  (iii) Int ( ea.(v-xdot).Vv dh3/dmesh |Jv|   )
                         *
                         * For unsteady problems, we have an
                         * additional term
                         *
                         *    advection_c =
                         *    	Int ( ea.d(v-xdot)/dmesh.Vv h3 |Jv| )
                         */

                        advection_a = R_advection;

                        advection_a *= (d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj);

                        d_vdotdels_dm = 0.;

                        advection_b = d_vdotdels_dm;
                        advection_b *= det_J * h3;

                        advection_c = 0.;
                        if (pd->TimeIntegration != STEADY) {
                          if (pd->e[pg->imtrx][eqn] & T_MASS) {
                            d_xdotdels_dm = (1. + 2. * tt) * phi_j / dt * grad_b[p][ii][jj];

                            advection_c -= d_xdotdels_dm;

                            advection_c *= h3 * det_J;
                          }
                        }

                        advection_d = 0.;

                        advection = advection_a + advection_b + advection_c + advection_d;

                        advection *= wt * at * lambda * pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                      }
                    }

                    diffusion = 0.;
                    if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                      diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                    }

                    /*
                     * Source term...
                     */

                    source = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                      source_a = R_source;
                      source_b = -at * (g[ii][jj] + gt[ii][jj]);

                      if (DOUBLE_NONZERO(alpha)) {
                        source_b += -s_dot_s[ii][jj] / (mup * mup) * alpha * lambda * saramitoCoeff;
                      }

                      source_a *= (d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj);

                      source_b *= det_J * h3 * d_mup_dmesh_pj;

                      source_c = 0.;

                      source = source_a + source_b + source_c;

                      source *= wt * pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                    }

                    d_func[mode][k][var][j] += mass + advection + diffusion + source;
                  }
                }
              }

              /*
               * J_S_G
               */
              if (evss_gradv == 0) {
                for (p = 0; p < VIM; p++) {
                  for (q = 0; q < VIM; q++) {
                    var = v_g[p][q];

                    if (pd->v[pg->imtrx][var]) {
                      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                        phi_j = bf[var]->phi[j];
                        advection = 0.;
                        advection_a = 0.;
                        if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                          if (DOUBLE_NONZERO(lambda)) {
                            for (int k = 0; k < VIM; k++) {
                              advection += -b[ii][k] * delta(p, k) * delta(jj, q);
                            }
                            advection += -d_a_dot_b_dG[p][q][ii][jj];
                            advection *= phi_j * h3 * det_J;

                            advection *=
                                wt * at * lambda * pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                          }
                        }

                        /*
                         * Diffusion...
                         */

                        diffusion = 0.;

                        if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                          diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                        }

                        /*
                         * Source term...
                         */

                        source = 0.;

                        d_func[mode][k][var][j] += advection + diffusion + source;
                      }
                    }
                  }
                }
              }

              /*
               * J_S_F
               */
              var = FILL;
              if (pd->v[pg->imtrx][var]) {
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  mass = 0.;

                  if (pd->TimeIntegration != STEADY) {
                    if (pd->e[pg->imtrx][eqn] & T_MASS) {

                      mass = b_dot[ii][jj];
                      mass *= d_lambda_dF[j];
                      mass *= at * det_J * wt;
                      mass *= h3;
                      mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                    }
                  }

                  advection = 0.;

                  if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                    if (d_lambda_dF[j] != 0.) {

                      advection += v_dot_del_b[ii][jj] - x_dot_del_b[ii][jj];

                      advection *= d_lambda_dF[j];
                      advection *= at * det_J * wt * h3;
                      advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                    }
                  }

                  diffusion = 0.;

                  if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                    /* add SU term in here when appropriate */

                    diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                  }

                  source = 0.;

                  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {

                    double invmup = 1 / mup;
                    // PTT

                    // Giesekus
                    if (alpha != 0.) {
                      source += s_dot_s[ii][jj] *
                                (-alpha * lambda * d_mup->F[j] * invmup * invmup +
                                 d_alpha_dF[j] * lambda * invmup + alpha * d_lambda_dF[j] * invmup);
                    }

                    source *= det_J * h3 * wt;

                    source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                  }

                  d_func[mode][k][var][j] += mass + advection + diffusion + source;
                }
              }

              /*
               * J_S_S
               */
              for (p = 0; p < VIM; p++) {
                for (q = 0; q < VIM; q++) {
                  if (q >= p) {
                    var = v_s[mode][p][q];

                    if (pd->v[pg->imtrx][var]) {
                      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                        phi_j = bf[var]->phi[j];
                        mass = 0.;
                        if (pd->TimeIntegration != STEADY) {
                          if (pd->e[pg->imtrx][eqn] & T_MASS) {
                            mass = (1. + 2. * tt) * phi_j / dt * (double)delta(ii, p) *
                                   (double)delta(jj, q);
                            mass *= h3 * det_J;
                            mass *= at * lambda * wt * pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                          }
                        }

                        advection = 0.;

                        if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                          if (DOUBLE_NONZERO(lambda)) {

                            for (int k = 0; k < VIM; k++) {
                              advection -=
                                  phi_j *
                                  (delta(ii, q) * delta(k, p) | delta(ii, p) * delta(k, q)) *
                                  g[k][jj];
                            }
                            advection -= phi_j * d_a_dot_b_db[p][q][ii][jj];

                            advection *= h3 * det_J;

                            advection *=
                                wt * at * lambda * pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                          }
                        }

                        /*
                         * Diffusion...
                         */

                        diffusion = 0.;

                        if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                          diffusion *= det_J * wt * h3;
                          diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                        }

                        /*
                         * Source term...
                         */

                        source = 0.;

                        if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                          source = d_source_term_db[ii][jj][p][q];
                          source *=
                              phi_j * det_J * h3 * wt * pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                        }

                        d_func[mode][k][var][j] += mass + advection + diffusion + source;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    } /* End Assemble Jacobian */
  }   /* End loop over modes */
}

void stress_no_v_dot_gradS_logc_transient(
    double func[MAX_MODES][6],
    double d_func[MAX_MODES][6][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
    const double dt,
    const double tt)

/****************************************************************************
 *
 *  Function which imposes fully-developed fluid stresses
 *  in the velocity direction
 *
 *  It's a clone for assemble_stress_fortin with the v dot grad S term removed
 *
 *  For more details see:
 *
 *    Xueying Xie and Matteo Pasquali J. Non-Newtonian Fluid Mech. 122 (2004) 159 - 176
 *
 *  Kristianto Tjiptowidjojo (tjiptowi@unm.edu)
 *
 *  This routine is adjusted for the log-conformation tensor
 *****************************************************************************/
{

  int i, j, eqn, a, b, mode, w, k;
  int siz;
  int R_s[MAX_MODES][DIM][DIM];
  int v_s[MAX_MODES][DIM][DIM];
  int inv_v_s[DIM][DIM];
  int logc_gradv = 0;

  dbl grad_v[DIM][DIM]; /* Velocity gradient based on velocity - discontinuous across element */
  dbl gamma[DIM][DIM];  /* Shear-rate tensor based on velocity */

  dbl exp_s[DIM][DIM]; // Exponential of log_conf
  dbl trace = 0.0;     /* trace of the stress tensor */
  dbl s_dot[DIM][DIM]; /* stress tensor from last time step */
  dbl gt[DIM][DIM];    /* transpose of velocity gradient tensor */

  dbl source_term1[DIM][DIM];
  dbl advection_term1[DIM][DIM];

  /* polymer viscosity and derivatives */
  dbl mup;
  VISCOSITY_DEPENDENCE_STRUCT d_mup_struct;
  VISCOSITY_DEPENDENCE_STRUCT *d_mup = &d_mup_struct;

  /*  shift function */
  dbl at = 0.0;
  dbl wlf_denom;

  // Decomposition of velocity vector
  dbl eig_values[DIM];
  dbl R1[DIM][DIM];
  dbl R1_T[DIM][DIM];
  dbl Rt_dot_gradv[DIM][DIM];
  dbl D[DIM][DIM];
  dbl D_dot_D[DIM][DIM];
  dbl M1[DIM][DIM];

  dbl tmp1[DIM][DIM], tmp2[DIM][DIM], tmp3[DIM][DIM];

  /* constitutive equation parameters */
  dbl alpha;      /* This is the Geisekus mobility parameter */
  dbl lambda = 0; /* polymer relaxation constant */
  dbl eps;        /* This is the PTT elongation parameter */
  dbl Z = 1.0;    /* This is the factor appearing in front of the stress tensor in PTT */

  /* ETMs*/
  dbl mass, advection, source;

  dbl d_lambda;

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  /* A bookeeping array from stress[i][j] into its position in func */
  inv_v_s[0][0] = 0; /* S11 */
  inv_v_s[0][1] = 1; /* S12 */
  inv_v_s[1][1] = 2; /* S22 */
  inv_v_s[0][2] = 3; /* S13 */
  inv_v_s[1][2] = 4; /* S23 */
  inv_v_s[2][2] = 5; /* S33 */

  (void)stress_eqn_pointer(v_s);
  (void)stress_eqn_pointer(R_s);

  memset(exp_s, 0, sizeof(double) * DIM * DIM);
  memset(R1, 0, sizeof(double) * DIM * DIM);
  memset(eig_values, 0, sizeof(double) * DIM);

  if (vn->evssModel == LOG_CONF_GRADV) {
    logc_gradv = 1;
  }

  /*
   * Load up Field variables...
   */

  /****  Velocity gradient, and shear rate ****/
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      grad_v[a][b] = fv->grad_v[a][b];
    }
  }
  /* load up shear rate tensor based on velocity */
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      gamma[a][b] = grad_v[a][b] + grad_v[b][a];
    }
  }
  /* Continuous velocity gradient field and its transpose */
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      gt[b][a] = fv->G[a][b];
    }
  }

  /*  Time-temperature shift factor  */
  if (pd->e[pg->imtrx][TEMPERATURE]) {
    if (vn->shiftModel == CONSTANT) {
      at = vn->shift[0];
    } else if (vn->shiftModel == MODIFIED_WLF) {
      wlf_denom = vn->shift[1] + fv->T - mp->reference[TEMPERATURE];
      if (wlf_denom != 0.) {
        at = exp(vn->shift[0] * (mp->reference[TEMPERATURE] - fv->T) / wlf_denom);
      } else {
        at = 1.;
      }
    }
  } else {
    at = 1.;
  }

  /* Begin loop over modes */
  for (mode = 0; mode < vn->modes; mode++) {

    /*
     * Load polymeric stress tensor for each mode
     */
    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        if (pd->TimeIntegration != STEADY) {
          s_dot[a][b] = fv_dot->S[mode][a][b];
        } else {
          s_dot[a][b] = 0.;
        }
      }
    }

    // Polymer viscosity
    mup = viscosity(ve[mode]->gn, gamma, d_mup);

    // Giesekus mobility parameter
    alpha = ve[mode]->alpha;

    // Polymer time constant
    if (ve[mode]->time_constModel == CONSTANT) {
      lambda = ve[mode]->time_const;
    } else if (ve[mode]->time_constModel == CARREAU || ve[mode]->time_constModel == POWER_LAW) {
      lambda = mup / ve[mode]->time_const;
    }
    const bool saramitoEnabled =
        (vn->ConstitutiveEquation == SARAMITO_OLDROYDB ||
         vn->ConstitutiveEquation == SARAMITO_PTT || vn->ConstitutiveEquation == SARAMITO_GIESEKUS);

    dbl saramitoCoeff = 1.;

#ifdef ANALEIG_PLEASE
    analytical_exp_s(fv_old->S[mode], exp_s, eig_values, R1, NULL);
#else
    compute_exp_s(fv_old->S[mode], exp_s, eig_values, R1);
#endif

    // Decompose velocity gradient

    memset(D, 0, sizeof(double) * DIM * DIM);
    D[0][0] = eig_values[0];
    D[1][1] = eig_values[1];
    if (VIM > 2) {
      D[2][2] = eig_values[2];
    }
    (void)tensor_dot(D, D, D_dot_D, VIM);

    // Decompose velocity gradient

    memset(M1, 0, sizeof(double) * DIM * DIM);
    memset(R1_T, 0, sizeof(double) * DIM * DIM);

    for (i = 0; i < VIM; i++) {
      for (j = 0; j < VIM; j++) {
        R1_T[i][j] = R1[j][i];
      }
    }

    for (i = 0; i < VIM; i++) {
      for (j = 0; j < VIM; j++) {
        Rt_dot_gradv[i][j] = 0.;
        for (w = 0; w < VIM; w++) {
          if (logc_gradv) {
            Rt_dot_gradv[i][j] += R1_T[i][w] * grad_v[j][w];
          } else {
            Rt_dot_gradv[i][j] += R1_T[i][w] * gt[w][j];
          }
        }
      }
    }

    for (i = 0; i < VIM; i++) {
      for (j = 0; j < VIM; j++) {
        M1[i][j] = 0.;
        for (w = 0; w < VIM; w++) {
          M1[i][j] += Rt_dot_gradv[i][w] * R1[w][j];
        }
      }
    }

    dbl tau[DIM][DIM] = {{0.0}};
    if (saramitoEnabled == TRUE) {
      for (int i = 0; i < VIM; i++) {
        for (int j = 0; j < VIM; j++) {
          tau[i][j] = mup / lambda * (exp_s[i][j] - delta(i, j));
        }
      }
      compute_saramito_model_terms(&saramitoCoeff, NULL, tau, ve[mode]->gn, FALSE);
    } else {
      saramitoCoeff = 1.;
    }

    // Predetermine advective terms
    trace = eig_values[0] + eig_values[1];
    if (VIM > 2) {
      trace += eig_values[2];
    }

    // PTT exponent
    eps = ve[mode]->eps;

    // PTT
    Z = 1;
    if (vn->ConstitutiveEquation == PTT) {
      if (vn->ptt_type == PTT_LINEAR) {
        Z = 1 + eps * (trace - (double)VIM);
      } else if (vn->ptt_type == PTT_EXPONENTIAL) {
        Z = exp(eps * (trace - (double)VIM));
      } else {
        GOMA_EH(GOMA_ERROR, "Unrecognized PTT Form %d", vn->ptt_type);
      }
    }

    siz = sizeof(double) * DIM * DIM;
    memset(tmp1, 0, siz);
    memset(tmp2, 0, siz);
    memset(tmp3, 0, siz);
    memset(advection_term1, 0, siz);
    memset(source_term1, 0, siz);

    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        if (a != b) {
          d_lambda = eig_values[b] - eig_values[a];
          if (fabs(d_lambda) > 1.e-8) {
            tmp1[a][b] += (log(eig_values[b]) - log(eig_values[a])) / d_lambda;
            tmp1[a][b] *= (eig_values[a] * M1[b][a] + eig_values[b] * M1[a][b]);
          } else {
            tmp1[a][b] += M1[a][b] + M1[b][a];
          }
        }
        if (a == b) {
          source_term1[a][b] += saramitoCoeff * Z * (1.0 - D[a][a]) / lambda;
          if (alpha != 0) {
            source_term1[a][b] += alpha * (2.0 * D[a][a] - 1.0 - D_dot_D[a][a]) / lambda;
          }
          source_term1[a][b] /= eig_values[a];
          source_term1[a][b] += 2.0 * M1[a][a];
        }
      }
    }

    (void)tensor_dot(R1, tmp1, tmp2, VIM);
    (void)tensor_dot(tmp2, R1_T, advection_term1, VIM);
    (void)tensor_dot(R1, source_term1, tmp3, VIM);
    (void)tensor_dot(tmp3, R1_T, source_term1, VIM);

    /**** Assemble func *****/

    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {

        if (a <= b) /* since the stress tensor is symmetric, only assemble the upper half */
        {
          eqn = R_s[mode][a][b];
          k = inv_v_s[a][b];

          /* Mass term */
          mass = 0.;
          if (pd->TimeIntegration != STEADY) {
            mass = at * s_dot[a][b];
            mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
          }

          /* Advection term - minus v_dot_grad_S and x_dot_del_S terms */
          advection = 0.;
          if (lambda != 0.) {
            advection -= advection_term1[a][b];
            advection *= at * pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
          }

          /* Source term */
          source = 0.;
          source -= source_term1[a][b];
          source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

          func[mode][k] += mass + advection + source;
        }
      }
    }
    if (af->Assemble_Jacobian) {
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          if (a <= b) /* since the stress tensor is symmetric, only assemble the upper half */
          {
            eqn = R_s[mode][a][b];
            k = inv_v_s[a][b];

            /* Sensitivities w.r.t. polymeric stress - J_S_S */
            int var = POLYMER_STRESS11;
            if (pd->v[pg->imtrx][var]) {
              for (int p = 0; p < VIM; p++) {
                for (int q = 0; q < VIM; q++) {
                  var = v_s[mode][p][q];

                  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                    double phi_j = bf[var]->phi[j];

                    /* mass term */
                    mass = 0.0;
                    if (pd->TimeIntegration != STEADY) {
                      mass =
                          (1. + 2. * tt) * phi_j / dt * (double)delta(a, p) * (double)delta(b, q);
                      mass *= at;
                    }
                    mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];

                    /* Load them up */
                    d_func[mode][k][var][j] += mass;
                  }
                }
              }
            } /* End of J_S_S */
          }   /* End of if a == b */
        }     /* End of loop over dimension "b" */
      }       /* End of loop over dimension "a" */
    }         /*  end of if Assemble_Jacobian  */
  }           /* End of loop over modes */
} /* END of routine stress_no_v_dot_gradS_logc_transient */

void PSPG_consistency_bc(double *func,
                         double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                         const dbl x_dot[DIM],
                         const dbl time, /* current time  */
                         const dbl dt,   /* time step size */
                         const dbl tt,   /* time step parameter */
                         const PG_DATA *pg_data)
/******************************************************************************
 *
 *  Function which calculates the missing pressure stabilization terms
 *    for PSPG and adds it to the continuity equation
 *
 ******************************************************************************/

{
  /* Local variables */

  int j, var, a, b, q, p;

  /*
   * Variables for Pressure Stabilization Petrov-Galerkin...
   */

  const int dim = pd->Num_Dim;
  int mode;
  int meqn;
  int var1;
  int r, s, t;
  int w = 0;
  int v_s[MAX_MODES][DIM][DIM];
  const int v_g[DIM][DIM] = {{VELOCITY_GRADIENT11, VELOCITY_GRADIENT12, VELOCITY_GRADIENT13},
                             {VELOCITY_GRADIENT21, VELOCITY_GRADIENT22, VELOCITY_GRADIENT23},
                             {VELOCITY_GRADIENT31, VELOCITY_GRADIENT32, VELOCITY_GRADIENT33}};

  dbl mass;
  dbl source_a;
  dbl diffusion;
  dbl advection_a;
  dbl stress;
  dbl pressure;
  dbl velocity_gradient;
  dbl pressure_stabilization;
  dbl momentum_residual[DIM]; /* momentum residual for PSPG */
  dbl v_dot[DIM];
  dbl div_s[DIM];
  dbl div_G[DIM];

  dbl h_elem;
  dbl rho;
  int pspg_local = 0;
  int pspg_global = 0;
  int pspg_shakib = 0;
  dbl tau_pspg = 0;
  dbl d_tau_pspg_dv[DIM][MDE];
  dbl d_tau_pspg_dX[DIM][MDE];

  dbl f[DIM];                                  /* Body force. */
  MOMENTUM_SOURCE_DEPENDENCE_STRUCT df_struct; /* Body force dependence */
  MOMENTUM_SOURCE_DEPENDENCE_STRUCT *df = &df_struct;

  /*
   * Variables for vicosity and derivative
   */
  dbl mu;
  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct; /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;

  dbl gamma[DIM][DIM]; /* shrearrate tensor based on velocity */

  dbl div_tau_p[DIM];                     /* divergence of the particle stress*/
  dbl d_div_tau_p_dgd[DIM][MDE];          /* derivative wrt shear_rate_invariant */
  dbl d_div_tau_p_dy[DIM][MAX_CONC][MDE]; /* derivative wrt concentration */
  dbl d_div_tau_p_dv[DIM][DIM][MDE];      /* derivative wrt velocity */
  dbl d_div_tau_p_dmesh[DIM][DIM][MDE];   /* derivative wrt mesh */
  dbl d_div_tau_p_dp[DIM][MDE];           /* derivative wrt pressure dir */

  stress_eqn_pointer(v_s);

  double h_elem_avg = pg_data->h_elem_avg;
  double mu_avg = pg_data->mu_avg;
  double U_norm = pg_data->U_norm;

  /* This is the flag for the standard global PSPG */
  if (PSPG == 1) {
    pspg_global = TRUE;
    pspg_local = FALSE;
  }
  /* This is the flag for the standard local PSPG */
  else if (PSPG == 2) {
    pspg_global = FALSE;
    pspg_local = TRUE;
  } else if (PSPG == 3) {
    pspg_global = FALSE;
    pspg_local = FALSE;
    pspg_shakib = TRUE;
  } else {
    return;
  }

  /*  h_elem = 0.;
    for ( p=0; p<dim; p++)
      {
        h_elem += h[p];
      }
    h_elem = sqrt(h_elem)/2.; */

  h_elem = h_elem_avg;

  /* Calculate a simple arithmetic average viscosity and density
     in the element                                              */

  rho = mp->density;

  /* Now calculate the element Reynolds number based on a global
     norm of the velocity */

  memset(d_tau_pspg_dv, 0, sizeof(double) * DIM * MDE);
  memset(d_tau_pspg_dX, 0, sizeof(double) * DIM * MDE);

  if (pspg_global) {

    /* Now calculate the element Reynolds number based on a global
     * norm of the velocity and determine tau_pspg discretely from Re
     * The global version has no Jacobian dependencies
     */
    double Re = rho * U_norm * h_elem / (2.0 * mu_avg);

    if (Re <= 3.0) {
      tau_pspg = PS_scaling * h_elem * h_elem / (12.0 * mu_avg);
    } else if (Re > 3.0) {
      tau_pspg = PS_scaling * h_elem / (2.0 * rho * U_norm);
    }
  } else if (pspg_local) {
    double hh_siz = 0.;
    for (p = 0; p < dim; p++) {
      hh_siz += pg_data->hsquared[p];
    }
    // Average value of h**2 in the element
    hh_siz = hh_siz / ((double)dim);

    // Average value of v**2 in the element
    double vv_speed = 0.0;
    for (a = 0; a < WIM; a++) {
      vv_speed += pg_data->v_avg[a] * pg_data->v_avg[a];
    }

    double rho_avg = pg_data->rho_avg;
    double mu_avg = pg_data->mu_avg;

    // Use vv_speed and hh_siz for tau_pspg, note it has a continuous dependence on Re
    double tau_pspg1 =
        rho_avg * rho_avg * vv_speed / hh_siz + (9.0 * mu_avg * mu_avg) / (hh_siz * hh_siz);
    if (pd->TimeIntegration != STEADY) {
      tau_pspg1 += 4.0 / (dt * dt);
    }
    tau_pspg = PS_scaling / sqrt(tau_pspg1);

    // tau_pspg derivatives wrt v from vv_speed
    if (pd->v[pg->imtrx][VELOCITY1]) {
      for (b = 0; b < dim; b++) {
        var = VELOCITY1 + b;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_tau_pspg_dv[b][j] = -tau_pspg / tau_pspg1;
            d_tau_pspg_dv[b][j] *=
                rho_avg * rho_avg / hh_siz * pg_data->v_avg[b] * pg_data->dv_dnode[b][j];
          }
        }
      }
    }

    // tau_pspg derivatives wrt mesh from hh_siz
    if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
      for (b = 0; b < dim; b++) {
        var = MESH_DISPLACEMENT1 + b;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_tau_pspg_dX[b][j] = tau_pspg / tau_pspg1;
            d_tau_pspg_dX[b][j] *=
                (rho_avg * rho_avg * vv_speed + 18.0 * mu_avg * mu_avg / hh_siz) /
                (hh_siz * hh_siz);
            d_tau_pspg_dX[b][j] *= pg_data->hhv[b][b] * pg_data->dhv_dxnode[b][j] / ((double)dim);
          }
        }
      }
    }
  } else if (pspg_shakib) {

    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
      }
    }

    mu = viscosity(gn, gamma, d_mu);

    dbl mup[MAX_MODES];
    VISCOSITY_DEPENDENCE_STRUCT d_mup[MAX_MODES]; /* viscosity dependence */
    if (pd->gv[POLYMER_STRESS11]) {
      for (int mode = 0; mode < vn->modes; mode++) {
        mup[mode] = viscosity(ve[mode]->gn, gamma, &d_mup[mode]);
      }
    }

    dbl G[DIM][DIM];
    get_metric_tensor(bf[pd->ShapeVar]->B, pd->Num_Dim, ei[pg->imtrx]->ielem_type, G);

    dbl tau_time = 0;
    // time term
    if (pd->TimeIntegration != STEADY) {
      tau_time += 4 * rho * rho / (dt * dt);
    }

    // advection
    dbl tau_adv = 0;
    for (int i = 0; i < dim; i++) {
      for (int j = 0; j < dim; j++) {
        tau_adv += rho * rho * fv->v[i] * G[i][j] * fv->v[j];
      }
    }

    // diffusion
    dbl tau_diff = 0;
    dbl mu_total = mu;
    for (int mode = 0; mode < vn->modes; mode++) {
      mu_total += mup[mode];
    }
    dbl coeff = 12 * (mu_total * mu_total);
    for (int i = 0; i < dim; i++) {
      for (int j = 0; j < dim; j++) {
        tau_diff += coeff * G[i][j] * G[i][j];
      }
    }

    tau_pspg = PS_scaling / sqrt(tau_time + tau_adv + tau_diff);

    // d/dx 1/sqrt(f(x)) => - f'(x) / (2 * f(x)^(3/2))
    if (pd->v[pg->imtrx][VELOCITY1]) {
      for (b = 0; b < dim; b++) {
        var = VELOCITY1 + b;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            dbl tau_adv_dv = 0;
            for (int a = 0; a < dim; a++) {
              tau_adv_dv += rho * rho * bf[var]->phi[j] * G[b][a] * fv->v[a];
              tau_adv_dv += rho * rho * bf[var]->phi[j] * G[a][b] * fv->v[a];
            }

            d_tau_pspg_dv[b][j] = -rho * tau_pspg * tau_pspg * tau_pspg * tau_adv_dv;
          }
        }
      }
    }
  }

  /* load up shearrate tensor based on velocity */
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
    }
  }

  /* get viscosity for velocity second derivative/diffusion term in PSPG stuff */
  mu = viscosity(gn, gamma, d_mu);

  for (a = 0; a < dim; a++) {
    if (pd->TimeIntegration != STEADY) {
      v_dot[a] = fv_dot->v[a];
    } else {
      v_dot[a] = 0.0;
    }
  }

  for (p = 0; p < WIM; p++) {
    div_s[p] = 0.;
  }

  if (pd->gv[POLYMER_STRESS11]) {
    for (p = 0; p < WIM; p++) {
      for (mode = 0; mode < vn->modes; mode++) {
        div_s[p] += fv->div_S[mode][p];
      }
    }
  }

  if (pd->gv[VELOCITY_GRADIENT11]) {
    for (p = 0; p < dim; p++) {
      div_G[p] = fv->div_G[p];
    }
  } else {
    for (p = 0; p < dim; p++) {
      div_G[p] = 0.;
    }
  }

  memset(div_tau_p, 0, sizeof(double) * DIM);
  memset(d_div_tau_p_dgd, 0, sizeof(double) * DIM * MDE);
  memset(d_div_tau_p_dy, 0, sizeof(double) * DIM * MAX_CONC * MDE);
  memset(d_div_tau_p_dv, 0, sizeof(double) * DIM * DIM * MDE);
  memset(d_div_tau_p_dmesh, 0, sizeof(double) * DIM * DIM * MDE);
  memset(d_div_tau_p_dp, 0, sizeof(double) * DIM * MDE);
  if (cr->MassFluxModel == DM_SUSPENSION_BALANCE && PSPG) {
    /* This is the divergence of the particle stress  */
    divergence_particle_stress(div_tau_p, d_div_tau_p_dgd, d_div_tau_p_dy, d_div_tau_p_dv,
                               d_div_tau_p_dmesh, d_div_tau_p_dp, w);
  }

  /* get momentum source term */
  momentum_source_term(f, df, time);

  for (a = 0; a < dim; a++) {
    meqn = R_MOMENTUM1 + a;
    momentum_residual[a] = rho * v_dot[a] * pd->etm[pg->imtrx][meqn][(LOG2_ADVECTION)] +
                           fv->grad_P[a] * pd->etm[pg->imtrx][meqn][(LOG2_DIFFUSION)] -
                           div_s[a] * pd->etm[pg->imtrx][meqn][(LOG2_DIFFUSION)] -
                           div_tau_p[a] * pd->etm[pg->imtrx][meqn][(LOG2_DIFFUSION)] -
                           mu * div_G[a] * pd->etm[pg->imtrx][meqn][(LOG2_DIFFUSION)] -
                           f[a] * pd->etm[pg->imtrx][meqn][(LOG2_SOURCE)];
    for (b = 0; b < dim; b++) {
      momentum_residual[a] += rho * (fv->v[b] - x_dot[b]) * fv->grad_v[b][a] *
                              pd->etm[pg->imtrx][meqn][(LOG2_ADVECTION)];
    }
  }

  /*  calculate pressure at current gauss point */
  for (p = 0; p < pd->Num_Dim; p++) {
    *func += momentum_residual[p] * fv->snormal[p];
  }
  *func *= tau_pspg;

  if (af->Assemble_LSA_Mass_Matrix) {
    for (a = 0; a < dim; a++) {
      var = VELOCITY1 + a;
      meqn = R_MOMENTUM1 + a;
      if (pd->v[pg->imtrx][var]) {
        if (pd->e[pg->imtrx][meqn] & T_MASS) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_func[0][var][j] += tau_pspg * bf[var]->phi[j] * rho *
                                 pd->etm[pg->imtrx][meqn][(LOG2_MASS)] * fv->snormal[a];
          }
        }
      }

      var = MESH_DISPLACEMENT1 + a;
      if (pd->v[pg->imtrx][var])
        if (pd->e[pg->imtrx][meqn] & T_ADVECTION)
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
            for (b = 0; b < pd->Num_Dim; b++) {
              meqn = R_MOMENTUM1 + b;
              d_func[0][var][j] -= rho * bf[var]->phi[j] * fv->grad_v[a][b] *
                                   pd->etm[pg->imtrx][meqn][(LOG2_ADVECTION)] * fv->snormal[b] *
                                   tau_pspg;
            }
    }
    return;
  }

  if (af->Assemble_Jacobian) {
    for (b = 0; b < dim; b++) {
      var = VELOCITY1 + b;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          pressure_stabilization = 0.;
          for (a = 0; a < dim; a++) {
            meqn = R_MOMENTUM1 + a;

            mass = 0.;
            if (pd->TimeIntegration != STEADY) {
              if (pd->e[pg->imtrx][meqn] & T_MASS) {
                mass += (1. + 2. * tt) * bf[var]->phi[j] / dt * delta(a, b);
                mass *= rho * pd->etm[pg->imtrx][meqn][(LOG2_MASS)];
              }
            }

            diffusion = 0.;
            if (pd->e[pg->imtrx][meqn] & T_DIFFUSION) {
              diffusion -= d_mu->v[b][j] * div_G[a];
              diffusion -= d_div_tau_p_dv[a][b][j];
              diffusion *= pd->etm[pg->imtrx][meqn][(LOG2_DIFFUSION)];
            }

            advection_a = 0.;
            if (pd->e[pg->imtrx][meqn] & T_ADVECTION) {
              advection_a += bf[var]->phi[j] * fv->grad_v[b][a];
              for (p = 0; p < dim; p++) {
                advection_a += (fv->v[p] - x_dot[p]) * bf[var]->grad_phi_e[j][b][p][a];
              }
              advection_a *= rho * pd->etm[pg->imtrx][meqn][(LOG2_ADVECTION)];
            }

            source_a = 0.;
            if (pd->e[pg->imtrx][meqn] & T_SOURCE) {
              source_a -= df->v[a][b][j] * pd->etm[pg->imtrx][meqn][(LOG2_SOURCE)];
            }

            pressure_stabilization +=
                tau_pspg * (mass + diffusion + advection_a + source_a) * fv->snormal[a];
          }

          d_func[0][var][j] +=
              pressure_stabilization + d_tau_pspg_dv[b][j] * momentum_residual[a] * fv->snormal[a];
        }
      }
    }
    var = MASS_FRACTION;
    if (pd->v[pg->imtrx][var]) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          pressure_stabilization = 0.;
          for (a = 0; a < dim; a++) {
            meqn = R_MOMENTUM1 + a;
            if (pd->e[pg->imtrx][meqn] & T_SOURCE) {
              pressure_stabilization = -tau_pspg * fv->snormal[a] * df->C[a][w][j] *
                                       pd->etm[pg->imtrx][meqn][(LOG2_SOURCE)];
            }

            if (pd->e[pg->imtrx][meqn] & T_DIFFUSION) {
              pressure_stabilization -= tau_pspg * fv->snormal[a] *
                                        (d_mu->C[w][j] * div_G[a] + d_div_tau_p_dy[a][w][j]) *
                                        pd->etm[pg->imtrx][meqn][(LOG2_DIFFUSION)];
            }
          }

          d_func[0][MAX_VARIABLE_TYPES + w][j] += pressure_stabilization;
        }
      }
    }

    var = PRESSURE;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        pressure_stabilization = 0.;
        for (a = 0; a < dim; a++) {
          meqn = R_MOMENTUM1 + a;
          if (pd->e[pg->imtrx][meqn] & T_DIFFUSION) {
            pressure_stabilization += tau_pspg * fv->snormal[a] * bf[var]->grad_phi[j][a] *
                                      pd->etm[pg->imtrx][meqn][(LOG2_DIFFUSION)];
          }
        }
        d_func[0][var][j] += pressure_stabilization;
      }
    }

    var = SHEAR_RATE;
    if (cr->MassFluxModel == DM_SUSPENSION_BALANCE) {
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          pressure_stabilization = 0.;
          for (a = 0; a < dim; a++) {
            meqn = R_MOMENTUM1 + a;
            if (pd->e[pg->imtrx][meqn] & T_DIFFUSION) {
              pressure_stabilization += tau_pspg * fv->snormal[a] * d_div_tau_p_dgd[a][j] *
                                        pd->etm[pg->imtrx][meqn][(LOG2_DIFFUSION)];
            }
          }
          d_func[0][var][j] += pressure_stabilization;
        }
      }
    }

    var = POLYMER_STRESS11;
    if (pd->v[pg->imtrx][var]) {
      for (mode = 0; mode < vn->modes; mode++) {
        for (p = 0; p < VIM; p++) {
          for (q = 0; q < VIM; q++) {
            var = v_s[mode][p][q];
            if (pd->v[pg->imtrx][var]) {
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                pressure_stabilization = 0.;
                for (a = 0; a < dim; a++) {
                  meqn = R_MOMENTUM1 + a;
                  if (pd->e[pg->imtrx][meqn] & T_DIFFUSION) {
                    for (r = 0; r < dim; r++) {
                      pressure_stabilization -= (delta(p, r) * delta(a, q)) * fv->snormal[a] *
                                                bf[var]->grad_phi[j][r] *
                                                pd->etm[pg->imtrx][meqn][(LOG2_DIFFUSION)];
                    }
                    if (pd->CoordinateSystem != CARTESIAN) {
                      for (r = 0; r < VIM; r++) {
                        for (s = 0; s < VIM; s++) {
                          for (t = 0; t < VIM; t++) {
                            pressure_stabilization -= fv->snormal[a] * bf[var]->phi[j] *
                                                      delta(p, s) * delta(q, t) *
                                                      (fv->grad_e[s][r][t] * delta(a, t) +
                                                       fv->grad_e[t][s][a] * delta(r, s));
                          }
                        }
                      }
                    }
                  }
                }
                pressure_stabilization *= tau_pspg;

                d_func[0][var][j] += pressure_stabilization;
              }
            }
          }
        }
      }
    }

    /*
     * J_c_G this term is only present for PSPG
     */
    var = VELOCITY_GRADIENT11;
    if (pd->v[pg->imtrx][var]) {
      for (p = 0; p < VIM; p++) {
        for (q = 0; q < VIM; q++) {
          var = v_g[p][q];
          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              pressure_stabilization = -fv->snormal[q] * bf[var]->grad_phi[j][p] *
                                       pd->etm[pg->imtrx][R_MOMENTUM1 + q][(LOG2_DIFFUSION)];

              if (pd->CoordinateSystem != CARTESIAN) {
                for (r = 0; r < VIM; r++) {
                  pressure_stabilization -= fv->snormal[q] * bf[var]->phi[j] * fv->grad_e[p][r][q] *
                                            pd->etm[pg->imtrx][R_MOMENTUM1 + a][(LOG2_DIFFUSION)];
                }
                for (a = 0; a < dim; a++) {
                  pressure_stabilization -= fv->snormal[a] * bf[var]->phi[j] * fv->grad_e[q][p][a] *
                                            pd->etm[pg->imtrx][R_MOMENTUM1 + a][(LOG2_DIFFUSION)];
                }
              }

              pressure_stabilization *= mu * tau_pspg;

              d_func[0][var][j] += pressure_stabilization;
            }
          }
        }
      }
    }

    for (b = 0; b < ei[pg->imtrx]->ielem_dim; b++) {
      var = MESH_DISPLACEMENT1 + b;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          pressure_stabilization = 0.;
          for (a = 0; a < pd->Num_Dim; a++) {
            meqn = R_MOMENTUM1 + a;
            if (pd->e[pg->imtrx][meqn]) {
              advection_a = 0.;
              if (pd->e[pg->imtrx][meqn] & T_ADVECTION) {
                if (pd->TimeIntegration != STEADY) {
                  advection_a = -rho * (1. + 2. * tt) * bf[var]->phi[j] / dt * fv->grad_v[b][a] *
                                pd->etm[pg->imtrx][meqn][(LOG2_ADVECTION)];
                }
                for (p = 0; p < dim; p++) {
                  advection_a += rho * (fv->v[p] - x_dot[p]) * fv->d_grad_v_dmesh[p][a][b][j] *
                                 pd->etm[pg->imtrx][meqn][(LOG2_ADVECTION)];
                }
              }

              diffusion = 0.;
              if (pd->e[pg->imtrx][meqn] & T_DIFFUSION) {
                diffusion -= (d_mu->X[b][j] * div_G[a] + d_div_tau_p_dmesh[a][b][j]) *
                             pd->etm[pg->imtrx][meqn][(LOG2_DIFFUSION)];
              }

              pressure = 0.;
              var1 = PRESSURE;
              if (pd->v[pg->imtrx][var1]) {
                if (pd->e[pg->imtrx][meqn] & T_DIFFUSION) {
                  pressure =
                      fv->d_grad_P_dmesh[a][b][j] * pd->etm[pg->imtrx][meqn][(LOG2_DIFFUSION)];
                }
              }

              stress = 0.;
              var1 = POLYMER_STRESS11;
              if (pd->v[pg->imtrx][var1]) {
                if (pd->e[pg->imtrx][meqn] & T_DIFFUSION) {
                  for (mode = 0; mode < vn->modes; mode++) {
                    stress -= fv->d_div_S_dmesh[mode][a][b][j] *
                              pd->etm[pg->imtrx][meqn][(LOG2_DIFFUSION)];
                  }
                }
              }

              velocity_gradient = 0.;
              var1 = VELOCITY_GRADIENT11;
              if (pd->v[pg->imtrx][var1]) {
                if (pd->e[pg->imtrx][meqn] & T_DIFFUSION) {
                  velocity_gradient -=
                      fv->d_div_G_dmesh[a][b][j] * pd->etm[pg->imtrx][meqn][(LOG2_DIFFUSION)];
                }
              }

              source_a = 0.;
              if (pd->e[pg->imtrx][meqn] & T_SOURCE) {
                source_a -= df->X[a][b][j] * pd->etm[pg->imtrx][meqn][(LOG2_SOURCE)];
              }

              pressure_stabilization +=
                  momentum_residual[a] * fv->dsnormal_dx[a][b][j] +
                  (advection_a + source_a + pressure + stress + velocity_gradient) * fv->snormal[a];
            }
          }
          d_func[0][var][j] += tau_pspg * pressure_stabilization +
                               d_tau_pspg_dX[b][j] * momentum_residual[a] * fv->snormal[a];
        }
      }
    }
  }

} /* END of routine PSPG_consistency_bc                                      */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void fapply_CA(double *func, /* Value of the Residual Function */
               double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
               double d_func_ss[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
               const double fsnormal[DIM],         /* free surface normal components      */
               double dfsnormal_dx[DIM][DIM][MDE], /* free surface normal
                                                    * derivatives ([i][j][k])
                                                    * component i wrt
                                                    * displacement j at
                                                    * node k                 */
               const double ssnormal[DIM],         /* solid surface normal components     */
               double dssnormal_dx[DIM][DIM][MDE], /* solid surface normal
                                                    * derivatives ([i][j][k])
                                                    * component i wrt
                                                    * displacement j at
                                                    * node k                 */
               const double contact_angle)         /*  Static or dynamic contact angle        */

/******************************************************************************
 *
 *  Function which applies a contact angle boundary condition in place of n.(v -vs)=0
 *
 *            Author: P. R. Schunk    (4/28/94)
 *            Revised: 6/1/95 RAC
 *            Revised: 6/7/95 RRR
 ******************************************************************************/

{
  int j, var, p, q;

  /***************************** EXECUTION BEGINS ******************************/

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  /* calculate the jacobian contributions */
  if (af->Assemble_Jacobian) {
    /*   free surface derivatives */
    for (p = 0; p < pd->Num_Dim; p++) {
      var = MESH_DISPLACEMENT1 + p;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          if (pd->Num_Dim == 2 && contact_angle > 0.75 * M_PIE) {
            d_func[0][var][j] +=
                dfsnormal_dx[0][p][j] * ssnormal[1] - dfsnormal_dx[1][p][j] * ssnormal[0];
          } else {
            for (q = 0; q < pd->Num_Dim; q++) {
              d_func[0][var][j] += dfsnormal_dx[q][p][j] * ssnormal[q];
            }
          }
          d_func[0][var][j] *= BIG_PENALTY;
        }
      }
    }
    /**   solid surface derivatives         */
    for (p = 0; p < pd->Num_Dim; p++) {
      var = MESH_DISPLACEMENT1 + p;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          if (pd->Num_Dim == 2 && contact_angle > 0.75 * M_PIE) {
            d_func_ss[0][var][j] +=
                fsnormal[0] * dssnormal_dx[1][p][j] - fsnormal[1] * dssnormal_dx[0][p][j];
          } else {
            for (q = 0; q < pd->Num_Dim; q++) {
              d_func_ss[0][var][j] += fsnormal[q] * dssnormal_dx[q][p][j];
            }
          }
          d_func_ss[0][var][j] *= BIG_PENALTY;
        }
      }
    }
  }

  /* Calculate the residual contribution	*/
  if (pd->Num_Dim != 1) {
    *func = -cos(contact_angle);
    for (p = 0; p < pd->Num_Dim; p++) {
      *func += fsnormal[p] * ssnormal[p];
    }

    /*  modify for big angles (use sine and cross-product)  */
    if (pd->Num_Dim == 2 && contact_angle > 0.75 * M_PIE) {
      *func = sin(contact_angle);
      *func += fsnormal[0] * ssnormal[1] - ssnormal[0] * fsnormal[1];
    }

    *func *= BIG_PENALTY;
  }
  /* NOTE this BIG_PENALTY results in CA getting much more weight than KINEMATIC
     because BIG_PENALTY is then multiplied twice times this bc */

} /* END of routine fapply_CA  */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void fapply_var_CA(double *func,
                   double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                   double d_func_ss[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                   const double fsnormal[MAX_PDIM], /* free surface normal components      */
                   double dfsnormal_dx[MAX_PDIM][MAX_PDIM][MDE], /* free surface normal
                                                                  * derivatives ([i][j][k])
                                                                  * component i wrt
                                                                  * displacement j at
                                                                  * node k                 */
                   const double ssnormal[MAX_PDIM], /* solid surface normal components     */
                   double dssnormal_dx[MAX_PDIM][MAX_PDIM][MDE], /* solid surface normal
                                                                  * derivatives ([i][j][k])
                                                                  * component i wrt
                                                                  * displacement j at
                                                                  * node k                 */
                   const double clnormal[MAX_PDIM], /* contact line normal components      */
                   double dclnormal_dx[MAX_PDIM][MAX_PDIM][MDE], /* contact line normal
                                                                  * derivatives ([i][j][k])
                                                                  * component i wrt
                                                                  * displacement j at
                                                                  * node k                 */
                   const dbl BC_Data_Float[8],  /* static contact angle, response slope,
                                                 * components of web speed vector         */
                   const double xdot[MAX_PDIM], /* Current mesh velocity vector */
                   const double tt,             /* parameter to vary time integration method
                                                 * from explicit (tt = 1) to
                                                 * implicit (tt = 0)                         */
                   const double dt)             /* time step size */

/*****************************************************************************
 *
 ******************************************************************************/

{
  int j, var, p, q;
  double vw[3], theta_static, Ca_local, cT, cos_CA, mu, sigma;
  double d_cos_CA_dx[MAX_PDIM][MDE];

  /***************************** EXECUTION BEGINS ******************************/

  mu = gn->mu0;
  sigma = mp->surface_tension;

  theta_static = BC_Data_Float[0] * M_PIE / 180.0;

  cT = BC_Data_Float[1];

  if (TimeIntegration == TRANSIENT) {
    vw[0] = BC_Data_Float[2] - xdot[0];
    vw[1] = BC_Data_Float[3] - xdot[1];
    vw[2] = BC_Data_Float[4] - xdot[2];
  } else {
    vw[0] = BC_Data_Float[2];
    vw[1] = BC_Data_Float[3];
    vw[2] = BC_Data_Float[4];
  }

  Ca_local = 0.0;
  memset(d_cos_CA_dx, 0, MAX_PDIM * MDE * sizeof(double));

  for (p = 0; p < pd->Num_Dim; p++) {
    Ca_local += mu * (clnormal[p] * vw[p]) / sigma;

    for (q = 0; q < pd->Num_Dim; q++) {
      var = MESH_DISPLACEMENT1 + q;
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_cos_CA_dx[q][j] -= cT * mu * dclnormal_dx[p][q][j] * vw[p] / sigma;

        if (TimeIntegration == TRANSIENT) {
          d_cos_CA_dx[q][j] +=
              cT * mu * clnormal[p] * ((1. + 2. * tt) * bf[var]->phi[j] / dt) * delta(p, q) / sigma;
        }
      }
    }
  }

  cos_CA = cos(theta_static) - cT * Ca_local;

  /* Set upper and lower limit on contact angle
   * to be 1 or 179 degrees, respectively
   */

  if (cos_CA < -1.0) {
    cos_CA = cos(179. * M_PIE / 180.);
    memset(d_cos_CA_dx, 0, MAX_PDIM * MDE);
  } else if (cos_CA > 1.0) {
    cos_CA = cos(1. * M_PIE / 180.);
    memset(d_cos_CA_dx, 0, MAX_PDIM * MDE);
  }

  if (af->Assemble_LSA_Mass_Matrix) {
    for (p = 0; p < pd->Num_Dim; p++) {
      var = MESH_DISPLACEMENT1 + p;
      if (pd->v[pg->imtrx][var])
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
          d_func[0][var][j] -= cT * mu * clnormal[p] * bf[var]->phi[j] / sigma * BIG_PENALTY;
    }
    return;
  }

  /* calculate the jacobian contributions */
  if (af->Assemble_Jacobian) {
    /*  free surface derivatives  */
    for (p = 0; p < pd->Num_Dim; p++) {
      var = MESH_DISPLACEMENT1 + p;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_func[0][var][j] = -d_cos_CA_dx[p][j];
          for (q = 0; q < pd->Num_Dim; q++) {
            d_func[0][var][j] += dfsnormal_dx[q][p][j] * ssnormal[q];
          }

          d_func[0][var][j] *= BIG_PENALTY;
        }
      }
    }
    /*  solid surface derivatives  */
    for (p = 0; p < pd->Num_Dim; p++) {
      var = MESH_DISPLACEMENT1 + p;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_func_ss[0][var][j] = 0.0;
          for (q = 0; q < pd->Num_Dim; q++) {
            d_func_ss[0][var][j] += fsnormal[q] * dssnormal_dx[q][p][j];
          }

          d_func_ss[0][var][j] *= BIG_PENALTY;
        }
      }
    }
  }

  /* Calculate the residual contribution	*/
  if (pd->Num_Dim != 1) {
    *func = -cos_CA;
    for (p = 0; p < pd->Num_Dim; p++) {
      *func += (fsnormal[p] * ssnormal[p]);
    }
    *func *= BIG_PENALTY;
  }
  /* NOTE this BIG_PENALTY results in CA getting much more weight than KINEMATIC
     because  BIG_PENALTY is then multiplied twice times this bc */

} /* END of routine fapply_var_CA  */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void fapply_var_CA_user(double *func,
                        double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                        double d_func_ss[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                        const double fsnormal[MAX_PDIM], /* free surface normal components      */
                        double dfsnormal_dx[MAX_PDIM][MAX_PDIM][MDE], /* free surface normal
                                                                       * derivatives ([i][j][k])
                                                                       * component i wrt
                                                                       * displacement j at
                                                                       * node k                 */
                        const double ssnormal[MAX_PDIM], /* solid surface normal components     */
                        double dssnormal_dx[MAX_PDIM][MAX_PDIM][MDE], /* solid surface normal
                                                                       * derivatives ([i][j][k])
                                                                       * component i wrt
                                                                       * displacement j at
                                                                       * node k                 */
                        const double clnormal[MAX_PDIM], /* contact line normal components      */
                        double dclnormal_dx[MAX_PDIM][MAX_PDIM][MDE], /* contact line normal
                                                                       * derivatives ([i][j][k])
                                                                       * component i wrt
                                                                       * displacement j at
                                                                       * node k                 */
                        const int num_user_const,    /* number of user constants in user model    */
                        const double *user_const,    /* pointer to array of user constants        */
                        const double xdot[MAX_PDIM], /* Current mesh velocity vector            */
                        const double tt,             /* vary time integration method              *
                                                      * explicit (tt = 1) to implicit (tt = 0)    */
                        const double dt)             /* time step size                            */

/************************************************************************
 *
 *
 ************************************************************************/
{
  int j, var, p, q;
  double vw[3], Ca_local, cos_CA, d_cos_CA_Ca_local;
  double d_cos_CA_dx[MAX_PDIM][MDE];

  /***************************** EXECUTION BEGINS ******************************/

  if (TimeIntegration == TRANSIENT) {
    vw[0] = user_const[0] - xdot[0];
    vw[1] = user_const[1] - xdot[1];
    vw[2] = user_const[2] - xdot[2];
  } else {
    vw[0] = user_const[0];
    vw[1] = user_const[1];
    vw[2] = user_const[2];
  }

  Ca_local = 0.0;
  memset(d_cos_CA_dx, 0, MAX_PDIM * MDE * sizeof(double));

  for (p = 0; p < pd->Num_Dim; p++) {
    Ca_local += mp->viscosity * (clnormal[p] * vw[p]) / mp->surface_tension;
  }

  cos_CA = var_CA_user(Ca_local, num_user_const, &user_const[3], &d_cos_CA_Ca_local);

  for (p = 0; p < pd->Num_Dim; p++) {

    for (q = 0; q < pd->Num_Dim; q++) {
      var = MESH_DISPLACEMENT1 + q;
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_cos_CA_dx[q][j] +=
            d_cos_CA_Ca_local * mp->viscosity * dclnormal_dx[p][q][j] * vw[p] / mp->surface_tension;

        if (TimeIntegration == TRANSIENT) {
          d_cos_CA_dx[q][j] -= d_cos_CA_Ca_local * mp->viscosity * clnormal[p] *
                               ((1. + 2. * tt) * bf[var]->phi[j] / dt) * delta(p, q) /
                               mp->surface_tension;
        }
      }
    }
  }

  /* calculate the jacobian contributions */
  if (af->Assemble_Jacobian) {
    for (p = 0; p < pd->Num_Dim; p++) {
      var = MESH_DISPLACEMENT1 + p;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_func[0][var][j] = d_cos_CA_dx[p][j];
          for (q = 0; q < pd->Num_Dim; q++) {
            d_func[0][var][j] += dfsnormal_dx[q][p][j] * ssnormal[q];
          }

          d_func[0][var][j] *= BIG_PENALTY;
        }
      }
    }
    /*   solid surface contribution  */
    for (p = 0; p < pd->Num_Dim; p++) {
      var = MESH_DISPLACEMENT1 + p;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_func_ss[0][var][j] = 0.0;
          for (q = 0; q < pd->Num_Dim; q++) {
            d_func_ss[0][var][j] += fsnormal[q] * dssnormal_dx[q][p][j];
          }

          d_func_ss[0][var][j] *= BIG_PENALTY;
        }
      }
    }
  }

  /* Calculate the residual contribution	*/
  /* NOTE : the residual equation is (fsnormal)dot(-ssnormal) = cos_CA
   * so that our contact is in agreement with the rest of the world
   */

  if (pd->Num_Dim != 1) {
    *func = cos_CA;
    for (p = 0; p < pd->Num_Dim; p++) {
      *func += (fsnormal[p] * ssnormal[p]);
    }
    *func *= BIG_PENALTY;
  }
  /* NOTE this BIG_PENALTY results in CA getting much more weight than KINEMATIC
     because  BIG_PENALTY is then multiplied twice times this bc */

} /* END of routine fapply_var_CA_user  */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int evaluate_gibbs_criterion(const double fsnormal[MAX_PDIM], /* free surface normal components */
                             const double ssnormal[MAX_PDIM], /* solid surface normal components */
                             int *ipin, /* Flag which tracks whether pinned or not   */
                             const double contact_angle, /* Static or dynamic contact angle */
                             const double x_pos, /* x-coordinate of sharp edge                */
                             const double y_pos, /* y-coordinate of sharp edge                */
                             const double z_pos, /* z-coordinate of sharp edge                */
                             const double sign_origx, /* original relative sign on x-position of   *
                                                       * contact line with specified sharp edge */
                             const double sign_origy, /* original relative sign on y-position of   *
                                                       * contact line with specified sharp edge */
                             const double sign_origz) /* original relative sign on z-position of   *
                                                       * contact line with specified sharp edge */

/*******************************************************************************
 *
 *  Function which evaluates the Gibb's inequality criterion for pinning
 *  and releasing a contact line, viz.,
 *    (theta - theta_s)*(dist_from_salient_point) = 0 with
 *        theta_s-theta >= 0, hs >= 0
 *
 *    If contact line should be or stay released (return 1)
 *    If contact line should be fixed (return 0)
 *
 *            Author: P. R. Schunk    (12/16/97)
 *            Author: P. R. Schunk    (12/16/97)
 *
 * The fact the Randy has "Author"'ed this twice opens the door
 * for a slew of jokes...
 ******************************************************************************/

{
  int p;
  dbl actual_angle, dot_prod, pos[3];

  /***************************** EXECUTION BEGINS ******************************/

  /* Compute distance from sharpe edge */

  pos[0] = x_pos - fv->x[0];
  pos[1] = y_pos - fv->x[1];
  if (pd->Num_Dim == 3) {
    pos[2] = x_pos - fv->x[2];
  }

  /* 2D only for now */

  if (pd->Num_Dim == 3)
    GOMA_EH(GOMA_ERROR, "CA_OR_FIX only in 2D now");

  if (((sign_of(pos[0]) == sign_of(sign_origx) && fabs(pos[0]) > 1.e-6)) ||
      ((sign_of(pos[1]) == sign_of(sign_origy)) && fabs(pos[1]) > 1.e-6)) {
    *ipin = 0;
    return (1);

  }

  else {
    if (!*ipin) {
      *ipin = 1;
      return (0);
    } else {
      *ipin = 1;
    }

    /* if dist is basically zero, or the line has gone past the feature,
       the contact line should be fixed from the previous iteration or should
       be fixed now.  Evaluate Gibbs Criterion */

    dot_prod = 0.;
    for (p = 0; p < pd->Num_Dim; p++) {
      dot_prod += fsnormal[p] * ssnormal[p];
    }
    actual_angle = acos(dot_prod);
    /* evaluate gibbs criterion here */
    if (actual_angle >= (contact_angle + 1.e-1)) {
      *ipin = 0;
      return (1);
    } else {
      *ipin = 1;
      return (0);
    }
  }

} /* END of routine evaluate_gibbs_criterion  */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void fapply_moving_CA(double *func,
                      double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                      double d_func_ss[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                      const double fsnormal[MAX_PDIM], /* free surface normal components          */
                      double dfsnormal_dx[MAX_PDIM][MAX_PDIM][MDE], /* free surface normal        *
                                                                     * derivatives ([i][j][k])    *
                                                                     * component i wrt            *
                                                                     * displacement j at node k   */
                      const double ssnormal[MAX_PDIM], /* solid surface normal components         */
                      double dssnormal_dx[MAX_PDIM][MAX_PDIM][MDE], /* solid surface normal       *
                                                                     * derivatives ([i][j][k])    *
                                                                     * component i wrt            *
                                                                     * displacement j at node k   */
                      const double stat_ca,         /* Static or dynamic contact angle           */
                      const double advancing_ca,    /* Static or dynamic contact angle           */
                      const double receding_ca,     /* Static or dynamic contact angle           */
                      const double scaling,         /* Static or dynamic contact angle           */
                      const double vwx,             /* x-Velocity of wall                        */
                      const double vwy,             /* y-Velocity of wall                        */
                      const double vwz,             /* z-Velocity of wall                        */
                      const double x_dot[MAX_PDIM], /* mesh velocity vector                      */
                      const dbl dt,                 /* current value of the time step            */
                      const dbl tt)                 /* parameter to vary time integration:       *
                                                     * explicit (tt = 1) -- implicit (tt = 0)    */

/******************************************************************************
 *
 *  Function which applies a contact angle boundary condition in place of n.(v -vs)=0
 *
 *            Author: P. R. Schunk    (4/28/94)
 *            Revised: 6/1/95 RAC
 *            Revised: 6/7/95 RRR
 *            Revised: 11/19/95 RAC
 ******************************************************************************/

{
  int j, var, p, q;
  double v_rel, ca;
  double vwall[MAX_PDIM];

  /***************************** EXECUTION BEGINS ******************************/
  /* find contact angle as a function of the contact line speed relative to the
   * wall speed - currently assume wall is stationary*/

  /* note that the angles in this routine are in degrees and converted to radians
   * below */

  vwall[0] = vwx;
  vwall[1] = vwy;
  vwall[2] = vwz;

  /* calculate the speed of the free surface relative to the solid surface
   * currently assume that the solid surface is not moving, and use a
   * trick of dotting the free surface velocity into the free surface normal
   * to find out if the contact line is advancing (+) or receding (-) */
  v_rel = 0.;
  for (p = 0; p < Num_Dim; p++) {
    v_rel += (x_dot[p] - vwall[p]) * fsnormal[p];
  }
  ca = (stat_ca + (advancing_ca - stat_ca) * tanh(scaling * v_rel)) * M_PIE / 180.;

  if (af->Assemble_LSA_Mass_Matrix) {
    for (p = 0; p < pd->Num_Dim; p++) {
      var = MESH_DISPLACEMENT1 + p;
      if (pd->v[pg->imtrx][var])
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
          d_func[0][var][j] +=
              BIG_PENALTY *
              (sin(ca) * (advancing_ca - stat_ca) * scaling / cosh(scaling * v_rel) /
               cosh(scaling * v_rel) * fsnormal[p] * bf[var]->phi[j]) *
              M_PIE / 180.0;
    }
    return;
  }

  /* calculate the jacobian contributions */
  if (af->Assemble_Jacobian) {
    /*   free surface contributions  */
    for (p = 0; p < pd->Num_Dim; p++) {
      var = MESH_DISPLACEMENT1 + p;
      if (pd->v[pg->imtrx][var]) {
        for (q = 0; q < pd->Num_Dim; q++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_func[0][var][j] +=
                BIG_PENALTY * (dfsnormal_dx[q][p][j] * ssnormal[q] +
                               (sin(ca) * (advancing_ca - stat_ca) * scaling /
                                cosh(scaling * v_rel) / cosh(scaling * v_rel) *
                                dfsnormal_dx[q][p][j] * (x_dot[q] - vwall[q]) * M_PIE / 180.));
          }
        }
        if (TimeIntegration != 0) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_func[0][var][j] +=
                BIG_PENALTY *
                (sin(ca) * (advancing_ca - stat_ca) * scaling / cosh(scaling * v_rel) /
                 cosh(scaling * v_rel) * fsnormal[p] * (1. + 2. * tt) * bf[var]->phi[j] / dt) *
                M_PIE / 180.;
          }
        }
      }
    }
    /*   solid surface contributions  */
    for (p = 0; p < pd->Num_Dim; p++) {
      var = MESH_DISPLACEMENT1 + p;
      if (pd->v[pg->imtrx][var]) {
        for (q = 0; q < pd->Num_Dim; q++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_func_ss[0][var][j] += BIG_PENALTY * fsnormal[q] * dssnormal_dx[q][p][j];
          }
        }
      }
    }
  }

  /* Calculate the residual contribution	*/
  if (pd->Num_Dim != 1) {
    *func = -cos(ca) * BIG_PENALTY;
    for (p = 0; p < pd->Num_Dim; p++) {
      *func += BIG_PENALTY * fsnormal[p] * ssnormal[p];
    }
  }
  /* NOTE this BIG_PENALTY results in CA getting much more weight than KINEMATIC
     because  BIG_PENALTY is then multiplied twice times this bc */

} /* END of routine fapply_moving_CA  */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/*
 *			      [0]    [1]  [2]  [3]  [4] [5] [6]  [7]
 *			int   flt    flt  flt  flt  flt flt flt  flt
 *  BC = VELO_THETA_TPL nsid theta_0 nssx nssy nssz v0 g t_relax v_old
 *
 * Constitutive relation for apparent contact line speed as a function of
 * apparent instantaneous contact angle.
 *
 * Blake model, but with relaxation time.
 *
 *        v = v_old + ( v_new - v_old ) * ( 1 - exp(-t/t_relax) )
 *
 * where
 *
 *        v_new = V_0 * sinh ( cos(theta_eq) - cos(theta ) * g )
 *
 * where
 *        V_0      is a pre-exponential velocity,
 *        theta_eq is an equilibrium contact angle, where v=0
 *        g        is a thermally-scaled surface energy (tension)
 *                 = gamma/(2nkT) where gamma is the surface
 *                   energy, n is number sites per unit area, kT
 *                   is thermal energy.
 *
 *
 * The actual residual equation is:
 *
 *		R = 0 = v_blake - t.dxdt
 *
 * where we've implicitly assumed the other constraints at the trijunction
 * will insure that v_fluid = dxdt.
 *
 * It may be better to set t.v_fluid = v_blake because that would enable
 * us to solve problems where dxdt=0 in the lab reference frame
 * (i.e., steady coating flows) as well as transient spreading.
 *
 *
 * Adapted from fapply_moving_CA() see above.
 *
 * Revised: 2004/04/14 pasacki@sandia.gov
 */

void fapply_moving_CA_sinh(
    double *func,
    double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
    double d_func_ss[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
    const double fsnormal[MAX_PDIM],              /* free surface normal components          */
    double dfsnormal_dx[MAX_PDIM][MAX_PDIM][MDE], /* free surface normal  *
                                                   * derivatives ([i][j][k])    *
                                                   * component i wrt            *
                                                   * displacement j at node k   */
    const double ssnormal[MAX_PDIM],              /* solid surface normal components         */
    double dssnormal_dx[MAX_PDIM][MAX_PDIM][MDE], /* solid surface normal *
                                                   * derivatives ([i][j][k])    *
                                                   * component i wrt            *
                                                   * displacement j at node k   */
    const double equilibrium_contact_angle,
    const double velocity_pre_exponential,
    const double modified_surface_tension,
    const double relaxation_time,
    const double old_tpl_velocity,
    const double x_dot[MAX_PDIM], /* mesh velocity vector                      */
    const dbl dt,                 /* current value of the time step            */
    const dbl tt,                 /* explicit (tt = 1) -- implicit (tt = 0)    */
    const double time,            /* current time */
    const double wall_velocity,
    const double theta_max_degrees,
    const double dewet_input,
    const double dcl_shearrate,
    const int bc_type, /*  bc identifier	*/
    double dwall_velo_dx[MAX_PDIM][MDE],
    const int local_node) {
  int j, var, p, q, w;

  double sstangent[MAX_PDIM];
  double dsstangent_qpj = 0.0;

  double sign;
  double dnnddpj, dnn_ss_ddpj;
  double dvmesh_ddpj, dvmesh_ss_ddpj;
  double dv_ddpj = 0.0, dv_ss_ddpj = 0.0;
  double costheta;

  const double t_relax = relaxation_time;
  double g; /* scaled surface tension	*/
  const double costhetaeq = cos(equilibrium_contact_angle * (M_PIE / 180));
  const double v0 = velocity_pre_exponential;
  const double t = time;
  double factor;
  double v;
  double v_new = 0.0;
  const double v_old = old_tpl_velocity;
  double v_mesh, v_mesh_dt;

  /*  Hoffman correlation variables	*/
  double ca_no = 0.0, g_sca = 0.0, g_dca = 0.0;
#ifdef NEW_HOFFMAN_FCN_PLEASE
  double g_deriv, g_deriv_ss;
  double hoff_C = 0.012874005, hoff_N = 2.80906762, hoff_F = 0.7093681;
  double hoff_M = 1.253351327, hoff_R = 9.614608063;
  double hoff_D = velocity_pre_exponential * M_PIE / 180.0;
#else
  double A_sca = 0.0, A_dca = 0.0;
  int iter, iter_max = 20;
  double eps_tol = 1.0e-12;
#endif
  double liq_visc = 0.0, gamma[DIM][DIM];
  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct; /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;

  /*  Cox wetting parameters	*/
  double theta = 0.0, thetaeq = 0.0, sintheta = 0.0, sinthetaeq = 0.0;
  double lambda = 0;
  double th, sinth, costh;
  double f_sca = 0.0, f_dca = 0.0, f_num = 0.0, f_den = 0.0, g_integral = 0.0;
  double reciprocal_slip = 0.0;
  const double q_inner = 0;
  const double q_outer = 0;
  /* use 10 point Gaussian quadrature for now	*/
  const int num_gauss_pts = 10;
  const double gpt[10] = {-0.973906528517172, -0.865063366688985, -0.679409568299024,
                          -0.433395394129247, -0.148874338981631, 0.148874338981631,
                          0.433395394129247,  0.679409568299024,  0.865063366688985,
                          0.973906528517172};
  const double wt[10] = {0.066671344308688, 0.149451349150581, 0.219086362515982, 0.269266719309996,
                         0.295524224714753, 0.295524224714753, 0.269266719309996, 0.219086362515982,
                         0.149451349150581, 0.066671344308688};
  static const char yo[] = "fapply_moving_CA_sinh";
  /* Shikhmurzaev wetting parameters	*/
  double rhs = 0.0, rhs_den = 0.0, rhs_num = 0.0, veloc0 = 0.0, veloc0max;
  double drhs_ddpj, drhs_den_ddpj, drhs_num_ddpj, dveloc0_ddpj;
  double drhs_ss_ddpj, drhs_den_ss_ddpj, drhs_num_ss_ddpj, dveloc0_ss_ddpj;
  double theta_max = 0.0, costhetamax = 0.0, sinthetamax = 0.0, dewet = 0.0;
  const double shik_max_factor = 1.01;
  /*  const double wall_sign = (TimeIntegration == STEADY) ? 1 : -1;  */
  /* disabling this sign change for now - doesn't seem necessary*/
  const double wall_sign = (TimeIntegration == STEADY) ? 1 : 1;

  /*
   * What's the cosine of the current instantaneous contact angle, the
   * angle between the free surface normal and solid surface normal?
   */
  costheta = (fsnormal[0] * ssnormal[0] + fsnormal[1] * ssnormal[1] + fsnormal[2] * ssnormal[2]);
  sintheta = -(fsnormal[0] * ssnormal[1] - ssnormal[0] * fsnormal[1]);
  /*  if(sintheta < 0 ) costheta = - 0.99; */

  switch (bc_type) {
  case VELO_THETA_COX_BC:
    reciprocal_slip = 1. / velocity_pre_exponential;
  /* fall through */
  case VELO_THETA_SHIK_BC:
    theta_max = M_PIE * theta_max_degrees / 180.;
    thetaeq = equilibrium_contact_angle * (M_PIE / 180);
    theta = acos(costheta);
    sintheta = sin(theta);
    sinthetaeq = sin(thetaeq);
    if (theta > theta_max) {
      theta_max = MIN(shik_max_factor * theta, M_PIE);
    }
    costhetamax = cos(theta_max);
    sinthetamax = sin(theta_max);
    break;
  case VELO_THETA_HOFFMAN_BC:
    theta_max = M_PIE * theta_max_degrees / 180.;
    costhetamax = cos(theta_max);
#ifdef NEW_HOFFMAN_FCN_PLEASE
    hoff_R = pow(hoff_F, hoff_N) * pow(1. - hoff_F, hoff_M) * pow(theta_max, hoff_N + hoff_M);
    theta = acos(costheta);
    sintheta = sin(theta);
    thetaeq = equilibrium_contact_angle * (M_PIE / 180);
#else
    costheta = MAX(costheta, costhetamax);
#endif
    break;
  }

  if (bc_type == VELO_THETA_COX_BC || bc_type == VELO_THETA_HOFFMAN_BC) {
    if (dcl_shearrate > 0) {
      memset(gamma, 0, sizeof(double) * DIM * DIM);
      gamma[0][1] = gamma[1][0] = dcl_shearrate;
      fv->T = mp->reference[TEMPERATURE];
      liq_visc = viscosity(gn, gamma, NULL);
    } else {
      for (p = 0; p < VIM; p++) {
        for (q = 0; q < VIM; q++) {
          gamma[p][q] = fv->grad_v[p][q] + fv->grad_v[q][p];
        }
      }
      liq_visc = viscosity(gn, gamma, d_mu);
    }
  }

  /*
   * Based on these input properties and the current contact angle, we
   * can compute what the velocity ought to be...
   */

  g = modified_surface_tension * mp->surface_tension;
  switch (bc_type) {
  case VELO_THETA_TPL_BC:
    v_new = v0 * sinh((costhetaeq - costheta) * g);
    dewet = (v_new < 0) ? dewet_input : 1.0;
    v_new *= dewet;
    break;
  case VELO_THETA_HOFFMAN_BC:
#ifdef NEW_HOFFMAN_FCN_PLEASE
    if (thetaeq < hoff_F * theta_max) {
      g_sca = hoff_C * pow(thetaeq, hoff_N);
    } else {
      g_sca = hoff_C * hoff_R / pow(theta_max - thetaeq, hoff_M);
    }
    if (theta < hoff_F * theta_max) {
      g_dca = hoff_C * pow(theta, hoff_N);
    } else if (theta < theta_max - hoff_D) {
      g_dca = hoff_C * hoff_R / pow(theta_max - theta, hoff_M);
    } else {
      g_dca = hoff_C * hoff_R / pow(hoff_D, hoff_M) *
              (1.0 + hoff_M / hoff_D * (theta - theta_max + hoff_D));
    }
    if (!isfinite(g_dca)) {
      g_dca = SGN(g_dca) * BIG_PENALTY;
    }
#else
    ca_no = 1.0E+06;
    iter = 0;
    eps = 10. * eps_tol;
    while (iter <= iter_max && fabs(eps) > eps_tol) {
      g_sca = log((3. - costhetaeq) / (1. + costhetaeq)) / (2. * 5.16);
      A_sca = pow(g_sca, 1. / 0.706);
      eps = -(ca_no - 1.31 * pow(ca_no, 0.99) * A_sca - A_sca) /
            (1. - 1.31 * 0.99 * A_sca / pow(ca_no, 0.01));
      ca_no += eps;
      iter++;
    }
    if (fabs(eps) > eps_tol)
      GOMA_EH(GOMA_ERROR, "Hoffman iteration not converged");
    g_sca = ca_no;
    ca_no = 1.0E+06;
    iter = 0;
    eps = 10. * eps_tol;
    while (iter <= iter_max && fabs(eps) > eps_tol) {
      g_dca = log((3. - costheta) / (1. + costheta)) / (2. * 5.16);
      A_dca = pow(g_dca, 1. / 0.706);
      eps = -(ca_no - 1.31 * pow(ca_no, 0.99) * A_dca - A_dca) /
            (1. - 1.31 * 0.99 * A_dca / pow(ca_no, 0.01));
      ca_no += eps;
      if (!isfinite(ca_no)) {
        ca_no = eps = DBL_MAX / 10.;
      }
      iter++;
    }
    if (fabs(eps) > eps_tol)
      fprintf(stderr, "Hoffman not converged ... %d %g\n", iter, eps);
    g_dca = ca_no;
#endif
    ca_no = g_dca - g_sca;
    dewet = (ca_no < 0) ? dewet_input : 1.0;
    v_new = dewet * ca_no * g / liq_visc;
    break;
  case VELO_THETA_COX_BC:
    /*  Cox analysis integral	*/
    g_integral = 0;
    for (j = 0; j < num_gauss_pts; j++) {
      th = thetaeq + (theta - thetaeq) * (gpt[j] + 1.) / 2.;
      sinth = sin(th);
      costh = cos(th);
      f_num = 2. * sinth *
              (SQUARE(lambda) * (SQUARE(th) - SQUARE(sinth)) +
               2. * lambda * (th * (M_PIE - th) + SQUARE(sinth)) +
               (SQUARE(M_PIE - th) - SQUARE(sinth)));
      f_den = lambda * (SQUARE(th) - SQUARE(sinth)) * (M_PIE - th + sinth * costh) +
              (SQUARE(M_PIE - th) - SQUARE(sinth)) * (th - sinth * costh);
      g_integral += wt[j] * f_den / f_num;
    }
    g_integral *= 0.5 * (theta - thetaeq);
    /*  inner and outer solution terms	*/
    f_num = 2. * sinthetaeq *
            (SQUARE(lambda) * (SQUARE(thetaeq) - SQUARE(sinthetaeq)) +
             2. * lambda * (thetaeq * (M_PIE - thetaeq) + SQUARE(sinthetaeq)) +
             (SQUARE(M_PIE - thetaeq) - SQUARE(sinthetaeq)));
    f_den = lambda * (SQUARE(thetaeq) - SQUARE(sinthetaeq)) *
                (M_PIE - thetaeq + sinthetaeq * costhetaeq) +
            (SQUARE(M_PIE - thetaeq) - SQUARE(sinthetaeq)) * (thetaeq - sinthetaeq * costhetaeq);
    f_sca = f_num / f_den;
    f_num = 2. * sintheta *
            (SQUARE(lambda) * (SQUARE(theta) - SQUARE(sintheta)) +
             2. * lambda * (theta * (M_PIE - theta) + SQUARE(sintheta)) +
             (SQUARE(M_PIE - theta) - SQUARE(sintheta)));
    f_den = lambda * (SQUARE(theta) - SQUARE(sintheta)) * (M_PIE - theta + sintheta * costheta) +
            (SQUARE(M_PIE - theta) - SQUARE(sintheta)) * (theta - sintheta * costheta);
    f_dca = f_num / f_den;
    /* solve for wetting speed	*/
    ca_no = g_integral / (log(reciprocal_slip) + (q_inner / f_sca - q_outer / f_dca));
    dewet = (ca_no < 0) ? dewet_input : 1.0;
    v_new = dewet * ca_no * g / liq_visc;
    break;
  case VELO_THETA_SHIK_BC:
    veloc0 = (sintheta - theta * costheta) / (sintheta * costheta - theta);
    veloc0max = (sinthetamax - theta_max * costhetamax) / (sinthetamax * costhetamax - theta_max);
    rhs_num = costhetaeq - costheta;
    rhs_den = (v0 - 1.) * (veloc0 - veloc0max) + costheta - costhetamax;
    rhs = rhs_num / rhs_den;
    v_new = sqrt(g * v0) * rhs / (2. * sqrt(1. + rhs));
    break;
  default:
    GOMA_EH(GOMA_ERROR, "bad DCA bc name\n");
  }

  /*
   * If time dependent, possibly relax any abrupt changes...
   */

  if (t_relax <= 0.) {
    factor = 1;
  } else if (t / t_relax > 36) {
    factor = 1;
  } else {
    factor = 1 - exp(-t / t_relax);
  }

  if (TimeIntegration != 0) {
    v = v_old + (v_new - v_old) * factor;
  } else {
    v = v_new;
  }

  /*
   * Velocity is considered positive when the free surface moves tangentially
   * across the solid surface such that v.n_LV > 0. Otherwise, it's negative.
   *
   * Get a solid surface tangent vector (sorry 2D only for now)...
   */

  /*
   * Which way is the outward normal of the free surface relative to the
   * solid surface tangent vector?
   *
   * Note that the derivative of sgn(x) is delta(x)...
   */

  if (pd->Num_Dim > 2) {
    log_err("BC VELO_THETA_TPL relies on 2D for ez tangents.");
  }

  sstangent[0] = -ssnormal[1];
  sstangent[1] = ssnormal[0];
  sstangent[2] = 0;

  sign = (sstangent[0] * fsnormal[0] + sstangent[1] * fsnormal[1] + sstangent[2] * fsnormal[2]) > 0
             ? 1
             : -1;

  /*sign = 1.0;*/

  v_mesh = sign * (wall_sign * wall_velocity + sstangent[0] * x_dot[0] + sstangent[1] * x_dot[1]);
  /*	+ sstangent[2] * x_dot[2] );	*/
  v_mesh_dt = sign * (sstangent[0] * fv_dot->v[0] + sstangent[1] * fv_dot->v[1]);

  /*
   * Residual equation
   *
   *          v - t.dxdt = 0
   */

#ifdef ALE_DCA_INFO_PLEASE
  if (af->Assemble_Jacobian)
    fprintf(stderr, "         ALE_DCA INFO - v_wetting: %g v_mesh: %g dvmdt: %g DCA: %g\n", v,
            v_mesh, v_mesh_dt, acos(costheta) * 180 / M_PIE);
#endif
#if 0
fprintf(stderr,"fs_normal: %g %g ss_normal: %g %g\n",fsnormal[0],fsnormal[1],ssnormal[0],ssnormal[1]);
fprintf(stderr,"\nwall_v: %g  x_dot: %g %g \n",wall_sign*wall_velocity,x_dot[0],x_dot[1]);
fprintf(stderr,"cos: %g  sin: %g  CA#:  %g \n",costheta, sintheta, ca_no);
fprintf(stderr,"dewet %g sign  %g wall_sign %g\n",dewet, sign, wall_sign);
fprintf(stderr,"velocity  %g v_mesh %g dvmdt: %g node: %d\n",v,v_mesh,v_mesh_dt,local_node);
#endif

  if (pd->Num_Dim != 1) {
    *func = (v - v_mesh);
    if (TimeIntegration == TRANSIENT) {
      *func += t_relax * v_mesh_dt;
    }
    *func *= BIG_PENALTY;
  }

  /* NOTE this BIG_PENALTY results in CA getting much more weight
   * than KINEMATIC because BIG_PENALTY is then multiplied twice times
   * this bc
   */

  if (!af->Assemble_Jacobian)
    return;

  /*
   * Compute d(Residual)/d_d[p][j]
   *
   * Residual depends only on mesh geometry... not any more....
   */

  for (p = 0; p < pd->Num_Dim; p++) {
    var = MESH_DISPLACEMENT1 + p;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        dnnddpj = 0;
        dnn_ss_ddpj = 0;
        dvmesh_ddpj = 0;
        dvmesh_ss_ddpj = 0;
        for (q = 0; q < pd->Num_Dim; q++) {
          dnnddpj += dfsnormal_dx[q][p][j] * ssnormal[q];
          dnn_ss_ddpj += dssnormal_dx[q][p][j] * fsnormal[q];

          /*
           * d(sstangent[q])/d_d[p][j] = +/- d(ssnormal[D-q])/d_d[p][j]
           *
           * Could be extended to 3D with ROT seed vectors, etc.
           */

          switch (q) {
          case 0:
            dsstangent_qpj = -dssnormal_dx[1][p][j];
            break;
          case 1:
            dsstangent_qpj = dssnormal_dx[0][p][j];
            break;
          case 2:
            /*		      dsstangent_qpj = dssnormal_dx[2][p][j];		      */
            dsstangent_qpj = 0.0;
            break;
          }

          if (TimeIntegration == TRANSIENT) {
            if (j == local_node) {
              dvmesh_ddpj +=
                  sign * sstangent[q] * (1. + 2. * tt) * delta(p, q) * (bf[var]->phi[j]) / dt;
            }
            dvmesh_ss_ddpj += sign * dsstangent_qpj * x_dot[q];
            dvmesh_ss_ddpj -= t_relax * sign * dsstangent_qpj * fv_dot->v[q];
          } else {
            dvmesh_ddpj += sign * (dsstangent_qpj * wall_velocity);
          }
        }
        dvmesh_ddpj += sign * wall_sign * dwall_velo_dx[p][j];
        switch (bc_type) {
        case VELO_THETA_TPL_BC:
          dv_ddpj = dewet * v0 * cosh((costhetaeq - costheta) * g) * g * (-dnnddpj) * factor;
          dv_ss_ddpj = dewet * v0 * cosh((costhetaeq - costheta) * g) * g * (-dnn_ss_ddpj) * factor;
          break;
        case VELO_THETA_HOFFMAN_BC:
#ifdef NEW_HOFFMAN_FCN_PLEASE
          dv_ddpj = factor * dewet * g / liq_visc;
          dv_ss_ddpj = dv_ddpj;
          if (theta < hoff_F * theta_max) {
            g_deriv = hoff_C * hoff_N * pow(theta, hoff_N - 1.) / (-sintheta) * (dnnddpj);
            g_deriv_ss = hoff_C * hoff_N * pow(theta, hoff_N - 1.) / (-sintheta) * (dnn_ss_ddpj);
          } else if (theta < theta_max - hoff_D) {
            g_deriv = hoff_C * hoff_R * hoff_M / pow(theta_max - theta, hoff_M + 1.) / (-sintheta) *
                      (dnnddpj);
            g_deriv_ss = hoff_C * hoff_R * hoff_M / pow(theta_max - theta, hoff_M + 1.) /
                         (-sintheta) * (dnn_ss_ddpj);
          } else {
            g_deriv = hoff_C * hoff_R * hoff_M / pow(hoff_D, hoff_M + 1.0) / (-sintheta) * dnnddpj;
            g_deriv_ss =
                hoff_C * hoff_R * hoff_M / pow(hoff_D, hoff_M + 1.0) / (-sintheta) * dnn_ss_ddpj;
          }
          if (!isfinite(g_deriv)) {
            g_deriv = SGN(g_deriv) * BIG_PENALTY;
          }
          if (!isfinite(g_deriv_ss)) {
            g_deriv_ss = SGN(g_deriv_ss) * BIG_PENALTY;
          }
          dv_ddpj *= g_deriv;
          dv_ss_ddpj *= g_deriv_ss;
#else
          dv_ddpj = factor * dewet * g / liq_visc * (1. + 1.31 * pow(g_dca, 0.99)) /
                    (1. - 1.31 * 0.99 * A_dca / pow(g_dca, 0.01)) * pow(A_dca, 0.294) * (-4.) *
                    dnnddpj / (0.706 * 2 * 5.16 * (3. - costheta) * (1. + costheta));
          dv_ss_ddpj = factor * dewet * g / liq_visc * (1. + 1.31 * pow(g_dca, 0.99)) /
                       (1. - 1.31 * 0.99 * A_dca / pow(g_dca, 0.01)) * pow(A_dca, 0.294) * (-4.) *
                       dnn_ss_ddpj / (0.706 * 2 * 5.16 * (3. - costheta) * (1. + costheta));
#endif
          dv_ddpj += factor * dewet * (g_dca - g_sca) * g * (-d_mu->X[p][j] / SQUARE(liq_visc));
          dv_ss_ddpj += factor * dewet * (g_dca - g_sca) * g * (-d_mu->X[p][j] / SQUARE(liq_visc));
          break;
        case VELO_THETA_COX_BC:
          /* sensitivity wrt to q_outer has not been included	*/
          f_den = log(reciprocal_slip) + (q_inner / f_sca - q_outer / f_dca);
          dv_ddpj = factor * dewet * g / liq_visc * (1. / (f_dca * (-sintheta)) * dnnddpj) / f_den;
          dv_ss_ddpj =
              factor * dewet * g / liq_visc * (1. / (f_dca * (-sintheta)) * dnn_ss_ddpj) / f_den;
          break;
        case VELO_THETA_SHIK_BC:
          dveloc0_ddpj =
              dnnddpj * (SQUARE(theta) + theta * sintheta * costheta - 2 * SQUARE(sintheta)) /
              (SQUARE(sintheta * costheta) - 2 * theta * sintheta * costheta + SQUARE(theta));
          dveloc0_ss_ddpj =
              dnn_ss_ddpj * (SQUARE(theta) + theta * sintheta * costheta - 2 * SQUARE(sintheta)) /
              (SQUARE(sintheta * costheta) - 2 * theta * sintheta * costheta + SQUARE(theta));
          drhs_num_ddpj = -dnnddpj;
          drhs_num_ss_ddpj = -dnn_ss_ddpj;
          drhs_den_ddpj = (v0 - 1) * dveloc0_ddpj + dnnddpj;
          drhs_den_ss_ddpj = (v0 - 1) * dveloc0_ss_ddpj + dnn_ss_ddpj;
          if (theta_max > theta_max_degrees * M_PIE / 180.) {
            drhs_den_ddpj *= (1. - 1. / shik_max_factor);
            drhs_den_ss_ddpj *= (1. - 1. / shik_max_factor);
          }
          drhs_ddpj = (rhs_den * drhs_num_ddpj - rhs_num * drhs_den_ddpj) / SQUARE(rhs_den);
          drhs_ss_ddpj =
              (rhs_den * drhs_num_ss_ddpj - rhs_num * drhs_den_ss_ddpj) / SQUARE(rhs_den);
          dv_ddpj = factor * 0.5 * sqrt(g * v0) * (1 - 0.5 * rhs) * drhs_ddpj / pow(1 + rhs, 1.5);
          dv_ss_ddpj =
              factor * 0.5 * sqrt(g * v0) * (1 - 0.5 * rhs) * drhs_ss_ddpj / pow(1 + rhs, 1.5);
          break;
        }

        d_func[0][var][j] += BIG_PENALTY * (dv_ddpj - dvmesh_ddpj);
        d_func_ss[0][var][j] += BIG_PENALTY * (dv_ss_ddpj - dvmesh_ss_ddpj);
      }
    }
  }

  /***	velocity sensitivities thru surface tension or viscosity	**/

  for (p = 0; p < pd->Num_Dim; p++) {
    var = VELOCITY1 + p;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        switch (bc_type) {
        case VELO_THETA_TPL_BC:
          dv_ddpj = 0;
          break;
        case VELO_THETA_HOFFMAN_BC:
          dv_ddpj =
              factor * dewet * (g_dca - g_sca) *
              (liq_visc * modified_surface_tension * mp->d_surface_tension[var] * bf[var]->phi[j] -
               g * d_mu->v[p][j]) /
              SQUARE(liq_visc);
          break;
        case VELO_THETA_COX_BC:
          f_den = log(reciprocal_slip) + (q_inner / f_sca - q_outer / f_dca);
          dv_ddpj =
              factor * dewet * g_integral / f_den *
              (liq_visc * modified_surface_tension * mp->d_surface_tension[var] * bf[var]->phi[j] -
               g * d_mu->v[p][j]) /
              SQUARE(liq_visc);
          break;
        case VELO_THETA_SHIK_BC:
          dv_ddpj = factor * 0.5 * v_new / g * modified_surface_tension *
                    mp->d_surface_tension[var] * bf[var]->phi[j];
          break;
        }
        d_func[0][var][j] += BIG_PENALTY * dv_ddpj;
        d_func_ss[0][var][j] += BIG_PENALTY * dv_ddpj;
        if (TimeIntegration != 0) {
          d_func[0][var][j] +=
              BIG_PENALTY * t_relax * sign * (sstangent[p] * (1. + 2. * tt) * bf[var]->phi[j] / dt);
          d_func_ss[0][var][j] +=
              BIG_PENALTY * t_relax * sign * (sstangent[p] * (1. + 2. * tt) * bf[var]->phi[j] / dt);
        }
      }
    }
  }

  /*** temperature sensitivities	**/

  var = TEMPERATURE;
  if (pd->v[pg->imtrx][var]) {
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      switch (bc_type) {
      case VELO_THETA_TPL_BC:
        dv_ddpj = dewet * v0 * cosh((costhetaeq - costheta) * g) * (costhetaeq - costheta) *
                  modified_surface_tension * mp->d_surface_tension[var] * bf[var]->phi[j] * factor;
        break;
      case VELO_THETA_HOFFMAN_BC:
        dv_ddpj =
            factor * dewet * (g_dca - g_sca) *
            (liq_visc * modified_surface_tension * mp->d_surface_tension[var] * bf[var]->phi[j] -
             g * d_mu->T[j]) /
            SQUARE(liq_visc);
        break;
      case VELO_THETA_COX_BC:
        f_den = log(reciprocal_slip) + (q_inner / f_sca - q_outer / f_dca);
        dv_ddpj =
            factor * dewet * g_integral / f_den *
            (liq_visc * modified_surface_tension * mp->d_surface_tension[var] * bf[var]->phi[j] -
             g * d_mu->T[j]) /
            SQUARE(liq_visc);
        break;
      case VELO_THETA_SHIK_BC:
        dv_ddpj = factor * 0.5 * v_new / g * modified_surface_tension * mp->d_surface_tension[var] *
                  bf[var]->phi[j];
        break;
      }
      d_func[0][var][j] += BIG_PENALTY * dv_ddpj;
      d_func_ss[0][var][j] += BIG_PENALTY * dv_ddpj;
    }
  }

  /*** concentration sensitivities	**/

  var = MASS_FRACTION;
  if (pd->v[pg->imtrx][var]) {
    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        switch (bc_type) {
        case VELO_THETA_TPL_BC:
          dv_ddpj = dewet * v0 * cosh((costhetaeq - costheta) * g) * (costhetaeq - costheta) *
                    modified_surface_tension * mp->d_surface_tension[MAX_VARIABLE_TYPES + w] *
                    bf[var]->phi[j] * factor;
          break;
        case VELO_THETA_HOFFMAN_BC:
          dv_ddpj = factor * dewet * (g_dca - g_sca) *
                    (liq_visc * modified_surface_tension *
                         mp->d_surface_tension[MAX_VARIABLE_TYPES + w] * bf[var]->phi[j] -
                     g * d_mu->C[w][j]) /
                    SQUARE(liq_visc);
          break;
        case VELO_THETA_COX_BC:
          f_den = log(reciprocal_slip) + (q_inner / f_sca - q_outer / f_dca);
          dv_ddpj = factor * dewet * g_integral / f_den *
                    (liq_visc * modified_surface_tension *
                         mp->d_surface_tension[MAX_VARIABLE_TYPES + w] * bf[var]->phi[j] -
                     g * d_mu->C[w][j]) /
                    SQUARE(liq_visc);
          break;
        case VELO_THETA_SHIK_BC:
          dv_ddpj = factor * 0.5 * v_new / g * modified_surface_tension *
                    mp->d_surface_tension[MAX_VARIABLE_TYPES + w] * bf[var]->phi[j];
          break;
        }
        d_func[0][MAX_VARIABLE_TYPES + w][j] += BIG_PENALTY * dv_ddpj;
        d_func_ss[0][MAX_VARIABLE_TYPES + w][j] += BIG_PENALTY * dv_ddpj;
      }
    }
  }

  return;

} /* END of routine fapply_moving_CA_sinh  */

void fapply_ST(double func[MAX_PDIM],
               double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
               const double tx,     /* x-component of surface tangent at outflow */
               const double ty,     /* y-component of surface tangent at outflow */
               const double tz,     /* z-component of surface tangent at outflow */
               const double sigma,  /*  ST                                       */
               const double rcoord, /* r coord in axisym problems: in r dr dz    */
               const int id_side)

/******************************************************************************
 *
 *  Function which applies a tangent at outflow boundary condition on the
 *  surface stress
 *
 *            Author: P. R. Schunk    (4/28/94)
 *            Revised: 6/1/95 RAC
 *
 ******************************************************************************/

{
  /*    TAB certifies that this function conforms to the exo/patran side numbering convention
   * 11/10/98. */
  int j, var, w;
  int p; /* Degree of freedom counter                   */
  double sign;
  double t[MAX_PDIM];

  /* Function and externals definitions */

  /***************************** EXECUTION BEGINS ******************************/

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  t[0] = tx;
  t[1] = ty;
  t[2] = tz;

  /*First id side and appropriate sign */

  sign = (id_side < 3) ? 1. : -1.;

  /* Calculate the residual contribution					     */

  for (p = 0; p < pd->Num_Dim; p++) {
    func[p] += sign * mp->surface_tension * sigma * t[p] * fv->h3;
    /* does this h3 need to be here?? The radial weighting is already incorporated
     * into the addition into the residual in bc_special
     * -- oops - no its not, so we need it here!!*/
    if (af->Assemble_Jacobian) {
      var = TEMPERATURE;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_func[p][var][j] +=
              sign * mp->d_surface_tension[TEMPERATURE] * sigma * t[p] * bf[var]->phi[j] * fv->h3;
        }
      }
      var = MASS_FRACTION;
      if (pd->v[pg->imtrx][var]) {
        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_func[p][MAX_VARIABLE_TYPES + w][j] += sign *
                                                    mp->d_surface_tension[MAX_VARIABLE_TYPES + w] *
                                                    sigma * t[p] * bf[var]->phi[j] * fv->h3;
          }
        }
      }

      if (pd->CoordinateSystem == CYLINDRICAL || pd->CoordinateSystem == SWIRLING ||
          pd->CoordinateSystem == PROJECTED_CARTESIAN || pd->CoordinateSystem == CARTESIAN_2pt5D) {
        var = MESH_DISPLACEMENT2;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_func[p][var][j] +=
                sign * mp->surface_tension * sigma * t[p] * bf[var]->phi[j] * fv->dh3dq[1];
          }
        }
      }
    }
  }

} /* END of routine fapply_ST  */
/*****************************************************************************/
/*****************************************************************************/

void apply_ST_scalar(double func[MAX_PDIM],
                     double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                     const int id_side,  /* ID of the side of the element        */
                     const double sigma, /* value of the surface tension        */
                     const int userSign) /* User added sign -> see Bug report below   */

/******************************************************************************
 *
 *  Function which applies a tangent at outflow boundary condition on
 *  the surface stress
 *
 *            Author: P. R. Schunk    (4/28/94)
 *            Revised: R. A. Cairncross (6/13/95)
 *
 ******************************************************************************/
/*
 *  HKM Bug report 12/12/2007
 *      I believe this function has no chance of being correct. The sign can not
 * be determined. We call this routine from the element side corresponding to the
 * free surface, because we call CA_ENDFORCE_SCALAR from the CAPILLARY condition logic
 * to apply this where the free surface hits the wall. Therefore, id_side refers to
 * that side. We have no idea in this function where the wall is. Is it to the left or
 * is it to the right? Who knows. We don't have enough information. However, in
 * order to apply the surface tangent condition, we need to pick the surface tangent
 * direction that points out of the domain. (there are two and we need to pick one of
 * the two). Therefore, we need to take the dot product of the wall normal with the
 * surface tangent in order to pick the sign. Currently, we don't have enough
 * information. Therefore, I added a sign input on the card as a starter
 * for handling this.
 */
{
  /*    TAB certifies that this function conforms to the exo/patran side numbering convention
   * 11/10/98. */
  int j, j_id, w;
  int p; /* Degree of freedom counter                   */
  int jvar;
  int var;
  double sign;

  /***************************** EXECUTION BEGINS ******************************/

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  /*First id side and appropriate sign */

  sign = (id_side < 3) ? 1. : -1.;

  if (userSign == -1) {
    sign = -1.0 * sign;
  }

  GOMA_WH(1, "Sign on Surface Tangent Scalar is be incorrect in certain cases");
  if (ei[pg->imtrx]->ielem_dim == 3)
    GOMA_EH(GOMA_ERROR, "Need to update surface tangent for 3D");

  /*
   *  Calculate the residual contribution
   *   -  multiplying by h3 accounts for radial weighting of this
   *      endpoint bc
   */
  for (p = 0; p < ei[pg->imtrx]->ielem_dim; p++) {
    func[p] += sign * mp->surface_tension * sigma * fv->stangent[0][p] * fv->h3;
  }

  /*
   *   Calculate the Jacobian contribution
   */

  if (af->Assemble_Jacobian) {
    for (p = 0; p < ei[pg->imtrx]->ielem_dim; p++) {
      for (jvar = 0; jvar < ei[pg->imtrx]->ielem_dim; jvar++) {
        var = MESH_DISPLACEMENT1 + jvar;
        if (pd->v[pg->imtrx][var]) {
          for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
            d_func[p][var][j_id] += sign * mp->surface_tension * sigma *
                                    (fv->dstangent_dx[0][p][jvar][j_id] * fv->h3 +
                                     fv->stangent[0][p] * fv->dh3dq[jvar] * bf[var]->phi[j_id]);
          }
        }
      }
      var = TEMPERATURE;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_func[p][var][j] += sign * mp->d_surface_tension[TEMPERATURE] * sigma *
                               fv->stangent[0][p] * bf[var]->phi[j] * fv->h3;
        }
      }
      var = MASS_FRACTION;
      if (pd->v[pg->imtrx][var]) {
        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_func[p][MAX_VARIABLE_TYPES + w][j] +=
                sign * mp->d_surface_tension[MAX_VARIABLE_TYPES + w] * sigma * fv->stangent[0][p] *
                bf[var]->phi[j] * fv->h3;
          }
        }
      }
    }
  }
} /* END of routine apply_ST_scalar  */

/*****************************************************************************/
/*****************************************************************************/
void apply_ST_3D(double func[MAX_PDIM],
                 double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                 const double tx,
                 const double ty,
                 const double tz,
                 const double sigma) /*  ST - an equation term multiplier for the *
                                      * capillary stress normally this should     *
                                      * either be 0 or 1                          */

/******************************************************************************
 *
 *  Function which applies a tangent at outflow boundary condition on the surface stress
 *   for 3D edges of interfaces
 *
 *            Author: P. R. Schunk    (4/28/94)
 *            Revised: R. A. Cairncross (6/13/95)
 *            Revised: R. A. Cairncross (7/29/97) - extension to 3D
 *
 ******************************************************************************/

{
  int j, w;
  int p; /* Degree of freedom counter                   */
  int var;
  double t[MAX_PDIM];

  /***************************** EXECUTION BEGINS ******************************/
  if (af->Assemble_LSA_Mass_Matrix)
    return;

  t[0] = tx;
  t[1] = ty;
  t[2] = tz;

  if (ei[pg->imtrx]->ielem_dim != 3)
    GOMA_EH(GOMA_ERROR, "Using 3D surface tangent in 2D");
  /* Calculate the residual contribution					     */

  for (p = 0; p < ei[pg->imtrx]->ielem_dim; p++) {
    func[p] += mp->surface_tension * sigma * t[p];
    /* h3 weighting is accounted for by integral in bc_curve */
  }

  /* Calculate the Jacobian contribution					     */

  if (af->Assemble_Jacobian) {
    for (p = 0; p < ei[pg->imtrx]->ielem_dim; p++) {
      var = TEMPERATURE;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_func[p][var][j] += mp->d_surface_tension[TEMPERATURE] * sigma * t[p] * bf[var]->phi[j];
        }
      }
      var = MASS_FRACTION;
      if (pd->v[pg->imtrx][var]) {
        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_func[p][MAX_VARIABLE_TYPES + w][j] +=
                mp->d_surface_tension[MAX_VARIABLE_TYPES + w] * sigma * t[p] * bf[var]->phi[j];
          }
        }
      }
    }
  }

} /* END of routine apply_ST_3D  */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void apply_ST_scalar_3D(double func[MAX_PDIM],
                        double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                        const double sigma) /* ST - an equation term multiplier   *
                                             * for the capillary stress -         *
                                             * normally this should               *
                                             * either be 0 or 1                   */

/******************************************************************************
 *
 *  Function which applies a tangent at outflow boundary condition on the surface stress
 *   for 3D edges of interfaces
 *
 *            Author: P. R. Schunk    (4/28/94)
 *            Revised: R. A. Cairncross (6/13/95)
 *            Revised: R. A. Cairncross (7/29/97) - extension to 3D
 *
 ******************************************************************************/

{
  int j, j_id, w;
  int p; /* Degree of freedom counter                   */
  int jvar;
  int var;

  /***************************** EXECUTION BEGINS ******************************/

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  if (ei[pg->imtrx]->ielem_dim != 3)
    GOMA_EH(GOMA_ERROR, "Using 3D surface tangent in 2D");
  /* Calculate the residual contribution					     */

  for (p = 0; p < ei[pg->imtrx]->ielem_dim; p++) {
    func[p] += mp->surface_tension * sigma * fv->stangent[0][p];
    /* h3 weighting is accounted for by integral in bc_curve */
  }

  /* Calculate the Jacobian contribution					     */

  if (af->Assemble_Jacobian) {
    for (p = 0; p < ei[pg->imtrx]->ielem_dim; p++) {
      for (jvar = 0; jvar < ei[pg->imtrx]->ielem_dim; jvar++) {
        var = MESH_DISPLACEMENT1 + jvar;
        for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
          d_func[p][var][j_id] += mp->surface_tension * sigma * fv->dstangent_dx[0][p][jvar][j_id];
        }
      }

      var = TEMPERATURE;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_func[p][var][j] +=
              mp->d_surface_tension[TEMPERATURE] * sigma * fv->stangent[0][p] * bf[var]->phi[j];
        }
      }
      var = MASS_FRACTION;
      if (pd->v[pg->imtrx][var]) {
        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_func[p][MAX_VARIABLE_TYPES + w][j] += mp->d_surface_tension[MAX_VARIABLE_TYPES + w] *
                                                    sigma * fv->stangent[0][p] * bf[var]->phi[j];
          }
        }
      }
    }
  }

} /* END of routine apply_ST_scalar_3D  */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void apply_CA_FILL(double func[MAX_PDIM],
                   double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                   const double ca) /* contact angle                             */

/******************************************************************************
 *
 *  Function which applies bc on momentum for fill ca
 *
 *            Author: David Noble (4/27/2000)
 *
 ******************************************************************************/

{
  int j, var;
  int a, b, p; /* Degree of freedom counter                   */
  double cos_ca, sin_ca, sign;
  double nw[DIM], nf[DIM];
  double tmag;
  int dim = pd->Num_Dim;
  double t[DIM], dot_prod;
  double dt_dmesh[DIM][DIM][MDE];
  double dp;
  double ddp_dFj, dtmag_dFj, dt_dFj[DIM];
  double ddp_dmeshbj, dtmag_dmeshbj, dt_dmeshbj[DIM];

  /***************************** EXECUTION BEGINS ******************************/

  if (af->Assemble_LSA_Mass_Matrix)
    return;

#if 1

  /* Fetch the level set interface quantities. Bail if we're not in the mushy zone.*/
  load_lsi(ls->Length_Scale);
  if (!lsi->near)
    return;

  /* Get the derivatives if necessary. */
  if (af->Assemble_Jacobian)
    load_lsi_derivs();

  /* Some shortcuts. */
  for (p = 0; p < dim; p++) {
    nw[p] = fv->snormal[p]; /* Normal to the "wall" */
    nf[p] = lsi->normal[p]; /* Normal to the level set (F=0) */
  }
  cos_ca = cos(M_PIE * ca / 180.0);
  sin_ca = sin(M_PIE * ca / 180.0);

  /*
   * PKN: I think the 3D case is general enough that it should work in
   * 2D as well.  It might require slightly more CPU, but could make
   * maintenance easier.
   */

  if (dim == 2) {

    sign = 1.0;
    if ((lsi->normal[0] * fv->snormal[1] - lsi->normal[1] * fv->snormal[0]) < 0.0)
      sign = -1.0;

    t[0] = sign * fv->snormal[1];
    t[1] = -sign * fv->snormal[0];

    dot_prod = dot_product(dim, t, lsi->normal);

    for (p = 0; p < dim; p++) {
      func[p] =
          lsi->delta * mp->surface_tension * dot_prod * (sin_ca * fv->snormal[p] + cos_ca * t[p]);
    }
    if (!af->Assemble_Jacobian)
      return;
    /* Compute the jacobian contributions. */

    /* Derivatives w.r.t FILL */
    var = FILL;
    for (p = 0; p < dim; p++) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        /* Fetch the basis functions and derivatives. */

        /* Compute the derivatives term-by-term; */
        d_func[p][var][j] = lsi->d_delta_dF[j] * mp->surface_tension * dot_prod *
                            (sin_ca * fv->snormal[p] + cos_ca * t[p]);

        d_func[p][var][j] += lsi->delta * mp->surface_tension *
                             (t[0] * lsi->d_normal_dF[0][j] + t[1] * lsi->d_normal_dF[1][j]) *
                             (sin_ca * fv->snormal[p] + cos_ca * t[p]);

      } /* for: j=0,...,ei[pg->imtrx]->dof[FILL] */
    }   /* for: p=0,...,dim */

    /* Derivatives w.r.t MESH */
    if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
      /* Precompute dt_dmesh[][][]. */
      for (b = 0; b < dim; b++) {
        var = MESH_DISPLACEMENT1 + b;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          dt_dmesh[0][b][j] = +sign * fv->dsnormal_dx[1][b][j];
          dt_dmesh[1][b][j] = -sign * fv->dsnormal_dx[0][b][j];
        }
      }

      for (p = 0; p < dim; p++) {
        for (b = 0; b < dim; b++) {
          var = MESH_DISPLACEMENT1 + b;
          /* dsnormal_dx[DIM][DIM][MDE] */
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_func[p][var][j] = lsi->d_delta_dmesh[b][j] * mp->surface_tension * dot_prod *
                                (sin_ca * fv->snormal[p] + cos_ca * t[p]);

            for (a = 0; a < VIM; a++) {
              d_func[p][var][j] +=
                  lsi->delta * mp->surface_tension *
                  (dt_dmesh[a][b][j] * lsi->normal[a] + t[a] * lsi->d_normal_dmesh[a][b][j]) *
                  (sin_ca * fv->snormal[p] + cos_ca * t[p]);
            }
            d_func[p][var][j] += lsi->delta * mp->surface_tension * dot_prod *
                                 (sin_ca * fv->dsnormal_dx[p][b][j] + cos_ca * dt_dmesh[p][b][j]);

          } /* for: j */

        } /* for: b */

      } /* for: p */

    } /* if: Moving Mesh */

  } /* if: 2D */
  else {

    /* We can simplify the vector math ahead of time.
     *
     * l = (nf x nw) / | nf x nw |
     *
     * t = nw x l / | nw x l |
     *   = nw x ( nf x nw ) / | nw x ( nf x nw ) |
     *
     * nw x ( nf x nw ) = nf - (nw.nf) nw    (using the "ed" rule)
     *
     * t = (nf - nw nw.nf) / | nf - nw nw.nf |
     *
     * | nf - nw nw.nf | = [ 1  - 2 (nw.nf)^2 ]^1/2
     *
     * t_i  = [ 1 + nw.nf - 2 (nw.nf)^2 ]^-1/2 * (nf_i - nw_i nw_j nf_j)
     *
     * Since we only need "t" and "nf", I won't actually compute "l".
     */

    dp = dot_product(dim, nf, nw);

    tmag = sqrt(1.0 - pow(dp, 2.0));
    for (p = 0; p < dim; p++) {
      t[p] = (nf[p] - dp * nw[p]) / tmag;
    }
    dot_prod = dot_product(dim, t, nf);

    for (p = 0; p < dim; p++) {
      func[p] =
          lsi->delta * mp->surface_tension * dot_prod * (sin_ca * fv->snormal[p] + cos_ca * t[p]);
    }

    if (!af->Assemble_Jacobian)
      return;

    /* Derivatives w.r.t. FILL */
    var = FILL;
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      /* Precompute ddp_dF */
      for (ddp_dFj = 0.0, a = 0; a < dim; a++)
        ddp_dFj += lsi->d_normal_dF[a][j] * nw[a];

      /* Precompute dtmag_dF */
      dtmag_dFj = (0.5 / tmag) * (-2.0 * dp * ddp_dFj);

      /* Precompute dt_dF[][] */
      for (a = 0; a < dim; a++)
        dt_dFj[a] = (lsi->d_normal_dF[a][j] - ddp_dFj * nw[a]) / tmag +
                    (nf[a] - dp * nw[a]) * (-pow(tmag, -2)) * dtmag_dFj;

      for (p = 0; p < dim; p++) {

        d_func[p][var][j] = lsi->d_delta_dF[j] * mp->surface_tension * dot_prod *
                            (sin_ca * fv->snormal[p] + cos_ca * t[p]);

        for (a = 0; a < dim; a++) {
          d_func[p][var][j] = lsi->delta * mp->surface_tension *
                              (dt_dFj[a] * nf[a] + t[a] * lsi->d_normal_dF[a][j]) *
                              (sin_ca * fv->snormal[p] + cos_ca * t[p]);
        } /* for: a */

        d_func[p][var][j] = lsi->delta * mp->surface_tension * dot_prod * (cos_ca * dt_dFj[p]);

      } /* for: p */

    } /* for: j */

    /* Derivatives w.r.t. MESH_DISPLACEMENTs */
    if (!pd->v[pg->imtrx][MESH_DISPLACEMENT1])
      return;
    for (b = 0; b < dim; b++) {
      var = MESH_DISPLACEMENT1 + b;
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        /* Precompute ddp_dmeshbj */
        for (a = 0; a < dim; a++) {
          ddp_dmeshbj = lsi->d_normal_dmesh[a][b][j] * nw[a] + nf[a] * fv->dsnormal_dx[a][b][j];
        }

        /* Precompute dtmag_dmeshbj */
        for (a = 0; a < dim; a++) {
          dtmag_dmeshbj = (0.5 / tmag) * (-ddp_dmeshbj - 4.0 * dp * ddp_dmeshbj);
        }

        /* Precompute dt_dmeshbj[] */
        for (a = 0; a < dim; a++) {
          dt_dmeshbj[a] =
              (lsi->d_normal_dmesh[a][b][j] - ddp_dmeshbj * nw[a] - dp * fv->dsnormal_dx[a][b][j]) /
                  tmag +
              (nf[a] - dp * nw[a]) * (-pow(tmag, -2.0)) * dtmag_dmeshbj;
        }

        for (p = 0; p < dim; p++) {
          d_func[p][var][j] = lsi->d_delta_dmesh[b][j] * mp->surface_tension * dot_prod *
                              (sin_ca * fv->snormal[p] + cos_ca * t[p]);

          for (a = 0; a < dim; a++) {
            d_func[p][var][j] = lsi->delta * mp->surface_tension *
                                (dt_dmeshbj[a] * nf[a] + t[a] * lsi->d_normal_dmesh[a][b][j]) *
                                (sin_ca * fv->snormal[p] + cos_ca * t[p]);
          }

          d_func[p][var][j] = lsi->delta * mp->surface_tension * dot_prod *
                              (sin_ca * fv->dsnormal_dx[p][b][j] + cos_ca * dt_dmeshbj[p]);

        } /* for: p */

      } /* for: j */

    } /* for: b */

  } /* else: 3D */

#endif /* 1 */

} /* END of routine apply_CA_FILL  */

void apply_sharp_ca(double func[MAX_PDIM],
                    double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                    const double ca) /* contact angle                             */
                                     /*
                                      *
                                      *   Apply wetting line forces as "Sharp" conditions *
                                      *
                                      */

{
  int j, var;
  int p; /* Degree of freedom counter                   */
  double cos_ca, sin_ca, sign;
  int dim = pd->Num_Dim;
  double t[DIM], dot_prod;

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  load_lsi(ls->Length_Scale);

  /* Get the derivatives if necessary. */
  if (af->Assemble_Jacobian)
    load_lsi_derivs();

  cos_ca = cos(M_PIE * ca / 180.0);
  sin_ca = sin(M_PIE * ca / 180.0);

  /*
   * PKN: I think the 3D case is general enough that it should work in
   * 2D as well.  It might require slightly more CPU, but could make
   * maintenance easier.
   */

  if (dim == 2) {
    double d_dot_prod_dF[MDE];

    sign = 1.0;
    if ((lsi->normal[0] * fv->snormal[1] - lsi->normal[1] * fv->snormal[0]) < 0.0)
      sign = -1.0;

    t[0] = sign * fv->snormal[1];
    t[1] = -sign * fv->snormal[0];

    dot_prod = dot_product(dim, t, lsi->normal);

    for (p = 0; p < dim; p++) {
      func[p] = mp->surface_tension * dot_prod * (sin_ca * fv->snormal[p] + cos_ca * t[p]);
    }
    if (af->Assemble_Jacobian) {
      for (j = 0; j < ei[pg->imtrx]->dof[FILL]; j++) {
        d_dot_prod_dF[j] = 0.0;

        for (p = 0; p < dim; p++) {
          d_dot_prod_dF[j] += t[p] * lsi->d_normal_dF[p][j];
        }
      }

      for (p = 0; p < dim; p++) {

        var = FILL;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_func[p][var][j] =
              mp->surface_tension * d_dot_prod_dF[j] * (sin_ca * fv->snormal[p] + cos_ca * t[p]);
        }
      }
    }
  } else {
    GOMA_EH(
        GOMA_ERROR,
        "Hello!....SHARP_CA_2D.   2D!  2D!      Doesn't work for three dimensional problems.  \n");
  }
  return;
}

void apply_wetting_velocity(double func[MAX_PDIM],
                            double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                            int id_side,
                            int element_type,
                            double depth,
                            double vs[MAX_PDIM],
                            double wetting_length,
                            double theta_s,
                            double cT) {
  int a, b, j, dim = pd->Num_Dim, var;
  double betainv = 1.e12;
  double wet_speed, wet_vector[MAX_PDIM], d_wet_vector_dF[MAX_PDIM][MDE];
  double d_wet_vector_dX[MAX_PDIM][MAX_PDIM][MDE];

  if (ls == NULL) {
    GOMA_EH(GOMA_ERROR,
            "Boundary condition WETTING_SPEED_LINEAR requires active level set dofs.\n");
  }

  /*	load_lsi(ls->Length_Scale);	*/

  if (wetting_length <= 0.0)
    wetting_length = ls->Length_Scale;

  load_lsi(wetting_length);

  for (a = 0; a < MAX_PDIM; a++)
    wet_vector[a] = 0.0;
  memset(d_wet_vector_dF, 0, sizeof(double) * MAX_PDIM * MDE);
  memset(d_wet_vector_dX, 0, sizeof(double) * MAX_PDIM * MAX_PDIM * MDE);

  betainv = 1.0 / depth;

  if (lsi->near) {
    double nw[MAX_PDIM], nf[MAX_PDIM], t[MAX_PDIM];
    double d_t_dFj[MAX_PDIM];
    double dp = 0.0, tmag = 0.0, tmaginv = 0.0;
    double d_dp_dFj, d_tmag_dFj;
    double cos_ca, cos_ca_static;
    double delta = 0, d_delta_dFj = 0, d_delta_dXj[DIM];
    double d_dp_dXj[DIM], d_t_dXj[DIM][DIM], d_tmag_dXj[DIM];

    memset(t, 0, sizeof(double) * MAX_PDIM);
    memset(d_t_dFj, 0, sizeof(double) * MAX_PDIM);
    memset(d_t_dXj, 0, sizeof(double) * DIM * DIM);

    cos_ca_static = cos(M_PIE * (180.0 - theta_s) / 180.0);

    if (af->Assemble_Jacobian)
      load_lsi_derivs();

    for (a = 0; a < dim; a++) {
      nw[a] = fv->snormal[a];
      nf[a] = lsi->normal[a];

      dp += nw[a] * nf[a];
    }

    tmag = sqrt(1.0 - pow(dp, 2.0));
    if (tmag <= DBL_SMALL) {
      printf("Underflow tmag detected.\n");
      tmaginv = 0.0;
    } else
      tmaginv = 1.0 / tmag;

    for (a = 0; a < dim; a++) {
      t[a] = (nf[a] - dp * nw[a]) * tmaginv;
    }

    cos_ca = dp;

    wet_speed = (cos_ca - cos_ca_static) / cT;

    if (wetting_length < 1.e-12)
      wetting_length = 1.;

    delta = lsi->delta;

    for (a = 0; a < dim; a++) {
      wet_vector[a] = t[a] * wet_speed * delta;
    }

    betainv = 1.0 / depth;
    var = ls->var;

    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      for (d_dp_dFj = 0, a = 0; a < dim; a++)
        d_dp_dFj += lsi->d_normal_dF[a][j] * nw[a];

      d_tmag_dFj = -tmaginv * dp * d_dp_dFj;

      for (a = 0, d_t_dFj[a] = 0; a < dim; a++)
        d_t_dFj[a] = (nf[a] - dp * nw[a]) * (-pow(tmaginv, 2)) * d_tmag_dFj +
                     (lsi->d_normal_dF[a][j] - d_dp_dFj * nw[a]) * tmaginv;

      d_delta_dFj = lsi->d_delta_dF[j];

      for (a = 0, d_wet_vector_dF[a][j] = 0.0; a < dim; a++)
        d_wet_vector_dF[a][j] = d_t_dFj[a] * wet_speed * delta + t[a] * d_dp_dFj * delta / cT +
                                t[a] * wet_speed * d_delta_dFj;
    }
    var = MESH_DISPLACEMENT1;
    if (pd->v[pg->imtrx][var]) {
      for (b = 0; b < dim; b++) {
        var = MESH_DISPLACEMENT1 + b;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (d_dp_dXj[b] = 0.0, a = 0; a < dim; a++) {
            d_dp_dXj[b] += lsi->d_normal_dmesh[a][b][j] * nw[a];
          }

          d_tmag_dXj[b] = -tmaginv * dp * d_dp_dXj[b];

          for (a = 0, d_t_dXj[a][b] = 0.0; a < dim; a++) {
            d_t_dXj[a][b] = (nf[a] - dp * nw[a]) * (-pow(tmaginv, 2)) * d_tmag_dXj[b] +
                            (lsi->d_normal_dmesh[a][b][j] - d_dp_dXj[b] * nw[a]) * tmaginv;
          }

          d_delta_dXj[b] = lsi->d_delta_dmesh[b][j];

          for (a = 0; a < dim; a++) {
            d_wet_vector_dX[a][b][j] = d_t_dXj[a][b] * wet_speed * delta +
                                       t[a] * d_dp_dXj[b] * delta / cT +
                                       t[a] * wet_speed * d_delta_dXj[b];
          }
        }
      }
    }
  }

  for (a = 0; a < pd->Num_Dim; a++) {
    func[a] = betainv * (wet_vector[a]);
  }

  if (af->Assemble_Jacobian) {
    var = ls->var;

    for (a = 0; a < dim; a++) {
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_func[a][var][j] += (betainv) * (d_wet_vector_dF[a][j]);
        }
      }
    }
    var = MESH_DISPLACEMENT1;

    if (pd->v[pg->imtrx][var]) {
      for (b = 0; b < dim; b++) {
        var = MESH_DISPLACEMENT1 + b;

        for (a = 0; a < dim; a++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_func[a][var][j] += betainv * d_wet_vector_dX[a][b][j];
          }
        }
      }
    }
  }

  return;
}

void apply_linear_wetting_sic(double func[MAX_PDIM],
                              double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                              const double dt,
                              const double tt,
                              int id_side,
                              int element_type,
                              double beta,
                              double vs[MAX_PDIM],
                              double wetting_length,
                              double theta_s,
                              double cT,
                              double tau) {
  int a, b, j, dim = pd->Num_Dim, var;
  double betainv = 1.e12;
  double wet_speed, wet_vector[MAX_PDIM], d_wet_vector_dF[MAX_PDIM][MDE];
  double d_wet_vector_dX[MAX_PDIM][MAX_PDIM][MDE];

  if (ls == NULL) {
    GOMA_EH(GOMA_ERROR,
            "Boundary condition WETTING_SPEED_LINEAR requires active level set dofs.\n");
  }

  /*	load_lsi(ls->Length_Scale);	*/

  if (wetting_length <= 0.0)
    wetting_length = ls->Length_Scale;

  load_lsi(wetting_length);

  for (a = 0; a < MAX_PDIM; a++)
    wet_vector[a] = 0.0;
  memset(d_wet_vector_dF, 0, sizeof(double) * MAX_PDIM * MDE);
  memset(d_wet_vector_dX, 0, sizeof(double) * MAX_PDIM * MAX_PDIM * MDE);

  betainv = 1.0 / beta;

  if (lsi->near) {
    double nw[MAX_PDIM], nf[MAX_PDIM], t[MAX_PDIM];
    double d_t_dFj[MAX_PDIM];
    double dp = 0.0, tmag = 0.0, tmaginv = 0.0;
    double d_dp_dFj, d_tmag_dFj;
    double cos_ca, cos_ca_static;
    double delta = 0, d_delta_dFj = 0, d_delta_dXj[DIM];
    double d_dp_dXj[DIM], d_t_dXj[DIM][DIM], d_tmag_dXj[DIM];

    memset(t, 0, sizeof(double) * MAX_PDIM);
    memset(d_t_dFj, 0, sizeof(double) * MAX_PDIM);
    memset(d_t_dXj, 0, sizeof(double) * DIM * DIM);

    cos_ca_static = cos(M_PIE * (180.0 - theta_s) / 180.0);

    if (af->Assemble_Jacobian)
      load_lsi_derivs();

    for (a = 0; a < dim; a++) {
      nw[a] = fv->snormal[a];
      nf[a] = lsi->normal[a];

      dp += nw[a] * nf[a];
    }

    tmag = sqrt(1.0 - pow(dp, 2.0));
    if (tmag <= DBL_SMALL) {
      printf("Underflow tmag detected.\n");
      tmaginv = 0.0;
    } else
      tmaginv = 1.0 / tmag;

    for (a = 0; a < dim; a++) {
      t[a] = (nf[a] - dp * nw[a]) * tmaginv;
    }

    cos_ca = dp;

    wet_speed = (cos_ca - cos_ca_static) / cT;

    if (wetting_length < 1.e-12)
      wetting_length = 1.;

    delta = 0.5 * lsi->delta * wetting_length;

    for (a = 0; a < dim; a++) {
      wet_vector[a] = t[a] * wet_speed * delta;
    }

    betainv = 1.0 / beta;
    var = ls->var;

    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      for (d_dp_dFj = 0, a = 0; a < dim; a++)
        d_dp_dFj += lsi->d_normal_dF[a][j] * nw[a];

      d_tmag_dFj = -tmaginv * dp * d_dp_dFj;

      for (a = 0, d_t_dFj[a] = 0; a < dim; a++)
        d_t_dFj[a] = (nf[a] - dp * nw[a]) * (-pow(tmaginv, 2)) * d_tmag_dFj +
                     (lsi->d_normal_dF[a][j] - d_dp_dFj * nw[a]) * tmaginv;

      d_delta_dFj = 0.5 * lsi->d_delta_dF[j] * wetting_length;

      for (a = 0, d_wet_vector_dF[a][j] = 0.0; a < dim; a++)
        d_wet_vector_dF[a][j] = d_t_dFj[a] * wet_speed * delta + t[a] * d_dp_dFj * delta / cT +
                                t[a] * wet_speed * d_delta_dFj;
    }
    var = MESH_DISPLACEMENT1;
    if (pd->v[pg->imtrx][var]) {
      for (b = 0; b < dim; b++) {
        var = MESH_DISPLACEMENT1 + b;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (d_dp_dXj[b] = 0.0, a = 0; a < dim; a++) {
            d_dp_dXj[b] += lsi->d_normal_dmesh[a][b][j] * nw[a];
          }

          d_tmag_dXj[b] = -tmaginv * dp * d_dp_dXj[b];

          for (a = 0, d_t_dXj[a][b] = 0.0; a < dim; a++) {
            d_t_dXj[a][b] = (nf[a] - dp * nw[a]) * (-pow(tmaginv, 2)) * d_tmag_dXj[b] +
                            (lsi->d_normal_dmesh[a][b][j] - d_dp_dXj[b] * nw[a]) * tmaginv;
          }

          d_delta_dXj[b] = 0.5 * lsi->d_delta_dmesh[b][j] * wetting_length;

          for (a = 0; a < dim; a++) {
            d_wet_vector_dX[a][b][j] = d_t_dXj[a][b] * wet_speed * delta +
                                       t[a] * d_dp_dXj[b] * delta / cT +
                                       t[a] * wet_speed * d_delta_dXj[b];
          }
        }
      }
    }
  }

  for (a = 0; a < pd->Num_Dim; a++) {
    func[a] = betainv * (-tau * fv_dot->v[a] + wet_vector[a] + vs[a] - fv->v[a]);
  }

  if (af->Assemble_Jacobian) {
    var = ls->var;

    for (a = 0; a < dim; a++) {
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_func[a][var][j] += (betainv) * (d_wet_vector_dF[a][j]);
        }
      }
    }
    var = MESH_DISPLACEMENT1;

    if (pd->v[pg->imtrx][var]) {
      for (b = 0; b < dim; b++) {
        var = MESH_DISPLACEMENT1 + b;

        for (a = 0; a < dim; a++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_func[a][var][j] += betainv * d_wet_vector_dX[a][b][j];
          }
        }
      }
    }

    var = VELOCITY1;
    if (pd->v[pg->imtrx][var]) {
      double phi_j;

      for (b = 0; b < dim; b++) {
        var = VELOCITY1 + b;

        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          d_func[b][var][j] += betainv * (-tau * (1. + 2. * tt) * phi_j / dt - phi_j);
        }
      }
    }
  }
  return;
}

void apply_sharp_wetting_velocity(double func[MAX_PDIM],
                                  double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                                  const int bc_type,
                                  const double depth,
                                  const double theta_s,
                                  const double cT,
                                  const double v0,
                                  const double t_relax,
                                  const double v_old) {
  int a, b, j, dim = pd->Num_Dim, var, q, w;
  double betainv = 1.e12;
  double wet_speed = 0.0, wet_vector[MAX_PDIM], d_wet_vector_dF[MAX_PDIM][MDE];
  double nw[MAX_PDIM], nf[MAX_PDIM], t[MAX_PDIM];
  double d_wet_vector_dX[MAX_PDIM][MAX_PDIM][MDE];
  double d_wet_vector_dpj = 0.0;
  double d_t_dFj[MAX_PDIM];
  double d_wet_speed_dFj = 0.0, d_wet_speed_dXj = 0.0;
  double dp = 0.0, tmag = 0.0, tmaginv = 0.0;
  double d_dp_dFj, d_tmag_dFj;
  double d_dp_dXj[DIM], d_t_dXj[DIM][DIM], d_tmag_dXj[DIM];
  double cos_ca, cos_ca_static;
  double factor;
  int elem_sign_org;

  int include_stress = FALSE;
  double Pi[DIM][DIM];
  STRESS_DEPENDENCE_STRUCT d_Pi_struct;
  STRESS_DEPENDENCE_STRUCT *d_Pi = &d_Pi_struct;

  /*  Blake condition parameters	*/
  double g;

  /*  Hoffman correlation variables       */
  double ca_no, g_sca = 0.0, g_dca = 0.0, A_sca = 0.0, A_dca = 0.0;
  int iter, iter_max = 20;
  double eps_tol = 1.0e-12;
  double liq_visc = 0.0, gamma[DIM][DIM];
  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct; /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;

  /*  Cox wetting parameters      */
  double theta = 0.0, thetaeq = 0.0, sintheta = 0.0, sinthetaeq = 0.0;
  double costheta = 0.0, costhetaeq = 0.0;
  double lambda = 0;
  double th, sinth, costh;
  double f_sca = 0.0, f_dca = 0.0, f_num = 0.0, f_den = 0.0, g_integral = 0.0;
  double reciprocal_slip = 0.0;
  const double q_inner = 0;
  const double q_outer = 0;
  /* use 10 point Gaussian quadrature for now     */
  const int num_gauss_pts = 10;
  const double gpt[10] = {-0.973906528517172, -0.865063366688985, -0.679409568299024,
                          -0.433395394129247, -0.148874338981631, 0.148874338981631,
                          0.433395394129247,  0.679409568299024,  0.865063366688985,
                          0.973906528517172};
  const double wt[10] = {0.066671344308688, 0.149451349150581, 0.219086362515982, 0.269266719309996,
                         0.295524224714753, 0.295524224714753, 0.269266719309996, 0.219086362515982,
                         0.149451349150581, 0.066671344308688};

  /* Shikhmurzaev wetting parameters	*/
  double rhs = 0.0, rhs_den = 0.0, rhs_num = 0.0, veloc0 = 0.0;
  double drhs_ddpj, drhs_den_ddpj, drhs_num_ddpj, dveloc0_ddpj;
  const double psg_surf = 0;

  if (!af->Assemble_Jacobian)
    d_Pi = NULL;

  memset(t, 0, sizeof(double) * MAX_PDIM);
  memset(d_t_dFj, 0, sizeof(double) * MAX_PDIM);

  cos_ca_static = cos(M_PIE * (180.0 - theta_s) / 180.0);
  costhetaeq = -cos_ca_static;

  if (ls == NULL) {
    GOMA_EH(GOMA_ERROR, "Sharp LS wetting Boundary conditions requires active level set dofs.\n");
  }

  for (a = 0; a < MAX_PDIM; a++)
    wet_vector[a] = 0.0;
  memset(d_wet_vector_dF, 0, sizeof(double) * MAX_PDIM * MDE);
  memset(d_wet_vector_dX, 0, sizeof(double) * MAX_PDIM * MAX_PDIM * MDE);

  g = cT * mp->surface_tension;

  if (bc_type == SHARP_HOFFMAN_VELOCITY_BC || bc_type == SHARP_COX_VELOCITY_BC) {
    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
      }
    }
    elem_sign_org = ls->Elem_Sign;
    switch (mp->mp2nd->viscositymask[1] - mp->mp2nd->viscositymask[0]) {
    case 1:
      ls->Elem_Sign = -1;
      break;
    case -1:
      ls->Elem_Sign = 1;
      break;
    }
    liq_visc = viscosity(gn, gamma, d_mu);
#if 0
      printf("sharp %d %g %g %g %g\n",ls->Elem_Sign, lsi->H, liq_visc, gn->mu0, mp->mp2nd->viscosity);
#endif
    ls->Elem_Sign = elem_sign_org;
  }

  for (a = 0; a < dim; a++) {
    nw[a] = fv->snormal[a];
    nf[a] = lsi->normal[a];

    dp += nw[a] * nf[a];
  }

  tmag = sqrt(1.0 - pow(dp, 2.0));
  if (tmag <= DBL_SMALL) {
    printf("Underflow tmag detected.\n");
    tmaginv = 0.0;
  } else
    tmaginv = 1.0 / tmag;

  for (a = 0; a < dim; a++) {
    t[a] = (nf[a] - dp * nw[a]) * tmaginv;
  }

  cos_ca = dp;
  costheta = -cos_ca;

  switch (bc_type) {
  case SHARP_WETLIN_VELOCITY_BC:
    wet_speed = (cos_ca - cos_ca_static) / g;
    break;
  case SHARP_BLAKE_VELOCITY_BC:
    wet_speed = v0 * sinh(g * (cos_ca - cos_ca_static));
    break;
  case SHARP_HOFFMAN_VELOCITY_BC:
    ca_no = 1.0E+06;
    iter = 0;
    eps = 10. * eps_tol;
    while (iter <= iter_max && fabs(eps) > eps_tol) {
      g_sca = log((3. - costhetaeq) / (1. + costhetaeq)) / (2. * 5.16);
      A_sca = pow(g_sca, 1. / 0.706);
      eps = -(ca_no - 1.31 * pow(ca_no, 0.99) * A_sca - A_sca) /
            (1. - 1.31 * 0.99 * A_sca / pow(ca_no, 0.01));
      ca_no += eps;
      iter++;
    }
    if (fabs(eps) > eps_tol)
      GOMA_EH(GOMA_ERROR, "Hoffman iteration not converged");
    g_sca = ca_no;
    ca_no = 1.0E+06;
    iter = 0;
    eps = 10. * eps_tol;
    while (iter <= iter_max && fabs(eps) > eps_tol) {
      g_dca = log((3. - costheta) / (1. + costheta)) / (2. * 5.16);
      A_dca = pow(g_dca, 1. / 0.706);
      eps = -(ca_no - 1.31 * pow(ca_no, 0.99) * A_dca - A_dca) /
            (1. - 1.31 * 0.99 * A_dca / pow(ca_no, 0.01));
      ca_no += eps;
      iter++;
    }
    if (fabs(eps) > eps_tol)
      GOMA_EH(GOMA_ERROR, "Hoffman iteration not converged");
    g_dca = ca_no;
    ca_no = g_dca - g_sca;
    wet_speed = ca_no * g / liq_visc;
    break;
  case SHARP_COX_VELOCITY_BC:
    thetaeq = M_PIE * theta_s / 180.0;
    theta = acos(costheta);
    sintheta = sin(theta);
    sinthetaeq = sin(thetaeq);
    reciprocal_slip = 1. / v0;

    /*  Cox analysis integral       */
    g_integral = 0;
    for (j = 0; j < num_gauss_pts; j++) {
      th = thetaeq + (theta - thetaeq) * (gpt[j] + 1.) / 2.;
      sinth = sin(th);
      costh = cos(th);
      f_num = 2. * sinth *
              (SQUARE(lambda) * (SQUARE(th) - SQUARE(sinth)) +
               2. * lambda * (th * (M_PIE - th) + SQUARE(sinth)) +
               (SQUARE(M_PIE - th) - SQUARE(sinth)));
      f_den = lambda * (SQUARE(th) - SQUARE(sinth)) * (M_PIE - th + sinth * costh) +
              (SQUARE(M_PIE - th) - SQUARE(sinth)) * (th - sinth * costh);
      g_integral += wt[j] * f_den / f_num;
    }
    g_integral *= 0.5 * (theta - thetaeq);
    /*  inner and outer solution terms      */
    f_num = 2. * sinthetaeq *
            (SQUARE(lambda) * (SQUARE(thetaeq) - SQUARE(sinthetaeq)) +
             2. * lambda * (thetaeq * (M_PIE - thetaeq) + SQUARE(sinthetaeq)) +
             (SQUARE(M_PIE - thetaeq) - SQUARE(sinthetaeq)));
    f_den = lambda * (SQUARE(thetaeq) - SQUARE(sinthetaeq)) *
                (M_PIE - thetaeq + sinthetaeq * costhetaeq) +
            (SQUARE(M_PIE - thetaeq) - SQUARE(sinthetaeq)) * (thetaeq - sinthetaeq * costhetaeq);
    f_sca = f_num / f_den;
    f_num = 2. * sintheta *
            (SQUARE(lambda) * (SQUARE(theta) - SQUARE(sintheta)) +
             2. * lambda * (theta * (M_PIE - theta) + SQUARE(sintheta)) +
             (SQUARE(M_PIE - theta) - SQUARE(sintheta)));
    f_den = lambda * (SQUARE(theta) - SQUARE(sintheta)) * (M_PIE - theta + sintheta * costheta) +
            (SQUARE(M_PIE - theta) - SQUARE(sintheta)) * (theta - sintheta * costheta);
    f_dca = f_num / f_den;
    /* solve for wetting speed      */
    ca_no = g_integral / (log(reciprocal_slip) + (q_inner / f_sca - q_outer / f_dca));
    wet_speed = ca_no * g / liq_visc;
    break;
  case SHARP_SHIK_VELOCITY_BC:
    theta = acos(costheta);
    sintheta = sin(theta);
    veloc0 = (sintheta - theta * costheta) / (sintheta * costheta - theta);
    rhs_num = costhetaeq - costheta;
    rhs_den = psg_surf + v0 + (v0 - 1.) * veloc0 + costheta;
    rhs = rhs_num / rhs_den;
    wet_speed = sqrt(g * v0) * rhs / (2. * sqrt(1. + rhs));
    break;
  default:
    GOMA_EH(GOMA_ERROR, "bad DCA bc name\n");
    break;
  }
#if 0
  printf("Ca angle vnew %g %g %g\n",ca_no, acos(costheta)*180/M_PIE, wet_speed);
#endif

  if (t_relax <= 0.) {
    factor = 1;
  } else if (tran->time_value / t_relax > 36) {
    factor = 1;
  } else {
    factor = 1 - exp(-tran->time_value / t_relax);
  }

  if (TimeIntegration != 0) {
    wet_speed = v_old + (wet_speed - v_old) * factor;
  }

  for (a = 0; a < dim; a++) {
    wet_vector[a] = t[a] * wet_speed;
  }
  betainv = 1.0 / depth;

  var = ls->var;

  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
    for (d_dp_dFj = 0, a = 0; a < dim; a++)
      d_dp_dFj += lsi->d_normal_dF[a][j] * nw[a];

    d_tmag_dFj = -tmaginv * dp * d_dp_dFj;

    for (a = 0, d_t_dFj[a] = 0; a < dim; a++)
      d_t_dFj[a] = (nf[a] - dp * nw[a]) * (-pow(tmaginv, 2)) * d_tmag_dFj +
                   (lsi->d_normal_dF[a][j] - d_dp_dFj * nw[a]) * tmaginv;

    switch (bc_type) {
    case SHARP_WETLIN_VELOCITY_BC:
      d_wet_speed_dFj = d_dp_dFj / cT;
      break;
    case SHARP_BLAKE_VELOCITY_BC:
      d_wet_speed_dFj = v0 * cosh(g * (cos_ca - cos_ca_static)) * g * d_dp_dFj;
      break;
    case SHARP_HOFFMAN_VELOCITY_BC:
      d_wet_speed_dFj = g / liq_visc * (1. + 1.31 * pow(g_dca, 0.99)) /
                        (1. - 1.31 * 0.99 * A_dca / pow(g_dca, 0.01)) * pow(A_dca, 0.294) * (-4.) *
                        (-d_dp_dFj) / (0.706 * 2 * 5.16 * (3. - costheta) * (1. + costheta));
      break;
    case SHARP_COX_VELOCITY_BC:
      /* sensitivity wrt to q_outer has not been included	*/
      f_den = log(reciprocal_slip) + (q_inner / f_sca - q_outer / f_dca);
      d_wet_speed_dFj = g / liq_visc * (1. / (f_dca * (-sintheta)) * (-d_dp_dFj)) / f_den;
      break;
    case SHARP_SHIK_VELOCITY_BC:
      dveloc0_ddpj =
          (-d_dp_dFj) * (SQUARE(theta) + theta * sintheta * costheta - 2 * SQUARE(sintheta)) /
          (SQUARE(sintheta * costheta) - 2 * theta * sintheta * costheta + SQUARE(theta));
      drhs_num_ddpj = d_dp_dFj;
      drhs_den_ddpj = (v0 - 1) * dveloc0_ddpj - d_dp_dFj;
      drhs_ddpj = (rhs_den * drhs_num_ddpj - rhs_num * drhs_den_ddpj) / SQUARE(rhs_den);
      d_wet_speed_dFj = 0.5 * sqrt(g * v0) * (1 - 0.5 * rhs) * drhs_ddpj / pow(1 + rhs, 1.5);
      break;
    } /* end of switch	*/
    d_wet_speed_dFj *= factor;
    for (a = 0, d_wet_vector_dF[a][j] = 0.0; a < dim; a++)
      d_wet_vector_dF[a][j] = d_t_dFj[a] * wet_speed + t[a] * d_wet_speed_dFj;
  }

  var = MESH_DISPLACEMENT1;
  if (pd->v[pg->imtrx][var]) {
    for (b = 0; b < dim; b++) {
      var = MESH_DISPLACEMENT1 + b;
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        for (d_dp_dXj[b] = 0.0, a = 0; a < dim; a++) {
          d_dp_dXj[b] += lsi->d_normal_dmesh[a][b][j] * nw[a];
        }

        d_tmag_dXj[b] = -tmaginv * dp * d_dp_dXj[b];

        for (a = 0, d_t_dXj[a][b] = 0.0; a < dim; a++) {
          d_t_dXj[a][b] = (nf[a] - dp * nw[a]) * (-pow(tmaginv, 2)) * d_tmag_dXj[b] +
                          (lsi->d_normal_dmesh[a][b][j] - d_dp_dXj[b] * nw[a]) * tmaginv;
        }

        for (a = 0, d_wet_vector_dX[a][b][j] = 0.0; a < dim; a++) {
          switch (bc_type) {
          case SHARP_WETLIN_VELOCITY_BC:
            d_wet_speed_dXj = d_dp_dXj[b] / cT;
            break;
          case SHARP_BLAKE_VELOCITY_BC:
            d_wet_speed_dXj = v0 * cosh(g * (cos_ca - cos_ca_static)) * g * d_dp_dXj[b];
            break;
          case SHARP_HOFFMAN_VELOCITY_BC:
            d_wet_speed_dXj = g / liq_visc * (1. + 1.31 * pow(g_dca, 0.99)) /
                              (1. - 1.31 * 0.99 * A_dca / pow(g_dca, 0.01)) * pow(A_dca, 0.294) *
                              (4.) / (0.706 * 2 * 5.16 * (3. - costheta) * (1. + costheta));
            d_wet_speed_dXj *= -d_dp_dXj[b];
            d_wet_speed_dXj += wet_speed * (-1. / liq_visc) * d_mu->X[b][j];
            break;
          case SHARP_COX_VELOCITY_BC:
            /* sensitivity wrt to q_outer has not been included	*/
            f_den = log(reciprocal_slip) + (q_inner / f_sca - q_outer / f_dca);
            d_wet_speed_dXj = g / liq_visc * (-1. / (f_dca * (-sintheta))) / f_den;
            d_wet_speed_dXj *= -d_dp_dXj[b];
            d_wet_speed_dXj += wet_speed * (-1. / liq_visc) * d_mu->X[b][j];
            break;
          case SHARP_SHIK_VELOCITY_BC:
            dveloc0_ddpj =
                (-d_dp_dXj[b]) *
                (SQUARE(theta) + theta * sintheta * costheta - 2 * SQUARE(sintheta)) /
                (SQUARE(sintheta * costheta) - 2 * theta * sintheta * costheta + SQUARE(theta));
            drhs_num_ddpj = d_dp_dXj[b];
            drhs_den_ddpj = (v0 - 1) * dveloc0_ddpj - d_dp_dXj[b];
            drhs_ddpj = (rhs_den * drhs_num_ddpj - rhs_num * drhs_den_ddpj) / SQUARE(rhs_den);
            d_wet_speed_dXj = 0.5 * sqrt(g * v0) * (1 - 0.5 * rhs) * drhs_ddpj / pow(1 + rhs, 1.5);
            break;
          } /* end of switch	*/
          d_wet_speed_dXj *= factor;
          d_wet_vector_dX[a][b][j] = d_t_dXj[a][b] * wet_speed + t[a] * d_wet_speed_dXj;
        }
      }
    }
  }
  /* compute stress tensor and its derivatives */
  if (include_stress) {
    fluid_stress(Pi, d_Pi);
  }

  for (a = 0; a < dim; a++) {
    func[a] = betainv * wet_vector[a];
    if (include_stress) {
      for (b = 0; b < dim; b++) {
        func[a] += nw[b] * Pi[a][b];
      }
    }
  }

  if (af->Assemble_Jacobian) {
    var = ls->var;

    for (a = 0; a < dim; a++) {
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_func[a][var][j] += betainv * d_wet_vector_dF[a][j];
          if (include_stress) {
            for (q = 0; q < dim; q++) {
              d_func[a][var][j] += nw[q] * d_Pi->F[a][q][j];
            }
          }
        }
      }
    }
    var = MESH_DISPLACEMENT1;

    if (pd->v[pg->imtrx][var]) {
      for (b = 0; b < dim; b++) {
        var = MESH_DISPLACEMENT1 + b;

        for (a = 0; a < dim; a++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_func[a][var][j] += betainv * d_wet_vector_dX[a][b][j];
          }
        }
      }
    }
    /* velocity senstivities		*/
    var = VELOCITY1;

    if (pd->v[pg->imtrx][var]) {
      for (b = 0; b < dim; b++) {
        var = VELOCITY1 + b;

        for (a = 0; a < dim; a++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            switch (bc_type) {
            case SHARP_WETLIN_VELOCITY_BC:
            case SHARP_BLAKE_VELOCITY_BC:
              d_wet_vector_dpj = 0;
              break;
            case SHARP_HOFFMAN_VELOCITY_BC:
              d_wet_vector_dpj = (g_dca - g_sca) *
                                 (liq_visc * cT * mp->d_surface_tension[var] * bf[var]->phi[j] -
                                  g * d_mu->v[b][j]) /
                                 SQUARE(liq_visc);
              break;
            case SHARP_COX_VELOCITY_BC:
              f_den = log(reciprocal_slip) + (q_inner / f_sca - q_outer / f_dca);
              d_wet_vector_dpj = g_integral / f_den *
                                 (liq_visc * cT * mp->d_surface_tension[var] * bf[var]->phi[j] -
                                  g * d_mu->v[b][j]) /
                                 SQUARE(liq_visc);
              break;
            case SHARP_SHIK_VELOCITY_BC:
              d_wet_vector_dpj =
                  0.5 * wet_speed / g * cT * mp->d_surface_tension[var] * bf[var]->phi[j];
              break;
            }
            d_func[a][var][j] += betainv * t[a] * d_wet_vector_dpj;
          }
        }
      }
    }
    /* Temperature senstivities		*/
    var = TEMPERATURE;

    if (pd->v[pg->imtrx][var]) {
      for (a = 0; a < dim; a++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          switch (bc_type) {
          case SHARP_WETLIN_VELOCITY_BC:
            d_wet_vector_dpj = 0;
            break;
          case SHARP_BLAKE_VELOCITY_BC:
            d_wet_vector_dpj = v0 * cosh((costhetaeq - costheta) * g) * (costhetaeq - costheta) *
                               cT * mp->d_surface_tension[var] * bf[var]->phi[j];
            break;
          case SHARP_HOFFMAN_VELOCITY_BC:
            d_wet_vector_dpj =
                (g_dca - g_sca) *
                (liq_visc * cT * mp->d_surface_tension[var] * bf[var]->phi[j] - g * d_mu->T[j]) /
                SQUARE(liq_visc);
            break;
          case SHARP_COX_VELOCITY_BC:
            f_den = log(reciprocal_slip) + (q_inner / f_sca - q_outer / f_dca);
            d_wet_vector_dpj =
                g_integral / f_den *
                (liq_visc * cT * mp->d_surface_tension[var] * bf[var]->phi[j] - g * d_mu->T[j]) /
                SQUARE(liq_visc);
            break;
          case SHARP_SHIK_VELOCITY_BC:
            d_wet_vector_dpj =
                0.5 * wet_speed / g * cT * mp->d_surface_tension[var] * bf[var]->phi[j];
            break;
          }
          d_func[a][var][j] += betainv * t[a] * d_wet_vector_dpj;
        }
      }
    }
    /*** concentration sensitivities	**/

    var = MASS_FRACTION;
    if (pd->v[pg->imtrx][var]) {
      for (a = 0; a < dim; a++) {
        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            switch (bc_type) {
            case SHARP_WETLIN_VELOCITY_BC:
              d_wet_vector_dpj = 0;
              break;
            case SHARP_BLAKE_VELOCITY_BC:
              d_wet_vector_dpj = v0 * cosh((costhetaeq - costheta) * g) * (costhetaeq - costheta) *
                                 cT * mp->d_surface_tension[MAX_VARIABLE_TYPES + w] *
                                 bf[var]->phi[j];
              break;
            case SHARP_HOFFMAN_VELOCITY_BC:
              d_wet_vector_dpj =
                  (g_dca - g_sca) *
                  (liq_visc * cT * mp->d_surface_tension[MAX_VARIABLE_TYPES + w] * bf[var]->phi[j] -
                   g * d_mu->C[w][j]) /
                  SQUARE(liq_visc);
              break;
            case SHARP_COX_VELOCITY_BC:
              f_den = log(reciprocal_slip) + (q_inner / f_sca - q_outer / f_dca);
              d_wet_vector_dpj =
                  g_integral / f_den *
                  (liq_visc * cT * mp->d_surface_tension[MAX_VARIABLE_TYPES + w] * bf[var]->phi[j] -
                   g * d_mu->C[w][j]) /
                  SQUARE(liq_visc);
              break;
            case SHARP_SHIK_VELOCITY_BC:
              d_wet_vector_dpj = 0.5 * wet_speed / g * cT *
                                 mp->d_surface_tension[MAX_VARIABLE_TYPES + w] * bf[var]->phi[j];
              break;
            }
            d_func[a][MAX_VARIABLE_TYPES + w][j] += betainv * t[a] * d_wet_vector_dpj;
          }
        }
      }
    } /*  end of concentration sensitivities	*/
    if (include_stress) {
      if (pd->v[pg->imtrx][VELOCITY1]) {
        for (a = 0; a < dim; a++) {
          for (q = 0; q < dim; q++) {
            for (b = 0; b < VIM; b++) {
              var = VELOCITY1 + b;
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                d_func[a][var][j] += nw[q] * d_Pi->v[a][q][b][j];
              }
            }
          }
        }
      }
      var = PRESSURE;
      if (pd->v[pg->imtrx][var]) {
        for (a = 0; a < dim; a++) {
          for (q = 0; q < dim; q++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_func[a][var][j] += fv->snormal[q] * d_Pi->P[a][q][j];
            }
          }
        }
      }
    }

  } /* if Jacobian	*/

  return;
} /* end of apply_sharp_wetting_velocity */

void apply_blake_wetting_velocity(double func[MAX_PDIM],
                                  double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                                  const int bc_type,
                                  int id_side,
                                  int element_type,
                                  double depth,
                                  double wetting_length,
                                  double theta_s,
                                  double v0,
                                  double surf_tens) {
  int a, b, j, dim = pd->Num_Dim, var, w;
  double betainv = 1.e12;
  double g = surf_tens;
  double wet_speed = 0.0, wet_vector[MAX_PDIM], d_wet_vector_dF[MAX_PDIM][MDE];
  double d_wet_vector_dX[MAX_PDIM][MAX_PDIM][MDE];
  double d_wet_vector_dpj = 0.0;
  double nw[MAX_PDIM], nf[MAX_PDIM], t[MAX_PDIM];
  double delta = 0, d_delta_dFj = 0, d_delta_dXj[DIM];
  /*  Hoffman correlation variables       */
  double ca_no, g_sca = 0.0, g_dca = 0.0, A_sca = 0.0, A_dca = 0.0;
  int iter, iter_max = 20;
  double eps_tol = 1.0e-12;
  double liq_visc = 0.0, gamma[DIM][DIM], visc_factor = 0.0;
  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct; /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;
  /*  Cox wetting parameters      */
  double theta = 0.0, thetaeq = 0.0, sintheta = 0.0, sinthetaeq = 0.0;
  double costhetaeq, costheta;
  double lambda = 0;
  double th, sinth, costh;
  double f_sca = 0.0, f_dca = 0.0, f_num, f_den, g_integral = 0.0;
  double reciprocal_slip = 0.0;
  const double q_inner = 0;
  const double q_outer = 0;
  /* use 10 point Gaussian quadrature for now     */
  const int num_gauss_pts = 10;
  const double gpt[10] = {-0.973906528517172, -0.865063366688985, -0.679409568299024,
                          -0.433395394129247, -0.148874338981631, 0.148874338981631,
                          0.433395394129247,  0.679409568299024,  0.865063366688985,
                          0.973906528517172};
  const double wt[10] = {0.066671344308688, 0.149451349150581, 0.219086362515982, 0.269266719309996,
                         0.295524224714753, 0.295524224714753, 0.269266719309996, 0.219086362515982,
                         0.149451349150581, 0.066671344308688};
  /* Shikhmurzaev wetting parameters	*/
  double rhs = 0.0, rhs_den = 0.0, rhs_num = 0.0, veloc0;
  double drhs_ddpj, drhs_den_ddpj, drhs_num_ddpj, dveloc0_ddpj;
  const double psg_surf = 0;

  if (ls == NULL) {
    GOMA_EH(GOMA_ERROR, "Boundary condition WETTING_SPEED bcs requires active level set dofs.\n");
  }
  /*  Should change surface_tension to come from mp->  */
  g = surf_tens * mp->surface_tension;

  if (wetting_length <= 0.0)
    wetting_length = ls->Length_Scale;

  load_lsi(wetting_length);

  for (a = 0; a < MAX_PDIM; a++)
    wet_vector[a] = 0.0;
  memset(d_wet_vector_dF, 0, sizeof(double) * MAX_PDIM * MDE);
  memset(d_wet_vector_dX, 0, sizeof(double) * MAX_PDIM * MAX_PDIM * MDE);

  betainv = 1.0 / depth;

  if (lsi->near) {
    double d_t_dFj[MAX_PDIM];
    double dp = 0.0, tmag = 0.0, tmaginv = 0.0;
    double d_dp_dFj, d_tmag_dFj;
    double d_wet_speed_dFj = 0.0, d_wet_speed_dXj = 0.0;
    double cos_ca, cos_ca_static;
    double d_dp_dXj[DIM], d_t_dXj[DIM][DIM], d_tmag_dXj[DIM];

    memset(t, 0, sizeof(double) * MAX_PDIM);
    memset(d_t_dFj, 0, sizeof(double) * MAX_PDIM);
    memset(d_t_dXj, 0, sizeof(double) * DIM * DIM);

    cos_ca_static = cos(M_PIE * (180.0 - theta_s) / 180.0);
    costhetaeq = -cos_ca_static;

    if (af->Assemble_Jacobian)
      load_lsi_derivs();
    if (bc_type == WETTING_SPEED_HOFFMAN_BC || bc_type == WETTING_SPEED_COX_BC) {
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
        }
      }
      liq_visc = viscosity(gn, gamma, d_mu);
      switch (mp->mp2nd->viscositymask[1] - mp->mp2nd->viscositymask[0]) {
      case 1:
        liq_visc = (liq_visc - mp->mp2nd->viscosity * lsi->H) / (1. - lsi->H);
        visc_factor = 1. / (1. - lsi->H);
        break;
      case -1:
        liq_visc = (liq_visc - mp->mp2nd->viscosity * (1. - lsi->H)) / lsi->H;
        visc_factor = 1. / lsi->H;
        break;
      }
    }

    for (a = 0; a < dim; a++) {
      nw[a] = fv->snormal[a];
      nf[a] = lsi->normal[a];

      dp += nw[a] * nf[a];
    }

    tmag = sqrt(1.0 - pow(dp, 2.0));
    if (tmag <= DBL_SMALL) {
      printf("Underflow tmag detected.\n");
      tmaginv = 0.0;
    } else
      tmaginv = 1.0 / tmag;

    for (a = 0; a < dim; a++) {
      t[a] = (nf[a] - dp * nw[a]) * tmaginv;
    }

    cos_ca = dp;
    costheta = -cos_ca;

    switch (bc_type) {
    case WETTING_SPEED_BLAKE_BC:
      wet_speed = v0 * sinh(g * (cos_ca_static - cos_ca));
      break;
    case WETTING_SPEED_HOFFMAN_BC:
      ca_no = 1.0E+06;
      iter = 0;
      eps = 10. * eps_tol;
      while (iter <= iter_max && fabs(eps) > eps_tol) {
        g_sca = log((3. - costhetaeq) / (1. + costhetaeq)) / (2. * 5.16);
        A_sca = pow(g_sca, 1. / 0.706);
        eps = -(ca_no - 1.31 * pow(ca_no, 0.99) * A_sca - A_sca) /
              (1. - 1.31 * 0.99 * A_sca / pow(ca_no, 0.01));
        ca_no += eps;
        iter++;
      }
      if (fabs(eps) > eps_tol)
        GOMA_EH(GOMA_ERROR, "Hoffman iteration not converged");
      g_sca = ca_no;
      ca_no = 1.0E+06;
      iter = 0;
      eps = 10. * eps_tol;
      while (iter <= iter_max && fabs(eps) > eps_tol) {
        g_dca = log((3. - costheta) / (1. + costheta)) / (2. * 5.16);
        A_dca = pow(g_dca, 1. / 0.706);
        eps = -(ca_no - 1.31 * pow(ca_no, 0.99) * A_dca - A_dca) /
              (1. - 1.31 * 0.99 * A_dca / pow(ca_no, 0.01));
        ca_no += eps;
        iter++;
      }
      if (fabs(eps) > eps_tol)
        GOMA_EH(GOMA_ERROR, "Hoffman iteration not converged");
      g_dca = ca_no;
      ca_no = g_dca - g_sca;
      wet_speed = ca_no * g / liq_visc;
      break;
    case WETTING_SPEED_COX_BC:
      thetaeq = M_PIE * theta_s / 180.0;
      theta = acos(costheta);
      sintheta = sin(theta);
      sinthetaeq = sin(thetaeq);
      reciprocal_slip = 1. / v0;

      /*  Cox analysis integral       */
      g_integral = 0;
      for (j = 0; j < num_gauss_pts; j++) {
        th = thetaeq + (theta - thetaeq) * (gpt[j] + 1.) / 2.;
        sinth = sin(th);
        costh = cos(th);
        f_num = 2. * sinth *
                (SQUARE(lambda) * (SQUARE(th) - SQUARE(sinth)) +
                 2. * lambda * (th * (M_PIE - th) + SQUARE(sinth)) +
                 (SQUARE(M_PIE - th) - SQUARE(sinth)));
        f_den = lambda * (SQUARE(th) - SQUARE(sinth)) * (M_PIE - th + sinth * costh) +
                (SQUARE(M_PIE - th) - SQUARE(sinth)) * (th - sinth * costh);
        g_integral += wt[j] * f_den / f_num;
      }
      g_integral *= 0.5 * (theta - thetaeq);
      /*  inner and outer solution terms      */
      f_num = 2. * sinthetaeq *
              (SQUARE(lambda) * (SQUARE(thetaeq) - SQUARE(sinthetaeq)) +
               2. * lambda * (thetaeq * (M_PIE - thetaeq) + SQUARE(sinthetaeq)) +
               (SQUARE(M_PIE - thetaeq) - SQUARE(sinthetaeq)));
      f_den = lambda * (SQUARE(thetaeq) - SQUARE(sinthetaeq)) *
                  (M_PIE - thetaeq + sinthetaeq * costhetaeq) +
              (SQUARE(M_PIE - thetaeq) - SQUARE(sinthetaeq)) * (thetaeq - sinthetaeq * costhetaeq);
      f_sca = f_num / f_den;
      f_num = 2. * sintheta *
              (SQUARE(lambda) * (SQUARE(theta) - SQUARE(sintheta)) +
               2. * lambda * (theta * (M_PIE - theta) + SQUARE(sintheta)) +
               (SQUARE(M_PIE - theta) - SQUARE(sintheta)));
      f_den = lambda * (SQUARE(theta) - SQUARE(sintheta)) * (M_PIE - theta + sintheta * costheta) +
              (SQUARE(M_PIE - theta) - SQUARE(sintheta)) * (theta - sintheta * costheta);
      f_dca = f_num / f_den;
      /* solve for wetting speed      */
      ca_no = g_integral / (log(reciprocal_slip) + (q_inner / f_sca - q_outer / f_dca));
      wet_speed = ca_no * g / liq_visc;
      break;
    case WETTING_SPEED_SHIK_BC:
      theta = acos(costheta);
      sintheta = sin(theta);
      veloc0 = (sintheta - theta * costheta) / (sintheta * costheta - theta);
      rhs_num = costhetaeq - costheta;
      rhs_den = psg_surf + v0 + (v0 - 1.) * veloc0 + costheta;
      rhs = rhs_num / rhs_den;
      wet_speed = sqrt(g * v0) * rhs / (2. * sqrt(1. + rhs));
      break;
    default:
      GOMA_EH(GOMA_ERROR, "bad DCA bc name\n");
      break;
    }

    delta = lsi->delta;

    for (a = 0; a < dim; a++) {
      wet_vector[a] = t[a] * wet_speed * delta;
    }

    betainv = 1.0 / depth;
    var = ls->var;

    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      for (d_dp_dFj = 0, a = 0; a < dim; a++)
        d_dp_dFj += lsi->d_normal_dF[a][j] * nw[a];

      switch (bc_type) {
      case WETTING_SPEED_BLAKE_BC:
        d_wet_speed_dFj = v0 * cosh(g * (cos_ca - cos_ca_static)) * g * d_dp_dFj;
        break;
      case WETTING_SPEED_HOFFMAN_BC:
        d_wet_speed_dFj = g / liq_visc * (1. + 1.31 * pow(g_dca, 0.99)) /
                          (1. - 1.31 * 0.99 * A_dca / pow(g_dca, 0.01)) * pow(A_dca, 0.294) *
                          (-4.) * (-d_dp_dFj) /
                          (0.706 * 2 * 5.16 * (3. - costheta) * (1. + costheta));
        break;
      case WETTING_SPEED_COX_BC:
        /* sensitivity wrt to q_outer has not been included	*/
        f_den = log(reciprocal_slip) + (q_inner / f_sca - q_outer / f_dca);
        d_wet_speed_dFj = g / liq_visc * (1. / (f_dca * (-sintheta)) * (-d_dp_dFj)) / f_den;
        break;
      case WETTING_SPEED_SHIK_BC:
        dveloc0_ddpj =
            (-d_dp_dFj) * (SQUARE(theta) + theta * sintheta * costheta - 2 * SQUARE(sintheta)) /
            (SQUARE(sintheta * costheta) - 2 * theta * sintheta * costheta + SQUARE(theta));
        drhs_num_ddpj = d_dp_dFj;
        drhs_den_ddpj = (v0 - 1) * dveloc0_ddpj - d_dp_dFj;
        drhs_ddpj = (rhs_den * drhs_num_ddpj - rhs_num * drhs_den_ddpj) / SQUARE(rhs_den);
        d_wet_speed_dFj = 0.5 * sqrt(g * v0) * (1 - 0.5 * rhs) * drhs_ddpj / pow(1 + rhs, 1.5);
        break;
      default:
        GOMA_EH(GOMA_ERROR, "bad DCA bc name\n");
        break;
      }

      d_tmag_dFj = -tmaginv * dp * d_dp_dFj;

      for (a = 0, d_t_dFj[a] = 0; a < dim; a++)
        d_t_dFj[a] = (nf[a] - dp * nw[a]) * (-pow(tmaginv, 2)) * d_tmag_dFj +
                     (lsi->d_normal_dF[a][j] - d_dp_dFj * nw[a]) * tmaginv;

      d_delta_dFj = lsi->d_delta_dF[j];

      for (a = 0, d_wet_vector_dF[a][j] = 0.0; a < dim; a++)
        d_wet_vector_dF[a][j] = d_t_dFj[a] * wet_speed * delta + t[a] * d_wet_speed_dFj * delta +
                                t[a] * wet_speed * d_delta_dFj;
    }

    var = MESH_DISPLACEMENT1;
    if (pd->v[pg->imtrx][var]) {
      for (b = 0; b < dim; b++) {
        var = MESH_DISPLACEMENT1 + b;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (d_dp_dXj[b] = 0.0, a = 0; a < dim; a++) {
            d_dp_dXj[b] += lsi->d_normal_dmesh[a][b][j] * nw[a];
          }

          d_tmag_dXj[b] = -tmaginv * dp * d_dp_dXj[b];

          for (a = 0, d_t_dXj[a][b] = 0.0; a < dim; a++) {
            d_t_dXj[a][b] = (nf[a] - dp * nw[a]) * (-pow(tmaginv, 2)) * d_tmag_dXj[b] +
                            (lsi->d_normal_dmesh[a][b][j] - d_dp_dXj[b] * nw[a]) * tmaginv;
          }

          d_delta_dXj[b] = lsi->d_delta_dmesh[b][j];

          for (a = 0, d_wet_vector_dX[a][b][j] = 0.0; a < dim; a++) {
            switch (bc_type) {
            case WETTING_SPEED_BLAKE_BC:
              d_wet_speed_dXj = v0 * cosh(g * (cos_ca - cos_ca_static)) * g * d_dp_dXj[b];
              break;
            case WETTING_SPEED_HOFFMAN_BC:
              d_wet_speed_dXj = g / liq_visc * (1. + 1.31 * pow(g_dca, 0.99)) /
                                (1. - 1.31 * 0.99 * A_dca / pow(g_dca, 0.01)) * pow(A_dca, 0.294) *
                                (4.) / (0.706 * 2 * 5.16 * (3. - costheta) * (1. + costheta));
              d_wet_speed_dXj *= -d_dp_dXj[b];
              d_wet_speed_dXj += wet_speed * (-1. / liq_visc) * d_mu->X[b][j] * visc_factor;
              break;
            case WETTING_SPEED_COX_BC:
              /* sensitivity wrt to q_outer has not been included	*/
              f_den = log(reciprocal_slip) + (q_inner / f_sca - q_outer / f_dca);
              d_wet_speed_dXj = g / liq_visc * (-1. / (f_dca * (-sintheta))) / f_den;
              d_wet_speed_dXj *= -d_dp_dXj[b];
              d_wet_speed_dXj += wet_speed * (-1. / liq_visc) * d_mu->X[b][j] * visc_factor;
              break;
            case WETTING_SPEED_SHIK_BC:
              dveloc0_ddpj =
                  (-d_dp_dXj[b]) *
                  (SQUARE(theta) + theta * sintheta * costheta - 2 * SQUARE(sintheta)) /
                  (SQUARE(sintheta * costheta) - 2 * theta * sintheta * costheta + SQUARE(theta));
              drhs_num_ddpj = d_dp_dXj[b];
              drhs_den_ddpj = (v0 - 1) * dveloc0_ddpj - d_dp_dXj[b];
              drhs_ddpj = (rhs_den * drhs_num_ddpj - rhs_num * drhs_den_ddpj) / SQUARE(rhs_den);
              d_wet_speed_dXj =
                  0.5 * sqrt(g * v0) * (1 - 0.5 * rhs) * drhs_ddpj / pow(1 + rhs, 1.5);
              break;
            } /* end of switch	*/
            d_wet_vector_dX[a][b][j] = d_t_dXj[a][b] * wet_speed * delta +
                                       t[a] * d_wet_speed_dXj * delta +
                                       t[a] * wet_speed * d_delta_dXj[b];
          }
        }
      }
    }

    /*}*/

    for (a = 0; a < pd->Num_Dim; a++) {
      func[a] = betainv * wet_vector[a];
    }

    if (af->Assemble_Jacobian) {
      var = ls->var;

      for (a = 0; a < dim; a++) {
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_func[a][var][j] += (betainv) * (d_wet_vector_dF[a][j]);
          }
        }
      }
      var = MESH_DISPLACEMENT1;

      if (pd->v[pg->imtrx][var]) {
        for (b = 0; b < dim; b++) {
          var = MESH_DISPLACEMENT1 + b;

          for (a = 0; a < dim; a++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_func[a][var][j] += betainv * d_wet_vector_dX[a][b][j];
            }
          }
        }
      }
      /* velocity senstivities		*/
      var = VELOCITY1;

      if (pd->v[pg->imtrx][var]) {
        for (b = 0; b < dim; b++) {
          var = VELOCITY1 + b;

          for (a = 0; a < dim; a++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              switch (bc_type) {
              case WETTING_SPEED_BLAKE_BC:
                d_wet_vector_dpj = 0;
                break;
              case WETTING_SPEED_HOFFMAN_BC:
                d_wet_vector_dpj =
                    (g_dca - g_sca) *
                    (liq_visc * surf_tens * mp->d_surface_tension[var] * bf[var]->phi[j] -
                     g * d_mu->v[b][j] * visc_factor) /
                    SQUARE(liq_visc);
                break;
              case WETTING_SPEED_COX_BC:
                f_den = log(reciprocal_slip) + (q_inner / f_sca - q_outer / f_dca);
                d_wet_vector_dpj =
                    g_integral / f_den *
                    (liq_visc * surf_tens * mp->d_surface_tension[var] * bf[var]->phi[j] -
                     g * d_mu->v[b][j] * visc_factor) /
                    SQUARE(liq_visc);
                break;
              case WETTING_SPEED_SHIK_BC:
                d_wet_vector_dpj =
                    0.5 * wet_speed / g * surf_tens * mp->d_surface_tension[var] * bf[var]->phi[j];
                break;
              }
              d_func[a][var][j] += betainv * t[a] * d_wet_vector_dpj * delta;
            }
          }
        }
      }
      /* Temperature senstivities		*/
      var = TEMPERATURE;

      if (pd->v[pg->imtrx][var]) {
        for (a = 0; a < dim; a++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            switch (bc_type) {
            case WETTING_SPEED_BLAKE_BC:
              d_wet_vector_dpj = v0 * cosh((costhetaeq - costheta) * g) * (costhetaeq - costheta) *
                                 surf_tens * mp->d_surface_tension[var] * bf[var]->phi[j];
              break;
            case WETTING_SPEED_HOFFMAN_BC:
              d_wet_vector_dpj =
                  (g_dca - g_sca) *
                  (liq_visc * surf_tens * mp->d_surface_tension[var] * bf[var]->phi[j] -
                   g * d_mu->T[j] * visc_factor) /
                  SQUARE(liq_visc);
              break;
            case WETTING_SPEED_COX_BC:
              f_den = log(reciprocal_slip) + (q_inner / f_sca - q_outer / f_dca);
              d_wet_vector_dpj =
                  g_integral / f_den *
                  (liq_visc * surf_tens * mp->d_surface_tension[var] * bf[var]->phi[j] -
                   g * d_mu->T[j] * visc_factor) /
                  SQUARE(liq_visc);
              break;
            case WETTING_SPEED_SHIK_BC:
              d_wet_vector_dpj =
                  0.5 * wet_speed / g * surf_tens * mp->d_surface_tension[var] * bf[var]->phi[j];
              break;
            }
            d_func[a][var][j] += betainv * t[a] * d_wet_vector_dpj * delta;
          }
        }
      }
      /*** concentration sensitivities	**/

      var = MASS_FRACTION;
      if (pd->v[pg->imtrx][var]) {
        for (a = 0; a < dim; a++) {
          for (w = 0; w < pd->Num_Species_Eqn; w++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              switch (bc_type) {
              case WETTING_SPEED_BLAKE_BC:
                d_wet_vector_dpj = v0 * cosh((costhetaeq - costheta) * g) *
                                   (costhetaeq - costheta) * surf_tens *
                                   mp->d_surface_tension[MAX_VARIABLE_TYPES + w] * bf[var]->phi[j];
                break;
              case WETTING_SPEED_HOFFMAN_BC:
                d_wet_vector_dpj =
                    (g_dca - g_sca) *
                    (liq_visc * surf_tens * mp->d_surface_tension[MAX_VARIABLE_TYPES + w] *
                         bf[var]->phi[j] -
                     g * d_mu->C[w][j] * visc_factor) /
                    SQUARE(liq_visc);
                break;
              case WETTING_SPEED_COX_BC:
                f_den = log(reciprocal_slip) + (q_inner / f_sca - q_outer / f_dca);
                d_wet_vector_dpj =
                    g_integral / f_den *
                    (liq_visc * surf_tens * mp->d_surface_tension[MAX_VARIABLE_TYPES + w] *
                         bf[var]->phi[j] -
                     g * d_mu->C[w][j] * visc_factor) /
                    SQUARE(liq_visc);
                break;
              case WETTING_SPEED_SHIK_BC:
                d_wet_vector_dpj = 0.5 * wet_speed / g * surf_tens *
                                   mp->d_surface_tension[MAX_VARIABLE_TYPES + w] * bf[var]->phi[j];
                break;
              }
              d_func[a][MAX_VARIABLE_TYPES + w][j] += betainv * t[a] * d_wet_vector_dpj * delta;
            }
          }
        }
      } /*  end of concentration sensitivities	*/
    }

  } /*  lsi->near  */
  return;
}

void apply_blake_wetting_velocity_sic(double func[MAX_PDIM],
                                      double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                                      const double dt,
                                      const double tt,
                                      int bc_type,
                                      const int id_side,
                                      const int element_type,
                                      double wetting_length,
                                      const double theta_s,
                                      const double v0,
                                      const double surf_tens,
                                      const double tau,
                                      double vs[MAX_PDIM],
                                      const double theta_max_degrees,
                                      const int rolling_substrate,
                                      const double beta_outside,
                                      const double gas_phase_factor,
                                      const double contact_fraction) {
  double wet_speed = 0.0, g;
  double cos_ca, cos_ca_static;
  double triangle;
  int triangle_form = TRUE;
  int var, a, b, j, w, dim = pd->Num_Dim;
  double phi_j;
  double t[DIM], tmaginv, dp, sin2;
  double betainv; /*penalty = pow(BIG_PENALTY*LITTLE_PENALTY, 0.5 );  */
  double nw[MAX_PDIM], nf[MAX_PDIM];
  double d_wet_vector_dF[MAX_PDIM][MDE];
  double d_wet_vector_dX[MAX_PDIM][MAX_PDIM][MDE];
  double d_wet_vector_dpj = 0.0, d_triangle_dF;
  double d_t_dFj[MAX_PDIM];
  double d_dp_dFj, d_tmag_dFj;
  double d_wet_speed_dFj = 0.0, d_wet_speed_dXj = 0.0, d_beta_dF[MDE];
  double d_dp_dXj[DIM], d_t_dXj[DIM][DIM], d_tmag_dXj[DIM];
  double beta, beta_inside = 1. / sqrt(BIG_PENALTY * LITTLE_PENALTY);
  /*  Hoffman correlation variables       */
  double ca_no = 0.0, g_sca = 0.0, g_dca = 0.0, A_sca = 0.0, A_dca = 0.0;
  int iter, iter_max = 20;
  double eps_tol = 1.0e-12;
  double liq_visc = 0.0, gamma[DIM][DIM], visc_factor = 0.0, lsiF_org;
  int elem_sign_org;
  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct; /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;
  /*  Cox wetting parameters      */
  double theta = 0.0, thetaeq, sintheta = 0.0, sinthetaeq;
  double costhetaeq, costheta;
  double lambda = 0;
  double th, sinth, costh;
  double f_sca = 0.0, f_dca = 0.0, f_num = 0.0, f_den = 0.0, g_integral = 0.0;
  double reciprocal_slip = 0.0;
  const double q_inner = 0.0;
  const double q_outer = 0.0;
  /* use 10 point Gaussian quadrature for now     */
  const int num_gauss_pts = 10;
  const double gpt[10] = {-0.973906528517172, -0.865063366688985, -0.679409568299024,
                          -0.433395394129247, -0.148874338981631, 0.148874338981631,
                          0.433395394129247,  0.679409568299024,  0.865063366688985,
                          0.973906528517172};
  const double wt[10] = {0.066671344308688, 0.149451349150581, 0.219086362515982, 0.269266719309996,
                         0.295524224714753, 0.295524224714753, 0.269266719309996, 0.219086362515982,
                         0.149451349150581, 0.066671344308688};
  /* Shikhmurzaev wetting parameters	*/
  double rhs = 0.0, rhs_den = 0.0, rhs_num = 0.0, veloc0, veloc0max;
  double drhs_ddpj, drhs_den_ddpj, drhs_num_ddpj, dveloc0_ddpj;
  double theta_max = M_PIE * theta_max_degrees / 180., costhetamax, sinthetamax;
  const double shik_max_factor = 1.01;

  if (ls == NULL) {
    GOMA_EH(GOMA_ERROR, "Boundary condition BLAKE_DIRICHLET requires active level set dofs.\n");
  }

  if (rolling_substrate == TRUE) {
    double rpm, X_0, X_1;
    /* For roll conditions compute variable wall velocity */
    rpm = vs[0]; /* revs / time */
    X_0 = vs[1]; /* Roll center , x */
    X_1 = vs[2]; /* Roll center, y */

    vs[0] = 2.0 * M_PIE * rpm * (fv->x[1] - X_1);
    vs[1] = -2.0 * M_PIE * rpm * (fv->x[0] - X_0);
    vs[2] = 0;

    /* Also map wetting line models from BC type */
    switch (bc_type) {
    case BLAKE_DIRICH_ROLL_BC:
      bc_type = BLAKE_DIRICHLET_BC;
      break;
    case HOFFMAN_DIRICH_ROLL_BC:
      bc_type = HOFFMAN_DIRICHLET_BC;
      break;
    case COX_DIRICH_ROLL_BC:
      bc_type = COX_DIRICHLET_BC;
      break;
    case SHIK_DIRICH_ROLL_BC:
      bc_type = SHIK_DIRICHLET_BC;
      break;
    default:
      GOMA_EH(GOMA_ERROR, "Tried to imposing rolling substrate on non-rolling BC.\n");
      break;
    }
  }

  if (wetting_length <= 0.0)
    wetting_length = ls->Length_Scale;
  if (wetting_length <= 0.0) {
    GOMA_EH(
        GOMA_ERROR,
        "Boundary condition BLAKE_DIRICHLET positive level set length scale or wetting length.\n");
  }

  load_lsi(wetting_length);
  if (af->Assemble_Jacobian) {
    load_lsi_derivs();
    memset(d_beta_dF, 0, MDE * sizeof(double));
  }

  beta = slip_coefficient(beta_inside, beta_outside, fv->F, wetting_length, d_beta_dF,
                          gas_phase_factor, contact_fraction);

  betainv = 1.0 / beta;

  if (lsi->near) {
    cos_ca_static = cos(M_PIE * theta_s / 180.0);
    costhetaeq = cos_ca_static;

    /*  Should change surface_tension to come from mp->  */
    g = surf_tens * mp->surface_tension;

    memset(d_wet_vector_dF, 0, sizeof(double) * MAX_PDIM * MDE);
    memset(d_wet_vector_dX, 0, sizeof(double) * MAX_PDIM * MAX_PDIM * MDE);
    memset(d_t_dFj, 0, sizeof(double) * MAX_PDIM);
    memset(d_t_dXj, 0, sizeof(double) * DIM * DIM);

    dp = 0.;
    for (a = 0; a < pd->Num_Dim; a++) {
      nw[a] = fv->snormal[a];
      nf[a] = lsi->normal[a];
      dp += nw[a] * nf[a];
    }

    sin2 = 1.0 - pow(dp, 2.0);
    if (sin2 <= 1. / BIG_PENALTY) {
      tmaginv = sqrt(BIG_PENALTY);
    } else {
      tmaginv = 1.0 / sqrt(sin2);
    }
    for (a = 0; a < pd->Num_Dim; a++) {
      t[a] = (lsi->normal[a] - dp * fv->snormal[a]) * tmaginv;
    }

    if (bc_type == HOFFMAN_DIRICHLET_BC || bc_type == COX_DIRICHLET_BC) {
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
        }
      }

      if (ls->SubElemIntegration) {
        elem_sign_org = ls->Elem_Sign;
        switch (mp->mp2nd->viscositymask[1] - mp->mp2nd->viscositymask[0]) {
        case 1:
          ls->Elem_Sign = -1;
          break;
        case -1:
          ls->Elem_Sign = 1;
          break;
        }
        liq_visc = viscosity(gn, gamma, d_mu);
        ls->Elem_Sign = elem_sign_org;
        visc_factor = 1.;
      } else {
        lsiF_org = fv->F;
        fv->F = 0.;
        liq_visc = viscosity(gn, gamma, d_mu);
        switch (mp->mp2nd->viscositymask[1] - mp->mp2nd->viscositymask[0]) {
        case 1:
          liq_visc = (liq_visc - mp->mp2nd->viscosity * lsi->H) / (1. - lsi->H);
          visc_factor = 1. / (1. - lsi->H);
          break;
        case -1:
          liq_visc = (liq_visc - mp->mp2nd->viscosity * (1. - lsi->H)) / lsi->H;
          visc_factor = 1. / lsi->H;
          break;
        }
        fv->F = lsiF_org;
      }
    }

    cos_ca = -dp;
    costheta = cos_ca;
    switch (bc_type) {
    case BLAKE_DIRICHLET_BC:
      wet_speed = v0 * sinh(g * (cos_ca_static - cos_ca));
      break;
    case HOFFMAN_DIRICHLET_BC:
      ca_no = 1.0E+06;
      iter = 0;
      eps = 10. * eps_tol;
      while (iter <= iter_max && fabs(eps) > eps_tol) {
        g_sca = log((3. - costhetaeq) / (1. + costhetaeq)) / (2. * 5.16);
        A_sca = pow(g_sca, 1. / 0.706);
        eps = -(ca_no - 1.31 * pow(ca_no, 0.99) * A_sca - A_sca) /
              (1. - 1.31 * 0.99 * A_sca / pow(ca_no, 0.01));
        ca_no += eps;
        iter++;
      }
      if (fabs(eps) > eps_tol)
        GOMA_EH(GOMA_ERROR, "Hoffman iteration not converged");
      g_sca = ca_no;
      ca_no = 1.0E+06;
      iter = 0;
      eps = 10. * eps_tol;
      while (iter <= iter_max && fabs(eps) > eps_tol) {
        g_dca = log((3. - costheta) / (1. + costheta)) / (2. * 5.16);
        A_dca = pow(g_dca, 1. / 0.706);
        eps = -(ca_no - 1.31 * pow(ca_no, 0.99) * A_dca - A_dca) /
              (1. - 1.31 * 0.99 * A_dca / pow(ca_no, 0.01));
        ca_no += eps;
        iter++;
      }
      if (fabs(eps) > eps_tol)
        GOMA_WH(-1, "Hoffman iteration not converged");
      g_dca = ca_no;
      ca_no = g_dca - g_sca;
      wet_speed = ca_no * g / liq_visc;
      break;
    case COX_DIRICHLET_BC:
      thetaeq = M_PIE * theta_s / 180.0;
      theta = acos(costheta);
      sintheta = sin(theta);
      sinthetaeq = sin(thetaeq);
      reciprocal_slip = 1. / v0;

      /*  Cox analysis integral       */
      g_integral = 0;
      for (j = 0; j < num_gauss_pts; j++) {
        th = thetaeq + (theta - thetaeq) * (gpt[j] + 1.) / 2.;
        sinth = sin(th);
        costh = cos(th);
        f_num = 2. * sinth *
                (SQUARE(lambda) * (SQUARE(th) - SQUARE(sinth)) +
                 2. * lambda * (th * (M_PIE - th) + SQUARE(sinth)) +
                 (SQUARE(M_PIE - th) - SQUARE(sinth)));
        f_den = lambda * (SQUARE(th) - SQUARE(sinth)) * (M_PIE - th + sinth * costh) +
                (SQUARE(M_PIE - th) - SQUARE(sinth)) * (th - sinth * costh);
        g_integral += wt[j] * f_den / f_num;
      }
      g_integral *= 0.5 * (theta - thetaeq);
      /*  inner and outer solution terms      */
      f_num = 2. * sinthetaeq *
              (SQUARE(lambda) * (SQUARE(thetaeq) - SQUARE(sinthetaeq)) +
               2. * lambda * (thetaeq * (M_PIE - thetaeq) + SQUARE(sinthetaeq)) +
               (SQUARE(M_PIE - thetaeq) - SQUARE(sinthetaeq)));
      f_den = lambda * (SQUARE(thetaeq) - SQUARE(sinthetaeq)) *
                  (M_PIE - thetaeq + sinthetaeq * costhetaeq) +
              (SQUARE(M_PIE - thetaeq) - SQUARE(sinthetaeq)) * (thetaeq - sinthetaeq * costhetaeq);
      f_sca = f_num / f_den;
      f_num = 2. * sintheta *
              (SQUARE(lambda) * (SQUARE(theta) - SQUARE(sintheta)) +
               2. * lambda * (theta * (M_PIE - theta) + SQUARE(sintheta)) +
               (SQUARE(M_PIE - theta) - SQUARE(sintheta)));
      f_den = lambda * (SQUARE(theta) - SQUARE(sintheta)) * (M_PIE - theta + sintheta * costheta) +
              (SQUARE(M_PIE - theta) - SQUARE(sintheta)) * (theta - sintheta * costheta);
      f_dca = f_num / f_den;
      /* solve for wetting speed      */
      ca_no = g_integral / (log(reciprocal_slip) + (q_inner / f_sca - q_outer / f_dca));
      wet_speed = ca_no * g / liq_visc;
      break;
    case SHIK_DIRICHLET_BC:
      theta = acos(costheta);
      sintheta = sin(theta);
      if (theta > theta_max) {
        theta_max = MIN(shik_max_factor * theta, M_PIE);
      }
      costhetamax = cos(theta_max);
      sinthetamax = sin(theta_max);
      veloc0 = (sintheta - theta * costheta) / (sintheta * costheta - theta);
      veloc0max = (sinthetamax - theta_max * costhetamax) / (sinthetamax * costhetamax - theta_max);
      rhs_den = (v0 - 1.) * (veloc0 - veloc0max) + costheta - costhetamax;
      rhs_num = costhetaeq - costheta;
      rhs = rhs_num / rhs_den;
      wet_speed = sqrt(g * v0) * rhs / (2. * sqrt(1. + rhs));
      break;
    default:
      GOMA_EH(GOMA_ERROR, "bad DCA bc name\n");
      break;
    }

    if (triangle_form) {
      triangle = 0.;
      d_triangle_dF = 0.;
      if (fv->F < -wetting_length) {
        triangle = 0.;
        d_triangle_dF = 0.;
      } else if (fv->F < 0.) {
        triangle = 1. + fv->F / wetting_length;
        d_triangle_dF = 1. / wetting_length;
      } else if (fv->F < wetting_length) {
        triangle = 1. - fv->F / wetting_length;
        d_triangle_dF = -1. / wetting_length;
      } else {
        triangle = 0.;
        d_triangle_dF = 0.;
      }
    } else {
      if (wetting_length < 1.e-12)
        wetting_length = 1.;
      triangle = 0.5 * lsi->delta * wetting_length;
    }

#if 0
      printf("Ca angle vnew %g %g %g %g\n",ca_no, acos(costheta)*180/M_PIE, wet_speed, triangle);
#endif

    var = ls->var;

    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      for (d_dp_dFj = 0, a = 0; a < dim; a++)
        d_dp_dFj += lsi->d_normal_dF[a][j] * nw[a];

      switch (bc_type) {
      case BLAKE_DIRICHLET_BC:
        d_wet_speed_dFj = v0 * cosh(g * (cos_ca - cos_ca_static)) * g * d_dp_dFj;
        break;
      case HOFFMAN_DIRICHLET_BC:
        d_wet_speed_dFj = g / liq_visc * (1. + 1.31 * pow(g_dca, 0.99)) /
                          (1. - 1.31 * 0.99 * A_dca / pow(g_dca, 0.01)) * pow(A_dca, 0.294) *
                          (-4.) * (-d_dp_dFj) /
                          (0.706 * 2 * 5.16 * (3. - costheta) * (1. + costheta));
        break;
      case COX_DIRICHLET_BC:
        /* sensitivity wrt to q_outer has not been included	*/
        f_den = log(reciprocal_slip) + (q_inner / f_sca - q_outer / f_dca);
        d_wet_speed_dFj = g / liq_visc * (1. / (f_dca * (-sintheta)) * (-d_dp_dFj)) / f_den;
        break;
      case SHIK_DIRICHLET_BC:
        dveloc0_ddpj =
            (-d_dp_dFj) * (SQUARE(theta) + theta * sintheta * costheta - 2 * SQUARE(sintheta)) /
            (SQUARE(sintheta * costheta) - 2 * theta * sintheta * costheta + SQUARE(theta));
        drhs_num_ddpj = d_dp_dFj;
        drhs_den_ddpj = (v0 - 1) * dveloc0_ddpj - d_dp_dFj;
        if (theta_max > theta_max_degrees * M_PIE / 180.)
          drhs_den_ddpj *= (1. - 1. / shik_max_factor);
        drhs_ddpj = (rhs_den * drhs_num_ddpj - rhs_num * drhs_den_ddpj) / SQUARE(rhs_den);
        d_wet_speed_dFj = 0.5 * sqrt(g * v0) * (1 - 0.5 * rhs) * drhs_ddpj / pow(1 + rhs, 1.5);
        break;
      default:
        GOMA_EH(GOMA_ERROR, "bad DCA bc name\n");
        break;
      }

      d_tmag_dFj = -tmaginv * dp * d_dp_dFj;

      for (a = 0, d_t_dFj[a] = 0; a < dim; a++)
        d_t_dFj[a] = (nf[a] - dp * nw[a]) * (-pow(tmaginv, 2)) * d_tmag_dFj +
                     (lsi->d_normal_dF[a][j] - d_dp_dFj * nw[a]) * tmaginv;

      for (a = 0, d_wet_vector_dF[a][j] = 0.0; a < dim; a++)
        d_wet_vector_dF[a][j] = d_t_dFj[a] * wet_speed + t[a] * d_wet_speed_dFj;
    }

    var = MESH_DISPLACEMENT1;
    if (pd->v[pg->imtrx][var]) {
      for (b = 0; b < dim; b++) {
        var = MESH_DISPLACEMENT1 + b;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (d_dp_dXj[b] = 0.0, a = 0; a < dim; a++) {
            d_dp_dXj[b] += lsi->d_normal_dmesh[a][b][j] * nw[a];
          }

          d_tmag_dXj[b] = -tmaginv * dp * d_dp_dXj[b];

          for (a = 0, d_t_dXj[a][b] = 0.0; a < dim; a++) {
            d_t_dXj[a][b] = (nf[a] - dp * nw[a]) * (-pow(tmaginv, 2)) * d_tmag_dXj[b] +
                            (lsi->d_normal_dmesh[a][b][j] - d_dp_dXj[b] * nw[a]) * tmaginv;
          }

          for (a = 0, d_wet_vector_dX[a][b][j] = 0.0; a < dim; a++) {
            switch (bc_type) {
            case BLAKE_DIRICHLET_BC:
              d_wet_speed_dXj = v0 * cosh(g * (cos_ca - cos_ca_static)) * g * d_dp_dXj[b];
              break;
            case HOFFMAN_DIRICHLET_BC:
              d_wet_speed_dXj = g / liq_visc * (1. + 1.31 * pow(g_dca, 0.99)) /
                                (1. - 1.31 * 0.99 * A_dca / pow(g_dca, 0.01)) * pow(A_dca, 0.294) *
                                (4.) / (0.706 * 2 * 5.16 * (3. - costheta) * (1. + costheta));
              d_wet_speed_dXj *= -d_dp_dXj[b];
              d_wet_speed_dXj += wet_speed * (-1. / liq_visc) * d_mu->X[b][j] * visc_factor;
              break;
            case COX_DIRICHLET_BC:
              /* sensitivity wrt to q_outer has not been included	*/
              f_den = log(reciprocal_slip) + (q_inner / f_sca - q_outer / f_dca);
              d_wet_speed_dXj = g / liq_visc * (-1. / (f_dca * (-sintheta))) / f_den;
              d_wet_speed_dXj *= -d_dp_dXj[b];
              d_wet_speed_dXj += wet_speed * (-1. / liq_visc) * d_mu->X[b][j] * visc_factor;
              break;
            case SHIK_DIRICHLET_BC:
              dveloc0_ddpj =
                  (-d_dp_dXj[b]) *
                  (SQUARE(theta) + theta * sintheta * costheta - 2 * SQUARE(sintheta)) /
                  (SQUARE(sintheta * costheta) - 2 * theta * sintheta * costheta + SQUARE(theta));
              drhs_num_ddpj = d_dp_dXj[b];
              drhs_den_ddpj = (v0 - 1) * dveloc0_ddpj - d_dp_dXj[b];
              if (theta_max > theta_max_degrees * M_PIE / 180.)
                drhs_den_ddpj *= (1. - 1. / shik_max_factor);
              drhs_ddpj = (rhs_den * drhs_num_ddpj - rhs_num * drhs_den_ddpj) / SQUARE(rhs_den);
              d_wet_speed_dXj =
                  0.5 * sqrt(g * v0) * (1 - 0.5 * rhs) * drhs_ddpj / pow(1 + rhs, 1.5);
              break;
            } /* end of switch	*/
            d_wet_vector_dX[a][b][j] = d_t_dXj[a][b] * wet_speed + t[a] * d_wet_speed_dXj;
          }
        }
      }
    }
    for (a = 0; a < pd->Num_Dim; a++) {
      func[a] = betainv * (wet_speed * triangle * t[a]);
    }

    if (af->Assemble_Jacobian) {

      var = ls->var;
      for (a = 0; a < dim; a++) {
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            if (triangle_form) {
              d_func[a][var][j] += betainv * (d_wet_vector_dF[a][j] * triangle +
                                              wet_speed * t[a] * d_triangle_dF * bf[var]->phi[j]);
            } else {
              d_func[a][var][j] +=
                  betainv * (d_wet_vector_dF[a][j] * triangle +
                             wet_speed * t[a] * 0.5 * lsi->d_delta_dF[j] * wetting_length);
            }
            d_func[a][var][j] += -SQUARE(betainv) * d_beta_dF[j] * wet_speed * triangle * t[a];
          }
        }
      }
      var = MESH_DISPLACEMENT1;
      if (pd->v[pg->imtrx][var]) {
        for (b = 0; b < dim; b++) {
          var = MESH_DISPLACEMENT1 + b;

          for (a = 0; a < dim; a++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_func[a][var][j] += betainv * triangle * d_wet_vector_dX[a][b][j];
            }
          }
        }
      }
      /* velocity senstivities		*/
      var = VELOCITY1;
      if (pd->v[pg->imtrx][var]) {
        for (b = 0; b < dim; b++) {
          var = VELOCITY1 + b;
          for (a = 0; a < dim; a++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              switch (bc_type) {
              case BLAKE_DIRICHLET_BC:
                d_wet_vector_dpj = 0;
                break;
              case HOFFMAN_DIRICHLET_BC:
                d_wet_vector_dpj =
                    (g_dca - g_sca) *
                    (liq_visc * surf_tens * mp->d_surface_tension[var] * bf[var]->phi[j] -
                     g * d_mu->v[b][j] * visc_factor) /
                    SQUARE(liq_visc);
                break;
              case COX_DIRICHLET_BC:
                f_den = log(reciprocal_slip) + (q_inner / f_sca - q_outer / f_dca);
                d_wet_vector_dpj =
                    g_integral / f_den *
                    (liq_visc * surf_tens * mp->d_surface_tension[var] * bf[var]->phi[j] -
                     g * d_mu->v[b][j] * visc_factor) /
                    SQUARE(liq_visc);
                break;
              case SHIK_DIRICHLET_BC:
                d_wet_vector_dpj =
                    0.5 * wet_speed / g * surf_tens * mp->d_surface_tension[var] * bf[var]->phi[j];
                break;
              }
              d_func[a][var][j] += betainv * triangle * t[a] * d_wet_vector_dpj;
            }
          }
        }
      }
      /* Temperature senstivities		*/
      var = TEMPERATURE;
      if (pd->v[pg->imtrx][var]) {
        for (a = 0; a < dim; a++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            switch (bc_type) {
            case BLAKE_DIRICHLET_BC:
              d_wet_vector_dpj = v0 * cosh((costhetaeq - costheta) * g) * (costhetaeq - costheta) *
                                 surf_tens * mp->d_surface_tension[var] * bf[var]->phi[j];
              break;
            case HOFFMAN_DIRICHLET_BC:
              d_wet_vector_dpj =
                  (g_dca - g_sca) *
                  (liq_visc * surf_tens * mp->d_surface_tension[var] * bf[var]->phi[j] -
                   g * d_mu->T[j] * visc_factor) /
                  SQUARE(liq_visc);
              break;
            case COX_DIRICHLET_BC:
              f_den = log(reciprocal_slip) + (q_inner / f_sca - q_outer / f_dca);
              d_wet_vector_dpj =
                  g_integral / f_den *
                  (liq_visc * surf_tens * mp->d_surface_tension[var] * bf[var]->phi[j] -
                   g * d_mu->T[j] * visc_factor) /
                  SQUARE(liq_visc);
              break;
            case SHIK_DIRICHLET_BC:
              d_wet_vector_dpj =
                  0.5 * wet_speed / g * surf_tens * mp->d_surface_tension[var] * bf[var]->phi[j];
              break;
            }
            d_func[a][var][j] += betainv * triangle * t[a] * d_wet_vector_dpj;
          }
        }
      }
      /*** concentration sensitivities	**/

      var = MASS_FRACTION;
      if (pd->v[pg->imtrx][var]) {
        for (a = 0; a < dim; a++) {
          for (w = 0; w < pd->Num_Species_Eqn; w++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              switch (bc_type) {
              case BLAKE_DIRICHLET_BC:
                d_wet_vector_dpj = v0 * cosh((costhetaeq - costheta) * g) *
                                   (costhetaeq - costheta) * surf_tens *
                                   mp->d_surface_tension[MAX_VARIABLE_TYPES + w] * bf[var]->phi[j];
                break;
              case HOFFMAN_DIRICHLET_BC:
                d_wet_vector_dpj =
                    (g_dca - g_sca) *
                    (liq_visc * surf_tens * mp->d_surface_tension[MAX_VARIABLE_TYPES + w] *
                         bf[var]->phi[j] -
                     g * d_mu->C[w][j] * visc_factor) /
                    SQUARE(liq_visc);
                break;
              case COX_DIRICHLET_BC:
                f_den = log(reciprocal_slip) + (q_inner / f_sca - q_outer / f_dca);
                d_wet_vector_dpj =
                    g_integral / f_den *
                    (liq_visc * surf_tens * mp->d_surface_tension[MAX_VARIABLE_TYPES + w] *
                         bf[var]->phi[j] -
                     g * d_mu->C[w][j] * visc_factor) /
                    SQUARE(liq_visc);
                break;
              case SHIK_DIRICHLET_BC:
                d_wet_vector_dpj = 0.5 * wet_speed / g * surf_tens *
                                   mp->d_surface_tension[MAX_VARIABLE_TYPES + w] * bf[var]->phi[j];
                break;
              }
              d_func[a][MAX_VARIABLE_TYPES + w][j] += betainv * triangle * t[a] * d_wet_vector_dpj;
            }
          }
        }
      } /*  end of concentration sensitivities	*/
    }   /* af->Assemble_Jacobian */
  }     /*  lsi->near		*/

  /* Elsewhere on the boundary we apply the slip velocity as a Dirichlet condition */

  for (a = 0; a < pd->Num_Dim; a++) {
    func[a] += betainv * (-tau * fv_dot->v[a] + vs[a] - fv->v[a]);
  }

  if (af->Assemble_Jacobian) {
    if (pd->v[pg->imtrx][VELOCITY1]) {
      for (a = 0; a < pd->Num_Dim; a++) {
        var = VELOCITY1 + a;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          d_func[a][var][j] += betainv * (-tau * (1. + 2. * tt) * phi_j / dt - phi_j);
        }
      }
    }
    var = ls->var;
    for (a = 0; a < dim; a++) {
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_func[a][var][j] +=
              -SQUARE(betainv) * d_beta_dF[j] * (-tau * fv_dot->v[a] + vs[a] - fv->v[a]);
        }
      }
    }
  }
}

void apply_wetting_tension(double func[MAX_PDIM],
                           double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                           const double wetting_tension) /* wetting tension */

/******************************************************************************
 *
 *  Function which applies wetting tension to fluid momentum equations in vicinity of
 *   zero level set.
 *
 *            Author: Tom Baer ( 6/2001)
 *
 ******************************************************************************/

{
  int p;

  int near_ls;
  double H_ls, delta_ls, normal_ls[MAX_PDIM] = {0., 0., 0.};

  /***************************** EXECUTION BEGINS ******************************/

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  level_set_interface(fv->F, fv->grad_F, ls->Length_Scale, 0, &near_ls, &H_ls, NULL, NULL,
                      &delta_ls, NULL, NULL, normal_ls, NULL, NULL);
  if (near_ls) {
    if (pd->Num_Dim == 2) {
      double t[2], dot_prod = 0., sign;

      sign = 1.0;
      if ((normal_ls[0] * fv->snormal[1] - normal_ls[1] * fv->snormal[0]) < 0.0)
        sign = -1.0;

      t[0] = sign * fv->snormal[1];
      t[1] = -sign * fv->snormal[0];

      dot_prod = dot_product(pd->Num_Dim, t, normal_ls);

      for (p = 0; p < pd->Num_Dim; p++) {
        func[p] = delta_ls * wetting_tension * dot_prod * t[p];
      }
    } else {
      double t[3], l[3];

      cross_really_simple_vectors(normal_ls, fv->snormal, l);

      normalize_really_simple_vector(l, pd->Num_Dim);

      cross_really_simple_vectors(fv->snormal, l, t);

      normalize_really_simple_vector(t, pd->Num_Dim);

      for (p = 0; p < pd->Num_Dim; p++) {
        func[p] = delta_ls * wetting_tension * t[p];
      }
      /* DPRINTF(stderr," F = %6.3f ,    func[0,1,2] = %6.3f %6.3f %6.3f \n", fv->F, func[0],
       * func[1], func[2] ); */
    }
  }

  return;
}

/*****************************************************************************/

void ftmelt_bc(double func[DIM],
               double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
               const double a1,              /*  function parameters from data card       */
               const double x_dot[MAX_PDIM], /* mesh velocity vector             */
               const dbl tt,                 /* parameter to vary time integration from   *
                                              * explicit (tt = 1) to implicit (tt = 0)    */
               const dbl dt)                 /* current value of the time step            */

/******************************************************************************
 *
 *  Function which evaluates the melting point distinguishing condition:
 *  n.grad(phi)*(t - tmp)=0
 *
 *            Author: P. R. Schunk    (9/6/96)
 *
 ******************************************************************************/

{

  /* Local variables */

  int j, kdir, var, p;
  double phi_j;
  /***************************** EXECUTION BEGINS *******************************/

  if (af->Assemble_LSA_Mass_Matrix) {
    for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
      var = MESH_DISPLACEMENT1 + kdir;
      if (pd->v[pg->imtrx][var])
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          d_func[0][var][j] -= phi_j * fv->snormal[kdir];
        }
    }
    return;
  }

  if (af->Assemble_Jacobian) {
    for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
      for (p = 0; p < pd->Num_Dim; p++) {
        var = MESH_DISPLACEMENT1 + p;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            d_func[0][var][j] += (fv->v[kdir] - x_dot[kdir]) * fv->dsnormal_dx[kdir][p][j];
            if (TimeIntegration != 0 && p == kdir) {
              d_func[0][var][j] +=
                  (-(1. + 2. * tt) * phi_j / dt) * fv->snormal[kdir] * delta(p, kdir);
            }
          }
        }
      }

      /*
       * Did not include particle stuff here.
       */
      var = VELOCITY1 + kdir;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          d_func[0][var][j] += phi_j * fv->snormal[kdir];
        }
      }
    }
  } /* end of if Assemble_Jacobian */

  /* Calculate the residual contribution	*/

  for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
    *func += (fv->T - a1) * fv->snormal[kdir];
  }

} /* END of routine tmelt_bc  */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void continuous_tangent_velocity(double func[DIM],
                                 double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                                 const int ielem_dim) /* dimension of element     */
{
  int j, kdir, var, p;
  double phi_j; /* Basis Function			      */

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  /*  initialize variables */

  if (pd->Num_Dim == 2) {
    if (af->Assemble_Jacobian) {
      for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
        for (p = 0; p < pd->Num_Dim; p++) {
          var = MESH_DISPLACEMENT1 + p;
          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_func[0][var][j] += fv->v[kdir] * fv->dstangent_dx[0][kdir][p][j];
            }
          }
        }

        /*
         * Didn't add particle stuff here either.
         */

        var = VELOCITY1 + kdir;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            d_func[0][var][j] += phi_j * fv->stangent[0][kdir];
          }
        }
      }
    } /* end of if Assemble_Jacobian */
  }   /* end if Num_Dim */

  /* Calculate the residual contribution	*/
  if (pd->Num_Dim == 2) {
    for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
      *func += fv->v[kdir] * fv->stangent[0][kdir];
    }

  } else if (pd->Num_Dim == 3) {
    GOMA_EH(GOMA_ERROR, " CONT_TANG_VEL in 3D is ill-defined");
  }
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void continuous_normal_velocity(double func[DIM],
                                double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                                const int ielem_dim) /* dimension of element     */
{
  int j, kdir, var, p;
  double phi_j; /* Basis Function			      */

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  /*  initialize variables */

  if (pd->Num_Dim == 2) {
    if (af->Assemble_Jacobian) {
      for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
        for (p = 0; p < pd->Num_Dim; p++) {
          var = MESH_DISPLACEMENT1 + p;
          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_func[0][var][j] += fv->v[kdir] * fv->dsnormal_dx[kdir][p][j];
            }
          }
        }

        /*
         * I didn't add particle stuff here.
         */

        var = VELOCITY1 + kdir;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            d_func[0][var][j] += phi_j * fv->snormal[kdir];
          }
        }
      }
    } /* end of if Assemble_Jacobian */
  }   /* end if Num_Dim */

  /* Calculate the residual contribution	*/
  if (pd->Num_Dim == 2) {
    for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
      *func += fv->v[kdir] * fv->snormal[kdir];
    }

  } else if (pd->Num_Dim == 3) {
    GOMA_EH(GOMA_ERROR, " CONT_NORM_VEL in 3D is ill-defined");
  }
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void discontinuous_velocity(double func[DIM],
                            double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                            double x_dot[MAX_PDIM],  /* mesh velocity vector                      */
                            const int mode,          /* Evporation or dissolution                 */
                            const int eb_mat_liquid, /* elem blk id for the liquid phase          */
                            const int eb_mat_gas,    /* elem blk id for the gas phase             */
                            const dbl tt,            /* parameter to vary time integration from   *
                                                      * explicit (tt = 1) to implicit (tt = 0)    */
                            const dbl dt)            /* current value of the time step            */
{
  int idblock = 0, j, kdir, var, p, w;
  double phi_j;                      /* Basis Function			     */
  double conc_insoluble_species = 0; /* concentration of insoluble species   *
                                      * like air                             */
  double diffusion_flux[MAX_CONC];   /* diffusion flux normal to interface of volatile species */
  double d_diffusion_flux_dc;        /* diffusion flux normal to interface of volatile species */
  double d_diffusion_flux_dx[MAX_CONC][DIM][MDE];
  /* derivative of diffusion flux normal to interface of volatile species */
  double vnormal = 0;

  /*  initialize variables */

  /* Calculate the residual contribution	*/
  if (mode == EVAPORATION)
    idblock = eb_mat_gas;
  else if (mode == DISSOLUTION)
    idblock = eb_mat_liquid;
  else
    GOMA_EH(GOMA_ERROR, "Don't recognize your mode on the DISCONTINUOUS_VELO BC");

  /* only apply if within gas phase */
  if (Current_EB_ptr->Elem_Blk_Id == idblock) {
    vnormal = 0.;
    for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
      vnormal += mp->density * (fv->v[kdir] - x_dot[kdir]) * fv->snormal[kdir];
    }
    *func = vnormal;
    if (mode == EVAPORATION) {
      conc_insoluble_species = 1.;
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        conc_insoluble_species -= fv->c[w];
      }

      for (w = 0; w < pd->Num_Species_Eqn; w++)
        diffusion_flux[w] = 0.;

      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
          diffusion_flux[w] += fv->grad_c[w][kdir] * fv->snormal[kdir];
        }
      }
    }
    if (mode == DISSOLUTION) {
      conc_insoluble_species = 0.;
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        conc_insoluble_species += fv->c[w];
      }

      for (w = 0; w < pd->Num_Species_Eqn; w++)
        diffusion_flux[w] = 0.;

      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
          diffusion_flux[w] -= fv->grad_c[w][kdir] * fv->snormal[kdir];
        }
      }
    }

    *func *= conc_insoluble_species;

    /* for now assume constant diffusivity */
    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      *func += mp->density * mp->diffusivity[w] * diffusion_flux[w];
    }
  }

  else {
    /* if not in the correct phase for this mode, do nothing for now */
  }

  if (af->Assemble_LSA_Mass_Matrix) {
    if (Current_EB_ptr->Elem_Blk_Id == idblock)
      for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
        var = MESH_DISPLACEMENT1 + kdir;
        if (pd->v[pg->imtrx][var])
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            d_func[0][var][j] -= conc_insoluble_species * mp->density * phi_j * fv->snormal[kdir];
          }
      }
    return;
  }

  if (af->Assemble_Jacobian) {
    if (Current_EB_ptr->Elem_Blk_Id == idblock) {
      for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
        for (w = 0; w < pd->Num_Species_Eqn; w++)
          for (p = 0; p < pd->Num_Dim; p++)
            for (j = 0; j < MDE; j++) {
              {
                { d_diffusion_flux_dx[w][p][j] = 0.; }
              }
            }
        for (p = 0; p < pd->Num_Dim; p++) {
          var = MESH_DISPLACEMENT1 + p;
          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];
              d_func[0][var][j] += conc_insoluble_species * mp->density *
                                   (fv->v[kdir] - x_dot[kdir]) * fv->dsnormal_dx[kdir][p][j];
              if (TimeIntegration != 0 && p == kdir) {
                d_func[0][var][j] += conc_insoluble_species * mp->density *
                                     (-(1. + 2. * tt) * phi_j / dt) * fv->snormal[kdir] *
                                     delta(p, kdir);
              }
              if (mode == EVAPORATION) {
                for (w = 0; w < pd->Num_Species_Eqn; w++) {
                  d_diffusion_flux_dx[w][p][j] +=
                      fv->d_grad_c_dmesh[kdir][w][p][j] * fv->snormal[kdir] +
                      fv->grad_c[w][kdir] * fv->dsnormal_dx[kdir][p][j];
                }
              } else if (mode == DISSOLUTION) {
                for (w = 0; w < pd->Num_Species_Eqn; w++) {
                  d_diffusion_flux_dx[w][p][j] -=
                      fv->d_grad_c_dmesh[kdir][w][p][j] * fv->snormal[kdir] +
                      fv->grad_c[w][kdir] * fv->dsnormal_dx[kdir][p][j];
                }
              }
              for (w = 0; w < pd->Num_Species_Eqn; w++) {
                d_func[0][var][j] +=
                    mp->density * mp->diffusivity[w] * d_diffusion_flux_dx[w][p][j];
              }
            }
          }
        }

        /*
         * Didn't put particle stuff here.
         */

        var = VELOCITY1 + kdir;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            d_func[0][var][j] += conc_insoluble_species * mp->density * phi_j * fv->snormal[kdir];
          }
        }
      }
      var = MASS_FRACTION;
      if (pd->v[pg->imtrx][var]) {

        for (w = 0; w < pd->Num_Species_Eqn; w++) {

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            if (mode == EVAPORATION) {
              d_diffusion_flux_dc = 0.;
              for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
                d_diffusion_flux_dc += mp->density * mp->diffusivity[w] *
                                       bf[var]->grad_phi[j][kdir] * fv->snormal[kdir];
              }
              d_func[0][MAX_VARIABLE_TYPES + w][j] += (-vnormal * phi_j + d_diffusion_flux_dc);
            }
            if (mode == DISSOLUTION) {
              d_diffusion_flux_dc = 0.;
              for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
                d_diffusion_flux_dc -= mp->density * mp->diffusivity[w] *
                                       bf[var]->grad_phi[j][kdir] * fv->snormal[kdir];
              }
              d_func[0][MAX_VARIABLE_TYPES + w][j] += (vnormal * phi_j + d_diffusion_flux_dc);
            }
          }
        }
      }
    }

    else {
      /* Do nothing for now */
    }
  }
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
#define SIGNN 1.

void fnormal_stress_bc(double func[DIM],
                       double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                       const double stress_normal, /* normal stress               */
                       const dbl relax)            /* relaxation parameters                  */

/****************************************************************************
 *
 *   Computes the residual and sensitivities of the constraint
 *
 *                 nn:T = stress_normal
 *
 *   where T is the fluid stress and n the normal to the surface
 *
 *
 *   Author: TA Baer ( 12/97 )
 *
 ******************************************************************************/

{
  int j, var, p, q, a;

  double T[DIM][DIM];
  double dT_dv[DIM][DIM][DIM][MDE];
  double dT_dmesh[DIM][DIM][DIM][MDE];

  /* compute fluid stress and its sensitivities */

  memset(T, 0, sizeof(double) * DIM * DIM);

  *func = -stress_normal;

  for (p = 0; p < pd->Num_Dim; p++) {
    T[p][p] = -fv->P;
    for (q = 0; q < pd->Num_Dim; q++) {
      T[p][q] += gn->mu0 * (fv->grad_v[p][q] + fv->grad_v[q][p]);

      for (a = 0; a < pd->Num_Dim; a++) {
        var = MESH_DISPLACEMENT1 + a;

        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          dT_dmesh[p][q][a][j] = fv->d_grad_v_dmesh[p][q][a][j] + fv->d_grad_v_dmesh[q][p][a][j];
          dT_dmesh[p][q][a][j] *= gn->mu0;
        }
      }

      /*
       * Didn't add particle stuff here.
       */

      for (a = 0; a < pd->Num_Dim; a++) {
        var = VELOCITY1 + a;

        for (j = 0; pd->v[pg->imtrx][var] && j < ei[pg->imtrx]->dof[var]; j++) {
          dT_dv[p][q][a][j] = bf[var]->grad_phi_e[j][a][p][q] + bf[var]->grad_phi_e[j][a][q][p];
          dT_dv[p][q][a][j] *= gn->mu0;
        }
      }

      /* compute residual while were at it */

      *func += fv->snormal[p] * fv->snormal[q] * T[p][q];
    }
  }

  *func *= SIGNN;

  if (af->Assemble_Jacobian) {
    for (p = 0; p < pd->Num_Dim; p++) {
      for (q = 0; q < pd->Num_Dim; q++) {

        for (a = 0; a < pd->Num_Dim; a++) {
          var = MESH_DISPLACEMENT1 + a;

          for (j = 0; pd->v[pg->imtrx][var] && j < ei[pg->imtrx]->dof[var]; j++) {

            d_func[0][var][j] += fv->dsnormal_dx[p][a][j] * fv->snormal[q] * T[p][q];
            d_func[0][var][j] += fv->dsnormal_dx[q][a][j] * fv->snormal[p] * T[p][q];
            d_func[0][var][j] += fv->snormal[p] * fv->snormal[q] * dT_dmesh[p][q][a][j];

            d_func[0][var][j] *= SIGNN;
          }
        }

        for (a = 0; a < pd->Num_Dim; a++) {
          var = VELOCITY1 + a;

          for (j = 0; pd->v[pg->imtrx][var] && j < ei[pg->imtrx]->dof[var]; j++) {
            d_func[0][var][j] += fv->snormal[p] * fv->snormal[q] * dT_dv[p][q][a][j];

            d_func[0][var][j] *= SIGNN;
          }
        }
      }

      var = PRESSURE;

      for (j = 0; pd->v[pg->imtrx][var] && j < ei[pg->imtrx]->dof[var]; j++) {
        d_func[0][var][j] += -fv->snormal[p] * fv->snormal[p] * bf[var]->phi[j];
        d_func[0][var][j] *= SIGNN;
      }
    }
  }
}
/***************************************************************************/
/***************************************************************************/
/*****************************************************************************/

void qside_directional(double func[DIM],
                       double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                       const double qx, /* prescribed heat flux in x direction */
                       const double qy, /* prescribed heat flux in y direction */
                       const double qz) /* prescribed heat flux in z direction */

/***********************************************************************
 *
 * qside_directonal():
 *
 *  Function which evaluates the expression specifying the
 *  heat flux  boundary condition at a quadrature point on a side
 *  of an element.
 *
 *         func =     n . qside_of_single_direction
 *
 *  Note: this bc is similar to qside except it is
 *        directional dependent.  qside replaces
 *        n dot q, but qside_directional does not.
 *
 * Input: qx, qy, and qz
 *
 *  vnormal = specified on the bc card as the first float
 *
 * Output:
 *
 *  func[0] = value of the function mentioned above
 *  d_func[0][varType][lvardof] =
 *              Derivate of func[0] wrt
 *              the variable type, varType, and the local variable
 *              degree of freedom, lvardof, corresponding to that
 *              variable type.
 *
 *   Author: ACSUN    (4/23/03)
 ********************************************************************/
{
  int j, kdir, var, p;
  double qside[DIM];
  double phi_j;

  qside[0] = qx;
  qside[1] = qy;
  qside[2] = qz;

  if (af->Assemble_LSA_Mass_Matrix) {
    for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
      var = MESH_DISPLACEMENT1 + kdir;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          d_func[0][var][j] -= phi_j * fv->snormal[kdir];
        }
      }
    }
    return;
  }

  if (af->Assemble_Jacobian) {

    for (kdir = 0; kdir < pd->Num_Dim; kdir++) {

      for (p = 0; p < pd->Num_Dim; p++) {
        var = MESH_DISPLACEMENT1 + p;

        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

            phi_j = bf[var]->phi[j];
            d_func[0][var][j] += qside[kdir] * fv->dsnormal_dx[kdir][p][j];
          }
        }
      }
    }

  } /* end of if Assemble_Jacobian */

  /* Calculate the residual contribution	*/

  func[0] = 0.;
  for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
    func[0] += qside[kdir] * fv->snormal[kdir];
  }
} /* END of routine qside_directional_bc  */
/*****************************************************************************/
/*****************************************************************************/

void q_velo_slip_bc(double func[MAX_PDIM],
                    double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                    int ibc,
                    double x[MAX_PDIM],
                    const double xsurf[MAX_PDIM], /* coordinates of surface Gauss  *
                                                   * point, i.e. current position  */
                    const double tt,
                    const double dt)

/*****************************************************************************
 *  Function which calculates the surface integral for viscous heating due
 *  to navier-slip.
 *
 ******************************************************************************/
{
  int var, p, j;
  double phi_j;
  double slip_stress[MAX_PDIM];
  double d_slip_stress[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE];
  double omega, X_0[2], vs[3], vslip[3];

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  /* if based on VELO_SLIP_ROT, we need to compute surface velocity */
  /* Note: positive omega is CLOCKWISE */
  if (BC_Types[ibc].BC_Name == VELO_SLIP_ROT_BC) {
    omega = BC_Types[ibc].BC_Data_Float[1];
    X_0[0] = BC_Types[ibc].BC_Data_Float[2];
    X_0[1] = BC_Types[ibc].BC_Data_Float[3];
    vs[0] = omega * (fv->x[1] - X_0[1]);
    vs[1] = -omega * (fv->x[0] - X_0[0]);
    vs[2] = 0.;
#define THREADED_SLIP FALSE
#if THREADED_SLIP
    pitch = BC_Types[ibc].BC_Data_Float[4];
    for (p = 0; p < 2; p++) {
      vslip[p] = (vs[p] - fv->v[p]);
    }
    vmag = normalize_really_simple_vector(vslip, 2);
    radius =
        sqrt((fv->x[0] - X_0[0]) * (fv->x[0] - X_0[0]) + (fv->x[1] - X_0[1]) * (fv->x[1] - X_0[1]));
    vs[2] = vmag / radius / (2. * M_PIE) * pitch;
#endif
  } else {
    vs[0] = BC_Types[ibc].BC_Data_Float[1];
    vs[1] = BC_Types[ibc].BC_Data_Float[2];
    vs[2] = BC_Types[ibc].BC_Data_Float[3];
  }

  /* slip velocity */
  for (p = 0; p < pd->Num_Dim; p++) {
    vslip[p] = (vs[p] - fv->v[p]);
  }

  /* zero out slip sensitivities */
  memset(slip_stress, 0, sizeof(double) * MAX_PDIM);
  memset(d_slip_stress, 0, sizeof(double) * MAX_PDIM * (MAX_VARIABLE_TYPES + MAX_CONC) * MDE);

  /* use fvelo_slip to evaluate slip and derivatives */
  fvelo_slip_bc(slip_stress, d_slip_stress, x, (int)BC_Types[ibc].BC_Name,
                (int)BC_Types[ibc].max_DFlt, BC_Types[ibc].BC_Data_Float,
                (int)BC_Types[ibc].BC_Data_Int[0], xsurf, tt, dt);

  if (af->Assemble_Jacobian) {
    for (p = 0; p < VIM; p++) {
      var = VELOCITY1 + p;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          d_func[0][var][j] += d_slip_stress[p][var][j] * vslip[p] - slip_stress[p] * phi_j;
        }
      }
      var = PVELOCITY1 + p;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          d_func[0][var][j] += d_slip_stress[p][var][j] * vslip[p] - slip_stress[p] * phi_j;
        }
      }
    }
  }

  /* Calculate the residual contribution					     */

  for (p = 0; p < VIM; p++) {
    *func += slip_stress[p] * vslip[p];
  }

  return;
} /* END of routine q_velo_slip_bc   */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void qnobc_surf(double func[DIM],
                double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                double time)

/******************************************************************************
 *
 *  Function which calculates the surface integral for the "no bc" heat transfer boundary condition
 *
 ******************************************************************************/
{

  /* Local variables */

  int j, b, w, p;
  int var;

  double q[DIM];
  HEAT_FLUX_DEPENDENCE_STRUCT d_q_struct;
  HEAT_FLUX_DEPENDENCE_STRUCT *d_q = &d_q_struct;

  /***************************** EXECUTION BEGINS *******************************/

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  heat_flux(q, d_q, time);

  if (af->Assemble_Jacobian) {

    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        for (p = 0; p < VIM; p++) {
          d_func[0][var][j] += fv->snormal[p] * d_q->T[p][j];
        }
      }
    }

    var = FILL;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        for (p = 0; p < VIM; p++) {
          d_func[0][var][j] += fv->snormal[p] * d_q->F[p][j];
        }
      }
    }

    if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
      for (b = 0; b < VIM; b++) {
        var = MESH_DISPLACEMENT1 + b;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (p = 0; p < VIM; p++) {
            d_func[0][var][j] += fv->snormal[p] * d_q->X[p][b][j] + q[p] * fv->dsnormal_dx[p][b][j];
          }
        }
      }
    }

    var = MASS_FRACTION;
    if (pd->v[pg->imtrx][var]) {
      for (p = 0; p < VIM; p++) {
        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_func[0][var][j] += fv->snormal[p] * d_q->C[p][w][j];
          }
        }
      }
    }
  }

  /* Calculate the residual contribution	     			     */
  for (p = 0; p < VIM; p++) {
    *func += fv->snormal[p] * q[p];
  }

  return;
} /* END of routine qnobc_surf                                               */
/*****************************************************************************/

void potential_nobc_surf(double func[DIM],
                         double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                         double time)

/************************************************************************
 *
 *  Function which calculates the surface integral for the "no bc"
 *  heat transfer boundary condition
 *
 ************************************************************************/
{
  int dim;
  int j, a, b, p, q;
  int var;
  double grad_phi_j;
  dbl k; /* permittivity or conductivity */

  /***************************** EXECUTION BEGINS ****************************/

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  dim = pd->Num_Dim;
  if (mp->VoltageFormulation == V_PERMITTIVITY) {
    k = mp->permittivity;
  } else {
    k = mp->electrical_conductivity;
  }

  if (af->Assemble_Jacobian) {

    var = VOLTAGE;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        for (p = 0; p < VIM; p++) {
          grad_phi_j = bf[var]->grad_phi[j][p];
          d_func[0][var][j] += fv->snormal[p] * (-k * grad_phi_j);
        }
      }
    }

    if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
      for (b = 0; b < VIM; b++) {
        var = MESH_DISPLACEMENT1 + b;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (q = 0; q < dim; q++) {
            d_func[0][var][j] += fv->snormal[q] * (-k * fv->d_grad_V_dmesh[q][b][j]) +
                                 fv->dsnormal_dx[q][b][j] * (-k * fv->grad_V[q]);
          }
        }
      }
    }
  }

  /* Calculate the residual contribution					     */

  for (a = 0; a < VIM; a++) {
    *func += fv->snormal[a] * (-k * fv->grad_V[a]);
  }

} /* END of routine potential_nobc_surf                          */
/*ARGSUSED*/

/*****************************************************************************/

void acoustic_plane_transmission(double func[DIM],
                                 double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                                 double time,
                                 const int bc_type,
                                 const double bdy_impedance,
                                 const double bdy_absorption,
                                 const double bdy_incident_real,
                                 const double bdy_incident_imag,
                                 const int blk_id)
/******************************************************************************
 *  Function which calculates the surface integral for the "plane transmission"
 *	acoustic wave equation
 ******************************************************************************/
{
  /* Local variables */

  int j;
  int var, conj_var;

  double imped_inv = 0.0;
  double normal_velo = 0.0;
  /***************************** EXECUTION BEGINS *******************************/
  if (af->Assemble_LSA_Mass_Matrix)
    return;

  if (bc_type == APR_PLANE_TRANS_BC || bc_type == API_PLANE_TRANS_BC)
    imped_inv = 1. / bdy_impedance;
  if (bc_type == APR_VELOCITY_BC || bc_type == API_VELOCITY_BC)
    normal_velo = bdy_impedance;

  if (bc_type == APR_PLANE_TRANS_BC) {
    var = ACOUS_PIMAG;
    conj_var = ACOUS_PREAL;
  } else {
    var = ACOUS_PREAL;
    conj_var = ACOUS_PIMAG;
  }

  if (af->Assemble_Jacobian) {
    if (bc_type == APR_PLANE_TRANS_BC || bc_type == API_PLANE_TRANS_BC) {
      switch (var) {
      case ACOUS_PIMAG:
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_func[0][var][j] = imped_inv * bf[var]->phi[j];
          }
        }
        if (pd->v[pg->imtrx][conj_var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[conj_var]; j++) {
            d_func[0][conj_var][j] = imped_inv * bf[var]->phi[j] * bdy_absorption;
          }
        }
        break;
      case ACOUS_PREAL:
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_func[0][var][j] = -imped_inv * bf[var]->phi[j];
          }
        }
        if (pd->v[pg->imtrx][conj_var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[conj_var]; j++) {
            d_func[0][conj_var][j] = -imped_inv * bf[var]->phi[j] * bdy_absorption;
          }
        }
        break;
      }
    } /*  plane_trans if	*/
  }

  /* Calculate the residual contribution					     */
  switch (bc_type) {
  case APR_PLANE_TRANS_BC:
    *func = imped_inv * ((fv->api - 2. * bdy_incident_imag) -
                         bdy_absorption * (fv->apr - 2. * bdy_incident_real));
    break;
  case API_PLANE_TRANS_BC:
    *func = imped_inv * (-(fv->apr - 2. * bdy_incident_real) -
                         bdy_absorption * (fv->api - 2. * bdy_incident_imag));
    break;
  case APR_VELOCITY_BC:
    if (blk_id == -1 || blk_id == ei[pg->imtrx]->elem_blk_id)
      *func = -normal_velo;
    break;
  case API_VELOCITY_BC:
    if (blk_id == -1 || blk_id == ei[pg->imtrx]->elem_blk_id)
      *func = normal_velo;
    break;
  }
  return;
} /* END of routine acoustic_plane_transmission                             */
/****************************************************************************/
/*****************************************************************************/
void acoustic_nobc_surf(double func[DIM],
                        double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                        double time,
                        const int bc_type)

/******************************************************************************
 *  Function which calculates the surface integral for the "no bc" acoustic bc
 ******************************************************************************/
{
  /* Local variables */
  int j, b, w, p;
  int var, ac_var, ac_eqn;

  double q[DIM];
  ACOUSTIC_FLUX_DEPENDENCE_STRUCT d_q_struct;
  ACOUSTIC_FLUX_DEPENDENCE_STRUCT *d_q = &d_q_struct;

  /***************************** EXECUTION BEGINS ******************************/

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  if (bc_type == APR_NOBC_BC) {
    ac_eqn = R_ACOUS_PREAL;
    ac_var = ACOUS_PREAL;
  } else {
    ac_eqn = R_ACOUS_PIMAG;
    ac_var = ACOUS_PIMAG;
  }

  acoustic_flux(q, d_q, time, ac_eqn, ac_var);

  if (af->Assemble_Jacobian) {

    var = ac_var;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        for (p = 0; p < VIM; p++) {
          d_func[0][var][j] += fv->snormal[p] * d_q->P[p][j];
        }
      }
    }

    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        for (p = 0; p < VIM; p++) {
          d_func[0][var][j] += fv->snormal[p] * d_q->T[p][j];
        }
      }
    }

    var = FILL;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        for (p = 0; p < VIM; p++) {
          d_func[0][var][j] += fv->snormal[p] * d_q->F[p][j];
        }
      }
    }

    if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
      for (b = 0; b < VIM; b++) {
        var = MESH_DISPLACEMENT1 + b;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (p = 0; p < VIM; p++) {
            d_func[0][var][j] += fv->snormal[p] * d_q->X[p][b][j] + q[p] * fv->dsnormal_dx[p][b][j];
          }
        }
      }
    }

    var = MASS_FRACTION;
    if (pd->v[pg->imtrx][var]) {
      for (p = 0; p < VIM; p++) {
        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_func[0][var][j] += fv->snormal[p] * d_q->C[p][w][j];
          }
        }
      }
    }
  }

  /* Calculate the residual contribution					     */
  for (p = 0; p < VIM; p++) {
    *func += fv->snormal[p] * q[p];
  }

  return;
} /* END of routine acoustic_nobc_surf                                                */
/****************************************************************************/
/******************************************************************************/
double calculate_vapor_cool(const double p[], double *d_evap_loss, double time)

/*******************************************************************************
 *
 *  Function which calculates the kernal to the surface integral for user-define
 *          keyhole welding specific heat transfer model.
 *******************************************************************************/

{
  double T_scale, q_scale, T_boil, theta, vapor_c0 = 0.0, vapor_c1 = 0.0, vapor_c2 = 0.0,
                                          vapor_c3 = 0.0;
  double evap_loss;

  /***************************** EXECUTION BEGINS *******************************/
  T_scale = p[0];
  q_scale = p[1];
  T_boil = mp->melting_point_solidus; /* T_boil stored as T_solidus in mat-file */
  theta = fv->T - T_boil;

  /* Mike Kanouff's curve-fit for evaporation energy loss in W/m/m          */
  /*                                              assume iron if T_boil>2000 */
  /*                                                      ice if T_boil<2000 */
  /* New Version 3/13/03                          */
  if (T_boil > 2000.0 * T_scale) {
    if (theta > 0.0 && theta <= 170.0 * T_scale) {
      vapor_c0 = 0.; /* iron coefficients */
      vapor_c1 = 8.14373e5 * (1.0 / T_scale);
      vapor_c2 = -2.24831e3 * (1.0 / T_scale) * (1.0 / T_scale);
      vapor_c3 = 2.71683e1 * (1.0 / T_scale) * (1.0 / T_scale) * (1.0 / T_scale);
    }
    if (theta > 170.0 * T_scale) {
      vapor_c0 = -3.1036e8; /* iron coefficients */
      vapor_c1 = 3.2724e6 * (1.0 / T_scale);
      vapor_c2 = -1.8084e3 * (1.0 / T_scale) * (1.0 / T_scale);
      vapor_c3 = 2.7284e0 * (1.0 / T_scale) * (1.0 / T_scale) * (1.0 / T_scale);
    }
    /* vapor_c0 = 0.;       */ /* iron coefficients */
    /*vapor_c1 = 76.452*1.0e4*(1.0/T_scale);*/
    /*vapor_c2 = 0.2246*1.0e4*(1.0/T_scale)*(1.0/T_scale);*/
    /*vapor_c3 = 1.7816e-4*1.0e4*(1.0/T_scale)*(1.0/T_scale)*(1.0/T_scale);*/
  } else {
    vapor_c0 = 0.;
    vapor_c1 = 3.442e3 * 1.0e4 * (1.0 / T_scale); /* ice coefficients */
    vapor_c2 = 7.7214 * 1.0e4 * (1.0 / T_scale) * (1.0 / T_scale);
    vapor_c3 = 0.34523 * 1.0e4 * (1.0 / T_scale) * (1.0 / T_scale) * (1.0 / T_scale);
  }

  /* evaporization related quantities */
  if (theta < 0.) {
    evap_loss = 0.;
    *d_evap_loss = 0.;
  } else {
    evap_loss = q_scale * (vapor_c0 + vapor_c1 * theta + vapor_c2 * pow(theta, 2.0) +
                           vapor_c3 * pow(theta, 3.0));
    *d_evap_loss = q_scale * (vapor_c1 + 2. * vapor_c2 * theta + 3. * vapor_c3 * pow(theta, 2.0));
  }

  return (evap_loss);
}

/******************************************************************************/
void q_vapor(double func[DIM],
             double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
             const double p[],
             const double x[],
             const double time)

/*******************************************************************************
 *
 *  Function which calculates the kernal to the surface integral for user-define
 *          keyhole welding specific heat transfer model.
 *******************************************************************************/

{
  int j, var;
  double qvaporloss, d_evap_loss = 0., theta, T_boil;

  /***************************** EXECUTION BEGINS *******************************/
  T_boil = mp->melting_point_solidus; /* T_boil stored as T_solidus in mat-file */
  theta = fv->T - T_boil;

  /* Lets CALCULATE LASER VAPOR HEAT LOSS  RARR */
  qvaporloss = calculate_vapor_cool(p, &d_evap_loss, time);

  /* Evaluate sensitivity to temperature  */
  var = TEMPERATURE;
  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
    if (theta > 0.0) {
      d_func[0][var][j] -= (d_evap_loss)*bf[var]->phi[j];
    }
  }

  *func -= qvaporloss;
}

/******************************************************************************/
void qlaser_surf(double func[DIM],
                 double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                 const double p[],
                 const double x[],
                 const double time)

/*******************************************************************************
 *
 *  Function which calculates the kernal to the surface integral for user-define
 *          keyhole welding specific heat transfer model.
 *
 *    SEE calculate_laser_flux for description of input parameters
 *******************************************************************************/

{
  int var, j, jvar;
  double absorp, absorp_base, qlaser, d_laser_dx[MAX_PDIM];

  /***************************** EXECUTION BEGINS *******************************/

  memset(d_laser_dx, 0, sizeof(double) * MAX_PDIM);

  absorp_base = p[2];
  absorp = absorp_base;

  /* Lets CALCULATE LASER FLUX  RARR */
  /* WARNING! sending p as x until we find what method can be used for use_pth and level sets */
  qlaser = calculate_laser_flux(p, time, p, fv->snormal, 1, 0, d_laser_dx);

  /* Evaluate sensitivity to displacements d()/dx  */
  for (jvar = 0; jvar < ei[pg->imtrx]->ielem_dim; jvar++) {
    var = MESH_DISPLACEMENT1 + jvar;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        /* input laser flux only a function of r  */

        /* fix this */

        /*d_func[0][var][j] += absorp*qlaser*fv->dsnormal_dx[0][jvar][j]+
          absorp*d_laser_dx[jvar]*fv->snormal[0]*bf[var]->phi[j]*delta(jvar,1);*/
        d_func[0][var][j] += absorp * d_laser_dx[jvar] * bf[var]->phi[j] * delta(jvar, 1);
      }
    }
  }

  *func += absorp * qlaser;

  return;

} /* End of routine qlaser_surf                                                */
/*******************************************************************************/

double calculate_laser_flux(const double p[],
                            double time,
                            const double x[],
                            double normal[],
                            int sw_ale,
                            int sw_ls,
                            double d_laser_dx[])

/*******************************************************************************
 *
 *  Setup here is for DIMENSIONLESS keyhole weld simulation:
 *
 *  >>> SEE DOCUMENTATION (HARDCOPY) FOR EXPLANATION OF PARAMETERS  <<<
 *
 *  Since we now have Mike's ICE recoil pressure correlations I included both
 *  in this version of QUSER. Basically the value of the evaporation temperature
 *  is used to select either the iron correlation (iff T_boil>2000) or the ice
 *  correlation for T_boil<=2000. Thought that this might reduce recompiling the
 *  code to switch between two materials.
 *
 *  CAUTION - the 3D capability calculates NO view factors (yet) and assumes that
 *            the laser beam will be normal to one of the primary coordinate axis!
 *
 *  Input parameters:
 *
 *
 *      p[0] -- Laser Power                                          [q_laserpow]
 *              representing the nominal laser radiation.
 *
 *      p[1] -- Base value of Laser Power                            [base_qlaser]
 *
 *      p[2] -- Base value of surface absorptivity                   [absorp_base]
 *
 *      p[3] -- Switch to use velo_normal with laser radiation       [sw_vn]
 *              [0] -off, [1] -on
 *
 *      p[4] -- Cutoff time for laser pulse radiation                [t_cutoff]
 *
 *      p[5] -- Time at which laser radiation drops to 1/e           [t_tapper]
 *
 *      p[6] -- Overshoot value (%) of Laser Power at time 0         [qlaser0]
 *
 *      p[7] -- Radius of laser beam, used for non-dimensional calcs [radius_r]
 *
 *      p[8]-- Time for laser pulse to reach peak                   [t_deltpk]
 *
 *      p[9]-- Time for laser pulse to reach steady state           [t_deltst]
 *
 *      p[10]-- Switch to activate Exponential Function with         [sw_th_exp]
 *              Radial Position only in the x-y plane [1] or
 *              Exponential Function with Radial Position based on
 *              Absolute Distance from Beam Center[0]
 *
 *      p[11]-- Laser center location (x)                            [laser_x]
 *
 *      p[12]-- Laser center location (x)                            [laser_y]
 *
 *      p[13]-- Laser center location (x)                            [laser_z]
 *
 *      p[14]-- Laser axis card (x)                                  [laser_x_dir]
 *
 *      p[15]-- Laser axis card (x)                                  [laser_y_dir]
 *
 *      p[16]-- Laser axis card (x)                                  [laser_z_dir]
 *
 *      p[17]-- Overlapping Spot Weld Frequency                      [1.0/t_spot]
 *
 *      p[18]-- Overlapping Spot Weld Total Number                   [spot_total]
 *
 *      p[19]-- Set 1= travel weld, 0=spot weld, -1 for pseudo-spot  [sw_trv_spt]
 *
 *      p[20]-- Spot Weld Spacing                                    [spot_space]
 *
 *      p[21]-- Traveling Beam Radius Linear <0                      [trav_rad]
 *
 *      p[22]-- Use laser pass-thru model for gap simulation         [use_pth]
 *
 *      p[23]-- Initial time for start of beam for restarts          [t_travrad0]
 *
 *      p[24]-- Laser Beam Travel Speed in x                         [beam_u]
 *
 *      p[25]-- Laser Beam Travel Speed in y                         [beam_v]
 *
 *      p[26]-- Laser Beam Travel Speed in z                         [beam_w]
 *
 *
 *    Laser flux distribution of the form q=(input Watts/m^2) * f(r), where
 *       f(r)= 2. * R_eff * exp ( - R_eff * r^2)
 *    R_eff set by the energy_concentration input card as
 *       R_eff = - ln( 1 - e_concen)  where e_concen is set as 0.4522
 *
 ********************************************************************************/
{
  int ispot, ispot_total, sw_trv_spt, onnegx, ns_id1, idx1, use_pth, node_1;
  double q_laserpow, sw_vn, t_cutoff, t_tapper, qlaser0;
  double e_concen, radius_r, t_deltpk, t_deltst, sw_th_exp, laser_x, laser_y;
  double laser_z, laser_x_dir, laser_y_dir, laser_z_dir, t_spot, spot_space;
  double trav_rad, t_travrad0, beam_u, beam_v, beam_w, base_qlaser, initial_pos1;
  double delta_pos1, actual_pos1 = 0.0, t_0, rispot = 0.0, R_eff = 0.0;
  double x0c, x1c, x2c, center_x = 0.0, center_y = 0.0, center_z = 0.0, radial_pos, beam_vx,
                        beam_vy, beam_vz;
  double traveldist, darc, ddx, ddy, qlaser, dldx, dldy, dldz;

  /***************************** EXECUTION BEGINS *******************************/

  q_laserpow = p[0];
  base_qlaser = p[1];
  sw_vn = p[3];
  t_cutoff = p[4];
  t_tapper = p[5];
  /* Next card disabled until fix becomes available */
  /*  e_concen    = p[6];     */
  qlaser0 = p[6];
  e_concen = 0.4522;
  radius_r = p[7];
  t_deltpk = p[8];
  t_deltst = p[9];
  sw_th_exp = p[10];
  laser_x = p[11];
  laser_y = p[12];
  laser_z = p[13];
  laser_x_dir = p[14];
  laser_y_dir = p[15];
  laser_z_dir = p[16];
  t_spot = 1.0 / p[17];
  ispot_total = p[18];
  sw_trv_spt = p[19];
  spot_space = p[20];
  trav_rad = p[21];
  use_pth = p[22];
  t_travrad0 = p[23];
  beam_u = p[24];
  beam_v = p[25];
  beam_w = p[26];

  if (sw_trv_spt < -0.2) /*if pseudo-spot weld*/
  {
    spot_space = 0.0;
  } else if (sw_trv_spt > -0.1 && sw_trv_spt < 0.1) /*if spot weld*/
  {
    beam_u = 0.0;
    beam_v = 0.0;
    beam_w = 0.0;
  } else if (sw_trv_spt > 0.9 && sw_trv_spt < 1.1) /* if linear traveling weld */
  {
    spot_space = 0.0;
    t_spot = 0.0;
    ispot_total = 1;
  } else if (sw_trv_spt > 1.9 && sw_trv_spt < 2.1) /* if sinusoidal traveling weld */
  {
    spot_space = 0.0;
    ispot_total = 1;
  }

  /* This is a very-rough first cut at calculating a view factor (or laser pass-thru if */
  /*      a welding model with a gap is used with a top piece and a bottom piece */
  /*      seperated by DY */

  if (use_pth == 1 && sw_ale == 1) {
    ns_id1 = 1001;
    /* find the x location of 1001 */
    node_1 = psid2nn(ns_id1);
    GOMA_EH(node_1, "Could not find the nsid needed for keyhole calculation.");
    initial_pos1 = Coor[0][node_1];
    idx1 = Index_Solution(node_1, MESH_DISPLACEMENT1, 0, 0, -2, pg->imtrx);
    delta_pos1 = x[idx1];
    actual_pos1 = (initial_pos1 + delta_pos1);
  }

  /* --------------------------------------------------------------------*/
  /* set laser flux to appropriate power  */

  if (sw_trv_spt < 0.1) { /* spot weld */

    t_0 = 0.0;
    ispot = 0;
    rispot = 0.0;
    while (ispot <= ispot_total - 1) {
      if (time >= (t_0 + ispot * t_spot) && time < (t_deltpk + ispot * t_spot)) {
        q_laserpow = base_qlaser +
                     (qlaser0 - base_qlaser) * (time - (t_0 + ispot * t_spot)) / (t_deltpk - t_0);
      }
      if (time >= (t_deltpk + ispot * t_spot) && time < (t_deltst + ispot * t_spot)) {
        q_laserpow = qlaser0 + (q_laserpow - qlaser0) * (time - (t_deltpk + ispot * t_spot)) /
                                   (t_deltst - t_deltpk);
      }
      /*if (time >= (t_deltst+ispot*t_spot) && time < (t_cutoff+ispot*t_spot))
        {
          q_laserpow = q_laserpow;
        }*/
      if (time >= (t_cutoff + ispot * t_spot) && time < (t_tapper + ispot * t_spot)) {
        q_laserpow = q_laserpow + (base_qlaser - q_laserpow) *
                                      (time - (t_cutoff + ispot * t_spot)) / (t_tapper - t_cutoff);
      }
      if (time >= (t_tapper + ispot * t_spot) && time < (ispot + 1.0) * t_spot) {
        q_laserpow = base_qlaser;
      }
      if (time >= (ispot + 1) * t_spot) {
        rispot = ispot + 1.0;
      }
      ispot = ispot + 1;
    }
  } else if (sw_trv_spt > 0.9 && sw_trv_spt < 1.1) { /* linear traveling weld */
    t_0 = 0.0;
    if (time >= t_0 && time < t_deltpk) {
      q_laserpow = base_qlaser + (qlaser0 - base_qlaser) * (time - t_0) / (t_deltpk - t_0);
    }
    if (time >= t_deltpk && time < t_deltst) {
      q_laserpow = qlaser0 + (q_laserpow - qlaser0) * (time - t_deltpk) / (t_deltst - t_deltpk);
    }
    /*if (time >= t_deltst )
      {
        q_laserpow = q_laserpow;
      }*/
  } else if (sw_trv_spt > 1.9 && sw_trv_spt < 2.1) { /* sinusoidal traveling weld */
    t_0 = 0.0;
    /* for first cut q_laserpow is mean power */
    /* for first cut use t_spot to set frequency of sinusoid */
    /* so start at mean power */
    /*     but need a way to ramp laser power - so hardwire a ramp */
    /*     for first period use linear ramp to mean power          */
    if (time < t_spot) {
      q_laserpow = q_laserpow * time / t_spot;
    } else {
      q_laserpow = q_laserpow * (1.0 + sin(2.0 * 3.14159 * time / t_spot));
    }
    /* for first cut min power is base_qlaser */
    if (q_laserpow < base_qlaser) {
      q_laserpow = base_qlaser;
    }
  }

  if (use_pth == 1 && sw_ale == 1) {
    if (fv->x[0] < actual_pos1) {
      q_laserpow = 0.0;
    }
  }

  /* --------------------------------------------------------------------*/
  /* calculate exponent for laser flux distribution */
  R_eff = 3.0;
  if (e_concen > 0.10 && e_concen < 0.999) {
    R_eff = -log(1.00 - e_concen);
  }

  /* --------------------------------------------------------------------*/
  /* Now lets calculate center of beam and current position */

  if (ei[pg->imtrx]->ielem_dim == 2) {
    /* Current position */
    x0c = fv->x[0];
    x1c = fv->x[1];
    x2c = 0.0;
    /* calculate new center of laser */
    if (laser_z < 0) { /* danger - assumes that simulation is x-y and laser moves in z for 2D */
      center_x = sqrt(pow(laser_x_dir, 2)) * (laser_x);
      center_y = sqrt(pow(laser_x_dir, 2)) * (laser_y);
      center_z = sqrt(pow(laser_x_dir, 2)) *
                 (laser_z + beam_w * (time - t_travrad0) + spot_space * rispot);
      radial_pos =
          sqrt(0.0 * pow(x0c - center_x, 2) + pow(x1c - center_y, 2) + pow(x2c - center_z, 2)) /
          radius_r;
    } else { /* standard x-y */
      center_x = sqrt(pow(laser_y_dir, 2)) *
                 (laser_x + beam_w * (time - t_travrad0) + spot_space * rispot);
      center_y = sqrt(pow(laser_x_dir, 2)) *
                 (laser_y + beam_w * (time - t_travrad0) + spot_space * rispot);
      center_z = 0.;
      radial_pos = sqrt(pow(laser_x_dir, 2)) * sqrt(pow(x1c - center_y, 2)) / radius_r +
                   sqrt(pow(laser_y_dir, 2)) * sqrt(pow(x0c - center_x, 2)) / radius_r;
    }
  } else /* 3D! */
  {
    /* Current position */
    x0c = fv->x[0];
    x1c = fv->x[1];
    x2c = fv->x[2];
    /* 3D fix this */
    if (trav_rad < 0.0) {
      beam_vx = beam_u;
      beam_vy = beam_v;
      beam_vz = beam_w;
      center_x = (laser_x + beam_vx * (time - t_travrad0) + spot_space * rispot);
      center_y = (laser_y + beam_vy * (time - t_travrad0) + 0.0);
      center_z = (laser_z + beam_vz * (time - t_travrad0) + 0.0);
    }
    if (sw_th_exp > 0.5) { /* DANGER Assuming Beam Variation in only x-y plane */
      radial_pos =
          sqrt(pow(x0c - center_x, 2) + pow(x1c - center_y, 2) + 0.0 * pow(x2c - center_z, 2)) /
          radius_r;
    } else {
      radial_pos =
          sqrt(pow(x0c - center_x, 2) + pow(x1c - center_y, 2) + pow(x2c - center_z, 2)) / radius_r;
    }
  }

  /* For a traveling beam along a radii calculate current position */
  if (trav_rad > 0.0) {
    /* Assume beam center travels in x-y coordinates for now */
    /* must start at (0,travrad) with center at (0,0)!! */
    onnegx = 0;
    traveldist = beam_v * (time - t_travrad0) + spot_space * rispot;
    if (traveldist > 3.14159 * trav_rad && traveldist < 2.0 * 3.14159 * trav_rad) {
      traveldist = 2.0 * 3.14159 * trav_rad - beam_v * (time - t_travrad0) + spot_space * rispot;
      onnegx = 1;
    }
    darc = 2.0 * trav_rad * sin(traveldist / 2.0 / trav_rad);
    ddy = darc * darc / 2.0 / trav_rad;
    ddx = pow(darc * darc - ddy * ddy, 0.5);
    if (onnegx == 1) {
      center_x = laser_x - ddx;
    } else {
      center_x = laser_x + ddx;
    }
    center_y = laser_y - ddy;
    center_z = laser_z;

    if (sw_th_exp > 0.5) { /* DANGER Assuming Beam Variation in only x-y plane */
      radial_pos =
          sqrt(pow(x0c - center_x, 2) + pow(x1c - center_y, 2) + 0.0 * pow(x2c - center_z, 2)) /
          radius_r;
    } else {
      radial_pos =
          sqrt(pow(x0c - center_x, 2) + pow(x1c - center_y, 2) + pow(x2c - center_z, 2)) / radius_r;
    }
  }

  /* --------------------------------------------------------------------*/
  /* Now lets calculate distribution of qlaser */

  if (radial_pos <= 1.0) {
    /* use exponential */
    qlaser = 2. * R_eff * q_laserpow * exp(-R_eff * pow(radial_pos, 2));

    if (sw_th_exp > 0.5) {                   /* DANGER Assuming Beam Variation in-plane only */
      if (sqrt(pow(laser_x_dir, 2)) > 0.5) { /* DANGER assume x normal with this switch */
        if (sw_vn > 0.5) {
          qlaser = laser_x_dir * normal[0] * qlaser;
        } else {
          qlaser = laser_x_dir * qlaser;
        }
      }
      if (sqrt(pow(laser_y_dir, 2)) > 0.5) { /* DANGER assume y normal with this switch */
        if (sw_vn > 0.5) {
          qlaser = laser_y_dir * normal[1] * qlaser;
        } else {
          qlaser = laser_y_dir * qlaser;
        }
      }
      if (sqrt(pow(laser_z_dir, 2)) > 0.5) { /* DANGER assume z normal with this switch */
        if (sw_vn > 0.5) {
          qlaser = laser_z_dir * normal[2] * qlaser;
        } else {
          qlaser = laser_z_dir * qlaser;
        }
      }
    }
  } else {
    qlaser = 0.0;
  }

  /* To help handle normals to surface around an edge - use safety net */
  if (qlaser < 0)
    qlaser = 0.0;

  /* --------------------------------------------------------------------*/
  /* Now lets calculate derivatives*/

  dldx = 0.0;
  dldy = 0.0;
  dldz = 0.0;
  if (sqrt(pow(laser_x_dir, 2)) > 0.5) { /* have x normal component of beam so varies in y and z */
    if (x1c - center_y > 0.0) {
      dldy = -1.0;
    }
    if (x1c - center_y < 0.0) {
      dldy = 1.0;
    }
    if (x2c - center_z > 0.0) {
      dldz = -1.0;
    }
    if (x2c - center_z < 0.0) {
      dldz = 1.0;
    }
  }
  if (sqrt(pow(laser_y_dir, 2)) > 0.5) { /* have y normal component of beam so varies in x and z */
    if (x0c - center_x > 0.0) {
      dldx = -1.0;
    }
    if (x0c - center_x < 0.0) {
      dldx = 1.0;
    }
    if (x2c - center_z > 0.0) {
      dldz = -1.0;
    }
    if (x2c - center_z < 0.0) {
      dldz = 1.0;
    }
  }
  if (sqrt(pow(laser_z_dir, 2)) > 0.5) { /* have z normal component of beam so varies in x and y */
    if (x0c - center_x > 0.0) {
      dldx = -1.0;
    }
    if (x0c - center_x < 0.0) {
      dldx = 1.0;
    }
    if (x1c - center_y > 0.0) {
      dldy = -1.0;
    }
    if (x1c - center_y < 0.0) {
      dldy = 1.0;
    }
  }
  d_laser_dx[0] = dldx * (-2. * R_eff * radial_pos * qlaser);
  d_laser_dx[1] = dldy * (-2. * R_eff * radial_pos * qlaser);
  if (ei[pg->imtrx]->ielem_dim == 3)
    d_laser_dx[2] = dldz * (-2. * R_eff * radial_pos * qlaser);

  return (qlaser);
} /* END of routine calculate_laser_flux */
/*****************************************************************************/

void qrad_surf(double func[DIM],
               double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
               double heat_tran_coeff, /* Heat transfer coefficient (cgs units)   */
               double T_c,             /* bath temperature (Kelvin)	             */
               double epsilon,         /* emissivity                              */
               double sigma)           /* Boltzmann's constant                    */

/******************************************************************************
 *
 *  Function which calculates the surface integral for radiative - convective
 *  heat transfer.
 *
 ******************************************************************************/

{

  /* Local variables */

  int j_id;
  int var;
  double phi_j;

  /***************************** EXECUTION BEGINS *******************************/

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  if (af->Assemble_Jacobian) {

    /* sum the contributions to the global stiffness matrix */

    var = TEMPERATURE;
    for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
      phi_j = bf[var]->phi[j_id];
      d_func[0][var][j_id] -= heat_tran_coeff * phi_j;
      d_func[0][var][j_id] -= 4. * epsilon * sigma * pow(fv->T, 3.0) * phi_j;
    }
  }

  /* Calculate the residual contribution					     */

  *func += heat_tran_coeff * (T_c - fv->T) + epsilon * sigma * (pow(T_c, 4.0) - pow(fv->T, 4.0));

  return;
} /* END of routine qrad_surf                                                */

/*****************************************************************************/
void qside_contact_resis(double func[DIM],
                         double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                         int ID_mat_1, /* block ID material 1   */
                         int ID_mat_2, /* block ID material 2   */
                         double R_inv) /* resistivity */
/******************************************************************************
 *
 *  Function which applies thermal contact resistance at a side set interface
 *  Author: P. R. Schunk (3/1/2013)
 *
 ******************************************************************************/

{

  /* Local variables */

  int j_id;
  int var;
  double phi_j;
  double sign_int = 0;

  /***************************** EXECUTION BEGINS *******************************/

  if (af->Assemble_LSA_Mass_Matrix)
    return;
  if (Current_EB_ptr->Elem_Blk_Id == ID_mat_1) {
    sign_int = -1.0;
  } else if (Current_EB_ptr->Elem_Blk_Id == ID_mat_2) {
    sign_int = 1.0;
  } else {
    GOMA_EH(GOMA_ERROR, "T_CONTACT_RESIS has incorrect material ids");
  }

  if (af->Assemble_Jacobian) {

    /* sum the contributions to the global stiffness matrix */

    var = TEMPERATURE;
    for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
      phi_j = bf[var]->phi[j_id];
      d_func[0][var][j_id] += sign_int * R_inv * phi_j;
    }
  }

  /* Calculate the residual contribution					     */

  *func += sign_int * R_inv * (fv->T);

  return;
} /* END of routine qrad_surf                                                */

/****************************************************************************/
void qrad_surf_repulse(double func[DIM],
                       double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                       const double heat_tran_coeff, /* Heat transfer coefficient (cgs units)   */
                       const double T_c,             /* bath temperature (Kelvin)	             */
                       const double epsilon,         /* emissivity                              */
                       const double sigma,           /* Boltzmann's constant                    */
                       const double roll_rad,        /* Roll Radius                  */
                       const double origin[3],       /* roll axis origin (x,y,z) */
                       const double dir_angle[3],    /* axis direction angles */
                       const double hscale,          /* repulsion length scale */
                       const double repexp,          /* repulsive force exponent */
                       const double htc_roll,        /* repulsion coefficient */
                       const double T_roll           /* Roll temperature  */
                       )

/******************************************************************************
 *
 *  Function which calculates the surface integral for radiative - convective
 *  heat transfer.
 *
 ******************************************************************************/

{

  /* Local variables */

  int j_id;
  int var;
  double phi_j;
  int dim, jvar;      /* Degree of freedom counter */
  double d_dist[DIM]; /* distance derivatives  */
  double dist = 1e12; /* squared distance from surface to wall     */
  double factor;
  double coord[3], axis_pt[3], rad_dir[3], t, R;
  double force = 0.0, d_force = 0.0;
  int j;

  /***************************** EXECUTION BEGINS *******************************/

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  /*  initialize variables */
  dim = pd->Num_Dim;

  /* calculate distance from free surface to solid surface for repulsion calculations */

  coord[0] = fv->x[0];
  coord[1] = fv->x[1];
  if (dim == 3) {
    coord[2] = fv->x[2];
  } else {
    coord[2] = 0.0;
  }

  /*  find intersection of axis with normal plane - i.e., locate point on
   *          axis that intersects plane normal to axis that contains local point. */

  factor = SQUARE(dir_angle[0]) + SQUARE(dir_angle[1]) + SQUARE(dir_angle[2]);
  t = (dir_angle[0] * (coord[0] - origin[0]) + dir_angle[1] * (coord[1] - origin[1]) +
       dir_angle[2] * (coord[2] - origin[2])) /
      factor;
  axis_pt[0] = origin[0] + dir_angle[0] * t;
  axis_pt[1] = origin[1] + dir_angle[1] * t;
  axis_pt[2] = origin[2] + dir_angle[2] * t;

  /*  compute radius and radial direction */

  R = sqrt(SQUARE(coord[0] - axis_pt[0]) + SQUARE(coord[1] - axis_pt[1]) +
           SQUARE(coord[2] - axis_pt[2]));
  rad_dir[0] = (coord[0] - axis_pt[0]) / R;
  rad_dir[1] = (coord[1] - axis_pt[1]) / R;
  rad_dir[2] = (coord[2] - axis_pt[2]) / R;
  dist = R - roll_rad;
  d_dist[0] = rad_dir[0] * (1. - SQUARE(dir_angle[0]) / factor) +
              rad_dir[1] * (-dir_angle[1] * dir_angle[0] / factor) +
              rad_dir[2] * (-dir_angle[2] * dir_angle[0] / factor);
  d_dist[1] = rad_dir[1] * (1. - SQUARE(dir_angle[1]) / factor) +
              rad_dir[0] * (-dir_angle[0] * dir_angle[1] / factor) +
              rad_dir[2] * (-dir_angle[2] * dir_angle[1] / factor);
  d_dist[2] = rad_dir[2] * (1. - SQUARE(dir_angle[2]) / factor) +
              rad_dir[0] * (-dir_angle[0] * dir_angle[2] / factor) +
              rad_dir[1] * (-dir_angle[1] * dir_angle[2] / factor);

  /*  repulsion function  */
  force = htc_roll / pow(dist / hscale, repexp);
  d_force = -htc_roll * repexp / pow(dist / hscale, repexp + 1) / hscale;
  if (af->Assemble_Jacobian) {

    /* sum the contributions to the global stiffness matrix */

    var = TEMPERATURE;
    for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
      phi_j = bf[var]->phi[j_id];
      d_func[0][var][j_id] -= heat_tran_coeff * phi_j;
      d_func[0][var][j_id] -= force * phi_j;
      d_func[0][var][j_id] -= 4. * epsilon * sigma * pow(fv->T, 3.0) * phi_j;
    }
    /*
     *  Evaluate sensitivity to displacements d()/dx
     */
    for (jvar = 0; jvar < ei[pg->imtrx]->ielem_dim; jvar++) {
      var = MESH_DISPLACEMENT1 + jvar;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

          d_func[0][var][j] += d_force * d_dist[jvar] * bf[var]->phi[j] * (T_roll - fv->T);
        }
      }
    }
  }

  /* Calculate the residual contribution					     */

  *func += heat_tran_coeff * (T_c - fv->T) + force * (T_roll - fv->T) +
           epsilon * sigma * (pow(T_c, 4.0) - pow(fv->T, 4.0));

  return;
} /* END of routine qrad_surf                                                */

/***********************************************/

void qside_ls(double func[DIM],
              double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
              int no_LS_block_id,
              int LS_block_id,
              double Q_neg_side,
              double Q_pos_side) {
  int j;
  int var = LS;

  if (Current_EB_ptr->Elem_Blk_Id == LS_block_id) {
    /* do nothing */
  } else if (Current_EB_ptr->Elem_Blk_Id == no_LS_block_id) {
    /* Add contributions from level set side of boundary to flux */

    load_lsi_adjmatr(ls->Length_Scale);

    func[0] += lsi->H * (Q_pos_side - Q_neg_side) + Q_neg_side;

    if (af->Assemble_Jacobian) {

      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_func[0][var][j] = lsi->d_H_dF[j] * (Q_pos_side - Q_neg_side);
      }
    }
  } else {
    GOMA_EH(GOMA_ERROR, "Cannot find matching block id in problem.");
  }
}

void shear_to_shell(double cfunc[MDE][DIM],
                    double d_cfunc[MDE][DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                    const int id_side, /* ID of the side of the element */
                    double scale,
                    struct elem_side_bc_struct *elem_side_bc,
                    const int iconnect_ptr,
                    double xi[DIM],
                    const Exo_DB *exo)

{
  int b, i, j, id, ldof;
  int eqn, var;
  int p, q;
  int *n_dof = NULL;
  int nf, el1;
  int n_dofptr[MAX_VARIABLE_TYPES][MDE], dof_map[MDE];
  double detJ;
  double Pi[DIM][DIM];
  STRESS_DEPENDENCE_STRUCT d_Pi;
  double TL; /* shear stress loading */
  double dTL_dX[DIM][MDE];
  double dTL_dv[DIM][MDE];
  double dTL_dP[MDE];

  /* Determine if sh_tens is active on neighboring shell block */
  if (num_shell_blocks != 0)
    nf = num_elem_friends[ei[pg->imtrx]->ielem];
  else
    nf = 0;

  eqn = SHELL_TENSION;

  /* Determine if sh_tens is active on neighboring shell block */
  if (num_shell_blocks != 0)
    nf = num_elem_friends[ei[pg->imtrx]->ielem];
  else
    nf = 0;

  /* If so, set up assembly to include variable shell tension */
  if (nf == 1 && upd->vp[pg->imtrx][SHELL_TENSION]) {
    el1 = elem_friends[ei[pg->imtrx]->ielem][0];
    n_dof = (int *)array_alloc(1, MAX_VARIABLE_TYPES, sizeof(int));
    load_neighbor_var_data(ei[pg->imtrx]->ielem, el1, n_dof, dof_map, n_dofptr, id_side, xi, exo);
  }

  memset(Pi, 0, DIM * DIM * sizeof(double));
  memset(&d_Pi, 0, sizeof(STRESS_DEPENDENCE_STRUCT));

  detJ = 1.0;

  fluid_stress(Pi, &d_Pi);

  TL = 0.0;
  memset(dTL_dv, 0, sizeof(double) * DIM * MDE);
  memset(dTL_dX, 0, sizeof(double) * DIM * MDE);
  memset(dTL_dP, 0, sizeof(double) * MDE);

  for (p = 0; p < ei[pg->imtrx]->ielem_dim; p++) {
    for (q = 0; q < ei[pg->imtrx]->ielem_dim; q++) {
      TL += -scale * fv->snormal[p] * fv->stangent[0][q] * Pi[q][p];

      for (j = 0; j < ei[pg->imtrx]->dof[VELOCITY1]; j++) {
        for (b = 0; b < ei[pg->imtrx]->ielem_dim; b++) {

          dTL_dv[b][j] += -scale * fv->snormal[p] * fv->stangent[0][q] * d_Pi.v[q][p][b][j];
        }
      }

      for (j = 0; j < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; j++) {
        for (b = 0; b < ei[pg->imtrx]->ielem_dim; b++) {
          dTL_dX[b][j] += -scale * fv->dsnormal_dx[p][b][j] * fv->stangent[0][q] * Pi[q][p];
          dTL_dX[b][j] += -scale * fv->snormal[p] * fv->dstangent_dx[0][q][b][j] * Pi[q][p];
          dTL_dX[b][j] += -scale * fv->snormal[p] * fv->stangent[0][q] * d_Pi.X[q][p][b][j];
        }
      }

      for (j = 0; j < ei[pg->imtrx]->dof[PRESSURE]; j++) {
        dTL_dP[j] += -scale * fv->snormal[p] * fv->stangent[0][q] * d_Pi.P[q][p][j];
      }
    }
  }

  for (i = 0; i < (int)elem_side_bc->num_nodes_on_side; i++) {
    id = (int)elem_side_bc->local_elem_node_id[i];
    ldof = ei[pg->imtrx]->ln_to_first_dof[eqn][id];

    if (ldof >= 0) {
      cfunc[ldof][0] = bf[eqn]->phi[ldof] * (TL)*detJ;
    }
  }

  if (af->Assemble_Jacobian) {

    for (i = 0; i < (int)elem_side_bc->num_nodes_on_side; i++) {
      id = (int)elem_side_bc->local_elem_node_id[i];
      ldof = ei[pg->imtrx]->ln_to_first_dof[eqn][id];

      if (ldof >= 0) {

        for (b = 0; b < ei[pg->imtrx]->ielem_dim; b++) {
          var = MESH_DISPLACEMENT1 + b;

          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

              d_cfunc[ldof][0][var][j] = bf[eqn]->phi[ldof] * (dTL_dX[b][j]) * detJ;
            }
          }
        }

        var = PRESSURE;

        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_cfunc[ldof][0][var][j] = bf[eqn]->phi[ldof] * dTL_dP[j] * detJ;
          }
        }

        var = VELOCITY1;
        if (pd->v[pg->imtrx][var]) {

          for (b = 0; b < ei[pg->imtrx]->ielem_dim; b++) {
            var = VELOCITY1 + b;

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_cfunc[ldof][0][var][j] = bf[eqn]->phi[ldof] * dTL_dv[b][j] * detJ;
            }
          }
        }
      }
    }
  }
  return;
}

#ifndef HAVE_HYSTERESIS_WETTING

void apply_hysteresis_wetting_sic(double *func,
                                  double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                                  double delta_t,
                                  double theta,
                                  double *float_data) {

  GOMA_EH(GOMA_ERROR, "Sorry, this model has not been included with this distribution \n");
  return;
}
#else
#include "func_hysteresis_wet.h"
#endif

/************************************************************************************************************/

void fgamma1_deriv_bc(double func[DIM],
                      double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                      const double normalDerivative)

/***********************************************************************
 *
 * fgamma1_deriv_bc:
 *   Function with evalulates the derivative of gamma1 to zero on
 *   a surface.
 */
{
  int dim, i, j, p, b;
  int eqn, var;
  int dofs;
  double grad_gamma1[DIM];
  if (af->Assemble_LSA_Mass_Matrix)
    return;
  BASIS_FUNCTIONS_STRUCT *bf_ptr;
  dim = pd->Num_Dim;

  eqn = R_SHELL_SURF_DIV_V;
  bf_ptr = bf[eqn];

  if (pd->v[pg->imtrx][R_SHELL_SURF_DIV_V]) {

    for (p = 0; p < dim; p++) {
      grad_gamma1[dim] = 0.0;
    }
    dofs = ei[pg->imtrx]->dof[eqn];
    for (p = 0; p < dim; p++) {
      for (i = 0; i < dofs; i++) {
        grad_gamma1[p] += *esp->div_s_v[i] * bf_ptr->grad_phi[i][p];
      }
    }
    func[0] = normalDerivative;
    for (p = 0; p < dim; p++) {
      func[0] += fv->snormal[p] * grad_gamma1[p];
    }

    if (af->Assemble_Jacobian) {
      var = R_SHELL_SURF_DIV_V;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (p = 0; p < VIM; p++) {
            d_func[0][var][j] += fv->snormal[p] * bf_ptr->grad_phi[j][p];
          }
        }
      }

      if (pd->v[pg->imtrx][R_MESH1]) {
        for (b = 0; b < dim; b++) {
          var = R_MESH1 + b;
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            for (i = 0; i < dofs; i++) {
              for (p = 0; p < VIM; p++) {
                d_func[0][var][j] +=
                    fv->snormal[p] * (*esp->div_s_v[i]) * bf_ptr->d_grad_phi_dmesh[i][p][b][j];
                d_func[0][var][j] += fv->dsnormal_dx[p][b][j] * grad_gamma1[p];
              }
            }
          }
        }
      }
    }
  }
  return;
} /* END of routine fgamma1_deriv_bc                                         */

/************************************************************************************************************/

void fgamma2_deriv_bc(double func[DIM],
                      double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                      const double normalDerivative)

/***********************************************************************
 *
 * fgamma2_deriv_bc:
 *   Function with evalulates the derivative of gamma2 to zero on
 *   a surface of a shell, which is actually an edge.
 */
{
  int dim, i, j, p, b;
  int eqn, var;
  int dofs;
  double grad_gamma2[DIM];
  if (af->Assemble_LSA_Mass_Matrix)
    return;
  BASIS_FUNCTIONS_STRUCT *bf_ptr;
  dim = pd->Num_Dim;

  eqn = R_SHELL_SURF_CURV;
  bf_ptr = bf[eqn];

  if (pd->v[pg->imtrx][R_SHELL_SURF_CURV]) {

    for (p = 0; p < dim; p++) {
      grad_gamma2[dim] = 0.0;
    }
    dofs = ei[pg->imtrx]->dof[eqn];
    for (p = 0; p < dim; p++) {
      for (i = 0; i < dofs; i++) {
        grad_gamma2[p] += *esp->curv[i] * bf_ptr->grad_phi[i][p];
      }
    }
    func[0] = normalDerivative;
    for (p = 0; p < dim; p++) {
      func[0] += fv->snormal[p] * grad_gamma2[p];
    }

    if (af->Assemble_Jacobian) {
      var = R_SHELL_SURF_CURV;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (p = 0; p < VIM; p++) {
            d_func[0][var][j] += fv->snormal[p] * bf_ptr->grad_phi[j][p];
          }
        }
      }

      if (pd->v[pg->imtrx][R_MESH1]) {
        for (b = 0; b < dim; b++) {
          var = R_MESH1 + b;
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            for (i = 0; i < dofs; i++) {
              for (p = 0; p < dim; p++) {
                // For shell equations snormal is in jShell numbering
                // for shell equations d_grad_phi_dmesh is
                d_func[0][var][j] +=
                    fv->snormal[p] * (*esp->curv[i]) * bf_ptr->d_grad_phi_dmesh[i][p][b][j];
                d_func[0][var][j] += fv->dsnormal_dx[p][b][j] * grad_gamma2[p];
              }
            }
          }
        }
      }
    }
  }
  return;
} /* END of routine fgamma2_deriv_bc                                         */

/************************************************************************************************************/

void dvzdr_zero_deriv_bc(double func[DIM],
                         double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                         const double nwall[DIM],
                         const double normalDerivative)

/***********************************************************************
 *
 *  dvzdr_zero_deriv_bc
 *   Function with evalulates the derivative of dvaxial to zero on
 *   the r = 0 surface
 */
{
  int dim, i, j, p, b;
  int eqn, var;
  int dofs;
  double grad_vz;
  if (af->Assemble_LSA_Mass_Matrix)
    return;
  BASIS_FUNCTIONS_STRUCT *bf_ptr;
  dim = pd->Num_Dim;

  eqn = R_MOMENTUM1;
  bf_ptr = bf[eqn];

  if (pd->v[pg->imtrx][R_MOMENTUM1]) {

    grad_vz = 0.0;

    dofs = ei[pg->imtrx]->dof[eqn];
    // grad_v[dir][xdir]
    for (p = 0; p < dim; p++) {
      for (i = 0; i < dofs; i++) {
        grad_vz += *esp->v[0][i] * bf_ptr->grad_phi[i][p] * nwall[p];
      }
    }
    func[0] = grad_vz - normalDerivative;

    if (af->Assemble_Jacobian) {
      var = R_MOMENTUM1;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (p = 0; p < VIM; p++) {
            d_func[0][var][j] += bf_ptr->grad_phi[j][p] * nwall[p];
          }
        }
      }

      if (pd->v[pg->imtrx][R_MESH1]) {
        for (b = 0; b < dim; b++) {
          var = R_MESH1 + b;
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            for (i = 0; i < dofs; i++) {
              for (p = 0; p < VIM; p++) {
                d_func[0][var][j] +=
                    *esp->v[0][i] * bf_ptr->d_grad_phi_dmesh[i][p][b][j] * nwall[p];
              }
            }
          }
        }
      }
    }
  }
  return;
} /* END of routine fgamma2_deriv_bc                                         */
/*****************************************************************************/

void light_transmission(double func[DIM],
                        double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                        double time,
                        const int bc_type,
                        const double bdy_refindex,
                        UNUSED const double bdy_absorption,
                        const double bdy_incident,
                        UNUSED const int blk_id)
/******************************************************************************
 *  Function which calculates the surface integral for the "plane transmission"
 *	acoustic wave equation
 ******************************************************************************/
{
  /* Local variables */

  int j, b, w, dim;
  int eqn = 0, eqn_alt = 0, var;

  /*
   *    Radiative transfer equation variables - connect to input file someday
   */
  double svect[3] = {0., -1., 0.};
  double mucos = 1.0;
  double mucos_tran, mu_crit;
  double Grefl, Rrefl, Xrefl, Yrefl;
  CONDUCTIVITY_DEPENDENCE_STRUCT d_alpha_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_alpha = &d_alpha_struct;
  double refindex, refratio, direction; /* Refractive Index */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_n_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_n = &d_n_struct;

  double mucos_tran_dn, Grefl_dn, Rrefl_dn, Xrefl_dn, Yrefl_dn;

  /***************************** EXECUTION BEGINS *******************************/
  dim = pd->Num_Dim;
  if (af->Assemble_LSA_Mass_Matrix)
    return;

  light_absorption(d_alpha, time);
  refindex = refractive_index(d_n, time);

  direction = svect[0] * fv->snormal[0] + svect[1] * fv->snormal[1] + svect[2] * fv->snormal[2];
  refratio = bdy_refindex / refindex;
  if (refratio >= 1.0) {
    mu_crit = sqrt(1.0 - 1.0 / SQUARE(refratio));
  } else {
    mu_crit = 0.0;
  }

  if (direction < 0.0) {
    mucos_tran = sqrt(1.0 - SQUARE(bdy_refindex / refindex) * (1.0 - SQUARE(mucos)));
    mucos_tran_dn = SQUARE(bdy_refindex) / CUBE(refindex) * (1.0 - SQUARE(mucos)) / mucos_tran;
  } else {
    mucos_tran = sqrt(1.0 - SQUARE(refindex / bdy_refindex) * (1.0 - SQUARE(mucos)));
    mucos_tran_dn = -refindex / SQUARE(bdy_refindex) * (1.0 - SQUARE(mucos)) / mucos_tran;
  }
  Grefl = 0.5 * (SQUARE((refindex * mucos - bdy_refindex * mucos_tran) /
                        (refindex * mucos + bdy_refindex * mucos_tran)) +
                 SQUARE((bdy_refindex * mucos - refindex * mucos_tran) /
                        (bdy_refindex * mucos + refindex * mucos_tran)));
  Grefl_dn =
      2 * bdy_refindex * mucos *
      ((mucos_tran - refindex * mucos_tran_dn) * (refindex * mucos - bdy_refindex * mucos_tran) /
           CUBE(refindex * mucos + bdy_refindex * mucos_tran) -
       (mucos_tran + refindex * mucos_tran_dn) * (bdy_refindex * mucos - refindex * mucos_tran) /
           CUBE(bdy_refindex * mucos + refindex * mucos_tran));
  if (mucos >= mu_crit) {
    Rrefl = Grefl;
    Rrefl_dn = Grefl_dn;
  } else {
    Rrefl = 1.0;
    Rrefl_dn = 0.;
  }
  if (bdy_refindex <= refindex) {
    Xrefl = Grefl;
    Xrefl_dn = Grefl_dn;
  } else {
    Xrefl = Rrefl;
    Xrefl_dn = Rrefl_dn;
  }

  Yrefl = SQUARE(bdy_refindex / refindex) * (1.0 - Xrefl);
  Yrefl_dn = SQUARE(bdy_refindex / refindex) * (-Xrefl_dn - 2 * (1.0 - Xrefl) / refindex);

#if 0
fprintf(stderr,"refl n nbdy X Y mu mut dir %g %g %g %g %g %g %g\n",refindex,bdy_refindex,Xrefl,Yrefl,mucos,mucos_tran, direction);
#endif

  if (bc_type == LIGHTP_TRANS_BC) {
    eqn = LIGHT_INTP;
    eqn_alt = LIGHT_INTM;
    *func = fv->poynt[0] - Xrefl * fv->poynt[1] - Yrefl * bdy_incident;
  } else if (bc_type == LIGHTM_TRANS_BC) {
    eqn = LIGHT_INTM;
    eqn_alt = LIGHT_INTP;
    *func = fv->poynt[1] - Xrefl * fv->poynt[0] - Yrefl * bdy_incident;
  } else {
    GOMA_EH(GOMA_ERROR, "invalid light transmission bc\n");
  }

  if (af->Assemble_Jacobian) {
    if (pd->v[pg->imtrx][eqn]) {
      for (j = 0; j < ei[pg->imtrx]->dof[eqn]; j++) {
        d_func[0][eqn][j] = bf[eqn]->phi[j];
      }
    }
    if (pd->v[pg->imtrx][eqn_alt]) {
      for (j = 0; j < ei[pg->imtrx]->dof[eqn_alt]; j++) {
        d_func[0][eqn_alt][j] = -Xrefl * bf[eqn_alt]->phi[j];
      }
    }
    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_func[0][var][j] =
            d_n->T[j] * (-Xrefl_dn * fv->poynt[LIGHT_INTM - eqn] - Yrefl_dn * bdy_incident);
      }
    }
    for (b = 0; b < dim; b++) {
      var = MESH_DISPLACEMENT1 + b;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_func[0][var][j] =
              d_n->X[b][j] * (-Xrefl_dn * fv->poynt[LIGHT_INTM - eqn] - Yrefl_dn * bdy_incident);
        }
      }
    }
    var = MASS_FRACTION;
    if (pd->e[pg->imtrx][eqn] && pd->v[pg->imtrx][var]) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_func[0][MAX_PROB_VAR + w][j] =
              d_n->C[w][j] * (-Xrefl_dn * fv->poynt[LIGHT_INTM - eqn] - Yrefl_dn * bdy_incident);
        }
      }
    }
  }
  return;
} /* END of routine acoustic_plane_transmission                             */
/****************************************************************************/
/*****************************************************************************/
void qside_light_jump(double func[DIM],
                      double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                      double time,
                      const int bc_type,
                      int ID_mat_1, /* block ID material 1   */
                      int ID_mat_2  /* block ID material 2   */
                      )
/******************************************************************************
 *
 *  Function which applies light intensity jump at a side set interface
 *  Author: Robert Secor (4/29/2015)
 *
 ******************************************************************************/

{

  /* Local variables */

  int sign_int = 0;
  int j, b, w, dim;
  int eqn = 0, eqn_alt = 0, var;

  /*
   *    Radiative transfer equation variables - connect to input file someday
   */
  double svect[3] = {0., -1., 0.};
  double mucos = 1.0;
  double mucos_tran, mu_crit;
  double Grefl, Rrefl, Xrefl, Yrefl;
  CONDUCTIVITY_DEPENDENCE_STRUCT d_alpha_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_alpha = &d_alpha_struct;
  double refindex, refratio, direction; /* Refractive Index */
  double other_refindex;                /* Refractive Index */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_n_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_n = &d_n_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_n_other = &d_n_struct;

  double mucos_tran_dn, Grefl_dn, Rrefl_dn, Xrefl_dn, Yrefl_dn;
  MATRL_PROP_STRUCT *mp_2;
  MATRL_PROP_STRUCT *mp_1;

  /***************************** EXECUTION BEGINS *******************************/
  dim = pd->Num_Dim;
  if (af->Assemble_LSA_Mass_Matrix)
    return;
  if (Current_EB_ptr->Elem_Blk_Id == ID_mat_1) {
    sign_int = -1;
  } else if (Current_EB_ptr->Elem_Blk_Id == ID_mat_2) {
    sign_int = 1;
  } else {
    GOMA_EH(GOMA_ERROR, "LIGHT_JUMP has incorrect material ids");
  }

  mp_1 = mp_glob[Current_EB_ptr->Elem_Blk_Id - 1];
  mp = mp_1;
  load_matrl_statevector(mp_1);
  load_properties(mp_1, time);
  light_absorption(d_alpha, time);
  refindex = refractive_index(d_n, time);
  /*other mp */
  if (sign_int == -1) {
    mp_2 = mp_glob[ID_mat_2 - 1];
  } else {
    mp_2 = mp_glob[ID_mat_1 - 1];
  }
  if (mp_2 != mp_1) {
    mp = mp_2;
    load_matrl_statevector(mp_2);
    load_properties(mp_2, time);
  }
  other_refindex = refractive_index(d_n_other, time);
  fprintf(stderr, "refractive index %g %g %d %d %d %d\n", refindex, other_refindex, bc_type,
          ID_mat_1, ID_mat_2, sign_int);
  /*reset material*/
  mp = mp_1;
  load_matrl_statevector(mp_1);
  load_properties(mp_1, time);

  direction = svect[0] * fv->snormal[0] + svect[1] * fv->snormal[1] + svect[2] * fv->snormal[2];
  refratio = other_refindex / refindex;
  if (refratio >= 1.0) {
    mu_crit = sqrt(1.0 - 1.0 / SQUARE(refratio));
  } else {
    mu_crit = 0.0;
  }

  if (direction < 0.0) {
    mucos_tran = sqrt(1.0 - SQUARE(other_refindex / refindex) * (1.0 - SQUARE(mucos)));
    mucos_tran_dn = SQUARE(other_refindex) / CUBE(refindex) * (1.0 - SQUARE(mucos)) / mucos_tran;
  } else {
    mucos_tran = sqrt(1.0 - SQUARE(refindex / other_refindex) * (1.0 - SQUARE(mucos)));
    mucos_tran_dn = -refindex / SQUARE(other_refindex) * (1.0 - SQUARE(mucos)) / mucos_tran;
  }
  Grefl = 0.5 * (SQUARE((refindex * mucos - other_refindex * mucos_tran) /
                        (refindex * mucos + other_refindex * mucos_tran)) +
                 SQUARE((other_refindex * mucos - refindex * mucos_tran) /
                        (other_refindex * mucos + refindex * mucos_tran)));
  Grefl_dn =
      2 * other_refindex * mucos *
      ((mucos_tran - refindex * mucos_tran_dn) * (refindex * mucos - other_refindex * mucos_tran) /
           CUBE(refindex * mucos + other_refindex * mucos_tran) -
       (mucos_tran + refindex * mucos_tran_dn) * (other_refindex * mucos - refindex * mucos_tran) /
           CUBE(other_refindex * mucos + refindex * mucos_tran));
  if (mucos >= mu_crit) {
    Rrefl = Grefl;
    Rrefl_dn = Grefl_dn;
  } else {
    Rrefl = 1.0;
    Rrefl_dn = 0.;
  }
  if (other_refindex <= refindex) {
    Xrefl = Grefl;
    Xrefl_dn = Grefl_dn;
  } else {
    Xrefl = Rrefl;
    Xrefl_dn = Rrefl_dn;
  }

  Yrefl = SQUARE(other_refindex / refindex) * (1.0 - Xrefl);
  Yrefl_dn = SQUARE(other_refindex / refindex) * (-Xrefl_dn - 2 * (1.0 - Xrefl) / refindex);

  fprintf(stderr, "mu mut Xrefl Yrefl %g %g %g %g \n", mucos, mucos_tran, Xrefl, Yrefl);
  fprintf(stderr, "coords %g %g  \n", fv->x[0], fv->x[1]);
  fprintf(stderr, "dofs elem %d %d %d\n", ei[pg->imtrx]->dof[LIGHT_INTP],
          ei[pg->imtrx]->dof[LIGHT_INTM], ei[pg->imtrx]->ielem);
  for (j = 0; j < ei[pg->imtrx]->dof[LIGHT_INTP]; j++) {
    fprintf(stderr, "vars %d %g %g %g\n", j, *esp->poynt[0][j], *esp->poynt[1][j],
            bf[LIGHT_INTP]->phi[j]);
  }
  fprintf(stderr, "light intensity %g %g %g \n", fv->poynt[0], fv->poynt[1], *func);
  /* Calculate the residual contribution					     */

  if (bc_type == LIGHTP_JUMP_BC || bc_type == LIGHTP_JUMP_2_BC) {
    eqn = LIGHT_INTP;
    eqn_alt = LIGHT_INTM;
    if (pd->v[pg->imtrx][eqn]) {
      if (sign_int == -1) {
        *func += fv->poynt[0] - Xrefl * fv->poynt[1];
      } else {
        *func -= Yrefl * fv->poynt[0];
      }
    }
  } else if (bc_type == LIGHTM_JUMP_BC || bc_type == LIGHTM_JUMP_2_BC) {
    eqn = LIGHT_INTM;
    eqn_alt = LIGHT_INTP;
    if (pd->v[pg->imtrx][eqn]) {
      if (sign_int == -1) {
        *func += fv->poynt[1] - Xrefl * fv->poynt[0];
      } else {
        *func -= Yrefl * fv->poynt[1];
      }
    }
  } else {
    GOMA_EH(GOMA_ERROR, "invalid light transmission bc\n");
  }
  fprintf(stderr, "light intensity %g %g %g \n", fv->poynt[0], fv->poynt[1], *func);

  if (af->Assemble_Jacobian) {
    if (pd->v[pg->imtrx][eqn]) {
      if (sign_int == -1) {
        for (j = 0; j < ei[pg->imtrx]->dof[eqn]; j++) {
          d_func[0][eqn][j] += bf[eqn]->phi[j];
        }
      } else {
        for (j = 0; j < ei[pg->imtrx]->dof[eqn]; j++) {
          d_func[0][eqn][j] += -Yrefl * bf[eqn]->phi[j];
        }
      }
    }
    if (pd->v[pg->imtrx][eqn_alt] && sign_int == -1) {
      for (j = 0; j < ei[pg->imtrx]->dof[eqn_alt]; j++) {
        d_func[0][eqn_alt][j] += -Xrefl * bf[eqn_alt]->phi[j];
      }
    }
    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_func[0][var][j] = d_n->T[j] * (-Xrefl_dn * fv->poynt[LIGHT_INTM - eqn] - Yrefl_dn);
      }
    }
    for (b = 0; b < dim; b++) {
      var = MESH_DISPLACEMENT1 + b;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_func[0][var][j] = d_n->X[b][j] * (-Xrefl_dn * fv->poynt[LIGHT_INTM - eqn] - Yrefl_dn);
        }
      }
    }
    var = MASS_FRACTION;
    if (pd->e[pg->imtrx][eqn] && pd->v[pg->imtrx][var]) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_func[0][MAX_PROB_VAR + w][j] =
              d_n->C[w][j] * (-Xrefl_dn * fv->poynt[LIGHT_INTM - eqn] - Yrefl_dn);
        }
      }
    }
  }
}
/****************************************************************************/
void fvelo_slip_ls_heaviside(double func[DIM],
                             double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                             double width,
                             double beta_negative,
                             double beta_positive,
                             const double vsx,
                             const double vsy,
                             const double vsz) {
  int j, var, jvar, p;
  double phi_j, vs[MAX_PDIM];
  double beta, betainv;
  double d_beta_dF[MDE];
  /************************* EXECUTION BEGINS *******************************/

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  load_lsi(width);
  if (af->Assemble_Jacobian) {
    load_lsi_derivs();
    memset(d_beta_dF, 0, MDE * sizeof(double));
  }

  level_set_property(beta_negative, beta_positive, width, &beta, d_beta_dF);
  // betainv = mu/beta;
  betainv = 1 / beta;

  vs[0] = vsx;
  vs[1] = vsy;
  vs[2] = vsz;

  if (af->Assemble_Jacobian) {
    for (jvar = 0; jvar < pd->Num_Dim; jvar++) {
      var = VELOCITY1 + jvar;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          d_func[jvar][var][j] += (-betainv) * (phi_j);
        }
      }
    }

    var = LS;
    if (pd->v[pg->imtrx][var]) {
      for (p = 0; p < pd->Num_Dim; p++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];

          d_func[p][var][j] += (d_beta_dF[j] / beta / beta) * (fv->v[p] - vs[p]);
          // d_func[p][var][j] += (d_beta_dF[j])*( fv->v[p] - vs[p] );
        }
      }
    }
  } /* end of of Assemble Jacobian		*/

  /* Calculate the residual contribution	*/
  for (p = 0; p < pd->Num_Dim; p++) {
    func[p] += (-betainv) * (fv->v[p] - vs[p]);
  }
}

void fvelo_slip_ls_oriented(double func[MAX_PDIM],
                            double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                            double width,
                            double beta_negative,
                            double beta_positive,
                            double gamma_negative,
                            double gamma_positive,
                            const double vsx, /* velocity components of solid  */
                            const double vsy, /* surface on which slip condition   */
                            const double vsz) /* is applied           */
{
  int j, var, jvar, p;
  double phi_j, vs[MAX_PDIM];
  double beta, betainv;
  double d_beta_dF[MDE];
  double gamma, gammainv;
  double d_gamma_dF[MDE];
  /************************* EXECUTION BEGINS *******************************/

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  load_lsi(width);
  if (af->Assemble_Jacobian) {
    load_lsi_derivs();
    memset(d_beta_dF, 0, MDE * sizeof(double));
  }
  level_set_property(beta_negative, beta_positive, width, &beta, d_beta_dF);
  level_set_property(gamma_negative, gamma_positive, width, &gamma, d_gamma_dF);
  betainv = 1 / beta;
  gammainv = 1 / gamma;

  vs[0] = vsx;
  vs[1] = vsy;
  vs[2] = vsz;

  if (af->Assemble_Jacobian) {
    for (jvar = 0; jvar < pd->Num_Dim; jvar++) {
      var = VELOCITY1 + jvar;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          d_func[jvar][var][j] += (-betainv) * (phi_j);

          for (int q = 0; q < pd->Num_Dim; q++) {
            d_func[jvar][var][j] +=
                -((gammainv - betainv) * (fv->snormal[jvar] * fv->snormal[q]) * phi_j);
          }
        }
      }

      var = LS;
      if (pd->v[pg->imtrx][var]) {
        for (p = 0; p < pd->Num_Dim; p++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_func[p][var][j] += (d_beta_dF[j] * betainv * betainv) * (fv->v[p] - vs[p]);
            for (int q = 0; q < pd->Num_Dim; q++) {
              d_func[p][var][j] +=
                  (d_gamma_dF[j] * gammainv * gammainv - d_beta_dF[j] * betainv * betainv) *
                  (fv->v[p] - vs[p]);
            }
          }
        }
      }
    }
  } /* end of of Assemble Jacobian		*/

  /* Calculate the residual contribution	*/
  for (p = 0; p < pd->Num_Dim; p++) {
    double p_vel = (fv->v[p] - vs[p]);
    func[p] += -(betainv)*p_vel;
    for (int q = 0; q < pd->Num_Dim; q++) {
      func[p] += -((gammainv - betainv) * (fv->snormal[p] * fv->snormal[q]) * p_vel);
    }
  }
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void ls_wall_angle_bc(double func[DIM],
                      double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                      const double angle) /* angle in radians */
{
  int j, kdir, var, p;
  const double cos_angle = cos(angle);

  /* Calculate the residual contribution	*/
  func[0] = -cos_angle;
  for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
    func[0] += fv->grad_F[kdir] * fv->snormal[kdir];
  }

  if (af->Assemble_Jacobian) {

    for (kdir = 0; kdir < pd->Num_Dim; kdir++) {

      for (p = 0; p < pd->Num_Dim; p++) {
        var = MESH_DISPLACEMENT1 + p;
        if (pd->gv[var]) {
          for (j = 0; j < ei[pd->mi[var]]->dof[var]; j++) {
            d_func[0][var][j] += fv->grad_F[kdir] * fv->dsnormal_dx[kdir][p][j] +
                                 fv->d_grad_F_dmesh[kdir][p][j] * fv->snormal[kdir];
          }
        }
      }

      var = FILL;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pd->mi[var]]->dof[var]; j++) {
          d_func[0][var][j] += bf[var]->grad_phi[j][kdir] * fv->snormal[kdir];
        }
      }
    } /* for: kdir */
  }   /* end of if Assemble_Jacobian */
}
