/* Copyright 2019 Andrew Cochrane
 * This file released under MIT license
 * */

#include "shell_tfmp_util.h"

#include <string.h>
#include <math.h>

#include "std.h"
#include "rf_fem_const.h"
#include "el_geom.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_eh.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "mm_shell_util.h"
#include "mm_std_models_shell.h"
#include "el_elm.h"
#include "shell_tfmp_struct.h"

enum clipping_kind my_clipping_kind = restorative;

void load_tfmp_viscosity_model
( double *mu_l,
  double *mu_g
) {
  switch(mp->tfmp_viscosity_model){
    case CONSTANT:
      *mu_g = mp->tfmp_viscosity_const[0];
      *mu_l = mp->tfmp_viscosity_const[1];
      break;
    default:
      WH(-1, "There is no tfmp viscosity");
      return;
      break;
  }
}

void load_gas_density_model
( double *Patm,
  double *rho_g,
  double *drho_g_dP
) {
 //temporary variables
 double Mg, R, T;
 switch(mp->tfmp_density_model) {
  case IDEAL_GAS:
    Mg =   mp->tfmp_density_const[0]; // grams per mole
    R =    mp->tfmp_density_const[1];
    T =    mp->tfmp_density_const[2];
    *Patm = mp->tfmp_density_const[3];
    *rho_g = fv->tfmp_pres*Mg/R/T;
    *drho_g_dP = Mg/R/T;
    break;

  case CONSTANT:
    *Patm = mp->tfmp_density_const[3];
    *rho_g = mp->tfmp_density_const[0]; // grams per cubic centimeter
    *drho_g_dP = 0.0;
    break;

  default:
    *Patm = 0.0;
    *rho_g = 1.0; // grams per cubic centimeter
    *drho_g_dP = 0.0;
    break;
  }
}

void load_molecular_diffusion_model
( double S,           // saturation
  double *D,          // Diffusion constant L^2/T
  double *Krd,        // molecular diffusion correction factor
  double *dKrd_dS     // partial wrt S
  ) {
  // temporary variables
  double Scd, betad, md, cd;

  switch (mp->tfmp_diff_model) {
    case CONSTANT:
      *D = mp->tfmp_diff_const[0];
      *Krd = 1.0;
      *dKrd_dS = 0.0;
    break;
    case PIECEWISE:
      *D = mp->tfmp_diff_const[0];
      // diffusion transition
      Scd = mp->tfmp_diff_const[1];
      betad = mp->tfmp_diff_const[2];
      md = 1.f/2.f/betad;
      cd = -md*(Scd - betad);

      if ( S < Scd - betad) {
        *Krd = 0.0;
        *dKrd_dS = 0.0;
      } else {
        *Krd = md*S + cd;
        *dKrd_dS = md;
      }
    break;

    default:
      *D = 0.0;
      *Krd = 0.0;
      *dKrd_dS = 0.0;
    break;
  }
}

void load_relative_permeability_model
( double S,         // saturation
  double *Krl,      // liquid rel perm
  double *dKrl_dS,  // partial wrt S
  double *Krg,      // gas rel perm
  double *dKrg_dS   // partial wrt S
) {
// temporary variables
double Scg, alphag, Scl, alphal, mg, cg, ml, cl;
// set model initiators
switch (mp->tfmp_rel_perm_model) {
  case PIECEWISE:
    Scg = mp->tfmp_rel_perm_const[0];
    alphag = mp->tfmp_rel_perm_const[1];
    Scl = mp->tfmp_rel_perm_const[2];
    alphal = mp->tfmp_rel_perm_const[3];
    break;
  case LEVER:
    Scg = 0.5;
    alphag = 0.5;
    Scl = 0.5;
    alphal = 0.5;
    break;
  case SATURATION:
    *Krl = S;
    *dKrl_dS = 1.0;
    *Krg = 1.0 - S;
    *dKrg_dS = -1.0;
    return;
    break;
  default:
    WH(-1, "relative permeability model not set");
    return;
    break;
  }

  //compute slopes and offsets
  mg = -1./2./alphag;
  cg = -mg*(Scg + alphag);
  ml = 1./2./alphal;
  cl = -ml*(Scl - alphal);

  // calculate perms based on saturation
  // gas
  if ( S <= Scg - alphag) {
    *Krg = 1.0;
    *dKrg_dS = 0.0;
  } else if ( S > Scg - alphag && S < Scg + alphag ) {
    *Krg = mg*S + cg;
    *dKrg_dS = mg;
  } else {
    *Krg = 0.0;
    *dKrg_dS = 0.0;
  }

  // liquid
  if ( S <= Scl - alphal) {
    *Krl = 0.0;
    *dKrl_dS = 0.0;
  } else if ( S > Scl - alphal && S < Scl + alphal ) {
    *Krl = ml*S + cl;
    *dKrl_dS = ml;

  } else { // S > 1.0
    *Krl = 1.0 ;
    *dKrl_dS = 0.0;
  }
}

void load_gas_dissolution_model(
  double h,
  double Patm,
  double *J,
  double *dJ_dP,
  double *dJ_dS,
  double *dJ_dh
) {
  // temporary variables
  double Dgl, Hgls, Mg,trans_diss,sqrtPI, S, bg, dbg_dP, fS, dfS_dS,
         lambda, Vd, L, dL_dh, pi;

  pi = 3.1415926;

  S = fv->tfmp_sat;

  switch(mp->tfmp_drop_lattice_model) {
  case TFMP_SQUARE:
    lambda = mp->tfmp_drop_lattice_const[0]; // cm
    Vd = mp->tfmp_drop_lattice_const[1];// cm^3
    break;
  default:
    lambda = Vd = 1.0;
    if (mp->tfmp_dissolution_model != NO_MODEL) {
      WH(-1, "The lattice model is not set and the dissolution model is not NO_MODEL");
      return;
    }
    break;
  }

  switch(mp->tfmp_dissolution_model) {
  case TFMP_SQUARE:
    lambda = mp->tfmp_drop_lattice_const[0]; // cm
    Vd = mp->tfmp_drop_lattice_const[1];// cm^3
    // Diffusion Coefficient of gas species in liquid solvent
    Dgl = mp->tfmp_dissolution_const[0];
    // Henry's Law Constant
    Hgls = mp->tfmp_dissolution_const[1];
    // Molecular weight of gas
    Mg = mp->tfmp_dissolution_const[2];

    trans_diss = pi/4.0;
    sqrtPI = sqrt(pi);

    if (S < 1.0 && S > 0.0 ) {

      // dependence on saturation
      if (S > trans_diss) {        
        fS = 2.0*sqrtPI/lambda*sqrt(S);
        dfS_dS = sqrtPI/lambda/sqrt(S);
        L = sqrt(Vd/pi/h);
        dL_dh = -sqrt(Vd/pi/h*h*h)/2.0;
      } else {

        fS = 4.0/lambda*sqrt(1.0-S);
        dfS_dS = -4.0/lambda/2.0/sqrt(1.0-S);

        if (Vd/lambda/lambda/h >= 1.0) { // Length is maxed for lattice defined by
          // Vd and lambda
          L = lambda/sqrt(2.0);
          dL_dh = 0.0;
        } else {
          L = lambda/sqrt(2.0)*(1.0 - sqrt((1.0-Vd/lambda/lambda/h)/2.0));
          dL_dh = lambda/2.0/sqrt(2.0*(1.0-Vd/lambda/lambda/h))*Vd/lambda/lambda/h/h;

        }
      }

      //Concentration of gas at interface moles/cm^3
      bg = (fv->tfmp_pres - Patm)*Hgls;
      dbg_dP = Hgls;

      // Dissolution Term
      *J = h*Dgl*Mg*bg*fS/L;
      *dJ_dP = h*Dgl*Mg*dbg_dP*fS/L;
      *dJ_dS = h*Dgl*Mg*bg*dfS_dS/L;
      *dJ_dh = fS*Dgl*Mg*bg*(-h/L/L*dL_dh + 1.0/L);


    } else {
      fS = 0.0;
      dfS_dS = 0.0;

      L = 1.0;
      dL_dh = 0.0;

      bg = 0.0;
      dbg_dP= 0.0;

      *J = 0.0;
      *dJ_dP = 0.0;
      *dJ_dS = 0.0;
      *dJ_dh = 0.0;
    }
    break;
  default:
    fS = 0.0;
    dfS_dS = 0.0;

    L = 1.0;
    dL_dh = 0.0;

    bg = 0.0;
    dbg_dP= 0.0;

    *J = 0.0;
    *dJ_dP = 0.0;
    *dJ_dS = 0.0;
    *dJ_dh = 0.0;

  }
}

void load_displacement_coupling_model
( double tt,
  double delta_t,
  double *h,
  double *dh_dtime,
  double *gradII_h,
  double dh_dmesh[][MDE],
  double dh_dnormal[][MDE],
  double d2h_dtime_dmesh[][MDE],
  double d2h_dtime_dnormal[][MDE],
  double d_gradIIh_dmesh[][DIM][MDE],
  double d_gradIIh_dnormal[][DIM][MDE],
  int *n_dof,
  int *dof_map
) {
  switch (mp->ehl_gap_model){
    case GM_NDOTD:
      h0_minus_ndotd(tt, delta_t, h, dh_dtime,
                     gradII_h, dh_dmesh, dh_dnormal,
                     d2h_dtime_dmesh, d2h_dtime_dnormal,
                     d_gradIIh_dmesh, d_gradIIh_dnormal,
                     n_dof, dof_map);
      break;
    case GM_RADIAL:
      rmesh_minus_rroller(tt, delta_t, h, dh_dtime,
                          gradII_h, dh_dmesh, dh_dnormal,
                          d2h_dtime_dmesh, d2h_dtime_dnormal,
                          d_gradIIh_dmesh, d_gradIIh_dnormal,
                          n_dof, dof_map);
      break;
  }


}

void h0_minus_ndotd (
  double tt,
  double delta_t,
  double *h,
  double *dh_dtime,
  double *gradII_h,
  double dh_dmesh[][MDE],
  double dh_dnormal[][MDE],
  double d2h_dtime_dmesh[][MDE],
  double d2h_dtime_dnormal[][MDE],
  double d_gradIIh_dmesh[][DIM][MDE],
  double d_gradIIh_dnormal[][DIM][MDE],
  int *n_dof,
  int *dof_map
) {
  // temporary variables
  int i, k, l, m, eqn;

  // assume the web is always below the gap

  double h_sign = -1.0;
  double phi_i, grad_phi_i[DIM], gradII_phi_i[DIM]; // Basis funcitons (i)
  double d_gradII_phi_i_dmesh[DIM][DIM][MDE];


  if (pd->TimeIntegration == STEADY) {
    // don't divide by 0
        delta_t = 1.0;
  }


  // lots of stuff needed for gradh when h += d_dot_n
  double grad_n[DIM][DIM], gradII_n[DIM][DIM];
  double dgrad_n_dmesh[DIM][DIM][DIM][MDE], dgradII_n_dmesh[DIM][DIM][DIM][MDE];
  double grad_disp[DIM][DIM], gradII_disp[DIM][DIM];
  double dgrad_disp_dmesh[DIM][DIM][DIM][MDE], dgradII_disp_dmesh[DIM][DIM][DIM][MDE];

  double normal[DIM];
  double grad_normal[DIM][DIM];
  double d_grad_n_dmesh[DIM][DIM][DIM][MDE];

  // need to make ds_dcsi = det_J;
  double det_J;
  double d_det_J_dmeshkj[DIM][MDE];

  memset(d_det_J_dmeshkj, 0.0, sizeof(double)*DIM*MDE);

  memset(grad_normal, 0.0, sizeof(double)*DIM*DIM);
  memset(d_grad_n_dmesh, 0.0, sizeof(double)*DIM*DIM*DIM*MDE);

  // construct the normal - there are probably better names for these enums
  switch(mp->ehl_normal_method) {
    case NCM_PRIMITIVE_XY:
    for (k=0; k<pd->Num_Dim; k++) { // vector element k
      normal[k] = fv->n[k];
      for (l=0; l<pd->Num_Dim; l++) { // derivative direction l
        grad_normal[k][l] = fv->grad_n[k][l];
        for (m=0; m<pd->Num_Dim; m++){
          for (i=0; i<ei[pg->imtrx]->dof[SHELL_NORMAL1]; i++) {
            d_grad_n_dmesh[k][l][m][i] = fv->d_grad_n_dmesh[k][l][m][i];
          }
        }
      }
    }
    break;

    case NCM_PRIMITIVE_S_ROLLER:
      // this has not been implemented
    case NCM_PRIMITIVE_S_WEB:
      detJ_2d_bar(&det_J, d_det_J_dmeshkj);
      //dn_dcsi
      double dn_dcsi[DIM];

      memset(dn_dcsi, 0.0, sizeof(double)*DIM);

      for (int i=0; i<ei[pg->imtrx]->dof[SHELL_NORMAL1]; i++){
          dn_dcsi[0] += *esp->n[0][i]*bf[SHELL_NORMAL1]->dphidxi[i][0];
          dn_dcsi[1] += *esp->n[1][i]*bf[SHELL_NORMAL2]->dphidxi[i][0];
      }

      for (int k=0; k<DIM; k++) { // vector element k
        normal[k] = fv->n[k];
        grad_normal[k][0] = dn_dcsi[k]/det_J;
        for (int m=0; m<DIM; m++) {
          for (int j=0; j<ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; j++) {
            d_grad_n_dmesh[k][0][m][j] = dn_dcsi[k]*(-1.0f)/det_J/det_J*d_det_J_dmeshkj[m][j];
          }
        }
      }
    break;

    case NCM_MAPPING:
    default:
    for (k=0; k<pd->Num_Dim; k++) {
      normal[k] = fv->snormal[k];
      for (l=0; l<pd->Num_Dim; l++) {
        grad_normal[k][l] = 0.0;
        for (m=0; m<pd->Num_Dim; m++) {
          for (i=0; i<ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; i++) {
            d_grad_n_dmesh[k][l][m][i] = 0.0;
          }
        }
      }
    }
    break;
  }

  double d_grad_d_dmesh[DIM][DIM][DIM][MDE];
  memset(d_grad_d_dmesh, 0.0, sizeof(double)*DIM*DIM*DIM*MDE);

  switch(mp->ehl_integration_kind) {
    case SIK_S:
      // in domain_s the basis direction 0 (l=0) is parallel to the element
      // need to make ds/dcsi = det_J
      detJ_2d_bar(&det_J, d_det_J_dmeshkj);

      //ddisp_dcsi
      double csigrad[DIM];
      memset(csigrad, 0.0, sizeof(double)*DIM);

      for (int i=0; i<ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; i++){ // elemental variable degree of freedom i
          csigrad[0] += *esp->d[0][i]*bf[MESH_DISPLACEMENT1]->dphidxi[i][0];
          csigrad[1] += *esp->d[1][i]*bf[MESH_DISPLACEMENT2]->dphidxi[i][0];
      }
      for (int k=0; k<DIM; k++) { // displacement index k
        grad_disp[k][0] = csigrad[k]/det_J;
        for (int m=0; m<DIM; m++) {
          for(int j=0; j<ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; j++){ // sensitivity to displacement k degree of freedom j
                         //   gradient direction 'csi'
                         //   |  sensitivity to this displacement
                         //   |  |
                         //   v  v
            d_grad_d_dmesh[k][0][m][j] = csigrad[k]*(-1.0f)/det_J/det_J*d_det_J_dmeshkj[m][j]
                                       + delta(k,m)*bf[MESH_DISPLACEMENT1]->dphidxi[j][0]/det_J;
          }
        }
      }
      // switch indices and load into II vars because we don't need another set of intermediate variables
      for (int k = 0; k<DIM; k++) { // vector element k
          for (int l = 0; l<DIM; l++) { // derivative direction l
            gradII_n[l][k] = grad_normal[k][l];
            gradII_disp[l][k] = grad_disp[k][l];
            for (int i=0; i<ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; i++) { // element degree of freedom i
              for (int m=0; m<DIM; m++) { // sensitivity to mesh displacement m
                dgradII_n_dmesh[m][k][l][i]    = d_grad_n_dmesh[k][l][m][i];
                dgradII_disp_dmesh[m][k][l][i] = d_grad_d_dmesh[k][l][m][i];
              }
            }
          }
        }

    break;
  case SIK_XY: // not really implemented yet but this kinda works
      for (k = 0; k<DIM; k++) { // vector element k
        for (l = 0; l<DIM; l++) { // derivative direction l
          grad_n[l][k] = grad_normal[k][l];
          grad_disp[l][k] = fv->grad_d[k][l];
          for (i=0; i<ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; i++) { // element degree of freedom i
            for (m=0; m<DIM; m++) { // sensitivity to mesh displacement m
              dgrad_n_dmesh[m][k][l][i] = fv->d_grad_n_dmesh[k][l][m][i];
              dgrad_disp_dmesh[m][k][l][i] = fv->d_grad_d_dmesh[k][l][m][i];
            }
          }
        }
      }
      for (k = 0; k<DIM; k++) {
          ShellRotate(grad_n[k], dgrad_n_dmesh[k], gradII_n[k], dgradII_n_dmesh[k], n_dof[MESH_DISPLACEMENT1]);
          ShellRotate(grad_disp[k], dgrad_disp_dmesh[k], gradII_disp[k], dgradII_disp_dmesh[k], n_dof[MESH_DISPLACEMENT1]);
      }
    break;
  }

  // assume the web is always the lower member of the gap boundary
  double n_dot_up = 0.0;

  double up[DIM];
  if (pd->Num_Dim == 2) {
    up[0] = up[2] = 0.0;
    up[1] = 1.0;
  } else if (pd->Num_Dim == 3){
    up[0] = up[1] = 0.0;
    up[2] = 1.0;
  }
  for (int k=0; k<DIM; k++) {
    n_dot_up += normal[k]*up[k];
  }

  if (n_dot_up > 0.0) {
    h_sign = 1.0;
  } else if (n_dot_up < 0.0) {
    h_sign = -1.0;
  } else {
    EH(-1, "Something is wrong! n_dot_up == 0.0");
    return;
  }

  // zero the arrays
  //memset (dh_dmesh, 0.0, sizeof(double)*DIM*MDE);
  memset (dh_dnormal, 0.0, sizeof(double)*DIM*MDE);
  memset (d2h_dtime_dmesh, 0.0, sizeof(double)*DIM*MDE);
  memset (d2h_dtime_dnormal, 0.0, sizeof(double)*DIM*MDE);
  //memset (d_gradIIh_dmesh, 0.0, sizeof(double)*DIM*DIM*MDE);
  memset (d_gradIIh_dnormal, 0.0, sizeof(double)*DIM*DIM*MDE);

  // Setup Height function model and sensitivities to mesh motion, and normal
  switch ( mp->FSIModel ) {
    case FSI_SHELL_ONLY_MESH:
      for (k=0; k<DIM; k++) { // vector direction k
        *h -= h_sign*(normal[k])*(fv->d[k]);
        if(*h < 0.0) {
          WH(-1, "negative h in load_displacement_coupling_model()");
        }
        for (l = 0; l<DIM; l++) { // gradient direction l

          gradII_h[l] -= h_sign*grad_normal[k][l]*fv->d[k];
          gradII_h[l] -= h_sign*fv->n[k]*grad_disp[k][l];

        }
        if (pd->TimeIntegration == TRANSIENT ) {
          *(dh_dtime) -= h_sign*(fv->n[k] * fv_dot->d[k] + fv_dot->n[k] * fv->d[k]);

          for (l = 0; l<DIM; l++) { // gradient direction l
            for (i = 0; i<ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; i++) {
              d2h_dtime_dmesh[k][i] -= h_sign*fv->n[l]*delta(l,k)*bf[MESH_DISPLACEMENT1]->phi[i]*(1.0f+2.0f*tt)/delta_t;
              d2h_dtime_dmesh[k][i] -= h_sign*fv_dot->n[l]*delta(k,l)*bf[MESH_DISPLACEMENT1]->phi[i];


            }
          }

          if (pd->e[pg->imtrx][R_SHELL_NORMAL1]) {
            for (i = 0; i<ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; i++) {
              d2h_dtime_dnormal[k][i] -= h_sign*fv_dot->d[k]*bf[SHELL_NORMAL1]->phi[i];
              d2h_dtime_dnormal[k][i] -= h_sign*fv->d[k]*bf[SHELL_NORMAL1]->phi[i]*(1.0+2.0*tt)/delta_t;
            }
          }
        }
        for(i = 0; i<ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; i++) {
          dh_dmesh[k][i] -= h_sign*normal[k]*bf[MESH_DISPLACEMENT1]->phi[i];
          if (pd->v[SHELL_NORMAL1 + k]) {
            dh_dnormal[k][i] -= h_sign*fv->d[k]*bf[SHELL_NORMAL1]->phi[i];
          }
        }

        for (l = 0; l<DIM; l++){ // gradient direction l
          for (i = 0; i<ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; i++){ // element degree of freedom i
            for (int m = 0; m<DIM; m++) { // sensitivity to displacement m

              d_gradIIh_dmesh[l][m][i] -= h_sign*d_grad_n_dmesh[k][l][m][i]*fv->d[k]
                                        + h_sign*grad_normal[k][l]*delta(m,k)*bf[MESH_DISPLACEMENT1]->phi[i];
              d_gradIIh_dmesh[l][m][i] -= h_sign*fv->n[k]*d_grad_d_dmesh[k][l][m][i];


              eqn = R_MESH1;

              ShellBF(eqn, i, &phi_i, grad_phi_i, gradII_phi_i, d_gradII_phi_i_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map);


              d_gradIIh_dnormal[l][m][i] -= h_sign*gradII_phi_i[l]*delta(m,k)*fv->d[k];
              d_gradIIh_dnormal[l][m][i] -= h_sign*phi_i*grad_disp[k][l]*delta(k,m);
            }
          }
        }
      }
      break;
    default:
      break;
  }
}

void rmesh_minus_rroller (
  double tt,
  double delta_t,
  double *h,
  double *dh_dtime,
  double *gradII_h,
  double dh_dmesh[][MDE],
  double dh_dnormal[][MDE],
  double d2h_dtime_dmesh[][MDE],
  double d2h_dtime_dnormal[][MDE],
  double d_gradIIh_dmesh[][DIM][MDE],
  double d_gradIIh_dnormal[][DIM][MDE],
  int *n_dof,
  int *dof_map
) {
  // Check to see if Height_UFunctionModel is ROLLER, otherwise break
  if( mp->HeightUFunctionModel != ROLLER){
    EH(-1,"I cannot roll without an 'Upper Height Function Model = ROLLER'  !!");
  }

  //notice there is no more need of normal
  memset(dh_dnormal, 0.0, sizeof(double)*DIM*MDE);
  memset(d2h_dtime_dnormal, 0.0, sizeof(double)*DIM*MDE);
  memset(d_gradIIh_dnormal, 0.0, sizeof(double)*DIM*DIM*MDE);

  // actually no other source of h matters, yeah?
  memset (dh_dmesh, 0.0, sizeof(double)*DIM*MDE);
  memset (d2h_dtime_dmesh, 0.0, sizeof(double)*DIM*MDE);
  memset (d_gradIIh_dmesh, 0.0, sizeof(double)*DIM*DIM*MDE);

  // need the position of the roller center, roller radius and position of the Gauss point.
  double hmin = mp->u_heightU_function_constants[0];
  double rroller = mp->u_heightU_function_constants[1];
  double xrc = mp->u_heightU_function_constants[2];
  //double exfield_mult = mp->u_heightU_function_constants[3]; // should not be used here.

  double yrc = hmin + rroller;
  double xc[DIM];
  xc[0] = xrc;
  xc[1] = yrc;
  xc[2] = 0.0;

  //position
  double x[DIM];
  x[0] = fv->x[0];
  x[1] = fv->x[1];

  //the radial distance from roller center to web
  double rmesh = sqrt((x[0] - xrc)*(x[0] - xrc) + (x[1] - yrc)*(x[1] - yrc));

  //need the det_J_2d_bar
  double det_J;
  double d_det_J_dmeshkj[DIM][MDE];
  detJ_2d_bar(&det_J, d_det_J_dmeshkj);
  //get dh/dcsi

  // but first need dxw_dcsi

  double dwx_dcsi[DIM], d_dwx_dcsi_dmesh[DIM][DIM][MDE];

  dpos_dcsi(dwx_dcsi, d_dwx_dcsi_dmesh);

  double dh_dcsi = ( (x[0]-xc[0])*dwx_dcsi[0] + (x[1] - xc[1])*dwx_dcsi[1] )/rmesh;

  // set h
  *h = rmesh - rroller;

  // set dh_dtime
  if( pd->TimeIntegration == TRANSIENT) {
    *dh_dtime = ( (x[0]-xrc)*fv_dot->x[0] + (x[1] - yrc)*fv_dot->x[1] )/rmesh;
  }

  // set gradII_h
  if( mp->ehl_integration_kind == SIK_S) {
    gradII_h[0] = dh_dcsi/det_J;
    gradII_h[1] = 0.0;
    gradII_h[2] = 0.0;

  } else {
    // don't know how to do anything else yet
    EH(-1, "I don't know how to set gradII_h unless 'Elastohydrodynamic Lubrication Shell Integration Kind' is 'S'.\n"
           "Check your material file or implement new features.");
  }
  // set dh_dmesh
  for (int i=0; i<ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; i++){
    dh_dmesh[0][i] = (x[0] - xrc)/rmesh*bf[MESH_DISPLACEMENT1]->phi[i];
  }
  for (int i=0; i<ei[pg->imtrx]->dof[MESH_DISPLACEMENT2]; i++){
    dh_dmesh[1][i] = (x[1] - yrc)/rmesh*bf[MESH_DISPLACEMENT2]->phi[i];
  }

  // set d2h_dtime_dmesh
  if( pd->TimeIntegration == TRANSIENT) {
    for (int k=0; k<DIM-1; k++) {
      for (int j=0; j<ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; j++){
        d2h_dtime_dmesh[k][j] = bf[MESH_DISPLACEMENT1+k]->phi[j]/rmesh
                                *((x[k]-xc[k])*(1.0 + 2.0*tt)/delta_t + fv_dot->x[k])
            -1.0/rmesh/rmesh*(*dh_dtime)*(x[k] - xc[k])
            *bf[MESH_DISPLACEMENT1+k]->phi[j];
      }
    }
  }
  // set d_gradIIh_dmesh
  for (int m=0; m<DIM-1; m++) {
    for (int j=0; j<ei[pg->imtrx]->dof[MESH_DISPLACEMENT1+m]; j++) {
      double d_dh_dcsi_dmesh[DIM][MDE];
      d_dh_dcsi_dmesh[m][j] =
          (
            (x[m] - xc[m])*d_dwx_dcsi_dmesh[m][m][j]
            + bf[MESH_DISPLACEMENT1]->phi[j]*dwx_dcsi[m]
          )/rmesh
          +
          (
            ((x[0]-xc[0])*dwx_dcsi[0] + (x[1] - xc[1])*dwx_dcsi[1])
            *-1.0/rmesh/rmesh/rmesh
            *bf[MESH_DISPLACEMENT1]->phi[j]
            *(x[m]-xc[m])
          );

      d_gradIIh_dmesh[0][m][j] = d_dh_dcsi_dmesh[m][j]/det_J
                               + -dh_dcsi/det_J/det_J
                                 *d_det_J_dmeshkj[m][j];
    }
  }
}

void dpos_dcsi(double dxw_dcsi[], double d_dxw_dcsi_dmesh[][DIM][MDE]) {
  double d_phi_dxi[MDE], d_sh_x_dxi, d_sh_y_dxi;
  int eqn, var, index, node;
  double displacement[DIM][MDE];
  memset(displacement, 0.0f, sizeof(double)*DIM*MDE);

  eqn = R_MESH1;
  var = MESH_DISPLACEMENT1;

  switch(mp->FSIModel) {
    case FSI_SHELL_ONLY_MESH:

      eqn = R_MESH1;
      var = MESH_DISPLACEMENT1;
      for (int i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        displacement[0][i] += *esp->d[0][i];
        displacement[1][i] += *esp->d[1][i];
      }

    break;

    case FSI_SHELL_ONLY:
      // there WAS a reason this is here, but it looks like
      // it should be broken now
      eqn = R_TFMP_BOUND;
      var = TFMP_PRES;
    break;
  }
  for (int i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
    d_phi_dxi[i] = bf[eqn]->dphidxi[i][0];
  }

  d_sh_x_dxi = d_sh_y_dxi = 0.;
  for (int i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
    node = ei[pg->imtrx]->dof_list[eqn][i];
    index = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr +node];

    d_sh_x_dxi +=    (Coor[0][index] + displacement[0][i]) * d_phi_dxi[i];
    d_sh_y_dxi +=    (Coor[1][index] + displacement[1][i]) * d_phi_dxi[i];
  }
  dxw_dcsi[0] = d_sh_x_dxi;
  dxw_dcsi[1] = d_sh_y_dxi;
  dxw_dcsi[2] = 0.0;

  for(int k=0; k<DIM-1; k++) {
    for (int m=0; m<DIM-1; m++){
      for (int j=0; j<ei[pg->imtrx]->dof[var]; j++) {
        d_dxw_dcsi_dmesh[k][m][j] = delta(m,k)*bf[var]->dphidxi[j][0];
      }
    }
  }

}

double dxdcsi(double *x_node, int x_var) {
  double dxdcsi = 0.0;
  int dof = ei[pg->imtrx]->dof[x_var];
  for(int i=0;i<dof;i++){
    dxdcsi += x_node[i]*bf[x_var]->dphidxi[i][0];
  }
  return dxdcsi;
}

// get s-based mapping determinates and mesh derivatives for 2d bar elements

void detJ_2d_bar(
  double *det_J,
  double d_det_J_dmesh[][MDE]
) {
  int eqn, var, i, node, index;
  double d_phi_dxi[MDE], d_sh_x_dxi, d_sh_y_dxi;

  double displacement[DIM][MDE];
  memset(displacement, 0.0f, sizeof(double)*DIM*MDE);

  eqn = R_MESH1;
  var = MESH_DISPLACEMENT1;

  switch(mp->FSIModel) {
    case FSI_SHELL_ONLY_MESH:

      eqn = R_MESH1;
      var = MESH_DISPLACEMENT1;
      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        displacement[0][i] += *esp->d[0][i];
        displacement[1][i] += *esp->d[1][i];
      }

    break;
    case FSI_SHELL_ONLY:
      eqn = R_TFMP_BOUND;
      var = TFMP_PRES;
    break;
  }
  for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
    d_phi_dxi[i] = bf[eqn]->dphidxi[i][0];
  }

  d_sh_x_dxi = d_sh_y_dxi = 0.;
  for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
    node = ei[pg->imtrx]->dof_list[eqn][i];
    index = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr +node];

    d_sh_x_dxi +=    (Coor[0][index] + displacement[0][i]) * d_phi_dxi[i];
    d_sh_y_dxi +=    (Coor[1][index] + displacement[1][i]) * d_phi_dxi[i];
  }

  *det_J = sqrt(d_sh_x_dxi*d_sh_x_dxi + d_sh_y_dxi*d_sh_y_dxi);
  if (mp->FSIModel == FSI_SHELL_ONLY_MESH) {
    for (int j = 0; j<ei[pg->imtrx]->dof[var]; j++) {
      d_det_J_dmesh[0][j] = d_sh_x_dxi*d_phi_dxi[j]/(*det_J);
      d_det_J_dmesh[1][j] = d_sh_y_dxi*d_phi_dxi[j]/(*det_J);
    }
  }
}

void ShellBF_2d_bar(
  int eq_var,
  int ii,
  double gradII_phi_i[DIM],
  double d_gradII_phi_i_dx[DIM][DIM][MDE]
) {
  memset(gradII_phi_i, 0.0, sizeof(double)*DIM);
  memset(d_gradII_phi_i_dx, 0.0, sizeof(double)*DIM*DIM*MDE);
  double det_J;
  double d_det_J_dmesh[DIM][MDE];

  detJ_2d_bar(&det_J, d_det_J_dmesh);

  gradII_phi_i[0] = bf[eq_var]->dphidxi[ii][0]/det_J;
  for (int j=0; j<ei[pg->imtrx]->dof[eq_var]; j++) {
    d_gradII_phi_i_dx[0][0][j] = bf[eq_var]->dphidxi[ii][0]*(-1.0)/det_J/det_J*d_det_J_dmesh[0][j];
    d_gradII_phi_i_dx[0][1][j] = bf[eq_var]->dphidxi[ii][0]*(-1.0)/det_J/det_J*d_det_J_dmesh[1][j];
  }
}

void load_gap_model(GAP_STRUCT *gap) {
  int var;
  int fp_type = FP_NORMAL;
  double det_J, d_det_J_dmeshkj[DIM][MDE];
  // set zeros
  memset(gap->dh_dmesh, 0.0, sizeof(double)*DIM*MDE);
  memset(gap->dh_dnormal, 0.0, sizeof(double)*DIM*MDE);

  memset(gap->d2h_dtime_dmesh, 0.0, sizeof(double)*DIM*MDE);
  memset(gap->d2h_dtime_dnormal, 0.0, sizeof(double)*DIM*MDE);

  memset(gap->gradII_h, 0.0, sizeof(double)*DIM);
  memset(gap->d_gradIIh_dmesh, 0.0, sizeof(double)*DIM*DIM*MDE);
  memset(gap->d_gradIIh_dnormal, 0.0, sizeof(double)*DIM*DIM*MDE);

  double gradII_x[DIM];
  double dgradII_x_dmesh[DIM][DIM][MDE];
  double csigrad[DIM];
  // load up the dx_ds and mesh sensitivities
  if (mp->ehl_integration_kind == SIK_S) {
    detJ_2d_bar(&det_J, d_det_J_dmeshkj);
    double* grad;

    if (pd->Num_Dim == 2 && ei[pg->imtrx]->ielem_type == LINEAR_BAR) {
      var = MESH_DISPLACEMENT1;
      grad = gradII_x;
      memset(grad, 0.0, sizeof(double)*DIM);
      memset (csigrad, 0.0, sizeof(double)*DIM);
      memset (dgradII_x_dmesh, 0.0, sizeof(double)*DIM*DIM*MDE);
      for (int i=0; i<ei[pg->imtrx]->dof[var]; i++) {
        int node = ei[pg->imtrx]->dof_list[R_MESH1][i];
        int index = Proc_Elem_Connect[Proc_Connect_Ptr[ei[pg->imtrx]->ielem] +node];
        csigrad[0] += (Coor[0][index] + *esp->d[0][i])*bf[var]->dphidxi[i][0];
      }

      grad[0] = csigrad[0]/det_J;

      for (int i=0; i<ei[pg->imtrx]->dof[var]; i++) {
        for (int k=0; k<DIM; k++) {
          dgradII_x_dmesh[0][k][i] = csigrad[0]*(-1.0)/det_J/det_J*d_det_J_dmeshkj[k][i]
                                    + delta(0,k)*bf[var]->dphidxi[i][0]/det_J;
        }
      }
    }
  }

  /* Use the height_function_model */
  //double h, H_U, dH_U_dtime, H_L, dH_L_dtime;
  double H_U, dH_U_dtime, H_L, dH_L_dtime;
  double dH_U_dX[DIM],dH_L_dX[DIM], dH_U_dp, dH_U_ddh;

  gap->h = height_function_model(&H_U, &dH_U_dtime, &H_L, &dH_L_dtime,
         dH_U_dX, dH_L_dX, &dH_U_dp, &dH_U_ddh, gap->time, gap->delta_t);

  if (fpclassify(gap->h)!= fp_type && gap->h != 0.0) {
    EH(-1, "mass term is not normal after retrieving initial value from height_function_model()");
  }

  for (int l=0; l<DIM; l++)  {
    for (int j=0; j<ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; j++){
      gap->dh_dmesh[l][j] += dH_U_dX[l]*bf[MESH_DISPLACEMENT1]->phi[j];
    }
  }

  gap->dh_dtime = dH_U_dtime - dH_L_dtime;



  if (mp->ehl_integration_kind == SIK_S) {
    if (mp->HeightUFunctionModel == ROLLER) {
      for (int k=0; k<DIM; k++) {
        gap->gradII_h[k] = dH_U_dX[k]*gradII_x[k] - dH_L_dX[k];
      }
      double r = mp->u_heightU_function_constants[1];
      double xc = mp->u_heightU_function_constants[2];
      double x = fv->x[0];

      for (int j=0; j<ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; j++){
        gap->d_gradIIh_dmesh[0][0][j] += 1.0/sqrt(CUBE(SQUARE(r) - SQUARE(x - xc)))
                                    *(x - xc)
                                    *bf[MESH_DISPLACEMENT1]->phi[j]
                                    *(x - xc)
                                    *gradII_x[0];
        gap->d_gradIIh_dmesh[0][0][j] += 1.0/sqrt(SQUARE(r) - SQUARE(x - xc))
                                    *bf[MESH_DISPLACEMENT1]->phi[j]
                                    *gradII_x[0];
        gap->d_gradIIh_dmesh[0][0][j] += 1.0/sqrt(SQUARE(r) - SQUARE(x - xc))
                                    *(x - xc)
                                    *dgradII_x_dmesh[0][0][j];
        gap->d_gradIIh_dmesh[0][1][j] += 1.0/sqrt(SQUARE(r) - SQUARE(x - xc))
                                    *(x - xc)
                                    *dgradII_x_dmesh[0][1][j];
      }
    } else {
      for (int k = 0; k<DIM; k++) {
        gap->gradII_h[k] = dH_U_dX[k] - dH_L_dX[k];
      }
    }
  }

  if (mp->FSIModel == FSI_SHELL_ONLY_MESH) {
    load_displacement_coupling_model(
          gap->tt,
          gap->delta_t,
          &(gap->h),
          &(gap->dh_dtime),
          gap->gradII_h,
          gap->dh_dmesh,
          gap->dh_dnormal,
          gap->d2h_dtime_dmesh,
          gap->d2h_dtime_dnormal,
          gap->d_gradIIh_dmesh,
          gap->d_gradIIh_dnormal,
          gap->n_dof,
          gap->dof_map
    );
  }

}

void load_roller_normal_into_fv(void) {
  // I think this function relies on shell_determinant_and_normal
  // to init fv->snormal

  double r = mp->u_heightU_function_constants[1];
  double x  = fv->x[0];
  double xc = mp->u_heightU_function_constants[2];
  double unnormalized_n[DIM], mag_unnormalized_n;
  unnormalized_n[0] = x - xc;
  unnormalized_n[1] = sqrt(r*r - unnormalized_n[0]*unnormalized_n[0]);
  mag_unnormalized_n = sqrt(unnormalized_n[0]*unnormalized_n[0]
                          + unnormalized_n[1]*unnormalized_n[1]);
  fv->snormal[0] = unnormalized_n[0]/mag_unnormalized_n;
  fv->snormal[1] = unnormalized_n[1]/mag_unnormalized_n;

  for (int j=0; j<ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; j++) {
    fv->dsnormal_dx[0][0][j] = 1./mag_unnormalized_n
                               *bf[MESH_DISPLACEMENT1]->phi[j];
    fv->dsnormal_dx[1][0][j] = -unnormalized_n[0]
                               /unnormalized_n[1]/mag_unnormalized_n
                               *bf[MESH_DISPLACEMENT1]->phi[j];
    fv->dsnormal_dx[0][1][j] = 0.0;
    fv->dsnormal_dx[1][1][j] = 0.0;
  }


  return;
}
