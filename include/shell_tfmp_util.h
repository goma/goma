/* Copyright 2019 Andrew Cochrane
 * This file released under MIT license
 * */

#include "mm_shell_bc.h"
#include "shell_tfmp_struct.h"
#include "el_elm.h"

#ifndef _SHELL_TFMP_UTIL
#define _SHELL_TFMP_UTIL


//Still Experimental
enum clipping_kind{none, continuity, restorative, constant_sat, var_swap};

// Make global variables (really? - eesh!)
EXTERN enum clipping_kind my_clipping_kind;

// load_viscosity_model(&mu_l, &mu_g)
// OUTPUTS:
//  mu_l = constant liquid viscosity
//  mu_g = constant gas viscosity
EXTERN void load_tfmp_viscosity_model
(double*,
       double*
);

// load_gas_density_model(&Patm, &rho_g, &drho_g_dP)
// OUTPUTS:
//  Patm      = atmospheric pressure (important for
//              ideal gas model)
//  rho_g     = computed gas density
//  drho_g_dP = partial wrt pressure
EXTERN void load_gas_density_model
(double*,
       double*,
       double*
);

// load_molecular_diffusion_model(S, &D, &Krd, &dKrd_dS)
// OUTPUTS:
//  D       = Diffusivity 
//  Krd     = molecular diffusion correction factor
//  dKrd_dS = partial wrt to saturation
// INPUTS:
//  S       = saturation
EXTERN void load_molecular_diffusion_model
(double,
       double*,
       double*,
       double*
);

// load_relative_permeability_model(S, &Krl, &dKrl_dS, 
//                                  &Krg, &dKrg_dS)
// OUTPUTS:
//  Krl     = liquid relative permeability
//  dKrl_dS = partial wrt saturation
//  Krg     = gas realtive permeability
//  dKrg_dS = partial wrt saturation
// INPUTS:
//  S       = saturation
EXTERN void load_relative_permeability_model
(double,
       double*,
       double*,
       double*,
       double*
);

// load_gas_dissolution_model(h, Patm, &J, &dJ_dP, &dJ_dS,
//                            &dJ_dh)
// OUTPUTS:
//  J     = volume average diffusive flux
//  dJ_dP = partial wrt pressure
//  dJ_dS = partial wrt saturation
//  dJ_dh = partial wrt to gap thickness
// INPUTS:
//  h     = gap thickness
//  Patm  = atmospheric pressure
EXTERN void load_gas_dissolution_model
(double,
       double,
       double*,
       double*,
       double*,
       double*
);

// load_displacement_coupling_model(tt, delta_t, &h,
//                                  &dh_dtime, gradII_h,
//                                  dh_dmesh, dh_dnormal,
//                                  d2h_dtime_dmesh,
//                                  d2h_dtime_dnormal,
//                                  d_gradIIh_dmesh,
//                                  d_gradIIh_dnormal,
//                                  n_dof, dof_map)
// OUTPUTS:
//  h                 = gap thickness
//  dh_dtime          = partial wrt time
//  gradII_h          = gap thickness gradient
//  dh_dmesh          = partial wrt to mesh motion
//  dh_dnormal        = partial wrt the normal
//  d2h_dtime_dmesh   = 2nd partial wrt to time and 
//                      mesh motion
//  d2h_dtime_dnormal = 2nd partial wrt to time and
//                      the normal
//  d_gradIIh_dmesh   = partial wrt mesh motion
//  d_gradIIh_dnormal = partial wrt to normal
// INPUTS:
//  tt                = time step parameter (BE or CN)
//  delta_t           = period of time step
//  n_dof             = number of degrees of freedom of
//                      mesh eq
//  dof_map           = used in ShellBF for mapping
//                      between basis function gradient
//                      sensitivity wrt to mesh
//                      dispalcement and shell basis
//                      function gradient sensitivity
//                      wrt to mesh displacement
EXTERN void load_displacement_coupling_model
(double,
       double,
       double*,
       double*,
       double*,             // gradII_h
       double[][MDE],       // dh_dmesh
       double[][MDE],       // dh_dnormal
       double[][MDE],       // d2h_dtime_dmesh
       double[][MDE],       // d2h_dtime_dnormal
       double[][DIM][MDE],  // d_gradIIh_dmesh
       double[][DIM][MDE],  // d_gradIIh_dnormal
       int*,
       int*
);

EXTERN void h0_minus_ndotd
(double,
       double,
       double*,
       double*,
       double*,             // gradII_h
       double[][MDE],       // dh_dmesh
       double[][MDE],       // dh_dnormal
       double[][MDE],       // d2h_dtime_dmesh
       double[][MDE],       // d2h_dtime_dnormal
       double[][DIM][MDE],  // d_gradIIh_dmesh
       double[][DIM][MDE],  // d_gradIIh_dnormal
       int*,
       int*
);

EXTERN void rmesh_minus_rroller
(double,
       double,
       double*,
       double*,
       double*,             // gradII_h
       double[][MDE],       // dh_dmesh
       double[][MDE],       // dh_dnormal
       double[][MDE],       // d2h_dtime_dmesh
       double[][MDE],       // d2h_dtime_dnormal
       double[][DIM][MDE],  // d_gradIIh_dmesh
       double[][DIM][MDE],  // d_gradIIh_dnormal
       int*,
       int*
);

EXTERN void dpos_dcsi
(double[], double[][DIM][MDE]);

EXTERN double dxdcsi
(double*, int);

EXTERN void detJ_2d_bar
(double*, double[][MDE]);

EXTERN void load_gap_model
(GAP_STRUCT*);

EXTERN void load_roller_normal_into_fv
(void);

#endif // _SHELL_TFMP_UTIL
