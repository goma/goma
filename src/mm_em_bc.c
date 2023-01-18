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

#include <complex.h>
#undef I

#include "mm_as.h"
#include "mm_eh.h"
#include "mm_em_bc.h"
#include "mm_fill_em.h"
#include "mm_mp.h"
#include "mm_qtensor_model.h"
#include "std.h"

int apply_em_farfield_direct_vec(double func[DIM],
                                 double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                                 double xi[DIM], /* Local stu coordinates */
                                 const int bc_name,
                                 double *bc_data) {
  /***********************************************************************
   * TODO AMC: rewrite this description
   * apply_em_farfield_direct():
   *
   *  Function which evaluates the expression specifying the
   *  a plane-wave directly incident (parallel to normal) on
   *  the boundary.
   *
   *  n_bound CROSS E = -n_bound CROSS E
   *                  * ( eta_2 * kappa_2)
   *                  / ( omega * mu_2 )
   *                  * ( 1 - 2*I)
   *
   *   func = -
   *
   *
   *  The boundary condition EM_DIRECT_BC employs this
   *  function.
   *
   *
   * Input:
   *
   *  em_eqn   = which equation is this applied to
   *  em_var   = which variable is this value sensitive to
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
   *   Author: Andrew Cochrane (10/9/2019)
   *
   ********************************************************************/

  int var;
  dbl mag_permeability = 12.57e-07; // H/m
  double n1, n2;                    /* Refractive index */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_n1_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_n1 = &d_n1_struct;

  double k1, k2; /* Extinction Coefficient */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k1_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k1 = &d_k1_struct;

  double normal[DIM]; // surface normal vector

  // double complex E1[DIM], H1[DIM]; // complex fields inside domain

  // Need material properties for both sides of interface

  // use mp for inside (subscript 1) ..

  double omega = upd->EM_Frequency;
  n1 = refractive_index(d_n1, 0.0);

  k1 = extinction_index(d_k1, 0.0);
  // Compute complex impedance
  complex cpx_refractive_index1, cpx_rel_permittivity1, cpx_permittivity1, impedance1, kappa1;

  cpx_refractive_index1 = n1 + _Complex_I * k1;
  cpx_rel_permittivity1 = SQUARE(cpx_refractive_index1);
  cpx_permittivity1 = cpx_rel_permittivity1 * mp->permittivity;

  impedance1 = csqrt(mag_permeability / cpx_permittivity1);
  kappa1 = omega * impedance1 * cpx_permittivity1;

  // use BC input for outside (subscript 2)
  n2 = bc_data[0];
  k2 = bc_data[1];
  // n2 = 1.000293; // air (wikipedia 2019)
  // k2 = 0;
  //  Compute complex impedance
  complex cpx_refractive_index2, cpx_rel_permittivity2, cpx_permittivity2, impedance2;

  cpx_refractive_index2 = n2 + _Complex_I * k2;
  cpx_rel_permittivity2 = SQUARE(cpx_refractive_index2);
  cpx_permittivity2 = cpx_rel_permittivity2 * mp->permittivity;

  impedance2 = csqrt(mag_permeability / cpx_permittivity2);

  // need Surface Normal vector
  for (int p = 0; p < DIM; p++) {
    normal[p] = fv->snormal[p];
  }

  complex Gamma, tau, incidentE[DIM];
  complex incidentH[DIM] = {0.0};
  Gamma = (impedance2 - impedance1) / (impedance2 + impedance1);
  tau = (2.0 * impedance2) / (impedance2 + impedance1);

  double complex reduction_factor;

  switch (bc_name) {
  case EM_ER_FARFIELD_DIRECT_BC:
  case EM_EI_FARFIELD_DIRECT_BC:
    reduction_factor = -_Complex_I * tau / kappa1 / (1 + Gamma);
    break;
  case EM_HR_FARFIELD_DIRECT_BC:
  case EM_HI_FARFIELD_DIRECT_BC:
    reduction_factor = -_Complex_I * tau / kappa1 / (1 - Gamma);
    break;
  default:
    reduction_factor = 0.0;
    break;
  }

  incidentE[0] = bc_data[2] + _Complex_I * bc_data[5];
  incidentE[1] = bc_data[3] + _Complex_I * bc_data[6];
  incidentE[2] = bc_data[4] + _Complex_I * bc_data[7];

  for (int p = 0; p < DIM; p++) {
    for (int q = 0; q < DIM; q++) {
      for (int r = 0; r < DIM; r++) {
        incidentH[p] += permute(p, q, r) * normal[q] * incidentE[r] / impedance2;
      }
    }
  }

  // construct curls and sensitivities
  // assuming normal is purely real.
  double Re_curl_E[DIM] = {0.0};
  double Im_curl_E[DIM] = {0.0};
  double Re_curl_H[DIM] = {0.0};
  double Im_curl_H[DIM] = {0.0};
  double d_dERb_Re_curl_E[DIM][DIM][MDE] = {{{0.0}}};
  double d_dEIb_Im_curl_E[DIM][DIM][MDE] = {{{0.0}}};
  double d_dHRb_Re_curl_H[DIM][DIM][MDE] = {{{0.0}}};
  double d_dHIb_Im_curl_H[DIM][DIM][MDE] = {{{0.0}}};

  for (int p = 0; p < pd->Num_Dim; p++) {
    for (int q = 0; q < pd->Num_Dim; q++) {
      for (int r = 0; r < pd->Num_Dim; r++) {
        Re_curl_E[p] += permute(p, q, r) * fv->grad_em_er[r][q];
        Im_curl_E[p] += permute(p, q, r) * fv->grad_em_ei[r][q];
        Re_curl_H[p] += permute(p, q, r) * fv->grad_em_hr[r][q];
        Im_curl_H[p] += permute(p, q, r) * fv->grad_em_hi[r][q];
        // assuming all variables have same degrees of freedom

        for (int b = 0; b < ei[pg->imtrx]->dof[EM_E1_REAL]; b++) {
          d_dERb_Re_curl_E[p][r][b] += permute(p, q, r) * bf[EM_E1_REAL + r]->grad_phi[b][q];
        }
        for (int b = 0; b < ei[pg->imtrx]->dof[EM_E1_IMAG]; b++) {
          d_dEIb_Im_curl_E[p][r][b] += permute(p, q, r) * bf[EM_E1_IMAG + r]->grad_phi[b][q];
        }
        for (int b = 0; b < ei[pg->imtrx]->dof[EM_H1_REAL]; b++) {
          d_dHRb_Re_curl_H[p][r][b] += permute(p, q, r) * bf[EM_H1_REAL + r]->grad_phi[b][q];
        }
        for (int b = 0; b < ei[pg->imtrx]->dof[EM_H1_IMAG]; b++) {
          d_dHIb_Im_curl_H[p][r][b] += permute(p, q, r) * bf[EM_H1_IMAG + r]->grad_phi[b][q];
        }
      }
    }
  }

  complex cpx_func[DIM] = {0.0};

  // double real, imag;
  switch (bc_name) {
  case EM_ER_FARFIELD_DIRECT_BC:
  case EM_EI_FARFIELD_DIRECT_BC:
    for (int p = 0; p < DIM; p++) {
      cpx_func[p] += (Re_curl_E[p] + _Complex_I * Im_curl_E[p]) * reduction_factor;
    }
    for (int p = 0; p < pd->Num_Dim; p++) {
      for (int q = 0; q < DIM; q++) {
        for (int r = 0; r < DIM; r++) {
          cpx_func[p] += permute(p, q, r) * normal[q] * incidentE[r];
        }
      }
    }
    break;
  case EM_HR_FARFIELD_DIRECT_BC:
  case EM_HI_FARFIELD_DIRECT_BC:
    for (int p = 0; p < DIM; p++) {
      cpx_func[p] += (Re_curl_H[p] + _Complex_I * Im_curl_H[p]) * reduction_factor;
    }
    for (int p = 0; p < pd->Num_Dim; p++) {
      for (int q = 0; q < DIM; q++) {
        for (int r = 0; r < DIM; r++) {
          cpx_func[p] -= permute(p, q, r) * normal[q] * incidentH[r];
        }
      }
    }
    break;
  }

  switch (bc_name) {
  case EM_ER_FARFIELD_DIRECT_BC:

    for (int p = 0; p < DIM; p++) {
      func[p] = creal(cpx_func[p]);
      // func[0] = fv->grad_em_er[1][2] - fv->grad_em_er[2][1];
      /*
      for(int q=0; q<DIM; q++) {
        for (int r=0; r<DIM; r++){
          func[p] += permute(p,q,r)*fv->grad_em_er[r][q];
        }
      }*/
    }
    // eqn = R_EM_H*_REAL;
    var = EM_E1_REAL;
    //  real = 1.0;
    //  imag = 0.0;
    break;
  case EM_EI_FARFIELD_DIRECT_BC:
    for (int p = 0; p < DIM; p++) {
      func[p] = cimag(cpx_func[p]);
    }
    // eqn = R_EM_H*_IMAG;
    var = EM_E1_REAL;
    // real = 0.0;
    // imag = 1.0;
    break;
  case EM_HR_FARFIELD_DIRECT_BC:

    for (int p = 0; p < DIM; p++) {
      func[p] = creal(cpx_func[p]);
    }

    /*
          for (int p=0; p<DIM; p++) {
            func[p] = Im_n_x_H[p];
          }
          */
    // eqn = R_EM_E*_REAL;
    var = EM_H1_REAL;
    // real = 1.0;
    // imag = 0.0;
    break;
  case EM_HI_FARFIELD_DIRECT_BC:
    for (int p = 0; p < DIM; p++) {
      func[p] = cimag(cpx_func[p]);
    }
    // eqn = R_EM_E*_IMAG;
    var = EM_H1_REAL;
    // real = 0.0;
    // imag = 1.0;
    break;
  default:
    var = 0;
    // real = 0;
    // imag = 0;
    reduction_factor = 0;
    GOMA_EH(GOMA_ERROR, "Must call apply_em_farfield_direct with an applicable BC_NAME");
    return -1;
    break;
  }

  if (af->Assemble_Jacobian) {
    /*
        for (int j=0; j< ei[pg->imtrx]->dof[EM_E2_REAL]; j++){
              d_func[0][EM_E2_REAL][j] = -bf[EM_E2_REAL]->grad_phi[j][2];
              d_func[0][EM_E3_REAL][j] =  bf[EM_E3_REAL]->grad_phi[j][1];
                     0        1           -        1                  2
                     0        2           +        2                  1
                     p        r           +        r                  q
        }
      */
    /*
    for (int p=0; p<DIM; p++) {

      //for (int q=0; q<DIM; q++) {
        for (int r=0; r<DIM; r++) {
          for (int j=0; j<ei[pg->imtrx]->dof[EM_E1_REAL + r]; j++) {
            d_func[p][EM_E1_REAL + r][j] += -d_dERb_Re_curl_E[p][r][j];
          }
        }
      //}
    }
//    d_dERb_Re_curl_E[p][r][b]
    */
    for (int p = 0; p < pd->Num_Dim; p++) {
      // for (int q=0; q<pd->Num_Dim; q++) {
      for (int g = 0; g < pd->Num_Dim; g++) {
        int gvar = var + g;
        for (int j = 0; j < ei[pg->imtrx]->dof[gvar]; j++) {
          switch (bc_name) {
          case EM_ER_FARFIELD_DIRECT_BC:
            // d_func[p][gvar][j] += d_dERb_Re_curl_E[p][g][j]*creal(reduction_factor);
            // d_func[p][gvar+3][j] -= d_dEIb_Im_curl_E[p][g][j]*cimag(reduction_factor);
            d_func[p][gvar][j] += -d_dERb_Re_curl_E[p][g][j] * creal(reduction_factor);
            d_func[p][gvar + 3][j] += d_dERb_Re_curl_E[p][g][j] * cimag(reduction_factor);
            break;

          case EM_EI_FARFIELD_DIRECT_BC:
            d_func[p][gvar][j] += -d_dERb_Re_curl_E[p][g][j] * cimag(reduction_factor);
            d_func[p][gvar + 3][j] += -d_dEIb_Im_curl_E[p][g][j] * creal(reduction_factor);
            break;

          case EM_HR_FARFIELD_DIRECT_BC:
            d_func[p][gvar][j] += -d_dHRb_Re_curl_H[p][g][j] * creal(reduction_factor);
            d_func[p][gvar + 3][j] += d_dHIb_Im_curl_H[p][g][j] * cimag(reduction_factor);
            // d_func[p][gvar][j] += d_dHRb_Re_n_x_H[p][g][j];
            // d_func[p][gvar+3][j] += d_dERb_Re_n_x_E[p][g][j];
            break;

          case EM_HI_FARFIELD_DIRECT_BC:
            d_func[p][gvar][j] += -d_dHRb_Re_curl_H[p][g][j] * cimag(reduction_factor);
            d_func[p][gvar + 3][j] += -d_dHIb_Im_curl_H[p][g][j] * creal(reduction_factor);
            break;
          }
        }
        //}
      }
    }
  }
  return 0;
} // end of apply_em_direct_vec

/* AMC TODO:  Here's the working cross product w/sensitivity
  for (int p=0; p<pd->Num_Dim; p++) {
    for (int q=0; q<pd->Num_Dim; q++) {
      for (int r=0; r<pd->Num_Dim; r++) {
        //func[p] += fv->em_er[q] + fv->em_ei[q] + fv->em_hr[q] + fv->em_hi[q];
        func[p] += permute(p,q,r)* creal(normal[q]*E1[r]);
      }
    }
  }

  if(af->Assemble_Jacobian) {

    for (int p=0; p<pd->Num_Dim; p++) {
      for (int q=0; q<pd->Num_Dim; q++) {
        for (int g=0; g<3; g++) {
          int gvar = EM_E1_REAL + g;
          for (int j=0; j<ei[pg->imtrx]->dof[gvar]; j++) {


            //d_func[p][gvar][j] = phi_j;
            d_func[p][gvar][j] += permute(q,p,g)*creal(normal[q])*phi_j;

          }
        }
      }
    }
*/

/*
  switch(bc_name) {
    case EM_ER_FARFIELD_DIRECT_BC:
    case EM_EI_FARFIELD_DIRECT_BC:
    case EM_HR_FARFIELD_DIRECT_BC:
    case EM_HI_FARFIELD_DIRECT_BC:
      for (int g=0; g<pd->Num_Dim; g++){
        int gvar = var + g;
        for (int j=0; j<ei[pg->imtrx]->dof[gvar]; j++) {
          double phi_j = bf[gvar]->phi[j];
          for (int p=0; p<pd->Num_Dim; p++) {
            for (int q=0; q<pd->Num_Dim; q++) {
              //for (int r=0; r<pd->Num_Dim; r++) {
                d_func[p][gvar][j] += real*creal(permute(p,q,g)
                                                 *reduction_factor
                                                 *normal[q]
                                                 *phi_j)
                                   +  imag*cimag(permute(p,q,g)
                                                 *reduction_factor
                                                 *normal[q]
                                                 *phi_j);
              //}
            }
          }
        }
      }
      break;

  }

}*/

int apply_em_free_vec(double func[DIM],
                      double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                      double xi[DIM], /* Local stu coordinates */
                      const int bc_name) {

  double normal[DIM]; // surface normal vector

  // need Surface Normal vector
  for (int p = 0; p < DIM; p++) {
    normal[p] = fv->snormal[p];
  }

  // Evaluate n cross E or n cross H
  switch (bc_name) {
  case EM_ER_FREE_BC:
    for (int p = 0; p < DIM; p++) {
      for (int q = 0; q < DIM; q++) {
        for (int r = 0; r < DIM; r++) {
          func[p] += permute(p, q, r) * normal[q] * fv->em_er[r];
        }
      }
    }
    break;
  case EM_EI_FREE_BC:
    for (int p = 0; p < DIM; p++) {
      for (int q = 0; q < DIM; q++) {
        for (int r = 0; r < DIM; r++) {
          func[p] += permute(p, q, r) * normal[q] * fv->em_ei[r];
        }
      }
    }
    break;
  case EM_HR_FREE_BC:
    for (int p = 0; p < DIM; p++) {
      for (int q = 0; q < DIM; q++) {
        for (int r = 0; r < DIM; r++) {
          func[p] += permute(p, q, r) * normal[q] * fv->em_hr[r];
        }
      }
    }
    break;
  case EM_HI_FREE_BC:
    for (int p = 0; p < DIM; p++) {
      for (int q = 0; q < DIM; q++) {
        for (int r = 0; r < DIM; r++) {
          func[p] += permute(p, q, r) * normal[q] * fv->em_hi[r];
        }
      }
    }
    break;
  }
  return 0;
}

int apply_em_sommerfeld_vec(double func[DIM],
                            double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                            double xi[DIM], /* Local stu coordinates */
                            const int bc_name,
                            double *bc_data) {
  /***********************************************************************
   * TODO AMC: rewrite this description
   * apply_em_sommerfeld():
   *
   *  Function for specifying incoming plane wave and outgoing scattered
   *  energy using a radiative boundary condition sometimes called
   *  Sommerfeld radiation condition
   *
   *  n_bound CROSS E = j/kappa * Curl(E - E_i)
   *                  + n_bound CROSS E_i
   *
   *   func =
   *
   *
   *  The boundary condition EM_ employs this
   *  function.
   *
   *
   * Input:
   *
   *  em_eqn   = which equation is this applied to
   *  em_var   = which variable is this value sensitive to
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
   *   Author: Andrew Cochrane (10/9/2019)
   *
   ********************************************************************/
  dbl mag_permeability = 12.57e-07; // H/m

  double impedance = sqrt(mag_permeability / mp->permittivity);
  double omega = upd->EM_Frequency;
  double n[DIM]; // surface normal vector

  // double complex E1[DIM], H1[DIM]; // complex fields inside domain

  // This BC assumes that the boundary is far from the subject and
  // the material properties are the same on both sides
  // Need the wave number

  double kappa = omega * sqrt(mp->permittivity * mag_permeability);

  // need Surface Normal vector
  for (int p = 0; p < DIM; p++) {
    n[p] = fv->snormal[p];
  }

  // polarization P and propagation direction k of incident plane wave
  // from input deck
  complex double P[DIM];
  double k[DIM];

  P[0] = bc_data[0] + _Complex_I * bc_data[3];
  P[1] = bc_data[1] + _Complex_I * bc_data[4];
  P[2] = bc_data[2] + _Complex_I * bc_data[5];

  k[0] = bc_data[6];
  k[1] = bc_data[7];
  k[2] = bc_data[8];

  // normalize k and use wavenumber kappa
  double k_mag = sqrt(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]);

  k[0] = k[0] / k_mag;
  k[1] = k[1] / k_mag;
  k[2] = k[2] / k_mag;

  // Compute E_i
  // E_i = P*exp(i(k DOT x - wt)
  double x[DIM];
  x[0] = fv->x[0];
  x[1] = fv->x[1];
  x[2] = fv->x[2];

  double kappa_dot_x = 0.0;
  for (int q = 0; q < DIM; q++) {
    kappa_dot_x += k[q] * x[q];
  }
  kappa_dot_x *= kappa;

  complex double E_i[DIM] = {0.0};
  complex double CurlE_i[DIM] = {0.0};
  complex double nCrossE_i[DIM] = {0.0};
  complex double CurlE[DIM] = {0.0};
  complex double H_i[DIM] = {0.0};
  complex double CurlH_i[DIM] = {0.0};
  complex double nCrossH_i[DIM] = {0.0};
  complex double CurlH[DIM] = {0.0};

  for (int p = 0; p < DIM; p++) {
    E_i[p] += P[p] * cexp(_Complex_I * kappa_dot_x);
  }

  switch (bc_name) {
  case EM_ER_SOMMERFELD_BC:
  case EM_EI_SOMMERFELD_BC:

    // Compute Curl(E_i)
    // [Curl(E_i)]_e = sum_{f,g} [permute(e,f,g)*[d_dx]_f([E_i]_g)]
    for (int e = 0; e < DIM; e++) {
      for (int f = 0; f < DIM; f++) {
        for (int g = 0; g < DIM; g++) {
          CurlE_i[e] +=
              permute(e, f, g) * _Complex_I * kappa * k[f] * P[g] * cexp(_Complex_I * kappa_dot_x);
        }
      }
    }

    // Compute n CROSS E_i
    for (int e = 0; e < DIM; e++) {
      for (int f = 0; f < DIM; f++) {
        for (int g = 0; g < DIM; g++) {
          nCrossE_i[e] += permute(e, f, g) * n[f] * E_i[g];
        }
      }
    }

    // Compute Curl E
    for (int e = 0; e < DIM; e++) {
      for (int f = 0; f < DIM; f++) {
        for (int g = 0; g < DIM; g++) {
          CurlE[e] += permute(e, f, g) * (fv->grad_em_er[g][f] + _Complex_I * fv->grad_em_ei[g][f]);
        }
      }
    }
    break;

  case EM_HR_SOMMERFELD_BC:
  case EM_HI_SOMMERFELD_BC:
    for (int e = 0; e < DIM; e++) {
      for (int f = 0; f < DIM; f++) {
        for (int g = 0; g < DIM; g++) {
          H_i[e] += permute(e, f, g) * k[f] * E_i[g] / impedance;
        }
      }
    }

    // Compute Curl(H_i)
    // [Curl(H_i)]_e = sum_{f,g} [permute(e,f,g)*[d_dx]_f([H_i]_g)]
    for (int e = 0; e < DIM; e++) {
      for (int f = 0; f < DIM; f++) {
        for (int g = 0; g < DIM; g++) {
          CurlH_i[e] += permute(e, f, g) * _Complex_I * kappa * k[f] * H_i[g];
        }
      }
    }

    // Compute n CROSS H_i
    for (int e = 0; e < DIM; e++) {
      for (int f = 0; f < DIM; f++) {
        for (int g = 0; g < DIM; g++) {
          nCrossH_i[e] += permute(e, f, g) * n[f] * H_i[g];
        }
      }
    }

    // Compute Curl H
    for (int e = 0; e < DIM; e++) {
      for (int f = 0; f < DIM; f++) {
        for (int g = 0; g < DIM; g++) {
          CurlE[e] += permute(e, f, g) * (fv->grad_em_hr[g][f] + _Complex_I * fv->grad_em_hi[g][f]);
        }
      }
    }
    break;
  }

  // Residual Components

  switch (bc_name) {
  case EM_ER_SOMMERFELD_BC:
    for (int p = 0; p < DIM; p++) {
      func[p] -= creal(_Complex_I / kappa * (CurlE[p] - CurlE_i[p]) + nCrossE_i[p]);
    }
    break;
  case EM_EI_SOMMERFELD_BC:
    for (int p = 0; p < DIM; p++) {
      func[p] -= cimag(_Complex_I / kappa * (CurlE[p] - CurlE_i[p]) + nCrossE_i[p]);
    }
    break;
  case EM_HR_SOMMERFELD_BC:
    for (int p = 0; p < DIM; p++) {
      func[p] -= creal(_Complex_I / kappa * (CurlH[p] - CurlH_i[p]) + nCrossH_i[p]);
    }
    break;
  case EM_HI_SOMMERFELD_BC:
    for (int p = 0; p < DIM; p++) {
      func[p] -= cimag(_Complex_I / kappa * (CurlH[p] - CurlH_i[p]) + nCrossH_i[p]);
    }
    break;
  }
  return 0;
} // end of apply_em_sommerfeld_vec

int apply_ewave_planewave_vec(double func[DIM],
                              double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                              double xi[DIM], /* Local stu coordinates */
                              const int bc_name,
                              double *bc_data) {
  /***********************************************************************
   * TODO AMC: rewrite this description + add Jacobians
   * apply_ewave_planewave_vec():
   *
   *  Function for specifying surface integral terms for plane wave test.
   *
   *  Surface Integral = n_bound CROSS Curl E_i dS
   *
   *   func =
   *
   *
   *  The boundary condition E_PLANEWAVE employs this
   *  function.
   *
   *
   * Input:
   *
   *  em_eqn   = which equation is this applied to
   *  em_var   = which variable is this value sensitive to
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
   *   Author: Andrew Cochrane (5/19/2020)
   *
   ********************************************************************/

  // double impedance = sqrt(mp->magnetic_permeability/mp->permittivity);
  double omega = upd->EM_Frequency;
  double n[DIM]; // surface normal vector

  // double complex E1[DIM], H1[DIM]; // complex fields inside domain

  // This BC assumes that the boundary is far from the subject and
  // the material properties are the same on both sides
  // Need the wave number

  double kappa = omega * sqrt(mp->permittivity * mp->magnetic_permeability);

  // need Surface Normal vector
  for (int p = 0; p < DIM; p++) {
    n[p] = fv->snormal[p];
  }

  // polarization P, propagation direction k of incident plane wave
  // as well as origin offset from input deck
  complex double P[DIM];
  double k[DIM];
  double x_orig[DIM];

  P[0] = bc_data[0] + _Complex_I * bc_data[3];
  P[1] = bc_data[1] + _Complex_I * bc_data[4];
  P[2] = bc_data[2] + _Complex_I * bc_data[5];

  k[0] = bc_data[6];
  k[1] = bc_data[7];
  k[2] = bc_data[8];

  x_orig[0] = bc_data[9];
  x_orig[1] = bc_data[10];
  x_orig[2] = bc_data[11];

  // normalize k and use wavenumber kappa
  double k_mag = sqrt(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]);

  k[0] = k[0] / k_mag;
  k[1] = k[1] / k_mag;
  k[2] = k[2] / k_mag;

  double x[DIM];
  x[0] = fv->x[0] - x_orig[0];
  x[1] = fv->x[1] - x_orig[1];
  x[2] = fv->x[2] - x_orig[2];

  double kappa_dot_x = 0.0;
  for (int q = 0; q < DIM; q++) {
    kappa_dot_x += k[q] * x[q];
  }
  kappa_dot_x *= kappa;

  complex double E[DIM] = {0.0};
  complex double CurlE[DIM] = {0.0};
  complex double nCrossCurlE[DIM] = {0.0};

  for (int p = 0; p < DIM; p++) {
    E[p] += P[p] * cexp(_Complex_I * kappa_dot_x);
  }

  // Compute Curl(E)
  // [Curl(E)]_e = sum_{f,g} [permute(e,f,g)*[d_dx]_f([E]_g)]
  for (int e = 0; e < DIM; e++) {
    for (int f = 0; f < DIM; f++) {
      for (int g = 0; g < DIM; g++) {
        CurlE[e] += permute(e, f, g) * kappa * k[f] * P[g] * cexp(_Complex_I * kappa_dot_x);
      }
    }
  }

  // Compute n CROSS Curl(E)
  for (int e = 0; e < DIM; e++) {
    for (int f = 0; f < DIM; f++) {
      for (int g = 0; g < DIM; g++) {
        nCrossCurlE[e] += permute(e, f, g) * n[f] * CurlE[g];
      }
    }
  }

  // Residual Components

  switch (bc_name) {
  case E_ER_PLANEWAVE_BC:
    for (int p = 0; p < DIM; p++) {
      func[p] += creal(nCrossCurlE[p]);
    }
    break;
  case E_EI_PLANEWAVE_BC:
    for (int p = 0; p < DIM; p++) {
      func[p] += cimag(nCrossCurlE[p]);
    }
    break;
  }
  return 0;
} // end of apply_ewave_planewave_vec

int apply_ewave_curlcurl_farfield_vec(double func[DIM],
                                      double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                                      double xi[DIM], /* Local stu coordinates */
                                      double time,    // present time
                                      const int bc_name,
                                      double *bc_data) {
  /***********************************************************************
   * TODO AMC: rewrite this description + add Jacobians
   * apply_ewave_curlcurl_farfield_vec():
   *
   *  Function for specifying surface integral terms for curlcurlE eqn.
   *
   *  Surface Integral = phi dot n_bound CROSS Curl E_1 dS
   *
   *   Derived Expression hasn't accounted for difference in polarization
   *   between E1 and thefarfield incident wave - for now assume they're
   *   the same and in the direction of the farfield incident wave
   *
   *
   *
   *   func = psi dot [ -i kappa_1 (E_1 - 2*E_f) eta_1/eta_2 ]
   *
   *
   *  The boundary condition E_FARFIELD employs this
   *  function.
   *
   *
   * Input:
   *
   *  em_eqn   = which equation is this applied to
   *  em_var   = which variable is this value sensitive to
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
   *   Author: Andrew Cochrane (9/21/2020)
   *
   ********************************************************************/

  // need Surface Normal vector
  /* maybe not
    double n_surf[DIM]; // surface normal vector
    for (int p=0; p<DIM; p++) {
      n_surf[p] = fv->snormal[p];
    }
  */
  // setup wave and material properties
  double omega = upd->EM_Frequency;
  dbl n_1; /* Refractive index. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_n_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_n = &d_n_struct;

  dbl k_1; /* extinction coefficient. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  n_1 = refractive_index(d_n, time);
  k_1 = extinction_index(d_k, time);

  // Compute complex material properties
  complex double cpx_refractive_index_1, cpx_rel_permittivity_1,
      cpx_permittivity_1; //, impedance;

  cpx_refractive_index_1 = n_1 + _Complex_I * k_1; // k > 0 is extinction
  cpx_rel_permittivity_1 = SQUARE(cpx_refractive_index_1);
  cpx_permittivity_1 =
      cpx_rel_permittivity_1 * mp->permittivity; // better set permittivity to vacuum in input deck?

  double complex E_1[DIM]; // complex field inside domain

  for (int p = 0; p < DIM; p++) {
    E_1[p] = fv->em_er[p] + _Complex_I * fv->em_ei[p];
  }

  // Need the wave number
  complex double kappa_1 = omega * csqrt(cpx_permittivity_1 * mp->magnetic_permeability);
  // Need the impedance
  complex double eta_1 = csqrt(mp->magnetic_permeability / cpx_permittivity_1);

  // polarization P of incident plane wave
  complex double P[DIM];

  // permittivity of the far-field medium
  double n_2, k_2;

  P[0] = bc_data[0] + _Complex_I * bc_data[3];
  P[1] = bc_data[1] + _Complex_I * bc_data[4];
  P[2] = bc_data[2] + _Complex_I * bc_data[5];

  // Need the impedance of the far-field material (superstrate)
  // use refractive index and extinction coefficient
  n_2 = bc_data[6];
  k_2 = bc_data[7];

  // Compute complex material properties
  complex double cpx_rel_permittivity_2, cpx_refractive_index_2, cpx_permittivity_2;

  cpx_refractive_index_2 = n_2 + _Complex_I * k_2; // k > 0 is extinction
  cpx_rel_permittivity_2 = SQUARE(cpx_refractive_index_2);
  cpx_permittivity_2 =
      cpx_rel_permittivity_2 * mp->permittivity; // better set permittivity to vacuum in input deck?
  // Need the impedance
  complex double eta_2 = csqrt(mp->magnetic_permeability / cpx_permittivity_2);

  complex double cpx_coeff = -_Complex_I * kappa_1 * eta_1 / eta_2;

  // Residual Components

  switch (bc_name) {
  case E_ER_FARFIELD_BC:

    for (int p = 0; p < pd->Num_Dim; p++) {
      func[p] = creal(cpx_coeff * (E_1[p] - 2 * P[p]));
      for (int q = 0; q < pd->Num_Dim; q++) {
        int qRvar = EM_E1_REAL + q;
        int qIvar = EM_E1_IMAG + q;
        for (int j = 0; j < ei[pg->imtrx]->dof[EM_E1_REAL]; j++) {
          double Rphi_j = bf[EM_E1_REAL + q]->phi[j];
          double Iphi_j = bf[EM_E1_IMAG + q]->phi[j];
          // d(Re[f{z}])/d(Re[z]) = Re(f'(z))
          d_func[p][qRvar][j] = delta(p, q) * creal(cpx_coeff * Rphi_j);

          // d(Re[f{z}])/d(Im[z]) = -Im(f'(z))
          d_func[p][qIvar][j] = -delta(p, q) * cimag(cpx_coeff * Iphi_j);
        }
      }
    }
    /*/
  for (int p=0; p<DIM; p++) {
      //func[p] += creal((E_1[p]));
    //func[p] = fv->em_er[p];
    func[p] = 1.0;
      for (int q=0; q<pd->Num_Dim; q++) {
        int qRvar = EM_E1_REAL + q;
        int qIvar = EM_E1_IMAG + q;
        for (int j=0; j<ei[pg->imtrx]->dof[EM_E1R_BC]; j++) {
          double Rphi_j = bf[EM_E1_REAL+q]->phi[j];
          double Iphi_j = bf[EM_E1_IMAG+q]->phi[j];
          // d(Re[f{z}])/d(Re[z]) = Re(f'(z))
          //d_func[p][qRvar][j] = delta(p,q)*Rphi_j;
          d_func[p][qRvar][j] = 0.0;
              // *creal(Rphi_j);

          // d(Re[f{z}])/d(Im[z]) = -Im(f'(z))
          d_func[p][qIvar][j] += -delta(p,q)
              *0.0;
              // *cimag(Iphi_j);
        }
      }
    }
  */
    break;
  case E_EI_FARFIELD_BC:

    for (int p = 0; p < DIM; p++) {
      func[p] += cimag(cpx_coeff * (E_1[p] - 2 * P[p]));
      for (int q = 0; q < pd->Num_Dim; q++) {
        int qRvar = EM_E1_REAL + q;
        int qIvar = EM_E1_IMAG + q;
        for (int j = 0; j < ei[pg->imtrx]->dof[EM_E1_IMAG]; j++) {
          double Rphi_j = bf[EM_E1_REAL + q]->phi[j];
          double Iphi_j = bf[EM_E1_IMAG + q]->phi[j];

          // d(Im[f{z}])/d(Re[z]) = Im(f'(z))
          d_func[p][qRvar][j] = delta(p, q) * cimag(cpx_coeff * Rphi_j);

          // d(Im[f{z}])/d(Im[z]) = Re(f'(z))
          d_func[p][qIvar][j] = delta(p, q) * creal(cpx_coeff * Iphi_j);
        }
      }
    }
    /*/
    for (int p=0; p<pd->Num_Dim; p++) {
      func[p] += cimag(cpx_coeff*(E_1[p]));
      //func[p] = fv->em_ei[p];
        for (int q=0; q<pd->Num_Dim; q++) {
          int qRvar = EM_E1_REAL + q;
          int qIvar = EM_E1_IMAG + q;
          for (int j=0; j<ei[pg->imtrx]->dof[EM_E1_IMAG]; j++) {
            double Rphi_j = bf[EM_E1_REAL+q]->phi[j];
            double Iphi_j = bf[EM_E1_IMAG+q]->phi[j];
            // d(Im[f{z}])/d(Re[z]) = Im(f'(z))
            d_func[p][qRvar][j] = delta(p,q)
                *cimag(cpx_coeff*Rphi_j);

            // d(Im[f{z}])/d(Im[z]) = Re(f'(z))
            d_func[p][qIvar][j] = delta(p,q)
                *creal(cpx_coeff*Iphi_j);
          }
        }
      }
*/
    break;
  }

  return 0;
} // end of apply_ewave_farfield_vec

int apply_ewave_mms_vec(double func[DIM],
                        double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                        double xi[DIM], /* Local stu coordinates */
                        const int bc_name,
                        double *bc_data) {
  /***********************************************************************
   * TODO AMC: rewrite this description + add Jacobians
   * apply_ewave_planewave_vec():
   *
   *  Function for specifying surface integral terms for plane wave test.
   *
   *  Surface Integral = n_bound CROSS Curl E_mms dS
   *
   *   func =
   *
   *
   *  The boundary condition E_MMS employs this
   *  function.
   *
   *
   * Input:
   *
   *  em_eqn   = which equation is this applied to
   *  em_var   = which variable is this value sensitive to
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
   *   Author: Andrew Cochrane (7/15/2020)
   *
   ********************************************************************/

  double omega = upd->EM_Frequency;
  double n[DIM]; // surface normal vector

  // double complex E1[DIM], H1[DIM]; // complex fields inside domain

  // This BC assumes that the boundary is far from the subject and
  // the material properties are the same on both sides
  // Need the wave number

  double kappa = omega * sqrt(mp->permittivity * mp->magnetic_permeability);

  // need Surface Normal vector
  for (int p = 0; p < DIM; p++) {
    n[p] = fv->snormal[p];
  }

  // polarization P, propagation direction k of incident plane wave
  // as well as origin offset from input deck
  complex double P[DIM];
  double k[DIM];
  double x_orig[DIM];

  P[0] = bc_data[0] + _Complex_I * bc_data[3];
  P[1] = bc_data[1] + _Complex_I * bc_data[4];
  P[2] = bc_data[2] + _Complex_I * bc_data[5];

  k[0] = bc_data[6];
  k[1] = bc_data[7];
  k[2] = bc_data[8];

  x_orig[0] = bc_data[9];
  x_orig[1] = bc_data[10];
  x_orig[2] = bc_data[11];

  // normalize k and use wavenumber kappa
  double k_mag = sqrt(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]);

  k[0] = k[0] / k_mag;
  k[1] = k[1] / k_mag;
  k[2] = k[2] / k_mag;

  double x[DIM];
  x[0] = fv->x[0] - x_orig[0];
  x[1] = fv->x[1] - x_orig[1];
  x[2] = fv->x[2] - x_orig[2];

  x[0] = fv->x[0];
  x[1] = fv->x[1];
  x[2] = fv->x[2];

  double kappa_dot_x = 0.0;
  for (int q = 0; q < DIM; q++) {
    kappa_dot_x += k[q] * x[q];
  }
  kappa_dot_x *= kappa;

  complex double E[DIM] = {0.0};
  complex double CurlE[DIM] = {0.0};
  complex double nCrossCurlE[DIM] = {0.0};

  for (int p = 0; p < DIM; p++) {
    E[p] += P[p] * cexp(_Complex_I * kappa_dot_x);
  }

  // Compute Curl(E)
  // [Curl(E)]_e = sum_{f,g} [permute(e,f,g)*[d_dx]_f([E]_g)]
  for (int e = 0; e < DIM; e++) {
    for (int f = 0; f < DIM; f++) {
      for (int g = 0; g < DIM; g++) {
        CurlE[e] += permute(e, f, g) * kappa * k[f] * P[g] * cexp(_Complex_I * kappa_dot_x);
      }
    }
  }

  // Compute n CROSS Curl(E)
  for (int e = 0; e < DIM; e++) {
    for (int f = 0; f < DIM; f++) {
      for (int g = 0; g < DIM; g++) {
        nCrossCurlE[e] += permute(e, f, g) * n[f] * CurlE[g];
      }
    }
  }

  // Residual Components

  switch (bc_name) {
  case E_ER_PLANEWAVE_BC:
    for (int p = 0; p < DIM; p++) {
      func[p] += creal(nCrossCurlE[p]);
    }
    break;
  case E_EI_PLANEWAVE_BC:
    for (int p = 0; p < DIM; p++) {
      func[p] += cimag(nCrossCurlE[p]);
    }
    break;
  }
  return 0;
} // end of apply_ewave_mms_vec

int apply_ewave_2D(double func[DIM],
                   double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                   double xi[DIM], /* Local stu coordinates */
                   const int bc_name) {
  /***********************************************************************
   * TODO AMC: rewrite this description + add Jacobians
   * apply_ewave_2D():
   *
   *  strongly integrate the dEz_dz = 0
   *
   *  Surface Integral = d/dx_3( x_3 dot E) = 0 dS
   *
   *   func_r[0,1,2] = fv->grad_em_er[0,1,2][2]
   *   func_i[0,1,2] = fv->grad_em_ei[0,1,2][2]
   *
   *  The boundary condition E_PLANEWAVE employs this
   *  function.
   *
   *
   * Input:
   *
   *
   *
   * Output:
   *
   *  func[0,1,2] = value of the function mentioned above
   *  d_func[0][varType][lvardof] =
   *              Derivate of func[0] wrt
   *              the variable type, varType, and the local variable
   *              degree of freedom, lvardof, corresponding to that
   *              variable type.
   *
   *   Author: Andrew Cochrane (6/22/2020)
   *
   ********************************************************************/
  int var;
  double dphi_j_dx3;
  // Residual Components

  switch (bc_name) {
  case E_ER_2D_BC:
    // func[0] += fv->grad_em_er[2][0];
    // func[1] += fv->grad_em_er[2][1];
    func[2] += fv->grad_em_er[2][2];
    break;
  case E_EI_2D_BC:
    // func[0] += fv->grad_em_ei[2][0];
    // func[1] += fv->grad_em_ei[2][1];
    func[2] += fv->grad_em_ei[2][2];
    break;
  }

  // Jacobian Components
  if (af->Assemble_Jacobian) {
    switch (bc_name) {
    case E_ER_2D_BC:
      for (int k = 2; k < DIM; k++) {
        var = EM_E1_REAL + k;
        for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          dphi_j_dx3 = bf[var]->grad_phi[j][2];
          d_func[k][var][j] += dphi_j_dx3;
        }
      }
      break;
    case E_EI_2D_BC:
      for (int k = 2; k < DIM; k++) {
        var = EM_E1_IMAG + k;
        for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          dphi_j_dx3 = bf[var]->grad_phi[j][2];
          d_func[k][var][j] += dphi_j_dx3;
        }
      }
      break;
    }
  }
  return 0;
} // end of apply_ewave_2D
int apply_ewave_nedelec_farfield(double func[DIM],
                                 double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                                 double xi[DIM], /* Local stu coordinates */
                                 double time,    // present time
                                 const int bc_name,
                                 double *bc_data) {
  /***********************************************************************
   * TODO AMC: rewrite this description + add Jacobians
   * apply_ewave_curlcurl_farfield_vec():
   *
   *  Function for specifying surface integral terms for curlcurlE eqn.
   *
   *  Surface Integral = phi dot n_bound CROSS Curl E_1 dS
   *
   *   Derived Expression hasn't accounted for difference in polarization
   *   between E1 and thefarfield incident wave - for now assume they're
   *   the same and in the direction of the farfield incident wave
   *
   *
   *
   *   func = psi dot [ -i kappa_1 (E_1 - 2*E_f) eta_1/eta_2 ]
   *
   *
   *  The boundary condition E_FARFIELD employs this
   *  function.
   *
   *
   * Input:
   *
   *  em_eqn   = which equation is this applied to
   *  em_var   = which variable is this value sensitive to
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
   *   Author: Andrew Cochrane (9/21/2020)
   *
   ********************************************************************/

  // need Surface Normal vector
  /* maybe not
    double n_surf[DIM]; // surface normal vector
    for (int p=0; p<DIM; p++) {
      n_surf[p] = fv->snormal[p];
    }
  */
  // setup wave and material properties
  double omega = upd->EM_Frequency;
  dbl n_1; /* Refractive index. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_n_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_n = &d_n_struct;

  dbl k_1; /* extinction coefficient. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  n_1 = refractive_index(d_n, time);
  k_1 = extinction_index(d_k, time);

  // Compute complex material properties
  complex double cpx_refractive_index_1, cpx_rel_permittivity_1,
      cpx_permittivity_1; //, impedance;

  cpx_refractive_index_1 = n_1 + _Complex_I * k_1; // k > 0 is extinction
  cpx_rel_permittivity_1 = SQUARE(cpx_refractive_index_1);
  cpx_permittivity_1 =
      cpx_rel_permittivity_1 * mp->permittivity; // better set permittivity to vacuum in input deck?

  double complex E_1[DIM]; // complex field inside domain

  for (int p = 0; p < DIM; p++) {
    E_1[p] = fv->em_er[p] + _Complex_I * fv->em_ei[p];
  }

  // Need the wave number
  complex double kappa_1 = omega * csqrt(cpx_permittivity_1 * mp->magnetic_permeability);
  // Need the impedance
  complex double eta_1 = csqrt(mp->magnetic_permeability / cpx_permittivity_1);

  // polarization P of incident plane wave
  complex double P[DIM];

  // permittivity of the far-field medium
  double n_2, k_2;

  P[0] = bc_data[0] + _Complex_I * bc_data[3];
  P[1] = bc_data[1] + _Complex_I * bc_data[4];
  P[2] = bc_data[2] + _Complex_I * bc_data[5];

  // Need the impedance of the far-field material (superstrate)
  // use refractive index and extinction coefficient
  n_2 = bc_data[6];
  k_2 = bc_data[7];

  // Compute complex material properties
  complex double cpx_rel_permittivity_2, cpx_refractive_index_2, cpx_permittivity_2;

  cpx_refractive_index_2 = n_2 + _Complex_I * k_2; // k > 0 is extinction
  cpx_rel_permittivity_2 = SQUARE(cpx_refractive_index_2);
  cpx_permittivity_2 =
      cpx_rel_permittivity_2 * mp->permittivity; // better set permittivity to vacuum in input deck?
  // Need the impedance
  complex double eta_2 = csqrt(mp->magnetic_permeability / cpx_permittivity_2);

  complex double cpx_coeff = -_Complex_I * kappa_1 * eta_1 / eta_2;

  // Residual Components

  switch (bc_name) {
  case EM_FARFIELD_REAL_NED_BC:

    for (int p = 0; p < pd->Num_Dim; p++) {
      func[p] = creal(cpx_coeff * (E_1[p] - 2 * P[p]));
      for (int j = 0; j < ei[pg->imtrx]->dof[EM_E1_REAL]; j++) {
        double Rphi_j = bf[EM_E1_REAL]->phi_e[j][p];
        double Iphi_j = bf[EM_E1_IMAG]->phi_e[j][p];
        // d(Re[f{z}])/d(Re[z]) = Re(f'(z))
        d_func[p][EM_E1_REAL][j] = creal(cpx_coeff * Rphi_j);

        // d(Re[f{z}])/d(Im[z]) = -Im(f'(z))
        d_func[p][EM_E1_IMAG][j] = -cimag(cpx_coeff * Iphi_j);
      }
    }
    break;
  case EM_FARFIELD_IMAG_NED_BC:

    for (int p = 0; p < DIM; p++) {
      func[p] += cimag(cpx_coeff * (E_1[p] - 2 * P[p]));
      for (int j = 0; j < ei[pg->imtrx]->dof[EM_E1_IMAG]; j++) {
        double Rphi_j = bf[EM_E1_REAL]->phi_e[j][p];
        double Iphi_j = bf[EM_E1_IMAG]->phi_e[j][p];
        d_func[p][EM_E1_REAL][j] = cimag(cpx_coeff * Rphi_j);
        d_func[p][EM_E1_IMAG][j] = creal(cpx_coeff * Iphi_j);
      }
    }
    break;
  }

  return 0;
} // end of apply_ewave_farfield_vec
void em_absorbing_bc_nedelec(int bc_name,
                             dbl *func,
                             double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE]) {
  const double c0 = 1.0 / (sqrt(upd->Free_Space_Permittivity * upd->Free_Space_Permeability));
  ;
  dbl x = fv->x[0];
  dbl y = fv->x[1];
  dbl z = fv->x[2];
  dbl freq = upd->EM_Frequency;
  dbl lambda0 = c0 / freq;
  dbl k0 = 2 * M_PI / lambda0;
  complex double wave[DIM];
  complex double curl_wave[DIM];
  incident_wave(x, y, z, k0, wave, curl_wave);
  complex double j = _Complex_I;
  complex double E[DIM] = {fv->em_er[0] + j * fv->em_ei[0], fv->em_er[1] + j * fv->em_ei[1],
                           fv->em_er[2] + j * fv->em_ei[2]};
  switch (bc_name) {
  case EM_ABSORBING_REAL_BC: {
    double wave_r[DIM] = {creal(wave[0]), creal(wave[1]), creal(wave[2])};
    double wave_i[DIM] = {cimag(wave[0]), cimag(wave[1]), cimag(wave[2])};

    double curl_wave_r[DIM] = {creal(curl_wave[0]), creal(curl_wave[1]), creal(curl_wave[2])};
    double curl_wave_i[DIM] = {cimag(curl_wave[0]), cimag(curl_wave[1]), cimag(curl_wave[2])};

    double nxWave_r[DIM];
    double nxWave_i[DIM];

    cross_really_simple_vectors(fv->snormal, wave_r, nxWave_r);
    cross_really_simple_vectors(fv->snormal, wave_i, nxWave_i);

    double nxnxWave_r[DIM];
    double nxnxWave_i[DIM];

    cross_really_simple_vectors(fv->snormal, nxWave_r, nxnxWave_r);
    cross_really_simple_vectors(fv->snormal, nxWave_i, nxnxWave_i);

    dbl ke = k0;

    double nxCurlWave_r[DIM];
    double nxCurlWave_i[DIM];
    cross_really_simple_vectors(fv->snormal, curl_wave_r, nxCurlWave_r);
    cross_really_simple_vectors(fv->snormal, curl_wave_i, nxCurlWave_i);

    complex double Uinc[DIM];
    for (int i = 0; i < DIM; i++) {
      Uinc[i] =
          nxCurlWave_r[i] + nxCurlWave_i[i] * j + j * k0 * (nxnxWave_r[i] + nxnxWave_i[i] * j);
    }
    double ABC[DIM];

    for (int i = 0; i < DIM; i++) {
      ABC[i] = creal(j * ke * E[i]);
    }

    dbl nxABC[DIM] = {0.0};
    cross_really_simple_vectors(fv->snormal, ABC, nxABC);

    for (int i = 0; i < DIM; i++) {
      func[i] += +(ABC[i] + creal(Uinc[i]));
    }

    for (int j = 0; j < ei[pg->imtrx]->dof[EM_E1_IMAG]; j++) {
      dbl Ei_imag[DIM] = {-ke * bf[EM_E1_IMAG]->phi_e[j][0], -ke * bf[EM_E1_IMAG]->phi_e[j][1],
                          -ke * bf[EM_E1_IMAG]->phi_e[j][2]};

      d_func[0][EM_E1_IMAG][j] += Ei_imag[0];
      d_func[1][EM_E1_IMAG][j] += Ei_imag[1];
      d_func[2][EM_E1_IMAG][j] += Ei_imag[2];
    }
  } break;
  case EM_ABSORBING_IMAG_BC: {
    double wave_r[DIM] = {creal(wave[0]), creal(wave[1]), creal(wave[2])};
    double wave_i[DIM] = {cimag(wave[0]), cimag(wave[1]), cimag(wave[2])};

    double curl_wave_r[DIM] = {creal(curl_wave[0]), creal(curl_wave[1]), creal(curl_wave[2])};
    double curl_wave_i[DIM] = {cimag(curl_wave[0]), cimag(curl_wave[1]), cimag(curl_wave[2])};

    double nxWave_r[DIM];
    double nxWave_i[DIM];

    cross_really_simple_vectors(fv->snormal, wave_r, nxWave_r);
    cross_really_simple_vectors(fv->snormal, wave_i, nxWave_i);

    double nxnxWave_r[DIM];
    double nxnxWave_i[DIM];

    cross_really_simple_vectors(fv->snormal, nxWave_r, nxnxWave_r);
    cross_really_simple_vectors(fv->snormal, nxWave_i, nxnxWave_i);

    dbl ke = k0;

    double nxCurlWave_r[DIM];
    double nxCurlWave_i[DIM];
    cross_really_simple_vectors(fv->snormal, curl_wave_r, nxCurlWave_r);
    cross_really_simple_vectors(fv->snormal, curl_wave_i, nxCurlWave_i);

    complex double Uinc[DIM];
    for (int i = 0; i < DIM; i++) {
      Uinc[i] =
          nxCurlWave_r[i] + nxCurlWave_i[i] * j + j * k0 * (nxnxWave_r[i] + nxnxWave_i[i] * j);
    }
    dbl ABC[DIM];

    for (int i = 0; i < DIM; i++) {
      ABC[i] = cimag(j * ke * E[i]);
    }
    dbl nxABC[DIM] = {0.0};
    cross_really_simple_vectors(fv->snormal, ABC, nxABC);

    for (int i = 0; i < DIM; i++) {
      func[i] = +(ABC[i] + cimag(Uinc[i]));
    }

    for (int j = 0; j < ei[pg->imtrx]->dof[EM_E1_REAL]; j++) {
      dbl Ei_real[DIM] = {ke * bf[EM_E1_REAL]->phi_e[j][0], ke * bf[EM_E1_REAL]->phi_e[j][1],
                          ke * bf[EM_E1_REAL]->phi_e[j][2]};
      dbl deriv[DIM] = {0.0};
      cross_really_simple_vectors(fv->snormal, Ei_real, deriv);
      for (int i = 0; i < pd->Num_Dim; i++) {
        d_func[i][EM_E1_REAL][j] += Ei_real[i];
      }
    }
  } break;
  default:
    GOMA_EH(GOMA_ERROR, "Unknown bc type %d for em_absorbing_bc_nedelec", bc_name);
    break;
  }
}

// This was used for testing but doesn't work at the moment
// Appraoch for applying a Dirichlet condition for Nedelec EM solves
void em_mms_nedelec_bc(int bc_name,
                       dbl *func,
                       double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE]) {
  switch (bc_name) {
  case EM_MMS_SIDE_BC: {
    complex double exact[DIM];
    em_mms_exact(fv->x[0], fv->x[1], fv->x[2], exact);
    dbl omega = upd->EM_Frequency;
    dbl x = fv->x[0];
    dbl y = fv->x[1];
    dbl z = fv->x[2];
    complex double wave[DIM];
    complex double curl_wave[DIM];
    incident_wave(x, y, z, omega, wave, curl_wave);
    dbl exact_real[DIM];
    for (int i = 0; i < pd->Num_Dim; i++) {
      exact_real[i] = 0 * creal(wave[i]);
    }

    for (int i = 0; i < pd->Num_Dim; i++) {
      exact_real[i] += fv->em_er[i];
    }
    cross_really_simple_vectors(fv->snormal, exact_real, func);
    for (int j = 0; j < ei[pg->imtrx]->dof[EM_E1_REAL]; j++) {
      for (int i = 0; i < pd->Num_Dim; i++) {
        exact_real[i] = bf[EM_E1_REAL]->phi_e[j][i];
      }
      dbl d_exact[DIM];
      cross_really_simple_vectors(fv->snormal, exact_real, d_exact);
      for (int i = 0; i < pd->Num_Dim; i++) {
        d_func[i][EM_E1_REAL][j] = d_exact[i];
      }
    }
  } break;
  case EM_MMS_SIDE_IMAG_BC: {
    dbl omega = upd->EM_Frequency;
    dbl x = fv->x[0];
    dbl y = fv->x[1];
    dbl z = fv->x[2];
    complex double wave[DIM];
    complex double curl_wave[DIM];
    incident_wave(x, y, z, omega, wave, curl_wave);
    complex double exact[DIM];
    em_mms_exact(fv->x[0], fv->x[1], fv->x[2], exact);
    dbl exact_imag[DIM];
    for (int i = 0; i < pd->Num_Dim; i++) {
      exact_imag[i] = 0 * cimag(wave[i]);
    }

    for (int i = 0; i < pd->Num_Dim; i++) {
      exact_imag[i] += fv->em_ei[i];
    }
    cross_really_simple_vectors(fv->snormal, exact_imag, func);
    for (int j = 0; j < ei[pg->imtrx]->dof[EM_E1_IMAG]; j++) {
      for (int i = 0; i < pd->Num_Dim; i++) {
        exact_imag[i] = bf[EM_E1_IMAG]->phi_e[j][i];
      }
      dbl d_exact[DIM] = {0.0};
      cross_really_simple_vectors(fv->snormal, exact_imag, d_exact);
      for (int i = 0; i < pd->Num_Dim; i++) {
        d_func[i][EM_E1_IMAG][j] = d_exact[i];
      }
    }
  } break;
  }
}