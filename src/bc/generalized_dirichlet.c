/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2025 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/

#include "bc/generalized_dirichlet.h"
#include "el_elm.h"
#include "load_variable.h"
#include "mm_as.h"
#include "rf_bc.h"
#include "std.h"
#include "table.h"

/******************************************************************************

 fgeneralized_dirichlet() - simple pointwise collocation contribution for BCs

     Function which adds simple contributions to boundary conditions.
     This is used as a pointwise collocation boundary condition, and has
     flexibility built in so that the condition can be applied to any equation
     and be sensitive with respect to any unknown

     Author:          Rich Cairncross (1511)
     Date:            22 December 1994
     Revised:
******************************************************************************/

/*
 * There are inconsistencies in the prototype declaration and how this function
 * is invoked in bc_curve.c (circa line 1170).
 *
 * Based on other boundary conditions, I believe the declaration needs to be
 * more like:
 *
 *	double func[],
 *      double d_func[][MAX_VARIABLE_TYPES+MAX_CONC][MDE],
 *      ...
 * Resolve this at some point. -PAS
 */
int fgeneralized_dirichlet(double *func,
                           double d_func[],        /* MAX_VARIABLE_TYPES + MAX_CONC */
                           const int gd_condition, /* denoting which condition
                                                    * applied */
                           const int bc_input_id,
                           const double tt, /* parameter to vary time integration
                                             * from explicit (tt = 1) to
                                             * implicit (tt = 0) */
                           const double dt) /* current time step size          */
{
  int jvar, wspec, vector_sens, b;
  int index_var;       /* Column index into the global stiffness matrix*/
  dbl x_var;           /* value of variable at this node */
  dbl d_x_var;         /* sensitivity of variable to nodal unknown */
  dbl d_vect_var[DIM]; /* sensitivity of vector variable to nodal unknown */
  dbl slope;           /* slope of interpolated function in table */
  dbl x_var_mp[1];     /* dummy variable for table lookup subroutines */

  if (af->Assemble_LSA_Mass_Matrix)
    return 0;

  /* ---- Find variable number and species number */

  jvar = BC_Types[bc_input_id].BC_Data_Int[2];

  wspec = BC_Types[bc_input_id].BC_Data_Int[3];

  /* put value of variable in GD Condition into x_var and put it's sensitivity in d_x_var */
  index_var = load_variable(&x_var, &d_x_var, jvar, wspec, tt, dt, d_vect_var);

  if (jvar == SPEED) {
    vector_sens = 1;
  } else {
    vector_sens = 0;
  }
  /* Now add in contributions to residual vector and jacobian matrix */

  switch (gd_condition) {
  case (GD_CONST_BC): /* x - c0 */

    *func = (x_var - BC_Types[bc_input_id].BC_Data_Float[0]);

    if (af->Assemble_Jacobian) {
      if (vector_sens) {
        for (b = 0; b < DIM; b++) {
          d_func[index_var + b] = d_vect_var[b];
        }
      } else {
        d_func[index_var] = d_x_var;
      }
    }
    break;

  case (GD_LINEAR_BC): /* C1 x + c0 */

    *func =
        (x_var * BC_Types[bc_input_id].BC_Data_Float[1] + BC_Types[bc_input_id].BC_Data_Float[0]);

    if (af->Assemble_Jacobian) {
      if (vector_sens) {
        for (b = 0; b < DIM; b++) {
          d_func[index_var + b] = d_vect_var[b] * BC_Types[bc_input_id].BC_Data_Float[1];
        }
      } else {
        d_func[index_var] = d_x_var * BC_Types[bc_input_id].BC_Data_Float[1];
      }
    }
    break;

  case (GD_INVERSE_BC): /* C1/x + c0 */

    *func =
        (BC_Types[bc_input_id].BC_Data_Float[1] / x_var + BC_Types[bc_input_id].BC_Data_Float[0]);

    if (af->Assemble_Jacobian) {
      if (vector_sens) {
        for (b = 0; b < DIM; b++) {
          d_func[index_var + b] =
              -d_vect_var[b] * BC_Types[bc_input_id].BC_Data_Float[1] / (x_var * x_var);
        }
      } else {
        d_func[index_var] = -d_x_var * BC_Types[bc_input_id].BC_Data_Float[1] / (x_var * x_var);
      }
    }
    break;

  case (GD_PARAB_BC): /* C2 x^2 + C1 x + c0 */

    *func =
        (x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[2] +
         x_var * BC_Types[bc_input_id].BC_Data_Float[1] + BC_Types[bc_input_id].BC_Data_Float[0]);

    if (af->Assemble_Jacobian) {
      if (vector_sens) {
        for (b = 0; b < DIM; b++) {
          d_func[index_var + b] =
              d_vect_var[b] * (BC_Types[bc_input_id].BC_Data_Float[1] +
                               2. * x_var * BC_Types[bc_input_id].BC_Data_Float[2]);
        }
      } else {
        d_func[index_var] = d_x_var * (BC_Types[bc_input_id].BC_Data_Float[1] +
                                       2. * x_var * BC_Types[bc_input_id].BC_Data_Float[2]);
      }
    }
    break;
  case (GD_PARAB_OFFSET_BC): /* C2 (x - C3)^2 + C1 (x - C3) + c0 */

    *func =
        ((x_var - BC_Types[bc_input_id].BC_Data_Float[3]) *
             (x_var - BC_Types[bc_input_id].BC_Data_Float[3]) *
             BC_Types[bc_input_id].BC_Data_Float[2] +
         (x_var - BC_Types[bc_input_id].BC_Data_Float[3]) * BC_Types[bc_input_id].BC_Data_Float[1] +
         BC_Types[bc_input_id].BC_Data_Float[0]);
    if (af->Assemble_Jacobian) {
      if (vector_sens) {
        for (b = 0; b < DIM; b++) {
          d_func[index_var + b] =
              d_vect_var[b] * (BC_Types[bc_input_id].BC_Data_Float[1] +
                               2. * (x_var - BC_Types[bc_input_id].BC_Data_Float[3]) *
                                   BC_Types[bc_input_id].BC_Data_Float[2]);
        }
      } else {
        d_func[index_var] = d_x_var * (BC_Types[bc_input_id].BC_Data_Float[1] +
                                       2. * (x_var - BC_Types[bc_input_id].BC_Data_Float[3]) *
                                           BC_Types[bc_input_id].BC_Data_Float[2]);
      }
    }
    break;
  case (GD_CIRC_BC): /* C2 ( x - C1 )^2  - c0^2 */
    /* C2 represents ellipticity and C1 represents origin */
    /* C0 is radius and should enter only one of the BC's */

    *func =
        (BC_Types[bc_input_id].BC_Data_Float[2] * (x_var - BC_Types[bc_input_id].BC_Data_Float[1]) *
             (x_var - BC_Types[bc_input_id].BC_Data_Float[1]) -
         BC_Types[bc_input_id].BC_Data_Float[0] * BC_Types[bc_input_id].BC_Data_Float[0]);

    if (af->Assemble_Jacobian) {
      if (vector_sens) {
        for (b = 0; b < DIM; b++) {
          d_func[index_var + b] =
              d_vect_var[b] * (2. * BC_Types[bc_input_id].BC_Data_Float[2] *
                               (x_var - BC_Types[bc_input_id].BC_Data_Float[1]));
        }
      } else {
        d_func[index_var] = d_x_var * (2. * BC_Types[bc_input_id].BC_Data_Float[2] *
                                       (x_var - BC_Types[bc_input_id].BC_Data_Float[1]));
      }
    }
    break;

  case (GD_POLYN_BC): /* up to 6th order polynomial */
                      /* C6 x^6 + C5 x^5 + C4 x^4 + C3 x^3 + C2 x^2 + C1 x + C0 */

    *func =
        (x_var * x_var * x_var * x_var * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[6] +
         x_var * x_var * x_var * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[5] +
         x_var * x_var * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[4] +
         x_var * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[3] +
         x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[2] +
         x_var * BC_Types[bc_input_id].BC_Data_Float[1] + BC_Types[bc_input_id].BC_Data_Float[0]);

    /* printf("POLYN fit X,F = %f %f\n", x_var, *func); */

    if (af->Assemble_Jacobian) {
      if (vector_sens) {
        for (b = 0; b < DIM; b++) {
          d_func[index_var + b] =
              d_vect_var[b] *
              (BC_Types[bc_input_id].BC_Data_Float[1] +
               2. * x_var * BC_Types[bc_input_id].BC_Data_Float[2] +
               3. * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[3] +
               4. * x_var * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[4] +
               5. * x_var * x_var * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[5] +
               6. * x_var * x_var * x_var * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[6]);
        }
      } else {
        d_func[index_var] =
            d_x_var *
            (BC_Types[bc_input_id].BC_Data_Float[1] +
             2. * x_var * BC_Types[bc_input_id].BC_Data_Float[2] +
             3. * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[3] +
             4. * x_var * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[4] +
             5. * x_var * x_var * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[5] +
             6. * x_var * x_var * x_var * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[6]);
      }
    }
    break;

  case (GD_TABLE_BC):

    *func = BC_Types[bc_input_id].BC_Data_Float[0];

    x_var_mp[0] = x_var;
    *func *= interpolate_table(BC_Types[bc_input_id].table, x_var_mp, &slope, NULL);

    if (af->Assemble_Jacobian) {
      if (vector_sens) {
        for (b = 0; b < DIM; b++) {
          d_func[index_var + b] = BC_Types[bc_input_id].BC_Data_Float[0] * slope * d_vect_var[b];
        }
      } else {
        d_func[index_var] = BC_Types[bc_input_id].BC_Data_Float[0] * slope * d_x_var;
      }
    }

    break;

  default:
    return (-1);
  }

  return (0);
} /* END of routine fgeneralized_dirichlet                                   */
/*****************************************************************************/

int evaluate_time_func(const double current_time,
                       double *f_time, /* computed time function */
                       const int bc_input_id)

/************************************************************************
 *
 * Function which multiplies a time function by previously
 * loaded GD conditions.
 ************************************************************************/
{
  double time = current_time;
  int time_function;
  /* Check if max time was specified and reset time if greater than max time */
  if (BC_Types[bc_input_id].BC_Data_Int[4] == GD_TIME_MAX) {
    if (time > BC_Types[bc_input_id].BC_Data_Float[2]) {
      time = BC_Types[bc_input_id].BC_Data_Float[2];
    }
  }

  /* Check if max time was specified  and reset time if so */
  if (BC_Types[bc_input_id].BC_Data_Int[3] == GD_TIME_MAX) {
    if (time > BC_Types[bc_input_id].BC_Data_Float[2]) {
      time = BC_Types[bc_input_id].BC_Data_Float[2];
    }
  }

  /* ---- Find variable number and species number */
  time_function = BC_Types[bc_input_id].BC_Data_Int[2];

  /* Now compute time function */

  switch (time_function) {
  case (GD_TIME_LIN): /* c0 +c1*t */

    *f_time =
        BC_Types[bc_input_id].BC_Data_Float[0] + BC_Types[bc_input_id].BC_Data_Float[1] * time;

    break;

  case (GD_TIME_EXP): /* exp(c0 +c1*t) */

    *f_time =
        exp(BC_Types[bc_input_id].BC_Data_Float[0] + BC_Types[bc_input_id].BC_Data_Float[1] * time);

    break;

  case (GD_TIME_SIN): /* sin(c0 +c1*t) */

    *f_time =
        sin(BC_Types[bc_input_id].BC_Data_Float[0] + BC_Types[bc_input_id].BC_Data_Float[1] * time);

    break;
  case (GD_TIME_TABLE): {
    double slope, _time[1];
    _time[0] = time;
    *f_time = interpolate_table(BC_Types[bc_input_id].table, _time, &slope, NULL);
  } break;
  default:
    return (-1);
  }

  return (0);
} /* END of routine evaluate_time_function                                   */