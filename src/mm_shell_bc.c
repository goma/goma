/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2014 Sandia Corporation.                                  *
*                                                                         *
* Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,  *
* the U.S. Government retains certain rights in this software.            *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
\************************************************************************/
 


/* Standard include files */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* GOMA include files */

#include "std.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_masks.h"
#include "rf_io_const.h"
#include "rf_io_structs.h"
#include "rf_io.h"
#include "rf_mp.h"
#include "el_elm.h"
#include "el_geom.h"
#include "rf_bc_const.h"
#include "rf_solver_const.h"
#include "rf_solver.h"
#include "rf_fill_const.h"
#include "rf_vars_const.h"
#include "mm_mp_const.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"

#include "mm_mp_structs.h"
#include "mm_mp.h"
#include "mm_std_models_shell.h"
#include "mm_qp_storage.h"

#include "mm_eh.h"

#define eps(i,j,k)   ( (i-j)*(j-k)*(k-i)/2 )

#define _MM_SHELL_BC_C
#include "goma.h"

/*
 * Global variables defined here. Declared frequently via rf_bc.h
 */



/*********** R O U T I N E S   I N   T H I S   F I L E ***********************
*
*						
*				NOTE:		 
*						 
*						 
*						 
*
*       NAME			        TYPE			CALLED BY
*  -----------------------------------------------------------------
*
*  shell_n_dot_flow_bc_confined         void                   
*
*  shell_n_dot_flow_bc_film             void                   
*
*  shell_t_dot_flow_bc_film             void
*
*  shell_n_dot_gradh_bc                 void
*
*  shell_n_dot_pflux_bc                 void
*
*
******************************************************************************/


void
shell_n_dot_flow_bc_confined(double func[DIM],
                             double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                             const double flowrate, /* imposed flow rate */
                             const double time,     /* current time */
                             const double dt,       /* current time step size */
                             double xi[DIM],        /* Local stu coordinates */
                             const Exo_DB *exo)

     /***********************************************************************
      *
      * shell_n_dot_flow_bc_confined():
      *
      *  Function which evaluates the expression specifying the
      *  pressure gradient (or flow rate) at a quadrature point normal to the side
      *  of an element.
      *
      *         func =   - flowrate + n .(  H^3 /(12*mu) * (-grad LUB_P)
      *                                   + 0.5 * (U_bot + U_top) * H )
      *
      *  The boundary condition GRAD_LUB_PRESS_BC employs this function.
      *
      *
      * Input:
      *
      *  flowrate      = specified on the bc card as the first float
      *  grad LUB_P    = Lubrication pressure gradient
      *  U_top         = Velocity of the top wall
      *  U_bot         = Velocity of the bottom wall
      *  H             = Distance between top and bottom wall
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
      *
      *   Author: K. Tjiptowidjojo    (12/13/2010)
      *
      ********************************************************************/
{
  int j, ii, var;
  int *n_dof = NULL;
  int dof_map[MDE];
  double grad_phi_j[DIM], grad_II_phi_j[DIM];
  double bound_normal[DIM];


/* Save the boundary normal vector */

  for(ii = 0; ii < pd->Num_Dim; ii++)
    { 
      bound_normal[ii] = fv->snormal[ii];
    }

 /*
  * Prepare geometry
  */
  n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
  lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);


/* Calculate the flow rate and its sensitivties */

  calculate_lub_q_v(R_LUBP, time, dt, xi, exo);



  if (af->Assemble_LSA_Mass_Matrix)
    {
      return;
    }

  if (af->Assemble_Jacobian)
    {
      var = LUBP;
      if (pd->v[var])
      {
        for ( j=0; j<ei->dof[var]; j++)
           {
            for (ii = 0; ii < pd->Num_Dim; ii++)
              {
                 grad_phi_j[ii] = bf[var]->grad_phi[j][ii];
                 grad_II_phi_j[ii] = 0.0;
              }

            Inn(grad_phi_j, grad_II_phi_j);

            for (ii=0; ii<pd->Num_Dim; ii++)
              {
               d_func[0][var][j] += LubAux->dq_dp1[ii][j] * grad_II_phi_j[ii] * bound_normal[ii];
              }
           }
      }


    } /* end of if Assemble_Jacobian */


  /* Calculate the residual contribution        */

  func[0] = - flowrate;
  for (ii = 0; ii < pd->Num_Dim; ii++)
    {
      func[0] +=  LubAux->q[ii] * bound_normal[ii];
    }

  /* clean-up */
  safe_free((void *) n_dof);

} /* END of routine shell_n_dot_flow_bc_confined  */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void
lub_static_pressure(double func[DIM],
                    double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                    const double P_atm,    /* imposed atmospheric pressure */
                    const double time,     /* current time */
                    const double dt,	   /* current time step size */
                    double xi[DIM],        /* Local stu coordinates */
                    const Exo_DB *exo)

     /***********************************************************************
      *
      * lub_static_pressure():
      *
      *  Function which evaluates the expression specifying the
      *  pressure at a quadrature point normal to the side
      *  of an element to be in static equilibrium.
      *
      *         func =   LUB_P - (P_atm - HEAVISIDE * CURV * sigma/H)
      *
      *  The boundary condition LUB_STATIC employs this function.
      *
      *
      * Input:
      *
      *
      *  LUB_P         = Lubrication pressure
      *  P_atm         = specified on the bc card as the first float
      *  HEAVISIDE     = Heaviside function -- active only when F < 0
      *  CURV          = Analytical curvature in direction normal to the substrate
      *  sigma         = Surface tension
      *  H             = Distance between top and bottom wall
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
      *
      *   Author: K. Tjiptowidjojo    (08/05/2014)
      *
      ********************************************************************/
{
  int i, j, k, jk, var;
  int dim;
  int *n_dof = NULL;
  int dof_map[MDE];

  dbl phi_j;

 /*
  * Prepare geometry
  */

  dim = pd->Num_Dim;

  n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
  lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);


  /* Extract wall heights */
  dbl H, H_U, dH_U_dtime, H_L, dH_L_dtime;
  dbl dH_U_dX[DIM],dH_L_dX[DIM], dH_U_dp, dH_U_ddh;
  H = height_function_model(&H_U, &dH_U_dtime, &H_L, &dH_L_dtime, dH_U_dX, dH_L_dX, &dH_U_dp, &dH_U_ddh, time, dt);

  /***** DEFORM HEIGHT AND CALCULATE SENSITIVITIES *****/

  /* Define variables */
  dbl D_H_DX[DIM][MDE], D_H_DP[MDE];
  dbl D_H_DRS[DIM][MDE];
  dbl D_H_DNORMAL[DIM][MDE];
  memset(D_H_DX,  0.0, sizeof(double)*DIM*MDE);
  memset(D_H_DRS,  0.0, sizeof(double)*DIM*MDE);
  memset(D_H_DNORMAL,  0.0, sizeof(double)*DIM*MDE);
  memset(D_H_DP,  0.0, sizeof(double)*MDE);

  /* Deform height */
  switch ( mp->FSIModel ) {
  case FSI_MESH_CONTINUUM:
  case FSI_MESH_UNDEF:
  case FSI_SHELL_ONLY_UNDEF:
  for ( i = 0; i < dim; i++)
     {
      H -= fv->snormal[i] * fv->d[i];
     }
  break;
  case FSI_SHELL_ONLY_MESH:
    if (pd->e[R_SHELL_NORMAL1] && pd->e[R_SHELL_NORMAL2] && pd->e[R_SHELL_NORMAL3] )
      {
       for ( i = 0; i < dim; i++)
          {
           H -= fv->n[i] * fv->d[i];
          }
      }
    else
      {
       for ( i = 0; i < dim; i++)
          {
           H -= fv->snormal[i] * fv->d[i];
          }
      }
  break;
  case FSI_REALSOLID_CONTINUUM:
  for ( i = 0; i < dim; i++) 
     {
      H -= fv->snormal[i] * fv->d_rs[i];
     }
  break;
 }

 /* Calculate height sensitivity to mesh */
 switch ( mp->FSIModel ) {
 case FSI_MESH_CONTINUUM:
 case FSI_MESH_UNDEF:
 case FSI_SHELL_ONLY_UNDEF:
 for ( i = 0; i < dim; i++)
    {
     for ( j = 0; j < dim; j++)
        {
         for ( k = 0; k < ei->dof[MESH_DISPLACEMENT1]; k++)
            {
              jk = dof_map[k];
              D_H_DX[j][jk] += delta(i,j)*(dH_U_dX[i]-dH_L_dX[i])*bf[MESH_DISPLACEMENT1]->phi[k];
              D_H_DX[j][jk] -= fv->dsnormal_dx[i][j][jk] * fv->d[i];
              D_H_DX[j][jk] -= fv->snormal[i] * delta(i,j) * bf[MESH_DISPLACEMENT1]->phi[k];
            }
	}
    }
  break;
 case FSI_SHELL_ONLY_MESH:
   if ( (pd->e[R_SHELL_NORMAL1]) && (pd->e[R_SHELL_NORMAL2]) && (pd->e[R_SHELL_NORMAL3]) )
     {
      for ( i = 0; i < dim; i++)
         {
          for ( j = 0; j < dim; j++)
             {
              for ( k = 0; k < ei->dof[MESH_DISPLACEMENT1]; k++)
                 {
                  jk = dof_map[k];
                  D_H_DX[j][jk] += delta(i,j)*(dH_U_dX[i]-dH_L_dX[i])*bf[MESH_DISPLACEMENT1]->phi[k];
                  D_H_DX[j][jk] -= fv->n[i] * delta(i,j) * bf[MESH_DISPLACEMENT1]->phi[k];
                 }
             }
         }
     }
   else
     {
      for ( i = 0; i < dim; i++)
         {
          for ( j = 0; j < dim; j++)
             {
              for ( k = 0; k < ei->dof[MESH_DISPLACEMENT1]; k++)
                 {
                  jk = dof_map[k];
                  D_H_DX[j][jk] += delta(i,j)*(dH_U_dX[i]-dH_L_dX[i])*bf[MESH_DISPLACEMENT1]->phi[k];
                  D_H_DX[j][jk] -= fv->dsnormal_dx[i][j][jk] * fv->d[i];
                  D_H_DX[j][jk] -= fv->snormal[i] * delta(i,j) * bf[MESH_DISPLACEMENT1]->phi[k];
                 }
             }
         }
     }
  break;
 case FSI_REALSOLID_CONTINUUM:
 for ( i = 0; i < dim; i++) 
    {
     for ( j = 0; j < dim; j++) 
        {
         for ( k = 0; k < ei->dof[MESH_DISPLACEMENT1]; k++) 
            {
             jk = dof_map[k];
             D_H_DX[j][jk] += delta(i,j)*(dH_U_dX[i]-dH_L_dX[i])*bf[MESH_DISPLACEMENT1]->phi[k];
             D_H_DX[j][jk] -= fv->dsnormal_dx[i][j][jk] * fv->d_rs[i];
	    }
         for ( k = 0; k < ei->dof[SOLID_DISPLACEMENT1]; k++) 
            {
             jk = dof_map[k];
             D_H_DRS[j][jk] -= fv->snormal[i] * delta(i,j) * bf[SOLID_DISPLACEMENT1]->phi[jk];
            }
        }
    }
  break;
 }

 /* Calculate height sensitivity to shell normal */
 switch ( mp->FSIModel ) {

 case FSI_SHELL_ONLY_MESH:
 if ( (pd->e[R_SHELL_NORMAL1]) && (pd->e[R_SHELL_NORMAL2]) && (pd->e[R_SHELL_NORMAL3]) )
   {
    for ( i = 0; i < dim; i++)
       {
        for ( j = 0; j < dim; j++)
           {
            for ( k = 0; k < ei->dof[SHELL_NORMAL1]; k++)
               {
                D_H_DNORMAL[j][k] -= delta(i,j) * bf[SHELL_NORMAL1]->phi[k] * fv->d[i];
               }
           }
       }
   }
  break;
 }

  /***** CALCULATE HEAVISIDE AND SENSITIVITIES *****/

  dbl Hn = 0.0;
  dbl D_Hn_DF[MDE], D_Hn_DX[DIM][MDE];
  memset(D_Hn_DF, 0.0, sizeof(double)*MDE);
  memset(D_Hn_DX, 0.0, sizeof(double)*DIM*MDE);


  if ( pd->v[FILL] )
    {
     load_lsi( ls->Length_Scale );
     load_lsi_derivs();
     Hn = 1.0 - lsi->Hn;
    }

  /* Calculate F sensitivity */
  if ( pd->v[FILL] )
    {
     for ( i = 0; i < ei->dof[FILL]; i++)
        {
         D_Hn_DF[i] = - lsi->d_Hn_dF[i];
        }
    }

  /***** CALCULATE CURVATURE AND SENSITIVITIES *****/

  dbl CURV = 0.0;
  dbl D_CURV_DH = 0.0;
  dbl D_CURV_DF[MDE], D_CURV_DX[DIM][MDE];
  memset(D_CURV_DF, 0.0, sizeof(double)*MDE);
  memset(D_CURV_DX, 0.0, sizeof(double)*DIM*MDE);

  /* Curvature - analytic in the "z" direction  */
  dbl dcaU, dcaL, slopeU, slopeL;
  dcaU = dcaL = slopeU = slopeL = 0;
  if ( pd->v[FILL] )
    {
     dcaU = mp->dcaU*M_PIE/180.0;
     dcaL = mp->dcaL*M_PIE/180.0;
     slopeU = slopeL = 0.;
     for ( i = 0; i < dim; i++)
        {
          slopeU += dH_U_dX[i]*lsi->normal[i];
          slopeL += dH_L_dX[i]*lsi->normal[i];
        }
     CURV += (cos(M_PIE-dcaU-atan(slopeU)) + cos(M_PIE-dcaL-atan(-slopeL)))/H ;
    }

  /* Sensitivity to height */
  if ( pd->v[FILL] )  D_CURV_DH = -(cos(M_PIE-dcaU-atan(slopeU)) + cos(M_PIE-dcaL-atan(-slopeL)))/(H*H);


  /* Sensitivity to level set F */
  if ( pd->v[FILL] )
    {
     for ( i = 0; i < ei->dof[FILL]; i++)
        {
         for ( j = 0; j < dim; j++)
            {
             D_CURV_DF[i] += sin(dcaU+atan(slopeU))/(H*(1+slopeU*slopeU))*dH_U_dX[j]**lsi->d_normal_dF[j];
             D_CURV_DF[i] += sin(dcaL+atan(slopeL))/(H*(1+slopeL*slopeL))*dH_L_dX[j]**lsi->d_normal_dF[j];
            }
	}
    }

  /* Sensitivity to mesh */
  if ( pd->v[FILL] )
    {
     for ( i = 0; i < dim; i++)
        {
         for ( j = 0; j < n_dof[MESH_DISPLACEMENT1]; j++)
            {
             D_CURV_DX[i][j] += D_CURV_DH * D_H_DX[i][j];
            }
	}
    }


  /***** CALCULATE RESIDUAL CONTRIBUTION *****/

  func[0] = fv->lubp - (P_atm + mp->surface_tension * Hn * CURV);


  /***** CALCULATE JACOBIAN CONTRIBUTION *****/

  if (af->Assemble_Jacobian)
    {

      /* Sensitivity w.r.t. pressure */
      var = LUBP;
      if (pd->v[var])
        {
         for (j = 0; j < ei->dof[var]; j++)
            {
             phi_j = bf[var]->phi[j];

             d_func[0][var][j] = phi_j - mp->surface_tension * Hn * D_CURV_DH * D_H_DP[j];
            }
	}


      /* Sensitivity w.r.t. level set */
      var = FILL;
      if (pd->v[var])
        {
         for (j = 0; j < ei->dof[var]; j++)
            {
             d_func[0][var][j] = - mp->surface_tension * (D_Hn_DF[j] * CURV + Hn * D_CURV_DF[j]);
            }
	}

      /* Sensitivity w.r.t. mesh */
      if (pd->v[MESH_DISPLACEMENT1])
        {
         for (i = 0; i < dim; i++)
            {
             var = MESH_DISPLACEMENT1 + i;

             for (j = 0; j < ei->dof[var]; j++)
                {
                 d_func[0][var][j] = - mp->surface_tension * (D_Hn_DX[i][j] * CURV + Hn * D_CURV_DX[i][j]);
                }
            }
	}

    } /* end of if Assemble_Jacobian */



  /* clean-up */
  safe_free((void *) n_dof);

} /* END of routine lub_static_pressure  */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void 
shell_n_dot_flow_bc_film(double func[DIM],
		         double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
		         const double flowrate, /* imposed flow rate */
                         const double time,     /* current time */
                         const double dt,       /* current time step size */
                         double xi[DIM],        /* Local stu coordinates */
                         const Exo_DB *exo)

     /***********************************************************************
      *
      * shell_n_dot_flow_bc_film():
      *
      *  Function which evaluates the expression specifying the
      *  pressure gradient (or flow rate) at a quadrature point normal to the side
      *  of an element.
      *
      *         func =   - flowrate + n .(  SH_FH^3 /(3*mu) * (-grad SH_FP)
      *                                   + U_bot * SH_FH )
      *
      *  The boundary condition SHELL_GRAD_FP_BC employs this function.
      *
      *
      * Input:
      *
      *  flowrate      = specified on the bc card as the first float
      *  grad SH_FP    = Lubrication pressure gradient
      *  U_bot         = Velocity of the substrate
      *  SH_FH         = Film thickness
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
      *   Author: K. Tjiptowidjojo    (3/1/2010)
      * 
      ********************************************************************/
{
  int j, ii, var;
  int *n_dof = NULL;
  int dof_map[MDE];
  double phi_j;
  double grad_phi_j[DIM], grad_II_phi_j[DIM];
  double bound_normal[DIM];



/* Save the boundary normal vector */

  for(ii = 0; ii < pd->Num_Dim; ii++)
    {
      bound_normal[ii] = fv->snormal[ii];
    }

/*
  * Prepare geometry
  */
  n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
  lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);


/* Calculate the flow rate and its sensitivties */

  calculate_lub_q_v(R_SHELL_FILMP, time, dt, xi, exo);


  if (af->Assemble_LSA_Mass_Matrix)
    {
      return;
    }
  

  if (af->Assemble_Jacobian) 
    {
      var = SHELL_FILMP;
      if (pd->v[var])
      {
        for ( j=0; j<ei->dof[var]; j++)
           {
            for (ii = 0; ii < pd->Num_Dim; ii++)
              {
                 grad_phi_j[ii] = bf[var]->grad_phi[j][ii];
                 grad_II_phi_j[ii] = 0.0;
              }

            Inn(grad_phi_j, grad_II_phi_j);
                       	              
            for (ii=0; ii<pd->Num_Dim; ii++)
              {
	       d_func[0][var][j] += LubAux->dq_dp1[ii][j] * grad_II_phi_j[ii] * bound_normal[ii];
              }
           }
      }

      var = SHELL_FILMH;
      if (pd->v[var])
      {
        for ( j=0; j<ei->dof[var]; j++)
           {
            phi_j = bf[var]->phi[j];

            for (ii = 0; ii < pd->Num_Dim; ii++)
              {
                 grad_phi_j[ii] = bf[var]->grad_phi[j][ii];
                 grad_II_phi_j[ii] = 0.0;
              }

            Inn(grad_phi_j, grad_II_phi_j);

                       	              
            for (ii=0; ii<pd->Num_Dim; ii++)
              {
	       d_func[0][var][j] += LubAux->dq_dh1[ii][j] * grad_II_phi_j[ii] * bound_normal[ii];
	       d_func[0][var][j] += LubAux->dq_dh2[ii][j] * phi_j * bound_normal[ii];
              }
           }
      }

      var = SHELL_PARTC;
      if (pd->v[var])
      {
        for ( j=0; j<ei->dof[var]; j++)
           {
            phi_j = bf[var]->phi[j];
                       	              
            for (ii=0; ii<pd->Num_Dim; ii++)
              {
	       d_func[0][var][j] += LubAux->dq_dc[ii][j] * phi_j * bound_normal[ii];
              }
           }
      }


    } /* end of if Assemble_Jacobian */
  

  /* Calculate the residual contribution	*/

  func[0] = - flowrate;
  for (ii = 0; ii < pd->Num_Dim; ii++)
    {
      func[0] += LubAux->q[ii] * bound_normal[ii]; 
    }
   
  /* clean-up */
  safe_free((void *) n_dof);

} /* END of routine shell_n_dot_flow_bc_film  */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/


void 
shell_n_dot_gradp_bc(double func[DIM],
		     double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                     const double time,     /* current time */
                     const double dt,       /* current time step size */
                     double xi[DIM],        /* Local stu coordinates */
                     const Exo_DB *exo)

     /***********************************************************************
      *
      * shell_n_dot_gradp_bc():
      *
      *  Function which sets the pressure gradient at the side to be zero
      *
      *         func =   n .(-grad SH_FP + grad DisjPress)
      *                         
      *
      *  The boundary condition SHELL_FLOW_DEVELOPED_BC employs this function.
      *
      *
      * Input:
      *
      *  grad SH_FP    = Lubrication pressure gradient
      *  grad_DisjPress = Disjoining pressure gradient
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
      *   Author: K. Tjiptowidjojo    (4/27/2011)
      * 
      ********************************************************************/
{
  int j, ii, var;
  int *n_dof = NULL;
  int dof_map[MDE];
  double grad_P[DIM];
  double grad_DisjPress[DIM], dgrad_DisjPress_dH1[DIM][MDE], dgrad_DisjPress_dH2[DIM][MDE]; 
  double phi_j;
  double grad_phi_j[DIM], grad_II_phi_j[DIM];
  double bound_normal[DIM];



/* Save the boundary normal vector */

  for(ii = 0; ii < pd->Num_Dim; ii++)
    {
      bound_normal[ii] = fv->snormal[ii];
    }

/*
  * Prepare geometry
  */
  n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
  lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);


/* Calculate pressure gradient  */

   Inn( fv->grad_sh_fp, grad_P );

/* Calculate disjoining pressure gradient and its sensitivities */
  
   disjoining_pressure_model(fv->sh_fh, fv->grad_sh_fh, grad_DisjPress, dgrad_DisjPress_dH1, dgrad_DisjPress_dH2);

  if (af->Assemble_LSA_Mass_Matrix)
    {
      return;
    }
  

  if (af->Assemble_Jacobian) 
    {
      var = SHELL_FILMP;
      if (pd->v[var])
      {
        for ( j=0; j<ei->dof[var]; j++)
           {
            for (ii = 0; ii < pd->Num_Dim; ii++)
              {
                 grad_phi_j[ii] = bf[var]->grad_phi[j][ii];
                 grad_II_phi_j[ii] = 0.0;
              }

            Inn(grad_phi_j, grad_II_phi_j);
                       	              
            for (ii=0; ii<pd->Num_Dim; ii++)
              {
	       d_func[0][var][j] += - grad_II_phi_j[ii] * bound_normal[ii];
              }
           }
      }

      var = SHELL_FILMH;
      if (pd->v[var])
      {
        for ( j=0; j<ei->dof[var]; j++)
           {
            phi_j = bf[var]->phi[j];

            for (ii = 0; ii < pd->Num_Dim; ii++)
              {
                 grad_phi_j[ii] = bf[var]->grad_phi[j][ii];
                 grad_II_phi_j[ii] = 0.0;
              }

            Inn(grad_phi_j, grad_II_phi_j);

                       	              
            for (ii=0; ii<pd->Num_Dim; ii++)
              {
	       d_func[0][var][j] += dgrad_DisjPress_dH1[ii][j] * grad_II_phi_j[ii] * bound_normal[ii];
	       d_func[0][var][j] += dgrad_DisjPress_dH2[ii][j] * phi_j * bound_normal[ii];
              }
           }
      }

    } /* end of if Assemble_Jacobian */
  

  /* Calculate the residual contribution	*/

  for (ii = 0; ii < pd->Num_Dim; ii++)
    {
      func[0] += (- grad_P[ii] + grad_DisjPress[ii]) * bound_normal[ii]; 
    }
   
  /* clean-up */
  safe_free((void *) n_dof);

} /* END of routine shell_n_dot_gradp_bc  */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/



void 
shell_n_dot_gradh_bc(double func[DIM],
		     double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
		     const double slope) /* imposed slope */

     /***********************************************************************
      *
      * shell_n_dot_gradh_bc():
      *
      *  Function which evaluates the expression specifying the
      *  film slope  at a quadrature point normal to the side
      *  of an element.
      *
      *         func =   - slope + n . grad SH_FH
      *
      *  The boundary condition SHELL_GRAD_FH_BC employs this function.
      *
      *
      * Input:
      *
      *  slope         = specified on the bc card as the first float
      *  grad SH_FH    = Film slope at the Gauss point
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
      *   Author: K. Tjiptowidjojo    (3/2/2010)
      * 
      ********************************************************************/
{
  int j, ii, jj, var;
  double grad_phi_j[DIM], grad_II_phi_j[DIM];
  double grad_II_H[DIM];
  double bound_normal[DIM], shell_normal[DIM];

/* Save the boundary normal vector */

  for(ii = 0; ii < VIM; ii++)
    {
      grad_II_H[ii] = 0.;
      bound_normal[ii] = fv->snormal[ii];
    }


/* Get the vector normal to the shell plane */

   shell_determinant_and_normal(ei->ielem, ei->iconnect_ptr, ei->num_local_nodes,
                                ei->ielem_dim, 1);

  for(ii = 0; ii < VIM; ii++)
    {
      shell_normal[ii] = fv->snormal[ii];
    }

/* grad H = (I - nn) dot grad sh fh */

  for (ii = 0; ii < VIM; ii++)
    {
      for (jj = 0; jj < VIM; jj++)
        {
          grad_II_H[ii] += (  delta(ii,jj) * fv->grad_sh_fh[jj]
                         - shell_normal[ii] * shell_normal[jj] * fv->grad_sh_fh[jj] );
        }
    }


  if (af->Assemble_LSA_Mass_Matrix)
    {
      return;
    }
  

  if (af->Assemble_Jacobian) 
    {
      var = SHELL_FILMH;
      if (pd->v[var])
      {
        for ( j=0; j<ei->dof[var]; j++)
           {	              
            for (ii=0; ii<VIM; ii++)
              {
               grad_phi_j[ii] = bf[var]->grad_phi[j][ii];
	       grad_II_phi_j[ii] = 0.0;
              }
            for (ii = 0; ii<VIM; ii++)
              {
               for (jj = 0; jj<VIM; jj++)
                 {
                   grad_II_phi_j[ii] += (delta(ii,jj)*grad_phi_j[jj]
                                        - shell_normal[ii]*shell_normal[jj]*grad_phi_j[jj]);
                 }
              }
            for (ii=0; ii<VIM; ii++)
              {
	       d_func[0][var][j] += grad_II_phi_j[ii] * bound_normal[ii];
              }

           }
      }
    } /* end of if Assemble_Jacobian */
  
  /* Calculate the residual contribution	*/
  func[0] = - slope;
  for (ii = 0; ii < VIM; ii++)
    {
      func[0] += bound_normal[ii] * grad_II_H[ii];
    }

} /* END of routine shell_n_dot_gradh_bc  */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/


void 
shell_n_dot_pflux_bc(double func[DIM],
		     double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
		     const double flux,     /* imposed particles flux */
                     const double time,     /* current time */
                     const double dt)       /* current time step size */
                     


     /***********************************************************************
      *
      * shell_n_dot_pflux_bc():
      *
      *  Function which evaluates the expression specifying the
      *  particles flux  at a quadrature point normal to the side
      *  of an element.
      *
      *         func =   - flux + n . (   SH_PC * SH_FH^3/(3*mu) * (- grad_SH_FP)    
      *                                 + SH_PC * Ubot * SH_FH 
      *                                 - SH_FH * Diff * grad SH_PC )
      *
      *  The boundary condition SHELL_GRAD_PC_BC employs this function.
      *
      *
      * Input:
      *
      *  flux         = specified on the bc card as the first float
      *  SH_PC        = particles concentration
      *  SH_FH        = film thickness
      *  grad_SH_FP   = film lubrication pressure gradient
      *  grad SH_FH   = film slope
      *  grad_SH_PC   = particles concentration gradient
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
      *   Author: K. Tjiptowidjojo    (3/10/2010)
      * 
      ********************************************************************/
{
  int j, ii, jj, var;
  double phi_j;
  double grad_phi_j[DIM], grad_phi_j_corrected[DIM];
  double veloU[DIM], veloL[DIM];
  double grad_C[DIM];
  double bound_normal[DIM], shell_normal[DIM];
  double H;
  double mu, dmu_dc, diff_coeff, ddiff_dmu, ddiff_dc;
  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct;  /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;

  H = fv->sh_fh;
  mu = viscosity(gn, NULL, d_mu);
  diff_coeff = diffusion_coefficient_model(mu, &ddiff_dmu);
  dmu_dc = mp->d_viscosity[SHELL_PARTC];
  ddiff_dc = ddiff_dmu * dmu_dc;
  velocity_function_model(veloU, veloL, time, dt);


/* Save the boundary normal vector */

  for(ii = 0; ii < VIM; ii++)
    {
      grad_C[ii] = 0.0;
      bound_normal[ii] = fv->snormal[ii];
    }

  shell_determinant_and_normal(ei->ielem, ei->iconnect_ptr, ei->num_local_nodes,
                               ei->ielem_dim, 1);

  for(ii = 0; ii < VIM; ii++)
    {
      shell_normal[ii] = fv->snormal[ii];
    }

/* Peform I- nn to gradient operator */

  for (ii = 0; ii < VIM; ii++)
    {
      for (jj = 0; jj < VIM; jj++)
        {
          grad_C[ii] += (  delta(ii,jj) * fv->grad_sh_pc[jj] 
                         - shell_normal[ii] * shell_normal[jj] * fv->grad_sh_pc[jj] );
        }
    }


  if (af->Assemble_LSA_Mass_Matrix)
    {
      return;
    }


  if (af->Assemble_Jacobian) 
    {

      var = SHELL_FILMH;
      if (pd->v[var])
      {
        for ( j=0; j<ei->dof[var]; j++)
           {
            phi_j = bf[var]->phi[j];
     
            if ( fabs(flux) < 1.0e-10 )
              {
               d_func[0][var][j] += 0.0;
              }
            else
              {
               for (ii=0; ii<VIM; ii++)
                  {
	           d_func[0][var][j] += phi_j * diff_coeff * grad_C[ii] * bound_normal[ii];
                  }
              }
           }
      }


      var = SHELL_PARTC;
      if (pd->v[var])
      {
        for ( j=0; j<ei->dof[var]; j++)
           {
            phi_j = bf[var]->phi[j];

            for (ii = 0; ii < VIM; ii++)
              {
                 grad_phi_j[ii] = bf[var]->grad_phi[j][ii];
                 grad_phi_j_corrected[ii] = 0.0;
              }
            for (ii = 0; ii < VIM; ii++)
              {
                 for (jj = 0; jj < VIM; jj++)
                    {
                     grad_phi_j_corrected[ii] += (delta(ii,jj)*grad_phi_j[jj]
                                                 - shell_normal[ii]*shell_normal[jj]*grad_phi_j[jj]);
                    }
              }
            if ( fabs(flux) == 0.0 )
              {
               for (ii=0; ii<VIM; ii++)
                  {
	           d_func[0][var][j] += grad_phi_j_corrected[ii] * bound_normal[ii];
                  }
              }
            else
              {
               for (ii=0; ii<VIM; ii++)
                  {
	           d_func[0][var][j] += H * ddiff_dc * grad_C[ii] * bound_normal[ii];
	           d_func[0][var][j] += H * diff_coeff * grad_phi_j_corrected[ii] * bound_normal[ii];
                  }
              }
           }
      }


    } /* end of if Assemble_Jacobian */

  /* Calculate the residual contribution	*/


  func[0] = - flux;
  if ( fabs(flux) == 0. )
    {
     for (ii = 0; ii < VIM; ii++)
        {
         func[0] += bound_normal[ii] * grad_C[ii];
        }
    }
  else
    {
     for (ii = 0; ii < VIM; ii++)
        {
         func[0] += bound_normal[ii] * H * diff_coeff * grad_C[ii];
        }
    }

} /* END of routine shell_n_dot_pflux_bc  */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/****************************************************************************/

void
apply_shell_traction_bc(double func[DIM],
                       	double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                        const int BC_name, /* BC name identifier */
                        const dbl tx,   /* Traction in x-direction */
                        const dbl ty,   /* Traction in y-direction */
                        const dbl tz)	   /* Traction in z-direction */

     /***********************************************************************
      *
      * apply_shell_traction_bc():
      *
      *  Function which which evaluates traction  at a quadrature point
      *  normal to the side of an element.
      *
      *   func[0] =   e . (i tx + j ty + k tz)
      *
      *  e is either e1 or e2 depending on the type of BC
      *
      *  The boundary condition SH_S11_WEAK_BC and SH_S22_WEAK_BC
      *  employ this function.
      *
      *
      * Input:
      *
      * tx, ty, tz   = Traction components
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
      *   Author: K. Tjiptowidjojo    (5/22/2014)
      *
      ********************************************************************/
{
  int var, dim, a, b;
  int j;

  dbl t0[DIM];
  dbl t1[DIM];
  dbl dt0_dx[DIM][DIM][MDE];
  dbl dt1_dx[DIM][DIM][MDE];
  dbl dt0_dnormal[DIM][DIM][MDE];
  dbl dt1_dnormal[DIM][DIM][MDE];

  dbl e[DIM];
  dbl de_dx[DIM][DIM][MDE];
  dbl de_dnormal[DIM][DIM][MDE];

/************** PRECALCULATION ***********************/

  /* Unpack variables from structures for local convenience... */
  dim = pd->Num_Dim;


  /******* TANGENTS AND THEIR SENSITIVITIES **********/

  memset(dt0_dx, 0.0, sizeof(double)*DIM*DIM*MDE);
  memset(dt1_dx, 0.0, sizeof(double)*DIM*DIM*MDE);
  memset(dt0_dnormal, 0.0, sizeof(double)*DIM*DIM*MDE);
  memset(dt1_dnormal, 0.0, sizeof(double)*DIM*DIM*MDE);

  shell_tangents(t0, t1, dt0_dx, dt1_dx);

//  shell_tangents_seeded(t0, t1, dt0_dnormal, dt1_dnormal);


  for (a = 0; a < dim; a++)
     {
      if (BC_name == SH_S11_WEAK_BC)
        {
         e[a] = t0[a];
        }
      else if (BC_name == SH_S22_WEAK_BC)
        {
         e[a] = t1[a];
        }
     }

  memset( de_dx, 0.0, sizeof(double)*DIM*DIM*MDE);
  memset( de_dnormal, 0.0, sizeof(double)*DIM*DIM*MDE);

  for (b = 0; b < dim; b++)
     {
      var = MESH_DISPLACEMENT1 + b;

      for (j = 0; j < ei->dof[var]; j++)
         {
          for (a = 0; a < dim; a++)
             {
              if (BC_name == SH_S11_WEAK_BC)
                {
                 de_dx[a][b][j] = dt0_dx[a][b][j];
                }
              else if (BC_name == SH_S22_WEAK_BC)
                {
                 de_dx[a][b][j] = dt1_dx[a][b][j];
                }
             }
         }
     }

  for (b = 0; b < dim; b++)
     {
      var = SHELL_NORMAL1 + b;

      for (j = 0; j < ei->dof[var]; j++)
         {
          for (a = 0; a < dim; a++)
             {
              if (BC_name == SH_S11_WEAK_BC)
                {
                 de_dnormal[a][b][j] = dt0_dnormal[a][b][j];
                }
              else if (BC_name == SH_S22_WEAK_BC)
                {
                 de_dnormal[a][b][j] = dt1_dnormal[a][b][j];
                }
             }
         }
     }

/************** APPLY TRACTION ***********************/

  func[0] = e[0] * tx + e[1] * ty + e[2] * tz;

  if (af->Assemble_Jacobian)
    {

     var = MESH_DISPLACEMENT1;
     if ( pd->v[var] )
       {

        /*** Loop over dimensions of mesh displacement ***/
        for ( b = 0; b < dim; b++)
           {
            var = MESH_DISPLACEMENT1 + b;

            for (j = 0; j < ei->dof[var]; j++)
               {
                d_func[0][var][j] = de_dx[0][b][j] * tx + de_dx[1][b][j] * ty + de_dx[2][b][j] * tz;
               }
           }
       }

     var = SHELL_NORMAL1;
     if ( pd->v[var] )
       {

        /*** Loop over dimensions of shell normal ***/
        for ( b = 0; b < dim; b++)
           {
            var = SHELL_NORMAL1 + b;

            for (j = 0; j < ei->dof[var]; j++)
               {
                d_func[0][var][j] = de_dnormal[0][b][j] * tx + de_dnormal[1][b][j] * ty + de_dnormal[2][b][j] * tz;
               }
           }
       }
    }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/****************************************************************************/



void 
 match_lubrication_film_pressure(double func[],
				 double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
				 int eb_mat_lubp,
				 int eb_mat_filmp)
/******************************************************************************
*
*  Function which matches the lubrication pressure from a confined region to
*  a film pressure from a lubrication film region. 
*			 Author: A. Venkatakrishnan, Randy Schunk (11/12/2010)
******************************************************************************/
{
  int j_id;
  int var = -1;
  double phi_j;

  /* local contributions of boundary condition to residual and jacobian */

/***************************** EXECUTION BEGINS *******************************/


/***************************** Confined Lubrication SIDE *******************************/

  /*
   *  If current material is the confined lubrication region (lubp the variable)
   *  simply grab the lubrication pressue from the element and add it to the residual. 
   */
  if (Current_EB_ptr->Elem_Blk_Id == eb_mat_lubp)
    {    
      func[0] += fv->lubp;
	  
      if (af->Assemble_Jacobian) {

	var = LUBP;
	for ( j_id=0; j_id<ei->dof[var]; j_id++)
	  {
	    if (pd->v[var])
	      {
		phi_j = bf[var]->phi[j_id];
		d_func[0][var][j_id] += phi_j;
	      }
	  }
      }

    }

/***************************** FILM SIDE *******************************/
  /*
   *  If current material is the open film phase, grab the film pressure and 
   * add it to the residual. 
   */
  else if (Current_EB_ptr->Elem_Blk_Id == eb_mat_filmp)
    {
      func[0] -= fv->sh_fp;
	  
      if (af->Assemble_Jacobian) {

	var = SHELL_FILMP;
	for ( j_id=0; j_id<ei->dof[var]; j_id++)
	  {
	    if (pd->v[var])
	      {
		phi_j = bf[var]->phi[j_id];
		d_func[0][var][j_id] -= phi_j;
	      }
	  }
      }

  
    }
  else EH(-1,"LUBP_SH_FP_MATCH BC called with incorrect block id");
  /* note, we may not want to quit at this point */
  return;
}

/****************************************************************************/
/*
 * This function equates the mass flux on the boundary between a lubrication 
 * film region and a lubrication confined flow region (like the exit of a slot
 * coating flow) by adding the local contribution of the
 * volume integrals in the confined flow region to the total contribution in
 * the film flow region, and decrementing the volume contribution in that phase
 *
 */

void 
put_lub_flux_in_film(int id, /* local element node number for the 
			      * current node whose residual contribution
			      * is being sought                        */
		     int I, /* Global node number                      */
		     int ielem_dim, /* physical dimension of the elem  */
		     double resid_vector[], /* Residual vector         */
		     int i_mat_lubp,  /* elem block id's of confined film */
		     int i_mat_filmp, /* elem block id's of free film     */
		     int local_node_list_fs[]) /* MDE list to keep track
						* of nodes at which 
						* solid contributions 
						* have been transfered
						* to liquid (fluid-solid
						* boundaries)          */

{
    int j_id, dim, var, pvar, q, id_doflubp, id_doffilmp;
    int peqn_lubp, peqn_filmp;
    int ieqn_lubp, ieqn_filmp;

    dim = pd->Num_Dim;

    /*
     * if you are in the film phase, return without doing anything
     * In the solid phase, there are no fluid momentum equations.
     */
    if (!pd->e[R_LUBP]) return;   

    
    id_doflubp = ei->ln_to_dof[R_LUBP][id];
    id_doffilmp = ei->ln_to_dof[R_SHELL_FILMP][id];

    if (Current_EB_ptr->Elem_Blk_Id != i_mat_lubp)
      { 
	if (Current_EB_ptr->Elem_Blk_Id != i_mat_filmp) {
	  EH(-1, "put_lub_flux_in_film: Improper lubp and filmp block ids");
	}
	/* note, this may not account for all situations - we may not want
	   to quit here, but we'll do it for now */
	return;
      }

    /*
     * if this nodal contribution has already been added to the filmp
     * equation (i.e. we are at a corner on the second side) return
     * without doing anything
     */
    if (local_node_list_fs[id] == -1) {
      local_node_list_fs[id] = 1;
    } else {
      return;
    }

    /*
     * check to make sure that both lubp dof and sh_fp
     * dof exist at this node
     */
    if (Dolphin[I][R_LUBP] <= 0)     return;
    if (Dolphin[I][R_SHELL_FILMP] <= 0) return;

    /*
     * add local contribution to lubp equation
     * into local contribution for filmp equation while
     * decrementing the filmp equation
     */
    if (af->Assemble_Residual) {
	  ieqn_lubp = R_LUBP;
	  ieqn_filmp = R_SHELL_FILMP;
	  id_doflubp = ei->ln_to_dof[ieqn_lubp][id];
	  id_doffilmp = ei->ln_to_dof[ieqn_filmp][id];
	  lec->R[upd->ep[ieqn_filmp]][id_doffilmp] =
	       -lec->R[upd->ep[ieqn_lubp]][id_doflubp];
    }
    
    /*
     * loop over directions and add local contribution to lubp
     * eqution into local contribution for filmp equation
     */
    if (af->Assemble_Jacobian)
      {
	ieqn_lubp = R_LUBP;
	peqn_lubp = upd->ep[ieqn_lubp];
	ieqn_filmp = R_SHELL_FILMP;
	peqn_filmp = upd->ep[ieqn_filmp];
	id_doflubp = ei->ln_to_dof[ieqn_lubp][id];
	id_doffilmp = ei->ln_to_dof[ieqn_filmp][id];

	/* Add contributions due to all nodal sensitivities in filmp element */

	/*
	 * local J_lubp_d -> J_filmp_d
	 */
	for ( q=0; q<dim; q++)
	  {
	    var = MESH_DISPLACEMENT1+q;
	    if ( pd->v[var] )
	      {
		pvar = upd->vp[var];
		for ( j_id=0; j_id<ei->dof[var]; j_id++)
		  {
		    lec->J[peqn_filmp][pvar][id_doffilmp][j_id] =
		      -lec->J[peqn_lubp][pvar][id_doflubp][j_id];
		  }
	      }
	  }

	/*
	 * local J_lubp_lubp -> J_filmp_lubp
	 */
	var = LUBP;
	if ( pd->v[var] )
	  {
	    pvar = upd->vp[var];
	    for ( j_id=0; j_id<ei->dof[var]; j_id++)
	      {			
		lec->J[peqn_filmp][pvar][id_doffilmp][j_id] =
		  -lec->J[peqn_lubp][pvar][id_doflubp][j_id];
	      }
	  }
				
      } /* end of Jacobian entries */
    
} /* end of routine put_lub_flux_in_film */

/*****************************************************************************/

void 
shell_n_dot_liq_velo_bc_tfmp(double func[DIM],
			     double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
			     const double flowrate, /* imposed flow rate */
			     const double time,     /* current time */
			     const double delta_t,       /* current time step size */
			     double xi[DIM],        /* Local stu coordinates */
			     const Exo_DB *exo)

     /***********************************************************************
      *
      * shell_n_dot_liq_velo_bc_tfmp():
      *
      *  Function which evaluates the expression specifying the
      *  gas velocity at a quadrature point normal to the side
      *  of an element.
      *
      *         func =   - velocity + n .( - h^2 /(12*mu_l)(krl) (grad tfmp_pres) )

      *
      *  The boundary condition SHELL_TFMP_FREE_LIQ_BC employs this function.
      *
      *
      * Input:
      *
      *  velocity       = 0.0 for now (could be specified on the bc card as the first float)
      *  grad tfmp_pres = Lubrication pressure gradient
      *  h              = Film thickness
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
      *   Author: Andrew Cochrane (7/20/2016)
      *   Borrowed heavily from shell_n_dot_flow_bc_film() (KT)
      * 
      ********************************************************************/
{
  int j, k, var;
  int *n_dof = NULL;
  int dof_map[MDE];
  double phi_j;
  double grad_phi_j[DIM], gradII_phi_j[DIM];
  double bound_normal[DIM];

  double S, h, grad_P[DIM], gradII_P[DIM], v_l[DIM];


/* Save the boundary normal vector */

  for(k = 0; k < pd->Num_Dim; k++) {
    bound_normal[k] = fv->snormal[k];
  }

/*
  * Prepare geometry
  */
  n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
  lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

  /* Gather necessary values (S, h, Krg, gradII_P)*/

  // need pure phase viscosity
  dbl mu_l;
  switch(mp->tfmp_viscosity_model){
  case CONSTANT:
    mu_l = mp->tfmp_viscosity_const[1];
    break;
  default:
    // if mp is not set, assemble function will bomb out.
    mu_l = 0.0; //otherwise setup a nofinite
    break;
  }
  
  S = fv->tfmp_sat;
  /* Use the height_function_model */
  double H_U, dH_U_dtime, H_L, dH_L_dtime;
  double dH_U_dX[DIM],dH_L_dX[DIM], dH_U_dp, dH_U_ddh;
  h = height_function_model(&H_U, &dH_U_dtime, &H_L, &dH_L_dtime,
			    dH_U_dX, dH_L_dX, &dH_U_dp, &dH_U_ddh, time, delta_t);

  
  // try shifted-scaled rel perms
  dbl Krl, dKrl_dS;
  dbl Scl, alphal;
  switch (mp->tfmp_rel_perm_model) {
  case PIECEWISE:
    Scl = mp->tfmp_rel_perm_const[2];
    alphal = mp->tfmp_rel_perm_const[3];
    break;
  default:
    Scl = 0.78;
    alphal = 0.1;
    break;
  }
    
  dbl ml = 1./2./alphal;
  dbl cl = -ml*(Scl - alphal);
  
  // liquid transition
  
  if ( S <= Scl - alphal) {
    Krl = 0.0;
    dKrl_dS = 0.0;
  } else if ( S > Scl - alphal && S < Scl + alphal ) {
    Krl = ml*S + cl;
    dKrl_dS = ml;
    
  } else { // S > 1.0
    Krl = 1.0 ;
    dKrl_dS = 0.0;
  }

  for (k = 0; k<DIM; k++) {
    grad_P[k] = fv->grad_tfmp_pres[k];
  }
  Inn(grad_P, gradII_P);
  
  /* Calculate Velocity */
  for (k = 0; k<DIM; k++) {
    v_l[k] = -h*h/12.0/mu_l*Krl*gradII_P[k];
  }
  
  /* Calculate residual contribution */
  func[0] = 0.0;
  for (k = 0; k<pd->Num_Dim; k++) {
    func[0] += v_l[k]*bound_normal[k];
  }
  func[0] *= h;

  /* Calculate Jacobian contributions */
  if (af->Assemble_Jacobian) {
    var = TFMP_PRES;
    if (pd->v[var]) {
      for (j=0; j<ei->dof[var]; j++) {
	phi_j = bf[var]->phi[j];
	for (k=0; k<pd->Num_Dim; k++) {
	  grad_phi_j[k] = bf[var]->grad_phi[j][k];
	  gradII_phi_j[k] = 0.0;
	}
	Inn(grad_phi_j, gradII_phi_j);

	for (k=0; k<pd->Num_Dim; k++) {
	  d_func[0][var][j] += bound_normal[k]*gradII_phi_j[k];
	}
	d_func[0][var][j] *= -h*h/12.0/mu_l*Krl * h;
      }
    }

    var = TFMP_SAT;
    if (pd->v[var]) {
      for (j=0; j<ei->dof[var]; j++) {
	//	d_func[0][var][j] = 0.0;
	phi_j = bf[var]->phi[j];
	for (k=0; k<pd->Num_Dim; k++) {
	  d_func[0][var][j] += bound_normal[k]*gradII_P[k];
	}
	d_func[0][var][j] *= -h*h/12.0/mu_l*dKrl_dS*phi_j * h;
      }
    }
  }
  /* Cleanup */
  safe_free((void *) n_dof);
} /* end of routine shell_n_dot_liq_velo_bc_tfmp */
/*****************************************************************************/

void 
shell_num_diff_bc_tfmp(double func[DIM],
		       double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
		       const double time,     /* current time */
		       const double delta_t,       /* current time step size */
		       double xi[DIM],        /* Local stu coordinates */
		       const Exo_DB *exo)

     /***********************************************************************
      *
      * shell_num_diff_bc_tfmp():
      *
      *  Function which evaluates the expression specifying the
      *  numerical saturation diffusion at a quadrature point normal to 
      *  the side of an element.
      *
      *         func =  + n dot ( (grad S)*phi_i *D*Krd )

      *
      *  The boundary condition SHELL_TFMP_NUM_DIFF_BC employs this function.
      *
      *
      * Input:
      *
      *  grad tfmp_sat = saturation gradient
      *  D             = numerical diffusion coefficient
      *  Krd           = factor describing diffusivity in terms of saturation
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
      *   Author: Andrew Cochrane (7/20/2016)
      *   Borrowed heavily from shell_n_dot_flow_bc_film() (KT)
      *
      *
      *
      ********************************************************************/
{
  int j, k, var;
  int *n_dof = NULL;
  int dof_map[MDE];
  double phi_j;
  double grad_phi_j[DIM], gradII_phi_j[DIM];
  double bound_normal[DIM], bound_normalII[DIM];

  double S, grad_S[DIM], gradII_S[DIM];
  double gradS_dot_n, gradphi_j_dot_n;

/* Save the boundary normal vector */

  for(k = 0; k < pd->Num_Dim; k++) {
    bound_normal[k] = fv->snormal[k];
  }
  Inn(bound_normal, bound_normalII);
/*
  * Prepare geometry
  */
  n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
  lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

  /* Gather necessary values (S, D, Krd)*/

  S = fv->tfmp_sat;
  
  //Artificial diffusion constant
  dbl D, Scd, betad, md, cd, Krd, dKrd_dS;
  //D = .00001;
  switch (mp->tfmp_diff_model) {
  case CONSTANT:
    D = mp->tfmp_diff_const[0];
    Krd = 0.0;
    dKrd_dS = 0.0;
    break;
  case PIECEWISE:
    D = mp->tfmp_diff_const[0];
    // diffusion transition
    Scd = mp->tfmp_diff_const[1];
    betad = mp->tfmp_diff_const[2];
    md = 1.f/2.f/betad;
    cd = -md*(Scd - betad);
  
    if ( S < Scd - betad) {
      Krd = 0.0;
      dKrd_dS = 0.0;
    } else { // ( S >= Scd - betad) { // && S <= Scd + betad ) {
      Krd = md*S + cd;
      dKrd_dS = md;
    }
    break;
    
  default:
    D = 0.0;
    Krd = 0.0;
    dKrd_dS = 0.0;
    break;
  }
   
  for (k = 0; k<DIM; k++) {
    grad_S[k] = fv->grad_tfmp_sat[k];
  }
  Inn(grad_S, gradII_S);
  
  /* Calculate residual contribution */
  func[0] = 0.0;
  for (k = 0; k<pd->Num_Dim; k++) {
    func[0] += bound_normalII[k]*gradII_S[k]*D*Krd;
  }

  /* Calculate Jacobian contributions */
  if (af->Assemble_Jacobian) {
    var = TFMP_SAT;
    if (pd->v[var]) {
      for (j=0; j<ei->dof[var]; j++) {
	phi_j = bf[var]->phi[j];
	gradS_dot_n = 0.0;
	gradphi_j_dot_n= 0.0;
	for (k=0; k<pd->Num_Dim; k++) {
	  grad_phi_j[k] = bf[var]->grad_phi[j][k];
	}
	Inn(grad_phi_j, gradII_phi_j);
	for (k=0; k<pd->Num_Dim; k++) {
	  gradS_dot_n += gradII_S[k]*bound_normalII[k];
	  gradphi_j_dot_n += gradII_phi_j[k]*bound_normalII[k];
	}

	d_func[0][var][j] += -(gradphi_j_dot_n*D*Krd
			      + gradS_dot_n*D*dKrd_dS*phi_j);
      }
    }
  }
  /* Cleanup */
  safe_free((void *) n_dof);
} /* end of routine shell_num_diff_bc_tfmp */
/*****************************************************************************/

