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

#include <stdio.h>
#include <string.h>
#include <math.h>

/* GOMA include files */

#include "std.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "el_elm.h"
#include "el_geom.h"
#include "rf_bc_const.h"
#include "mm_mp_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_mp_structs.h"
#include "mm_mp.h"
#include "mm_std_models_shell.h"
#include "mm_eh.h"
#include "mm_elem_block_structs.h"
#include "mm_fill_ls.h"
#include "mm_fill_util.h"
#include "mm_shell_util.h"
#include "mm_viscosity.h"
#include "rf_allo.h"
#include "shell_tfmp_struct.h"
#include "exo_struct.h"

#define GOMA_MM_SHELL_BC_C
#include "mm_shell_bc.h"

#include "shell_tfmp_util.h"

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
      if (pd->v[pg->imtrx][var])
      {
        for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
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
    if (pd->e[pg->imtrx][R_SHELL_NORMAL1] && pd->e[pg->imtrx][R_SHELL_NORMAL2] && pd->e[pg->imtrx][R_SHELL_NORMAL3] )
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
         for ( k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++)
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
   if ( (pd->e[pg->imtrx][R_SHELL_NORMAL1]) && (pd->e[pg->imtrx][R_SHELL_NORMAL2]) && (pd->e[pg->imtrx][R_SHELL_NORMAL3]) )
     {
      for ( i = 0; i < dim; i++)
         {
          for ( j = 0; j < dim; j++)
             {
              for ( k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++)
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
              for ( k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++)
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
         for ( k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) 
            {
             jk = dof_map[k];
             D_H_DX[j][jk] += delta(i,j)*(dH_U_dX[i]-dH_L_dX[i])*bf[MESH_DISPLACEMENT1]->phi[k];
             D_H_DX[j][jk] -= fv->dsnormal_dx[i][j][jk] * fv->d_rs[i];
	    }
         for ( k = 0; k < ei[pg->imtrx]->dof[SOLID_DISPLACEMENT1]; k++) 
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
 if ( (pd->e[pg->imtrx][R_SHELL_NORMAL1]) && (pd->e[pg->imtrx][R_SHELL_NORMAL2]) && (pd->e[pg->imtrx][R_SHELL_NORMAL3]) )
   {
    for ( i = 0; i < dim; i++)
       {
        for ( j = 0; j < dim; j++)
           {
            for ( k = 0; k < ei[pg->imtrx]->dof[SHELL_NORMAL1]; k++)
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


  if ( pd->v[pg->imtrx][FILL] )
    {
     load_lsi( ls->Length_Scale );
     load_lsi_derivs();
     Hn = 1.0 - lsi->Hn;
    }

  /* Calculate F sensitivity */
  if ( pd->v[pg->imtrx][FILL] )
    {
     for ( i = 0; i < ei[pg->imtrx]->dof[FILL]; i++)
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
  if ( pd->v[pg->imtrx][FILL] )
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
  if ( pd->v[pg->imtrx][FILL] )  D_CURV_DH = -(cos(M_PIE-dcaU-atan(slopeU)) + cos(M_PIE-dcaL-atan(-slopeL)))/(H*H);


  /* Sensitivity to level set F */
  if ( pd->v[pg->imtrx][FILL] )
    {
     for ( i = 0; i < ei[pg->imtrx]->dof[FILL]; i++)
        {
         for ( j = 0; j < dim; j++)
            {
             D_CURV_DF[i] += sin(dcaU+atan(slopeU))/(H*(1+slopeU*slopeU))*dH_U_dX[j]**lsi->d_normal_dF[j];
             D_CURV_DF[i] += sin(dcaL+atan(slopeL))/(H*(1+slopeL*slopeL))*dH_L_dX[j]**lsi->d_normal_dF[j];
            }
	}
    }

  /* Sensitivity to mesh */
  if ( pd->v[pg->imtrx][FILL] )
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
      if (pd->v[pg->imtrx][var])
        {
         for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
            {
             phi_j = bf[var]->phi[j];

             d_func[0][var][j] = phi_j - mp->surface_tension * Hn * D_CURV_DH * D_H_DP[j];
            }
	}


      /* Sensitivity w.r.t. level set */
      var = FILL;
      if (pd->v[pg->imtrx][var])
        {
         for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
            {
             d_func[0][var][j] = - mp->surface_tension * (D_Hn_DF[j] * CURV + Hn * D_CURV_DF[j]);
            }
	}

      /* Sensitivity w.r.t. mesh */
      if (pd->v[pg->imtrx][MESH_DISPLACEMENT1])
        {
         for (i = 0; i < dim; i++)
            {
             var = MESH_DISPLACEMENT1 + i;

             for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
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
      if (pd->v[pg->imtrx][var])
      {
        for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
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
      if (pd->v[pg->imtrx][var])
      {
        for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
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
      if (pd->v[pg->imtrx][var])
      {
        for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
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
      if (pd->v[pg->imtrx][var])
      {
        for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
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
      if (pd->v[pg->imtrx][var])
      {
        for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
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

   shell_determinant_and_normal(ei[pg->imtrx]->ielem, ei[pg->imtrx]->iconnect_ptr, ei[pg->imtrx]->num_local_nodes,
                                ei[pg->imtrx]->ielem_dim, 1);

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
      if (pd->v[pg->imtrx][var])
      {
        for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
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

  shell_determinant_and_normal(ei[pg->imtrx]->ielem, ei[pg->imtrx]->iconnect_ptr, ei[pg->imtrx]->num_local_nodes,
                               ei[pg->imtrx]->ielem_dim, 1);

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
      if (pd->v[pg->imtrx][var])
      {
        for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
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
      if (pd->v[pg->imtrx][var])
      {
        for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
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

  shell_tangents(t0, t1, dt0_dx, dt1_dx, dt0_dnormal, dt1_dnormal);


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

      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
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

      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
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
     if ( pd->v[pg->imtrx][var] )
       {

        /*** Loop over dimensions of mesh displacement ***/
        for ( b = 0; b < dim; b++)
           {
            var = MESH_DISPLACEMENT1 + b;

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
               {
                d_func[0][var][j] = de_dx[0][b][j] * tx + de_dx[1][b][j] * ty + de_dx[2][b][j] * tz;
               }
           }
       }

     var = SHELL_NORMAL1;
     if ( pd->v[pg->imtrx][var] )
       {

        /*** Loop over dimensions of shell normal ***/
        for ( b = 0; b < dim; b++)
           {
            var = SHELL_NORMAL1 + b;

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
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
	for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
	  {
	    if (pd->v[pg->imtrx][var])
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
	for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
	  {
	    if (pd->v[pg->imtrx][var])
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
    if (!pd->e[pg->imtrx][R_LUBP]) return;   

    
    id_doflubp = ei[pg->imtrx]->ln_to_dof[R_LUBP][id];
    id_doffilmp = ei[pg->imtrx]->ln_to_dof[R_SHELL_FILMP][id];

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
    if (Dolphin[pg->imtrx][I][R_LUBP] <= 0)     return;
    if (Dolphin[pg->imtrx][I][R_SHELL_FILMP] <= 0) return;

    /*
     * add local contribution to lubp equation
     * into local contribution for filmp equation while
     * decrementing the filmp equation
     */
    if (af->Assemble_Residual) {
	  ieqn_lubp = R_LUBP;
	  ieqn_filmp = R_SHELL_FILMP;
	  id_doflubp = ei[pg->imtrx]->ln_to_dof[ieqn_lubp][id];
	  id_doffilmp = ei[pg->imtrx]->ln_to_dof[ieqn_filmp][id];
	  lec->R[upd->ep[pg->imtrx][ieqn_filmp]][id_doffilmp] =
	       -lec->R[upd->ep[pg->imtrx][ieqn_lubp]][id_doflubp];
    }
    
    /*
     * loop over directions and add local contribution to lubp
     * eqution into local contribution for filmp equation
     */
    if (af->Assemble_Jacobian)
      {
	ieqn_lubp = R_LUBP;
	peqn_lubp = upd->ep[pg->imtrx][ieqn_lubp];
	ieqn_filmp = R_SHELL_FILMP;
	peqn_filmp = upd->ep[pg->imtrx][ieqn_filmp];
	id_doflubp = ei[pg->imtrx]->ln_to_dof[ieqn_lubp][id];
	id_doffilmp = ei[pg->imtrx]->ln_to_dof[ieqn_filmp][id];

	/* Add contributions due to all nodal sensitivities in filmp element */

	/*
	 * local J_lubp_d -> J_filmp_d
	 */
	for ( q=0; q<dim; q++)
	  {
	    var = MESH_DISPLACEMENT1+q;
	    if ( pd->v[pg->imtrx][var] )
	      {
		pvar = upd->vp[pg->imtrx][var];
		for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
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
	if ( pd->v[pg->imtrx][var] )
	  {
	    pvar = upd->vp[pg->imtrx][var];
	    for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
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
      *  liquid velocity at a quadrature point normal to the side
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
  int j, k, l, var;
  int *n_dof = NULL;
  int dof_map[MDE];
  double phi_j;
  double grad_phi_j[DIM], gradII_phi_j[DIM];
  double bound_normal[DIM], dbound_normal_dx[DIM][DIM][MDE];
  double S;
  dbl gradII_P[DIM];
  dbl dgradII_P_dmesh[DIM][DIM][MDE];

  // mesh sensitivity dot products
  dbl gradIIP_dot_bound_normal, dgradIIP_dmesh_dot_bound_normal, gradIIP_dot_dbound_normal_dmesh;
  
  /* Save the boundary normal vector */

  memset(bound_normal,0.0, sizeof(double)*DIM);
  memset(dbound_normal_dx, 0.0, sizeof(double)*DIM*DIM*MDE);

  switch(mp->ehl_integration_kind){
    case SIK_XY:;
      for(k = 0; k < pd->Num_Dim; k++) {
        bound_normal[k] = fv->snormal[k];
        for (l = 0; l<DIM; l++) {
          for (j = 0; j <ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; j++) {
            dbound_normal_dx[k][l][j] = fv->dsnormal_dx[k][l][j];
          }
        }
      }
      break;
    default:
    case SIK_S:
      bound_normal[0] = 1.0;
      bound_normal[1] = 0.0;
      bound_normal[2] = 0.0;
      break;
  }

  /*
  * Prepare geometry
  */
  n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
  lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

  /* Gather necessary values (S, h, Krg, gradII_P)*/
  S = fv->tfmp_sat;

  // need pure phase viscosities
  double mu_l, mu_g;

  load_tfmp_viscosity_model(
    &mu_l, &mu_g
  );

  // here dh_dtime is not used, so it doesn't matter what tt is.
  double tt = 1.0;

  GAP_STRUCT gap_v;
  GAP_STRUCT *gap = &gap_v;
  gap->time = time;
  gap->tt = tt;
  gap->delta_t = delta_t;
  gap->n_dof = n_dof;
  gap->dof_map = dof_map;
  load_gap_model(gap);

  double h = gap->h;
  double dh_dmesh[DIM][MDE];
  double dh_dnormal[DIM][MDE];

  for (int k=0; k<DIM; k++) {
    for (int i=0; i<MDE; i++) {
      dh_dmesh[k][i] = gap->dh_dmesh[k][i];
      dh_dnormal[k][i] = gap->dh_dnormal[k][i];
    }
  }

  /* Use the velocity function model */
  double veloU[DIM], veloL[DIM], veloAVG[DIM];
  velocity_function_model(veloU, veloL, time, delta_t);

  for (k=0; k<DIM; k++) {
    veloAVG[k] = (veloU[k] + veloL[k])/2.;
  }
  veloAVG[2] = 0.0;
  double n_dot_v_avg = 0.0;
  for (k=0; k<DIM; k++) {
    n_dot_v_avg += bound_normal[k]*veloAVG[k];
  }

  //  rel perms
  double Krl, dKrl_dS, Krg, dKrg_dS;
  load_relative_permeability_model(S, &Krl, &dKrl_dS, &Krg, &dKrg_dS);
    
  ShellRotate(fv->grad_tfmp_pres, fv->d_grad_tfmp_pres_dmesh, gradII_P, dgradII_P_dmesh, n_dof[MESH_DISPLACEMENT1]);

  // because I didn't fix ShellRotate to work with domain_s integration
  double csigrad[DIM];
  double det_J;
  double d_det_J_dmeshkj[DIM][MDE];

  if (mp->ehl_integration_kind == SIK_S){
    // fill mapping determinate and sensitivity to mesh motion
    detJ_2d_bar(&det_J, d_det_J_dmeshkj);

    double *grad;

    if (pd->Num_Dim == 2 && ei[pg->imtrx]->ielem_type == LINEAR_BAR) {
      // only one dimension to integrate over, s.
      var = TFMP_PRES;
      grad = gradII_P;

      memset (grad, 0.0, sizeof(double)*DIM);
      memset (csigrad, 0.0, sizeof(double)*DIM);
      for (int i=0; i<ei[pg->imtrx]->dof[var]; i++) {
        csigrad[0] += *esp->tfmp_pres[i]*bf[var]->dphidxi[i][0];
      }

      grad[0] = csigrad[0]/det_J;

      for (int k=0; k<DIM; k++) {
        for (int i=0; i<ei[pg->imtrx]->dof[var]; i++){
          dgradII_P_dmesh[0][k][i] = csigrad[0]*(-1.0)/det_J/det_J*d_det_J_dmeshkj[k][i];
        }
      }
    }
  }


  /* Calculate residual contribution */
  for (k = 0; k<DIM; k++) {
    // Pressure gradient driven
    func[0] += -h*h*h/12.0/mu_l*Krl*gradII_P[k]*bound_normal[k];
  }

  double n_dot_grad_P, n_dot_phi_j;
  double d_gradII_phi_j_dmesh[DIM][DIM][MDE];

  /* Calculate Jacobian contributions */
  if (af->Assemble_Jacobian) {
    var = TFMP_PRES;
    if (pd->v[pg->imtrx][var]) {
      for (j=0; j<ei[pg->imtrx]->dof[var]; j++) {
        ShellBF(var, j, &phi_j, grad_phi_j, gradII_phi_j, d_gradII_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map);
	
        n_dot_phi_j = 0.0;
      	for (k=0; k<DIM; k++) {
	        n_dot_phi_j += bound_normal[k]*gradII_phi_j[k];
	      }
        if (S <= 1.0 ){
          // Pressure gradient driven
          d_func[0][var][j] += n_dot_phi_j*(-h*h*h/12.0/mu_l*Krl);
        }
      }
    }

    var = TFMP_SAT;
    if (pd->v[pg->imtrx][var]) {
      for (j=0; j<ei[pg->imtrx]->dof[var]; j++) {
        ShellBF(R_TFMP_MASS, j, &phi_j, grad_phi_j, gradII_phi_j, d_gradII_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map);
        n_dot_grad_P = 0.0;
        for (k=0; k<DIM; k++) {
          n_dot_grad_P += bound_normal[k]*gradII_P[k];
        }
        if (S <= 1.0 ){
          // Pressure gradient driven
          d_func[0][var][j] += n_dot_grad_P	*(-h*h/12.0/mu_l*dKrl_dS*phi_j*h);
        }
      }
    }

    for (l = 0; l<DIM; l++) {
      var = MESH_DISPLACEMENT1 + l;
      if (pd->v[pg->imtrx][var]) {
	      for (j=0; j<ei[pg->imtrx]->dof[var]; j++) {
          gradIIP_dot_bound_normal = 0.0;
          dgradIIP_dmesh_dot_bound_normal = 0.0;
          gradIIP_dot_dbound_normal_dmesh = 0.0;
	  
          for (k = 0; k<DIM; k++) {
            gradIIP_dot_bound_normal += gradII_P[k]*bound_normal[k];
            dgradIIP_dmesh_dot_bound_normal += dgradII_P_dmesh[k][l][j]*bound_normal[k];
            gradIIP_dot_dbound_normal_dmesh += gradII_P[k]*dbound_normal_dx[k][l][j];
          }
          if (S <= 1.0 ){
            // Pressure gradient driven
            d_func[0][var][j] += -3.0*h*h*dh_dmesh[l][j]/12.0/mu_l*Krl*gradIIP_dot_bound_normal;
            d_func[0][var][j] += -h*h*h/12.0/mu_l*Krl*dgradIIP_dmesh_dot_bound_normal;
            d_func[0][var][j] += -h*h*h/12.0/mu_l*Krl*gradIIP_dot_dbound_normal_dmesh;
          }
        }
      }
    }

    for (l = 0; l<DIM; l++) {
      var = SHELL_NORMAL1 + l;
      if (pd->v[pg->imtrx][var]) {
        for (j=0; j<ei[pg->imtrx]->dof[var]; j++) {
          gradIIP_dot_bound_normal = 0.0;
          // AMCTODO - is bound_normal sensitive to normal?
          // for ds: no
          for (k = 0; k<DIM; k++) {
            gradIIP_dot_bound_normal += gradII_P[k]*bound_normal[k];
          }
          if (S <= 1.0 ){
            // Pressure gradient driven
            d_func[0][var][j] += -3.0*h*h*dh_dnormal[l][j]/12.0/mu_l*Krl*gradIIP_dot_bound_normal;
          }
        }
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

  double S, gradII_S[DIM];
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
  double D, Krd, dKrd_dS;
  load_molecular_diffusion_model(S, &D, &Krd, &dKrd_dS);

  Inn(fv->grad_tfmp_sat, gradII_S);
  
  /* Calculate residual contribution */
  func[0] = 0.0;
  for (k = 0; k<pd->Num_Dim; k++) {
    func[0] += bound_normalII[k]*gradII_S[k]*D*Krd;
  }

  double d_gradII_phi_j_dmesh[DIM][DIM][MDE];

  /* Calculate Jacobian contributions */
  if (af->Assemble_Jacobian) {
    var = TFMP_SAT;
    if (pd->v[pg->imtrx][var]) {
      for (j=0; j<ei[pg->imtrx]->dof[var]; j++) {
        ShellBF( var, j, &phi_j, grad_phi_j, gradII_phi_j, d_gradII_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
        gradS_dot_n = 0.0;
        gradphi_j_dot_n= 0.0;
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

void 
shell_tfmp_avg_plate_velo_liq(double func[DIM],
			     double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
			     const double time,     /* current time */
			     const double delta_t,       /* current time step size */
			     double xi[DIM],        /* Local stu coordinates */
			     const Exo_DB *exo) 

     /***********************************************************************
      *
      * shell_tfmp_avg_plate_velo():
      *
      *  Function which evaluates the expression specifying the
      *  liquid velocity at a quadrature point normal to the side
      *  of an element.
      *
      *   func = h*S*n dot v_avg
      *
      *
      *  The boundary condition SHELL_TFMP_AVG_PLATE_VELO_BC employs this
      *  function.
      *
      *
      * Input:
      *
      *  v_avg   = comes from mp?
      *  S       = saturation (fv->tfmp_sat)
      *  h       = Film thickness
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
      *   Author: Andrew Cochrane (2/20/2018)
      *   
      ********************************************************************/
{
  int j, k, l, var;
  int *n_dof = NULL;
  int dof_map[MDE];
  double phi_j;
  double bound_normal[DIM], dbound_normal_dx[DIM][DIM][MDE];

  double S, h;
  double n_dot_v_avg;
  /* Save the boundary normal vector */

  memset(bound_normal,0.0, sizeof(double)*DIM);
  memset(dbound_normal_dx, 0.0, sizeof(double)*DIM*DIM*MDE);

  switch(mp->ehl_integration_kind){
    case SIK_XY:;
      for(k = 0; k < pd->Num_Dim; k++) {
        bound_normal[k] = fv->snormal[k];
        for (l = 0; l<DIM; l++) {
          for (j = 0; j <ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; j++) {
            dbound_normal_dx[k][l][j] = fv->dsnormal_dx[k][l][j];
          }
        }
      }
      break;
    default:
    case SIK_S:
      bound_normal[0] = 1.0;
      bound_normal[1] = 0.0;
      bound_normal[2] = 0.0;
      break;
  }

  /*
  * Prepare geometry
  */
  n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
  lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

  /* Gather necessary values (S, h, v_avg)*/
  
  S = fv->tfmp_sat;
  /* Use the height_function_model */
  double H_U, dH_U_dtime, H_L, dH_L_dtime;
  double dH_U_dX[DIM],dH_L_dX[DIM], dH_U_dp, dH_U_ddh;

  h = height_function_model(&H_U, &dH_U_dtime, &H_L, &dH_L_dtime,
			    dH_U_dX, dH_L_dX, &dH_U_dp, &dH_U_ddh, time, delta_t);

  double dh_dmesh[DIM][MDE];
  double dh_dnormal[DIM][MDE];
  double d2h_dtime_dmesh[DIM][MDE];
  double d2h_dtime_dnormal[DIM][MDE];
	double d_gradIIh_dmesh[DIM][DIM][MDE];
	double d_gradIIh_dnormal[DIM][DIM][MDE];

  double gradII_h[DIM];

	for (k = 0; k<DIM; k++) {
		gradII_h[k] = dH_U_dX[k] - dH_L_dX[k];
	}


  // here dh_dtime is not used, so it doesn't matter what tt is.
  double tt = 1.0;
  double dh_dtime;
  load_displacement_coupling_model(
    tt,
    delta_t,
    &h,
    &dh_dtime,
    gradII_h,
    dh_dmesh,
    dh_dnormal,
    d2h_dtime_dmesh,
    d2h_dtime_dnormal,
    d_gradIIh_dmesh,
    d_gradIIh_dnormal,
    n_dof,
    dof_map
  );

  double v_avg[DIM], veloU[DIM], veloL[DIM];
  //double liquid_flux_density[DIM];

  velocity_function_model(veloU, veloL, time, delta_t);

  for (k=0; k<ei[pg->imtrx]->ielem_dim; k++) {
    v_avg[k] = (veloU[k] + veloL[k])/2.;
  }
  v_avg[2] = 0.0;

  n_dot_v_avg = 0.0;
  for (k=0; k<DIM; k++) {
    //liquid_flux_density[k] = h*v_avg[k]*S;
    n_dot_v_avg += bound_normal[k]*v_avg[k];
  }

  func[0] += h*S*n_dot_v_avg;


  if (af->Assemble_Jacobian) {
    
    var = TFMP_SAT;
    if (pd->v[pg->imtrx][var]) {
      for (j=0; j<ei[pg->imtrx]->dof[var]; j++) {
        phi_j = bf[var]->phi[j];
        for (k=0; k<DIM; k++) {
          d_func[0][var][j] += bound_normal[k]*v_avg[k]*h*phi_j;
        }
      }
    }

    for (l=0; l<DIM; l++) {
      var = MESH_DISPLACEMENT1 + l;
      if (pd->v[pg->imtrx][var]) {

        for (j=0; j<ei[pg->imtrx]->dof[var]; j++) {
          d_func[0][var][j] += n_dot_v_avg*S*dh_dmesh[l][j];
        }
      }
      var = SHELL_NORMAL1 + l;
      if (pd->v[pg->imtrx][var]) {
        for (j=0; j<ei[pg->imtrx]->dof[var]; j++) {
          d_func[0][var][j] += n_dot_v_avg*S*dh_dnormal[l][j];
        }
      }
    }  
  }
  /* Cleanup */
  safe_free((void *) n_dof);
}/* end of routine shell_tfmp_avg_plate_velo_liq */

void 
shell_tfmp_n_dot_grad_s(double func[DIM],
			     double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
			     const double time,     /* current time */
			     const double delta_t,       /* current time step size */
			     double xi[DIM],        /* Local stu coordinates */
			     const Exo_DB *exo) 

     /***********************************************************************
      * AMCTODO - rewrite this description
      * shell_tfmp_avg_plate_velo():
      *
      *  Function for strong integration of n_dot_grad_s
      *  
      *   func = n dot grad(S)
      *
      *
      *  The boundary condition SHELL_TFMP_AVG_PLATE_VELO_BC employs this
      *  function.
      *
      *
      * Input:
      *
      *  S       = saturation (fv->tfmp_sat)
      *  n       = outflow normal vector
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
      *   Author: Andrew Cochrane (2/27/2018)
      *   
      ********************************************************************/
{
  int j, k, l, var;
  int *n_dof = NULL;
  int dof_map[MDE];
  double phi_j;
  double grad_phi_j[DIM], gradII_phi_j[DIM];
  double bound_normal[DIM], bound_normalII[DIM];

  double gradII_S[DIM];
  double gradphi_j_dot_n;

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

  Inn(fv->grad_tfmp_sat, gradII_S);
   /* Calculate residual contribution */

  for (k = 0; k<pd->Num_Dim; k++) {
    func[0] += bound_normalII[k]*gradII_S[k];
  }
  double d_gradII_phi_j_dmesh[DIM][DIM][MDE];
  /* Calculate Jacobian contributions */
  if (af->Assemble_Jacobian) {

    var = TFMP_SAT;
    if (pd->v[pg->imtrx][var]) {
      for (j=0; j<ei[pg->imtrx]->dof[var]; j++) {
        ShellBF( var, j, &phi_j, grad_phi_j, gradII_phi_j, d_gradII_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
        gradphi_j_dot_n= 0.0;
        for (k=0; k<DIM; k++) {
          gradphi_j_dot_n += gradII_phi_j[k]*bound_normalII[k];
        }
        d_func[0][var][j] += -(gradphi_j_dot_n);
      }
    }

    for (l=0; l<DIM; l++) {
      var = SHELL_NORMAL1 + l;
      if (pd->v[pg->imtrx][var]) {
        for (j=0; j<ei[pg->imtrx]->dof[var]; j++) {
          ShellBF( var, j, &phi_j, grad_phi_j, gradII_phi_j, d_gradII_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
          for (k=0; k<DIM; k++) {
            d_func[0][var][j] += phi_j*gradII_S[k];
          }
        }
      }
    } 
  }
 /* Cleanup */
  safe_free((void *) n_dof);
}

void 
shell_n_dot_gas_velo_bc_tfmp(double func[DIM],
			     double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
			     const double flowrate, /* imposed flow rate */
			     const double time,     /* current time */
			     const double delta_t,       /* current time step size */
			     double xi[DIM],        /* Local stu coordinates */
			     const Exo_DB *exo)

     /***********************************************************************
      *
      * shell_n_dot_gas_velo_bc_tfmp():
      *
      *  Function which evaluates the expression specifying the
      *  liquid velocity at a quadrature point normal to the side
      *  of an element.
      *
      *         func =   - velocity + n .( - h^2 /(12*mu_g)(krg) (grad tfmp_pres) )

      *
      *  The boundary condition SHELL_TFMP_FREE_GAS_BC employs this function.
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
      *   Author: Andrew Cochrane (3/01/2018)
      * 
      ********************************************************************/
{
  int j, k, l, var;
  int *n_dof = NULL;
  int dof_map[MDE];
  double phi_j;
  double grad_phi_j[DIM], gradII_phi_j[DIM];
  double bound_normal[DIM], dbound_normal_dx[DIM][DIM][MDE];
  double S;
  dbl gradII_P[DIM];
  dbl dgradII_P_dmesh[DIM][DIM][MDE];

  // mesh sensitivity dot products
  dbl gradIIP_dot_bound_normal, dgradIIP_dmesh_dot_bound_normal, gradIIP_dot_dbound_normal_dmesh;
  
  /* Save the boundary normal vector */

  memset(bound_normal,0.0, sizeof(double)*DIM);
  memset(dbound_normal_dx, 0.0, sizeof(double)*DIM*DIM*MDE);

  switch(mp->ehl_integration_kind){
    case SIK_XY:;
      for(k = 0; k < pd->Num_Dim; k++) {
        bound_normal[k] = fv->snormal[k];
        for (l = 0; l<DIM; l++) {
          for (j = 0; j <ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; j++) {
            dbound_normal_dx[k][l][j] = fv->dsnormal_dx[k][l][j];
          }
        }
      }
      break;
    default:
    case SIK_S:
      bound_normal[0] = 1.0;
      bound_normal[1] = 0.0;
      bound_normal[2] = 0.0;
      break;
  }

  /*
  * Prepare geometry
  */
  n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
  lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

  /* Gather necessary values (S, h, Krg, gradII_P)*/
  
  // need pure phase viscosities
  double mu_l, mu_g;

  load_tfmp_viscosity_model(
    &mu_l, &mu_g
  );
  
  S = fv->tfmp_sat;

  // here dh_dtime is not used, so it doesn't matter what tt is.
  double tt = 1.0;

  GAP_STRUCT gap_v;
  GAP_STRUCT *gap = &gap_v;
  gap->time = time;
  gap->tt = tt;
  gap->delta_t = delta_t;
  gap->n_dof = n_dof;
  gap->dof_map = dof_map;
  load_gap_model(gap);

  double h = gap->h;
  double dh_dmesh[DIM][MDE];
  double dh_dnormal[DIM][MDE];
  //double dh_dtime = gap->dh_dtime;
  //double d2h_dtime_dmesh[DIM][MDE];
  //double d2h_dtime_dnormal[DIM][MDE];
  //double gradII_h[DIM];
  //double d_gradIIh_dmesh[DIM][DIM][MDE];
  //double d_gradIIh_dnormal[DIM][DIM][MDE];

  for (int k=0; k<DIM; k++) {
    //gradII_h[k] = gap->gradII_h[k];
    for (int i=0; i<MDE; i++) {
      dh_dmesh[k][i] = gap->dh_dmesh[k][i];
      dh_dnormal[k][i] = gap->dh_dnormal[k][i];
      //d2h_dtime_dmesh[k][i] = gap->d2h_dtime_dmesh[k][i];
      //d2h_dtime_dnormal[k][i] = gap->d2h_dtime_dnormal[k][i];
    }
    /*
    for (int l=0; l<DIM; l++){
      for (int i=0; i<MDE; i++){
        d_gradIIh_dmesh[k][l][i] = gap->d_gradIIh_dmesh[k][l][i];
        d_gradIIh_dnormal[k][l][i] = gap->d_gradIIh_dnormal[k][l][i];
      }
    }
    */
  }

  /* Use the velocity function model */
  double veloU[DIM], veloL[DIM], veloAVG[DIM];
  velocity_function_model(veloU, veloL, time, delta_t);
  for (k=0; k<DIM; k++) {
    veloAVG[k] = (veloU[k] + veloL[k])/2.;
  }
  veloAVG[2] = 0.0;
  double n_dot_v_avg = 0.0;
  for (k=0; k<DIM; k++) {
    n_dot_v_avg += bound_normal[k]*veloAVG[k];
  }

  //  rel perms
  double Krl, dKrl_dS, Krg, dKrg_dS;
  load_relative_permeability_model(S, &Krl, &dKrl_dS, &Krg, &dKrg_dS);

  ShellRotate(fv->grad_tfmp_pres, fv->d_grad_tfmp_pres_dmesh, gradII_P, dgradII_P_dmesh, n_dof[MESH_DISPLACEMENT1]);

  // because I didn't fix ShellRotate to work with domain_s integration
  double csigrad[DIM];
  double det_J;
  double d_det_J_dmeshkj[DIM][MDE];

  if (mp->ehl_integration_kind == SIK_S){
    // fill mapping determinate and sensitivity to mesh motion
    detJ_2d_bar(&det_J, d_det_J_dmeshkj);

    double *grad;

    if (pd->Num_Dim == 2 && ei[pg->imtrx]->ielem_type == LINEAR_BAR) {
      // only one dimension to integrate over, s.
      var = TFMP_PRES;
      grad = gradII_P;

      memset (grad, 0.0, sizeof(double)*DIM);
      memset (csigrad, 0.0, sizeof(double)*DIM);
      for (int i=0; i<ei[pg->imtrx]->dof[var]; i++) {
        csigrad[0] += *esp->tfmp_pres[i]*bf[var]->dphidxi[i][0];
      }

      grad[0] = csigrad[0]/det_J;

      for (int k=0; k<DIM; k++) {
        for (int i=0; i<ei[pg->imtrx]->dof[var]; i++){
          dgradII_P_dmesh[0][k][i] = csigrad[0]*(-1.0)/det_J/det_J*d_det_J_dmeshkj[k][i];
        }
      }
    }
  }

  /* Calculate residual contribution */
  for (k = 0; k<pd->Num_Dim; k++) {
    func[0] += -h*h*h/12.0/mu_g*Krg*gradII_P[k]*bound_normal[k];
  }

  // boundary velocity term
  func[0] += h*(1.0 - S)*n_dot_v_avg;

  // res = -h*h*h/12.0/mu_g*Krg*gradIIP_dot_grad
  double n_dot_grad_P, n_dot_grad_phi_j;
  double d_gradII_phi_j_dmesh[DIM][DIM][MDE];

  /* Calculate Jacobian contributions */
  if (af->Assemble_Jacobian) {
    var = TFMP_PRES;
    if (pd->v[pg->imtrx][var]) {
      for (j=0; j<ei[pg->imtrx]->dof[var]; j++) {
	      ShellBF(var, j, &phi_j, grad_phi_j, gradII_phi_j, d_gradII_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map);
        n_dot_grad_phi_j = 0.0;
        for (k=0; k<DIM; k++) {
          n_dot_grad_phi_j += bound_normal[k]*gradII_phi_j[k];
        }
        // pressure gradient term
        d_func[0][var][j] += n_dot_grad_phi_j*(-h*h*h/12.0/mu_g*Krg);
      }
    }

    var = TFMP_SAT;
    if (pd->v[pg->imtrx][var]) {
      for (j=0; j<ei[pg->imtrx]->dof[var]; j++) {
        ShellBF(var, j, &phi_j, grad_phi_j, gradII_phi_j, d_gradII_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map);
        n_dot_grad_P = 0.0;
        for (k=0; k<DIM; k++) {
          n_dot_grad_P += bound_normal[k]*gradII_P[k];
        }
        // pressure gradient term
        d_func[0][var][j] += n_dot_grad_P*(-h*h*h/12.0/mu_g*dKrg_dS*phi_j);
        // Boundary velocity terms
        d_func[0][var][j] += h*(-1.0)*n_dot_v_avg;
      }
    }

    for (l = 0; l<DIM; l++) {
      var = MESH_DISPLACEMENT1 + l;
      if (pd->v[pg->imtrx][var]) {
        for (j=0; j<ei[pg->imtrx]->dof[var]; j++) {
          gradIIP_dot_bound_normal = 0.0;
          dgradIIP_dmesh_dot_bound_normal = 0.0;
          gradIIP_dot_dbound_normal_dmesh = 0.0;
          
          for (k = 0; k<DIM; k++) {
            gradIIP_dot_bound_normal += gradII_P[k]*bound_normal[k];
            dgradIIP_dmesh_dot_bound_normal += dgradII_P_dmesh[k][l][j]*bound_normal[k];
            gradIIP_dot_dbound_normal_dmesh += gradII_P[k]*dbound_normal_dx[k][l][j];
          }
          // pressure gradient term
          d_func[0][var][j] += -3.0*h*h*dh_dmesh[l][j]/12.0/mu_g*Krg*gradIIP_dot_bound_normal;
          d_func[0][var][j] += -h*h*h/12.0/mu_g*Krg*dgradIIP_dmesh_dot_bound_normal;
          d_func[0][var][j] += -h*h*h/12.0/mu_g*Krg*gradIIP_dot_dbound_normal_dmesh;
          // Boundary velocity terms
          d_func[0][var][j] += dh_dmesh[l][j]*(1.0 - S)*n_dot_v_avg;
        }
      }
	
    }

    for (l = 0; l<DIM; l++) {
      var = SHELL_NORMAL1 + l;
      if (pd->v[pg->imtrx][var]) {
        for (j=0; j<ei[pg->imtrx]->dof[var]; j++) {
          gradIIP_dot_bound_normal = 0.0;
          // AMCTODO - is bound_normal sensitive to normal?
          for (k = 0; k<DIM; k++) {
            gradIIP_dot_bound_normal += gradII_P[k]*bound_normal[k];
          }
          // pressure gradient term
          d_func[0][var][j] += -3.0*h*h*dh_dnormal[l][j]/12.0/mu_g*Krg*gradIIP_dot_bound_normal;
          // Boundary velocity terms
          d_func[0][var][j] += dh_dnormal[l][j]*(1.0 - S)*n_dot_v_avg;
        }
      }
	
    }
    
  }
  /* Cleanup */
  safe_free((void *) n_dof);
} /* end of routine shell_n_dot_gas_velo_bc_tfmp */
/*****************************************************************************/

void
shell_lubrication_outflow(
  double func[DIM],
  double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
  const double time,     /* current time */
  const double delta_t,       /* current time step size */
  double xi[DIM],        /* Local stu coordinates */
  const Exo_DB *exo)
{
  int var;
  // 1d case only right now gradII_P[0] === dP_ds
  // need pure phase viscosities
  double mu_l, mu_g;

  load_tfmp_viscosity_model(
    &mu_l, &mu_g
  );
  /*
  * Prepare geometry
  */
  int *n_dof = NULL;
  int dof_map[MDE];
  n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
  lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);
  /* Use the height_function_model */
  double H_U, dH_U_dtime, H_L, dH_L_dtime;
  double dH_U_dX[DIM],dH_L_dX[DIM], dH_U_dp, dH_U_ddh;

  double h;

  h = height_function_model(&H_U, &dH_U_dtime, &H_L, &dH_L_dtime,
          dH_U_dX, dH_L_dX, &dH_U_dp, &dH_U_ddh, time, delta_t);

  double dh_dmesh[DIM][MDE];
  double dh_dnormal[DIM][MDE];
  double d2h_dtime_dmesh[DIM][MDE];
  double d2h_dtime_dnormal[DIM][MDE];
  double d_gradIIh_dmesh[DIM][DIM][MDE];
  double d_gradIIh_dnormal[DIM][DIM][MDE];

  double gradII_h[DIM];

  for (int k = 0; k<DIM; k++) {
    gradII_h[k] = dH_U_dX[k] - dH_L_dX[k];
  }
  memset (dh_dmesh, 0.0, sizeof(double)*DIM*MDE);
  memset (dh_dnormal, 0.0, sizeof(double)*DIM*MDE);
  memset (d2h_dtime_dmesh, 0.0, sizeof(double)*DIM*MDE);
  memset (d2h_dtime_dnormal, 0.0, sizeof(double)*DIM*MDE);
  memset (d_gradIIh_dmesh, 0.0, sizeof(double)*DIM*DIM*MDE);
  memset (d_gradIIh_dnormal, 0.0, sizeof(double)*DIM*DIM*MDE);
  // here dh_dtime is not used, so it doesn't matter what tt is.
  double tt = 1.0;
  double dh_dtime;
  load_displacement_coupling_model(
    tt,
    delta_t,
    &h,
    &dh_dtime,
    gradII_h,
    dh_dmesh,
    dh_dnormal,
    d2h_dtime_dmesh,
    d2h_dtime_dnormal,
    d_gradIIh_dmesh,
    d_gradIIh_dnormal,
    n_dof,
    dof_map
  );

  double gradII_P[DIM];
  double dgradII_P_dmesh[DIM][DIM][MDE];
  ShellRotate(fv->grad_tfmp_pres, fv->d_grad_tfmp_pres_dmesh, gradII_P, dgradII_P_dmesh, n_dof[MESH_DISPLACEMENT1]);

  double phi_j, grad_phi_j[DIM], gradII_phi_j[DIM], d_gradII_phi_j_dmesh[DIM][DIM][MDE];

  // because I didn't fix ShellRotate to work with domain_s integration
  double csigrad[DIM];
  double det_J;
  double d_det_J_dmeshkj[DIM][MDE];

  //double gradII_curv[DIM];
  //double dgradII_curv_dmesh[DIM][DIM][MDE];

  if (mp->ehl_integration_kind == SIK_S) {
    // fill mapping determinate and sensitivity to mesh motion
    detJ_2d_bar(&det_J, d_det_J_dmeshkj);

    double *grad;

    if (pd->Num_Dim == 2 && ei[pg->imtrx]->ielem_type == LINEAR_BAR) {
      // only one dimension to integrate over, s.
      var = TFMP_PRES;
      grad = gradII_P;

      memset (grad, 0.0, sizeof(double)*DIM);
      memset (csigrad, 0.0, sizeof(double)*DIM);
      for (int i=0; i<ei[pg->imtrx]->dof[var]; i++) {
        csigrad[0] += *esp->tfmp_pres[i]*bf[var]->dphidxi[i][0];
      }

      grad[0] = csigrad[0]/det_J;

      for (int k=0; k<DIM; k++) {
        for (int i=0; i<ei[pg->imtrx]->dof[var]; i++){
          dgradII_P_dmesh[0][k][i] = csigrad[0]*(-1.0)/det_J/det_J*d_det_J_dmeshkj[k][i];
        }
      }

      /*
      var = SHELL_CURVATURE;
      grad = gradII_curv;
      memset (grad, 0.0, sizeof(double)*DIM);
      memset (csigrad, 0.0, sizeof(double)*DIM);
      for (int i=0; i<ei[pg->imtrx]->dof[var]; i++) {
        csigrad[0] += *esp->sh_K[i]*bf[var]->dphidxi[i][0];
      }

      grad[0] = csigrad[0]/det_J;

      for (int k=0; k<DIM; k++) {
        for (int i=0; i<ei[pg->imtrx]->dof[var]; i++){
          dgradII_curv_dmesh[0][k][i] = csigrad[0]*(-1.0)/det_J/det_J*d_det_J_dmeshkj[k][i];
        }
      }
      */
    }
  }

  double scale_factor = 1.0;

  if (af->Assemble_Residual) {
    func[0] +=  scale_factor*h*(-h*h/12.0/mu_l*gradII_P[0]);
    //func[0] = fv->tfmp_pres;
    //func[0] = gradII_P[0];
    func[0] = -h*h*h/12.0/mu_l*gradII_P[0];
  }
  if (af->Assemble_Jacobian) {
    var = TFMP_PRES;

    for (int j=0; j<ei[pg->imtrx]->dof[var]; j++) {
      ShellBF(var, j, &phi_j, grad_phi_j, gradII_phi_j, d_gradII_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map);
      //d_func[0][var][j] += scale_factor*h*(-h*h/12.0f/mu_l*gradII_phi_j[0]);
      d_func[0][var][j] = -h*h*h/12.0/mu_l*gradII_phi_j[0];
    }

    var = SHELL_CURVATURE;
    for (int j=0; j<ei[pg->imtrx]->dof[var]; j++) {
      ShellBF(var, j, &phi_j, grad_phi_j, gradII_phi_j, d_gradII_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map);
      //d_func[0][var][j] += scale_factor*h*(-h*h/12.0/mu_l*fv->sh_tens*gradII_phi_j[0]);
    }
    var = SHELL_TENSION;
    for (int j=0; j<ei[pg->imtrx]->dof[var]; j++) {
      //ShellBF(var, j, &phi_j, grad_phi_j, gradII_phi_j, d_gradII_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map);
      //d_func[0][var][j] += scale_factor*h*(-h*h/12.0/mu_l*phi_j*gradII_curv[0]);
    }
    for (int l=0; l<pd->Num_Dim; l++) {
      var = MESH_DISPLACEMENT1 + l;
      for (int j=0; j<ei[pg->imtrx]->dof[var]; j++) {
        ShellBF(var, j, &phi_j, grad_phi_j, gradII_phi_j, d_gradII_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map);
        //d_func[0][var][j] += -scale_factor*3.0*h*h*dh_dmesh[l][j]/12.0/mu_l*gradII_P[0];
        d_func[0][var][j] += -h*h*dh_dmesh[l][j]/4.0/mu_l*gradII_P[0];
        d_func[0][var][j] += -h*h*h/12.0/mu_l*dgradII_P_dmesh[0][l][j];
      }
      var = SHELL_NORMAL1 + l;
      for (int j=0; j<ei[pg->imtrx]->dof[var]; j++) {
        ShellBF(var, j, &phi_j, grad_phi_j, gradII_phi_j, d_gradII_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map);
        //d_func[0][var][j] += -scale_factor*3.0*h*h*dh_dnormal[l][j]/12.0/mu_l*gradII_P[0];
        d_func[0][var][j] += -h*h*dh_dnormal[l][j]/4.0/mu_l*gradII_P[0];
      }
    }
  }

}/* end of routine shell_lubrication_outflow */


void 
shell_tfmp_avg_plate_velo_gas(double func[DIM],
			     double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
			     const double time,     /* current time */
			     const double delta_t,       /* current time step size */
			     double xi[DIM],        /* Local stu coordinates */
			     const Exo_DB *exo) 

     /***********************************************************************
      *
      * shell_tfmp_avg_plate_velo_gas():
      *
      *  Function which evaluates the expression specifying the
      *  gas flux at a quadrature point normal to the side
      *  of an element.
      *
      *   func = h*(1-S)*n dot v_avg
      *
      *
      *  The boundary condition SHELL_TFMP_AVG_PLATE_VELO_BC employs this
      *  function.
      *
      *
      * Input:
      *
      *  v_avg   = comes from mp?
      *  S       = saturation (fv->tfmp_sat)
      *  h       = Film thickness
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
      *   Author: Andrew Cochrane (2/20/2018)
      *   
      ********************************************************************/
{
  int j, k, l, var;
  int *n_dof = NULL;
  int dof_map[MDE];
  double phi_j;
  double grad_phi_j[DIM], gradII_phi_j[DIM];
  double bound_normal[DIM];//, dbound_normal_dx[DIM][DIM][MDE];

  double S, h, v_avg[DIM], veloU[DIM], veloL[DIM];
  double gas_flux_density[DIM];
  /* Save the boundary normal vector */

  for(k = 0; k < pd->Num_Dim; k++) {
    bound_normal[k] = fv->snormal[k];
    /*
    for (l = 0; l<DIM; l++) {
      for (j = 0; j <ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; j++) {
	      dbound_normal_dx[k][l][j] = fv->dsnormal_dx[k][l][j];
      }
    }
    */
  }

  /*
  * Prepare geometry
  */
  n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
  lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

  /* Gather necessary values (S, h, v_avg)*/
  
  S = fv->tfmp_sat;
  /* Use the height_function_model */
  double H_U, dH_U_dtime, H_L, dH_L_dtime;
  double dH_U_dX[DIM],dH_L_dX[DIM], dH_U_dp, dH_U_ddh;

  h = height_function_model(&H_U, &dH_U_dtime, &H_L, &dH_L_dtime,
			    dH_U_dX, dH_L_dX, &dH_U_dp, &dH_U_ddh, time, delta_t);

  double dh_dmesh[DIM][MDE];
  double dh_dnormal[DIM][MDE];
  double d2h_dtime_dmesh[DIM][MDE];
  double d2h_dtime_dnormal[DIM][MDE];
	double d_gradIIh_dmesh[DIM][DIM][MDE];
	double d_gradIIh_dnormal[DIM][DIM][MDE];

  double gradII_h[DIM];

	for (k = 0; k<DIM; k++) {
		gradII_h[k] = dH_U_dX[k] - dH_L_dX[k];
	}

  // here dh_dtime is not used, so it doesn't matter what tt is.
  double tt = 1.0;
  double dh_dtime;
  load_displacement_coupling_model(
    tt,
    delta_t,
    &h,
    &dh_dtime,
    gradII_h,
    dh_dmesh,
    dh_dnormal,
    d2h_dtime_dmesh,
    d2h_dtime_dnormal,
    d_gradIIh_dmesh,
    d_gradIIh_dnormal,
    n_dof,
    dof_map
  );

  velocity_function_model(veloU, veloL, time, delta_t);

  for (k=0; k<DIM; k++) {
    v_avg[k] = (veloU[k] + veloL[k])/2.;
  }
  v_avg[2] = 0.0;


  for (k=0; k<DIM; k++) {
    gas_flux_density[k] = h*v_avg[k]*(1.0-S);
  }
  
  for (k=0; k<DIM; k++){
    func[0] += gas_flux_density[k]*bound_normal[k];
  }
  
  double d_gradII_phi_j_dmesh[DIM][DIM][MDE];
  if (af->Assemble_Jacobian) {
    var = TFMP_SAT;
    if (pd->v[pg->imtrx][var]) {
      for (j=0; j<ei[pg->imtrx]->dof[var]; j++) {
        ShellBF(var, j, &phi_j, grad_phi_j, gradII_phi_j, d_gradII_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map);
        for (k=0; k<DIM; k++) {
          d_func[0][var][j] += -bound_normal[k]*v_avg[k]*h*phi_j;
        }
      }
    }

    for (l=0; l<DIM; l++) {
      var = MESH_DISPLACEMENT1 + l;
      if (pd->v[pg->imtrx][var]) {
        for (j=0; j<ei[pg->imtrx]->dof[var]; j++) {
          for (k=0; k<DIM; k++) {
            d_func[0][var][j] += bound_normal[k]*v_avg[k]*(1.0 - S)*dh_dmesh[k][j];
          }
        }
      }
      var = SHELL_NORMAL1 + l;
      if (pd->v[pg->imtrx][var]) {
        for (j=0; j<ei[pg->imtrx]->dof[var]; j++) {
          for (k=0; k<DIM; k++) {
            d_func[0][var][j] += bound_normal[k]*v_avg[k]*(1.0 - S)*dh_dnormal[k][j];
          }
        }
      }
    }
  }
  /* Cleanup */
  safe_free((void *) n_dof);
}/* end of routine shell_tfmp_avg_plate_velo_gas */

void 
apply_sdet(double func[DIM],
			     double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
           double xi[DIM],        /* Local stu coordinates */
           const Exo_DB *exo)
     /***********************************************************************
      *
      * apply_sdet():
      *
      *  Function which evaluates the expression specifying the
      *  surface determinant at a quadrature point. sign TBD
      *  
      *
      *   func = 1/2*phi_i*sdet*sdet
      *
      *
      *  The boundary condition SH_SDET_BC employs this
      *  function.
      *
      *
      * Input:
      *
      *  sdet   = surface determinant
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
      *   Author: Andrew Cochrane (4/5/2018)
      *   
      ********************************************************************/
           {
            int j, var;

            int *n_dof = NULL;
            int dof_map[MDE];

            n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
            lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);
            double det_J;
            double d_det_J_dmeshkj[DIM][MDE];
            memset(d_det_J_dmeshkj, 0.0, sizeof(double)*DIM*MDE);
            detJ_2d_bar(&det_J, d_det_J_dmeshkj);


            double factor = 1.0;
            func[0] += factor*0.5*det_J*det_J;

            if (af->Assemble_Jacobian) {
              var = MESH_DISPLACEMENT1;
              if (pd->v[pg->imtrx][var]) {
                for (j=0; j<ei[pg->imtrx]->dof[var]; j++) {
                  d_func[0][var][j] += factor*det_J*d_det_J_dmeshkj[0][j];
                }
              }
              var = MESH_DISPLACEMENT2;
              if (pd->v[pg->imtrx][var]) {
                for (j=0; j<ei[pg->imtrx]->dof[var]; j++) {
                  d_func[0][var][j] += factor*det_J*d_det_J_dmeshkj[1][j];
                }
              }
            }

           }/* end of routine apply_sdet */

void 
apply_sh_weak(double func[DIM],
			     double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
           double xi[DIM],        /* Local stu coordinates */
			     const Exo_DB *exo,
           const double dy_ds)
     /***********************************************************************
      *
      * apply_sh_weak():
      *
      *  Function which evaluates the expression specifying the
      *  surface determinant at a quadrature point. sign TBD
      *  
      *
      *   func = -phi_i*dy_ds
      *
      *
      *  The boundary condition SH_MESH2_WEAK_BC employs this
      *  function.
      *
      *
      * Input:
      *
      *  y   = y(global)-position of mesh
      *  s   = position in mesh basis
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
      *   Author: Andrew Cochrane (4/10/2018)
      *   
      ********************************************************************/
{
  int j, var;

  int i, eqn, node, index;

  int *n_dof = NULL;
  int dof_map[MDE];

  // if value of dy_ds isn't defined on BC input, it is set to 0

  n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
  lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

  // get the gradient of mesh position
  eqn = R_MESH2;
  dbl d_sh_y_dcsi = 0.;
  dbl d_phi_dcsi[MDE];
  double d_sh_y_dcsi_dmesh[DIM][MDE];

  memset(d_sh_y_dcsi_dmesh, 0.0, sizeof(double)*DIM*MDE);

  for (i=0; i< ei[pg->imtrx]->dof[eqn]; i++) {
    node = ei[pg->imtrx]->dof_list[R_MESH2][i];
    index = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + node];
    d_phi_dcsi[i] = bf[eqn]->dphidxi[i][0];
    d_sh_y_dcsi += (Coor[1][index] + *esp->d[1][i]) * d_phi_dcsi[i];

    // don't need to loop over k, it's just sensitive to d[1];
    d_sh_y_dcsi_dmesh[1][i] += d_phi_dcsi[i];
  }

  // get the linear bar determinate of mapping
  //(need to modify for elements with more dim than BAR element)
  double det_J;
  double d_det_J_dmeshkj[DIM][MDE];
  detJ_2d_bar(&det_J, d_det_J_dmeshkj);

  if (af->Assemble_Residual) {
    if (dy_ds == 0.0) {
      func[0] += -d_sh_y_dcsi/det_J;
    } else { // use input value
      func[0] += dy_ds/det_J;
    }
  }
  
  if (af->Assemble_Jacobian) {
    var = MESH_DISPLACEMENT1;
    if (pd->v[pg->imtrx][var]) {
      for (j=0; j<ei[pg->imtrx]->dof[var]; j++) {
        if (dy_ds == 0.0) {
          d_func[0][var][j] += d_sh_y_dcsi/det_J/det_J*d_det_J_dmeshkj[0][j];
        } else { // use input value (not sensitive to anything)
          d_func[0][var][j] += -1.0*dy_ds/det_J/det_J*d_det_J_dmeshkj[0][j];
        }

      }
    }
    var = MESH_DISPLACEMENT2;
    if (pd->v[pg->imtrx][var]) {
      for (j=0; j<ei[pg->imtrx]->dof[var]; j++) {
        if (dy_ds == 0.0) {
          d_func[0][var][j] +=  d_sh_y_dcsi/det_J/det_J*d_det_J_dmeshkj[1][j];
          d_func[0][var][j] +=  -d_phi_dcsi[j]/det_J;
        } else { // use input value (not sensitive to anything)
          d_func[0][var][j] += -1.0*dy_ds/det_J/det_J*d_det_J_dmeshkj[1][j];
        }
      }
    }
  }

}/* end of routine apply_sdet */
