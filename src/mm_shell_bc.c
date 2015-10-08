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
      if (pd->v[pg->imtrx][var])
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

  calculate_lub_q_v(R_LUBP, time, dt, xi, exo);


  if (af->Assemble_LSA_Mass_Matrix)
    {
      return;
    }
  

  if (af->Assemble_Jacobian) 
    {
      var = SHELL_FILMP;
      if (pd->v[pg->imtrx][var])
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
      if (pd->v[pg->imtrx][var])
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
      if (pd->v[pg->imtrx][var])
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
      if (pd->v[pg->imtrx][var])
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
      if (pd->v[pg->imtrx][var])
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
      if (pd->v[pg->imtrx][var])
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
      if (pd->v[pg->imtrx][var])
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
      if (pd->v[pg->imtrx][var])
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
	for ( j_id=0; j_id<ei->dof[var]; j_id++)
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
	  id_doflubp = ei->ln_to_dof[ieqn_lubp][id];
	  id_doffilmp = ei->ln_to_dof[ieqn_filmp][id];
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
	id_doflubp = ei->ln_to_dof[ieqn_lubp][id];
	id_doffilmp = ei->ln_to_dof[ieqn_filmp][id];

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
	if ( pd->v[pg->imtrx][var] )
	  {
	    pvar = upd->vp[pg->imtrx][var];
	    for ( j_id=0; j_id<ei->dof[var]; j_id++)
	      {			
		lec->J[peqn_filmp][pvar][id_doffilmp][j_id] =
		  -lec->J[peqn_lubp][pvar][id_doflubp][j_id];
	      }
	  }
				
      } /* end of Jacobian entries */
    
} /* end of routine put_lub_flux_in_film */
/*****************************************************************************/
