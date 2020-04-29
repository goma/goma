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

/* 
 *$Id: user_mp_gen.c,v 5.2 2010-04-05 16:49:21 prschun Exp $
 */


/* Standard include files */

#include <math.h>

/* GOMA include files */

#include "std.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "el_elm.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_eh.h"
#include "mm_fill_util.h"
#include "user_mp_gen.h"

#define GOMA_USER_MP_GEN_C


/*********** R O U T I N E S   I N   T H I S   F I L E ************************
*
*       NAME            TYPE            CALLED_BY
*    ------------             ---------               --------------
* user_heat_source_gen    int          assemble_energy
* user_viscosity_gen      int          viscosity
*    
******************************************************************************/
/*********** R E C I P E   F O R   U S A G E **********************************
* Each routine in this file corresponds to one material property or source term.
* Each of the routines is responsible for:
*
*     (1)Calculating the value of that material property or source term at the current
*        current gauss integration point.    
*     (2)Calculating the value of all derivatives with respect to all degrees
*        of freedom REQUESTED by the subroutine at each gauss point.
*
* If these functions are not performed, then the routine should stop and return
* an error.  If the routine is used, then the call to  EH (error handler routine) routine
* should be removed.  The following is a well-documented example for the heat source 
* card.  The heat source in this case can depend on anything, including gradients.  
* A complimentary much simpler set of routines for variations without dependencies on 
* gradients can be found in the file "user_mp.c".  
*
*******************************************************************************/
/*
 * Heat Source Model 
 */

/*
 * int usr_heat_source_gen (h, dhdT, dhdX, dhdC, dhdV, param)
 *
 * ------------------------------------------------------------------------------
 * This routine is responsible for filling up the following forces and sensitivities
 * at the current gauss point:
 *     intput:  
 *
 *     output:  h             - heat source
 *              dhdT[j]    - derivative wrt temperature at node j.
 *              dhdC[i][j] - derivative wrt mass frac species i at node j
 *              dhdV[i][j] - derivative wrt velocity component i at node j
 *              dhdX[0][j] - derivative wrt mesh displacement components i at node j
 *
 *   NB: The user need only supply f, dfdT, dfdC, etc....mp struct is loaded up for you
 * ---------------------------------------------------------------------------
 */

int
usr_heat_source_gen(dbl *h,	/* volumetric heat source */
		    dbl dhdT[MDE], /* temperature dependence. */
		    dbl dhdX[DIM][MDE],	/* spatial dependence. */
		    dbl dhdV[DIM][MDE],	/* velocity dependence. */
		    dbl dhdC[MAX_CONC][MDE], /* concentration dependence. */
		    dbl dhdVolt[MDE], /* voltage dependence. */
		    dbl *param,	/* ptr to the user-defined parameter list */
		    dbl time)
{
  int var;
  
  int dim;
  int a, b;

  /*  int eqn */
  /*  int err; */
  /*  dbl X[DIM], T, C[MAX_CONC]; */ /* Convenient local variables */
  /*  dbl rho, Cp, K; */            /* Convenient property names  */

  /*  dbl dCpdT[MDE];*/

  /*  int i; */
  int j;

  /* Begin Execution */
 /**********************************************************/

 /* Comment out our remove this line if using this routine */

    EH(GOMA_ERROR,"No usr_heat_source_gen model defined.");

 /**********************************************************/

  dim   = pd->Num_Dim;
  /*  eqn   = R_ENERGY;	*/

  /**********************************************************/
  
  /***Load up convenient local variables and properties******/
  /*NB This ought to be done once for all fields at gauss pt*/

  /*  
  T = fv->T;                                       
  for(a=0; a<DIM; a++)X[a] = fv->x[a];		   
  for(i=0; i<pd->Num_Species_Eqn; i++) C[i] = fv->c[i];

  rho  = mp->density;
  K   = mp->thermal_conductivity;
  */
  /* Accounting here must be made of potentially variable heat capacity */
  /*
  if (mp->HeatCapacityModel == CONSTANT)
    {  
      Cp = mp->heat_capacity;
    } 
  else if (mp->HeatCapacityModel == USER) 
    {
      err = usr_heat_capacity(mp->u_heat_capacity, time);
      Cp = mp->heat_capacity;
      var = TEMPERATURE;

      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	{
	  dCpdT[j]= mp->d_heat_capacity[var]*bf[var]->phi[j];
	}
    }
  */
  /**********************************************************/

  /*
   * Example:
   *	Source Term is 
   *
   *		vz.dT/dz + J.J/sigma 
   *
   * for an axisymmetric advection-diffusion problem with Joule heating.
   *
   *   if (!strcmp(pd->MaterialName, "alloy625_e") )  
   *     { 
   *       *h +=  -rho*Cp*param[0]*fv->grad_T[0];
   *     }
   *   if (!strcmp(pd->MaterialName, "60_20_20") )  
   *     { 
   *       *h +=  -rho*Cp*param[0]*fv->grad_T[0];
   *     }
   *   if (!strcmp(pd->MaterialName, "alloy625_i") )  
   *     { 
   *       *h +=  -rho*Cp*param[0]*fv->grad_T[0];
   *     }
   */

  /* Now do sensitivies */

  if (pd->v[pg->imtrx][MASS_FRACTION] )
    {
      var = MASS_FRACTION;
      for(a = 0; a<DIM; a++)
	{
	  for (j=0; j<ei[pg->imtrx]->dof[var]; j++)
	    {
	      dhdC[0][j] += 0.;
	    }
	}
    }

  if (pd->v[pg->imtrx][TEMPERATURE] )
    {
      var = TEMPERATURE;
      for (j=0; j<ei[pg->imtrx]->dof[var]; j++)
	{
	  /*
	   * Example: (contd)
	   *
	   *	  if (!strcmp(pd->MaterialName, "alloy625_e") )  
	   *	    {
	   *	      dhdT[j] += -rho* Cp*      param[0]*bf[var]->grad_phi[j][0] 
	   *		- rho*dCpdT[j]*param[0]*fv->grad_T[0];
	   *	    }
	   *	  if (!strcmp(pd->MaterialName, "60_20_20") )  
	   *	    {
	   *	      dhdT[j] += -rho* Cp*      param[0]*bf[var]->grad_phi[j][0]; 
	   *	    }
	   *	  if (!strcmp(pd->MaterialName, "alloy625_i") )  
	   *	    {
	   *	      dhdT[j] += -rho* Cp*      param[0]*bf[var]->grad_phi[j][0] 
	   *		- rho*dCpdT[j]*param[0]*fv->grad_T[0];
	   *	    }
	   */
	}
    }

      if (pd->v[pg->imtrx][MESH_DISPLACEMENT1] )
	{
	  for ( b=0; b<dim; b++)
	    {
	      var = MESH_DISPLACEMENT1+b;
	      for (j=0; j<ei[pg->imtrx]->dof[var]; j++)
		{
		  /*
		   * Example: (contd)
		   *
		   *		  if (!strcmp(pd->MaterialName, "alloy625_e") )
		   *		    {  
		   *		      dhdX[0][j] +=  -rho*Cp*param[0]* 
		   *			fv->d_grad_T_dmesh[0][b][j];
		   *		    }
		   *		  if (!strcmp(pd->MaterialName, "60_20_20") )
		   *		    {  
		   *		      dhdX[0][j] +=  -rho*Cp*param[0]* 
		   *			fv->d_grad_T_dmesh[0][b][j];
		   *		    }
		   *		  if (!strcmp(pd->MaterialName, "alloy625_i") )
		   *		    {  
		   *		      dhdX[0][j] +=  -rho*Cp*param[0]* 
		   *			fv->d_grad_T_dmesh[0][b][j];
		   *		    }
		   */
		}
	    }
	}

    /* Now add source term for Joule heating. (Recall J=-sigma*grad_Voltage
     * and we're abusing c[0] as an electric potential function.)
     *
     *
     *    if (!strcmp(pd->MaterialName, "60_20_20") ) 
     *      { 
     *	for(a = 0; a<dim; a++)
     *	  {
     *	    *h +=  mp->diffusivity[0] * fv->grad_c[0][a] * fv->grad_c[0][a];
     *	  }
     *
     * Sensitivities of this part of the term.
     *
     * For right now you must have constant mp->diffusivity model 
     *
     *     if(mp->DiffusivityModel[0] != CONSTANT)
     *       {
     *	 EH(GOMA_ERROR,"Cannot use Joule heating term with variable conductivity yet");
     *       }
     *
     *	if (pd->v[pg->imtrx][MASS_FRACTION] )
     *	  {
     *	    var = MASS_FRACTION;
     *	    for(a = 0; a<DIM; a++)
     *	      {
     *		for (j=0; j<ei[pg->imtrx]->dof[var]; j++)
     *		  {
     *		    dhdC[0][j] += mp->diffusivity[0]* 
     *		      2.* fv->grad_c[0][a] * bf[var]->grad_phi[j][a];
     *		  }
     *	      }
     *	  }
     *
     *	if (pd->v[pg->imtrx][MESH_DISPLACEMENT1] )
     *	  {
     *	    for ( a=0; a<dim; a++)
     *	      {
     *		for ( b=0; b<dim; b++)
     *		  {
     *		    var = MESH_DISPLACEMENT1+b;
     *		    for (j=0; j<ei[pg->imtrx]->dof[var]; j++)
     *		      {
     *			dhdX[a][j] +=  mp->diffusivity[0]*2.* fv->grad_c[0][a] * 
     *			fv->d_grad_c_dmesh[a][0][b][j];
     *		      }
     *		  }
     *	      }
     *	  }
     *      }
     */

      return(0);		/* everything is "fine" */
} /* End of routine usr_heat_source_gen */
/*****************************************************************************/

/*
 * VISCOSITY 
 */
/*
 * int usr_viscosity_gen ()
 *
 * ------------------------------------------------------------------------------
 * This routine is responsible for filling up the following portions of the mp structure
 * at the current gauss point:
 *     intput:  param - array of constants input on the property card.  
 *
 *     output:  mu       => mp->viscosity - viscosity
 *              dmudT    => mp->d_viscosity[TEMPERATURE] 
 *                                         - derivative wrt temperature.
 *              dmudC    => mp->d_viscosity[MASS_FRACTION][i]
 *                                         - derivative wrt mass frac species i
 *              dmudV[0] => mp->d_viscosity[VELOCITY1]
 *                          mp->d_viscosity[VELOCITY2]
 *                          mp->d_viscosity[VELOCITY3]
 *                                         - derivative wrt velocities
 *              dmudX[0] => mp->d_viscosity[MESH_DISPLACEMENT1]
 *                          mp->d_viscosity[MESH_DISPLACEMENT2]
 *                          mp->d_viscosity[MESH_DISPLACEMENT3]
 *                                         - derivative wrt mesh displacements
 *
 *   NB: The user need only supply mu, dmudT, dmudC, etc....mp struct is loaded up for you
 * --------------------------------------------------------------------------------
 * Example of an everything-dependent viscosity for material "sample"
 * -----------Add these lines just before last section.-------------
 *
 * Simple heat source function for one material 
 *
 *  if (!strcmp(pd->MaterialName, "sample") )   
 *    {
 *	mu = param[0]*T + param[1]*X[0] + param[2]*X[1] + param[3]*C[0] + param[4];
 *     dmudT    = param[0];
 *     dmudX[0] = param[1];
 *     dmudX[1] = param[2];
 *     dmudC[0] = param[3];
 *    } 
 *  else
 *    {
 *	EH(GOMA_ERROR,"No user-defined function for this material");
 *    }
 * ----------------------------------------------------------------------------
 */

int 
usr_viscosity_gen(dbl *mu,
		  dbl gamma_dot[DIM][DIM], /* strain rate tensor             */
		  dbl *d_mu_dgd, /* deriv of viscosity wrt strain rate inv   */
		  dbl d_mu_dv[DIM][MDE],
		  dbl d_mu_dmesh[DIM][MDE],
		  dbl d_mu_dT[MDE],
		  dbl d_mu_dp[MDE],
		  dbl d_mu_dC[MAX_CONC][MDE],
		  dbl *param)	/* user-defined parameter list               */
{
  int a, b;

  int var;
  /* int var_offset; */
  int w;

  int mdofs=0,vdofs;

  int i, j;

  dbl dmudC, dmudT;

  dbl d_gd_dv[DIM][MDE];        /* derivative of strain rate invariant 
				   wrt velocity */ 
  dbl d_gd_dmesh[DIM][MDE];     /* derivative of strain rate invariant 
				   wrt mesh */ 
  dbl gammadot = 0.;	                /* strain rate invariant */ 

  dbl val, val1, val2;
  dbl mu0;
  dbl muinf;
  dbl nexp;
  dbl aexp;
  dbl lambda;

  /* Begin Execution */

 /**********************************************************/

 /* Comment out our remove this line if using this routine */

     EH(GOMA_ERROR,"No user_viscosity_gen model implemented");
 /**********************************************************/

 /*******Add property function and sensitivities here*******/


  /* Example: */
  mu0 = param[0];
  nexp = param[1];
  muinf = param[2];
  lambda = param[3];
  aexp = param[4];
  
  vdofs = ei[pg->imtrx]->dof[VELOCITY1];
  
  if ( pd->e[pg->imtrx][R_MESH1] )
    {
      mdofs = ei[pg->imtrx]->dof[R_MESH1];
    }
  
  calc_shearrate(&gammadot, gamma_dot, d_gd_dv, d_gd_dmesh);
  
  
  if(gammadot != 0.)
    {
      val2 = pow( lambda*gammadot, aexp);
    }
  else
    {
      val2 = 0.;
    }
  val = pow(1. + val2,(nexp-1.)/aexp);
  *mu = muinf + (mu0 - muinf)* val; 
  
  if(gammadot != 0.)
    {
      val = pow( lambda*gammadot, aexp-1.);
    }
  else
    {
      val = 0.;
    }
  val1 = pow(1. + val2,(nexp-1.-aexp)/aexp);
  
  *d_mu_dgd = (mu0 - muinf)* (nexp-1.) * lambda * val * val1 ;
  
  /*
   * d( mu )/dmesh
   */
  if ( pd->e[pg->imtrx][R_MESH1] )
    {
      for ( b=0; b<VIM; b++)
	{
	  for ( j=0; j<mdofs; j++)
	    {
	      if(gammadot != 0.0)
		{
		  d_mu_dmesh [b][j] =
		    *d_mu_dgd * d_gd_dmesh [b][j] ;
		}
	    }
	}
    }

  /*
   * d( mu )/dv
   */
  
  for ( a=0; a<VIM; a++)
    {
      for ( i=0; i<vdofs; i++)
	{
	  if(gammadot != 0.0)
	    {
	      d_mu_dv[a][i] = 
		*d_mu_dgd *  d_gd_dv [a][i] ;
	    }
	}
    }

  /*
   * d( mu )/dT
   */
  var = TEMPERATURE;
  if ( pd->e[pg->imtrx][var] )
    {
      dmudT = 0.; /* change this line for your function */
      
      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	{
	  d_mu_dT[j]= dmudT * bf[var]->phi[j];
	}
    }
  
  
  /*
   * d( mu )/dC
   */
  var = MASS_FRACTION;
  if (pd->v[pg->imtrx][var] )
    {
      dmudC = 0.;   /* change this line for your function */

      for ( w=0; w<pd->Num_Species_Eqn; w++)
	{
	  var = MASS_FRACTION;
	  /*	  var_offset = MAX_VARIABLE_TYPES + w; */
	  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	    {
	      d_mu_dC[w][j] = dmudC * bf[var]->phi[j];
	    }
	}
    }   
  
  return(0);
} /* End of routine usr_viscosity_gen */
/*****************************************************************************/
/* End of routine user_mp_gen.c  */
/*****************************************************************************/
