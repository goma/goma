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
 * Routines for calculating Boundary Conditions and Adding them to matrix_fill
 */

/*
 *$Id: bc_surfacedomain.c,v 5.1 2007-09-18 18:53:41 prschun Exp $
 */

/* Standard include files */
 
/* GOMA include files */
 
#include "std.h"
#include "rf_fem_const.h"
#include "el_elm.h"
#include "mm_mp_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_mp.h"
#include "mm_mp_structs.h"
#include "mm_eh.h"
#include "bc_surfacedomain.h"
#include "mm_fill_solid.h"
#include "mm_fill_species.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void 
mass_flux_sd_bc(double func[],
		double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
		int wspec, double mass_tran_coeff, double Y_c, double dt,       
		double tt)
    /*************************************************************************
     *
     *  Function which calculates the surface integral mass transfer
     *  between volumetric domains.
     *
     *  ----------------------------------------------------------------------
     *
     *  Functions called:
     *  mass_flux_surf_mtc  -- calculates mass flux for one component
     *
     *  ----------------------------------------------------------------------
     *
     *************************************************************************/
{  
  int j, j_id,w1,dim,kdir;
  int var,jvar;
  double phi_j;
  double Y_w; /* local concentration of current species */

  double vconv[MAX_PDIM]; /*Calculated convection velocity */
  double vconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;

  int err=0;


  /***************************** EXECUTION BEGINS ****************************/

  /* call routine to calculate surface flux of this component and it's
   * sensitivity to all variable types
   */
  
  mass_flux_surf_mtc(mp->mass_flux,mp->d_mass_flux, fv->T, 
		      fv->c, wspec, mass_tran_coeff, Y_c);

  dim   = pd->Num_Dim;
  Y_w = fv->c[wspec];

  /* get the convection velocity (it's different for arbitrary and
     lagrangian meshes) */

  /* get deformation gradient, etc. if lagrangian mesh with inertia */
  if ((pd->MeshMotion == LAGRANGIAN ||
       pd->MeshMotion == DYNAMIC_LAGRANGIAN ) && pd->MeshInertia == 1)
    {
      err = belly_flop(elc->lame_mu);
      EH(err, "error in belly flop");
      if( neg_elem_volume ) return;
    }
  if (mp->PorousMediaType == CONTINUOUS)
    err = get_convection_velocity(vconv, vconv_old, d_vconv, dt, tt);
  EH(err, "Error in calculating effective convection velocity");

  if (af->Assemble_Residual) {
      *func -= mp->mass_flux[wspec];
      
      /* Calculate the residual contribution from convective flux	*/
      for(kdir=0; kdir<dim; kdir++)  {
	  *func += Y_w*vconv[kdir] * fv->snormal[kdir];
      }
  }
  
  if (af->Assemble_Jacobian )  {
      /* sum the contributions to the global stiffness matrix  for Species*/
      
      /*
       * J_s_c
       */
      var=MASS_FRACTION;
      for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
	phi_j = bf[var]->phi[j_id];
	
	for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++ )
	  {
	    d_func[0][MAX_VARIABLE_TYPES + w1][j_id] -= 
	      mp->d_mass_flux[wspec][MAX_VARIABLE_TYPES + w1] 
		* phi_j;
	  }
	
	for(kdir=0; kdir<dim; kdir++) 
	  {
	    d_func[0][MAX_VARIABLE_TYPES + wspec][j_id] += 
	      phi_j*vconv[kdir] * fv->snormal[kdir];

	      for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++ )
		{
		  d_func[0][MAX_VARIABLE_TYPES + w1][j_id] += 
		    Y_w*d_vconv->C[kdir][w1][j_id] * fv->snormal[kdir];
		}
	  }
      }
      
      /*
       * J_s_T
       */
      var=TEMPERATURE;
      if (pd->v[pg->imtrx][var]){
	for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
	  j_id = j;
	  phi_j = bf[var]->phi[j_id];
	  
	  d_func[0][var][j_id] -= mp->d_mass_flux[wspec][TEMPERATURE] 
	    * phi_j;
	  
	  for(kdir=0; kdir<dim; kdir++) 
	    {
	      d_func[0][var][j_id]  += 
		Y_w*d_vconv->T[kdir][j_id] * fv->snormal[kdir];
	    }
	}
      }
      
      /*
       * J_s_d
       */
      for (jvar = 0; jvar < dim; jvar++) 
	{
	  var=MESH_DISPLACEMENT1+jvar;
	  if (pd->v[pg->imtrx][var])
	    {
	      for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) 
		{
		  /*     d( )/dx        */
		  /* additional terms due to convective flux */
		  
		  phi_j = bf[var]->phi[j_id];
		  for(kdir=0; kdir<dim; kdir++) 
		    {
		      d_func[0][var][j_id] += 
			Y_w*vconv[kdir]*fv->dsnormal_dx[kdir][jvar][j_id] +
			Y_w*d_vconv->X[kdir][jvar][j_id]*fv->snormal[kdir];
		    }
		}
	    }
	}
      
      for(jvar=0; jvar<dim; jvar++) 
	{
	  var = VELOCITY1 + jvar;
	  if (pd->v[pg->imtrx][var])
	    {
	      for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) 
		{
		  phi_j = bf[var]->phi[j_id];
		  d_func[0][var][j_id] += 
		    Y_w*d_vconv->v[jvar][jvar][j_id]*fv->snormal[jvar];
		}
	    }
	}
      
    } /* End of if Assemble_Jacobian */
  
} /* END of routine mass_flux_surf_bc         */
