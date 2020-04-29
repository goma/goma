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

/* auxiliary routines for calculating tensors, vectors, and scalars that and 
 * sundry boundary conditions that are
 * used in solid mechanics for mesh deformation in the global Jacobian & resid */


#include <math.h>
#include <stdio.h>
#include <string.h>

#include "std.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_mp.h"
#include "rf_io.h"
#include "el_elm.h"
#include "el_geom.h"
#include "rf_bc_const.h"
#include "rf_vars_const.h"
#include "mm_mp_const.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_mp_structs.h"
#include "mm_mp.h"
#include "mm_eh.h"
#include "bc_colloc.h"
#include "el_elm_info.h"
#include "exo_struct.h"
#include "mm_elem_block_structs.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_rs.h"
#include "mm_fill_solid.h"
#include "mm_fill_species.h"
#include "mm_fill_stress.h"
#include "mm_fill_terms.h"
#include "mm_fill_util.h"
#include "mm_qtensor_model.h"
#include "rf_element_storage_struct.h"
#include "rf_node_const.h"
#include "user_mp.h"

#define GOMA_MM_FILL_SOLID_C


/*********** R O U T I N E S   I N   T H I S   F I L E *************************
*
*       NAME			TYPE			CALLED_BY
*    ------------             ---------               --------------
*
*  belly_flop()		int		matrix_fill_terms 
*                                        (assemble_mesh and assemble_continuity)
*
*******************************************************************************/

/* belly_flop() -- find strain tensors for all psuedo-solid and solid models
 *
 * Created:	Tue Mar 15 12:29:51 MST 1994 pasacki@sandia.gov
 * Revised: 12/8/94 from beer_belly by RAC
 * Revised: 7/8/96 by RAC to use deformation gradients in Eulerian Frame
 * Revised: 2/20/99 by PRS for EVP implementation
 *
 */

int
belly_flop(dbl mu)		/* elastic modulus (plane stress case) */
{
  int status;
  int i, j, k, v, p, q, a, b;
  int dim;
  int mdof, dofs;

  dbl factor;			/* coefficients in deviatoric strain matrix */
  dbl deform_grad[DIM][DIM];  /* deformation gradient for nonlinear elasticity 
                       is d(deformation)/d(initial mesh coord) + identity */
  dbl deform_grad_old[DIM][DIM];
  dbl invdeform_grad[DIM][DIM]; 
  dbl d_invdeform_grad_dx[DIM][DIM] [DIM][MDE]; 
  dbl grad_d[DIM][DIM];  /* displacement gradient*/
  dbl d_grad_d[DIM][DIM] [DIM][MDE];  /* displacement gradient*/
  dbl grad_d_old[DIM][DIM];
  dbl det2d;   /* determinant of 2D deformation gradient tensor */
  dbl det2d_old;
  dbl ddet2d_dx[DIM][MDE];  /* sensitivity */
  dbl cauchy_green[DIM][DIM];  /* strain tensor without division by determinant, etc. */
  dbl d_cauchy_green_dx[DIM][DIM][DIM][MDE];  /* sensitivity */
  dbl cauchy_green_old[DIM][DIM];
  static int is_initialized=FALSE;
  
  struct Basis_Functions  *bfv;

  status = 0;

  dim = ei[pg->imtrx]->ielem_dim;
  mdof = ei[pg->imtrx]->dof[MESH_DISPLACEMENT1];
  /*******************************************************************************/
  /* load up the displacement gradient and it's sensitivities or calculate it in 
   * psuedo-cartesian coordinates if using arbitrary mesh */
  /*******************************************************************************/
  /* Initialize arrays used here */
  memset( deform_grad, 0, sizeof(double)*DIM*DIM);
  memset( deform_grad_old, 0, sizeof(double)*DIM*DIM);
  memset( grad_d, 0, sizeof(double)*DIM*DIM);
  memset( grad_d_old, 0, sizeof(double)*DIM*DIM);
  memset( cauchy_green, 0, sizeof(double)*DIM*DIM);
  memset( cauchy_green_old, 0, sizeof(double)*DIM*DIM);
  memset( fv->strain, 0, sizeof(double)*DIM*DIM);
  memset( fv_old->strain, 0, sizeof(double)*DIM*DIM);
  memset( fv->deform_grad, 0, sizeof(double)*DIM*DIM);
  if (af->Assemble_Jacobian) {
    memset(d_grad_d, 0, sizeof(double)*DIM*DIM*DIM*MDE);
    memset(d_cauchy_green_dx, 0, sizeof(double)*DIM*DIM*DIM*MDE);
    memset(ddet2d_dx, 0, sizeof(double)*DIM*MDE);
	
	if( !is_initialized ) {
		memset(fv->d_volume_change_dx, 0, sizeof(double)*DIM*MDE);
		memset(fv->d_volume_strain_dx, 0, sizeof(double)*DIM*MDE);
		memset(fv->d_volume_change_dp, 0, sizeof(double)*MDE);
		memset(fv->d_volume_strain_dp, 0, sizeof(double)*MDE);
		memset(fv->d_strain_dx, 0, sizeof(double)*DIM*DIM*DIM*MDE);
		memset(fv->d_strain_dp, 0, sizeof(double)*DIM*DIM*MDE);
		memset(fv->d_deform_grad_dx, 0, sizeof(double)*DIM*DIM*DIM*MDE); 
		memset(d_invdeform_grad_dx, 0, sizeof(double)*DIM*DIM*DIM*MDE);
		is_initialized = TRUE;
	}
  }

  /*
   * grad(d)
   */
  /* For ARBITRARY MESH and non-cartesian coordinates
   * use Psuedo-Cartesian Coordinates for mesh motion */
  if (cr->MeshMotion == ARBITRARY && pd->CoordinateSystem != CARTESIAN)
    {
      for ( p=0; p<dim; p++)
	{
	  for ( q=0; q<dim; q++)
	    {
	      v = MESH_DISPLACEMENT1 + p;
	      if ( pd->v[pg->imtrx][v] )
		{
		  dofs     = ei[pg->imtrx]->dof[v];
		  for ( i=0; i<dofs; i++)
		    {
		      grad_d[p][q] += 
			*esp->d[q][i] * bf[v]->d_phi[i][p];
		      grad_d_old[p][q] += 
			*esp_old->d[q][i] * bf[v]->d_phi[i][p];
		      
		    }
		} else EH(GOMA_ERROR,"Cant get deformation gradient without mesh!");
	    }
	}

      if (af->Assemble_Jacobian) {
		
	/* calculate sensitivity of displacement gradient */
#ifdef DO_NO_UNROLL
	for ( b=0; b<dim; b++)
	  {
	    v = MESH_DISPLACEMENT1 + b;
	    if ( pd->v[pg->imtrx][v] )
	      {
		for ( j=0; j<ei[pg->imtrx]->dof[MESH_DISPLACEMENT1+b]; j++)
		  {
		    for ( p=0; p<dim; p++)
		      {
			d_grad_d[p][b][b][j] += bf[v]->d_phi[j][p];
		      }
		    for ( i=0; i<ei[pg->imtrx]->dof[v]; i++)
		      {
			for ( p=0; p<dim; p++)
			  {
			    for ( q=0; q<dim; q++)
			      {
				d_grad_d[p][q][b][j] += *esp->d[q][i] * bf[v]->d_d_phi_dmesh[i][p][b][j];
			      }
			  }
		      }
		  }
	      }
	  }
#else
	for ( b=0; b<dim; b++)
	  {
	    v = MESH_DISPLACEMENT1 + b;
	    if ( pd->v[pg->imtrx][v] )
	      {
			bfv = bf[v];
		for ( j=0; j<ei[pg->imtrx]->dof[MESH_DISPLACEMENT1+b]; j++)
		  {
			d_grad_d[0][b][b][j] = bfv->d_phi[j][0];
			d_grad_d[1][b][b][j] = bfv->d_phi[j][1];
			
			if ( dim == 3 )
			 d_grad_d[2][b][b][j] = bfv->d_phi[j][2];

		    for ( i=0; i<ei[pg->imtrx]->dof[v]; i++)
			{
				d_grad_d[0][0][b][j] += *esp->d[0][i] * bfv->d_d_phi_dmesh[i][0][b][j];
				d_grad_d[1][1][b][j] += *esp->d[1][i] * bfv->d_d_phi_dmesh[i][1][b][j];
				d_grad_d[1][0][b][j] += *esp->d[0][i] * bfv->d_d_phi_dmesh[i][1][b][j];
				d_grad_d[0][1][b][j] += *esp->d[1][i] * bfv->d_d_phi_dmesh[i][0][b][j];
				
				if ( dim == 3 )
				{
					d_grad_d[2][2][b][j] += *esp->d[2][i] * bfv->d_d_phi_dmesh[i][2][b][j];
					d_grad_d[2][1][b][j] += *esp->d[1][i] * bfv->d_d_phi_dmesh[i][2][b][j];
					d_grad_d[2][0][b][j] += *esp->d[0][i] * bfv->d_d_phi_dmesh[i][2][b][j];
					d_grad_d[1][2][b][j] += *esp->d[2][i] * bfv->d_d_phi_dmesh[i][1][b][j];
					d_grad_d[0][2][b][j] += *esp->d[2][i] * bfv->d_d_phi_dmesh[i][0][b][j];
				}
			}
		  }
	      }
	  }
#endif


      } /* end of if Assemble_Jacobian */
    } else {  /* Lagrangian Mesh!! */
      /* Use Coordinates of Main Problem */
      for ( p=0; p<VIM; p++)
	{
	  for ( q=0; q<VIM; q++)
	    {
	      grad_d[p][q] = fv->grad_d[p][q];
	      grad_d_old[p][q] = fv_old->grad_d[p][q];
	    }
	}
      if (af->Assemble_Jacobian) {
	for ( p=0; p<VIM; p++)
	  {
	    for ( q=0; q<VIM; q++)
	      {
		for ( b=0; b<dim; b++)
		  {
		    for ( j=0; j<ei[pg->imtrx]->dof[MESH_DISPLACEMENT1+b]; j++)
		      {
			d_grad_d[p][q][b][j] = fv->d_grad_d_dmesh[p][q][b][j];			
		      }
		  }
	      }
	  }
      } /* end of Assemble_Jacobian */
    } /* end of if Arbitrary or Lagrangian */

  /*******************************************************************************/
  /* calculate basic Cauchy-Green strain tensor (grad_d + grad_d_transpose - grad_d*grad_d_transpose)
   * and it's sensitivity */
  /*******************************************************************************/

  /* Start with contribution from linear elasticity */
  for ( p=0; p<VIM; p++)
    {
      for ( q=0; q<VIM; q++)
	{
	  cauchy_green[p][q] = 0.5 * (grad_d[p][q] + grad_d[q][p]);
	  cauchy_green_old[p][q] = 0.5 * (grad_d_old[p][q] + grad_d_old[q][p]);
	}
    }
  /* add on nonlinear term to Eulerian Strain Tensor */
  if (cr->MeshFluxModel != LINEAR) {
    for ( p=0; p<VIM; p++)
      {
	for ( q=0; q<VIM; q++)
	  {
	    for ( a=0; a<VIM; a++)
	      {
		/* dot product between displacement gradient and 
		 * it's transpose - note that there are two ways to 
		 * do this
		 * the right way sums d(d_j)/d(x_m) d(d_j)/d(x_k) 
		 * over j to get B_m_k
		 */
		cauchy_green[p][q] -= 0.5 * grad_d[p][a] * grad_d[q][a];
		cauchy_green_old[p][q] -= 0.5 * grad_d_old[p][a] * grad_d_old[q][a];
	      }
	  }
      }
  }

  if (af->Assemble_Jacobian) {
#if DO_NO_UNROLL
    for ( p=0; p<VIM; p++)
      {
	for ( q=0; q<VIM; q++)
	  {
	    for ( b=0; b<dim; b++)
	      {
		for ( j=0; j<ei[pg->imtrx]->dof[MESH_DISPLACEMENT1+b]; j++)
		  {
		    d_cauchy_green_dx[p][q][b][j] = 0.5 * 
		      (d_grad_d[p][q][b][j] + d_grad_d[q][p][b][j]);
		  }
	      }
	    if (cr->MeshFluxModel != LINEAR) {
	      for ( a=0; a<VIM; a++)
		{
		  for ( b=0; b<dim; b++)
		    {
		      for ( j=0; j<ei[pg->imtrx]->dof[MESH_DISPLACEMENT1+b]; j++)
			{
			  /* dot product between displacement gradient and 
			   * it's transpose - note that there are two ways to 
			   * do this
			   * the right way sums d(d_j)/d(x_m) d(d_j)/d(x_k) 
			   * over j to get B_m_k
			   */
			  d_cauchy_green_dx[p][q][b][j] -= 0.5 * 
			    ( d_grad_d[p][a][b][j] * grad_d[q][a] 
			      + grad_d[p][a] * d_grad_d[q][a][b][j] ); 
			}
		    }
		}
	    }
	  }
      }
#else
	  for ( b=0; b<dim; b++)
	  {
		  for ( j=0; j<ei[pg->imtrx]->dof[MESH_DISPLACEMENT1+b]; j++)
		  {
			  /*d_cauchy_green_dx[p][q][b][j] = 0.5 *  (d_grad_d[p][q][b][j] + d_grad_d[q][p][b][j]);*/

			  d_cauchy_green_dx[0][0][b][j] = 0.5 *  (d_grad_d[0][0][b][j] + d_grad_d[0][0][b][j]);
			  d_cauchy_green_dx[1][1][b][j] = 0.5 *  (d_grad_d[1][1][b][j] + d_grad_d[1][1][b][j]);
			  d_cauchy_green_dx[0][1][b][j] = 0.5 *  (d_grad_d[0][1][b][j] + d_grad_d[1][0][b][j]);
			  d_cauchy_green_dx[1][0][b][j] = 0.5 *  (d_grad_d[1][0][b][j] + d_grad_d[0][1][b][j]);
			  
			  if( VIM == 3 ) 
			  {
				d_cauchy_green_dx[2][2][b][j] = 0.5 *  (d_grad_d[2][2][b][j] + d_grad_d[2][2][b][j]);
				d_cauchy_green_dx[2][1][b][j] = 0.5 *  (d_grad_d[2][1][b][j] + d_grad_d[1][2][b][j]);
				d_cauchy_green_dx[2][0][b][j] = 0.5 *  (d_grad_d[2][0][b][j] + d_grad_d[0][2][b][j]);
				d_cauchy_green_dx[1][2][b][j] = 0.5 *  (d_grad_d[1][2][b][j] + d_grad_d[2][1][b][j]);
				d_cauchy_green_dx[0][2][b][j] = 0.5 *  (d_grad_d[0][2][b][j] + d_grad_d[2][0][b][j]);
			  }
		  }
	  }
	  if (cr->MeshFluxModel != LINEAR) {
	      for ( a=0; a<VIM; a++)
		  {
			  for ( b=0; b<dim; b++)
			  {
				  for ( j=0; j<ei[pg->imtrx]->dof[MESH_DISPLACEMENT1+b]; j++)
				  {
					  /*d_cauchy_green_dx[p][q][b][j] -= 0.5 * ( d_grad_d[p][a][b][j] * grad_d[q][a] + grad_d[p][a] * d_grad_d[q][a][b][j] ); */
					  
					  d_cauchy_green_dx[0][0][b][j] -= 0.5 * ( d_grad_d[0][a][b][j] * grad_d[0][a] + grad_d[0][a] * d_grad_d[0][a][b][j] ); 
					  d_cauchy_green_dx[1][1][b][j] -= 0.5 * ( d_grad_d[1][a][b][j] * grad_d[1][a] + grad_d[1][a] * d_grad_d[1][a][b][j] ); 
					  d_cauchy_green_dx[0][1][b][j] -= 0.5 * ( d_grad_d[0][a][b][j] * grad_d[1][a] + grad_d[0][a] * d_grad_d[1][a][b][j] ); 
					  d_cauchy_green_dx[1][0][b][j] -= 0.5 * ( d_grad_d[1][a][b][j] * grad_d[0][a] + grad_d[1][a] * d_grad_d[0][a][b][j] ); 
					  
					  if( VIM == 3 ) 
					  {
						  d_cauchy_green_dx[2][2][b][j] -= 0.5 * ( d_grad_d[2][a][b][j] * grad_d[2][a] + grad_d[2][a] * d_grad_d[2][a][b][j] ); 
						  d_cauchy_green_dx[2][1][b][j] -= 0.5 * ( d_grad_d[2][a][b][j] * grad_d[1][a] + grad_d[2][a] * d_grad_d[1][a][b][j] ); 
						  d_cauchy_green_dx[2][0][b][j] -= 0.5 * ( d_grad_d[2][a][b][j] * grad_d[0][a] + grad_d[2][a] * d_grad_d[0][a][b][j] ); 
						  d_cauchy_green_dx[1][2][b][j] -= 0.5 * ( d_grad_d[1][a][b][j] * grad_d[2][a] + grad_d[1][a] * d_grad_d[2][a][b][j] ); 
						  d_cauchy_green_dx[0][2][b][j] -= 0.5 * ( d_grad_d[0][a][b][j] * grad_d[2][a] + grad_d[0][a] * d_grad_d[2][a][b][j] );				
					  }
				  }
			  }
		  }
	  }
#endif

  
  }
  /*******************************************************************************/
  /* calculate Lagrangian deformation gradient and it's invariants */
      /* First invariant - trace */
      /* Second invariant - shearing */
      /* Third invariant - volume change */
  /*******************************************************************************/
  if ( pd->MeshMotion == LAGRANGIAN ||
       pd->MeshMotion == DYNAMIC_LAGRANGIAN || 
       pd->MeshMotion == TOTAL_ALE)
    {
      if (cr->MeshFluxModel == LINEAR)
	{
	  /* 
	   * Lagrangian and Eulerian displacement Gradients are the same
	   */
	  for ( p=0; p<VIM; p++)
	    {
	      for ( q=0; q<VIM; q++)
		{
		  fv->deform_grad[p][q] = delta(p,q) + grad_d[p][q];
		}
	    }

	  if (af->Assemble_Jacobian) {
	    for ( p=0; p<VIM; p++)
	      {
		for ( q=0; q<VIM; q++)
		  {
		    for ( b=0; b<dim; b++)
		      {
			for ( j=0; j<ei[pg->imtrx]->dof[MESH_DISPLACEMENT1+b]; j++)
			  {
			    fv->d_deform_grad_dx[p][q][b][j] = d_grad_d[p][q][b][j];
			  }
		      }
		  }
	      }
	  }
	} else {  /* NONLINEAR MESH */
	  /* 
	   * Lagrangian and Eulerian Deformation Gradients are inverses
	   */
	  for ( p=0; p<VIM; p++)
	    {
	      for ( q=0; q<VIM; q++)
		{
		  invdeform_grad[p][q] = delta(p,q) - grad_d[p][q];
		}
	    }

	  if (af->Assemble_Jacobian) {
	    for ( p=0; p<VIM; p++)
	      {
		for ( q=0; q<VIM; q++)
		  {
		    for ( b=0; b<dim; b++)
		      {
			for ( j=0; j<ei[pg->imtrx]->dof[MESH_DISPLACEMENT1+b]; j++)
			  {
			    d_invdeform_grad_dx[p][q][b][j] = - d_grad_d[p][q][b][j];
			  }
		      }
		  }
	      }
	  }
	  invert_tensor(invdeform_grad, fv->deform_grad, VIM,
			d_invdeform_grad_dx, fv->d_deform_grad_dx, 
			ei[pg->imtrx]->dof[MESH_DISPLACEMENT1], af->Assemble_Jacobian);
	}
    }

  /*******************************************************************************/
  /* calculate the volume change needed in Neo-Hookean Constitutive Law 
   * (V - Vo)/Vo */
  /*******************************************************************************/
  fv->volume_change = 0.;
  fv_old->volume_change = 0.;
  fv->volume_strain = 0.;

  if (cr->MeshFluxModel == LINEAR)
    {
      /*
       * Find the trace of the deformation gradient, which in linear elasticity is the 
       * volume change
       */
      fv->volume_change = 1.;
      fv_old->volume_change = 1.;
      fv->volume_strain = 0.;
      for (p=0; p<VIM; p++) 
	{
	  fv->volume_change += cauchy_green[p][p];
	  fv_old->volume_change += cauchy_green_old[p][p];
	  fv->volume_strain += cauchy_green[p][p];
	}

      if (af->Assemble_Jacobian) {


	    for ( b=0; b<dim; b++)
	      {
		for ( j=0; j<ei[pg->imtrx]->dof[MESH_DISPLACEMENT1+b]; j++)
		  {
		  	fv->d_volume_change_dx[b][j] = d_cauchy_green_dx[0][0] [b][j];		    
			fv->d_volume_strain_dx[b][j] = d_cauchy_green_dx[0][0] [b][j];		    

		    fv->d_volume_change_dx[b][j] += d_cauchy_green_dx[1][1] [b][j];		    
		    fv->d_volume_strain_dx[b][j] += d_cauchy_green_dx[1][1] [b][j];	
			if ( VIM == 3) {	    
		    fv->d_volume_change_dx[b][j] += d_cauchy_green_dx[2][2] [b][j];		    
		    fv->d_volume_strain_dx[b][j] += d_cauchy_green_dx[2][2] [b][j];		    
			}
		  }
		
	  }
      }
	  
    } else { /* Non-Linear Mesh Motion */
      /*
       * volume change is the determinant of the Lagrangian deformation gradient tensor
       * or reciprocal of determinant of Eulerian deformation gradient tensor
       */

      /* Eulerian Deformation Gradient is Identity MINUS displacement gradient !! */
      for (p=0; p<VIM; p++) {
	for (q=0; q<VIM; q++) {
	  deform_grad[p][q] = delta(p,q) - grad_d[p][q];
	  deform_grad_old[p][q] = delta(p,q) - grad_d_old[p][q];
        /* fv->deform_grad[p][q] = deform_grad[p][q]; */  /*Uncomment for ST DILATATION MODEL*/
	}
      }

      switch ( dim )
	{
	case 1:
	  fv->volume_change    = 1. / deform_grad[0][0];
	  fv_old->volume_change= 1. / deform_grad_old[0][0];
	  fv->volume_strain    = fv->volume_change - 1.;
	  if ( af->Assemble_Jacobian ) {
	    for (i=0; i<dim; i++) {
	      for (k=0; k<mdof; k++) {
		fv->d_volume_change_dx[i][k] = 1. / deform_grad[0][0] / deform_grad[0][0]
		  * d_grad_d[0][0] [i][k];
		fv->d_volume_strain_dx[i][k] = fv->d_volume_change_dx[i][k];
	      }
	    }
	  }
	  break;

	case 2:
	  /* find determinant of 2-d deformation gradient (note this is not the volume change, that
	     is the determinant of the 3-d deformation gradient which is here approximated by plane
	     strain or plane stress for 2-d) */
	  det2d    = 1. / (
			   deform_grad[0][0] * deform_grad[1][1]
			   - deform_grad[0][1] * deform_grad[1][0]);
	  det2d_old= 1. / (
			   deform_grad_old[0][0] * deform_grad_old[1][1]
			   - deform_grad_old[0][1] * deform_grad_old[1][0]);

	  /* escape if element has inverted */
	  if ( (det2d <= 0.) && (Debug_Flag >= 0 ) ) 
	    {
#ifdef PARALLEL
              fprintf(stderr,"\nP_%d: Volume change  %f\n",ProcID,det2d);
#else
              fprintf(stderr,"\nVolume change  %f\n",det2d);
#endif
              for (i=0; i<dim; i++)
                 {
#ifdef PARALLEL
                  fprintf(stderr,
                     "P_%d: Deformation Gradient %f %f in elem [%d]\n",
                      ProcID, deform_grad[i][0],
                      deform_grad[i][1], ei[pg->imtrx]->ielem );
#else
                  fprintf(stderr,
                     "Deformation Gradient %f %f in elem [%d] %g %g\n",
                      deform_grad[i][0], deform_grad[i][1], ei[pg->imtrx]->ielem,fv->x[0],fv->x[1] );
#endif
                 }
              neg_elem_volume = TRUE;
              return(2);
            }
	  if (det2d <= 0.) 
            {
              neg_elem_volume = TRUE;
              return(2);
            }
	  
	  if ( af->Assemble_Jacobian )
	    {
	      for (i=0; i<dim; i++)
		{
		  for (k=0; k<mdof; k++)
		    {
		      /* N.B. the puzzling 
			 sign change here is absorbed in d_grad_d */
		      ddet2d_dx[i][k] = 
			(
			 deform_grad[0][0] * d_grad_d[1][1] [i][k]
			 - deform_grad[0][1] * d_grad_d[1][0] [i][k]
			 + d_grad_d[0][0] [i][k] * deform_grad[1][1]
			 - d_grad_d[0][1] [i][k] * deform_grad[1][0]
			 ) * det2d * det2d;
		    }
		}
	    }

	  /* PLANE STRAIN CASES */
	  if (cr->MeshFluxModel == NONLINEAR || cr->MeshFluxModel == INCOMP_PSTRAIN 
	      || cr->MeshFluxModel == HOOKEAN_PSTRAIN)
	    {
	      fv->volume_change    = det2d;
	      fv_old->volume_change= det2d_old;
	      fv->volume_strain    = 3. * (pow(det2d, 1./3.) - 1.);

	      if ( af->Assemble_Jacobian )
		{
		  for (i=0; i<dim; i++)
		    {
		      for (k=0; k<mdof; k++)
			{
			  fv->d_volume_change_dx[i][k] = ddet2d_dx[i][k]; 
			  fv->d_volume_strain_dx[i][k] = ddet2d_dx[i][k] * pow(det2d, -2./3.) ;
			}
		    }
		}
	    }
	  /* PLANE STRESS CASES */
	  else if (cr->MeshFluxModel == INCOMP_PSTRESS)
	    {
	      EH(GOMA_ERROR, "need to fix Plane Stress cases");
/* 	      if ((fv->P / mu) >= 1.) EH(GOMA_ERROR, "Zero or Negative Denominator in PLANE STRESS"); */
/* 	      if ( cr->MeshMotion == ARBITRARY ) EH(GOMA_ERROR,"Can't have ARBITRARY mesh with PLANE_STRESS"); */
/* 	      fv->volume_change    =  */
/* 	        pow(det2d, 3./2.) */
/* 		  / pow((1. - fv->P / mu ), 3./4.); */

/* 	      if ( af->Assemble_Jacobian ) */
/* 		{ */
/* 		  for (i=0; i<dim; i++) */
/* 		    { */
/* 		      for (k=0; k<mdof; k++) */
/* 			{ */
/* 			  fv->d_volume_change_dx[i][k] = 3./2.  */
/* 			     * pow(det2d, 1./2.) * ddet2d_dx[i][k] */
/* 			     / pow((1. - fv->P / mu ), 3./4.) ; */
/* 			} */
/* 		    }      */
/* 		  mdofp = ei[pg->imtrx]->dof[R_PRESSURE]; */
/* 		      for (k=0; k<mdofp; k++) */
/* 			{ */
/* 			  fv->d_volume_change_dp[k] =  3./4./mu * bf[PRESSURE]->phi[k] */
/* 			    * pow(det2d, 3./2.) */
/* 			    / pow((1. - fv->P / mu ), 7./4.); */
/* 			} */
/* 		} */
	    }
	  break;
	case 3:
	  fv->volume_change    = 1. / (
				       deform_grad[0][0] * 
				       ( deform_grad[1][1] * deform_grad[2][2] 
					 -deform_grad[1][2] * deform_grad[2][1])
				       - deform_grad[0][1] * 
				       ( deform_grad[1][0] * deform_grad[2][2] 
					 -deform_grad[2][0] * deform_grad[1][2])
				       + deform_grad[0][2] * 
				       ( deform_grad[1][0] * deform_grad[2][1] 
					 -deform_grad[2][0] * deform_grad[1][1]) );
	  fv_old->volume_change= 1. / (
				       deform_grad_old[0][0] * 
				       ( deform_grad_old[1][1] * deform_grad_old[2][2] 
					 -deform_grad_old[1][2] * deform_grad_old[2][1])
				       - deform_grad_old[0][1] * 
				       ( deform_grad_old[1][0] * deform_grad_old[2][2] 
					 -deform_grad_old[2][0] * deform_grad_old[1][2])
				       + deform_grad_old[0][2] * 
				       ( deform_grad_old[1][0] * deform_grad_old[2][1] 
					 -deform_grad_old[2][0] * deform_grad_old[1][1]) );
	  /* Check to make sure element hasn't inverted */
          if ((fv->volume_change <= 0.) && (Debug_Flag >= 0 )) 
            {
#ifdef PARALLEL
              fprintf(stderr,"\nP_%d: Volume change  %f\n",
                      ProcID,fv->volume_change);
#else
              fprintf(stderr,"\nVolume change  %f\n",fv->volume_change);
#endif

              for (i=0; i<dim; i++)
                 {
#ifdef PARALLEL
                  fprintf(stderr,
                     "P_%d: Deformation Gradient %f %f %f in elem [%d]\n",
                      ProcID, deform_grad[i][0], deform_grad[i][1],
                      deform_grad[i][2], ei[pg->imtrx]->ielem );
#else
                  fprintf(stderr,
                     "Deformation Gradient %f %f %f in elem [%d]\n",
                      deform_grad[i][0], deform_grad[i][1],
                      deform_grad[i][2], ei[pg->imtrx]->ielem );
#endif
                 }
                 neg_elem_volume = TRUE;
                 return(2);
            }
	  if (fv->volume_change <= 0.)
            {
              neg_elem_volume = TRUE;
              return(2);
            }

	  fv->volume_strain    = 3. * (pow(fv->volume_change, 1./3.) - 1.);

	  if ( af->Assemble_Jacobian )
	    {
	      for (i=0; i<dim; i++)
		{
		  for (k=0; k<mdof; k++)
		    {
		      fv->d_volume_change_dx[i][k] = 
			(
			 d_grad_d[0][0] [i][k] * 
			   ( deform_grad[1][1] * deform_grad[2][2] 
			     -deform_grad[1][2] * deform_grad[2][1])
			 - d_grad_d[0][1] [i][k] * 
			   ( deform_grad[1][0] * deform_grad[2][2] 
			     -deform_grad[2][0] * deform_grad[1][2])
			 + d_grad_d[0][2] [i][k] * 
			   ( deform_grad[1][0] * deform_grad[2][1] 
			     -deform_grad[2][0] * deform_grad[1][1])
			 + deform_grad[0][0] * 
			   ( d_grad_d[1][1] [i][k] * deform_grad[2][2] 
			     -d_grad_d[1][2] [i][k] * deform_grad[2][1])
			 - deform_grad[0][1] * 
			   ( d_grad_d[1][0] [i][k] * deform_grad[2][2] 
			     -d_grad_d[2][0] [i][k] * deform_grad[1][2])
			 + deform_grad[0][2] * 
			   ( d_grad_d[1][0] [i][k] * deform_grad[2][1] 
			     -d_grad_d[2][0] [i][k] * deform_grad[1][1])
			 + deform_grad[0][0] * 
			   ( deform_grad[1][1] * d_grad_d[2][2] [i][k] 
			     -deform_grad[1][2] * d_grad_d[2][1] [i][k])
			 - deform_grad[0][1] * 
			   ( deform_grad[1][0] * d_grad_d[2][2] [i][k] 
			     -deform_grad[2][0] * d_grad_d[1][2] [i][k])
			 + deform_grad[0][2] * 
			   ( deform_grad[1][0] * d_grad_d[2][1] [i][k] 
			     -deform_grad[2][0] * d_grad_d[1][1] [i][k]) 
			 ) * fv->volume_change * fv->volume_change;
		      fv->d_volume_strain_dx[i][k] = fv->d_volume_change_dx[i][k]
			* pow(fv->volume_change, -2./3.) ;
		    }
		}
	    }
	  break;
	default:
	  EH( -1, "Bad dim.");
	}      
    } /* end of Non-Linear volume change */

  /*******************************************************************************/
  /* convert Cauchy-Green strain tensor to strain tensor for desired consitutive laws */
  /*******************************************************************************/

  if (cr->MeshFluxModel == LINEAR || cr->MeshFluxModel == NONLINEAR 
      || cr->MeshFluxModel == HOOKEAN_PSTRAIN)
    {
#ifdef DO_NO_UNROLL
      for (p=0; p<VIM; p++) 
	{
	  for (q=0; q<VIM; q++) 
	    {
	      fv->strain[p][q] = cauchy_green[p][q];
	      
	      if (af->Assemble_Jacobian) {
		for (i=0; i<dim; i++) 
		  {
		    for (k=0; k<mdof; k++)
		      {
			fv->d_strain_dx[p][q] [i][k] = d_cauchy_green_dx[p][q] [i][k];
		      }
		  }
	      }
	    }
	}
#else
	  fv->strain[0][0] = cauchy_green[0][0];
	  fv->strain[1][1] = cauchy_green[1][1];
	  fv->strain[1][0] = cauchy_green[1][0];
	  fv->strain[0][1] = cauchy_green[0][1];
          
	  fv_old->strain[0][0] = cauchy_green_old[0][0];
	  fv_old->strain[1][1] = cauchy_green_old[1][1];
	  fv_old->strain[1][0] = cauchy_green_old[1][0];
	  fv_old->strain[0][1] = cauchy_green_old[0][1];
	  
	  if( VIM == 3 ) 
	  {
		  fv->strain[2][2] = cauchy_green[2][2];
		  fv->strain[2][1] = cauchy_green[2][1];
		  fv->strain[2][0] = cauchy_green[2][0];
		  fv->strain[1][2] = cauchy_green[1][2];
		  fv->strain[0][2] = cauchy_green[0][2];
                  
		  fv_old->strain[2][2] = cauchy_green_old[2][2];
		  fv_old->strain[2][1] = cauchy_green_old[2][1];
		  fv_old->strain[2][0] = cauchy_green_old[2][0];
		  fv_old->strain[1][2] = cauchy_green_old[1][2];
		  fv_old->strain[0][2] = cauchy_green_old[0][2];
	  }
	  
	  if (af->Assemble_Jacobian) {
		  for (i=0; i<dim; i++) 
		  {
			  for (k=0; k<mdof; k++)
		      {
				  fv->d_strain_dx[0][0] [i][k] = d_cauchy_green_dx[0][0] [i][k];
				  fv->d_strain_dx[1][1] [i][k] = d_cauchy_green_dx[1][1] [i][k];
				  fv->d_strain_dx[0][1] [i][k] = d_cauchy_green_dx[0][1] [i][k];
				  fv->d_strain_dx[1][0] [i][k] = d_cauchy_green_dx[1][0] [i][k];
				  
				  if( VIM == 3)
				  {
					  fv->d_strain_dx[2][2] [i][k] = d_cauchy_green_dx[2][2] [i][k];
					  fv->d_strain_dx[2][1] [i][k] = d_cauchy_green_dx[2][1] [i][k];
					  fv->d_strain_dx[2][0] [i][k] = d_cauchy_green_dx[2][0] [i][k];
					  fv->d_strain_dx[0][2] [i][k] = d_cauchy_green_dx[0][2] [i][k];
					  fv->d_strain_dx[1][2] [i][k] = d_cauchy_green_dx[1][2] [i][k];				
				  }
				  
		      }
		  }
	  }
#endif	  
    } else if (cr->MeshFluxModel == INCOMP_PSTRAIN || cr->MeshFluxModel == INCOMP_3D) {

      /* this is the model of Segalman which removes dilation from the strain tensor */
      factor = 0.5;
      for (p=0; p<dim; p++)
	{
	  for (q=0; q<dim; q++)
	    {
 	      fv->strain[p][q] = factor * delta(p,q) * (1. - pow(fv->volume_change,2./3.))
		+ cauchy_green[p][q] * pow(fv->volume_change,2./3.); 
 	      fv_old->strain[p][q] = factor * delta(p,q) * (1. - pow(fv_old->volume_change,2./3.))
		+ cauchy_green_old[p][q] * pow(fv_old->volume_change,2./3.); 
	    }
	}
      if (af->Assemble_Jacobian) {
	/*
	 * Now find sensitivity with respect to nodal displacements
	 */
	for (i=0; i<dim; i++)
	  {
	    for (k=0; k<mdof; k++)
	      {
		for (p=0; p<dim; p++)
		  {
		    for (q=0; q<dim; q++)
		      {
 			fv->d_strain_dx[p][q] [i][k] = 
			  -2./3. * factor * delta(p,q) * pow(fv->volume_change, -1./3.) 
			    * fv->d_volume_change_dx[i][k] 
 			  + ( d_cauchy_green_dx[p][q][i][k] * pow(fv->volume_change, 2./3.) 
 			      + 2./3. * cauchy_green[p][q] * pow(fv->volume_change, -1./3.)  
			      * fv->d_volume_change_dx[i][k] ); 
		      }
		  }
	      }
	  }
      } 
      /* end of INCOMPRESSIBLE MODELS */
    } else if (cr->MeshFluxModel == INCOMP_PSTRESS || cr->MeshFluxModel == HOOKEAN_PSTRESS) {
      EH(GOMA_ERROR,"Need to fix PSTRESS implementation");
/* 	    factor = 0.5; */
/* 	    for (p=0; p<dim; p++) */
/* 	      { */
/* 		for (q=0; q<dim; q++) */
/* 		  { */
/* 		    fv->strain[p][q] = factor * (delta(p,q) - cauchy_green[p][q] / */
/* 		      det2d / pow((1. - fv->P / mu ), 3./4.)); */
/* 		  } */
/* 	      } */
/* */
/*  * Now find sensitivity with respect to nodal displacements */
/*  */
/* 	    if ( af->Assemble_Jacobian ) */
/* 	      { */
/* 		for (i=0; i<dim; i++) */
/* 		  { */
/* 		    for (k=0; k<mdof; k++) */
/* 		      { */
/* 			for (p=0; p<dim; p++) */
/* 			  { */
/* 			    for (q=0; q<dim; q++) */
/* 			      { */
/* 				fv->d_strain_dx[p][q] [i][k] =  -factor * (dcauchy_green_dx[p][q][i][k] */
/* 				    - cauchy_green[p][q] / det2d * ddet2d_dx[i][k]) */
/* 				     / det2d / pow((1. - fv->P / mu), 3./4.); */
/* 			      } */
/* 			  } */
/* 		      } */
/* 		  } */
/* 		for (k=0; k<mdofp; k++) */
/* 		  { */
/* 		    for (p=0; p<dim; p++) */
/* 		      { */
/* 			for (q=0; q<dim; q++) */
/* 			  { */
/* 			    fv->d_strain_dp[p][q] [k] = - factor * 3./4. * cauchy_green[p][q]  */
/* 				  / det2d / pow((1. - fv->P / mu), 7./4.) / mu * bf[PRESSURE]->phi[k]; */
/* 			  } */
/* 		      } */
/* 		  } */
/* 	      } */
/* 	  }  end of PLANE STRESS  */
    } else EH(GOMA_ERROR,"Illegal Mesh Constitutive Equation");

  return(status);
}     
/*
 * This function inverts a matrix and evaluates it's sensitivities 
 * Written by RAC 1- July 1996
 */
void 
invert_tensor(double A[DIM][DIM], /* tensor to be inverted */
	      double B[DIM][DIM], /* inverted tensor */
	      int dim,		/* dimensions of tensor */
	      double dA[DIM][DIM][DIM][MDE], /* sensitivities of tensor to be 
					      * inverted */
	      double dB[DIM][DIM][DIM][MDE], /* sensitivities of inverted 
					      * tensor */
	      int dof,		/* number of dofs of variable for 
				 * sensitivities */
	      int sense)	/* flag to calculate sensitivities */
{
  double detA, detA_i, ddetA_i;
  int j, p;

  switch ( dim )
    {
    case 1:
      detA    = A[0][0];
      B[0][0] = 1. / detA;

      if (sense && dA != NULL) {
	for (j=0; j<dof; j++) {
	  for (p=0; p<dim; p++) {
	    dB[0][0][p][j] = - B[0][0] * B[0][0] * dA[0][0][p][j];
	  }
	}
      }
      break;
      
    case 2:
      detA    =  A[0][0] * A[1][1]
	- A[0][1] * A[1][0];
      detA_i = 1. / detA;
      
      B[0][0] =  A[1][1] * detA_i;
      B[0][1] = -A[0][1] * detA_i;
      B[1][0] = -A[1][0] * detA_i;
      B[1][1] =  A[0][0] * detA_i;

      if (sense && dA != NULL) {
	for (j=0; j<dof; j++) {
	  for (p=0; p<dim; p++) {
	    ddetA_i = - detA_i * detA_i * 
	      (dA[0][0][p][j] * A[1][1] - dA[0][1][p][j] * A[1][0]
	       + A[0][0] * dA[1][1][p][j] - A[0][1] * dA[1][0][p][j]);
	    dB[0][0][p][j] =   dA[1][1][p][j] * detA_i + A[1][1] * ddetA_i;
	    dB[0][1][p][j] = - dA[0][1][p][j] * detA_i - A[0][1] * ddetA_i;
	    dB[1][0][p][j] = - dA[1][0][p][j] * detA_i - A[1][0] * ddetA_i;
	    dB[1][1][p][j] =   dA[0][0][p][j] * detA_i + A[0][0] * ddetA_i;
	  }
	}
      }
      break;
      
    case 3:
      detA = A[0][0] * ( A[1][1] * A[2][2] -A[1][2] * A[2][1])
	- A[0][1] * ( A[1][0] * A[2][2] -A[2][0] * A[1][2])
	+ A[0][2] * ( A[1][0] * A[2][1] -A[2][0] * A[1][1]);
      detA_i = 1. / detA;
      
      B[0][0] = ( A[1][1] * A[2][2] -A[2][1] * A[1][2] )
	* detA_i;
      
      B[0][1] =-( A[0][1] * A[2][2] -A[2][1] * A[0][2])
	* detA_i;
      
      B[0][2] = ( A[0][1] * A[1][2] -A[1][1] * A[0][2])
	* detA_i;
      
      B[1][0] =-( A[1][0] * A[2][2] -A[2][0] * A[1][2])
	* detA_i;
      
      B[1][1] = ( A[0][0] * A[2][2] -A[2][0] * A[0][2])
	* detA_i;
      
      B[1][2] =-( A[0][0] * A[1][2] -A[1][0] * A[0][2])
	* detA_i;
      
      B[2][0] = ( A[1][0] * A[2][1] -A[1][1] * A[2][0])
	* detA_i;
      
      B[2][1] =-( A[0][0] * A[2][1] -A[2][0] * A[0][1])
	* detA_i;
      
      B[2][2] = ( A[0][0] * A[1][1] -A[1][0] * A[0][1])
	* detA_i;
      
      if (sense && dA != NULL) {
	for (j=0; j<dof; j++) {
	  for (p=0; p<dim; p++) {
	    ddetA_i = - detA_i * detA_i 
	      * (
		   dA[0][0][p][j] * ( A[1][1] * A[2][2] - A[1][2] * A[2][1])
		 - dA[0][1][p][j] * ( A[1][0] * A[2][2] - A[2][0] * A[1][2])
		 + dA[0][2][p][j] * ( A[1][0] * A[2][1] - A[2][0] * A[1][1])
		 + A[0][0] * ( dA[1][1][p][j] * A[2][2] - dA[1][2][p][j] * A[2][1])
		 - A[0][1] * ( dA[1][0][p][j] * A[2][2] - dA[2][0][p][j] * A[1][2])
		 + A[0][2] * ( dA[1][0][p][j] * A[2][1] - dA[2][0][p][j] * A[1][1])
		 + A[0][0] * ( A[1][1] * dA[2][2][p][j] - A[1][2] * dA[2][1][p][j])
		 - A[0][1] * ( A[1][0] * dA[2][2][p][j] - A[2][0] * dA[1][2][p][j])
		 + A[0][2] * ( A[1][0] * dA[2][1][p][j] - A[2][0] * dA[1][1][p][j])
		 );
	    dB[0][0][p][j] = ( A[1][1] * A[2][2] -A[2][1] * A[1][2] )
	      * ddetA_i +  
	      ( dA[1][1][p][j] * A[2][2] -dA[2][1][p][j] * A[1][2] 
		+ A[1][1] * dA[2][2][p][j] -A[2][1] * dA[1][2][p][j] )
	      * detA_i;
      
	    dB[0][1][p][j] = -( A[0][1] * A[2][2] -A[2][1] * A[0][2])
	      * ddetA_i - 
	      ( dA[0][1][p][j] * A[2][2] -dA[2][1][p][j] * A[0][2]
		+ A[0][1] * dA[2][2][p][j] -A[2][1] * dA[0][2][p][j])
	      * detA_i;
	    
	    dB[0][2][p][j] = ( A[0][1] * A[1][2] -A[1][1] * A[0][2])
	      * ddetA_i + 
	      ( dA[0][1][p][j] * A[1][2] -dA[1][1][p][j] * A[0][2] 
		+ A[0][1] * dA[1][2][p][j] -A[1][1] * dA[0][2][p][j])
	      * detA_i;
	    
	    dB[1][0][p][j] = -( A[1][0] * A[2][2] -A[2][0] * A[1][2])
	      * ddetA_i - 
	      ( dA[1][0][p][j] * A[2][2] -dA[2][0][p][j] * A[1][2] 
		+ A[1][0] * dA[2][2][p][j] -A[2][0] * dA[1][2][p][j])
	      * detA_i;
	    
	    dB[1][1][p][j] = ( A[0][0] * A[2][2] -A[2][0] * A[0][2])
	      * ddetA_i + 
	      ( dA[0][0][p][j] * A[2][2] -dA[2][0][p][j] * A[0][2] 
		+ A[0][0] * dA[2][2][p][j] -A[2][0] * dA[0][2][p][j])
	      * detA_i;
	    
	    dB[1][2][p][j] = -( A[0][0] * A[1][2] -A[1][0] * A[0][2])
	      * ddetA_i - 
	      ( dA[0][0][p][j] * A[1][2] -dA[1][0][p][j] * A[0][2] 
		+ A[0][0] * dA[1][2][p][j] -A[1][0] * dA[0][2][p][j])
	      * detA_i;
	    
	    dB[2][0][p][j] = ( A[1][0] * A[2][1] -A[1][1] * A[2][0])
	      * ddetA_i + 
	      ( dA[1][0][p][j] * A[2][1] -dA[1][1][p][j] * A[2][0] 
		+ A[1][0] * dA[2][1][p][j] -A[1][1] * dA[2][0][p][j])
	      * detA_i;
	    
	    dB[2][1][p][j] = -( A[0][0] * A[2][1] -A[2][0] * A[0][1])
	      * ddetA_i - 
	      ( dA[0][0][p][j] * A[2][1] -dA[2][0][p][j] * A[0][1]
		+ A[0][0] * dA[2][1][p][j] -A[2][0] * dA[0][1][p][j])
	      * detA_i;
	    
	    dB[2][2][p][j] = ( A[0][0] * A[1][1] -A[1][0] * A[0][1])
	      * ddetA_i + 
	      ( dA[0][0][p][j] * A[1][1] -dA[1][0][p][j] * A[0][1]
		+ A[0][0] * dA[1][1][p][j] -A[1][0] * dA[0][1][p][j])
	      * detA_i;
	  }
	}
      }
     break;
      
    default:
      EH(GOMA_ERROR,"can't invert a tensor with more than 3 dimensions");
      break;
    }
  return;
}

/*
 * This function sets a boundary condition on the slope of the boundary
 * by setting the residual equal to the normal to the boundary dotted
 * into the unit vector of the desired slope
 *
 */

void 
slope_n_dot_n0_bc(double func[DIM],
		  double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
		  double slopex, 
		  double slopey, 
		  double slopez)
{
    int j_id, jvar, kdir, dim, var;
    double slope[3];

    if(af->Assemble_LSA_Mass_Matrix)
      return;

 /* load vector of the desired slope */
    slope[0] = slopex;
    slope[1] = slopey;
    slope[2] = slopez;

    dim = pd->Num_Dim;
 
  if (af->Assemble_Jacobian) {
    for(jvar=0; jvar<dim; jvar++)
      {
	var = MESH_DISPLACEMENT1 + jvar;
	if (pd->v[pg->imtrx][var])
	  {
	    for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++)
	      {
		for(kdir=0; kdir<dim; kdir++)
		  {
		    var = MESH_DISPLACEMENT1 + jvar;
		    d_func[0][var][j_id] += 
		      slope[kdir] * fv->dsnormal_dx[kdir][jvar][j_id] ;
		  }
	      }
	  }
      }
  }

/* Calculate the residual contribution  */
 
    for(kdir=0; kdir<dim; kdir++)
    {
        *func += slope[kdir] * fv->snormal[kdir];
    }

}
/* end of routine slope_n_dot_n0_bc */

/*
 * This function specifies the external force imposed on a boundary
 * which is balanced by the normal elastic force within the material
 * this force is applied weakly to the boundary term in the integral
 * (by parts) from the Galerkin form of the Lagrangian Mesh Momentum
 * equation
 *
 */

void 
force_n_dot_f_bc(double func[DIM],
		 double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
		 const double forcex, 
		 const double forcey,
		 const double forcez,
		 const int sic_flag,
		 const double delta_t,
		 const double theta,
		 const int ip,
		 const int ip_total,
		 const double time_value)
{
    int kdir, dim;
    double force[3], phi_j;
  int j_id, jvar, var;
  int err, a, b, c, j;
  int eqn;
  double TT[MAX_PDIM][MAX_PDIM];   /**  solid stresses  **/
  double dTT_drs[DIM][DIM][DIM][MDE];
  double dTT_dx[MAX_PDIM][MAX_PDIM][MAX_PDIM][MDE];
  double dTT_dp[MAX_PDIM][MAX_PDIM][MDE];
  double dTT_dc[MAX_PDIM][MAX_PDIM][MAX_CONC][MDE];
  double dTT_dp_liq[DIM][DIM][MDE];/* Sensitivity of stress tensor...
                                    to nodal porous liquid pressure*/
  double dTT_dp_gas[DIM][DIM][MDE];/* Sensitivity of stress tensor...
                                    to nodal porous gas pressure*/
  double dTT_dporosity[DIM][DIM][MDE];/* Sensitivity of stress tensor...
                                    to nodal porosity*/
  double dTT_dsink_mass[DIM][DIM][MDE];/* Sensitivity of stress tensor...
                                    to sink_mass*/
  double dTT_dT[MAX_PDIM][MAX_PDIM][MDE];
  double dTT_dmax_strain[MAX_PDIM][MAX_PDIM][MDE];
  double dTT_dcur_strain[MAX_PDIM][MAX_PDIM][MDE];
  double elast_modulus;

  double vconv[MAX_PDIM], vconv_old[MAX_PDIM];/*Calculated convection velocity*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;

    if(af->Assemble_LSA_Mass_Matrix)
      return;

 /* load vector of the desired force */
    force[0] = forcex;
    force[1] = forcey;
    force[2] = forcez;

/*  initialize variables */
    dim = pd->Num_Dim;

/* Calculate the residual contribution  */
 
if(sic_flag != 2)
  {
    for(kdir=0; kdir<dim; kdir++)
    {
          func[kdir] += force[kdir];
    }
  }

/*  If the strong integrated form add the traction force	*/
if(sic_flag)
 {
  /* initialize some arrays */
  memset( TT, 0, sizeof(double)*DIM*DIM);
  if (af->Assemble_Jacobian) 
   {
    memset( dTT_dx,         0, sizeof(double)*DIM*DIM*DIM*MDE);
    memset( dTT_dp,         0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_drs,        0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dc,         0, sizeof(double)*DIM*DIM*MAX_CONC*MDE);
    memset( dTT_dp_liq,     0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dp_gas,     0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dporosity,  0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dsink_mass, 0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dT,         0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dmax_strain,0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dcur_strain,0, sizeof(double)*DIM*DIM*MDE);
   }
                                                                                       
if(pd->e[pg->imtrx][R_MESH1] && cr->MeshMotion != ARBITRARY)
  {
  err = belly_flop(elc->lame_mu);
  EH(err, "error in belly flop");

  err = mesh_stress_tensor(TT, dTT_dx, dTT_dp, dTT_dc, dTT_dp_liq,
                           dTT_dp_gas, dTT_dporosity, dTT_dsink_mass, dTT_dT, dTT_dmax_strain, dTT_dcur_strain,
                           elc->lame_mu, elc->lame_lambda, delta_t, ei[pg->imtrx]->ielem,
			   ip, ip_total);
  if ( cr->MeshMotion == LAGRANGIAN ||
       cr->MeshMotion == DYNAMIC_LAGRANGIAN)
    {
      err = get_convection_velocity(vconv, vconv_old, d_vconv, delta_t, theta);
    }
  else /* No inertia in an Arbitrary Mesh */
    {
      memset( vconv, 0, sizeof(double)*MAX_PDIM );
      if (pd->v[pg->imtrx][MESH_DISPLACEMENT1])
        memset(d_vconv->X, 0, DIM*DIM*MDE*sizeof(dbl));
      if (pd->v[pg->imtrx][VELOCITY1] || pd->v[pg->imtrx][POR_LIQ_PRES])
        memset(d_vconv->v, 0, DIM*DIM*MDE*sizeof(dbl));
      if (pd->v[pg->imtrx][MASS_FRACTION] || pd->v[pg->imtrx][POR_LIQ_PRES])
        memset(d_vconv->C, 0, DIM*MAX_CONC*MDE*sizeof(dbl));
      if (pd->v[pg->imtrx][TEMPERATURE])
        memset(d_vconv->T, 0, DIM*MDE*sizeof(dbl) );
    }
/* For LINEAR ELASTICITY */
      if (cr->MeshFluxModel == LINEAR)
         {
          if (dim == 2)
		{
                TT[2][2] = 1.;
                TT[1][2] = 0.;
                TT[0][2] = 0.;
                }
         }
/*  For Hookian Elasticity and shrinkage */
       else
         {
          if (dim == 2)
		{
                 elast_modulus = elc->lame_mu;
                 if (cr->MeshMotion == ARBITRARY)
                    {
                     TT[2][2] = (1. - fv->volume_change) * elast_modulus;
                    }
                 else
                    {
                     if (cr->MeshFluxModel == NONLINEAR ||
                         cr->MeshFluxModel == HOOKEAN_PSTRAIN ||
                         cr->MeshFluxModel == INCOMP_PSTRAIN )
        			TT[2][2] = (1. - pow(fv->volume_change,2./3.))
					 * elast_modulus - fv->P;
/*              if (cr->MeshFluxModel == INCOMP_PSTRESS) */
                     else  TT[2][2] = 0.;
                    }
                TT[1][2] = 0.;
                TT[0][2] = 0.;
                }
         }
  } /* end of STRESS_TENSOR */
/* calculate real-solid stress here !!*/
if(pd->e[pg->imtrx][R_SOLID1] && cr->MeshMotion != ARBITRARY)
  {
   eqn = R_SOLID1;
   err = belly_flop_rs(elc_rs->lame_mu);
   EH(err, "error in belly flop");

   err = solid_stress_tensor(TT, dTT_dx, dTT_drs, dTT_dp,
                             dTT_dc,  dTT_dp_liq, dTT_dp_gas, dTT_dporosity,
                             dTT_dT, dTT_dmax_strain, elc_rs->lame_mu, elc_rs->lame_lambda);
  if (  cr->MeshMotion == LAGRANGIAN ||
        cr->MeshMotion == DYNAMIC_LAGRANGIAN ||
        cr->MeshMotion == TOTAL_ALE)
    {
      err = get_convection_velocity_rs(vconv, vconv_old, d_vconv, delta_t, theta);
      if ( pd->TimeIntegration != STEADY &&
           pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)] )
        {
          for ( a=0; a<dim; a++)
            {
              vconv[a] = fv_dot->d_rs[a];
              for(b=0; b < VIM; b++)
                {
                  vconv[a] -= fv_dot->d_rs[b]*fv->grad_d_rs[b][a];
                }
            }
          for (a = 0; a < VIM; a++)
            {
              for(b=0; b < VIM; b++)
                {
                  var = SOLID_DISPLACEMENT1 + b;
                  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];
                      d_vconv->rs[a][b][j] += phi_j * (1.+2.*theta)*delta_t*delta(a,b) ;
                      d_vconv->rs[a][b][j] -= (1.+2.*theta)*phi_j*delta_t*fv->grad_d_rs[b][a] ;
                      for(c=0; c < VIM; c++)
                        {
                          /*grad_phi_e[i][j][k][l] = e_k e_l : grad(phi_i e_j )*/
                          d_vconv->rs[a][b][j] -= fv_dot->d_rs[c]* bf[var]->grad_phi_e[j][b][c][a];
                        }
                    }
                }
            }
        }
    }
  else /* No inertia in an Arbitrary Mesh */
    {
      for ( a=0; a<dim; a++)
        {
          vconv[a] = 0.;
        }
    }

   if (dim == 2)
      {
       elast_modulus = elc_rs->lame_mu;
       if (cr->RealSolidFluxModel == NONLINEAR ||
           cr->RealSolidFluxModel == HOOKEAN_PSTRAIN ||
           cr->RealSolidFluxModel == INCOMP_PSTRAIN )
                TT[2][2] = (1.-pow(fv->volume_change,2./3.))*elast_modulus-fv->P;
/*              if (cr->MeshFluxModel == INCOMP_PSTRESS) */
       else  TT[2][2] = 0.;
       TT[1][2] = 0.;
       TT[0][2] = 0.;
      }
  } /* end of REAL_STRESS_TENSOR */
  /* Calculate the residual contribution  */
  
  for ( a=0; a<dim; a++)
     {
        for ( b=0; b<dim; b++)
      	  {
	   func[b] -= fv->snormal[a]*TT[a][b];
	  }
     }
  /*  Jacobian contributions		*/
  if (af->Assemble_Jacobian) 
	{
	for(jvar=0; jvar<dim; jvar++)
	  {
	    var = MESH_DISPLACEMENT1 + jvar;
	    if (pd->v[pg->imtrx][var])
	      {
		for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++)
		  {
		    for(a=0; a<dim; a++)
	  		{
		    	for(b=0; b<dim; b++)
	  		    {
		    		d_func[b][var][j_id] -= 
				   (fv->snormal[a]*dTT_dx[a][b][jvar][j_id] + 
				   fv->dsnormal_dx[a][jvar][j_id]*TT[a][b]);
			    }
			}
		  }
	      }
	  }
  	}
 }	/*  end of if sic_flag	*/

}
/* end of routine force_n_dot_f_bc */

/*
 * This function specifies the external spring resisting force 
 * imposed on a boundary determined by the displacement of a 
 * specified single node nodeset which is balanced by the normal 
 * elastic force within the material this force is applied weakly 
 * to the boundary term in the integral (by parts) from the 
 * Galerkin form of the Lagrangian Mesh Momentum equation
 *
 */
#ifdef NOT_USED
void 
force_n_dot_f_spring_bc(double func[DIM],
		 double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
		 double forcex, 
		 double forcey,
		 double forcez,
		 const int NS1,
		 double x[])
{
    int kdir, dim;
    double force[3];

    int ns_id1=NS1;               /* nodeset ID to track repulsive plane */
    int node_1;                   /* node number of nodeset */
    double initial_xpos1;
    double initial_ypos1;
    double initial_zpos1;
    int idx1;
    int idy1;
    int idz1;
    double delta_xpos1;
    double delta_ypos1;
    double delta_zpos1;
    double actual_xpos1;
    double actual_ypos1;
    double actual_zpos1;
  
    node_1 = psid2nn(ns_id1);                       
    EH(node_1, "Could not find the nsid needed for calculation");
    initial_xpos1 = Coor[0][node_1];
    initial_ypos1 = Coor[1][node_1];
    initial_zpos1 = Coor[2][node_1];
    idx1 = Index_Solution(node_1, MESH_DISPLACEMENT1, 0, 0, -2);
    idy1 = Index_Solution(node_1, MESH_DISPLACEMENT2, 0, 0, -2);
    idz1 = Index_Solution(node_1, MESH_DISPLACEMENT3, 0, 0, -2);
    delta_xpos1 = x[idx1];
    delta_ypos1 = x[idy1];
    delta_zpos1 = x[idz1];
    actual_xpos1 = (initial_xpos1 + delta_xpos1);
    actual_ypos1 = (initial_ypos1 + delta_ypos1);
    actual_zpos1 = (initial_zpos1 + delta_zpos1);
  
    if(af->Assemble_LSA_Mass_Matrix)
      return;

 /* load vector of the desired force */
    force[0] = forcex*delta_xpos1;
    force[1] = forcey*delta_ypos1;
    force[2] = forcez*delta_zpos1;

/*  initialize variables */
    dim = pd->Num_Dim;

/* Calculate the residual contribution  */
 
    for(kdir=0; kdir<dim; kdir++)
    {
          func[kdir] += force[kdir];
    }

}
/* end of routine force_n_dot_f_spring_bc */
#endif

/*
 * This function specifies the external force imposed on a boundary
 * which is apportioned with the distance the solid boundary is with
 * a plane.  The force is balanced by the normal elastic force
 *  within the material
 * this force is applied weakly to the boundary term in the integral
 * (by parts) from the Galerkin form of the Lagrangian Mesh Momentum
 * equation
 *
 */

void 
rep_force_n_dot_f_bc(double func[DIM],
		     double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
	const double pr,	/* coefficient for repulsion force to ensure
				   no penetration of the solid boundary by the
				   free surface                              */
	const double ap,	/* a coefficient in plane equation */
	const double bp,	/* b coefficient in plane equation */
	const double cp,	/* c coefficient in plane equation */
	const double dp,	/* d coefficient in plane equation */
	const double repexp,	/* repulsive force exponent */
	const double friction,	/* friction coefficient */
	const int bc_type)       /*  bc id  */
{
    int j, var, dim, a;
    int jvar;			/* Degree of freedom counter */

    double d_dist[DIM];                /* distance derivatives  */  
    double dist=1e12;		/* squared distance from surface to wall     */
    double xpl, ypl, zpl;   /* coordinates of plane  */
    double factor, denom, dist_sign = 0.0;
    double force = 0.0, d_force = 0.0;

    if(af->Assemble_LSA_Mass_Matrix)
      return;

/*  initialize variables */
    dim = pd->Num_Dim;

/* if pr is 0. => we don't want free surface/wall repulsion and just
   ensure that dist and dist2 are nonzero so nothing bad happens */
  if (pr == 0) return;
/* calculate distance from free surface to solid surface for repulsion calculations */
#if 0
/* first find coordinates of point on plane closest to current gauss point */
  if (cp != 0) EH(GOMA_ERROR, "Can't do repulsion in 3D yet");

/* the following setup is for a plane in 2D (e.g. a line) */
  d_dist[0] = 0;
  d_dist[1] = 0;
  d_dist[2] = 0;
  if (ap == 0)
    {
      xpl = fv->x[0];
      ypl = - dp / bp;
      dist = pow(pow(ypl - fv->x[1],2.0),0.5);
      d_dist[1] = -(ypl - fv->x[1])/dist;
    }
  if (bp == 0)
    {
      ypl = fv->x[1];
      xpl = - dp / ap;
      dist = pow(pow(xpl - fv->x[0],2.0),0.5);
      d_dist[0] = -(xpl - fv->x[0])/dist;
    }
  if (ap != 0 && bp != 0)
    {
      factor = bp/ap + ap /bp;
      xpl = (bp/ap * fv->x[0] - fv->x[1] - dp/bp) / factor;
      ypl = (ap/bp * fv->x[1] - fv->x[0] - dp/ap) / factor;
      dist = pow(pow(xpl - fv->x[0],2.0) + pow(ypl - fv->x[1],2.0),0.5);
      d_dist[0] = (-(xpl - fv->x[0])* (1 - bp/ap/factor) 
		   - (ypl - fv->x[1])/factor 
		   )/dist;
      d_dist[1] = (-(ypl - fv->x[1])* (1 - ap/bp/factor) 
		   - (xpl - fv->x[0])/factor 
		   )/dist;
    }
#else
      xpl = fv->x[0];
      ypl = fv->x[1];
      if( dim == 3)
	{ zpl = fv->x[2];}
      else
	{ zpl = 0.0;}
      denom = sqrt(ap*ap+bp*bp+cp*cp);
      factor = ap*xpl+bp*ypl+cp*zpl+dp;
      dist = fabs(factor)/denom;
      d_dist[0] = SGN(factor)*ap/denom;
      d_dist[1] = SGN(factor)*bp/denom;
      d_dist[2] = SGN(factor)*cp/denom;
      switch ( bc_type )
	{
	case REP_FORCE_BC:
	case REP_FORCE_RS_BC:
           force = -pr/pow(dist, repexp); 
           d_force = pr*repexp/pow(dist, repexp+1);
        break;
        case ATTR_FORCE_BC:
        case ATTR_FORCE_RS_BC:
           if(DOUBLE_NONZERO(dp))
             { dist_sign = -SGN(dp); }
           else if (DOUBLE_NONZERO(cp))
             { dist_sign = SGN(cp); }
           else if (DOUBLE_NONZERO(bp))
             { dist_sign = SGN(bp); }
           dist_sign *= SGN(factor);
           force = pr*dist_sign*pow(dist, repexp); 
           d_force = pr*dist_sign*repexp*pow(dist, repexp-1);
        break;
        default:
           EH(GOMA_ERROR, "Invalid BC id in force_repulsion");
        break;
         }
#endif

  if (af->Assemble_Jacobian)
    {
      /* 
       *  Evaluate sensitivity to displacements d()/dx 
       */
      for (jvar=0; jvar<dim; jvar++)
	{
	  var = MESH_DISPLACEMENT1 + jvar;
	  if (pd->v[pg->imtrx][var]) 
	    {
	      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		{
		  for (a=0; a<ei[pg->imtrx]->ielem_dim; a++)
		    {
			  
		     d_func[a][var][j] += force*fv->dsnormal_dx[a][jvar][j]  
		       + d_force * fv->snormal[a]*d_dist[jvar]*bf[var]->phi[j]; 
                     d_func[a][var][j] += friction*
                            (force*fv->dstangent_dx[0][a][jvar][j]  
		                    + d_force * fv->stangent[0][a]
			     * d_dist[jvar] * bf[var]->phi[j]); 
                      if (dim == 3)
	                {
                     d_func[a][var][j] += friction*
                            (force*fv->dstangent_dx[1][a][jvar][j]  
		                    + d_force * fv->stangent[1][a]
			     * d_dist[jvar] * bf[var]->phi[j]); 
	                }
		    }
		}
	    }
	}
    }

  for (a=0; a<ei[pg->imtrx]->ielem_dim; a++)
    {  
      func[a] += force * fv->snormal[a]; 
      func[a] += force * fv->stangent[0][a]*friction; 
      if (dim == 3)
	{
          func[a] += force * fv->stangent[1][a]*friction; 
	}
    }
    

}
/* end of routine rep_force_n_dot_f_bc */
/********************************************************************************/

/*
 * This function specifies the external force imposed on a boundary
 * which is apportioned with the distance the solid boundary is from
 * a roll.  The force is balanced by the normal elastic force
 *  within the material
 * this force is applied weakly to the boundary term in the integral
 * (by parts) from the Galerkin form of the Lagrangian Mesh Momentum
 * equation
 *
 */

void 
rep_force_roll_n_dot_f_bc(double func[DIM],
		     double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
	const double pr,	/* coefficient for repulsion force to ensure
				   no penetration of the solid boundary by the
				   free surface                              */
	const double origin[3],	/* roll axis origin (x,y,z) */
	const double dir_angle[3],	/* axis direction angles */
	const double roll_rad,	/* roll radius */
	const double repexp,	/* repulsive force exponent */
	const double friction,	/* friction coefficient */
	const int bc_type)       /*  bc id  */
{
    int j, var, dim, a;
    int jvar;			/* Degree of freedom counter */

    double d_dist[DIM];                /* distance derivatives  */  
    double dist=1e12;		/* squared distance from surface to wall     */
    double factor;
    double coord[3], axis_pt[3], t;
    double force = 0.0, d_force = 0.0;

    if(af->Assemble_LSA_Mass_Matrix)
      return;

/*  initialize variables */
    dim = pd->Num_Dim;

/* if pr is 0. => we don't want free surface/wall repulsion and just
   ensure that dist and dist2 are nonzero so nothing bad happens */
  if (pr == 0) return;
/* calculate distance from free surface to solid surface for repulsion calculations */

      coord[0] = fv->x[0];
      coord[1] = fv->x[1];
      if( dim == 3)
	{ coord[2] = fv->x[2];}
      else
	{ coord[2] = 0.0;}

/*  find intersection of axis with normal plane - i.e., locate point on
 *          axis that intersects plane normal to axis that contains local point. */

    factor = SQUARE(dir_angle[0]) + SQUARE(dir_angle[1]) + SQUARE(dir_angle[2]);
    t = (dir_angle[0]*(coord[0]-origin[0]) + dir_angle[1]*(coord[1]-origin[1])
        + dir_angle[2]*(coord[2]-origin[2]))/factor;
    axis_pt[0] = origin[0]+dir_angle[0]*t;
    axis_pt[1] = origin[1]+dir_angle[1]*t;
    axis_pt[2] = origin[2]+dir_angle[2]*t;

/*  compute radius and radial direction */

    dist = sqrt( SQUARE(coord[0]-axis_pt[0]) + SQUARE(coord[1]-axis_pt[1]) +
                SQUARE(coord[2]-axis_pt[2]) ) - roll_rad;
    d_dist[0] = (coord[0]-axis_pt[0])/(dist+roll_rad)*(1.-SQUARE(dir_angle[0])/factor)
          +(coord[1]-axis_pt[1])/(dist+roll_rad)*(-dir_angle[1]*dir_angle[0]/factor)
          +(coord[2]-axis_pt[2])/(dist+roll_rad)*(-dir_angle[2]*dir_angle[0]/factor);
    d_dist[1] = (coord[1]-axis_pt[1])/(dist+roll_rad)*(1.-SQUARE(dir_angle[1])/factor)
          +(coord[0]-axis_pt[0])/(dist+roll_rad)*(-dir_angle[0]*dir_angle[1]/factor)
          +(coord[2]-axis_pt[2])/(dist+roll_rad)*(-dir_angle[2]*dir_angle[1]/factor);
    d_dist[2] = (coord[2]-axis_pt[2])/(dist+roll_rad)*(1.-SQUARE(dir_angle[2])/factor)
          +(coord[0]-axis_pt[0])/(dist+roll_rad)*(-dir_angle[0]*dir_angle[2]/factor)
          +(coord[1]-axis_pt[1])/(dist+roll_rad)*(-dir_angle[1]*dir_angle[2]/factor);
      switch ( bc_type )
	{
	case REP_FORCE_ROLL_BC:
	case REP_FORCE_ROLL_RS_BC:
           force = -pr/pow(dist, repexp); 
           d_force = pr*repexp/pow(dist, repexp+1);
        break;
        default:
           EH(GOMA_ERROR, "Invalid BC id in force_repulsion roll");
        break;
         }

  if (af->Assemble_Jacobian)
    {
      /* 
       *  Evaluate sensitivity to displacements d()/dx 
       */
      for (jvar=0; jvar<dim; jvar++)
	{
	  var = MESH_DISPLACEMENT1 + jvar;
	  if (pd->v[pg->imtrx][var]) 
	    {
	      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		{
		  for (a=0; a<ei[pg->imtrx]->ielem_dim; a++)
		    {
			  
		     d_func[a][var][j] += force*fv->dsnormal_dx[a][jvar][j]  
		       + d_force * fv->snormal[a]*d_dist[jvar]*bf[var]->phi[j]; 
                     d_func[a][var][j] += friction*
                            (force*fv->dstangent_dx[0][a][jvar][j]  
		                    + d_force * fv->stangent[0][a]
			     * d_dist[jvar] * bf[var]->phi[j]); 
                      if (dim == 3)
	                {
                     d_func[a][var][j] += friction*
                            (force*fv->dstangent_dx[1][a][jvar][j]  
		                    + d_force * fv->stangent[1][a]
			     * d_dist[jvar] * bf[var]->phi[j]); 
	                }
		    }
		}
	    }
	}
    }

  for (a=0; a<ei[pg->imtrx]->ielem_dim; a++)
    {  
      func[a] += force * fv->snormal[a]; 
      func[a] += force * fv->stangent[0][a]*friction; 
      if (dim == 3)
	{
          func[a] += force * fv->stangent[1][a]*friction; 
	}
    }
    

}
/* end of routine rep_force_roll_n_dot_f_bc */
/********************************************************************************/
/*
 * This function specifies the external force imposed on a boundary
 * which is balanced by the normal elastic force within the material
 * this force is applied weakly to the boundary term in the integral
 * (by parts) from the Galerkin form of the Lagrangian Mesh Momentum
 * equation
 *
 */

void
norm_force_n_dot_f_bc(double func[DIM],
		      double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
		      double forcex,
		      double forcey,
		      double forcez)
{
  int j_id, jvar, kdir, dim, var;
  double force[3];
  
  if(af->Assemble_LSA_Mass_Matrix)
    return;

  /* load vector of the desired force */
  force[0] = forcex; /* normal for applied to surface */
  force[1] = forcey; /* tangential force applied to surface */
  force[2] = forcez;
  
  /*  initialize variables */
  dim = pd->Num_Dim;
  
  if (af->Assemble_Jacobian) {
    for(kdir=0; kdir<dim; kdir++)
      {
	for(jvar=0; jvar<dim; jvar++)
	  {
	    var = MESH_DISPLACEMENT1 + jvar;
	    if (pd->v[pg->imtrx][var])
	      {
		for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++)
		  {
		    
		    d_func[kdir][var][j_id] += (
						fv->dsnormal_dx[kdir][jvar][j_id]*force[0] +
						fv->dstangent_dx[0][kdir][jvar][j_id]*force[1]);
		    if (dim == 3)
		      {
			d_func[kdir][var][j_id] += fv->dstangent_dx[1][kdir][jvar][j_id]*force[2];
		      }
		  }
	      }
	  }
      }
  }
  
  /* Calculate the residual contribution  */
  
  for(kdir=0; kdir<dim; kdir++)
    {
      /* add on normal and tangential components of surface force */
      func[kdir] += fv->snormal[kdir]*force[0];
      func[kdir] += fv->stangent[0][kdir]*force[1];
      if (dim == 3)
	{
	  func[kdir] += fv->stangent[1][kdir]*force[2];
	}
    }
  
}
/* end of routine norm_force_n_dot_f_bc */
/********************************************************************************/

/*
 * This function specifies the external force imposed on a boundary
 * which is balanced by the normal elastic force within the material
 * this force is applied weakly to the boundary term in the integral
 * (by parts) from the Galerkin form of the Lagrangian Mesh Momentum
 * equation
 *
 */

void
friction_n_dot_f_bc(double func[DIM],
		      double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
		      double frict_coeff,
		      const int blk_id,
		      const double delta_t,
		      const double theta,
		      const int ip,
		      const int ip_total,
		      const int bc_type,
		      const double time_value,
		      const double u_par[],
		      const int n_par)
{
  int j_id, jvar, dim, var;
  int err, a, b, c, j;
  int eqn;
  double normal_force, velo_mag, tang_velo[2], tang_force;
  double d_nforce, d_tforce, phi_j;
  
  double TT[MAX_PDIM][MAX_PDIM];   /**  solid stresses  **/
  double dTT_drs[DIM][DIM][DIM][MDE];
  double dTT_dx[MAX_PDIM][MAX_PDIM][MAX_PDIM][MDE];
  double dTT_dp[MAX_PDIM][MAX_PDIM][MDE];
  double dTT_dc[MAX_PDIM][MAX_PDIM][MAX_CONC][MDE];
  double dTT_dp_liq[DIM][DIM][MDE];/* Sensitivity of stress tensor...
                                    to nodal porous liquid pressure*/
  double dTT_dp_gas[DIM][DIM][MDE];/* Sensitivity of stress tensor...
                                    to nodal porous gas pressure*/
  double dTT_dporosity[DIM][DIM][MDE];/* Sensitivity of stress tensor...
                                    to nodal porosity*/
  double dTT_dsink_mass[DIM][DIM][MDE];/* Sensitivity of stress tensor...
                                    to sink_mass*/
  double dTT_dT[MAX_PDIM][MAX_PDIM][MDE];
  double dTT_dmax_strain[MAX_PDIM][MAX_PDIM][MDE];
  double dTT_dcur_strain[MAX_PDIM][MAX_PDIM][MDE];
  double elast_modulus;

  double vconv[MAX_PDIM], vconv_old[MAX_PDIM];/*Calculated convection velocity*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;

  double velo_acous, frict_acous;    
  double factor, dnum, dden;
  const double *pnum, *pden;
  int n_pade;


  if(af->Assemble_LSA_Mass_Matrix)
    return;

  if(blk_id != -1 && blk_id != ei[pg->imtrx]->elem_blk_id)
	return;

  /*  initialize variables */
  dim = pd->Num_Dim;
  
  /* initialize some arrays */
  memset( TT, 0, sizeof(double)*DIM*DIM);
  if (af->Assemble_Jacobian) {
    memset( dTT_dx,         0, sizeof(double)*DIM*DIM*DIM*MDE);
    memset( dTT_dp,         0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_drs,        0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dc,         0, sizeof(double)*DIM*DIM*MAX_CONC*MDE);
    memset( dTT_dp_liq,     0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dp_gas,     0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dporosity,  0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dsink_mass, 0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dT,         0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dmax_strain,0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dcur_strain,0, sizeof(double)*DIM*DIM*MDE);
  }
                                                                                       
if(pd->e[pg->imtrx][R_MESH1] && cr->MeshMotion != ARBITRARY)
  {
  err = belly_flop(elc->lame_mu);
  EH(err, "error in belly flop");

  err = mesh_stress_tensor(TT, dTT_dx, dTT_dp, dTT_dc, dTT_dp_liq,
                           dTT_dp_gas, dTT_dporosity, dTT_dsink_mass, dTT_dT, dTT_dmax_strain, dTT_dcur_strain,
                           elc->lame_mu, elc->lame_lambda, delta_t, ei[pg->imtrx]->ielem,
			   ip, ip_total);
  if ( cr->MeshMotion == LAGRANGIAN ||
       cr->MeshMotion == DYNAMIC_LAGRANGIAN)
    {
      err = get_convection_velocity(vconv, vconv_old, d_vconv, delta_t, theta);
    }
  else /* No inertia in an Arbitrary Mesh */
    {
      memset( vconv, 0, sizeof(double)*MAX_PDIM );
      if (pd->v[pg->imtrx][MESH_DISPLACEMENT1])
        memset(d_vconv->X, 0, DIM*DIM*MDE*sizeof(dbl));
      if (pd->v[pg->imtrx][VELOCITY1] || pd->v[pg->imtrx][POR_LIQ_PRES])
        memset(d_vconv->v, 0, DIM*DIM*MDE*sizeof(dbl));
      if (pd->v[pg->imtrx][MASS_FRACTION] || pd->v[pg->imtrx][POR_LIQ_PRES])
        memset(d_vconv->C, 0, DIM*MAX_CONC*MDE*sizeof(dbl));
      if (pd->v[pg->imtrx][TEMPERATURE])
        memset(d_vconv->T, 0, DIM*MDE*sizeof(dbl) );
    }
/* For LINEAR ELASTICITY */
      if (cr->MeshFluxModel == LINEAR)
         {
          if (dim == 2)
		{
                TT[2][2] = 1.;
                TT[1][2] = 0.;
                TT[0][2] = 0.;
                }
         }
/*  For Hookian Elasticity and shrinkage */
       else
         {
          if (dim == 2)
		{
                 elast_modulus = elc->lame_mu;
                 if (cr->MeshMotion == ARBITRARY)
                    {
                     TT[2][2] = (1. - fv->volume_change) * elast_modulus;
                    }
                 else
                    {
                     if (cr->MeshFluxModel == NONLINEAR ||
                         cr->MeshFluxModel == HOOKEAN_PSTRAIN ||
                         cr->MeshFluxModel == INCOMP_PSTRAIN )
        			TT[2][2] = (1. - pow(fv->volume_change,2./3.))
					 * elast_modulus - fv->P;
/*              if (cr->MeshFluxModel == INCOMP_PSTRESS) */
                     else  TT[2][2] = 0.;
                    }
                TT[1][2] = 0.;
                TT[0][2] = 0.;
                }
         }
  } /* end of STRESS_TENSOR */
/* calculate real-solid stress here !!*/
if(pd->e[pg->imtrx][R_SOLID1] && cr->MeshMotion != ARBITRARY)
  {
   eqn = R_SOLID1;
   err = belly_flop_rs(elc_rs->lame_mu);
   EH(err, "error in belly flop");

   err = solid_stress_tensor(TT, dTT_dx, dTT_drs, dTT_dp,
                             dTT_dc,  dTT_dp_liq, dTT_dp_gas, dTT_dporosity,
                             dTT_dT, dTT_dmax_strain, elc_rs->lame_mu, elc_rs->lame_lambda);
  if (  cr->MeshMotion == LAGRANGIAN ||
        cr->MeshMotion == DYNAMIC_LAGRANGIAN ||
        cr->MeshMotion == TOTAL_ALE)
    {
      err = get_convection_velocity_rs(vconv, vconv_old, d_vconv, delta_t, theta);
      if ( pd->TimeIntegration != STEADY &&
           pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)] )
        {
          for ( a=0; a<dim; a++)
            {
              vconv[a] = fv_dot->d_rs[a];
              for(b=0; b < VIM; b++)
                {
                  vconv[a] -= fv_dot->d_rs[b]*fv->grad_d_rs[b][a];
                }
            }
          for (a = 0; a < VIM; a++)
            {
              for(b=0; b < VIM; b++)
                {
                  var = SOLID_DISPLACEMENT1 + b;
                  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];
                      d_vconv->rs[a][b][j] += phi_j * (1.+2.*theta)*delta_t*delta(a,b) ;
                      d_vconv->rs[a][b][j] -= (1.+2.*theta)*phi_j*delta_t*fv->grad_d_rs[b][a] ;
                      for(c=0; c < VIM; c++)
                        {
                          /*grad_phi_e[i][j][k][l] = e_k e_l : grad(phi_i e_j )*/
                          d_vconv->rs[a][b][j] -= fv_dot->d_rs[c]* bf[var]->grad_phi_e[j][b][c][a];
                        }
                    }
                }
            }
        }
    }
  else /* No inertia in an Arbitrary Mesh */
    {
      for ( a=0; a<dim; a++)
        {
          vconv[a] = 0.;
        }
    }

   if (dim == 2)
      {
       elast_modulus = elc_rs->lame_mu;
       if (cr->RealSolidFluxModel == NONLINEAR ||
           cr->RealSolidFluxModel == HOOKEAN_PSTRAIN ||
           cr->RealSolidFluxModel == INCOMP_PSTRAIN )
                TT[2][2] = (1.-pow(fv->volume_change,2./3.))*elast_modulus-fv->P;
/*              if (cr->MeshFluxModel == INCOMP_PSTRESS) */
       else  TT[2][2] = 0.;
       TT[1][2] = 0.;
       TT[0][2] = 0.;
      }
  } /* end of REAL_STRESS_TENSOR */
                                                                                       
/* compute normal force and velocity magnitude	*/
  normal_force = 0;
  tang_force = 0;
  for ( a=0; a<dim; a++)
     {
        for ( b=0; b<dim; b++)
      	  {
	   normal_force += fv->snormal[a]*TT[a][b]*fv->snormal[b];
	   tang_force += fv->snormal[a]*TT[a][b]*fv->stangent[0][b];
	  }
     }
  velo_mag = 0;
  tang_velo[0] = 0;
  tang_velo[1] = 0;
  for ( a=0; a<dim; a++)
     {
	velo_mag += vconv[a]*vconv[a];
	tang_velo[0] += fv->stangent[0][a]*vconv[a];
	tang_velo[1] += fv->stangent[1][a]*vconv[a];
     }
  velo_mag = sqrt(velo_mag);
  if( DOUBLE_ZERO(velo_mag))
       { EH(GOMA_ERROR,"Trouble with sliding friction bc - zero relative velocity.");}

  frict_acous = 1.0;
  if( bc_type == FRICTION_ACOUSTIC_BC)
	{
	 factor = u_par[0];
	 n_pade = n_par/2;
	 pnum = &u_par[1];
	 pden = &u_par[n_pade+1];
	 dnum = pnum[n_pade-1];
	 dden = pden[n_pade-2];
	 for(a=n_pade-2 ; a>=0 ; a--)
		{dnum = pnum[a]+fv->x[0]*dnum;  }
	 for (a=n_pade-3 ; a>=0 ; a--)
		{dden = pden[a] + fv->x[0]*dden;  }
	 velo_acous = factor*dnum/dden;
	 if( velo_acous > velo_mag)
		{
		frict_acous = 2.0/M_PIE*asin(velo_mag/velo_acous);
		}
	}

  if (af->Assemble_Jacobian) {
	for(jvar=0; jvar<dim; jvar++)
	  {
	    var = MESH_DISPLACEMENT1 + jvar;
	    if (pd->v[pg->imtrx][var])
	      {
		for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++)
		  {
		    
		    d_nforce = 0;
		    d_tforce = 0;
		    for(a=0; a<dim; a++)
	  		{
		    	for(b=0; b<dim; b++)
	  		    {
	   			d_nforce += 
			fv->dsnormal_dx[a][jvar][j_id]*TT[a][b]*fv->snormal[b] +
			fv->snormal[a]*dTT_dx[a][b][jvar][j_id]*fv->snormal[b] +
			fv->snormal[a]*TT[a][b]*fv->dsnormal_dx[b][jvar][j_id];
	   			d_tforce += 
			fv->dsnormal_dx[a][jvar][j_id]*TT[a][b]*fv->stangent[0][b] +
			fv->snormal[a]*dTT_dx[a][b][jvar][j_id]*fv->stangent[0][b] +
			fv->snormal[a]*TT[a][b]*fv->dstangent_dx[0][b][jvar][j_id];
			    }
			}

		    d_func[0][var][j_id] += (d_tforce - frict_coeff*frict_acous
				*d_nforce*tang_velo[0]/velo_mag);
		  }
	      }
	  }
  }
  
  /* Calculate the residual contribution  */
  
  *func += tang_force - frict_coeff*frict_acous*normal_force
		*tang_velo[0]/velo_mag;
  
}
/* end of routine friction_n_dot_f_bc */
/********************************************************************************/

/*
 * This function equates the normal (natural) force on the boundary
 * in the solid phase to the normal (natural) force on the boundary
 * in the liquid phase by adding the local contribution of the
 * volume integrals in the liquid into the total contribution in
 * the solid, and decrementing the arbitrary mesh contribution
 * when in the fluid phase.
 *
 * MMH: Since this constructs Jacobian entries based on
 * already-existing Jacobian entries, I believe the LSA details have
 * already been taken care of by the time this is called.
 */

void 
put_liquid_stress_in_solid(int id, /* local element node number for the 
				    * current node whose residual contribution
				    * is being sought                        */
			   int I, /* Global node number                      */
			   int ielem_dim, /* physical dimension of the elem  */
			   double resid_vector[], /* Residual vector         */
			   int i_mat_solid, /* elem block id's of solid       */
			   int i_mat_fluid, /* elem block id's of liquid      */
			   int local_node_list_fs[], /* MDE list to keep track
						      * of nodes at which 
						      * solid contributions 
						      * have been transfered
						      * to liquid (fluid-solid
						      * boundaries)          */
			   double scale) /* Scale factor, nondimension       */
{
    int j_id, dim, wim, var, pvar,  p, q, p2, w, id_dofmom, id_dofmesh, offset, mode;
    int peqn_mom, peqn_solid;
    int ieqn_mom, ieqn_solid;
    int xfixed[DIM];

    int v_s[MAX_MODES][DIM][DIM]; 
    int v_g[DIM][DIM]; 

    NODE_INFO_STRUCT *node = Nodes[I];
    NODAL_VARS_STRUCT *nv = node->Nodal_Vars_Info[pg->imtrx];

    dim = pd->Num_Dim;
    wim   = dim;
    if(pd->CoordinateSystem == SWIRLING) wim = wim+1;

    if (node->DBC[pg->imtrx]) {
     
      offset = get_nodal_unknown_offset(nv, R_MESH1, -2, 0, NULL);
      xfixed[0] = (offset >= 0) && (node->DBC[pg->imtrx][offset] != -1);
      offset = get_nodal_unknown_offset(nv, R_MESH2, -2, 0, NULL);
      xfixed[1] = (offset >= 0) && (node->DBC[pg->imtrx][offset] != -1);	
      offset = get_nodal_unknown_offset(nv, R_MESH3, -2, 0, NULL);
      xfixed[2] = (offset >= 0) && (node->DBC[pg->imtrx][offset] != -1);
    } else {
      xfixed[0] = xfixed[1] = xfixed[2] = 0;
    }


    /*
     * if you are in the solid phase, return without doing anything
     * In the solid phase, there are no fluid momentum equations.
     */
    if (!pd->e[pg->imtrx][R_MOMENTUM1]) return;

    
    id_dofmom = ei[pg->imtrx]->ln_to_dof[R_MOMENTUM1][id];
    id_dofmesh = ei[pg->imtrx]->ln_to_dof[R_MESH1][id];

    if (Current_EB_ptr->Elem_Blk_Id != i_mat_fluid)
      { 
	if (Current_EB_ptr->Elem_Blk_Id != i_mat_solid) {
	  EH(GOMA_ERROR, "Improper solid-fluid block ids");
	}
	/* note, this may not account for all situations - we may not want
	   to quit here, but we'll do it for now */
	return;
      }

    /*
     * if this nodal contribution has already been added to fluid momentum
     * equation (i.e. we are at a corner on the second side) return
     * without doing anything
     */
    if (local_node_list_fs[id] == -1) {
      local_node_list_fs[id] = 1;
    } else {
      return;
    }

    /*
     * check to make sure that both mesh dof and velocity
     * dof exist at this node
     */
    if (Dolphin[pg->imtrx][I][R_MESH1] <= 0)     return;
    if (Dolphin[pg->imtrx][I][R_MOMENTUM1] <= 0) return;

    /*
     * loop over directions and add local contribution to liquid
     * momentum into local contribution for solid momentum while
     * decrementing the psuedo-solid momentum
     */
    if (af->Assemble_Residual) {
      for (p = 0; p < wim; p++) {
	/* don't do this if the displacement is fixed for this node */
	if (! xfixed[p]) {
	  ieqn_mom = R_MOMENTUM1 + p;
	  ieqn_solid = R_MESH1 + p;
	  id_dofmom = ei[pg->imtrx]->ln_to_dof[ieqn_mom][id];
	  id_dofmesh = ei[pg->imtrx]->ln_to_dof[ieqn_solid][id];
	  lec->R[upd->ep[pg->imtrx][ieqn_solid]][id_dofmesh] =
	      scale * lec->R[upd->ep[pg->imtrx][ieqn_mom]][id_dofmom];
	}
      }
    }
    
    /*
     * loop over directions and add local contribution to solid
     * momentum into local contribution for liquid momentum
     */
    if (af->Assemble_Jacobian)
      {
	for (p = 0; p < wim; p++) {
	    /* 
	     * don't do this if the mesh is fixed for this node 
	     *  HKM -> I assume this is the mesh since xfixed[] is used
	     *         in the statement below, even though the original
	     *         comment stated something about velocity being fixed.
	     */
	    if (! xfixed[p])
	      {
		ieqn_mom = R_MOMENTUM1 + p;
		peqn_mom = upd->ep[pg->imtrx][ieqn_mom];
		ieqn_solid = R_MESH1 + p;
		peqn_solid = upd->ep[pg->imtrx][ieqn_solid];
		id_dofmom = ei[pg->imtrx]->ln_to_dof[ieqn_mom][id];
		id_dofmesh = ei[pg->imtrx]->ln_to_dof[ieqn_solid][id];

		/* Add contributions due to all nodal sensitivities in solid element */

		/*
		 * local J_m_d -> J_d_d
		 */
		for ( q=0; q<dim; q++)
		  {
		    var = MESH_DISPLACEMENT1+q;
		    if ( pd->v[pg->imtrx][var] )
		      {
			pvar = upd->vp[pg->imtrx][var];
			for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
			  {
			    lec->J[peqn_solid][pvar][id_dofmesh][j_id] =
				scale*lec->J[peqn_mom][pvar][id_dofmom][j_id];
			  }
		      }
		  }

		/*
		 * local J_m_P -> J_d_P
		 */
		var = PRESSURE;
		if ( pd->v[pg->imtrx][var] )
		  {
		    pvar = upd->vp[pg->imtrx][var];
		    for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
		      {			
			lec->J[peqn_solid][pvar][id_dofmesh][j_id] =
			    scale*lec->J[peqn_mom][pvar][id_dofmom][j_id];
		      }
		  }
		

		/*
		 * local J_m_T -> J_d_T
		 */
		var = TEMPERATURE;
		if ( pd->v[pg->imtrx][var] )
		  {
		    pvar = upd->vp[pg->imtrx][var];
		    for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
		      {
			lec->J[peqn_solid][pvar][id_dofmesh][j_id] =
			    scale*lec->J[peqn_mom][pvar][id_dofmom][j_id];
		      }
		  }
		
		/*
		 * local J_m_c -> J_d_c
		 */
		var = MASS_FRACTION;
		if ( pd->v[pg->imtrx][var] )
		  {
		    for ( w=0; w<pd->Num_Species_Eqn; w++)
		      {
			pvar = MAX_PROB_VAR + w;
			for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
			  {
			    if (fabs(lec->J[peqn_solid][pvar][id_dofmesh][j_id]) > 1.e-8)
			      {
				lec->J[peqn_solid][pvar][id_dofmesh][j_id] =
				    scale*lec->J[peqn_mom][pvar][id_dofmom][j_id];
			      }
			  }
		      }
		  }
		
		/*
		 * local J_m_p_liq -> J_d_p_liq
		 */
                    var = POR_LIQ_PRES;
                    pvar = upd->vp[pg->imtrx][var];
		    if ( pd->v[pg->imtrx][var] )
		      {
			for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
			  {
			    if (fabs(lec->J[peqn_solid][pvar][id_dofmesh][j_id]) > 1.e-8)
			      {
				lec->J[peqn_solid][pvar][id_dofmesh][j_id] =
				    scale*lec->J[peqn_mom][pvar][id_dofmom][j_id];
			      }
			  }
		      }
		/*
		 * local J_m_p_gas -> J_d_p_gas
		 */
                    var = POR_GAS_PRES;
                    pvar = upd->vp[pg->imtrx][var];
		    if ( pd->v[pg->imtrx][var] )
		      {
			for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
			  {
			    if (fabs(lec->J[peqn_solid][pvar][id_dofmesh][j_id]) > 1.e-8)
			      {
				lec->J[peqn_solid][pvar][id_dofmesh][j_id] =
				    scale*lec->J[peqn_mom][pvar][id_dofmom][j_id];
			      }
			  }
		      }
		/*
		 * local J_m_porosity -> J_d_porosity
		 */
                    var = POR_POROSITY;
                    pvar = upd->vp[pg->imtrx][var];
		    if ( pd->v[pg->imtrx][var] )
		      {
			for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
			  {
			    if (fabs(lec->J[peqn_solid][pvar][id_dofmesh][j_id]) > 1.e-8)
			      {
				lec->J[peqn_solid][pvar][id_dofmesh][j_id] =
				    scale*lec->J[peqn_mom][pvar][id_dofmom][j_id];
			      }
			  }
		      }
		
		/*
		 * local J_m_v -> J_d_v  
		 */
		for ( q=0; q<wim; q++)
		  {
		    var = VELOCITY1+q;
		    if ( pd->v[pg->imtrx][var] )
		      {
			pvar = upd->vp[pg->imtrx][var];
			for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
			  {
			    lec->J[peqn_solid][pvar][id_dofmesh][j_id] =
				scale*lec->J[peqn_mom][pvar][id_dofmom][j_id];
			  }
		      }
		  }

		/*
		 * local J_m_S -> J_d_S 
		 */
		var = POLYMER_STRESS11;
		if ( pd->v[pg->imtrx][var] )
		  {
		    (void) stress_eqn_pointer(v_s);
		    for ( mode=0; mode<vn->modes; mode++)
		      {
			for ( p2=0; p2<VIM; p2++)
			  {
			    for ( q=0; q<VIM; q++)
			      {
				var =  v_s[mode][p2][q];
				if ( pd->v[pg->imtrx][var] )
				  {
				    pvar = upd->vp[pg->imtrx][var];
				    for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
				      {
					lec->J[peqn_solid][pvar][id_dofmesh][j_id] =
					  scale*lec->J[peqn_mom][pvar][id_dofmom][j_id];
				      }
				  }
			      }
			  }
		      }
		  }

		/*
		 * local J_m_G -> J_d_G
		 */
		var = VELOCITY_GRADIENT11;
		if ( pd->v[pg->imtrx][var] )
		  {
		    v_g[0][0] = VELOCITY_GRADIENT11;
		    v_g[0][1] = VELOCITY_GRADIENT12;
		    v_g[1][0] = VELOCITY_GRADIENT21;
		    v_g[1][1] = VELOCITY_GRADIENT22;
		    v_g[0][2] = VELOCITY_GRADIENT13;
		    v_g[1][2] = VELOCITY_GRADIENT23;
		    v_g[2][0] = VELOCITY_GRADIENT31;
		    v_g[2][1] = VELOCITY_GRADIENT32; 
		    v_g[2][2] = VELOCITY_GRADIENT33; 
		    for ( p2=0; p2<VIM; p2++)
		      {
			for ( q=0; q<VIM; q++)
			  {
			    var =  v_g[p2][q];
			    if ( pd->v[pg->imtrx][var] )
			      {
				pvar = upd->vp[pg->imtrx][var];
				for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
				  {
				    lec->J[peqn_solid][pvar][id_dofmesh][j_id] =
				      scale*lec->J[peqn_mom][pvar][id_dofmom][j_id];
				    
				  }
			      }
			  }
		      }
		  }

	      } /* end of check for velocity dirichlet condition */
	}
      } /* end of Jacobian entries */
    
} /* end of routine put_liquid_stress_in_solid */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/*
 * This function equates the normal (natural) force on the boundary in the GENERALIZED ALE
 * solid phase to the normal (natural) force on the boundary in the
 * liquid phase by adding the local
 * contribution of the volume integrals in the liquid into the total contribution in the 
 * REAL solid, and decrementing the arbitrary mesh contribution when in the fluid phase
 * This routine DIFFERS from "put_liquid_stressin_solid" by the equation to which it is
 * applied, i.e., that one is lagrangian, this one applies it to the real ALE solid. 
 *
 * MMH: Since this constructs Jacobian entries based on
 * already-existing Jacobian entries, I believe the LSA details have
 * already been taken care of by the time this is called.
 */

void 
put_liquid_stress_in_solid_ALE(int id, /* local element node number for the 
					* current node whose residual 
					* contribution is being sought       */
			   int I, /* Global node number                      */
			   int ielem_dim, /* physical dimension of the elem  */
			   double resid_vector[], /* Residual vector         */
			   int i_mat_solid, /* elem block id's of solid       */
			   int i_mat_fluid, /* elem block id's of liquid      */
			   int local_node_list_fs[], /* MDE list to keep track
						      * of nodes at which 
						      * solid contributions 
						      * have been transfered
						      * to liquid (fluid-solid
						      * boundaries)          */
			   double scale) /* Scale factor, nondimension       */
{
    int j_id, dim, var, p, q, w, id_dofsol, id_dofmom, offset;
    int peqn_mom, peqn_solid, pvar;
    int ieqn_mom, ieqn_solid;
    int xfixed[DIM];
    NODE_INFO_STRUCT *node = Nodes[I];
    NODAL_VARS_STRUCT *nv = node->Nodal_Vars_Info[pg->imtrx];

    dim = pd->Num_Dim;

    if (node->DBC[pg->imtrx]) {
      offset = get_nodal_unknown_offset(nv, R_SOLID1, -2, 0, NULL);
      xfixed[0] = (offset >= 0) && (node->DBC[pg->imtrx][offset] != -1);
      offset = get_nodal_unknown_offset(nv, R_SOLID2, -2, 0, NULL);
      xfixed[1] = (offset >= 0) && (node->DBC[pg->imtrx][offset] != -1);	
      offset = get_nodal_unknown_offset(nv, R_SOLID3, -2, 0, NULL);
      xfixed[2] = (offset >= 0) && (node->DBC[pg->imtrx][offset] != -1);
    } else {
      xfixed[0] = xfixed[1] = xfixed[2] = 0;
    } 

    id_dofmom = ei[pg->imtrx]->ln_to_dof[R_MOMENTUM1][id];
    id_dofsol = ei[pg->imtrx]->ln_to_dof[R_SOLID1][id];
    /* 
     * if you are in the solid phase, return without doing anything.
     * In the solid phase, there are no fluid momentum equations
     */
    if (!pd->e[pg->imtrx][R_MOMENTUM1]) return;
    if (Current_EB_ptr->Elem_Blk_Id != i_mat_fluid)
      { 
	if (Current_EB_ptr->Elem_Blk_Id != i_mat_solid) EH(GOMA_ERROR, "Improper solid-fluid block ids");
	/* note, this may not account for all situations - we may not want
	   to quit here, but we'll do it for now */
	return;
      }

/* if this nodal contribution has already been added to fluid momentum equation (i.e.
   we are at a corner on the second side) return without doing anything */
    if (local_node_list_fs[id] == -1)
      {
	local_node_list_fs[id] = 1;
      }
    else
      {
	return;
      }

    /* check to make sure that both SOLID dof and velocity dof exist at this node */
    if (Dolphin[pg->imtrx][I][R_SOLID1] <= 0) return;
    if (Dolphin[pg->imtrx][I][R_MOMENTUM1] <= 0) return;

    /* loop over directions and add local contribution to liquid momentum into
     *  local contribution  for solid momentum while decrementing the liquid momentum
     */
    if ( af->Assemble_Residual )
      {
	for(p=0; p<dim; p++) {
	    /* don't do this if the SOLID displacement is fixed for this node */
	    if (! xfixed[p])
	      {
		ieqn_mom = R_MOMENTUM1 + p;
		ieqn_solid = R_SOLID1 + p;
		lec->R[upd->ep[pg->imtrx][ieqn_solid]][id_dofsol] += 
		    scale*lec->R[upd->ep[pg->imtrx][ieqn_mom]][id_dofmom];
/*fprintf(stderr,"solid_fluid %g %g %g\n",scale,fv->x[0],fv->P);*/
	      }
	}
      }
    
/* loop over directions and add local contribution to solid momentum into local contribution
 * for liquid momentum */
    if ( af->Assemble_Jacobian )
      {
	for(p=0; p<dim; p++) {
	    /* don't do this if the solid displacement is fixed for this node */
	    if (! xfixed[p])
	      {
		ieqn_mom = R_MOMENTUM1 + p;
		peqn_mom = upd->ep[pg->imtrx][ieqn_mom];
		ieqn_solid = R_SOLID1 + p;
		peqn_solid = upd->ep[pg->imtrx][ieqn_solid];	

		/* Add contributions due to all nodal sensitivities in solid element */

		/*
		 * local J_m_r -> J_r_r
		 */
		for ( q=0; q<dim; q++)
		  {
		    var = SOLID_DISPLACEMENT1+q;
                    pvar = upd->vp[pg->imtrx][var];
		    if ( pd->v[pg->imtrx][var] )
		      {
			for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
			  {
			    lec->J[peqn_solid][pvar][id_dofsol][j_id] +=
				scale*lec->J[peqn_mom][pvar][id_dofmom][j_id];
			  }
		      }
		  }

		/*
		 * local J_m_d -> J_r_d
		 */
		for ( q=0; q<dim; q++)
		  {
		    var = MESH_DISPLACEMENT1+q;
		    pvar = upd->vp[pg->imtrx][var];
		    if ( pd->v[pg->imtrx][var] )
		      {
			for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
			  {
			    lec->J[peqn_solid][pvar][id_dofsol][j_id] += 
				scale*lec->J[peqn_mom][pvar][id_dofmom][j_id];
			  }
		      }
		  }


		/*
		 * local J_m_P -> J_r_P
		 */
		var = PRESSURE;
		pvar = upd->vp[pg->imtrx][var];
		if ( pd->v[pg->imtrx][var] )
		  {
		    for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
		      {			
			lec->J[peqn_solid][pvar][id_dofsol][j_id] += 
			    scale*lec->J[peqn_mom][pvar][id_dofmom][j_id];
		      }
		  }
		

		/*
		 * local J_m_T -> J_r_T
		 */
		var = TEMPERATURE;
		pvar = upd->vp[pg->imtrx][var];
		if ( pd->v[pg->imtrx][var] )
		  {
		    for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
		      {
			lec->J[peqn_solid][pvar][id_dofsol][j_id] += 
			    scale*lec->J[peqn_mom][pvar][id_dofmom][j_id];
		      }
		  }
		
		/*
		 * local J_m_c -> J_r_c
		 */
		var = MASS_FRACTION;
		if ( pd->v[pg->imtrx][var] )
		  {
		    for ( w=0; w<pd->Num_Species_Eqn; w++)
		      {
		        pvar = upd->vp[pg->imtrx][var] + w;
			for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
			  {
			    lec->J[peqn_solid][pvar][id_dofsol][j_id] += 
                                scale*lec->J[peqn_mom][pvar][id_dofmom][j_id];
			  }
		      }
		  }
		
		/*
		 * local J_m_pmv -> J_r_pmv  (PRS: Don't think you have to split this
		 *                            loop up, even for new system. 
		 *                            cf. Lagrangian case above)
		 */
		for ( w=0; w<pd->Num_Porous_Eqn; w++)
		  {
                    var = POR_LIQ_PRES + w;
		    pvar = upd->vp[pg->imtrx][var];
		    if ( pd->v[pg->imtrx][var] )
		      {
			for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
			  {
			    lec->J[peqn_solid][pvar][id_dofsol][j_id] +=
				scale*lec->J[peqn_mom][pvar][id_dofmom][j_id];
			  }
		      }
		  }
		
		/*
		 * local J_m_v -> J_r_v  
		 */
		for ( q=0; q<dim; q++)
		  {
		    var = VELOCITY1+q;
		    if ( pd->v[pg->imtrx][var] )
		      {
			pvar = upd->vp[pg->imtrx][var];
			for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
			  {
			    lec->J[peqn_solid][pvar][id_dofsol][j_id] += 
				scale*lec->J[peqn_mom][pvar][id_dofmom][j_id];
			  }
		      }
		  }

	      } /* end of check for velocity dirichlet condition */
	}
      } /* End of Jacobian entries */
}
/* end of routine put_liquid_stress_in_solid_ALE */

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void 
penetration(double func[],
	    double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
	    double x_dot[DIM],  /* mesh velocity vector                      */
	    dbl tt,		/* parameter to vary time integration from 
				 * explicit (tt = 1) to implicit (tt = 0)    */
	    dbl dt,		/* current value of the time step            */
	    int bc_input_id,
	    struct Boundary_Condition *BC_Types,
	    int i_mat_solid,
	    int i_mat_fluid)

/*******************************************************************************
*
*  Function which evaluates the kinematic boundary condition 
*  n.(vfluid)=n.(vsolid) + (diffusive or porous flux) 
*  for adjacent fluid and solid (continuous or porous) phases
*			 Author: R. A. Cairncross (8/28/95)
*******************************************************************************/
{
  
/* Local variables */
  
  int j_id;
  int var,jvar,kdir=-1, a;
  int w, dim;
  double phi_j;

  double vconv[MAX_PDIM]; /*Calculated convection velocity */
  double vconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;

  int err;

  /* local contributions of boundary condition to residual and jacobian */

/***************************** EXECUTION BEGINS *******************************/
  if(af->Assemble_LSA_Mass_Matrix)
    return;

/*  initialize variables */
    dim = pd->Num_Dim;

/***************************** SOLID SIDE *******************************/
  /* 
   *  If current element block is the solid phase, calculate normal fluid flux
   *  in solid phase  - this is equal to the normal flux of the solid material
   *  plus the diffusive or porous flux of the solvents
   */
  if (Current_EB_ptr->Elem_Blk_Id == i_mat_solid)
    {
      /* get convection velocity in this phase */
      if (pd->MeshMotion == LAGRANGIAN ||
	  pd->MeshMotion == DYNAMIC_LAGRANGIAN)
	{
	  err = belly_flop(elc->lame_mu);
	  EH(err, "error in belly flop");
	}
	  err = get_convection_velocity(vconv, vconv_old, d_vconv, dt, tt);
      if ( mp->PorousMediaType == CONTINUOUS )
	{
	  for ( a=0; a<VIM; a++)
	    {
	      *func += fv->snormal[a] * vconv[a];
	    }

	  if (af->Assemble_Jacobian) {
	    
	    /* sum the contributions to the global stiffness matrix */
	    
	    var = TEMPERATURE;
	    for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
	      {
		if (pd->v[pg->imtrx][var])
		  {
		    phi_j = bf[var]->phi[j_id];
		    for(kdir=0; kdir<dim; kdir++) 
		      {
			/*     d( )/dx        */
			d_func[0][var][j_id] += d_vconv->T[kdir][j_id]*fv->snormal[kdir];
		      }
		  }
	      }
	    
	    var = MASS_FRACTION;
	    for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
	      {
		if (pd->v[pg->imtrx][var])
		  {
		    phi_j = bf[var]->phi[j_id];
		    for(kdir=0; kdir<dim; kdir++) 
		      {
			for(w=0; w<pd->Num_Species_Eqn; w++) 
			  {
			    /*     d( )/dx        */
			    d_func[0][MAX_VARIABLE_TYPES + w][j_id] += d_vconv->C[kdir][w][j_id]*fv->snormal[kdir];
			  }
		      }
		  }
	      }
	    
	    for(jvar=0; jvar<dim; jvar++) 
	      {
		var = MESH_DISPLACEMENT1 + jvar;
		for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
		  {
		    if (pd->v[pg->imtrx][var])
		      {
			phi_j = bf[var]->phi[j_id];
			for(kdir=0; kdir<dim; kdir++) 
			  {
			    /*     d( )/dx        */
			    d_func[0][var][j_id] += vconv[kdir]*fv->dsnormal_dx[kdir][jvar][j_id]
			      + d_vconv->X[kdir][jvar][j_id]*fv->snormal[kdir];
			  }
		      }
		  }
	      }
	    
	    for(jvar=0; jvar<dim; jvar++) 
	      {
		var = VELOCITY1 + jvar;
		for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
		  {
		    if (pd->v[pg->imtrx][var])
		      {
			phi_j = bf[var]->phi[j_id];
			d_func[0][var][j_id] += d_vconv->v[kdir][jvar][j_id] * fv->snormal[jvar];
		      }
		  }
	      }
	  }
	}
      else if (mp->PorousMediaType == POROUS_UNSATURATED || mp->PorousMediaType == POROUS_SATURATED)
	{
	  EH(GOMA_ERROR,"Not ready for porous penetration");
	}
      else  EH(GOMA_ERROR,"bad media type in penetration");
      
    }

/***************************** FLUID SIDE *******************************/
  /*
   *  If current material is the fluid phase, calculate normal fluid velocity
   *  in liquid at the interface
   */
  else if (Current_EB_ptr->Elem_Blk_Id == i_mat_fluid)
    {
      /* Calculate the residual contribution from the normal component of velocity */
      for(kdir=0; kdir<dim; kdir++) 
	{
	  *func += fv->v[kdir] * fv->snormal[kdir];
	}

      if (af->Assemble_Jacobian) {
	
	/* sum the contributions to the global stiffness matrix */
	
	for(jvar=0; jvar<dim; jvar++) 
	  {
	    var = MESH_DISPLACEMENT1 + jvar;
	    for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
	      {
		if (pd->v[pg->imtrx][var])
		  {
		    phi_j = bf[var]->phi[j_id];
		    for(kdir=0; kdir<dim; kdir++) 
		      {
			/*     d( )/dx        */
			d_func[0][var][j_id] += fv->v[kdir]*fv->dsnormal_dx[kdir][jvar][j_id];
		      }
		  }
	      }
	  }
	
	for(jvar=0; jvar<dim; jvar++) 
	  {
	    var = VELOCITY1 + jvar;
	    for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
	      {
		if (pd->v[pg->imtrx][var])
		  {
		    phi_j = bf[var]->phi[j_id];
		    d_func[0][var][j_id] += phi_j * fv->snormal[jvar];
		  }
	      }
	  }
      }

    }
  else EH(GOMA_ERROR,"Penetration called with incorrect block id");
  /* note, we may not want to quit at this point */
  return;
}

/****************************************************************************/

void 
no_slip(double func[],
	double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
	double x_dot[DIM],     /* mesh velocity vector                      */
        double x_rs_dot[DIM],  /* real solid velocity vector                     */
	dbl tt,		      /* parameter to vary time integration from 
			       * explicit (tt = 1) to implicit (tt = 0)    */
	dbl dt,		       /* current value of the time step            */
	int bc_input_id,
	struct Boundary_Condition *BC_Types,
	int eb_mat_solid,
	int eb_mat_fluid)
/******************************************************************************
*
*  Function which evaluates the kinematic boundary condition 
*  n.(vfluid)=n.(vsolid) + (diffusive or porous flux) 
*  for adjacent fluid and solid (continuous or porous) phases
*			 Author: R. A. Cairncross (8/28/95)
******************************************************************************/
{
  int j, j_id;
  int var = -1, jvar = -1, a;
  int w, dim;
  double phi_j;

  double vconv[MAX_PDIM]; /*Calculated convection velocity */
  double vconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;

  double v_solid_mesh[DIM];
  int err;

  /* local contributions of boundary condition to residual and jacobian */

/***************************** EXECUTION BEGINS *******************************/
/*  initialize variables */
    dim = pd->Num_Dim;

    for(a=0; a<DIM; a++) v_solid_mesh[a]=0.;

/***************************** SOLID SIDE *******************************/
    if (af->Assemble_LSA_Mass_Matrix)
      {
	if (Current_EB_ptr->Elem_Blk_Id == eb_mat_solid)
	  {
	    if ( pd->MeshMotion != LAGRANGIAN &&
		 pd->MeshMotion != DYNAMIC_LAGRANGIAN &&
	         pd->MeshMotion != TOTAL_ALE)
	      EH(GOMA_ERROR, "Shouldn't be in this section of no-slip with an arbitrary solid");
	    if ( mp->PorousMediaType == CONTINUOUS )
	      for (jvar = 0; jvar < dim; jvar++)
		{
		  if(pd->MeshMotion == LAGRANGIAN ||
		     pd->MeshMotion == DYNAMIC_LAGRANGIAN) 
		    var = MESH_DISPLACEMENT1 + jvar; 
		  else if (pd->MeshMotion == TOTAL_ALE)
		    var = SOLID_DISPLACEMENT1 + jvar;
		  else
		    EH(GOMA_ERROR, "Bad pd->MeshMotion");
		  if (pd->v[pg->imtrx][var])
		    for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		      {
			phi_j = bf[var]->phi[j];
			d_func[jvar][var][j] -= phi_j;
		      }
		}
	  }
	return;
      }

  /*
   *  If current material is the solid phase, calculate normal fluid flux
   *  in solid phase  - this is equal to the normal flux of the solid material
   *  plus the diffusive or porous flux of the solvents
   */
  if (Current_EB_ptr->Elem_Blk_Id == eb_mat_solid)
    {
      /*first load up the correct reference velocity and do this in the
       solid phase because you may not have it in the liquid phase, as in
       TALE */
      if(pd->MeshMotion == LAGRANGIAN ||
	  pd->MeshMotion == DYNAMIC_LAGRANGIAN)
	{
	  for(a=0; a<DIM; a++) v_solid_mesh[a]=x_dot[a];
	  err = belly_flop(elc->lame_mu);
	  EH(err, "error in belly flop");
	  err = get_convection_velocity(vconv, vconv_old, d_vconv, dt, tt);

	}
      else if (pd->MeshMotion == TOTAL_ALE)
	{
	  for(a=0; a<DIM; a++) v_solid_mesh[a]=x_rs_dot[a];
	  err = belly_flop_rs(elc_rs->lame_mu);
	  EH(err, "error in belly flop");
	  err = get_convection_velocity_rs(vconv, vconv_old, d_vconv, dt, tt);

	}
      else
	{
	  EH(GOMA_ERROR,"Shouldn't be in this section of no-slip with an arbitrary solid");
	}


      if ( mp->PorousMediaType == CONTINUOUS )
	{
	  for ( a=0; a<VIM; a++)
	    {
	      func[a] -= vconv[a] + v_solid_mesh[a];
	    }

	  if (af->Assemble_Jacobian) {

	    /* sum the contributions to the global stiffness matrix */

	    var = TEMPERATURE;
	    for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
	      {
		if (pd->v[pg->imtrx][var])
		  {
		    phi_j = bf[var]->phi[j_id];
			/*     d( )/dx        */
		    for ( a=0; a<VIM; a++)
		      {
			d_func[a][var][j_id] -= d_vconv->T[a][j_id];
		      }
		  }
	      }

	    var = MASS_FRACTION;
	    for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
	      {
		if (pd->v[pg->imtrx][var])
		  {
		    phi_j = bf[var]->phi[j_id];
		    for ( a=0; a<VIM; a++)
		      {
			for(w=0; w<pd->Num_Species_Eqn; w++)
			  {
			    /*     d( )/dx        */
			    d_func[a][MAX_VARIABLE_TYPES + w][j_id] -= d_vconv->C[a][w][j_id];
			  }
		      }
		  }
	      }

	    for(jvar=0; jvar<dim; jvar++)
	      {
		var = MESH_DISPLACEMENT1 + jvar;
		for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
		  {
		    if (pd->v[pg->imtrx][var])
		      {
			phi_j = bf[var]->phi[j_id];
			for ( a=0; a<VIM; a++)
			  {
			    /*     d( )/dx        */
			    d_func[a][var][j_id] -= d_vconv->X[a][jvar][j_id];
			  }
		      }
		  }
	      }

	    if (pd->MeshMotion == TOTAL_ALE)
	      {
		for(jvar=0; jvar<dim; jvar++)
		  {
		    var = SOLID_DISPLACEMENT1 + jvar;
		    for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
		      {
			if (pd->v[pg->imtrx][var])
			  {
			    phi_j = bf[var]->phi[j_id];
			    for ( a=0; a<VIM; a++)
			      {
				/*     d( )/dx        */
				d_func[a][var][j_id] -= d_vconv->rs[a][jvar][j_id];
			      }
			 }
		      }
		  }
	      }

	    for(jvar=0; jvar<dim; jvar++)
	      {
		var = VELOCITY1 + jvar;
		for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
		  {
		    if (pd->v[pg->imtrx][var])
		      {
			phi_j = bf[var]->phi[j_id];
			d_func[jvar][var][j_id] -= d_vconv->v[jvar][jvar][j_id];
		      }
		  }
	      }

	    for (jvar=0; jvar<dim; jvar++)
	      {
		if(pd->MeshMotion == LAGRANGIAN ||
		   pd->MeshMotion == DYNAMIC_LAGRANGIAN)
		  {
		    var = MESH_DISPLACEMENT1 + jvar;
		  }
		else if (pd->MeshMotion == TOTAL_ALE)
		  {
		    var = SOLID_DISPLACEMENT1 + jvar;
		  }
		if (pd->v[pg->imtrx][var])
		  {
		    for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		      {
			phi_j = bf[var]->phi[j];
			if(TimeIntegration != 0)
			  {
			    d_func[jvar][var][j] += ( -(1.+2.*tt)* phi_j/dt);
			  }
		      }
		  }
	      }
	  }
	}
      else if (mp->PorousMediaType == POROUS_UNSATURATED || mp->PorousMediaType == POROUS_SATURATED)
	{
	  EH(GOMA_ERROR,"Not ready for porous penetration");
	}
      else  EH(GOMA_ERROR,"bad media type in penetration");

    }

/***************************** FLUID SIDE *******************************/
  /*
   *  If current material is the fluid phase, calculate normal fluid velocity
   *  in liquid at the interface
   */
  else if (Current_EB_ptr->Elem_Blk_Id == eb_mat_fluid)
    {


      /* Calculate the residual contribution from the velocity */
      for ( a=0; a<VIM; a++)
	{
	  func[a] += fv->v[a];
	}

      if (af->Assemble_Jacobian) {

	/* sum the contributions to the global stiffness matrix */

	for(jvar=0; jvar<dim; jvar++)
	  {
	    var = VELOCITY1 + jvar;
	    for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
	      {
		if (pd->v[pg->imtrx][var])
		  {
		    phi_j = bf[var]->phi[j_id];
		    d_func[jvar][var][j_id] += phi_j;
		  }
	      }
	  }
      }
    }
  else EH(GOMA_ERROR,"No Slip called with incorrect block id");
  /* note, we may not want to quit at this point */
  return;
}

/*------------------------------------------------------------*/
/* mesh_stress_tensor() -- find TT (stress tensor) for LAGRANGIAN
 * or arbitrary mesh motion
 *
 * Created: 1/16/95 by RAC
 *
 */

int

mesh_stress_tensor(dbl TT[DIM][DIM],
		   dbl dTT_dx[DIM][DIM][DIM][MDE],
		   dbl dTT_dp[DIM][DIM][MDE],
		   dbl dTT_dc[DIM][DIM][MAX_CONC][MDE],
		   dbl dTT_dp_liq[DIM][DIM][MDE],
		   dbl dTT_dp_gas[DIM][DIM][MDE],
		   dbl dTT_dporosity[DIM][DIM][MDE],
		   dbl dTT_dsink_mass[DIM][DIM][MDE],
 		   dbl dTT_dT[DIM][DIM][MDE],
 		   dbl dTT_dmax_strain[DIM][DIM][MDE],
                   dbl dTT_dcur_strain[DIM][DIM][MDE],
		   dbl mu,
		   dbl lambda,
		   dbl delta_t,
		   int ielem,	/* current element number */
		   int ip,	/* current integration point */
		   int ip_total) /* number of gauss points in this element */
{
 dbl factor = 0.0;
 dbl d_mu_dx[DIM][MDE], d_lambda_dx[DIM][MDE];
 int a=0, b, j, p, q, dim, var, w, v, dofs, err, mat_ielem;
 int SPECIES = MAX_VARIABLE_TYPES;
 dbl p_gas_star = 0.0;

  dbl thermexp;
  dbl speciesexp[MAX_CONC];
  dbl d_thermexp_dx[MAX_VARIABLE_TYPES+MAX_CONC];
  dbl d_speciesexp_dx[MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC];

  dim = ei[pg->imtrx]->ielem_dim;
  mat_ielem = PRS_mat_ielem;

  
  memset(d_mu_dx,0,sizeof(double)*DIM*MDE);
  memset(d_lambda_dx,0,sizeof(double)*DIM*MDE);
  memset(d_thermexp_dx,0,sizeof(double)*(MAX_VARIABLE_TYPES+MAX_CONC));
  memset(d_speciesexp_dx,0,sizeof(double)*MAX_CONC*(MAX_VARIABLE_TYPES+MAX_CONC));
  memset(speciesexp,0,sizeof(double)*MAX_CONC);
  memset(TT,0,sizeof(dbl)*DIM*DIM);

/* 
 * Calculate the lame coefficients if they are not constant. Note, some of the
 * sensitivities of mu and lambda to porous media vars come out through the elc->d_lame
 * structure elements, and hence are not args. 
 */

  err = load_elastic_properties(elc, &mu, &lambda, &thermexp, speciesexp, d_mu_dx, d_lambda_dx, d_thermexp_dx, d_speciesexp_dx);
  EH(err," Problem in loading up elastic constants");


  /* Here we will simple use our cadre of Elastic models if no Viscoplastic
   * strain is allowed, otherwise we will call the EVP routine get_evp_stress_tensor
   */

  if (evpl->ConstitutiveEquation == NO_MODEL)
    {
      for ( p=0; p<VIM; p++)
	{
	  for ( q=0; q<VIM; q++)
	    {
	      TT[p][q] = lambda * fv->volume_strain * delta(p,q) + 2. * mu * fv->strain[p][q];
	    }
	}

      /* add shrinkage stress, if called for */
      if(elc->thermal_expansion_model == SHRINKAGE)
	{
	  if((fv->external_field[0] >= 1.63 && fv->external_field[0] <= 1.7) ||
	     Element_Blocks[ei[pg->imtrx]->elem_blk_index].ElemStorage[mat_ielem].solidified[ip])
	    {
	      for ( p=0; p<VIM; p++)
		{
		  for ( q=0; q<VIM; q++)
		    {
		      TT[p][q] -=  (2.* mu + 3.*lambda) * (-0.04) * delta(p,q);
		    }
		}
	      Element_Blocks[ei[pg->imtrx]->elem_blk_index].ElemStorage[mat_ielem].solidified[ip] = 1.0;
	    }
	}
 
      /*  add thermo-elasticity  */
      if( pd->e[pg->imtrx][R_ENERGY] )
 	{
	  if( elc->thermal_expansion_model == CONSTANT)
	    {
	      for ( p=0; p<VIM; p++)
		{
		  for ( q=0; q<VIM; q++)
		    {
		      TT[p][q] -=  (2.* mu + 3.*lambda) * thermexp * 
			(fv->T - elc->solid_reference_temp) * delta(p,q);
		    }
		}
	    }
	  if( elc->thermal_expansion_model == IDEAL_GAS)
	    {
	      for ( p=0; p<VIM; p++)
		{
		  for ( q=0; q<VIM; q++)
		    {
		      TT[p][q] -=  (2.* mu + 3.*lambda) / (thermexp + fv->T) * 
			(fv->T - elc->solid_reference_temp) * delta(p,q);
		    }
		}
	    }
	  if( elc->thermal_expansion_model == USER)
	    {
	      for ( p=0; p<VIM; p++)
		{
		  for ( q=0; q<VIM; q++)
		    {
		      TT[p][q] -=  (2.* mu + 3.*lambda) * thermexp * delta(p,q);
		    }
		}
	    }
	}
      
/*   add species expansion/shrinkage	*/
      if (pd->e[pg->imtrx][R_MASS])
	{
	     for (w=0; w<pd->Num_Species_Eqn; w++) 
	       {
		 if(mp->SpecVolExpModel[w] == CONSTANT ||
                     mp->SpecVolExpModel[w] == PHOTO_CURING)
		   {
		     for ( p=0; p<VIM; p++)
		       {
			 for ( q=0; q<VIM; q++)
			   {
			     TT[p][q] -=  (2.* mu + 3.*lambda) * speciesexp[w] * 
			       (fv->c[w] - mp->reference_concn[w]) * delta(p,q);
			   }
		       }
		   }
	       }
	   }

   if ( af->Assemble_Jacobian )
     {
       for ( p=0; p<VIM; p++)
	 {
	   for ( q=0; q<VIM; q++)
	     {
	       for ( b=0; b<dim; b++)
		 {
		   v = MESH_DISPLACEMENT1 + b;
		   if ( pd->v[pg->imtrx][v] )
		     {
		       dofs     = ei[pg->imtrx]->dof[v];
		       for ( j=0; j<dofs; j++)
			 {
			   dTT_dx[p][q][b][j] =
			     lambda * fv->d_volume_strain_dx[b][j] * delta(p, q)
			     + 2. * mu * fv->d_strain_dx[p][q][b][j];
			   
			   dTT_dx[p][q][b][j] += d_lambda_dx[b][j] * fv->volume_strain * delta(p,q) 
			     + 2. * d_mu_dx[b][j] * fv->strain[p][q];
			   if( pd->e[pg->imtrx][R_ENERGY] )
			     {
	  			if( elc->thermal_expansion_model == CONSTANT || 
                                    elc->thermal_expansion_model == IDEAL_GAS )
	  			{
			       dTT_dx[p][q][b][j] -=  (2. * d_mu_dx[b][j] + 3.*d_lambda_dx[b][j])
				 * thermexp * 
				 (fv->T - elc->solid_reference_temp) * delta(p,q);
	  			}
	  			if( elc->thermal_expansion_model == USER)
	  			{
				  dTT_dx[p][q][b][j] -=  ((2.*d_mu_dx[b][j]+3.*d_lambda_dx[b][j])*thermexp 
							  + (2.*mu+3.*lambda)*d_thermexp_dx[v]*bf[v]->phi[j])
				    * delta(p,q);
	  			}
			     }
			   if( pd->e[pg->imtrx][R_MASS] )
			     {
			       for (w=0; w<pd->Num_Species_Eqn; w++) 
				 {
				   if(mp->SpecVolExpModel[w] == CONSTANT ||
                                        mp->SpecVolExpModel[w] == PHOTO_CURING)
				     {
				       dTT_dx[p][q][b][j] -=  (2. * d_mu_dx[b][j] + 3.*d_lambda_dx[b][j])
					 * speciesexp[w] * (fv->c[w] - mp->reference_concn[w]) 
					 * delta(p,q);
				     }
				 }
			     }
			 }
		     }
		 }
	       /*  Temperature sensitivities here */
	       v = TEMPERATURE;
	       if ( pd->v[pg->imtrx][v] )
		 {
		   dofs     = ei[pg->imtrx]->dof[v];

		   for (j=0; j<dofs; j++)
		     {
		       /*if no temperature dependence d_lame_mu is zero because we initialize in set_mp_to_unity */
		       dTT_dT[p][q][j] += ( 2.*elc->d_lame_mu[TEMPERATURE]*fv->strain[p][q]      )*bf[v]->phi[j];
		     }

		   if( elc->thermal_expansion_model == CONSTANT)
		     {
		       for ( j=0; j<dofs; j++)
			 {
			   dTT_dT[p][q][j] -=  (2.* mu + 3.*lambda) * thermexp * bf[v]->phi[j] 
			     * delta(p,q); 
			 }
		     }
		   if( elc->thermal_expansion_model == IDEAL_GAS)
		     {
		       for ( j=0; j<dofs; j++)
			 {
			   dTT_dT[p][q][j] -=  (2.* mu + 3.*lambda) * 
                                    (thermexp+elc->solid_reference_temp)/SQUARE(fv->T+thermexp)
                                      * bf[v]->phi[j] * delta(p,q); 
			 }
		     }
		   if( elc->thermal_expansion_model == USER)
		     {
		       for ( j=0; j<dofs; j++)
			 {
			   dTT_dT[p][q][j] -=  (2.* mu + 3.*lambda) * d_thermexp_dx[v]*bf[v]->phi[j] 
			     * delta(p,q); 
			 }
		     }

		   if(elc->lame_mu_model == USER)
		     {
		       dTT_dT[p][q][j] += ( 2.*elc->d_lame_mu[TEMPERATURE]*fv->strain[p][q]      )*bf[v]->phi[j];
		     }
		   if(elc->lame_lambda_model == USER)
		     {
		       dTT_dT[p][q][j] += (elc->d_lame_lambda[TEMPERATURE]*fv->volume_strain*delta(p,q))*bf[v]->phi[j];
		     }
		 }
	       /*  max_strain sensitivities here */
	       v = MAX_STRAIN;
	       if ( pd->v[pg->imtrx][v] )
		 {
		   dofs     = ei[pg->imtrx]->dof[v];

		   if(elc->lame_mu_model == TABLE)
		     {
		       for ( j=0; j<dofs; j++)
			 {
			   dTT_dmax_strain[p][q][j] += fv->volume_strain * delta(p,q) * elc->d_lame_lambda[v] * bf[v]->phi[j];
			   dTT_dmax_strain[p][q][j] += 2.0 * fv->strain[p][q] * elc->d_lame_mu[v] * bf[v]->phi[j];
			 }
		     }
		 }
	       /*  cur_strain sensitivities here */
	       v = CUR_STRAIN;
	       if ( pd->v[pg->imtrx][v] )
		 {
		   dofs     = ei[pg->imtrx]->dof[v];

		   if(elc->lame_mu_model == TABLE)
		     {
		       for ( j=0; j<dofs; j++)
			 {
			   dTT_dcur_strain[p][q][j] += fv->volume_strain * delta(p,q) * elc->d_lame_lambda[v] * bf[v]->phi[j];
			   dTT_dcur_strain[p][q][j] += 2.0 * fv->strain[p][q] * elc->d_lame_mu[v] * bf[v]->phi[j];
			 }
		     }
		 }
	       /*  Species expansion sensitivities here */
 	       v = MASS_FRACTION;
 	       if ( pd->v[pg->imtrx][v] )
 		 {
 		   dofs     = ei[pg->imtrx]->dof[v];
 		   for ( j=0; j<dofs; j++)
 		     {
		       for (w=0; w<pd->Num_Species_Eqn; w++) 
			 {
			   if(mp->SpecVolExpModel[w] == CONSTANT ||
                                 mp->SpecVolExpModel[w] == PHOTO_CURING)
			     {
			       dTT_dc[p][q][w][j] -=  (2.* mu + 3.*lambda)*
				 speciesexp[w]*bf[v]->phi[j]*delta(p,q); 
			     }
			 }
 		     }
		   if( elc->thermal_expansion_model == USER)
		     {
		       for ( j=0; j<dofs; j++)
			 {
			   for (w=0; w<pd->Num_Species_Eqn; w++) 
			     {
			       dTT_dc[p][q][w][j] -=  (2.* mu + 3.*lambda)*
				 d_thermexp_dx[MAX_VARIABLE_TYPES+w]*bf[v]->phi[j]
				 *delta(p,q); 
			     }
			 }
		     }
 		 }
	     }
	 }
       /* PRESSURE and MASS_FRACTION sensitivities are taken care of later */
     }
    }
  else if(evpl->ConstitutiveEquation == EVP_HYPER)
    {
      err =get_evp_stress_tensor(TT, 
				 dTT_dx, 
				 dTT_dp, 
				 dTT_dc,
				 mu,
				 lambda, 
				 delta_t,  
				 ielem, 
				 ip, 
				 ip_total);
      /*If your elastic model calls for an compressible formulation, add in volume strain term */
      if (cr->MeshFluxModel != INCOMP_PSTRAIN 
	  && cr->MeshFluxModel != INCOMP_PSTRESS && cr->MeshFluxModel != INCOMP_3D)
	{
	  for ( p=0; p<VIM; p++)
	    {
	      for ( q=0; q<VIM; q++)
		{
		  TT[p][q] += lambda * fv->volume_strain * delta(p,q);
		}
	    }
	
	  if ( af->Assemble_Jacobian )
	    {
	      for ( p=0; p<VIM; p++)
		{
		  for ( q=0; q<VIM; q++)
		    {
		      for ( b=0; b<dim; b++)
			{
			  v = MESH_DISPLACEMENT1 + b;
			  if ( pd->v[pg->imtrx][v] )
			    {
			      dofs     = ei[pg->imtrx]->dof[v];
			      for ( j=0; j<dofs; j++)
				{
				  dTT_dx[p][q][b][j] +=
				    lambda * fv->d_volume_strain_dx[b][j] * delta(p, q);
			      
				  dTT_dx[p][q][b][j] += d_lambda_dx[b][j] * fv->volume_strain * delta(p,q);
				}
			    }
			}
		    }
		}
	      /* PRESSURE and MASS_FRACTION sensitivities are taken care of later */
	    }
	}
    }
  else
    {
      EH(GOMA_ERROR,"Unrecognizable Plastic Constitutive Equation specification");
    }

  if (cr->MeshMotion == LAGRANGIAN ||
      cr->MeshMotion == DYNAMIC_LAGRANGIAN)
   {
     /*
      * add on the pressure contribution to the stress tensor for the porous or 
      * incompressible formulations
      */
     /* now, add in pressure due to swelling for incompressible lagrangian
      * mesh motion */
     var = PRESSURE;
     if (pd->v[pg->imtrx][var] && ( mp->PorousMediaType == CONTINUOUS ) ) 
       {
	 if (cr->MeshFluxModel == INCOMP_PSTRAIN ||
	     cr->MeshFluxModel == INCOMP_PSTRESS || 
	     cr->MeshFluxModel == INCOMP_3D      ||
	     cr->MeshFluxModel == LINEAR)
	   {
	     for ( a=0; a<DIM; a++)
	       {
		 TT[a][a] -=  fv->P; 
	       }
	   }
       }

     /* pressure force is the pressure of each phase times its volume fraction */
     if ((mp->PorousMediaType == POROUS_UNSATURATED || 
	 mp->PorousMediaType == POROUS_SATURATED || 
	 mp->PorousMediaType == POROUS_TWO_PHASE) && pd->e[pg->imtrx][R_POR_POROSITY]) 
       {
	 for ( a=0; a<VIM; a++)
	   {
	     if (mp->CapStress == NO_CAP_STRESS)
	       {
		 /*do nothing*/
	       }
	     else if (mp->CapStress == COMPRESSIBLE && mp->PorousMediaType != POROUS_SATURATED) 
	       {
		 EH(GOMA_ERROR,"Need to reconcile the COMPRESSIBLE model pressures.  Not ready yet");
		 /* Compressible model of effective stress law with
		  *  partially saturated media from Zienkeivicz and Garg and Nur */
		 if (elc->lame_lambda_model != POWER_LAW) 
		   EH(GOMA_ERROR,"Effective stress law may be missing constant");

		 TT[a][a] -=  (1 - (1 - mp->porosity) * elc->lame_lambda / elc->u_lambda[0]) * 
		   mp->saturation * fv->p_liq; 
		 if (pd->v[pg->imtrx][POR_GAS_PRES])
		   {
		     TT[a][a] -=  (1 - (1 - mp->porosity) 
				   * elc->lame_lambda / elc->u_lambda[0]) * 
		       (1. - mp->saturation) * fv->p_gas; 
	       
		   }
	       } 
	     else if (mp->CapStress == PARTIALLY_WETTING && mp->PorousMediaType != POROUS_SATURATED) 
	       {
		 p_gas_star =  mp->u_porous_gas_constants[3];
		 TT[a][a] -= (1. - mp->saturation)*p_gas_star + mp->saturation * fv->p_liq;

		 if (pd->v[pg->imtrx][POR_GAS_PRES]) TT[a][a] -= (1. - mp->saturation) * fv->p_gas; 
	       }
	       
	     else if (mp->CapStress == WETTING && mp->PorousMediaType != POROUS_SATURATED) 
	       {
		 /* If liquid is wetting, so that all surfaces are covered
		    by a thin layer of liquid */
		 EH(GOMA_ERROR,"Need to reconcile the WETTING model pressures.  Not ready yet");
		 TT[a][a] -= (1 - mp->porosity * (1. - mp->saturation)) * fv->p_liq;  
		 if (pd->v[pg->imtrx][POR_GAS_PRES]) TT[a][a] -= mp->porosity * (1. - mp->saturation) * fv->p_gas;  
	       }
	     else if (mp->PorousMediaType == POROUS_SATURATED && pd->e[pg->imtrx][POR_POROSITY])
	       { 
		 TT[a][a] -= fv->p_liq; 
	       
	       }
	     else 
	       {
		 WH(-1,"No way to put liquid stress into porous matrix because you have const porosity and/or no pore eqn");
	       }

	     if(pd->e[pg->imtrx][R_POR_SINK_MASS])
	       {

		 /*This factor is used to partition the strain of an agm particle between the network
		  *and the pore-space.  factor=0 implies porespace only. 
		  */
		 factor = mp->u_porous_sink_constants[7];  

		 TT[a][a] -= factor*2*mu*fv->sink_mass*mp->u_porous_sink_constants[5]/mp->density; 
	       }

	   }
       }
 
     if ( af->Assemble_Jacobian )
       {
	 var = PRESSURE; 
	 if (pd->v[pg->imtrx][var])
	   {
	     if (cr->MeshFluxModel == INCOMP_PSTRESS || cr->MeshFluxModel == HOOKEAN_PSTRESS )
	       { 
		 for ( p=0; p<VIM; p++)
		   {
		     for ( q=0; q<VIM; q++)
		       {
			 for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
			   {
			     /* didn't add this sensitivity in before - add it now */
			     dTT_dp[p][q][j] = 
			       2. * mu * fv->d_strain_dp[p][q][j]
			       + lambda * fv->d_volume_strain_dp[j] * delta(p,q);
			   }
		       }
		   }
	       }

	     if ((cr->MeshFluxModel == INCOMP_PSTRAIN 
		  || cr->MeshFluxModel == INCOMP_PSTRESS || cr->MeshFluxModel == INCOMP_3D || cr->MeshFluxModel == LINEAR) 
		 && ( mp->PorousMediaType == CONTINUOUS ) )
	       {
		 for ( a=0; a<VIM; a++)
		   {
		     for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		       {
			 dTT_dp[a][a][j] -= bf[var]->phi[j]; 
		       }
		   }
	       }
	   }
	 
     
 /* calculate sensitivity to concentration due to concentration dependent
    elastic and bulk modulus */

	 var = MASS_FRACTION;
	 if (pd->v[pg->imtrx][var])
	   {

	     /* sensitivity of total stress to species concentrations
		- first calculate contribution from sensitivity of elastic moduli */
	     for (w=0; w<pd->Num_Species_Eqn; w++) 
	       {
		 for ( a=0; a<VIM; a++)
		   {
		     for ( b=0; b<VIM; b++)
		       {
			 for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
			   {
			     /* didn't add this in before - add it in now */
			     dTT_dc[a][b][w][j] += bf[var]->phi[j] 
			       * (
				  elc->d_lame_lambda[SPECIES + w] * fv->volume_strain * delta(a,b)
				  + 2. * elc->d_lame_mu[SPECIES + w] * fv->strain[a][b] );
			   }
		       }
		   }
	       }
	   }
 
	 /* sensitivity of total stress to porous media variables
	    - first calculate contribution from sensitivity of elastic moduli */
	 var = POR_POROSITY;
	 if (pd->v[pg->imtrx][var])
	   for ( a=0; a<VIM; a++)
	     {
	       for ( b=0; b<VIM; b++)
		 {
		   for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		     {
		       /* didn't add this in before - add it in now */
		       dTT_dporosity[a][b][j] += bf[var]->phi[j] 
			 * (
			    elc->d_lame_lambda[POR_POROSITY] * fv->volume_strain * delta(a,b)
			    + 2. * elc->d_lame_mu[POR_POROSITY] * fv->strain[a][b] );
		     }
		 }
	     }
 
	 /* in porous media, stress has contribution from phase pressures */
	 /*Let's deal with liquid phase pressure first */
	 var = POR_LIQ_PRES;
	 if (pd->v[pg->imtrx][var] && pd->e[pg->imtrx][R_POR_POROSITY])
	   {
		 for ( p=0; p<VIM; p++)
		   {
		     for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		       {
			 if (mp->CapStress == COMPRESSIBLE) 
			   {
			     dTT_dp_liq[p][p][j] -= bf[var]->phi[j] *
			       (1 - (1 - mp->porosity) * elc->lame_lambda / elc->u_lambda[0]) * 
			       ( mp->saturation + mp->d_saturation[POR_LIQ_PRES] * fv->p_liq);
 
			     if (pd->v[pg->imtrx][POR_GAS_PRES]) 
			       {
				 dTT_dp_gas[p][p][j] -= bf[var]->phi[j] * 
				   (1 - (1 - mp->porosity) * elc->lame_lambda / elc->u_lambda[0]) * 
				   ((1. - mp->saturation) - mp->d_saturation[POR_GAS_PRES] * fv->p_gas); 
			       }
			   if (pd->v[pg->imtrx][POR_POROSITY])  
			     {
			       dTT_dporosity[p][p][j] -= bf[var]->phi[j] 
				 * (
				    mp->d_porosity[POR_POROSITY] * elc->lame_lambda 
				    / elc->u_lambda[0] * mp->saturation
				    - (1 - mp->porosity) * elc->d_lame_lambda[POR_POROSITY] 
				    / elc->u_lambda[0] * mp->saturation
				    + (1 - (1 - mp->porosity) * elc->lame_lambda / elc->u_lambda[0]) 
				    * mp->d_saturation[POR_POROSITY]
				    )  * fv->p_liq; 
			     }
			   if (pd->v[pg->imtrx][POR_GAS_PRES] && pd->v[pg->imtrx][POR_POROSITY] )  
			     {
			       dTT_dporosity[p][p][j] -= bf[var]->phi[j]
				 * ( 
				    mp->d_porosity[POR_POROSITY] 
				    * elc->lame_lambda / elc->u_lambda[0] * (1 - mp->saturation)
				    - (1 - (1 - mp->porosity) * elc->lame_lambda / elc->u_lambda[0]) 
				    * mp->d_saturation[POR_POROSITY]
				    ) * fv->p_gas;
			     }

			   }
			 else if (mp->CapStress == PARTIALLY_WETTING)
			   {
			     dTT_dp_liq[p][p][j] -= 
			       -bf[var]->phi[j] *
			       p_gas_star*mp->d_saturation[POR_LIQ_PRES]
			       +bf[var]->phi[j] *
			       (mp->saturation + mp->d_saturation[POR_LIQ_PRES] *  fv->p_liq ); 

			     if (pd->v[pg->imtrx][POR_GAS_PRES]) 
			       {
				 dTT_dp_gas[p][p][j] -= bf[var]->phi[j] * 
				   ( 
				    (1. - mp->saturation) + mp->d_saturation[POR_GAS_PRES] * 
				    (fv->p_liq - fv->p_gas) );
			       }
			   if (pd->v[pg->imtrx][POR_POROSITY])  
			     {
			       dTT_dporosity[p][p][j] -= bf[var]->phi[j]
						* mp->d_saturation[POR_POROSITY] * fv->p_liq;
			     }
			   if (pd->v[pg->imtrx][POR_POROSITY] && pd->v[pg->imtrx][POR_GAS_PRES] ) 
			     {
			       dTT_dporosity[p][p][j] += bf[var]->phi[j]*
				 mp->d_saturation[POR_POROSITY] * fv->p_gas;
			     }

			   }
			 else if (mp->CapStress == WETTING)
			   {
			   /* If liquid is wetting, so that all surfaces are covered 
			      by a thin layer of liquid */
			     dTT_dp_liq[p][p][j] -= bf[var]->phi[j]
			       * ((1 -  mp->porosity * (1. - mp->saturation))
				  + mp->porosity * mp->d_saturation[POR_LIQ_PRES] * fv->p_liq );
			     if (pd->v[pg->imtrx][POR_GAS_PRES]) 
			       {
				 dTT_dp_gas[p][p][j] -= bf[var]->phi[j] * mp->porosity * 
				   (
				    (1. - mp->saturation) - mp->d_saturation[POR_GAS_PRES] *fv->p_gas );
			       }
			     if (pd->v[pg->imtrx][POR_POROSITY])  
			       {
				 dTT_dporosity[p][p][j] -= bf[var]->phi[j] 
				   * (
				      - mp->d_porosity[POR_POROSITY] * (1. - mp->saturation) 
				      + mp->porosity * mp->d_saturation[POR_POROSITY] 
				      ) * fv->p_liq;
			       }
			     if (pd->v[pg->imtrx][POR_GAS_PRES] && pd->v[pg->imtrx][POR_POROSITY] ) 
			       {
				 dTT_dporosity[p][p][j] -= bf[var]->phi[j]
				   * (
				      mp->d_porosity[POR_POROSITY] * (1. - mp->saturation) 
				      - mp->porosity * mp->d_saturation[POR_POROSITY] 
				      ) * fv->p_gas;
			       }

			   }
			 else if (mp->PorousMediaType == POROUS_SATURATED && pd->e[pg->imtrx][POR_POROSITY]) 
			   { 
			     dTT_dp_liq[p][p][j] -= bf[var]->phi[j]; 
			   }
		       }

		   }
	   }

	 var = POR_SINK_MASS;
	 if (pd->v[pg->imtrx][var] && pd->e[pg->imtrx][R_POR_POROSITY])
	   {
	     for ( p=0; p<VIM; p++)
	       {
		 for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		   {
		     /* dTT_dsink_mass[p][p][j] -= factor*2*mu*bf[var]->phi[j]*mp->u_porous_sink_constants[5]/mp->density; */
		     dTT_dsink_mass[p][p][j] -= factor*2*mu*bf[var]->phi[j] *
		       (mp->u_porous_sink_constants[2]/mp->u_porous_sink_constants[5])/(1.0 - fv->porosity);
		   }
	       }

	     var = POR_POROSITY;
	     if(pd->v[pg->imtrx][var])
	       {
		 for ( p=0; p<VIM; p++)
		   {
		     for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		       {
		 
			 dTT_dporosity[p][p][j] -= factor*2*mu*bf[var]->phi[j] *
			   ((1.0 - mp->u_porosity[0]) + fv->sink_mass*mp->u_porous_sink_constants[2]/mp->u_porous_sink_constants[5])
			   /(1.0 - fv->porosity)/(1.0 - fv->porosity);
		       }
		   }
	       }
	   }
		     

       } /* end of if Jacobian */
   } /* end of if LAGRANGIAN */
 
 return(0);
} /* end of mesh_stress_tensor()*/

/*------------------------------------------------------------*/
/* get_evp_stress_tensor() -- find deviatoric part of TT 
 * (stress tensor) for LAGRANGIAN, elasto-viscoplastic contribution
 *
 * Created: 2/20/99 by PRS/SYT
 *
 */

int 
get_evp_stress_tensor(double TT[DIM][DIM],
		      double dTT_dx[DIM][DIM][DIM][MDE],
		      double dTT_dp[DIM][DIM][MDE],
		      double dTT_dc[DIM][DIM][MAX_CONC][MDE],
		      double mu,
		      double lambda,
		      double delta_t,
		      int ielem, /* current element number */
		      int ip,	/* current integration point */
		      int ip_total) /* number of gauss points in element    */
{

  double F[DIM][DIM];  /* deformation gradient tensor with the gradient operator
                       based on the undeformed coordinates */
  double F_inv[DIM][DIM];  /* the inverse of the deformation gradient tensor */
  double F_vp[DIM][DIM];  /* the viscoplastic deformation gradient tensor */
  double F_vp_old[DIM][DIM];  /* the viscoplastic deformation gradient tensor of
                              previous converged time step */
  double F_vp_inv[DIM][DIM];  /* inverse of viscoplastic deformation gradient */
  double dF_dx[DIM][DIM][DIM][MDE];  
  double dF_vp_inv_dx[DIM][DIM][DIM][MDE];  
  double dF_vp_inv_dc[DIM][DIM][MAX_CONC][MDE];  
  double dF_vp_dx[DIM][DIM][DIM][MDE];  
  double dF_vp_dc[DIM][DIM][MAX_CONC][MDE];  
  double dF_e_dx[DIM][DIM][DIM][MDE];  
  double dF_e_dc[DIM][DIM][MAX_CONC][MDE];
  double F_e[DIM][DIM];  /* the elastic deformation gradient tensor */
  double C_e[DIM][DIM];  /* product of (transpose F_e) by F_e */
  double E_e[DIM][DIM];  /* the elastic strain tensor */
  double dE_e_dx[DIM][DIM][DIM][MDE];
  double dE_e_dc[DIM][DIM][MAX_CONC][MDE];
  double S[DIM][DIM];  /* the second Piola-Kirchoff stress tensor, i.e.
                       the stress tensor based on the stress-free
                       frame of reference */
  double alpha;          /* the linear shrinkage parameter that depends on solvent
                         concentration */
  int k;

 double grad_d[DIM][DIM]; /* displacement gradient */

 double d_mu_dx[DIM][MDE], d_lambda_dx[DIM][MDE];
 double d_plastic_mu_dc[MAX_CONC][MDE];
 double d_yield_dc[MAX_CONC][MDE];
 int a, b, i, j = -1, p, q, m, n, dim, v, v1, var, dofs, dofs1, err, F_vp_flag;
 double thermexp;
 double speciesexp[MAX_CONC];
 double d_thermexp_dx[MAX_VARIABLE_TYPES+MAX_CONC];
 double d_speciesexp_dx[MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC];

 dim = ei[pg->imtrx]->ielem_dim;
 if(dim > 2) EH(GOMA_ERROR,"EVP models only implemented for plane strain case just now");

 /******************************************************/
 /* Update glossary */
 /* evpl_glob[0]->update_flag=0 -- first time in here */
 /* evpl_glob[0]->update_flag=1 -- all other times in here except...*/
 /* evpl_glob[0]->update_flag=2 -- in here but failed on the last time step to converge*/
 /*****************************************************/

/*Initialize this element for SFS case, viz F_vp=I  and T=0, or with restart
information */

 if (**evpl->F_vp_old_glob[ielem][ip] == 0.) 
   {
     for (a = 0; a < efv->Num_external_field; a++)
       {
	 if(!strcmp(efv->name[a], "FVP11"))
	   {
	     evpl->F_vp_old_glob[ielem][ip][0][0] = fv->external_field[a];
		   
	   }
	 if(!strcmp(efv->name[a], "FVP22"))
	   {
	     evpl->F_vp_old_glob[ielem][ip][1][1] = fv->external_field[a];
	   }
	 if(!strcmp(efv->name[a], "FVP12"))
	   {
	     evpl->F_vp_old_glob[ielem][ip][0][1] = fv->external_field[a];
		   
	   }
	 if(!strcmp(efv->name[a], "FVP21"))
	   {
	     evpl->F_vp_old_glob[ielem][ip][1][0] = fv->external_field[a];
		   
	   }
	 if(!strcmp(efv->name[a], "FVP13"))
	   {
	     evpl->F_vp_old_glob[ielem][ip][0][2] = fv->external_field[a];
	     
       }
	 if(!strcmp(efv->name[a], "FVP31"))
	   {
	     evpl->F_vp_old_glob[ielem][ip][2][0] = fv->external_field[a];
		   
	   }
	 if(!strcmp(efv->name[a], "FVP23"))
	   {
	     evpl->F_vp_old_glob[ielem][ip][1][2] = fv->external_field[a];
		   
	   }
	 if(!strcmp(efv->name[a], "FVP32"))
       {
	 evpl->F_vp_old_glob[ielem][ip][2][1] = fv->external_field[a];
		   
       }
	 if(!strcmp(efv->name[a], "FVP33"))
	   {
	     evpl->F_vp_old_glob[ielem][ip][2][2] = fv->external_field[a];
		   
	   }
	 if(!strcmp(efv->name[a], "TVP11"))
	   {
	     evpl->TT_glob[ielem][ip][0][0] = fv->external_field[a];
	     evpl->TT_old_glob[ielem][ip][0][0] = fv->external_field[a];
	     
	   }
	 if(!strcmp(efv->name[a], "TVP22"))
	   {
	     evpl->TT_glob[ielem][ip][1][1] = fv->external_field[a];
	     evpl->TT_old_glob[ielem][ip][1][1] = fv->external_field[a];
		   
	   }
	 if(!strcmp(efv->name[a], "TVP12"))
	   {
	     evpl->TT_glob[ielem][ip][0][1] = fv->external_field[a];
	     evpl->TT_old_glob[ielem][ip][0][1] = fv->external_field[a];
		   
	   }
	 if(!strcmp(efv->name[a], "TVP21"))
	   {
	     evpl->TT_glob[ielem][ip][1][0] = fv->external_field[a];
	     evpl->TT_old_glob[ielem][ip][1][0] = fv->external_field[a];
		   
	   }
	 if(!strcmp(efv->name[a], "TVP13"))
	   {
	     evpl->TT_glob[ielem][ip][0][2] = fv->external_field[a];
	     evpl->TT_old_glob[ielem][ip][0][2] = fv->external_field[a];
	   }
	 if(!strcmp(efv->name[a], "TVP31"))
	   {
	     evpl->TT_glob[ielem][ip][2][0] = fv->external_field[a];
	     evpl->TT_old_glob[ielem][ip][2][0] = fv->external_field[a];
	     
	   }
	 if(!strcmp(efv->name[a], "TVP23"))
	   {
	     evpl->TT_glob[ielem][ip][1][2] = fv->external_field[a];
	     evpl->TT_old_glob[ielem][ip][1][2] = fv->external_field[a];
	     
	   }
	 if(!strcmp(efv->name[a], "TVP32"))
	   {
	     evpl->TT_glob[ielem][ip][2][1] = fv->external_field[a];
	     evpl->TT_old_glob[ielem][ip][2][1] = fv->external_field[a];
	     
	   }
	 if(!strcmp(efv->name[a], "TVP33"))
	   {
	     evpl->TT_glob[ielem][ip][2][2] = fv->external_field[a];
	     evpl->TT_old_glob[ielem][ip][2][2] = fv->external_field[a];
	     
	   }
       }
   }
 /*if we have gotten to here and still have no vp strain, then just set to SFS */
 if (**evpl->F_vp_old_glob[ielem][ip] == 0.) 
   {
     for (p=0; p<DIM; p++)
       {
	 for (q=0; q<DIM; q++)
	   {
	     evpl->F_vp_old_glob[ielem][ip][p][q] = delta(p,q);
	     evpl->TT_glob[ielem][ip][p][q] = 0.0;
	     evpl->TT_old_glob[ielem][ip][p][q] = 0.0;

	     for ( b=0; b<DIM; b++)
	       {
			 
		 v1 = MESH_DISPLACEMENT1 + b;
		 if ( pd->v[pg->imtrx][v1] )
		   {
		     dofs1    = ei[pg->imtrx]->dof[v1];
		     for ( m=0; m<dofs1; m++)
		       {
			 evpl->dTT_dx_glob[ielem][ip][p][q][b][m] = 0.;
			 evpl->dTT_dx_old_glob[ielem][ip][p][q][b][m] = 0.;
		       }
		   }
	       }
	   }
       }
   }


  memset(d_mu_dx,0,sizeof(double)*DIM*MDE);
  memset(d_lambda_dx,0,sizeof(double)*DIM*MDE);
  memset(d_plastic_mu_dc,0,sizeof(double)*MAX_CONC*MDE);
  memset(d_yield_dc,0,sizeof(double)*MAX_CONC*MDE);
  memset(grad_d,0,sizeof(double)*DIM*DIM);
  memset(F_e,0,sizeof(double)*DIM*DIM);
  memset(S, 0, sizeof(double)*DIM*DIM);
  memset(C_e,0,sizeof(double)*DIM*DIM);
  memset(F_vp,0,sizeof(double)*DIM*DIM);
  memset( F_vp_inv, 0, sizeof(double)*DIM*DIM); 
  memset( dF_dx, 0, sizeof(double)*DIM*DIM*DIM*MDE);  
  memset( dF_vp_inv_dx,0, sizeof(double)*DIM*DIM*DIM*MDE);  
  memset( dF_vp_inv_dc, 0,  sizeof(double)*DIM*DIM*MAX_CONC*MDE);  
  memset( dF_vp_dx, 0, sizeof(double)*DIM*DIM*DIM*MDE); 
  memset( dF_vp_dc, 0, sizeof(double)*DIM*DIM*MAX_CONC*MDE); 
  memset(d_thermexp_dx,0,sizeof(double)*(MAX_VARIABLE_TYPES+MAX_CONC));
  memset(d_speciesexp_dx,0,sizeof(double)*MAX_CONC*(MAX_VARIABLE_TYPES+MAX_CONC));
  memset(speciesexp,0,sizeof(double)*MAX_CONC);

/* 
 * Calculate the lame coefficients if they are not constant
 * Probably redundant here.  Already called up in mesh_stress_tensor
 */

  err = load_elastic_properties(elc, &mu, &lambda, &thermexp, speciesexp, d_mu_dx, d_lambda_dx, d_thermexp_dx, d_speciesexp_dx);
  EH(err," Problem in loading up elastic constants");

  /* will not need plastic_mu nor yield from the routine;
     otherwise, UMR results */

  err = load_plastic_properties(d_plastic_mu_dc, d_yield_dc);
  EH(err," Problem in loading up plastic constants");

	 
/* Update TT, and F_vp */
   if (evpl_glob[0]->update_flag == 1)
     {
       for ( i=0; i<Num_Elem; i++)
	 {
         for ( j=0; j<ip_total; j++)
	   {
           for ( p=0; p<DIM; p++)
             for ( q=0; q<DIM; q++)
               {
                 evpl->TT_old_glob[i][j][p][q] = evpl->TT_glob[i][j][p][q];
                 evpl->F_vp_old_glob[i][j][p][q] = evpl->F_vp_glob[i][j][p][q];
                 for ( b=0; b<DIM; b++)
                   {
                     v1 = MESH_DISPLACEMENT1 + b;
                     if ( pd->v[pg->imtrx][v1] )
                       {
                         dofs1    = ei[pg->imtrx]->dof[v1];
                         for ( m=0; m<dofs1; m++)
                           {
		 	     evpl->dTT_dx_old_glob[i][j][p][q][b][m] = evpl->dTT_dx_glob[i][j][p][q][b][m];
			   }
		       }
		   } 
                 var = MASS_FRACTION;
	         for (b=0; b<pd->Num_Species_Eqn; b++) 
	           {
	             for ( m=0; m<ei[pg->imtrx]->dof[var]; m++)
	               {
		         evpl->dTT_dc_old_glob[i][j][p][q][b][m] = evpl->dTT_dc_glob[i][j][p][q][b][m];
	               }
	           }
               }
	   }
	 }
       evpl_glob[0]->update_flag = 0;
     }

   /* Get the true-stress tensor in elasto-viscoplastic deformation modeling
      to replace that calculated by RAC. KSC on 6-5-98  */

   /* First calculate the displacement gradient, which is based on the current
      configuration */

      for ( p=0; p<VIM; p++)
        {
          for ( q=0; q<VIM; q++)
            {
              grad_d[p][q] = fv->grad_d[p][q];
            }
        }

  /* Calculate the inverse of the deformation gradient tensor */
  for ( p=0; p<VIM; p++)
    {
      for ( q=0; q<VIM; q++)
        {
          F_inv[p][q] = delta(p,q) - grad_d[p][q];
        }
    }

/* determine the deformation gradient tensor from its inverse */
   invert_tensor(F_inv, F, VIM, NULL, NULL, 0, 0);
   
/* Because we are dealing with plane stress case */
   F[0][2] = 0;
   F[1][2] = 0;
   F[2][0] = 0;
   F[2][1] = 0;
   F[2][2] = 1;

/* Unpack global variables into local variables */
   for ( p=0; p<DIM; p++)
     for ( q=0; q<DIM; q++)
       {
         if (evpl_glob[0]->update_flag != 2)
           {
             TT[p][q] = evpl->TT_glob[ielem][ip][p][q];
             for ( b=0; b<DIM; b++)
               {
                 v1 = MESH_DISPLACEMENT1 + b;
                 if ( pd->v[pg->imtrx][v1] )
                   {
                     dofs1    = ei[pg->imtrx]->dof[v1];
                     for ( m=0; m<dofs1; m++)
                       {
		         dTT_dx[p][q][b][m] = evpl->dTT_dx_glob[ielem][ip][p][q][b][m];
	               }
	           }
	       } 
             var = MASS_FRACTION;
	     for (b=0; b<pd->Num_Species_Eqn; b++) 
	       {
	         for ( m=0; m<ei[pg->imtrx]->dof[var]; m++)
	           {
		     dTT_dc[p][q][b][m] = evpl->dTT_dc_glob[ielem][ip][p][q][b][m];
	           }
	       }
           }
         else
           {
             TT[p][q] = evpl->TT_old_glob[ielem][ip][p][q];
             for ( b=0; b<DIM; b++)
               {
                 v1 = MESH_DISPLACEMENT1 + b;
                 if ( pd->v[pg->imtrx][v1] )
                   {
                     dofs1    = ei[pg->imtrx]->dof[v1];
                     for ( m=0; m<dofs1; m++)
                       {
		         dTT_dx[p][q][b][m] = evpl->dTT_dx_old_glob[ielem][ip][p][q][b][m];
		       }
		   }
	       } 
             var = MASS_FRACTION;
	     for (b=0; b<pd->Num_Species_Eqn; b++) 
	       {
	         for ( m=0; m<ei[pg->imtrx]->dof[var]; m++)
	           {
		     dTT_dc[p][q][b][m] = evpl->dTT_dc_old_glob[ielem][ip][p][q][b][m];
	           }
	       }
           }
         F_vp_old[p][q] = evpl->F_vp_old_glob[ielem][ip][p][q];
       }


/* calculate the viscoplastic deformation gradient tensor */
  F_vp_flag=get_F_vp(F_vp, F_vp_old, TT, dTT_dx, dTT_dc, dF_vp_dx, dF_vp_dc, d_plastic_mu_dc, d_yield_dc, delta_t);

/* calculate the isotropic shrinkage parameter */
  alpha = pow(fv->volume_change, 1./3.); 

/* invert the viscoplastic deformation gradient tensor only if material yields */
  if (F_vp_flag ==1)
    {
       invert_tensor(F_vp, F_vp_inv, VIM, NULL, NULL, 0, 0);
    }

/* calculate the elastic deformation gradient tensor */
  for ( p=0; p<VIM; p++)
    {
      for ( q=0; q<VIM; q++)
        {
	  if (F_vp_flag==1)
	    {
	      F_e[p][q] = 0.0;
	      for ( k=0; k<VIM; k++)
		{
		  F_e[p][q] += (1.0/alpha)*F[p][k]*F_vp_inv[k][q];
		}
	    }
	  else
	    {
               F_e[p][q] = (1.0/alpha)*F[p][q];
	    }
	}
    }

/* calculate the strain tensor based on the stress-free reference frame */
  for ( p=0; p<VIM; p++)
    {
      for ( q=0; q<VIM; q++)
        {
          C_e[p][q] = 0.0;
          for ( k=0; k<VIM; k++)
            {
              C_e[p][q] += F_e[k][p]*F_e[k][q];
            }
          E_e[p][q]=0.5*(C_e[p][q] - delta(p,q));
        }
    }

/* calculate the stress tensor */

  for ( p=0; p<VIM; p++)
    {
      for ( q=0; q<VIM; q++)
        {
         S[p][q] = 2.0*elc->lame_mu*E_e[p][q];  
        }
    }
/* calculate the true stress tensor, i.e. the stress tensor based on
   the current configuration */

  for ( p=0; p<VIM; p++)
    {
      for ( q=0; q<VIM; q++)
        {
         TT[p][q]=S[p][q];
        }
    }

  if ( af->Assemble_Jacobian )
 
    {
      for (p=0; p<VIM; p++)
        {
           for (q=0; q<VIM; q++)
             {
	        for ( b=0; b<dim; b++)
	          {
	            v = MESH_DISPLACEMENT1 + b;
	            if ( pd->v[pg->imtrx][v] )
		      {
		        dofs     = ei[pg->imtrx]->dof[v];
		        for ( j=0; j<dofs; j++)
		          {
			    dF_dx[p][q][b][j] = 0.0;
                            for (m=0; m<VIM; m++)
                   	      {  
                                for (n=0; n<VIM; n++)
                   	          {  
                                    dF_dx[p][q][b][j] += F[p][m] * fv->d_grad_d_dmesh[m][n][b][j] *  F[n][q]; 
			          }
			      }
		          }
		      }
	          }
              }
         }

      for ( b=0; b<dim; b++)
        {
          v = MESH_DISPLACEMENT1 + b;
          if ( pd->v[pg->imtrx][v] )
	    {
	      dofs     = ei[pg->imtrx]->dof[v];
	      for ( j=0; j<dofs; j++)
	        {
                  dF_dx[0][2][b][j] = 0.;
                  dF_dx[1][2][b][j] = 0.;
                  dF_dx[2][0][b][j] = 0.;
                  dF_dx[2][1][b][j] = 0.;
                  dF_dx[2][2][b][j] = 0.;
		}
            }
        }

      if (F_vp_flag == 1)
	{
          for (p=0; p<VIM; p++)
            {
              for (q=0; q<VIM; q++)
                {
	          for ( b=0; b<dim; b++)
	            {
	              v = MESH_DISPLACEMENT1 + b;
	              if ( pd->v[pg->imtrx][v] )
		        {
		          dofs     = ei[pg->imtrx]->dof[v];
		          for ( j=0; j<dofs; j++)
			    {
			      dF_vp_inv_dx[p][q][b][j] = 0.;
	      		      for ( m=0; m<VIM; m++)
	        	        {
	      		          for ( n=0; n<VIM; n++)
	        	            {
                                      dF_vp_inv_dx[p][q][b][j] -= F_vp_inv[p][m] * dF_vp_dx[m][n][b][j] * F_vp_inv[n][q];
	        	            }
		   	        }
		    	    }
		        }
		    }
	        }
            }
          var=MASS_FRACTION;
          for (p=0; p<VIM; p++)
            {
              for (q=0; q<VIM; q++)
                {
	          for (b=0; b<pd->Num_Species_Eqn; b++) 
	            {
	              for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		        {
			  dF_vp_inv_dc[p][q][b][j] = 0.;
	      		  for ( m=0; m<VIM; m++)
	        	    {
	      		      for ( n=0; n<VIM; n++)
	        	        {
                                  dF_vp_inv_dc[p][q][b][j] -= F_vp_inv[p][m] * dF_vp_dc[m][n][b][j] * F_vp_inv[n][q];
		   	        }
		    	    }
		        }
		    }
	        }
            }
        }

      for (p=0; p<VIM; p++)
        {
          for (q=0; q<VIM; q++)
            {
	      for ( b=0; b<dim; b++)
	        {
	          v = MESH_DISPLACEMENT1 + b;
	          if ( pd->v[pg->imtrx][v] )
		    {
		      dofs     = ei[pg->imtrx]->dof[v];
		      for ( j=0; j<dofs; j++)
			{
		          if (F_vp_flag == 1)
			    {
			      dF_e_dx[p][q][b][j]=0.;
	      		      for ( m=0; m<VIM; m++)
	        	        {
                                  dF_e_dx[p][q][b][j] += -1./3.*pow(fv->volume_change,-4./3.)*fv->d_volume_change_dx[b][j] * F[p][m] * F_vp_inv[m][q] + 
			          1./alpha * dF_dx[p][m][b][j] * F_vp_inv[m][q] + 
			          1./alpha * F[p][m] * dF_vp_inv_dx[m][q][b][j]; 
				}
	        	    }
		          else
		            {
                          dF_e_dx[p][q][b][j] = -1./3.*pow(fv->volume_change,-4./3.)*fv->d_volume_change_dx[b][j] * F[p][q] + 1./alpha * dF_dx[p][q][b][j]; 
			    }
		        }
		    }
	        }
            }
        }

      if (F_vp_flag == 1)
	{
          for (p=0; p<VIM; p++)
            {
              for (q=0; q<VIM; q++)
                {
	          for (b=0; b<pd->Num_Species_Eqn; b++) 
	            {
	              for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	                {
			  dF_e_dc[p][q][b][j]=0.;
	      		  for ( m=0; m<VIM; m++)
	        	    {
                              dF_e_dc[p][q][b][j] += -1./3.*pow(fv->volume_change,-4./3.)*fv->d_volume_change_dx[b][j] * F[p][m] * F_vp_inv[m][q] + 
			      1./alpha * F[p][m] * dF_vp_inv_dc[m][q][b][j]; 
	        	    }
			}
		    }
	        }
            }
        }
      for (p=0; p<VIM; p++)
        {
          for (q=0; q<VIM; q++)
             {
	       for ( b=0; b<dim; b++)
	         {
	           v = MESH_DISPLACEMENT1 + b;
	           if ( pd->v[pg->imtrx][v] )
		     {
		       dofs     = ei[pg->imtrx]->dof[v];
		       for ( j=0; j<dofs; j++)
		         {
			   dE_e_dx[p][q][b][j] = 0.0;
                           for (i=0; i<VIM; i++)
                   	     {  
                               dE_e_dx[p][q][b][j] += 0.5 * (dF_e_dx[i][p][b][j] * F_e[i][q] + F_e[i][p] * dF_e_dx[i][q][b][j]);
			     }
		         }
		     }
	         }
             }
         }

      if (F_vp_flag == 1)
	{
          for (p=0; p<VIM; p++)
            {
              for (q=0; q<VIM; q++)
                {
	           for (b=0; b<pd->Num_Species_Eqn; b++) 
	             {
	               for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	                 {
		           dE_e_dc[p][q][b][j] = 0.0;
                           for (i=0; i<VIM; i++)
                   	     {  
                               dE_e_dc[p][q][b][j] += 0.5 * (dF_e_dc[i][p][b][j] * F_e[i][q] + F_e[i][p] * dF_e_dc[i][q][b][j]);
			     }
		         }
		     } 
	         }
             }
         }
       for (p=0; p<VIM; p++)
         {
           for (q=0; q<VIM; q++)
             {
	       for ( b=0; b<dim; b++)
	         {
	           v = MESH_DISPLACEMENT1 + b;
	           if ( pd->v[pg->imtrx][v] )
		     {
		       dofs     = ei[pg->imtrx]->dof[v];
		       for ( j=0; j<dofs; j++)
		         {
                           dTT_dx[p][q][b][j] = 2.0 * elc->lame_mu * dE_e_dx[p][q][b][j] + 2.0 * d_mu_dx[b][j] * E_e[p][q]; 
		         }
		     } 
	         }
             }
         }

       for (p=0; p<VIM; p++)
         {
           for (q=0; q<VIM; q++)
             {
	       for (b=0; b<pd->Num_Species_Eqn; b++) 
	         {
	           for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	             {
		       if (F_vp_flag ==1)
			 {
                           dTT_dc[p][q][b][j] = 2.0 * elc->lame_mu * dE_e_dc[p][q][b][j];
			 }
		       else
			 {
                           dTT_dc[p][q][b][j] = 0.; 
			 }
		     } 
	         }
             }
         }

    }  /* End if Assemble Jacobian */

 
/* Load local variables into global storage */
   for ( p=0; p<VIM; p++)
     for ( q=0; q<VIM; q++)
       {
         evpl->TT_glob[ielem][ip][p][q] = TT[p][q];
         evpl->F_vp_glob[ielem][ip][p][q] = F_vp[p][q];
       }
   for (p=0; p<VIM; p++)
     {
       for (q=0; q<VIM; q++)
         {
           for ( b=0; b<dim; b++)
	     {
	       v = MESH_DISPLACEMENT1 + b;
	       if ( pd->v[pg->imtrx][v] )
	         {
		   dofs     = ei[pg->imtrx]->dof[v];
		   for ( j=0; j<dofs; j++)
		     {
                       evpl->dTT_dx_glob[ielem][ip][p][q][b][j] = dTT_dx[p][q][b][j]; 
		     }
		 } 
	     }
           var = MASS_FRACTION;
	   for (b=0; b<pd->Num_Species_Eqn; b++) 
	     {
	       for ( m=0; m<ei[pg->imtrx]->dof[var]; m++)
	         {
                   evpl->dTT_dc_glob[ielem][ip][p][q][b][j] = dTT_dc[p][q][b][j]; 
	         }
	     }
         }
     }

 return(0);
} /* end of get_evp_stress_tensor()*/

/******************************************************************************/

/*
 *  This function calculates the viscoplastic deformation gradient.
 *  If the viscoplastic deformation gradient is nonzero,
 *  it returns 1, otherwise, it returns 0.
 *
 *  Created by S. Y. Tam 07/23/1998
 *  Modified to include sensitibities SYT 4/2000
 */

int 
get_F_vp(double F_vp[DIM][DIM],
	 double F_vp_old[DIM][DIM],
	 double TT[DIM][DIM],
	 double dTT_dx[DIM][DIM][DIM][MDE],
	 double dTT_dc[DIM][DIM][MAX_CONC][MDE],
	 double dF_vp_dx[DIM][DIM][DIM][MDE],
	 double dF_vp_dc[DIM][DIM][MAX_CONC][MDE],
  	 double d_plastic_mu_dc[MAX_CONC][MDE],
  	 double d_yield_dc[MAX_CONC][MDE],
	 double delta_t)
{
  double TTddC, phi, phe;
  int i, j, p, q, F_vp_flag;
  int b, dim, v, dofs, var;
  double T[DIM][DIM];
  double d_phi_dTT[DIM][DIM];
  double d_phi_dTT_norm[DIM][DIM];
  double D_vp[DIM][DIM];
  double norm = 0.0;
  double dT2_dTT_dx[DIM][DIM][DIM][MDE];
  double dT2_dTT_dc[DIM][DIM][MAX_CONC][MDE];
  double d_phi_dTT_norm_dx[DIM][DIM][DIM][MDE];
  double d_phi_dTT_norm_dc[DIM][DIM][MAX_CONC][MDE];
  double d_phi_dTT_dx[DIM][DIM][DIM][MDE];
  double d_phi_dTT_dc[DIM][DIM][MAX_CONC][MDE];
  double dT2_dTT[DIM][DIM];
  double d_phi_dx[DIM][MDE];
  double d_phi_dc[MAX_CONC][MDE];
  double d_norm2_dx[DIM][MDE];
  double d_norm2_dc[MAX_CONC][MDE];

  F_vp_flag = 0;

  dim = ei[pg->imtrx]->ielem_dim;

  TTddC = 0.;
  for ( p=0; p<VIM; p++)
    {
      TTddC += TT[p][p];
    }

  for ( p=0; p<VIM; p++)
    {
      for ( q=0; q<VIM; q++)
        {
          T[p][q] = TT[p][q] - 1./3. * TTddC * delta(p,q);
        }
    }

  phe = 0.;
  for ( p=0; p<VIM; p++)
    {
      for ( q=0; q<VIM; q++)
        {
          phe += T[q][p]*T[p][q];
        }
    }

  phe = sqrt(1.5*phe);
  phi = phe - evpl->yield;

/* calculate only if the material has yielded */
  if (phi>0.) 
   {
     for ( p=0; p<VIM; p++)
       {
         for ( q=0; q<VIM; q++)
           {
             d_phi_dTT[p][q] = 0.;
             for ( i=0; i<VIM; i++)
               {
                 for ( j=0; j<VIM; j++)
		   {
                     d_phi_dTT[p][q] += T[i][j]*(delta(p,i)*delta(q,j) - 1./3.*delta(p,q)*delta(i,j)); 
                   }
               }
             d_phi_dTT[p][q] = 1.5/phe * d_phi_dTT[p][q];
           }
       }

/* normalize d_phi_dTT */
     norm = 0.;
     for ( p=0; p<VIM; p++)
       {
         for ( q=0; q<VIM; q++)
           {  
	     norm += d_phi_dTT[p][q] * d_phi_dTT[p][q];
           } 
       }
     norm = sqrt(norm);
     if (norm < 1.e-12) norm = 1.;

     for ( p=0; p<VIM; p++)
       {
         for ( q=0; q<VIM; q++)
           {
	      d_phi_dTT_norm[p][q] = d_phi_dTT[p][q]/norm;
           }
       }

/* calculate the viscoplastic strain rate tensor */
     for ( p=0; p<VIM; p++)
       {
         for ( q=0; q<VIM; q++)
           {
             D_vp[p][q] = 1/evpl->plastic_mu * phi * d_phi_dTT_norm[p][q];
             F_vp_flag = 1;
           }
       }

     for ( p=0; p<VIM; p++)
       {
         for ( q=0; q<VIM; q++)
           {
             F_vp[p][q] = 0;
             F_vp[p][q] = D_vp[p][q] * delta_t + F_vp_old[p][q];
           }
       } 
    }
  else
    {
      for ( p=0; p<VIM; p++)
	{
	  for ( q=0; q<VIM; q++)
	    {
	      D_vp[p][q] = 0.0;
              F_vp[p][q] = F_vp_old[p][q];
            }
        }
    }


  if ( af->Assemble_Jacobian )
    {
      var = MASS_FRACTION;
      if (F_vp_flag == 1)
	{

	  for ( b=0; b<dim; b++)
	    {
	      v = MESH_DISPLACEMENT1 + b;
	      if ( pd->v[pg->imtrx][v] )
	        {
	          dofs     = ei[pg->imtrx]->dof[v];
	          for ( j=0; j<dofs; j++)
	            {
	              dT2_dTT_dx[0][0][b][j] = 2./3. * (2. * dTT_dx[0][0][b][j] 
					       		- dTT_dx[1][1][b][j] 
							- dTT_dx[2][2][b][j]);
	              dT2_dTT_dx[0][1][b][j] = 2. * dTT_dx[0][1][b][j]; 
	              dT2_dTT_dx[0][2][b][j] = 2. * dTT_dx[0][2][b][j];
	              dT2_dTT_dx[1][0][b][j] = 2. * dTT_dx[1][0][b][j];
	              dT2_dTT_dx[1][1][b][j] = 2./3. * (-dTT_dx[0][0][b][j] + 
					       		2. * dTT_dx[1][1][b][j] 
					      		- dTT_dx[2][2][b][j]);
	              dT2_dTT_dx[1][2][b][j] = 2. * dTT_dx[1][2][b][j];
	              dT2_dTT_dx[2][0][b][j] = 2. * dTT_dx[2][0][b][j];
	              dT2_dTT_dx[2][1][b][j] = 2. * dTT_dx[2][1][b][j];
	              dT2_dTT_dx[2][2][b][j] = 2./3. * (-dTT_dx[0][0][b][j] - 
					       		dTT_dx[1][1][b][j] +
					      		2. * dTT_dx[2][2][b][j]);
		    }
	        }
            }
	 
	 if (pd->v[pg->imtrx][var]) 
	   {
	     for (b=0; b<pd->Num_Species_Eqn; b++) 
	       {
	         for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		   {
	             dT2_dTT_dc[0][0][b][j] = 2./3. * (2. * dTT_dc[0][0][b][j] 
				       		- dTT_dc[1][1][b][j] 
						- dTT_dc[2][2][b][j]);
	             dT2_dTT_dc[0][1][b][j] = 2. * dTT_dc[0][1][b][j]; 
	             dT2_dTT_dc[0][2][b][j] = 2. * dTT_dc[0][2][b][j];
	             dT2_dTT_dc[1][0][b][j] = 2. * dTT_dc[1][0][b][j];
	             dT2_dTT_dc[1][1][b][j] = 2./3. * (-dTT_dc[0][0][b][j] + 
				       		2. * dTT_dc[1][1][b][j] 
				      		- dTT_dc[2][2][b][j]);
	             dT2_dTT_dc[1][2][b][j] = 2. * dTT_dc[1][2][b][j];
	             dT2_dTT_dc[2][0][b][j] = 2. * dTT_dc[2][0][b][j];
	             dT2_dTT_dc[2][1][b][j] = 2. * dTT_dc[2][1][b][j];
	             dT2_dTT_dc[2][2][b][j] = 2./3. * (-dTT_dc[0][0][b][j] - 
				       		dTT_dc[1][1][b][j] +
				      		2. * dTT_dc[2][2][b][j]);
	           }
               }
	   }

          dT2_dTT[0][0] = 2./3. * (2. * TT[0][0] - TT[1][1] - TT[2][2]);
          dT2_dTT[0][1] = 2. * TT[0][1]; 
          dT2_dTT[0][2] = 2. * TT[0][2];
          dT2_dTT[1][0] = 2. * TT[1][0];
          dT2_dTT[1][1] = 2./3. * (-TT[0][0] + 2. * TT[1][1] - TT[2][2]);
          dT2_dTT[1][2] = 2. * TT[1][2];
          dT2_dTT[2][0] = 2. * TT[2][0];
          dT2_dTT[2][1] = 2. * TT[2][1];
          dT2_dTT[2][2] = 2./3. * (-TT[0][0] - TT[1][1] + 2. * TT[2][2]);

	  for ( b=0; b<dim; b++)
	    {
	      v = MESH_DISPLACEMENT1 + b;
	      if ( pd->v[pg->imtrx][v] )
	        {
	          dofs     = ei[pg->imtrx]->dof[v];
	          for ( j=0; j<dofs; j++)
	            {
		      d_phi_dx[b][j] = 0.5/phe * 
			((2.*TT[0][0]-TT[1][1]-TT[2][2]) * dTT_dx[0][0][b][j] +
			(-TT[0][0]+2.*TT[1][1]-TT[2][2]) * dTT_dx[1][1][b][j] +
			(-TT[0][0]-TT[1][1]+2.*TT[2][2]) * dTT_dx[2][2][b][j]) +
			3./phe * 
			(TT[0][1] * dTT_dx[0][1][b][j] + 
			TT[0][2] * dTT_dx[0][2][b][j] + 
			TT[1][2] * dTT_dx[1][2][b][j]); 
		    }
	        }
            }
	 
	   if (pd->v[pg->imtrx][var])
	     {
	       for (b=0; b<pd->Num_Species_Eqn; b++) 
	         {
	           for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		     {
		       d_phi_dc[b][j] = 0.5/phe * 
			((2.*TT[0][0]-TT[1][1]-TT[2][2]) * dTT_dc[0][0][b][j] +
			(-TT[0][0]+2.*TT[1][1]-TT[2][2]) * dTT_dc[1][1][b][j] +
			(-TT[0][0]-TT[1][1]+2.*TT[2][2]) * dTT_dc[2][2][b][j]) +
			3./phe * 
			(TT[0][1] * dTT_dc[0][1][b][j] + 
			TT[0][2] * dTT_dc[0][2][b][j] + 
			TT[1][2] * dTT_dc[1][2][b][j]); 
		     }
	         }
             }
	 
	  for ( b=0; b<dim; b++)
	    {
	      v = MESH_DISPLACEMENT1 + b;
	      if ( pd->v[pg->imtrx][v] )
	        {
	          dofs     = ei[pg->imtrx]->dof[v];
	          for ( j=0; j<dofs; j++)
	            {
		      d_norm2_dx[b][j] = -3. * pow(phe,-3.) * d_phi_dx[b][j]*
			(TT[0][0]*TT[0][0] + TT[1][1]*TT[1][1] + TT[2][2]*TT[2][2]
			-TT[0][0]*TT[1][1]-TT[0][0]*TT[2][2]-TT[1][1]*TT[2][2]
			+ 3.*(TT[0][1]*TT[0][1]+TT[0][2]*TT[0][2]+TT[1][2]*TT[1][2]))  +
			1.5*pow(phe,-2.)*(2.*TT[0][0]*dTT_dx[0][0][b][j] + 
			2.*TT[1][1]*dTT_dx[1][1][b][j] +  
			2.*TT[2][2]*dTT_dx[2][2][b][j] +
			-dTT_dx[0][0][b][j]*TT[1][1] - TT[0][0]*dTT_dx[1][1][b][j] +
			-dTT_dx[1][1][b][j]*TT[2][2] - TT[1][1]*dTT_dx[2][2][b][j] +
			-dTT_dx[0][0][b][j]*TT[2][2] - TT[0][0]*dTT_dx[2][2][b][j] +
			6.*TT[0][1]*dTT_dx[0][1][b][j] +  
			6.*TT[1][2]*dTT_dx[1][2][b][j] +  
			6.*TT[0][2]*dTT_dx[0][2][b][j]);   
		    }
	        }
            }

	   if (pd->v[pg->imtrx][var])
	     {
	       for (b=0; b<pd->Num_Species_Eqn; b++) 
	         {
	           for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		     {
		       d_norm2_dc[b][j] = -3. * pow(phe,-3.) * d_phi_dc[b][j]*
		       (TT[0][0]*TT[0][0] + TT[1][1]*TT[1][1] + TT[2][2]*TT[2][2]
		        -TT[0][0]*TT[1][1]-TT[0][0]*TT[2][2]-TT[1][1]*TT[2][2]
		        + 3.*(TT[0][1]*TT[0][1]+TT[0][2]*TT[0][2]+TT[1][2]*TT[1][2]))  +
		        1.5*pow(phe,-2.)*(2.*TT[0][0]*dTT_dc[0][0][b][j] + 
		        2.*TT[1][1]*dTT_dc[1][1][b][j] +  
		        2.*TT[2][2]*dTT_dc[2][2][b][j] +
	       	        -dTT_dc[0][0][b][j]*TT[1][1] - TT[0][0]*dTT_dc[1][1][b][j] +
		        -dTT_dc[1][1][b][j]*TT[2][2] - TT[1][1]*dTT_dc[2][2][b][j] +
		        -dTT_dc[0][0][b][j]*TT[2][2] - TT[0][0]*dTT_dc[2][2][b][j] +
		        6.*TT[0][1]*dTT_dc[0][1][b][j] +  
		        6.*TT[1][2]*dTT_dc[1][2][b][j] +  
		        6.*TT[0][2]*dTT_dc[0][2][b][j]);   
		      }
		  }
	      }

          for (p=0; p<VIM; p++)
            {
              for (q=0; q<VIM; q++)
                {
	          for ( b=0; b<dim; b++)
	            {
	              v = MESH_DISPLACEMENT1 + b;
	              if ( pd->v[pg->imtrx][v] )
	                {
		          dofs     = ei[pg->imtrx]->dof[v];
		          for ( j=0; j<dofs; j++)
		            {
			      d_phi_dTT_dx[p][q][b][j] = -0.75*pow(phe,-2.) * d_phi_dx[b][j] * dT2_dTT[p][q] +
				0.75/phe * dT2_dTT_dx[p][q][b][j];
		            }
		        }
		    }
	        }
            }

	  if (pd->v[pg->imtrx][var])
	    {
              for (p=0; p<VIM; p++)
                {
                  for (q=0; q<VIM; q++)
                    {
	              for (b=0; b<pd->Num_Species_Eqn; b++) 
	                {
	                  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		            {
			      d_phi_dTT_dc[p][q][b][j] = -0.75*pow(phe,-2.) * d_phi_dc[b][j] * dT2_dTT[p][q] +
			      0.75/phe * dT2_dTT_dc[p][q][b][j];
		            }
		        }
		    }
	        }
            }

          for (p=0; p<VIM; p++)
            {
              for (q=0; q<VIM; q++)
                {
	          for ( b=0; b<dim; b++)
	            {
	              v = MESH_DISPLACEMENT1 + b;
	              if ( pd->v[pg->imtrx][v] )
	                {
		          dofs     = ei[pg->imtrx]->dof[v];
		          for ( j=0; j<dofs; j++)
		            {
			      d_phi_dTT_norm_dx[p][q][b][j] = (d_phi_dTT_dx[p][q][b][j] * norm - d_phi_dTT[p][q] * 0.5 * d_norm2_dx[b][j]/norm)/pow(norm,2.);  
		            }
		        }
		    }
	        }
            }

	  if (pd->v[pg->imtrx][var])
	    {
              for (p=0; p<VIM; p++)
                {
                  for (q=0; q<VIM; q++)
                    {
	              for (b=0; b<pd->Num_Species_Eqn; b++) 
	                {
	                  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		            {
			      d_phi_dTT_norm_dc[p][q][b][j] = (d_phi_dTT_dc[p][q][b][j] * norm - d_phi_dTT[p][q] * 0.5 * d_norm2_dc[b][j]/norm)/pow(norm,2.);  
		            }
		        }
		    }
	        }
            }

          for (p=0; p<VIM; p++)
            {
              for (q=0; q<VIM; q++)
                {
	          for ( b=0; b<dim; b++)
	            {
	              v = MESH_DISPLACEMENT1 + b;
	              if ( pd->v[pg->imtrx][v] )
	                {
		          dofs     = ei[pg->imtrx]->dof[v];
		          for ( j=0; j<dofs; j++)
		            {
                              dF_vp_dx[p][q][b][j] = delta_t *  
		              (1./evpl->plastic_mu * d_phi_dx[b][j] * d_phi_dTT_norm[p][q] +
		              1./evpl->plastic_mu * phi * d_phi_dTT_norm_dx[p][q][b][j]); 
		            }
		        }
		    }
	        }
            }

	  if (pd->v[pg->imtrx][var])
	    {
              for (p=0; p<VIM; p++)
                {
                  for (q=0; q<VIM; q++)
                    {
	              for (b=0; b<pd->Num_Species_Eqn; b++) 
	                {
	                  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		            {
                              dF_vp_dc[p][q][b][j] = delta_t *  
		              (-1.*pow(evpl->plastic_mu, -2.)* 
			       d_plastic_mu_dc[b][j] *  phi * d_phi_dTT_norm[p][q] +
			       1./evpl->plastic_mu * 
			      (d_phi_dc[b][j] - d_yield_dc[b][j]) * 
			       d_phi_dTT_norm[p][q] +
		              1./evpl->plastic_mu * phi * d_phi_dTT_norm_dc[p][q][b][j]);
		            }
		        }
		    }
	        }
            }
	}
      else
	{
          for (p=0; p<VIM; p++)
            {
              for (q=0; q<VIM; q++)
                {
	          for ( b=0; b<dim; b++)
	            {
	              v = MESH_DISPLACEMENT1 + b;
	              if ( pd->v[pg->imtrx][v] )
	                {
		          dofs     = ei[pg->imtrx]->dof[v];
		          for ( j=0; j<dofs; j++)
		            {
                              dF_vp_dx[p][q][b][j] = 0.; 
		            }
		        }
		    }
	        }
            }

	  if (pd->v[pg->imtrx][var])
	    {
              for (p=0; p<VIM; p++)
                {
                  for (q=0; q<VIM; q++)
                    {
	              for (b=0; b<pd->Num_Species_Eqn; b++) 
	                {
	                  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		            {
                              dF_vp_dc[p][q][b][j] = 0.;  
		            }
		        }
		    }
	        }
            }
	} /* end of if F_vp_flag loop */
    } /* end of Jacobian Assembly */

  return(F_vp_flag);
} /* end of  get_F_vp */

/*****************************************************************************/
/*  load_elastic_properties                                                  */
/*****************************************************************************/

int 
load_elastic_properties(struct Elastic_Constitutive *elcp,
			double *mu, 
			double *lambda,
 			double *thermexp,
			double speciesexp[MAX_CONC],
			double d_mu_dx[DIM][MDE], 
			double d_lambda_dx[DIM][MDE],
			double d_thermexp_dx[MAX_VARIABLE_TYPES+MAX_CONC],
			double d_speciesexp_dx[MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC])

/*
 *  This function calculates the the elastic properties
 *  for all linear and nonlinear elastic material behavior.
 *  i.e., it is called from both EVP and nonlinear elastic
 *  models as well as linear elastic models.  
 *  The routine supplants evaluations that were done in
 *  mesh_stress, and can eventually be called up above in
 *  matrix_fill
 *
 *  Created by P. R. Schunk, 2/21/99
 */
{
  /*local variables */
  int err;
  int nnn, ins, inode, dim, dofs;
  int a, b, q, v, i, j, w;

  ELASTIC_CONST_STRUCT *elc_ptr;

  memset(d_mu_dx, 0.0, sizeof(double)*DIM*MDE);

  elc_ptr = elcp;  

  double dist, iie, value;

  dim = ei[pg->imtrx]->ielem_dim;

  if(elc_ptr->lame_mu_model == USER )
    {
      err = usr_lame_mu(elc_ptr, elc_ptr->u_mu);
      EH(err,"bad user mu model");
      *mu = elc_ptr->lame_mu;
    }
  else if (elc_ptr->lame_mu_model == CONSTANT )
    {
      *mu = elc_ptr->lame_mu;
    }
  else if (elc_ptr->lame_mu_model == CONTACT_LINE )
    {
      /* make mu (shear modulus) become very large near a critical point,
       * and decay to lower value radially from that point
       * this helps keep elements near critical point from shearing excessively
       * while elements far away won't dilate much
       *  u_mu[0]  is the node set number for the point
       *  u_mu[1]  is minimum value of the shear modulus (comparable to lambda, e.g. 0.5)
       *  u_mu[2]  is the maximum value of the shear modulus (large, e.g. 1e4)
       *  u_mu[3]  is the decay length (problem dependent)
       */

      nnn = -1;
      for (ins = 0; ins < Proc_Num_Node_Sets; ins++) {
	
	/* Check for a match between the ID of the current node set and the node set
	   ID specified in the input file - continue if a match is found */
	
	if (Proc_NS_Ids[ins] ==  (int) elc_ptr->u_mu[0]) nnn = ins;
      }
 
      if( nnn == -1 ) 
       {
         if( Num_Proc == 1 )
             EH(GOMA_ERROR,"Cannot find NS number on Lame MU CONTACT_LINE model");
       }
      else
       {
         if( Proc_NS_Count[ nnn] != 1 ) 
          {
            if( Num_Proc == 1 ) 
                  EH(GOMA_ERROR,"Illegal NS number on Lame MU CONTACT_LINE model");
          }
         else
          {
           inode = Proc_NS_List[Proc_NS_Pointers[ nnn]];

           dist = pow( (  pow(fv->x0[0] - Coor[0][inode],2.0) 
 		        + pow(fv->x0[1] - Coor[1][inode],2.0) ), 0.5) 
                                      / elc_ptr->u_mu[3];
           *mu = elc_ptr->lame_mu = elc_ptr->u_mu[1] + 0.1 / (pow(dist,3.0) 
                                             + 0.1 / elc_ptr->u_mu[2]);
          }
       }
    }
  else if (elc_ptr->lame_mu_model == POWER_LAW )
    {
      /* Shear modulus follows a power law from Scherer
       *  u_mu[0]  is lame_mu (shear modulus) at the reference porosity
       *  u_mu[1]  is the reference porosity
       *  u_mu[2]  is the power law exponent
       */
      *mu = elc_ptr->lame_mu = elc_ptr->u_mu[0] * pow((1. - mp->porosity)/(1. - elc_ptr->u_mu[1]), elc_ptr->u_mu[2]);
      for (i=0; i<MAX_VARIABLE_TYPES + MAX_CONC; i++)
	{
	  elc_ptr->d_lame_mu[i] = - elc_ptr->u_mu[0] * elc_ptr->u_mu[2] *
	    pow((1. - mp->porosity), elc_ptr->u_mu[2] - 1.) /
	    pow((1. - elc_ptr->u_mu[1]), elc_ptr->u_mu[2]) * mp->d_porosity[i];
	}
    }
  else if (elc_ptr->lame_mu_model == SHEAR_HARDEN )
    {
      if(pd->MeshMotion == TOTAL_ALE) EH(GOMA_ERROR,"No TALE Real-solid jacobian entries for SHEAR_HARDEN");
      /* Shear modulus increases when the second invariant of the strain is non-zero
       *  u_mu[0]  is lame_mu (shear modulus) at zero shear
       *  u_mu[1]  is the coefficient of variation
       */
      iie = 0;
      for ( a=0; a<dim; a++)
	{
	  for ( b=0; b<dim; b++)
	    {
	      iie += 0.5 * (fv->strain[a][b] * fv->strain[a][b] - fv->strain[a][a] * fv->strain[b][b]);
	    }
	}
       *mu = elc_ptr->lame_mu = elc_ptr->u_mu[0] + elc_ptr->u_mu[1]*iie*iie;

     for ( a=0; a<dim; a++)
	{
	  for ( b=0; b<dim; b++)
	    {
	      for ( q=0; q<dim; q++)
		{
		  v = MESH_DISPLACEMENT1 + q;
		  if ( pd->v[pg->imtrx][v] )
		    {
		      dofs     = ei[pg->imtrx]->dof[v];
		      for ( j=0; j<dofs; j++)
			{
			  d_mu_dx[q][j] += iie * elc_ptr->u_mu[1]*(2. * fv->strain[a][b] * fv->d_strain_dx[a][b][q][j] 
							       - fv->d_strain_dx[a][a][q][j] * fv->strain[b][b]
							       - fv->strain[a][a] * fv->d_strain_dx[b][b][q][j]);
			}
		    }
		}
	    }
	}
      
    }
  else if (elc_ptr->lame_mu_model == EXPONENTIAL)
    {
    /* Shear modulus decreases exponential with porosity 
     * u_mu[0] = Reference modulus 
     * u_mu[1] = rate of decay
     * u_mu[2] = reference porosity
     */

    *mu = elc_ptr->lame_mu = elc_ptr->u_mu[0]*exp( elc_ptr->u_mu[1]*( elc_ptr->u_mu[2] - mp->porosity ) );
 
    for (i=0; i<MAX_VARIABLE_TYPES + MAX_CONC; i++)
      {
	elc_ptr->d_lame_mu[i] = -  elc_ptr->u_mu[1] * (*mu) *  mp->d_porosity[i];
      } 
    }
  else if (elc_ptr->lame_mu_model == DENSE_POWER_LAW )
    {
      if(pd->MeshMotion == TOTAL_ALE) EH(GOMA_ERROR,"No TALE Real-solid jacobian entries for DENSE_POWER_LAW");
      /* Shear modulus follows a power law 
       *  u_mu[0]  is ultimate lame_mu (shear modulus) when the film is dried
       *  u_mu[1]  is the power law exponent
       */
      /*mu = elc_ptr->lame_mu = elc_ptr->u_mu[0] * pow((1.0000 - fv->volume_change)/(1. - elc_ptr->Strss_fr_sol_vol_frac), elc_ptr->u_mu[1]); */
      *mu =  elc_ptr->u_mu[0] * pow((1.0001 - fv->volume_change)/(1. - elc_ptr->Strss_fr_sol_vol_frac), elc_ptr->u_mu[1]);
      if (*mu < 1.e-12)
        {
	  elc_ptr->lame_mu=1.e-12;
        }
      else
        {
          elc_ptr->lame_mu=*mu;
        }
      for ( q=0; q<dim; q++)
	{
	  v = MESH_DISPLACEMENT1 + q;
	  if ( pd->v[pg->imtrx][v] )
	    {
	      dofs     = ei[pg->imtrx]->dof[v];
	      for ( j=0; j<dofs; j++)
		{
                  if (*mu < 1.e-12)
		    {
		      d_mu_dx[q][j] = 0.;
		    }
		  else
		    {
		  d_mu_dx[q][j] = - elc_ptr->u_mu[0] * elc_ptr->u_mu[1] *
	    	  pow(((1.0000 - fv->volume_change)/(1.-elc_ptr->Strss_fr_sol_vol_frac)), elc_ptr->u_mu[1] - 1.) *
	    	 1./(1.-elc_ptr->Strss_fr_sol_vol_frac) * fv->d_volume_change_dx[q][j]; 
		    }
		}
	    }
	}
    }
  else if (elc_ptr->lame_mu_model == TABLE)
    {
      /* Shear modulus is a table lookup with linear interpolation */
      struct  Data_Table *table_local;
      int var;
      table_local = MP_Tables[elc_ptr->lame_mu_tableid];

      if ( !strcmp(table_local->t_name[0], "FAUX_PLASTIC") ) {
	double nu, MP, slope, var1[1], Enew, Emin;
        double strainI, strainE, strainT, strain;
        double strainIc, strainEc, strainEFF;
        double stressE, stressT, stressEFF;
        double strain_sm, stress_sm, E;
        double d_Enew, d_Enew_dc;
        double d_strainI, d_strainE, d_strainT, d_strain;
        double d_strainIc_dc, d_strainEc_dc, d_strainEFF, d_strainEFF_dc;
        double d_stressE, d_stressT, d_stressEFF, d_stressEFF_dc;
        
        // Internal parameters
        int use_new_behavior = 1;
        int lag_strain = 0;
        int true_stress_strain = 1;
        double small_strain = 1e-6;
        double compare_factor = 0.99999;
        double downslope_scale = 0.1;

	// Check for Poisson ratio
	if ( elc_ptr->lame_lambda_model != POISSON_RATIO ) {
	  EH(GOMA_ERROR, "Must specify POISSON_RATIO for Lame Lambda in order to use FAUX_PLASTIC");
	}
        
        // If Emin is not specified, don't use the new behavior
        Emin = table_local->Emin;
        if ( Emin < 0.0 ) use_new_behavior = 0;

	// Check that the max strain equation is active
	var = MAX_STRAIN;
	if ( !pd->v[pg->imtrx][var] ) {
	  EH(GOMA_ERROR, "The max_strain equation must be activated to use FAUX_PLASTIC");	  
	}
        
        // If using the new behavior, then the cur_strain equation must be on.
	var = CUR_STRAIN;
        if ( use_new_behavior ) {
          if ( !pd->v[pg->imtrx][var] ) {
            EH(GOMA_ERROR, "The cur_strain equation must be activated to use 'new' FAUX_PLASTIC");	  
          }
        }
	
	// Calculate poisson's ratio and multiplier, mu=MP*E
	nu = elc_ptr->u_lambda[0];
	MP = 1.0/(2*(1+nu));

	// Calculate MAXIMUM von Mises strain and engineering strain
        if( lag_strain ) {
          strainI = fv_old->max_strain;
          d_strainI = 0.0;
        } else {
          strainI = fv->max_strain;
          d_strainI = 1.0;
        }
	strainE = 3*strainI/(2.0*(1+nu));
        d_strainE = 3.0/(2.0*(1+nu))*d_strainI;
        strainT = log(1+strainE);
        d_strainT = 1.0/(1+strainE)*d_strainE;

	// Calculate CURRENT (old) von Mises strain and engineering strain 
        if( lag_strain ) {
          strainIc = fv_old->cur_strain;
          d_strainIc_dc = 0.0;
        } else {
          strainIc = fv->cur_strain;
          d_strainIc_dc = 1.0;
        }
	strainEc = 3*strainIc/(2.0*(1+nu));
        d_strainEc_dc = 3.0/(2.0*(1+nu))*d_strainIc_dc;
        
        // Select value of strain for calculations
        if (true_stress_strain) {
          strain = strainT;	
          d_strain = d_strainT;
        } else {
          strain = strainE;
          d_strain = d_strainE;
        }
        
        // Modulate strain for small values
        if ( strainE  < small_strain ) {
          strainE  = small_strain;
          d_strainE = 0.0;
        }
        if ( strainEc < small_strain ) {
          strainEc = small_strain;
          d_strainEc_dc = 0.0;
        }
        if ( strain   < small_strain ) {
          strain   = small_strain;
          d_strain = 0.0;
        }

	// Calculate engineering stress from table, based on format
        var1[0] = strain;
        if (true_stress_strain) {
          stressT = interpolate_table( table_local, var1, &slope, NULL);
          d_stressT = slope*d_strain;
          stressE = stressT/(1.0+strainE);
          d_stressE = 1.0/(1+strainE)*d_stressT - stressT/pow(1+strainE,2)*d_strainE;
        } else {
          stressE = interpolate_table( table_local, var1, &slope, NULL);
          d_stressE = slope*d_strain;
        }
        
        // Decide if we are on the downslope, follow a different path
        if ( (strainEc < compare_factor*strainE) && use_new_behavior ) {
          
          // Calculate elastic modulus of the small-strain elastic regime
          strain_sm = small_strain;
          stress_sm = interpolate_table( table_local, &strain_sm, &slope, NULL);
          E = stress_sm/strain_sm*table_local->yscale*downslope_scale;
          
          // Calculate stress from an elastic response from max strain
          stressEFF = E*(strainEc-strainE) + stressE;
          d_stressEFF = -E*d_strainE + d_stressE;
          d_stressEFF_dc = E*d_strainEc_dc;
          
          // The effective strain is the current engineering strain
          strainEFF = strainEc;
          d_strainEFF = 0.0;
          d_strainEFF_dc = d_strainEc_dc;
          
        } else {
          
          // Use the max engineering stress
          stressEFF = stressE;
          d_stressEFF = d_stressE;
          d_stressEFF_dc = 0.0;
          
          // Use the max engineering strain
          strainEFF = strainE;
          d_strainEFF = d_strainE;
          d_strainEFF_dc = 0.0;
          
        }
        
        // Scale the engineering stress based on the input scaling parameter
        stressEFF *= table_local->yscale;
        d_stressEFF *= table_local->yscale;
        d_stressEFF_dc *= table_local->yscale;
        
        // Calculate modulus and compare to minimum value
        Enew = stressEFF/strainEFF;
        d_Enew = 1.0/strainEFF*d_stressEFF - stressEFF/pow(strainEFF,2)*d_strainEFF;
        d_Enew_dc = 1.0/strainEFF*d_stressEFF_dc - stressEFF/pow(strainEFF,2)*d_strainEFF_dc;
        if ( Enew < Emin ) {
          Enew = Emin;
          d_Enew = 0.0;
          d_Enew_dc = 0.0;
        }
        
	// Store final Lame MU
	*mu = elc_ptr->lame_mu = MP*Enew;
        elc_ptr->d_lame_mu[MAX_STRAIN] = MP*d_Enew;
        elc_ptr->d_lame_mu[CUR_STRAIN] = MP*d_Enew_dc;

      } else {
	apply_table_mp(&elc_ptr->lame_mu, table_local);
	*mu = elc_ptr->lame_mu;
	
	for ( i = 0; i < table_local->columns-1; i++) {
	  var = table_local->t_index[i];
	  elc_ptr->d_lame_mu[var] = table_local->slope[i];
	}
      }
    }
  else
    {
      EH(GOMA_ERROR,"Unrecognized mu model");
    }

  /*If needed, hit lame_me with temperature shift */
  if (elc_ptr->lameTempShiftModel == CONSTANT)
    {
      *mu *= elc_ptr->lame_TempShift;
      elc_ptr->lame_mu *= elc_ptr->lame_TempShift;
      elc_ptr->d_lame_mu[TEMPERATURE] = 0.;
    }
  else if (elc_ptr->lameTempShiftModel == POWER_LAW)
    {
      double Troom = elc_ptr->u_lame_TempShift[0];
      double Tmelt = elc_ptr->u_lame_TempShift[1];
      double exponent = elc_ptr->u_lame_TempShift[2];
      if(fv->T > Troom)
	{		      
	  *mu *= 1.-pow(((fv->T - Troom)/(Tmelt-Troom)), exponent);
          for (i=0; i<MAX_VARIABLE_TYPES + MAX_CONC; i++) {
            elc_ptr->d_lame_lambda[i] *= 1.-pow(((fv->T - Troom)/(Tmelt-Troom)), exponent);
          }
	  elc_ptr->d_lame_mu[TEMPERATURE] = -exponent*pow(((fv->T - Troom)/(Tmelt-Troom)), exponent - 1.0)*elc_ptr->lame_mu/(Tmelt-Troom);
	}
      else
	{
	  *mu *=1.;
	  elc_ptr->d_lame_mu[TEMPERATURE] = 0.;
	}

      elc_ptr->lame_mu = *mu;

    }
  else if (elc_ptr->lameTempShiftModel == TABLE)
    {
      struct  Data_Table *table_local;
      int var;
      table_local = MP_Tables[elc_ptr->lame_TempShift_tableid];
      apply_table_mp(&elc_ptr->lame_TempShift, table_local);

      *mu *= elc_ptr->lame_TempShift;
      elc_ptr->lame_mu *= elc_ptr->lame_TempShift;
      
      for ( i = 0; i < table_local->columns-1; i++) {
	var = table_local->t_index[i];
	elc_ptr->d_lame_mu[var] = elc_ptr->lame_mu/elc_ptr->lame_TempShift*table_local->slope[i];
      }
    }

  if(elc_ptr->lame_lambda_model == USER )
    {
      err = usr_lame_lambda(elc_ptr, elc_ptr->u_lambda);
      EH(err,"bad user lambda model");
      *lambda = elc_ptr->lame_lambda;
    }
  else if (elc_ptr->lame_lambda_model == CONSTANT )
    {
      *lambda = elc_ptr->lame_lambda;
    }
  else if (elc_ptr->lame_lambda_model == POWER_LAW )
    {
      /* Shear modulus follows a power law from Scherer
       *  u_lambda[0]  is lame_lambda at the reference porosity
       *  u_lambda[1]  is the reference porosity
       *  u_lambda[2]  is the power law exponent
       */
      *lambda = elc_ptr->lame_lambda = elc_ptr->u_lambda[0] * pow((1. - mp->porosity)/
							 (1. - elc_ptr->u_lambda[1]), elc_ptr->u_lambda[2]);
      for (i=0; i<MAX_VARIABLE_TYPES + MAX_CONC; i++)
	{
	  elc_ptr->d_lame_lambda[i] = - elc_ptr->u_lambda[0] * elc_ptr->u_lambda[2] *
	    pow((1. - mp->porosity), elc_ptr->u_lambda[2] - 1.) /
	    pow((1. - elc_ptr->u_lambda[1]), elc_ptr->u_lambda[2]) * mp->d_porosity[i];
	}
    }
  else if (elc_ptr->lame_lambda_model == EXPONENTIAL)
    {
    /* Shear modulus decreases exponential with porosity 
     * u_lambda[0] = Reference modulus 
     * u_lambda[1] = rate of decay
     * u_lambda[2] = reference porosity
     */

      *lambda = elc_ptr->lame_lambda = elc_ptr->u_lambda[0]*exp( elc_ptr->u_lambda[1]*( elc_ptr->u_lambda[2] - mp->porosity ) );
 
      for (i=0; i<MAX_VARIABLE_TYPES + MAX_CONC; i++)
	{
	  elc_ptr->d_lame_lambda[i] = -elc_ptr->u_lambda[1] * (*lambda) *  mp->d_porosity[i];
	} 
    }
  else if (elc_ptr->lame_lambda_model == POISSON_RATIO)
    {
      /* 
	 Calculate lambda using mu and poisson's ratio using the formula:
	 lambda = 2*mu*nu/(1-2*nu)
      */
      dbl nu = elc_ptr->u_lambda[0];

      *lambda = elc_ptr->lame_lambda = 2*elc_ptr->lame_mu*nu/(1.0-2*nu);
      
      for (i=0; i<MAX_VARIABLE_TYPES + MAX_CONC; i++) {
	elc_ptr->d_lame_lambda[i] = 2*elc_ptr->d_lame_mu[i]*nu/(1.0-2*nu);
      }

      // For FAUX_PLASTIC model, we need to incorporate mesh sensitivities
      if ( elc_ptr->lame_mu_model == TABLE ) {
	for ( i = 0; i < DIM; i++) {
	  for ( j = 0; j < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; j++) {
	    d_lambda_dx[i][j] = d_mu_dx[i][j] * 2*nu/(1.0-2*nu);
	  }
	}
      }

    }
  else
    {
      EH(GOMA_ERROR,"Unrecognized lambda model");
    }

/*  thermal expansion	*/
   if(elc_ptr->thermal_expansion_model == CONSTANT )
     {
       *thermexp = elc_ptr->thermal_expansion;
     }
   else if(elc_ptr->thermal_expansion_model == SHRINKAGE)
     {
	*thermexp = elc_ptr->u_thermal_expansion[0];
     }
   else if(elc_ptr->thermal_expansion_model == IDEAL_GAS)
     {
	*thermexp = elc_ptr->u_thermal_expansion[0];
     }
   else if(elc_ptr->thermal_expansion_model == USER )
     {
      if(pd->MeshMotion == TOTAL_ALE) EH(GOMA_ERROR,"No TALE Real-solid jacobian entries for USER");
      err = usr_expansion(elc_ptr->u_thermal_expansion, &value, d_thermexp_dx);
	*thermexp = value;
     }
   else
     {
       EH(GOMA_ERROR,"Unrecognized thermal expansion model");
     }

/*  species expansion	*/
   if (pd->e[pg->imtrx][R_MASS] && 
        (pd->MeshMotion != ARBITRARY || mp->SpecVolExpModel[0] == PHOTO_CURING)) 
   {
	for(w=0 ; w<pd->Num_Species_Eqn ; w++)
	   {
   	    if(mp->SpecVolExpModel[w] == CONSTANT  ||
                mp->SpecVolExpModel[w] == PHOTO_CURING)
     		{
       	    	speciesexp[w] = mp->species_vol_expansion[w];
     		}
   	    else
     		{
       		EH(GOMA_ERROR,"Unrecognized species expansion model");
     		}
	   }
   }
  return(1);
} /*End of load_elastic_properties*/

/*****************************************************************************/
/*  load_plastic_properties                                                  */
/*****************************************************************************/

int 
load_plastic_properties(double d_plastic_mu_dc[MAX_CONC][MDE], 
			double d_yield_dc[MAX_CONC][MDE])

/*
 *  This function calculates the the plastic properties
 *  that are used in EVP. The properties are functions of
 *  solvent concentration
 *
 *  Created by S.Y Tam 5/30/2000
 */
{
  /*local variables */
  int err = 0;
  double min;
  int b, j, var;

  if(evpl->plastic_mu_model == USER )
    {
/*  Take care of this later
      err = usr_plastic_mu(evpl->u_plastic_mu); */
      EH(err,"bad user plastic_mu model");
    }
  else if (evpl->plastic_mu_model == LINEAR )
    {
      /*  Plastic viscosity increases linearly from the smaller value to
       *  the larger value with decrease in solvent concentration. 
       *  The 2 terminal values can be entered in random order. 
       *  u_mu[0]  is one of the two terminal plastic viscosities 
       *  u_mu[1]  is the other terminal plastic viscosity
       *  fv->c[0] is the solvent concentration for binary system
       *  if more than 1 solvent is present, then fv->c[0] should be replaced.
       */
        if (evpl->u_plastic_mu[1] > evpl->u_plastic_mu[0])
        {
          min= evpl->u_plastic_mu[0];
        }
      else
        {
          min= evpl->u_plastic_mu[1];
	}   
      evpl->plastic_mu = min + (elc->Strss_fr_sol_vol_frac - fv->c[0])/elc->Strss_fr_sol_vol_frac *
		   fabs(evpl->u_plastic_mu[1] - evpl->u_plastic_mu[0]); 
      var=MASS_FRACTION;
      for (b=0; b<pd->Num_Species_Eqn; b++) 
        {
          for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
            {
	  d_plastic_mu_dc[b][j] = -fv->c[b]/elc->Strss_fr_sol_vol_frac *
		   fabs(evpl->u_plastic_mu[1] - evpl->u_plastic_mu[0]);
            }
	} 
    }
  else
    {
      EH(GOMA_ERROR,"Unrecognized plastic_mu model");
    }

  if(evpl->yield_model == USER )
    {
/*   take care of this later 
      err = usr_plastic_mu(evpl->u_yield); */
      EH(err,"bad user yield model");
    }
  else if (evpl->yield_model == LINEAR )
    {
      /*  Yield stress increases linearly from the smaller value to 
       *  the larger value with decrease in solvent concentration. 
       *  The 2 terminal values can be entered in random order. 
       *  u_mu[0]  is one of the two terminal yield stresses 
       *  u_mu[1]  is the other terminal yield stress
       *  fv->c[0] is the solvent concentration for binary system
       *  if more than 1 solvent is present, then fv->c[0] should be replaced.
       */
      if (evpl->u_yield[1] > evpl->u_yield[0])
        {
          min= evpl->u_yield[0];
        }
      else
        {
          min= evpl->u_yield[1];
	} 
      evpl->yield = min + (elc->Strss_fr_sol_vol_frac - fv->c[0])/elc->Strss_fr_sol_vol_frac *
		   fabs(evpl->u_yield[1] - evpl->u_yield[0]); 
      var=MASS_FRACTION;
      for (b=0; b<pd->Num_Species_Eqn; b++) 
        {
          for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
            {
	  d_yield_dc[b][j] = - fv->c[b]/elc->Strss_fr_sol_vol_frac *
		   fabs(evpl->u_yield[1] - evpl->u_yield[0]);
            }
	} 
    }
  else
    {
      EH(GOMA_ERROR,"Unrecognized yield model");
    }

  return(1);
} /*End of load_plastic_properties*/


int check_for_neg_elem_volume( int eb,
                               double x[],
                               double resid_vector[],
                               double x_old[],
                               double x_older[],
                               double xdot[],
                               double xdot_old[],
                               double x_update[],
                               double *ptr_delta_t,
                               double *ptr_theta,
                               double *ptr_time_value,
                               Exo_DB *exo )
{
  /* check all elems of the given element block for a neg_volume */
  int e_start, e_end, e;
  int ib;
  double xi[3];
  int err, ip_total, ip;

  /* Locate block in exo structure */
  for (ib = 0; ib < exo->num_elem_blocks; ib++)
    {
      if (exo->eb_id[ib] == eb) break;
    }
  if ( exo->eb_id[ib] != eb ) EH(GOMA_ERROR,"Couldn't locate element block.");
    
  e_start = exo->eb_ptr[ib];
  e_end   = exo->eb_ptr[ib+1];

  /* now populate augc with elem number */
  for (e = e_start; e < e_end && !neg_elem_volume; e++)
    {
      err = load_elem_dofptr(e, exo, x, x_old, xdot, xdot_old, 0);

      ip_total = elem_info(NQUAD, ei[pg->imtrx]->ielem_type);

      /* Loop over all the Volume Quadrature integration points */

      for (ip = 0; ip < ip_total; ip++)
        {

          find_stu(ip, ei[pg->imtrx]->ielem_type, xi, xi+1, xi+2); /* find quadrature point */
          fv->wt = Gq_weight (ip, ei[pg->imtrx]->ielem_type); /* find quadrature weights for */

          err = load_basis_functions(xi, bfd);
          EH( err, "problem from load_basis_functions");

          err = beer_belly();
          EH( err, "beer_belly");

          err = load_fv();
          EH( err, "load_fv");

          err = load_bf_grad();
          EH( err, "load_bf_grad");

          err = load_bf_mesh_derivs();
          EH( err, "load_bf_mesh_derivs");

          err = load_fv_grads();
          EH( err, "load_fv_grads");

          err = load_fv_mesh_derivs(1);
          EH( err, "load_fv_mesh_derivs");
    
          err = belly_flop(elc->lame_mu);
          EH(err, "error in belly flop");
          
          if (neg_elem_volume) break;
        }
    }                      
  return neg_elem_volume;
}
                                        
/*****************************************************************************/
/* end of file mm_fill_solid.c */
/*****************************************************************************/
