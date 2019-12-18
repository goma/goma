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

/* mm_fill_species -- auxiliary routines helpful in matrix & residual assembly 
 * for species equations and boundary conditions
 */

/*
 *$Id: mm_fill_species.c,v 5.15 2009-03-13 01:39:03 hkmoffa Exp $
 */

/* Standard include files */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

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
#include "rf_fill_const.h"
#include "rf_vars_const.h"
#include "mm_mp_const.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"

#include "mm_eh.h"

#ifdef USE_CHEMKIN
#include "ck_chemkin_const.h"
#endif

#include "mm_mp_structs.h"
#include "mm_mp_const.h"
#include "mm_mp.h"

#include "mm_fill_terms.h"
#include "mm_fill_population.h"

#define GOMA_MM_FILL_SPECIES_C
#include "mm_fill_species.h"
#include "sl_aux.h"
#include "goma.h"

/*********** R O U T I N E S   I N   T H I S   F I L E *************************
*
*						-All routines in this
*				NOTE:		 file are called only by
*						 mm_fill.c: matrix_fill
*						 except possibly for static
*						 functions.
*
*       NAME			TYPE			CALLED BY
*  -----------------------------------------------------------------
*  assemble_mass_transport       int  
*  mass_flux_surf_mtc		 void
*  mass_flux_surf_BV2            void  (RSL 6/22/02)
*  mass_flux_surf_NI             void  (RSL 6/22/02)
*  mass_flux_surf_const		 void  
*  raoults_law                   void
*  flory_huggins                 void
*  kinematic_species_bc          void
*  mass_flux_surf_bc		 void
*  const_mass_flux_surf_bc	 void  
*  mass_flux_BV2_surf_bc         void  (RSL 6/22/02)
*  mass_flux_NI_surf_bc          void  (RSL 6/22/02)
*  current_BV2_surf_bc           void  (RSL 6/22/02)
*  current_NI_surf_bc            void  (RSL 6/22/02)
*  mass_flux_equil_mtc		 void
*  mtc_chilton_coburn		 void  (RBS 4/10/03)
*  act_coeff                     void
*  get_equil_surf_bc             void
*  sus_mass_flux_surf_bc         void
*  kin_bc_leak                   void
*  kin_bc_electrodeposition      void  (RSL 5/27/02)
*  vnorm_bc_electrodeposition    void  (RSL 5/30/02)
*  lat_heat_bc                   void
*  lat_heat_internal_bc          void
*  mass_flux_surf		 void
*  get_convection_velocity       int
*  get_continuous_species_terms  int
*  ludcmp                        void
*  lubksb                        void
*  Stefan_Maxwell_diff_flux      int
*  fickian_flux                  int
*  fickian_charged_flux          int
*  fickian_charged_flux_x        int  (RSL 9/18/00)
*  generalized_fickian_flux      int
*  assemble_invariant            int
*  get_particle_convection_velocity   int       get_continuous_species_terms.
*
*******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/* assemble_mass_transport() -- assemble terms (Residual &| Jacobian) for mass
 *                              conservation eqns
 *
 * in:
 * 	ei -- pointer to Element Indices	structure
 *	pd -- pointer to Problem Description	structure
 *	af -- pointer to Action Flag		structure
 *	bf -- pointer to Basis Function		structure
 *	fv -- pointer to Field Variable		structure
 *  fv_old -- pointer to old Field Variable	structure
 *  fv_dot -- pointer to dot Field Variable	structure
 *	cr -- pointer to Constitutive Relation	structure
 *	md -- pointer to Mesh Derivative	structure
 *	me -- pointer to Material Entity	structure
 *
 * out:
 *	a   -- gets loaded up with proper contribution 
 *	lec -- gets loaded up with local contributions to resid, Jacobian
 * 	r   -- residual RHS vector
 *
 * Created:	Thu Mar  3 07:48:01 MST 1994 pasacki@sandia.gov
 *
 * Revised:	9/24/94 by RRR
 * Revised:	10/21/94 by RAC (adapted energy -> species)
 * Revised:     8/98 by KSC to enable modeling multicomponent transport of
 *                      neutral or charged species  via Stefan_Maxwell flux
 *                      equations
 * Revised:	4/00 by RSL to allow for sensitivity of diffusive flux to
 *         	voltage
 * Revised:	6/00 by RSL to implement constant total flux boundary condition
 * Revised:	7/00 by RSL to generalize species conservation equation for a
 *         	porous electrode
 * Revised:	8/00 by RSL to supply Jacobian terms needed in modeling electrode
 *         	processes
 * Revised:     9/00 by KSC to implement the fickian_charged routine and 
 *              the J_s_V sensitivity for modeling transport of charged species
 *              in dilute electrolyte solutions. 
 * Revised:     11/00 by KSC to implement Butler-Volmer kinetics routine for
 *              handling surface electrochemical reactions and to modify the 
 *              kin_bc_leak routine to account for the effects of molar-volume
 *              difference between a liquid reactant and a solid product on 
 *              moving boundary velocity (as in LIGA electrodeposition).  
 * Revised:     12/2001 by KSC (to clean up call to the electrode_species_source routine) 
 *
 * Note: currently we do a "double load" into the addresses in the global
 *       "a" matrix and resid vector held in "esp", as well as into the
 *       local accumulators in "lec".
 *
 * MMH
 * This required extensive modification for the SUSPENSION_PM model.  For one
 * of the species, the velocities are not the normal Navier-Stokes velocities,
 * but the particle phase velocities...  For SUSPENSION_PM, I am assuming:
 *
 *   o CONTINUOUS PorousMediaType.
 *   o FICKIAN MassFluxModel (although it shouldn't matter b/c of
 *     zero diffusion).
 *   o Only mass and advection terms are turned on.
 *   o All weight functions are GALERKIN (supg=0.0).
 *
 * HKM Note about the sign convention:
 *
 *    A note about the sign convention of all of the volumetric and
 *    surface integral terms in these equations. The time derivative
 *    in the residual equation enters with a NEGATIVE sign. Thus, the
 *    term in the residual is
 *
 *               vol_integral( - dc/dt * phi_i)
 *
 *    All terms, volume and surface integral, must reflect this
 *    sign convention. In particular, surface integrals added into
 *    the residual must equate to the negative of the outward facing
 *    flux.
 */

int 
assemble_mass_transport(double time, /* present time valuel; KSC             */
			double tt, /* parameter to vary time integration from 
				    * explicit (tt = 1) to implicit (tt = 0) */
			double dt, /* current time step size */
                        PG_DATA *pg_data)
{
  int var, ii,  pvar, ledof;
  const int eqn = R_MASS;
  const int dim = pd->Num_Dim;
  int p, b, q, w, w1, i, j, status, c;
  /*
   *    species_eqn_type:
   *        This is a temp variable that gets set at the top
   *        from Species_Var_Type.
   *        Specifies what identity of the independent variable
   *        The form of the equation and the units for the equation
   *        are dependent on the this determination.
   */
  int species_eqn_type;
  /*
   *    coeff_rho:
   *        If coeff_rho_nonunity is TRUE, it means that the time
   *        derivative and advection term have variable density or
   *        concentration terms in them that are not part of the
   *        independent variable itself.
   *
   *              Species_Var_type   Time_Deriv_Term  coeff_rho_nonunity
   *             ---------------------------------------------------------
   *             SPECIES_MASS_FRACTION     rho * d(Y_k)/dt       TRUE
   *             SPECIES_MOLE_FRACTION     conc * d(X_k)/dt      TRUE
   *             SPECIES_CONCENTRATION     d(conc_k)/dt         FALSE
   *                DEFAULT                d(????)/dt           FALSE
   *           ---------------------------------------------------------
   *
   *        coeff_rho is equal to one for coeff_rho_nonunity FALSE
   *        cases. It is non-unity and variable (with units) for
   *        for coeff_rho_nonunity TRUE cases.
   */
  int coeff_rho_nonunity = FALSE;
  dbl coeff_rho;
  dbl rho;				/* Density. */
  dbl sumxm = 0.0, sumrm, epsilon = 0.0, small_c = 0.0, sumxdot, sumxdotm;  
  dbl sumdxdot, sumdxdotm, sumxms;  
  dbl psum, sumflux, sumfluxm;  
  dbl psum1, psum2, psum3, psum4, sumdflux, sumdfluxm, dmdx;  
  dbl sumdrm, factor;  
  dbl sumx;
  int num_species;  
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;  /* density dependence */
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;

  dbl M[MAX_CONC];            /* species molecular weight */ 
  dbl x[MAX_CONC];            /* mole fraction */

  struct Species_Conservation_Terms s_terms;
  memset(&s_terms, 0, sizeof(struct Species_Conservation_Terms));

  dbl mass;		         	/* For terms and their derivatives */

  dbl advection;			/* For terms and their derivatives */
  dbl advection_a, advection_b, advection_c, advection_d, advection_e,
      advection_f;
  dbl diffusion;
  dbl diff_a, diff_b, diff_c, source;

  /*
   * Galerkin weighting functions for i-th energy residuals 
   * and some of their derivatives...
   */

  dbl phi_i;
  dbl grad_phi_i[DIM];

  /*
   * Petrov-Galerkin weighting functions for i-th residuals 
   * and some of their derivatives...
   */

  dbl wt_func = 0;

  /* SUPG variables */
  dbl supg=0, d_wt_func = 0;

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;

  dbl h3;			/* Product of all 3 scale factors used for */
				/* volume integrals in this coordinate system*/
  dbl dh3dmesh_bj;		/* Mesh derivative of same. */

  dbl det_J;

  dbl d_det_J_dmeshbj;			/* for specified (b,j) mesh dof */
  dbl dgrad_phi_i_dmesh[DIM];		/* ditto.  */
  dbl wt;
  int err;
  int v_g[DIM][DIM];
  int taylor_galerkin[MAX_CONC];
  int species = -1;             /* Which "w" is the particle phase ... */
  status = 0;

  /*initialize grad_phi_i */
  for (i=0; i<DIM; i++) {
    grad_phi_i[i] = 0;
  }

  /*
   * Unpack variables from structures for local convenience...
   */



  /*
   * Bail out fast if there's nothing to do...
   */

  if ( ! pd->e[pg->imtrx][eqn] )
    {
      return(status);
    }

  for ( w=0; w<pd->Num_Species_Eqn; w++)
    {
      taylor_galerkin[w] =mp->SpeciesTimeIntegration[w]; 
    }

  wt = fv->wt;				/* Gauss point weight. */

  det_J = bf[eqn]->detJ;		/* Really, ought to be mesh eqn. */

  h3 = fv->h3;			        /* Differential volume element. */

  num_species = pd->Num_Species;  

  for (j=0; j<num_species; j++)  
     {
      M[j] = mp->molecular_weight[j];
     }

  sumx = 0.;
  for (j=0; j<num_species-1; j++) 
     {
      x[j] = fv->c[j];
      sumx += x[j];
     }
  x[num_species-1] = 1. - sumx;

  /* MMH
   * For species transport in the SUSPENSION_PM model, there is one 
   * species that represents the particles.  There may be other species
   * involved, though, and they care about the "averaged" density, which
   * is the fluid density here.
   *
   * I don't see where this is used any where.  The functions that
   * this function calls seem to be able to figure out the correct
   * density on their own...
   */
  rho  = density(d_rho, time);


  SUPG_terms supg_terms;

  if( mp->Spwt_funcModel == GALERKIN)
    {
      supg = 0.;
    }
  else if( mp->Spwt_funcModel == SUPG)
    {
      supg = mp->Spwt_func;
    }

  /************************************************************************/
  /*                       START OF SPECIES ASSEMBLE                      */
  /************************************************************************/
  /*
   *  Initialize the Species_Conservation_Terms temporary structure
   *  before filling it up
   */
  zero_structure(&s_terms, sizeof(struct Species_Conservation_Terms), 1);

/*  if (mp->PorousMediaType == CONTINUOUS)
    { */
      /*
       * ---------------------------------------------------------------------
       * ----- Calculate the capacity, flux, convection, and source terms for
       *       a continuous medium 
       * ---------------------------------------------------------------------
       */
      err = get_continuous_species_terms(&s_terms, time, tt, dt, pg_data->hsquared);
      EH(err,"problem in getting the species terms");
      
/*    } */   /* end of if CONTINUOUS */

  for ( w=0; w<pd->Num_Species_Eqn; w++)
    {

      if (supg != 0.) {
        dbl D = 1e-6;
        if (mp->DiffusivityModel[w] == CONSTANT) {
          D = mp->diffusivity[w];
        }
        get_supg_tau(&supg_terms, dim, D, pg_data);
      }

      double Heaviside = 1;

      if (ls != NULL) {
        if (mp->SpeciesOnlyDiffusion[w] == DIFF_POSITIVE) {
          load_lsi(ls->Length_Scale);
          Heaviside = 1 - lsi->H;
        } else if (mp->SpeciesOnlyDiffusion[w] == DIFF_NEGATIVE) {
          load_lsi(ls->Length_Scale);
          Heaviside = lsi->H;
        }
      }


      /*
       * Residuals_________________________________________________________________
       */
      if ( af->Assemble_Residual )
        {
          var = MASS_FRACTION;
          /*
           *  Store the species eqn type (which is keyed to the Variable
           *  type in a temporary variable).
           */
          species_eqn_type = mp->Species_Var_Type;

          /*
           *   START loop over species equations. The outer loop is over
           *   the species number
           */

	  /*
	   *  Calculate the coef_rho term based upon the value of
	   *  species_eqn_type.
	   *     This coefficient indicates whether the time and advection
	   *     term should be multipled by the density or maybe the
	   *     concentration. 
	   */
	  coeff_rho = 1.0;
	  coeff_rho_nonunity = FALSE;
	  if ((species_eqn_type == SPECIES_MASS_FRACTION) ||
	      (species_eqn_type == SPECIES_MOLE_FRACTION)) {
	    coeff_rho = rho;
	    coeff_rho_nonunity = TRUE;
	  }

          if (mp->SpeciesSourceModel[w]  == ELECTRODE_KINETICS ||
              mp->SpeciesSourceModel[w]  == ION_REACTIONS) /*  RSL 3/19/01  */
	     {
	      sumxm = 0.;
	      for (j=0; j<num_species; j++)
	         {
	          sumxm += x[j] * M[j];
	         }
	      if (mp->PorosityModel == CONSTANT)
	         {
	          epsilon = mp->porosity;
	         }
	      else if (mp->PorosityModel == THERMAL_BATTERY)
	         {
	          epsilon = mp->u_porosity[0];
	         }
	      else
	         {
	          EH(-1, "invalid porosity model");
	         }
	      /*  rho = density(d_rho); /\*  RSL 6/22/02  *\/ */ /* not sure why we are calling density again, so I commented this out -RRR*/
	      small_c = rho/sumxm;
              coeff_rho = small_c; /*  RSL 9/27/01  */
	     }
		  
	  /*
	   *  Loop over the degrees of freedom in the element corresponding
	   *  to species unknowns. These are at different nodes in the
	   *  element, usually. However, there can be more than one set
	   *  of degrees of freedom at each node. Note, this
	   *  step doesn't depend upon the species equation number, so
	   *  we might think about exchanging the order of the loops!
	   */
	  for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
	    ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];
	    if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
	      /*
	       *  Here is where we figure out whether the row is to placed in
	       *  the normal spot (e.g., ii = i), or whether a boundary condition
	       *  require that the volumetric contribution be stuck in another
	       *  ldof pertaining to the same variable type.
	       */
	      ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

		  phi_i = bf[eqn]->phi[i];

		  /* only use Petrov Galerkin on advective term - if required */
		  wt_func = phi_i;
		  /* add Petrov-Galerkin terms as necessary */
		  if (supg != 0.0)
		    {
		      for(p=0; p<dim; p++)
			{
                          wt_func += supg * supg_terms.supg_tau * fv->v[p] * bf[eqn]->grad_phi[i][p];
			}
		    }

		  mass = 0.;
		  if ( pd->TimeIntegration != STEADY )
		    {
		      if ( pd->e[pg->imtrx][eqn] & T_MASS ) {
			  mass = coeff_rho * s_terms.Y_dot[w];

                          if (mp->SpeciesSourceModel[w]  == ELECTRODE_KINETICS ||
                              mp->SpeciesSourceModel[w]  == ION_REACTIONS) /*  RSL 3/19/01  */
			  {
			    mass = s_terms.Y_dot[w];
			    sumxdot = 0.0;
			    sumxdotm = 0.0;
			    for (j = 0; j < num_species - 1; j++) {
			      sumxdotm += s_terms.Y_dot[j] * M[j];
			      sumxdot  += s_terms.Y_dot[j];
			    }
			    sumxdotm -= sumxdot * M[num_species-1];
			    mass -= x[w]*sumxdotm/sumxm;
			    mass *= (epsilon*small_c);
			  }

			  mass *= - wt_func * h3 * det_J * wt;
			  mass *= pd->etm[pg->imtrx][eqn][LOG2_MASS];
			}
		    }
		  

		  
		  /*
		   *   Advection is velocity times gradient of the species unknown
		   *   variable, possibly multiplied by a density or total
		   *   concentration, depending upon the species variable
		   *   type.
		   */
		  advection = 0.0;
		  if ( pd->e[pg->imtrx][eqn] & T_ADVECTION )
		    {
		      advection_a = 0.0;
		      for ( p=0; p<VIM; p++)
			{
			  advection_a += s_terms.conv_flux[w][p];
			}

                      if (mp->SpeciesSourceModel[w]  == ELECTRODE_KINETICS ||
                          mp->SpeciesSourceModel[w]  == ION_REACTIONS) /*  RSL 3/19/01  */
		         {
		          psum = 0.;
		          for ( p=0; p<VIM; p++)
		             {
		              sumflux = 0.;
		              sumfluxm = 0.;
		              for (j=0; j<num_species-1; j++)
		                 {
		                  sumfluxm += s_terms.conv_flux[j][p] * M[j];
		                  sumflux  += s_terms.conv_flux[j][p];
		                 }
		              sumfluxm -= sumflux * M[num_species-1];
		              psum += sumfluxm;
		             }
		          advection_a -= x[w]*psum/sumxm;
		          advection_a *= small_c;
		          advection_a *= (-wt_func);
		         }
		      else
		         {
		          advection_a *= - coeff_rho * wt_func;
		         }

		      advection_b = 0.;
		      if(taylor_galerkin[w])
			{
			  for ( p=0; p<VIM; p++)
			    {
			      advection_b += s_terms.taylor_flux[w][p];
			    }
			  advection_b *= - coeff_rho * s_terms.taylor_flux_wt[i]*dt/2.;
			}
		      advection = advection_a + advection_b;
		      advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)] * det_J * h3 * wt 
		       	           * mp->AdvectiveScaling[w];
		    }
		  
		  /*
		   * the diffusion term contains all the fluxes that are in
		   * a divergence in the mass conservation equation and are 
		   * integrated by parts in finite elements
		   *   Units of s_terms.diff_flux:
		   *      SPECIES_MASS_FRACTION:    gm cm-2 sec-1
		   *      SPECIES_MOLE_FRACTION:    mol cm-2 sec-1
		   */
		  diffusion = 0.;
		  if ( pd->e[pg->imtrx][eqn] & T_DIFFUSION )
		    {
		      for ( p=0; p<VIM; p++)
			{
			  grad_phi_i[p] = bf[eqn]->grad_phi[i] [p];
			}
		      
		      for ( p=0; p<VIM; p++)
			{
			  diffusion += grad_phi_i[p] * s_terms.diff_flux[w][p];
			}

                      if ((mp->SpeciesSourceModel[w]  == ELECTRODE_KINETICS ||
                           mp->SpeciesSourceModel[w]  == ION_REACTIONS)  &&
                          cr->MassFluxModel == STEFAN_MAXWELL_VOLUME) /*  RSL 3/19/01  */
		         {
		          diffusion *= small_c;
		         }

		      diffusion *= h3 * det_J * wt;
		      diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
		    }
		  /*
		   * HKM -> Note the addition of a species molecular weight
		   *        term is currently done in the source term
		   *        calculation section of the code
		   *   Units ofs_terms.MassSource:
		   *      SPECIES_MASS_FRACTION:    gm cm-3 sec-1
		   *      SPECIES_MOLE_FRACTION:    mol cm-3 sec-1
		   */
		  source = 0.;
		  if ( pd->e[pg->imtrx][eqn] & T_SOURCE )
		    {
		      source = s_terms.MassSource[w];

                      if (mp->SpeciesSourceModel[w]  == ELECTRODE_KINETICS ||
                          mp->SpeciesSourceModel[w]  == ION_REACTIONS) /*  RSL 3/19/01  */
		         {
		          sumrm = 0.;
		          for (j=0; j<num_species; j++)
		             {
		              sumrm += s_terms.MassSource[j] * M[j];
		             }
		          source -= x[w]*sumrm/sumxm;
		         }

		      source *= wt_func * h3 * det_J * wt;
		      source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
		    }

                   /*
		   *  Sum up all of the individual contributions and store it
		   *  in the local element residual vector.
		   */
                  lec->R[MAX_PROB_VAR + w][ii] +=
                        Heaviside*(mass + advection) +  diffusion + source;
		  
		}   /* if active_dofs */
	      
	    } /* end of loop over equations */
        } /* end of assemble residuals */
  

      /*
       * Jacobian terms...
       */

      if ( af->Assemble_Jacobian )
        {
          /*
           *         START loop over the rows corresponding to difference
           *	 species conservation equations
           *             w = row Species ktype
           *             i = node (i.e., dof) where species conservation
           *                 equation is located.
           */
          /*
           *  Store the species eqn type (which is keyed to the Variable
           *  type in a temporary variable).
           */
          species_eqn_type = mp->Species_Var_Type;
      

	  /*
	   *  Calculate the coef_rho term based upon the value of
	   *  species_eqn_type.
	   *     This coefficient indicates whether the time and advection
	   *     term should be multipled by the density or maybe the
	   *     concentration. 
	   */
	  coeff_rho = 1.0;
	  coeff_rho_nonunity = FALSE;
	  if ((species_eqn_type == SPECIES_MASS_FRACTION) ||
	      (species_eqn_type == SPECIES_MOLE_FRACTION)) {
	    coeff_rho = rho;
	    coeff_rho_nonunity = TRUE;	    
	  }

          if (mp->SpeciesSourceModel[w]  == ELECTRODE_KINETICS ||
              mp->SpeciesSourceModel[w]  == ION_REACTIONS) /*  RSL 3/19/01  */
	     {
	      sumxm = 0.0;
	      for (j=0; j<num_species; j++)
	         {
	          sumxm += x[j] * M[j];
	         }
	      if (mp->PorosityModel == CONSTANT)
	         {
	          epsilon = mp->porosity;
	         }
	      else if (mp->PorosityModel == THERMAL_BATTERY)
	         {
	          epsilon = mp->u_porosity[0];
	         }
	      else
	         {  
	          EH(-1, "invalid porosity model");
	         }
       /*        rho = density(d_rho); /\*  RSL 6/22/02  *\/ */
              small_c = rho/sumxm;
              coeff_rho = small_c; /*  RSL 9/27/01  */
	     }
	  
	  for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
	    ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];
	    if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
	      /*
	       *  Here is where we figure out whether the row is to placed in
	       *  the normal spot (e.g., ii = i), or whether a boundary condition
	       *  require that the volumetric contribution be stuck in another
	       *  ldof pertaining to the same variable type.
	       */
	      ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];		  
	      phi_i = bf[eqn]->phi[i];
	
		  wt_func = bf[eqn]->phi[i];
		  /* add Petrov-Galerkin terms as necessary */
		  if(supg!=0.0)
		    {
		      for(p=0; p<dim; p++)
			{
                          wt_func += supg * supg_terms.supg_tau * fv->v[p] * bf[eqn]->grad_phi[i][p];
			}
		    }
		  /*
		   * Set up some preliminaries that are needed for the (a,i)
		   * equation for bunches of (b,j) column variables...
		   */
		  
		  for ( p=0; p<VIM; p++)
		    {
		      grad_phi_i[p] = bf[eqn]->grad_phi[i][p];
		    }
		  
		  /* 
		   * J_s_c derivative of residual pieces w.r.t. the species
		   *       unknowns
		   */
		  var = MASS_FRACTION;
		  if ( pd->e[pg->imtrx][eqn] && pd->v[pg->imtrx][var] )
		    {
		      for ( w1=0; w1<pd->Num_Species_Eqn; w1++)
			{
			  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
			    {
			      phi_j = bf[var]->phi[j];	      
			      
			      mass = 0.;
			      if ( pd->TimeIntegration != STEADY )
				{
				  if ( pd->e[pg->imtrx][eqn] & T_MASS )
				  {
				    /*
				     * HKM -
				     *  Added dependence of rho on the species unknowns here.
				     */
				    mass  = coeff_rho * s_terms.d_Y_dot_dc[w][w1][j];
				    if (coeff_rho_nonunity) {
				      mass += d_rho->C[w1][j] * s_terms.Y_dot[w];
				    }

                                    if (mp->SpeciesSourceModel[w]  == ELECTRODE_KINETICS ||
                                        mp->SpeciesSourceModel[w]  == ION_REACTIONS) /*  RSL 3/19/01  */
				       {
				        sumxdot = 0.;
				        sumxdotm = 0.;
				        sumdxdot = 0.;
				        sumdxdotm = 0.;
				        for (b=0; b<num_species-1; b++)
				           {
				            sumxdot  += s_terms.Y_dot[b];
				            sumxdotm += s_terms.Y_dot[b] * M[b];
				            sumdxdot  += s_terms.d_Y_dot_dc[b][w1][j];
				            sumdxdotm += s_terms.d_Y_dot_dc[b][w1][j] * M[b];
				           }
				        sumxdotm -= sumxdot * M[num_species-1];
				        sumdxdotm -= sumdxdot * M[num_species-1];
				        sumxms = sumxm*sumxm;
				        mass = x[w]*d_rho->C[w1][j]/sumxms;
				        mass -= 2.*(M[w1] - M[num_species-1])*phi_j*
				                rho*x[w]/(sumxms*sumxm);
				        if (w == w1) mass += rho*phi_j/sumxms;
				        mass *= (-sumxdotm);
				        mass += s_terms.Y_dot[w]*d_rho->C[w1][j]/sumxm;
				        mass += rho*s_terms.d_Y_dot_dc[w][w1][j]/sumxm;
				        mass -= rho*s_terms.Y_dot[w]*
				                (M[w1] - M[num_species-1])*phi_j/sumxms;
				        mass -= rho*x[w]*sumdxdotm/sumxms;
				        mass *= epsilon;
				       }
				    mass *= - wt_func * h3 * det_J * wt;
				    mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
				  }
				}
			      
			      advection = 0.;
			      if ( pd->e[pg->imtrx][eqn] & T_ADVECTION )
				{
				  advection_a = 0.;
				  for ( p=0; p<VIM; p++)
				    {
				      advection_a += coeff_rho * s_terms.d_conv_flux_dc[w][p][w1][j];
				      if (coeff_rho_nonunity) {
				        advection_a += d_rho->C[w1][j] * s_terms.conv_flux[w][p];
				      }
				    }
				  
                                  if (mp->SpeciesSourceModel[w]  == ELECTRODE_KINETICS ||
                                      mp->SpeciesSourceModel[w]  == ION_REACTIONS) /*  RSL 3/19/01  */
				     {
				      psum1 = 0.;
				      psum2 = 0.;
				      psum3 = 0.;
				      psum4 = 0.;
				      for ( p=0; p<VIM; p++)
				         {
				          psum1 += s_terms.conv_flux[w][p];
				          psum2 += s_terms.d_conv_flux_dc[w][p][w1][j];
				          sumflux = 0.;
				          sumfluxm = 0.;
				          sumdflux = 0.;
				          sumdfluxm = 0.;
				          for (b=0; b<num_species-1; b++)
				             {
				              sumflux += s_terms.conv_flux[b][p];
				              sumfluxm += s_terms.conv_flux[b][p] * M[b];
				              sumdflux += s_terms.d_conv_flux_dc[b][p][w1][j];
				              sumdfluxm += s_terms.d_conv_flux_dc[b][p][w1][j] * M[b];
				             }
				          sumfluxm -= sumflux * M[num_species-1];
				          sumdfluxm -= sumdflux * M[num_species-1];
				          psum3 += sumdfluxm;
				          psum4 += sumfluxm;
				         }
				      sumxms = sumxm*sumxm;
				      advection_a = x[w]*d_rho->C[w1][j]/sumxms;
				      dmdx = (M[w1] - M[num_species-1])*phi_j;
				      advection_a -= 2.*dmdx*rho*x[w]/(sumxms*sumxm);
				      if (w == w1) advection_a += rho*phi_j/sumxms;
				      advection_a *= (-psum4);
				      advection_a += psum1*d_rho->C[w1][j]/sumxm;
				      advection_a += rho*psum2/sumxm;
				      advection_a -= rho*psum1*dmdx/sumxms;
				      advection_a -= rho*x[w]*psum3/sumxms;
				     }
				  
				  advection_a *= - wt_func;

				  advection_b = 0.;
				  if (taylor_galerkin[w])
				    {
				      for ( p=0; p<VIM; p++) {
					advection_b += coeff_rho * s_terms.d_taylor_flux_dc[w][p][w1][j];
					if (coeff_rho_nonunity) {
					  advection_b += d_rho->C[w1][j] * s_terms.conv_flux[w][p];
					}
				      }
				      advection_b *= -  s_terms.taylor_flux_wt[i]*dt/2.;
				    }
				  
				  advection = advection_a + advection_b;  
				  advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)]* h3 * det_J * wt
				    * mp->AdvectiveScaling[w];
				}
			      
			      diffusion = 0.;
			      if ( pd->e[pg->imtrx][eqn] & T_DIFFUSION )
				{
				  for ( p=0; p<VIM; p++)
				    {
				      diffusion += grad_phi_i[p]
					* s_terms.d_diff_flux_dc[w][p] [w1][j];
				    }

                                  if ((mp->SpeciesSourceModel[w]  == ELECTRODE_KINETICS ||
                                       mp->SpeciesSourceModel[w]  == ION_REACTIONS)  &&
                                      cr->MassFluxModel == STEFAN_MAXWELL_VOLUME) /*  RSL 3/19/01  */
				     {
				      diffusion *= small_c;
				     }
				  
				  diffusion *= h3 * det_J * wt;
				  diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
				}	  
			      
			      source = 0.;
			      if ( pd->e[pg->imtrx][eqn] & T_SOURCE )
				{
				  source += s_terms.d_MassSource_dc[w][w1][j];

                                  if (mp->SpeciesSourceModel[w] == ELECTRODE_KINETICS ||
                                      mp->SpeciesSourceModel[w]  == ION_REACTIONS) /*  RSL 3/19/01  */
				     {
				      sumrm = 0.;
				      sumdrm = 0.;
				      for (b=0; b<num_species; b++)
				         {
				          sumrm += s_terms.MassSource[b] * M[b];
				          sumdrm += s_terms.d_MassSource_dc[b][w1][j] * M[b];
				         }
				      source -= x[w]*sumdrm/sumxm;
				      factor = -x[w]*(M[w1] - M[num_species-1])/sumxm;
				      if (w == w1) factor += 1.;
				      source -= sumrm*phi_j*factor/sumxm;
				     }

				  source *= wt_func;
				  source *= h3 * det_J * wt;
				  source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
				}
			      
                              lec->J[MAX_PROB_VAR + w][MAX_PROB_VAR + w1][ii][j] +=
                                    Heaviside*(mass + advection) + diffusion + source;
			    }
			}
		    }
		  
		  /* MMH
		   * J_s_v  sensitivity of species equation w.r.t. velocity
		   *
		   * NOTE: If DensityModel == SUSPENSION_PM && w == species 
		   * then this is really J_s_pv.  
		   */
		  for ( b=0; b<dim; b++)
		    {
		      /* MMH */
		      if(mp->DensityModel == SUSPENSION_PM &&
			 w==species)
			var = PVELOCITY1+b;
		      else
			var = VELOCITY1+b;
		      if ( pd->v[pg->imtrx][var] )
			{
			  pvar = upd->vp[pg->imtrx][var];
			  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
			    {
			      phi_j = bf[var]->phi[j];	

			      d_wt_func = 0.;

                              if(supg != 0.)
                                {
                                  d_wt_func = supg * supg_terms.supg_tau * phi_j * bf[eqn]->grad_phi[i][b];
                                  for( p=0; p<dim; p++ )
                                    {
                                      d_wt_func += supg * supg_terms.d_supg_tau_dv[b][j] * fv->v[p]* bf[eqn]->grad_phi[i][p];
                                    }
                                }
			      
			      mass = 0.;

			      if ( (supg != 0.0) && (pd->e[pg->imtrx][eqn] & T_MASS) ) {
				  mass = coeff_rho * s_terms.Y_dot[w];

				  if (mp->SpeciesSourceModel[w]  == ELECTRODE_KINETICS ||
				      mp->SpeciesSourceModel[w]  == ION_REACTIONS) /*  RSL 3/19/01  */
				  {
				    mass = s_terms.Y_dot[w];
				    sumxdot = 0.0;
				    sumxdotm = 0.0;
				    for (j = 0; j < num_species - 1; j++) {
				      sumxdotm += s_terms.Y_dot[j] * M[j];
				      sumxdot  += s_terms.Y_dot[j];
				    }
				    sumxdotm -= sumxdot * M[num_species-1];
				    mass -= x[w]*sumxdotm/sumxm;
				    mass *= (epsilon*small_c);
				  }

				  mass *= - d_wt_func * h3 * det_J * wt;
				  mass *= pd->etm[pg->imtrx][eqn][LOG2_MASS];
				}
			      
			      advection = 0.;
			      if ( pd->e[pg->imtrx][eqn] & T_ADVECTION )
				{
				  advection_a =  -wt_func * coeff_rho * s_terms.d_conv_flux_dv[w][b] [b][j];

                                  if (mp->SpeciesSourceModel[w]  == ELECTRODE_KINETICS ||
                                      mp->SpeciesSourceModel[w]  == ION_REACTIONS) /*  RSL 3/19/01  */
				     {
				      sumdflux = 0.;
				      sumdfluxm = 0.;
				      for (q=0; q<num_species-1; q++)
				         {
				          sumdflux += s_terms.d_conv_flux_dv[q][b][b][j];
				          sumdfluxm += s_terms.d_conv_flux_dv[q][b][b][j] * M[q];
				         }
				      sumdfluxm -= sumdflux * M[num_species-1];
				      advection_a = -x[w]*sumdfluxm/sumxm;
				      advection_a += s_terms.d_conv_flux_dv[w][b][b][j];
				      advection_a *= small_c;
				      advection_a *= (-wt_func);
				     }

				  advection_b = 0.;
				  if(supg!=0.)
				    {				      
				      for(p=0;p<dim;p++)
					{
					  advection_b +=  coeff_rho * s_terms.conv_flux[w][p];
					}
				      
				      advection_b *=  -d_wt_func;
				    }
				  
				  advection_c = 0.;
				  if(taylor_galerkin[w])
				    {
				      for ( p=0; p<VIM; p++)
					{
					  advection_c += s_terms.taylor_flux[w][p];
					}
				      advection_c *= - coeff_rho * s_terms.d_taylor_flux_wt_dv[i] [b][j]*dt/2.
					-s_terms.taylor_flux_wt[i]*dt/2. 
					* s_terms.d_taylor_flux_dv[w][b] [b][j];
				    }
				  
				  advection = advection_a + advection_b + advection_c;
				  advection *= det_J * h3 * wt;
				  advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)]
				    * mp->AdvectiveScaling[w];
				}
			      
			      diffusion = 0.;
			      
			      if ( pd->e[pg->imtrx][eqn] & T_DIFFUSION )
				{
				  for ( p=0; p<VIM; p++)
				    {
				      diffusion += grad_phi_i[p] * s_terms.d_diff_flux_dv[w][p][b][j];
				    }
				  diffusion *= h3 * det_J * wt * pd->etm[pg->imtrx][eqn][LOG2_DIFFUSION]; 
				}
			      
			      source = 0.;
			      if ( pd->e[pg->imtrx][eqn] & T_SOURCE )
				{
				  source += wt_func * s_terms.d_MassSource_dv[w][b][j];
				  source *= h3 * det_J * wt;
				  source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

				  double source_b = 0;
				  if (supg != 0)
				    {
				      source_b = s_terms.MassSource[w];

				      if (mp->SpeciesSourceModel[w]  == ELECTRODE_KINETICS ||
					  mp->SpeciesSourceModel[w]  == ION_REACTIONS) /*  RSL 3/19/01  */
					 {
					  sumrm = 0.;
					  for (j=0; j<num_species; j++)
					     {
					      sumrm += s_terms.MassSource[j] * M[j];
					     }
					  source_b -= x[w]*sumrm/sumxm;
					 }

				      source_b *= d_wt_func * h3 * det_J * wt;
				      source_b *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
				    }
				}
			      
			      lec->J[MAX_PROB_VAR + w][pvar][ii][j] += 
                                Heaviside*advection+ diffusion + source;
			      
			    }
			}
		    }
		  
		  /* 
		   * J_s_vd  sensitivity of species equation w.r.t. vorticity direction
		   *
		   */
		  for ( b=0; b<dim; b++)
		    {
		      var = VORT_DIR1+b;
		      if ( pd->v[pg->imtrx][var] )
			{
			  pvar = upd->vp[pg->imtrx][var];
			  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
			    {
			      phi_j = bf[var]->phi[j];	
			      
			      diffusion = 0.;
			      
			      if ( pd->e[pg->imtrx][eqn] & T_DIFFUSION )
				{
				  for ( p=0; p<VIM; p++)
				    {
				      diffusion += grad_phi_i[p] * s_terms.d_diff_flux_dvd[w][p][b][j];
				    }
				  diffusion *= h3 * det_J * wt * pd->etm[pg->imtrx][eqn][LOG2_DIFFUSION]; 
				}
			      
			      lec->J[MAX_PROB_VAR + w][pvar][ii][j] += 
				diffusion;
			      
			    }
			}
		    }

		  /*
		   * J_s_d  sensitivity of species equation w.r.t. mesh displacement
		   */
		  for ( b=0; b<dim; b++)
		    {
		      var = MESH_DISPLACEMENT1+b;
		      if ( pd->v[pg->imtrx][var] )
			{
			  pvar = upd->vp[pg->imtrx][var];
			  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
			    {
			      phi_j = bf[var]->phi[j];	      
			      
			      d_det_J_dmeshbj = bf[var]->d_det_J_dm[b][j];
			      
			      dh3dmesh_bj = fv->dh3dq[b] * bf[var]->phi[j];
			      

                              d_wt_func = 0.;

                              if(supg != 0.)
                                {
                                  for( p=0; p<dim; p++ )
                                    {
                                      d_wt_func += supg *
                                        (supg_terms.d_supg_tau_dX[b][j] * fv->v[p] * bf[eqn]->grad_phi[i][p]
                                         + supg_terms.supg_tau * fv->v[p] * bf[eqn]->d_grad_phi_dmesh[i][p][b][j]);
                                    }
                                }
			      
			      mass = 0.;
			      if ( pd->TimeIntegration != STEADY )
				{
				  if ( pd->e[pg->imtrx][eqn] & T_MASS )
				    {
				      mass  = coeff_rho * s_terms.Y_dot[w];

                                      if (mp->SpeciesSourceModel[w]  == ELECTRODE_KINETICS
 ||
                                          mp->SpeciesSourceModel[w]  == ION_REACTIONS) /*  RSL 9/27/01  */
                                         {
                                          mass = s_terms.Y_dot[w];
                                          sumxdot = 0.;
                                          sumxdotm = 0.;
                                          for (q=0; q<num_species-1; q++)
                                             {
                                              sumxdotm += s_terms.Y_dot[q] * M[q];
                                              sumxdot  += s_terms.Y_dot[q];
                                             }
                                          sumxdotm -= sumxdot * M[num_species-1];
                                          mass -= x[w]*sumxdotm/sumxm;
                                          mass *= (epsilon*small_c);
                                         }

                                      mass *= - (d_wt_func * h3 * det_J + wt_func *
					( h3 * d_det_J_dmeshbj 
                                          + dh3dmesh_bj * det_J )) * wt;
				      mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
				    }
				}
			      
			      advection = 0.;
			      
			      if ( pd->e[pg->imtrx][eqn] & T_ADVECTION )
				{
				  /*
				   * Two parts:
				   *
				   *	Int( (v- xdot).d(Vc)/dmesh h3 |Jv| )
				   *	Int( v.Vc d(h3|Jv|)/dmesh )
				   */
				  advection_a = 0.;
				  
				  for ( p=0; p<VIM; p++)
				    {
				      advection_a += s_terms.d_conv_flux_dmesh[w][p] [b][j];
				    }

                                  if (mp->SpeciesSourceModel[w]  == ELECTRODE_KINETICS ||
                                      mp->SpeciesSourceModel[w]  == ION_REACTIONS) /*  RSL 9/27/01  */
                                     {
                                      psum = 0.;
                                      for ( p=0; p<VIM; p++)
                                         {
                                          sumdflux = 0.;
                                          sumdfluxm = 0.;
                                          for (q=0; q<num_species-1; q++)
                                             {
                                              sumdfluxm += s_terms.d_conv_flux_dmesh[q][p]
 [b][j] * M[q];
                                              sumdflux  += s_terms.d_conv_flux_dmesh[q][p]
 [b][j];
                                             }
                                          sumdfluxm -= sumdflux * M[num_species-1];
                                          psum += sumdfluxm;
                                         }
                                      advection_a -= x[w]*psum/sumxm;
                                      advection_a *= - small_c * wt_func * h3 * det_J * wt
;
                                     }
                                  else
                                     {
                                      advection_a *= - coeff_rho * wt_func * h3 * det_J *wt;
                                     }
				  
				  advection_b = 0.;
				  
				  for ( p=0; p<VIM; p++)
				    {
				      advection_b += s_terms.conv_flux[w][p];
				    }

                                  if (mp->SpeciesSourceModel[w]  == ELECTRODE_KINETICS ||
                                      mp->SpeciesSourceModel[w]  == ION_REACTIONS) /*  RSL 9/27/01  */
                                     {
                                      psum = 0.;
                                      for ( p=0; p<VIM; p++)
                                         {
                                          sumflux = 0.;
                                          sumfluxm = 0.;
                                          for (q=0; q<num_species-1; q++)
                                             {
                                              sumfluxm += s_terms.conv_flux[q][p] * M[q];
                                              sumflux  += s_terms.conv_flux[q][p];
                                             }
                                          sumfluxm -= sumflux * M[num_species-1];
                                          psum += sumfluxm;
                                         }
                                      advection_b -= x[w]*psum/sumxm;
                                      advection_b *= - small_c * wt_func;
                                     }
                                  else
                                     {
                                      advection_b *= - coeff_rho * wt_func;
                                     }

				  if(taylor_galerkin[w])
				    {
				      advection_c = 0.;
				      for ( p=0; p<VIM; p++)
					{
                                          advection_c += s_terms.taylor_flux[w][p];
					}
				      advection_c *= - coeff_rho * s_terms.taylor_flux_wt[i]*dt/2.;
				      advection_b += advection_c;
				    }
				  
				  advection_b *= ( h3 * d_det_J_dmeshbj 
						   +dh3dmesh_bj * det_J) * wt;
				  
				  advection_f = 0.;
				  if(supg != 0.)
				    {
				      for( p=0; p<dim; p++ )
					{ 
					  advection_f += s_terms.conv_flux[w][p];
					}
				      advection_f *= -d_wt_func * coeff_rho * h3 * det_J * wt;
				    }
				  
				  advection_d = 0.;
				  if(taylor_galerkin[w])
				    {
				      advection_e = 0.;
				      for ( p=0; p<VIM; p++)
					{
					  advection_d += s_terms.d_taylor_flux_dmesh[w][p] [b][j];
					  advection_e += s_terms.taylor_flux[w][p];
					}
				      
				      advection_d *= - coeff_rho * s_terms.taylor_flux_wt[i]*dt/2.;
				      advection_d += - coeff_rho * advection_e *
					               s_terms.d_taylor_flux_wt_dmesh[i] [b][j]*dt/2.;
				      
				      advection_d *= h3 * det_J * wt;
				    }
				  
				  advection = advection_a + advection_f
				            + advection_b + advection_d;
				  
				  advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)]
				    * mp->AdvectiveScaling[w];
				}
			      
			      
			      diffusion = 0.;
			      
			      if ( pd->e[pg->imtrx][eqn] & T_DIFFUSION )
				{
				  /*
				   * Three parts:
				   *  diff_a = Int(d(grad_phi_i)/dmesh.q h3 |Jv|)
				   *  diff_b = Int(grad_phi_i.d(q)/dmesh h3 |Jv|)
				   *  diff_c = Int(grad_phi_i.q d(h3 |Jv|)/dmesh)
				   */
				  
				  diff_a = 0.;
				  
				  for ( p=0; p<VIM; p++)
				    {
				      dgrad_phi_i_dmesh[p]
					= bf[eqn]->d_grad_phi_dmesh[i][p] [b][j];
				      
				      diff_a += dgrad_phi_i_dmesh[p] * s_terms.diff_flux[w][p];
				    }
				  
				  diff_a *= h3 * det_J * wt;
				  
				  diff_b = 0.;
				  for ( p=0; p<VIM; p++)
				    {
				      diff_b += bf[eqn]->grad_phi[i] [p]* 
					s_terms.d_diff_flux_dmesh[w][p] [b][j];
				    }
				  diff_b *= h3 * det_J * wt;
				  
				  diff_c = 0.;
				  for ( p=0; p<VIM; p++)
				    {
				      diff_c +=
					bf[eqn]->grad_phi[i] [p] 
					* s_terms.diff_flux[w][p];
				    }
				  diff_c *= ( h3 * d_det_J_dmeshbj 
					      + dh3dmesh_bj * det_J ) * wt;
				  
				  
				  diffusion = diff_a + diff_b + diff_c;

                                  if ((mp->SpeciesSourceModel[w]  == ELECTRODE_KINETICS ||
                                       mp->SpeciesSourceModel[w]  == ION_REACTIONS)  &&
                                      cr->MassFluxModel == STEFAN_MAXWELL_VOLUME) /*  RSL 9/27/01  */
                                     {
                                      diffusion *= small_c;
                                     }

				  diffusion *= pd->etm[pg->imtrx][eqn][LOG2_DIFFUSION];
				}
			      
			      source = 0.;
			      if ( pd->e[pg->imtrx][eqn] & T_SOURCE )
				{
				  source += s_terms.MassSource[w];

                                  if (mp->SpeciesSourceModel[w]  == ELECTRODE_KINETICS ||
                                      mp->SpeciesSourceModel[w]  == ION_REACTIONS) /*  RSL 9/27/01  */
                                     {
                                      sumrm = 0.;
                                      for (q=0; q<num_species; q++)
                                         {
                                          sumrm += s_terms.MassSource[q] * M[q];
                                         }
                                      source -= x[w]*sumrm/sumxm;
                                     }

                                  source *= (d_wt_func * h3 * det_J + wt_func * ( h3 * d_det_J_dmeshbj + dh3dmesh_bj * det_J ))
                                    * wt;
				  source += s_terms.d_MassSource_dmesh[w][b][j]*det_J*h3*wt*wt_func;

                                  if (mp->SpeciesSourceModel[w]  == ELECTRODE_KINETICS ||
                                      mp->SpeciesSourceModel[w]  == ION_REACTIONS) /*  RSL 9/27/01  */
                                     {
                                      sumrm = 0.;
                                      for (q=0; q<num_species; q++)
                                         {
                                          sumrm += s_terms.d_MassSource_dmesh[q][b][j] * M[q];
                                         }
                                      sumrm *= (-x[w]/sumxm);
                                      source += sumrm*det_J*h3*wt*wt_func;
                                     }

				  source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];
				}
			      
                              lec->J[MAX_PROB_VAR + w][pvar][ii][j] +=
                                    Heaviside*(mass + advection) + diffusion + source;
			    }
			}
		    }
		  
		  
		  /*
		   * J_s_T  sensitivity of species equation w.r.t. temperature
		   */
		  var = TEMPERATURE;
		  if ( pd->e[pg->imtrx][eqn] && pd->v[pg->imtrx][var] )
		    {
		      pvar = upd->vp[pg->imtrx][var];
		      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
			{
			  phi_j = bf[var]->phi[j];	      

			  /*
			   * HKM -
			   *    Added in the density dependence on temperature into
			   *    the time derivative term
			   */
			  mass = 0.0;
			  if (pd->TimeIntegration != STEADY) {
			    if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
			      if (coeff_rho_nonunity) {
				mass = d_rho->T[j] * s_terms.Y_dot[w];
				mass *= - wt_func * h3 * det_J * wt;
				mass *= pd->etm[pg->imtrx][eqn] [(LOG2_MASS)];
			      }
			    }
			  }
			  
			  advection = 0.;
			  if ( pd->e[pg->imtrx][eqn] & T_ADVECTION )
			    {
			      for ( p=0; p<VIM; p++)
				{
				  advection += coeff_rho * s_terms.d_conv_flux_dT[w][p] [j];
				  if (coeff_rho_nonunity) {
				    advection += d_rho->T[j] * s_terms.conv_flux[w][p];
				  }
				}
			      advection *= - wt_func;
			      
			      if(taylor_galerkin[w])
				{
				  advection_a = 0.;
				  for ( p=0; p<VIM; p++) {
				    advection_a += coeff_rho * s_terms.taylor_flux[w][p];
				    if (coeff_rho_nonunity) {
				      advection += d_rho->T[j] * s_terms.conv_flux[w][p];
				    }
				  }
				  advection *= -advection_a *s_terms.d_taylor_flux_wt_dT[i] [j]*dt/2.
				    -s_terms.taylor_flux_wt[i]*dt/2. * s_terms.d_taylor_flux_dT[w][b] [j];
				}
			      
			      advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)] * h3 * det_J * wt
				* mp->AdvectiveScaling[w];
			    }
			  
			  diffusion = 0.;
			  if ( pd->e[pg->imtrx][eqn] & T_DIFFUSION )
			    {
			      for ( p=0; p<VIM; p++)
				{
				  diffusion += grad_phi_i[p]
				    * s_terms.d_diff_flux_dT[w][p] [j];
				}
			      diffusion *= h3 * det_J * wt;
			      diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
			    }
			  
			  source = 0.;
			  if ( pd->e[pg->imtrx][eqn] & T_SOURCE )
			    {
			      source += s_terms.d_MassSource_dT[w][j];

                              if (mp->SpeciesSourceModel[w] == ELECTRODE_KINETICS ||
                                  mp->SpeciesSourceModel[w]  == ION_REACTIONS) /*  RSL 3/19/01  */
			         {
			          sumdrm = 0.;
			          for (w1=0; w1<num_species; w1++)
			             {
			              sumdrm += s_terms.d_MassSource_dT[w1][j] * M[w1];
			             }
			          source -= x[w]*sumdrm/sumxm;
			         }

			      source *= det_J*h3*wt*wt_func;
			      source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];
                            }

                            lec->J[MAX_PROB_VAR + w][pvar][ii][j] += Heaviside*(mass + advection) + diffusion + source;
			}		/* for(j) .... */
		    }			/* if ( e[eqn], v[var]) .... */	      
		  var = LIGHT_INTP;
		  if ( pd->e[pg->imtrx][eqn] && pd->v[pg->imtrx][var] )
		    {
		      pvar = upd->vp[pg->imtrx][var];
		      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
			{
			  phi_j = bf[var]->phi[j];	      
			  source = 0.;
			  if ( pd->e[pg->imtrx][eqn] & T_SOURCE )
			    {
			      source += s_terms.d_MassSource_dI[w][j];
			      source *= det_J*h3*wt*wt_func;
			      source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];
			    }
			  
			  lec->J[MAX_PROB_VAR + w][pvar][ii][j] += source;
			}		/* for(j) .... */
		    }			/* if ( e[eqn], v[var]) .... */	      
		  var = LIGHT_INTM;
		  if ( pd->e[pg->imtrx][eqn] && pd->v[pg->imtrx][var] )
		    {
		      pvar = upd->vp[pg->imtrx][var];
		      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
			{
			  phi_j = bf[var]->phi[j];	      
			  source = 0.;
			  if ( pd->e[pg->imtrx][eqn] & T_SOURCE )
			    {
			      source += s_terms.d_MassSource_dI[w][j];
			      source *= det_J*h3*wt*wt_func;
			      source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];
			    }
			  
			  lec->J[MAX_PROB_VAR + w][pvar][ii][j] += source;
			}		/* for(j) .... */
		    }			/* if ( e[eqn], v[var]) .... */	      
		  var = LIGHT_INTD;
		  if ( pd->e[pg->imtrx][eqn] && pd->v[pg->imtrx][var] )
		    {
		      pvar = upd->vp[pg->imtrx][var];
		      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
			{
			  phi_j = bf[var]->phi[j];	      
			  source = 0.;
			  if ( pd->e[pg->imtrx][eqn] & T_SOURCE )
			    {
			      source += s_terms.d_MassSource_dI[w][j];
			      source *= det_J*h3*wt*wt_func;
			      source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];
			    }
			  
			  lec->J[MAX_PROB_VAR + w][pvar][ii][j] += source;
			}		/* for(j) .... */
		    }			/* if ( e[eqn], v[var]) .... */	      
		  
		  
		  /*
		   * J_s_V  sensitivity of species equation w.r.t. voltage -- RSL 4/4/00
		   */
               if ( cr->MassFluxModel != FICKIAN_CHARGED) 
                {
		  var = VOLTAGE;
		  if ( pd->e[pg->imtrx][eqn] && pd->v[pg->imtrx][var] )
		    {
		      pvar = upd->vp[pg->imtrx][var];
		      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
			{
			  phi_j = bf[var]->phi[j];	      

			  mass = 0.0;
			  
			  advection = 0.;
			  
			  diffusion = 0.;
			  if ( pd->e[pg->imtrx][eqn] & T_DIFFUSION )
			    {
			      for ( p=0; p<VIM; p++)
				{
				  diffusion += grad_phi_i[p]
				    * s_terms.d_diff_flux_dV[w][p] [j];
				}
		  				  
                                  if ((mp->SpeciesSourceModel[w]  == ELECTRODE_KINETICS ||
                                       mp->SpeciesSourceModel[w]  == ION_REACTIONS)  &&
                                      cr->MassFluxModel == STEFAN_MAXWELL_VOLUME) /*  RSL 3/19/01  */
				     {
				      diffusion *= small_c;
				     }
				  
			      diffusion *= h3 * det_J * wt;
			      diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
			    }
			  
			  source = 0.;
			  if ( pd->e[pg->imtrx][eqn] & T_SOURCE )
			    {
			      source += s_terms.d_MassSource_dV[w][j];

                              if (mp->SpeciesSourceModel[w] == ELECTRODE_KINETICS ||
                                  mp->SpeciesSourceModel[w] == ION_REACTIONS) /*  RSL 3/19/01  */
			         {
			          sumdrm = 0.;
			          for (w1=0; w1<num_species; w1++)
			             {
			              sumdrm += s_terms.d_MassSource_dV[w1][j] * M[w1];
			             }
			          source -= x[w]*sumdrm/sumxm;
			         }

			      source *= det_J*h3*wt*wt_func;
			      source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];
                            }

                            lec->J[MAX_PROB_VAR + w][pvar][ii][j] += Heaviside * (mass + advection) + diffusion + source;
			}  /* end of loop over j */
		    }  /* end of var = VOLTAGE */
                }

	       /*  KSC: 9/9/00
	        * J_s_V  sensitivity of species equation w.r.t. electrical potential for 
                * the case of Fickian diffusion of charged species  
	        */
               else if ( cr->MassFluxModel == FICKIAN_CHARGED) 
                {
		  var = VOLTAGE;
		  if ( pd->e[pg->imtrx][eqn] && pd->v[pg->imtrx][var] )
		    {
		      pvar = upd->vp[pg->imtrx][var];
		      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
			{
			  diffusion = 0.;
			  if ( pd->e[pg->imtrx][eqn] & T_DIFFUSION )
			    {
			      for ( p=0; p<VIM; p++)
				{
				  diffusion += grad_phi_i[p]
				    * s_terms.d_diff_flux_dV[w][p] [j];
				}
			      diffusion *= h3 * det_J * wt;
			      diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
			    }
			  
			  lec->J[MAX_PROB_VAR + w][pvar][ii][j] += diffusion;
			}		/* for(j) .... */
		    }			/* if ( e[eqn], v[var]) .... */	      
                }                       /* if cr->MassFluxModel == FICKIAN_CHARGED ...  */   

		  /*
		   * Jacobian with respect to pressure
		   */
		  var = PRESSURE; 
		  if ( pd->e[pg->imtrx][eqn] && pd->v[pg->imtrx][var] )
		    {
		      pvar = upd->vp[pg->imtrx][var];
		      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
			{
			  phi_j = bf[var]->phi[j];	      
			  /*
			   *  HKM - 
			   *   Right now, the density does not depend on pressure. We could
			   *   easily add it in here
			   */
			  mass = 0.;
			  if ( pd->e[pg->imtrx][eqn] & T_MASS )
			    {
			      mass += coeff_rho * s_terms.d_Y_dot_dP[w] [j];
			      mass *= - wt_func * det_J * h3 * wt;
			      mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
			    }
			  
			  diffusion = 0.;
			  if ( pd->e[pg->imtrx][eqn] & T_DIFFUSION )
			    {
			      for ( p=0; p<VIM; p++)
				{
				  diffusion += grad_phi_i[p]
				    * s_terms.d_diff_flux_dP[w][p] [j];
				}
			      diffusion *= h3 * det_J * wt;
			      diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
			    }
			  
			  /* Advection is velocity times gradient of 
			     mass (or volume) fraction */
			  advection = 0.;
			  if ( pd->e[pg->imtrx][eqn] & T_ADVECTION )
			    {
			      for ( p=0; p<VIM; p++)
				{
				  advection += rho * s_terms.d_conv_flux_dP[w][p] [j];
				}
			      advection *= - wt_func ;
			      
			      advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)]* det_J * h3 * wt
				* mp->AdvectiveScaling[w];
			    }
			  
                          lec->J[MAX_PROB_VAR + w][pvar][ii][j] += Heaviside * advection + diffusion + mass;
			}		/* for(j) .... */
		    }			/* if ( e[eqn], v[var]) .... */
		  
		  /*
		   * J_s_SH:
		   *      Dependence of the species continuity equation on the Shear Rate
		   *      inpendent unknown.
		   */
		  
		  if ( cr->MassFluxModel == HYDRODYNAMIC || cr->MassFluxModel == DM_SUSPENSION_BALANCE ) /* These terms only appear for this model */
		    {
		      var = SHEAR_RATE;
		      if ( pd->v[pg->imtrx][var])
			{
			  pvar = upd->vp[pg->imtrx][var];
			  for( j=0;  j<ei[pg->imtrx]->dof[var] ; j++)
			    {
			      diffusion = 0.0;
			      if ( pd->e[pg->imtrx][eqn] & T_DIFFUSION )
				for ( p=0; p<VIM; p++)
				  {
				    diffusion += grad_phi_i[p]
				      * s_terms.d_diff_flux_dSH[w][p][j];
				  }
			      diffusion *= h3 * det_J * wt;
			      diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
			      
			      lec->J[MAX_PROB_VAR + w][pvar][ii][j] += diffusion;
			      
			    } /* for (j) .. J_s_SH */
			}/* if (pd) */
		    } /* if( cr) */ 

		  /*
		   * J_s_G:
		   *      Dependence of the species continuity equation on the velocity gradient
		   *      tensor unknowns.
		   */
		   /* These terms only appear for this model */
		  if ( cr->MassFluxModel == DM_SUSPENSION_BALANCE )
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

		      for ( b=0; b<VIM; b++)
			{
			  for ( c=0; c<VIM; c++)
			    {
			      var = v_g[b][c];		      
			      if ( pd->v[pg->imtrx][var])
				{
				  pvar = upd->vp[pg->imtrx][var];
				  for( j=0;  j<ei[pg->imtrx]->dof[var] ; j++)
				    {
				      diffusion = 0.0;
				      if ( pd->e[pg->imtrx][eqn] & T_DIFFUSION )
					for ( p=0; p<VIM; p++)
					  {
					    diffusion += grad_phi_i[p]
					      * s_terms.d_diff_flux_dG[w][p][b][c][j];
					  }
				      diffusion *= h3 * det_J * wt;
				      diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
				      
				      lec->J[MAX_PROB_VAR + w][pvar][ii][j] += diffusion;
				    }
				}
			    }
			}
		    } /*if ( cr->MassFluxModel == DM_SUSPENSION_BALANCE ) */
			      

		}     /* if active_dofs */
            }				/* for (i) .... */
        }					/* if ( assemble Jacobian ) */
    }				/* for (w) ... */
  
  return(status);
  
} /* end of assemble_mass_transport */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int 
assemble_mass_transport_path_dependence
                       (double time, /* present time valuel; KSC             */
			double tt, /* parameter to vary time integration from 
				    * explicit (tt = 1) to implicit (tt = 0) */
			double dt, /* current time step size */
			const dbl h[DIM], /* element sizes, not scale factors.     */
			const dbl hh[DIM][DIM],
			const dbl dh_dxnode[DIM][MDE],
			const dbl vcent[DIM],	/* average element velocity, which is
					 * the centroid velocity for Q2 and 
					 * the average of the vertices for Q1.
					 * From routine "element_velocity."  */
			const dbl dvc_dnode[DIM][MDE])
{
  int var, ii,  pvar, ledof;
  const int eqn = R_MASS;
  const int dim = pd->Num_Dim;
  int p, w,  i, j, status;
  /*
   *    species_eqn_type:
   *        This is a temp variable that gets set at the top
   *        from Species_Var_Type.
   *        Specifies what identity of the independent variable
   *        The form of the equation and the units for the equation
   *        are dependent on the this determination.
   */
  int species_eqn_type;
  /*
   *    coeff_rho:
   *        If coeff_rho_nonunity is TRUE, it means that the time
   *        derivative and advection term have variable density or
   *        concentration terms in them that are not part of the
   *        independent variable itself.
   *
   *              Species_Var_type   Time_Deriv_Term  coeff_rho_nonunity
   *             ---------------------------------------------------------
   *             SPECIES_MASS_FRACTION     rho * d(Y_k)/dt       TRUE
   *             SPECIES_MOLE_FRACTION     conc * d(X_k)/dt      TRUE
   *             SPECIES_CONCENTRATION     d(conc_k)/dt         FALSE
   *                DEFAULT                d(????)/dt           FALSE
   *           ---------------------------------------------------------
   *
   *        coeff_rho is equal to one for coeff_rho_nonunity FALSE
   *        cases. It is non-unity and variable (with units) for
   *        for coeff_rho_nonunity TRUE cases.
   */
  dbl coeff_rho;
  dbl rho;				/* Density. */
  dbl sumxm = 0.0, sumrm, epsilon = 0.0, small_c = 0.0, sumxdot, sumxdotm; 
  dbl psum, sumflux, sumfluxm;  
  dbl sumx;
  int num_species;  
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;  /* density dependence */
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;

  dbl M[MAX_CONC];            /* species molecular weight */ 
  dbl x[MAX_CONC];            /* mole fraction */

  struct Species_Conservation_Terms s_terms;

  dbl mass;		         	/* For terms and their derivatives */

  dbl advection;			/* For terms and their derivatives */
  dbl advection_a, advection_b;
  dbl diffusion;
  dbl source;

  /*
   * Galerkin weighting functions for i-th energy residuals 
   * and some of their derivatives...
   */

  dbl phi_i;
  dbl grad_phi_i[DIM];

  /*
   * Petrov-Galerkin weighting functions for i-th residuals 
   * and some of their derivatives...
   */

  dbl wt_func;

  /* SUPG variables */
  dbl h_elem=0;
  dbl supg=0;

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl h3;			/* Product of all 3 scale factors used for */
				/* volume integrals in this coordinate system*/
  dbl det_J;
  dbl wt;
  int sign = ls->Elem_Sign;
  int err;
  int taylor_galerkin[MAX_CONC];
  dbl residual;
  status = 0;


  /*
   * Unpack variables from structures for local convenience...
   */



  /*
   * Bail out fast if there's nothing to do...
   */

  if ( ! pd->e[pg->imtrx][eqn] )
    {
      return(status);
    }
    
  wt = fv->wt;				/* Gauss point weight. */
  h3 = fv->h3;			        /* Differential volume element. */
  
  if ( ls->on_sharp_surf ) /* sharp interface */
	{det_J = fv->sdet;}
  else
    	{det_J = bf[eqn]->detJ; }

  for ( w=0; w<pd->Num_Species_Eqn; w++)
    {
      taylor_galerkin[w] =mp->SpeciesTimeIntegration[w]; 
    }

  num_species = pd->Num_Species;  

  for (j=0; j<num_species; j++)  
     {
      M[j] = mp->molecular_weight[j];
     }

  sumx = 0.;
  for (j=0; j<num_species-1; j++) 
     {
      x[j] = fv->c[j];
      sumx += x[j];
     }
  x[num_species-1] = 1. - sumx;

  /* MMH
   * For species transport in the SUSPENSION_PM model, there is one 
   * species that represents the particles.  There may be other species
   * involved, though, and they care about the "averaged" density, which
   * is the fluid density here.
   *
   * I don't see where this is used any where.  The functions that
   * this function calls seem to be able to figure out the correct
   * density on their own...
   */
  rho  = density(d_rho, time);

  if( mp->Spwt_funcModel == GALERKIN)
    {
      supg = 0.;
    }
  else if( mp->Spwt_funcModel == SUPG)
    {
      supg = mp->Spwt_func;
    }


  if(supg!=0.)
    {
      h_elem = 0.;
      for ( p=0; p<dim; p++)
	{
	  h_elem += vcent[p]*vcent[p]*h[p];
	}
      h_elem = sqrt(h_elem)/2.;
    }
  /* end Petrov-Galerkin addition */

  /************************************************************************/
  /*                       START OF SPECIES ASSEMBLE                      */
  /************************************************************************/
  /*
   *  Initialize the Species_Conservation_Terms temporary structure
   *  before filling it up
   */
  zero_structure(&s_terms, sizeof(struct Species_Conservation_Terms), 1);

/*  if (mp->PorousMediaType == CONTINUOUS)
    { */
      /*
       * ---------------------------------------------------------------------
       * ----- Calculate the capacity, flux, convection, and source terms for
       *       a continuous medium 
       * ---------------------------------------------------------------------
       */
      err = get_continuous_species_terms(&s_terms, time, tt, dt, h);
      EH(err,"problem in getting the species terms");
      
      /*    } */   /* end of if CONTINUOUS */

  /*
   * Residuals_________________________________________________________________
   */
  if (af->Assemble_Residual)
    {
      var = MASS_FRACTION;
      /* 
       *  Store the species eqn type (which is keyed to the Variable
       *  type in a temporary variable).
       */
      species_eqn_type = mp->Species_Var_Type;

      /*
       *   START loop over species equations. The outer loop is over
       *   the species number
       */
      for ( w=0; w<pd->Num_Species_Eqn; w++)
	{
	  /*
	   *  Calculate the coef_rho term based upon the value of
	   *  species_eqn_type.
	   *     This coefficient indicates whether the time and advection
	   *     term should be multipled by the density or maybe the
	   *     concentration. 
	   */
	  coeff_rho = 1.0;
	  if ((species_eqn_type == SPECIES_MASS_FRACTION) ||
	      (species_eqn_type == SPECIES_MOLE_FRACTION)) {
	    coeff_rho = rho;
	  }

          if (mp->SpeciesSourceModel[w]  == ELECTRODE_KINETICS ||
              mp->SpeciesSourceModel[w]  == ION_REACTIONS) /*  RSL 3/19/01  */
	     {
	      sumxm = 0.;
	      for (j=0; j<num_species; j++)
	         {
	          sumxm += x[j] * M[j];
	         }
	      if (mp->PorosityModel == CONSTANT)
	         {
	          epsilon = mp->porosity;
	         }
	      else if (mp->PorosityModel == THERMAL_BATTERY)
	         {
	          epsilon = mp->u_porosity[0];
	         }
	      else
	         {
	          EH(-1, "invalid porosity model");
	         }
              /* rho = density(d_rho); /\*  RSL 6/22/02  *\/ */
	      small_c = rho/sumxm;
              coeff_rho = small_c; /*  RSL 9/27/01  */
	     }
		  
	  /*
	   *  Loop over the degrees of freedom in the element corresponding
	   *  to species unknowns. These are at different nodes in the
	   *  element, usually. However, there can be more than one set
	   *  of degrees of freedom at each node. Note, this
	   *  step doesn't depend upon the species equation number, so
	   *  we might think about exchanging the order of the loops!
	   */
	  for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
	    ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];
	    if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
	      /*
	       *  Here is where we figure out whether the row is to placed in
	       *  the normal spot (e.g., ii = i), or whether a boundary condition
	       *  require that the volumetric contribution be stuck in another
	       *  ldof pertaining to the same variable type.
	       */
	      ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

		  phi_i = bf[eqn]->phi[i];
		  mass = 0.;
		  if ( pd->TimeIntegration != STEADY )
		    {
		      if ( pd->e[pg->imtrx][eqn] & T_MASS ) {
			  mass  = coeff_rho * s_terms.Y_dot[w];

                          if (mp->SpeciesSourceModel[w]  == ELECTRODE_KINETICS ||
                              mp->SpeciesSourceModel[w]  == ION_REACTIONS) /*  RSL 3/19/01  */
			  {
			    mass = s_terms.Y_dot[w];
			    sumxdot = 0.0;
			    sumxdotm = 0.0;
			    for (j = 0; j < num_species - 1; j++) {
			      sumxdotm += s_terms.Y_dot[j] * M[j];
			      sumxdot  += s_terms.Y_dot[j];
			    }
			    sumxdotm -= sumxdot * M[num_species-1];
			    mass -= x[w]*sumxdotm/sumxm;
			    mass *= (epsilon*small_c);
			  }

			  mass *= - phi_i * h3 * det_J * wt;
			  mass *= pd->etm[pg->imtrx][eqn][LOG2_MASS];
			}
		    }
		  
		  /* only use Petrov Galerkin on advective term - if required */
		  wt_func = phi_i;
		  /* add Petrov-Galerkin terms as necessary */
		  if (supg != 0.0)
		    {
		      for(p=0; p<dim; p++)
			{
			  wt_func += supg * h_elem * fv->v[p] * bf[eqn]->grad_phi[i][p];
			}
		    }
		  
		  /*
		   *   Advection is velocity times gradient of the species unknown
		   *   variable, possibly multiplied by a density or total
		   *   concentration, depending upon the species variable
		   *   type.
		   */
		  advection = 0.0;
		  if ( pd->e[pg->imtrx][eqn] & T_ADVECTION )
		    {
		      advection_a = 0.0;
		      for ( p=0; p<VIM; p++)
			{
			  advection_a += s_terms.conv_flux[w][p];
			}

                      if (mp->SpeciesSourceModel[w]  == ELECTRODE_KINETICS ||
                          mp->SpeciesSourceModel[w]  == ION_REACTIONS) /*  RSL 3/19/01  */
		         {
		          psum = 0.;
		          for ( p=0; p<VIM; p++)
		             {
		              sumflux = 0.;
		              sumfluxm = 0.;
		              for (j=0; j<num_species-1; j++)
		                 {
		                  sumfluxm += s_terms.conv_flux[j][p] * M[j];
		                  sumflux  += s_terms.conv_flux[j][p];
		                 }
		              sumfluxm -= sumflux * M[num_species-1];
		              psum += sumfluxm;
		             }
		          advection_a -= x[w]*psum/sumxm;
		          advection_a *= small_c;
		          advection_a *= (-wt_func);
		         }
		      else
		         {
		          advection_a *= - coeff_rho * wt_func;
		         }

		      advection_b = 0.;
		      if(taylor_galerkin[w])
			{
			  for ( p=0; p<VIM; p++)
			    {
			      advection_b += s_terms.taylor_flux[w][p];
			    }
			  advection_b *= - coeff_rho * s_terms.taylor_flux_wt[i]*dt/2.;
			}
		      advection = advection_a + advection_b;
		      advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)] * det_J * h3 * wt 
		       	           * mp->AdvectiveScaling[w];
		    }
		  
		  /*
		   * the diffusion term contains all the fluxes that are in
		   * a divergence in the mass conservation equation and are 
		   * integrated by parts in finite elements
		   *   Units of s_terms.diff_flux:
		   *      SPECIES_MASS_FRACTION:    gm cm-2 sec-1
		   *      SPECIES_MOLE_FRACTION:    mol cm-2 sec-1
		   */
		  diffusion = 0.;
		  if ( pd->e[pg->imtrx][eqn] & T_DIFFUSION )
		    {
		      for ( p=0; p<VIM; p++)
			{
			  grad_phi_i[p] = bf[eqn]->grad_phi[i] [p];
			}
		      
		      for ( p=0; p<VIM; p++)
			{
			  diffusion += grad_phi_i[p] * s_terms.diff_flux[w][p];
			}

                      if ((mp->SpeciesSourceModel[w]  == ELECTRODE_KINETICS ||
                           mp->SpeciesSourceModel[w]  == ION_REACTIONS)  &&
                          cr->MassFluxModel == STEFAN_MAXWELL_VOLUME) /*  RSL 3/19/01  */
		         {
		          diffusion *= small_c;
		         }

		      diffusion *= h3 * det_J * wt;
		      diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
		    }
		  /*
		   * HKM -> Note the addition of a species molecular weight
		   *        term is currently done in the source term
		   *        calculation section of the code
		   *   Units ofs_terms.MassSource:
		   *      SPECIES_MASS_FRACTION:    gm cm-3 sec-1
		   *      SPECIES_MOLE_FRACTION:    mol cm-3 sec-1
		   */
		  source = 0.;
		  if ( pd->e[pg->imtrx][eqn] & T_SOURCE )
		    {
		      source = s_terms.MassSource[w];

                      if (mp->SpeciesSourceModel[w]  == ELECTRODE_KINETICS ||
                          mp->SpeciesSourceModel[w]  == ION_REACTIONS) /*  RSL 3/19/01  */
		         {
		          sumrm = 0.;
		          for (j=0; j<num_species; j++)
		             {
		              sumrm += s_terms.MassSource[j] * M[j];
		             }
		          source -= x[w]*sumrm/sumxm;
		         }

		      source *= phi_i * h3 * det_J * wt; 
		      source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
		    }
		    
                  residual = mass + advection +  diffusion + source;
		  
		  /*
		   * Path dependence is residual times derivative of H wrt to F
		   */
		  var = FILL;
                  pvar = upd->vp[pg->imtrx][var];
                  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++ )
                    {
                      lec->J[MAX_PROB_VAR + w][pvar][ii][j] += lsi->d_H_dF[j] * residual * sign;
                    }
		}   /* if active_dofs */
	      
	    } /* end of loop over equations */
	} /* end of loop over species */
    } /* end of assemble residuals */
  return 0;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void 
mass_flux_surf_mtc(dbl mass_flux[MAX_CONC],
		   dbl d_mass_flux[MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC],
		   double Tsurf_quad, /* Value of temperature at the surface */
		   double Ysurf_quad[MAX_CONC], /* concentration at surface 
						 * quadrature point          */
		   int wspec,              /* subvariable species number     */
		   double mass_tran_coeff, /* Mass transfer coefficient 
					    * (cgs?? MKS units)              */
		   double Y_c)	           /* bath concentration             */

/******************************************************************************
*
*  Function which calculates the mass flux rate for convective
*  mass transfer using a mass transfer coefficient.
*
******************************************************************************/
{
  int w;
  /* HARDWIRE a linear increase in MTC from zero to mass_tran_coeff
   * along free surface boundary */ 
  /*      mtc = mass_tran_coeff * (fv->x[0] - 5)/20.; */
  if (af->Assemble_Jacobian) {
    /*
     * Currently there is no sensitivity to temperature 
     * or to other concentrations
     */
    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      d_mass_flux[wspec][MAX_VARIABLE_TYPES + w] = 0.;
    }
    d_mass_flux[wspec][MAX_VARIABLE_TYPES + wspec] = mass_tran_coeff;
    d_mass_flux[wspec][TEMPERATURE] = 0.;
  }
  mass_flux[wspec] = mass_tran_coeff * (Ysurf_quad[wspec] - Y_c);
  return;
} /* END of routine mass_flux_surf_mtc                           */ 
/****************************************************************************/

void 
mass_flux_surf_BV(dbl mass_flux[MAX_CONC],
		  dbl d_mass_flux[MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC],
		  int wspec,                    /* species number                     */
		  double nu,                    /* stoichiometric coefficient         */  
		  double k,                     /* kinetic rate constant              */  
		  double beta,                  /* reaction order                     */
		  double alphaa,                /* anodic direction transfer coef.    */
		  double alphac,                /* cathodic direction transfer coef.  */
		  double V,                     /* electrode potential                */
		  double U0,                    /* open circuit electrolyte potential */
		  double T)	                /* electrolyte solution temperature   */

/******************************************************************************
*
*  A function that calculates the mass flux rate of a given species, which is equal to   
*  the heterogeneous or surface reaction rate given by Butler-Volmer kinetics. 
*  This routine was created by cloning after mass_flux_surf_mtc.
*
*  Ken S. Chen (11/2000)
*
*  Revised by K. S. Chen on 10/1/2001. 
*
******************************************************************************/
{
  int w;

  const double R = 8.314;         /* Universal gas constant in units of J/mole K */
  const double F = 96487.0;       /* Faraday's constant in units of C/equiv. */
  double FRT;                     /* product of F/R/T */
  double PHI;                     /* electrical potential in electrolyte phase */ 
  double conc, conc1, grpa, grpc;

  if(nAC) 
   {                              /* if augmenting condition card is active, then  */
     mp->electrode_potential = V; /* set electrode potential to be used in current_BV_surf */
   }                      
  
  FRT = F/R/T; 
  PHI = fv->V;
  conc = pow(fv->c[wspec], beta); 
  conc1 = beta*pow(fv->c[wspec], beta-1.0);
  grpa = alphaa*FRT*(V-PHI-U0); 
  grpc = alphac*FRT*(V-PHI-U0);   

  mass_flux[wspec] = nu*k*conc*(exp(grpa) - exp(-grpc));

  if (af->Assemble_Jacobian) 
   {
     for (w = 0; w < pd->Num_Species_Eqn; w++) 
       {
        d_mass_flux[wspec][MAX_VARIABLE_TYPES + w] = 0.;   /* no dependence on other species */ 
       }
     d_mass_flux[wspec][MAX_VARIABLE_TYPES + wspec] = nu*k*conc1*(exp(grpa) - exp(-grpc));
     d_mass_flux[wspec][TEMPERATURE] = nu*k*conc*(FRT/T)*(-V+PHI+U0)*(alphaa*exp(grpa)+alphac*exp(-grpc)); 
     d_mass_flux[wspec][VOLTAGE] = -nu*k*conc*FRT*(alphaa*exp(grpa)+alphac*exp(-grpc));  
   } 


  return;
} /* END of routine mass_flux_surf_BV                           */ 
/****************************************************************************/

void
mass_flux_surf_HOR(dbl mass_flux[MAX_CONC],
                   dbl d_mass_flux[MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC],
                   int wspec,      /* species number                          */
                   double ai0,     /* product of interfacial area by          */
                                   /* exchange current density (A/cm^3)       */
                   double H,       /* thickness of catalyst layer (cm)        */
                   double cref,    /* ref. concentration (moles/cm^3)         */
                   double alphaa,  /* anodic direction transfer coefficient   */
                   double alphac,  /* cathodic direction transfer coefficient */
                   double T,       /* cell temperature (K)                    */
                   double U0,      /* open-circuit potential (V)              */
                   double beta,    /* reaction order                          */
                   double n,       /* number of electrons involved in rxn     */
                   double V)       /* electrode potential                     */

/******************************************************************************
*
*  A function that calculates the mass flux rate of a given species,
*  which is equal to the heterogeneous or surface reaction rate
*  given by a linearized kinetic model for an electrochemical reaction
*  such as the hydrogen oxidation reaction in PEM fuel cells
*  (cf. Chen and Hickner 2006). This routine was created
*  by cloning after mass_flux_surf_BV.
*
*  Ken S. Chen (1/2006)
*  Modified: 5/22/2006 by KSC.
*
******************************************************************************/
{
  int w;                          /* species index                               */
  const double R = 8.314;         /* Universal gas constant in units of J/mole K */
  double nRT;                     /* product of n*R*T                            */
  double PHI;                     /* electrical potential in electrolyte phase   */
  double c;                       /* species molar concentration (moles/cm^3)    */
  double cratio;

  nRT = n*R*T;
  PHI = fv->V;
  c = fv->c[wspec];
  if(c < 0.0) c = 1.0e-10;
  if(c == 0) c = cref;
  cratio = pow(c/cref, beta);

  mass_flux[wspec] = (ai0*H/nRT)*cratio*(alphaa+alphac)*(V-PHI-U0);

  if (af->Assemble_Jacobian)
   {
     for (w = 0; w < pd->Num_Species_Eqn; w++)
       {
        d_mass_flux[wspec][MAX_VARIABLE_TYPES + w] = 0.;
       }

     d_mass_flux[wspec][MAX_VARIABLE_TYPES + wspec] =
                                           mass_flux[wspec]*(beta/c);
     d_mass_flux[wspec][TEMPERATURE] = -mass_flux[wspec]/T;
     d_mass_flux[wspec][VOLTAGE] =
                                 -(ai0*H/nRT)*cratio*(alphaa+alphac);
   }

  return;
} /* END of routine mass_flux_surf_HOR                           */ 
/****************************************************************************/

void
mass_flux_surf_ORR(dbl mass_flux[MAX_CONC],
                   dbl d_mass_flux[MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC],
                   int wspec,      /* species number                            */
                   double ai0,     /* product of interfacial area by            */
                                   /* exchange current density (A/cm^3)         */
                   double H,       /* thickness of catalyst layer (cm)          */
                   double cref,    /* species ref. concentration (moles/cm^3)   */
                   double alphac,  /* cathodic direction transfer coefficient   */
                   double T,       /* cell temperature (K)                      */
                   double V,       /* cel voltage (V)                           */
                   double U0,      /* open-circuit potential (V)                */
                   double beta,    /* reaction order                            */
                   double n)       /* number of electrons involved in rxn       */

/******************************************************************************
*
*  A function that calculates the mass flux rate of a given species,
*  which is equal to the heterogeneous or surface reaction rate
*  given by a Tafel kinetic model for an electrochemical reaction
*  such as the oxygen reduction reaction (Chen and Hickner 2006).
*  reaction (Chen and Hickner 2006). This routine was created
*  by cloning after mass_flux_surf_HOR.
*
*  Ken S. Chen (1/2006)
*  Modified: 5/22/2006 by KSC.
*
******************************************************************************/
{
  int w;                          /* species index                               */
  const double F = 96487;         /* Faraday's constant in units of C/mole       */
  const double R = 8.314;         /* Universal gas constant in units of J/mole-K */
  double FRT;                     /* F/R/T                                       */
  double PHI;                     /* electrical potential in electrolyte phase   */
  double c;                       /* molar concentration of O2 (moles/cm^3)      */
  double cratio;
  double grp;

  FRT = F/R/T;
  PHI = fv->V;
  c = fv->c[wspec];
  if(c == 0) c = cref;
  if(c < 0.0) c = 1.0e-10;
  cratio = pow(c/cref, beta);
  grp = alphac*FRT*(V-PHI-U0);

  mass_flux[wspec] = (ai0*H/(n*F))*cratio*exp(-grp);

  if (af->Assemble_Jacobian)
   {
     for (w = 0; w < pd->Num_Species_Eqn; w++)
       {
        d_mass_flux[wspec][MAX_VARIABLE_TYPES + w] = 0.;
       }

     d_mass_flux[wspec][MAX_VARIABLE_TYPES + wspec] =
                                      mass_flux[wspec]*(beta/c);
     d_mass_flux[wspec][TEMPERATURE] = mass_flux[wspec]*(grp/T);
     d_mass_flux[wspec][VOLTAGE] = mass_flux[wspec]*alphac*FRT;
   }

  return;

} /* END of routine mass_flux_surf_ORR                           */ 
/****************************************************************************/

void
mass_flux_surf_H2O_ANODE(dbl mass_flux[MAX_CONC],
                         dbl d_mass_flux[MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC],
                         int wspec,      /* species number                          */
                         double ai0,     /* product of interfacial area by          */
                                         /* exchange current density (A/cm^3)       */
                         double Ha,      /* thickness of anode catalyst layer (cm)  */
                         double cH2ref,  /* ref. concentration for H2 (moles/cm^3)  */
                         double alphaa,  /* anodic direction transfer coefficient   */
                         double alphac,  /* cathodic direction transfer coefficient */
                         double T,       /* cell temperature (K)                    */
                         double U0a,     /* Open-circuit potential for HOR (V)      */
                         double nd)      /* electro-osmatic drag coefficient        */

/******************************************************************************
*
*  A function that calculates the mass flux rate of the H2O species,
*  which is due to electro-osmatic drag. Briefly,
*  
*  rH2O = nd (i/F) = nd rH2   where rH2 is the surface rate of H2
*  consumption due to the hydrogen oxidation reaction
*
*  Ken S. Chen (1/2006)
*
******************************************************************************/
{
  int w;                          /* species index                               */
  const double R = 8.314;         /* Universal gas constant in units of J/mole K */
  double RT;                      /* product of R*T                              */
  double PHI;                     /* electrical potential in electrolyte phase   */ 
  double cH2;                     /* molar concentration of H2 (moles/cm^3)      */
  double cratio;

  RT = R*T;
  PHI = fv->V;
  cH2 = fv->c[wspec];
  if(cH2 == 0) cH2 = cH2ref;
  cratio = pow(cH2/cH2ref, 0.5);

  mass_flux[wspec] = -ai0*Ha*cratio*0.5*(alphaa+alphac)*((PHI+U0a)/RT)*nd;

  if (af->Assemble_Jacobian) 
   {
     for (w = 0; w < pd->Num_Species_Eqn; w++) 
       {
        d_mass_flux[wspec][MAX_VARIABLE_TYPES + w] = 0.;
       }

     d_mass_flux[wspec][MAX_VARIABLE_TYPES + wspec] =
                                           mass_flux[wspec]/(2.0*cH2);
     d_mass_flux[wspec][TEMPERATURE] = -mass_flux[wspec]/T;
     d_mass_flux[wspec][VOLTAGE] =
                                 -ai0*Ha*cratio*0.5*(alphaa+alphac)*nd/RT;

   }

  return;
} /* END of routine mass_flux_surf_H2O_ANODE                           */ 
/****************************************************************************/

void
mass_flux_surf_H2O_CATHODE(dbl mass_flux[MAX_CONC],
                           dbl d_mass_flux[MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC],
                           int wspec,      /* species number                           */
                           double ai0,     /* product of interfacial area by           */
                                           /* exchange current density (A/cm^3)        */
                           double Hc,      /* thickness of cathode catalyst layer (cm) */
                           double cO2ref,  /* ref. concentration for O2 (moles/cm^3)   */
                           double alphac,  /* cathodic direction transfer coefficient  */
                           double T,       /* cell temperature (K)                     */
                           double V,       /* cel voltage (V)                          */
                           double U0c,     /* Open-circuit potential for ORR (V)       */
                           double nd)      /* electro-osmatic drag coefficient         */

/******************************************************************************
*
*  A function that calculates the mass flux rate of the H2O species,
*  which is due to electro-osmatic drag and oxygen reduction reaction:
*
*  rH2O = nd (i/F) + 2 rO2 = (nd+2) rO2   where rO2 is the surface rate of O2
*  consumption due to the oxygen reduction reaction and nd is the
*  electro-osmatic drag coefficient
*
*  Ken S. Chen (1/2006)
*
******************************************************************************/
{
  int w;                          /* species index                               */
  const double F = 96487;         /* Faraday's constant in units of C/mole       */
  const double R = 8.314;         /* Universal gas constant in units of J/mole-K */
  double FRT;                     /* F/R/T                                       */
  double PHI;                     /* electrical potential in electrolyte phase   */
  double cO2;                     /* molar concentration of O2 (moles/cm^3)      */
  double cratio;
  double grp;

  FRT = F/R/T;
  PHI = fv->V;
  cO2 = fv->c[wspec];
  if(cO2 < 0.0) cO2 = 1.0e-10;
  if(cO2 == 0) cO2 = cO2ref;
  cratio = cO2/cO2ref;
  grp = alphac*FRT*(V-PHI-U0c);

  mass_flux[wspec] = -(nd+2.0)*(ai0*Hc/(4.0*F))*cratio*exp(-grp);

  if (af->Assemble_Jacobian)
   {
     for (w = 0; w < pd->Num_Species_Eqn; w++)
       {
        d_mass_flux[wspec][MAX_VARIABLE_TYPES + w] = 0.;
       }

     d_mass_flux[wspec][MAX_VARIABLE_TYPES + wspec] = mass_flux[wspec]/cO2;
     d_mass_flux[wspec][TEMPERATURE] = mass_flux[wspec]*(grp/T);
     d_mass_flux[wspec][VOLTAGE] = mass_flux[wspec]*alphac*FRT;
   }

  return;

} /* END of routine mass_flux_surf_H2O_CATHODE                           */
/****************************************************************************/

void 
mass_flux_surf_SULFIDATION (dbl mass_flux[MAX_CONC],
		            dbl d_mass_flux[MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC],
                            int mode,      /* key word for sulfidation kinetic model */
		            int wspec,     /* species number             */
		            double nu,     /* stoichiometric coefficient */  
		            double k1,     /* forward rate constant      */  
		            double E1,     /* forward activation energy  */
		            double kn1,    /* backward rate constant     */  
		            double En1,    /* backward activation energy */
		            double T,      /* Temperature                */
		            double c_H2S,  /* bulk concentration of H2S  */
		            double c_O2)   /* bulk concentration of O2   */

/******************************************************************************
*
*  A function that calculates the mass flux rate of a given species, which is    
*  equl to the heterogeneous or surface reaction rate given by the various 
*  sulfidation kinetic models as selected by the user.
*
*  SOLID_DIFFUSION_SIMPLIFIED:
*    Simplified solid-diffusion-controlled with Cu as the diffusing species:
*   
*    r = k exp(-E/R/T) c_H2S c_Cu       ( c_H2S is fixed at its bulk value ) 
*
*  SOLID_DIFFUSION_ELECTRONEUTRALITY:
*    Solid-diffusion-controlled with Cu vacancies as diffusing species and
*    the approximation of electroneutrality ( i.e., c_V = c_h ):
*
*    r = k1 exp(-E1/R/T) c_H2S c_O2**0.5 - kn1 exp(-En1/R/T) c_V**4 
*
*    ( c_H2S and c_O2 are taken to be fixed at its bulk values )
*
*  SOLID_DIFFUSION_ELECTRONEUTRALITY_LINEAR:
*    Solid-diffusion-controlled with Cu vacancies as diffusing species and
*    the approximation of electroneutrality ( i.e., c_V = c_h )
*    using a reduced reverse reaction rate dependence on the vacancies
*
*    r = k1 exp(-E1/R/T) c_H2S c_O2**0.5 - kn1 exp(-En1/R/T) c_V**2
*
*    ( c_H2S and c_O2 are taken to be fixed at its bulk values )
*
*  SOLID_DIFFUSION:
*    Solid-diffusion-controlled with Cu vacancies as diffusing species and
*    c_V being different from c_h near the Cu/Cu2S and gas/Cu2S interfaces:
*
*    r = k1 exp(-E1/R/T) c_H2S c_O2**0.5 - kn1 exp(-En1/R/T) c_V**2 c_h**2 
*
*    ( c_H2S and c_O2 are taken to be fixed at its bulk values )
*
*  GAS_DIFFUSION:
*    Gas-diffusion-controlled with H2S bnd O2 being the diffusing species:
*
*    r = k1 exp(-E1/R/T) c_H2S c_O2**0.5     (c_V is taken to be zero)
*
*  ANNIHILATION_ELECTRONEUTRALITY: 
*    Annihilation reaction at the Cu/Cu2S interface:
*
*    r = k2 exp(-E2/R/T) c_V c_V     (with electroneutrality approximation)           
*
*  ANNIHILATION:
*    Annihilation reaction at the Cu/Cu2S interface:
*    r = k2 exp(-E2/R/T) c_V c_h     (general case: c_h != c_V)           
*
*  This routine was created by cloning after mass_flux_surf_BV.
*  Ken S. Chen (3/2002)
******************************************************************************/
{
  int w;
  const double R = 1.987; /* Universal gas constant in units of cal/mole K */
  double c_Cu;            /* concentration of Cu atoms       */
  double c_V;             /* concentration of Cu vacancies   */
  double c_h;             /* concentration of electron holes */

  if(mode == SOLID_DIFFUSION_SIMPLIFIED)
    {
      c_Cu = fv->c[0];    /* Cu is the diffusing species */
      mass_flux[wspec] = nu*k1*exp(-E1/R/T)*c_H2S*c_Cu;
    } 
  else if(mode == SOLID_DIFFUSION_ELECTRONEUTRALITY)
    {
      c_V = fv->c[wspec];  /* Cu vacancies and electron holes are diffusing species */
                           /* c_h = c_V due to electroneutrality approximation */
      mass_flux[wspec] = k1*exp(-E1/R/T)*c_H2S*sqrt(c_O2) - 
                         kn1*exp(-En1/R/T)*c_V*c_V*c_V*c_V; 
      mass_flux[wspec] *= nu;
    }
  else if(mode == SOLID_DIFFUSION_ELECTRONEUTRALITY_LINEAR)
    {
      c_V = fv->c[wspec];  /* Cu vacancies and electron holes are diffusing species */
                           /* c_h = c_V due to electroneutrality approximation */
      mass_flux[wspec] = k1*exp(-E1/R/T)*c_H2S*sqrt(c_O2) - 
                         kn1*exp(-En1/R/T)*c_V*c_V; 
      mass_flux[wspec] *= nu;
    }
  else if(mode == SOLID_DIFFUSION)
    {
      c_V = fv->c[0];     /* Cu vacancy is the 1st diffusing species */
      c_h = fv->c[1];     /* Electron hole is the 2nd diffusion species */
      mass_flux[wspec] = k1*exp(-E1/R/T)*c_H2S*sqrt(c_O2) - 
                         kn1*exp(-En1/R/T)*c_V*c_V*c_h*c_h; 
      mass_flux[wspec] *= nu;
    }
  else if(mode == GAS_DIFFUSION)
    { 
      c_H2S = fv->c[0];   /* H2S is the 1st diffusing species  */
      c_O2  = fv->c[1];   /* O2 is the 2nd diffusion species   */
      mass_flux[wspec] = nu*k1*exp(-E1/R/T)*c_H2S*sqrt(c_O2);
    }
  else if(mode == FULL)
    {
      fprintf(stderr, "The full model has not yet implemented - awaits future efforts\n");
      exit(1);
    }
  else if(mode == ANNIHILATION_ELECTRONEUTRALITY)
    {
      /* c_V = fv->c[0]; */    /* Cu is the diffusing species */
      c_V = fv->c[wspec];      
      mass_flux[wspec] = nu*k1*exp(-E1/R/T)*c_V*c_V;
    } 
  else if(mode == ANNIHILATION)
    {
      fprintf(stderr, "The ANNIHILATION model without electroneutrality has not yet implemented\n");
      exit(1);
    }

  if (af->Assemble_Jacobian) 
   {
     /*
      *  no dependency on other species 
      */ 
     for (w = 0; w < pd->Num_Species_Eqn; w++)  {
       d_mass_flux[wspec][MAX_VARIABLE_TYPES + w] = 0.0;
     }
     /*
      * There is no explicit dependence on the voltage drop across the interface
      */
     d_mass_flux[wspec][VOLTAGE] = 0.0;
     if (mode == SOLID_DIFFUSION_SIMPLIFIED)
       {
        d_mass_flux[wspec][MAX_VARIABLE_TYPES + wspec] = nu*k1*exp(-E1/R/T)*c_H2S;
        c_Cu = fv->c[0];    
        d_mass_flux[wspec][TEMPERATURE] = nu*k1*(E1/R/T/T)*exp(-E1/R/T)*c_H2S*c_Cu; 
       }
     else if (mode == SOLID_DIFFUSION_ELECTRONEUTRALITY)
       {
        c_V = fv->c[wspec];  
        d_mass_flux[wspec][MAX_VARIABLE_TYPES + wspec] = -kn1*exp(-En1/R/T)*4.0*c_V*c_V*c_V;
        d_mass_flux[wspec][MAX_VARIABLE_TYPES + wspec] *= nu;
        d_mass_flux[wspec][TEMPERATURE] = k1*(E1/R/T/T)*exp(-E1/R/T)*c_H2S*sqrt(c_O2) -
                                          kn1*(En1/R/T/T)*exp(-En1/R/T)*c_V*c_V*c_V*c_V;
        d_mass_flux[wspec][TEMPERATURE] *= nu; 
       }
     else if (mode == SOLID_DIFFUSION_ELECTRONEUTRALITY_LINEAR)
       {
        c_V = fv->c[wspec];  
        d_mass_flux[wspec][MAX_VARIABLE_TYPES + wspec] = -kn1*exp(-En1/R/T)*2.0*c_V;
        d_mass_flux[wspec][MAX_VARIABLE_TYPES + wspec] *= nu;
        d_mass_flux[wspec][TEMPERATURE] = k1*(E1/R/T/T)*exp(-E1/R/T)*c_H2S*sqrt(c_O2) -
                                          kn1*(En1/R/T/T)*exp(-En1/R/T)*c_V*c_V;
        d_mass_flux[wspec][TEMPERATURE] *= nu; 
       }
     else if (mode == SOLID_DIFFUSION)
       {
         c_V = fv->c[0];
         c_h = fv->c[1];
         d_mass_flux[wspec][MAX_VARIABLE_TYPES + wspec] = -kn1*exp(-En1/R/T)*2.0*fv->c[wspec];
         if (wspec == 0)
           {
             d_mass_flux[wspec][MAX_VARIABLE_TYPES + wspec] *= c_h*c_h;
           }
         else if (wspec == 1)
           {
             d_mass_flux[wspec][MAX_VARIABLE_TYPES + wspec] *= c_V*c_V;
           }
         d_mass_flux[wspec][MAX_VARIABLE_TYPES + wspec] *= nu;
         d_mass_flux[wspec][TEMPERATURE] = k1*(E1/R/T/T)*exp(-E1/R/T)*c_H2S*sqrt(c_O2) -
                                           kn1*(En1/R/T/T)*exp(-En1/R/T)*c_V*c_V*c_h*c_h;
         d_mass_flux[wspec][TEMPERATURE] *= nu; 
       }
     else if (mode == GAS_DIFFUSION)
       {
         c_H2S = fv->c[0];
         c_O2  = fv->c[1];
         d_mass_flux[wspec][MAX_VARIABLE_TYPES + wspec]  = nu*k1*exp(-E1/R/T)*fv->c[wspec];
         if (wspec == 0)
           {
             d_mass_flux[wspec][MAX_VARIABLE_TYPES + wspec] *= c_O2;
           }
         else if (wspec == 1)
           {
             d_mass_flux[wspec][MAX_VARIABLE_TYPES + wspec] *= c_H2S;
           }
         d_mass_flux[wspec][TEMPERATURE] = nu*k1*(E1/R/T/T)*exp(-E1/R/T)*c_H2S*c_O2; 
       }
     else if (mode == FULL)
       {
         fprintf(stderr, "The full model has not yet implemented - awaits future efforts\n");
         exit(1);
       }
     else if (mode == ANNIHILATION_ELECTRONEUTRALITY)
       {
        c_V = fv->c[wspec];      
        d_mass_flux[wspec][MAX_VARIABLE_TYPES + wspec] = nu*k1*exp(-E1/R/T)*2.0*c_V;
        d_mass_flux[wspec][TEMPERATURE] = nu*k1*(E1/R/T/T)*exp(-E1/R/T)*c_V*c_V; 
       }
     else if (mode == ANNIHILATION)
       {
         fprintf(stderr, 
		 "The ANNIHILATION model without electroneutrality"
		 " has not yet implemented\n");
         exit(1);
       }
   }
  return;
} /* END of routine mass_flux_surf_SULFIDATION                           */ 
/****************************************************************************/

void
mass_flux_surf_BV2(dbl time,        /* current time value                    */
		   dbl mass_flux[MAX_CONC],
                   dbl d_mass_flux[MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC],
                   int wspec,       /* species number                        */
                   int flag,        /* output option -- mass flux or current */
                   double k,        /* rate constant                         */
                   double PHI_E,    /* electrode potential                   */
                   double alphaa,   /* anodic transfer coefficient           */
                   double alphac,   /* cathodic transfer coefficient         */
                   double U0,       /* standard state open circuit potential */
                   double T)        /* electrolyte temperature               */

/**********************************************************************************
*
*  An alternate to the routine mass_flux_surf_BV, using a less conventional but
*  more convenient form of the Butler-Volmer rate law.  The reaction order beta
*  used in mass_flux_surf_BV is not needed here.  In addition, the U0 in this
*  version is the open circuit potential in the standard state, not at the
*  prevailing concentration, so it is at most a function of temperature.  The
*  standard state concentration is taken to be 0.001 mol/cm^3, i.e., unit molarity.
*  The computed mass flux is the rate of the electrode reaction in the anodic
*  direction, so a negative rate constant must be used if the rate in the cathodic
*  direction is needed.  Since the activation energy of the rate constant (or
*  exchange current density) is not specified, temperature derivatives of the mass
*  flux are not calculated here.
*
*  In this version, the ionic species is assumed to be consumed by reduction or
*  produced by oxidation, i.e., it is a reactant in the cathodic direction.  This
*  covers the usual thermal battery and electroplating applications.
*
*  RSL 1/10/01
*
*  Revised 10/1/01 to fine-tune calculation of Jacobian terms
*
**********************************************************************************/
{
  int w, j, store;

  const double R = 8.314;         /* universal gas constant in units of J/mole K */
  const double F = 96487.;        /* Faraday's constant in units of C/equiv      */
  double FRT;                     /* F/(R*T)                                     */
  double PHI_S;                   /* electrical potential in solution phase      */
  double conc, grpa, grpc, diff_phi, kgc, derivative;
  dbl c, rho, M_mix;
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;  /* density dependence */
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;

  FRT = F/R/T;
  PHI_S = fv->V;
  M_mix = 0.;
  for (j=0; j<pd->Num_Species; j++)
     {
      M_mix += fv->c[j] * mp->molecular_weight[j];
     }
  rho = density(d_rho, time); /*  RSL 6/22/02  */
  c = rho/M_mix;
  conc = c*fv->c[wspec];
  diff_phi = FRT*(PHI_E - PHI_S - U0);
  grpa = exp( alphaa*diff_phi);
  grpc = exp(-alphac*diff_phi);

  if (flag == 0)
     store = wspec;
  else
     store = MAX_CONC - 1;

  mass_flux[store] = k*(grpa - (conc/0.001)*grpc);

  if (af->Assemble_Jacobian)
    {
     kgc = k*grpc/0.001;
     for (w = 0; w < pd->Num_Species_Eqn; w++)
       {
        for ( j=0; j<ei[pg->imtrx]->dof[MASS_FRACTION]; j++)
           {
            if ( bf[MASS_FRACTION]->phi[j] > 0.0 ) break;
           }
        derivative = d_rho->C[w][j]/bf[MASS_FRACTION]->phi[j];
        d_mass_flux[store][MAX_VARIABLE_TYPES + w] = -kgc*fv->c[wspec]*(derivative -
            c*(mp->molecular_weight[w] - mp->molecular_weight[pd->Num_Species - 1]))/M_mix
;
       }
     d_mass_flux[store][MAX_VARIABLE_TYPES + wspec] -= kgc*c;
     d_mass_flux[store][TEMPERATURE] = 0.;
     d_mass_flux[store][VOLTAGE] = -k*FRT*(alphaa*grpa+(conc/0.001)*alphac*grpc);
    }

  return;
} /* END of routine mass_flux_surf_BV2                                      */
/****************************************************************************/
/****************************************************************************/

void
mass_flux_surf_NI (dbl mass_flux[MAX_CONC],
                   dbl d_mass_flux[MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC],
                   dbl time,        /* current time                             */
		   int wspec,       /* species number (zero for current)        */
                   int flag,        /* output option -- species flux or current */
                   double PHI_E,    /* electrode potential                      */
                   double T)        /* electrolyte temperature                  */

/**********************************************************************************
*
*  Calculation of either (a) the net consumption rate (equal to the outwardly
*  directed flux) at the electrode surface for an ionic species involved in nickel
*  electroplating, or (b) the total outwardly directed current due to all reactions
*  at the electrode.
*
*  RSL 3/13/01
*
**********************************************************************************/
{
  int w, j, store, n;
  int four;

  const double R = 8.314;         /* universal gas constant in units of J/mole K */
  const double F = 96487.;        /* Faraday's constant in units of C/equiv      */
  double FRT;                     /* F/(R*T)                                     */
  double PHI_S;                   /* electrical potential in solution phase      */
  double conc, grpa, grpc, diff_phi, i00, alphaa, alphac, U00, Q, dQdV, dQdx, dQdy;
  dbl c, rho, M_mix;
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;  /* density dependence */
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;

  if (MAX_CONC < 5) {
    EH(-1, "mass_flux_surf_NI expects MAX_CONC >= 5");
    return;
  }

  four = 4;

  if (flag == 0)
     store = wspec;
  else
     store = MAX_CONC - 1;

  mass_flux[store] = 0.;
  if (af->Assemble_Jacobian)
    {
     for (w = 0; w < pd->Num_Species_Eqn; w++)
       {
        d_mass_flux[store][MAX_VARIABLE_TYPES + w] = 0.;
       }
     d_mass_flux[store][TEMPERATURE] = 0.;
     d_mass_flux[store][VOLTAGE] = 0.;
    }

  if (wspec == 1 || wspec == 5) return;

  FRT = F/R/T;
  PHI_S = fv->V;
  M_mix = 0.;
  for (j=0; j<pd->Num_Species; j++)
     {
      M_mix += fv->c[j] * mp->molecular_weight[j];
     }
  rho = density(d_rho, time); /*  RSL 6/22/02  */
  c = rho/M_mix;

/***  Reaction  H+ + e- = (1/2)H2  ***/

  if (wspec == 2 || flag == 1)
     {
      alphac = 0.11;
      n = 1;
      i00 = 0.289;
      U00 = 0.;
      diff_phi = n*FRT*(PHI_E - PHI_S - U00);
      grpc = exp(-alphac*diff_phi);
      conc = c*fv->c[2];
      Q = (i00/n/F)*grpc*(conc/0.001);
      mass_flux[store] += Q;
      if (af->Assemble_Jacobian)
        {
         dQdx = Q/fv->c[2];
         dQdV = Q*alphac*n*FRT;
         d_mass_flux[store][MAX_VARIABLE_TYPES + 2] += dQdx;
         d_mass_flux[store][VOLTAGE] += dQdV;
        }
     }

/***  Reaction  H2O + e- = (1/2)H2 + OH-  ***/

  if (wspec == 3 || wspec == 6 || flag == 1)
     {
      alphac = 0.4;
      n = 1;
      i00 = 8.80e-09;
      U00 = -0.8280;
      diff_phi = n*FRT*(PHI_E - PHI_S - U00);
      grpc = exp(-alphac*diff_phi);
      Q = (i00/n/F)*grpc;
      if (wspec == 3)
         {
          mass_flux[store] -= Q;
         }
      else
         {
          mass_flux[store] += Q;
         }
      if (af->Assemble_Jacobian)
        {
         dQdV = Q*alphac*n*FRT;
         if (wspec == 3)
            {
             d_mass_flux[store][VOLTAGE] -= dQdV;
            }
         else
            {
             d_mass_flux[store][VOLTAGE] += dQdV;
            }
        }
     }

/***  Reaction  NiOH+ + 2e- = Ni + OH-  ***/

  if (wspec == 3 || wspec == 4 || flag == 1)
     {
      alphac = 0.33;
      n = 2;
      i00 = 1.34;
      U00 = -0.3648;
      diff_phi = n*FRT*(PHI_E - PHI_S - U00);
      grpc = exp(-alphac*diff_phi);
      alphaa = 1 - alphac;
      grpa = exp(alphaa*diff_phi);
      Q = (i00/n/F)*(c/0.001)*(grpc*fv->c[four] - grpa*fv->c[3]);
      if (wspec == 3)
         {
          mass_flux[store] -= Q;
         }
      else if (wspec == 4)
         {
          mass_flux[store] += Q;
         }
      else
         {
          mass_flux[store] += 2.*Q;
         }
      if (af->Assemble_Jacobian)
        {
         dQdx = (i00/n/F)*(c/0.001)*grpc;
         dQdy = -(i00/n/F)*(c/0.001)*grpa;
         dQdV = (i00/F)*(c/0.001)*FRT*(grpc*alphac*fv->c[four] + grpa*alphaa*fv->c[3]);
         if (wspec == 3)
            {
             d_mass_flux[store][MAX_VARIABLE_TYPES + four] -= dQdx;
             d_mass_flux[store][MAX_VARIABLE_TYPES + 3] -= dQdy;
             d_mass_flux[store][VOLTAGE] -= dQdV;
            }
         else if (wspec == 4)
            {
             d_mass_flux[store][MAX_VARIABLE_TYPES + four] += dQdx;
             d_mass_flux[store][MAX_VARIABLE_TYPES + 3] += dQdy;
             d_mass_flux[store][VOLTAGE] += dQdV;
            }
         else
            {
             d_mass_flux[store][MAX_VARIABLE_TYPES + four] += 2.*dQdx;
             d_mass_flux[store][MAX_VARIABLE_TYPES + 3] += 2.*dQdy;
             d_mass_flux[store][VOLTAGE] += 2.*dQdV;
            }
        }
     }

/***  Reaction  Ni+2 + 2e- = Ni  ***/

  if (wspec == 0 || flag == 1)
     {
      alphac = 0.21;
      n = 2;
      i00 = 1.07e-05;
      U00 = -0.2363;
      diff_phi = n*FRT*(PHI_E - PHI_S - U00);
      grpc = exp(-alphac*diff_phi);
      alphaa = 1 - alphac;
      grpa = exp(alphaa*diff_phi);
      conc = c*fv->c[0];
      Q = (i00/n/F)*(grpc*conc/0.001 - grpa);
      if (flag == 1)
         {
          mass_flux[store] += 2.*Q;
         }
      else
         {
          mass_flux[store] += Q;
         }
      if (af->Assemble_Jacobian)
        {
         dQdx = (i00/n/F)*(c/0.001)*grpc;
         dQdV = (i00/F)*FRT*(grpc*alphac*conc/0.001 + grpa*alphaa);
         if (flag == 1)
            {
             d_mass_flux[store][MAX_VARIABLE_TYPES + 0] += 2.*dQdx;
             d_mass_flux[store][VOLTAGE] += 2.*dQdV;
            }
         else
            {
             d_mass_flux[store][MAX_VARIABLE_TYPES + 0] += dQdx;
             d_mass_flux[store][VOLTAGE] += dQdV;
            }
        }
     }

  return;
} /* END of routine mass_flux_surf_NI                                       */
/****************************************************************************/
/****************************************************************************/

void 
mass_flux_surf_const(dbl mass_flux[MAX_CONC],
		     dbl d_mass_flux[MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC],
		     double Tsurf_quad, /* Surface Temperature              */
		     double Ysurf_quad[MAX_CONC], /* concentration at surface 
						 * quadrature point         */
		     int wspec,              /* subvar species number       */
		     double const_mass_flux) /* specified mass flux         */

/****************************************************************************
*
*  Function which returns the specified constant total mass flux and its
*  sensitivities (all zero).  This routine cloned from mass_flux_surf_mtc
*  by RSL on 6/8/00.
*
*****************************************************************************/
{
  int w;
  if (af->Assemble_Jacobian) {
    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      d_mass_flux[wspec][MAX_VARIABLE_TYPES + w] = 0.;
    }
    d_mass_flux[wspec][TEMPERATURE] = 0.;
  }
  mass_flux[wspec] = const_mass_flux;
  return;
} /* END of routine mass_flux_surf_const                           */ 
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void 
raoults_law(double func[],
	    double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
	    int wspec,
	    int eb_mat_liq,
	    int eb_mat_gas,      /* elem block id's of liq and gas */
	    double amb_pres,	/* ambient pressure */
	    double M1,
	    double M2,
	    double M3,
	    double M4)		/* molecular weights of four species. M1 must
				 * be the first volatile species, M2 is the
                                 * second volatile species, M3 is the condensed
                                 * phase in the liquid, and M4 is the insoluble
                                 * gas phase, e.g. air                       */


/*******************************************************************************
*  Added Riedel Vapor Pressure Model to handle saturation pressure
*  of volatile liquid components. Author: acsun (6/23/98)
*
*  Added Antoine Vapor Pressure Model to handle saturation pressure
*  of volatile liquid components. Author: acsun (5/20/98)
*
*  Function which evaluates vapor-liquid equilibrium 
* for adjacent phases using a discontinuous variable approach
*			 Author: p. r. schunk (9/14/97)
*
*  Implements the following equation:
*
*           func =   p_v[w] / p_total * X_liq_w  - X_gas_w = 0
*
*  where  X_liq_w is the liquid mole fraction of species w
*         X_gas_w is the gas mole fraction of species w
*         p_v[w] is the vapor pressure of species w in the gas above a
*                liquid of pure w
*         p_total is the total ambient pressure.
*
*  Note: It is assumed that the species dependent variable types are
*        mass fractions.
*******************************************************************************/
{
  int j_id;
  int var;
  int w;
  double phi_j;
  double C[MAX_CONC];
  double A, B_L, B_G, a, b, c, d, e;
  double psat[MAX_CONC], dpsatdt[MAX_CONC];

  if(af->Assemble_LSA_Mass_Matrix)
    return;

  /* Nonideal VP Calculations based on either ANTOINE or RIEDEL models */

   if(mp->VaporPressureModel[wspec] == ANTOINE )
     {
       antoine_psat(wspec, mp->u_vapor_pressure[wspec],
                    &psat[wspec], &dpsatdt[wspec]);
       mp->vapor_pressure[wspec] = psat[wspec];
     }

   else if(mp->VaporPressureModel[wspec] == RIEDEL )
     {
       riedel_psat(wspec, mp->u_vapor_pressure[wspec],
                   &psat[wspec], &dpsatdt[wspec]);
       mp->vapor_pressure[wspec] = psat[wspec];
     }

  /* Define some convenient/repetitive chunks to make eqns more compact */

   C[0]=0.; C[1] = 0.;
   for (w = 0; w < pd->Num_Species_Eqn; w++) C[w] = fv->c[w];
   a = M2*M3;
   b = M1*M2;
   c = M3*M1;
   d = M2*M4;
   e = M4*M1;
   A = (mp->vapor_pressure[wspec]/amb_pres);
   B_L = SQUARE(b*(1. - C[0] - C[1]) + c*C[1] + a*C[0]);
   B_G = SQUARE(b*(1. - C[0] - C[1]) + e*C[1] + d*C[0]);

  /* local contributions of boundary condition to residual and jacobian */
  
/***************************** EXECUTION BEGINS *******************************/

     if (af->Assemble_Jacobian) {
	
	/* sum the contributions to the global stiffness matrix */
	    var = MASS_FRACTION;
	    for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
	      {
		if (pd->v[pg->imtrx][var])
		  {
		    phi_j = bf[var]->phi[j_id];

		    if(pd->Num_Species_Eqn == 1)
		      {
			if (Current_EB_ptr->Elem_Blk_Id == eb_mat_liq)
			  {
			    d_func[0][MAX_VARIABLE_TYPES+0][j_id] +=  phi_j*A*M3*
			                  ((M1*(1.- C[0]) + M3*C[0]) - C[0]*(-M1+M3))/
                                          SQUARE(M1*(1.- C[0]) + M3*C[0]);
			  }
			else
			  {
			    d_func[0][MAX_VARIABLE_TYPES+0][j_id] +=  -phi_j*M4*
                                          ((M1*(1.- C[0]) + M4*C[0])- C[0]*(-M1+M4))/
                                          SQUARE(M1*(1.- C[0]) + M4*C[0]);
			  }
		      }
		    else if (pd->Num_Species_Eqn == 2)
		      {
			if (Current_EB_ptr->Elem_Blk_Id == eb_mat_liq)
			  {
			    if(wspec==0)
			      {
				d_func[0][MAX_VARIABLE_TYPES+0][j_id] +=  phi_j*A*
				   (a*b*(1.-C[1]) + a*c*C[1])/B_L;
				d_func[0][MAX_VARIABLE_TYPES+1][j_id] +=  phi_j*A*
				   (a*b*C[0] - a*c*C[0])/B_L;
			      }
			    if(wspec==1)
			      {
                                d_func[0][MAX_VARIABLE_TYPES+0][j_id] +=  phi_j*A*
				  (b*c*C[1] - a*c*C[1])/B_L;
                                d_func[0][MAX_VARIABLE_TYPES+1][j_id] +=  phi_j*A*
				  (b*c*(1-C[0]) + a*c*C[0])/B_L;
			      }
			  }
			else
			  {
			    if(wspec==0)
			      {
				d_func[0][MAX_VARIABLE_TYPES+0][j_id] +=  -phi_j*
				   (d*b*(1.-C[1]) + d*e*C[1])/B_G;
				d_func[0][MAX_VARIABLE_TYPES+1][j_id] +=  -phi_j*
				   (d*b*C[0] - d*e*C[0])/B_G;
				
			      }
			    if(wspec==1)
			      {
				d_func[0][MAX_VARIABLE_TYPES+0][j_id] +=  -phi_j*
				   (b*e*C[1] - d*e*C[1])/B_G;
				d_func[0][MAX_VARIABLE_TYPES+1][j_id] +=  -phi_j*
				   (b*e*(1.-C[0]) + d*e*C[0])/B_G;
			      }
			  }
		      }
		    else
		      {
			EH(-1,"For the Raoult's law model we only allow 1 or 2 species");
		      }
		  }
	      }
     }

 /* Since nonideal VP equation is temperature dependent, d_func is temperature
      dependent */

if(mp->VaporPressureModel[wspec] == ANTOINE 
   || mp->VaporPressureModel[wspec] == RIEDEL )
  {
    var = TEMPERATURE ;
    for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
      {
	if (pd->v[pg->imtrx][var])
	  {
	    phi_j = bf[var]->phi[j_id];

	    if (Current_EB_ptr->Elem_Blk_Id == eb_mat_liq)
	      {
	      if(pd->Num_Species_Eqn == 1)
		{
		  d_func[0][var][j_id] +=  phi_j*dpsatdt[wspec]
		         *M3*C[0]/(M1*(1.- C[0]) + M3*C[0])/amb_pres;
		}
              else if(pd->Num_Species_Eqn == 2)
		{
		  if(wspec==0)
		    {
		      d_func[0][var][j_id] +=  phi_j*dpsatdt[wspec]
			*a*C[0]/(b*(1.- C[0] - C[1]) + c*C[1] + a*C[0])/amb_pres;
		    }
		  if(wspec==1)
		    {
		      d_func[0][var][j_id] +=  phi_j*dpsatdt[wspec]
			*c*C[1]/(b*(1. -C[0] - C[1]) + c*C[1] + a*C[0])/amb_pres;
		    }
		}
	      }
	  }
      }
  }
 
if(pd->Num_Species_Eqn == 1)
  {
    if (Current_EB_ptr->Elem_Blk_Id == eb_mat_liq)
      {
	/*	*func += (mp->vapor_pressure[wspec]/amb_pres)*fv->c[wspec]; */
	*func +=     A*M3*C[0]/(M1*(1.- C[0]) + M3*C[0]);
      }
    else if (Current_EB_ptr->Elem_Blk_Id == eb_mat_gas)
      {
	/* *func += -fv->c[wspec]; */
	*func +=       -M4*C[0]/(M1*(1.- C[0]) + M4*C[0]);
      }
    else
      {
	EH(-1,"Something screwy in your material block specs on the VL_EQUIL card");
      }
  }
else if (pd->Num_Species_Eqn == 2)
  {
    if (Current_EB_ptr->Elem_Blk_Id == eb_mat_liq)
      {
	if(wspec==0)
	  {
	    *func=A*a*C[0]/(b*(1.- C[0] - C[1]) + c*C[1] + a*C[0]);
	  }
	if(wspec==1)
	  {
	    *func=A*c*C[1]/(b*(1. -C[0] - C[1]) + c*C[1] + a*C[0]);
	  }
      }
    else
      {
	if(wspec==0)
	  {
	    *func= -d*C[0]/(b*(1. - C[0] - C[1]) + e*C[1] + d*C[0]);
	  }
	if(wspec==1)
	  {
	    *func= -e*C[1]/(b*(1. - C[0] - C[1]) + e*C[1] + d*C[0]);
	  }
      }
  }
else
  {
    EH(-1,"For the Raoult's law model we only allow 1 or 2 species");
  }
  return;
} /* END of routine raoults_law   */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void 
raoults_law_new(double func[],
	        double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
		BOUNDARY_CONDITION_STRUCT *bc)

    /************************************************************************
     *
     * raoults_law_new():
     *
     *  Function which evaluates vapor-liquid equilibrium 
     * for adjacent phases using a discontinuous variable approach.
     * It gets called from both sides of the interface. 
     *
     *
     *  Implements the following equation:
     *
     *           func =   p_v[w] / p_total * X_liq_w  - X_gas_w = 0
     *
     *  where  X_liq_w is the liquid mole fraction of species w
     *         X_gas_w is the gas mole fraction of species w
     *         p_v[w] is the vapor pressure of species w in the gas above a
     *                liquid of pure w
     *         p_total is the total ambient pressure.
     *
     *  Note: It is assumed that the species dependent variable types are
     *        mass fractions.
     *		
     * Input
     * ----------
     *  wspec -> Species number on both sides of the interface whose
     *           concentrations are linked via a raoult's law-type 
     *           expression.
     *  eb_mat_liq = Element block id for the first phase -> this is 
     *               identified with the liquid phase
     *  eb_mat_gas = Element block id for the second phase -> this is 
     *               identified with the gas phase
     ***********************************************************************/
{
  int j_id, var, k;
  int wspec = bc->BC_Data_Int[0];
  int eb_mat_liq = bc->BC_Data_Int[1];
  int liquidSide = (Current_EB_ptr->Elem_Blk_Id == eb_mat_liq);
  double phi_j, psat_w, dpsatdT_w = 0.0, Mbar, A, amb_pres, Xk_wspec;
  double tmp1, tmp2, tmp3;
  double *mw = mp->molecular_weight;
  double *Yk =  fv->c;
  

  if (af->Assemble_LSA_Mass_Matrix) return;

  /*
   * Obtain the ambient thermodynamic pressure
   * -> Currently, we do not consider any compressibility effects, nor
   *    do we link the thermodynamic pressure with the hydrodynamic
   *    pressure. However, this could be done in the future.
   */
  amb_pres = upd->Pressure_Datum;

  /*
   * Calculate the average molecular weight of the mixture 
   */
  Mbar = wt_from_Yk(mp->Num_Species, Yk, mw);

  /*
   * Calculate the mole fraction of the pertinent species
   */
  Xk_wspec = Yk[wspec] * Mbar / mw[wspec];

  /*
   * Residual Contributions
   */
  if (liquidSide) {
    if (mp->VaporPressureModel[wspec] == ANTOINE) {
      antoine_psat(wspec, mp->u_vapor_pressure[wspec],
		   &psat_w, &dpsatdT_w);
      mp->vapor_pressure[wspec] = psat_w;
    } else if (mp->VaporPressureModel[wspec] == RIEDEL) {
      riedel_psat(wspec, mp->u_vapor_pressure[wspec],
		  &psat_w, &dpsatdT_w);
      mp->vapor_pressure[wspec] = psat_w;
    } else {
      psat_w = mp->vapor_pressure[wspec];
    }
    A = psat_w/amb_pres;
  } else {
    A = -1;
  }
  func[0] = A * Xk_wspec;

  if (af->Assemble_Jacobian) {
    var = MASS_FRACTION;
    if (pd->v[pg->imtrx][var]) {
      tmp2 = A * Mbar * Xk_wspec;
      if (mp->Num_Species == mp->Num_Species_Eqn) {
	for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
	  phi_j = bf[var]->phi[j_id];
	  tmp3 = phi_j * tmp2;
	  for (k = 0; k < mp->Num_Species; k++) {
	    d_func[0][MAX_VARIABLE_TYPES + k][j_id] = - tmp3 / mw[k];
	  }
	  d_func[0][MAX_VARIABLE_TYPES + wspec][j_id] +=
	      phi_j * A * Mbar / mw[wspec];
	}
      } else {
        tmp1 = 1.0/ mw[mp->Num_Species_Eqn];
	for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
	  phi_j = bf[var]->phi[j_id];
          tmp3 = phi_j * tmp2;
	  for (k = 0; k < mp->Num_Species_Eqn; k++) {
	    d_func[0][MAX_VARIABLE_TYPES + k][j_id] = 
                tmp3 * (tmp1 - 1.0 / mw[k]);
	  }
	  d_func[0][MAX_VARIABLE_TYPES + wspec][j_id] +=
	      phi_j * A * Mbar / mw[wspec];
	}

      }
    }
  
    if (dpsatdT_w != 0.0) {
      var = TEMPERATURE;
      if (pd->v[pg->imtrx][var]) {
	tmp1 =  dpsatdT_w * Xk_wspec / amb_pres;
	for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
	  phi_j = bf[var]->phi[j_id];
	  d_func[0][var][j_id] = phi_j * tmp1;
	}
      }
    }
  }
  return;
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void 
flory_huggins(double func[],
	      double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
	      int wspec,
	      int eb_mat_gas,
	      int eb_mat_liq,    /* element block id's for gas and liquid */
	      int mode,		/* VOLUME or MASS based formulation */
	      double amb_pres)         /* ambient pressure */

/*******************************************************************************
*
*  Function which evaluates vapor_polymer equilibrium based on 
*  Flory_Huggins solution theory.  mode=VOLUME is a volume fraction based formulation,
*  which is common for polymer-solvent solutions. i.e. C[] = volume fraction
*  mode=MASS is a mass fraction based formulation.
*  The definition for VLE in these solutions were based on Prausnitz (1986):
*
*  y(i) = Ptot*activity(i)/psat
*  acti_coef(i) = activity(i)/mass_frac(i) for polymer-solvent systems.
*  where y(i) mole_fraction of solvent i in the vapor.
*        Ptot = total pressure in the vapor.
*        activity(i) = departure from ideality for the solvent i in the liquid.
*        psat(i) = saturation pressure for pure liquid i. Note this is the
*                  REFERENCE STATE.
*
*			 Author: a.c.sun (8/98)
*******************************************************************************/
{
  
/* Local variables */
  
  int i, j, jac, j_id, k, l, mn, var;
  int Num_S1, Num_S2;
  double flory1[MAX_CONC],flory2[MAX_CONC],flory3[MAX_CONC],flory[MAX_CONC];
  double df1_dc[MAX_CONC][MAX_CONC],df2_dc[MAX_CONC][MAX_CONC];
  double df3_dc[MAX_CONC][MAX_CONC],df_dc[MAX_CONC][MAX_CONC];
  double truedf_dc[MAX_CONC][MAX_CONC];
  double phi_j, A;
  double C[MAX_CONC], vol[MAX_CONC], sv[MAX_CONC];
  double y_mol[MAX_CONC], y_mass[MAX_CONC];
  double mw[MAX_CONC], prod[MAX_CONC]; 
  double activity[MAX_CONC];
  double dact_dC[MAX_CONC][MAX_CONC],dv_dw[MAX_CONC][MAX_CONC];
  double dmol_dC[MAX_CONC][MAX_CONC];
  double psat[MAX_CONC], dpsatdt[MAX_CONC];
  double bottom, prod2, sum_C;
  double chi[MAX_CONC][MAX_CONC]; /* chi is the binary interaction parameter*/
  double mw_last=0, dmdv=0; /* Molecular weight of non-condensable and conversion factor */

  memset(y_mass, 0, sizeof(double)*MAX_CONC);
  memset(prod, 0, sizeof(double)*MAX_CONC);

  if(af->Assemble_LSA_Mass_Matrix)
    return;

  /* Nonideal VP Calculations based on either ANTOINE or RIEDEL models */

   if(mp->VaporPressureModel[wspec] == ANTOINE )
     {
       antoine_psat(wspec, mp->u_vapor_pressure[wspec],
                    &psat[wspec], &dpsatdt[wspec]);
       mp-> vapor_pressure[wspec] = psat[wspec];
     }

   else if(mp->VaporPressureModel[wspec] == RIEDEL )
     {
       riedel_psat(wspec, mp->u_vapor_pressure[wspec],
                   &psat[wspec], &dpsatdt[wspec]);
       mp-> vapor_pressure[wspec] = psat[wspec];
     }

  /* Define some convenient/repetitive chunks to make eqns more compact */

   Num_S1 = pd->Num_Species_Eqn + 1;
   Num_S2 = pd->Num_Species_Eqn - 1;

   memset(flory1, 0,sizeof(double)*MAX_CONC);
   memset(flory2, 0,sizeof(double)*MAX_CONC);
   memset(flory3, 0,sizeof(double)*MAX_CONC);
   memset(flory, 0,sizeof(double)*MAX_CONC);
   memset(truedf_dc, 0,sizeof(double)*MAX_CONC*MAX_CONC);
   memset(df_dc, 0,sizeof(double)*MAX_CONC*MAX_CONC);
   memset(df1_dc, 0,sizeof(double)*MAX_CONC*MAX_CONC);
   memset(df2_dc, 0,sizeof(double)*MAX_CONC*MAX_CONC);
   memset(df3_dc, 0,sizeof(double)*MAX_CONC*MAX_CONC);
 
   A = mp->vapor_pressure[wspec]/amb_pres;

   if (Current_EB_ptr->Elem_Blk_Id == eb_mat_liq)
     {
      /* calculate liquid activities of solvent components in a polymeric
         solution and their derivatives  
         use generalized activity expression for multicomponents  */

       /* mass based formulation must require specific volume for conversion
         to volume fractions */

      if(mode==MASS)
	 {
	   bottom = 0.;

	   for(i=0;i<Num_S1;i++)
	     {
	       if(mp->specific_volume[i] < 0.)
		 {
		   EH(-1, "Specific volume not specified in the material file.");
		 }
	       else
		 {
		   sv[i] = mp->specific_volume[i];
		 }
	     }
	   for(i=0;i<pd->Num_Species_Eqn;i++)
	     {
	       y_mass[i] = fv->c[i];
	       bottom += y_mass[i]*(sv[i]-sv[pd->Num_Species_Eqn]);
	     }
           bottom += sv[pd->Num_Species_Eqn];
	   
	   for(i=0;i<pd->Num_Species_Eqn; i++)
	     {
	       C[i] = y_mass[i]*sv[i]/bottom;

	       for(j=0;j<pd->Num_Species_Eqn; j++)
		 {
		   if (j==i)
		     {
		       dv_dw[i][j] = sv[j]/bottom 
			 -y_mass[i]*sv[i]*(sv[j]-sv[pd->Num_Species_Eqn])/(bottom*bottom);
		     }
		   else
		     {
		       dv_dw[i][j] = -y_mass[i]*sv[i]
			 *(sv[j]-sv[pd->Num_Species_Eqn])/(bottom*bottom);
		     }
		 }
	     }
	 }
      else
	{
	  for (i = 0; i<pd->Num_Species_Eqn; i++) 
	    {
	      C[i] = fv->c[i] ;
	      for(j=0;j<pd->Num_Species_Eqn; j++)
		{
		  dv_dw[i][j]= delta(i,j);
		}
	    }
	}

       sum_C = 0.;
       for (i = 0; i<pd->Num_Species_Eqn; i++) 
	 {
	   sum_C += C[i];
	 }
       for (i = 0; i<Num_S1; i++) 
	 {
	   if(mp->molar_volume[i] < 0.)
	     {
	       EH(-1, "Molar volume not specified in the material file.");
	     }
	   else
	     {
	       vol[i] = mp->molar_volume[i];
	     }
	   for(k=0;k<Num_S1;k++)
	     {
	     if(mp->flory_param[i][k] < 0.)
	       {
		 EH(-1, "Flory-Huggins binary parameters not specified in the material file.");
	       }
	     else
	       {
		 chi[i][k] = mp->flory_param[i][k];
	       }
	     }
	 }
       
       for(k=0;k<pd->Num_Species_Eqn;k++)
	 {
	   flory1[wspec] += vol[wspec]*(delta(k,wspec)-C[k])/vol[k]
	             +vol[wspec]*C[k]/vol[pd->Num_Species_Eqn];

	   df1_dc[wspec][k]= (vol[wspec]/vol[pd->Num_Species_Eqn])
	              - (vol[wspec]/vol[k]) ;
	 }
       flory1[wspec] += -vol[wspec]/vol[pd->Num_Species_Eqn];

       if (pd->Num_Species_Eqn > 1)
	 {
	   for(k=1;k<pd->Num_Species_Eqn; k++)
	     {
	       for(jac=0;jac<k;jac++)
		 {   
		   flory2[wspec] += (delta(wspec,jac)*C[k]+(vol[wspec]/vol[jac])
			       *C[jac]*(delta(k,wspec)-C[k]))*chi[jac][k];
		 }
	       /* derivative for the 1st component */
	       df2_dc[wspec][0] += (vol[wspec]*(delta(k,wspec)-C[k])/vol[0])
		 *chi[0][k];
	     }

	   /* derivative for the last component */
	   for(jac=0;jac<Num_S2;jac++)
	     {
	       df2_dc[wspec][Num_S2] += (delta(wspec,jac)-vol[wspec]*C[jac]
				  /vol[jac])*chi[jac][Num_S2];
	     }
	   /* derivative for the components in between */
	   for(l=1;l<Num_S2 ;l++)
	     {
	       for(jac=0;jac<l;jac++)
		 {
		  df2_dc[wspec][l] += (delta(wspec,jac)-(vol[wspec]*C[jac]
				  /vol[jac]))*chi[jac][l];
		 }
	       for(k=l+1;k<pd->Num_Species_Eqn; k++)
		 {
		  df2_dc[wspec][l] +=(vol[wspec]*(delta(k,wspec)-C[k])/vol[l])
		                 *chi[l][k];
		 }
	     }
	 }
       else
	 {
	   flory2[wspec] = 0.;
           for(k=0;k<pd->Num_Species_Eqn; k++)  df2_dc[wspec][k] = 0.;
	 }

       for(jac=0;jac<pd->Num_Species_Eqn; jac++)
	 {
	   flory3[wspec] +=(delta(wspec,jac)-(vol[wspec]/vol[jac])*C[jac])
	           *(1.-sum_C)*chi[jac][pd->Num_Species_Eqn];
           df3_dc[wspec][jac] -= (vol[wspec]/vol[jac])*(1.-sum_C)
                            *chi[jac][pd->Num_Species_Eqn];

	   for(mn=0;mn<pd->Num_Species_Eqn;mn++)
	     {
	       df3_dc[wspec][jac]
		 -= (delta(wspec,mn)-(vol[wspec]/vol[mn])*C[mn])
		       *chi[mn][pd->Num_Species_Eqn];
	     }
	 }
       flory[wspec]= flory3[wspec];

       /* check the simplest case: 1solvent, 1polymer
       check=(1-sum_C)+chi[0][1]*(1.-sum_C)*(1.-sum_C);
       printf("flory = %e, check =%e\n", flory[i],check); */

       for(k=0;k<pd->Num_Species_Eqn;k++)
         {
	   df_dc[wspec][k]= df3_dc[wspec][k]; 
	 }
      

       for(i=0;i<pd->Num_Species_Eqn;i++)
	 {
	   for(j=0;j<pd->Num_Species_Eqn;j++)
	     {
	       truedf_dc[wspec][i] += df_dc[wspec][j]*dv_dw[j][i];
	     }
	 }

      activity[wspec] = C[wspec]*exp(flory[wspec]);

      *func +=  activity[wspec]*A;

      for (i=0; i<pd->Num_Species_Eqn; i++)
       {
       	 if (mode==MASS)
	   {
	     dact_dC[wspec][i] = (dv_dw[wspec][i]*exp(flory[wspec])
			      +activity[wspec]*truedf_dc[wspec][i]);
	   }
	 else
	   {
	     dact_dC[wspec][i] = (dv_dw[wspec][i]*exp(flory[wspec])
			      +activity[wspec]*df_dc[wspec][i]);
	   }
       }
     
      if (af->Assemble_Jacobian) {
	
	/* sum the contributions to the global stiffness matrix */
	var = MASS_FRACTION;
	for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
	  {
	    if (pd->v[pg->imtrx][var])
	      {
		phi_j = bf[var]->phi[j_id];
		for (i=0; i<pd->Num_Species_Eqn; i++)
		  {
		    d_func[0][MAX_VARIABLE_TYPES+i][j_id] +=  
		      phi_j*A*dact_dC[wspec][i]; 
		  }
	      }
	  }
      }
     }

   else if (Current_EB_ptr->Elem_Blk_Id == eb_mat_gas)
     {
       /* from the activity, vapor phase composition is determined,
      we write the vapor mole fraction in terms of mass fraction */


       if(mode==VOLUME)
	 {
	   dmdv = mp->density;
	 }
       else if(mode==MASS)
	 {
	   dmdv = 1.;
	 }

       for (i = 0; i<pd->Num_Species_Eqn; i++) C[i] = fv->c[i]*dmdv;

       bottom=0.;
       prod2=1.;

       for (i = 0; i<pd->Num_Species_Eqn; i++) 
	 {
	   if(mp->molecular_weight[i] < 0.)
	   {
	     EH(-1, "Molecular weight of the species not specified in the material file");
	   }
	 else
	   {
	     mw[i]  = mp->molecular_weight[i];
	   }
	 }

       if(mp->molecular_weight[pd->Num_Species_Eqn] < 0.)
	 {
	 EH(-1, "M.W. of a non-condensable species not specified in the material file");
	 }
       else
	 {
	   mw_last = mp -> molecular_weight[pd->Num_Species_Eqn];
	 }

       for(i=0;i<pd->Num_Species_Eqn;i++)
	 {
	   y_mass[i] = C[i];
	   prod[i]=1.;
	   for(j=0;j<pd->Num_Species_Eqn;j++)
	     {
	       prod[i] *= (i==j? 1. : mw[j]);
	     } 
	   prod2 *= mw[i];
	 }

       for(i=0;i<pd->Num_Species_Eqn;i++)
	 {
	   bottom += (mw_last*prod[i]-prod2)*y_mass[i];
	 }
       bottom += prod2;
   
       y_mol[wspec] = mw_last*prod[wspec]*y_mass[wspec]/bottom;
   
       *func -= y_mol[wspec];

       /* write out the derivatives of mol_frac w.r.t. mass_frac 
	  then from mass_frac to vol_frac in the VAPOR */
       
       for(j=0;j<pd->Num_Species_Eqn;j++)
	 {
	   if (j==wspec)
	     {
	       dmol_dC[wspec][j] = dmdv*(mw_last*prod[wspec]/bottom
                -mw_last*prod[wspec]*y_mass[wspec]*(mw_last*prod[j]-prod2)
                /(bottom*bottom));

	     }
	   else
	     {
	       dmol_dC[wspec][j] = dmdv*(-mw_last*prod[wspec]
		 *y_mass[wspec]*(mw_last*prod[j]-prod2)/(bottom*bottom)) ;
	     }
	 }

       if (af->Assemble_Jacobian) {
	 
	 /* sum the contributions to the global stiffness matrix */
	 var = MASS_FRACTION;
	 for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
	   {
	     if (pd->v[pg->imtrx][var])
	       {
		 phi_j = bf[var]->phi[j_id];
		 
		 for (i=0; i<pd->Num_Species_Eqn; i++)
		   {
		     d_func[0][MAX_VARIABLE_TYPES+i][j_id] -=  
		       phi_j*dmol_dC[wspec][i]; 
		   }
	       }
	   }
       }
     }

  else
    {
     EH(-1,"Material file id for gas and liquid in VL_POLY card must be specified.");
     return;
    }



/* Since nonideal VP equation is temperature dependent, d_func is temperature
      dependent */

     if(mp->VaporPressureModel[wspec] == ANTOINE 
	|| mp->VaporPressureModel[wspec] == RIEDEL )
       {
	 var = TEMPERATURE ;
	 for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
	   {
	     if (pd->v[pg->imtrx][var])
	       {
		 phi_j = bf[var]->phi[j_id];
		 
		 d_func[0][var][j_id] +=  phi_j*dpsatdt[wspec]
		   *activity[wspec]/amb_pres;
	       }
	   }
       }
} /* END of routine flory_huggins   */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void 
kinematic_species_bc(double func[DIM],
		     double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
		     int wspec,	/* current species no. */
		     double vnormal, /* normal velocity */
		     double x_dot[MAX_PDIM], /* mesh velocity vector         */
		     dbl tt,	/* parameter to vary time integration from 
				   explicit (tt = 1) to implicit (tt = 0)    */
		     dbl dt)	/* current value of the time step            */

    /*************************************************************************
     *
     * kinematic_species_bc():
     *
     *  Function which evaluates the species kinematic boundary condition.
     *  This function is applied to replace the negative of the diffusive
     *  flux out of the liquid side of the interface:
     *     (actually, in this case, it replaces the difference between
     *      this flux and the diffusive flux in the "gas" for species, i)
     *
     *   - n_l dot j_k_l + n_l dot j_k_g) 
     *
     *     = n_g dot (rho_g * Y_k_g * (v_g  - vs)) - vnormal +
     *       n_l dot (rho_l * Y_k_l * (v_l  - vs)) - vnormal
     *
     * where
     *       n_g = outward facing normal from the gas side
     *       n_l = outward facing normal from the liquid side (n_g = - n_l)
     *       vs = xdot = velocity of the moving mesh
     *       rho_g = density in the gas
     *       rho_l = density in the liquid
     *       Y_k_g = mass fraction of species k in the gas
     *       Y_k_l = mass fraction of species k in the liquid
     *       vnormal = velocity in a moving reference frame
     *                 (first single float on the bc data card,
     *                  almost always equal to zero)
     *       v_g = mass average velocity in the gas
     *       v_l = mass average velocity in the liquid.
     *
     *
     *            Author: P. R. Schunk    (9/24/97)
     *
     *  It replaces the negative of the outward diffusive flux because
     *  the time derivative in the species residual equation has a
     *  negative sign.
     *
     *  This function is called from both the gas and liquid sides of the
     *  interface, forming each of the terms in the expression
     *  above one at a time.
     *************************************************************************/
{
  int j, kdir, var, p;
  double phi_j;
  double mass_c=fv->c[wspec], mass_dc=1.0;

  if(mp->Species_Var_Type == SPECIES_DENSITY)
	{
         mass_c = fv->c[wspec];
         mass_dc = 1.0;
	}
  else if(mp->Species_Var_Type == SPECIES_CONCENTRATION)
	{
         mass_c = fv->c[wspec];
         mass_dc = 1.0;
	}
  else if(mp->Species_Var_Type == SPECIES_MASS_FRACTION
           || mp->Species_Var_Type == SPECIES_UNDEFINED_FORM)
	{
         mass_c = mp->density*fv->c[wspec];
         mass_dc = mp->density;
	}
/*  else if(mp->Species_Var_Type == SPECIES_MOLE_FRACTION)
	{
         mass_c = mp->molecular_weight[wspec]*fv->c[wspec];
         mass_dc = mp->density;
	}*/
  else
	{ EH(-1,"That species formulation not done in kinematic_species\n"); }

  if (af->Assemble_LSA_Mass_Matrix) {
    for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
	for (p = 0; p < pd->Num_Dim; p++) {
	  var = MESH_DISPLACEMENT1 + p;
	  if (pd->v[pg->imtrx][var]) {
	    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
	      phi_j = bf[var]->phi[j];
	      d_func[0][var][j] -= (mass_c *
				    phi_j * fv->snormal[kdir]);
	    }
	  }
	}
    }
    return;
  }

  if (af->Assemble_Jacobian) {
    for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
      for (p = 0; p < pd->Num_Dim; p++) {
	var = MESH_DISPLACEMENT1 + p;
	if (pd->v[pg->imtrx][var]) {
	  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
	    phi_j = bf[var]->phi[j];
	    d_func[0][var][j] += (mass_c *
				  (fv->v[kdir] - x_dot[kdir]) *
				  fv->dsnormal_dx[kdir][p][j]);
	    if (TimeIntegration != 0 && p == kdir) {
              d_func[0][var][j] += (mass_c *
				    ( -(1.+2.*tt)* phi_j/dt) *
				    fv->snormal[kdir] * delta(p, kdir));
	    }
	  }
	}
      }
	      
      var = VELOCITY1 + kdir;
      if (pd->v[pg->imtrx][var]) {
	for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
	  phi_j = bf[var]->phi[j];
          d_func[0][var][j] += (mass_c *
	                        phi_j * fv->snormal[kdir]);
	}
      }

      var = MASS_FRACTION;
      if (pd->v[pg->imtrx][var]) {
	for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)	{
	  phi_j = bf[var]->phi[j];
	  d_func[0][MAX_VARIABLE_TYPES+wspec][j] += (mass_dc * phi_j *
		     (fv->v[kdir] - x_dot[kdir]) * fv->snormal[kdir]);
	}
      }
    }
  } /* end of if Assemble_Jacobian */
  
  /* Calculate the residual contribution	*/
  *func = -vnormal;
  for (kdir = 0; kdir < pd->Num_Dim; kdir++) {
    *func += mass_c * (fv->v[kdir] - x_dot[kdir]) * fv->snormal[kdir];
  }
} /* END of routine kinematic_species_bc  */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void 
mass_flux_surf_bc(double func[],
		  double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
		  int wspec,	/* species number of this boundary condition */
		  double mass_tran_coeff, /* Mass transfer coefficient       */
		  double Y_c,	/* bath concentration 	                     */
		  dbl dt,	/* current value of the time step            */
		  dbl tt)	/* parameter to vary time integration        */

     /*******************************************************************************
      *
      *  Function which calculates the surface integral for convective
      *  mass transfer.
      *
      *  ----------------------------------------------------------------------------
      *
      *  Functions called:
      *  mass_flux_surf_mtc  -- calculates mass flux for one component
      *
      *  ----------------------------------------------------------------------------
      *
      *******************************************************************************/
{  
  int j, j_id, w1, dim, kdir, var, jvar;
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

  mass_flux_surf_mtc(mp->mass_flux, mp->d_mass_flux, fv->T,
		     fv->c, wspec, mass_tran_coeff, Y_c);

  dim   = pd->Num_Dim;
  Y_w = fv->c[wspec];

  /* get the convection velocity (it's different for arbitrary and
     lagrangian meshes) */

  /* get deformation gradient, etc. if lagrangian mesh with inertia */
  if ((pd->MeshMotion == LAGRANGIAN ||
       pd->MeshMotion == DYNAMIC_LAGRANGIAN) && pd->MeshInertia == 1) {
      err = belly_flop(elc->lame_mu);
      EH(err, "error in belly flop");
      if( neg_elem_volume ) return;
    }

/*  if (mp->PorousMediaType == CONTINUOUS) */

    err = get_convection_velocity(vconv, vconv_old, d_vconv, dt, tt);
  EH(err, "Error in calculating effective convection velocity");

  /*
   *  Calculate the Residual contribution for this boundary condition
   */
  if (af->Assemble_Residual) {
    *func -= mp->mass_flux[wspec];
    for (kdir = 0; kdir < dim; kdir++) {
      *func += Y_w * vconv[kdir] * fv->snormal[kdir];
    }
  }
  
  if (af->Assemble_Jacobian ) 
    {
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
	
	for (kdir=0; kdir<dim; kdir++) 
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
       * J_s_V
       */
      var=VOLTAGE;
      if (pd->v[pg->imtrx][var]){
	for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
	  j_id = j;
	  phi_j = bf[var]->phi[j_id];
	  
	  d_func[0][var][j_id] -= mp->d_mass_flux[wspec][VOLTAGE] 
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

/****************************************************************************/

void 
yflux_disc_rxn_bc(double func[],
		  double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
		  int wspec_a,	/* species number of this boundary condition */
		  int mat_a,	/* material block for side 1, reactant */
		  int mat_b,	/* material block for side 2, product  */
		  double kf,    /* forward rate constant      */
		  double kr,    /* reverse rate constant      */
		  dbl dt,	/* current value of the time step            */
		  dbl tt)	/* parameter to vary time integration        */

     /*******************************************************************************
      *
      *  Function which calculates the surface integral for convective
      *  mass transfer and adds in the interfacial reaction: kf*cf-kr*cr.
      *  This routine must be used with discontinuous concentration.
      *
      *  ----------------------------------------------------------------------------
      *
      *  Functions called:
      *  mass_flux_surf_mtc  -- calculates mass flux for one component
      *
      *  ----------------------------------------------------------------------------
      *
      *******************************************************************************/
{  
  int j, j_id, w1, dim, kdir, var, jvar, wspec_b, wspec;
  double phi_j;
  double Y_w; /* local concentration of current species */

  double vconv[MAX_PDIM]; /*Calculated convection velocity */
  double vconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;

  /* MATRL_PROP_STRUCT *mp_a = 0, *mp_b = 0; */

  int err=0;

  wspec= wspec_b = wspec_a; /* will this work? */

  /***************************** EXECUTION BEGINS ****************************/

  /* call routine to calculate surface flux of this component and it's
   * sensitivity to all variable types
   */

  /* need to pass down this struct ... ? */

  /* ebIndex_a = ebID_to_ebIndex(mat_a); */
  /* ebIndex_b = ebID_to_ebIndex(mat_b); */
  /* mp_a = mp_glob[Matilda[ebIndex_a]]; */
  /* mp_b = mp_glob[Matilda[ebIndex_b]]; */
  
  /* how do we get these concentrations */

  /* if (Current_EB_ptr->Elem_Blk_Id == mat_a) */
  /*     { */
  /* 	src_rxn = kf*Y_w-kr*Y_otherside; */
  /*     } */
  /* else */
  /*     { */
  /* 	src_rxn = kf*Y_otherside -kr*Y_w; */
  /*     } */



  /* interface_rxn(mp->mass_flux, mp->d_mass_flux, fv->T,fv->c, wspec, kf, kr); */

  dim   = pd->Num_Dim;
  Y_w = fv->c[wspec_a];

  /* get the convection velocity (it's different for arbitrary and
     lagrangian meshes) */

  /* get deformation gradient, etc. if lagrangian mesh with inertia */
  if ((pd->MeshMotion == LAGRANGIAN ||
       pd->MeshMotion == DYNAMIC_LAGRANGIAN) && pd->MeshInertia == 1) {
      err = belly_flop(elc->lame_mu);
      EH(err, "error in belly flop");
      if( neg_elem_volume ) return;
    }

/*  if (mp->PorousMediaType == CONTINUOUS) */

    err = get_convection_velocity(vconv, vconv_old, d_vconv, dt, tt);
  EH(err, "Error in calculating effective convection velocity");

  /*
   *  Calculate the Residual contribution for this boundary condition
   */
  if (af->Assemble_Residual) {
    *func -= mp->mass_flux[wspec];
    for (kdir = 0; kdir < dim; kdir++) {
      *func += Y_w * vconv[kdir] * fv->snormal[kdir];
    }
  }
  
  if (af->Assemble_Jacobian ) 
    {
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
	
	for (kdir=0; kdir<dim; kdir++) 
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
       * J_s_V
       */
      var=VOLTAGE;
      if (pd->v[pg->imtrx][var]){
	for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
	  j_id = j;
	  phi_j = bf[var]->phi[j_id];
	  
	  d_func[0][var][j_id] -= mp->d_mass_flux[wspec][VOLTAGE] 
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
  
} /* END of routine yflux_disc         */

/****************************************************************************/

void
mass_flux_surf_user_bc(double func[DIM],
                       double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                       const int wspec,
                       const double p[],  /* parameters to parameterize heat transfer model*/
                       const dbl time)

/******************************************************************************
*
*  Function which calculates the surface integral for user-defined
*  mass flux.
*
*  Functions called:
*  mass_flux_user_surf (user_bc.c)  -- calculates mass flux for one component
*
*  ----------------------------------------------------------------------------
*
******************************************************************************/
{
  
  int j, j_id, w1;
  int var;
  double phi_j;
  
  /* Execution begins */

  /* Grab your user-defined mass flux from routine in user_bc.c */
  mass_flux_user_surf (mp->mass_flux, mp->d_mass_flux, wspec, p, time);

  /* Load up func and d_func and return to apply_integrated BC */  
  
  if (af->Assemble_Residual ) 
    {
      *func -= mp->mass_flux[wspec];
    }
  
  if (af->Assemble_Jacobian ) 
    {
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
	}
      }

      /*
       * J_s_V
       */
      var=VOLTAGE;
      if (pd->v[pg->imtrx][var]){
	for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
	  j_id = j;
	  phi_j = bf[var]->phi[j_id];
	  
	  d_func[0][var][j_id] -= mp->d_mass_flux[wspec][VOLTAGE] 
	    * phi_j; 
	}
      }
 
    } /* End of if Assemble_Jacobian */
  return;
} /* END of routine mass_flux_surf_user_bc                                   */

/*****************************************************************************/
/****************************************************************************/

void
mass_flux_surf_etch(double func[DIM],
                    double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                    const int wspec,
                    const int etch_plane,
                    const double time,
                    const double dt,
                    const double tt )

/******************************************************************************
*
*  Function which calculates the surface integral for mass flux resulting from
*  etching reaction
*
*  Right now, it only handles KOH wet etching on crystalline silicon surface plane 100
*  It also assumes a species ordering as follows:
*                   0: H2O - water
*                   1: KOH - potassium hydroxide
*                   2: H2 - hydrogen
*                   3: Silicon hydroxyl byproducts
*
*  Flux units are given in CGS
*
*  Functions called:
*  calc_KOH_Si_etch_rate_100  -- calculates etch rate based on local concentration
*
*  Kristianto Tjiptowidjojo (5/2017)
*  ----------------------------------------------------------------------------
*
******************************************************************************/
{

  int j_id, w1;
  int kdir, var;
  int err = 0;
  int dim = pd->Num_Dim;
  double phi_j;

  double mass_flux[MAX_CONC] = {0.0};
  double d_mass_flux[MAX_CONC][MAX_VARIABLE_TYPES + MAX_CONC] = {{0.0}};

  struct Species_Conservation_Terms s_terms;
  /* Use fake values for this since we do not need the SUPG term */
  double h[DIM] = {1.0, 1.0, 1.0};

  double etch_rate = 0.0;
  double d_etch_rate_d_C[MAX_CONC] = {0.0};

  /* Bulk density of crystalline silicon (g/cm^3) */
  double rho_bulk_Si = 2.3290;

  /* Molecular weight in mole/g */
  double MW_H2O = 18.01528;
  double MW_OH = 17.008;
  double MW_Si = 28.0855;
  double MW_H2 = (2.0 * 1.00794);
  double MW_SiO2OH2 = (28.0855 + 2.0*15.9994 + 2.0*17.008);

  /*
   *  Initialize the Species_Conservation_Terms temporary structure
   *  before filling it up
   */
  zero_structure(&s_terms, sizeof(struct Species_Conservation_Terms), 1);
  err = get_continuous_species_terms(&s_terms, time, tt, dt, h);
  EH(err, "get_continuous_species_terms");

 /* Right now it only handles KOH wet etching on plane 100 of crystalline silicon*/
  if (etch_plane == 100)
    {
     /* Get etch rate */
     etch_rate = calc_KOH_Si_etch_rate_100(d_etch_rate_d_C);

     /* Export it to mass_flux array, depending on their stochiometric coefficient */
     switch (wspec)
       {
        case 0: /* Water */
           mass_flux[wspec] = 2.0 * rho_bulk_Si/MW_Si * MW_H2O * etch_rate;
           break;

        case 1: /* OH */
           mass_flux[wspec] = 2.0 * rho_bulk_Si/MW_Si * MW_OH * etch_rate;
           break;

        case 2: /* H2 */
           mass_flux[wspec] = -2.0 * rho_bulk_Si/MW_Si * MW_H2 * etch_rate;
           break;

        case 3: /* SiO2OH2 */
           mass_flux[wspec] = -1.0 * rho_bulk_Si/MW_Si * MW_SiO2OH2 * etch_rate;
           break;
       }

     /* Export sensitivity of mass flux array */
     switch (wspec)
       {
        case 0: /* Water */
           d_mass_flux[wspec][MAX_VARIABLE_TYPES+0] = 2.0 * rho_bulk_Si/MW_Si * MW_H2O * d_etch_rate_d_C[0];
           d_mass_flux[wspec][MAX_VARIABLE_TYPES+1] = 2.0 * rho_bulk_Si/MW_Si * MW_H2O * d_etch_rate_d_C[1];
           break;

        case 1: /* KOH */
           d_mass_flux[wspec][MAX_VARIABLE_TYPES+0] = 2.0 * rho_bulk_Si/MW_Si * MW_OH * d_etch_rate_d_C[0];
           d_mass_flux[wspec][MAX_VARIABLE_TYPES+1] = 2.0 * rho_bulk_Si/MW_Si * MW_OH * d_etch_rate_d_C[1];
           break;

        case 2: /* H2 */
           d_mass_flux[wspec][MAX_VARIABLE_TYPES+0] = -2.0 * rho_bulk_Si/MW_Si * MW_H2 * d_etch_rate_d_C[0];
           d_mass_flux[wspec][MAX_VARIABLE_TYPES+1] = -2.0 * rho_bulk_Si/MW_Si * MW_H2 * d_etch_rate_d_C[1];
           break;

        case 3: /* SiO2OH2 */
           d_mass_flux[wspec][MAX_VARIABLE_TYPES+0] = -1.0 * rho_bulk_Si/MW_Si * MW_SiO2OH2 * d_etch_rate_d_C[0];
           d_mass_flux[wspec][MAX_VARIABLE_TYPES+1] = -1.0 * rho_bulk_Si/MW_Si * MW_SiO2OH2 * d_etch_rate_d_C[1];
           break;
       }

     d_mass_flux[wspec][MAX_VARIABLE_TYPES+2] = 0.0;
     d_mass_flux[wspec][MAX_VARIABLE_TYPES+3] = 0.0;
    }


  /* Load up func and d_func and return to apply_integrated BC */

  if (af->Assemble_Residual )
    {
     for (kdir = 0; kdir < dim; kdir++)
        {
          *func += s_terms.diff_flux[wspec][kdir] * fv->snormal[kdir];
        }

     *func -= mass_flux[wspec];
    }

  if (af->Assemble_Jacobian )
    {
      /* sum the contributions to the global stiffness matrix  for Species*/

      /*
       * J_s_c
       */
      var=MASS_FRACTION;
      for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
	phi_j = bf[var]->phi[j_id];
	for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++ )
	  {
           for (kdir = 0; kdir < dim; kdir++)
              {
               d_func[0][MAX_VARIABLE_TYPES + w1][j_id] +=
               s_terms.d_diff_flux_dc[wspec][kdir] [w1][j_id] * fv->snormal[kdir];
              }
	   d_func[0][MAX_VARIABLE_TYPES + w1][j_id] -=
	   d_mass_flux[wspec][MAX_VARIABLE_TYPES + w1] * phi_j;
	  }
      }
 
    } /* End of if Assemble_Jacobian */
  return;
} /* END of routine mass_flux_surf_etch_bc                                   */

/*****************************************************************************/
/****************************************************************************/

void
mass_flux_alloy_surf(double func[DIM],
		     double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
		     const int wspec,
		     const double p[],  /* parameters to parameterize heat transfer model*/
		     const dbl time)

/******************************************************************************
*
*  Function which calculates the surface integral for a cubic mass-flux transfer
*  model which fits metal evaporation rate curves well.  Basically the difference
*  between this model and the simple convective mass transfer coefficient is that the
*   coefficient is a cubic dependence on temperature.
*  Added by P. R. Schunk (2/9/01)
*
*  Functions called: none
*
*  Form of flux:  n.j_i = exp[c_0 + c_1*(T-Tm) - c_2*(T-Tm)^2 + c3*(T-Tm)^3]*(y_i-y_inf)
*
*  ----------------------------------------------------------------------------
*
******************************************************************************/
{
  
  int j, j_id;
  int var, w1;
  double phi_j;

  /* Execution begins */

  /* Now transfer into func, d_func */
  if (af->Assemble_Residual ) 
    {
      *func -= mp->mass_flux[wspec];
    }
  
  if (af->Assemble_Jacobian ) 
    {
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
	}
      }

    } /* End of if Assemble_Jacobian */
  return;
} /* END of routine mass_flux_alloy_surf                                   */

/*****************************************************************************/
void 
mass_flux_BV_surf_bc(
                double func[],
	        double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
		int wspec,     /* species number of this boundary condition */
		double nu,     /* stoichiometric coefficient                */
		double k,      /* kinetic rate constant                     */
		double beta,   /* reaction order                            */
		double alphaa, /* anodic direction transfer coefficient     */
		double alphac, /* cathodic direction transfer coefficient   */
		double V,      /* electrode potential                       */
		double U0,     /* theoretical open-circuit potential        */
		double dt,     /* current value of the time step            */
		double tt)     /* parameter to vary time integration        */

/*****************************************************************************
*
*  A function that calculates the mass flux rate (and associated Jacobian
*  entries) given by Butler-Volmer kinetics. This routine was cloned after
*  mass_flux_surf_bc. 
*
*  Ken S. Chen (11/2000)
*
*  Revised by KSC on 10/1/2001. 
*
*  Revised by KSC on 12/17/2001.
*  ---------------------------------------------------------------------------
*  Functions called:
*  mass_flux_surf_BV  -- calculates the total mass flux and associated
*                        sensititives for one component. 
*
*  ---------------------------------------------------------------------------
*
******************************************************************************/
{
  
  /* Local variables */
  
  int j, j_id,w1,dim,kdir;
  int var,jvar;
  double phi_j;
  double Y_w; /* local concentration of current species */
  double T = -1.0;   /* electrolyte solution temperature */
  double mass_flux[MAX_CONC];
  double d_mass_flux[MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC];

  double vconv[MAX_PDIM]; /*Calculated convection velocity */
  double vconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;

  int err=0;


  /***************************** EXECUTION BEGINS ****************************/

  /* first set the electrolyte temperature */
  if (mp->SolutionTemperatureModel == CONSTANT)
    {
      T = mp->solution_temperature;
    }
  else 
    {
      EH(-1, "Solution-temperature model not yet implemented");
    }
  if(pd->e[pg->imtrx][R_ENERGY]) /* if energy equation is active, re-set electrolyte temperature */
    {
      T = fv->T;   
    }
  /* set solution temperature to 298 K if it is zero - safety feature */
  if (T == 0.0) 
    {
      T = 298.0;  
      fprintf(stderr, "Warning!: a default electrolyte temperature of 298 K is being used!");
    }
  
  /* call routine to calculate surface flux of this component and its 
   * sensitivity to all variable types 
   */
  mass_flux_surf_BV (mass_flux, d_mass_flux, wspec, nu, k, beta, alphaa, alphac, V, U0, T); 

  dim   = pd->Num_Dim;
  Y_w = fv->c[wspec];

  /* get the convection velocity (it's different for arbitrary and
     lagrangian meshes) */

  /* get deformation gradient, etc. if lagrangian mesh with inertia */
  if ((pd->MeshMotion == LAGRANGIAN ||
       pd->MeshMotion == DYNAMIC_LAGRANGIAN) && pd->MeshInertia == 1)
    {
      err = belly_flop(elc->lame_mu);
      EH(err, "error in belly flop");
      if( neg_elem_volume ) return;
    }
  if (mp->PorousMediaType == CONTINUOUS) {
    err = get_convection_velocity(vconv, vconv_old, d_vconv, dt, tt);
  }
  EH(err, "Error in calculating effective convection velocity");

  if (af->Assemble_Residual )
    {
      *func -= mass_flux[wspec];
      
      /* Calculate the residual contribution from convective flux	*/
      for(kdir=0; kdir<dim; kdir++) 
	{
	  *func += Y_w*vconv[kdir] * fv->snormal[kdir];
	}
    }
  
  if (af->Assemble_Jacobian ) 
    {
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
	                    d_mass_flux[wspec][MAX_VARIABLE_TYPES + w1]*phi_j; 
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
	  
	  d_func[0][var][j_id] -= d_mass_flux[wspec][TEMPERATURE]*phi_j; 
	  
	  for(kdir=0; kdir<dim; kdir++) 
	    {
	      d_func[0][var][j_id]  += 
		          Y_w*d_vconv->T[kdir][j_id] * fv->snormal[kdir];
	    }
	}
      }

      /*
       * J_s_V
       */
      var=VOLTAGE;
      if (pd->v[pg->imtrx][var]){
	for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
	  j_id = j;
	  phi_j = bf[var]->phi[j_id];
	  
	  d_func[0][var][j_id] -= d_mass_flux[wspec][VOLTAGE]*phi_j; 
	  
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

return;
  
} /* END of routine mass_flux_BV_surf_bc                                    */
/****************************************************************************/

void
mass_flux_HOR_surf_bc(
                double func[],
                double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                int wspec,     /* species number of this boundary condition */
                double ai0,    /* product of interfacial area by            */
                               /* exchange current density (A/cm^3)         */
                double H,      /* thickness of catalyst layer (cm)          */
                double cref,   /* species ref. concentration (moles/cm^3)   */
                double alphaa, /* anodic direction transfer coefficient     */
                double alphac, /* cathodic direction transfer coefficient   */
                double T,      /* cell temperature (K)                      */
                double U0,     /* open-circuit potential (V)                */
                double beta,   /* reaction order                            */
                double n,      /* number of electrons involved in rxn       */
                double V,      /* electrode potential (V)                   */
                double dt,     /* current value of the time step            */
                double tt)     /* parameter to vary time integration        */

/*****************************************************************************
*
*  A function that calculates the surface mass flux rate (and associated Jacobian
*  entries) given by the linearized kinetic model for an electrochemical
* reaction such as the hydrogen oxidation reaction in PEM fuel cells.
*  (see Chen and Hickner 2006)
*
*  Ken S. Chen (1/2006)
*
*  ---------------------------------------------------------------------------
*  Functions called:
*  mass_flux_surf_HOR  -- calculates the total mass flux and associated
*                         sensititives for one component. 
*
*  ---------------------------------------------------------------------------
*
******************************************************************************/
{

  /* Local variables */

  int j, j_id, w1;
  int var;
  double phi_j;
  double mass_flux[MAX_CONC];
  double d_mass_flux[MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC];

  /***************************** EXECUTION BEGINS ****************************/

  /* call routine to calculate surface flux of this component and its
   * sensitivity to all variable types
   */

  mass_flux_surf_HOR (mass_flux, d_mass_flux, wspec, ai0, H, cref,
                      alphaa, alphac, T, U0, beta, n, V);

  /* residual */
  if (af->Assemble_Residual )
    {
      *func -= mass_flux[wspec];
    }

  /* sensitivity entries */
  if (af->Assemble_Jacobian )
    {

      /* J_s_c --- sensitivity wrt species concentrations */
      var=MASS_FRACTION;
      for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
        phi_j = bf[var]->phi[j_id];

       for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++ )
          {
           d_func[0][MAX_VARIABLE_TYPES + w1][j_id] = 0.0;
          }

          d_func[0][MAX_VARIABLE_TYPES][j_id] =
                    -d_mass_flux[wspec][MAX_VARIABLE_TYPES + wspec]*phi_j;
      }

      /* J_s_T --- sensitivity wrt temperature */
      var=TEMPERATURE;
      if (pd->v[pg->imtrx][var]){
       for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          j_id = j;
          phi_j = bf[var]->phi[j_id];

          d_func[0][var][j_id] = -d_mass_flux[wspec][TEMPERATURE]*phi_j;
       }
      }

      /* J_s_V --- sensitivity wrt voltage or potential */
      var=VOLTAGE;
      if (pd->v[pg->imtrx][var]){
       for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          j_id = j;
          phi_j = bf[var]->phi[j_id];

          d_func[0][var][j_id] = -d_mass_flux[wspec][VOLTAGE]*phi_j;
       }
      }

    } /* End of if Assemble_Jacobian */

  return;

} /* END of routine mass_flux_HOR_surf_bc                                    */
/****************************************************************************/

void
mass_flux_ORR_surf_bc(
                double func[],
                double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                int wspec,     /* species number of this boundary condition */
                double ai0,    /* product of interfacial area by            */
                               /* exchange current density (A/cm^3)         */
                double H,      /* thickness of catalyst layer (cm)          */
                double cref,   /* species ref. concentration (moles/cm^3)   */
                double alphac, /* cathodic direction transfer coefficient   */
                double T,      /* cell temperature (K)                      */
                double V,      /* cell voltage (K)                          */
                double U0,     /* open-circuit potential (V)                */
                double beta,   /* reaction order                            */
                double n,      /* number of electrons involved in rxn       */
                double dt,     /* current value of the time step            */
                double tt)     /* parameter to vary time integration        */

/*****************************************************************************
*
*  A function that calculates the surface mass flux rate (and associated Jacobian
*  entries) given by the Tafel kinetic model for an electrochemical
*  reaction such as the oxygen reduction reaction in PEM fuel cells
*  (see Chen and Hickner 2006).
*
*  Ken S. Chen (1/2006)
*  Modified: 5/22/2006 by KSC.
*
*  ---------------------------------------------------------------------------
*  Functions called:
*  mass_flux_surf_ORR  -- calculates the total mass flux and associated
*                         sensititives for one component.
*
*  ---------------------------------------------------------------------------
*
******************************************************************************/
{
  int j, j_id, w1;
  int var;
  double phi_j;
  double mass_flux[MAX_CONC];
  double d_mass_flux[MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC];

  /* call routine to calculate surface flux of this component and its 
   * sensitivity to all variable types 
   */

  mass_flux_surf_ORR (mass_flux, d_mass_flux, wspec, ai0, H,
                      cref, alphac, T, V, U0, beta, n);

  /* residual */
  if (af->Assemble_Residual )
    {
      *func -= mass_flux[wspec];
    }

  /* sensitivity entries */
  if (af->Assemble_Jacobian )
    {

      /* J_s_c --- sensitivity wrt species concentrations */
      var=MASS_FRACTION;
      for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
        phi_j = bf[var]->phi[j_id];

       for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++ )
          {
           d_func[0][MAX_VARIABLE_TYPES + w1][j_id] = 0.0;
          }

          d_func[0][MAX_VARIABLE_TYPES][j_id] =
                    -d_mass_flux[wspec][MAX_VARIABLE_TYPES + wspec]*phi_j;
      }

      /* J_s_T --- sensitivity wrt temperature */
      var=TEMPERATURE;
      if (pd->v[pg->imtrx][var]){
       for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          j_id = j;
          phi_j = bf[var]->phi[j_id];

          d_func[0][var][j_id] = -d_mass_flux[wspec][TEMPERATURE]*phi_j;
       }
      }

      /* J_s_V --- sensitivity wrt voltage or potential */
      var=VOLTAGE;
      if (pd->v[pg->imtrx][var]){
       for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          j_id = j;
          phi_j = bf[var]->phi[j_id];

          d_func[0][var][j_id] = -d_mass_flux[wspec][VOLTAGE]*phi_j;
       }
      }

    } /* End of if Assemble_Jacobian */

  return;

} /* END of routine mass_flux_ORR_surf_bc                                    */
/****************************************************************************/

void
mass_flux_H2O_ANODE_surf_bc(
                double func[],
                double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                int wspec,     /* species number of this boundary condition */
                double ai0,    /* product of interfacial area by            */
                               /* exchange current density (A/cm^3)         */
                double Ha,     /* thickness of anode catalyst layer (cm)    */
                double cH2ref, /* ref. concentration for H2 (moles/cm^3)    */
                double alphaa, /* anodic direction transfer coefficient     */
                double alphac, /* cathodic direction transfer coefficient   */
                double T,      /* cell temperature (K)                      */
                double U0,     /* Open-circuit potential for HOR (V)        */
                double nd,     /* electro-osmatic drag coefficient          */
                double dt,     /* current value of the time step            */
                double tt)     /* parameter to vary time integration        */

/*****************************************************************************
*
*  A function that calculates the surface mass flux rate (and associated Jacobian
*  entries) of H2O due to electro-osmatic drag.
*
*  Ken S. Chen (1/2006)
*
*  ---------------------------------------------------------------------------
*  Functions called:
*  mass_flux_surf_H2O_ANODE  -- calculates the total mass flux and associated
*                               sensititives for one component. 
*
*  ---------------------------------------------------------------------------
*
******************************************************************************/
{
  int j, j_id, w1;
  int var;
  double phi_j;
  double mass_flux[MAX_CONC];
  double d_mass_flux[MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC];

  /***************************** EXECUTION BEGINS ****************************/

  /* call routine to calculate surface flux of this component and its 
   * sensitivity to all variable types 
   */

  mass_flux_surf_H2O_ANODE (mass_flux, d_mass_flux, wspec, ai0, Ha, cH2ref, alphaa, alphac, T, U0, nd);

  if (af->Assemble_Residual )
    {
      *func -= mass_flux[wspec];
    }

  /* sensitivity entries */
  if (af->Assemble_Jacobian )
    {

      /* J_s_c --- sensitivity wrt species concentrations */
      var=MASS_FRACTION;
      for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
        phi_j = bf[var]->phi[j_id];

       for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++ )
          {
           d_func[0][MAX_VARIABLE_TYPES + w1][j_id] = 0.0;
          }

          d_func[0][MAX_VARIABLE_TYPES][j_id] =
                    -d_mass_flux[wspec][MAX_VARIABLE_TYPES + wspec]*phi_j;
      }

      /* J_s_T --- sensitivity wrt temperature */
      var=TEMPERATURE;
      if (pd->v[pg->imtrx][var]){
       for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          j_id = j;
          phi_j = bf[var]->phi[j_id];

          d_func[0][var][j_id] = -d_mass_flux[wspec][TEMPERATURE]*phi_j;
       }
      }

      /* J_s_V --- sensitivity wrt voltage or potential */
      var=VOLTAGE;
      if (pd->v[pg->imtrx][var]){
       for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          j_id = j;
          phi_j = bf[var]->phi[j_id];

          d_func[0][var][j_id] = -d_mass_flux[wspec][VOLTAGE]*phi_j;
       }
      }

    } /* End of if Assemble_Jacobian */


  return;

} /* END of routine mass_flux_H2O_ANODE_surf_bc                                    */
/****************************************************************************/

void
mass_flux_H2O_CATHODE_surf_bc(
                double func[],
                double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                int wspec,     /* species number of this boundary condition */
                double ai0,    /* product of interfacial area by            */
                               /* exchange current density (A/cm^3)         */
                double Hc,     /* thickness of cathode catalyst layer (cm)  */
                double cO2ref, /* ref. concentration for O2 (moles/cm^3)    */
                double alphac, /* cathodic direction transfer coefficient   */
                double T,      /* cell temperature (K)                      */
                double V,      /* cell voltage (K)                          */
                double U0,     /* Open-circuit potential for HOR (V)        */
                double nd,     /* electro-osmatic drag coefficient          */
                double dt,     /* current value of the time step            */
                double tt)     /* parameter to vary time integration        */

/*****************************************************************************
*
*  A function that calculates the surface mass flux rate (and associated Jacobian
*  entries) of H2O due to electro-osmatic drag and oxygen reduction
*  reaction.
*
*  Ken S. Chen (1/2006)
*
*  ---------------------------------------------------------------------------
*  Functions called:
*  mass_flux_surf_ORR  -- calculates the total mass flux and associated
*                         sensititives for one component.
*
*  ---------------------------------------------------------------------------
*
******************************************************************************/
{

  /* Local variables */

  int j, j_id, w1;
  int var;
  double phi_j;
  double mass_flux[MAX_CONC];
  double d_mass_flux[MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC];


  /***************************** EXECUTION BEGINS ****************************/

  /* call routine to calculate surface flux of this component and its 
   * sensitivity to all variable types 
   */

  mass_flux_surf_H2O_CATHODE (mass_flux, d_mass_flux, wspec, ai0, Hc, cO2ref, alphac, T, V, U0, nd);

  if (af->Assemble_Residual )
    {
      *func -= mass_flux[wspec];
    }

  /* sensitivity entries */
  if (af->Assemble_Jacobian )
    {

      /* J_s_c --- sensitivity wrt species concentrations */
      var=MASS_FRACTION;
      for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
        phi_j = bf[var]->phi[j_id];

       for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++ )
          {
           d_func[0][MAX_VARIABLE_TYPES + w1][j_id] = 0.0;
          }

          d_func[0][MAX_VARIABLE_TYPES][j_id] =
                    -d_mass_flux[wspec][MAX_VARIABLE_TYPES + wspec]*phi_j;
      }

      /* J_s_T --- sensitivity wrt temperature */
      var=TEMPERATURE;
      if (pd->v[pg->imtrx][var]){
       for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          j_id = j;
          phi_j = bf[var]->phi[j_id];

          d_func[0][var][j_id] = -d_mass_flux[wspec][TEMPERATURE]*phi_j;
       }
      }

      /* J_s_V --- sensitivity wrt voltage or potential */
      var=VOLTAGE;
      if (pd->v[pg->imtrx][var]){
       for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          j_id = j;
          phi_j = bf[var]->phi[j_id];

          d_func[0][var][j_id] = -d_mass_flux[wspec][VOLTAGE]*phi_j;
       }
      }

    } /* End of if Assemble_Jacobian */


  return;

} /* END of routine mass_flux_H2O_CATHODE_surf_bc                                    */
/****************************************************************************/

void 
mass_flux_SULFIDATION_surf_bc(
                double func[],
	        double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                int mode,         /* key word for sulfidation kinetic model    */
		int wspec,        /* species number of this boundary condition */
		double nu,        /* stoichiometric coefficient                */
		double k1,        /* forward kinetic rate constant             */
		double E1,        /* forward activation energy                 */
		double kn1,       /* backward kinetic rate constant            */
		double En1,       /* backward activation energy                */
		double T,         /* Temperature                               */
		double c_H2S,     /* bulk concentration of H2S                 */
		double c_O2,      /* bulk concentration of O2                  */
		double M_solid,   /* molecular weight of Cu2S                  */
		double rho_solid, /* theoretical open-circuit potential        */
		double dt,        /* current value of the time step            */
		double tt)        /* parameter to vary time integration        */

/*****************************************************************************
*
*  A function that calculates the mass flux rate (and associated Jacobian
*  entries) given by the various sulfidation kinetic models as selected by 
*  the user:
*
*
*  Simplified solid-diffusion-controlled with Cu as the diffusing species:
*
*  r = k exp(-E/R/T) c_H2S c_Cu         ( c_H2S is fixed at its bulk value ) 
*
*
*  Solid-diffusion-controlled with Cu vacancies as diffusing species and
*  the approximation of electroneutrality ( i.e., c_V = c_h ):
*
*  r = k1 exp(-E1/R/T) c_H2S c_O2**0.5 - kn1 exp(-En1/R/T) c_V**4 
*
*  ( c_H2S and c_O2 are taken to be fixed at its bulk values )
*
*
*  Solid-diffusion-controlled with Cu vacancies as diffusing species and
*  c_V being different from c_h near the Cu/Cu2S and gas/Cu2S interfaces:
*
*  r = k1 exp(-E1/R/T) c_H2S c_O2**0.5 - kn1 exp(-En1/R/T) c_V**2 c_h**2 
*
*  ( c_H2S and c_O2 are taken to be fixed at its bulk values )
*
*
*  Gas-diffusion-controlled with H2S bnd O2 being the diffusing species:
*
*  r = k1 exp(-E1/R/T) c_H2S c_O2**0.5     (c_V is taken to be zero)
*   
* 
*  This routine was cloned after mass_flux_BV_surf_bc.
*
*  Ken S. Chen (3/2002)
*
*  ---------------------------------------------------------------------------
*  Functions called:
*  mass_flux_surf_SULFIDATION  -- calculates the total mass flux and associated
*                                 sensititives for one component. 
*
*  ---------------------------------------------------------------------------
*
******************************************************************************/
{
  int j, j_id, w1, dim, kdir, var, jvar;
  double phi_j;
  double Y_w; /* local concentration of current species */
  double mass_flux[MAX_CONC];
  double d_mass_flux[MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC];

  double vconv[MAX_PDIM]; /*Calculated convection velocity */
  double vconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;

  int err=0;


  /***************************** EXECUTION BEGINS ****************************/

  /* call routine to calculate surface flux of this component and its
   * sensitivity to all variable types
   */
  mass_flux_surf_SULFIDATION (mass_flux, d_mass_flux, mode, wspec,
                              nu, k1, E1, kn1, En1, T, c_H2S, c_O2); 

  dim   = pd->Num_Dim;
  Y_w = fv->c[wspec];

  /* get the convection velocity (it's different for arbitrary and
     lagrangian meshes) */

  /* get deformation gradient, etc. if lagrangian mesh with inertia */
  if ((pd->MeshMotion == LAGRANGIAN ||
       pd->MeshMotion == DYNAMIC_LAGRANGIAN) && pd->MeshInertia == 1)
    {
      err = belly_flop(elc->lame_mu);
      EH(err, "error in belly flop");
      if( neg_elem_volume ) return;
    }
  if (mp->PorousMediaType == CONTINUOUS) {
    err = get_convection_velocity(vconv, vconv_old, d_vconv, dt, tt);
    EH(err, "Error in calculating effective convection velocity");
  }

  if (af->Assemble_Residual ) 
    {
      *func -= mass_flux[wspec];
      
      /* Calculate the residual contribution from convective flux	*/
      for(kdir=0; kdir<dim; kdir++) 
	{
	  *func += Y_w*vconv[kdir] * fv->snormal[kdir];
	}
    }
  
  if (af->Assemble_Jacobian ) 
    {
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
	                    d_mass_flux[wspec][MAX_VARIABLE_TYPES + w1]*phi_j; 
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
	  
	  d_func[0][var][j_id] -= d_mass_flux[wspec][TEMPERATURE]*phi_j; 
	  
	  for(kdir=0; kdir<dim; kdir++) 
	    {
	      d_func[0][var][j_id]  += 
		          Y_w*d_vconv->T[kdir][j_id] * fv->snormal[kdir];
	    }
	}
      }

      /*
       * J_s_V
       */
      var=VOLTAGE;
      if (pd->v[pg->imtrx][var]) {
	for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
	  phi_j = bf[var]->phi[j];
	  d_func[0][var][j] -= d_mass_flux[wspec][VOLTAGE] * phi_j;
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

return;
  
} /* END of routine mass_flux_SULFIDATION_surf_bc                           */
/******************************************************************************/

void
mass_flux_BV2_surf_bc(double func[],
                      double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                      int wspec,       /* species number of this boundary condition */
                      double k,        /* rate constant                             */
                      double PHI_E,    /* electrode potential                       */
                      double alphaa,   /* anodic transfer coefficient               */
                      double alphac,   /* cathodic transfer coefficient             */
                      double U0,       /* standard state open circuit potential     */
                      double T,        /* electrolyte temperature                   */
		      dbl time,        /* current value of the time                 */
                      dbl dt,          /* current value of the time step            */
                      dbl tt)          /* parameter to vary time integration        */

     /*******************************************************************************
      *
      *  Function which calculates the surface integral for a total mass flux
      *  given by Butler-Volmer kinetics.  This routine cloned from
      *  const_mass_flux_surf_bc by RSL on 1/10/01.
      *
      *  ----------------------------------------------------------------------------
      *
      *  Functions called:
      *  mass_flux_surf_BV2  -- returns mass flux for one component
      *
      *  ----------------------------------------------------------------------------
      *
      *******************************************************************************/
{ 
  int j, j_id, w1, dim, kdir, var, jvar, flag;
  double phi_j;
  double Y_w; /* local concentration of current species */
  double vconv[MAX_PDIM]; /*Calculated convection velocity */
  double vconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;
  dbl c, rho, M_mix;
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;  /* density dependence */
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;
  int err = 0;

  /***************************** EXECUTION BEGINS ****************************/

  /* calculate concentration to be used in the convective flux; this is the
   * actual total molar concentration if molar fluxes are being used, but
   * simply unity if volume fluxes are being used.
   */

  if ((mp->SpeciesSourceModel[wspec]  == ELECTRODE_KINETICS ||
       mp->SpeciesSourceModel[wspec]  == ION_REACTIONS) &&
      cr->MassFluxModel != STEFAN_MAXWELL_VOLUME) /*  RSL 3/19/01  */
     {
      M_mix = 0.;
      for (j=0; j<pd->Num_Species; j++)
         {
          M_mix += fv->c[j] * mp->molecular_weight[j];
         }
      rho = density(d_rho, time); /*  RSL 6/22/02  */
      c = rho/M_mix;
     }
  else
     {
      c = 1.;
     }

  /* call routine to return surface flux of this component and its
   * sensitivity to all variable types
   */
 
  flag = 0;
  mass_flux_surf_BV2(time, mp->mass_flux, mp->d_mass_flux, wspec, flag, k, PHI_E, alphaa, alphac
, U0, T);

  dim   = pd->Num_Dim;
  Y_w = fv->c[wspec];

  /* get deformation gradient, etc. if lagrangian mesh with inertia */

  if (pd->MeshMotion == LAGRANGIAN && pd->MeshInertia == 1) {
    err = belly_flop(elc->lame_mu);
    EH(err, "error in belly flop");
    if( neg_elem_volume ) return;
  }

  /* get the convection velocity (it's different for arbitrary and
     lagrangian meshes) */

  err = get_convection_velocity(vconv, vconv_old, d_vconv, dt, tt);
  EH(err, "Error in calculating effective convection velocity");

  if (af->Assemble_Residual) {
    *func -= mp->mass_flux[wspec];

    /* Calculate the residual contribution from convective flux */
    for (kdir = 0; kdir < dim; kdir++) {
      *func += c * Y_w * vconv[kdir] * fv->snormal[kdir];
    }
  }

  if (af->Assemble_Jacobian) {
    /* sum the contributions to the global stiffness matrix  for Species*/

    /*
     * J_s_c
     */
    var = MASS_FRACTION;
    for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
      phi_j = bf[var]->phi[j_id];
      for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++ ) {
        d_func[0][MAX_VARIABLE_TYPES + w1][j_id] -=
            mp->d_mass_flux[wspec][MAX_VARIABLE_TYPES + w1]
            * phi_j;
      }
      for (kdir = 0; kdir < dim; kdir++) {
        d_func[0][MAX_VARIABLE_TYPES + wspec][j_id] +=
            c*phi_j*vconv[kdir] * fv->snormal[kdir];

        for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++ ) {
          d_func[0][MAX_VARIABLE_TYPES + w1][j_id] +=
              c*Y_w*d_vconv->C[kdir][w1][j_id] * fv->snormal[kdir];
        }
      }
    }
     
      /*
       * J_s_T
       */
    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        j_id = j;
        phi_j = bf[var]->phi[j_id];
        d_func[0][var][j_id] -= (mp->d_mass_flux[wspec][TEMPERATURE]
                                 * phi_j);
        for(kdir=0; kdir<dim; kdir++)  {
          d_func[0][var][j_id]  +=
              c*Y_w*d_vconv->T[kdir][j_id] * fv->snormal[kdir];
        }
      }
    }
     
      /*
       * J_s_V
       */
    var = VOLTAGE;
    if (pd->v[pg->imtrx][var])
       {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
           {
            j_id = j;
            phi_j = bf[var]->phi[j_id];
            d_func[0][var][j_id] -= (mp->d_mass_flux[wspec][VOLTAGE] * phi_j);
           }
       }

      /*
       * J_s_d
       */
    for (jvar = 0; jvar < dim; jvar++)  {
      var = MESH_DISPLACEMENT1 + jvar;
      if (pd->v[pg->imtrx][var]) {
        for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
          /*     d( )/dx        */
          /* additional terms due to convective flux */

          phi_j = bf[var]->phi[j_id];
          for (kdir = 0; kdir < dim; kdir++) {
            d_func[0][var][j_id] += c*(
                Y_w*vconv[kdir]*fv->dsnormal_dx[kdir][jvar][j_id] +
                Y_w*d_vconv->X[kdir][jvar][j_id]*fv->snormal[kdir]);
          }
        }
      }
    }
     
    for (jvar = 0; jvar < dim; jvar++) {
      var = VELOCITY1 + jvar;
      if (pd->v[pg->imtrx][var]) {
        for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++)  {
          phi_j = bf[var]->phi[j_id];
          d_func[0][var][j_id] +=
              c*Y_w*d_vconv->v[jvar][jvar][j_id]*fv->snormal[jvar];
        }
      }
    }
  } /* End of if Assemble_Jacobian */
} /* END of routine mass_flux_BV2_surf_bc         */

/******************************************************************************/
/******************************************************************************/

void
mass_flux_NI_surf_bc(double func[],
                     double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                     int wspec,       /* species number of this boundary condition */
                     double PHI_E,    /* electrode potential                       */
                     double T,        /* electrolyte temperature                   */
		     double time,     /* current time                              */
                     dbl dt,          /* current value of the time step            */
                     dbl tt)          /* parameter to vary time integration        */

     /*******************************************************************************
      *
      *  Function which calculates the surface integral for a total mass flux
      *  in Ni electroplating.  This routine cloned from
      *  mass_flux_BV2_surf_bc by RSL on 3/9/01.
      *
      *  ----------------------------------------------------------------------------
      *
      *  Functions called:
      *  mass_flux_surf_NI  -- returns mass flux for one component
      *
      *  ----------------------------------------------------------------------------
      *
      *******************************************************************************/
{ 
  int j, j_id, w1, dim, kdir, var, jvar, flag;
  double phi_j;
  double Y_w; /* local concentration of current species */
  double vconv[MAX_PDIM]; /*Calculated convection velocity */
  double vconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;
  dbl c, rho, M_mix;
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;  /* density dependence */
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;
  int err = 0;

  /***************************** EXECUTION BEGINS ****************************/

  /* calculate total molar concentration to be used in the convective flux */

  M_mix = 0.;
  for (j=0; j<pd->Num_Species; j++)
     {
      M_mix += fv->c[j] * mp->molecular_weight[j];
     }
  rho = density(d_rho, time); /*  RSL 6/22/02  */
  c = rho/M_mix;

  /* call routine to return surface flux of this component and its
   * sensitivity to all variable types
   */
 
  flag = 0;
  mass_flux_surf_NI(mp->mass_flux, mp->d_mass_flux, time, wspec, flag, PHI_E, T);

  dim   = pd->Num_Dim;
  Y_w = fv->c[wspec];

  /* get deformation gradient, etc. if lagrangian mesh with inertia */

  if (pd->MeshMotion == LAGRANGIAN && pd->MeshInertia == 1) {
    err = belly_flop(elc->lame_mu);
    EH(err, "error in belly flop");
    if( neg_elem_volume ) return;
  }

  /* get the convection velocity (it's different for arbitrary and
     lagrangian meshes) */

  err = get_convection_velocity(vconv, vconv_old, d_vconv, dt, tt);
  EH(err, "Error in calculating effective convection velocity");

  if (af->Assemble_Residual) {
    *func -= mp->mass_flux[wspec];

    /* Calculate the residual contribution from convective flux */
    for (kdir = 0; kdir < dim; kdir++) {
      *func += c * Y_w * vconv[kdir] * fv->snormal[kdir];
    }
  }

  if (af->Assemble_Jacobian) {
    /* sum the contributions to the global stiffness matrix  for Species*/

    /*
     * J_s_c
     */
    var = MASS_FRACTION;
    for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
      phi_j = bf[var]->phi[j_id];
      for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++ ) {
        d_func[0][MAX_VARIABLE_TYPES + w1][j_id] -=
            mp->d_mass_flux[wspec][MAX_VARIABLE_TYPES + w1]
            * phi_j;
      }
      for (kdir = 0; kdir < dim; kdir++) {
        d_func[0][MAX_VARIABLE_TYPES + wspec][j_id] +=
            c*phi_j*vconv[kdir] * fv->snormal[kdir];

        for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++ ) {
          d_func[0][MAX_VARIABLE_TYPES + w1][j_id] +=
              c*Y_w*d_vconv->C[kdir][w1][j_id] * fv->snormal[kdir];
        }
      }
    }
     
      /*
       * J_s_T
       */
    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        j_id = j;
        phi_j = bf[var]->phi[j_id];
        d_func[0][var][j_id] -= (mp->d_mass_flux[wspec][TEMPERATURE]
                                 * phi_j);
        for(kdir=0; kdir<dim; kdir++)  {
          d_func[0][var][j_id]  +=
              c*Y_w*d_vconv->T[kdir][j_id] * fv->snormal[kdir];
        }
      }
    }
     
      /*
       * J_s_V
       */
    var = VOLTAGE;
    if (pd->v[pg->imtrx][var])
       {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
           {
            j_id = j;
            phi_j = bf[var]->phi[j_id];
            d_func[0][var][j_id] -= (mp->d_mass_flux[wspec][VOLTAGE] * phi_j);
           }
       }

      /*
       * J_s_d
       */
    for (jvar = 0; jvar < dim; jvar++)  {
      var = MESH_DISPLACEMENT1 + jvar;
      if (pd->v[pg->imtrx][var]) {
        for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
          /*     d( )/dx        */
          /* additional terms due to convective flux */

          phi_j = bf[var]->phi[j_id];
          for (kdir = 0; kdir < dim; kdir++) {
            d_func[0][var][j_id] += c*(
                Y_w*vconv[kdir]*fv->dsnormal_dx[kdir][jvar][j_id] +
                Y_w*d_vconv->X[kdir][jvar][j_id]*fv->snormal[kdir]);
          }
        }
      }
    }

    for (jvar = 0; jvar < dim; jvar++) {
      var = VELOCITY1 + jvar;
      if (pd->v[pg->imtrx][var]) {
        for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++)  {
          phi_j = bf[var]->phi[j_id];
          d_func[0][var][j_id] +=
              c*Y_w*d_vconv->v[jvar][jvar][j_id]*fv->snormal[jvar];
        }
      }
    }
  } /* End of if Assemble_Jacobian */
} /* END of routine mass_flux_NI_surf_bc         */

/******************************************************************************/
/******************************************************************************/

void
current_BV2_surf_bc(double func[],
                    double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                    int wspec,       /* species number of this boundary condition */
                    double k,        /* rate constant                             */
                    double PHI_E,    /* electrode potential                       */
                    double alphaa,   /* anodic transfer coefficient               */
                    double alphac,   /* cathodic transfer coefficient             */
                    double U0,       /* standard state open circuit potential     */
                    double T,        /* electrolyte temperature                   */
                    dbl time,        /* current time                              */
                    dbl dt,          /* current value of the time step            */
                    dbl tt)          /* parameter to vary time integration        */

     /*******************************************************************************
      *
      *  Function which calculates the surface integral for an electric current
      *  given by Butler-Volmer kinetics.  This routine cloned from
      *  mass_flux_BV2_surf_bc by RSL on 1/22/01.  However, note that convective
      *  terms are not included here, because convection cannot contribute to the
      *  current in an electrically neutral phase.
      *
      *  ----------------------------------------------------------------------------
      *
      *  Functions called:
      *  mass_flux_surf_BV2  -- returns current for the component in question
      *
      *  ----------------------------------------------------------------------------
      *
      *******************************************************************************/
{
  int j, j_id, w1, var, flag, store;
  double phi_j;

  /***************************** EXECUTION BEGINS ****************************/

  /* call routine to return current for this component and its
   * sensitivity to all variable types
   */

  flag = 1;
  mass_flux_surf_BV2(time, mp->mass_flux, mp->d_mass_flux, wspec, flag, k, PHI_E, alphaa, alphac, U0, T);
  store = MAX_CONC - 1;

  if (af->Assemble_Residual) {
    *func -= mp->mass_flux[store];
  }

  if (af->Assemble_Jacobian) {
    /* sum the contributions to the global stiffness matrix for species */

    /*
     * J_s_c
     */
    var = MASS_FRACTION;
    for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++)
       {
        phi_j = bf[var]->phi[j_id];
        for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++ )
           {
            d_func[0][MAX_VARIABLE_TYPES + w1][j_id] -=
            mp->d_mass_flux[store][MAX_VARIABLE_TYPES + w1] * phi_j;
           }
       }

    /*
     * J_s_T
     */
    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var])
       {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
           {
            j_id = j;
            phi_j = bf[var]->phi[j_id];
            d_func[0][var][j_id] -= (mp->d_mass_flux[store][TEMPERATURE] * phi_j);
           }
       }

    /*
     * J_s_V
     */
    var = VOLTAGE;
    if (pd->v[pg->imtrx][var])
       {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
           {
            j_id = j;
            phi_j = bf[var]->phi[j_id];
            d_func[0][var][j_id] -= (mp->d_mass_flux[store][VOLTAGE] * phi_j);
           }
       }

  } /* End of if Assemble_Jacobian */
} /* END of routine current_BV2_surf_bc         */

/******************************************************************************/
/******************************************************************************/

void
current_NI_surf_bc(double func[],
                   double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                   int wspec,       /* dummy variable -- not used here           */
                   double PHI_E,    /* electrode potential                       */
                   double T,        /* electrolyte temperature                   */
                   dbl time,        /* current value of the time                 */
                   dbl dt,          /* current value of the time step            */
                   dbl tt)          /* parameter to vary time integration        */

     /*******************************************************************************
      *
      *  Function which calculates the surface integral for an electric current
      *  in Ni electroplating.  This routine cloned from
      *  current_BV2_surf_bc by RSL on 3/9/01.
      *
      *  ----------------------------------------------------------------------------
      *
      *  Functions called:
      *  mass_flux_surf_NI  -- returns total current from all contributing reactions
      *
      *  ----------------------------------------------------------------------------
      *
      *******************************************************************************/
{
  int j, j_id, w1, var, flag, store;
  double phi_j;

  /***************************** EXECUTION BEGINS ****************************/

  /* call routine to return the total current and its sensitivity to all variable types */

  flag = 1;
  mass_flux_surf_NI(mp->mass_flux, mp->d_mass_flux, time, wspec, flag, PHI_E, T);
  store = MAX_CONC - 1;

  if (af->Assemble_Residual) {
    *func -= mp->mass_flux[store];
  }

  if (af->Assemble_Jacobian) {
    /* sum the contributions to the global stiffness matrix for species */

    /*
     * J_s_c
     */
    var = MASS_FRACTION;
    for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++)
       {
        phi_j = bf[var]->phi[j_id];
        for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++ )
           {
            d_func[0][MAX_VARIABLE_TYPES + w1][j_id] -=
            mp->d_mass_flux[store][MAX_VARIABLE_TYPES + w1] * phi_j;
           }
       }

    /*
     * J_s_T
     */
    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var])
       {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
           {
            j_id = j;
            phi_j = bf[var]->phi[j_id];
            d_func[0][var][j_id] -= (mp->d_mass_flux[store][TEMPERATURE] * phi_j);
           }
       }

    /*
     * J_s_V
     */
    var = VOLTAGE;
    if (pd->v[pg->imtrx][var])
       {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
           {
            j_id = j;
            phi_j = bf[var]->phi[j_id];
            d_func[0][var][j_id] -= (mp->d_mass_flux[store][VOLTAGE] * phi_j);
           }
       }

  } /* End of if Assemble_Jacobian */
} /* END of routine current_NI_surf_bc         */

/******************************************************************************/
/****************************************************************************/

void
const_mass_flux_surf_bc(
                  double func[],
		  double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
		  int wspec,	/* species number of this boundary condition */
		  double const_mass_flux, /* specified flux       */
		  double time,  /* current value of the time                 */
		  dbl dt,	/* current value of the time step            */
		  dbl tt)	/* parameter to vary time integration        */

/*****************************************************************************
*
*  Function which calculates the surface integral for a constant
*  total mass flux.  This routine cloned from mass_flux_surf_bc
*  by RSL on 6/8/00.  Admittedly this involves some inefficiency.
*
*  Revised:  RSL 9/20/00 to allow for molar (rather than volume) flux.
*
*  ---------------------------------------------------------------------------
*
*  Functions called:
*  mass_flux_surf_const  -- returns mass flux for one component
*
*  ---------------------------------------------------------------------------
*
******************************************************************************/
{
  int j, j_id, w1, dim, kdir, var, jvar;
  double phi_j;
  double Y_w; /* local concentration of current species */
  double vconv[MAX_PDIM]; /*Calculated convection velocity */
  double vconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;
  dbl c, rho, M_mix;
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;  /* density dependence */
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;
  int err = 0;
 
  /***************************** EXECUTION BEGINS ****************************/
 
  /* calculate concentration to be used in the convective flux; this is the
   * actual total molar concentration if molar fluxes are being used, but
   * simply unity if volume fluxes are being used.
   */

  if ((mp->SpeciesSourceModel[wspec]  == ELECTRODE_KINETICS ||
       mp->SpeciesSourceModel[wspec]  == ION_REACTIONS) &&
      cr->MassFluxModel != STEFAN_MAXWELL_VOLUME) /*  RSL 3/19/01  */
     {
      M_mix = 0.;
      for (j=0; j<pd->Num_Species; j++)
         {
          M_mix += fv->c[j] * mp->molecular_weight[j];
         }
      rho = density(d_rho, time); /*  RSL 6/22/02  */
      c = rho/M_mix;
     }
  else
     {
      c = 1.;
     }

  /* call routine to return surface flux of this component and its
   * sensitivity to all variable types
   */
 
  mass_flux_surf_const(mp->mass_flux, mp->d_mass_flux, fv->T,
                       fv->c, wspec, const_mass_flux);

  dim   = pd->Num_Dim;
  Y_w = fv->c[wspec];

  /* get deformation gradient, etc. if lagrangian mesh with inertia */

  if (pd->MeshMotion == LAGRANGIAN && pd->MeshInertia == 1) {
    err = belly_flop(elc->lame_mu);
    EH(err, "error in belly flop");
    if( neg_elem_volume ) return;
  }

  /* get the convection velocity (it's different for arbitrary and
     lagrangian meshes) */

  err = get_convection_velocity(vconv, vconv_old, d_vconv, dt, tt);
  EH(err, "Error in calculating effective convection velocity");

  if (af->Assemble_Residual) {
    *func -= mp->mass_flux[wspec];

    /* Calculate the residual contribution from convective flux */
    for (kdir = 0; kdir < dim; kdir++) {
      *func += c * Y_w * vconv[kdir] * fv->snormal[kdir];
    }
  }

  if (af->Assemble_Jacobian) {
    /* sum the contributions to the global stiffness matrix  for Species*/

    /*
     * J_s_c
     */
    var = MASS_FRACTION;
    for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
      phi_j = bf[var]->phi[j_id];
      for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++ ) {
        d_func[0][MAX_VARIABLE_TYPES + w1][j_id] -=
            mp->d_mass_flux[wspec][MAX_VARIABLE_TYPES + w1]
            * phi_j;
      }
      for (kdir = 0; kdir < dim; kdir++) {
        d_func[0][MAX_VARIABLE_TYPES + wspec][j_id] +=
            c*phi_j*vconv[kdir] * fv->snormal[kdir];

        for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++ ) {
          d_func[0][MAX_VARIABLE_TYPES + w1][j_id] +=
              c*Y_w*d_vconv->C[kdir][w1][j_id] * fv->snormal[kdir];
        }
      }
    }
     
    /*
       * J_s_T
       */
    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        j_id = j;
        phi_j = bf[var]->phi[j_id];
        d_func[0][var][j_id] -= (mp->d_mass_flux[wspec][TEMPERATURE]
                                 * phi_j);
        for(kdir=0; kdir<dim; kdir++)  {
          d_func[0][var][j_id]  +=
              c*Y_w*d_vconv->T[kdir][j_id] * fv->snormal[kdir];
        }
      }
    }
     
    /*
       * J_s_d
       */
    for (jvar = 0; jvar < dim; jvar++)  {
      var = MESH_DISPLACEMENT1 + jvar;
      if (pd->v[pg->imtrx][var]) {
        for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
          /*     d( )/dx        */
          /* additional terms due to convective flux */
               
          phi_j = bf[var]->phi[j_id];
          for (kdir = 0; kdir < dim; kdir++) {
            d_func[0][var][j_id] += c*(
                Y_w*vconv[kdir]*fv->dsnormal_dx[kdir][jvar][j_id] +
                Y_w*d_vconv->X[kdir][jvar][j_id]*fv->snormal[kdir]);
          }
        }
      }
    }
     
    for (jvar = 0; jvar < dim; jvar++) {
      var = VELOCITY1 + jvar;
      if (pd->v[pg->imtrx][var]) {
        for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++)  {
          phi_j = bf[var]->phi[j_id];
          d_func[0][var][j_id] +=
              c*Y_w*d_vconv->v[jvar][jvar][j_id]*fv->snormal[jvar];
        }
      }
    }
  } /* End of if Assemble_Jacobian */
} /* END of routine const_mass_flux_surf_bc         */

/******************************************************************************/
/******************************************************************************
*  Function which calculates the mass flux rate for convective
*  mass transfer using a mass transfer coefficient. Equilibrium relations
*  are utilized to obtained the concentration at the surface.  Thus far,
*  two options are allowed for calculating activities of liquid components.
*  Raoult's law and Flory Huggins for solvent-polymermixture.
*
*           Author: A.C. Sun (9/98)
*
******************************************************************************/

void 
mass_flux_equil_mtc(dbl mass_flux[MAX_CONC],
		    dbl d_mass_flux[MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC],
		    double activity[MAX_CONC],
		    double dact_dC[MAX_CONC][MAX_CONC],
		    double y_mass[MAX_CONC], /* conc at boundary             */
		    int mode,	/* model to which the VLE is based on        */
		    double amb_pres, /* ambient pressure                     */
		    int wspec,	/* species no.                               */
		    double mass_tran_coeff, /* MASS transfer coeff           */
 		    double d_mtc[MAX_VARIABLE_TYPES+MAX_CONC],
		    double Y_c)	/* bulk concentration 	                     */
/*****************************************************************************/
{
   /* Local variables */
     int w;
     double mtc;
     double psat[MAX_CONC], dpsatdt[MAX_CONC];
     double A ;

     double P_diff_1, P_diff_2, P_diff_log;
     double P_diff_1_dT, P_diff_log_dT, P_diff_log_dw[MAX_CONC], P_diff_1_dw;
     double P_bulk=0.;

     int i, j, jac, k, mn;
     int Num_S1;
     double flory1,flory2,flory3,flory;
     double df1_dc[MAX_CONC],df2_dc[MAX_CONC];
     double df3_dc[MAX_CONC],df_dc[MAX_CONC];
     double truedf_dc[MAX_CONC],dv_dw[MAX_CONC][MAX_CONC];
     double C[MAX_CONC], vol[MAX_CONC], sv[MAX_CONC];
     double mw[MAX_CONC], prod[MAX_CONC]; 
     double bottom, prod2, sum_C;
     double chi[MAX_CONC][MAX_CONC]; /* chi is the binary interaction parameter*/
     double mw_last=0; /* Molecular weight of non-condensable and conversion factor */
     
     if (MAX_CONC < 3) {
       EH(-1, "mass_flux_equil_mtc expects MAX_CONC >= 3");
       return;
     }
     memset(prod, 0, sizeof(double)*MAX_CONC);

/***************************** EXECUTION BEGINS *******************************/

  /* Insure concentration are positive */
     for( j=0; j<pd->Num_Species_Eqn; j++  ) y_mass[j] = MAX(DBL_SMALL,y_mass[j]);
     for( j=0; j<pd->Num_Species_Eqn; j++  ) C[j] = y_mass[j];

  /* Nonideal VP Calculations based on either ANTOINE or RIEDEL models */

  if(mp->VaporPressureModel[wspec] == ANTOINE )
     {
       antoine_psat(wspec, mp->u_vapor_pressure[wspec],
                    &psat[wspec], &dpsatdt[wspec]);
       mp-> vapor_pressure[wspec] = psat[wspec];
     }

  else if(mp->VaporPressureModel[wspec] == RIEDEL )
     {
       riedel_psat(wspec, mp->u_vapor_pressure[wspec],
                   &psat[wspec], &dpsatdt[wspec]);
       mp-> vapor_pressure[wspec] = psat[wspec];
     }

  A = mp->vapor_pressure[wspec]/amb_pres;

  if(mode==RAOULT)
    {
      bottom=0.;
      prod2=1.;

      for (i = 0; i<pd->Num_Species_Eqn; i++) 
	{
	  if(mp->molecular_weight[i] < 0.)
	    {
	      EH(-1, "Molecular weight of the species not specified in the material file");
	    }
	  else
	    {
	      mw[i]  = mp->molecular_weight[i];
	    }
	}

      if(mp->molecular_weight[pd->Num_Species_Eqn] < 0.)
	{
	  EH(-1, "M.W. of a non-condensable species not specified in the material file");
	}
      else
	{
	  mw_last = mp -> molecular_weight[pd->Num_Species_Eqn];
	}
      
      for(i=0;i<pd->Num_Species_Eqn;i++)
	 {
	   prod[i]=1.;
	   for(j=0;j<pd->Num_Species_Eqn;j++)
	     {
	       prod[i] *= (i==j? 1. :mw[j]);
	     } 
	   prod2 *= mw[i];
	 }

      for(i=0;i<pd->Num_Species_Eqn ;i++)
	 {
	   bottom += (mw_last*prod[i]-prod2)*C[i];
	 }
      bottom += prod2;
   
      activity[wspec] = mw_last*prod[wspec]*C[wspec]/bottom;

      /* write out the derivatives of mol_frac w.r.t. mass_frac*/
      
      for(j=0;j<pd->Num_Species_Eqn;j++)
	 {
	   if (j==wspec)
	     {
	       dact_dC[wspec][j] = mw_last*prod[wspec]/bottom
                -mw_last*prod[wspec]*C[wspec]*(mw_last*prod[j]-prod2)
                /(bottom*bottom);

	     }
	   else
	     {
	       dact_dC[wspec][j] = -mw_last*prod[wspec]
		 *C[wspec]*(mw_last*prod[j]-prod2)/(bottom*bottom) ;
	     }
	 }
    }

  else if(mode==FLORY || mode==FLORY_CC)

    {
  /* Define some convenient/repetitive chunks to make eqns more compact */

      Num_S1 = pd->Num_Species_Eqn + 1;
      
      memset(df_dc, 0,sizeof(double)*MAX_CONC);
      memset(truedf_dc, 0,sizeof(double)*MAX_CONC);
      memset(df1_dc, 0,sizeof(double)*MAX_CONC);
      memset(df2_dc, 0,sizeof(double)*MAX_CONC);
      memset(df3_dc, 0,sizeof(double)*MAX_CONC);
      memset(dv_dw, 0,sizeof(double)*MAX_CONC*MAX_CONC);
      flory = flory1 = flory2 = flory3 = 0;
  

       for (i = 0; i<Num_S1; i++) 
	 {
	   if(mp->specific_volume[i] < 0.)
	     {
	       EH(-1, "Specific volume not specified in the material file.");
	     }
	   else
	     {
	       sv[i] = mp->specific_volume[i];
	     }
	   if(mp->molar_volume[i] < 0.)
	     {
	       EH(-1, "Molar volume not specified in the material file");
	     }
	   else
	     {
	       vol[i] = mp->molar_volume[i];
	     }
	   for(k=0;k<Num_S1;k++)
	     {
	     if(mp->flory_param[i][k] < 0.)
	       {
		 EH(-1, "Flory-Huggins binary parameters not specified in the material file");
	       }
	     else
	       {
		 chi[i][k] = mp->flory_param[i][k];
	       }
	     }
	 }

       bottom = 0.;
       sum_C = 0.;
/*  denominators for fraction-based formulations */
	if(mp->Species_Var_Type == SPECIES_MASS_FRACTION
             || mp->Species_Var_Type == SPECIES_UNDEFINED_FORM)
	{
         for (i = 0; i<pd->Num_Species_Eqn; i++) 
	    {	   
	     bottom += y_mass[i]*(sv[i]-sv[pd->Num_Species_Eqn]);
	    }
         bottom += sv[pd->Num_Species_Eqn];
	}
	if(mp->Species_Var_Type == SPECIES_MOLE_FRACTION)
	{
         for (i = 0; i<pd->Num_Species_Eqn; i++) 
	    {	   
	     bottom += y_mass[i]*(vol[i]-vol[pd->Num_Species_Eqn]);
	    }
         bottom += vol[pd->Num_Species_Eqn];
	}
	   
	if(mp->Species_Var_Type == SPECIES_DENSITY)
		{
                 for(i=0;i<pd->Num_Species_Eqn;i++)
	           {
	            C[i] = y_mass[i]*sv[i];
	            sum_C += C[i];
		    dv_dw[i][i]= sv[i];
	           }
		}
	 else if(mp->Species_Var_Type == SPECIES_CONCENTRATION)
		{
                 for(i=0;i<pd->Num_Species_Eqn;i++)
	           {
	            C[i] = y_mass[i]*vol[i];
	            sum_C += C[i];
	            dv_dw[i][i]= vol[i];
	           }
		}
	else if(mp->Species_Var_Type == SPECIES_MASS_FRACTION
             || mp->Species_Var_Type == SPECIES_UNDEFINED_FORM)
		{
                 for(i=0;i<pd->Num_Species_Eqn;i++)
	           {
	            C[i] = y_mass[i]*sv[i]/bottom;
	            sum_C += C[i];
	            for(j=0;j<pd->Num_Species_Eqn;j++)
	               {
                        dv_dw[i][j] = -y_mass[i]*sv[i]
		              *(sv[j]-sv[pd->Num_Species_Eqn])/SQUARE(bottom); 
                       }
		    dv_dw[i][i] += sv[i]/bottom; 
                   }
		}
	else if(mp->Species_Var_Type == SPECIES_MOLE_FRACTION)
		{
                 for(i=0;i<pd->Num_Species_Eqn;i++)
	           {
	            C[i] = y_mass[i]*vol[i]/bottom;
	            sum_C += C[i];
	            for(j=0;j<pd->Num_Species_Eqn;j++)
	               {
                        dv_dw[i][j] = -y_mass[i]*vol[i]
		              *(vol[j]-vol[pd->Num_Species_Eqn])/SQUARE(bottom); 
                       }
		    dv_dw[i][i] += vol[i]/bottom; 
		   }
		}
	else
		{ EH(-1,"That species formulation not done in mtc_flory\n"); }
	   
       
       for(k=0;k<pd->Num_Species_Eqn; k++)
	 {
	   flory1 += vol[wspec]*(delta(k,wspec)-C[k])/vol[k]
	             +vol[wspec]*C[k]/vol[pd->Num_Species_Eqn];

	   df1_dc[k]= (vol[wspec]/vol[pd->Num_Species_Eqn])
	              - (vol[wspec]/vol[k]) ;
	 }
       flory1 += -vol[wspec]/vol[pd->Num_Species_Eqn];

       if (pd->Num_Species_Eqn > 1)
	 {
	   for(k=1;k<pd->Num_Species_Eqn; k++)
	     {
	       for(jac=0;jac<k;jac++)
		 {   
		   flory2 += (delta(wspec,jac)*C[k]+(vol[wspec]/vol[jac])
			       *C[jac]*(delta(k,wspec)-C[k]))*chi[jac][k];
	       df2_dc[k] += (delta(wspec,jac) 
                                 - vol[wspec]/vol[jac]*C[jac])*chi[jac][k];
	       df2_dc[jac] += (vol[wspec]/vol[jac])*(delta(k,wspec)-C[k])
                                      *chi[jac][k];
		 }
	     }

	 }

       for(jac=0;jac<pd->Num_Species_Eqn;jac++)
	 {
	   flory3 +=(delta(wspec,jac)-(vol[wspec]/vol[jac])*C[jac])
	           *(1.-sum_C)*chi[jac][pd->Num_Species_Eqn];
           df3_dc[jac] -= (vol[wspec]/vol[jac])*(1.-sum_C)
                            *chi[jac][pd->Num_Species_Eqn];

	   for(mn=0;mn<pd->Num_Species_Eqn; mn++)
	     {
	       df3_dc[mn]
		 -= (delta(wspec,jac)-(vol[wspec]/vol[jac])*C[jac])
		       *chi[mn][pd->Num_Species_Eqn];
	     }
	 }
       flory= flory1+flory2+flory3;

       /* check the simplest case: 1solvent, 1polymer 
       check=(1-sum_C)+chi[0][1]*(1.-sum_C)*(1.-sum_C);
       printf("flory = %e, check =%e\n", flory[i],check); */

       for(k=0;k<pd->Num_Species_Eqn;k++)
         {
           df_dc[k]= df1_dc[k]+df2_dc[k]+df3_dc[k];
	 }
      
       for(i=0;i<pd->Num_Species_Eqn;i++)
	 {
	   for(j=0;j<pd->Num_Species_Eqn;j++)
	     {
	       truedf_dc[i] += df_dc[j]*dv_dw[j][i];
	     }
	 }

      activity[wspec] = C[wspec]*exp(flory);

      for (i = 0; i<pd->Num_Species_Eqn; i++)
       {
	dact_dC[wspec][i] = activity[wspec]*truedf_dc[i]; 
       }
	dact_dC[wspec][wspec] += dv_dw[wspec][wspec]*exp(flory);

/**   log_mean pressure difference for Chilton_Coburn  **/

  if(mode == FLORY_CC)
        {
        P_diff_1 = amb_pres - amb_pres*A*activity[wspec];
        P_diff_2 = amb_pres - P_bulk;
	P_diff_1_dT = -dpsatdt[wspec]*activity[wspec];
                if( P_diff_1 >= P_diff_2 || P_diff_1 <= 0.0)
                        {
			 P_diff_log = 1.0;
			 P_diff_log_dT = 0.0;
    			 for (w=0; w<pd->Num_Species_Eqn; w++) 
				{
				 P_diff_log_dw[w]=0.;
				}
			}
                        else
                        {
			 P_diff_log = (P_diff_2 - P_diff_1)/log(P_diff_2/P_diff_1);
			 P_diff_log_dT = 
				(P_diff_log/P_diff_1-1.)*P_diff_1_dT/
				log(P_diff_2/P_diff_1);
    			 for (w=0; w<pd->Num_Species_Eqn; w++) 
				{
			 	P_diff_1_dw = -amb_pres*A*dact_dC[wspec][w];
			 	P_diff_log_dw[w] = 
					(P_diff_log/P_diff_1-1.)*P_diff_1_dw/
					log(P_diff_2/P_diff_1);
				}
			}
        d_mtc[TEMPERATURE] = amb_pres*
		(P_diff_log*d_mtc[TEMPERATURE] - mass_tran_coeff*P_diff_log_dT)
			/SQUARE(P_diff_log);
    	for (w=0; w<pd->Num_Species_Eqn; w++) 
		{
        	d_mtc[MAX_VARIABLE_TYPES+w] = amb_pres*
			(P_diff_log*d_mtc[MAX_VARIABLE_TYPES+w] 
			- mass_tran_coeff*P_diff_log_dw[w])/SQUARE(P_diff_log);
		}
        mass_tran_coeff *= amb_pres/P_diff_log;
        }
    }

 /* HARDWIRE a linear increase in MTC from zero to mass_tran_coeff
   along free surface boundary */

      mtc = mass_tran_coeff ;
  
 /* there is sensitivity of flux with respect to concentration
     and temperature */

  if (af->Assemble_Jacobian) {

    for (w=0; w<pd->Num_Species_Eqn; w++) {
      d_mass_flux[wspec][MAX_VARIABLE_TYPES + w] =  mtc*A*dact_dC[wspec][w]
		+ d_mtc[MAX_VARIABLE_TYPES+w]*(A*activity[wspec] - Y_c);
    		}

    d_mass_flux[wspec][TEMPERATURE] = 
		mtc*activity[wspec]*dpsatdt[wspec]/amb_pres
 		+d_mtc[TEMPERATURE]*(A*activity[wspec] - Y_c);
  }

  mass_flux[wspec] = mtc*(A*activity[wspec] - Y_c);

return;

} /* END of routine mass_flux_equil_mtc      */ 

/*****************************************************************************/
/******************************************************************************
*  Function which calculates the mass transfer coefficent according to the 
*  Chilton-Coburn correlation
*
*           Author: R. B. Secor (4/03)
*
******************************************************************************/
 
void 
mtc_chilton_coburn(dbl *mtc,
 		    dbl d_mtc[MAX_VARIABLE_TYPES+MAX_CONC],
 		    int wspec,		/* species number */
 		    double htc,		/* heat transfer coefficient */
 		    double T_gas,	/* gas Temperature */
 		    double P_total,	/* total pressure */
 		    double diff_gas_25)	/* solvent diff. at 25 C.    */
 
/*****************************************************************************/
{
   /* Local variables */
     double mw_gas=28.951;		/* Molecular wt. of gas (g/mole) */
     double pr_gas=0.71;		/* Prandtl Number		 */
     double rho_gas;			/* gas density	(g/cc)		 */
     double cp_gas;			/* heat capacity (cal/g/deg K)   */
     double diff_gas;			/* solvent diffusivity in gas    */
     double visc_gas;			/* gas viscosity (p.)    */
     double T_film,T_liquid;
     double temp1, temp2;
     int w;
     double convF, convF_dc[MAX_CONC],convF_dT;

/*  Assuming temperature is in degrees K  */
 
	T_liquid = MAX(DBL_SMALL,fv->T);
	T_film = 0.5*(T_liquid + T_gas);	
#if 0
if(T_gas < 0 || T_gas > 1000)
	{
	fprintf(stderr,"chilton-coburn %g\n",T_gas);
	T_film = fv->T;	
	}
#endif

 
/* gas density - ideal gas law */
 
	rho_gas = mw_gas*P_total/(82.05*T_film);
 
/* gas heat capacity	*/
 	if(T_film < 200.)	
 		{  cp_gas = 0.23901; }
 	else
 		{
 		cp_gas = 0.00023901*(28958. 
 			+9390.*SQUARE((3012./T_film)/sinh(3012./T_film))
 			+7580.*SQUARE((1484./T_film)/cosh(1484./T_film))
 			)/mw_gas;
 		}
 
/*  solvent diffusivity	*/
 
 	diff_gas = diff_gas_25* pow(T_film/298.15,1.5);
 
/*  gas viscosity*/
 
 	visc_gas = (1.425E-05*pow(T_film, 0.5039))/(1.0 + 108.3/T_film);
 
/*	Factor for converting back to volume fraction basis  */
        memset(convF_dc,0,sizeof(dbl)*MAX_CONC);
        convF = 1.0;   convF_dT = 0.0;
        switch(mp->Species_Var_Type)   {
          case SPECIES_MASS_FRACTION:
          case SPECIES_UNDEFINED_FORM:
              break;
          case SPECIES_CONCENTRATION:
              convF = 0.0;
              for( w=0 ; w<pd->Num_Species_Eqn ; w++)
		{  
                    convF += fv->c[w]*(1.0-mp->molecular_weight[w]/mp->molecular_weight[pd->Num_Species_Eqn]);
                    convF_dc[w] += (1.0-mp->molecular_weight[w]/mp->molecular_weight[pd->Num_Species_Eqn]);
                }
               convF += 1.0/mp->molecular_weight[pd->Num_Species_Eqn];
              break;
          case SPECIES_DENSITY:
          default:
              convF = mp->density;
              for( w=0 ; w<pd->Num_Species_Eqn ; w++)
		{  convF_dc[w] = mp->d_density[MAX_VARIABLE_TYPES+w];}
              convF_dT = mp->d_density[TEMPERATURE];
              break;
           }

/*  mass transfer coefficient	*/
 
/*	mass transfer coefficient based on Pr no. = constant
		- favored method by Pete Price		*/
	temp1 = pow(pr_gas*rho_gas*diff_gas/visc_gas,0.67);
	temp2 = rho_gas*cp_gas*82.05*T_film;
 	*mtc = htc*mp->molecular_weight[wspec]*temp1/temp2;
 	*mtc = *mtc/convF;
 
/*  Due to sensitivity problems only compute explicit temperature derivs  */
 	d_mtc[TEMPERATURE] = -0.5*htc*mp->molecular_weight[wspec]*
		(temp1/temp2/T_film);
 	d_mtc[TEMPERATURE] = d_mtc[TEMPERATURE]/convF;
	d_mtc[TEMPERATURE] += - (*mtc)*convF_dT/convF;
	for( w=0 ; w<pd->Num_Species_Eqn ; w++)
		{
		d_mtc[MAX_VARIABLE_TYPES+w] = -(*mtc)*convF_dc[w]/convF;
		}
 	return;
}  /*  end of mtc_chilton_coburn  */
/******************************************************************************
*  Function which calculates the first and second derivatives of
*  the natural log of activity coefficient.
*  Only implemented for FLORY-HUGGINS activity for polymer-solvent mixtures.
*
*           Author: A.C. Sun (5/00)
*
******************************************************************************/

void 
act_coeff(dbl lngamma[MAX_CONC], dbl dlngamma_dC[MAX_CONC][MAX_CONC],
		    dbl d2lngamma_dC2[MAX_CONC][MAX_CONC][MAX_CONC],
		    double y_mass[MAX_CONC], /* conc at boundary             */
		    int mode,	/* model to which the VLE is based on        */
		    int wspec)	/* species no.                               */
/*****************************************************************************/
{
     int i, j, jac, k, l, mn;
     int Num_S1, Num_S2;
     double flory1[MAX_CONC],flory2[MAX_CONC],flory3[MAX_CONC];
     double df1_dc[MAX_CONC][MAX_CONC],df2_dc[MAX_CONC][MAX_CONC];
     double df3_dc[MAX_CONC][MAX_CONC],df_dc[MAX_CONC][MAX_CONC];
     double d2f1_dc2[MAX_CONC][MAX_CONC][MAX_CONC], 
       d2f2_dc2[MAX_CONC][MAX_CONC][MAX_CONC], 
       d2f3_dc2[MAX_CONC][MAX_CONC][MAX_CONC],
       d2f_dc2[MAX_CONC][MAX_CONC][MAX_CONC];
     double dv_dw[MAX_CONC][MAX_CONC];
     double C[MAX_CONC], vol[MAX_CONC], sv[MAX_CONC];
     double mw[MAX_CONC], prod[MAX_CONC]; 
     double bottom, prod2, sum_C;
     double chi[MAX_CONC][MAX_CONC]; /* chi is the binary interaction parameter*/
     double mw_last=0; /* Molecular weight of non-condensable and conversion factor */

/***************************** EXECUTION BEGINS *******************************/

      memset(lngamma, 0,sizeof(double)*MAX_CONC);
      memset(dlngamma_dC, 0,sizeof(double)*MAX_CONC*MAX_CONC);
      memset(C, 0, sizeof(double)*MAX_CONC);
      memset(prod, 0, sizeof(double)*MAX_CONC);

  if(mode==RAOULT)
    {
      bottom=0.;
      prod2=1.;

      for (i = 0; i<pd->Num_Species_Eqn; i++) 
	{
	  if(mp->molecular_weight[i] < 0.)
	    {
	      EH(-1, "Molecular weight of the species not specified in the material file");
	    }
	  else
	    {
	      mw[i]  = mp->molecular_weight[i];
	    }
	}

      if(mp->molecular_weight[pd->Num_Species_Eqn] < 0.)
	{
	  EH(-1, "M.W. of a non-condensable species not specified in the material file");
	}
      else
	{
	  mw_last = mp -> molecular_weight[pd->Num_Species_Eqn];
	}
      
      for(i=0;i<pd->Num_Species_Eqn;i++)
	 {
	   prod[i]=1.;
	   for(j=0;j<pd->Num_Species_Eqn;j++)
	     {
	       prod[i] *= (i==j? 1.0 :mw[j]);
	     } 
	   prod2 *= mw[i];
	 }

      for(i=0;i<pd->Num_Species_Eqn ;i++)
	 {
	   bottom += (mw_last*prod[i]-prod2)*C[i];
	 }
      bottom += prod2;
   
      lngamma[wspec] = mw_last*prod[wspec]/bottom;

      /* write out the derivatives of mol_frac w.r.t. mass_frac*/
      
      for(j=0;j<pd->Num_Species_Eqn;j++)
	{
	  dlngamma_dC[wspec][j] = -mw_last*prod[wspec]
		 *(mw_last*prod[j]-prod2)/(bottom*bottom) ;
	}
    }

  else if(mode==FLORY || mode==FLORY_CC)

    {
  /* Define some convenient/repetitive chunks to make eqns more compact */

      Num_S1 = pd->Num_Species_Eqn + 1;
      Num_S2 = pd->Num_Species_Eqn - 1;
      
      memset(flory1, 0,sizeof(double)*MAX_CONC);
      memset(flory2, 0,sizeof(double)*MAX_CONC);
      memset(flory3, 0,sizeof(double)*MAX_CONC);
      memset(df_dc, 0,sizeof(double)*MAX_CONC*MAX_CONC);
      memset(df1_dc, 0,sizeof(double)*MAX_CONC*MAX_CONC);
      memset(df2_dc, 0,sizeof(double)*MAX_CONC*MAX_CONC);
      memset(df3_dc, 0,sizeof(double)*MAX_CONC*MAX_CONC);
      memset(dv_dw, 0,sizeof(double)*MAX_CONC*MAX_CONC);
      memset(d2f1_dc2, 0,sizeof(double)*MAX_CONC*MAX_CONC*MAX_CONC);
      memset(d2f2_dc2, 0,sizeof(double)*MAX_CONC*MAX_CONC*MAX_CONC);
      memset(d2f3_dc2, 0,sizeof(double)*MAX_CONC*MAX_CONC*MAX_CONC);
      memset(d2f_dc2, 0,sizeof(double)*MAX_CONC*MAX_CONC*MAX_CONC);
  
/* Seems this error checking should be done at input file reading time...*/
       for (i = 0; i<Num_S1; i++) 
	 {
	   if(mp->specific_volume[i] < 0.)
	     { EH(-1, "Specific volume not specified in the material file."); }
	   else
	     { sv[i] = mp->specific_volume[i]; }

	   if(mp->molar_volume[i] < 0.)
	     { EH(-1, "Molar volume not specified in the material file"); }
	   else
	     { vol[i] = mp->molar_volume[i]; }

	   for(k=0;k<Num_S1;k++)
	     {
	     if(mp->flory_param[i][k] < 0.)
	       { EH(-1, "Flory-Huggins binary parameters not specified in the material file"); }
	     else
	       { chi[i][k] = mp->flory_param[i][k]; }
	     }
	 }

       bottom = 0.;
       sum_C = 0.;
/*  denominators for fraction-based formulations */
	if(mp->Species_Var_Type == SPECIES_MASS_FRACTION
           || mp->Species_Var_Type == SPECIES_UNDEFINED_FORM)
	{
         for (i = 0; i<pd->Num_Species_Eqn; i++) 
	    {	   
	     bottom += y_mass[i]*(sv[i]-sv[pd->Num_Species_Eqn]);
	    }
         bottom += sv[pd->Num_Species_Eqn];
	}
	if(mp->Species_Var_Type == SPECIES_MOLE_FRACTION)
	{
         for (i = 0; i<pd->Num_Species_Eqn; i++) 
	    {	   
	     bottom += y_mass[i]*(vol[i]-vol[pd->Num_Species_Eqn]);
	    }
         bottom += vol[pd->Num_Species_Eqn];
	}
	   
	if(mp->Species_Var_Type == SPECIES_DENSITY)
		{
                 for(i=0;i<pd->Num_Species_Eqn;i++)
	           {
	            C[i] = y_mass[i]*sv[i];
	            sum_C += C[i];
		    dv_dw[i][i]= sv[i];
	           }
		}
	 else if(mp->Species_Var_Type == SPECIES_CONCENTRATION)
		{
                 for(i=0;i<pd->Num_Species_Eqn;i++)
	           {
	            C[i] = y_mass[i]*vol[i];
	            sum_C += C[i];
	            dv_dw[i][i]= vol[i];
	           }
		}
	else if(mp->Species_Var_Type == SPECIES_MASS_FRACTION
           || mp->Species_Var_Type == SPECIES_UNDEFINED_FORM)
		{
                 for(i=0;i<pd->Num_Species_Eqn;i++)
	           {
/*	            C[i] = y_mass[i]*sv[i]/bottom;  */
	            C[i] = y_mass[i]*sv[i];
	            sum_C += C[i];
	            for(j=0;j<pd->Num_Species_Eqn;j++)
	               {
                   /*     dv_dw[i][j] = -y_mass[i]*sv[i]
		              *(sv[j]-sv[pd->Num_Species_Eqn])/SQUARE(bottom); 
*/
                       }
/*		    dv_dw[i][i] += sv[i]/bottom;    */
		    dv_dw[i][i] += sv[i]; 
                   }
		}
	else if(mp->Species_Var_Type == SPECIES_MOLE_FRACTION)
		{
                 for(i=0;i<pd->Num_Species_Eqn;i++)
	           {
	            C[i] = y_mass[i]*vol[i]/bottom;
	            sum_C += C[i];
	            for(j=0;j<pd->Num_Species_Eqn;j++)
	               {
                        dv_dw[i][j] = -y_mass[i]*vol[i]
		              *(vol[j]-vol[pd->Num_Species_Eqn])/SQUARE(bottom); 
                       }
		    dv_dw[i][i] += vol[i]/bottom; 
		   }
		}
	else
		{ EH(-1,"That species formulation not done in mtc_flory\n"); }
       for(k=0;k<pd->Num_Species_Eqn; k++)
	 {
	   flory1[wspec] += vol[wspec]*(delta(k,wspec)-C[k])/vol[k]
	             +vol[wspec]*C[k]/vol[pd->Num_Species_Eqn];

	   df1_dc[wspec][k]= (vol[wspec]/vol[pd->Num_Species_Eqn])
	              - (vol[wspec]/vol[k]) ;
	 }
       flory1[wspec] += -vol[wspec]/vol[pd->Num_Species_Eqn];

       if (pd->Num_Species_Eqn > 1)
	 {
	   for(k=1;k<pd->Num_Species_Eqn; k++)
	     {
	       for(jac=0;jac<k;jac++)
		 {   
		   flory2[wspec] += (delta(wspec,jac)*C[k]+(vol[wspec]/vol[jac])
			       *C[jac]*(delta(k,wspec)-C[k]))*chi[jac][k];
		 }
	       /* derivative for the 1st component */
	       df2_dc[wspec][0] += (vol[wspec]*(delta(k,wspec)-C[k])/vol[0])
		 *chi[0][k];
	       d2f2_dc2[wspec][0][k] += (-vol[wspec]/vol[0])*chi[0][k];
	     }

	   /* derivative for the last component */
	   for(jac=0;jac<Num_S2;jac++)
	     {
	       df2_dc[wspec][Num_S2] += (delta(wspec,jac)-vol[wspec]*C[jac]
				  /vol[jac])*chi[jac][Num_S2];
	       d2f2_dc2[wspec][Num_S2][jac] += (-vol[wspec]/vol[jac])*chi[jac][Num_S2];
	     }
	   /* derivative for the components in between */
	   for(l=1;l<Num_S2 ;l++)
	     {
	       for(jac=0;jac<l;jac++)
		 {
		  df2_dc[wspec][l] += (delta(wspec,jac)-(vol[wspec]*C[jac]
				  /vol[jac]))*chi[jac][l];
		  d2f2_dc2[wspec][l][jac] += (-vol[wspec]/vol[jac])*chi[jac][l];
		 }
	       for(k=l+1;k<pd->Num_Species_Eqn;k++)
		 {
		  df2_dc[wspec][l] +=(vol[wspec]*(delta(k,wspec)-C[k])/vol[l])
		                 *chi[l][k];
		  d2f2_dc2[wspec][l][k] += (-vol[wspec]/vol[l])*chi[l][k];

		 }
	     }
	 }
       else
	 {
	   flory2[wspec] = 0.;
           for(k=0;k<pd->Num_Species_Eqn;k++)  df2_dc[wspec][k] = 0.;
	 }

       for(jac=0;jac<pd->Num_Species_Eqn;jac++)
	 {
	   flory3[wspec] +=(delta(wspec,jac)-(vol[wspec]/vol[jac])*C[jac])
	           *(1.-sum_C)*chi[jac][pd->Num_Species_Eqn];
           df3_dc[wspec][jac] -= (vol[wspec]/vol[jac])*(1.-sum_C)
                            *chi[jac][pd->Num_Species_Eqn];

	   for(mn=0;mn<pd->Num_Species_Eqn; mn++)
	     {
	       df3_dc[wspec][jac]
		 -= (delta(wspec,mn)-(vol[wspec]/vol[mn])*C[mn])
		       *chi[mn][pd->Num_Species_Eqn];

	       d2f3_dc2[wspec][jac][mn] += 
		 (vol[wspec]/vol[jac])*chi[jac][pd->Num_Species_Eqn]
		 + (vol[wspec]/vol[mn])*chi[mn][pd->Num_Species_Eqn];
	     }
	 }
       lngamma[wspec]= flory1[wspec]+flory2[wspec]+flory3[wspec];

       /* check the simplest case: 1solvent, 1polymer 
       check=(1-sum_C)+chi[0][1]*(1.-sum_C)*(1.-sum_C);
       printf("flory = %e, check =%e\n", flory[i],check); */

       for(k=0;k<pd->Num_Species_Eqn;k++)
         {
           df_dc[wspec][k]= df1_dc[wspec][k]+df2_dc[wspec][k]+df3_dc[wspec][k];

	   for(l=0;l<pd->Num_Species_Eqn;l++)
	     {
	      d2f_dc2[wspec][k][l] = 
	      d2f1_dc2[wspec][k][l]+d2f2_dc2[wspec][k][l]+d2f3_dc2[wspec][k][l];
	     }
	 }
      
       for(i=0;i<pd->Num_Species_Eqn;i++)
	 {
	   for(j=0;j<pd->Num_Species_Eqn;j++)
	     {
	       dlngamma_dC[wspec][i] += df_dc[wspec][j]*dv_dw[j][i];

	       for(k=0;k<pd->Num_Species_Eqn;k++)
		 {
		   for(l=0;l<pd->Num_Species_Eqn;l++)
		     {
		       d2lngamma_dC2[wspec][i][j] += d2f_dc2[wspec][k][l]
			 *dv_dw[l][i]*dv_dw[k][j];
		     }
		 }
	     }
	 }
    }

return;

} /* END of routine act_coeff */ 

/*****************************************************************************/
/*****************************************************************************
 *  Function which calculates the surface integral for convective
 *  mass transfer.
 *
 *  ----------------------------------------------------------------------------
 *
 * To account for VLE at the external boundary using a YFLUX_EQUIL, 
 * the concentration at the surface is obtained from applying Raoult or
 * Flory-Huggins VLE model.
 *          Author: A.C. Sun 9/98                                            */
/*****************************************************************************/

void 
get_equil_surf_bc(double func[],
		  double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
		  int mode,	/* model on which the VLE is based           */
		  int wspec,	/* species number of this boundary condition */
		  double amb_pres,
		  double mass_tran_coeff,
		  double Y_c,	/* bath concentration 	                     */
 		  double T_gas,
 		  double diff_gas_25,	/* Chilton-Coburn parameters	     */
		  dbl dt,	/* current value of the time step            */
		  dbl tt)	/* parameter to vary time integration        */

/********************************************************************************
 *  Functions called:
 *  mass_flux_equil_mtc  -- calculates mass flux for one component
 *
 *******************************************************************************/
{
  int j, j_id, w1, dim, kdir, var, jvar;
  double phi_j;
  double Y_w; /* local concentration of current species */
  
  double activity[MAX_CONC]; /* nonideal activity of species */
  double dact_dC[MAX_CONC][MAX_CONC];
  double vconv[MAX_PDIM]; /*Calculated convection velocity */
  double vconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;
  double mtc, d_mtc[MAX_VARIABLE_TYPES+MAX_CONC];
  int err=0;  

  /***************************** EXECUTION BEGINS ****************************/
  
  if (af->Assemble_LSA_Mass_Matrix) return;

  /* 
   *  call routine to calculate surface flux of this component and it's 
   *  sensitivity to all variable types 
   */
 		/*  Chilton-Coburn correlation if FLORY_CC */
  memset(d_mtc, 0, sizeof(dbl)*(MAX_VARIABLE_TYPES+MAX_CONC))    ;
		  if( mode == FLORY_CC)
 			{
 			mtc_chilton_coburn(&mtc, d_mtc, wspec, 
 				mass_tran_coeff,
 				T_gas,
 				amb_pres,
 				diff_gas_25);
 			}
 		  	else
 			{
 			mtc=mass_tran_coeff;
 			}
  mass_flux_equil_mtc (mp->mass_flux,mp->d_mass_flux, activity, dact_dC,
                      fv->c, mode, amb_pres, wspec, mtc, d_mtc, Y_c);

  dim   = pd->Num_Dim;
  Y_w = fv->c[wspec];

  /* get the convection velocity (it's different for arbitrary and
     lagrangian meshes) */

  /* get deformation gradient, etc. if lagrangian mesh with inertia */
  if ((pd->MeshMotion == LAGRANGIAN ||
       pd->MeshMotion == DYNAMIC_LAGRANGIAN) && pd->MeshInertia == 1)
    {
      err = belly_flop(elc->lame_mu);
      EH(err, "error in belly flop");
      if( neg_elem_volume ) return;
    }

/*  if (mp->PorousMediaType == CONTINUOUS) */
  err = get_convection_velocity(vconv, vconv_old, d_vconv, dt, tt);
  EH(err, "Error in calculating effective convection velocity");


  /*
   *   Calculate the residual contribution for this boundary condition
   * 
   *     flux = int ( Y_k * (v dot snormal) - mass_flux)
   *
   *  HKM Note:
   *   The density should multiply the first term for the units to work out.
   *   The interfacial velocity is subtracted out from the mass-averaged
   *   velocity in the routine get_convection_velocity() first, so it has
   *   not been neglected. 
   */
  if (af->Assemble_Residual) {
    *func -= mp->mass_flux[wspec];
    for (kdir = 0; kdir < dim; kdir++) {
      *func += Y_w * vconv[kdir] * fv->snormal[kdir];
    }
  }
  
  if (af->Assemble_Jacobian ) 
    {
      /* sum the contributions to the global stiffness matrix  for Species*/
      
      /*
       * J_s_c
       */
      var = MASS_FRACTION;
      for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
	phi_j = bf[var]->phi[j_id];
	
	for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++ )
	  {
	    d_func[0][MAX_VARIABLE_TYPES + w1][j_id] -= 
	      mp->d_mass_flux[wspec][MAX_VARIABLE_TYPES + w1] 
		* phi_j;
	  }
	
	for (kdir = 0; kdir < dim; kdir++) 
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
		  for (kdir = 0; kdir < dim; kdir++) 
		    {
		      d_func[0][var][j_id] += 
			Y_w*vconv[kdir]*fv->dsnormal_dx[kdir][jvar][j_id] +
			Y_w*d_vconv->X[kdir][jvar][j_id]*fv->snormal[kdir];
		    }
		}
	    }
	}
      
      for (jvar=0; jvar<dim; jvar++) 
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
  
} /* END of routine get_equil_surf_bc       */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/****************************************************************************
 *
 *  Function which calculates the surface integral for convective
 *  mass transfer.
 *
 *  -------------------------------------------------------------------------
 *
 *  Functions called:
 *  mass_flux_surf_mtc  -- calculates mass flux for one component
 *
 *  -------------------------------------------------------------------------
 *
 ****************************************************************************/

void 
sus_mass_flux_surf_bc (double func[],
		       double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
		       const int sus_species, /* species number of this 
					       * boundary condition */
		       const dbl time, /* current time */
		       const dbl dt, /* current value of the time step */
		       const dbl tt, /* parameter to vary time integration 
				      * from BE(0) to CN(1/2) to FE(1) */
		       const dbl hsquared[DIM])
{
  
  /* Local variables */
  
  int j, j_id,w1,dim, a;
  int var, jvar;
  
  dbl rho;
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;  /* density dependence */
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;

  struct Species_Conservation_Terms s_terms; 
  
  /***************************** EXECUTION BEGINS *******************************/
  
  if(af->Assemble_LSA_Mass_Matrix)
    return;

  dim   = pd->Num_Dim;

  if( mp->DensityModel == SUSPENSION  )
    { 
      if(sus_species != (int) mp->u_density[0])
	{
	  EH( -1, "YFLUX_SUS_BC and density model species # must be consistent");
	}

      /*** Density ***/
      rho = density(d_rho, time);
      /* MMH
       * "Fix" the density call if we are in the SUSPENSION_PM model, and
       * we are actually in the particle phase... 
       */
      if(mp->DensityModel == SUSPENSION_PM &&
	 sus_species == mp->u_density[0])
	rho = mp->u_density[2];

      hydro_flux( &s_terms, sus_species, tt, dt, hsquared); 
    }
  else
    {
      return;
    }

  if (af->Assemble_Residual ) 
    {
      /* Calculate the residual contribution from convective flux	*/
      for(a=0; a<dim; a++) 
	{
	  *func += fv->snormal[a] * rho
	    * s_terms.diff_flux[sus_species][a]; 	
	}
    }
  
  if (af->Assemble_Jacobian ) 
    {
      /* sum the contributions to the global stiffness matrix  for Species*/
      
      /*
       * J_s_c
       */
      var=MASS_FRACTION;
      for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) 
	{
	  for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++ )
	    {
	      for(a=0; a<dim; a++) 
		{
		  d_func[0][MAX_VARIABLE_TYPES + w1][j_id] += fv->snormal[a] *
		    s_terms.d_diff_flux_dc[sus_species][a][w1][j_id];
		}
	    }
	}
      
      /*
       * J_s_T
       */
      var=TEMPERATURE;
      if (pd->v[pg->imtrx][var])
	{
	  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) 
	    {
	      j_id = j;
	      for(a=0; a<dim; a++) 
		{
		  d_func[0][var][j_id]  += fv->snormal[a] 
		    * ( rho * s_terms.d_diff_flux_dT[sus_species][a][j_id] +
			d_rho->T[j] * s_terms.diff_flux[sus_species][a] );
		}
	    }
	}

     /*
       * J_s_SH
       */
      var= SHEAR_RATE;
      if (pd->v[pg->imtrx][var])
	{
	  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) 
	    {
	      j_id = j;
	      for(a=0; a<dim; a++) 
		{
		  d_func[0][var][j_id]  += fv->snormal[a]
		    * (s_terms.d_diff_flux_dSH[sus_species][a][j_id]);
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
		  for(a=0; a<dim; a++) 
		    {
		      d_func[0][var][j_id] += fv->dsnormal_dx[a][jvar][j_id] * rho
			* s_terms.diff_flux[sus_species][a]
			+ fv->snormal[a] * rho 
			* s_terms.d_diff_flux_dmesh[sus_species][a][jvar][j_id]; 
		    }
		}
	    }
	}
      
    } /* End of if Assemble_Jacobian */
  
} /* END of routine sus_mass_flux_surf_bc  */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void 
compute_leak_velocity(double *vnorm,
	              NORMAL_VELOCITY_DEPENDENCE_STRUCT *d_vnorm,
	              dbl tt,		/* parameter to vary time integration from 
		        		 * explicit (tt = 1) to implicit (tt = 0)    */
	              dbl dt,		/* current value of the time step            */
	              struct Boundary_Condition *bc,
		      struct Boundary_Condition *fluxbc)

/******************************************************************************
*
*  Function which evaluates the normal leak velocity used in the
*  kinematic boundary condition n.(v -vs)=0
*/
{
  int i, j, p, q;
  int var;
  int w, wspec, w1, w2;
  int num_comp;
  int add_fluxes = FALSE;
  double StoiCoef[MAX_CONC];
  double mass_tran_coeff,Y_c, density_tot;
  double nu, k, beta, alphaa, alphac, V = 0.0, U0, T = 0.0, nd;
  double M_solid, rho_solid, molar_volume;
  double ai0, H, cref, n;
  double k1, E1, kn1, En1, c_H2S, c_O2;
  double xbulk, d_xbulk_dC[MAX_CONC];
  double vnormal, phi_j;

  int mode;
  double amb_pres, A, mtc, Y_inf, driving_force;
  double d_mtc[MAX_VARIABLE_TYPES+MAX_CONC];
  double activity[MAX_CONC];
  double dact_dC[MAX_CONC][MAX_CONC];
  double psat[MAX_CONC], dpsatdt[MAX_CONC];
  PROPERTYJAC_STRUCT *densityJac = NULL;
  propertyJac_realloc(&densityJac, mp->Num_Species+1);


  /* local contributions of boundary condition to residual and jacobian */

/***************************** EXECUTION BEGINS *******************************/

  /* Need to calculate vnorm and it's derivatives here.
     this depends on whether YFLUX_BC or YFLUX_BV_BC or YFLUX_EQUIL_BC is chosen */

  vnormal = 0.;
  
  memset(d_vnorm->v, 0, sizeof(dbl)*DIM*MDE);
  memset(d_vnorm->T, 0, sizeof(dbl)*MDE);
  memset(d_vnorm->C, 0, sizeof(dbl)*MAX_CONC*MDE);
  memset(d_vnorm->V, 0, sizeof(dbl)*MDE);
  memset(d_vnorm->F, 0, sizeof(dbl)*MDE);
  memset(d_vnorm->X, 0, sizeof(dbl)*DIM*MDE);

  

  for (i=0; i < pd->Num_Species_Eqn; i++) 
      { StoiCoef[i] = 1.0; }
  /* Expand function use to include KIN_CHEM BC */
  if ( bc != NULL && bc->BC_Name == KIN_CHEM_BC ) 
    {
    num_comp = bc->len_u_BC;
    for (i=0; i < num_comp; i++) 
      { StoiCoef[i] = bc->u_BC[i]; }
    }

/* convert mass flux to volume flux depending on Species Formulation */
 density_tot = calc_density(mp, TRUE, densityJac, 0.0);
 switch(mp->Species_Var_Type)   {
    case SPECIES_UNDEFINED_FORM:
    case SPECIES_MASS_FRACTION:
         for (i=0; i < MAX_CONC; i++) 
              { StoiCoef[i] *= density_tot*mp->specific_volume[i]; }
/* Probably need a d_StoiCoef_dC[][] for FRACTION formulations  */
          break;
    case SPECIES_MOLE_FRACTION:
    case SPECIES_VOL_FRACTION:
          EH(-1, "Volume conversion not done for that Species Formulation");
          break;
    case SPECIES_CONCENTRATION:
         for (i=0; i < MAX_CONC; i++) 
              { StoiCoef[i] *= mp->molar_volume[i]; }
         break;   
    case SPECIES_DENSITY:
    default:
         for (i=0; i < MAX_CONC; i++) 
              { StoiCoef[i] *= mp->specific_volume[i]; }
         break;
    }

  /* call routine to calculate surface flux of this component and it's 
   * sensitivity to all variable types 
   */
  
  /* Calculate flux contribution of bulk component */
  xbulk = 1.;
  if (pd->v[pg->imtrx][MASS_FRACTION])
    {
        switch(mp->Species_Var_Type)   {
        case SPECIES_UNDEFINED_FORM:
        case SPECIES_MASS_FRACTION:
           for (w=0; w<pd->Num_Species_Eqn; w++) 
              { xbulk -= fv->c[w]; d_xbulk_dC[w] = -1.0; }
           break;
        case SPECIES_MOLE_FRACTION:
        case SPECIES_VOL_FRACTION:
           EH(-1, "BC mass fraction conversion not done for that Species Formulation");
           break;
        case SPECIES_CONCENTRATION:
           for (w=0; w<pd->Num_Species_Eqn; w++) 
              {
               xbulk -= fv->c[w]*mp->molecular_weight[w]/density_tot;
               d_xbulk_dC[w] = -mp->molecular_weight[w]*
        (1.0/density_tot+fv->c[w]*mp->d_density[MAX_VARIABLE_TYPES+w]/SQUARE(density_tot));
              }
           break;   
        case SPECIES_DENSITY:
        default:
           for (w=0; w<pd->Num_Species_Eqn; w++) 
              {
               xbulk -= fv->c[w]/density_tot;
               d_xbulk_dC[w] = -(1.0/density_tot
                      +fv->c[w]*mp->d_density[MAX_VARIABLE_TYPES+w]/SQUARE(density_tot));
              }
           break;
         }
    }
  
  /* Calculate volume flux of bulk component through surface */
  mass_tran_coeff = bc->BC_Data_Float[0];
  Y_c             = bc->BC_Data_Float[1];
  
  if (af->Assemble_Jacobian) {
    var = MASS_FRACTION;
    if (pd->v[pg->imtrx][var])
      {
	for (w=0; w<pd->Num_Species_Eqn; w++) {
	  for (j=0; j<ei[pg->imtrx]->dof[var]; j++) {
	    phi_j = bf[var]->phi[j];
	    d_vnorm->C[w][j] = -mp->specific_volume[pd->Num_Species_Eqn]*
                  mass_tran_coeff*(density_tot * d_xbulk_dC[w]
                    +mp->d_density[MAX_VARIABLE_TYPES+w]*xbulk)*phi_j;
	  }
	}
      }
  }
  vnormal += density_tot*mp->specific_volume[pd->Num_Species_Eqn]*
                  mass_tran_coeff * ( xbulk - Y_c );

  add_fluxes = (fluxbc == NULL);
  while ( fluxbc != NULL )
    {
      
      if(!strcmp(fluxbc->desc->name1,"YFLUX_EQUIL"))  /* this is limited to multicomponent */
	{
	  
	  /* call routine to calculate surface flux of bulk component and it's 
	   * sensitivity to all variable types    
           * ACS: modified 10/99 to accommodate mass conc. formulation for YFLUX_EQUIL  */
	  
	  driving_force = 1.;	  
	  
	  if (pd->v[pg->imtrx][MASS_FRACTION])
	    {
	      wspec    = fluxbc->BC_Data_Int[0];
	      mode     = fluxbc->BC_Data_Int[2];
	      amb_pres = fluxbc->BC_Data_Float[0];
	      mtc      = fluxbc->BC_Data_Float[1];
	      Y_inf    = fluxbc->BC_Data_Float[2];

              /*  Chilton-Coburn correlation if FLORY_CC */
              memset(d_mtc, 0, sizeof(dbl)*(MAX_VARIABLE_TYPES+MAX_CONC));
	      if( mode == FLORY_CC)
 		{
		  mtc_chilton_coburn(&mtc, d_mtc, wspec, 
 				fluxbc->BC_Data_Float[1],
 				fluxbc->BC_Data_Float[3],
 				fluxbc->BC_Data_Float[0],
 				fluxbc->BC_Data_Float[4]);
 		}
		  
	      /* Nonideal VP Calculations based on either ANTOINE or RIEDEL models */
	      
	      if(mp->VaporPressureModel[wspec] == ANTOINE )
		{
		  antoine_psat(wspec, mp->u_vapor_pressure[wspec],
			       &psat[wspec], &dpsatdt[wspec]);
		  mp-> vapor_pressure[wspec] = psat[wspec];
		}
	      
	      else if(mp->VaporPressureModel[wspec] == RIEDEL )
		{
		  riedel_psat(wspec, mp->u_vapor_pressure[wspec],
			      &psat[wspec], &dpsatdt[wspec]);
		  mp-> vapor_pressure[wspec] = psat[wspec];
		}
		  
	      /* Calculate mole flux of other components through surface */
	      mtc *= StoiCoef[wspec];
	      d_mtc[TEMPERATURE] *= StoiCoef[wspec];
    	      for (w=0; w<pd->Num_Species_Eqn; w++) 
	          d_mtc[MAX_VARIABLE_TYPES+w] *= StoiCoef[wspec];
	      mass_flux_equil_mtc (mp->mass_flux,mp->d_mass_flux, activity, dact_dC, 
				   fv->c, mode, amb_pres, wspec, mtc, d_mtc, Y_inf);
		  
	      A = mp->vapor_pressure[wspec]/amb_pres;
	      driving_force -=  A*activity[wspec];
	      vnormal += mp->mass_flux[wspec] ;
/* This was causing compiler warnings: probably due to C[w] which should have been C[w][j]  */
    	      var = MASS_FRACTION;
    	      if (pd->v[pg->imtrx][var])
      		{
	      	for (w=0; w<pd->Num_Species_Eqn; w++) 
		   {
	  	    for (j=0; j<ei[pg->imtrx]->dof[var]; j++) 
			{
	    		phi_j = bf[var]->phi[j];
			d_vnorm->C[w][j] += 
				mp->d_mass_flux[wspec][MAX_VARIABLE_TYPES+w]*phi_j;
			}
	      	   }
		}

	      var = TEMPERATURE;
	      if (pd->v[pg->imtrx][var])
	        {
                  for (j=0; j<ei[pg->imtrx]->dof[var]; j++) 
		  {
                    phi_j = bf[var]->phi[j];
                    d_vnorm->T[j] += mp->d_mass_flux[wspec][TEMPERATURE] * phi_j;
                  }
                }
	    }
	  	  
	}
      else if (!strcmp(fluxbc->desc->name1,"YFLUX_CONST"))  
	{
	  
	  if (pd->v[pg->imtrx][MASS_FRACTION])
	    {
	      wspec    = fluxbc->BC_Data_Int[0];
	      vnormal -= fluxbc->BC_Data_Float[0] * StoiCoef[wspec];
	    }
	}
      else if (!strcmp(fluxbc->desc->name1,"YFLUX_BV"))  
	{                            /* mass flux given by Butler-Volmer kinectics, KSC: 11/2000 */
	  
	  if (pd->v[pg->imtrx][MASS_FRACTION])
	    {
	      wspec         = fluxbc->BC_Data_Int[0];   /* species number */   
	      nu            = fluxbc->BC_Data_Float[0]; /* stoichiometric coeff. */
	      k             = fluxbc->BC_Data_Float[1]; /* rate constant  */
	      beta          = fluxbc->BC_Data_Float[2]; /* reaction order */
	      alphaa        = fluxbc->BC_Data_Float[3]; /* anodic transfer coeff. */
	      alphac        = fluxbc->BC_Data_Float[4]; /* cathodic transfer coeff. */ 
	      V             = fluxbc->BC_Data_Float[5]; /* electrode potential */
	      U0            = fluxbc->BC_Data_Float[6]; /* electrolyte open-circuit potential */
	      M_solid       = fluxbc->BC_Data_Float[7]; /* molecular weight of solid deposit  */ 
	      rho_solid     = fluxbc->BC_Data_Float[8]; /* density of solid deposit */ 
	      molar_volume  = M_solid/rho_solid;              /* molar volume of solid deposit */ 
	      if (mp->SolutionTemperatureModel == CONSTANT)  
		{
		  T = mp->solution_temperature;
		} 
	      else 
		{
		  EH(-1, "Solution-temperature model not yet implemented");
		}
	      if(pd->e[pg->imtrx][R_ENERGY]) /* if energy equation is active, re-set electrolyte temperature */
		{
		  T = fv->T;   
		}
	      /* set solution temperature to 298 K if it is zero - safety feature */
	      if (T == 0.0) 
		{
		  T = 298.0;  
		  fprintf(stderr, "Warning!: a default electrolyte temperature of 298 K is being used!");
		}

	      mass_flux_surf_BV (mp->mass_flux, mp->d_mass_flux, wspec,  
				 nu, k, beta, alphaa, alphac, V, U0, T);
	      vnormal += molar_volume*mp->mass_flux[wspec] * StoiCoef[wspec];   

	      var = MASS_FRACTION;
              for (w = 0; w < pd->Num_Species_Eqn; w++) 
                {
		  for (j=0; j<ei[pg->imtrx]->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];
                      d_vnorm->C[w][j] += molar_volume*mp->d_mass_flux[wspec][MAX_VARIABLE_TYPES +w] * StoiCoef[wspec] * phi_j; 
                    }
                }
              var = TEMPERATURE;
	      if (pd->v[pg->imtrx][var]) 
	        {
                  for (j=0; j<ei[pg->imtrx]->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];
                      d_vnorm->T[j] += molar_volume*mp->d_mass_flux[wspec][TEMPERATURE] * StoiCoef[wspec] * phi_j; 
                    }
                }
              var = VOLTAGE;
	      if (pd->v[pg->imtrx][var]) 
	        {
                  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];
                      d_vnorm->V[j] += molar_volume*mp->d_mass_flux[wspec][VOLTAGE] * StoiCoef[wspec] * phi_j; 
                    }
                }
	    }
	}  /* end of the else on the if(YFLUX_BV) */

      else if (!strcmp(fluxbc->desc->name1,"YFLUX_HOR"))
            {             /* mass flux given by linearized kinetic model */
             if (pd->v[pg->imtrx][MASS_FRACTION])
              {
                wspec     = fluxbc->BC_Data_Int[0];
                ai0       = fluxbc->BC_Data_Float[0];
                H         = fluxbc->BC_Data_Float[1];
                cref      = fluxbc->BC_Data_Float[2];
                alphaa    = fluxbc->BC_Data_Float[3];
                alphac    = fluxbc->BC_Data_Float[4];
                T         = fluxbc->BC_Data_Float[5];
                U0        = fluxbc->BC_Data_Float[6];
                beta      = fluxbc->BC_Data_Float[7];
                n         = fluxbc->BC_Data_Float[8];
                V         = fluxbc->BC_Data_Float[9];

                mass_flux_surf_HOR (mp->mass_flux, mp->d_mass_flux, wspec,
                                    ai0, H, cref, alphaa, alphac, T, U0,
                                    beta, n, V);
                vnormal += mp->mass_flux[wspec] * StoiCoef[wspec];

                var = MASS_FRACTION;
                for (w = 0; w < pd->Num_Species_Eqn; w++) 
                 {
                  for (j=0; j<ei[pg->imtrx]->dof[var]; j++)
                   {
                     phi_j = bf[var]->phi[j];
                     d_vnorm->C[w][j] += mp->d_mass_flux[wspec][MAX_VARIABLE_TYPES +w] * StoiCoef[wspec] * phi_j;
                   }
                 }
                var = TEMPERATURE;
                if (pd->v[pg->imtrx][var]) 
                 {
                  for (j=0; j<ei[pg->imtrx]->dof[var]; j++)
                   {
                     phi_j = bf[var]->phi[j];
                     d_vnorm->T[j] += mp->d_mass_flux[wspec][TEMPERATURE] * StoiCoef[wspec] * phi_j;
                   }
                 }
                var = VOLTAGE;
                if (pd->v[pg->imtrx][var])
                 {
                   for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];
                      d_vnorm->V[j] += mp->d_mass_flux[wspec][VOLTAGE] * StoiCoef[wspec] * phi_j;
                    }
                 }
              }
           }  /* end of the else on the if(YFLUX_HOR) */

      else if (!strcmp(fluxbc->desc->name1,"YFLUX_ORR"))
            {  /* mass flux given by the Tafel kinetic model for ORR */
             if (pd->v[pg->imtrx][MASS_FRACTION])
              {
                wspec     = fluxbc->BC_Data_Int[0];
                ai0       = fluxbc->BC_Data_Float[0];
                H         = fluxbc->BC_Data_Float[1];
                cref      = fluxbc->BC_Data_Float[2];
                alphac    = fluxbc->BC_Data_Float[3];
                T         = fluxbc->BC_Data_Float[4];
                V         = fluxbc->BC_Data_Float[5];
                U0        = fluxbc->BC_Data_Float[6];
                beta      = fluxbc->BC_Data_Float[7];
                n         = fluxbc->BC_Data_Float[8];

                mass_flux_surf_ORR (mp->mass_flux, mp->d_mass_flux, wspec,
                                    ai0, H, cref, alphac, T, V, U0, beta, n);
                vnormal += mp->mass_flux[wspec] * StoiCoef[wspec];

                var = MASS_FRACTION;
                for (w = 0; w < pd->Num_Species_Eqn; w++) 
                 {
                  for (j=0; j<ei[pg->imtrx]->dof[var]; j++)
                   {
                     phi_j = bf[var]->phi[j];
                     d_vnorm->C[w][j] += mp->d_mass_flux[wspec][MAX_VARIABLE_TYPES +w] * StoiCoef[wspec] * phi_j;
                   }
                 }
                var = TEMPERATURE;
                if (pd->v[pg->imtrx][var]) 
                 {
                  for (j=0; j<ei[pg->imtrx]->dof[var]; j++)
                   {
                     phi_j = bf[var]->phi[j];
                     d_vnorm->T[j] += mp->d_mass_flux[wspec][TEMPERATURE] * StoiCoef[wspec] * phi_j;
                   }
                 }
                var = VOLTAGE;
                if (pd->v[pg->imtrx][var])
                 {
                   for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];
                      d_vnorm->V[j] += mp->d_mass_flux[wspec][VOLTAGE] * StoiCoef[wspec] * phi_j;
                    }
                 }
              }
           }  /* end of the else on the if(YFLUX_ORR) */

      else if (!strcmp(fluxbc->desc->name1,"YFLUX_H2O_ANODE"))
            {  /* H2O mass flux electro-osmatic drag */
             if (pd->v[pg->imtrx][MASS_FRACTION])
              {
                wspec     = fluxbc->BC_Data_Int[0];
                ai0       = fluxbc->BC_Data_Float[0];
                H         = fluxbc->BC_Data_Float[1];
                cref      = fluxbc->BC_Data_Float[2];
                alphaa    = fluxbc->BC_Data_Float[3];
                alphac    = fluxbc->BC_Data_Float[4];
                T         = fluxbc->BC_Data_Float[5];
                U0        = fluxbc->BC_Data_Float[6];
                nd        = fluxbc->BC_Data_Float[7];

                mass_flux_surf_H2O_CATHODE (mp->mass_flux, mp->d_mass_flux, wspec,
                                            ai0, H, cref, alphac, T, V, U0, nd);
                vnormal += mp->mass_flux[wspec] * StoiCoef[wspec];

                var = MASS_FRACTION;
                for (w = 0; w < pd->Num_Species_Eqn; w++) 
                 {
                  for (j=0; j<ei[pg->imtrx]->dof[var]; j++)
                   {
                     phi_j = bf[var]->phi[j];
                     d_vnorm->C[w][j] += mp->d_mass_flux[wspec][MAX_VARIABLE_TYPES +w] * StoiCoef[wspec] * phi_j;
                   }
                 }
                var = TEMPERATURE;
                if (pd->v[pg->imtrx][var]) 
                 {
                  for (j=0; j<ei[pg->imtrx]->dof[var]; j++)
                   {
                     phi_j = bf[var]->phi[j];
                     d_vnorm->T[j] += mp->d_mass_flux[wspec][TEMPERATURE] * StoiCoef[wspec] * phi_j;
                   }
                 }
                var = VOLTAGE;
                if (pd->v[pg->imtrx][var])
                 {
                   for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];
                      d_vnorm->V[j] += mp->d_mass_flux[wspec][VOLTAGE] * StoiCoef[wspec] * phi_j;
                    }
                 }
              }
           }  /* end of the else on the if(YFLUX_H2O_ANODE) */

      else if (!strcmp(fluxbc->desc->name1,"YFLUX_H2O_CATHODE"))
            {  /* H2O mass flux electro-osmatic drag */
             if (pd->v[pg->imtrx][MASS_FRACTION])
              {
                wspec     = fluxbc->BC_Data_Int[0];
                ai0       = fluxbc->BC_Data_Float[0];
                H         = fluxbc->BC_Data_Float[1];
                cref      = fluxbc->BC_Data_Float[2];
                alphaa    = fluxbc->BC_Data_Float[3];
                alphac    = fluxbc->BC_Data_Float[4];
                T         = fluxbc->BC_Data_Float[5];
                U0        = fluxbc->BC_Data_Float[6];
                nd        = fluxbc->BC_Data_Float[7];

                mass_flux_surf_H2O_ANODE (mp->mass_flux, mp->d_mass_flux, wspec,
                                          ai0, H, cref, alphaa, alphac, T, U0, nd);
                vnormal += mp->mass_flux[wspec] * StoiCoef[wspec];

                var = MASS_FRACTION;
                for (w = 0; w < pd->Num_Species_Eqn; w++) 
                 {
                  for (j=0; j<ei[pg->imtrx]->dof[var]; j++)
                   {
                     phi_j = bf[var]->phi[j];
                     d_vnorm->C[w][j] += mp->d_mass_flux[wspec][MAX_VARIABLE_TYPES +w] * StoiCoef[wspec] * phi_j;
                   }
                 }
                var = TEMPERATURE;
                if (pd->v[pg->imtrx][var]) 
                 {
                  for (j=0; j<ei[pg->imtrx]->dof[var]; j++)
                   {
                     phi_j = bf[var]->phi[j];
                     d_vnorm->T[j] += mp->d_mass_flux[wspec][TEMPERATURE] * StoiCoef[wspec] * phi_j;
                   }
                 }
                var = VOLTAGE;
                if (pd->v[pg->imtrx][var])
                 {
                   for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];
                      d_vnorm->V[j] += mp->d_mass_flux[wspec][VOLTAGE] * StoiCoef[wspec] * phi_j;
                    }
                 }
              }
           }  /* end of the else on the if(YFLUX_H2O_CATHODE) */

      else if (!strcmp(fluxbc->desc->name1,"YFLUX_SULFIDATION")) 
	{
         if(fluxbc->BC_Data_Int[2] != ANNIHILATION_ELECTRONEUTRALITY)
          {
	   if (pd->v[pg->imtrx][MASS_FRACTION])
	    {
	      mode         = fluxbc->BC_Data_Int[2];
	      wspec        = fluxbc->BC_Data_Int[0];
	      nu           = fluxbc->BC_Data_Float[0]; /* stoichiometry coefficient  */
	      k1           = fluxbc->BC_Data_Float[1]; /* forward rate constant      */
	      E1           = fluxbc->BC_Data_Float[2]; /* forward activation energy  */
	      kn1          = fluxbc->BC_Data_Float[3]; /* backward rate constant     */
	      En1          = fluxbc->BC_Data_Float[4]; /* backward activation energy */
	      T            = fluxbc->BC_Data_Float[5]; /* Temperature                */
	      c_H2S        = fluxbc->BC_Data_Float[6]; /* bulk concentration of H2S  */
	      c_O2         = fluxbc->BC_Data_Float[7]; /* bulk concentration of O2   */
	      M_solid      = fluxbc->BC_Data_Float[8]; /* molecular weight of Cu2S   */
	      rho_solid    = fluxbc->BC_Data_Float[9]; /* density of Cu2S            */
	      molar_volume = M_solid/rho_solid;              /* molar volume of Cu2S       */

	      mass_flux_surf_SULFIDATION (mp->mass_flux, mp->d_mass_flux, mode, wspec,   
					  nu, k1, E1, kn1, En1, T, c_H2S, c_O2); 
	      if(fv->c[wspec] <= 0.0)
		{
		  fv->c[wspec] = 1.0e-10;
		}

	      if(mode == SOLID_DIFFUSION_SIMPLIFIED)
		{
		  StoiCoef[wspec] =  -0.5;   /* 0.5 mole of Cu2S produced  */
		  /* per mole of Cu consumped   */    
		}
	      else if(mode == SOLID_DIFFUSION_ELECTRONEUTRALITY || 
		      mode == SOLID_DIFFUSION || 
		      mode == SOLID_DIFFUSION_ELECTRONEUTRALITY_LINEAR)
		{
		  StoiCoef[wspec] =  0.5;    /* 0.5 mole of Cu2S produced          */
		  /* per mole of Cu vacancies generated */    
		}
	      else if(mode == GAS_DIFFUSION)
		{
		  StoiCoef[wspec] =  -1.0;   /* 1 mole of Cu2S produced    */
		  /* per mole of H2S consumped  */    
		}
	      else if(mode == FULL)
		{
		  fprintf(stderr, "The full model has not yet implemented - awaits future efforts\n");
		  exit(1);
		}

	      vnormal += molar_volume*mp->mass_flux[wspec] * StoiCoef[wspec];
              
              var = MASS_FRACTION;
              for (w = 0; w < pd->Num_Species_Eqn; w++) 
                {
		  for (j=0; j<ei[pg->imtrx]->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];
                      d_vnorm->C[w][j] += molar_volume*mp->d_mass_flux[wspec][MAX_VARIABLE_TYPES +w] * StoiCoef[wspec] * phi_j; 
                    }
                }
              var = TEMPERATURE;
	      if (pd->v[pg->imtrx][var]) 
	        {
                  for (j=0; j<ei[pg->imtrx]->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];
                      d_vnorm->T[j] += molar_volume*mp->d_mass_flux[wspec][TEMPERATURE] * StoiCoef[wspec] * phi_j; 
                    }
                }
	    }
	    
          }
	}  /* end of if else on YFLUX_SULFIDATION */ 
      else if (!strcmp(fluxbc->desc->name1,"YFLUX_USER"))  
	{
	  
	  if (pd->v[pg->imtrx][MASS_FRACTION])
	    {
	      wspec           = fluxbc->BC_Data_Int[0];
	      mass_flux_user_surf (mp->mass_flux, mp->d_mass_flux,
				   wspec,
				   fluxbc->u_BC, 
				   tt); 

	      vnormal += mp->mass_flux[wspec] * StoiCoef[wspec];
	      var = MASS_FRACTION;
              for (w = 0; w < pd->Num_Species_Eqn; w++) 
                {
		  for (j=0; j<ei[pg->imtrx]->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];
                      d_vnorm->C[w][j] += mp->d_mass_flux[wspec][MAX_VARIABLE_TYPES +w] * StoiCoef[wspec] * phi_j; 
                    }
                }
              var = TEMPERATURE;
	      if (pd->v[pg->imtrx][var]) 
	        {
                  for (j=0; j<ei[pg->imtrx]->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];
                      d_vnorm->T[j] += mp->d_mass_flux[wspec][TEMPERATURE] * StoiCoef[wspec] * phi_j; 
                    }
                }
	    }
	
	}  /*else on the if(YFLUX_USER***) */
      else if (!strcmp(fluxbc->desc->name1,"YFLUX_ALLOY"))  
	{
	  
	  if (pd->v[pg->imtrx][MASS_FRACTION])
	    {
	      wspec           = fluxbc->BC_Data_Int[0];
	      EH(-1, "KIN_LEAK: no yflux_alloy implemented. See source.");
	      
	      /*If you want yflux_alloy to participate in kin_leak 
	       *mass loss, you need to shore up mass_flux_alloy_surf
	       *with func level, ala YFLUX and YFLUX_USER */
	    }
	}  /*else on the if(YFLUX_ALLLOY***) */
      else  /* This is the YFLUX default in case you weren't paying attention */
	{
	  
	  if (pd->v[pg->imtrx][MASS_FRACTION])
	    {
	      wspec           = fluxbc->BC_Data_Int[0];
	      mass_tran_coeff = fluxbc->BC_Data_Float[0];
	      Y_c             = fluxbc->BC_Data_Float[1];
	      mass_flux_surf_mtc (mp->mass_flux,mp->d_mass_flux,fv->T, 
				  fv->c, wspec, mass_tran_coeff, Y_c);
	      vnormal += mp->mass_flux[wspec] * StoiCoef[wspec];
/*
	      fprintf(stderr,"vnormal=%g, mass_trans_coeff=%g, xbulk=%g, Y_c=%g, StoiCoef=%g\n",
	        vnormal,mass_tran_coeff,fv->c[wspec],Y_c,StoiCoef[wspec]);
*/

              var = MASS_FRACTION;
              for (w = 0; w < pd->Num_Species_Eqn; w++) 
                {
		  for (j=0; j<ei[pg->imtrx]->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];
                      d_vnorm->C[w][j] += mp->d_mass_flux[wspec][MAX_VARIABLE_TYPES +w] * StoiCoef[wspec] * phi_j; 
                    }
                }
              var = TEMPERATURE;
	      if (pd->v[pg->imtrx][var]) 
	        {
                  for (j=0; j<ei[pg->imtrx]->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];
                      d_vnorm->T[j] += mp->d_mass_flux[wspec][TEMPERATURE] * StoiCoef[wspec] * phi_j; 
                    }
                }
	    }
	}  /*else on the if(YFLUX_EQUIL) */

      if ( fluxbc->BC_Data_Int[1] == -1 ) fluxbc = NULL;
      else fluxbc = BC_Types + fluxbc->BC_Data_Int[1];

    } /*while ( fluxbc != NULL ) */

    if(add_fluxes)
      {
      for (p = 0; p < VIM; p++) {
	if ( cr->MassFluxModel == FICKIAN  ||
	     cr->MassFluxModel == STEFAN_MAXWELL  ||
	     cr->MassFluxModel == STEFAN_MAXWELL_CHARGED  ||
	     cr->MassFluxModel == STEFAN_MAXWELL_VOLUME ) 
	{
	  if ( Diffusivity() )  EH( -1, "Error in Diffusivity.");
	  for (w=0; w<pd->Num_Species_Eqn; w++)
	  {
	    vnormal -= fv->snormal[p]*mp->diffusivity[w]*fv->grad_c[w][p]*StoiCoef[w];
	  }
	}
	else if ( cr->MassFluxModel == GENERALIZED_FICKIAN)
	{
	  if ( Generalized_Diffusivity() )  EH( -1, "Error in Diffusivity.");
	  for (w=0; w<pd->Num_Species_Eqn; w++)
	  {
	    for (w1=0; w1<pd->Num_Species_Eqn; w1++)
	    {
	      vnormal -= fv->snormal[p]*mp->diffusivity_gen_fick[w][w1]
		  * fv->grad_c[w1][p]*StoiCoef[w];
	    }
	  }
	}
	else  if ( cr->MassFluxModel == DARCY )
	{ /* diffusion induced convection is zero */
	}
	else
	{
	  EH( -1, "Unimplemented mass flux constitutive relation.");
	}
      }
      if ( af->Assemble_Jacobian && d_vnorm != NULL )
      {
	for (p=0; p<VIM; p++)
	{
	  var = MESH_DISPLACEMENT1;
	  if (pd->v[pg->imtrx][var])
	  {
	    for (q=0; q<VIM; q++)
	    {
	      for ( i=0; i<ei[pg->imtrx]->dof[var]; i++)
	      {
		if ( cr->MassFluxModel == FICKIAN )
		{
		  if ( Diffusivity() )  EH( -1, "Error in Diffusivity.");
				      
		  for (w=0; w<pd->Num_Species_Eqn; w++)
		  {
		    d_vnorm->X[q][i] -= fv->snormal[p]*mp->diffusivity[w] *
			fv->d_grad_c_dmesh[p][w] [q][i]*StoiCoef[w];
		    d_vnorm->X[q][i] -= fv->dsnormal_dx[p][q][i]*mp->diffusivity[w] *
			fv->grad_c[w][p]*StoiCoef[w];
		  }
		}
	        else if ( cr->MassFluxModel == GENERALIZED_FICKIAN)
	        {
	          if ( Generalized_Diffusivity() )  EH( -1, "Error in Diffusivity.");
	          for (w=0; w<pd->Num_Species_Eqn; w++)
	          {
	           for (w1=0; w1<pd->Num_Species_Eqn; w1++)
	             {
		         d_vnorm->X[q][i] -= fv->snormal[p]*
                               mp->diffusivity_gen_fick[w][w1] *
                                fv->d_grad_c_dmesh[p][w1] [q][i]*StoiCoef[w];
		         d_vnorm->X[q][i] -= fv->dsnormal_dx[p][q][i]*
                               mp->diffusivity_gen_fick[w][w1] *
                                fv->grad_c[w1][p]*StoiCoef[w];
	             }
	          }
	        }
	      }
	    }
	  }
	  
	  var = MASS_FRACTION;
	  if (pd->v[pg->imtrx][var])
	  {
	    if ( cr->MassFluxModel == FICKIAN ||
		 cr->MassFluxModel == STEFAN_MAXWELL ||
		 cr->MassFluxModel == STEFAN_MAXWELL_CHARGED ||
		 cr->MassFluxModel == STEFAN_MAXWELL_VOLUME )
	    {
	      if ( Diffusivity() )  EH( -1, "Error in Diffusivity.");

	      for (w=0; w<pd->Num_Species_Eqn; w++)
	      {
		for ( i=0; i<ei[pg->imtrx]->dof[var]; i++)
		{
		  d_vnorm->C[w][i] -=  fv->snormal[p]*StoiCoef[w]*
		      mp->diffusivity[w] * bf[var]->grad_phi[i][p];
		  for (w1=0; w1<pd->Num_Species_Eqn; w1++)
		  {
		    d_vnorm->C[w][i] -= fv->snormal[p]*
			mp->d_diffusivity[w][MAX_VARIABLE_TYPES + w1] 
			* fv->grad_c[w][p]*StoiCoef[w] ;
		  }
		}
	      }
	    }
	    else if ( 0 && cr->MassFluxModel == GENERALIZED_FICKIAN)
	    {
	     if ( Generalized_Diffusivity() )  EH( -1, "Error in Diffusivity.");
	     for (w=0; w<pd->Num_Species_Eqn; w++)
	     {
	      for (w1=0; w1<pd->Num_Species_Eqn; w1++)
	       {
		for ( i=0; i<ei[pg->imtrx]->dof[var]; i++)
		{
		  d_vnorm->C[w][i] -=  fv->snormal[p]*StoiCoef[w]*
		      mp->diffusivity_gen_fick[w][w1] * bf[var]->grad_phi[i][p];
	        for (w2=0; w2<pd->Num_Species_Eqn; w2++)
	         {
		 d_vnorm->C[w][i] -= fv->snormal[p]*StoiCoef[w]*
			mp->d_diffusivity_gf[w][w2][MAX_VARIABLE_TYPES + w1] 
			* fv->grad_c[w2][p] ;  
	         }
	        }
	       }
	     }
	    }
	  }
	} /* end of loop over vconv directions */
      } /* end of if Assemble Jacobian */
    }
  *vnorm = vnormal;

}
/*****************************************************************************/
/*****************************************************************************/

void 
compute_leak_velocity_heat(double *vnorm,
			   NORMAL_VELOCITY_DEPENDENCE_STRUCT *d_vnorm,
			   dbl tt,		/* parameter to vary time integration from 
						 * explicit (tt = 1) to implicit (tt = 0)    */
			   dbl dt,		/* current value of the time step            */
			   struct Boundary_Condition *bc,
			   struct Boundary_Condition *fluxbc)

/******************************************************************************
*
*  Function which evaluates the normal leak velocity used in the
*  kinematic boundary condition n.(v -vs)=0
*/
{
  int a, b, j;
  int var;
  int dim;
  double latent_heat_vap, local_q;
  double vnormal;
  dbl k;                                /* Thermal conductivity. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct; /* Thermal conductivity dependence. */
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  dim   = pd->Num_Dim;
  
  /* Need to calculate vnorm and it's derivatives here.
     this depends on whether YFLUX_BC or YFLUX_BV_BC or YFLUX_EQUIL_BC is chosen */

  vnormal = 0.;
  
  memset(d_vnorm->v, 0, sizeof(dbl)*DIM*MDE);
  memset(d_vnorm->T, 0, sizeof(dbl)*MDE);
  memset(d_vnorm->C, 0, sizeof(dbl)*MAX_CONC*MDE);
  memset(d_vnorm->V, 0, sizeof(dbl)*MDE);
  memset(d_vnorm->F, 0, sizeof(dbl)*MDE);
  memset(d_vnorm->X, 0, sizeof(dbl)*MDE*DIM);

  /* Calculate volume flux of bulk component through surface */
  latent_heat_vap= fluxbc->BC_Data_Float[0];
 
  /*Add heat_residual here divided by latent_heat_vap */

  k   = conductivity( d_k, tran->time_value);
  local_q = 0.;

  for (a=0; a<VIM; a++)
    {
      local_q +=  -k * fv->snormal[a] * fv->grad_T[a] ;
    }  
  vnormal += -local_q/latent_heat_vap;

  if (af->Assemble_Jacobian) {
    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var])
      {
	for (j=0; j<ei[pg->imtrx]->dof[var]; j++) {
	  for (a=0; a<VIM; a++)
	    {
	      d_vnorm->T[j] +=  k * fv->snormal[a] * bf[var]->grad_phi[j][a]/latent_heat_vap ;
	    }  
	}
      }
    for (b = 0; b < dim; b++) {
      var = MESH_DISPLACEMENT1 + b;
      if (pd->v[pg->imtrx][var] && fluxbc->BC_Name != LS_LATENT_HEAT_BC ) {
	for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
		      
	  for (a = 0; a < dim; a++) {
	    d_vnorm->X[b][j] +=
	      k * fv->dsnormal_dx[a][b][j] * fv->grad_T[a]/latent_heat_vap;
	  }
	}
      }
    }
  }

  
 
    
  *vnorm = vnormal;

}
/*****************************************************************************/

void 
compute_leak_energy(double *enorm,
	              NORMAL_ENERGY_DEPENDENCE_STRUCT *d_enorm,
	              dbl tt,		/* parameter to vary time integration from 
		        		 * explicit (tt = 1) to implicit (tt = 0)    */
	              dbl dt,		/* current value of the time step            */
	              struct Boundary_Condition *bc,
		      struct Boundary_Condition *fluxbc)

/******************************************************************************
*
*  Function which evaluates the normal leak energy used in the
*  kinematic boundary condition n_dot_q.  This is a direct clone of
*  compute_leak_velocity where the energy contribution is calculated.
*/
{
  int i, j;
  int var;
  int w, wspec;
  int num_comp;
  double StoiCoef[MAX_CONC];
  double mass_tran_coeff,Y_c;
  double xbulk, lat_heat_bulk;
  double enormal, phi_j;


  int mode;
  double amb_pres, A, mtc, Y_inf, driving_force;
  double d_mtc[MAX_VARIABLE_TYPES+MAX_CONC];
  double activity[MAX_CONC];
  double dact_dC[MAX_CONC][MAX_CONC];
  double psat[MAX_CONC], dpsatdt[MAX_CONC];

  /* local contributions of boundary condition to residual and jacobian */

/***************************** EXECUTION BEGINS *******************************/

  /* Need to calculate vnorm and it's derivatives here.
     this depends on whether YFLUX_BC or YFLUX_BV_BC or YFLUX_EQUIL_BC is chosen */

  enormal = 0.;
  
  memset(d_enorm->v, 0, sizeof(dbl)*DIM*MDE);
  memset(d_enorm->T, 0, sizeof(dbl)*MDE);
  memset(d_enorm->C, 0, sizeof(dbl)*MAX_CONC*MDE);
  memset(d_enorm->V, 0, sizeof(dbl)*MDE);
  memset(d_enorm->F, 0, sizeof(dbl)*MDE);

  /* Expand function use to include KIN_CHEM BC */

  if ( bc != NULL && bc->BC_Name == KIN_CHEM_BC ) {
    num_comp = bc->len_u_BC;
    for (i=0; i < num_comp; i++) {
      StoiCoef[i] = bc->u_BC[i]; }
  }
  else {
    for (i=0; i < MAX_CONC; i++) {
      StoiCoef[i] = 1.0; }
  }

  /* Calculate energy flux of bulk component */
  xbulk = 1.;
  if (pd->v[pg->imtrx][MASS_FRACTION])
    for (w=0; w<pd->Num_Species_Eqn; w++) xbulk -= fv->c[w];
  
  /* Calculate volume flux of bulk component through surface */
  lat_heat_bulk   = bc->BC_Data_Float[0];
  mass_tran_coeff = bc->BC_Data_Float[1];
  Y_c             = bc->BC_Data_Float[2];
  
  if (af->Assemble_Jacobian) {
    var = MASS_FRACTION;
    if (pd->v[pg->imtrx][var])
      {
	for (w=0; w<pd->Num_Species_Eqn; w++) {
	  for (j=0; j<ei[pg->imtrx]->dof[var]; j++) {
	    phi_j = bf[var]->phi[j];
	    d_enorm->C[w][j] = - mp->density * lat_heat_bulk * mass_tran_coeff * phi_j;
	  }
	}
      }
  }
  enormal +=  mp->density * lat_heat_bulk * mass_tran_coeff * ( xbulk - Y_c );
  
  while ( fluxbc != NULL )
    {
      
      if(!strcmp(fluxbc->desc->name1,"YFLUX_EQUIL"))  /* this is limited to multicomponent */
	{
	  
	  /* call routine to calculate surface flux of bulk component and it's 
	   * sensitivity to all variable types    
           * ACS: modified 10/99 to accommodate mass conc. formulation for YFLUX_EQUIL  */
	  
	  driving_force = 1.;	  
	  
	  if (pd->v[pg->imtrx][MASS_FRACTION])
	    {
	      wspec    = fluxbc->BC_Data_Int[0];
	      mode     = fluxbc->BC_Data_Int[2];
	      amb_pres = fluxbc->BC_Data_Float[0];
	      mtc      = fluxbc->BC_Data_Float[1];
	      Y_inf    = fluxbc->BC_Data_Float[2];

              /*  Chilton-Coburn correlation if FLORY_CC */
              memset(d_mtc, 0, sizeof(dbl)*(MAX_VARIABLE_TYPES+MAX_CONC));
	      if( mode == FLORY_CC)
 		{
		  mtc_chilton_coburn(&mtc, d_mtc, wspec, 
 				fluxbc->BC_Data_Float[1],
 				fluxbc->BC_Data_Float[3],
 				fluxbc->BC_Data_Float[0],
 				fluxbc->BC_Data_Float[4]);
 		}
	      /* Shouldn't this be here, PRS 3/19/01. cf. KIN_LEAK
	       *
	       * mtc = mtc*mp->specific_volume[wspec];   
	       * 
	       * Adding it changes
	       * Amy's benchmark, which did agree with experiment. 
	       *   It's required in KIN_LEAK because a volume flux is required.
	       *   Elsewhere, a mass flux is required. - RBS (5/6/03)
	       */
		  
	      /* Nonideal VP Calculations based on either ANTOINE or RIEDEL models */
	      
	      if(mp->VaporPressureModel[wspec] == ANTOINE )
		{
		  antoine_psat(wspec, mp->u_vapor_pressure[wspec],
			       &psat[wspec], &dpsatdt[wspec]);
		  mp-> vapor_pressure[wspec] = psat[wspec];
		}
	      
	      else if(mp->VaporPressureModel[wspec] == RIEDEL )
		{
		  riedel_psat(wspec, mp->u_vapor_pressure[wspec],
			      &psat[wspec], &dpsatdt[wspec]);
		  mp-> vapor_pressure[wspec] = psat[wspec];
		}
		  
	      /* Calculate mole flux of other components through surface */
	      mass_flux_equil_mtc (mp->mass_flux,mp->d_mass_flux, activity, dact_dC, 
				   fv->c, mode, amb_pres, wspec, mtc, d_mtc, Y_inf);
		  
	      A = mp->vapor_pressure[wspec]/amb_pres;
	      driving_force -=  A*activity[wspec];
	      enormal += mp->density*mp->latent_heat_vap[wspec]*mp->mass_flux[wspec]
                        *StoiCoef[wspec];


	      var = MASS_FRACTION;
	      if (pd->v[pg->imtrx][var])
		{
		  for (w=0; w<pd->Num_Species_Eqn; w++) {
		    for (j=0; j<ei[pg->imtrx]->dof[var]; j++) {
		      phi_j = bf[var]->phi[j];
		      d_enorm->C[w][j] = - mp->density*mp->latent_heat_vap[wspec]*StoiCoef[wspec]
			                   *mp->d_mass_flux[wspec][MAX_VARIABLE_TYPES+w]*phi_j ;
		    }
		  }
		}

	      var = TEMPERATURE;
	      if (pd->v[pg->imtrx][var])
	        {
                  for (j=0; j<ei[pg->imtrx]->dof[var]; j++) {
                    phi_j = bf[var]->phi[j];
                    d_enorm->T[j] += mp->density*mp->latent_heat_vap[wspec]*StoiCoef[wspec]
		                     *mp->d_mass_flux[wspec][TEMPERATURE] * phi_j;
                  }
                }
	    }
	  	  
	} /*  if (!strcmp(fluxbc->desc->name1,"YFLUX_EQUIL") */
      else if (!strcmp(fluxbc->desc->name1,"YFLUX_CONST"))  
	{	  
	  if (pd->v[pg->imtrx][MASS_FRACTION])
	    {
	      wspec    = fluxbc->BC_Data_Int[0];
	      enormal -= mp->density*mp->latent_heat_vap[wspec]*fluxbc->BC_Data_Float[0] * StoiCoef[wspec];
	    }
	} /*  if (!strcmp(fluxbc->desc->name1,"YFLUX_CONST") */
      else if (!strcmp(fluxbc->desc->name1,"YFLUX_USER"))  
	{
	  if (pd->v[pg->imtrx][MASS_FRACTION])
	    {
	      wspec           = fluxbc->BC_Data_Int[0];
	      mass_flux_user_surf (mp->mass_flux, mp->d_mass_flux,
				   wspec,
				   fluxbc->u_BC, 
				   tt); 

	      enormal +=  mp->density*mp->latent_heat_vap[wspec]*mp->mass_flux[wspec] * StoiCoef[wspec];

	      var = MASS_FRACTION;
              for (w = 0; w < pd->Num_Species_Eqn; w++) 
                {
		  for (j=0; j<ei[pg->imtrx]->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];
                      d_enorm->C[w][j] +=  mp->density*mp->latent_heat_vap[wspec]
			*mp->d_mass_flux[wspec][MAX_VARIABLE_TYPES +w] * StoiCoef[wspec] * phi_j; 
                    }
                }
              var = TEMPERATURE;
	      if (pd->v[pg->imtrx][var]) 
	        {
                  for (j=0; j<ei[pg->imtrx]->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];
                      d_enorm->T[j] += mp->density*mp->latent_heat_vap[wspec]
			               *mp->d_mass_flux[wspec][TEMPERATURE] * StoiCoef[wspec] * phi_j; 
                    }
                }
	    }
	
	}  /*else on the if(YFLUX_USER***) */
      else  /* This is the YFLUX default in case you weren't paying attention */
	{
	  
	  if (pd->v[pg->imtrx][MASS_FRACTION])
	    {
	      wspec           = fluxbc->BC_Data_Int[0];
	      mass_tran_coeff = fluxbc->BC_Data_Float[0];
	      Y_c             = fluxbc->BC_Data_Float[1];
	      mass_flux_surf_mtc (mp->mass_flux,mp->d_mass_flux,fv->T, 
				  fv->c, wspec, mass_tran_coeff, Y_c);

	      enormal += mp->density*mp->latent_heat_vap[wspec]*mp->mass_flux[wspec] * StoiCoef[wspec];
/*
	      fprintf(stderr,"vnormal=%g, mass_trans_coeff=%g, xbulk=%g, Y_c=%g, StoiCoef=%g\n",
	        vnormal,mass_tran_coeff,fv->c[wspec],Y_c,StoiCoef[wspec]);
*/
              var = MASS_FRACTION;
              for (w = 0; w < pd->Num_Species_Eqn; w++) 
                {
		  for (j=0; j<ei[pg->imtrx]->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];
                      d_enorm->C[w][j] += mp->density*mp->latent_heat_vap[wspec]
			                  *mp->d_mass_flux[wspec][MAX_VARIABLE_TYPES +w] * StoiCoef[wspec] * phi_j; 
                    }
                }
              var = TEMPERATURE;
	      if (pd->v[pg->imtrx][var]) 
	        {
                  for (j=0; j<ei[pg->imtrx]->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];
                      d_enorm->T[j] += mp->density*mp->latent_heat_vap[wspec]
			               *mp->d_mass_flux[wspec][TEMPERATURE] * StoiCoef[wspec] * phi_j; 
                    }
                }
	    }
	}  /*else default YFLUX */

      if (fluxbc->BC_Data_Int[1] == -1 )fluxbc = NULL;
      else fluxbc = BC_Types + fluxbc->BC_Data_Int[1];

    } /*while ( fluxbc != NULL ) */
    
  *enorm = enormal;

}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void 
kin_bc_leak(double func[],
	    double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
	    double x_dot[DIM],  /* mesh velocity                             */
	    dbl tt,		/* parameter to vary time integration from 
				 * explicit (tt = 1) to implicit (tt = 0)    */
	    dbl dt,		/* current value of the time step            */
	    int bc_input_id,
	    struct Boundary_Condition *BC_Types)

/******************************************************************************
*
*  Function which evaluates the kinematic boundary condition n.(v -vs)=0
*  ALSO Function which evaluates the velo-normal boundary condition n.(v -vs)=0
*
*            Author: P. R. Schunk    (3/24/94)
*
*  With leak caused by flux of species through surface n.(v-vs)=sum(evap vol flux)
*			 Reviser: R. A. Cairncross (11/1/94)
*
*  With ties to equilibrium mole fraction of species by calling Raoult's law
*  or Flory-Huggins VLE relations.
*			 Reviser: A. C. Sun (9/20/98)
*
*  With leak caused by species mass flux that is given by Butler-Volmer kinetics.
*  Also introduced M_solid (molecular weight of solid deposit) and
*  rho_solid (density of solid deposit) to modify vnormal for processes
*  involving a moving boundary due to formation of a solid material 
*  (as in LIGA electrodeposition).  
*
*			 Reviser: K. S. Chen (11/16/2000, 11/21/2000, 10/1/2001)
*
*  With leak caused by copper sulfidation reaction with surface rate or 
*  species mass flux given by 
*
*
*  Simplified kinetic model:
*
*  Surface Rate or species mass flux = n.(v - vs) = k exp(-E/R/T) c_H2S c_Cu 
*
*
*  General model:
*
*  Species mass flux = k1 exp(-E1/R/T) c_H2S c_O2**0.5 - kn1 exp(-En1/R/T) c_V c_h  
*
*                        Reviser: K. S. Chen (3/20/2002)
*
*  ----------------------------------------------------------------------------
*  Take out all of the flux calculations and move them to compute_leak_velocity
*  to facilitate LS implementation of kin_leak.
*                        Reviser: D. Noble (5/2005) 
*
*  Functions Called
*
*  compute_leak_velocity - compute summation of species fluxes.
*            sum_over_i(k_i*(y_i - y_inf_i))+k_n*((1-sum_over_i(y_i))-y_inf_n)
*
*  get_convection_velocity - compute the relative convective velocity.
*  ----------------------------------------------------------------------------
*******************************************************************************/
{
  int j, j_id;
  int var,jvar,kdir;
  int w, dim;
  double phi_j;

  double vnorm; /*Calculated normal velocity */
  NORMAL_VELOCITY_DEPENDENCE_STRUCT d_vnorm_struct;
  NORMAL_VELOCITY_DEPENDENCE_STRUCT *d_vnorm = &d_vnorm_struct;
  
  double vconv[MAX_PDIM]; /*Calculated convection velocity */
  double vconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;
  int err;
  struct Boundary_Condition *bc = NULL;
  struct Boundary_Condition *fluxbc = NULL;

  /* local contributions of boundary condition to residual and jacobian */

/***************************** EXECUTION BEGINS *******************************/

  if(af->Assemble_LSA_Mass_Matrix)
    {
      for(kdir = 0; kdir < pd->Num_Dim; kdir++)
	{
	  var = MESH_DISPLACEMENT1 + kdir;
	  if(pd->v[pg->imtrx][var])
	    for(j = 0; j < ei[pg->imtrx]->dof[var]; j++)
	      {
		phi_j = bf[var]->phi[j];
		d_func[0][var][j] -= phi_j * fv->snormal[kdir];
	      }
	}
      return;
    }

  dim   = pd->Num_Dim;
  
  bc = BC_Types + bc_input_id;
  if ( bc->BC_Data_Int[1] == -1 ) fluxbc = NULL;
  else fluxbc = BC_Types + bc->BC_Data_Int[1];

  /* Get leak velocity */
  if (bc->BC_Name == KIN_LEAK_HEAT_BC)
    {
      compute_leak_velocity_heat(&vnorm, d_vnorm, tt, dt, bc, fluxbc);
    }
  else
    {
      compute_leak_velocity(&vnorm, d_vnorm, tt, dt, bc, fluxbc);
    }
  
  /* get the convection velocity (it's different for arbitrary and
     lagrangian meshes) */
  /* get deformation gradient, etc. if lagrangian mesh with inertia */
  if ((pd->MeshMotion == LAGRANGIAN ||
       pd->MeshMotion == DYNAMIC_LAGRANGIAN) && pd->MeshInertia == 1)
    {
      err = belly_flop(elc->lame_mu);
      EH(err, "error in belly flop");
      if( neg_elem_volume ) return;
    }
  err =   get_convection_velocity(vconv, vconv_old, d_vconv, dt, tt);
  EH(err, "Error in calculating effective convection velocity");

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
		    d_func[0][var][j_id] += vconv[kdir]*fv->dsnormal_dx[kdir][jvar][j_id]
		      + d_vconv->X[kdir][jvar][j_id]*fv->snormal[kdir];
		    
		  } 
	       	d_func[0][var][j_id] += -d_vnorm->X[jvar][j_id]; 
	      }
	  }
      }
    
    for(jvar=0; jvar<dim; jvar++) 
      {
	var = VELOCITY1 + jvar;
        if (pd->v[pg->imtrx][var])
          {
	    for (j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
	      {
                for(kdir=0; kdir<dim; kdir++) 
	          {
                    d_func[0][var][j_id] += d_vconv->v[kdir][jvar][j_id] * fv->snormal[kdir];
                  }
	      }
	  }
      }
    
    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var]) 
      {
	for (j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
	  {
	    d_func[0][var][j_id] -= d_vnorm->T[j_id];
	    for(kdir=0; kdir<dim; kdir++) 
	      {
		d_func[0][var][j_id] += d_vconv->T[kdir][j_id] * fv->snormal[kdir];
	      }
	  }
      }

    var = VOLTAGE;
    if (pd->v[pg->imtrx][var]) {
      for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
	d_func[0][var][j_id] -= d_vnorm->V[j_id];
      }
    }
    
    var = MASS_FRACTION;
    if (pd->v[pg->imtrx][var]) 
      {
	for (w = 0; w < pd->Num_Species_Eqn; w++) 
	  {
	    for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
	      {
		d_func[0][MAX_VARIABLE_TYPES + w][j_id] -= d_vnorm->C[w][j_id];
		for(kdir=0; kdir<dim; kdir++) 
		  {
		    d_func[0][MAX_VARIABLE_TYPES + w][j_id] += 
		      d_vconv->C[kdir][w][j_id] * fv->snormal[kdir];
		  }
	      }
	  }
      }
  }
  /* Calculate the residual contribution	*/
  
  *func = -vnorm;
  for (kdir = 0; kdir < dim; kdir++) {
    *func += vconv[kdir] * fv->snormal[kdir];
  }

} /* END of routine kin_bc_leak  */
/******************************************************************************/

/*****************************************************************************/

void
kin_bc_electrodeposition (double func[],
			  double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
			  double time, /* current time */
			  dbl tt,  /* parameter to vary time integration from explicit (tt = 1)
				      to implicit (tt = 0) */
			  dbl dt,  /* current value of the time step */
			  int bc_input_id,
			  struct Boundary_Condition *BC_Types)
     
/******************************************************************************
*
*  Function which evaluates the boundary condition for mesh motion arising
*  from electrodeposition.  This routine adapted from kin_bc_leak and
*  get_convection_velocity.
*
*            Author: R. S. Larson    (5/27/02)
*
*  ----------------------------------------------------------------------------
*
*  Functions Called
*
*  mass_flux_surf_BV2 - calculates the flux of a given species to the surface
*                       using general Butler-Volmer kinetics.
*
*  mass_flux_surf_NI  - calculates the flux of a given species to the surface
*                       for the special case of nickel electrodeposition.
*
*******************************************************************************/

{

/* Local variables */

  int j, j_id, p, i;
  int var, jvar, kdir;
  int w, wspec, ibc, dim;
  double d_vsurface[MAX_VARIABLE_TYPES][MAX_CONC];
  double vsurface;
  double phi_j;

  double nxdot[DIM]; /* negative of mesh velocity */
  double dnxdot_dx[DIM][DIM][MDE]; /* sensitivity of nxdot to nodal x */

  int flag;
  double PHI_E, T, Vmolar;
  double k, alphaa, alphac, U0;

  /* local contributions of boundary condition to residual and jacobian */
 
/***************************** EXECUTION BEGINS *******************************/

  if (af->Assemble_LSA_Mass_Matrix)
     {
      for (kdir = 0; kdir < pd->Num_Dim; kdir++)
         {
          var = MESH_DISPLACEMENT1 + kdir;
          if (pd->v[pg->imtrx][var])
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
                 {
                  phi_j = bf[var]->phi[j];
                  d_func[0][var][j] -= phi_j * fv->snormal[kdir];
                 }
         }
      return;
     }

  dim   = pd->Num_Dim;
 
  /* calculate vsurface and its derivatives. */

  vsurface = 0.;
  memset(d_vsurface, 0, sizeof(dbl)*MAX_VARIABLE_TYPES*MAX_CONC)    ;

  /* using pointer to boundary condition types in BC_Data_Int  */
  ibc = BC_Types[bc_input_id].BC_Data_Int[1];
 
  if (!strcmp(BC_Types[ibc].desc->name1,"YFLUX_NI"))
     {
      if (pd->v[pg->imtrx][MASS_FRACTION])
         {
          Vmolar = 6.596;  /*  molar volume of solid nickel  */
          while (ibc != -1)
             {
              wspec = BC_Types[ibc].BC_Data_Int[0];    /*  species number           */
              flag  = 0;                               /*  need species flux        */
              PHI_E = BC_Types[ibc].BC_Data_Float[0];  /*  electrode potential      */
              T     = BC_Types[ibc].BC_Data_Float[1];  /*  electrolyte temperature  */
              if (wspec == 0 || wspec == 4)
                 {
		   mass_flux_surf_NI (mp->mass_flux, mp->d_mass_flux, time,
				      wspec, flag, PHI_E, T);
                  vsurface += mp->mass_flux[wspec]*Vmolar;
                  for (w=0; w<pd->Num_Species_Eqn; w++)
                     {
                      d_vsurface[MASS_FRACTION][w] +=
                         mp->d_mass_flux[wspec][MAX_VARIABLE_TYPES +w]*Vmolar;
                     }
                  d_vsurface[TEMPERATURE][0] +=
                     mp->d_mass_flux[wspec][TEMPERATURE]*Vmolar;
                  d_vsurface[VOLTAGE][0] += mp->d_mass_flux[wspec][VOLTAGE]*Vmolar;
                 }
              ibc = BC_Types[ibc].BC_Data_Int[1];
             }
         }
     }  /*  end of YFLUX_NI  */

  else if (!strcmp(BC_Types[ibc].desc->name1,"YFLUX_BV2"))
     {
      if (pd->v[pg->imtrx][MASS_FRACTION])
         {
          while (ibc != -1)
             {
              wspec  = BC_Types[ibc].BC_Data_Int[0];    /*  species number
   */
              flag   = 0;                               /*  need species flux
   */
              k      = BC_Types[ibc].BC_Data_Float[0];  /*  rate constant
   */
              PHI_E  = BC_Types[ibc].BC_Data_Float[1];  /*  electrode potential
   */
              alphaa = BC_Types[ibc].BC_Data_Float[2];  /*  anodic transfer coefficient
   */
              alphac = BC_Types[ibc].BC_Data_Float[3];  /*  cathodic transfer coefficient
   */
              U0     = BC_Types[ibc].BC_Data_Float[4];  /*  standard open circuit potential  */
              T      = BC_Types[ibc].BC_Data_Float[5];  /*  electrolyte temperature
   */
              Vmolar = BC_Types[ibc].BC_Data_Float[6];  /*  molar volume of solid deposit
   */
              mass_flux_surf_BV2 (time, mp->mass_flux, mp->d_mass_flux, 
                                  wspec, flag, k, PHI_E, alphaa, alphac, U0, T);
              vsurface += mp->mass_flux[wspec]*Vmolar;
              for (w=0; w<pd->Num_Species_Eqn; w++)
                 {
                  d_vsurface[MASS_FRACTION][w] +=
                     mp->d_mass_flux[wspec][MAX_VARIABLE_TYPES +w]*Vmolar;
                 }
              d_vsurface[TEMPERATURE][0] +=
                 mp->d_mass_flux[wspec][TEMPERATURE]*Vmolar;
              d_vsurface[VOLTAGE][0] += mp->d_mass_flux[wspec][VOLTAGE]*Vmolar;
              ibc = BC_Types[ibc].BC_Data_Int[1];
             }
         }
     }  /*  end of YFLUX_BV2  */

  /* initialize sensitivity arrays if necessary */

  if ( af->Assemble_Jacobian )
     {
      if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) memset(dnxdot_dx, 0, DIM*DIM*MDE*sizeof(dbl));
     }
 
  /* calculate -xdot and its derivatives */
 
  for (p=0; p<VIM; p++)
     {
      nxdot[p] = 0.;
      if ( pd->TimeIntegration != STEADY )
         {
          if (pd->v[pg->imtrx][R_MESH1])
             {
              nxdot[p] -= (1. + 2.*tt) * (fv->x[p] - fv_old->x[p])/dt
                          - 2. * tt * fv_dot->x[p];
             }
         }
     }
 
  if ( af->Assemble_Jacobian )
     {
      for (p=0; p<VIM; p++)
         {
          if ( pd->TimeIntegration != STEADY )
             {
              var = MESH_DISPLACEMENT1+p;
              if (pd->v[pg->imtrx][var])
                 {
                  for ( i=0; i<ei[pg->imtrx]->dof[var]; i++)
                     {
                      dnxdot_dx[p][p][i] -= bf[var]->phi[i] * (1. + 2.* tt)/dt;
                     }
                 }
             }
         }
     }

  /* sum the contributions to the global stiffness matrix */
 
  if (af->Assemble_Jacobian)
     {
      for (jvar=0; jvar<dim; jvar++)
         {
          var = MESH_DISPLACEMENT1 + jvar;
          if (pd->v[pg->imtrx][var])
             {
              for (j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
                 {
                  phi_j = bf[var]->phi[j_id];
                  for (kdir=0; kdir<dim; kdir++)
                     {
                      d_func[0][var][j_id] += nxdot[kdir]*fv->dsnormal_dx[kdir][jvar][j_id]
                                              + dnxdot_dx[kdir][jvar][j_id]*fv->snormal[kdir];
                     }
                 }
             }
         }

      var = TEMPERATURE;
      if (pd->v[pg->imtrx][var])
         {
          for (j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
             {
              phi_j = bf[var]->phi[j_id];
              d_func[0][var][j_id] -= d_vsurface[TEMPERATURE][0] * phi_j;
             }
         }

      var = VOLTAGE;
      if (pd->v[pg->imtrx][var])
         {
          for (j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
             {
              phi_j = bf[var]->phi[j_id];
              d_func[0][var][j_id] -= d_vsurface[VOLTAGE][0] * phi_j;
             }
         }

      var = MASS_FRACTION;
      if (pd->v[pg->imtrx][var])
         {
          for (w = 0; w < pd->Num_Species_Eqn; w++)
             {
              for (j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
                 {
                  phi_j = bf[var]->phi[j_id];
                  d_func[0][MAX_VARIABLE_TYPES + w][j_id] -= d_vsurface[MASS_FRACTION][w]*phi_j;
                 }
             }
         }
     }

  /* calculate the residual contribution        */
 
  *func = -vsurface;
 
  for (kdir=0; kdir<dim; kdir++)
     {
      *func += nxdot[kdir] * fv->snormal[kdir];
     }

} /* END of routine kin_bc_electrodeposition  */

/*****************************************************************************/

void
vnorm_bc_electrodeposition (double func[],
			    double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
			    double time,
			    dbl tt,  /* parameter to vary time integration from explicit (tt = 1)
					to implicit (tt = 0) */
			    dbl dt,  /* current value of the time step */
			    int bc_input_id,
			    struct Boundary_Condition *BC_Types)

/******************************************************************************
*
*  Function which evaluates the boundary condition for normal velocity arising
*  from electrodeposition.  This routine adapted from kin_bc_leak.
*
*            Author: R. S. Larson    (5/30/02)
*
*  ----------------------------------------------------------------------------
*
*  Functions Called
*
*  mass_flux_surf_BV2 - calculates the flux of a given species to the surface
*                       using general Butler-Volmer kinetics.
*
*  mass_flux_surf_NI  - calculates the flux of a given species to the surface
*                       for the special case of nickel electrodeposition.
*
*******************************************************************************/

{
 
/* Local variables */
 
  int j, j_id;
  int var, jvar, kdir;
  int w, wspec, ibc, dim;
  double d_surface_flux[MAX_VARIABLE_TYPES][MAX_CONC];
  double surface_flux;
  double phi_j;

  double vconv[MAX_PDIM]; /*Calculated convection velocity */
  double vconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;
  int err;

  int flag;
  double PHI_E, T, PHI_E_save = 0.0, T_save = -1.0, mw;
  double k, alphaa, alphac, U0;

  dbl rho;
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;  /* density dependence */
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;

  if (MAX_CONC < 7) {
    EH(-1, "vnorm_bc_electrodeposition expects MAX_CONC >= 7");
    return;
  }

  /* local contributions of boundary condition to residual and jacobian */
 
/***************************** EXECUTION BEGINS *******************************/

  if (af->Assemble_LSA_Mass_Matrix)
     {
      for (kdir = 0; kdir < pd->Num_Dim; kdir++)
         {
          var = MESH_DISPLACEMENT1 + kdir;
          if (pd->v[pg->imtrx][var])
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
                 {
                  phi_j = bf[var]->phi[j];
                  d_func[0][var][j] -= phi_j * fv->snormal[kdir];
                 }
         }
      return;
     }

  dim   = pd->Num_Dim;
 
  /* calculate surface_flux and its derivatives. */

  surface_flux = 0.;
  memset(d_surface_flux, 0, sizeof(dbl)*MAX_VARIABLE_TYPES*MAX_CONC)    ;

  /* using pointer to boundary condition types in BC_Data_Int  */
  ibc = BC_Types[bc_input_id].BC_Data_Int[1];
 
  if (!strcmp(BC_Types[ibc].desc->name1,"YFLUX_NI"))
     {
      if (pd->v[pg->imtrx][MASS_FRACTION])
         {
          while (ibc != -1)
             {
              wspec = BC_Types[ibc].BC_Data_Int[0];    /*  species number           */
              flag  = 0;                               /*  need species flux        */
              PHI_E = BC_Types[ibc].BC_Data_Float[0];  /*  electrode potential      */
              T     = BC_Types[ibc].BC_Data_Float[1];  /*  electrolyte temperature  */
              if (wspec == 3)
                 {
                  PHI_E_save = PHI_E;
                  T_save = T;
                 }
              mass_flux_surf_NI (mp->mass_flux, mp->d_mass_flux, time, wspec, flag, PHI_E, T);
              mw = mp->molecular_weight[wspec];
              surface_flux += mp->mass_flux[wspec]*mw;
              for (w=0; w<pd->Num_Species_Eqn; w++)
                 {
                  d_surface_flux[MASS_FRACTION][w] +=
                     mp->d_mass_flux[wspec][MAX_VARIABLE_TYPES +w]*mw;
                 }
              d_surface_flux[TEMPERATURE][0] += mp->d_mass_flux[wspec][TEMPERATURE]*mw;
              d_surface_flux[VOLTAGE][0] += mp->d_mass_flux[wspec][VOLTAGE]*mw;
              ibc = BC_Types[ibc].BC_Data_Int[1];
             }
          /* add in contribution from solvent species (water) */
          wspec = 6;
          flag  = 0;
          mass_flux_surf_NI (mp->mass_flux, mp->d_mass_flux, time, wspec, flag, PHI_E_save, T_save);
          mw = 18.01534;
          surface_flux += mp->mass_flux[wspec]*mw;
          for (w=0; w<pd->Num_Species_Eqn; w++)
             {
              d_surface_flux[MASS_FRACTION][w] +=
                 mp->d_mass_flux[wspec][MAX_VARIABLE_TYPES +w]*mw;
             }
          d_surface_flux[TEMPERATURE][0] += mp->d_mass_flux[wspec][TEMPERATURE]*mw;
          d_surface_flux[VOLTAGE][0] += mp->d_mass_flux[wspec][VOLTAGE]*mw;
         }
     }  /*  end of YFLUX_NI  */

  else if (!strcmp(BC_Types[ibc].desc->name1,"YFLUX_BV2"))
     {
      if (pd->v[pg->imtrx][MASS_FRACTION])
         {
          while (ibc != -1)
             {
              wspec  = BC_Types[ibc].BC_Data_Int[0];    /*  species number
   */
              flag   = 0;                               /*  need species flux
   */
              k      = BC_Types[ibc].BC_Data_Float[0];  /*  rate constant
   */
              PHI_E  = BC_Types[ibc].BC_Data_Float[1];  /*  electrode potential
   */
              alphaa = BC_Types[ibc].BC_Data_Float[2];  /*  anodic transfer coefficient
   */
              alphac = BC_Types[ibc].BC_Data_Float[3];  /*  cathodic transfer coefficient
   */
              U0     = BC_Types[ibc].BC_Data_Float[4];  /*  standard open circuit potential  */
              T      = BC_Types[ibc].BC_Data_Float[5];  /*  electrolyte temperature
   */
              mass_flux_surf_BV2 (time, mp->mass_flux, mp->d_mass_flux,
                                  wspec, flag, k, PHI_E, alphaa, alphac, U0, T);
              mw = mp->molecular_weight[wspec];
              surface_flux += mp->mass_flux[wspec]*mw;
              for (w=0; w<pd->Num_Species_Eqn; w++)
                 {
                  d_surface_flux[MASS_FRACTION][w] +=
                     mp->d_mass_flux[wspec][MAX_VARIABLE_TYPES +w]*mw;
                 }
              d_surface_flux[TEMPERATURE][0] += mp->d_mass_flux[wspec][TEMPERATURE]*mw;
              d_surface_flux[VOLTAGE][0] += mp->d_mass_flux[wspec][VOLTAGE]*mw;
              ibc = BC_Types[ibc].BC_Data_Int[1];
             }
         }
     }  /*  end of YFLUX_BV2  */

  /* get the convection velocity and its derivatives */

  err = get_convection_velocity(vconv, vconv_old, d_vconv, dt, tt);
  EH(err, "Error in calculating effective convection velocity");

  /* get the fluid density and its derivatives */

  rho = density(d_rho, time); /*  RSL 6/22/02  */

  /* sum the contributions to the global stiffness matrix */
 
  if (af->Assemble_Jacobian)
     {
      for (jvar=0; jvar<dim; jvar++)
         {
          var = MESH_DISPLACEMENT1 + jvar;
          for (j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
             {
              if (pd->v[pg->imtrx][var])
                 {
                  phi_j = bf[var]->phi[j_id];
                  for (kdir=0; kdir<dim; kdir++)
                     {
                      d_func[0][var][j_id] += vconv[kdir]*fv->dsnormal_dx[kdir][jvar][j_id]
                                              + d_vconv->X[kdir][jvar][j_id]*fv->snormal[kdir];
                     }
                 }
             }
         }

      for (jvar=0; jvar<dim; jvar++)
         {
          var = VELOCITY1 + jvar;
          for (j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
             {
              if (pd->v[pg->imtrx][var])
                 {
                  phi_j = bf[var]->phi[j_id];
                  d_func[0][var][j_id] += d_vconv->v[jvar][jvar][j_id] * fv->snormal[jvar];
                 }
             }
         }

      var = TEMPERATURE;
      if (pd->v[pg->imtrx][var])
         {
          for (j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
             {
              phi_j = bf[var]->phi[j_id];
              d_func[0][var][j_id] -= d_surface_flux[TEMPERATURE][0] * phi_j / rho;
              d_func[0][var][j_id] += surface_flux * d_rho->T[j_id] / (rho*rho);
              for (kdir=0; kdir<dim; kdir++)
                 {
                  d_func[0][var][j_id] += d_vconv->T[kdir][j_id] * fv->snormal[kdir];
                 }
             }
         }

      var = VOLTAGE;
      if (pd->v[pg->imtrx][var])
         {
          for (j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
             {
              phi_j = bf[var]->phi[j_id];
              d_func[0][var][j_id] -= d_surface_flux[VOLTAGE][0] * phi_j / rho;
             }
         }

      var = MASS_FRACTION;
      if (pd->v[pg->imtrx][var])
         {
          for (w = 0; w < pd->Num_Species_Eqn; w++)
             {
              for (j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
                 {
                  phi_j = bf[var]->phi[j_id];
                  d_func[0][MAX_VARIABLE_TYPES + w][j_id] -=
                     d_surface_flux[MASS_FRACTION][w] * phi_j / rho;
                  d_func[0][MAX_VARIABLE_TYPES + w][j_id] +=
                     surface_flux * d_rho->C[w][j_id] / (rho*rho);
                  for (kdir=0; kdir<dim; kdir++)
                     {
                      d_func[0][MAX_VARIABLE_TYPES + w][j_id] +=
                         d_vconv->C[kdir][w][j_id] * fv->snormal[kdir];
                     }
                 }
             }
         }
     }

  /* calculate the residual contribution        */
 
  *func = -surface_flux/rho;
 
  for (kdir=0; kdir<dim; kdir++)
     {
      *func += vconv[kdir] * fv->snormal[kdir];
     }
 
} /* END of routine vnorm_bc_electrodeposition  */

/******************************************************************************/
/******************************************************************************/

void 
lat_heat_bc(double func[],
	    double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
	    double x_dot[DIM],  /* mesh velocity                             */
	    dbl tt,		/* parameter to vary time integration from 
				 * explicit (tt = 1) to implicit (tt = 0)    */
	    dbl dt,		/* current value of the time step            */
	    int bc_input_id,
	    struct Boundary_Condition *BC_Types)

/******************************************************************************
*
*  Function which evaluates the latent heat boundary condition with leak 
*  caused by flux of species through surface n.(v-vs)=sum(evap vol flux)
*  It calculates the local contributions of the boundary condition to the
*  residual and jacobian
*
*      Adapted by D.A. Labreche (12/97) from kin_bc_leak by prs & rac
*
*  With ties to equilibrium mole fraction of species by calling Raoult's law
*  or Flory-Huggins VLE relations.
*			 Reviser: A. C. Sun (9/20/98)
*
*  Hopefully finished off with all generality by PRS (2/15/01)
* 
*  ----------------------------------------------------------------------------
*
*  Move all of the energy flux calculations to compute_leak_energy
*   by A.C.Sun (01/06)
*
*  Functions Called
*
*  compute_leak_energy - compute all of enthalpy contribution due to 
*                        mass flux changes across the interface
*
*  ----------------------------------------------------------------------------
*
*******************************************************************************/
     
{
  int j_id;
  int var;
  int w;
  double enorm; /*Calculated normal energy flux */
  NORMAL_ENERGY_DEPENDENCE_STRUCT d_enorm_struct;
  NORMAL_ENERGY_DEPENDENCE_STRUCT *d_enorm = &d_enorm_struct;
  struct Boundary_Condition *bc = NULL;
  struct Boundary_Condition *fluxbc = NULL;
  
  if(af->Assemble_LSA_Mass_Matrix)
    return;
  
  /* Calculate volume flux of bulk component through surface. viz. 
  *  this is the component that has no corresponding convective diffusion
  *  equation because it respresents the balance of the species accounted
  *  for by the overall mass balance.  Its evaporation is controlled through
  *  the input values on the latent_heat card directly, and is assumed to
  *  be that of a convective mass transfer model */
 
  
  bc = BC_Types + bc_input_id;
  if ( bc->BC_Data_Int[1] == -1 ) fluxbc = NULL;
 
  else fluxbc = BC_Types + bc->BC_Data_Int[1];
  
  compute_leak_energy(&enorm, d_enorm, tt, dt, bc, fluxbc);
  
  if (af->Assemble_Jacobian) {

 /* sum the contributions to the global stiffness matrix */
  
    var = TEMPERATURE;
        if (pd->v[pg->imtrx][var]) 
          {
            for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
              {
	        d_func[0][var][j_id] -= d_enorm->T[j_id];

              }
          }

    var = MASS_FRACTION;
    if (pd->v[pg->imtrx][var]) 
      {
	for (w = 0; w < pd->Num_Species_Eqn; w++) 
	  {
	    for ( j_id=0; j_id<ei[pg->imtrx]->dof[var]; j_id++)
	      {
		d_func[0][MAX_VARIABLE_TYPES + w][j_id] -= d_enorm->C[w][j_id];
	      }
	  }
      }
  }
  /* Calculate the residual contribution	*/
  
  *func = -enorm; 

} /* END of routine lat_heat_bc */
/******************************************************************************/
/*****************************************************************************/

void
lat_heat_internal_bc(double func[],
		     double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
		     const double tt,		/* parameter to vary time integration from 
						   explicit (tt = 1) to implicit (tt = 0)*/
		     const double dt,		/* current value of the time step */
		     const int bc_input_id,
		     const struct Boundary_Condition *BC_Types) 
     /******************************************************************************
      *
      *  Function which evaluates the latent heat boundary condition with leak 
      *  caused by convection in the absence of species transfer (lazy man's way).
      *  For example, phase transition via liquidus temperature.
      *  n.(v-vs) drives the enthalpy transfer weighted by latent_heat.
      *  The old routine was incorrect that it had called vconv then the quantity
      *  was then multiplied  by c[w].
      *  
      *  This is rewritten such that the species variable does not enter in the
      *  enthalpy calculations at all.  One can think of this as a lumped enthalpy
      *  contribution at the interface, whether or not it is an internal or external
      *  boundary.  If species are involved, one should be using latent_heat_bc
      *  to properly account for the component contribution and not this routine.
      *  
      *  The meaning of SOLID_LIQUID or LIQUID_VAPOR is lost as well due to
      *  this revision.  The input order is however preserved for backward compatibility.
      *
      *  It calculates the local contributions of the boundary condition to the
      *  residual and jacobian
      *
      *     Rewritten by A.C.Sun (01/06)
      *
      *  Functions Called
      *
      *  get_convection_velocity
      *
      *******************************************************************************/
{
  int j, j_id;
  int var,jvar,kdir;
  int dim;
  double phi_j;

  double vconv[MAX_PDIM]; /*Calculated convection velocity */
  double vconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;
  int err;

  /* local contributions of boundary condition to residual and jacobian */

  double  lhvflux, lat_heat_lump;

  if (af->Assemble_LSA_Mass_Matrix)
    {
      for (kdir = 0; kdir < pd->Num_Dim; kdir++)
	{
	  var = MESH_DISPLACEMENT1 + kdir;
	  if (pd->v[pg->imtrx][var])
	    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
	      {
		phi_j = bf[var]->phi[j];
		d_func[0][var][j] -= phi_j * fv->snormal[kdir];
	      }
	}
      return;
    }

  dim   = pd->Num_Dim;

/* Need to calculate lhvflux  and it's derivatives here */

/* call routine to calculate surface flux of this component and it's 
 * sensitivity to all variable types 
 */
  lhvflux = 0.;

/* get the convection velocity (it's different for arbitrary and
   lagrangian meshes) */

/* get deformation gradient, etc. if lagrangian mesh with inertia */

    if ((pd->MeshMotion == LAGRANGIAN ||
	 pd->MeshMotion == DYNAMIC_LAGRANGIAN) && pd->MeshInertia == 1)
    {
      err = belly_flop(elc->lame_mu);
      EH(err, "error in belly flop");
      if( neg_elem_volume ) return;
    }
    err = get_convection_velocity(vconv, vconv_old, d_vconv, dt, tt);
   EH(err, "Error in calculating effective convection velocity");


/* Residual calculation, note no contribution from species at ALL. */
  lat_heat_lump = BC_Types[bc_input_id].BC_Data_Float[0];

  for(kdir=0; kdir<dim; kdir++) 
    {    
      lhvflux +=  mp->density * lat_heat_lump * vconv[kdir]* fv->snormal[kdir];
    }

  if (af->Assemble_Jacobian) {

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
		  
                    d_func[0][var][j_id] -=  mp->density * lat_heat_lump*
			  (vconv[kdir]*fv->dsnormal_dx[kdir][jvar][j_id]
			  + d_vconv->X[kdir][jvar][j_id]*fv->snormal[kdir]);
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
		d_func[0][var][j_id] -=  mp->density * lat_heat_lump*
		  d_vconv->v[jvar][jvar][j_id] * fv->snormal[jvar];

	      }
	  }
      }

  }
  /* Calculate the residual contribution	*/
  
  *func = -lhvflux; 
  
} /* END of routine lat_heat_internal */
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
/* get_convection_velocity
 * This routine calculates the effective convection velocity to be used
 * in the species transport and energy transport equations.
 * If the mesh motion is arbitrary, then the convection velocity is
 * v - xdot.  If the mesh motion is Lagrangian (moving with solid), then
 * the convection velo is the sum of diffusion volume fluxes of 
 * non-solid species divided by the volume fraction of solids, plus an
 * advected Lagrangian solid velocity of the stress-free-state. 
 * 
 * Major change as of 5/17/01: This routine no longer adds an x_dot to the
 * porous media convective velocity, as we found errors in the equation.  Where
 * those equations need x_dot individually (as called from mm_fill_porous) we 
 * get it from another source. 
 */

int 
get_convection_velocity(double vconv[DIM], /*Calculated convection velocity */
			double vconv_old[DIM], /*Calculated convection velocity at previous time*/
                        CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv,
			double dt,
			double tt)
{
  int dim, p, q, b, var, i, w, w1;
  double volsolid, volsolid_old;
  
  dim = pd->Num_Dim;
  memset(vconv, 0, sizeof(double)*MAX_PDIM);
  if (af->Assemble_Jacobian && d_vconv != NULL ) {
    if (pd->v[pg->imtrx][MESH_DISPLACEMENT1])
      memset(d_vconv->X, 0, DIM*DIM*MDE*sizeof(dbl));

    if (pd->v[pg->imtrx][VELOCITY1] || pd->v[pg->imtrx][POR_LIQ_PRES])
      memset(d_vconv->v, 0, DIM*DIM*MDE*sizeof(dbl));

    if (pd->v[pg->imtrx][MASS_FRACTION] || pd->v[pg->imtrx][POR_LIQ_PRES])
      memset(d_vconv->C, 0, DIM*MAX_CONC*MDE*sizeof(dbl));

    if (pd->v[pg->imtrx][TEMPERATURE])
      memset(d_vconv->T, 0, DIM*MDE*sizeof(dbl));
  }


  if (cr->MeshMotion == ARBITRARY) {
    /* calculate v - xdot and it's derivatives */

    for (p = 0; p < VIM; p++) {
      vconv[p] = fv->v[p];
      if (pd->TimeIntegration != STEADY) {
	vconv_old[p] = fv_old->v[p];
	if (pd->v[pg->imtrx][R_MESH1]) { 
	  vconv[p]     -= fv_dot->x[p];
	  vconv_old[p] -= fv_dot_old->x[p];
	}
      } else {
        vconv_old[p] = 0.;
      }
    }
       
    if (af->Assemble_Jacobian && d_vconv != NULL) {
      for (p = 0; p < VIM; p++) {
	var = VELOCITY1 + p;
	if (pd->v[pg->imtrx][var]) {
	  for (i = 0; i < ei[pg->imtrx]->dof[var]; i++) {
	    d_vconv->v[p][p][i] = bf[var]->phi[i];
	  }
	}
	  
	if (pd->TimeIntegration != STEADY ) {
	  var = MESH_DISPLACEMENT1 + p;
	  if (pd->v[pg->imtrx][var]) {
	    for (i = 0; i < ei[pg->imtrx]->dof[var]; i++) {
	      d_vconv->X[p][p][i] -= bf[var]->phi[i] * (1 +2.* tt) /dt;
	    }
	  }
	}
      }
    }
  } /* end of if ARBITRARY */
   
   if ( cr->MeshMotion == LAGRANGIAN ||
	cr->MeshMotion == DYNAMIC_LAGRANGIAN || 
	cr->MeshMotion == TOTAL_ALE )
     {
       /*
	* In a lagrangian, dynamic_lagrangian, or total_ale mesh,
	* the convective velocity has two contributions
        * which are different in continuous and porous systems.  If the material
        * is continuous than the convective velocity is:
        *    v_conv  =  v_solvent_flux  +  v_solid_lagrangian_convection
	* PRS: I think this needs to be carefully checked, as the analogy is
	*      NOT what you want for porous materials. 
	*
        * If the material is porous, than the convective velocity is:
        *    v_conv  =   v_solid_lagrangian_convection
	* in some cases, and just
	*    v_conv  =  0. 
	*
	* In both of these cases, v_solid_lagrangian_convection is the velocity of
        * a moving solid in a fixed frame of reference where the velocity of the
        * material is know if it is unstressed - then v_solid_lagrangian_convection
	* is a mapping of the know velocity of the stress free state to the current
	* state
	*
	* The following sections add these contributions into the convective velocity
	*/


       /* calculate (sum diffusion flux)/solid volume frac. and it's derivatives 
	*  v_solvent_flux  */

       if (pd->v[pg->imtrx][MASS_FRACTION])   
	 /* if no solvents - then no convective velocity */
	 {
	   /* calculate volume fraction of solids and total diffusion flux */
        switch (mp->Species_Var_Type)
                {
                case SPECIES_DENSITY:
	                volsolid = 1.;
	                volsolid_old = 1.;
                        for(w=0 ; w<pd->Num_Species_Eqn ; w++)
	                   {
	                    volsolid -= MAX(fv->c[w],0)*mp->specific_volume[w];
	                    if ( pd->TimeIntegration != STEADY )
		             volsolid_old -= fv_old->c[w]*mp->specific_volume[w];
	                   }
                        break;
                case SPECIES_CONCENTRATION:
	                volsolid = 1.;
	                volsolid_old = 1.;
                        for(w=0 ; w<pd->Num_Species_Eqn ; w++)
	                   {
	                    volsolid -= fv->c[w]*mp->molar_volume[w];
	                    if ( pd->TimeIntegration != STEADY )
		             volsolid_old -= fv_old->c[w]*mp->molar_volume[w];
	                   }
                        break;
                default:
	                volsolid = 1.;
	                volsolid_old = 1.;
                        for(w=0 ; w<pd->Num_Species_Eqn ; w++)
	                   {
	                    volsolid -= fv->c[w];
	                    if ( pd->TimeIntegration != STEADY )
		                     volsolid_old -= fv_old->c[w];
	                   }
                        break;
                }
   if( volsolid <= 0)
        {
	double volso=1.0;
        fprintf(stderr,"nonvolatile fraction %g %g %g %g\n",
                        volsolid,fv->x[0],fv->x[1],fv->x[2]);
        for(w=0 ; w<pd->Num_Species_Eqn ; w++)
	     {
	      volso -= fv->c[w]*mp->specific_volume[w];
              fprintf(stderr,"spec %d %g %g %g\n",w,fv->c[w],
                                mp->specific_volume[w],volso);
             }
        WH(-1,"negative nonvolatile volume fraction");
        }

      for (p = 0; p < VIM; p++) {
	vconv[p] = 0.;
	vconv_old[p] = 0.;
	if ( cr->MassFluxModel == FICKIAN  ||
	     cr->MassFluxModel == STEFAN_MAXWELL  ||
	     cr->MassFluxModel == STEFAN_MAXWELL_CHARGED  ||
	     cr->MassFluxModel == STEFAN_MAXWELL_VOLUME ) /* Last modified; KSC: 9/98  and  RSL 6/29/00  */
	{
	  if ( Diffusivity() )  EH( -1, "Error in Diffusivity.");
	  /* Get diffusivity and derivatives */

	  for (w=0; w<pd->Num_Species_Eqn; w++)
	  {
	    vconv[p] -= mp->diffusivity[w] * fv->grad_c[w][p];
	    if ( pd->TimeIntegration != STEADY )
		vconv_old[p] -= mp->diffusivity[w] * fv_old->grad_c[w][p];
	  }
	}
	else if ( cr->MassFluxModel == GENERALIZED_FICKIAN)
	{
	  /* diffusion induced convection */
	  if ( Generalized_Diffusivity() )  EH( -1, "Error in Diffusivity.");
	  for (w=0; w<pd->Num_Species_Eqn; w++)
	  {
	    for (w1=0; w1<pd->Num_Species_Eqn; w1++)
	    {
	      vconv[p] -= mp->diffusivity_gen_fick[w][w1]
		  * fv->grad_c[w1][p];
	    }
	    if ( pd->TimeIntegration != STEADY ) {
	      for (w1=0; w1<pd->Num_Species_Eqn; w1++)
	      {
		vconv_old[p] -= mp->diffusivity_gen_fick[w][w1]
		    * fv_old->grad_c[w1][p];
	      }
	    }
	  }
	}
	else  if ( cr->MassFluxModel == DARCY )
	{ /* diffusion induced convection is zero */
	}
	else
	{
	  EH( -1, "Unimplemented mass flux constitutive relation.");
	}
      }
	   
      for (p=0; p<VIM; p++)
      {
	vconv[p] /= volsolid;
	if ( pd->TimeIntegration != STEADY )
	    vconv_old[p] /= volsolid_old;
      }
      if ( af->Assemble_Jacobian && d_vconv != NULL )
      {
	for (p=0; p<dim; p++)
	{
	  var = MESH_DISPLACEMENT1;
	  if (pd->v[pg->imtrx][var])
	  {
	    for (q=0; q<dim; q++)
	    {
	      for ( i=0; i<ei[pg->imtrx]->dof[var]; i++)
	      {
		if ( cr->MassFluxModel == FICKIAN )
		{
		  if ( Diffusivity() )  EH( -1, "Error in Diffusivity.");
				      
		  for (w=0; w<pd->Num_Species_Eqn; w++)
		  {
		    d_vconv->X[p][q][i] -= mp->diffusivity[w] *
			fv->d_grad_c_dmesh[p][w] [q][i] / volsolid;
		  }
		}
	      }
	    }
	  }
	  
	  var = MASS_FRACTION;
	  if (pd->v[pg->imtrx][var])
	  {
	    if ( cr->MassFluxModel == FICKIAN ||
		 cr->MassFluxModel == STEFAN_MAXWELL ||
		 cr->MassFluxModel == STEFAN_MAXWELL_CHARGED ||
		 cr->MassFluxModel == STEFAN_MAXWELL_VOLUME )
	    {
	      if ( Diffusivity() )  EH( -1, "Error in Diffusivity.");

	      for (w=0; w<pd->Num_Species_Eqn; w++)
	      {
		for ( i=0; i<ei[pg->imtrx]->dof[var]; i++)
		{
		  d_vconv->C[p][w][i] -=
		      mp->diffusivity[w] * bf[var]->grad_phi[i][p]
		      / volsolid;
		  for (w1=0; w1<pd->Num_Species_Eqn; w1++)
		  {
		    d_vconv->C[p][w][i] -= (
			mp->d_diffusivity[w][MAX_VARIABLE_TYPES + w1] 
			* fv->grad_c[w][p] 
			+ mp->diffusivity[w] *fv->grad_c[w][p] / volsolid )
			* bf[var]->phi[i] / volsolid;
		  }
		}
	      }
	    }
	  }
	} /* end of loop over vconv directions */
      } /* end of if Assemble Jacobian */
    } /* end of if MASS_FRACTION */
       
       /*
	* Add in convection due to motion of Stress Free State - Pseudo Lagrangian Convection 
	*  NOTE: this formulation still assumes mesh is quasi -   static, i.e. that mesh motion
	*        in transient calculations contributes negligibly to the momentum equation
        *    v_solid_lagrangian_convection
	*/

       /* First test to see what type of prescribed kinematics model */

    if (elc->v_mesh_sfs_model == CONSTANT)
    {
      /*do nothing??*/
    }
    else if (elc->v_mesh_sfs_model == ROTATIONAL ||
		elc->v_mesh_sfs_model == ROTATIONAL_3D )
    {
      (void) V_mesh_sfs_model(elc->u_v_mesh_sfs, elc->v_mesh_sfs, 
				elc->v_mesh_sfs_model, -1);
    }

       if ( pd->MeshInertia == 1)
	 {
	   if ( pd->TimeIntegration != STEADY )
	     WH(-1, "Can't have Unsteady Mesh Inertia. Maybe you don't want Convective Lagrangian Velocity? ");
	   /*
	    * Velocity of solid in lab coordinates is the velocity of the stress free state
	    * dotted into the deformation gradient tensor
	    */
	   if (! pd->v[pg->imtrx][MESH_DISPLACEMENT1])
	     {
	       for (p=0; p < dim; p++)
		 {
		   vconv[p] += elc->v_mesh_sfs[p];
		 }
	     }
	   else 
	     {
	       for (p=0; p < dim; p++)
		 {
		   for (q=0; q < dim; q++)
		     {
		       vconv[p] += elc->v_mesh_sfs[q] * fv->deform_grad[q][p];
		     }
		 }
	       if ( af->Assemble_Jacobian )
		 {
		   for (b=0; b < dim; b++)
		     {
		       var = MESH_DISPLACEMENT1 + b;
		       if (pd->v[pg->imtrx][var])
			 {
			   for ( i=0; i<ei[pg->imtrx]->dof[var]; i++)
			     {
			       for (p=0; p < dim; p++)
				 {
				   for (q=0; q < dim; q++)
				     {
				       d_vconv->X[p][b][i] += elc->v_mesh_sfs[q]
					 * fv->d_deform_grad_dx[q][p] [b][i];
				     }
				 }
			     }
			 }
		     }
		 }
	     }
	 }
       
     } /* end of if LAGRANGIAN, DYNAMIC_LAGRANGIAN or TOTAL_ALE */
   
  return 0;
} /* END of routine get_convection_velocity */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int
get_continuous_species_terms(struct Species_Conservation_Terms *st,
			     double time,             /* present time value; KSC: 10/98 */
			     double tt,               /* parameter to vary time integration
							 explicit (tt = 1), implicit (tt = 0) */
			     double dt,               /* current time step size */
			     const double hsquared[DIM])    /* element size measure */
     /**************************************************************************
      *   get_continuous_species_terms
      *
      *  This routine calculates source terms that are needed in species equations
      *  for continuous material types. It is called by assemble_mass_transport()
      *  at the quadrature point level once for all species degrees of freedom in
      *  the element.
      *
      * --------------------------------------------------------------------------
      *       Calculate the capacity, flux, convection, and source terms for a
      *       continuous medium
      * --------------------------------------------------------------------------
      *  Written by: Richard Cairncross  5/6/95 
      *  Modified by Ken S. Chen to enable diffusive-flux calculations using
      *                          Stefan-Maxwell model: 7/98   
      *
      **************************************************************************/
{
  double vconv[MAX_PDIM]; /*Calculated convection velocity */
  double vconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;

  /* MMH
   * For particle phase species.  Notice the pv's instead of the
   * v's everywhere.
   */
  double pvconv[MAX_PDIM]; /*Calculated convection velocity */
  double pvconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_pvconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_pvconv = &d_pvconv_struct;

  /* MMH
   * So we only have to set one ptr to take care of the separate species.
   * These replace the vconv[] stuff from before in the code below.
   */
  double *conv; /*Calculated convection velocity */
  double *conv_old; /*Calculated convection velocity at previous time*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_conv;

  dbl c_new = 0.0, c_old = 0.0;                 /* velocity averages */
  int w, w1, i, j, a, b, err, var, eqn, var_offset;
  int taylor_galerkin[MAX_CONC];
  int explicit[MAX_CONC]; /* move to input deck asap */

  int species;			/* Species number for the particle phase. */
  int wim   = pd->Num_Dim;    /* wim is the number of velocity unknowns */
  if(pd->CoordinateSystem == SWIRLING ||
     pd->CoordinateSystem == PROJECTED_CARTESIAN ||
     pd->CoordinateSystem == CARTESIAN_2pt5D)
    wim = wim+1;

  if (mp->DensityModel == SUSPENSION_PM)
    species = (int) mp->u_density[0];
  else
    species = -1;

  for ( w=0; w<pd->Num_Species_Eqn; w++)
    {
      taylor_galerkin[w] =mp->SpeciesTimeIntegration[w];
      if(mp->SpeciesTimeIntegration[w] == TAYLOR_GALERKIN_EXP)
	{
	  explicit[w]=1;
	}
      else
	{
	  explicit[w]=0;
	}
    }

  /*
   * CAPACITY TERM - concentration, volume fraction, mass fraction, etc.
   *                 (Calculate the derivative of the mass fraction
   *                  or volume fraction wrt time at the current gauss point)
   */      
  for ( w=0; w<pd->Num_Species_Eqn; w++) {
    st->Y[w]         = fv->c[w];
    st->Y_old[w]     = fv_old->c[w];
    st->Y_dot_old[w] = fv_dot_old->c[w];
    if (pd->TimeIntegration != STEADY) {
      st->Y_dot[w]    = fv_dot->c[w];
    } else {
      st->Y_dot[w]    = 0.0;
    }
  }
  for ( w=0; w<pd->Num_Species_Eqn; w++)
    {
      for ( a=0; a<VIM; a++)
	{
	  st->grad_Y[w][a] = fv->grad_c[w][a];
	  /* note, this is zero for steady calculations */
	  /* keep this only for VOF/Taylor-Galerkin stuff */
	  st->grad_Y_old[w][a] = fv_old->grad_c[w][a];
	}
    }
      
  /*
   * DIFFUSIVE FLUX TERM
   */      
  if ( cr->MassFluxModel == FICKIAN )
    {
      if ( Diffusivity() )  EH( -1, "Error in Diffusivity.");
	  
      for ( w=0; w<pd->Num_Species_Eqn; w++)
	{
	  err = fickian_flux(st,w);
	}
    }
  else if ( cr->MassFluxModel == GENERALIZED_FICKIAN )
    {
	  
      for ( w=0; w<pd->Num_Species_Eqn; w++)
	{
	  err = generalized_fickian_flux(st,w);
	}
    }
  else if ( cr->MassFluxModel == FICKIAN_CHARGED)     /* Fickian diffusion of charged species, KSC: 9/2000 */ 
    {
      if ( Diffusivity() )  EH( -1, "Error in Diffusivity.");
	  
      for ( w=0; w<pd->Num_Species_Eqn; w++)
	{
	  err = fickian_charged_flux(st,w);
	}
    }
  else if ( cr->MassFluxModel == FICKIAN_CHARGED_X )  /* RSL 9/18/00 */
     {
      err = fickian_charged_flux_x(st, time, dt);
     }
  else if ( cr->MassFluxModel == STEFAN_MAXWELL || 
            cr->MassFluxModel == STEFAN_MAXWELL_CHARGED ||
            cr->MassFluxModel == STEFAN_MAXWELL_VOLUME)  
    { 
      err = Stefan_Maxwell_diff_flux(st, time, dt);
    } 
  else if ( cr->MassFluxModel == HYDRODYNAMIC 
	    || cr->MassFluxModel == HYDRODYNAMIC_QTENSOR 
	    || cr->MassFluxModel == HYDRODYNAMIC_QTENSOR_OLD)
    {
      if ( Diffusivity() )  EH( -1, "Error in Diffusivity.");
      
      for ( w=0; w<pd->Num_Species_Eqn; w++)
	{
	  switch(mp->DiffusivityModel[w])
	    {
	    case CONSTANT:
	    case USER:
	    case POROUS:
	      fickian_flux(st, w);
	      break;
	      
	    case HYDRO:
	      hydro_flux(st, w, tt, dt, hsquared);
	      break;
	 
	    default:
	      EH( -1, "Unknown Diffusivity Model.");
	      break;

	    } /* end of switch(DiffusivityModel) */

	}
    }
  else if ( cr->MassFluxModel ==  DM_SUSPENSION_BALANCE )
    {
      if ( Diffusivity() )  EH( -1, "Error in Diffusivity.");
      
      for ( w=0; w<pd->Num_Species_Eqn; w++)
	{
	  switch(mp->DiffusivityModel[w])
	    {
	    case SUSP_BAL:
	      suspension_balance(st,w);
	      break;

	    default:
	      EH( -1, "Unknown Diffusivity Model for Suspension Balance.");
	      break;

	    } /* end of switch(DiffusivityModel) */

	}
    }
  else
    {
      EH( -1, "Unimplemented mass flux constitutive relation in continuous media.");
    }

  /* 
   * CONVECTIVE FLUX TERM
   */
  /* get the convection velocity (it's different for arbitrary and
     lagrangian meshes) */
  if ((pd->MeshMotion == LAGRANGIAN ||
       pd->MeshMotion == DYNAMIC_LAGRANGIAN) && pd->MeshInertia == 1)
    {
      err = belly_flop(elc->lame_mu);
      EH(err, "error in belly flop");
      if( neg_elem_volume ) return(err);
    }
  err = get_convection_velocity(vconv, vconv_old, d_vconv, dt, tt);
  EH(err, "Error in calculating effective convection velocity");

  /* MMH
   * If we are doing a SUSPENSION_PM model, get all the same information
   * for the particle velocity phase.  Use it for the proper species.
   */
  if(mp->DensityModel == SUSPENSION_PM)
    {
      err = get_particle_convection_velocity(pvconv, pvconv_old, d_pvconv, dt, tt);
      EH(err, "Error in calculating effective convection velocity for particle phase");
    }
      
  for ( w=0; w<pd->Num_Species_Eqn; w++)
    {
      if(mp->DensityModel == SUSPENSION_PM &&
	 w == species)
	{
	  conv = pvconv;
	  conv_old = pvconv_old;
          d_conv = d_pvconv;
	}
      else
	{
	  conv = vconv;
	  conv_old = vconv_old;
	  d_conv = d_vconv;
	}
      if(taylor_galerkin[w])
	{
	  if(tt == 0.0)
	    {
	      c_new=3./2.;
	      c_old=-1./2.;
	    }
	  else if(tt == 0.5)
	    {
	      c_new=1./2.;
	      c_old=1./2.;
	    }
	  eqn = R_MASS;
	  if(explicit[w])
	    {
	      for ( j=0; j<ei[pg->imtrx]->dof[eqn]; j++)
		{
		  st->taylor_flux_wt[j] = 0.;
		  for ( a=0; a<VIM; a++)
		    {
		      st->taylor_flux_wt[j]+= conv_old[a] * bf[eqn]->grad_phi[j] [a] ;
		    }
		}
	    }
	  else
	    {
	      for ( j=0; j<ei[pg->imtrx]->dof[eqn]; j++)
		{
		  st->taylor_flux_wt[j] = 0.;
		  for ( a=0; a<VIM; a++)
		    {
		      st->taylor_flux_wt[j]+= conv[a] * bf[eqn]->grad_phi[j] [a] ;
		    }
		}
	    }
	}
    }

  /*
   *  Loop over each species equation, entering the source term
   *  into a temporary vector
   */
  for ( w=0; w<pd->Num_Species_Eqn; w++)
    {
      /* MMH
       * Pick the correct set of convection velocities.
       * Fluid phase vs. Particle phase.
       */
      if( mp->DensityModel == SUSPENSION_PM &&
	  w == species )
	{
	  conv = pvconv;
	  conv_old = pvconv_old;
	  d_conv = d_pvconv;
	}
      else
	{
	  conv = vconv;
	  conv_old = vconv_old;
	  d_conv = d_vconv;
	}
      for (a = 0; a < VIM; a++)
	{
	  if (taylor_galerkin[w] && !(explicit[w]) )
	    {
	      st->conv_flux[w][a] = (c_new*conv[a] + c_old*conv_old[a]) * st->grad_Y[w][a];
	      st->taylor_flux[w][a] = conv[a] * st->grad_Y[w][a];
	    }
	  else if (taylor_galerkin[w] && explicit[w])
	    {
	      st->conv_flux[w][a] = (conv[a] + conv_old[a])/2. * st->grad_Y_old[w][a];
	      st->taylor_flux[w][a] = conv_old[a] * st->grad_Y_old[w][a];
	    }
	  else
	    {
	      st->conv_flux[w][a] = conv[a] * st->grad_Y[w][a];
	      st->taylor_flux[w][a] =0.;
	    }

	  if (mp->ExtrinsicIndependentSpeciesVar[w])
	    {
	      /* for compressible materials, add c(del dot v) portion of advection term */
	      st->conv_flux[w][a] += st->Y[w] * fv->grad_v[a][a];
	    }
	}
    }

  /*
   * SOURCE TERM
   */

  /*** Species Source ****/
  /*
   *  Do source terms that are best done with the species index
   *  as the inner loop.
   */
  if (mp->SpeciesSourceModel[0] == SSM_CHEMKIN_GAS) {
#ifdef USE_CHEMKIN
    /*
     *  Right now, for the material ID, we will use the material number.
     *  This really should be identified with the material number specified
     *  in the exodus data file. However, currently GOMA doesn't read
     *  that in .((.
     */
    /*
     *  Calculate the thermodynamic pressure in cgs units
     *  HKM -> Right now, we take the pressure from the Pressure Datum
     *         field in th upd structure. We do not consider the pressure field
     *         from the calculation itself, yet.
     */
    pressureCGS = upd->Pressure_Datum;
    if (af->Assemble_Jacobian) {
      jac_Species_Source = mp->Jac_Species_Source;
      d_species_source_T = mp->d_species_source;
      iopt[0] = 1;
      if (mp->Species_Var_Type == SPECIES_MASS_FRACTION) {
	/*
	 * Calculate the Jacobian d_source_i / d_MF_j under the
	 * conditions of constant other MF_l, l = 1, ..., Num_species
	 */
	err = ck_VD_dsdy(ei[pg->imtrx]->mn, &pressureCGS, &(fv->T), fv->c, NULL,
			 jac_Species_Source, num_species,
			 d_species_source_T, iopt, st->MassSource);
	/*
         * We need to fix the conditions up to account for Goma's
         * usage of n-1 equations. Thus, Goma requires the following
         * Jacobian entries,
         * d_source_i / d_MF_j under the
         * conditions of constant other MF_l, l = 1, ..., Num_species -1
         * and the additional constraing of sum(MF_l) = 1, where the
         * sum is over l = 1, ..., Num_Species.
         *
         * Also, let's convert the source terms and its Jacobian to g/cm**3sec
         * from mol/cm**3*sec.
         */
	if (mp->Dropped_Last_Species_Eqn) {
	  for (w = 0; w < pd->Num_Species_Eqn; w++) {
	    dtmp =  jac_Species_Source[w + pd->Num_Species * pd->Num_Species_Eqn];
	    mw = mp->molecular_weight[w];
	    for (j = 0; j < pd->Num_Species_Eqn; j++) {
	      jac_Species_Source[w + j * pd->Num_Species] -= dtmp;
	      jac_Species_Source[w + j * pd->Num_Species] *= mw;
	    }
	    st->MassSource[w] *= mw;
	  }
	} else {
	  for (w = 0; w < pd->Num_Species_Eqn; w++) {
	    mw = mp->molecular_weight[w];
	    for (j = 0; j < pd->Num_Species_Eqn; j++) {
	      jac_Species_Source[w + j * pd->Num_Species] *= mw;
	    }
	    st->MassSource[w] *= mw;
	  }
	}
      } else {
        EH(-1,"SPECIES_MASS_FRACTION only implementation so far");
      }

    } else {
    err = ck_VD_wyp(ei[pg->imtrx]->mn, &pressureCGS, &(fv->T), fv->c, NULL,
		    st->MassSource);
    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      st->MassSource[w] *= mp->molecular_weight[w];
    }
    }
    if (err != CPC_SUCCESS) {
      printf("failure\n");
      exit(-1);
    }
    if (af->Assemble_Jacobian) {
      var = TEMPERATURE;
      if (pd->v[pg->imtrx][var])   {
	for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
	  dtmp = bf[var]->phi[j];
	  for (w = 0; w < pd->Num_Species_Eqn; w++)  {
	    st->d_MassSource_dT[w][j] = d_species_source_T[w] * dtmp;
	  }
	}
      }
      if (pd->v[pg->imtrx][MASS_FRACTION] ) {
	for (w = 0; w < pd->Num_Species_Eqn; w++) {
	  var = MASS_FRACTION;
	  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
	    dtmp = bf[var]->phi[j];
	    for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++) {
	      st->d_MassSource_dc[w][w1][j] = jac_Species_Source[w + num_species*w1] * dtmp;
	    }
	  }
	}
      }
    }
#else
    chemkin_not_linked("SSM_CHEMKIN_GAS");
#endif
  } else if (mp->SpeciesSourceModel[0] == FOAM) {

    /* this is called in the interior nodes */

	err = foam_species_source( mp->u_species_source[0]);

    for ( w=0; w<pd->Num_Species_Eqn; w++)  {
	st->MassSource[w]= mp->species_source[w];

	  if ( af->Assemble_Jacobian ) {
	    var_offset = MAX_VARIABLE_TYPES + w;

	    var = TEMPERATURE;
	    if(pd->v[pg->imtrx][TEMPERATURE])
	      {
		for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		  {
		    st->d_MassSource_dT[w][j]=
		      mp->d_species_source[var_offset]*bf[var]->phi[j];
		  }
	      }

	    var = MASS_FRACTION;
	    if (pd->v[pg->imtrx][MASS_FRACTION] )
	      {
		for ( w1=0; w1<pd->Num_Species; w1++)
		  {
		    for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		      {
			st->d_MassSource_dc[w][w1][j]=
			  mp->Jac_Species_Source[w1+w*pd->Num_Species]*bf[var]->phi[j];
		      }
		  }
	      }
	  }
    }
  } else {
    /*
     *  Do source terms that are best done with the species index as the
     *  outer loop
     */
    for ( w=0; w<pd->Num_Species_Eqn; w++)  {
      if(mp->SpeciesSourceModel[w] == USER )
      {
	err = usr_species_source(w, mp->u_species_source[w]);
	st->MassSource[w]= mp->species_source[w];

	if ( af->Assemble_Jacobian )
	{
	  var = TEMPERATURE;
	  if(pd->v[pg->imtrx][var])
	  {
	    for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	    {
	      st->d_MassSource_dT[w][j]= mp->d_species_source[var]*bf[var]->phi[j];
	    }
	  }

	  var = VOLTAGE;
	  if(pd->v[pg->imtrx][var])
	  {
	    for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	    {
	      st->d_MassSource_dV[w][j]= mp->d_species_source[var]*bf[var]->phi[j];
	    }
	  }

	  if(pd->v[pg->imtrx][VELOCITY1])
	  {
	    for ( a=0; a<DIM; a++)
	    {
	      var = VELOCITY1 + a;
	      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	      {
		st->d_MassSource_dv[w][a][j]=mp->d_species_source[var]*bf[var]->phi[j];
	      }
	    }
	  }

	  if(pd->v[pg->imtrx][MESH_DISPLACEMENT1])
	  {
	    for ( a=0; a<DIM; a++)
	    {
	      var = MESH_DISPLACEMENT1 + a;
	      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	      {
		st->d_MassSource_dmesh[w][a][j] =mp->d_species_source[var]*bf[var]->phi[j];
	      }
	    }
	  }

	  if (pd->v[pg->imtrx][MASS_FRACTION] )
	  {
	    for ( w1=0; w1<pd->Num_Species_Eqn; w1++)
	    {
	      var = MASS_FRACTION;
	      var_offset = MAX_VARIABLE_TYPES + w1;
	      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	      {
		st->d_MassSource_dc[w][w1] [j]=mp->d_species_source[var_offset]*bf[var]->phi[j];
	      }
	    }
	  }

	}
      }
      else if(mp->SpeciesSourceModel[w] == PHOTO_CURING )
      {
       double intensity=0;
       double k_prop=0, k_term=0, k_inh=0, free_rad, d_free_rad_dI;
       double *param,s,dsdC[MAX_CONC],dsdT,dsdI,Conc[MAX_CONC];
       int model_bit, num_mon, O2_spec=-1, rad_spec=-1, init_spec = 0;
       double k_propX=1, k_propT=0, k_propX_num=0, k_propX_den=0;
       double intensity_cgs = 2.998e+10*8.85e-12/200.0;
       double dbl_small = 1.0e-15, Xconv_denom=0, sum_init=0;
       double Xconv=0.0, Xconv_init=0.0, sum_mon=0;
       double dXdC[MAX_CONC] = {0.0};

       param = mp->u_species_source[w];
       model_bit = ((int)param[0]);
       s = 0;       
       dsdT = 0;  dsdI = 0;  d_free_rad_dI = 0;
       for(a=0; a<MAX_CONC; a++) dsdC[a]=0.;  

       if(pd->e[pg->imtrx][R_LIGHT_INTP])
         {
         intensity = fv->poynt[0];
         if(pd->e[pg->imtrx][R_LIGHT_INTM])
          { intensity += fv->poynt[1];}
         if(pd->e[pg->imtrx][R_LIGHT_INTD])
          { intensity += fv->poynt[2];}
         intensity *= mp->u_species_source[init_spec][1];
         intensity = MAX(intensity,0.0);
         }   
       else if(pd->e[pg->imtrx][R_ACOUS_PREAL])
         {
         intensity = mp->u_species_source[init_spec][1]*
                     intensity_cgs*
                     (SQUARE(fv->apr)+SQUARE(fv->api));
         } 
       else
        { WH(-1,"No Intensity field found in PHOTO_CURING\n"); }

      /* insure concentrations are positive  */
	for ( w1=0; w1<pd->Num_Species_Eqn; w1++)
	    {
	     Conc[w1] = MAX(dbl_small,fv->c[w1]);
	    }

      /*  free radical concentration  */
	num_mon = model_bit>>2;
       if( (model_bit & 1)  && (model_bit & 2))
            {
	    O2_spec = init_spec + num_mon +2;
	    rad_spec = O2_spec + 1;
            free_rad = Conc[rad_spec];
            k_term = mp->u_species_source[rad_spec][1]*
                     exp(-mp->u_species_source[rad_spec][2]*
                     (1./fv->T - 1./mp->u_species_source[rad_spec][3]));
            k_inh = mp->u_species_source[O2_spec][1]*
                    exp(-mp->u_species_source[O2_spec][2]*
                    (1./fv->T - 1./mp->u_species_source[O2_spec][3]));
            } 
       else if( model_bit & 1)
            {
	    O2_spec = init_spec + num_mon +2;
            k_inh = mp->u_species_source[O2_spec][1]*
                    exp(-mp->u_species_source[O2_spec][2]*
                    (1./fv->T - 1./mp->u_species_source[O2_spec][3]));
            free_rad = sqrt(SQUARE(k_inh*Conc[O2_spec])/4.+
                mp->u_species_source[init_spec+1][2]*intensity*Conc[init_spec])
                - k_inh*Conc[O2_spec]/2.;
            d_free_rad_dI += 0.5*mp->u_species_source[init_spec+1][2]
                             *Conc[init_spec]
                  /sqrt(SQUARE(k_inh*Conc[O2_spec])/4.+
                mp->u_species_source[init_spec+1][2]*intensity*Conc[init_spec]);
            }
       else if( model_bit & 2)
            {
	    rad_spec = init_spec + num_mon +2;
            k_term = mp->u_species_source[rad_spec][1]*
                     exp(-mp->u_species_source[rad_spec][2]*
                     (1./fv->T - 1./mp->u_species_source[rad_spec][3]));
            free_rad = Conc[rad_spec];
            }
       else
            {
            free_rad = sqrt(mp->u_species_source[init_spec+1][2]*
                       intensity*Conc[init_spec]);
            if(free_rad > 0){
            d_free_rad_dI += mp->u_species_source[init_spec+1][2]*
                            Conc[init_spec]/free_rad;
                  }
            }

       switch(mp->Species_Var_Type)   {
           case SPECIES_DENSITY:
                for ( w1=init_spec+2; w1<init_spec+2+num_mon; w1++)
	            { 
                     Xconv += fv->external_field[w1]/mp->molecular_weight[w1]/mp->specific_volume[w1];
                     sum_mon += Conc[w1]/mp->molecular_weight[w1];
                     dXdC[w1] = -1.0/mp->molecular_weight[w1];
                     Xconv_init +=  mp->u_reference_concn[w1][1]/mp->molecular_weight[w1]/mp->specific_volume[w1];
                     sum_init += mp->u_reference_concn[w1][0]/mp->molecular_weight[w1];
                    }
                break;
           case SPECIES_CONCENTRATION:
                for ( w1=init_spec+2; w1<init_spec+2+num_mon; w1++)
	            { 
                     Xconv += fv->external_field[w1]/mp->specific_volume[w1];
                     sum_mon += Conc[w1];
                     dXdC[w1] = -1.0;
                     Xconv_init += mp->u_reference_concn[w1][1]/mp->specific_volume[w1];
                     sum_init += mp->u_reference_concn[w1][0];
                    }
                break;
           default:
                EH(-1,"invalid Species Type for PHOTO_CURING\n");
           }
       Xconv *= mp->specific_volume[pd->Num_Species_Eqn];
       Xconv_init *= mp->specific_volume[pd->Num_Species_Eqn];
       Xconv_denom = Xconv + sum_mon;
       Xconv /= Xconv_denom;
       Xconv = MAX(dbl_small,Xconv);
       Xconv = MIN(1.0-dbl_small,Xconv);
       Xconv_init /= (Xconv_init + sum_init);
       for ( w1=init_spec+2; w1<init_spec+2+num_mon; w1++)
            { dXdC[w1] *= Xconv/Xconv_denom; }
        if(Xconv <= dbl_small || Xconv >= (1.0-dbl_small) )
            { memset( dXdC, 0, sizeof(double) * MAX_CONC); }

       if(w == init_spec) 
         {
            s = -Conc[w]*intensity;
            dsdC[w] = -intensity;
            dsdI = -Conc[w];
         }
       else if(w == init_spec + 1)
         {
            s = param[2]*Conc[init_spec]*intensity;
            dsdC[init_spec] = param[2]*intensity;
            dsdI = param[2]*Conc[init_spec];
          }
       else if(w > init_spec+1 && w <= init_spec+num_mon+1)
         {
            k_prop = param[1]*exp(-param[2]*(1./fv->T - 1./param[3]));
            k_propX_num = (1.0-param[4])*(1.-Xconv)+param[4]*(1.0-Xconv_init);
            k_propX_den = k_propX_num - (1.0-param[4])*(1.-Xconv)*log((1.-Xconv)/(1.0-Xconv_init));
            k_propX = SQUARE(k_propX_num)/k_propX_den;
            k_propT = k_prop*k_propX;
            s = -k_propT*Conc[w]*free_rad;
            dsdC[w] = dXdC[w]*(param[4]-1.0)*k_propX*(2./k_propX_num
                          +log((1.-Xconv)/(1.0-Xconv_init))/k_propX_den);
            dsdC[w] *= k_prop*Conc[w]*free_rad;
            dsdC[w] += -k_propT*free_rad;
            dsdT = -Conc[w]*free_rad*k_propT*param[2]/SQUARE(fv->T);
            dsdI += -k_propT*Conc[w]*d_free_rad_dI;
	  if( (model_bit & 1)  && (model_bit & 2))
            {
            dsdC[rad_spec] = -k_propT*Conc[w];
            } 
          else if( model_bit & 1)
            {
            dsdC[O2_spec] = -k_propT*Conc[w]*(SQUARE(k_inh/2.)*Conc[O2_spec]/
                       sqrt(SQUARE(k_inh*Conc[O2_spec])/4.+
                mp->u_species_source[init_spec+1][2]*intensity*Conc[init_spec])
                       -k_inh/2.);
            dsdC[init_spec] = -k_propT*Conc[w]*0.5*
                       (mp->u_species_source[init_spec+1][2]*intensity/
                       sqrt(SQUARE(k_inh*Conc[O2_spec])/4.+
                mp->u_species_source[init_spec+1][2]*intensity*Conc[init_spec]));
            }
          else if( model_bit & 2)
            {
            dsdC[rad_spec] = -k_propT*Conc[w];
            }
          else
            {
            dsdC[init_spec] = -k_propT*Conc[w]*
                       (sqrt(mp->u_species_source[init_spec+1][2]*intensity/
                       Conc[init_spec])/2.);
            }
          }
       else if(w == O2_spec && (model_bit & 1))
         {
            s = -k_inh*Conc[w]*free_rad;
            dsdT = -Conc[w]*free_rad*k_inh*param[2]/SQUARE(fv->T);
            dsdI += -k_inh*Conc[w]*d_free_rad_dI;
	  if( model_bit & 2)
            {
            dsdC[w] = -k_inh*free_rad;
            dsdC[rad_spec] = -k_inh*Conc[w];
            } 
          else 
            {
            dsdC[w] = -k_inh*free_rad-k_inh*Conc[w]*(SQUARE(k_inh/2.)*Conc[w]/
                       sqrt(SQUARE(k_inh*Conc[w])/4.+
                mp->u_species_source[init_spec+1][2]*intensity*Conc[init_spec])
                       -k_inh/2.);
            dsdC[init_spec] = -k_inh*Conc[w]*0.5*intensity*
                       mp->u_species_source[init_spec+1][2]/
                       sqrt(SQUARE(k_inh*Conc[O2_spec])/4.+
           mp->u_species_source[init_spec+1][2]*intensity*Conc[init_spec]);
            }
          }
       else if(w == rad_spec && (model_bit & 2))
         {
	  if( model_bit & 1 )
            {
	    s = mp->u_species_source[init_spec+1][2]*Conc[init_spec]*intensity
                  -k_term*SQUARE(Conc[w])
                  -k_inh*Conc[O2_spec]*Conc[w];
            dsdC[init_spec] = mp->u_species_source[init_spec+1][2]*intensity;
            dsdC[O2_spec] = -k_inh*Conc[w];
            dsdC[w] = -k_term*2*Conc[w]-k_inh*Conc[O2_spec];
            dsdT = - SQUARE(Conc[w])*k_term*param[2]/SQUARE(fv->T)
                           -Conc[O2_spec]*Conc[w]*k_inh*
                            mp->u_species_source[O2_spec][2]/SQUARE(fv->T);
	    dsdI = mp->u_species_source[init_spec+1][2]*Conc[init_spec];
            } 
          else 
            {
	    s = mp->u_species_source[init_spec+1][2]*Conc[init_spec]*intensity
                  -k_term*SQUARE(Conc[w]);
            dsdC[init_spec] = mp->u_species_source[init_spec+1][2]*intensity;
            dsdC[rad_spec] = -k_term*free_rad*2.;
            dsdT = - SQUARE(Conc[w])*k_term*param[2]/SQUARE(fv->T);
	    dsdI = mp->u_species_source[init_spec+1][2]*Conc[init_spec];
            }
          }
        else
          {
            s = 0.;
          }

	st->MassSource[w]= s;

	if ( af->Assemble_Jacobian )
	{
	  var = TEMPERATURE;
	  if(pd->v[pg->imtrx][var])
	  {
	    for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	    {
	      st->d_MassSource_dT[w][j]= dsdT*bf[var]->phi[j];
	    }
	  }

	  if (pd->v[pg->imtrx][MASS_FRACTION] )
	  {
	    for ( w1=0; w1<pd->Num_Species_Eqn; w1++)
	    {
	      var = MASS_FRACTION;
	      var_offset = MAX_VARIABLE_TYPES + w1;
	      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	      {
		st->d_MassSource_dc[w][w1] [j]=dsdC[w1]*bf[var]->phi[j];
	      }
	    }
	  }
	  var = LIGHT_INTP;
	  if(pd->v[pg->imtrx][var])
	  {
	    for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	    {
	      st->d_MassSource_dI[w][j]= dsdI*bf[var]->phi[j]
                         *mp->u_species_source[init_spec][1];
	    }
	  }
	  var = LIGHT_INTM;
	  if(pd->v[pg->imtrx][var])
	  {
	    for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	    {
	      st->d_MassSource_dI[w][j]= dsdI*bf[var]->phi[j]
                         *mp->u_species_source[init_spec][1];
	    }
	  }
	  var = LIGHT_INTD;
	  if(pd->v[pg->imtrx][var])
	  {
	    for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	    {
	      st->d_MassSource_dI[w][j]= dsdI*bf[var]->phi[j]
                         *mp->u_species_source[init_spec][1];
	    }
	  }
	  var = ACOUS_PREAL;
	  if(pd->v[pg->imtrx][var])
	  {
	    for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	    {
	      st->d_MassSource_dI[w][j]= dsdI*bf[var]->phi[j]
                    *intensity_cgs*mp->u_species_source[init_spec][1]
                    *2.*fv->apr;
	    }
	  }
	  var = ACOUS_PIMAG;
	  if(pd->v[pg->imtrx][var])
	  {
	    for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	    {
	      st->d_MassSource_dI[w][j]= dsdI*bf[var]->phi[j]
                    *intensity_cgs*mp->u_species_source[init_spec][1]
                    *2.*fv->api;
	    }
	  }

	}
      }
      else if (mp->SpeciesSourceModel[w]  == EPOXY )
      {
	err = epoxy_species_source(w, mp->u_species_source[w]);
	st->MassSource[w]    =  mp->species_source[w];

	if ( af->Assemble_Jacobian )
	{
	  var = TEMPERATURE;
	  if(pd->v[pg->imtrx][var])
	  {
	    for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	    {
	      st->d_MassSource_dT[w][j]= mp->d_species_source[var]*bf[var]->phi[j];
	    }
	  }

	  var = MASS_FRACTION;
	  if (pd->v[pg->imtrx][MASS_FRACTION] )
	  {
	    var_offset = MAX_VARIABLE_TYPES + w;
	    for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	    {
	      st->d_MassSource_dc[w][w] [j]=mp->d_species_source[var_offset]
		  *bf[var]->phi[j];
	    }
	  }
	}

      }
      else if (mp->SpeciesSourceModel[w]  == FOAM_EPOXY)
      {
	err = foam_epoxy_species_source(w, mp->u_species_source[w], tt, dt);
	st->MassSource[w] =  mp->species_source[w];
	
	if ( af->Assemble_Jacobian )
	  {

	    var = MASS_FRACTION;
	    if (pd->v[pg->imtrx][MASS_FRACTION] )
	      {
                // Just include the diagonal contribution here
		var_offset = MAX_VARIABLE_TYPES + w;
		for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
		  {
		    st->d_MassSource_dc[w][w][j] = mp->d_species_source[var_offset]
		      *bf[var]->phi[j];
		  }
	      }

	    var = TEMPERATURE;
	    if (pd->v[pg->imtrx][var])
	      {
		for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
		  {
		    st->d_MassSource_dT[w][j] = mp->d_species_source[var] * bf[var]->phi[j];
		  }
	      }
	    
	  }
      }
      else if (mp->SpeciesSourceModel[w]  == EPOXY_DEA )
      {
	err = epoxy_dea_species_source(w, mp->u_species_source[w]);
	st->MassSource[w]    =  mp->species_source[w];
	if ( af->Assemble_Jacobian )
	{
	  var = TEMPERATURE;
	  if(pd->v[pg->imtrx][var])
	  {
	    for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	    {
	      st->d_MassSource_dT[w][j]= mp->d_species_source[var]*bf[var]->phi[j];
	    }
	  }

	  var = MASS_FRACTION;
	  if (pd->v[pg->imtrx][MASS_FRACTION] )
	  {
	    var_offset = MAX_VARIABLE_TYPES + w;
	    for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	    {
	      st->d_MassSource_dc[w][w] [j]=mp->d_species_source[var_offset]
		  *bf[var]->phi[j];
	    }
	  }
	}
      }
      else if (mp->SpeciesSourceModel[w]  == BUTLER_VOLMER)     /* added by KSC: 05/15/06 */
      {
        dbl dh[3], p[10];
        p[0] = w;
        for (j=1; j<10; j++)
         {
           p[j]= mp->u_species_source[w][j-1];
         }

        /* Computing current source by calling butler_volmer_source */
        st->MassSource[w] = butler_volmer_source(p, 2, dh);

        if ( af->Assemble_Jacobian )
         {
           var = TEMPERATURE;
           if(pd->v[pg->imtrx][var])
            {
              for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
               {
                 st->d_MassSource_dT[w][j]= dh[0]*bf[var]->phi[j];
               }
            }

           var = VOLTAGE;
           if(pd->v[pg->imtrx][var])
            {
              for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
               {
                 st->d_MassSource_dV[w][j]= dh[1]*bf[var]->phi[j];
               }
            }

           var = MASS_FRACTION;
           if (pd->v[pg->imtrx][MASS_FRACTION] )
            {
              for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
               {
                 st->d_MassSource_dc[w][w][j]=dh[2]*bf[var]->phi[j];
               }
            }
         }
      }
      else if (mp->SpeciesSourceModel[w]  == ELECTROOSMOTIC)     /* added by KSC: 06/09/06 */
      {
        dbl dh[3], n, nd;

        /* Computing species sink due to electro-osmotic drag  */
        n = mp->u_species_source[w][9];  // number of electrons involved
        nd = mp->u_species_source[w][10];  // electro-osmotic drag coef.
        st->MassSource[w] = n * nd * butler_volmer_source(mp->u_species_source[w], 2, dh);

        if ( af->Assemble_Jacobian )
         {
           var = TEMPERATURE;
           if(pd->v[pg->imtrx][var])
            {
              for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
               {
                 st->d_MassSource_dT[w][j] = n*nd*dh[0]*bf[var]->phi[j];
               }
            }

           var = VOLTAGE;
           if(pd->v[pg->imtrx][var])
            {
              for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
               {
                 st->d_MassSource_dV[w][j] = n*nd*dh[1]*bf[var]->phi[j];
               }
            }

           var = MASS_FRACTION;
           if (pd->v[pg->imtrx][MASS_FRACTION] )
            {
              for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
               {
                 st->d_MassSource_dc[w][w][j] = n*nd*dh[2]*bf[var]->phi[j];
               }
            }
         }
      }

      else if (mp->SpeciesSourceModel[w]  == ELECTRODE_KINETICS ) /* KSC: 10/98 */
      {
	electrode_species_source(w, time, dt);
	st->MassSource[w]    =  mp->species_source[w];
	if ( af->Assemble_Jacobian )
	{
	  var = TEMPERATURE;
	  if(pd->v[pg->imtrx][var])
	  {
	    for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	    {
	      st->d_MassSource_dT[w][j]= mp->d_species_source[var]*bf[var]->phi[j];
	      /*  above line activated by RSL 8/4/00; line below deactivated
	          st->d_MassSource_dT[w][j]= 0.0; */
	    }
	  }

	  var = VOLTAGE;  /*  RSL 8/4/00  */
	  if(pd->v[pg->imtrx][var])
	  {
	    for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	    {
	      st->d_MassSource_dV[w][j]= mp->d_species_source[var]*bf[var]->phi[j];
	    }
	  }

	  var = MASS_FRACTION;
	  if (pd->v[pg->imtrx][MASS_FRACTION] )
/* The following block generalized by RSL 8/8/00 */
	  {
	    for ( w1=0; w1<pd->Num_Species_Eqn; w1++)
	    {
	      var_offset = MAX_VARIABLE_TYPES + w1;
	      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	      {
		  st->d_MassSource_dc[w][w1][j]=mp->d_species_source[var_offset]*bf[var]->phi[j];
	      }
	    }
	  }
	}
      }

      else if (mp->SpeciesSourceModel[w] == FOAM_PBE_WATER) {
        foam_pbe_conversion_water(st, time, tt, dt);
      } else if (mp->SpeciesSourceModel[w] == FOAM_PBE_OH) {
        foam_pbe_conversion_OH(st, time, tt, dt);
      } else if (mp->SpeciesSourceModel[w] == FOAM_PBE_BA_G) {
        foam_pbe_ba_gas_source(st, time, tt, dt);
      } else if (mp->SpeciesSourceModel[w] == FOAM_PBE_BA_L) {
        foam_pbe_ba_liquid_source(st, time, tt, dt);
      } else if (mp->SpeciesSourceModel[w] == FOAM_PBE_CO2_G) {
        foam_pbe_co2_gas_source(st, time, tt, dt);
      } else if (mp->SpeciesSourceModel[w] == FOAM_PBE_CO2_L) {
        foam_pbe_co2_liquid_source(st, time, tt, dt);
      }
      else if (mp->SpeciesSourceModel[w]  == FOAM_PMDI_10_RXN)
      {
        err = foam_pmdi10_rxn_species_source(w, mp->u_species_source[w], tt, dt);
        st->MassSource[w] =  mp->species_source[w];

        if ( af->Assemble_Jacobian )
        {

          var = MASS_FRACTION;
          if (pd->v[pg->imtrx][MASS_FRACTION] )
          {
            // Just include the diagonal contribution here
            var_offset = MAX_VARIABLE_TYPES + w;
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
            {
              st->d_MassSource_dc[w][w][j] = mp->d_species_source[var_offset]
                                             *bf[var]->phi[j];
            }
          }

          var = TEMPERATURE;
          if (pd->v[pg->imtrx][var])
          {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
            {
              st->d_MassSource_dT[w][j] = mp->d_species_source[var] * bf[var]->phi[j];
            }
          }
        }
      }
      else if (mp->SpeciesSourceModel[w]  == FOAM_PMDI_10_H2O)
      {
        for ( w1=0; w1<pd->Num_Species_Eqn; w1++)
        {
          var_offset = MAX_VARIABLE_TYPES + w1;
          mp->d_species_source[var_offset] = 0.0;
        }

        err = foam_pmdi10_h2o_species_source(w, mp->u_species_source[w], time, tt, dt);
        st->MassSource[w] =  mp->species_source[w];

        if ( af->Assemble_Jacobian )
        {

          var = MASS_FRACTION;
          if (pd->v[pg->imtrx][MASS_FRACTION] )
          {
            for ( w1=0; w1<pd->Num_Species_Eqn; w1++)
            {
              var_offset = MAX_VARIABLE_TYPES + w1;
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
              {
                st->d_MassSource_dc[w][w1][j] = mp->d_species_source[var_offset]
                                                *bf[var]->phi[j];
              }
            }
          }

          var = TEMPERATURE;
          if (pd->v[pg->imtrx][var])
          {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
            {
              st->d_MassSource_dT[w][j] = mp->d_species_source[var] * bf[var]->phi[j];
            }
          }
        }
      }
      else if (mp->SpeciesSourceModel[w]  == FOAM_PMDI_10_CO2)
      {
        err = foam_pmdi10_co2_species_source(w, mp->u_species_source[w], time, tt, dt);
        st->MassSource[w] =  mp->species_source[w];

        if ( af->Assemble_Jacobian )
        {

          var = MASS_FRACTION;
          if (pd->v[pg->imtrx][MASS_FRACTION] )
          {
            for ( w1=0; w1<pd->Num_Species_Eqn; w1++)
            {
              var_offset = MAX_VARIABLE_TYPES + w1;
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
              {
                st->d_MassSource_dc[w][w1][j] = mp->d_species_source[var_offset]
                                                *bf[var]->phi[j];
              }
            }
          }

          var = TEMPERATURE;
          if (pd->v[pg->imtrx][var])
          {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
            {
              st->d_MassSource_dT[w][j] = mp->d_species_source[var] * bf[var]->phi[j];
            }
          }
        }
      }
      else if (mp->SpeciesSourceModel[w]  == FOAM_PMDI_10_CO2_LIQ)
      {
        err = foam_pmdi10_co2_liq_species_source(w, st, mp->u_species_source[w], time, tt, dt);

        st->MassSource[w] =  mp->species_source[w];

        if ( af->Assemble_Jacobian )
        {

          var = MASS_FRACTION;
          if (pd->v[pg->imtrx][MASS_FRACTION] )
          {
            for ( w1=0; w1<pd->Num_Species_Eqn; w1++)
            {
              var_offset = MAX_VARIABLE_TYPES + w1;
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
              {
                st->d_MassSource_dc[w][w1][j] = mp->d_species_source[var_offset]
                                                *bf[var]->phi[j];
              }
            }
          }

          var = TEMPERATURE;
          if (pd->v[pg->imtrx][var])
          {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
            {
              st->d_MassSource_dT[w][j] = mp->d_species_source[var] * bf[var]->phi[j];
            }
          }
        }
      }
      else if (mp->SpeciesSourceModel[w]  == FOAM_PMDI_10_CO2_GAS)
      {
        err = foam_pmdi10_co2_gas_species_source(w, st, mp->u_species_source[w], time, tt, dt);

        st->MassSource[w] =  mp->species_source[w];

        if ( af->Assemble_Jacobian )
        {

          var = MASS_FRACTION;
          if (pd->v[pg->imtrx][MASS_FRACTION] )
          {
            for ( w1=0; w1<pd->Num_Species_Eqn; w1++)
            {
              var_offset = MAX_VARIABLE_TYPES + w1;
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
              {
                st->d_MassSource_dc[w][w1][j] = mp->d_species_source[var_offset]
                                                *bf[var]->phi[j];
              }
            }
          }

          var = TEMPERATURE;
          if (pd->v[pg->imtrx][var])
          {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
            {
              st->d_MassSource_dT[w][j] = mp->d_species_source[var] * bf[var]->phi[j];
            }
          }
        }
      }


      else if (mp->SpeciesSourceModel[w]  == CONSTANT )
      {
	st->MassSource[w]    = mp->species_source[w];
	/* Jacobians should default to zero without us doing anything b/c of
	 *  the initial memset of st - hopefully!
	 */
      }

      else if (mp->SpeciesSourceModel[w]  == ION_REACTIONS) /*  RSL 6/6/02  */
         ; /*  this just prevents an error message from being printed  */

      /*
       *  HKM -> Chemkin Homogeneneous kinetics goes here
       *               SSM_CHEMKIN_GAS  -> Gas phase package call
       *               SSM_CHEMKIN_LIQ  -> Liquid phase package call
       *               SSM_CHEMKIN_CPC  -> Condensed phase package call
       */
      else if (mp->SpeciesSourceModel[w] == SSM_CHEMKIN_GAS) {
	EH(-1,"Chemkin Gas source model is implemented for some species and not others");
      }
      else if (mp->SpeciesSourceModel[w] == SSM_CHEMKIN_LIQ) {
	EH(-1,"Chem Liq source model is implemented for some species and not others");
      }
      else if (mp->SpeciesSourceModel[w] == SSM_CHEMKIN_CPC) {
	EH(-1,"CPC source model is implemented for some species and not others");
      } else {
	printf("species source model = %d\n", mp->SpeciesSourceModel[w]);
	EH(-1,"Unrecognized species source model");
      }

      if( ls != NULL ) ls_modulate_speciessource ( w,  st );

    }

    for (w=0; w<pd->Num_Species; w++) /*  RSL 3/20/01 and 6/6/02 -- need all species  */
       {
        if (mp->SpeciesSourceModel[w]  == ION_REACTIONS)
           {
            ion_reaction_source(w);
            st->MassSource[w] = mp->species_source[w];
            if (af->Assemble_Jacobian)
               {
                var = TEMPERATURE;
                if (pd->v[pg->imtrx][var])
                   {
                    for (j=0; j<ei[pg->imtrx]->dof[var]; j++)
                       {
                        st->d_MassSource_dT[w][j] = mp->d_species_source[var]*bf[var]->phi[j];
                       }
                    }

                var = MASS_FRACTION;
                if (pd->v[pg->imtrx][var])
                   {
                    for (w1=0; w1<pd->Num_Species_Eqn; w1++)
                       {
                        var_offset = MAX_VARIABLE_TYPES + w1;
                        for (j=0; j<ei[pg->imtrx]->dof[var]; j++)
                           {
                            st->d_MassSource_dc[w][w1][j] = mp->d_species_source[var_offset]*
                                                            bf[var]->phi[j];
                           }
                       }
                   }
               }
               if( ls != NULL ) ls_modulate_speciessource ( w,  st );
           }
       }
  }
    

    /*
   * NOW, CALCULATE SENSITIVITIES for the Jacobian, if needed
   */
  if ( af->Assemble_Jacobian )
    {
      /*
       * sensitivity of CAPACITY TERM - concentration (volume fraction)
       */
      var = MASS_FRACTION;
      if (pd->TimeIntegration != STEADY) {
	for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
	  for (w = 0; w < pd->Num_Species_Eqn; w++) {
	    st->d_Y_dot_dc[w][w] [j]= (1 + 2. * tt) * bf[var]->phi[j]/dt;
	  }
	}
      }

      /*
       * CONVECTIVE FLUX TERM
       */
      for ( w=0; w<pd->Num_Species_Eqn; w++)
	{
          var = MASS_FRACTION;
	  /* MMH
	   * Pick the correct set of convection velocities.
	   * Fluid phase vs. Particle phase.
	   */
	  if(mp->DensityModel == SUSPENSION_PM &&
	     w == species)
	    {
	      conv = pvconv;
	      conv_old = pvconv_old;
	      d_conv = d_pvconv;
	    }
	  else
	    {
	      conv = vconv;
	      conv_old = vconv_old;
	      d_conv = d_vconv;
	    }
	  for ( a=0; a<VIM; a++)
	    {
	      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		{
		  if(taylor_galerkin[w] && !(explicit[w]))
		    {
		      for ( w1=0; w1<pd->Num_Species_Eqn; w1++)
			{
			  st->d_conv_flux_dc[w][a] [w1][j] =  c_new*d_conv->C[a][w1][j] * st->grad_Y[w][a];
			  st->d_taylor_flux_dc[w][a] [w1][j] =  d_conv->C[a][w1][j] * st->grad_Y[w][a];
			}
		      st->d_conv_flux_dc[w][a] [w][j] += (c_new*conv[a] + c_old*conv_old[a]) * bf[var]->grad_phi[j][a];
		      st->d_taylor_flux_dc[w][a] [w][j] += conv[a] * bf[var]->grad_phi[j][a];
		    }
		  else if(taylor_galerkin[w] && explicit[w])
		    {
		      for ( w1=0; w1<pd->Num_Species_Eqn; w1++)
			{
			  st->d_conv_flux_dc[w][a] [w1][j] =  0.;
			  st->d_taylor_flux_dc[w][a] [w1][j] =  0.;
			}
		    }
		  else
		    {
		      for ( w1=0; w1<pd->Num_Species_Eqn; w1++)
			{
			  st->d_conv_flux_dc[w][a] [w1][j] =  d_conv->C[a][w1][j] * st->grad_Y[w][a];
			  st->d_taylor_flux_dc[w][a] [w1][j] =  0.;
			}
		      st->d_conv_flux_dc[w][a] [w][j] += conv[a] * bf[var]->grad_phi[j][a];
		    }
	 
		  if (mp->ExtrinsicIndependentSpeciesVar[w]) 
		    {
		      /* for compressible materials, add c(div_v) portion of advection term */
		      
		      st->d_conv_flux_dc[w][a] [w][j] += bf[var]->phi[j] * fv->grad_v[a][a];
		    }
		}
	    }
	  eqn = R_MASS;
	  for (b=0; b < pd->Num_Dim; b++)
	    {
	      var = MESH_DISPLACEMENT1 + b;
	      if ( pd->v[pg->imtrx][var] )
		{
		  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		    {
		      if(taylor_galerkin[w] && !(explicit[w]))
			{
			  for (a=0; a < VIM; a++)
			    {
			      st->d_conv_flux_dmesh[w][a] [b][j] =  (c_new*conv[a] + c_old*conv_old[a])
				* fv->d_grad_c_dmesh[a][w] [b][j]
				+ c_new*d_conv->X[a][b][j] * st->grad_Y[w][a];
			      st->d_taylor_flux_dmesh[w][a] [b][j] = conv[a]
				* fv->d_grad_c_dmesh[a][w] [b][j]
				+ d_conv->X[a][b][j] * st->grad_Y[w][a];
			    }
			  for ( i=0; i<ei[pg->imtrx]->dof[eqn]; i++)
			    {
			      st->d_taylor_flux_wt_dmesh[i] [b][j] = 0.;
			      for (a=0; a < VIM; a++)
				{
				  st->d_taylor_flux_wt_dmesh[i] [b][j]+= conv[a]
				    * bf[eqn]->d_grad_phi_dmesh[i][a] [b][j]
				    + d_conv->X[a][b][j] * st->grad_Y[w][a];
				}
			    }
			}
		      else if(taylor_galerkin[w] && explicit[w])
			{
			  for (a=0; a < VIM; a++)
			    {
			      st->d_conv_flux_dmesh[w][a] [b][j] =  0.;
			      st->d_taylor_flux_dmesh[w][a] [b][j] = 0.;
			    }
			  for ( i=0; i<ei[pg->imtrx]->dof[eqn]; i++)
			    {
			      st->d_taylor_flux_wt_dmesh[i] [b][j] = 0.;
			    }
			}
		      else
			{
			  for (a=0; a < VIM; a++)
			    {
			      st->d_conv_flux_dmesh[w][a] [b][j] = conv[a]
				* fv->d_grad_c_dmesh[a][w] [b][j]
				+ d_conv->X[a][b][j] * st->grad_Y[w][a];
			      st->d_taylor_flux_dmesh[w][a] [b][j] = 0.;

			    }
			  for ( i=0; i<ei[pg->imtrx]->dof[eqn]; i++)
			    {
			      st->d_taylor_flux_wt_dmesh[i] [b][j] = 0.;
			    }
			}

		      if (mp->ExtrinsicIndependentSpeciesVar[w]) 
			{
			  /* for compressible materials, add c(div_v) portion of advection term */
			  for (a = 0; a < VIM; a++)
			    { 
			      st->d_conv_flux_dmesh[w][a] [b][j] += st->Y[w] * fv->d_grad_v_dmesh[a][a][b][j];
			    }
			}
		    }
		}
	    }

	  for (a = 0; a < wim; a++)
	    {
	      var = VELOCITY1 + a;
	      if ( pd->v[pg->imtrx][var] )
		{
		  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		    {
		      if(taylor_galerkin[w] && !(explicit[w]))
			{
			  for (b=0; b < pd->Num_Dim; b++)
			    {
			      st->d_conv_flux_dv[w][b] [a][j] = c_new*d_conv->v[b][a][j]
				* st->grad_Y[w][b];
			      st->d_taylor_flux_dv[w][b] [a][j] = d_conv->v[b][a][j]
				* st->grad_Y[w][b];
			    }

			  for ( i=0; i<ei[pg->imtrx]->dof[var]; i++)
			    {
			      st->d_taylor_flux_wt_dv[i][a][j] = 0.;
			      for (b=0; b < pd->Num_Dim; b++)
				{
				  st->d_taylor_flux_wt_dv[i][a][j] += d_conv->v[b][a][j]  * bf[eqn]->grad_phi[i] [b] ;
				}
			    }
			}
		      else if(taylor_galerkin[w] && explicit[w])
			{
			  for (b=0; b < pd->Num_Dim; b++)
			    {
			      st->d_conv_flux_dv[w][b] [a][j] = 0.;
			      st->d_taylor_flux_dv[w][b] [a][j] = 0.;
			    }
			  for ( i=0; i<ei[pg->imtrx]->dof[var]; i++)
			    {
			      st->d_taylor_flux_wt_dv[i][a][j] = 0.;
			    }
			}
		      else
			{
			  for (b=0; b < pd->Num_Dim; b++)
			    {
			      st->d_conv_flux_dv[w][b] [a][j] = d_conv->v[b][a][j]
				* st->grad_Y[w][b];
			      st->d_taylor_flux_dv[w][b] [a][j] = 0.;

			    }
			  for ( i=0; i<ei[pg->imtrx]->dof[var]; i++)
			    {
			      st->d_taylor_flux_wt_dv[i][a][j] = 0.;
			    }
			}

		      if (mp->ExtrinsicIndependentSpeciesVar[w]) 
			{
			  /* for compressible materials, add c(div_v) portion of advection term */
			  for (b = 0; b < VIM; b++)
			    {
			      st->d_conv_flux_dv[w][b][a][j] += st->Y[w] * bf[VELOCITY1+a]->grad_phi_e[j][a][b][b];
			    }
			}
		    }
		}
	    }

	  var = TEMPERATURE;
	  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	    {
	      if(taylor_galerkin[w] && !(explicit[w]))
		{
		  for ( a=0; a<VIM; a++)
		    {
		      st->d_conv_flux_dT[w][a] [j] = c_new*d_conv->T[a][j] * st->grad_Y[w][a];
		      st->d_taylor_flux_dT[w][a] [j] = d_conv->T[a][j] * st->grad_Y[w][a];

		    }
		  for ( i=0; i<ei[pg->imtrx]->dof[eqn]; i++)
		    {
		      st->d_taylor_flux_wt_dT[i] [j] = 0.;
		      for (a=0; a < VIM; a++)
			{
			  st->d_taylor_flux_wt_dT[i] [j]+= d_conv->T[a][j]* bf[eqn]->grad_phi[j] [a] ;
			}
		    }
		}
	      else if(taylor_galerkin[w] && explicit[w])
		{
		  for ( a=0; a<VIM; a++)
		    {
		      st->d_conv_flux_dT[w][a] [j] = 0.;
		      st->d_taylor_flux_dT[w][a] [j] = 0.;

		    }
		  for ( i=0; i<ei[pg->imtrx]->dof[eqn]; i++)
		    {
		      st->d_taylor_flux_wt_dT[i] [j] = 0.;
		    }
		}
	      else
		{
		  for ( a=0; a<VIM; a++)
		    {
		      st->d_conv_flux_dT[w][a] [j] = d_conv->T[a][j] * st->grad_Y[w][a]; /* UMR FLAGGED on tbc */
		      st->d_taylor_flux_dT[w][a] [j] = 0.;
		    }
		  for ( i=0; i<ei[pg->imtrx]->dof[eqn]; i++)
		    {
		      st->d_taylor_flux_wt_dT[i] [j] = 0.;
		    }
		}
	    }
	}
    } /* end of if Jacobian */
  return 0;
}
/******************************************************************************
*
*  Function that computes the Residuals and Jacobian sensitivities
*  for the Stefan-Maxwell diffusion terms and loads them into st
*
*  Author: Ken S. Chen, 8/98
*  Revised: R. S. Larson, 4/00, mostly to supply missing Jacobian terms
*
******************************************************************************/
int
Stefan_Maxwell_diff_flux( struct Species_Conservation_Terms *st,
			  double time,
			  double dt )
{
  int w1, a, j, k, l, q, var;
  int w, wi, wj;
  int status = 0;
  dbl dB[MAX_CONC*DIM];           /* dB/dx vector where x is mole fraction   */
  dbl dA[MAX_CONC*DIM][MAX_CONC*DIM]; /* dA/dx matrix                        */
  dbl dAJ[MAX_CONC*DIM];          /* dA/dx vector multiplied by the J vector */
  dbl scratch[MAX_CONC*DIM][MAX_CONC*DIM];

  const double R = 8.314;         /* Universal gas constant ( J/mole K ) */
  const double F = 96487.0;       /* Faraday's constant  ( C/equiv ) */
  dbl T=298.0;                    /* Temperature; set default value to room temperature */
  dbl rho;                        /* density */
  dbl c;                          /* total molar concentration */
  dbl e=1.0;                      /* porosity; set default value to unity */ 
  dbl x[MAX_CONC];                /* Mole Fraction */ 
  dbl M[MAX_CONC], M_mix;         /* molecular weight of individual species */
                                  /* and of mixture, respectively */ 
  dbl z[MAX_CONC];                /* charge number */
  dbl J[MAX_CONC*DIM];            /* Stefan_Maxwell flux vector */ 
  dbl dJ[MAX_CONC*DIM];            
  dbl temp_var;                   
  int volume_flag;               
  dbl correction;                

  static unsigned int A_allocated=FALSE; /* Boolean - allocation flag         */
  static double **A;              /* Coef. Matrix in Stefan_Maxwell Matrix
				     Equation                                 */
  static double **A_inv;          /*  inverse of Coef. Matrix                 */
  static double  *B;              /* driving force vector in Stefan-Maxwell
				     Equation                                 */
  static double  *bbb;            /*    reused for repeated                   */
  static    int  *indx;           /*  scaling array used in lu decomposition  */
  
  dbl D[MAX_CONC][MAX_CONC];            /* Stefan_Maxwell diffusivities       */
  dbl grad_mu[MAX_CONC][DIM];           /* electrochemical potential gradient */
  dbl grad_phi2[DIM];                   /* gradient of electrical potential   */
                                        /*       in electrolyte               */
  dbl grad_x[MAX_CONC][DIM];            /* mole fraction gradient */

  /*  double cc[MAX_CONC*DIM][MAX_CONC*DIM];  cc is a dummy matrix */
  int n_species, n;

  int mn, m, i, i1, i2, i3, j1, j2, j3;
  dbl sumx, sumdelx; 
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;  /* density dependence */
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;

  dbl T0, EE, alpha;  /* KSC on 9/24/04 */

  if (MAX_CONC < 3) {
    EH(-1, "Stefan_Maxwell_Diff_flux expects MAX_CONC >= 3");
    return -1;
  }

  var = MASS_FRACTION;
  n_species = pd->Num_Species_Eqn + 1; 
  n = n_species*VIM;
  mn = ei[pg->imtrx]->mn;

  if (!A_allocated) {
    A     = (double **) smalloc( MAX_CONC *DIM * sizeof(double *) );
    A_inv = (double **) smalloc( MAX_CONC *DIM * sizeof(double *) );
    bbb   =  (double *) smalloc( MAX_CONC *DIM * sizeof(double) );
    B     =  (double *) smalloc( MAX_CONC *DIM * sizeof(double) );
    indx  =     (int *) smalloc( MAX_CONC *DIM * sizeof(int) );
    for (i = 0; i < MAX_CONC * DIM; i++ )  {
      A[i]     = (double *) smalloc( MAX_CONC *DIM * sizeof( double ) );
      A_inv[i] = (double *) smalloc( MAX_CONC *DIM * sizeof( double ) );
    }
    A_allocated = TRUE;
  }

  for (i = 0; i < (MAX_CONC*DIM); i++) {
    memset(    A[i], 0, sizeof(double) * MAX_CONC * DIM); 
    memset(A_inv[i], 0, sizeof(double) * MAX_CONC * DIM); 
  }
  memset(indx, 0, sizeof(int)    * MAX_CONC * DIM); 
  memset(   B, 0, sizeof(double) * MAX_CONC * DIM); 
  memset( bbb, 0, sizeof(double) * MAX_CONC * DIM); 

  /*
   * Put in a warning about these species variable types until
   * the equations can be worked out
   */
  if (mp->Species_Var_Type == SPECIES_MASS_FRACTION ||
      mp->Species_Var_Type == SPECIES_MOLE_FRACTION ||
      mp->Species_Var_Type == SPECIES_VOL_FRACTION) {
    EH(-1,"Possible Conflict: Stefan Flux Expression hasn't been checked out for this species var type");
  }

  volume_flag = (pd_glob[mn]->MassFluxModel == STEFAN_MAXWELL_VOLUME); /* RSL 8/24/00 */

  if(pd_glob[mn]->e[R_ENERGY])  /* if the energy-transport equation is active, get temperature from fv */ 
    {
      T = fv->T;   
      if(pd_glob[mn]->MassFluxModel == STEFAN_MAXWELL_CHARGED)
        {
          if(mp->SolutionTemperatureModel == CONSTANT)
            {
              T = mp->solution_temperature;
            }
          else if(mp->SolutionTemperatureModel == THERMAL_BATTERY)
            {              /* need to get T from electrolyte_temperature() for the case of thermal battery */
              electrolyte_temperature(time, dt, 0);  /* calculate electrolyte temperature at present time */
              T = mp->electrolyte_temperature;
            }
          else
            {
              EH(-1, "Solution-temperature model other than THERMAL_BATTERY awaits future implementation");
            }
        }
      if(T == 0.0) T = 298.0;  /* set the solution temperature to the 298 K if it is zero - safety feature */

    }
  else                          /* when the energy-transport equation is NOT being solved */ 
    { 
      if(pd_glob[mn]->MassFluxModel == STEFAN_MAXWELL)
        { 
          if(mp->SolutionTemperatureModel == CONSTANT)
            { 
              T = mp->solution_temperature;
            }
          else
            {
              EH(-1, "User needs to define non-CONSTANT solution temperature model in GOMA");
            }
        }
      else if(pd_glob[mn]->MassFluxModel == STEFAN_MAXWELL_CHARGED || volume_flag)  /*  RSL 8/24/00  */
        {
          if(mp->SolutionTemperatureModel == CONSTANT)    
            {              
              T = mp->solution_temperature;
            }
          else if(mp->SolutionTemperatureModel == THERMAL_BATTERY)    
            {              /* need to get T from electrolyte_temperature() for the case of thermal battery */
              electrolyte_temperature(time, dt, 0);  /* calculate electrolyte temperature at present time */
              T = mp->electrolyte_temperature;
            }
          else
            {
              EH(-1, "Solution-temperature model other than THERMAL_BATTERY awaits future implementation"); 
            } 
        }
      if(T == 0.0) T = 298.0;  /* set the solution temperature to the 298 K if it is zero - safety feature */  
    }


  if (mp->PorosityModel == CONSTANT )                  /* constant porosity */ 
    {
      e = mp->porosity;                                         
    }
  else                                                /* non-CONSTANT porosity model */  
    {
      if (mp->PorosityModel == THERMAL_BATTERY)       /* thermal battery porosity model */ 
        {
          e = mp->u_porosity[0];
        }
      else
        {
	 /* need to implement non-constant porosity model */ 
         EH(-1, "other non-CONSTANT porosity model to be implemented");
        }
    }

  if (mp->DensityModel == CONSTANT )
    {
      rho   = mp->density;  /* get density from the xxx.mat file if it is constant */
      memset(d_rho->T, 0, sizeof(double)*MDE);
      memset(d_rho->F, 0, sizeof(double)*MDE);
      memset(d_rho->C, 0, sizeof(double)*MAX_CONC*MDE);
    }
  else 
    {
      /* get density from density() routine defined in mm_fill_terms.c */ 
      rho = density(d_rho, time);
    }

  for ( i=0; i<n_species; i++)
    {
      M[i] = mp->molecular_weight[i];
      for ( j=0; j<n_species; j++) 
	{
	  D[i][j] = mp->diffusivity_Stefan_Maxwell[i][j];
          EE = mp->u_diffusivity_Stefan_Maxwell[i][j][1];
          T0 = mp->u_diffusivity_Stefan_Maxwell[i][j][2];
          alpha = exp(-(EE/R)*(1.0/T-1.0/T0));  /* Arrhenius model */
	  D[i][j] *= alpha; /* Arrhenius S-M diffusivity model, KSC 9/04 */
	}
    }
  
  sumx = 0.0;
  for ( i=0; i<n_species-1; i++)
    {
      x[i] = fv->c[i];
#ifdef DEBUG_HKM
      if (x[i] <= 0.0) {
        fprintf(stderr,
		"Stefan_Maxwell_diff_flux WARNING P_%d: x[%d] = %g\n",
		ProcID, i, x[i]);
	EH(-1, "Stefan_Maxwell_Diff_flux zero or neg species");
      }
#endif
      sumx += x[i]; 
    }
  x[n_species-1] = 1.0 - sumx;

  M_mix = 0.0; 
  for ( i=0; i<n_species; i++)
    {
      M_mix += M[i]*x[i]; 
    }
  c = rho/M_mix; 

/* In order to obtain a volume flux for use in the species balance equation, we must divide the molar flux
   by c; the easiest way to do that is to set c = 1 here.  This will yield a flux that is analogous to the
   one computed by the Fick's law routine; see Equation (6.23) in the GOMA manual. */

  if (volume_flag)  /*  RSL 8/24/00  */
    {
     c = 1.0;
    }
  
  for ( i=0; i<n_species-1; i++)
    {
      for ( a=0; a<VIM; a++)
	{
	  grad_x[i][a] = fv->grad_c[i][a];         /* more-fraction or concentration gradient */ 
	}
    }
  
  for ( a=0; a<VIM; a++)
    {
      sumdelx = 0.0;
      for ( i=0; i<n_species-1; i++)
	{
	  sumdelx += grad_x[i][a]; 
	}
      grad_x[n_species-1][a] = - sumdelx;
    }

                         /* gradient of electrical potential in electrolyte */ 
  for ( a=0; a<VIM; a++)
    {
      if( pd_glob[mn]->MassFluxModel == STEFAN_MAXWELL_CHARGED || volume_flag)  /*  RSL 8/24/00  */
          {
             grad_phi2[a] = fv->grad_V[a];
          }
        else
          {
             grad_phi2[a] = 0.0;
          }
    }
  
  for ( i=0; i<n_species; i++)
    {
      z[i] = mp->charge_number[i];  
      for ( a=0; a<VIM; a++)
	{
           grad_mu[i][a] = z[i]*F*grad_phi2[a] + (R*T/x[i])*grad_x[i][a];    
                                      /* electrochemical potential gradient */
	}     
    }

  if ( af->Assemble_Residual )
    { 
      
      switch ( VIM )
	{
	case 1:  /* 1-D approximation */ 
	  B[0] = 0.0; 
	  for (j=0; j<n_species; j++)
	    {
	      A[0][j] = M[j];
	    }
	  for (i=1; i<n_species; i++)
	    {
	      B[i] = x[i]*grad_mu[i][0]/R/T; 
	      for ( j=0; j<n_species; j++)
		{
		  if(j != i) 
		    {
		      A[i][j] = x[i]/(c*e*D[i][j]); 
		    }
		  else if (j == i)
		    {
		      A[i][i] = 0.0;
		      for ( k=0; k<n_species; k++)
			{
			  if(k != i) A[i][i] -= x[k]/(c*e*D[i][k]); 
			}
		    }
		}
	    }   
	  break; 
	  
	case 2:  /* 2-D approximation */ 
	  B[0] = 0.0; 
	  B[1] = 0.0; 
	  for (k=0; k<n_species; k++)
	    {
	      j1=k*VIM;
	      j2=k*VIM+1;
	      A[0][j1] = M[k];
	      A[0][j2] = 0.0;
	      A[1][j1] = 0.0;
	      A[1][j2] = M[k];
	    }
	  
	  for (k=1; k<n_species; k++)
	    {
	      i1=k*VIM;
	      i2=k*VIM+1;
	      B[i1] = x[k]*grad_mu[k][0]/R/T;
	      B[i2] = x[k]*grad_mu[k][1]/R/T;
	      for (l=0; l<n_species; l++)
		{
		  j1=l*VIM;
		  j2=l*VIM+1;
		  if(l != k) 
		    {
		      A[i1][j1] = x[k]/(c*e*D[k][l]); 
		      A[i1][j2] = 0.0;
		      A[i2][j1] = 0.0;
		      A[i2][j2] = A[i1][j1]; 
		    }
		  else if (l == k)
		    {
		      A[i1][i1] = 0.0;
		      for ( m=0; m<n_species; m++)
			{
			  if(m != k) 
			    {
			      A[i1][i1] -= x[m]/(c*e*D[k][m]); 
			    }
			}
                A[i2][i2] = A[i1][i1];
                A[i1][i2] = 0.0;
                A[i2][i1] = 0.0;
		    }
		}
	    }   
	  break; 
	  
	case 3:  /* 3-D: full dimensions */ 
	  B[0] = 0.0; 
	  B[1] = 0.0; 
	  B[2] = 0.0; 
	  for (k=0; k<n_species; k++)
	    {
	      j1=k*VIM;
	      j2=k*VIM+1;
	      j3=k*VIM+2;
	      A[0][j1] = M[k];
	      A[0][j2] = 0.0;
	      A[0][j3] = 0.0;
	      A[1][j1] = 0.0;
	      A[1][j2] = M[k];
	      A[1][j3] = 0.0;
	      A[2][j1] = 0.0;
	      A[2][j2] = 0.0;
	      A[2][j3] = M[k];
	    }
	  
	  for (k=1; k<n_species; k++)
	    {
	      i1=k*VIM;
	      i2=k*VIM+1;
	      i3=k*VIM+2;
	      B[i1] = x[k]*grad_mu[k][0]/R/T;
	      B[i2] = x[k]*grad_mu[k][1]/R/T;
	      B[i3] = x[k]*grad_mu[k][2]/R/T;
	      for (l=0; l<n_species; l++)
		{
		  j1=l*VIM;
		  j2=l*VIM+1;
		  j3=l*VIM+2;
		  if(l != k) 
		    {
		      A[i1][j1] = x[k]/(c*e*D[k][l]); 
		      A[i1][j2] = 0.0;
		      A[i1][j3] = 0.0;
		      A[i2][j1] = 0.0;
		      A[i2][j2] = A[i1][j1]; 
		      A[i2][j3] = 0.0;
		      A[i3][j1] = 0.0;
		      A[i3][j2] = 0.0;
		      A[i3][j3] = A[i1][j1]; 
		    }
		  else if (l == k)
		    {
		      A[i1][i1] = 0.0;
		      for ( m=0; m<n_species; m++)
			{
			  if(m != k) 
			    {
			      A[i1][i1] -= x[m]/(c*e*D[k][m]); 
			    }
			}
		      A[i2][i2] = A[i1][i1];
		      A[i3][i3] = A[i1][i1];
		      A[i1][i2] = 0.0;
		      A[i1][i3] = 0.0;
		      A[i2][i1] = 0.0;
		      A[i2][i3] = 0.0;
		      A[i3][i1] = 0.0;
		      A[i3][i2] = 0.0;
		    }
		}
	    }   
	  break; 
	  
	}

                 /* LU decomposition of matrix A[1..n][1..n], then  */ 
                 /* forward and back substitution to get solution vector B */
      
      lu_decomp_backsub_driver ( A, B, indx, n, 1);
      
      for(i=0; i<n; i++)
	{
	  J[i] = B[i];  /* assign the solution vector bb to the flux vector */  
                   /* put the Stefan-Maxwell flux in the global fv structure */
          fv->SM_flux[i] = B[i];  
	} 
      
      for (j=0; j<n; j++)     /* calculate the inverse of aa; KSC: 10/98 */
        {
          for (i=0; i<n; i++) bbb[i] = 0.0;
          bbb[j] = 1.0;

          lu_decomp_backsub_driver ( A, bbb, indx, n, 0);

          for (i=0; i<n; i++) A_inv[i][j] = bbb[i];
        }

       /* printf("In mm_fill_species.c: the coeff. matrix A:\n");
          for (i=0; i<n; i++)
            {
               printf("%g %g %g %g %g %g\n", A[i][0], A[i][1], A[i][2],
                                            A[i][3], A[i][4], A[i][5]);
            }  */

       /* printf("In mm_fill_species.c: the inverse of matrix A:\n");   */

       for (i=0; i<n; i++)
         {
           for (j=0; j<n; j++)
             {
               fv->SM_matrix_inv[i][j] = A_inv[i][j];    
             }
    /* printf("%g %g %g %g %g %g\n", A_inv[i][0], A_inv[i][1], A_inv[i][2],
                                     A_inv[i][3], A_inv[i][4], A_inv[i][5]); */
         }

          /* check on the correctness of A_inv; KSC: 10/22/98 */
          /* for (i=0; i<n; i++)
            {
              for (j=0; j<n; j++)
                {
                  cc[i][j] = 0.0;
                  for (k=0; k<n; k++)
                    {
                      cc[i][j] += A[i][k]*A_inv[k][j];
                    }
                }
             }

          printf("In mm_fill_species.c: the product of coefficient matrix A and its inverse:\n");
          for (i=0; i<n; i++)
            {
              printf("%g %g %g %g %g %g\n", cc[i][0], cc[i][1], cc[i][2],
                                            cc[i][3], cc[i][4], cc[i][5]);
            }  */

      /* to check on the solution accuracy: sum of (MiJi) = 0 ? */ 
      /* printf("\n");      
	 printf("For the x-component fluxes:\n");    
	 sumMJ=0.0;
	 for(k=0; k<n_species; k++) 
	 {
	 sumMJ += M[k]*J[VIM*k];
	 printf("i=, M[i]=, J[i]= %d %e %e\n", k, M[k], J[VIM*k]);      
	 }
	 printf("sumMJ= %e\n", sumMJ);      
	 
	 if(VIM >1)
	 {
	 printf("\n");       
	 printf("For the y-component fluxes:\n");  
	 sumMJ=0.0;
	 for(k=0; k<n_species; k++) 
	 {
	 sumMJ += M[k]*J[VIM*k+1];
	 printf("i=, M[i]=, J[i]= %d %e %e\n", k, M[k], J[VIM*k+1]);    
	 }
	 printf("sumMJ= %e\n", sumMJ);     
	 }  
	 
	 if(VIM >2)
	 {
	 printf("\n");   
	 printf("For the z-component fluxes:\n");   
	 sumMJ=0.0;
	 for(k=0; k<n_species; k++) 
	 {
	 sumMJ += M[k]*J[VIM*k+2];
	 printf("i=, M[i]=, J[i]= %d %e %e\n", k, M[k], J[VIM*k+2]);     
	 }
	 printf("sumMJ= %e\n", sumMJ);     
	 } */   
      
      for ( w=0; w<n_species-1; w++)
	{
	  switch ( VIM )
	    {
	    case 1:
	      st->diff_flux[w][0] = B[w];
	      break;
	      
	    case 2:
	      st->diff_flux[w][0] = B[VIM*w  ];
	      st->diff_flux[w][1] = B[VIM*w+1];
	      break;
	      
	    case 3:
	      st->diff_flux[w][0] = B[VIM*w  ];
	      st->diff_flux[w][1] = B[VIM*w+1];
	      st->diff_flux[w][2] = B[VIM*w+2];
	      break;
	    }
	}
    } /* end if af->Assemble_Residual */
  
  if ( af->Assemble_Jacobian )
    { 
      /* first calculate dA/dx and dB/dx where A is the coefficent matrix and B is the
         driving-force vector in the AJ = B matrix equation system; x is mole fraction */
      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	{
	  for ( w1=0; w1<pd->Num_Species_Eqn; w1++) 
	    {
	      correction = d_rho->C[w1][j]/rho -
	                   bf[var]->phi[j]*(M[w1] - M[n_species-1])/M_mix; /* RSL 8/24/00 */
	      switch ( VIM )
		{
		case 1:  /* 1-D approximation */ 
		  for (wi=0; wi<n; wi++) /* first zero the dA matrix and dB vector */
		    {
		      dB[wi]=0.0;
		      for (wj=0; wj<n; wj++)
			{
			  dA[wi][wj]=0.0;
			}
		    }
		  for (wi=1; wi<n_species-1; wi++)
		    {
		      if(wi == w1) 
			{
			  dB[wi] = bf[var]->grad_phi[j][0];
			  dB[wi] += z[wi]*F*grad_phi2[0]*bf[var]->phi[j]/(R*T); /* RSL 3/23/00 */
			}
		      for ( wj=0; wj<n_species; wj++)
			{
			  if(wj != wi) 
			    {
			      if(w1 == wi) 
				{
				  dA[wi][wj] = bf[var]->phi[j]/(c*e*D[wi][wj]); 
				  if (!volume_flag)  /* RSL 8/24/00 */
				     {
				      dA[wi][wj] -= correction * A[wi][wj];
				     }
				}
			    }
			  else if (wj == wi && w1 != wi)
			    {
			      dA[wi][wi] = -bf[var]->phi[j]/(c*e*D[wi][w1]);
			      dA[wi][wi] += bf[var]->phi[j]/(c*e*D[wi][n_species-1]); /* RSL 8/24/00 */
			      if (!volume_flag)  /* RSL 8/24/00 */
			         {
			          dA[wi][wi] -= correction * A[wi][wi];
			         }
			    }
			  else /* RSL 8/24/00 */
			    {
			      dA[wi][wi] = bf[var]->phi[j]/(c*e*D[wi][n_species-1]);
			      if (!volume_flag)
			         {
			          dA[wi][wi] -= correction * A[wi][wi];
			         }
			    }
			}
		    }
		  dB[n-1] = -bf[var]->grad_phi[j][0];
		  dB[n-1] -= z[n-1]*F*grad_phi2[0]*bf[var]->phi[j]/(R*T); /* RSL 3/23/00 */
		  for (wj=0; wj<n_species-1; wj++)
		    {
		      dA[n-1][wj] = -bf[var]->phi[j]/(c*e*D[n_species-1][wj]);
		      if (!volume_flag)  /* RSL 8/24/00 */
		         {
		          dA[n-1][wj] -= correction * A[n-1][wj];
		         }
		    }
/*		  dA[n-1][n-1] = -bf[var]->phi[j]/(c*e*D[n_species-1][n_species-2]);       */
		  dA[n-1][n-1] = -bf[var]->phi[j]/(c*e*D[n_species-1][w1]); /* RSL 8/24/00 */
		  if (!volume_flag)  /* RSL 8/24/00 */
		     {
		      dA[n-1][n-1] -= correction * A[n-1][n-1];
		     }
		  break; 
		  
		case 2:  /* 2-D approximation */ 
		  for (wi=0; wi<n; wi++) /* first zero the dA matrix and dB vector */ 
		    {
		      dB[wi]=0.0;
		      for (wj=0; wj<n; wj++)
			{
			  dA[wi][wj]=0.0;
			}
		    }
		  
		  for (k=1; k<n_species-1; k++)  /* next fill in the non-zero entries */ 
		    {
		      i1=k*VIM;
		      i2=k*VIM+1;
		      if(k == w1) 
			{
			  dB[i1] = bf[var]->grad_phi[j][0];
			  temp_var = z[k]*F*bf[var]->phi[j]/(R*T); /* RSL 8/23/00 */
			  dB[i1] += temp_var*grad_phi2[0];         /* RSL 3/23/00 */
			  dB[i2] = bf[var]->grad_phi[j][1];
			  dB[i2] += temp_var*grad_phi2[1];         /* RSL 3/23/00 */
			}
		      for (l=0; l<n_species; l++)
			{
			  j1=l*VIM;
			  j2=l*VIM+1;
			  if(l != k) 
			    {
			      if (k == w1) 
				{ 
				  dA[i1][j1] = bf[var]->phi[j]/(c*e*D[k][l]);
				  if (!volume_flag)  /* RSL 8/24/00 */
				     {
				      dA[i1][j1] -= correction * A[i1][j1];
				     }
				  dA[i2][j2] = dA[i1][j1]; 
				}
			    }
			  else if (l == k  &&  k != w1)
			    {
			      dA[i1][i1] = -bf[var]->phi[j]/(c*e*D[k][w1]); 
			      dA[i1][i1] += bf[var]->phi[j]/(c*e*D[k][n_species-1]); /* RSL 8/24/00 */
			      if (!volume_flag)  /* RSL 8/24/00 */
			         {
			          dA[i1][i1] -= correction * A[i1][i1];
			         }
			      dA[i2][i2] = dA[i1][i1]; 
			    }
			  else /* RSL 8/24/00 */
			    {
			      dA[i1][i1] = bf[var]->phi[j]/(c*e*D[k][n_species-1]);
			      if (!volume_flag)
			         {
			          dA[i1][i1] -= correction * A[i1][i1];
			         }
			      dA[i2][i2] = dA[i1][i1];
			    }
			}
		    }   
		  dB[n-2] = -bf[var]->grad_phi[j][0];
		  temp_var = z[n_species-1]*F*bf[var]->phi[j]/(R*T); /* RSL 8/23/00 */
		  dB[n-2] -= temp_var*grad_phi2[0];                  /* RSL 3/23/00 */
		  dB[n-1] = -bf[var]->grad_phi[j][1];  
		  dB[n-1] -= temp_var*grad_phi2[1];                  /* RSL 3/23/00 */
		  for (l=0; l<n_species-1; l++)
		    {
		      j1=l*VIM;
		      j2=l*VIM+1;
		      dA[n-2][j1] = -bf[var]->phi[j]/(c*e*D[n_species-1][l]);
		      if (!volume_flag)  /* RSL 8/24/00 */
		         {
		          dA[n-2][j1] -= correction * A[n-2][j1];
		         }
		      dA[n-1][j2] = dA[n-2][j1];
		    }
/*		  dA[n-2][n-2] = -bf[var]->phi[j]/(c*e*D[n_species-1][n_species-2]);       */
		  dA[n-2][n-2] = -bf[var]->phi[j]/(c*e*D[n_species-1][w1]); /* RSL 8/24/00 */
		  if (!volume_flag)  /* RSL 8/24/00 */
		     {
		      dA[n-2][n-2] -= correction * A[n-2][n-2];
		     }
		  dA[n-1][n-1] = dA[n-2][n-2];
		  break; 
		  
		case 3:  /* 3-D: full dimensions */ 
		  for (wi=0; wi<n; wi++) /* first zero the dA matrix and dB vector */
		    {
		      dB[wi]=0.0;
		      for (wj=0; wj<n; wj++)
			{
			  dA[wi][wj]=0.0;
			}
		    }
		  for (k=1; k<n_species-1; k++)
		    {
		      i1=k*VIM;
		      i2=k*VIM+1;
		      i3=k*VIM+2;
		      if(k == w1)
			{
			  dB[i1] = bf[var]->grad_phi[j][0];
			  temp_var = z[k]*F*bf[var]->phi[j]/(R*T); /* RSL 8/23/00 */
			  dB[i1] += temp_var*grad_phi2[0];         /* RSL 3/23/00 */
			  dB[i2] = bf[var]->grad_phi[j][1];
			  dB[i2] += temp_var*grad_phi2[1];         /* RSL 3/23/00 */
			  dB[i3] = bf[var]->grad_phi[j][2];
			  dB[i3] += temp_var*grad_phi2[2];         /* RSL 3/23/00 */
			}
		      for (l=0; l<n_species; l++)
			{
			  j1=l*VIM;
			  j2=l*VIM+1;
			  j3=l*VIM+2;
			  if(l != k) 
			    {
			      if (k == w1)
				{
				  dA[i1][j1] = bf[var]->phi[j]/(c*e*D[k][l]);
				  if (!volume_flag)  /* RSL 8/24/00 */
				     {
				      dA[i1][j1] -= correction * A[i1][j1];
				     }
				  dA[i2][j2] = dA[i1][j1];
				  dA[i3][j3] = dA[i1][j1];
				}
			    }
			  else if (l == k  &&  k != w1)
			    {
			      dA[i1][i1] = -bf[var]->phi[j]/(c*e*D[k][w1]);
			      dA[i1][i1] += bf[var]->phi[j]/(c*e*D[k][n_species-1]); /* RSL 8/24/00 */
			      if (!volume_flag)  /* RSL 8/24/00 */
			         {
			          dA[i1][i1] -= correction * A[i1][i1];
			         }
			      dA[i2][i2] = dA[i1][i1];
			      dA[i3][i3] = dA[i1][i1];
			    }
			  else /* RSL 8/24/00 */
			    {
			      dA[i1][i1] = bf[var]->phi[j]/(c*e*D[k][n_species-1]);
			      if (!volume_flag)
			         {
			          dA[i1][i1] -= correction * A[i1][i1];
			         }
			      dA[i2][i2] = dA[i1][i1];
			      dA[i3][i3] = dA[i1][i1];
			    }
			}
		    }   
		  dB[n-3] = -bf[var]->grad_phi[j][0];
		  temp_var = z[n_species-1]*F*bf[var]->phi[j]/(R*T); /* RSL 8/23/00 */
		  dB[n-3] -= temp_var*grad_phi2[0];                  /* RSL 3/23/00 */
		  dB[n-2] = -bf[var]->grad_phi[j][1];
		  dB[n-2] -= temp_var*grad_phi2[1];                  /* RSL 3/23/00 */
		  dB[n-1] = -bf[var]->grad_phi[j][2];
		  dB[n-1] -= temp_var*grad_phi2[2];                  /* RSL 3/23/00 */
		  for (l=0; l<n_species-1; l++)
		    {
		      j1=l*VIM;
		      j2=l*VIM+1;
		      j3=l*VIM+2;
		      dA[n-3][j1] = -bf[var]->phi[j]/(c*e*D[n_species-1][l]);
		      if (!volume_flag)  /* RSL 8/24/00 */
		         {
		          dA[n-3][j1] -= correction * A[n-3][j1];
		         }
		      dA[n-2][j2] = dA[n-3][j1];
		      dA[n-1][j3] = dA[n-3][j1];
		    }
/*		  dA[n-3][n-3] = -bf[var]->phi[j]/(c*e*D[n_species-1][n_species-2]);       */
		  dA[n-3][n-3] = -bf[var]->phi[j]/(c*e*D[n_species-1][w1]); /* RSL 8/24/00 */
		  if (!volume_flag)  /* RSL 8/24/00 */
		     {
		      dA[n-3][n-3] -= correction * A[n-3][n-3];
		     }
		  dA[n-2][n-2] = dA[n-3][n-3];
		  dA[n-1][n-1] = dA[n-3][n-3];
		  break; 
		  
		} /* end switch; computation of dA and dB is complete for the current values of
		     w1 and j.  Now we need to compute the current contribution to the giant
		     C matrix and store it for use in mm_fill_potential.c, just as we stored
		     A_inv -- RSL 3/31/00 */

	      for (wi=0; wi<n; wi++)
	         {
	          for (wj=0; wj<n; wj++)
	             {
	              scratch[wi][wj]=0.0;
	              for (k=0; k<n; k++)
	                 {
	                  scratch[wi][wj] -= fv->SM_matrix_inv[wi][k] * dA[k][wj];
	                 }
	             }
	         }

	      for (wi=0; wi<n; wi++)
	         {
	          for (wj=0; wj<n; wj++)
	             {
	              fv->giant_C_matrix[w1][j][wi][wj]=0.0;
	              for (k=0; k<n; k++)
	                 {
	                  fv->giant_C_matrix[w1][j][wi][wj] += scratch[wi][k] *
	                                                       fv->SM_matrix_inv[k][wj];
	                 }
	             }
	         }

/* end of computation of giant C matrix -- RSL 3/31/00 */
	      
	      for(wi=0; wi<n; wi++)
		{
		  dAJ[wi] = 0.0;
		  for (k=0; k<n; k++) 
		    {
		      dAJ[wi] += dA[wi][k]*J[k];
		    }
		  B[wi] = dB[wi] - dAJ[wi];
		} 
	      
                /* forward and back substitution to get solution vector B */ 
              lu_decomp_backsub_driver ( A, B, indx, n, 0);
	      
	      switch ( VIM )
		{
		case 1:
		  for ( w=0; w<n_species-1; w++)
		    {
		      st->d_diff_flux_dc[w][0][w1][j] = B[w];
		    }
		  break;
		  
		case 2:
		  for ( w=0; w<n_species-1; w++)
		    {
		      st->d_diff_flux_dc[w][0][w1][j] = B[VIM*w  ];
		      st->d_diff_flux_dc[w][1][w1][j] = B[VIM*w+1];
		    }
		  break;
		  
		case 3:
		  for ( w=0; w<n_species-1; w++)
		    {
		      st->d_diff_flux_dc[w][0][w1][j] = B[VIM*w  ];
		      st->d_diff_flux_dc[w][1][w1][j] = B[VIM*w+1];
		      st->d_diff_flux_dc[w][2][w1][j] = B[VIM*w+2];
		    }
		  break;
		  
		}
	      
	    }
	}
      
      /* dJ/dmesh */
      for (q=0; q < pd->Num_Dim; q++)
	{
	  var = MESH_DISPLACEMENT1 + q;
	  if ( pd->v[pg->imtrx][var] )
	    {
	      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		{
		  switch ( VIM )
		    {
		    case 1:  /* 1-D approximation */ 
		      dB[0] = 0.0; 
		      for (wi=1; wi<n_species-1; wi++)
			{
			  dB[wi] = fv->d_grad_c_dmesh[0][wi][q][j];
			}   
		      dB[n_species-1] = 0.0 ;
		      for (wi=0; wi<n_species-1; wi++)
			{
			  dB[n_species-1] -= fv->d_grad_c_dmesh[0][wi][q][j];
			}   
		      break; 
		      
		    case 2:  /* 2-D approximation */ 
		      for (wi=0; wi<n; wi++) /* first zero the dB vector */
			{
			  dB[wi]=0.0;
			}
		      for (k=1; k<n_species-1; k++) /* next fill in the non-zero entries */ 
			{
			  i1=(k+1)*VIM-2;
			  i2=(k+1)*VIM-1;
			  dB[i1] = fv->d_grad_c_dmesh[0][k][q][j];
			  dB[i2] = fv->d_grad_c_dmesh[1][k][q][j];
			}   
		      dB[n_species*VIM-2] = 0.0;
		      dB[n_species*VIM-1] = 0.0;
		      for (k=0; k<n_species-1; k++)
			{
			  dB[n_species*VIM-2] -= fv->d_grad_c_dmesh[0][k][q][j];
			  dB[n_species*VIM-1] -= fv->d_grad_c_dmesh[1][k][q][j];
			}   
		      break; 
		      
		    case 3:  /* 3-D: full dimensions */ 
		      dB[0] = 0.0; 
		      dB[1] = 0.0; 
		      dB[2] = 0.0; 
		      for (k=1; k<n_species-1; k++)
			{
			  i1=(k+1)*VIM-3;
			  i2=(k+1)*VIM-2;
			  i3=(k+1)*VIM-1;
			  dB[i1] = fv->d_grad_c_dmesh[0][k][q][j];
			  dB[i2] = fv->d_grad_c_dmesh[1][k][q][j];
			  dB[i3] = fv->d_grad_c_dmesh[2][k][q][j];
			}   
		      dB[n_species*VIM-3] = 0.0 ;
		      dB[n_species*VIM-2] = 0.0 ;
		      dB[n_species*VIM-1] = 0.0 ;
		      for (k=0; k<n_species-1; k++)
			{
			  i1=(k+1)*VIM-3;
			  i2=(k+1)*VIM-2;
			  i3=(k+1)*VIM-1;
			  dB[n_species*VIM-3] -= fv->d_grad_c_dmesh[0][k][q][j];
			  dB[n_species*VIM-2] -= fv->d_grad_c_dmesh[1][k][q][j];
			  dB[n_species*VIM-1] -= fv->d_grad_c_dmesh[2][k][q][j];
			}   
		      break; 
		      
		    }
		  for(wi=0; wi<n; wi++)
		    {
		      B[wi] = dB[wi];
		    } 

		  /* forward and back substitution to get solution vector B */ 
                  lu_decomp_backsub_driver ( A, B, indx, n, 0);
		  
		  switch ( VIM )
		    {
		    case 1:
		      for ( w=0; w<n_species-1; w++)
			{
			  st->d_diff_flux_dmesh[w][0][q][j] = B[w];
			}
		      break;
		      
		    case 2:
		      for ( w=0; w<n_species-1; w++)
			{
			  st->d_diff_flux_dmesh[w][0][q][j] = B[VIM*w  ];
			  st->d_diff_flux_dmesh[w][1][q][j] = B[VIM*w+1];
			}
		      break;
		      
		    case 3:
		      for ( w=0; w<n_species-1; w++)
			{
			  st->d_diff_flux_dmesh[w][0][q][j] = B[VIM*w  ];
			  st->d_diff_flux_dmesh[w][1][q][j] = B[VIM*w+1];
			  st->d_diff_flux_dmesh[w][2][q][j] = B[VIM*w+2];
			}
		      break;
		    }
		  
		}
	    }
	}
      
      var = TEMPERATURE;
      for ( a=0; a<VIM; a++)
	{
	  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	    {
	      for ( w=0; w<pd->Num_Species_Eqn; w++)
		{
		  st->d_diff_flux_dT[w][a] [j] =  0.0;
		}
	    }
	}
      
      var = PRESSURE;
      
      for ( a=0; a<VIM; a++)
	{
	  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	    {
	      for ( w=0; w<pd->Num_Species_Eqn; w++)
		{
		  st->d_diff_flux_dP[w][a] [j] = 0.;
		}
	    }
	}

/* We need derivatives of the diffusive flux with respect to voltage.  Evidently the
   Jacobian terms associated with these derivatives have not been considered in GOMA
   before.  RSL 4/4/00 */
      
      var = VOLTAGE;

      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
         {
          for ( a=0; a<VIM; a++)
             {
              dB[a] = 0.0;
              for ( l=1; l<n_species; l++)
                 {
                  dB[l*VIM+a] = F*x[l]*z[l]*bf[var]->grad_phi[j][a]/(R*T);
                 }
             }
          for (wi=0; wi<n; wi++)
             {
              dJ[wi] = 0.0;
              for (k=0; k<n; k++)
                 {
                  dJ[wi] += fv->SM_matrix_inv[wi][k] * dB[k];
                 }
             }
          for ( a=0; a<VIM; a++)
             {
              for ( w=0; w<pd->Num_Species_Eqn; w++)
                 {
                  st->d_diff_flux_dV[w][a][j] = dJ[w*VIM+a];
                 }
             }
         }

    } /* end if ( af->Assemble_Jacobian) */

  return (status);  
} /* END of Stefan_Maxwell_diff_flux; KSC: 8/98 */
/********************************************************************************/
/********************************************************************************/
/********************************************************************************/

int
fickian_charged_flux (struct Species_Conservation_Terms *st, int w)

/*******************************************************************************
*
*  A routine for computing the total species fluxes (diffusive and migration) and
*  the associated Jacobian sensitivity entries. The total fluxes are loaded 
*  into st->diff_flux[][] whereas the Jacobian entries are loaded into
*  st->d_diff_flux[][][]. 
*
*
*   If species var type is, the flux has units of:
*        SPECIES_MOLE_FRACTION    -> mol cm-2 sec-1
*            coeff_rho = C_avg
*        SPECIES_MASS_FRACTION    -> gm  cm-2 sec-1
*            coeff_rho = rho
*        SPECIES_DENSITY          -> gm  cm-2 sec-1
*            coeff_rho = 1.0
*        SPECIES_UNDEFINED_FORM   -> nominally cm sec-1 (if diffusivity
*                                        is in cm^2 sec-1)
*            coeff_rho = 1.0
*           (However, units are left up to the user, in practice)
*        SPECIES_CONCENTRATION    -> mol cm-2 sec-1
*            coeff_rho = 1.0
*
*     Author: Ken S. Chen, 9/2000 (clone after fickian_flux and 
*                                  Stefan_Maxwell_diff_flux)
*
******************************************************************************/
{
  const double F = 96487.0;       /* Faraday's constant  ( C/equiv ) */
  const double R = 8.314;   /* Universal gas constant in units of J/mole K */
  dbl T=298.0;              /* default Electrolyte solution temperature */
  dbl FRT;                  /* product of F/R/T */ 
  int w1, a, j, q, var;
  double FRTzD, avg_molec_weight, coeff_rho;

  /*
   *  Add in rho or C depending upon species variable type
   */
  coeff_rho = 1.0;
  if (mp->Species_Var_Type == SPECIES_MASS_FRACTION) {
    coeff_rho = mp->density;
  } else if (mp->Species_Var_Type == SPECIES_MOLE_FRACTION) {
    for (w1 = 0, avg_molec_weight = 0.0; w1 < pd->Num_Species; w1++) {
      avg_molec_weight += mp->molecular_weight[w1] * fv->c[w];
    }
    coeff_rho = mp->density / avg_molec_weight;
  }


  if (mp->SolutionTemperatureModel == CONSTANT)  {
    T = mp->solution_temperature;
  } else {
    EH(-1, "Solution-temperature model other than CONSTANT awaits future implementation");
  }
  /* set solution temperature to 298 K if it is zero - safety feature */
  if (T == 0.0) T = 298.0;
  FRT = F/(R * T);
  FRTzD = FRT *  mp->charge_number[w] * mp->diffusivity[w];

  for (a = 0; a < VIM; a++) {
    st->diff_flux[w][a] = - coeff_rho *
	(mp->diffusivity[w] * fv->grad_c[w][a] +
	 FRTzD * fv->c[w] * fv->grad_V[a]);
  }

  if (af->Assemble_Jacobian) { 
    var = MASS_FRACTION;
    for (a = 0; a < VIM; a++) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
	for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++) {
	  st->d_diff_flux_dc[w][a][w1][j] = 0.0;
	}
	st->d_diff_flux_dc[w][a][w][j] = - coeff_rho * 
	    (mp->diffusivity[w] * bf[var]->grad_phi[j][a] +
	     FRTzD * bf[var]->phi[j] * fv->grad_V[a]);
      }
    }
      
    for (q = 0; q < pd->Num_Dim; q++) {
      var = MESH_DISPLACEMENT1 + q;
      if (pd->v[pg->imtrx][var]) {
	for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
	  for (a = 0; a < VIM; a++) {
	    st->d_diff_flux_dmesh[w][a] [q][j] = - coeff_rho *
		(mp->diffusivity[w] * fv->d_grad_c_dmesh[a][w] [q][j] +
		FRTzD * fv->c[w] * fv->d_grad_V_dmesh[a][q][j]);
	  }
	}
      }
    }
      
    var = TEMPERATURE;
    for (a = 0; a < VIM; a++) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
	st->d_diff_flux_dT[w][a][j] = coeff_rho * (FRTzD/T) * 
	    fv->c[w] * bf[var]->phi[j] * fv->grad_V[a];
      }
    }

    var = VOLTAGE;
    for (a = 0; a < VIM; a++) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
	st->d_diff_flux_dV[w][a][j] = - coeff_rho * FRTzD * 
	    fv->c[w] * bf[var]->grad_phi[j][a]; 
      }
    }
  }
  return (0);
} /* END of fickian_charged_flux */
/*******************************************************************************/

/******************************************************************************
*
*  Function that computes the fluxes and corresponding sensitivities for a
*  system of charged species obeying Fick's law -- mole fraction version
*
*  RSL 9/18/00
*
******************************************************************************/
int
fickian_charged_flux_x (struct Species_Conservation_Terms *st, double time,
                        double dt)
{
  int a, j, q, var;
  int w;
  int status = 0;

  const double R = 8.314;         /* Universal gas constant ( J/mole K ) */
  const double F = 96487.0;       /* Faraday's constant  ( C/equiv ) */
  dbl T = 298.0;                  /* temperature; set default value to room temperature */
  dbl rho;                        /* density */
  dbl c;                          /* total molar concentration */
  dbl x[MAX_CONC];                /* mole fraction */
  dbl M[MAX_CONC], M_mix;         /* molecular weight of individual species and of mixture
 */
  dbl z[MAX_CONC];                /* charge number */
  dbl D[MAX_CONC];                /* Fickian diffusivities */
  dbl grad_phi2[DIM];             /* gradient of electrical potential in electrolyte */
  dbl grad_x[MAX_CONC][DIM];      /* mole fraction gradient */

  int n_species;

  int mn, i;
  dbl sumx, frt = 0.0, derivative, save1, save2;
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;  /* density dependence */
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;

  n_species = pd->Num_Species_Eqn + 1;
  mn = ei[pg->imtrx]->mn;

  if (pd_glob[mn]->e[R_ENERGY])  /* if the energy equation is being solved */
     {
      T = fv->T;
     }
  else                           /* if the energy equation is NOT being solved */
     {
      if (mp->SolutionTemperatureModel == CONSTANT)
         {
          T = mp->solution_temperature;
         }
      else if (mp->SolutionTemperatureModel == THERMAL_BATTERY)
         {
          electrolyte_temperature(time, dt, 0);
          T = mp->electrolyte_temperature;
         }
      else
         {
          EH(-1, "Invalid solution temperature model");
         }
     }

  rho = density(d_rho, time); /*  RSL 6/22/02  */

  for ( i=0; i<n_species; i++)
     {
      M[i] = mp->molecular_weight[i];
      z[i] = mp->charge_number[i];
     }

  sumx = 0.0;
  for ( i=0; i<n_species-1; i++)
     {
      D[i] = mp->diffusivity[i];
      x[i] = fv->c[i];
      sumx += x[i];
     }
  x[n_species-1] = 1.0 - sumx;

  M_mix = 0.0;
  for ( i=0; i<n_species; i++)
     {
      M_mix += M[i]*x[i];
     }
  c = rho/M_mix;

  for ( i=0; i<n_species-1; i++)
     {
      for ( a=0; a<VIM; a++)
            {
             grad_x[i][a] = fv->grad_c[i][a];
            }
     }

  for ( a=0; a<VIM; a++)
     {
      grad_phi2[a] = fv->grad_V[a]; /* gradient of electrical potential in electrolyte */
     }

  if ( af->Assemble_Residual )
     {
      frt = F/(R*T);
      for ( i=0; i<n_species-1; i++)
          {
           for ( a=0; a<VIM; a++)
              {
               st->diff_flux[i][a] = -D[i]*c*(grad_x[i][a] + z[i]*frt*x[i]*grad_phi2[a]);
              }
          }
     }

  if ( af->Assemble_Jacobian )
     {
      var = MASS_FRACTION;
      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
         {
          for ( w=0; w<pd->Num_Species_Eqn; w++)
             {
              derivative = d_rho->C[w][j]/rho - bf[var]->phi[j]*(M[w] - M[n_species-1])/M_mix;
              for ( i=0; i<pd->Num_Species_Eqn; i++)
                 {
                  for ( a=0; a<VIM; a++)
                     {
                      st->d_diff_flux_dc[i][a][w][j] = st->diff_flux[i][a]*derivative;
                      if ( i == w )
                         {
                          save1 = bf[var]->grad_phi[j][a];
                          save2 = z[i]*frt*bf[var]->phi[j]*grad_phi2[a];
                          st->d_diff_flux_dc[i][a][w][j] -= D[i]*c*(save1 + save2);
                         }
                     }
                 }
             }
         }

      for ( q=0; q<pd->Num_Dim; q++ )
         {
          var = MESH_DISPLACEMENT1 + q;
          if ( pd->v[pg->imtrx][var] )
             {
              for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
                 {
                  for ( i=0; i<pd->Num_Species_Eqn; i++)
                     {
                      for ( a=0; a<VIM; a++)
                         {
                          st->d_diff_flux_dmesh[i][a][q][j] = -D[i]*c*(
                              fv->d_grad_c_dmesh[a][i] [q][j] +
                              z[i]*frt*x[i]*fv->d_grad_V_dmesh[a][q][j]);
                         }
                     }
                 }
             }
         }

      var = TEMPERATURE;
      for ( i=0; i<pd->Num_Species_Eqn; i++)
         {
          save2 = D[i]*c*z[i]*frt*x[i]/T;
          for ( a=0; a<VIM; a++)
             {
              save1 = save2*grad_phi2[a];
              for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
                 {
                  st->d_diff_flux_dT[i][a][j] =  save1*bf[var]->phi[j];
                 }
             }
         }

      var = PRESSURE;
      for ( i=0; i<pd->Num_Species_Eqn; i++)
         {
          for ( a=0; a<VIM; a++)
             {
              for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
                 {
                  st->d_diff_flux_dP[i][a][j] = 0.;
                 }
             }
         }

      var = VOLTAGE;
      for ( i=0; i<pd->Num_Species_Eqn; i++)
         {
          save2 = -D[i]*c*z[i]*frt*x[i];
          for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
             {
              for ( a=0; a<VIM; a++)
                 {
                  st->d_diff_flux_dV[i][a][j] = save2*bf[var]->grad_phi[j][a];
                 }
             }
         }

     } /* end if (af->Assemble_Jacobian) */

  return (status);

} /* END of fickian_charged_flux_x */

/*******************************************************************************/
/*******************************************************************************/

int
fickian_flux (struct Species_Conservation_Terms *st, int w)

    /***************************************************************************
     *
     *  Function that computes the diffusive flux of species w and its
     *  Jacobian terms.
     *  It loads the jacobian terms into st->d_diff_flux[][][]
     *
     *   The units of this term, as well as all of the other terms in the
     *   species conservation equation, depend upon the value of
     *   species_var_type.
     *
     *   The exact form of the flux is as follows:
     *   If species var type is:
     *        SPECIES_MOLE_FRACTION    -> mol cm-2 sec-1
     *            diff_flux_i = - C_avg * D_i del X_i
     *        SPECIES_MASS_FRACTION    -> gm  cm-2 sec-1
     *            diff_flux_i = - rho * D_i del Y_i
     *        SPECIES_DENSITY          -> gm  cm-2 sec-1
     *            diff_flux_i = -  D_i del Rho_i
     *        SPECIES_UNDEFINED_FORM   -> nominally cm sec-1 (if diffusivity
     *                                        is in cm^2 sec-1)
     *            diff_flux_i = - D_i del Y_i
     *           (However, units are left up to the user, in practice)
     *        SPECIES_CONCENTRATION    -> mol cm-2 sec-1
     *            diff_flux_i = - D_i del C_i
     *
     * Note, except for the SPECIES_MASS_FRACTION case, the form of the
     * fickian diffusive flux is not the form that produces the sum of
     * the diffusive fluxes equal zero for the case of binary diffusion and
     * the velocity being the mass-averaged velocity.
     **************************************************************************/
{
  int w1, a, j, q, var;
  double tmp, avg_molec_weight, coeff_rho, rhoD, *phi_ptr;

  /*
   *  Get diffusivity and Jacobian dependence on the diffusivity
   */
  if (Diffusivity()) EH(-1, "Error in Diffusivity.");

  /*
   *  Add in rho or C depending upon species variable type
   */
  coeff_rho = 1.0;
  if (mp->Species_Var_Type == SPECIES_MASS_FRACTION) {
    coeff_rho = mp->density;
  } else if (mp->Species_Var_Type == SPECIES_MOLE_FRACTION) {
    for (w1 = 0, avg_molec_weight = 0.0; w1 < pd->Num_Species; w1++) {
      avg_molec_weight += mp->molecular_weight[w1] * fv->c[w];
    }
    coeff_rho = mp->density / avg_molec_weight;
  }
  rhoD = coeff_rho * mp->diffusivity[w];

  for (a = 0; a < VIM; a++) {
    st->diff_flux[w][a] = - rhoD * st->grad_Y[w][a];
  }

  if (af->Assemble_Jacobian) {
    var = MASS_FRACTION;
    phi_ptr = bf[var]->phi;
    for (a = 0; a < VIM; a++) {
      for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++) {
        tmp = - coeff_rho * st->grad_Y[w][a]
            * mp->d_diffusivity[w][MAX_VARIABLE_TYPES + w1];
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          st->d_diff_flux_dc[w][a][w1][j] = tmp * phi_ptr[j];
        }
      }
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        st->d_diff_flux_dc[w][a][w][j] -= rhoD * bf[var]->grad_phi[j][a];
      }
    }

    for (q = 0; q < pd->Num_Dim; q++) {
      var = MESH_DISPLACEMENT1 + q;
      if (pd->v[pg->imtrx][var]) {
        for (a = 0; a < VIM; a++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            st->d_diff_flux_dmesh[w][a][q][j] =
                - rhoD * fv->d_grad_c_dmesh[a][w][q][j];
          }
        }
      }
    }

    var = TEMPERATURE;
    phi_ptr =  bf[var]->phi;
    for (a = 0; a < VIM; a++) {
      tmp = - coeff_rho * mp->d_diffusivity[w][var] * st->grad_Y[w][a];
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        st->d_diff_flux_dT[w][a][j] = tmp * phi_ptr[j];
      }
    }
  }
  return (0);
} /* END of routine fickian_flux */
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

int
generalized_fickian_flux (struct Species_Conservation_Terms *st, int w)

    /***************************************************************************
     *
     *  Function that computes the Residuals & Jacobian sensitivities of 
     *  Generalized Fickian diffusion terms
     *
     *   Flux - j(i) = coeff_rho * sum(D(ij) * grad_Y(j)) or
     *    {J} = coeff_rho [D] {grad_Y}
     *
     *   D(ij) is the mutual/binary diffusivity matrix.  It is known to
     *          have temperature and composition dependency
     *
     *   Generalized Fickian is a explicit formulation of diffusive flux.
     *   Two types of diffusivity calculator is available
     *    (1) constant D(ij)
     *    (2) D(ij) as a function of activity and free_volume diffusivity, 
     *        i.e.,
     *
     *        D(ij) = rho(i)*D(ii)*(d_mu_i/d_rho_j)/RT
     *        where D(ii) is the FV self-diffusivity of component i.
     *
     *     Author: A. C. Sun 4/00 (well...I did clone it)
     *
     *   If species var type is, the flux has units of:
     *        SPECIES_MOLE_FRACTION    -> mol cm-2 sec-1
     *            coeff_rho = C_avg
     *        SPECIES_MASS_FRACTION    -> gm  cm-2 sec-1
     *            coeff_rho = rho
     *        SPECIES_DENSITY          -> gm  cm-2 sec-1
     *            coeff_rho = 1.0
     *        SPECIES_UNDEFINED_FORM   -> nominally cm sec-1 (if diffusivity
     *                                        is in cm^2 sec-1)
     *            coeff_rho = 1.0
     *           (However, units are left up to the user, in practice)
     *        SPECIES_CONCENTRATION    -> mol cm-2 sec-1
     *            coeff_rho = 1.0
     **************************************************************************/
{
  int w1, w2, a, j, q, var, status = 0;
  double tmp, avg_molec_weight, coeff_rho, *phi_ptr;
  /*
   * Get diffusivity and derivatives 
   */
  if (Generalized_Diffusivity())  EH(-1, "Error in Diffusivity.");

  /*
   *  Add in rho or C depending upon species variable type
   */
  coeff_rho = 1.0;
  if (mp->Species_Var_Type == SPECIES_MASS_FRACTION) {
    coeff_rho = mp->density;
  } else if (mp->Species_Var_Type == SPECIES_MOLE_FRACTION) {
    for (w1 = 0, avg_molec_weight = 0.0; w1 < pd->Num_Species; w1++) {
      avg_molec_weight += mp->molecular_weight[w1] * fv->c[w];
    }
    coeff_rho = mp->density / avg_molec_weight;
  }

  for (a = 0; a < VIM; a++) {
    for (w1 = 0; w1 < pd->Num_Species; w1++) {
      st->diff_flux[w][a] -= 
	  coeff_rho * mp->diffusivity_gen_fick[w][w1] * st->grad_Y[w1][a];
    }
  }

  if (af->Assemble_Jacobian) { 
    var = MASS_FRACTION;
    phi_ptr =  bf[var]->phi;
    for (a = 0; a < VIM; a++) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
	tmp = coeff_rho * phi_ptr[j];
	for (w1 = 0; w1 < pd->Num_Species; w1++) {
	  for (w2 = 0; w2 < pd->Num_Species; w2++) {
	    st->d_diff_flux_dc[w][a] [w1][j] -= 
		mp->d_diffusivity_gf[w][w2][MAX_VARIABLE_TYPES + w1]
		* tmp * st->grad_Y[w2][a];
	  }
	  st->d_diff_flux_dc[w][a] [w1][j] -= 
	      mp->diffusivity_gen_fick[w][w1] * bf[var]->grad_phi[j][a];
	}
      }
    }
      
    for (q = 0; q < pd->Num_Dim; q++) {
      var = MESH_DISPLACEMENT1 + q;
      if (pd->v[pg->imtrx][var]) {
	for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
	  for (a = 0; a < VIM; a++) {
	    for (w1 = 0; w1 < pd->Num_Species; w1++) {
	      st->d_diff_flux_dmesh[w][a] [q][j] -=
		  coeff_rho * mp->diffusivity_gen_fick[w][w1]
		  * fv->d_grad_c_dmesh[a][w1][q][j];
	    }
	  }
	}
      }
    }
      
    var = TEMPERATURE;
    phi_ptr = bf[var]->phi;
    for (a = 0; a < VIM; a++) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
	tmp = coeff_rho * phi_ptr[j];
	st->d_diff_flux_dT[w][a] [j] = 0.;
	for (w1 = 0; w1 < pd->Num_Species; w1++) { 
	  st->d_diff_flux_dT[w][a][j] -=  
	      mp->d_diffusivity_gf[w][w1][var] 
	      * tmp * st->grad_Y[w1][a];
	}
      }
    }
  }
  return (status);
} /* END of routine generalized_fickian_flux */
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
/* assemble_invariant--
 *   
 *  assemble terms (Residual &| Jacobian) to generalized shear rate relation
 *
 *          shear_rate = ( 1/2 gamma_dot:gamma_dot)**1/2
 *
 *    where gamma_dot is the rate-of-deformation tensor
 *
 *    
 *
 * in:
 * 	ei -- pointer to Element Indices	structure
 *	pd -- pointer to Problem Description	structure
 *	af -- pointer to Action Flag		structure
 *	bf -- pointer to Basis Function		structure
 *	fv -- pointer to Field Variable		structure
 *	cr -- pointer to Constitutive Relation	structure
 *	md -- pointer to Mesh Derivative	structure
 *	me -- pointer to Material Entity	structure
 *
 * out:
 *	lec -- gets loaded up with local contributions to resid, Jacobian
 *
 * Created:	Apr 1, 1997 tabaer@sandia.gov
 *
 * Revised:	
 *
 *
 *
 */

int
assemble_invariant ( double tt,	/* parameter to vary time integration from 
				 * explicit (tt = 1) to implicit (tt = 0)    */
                     double dt ) /*  time step size                          */
{
  int dim;
  int p, a;
  
  int eqn, var;
  int peqn , pvar;
  int i, j, k, l;
  int status;
  
  dbl I2;
  dbl gd;
  dbl d_I2_dv[DIM][MDE];
  dbl d_I2_dmesh[DIM][MDE];
  
  dbl h3;		        	/* Volume element (scale factors). */
  dbl dh3dmesh_pj;	        	/* Sensitivity to (p,j) mesh dof. */
  
  dbl det_J;                            /* determinant of element Jacobian */
  
  dbl d_det_J_dmesh_pj;			/* for specific (p,j) mesh dof */
  
  dbl advection;	
  dbl source;
  dbl advection_a, advection_b;

  dbl diffusion=0;
  dbl diffusion_a, diffusion_b, diffusion_c;
  
 
  /*
   * Galerkin weighting functions for i-th and a-th momentum residuals 
   * and some of their derivatives...
   */
  
  dbl wt_func;
   
  
  /*
   * Interpolation functions for variables and some of their derivatives.
   */
  
  dbl phi_j;
  
  dbl wt;
    
  status = 0;
  
  /*
   * Unpack variables from structures for local convenience...
   */
  
  dim   = pd->Num_Dim;
  
 
  /*
   * Bail out fast if there's nothing to do...
   */
  
  if ( ! pd->e[pg->imtrx][eqn = R_SHEAR_RATE] )
    {
      return(status);
    }
  
  peqn = upd->ep[pg->imtrx][eqn];

  wt = fv->wt;                 /* Numerical integration weight */
  
  det_J = bf[eqn]->detJ;		/* Really, ought to be mesh eqn. */
  
  h3 = fv->h3;			/* Differential volume element (scales). */

  I2 = 0.0;

  /* 
   * Compute the 2nd invariant of the rate of strain tensor, I2 
   */

  for( k=0; k<VIM; k++)
    {
      for( l=0; l<VIM; l++)
	{
	  I2 += fv->grad_v[k][l]*fv->grad_v[k][l];
	  I2 += 2.0*fv->grad_v[k][l]*fv->grad_v[l][k];
	  I2 += fv->grad_v[l][k]*fv->grad_v[l][k];
	}
    }

  /* 
   * Typically, this function is called the "shear rate" and labeled gamma dot, 
   * Hence, gd for "gamma dot"
   */ 

  if(I2)
    {
      gd = pow(0.5*I2,0.5);
    }
  else
    {
      gd = 0.;
    }
  /*
   * Compute derivatives of 2nd invariant wrt velocities and displacements
   *
   */

  memset(d_I2_dv, 0 , sizeof(dbl)*DIM*MDE);
  memset(d_I2_dmesh, 0 , sizeof(dbl)*DIM*MDE);
  
  
  for(k=0; k<VIM; k++)
    {
      for(l=0; l<VIM; l++)
	{
	  for(a=0; a<VIM; a++)
	    
	    {
	      var = VELOCITY1+a;
	      
	      for(j=0; pd->v[pg->imtrx][var] && j<ei[pg->imtrx]->dof[var]; j++)
		{
		  d_I2_dv[a][j] += 2.0*(fv->grad_v[k][l]+fv->grad_v[l][k])* 
		    (bf[var]->grad_phi_e[j][a][k][l]+bf[var]->grad_phi_e[j][a][l][k]);
		}
	    }
	  
	  var = MESH_DISPLACEMENT1+a;
	  
	  for(j=0; pd->v[pg->imtrx][var] && j<ei[pg->imtrx]->dof[var]; j++)
	    {
	      d_I2_dmesh[a][j] += 2.0*(fv->grad_v[k][l]+fv->grad_v[l][k])*
		(fv->d_grad_v_dmesh[k][l][a][j] +
		 fv->d_grad_v_dmesh[l][k][a][j]);
	    }
	}
    }


  /*
   * Residuals_________________________________________________________________
   */
  
  if ( af->Assemble_Residual )
    {
      /*
       * Assemble the second_invariant equation
       */

	      for ( i=0; i<ei[pg->imtrx]->dof[eqn]; i++)
		{
		  
		  wt_func = bf[eqn]->phi[i];  
		  
		  advection = 0.;
		  
		  if ( pd->e[pg->imtrx][eqn] & T_ADVECTION )
		    {
		      advection = -gd;
		      advection *= wt_func * det_J * wt * h3;
		      advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
		    }

		  /*
		   * Diffusion Term..
		   */

		  /* OK this really isn't a diffusion term.  Its really a 
		   * filtering term.  But it looks like a diffusion operator.No?
		   */

		  diffusion = 0.;

		  if( pd->e[pg->imtrx][eqn] & T_DIFFUSION )
		    {
		      for( p=0; p < dim ; p++)
			{
			  diffusion += bf[eqn]->grad_phi[i][p]*fv->grad_SH[p];
			}

		      diffusion *= det_J*wt*h3;
		      diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
		    }
		  
		  /*
		   * Source term...
		   */
		  
		  source = 0;
		  
		  if ( pd->e[pg->imtrx][eqn] & T_SOURCE )
		    {
		      source += fv->SH;    
		      source *= wt_func * det_J * h3 * wt;
		      source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
		    }
		  
		  lec->R[peqn][i] += 
		    advection  + source +diffusion;      
  
		}
    }

  /*
   * Jacobian terms_________________________________________________________________
   */
  
  if ( af->Assemble_Jacobian )
    {

      for ( i=0; i<ei[pg->imtrx]->dof[eqn]; i++)
	{
	  wt_func = bf[eqn]->phi[i];


	  /*
	   * J_sh_v
	   */
	  
	  for( a=0; a<VIM; a++)
	    {
	      var = VELOCITY1+a;

	      if( pd->v[pg->imtrx][var])
		{
		  pvar = upd->vp[pg->imtrx][var];
		  for( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		    {

		      phi_j = bf[var]->phi[j];	 
     
		      advection = 0.;
			      
		      if ( pd->e[pg->imtrx][eqn] & T_ADVECTION )
			{
			  if(gd != 0.)
			    advection = -0.25*d_I2_dv[a][j]/gd;
			  else
			    advection = 0.;

			  advection *= wt_func * det_J * wt *h3;
			  advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
				  
			}
			      
		      lec->J[peqn][pvar][i][j] +=  advection;
		    }
		}
	    }

	  /*
	   * J_sh_SH
	   */
	  
	  var = SHEAR_RATE;
	  pvar = upd->vp[pg->imtrx][var];
	  for( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	    {

	      phi_j = bf[var]->phi[j];	 

	      advection = 0.0;
	      diffusion = 0.0;

	      if( pd->e[pg->imtrx][eqn] & T_DIFFUSION )
		{
		  for( p=0; p<dim; p++)
		    {
		      diffusion += bf[eqn]->grad_phi[i][p]*bf[var]->grad_phi[j][p];
		    }
		  
		  diffusion *= det_J*wt*h3;
		  diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
		}
	      
	      source = 0.0;

		if ( pd->e[pg->imtrx][eqn] & T_SOURCE )
		  {
		    source += phi_j;
		    source *= wt_func*det_J*wt*h3;
		    source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
		  }

	      lec->J[peqn][pvar][i][j] += advection + source + diffusion ;
	    }

	  /*
	   * J_sh_d
	   */
	  
	  for( p=0; p<dim; p++)
	    {
	      var = MESH_DISPLACEMENT1+p;

	      if( pd->v[pg->imtrx][var])
		{
		  pvar = upd->vp[pg->imtrx][var];
		  for( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];

		      d_det_J_dmesh_pj = bf[eqn]->d_det_J_dm[p][j];
			      
		      dh3dmesh_pj = fv->dh3dq[p] *phi_j;
     
		      advection = 0.;
		      advection_a = 0.;  /* d[gd]/dmeshb,j* |J|*h3 */
		      advection_b = 0.; /* gd*(dh3/dmesh*|J| + d|j|/dmesh*h3 ) */

		      if ( pd->e[pg->imtrx][eqn] & T_ADVECTION )
			{
			  if( gd != 0.)
			    advection_a = -0.25*d_I2_dmesh[p][j]*det_J*h3/gd;
			  else
			    advection_a = 0.;
			
			  advection_b -= gd;
			  advection_b *= d_det_J_dmesh_pj*h3 + det_J*dh3dmesh_pj;
			    
			  advection = advection_a + advection_b;
			
			  advection *= wt_func  * wt * pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

			}

		      diffusion_a = 0. ; /* grad_phi_i*d(grad_gd)/dmesh,bj *|J|*h3 */
		      diffusion_b = 0. ; /* d(grad_phi_i)/dmesh,bj*grad_gd *|J|*h3 */
		      diffusion_c = 0. ; /* grad_phi_i*grad_gd*( d|J|/dmesh,bj*h3 +|J|*dh3/dmesh,bj) */

		      if( pd->e[pg->imtrx][eqn] & T_DIFFUSION )
			{
			  for( a=0; a<dim ; a++)
			    {
			      diffusion_a += bf[eqn]->grad_phi[i][a]*fv->d_grad_SH_dmesh[a] [p][j];
			      diffusion_b += bf[eqn]->d_grad_phi_dmesh[i][a][p][j]*fv->grad_SH[a];
			      diffusion_c += bf[eqn]->grad_phi[i][a]*fv->grad_SH[a];
			    }

			  diffusion_a *= wt*det_J*h3;
			  diffusion_b *= wt*det_J*h3;
			  diffusion_c *= d_det_J_dmesh_pj*h3 + det_J*dh3dmesh_pj;

			  diffusion = diffusion_a + diffusion_b + diffusion_c;
			  diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
			  
			}
		      
			      
		      source    = 0.;

		      if ( pd->e[pg->imtrx][eqn] & T_SOURCE )
			{
			  source += fv->SH;
			  source *= d_det_J_dmesh_pj*h3 + det_J*dh3dmesh_pj;
			  source *=  wt_func*wt*pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
			}
			      
		      lec->J[peqn][pvar][i][j] += source + advection + diffusion;
		    }
		}
	    }
	} /* end of for(i.. loop */
    } /* end of if(af,, */

return(status);

 }  /* END of assemble_invariant */
/*******************************************************************************/
/* END of file mm_fill_species.c  */
/*******************************************************************************/


/* MMH
 *
 * I basically cloned this from get_convection_velocity, but used
 * the particle velocities instead of the fluid velocities. (fv->pv[]
 * instead of fv->v[]).
 * 
 * get_particle_convection_velocity()
 *
 * This routine calculates the effective convection velocity to be used
 * in the species transport and energy transport equations.
 * If the mesh motion is arbitrary, then the convection velocity is
 * v - xdot.  If the mesh motion is Lagrangian (moving with solid), then
 * the convection velockin_bc_leakity is the sum of diffusion volume fluxes of 
 * non-solid species divided by the volume fraction of solids
 */
int
get_particle_convection_velocity(double pvconv[DIM],
				 double pvconv_old[DIM], /* previous time    */
				 CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_pvconv,
				 double dt,
				 double tt)
{
  /*
   * Some local variables for convenience...
   */
  
  int dim, p, q, b, var, i, w, w1;
  int species;			/* MMH species number for particle phase. */
  double volsolid, volsolid_old;
  
  dim = pd->Num_Dim;
  species = (int) mp->u_density[0];

  for (p=0; p<VIM; p++)
    {
      pvconv[p] = 0.;
    }

   /* initialize sensitivity arrays if necessary */
   if ( af->Assemble_Jacobian )
     {
       var = MESH_DISPLACEMENT1;
       if (pd->v[pg->imtrx][var])
	 {
	   memset(d_pvconv->X, 0, DIM*DIM*MDE*sizeof(dbl));
	 }

       var = PVELOCITY1;
       if (pd->v[pg->imtrx][var])
	 {
	   memset(d_pvconv->v, 0, DIM*DIM*MDE*sizeof(dbl));
	 }

       var = MASS_FRACTION;
       if (pd->v[pg->imtrx][var])
	 {
	   memset(d_pvconv->C, 0, DIM*MAX_CONC*MDE*sizeof(dbl));
	 }

       var = TEMPERATURE;
       if (pd->v[pg->imtrx][var])
	 {
	   memset(d_pvconv->T, 0, DIM*MDE*sizeof(dbl));
	 }
     }

   if ( cr->MeshMotion == ARBITRARY )
     {
       /* calculate v - xdot and it's derivatives */

       for (p=0; p<VIM; p++)
	 {
	   pvconv[p] = fv->pv[p];
	   pvconv_old[p] = 0.;
	   if ( pd->TimeIntegration != STEADY )
	     {
	       pvconv_old[p] = fv_old->pv[p];
	       if (pd->v[pg->imtrx][R_MESH1]) {
		 pvconv[p] -= fv_dot->x[p];
		 pvconv_old[p] -= fv_dot_old->x[p];
	       }
	     }
	 }
       
       /* MMH
	* Sensitivities of particle phase wrt particle phase.
	*/
       if ( af->Assemble_Jacobian )
	 {
	   for (p=0; p<VIM; p++)
	     {
	       var = PVELOCITY1+p;
	       if (pd->v[pg->imtrx][var])
		 {
		   for ( i=0; i<ei[pg->imtrx]->dof[var]; i++)
		     {
		       d_pvconv->v[p][p][i] = bf[var]->phi[i];
		     }
		 }
	       
	       if ( pd->TimeIntegration != STEADY )
		 {
		   var = MESH_DISPLACEMENT1+p;
		   if (pd->v[pg->imtrx][var])
		     { 
		       for ( i=0; i<ei[pg->imtrx]->dof[var]; i++)
			 {
			   d_pvconv->X[p][p][i] -= bf[var]->phi[i] * (1 +2.* tt) /dt;
			 }
		     }
		 }
	     }
	 }
     } /* end of if ARBITRARY */
   
   if ( cr->MeshMotion == LAGRANGIAN ||
	cr->MeshMotion == DYNAMIC_LAGRANGIAN)
     {
       /*
	* in a lagrangian mesh, the convective velocity has two contributions
        * which are different in continuous and porous systems.  If the material
        * is continuous than the convective velocity is:
        *    v_conv  =  v_solvent_flux  +  v_solid_lagrangian_convection
        * If the material is porous, than the convective velocity is:
        *    v_conv  =  d_node_position/dt  +  v_solid_lagrangian_convection
	*
	* In both of these cases, v_solid_lagrangian_convection is the velocity of
        * a moving solid in a fixed frame of reference where the velocity of the
        * material is know if it is unstressed - then v_solid_lagrangian_convection
	* is a mapping of the know velocity of the stress free state to the current
	* state
	*
	* The following sections add these contributions into the convective velocity
	*/

       /* calculate (sum diffusion flux)/solid volume frac. and it's derivatives 
	*  v_solvent_flux  */
/*       if (pd->v[pg->imtrx][MASS_FRACTION] && mp->PorousMediaType == CONTINUOUS)   */
       if (pd->v[pg->imtrx][MASS_FRACTION] )  
	 /* if no solvents - then no convective velocity */
	 {
	   /* calculate volume fraction of solids and total diffusion flux */
	   volsolid = 1.;
	   volsolid_old = 1.;
	   for (w=0; w<pd->Num_Species_Eqn; w++)
	     {
/* MMH */
	       if(w != species)
		 {
		   volsolid -= fv->c[w];
		   if ( pd->TimeIntegration != STEADY )
		     volsolid_old -= fv_old->c[w];
		 }
	     }
	   
	   for (p=0; p<VIM; p++)
	     {
	       pvconv[p] = 0.;
	       pvconv_old[p] = 0.;
	       if ( cr->MassFluxModel == FICKIAN )
		 {
		   if ( Diffusivity() )  EH( -1, "Error in Diffusivity.");
		     /* Get diffusivity and derivatives */

/* MMH I don't think this matters for particle phase, there is no diffusion.
 * diffusivity[species] should be zero
 */
		   for (w=0; w<pd->Num_Species_Eqn; w++)
		     {
		       pvconv[p] -= mp->diffusivity[w] * fv->grad_c[w][p];
		       if ( pd->TimeIntegration != STEADY )
			 pvconv_old[p] -= mp->diffusivity[w] * fv_old->grad_c[w][p];
		     }
		 }
	       else  if ( cr->MassFluxModel == DARCY )
		 { /* diffusion induced convection is zero */
		 }
	       else
		 {
		   EH( -1, "Unimplemented mass flux constitutive relation.");
		 }
	     }
	   
	   for (p=0; p<VIM; p++)
	     {
	       pvconv[p] /= volsolid;
	       if ( pd->TimeIntegration != STEADY )
		 pvconv_old[p] /= volsolid_old;
	     }
	   if ( af->Assemble_Jacobian )
	     {
	       for (p=0; p<dim; p++)
		 {
		   var = MESH_DISPLACEMENT1;
		   if (pd->v[pg->imtrx][var])
		     {
		       for (q=0; q<dim; q++)
			 {
			   for ( i=0; i<ei[pg->imtrx]->dof[var]; i++)
			     {
			       if ( cr->MassFluxModel == FICKIAN )
				 {
				   if ( Diffusivity() )  EH( -1, "Error in Diffusivity.");
				      
				   for (w=0; w<pd->Num_Species_Eqn; w++)
				     {
				       d_pvconv->X[p][q][i] -= mp->diffusivity[w] *
					 fv->d_grad_c_dmesh[p][w] [q][i] / volsolid;
				     }
				 }
			     }
			 }
		     }
		   
		   /* Temperature dependence of diffusivity not yet implemented  */
		   /*
		   var = TEMPERATURE;
		   if (pd->v[pg->imtrx][var])
		     {
		       for ( i=0; i<ei[pg->imtrx]->dof[var]; i++)
			 {
			   if ( cr->MassFluxModel == FICKIAN )
			     {
			       if ( Diffusivity() )  EH( -1, "Error in Diffusivity.");
			       for (w=0; w<pd->Num_Species_Eqn; w++)
				 {
				   d_pvconv->T[p][i] -= mp->d_diffusivity[w][var][0] * bf[var]->phi[i] *
				     fv->grad_c[w][p] / volsolid;
				 }
			     }
			 }
		     }
		   */
		   
		   var = MASS_FRACTION;
		   if (pd->v[pg->imtrx][var])
		     {
		       if ( cr->MassFluxModel == FICKIAN )
			 {
			   if ( Diffusivity() )  EH( -1, "Error in Diffusivity.");

			   for (w=0; w<pd->Num_Species_Eqn; w++)
			     {
			       for ( i=0; i<ei[pg->imtrx]->dof[var]; i++)
				 {
				   d_pvconv->C[p][w][i] -=
				     mp->diffusivity[w] * bf[var]->grad_phi[i][p]
				       / volsolid;
				   for (w1=0; w1<pd->Num_Species_Eqn; w1++)
				     {
				       d_pvconv->C[p][w][i] -= (
							      mp->d_diffusivity[w][MAX_VARIABLE_TYPES + w1] 
							      * fv->grad_c[w][p] 
							      + mp->diffusivity[w] *fv->grad_c[w][p] / volsolid )
					 * bf[var]->phi[i] / volsolid;
				     }
				 }
			     }
			 }
		     }
		 } /* end of loop over pvconv directions */
	     } /* end of if Assemble Jacobian */
	 } /* end of if MASS_FRACTION */
       
        /*    d_node_position/dt */
       if (pd->TimeIntegration != STEADY && 
	   (mp->PorousMediaType == POROUS_UNSATURATED || 
	    mp->PorousMediaType == POROUS_SATURATED || 
	    mp->PorousMediaType == POROUS_TWO_PHASE))  
	 {
	   for (p=0; p<VIM; p++)
	     {
	       pvconv[p] = fv_dot->x[p];
	       pvconv_old[p] = fv_dot_old->x[p];
	     }
	   if ( af->Assemble_Jacobian )
	     {
	       for (p=0; p<VIM; p++)
		 {
		   
		   if ( pd->TimeIntegration != STEADY )
		     {
		       var = MESH_DISPLACEMENT1+p;
		       if (pd->v[pg->imtrx][var])
			 { 
			   for ( i=0; i<ei[pg->imtrx]->dof[var]; i++)
			     {
			       d_pvconv->X[p][p][i] += bf[var]->phi[i] * (1 +2.* tt) /dt;
			     }
			 }
		     }
		 }
	     }
	 } /* end of if POROUS media */
       /*
	* Add in convection due to motion of Stress Free State - Pseudo Lagrangian Convection 
	*  NOTE: this formulation still assumes mesh is quasi -   static, i.e. that mesh motion
	*        in transient calculations contributes negligibly to the momentum equation
        *    v_solid_lagrangian_convection
	*/

       /* First test to see what type of prescribed kinematics model */

       if (elc->v_mesh_sfs_model == CONSTANT)
	 {
           /*do nothing??*/
	 }
       else if (elc->v_mesh_sfs_model == ROTATIONAL ||
		elc->v_mesh_sfs_model == ROTATIONAL_3D )
	 {
           (void) V_mesh_sfs_model(elc->u_v_mesh_sfs, elc->v_mesh_sfs, 
				elc->v_mesh_sfs_model, -1);
	 }

       if ( pd->MeshInertia == 1)
	 {
	   if ( pd->TimeIntegration != STEADY ) EH(-1, "Can't have Unsteady Mesh Inertia ");
	   /*
	    * Velocity of solid in lab coordinates is the velocity of the stress free state
	    * dotted into the deformation gradient tensor
	    */
	   if (! pd->v[pg->imtrx][MESH_DISPLACEMENT1])
	     {
	       for (p=0; p < dim; p++)
		 {
		   pvconv[p] += elc->v_mesh_sfs[p];
		 }
	     }
	   else 
	     {
	       for (p=0; p < dim; p++)
		 {
		   for (q=0; q < dim; q++)
		     {
		       pvconv[p] += elc->v_mesh_sfs[q] * fv->deform_grad[q][p];
		     }
		 }
	       if ( af->Assemble_Jacobian )
		 {
		   for (b=0; b < dim; b++)
		     {
		       var = MESH_DISPLACEMENT1 + b;
		       if (pd->v[pg->imtrx][var])
			 {
			   for ( i=0; i<ei[pg->imtrx]->dof[var]; i++)
			     {
			       for (p=0; p < dim; p++)
				 {
				   for (q=0; q < dim; q++)
				     {
				       d_pvconv->X[p][b][i] += elc->v_mesh_sfs[q]
					 * fv->d_deform_grad_dx[q][p] [b][i];
				     }
				 }
			     }
			 }
		     }
		 }
	     }
	 }
       
     } /* end of if LAGRANGIAN */
   
   return 0;
 }
 
 
int
ls_modulate_speciessource(int w,
			  struct Species_Conservation_Terms *st)
{
  int i, b = -1, var;
  int dim = pd->Num_Dim;
  double factor;
  double pm_minus, pm_plus, width, f1, f2;
  struct Level_Set_Data *ls_old = ls;
	
 if(  mp->mp2nd == NULL ||
      mp->mp2nd->SpeciesSourceModel[w] != CONSTANT ||
      (!pd->gv[R_MASS] || !(pd->e[upd->matrix_index[R_MASS]][R_MASS] & T_SOURCE)) ) return(0);
	  
/* kludge for solidification tracking with phase function 0 */
 if(pfd != NULL && pd->gv[R_EXT_VELOCITY])
	{
	ls_old = ls;
	ls = pfd->ls[0];
 	width = ls->Length_Scale;
 	pm_minus = mp->mp2nd->speciessourcemask[0][w];	  
 	pm_plus =  mp->mp2nd->speciessourcemask[1][w];	  
 
 	f1 = st->MassSource[w];
 	f2 = mp->mp2nd->speciessource_phase[0][w];
 
 	if(f2 == 0.0 ) f2 = DBL_SMALL;
 	st->MassSource[w] = ls_modulate_property(  f1,
					f2,
					width,
					pm_minus,
					pm_plus,
					st->d_MassSource_dF[w],
					&factor);
	ls = ls_old;
	}
 width = ls->Length_Scale;
 pm_minus = mp->mp2nd->speciessourcemask[0][w];	  
 pm_plus =  mp->mp2nd->speciessourcemask[1][w];	  
 
 f1 = st->MassSource[w];
 f2 = mp->mp2nd->speciessource[w];
 
 if(f2 == 0.0 ) f2 = DBL_SMALL;
 
 st->MassSource[w] = ls_modulate_property(  f1,
					f2,
					width,
					pm_minus,
					pm_plus,
					st->d_MassSource_dF[w],
					&factor);

	  if( pd->v[pg->imtrx][var=TEMPERATURE] )
	    {
	      for ( i=0; i<ei[pg->imtrx]->dof[var]; i++)
		{
		  st->d_MassSource_dT[w][i] *=factor;
		}
	    }

	  if(  pd->v[pg->imtrx][var=MESH_DISPLACEMENT1] )
	    {
	      for( b=0; b<dim; b++ )
		{
		  for( i=0 ; i<ei[pg->imtrx]->dof[var+b]; i++)
		    {
		      st->d_MassSource_dmesh[w][b][i] *= factor;
		    }
		}
	    }

	  if( pd->v[pg->imtrx][var=VELOCITY1] )
	    {
	      for( b=0; b<dim; b++ )
		{
		  for( i=0 ; i<ei[pg->imtrx]->dof[var]; i++)
		    {
		      st->d_MassSource_dv[w][b][i] *= factor;
		    }
		}
	    }

	  if(  pd->v[pg->imtrx][var=MASS_FRACTION] )
	    {
	      for( b=0; b<pd->Num_Species; b++ )
		{
		  for( i=0 ; i<ei[pg->imtrx]->dof[var]; i++)
		    {
		      st->d_MassSource_dc[w][b][i] *= factor;
		    }
		}
	    }
	
	  if( pd->v[pg->imtrx][var=VOLTAGE] )
	  {
		  for( i=0 ; i<ei[pg->imtrx]->dof[var+b]; i++)
		  {
		      st->d_MassSource_dV[w][i] *= factor;
		  }
		  
	  }	  
	
return(0);
}
							
/* END of routine get_particle_convection_velocity */


void 
fickian_charged_gradient_bc(double func[],
		  	    double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
			    int wspec,	/* species number of this boundary condition */
		 	    double value) /* Value for normal surface gradient component       */

     /*******************************************************************************
      *
      *  Function which will set the normal component of concentraction gradient to the specified value
      *
      *  ----------------------------------------------------------------------------
      *
      *    calls fickian_charged_flux
      *  ----------------------------------------------------------------------------
      *
      *******************************************************************************/
{  
        int w, w1, a,q,j,var;
	double coeff_rho, avg_molec_weight;
	double sign = -1.0;

	w = wspec;

	struct Species_Conservation_Terms st;

	zero_structure(&st, sizeof(struct Species_Conservation_Terms), 1);
 	
	fickian_charged_flux(&st,wspec);


 	/*
 	  *  Add in rho or C depending upon species variable type
  	 */
  	coeff_rho = 1.0;
 	if (mp->Species_Var_Type == SPECIES_MASS_FRACTION) 
	{
	  	coeff_rho = mp->density;
	} 
	else if (mp->Species_Var_Type == SPECIES_MOLE_FRACTION) 
	{
    	     for (w1 = 0, avg_molec_weight = 0.0; w1 < pd->Num_Species; w1++) {
      		avg_molec_weight += mp->molecular_weight[w1] * fv->c[w1];
             }
    	coeff_rho = mp->density / avg_molec_weight;
  	}

	func[0] = value;

	/* Add on Fickian diffusivity portion so that only electrostatic flux is left */

 	for (a = 0; a < VIM; a++) 
	{
 	  	st.diff_flux[w][a] += coeff_rho * mp->diffusivity[w] * fv->grad_c[w][a];

		func[0] += sign*fv->snormal[a]*st.diff_flux[w][a];
  	}



  if (af->Assemble_Jacobian) { 
   	 var = MASS_FRACTION;
   	 for (a = 0; a < VIM; a++) {
     	 for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
		st.d_diff_flux_dc[w][a][w][j] += coeff_rho * (mp->diffusivity[w] * bf[var]->grad_phi[j][a]);
		d_func[0][MAX_VARIABLE_TYPES+w][j] += sign*fv->snormal[a]*st.d_diff_flux_dc[w][a][w][j];
		}
	}
      
   	 for (q = 0; q < pd->Num_Dim; q++) {
     	 var = MESH_DISPLACEMENT1 + q;
      		if (pd->v[pg->imtrx][var]) {
		for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
	 	 for (a = 0; a < VIM; a++) {
	   	 st.d_diff_flux_dmesh[w][a] [q][j] += coeff_rho * mp->diffusivity[w]*fv->d_grad_c_dmesh[a][w] [q][j];
		 
		 d_func[0][var][j] += sign*fv->snormal[a]*st.d_diff_flux_dmesh[w][a] [q][j];
		 d_func[0][var][j] += sign*fv->dsnormal_dx[a][q][j]*st.diff_flux[w][a];
		 
		 }
	     }
           }
    	}

	if(pd->v[pg->imtrx][var=TEMPERATURE] ){
  	 for (a = 0; a < VIM; a++) {
     	 for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
	   d_func[0][var][j] += sign*fv->snormal[a]*st.d_diff_flux_dT[w][a][j];
	  }
	 }	
	}

	if(pd->v[pg->imtrx][var=VOLTAGE] ){
  	 for (a = 0; a < VIM; a++) {
     	 for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
	   d_func[0][var][j] += sign*fv->snormal[a]*st.d_diff_flux_dV[w][a][j];
	  }
	 }	
	}

      }


  
  return;
}

/* End of fickian_charged_surf_gradient_bc */
