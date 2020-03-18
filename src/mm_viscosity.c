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

/* This routine contains all the viscosity models along with the routine 
 * that calls them. If you add a new viscosity model, it would go in here.
 */

/* Standard include files */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* GOMA include files */

#include "mm_viscosity.h"

#include "std.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "el_elm.h"
#include "rf_bc_const.h"
#include "rf_solver.h"
#include "mm_mp_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_mp_structs.h"
#include "mm_mp.h"
#include "mm_eh.h"
#include "mm_fill_ls.h"
#include "mm_fill_terms.h"
#include "mm_fill_util.h"
#include "rf_allo.h"
#include "user_mp.h"
#include "user_mp_gen.h"
#include "mm_qtensor_model.h"

#define GOMA_MM_VISCOSITY_C

/*********** R O U T I N E S   I N   T H I S   F I L E ***********************
 * 
 *
 *     viscosity()
 *      (then submodels for viscosity called by viscosity:)
 *        power_law_viscosity()
 *        herschel_buckley_viscosity()
 *        carreau_viscosity()
 *        bingham_viscosity()
 *        bingham_wlf_viscosity()
 *        carreau_wlf_viscosity()
 *        fill_viscosity()
 *        suspension_viscosity()
 *        carreau_suspension_viscosity(()
 *        powerlaw_suspension_viscosity()
 *        epoxy_viscosity()
 *        sylgard_viscosity(()
 *        filled_epoxy_viscosity()
 *        foam_epoxy_viscosity()
 *        thermal_viscosity()
 *        cure_viscosity()
 *        bond_viscosity()
 *        carreau_wlf_conc_viscosity()
 *        ls_modulate_viscosity()
 *     copy_pF_to_F()
 */
/*******************************************************************************
 * viscosity(): Calculate the viscosity and derivatives of viscosity
 *              with respect to solution unknowns at the Gauss point. Most 
 *              non-Newtonian purely viscous models depend on the shear rate
 *              invariant and may also depend on temperature and pressure.
 *
 * Input
 *----------
 *   gn_local	= the Generalized_Newtonian structure for solvent or polymer
 *   gamma_dot	= Strain rate tensor
 *
 * Output
 * -----
 *
 *   mu		= viscosity
 *   d_mu       = dependence of viscosity on the indendent unknowns in the
 *                local element stiffness matrix, where:
 *
 *   d_mu->gd	= derivative of viscosity wrt to strain rate invariant variables 
 *                (which one?)
 *   d_mu->v	= derivative of viscosity wrt to velocity variables
 *   d_mu->X	= derivative of viscosity wrt to mesh variables
 *   d_mu->T	= derivative of viscosity wrt to temperature variables
 *   d_mu->P	= derivative of viscosity wrt to pressure variables
 *   d_mu->C	= derivative of viscosity wrt to concentration/species variables
 *   d_mu->F	= derivative of viscosity wrt to FILL (Level Set / VoF) variables
 *   d_mu->nn	= derivative of viscosity wrt to bond concentration variables
 *
 *
 *******************************************************************************/
double
viscosity(struct Generalized_Newtonian *gn_local,
	  dbl gamma_dot[DIM][DIM],
          VISCOSITY_DEPENDENCE_STRUCT *d_mu)
{
  int err;
  int a;

  int var, var_offset;
  int v,w, w1;

  int vdofs;

  int i, j;

  int species;			/* Species number of particle phase. */
  dbl p_vol_frac;		/* local particle volume fraction. */
  dbl mu = 0.;
  int dim = ei[pg->imtrx]->ielem_dim;

  struct Level_Set_Data *ls_old;

  /* Zero out sensitivities */
  zeroStructures(d_mu, 1);
 
  /* this section is for all Newtonian models */

  if (gn_local->ConstitutiveEquation == NEWTONIAN) 
    {
      if (mp->ViscosityModel == USER )
	{
	  err = usr_viscosity(mp->u_viscosity);
	  mu = mp->viscosity;
	  
	  var = TEMPERATURE;
	  if (d_mu != NULL && pd->v[pg->imtrx][var])
	    {
	      for (j=0; j<ei[pg->imtrx]->dof[var]; j++)
		{
		  d_mu->T[j]= mp->d_viscosity[var]*bf[var]->phi[j];
		}
	    }
	  
#ifdef COUPLED_FILL
	  var = FILL;
	  if (d_mu != NULL && pd->v[pg->imtrx][var])
	    {
	      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
		{
		  d_mu->F[j]= mp->d_viscosity[var] * bf[var]->phi[j];
		}
	    }
#endif /* COUPLED_FILL */

	  if ( d_mu != NULL && pd->v[pg->imtrx][MESH_DISPLACEMENT1] )
	    {
	      for ( a=0; a<dim; a++)
		{
		  var = MESH_DISPLACEMENT1 + a;
		  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		    {
		      d_mu->X[a][j] =mp->d_viscosity[var]*bf[var]->phi[j];
		    }
		}
	    }

	  if ( d_mu != NULL && pd->v[pg->imtrx][MASS_FRACTION] )
	    {
	      for ( w=0; w<pd->Num_Species_Eqn; w++)
		{
		  var	     = MASS_FRACTION;
		  var_offset = MAX_VARIABLE_TYPES + w;
		  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		    {
		      d_mu->C[w][j] =mp->d_viscosity[var_offset]*bf[var]->phi[j];
		    }
		}
	    }
	    
	  v	= VELOCITY1;
	  vdofs	= ei[pg->imtrx]->dof[v];
	  if ( d_mu != NULL && pd->v[pg->imtrx][v] )
	    {
	      for ( a=0; a<dim; a++)
		{
		  for ( j=0; j<vdofs; j++)
		    {
		      d_mu->v[a][j] = mp->d_viscosity[v+a]*bf[v]->phi[j];
		    }
		}
	    }

	  v	= PRESSURE;
	  vdofs	= ei[pg->imtrx]->dof[v];
	  if ( d_mu != NULL && pd->v[pg->imtrx][v] )
	    {
	      for ( j=0; j<vdofs; j++)
		{
		  d_mu->P[j] = mp->d_viscosity[v]*bf[v]->phi[j];
		}
	    }
	}
      else if(mp->ViscosityModel == USER_GEN)
	{
	  err = usr_viscosity_gen(&mu, gamma_dot,
				  &d_mu->gd, d_mu->v, d_mu->X,
				  d_mu->T, d_mu->P, d_mu->C,
				  mp->u_viscosity);
	}
      else if (mp->ViscosityModel == FILL )
	{
	  /* let d_mu->F be zero here since mu is discontinuous */
	  err = fill_viscosity(mp->u_viscosity);
	  mu =  mp->viscosity;
	}
      else if (mp->ViscosityModel == LEVEL_SET )
	{
	  double mu0   = mp->u_viscosity[0];
	  double mu1   = mp->u_viscosity[1];
	  double width = mp->u_viscosity[2];
	  if ( d_mu != NULL )
            err = level_set_property(mu0, mu1, width, &mu, d_mu->F);
          else
            err = level_set_property(mu0, mu1, width, &mu, NULL);
	  EH(err, "level_set_property() failed for viscosity.");
	}
      else if (mp->ViscosityModel == LS_QUADRATIC )
	{
	  double mu0   = mp->u_viscosity[0];
	  double mu1   = mp->u_viscosity[1];
	  double width = mp->u_viscosity[2];
	  double rho0  = mp->u_density[0];
	  double rho1  = mp->u_density[1];
	  double nu0 = mu0/rho0;
	  double nu1 = mu1/rho1;

	  load_lsi(width);
	  
	  mu = 0.0;

	  mu += mu0;

	  mu += ( rho0*nu1 + nu0*rho1 - 2.0*mu0)*lsi->H;

	  mu += ( nu1 - nu0 ) * ( rho1 - rho0) * lsi->H * lsi->H;
	}
      else if (mp->ViscosityModel == CONST_PHASE_FUNCTION)
	{
	  // PRS: this is fubar'd    Needs fixing and documentation. This is cludged.  
	  double mu1, mu2, tmp_mu;
	  double width;
	  int a, num;
	  struct Level_Set_Data *ls_save = ls;
	  
	  num = pfd->num_phase_funcs;

	  width = mp->u_viscosity[num];

	  for (a = 0; a < num; a++)
	    {
	      mu1 = mp->u_viscosity[a]; 
	      mu2 = mp->u_viscosity[a+2];
	      ls = pfd->ls[a];
	      if( d_mu != NULL )
		err = level_set_property ( mu1, mu2, width, &tmp_mu, d_mu->pf[a] );
	      else
		err = level_set_property ( mu1, mu2, width, &tmp_mu, NULL );

	      mu += tmp_mu;
	    }
		
	 if ( fabs(mu) < DBL_SMALL) mu = mp->u_viscosity[num + 1 ];
	 ls = ls_save;
	}
	  
      /* MMH
       * Remember, this is only called for the fluid phase.
       */
      else if (mp->ViscosityModel == SUSPENSION_PM )
	{
	  species    = (int) mp->u_density[0];
	  p_vol_frac = fv->c[species];
	  mu	     = mp->viscosity* pow(1.0-p_vol_frac,-2.5);

          if ( d_mu != NULL )
            {
              for (j = 0; j < ei[pg->imtrx]->dof[MASS_FRACTION]; j++) {
                d_mu->C[species][j] = (mu*2.5/(1.0 - p_vol_frac)) *
                  bf[MASS_FRACTION]->phi[j];
	      }
            }
	}
      else if (mp->ViscosityModel == CONSTANT )
	{
	  /*  mu   = gn_local->mu0; corrected for auto continuation 3/01 */
	  if(gn_local->ConstitutiveEquation == CONSTANT)
	    {
	      mu = gn_local->mu0;
	    }
	  else
	    {
	      mu = mp->viscosity;
	    }
          mp_old->viscosity = mu;
	  /*Sensitivities were already set to zero */
	}
      else if (mp->ViscosityModel == TABLE)
	{
          struct  Data_Table *table_local;
          table_local = MP_Tables[mp->viscosity_tableid];
          apply_table_mp(&mp->viscosity, table_local);
	  mu = mp->viscosity;

          if ( d_mu != NULL )
            {
              for(i=0;i<table_local->columns-1;i++)
	        {
                  var = table_local->t_index[i];

	          switch (var)
		    {
		    case TEMPERATURE:
		      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
		        {
                          d_mu->T[j] = table_local->slope[i];

		        }
		      break;
		    case MASS_FRACTION:
		      if (pd->v[pg->imtrx][MASS_FRACTION] )
		        {
		          for ( w=0; w<pd->Num_Species_Eqn; w++)
			    {
			      var	     = MASS_FRACTION;
			      var_offset = MAX_VARIABLE_TYPES + w;
			      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
			        {
                                  d_mu->C[w][j] =table_local->slope[i];
			        }
			    }
		            }
		      break;
		    default:
		      EH(-1, "Variable function not yet implemented in material property table");
		    }
	        }
            }
	}
       else
	{
	  EH(-1,"Unrecognized viscosity model for Newtonian fluid");
	}

    } /* end Newtonian section */

  else if (gn_local->ConstitutiveEquation == CONSTANT)
    {
      mu   = gn_local->mu0;
      mp_old->viscosity = mu;
      /*Sensitivities were already set to zero */

    }
  /* this section is for all non-Newtonian models */
  else if (gn_local->ConstitutiveEquation == SUSPENSION)
    {
      if ( pd->v[pg->imtrx][SHELL_PARTC] )
        {
	  err = suspension_viscosity(gn_local->sus_species_no, gn_local->mu0,
				     gn_local->maxpack, gn_local->nexp, fv->sh_pc);
	  EH(err, "suspension_viscosity");
	  mu = mp->viscosity;
        }
      else
	{
	  err = suspension_viscosity(gn_local->sus_species_no, gn_local->mu0,
				     gn_local->maxpack, gn_local->nexp, fv->c[gn_local->sus_species_no]);
	
	  EH(err, "suspension_viscosity");
	  mu = mp->viscosity;
	}
      
      if ( d_mu != NULL && pd->v[pg->imtrx][MASS_FRACTION] )
	{
	  var	     = MASS_FRACTION;
	  var_offset = MAX_VARIABLE_TYPES + gn_local->sus_species_no;
	  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	    {
	      d_mu->C[gn_local->sus_species_no][j] =mp->d_viscosity[var_offset]*bf[var]->phi[j];
	    }
	}
    }
  else if (gn_local->ConstitutiveEquation == CURE)
    {
      err = cure_viscosity(gn_local->cure_species_no, gn_local->mu0,
			   gn_local->gelpoint, gn_local->cureaexp, gn_local->curebexp);
      EH(err, "cure_viscosity");
      
      mu = mp->viscosity;
      
      var = MASS_FRACTION;
      if (d_mu != NULL && pd->v[pg->imtrx][var] )
	{
	  w = gn_local->cure_species_no;
	  var_offset = MAX_VARIABLE_TYPES + w;
	  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	    {
	      d_mu->C[w][j] = mp->d_viscosity[var_offset]*bf[var]->phi[j];
	    }
	}
    }
  else if (gn_local->ConstitutiveEquation == THERMAL)
    {
      double mu0, atexp;
      mu0=gn_local->mu0; atexp= gn_local->atexp;
      if(gn_local->mu0Model == LEVEL_SET)
	{
	  if(d_mu != NULL)
	    {
	      err = level_set_property(gn_local->u_mu0[0], gn_local->u_mu0[1], 
				       gn_local->u_mu0[2], &mu0, d_mu->F);
	    }
	  else
	    {
	      err = level_set_property(gn_local->u_mu0[0], gn_local->u_mu0[1], 
				       gn_local->u_mu0[2], &mu0, NULL);
	    }
	}
	  
      if(gn_local->atexpModel == LEVEL_SET)
	{
	  if(d_mu != NULL)
	    {
	      err = level_set_property(gn_local->u_atexp[0], gn_local->u_atexp[1], 
				       gn_local->u_atexp[2], &atexp, d_mu->F);
	    }
	  else
	    {
	      err = level_set_property(gn_local->u_atexp[0], gn_local->u_atexp[1], 
				       gn_local->u_atexp[2], &atexp, NULL);
	    }
	}
	    
      err = thermal_viscosity(mu0, atexp);
      EH(err, "thermal_viscosity");
      
      mu = mp->viscosity;
      
      var = TEMPERATURE;
      if ( d_mu != NULL && pd->v[pg->imtrx][TEMPERATURE] )
	
	{
	  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	    {
	      d_mu->T[j]= mp->d_viscosity[var]*bf[var]->phi[j];
	    }
	}
    }
  else if (gn_local->ConstitutiveEquation == BOND)
    {
      err = bond_viscosity(gn_local->mu0, gn_local->muinf, gn_local->aexp);
      EH(err, "bond_viscosity");
      
      mu = mp->viscosity;
      
      var = BOND_EVOLUTION;
      if ( d_mu != NULL && pd->v[pg->imtrx][var] )
	
	{
	  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	    {
	      d_mu->nn[j]= mp->d_viscosity[var]*bf[var]->phi[j];
	    }
	}
    }
  else if (gn_local->ConstitutiveEquation == EPOXY)
    {
      err = epoxy_viscosity(gn_local->cure_species_no, gn_local->mu0,
			    gn_local->gelpoint, gn_local->cureaexp, 
			    gn_local->curebexp, gn_local->atexp);
      EH(err, "epoxy_viscosity");
      
      mu = mp->viscosity;
      
      var = TEMPERATURE;
      if ( d_mu != NULL && pd->v[pg->imtrx][TEMPERATURE] )
	{
	  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	    {
	      d_mu->T[j]= mp->d_viscosity[var]*bf[var]->phi[j];
	    }
	}
      
      var = MASS_FRACTION;
      if ( d_mu != NULL && pd->v[pg->imtrx][var] )
	{
	  w = gn_local->cure_species_no;
	  var_offset = MAX_VARIABLE_TYPES + w;
	  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	    {
	      d_mu->C[w][j] =mp->d_viscosity[var_offset]*bf[var]->phi[j];
	    }
	}
    }
   else if (gn_local->ConstitutiveEquation == FOAM_PMDI_10)
    {
      err = foam_pmdi10_viscosity(gn_local->cure_species_no, gn_local->mu0,
                            gn_local->gelpoint, gn_local->cureaexp,
                            gn_local->curebexp, gn_local->atexp);
      EH(err, "foam_pmdi10_viscosity");

      mu = mp->viscosity;

      var = TEMPERATURE;
      if ( d_mu != NULL && pd->v[pg->imtrx][TEMPERATURE] )
        {
          for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
            {
              d_mu->T[j]= mp->d_viscosity[var]*bf[var]->phi[j];
            }
        }

      var = MASS_FRACTION;
      if ( d_mu != NULL && pd->v[pg->imtrx][var] )
        {
          w = gn_local->cure_species_no;
          var_offset = MAX_VARIABLE_TYPES + w;
          for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
            {
              d_mu->C[w][j] =mp->d_viscosity[var_offset]*bf[var]->phi[j];
            }
        }
    }
  else if (gn_local->ConstitutiveEquation == SYLGARD)
    {
      err = sylgard_viscosity(gn_local->cure_species_no, gn_local->mu0,
			      gn_local->gelpoint, gn_local->cureaexp, 
			      gn_local->atexp);
      EH(err, "sylgard_viscosity");
      
      mu = mp->viscosity;
      
      var = TEMPERATURE;
      if ( d_mu != NULL && pd->v[pg->imtrx][TEMPERATURE] )
	{
	  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	    {
	      d_mu->T[j]= mp->d_viscosity[var]*bf[var]->phi[j];
	    }
	}
      
      var = MASS_FRACTION;
      if ( d_mu != NULL && pd->v[pg->imtrx][var] )
	{
	  w = gn_local->cure_species_no;
	  var_offset = MAX_VARIABLE_TYPES + w;
	  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	    {
	      d_mu->C[w][j] =mp->d_viscosity[var_offset]*bf[var]->phi[j];
	    }
	}
    }
  else if (gn_local->ConstitutiveEquation == FILLED_EPOXY )
    {
      err = filled_epoxy_viscosity(gn_local->sus_species_no, gn_local->cure_species_no, 
				   gn_local->mu0, gn_local->maxpack,
				   gn_local->nexp, gn_local->gelpoint, gn_local->cureaexp, 
				   gn_local->curebexp, gn_local->tgel0, gn_local->atexp);
      
      mu = mp->viscosity;
      if (d_mu != NULL) 
	{
	  if (pd->v[pg->imtrx][TEMPERATURE])
	    {
	      for (j = 0; j < ei[pg->imtrx]->dof[TEMPERATURE]; j++)
		{
		  d_mu->T[j] = mp->d_viscosity[TEMPERATURE] * bf[TEMPERATURE]->phi[j];
		}
	    }
      
	  var = MASS_FRACTION;
	  var_offset = MAX_VARIABLE_TYPES;
	  if (pd->v[pg->imtrx][var])
	    {
	      w  = gn_local->cure_species_no;
	      w1 = gn_local->sus_species_no;
	      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		{
		  d_mu->C[w][j]  = mp->d_viscosity[var_offset + w ]*bf[var]->phi[j];
		  d_mu->C[w1][j] = mp->d_viscosity[var_offset + w1]*bf[var]->phi[j];
		}
	    }
	}
    }
  else if (gn_local->ConstitutiveEquation == FOAM_EPOXY )
    {
      err = foam_epoxy_viscosity(gn_local->sus_species_no, gn_local->cure_species_no,
				 gn_local->mu0, gn_local->gelpoint, gn_local->atexp);
      mu = mp->viscosity;
      if (d_mu != NULL)
	{
	  if (pd->v[pg->imtrx][TEMPERATURE])
	    {
	      for (j = 0; j < ei[pg->imtrx]->dof[TEMPERATURE]; j++)
		{
		  d_mu->T[j] = mp->d_viscosity[TEMPERATURE] * bf[TEMPERATURE]->phi[j];
		}
	    }
	  if (pd->v[pg->imtrx][MASS_FRACTION])
	    {
	      for (w = 0; w < pd->Num_Species_Eqn; w++) 
		{
		  dbl tmp = mp->d_viscosity[MAX_VARIABLE_TYPES + w];
		  if (tmp != 0.0)
		    {
		      for (j = 0; j < ei[pg->imtrx]->dof[MASS_FRACTION]; j++)
			{
			  d_mu->C[w][j] = tmp * bf[MASS_FRACTION]->phi[j];
			}
		    }
		}
	    }
	}
    }
  else if (gn_local->ConstitutiveEquation == POWER_LAW) 
    {
      mu = power_law_viscosity(gn_local, gamma_dot, d_mu);
    }
  else if (gn_local->ConstitutiveEquation == CARREAU)
    {
      mu = carreau_viscosity(gn_local, gamma_dot, d_mu);
    }
  else if (gn_local->ConstitutiveEquation == BINGHAM)
    {
      mu = bingham_viscosity(gn_local, gamma_dot, d_mu);
    }
  else if (gn_local->ConstitutiveEquation == BINGHAM_WLF)
    {
      mu = bingham_wlf_viscosity(gn_local, gamma_dot, d_mu);
    }
  else if (gn_local->ConstitutiveEquation == CARREAU_WLF)
    {
      mu = carreau_wlf_viscosity(gn_local, gamma_dot, d_mu);
    }
  else if (gn_local->ConstitutiveEquation == CARREAU_SUSPENSION)
    {
      mu = carreau_suspension_viscosity(gn_local, gamma_dot, d_mu);
    }
  else if (gn_local->ConstitutiveEquation == POWERLAW_SUSPENSION)
    {
      mu = powerlaw_suspension_viscosity(gn_local, gamma_dot, d_mu);
    }

  else if (gn_local->ConstitutiveEquation == HERSCHEL_BULKLEY)
    {
      mu = herschel_buckley_viscosity(gn_local, gamma_dot, d_mu);
    }
  else if (gn_local->ConstitutiveEquation == CARREAU_WLF_CONC_PL ||
	   gn_local->ConstitutiveEquation == CARREAU_WLF_CONC_EXP)
    {
      mu = carreau_wlf_conc_viscosity(gn_local, gamma_dot, d_mu,
				      gn_local->ConstitutiveEquation);
    }
  else if (ls != NULL && gn_local->ConstitutiveEquation == VE_LEVEL_SET )
    {
      double pos_mup   = gn_local->pos_ls_mup;
      double neg_mup   = gn_local->mu0;
      double width     = ls->Length_Scale;
      if ( d_mu != NULL )
        err = level_set_property(neg_mup, pos_mup, width, &mu, d_mu->F);
      else
        err = level_set_property(neg_mup, pos_mup, width, &mu, NULL);
      EH(err, "level_set_property() failed for polymer viscosity.");
    }
  else
    {
      EH(-1,"Unrecognized viscosity model for non-Newtonian fluid");
    }
  
  if (ls != NULL && gn_local->ConstitutiveEquation != VE_LEVEL_SET &&
      mp->ViscosityModel != LEVEL_SET &&
      mp->ViscosityModel != LS_QUADRATIC &&
      mp->mp2nd != NULL &&
      (mp->mp2nd->ViscosityModel == CONSTANT || mp->mp2nd->ViscosityModel == RATIO)
    )
    {
      /* kludge for solidification tracking with phase function 0 */
      if (pfd != NULL && pd->e[pg->imtrx][R_EXT_VELOCITY])
	{
	  ls_old = ls;
	  ls = pfd->ls[0];
	  err = ls_modulate_viscosity(&mu, mp->mp2nd->viscosity_phase[0],
				      ls->Length_Scale,
				      (double) mp->mp2nd->viscositymask[0],
                                      (double) mp->mp2nd->viscositymask[1],
                                      d_mu, mp->mp2nd->ViscosityModel );
	  ls = ls_old;
	}
      err= ls_modulate_viscosity(&mu, mp->mp2nd->viscosity, ls->Length_Scale,
				 (double) mp->mp2nd->viscositymask[0],
                                 (double) mp->mp2nd->viscositymask[1],
                                 d_mu, mp->mp2nd->ViscosityModel );
      EH(err, "ls_modulate_viscosity");
    }
  if(DOUBLE_NONZERO(gn_local->thixo_factor))
    { 
     dbl thixotropic = 1;
     dbl thixo_time = 0;
     if(fv->restime > 0.)
	{ thixotropic += gn_local->thixo_factor*fv->restime; }
     if(tran->time_value > 0.)
	{ thixotropic += thixo_time*tran->time_value; }
     if ( d_mu != NULL )
       {
       d_mu->gd *= thixotropic;
       if ( pd->v[pg->imtrx][var=BOND_EVOLUTION] )
         {
          for( i=0 ; i<ei[pg->imtrx]->dof[var]; i++)
	    {
	     d_mu->nn[i] *= thixotropic;
	    }
         }
       var = RESTIME;
       if (pd->v[pg->imtrx][var] && fv->restime>0.)
        {
         for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	   {
	     d_mu->degrade[j] = mu*gn_local->thixo_factor*bf[var]->phi[j];
	   }
        }
       var = TEMPERATURE;
       if (pd->v[pg->imtrx][var] )
        {
         for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	   {
	     d_mu->T[j] *= thixotropic;
	   }
        }
       if (pd->e[pg->imtrx][R_MESH1] )
        {
         for ( i=0; i<VIM; i++)
	   {
	    for ( j=0; j< ei[pg->imtrx]->dof[R_MESH1]; j++)
	       {
		  d_mu->X [i][j] *= thixotropic;
	       }
	   }
        }
       if (pd->e[pg->imtrx][R_MOMENTUM1] )
        {
         for ( a=0; a<VIM; a++)
           {
            for ( i=0; i < ei[pg->imtrx]->dof[VELOCITY1]; i++)
	      {
	          d_mu->v[a][i] *= thixotropic;
	      }
	   }
	}
       if ( pd->v[pg->imtrx][var=MASS_FRACTION ] )
        {
         for ( w=0; w<pd->Num_Species_Eqn; w++)
	   {
	    for( i=0; i<ei[pg->imtrx]->dof[var]; i++) 
	      {
	       d_mu->C[w][i] *= thixotropic;
	      }
	   }
        }
       if( pd->v[pg->imtrx][var=PRESSURE] )
       {
        for( i=0; i<ei[pg->imtrx]->dof[var]; i++ )
	  {
	   d_mu->P[i] *= thixotropic;
	  }
       }
      }
     mu *= thixotropic;
    }
  return(mu);
}



double
power_law_viscosity(struct Generalized_Newtonian *gn_local,
		    dbl gamma_dot[DIM][DIM], /* strain rate tensor */
                    VISCOSITY_DEPENDENCE_STRUCT *d_mu )
{

  int a, b;
  int mdofs=0,vdofs;

  int i, j;

  dbl gammadot;	                /* strain rate invariant */

  dbl d_gd_dv[DIM][MDE];        /* derivative of strain rate invariant
				   wrt velocity */
  dbl d_gd_dmesh[DIM][MDE];     /* derivative of strain rate invariant
				   wrt mesh */

  dbl val;
  dbl mu0;
  dbl nexp;
  dbl offset;
  dbl mu = 0.;

  vdofs = ei[pg->imtrx]->dof[VELOCITY1];

  if ( pd->e[pg->imtrx][R_MESH1] )
    {
      mdofs = ei[pg->imtrx]->dof[R_MESH1];
    }


  calc_shearrate(&gammadot, gamma_dot, d_gd_dv, d_gd_dmesh);

  /* calculate power law viscosity
     mu = mu0 * (offset + gammadot)**(nexp-1))                 */
  mu0 = gn_local->mu0;
  nexp = gn_local->nexp;
  offset = 0.00001;

  val = pow( gammadot + offset, nexp-1.);
  mu = mu0 * val;
  /*
   * d( mu )/dmesh
   */

  val = pow( gammadot+offset, nexp-2.);


  if (d_mu != NULL)
    {
     d_mu->gd = mu0*(nexp - 1.0)*val;
    }

  if ( d_mu != NULL && pd->e[pg->imtrx][R_MESH1] )
    {
      for ( b=0; b<VIM; b++)
	{
	  for ( j=0; j<mdofs; j++)
	    {
	      if(gammadot != 0.0 && Include_Visc_Sens )
		{

		  d_mu->X [b][j] =
		    d_mu->gd * d_gd_dmesh [b][j] ;

		}
	      else
		{
		  /* printf("\ngammadot is zero in viscosity function");*/
		  d_mu->X [b][j] = 0.0;
		}
	    }
	}
    }
  
  /*
   * d( mu )/dv
   */
  if ( d_mu != NULL && pd->e[pg->imtrx][R_MOMENTUM1] )
    {
      for ( a=0; a<VIM; a++)
        {
          for ( i=0; i<vdofs; i++)
	    {
	      if(gammadot != 0.0 && Include_Visc_Sens )
	        {
	          d_mu->v[a][i] =
		    d_mu->gd * d_gd_dv[a][i] ;
	        }
	      else
	        {
	          d_mu->v[a][i] = 0.0 ;
	        }
	    }
        }
    }

  return(mu);
}

double
herschel_buckley_viscosity(struct Generalized_Newtonian *gn_local,
			   dbl gamma_dot[DIM][DIM], /* strain rate tensor */
			   VISCOSITY_DEPENDENCE_STRUCT *d_mu)
{
  int a, b;
  int mdofs=0,vdofs;

  int i, j;

  dbl gammadot;	                /* strain rate invariant */ 

  dbl d_gd_dv[DIM][MDE];        /* derivative of strain rate invariant 
				   wrt velocity */ 
  dbl d_gd_dmesh[DIM][MDE];     /* derivative of strain rate invariant 
				   wrt mesh */ 

  dbl val;
  dbl mu = 0.;
  dbl mu0;
  dbl nexp;
  dbl tau_y;
  dbl offset;

  vdofs = ei[pg->imtrx]->dof[VELOCITY1];
  
  if ( pd->e[pg->imtrx][R_MESH1] )
    {
      mdofs = ei[pg->imtrx]->dof[R_MESH1];
    }
  

  calc_shearrate(&gammadot, gamma_dot, d_gd_dv, d_gd_dmesh);

  /* calculate power law viscosity
     power law constants - make this migrate to input deck asap
     mu = mu0 * (offset + gammadot)**(nexp-1))                 */
  mu0 = gn_local->mu0;
  nexp = gn_local->nexp;
  tau_y = gn_local->tau_y;
  offset = 0.00001;
  
  val = pow( gammadot + offset, nexp-1.);
  mu = mu0 * val;
  mu += tau_y/(gammadot+offset);

  /*
   * d( mu )/dmesh
   */
  
  val = pow( gammadot+offset, nexp-2.);
  
  if ( d_mu != NULL ) d_mu->gd = mu0*(nexp - 1.0)*val;
/*   *d_mu_dgd -= tau_y/pow(gammadot+offset, 2.0); Disabling the sensitivities on this term
*  otherwise converges not 
*/
  
  if ( d_mu != NULL && pd->e[pg->imtrx][R_MESH1] )
    {
      for ( b=0; b<VIM; b++)
	{
	  for ( j=0; j<mdofs; j++)
	    {
	      if(gammadot != 0.0 && Include_Visc_Sens )
		{
		  
		  d_mu->X [b][j] =
		    d_mu->gd * d_gd_dmesh [b][j] ;
		  
		}
	      else
		{
		  /* printf("\ngammadot is zero in viscosity function");*/
		  d_mu->X [b][j] = 0.0;
		}
	    }
	}
    }
  
  /*
   * d( mu )/dv
   */
  /*
   * d( mu )/dv
   */
  if ( d_mu != NULL && pd->e[pg->imtrx][R_MOMENTUM1] )
    {
      for ( a=0; a<VIM; a++)
        {
          for ( i=0; i<vdofs; i++)
	    {
	      if(gammadot != 0.0 && Include_Visc_Sens )
	        {
	          d_mu->v[a][i] =
		    d_mu->gd * d_gd_dv[a][i] ;
	        }
	      else
	        {
	          d_mu->v[a][i] = 0.0 ;
	        }
	    }
        }
    }

  return(mu);
}


double
carreau_viscosity(struct Generalized_Newtonian *gn_local,
		  dbl gamma_dot[DIM][DIM], /* strain rate tensor */
		  VISCOSITY_DEPENDENCE_STRUCT *d_mu)
{

    int a, b;
  int mdofs=0,vdofs;


  int i, j;

  dbl gammadot;	                /* strain rate invariant */ 

  dbl d_gd_dv[DIM][MDE];        /* derivative of strain rate invariant 
				   wrt velocity */ 
  dbl d_gd_dmesh[DIM][MDE];     /* derivative of strain rate invariant 
				   wrt mesh */ 

  dbl val, val1, val2;
  dbl mu = 0.;
  dbl mu0;
  dbl muinf;
  dbl nexp;
  dbl aexp;
  dbl lambda;

  vdofs = ei[pg->imtrx]->dof[VELOCITY1];
  
  if ( pd->e[pg->imtrx][R_MESH1] )
    {
      mdofs = ei[pg->imtrx]->dof[R_MESH1];
    }
  
  calc_shearrate(&gammadot, gamma_dot, d_gd_dv, d_gd_dmesh);

  mu0 = gn_local->mu0;
  nexp = gn_local->nexp;
  muinf = gn_local->muinf;
  aexp = gn_local->aexp;
  lambda = gn_local->lam;

 if(gn_local->mu0Model == LEVEL_SET)
   {
     if(d_mu != NULL)
       {
	 level_set_property(gn_local->u_mu0[0], gn_local->u_mu0[1], 
			    gn_local->u_mu0[2], &mu0, d_mu->F);
       }
     else
       {
	 level_set_property(gn_local->u_mu0[0], gn_local->u_mu0[1], 
			    gn_local->u_mu0[2], &mu0, NULL);
       }
   }

 if(gn_local->nexpModel == LEVEL_SET)
   {
     if(d_mu != NULL)
       {
	 level_set_property(gn_local->u_nexp[0], gn_local->u_nexp[1], 
			    gn_local->u_nexp[2], &nexp, d_mu->F);
       }
     else
       {
	 level_set_property(gn_local->u_nexp[0], gn_local->u_nexp[1], 
			    gn_local->u_nexp[2], &nexp, NULL);
       }
   }

 if(gn_local->muinfModel == LEVEL_SET)
   {
     if(d_mu != NULL)
       {
	 level_set_property(gn_local->u_muinf[0], gn_local->u_muinf[1], 
			    gn_local->u_muinf[2], &muinf, d_mu->F);
       }
     else
       {
	 level_set_property(gn_local->u_muinf[0], gn_local->u_muinf[1], 
			    gn_local->u_muinf[2], &muinf, NULL);
       }
   }

 if(gn_local->aexpModel == LEVEL_SET)
   {
     if(d_mu != NULL)
       {
	 level_set_property(gn_local->u_aexp[0], gn_local->u_aexp[1], 
			    gn_local->u_aexp[2], &aexp, d_mu->F);
       }
     else
       {
	 level_set_property(gn_local->u_aexp[0], gn_local->u_aexp[1], 
			    gn_local->u_aexp[2], &aexp, NULL);
       }
   }
	  

 if(gn_local->lamModel == LEVEL_SET)
   {
     if(d_mu != NULL)
       {
	 level_set_property(gn_local->u_lam[0], gn_local->u_lam[1], 
			    gn_local->u_lam[2], &lambda, d_mu->F);
       }
     else
       {
	 level_set_property(gn_local->u_lam[0], gn_local->u_lam[1], 
			    gn_local->u_lam[2], &lambda, NULL);
       }
   }  

  
  if(gammadot != 0.)
    {
      val2 = pow( lambda*gammadot, aexp);
    }
  else
    {
      val2 = 0.;
    }
  val = pow(1. + val2,(nexp-1.)/aexp);
  mu = muinf + (mu0 - muinf)* val;
  
  /* gammadot = 0.0; */
  /* this effectively turns off the viscosity Jac terms */
  
  if(gammadot != 0.)
    {
      val = pow( lambda*gammadot, aexp-1.);
    }
  else
    {
      val = 0.;
    }
  val1 = pow(1. + val2,(nexp-1.-aexp)/aexp);
  
  if ( d_mu != NULL ) d_mu->gd = (mu0 - muinf)* (nexp-1.) * lambda * val * val1 ;
  
  /*
   * d( mu )/dmesh
   */
  if ( d_mu != NULL && pd->e[pg->imtrx][R_MESH1] )
    {
      for ( b=0; b<VIM; b++)
	{
	  for ( j=0; j<mdofs; j++)
	    {
	      if(gammadot != 0.0 && Include_Visc_Sens )
		{
		  d_mu->X [b][j] =
		    d_mu->gd * d_gd_dmesh [b][j] ;
		}
	      else
		{
		  /* printf("\ngammadot is zero in viscosity function");*/
		  d_mu->X [b][j] = 0.0;
		}
	    }
	}
    }
  
  /*
   * d( mu )/dv
   */
  if ( d_mu != NULL && pd->e[pg->imtrx][R_MOMENTUM1] )
    {
      for ( a=0; a<VIM; a++)
        {
          for ( i=0; i<vdofs; i++)
	    {
	      if(gammadot != 0.0 && Include_Visc_Sens )
	        {
	          d_mu->v[a][i] =
		    d_mu->gd * d_gd_dv[a][i] ;
	        }
	      else
	        {
	          d_mu->v[a][i] = 0.0 ;
	        }
	    }
        }
    }
  return(mu);
}

#define MELTING_BINGHAM FALSE

double
bingham_viscosity(struct Generalized_Newtonian *gn_local,
		  dbl gamma_dot[DIM][DIM], /* strain rate tensor */
		  VISCOSITY_DEPENDENCE_STRUCT *d_mu)
{

  int a, b;
  int var;
  int mdofs=0,vdofs;
  int i, j;

  dbl gammadot;	                /* strain rate invariant */ 

  dbl d_gd_dv[DIM][MDE];        /* derivative of strain rate invariant 
				   wrt velocity */ 
  dbl d_gd_dmesh[DIM][MDE];     /* derivative of strain rate invariant 
				   wrt mesh */ 

  dbl val1;
  dbl yield, shear, visc_cy;
  dbl d_yield = 0.0, d_shear, d_visc_cy, d_shear_d_at, d_yield_d_at = 0.0;
  dbl mu = 0.;
  dbl mu0;
  dbl muinf;
  dbl dmudT;
  dbl nexp;
  dbl atexp;
  dbl aexp;
  dbl at_shift;
  dbl lambda;
  dbl tau_y = 0.0;
  dbl fexp;
  dbl d_at_s;
  dbl temp;
#if MELTING_BINGHAM
  dbl tmelt;
#endif

  calc_shearrate(&gammadot, gamma_dot, d_gd_dv, d_gd_dmesh);

  mu0 = gn_local->mu0;
  nexp = gn_local->nexp;
  muinf = gn_local->muinf;
  aexp = gn_local->aexp;
  atexp = gn_local->atexp;
  lambda = gn_local->lam;
   if( gn_local->tau_yModel == CONSTANT )
	{
   	 tau_y = gn_local->tau_y;
 	}
   else if (gn_local->tau_yModel == USER )
 	{
 	  usr_yield_stress(gn_local->u_tau_y, tran->time_value);
 	  tau_y = gn_local->tau_y;
 	}
   else
 	{
 	EH(-1,"Invalid Yield Stress Model");
 	}
  fexp = gn_local->fexp;

  if ( pd->e[pg->imtrx][TEMPERATURE] )
       {temp = fv->T;}
  else
       {temp = upd->Process_Temperature;}

#if MELTING_BINGHAM
  tmelt = mp->melting_point_liquidus;
#endif

  vdofs = ei[pg->imtrx]->dof[VELOCITY1];
  
  if ( pd->e[pg->imtrx][R_MESH1] )
    {
      mdofs = ei[pg->imtrx]->dof[R_MESH1];
    }

  var = TEMPERATURE;
  if ( pd->e[pg->imtrx][var] && (temp != 0.) &&  (mp->reference[TEMPERATURE] != 0.))
    {
#if MELTING_BINGHAM
      /* melting version */
      if(temp <= tmelt)
        {
          at_shift = exp(-atexp*(mp->reference[TEMPERATURE]-temp)/
                 (tmelt-temp)/(mp->reference[TEMPERATURE]-tmelt));
	  if(!isfinite(at_shift)) { at_shift = DBL_MAX; }
                 
          d_at_s = - at_shift*atexp /(tmelt-temp)/(tmelt-temp);
        }
      else
        {
          at_shift = 1.;
          d_at_s = 0.;
        }
#else
      /* normal, non-melting version */
      at_shift = exp(atexp * (1./temp - 1./mp->reference[TEMPERATURE]));
      if(!isfinite(at_shift)) { at_shift = DBL_MAX; }
      d_at_s = -at_shift *atexp /(temp*temp) ;
#endif
    }
  else
    {
      at_shift = 1.;
      d_at_s = 0. ;
    }
  
  if((at_shift * lambda * gammadot) != 0.)
    {
      shear = pow( at_shift* lambda * gammadot, aexp);
      val1 = pow(at_shift * lambda * gammadot, aexp-1.);
      d_shear = aexp * at_shift * lambda * val1;
      d_shear_d_at = aexp * gammadot * lambda * val1; 
    }
  else
    {
      shear = 0.;
      d_shear = 0.;
      d_shear_d_at = 0.;
    }
  
  if((gammadot != 0.) && (at_shift != 0.))
    {
      yield = tau_y * (1. - exp(-at_shift*fexp*gammadot))/(at_shift*gammadot);
      d_yield = (-yield + tau_y*fexp * exp(-at_shift*fexp*gammadot) )/gammadot;
      d_yield_d_at = d_yield*gammadot/at_shift ;
    }
  else
    {
      yield = tau_y * fexp;
      d_yield = 0.;
    }
  
  visc_cy = pow(1. + shear,(nexp-1.)/aexp);
  d_visc_cy = (nexp-1.)/aexp * pow(1. + shear,(nexp-1.-aexp)/aexp);
  
  mu = at_shift * (muinf + (mu0 - muinf + yield) * visc_cy);
  
  
  if ( d_mu != NULL ) d_mu->gd =  at_shift * ( d_yield * visc_cy
			    + (mu0 - muinf + yield) * d_visc_cy *d_shear);
			   
#if MELTING_BINGHAM
  /* melting version */

  mu = at_shift * (muinf + (mu0 - muinf + yield*at_shift) * visc_cy);
  
  
  if ( d_mu != NULL ) d_mu->gd =  at_shift * ( d_yield*at_shift * visc_cy
		 + (mu0 - muinf + yield*at_shift) * d_visc_cy *d_shear);
			    
if(mu <= 1.)
    {
      mu = 1.;
      d_at_s = 0.;
      d_mu->gd = 0.;
      at_shift = 1.;
    }  			   
#endif		   
  
  /*
   * d( mu )/dT
   */
  
  var = TEMPERATURE;
  if ( d_mu != NULL && pd->e[pg->imtrx][var] )
    {
      if ( (temp != 0.) &&  (mp->reference[TEMPERATURE] != 0.) && (gammadot != 0.))
	{
#if MELTING_BINGHAM
          /* melting version */
	  dmudT = mu/at_shift -
	     at_shift*tau_y *fexp *exp(-at_shift*fexp*gammadot)* visc_cy  +
	     at_shift * (mu0 - muinf + yield*at_shift) * d_visc_cy * d_shear_d_at;
          dmudT *= d_at_s;
	  if(!isfinite(dmudT)) { dmudT = DBL_MAX; }
#else
 	  /* normal, non-melting version */
          dmudT = mu/at_shift +
	    at_shift* d_yield_d_at * visc_cy +
	    at_shift * (mu0 - muinf + yield) * d_visc_cy * d_shear_d_at;
 	  dmudT *= d_at_s;
	  if(!isfinite(dmudT)) { dmudT = DBL_MAX; }
#endif	
	}
      else
	{
	  dmudT = 0.;
	}
      
	   
      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	{
	  d_mu->T[j]= dmudT * bf[var]->phi[j];
	}
    }
  
  /*
   * d( mu )/dmesh
   */
  if ( d_mu != NULL && pd->e[pg->imtrx][R_MESH1] )
    {
      for ( b=0; b<VIM; b++)
	{
	  for ( j=0; j<mdofs; j++)
	    {
	      if(gammadot != 0.0 && Include_Visc_Sens )
		{
		  d_mu->X [b][j] =
		    d_mu->gd * d_gd_dmesh [b][j] ;
		}
	      else
		{
		  /* printf("\ngammadot is zero in viscosity function");*/
		  d_mu->X [b][j] = 0.0;
		}
	    }
	}
    }
  
  /*
   * d( mu )/dv
   */
  if ( d_mu != NULL && pd->e[pg->imtrx][R_MOMENTUM1] )
    {
      for ( a=0; a<VIM; a++)
        {
          for ( i=0; i<vdofs; i++)
	    {
	      if(gammadot != 0.0 && Include_Visc_Sens )
	        {
	          d_mu->v[a][i] =
		    d_mu->gd * d_gd_dv[a][i] ;
	        }
	      else
	        {
	          d_mu->v[a][i] = 0.0 ;
	        }
	    }
        }
    }

  return(mu);
}
double
bingham_wlf_viscosity(struct Generalized_Newtonian *gn_local,
		      dbl gamma_dot[DIM][DIM], /* strain rate tensor */
		      VISCOSITY_DEPENDENCE_STRUCT *d_mu)
{

  int a, b;
  int var;
  int mdofs=0,vdofs;
  int i, j;

  dbl gammadot;	                /* strain rate invariant */ 

  dbl d_gd_dv[DIM][MDE];        /* derivative of strain rate invariant 
				   wrt velocity */ 
  dbl d_gd_dmesh[DIM][MDE];     /* derivative of strain rate invariant 
				   wrt mesh */ 

  dbl val1;
  dbl yield, shear, visc_cy;
  dbl d_yield, d_shear, d_visc_cy, d_shear_d_at, d_yield_d_at;
  dbl mu = 0.;
  dbl mu0;
  dbl muinf;
  dbl dmudT;
  dbl nexp;
  dbl atexp;
  dbl aexp;
  dbl at_shift, d_at_dT = 0.0;
  dbl lambda;
  dbl tau_y = 0.0;
  dbl fexp;
  dbl Tref;
  dbl wlfc2, wlf_denom;
  dbl temp;

  calc_shearrate(&gammadot, gamma_dot, d_gd_dv, d_gd_dmesh);

  mu0 = gn_local->mu0;
  nexp = gn_local->nexp;
  muinf = gn_local->muinf;
  aexp = gn_local->aexp;
  atexp = gn_local->atexp;
  lambda = gn_local->lam;
  wlfc2 = gn_local->wlfc2;

  if ( pd->e[pg->imtrx][TEMPERATURE] )
       {temp = fv->T;}
  else
       {temp = upd->Process_Temperature;}

  if( gn_local->tau_yModel == CONSTANT )
    {
      tau_y = gn_local->tau_y;
    }
  else if (gn_local->tau_yModel == USER )
    {
      usr_yield_stress(gn_local->u_tau_y, tran->time_value);
      tau_y = gn_local->tau_y;
    }
  else
    {
      EH(-1,"Invalid Yield Stress Model");
    }
  fexp = gn_local->fexp;
  
  vdofs = ei[pg->imtrx]->dof[VELOCITY1];
  
  if ( pd->e[pg->imtrx][R_MESH1] )
    {
      mdofs = ei[pg->imtrx]->dof[R_MESH1];
    }

  var = TEMPERATURE;
  at_shift = 1.;
      Tref = mp->reference[TEMPERATURE];
      wlf_denom = wlfc2 + temp - Tref;
      if(wlf_denom != 0.)
	{
	  at_shift=exp(atexp*(Tref-temp)/wlf_denom);
          if(!isfinite(at_shift)) { at_shift = DBL_MAX; }
	  d_at_dT = at_shift*atexp/wlf_denom*(-1.-(Tref-temp)/wlf_denom);
	}
  
  if((at_shift * lambda * gammadot) != 0.)
    {
      shear = pow( at_shift* lambda * gammadot, aexp);
      val1 = pow(at_shift * lambda * gammadot, aexp-1.);
      d_shear = aexp * at_shift * lambda * val1;
      d_shear_d_at = aexp * gammadot * lambda * val1; 
    }
  else
    {
      shear = 0.;
      d_shear = 0.;
      d_shear_d_at = 0.;
    }
  
  if((gammadot != 0.) && (at_shift != 0.))
    {
      yield = tau_y * (1. - exp(-at_shift*fexp*gammadot))/(at_shift*gammadot);
      d_yield = (-yield + tau_y*fexp * exp(-at_shift*fexp*gammadot) )/gammadot;
      d_yield_d_at = d_yield*gammadot/at_shift ;
   }
  else
    {
      yield = tau_y * fexp;
      d_yield = 0.;
      d_yield_d_at = 0.;
    }
  
  visc_cy = pow(1. + shear,(nexp-1.)/aexp);
  d_visc_cy = (nexp-1.)/aexp * pow(1. + shear,(nexp-1.-aexp)/aexp);
  
  mu = at_shift * (muinf + (mu0 - muinf + yield) * visc_cy);
  
  
  if ( d_mu != NULL ) d_mu->gd =  at_shift * ( d_yield * visc_cy
			    + (mu0 - muinf + yield) * d_visc_cy *d_shear);
  
  /*
   * d( mu )/dT
   */
  
  var = TEMPERATURE;
  if ( d_mu != NULL && pd->e[pg->imtrx][var] )
    {
      if ( (temp != 0.) &&  (mp->reference[TEMPERATURE] != 0.) && (gammadot != 0.))
	{
          dmudT = mu/at_shift +
	    at_shift* d_yield_d_at * visc_cy +
	    at_shift * (mu0 - muinf + yield) * d_visc_cy * d_shear_d_at;
 	  dmudT *= d_at_dT;
	  if(!isfinite(dmudT)) { dmudT = DBL_MAX; }
	}
      else
	{
	  dmudT = 0.;
	}
      
	   
      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	{
	  d_mu->T[j]= dmudT * bf[var]->phi[j];
	}
    }
  
  /*
   * d( mu )/dmesh
   */
  if ( d_mu != NULL && pd->e[pg->imtrx][R_MESH1] )
    {
      for ( b=0; b<VIM; b++)
	{
	  for ( j=0; j<mdofs; j++)
	    {
	      if(gammadot != 0.0 && Include_Visc_Sens )
		{
		  d_mu->X [b][j] =
		    d_mu->gd * d_gd_dmesh [b][j] ;
		}
	      else
		{
		  /* printf("\ngammadot is zero in viscosity function");*/
		  d_mu->X [b][j] = 0.0;
		}
	    }
	}
    }
  
  /*
   * d( mu )/dv
   */
  if ( d_mu != NULL && pd->e[pg->imtrx][R_MOMENTUM1] )
    {
      for ( a=0; a<VIM; a++)
        {
          for ( i=0; i<vdofs; i++)
	    {
	      if(gammadot != 0.0 && Include_Visc_Sens )
	        {
	          d_mu->v[a][i] =
		    d_mu->gd * d_gd_dv[a][i] ;
	        }
	      else
	        {
	          d_mu->v[a][i] = 0.0 ;
	        }
	    }
        }
    }

  return(mu);
}





double
carreau_wlf_viscosity(struct Generalized_Newtonian *gn_local,
		  dbl gamma_dot[DIM][DIM], /* strain rate tensor */
		  VISCOSITY_DEPENDENCE_STRUCT *d_mu)
{

  int a, b;

  int var;
  int mdofs=0,vdofs;

  int i, j;

  dbl gammadot;	                /* strain rate invariant */ 

  dbl d_gd_dv[DIM][MDE];        /* derivative of strain rate invariant 
				   wrt velocity */ 
  dbl d_gd_dmesh[DIM][MDE];     /* derivative of strain rate invariant 
				   wrt mesh */ 

  dbl val1;
  dbl shear, visc_cy;
  dbl d_shear, d_visc_cy, d_shear_d_at;
  dbl mu = 0.;
  dbl mu0;
  dbl muinf;
  dbl dmudT;
  dbl nexp;
  dbl atexp;
  dbl aexp;
  dbl at_shift;
  dbl lambda;
  dbl wlf_denom = 0.0;
  dbl wlfc2;
  dbl temp;

  calc_shearrate(&gammadot, gamma_dot, d_gd_dv, d_gd_dmesh);

  mu0 = gn_local->mu0;
  nexp = gn_local->nexp;
  muinf = gn_local->muinf;
  aexp = gn_local->aexp;
  atexp = gn_local->atexp;
  wlfc2 = gn_local->wlfc2;
  lambda = gn_local->lam;

  if ( pd->e[pg->imtrx][TEMPERATURE] )
       {temp = fv->T;}
  else
       {temp = upd->Process_Temperature;}

 if(gn_local->mu0Model == LEVEL_SET)
   {
     if(d_mu != NULL)
       {
	 level_set_property(gn_local->u_mu0[0], gn_local->u_mu0[1], 
			    gn_local->u_mu0[2], &mu0, d_mu->F);
       }
     else
       {
	 level_set_property(gn_local->u_mu0[0], gn_local->u_mu0[1], 
			    gn_local->u_mu0[2], &mu0, NULL);
       }
   }

 if(gn_local->nexpModel == LEVEL_SET)
   {
     if(d_mu != NULL)
       {
	 level_set_property(gn_local->u_nexp[0], gn_local->u_nexp[1], 
			    gn_local->u_nexp[2], &nexp, d_mu->F);
       }
     else
       {
	 level_set_property(gn_local->u_nexp[0], gn_local->u_nexp[1], 
			    gn_local->u_nexp[2], &nexp, NULL);
       }
   }

 if(gn_local->muinfModel == LEVEL_SET)
   {
     if(d_mu != NULL)
       {
	 level_set_property(gn_local->u_muinf[0], gn_local->u_muinf[1], 
			    gn_local->u_muinf[2], &muinf, d_mu->F);
       }
     else
       {
	 level_set_property(gn_local->u_muinf[0], gn_local->u_muinf[1], 
			    gn_local->u_muinf[2], &muinf, NULL);
       }
   }

 if(gn_local->aexpModel == LEVEL_SET)
   {
     if(d_mu != NULL)
       {
	 level_set_property(gn_local->u_aexp[0], gn_local->u_aexp[1], 
			    gn_local->u_aexp[2], &aexp, d_mu->F);
       }
     else
       {
	 level_set_property(gn_local->u_aexp[0], gn_local->u_aexp[1], 
			    gn_local->u_aexp[2], &aexp, NULL);
       }
   }
	  
 if(gn_local->atexpModel == LEVEL_SET)
   {
     if(d_mu != NULL)
       {
	 level_set_property(gn_local->u_atexp[0], gn_local->u_atexp[1], 
			    gn_local->u_atexp[2], &atexp, d_mu->F);
       }
     else
       {
	 level_set_property(gn_local->u_atexp[0], gn_local->u_atexp[1], 
			    gn_local->u_atexp[2], &atexp, NULL);
       }
   }

 if(gn_local->wlfc2Model == LEVEL_SET)
   {
     if(d_mu != NULL)
       {
	 level_set_property(gn_local->u_wlfc2[0], gn_local->u_wlfc2[1], 
			    gn_local->u_wlfc2[2], &wlfc2, d_mu->F);
       }
     else
       {
	 level_set_property(gn_local->u_wlfc2[0], gn_local->u_wlfc2[1], 
			    gn_local->u_wlfc2[2], &wlfc2, NULL);
       }
   }

 if(gn_local->lamModel == LEVEL_SET)
   {
     if(d_mu != NULL)
       {
	 level_set_property(gn_local->u_lam[0], gn_local->u_lam[1], 
			    gn_local->u_lam[2], &lambda, d_mu->F);
       }
     else
       {
	 level_set_property(gn_local->u_lam[0], gn_local->u_lam[1], 
			    gn_local->u_lam[2], &lambda, NULL);
       }
   }  

  vdofs = ei[pg->imtrx]->dof[VELOCITY1];
  
  if ( pd->e[pg->imtrx][R_MESH1] )
    {
      mdofs = ei[pg->imtrx]->dof[R_MESH1];
    }

  at_shift = 1.;
      wlf_denom = wlfc2 + temp - mp->reference[TEMPERATURE];
      if(wlf_denom != 0.)
	{
	  at_shift=exp(atexp*(mp->reference[TEMPERATURE]-temp)/wlf_denom);
	  if(!isfinite(at_shift)) { at_shift = DBL_MAX; }
	}
  
  if(gammadot != 0.)
    {
      shear = pow( at_shift* lambda * gammadot, aexp);
      val1 = pow(at_shift * lambda * gammadot, aexp-1.);
      d_shear = aexp * at_shift * lambda * val1;
      d_shear_d_at = aexp * gammadot * lambda * val1; 
    }
  else
    {
      shear = 0.;
      d_shear = 0.;
      d_shear_d_at = 0.;
    }
  
  
  visc_cy = pow(1. + shear,(nexp-1.)/aexp);
  d_visc_cy = (nexp-1.)/aexp * pow(1. + shear,(nexp-1.-aexp)/aexp);
  
  mu = at_shift * (muinf + (mu0 - muinf) * visc_cy);
  
  if ( d_mu != NULL ) d_mu->gd =  at_shift * ((mu0 - muinf) * d_visc_cy *d_shear);
  
  /*
   * d( mu )/dT
   */
  
  var = TEMPERATURE;
  if ( d_mu != NULL && pd->e[pg->imtrx][var] )
    {
      if(wlf_denom != 0.)
	{
	  dmudT = mu +
	  at_shift * at_shift * (mu0 - muinf) * d_visc_cy * d_shear_d_at;
	  dmudT *= - atexp*wlfc2/(wlf_denom*wlf_denom);
	  if(!isfinite(dmudT)) { dmudT = DBL_MAX; }
	}
      else
	{
	  dmudT = 0.;
	}
      
      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	{
	  d_mu->T[j]= dmudT * bf[var]->phi[j];
	}
    }
  
  /*
   * d( mu )/dmesh
   */
  if ( d_mu != NULL && pd->e[pg->imtrx][R_MESH1] )
    {
      for ( b=0; b<VIM; b++)
	{
	  for ( j=0; j<mdofs; j++)
	    {
	      if(gammadot != 0.0 && Include_Visc_Sens )
		{
		  d_mu->X [b][j] =
		    d_mu->gd * d_gd_dmesh [b][j] ;
		}
	      else
		{
		  /* printf("\ngammadot is zero in viscosity function");*/
		  d_mu->X [b][j] = 0.0;
		}
	    }
	}
    }
  
  /*
   * d( mu )/dv
   */
  if ( d_mu != NULL && pd->e[pg->imtrx][R_MOMENTUM1] )
    {
      for ( a=0; a<VIM; a++)
        {
          for ( i=0; i<vdofs; i++)
	    {
	      if(gammadot != 0.0 && Include_Visc_Sens )
	        {
	          d_mu->v[a][i] =
		    d_mu->gd * d_gd_dv[a][i] ;
	        }
	      else
	        {
	          d_mu->v[a][i] = 0.0 ;
	        }
	    }
        }
    }

  return(mu);
}









/* 
 * Viscosity model for FILL equation with a
 * fluid of one viscosity displacing a fluid
 * of a different viscosity
 */


/*
 * int fill_viscosity ()
 *
 * ----------------------------------------------------------------------------
 * This routine is responsible for filling up the following portions of the mp structures
 * at the current gauss point:
 *     intput:  param - array of constants input on the property card.  
 *
 *     output:  mu       => mp->viscosity - viscosity
 */

int 
fill_viscosity(dbl *param)	/* ptr to the user-defined parameter list    */
{
  dbl mu; /* Newtonian viscosity*/
  dbl F; /* Convenient local variables */


 /***********Load up convenient local variables*************/
  F = fv->F;

 /**********************************************************/
  if ( ls != NULL )
    {
      EH(-1, "Use LEVEL_SET instead of FILL viscosity model.");
      return(0);
    }

   mu = (1.- F)*param[1] + F*param[0];
   if(mu <= 0. )
     {
       mu = param[1];
     }
   else if(mu > param[0] )
     {
       mu = param[0];
     }

   mp->viscosity = mu;
   return(0);
}

/*
 *
 *  Suspension viscosity model
 *
 */

/******************************************************************************
 *     Function that computes the viscosity of a solid suspension
 *     via the relation
 *                     mu = mu0*(1. - C/Cmaxpack)^nexp
 *       where
 *                      mu0 = solvent viscosity
 *                      Cmaxpack  = binding solid fraction, i.e. solid fraction for infinite viscosity
 *                      nexp   = exponent
 *
 *     Function sets the viscosity members of mp.
 *
 *    Author: Tom Baer
 *      Date: 11/26/96
 *   Revised: 3/24/97 RRR
 *            1/21/2010 KT -> Add capability for dealing with shell particle concentration
 *
 *
 *****************************************************************************/
int
suspension_viscosity(int species, /* species for solid volume fraction track */
		     dbl mu0,	/* carrier fluid viscosity                   */
		     dbl maxpack, /* maximum solid volume fraction           */
		     dbl nexp, 	/* exponent for constitutive equation        */
		     dbl C )   /* Concentration  */
{
  int a;
  int w;

  dbl mu;   /* viscosity */

  /* Concentration at previous time step */
  dbl C_old = fv_old->sh_pc;

  /* Viscosity at previous time step */
  dbl mu_old;
  
  int status = 1;

  /* initialize everything */

  mp->d_viscosity[TEMPERATURE]	= 0.0;
  mp->d2_viscosity[TEMPERATURE]	= 0.0;
  mp->d_viscosity[PRESSURE]	= 0.0;
  mp->d2_viscosity[PRESSURE]	= 0.0;
  mp->d_viscosity[SHEAR_RATE]	= 0.0;
  mp->d2_viscosity[SHEAR_RATE]	= 0.0;
  mp->d_viscosity[SHELL_PARTC] = 0.0;


  for ( a=0; a<DIM; a++)
    {
      mp->d_viscosity[VELOCITY1+a] = 0.0;
      mp->d2_viscosity[VELOCITY1+a] = 0.0;
      mp->d_viscosity[MESH_DISPLACEMENT1+a] = 0.0;
      mp->d2_viscosity[MESH_DISPLACEMENT1+a] = 0.0;
    }
  for(w=0; w<pd->Num_Species_Eqn; w++)
    {
      mp->d_viscosity[MAX_VARIABLE_TYPES + w] = 0.0;
      mp->d2_viscosity[MAX_VARIABLE_TYPES + w] = 0.0;
    }
  /* Implement a cutoff concentration that caps the viscosity at a lower
     value. This should make the viscosity more stable and allow the
     solidification Brinkman porous source term to slow the flow appropriately.
     Here I have chosen a cutoff value 10% below maximum packing. However,
     this may be too low and I will have to adjust it.
     */

  if ( nexp > 0.0 || ( C > 0.0 && C < 0.98 * maxpack)  )
    {    

      mu = mu0*pow( 1.0 - C/maxpack, nexp );
      mu_old = mu0*pow( 1.0 - C_old/maxpack, nexp );

      mp->viscosity = mu;
      mp_old->viscosity = mu_old;

      if (pd->v[pg->imtrx][SHELL_PARTC]  )
         {

           /* dmu/dc */
           mp->d_viscosity[SHELL_PARTC] = -mu*nexp/( maxpack - C) ;

           /* d2mu/dc2 .... */
           mp->d2_viscosity[SHELL_PARTC] = mu0*nexp*(nexp-1.0)*pow( 1.0 
                                           - C/maxpack,nexp-2. )/maxpack/maxpack;
         }

      else
         {

           /* dmu/dc */
           mp->d_viscosity[MAX_VARIABLE_TYPES+species] = -mu*nexp/( maxpack - C) ;

          /* d2mu/dc2 .... */
           mp->d2_viscosity[MAX_VARIABLE_TYPES+species]= mu0*nexp*(nexp-1.0)*pow( 1.0 
                                                         - C/maxpack,nexp-2. )/maxpack/maxpack;
         }
 
    }
  else if ( C <= 0. )
    {
      mu = mu0;
      mu_old = mu;
      mp->viscosity = mu;
      mp_old->viscosity = mu_old;
    }

  else if ( C >=  0.98 *maxpack )
    {
      mu = mu0*pow( 0.02, nexp );
      mu_old = mu;
      mp->viscosity = mu;
      mp_old->viscosity = mu_old;

      /* setting the derivatives to zero for this range seems to improve convergence. Odd.*/
    }

  return(status);
} /* end of suspension_viscosity */

/*
 *
 *  Carreau_suspension viscosity model
 *
 */

double
carreau_suspension_viscosity(struct Generalized_Newtonian *gn_local,
			     dbl gamma_dot[DIM][DIM], /* strain rate tensor  */
			     VISCOSITY_DEPENDENCE_STRUCT *d_mu)
{
  int a, b;
  int var, var_offset;
  int w;

  int mdofs=0,vdofs;

  dbl C[MAX_CONC]; /* Convenient local variables */
  dbl mu = 0.;
  dbl mu_car, mu_sus=0;

  dbl val, val1, val2;
  dbl muinf;
  dbl aexp;
  dbl lambda;
  dbl nexp;       

  dbl gammadot;	                /* strain rate invariant */ 

  dbl d_gd_dv[DIM][MDE];        /* derivative of strain rate invariant 
				   wrt velocity */ 
  dbl d_gd_dmesh[DIM][MDE];     /* derivative of strain rate invariant 
				   wrt mesh */ 

  /* parameters for suspension model */
  dbl maxpack;    /* maximum solid volume fraction */
  dbl nexp_species;  /* exponent for constitutive equation */
  dbl mu0;        /* zero shear-rate viscosity */

  int species;    /* species number for solid volume fraction tracking */

  int i, j;
  
  calc_shearrate(&gammadot, gamma_dot, d_gd_dv, d_gd_dmesh);

  vdofs = ei[pg->imtrx]->dof[VELOCITY1];
  
  if ( pd->e[pg->imtrx][R_MESH1] )
    {
      mdofs = ei[pg->imtrx]->dof[R_MESH1];
    }

  /* initialize everything */
                                       
  mp->d_viscosity[TEMPERATURE] = 0.0; 
  mp->d2_viscosity[TEMPERATURE] = 0.0;
  mp->d_viscosity[PRESSURE]    = 0.0;	
  mp->d2_viscosity[PRESSURE]    = 0.0;

  for ( a=0; a<DIM; a++)                                
    {							
      mp->d_viscosity[VELOCITY1+a] = 0.0;      
      mp->d2_viscosity[VELOCITY1+a] = 0.0;
      mp->d_viscosity[MESH_DISPLACEMENT1+a] = 0.0;    
      mp->d2_viscosity[MESH_DISPLACEMENT1+a] = 0.0; 		
    }                                                    
  for(w=0; w<pd->Num_Species_Eqn; w++)
    {					
      mp->d_viscosity[MAX_VARIABLE_TYPES + w] = 0.0;
      mp->d2_viscosity[MAX_VARIABLE_TYPES + w] = 0.0;
      C[w] = fv->c[w];
    }

  nexp = gn_local->nexp;
  muinf = gn_local->muinf;
  mu0 = gn_local->mu0;
  aexp = gn_local->aexp;
  lambda = gn_local->lam;
  maxpack = gn_local->maxpack;   
  nexp_species = gn_local->atexp;      
  species = gn_local->sus_species_no;    
  
  if(gammadot != 0.)
    {
      val2 = pow( lambda*gammadot, aexp);
    }
  else
    {
      val2 = 0.;
    }
  val = pow(1. + val2,(nexp-1.)/aexp);
  mu_car = muinf + (mu0 - muinf)* val;


  if ( nexp_species > 0.0 || (C[species] > 0.0 
			      && C[species] < (maxpack - .01))  )
    {
      mu_sus = pow( 1.0 - C[species]/maxpack, nexp_species );
      /* dmu/dc */
      mp->d_viscosity[MAX_VARIABLE_TYPES+species] 
	= -mu_car * mu_sus * nexp_species/( maxpack - C[species]);

      /* d2mu/dc2 .... */

      mp->d2_viscosity[MAX_VARIABLE_TYPES+species] 
	= mu_car*nexp_species*(nexp_species-1.0)
	*pow( 1.0 - C[species]/maxpack,nexp_species-2. )/maxpack/maxpack; 
    }
  else if ( C[species] <= 0. )
    {    
      mu_sus = 1.;
    }

  else if ( C[species] >=  (maxpack - .01) ) 
    {      
      mu_sus = pow( 0.01/maxpack, nexp_species);
    }

  mu = mu_car * mu_sus;
  mp->viscosity = mu;
							 

  if(gammadot != 0.)
    {
      val = pow( lambda*gammadot, aexp-1.);
    }
  else
    {
      val = 0.;
    }
  val1 = pow(1. + val2,(nexp-1.-aexp)/aexp);
  val1 = pow(1. + val2,(nexp-1.-aexp)/aexp);
  
  if ( d_mu != NULL )
    {
      d_mu->gd = (mu0 - muinf)* (nexp-1.) * lambda * val * val1 * mu_sus;

      if(C[species] < maxpack)
        {
          mp->d2_viscosity[MAX_VARIABLE_TYPES + MAX_CONC] =
            -d_mu->gd * nexp_species/( maxpack - C[species]);
        }

      mp->d2_viscosity[SHEAR_RATE]  = d_mu->gd *
        ((aexp-1)/(lambda*gammadot) + lambda*(nexp-1.-aexp)/(1+val2));
    }
  
  /*
   * d( mu )/dmesh
   */
  if ( d_mu != NULL && pd->e[pg->imtrx][R_MESH1] )
    {
      for ( b=0; b<VIM; b++)
	{
	  for ( j=0; j<mdofs; j++)
	    {
	      if(gammadot != 0.0 && Include_Visc_Sens )
		{
		  d_mu->X [b][j] =
		    d_mu->gd * d_gd_dmesh [b][j] ;
		}
	      else
		{
		  /* printf("\ngammadot is zero in viscosity function");*/
		  d_mu->X [b][j] = 0.0;
		}
	    }
	}
    }
  
  /*
   * d( mu )/dv
   */
  if ( d_mu != NULL && pd->e[pg->imtrx][R_MOMENTUM1] )
    {
      for ( a=0; a<VIM; a++)
        {
          for ( i=0; i<vdofs; i++)
	    {
	      if(gammadot != 0.0 && Include_Visc_Sens )
	        {
	          d_mu->v[a][i] =
		    d_mu->gd * d_gd_dv[a][i] ;
	        }
	      else
	        {
	          d_mu->v[a][i] = 0.0 ;
	        }
	    }
        }
    }
  
  /*
   * d( mu )/dc
   */
  var = MASS_FRACTION;
  if ( d_mu != NULL && pd->v[pg->imtrx][var] )
    {
      var_offset = MAX_VARIABLE_TYPES + species;
      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	{
	  d_mu->C[species][j] =mp->d_viscosity[var_offset]*bf[var]->phi[j];
	}
    }
  

  return(mu);
} /* end of carreau_suspension_viscosity */


/* powerlaw_suspension viscosity model */

double
powerlaw_suspension_viscosity(struct Generalized_Newtonian *gn_local,
			      dbl gamma_dot[DIM][DIM], /* strain rate tensor */
			      VISCOSITY_DEPENDENCE_STRUCT *d_mu)
{
  int a, b;
  int var, var_offset;
  int w;
  int mdofs=0,vdofs;

  dbl C[MAX_CONC];		/* Convenient local variables */
  dbl mu = 0.;
  dbl mu_car, mu_sus=0;

  dbl val;
  dbl offset;
  dbl nexp;       

  dbl gammadot;	                /* strain rate invariant */ 

  dbl d_gd_dv[DIM][MDE];        /* derivative of strain rate invariant 
				   wrt velocity */ 
  dbl d_gd_dmesh[DIM][MDE];     /* derivative of strain rate invariant 
				   wrt mesh */ 

  /* parameters for suspension model */
  dbl maxpack;    /* maximum solid volume fraction */
  dbl nexp_species;  /* exponent for constitutive equation */
  dbl mu0;        /* zero shear-rate viscosity */

  int species;    /* species number for solid volume fraction tracking */

  int i, j;

  calc_shearrate(&gammadot, gamma_dot, d_gd_dv, d_gd_dmesh);

  vdofs = ei[pg->imtrx]->dof[VELOCITY1];
  
  if ( pd->e[pg->imtrx][R_MESH1] )
    {
      mdofs = ei[pg->imtrx]->dof[R_MESH1];
    }

  /* initialize everything */
                                       
  mp->d_viscosity[TEMPERATURE] = 0.0;	
  mp->d2_viscosity[TEMPERATURE] = 0.0;
  mp->d_viscosity[PRESSURE]    = 0.0;	
  mp->d2_viscosity[PRESSURE]    = 0.0;

  for ( a=0; a<DIM; a++)                                
    {							
      mp->d_viscosity[VELOCITY1+a] = 0.0;      
      mp->d2_viscosity[VELOCITY1+a] = 0.0;
      mp->d_viscosity[MESH_DISPLACEMENT1+a] = 0.0;    
      mp->d2_viscosity[MESH_DISPLACEMENT1+a] = 0.0; 		
    }                                                    
  for(w=0; w<pd->Num_Species_Eqn; w++)
    {					
      mp->d_viscosity[MAX_VARIABLE_TYPES + w] = 0.0;
      mp->d2_viscosity[MAX_VARIABLE_TYPES + w] = 0.0;
      C[w] = fv->c[w];
    }

  nexp = gn_local->nexp;
  mu0 = gn_local->mu0;
  maxpack = gn_local->maxpack;   
  nexp_species = gn_local->atexp;      
  species = gn_local->sus_species_no;    
  offset = 0.00001;
  

  val = pow( gammadot + offset, nexp-1.);
  mu_car = mu0 * val;


  if ( nexp_species > 0.0 || (C[species] > 0.0 
			      && C[species] < (maxpack - .01))  )
    {
      mu_sus = pow( 1.0 - C[species]/maxpack, nexp_species );
      /* dmu/dc */
      mp->d_viscosity[MAX_VARIABLE_TYPES+species] 
	= -mu_car * mu_sus * nexp_species/( maxpack - C[species]);

      /* d2mu/dc2 .... */

      mp->d2_viscosity[MAX_VARIABLE_TYPES+species] 
	= mu_car*nexp_species*(nexp_species-1.0)
	*pow( 1.0 - C[species]/maxpack,nexp_species-2. )/maxpack/maxpack; 
    }
  else if ( C[species] <= 0. )
    {    
      mu_sus = 1.;
    }

  else if ( C[species] >=  (maxpack - .01) ) 
    {      
      mu_sus = pow( 0.01/maxpack, nexp_species);
    }

  mu = mu_car * mu_sus;
  mp->viscosity = mu;

  if ( d_mu != NULL ) d_mu->gd = mu *(nexp - 1.0)/(gammadot + offset);

  mp->d2_viscosity[SHEAR_RATE]  = mu *(nexp - 1.0)
    *(nexp - 2.0)/SQUARE(gammadot + offset);
  
  /*
   * d( mu )/dmesh
   */
  if ( d_mu != NULL && pd->e[pg->imtrx][R_MESH1] )
    {
      for ( b=0; b<VIM; b++)
	{
	  for ( j=0; j<mdofs; j++)
	    {
	      if(gammadot != 0.0 && Include_Visc_Sens )
		{
		  d_mu->X [b][j] =
		    d_mu->gd * d_gd_dmesh [b][j] ;
		}
	      else
		{
		  /* printf("\ngammadot is zero in viscosity function");*/
		  d_mu->X [b][j] = 0.0;
		}
	    }
	}
    }
  
  /*
   * d( mu )/dv
   */
  if ( d_mu != NULL && pd->e[pg->imtrx][R_MOMENTUM1] )
    {
      for ( a=0; a<VIM; a++)
        {
          for ( i=0; i<vdofs; i++)
	    {
	      if(gammadot != 0.0 && Include_Visc_Sens )
	        {
	          d_mu->v[a][i] =
		    d_mu->gd * d_gd_dv[a][i] ;
	        }
	      else
	        {
	          d_mu->v[a][i] = 0.0 ;
	        }
	    }
        }
    }
  
  /*
   * d( mu )/dc
   */
  var = MASS_FRACTION;
  if ( d_mu != NULL && pd->v[pg->imtrx][var] )
    {
      var_offset = MAX_VARIABLE_TYPES + species;
      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	{
	  d_mu->C[species][j] =mp->d_viscosity[var_offset]*bf[var]->phi[j];
	}
    }
  

  return(mu);
} /* end of powerlaw_suspension_viscosity */



/*
 *
 *  epoxy  viscosity model
 *
 */

/******************************************************************************
 *     Function that computes the viscosity of an epoxy that is polymerizing
 *     and whose behavior depends on the extent of reaction and the
 *     temperature of the system via the relation:
 *
 *                     mu = mu0* alpha_g/(alpha-alpha_g)^(A+ B*alpha)
 *                          * exp(Aexp/T)
 *       where
 *                      mu0       = solvent viscosity
 *                      alpha     = extent of reaction
 *                      alpha_g   = extent of reaction at the gel point
 *                      A, B      = exponent for cure behavior
 *                      Aexp      = exponent for temperature dependence of viscosity
 *
 *     Function sets the viscosity members of mp.
 *
 *    Author: RRR
 *      Date: 3/24/97
 *
 *
 *****************************************************************************/



int
epoxy_viscosity(int species,    /* species number for cure equation */
		dbl mu0,        /* monomer reference temperature viscosity */
		dbl alpha_g,    /* extent of reaction at the gel point */
		dbl A,          /* exponent for constitutive equation */
		dbl B,          /* exponent for constitutive equation */
		dbl Aexp)       /* exponent for thermal viscosity dependence */
{
  /* Local Variables */
  dbl mu;   /* viscosity */
  dbl alpha;  /* extent of reaction */
  dbl exponent;
  dbl ratio;
  dbl deriv;  /* stuff for the first derivative */
  dbl d_deriv; 
  dbl T;      /* Convenient local variables */
  int var;
  int status = 1;
  
  alpha = fv->c[species]; /* extent of reaction */
  
  if(alpha < alpha_g)
    {
      ratio = (alpha_g)/(alpha_g - alpha);
      exponent = A + B * alpha;
    }
  else /* do something special at the gel point */
    {
      ratio = 100000;
      exponent = A + B*alpha_g;
    }
    
  
  if ( pd->e[pg->imtrx][TEMPERATURE] )
       {T = fv->T;}
  else
       {T = upd->Process_Temperature;}

  if(T <= 0.)
    {
      mu = mu0  * pow ( ratio, exponent );
    }
  else
    {
      mu = mu0 * exp (Aexp/T) * pow ( ratio, exponent );
    }
    
  mp->viscosity = mu;
  
  /* dmu_dT */
  var = TEMPERATURE;
  if(pd->v[pg->imtrx][var])
    {
      if( T <= 0.)
	{
	  mp->d_viscosity[var] 	= 0.; 
	  mp->d2_viscosity[var] = 0.; 
	}
      else
	{
	  mp->d_viscosity[var] = -mu * Aexp/(T*T) ;
	  mp->d2_viscosity[var] 	= -mu* SQUARE( Aexp/(T*T))
	    +2.* mu * Aexp/(T*T*T);
	}
    }
  
  /* dmu/dc */
  var  = MASS_FRACTION;
  if(pd->v[pg->imtrx][var])
    {
      if(alpha < alpha_g)
	{
	  deriv = exponent/(alpha_g - alpha) + B*log(ratio);
	  d_deriv = 2.*B/(alpha_g - alpha) + exponent/SQUARE(alpha_g - alpha);
	  mp->d_viscosity[MAX_VARIABLE_TYPES+species] = mu * deriv;
	  mp->d2_viscosity[MAX_VARIABLE_TYPES+species] =  mu * deriv * deriv
	    + mu * d_deriv;
	}
      else
	{
	  mp->d_viscosity[MAX_VARIABLE_TYPES+species] = 0.;
	  mp->d2_viscosity[MAX_VARIABLE_TYPES+species] = 0.;
	}
	
    }

  return(status);
} /* end of epoxy_viscosity */


/******************************************************************************
 *     Function that computes the viscosity of foam pmdi10 model
 *
 *                     mu = mu0 * exp(E/RT) * ((alpha_g^A - alpha^A)/alpha^A)^(-B)
 *       where
 *                      mu0       = solvent viscosity
 *                      alpha     = extent of reaction
 *                      alpha_g   = extent of reaction at the gel point
 *                      A, B      = exponent for cure behavior
 *                      E/R       = Normalized activation energy
 *
 *     Function sets the viscosity members of mp.
 *
 *
 *****************************************************************************/



int
foam_pmdi10_viscosity(int species,    /* species number for cure equation */
		      dbl mu0,        /* monomer reference temperature viscosity */
		      dbl alpha_g,    /* extent of reaction at the gel point */
		      dbl A,          /* exponent for constitutive equation */
		      dbl B,          /* exponent for constitutive equation */
		      dbl norm_E)     /* Normalized activation energy */
{
  /* Local Variables */
  dbl mu;   /* viscosity */
  double muL;
  dbl alpha;  /* extent of reaction */
  dbl ratio;
  dbl deriv;  /* stuff for the first derivative */
  dbl T;      /* Convenient local variables */
  int var;
  int status = 1;
  double volF;
  int w;

  alpha = fv->c[species]; /* extent of reaction */

  double alpha_g_pow = pow(alpha_g, A);

  if(alpha < alpha_g)
    {
      if (fabs(alpha) > 1e-8)
	{
	  ratio = (alpha_g_pow - pow(alpha, A))/alpha_g_pow;
	}
      else
	{
	  ratio = 1.0;
	}
    }
  else /* do something special at the gel point */
    {
      ratio = 1e8;
    }

  T = upd->Process_Temperature;
  if ( pd->gv[TEMPERATURE] )
    {
      T = fv->T;
    }


  if(T <= 0.)
    {
      muL = mu0  * pow (ratio, -B );
    }
  else
    {
      muL = mu0 * exp(-norm_E/T) * pow(ratio, -B);
    }

  volF = mp->volumeFractionGas;

  if (volF > 0.98) {
    volF = 0.98;
  }

  mu = muL * exp(volF / (1- volF));
  mp->viscosity = mu;
  mp->FlowingLiquid_viscosity = muL;

  double partial = (1 / ((1-volF)*(1-volF))) * exp(volF / (1 - volF));

  /* dmu_dT */
  var = TEMPERATURE;
  if (T <= 0) {
    mp->d_viscosity[var] = 0;
    mp->d_FlowingLiquid_viscosity[var] = 0;
  } else {
    mp->d_FlowingLiquid_viscosity[var] = muL * norm_E/(T*T);
    mp->d_viscosity[var] = mu * norm_E/(T*T) +
      muL * mp->d_volumeFractionGas[var] * partial;
  }
  mp->d2_viscosity[var] = 0.0; // This doesn't seem to be used ever

  if(pd->gv[var])
    {
      if( T <= 0.)
	{
	  mp->d_viscosity[var] 	= 0.;
	  mp->d2_viscosity[var] = 0.;
	}
    }


  /* dmu/dc */
  var  = MASS_FRACTION;
  mp->d_viscosity[MAX_VARIABLE_TYPES+species] = 0.;
  mp->d2_viscosity[MAX_VARIABLE_TYPES+species] = 0.;
  for (w = 0; w < pd->Num_Species; w++) {
    mp->d_FlowingLiquid_viscosity[MAX_VARIABLE_TYPES+w] = 0;
  }
  if(pd->gv[var])
    {
      if(alpha < alpha_g && alpha > 1e-8)
	{
	  deriv = (-A*pow(alpha, A-1)/alpha_g_pow)*(-B)*(1/ratio) ;
	  for (w = 0; w < pd->Num_Species; w++) {
	    mp->d_viscosity[MAX_VARIABLE_TYPES+species] = mu * deriv +
	      muL * mp->d_volumeFractionGas[MAX_VARIABLE_TYPES+w] * partial;
	  }
	  mp->d_FlowingLiquid_viscosity[MAX_VARIABLE_TYPES+species] = muL * deriv;
	}
    }

  return(status);
}


int
sylgard_viscosity(int species,    /* species number for cure equation */
		dbl mu0,        /* monomer reference temperature viscosity */
		dbl alpha_g,    /* extent of reaction at the gel point */
		dbl A,          /* exponent for constitutive equation */
		dbl Aexp)       /* exponent for thermal viscosity dependence */
{
  /* Local Variables */
  dbl mu;   /* viscosity */
  dbl alpha;  /* extent of reaction */
  dbl exponent;
  dbl ratio;
  dbl deriv;  /* stuff for the first derivative */
  dbl d_deriv; 
  dbl T;      /* Convenient local variables */
  int var;
  int status = 1;
  
  alpha = fv->c[species]; /* extent of reaction */
  exponent = A;

  if(alpha < alpha_g )
    {
      if(alpha > 0.0)
	{
	  ratio = (alpha_g)/(alpha_g - alpha);
	}
      else 
	{
	  ratio = 1.;
	}
    }
  else /* do something special at the gel point */
    {
      ratio = 100000;
    }
    
  
  if ( pd->e[pg->imtrx][TEMPERATURE] )
       {T = fv->T;}
  else
       {T = upd->Process_Temperature;}

  if(T <= 0.)
    {
      mu = mu0  * pow ( ratio, exponent );
    }
  else
    {
      mu = mu0 * exp (Aexp/T) * pow ( ratio, exponent );
    }
    
  mp->viscosity = mu;
  
  /* dmu_dT */
  var = TEMPERATURE;
  if(pd->v[pg->imtrx][var])
    {
      if( T <= 0.)
	{
	  mp->d_viscosity[var] 	= 0.; 
	  mp->d2_viscosity[var] = 0.; 
	}
      else
	{
	  mp->d_viscosity[var] = -mu * Aexp/(T*T) ;
	  mp->d2_viscosity[var] 	= -mu* SQUARE( Aexp/(T*T))
	    +2.* mu * Aexp/(T*T*T);
	}
    }
  
  /* dmu/dc */
  var  = MASS_FRACTION;
  if(pd->v[pg->imtrx][var])
    {
      if(alpha < alpha_g && alpha > 0.)
	{
	  deriv = exponent/(alpha_g - alpha);
	  d_deriv = exponent/SQUARE(alpha_g - alpha);
	  mp->d_viscosity[MAX_VARIABLE_TYPES+species] = mu * deriv;
	  mp->d2_viscosity[MAX_VARIABLE_TYPES+species] =  mu * deriv * deriv
	    + mu * d_deriv;
	}
      else
	{
	  mp->d_viscosity[MAX_VARIABLE_TYPES+species] = 0.;
	  mp->d2_viscosity[MAX_VARIABLE_TYPES+species] = 0.;
	}
	
    }

  return(status);
} /* end of sylgard_viscosity */


/*
 *
 *  filled_epoxy_viscosity model
 *
 */

/******************************************************************************
 *     Function that computes the viscosity of a filled epoxy that is a solid suspension
 *     whose behavior depends on the solid volume fraction, the extent of reaction and the
 *     temperature of the system via the relation: 
 *     
 *                     mu = mu0*(1. - C/Cmaxpack)^nexp
 *                          * alpha_g/(alpha-alpha_g)^(A+ B*alpha)
 *                          * exp(Aexp/T)
 *       where
 *                      mu0       = solvent viscosity
 *                      Cmaxpack  = binding solid fraction, i.e. solid fraction for infinite viscosity
 *                      nexp      = exponent
 *                      alpha     = extent of reaction
 *                      alpha_g   = extent of reaction at the gel point
 *                      A, B      = exponent for cure behavior
 *                      Aexp      = exponent for temperature dependence of viscosity
 *
 *     Function sets the viscosity members of mp.
 *
 *    Author: RRR
 *      Date: 3/24/97 
 *
 *
 *****************************************************************************/

int 
filled_epoxy_viscosity(int species_sus,	/* species num, solid volume fraction*/
		       int species_cur,	/* species num, extent of rxn        */
		       dbl mu0,	/* monomer reference viscosity               */
		       dbl maxpack, /* maximum solid volume fraction         */
		       dbl nexp, /* exponent for constitutive equation       */
		       dbl alpha_g, /* extent of reaction at the gel point   */
		       dbl A,	/* exponent for cure constitutive equation   */
		       dbl B,	/* exponent for cure constitutive equation   */
		       dbl T_g0, /* gelation temp unreacted filled_epoxy     */
		       dbl Atexp) /* exponent for thermal viscosity function */
{
  /* Local Variables */
  dbl mu=0;			/* viscosity */
  dbl mu0_local;
  dbl T;			/* temperature */
  dbl vf;			/* volume fraction solid */
  dbl ratio;
  dbl exponent, exponent_t;
  dbl alpha, alpha2, alpha_g2;
  dbl deriv, d_deriv, deriv_t;
  dbl uu, duu_dT=0, d2uu_dT2, deriv2_t2 ;
  dbl c1, c2, T_g, ln_10;
  int var;
  int status = 1;
  static int gelled=FALSE;

  /* initialize everything */
                                       
  mp->d_viscosity[TEMPERATURE] = 0.0;	
  mp->d2_viscosity[TEMPERATURE] = 0.0;

  mp->d_viscosity[MAX_VARIABLE_TYPES + species_sus] = 0.0;
  mp->d2_viscosity[MAX_VARIABLE_TYPES + species_sus] = 0.0;  
  mp->d_viscosity[MAX_VARIABLE_TYPES + species_cur] = 0.0;
  mp->d2_viscosity[MAX_VARIABLE_TYPES + species_cur] = 0.0;  

  vf = fv->c[species_sus]; /* volume fraction solid */

  alpha = fv->c[species_cur]; /* extent of reaction */
  alpha2 = alpha * alpha ;   /* extent of reaction squared */
  alpha_g2 = alpha_g * alpha_g;

  if(alpha < alpha_g)
    {
      ratio = 1./(1. - alpha2/alpha_g2);
      exponent = 4./3.;
    }
  else /* do something special at the gel point */
    {
      ratio = 100000;
      exponent = 4./3.;
      if (!gelled)
        {
          fprintf(stderr, "Viscosity has gelled: Now divergent");
          fprintf(stderr," The viscosity has reached gel-point and diverged.  Calculation stopped.");
          gelled = TRUE;
        }
      exit(-1);
    }
                                       
  if ( pd->e[pg->imtrx][TEMPERATURE] )
       {T = fv->T;}
  else
       {T = upd->Process_Temperature;}

  c1 = Atexp;
  c2 = B;
  T_g = T_g0/(1.- A * alpha);
  ln_10 = log(10.0);


  if(T <= 0.)
    {
      mu0_local = mu0 * pow ( ratio, exponent ); /*  mu0 should be mu(0,Tg) here*/
    }
  else
    {
      exponent_t = -c1*(T-T_g)/(c2+T-T_g);  /*CAR*/
      if(exponent_t > 20.) exponent_t=20.;
      mu0_local = mu0 * pow(10.0, exponent_t) * pow ( ratio, exponent );
      /*in above should use mu(0,Tg) rather than mu0 CAR */
    }

  if ( nexp > 0.0 || ( vf > 0.0 && vf <= (maxpack - .01) ) )
    {
      mu = mu0_local*pow( 1.0 - vf/maxpack, nexp );
      mp->viscosity = mu;
      mp->d_viscosity[MAX_VARIABLE_TYPES+species_sus] 
	= -mu*nexp/( maxpack - vf );
      mp->d2_viscosity[MAX_VARIABLE_TYPES+species_sus] 
	= mu0_local*nexp*(nexp-1.0)*pow( 1.0 - vf/maxpack,nexp-2.0 )/maxpack/maxpack ;
    }
  else if ( vf >=  (maxpack - .01) )
    {
      mu = mu0_local*pow( 0.01/maxpack, nexp );
      mp->viscosity = mu;
      mp->d_viscosity[MAX_VARIABLE_TYPES+species_sus] = 0.;
      mp->d2_viscosity[MAX_VARIABLE_TYPES+species_sus] = 0.;

    }
  else if ( vf <= 0. )
    {      
      mu = mu0_local;
      mp->viscosity = mu;
      mp->d_viscosity[MAX_VARIABLE_TYPES+species_sus] = 0.;
      mp->d2_viscosity[MAX_VARIABLE_TYPES+species_sus] = 0.;
    }

  /* dmu_dT */
  var = TEMPERATURE;
  if(Num_Var_In_Type[pg->imtrx][var])
    {
      if( T <= 0.)
	{
	  mp->d_viscosity[var] 	= 0. ;
	  deriv_t = 0.;
	  mp->d2_viscosity[TEMPERATURE] = 0.;
	}
      else
	{
	  uu = -c1*c2/pow(c2+T-T_g,2.0) ;
          duu_dT = ln_10*uu ;
          d2uu_dT2 = ln_10*(2*c1*c2*pow(c2+T-T_g,-3.0) + duu_dT*uu) ;
          deriv_t = duu_dT*mu ;
          deriv2_t2 = d2uu_dT2*mu ;
          mp->d_viscosity[var] = deriv_t ;
          mp->d2_viscosity[var] =  deriv2_t2 ;
	}
    }

  /* dmu/dc */
  var = MASS_FRACTION;
  if(Num_Var_In_Type[pg->imtrx][var])
    {
      if(alpha < alpha_g)
	{
	  deriv = 8./3.* alpha /(alpha_g2 - alpha2) 
	    - duu_dT * (A* T_g/(1.-A*alpha)) ;

	  d_deriv = 8./3.*(1./(alpha_g2 - alpha2) -2.*alpha2 /SQUARE(alpha_g2 - alpha2)) 
	    - duu_dT * A*A* T_g/SQUARE(1.-A*alpha);

	  mp->d_viscosity[MAX_VARIABLE_TYPES+species_cur] = mu * deriv;
	  mp->d2_viscosity[MAX_VARIABLE_TYPES+species_cur] = mu * deriv * deriv
	    + mu * d_deriv;
	}
      else
	{
	  mp->d_viscosity[MAX_VARIABLE_TYPES+species_cur] = 0.;
	  mp->d2_viscosity[MAX_VARIABLE_TYPES+species_cur] = 0.;
	}
    }

  return(status);
} /* end of filled_epoxy_viscosity */



/******************************************************************************************
 *     Function that computes the viscosity of a filling epoxy that is a solid suspension
 *     whose behavior depends on the solid volume fraction, the extent of reaction, and the
 *     temperature of the system via the relation: 
 *     
 *                     mu_L = mu0  alpha_g2/(alpha2-alpha_g2)^(4/3)  exp(Emu / RT)
 *       where
 *                      mu0       = solvent viscosity
 *                      alpha     = extent of reaction
 *                      alpha2    = alpha * alpha
 *                      alpha_g   = extent of reaction at the gel point
 *                      alpha_g2  = alpha_g * alpha_g
 *                      Emu       = exponent for thermal viscosity function
 *
 *       and
 *                      mu        = mu_L * exp(volFract/(1 - volFract))
 * 
 *
 *     Function sets the viscosity members of mp. It also calculates the jacobian
 *     contributions and sticks them in mp->d_viscosity
 *
 *    Author: RRR
 *      Date: 3/24/97 
 *
 * Args:
 *    int species_fluor;    species number for fluorinert tracking 
 *    int species_cur;      species number for extent of reaction tracking 
 *    dbl mu0;              monomer reference viscosity 
 *    dbl alpha_g;          extent of reaction at the gel point 
 *    dbl Emu;              exponent for thermal viscosity function
 *    dbl volFract          volume fraction
 *    dbl *volFract         dependence of volume fraction on independent unknowns
 *
 ******************************************************************************************/
int 
foam_epoxy_viscosity(int species_fluor, int species_cur, dbl mu0, 
		     dbl alpha_g, dbl Aexp)
{
  dbl mu;   /* viscosity */
  dbl T;    /* temperature */
  dbl ratio;
  dbl volFac = 1.0;
  dbl muL;
  dbl exponent;
  dbl alpha, alpha2, alpha_g2;
  dbl deriv;
  int j;
  int status = 1;
  dbl cVolFrac = 0.0;

  /* initialize everything */
                                       
  mp->d_viscosity[TEMPERATURE] = 0.0;	
  mp->d2_viscosity[TEMPERATURE] = 0.0;
  mp->d_FlowingLiquid_viscosity[TEMPERATURE] = 0.0;

  mp->d_viscosity[MAX_VARIABLE_TYPES + species_fluor] = 0.0;
  mp->d2_viscosity[MAX_VARIABLE_TYPES + species_fluor] = 0.0;  
  mp->d_viscosity[MAX_VARIABLE_TYPES + species_cur] = 0.0;
  mp->d2_viscosity[MAX_VARIABLE_TYPES + species_cur] = 0.0;  


  alpha = fv->c[species_cur]; /* extent of reaction */
  alpha2 = alpha*alpha;
  alpha_g2 = alpha_g*alpha_g;

  if (alpha < alpha_g)
    {
      ratio = 1./(1. - (alpha2/alpha_g2));
      exponent = 4./.3;
    }
  else /* do something special after the gel point */
    {
      ratio = 100000;
      exponent = 4./.3;
    }

  if (pd->e[pg->imtrx][TEMPERATURE]) {
    T = fv->T;
  } else {
    T = upd->Process_Temperature;
  }

  // Calculate the pure liquid-phase viscosity
  muL = mu0 * pow ( ratio, exponent )* exp(Aexp/T);

  if (mp->volumeFractionGas >= 0.0) 
    {
      cVolFrac = mp->volumeFractionGas;
      if (cVolFrac > 0.99) cVolFrac = 0.99;
      volFac = exp(cVolFrac / (1.0 - cVolFrac));
    }

  mu = muL * volFac;

  /* dmu/dc */
  if (Num_Var_In_Type[pg->imtrx][MASS_FRACTION])
    {
      mp->d_viscosity[MAX_VARIABLE_TYPES+species_fluor] = 0.0;
      mp->d_FlowingLiquid_viscosity[MAX_VARIABLE_TYPES+species_fluor] = 0.0;
    }
    

  /* dmu_dT */
  if (Num_Var_In_Type[pg->imtrx][TEMPERATURE])
    {
      mp->d_viscosity[TEMPERATURE]  = -mu * Aexp/(T*T);
      mp->d_FlowingLiquid_viscosity[TEMPERATURE] = -muL * Aexp/(T*T);
    }

  /* dmu/dc */
  if (Num_Var_In_Type[pg->imtrx][MASS_FRACTION])
    {
      if (alpha < alpha_g)
	{
	  deriv = 2.*exponent*alpha/(alpha_g2 - alpha*alpha);
	  mp->d_viscosity[MAX_VARIABLE_TYPES+species_cur] = mu * deriv;
	  mp->d_FlowingLiquid_viscosity[MAX_VARIABLE_TYPES+species_cur]
	    = muL * deriv;
	}
      else
	{
	  mp->d_viscosity[MAX_VARIABLE_TYPES+species_cur] = 0.0;
	  mp->d_FlowingLiquid_viscosity[MAX_VARIABLE_TYPES+species_cur] = 0.0;
	}

      if (mp->volumeFractionGas >= 0.0 && 
	  mp->volumeFractionGas <= 0.99 ) 
	{
	  dbl volFacDeriv = mu / ((1-cVolFrac) * (1-cVolFrac));
	  for (j = 0; j < Num_Var_In_Type[pg->imtrx][MASS_FRACTION]; j++) 
	    {
	      mp->d_viscosity[MAX_VARIABLE_TYPES + j]
		+= volFacDeriv * mp->d_volumeFractionGas[MAX_VARIABLE_TYPES + j];
	    }
	}
    }

  mp->FlowingLiquid_viscosity = muL;
  mp->viscosity = mu;

  return(status);
} /* end of foam_epoxy_viscosity */

/******************************************************************************
 *     Function that computes the viscosity of a solution  whose behavior depends on the 
 *     temperature of the system via the relation: 
 *     
 *                     mu = mu0*exp(Aexp/T)
 *       where
 *                      mu0       = reference temperature viscosity
 *                      Aexp      = exponent for temperature dependence of viscosity
 *
 *     Function sets the viscosity members of mp.
 *
 *    Author: RRR
 *      Date: 3/24/97
 *
 *
 *****************************************************************************/
int
thermal_viscosity(dbl mu0,	/* reference temperature fluid viscosity */
		  dbl Aexp)	/* exponent for constitutive equation */
{
  /* Local Variables */
  dbl mu;   /* viscosity */

  dbl T; /* Convenient local variables */
  int status = 1;

  if (! pd->v[pg->imtrx][TEMPERATURE] )
    {
      return(0);
    }

                                       
  if ( pd->e[pg->imtrx][TEMPERATURE] )
       {T = fv->T;}
  else
       {T = upd->Process_Temperature;}

  if( T <= 0.3)
    {
      mu = mu0*exp(Aexp/0.3);
      mp->viscosity = mu;
      mp->d_viscosity[TEMPERATURE] 	= 0. ;
    }
  else
    {
      mu = mu0 * exp (Aexp/T);
      mp->viscosity = mu;
      mp->d_viscosity[TEMPERATURE] = -mu0 * Aexp/(T*T) ;
    }

  return(status);
} /* end of thermal_viscosity */


/*
 *
 *  cure viscosity model
 *
 */

/******************************************************************************
 *     Function that computes the viscosity as a function of the extent of reaction 
 *     system via the relation: 
 *     
 *                     mu = mu0 * alpha_g/(alpha-alpha_g)^(A+ B*alpha)
 *       where
 *                      mu0       = monomer viscosity
 *                      alpha     = extent of reaction
 *                      alpha_g   = extent of reaction at the gel point
 *                      A, B      = exponent for cure behavior
 *
 *     Function sets the viscosity members of mp.
 *
 *    Author: RRR
 *      Date: 3/24/97
 *
 *
 *****************************************************************************/



int
cure_viscosity(int species,	/* species num, solid volume fraction        */
	       dbl mu0,		/* carrier fluid viscosity                   */
	       dbl alpha_g,     /* extent of reaction at the gel point       */
	       dbl A,		/* exponent for constitutive equation        */
	       dbl B)		/* exponent for constitutive equation        */
{
  dbl mu;     /* viscosity */
  dbl alpha;  /* extent of reaction */
  dbl exponent;
  dbl ratio;
  dbl deriv;  /* stuff for the first derivative */
  int status = 1;

  if (! pd->v[pg->imtrx][MASS_FRACTION] )
    {
      return(0);
    }

  alpha = fv->c[species]; /* extent of reaction */
  if(alpha < alpha_g)
    {
      ratio = (alpha_g)/(alpha_g - alpha);
      exponent = A + B * alpha;
    }
  else /* do something special at the gel point */
    {
      ratio = 100000;
      exponent = A + B*alpha_g;
    }

  mu = mu0 * pow ( ratio, exponent );

  mp->viscosity = mu;
  
  /* dmu/dc */

  if(alpha < alpha_g)
    {
      deriv = exponent/(alpha_g - alpha) + B*log(ratio);
      mp->d_viscosity[MAX_VARIABLE_TYPES+species] = mu * deriv;
    }
  else
    {
      mp->d_viscosity[MAX_VARIABLE_TYPES+species] = 0.;
    }
  
  return(status);
} /* end of cure_viscosity */

/******************************************************************************************
 *     Function that computes the viscosity of a solution  whose behavior depends on the 
 *     structure formation of of the system via the relation: 
 *     
 *                     mu = mu_inf + mu0*nn**x
 *       where
 *                      mu_inf    = plateau viscosity at high shear
 *                      mu0       = reference  viscosity when x is zero
 *                      Aexp      = exponent for bond dependence of viscosity
 *
 *     Function sets the viscosity members of mp.
 *
 *    Author: RRR
 *      Date: 3/24/97
 *
 *
 ******************************************************************************************/


int bond_viscosity(dbl mu0,         /* reference zero shear rate fluid viscosity */
		   dbl mu_inf,     /* reference high shear rate fluid viscosity */
		   dbl Aexp)        /* exponent for constitutive equation */
{
  /* Local Variables */
  dbl mu;   /* viscosity */

  dbl nn; /* Convenient local variables */
  int status = 1;

  if (! pd->v[pg->imtrx][BOND_EVOLUTION] )
    {
      return(0);
    }

                                       
  nn= fv->nn;
  if( nn <= 0.0)
    {
      mu = mu_inf;
      mp->viscosity = mu;
      mp->d_viscosity[BOND_EVOLUTION] 	= 0. ;
    }
  else
    {
      mu = mu_inf + mu0 * pow(nn, Aexp);
      mp->viscosity = mu;
      mp->d_viscosity[BOND_EVOLUTION] = mu0 *Aexp*pow(nn, Aexp-1.);
    }

  return(status);
} /* end of bond_viscosity */


double
carreau_wlf_conc_viscosity(struct Generalized_Newtonian *gn_local,
 		  dbl gamma_dot[DIM][DIM], /* strain rate tensor */
 		  VISCOSITY_DEPENDENCE_STRUCT *d_mu,
		  const int const_eqn)
{
 
  int a, b;
 
  int var;
  int mdofs=0,vdofs;
 
  int i, j, w;
 
  dbl gammadot;	                /* strain rate invariant */ 
 
  dbl d_gd_dv[DIM][MDE];        /* derivative of strain rate invariant 
 				   wrt velocity */ 
  dbl d_gd_dmesh[DIM][MDE];     /* derivative of strain rate invariant 
 				   wrt mesh */ 
 
  dbl val1;
  dbl shear, visc_cy;
  dbl d_shear, d_visc_cy, d_shear_d_at,d_shear_d_ac;
  dbl mu = 0.;
  dbl mu0;
  dbl muinf;
  dbl dmudT = 0.0, dmudC = 0.0;
  dbl nexp;
  dbl atexp;
  dbl aexp;
  dbl at_shift,at_conc = 0.0;
  dbl lambda;
  dbl wlf_denom;
  dbl wlfc2;
  dbl ref_conc;
  int conc_exp;
  dbl nonvol_conc = 0.0;
  dbl temp;
 
  calc_shearrate(&gammadot, gamma_dot, d_gd_dv, d_gd_dmesh);
 
  mu0 = gn_local->mu0;
  nexp = gn_local->nexp;
  muinf = gn_local->muinf;
  aexp = gn_local->aexp;
  atexp = gn_local->atexp;
  wlfc2 = gn_local->wlfc2;
  lambda = gn_local->lam;
  ref_conc = gn_local->maxpack;
  conc_exp = gn_local->fexp;
   
 
  vdofs = ei[pg->imtrx]->dof[VELOCITY1];
   
  if ( pd->e[pg->imtrx][R_MESH1] )
     {
       mdofs = ei[pg->imtrx]->dof[R_MESH1];
     }
 
 /*  temperature shift factor  */
 
  if ( pd->e[pg->imtrx][TEMPERATURE] )
       {temp = fv->T;}
  else
       {temp = upd->Process_Temperature;}

   wlf_denom = wlfc2 + temp - mp->reference[TEMPERATURE];
   if(wlf_denom != 0.)
 	{
       at_shift=exp(atexp*(mp->reference[TEMPERATURE]-temp)/wlf_denom);
       if(!isfinite(at_shift)) { at_shift = DBL_MAX; }
 	}
   else
     {
       at_shift = 1.;
     }
 
 /*  concentration shift factor  */
 
	switch (mp->Species_Var_Type)
 		{
 		case SPECIES_MASS_FRACTION:
    			nonvol_conc = 1.;
    			for(w=0 ; w<pd->Num_Species_Eqn ; w++)
  				{ nonvol_conc -= fv->c[w]; }
 			break;
 		case SPECIES_DENSITY:
    			nonvol_conc = 1.;
    			for(w=0 ; w<pd->Num_Species_Eqn ; w++)
  				{ nonvol_conc -= fv->c[w]*mp->specific_volume[w]; }
 			nonvol_conc /= mp->u_density[0];
 			break;
 		default:
 			EH(-1,"That species type not completed yet.");
 			break;
 		}
   if( nonvol_conc < 0)
 	{
	fprintf(stderr,"nonvolatile conc %g %g \n",nonvol_conc,fv->c[0]);
 	if(const_eqn == CARREAU_WLF_CONC_PL) nonvol_conc = 0.;
  	WH(-1,"negative nonvolatile concentration");
 	}
   if(const_eqn == CARREAU_WLF_CONC_PL)
	{
    	at_conc = pow(nonvol_conc/ref_conc,conc_exp);
 	}
   else if( const_eqn == CARREAU_WLF_CONC_EXP)
	{
   	at_conc = exp(conc_exp*(nonvol_conc-ref_conc));
	}
   else
	{
  	EH(-1,"invalid constitutive model for WLF_CONC");
	}
   
   if(gammadot != 0.)
     {
       shear = pow( at_conc*at_shift*lambda*gammadot, aexp);
       val1 = pow( at_conc*at_shift*lambda*gammadot, aexp-1.);
       d_shear = aexp * at_conc * at_shift * lambda * val1;
       d_shear_d_at = aexp * at_conc * gammadot * lambda * val1; 
       d_shear_d_ac = aexp * at_shift * gammadot * lambda * val1; 
     }
   else
     {
       shear = 0.;
       d_shear = 0.;
       d_shear_d_at = 0.;
       d_shear_d_ac = 0.;
     }
   
   
   visc_cy = pow(1. + shear,(nexp-1.)/aexp);
   d_visc_cy = (nexp-1.)/aexp * pow(1. + shear,(nexp-1.-aexp)/aexp);
   
   mu = at_conc * at_shift * (muinf + (mu0 - muinf) * visc_cy);
   
   if ( d_mu != NULL ) d_mu->gd =  at_conc * at_shift * ((mu0 - muinf) * d_visc_cy *d_shear);
   
   /*
    * d( mu )/dT
    */
   
   var = TEMPERATURE;
   if ( d_mu != NULL && pd->e[pg->imtrx][var] )
     {
       if(wlf_denom != 0.)
 	{
 	  dmudT = mu +
 	  at_conc*at_conc*at_shift*at_shift*(mu0 - muinf)*d_visc_cy*d_shear_d_at;
 	  dmudT *= - atexp*wlfc2/(wlf_denom*wlf_denom);
	  if(!isfinite(dmudT)) { dmudT = DBL_MAX; }
 	}
       else
 	{
 	  dmudT = 0.;
 	}
       
       for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
 	{
 	  d_mu->T[j]= dmudT * bf[var]->phi[j];
 	}
     }
   /*
    * d( mu )/dC
    */
   
   var = MASS_FRACTION;
   if ( d_mu != NULL && pd->v[pg->imtrx][var] )
     {
   	if(const_eqn == CARREAU_WLF_CONC_PL)
		{
       		if(nonvol_conc > 0.)
 			{
 	  		dmudC = mu +
 	  			at_conc*at_conc*at_shift*at_shift*(mu0 - muinf)
				*d_visc_cy*d_shear_d_ac;
 	  		dmudC *= conc_exp/nonvol_conc;
 			}
       		else
 			{ dmudC = 0.; }
		}
   	else if( const_eqn == CARREAU_WLF_CONC_EXP)
		{
 	  	dmudC = mu +
 	  		at_conc*at_conc*at_shift*at_shift*(mu0 - muinf)
			*d_visc_cy*d_shear_d_ac;
 	  	dmudC *= conc_exp;
		}
   	else
		{
  		EH(-1,"invalid constitutive model for WLF_CONC");
		}
       
       for(w=0 ; w<pd->Num_Species_Eqn ; w++)
 	{
       	   for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
 		{
 	  	d_mu->C[w][j]= -dmudC * bf[var]->phi[j];
 		}
 	}
     }
   
   
   /*
    * d( mu )/dmesh
    */
   if ( d_mu != NULL && pd->e[pg->imtrx][R_MESH1] )
     {
       for ( b=0; b<VIM; b++)
 	{
 	  for ( j=0; j<mdofs; j++)
 	    {
 	      if(gammadot != 0.0 && Include_Visc_Sens )
 		{
 		  d_mu->X [b][j] =
 		    d_mu->gd * d_gd_dmesh [b][j] ;
 		}
 	      else
 		{
 		  /* printf("\ngammadot is zero in viscosity function");*/
 		  d_mu->X [b][j] = 0.0;
 		}
 	    }
 	}
     }
   
  /*
   * d( mu )/dv
   */
  if ( d_mu != NULL && pd->e[pg->imtrx][R_MOMENTUM1] )
    {
      for ( a=0; a<VIM; a++)
        {
          for ( i=0; i<vdofs; i++)
	    {
	      if(gammadot != 0.0 && Include_Visc_Sens )
	        {
	          d_mu->v[a][i] =
		    d_mu->gd * d_gd_dv[a][i] ;
	        }
	      else
	        {
	          d_mu->v[a][i] = 0.0 ;
	        }
	    }
        }
    }
 
   return(mu);
 }

int
ls_modulate_viscosity ( double *mu1,
			double  mu2,
			double width,
			double pm_minus,
			double pm_plus,
                        VISCOSITY_DEPENDENCE_STRUCT *d_mu,
                        const int model )
{
  double factor, ratio=0.0;
  int i,a, w, var;

  if (model == RATIO)
      {
      ratio = mu2;
      mu2 = *mu1*ratio;
      }
  if ( d_mu == NULL )
    {
      *mu1 = ls_modulate_property( *mu1, mu2, width, pm_minus, pm_plus, NULL, &factor);
      return(1);
    }

  *mu1 = ls_modulate_property( *mu1, mu2, width, pm_minus, pm_plus, d_mu->F, &factor);

  if (model == RATIO)
      {
      factor *= (1.-ratio);
      factor += ratio;
      }

  d_mu->gd *= factor;

  if ( pd->v[pg->imtrx][var=TEMPERATURE ] )
    {
      for(i=0; i<ei[pg->imtrx]->dof[var]; i++)
	{
	  d_mu->T[i] *= factor;
	}
    }
      
  if ( pd->v[pg->imtrx][var=MASS_FRACTION ] )
    {
      for ( w=0; w<pd->Num_Species_Eqn; w++)
	{
	  for( i=0; i<ei[pg->imtrx]->dof[var]; i++) 
	    {
	      d_mu->C[w][i] *= factor;
	    }
	}
    }

  if( pd->v[pg->imtrx][var=VELOCITY1] )
    {
      for( a=0; a<pd->Num_Dim; a++ )
	{
	  for( i=0; i<ei[pg->imtrx]->dof[var]; i++ )
	    {
	      d_mu->v[a][i] *= factor;
	    }
	}
    }

  if( pd->v[pg->imtrx][var=MESH_DISPLACEMENT1] )
    {
      for( a=0; a<pd->Num_Dim; a++ )
	{
	  for( i=0; i<ei[pg->imtrx]->dof[var]; i++ )
	    {
	      d_mu->X [a][i] *= factor;
	    }
	}
    }

  if( pd->v[pg->imtrx][var=PRESSURE] )
    {
      for( i=0; i<ei[pg->imtrx]->dof[var]; i++ )
	{
	  d_mu->P[i] *= factor;
	}
    }

  if ( pd->v[pg->imtrx][var=BOND_EVOLUTION] )
    {
      for( i=0 ; i<ei[pg->imtrx]->dof[var]; i++)
	{
	  d_mu->nn[i] *= factor;
	}
    }
  return ( 1 );
  
}
      

void
copy_pF_to_F ( int phase )
{
  int a;

  fv->F = fv->pF[phase];

  for( a=0; a<pd->Num_Dim; a++ ) fv->grad_F[a] = fv->grad_pF[phase][a];

}

/*******************************************************************************
 * flowing_liquid_viscosity(): Calculate the flowing liquid viscosity used in
 *                             Brinkman term of momentum equation together with
 *                             its sensitivities with respect to solution unknowns
 *                             at the Gauss point. It typically depends on species
 *                             concentration, temperature, etc
 *
 * Input
 *----------
 *
 * Output
 * -----
 *
 *   flow_vis    = flowing liquid viscosity
 *   d_flow_vis  = dependence of flowing liquid viscosity on the independent unknowns
 *                 in the local element stiffness matrix.
 *
 *
 *
 *******************************************************************************/
double
flowing_liquid_viscosity(VISCOSITY_DEPENDENCE_STRUCT *d_flow_vis)
{

  int err;
  int var, var_offset, vdofs;
  int j, a, w;
  int dim = ei[pg->imtrx]->ielem_dim;

  double flow_vis = 0.;

  /* Zero out sensitivities */
  if (d_flow_vis != NULL)
    {
     zeroStructures(d_flow_vis, 1);
    }
  if (mp->PorousMediaType != POROUS_BRINKMAN)
        WH(-1, "Set Porous term multiplier in continuous medium");


  /***** Evaluate FlowingLiquid Viscosity based on the specified model ****/

  if (mp->FlowingLiquidViscosityModel == CONSTANT)
     {
      flow_vis = mp->FlowingLiquid_viscosity;
      mp_old->FlowingLiquid_viscosity = flow_vis;
     }

  else if (mp->FlowingLiquidViscosityModel == MOLTEN_GLASS)
     {
      (void) molten_glass_viscosity(&flow_vis,
                                    d_flow_vis->T, mp->u_FlowingLiquid_viscosity);
     }

  else if (mp->FlowingLiquidViscosityModel == EPOXY)
     {
      (void) epoxy_flowing_liquid_viscosity(&flow_vis, d_flow_vis, mp->u_FlowingLiquid_viscosity);
     }

  else if (mp->FlowingLiquidViscosityModel == USER)
     {
      (void) usr_FlowingLiquidViscosity(mp->u_FlowingLiquid_viscosity);
      flow_vis = mp->FlowingLiquid_viscosity;
      mp_old->FlowingLiquid_viscosity = flow_vis;

      if (d_flow_vis != NULL)
        {
	 if (pd->v[pg->imtrx][TEMPERATURE] )
          {
           var = TEMPERATURE;
           for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
            {
             d_flow_vis->T[j]= mp->d_FlowingLiquid_viscosity[var]*bf[var]->phi[j];
            }
          }
	 if (pd->v[pg->imtrx][VELOCITY1] )
	   {
	    var = VELOCITY1;
	    vdofs = ei[pg->imtrx]->dof[var];
	    for ( a=0; a<dim; a++)
	      {
	       for ( j=0; j<vdofs; j++)
		  {
		   d_flow_vis->v[a][j] = mp->d_FlowingLiquid_viscosity[var+a]*bf[var]->phi[j];
		  }
	      }
	   }
         if (pd->v[pg->imtrx][MESH_DISPLACEMENT1])
          {
	   var = MESH_DISPLACEMENT1;
	   for ( a=0; a<dim; a++)
	      {
               for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
                {
                d_flow_vis->X[a][j]= mp->d_FlowingLiquid_viscosity[var+a]*bf[var]->phi[j];
                }
              }
          }
	 if (pd->v[pg->imtrx][MASS_FRACTION] )
	    {
	     for ( w=0; w<pd->Num_Species_Eqn; w++)
	       {
		var	     = MASS_FRACTION;
		var_offset = MAX_VARIABLE_TYPES + w;
		for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		  {
		   d_flow_vis->C[w][j] =mp->d_FlowingLiquid_viscosity[var_offset]*bf[var]->phi[j];
		  }
	       }
	    }
	}
     }

  else
     {
      EH(-1,"Don't recognize your FlowingLiquidViscosity model");
     }


  if (ls != NULL && 
      mp->mp2nd != NULL &&
      (mp->mp2nd->FlowingLiquidViscosityModel == CONSTANT )
    )
    {
      err= ls_modulate_viscosity(&flow_vis, mp->mp2nd->FlowingLiquid_viscosity, ls->Length_Scale,
                                 (double) mp->mp2nd->FlowingLiquid_viscositymask[0],
                                 (double) mp->mp2nd->FlowingLiquid_viscositymask[1],
                                 d_flow_vis, mp->mp2nd->FlowingLiquidViscosityModel );
      EH(err, "ls_modulate_viscosity");
    }


  return(flow_vis);

} /* End of flowing_liquid_viscosity*/

/* end of file mm_viscosity.c */
