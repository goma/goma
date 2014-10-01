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
 *$Id: user_mp.c,v 5.5 2009-01-12 16:41:49 hkmoffa Exp $
 */

#ifdef USE_RCSID
static char rcsid[] = "$Id: user_mp.c,v 5.5 2009-01-12 16:41:49 hkmoffa Exp $";
#endif

/* Standard include files */

#include <stdlib.h>
#include <stdio.h>
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
#include "rf_fill_const.h"
#include "rf_vars_const.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"

#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "mm_mp.h"

#include "mm_eh.h"

#define _USER_MP_C
#include "goma.h"

/*********** R O U T I N E S   I N   T H I S   F I L E *************************
*
*       NAME            TYPE            CALLED_BY
*    ------------             ---------               --------------
*
*    usr_thermal_conductivity    ()  int      assemble_energy
*    usr_electrical_conductivity ()  int      joule_heat_source
*    usr_density          ()         int      density
*    usr_heat_capacity    ()         int      assemble_energy
*    usr_heat_source      ()         int      assemble_energy
*    usr_species_source   ()         int      get_continuous_species_terms
*    usr_current_source   ()         int      assemble_potential
*    usr_viscosity        ()         int      viscosity
*    usr_surface_tension  ()         int      load_surface_tension
*    usr_momentum_source  ()         int      momentum_source_term
*    usr_lame_mu          ()         int      load_elastic_properties
*    usr_lame_lambda      ()         int      load_elastic_properties
*    usr_diffusivity      ()         int      Diffusivity
*    usr_FlowingLiquidViscosity ()   int      assemble_momentum & _continuity
*    usr_heat_flux        ()         int      assemble_energy
*
*******************************************************************************/
/*
*               ******SPECIAL RESTRICTION for these routines*******
* NB: For now (3/2/95) we will not allow any dependencies on gradients, a fairly
*     severe restriction.  I you are interested in that see the "viscosity" routine
*     in mm_fill_terms.c, where RRRao has taken such considerations into account.
*
*     There is a set of complimentary routines which allow for
*     dependencies on gradients. See file "user_mp_gen.c"  These correspond 
*     to the USER_GEN option for the property model. Use of those routines 
*     will undoubtedly be  more complex.
*/  

/*********** R E C I P E   F O R   U S A G E **********************************
* Each routine in this file correspond to one material property.  Each of the 
* routines is responsible for:
*
*     (1)Calculating the value of that material property at the current gauss 
*        integration point.    
*     (2)Calculating the value of all derivatives with respect to all degrees
*        of freedom REQUESTED by the subroutine at each gauss point.
*
* If these functions are not performed, then the routine should stop and return
* an error.  If the routine is used, then the call to  EH (error handler routine) routine
* should be removed.  The following is a well-documented example for the thermal
* conductivity property that presumably depends on everything but the phase of the
* moon and, by the way gradients of anything:
*
*******************************************************************************/
/*
 * THERMAL CONDUCTIVITY 
 */
/*
 * int usr_thermal_conductivity (param)
 *
 * ------------------------------------------------------------------------------
 * This routine is responsible for filling up the following portions of the mp structure
 * at the current gauss point:
 *     intput:    param - array of constants input on the property card.  
 *
 *     output:   k   => mp->thermal_conductivity - thermal conductivity at 
 *                                                  current gauss point
 *              dkdT => mp->d_thermal_conductivity[TEMPERATURE] 
 *                                         - derivative wrt temperature.
 *            dkdC[i]=> mp->d_thermal_conductivity[MASS_FRACTION][i]
 *                                         - derivative wrt mass frac species i
 *            dkdX[0]=> mp->d_thermal_conductivity[MESH_DISPLACEMENT1]
 *            dkdX[1]=> mp->d_thermal_conductivity[MESH_DISPLACEMENT2]
 *            dkdX[2]=> mp->d_thermal_conductivity[MESH_DISPLACEMENT3]
 *                                         - derivative wrt mesh displacements
 *
 *   NB: The user need only supply k, dkdT, dkdC, etc....mp struct is loaded up for you
 * --------------------------------------------------------------------------------
 * Example of an everything-dependent thermal conductivity for material "sample"
 * -----------Add these lines just before last section.-------------
 *
 * Simple conductivity function for one material 
 *
 *  if (!strcmp(pd->MaterialName, "sample") )   
 *    {
 *	k = param[0]*T + param[1]*X[0] + param[2]*X[1] + param[3]*C[0] + param[4];
 *     dkdT    = param[0];
 *     dkdX[0] = param[1];
 *     dkdX[1] = param[2];
 *     dkdC[0] = param[3];
 *    } 
 *  else
 *    {
 *	EH(-1,"No user-defined function for this material");
 *    }
 * ----------------------------------------------------------------------------
 */

int
usr_thermal_conductivity(dbl *param, dbl time) /* user-defined parameter list */
{
  int a;

  dbl k, dkdT;        /* thermal conductivity and its derivative wrt temperature*/
  dbl dkdX[DIM];      /* thermal conductivity derivative wrt displacement*/
  dbl dkdC[MAX_CONC]; /* thermal conductivity derivative wrt concentration*/

  dbl X[DIM], T, C[MAX_CONC]; /* Convenient local variables */

  int i;

  /* Begin Execution */
 /**********************************************************/

 /* Comment out or remove this line if using this routine */

 /*EH(-1,"No user_defined thermal conductivity model implemented");*/

 /**********************************************************/

 
 /************Initialize everything for safety**************/
  k = 0.;                                   /*Do not touch */
  dkdT = 0.;                                /*Do not touch */
  for (i=0; i<DIM; i++)                     /*Do not touch */
    {                                       /*Do not touch */
      dkdX[i]=0.;                           /*Do not touch */
    }                                       /*Do not touch */
  for (i=0; i<MAX_CONC; i++)                /*Do not touch */
    {                                       /*Do not touch */
      dkdC[i] = 0.;                         /*Do not touch */
    }                                       /*Do not touch */
 /**********************************************************/

 /***********Load up convenient local variables*************/
 /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T;                                             /*Do not touch */
  for(a=0; a<DIM; a++)X[a] = fv->x[a];		         /*Do not touch */
  for(i=0; i<pd->Num_Species_Eqn; i++) C[i] = fv->c[i];  /*Do not touch */

 /**********************************************************/

 /*******Add property function and sensitivities here*******/

   
 /****************Don't touch these lines***********************/ 
  mp->thermal_conductivity = k; 			      /*Do not touch */
  mp->d_thermal_conductivity[TEMPERATURE] = dkdT;	      /*Do not touch */
  for ( a=0; a<DIM; a++)				      /*Do not touch */
    {							      /*Do not touch */
      mp->d_thermal_conductivity[MESH_DISPLACEMENT1+a] = dkdX[a]; /*Do not touch */
    }							      /*Do not touch */
  for ( a=0; a<MAX_CONC; a++)				      /*Do not touch */
    {							      /*Do not touch */
      mp->d_thermal_conductivity[MAX_VARIABLE_TYPES+a] = dkdC[a]; /*Do not touch */
    }
 /**************************************************************/
  return(0);
} /* End of usr_thermal_conductivity */
/*****************************************************************************/

int
usr_electrical_conductivity(dbl *param, dbl time)	/* user-defined parameter list */
{
  int a;

  dbl k, dkdT;        /* electrical conductivity and its derivative wrt temperature*/
  dbl dkdV;           /* electrical conductivity derivative wrt voltage*/
  dbl dkdX[DIM];      /* electrical conductivity derivative wrt displacement*/
  dbl dkdC[MAX_CONC]; /* electrical conductivity derivative wrt concentration*/

  dbl X[DIM], T, C[MAX_CONC], V; /* Convenient local variables */

  int i;

  /* Begin Execution */
 /**********************************************************/

 /* Comment out our remove this line if using this routine */

    EH(-1,"No user_defined electrical conductivity model implemented"); 
 /**********************************************************/
 
 /************Initialize everything for safety**************/
  k = 0.;                                   /*Do not touch */
  dkdT = 0.;                                /*Do not touch */
  dkdV = 0.;                                /*Do not touch */
  for (i=0; i<DIM; i++)                     /*Do not touch */
    {                                       /*Do not touch */
      dkdX[i]=0.;                           /*Do not touch */
    }                                       /*Do not touch */
  for (i=0; i<MAX_CONC; i++)                /*Do not touch */
    {                                       /*Do not touch */
      dkdC[i] = 0.;                         /*Do not touch */
    }                                       /*Do not touch */
 /**********************************************************/

 /***********Load up convenient local variables*************/
 /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T;                                      /*Do not touch */
  V = fv->V;                                      /*Do not touch */
  for(a=0; a<DIM; a++)X[a] = fv->x[a];		  /*Do not touch */
  for(i=0; i<pd->Num_Species_Eqn; i++) C[i] = fv->c[i];  /*Do not touch */

 /**********************************************************/

 /*******Add property function and sensitivities here*******/


 /****************Don't touch these lines***********************/ 
  mp->electrical_conductivity = k; 			      /*Do not touch */
  mp->d_electrical_conductivity[TEMPERATURE] = dkdT;	      /*Do not touch */
  mp->d_electrical_conductivity[VOLTAGE] = dkdV;	      /*Do not touch */
  for ( a=0; a<DIM; a++)				      /*Do not touch */
    {							      /*Do not touch */
      mp->d_electrical_conductivity[MESH_DISPLACEMENT1+a] = dkdX[a];  /*Do not touch */
    }							      /*Do not touch */
  for ( a=0; a<MAX_CONC; a++)				      /*Do not touch */
    {							      /*Do not touch */
      mp->d_electrical_conductivity[MAX_VARIABLE_TYPES+a] = dkdC[a];  /*Do not touch */
    }
 /**************************************************************/

  return(0);
} /* End of usr_electrical_conductivity */
/*****************************************************************************/

int 
usr_density(dbl *param)         /* pointer to user-defined parameter list    */
{
  /* Local Variables */
  dbl rho, d_rho_dT;        /* density and its derivative wrt temperature*/
  dbl d_rho_dC[MAX_CONC];   /* density derivative wrt concentration*/
  dbl T, C[MAX_CONC]; /* Convenient local variables */
  int w;

  /* Begin Execution */
 /**********************************************************/

 /* Comment out or remove this line if using this routine */

   EH(-1,"No user_defined density model implemented"); 
 
 /**********************************************************/
 
 /************Initialize everything for safety**************/
  rho = 0.;                                   /*Do not touch */
  d_rho_dT = 0.;                            /*Do not touch */
  for (w=0; w<MAX_CONC; w++)                /*Do not touch */
    {                                       /*Do not touch */
      d_rho_dC[w] = 0.;                     /*Do not touch */
    }                                       /*Do not touch */
 /**********************************************************/

 /***********Load up convenient local variables*************/
 /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T;                                      /*Do not touch */
  for(w=0; w<pd->Num_Species_Eqn; w++) C[w] = fv->c[w];  /*Do not touch */

 /**********************************************************/

 /*******Add property function and sensitivities here*******/


 /****************Don't touch these lines***********************/ 
  mp->density = rho; 				              /*Do not touch */
  mp->d_density[TEMPERATURE] = d_rho_dT;	              /*Do not touch */
  for ( w=0; w<MAX_CONC; w++)				      /*Do not touch */
    {							      /*Do not touch */
      mp->d_density[MAX_VARIABLE_TYPES+w] = d_rho_dC[w];      /*Do not touch */
    }
 /**************************************************************/

  return(0);
} /* End of usr_density */
/*****************************************************************************/

/*
 * HEAT CAPACITY 
 */
/*
 * int usr_heat_capacity ()
 *
 * ----------------------------------------------------------------------------
 * This routine is responsible for filling up the following portions of the
 * mp structure at the current gauss point:
 *     intput:    param - array of constants input on the property card.  
 *
 *     output: Cp       => mp->heat_capacity - heat capacity at current gauss point
 *             dCpdT    => mp->d_heat_capacity[TEMPERATURE] 
 *                                         - derivative wrt temperature.
 *             dCpdC[i] => mp->d_heat_capacity[MASS_FRACTION][i]
 *                                         - derivative wrt mass frac species i
 *             dCpdX[0] => mp->d_heat_capacity[MESH_DISPLACEMENT1]
 *             dCpdX[1] => mp->d_heat_capacity[MESH_DISPLACEMENT2]
 *             dCpdX[2] => mp->d_heat_capacity[MESH_DISPLACEMENT3]
 *                                         - derivative wrt mesh displacements
 *
 * 
 *
 *   NB: The user need only supply Cp, dCpdT, dCpdC, etc....mp struct is loaded up for you
 * ----------------------------------------------------------------------------
 * Example of an everything-dependent heat capacity for material "sample"
 * -----------Add these lines just before last section.-------------
 *
 * Simple conductivity function for one material 
 *
 *  if (!strcmp(pd->MaterialName, "sample") )   
 *    {
 *	Cp = param[0]*T + param[1]*X[0] + param[2]*X[1] + param[3]*C[0] + param[4];
 *     dCpdT    = param[0];
 *     dCpdX[0] = param[1];
 *     dCpdX[1] = param[2];
 *     dCpdC[0] = param[3];
 *    } 
 *  else
 *    {
 *	EH(-1,"No user-defined function for this material");
 *    }
 * ----------------------------------------------------------------------------
 */
int
usr_heat_capacity(dbl *param, dbl time)	/* pt to user-defined parameter list */
{
  int a;

  dbl Cp, dCpdT;  /* heat capacity and its derivative wrt temperature*/
  dbl dCpdX[DIM]; /* heat capacity derivative wrt displacement */
  dbl dCpdC[MAX_CONC]; /* heat capacity derivative wrt concentration*/

  dbl X[DIM], T, C[MAX_CONC]; /* Convenient local variables */

  int i;

  /* Begin Execution */
 /**********************************************************/

 /* Comment out our remove this line if using this routine */

 /* EH(-1,"No user_defined heat capacity  model implemented");*/
 /**********************************************************/


 /************Initialize everything for saftey**************/
  Cp = 0;                                   /*Do not touch */
  dCpdT = 0;				    /*Do not touch */
  for ( a=0; a<DIM; a++)		    /*Do not touch */
    {					    /*Do not touch */
      dCpdX[a] = 0.;			    /*Do not touch */
    }					    /*Do not touch */
  for (i=0; i<MAX_CONC; i++)		    /*Do not touch */
    {					    /*Do not touch */
      dCpdC[i] = 0.;			    /*Do not touch */
    }					    /*Do not touch */
 /**********************************************************/

 /***********Load up convenient local variables*************/
 /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T;                                            /*Do not touch */
  for(a=0; a<DIM; a++)X[a] = fv->x[a];		        /*Do not touch */
  for(i=0; i<pd->Num_Species_Eqn; i++) C[i] = fv->c[i]; /*Do not touch */

 /**********************************************************/

 /*******Add property function and sensitivities here*******/


 /****************Don't touch these lines***********************/  
  /* NB nothing can depend on gradients of variables */      
  mp->heat_capacity = Cp;                                     /*Do not touch */
  mp->d_heat_capacity[TEMPERATURE] = dCpdT;		      /*Do not touch */
  for ( a=0; a<DIM; a++)				      /*Do not touch */
    {							      /*Do not touch */
      mp->d_heat_capacity[MESH_DISPLACEMENT1+a] = dCpdX[a];   /*Do not touch */
    }							      /*Do not touch */
  for ( a=0; a<MAX_CONC; a++)				      /*Do not touch */
    {							      /*Do not touch */
      mp->d_heat_capacity[MAX_VARIABLE_TYPES+a] = dCpdC[a];   /*Do not touch */
    }							      /*Do not touch */
 /**********************************************************/

  return(0);
} /* End of usr_heat_capacity */
/*****************************************************************************/

/*
 * HEAT SOURCE 
 */
/*
 * int usr_heat_source ()
 *
 * ------------------------------------------------------------------------------
 * This routine is responsible for filling up the following portions of the
 * mp structure at the current gauss point:
 *     intput:    param - array of constants input on the property card.  
 *
 *     output:  h       => mp->heat_source - thermal conductivity
 *              dhdT    => mp->d_heat_source[TEMPERATURE] 
 *                                         - derivative wrt temperature.
 *              dhdC    =>  mp->d_heat_source[MASS_FRACTION][i]
 *                                         - derivative wrt mass frac species i
 *              dhdV[0] => mp->d_heat_source[VELOCITY1]
 *                         mp->d_heat_source[VELOCITY2]
 *                         mp->d_heat_source[VELOCITY3]
 *                                         - derivative wrt velocities
 *              dhdX[0] => mp->d_heat_source[MESH_DISPLACEMENT1]
 *                         mp->d_heat_source[MESH_DISPLACEMENT2]
 *                         mp->d_heat_source[MESH_DISPLACEMENT3]
 *                                         - derivative wrt mesh displacements
 *
 *   NB: The user need only supply Cp, dCpdT, dCpdC, etc....mp struct is loaded up for you
 * ----------------------------------------------------------------------------
 * Example of an everything-dependent heat source for material "sample"
 * -----------Add these lines just before last section.-------------
 *
 * Simple heat source function for one material 
 *
 *  if (!strcmp(pd->MaterialName, "sample") )   
 *    {
 *	h = param[0]*T + param[1]*X[0] + param[2]*X[1] + param[3]*C[0] + param[4];
 *     dhdT    = param[0];
 *     dhdX[0] = param[1];
 *     dhdX[1] = param[2];
 *     dhdC[0] = param[3];
 *    } 
 *  else
 *    {
 *	EH(-1,"No user-defined function for this material");
 *    }
 * ----------------------------------------------------------------------------
 */
int
usr_heat_source(dbl *param, dbl time)	/* ptr to the user-defined parameter list */
{
  int a;

  dbl h, dhdT, dhdV;   /*heat sourceand its derivative wrt temperature, voltage*/
  dbl dhdv[DIM]; /* heat source derivative wrt velocity*/
  dbl dhdC[MAX_CONC]; /* heat source derivative wrt concentration*/
  dbl dhdX[DIM]; /* heat source derivative wrt displacement*/

  dbl X[DIM], T, V, C[MAX_CONC]; /* Convenient local variables */

  int i;
  int model_id;

  dbl intensity, intensity_dt=0;
  /* Begin Execution */

 /**********************************************************/

 /* Comment out our remove this line if using this routine 

    EH(-1,"No user_defined heat source model implemented");*/
 /**********************************************************/


 /************Initialize everything for saftey**************/
  h = 0;                                 /*Do not touch */
  dhdT = 0;				 /*Do not touch */
  dhdV = 0;				 /*Do not touch */
  for(a=0; a<DIM; a++) dhdv[a]=0.; 	 /*Do not touch */
  for(a=0; a<DIM; a++) dhdX[a]=0.;       /*Do not touch */
  for(a=0; a<MAX_CONC; a++) dhdC[a]=0.;  /*Do not touch */
 /**********************************************************/

 /***********Load up convenient local variables*************/
 /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T;                                             /*Do not touch */
  V = fv->V;                                             /*Do not touch */
  for(a=0; a<DIM; a++)X[a] = fv->x[a];		         /*Do not touch */
/*  for(i=0; i<pd->Num_Species_Eqn; i++) C[i] = fabs(fv->c[i]);  Do not touch */
  for(i=0; i<pd->Num_Species_Eqn; i++) C[i] = fv->c[i];  /*Do not touch */

 /**********************************************************/

 /*******Add property function and sensitivities here*******/

 model_id = ((int)param[0]);

 if (model_id == 1)
	{

/* Photo-polymerization reaction heat source            */
/*   if( efv->Num_external_field > 1) intensity_dt = fv->external_field[1];
 *          intensity = fv->external_field[0] + 
 *                                  intensity_dt*(tran->time_value - tran->init_time);
 *                                  */
      intensity = param[5]*(SQUARE(fv->apr)+SQUARE(fv->api));

/**  heat source from momomer heat of reaction     **/

 h = (param[1]*exp(-param[2]*(1./T - 1./param[3])))
                *C[2]*sqrt( intensity*C[0] );
 dhdC[0] = (param[1]*exp(-param[2]*(1./T - 1./param[3])))
                *C[2]*0.5*sqrt( intensity/C[0] );
 dhdC[2] = (param[1]*exp(-param[2]*(1./T - 1./param[3])))
                *sqrt( intensity*C[0] );
 dhdT = (param[1]*(param[2]/(T*T))*exp(-param[2]*(1./T - 1./param[3])))
                *C[2]*sqrt( intensity*C[0] );

/**  add heat generation from light absorption  **/

 h +=  param[4]*intensity*C[0];
 dhdC[0] += param[4]*intensity;
        }
 else
 if (model_id == 2)
       {

/*   Microwave Drying Heat Source		*/
	double amplitude, beta, y_surf;
	amplitude = param[1];
	beta = param[2];
	y_surf = param[3];

	h = amplitude*exp(-beta*(y_surf - X[1]));
	dhdX[1] = amplitude*beta*exp(-beta*(y_surf-X[1]));

       }
 else
 if (model_id == 3)
       {

/*   Microwave Drying Heat Source		*/
	double prefactor, k, R, alpha;
	prefactor = param[1];

	R = acoustic_impedance( NULL, time );
  	k = wave_number( NULL, time );
  	alpha = acoustic_absorption( NULL, time );

	h = prefactor*(alpha*k/R)*(SQUARE(fv->apr)+SQUARE(fv->api));
	mp->d_heat_source[ACOUS_PREAL] = 
		prefactor*(alpha*k/R)*2.*fv->apr;
	mp->d_heat_source[ACOUS_PIMAG] = 
		prefactor*(alpha*k/R)*2.*fv->api;
#if 0
        if(mp->wave_numberModel == TABLE)
           {
            struct  Data_Table *table_local;
            int var;
            table_local = MP_Tables[mp->wave_number_tableid];
            apply_table_mp( &mp->wave_number, table_local );
            for(i=0;i<table_local->columns-1;i++)
              {
              var = table_local->t_index[i];
              /* currently only set up to vary w.r.t. temperature */
              switch (var)
                {
                case TEMPERATURE:
                      dhdT += prefactor*(alpha/R)*table_local->slope[i]*
                               (SQUARE(fv->apr)+SQUARE(fv->api));
                  break;
                default:
                      EH(-1, "Variable function not yet implemented in material property table");
                }
               }
           }
        if(mp->Acoustic_ImpedanceModel == TABLE)
           {
            struct  Data_Table *table_local;
            int var;
            table_local = MP_Tables[mp->acoustic_impedance_tableid];
            apply_table_mp( &mp->acoustic_impedance, table_local );
            for(i=0;i<table_local->columns-1;i++)
              {
              var = table_local->t_index[i];
              /* currently only set up to vary w.r.t. temperature */
              switch (var)
                {
                case TEMPERATURE:
                      dhdT += prefactor*(-alpha*k*SQUARE(R))*table_local->slope[i]*
                               (SQUARE(fv->apr)+SQUARE(fv->api));
                  break;
                default:
                      EH(-1, "Variable function not yet implemented in material property table");
                }
               }
           }
        if(mp->Acoustic_AbsorptionModel == TABLE)
           {
            struct  Data_Table *table_local;
            int var;
            table_local = MP_Tables[mp->acoustic_absorption_tableid];
            apply_table_mp( &mp->acoustic_absorption, table_local );
            for(i=0;i<table_local->columns-1;i++)
              {
              var = table_local->t_index[i];
              /* currently only set up to vary w.r.t. temperature */
              switch (var)
                {
                case TEMPERATURE:
                      dhdT += prefactor*(k/R)*table_local->slope[i]*
                               (SQUARE(fv->apr)+SQUARE(fv->api));
                  break;
                default:
                      EH(-1, "Variable function not yet implemented in material property table");
                }
               }
           }
#endif
       }
 else 
 if (model_id == 4)
        {
  double k_prop, k_term, k_inh;

/* Photo-polymerization reaction heat source		*/
/*   if( efv->Num_external_field > 1) intensity_dt = fv->external_field[1];
       intensity = fv->external_field[0] + 
			intensity_dt*(tran->time_value - tran->init_time);
      intensity = param[5]*(SQUARE(fv->apr)+SQUARE(fv->api));
*/
      intensity = param[1]*(fv->poynt[1]);

	    k_prop = param[3]*exp(-param[4]*(1./T - 1./param[5]));
	    k_term = param[6]*exp(-param[7]*(1./T - 1./param[8]));
	    k_inh = param[9]*exp(-param[10]*(1./T - 1./param[11]));

/**  heat source from momomer heat of reaction     **/

 h = k_prop*C[2]*C[4] * param[12];
 dhdC[2] = k_prop*C[4] * param[12];
 dhdC[4] = k_prop*C[2] * param[12];
 dhdT = k_prop*param[4]/SQUARE(T)*C[2]*C[4]*param[12];

/**  add heat generation from light absorption  **/

 h +=  param[13]*intensity*C[0];
 dhdC[0] += param[13]*intensity;
        }
 else
       {
    EH(-1,"Invalid user heat source model number");
       }


 /****************Don't touch these lines***********************/  
  mp->heat_source = h;                                   /*Do not touch */
  mp->d_heat_source[TEMPERATURE] = dhdT;		 /*Do not touch */
  mp->d_heat_source[VOLTAGE] = dhdV;	        	 /*Do not touch */
							 /*Do not touch */
  for ( a=0; a<DIM; a++)                                 /*Do not touch */
    {							 /*Do not touch */
      mp->d_heat_source[VELOCITY1+a] = dhdv[a];		 /*Do not touch */
    }                                                    /*Do not touch */
  for ( a=0; a<DIM; a++)                                 /*Do not touch */
    {             					 /*Do not touch */
      mp->d_heat_source[MESH_DISPLACEMENT1+a] = dhdX[a]; /*Do not touch */
    }							 /*Do not touch */
  for ( a=0; a<MAX_CONC; a++)				 /*Do not touch */
    {                                                    /*Do not touch */
      mp->d_heat_source[MAX_VARIABLE_TYPES+a] = dhdC[a]; /*Do not touch */
    }							 /*Do not touch */

  return(0);
} /* End of usr_heat_source */
/*****************************************************************************/

/*
 * SPECIES SOURCE 
 */
/*
 * int usr_species_source ()
 *
 * ----------------------------------------------------------------------------
 * This routine is responsible for filling up the following portions of the mp structure
 * at the current gauss point:
 *     intput:    param - array of constants input on the property card.  
 *
 *     output:  s       => mp->species_source - thermal conductivity
 *              dsdT    => mp->d_species_source[TEMPERATURE] 
 *                                         - derivative wrt temperature.
 *              dsdC    =>  mp->d_species_source[MASS_FRACTION][i]
 *                                         - derivative wrt mass frac species i
 *              dsdV[0] => mp->d_species_source[VELOCITY1]
 *                         mp->d_species_source[VELOCITY2]
 *                         mp->d_species_source[VELOCITY3]
 *                                         - derivative wrt velocities
 *              dsdX[0] => mp->d_species_source[MESH_DISPLACEMENT1]
 *                         mp->d_species_source[MESH_DISPLACEMENT2]
 *                         mp->d_species_source[MESH_DISPLACEMENT3]
 *                                         - derivative wrt mesh displacements
 *
 *   NB: The user need only supply s, dsdT, dsdC, etc....mp struct is loaded up for you
 * ----------------------------------------------------------------------------
 * EXAMPLE of an everything-dependent spies source for material "sample"
 * -----------Add these lines just before last section.-------------
 *
 * Simple species source function for one material 
 *
 *  if (!strcmp(pd->MaterialName, "sample") )   
 *    {
 *      if(species_no == 0) 
 *        {
 *           s = param[0]*T + param[1]*X[0] + param[2]*C[0];
 *           dsdT = param[0];
 *           dsdX[0] = param[1];
 *           dsdC[0] = param[2];
 *        }
 *      if(species_no == 1)
 *        {
 *           s = param[2]*C[0]*C[0];
 *           dsdT = 0;
 *           dsdX[0] =0;
 *           dsdC[0] = 2.*param[2]*C[0];
 *         }
 *    } 
 *  else
 *    {
 *	EH(-1,"No user-defined function for this material");
 *    }                                                        
 * ----------------------------------------------------------------------------
 */

int 
usr_species_source(int species_no, /* Current species number                 */
		   dbl *param)	/* pointer to user-defined parameter list    */

{
  int a;

  dbl s, dsdT, dsdV;   /*species sourceand its derivative wrt temperature, voltage*/
  dbl dsdv[DIM]; /* species source derivative wrt velocity*/
  dbl dsdC[MAX_CONC]; /* species source derivative wrt concentration*/
  dbl dsdX[DIM]; /* species source derivative wrt displacement*/

  dbl X[DIM], T, V, C[MAX_CONC]; /* Convenient local variables */

  int i;

  dbl intensity, intensity_dt=0;
  double k_prop, k_term, k_inh;
  /* Begin Execution */

 /**********************************************************/

 /* Comment out our remove this line if using this routine */

 /*  EH(-1,"No user_defined species source model implemented");
 *********************************************************/


 /************Initialize everything for saftey**************/
  s = 0;                                 /*Do not touch */
  dsdT = 0;				 /*Do not touch */
  dsdV = 0;				 /*Do not touch */
  for(a=0; a<DIM; a++) dsdv[a]=0.; 	 /*Do not touch */
  for(a=0; a<DIM; a++) dsdX[a]=0.;       /*Do not touch */
  for(a=0; a<MAX_CONC; a++) dsdC[a]=0.;  /*Do not touch */
 /**********************************************************/

 /***********Load up convenient local variables*************/
 /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T;                                              /*Do not touch */
  V = fv->V;                                              /*Do not touch */
  for(a=0; a<DIM; a++)X[a] = fv->x[a];	                  /*Do not touch */
  for(i=0; i<pd->Num_Species_Eqn; i++) C[i] = fv->c[i];   /*Do not touch */

 /**********************************************************/

 /*******Add property function and sensitivities here*******/

 /**  Photolysis Reaction for Photo-modeling		**/
 /**  C[0] = Photo-initiator
      C[1] = Photolysis product
      C[2] = Monomer
      C[3] = Polymerization product?
  **/
/*   if( efv->Num_external_field > 1) intensity_dt = fv->external_field[1];
       intensity = fv->external_field[0] + 
			intensity_dt*(tran->time_value - tran->init_time);
       intensity = MAX(intensity,0);
*/
/*fprintf(stderr,"intensity %g %g %g\n",fv->x[0],fv->x[1],intensity);*/
/*      intensity = param[3]*(SQUARE(fv->apr)+SQUARE(fv->api));  */
#if 1
      intensity = param[3]*(SQUARE(fv->apr)+SQUARE(fv->api));

       if(species_no == 0)
         {
            s = -param[0]*C[0]*intensity;
            dsdC[0] = -param[0]*intensity;
         }
       else if(species_no == 1)
         {
            s = param[0]*C[0]*intensity;
            dsdC[0] = param[0]*intensity;
          }
       else if(species_no == 2)
         {
            s = -(param[0]*exp(-param[1]*(1./T - 1./param[2])))
                        *C[2]*sqrt( intensity*C[0] );
            dsdC[0] = -(param[0]*exp(-param[1]*(1./T - 1./param[2])))
                        *C[2]*0.5*sqrt( intensity/C[0] );
            dsdC[2] = -(param[0]*exp(-param[1]*(1./T - 1./param[2])))
                        *sqrt( intensity*C[0] );
            dsdT = -(param[0]*(param[1]/(T*T))*exp(-param[1]*(1./T-1./param[2])))
                        *C[2]*sqrt( intensity*C[0] );
          }
       else
          {
                EH(-1,"too many species for usr_species_source");
          }
#else
      intensity = param[0]*(fv->poynt[1]);

	    k_prop = param[2]*exp(-param[3]*(1./T - 1./param[4]));
	    k_term = param[5]*exp(-param[6]*(1./T - 1./param[7]));
	    k_inh = param[8]*exp(-param[9]*(1./T - 1./param[10]));
       if(species_no == 0) 
         {
            s = -C[0]*intensity;
            dsdC[0] = -intensity;
         }
       else if(species_no == 1)
         {
            s = param[1]*C[0]*intensity;
            dsdC[0] = param[1]*intensity;
          }
       else if(species_no == 2)
         {
            s = -k_prop*C[2]*C[4];
            dsdC[2] = -k_prop*C[4];
            dsdC[4] = -k_prop*C[2];
            dsdT = -C[2]*C[4]*k_prop*param[3]/SQUARE(T);
          }
       else if(species_no == 3)
         {
            s = -k_inh*C[3]*C[4];
            dsdC[3] = -k_inh*C[4];
            dsdC[4] = -k_inh*C[3];
            dsdT = -C[3]*C[4]*k_inh*param[9]/SQUARE(T);
          }
       else if(species_no == 4)
         {
		 s = param[1]*C[0]*intensity-k_term*SQUARE(C[4])-k_inh*C[3]*C[4];
                 dsdC[0] = param[1]*intensity;
                 dsdC[3] = -k_inh*C[4];
                 dsdC[4] = -k_term*2*C[4] - k_inh*C[3];
                 dsdT = - SQUARE(C[4])*k_term*param[6]/SQUARE(T)
                           -C[3]*C[4]*k_inh*param[9]/SQUARE(T);
          }
       else
	  {
		EH(-1,"too many species for usr_species_source");
	  }
#endif

 /****************Don't touch these lines***********************/  
  mp->species_source[species_no] = s;                    /*Do not touch */
  mp->d_species_source[TEMPERATURE] = dsdT;		 /*Do not touch */
  mp->d_species_source[VOLTAGE] = dsdV;	        	 /*Do not touch */
							 /*Do not touch */
  for ( a=0; a<DIM; a++)                                 /*Do not touch */
    {							 /*Do not touch */
      mp->d_species_source[VELOCITY1+a] = dsdv[a];	 /*Do not touch */
    }                                                    /*Do not touch */
  for ( a=0; a<DIM; a++)                                 /*Do not touch */
    {             					 /*Do not touch */
      mp->d_species_source[MESH_DISPLACEMENT1+a] = dsdX[a]; /*Do not touch */
    }							 /*Do not touch */
  for ( a=0; a<MAX_CONC; a++)				 /*Do not touch */
    {                                                    /*Do not touch */
      mp->d_species_source[MAX_VARIABLE_TYPES+a] = dsdC[a]; /*Do not touch */
    }							 /*Do not touch */

  return(0);
} /* End of usr_species_source */
/*****************************************************************************/

/*
 * CURRENT SOURCE 
 */
/*
 * int usr_current_source ()
 *
 * ----------------------------------------------------------------------------
 * This routine is responsible for filling up the following portions of the
 * mp structure at the current gauss point:
 *     intput:    param - array of constants input on the property card.  
 *
 *     output:  h       => mp->current_source - current source
 *              dhdT    => mp->d_current_source[TEMPERATURE] 
 *                                         - derivative wrt temperature.
 *              dhdV    => mp->d_current_source[VOLTAGE] 
 *                                         - derivative wrt voltage.
 *              dhdC    =>  mp->d_current_source[MASS_FRACTION][i]
 *                                         - derivative wrt mass frac species i
 *              dhdv[0] => mp->d_current_source[VELOCITY1]
 *                         mp->d_current_source[VELOCITY2]
 *                         mp->d_current_source[VELOCITY3]
 *                                         - derivative wrt velocities
 *              dhdX[0] => mp->d_current_source[MESH_DISPLACEMENT1]
 *                         mp->d_current_source[MESH_DISPLACEMENT2]
 *                         mp->d_current_source[MESH_DISPLACEMENT3]
 *                                         - derivative wrt mesh displacements
 *
 *   NB: The user need only supply h, dhdT, dhdC, etc....mp struct is loaded up for you
 * ----------------------------------------------------------------------------
 * Example of an everything-dependent heat source for material "sample"
 * -----------Add these lines just before last section.-------------
 *
 * Simple heat source function for one material 
 *
 *  if (!strcmp(pd->MaterialName, "sample") )   
 *    {
 *	h = param[0]*T + param[1]*X[0] + param[2]*X[1] + param[3]*C[0] + param[4];
 *     dhdT    = param[0];
 *     dhdX[0] = param[1];
 *     dhdX[1] = param[2];
 *     dhdC[0] = param[3];
 *    } 
 *  else
 *    {
 *	EH(-1,"No user-defined function for this material");
 *    }
 * ----------------------------------------------------------------------------
 */

int 
usr_current_source(dbl *param)	/* pointer to user-defined parameter list */
{
  int a;

  dbl h, dhdT;   /* current source and its derivative wrt temperature*/
  dbl dhdV;      /* current source derivative wrt voltage */
  dbl dhdv[DIM]; /* current source derivative wrt velocity*/
  dbl dhdC[MAX_CONC]; /* current source derivative wrt concentration*/
  dbl dhdX[DIM]; /* current source derivative wrt displacement*/

  dbl X[DIM], T, V, C[MAX_CONC]; /* Convenient local variables */

  int i;

  /* Begin Execution */

 /**********************************************************/

 /* Comment out our remove this line if using this routine */

    EH(-1,"No user_defined current source model implemented");  
 /**********************************************************/


 /************Initialize everything for saftey**************/
  h = 0;                                 /*Do not touch */
  dhdT = 0;				 /*Do not touch */
  dhdV = 0;				 /*Do not touch */
  for(a=0; a<DIM; a++) dhdv[a]=0.; 	 /*Do not touch */
  for(a=0; a<DIM; a++) dhdX[a]=0.;       /*Do not touch */
  for(a=0; a<MAX_CONC; a++) dhdC[a]=0.;  /*Do not touch */
 /**********************************************************/

 /***********Load up convenient local variables*************/
 /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T;                                             /*Do not touch */
  V = fv->V;                                             /*Do not touch */
  for(a=0; a<DIM; a++)X[a] = fv->x[a];		         /*Do not touch */
  for(i=0; i<pd->Num_Species_Eqn; i++) C[i] = fv->c[i];  /*Do not touch */

 /**********************************************************/

 /*******Add property function and sensitivities here*******/


 /****************Don't touch these lines***********************/  
  mp->current_source = h;                                   /*Do not touch */
  mp->d_current_source[TEMPERATURE] = dhdT;		    /*Do not touch */
  mp->d_current_source[VOLTAGE] = dhdV;		            /*Do not touch */
							    /*Do not touch */
  for ( a=0; a<DIM; a++)                                    /*Do not touch */
    {							    /*Do not touch */
      mp->d_current_source[VELOCITY1+a] = dhdv[a];	    /*Do not touch */
    }                                                       /*Do not touch */
  for ( a=0; a<DIM; a++)                                    /*Do not touch */
    {             					    /*Do not touch */
      mp->d_current_source[MESH_DISPLACEMENT1+a] = dhdX[a]; /*Do not touch */
    }							    /*Do not touch */
  for ( a=0; a<MAX_CONC; a++)				    /*Do not touch */
    {                                                       /*Do not touch */
      mp->d_current_source[MAX_VARIABLE_TYPES+a] = dhdC[a]; /*Do not touch */
    }							    /*Do not touch */

  return(0);
} /* End of usr_current_source */
/*****************************************************************************/

/*
 * VISCOSITY 
 */
/*
 * int usr_viscosity ()
 *
 * ----------------------------------------------------------------------------
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
 * ----------------------------------------------------------------------------
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
 *	EH(-1,"No user-defined function for this material");
 *    }
 * ----------------------------------------------------------------------------
 */

int 
usr_viscosity(dbl *param)	/* pointer to user-defined parameter list    */
{
  /* Local Variables */
  int a;
  dbl mu, dmudT;   /* thermal conductivity and its derivative wrt temperature*/
  dbl dmudV[DIM]; /* heat source derivative wrt velocity*/
  dbl dmudC[MAX_CONC]; /* heat source derivative wrt concentration*/
  dbl dmudX[DIM]; /* heat source derivative wrt displacement*/
  dbl X[DIM], F, T, C[MAX_CONC]; /* Convenient local variables */

  /* Begin Execution */

 /**********************************************************/

 /* Comment out our remove this line if using this routine */

 EH(-1,"No user_defined viscosity model implemented");

 /**********************************************************/


 /************Initialize everything for saftey**************/
  mu = 0;                                  /*Do not touch */
  dmudT = 0;				  /*Do not touch */
  for(a=0; a<DIM; a++) dmudV[a]=0.; 	  /*Do not touch */
  for(a=0; a<DIM; a++) dmudX[a]=0.;       /*Do not touch */
  for(a=0; a<MAX_CONC; a++) dmudC[a]=0.;  /*Do not touch */
 /**********************************************************/

 /***********Load up convenient local variables*************/
 /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T;                                        /*Do not touch */
  F = fv->F;                                        /*Do not touch */
  for(a=0; a<DIM; a++)X[a] = fv->x[a];		    /*Do not touch */
  for(a=0; a<pd->Num_Species_Eqn; a++) C[a] = fv->c[a];    /*Do not touch */

 /**********************************************************/

 /*******Add property function and sensitivities here*******/
 /* An example:
  *  mu = param[0] * exp(param[1] * (param[2] - C[0]));
  *
  * Add sensitivities here
  *  dmudC[0] = - mu * param[1];
  *
  */

  mu = param[0];
  if(fv->c[0] < 0.5)
    {
      mu=param[1];
    } 

 /****************Don't touch these lines***********************/  
  mp->viscosity = mu;                                    /*Do not touch */
  mp->d_viscosity[TEMPERATURE] = dmudT;		         /*Do not touch */
							 /*Do not touch */
  for ( a=0; a<DIM; a++)                                 /*Do not touch */
    {							 /*Do not touch */
      mp->d_viscosity[VELOCITY1+a] = dmudV[a];		 /*Do not touch */
    }                                                    /*Do not touch */
  for ( a=0; a<DIM; a++)                                 /*Do not touch */
    {             					 /*Do not touch */
      mp->d_viscosity[MESH_DISPLACEMENT1+a] = dmudX[a];  /*Do not touch */
    }							 /*Do not touch */
  for ( a=0; a<MAX_CONC; a++)				 /*Do not touch */
    {                                                    /*Do not touch */
      mp->d_viscosity[MAX_VARIABLE_TYPES+a] = dmudC[a];  /*Do not touch */
    }							 /*Do not touch */

  return(0);
} /* End of usr_viscosity */
/*****************************************************************************/

int 
usr_surface_tension(dbl *param)	/* ptr to user-defined parameter list        */
{
  int a;

  dbl sigma, dsigmadT;   /* thermal conductivity and its derivative wrt temperature*/
  dbl dsigmadV[DIM]; /* heat source derivative wrt velocity*/
  dbl dsigmadC[MAX_CONC]; /* heat source derivative wrt concentration*/
  dbl dsigmadX[DIM]; /* heat source derivative wrt displacement*/

  dbl X[DIM], T, C[MAX_CONC]; /* Convenient local variables */

  int i;

  /* Begin Execution */

 /**********************************************************/

 /* Comment out our remove this line if using this routine */

 /*  EH(-1,"No user_defined surface_tension model implemented");  */

 /**********************************************************/


 /************Initialize everything for saftey**************/
  sigma = 0;                                  /*Do not touch */
  dsigmadT = 0;				  /*Do not touch */
  for(a=0; a<DIM; a++) dsigmadV[a]=0.; 	  /*Do not touch */
  for(a=0; a<DIM; a++) dsigmadX[a]=0.;       /*Do not touch */
  for(a=0; a<MAX_CONC; a++) dsigmadC[a]=0.;  /*Do not touch */
 /**********************************************************/

 /***********Load up convenient local variables*************/
 /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T;                                        /*Do not touch */
  for(a=0; a<DIM; a++)X[a] = fv->x[a];		    /*Do not touch */
  for(i=0; i<pd->Num_Species_Eqn; i++) C[i] = fv->c[i];    /*Do not touch */

 /*******Add property function and sensitivities here*******/
 /* An example:
  *  sigma = param[0] + (param[1] - param[0]) * C[0];
  *
  * Add sensitivities here
  *
  *  dsigmadC[0] =  param[1] - param[0];
  *
  */

/*  surface tension as a fcn of temperature, concentration */


	sigma = param[0] + param[1]*(T - mp->reference[TEMPERATURE]);
	dsigmadT = param[1];
	for(i=2 ; i<mp->len_u_surface_tension ; i++)
		{
		sigma += param[i]*(C[i-2]-mp->reference_concn[i-2]);
		dsigmadC[i-2] = param[i];
		}

 /****************Don't touch these lines***********************/  
  mp->surface_tension = sigma;		         /*Do not touch */
  mp->d_surface_tension[TEMPERATURE] = dsigmadT;	 /*Do not touch */
							 /*Do not touch */
  for ( a=0; a<DIM; a++)                                 /*Do not touch */
    {							 /*Do not touch */
      mp->d_surface_tension[VELOCITY1+a] = dsigmadV[a];	 /*Do not touch */
    }                                                    /*Do not touch */
  for ( a=0; a<DIM; a++)                                 /*Do not touch */
    {             					 /*Do not touch */
      mp->d_surface_tension[MESH_DISPLACEMENT1+a] = dsigmadX[a]; /*Do not touch */
    }							 /*Do not touch */
  for ( a=0; a<MAX_CONC; a++)				 /*Do not touch */
    {                                                    /*Do not touch */
      mp->d_surface_tension[MAX_VARIABLE_TYPES+a] = dsigmadC[a]; /*Do not touch */
    }							 /*Do not touch */

  return(0);
} /* End of usr_surface_tension */
/*****************************************************************************/

/*
 *  MOMENTUM SOURCE 
 */
/*
 * int usr_momentum_source ()
 *
 * ----------------------------------------------------------------------------
 * This routine is responsible for filling up the following portions of the mp structure
 * at the current gauss point:
 *     intput:  param - array of constants input on the property card.  
 *
 *     output:  f[a]       => mp->momentum_source[a] - body force in direction [a]
 *              dfdT[a]    => mp->d_momentum_source[a][TEMPERATURE] 
 *                                         - derivative wrt temperature.
 *              dfdC[a][i] => mp->d_momentum_source[a][MAX_VARIABLE_TYPES+i]
 *                                         - derivative wrt mass frac species i
 *              dfdV[a][0] => mp->d_momentum_source[a][VELOCITY1]
 *                            mp->d_momentum_source[a][VELOCITY2]
 *                            mp->d_momentum_source[VELOCITY3]
 *                                         - derivative wrt velocities
 *              dfdX[a][0] => mp->d_momentum_source[a][MESH_DISPLACEMENT1]
 *                            mp->d_momentum_source[a][MESH_DISPLACEMENT2]
 *                            mp->d_momentum_source[a][MESH_DISPLACEMENT3]
 *                                         - derivative wrt mesh displacements
 *
 *   NB: The user need only supply f, dfdT, dfdC, etc....mp struct is loaded up for you
 * ----------------------------------------------------------------------------
 *
 * NB: This is the *entire* momentum source term -- if you want rho*g, then
 *     it's your responsibility to multiply by rho! If you only want JxB,
 *     then do not include the density! Either way, goma does not presuppose
 *     to multiply your source term by density.
 *
 *
 * Example of an everything-dependent momentum source in [0] direction for material "sample"
 * -----------Add these lines just before last section.-------------
 *
 * Simple heat source function for one material 
 *
 *  if (!strcmp(pd->MaterialName, "sample") )   
 *    {
 *     f[0] = param[0]*T + param[1]*X[0] + param[2]*X[1] + param[3]*C[0] + param[4];
 *     dfdT[0]    = param[0];
 *     dfdX[0][0] = param[1];
 *     dfdX[0][1] = param[2];
 *     dfdC[0][0] = param[3];
 *    } 
 *  else
 *    {
 *	EH(-1,"No user-defined function for this material");
 *    }
 * ----------------------------------------------------------------------------
 */

int
usr_momentum_source(dbl *param)	/* ptr to user-defined parameter list        */
{
  int a, b;
  int i;
  int w;

  dbl f[DIM], dfdT[DIM];   /* momentum sources and its derivative wrt temperature*/
  dbl dfdV[DIM][DIM];      /* momentum source derivative wrt velocity*/
  dbl dfdC[DIM][MAX_CONC]; /* momentum source derivative wrt concentration*/
  dbl dfdX[DIM][DIM];      /* momentum source derivative wrt displacement*/

  dbl X[DIM], T, C[MAX_CONC]; /* Convenient local variables */

  /* Begin Execution */

 /**********************************************************/

 /* Comment out our remove this line if using this routine */

    EH(-1,"No user_momentum_source model implemented");  
 /**********************************************************/


 /************Initialize everything for saftey**************/
  for(a=0; a<DIM; a++)
    {
      f[a] = 0;                                   /*Do not touch */
      dfdT[a] = 0;				  /*Do not touch */
      for(b=0; b<DIM; b++) dfdV[a][b]=0.; 	  /*Do not touch */
      for(b=0; b<DIM; b++) dfdX[a][b]=0.;         /*Do not touch */
      for(w=0; w<MAX_CONC; w++) dfdC[a][w]=0.;    /*Do not touch */
    }
 /**********************************************************/

 /***********Load up convenient local variables*************/
 /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T;                                             /*Do not touch */
  for(a=0; a<DIM; a++)X[a] = fv->x[a];		         /*Do not touch */
  for(i=0; i<pd->Num_Species_Eqn; i++) C[i] = fv->c[i];  /*Do not touch */

 /**********************************************************/

 /*******Add property function and sensitivities here*******/


 /****************Don't touch these lines***********************/ 
  for(a=0; a<DIM; a++)
    { 
      mp->momentum_source[a] = f[a];                          /*Do not touch */
      mp->d_momentum_source[a][TEMPERATURE] = dfdT[a];	      /*Do not touch */
                                                              /*Do not touch */
      for ( b=0; b<DIM; b++)                                  /*Do not touch */
	{						      /*Do not touch */
	  mp->d_momentum_source[a][VELOCITY1+b] = dfdV[a][b]; /*Do not touch */
	}                                                     /*Do not touch */
      for ( b=0; b<DIM; b++)                                  /*Do not touch */
	{             					      /*Do not touch */
	  mp->d_momentum_source[a][MESH_DISPLACEMENT1+a] = dfdX[a][b];    /*Do not touch */
	}						      /*Do not touch */
      for ( w=0; w<MAX_CONC; w++)			      /*Do not touch */
	{                                                     /*Do not touch */
	  mp->d_momentum_source[a][MAX_VARIABLE_TYPES+w] = dfdC[a][w];    /*Do not touch */
	}		
    }					                      /*Do not touch */

  return(0);
} /* End of usr_momentum_source */
/*****************************************************************************/

/*
 *  LAME MU - the shear modulus in elastic constitutive equations
 */
/*
 * int usr_lame_mu ()
 *
 * ----------------------------------------------------------------------------
 * This routine is responsible for filling up the following portions of the elc structure
 * at the current gauss point:
 *     intput:  param - array of constants input on the property card.  
 *
 *     output:  f       => elc->lame_mu - shear modulus
 *              dfdT    => elc->d_mu[a][TEMPERATURE] 
 *                                         - derivative wrt teelcerature.
 *              dfdC[i] => elc->d_mu[MAX_VARIABLE_TYPES+i]
 *                                         - derivative wrt mass frac species i
 *              dfdV[0] => elc->d_mu[VELOCITY1]
 *                            elc->d_mu[VELOCITY2]
 *                            elc->d_mu[VELOCITY3]
 *                                         - derivative wrt velocities
 *              dfdX[0] => elc->d_mu[MESH_DISPLACEMENT1]
 *                            elc->d_mu[MESH_DISPLACEMENT2]
 *                            elc->d_mu[MESH_DISPLACEMENT3]
 *                                         - derivative wrt mesh displacements
 *
 *   NB: The user need only supply f, dfdT, dfdC, etc....elc struct is loaded up for you
 */

int
usr_lame_mu(struct Elastic_Constitutive *ep, dbl *param)		/* ptr to user-defined parameter list        */
{
  int a, b;
  int w;

  dbl f, dfdT;   /* momentum sources and its derivative wrt temperature*/
  dbl dfdV[DIM];      /* momentum source derivative wrt velocity*/
  dbl dfdC[MAX_CONC]; /* momentum source derivative wrt concentration*/
  dbl dfdX[DIM];      /* momentum source derivative wrt displacement*/

  dbl X[DIM], T, C[MAX_CONC]; /* Convenient local variables */
  /*  dbl dist;*/
  int i;
  double Gzero, Geq, tl, C1, C2, Tref, at, dat_dT;

  /* Begin Execution */

 /**********************************************************/

 /* Comment out our remove this line if using this routine */

  /*EH(-1,"No user_lame_mu model implemented.");   */

 /**********************************************************/


 /************Initialize everything for saftey**************/
  f = 0;                                   /*Do not touch */
  dfdT = 0;				  /*Do not touch */
  for(b=0; b<DIM; b++) dfdV[b]=0.; 	  /*Do not touch */
  for(b=0; b<DIM; b++) dfdX[b]=0.;         /*Do not touch */
  for(w=0; w<MAX_CONC; w++) dfdC[w]=0.;    /*Do not touch */
  /**********************************************************/
  
 /***********Load up convenient local variables*************/
 /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T;                                             /*Do not touch */
  for(a=0; a<DIM; a++)X[a] = fv->x[a];		         /*Do not touch */
  for(i=0; i<pd->Num_Species_Eqn; i++) C[i] = fv->c[i];  /*Do not touch */

 /*******Add property function and sensitivities here*******/
  /* Example:
   *
   * This function makes mu large (10^4) near a contact point and
   * 0.5 far away from the contact point.
   * here we use the initial coordinates, so this value has no sensitivity 
   * with respect to position 
   *
   *
   *  if (mp->porosity < param[1])
   *    {
   *      f = param[0] * pow((1. - mp->porosity)/(1. - param[1]), param[2]);
   *    }
   *  else
   *    {
   *      f = param[0];
   *    }
   *
   */

  /* Example:
   *
   *  dist = pow(( fv->x0[0] * fv->x0[0] 
   *     + (0.009398 - fv->x0[1]) * (0.009398 - fv->x0[1])), 0.5);
   *  f = param[0] + param[1] / (pow(dist,param[3]) + param[2]);
   */ 

  /****************Add sensitivities here********************/

  /*  if (mp->porosity < param[1])
   *    {
   *      elc->d_lame_mu[POROSITY] = - param[2] * param[0] * 
   *	pow((1. - mp->porosity), param[2] - 1)/
   *	pow((1. - param[1]), param[2]);
   *    }
   *  else
   *    {
   *      elc->d_lame_mu[POROSITY] = 0;
   *    }
   */

/*  exponential decay wrt WLF time-temperature superposition	*/
  
	Gzero = param[0];
	Geq = param[1];
	tl = param[2];
	C1 = param[3];
	C2 = param[4];
	Tref = param[5];

	at = exp(-C1*(T-Tref)/(C2+T-Tref));
	f = Gzero*exp(-tl/at) + Geq;
	dat_dT = at*(-C1*C2)/SQUARE(C2+T-Tref);
	dfdT = Gzero*exp(-tl/at)*tl*dat_dT/SQUARE(at);
  /****************Don't touch these lines***********************/ 
  for(a=0; a<DIM; a++)
    { 
      ep->lame_mu = f;                                      /*Do not touch */
      ep->d_lame_mu[TEMPERATURE] = dfdT;		     /*Do not touch */
                                                             /*Do not touch */
      for ( b=0; b<DIM; b++)                                 /*Do not touch */
	{						     /*Do not touch */
	  ep->d_lame_mu[VELOCITY1+b] = dfdV[b];	     /*Do not touch */
	}                                                    /*Do not touch */
      for ( b=0; b<DIM; b++)                                 /*Do not touch */
	{             					     /*Do not touch */
	  ep->d_lame_mu[MESH_DISPLACEMENT1+a] = dfdX[b];    /*Do not touch */
	}						     /*Do not touch */
      for ( w=0; w<MAX_CONC; w++)			     /*Do not touch */
	{                                                    /*Do not touch */
	  ep->d_lame_mu[MAX_VARIABLE_TYPES+w] = dfdC[w];    /*Do not touch */
	}			                             /*Do not touch */
    }					                     /*Do not touch */

  return(0);
} /* End of usr_lame_mu */
/*****************************************************************************/

/*
 *  LAME LAMBDA - the shear modulus in elastic constitutive equations
 */
/*
 * int usr_lame_lambda ()
 *
 * ----------------------------------------------------------------------------
 * This routine is responsible for filling up the following portions of the elc structure
 * at the current gauss point:
 *     intput:  param - array of constants input on the property card.  
 *
 *     output:  f       => elc->lame_lambda - shear modulus
 *              dfdT    => elc->d_lambda[a][TEMPERATURE] 
 *                                         - derivative wrt teelcerature.
 *              dfdC[i] => elc->d_lambda[MAX_VARIABLE_TYPES+i]
 *                                         - derivative wrt mass frac species i
 *              dfdV[0] => elc->d_lambda[VELOCITY1]
 *                            elc->d_lambda[VELOCITY2]
 *                            elc->d_lambda[VELOCITY3]
 *                                         - derivative wrt velocities
 *              dfdX[0] => elc->d_lambda[MESH_DISPLACEMENT1]
 *                            elc->d_lambda[MESH_DISPLACEMENT2]
 *                            elc->d_lambda[MESH_DISPLACEMENT3]
 *                                         - derivative wrt mesh displacements
 *
 *   NB: The user need only supply f, dfdT, dfdC, etc....elc struct is loaded up for you
 */

int
usr_lame_lambda(struct Elastic_Constitutive *ep, dbl *param)	/* ptr to user-defined parameter list        */
{
  int a, b;
  int w;

  dbl f, dfdT;   /* momentum sources and its derivative wrt temperature*/
  dbl dfdV[DIM];      /* momentum source derivative wrt velocity*/
  dbl dfdC[MAX_CONC]; /* momentum source derivative wrt concentration*/
  dbl dfdX[DIM];      /* momentum source derivative wrt displacement*/

  dbl X[DIM], T, C[MAX_CONC]; /* Convenient local variables */

  int i;
  double Gzero, Geq, tl, C1, C2, Tref, at, dat_dT;

  /* Begin Execution */

 /**********************************************************/

 /* Comment out our remove this line if using this routine */

 /* EH(-1,"No user_lame_lambda model implemented.");     */

 /**********************************************************/


 /************Initialize everything for saftey**************/
  f = 0;                                   /*Do not touch */
  dfdT = 0;			           /*Do not touch */
  for(b=0; b<DIM; b++) dfdV[b]=0.;         /*Do not touch */
  for(b=0; b<DIM; b++) dfdX[b]=0.;         /*Do not touch */
  for(w=0; w<MAX_CONC; w++) dfdC[w]=0.;    /*Do not touch */
  /**********************************************************/
  
 /***********Load up convenient local variables*************/
 /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T;                                             /*Do not touch */
  for(a=0; a<DIM; a++)X[a] = fv->x[a];		         /*Do not touch */
  for(i=0; i<pd->Num_Species_Eqn; i++) C[i] = fv->c[i];  /*Do not touch */

 /**********************************************************/

 /*******Add property function and sensitivities here*******/
  /* Example:
   *
   * if (mp->porosity < param[1])
   *    {
   *      f = param[0] * pow((1. - mp->porosity)/(1. - param[1]), param[2]);
   *    }
   *  else
   *    {
   *      f = param[0];
   *    }
   *
   * Add sensitivities here.
   * 
   *  if (mp->porosity < param[1])
   *    {
   *      elc->d_lame_lambda[POROSITY] = - param[2] * param[0] * 
   *	pow((1. - mp->porosity), param[2] - 1)/
   *	pow((1. - param[1]), param[2]);
   *    }
   *  else
   *    {
   *      elc->d_lame_lambda[POROSITY] = 0;
   *    }
   */

/*  exponential decay wrt WLF time-temperature superposition	*/
   
	Gzero = param[0];
 	Geq = param[1];
 	tl = param[2];
 	C1 = param[3];
 	C2 = param[4];
 	Tref = param[5];
 
 	at = exp(-C1*(T-Tref)/(C2+T-Tref));
 	f = Gzero*exp(-tl/at) + Geq;
 	dat_dT = at*(-C1*C2)/SQUARE(C2+T-Tref);
 	dfdT = Gzero*exp(-tl/at)*tl*dat_dT/SQUARE(at);
 /****************Don't touch these lines***********************/ 
  for(a=0; a<DIM; a++)
    { 
      ep->lame_lambda = f;                                   /*Do not touch */
      ep->d_lame_lambda[TEMPERATURE] = dfdT;		      /*Do not touch */
                                                              /*Do not touch */
      for ( b=0; b<DIM; b++)                                  /*Do not touch */
	{						      /*Do not touch */
	  ep->d_lame_lambda[VELOCITY1+b] = dfdV[b];	      /*Do not touch */
	}                                                     /*Do not touch */
      for ( b=0; b<DIM; b++)                                  /*Do not touch */
	{             					      /*Do not touch */
	  ep->d_lame_lambda[MESH_DISPLACEMENT1+a] = dfdX[b]; /*Do not touch */
	}						      /*Do not touch */
      for ( w=0; w<MAX_CONC; w++)			      /*Do not touch */
	{                                                     /*Do not touch */
	  ep->d_lame_lambda[MAX_VARIABLE_TYPES+w] = dfdC[w]; /*Do not touch */
	}		                                      /*Do not touch */
    }					                      /*Do not touch */

  return(0);
} /* End of usr_lame_lambda */
/*****************************************************************************/

/*
 *  SOLID EXPANSION 
 */
/*
 * int usr_expansion ()
 *
 * ----------------------------------------------------------------------------
 * This routine is responsible for filling up the following portions of the elc structure
 * at the current gauss point:
 *     intput:  param - array of constants input on the property card.  
 *
 *     output:  f       => elc->lame_lambda - shear modulus
 *              dfdT    => elc->d_lambda[a][TEMPERATURE] 
 *                                         - derivative wrt teelcerature.
 *              dfdC[i] => elc->d_lambda[MAX_VARIABLE_TYPES+i]
 *                                         - derivative wrt mass frac species i
 *              dfdV[0] => elc->d_lambda[VELOCITY1]
 *                            elc->d_lambda[VELOCITY2]
 *                            elc->d_lambda[VELOCITY3]
 *                                         - derivative wrt velocities
 *              dfdX[0] => elc->d_lambda[MESH_DISPLACEMENT1]
 *                            elc->d_lambda[MESH_DISPLACEMENT2]
 *                            elc->d_lambda[MESH_DISPLACEMENT3]
 *                                         - derivative wrt mesh displacements
 *
 *   NB: The user need only supply f, dfdT, dfdC, etc....elc struct is loaded up for you
 */

int
usr_expansion(dbl *param, 	 /* ptr to user-defined parameter list        */
		double *thermexp, 
		double d_thermexp_dx[MAX_VARIABLE_TYPES+MAX_CONC])	
{
  int a, b;
  int w;

  dbl f, dfdT;   /* momentum sources and its derivative wrt temperature*/
  dbl dfdV[DIM];      /* momentum source derivative wrt velocity*/
  dbl dfdC[MAX_CONC]; /* momentum source derivative wrt concentration*/
  dbl dfdX[DIM];      /* momentum source derivative wrt displacement*/

  dbl X[DIM], T, C[MAX_CONC]; /* Convenient local variables */

  int i;

  /* Begin Execution */

 /**********************************************************/

 /* Comment out our remove this line if using this routine */

/*  EH(-1,"No user_expansion model implemented.");     */

 /**********************************************************/


 /************Initialize everything for saftey**************/
  f = 0;                                   /*Do not touch */
  dfdT = 0;			           /*Do not touch */
  for(b=0; b<DIM; b++) dfdV[b]=0.;         /*Do not touch */
  for(b=0; b<DIM; b++) dfdX[b]=0.;         /*Do not touch */
  for(w=0; w<MAX_CONC; w++) dfdC[w]=0.;    /*Do not touch */
  /**********************************************************/
  
 /***********Load up convenient local variables*************/
 /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T;                                             /*Do not touch */
  for(a=0; a<DIM; a++)X[a] = fv->x[a];		         /*Do not touch */
  for(i=0; i<pd->Num_Species_Eqn; i++) C[i] = fv->c[i];  /*Do not touch */

	f= param[0]*sin(param[1]*tran->time_value);

 /**********************************************************/

 /*******Add property function and sensitivities here*******/
  /* Example:
   *
   * if (mp->porosity < param[1])
   *    {
   *      f = param[0] * pow((1. - mp->porosity)/(1. - param[1]), param[2]);
   *    }
   *  else
   *    {
   *      f = param[0];
   *    }
   *
   * Add sensitivities here.
   * 
   *  if (mp->porosity < param[1])
   *    {
   *      elc->d_lame_lambda[POROSITY] = - param[2] * param[0] * 
   *	pow((1. - mp->porosity), param[2] - 1)/
   *	pow((1. - param[1]), param[2]);
   *    }
   *  else
   *    {
   *      elc->d_lame_lambda[POROSITY] = 0;
   *    }
   */

 /****************Don't touch these lines***********************/ 
  for(a=0; a<DIM; a++)
    { 
      *thermexp = f;                                   /*Do not touch */
      d_thermexp_dx[TEMPERATURE] = dfdT;		      /*Do not touch */
                                                              /*Do not touch */
      for ( b=0; b<DIM; b++)                                  /*Do not touch */
	{						      /*Do not touch */
	  d_thermexp_dx[VELOCITY1+b] = dfdV[b];	      /*Do not touch */
	}                                                     /*Do not touch */
      for ( b=0; b<DIM; b++)                                  /*Do not touch */
	{             					      /*Do not touch */
	  d_thermexp_dx[MESH_DISPLACEMENT1+a] = dfdX[b]; /*Do not touch */
	}						      /*Do not touch */
      for ( w=0; w<MAX_CONC; w++)			      /*Do not touch */
	{                                                     /*Do not touch */
	  d_thermexp_dx[MAX_VARIABLE_TYPES+w] = dfdC[w]; /*Do not touch */
	}		                                      /*Do not touch */
    }					                      /*Do not touch */

  return(0);
} /* End of usr_expansion */
/*****************************************************************************/

int 
usr_diffusivity(int species_no,	/* Species number of diffusivity etc. needed */
		dbl *param)	/* ptr to user-defined parameter list        */
{
  int a;
  int w;

  dbl D, dDdT;        /* diffusivity and its derivative wrt temperature*/
  dbl dDdX[DIM];      /* diffusivity derivative wrt displacement*/
  dbl dDdC[MAX_CONC]; /* diffusivity derivative wrt concentration*/

  dbl X[DIM], T, C[MAX_CONC]; /* Convenient local variables */

  int i;

  dbl a0, a1, a2, a3, a4;


  /* Begin Execution */
 /**********************************************************/

 /* Comment out our remove this line if using this routine */

   EH(-1,"No user_diffusivity model implemented");
 
 /**********************************************************/

  /************Initialize everything for safety**************/
  D = 0.;                                   /*Do not touch */
  dDdT = 0.;                                /*Do not touch */
  for (i=0; i<DIM; i++)                     /*Do not touch */
    {                                       /*Do not touch */
      dDdX[i]=0.;                           /*Do not touch */
    }                                       /*Do not touch */
  for (i=0; i<MAX_CONC; i++)                /*Do not touch */
    {                                       /*Do not touch */
      dDdC[i] = 0.;                         /*Do not touch */
    }                                       /*Do not touch */
 /**********************************************************/

 /***********Load up convenient local variables*************/
 /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T;                                             /*Do not touch */
  for(a=0; a<DIM; a++)X[a] = fv->x[a];		         /*Do not touch */
  for(i=0; i<pd->Num_Species_Eqn; i++) C[i] = fv->c[i];  /*Do not touch */

 /**********************************************************/

 /*******Add property function and sensitivities here*******/

 a0 = *param;
 a1 = *(param+1);
 a2 = *(param+2);
 a3 = *(param+3);
 a4 = *(param+4);

/* Constant model 
 D = a0;
 dDdT = 0.0;
 for ( a=0; a<DIM; a++)	dDdX[a]=0.;
 for ( i=0; i<pd->Num_Species_Eqn; i++) dDdC[i]=0.;
*/

/* Linear species model
 D = a0*(1.0 + a1*C[species_no]);
 dDdT = 0.0;
 for ( a=0; a<DIM; a++)	dDdX[a]=0.;
 for ( i=0; i<pd->Num_Species_Eqn; i++) dDdC[i]=0.; dDdC[species_no]=a0*a1;
*/

/* Linear model 
 D = a0*(1.0 + a1*C[species_no] + a2*T + a3*X[0] + a4*X[1]);
 dDdC[species_no]= a0*a1;
 dDdT = a0*a2;
 dDdX[0]= a0*a3;
 dDdX[1]= a0*a4;
*/

 /****************Don't touch these lines***********************/ 
  w = species_no;
  mp->diffusivity[w] = D; 				   /*Do not touch */
  mp->d_diffusivity[w][TEMPERATURE] = dDdT;		   /*Do not touch */
  for ( a=0; a<DIM; a++)				   /*Do not touch */
    {							   /*Do not touch */
      mp->d_diffusivity[w][MESH_DISPLACEMENT1+a] = dDdX[a];/*Do not touch */
    }							   /*Do not touch */
  for ( a=0; a<MAX_CONC; a++)				   /*Do not touch */
    {							   /*Do not touch */
      mp->d_diffusivity[w][MAX_VARIABLE_TYPES+a] = dDdC[a];/*Do not touch */
    }
 /**************************************************************/
  
  return 0;
} /* End of usr_diffusivity */
/*****************************************************************************/

int
usr_FlowingLiquidViscosity(dbl *param) /* ptr to user-defined parameter list */
{
  /* Local Variables */
  int a;
  dbl mu, dmudT;   /* thermal conductivity and its derivative wrt temperature*/
  dbl dmudV[DIM]; /* heat source derivative wrt velocity*/
  dbl dmudC[MAX_CONC]; /* heat source derivative wrt concentration*/
  dbl dmudX[DIM]; /* heat source derivative wrt displacement*/
  dbl X[DIM], F, T, C[MAX_CONC]; /* Convenient local variables */

  /* Begin Execution */

 /**********************************************************/

 /* Comment out our remove this line if using this routine */

   EH(-1,"No user_FlowingLiquidViscosity model implemented");

 /**********************************************************/


 /************Initialize everything for safety**************/
  mu = 0;                                 /*Do not touch */
  dmudT = 0;				  /*Do not touch */
  for(a=0; a<DIM; a++) dmudV[a]=0.; 	  /*Do not touch */
  for(a=0; a<DIM; a++) dmudX[a]=0.;       /*Do not touch */
  for(a=0; a<MAX_CONC; a++) dmudC[a]=0.;  /*Do not touch */
 /**********************************************************/

 /***********Load up convenient local variables*************/
 /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T;                                              /*Do not touch */
  F = fv->F;                                              /*Do not touch */
  for(a=0; a<DIM; a++)X[a] = fv->x[a];		          /*Do not touch */
  for(a=0; a<pd->Num_Species_Eqn; a++) C[a] = fv->c[a];   /*Do not touch */

 /**********************************************************/

 /*******Add property function and sensitivities here*******/
 /* An example:
  *  mu = param[0] * exp(param[1] * (param[2] - C[0]));
  *
  * Add sensitivities here
  *  dmudC[0] = - mu * param[1];
  *
  */
 /****************Don't touch these lines***********************/  
  mp->FlowingLiquid_viscosity = mu;                           /*Do not touch */
  mp->d_FlowingLiquid_viscosity[TEMPERATURE] = dmudT;	      /*Do not touch */
						              /*Do not touch */
  for ( a=0; a<DIM; a++)                                      /*Do not touch */
    {						              /*Do not touch */
      mp->d_FlowingLiquid_viscosity[VELOCITY1+a] = dmudV[a];  /*Do not touch */
    }                                                         /*Do not touch */
  for ( a=0; a<DIM; a++)                                      /*Do not touch */
    {             					      /*Do not touch */
      mp->d_FlowingLiquid_viscosity[MESH_DISPLACEMENT1+a] = dmudX[a];  /*Do not touch */
    }					               	      /*Do not touch */
  for ( a=0; a<MAX_CONC; a++)			              /*Do not touch */
    {                                                         /*Do not touch */
      mp->d_FlowingLiquid_viscosity[MAX_VARIABLE_TYPES+a] = dmudC[a];  /*Do not touch */
    }			                                      /*Do not touch */

  return(0);
} /* End of usr_FlowingLiquidViscosity */
/*****************************************************************************/

/********************************************
 *    general 2D lubrication approximation  *
 ********************************************/
#if defined SECOR_HEAT_FLUX
double  usr_heat_flux(const double gradP[],     /*   pressure gradient  */
                   double q[],                  /*   flow vector            */
                   double dq_gradP[DIM][DIM],   /*   flow sens wrt gradP    */
                   double dq_dX[DIM][DIM], /*   flow sens wrt coords   */
		   const double time,
		   const double h,
		   const double dh_dX[DIM],
		   const double vb[DIM],
		   const double vt[DIM],
                   double dq_dVb[DIM][DIM],
                   double dq_dVt[DIM][DIM]
                   )
{
        double vp0[10],vp1[10],vpm[10],vparm[10],tmp;
        double r[3],xj[3][3],u[3],rdpx[3],rdpy[3],rdvbx[3],rdvby[3],rdvtx[3],rdvty[3];
        static double udpx[3],udpy[3],gradPold[3],udh[3];
        double udvbx[3],udvby[3],udvtx[3],udvty[3];

        double eps, epstol=1.e-3, epsvisc=1.e-1; /*  integral iteration tolerance */
        double pgrad,pgrad2,pgrad3; /* pressure gradient magnitude */
        double shr;
        double visw;
        double e[2],f[2];  /* unit vectors */
        double x0[3],x1[3],xd[3],x[3];  /*  unknowns rate0,rate1,const */
        int i,l,idiv,jdiv,jdi;      /* loop counters */
        int iconv, ifin, ifail, icnt, iter;     /*  continuation parameters */
        double step,th,th1,step_beta;

        double rate0f,rate1f,shear,vis0,vis1,dvis0,dvis1;
        double drate, xint1,xint1old;
        double c0,c1,delx,rate,cee,vis,dvis;
        double ratea, rateb,sheara, shearb, visa, visb, dvisa, dvisb;
        double xint2,xint2old, xint3, xint3old;
        double xint1d,xint2d,xint3d;
        double dvis0dK,dvis1dK,dvisdK,dvisadK,dvisbdK;

        /*  gauss pts and weights  */
        /*double gp[3] = {0.5-sqrt(15.)/10.,0.5,0.5+sqrt(15.)/10.};
        double w[3] = {5./18.,4./9.,5./18.};  */
        double gp[3] = {0.112701665379258,0.5,0.887298334620742};
        double w[3] = {0.277777777777778,0.444444444444444,0.277777777777778};

        double deter, detinv;   /* determinant, inverse   */
        double nt_eps, nt_epsr, nt_epstol=1.e-9;        /*  newton tolerance */

        double term1, term2, term3, term4;
        double termv0, termv1, qxd0, qxd1, qyd0, qyd1;

 	double dqdp,dqdpdh;
        double qxdp, qydp, qxdpx, qxdpy, qydpx, qydpy;
        double wrate;
        double qxdK,qydK;

/*      double gamma_dot[DIM][DIM] = {0,0,0,0,0,0,0,0,0};
        double d_mu_dv[DIM][MDE], d_mu_dmesh[DIM][MDE], d_mu_dT[MDE];
        double d_mu_dp[MDE], d_mu_dC[MAX_CONC][MDE], d_mu_dgd;
        double mu;
        int err;  */

        static int first_call = TRUE;
        static double xold[3];

/*      recover viscosity parameters  */
        vparm[0]=gn->mu0;
        vparm[1]=gn->nexp;
        vparm[2]=gn->lam;
        vparm[3]=gn->aexp;
        vparm[4]=gn->muinf;
        vparm[5]=0.0;
        vparm[6]=10.0;
        vparm[7]=gn->atexp;
        vparm[8]=gn->wlfc2;
        vparm[9]=mp->reference[TEMPERATURE];

        tmp = mp->melting_point_liquidus;

/*  pressure gradient  */

        pgrad=sqrt(gradP[0]*gradP[0] + gradP[1]*gradP[1]);

/*      check if the flow is only drag flow  */

        shr=sqrt(SQUARE(vt[0]-vb[0])+SQUARE(vt[1]-vb[1]))/h ;
        visw=viscar(shr,vparm,tmp);
 	x[0]=shr-0.5*h*pgrad/visw;
 	x[1]=shr+0.5*h*pgrad/visw;
 	vis0=viscar(x[0],vparm,tmp);
 	vis1=viscar(x[1],vparm,tmp);
 	x[2]=visw/h*(-gradP[1]*(vt[0]-vb[0])+gradP[0]*(vt[1]-vb[1]))/pgrad;
 	eps=(fabs(vis1-visw)+fabs(vis0-visw))/visw;
 	if(eps <= epsvisc  )	{
 		dqdp = -h*h*h/(12.*visw);
                dq_gradP[0][0]=dq_gradP[1][1] = dqdp;
                dq_gradP[0][1]=dq_gradP[1][0] = 0. ;
                q[0]=h*(vt[0]+vb[0])/2.+ dqdp*gradP[0];
                q[1]=h*(vt[1]+vb[1])/2. + dqdp*gradP[1];
 		dqdpdh = -3.*h*h/(12.*visw);
         	dq_dX[0][0] = (0.5*(vt[0]+vb[0])+dqdpdh*gradP[0])*dh_dX[0];
         	dq_dX[0][1] = (0.5*(vt[0]+vb[0])+dqdpdh*gradP[0])*dh_dX[1];
         	dq_dX[1][0] = (0.5*(vt[1]+vb[1])+dqdpdh*gradP[1])*dh_dX[0];
         	dq_dX[1][1] = (0.5*(vt[1]+vb[1])+dqdpdh*gradP[1])*dh_dX[1];
                dq_dVb[0][0] = h/2.;  dq_dVb[0][1] = 0.;
                dq_dVb[1][0] = 0.;  dq_dVb[1][1] = h/2.;
                dq_dVt[0][0] = h/2.;  dq_dVt[0][1] = 0.;
                dq_dVt[1][0] = 0.;  dq_dVt[1][1] = h/2.;
         	for(i=0;i<3;i++){xold[i] = x[i];}
 	        for(i=0;i<2;i++){gradPold[i]=gradP[i];}
        		gradPold[2] = h;
 		wrate=MAX(fabs(x[0]),fabs(x[1]));
                return(wrate);
                }
        else
                {

/*      unit vector directions   */

                e[0]=gradP[0]/pgrad;
                e[1]=gradP[1]/pgrad;
                f[0] = -e[1];
                f[1] = e[0];
                }

        for(i=0;i<10;i++)       {
/*      printf("%d  %g\n",i,vparm[i]);  */
                vp0[i]=vp1[i]=vparm[i];}
        vp0[1]=1.0;
        vp0[3]=2.0;
        vp0[5]=0.0;
        vp0[7]=0.0;

/**
        solve for the shear rate at the top and bottom surface
        and crosswise shear stress value

        initial guesses based solutions for Newtonian liquids
**/

        x0[0] = (e[0]*(vt[0]-vb[0])+e[1]*(vt[1]-vb[1]))/h
                        -h*pgrad/(2.*vparm[0]);
        x0[1] = (e[0]*(vt[0]-vb[0])+e[1]*(vt[1]-vb[1]))/h
                        +h*pgrad/(2.*vparm[0]);
        x0[2] = vparm[0]/h*(f[0]*(vt[0]-vb[0])+f[1]*(vt[1]-vb[1]));

        if (first_call)
                {
                for(i=0;i<3;i++){xd[i] = 0.;}
                first_call = FALSE;
                }
                else
                {
                for(i=0;i<3;i++){xd[i]=xold[i]+udpx[i]*(gradP[0]-gradPold[0])
                                +udpy[i]*(gradP[1]-gradPold[1])
                                +udh[i]*(h-gradPold[2]) - x0[i];}
                }

/********************************************
      perform continuation from newtonian to viscosity parameters
*********************************************/

        iconv=1;
        step=1.0;
        ifail=0;
        th=0.;
        ifin=0;
        th1=1.0;

/**   inside continuation loop **/

        icnt=0;
        while ( icnt < 100 && ifin != 1)        {
                icnt++;

/**   check if past target value **/

          if(th+step >= th1){ step=th1-th; ifin=1;}
          th += step;

        for(i=0;i<10;i++){vpm[i]=vp0[i] + th*(vp1[i]-vp0[i]);}

        for(i=0;i<3;i++){x[i] = x0[i] + step*xd[i];}

        /*printf("%g   %g   %g\n",x[0],x[1],x[2]);*/
/**     iterate to determine shear rates **/

        iter=0; nt_eps=100.*nt_epstol;
        while ( iter < 40 && nt_eps > nt_epstol)        {
        iter++;

/*      evaluate viscosity  */

        rate0f=cross_rate(x[0],vpm,x[2],tmp);
        shear=sqrt(x[0]*x[0]+rate0f*rate0f);
        vis0=viscar(shear,vpm,tmp);
        dvis0=viscard(shear,vpm,tmp);
        dvis0dK=rate0f*dvis0/(vis0+rate0f*rate0f*dvis0);

/*      gamma_dot[0][1] = gamma_dot[1][0] = shear;
            err = viscosity(gn, &vis0, gamma_dot, &dvis0, d_mu_dv,
                d_mu_dmesh, d_mu_dT, d_mu_dp, d_mu_dC);
                dvis0 = dvis0/shear;*/

        rate1f=cross_rate(x[1],vpm,x[2],tmp);
        shear=sqrt(x[1]*x[1]+rate1f*rate1f);
        vis1=viscar(shear,vpm,tmp);
        dvis1=viscard(shear,vpm,tmp);
        dvis1dK=rate1f*dvis1/(vis1+rate1f*rate1f*dvis1);

/*      gamma_dot[0][1] = gamma_dot[1][0] = shear;
            err = viscosity(gn, &vis1, gamma_dot, &dvis1, d_mu_dv,
                d_mu_dmesh, d_mu_dT, d_mu_dp, d_mu_dC);
                dvis1 = dvis1/shear; */

/**       evaluate first viscosity integral  */

        if(x[0]*x[1] >= 0.0)    {

/*       rates of the same sign - do as one integral  */

        drate=x[1]-x[0];
        xint1old=0.;
        jdi=0; eps=10.*epstol;
        while (jdi < 8 && eps > epstol)      {
                jdi++;
		jdiv = pow(2,jdi-1);
                xint1=0.;xint1d=0.;
                for(idiv=1;idiv<=jdiv;idiv++)   {
                        c0=((double)(idiv-1))/((double)jdiv);
                        c1=((double)idiv)/((double)jdiv);
                        delx=c1-c0;
                        for(l=0;l<3;l++)        {
                                cee=c0+gp[l]*delx;
                                rate=drate*cee+x[0];

/*      determine cross rate which corresponds to this shear  */

                                shear=cross_rate(rate,vpm,x[2],tmp);
                                shr=sqrt(rate*rate+shear*shear);

/*      gamma_dot[0][1] = gamma_dot[1][0] = shr;
            err = viscosity(gn, &vis, gamma_dot, &d_mu_dgd, d_mu_dv,
                d_mu_dmesh, d_mu_dT, d_mu_dp, d_mu_dC);*/

                                vis=viscar(shr,vpm,tmp);
                                dvis=viscard(shr,vpm,tmp);
                                dvisdK=shear*dvis/(vis+shear*shear*dvis);
                                xint1 += vis*rate*drate*delx*w[l];
                                xint1d += dvisdK*rate*drate*delx*w[l];
                                }
                        }       /*   !  idiv loop  */
                        eps=fabs(xint1-xint1old)/(1.+fabs(xint1));
                        xint1old=xint1;
                }     /*        !   jdiv loop  */

        if(eps > epstol)        {
        fprintf(stderr,"xint1 integration diverged %g\n",eps);
                }
        }
        else    {

/*       rates of opposite sign - do as two integrals   */

        xint1old=0.;
        jdi=0; eps=10.*epstol;
        while (jdi < 8 && eps > epstol)      {
                jdi++;
 		jdiv = pow(2,jdi-1);
                xint1=0.; xint1d=0.;
                for(idiv=1;idiv <= jdiv; idiv++)        {
                        c0=((double)(idiv-1))/((double)jdiv);
                        c1=((double)idiv)/((double)jdiv);
                        delx=c1-c0;
                        for(l=0;l<3;l++)        {
                                cee=c0+gp[l]*delx;
                                ratea=cee*x[1];
                                rateb=cee*x[0];

/*      determine cross rate which corresponds to this shear   */

                                sheara=cross_rate(ratea,vpm,x[2],tmp);
                                shearb=cross_rate(rateb,vpm,x[2],tmp);
                                shr=sqrt(ratea*ratea+sheara*sheara);

/*      gamma_dot[0][1] = gamma_dot[1][0] = shr;
            err = viscosity(gn, &visa, gamma_dot, &d_mu_dgd, d_mu_dv,
                d_mu_dmesh, d_mu_dT, d_mu_dp, d_mu_dC); */

                                visa=viscar(shr,vpm,tmp);
                                dvisa=viscard(shr,vpm,tmp);
                                dvisadK=sheara*dvisa/(visa+sheara*sheara*dvisa);
                                shr=sqrt(rateb*rateb+shearb*shearb);

/*      gamma_dot[0][1] = gamma_dot[1][0] = shr;
            err = viscosity(gn, &visb, gamma_dot, &d_mu_dgd, d_mu_dv,
                d_mu_dmesh, d_mu_dT, d_mu_dp, d_mu_dC);*/

                                visb=viscar(shr,vpm,tmp);
                                dvisb=viscard(shr,vpm,tmp);
                                dvisbdK=shearb*dvisb/(visb+shearb*shearb*dvisb);
                xint1 += (visa*ratea*x[1]-visb*rateb*x[0])*delx*w[l];
                xint1d += (dvisadK*ratea*x[1] - dvisbdK*rateb*x[0])*delx*w[l];
                                }
                        }       /*   !  idiv loop   */
        eps=fabs(xint1-xint1old)/(1.+fabs(xint1));
        xint1old=xint1;
                }       /*      !  jdiv loop   */
        if(eps > epstol)        {
        fprintf(stderr,"xint1 integration diverged %g\n",eps);
                }
        }

/*      evaluate second viscosity integrals   */

        if(x[0]*x[1] >= 0.0)    {

/*       rates of the same sign - do as one integral   */

        drate=x[1]-x[0];
        xint2old=0.;
        jdi=0;  eps=10.*epstol;
        while ( jdi < 8 && eps > epstol )    {
                 jdi++;
 		jdiv = pow(2,jdi-1);
                xint2=0.; xint2d=0.;
                for(idiv=1;idiv<=jdiv;idiv++)   {
                        c0=((double)(idiv-1))/((double)jdiv);
                        c1=((double)idiv)/((double)jdiv);
                        delx=c1-c0;
                        for(l=0;l<3;l++)        {
                                cee=c0+gp[l]*delx;
                                rate=drate*cee+x[0];

/*      determine cross rate which corresponds to this shear  */

                                shear=cross_rate(rate,vpm,x[2],tmp);
                                shr=sqrt(rate*rate+shear*shear);

/*      gamma_dot[0][1] = gamma_dot[1][0] = shr;
            err = viscosity(gn, &vis, gamma_dot, &d_mu_dgd, d_mu_dv,
                d_mu_dmesh, d_mu_dT, d_mu_dp, d_mu_dC);*/

                                vis=viscar(shr,vpm,tmp);
                                dvis=viscard(shr,vpm,tmp);
                                dvisdK=shear*dvis/(vis+shear*shear*dvis);
                                xint2 +=  log(vis)*drate*delx*w[l];
                                xint2d += dvisdK/vis*drate*delx*w[l];
                                }
                        }       /*   !  idiv loop  */
        eps=fabs(xint2-xint2old)/(1.+fabs(xint2));
        xint2old=xint2;
                }    /*  !  jdiv loop    */
        if(eps > epstol)        {
        fprintf(stderr,"xint2 integration diverged %g\n",eps);
                }
        }
        else    {

/*       rates of opposite sign - do as two integrals   */

        xint2old=0.;
        jdi=0;  eps=10.*epstol;
        while (jdi < 8 && eps > epstol )     {
                jdi++;
		jdiv = pow(2,jdi-1);
                xint2=0.; xint2d=0.;
                for(idiv=1;idiv<=jdiv;idiv++)   {
                        c0=((double)(idiv-1))/((double)jdiv);
                        c1=((double)idiv)/((double)jdiv);
                        delx=c1-c0;
                        for(l=0;l<3;l++)        {
                                cee=c0+gp[l]*delx;
                                ratea=cee*x[1];
                                rateb=cee*x[0];

/*      determine cross rate which corresponds to this shear    */

                                sheara=cross_rate(ratea,vpm,x[2],tmp);
                                shearb=cross_rate(rateb,vpm,x[2],tmp);
                                shr=sqrt(ratea*ratea+sheara*sheara);

/*      gamma_dot[0][1] = gamma_dot[1][0] = shr;
            err = viscosity(gn, &visa, gamma_dot, &d_mu_dgd, d_mu_dv,
                d_mu_dmesh, d_mu_dT, d_mu_dp, d_mu_dC);*/

                                visa=viscar(shr,vpm,tmp);
                                dvisa=viscard(shr,vpm,tmp);
                                dvisadK=sheara*dvisa/(visa+sheara*sheara*dvisa);
                                shr=sqrt(rateb*rateb+shearb*shearb);

/*      gamma_dot[0][1] = gamma_dot[1][0] = shr;
            err = viscosity(gn, &visb, gamma_dot, &d_mu_dgd, d_mu_dv,
                d_mu_dmesh, d_mu_dT, d_mu_dp, d_mu_dC);*/

                                visb=viscar(shr,vpm,tmp);
                                dvisb=viscard(shr,vpm,tmp);
                                dvisbdK=shearb*dvisb/(visb+shearb*shearb*dvisb);
                xint2 += (log(visa)*x[1]-log(visb)*x[0])*delx*w[l];
                xint2d += (dvisadK/visa*x[1] - dvisbdK/visb*x[0])*delx*w[l];
                                }
                        }       /*   !   idiv loop   */
                eps=fabs(xint2-xint2old)/(1.+fabs(xint2));
                xint2old=xint2;
                }       /*   !   jdiv loop      */
        if(eps > epstol)        {
        fprintf(stderr,"xint2 integration diverged %g\n",eps);
                }
        }

/*      evaluate residuals, jacobian entries    */

        r[0]=h*pgrad-vis1*x[1]+vis0*x[0];
        r[1]=pgrad*(e[0]*(vt[0]-vb[0])+e[1]*(vt[1]-vb[1]))
                -(vis1*x[1]*x[1]-vis0*x[0]*x[0]-xint1);
        r[2]=pgrad*(f[0]*(vt[0]-vb[0])+f[1]*(vt[1]-vb[1]))
                -x[2]*(x[1]*(1.+log(vis1))-x[0]*(1.+log(vis0))
                -xint2);

        /*printf("resid %g %g %g\n",r[0],r[1],r[2]);*/
/*      jacobian entries   */

        xj[0][0]=vis0+x[0]*x[0]*dvis0;
        xj[0][1]=-vis1-x[1]*x[1]*dvis1;
        xj[0][2]=-x[1]*dvis1dK+x[0]*dvis0dK;
        xj[1][0]=vis0*x[0]+x[0]*x[0]*x[0]*dvis0;
        xj[1][1]=-vis1*x[1]-x[1]*x[1]*x[1]*dvis1;
        xj[1][2]=-x[1]*x[1]*dvis1dK+x[0]*x[0]*dvis0dK+xint1d;
        xj[2][0]=x[2]*(1.+x[0]*x[0]*dvis0/vis0);
        xj[2][1]=-x[2]*(1.+x[1]*x[1]*dvis1/vis1);
        xj[2][2]=-(x[1]*(1.+log(vis1))-x[0]*(1.+log(vis0))-xint2)
                        -x[2]*(x[1]*dvis1dK/vis1-x[0]*dvis0dK/vis0-xint2d);

/*      solve 3x3 linear system         */

        deter=xj[0][0]*(xj[1][1]*xj[2][2]-xj[2][1]*xj[1][2])
              - xj[1][0]*(xj[0][1]*xj[2][2]-xj[2][1]*xj[0][2])
              + xj[2][0]*(xj[1][2]*xj[0][1]-xj[1][1]*xj[0][2]);
        detinv=1./deter;
        u[0]=detinv*(-r[0]*(xj[2][2]*xj[1][1]-xj[2][1]*xj[1][2])
                - r[1]*(xj[2][1]*xj[0][2]-xj[2][2]*xj[0][1])
                - r[2]*(xj[1][2]*xj[0][1]-xj[1][1]*xj[0][2]));
        u[1]=detinv*(-r[0]*(xj[2][0]*xj[1][2]-xj[2][2]*xj[1][0])
                - r[1]*(xj[2][2]*xj[0][0]-xj[2][0]*xj[0][2])
                - r[2]*(xj[1][0]*xj[0][2]-xj[1][2]*xj[0][0]));
        u[2]=detinv*(-r[0]*(xj[2][1]*xj[1][0]-xj[2][0]*xj[1][1])
                -r[1]*(xj[2][0]*xj[0][1]-xj[2][1]*xj[0][0])
                -r[2]*(xj[1][1]*xj[0][0]-xj[1][0]*xj[0][1]));
        for(i=0;i<3;i++)        {x[i] += u[i];}
        nt_epsr=sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
        nt_eps=sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2])
                /(1.+sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]));

        if( (iter > 3 && nt_eps > 1.1) || (nt_epsr > 1.0e+10) ) { break; }
        if(icnt > 10)   {
                printf("usr_heat_flux:iter eps x[] %d %g  %g %g %g %g\n"
                        ,iter,nt_eps,nt_epsr,x[0],x[1],x[2]);
                        }

        }       /*   newton iteration loop */
        if(nt_eps < nt_epstol ) { iconv=1;}
                else    {iconv = 0;}

/*     if not converged : cut step, go back   */

         if(iconv == 0) {
                ifin=0;
                th -= step;
                step *=  0.1;
                ifail +=1;
                if(ifail >= 6)  {
                        fprintf(stderr,"6 consecutive failures  \n");
                fprintf(stderr,"icnt  ifail gradP %d %d %g %g %g %g\n",
                                icnt,ifail,th,th1,gradP[0],gradP[1]);
                        exit(0);
                        }
                }
         else   {
                ifail=0;

/*    compute continuation derivative   */

         if(step != 0.0)        {
                for(i=0;i<3;i++){xd[i] = (x[i]-x0[i])/step;}
                }
          else  {
                for(i=0;i<3;i++){xd[i] = 0.0;}
                }
        for(i=0;i<3;i++)        {x0[i] = x[i];}

/*      if converged; compute new step  */
           step_beta=pow(10.,((pow(2.,7-iter)-1.)/4.));
           step_beta=MAX(0.333,step_beta);
           step_beta=MIN(3.0,step_beta);
           step=step*step_beta;
                }  /*   end of iconv conditional  */
        if(icnt > 20)
        {printf("usr_heat_flux:icnt th step %d %g %g\n",icnt,th,step);}

        }       /*   end of icnt while loop       */
        if(ifin != 1)   {
                fprintf(stderr,"flow2d continuation loop failed \n");
                fprintf(stderr,"icnt  ifin gradP %d %d %g %g %g %g\n",
                                icnt,ifin,th,th1,gradP[0],gradP[1]);
                exit(0);
                }
                /*printf("cont loop %d %d\n",icnt,ifin);*/

/*      solve for flowrate, knowing both wall shear rates and
        the value of the crosswise shear stress  */


        if(x[0]*x[1] >= 0.0)    {

/*       rates of the same sign - do as one integral  */

        drate=x[1]-x[0];
        xint3old=0.;
        jdi=0;  eps=10.*epstol;
        while ( jdi < 8 && eps > epstol)     {
                jdi++;
		jdiv = pow(2,jdi-1);
                xint3=0.;  xint3d=0.;
                for(idiv=1;idiv<=jdiv;idiv++)   {
                        c0=((double)(idiv-1))/((double)jdiv);
                        c1=((double)idiv)/((double)jdiv);
                        delx=c1-c0;
                        for(l=0;l<3;l++)        {
                                cee=c0+gp[l]*delx;
                                rate=drate*cee+x[0];

/*      determine cross rate which corresponds to this shear   */

                                shear=cross_rate(rate,vpm,x[2],tmp);
                                shr=sqrt(rate*rate+shear*shear);

/*      gamma_dot[0][1] = gamma_dot[1][0] = shr;
            err = viscosity(gn, &vis, gamma_dot, &d_mu_dgd, d_mu_dv,
                d_mu_dmesh, d_mu_dT, d_mu_dp, d_mu_dC);*/

                                vis=viscar(shr,vpm,tmp);
                                dvis=viscard(shr,vpm,tmp);
                                dvisdK=shear*dvis/(vis+shear*shear*dvis);
                                xint3 += SQUARE(vis*rate)*drate*delx*w[l];
                                xint3d += 2.*vis*dvisdK*SQUARE(rate)*drate*delx*w[l
];
                                }
                        }       /*   !  idiv loop   */
        eps=fabs(xint3-xint3old)/(1.+fabs(xint3));
        xint3old=xint3;
                }       /*   !   jdiv loop      */
        if(eps > epstol)        {
        fprintf(stderr,"xint3 integration diverged %g\n",eps);
                }
        }
        else    {

/*       rates of opposite sign - do as two integrals   */

        xint3old=0.;
        jdi=0;  eps=10.*epstol;
        while (jdi < 8 && eps > epstol )     {
                jdi++;
		jdiv = pow(2,jdi-1);
                xint3=0.;  xint3d=0.;
                for(idiv=1;idiv<=jdiv;idiv++)   {
                        c0=((double)(idiv-1))/((double)jdiv);
                        c1=((double)idiv)/((double)jdiv);
                        delx=c1-c0;
                        for(l=0;l<3;l++)        {
                                cee=c0+gp[l]*delx;
                                ratea=cee*x[1];
                                rateb=cee*x[0];

/*      determine cross rate which corresponds to this shear   */

                                sheara=cross_rate(ratea,vpm,x[2],tmp);
                                shearb=cross_rate(rateb,vpm,x[2],tmp);
                                shr=sqrt(ratea*ratea+sheara*sheara);

/*      gamma_dot[0][1] = gamma_dot[1][0] = shr;
            err = viscosity(gn, &visa, gamma_dot, &d_mu_dgd, d_mu_dv,
                d_mu_dmesh, d_mu_dT, d_mu_dp, d_mu_dC); */

                                visa=viscar(shr,vpm,tmp);
                                dvisa=viscard(shr,vpm,tmp);
                                dvisadK=sheara*dvisa/(visa+sheara*sheara*dvisa);
                                shr=sqrt(rateb*rateb+shearb*shearb);

/*      gamma_dot[0][1] = gamma_dot[1][0] = shr;
            err = viscosity(gn, &visb, gamma_dot, &d_mu_dgd, d_mu_dv,
                d_mu_dmesh, d_mu_dT, d_mu_dp, d_mu_dC);*/

                                visb=viscar(shr,vpm,tmp);
                                dvisb=viscard(shr,vpm,tmp);
                                dvisbdK=shearb*dvisb/(visb+shearb*shearb*dvisb);
                                xint3 += (SQUARE(visa*ratea)*x[1]
                                        -SQUARE(visb*rateb)*x[0])*delx*w[l];
                                xint3d += (2.*visa*dvisadK*SQUARE(ratea)*x[1]
                                -2.*visb*dvisbdK*SQUARE(rateb)*x[0])*delx*w[l];
                                }
                        }       /*   !  idiv loop     */
        eps=fabs(xint3-xint3old)/(1.+fabs(xint3));
        xint3old=xint3;
                }       /*   !  jdiv loop       */
        }
        if(eps > epstol)        {
        fprintf(stderr,"xint3 integration diverged %g\n",eps);
                }



        term1 = vt[0]*vis1*x[1]-vb[0]*vis0*x[0];
        term2 = vis1*vis1*x[1]*x[1]*x[1]-vis0*vis0*x[0]*x[0]*x[0]-xint3;
        term3 = vis1*x[1]*x[1]-vis0*x[0]*x[0]-xint1;
        term4 = vt[1]*vis1*x[1]-vb[1]*vis0*x[0];
        pgrad2=pgrad*pgrad;
        pgrad3=pgrad2*pgrad;

/*      flow vector     */

        q[0] = term1/pgrad -e[0]*term2/2./pgrad2- f[0]*x[2]*term3/pgrad2;
        q[1] = term4/pgrad -e[1]*term2/2./pgrad2- f[1]*x[2]*term3/pgrad2;
        termv0 = vis0+x[0]*x[0]*dvis0;
        termv1 = vis1+x[1]*x[1]*dvis1;

/*      compute derivatives wrt rate0,rate1   */

        qxd0 = (-vb[0]*termv0)/pgrad
                        - e[0]*(-2.*vis0*x[0]*x[0]*termv0)/(2.*pgrad2)
                        - f[0]*x[2]*(-x[1]*termv0)/pgrad2;
        qxd1 = (vt[0]*termv1)/pgrad
                        - e[0]*(2.*vis1*x[1]*x[1]*termv1)/(2.*pgrad2)
                        - f[0]*x[2]*(x[1]*termv1)/pgrad2;
        qyd0 = (-vb[1]*termv0)/pgrad
                        - e[1]*(-2.*vis0*x[0]*x[0]*termv0)/(2.*pgrad2)
                        - f[1]*x[2]*(-x[0]*termv0)/pgrad2;
        qyd1 = (vt[1]*termv1)/pgrad
                        -e[1]*(2.*vis1*x[1]*x[1]*termv1)/(2.*pgrad2)
                        - f[1]*x[2]*(x[1]*termv1)/pgrad2;
        qxdK = (vt[0]*x[1]*dvis1dK-vb[0]*x[0]*dvis0dK)/pgrad
                -(vis1*dvis1dK*x[1]*x[1]*x[1]-vis0*dvis0dK*x[0]*x[0]*x[0]
                -0.5*xint3d)*e[0]/pgrad2
                -f[0]*term3/pgrad2
        -(vt[1]*x[1]*dvis1dK-vb[1]*x[0]*dvis0dK-xint1d)*f[0]*x[2]/pgrad2;
        qydK = (vt[1]*x[1]*dvis1dK-vb[1]*x[0]*dvis0dK)/pgrad
                -(vis1*dvis1dK*x[1]*x[1]*x[1]-vis0*dvis0dK*x[0]*x[0]*x[0]
                -0.5*xint3d)*e[1]/pgrad2
                -f[1]*term3/pgrad2
        -(vt[1]*x[1]*dvis1dK-vb[1]*x[0]*dvis0dK-xint1d)*f[1]*x[2]/pgrad2;


/*      derivatives wrt pgrad   */

        qxdp = -term1/pgrad2 + e[0]*term2/pgrad3
                + 2.*f[0]*x[2]*term3/pgrad3;
        qydp = -term4/pgrad2 + e[1]*term2/pgrad3
                + 2.*f[1]*x[2]*term3/pgrad3;

/*      derivatives wrt px,py through ex,ey,fx,fy       */

        qxdpx = -e[1]*e[1]*term2/2./pgrad3 -e[0]*e[1]*x[2]*term3/pgrad3;
        qxdpy = e[0]*e[1]*term2/2./pgrad3 +e[0]*e[0]*x[2]*term3/pgrad3;
        qydpx = e[0]*e[1]*term2/2./pgrad3 -e[1]*e[1]*x[2]*term3/pgrad3;
        qydpy = -e[0]*e[0]*term2/2./pgrad3 +e[0]*e[1]*x[2]*term3/pgrad3;

 /*   compute sensitivities of the variables to px,py  */

        rdpx[0]=h*e[0]; rdpy[0]=h*e[1];
        rdpx[1]=vt[0]-vb[0];
        rdpy[1]=vt[1]-vb[1];
        rdpx[2]=vt[1]-vb[1];
        rdpy[2]=-(vt[0]-vb[0]);
        rdvtx[1] = gradP[0];  rdvty[1] = gradP[1];
        rdvbx[1] = -gradP[0];  rdvby[1] = -gradP[1];
        rdvtx[2] = -gradP[1];  rdvty[2] = gradP[0];
        rdvbx[2] = gradP[1];  rdvby[2] = -gradP[0];
        udpx[0]=detinv*(-rdpx[0]*(xj[2][2]*xj[1][1]-xj[2][1]*xj[1][2])
                - rdpx[1]*(xj[2][1]*xj[0][2]-xj[2][2]*xj[0][1])
                - rdpx[2]*(xj[1][2]*xj[0][1]-xj[1][1]*xj[0][2]));
        udpx[1]=detinv*(-rdpx[0]*(xj[2][0]*xj[1][2]-xj[2][2]*xj[1][0])
                - rdpx[1]*(xj[2][2]*xj[0][0]-xj[2][0]*xj[0][2])
                - rdpx[2]*(xj[1][0]*xj[0][2]-xj[1][2]*xj[0][0]));
        udpx[2]=detinv*(-rdpx[0]*(xj[2][1]*xj[1][0]-xj[2][0]*xj[1][1])
                -rdpx[1]*(xj[2][0]*xj[0][1]-xj[2][1]*xj[0][0])
                -rdpx[2]*(xj[1][1]*xj[0][0]-xj[1][0]*xj[0][1]));
        udpy[0]=detinv*(-rdpy[0]*(xj[2][2]*xj[1][1]-xj[2][1]*xj[1][2])
                - rdpy[1]*(xj[2][1]*xj[0][2]-xj[2][2]*xj[0][1])
                - rdpy[2]*(xj[1][2]*xj[0][1]-xj[1][1]*xj[0][2]));
        udpy[1]=detinv*(-rdpy[0]*(xj[2][0]*xj[1][2]-xj[2][2]*xj[1][0])
                - rdpy[1]*(xj[2][2]*xj[0][0]-xj[2][0]*xj[0][2])
                - rdpy[2]*(xj[1][0]*xj[0][2]-xj[1][2]*xj[0][0]));
        udpy[2]=detinv*(-rdpy[0]*(xj[2][1]*xj[1][0]-xj[2][0]*xj[1][1])
                -rdpy[1]*(xj[2][0]*xj[0][1]-xj[2][1]*xj[0][0])
                -rdpy[2]*(xj[1][1]*xj[0][0]-xj[1][0]*xj[0][1]));
        udh[0]=detinv*(-pgrad*(xj[2][2]*xj[1][1]-xj[2][1]*xj[1][2]));
        udh[1]=detinv*(-pgrad*(xj[2][0]*xj[1][2]-xj[2][2]*xj[1][0]));
        udh[2]=detinv*(-pgrad*(xj[2][1]*xj[1][0]-xj[2][0]*xj[1][1]));
        udvtx[0]=detinv*(-rdvtx[1]*(xj[2][1]*xj[0][2]-xj[2][2]*xj[0][1])
                - rdvtx[2]*(xj[1][2]*xj[0][1]-xj[1][1]*xj[0][2]));
        udvtx[1]=detinv*(-rdvtx[1]*(xj[2][2]*xj[0][0]-xj[2][0]*xj[0][2])
                - rdvtx[2]*(xj[1][0]*xj[0][2]-xj[1][2]*xj[0][0]));
        udvtx[2]=detinv*(-rdvtx[1]*(xj[2][0]*xj[0][1]-xj[2][1]*xj[0][0])
                -rdvtx[2]*(xj[1][1]*xj[0][0]-xj[1][0]*xj[0][1]));
        udvty[0]=detinv*(-rdvty[1]*(xj[2][1]*xj[0][2]-xj[2][2]*xj[0][1])
                - rdvty[2]*(xj[1][2]*xj[0][1]-xj[1][1]*xj[0][2]));
        udvty[1]=detinv*(-rdvty[1]*(xj[2][2]*xj[0][0]-xj[2][0]*xj[0][2])
                - rdvty[2]*(xj[1][0]*xj[0][2]-xj[1][2]*xj[0][0]));
        udvty[2]=detinv*(-rdvty[1]*(xj[2][0]*xj[0][1]-xj[2][1]*xj[0][0])
                -rdvty[2]*(xj[1][1]*xj[0][0]-xj[1][0]*xj[0][1]));
        udvbx[0]=detinv*(-rdvbx[1]*(xj[2][1]*xj[0][2]-xj[2][2]*xj[0][1])
                - rdvbx[2]*(xj[1][2]*xj[0][1]-xj[1][1]*xj[0][2]));
        udvbx[1]=detinv*(-rdvbx[1]*(xj[2][2]*xj[0][0]-xj[2][0]*xj[0][2])
                - rdvbx[2]*(xj[1][0]*xj[0][2]-xj[1][2]*xj[0][0]));
        udvbx[2]=detinv*(-rdvbx[1]*(xj[2][0]*xj[0][1]-xj[2][1]*xj[0][0])
                -rdvbx[2]*(xj[1][1]*xj[0][0]-xj[1][0]*xj[0][1]));
        udvby[0]=detinv*(-rdvby[1]*(xj[2][1]*xj[0][2]-xj[2][2]*xj[0][1])
                - rdvby[2]*(xj[1][2]*xj[0][1]-xj[1][1]*xj[0][2]));
        udvby[1]=detinv*(-rdvby[1]*(xj[2][2]*xj[0][0]-xj[2][0]*xj[0][2])
                - rdvby[2]*(xj[1][0]*xj[0][2]-xj[1][2]*xj[0][0]));
        udvby[2]=detinv*(-rdvby[1]*(xj[2][0]*xj[0][1]-xj[2][1]*xj[0][0])
                -rdvby[2]*(xj[1][1]*xj[0][0]-xj[1][0]*xj[0][1]));

        dq_gradP[0][0] = qxd0*udpx[0]+qxd1*udpx[1]+qxdp*e[0]+qxdK*udpx[2]+qxdpx;
        dq_gradP[0][1] = qxd0*udpy[0]+qxd1*udpy[1]+qxdp*e[1]+qxdK*udpy[2]+qxdpy;
        dq_gradP[1][0] = qyd0*udpx[0]+qyd1*udpx[1]+qydp*e[0]+qydK*udpx[2]+qydpx;
        dq_gradP[1][1] = qyd0*udpy[0]+qyd1*udpy[1]+qydp*e[1]+qydK*udpy[2]+qydpy;

        dq_dX[0][0] = (qxd0*udh[0]+qxd1*udh[1]+qxdK*udh[2])*dh_dX[0];
        dq_dX[0][1] = (qxd0*udh[0]+qxd1*udh[1]+qxdK*udh[2])*dh_dX[1];
        dq_dX[1][0] = (qyd0*udh[0]+qyd1*udh[1]+qydK*udh[2])*dh_dX[0];
        dq_dX[1][1] = (qyd0*udh[0]+qyd1*udh[1]+qydK*udh[2])*dh_dX[1];

        dq_dVb[0][0] = -vis0*x[0]/pgrad + qxd0*udvbx[0]+qxd1*udvbx[1]+qxdK*udvbx[2];
        dq_dVb[0][1] = qxd0*udvby[0]+qxd1*udvby[1]+qxdK*udvby[2];
        dq_dVb[1][0] = qyd0*udvbx[0]+qyd1*udvbx[1]+qydK*udvbx[2];
        dq_dVb[1][1] = -vis0*x[0]/pgrad + qyd0*udvby[0]+qyd1*udvby[1]+qydK*udvby[2];

        dq_dVt[0][0] = vis1*x[1]/pgrad + qxd0*udvtx[0]+qxd1*udvtx[1]+qxdK*udvtx[2];
        dq_dVt[0][1] = qxd0*udvty[0]+qxd1*udvty[1]+qxdK*udvty[2];
        dq_dVt[1][0] = qyd0*udvtx[0]+qyd1*udvtx[1]+qydK*udvtx[2];
        dq_dVt[1][1] = vis1*x[1]/pgrad + qyd0*udvty[0]+qyd1*udvty[1]+qydK*udvty[2];

/**  max shear rate  **/

        sheara=cross_rate(x[0],vpm,x[2],tmp);
        shearb=cross_rate(x[1],vpm,x[2],tmp);
        wrate=sqrt(sheara*sheara+x[0]*x[0]);
        wrate=MAX(wrate,sqrt(shearb*shearb+x[1]*x[1]));
        /*printf("done %g\n",wrate);  */

        for(i=0;i<3;i++){xold[i] = x[i];}
        for(i=0;i<2;i++){gradPold[i]=gradP[i];}
        gradPold[2] = h;

        return(wrate);
}
/************************************
        iteration routine for finding cross shear rate
************************************/

double cross_rate(const double rate,  /* shear rate in main direction  */
                  const double vpm[], /* viscosity parameters */
                  const double strs,  /* shear stress */
                  const double tmp    /*  temperature */
                  )
{
        double epstol=1.0e-9;  /*  convergence tolerance  */
        double shr, vis, crsrate;
        double res,visd,xjac,delta,eps;
        int iter;

/**     guess initial value for cross shear rate   **/

        shr=fabs(rate);
        vis=viscar(shr,vpm,tmp);
        crsrate=strs/vis;

/** newton iteration loop **/

        iter=0;eps=10*epstol;
        while (iter <= 40 && eps > epstol)      {

        iter++;
        shr=sqrt(rate*rate+crsrate*crsrate);
        vis=viscar(shr,vpm,tmp);
        res=strs-vis*crsrate;
        visd=viscard(shr,vpm,tmp);
        xjac=-vis-crsrate*crsrate*visd;
        delta=-res/xjac;
        crsrate += delta;
        eps = fabs(delta)/(1.+fabs(crsrate));
                }

        if ( eps > epstol )     {
                fprintf(stderr,"iteration diverged - cross_rate\n");
        printf("crsrate rate strs vis %g %g %g %g\n",crsrate,rate,strs,vis);
                exit(1);
                }
                else
                {
                return(crsrate);
                }
}

/*************************************
  viscosity and derivative functions
     (yield stress plus thinning)
************************************/

double viscar(const double shr,  /*  shear rate  */
              const double vparm[],     /*  viscosity parameters  */
              const double tmp          /*  temperature  */
              )
{
        double eta0,xn,xlamb,xfact,etainf,yield,fparm;
        double tparm1,tparm2,tref,at,shr1,tp,tp1,tp2,tp3,tpe;
        double vis1;

        eta0=vparm[0];
        xn=vparm[1];
        xlamb=vparm[2];
        xfact=vparm[3];
        etainf=vparm[4];
        yield=vparm[5];
        fparm=vparm[6];

        tparm1=vparm[7];
        tparm2=vparm[8];
        tref=vparm[9];

/*      calculate temperature shift factor  */

        at=exp(tparm1*(tref-tmp)/(tparm2+tmp-tref));
        shr1=fabs(shr);

/*      calculate exponential term  */

        tp=at*shr1*fparm;
        if(tp > 0.1)    { tp1=(1.-exp(-tp))/(at*shr1);}
        else    { tp1=fparm*(1.+tp*(-0.5+tp*(1./6.-tp/24.)));}

        tpe=eta0-etainf+yield*tp1;
        tp2=at*xlamb*shr1;

        if(tp2 >= pow(10.0,8./xfact))   { tp3=pow(tp2,xn-1.);}
        else    { tp3=pow(1.+pow(tp2,xfact),(xn-1.)/xfact); }

        vis1 = at*(etainf+tpe*tp3);
        return(at*(etainf+tpe*tp3));
}



double viscard(const double shr,  /*  shear rate  */
              const double vparm[],     /*  viscosity parameters  */
              const double tmp          /*  temperature  */
              )
{
        double eta0,xn,xlamb,xfact,etainf,yield,fparm;
        double tparm1,tparm2,tref,at,shr1,tp,tp1,tp2,tp3,tp4,tp5,tpe;
        double return_value;

        eta0=vparm[0];
        xn=vparm[1];
        xlamb=vparm[2];
        xfact=vparm[3];
        etainf=vparm[4];
        yield=vparm[5];
        fparm=vparm[6];

        tparm1=vparm[7];
        tparm2=vparm[8];
        tref=vparm[9];

/*      calculate temperature shift factor  */

        at=exp(tparm1*(tref-tmp)/(tparm2+tmp-tref));
        shr1=fabs(shr);

/*      calculate exponential term  */

        tp=at*shr1*fparm;
        if(tp > 0.1)    {
        tp1=(1.-exp(-tp))/(at*shr1);
        tp2=exp(-tp)*(1.+tp)-1.;
                }
        else    {
        tp1=fparm*(1.+tp*(-0.5+tp*(1./6.-tp/24.)));
        tp2=tp*tp*(-0.5+tp*(1./3.+tp*(-1./8.+tp/24.)));
                }

        if(shr > 0.0)   {
                tpe=eta0-etainf+yield*tp1;
                tp3=at*xlamb*shr1;
                if(tp3 >= pow(10.0,8/xfact))    {
                        tp4=pow(tp3,xn-1.);
                        tp5=pow(tp3,xn-3.);
                        }
        else    {
                tp4=pow(1+pow(tp3,xfact),(xn-1.)/xfact);
                tp5=pow((1.+pow(tp3,xfact)),(xn-1.-xfact)/xfact)*pow(tp3,xfact-2.);
                }

        return_value = at*tpe*(xn-1.)*tp5*SQUARE(at*xlamb)
                                +tp4*yield*tp2/(shr1*shr1*shr1);
                }
        if(shr == 0.0)  {return_value=0.0;}
           return(return_value);
}
#else
double  usr_heat_flux(const double gradP[],     /*   pressure gradient  */
                   double q[],                  /*   flow vector            */
                   double dq_gradP[DIM][DIM],   /*   flow sens wrt gradP    */
                   double dq_dX[DIM][DIM], /*   flow sens wrt coords   */
                   const double time)
{
        EH(-1,"No usr_heat_flux model supplied");
        return(1.);
} /* End of usr_heat_flux */
#endif


/****************************************************
 *    user-defined porous media permeability model  *
 ****************************************************/
int
usr_permeability(dbl *param) /* user-defined parameter list */
{
 dbl k,phi,dkdphi;        /* permeability and its derivative wrt porosity*/

 /* Comment out or remove this line if using this routine */
                                                                                
 /*   EH(-1,"No user_defined  permeability model implemented");  */

 /**********************************************************/

 
 /************Initialize everything for safety**************/
  k = 0.;                                                  /*Do not touch */
  dkdphi = 0.;                                             /*Do not touch */
 /**********************************************************/

 /***********Load up convenient variables*************/

  phi = mp->porosity;                                      /*Do not touch */

 /**********************************************************/

 /*******Add property function and sensitivities here*******/

  if(phi<param[3])
    { 
      k=param[2];
      dkdphi=0.0; 
    }
  else 
    { 
      k = param[0]+(param[2]-param[0])/(param[3]-param[1])*(phi-param[1]);      
      dkdphi = (param[2]-param[0])/(param[3]-param[1]);  
    }
   
 /****************Don't touch these lines***********************/ 
  mp->permeability = k;         			      /*Do not touch */
  mp->d_permeability[POR_POROSITY] = dkdphi;	              /*Do not touch */
 /**************************************************************/

  return(0);
}
/*****************************************************************************/
/*****************************************************************************/
/*
 * Yield Stress value
 */
/*
 * int usr_yield_stress ()
 *
 * ----------------------------------------------------------------------------
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
 */

int 
usr_yield_stress(dbl *param, dbl time)	/* pointer to user-defined parameter list    */
{
  /* Local Variables */
  int a;
  dbl tau_y, dtau_ydT;   /* thermal conductivity and its derivative wrt temperature*/
  dbl dtau_ydV[DIM]; /* heat source derivative wrt velocity*/
  dbl dtau_ydC[MAX_CONC]; /* heat source derivative wrt concentration*/
  dbl dtau_ydX[DIM]; /* heat source derivative wrt displacement*/
  dbl X[DIM], F, T, C[MAX_CONC]; /* Convenient local variables */
  dbl yield_oo;
  dbl yield_0;
  dbl time_const;
  dbl x0, y0, radius, dist, dist_factor;

  /* Begin Execution */

 /**********************************************************/

 /* Comment out our remove this line if using this routine */

 /*EH(-1,"No user_defined yield stress model implemented");  */

 /**********************************************************/


 /************Initialize everything for saftey**************/
  tau_y = 0;                                  /*Do not touch */
  dtau_ydT = 0;				  /*Do not touch */
  for(a=0; a<DIM; a++) dtau_ydV[a]=0.; 	  /*Do not touch */
  for(a=0; a<DIM; a++) dtau_ydX[a]=0.;       /*Do not touch */
  for(a=0; a<MAX_CONC; a++) dtau_ydC[a]=0.;  /*Do not touch */
 /**********************************************************/

 /***********Load up convenient local variables*************/
 /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T;                                        /*Do not touch */
  F = fv->F;                                        /*Do not touch */
  for(a=0; a<DIM; a++)X[a] = fv->x[a];		    /*Do not touch */
  for(a=0; a<pd->Num_Species_Eqn; a++) C[a] = fv->c[a];    /*Do not touch */

 /**********************************************************/

 /*******Add property function and sensitivities here*******/
 /* An example:
  *  mu = param[0] * exp(param[1] * (param[2] - C[0]));
  *
  * Add sensitivities here
  *  dmudC[0] = - mu * param[1];
  *
  */

  yield_oo=param[0];
  yield_0=param[1];
  time_const=param[2];

/*  Position-dependent yield  referenced to the undeformed coordinates  */
 
  if(gn->len_u_tau_y > 3)	
	{
 	x0 = param[3];
 	y0 = param[4];
 	radius = param[5];
 	dist = 1-sqrt(SQUARE(fv->x0[0]-x0)+SQUARE(fv->x0[1]-y0))/radius;
 	dist_factor = param[gn->len_u_tau_y-1];
 	for(a=gn->len_u_tau_y-2 ; a>5 ; a--)	
 		{
 		dist_factor = dist_factor*dist + param[a];
 		}
 	dist_factor = dist_factor*dist;
 	}
  tau_y = yield_oo + (yield_0-yield_oo)*exp(-time/time_const)*exp(-dist_factor);

 /****************Don't touch these lines***********************/  
  gn->tau_y = tau_y;                                    /*Do not touch */
							 /*Do not touch */
  return(0);
} /* End of usr_viscosity */

/*****************************************************************************/
/* End of file user_mp.c */
/*****************************************************************************/
