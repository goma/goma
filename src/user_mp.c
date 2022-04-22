/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2022 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/

/*
 *$Id: user_mp.c,v 5.5 2009-01-12 16:41:49 hkmoffa Exp $
 */

/* Standard include files */

#include <math.h>

/* GOMA include files */

#include "el_elm.h"
#include "mm_as.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "mm_mp.h"
#include "mm_mp_structs.h"
#include "rf_fem_const.h"
#include "std.h"
#include "user_mp.h"

#define GOMA_USER_MP_C

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
 * an error.  If the routine is used, then the call to  GOMA_EH (error handler routine) routine
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
 *	GOMA_EH(GOMA_ERROR,"No user-defined function for this material");
 *    }
 * ----------------------------------------------------------------------------
 */
// ignore warnings for user function
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"

int usr_thermal_conductivity(dbl *param, dbl time) /* user-defined parameter list */
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

  GOMA_EH(GOMA_ERROR, "No user_defined thermal conductivity model implemented");
  /**********************************************************/

  /************Initialize everything for safety**************/
  k = 0.;                        /*Do not touch */
  dkdT = 0.;                     /*Do not touch */
  for (i = 0; i < DIM; i++)      /*Do not touch */
  {                              /*Do not touch */
    dkdX[i] = 0.;                /*Do not touch */
  }                              /*Do not touch */
  for (i = 0; i < MAX_CONC; i++) /*Do not touch */
  {                              /*Do not touch */
    dkdC[i] = 0.;                /*Do not touch */
  }                              /*Do not touch */
                                 /**********************************************************/

  /***********Load up convenient local variables*************/
  /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T; /*Do not touch */
  for (a = 0; a < DIM; a++)
    X[a] = fv->x[a]; /*Do not touch */
  for (i = 0; i < pd->Num_Species_Eqn; i++)
    C[i] = fv->c[i]; /*Do not touch */

  /**********************************************************/

  /*******Add property function and sensitivities here*******/

  /****************Don't touch these lines***********************/
  mp->thermal_conductivity = k;                                   /*Do not touch */
  mp->d_thermal_conductivity[TEMPERATURE] = dkdT;                 /*Do not touch */
  for (a = 0; a < DIM; a++)                                       /*Do not touch */
  {                                                               /*Do not touch */
    mp->d_thermal_conductivity[MESH_DISPLACEMENT1 + a] = dkdX[a]; /*Do not touch */
  }                                                               /*Do not touch */
  for (a = 0; a < MAX_CONC; a++)                                  /*Do not touch */
  {                                                               /*Do not touch */
    mp->d_thermal_conductivity[MAX_VARIABLE_TYPES + a] = dkdC[a]; /*Do not touch */
  }
  /**************************************************************/
  return (0);
} /* End of usr_thermal_conductivity */
/*****************************************************************************/

int usr_electrical_conductivity(dbl *param, dbl time) /* user-defined parameter list */
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

  GOMA_EH(GOMA_ERROR, "No user_defined electrical conductivity model implemented");
  /**********************************************************/

  /************Initialize everything for safety**************/
  k = 0.;                        /*Do not touch */
  dkdT = 0.;                     /*Do not touch */
  dkdV = 0.;                     /*Do not touch */
  for (i = 0; i < DIM; i++)      /*Do not touch */
  {                              /*Do not touch */
    dkdX[i] = 0.;                /*Do not touch */
  }                              /*Do not touch */
  for (i = 0; i < MAX_CONC; i++) /*Do not touch */
  {                              /*Do not touch */
    dkdC[i] = 0.;                /*Do not touch */
  }                              /*Do not touch */
                                 /**********************************************************/

  /***********Load up convenient local variables*************/
  /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T; /*Do not touch */
  V = fv->V; /*Do not touch */
  for (a = 0; a < DIM; a++)
    X[a] = fv->x[a]; /*Do not touch */
  for (i = 0; i < pd->Num_Species_Eqn; i++)
    C[i] = fv->c[i]; /*Do not touch */

  /**********************************************************/

  /*******Add property function and sensitivities here*******/

  /****************Don't touch these lines***********************/
  mp->electrical_conductivity = k;                                   /*Do not touch */
  mp->d_electrical_conductivity[TEMPERATURE] = dkdT;                 /*Do not touch */
  mp->d_electrical_conductivity[VOLTAGE] = dkdV;                     /*Do not touch */
  for (a = 0; a < DIM; a++)                                          /*Do not touch */
  {                                                                  /*Do not touch */
    mp->d_electrical_conductivity[MESH_DISPLACEMENT1 + a] = dkdX[a]; /*Do not touch */
  }                                                                  /*Do not touch */
  for (a = 0; a < MAX_CONC; a++)                                     /*Do not touch */
  {                                                                  /*Do not touch */
    mp->d_electrical_conductivity[MAX_VARIABLE_TYPES + a] = dkdC[a]; /*Do not touch */
  }
  /**************************************************************/

  return (0);
} /* End of usr_electrical_conductivity */
/*****************************************************************************/

int usr_density(dbl *param) /* pointer to user-defined parameter list    */
{
  /* Local Variables */
  dbl rho, d_rho_dT;      /* density and its derivative wrt temperature*/
  dbl d_rho_dC[MAX_CONC]; /* density derivative wrt concentration*/
  dbl T, C[MAX_CONC];     /* Convenient local variables */
  int w;

  /* Begin Execution */
  /**********************************************************/

  /* Comment out our remove this line if using this routine */

  GOMA_EH(GOMA_ERROR, "No user_defined density model implemented");
  /**********************************************************/

  /************Initialize everything for safety**************/
  rho = 0.;                      /*Do not touch */
  d_rho_dT = 0.;                 /*Do not touch */
  for (w = 0; w < MAX_CONC; w++) /*Do not touch */
  {                              /*Do not touch */
    d_rho_dC[w] = 0.;            /*Do not touch */
  }                              /*Do not touch */
                                 /**********************************************************/

  /***********Load up convenient local variables*************/
  /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T; /*Do not touch */
  for (w = 0; w < pd->Num_Species_Eqn; w++)
    C[w] = fv->c[w]; /*Do not touch */

  /**********************************************************/

  /*******Add property function and sensitivities here*******/

  /****************Don't touch these lines***********************/
  mp->density = rho;                                     /*Do not touch */
  mp->d_density[TEMPERATURE] = d_rho_dT;                 /*Do not touch */
  for (w = 0; w < MAX_CONC; w++)                         /*Do not touch */
  {                                                      /*Do not touch */
    mp->d_density[MAX_VARIABLE_TYPES + w] = d_rho_dC[w]; /*Do not touch */
  }
  /**************************************************************/

  return (0);
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
 *	GOMA_EH(GOMA_ERROR,"No user-defined function for this material");
 *    }
 * ----------------------------------------------------------------------------
 */
int usr_heat_capacity(dbl *param, dbl time) /* pt to user-defined parameter list */
{
  int a;

  dbl Cp, dCpdT;       /* heat capacity and its derivative wrt temperature*/
  dbl dCpdX[DIM];      /* heat capacity derivative wrt displacement */
  dbl dCpdC[MAX_CONC]; /* heat capacity derivative wrt concentration*/

  dbl X[DIM], T, C[MAX_CONC]; /* Convenient local variables */

  int i;

  /* Begin Execution */
  /**********************************************************/

  /* Comment out our remove this line if using this routine */

  GOMA_EH(GOMA_ERROR, "No user_defined heat capacity model implemented");
  /**********************************************************/

  /************Initialize everything for saftey**************/
  Cp = 0;                        /*Do not touch */
  dCpdT = 0;                     /*Do not touch */
  for (a = 0; a < DIM; a++)      /*Do not touch */
  {                              /*Do not touch */
    dCpdX[a] = 0.;               /*Do not touch */
  }                              /*Do not touch */
  for (i = 0; i < MAX_CONC; i++) /*Do not touch */
  {                              /*Do not touch */
    dCpdC[i] = 0.;               /*Do not touch */
  }                              /*Do not touch */
                                 /**********************************************************/

  /***********Load up convenient local variables*************/
  /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T; /*Do not touch */
  for (a = 0; a < DIM; a++)
    X[a] = fv->x[a]; /*Do not touch */
  for (i = 0; i < pd->Num_Species_Eqn; i++)
    C[i] = fv->c[i]; /*Do not touch */

  /**********************************************************/

  /*******Add property function and sensitivities here*******/

  /****************Don't touch these lines***********************/
  /* NB nothing can depend on gradients of variables */
  mp->heat_capacity = Cp;                                   /*Do not touch */
  mp->d_heat_capacity[TEMPERATURE] = dCpdT;                 /*Do not touch */
  for (a = 0; a < DIM; a++)                                 /*Do not touch */
  {                                                         /*Do not touch */
    mp->d_heat_capacity[MESH_DISPLACEMENT1 + a] = dCpdX[a]; /*Do not touch */
  }                                                         /*Do not touch */
  for (a = 0; a < MAX_CONC; a++)                            /*Do not touch */
  {                                                         /*Do not touch */
    mp->d_heat_capacity[MAX_VARIABLE_TYPES + a] = dCpdC[a]; /*Do not touch */
  }                                                         /*Do not touch */
  /**********************************************************/

  return (0);
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
 *	GOMA_EH(GOMA_ERROR,"No user-defined function for this material");
 *    }
 * ----------------------------------------------------------------------------
 */
int usr_heat_source(dbl *param, dbl time) /* ptr to the user-defined parameter list */
{
  int a;

  dbl h, dhdT, dhdV;  /*heat sourceand its derivative wrt temperature, voltage*/
  dbl dhdv[DIM];      /* heat source derivative wrt velocity*/
  dbl dhdC[MAX_CONC]; /* heat source derivative wrt concentration*/
  dbl dhdX[DIM];      /* heat source derivative wrt displacement*/

  int i;
  dbl X[DIM], T, V, C[MAX_CONC]; /* Convenient local variables */

  /* Begin Execution */

  /**********************************************************/

  /* Comment out our remove this line if using this routine */

  GOMA_EH(GOMA_ERROR, "No user_defined heat source model implemented");
  /**********************************************************/

  /************Initialize everything for saftey**************/
  h = 0;    /*Do not touch */
  dhdT = 0; /*Do not touch */
  dhdV = 0; /*Do not touch */
  for (a = 0; a < DIM; a++)
    dhdv[a] = 0.; /*Do not touch */
  for (a = 0; a < DIM; a++)
    dhdX[a] = 0.; /*Do not touch */
  for (a = 0; a < MAX_CONC; a++)
    dhdC[a] = 0.; /*Do not touch */
                  /**********************************************************/

  /***********Load up convenient local variables*************/
  /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T; /*Do not touch */
  V = fv->V; /*Do not touch */
  for (a = 0; a < DIM; a++)
    X[a] = fv->x[a]; /*Do not touch */
  for (i = 0; i < pd->Num_Species_Eqn; i++)
    C[i] = fv->c[i]; /*Do not touch */

  /**********************************************************/

  /*******Add property function and sensitivities here*******/

  /****************Don't touch these lines***********************/
  mp->heat_source = h;                                   /*Do not touch */
  mp->d_heat_source[TEMPERATURE] = dhdT;                 /*Do not touch */
  mp->d_heat_source[VOLTAGE] = dhdV;                     /*Do not touch */
                                                         /*Do not touch */
  for (a = 0; a < DIM; a++)                              /*Do not touch */
  {                                                      /*Do not touch */
    mp->d_heat_source[VELOCITY1 + a] = dhdv[a];          /*Do not touch */
  }                                                      /*Do not touch */
  for (a = 0; a < DIM; a++)                              /*Do not touch */
  {                                                      /*Do not touch */
    mp->d_heat_source[MESH_DISPLACEMENT1 + a] = dhdX[a]; /*Do not touch */
  }                                                      /*Do not touch */
  for (a = 0; a < MAX_CONC; a++)                         /*Do not touch */
  {                                                      /*Do not touch */
    mp->d_heat_source[MAX_VARIABLE_TYPES + a] = dhdC[a]; /*Do not touch */
  }                                                      /*Do not touch */

  return (0);
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
 *	GOMA_EH(GOMA_ERROR,"No user-defined function for this material");
 *    }
 * ----------------------------------------------------------------------------
 */

int usr_species_source(int species_no, /* Current species number                 */
                       dbl *param)     /* pointer to user-defined parameter list    */

{
  int a;

  dbl s, dsdT, dsdV;  /*species sourceand its derivative wrt temperature, voltage*/
  dbl dsdv[DIM];      /* species source derivative wrt velocity*/
  dbl dsdC[MAX_CONC]; /* species source derivative wrt concentration*/
  dbl dsdX[DIM];      /* species source derivative wrt displacement*/

  int i;
  dbl X[DIM], T, V, C[MAX_CONC]; /* Convenient local variables */

  /* Begin Execution */

  /**********************************************************/

  /* Comment out our remove this line if using this routine */

  GOMA_EH(GOMA_ERROR, "No user_defined species source model implemented");
  /**********************************************************/

  /************Initialize everything for saftey**************/
  s = 0;    /*Do not touch */
  dsdT = 0; /*Do not touch */
  dsdV = 0; /*Do not touch */
  for (a = 0; a < DIM; a++)
    dsdv[a] = 0.; /*Do not touch */
  for (a = 0; a < DIM; a++)
    dsdX[a] = 0.; /*Do not touch */
  for (a = 0; a < MAX_CONC; a++)
    dsdC[a] = 0.; /*Do not touch */
                  /**********************************************************/

  /***********Load up convenient local variables*************/
  /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T; /*Do not touch */
  V = fv->V; /*Do not touch */
  for (a = 0; a < DIM; a++)
    X[a] = fv->x[a]; /*Do not touch */
  for (i = 0; i < pd->Num_Species_Eqn; i++)
    C[i] = fv->c[i]; /*Do not touch */

  /**********************************************************/

  /*******Add property function and sensitivities here*******/

  /****************Don't touch these lines***********************/
  mp->species_source[species_no] = s;                       /*Do not touch */
  mp->d_species_source[TEMPERATURE] = dsdT;                 /*Do not touch */
  mp->d_species_source[VOLTAGE] = dsdV;                     /*Do not touch */
                                                            /*Do not touch */
  for (a = 0; a < DIM; a++)                                 /*Do not touch */
  {                                                         /*Do not touch */
    mp->d_species_source[VELOCITY1 + a] = dsdv[a];          /*Do not touch */
  }                                                         /*Do not touch */
  for (a = 0; a < DIM; a++)                                 /*Do not touch */
  {                                                         /*Do not touch */
    mp->d_species_source[MESH_DISPLACEMENT1 + a] = dsdX[a]; /*Do not touch */
  }                                                         /*Do not touch */
  for (a = 0; a < MAX_CONC; a++)                            /*Do not touch */
  {                                                         /*Do not touch */
    mp->d_species_source[MAX_VARIABLE_TYPES + a] = dsdC[a]; /*Do not touch */
  }                                                         /*Do not touch */

  return (0);
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
 *	GOMA_EH(GOMA_ERROR,"No user-defined function for this material");
 *    }
 * ----------------------------------------------------------------------------
 */

int usr_current_source(dbl *param) /* pointer to user-defined parameter list */
{
  int a;

  dbl h, dhdT;        /* current source and its derivative wrt temperature*/
  dbl dhdV;           /* current source derivative wrt voltage */
  dbl dhdv[DIM];      /* current source derivative wrt velocity*/
  dbl dhdC[MAX_CONC]; /* current source derivative wrt concentration*/
  dbl dhdX[DIM];      /* current source derivative wrt displacement*/

  int i;
  dbl X[DIM], T, V, C[MAX_CONC]; /* Convenient local variables */

  /* Begin Execution */

  /**********************************************************/

  /* Comment out our remove this line if using this routine */

  GOMA_EH(GOMA_ERROR, "No user_defined current source model implemented");
  /**********************************************************/

  /************Initialize everything for saftey**************/
  h = 0;    /*Do not touch */
  dhdT = 0; /*Do not touch */
  dhdV = 0; /*Do not touch */
  for (a = 0; a < DIM; a++)
    dhdv[a] = 0.; /*Do not touch */
  for (a = 0; a < DIM; a++)
    dhdX[a] = 0.; /*Do not touch */
  for (a = 0; a < MAX_CONC; a++)
    dhdC[a] = 0.; /*Do not touch */
                  /**********************************************************/

  /***********Load up convenient local variables*************/
  /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T; /*Do not touch */
  V = fv->V; /*Do not touch */
  for (a = 0; a < DIM; a++)
    X[a] = fv->x[a]; /*Do not touch */
  for (i = 0; i < pd->Num_Species_Eqn; i++)
    C[i] = fv->c[i]; /*Do not touch */

  /**********************************************************/

  /*******Add property function and sensitivities here*******/

  /****************Don't touch these lines***********************/
  mp->current_source = h;                                   /*Do not touch */
  mp->d_current_source[TEMPERATURE] = dhdT;                 /*Do not touch */
  mp->d_current_source[VOLTAGE] = dhdV;                     /*Do not touch */
                                                            /*Do not touch */
  for (a = 0; a < DIM; a++)                                 /*Do not touch */
  {                                                         /*Do not touch */
    mp->d_current_source[VELOCITY1 + a] = dhdv[a];          /*Do not touch */
  }                                                         /*Do not touch */
  for (a = 0; a < DIM; a++)                                 /*Do not touch */
  {                                                         /*Do not touch */
    mp->d_current_source[MESH_DISPLACEMENT1 + a] = dhdX[a]; /*Do not touch */
  }                                                         /*Do not touch */
  for (a = 0; a < MAX_CONC; a++)                            /*Do not touch */
  {                                                         /*Do not touch */
    mp->d_current_source[MAX_VARIABLE_TYPES + a] = dhdC[a]; /*Do not touch */
  }                                                         /*Do not touch */

  return (0);
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
 *	GOMA_EH(GOMA_ERROR,"No user-defined function for this material");
 *    }
 * ----------------------------------------------------------------------------
 */

int usr_viscosity(dbl *param) /* pointer to user-defined parameter list    */
{
  /* Local Variables */
  int a;
  dbl mu, dmudT;       /* thermal conductivity and its derivative wrt temperature*/
  dbl dmudV[DIM];      /* heat source derivative wrt velocity*/
  dbl dmudC[MAX_CONC]; /* heat source derivative wrt concentration*/
  dbl dmudX[DIM];      /* heat source derivative wrt displacement*/

  dbl X[DIM], F, T, C[MAX_CONC]; /* Convenient local variables */

  /* Begin Execution */

  /**********************************************************/

  /* Comment out our remove this line if using this routine */

  GOMA_EH(GOMA_ERROR, "No user_defined viscosity model implemented");

  /**********************************************************/

  /************Initialize everything for saftey**************/
  mu = 0;    /*Do not touch */
  dmudT = 0; /*Do not touch */
  for (a = 0; a < DIM; a++)
    dmudV[a] = 0.; /*Do not touch */
  for (a = 0; a < DIM; a++)
    dmudX[a] = 0.; /*Do not touch */
  for (a = 0; a < MAX_CONC; a++)
    dmudC[a] = 0.; /*Do not touch */
                   /**********************************************************/

  /***********Load up convenient local variables*************/
  /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T; /*Do not touch */
  F = fv->F; /*Do not touch */
  for (a = 0; a < DIM; a++)
    X[a] = fv->x[a]; /*Do not touch */
  for (a = 0; a < pd->Num_Species_Eqn; a++)
    C[a] = fv->c[a]; /*Do not touch */

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
  if (fv->c[0] < 0.5) {
    mu = param[1];
  }

  /****************Don't touch these lines***********************/
  mp->viscosity = mu;                                   /*Do not touch */
  mp->d_viscosity[TEMPERATURE] = dmudT;                 /*Do not touch */
                                                        /*Do not touch */
  for (a = 0; a < DIM; a++)                             /*Do not touch */
  {                                                     /*Do not touch */
    mp->d_viscosity[VELOCITY1 + a] = dmudV[a];          /*Do not touch */
  }                                                     /*Do not touch */
  for (a = 0; a < DIM; a++)                             /*Do not touch */
  {                                                     /*Do not touch */
    mp->d_viscosity[MESH_DISPLACEMENT1 + a] = dmudX[a]; /*Do not touch */
  }                                                     /*Do not touch */
  for (a = 0; a < MAX_CONC; a++)                        /*Do not touch */
  {                                                     /*Do not touch */
    mp->d_viscosity[MAX_VARIABLE_TYPES + a] = dmudC[a]; /*Do not touch */
  }                                                     /*Do not touch */

  return (0);
} /* End of usr_viscosity */
/*****************************************************************************/

int usr_surface_tension(dbl *param) /* ptr to user-defined parameter list        */
{
  int a;

  dbl sigma, dsigmadT;    /* surface tension and its derivative wrt temperature*/
  dbl dsigmadV[DIM];      /* surface tension derivative wrt velocity*/
  dbl dsigmadC[MAX_CONC]; /* surface tension derivative wrt concentration*/
  dbl dsigmadX[DIM];      /* surface tension derivative wrt displacement*/

  int i;
  dbl X[DIM], T, C[MAX_CONC]; /* Convenient local variables */

  /* Begin Execution */

  /**********************************************************/

  /* Comment out our remove this line if using this routine */

  GOMA_EH(GOMA_ERROR, "No user_defined surface_tension model implemented");

  /**********************************************************/

  /************Initialize everything for saftey**************/
  sigma = 0;    /*Do not touch */
  dsigmadT = 0; /*Do not touch */
  for (a = 0; a < DIM; a++)
    dsigmadV[a] = 0.; /*Do not touch */
  for (a = 0; a < DIM; a++)
    dsigmadX[a] = 0.; /*Do not touch */
  for (a = 0; a < MAX_CONC; a++)
    dsigmadC[a] = 0.; /*Do not touch */
                      /**********************************************************/

  /***********Load up convenient local variables*************/
  /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T; /*Do not touch */
  for (a = 0; a < DIM; a++)
    X[a] = fv->x[a]; /*Do not touch */
  for (i = 0; i < pd->Num_Species_Eqn; i++)
    C[i] = fv->c[i]; /*Do not touch */

  /*******Add property function and sensitivities here*******/
  /* An example:
   *  sigma = param[0] + (param[1] - param[0]) * C[0];
   *
   * Add sensitivities here
   *
   *  dsigmadC[0] =  param[1] - param[0];
   *
   */

  /****************Don't touch these lines***********************/
  mp->surface_tension = sigma;                                   /*Do not touch */
  mp->d_surface_tension[TEMPERATURE] = dsigmadT;                 /*Do not touch */
                                                                 /*Do not touch */
  for (a = 0; a < DIM; a++)                                      /*Do not touch */
  {                                                              /*Do not touch */
    mp->d_surface_tension[VELOCITY1 + a] = dsigmadV[a];          /*Do not touch */
  }                                                              /*Do not touch */
  for (a = 0; a < DIM; a++)                                      /*Do not touch */
  {                                                              /*Do not touch */
    mp->d_surface_tension[MESH_DISPLACEMENT1 + a] = dsigmadX[a]; /*Do not touch */
  }                                                              /*Do not touch */
  for (a = 0; a < MAX_CONC; a++)                                 /*Do not touch */
  {                                                              /*Do not touch */
    mp->d_surface_tension[MAX_VARIABLE_TYPES + a] = dsigmadC[a]; /*Do not touch */
  }                                                              /*Do not touch */

  return (0);
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
 *	GOMA_EH(GOMA_ERROR,"No user-defined function for this material");
 *    }
 * ----------------------------------------------------------------------------
 */

int usr_momentum_source(dbl *param) /* ptr to user-defined parameter list        */
{
  int a, b;
  int w;

  dbl f[DIM], dfdT[DIM];   /* momentum sources and its derivative wrt temperature*/
  dbl dfdV[DIM][DIM];      /* momentum source derivative wrt velocity*/
  dbl dfdC[DIM][MAX_CONC]; /* momentum source derivative wrt concentration*/
  dbl dfdX[DIM][DIM];      /* momentum source derivative wrt displacement*/

  int i;
  dbl X[DIM], T, C[MAX_CONC]; /* Convenient local variables */

  /* Begin Execution */

  /**********************************************************/

  /* Comment out our remove this line if using this routine */

  GOMA_EH(GOMA_ERROR, "No user_momentum_source model implemented");
  /**********************************************************/

  /************Initialize everything for saftey**************/
  for (a = 0; a < DIM; a++) {
    f[a] = 0;    /*Do not touch */
    dfdT[a] = 0; /*Do not touch */
    for (b = 0; b < DIM; b++)
      dfdV[a][b] = 0.; /*Do not touch */
    for (b = 0; b < DIM; b++)
      dfdX[a][b] = 0.; /*Do not touch */
    for (w = 0; w < MAX_CONC; w++)
      dfdC[a][w] = 0.; /*Do not touch */
  }
  /**********************************************************/

  /***********Load up convenient local variables*************/
  /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T; /*Do not touch */
  for (a = 0; a < DIM; a++)
    X[a] = fv->x[a]; /*Do not touch */
  for (i = 0; i < pd->Num_Species_Eqn; i++)
    C[i] = fv->c[i]; /*Do not touch */

  /**********************************************************/

  /*******Add property function and sensitivities here*******/

  /****************Don't touch these lines***********************/
  for (a = 0; a < DIM; a++) {
    mp->momentum_source[a] = f[a];                                   /*Do not touch */
    mp->d_momentum_source[a][TEMPERATURE] = dfdT[a];                 /*Do not touch */
                                                                     /*Do not touch */
    for (b = 0; b < DIM; b++)                                        /*Do not touch */
    {                                                                /*Do not touch */
      mp->d_momentum_source[a][VELOCITY1 + b] = dfdV[a][b];          /*Do not touch */
    }                                                                /*Do not touch */
    for (b = 0; b < DIM; b++)                                        /*Do not touch */
    {                                                                /*Do not touch */
      mp->d_momentum_source[a][MESH_DISPLACEMENT1 + a] = dfdX[a][b]; /*Do not touch */
    }                                                                /*Do not touch */
    for (w = 0; w < MAX_CONC; w++)                                   /*Do not touch */
    {                                                                /*Do not touch */
      mp->d_momentum_source[a][MAX_VARIABLE_TYPES + w] = dfdC[a][w]; /*Do not touch */
    }
  } /*Do not touch */

  return (0);
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

int usr_lame_mu(struct Elastic_Constitutive *ep,
                dbl *param) /* ptr to user-defined parameter list        */
{
  int a, b;
  int w;

  dbl f, dfdT;        /* momentum sources and its derivative wrt temperature*/
  dbl dfdV[DIM];      /* momentum source derivative wrt velocity*/
  dbl dfdC[MAX_CONC]; /* momentum source derivative wrt concentration*/
  dbl dfdX[DIM];      /* momentum source derivative wrt displacement*/

  int i;
  dbl X[DIM], T, C[MAX_CONC]; /* Convenient local variables */

  /* Begin Execution */

  /**********************************************************/

  /* Comment out our remove this line if using this routine */

  GOMA_EH(GOMA_ERROR, "No user_lame_mu model implemented.");

  /**********************************************************/

  /************Initialize everything for saftey**************/
  f = 0;    /*Do not touch */
  dfdT = 0; /*Do not touch */
  for (b = 0; b < DIM; b++)
    dfdV[b] = 0.; /*Do not touch */
  for (b = 0; b < DIM; b++)
    dfdX[b] = 0.; /*Do not touch */
  for (w = 0; w < MAX_CONC; w++)
    dfdC[w] = 0.; /*Do not touch */
  /**********************************************************/

  /***********Load up convenient local variables*************/
  /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T; /*Do not touch */
  for (a = 0; a < DIM; a++)
    X[a] = fv->x[a]; /*Do not touch */
  for (i = 0; i < pd->Num_Species_Eqn; i++)
    C[i] = fv->c[i]; /*Do not touch */

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
   *  dbl dist;
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

  /* Aluminum lame mu vs. T */
  // if(fv->T <= 500)
  //   {
  //    f = param[0]*(7.6e11 + fv->T*(5e11-9.2e11)/300.);
  //    dfdT = (5e11-7.6e11)/300.;
  //  }
  // else
  // {
  //    f = param[0]*(7.6e11 + 500*(5e11-9.2e11)/300.);
  //     dfdT = 0.;
  // }

  /****************Don't touch these lines***********************/
  for (a = 0; a < DIM; a++) {
    ep->lame_mu = f;                                   /*Do not touch */
    ep->d_lame_mu[TEMPERATURE] = dfdT;                 /*Do not touch */
                                                       /*Do not touch */
    for (b = 0; b < DIM; b++)                          /*Do not touch */
    {                                                  /*Do not touch */
      ep->d_lame_mu[VELOCITY1 + b] = dfdV[b];          /*Do not touch */
    }                                                  /*Do not touch */
    for (b = 0; b < DIM; b++)                          /*Do not touch */
    {                                                  /*Do not touch */
      ep->d_lame_mu[MESH_DISPLACEMENT1 + a] = dfdX[b]; /*Do not touch */
    }                                                  /*Do not touch */
    for (w = 0; w < MAX_CONC; w++)                     /*Do not touch */
    {                                                  /*Do not touch */
      ep->d_lame_mu[MAX_VARIABLE_TYPES + w] = dfdC[w]; /*Do not touch */
    }                                                  /*Do not touch */
  }                                                    /*Do not touch */

  return (0);
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

int usr_lame_lambda(struct Elastic_Constitutive *ep,
                    dbl *param) /* ptr to user-defined parameter list        */
{
  int a, b;
  int w;

  dbl f, dfdT;        /* momentum sources and its derivative wrt temperature*/
  dbl dfdV[DIM];      /* momentum source derivative wrt velocity*/
  dbl dfdC[MAX_CONC]; /* momentum source derivative wrt concentration*/
  dbl dfdX[DIM];      /* momentum source derivative wrt displacement*/

  int i;
  dbl X[DIM], T, C[MAX_CONC]; /* Convenient local variables */

  /* Begin Execution */

  /**********************************************************/

  /* Comment out our remove this line if using this routine */

  GOMA_EH(GOMA_ERROR, "No user_lame_lambda model implemented.");

  /**********************************************************/

  /************Initialize everything for saftey**************/
  f = 0;    /*Do not touch */
  dfdT = 0; /*Do not touch */
  for (b = 0; b < DIM; b++)
    dfdV[b] = 0.; /*Do not touch */
  for (b = 0; b < DIM; b++)
    dfdX[b] = 0.; /*Do not touch */
  for (w = 0; w < MAX_CONC; w++)
    dfdC[w] = 0.; /*Do not touch */
  /**********************************************************/

  /***********Load up convenient local variables*************/
  /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T; /*Do not touch */
  for (a = 0; a < DIM; a++)
    X[a] = fv->x[a]; /*Do not touch */
  for (i = 0; i < pd->Num_Species_Eqn; i++)
    C[i] = fv->c[i]; /*Do not touch */

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

  /* Aluminum lame mu vs. T */
  //  f = param[0]*(7.6e11 + fv->T*(5e11-7.6e11)/300.);
  //  dfdT = (5e11-7.6e11)/300.;

  /****************Don't touch these lines***********************/
  for (a = 0; a < DIM; a++) {
    ep->lame_lambda = f;                                   /*Do not touch */
    ep->d_lame_lambda[TEMPERATURE] = dfdT;                 /*Do not touch */
                                                           /*Do not touch */
    for (b = 0; b < DIM; b++)                              /*Do not touch */
    {                                                      /*Do not touch */
      ep->d_lame_lambda[VELOCITY1 + b] = dfdV[b];          /*Do not touch */
    }                                                      /*Do not touch */
    for (b = 0; b < DIM; b++)                              /*Do not touch */
    {                                                      /*Do not touch */
      ep->d_lame_lambda[MESH_DISPLACEMENT1 + a] = dfdX[b]; /*Do not touch */
    }                                                      /*Do not touch */
    for (w = 0; w < MAX_CONC; w++)                         /*Do not touch */
    {                                                      /*Do not touch */
      ep->d_lame_lambda[MAX_VARIABLE_TYPES + w] = dfdC[w]; /*Do not touch */
    }                                                      /*Do not touch */
  }                                                        /*Do not touch */

  return (0);
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

int usr_expansion(dbl *param, /* ptr to user-defined parameter list        */
                  double *thermexp,
                  double d_thermexp_dx[MAX_VARIABLE_TYPES + MAX_CONC]) {
  int a, b;
  int w;

  dbl f, dfdT;        /* momentum sources and its derivative wrt temperature*/
  dbl dfdV[DIM];      /* momentum source derivative wrt velocity*/
  dbl dfdC[MAX_CONC]; /* momentum source derivative wrt concentration*/
  dbl dfdX[DIM];      /* momentum source derivative wrt displacement*/

  int i;
  dbl X[DIM], T, C[MAX_CONC]; /* Convenient local variables */

  /* Begin Execution */

  /**********************************************************/

  /* Comment out our remove this line if using this routine */

  GOMA_EH(GOMA_ERROR, "No user_expansion model implemented.");

  /**********************************************************/

  /************Initialize everything for saftey**************/
  f = 0;    /*Do not touch */
  dfdT = 0; /*Do not touch */
  for (b = 0; b < DIM; b++)
    dfdV[b] = 0.; /*Do not touch */
  for (b = 0; b < DIM; b++)
    dfdX[b] = 0.; /*Do not touch */
  for (w = 0; w < MAX_CONC; w++)
    dfdC[w] = 0.; /*Do not touch */
  /**********************************************************/

  /***********Load up convenient local variables*************/
  /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T; /*Do not touch */
  for (a = 0; a < DIM; a++)
    X[a] = fv->x[a]; /*Do not touch */
  for (i = 0; i < pd->Num_Species_Eqn; i++)
    C[i] = fv->c[i]; /*Do not touch */

  f = param[0] * sin(param[1] * tran->time_value);

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
  for (a = 0; a < DIM; a++) {
    *thermexp = f;                                     /*Do not touch */
    d_thermexp_dx[TEMPERATURE] = dfdT;                 /*Do not touch */
                                                       /*Do not touch */
    for (b = 0; b < DIM; b++)                          /*Do not touch */
    {                                                  /*Do not touch */
      d_thermexp_dx[VELOCITY1 + b] = dfdV[b];          /*Do not touch */
    }                                                  /*Do not touch */
    for (b = 0; b < DIM; b++)                          /*Do not touch */
    {                                                  /*Do not touch */
      d_thermexp_dx[MESH_DISPLACEMENT1 + a] = dfdX[b]; /*Do not touch */
    }                                                  /*Do not touch */
    for (w = 0; w < MAX_CONC; w++)                     /*Do not touch */
    {                                                  /*Do not touch */
      d_thermexp_dx[MAX_VARIABLE_TYPES + w] = dfdC[w]; /*Do not touch */
    }                                                  /*Do not touch */
  }                                                    /*Do not touch */

  return (0);
} /* End of usr_expansion */
/*****************************************************************************/

int usr_diffusivity(int species_no, /* Species number of diffusivity etc. needed */
                    dbl *param)     /* ptr to user-defined parameter list        */
{
  int a;
  int w;

  dbl D, dDdT;        /* diffusivity and its derivative wrt temperature*/
  dbl dDdX[DIM];      /* diffusivity derivative wrt displacement*/
  dbl dDdC[MAX_CONC]; /* diffusivity derivative wrt concentration*/

  dbl X[DIM], T, C[MAX_CONC]; /* Convenient local variables */

  int i;

  /*  dbl a0, a1, a2, a3, a4; */

  /* Begin Execution */
  /**********************************************************/

  /* Comment out our remove this line if using this routine */

  GOMA_EH(GOMA_ERROR, "No user_diffusivity model implemented");

  /**********************************************************/

  /************Initialize everything for safety**************/
  D = 0.;                        /*Do not touch */
  dDdT = 0.;                     /*Do not touch */
  for (i = 0; i < DIM; i++)      /*Do not touch */
  {                              /*Do not touch */
    dDdX[i] = 0.;                /*Do not touch */
  }                              /*Do not touch */
  for (i = 0; i < MAX_CONC; i++) /*Do not touch */
  {                              /*Do not touch */
    dDdC[i] = 0.;                /*Do not touch */
  }                              /*Do not touch */
                                 /**********************************************************/

  /***********Load up convenient local variables*************/
  /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T; /*Do not touch */
  for (a = 0; a < DIM; a++)
    X[a] = fv->x[a]; /*Do not touch */
  for (i = 0; i < pd->Num_Species_Eqn; i++)
    C[i] = fv->c[i]; /*Do not touch */

  /**********************************************************/

  /*******Add property function and sensitivities here*******/

  /*
 a0 = *param;
 a1 = *(param+1);
 a2 = *(param+2);
 a3 = *(param+3);
 a4 = *(param+4);
  */
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
  mp->diffusivity[w] = D;                                   /*Do not touch */
  mp->d_diffusivity[w][TEMPERATURE] = dDdT;                 /*Do not touch */
  for (a = 0; a < DIM; a++)                                 /*Do not touch */
  {                                                         /*Do not touch */
    mp->d_diffusivity[w][MESH_DISPLACEMENT1 + a] = dDdX[a]; /*Do not touch */
  }                                                         /*Do not touch */
  for (a = 0; a < MAX_CONC; a++)                            /*Do not touch */
  {                                                         /*Do not touch */
    mp->d_diffusivity[w][MAX_VARIABLE_TYPES + a] = dDdC[a]; /*Do not touch */
  }
  /**************************************************************/

  return 0;
} /* End of usr_diffusivity */
/*****************************************************************************/

int usr_FlowingLiquidViscosity(dbl *param) /* ptr to user-defined parameter list */
{
  /* Local Variables */
  int a;
  dbl mu, dmudT;       /* thermal conductivity and its derivative wrt temperature*/
  dbl dmudV[DIM];      /* heat source derivative wrt velocity*/
  dbl dmudC[MAX_CONC]; /* heat source derivative wrt concentration*/
  dbl dmudX[DIM];      /* heat source derivative wrt displacement*/

  dbl X[DIM], F, T, C[MAX_CONC]; /* Convenient local variables */

  /* Begin Execution */

  /**********************************************************/

  /* Comment out our remove this line if using this routine */

  GOMA_EH(GOMA_ERROR, "No user_FlowingLiquidViscosity model implemented");

  /**********************************************************/

  /************Initialize everything for safety**************/
  mu = 0;    /*Do not touch */
  dmudT = 0; /*Do not touch */
  for (a = 0; a < DIM; a++)
    dmudV[a] = 0.; /*Do not touch */
  for (a = 0; a < DIM; a++)
    dmudX[a] = 0.; /*Do not touch */
  for (a = 0; a < MAX_CONC; a++)
    dmudC[a] = 0.; /*Do not touch */
                   /**********************************************************/

  /***********Load up convenient local variables*************/
  /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T; /*Do not touch */
  F = fv->F; /*Do not touch */
  for (a = 0; a < DIM; a++)
    X[a] = fv->x[a]; /*Do not touch */
  for (a = 0; a < pd->Num_Species_Eqn; a++)
    C[a] = fv->c[a]; /*Do not touch */

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
  mp->FlowingLiquid_viscosity = mu;                                   /*Do not touch */
  mp->d_FlowingLiquid_viscosity[TEMPERATURE] = dmudT;                 /*Do not touch */
                                                                      /*Do not touch */
  for (a = 0; a < DIM; a++)                                           /*Do not touch */
  {                                                                   /*Do not touch */
    mp->d_FlowingLiquid_viscosity[VELOCITY1 + a] = dmudV[a];          /*Do not touch */
  }                                                                   /*Do not touch */
  for (a = 0; a < DIM; a++)                                           /*Do not touch */
  {                                                                   /*Do not touch */
    mp->d_FlowingLiquid_viscosity[MESH_DISPLACEMENT1 + a] = dmudX[a]; /*Do not touch */
  }                                                                   /*Do not touch */
  for (a = 0; a < MAX_CONC; a++)                                      /*Do not touch */
  {                                                                   /*Do not touch */
    mp->d_FlowingLiquid_viscosity[MAX_VARIABLE_TYPES + a] = dmudC[a]; /*Do not touch */
  }                                                                   /*Do not touch */

  return (0);
} /* End of usr_FlowingLiquidViscosity */
/*****************************************************************************/

/*
 *  SOLID VISCOSITY
 */
/*
 * int usr_solid_viscosity ()
 *
 * ----------------------------------------------------------------------------
 * This routine is responsible for filling up the following portions of the elc structure
 * at the current gauss point:
 *     intput:  param - array of constants input on the property card.
 *
 *     output:  f       => elc->solid_viscosity
 *              dfdT    => elc->d_viscoslambda[a][TEMPERATURE]
 *                                         - derivative wrt temperature.
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

int usr_solid_viscosity(dbl *param, /* ptr to user-defined parameter list        */
                        double *viscos,
                        double d_viscos_dx[MAX_VARIABLE_TYPES + MAX_CONC]) {
  int a, b;
  int w;

  dbl f, dfdT;        /* momentum sources and its derivative wrt temperature*/
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
  f = 0;    /*Do not touch */
  dfdT = 0; /*Do not touch */
  for (b = 0; b < DIM; b++)
    dfdV[b] = 0.; /*Do not touch */
  for (b = 0; b < DIM; b++)
    dfdX[b] = 0.; /*Do not touch */
  for (w = 0; w < MAX_CONC; w++)
    dfdC[w] = 0.; /*Do not touch */
  /**********************************************************/

  /***********Load up convenient local variables*************/
  /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T; /*Do not touch */
  for (a = 0; a < DIM; a++)
    X[a] = fv->x[a]; /*Do not touch */
  for (i = 0; i < pd->Num_Species_Eqn; i++)
    C[i] = fv->c[i]; /*Do not touch */

  f = param[0] * sin(param[1] * tran->time_value);

  /**********************************************************/
  f = 1.;

  /****************Don't touch these lines***********************/
  for (a = 0; a < DIM; a++) {
    *viscos = f;                                     /*Do not touch */
    d_viscos_dx[TEMPERATURE] = dfdT;                 /*Do not touch */
                                                     /*Do not touch */
    for (b = 0; b < DIM; b++)                        /*Do not touch */
    {                                                /*Do not touch */
      d_viscos_dx[VELOCITY1 + b] = dfdV[b];          /*Do not touch */
    }                                                /*Do not touch */
    for (b = 0; b < DIM; b++)                        /*Do not touch */
    {                                                /*Do not touch */
      d_viscos_dx[MESH_DISPLACEMENT1 + a] = dfdX[b]; /*Do not touch */
    }                                                /*Do not touch */
    for (w = 0; w < MAX_CONC; w++)                   /*Do not touch */
    {                                                /*Do not touch */
      d_viscos_dx[MAX_VARIABLE_TYPES + w] = dfdC[w]; /*Do not touch */
    }                                                /*Do not touch */
  }                                                  /*Do not touch */

  return (0);
} /* End of usr_solid_viscosity */

/*****************************************************************************/
/********************************************
 *    user-defined heat flux model  *
 ********************************************/
#if defined SECOR_HEAT_FLUX
double usr_heat_flux(const double gradP[],      /*   pressure gradient  */
                     double q[],                /*   flow vector            */
                     double dq_gradP[DIM][DIM], /*   flow sens wrt gradP    */
                     double dq_dX[DIM][DIM],    /*   flow sens wrt coords   */
                     const double time,
                     const double h,
                     const double dh_dX[DIM],
                     const double vb[DIM],
                     const double vt[DIM],
                     double dq_dVb[DIM][DIM],
                     double dq_dVt[DIM][DIM])
#else
double usr_heat_flux(const double gradP[],      /*   pressure gradient  */
                     double q[],                /*   flow vector            */
                     double dq_gradP[DIM][DIM], /*   flow sens wrt gradP    */
                     double dq_dX[DIM][DIM],    /*   flow sens wrt coords   */
                     const double time)
#endif
{
  GOMA_EH(GOMA_ERROR, "No usr_heat_flux model supplied");
  return (1.);
} /* End of usr_heat_flux */
/****************************************************************************/

/****************************************************
 *    user-defined porous media permeability model  *
 ****************************************************/
int usr_permeability(dbl *param) /* user-defined parameter list */
{
  dbl k, phi, dkdphi; /* permeability and its derivative wrt porosity*/

  /* Comment out or remove this line if using this routine */

  /*   GOMA_EH(GOMA_ERROR,"No user_defined  permeability model implemented");  */

  /**********************************************************/

  /************Initialize everything for safety**************/
  k = 0.;      /*Do not touch */
  dkdphi = 0.; /*Do not touch */
               /**********************************************************/

  /***********Load up convenient variables*************/

  phi = mp->porosity; /*Do not touch */

  /**********************************************************/

  /*******Add property function and sensitivities here*******/

  if (phi < param[3]) {
    k = param[2];
    dkdphi = 0.0;
  } else {
    k = param[0] + (param[2] - param[0]) / (param[3] - param[1]) * (phi - param[1]);
    dkdphi = (param[2] - param[0]) / (param[3] - param[1]);
  }

  /****************Don't touch these lines***********************/
  mp->permeability = k;                      /*Do not touch */
  mp->d_permeability[POR_POROSITY] = dkdphi; /*Do not touch */
  /**************************************************************/

  return (0);
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
int usr_yield_stress(dbl *param, dbl time) /* pointer to user-defined parameter list    */
{
  /* Local Variables */
  int a;
  dbl tau_y; /* thermal conductivity and its derivative wrt temperature*/
  dbl dtau_ydT;
  dbl dtau_ydV[DIM];      /* heat source derivative wrt velocity*/
  dbl dtau_ydC[MAX_CONC]; /* heat source derivative wrt concentration*/
  dbl dtau_ydX[DIM];      /* heat source derivative wrt displacement*/

  dbl X[DIM], F, T, C[MAX_CONC]; /* Convenient local variables */

  /* Begin Execution */

  /**********************************************************/

  /* Comment out our remove this line if using this routine */

  GOMA_EH(GOMA_ERROR, "No user_defined yield stress model implemented");

  /**********************************************************/

  /************Initialize everything for saftey**************/
  tau_y = 0;    /*Do not touch */
  dtau_ydT = 0; /*Do not touch */
  for (a = 0; a < DIM; a++)
    dtau_ydV[a] = 0.; /*Do not touch */
  for (a = 0; a < DIM; a++)
    dtau_ydX[a] = 0.; /*Do not touch */
  for (a = 0; a < MAX_CONC; a++)
    dtau_ydC[a] = 0.; /*Do not touch */
  /**********************************************************/

  /***********Load up convenient local variables*************/
  /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T; /*Do not touch */
  F = fv->F; /*Do not touch */
  for (a = 0; a < DIM; a++)
    X[a] = fv->x[a]; /*Do not touch */
  for (a = 0; a < pd->Num_Species_Eqn; a++)
    C[a] = fv->c[a]; /*Do not touch */

  /**********************************************************/

  /*******Add property function and sensitivities here*******/
  /* An example:
   *  mu = param[0] * exp(param[1] * (param[2] - C[0]));
   *
   * Add sensitivities here
   *  dmudC[0] = - mu * param[1];
   *
   */

  tau_y = param[0] + (param[1] - param[0]) * exp(-time / param[2]);

  /****************Don't touch these lines***********************/
  gn->tau_y = tau_y; /*Do not touch */
                     /*Do not touch */
  return (0);
} /* End of usr_yield_stress */

// restore warnings
#pragma GCC diagnostic pop

/*****************************************************************************/
/* End of file user_mp.c */
/*****************************************************************************/