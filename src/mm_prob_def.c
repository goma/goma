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
 

#include "mm_prob_def.h"

#include <stdio.h>

#include "std.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "mm_input.h"

#ifndef lint
#endif

#include "el_geom.h"		/* Has info I'd like to replicate into the */
				/* Problem_Description structure... */

#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"

#define GOMA_MM_PROB_DEF_C

/* 
 * The following routine sets up parts of the Problem_Description structure
 * that are used during assembly. This structure is meant to contain information
 * similar to the numerous flags above, but in a handier format that can easily
 * be examined during assembly for a variety of equations and problem types.
 * Assume "pd" defined in one of the include files and allocated earlier....
 * 
 * in:
 * 	(none)
 *
 * out:
 * 	(none)
 * 
 * 
 * Return values: 0 == everything went OK fine
 *               -1 == something went wrong
 * 
 * Notes:	Remember that this routine is after mm_input.
 *
 */

int 
setup_pd(void)
{
  int	i;
  int   mn;                     /* Current material number */
  int   imtrx;                  /* Current matrix number */
  int	ce;			/* Current equation. */
  int   v;			/* Current variable. */
  int status;

  int kkSurf=10000;		/* Legacy variables made ridiculously large */
  int kkBulk=10000;		/* so you'll know if you attempt to use them.
				 * i.e., you should be implementing a better
				 * solution */
  int CoordinateSystem;		/* Used to test for single coordinate
				 * system over multiple materials. */


  static const char yo[] = "setup_pd";

  /*
   * WARNING: This hardwired coordinate flag is not general. - pas 94/03/22
   */
  status = 0.;

  CoordinateSystem = pd_glob[0]->CoordinateSystem;

  for(mn = 0; mn < upd->Num_Mat; mn++)
    {
      pd_glob[mn]->Num_Dim          = Num_Dim; 		/* from "el_geom.h" */
      pd_glob[mn]->TimeIntegration  = TimeIntegration; 	/* from "rf_fem.h" */
      if(pd_glob[mn]->CoordinateSystem != CoordinateSystem)
	EH(-1, "Not all materials have the same coordinate system!");
    }

  if(CoordinateSystem == CYLINDRICAL || CoordinateSystem == SWIRLING)
    {
      if (pd_glob[0]->Num_Dim == 3)
         EH(-1,"Whoa, Whoa.  3D mesh but CYLINDRICAL COORDINATE SYSTEM???");
      VIM = 3;
    }
  else if(CoordinateSystem == PROJECTED_CARTESIAN)
    {
      if(pd_glob[0]->Num_Dim == 3)
	EH(-1, "Achtung!  You cannot combine the PROJECTED_CARTESIAN coordinate system with a 3D mesh.");
      VIM = 3;
    }
  else if(CoordinateSystem == CARTESIAN_2pt5D)
    {
      if(pd_glob[0]->Num_Dim == 3)
	EH(-1, "Whoa!  3D mesh for 2-1/2D Calculation.");
      VIM = 3;
    }
  else
    {
      VIM = Num_Dim;
    }


   /*
    * Make one to one correspondence of each material in input file with each
    * block id in exodusII file, in numerical order.
    *
    * Too restrictive for distributed processing. See routine
    * check_elemblocks in rd_mesh.c for improved version of this intent.
    */

  /*
   * This section of code tries to "OR" together all the Booleans describing
   * the active terms in each equation. It does the ORing based on the floating
   * point values specified for the term multipliers. Those term multipliers
   * that are nonzero cause a flip in the appropriate Boolean variable for
   * this equation...
   */

  for (mn = 0; mn < upd->Num_Mat; mn++)    
    {
      for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++)
         {
          for ( i=0; i<pd_glob[mn]->Num_EQ[imtrx]; i++)
	     {
	      ce = pd_glob[mn]->m[imtrx][i];
	      if((ce == R_PRESSURE)   ||
	         (ce == R_GRADIENT11) ||
	         (ce == R_GRADIENT12) ||
	         (ce == R_GRADIENT13) ||
	         (ce == R_GRADIENT21) ||
	         (ce == R_GRADIENT22) ||
	         (ce == R_GRADIENT23) ||
	         (ce == R_GRADIENT31) ||
	         (ce == R_GRADIENT32) ||
	         (ce == R_GRADIENT33) ||
	         (ce == R_EFIELD1)    ||
	         (ce == R_EFIELD2)    ||
	         (ce == R_EFIELD3)    ||
	         (ce == R_ENORM)      ||
	         (ce == R_EXT_VELOCITY)||
	         (ce == R_VORT_DIR1) ||
	         (ce == R_VORT_DIR2) ||
	         (ce == R_VORT_DIR3) ||
	         (ce == R_NORMAL1)    ||
	         (ce == R_NORMAL2)    ||
	         (ce == R_NORMAL3)    ||
                 (ce == R_SHELL_SHEAR_TOP) ||
                 (ce == R_SHELL_SHEAR_BOT) ||
                 (ce == R_SHELL_CROSS_SHEAR))  
	        {
	         if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_ADVECTION)] != 0. )
		   {
		    pd_glob[mn]->e[imtrx][ce] |= T_ADVECTION;
		   }
	         if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_SOURCE)] != 0. )
		   {
		    pd_glob[mn]->e[imtrx][ce] |= T_SOURCE;
		   }
	        }
	      else if(
		  (ce == R_SHELL_SAT_CLOSED) ||
		  (ce == R_SHELL_DELTAH) ||
		  (ce == R_SHELL_LUB_CURV) ||
		  (ce == R_SHELL_LUB_CURV_2)
		  )
	        {
	         if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_MASS)] != 0. )
		   {
		    pd_glob[mn]->e[imtrx][ce] |= T_MASS;
		   }
	         if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_DIVERGENCE)] != 0. )
		   {
		    pd_glob[mn]->e[imtrx][ce] |= T_DIVERGENCE;
		   }
	        }
	      else if(
		  (ce == R_SHELL_SAT_GASN ||
		   ce == R_MAX_STRAIN ||
                   ce == R_CUR_STRAIN)
		  )
	        {
	         if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_MASS)] != 0. )
		   {
		    pd_glob[mn]->e[imtrx][ce] |= T_MASS;
		   }
	         if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_SOURCE)] != 0. )
		   {
		    pd_glob[mn]->e[imtrx][ce] |= T_SOURCE;
		   }
	        }
	      else if(
		  (ce == R_SHELL_SAT_OPEN) ||
		  (ce == R_SHELL_SAT_OPEN_2)
		  )
	        {
	         if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_MASS)] != 0. )
		   {
		    pd_glob[mn]->e[imtrx][ce] |= T_MASS;
		   }
	         if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_DIFFUSION)] != 0. )
		   {
		    pd_glob[mn]->e[imtrx][ce] |= T_DIFFUSION;
		   }
	         if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_SOURCE)] != 0. )
		   {
		    pd_glob[mn]->e[imtrx][ce] |= T_SOURCE;
		   }
	        }
	      else if (
                   (ce == R_FILL) || 
		   (ce == R_PHASE1) ||
		   (ce == R_PHASE2) ||
		   (ce == R_PHASE3) ||
		   (ce == R_PHASE4) ||
		   (ce == R_PHASE5) ||
		   (ce == R_ACOUS_REYN_STRESS) ||
		   (ce == R_POR_SINK_MASS) ||
                   (ce == R_SHELL_LUBP))
	    {
	      if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_MASS)] != 0. )
		{
		  pd_glob[mn]->e[imtrx][ce] |= T_MASS;
		}
	      if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_ADVECTION)] != 0. )
		{
		  pd_glob[mn]->e[imtrx][ce] |= T_ADVECTION;
		}
	      if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_SOURCE)] != 0. )
		{
		  pd_glob[mn]->e[imtrx][ce] |= T_SOURCE;
		}
	    }
          else if (ce == R_BOND_EVOLUTION ||
                         (ce == R_TFMP_MASS) ||
                         (ce == R_TFMP_BOUND))
	    {
	      if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_MASS)] != 0. )
		{
		  pd_glob[mn]->e[imtrx][ce] |= T_MASS;
		}
	      if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_ADVECTION)] != 0. )
		{
		  pd_glob[mn]->e[imtrx][ce] |= T_ADVECTION;
		}
	      if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_DIFFUSION)] != 0. )
		{
		  pd_glob[mn]->e[imtrx][ce] |= T_DIFFUSION;
		}
	      if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_SOURCE)] != 0. )
		{
		  pd_glob[mn]->e[imtrx][ce] |= T_SOURCE;
		}
	    }
	  else if(ce == R_POTENTIAL)
	    {
	      if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_MASS)] != 0. )   /* mass term added by KSC: 2/4/99 */
		{
		  pd_glob[mn]->e[imtrx][ce] |= T_MASS;
		}
	      if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_BOUNDARY)] != 0. )
		{
		  pd_glob[mn]->e[imtrx][ce] |= T_BOUNDARY;
		}
	      if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_DIFFUSION)] != 0. )
		{
		  pd_glob[mn]->e[imtrx][ce] |= T_DIFFUSION;
		}
	      if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_SOURCE)] != 0. )
		{
		  pd_glob[mn]->e[imtrx][ce] |= T_SOURCE;
		}
	    }
	  else if(ce == R_SHELL_CURVATURE)
	    {
	      if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_DIFFUSION)] != 0. )
		{
		  pd_glob[mn]->e[imtrx][ce] |= T_DIFFUSION;
		}
	    }
	  else if(ce == R_SHELL_CURVATURE2)
	    {
	      if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_DIFFUSION)] != 0. )
		{
		  pd_glob[mn]->e[imtrx][ce] |= T_DIFFUSION;
		}
	    }
	  else if(ce == R_SHELL_TENSION)
	    {
	      if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_DIFFUSION)] != 0. )
		{
		  pd_glob[mn]->e[imtrx][ce] |= T_DIFFUSION;
		}
	    }
	  else if(ce == R_SHELL_X)
	    {
	      if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_DIFFUSION)] != 0. )
		{
		  pd_glob[mn]->e[imtrx][ce] |= T_DIFFUSION;
		}
	    }
	  else if(ce == R_SHELL_Y)
	    {
	      if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_DIFFUSION)] != 0. )
		{
		  pd_glob[mn]->e[imtrx][ce] |= T_DIFFUSION;
		}
	    }
          else if(ce == R_SHELL_DIFF_FLUX)
            {
                 if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_DIFFUSION)] != 0. )
                   {
                    pd_glob[mn]->e[imtrx][ce] |= T_DIFFUSION;
                   }
                }
              else if(ce == R_SHELL_DIFF_CURVATURE)
                {
                 if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_DIFFUSION)] != 0. )
                   {
                    pd_glob[mn]->e[imtrx][ce] |= T_DIFFUSION;
                   }
                }
              else if(ce == R_SHELL_NORMAL1)
                {
                 if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_DIFFUSION)] != 0. )
                   {
                    pd_glob[mn]->e[imtrx][ce] |= T_DIFFUSION;
                   }
                }
              else if(ce == R_SHELL_NORMAL2)
                {
                 if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_DIFFUSION)] != 0. )
                   {
                    pd_glob[mn]->e[imtrx][ce] |= T_DIFFUSION;
                   }
                }
                else if(ce == R_SHELL_NORMAL3)
                {
                  if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_DIFFUSION)] != 0. )
                  {
                    pd_glob[mn]->e[imtrx][ce] |= T_DIFFUSION;
                  }
                }
	      else if(ce == R_SHEAR_RATE )
	        {
	         if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_ADVECTION)] != 0. )
		   {
		    pd_glob[mn]->e[imtrx][ce] |= T_ADVECTION;
		   }
	         if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_DIFFUSION)] != 0. )
		   {
		    pd_glob[mn]->e[imtrx][ce] |= T_DIFFUSION;
		   }
	         if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_SOURCE)] != 0. )
		   {
		    pd_glob[mn]->e[imtrx][ce] |= T_SOURCE;
		   }
	        }
	      else if(
                  (ce == R_CURVATURE ) ||
		  (ce == R_LUBP) ||
		  (ce == R_LUBP_2))
	        {
	         if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_BOUNDARY)] != 0. )
		   {
		    pd_glob[mn]->e[imtrx][ce] |= T_BOUNDARY;
		   }
	         if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_DIFFUSION)] != 0. )
		   {
		    pd_glob[mn]->e[imtrx][ce] |= T_DIFFUSION;
		   }
	         if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_SOURCE)] != 0. )
		   {
		    pd_glob[mn]->e[imtrx][ce] |= T_SOURCE;
		   }
	        }
	      else if(
                  (ce == R_ENERGY )||
		  (ce == R_MASS)||
		  (ce == R_MESH1)||
		  (ce == R_MESH2)||
		  (ce == R_MESH3)||
		  (ce == R_SOLID1)||
		  (ce == R_SOLID2)||
		  (ce == R_SOLID3)||
		  (ce == R_POR_LIQ_PRES)||
		  (ce == R_POR_GAS_PRES)||
		  (ce == R_POR_POROSITY)||
		  (ce == R_POR_ENERGY)||
		  (ce == R_POR_SATURATION) ||
                  (ce == R_SURF_CHARGE)  ||
 		  (ce == R_SHELL_BDYVELO)    ||
 		  (ce == R_SHELL_USER)    ||
		  (ce == R_ACOUS_PREAL)  ||
		  (ce == R_ACOUS_PIMAG)   ||
		  (ce == R_LIGHT_INTP)   ||
		  (ce == R_LIGHT_INTM)   ||
		  (ce == R_LIGHT_INTD)   ||
		  (ce == R_RESTIME)   ||  
		  (ce == R_EM_E1_REAL)  ||
		  (ce == R_EM_E2_REAL)  ||
		  (ce == R_EM_E3_REAL)  ||
		  (ce == R_EM_E1_IMAG)  ||
		  (ce == R_EM_E2_IMAG)  ||
		  (ce == R_EM_E3_IMAG)  ||
		  (ce == R_EM_H1_REAL)  ||
		  (ce == R_EM_H2_REAL)  ||
		  (ce == R_EM_H3_REAL)  ||
		  (ce == R_EM_H1_IMAG)  ||
		  (ce == R_EM_H2_IMAG)  ||
		  (ce == R_EM_H3_IMAG)  ||
		  (ce == R_SHELL_FILMP) ||
                  (ce == R_SHELL_FILMH) ||
                  (ce == R_SHELL_PARTC) || 
		  (ce == R_SHELL_ENERGY))
	        {
	         if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_MASS)] != 0. )
		   {
		    pd_glob[mn]->e[imtrx][ce] |= T_MASS;
		   }
	         if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_ADVECTION)] != 0. )
		   {
		    pd_glob[mn]->e[imtrx][ce] |= T_ADVECTION;
		   }
	         if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_BOUNDARY)] != 0. )
		   {
		    pd_glob[mn]->e[imtrx][ce] |= T_BOUNDARY;
		   }
	         if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_DIFFUSION)] != 0. )
		   {
		    pd_glob[mn]->e[imtrx][ce] |= T_DIFFUSION;
		   }
	         if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_SOURCE)] != 0. )
		   {
		    pd_glob[mn]->e[imtrx][ce] |= T_SOURCE;
		   }
	        }
	      else if(
                  (ce == R_MOMENTUM1)||
		  (ce == R_MOMENTUM2)||
		  (ce == R_MOMENTUM3)||
		  (ce == R_PMOMENTUM1)||
		  (ce == R_PMOMENTUM2)||
		  (ce == R_PMOMENTUM3)||
		  (ce == R_STRESS11)||
		  (ce == R_STRESS12)||
		  (ce == R_STRESS13)||
		  (ce == R_STRESS22)||
		  (ce == R_STRESS23)||
		  (ce == R_STRESS33)||

		  (ce == R_STRESS11_1)||
		  (ce == R_STRESS12_1)||
		  (ce == R_STRESS13_1)||
		  (ce == R_STRESS22_1)||
		  (ce == R_STRESS23_1)||
		  (ce == R_STRESS33_1)||		  

		  (ce == R_STRESS11_2)||
		  (ce == R_STRESS12_2)||
		  (ce == R_STRESS13_2)||
		  (ce == R_STRESS22_2)||
		  (ce == R_STRESS23_2)||
		  (ce == R_STRESS33_2)||

		  (ce == R_STRESS11_3)||
		  (ce == R_STRESS12_3)||
		  (ce == R_STRESS13_3)||
		  (ce == R_STRESS22_3)||
		  (ce == R_STRESS23_3)||
		  (ce == R_STRESS33_3)||

		  (ce == R_STRESS11_4)||
		  (ce == R_STRESS12_4)||
		  (ce == R_STRESS13_4)||
		  (ce == R_STRESS22_4)||
		  (ce == R_STRESS23_4)||
		  (ce == R_STRESS33_4)||

		  (ce == R_STRESS11_5)||
		  (ce == R_STRESS12_5)||
		  (ce == R_STRESS13_5)||
		  (ce == R_STRESS22_5)||
		  (ce == R_STRESS23_5)||
		  (ce == R_STRESS33_5)||

		  (ce == R_STRESS11_6)||
		  (ce == R_STRESS12_6)||
		  (ce == R_STRESS13_6)||
		  (ce == R_STRESS22_6)||
		  (ce == R_STRESS23_6)||
		  (ce == R_STRESS33_6)||

		  (ce == R_STRESS11_7)||
		  (ce == R_STRESS12_7)||
		  (ce == R_STRESS13_7)||
		  (ce == R_STRESS22_7)||
		  (ce == R_STRESS23_7)||
		  (ce == R_STRESS33_7))
	        {
	         if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_MASS)] != 0. )
		   {
		    pd_glob[mn]->e[imtrx][ce] |= T_MASS;
		   }
	         if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_ADVECTION)] != 0. )
		   {
		    pd_glob[mn]->e[imtrx][ce] |= T_ADVECTION;
		   }
                 if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_BOUNDARY)] != 0. )
		   {
		    pd_glob[mn]->e[imtrx][ce] |= T_BOUNDARY;
		   }
	         if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_DIFFUSION)] != 0. )
		   {
		    pd_glob[mn]->e[imtrx][ce] |= T_DIFFUSION;
		   }
	         if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_SOURCE)] != 0. )
		   {
		    pd_glob[mn]->e[imtrx][ce] |= T_SOURCE;
		   }
	         if ( pd_glob[mn]->etm[imtrx][ce][(LOG2_POROUS_BRINK)] != 0. )
		   {
		    pd_glob[mn]->e[imtrx][ce] |= T_POROUS_BRINK;
		   }
	        }
	      else if(
                  (ce == R_VORT_LAMBDA) ||
		  (ce == R_LAGR_MULT1)  ||
		  (ce == R_LAGR_MULT2)  ||
 		  (ce == R_LAGR_MULT3)  ||
		  (ce == R_SHELL_ANGLE1)  ||
		  (ce == R_SHELL_ANGLE2)  ||
		  (ce == R_SHELL_SURF_DIV_V) ||
		  (ce == R_SHELL_SURF_CURV) ||
		  (ce == R_N_DOT_CURL_V) ||
		  (ce == R_GRAD_S_V_DOT_N1) ||
		  (ce == R_GRAD_S_V_DOT_N2) ||
		  (ce == R_GRAD_S_V_DOT_N3) ||
		  (ce == R_DENSITY_EQN))
	        {
	      /* These equations have no term multipliers, but
	       * something needs to be set to make the equation
	       * "active". */
	         pd_glob[mn]->e[imtrx][ce] |= T_MASS;
	        }
	      else
	        {
	         fprintf(stderr, "%s: problem with unknown equation\n", yo);
	         status=-1;
	        }
             }

          if(pd_glob[mn]->e[imtrx][R_MOMENTUM3])
            {
	     if(CoordinateSystem == CARTESIAN && Num_Dim == 2)
	     EH(-1, "You have 3 velocity components, but only a 2D CARTESIAN mesh.\nDid you mean to use the PROJECTED_CARTESIAN coordinate system?");
            }
         }
    }
  /* 
   * Now define the number of variables of each type in the problem
   * (i.e. the number of species of Y and S. along with U1,U2,U3,T,P 
   *
   * Here this is a worst-case approach, used mainly for the post processing
   * We will have to revist this if we want to go on with parallel processing
   */

  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++)
     {
      for (v=V_FIRST; v<V_LAST; v++)
         {
          Num_Var_In_Type[imtrx][v] = 0;
         }
     }
  for( mn = 0; mn < upd->Num_Mat; mn++)
     {
    
      /* 
       * Assign local pointer pd to appropriate material
       */
      pd = pd_glob[mn];

      for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++)
         {

          for (v=V_FIRST; v<V_LAST; v++)
	     {
	      if(v == MASS_FRACTION)
	        {
	         Num_Var_In_Type[imtrx][MASS_FRACTION]      |= ( pd->e[imtrx][R_MASS] )      ? upd->Max_Num_Species_Eqn: 0;
	        }
	      else if (v == SURFACE )
	        {
	         Num_Var_In_Type[imtrx][SURFACE]            |= ( pd->e[imtrx][R_MASS_SURF])  ?  (kkSurf + kkBulk): 0;
	        }
	      else
	        {
	         Num_Var_In_Type[imtrx][v]          |= ( pd->e[imtrx][v] ) ?         1         : 0;
	        }
	     }

         }
     }


  return(status);

} /* end of routine setup_pd() */

/* end of file mm_prob_def.c */
