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
 *$Id: user_bc.c,v 5.1 2010-04-05 15:05:19 hkmoffa Exp $
 */

/* Standard include files */

#ifdef USE_RCSID
static char rcsid[] =
"$Id: user_bc.c,v 5.1 2010-04-05 15:05:19 hkmoffa Exp $";
#endif

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

#include "goma.h"

/*
 * Prototype declarations of functions defined in this file.
 */

#ifdef EXTERN
#undef EXTERN
#endif

#define EXTERN /* nothing */

#include "user_bc.h"

#ifdef EXTERN
#undef EXTERN
#endif

#define EXTERN extern

/*
 * Prototype declarations of functions not defined in this file.
 */

/*
 * Function definitions.
 */

/*****************************************************************************/
/*       Functions for user_defined velocity profiles                        */
/*****************************************************************************/

dbl velo_vary_fnc( velo_condition, x1, x2, x3, p, time)
     const int velo_condition;
     const dbl x1;
     const dbl x2;
     const dbl x3;
     const dbl p[];
     const dbl time;
{
  double f = 0.0;
  int model_id;
	  model_id = ((int)p[0]);
  /*  double a1,a2;

  a1 = p[0];
  a2 = p[1];
  */

  /*  radius = sqrt( x1*x1 + x2*x2 ); 
      
      if( velo_condition == UVARY_BC) 
      { 
      f = -a2*x2/radius; 
      } 
      else if ( velo_condition == VVARY_BC ) 
      { 
      
      f = a1*x1/radius; 
      }  */
  
  
  if( velo_condition == UVARY_BC)
    {
	if(model_id == 1)
	{
/*    parabolic velocity profile
 *      p[1] = coordinate1
 *      p[2] = coordinate2
 *      p[3] = flow in positive coordinate direction
 */
double coord1, coord2, qflow, gap;
	coord1 = p[1];
	coord2 = p[2];
	qflow = p[3];
	gap = fabs(coord2-coord1);
	if(gap != 0.0)	{
		f = 6.*qflow*(x2-coord1)*(coord2-x2)/(gap*gap*gap);
		}	else	{
		f=0.0;
		}
	}
	else if(model_id == 2)
	{
/*    parabolic velocity profile
 *      p[1] = coordinate1
 *      p[2] = coordinate2
 *      p[3] = flow in positive coordinate direction
 *	p[4] = velocity1
 *	p[5] = velocity2
 */
double coord1, coord2, qflow, gap, veloc1, veloc2;
	coord1 =  p[1];
	coord2 = p[2];
	qflow = p[3];
	veloc1 = p[4];
	veloc2 = p[5];

	gap = fabs(coord2-coord1);
	if(gap != 0.0)	{
		f = (6.*qflow-3.*gap*(veloc1+veloc2))
				*(x2-coord1)*(coord2-x2)/(gap*gap*gap)
				+ (veloc1*(coord2-x2)+veloc2*(x2-coord1))/gap;
		}	else	{
		f=0.5*(veloc1+veloc2);
		}
	}
	else if(model_id == 3)
	{
/*    upstream profile of forward roll coating nip
 *      p[1] = web_speed_lower
 *      p[2] = roll_radius_lower
 *      p[3] = web_speed_upper
 *	p[4] = roll_radius_upper
 *	p[5] = roll_separation
 *	p[6] = flowrate
 */
double coord1, coord2, qflow, gap, veloc1, veloc2;
	coord1 =  -p[2]+sqrt(SQUARE(p[2])-SQUARE(x1));
	coord2 = p[5]+p[4]-sqrt(SQUARE(p[4])-SQUARE(x1));
	qflow = p[6];
	veloc1 = p[1]*sqrt(1.0-SQUARE(x1/p[2]));
	veloc2 = p[3]*sqrt(1.0-SQUARE(x1/p[4]));

	gap = fabs(coord2-coord1);
	if(gap != 0.0)	{
		f = (6.*qflow-3.*gap*(veloc1+veloc2))
				*(x2-coord1)*(coord2-x2)/(gap*gap*gap)
				+ (veloc1*(coord2-x2)+veloc2*(x2-coord1))/gap;
		}	else	{
		f=0.5*(veloc1+veloc2);
		}
	}
	/* micro-flexo-printing bottom surface - omega, roll_rad, gap */
	else
	if (model_id == 15) {
	    
	  double omega, roll_rad, gap, t_offset, time1, yc;
	  omega = p[1];
	  roll_rad = p[2];
	  gap = p[3];
	  t_offset = p[4];

	  time1 = time + t_offset;
	  yc = gap + 2.*roll_rad*(1.-sin(omega*time1));
	  f = -omega*(x2 - yc);
		if(nAC > 0)	{
			f += augc[0].tmp2*cos(omega*time1);
			}
	  }
	/* micro-printing (axi-symmetric)	*/
	else 
	if (model_id == 19)
	{
	  dbl tdc_delay;
	  double omega, roll_rad, t_offset, time1;
	  omega = p[1];
	  roll_rad = p[2];
	  t_offset = p[5];
	  tdc_delay=p[8];

	  time1 = time + t_offset;
	  if (time <= 0.0)	
		{ f = 2.*roll_rad*(-omega*cos(omega*time1));}
	  else if(time <= tdc_delay)
		{ f = 0.0;}
	  else 
		{ time1 -= tdc_delay;
		 f = 2.*roll_rad*(-omega*cos(omega*time1));
		}
	}   /*  end of micro-printing (axi-symmetric) */
	/* micro-flexo-printing bottom surface - general */
	else
	if (model_id == 22) {
	    
	  double omega_l, omega_u, roll_rad_l, roll_rad_u, t_offset, gap;
	  double time1, angle_l, angle_u, yc;
	  omega_l = p[1];
	  omega_u = p[2];
	  roll_rad_l = p[3];
	  roll_rad_u = p[4];
	  t_offset = p[5];
	  gap = p[6];

	  time1 = time + t_offset;
	  angle_l = omega_l*time1;
	  angle_u = omega_u*time1 + 0.5*M_PIE*(1.-omega_u/omega_l);
	  yc = gap + roll_rad_l*(1.-sin(angle_l)) + roll_rad_u*(1.-sin(angle_u));
	  f = -omega_u*(x2 - yc);
		if(nAC > 0)	{
		     f += augc[0].tmp2*cos(angle_u) - augc[0].tmp1*omega_u*sin(angle_u);
		     }
	  }
	else
	{
	EH(-1,"invalid model_id in velo_user_bc");
	}
    }
  else if ( velo_condition == VVARY_BC )
    {
    	if(model_id == 1)
	{
/*    parabolic velocity profile
 *      p[1] = coordinate1
 *      p[2] = coordinate2
 *      p[3] = flow in positive coordinate direction
 */
double coord1, coord2, qflow,gap;
	coord1 = p[1];
	coord2 = p[2];
	qflow = p[3];
	gap = fabs(coord2-coord1);
	if(gap != 0.0)	{
		f = 6.*qflow*(x1-coord1)*(coord2-x1)/(gap*gap*gap);
		}	else	{
		f=0.0;
		}
	}
	else if(model_id == 2)
	{
/*    parabolic velocity profile
 *      p[1] = coordinate1
 *      p[2] = coordinate2
 *      p[3] = flow in positive coordinate direction
 *	p[4] = velocity1
 *	p[5] = velocity2
 */
double coord1, coord2, qflow, gap, veloc1, veloc2;
	coord1 = p[1];
	coord2 = p[2];
	qflow = p[3];
	veloc1 = p[4];
	veloc2 = p[5];
	gap = fabs(coord2-coord1);
	if(gap != 0.0)	{
		f = (6.*qflow-3.*gap*(veloc1+veloc2))
				*(x1-coord1)*(coord2-x1)/(gap*gap*gap)
				+ (veloc1*(coord2-x1)+veloc2*(x1-coord1))/gap;
		}	else	{
		f=0.5*(veloc1+veloc2);
		}
	}
	else if(model_id == 3)
	{
/*    upstream profile of forward roll coating nip
 *      p[1] = web_speed_lower
 *      p[2] = roll_radius_lower
 *      p[3] = web_speed_upper
 *	p[4] = roll_radius_upper
 *	p[5] = roll_separation
 *	p[6] = flowrate
 */
double coord1, coord2, gap, veloc1, veloc2;
	coord1 =  -p[2]+sqrt(SQUARE(p[2])-SQUARE(x1));
	coord2 = p[5]+p[4]-sqrt(SQUARE(p[4])-SQUARE(x1));
	veloc1 = -p[1]*x1/p[2];
	veloc2 = p[3]*x1/p[4];

	gap = fabs(coord2-coord1);
	if(gap != 0.0)	{
		f = (veloc1*(coord2-x2)+veloc2*(x2-coord1))/gap;
		}	else	{
		f=0.5*(veloc1+veloc2);
		}
	}
	/* micro-flexo-printing bottom surface - omega, roll_rad, gap */
	else
	if (model_id == 15) {
	    
	  double omega, roll_rad, t_offset, time1, ycdot;
	  omega = p[1];
	  roll_rad = p[2];
	  t_offset = p[4];

	  time1 = time + t_offset;
	  ycdot = -2.*omega*roll_rad*cos(omega*time1);
	  f = ycdot + omega*x1;
/* addition component due to feature compression	*/
		if(nAC > 0)	{
			f += augc[0].tmp2*sin(omega*time1);
			}
	}
	/* micro-printing (axi-symmetric)	*/
	else 
	if (model_id == 19)
	{
	  f = 0.0;
	}   /*  end of micro-printing (axi-symmetric) */
	/* micro-flexo-printing bottom surface - general */
	else
	if (model_id == 22) {
	    
	  double omega_l, omega_u, roll_rad_l, roll_rad_u, t_offset;
	  double time1, ycdot, angle_l, angle_u, xc;
	  omega_l = p[1];
	  omega_u = p[2];
	  roll_rad_l = p[3];
	  roll_rad_u = p[4];
	  t_offset = p[5];

	  time1 = time + t_offset;
	  angle_l = omega_l*time1;
	  angle_u = omega_u*time1 + 0.5*M_PIE*(1.-omega_u/omega_l);
	  xc = -roll_rad_u*cos(angle_u) + roll_rad_l*cos(angle_l);
	  ycdot = -(omega_l*roll_rad_l*cos(angle_l) + 
			omega_u*roll_rad_u*cos(angle_u));
	  f = ycdot + omega_u*(x1 - xc);
/* addition component due to feature compression	*/
		if(nAC > 0)	{
		    f += augc[0].tmp2*sin(angle_u) + augc[0].tmp1*omega_u*cos(angle_u);
		    }
	}
	else
	{
	EH(-1,"invalid model_id in velo_user_bc");
	}
    }
  else if ( velo_condition == WVARY_BC )
    {
    f=0.;
    }
  


  return(f);

/*  return((dbl)0);	*/	/* Here's a good default behavior! */

  /* Example code fragments:
   *
   * channel Poisuelle flow in negative y-direction
   * with p[0]=(flowrate/width) and p[1]=h/2
   *
   *  if ( velo_condition == UVARY_BC )
   *    {
   *       f = (0.1 - 0.3*(x2/p[0])*(x2/p[0]))*p[1]; 
   *
   *     for fountain flow: KSC on 8/1/94 
   *
   *  f = (1.0-(x2/p[0]))*p[1];  KSC on 8/17/94  for blade coating
   *
   *
   * For Couette Poisuelle flow across a gap p[0] thick 
   * with fractional flowrate Q=p[1] V p[0] (Couette if p[1]=0.5)
   * and  v=0 at y=p[0] and v=1 at y=0 
   * (change y/p[0] to 1-y/p[0] to invert this) RAC 10/17/94
   *
   * f=1.0 + (6.0*p[1]-4.0)*((x2 - 2)/p[0])
   *       + (3.0-6.0*p[1])*((x2 - 2)/p[0])*((x2 - 2)/p[0]);
   *
   * printf("in UVARY, x2=, p[0]=, p[1]=, f=%13.5f %13.5f %13.5f %13.5\n",  
   *	  x2,  p[0],  p[1],  f);
   *
   *    f = + 1.5*(1. - SQUARE(x2) ) -1.;
   *     f = 0.2 - 0.4 * SQUARE(x2);     
   *     f = 0.1 - 0.3 * SQUARE(x2);      
   *    f = +1.5*(1. - SQUARE(x2) );
   *
   * Inlet of the backward facing step:
   *      u(y) = 0    at y=step
   *             v_in at y=(step+channel)/2
   *	       0    at y= channel
   *   
   *	v_max   = p[0];
   *	gap     = p[1];
   *	channel = 1;
   *	midpt   = channel - gap/2.;
   *	u       = v_max - (SQUARE(x2-midpt))/(v_max*(SQUARE(gap/2.)));
   *	f       = u;
   *      }
   *  else if (velo_condition == VVARY_BC)
   *    { 
   *      f = -p[1] * (1 - (x1/p[0]) * (x1/p[0]));
   *    }
   *  else if (velo_condition == WVARY_BC)
   *    {
   *      f = 0.;
   *    }
   *  return (f);
   */
}
/*****************************************************************************/
dbl dvelo_vary_fnc_d1( velo_condition, x1, x2, x3, p, time)
     const int velo_condition;
     const dbl x1;
     const dbl x2;
     const dbl x3;
     const dbl p[];
     const dbl time;
{
  dbl f = 0.0;
  int model_id;
  model_id = ((int)p[0]);


  if( velo_condition == UVARY_BC)
     {
	if(model_id == 1)
	{
	f = 0.;
	}
	else
	if(model_id == 2)
	{
	f = 0.;
	}
	else if(model_id == 3)
	{
double coord1, coord2, qflow, gap, veloc1, veloc2;
double coord1dx, coord2dx, gapdx, veloc1dx, veloc2dx;
double tmp1, tmp2, tmp3;
	coord1 =  -p[2]+sqrt(SQUARE(p[2])-SQUARE(x1));
	coord2 = p[5]+p[4]-sqrt(SQUARE(p[4])-SQUARE(x1));
	qflow = p[6];
	veloc1 = p[1]*sqrt(1.0-SQUARE(x1/p[2]));
	veloc2 = p[3]*sqrt(1.0-SQUARE(x1/p[4]));
	coord1dx=-x1/sqrt(SQUARE(p[2])-SQUARE(x1));
	coord2dx=x1/sqrt(SQUARE(p[4])-SQUARE(x1));
	gapdx = coord2dx - coord1dx;
	veloc1dx = -p[1]*x1/sqrt(1.0-SQUARE(x1/p[2]));
	veloc2dx = -p[3]*x1/sqrt(1.0-SQUARE(x1/p[4]));
	gap = fabs(coord2-coord1);
	tmp1 = 2.*qflow - gap*(veloc1+veloc2);
	tmp2 = 3.*qflow - gap*(veloc1+veloc2);
	tmp3 = veloc1*(coord2-x2)+veloc2*(x2-coord1);

	if(gap != 0.0)	{
		f = (x2-coord1)*(coord2-x2)*(-6.*gapdx*tmp2/SQUARE(gap*gap)
			-3.*(veloc1dx+veloc2dx)/SQUARE(gap)) 
		+3.*tmp1*(coord2dx*(x2-coord1)-coord1dx*(coord2-x2))/(gap*gap*gap)
			-gapdx*tmp3/SQUARE(gap)+ (veloc1*coord2dx-veloc2*coord1dx
			+(coord2-x2)*veloc1dx+(x2-coord1)*veloc2dx)/gap;
		}	else	{
		f=0.5*(veloc1dx+veloc2dx);
		}
	}
	else
	{
	  f = 0;
	}
     }
  else if ( velo_condition == VVARY_BC )
     {
        if(model_id == 1)
        {
double coord1, coord2, qflow, gap;
	coord1 = p[1];
	coord2 = p[2];
	qflow = p[3];
	gap = fabs(coord2-coord1);
	if(gap != 0.0)	{
		f = 6.*qflow*(coord1+coord2 - 2.*x1)/(gap*gap*gap);
		}	else	{
		f=0.0;
		}
	}
	else
        if(model_id == 2)
        {
double coord1, coord2, qflow, gap, veloc1, veloc2;
	coord1 = p[1];
	coord2 = p[2];
	qflow = p[3];
	veloc1 = p[4];
	veloc2 = p[5];
	gap = fabs(coord2-coord1);
	if(gap != 0.0)	{
		f = (6.*qflow-3.*gap*(veloc1+veloc2))
			*(coord1+coord2 - 2.*x1)/(gap*gap*gap)
			+ (veloc2-veloc1)/gap;
		}	else	{
		f=0.0;
		}
	}
	else if(model_id == 3)
	{
double coord1, coord2, gap, veloc1, veloc2;
double coord1dx, coord2dx, gapdx, veloc1dx, veloc2dx, tmp3;
	coord1 =  -p[2]+sqrt(SQUARE(p[2])-SQUARE(x1));
	coord2 = p[5]+p[4]-sqrt(SQUARE(p[4])-SQUARE(x1));
	veloc1 = -p[1]*x1/p[2];
	veloc2 = p[3]*x1/p[4];
	coord1dx=-x1/sqrt(SQUARE(p[2])-SQUARE(x1));
	coord2dx=x1/sqrt(SQUARE(p[4])-SQUARE(x1));
	gapdx = coord2dx - coord1dx;
	veloc1dx = -p[1]/p[2];
	veloc2dx = p[3]/p[4];
	tmp3 = veloc1*(coord2-x2)+veloc2*(x2-coord1);

	gap = fabs(coord2-coord1);
	if(gap != 0.0)	{
		f = -gapdx*tmp3/SQUARE(gap) + 
			(veloc1*coord2dx - veloc2*coord1dx + 
			(coord2-x2)*veloc1dx + (x2-coord1)*veloc2dx)/gap;
		f = (veloc1*(coord2-x2)+veloc2*(x2-coord1))/gap;
		}	else	{
		f=0.5*(veloc1dx+veloc2dx);
		}
	}
	/* micro-flexo-printing bottom surface - omega, roll_rad, gap */
	else
	if (model_id == 15) {
	    
	  double omega;
	  omega = p[1];

	  f = omega;
	  }
	/* micro-flexo-printing bottom surface - general */
	else
	if (model_id == 22) {
	    
	  double omega_u;
	  omega_u = p[2];

	  f = omega_u;
	  }
	else
	{
	 f = 0.;
	}
      }
  else
    {
      f = 0.;
    }


  return(f);		/* Here's a good default behavior! */

  /* Example code fragments:
   *  if(velo_condition == UVARY_BC) {
   *    f = 0.; 
   *  } else if (velo_condition == VVARY_BC) { 
   *    f = 2 * p[1] * x1/p[0]/p[0];
   *  } else if (velo_condition == WVARY_BC) {
   *    f = 0;
   *  }
   *  return f;
   */
}
/*****************************************************************************/
dbl dvelo_vary_fnc_d2( velo_condition, x1, x2, x3, p, time)
     const int velo_condition;
     const dbl x1;
     const dbl x2;
     const dbl x3;
     const dbl p[];
     const dbl time;
{
/* for PI use M_PIE Constant from std.h include file. */

  dbl f=0.0;
  int model_id;
  model_id = ((int)p[0]);

if(af->Assemble_LSA_Mass_Matrix)  return 0 ;

  if( velo_condition == UVARY_BC)
	{
		if(model_id == 1)
		{
double coord1, coord2, qflow, gap;
		coord1 = p[1];
		coord2 = p[2];
		qflow = p[3];
		gap = fabs(coord2-coord1);
		if(gap != 0.0)	{
			f = 6.*qflow*(coord1+coord2-2.*x2)/(gap*gap*gap);
			}	else	{
			f=0.0;
			}
		}
 		else
		if(model_id == 2)
		{
double coord1, coord2, qflow, gap, veloc1, veloc2;
		coord1 = p[1];
		coord2 = p[2];
		qflow = p[3];
		veloc1 = p[4];
		veloc2 = p[5];
		gap = fabs(coord2-coord1);
		if(gap != 0.0)	{
			f = (6.*qflow-3.*gap*(veloc1+veloc2))
				*(coord1+coord2-2.*x2)/(gap*gap*gap)
				+ (veloc2-veloc1)/gap;
			}	else	{
			f=0.0;
			}
		}
		else if(model_id == 3)
		{
double coord1, coord2, qflow, gap, veloc1, veloc2;
	coord1 =  -p[2]+sqrt(SQUARE(p[2])-SQUARE(x1));
	coord2 = p[5]+p[4]-sqrt(SQUARE(p[4])-SQUARE(x1));
	qflow = p[6];
	veloc1 = p[1]*sqrt(1.0-SQUARE(x1/p[2]));
	veloc2 = p[3]*sqrt(1.0-SQUARE(x1/p[4]));

	gap = fabs(coord2-coord1);
	if(gap != 0.0)	{
		f = (6.*qflow-3.*gap*(veloc1+veloc2))
				*(coord1+coord2-2.*x2)/(gap*gap*gap)
				+ (veloc2-veloc1)/gap;
		}	else	{
		f=0.5*(veloc1+veloc2);
		}
	}
	/* micro-flexo-printing bottom surface - omega, roll_rad, gap */
 		else
 		if (model_id == 15) {
 	    
 	  	double omega;
 	  	omega = p[1];
 
 	  	f = -omega;
 	  	}
	/* micro-flexo-printing bottom surface - general */
 		else
 		if (model_id == 22) {
 	    
 	  	double omega_u;
 	  	omega_u = p[2];
 
 	  	f = -omega_u;
 	  	}
		else
	 	{
	  	f = 0.;
	 	}
	 }
  else if ( velo_condition == VVARY_BC )
	 {
		if(model_id == 3)
		{
		double coord1, coord2, gap, veloc1, veloc2;
		coord1 =  -p[2]+sqrt(SQUARE(p[2])-SQUARE(x1));
		coord2 = p[5]+p[4]-sqrt(SQUARE(p[4])-SQUARE(x1));
		veloc1 = -p[1]*x1/p[2];
		veloc2 = p[3]*x1/p[4];

		gap = fabs(coord2-coord1);
			if(gap != 0.0)	{
				f = (veloc2-veloc1)/gap;
				}	else	{
				f=0.0;
				}
		}
		else
	 	{
	  	f = 0.;
	 	}
	}
  else
    {
      f = 0.;
    }
  

  return(f);
/*  return((dbl)0);	*/	/* Here's a good default behavior! */

  /* Example code fragments:
   *
   *  if(velo_condition == UVARY_BC) {
   *  f = -0.3*2.0*(x2/p[0])*(1.0/p[0])*p[1];   for fountain flow: KSC on 8/1/94 
   *  f = (-1.0/p[0])*p[1];  for blade coating: KSC on 8/17/94 
   *
   * For Couette Poisuelle flow across a gap p[0] thick 
   *       with fractional flowrate Q=p[1] V p[0] (Couette if p[1]=0.5)
   *       and  v=1 at y=p[0] and v=1 at y=1 
   *       (change y/p[0] to 1-y/p[0] to invert this) RAC 10/17/94
   *	f=(6.0*p[1]-4.0)/p[0]+2.*(3.0-6.0*p[1])/p[0]*((x2-2)/p[0]);
   *
   *   f = -2.*0.4* x2; 
   *    f = -2.*0.3* x2; 
   *    f = -1.5*2.*x2; 
   *  } else if (velo_condition == VVARY_BC) { 
   *    f = 0.; 
   *  } else if (velo_condition == WVARY_BC) {
   *    f = 0;
   *  }
   *  return f;
   */
}

/*****************************************************************************/
dbl dvelo_vary_fnc_d3( velo_condition, x1, x2, x3, p, time)
     const int velo_condition;
     const dbl x1;
     const dbl x2;
     const dbl x3;
     const dbl p[];
     const dbl time;
{
/* for PI use M_PIE Constant from std.h include file. */
  dbl f;

  f = 0.;               /* enter your function here */

  return(f);		/* Here's a good default behavior! */

  /* Example code fragments:
   *
   *  if(velo_condition == UVARY_BC) {
   *    f = 0.; 
   *  } else if (velo_condition == VVARY_BC) { 
   *    f = 0.; 
   *  } else if (velo_condition == WVARY_BC) {
   *    f = 0;
   *  }
   *  return f;
   */
}

/*****************************************************************************/
/*       Functions for solid boundary description                            */
/*****************************************************************************/
/*
 *	user geometry model selection
 *
 *	1  - Gaussian bump on a roll
 *	2  - micro -replicated, gravure roll surface
 *	3  - radiused die lip corner in quadrant I&II
 *      4  - fitted function (must supply pade coefficients)
 *      5  - cell with hyperbolic tan functions
 *      6  - tapered surface for 3D manifold cavities
 *      7  - curve-driven cylindrical surface
 *      8  - single V cell 
 *      9  - Doug's bump
 *      10 - normal line to radiused die lip
 *      11 - 2 radiused die lip corners in quadrant I & II 
 *      12 - circular arc of changing radius
 *      13 - radiused die corner in combination with radiused die
 *      14 - needle collar with tanh function and planes
 *      15 - micro-flexo-printing bottom plane
 *      16 - micro-flexo-printing side plane
 *      17 - 2 radiused corners plus micro-flexo-printing motion
 *      18 - microprinting w/ force and torque
 *      19 - microprinting - axisymmetric
 *      20 - microprinting bottom roll surface
 *      21 - microprinting normal to bottom roll surface
 *      22 - microprinting - general roll, speed
 *      25 - fluid bearing die lip
 *      26 - syringe with curved lip
 *      27 - multiple plane, radiused corner combinations
 *      30 - falling circle
 *      31 - Radiused upstream bead
 */
dbl fnc( x1, x2,  x3, p,  time)
     const dbl x1;
     const dbl x2;
     const dbl x3;
     const dbl p[];
     const dbl time;
{
/* for PI use M_PIE Constant from std.h include file. */
	dbl f;
	int model_id;
	model_id = ((int)p[0]);


	/*   Gaussian bump on a roll surface
	 *      p[1] = rad = roll radius
	 *      p[2] = bump amplitude
	 *      p[3] = initial location of bump - radians
	 *      p[4] = width of bump
	 *      p[5] = web speed
  	 *	p[6] = xcenter
 	 *      p[7] = ycenter
	 */
	if (model_id == 1)
	{
	dbl rad,ampl,th0,wid,web_sp,arc_position,delta_rad;
        rad=p[1];
        ampl=p[2];
        web_sp=p[5];
        th0=p[3]+(web_sp/rad)*time;
        wid=p[4];
	arc_position = rad*(atan2(x2-p[7],x1-p[6])-th0)/wid;
	delta_rad = ampl*exp(-arc_position*arc_position);

        f=SQUARE(x1-p[6]) + SQUARE(x2-p[7]) - SQUARE(rad+delta_rad);

	}   /*  end of Gaussian bump if block */

	/*   cube corner and gravure geometries
	 *      - specific to flat geometries for right now
	 *      p[1] = cube_dep = depth of cube corner
	 *      p[2] = cube_lng = length across top of cube corner
	 *      p[3] = flat = length of flat between cubes
	 *      p[4] = flat_bot = length of flat at bottom of cubes
	 *      p[5] = radius = radius of corner fillet
	 *      p[6] = xveloc = web speed
	 *      p[7] = zangle = angle of groove in z-direction
	 *      p[8] = xoffset = starting offset in x-direction
	 *      p[9] = roll_rad
	 *	p[10] = xcenter = x-coord of roll center
	 *	p[11] = ycenter = y-coord of roll center
	 *
	 *      other variables:
	 *      xp = xprime between (0,cube_lng+flat+flat_bot)
	 *      slope = slope of cube walls
	 *	roll_rad = roll radius
	 *	arc_pos = arc_length position on roll from TDC
	 *	ysurf = deviation from smooth surface
	 */
	else 
	if (model_id == 2)
	{
	dbl xp, slope;
	dbl cube_dep, cube_lng, flat, radius, xveloc, zangle, xoffset;
	dbl xcenter, ycenter, roll_rad, arc_pos, ysurf=0, flat_bot, xref=0;

	double x_dist[4], theta;
        cube_dep = p[1];
        cube_lng = p[2];
        flat = p[3];
        flat_bot = p[4];
        radius = p[5];
        xveloc = p[6];
        zangle = p[7];
        xoffset = p[8];
        roll_rad = p[9];
	xcenter = p[10];
	ycenter = p[11];

	theta = atan2(cube_lng*0.5,cube_dep);
	x_dist[0] = radius/tan(0.5*theta+0.25*M_PIE);
	x_dist[1] = x_dist[0]*sin(theta);
	flat = MAX(flat,2*x_dist[0]);
	flat_bot = MAX(flat_bot,2*x_dist[0]);
        slope = 2.*cube_dep/cube_lng;
	arc_pos = roll_rad * (0.5*M_PIE - acos((x1-xcenter)
		/sqrt(SQUARE(x1-xcenter)+SQUARE(x2-ycenter))));

        xp = fmod(arc_pos-xveloc*time-zangle*x3-xoffset,cube_lng+flat+flat_bot);
        if(xp < 0.0){xp += cube_lng+flat+flat_bot;}

	if(xp >= (-0.5*flat+x_dist[0]) && xp <= (0.5*flat-x_dist[0]))
        	{ ysurf = 0.0; }

	else if(xp >= (0.5*flat-x_dist[0]) && xp <= (0.5*flat+x_dist[1]))
        	{ 
		xref = 0.5*flat-x_dist[0];
		ysurf = -radius + sqrt(SQUARE(radius)-SQUARE(xp-xref)); 
		}
        else if(xp >= (0.5*flat+x_dist[1]) && xp <= (0.5*flat+0.5*cube_lng-x_dist[1]))
        	{
		xref = 0.5*flat;
        	ysurf =  -slope*(xp - xref);
        	}

        else if(xp >= (0.5*flat+0.5*cube_lng-x_dist[1]) && 
		xp <= (0.5*flat+0.5*cube_lng+x_dist[0]))
        	{
		xref = 0.5*flat + 0.5*cube_lng + x_dist[0];
        	ysurf = -cube_dep + radius - sqrt(SQUARE(radius) - SQUARE(xp-xref));
        	}

        else if(xp >= (0.5*flat+0.5*cube_lng+x_dist[0]) && 
		xp <= (0.5*flat+0.5*cube_lng+flat_bot-x_dist[0]))
        	{
        	ysurf = -cube_dep;
        	}

        else if(xp >= (0.5*flat+0.5*cube_lng+flat_bot-x_dist[0]) && 
		xp <= (0.5*flat+0.5*cube_lng+flat_bot+x_dist[1]))
        	{
		xref = 0.5*flat + 0.5*cube_lng + flat_bot - x_dist[0];
        	ysurf = -cube_dep + radius - sqrt(SQUARE(radius) - SQUARE(xp - xref));
        	}

        else if(xp >= (0.5*flat+0.5*cube_lng+flat_bot+x_dist[1]) && 
		xp <= (0.5*flat+cube_lng+flat_bot-x_dist[1]))
        	{
		xref = 0.5*flat + 0.5*cube_lng + flat_bot;
        	ysurf = -cube_dep + slope*(xp-xref);
        	}
	else if(xp >= (0.5*flat+cube_lng+flat_bot-x_dist[1]) && 
		xp <= (0.5*flat+cube_lng+flat_bot+x_dist[0]))
        	{
		xref = 0.5*flat + cube_lng + flat_bot + x_dist[0];
		ysurf = -radius + sqrt(SQUARE(radius)-SQUARE(xp-xref));
        	}
	else if(xp >= (0.5*flat+cube_lng+flat_bot+x_dist[0]) && 
		xp <= (flat+cube_lng+flat_bot-x_dist[0]))
        	{
		ysurf = 0.0;
        	}
	else if(xp >= (flat+cube_lng+flat_bot-x_dist[0]) && 
		xp <= (flat+cube_lng+flat_bot+x_dist[1]))
        	{ 
		xref = flat+cube_lng+flat_bot-x_dist[0];
		ysurf = -radius + sqrt(SQUARE(radius)-SQUARE(xp-xref)); 
		}
        else
        	{
      		EH(-1,"invalid gravure user_bc position");
        	}

	f = SQUARE(x1-xcenter) + SQUARE(x2-ycenter) - SQUARE(roll_rad+ysurf);

	}	/**  end of gravure if block  */

	/*   radiused lip corner
	 *      p[1] = xpt - xcoord of die tip
	 *      p[2] = ypt - ycoord of die tip
	 *      p[3] = angle (radians) of side 1 (in CCW direction)
	 *      p[4] = angle (radians) of side 2 (in CCW direction)
	 *      p[5] = radius of corner
	 *	xcen = x-coord of radius center
	 *	ycen = y-coord of radius center
	 *	alpha = subtended half angle
	 */
	else 
	if (model_id == 3)
	{
	  dbl xpt, ypt, theta1, theta2, rad, xcen , ycen, alpha;
	  dbl theta;
	  xpt=p[1];
	  ypt=p[2];
	  theta1=p[3];
	  theta2=p[4];
	  rad=p[5];
	  
	  alpha = 0.5*(theta2-theta1);
	  xcen = xpt + (rad/sin(alpha))*cos(theta1+alpha);
	  ycen = ypt + (rad/sin(alpha))*sin(theta1+alpha);
	  
	  /**   compute angle of point on curve from arc center **/
	  
	  theta = atan2(x2-ycen,x1-xcen);
	  theta = theta > theta2-1.5*M_PIE ? theta : theta + 2*M_PIE;  
	  
	  /**  use different f depending on theta  **/
	  
	  if( (theta1-0.5*M_PIE) <= theta && theta <= (theta1+alpha))
	    {
	      f = (x2-ypt)*cos(theta1) - (x1-xpt)*sin(theta1);
	    }
	  else if ( (theta1+alpha) <= theta && (theta - 0.5*M_PIE) <= theta2)
	    {
	      f = (x2-ypt)*cos(theta2) - (x1-xpt)*sin(theta2);
	    }
	  else
	    {
	      f = SQUARE(x1-xcen)+SQUARE(x2-ycen)-SQUARE(rad);
	    }
      
	}   /*  end of radiused lip corner */
  
	/*
	   FITTED CURVE
	   SPECIFIED BY PADE COEFFICIENTS
	   
                                  2       3               N-1
                  p1 + p2 x + p3 x  + p4 x  + . . . + pN x
     f(x) = ------------------------------------------------------
                                  2         3                  N-2
            pN+1 + pN+2 x + pN+3 x  + pN+4 x  + . . . + p2N-1 x
     
	    HAVE TO PROVIDE COEFFICIENTS

	    p[0] = model_id
	    p[1] = N 
	    p[2] = 21: Y AS A FUNCTION OF X
                   31: Z AS A FUNCTION OF X
                   12: X AS A FUNCTION OF Y
                   32: Z AS A FUNCTION OF Y
                   13: X AS A FUNCTION OF Z
                   23: Y AS A FUNCTION OF Z
	    p[3] = coefficients

	    */
	else
	if (model_id == 4) {
	    
	  dbl w1=0, w2=0, dnum, dden; 
	  const dbl *pn, *pd;
	  int i, n, xyz;
	  
	  n      = abs(((int)p[1]));
	  xyz = ((int)p[2]);
	  
	  
	  pn = &p[3];
	  pd = &p[n+3];
	  
	switch (xyz)
		{
		case 21:
			w1=x1; w2=x2;break;
		case 31:
			w1=x1; w2=x3; break;
		case 12:
			w1=x2; w2=x1; break;
		case 32:
			w1=x2; w2=x3; break;
		case 13:
			w1=x3; w2=x1; break;
		case 23:
			w1=x3; w2=x2; break;
		default:
      EH(-1,"invalid xyz in Pade function user_bc");

		}

	  /*
	     
	  EVALUATE FUNCTION
	     
	  */
	  
	  dnum = pn[n-1];
	  dden = pd[n-2];
	  for (i=n-2;i>=0;i--)
	    { dnum = pn[i]+w1*dnum; }
	  for (i=n-3;i>=0;i--)
	    { dden = pd[i]+w1*dden; }
	  
	  f = dnum/dden-w2; } /* END OF FITTED CURVE */

	/*
	   CELL 
	   USES HYPERBOLIC TAN FUNCTIONS
	   
	   f(x) = y-(-d*0.5*(tanh(m*(x-s))+1.0)+d*0.5*tanh(m*(x-t))-1.0)+d)
	        = y-d*0.5*(-tanh(m*(x-s))+tanh(m*(x-t)))

	   HAVE TO PROVIDE COEFFICIENTS

	   p[0] = model_id
	   p[1] = d 
	   p[2] = m
	   p[3] = s
	   p[4] = width
	   p[5] = speed (transient analysis)
	   p[6] = coordinate direction (1 or 2)

	    */
	else
	if (model_id == 5) {
	    
	  dbl d, m, s, t, w1, w2;
          double coord, fcn;
          int dir;
	  
	  d = p[1];
	  m = p[2];
	  s = p[3]+p[5]*time;
	  t = s+p[4];
	  dir = ((int)p[6]);

	  if(dir == 1)	
		{coord = x1;fcn = x2;}
	  else
		{coord = x2; fcn = x1;}

	  w1 = tanh(m*(coord-s));
	  w2 = tanh(m*(coord-t));

	  f = fcn-d*0.5*(-w1+w2); 

	   } /* END OF CELL */

	    /*
		tapered surface for 3D manifolds
		driven by a Pade curve

	    p[0] = model_id
	    p[1] = theta, angle in xy-plane (degrees) 
	    p[2] = N, pade order
	    p[3 ... 2N+1] = pade coefficients

	    */
	else
	if (model_id == 6) {
	    
	  dbl x, dnum, dden, pade, th; 
	  const dbl *pn, *pd;
	  int i, n;
	  
	  th     = M_PIE * p[1]/180.;
	  n      = abs(((int)p[2]));
	  pn = &p[3];
	  pd = &p[n+3];
	  x = x3;
	  
	  /*
	  EVALUATE FUNCTION
	  */
	  
	  dnum = pn[n-1];
	  dden = pd[n-2];
	  for (i=n-2;i>=0;i--)
	    { dnum = pn[i]+x*dnum; }
	  for (i=n-3;i>=0;i--)
	    { dden = pd[i]+x*dden; }
	  
	  pade = dnum/dden;
	  
	  f = x2*cos(th) +sin(th)*(pade-x1);

	  } /* END OF tapered surface */

	/*   3D radiused lip corner
	 *      p[1] = angle (radians) of side 1 (in CCW direction)
	 *      p[2] = angle (radians) of side 2 (in CCW direction)
	 *      p[3] = radius of corner
	 *      p[4] = Nx, pade order for x0(z)
	 *      p[5 ... 2N+3] = pade coefficients
	 *      p[2Nx+4] = Ny, pade order for y0(z)
	 *      p[2Nx+5 ... 2Nx+2Ny+3] = pade coefficients
	 *	xcen = x-coord of radius center
	 *	ycen = y-coord of radius center
	 *	alpha = subtended half angle
	 */
	else
	if (model_id == 7)	{

	  dbl theta1, theta2, rad, xcen , ycen, alpha;
	  dbl theta;
	  dbl x, dnum, dden, padex, padey;
	  const dbl *pnx, *pdx, *pny, *pdy;
	  int i,nx, ny;
	  theta1=p[1];
	  theta2=p[2];
	  rad=p[3];
	  nx      = abs(((int)p[4]));
	  pnx = &p[5];
	  pdx = &p[nx+5];
	  ny      = abs(((int)p[nx+6]));
	  pny = &p[nx+7];
	  pdy = &p[nx+7+ny];
	  x = x3;
	  
	  /*
	  EVALUATE FUNCTION
	  */
	  
	  dnum = pnx[nx-1];
	  dden = pdx[nx-2];
	  for (i=nx-2;i>=0;i--)
	    { dnum = pnx[i]+x*dnum; }
	  for (i=nx-3;i>=0;i--)
	    { dden = pdx[i]+x*dden; }
	  
	  padex = dnum/dden;

	  dnum = pny[ny-1];
	  dden = pdy[ny-2];
	  for (i=ny-2;i>=0;i--)
	    { dnum = pny[i]+x*dnum; }
	  for (i=ny-3;i>=0;i--)
	    { dden = pdy[i]+x*dden; }
	  
	  padey = dnum/dden;

	  alpha = 0.5*(theta2-theta1);
	  xcen = padex + (rad/sin(alpha))*cos(theta1+alpha);
	  ycen = padey + (rad/sin(alpha))*sin(theta1+alpha);

	  /**   compute angle of point on curve from arc center **/
	  
	  theta = atan2(x2-ycen,x1-xcen);
	  theta = theta > theta2-1.5*M_PIE ? theta : theta + 2*M_PIE;  
	  
	  /**  use different f depending on theta  **/
	  
	  if( (theta1-0.5*M_PIE) <= theta && theta <= (theta1+alpha))
	    {
	      f = (x2-padey)*cos(theta1) - (x1-padex)*sin(theta1);
	    }
	  else if ( (theta1+alpha) <= theta && (theta - 0.5*M_PIE) <= theta2)
	    {
	      f = (x2-padey)*cos(theta2) - (x1-padex)*sin(theta2);
	    }
	  else
	    {
	      f = SQUARE(x1-xcen)+SQUARE(x2-ycen)-SQUARE(rad);
	    }
      
	  }

	    /*
		SINGLE V CELL

	    p[0] = model_id = 8
	    p[1] = x01
	    p[2] = y01 
	    p[3] = H
	    p[4] = L
	    p[5] = A
	    p[6] = R
	    p[7] = G

	    */
	else
	if (model_id == 8)	{

	  dbl 
	    x01, x02, x03, x04, x05, x06, x07, y01, y02, y03, y04, y05,
	    RcA, RsA, tA, xR, xC, xL, yR, yC, yL, H, A, R;

	x01 = p[1];
	y01 = p[2];
	H   = p[3];
	A   = p[5];
	R   = p[6];

	RcA = R*cos(A);
	RsA = R*sin(A);
	tA  = tan(A);

	y02 = y01+R-RsA;
	x02 = x01-RcA;

	y03 = y02;
	x03 = x01+RcA;
	y04 = y01+H-R+RsA;
	x04 = x02-(y04-y02)*tA;
	y05 = y04;
	x05 = x03+(y05-y03)*tA;
	x06 = x04-RcA;
	x07 = x05+RcA;

	xC = x01;
	yC = y01+R;

	xR = x05+RcA;
	yR = y05-RsA;

	xL = x04-RcA;
	yL = y04-RsA;

	  /*
	  EVALUATE FUNCTION
	  */
	  
	if (x2 <= x06)
	  {
	    f = x1-H-y01;
	  }
	else
	if ((x2 <= x04) && (x2 > x06))
	  {
	    f = (x2-xL)*(x2-xL)+(x1-yL)*(x1-yL)-R*R;
	  }
	else
	if ((x2 <= x02) && (x2 > x04))
	  {
	    f = (x04-x02)*(x1-y02)-(y04-y02)*(x2-x02);
	  }
	else
	if ((x2 <= x03) && (x2 > x02))
	  {
	    f = (x2-xC)*(x2-xC)+(x1-yC)*(x1-yC)-R*R;
	  }
	else
	if ((x2 <= x05) && (x2 > x03))
	  {
	    f = (x05-x03)*(x1-y03)-(y05-y03)*(x2-x03);
	  }
	else
	if ((x2 <= x07) && (x2 > x05))
	  {
	    f = (x2-xR)*(x2-xR)+(x1-yR)*(x1-yR)-R*R;
	  }
	else
	if (x2 > x07)
	  {
	    f = x1-H-y01;
	  } }
        else
        if (model_id == 9)
        {
        dbl x0,y0,ampl,speed,width,xpt;
          x0=p[1];
          y0=p[2];
          ampl=p[3];
          speed=p[4];
          width=p[5];
          xpt = x0+speed*time;
  
         f=x2 - (y0+ampl*exp(-SQUARE((x1-xpt)/width)));
  
        }   /*  end of Doug's Gaussian bump if block */
	/*   normal line to radiused lip corner
	 *      p[1] = xpt - xcoord of die tip
	 *      p[2] = ypt - ycoord of die tip
	 *      p[3] = angle (radians) of side 1 (in CCW direction)
	 *      p[4] = angle (radians) of side 2 (in CCW direction)
	 *      p[4] = radius of corner
	 *	xcen = x-coord of radius center
	 *	ycen = y-coord of radius center
	 *	alpha = subtended half angle
	 */
	else 
	if (model_id == 10)
	{
	  dbl xpt, ypt, theta1, theta2, rad, xcen , ycen, alpha;
	  dbl theta;
	  xpt=p[1];
	  ypt=p[2];
	  theta1=p[3];
	  theta2=p[4];
	  rad=p[5];
	  
	  alpha = 0.5*(theta2-theta1);
	  xcen = xpt + (rad/sin(alpha))*cos(theta1+alpha);
	  ycen = ypt + (rad/sin(alpha))*sin(theta1+alpha);
	  
	  /**   compute angle of point on curve from arc center **/
	  
	  theta = atan2(x2-ycen,x1-xcen);
	  theta = theta > theta2-1.5*M_PIE ? theta : theta + 2*M_PIE;  
	  
	  /**  use different f depending on theta  **/
	  
	  if( (theta1-0.5*M_PIE) <= theta && theta <= (theta1+alpha))
	    {
	      f = (x2-ypt)*sin(theta1) + (x1-xpt)*cos(theta1);
	    }
	  else if ( (theta1+alpha) <= theta && (theta - 0.5*M_PIE) <= theta2)
	    {
	      f = (x2-ypt)*sin(theta2) + (x1-xpt)*cos(theta2);
	    }
	  else
	    {
	      f = (x2-ypt)*cos(theta) - (x1-xpt)*sin(theta);
	    }
      
	}   /*  end of normal line to radiused lip corner */
  

	/*   two radiused lip corners
	 *      p[1] = xpt - xcoord of die tip 1
	 *      p[2] = ypt - ycoord of die tip 1
	 *      p[3] = angle (radians) of side 1 (in CCW direction)
	 *      p[4] = radius of corner 1
	 *      p[4] = xpt - xcoord of die tip 2
	 *      p[5] = ypt - ycoord of die tip 2
	 *      p[6] = angle (radians) of side 2 (in CCW direction)
	 *      p[7] = radius of corner 2
	 *	xcen = x-coord of radius center
	 *	ycen = y-coord of radius center
	 *	alpha = subtended half angle
	 */
	else 
	if (model_id == 11)
	{
	  dbl xpt1, ypt1, theta1, rad1, xcen1 , ycen1, alpha1;
	  dbl xpt2, ypt2, theta2, rad2, xcen2 , ycen2, alpha2;
	  dbl theta1m, theta2m, th1, th2, th2t;
	  xpt1=p[1];
	  ypt1=p[2];
	  theta1=p[3];
	  rad1=p[4];
	  xpt2=p[5];
	  ypt2=p[6];
	  theta2=p[7];
	  rad2=p[8];
	  
	/*  slope of middle line		*/

	  theta1m = atan2(ypt2-ypt1,xpt2-xpt1);
	  theta1m = theta1m > theta1 ? theta1m : theta1m + 2*M_PIE;  
	  theta2m = atan2(ypt1-ypt2,xpt1-xpt2);
	  alpha1 = 0.5*(theta1m-theta1);
	  alpha2 = 0.5*(theta2-theta2m);

	  xcen1 = xpt1 + (rad1/sin(alpha1))*cos(theta1+alpha1);
	  ycen1 = ypt1 + (rad1/sin(alpha1))*sin(theta1+alpha1);
	  xcen2 = xpt2 + (rad2/sin(alpha2))*cos(theta2m+alpha2);
	  ycen2 = ypt2 + (rad2/sin(alpha2))*sin(theta2m+alpha2);
	  
	  /**   compute angle of point on curve from arc center **/
	  
	  th1 = atan2(x2-ycen1,x1-xcen1);
	  th2 = atan2(x2-ycen2,x1-xcen2);
	  th2t = th2 > 0.0 ? th2 : th2 + 2*M_PIE;  

	  /**  use different f depending on theta  **/
	  
	  if( (theta1-0.5*M_PIE) <= th1 && th1 <= theta1)
	    {
	      f = (x2-ypt1)*cos(theta1) - (x1-xpt1)*sin(theta1);
	    }
	  else if ( theta2 <= th2t && (th2t - 0.5*M_PIE) <= theta2)
	    {
	      f = (x2-ypt2)*cos(theta2) - (x1-xpt2)*sin(theta2);
	    }
	  else if ( theta2m <= (th1+0.5*M_PIE) && th1 <= (theta1-0.5*M_PIE))
	    {
	      f = SQUARE(x1-xcen1)+SQUARE(x2-ycen1)-SQUARE(rad1);
	    }
	  else if ( (theta2m-0.5*M_PIE) <= th2 && (th1+0.5*M_PIE) <= theta2m)
	    {
	      f = (x2-ypt1)*cos(theta1m) - (x1-xpt1)*sin(theta1m);
	    }
	  else
	    {
	      f = SQUARE(x1-xcen2)+SQUARE(x2-ycen2)-SQUARE(rad2);
	    }
      
	}   /*  end of 2 radiused lip corners */
	else 
	if (model_id == 12)
	{
	  dbl xcen1 , ycen1, rad_init, alpha, radius;
	  xcen1=p[1];
	  ycen1=p[2];
	  rad_init=p[3];
	  alpha=p[4];
	  radius = rad_init + alpha*time;
	      f = SQUARE(x1-xcen1)+SQUARE(x2-ycen1)-SQUARE(radius);
	}
	else 
	if (model_id == 13)
	{
	  dbl xpt, ypt, theta1, theta2, rad, xcen , ycen, alpha;
	  dbl theta, xcirc, ycirc, dierad;
	  int iside;
	  xpt=p[1];
	  ypt=p[2];
	  theta1=p[3];
	  theta2=p[4];
	  rad=p[5];
	  dierad=p[6];
	  iside = ((int)p[7]);
	  
	  /**  find center of die face  **/

	if( iside == 1)
	  {
	  xcirc = xpt + dierad*sin(theta1);
	  ycirc = ypt - dierad*cos(theta1);
	  }	else	{
	  xcirc = xpt - dierad*sin(theta2);
	  ycirc = ypt + dierad*cos(theta2);
	  }

	  /**  find center of corner radius  **/

	if( iside == 1)
	  	{
		alpha = theta2 - asin((rad-dierad*cos(theta2-theta1))/(rad+dierad));
	  	}	else	{
		alpha = theta1 + asin((rad-dierad*cos(theta2-theta1))/(rad+dierad));
	  	}
		xcen = xcirc + (rad+dierad)*cos(alpha);
		ycen = ycirc + (rad+dierad)*sin(alpha);
	  
	  /**   compute angle of point on curve from arc center **/
	  
	  theta = atan2(x2-ycen,x1-xcen);
	  theta = theta > alpha-M_PIE ? theta : theta + 2*M_PIE;  
	  
	  /**  use different f depending on theta  **/
	  
	if( iside == 1)
	  {
	  if( (alpha-M_PIE) <= theta && theta <= theta1)
	    {
	      f = SQUARE(x1-xcirc)+SQUARE(x2-ycirc)-SQUARE(dierad);
	    }
	  else if ( theta2 <= theta && (theta - 0.5*M_PIE) <= theta2)
	    {
	      f = (x2-ypt)*cos(theta2) - (x1-xpt)*sin(theta2);
	    }
	  else
	    {
	      f = SQUARE(x1-xcen)+SQUARE(x2-ycen)-SQUARE(rad);
	    }
	  }	else	
	  {
	  if( theta2 <= theta && theta <= (alpha+M_PIE))
	    {
	      f = SQUARE(x1-xcirc)+SQUARE(x2-ycirc)-SQUARE(dierad);
	    }
	  else if ( (theta1-0.5*M_PIE) <= theta && theta <= theta1)
	    {
	      f = (x2-ypt)*cos(theta1) - (x1-xpt)*sin(theta1);
	    }
	  else
	    {
	      f = SQUARE(x1-xcen)+SQUARE(x2-ycen)-SQUARE(rad);
	    }
	  }
      
	}   /*  end of radiused lip corner */
	/*
	   Needle collar based on a tanh modified plane equation:

		Ax + By +Cz + D + E(tanh(beta*(x-x0))+1)/2 = 0
	   USES HYPERBOLIC TAN FUNCTIONS
	    */
	else
	if (model_id == 14) {
	    
	  double A, B, C, D, E, beta, x0;
	  
	  A = p[1];
	  B = p[2];
	  C = p[3];
	  D = p[4];
	  E = p[5];
	  beta = p[6];
	  x0 = p[7];

	  f = A*x1 + B*x2 + C*x3 + D + E*(tanh(beta*(x1-x0))+1)/2.;
	  }
	/* micro-flexo-printing bottom surface - omega, roll_rad, gap */
	else
	if (model_id == 15) {
	    
	  double omega, roll_rad, gap, t_offset, time1;
	  omega = p[1];
	  roll_rad = p[2];
	  gap = p[3];
	  t_offset = p[4];

	  time1 = time + t_offset;
	  f = sin(2.*omega*time1)*x1 - cos(2.*omega*time1)*
			(x2-gap-2.*roll_rad*(1.-sin(omega*time1)));
	  }
	/* micro-flexo-printing vertical surface - omega, roll_rad, gap, width */
	else
	if (model_id == 16) {
	    
	  double omega, roll_rad, gap, width, t_offset, time1;
	  omega = p[1];
	  roll_rad = p[2];
	  gap = p[3];
	  width = p[4];
	  t_offset = p[5];

	  time1 = time + t_offset;

	  f = - cos(2.*omega*time1)*x1 - sin(2.*omega*time1)*
			(x2-gap-2.*roll_rad*(1.-sin(omega*time1))) - 0.5*width;
	  }
	/*   two radiused lip corners plus micro-flexo-printing
	 *      p[1] = omega - roll angular velocity
	 *      p[2] = roll radius
	 *      p[3] = roll gap
	 *      p[4] = width of feature
	 *      p[5] = time offset (PI/2*omega)
	 *      p[6] = radius of corner 1
	 *      p[7] = radius of corner 2
	 *	xcen = x-coord of radius center
	 *	ycen = y-coord of radius center
	 *	alpha = subtended half angle
	 */
	else 
	if (model_id == 17)
	{
	  dbl xpt1, ypt1, theta1, rad1, xcen1 , ycen1, alpha1;
	  dbl xpt2, ypt2, theta2, rad2, xcen2 , ycen2, alpha2;
	  dbl theta1m, theta2m, th1, th2, th2t;
	  double omega, roll_rad, gap, width, t_offset, time1, yc;
	  omega = p[1];
	  roll_rad = p[2];
	  gap = p[3];
	  width = p[4];
	  t_offset = p[5];

	  time1 = time + t_offset;
	  yc = gap + 2.*roll_rad*(1.-sin(omega*time1));
	  xpt1= -0.5*width*cos(2.*omega*time1);
	  ypt1= yc - 0.5*width*sin(2.*omega*time1);
	  theta1= omega*time1;
	  rad1=p[6];
	  xpt2= 0.5*width*cos(2.*omega*time1);
	  ypt2= yc + 0.5*width*sin(2.*omega*time1);
	  theta2= omega*time1;
	  rad2=p[7];
/*printf("position %g %g %g %g %g %g\n",time,time1,xpt1,ypt1,xpt2,ypt2);*/
	  
	/*  slope of middle line		*/

	  theta1m = atan2(ypt2-ypt1,xpt2-xpt1);
	  theta1m = theta1m > theta1 ? theta1m : theta1m + 2*M_PIE;  
	  theta2m = atan2(ypt1-ypt2,xpt1-xpt2);
	  alpha1 = 0.5*(theta1m-theta1);
	  alpha2 = 0.5*(theta2-theta2m);

	  xcen1 = xpt1 + (rad1/sin(alpha1))*cos(theta1+alpha1);
	  ycen1 = ypt1 + (rad1/sin(alpha1))*sin(theta1+alpha1);
	  xcen2 = xpt2 + (rad2/sin(alpha2))*cos(theta2m+alpha2);
	  ycen2 = ypt2 + (rad2/sin(alpha2))*sin(theta2m+alpha2);
	  
	  /**   compute angle of point on curve from arc center **/
	  
	  th1 = atan2(x2-ycen1,x1-xcen1);
	  th2 = atan2(x2-ycen2,x1-xcen2);
	  th2t = th2 > 0.0 ? th2 : th2 + 2*M_PIE;  

	  /**  use different f depending on theta  **/
	  
	  if( (theta1-0.5*M_PIE) <= th1 && th1 <= theta1)
	    {
	      f = (x2-ypt1)*cos(theta1) - (x1-xpt1)*sin(theta1);
	    }
	  else if ( theta2 <= th2t && (th2t - 0.5*M_PIE) <= theta2)
	    {
	      f = (x2-ypt2)*cos(theta2) - (x1-xpt2)*sin(theta2);
	    }
	  else if ( theta2m <= (th1+0.5*M_PIE) && th1 <= (theta1-0.5*M_PIE))
	    {
	      f = SQUARE(x1-xcen1)+SQUARE(x2-ycen1)-SQUARE(rad1);
	    }
	  else if ( (theta2m-0.5*M_PIE) <= th2 && (th1+0.5*M_PIE) <= theta2m)
	    {
	      f = (x2-ypt1)*cos(theta1m) - (x1-xpt1)*sin(theta1m);
	    }
	  else
	    {
	      f = SQUARE(x1-xcen2)+SQUARE(x2-ycen2)-SQUARE(rad2);
	    }
      
	}   /*  end of micro-printing (original) */
	/*   two radiused lip corners plus micro-flexo-printing
	 *     with accounting for displacement and twist
	 *      p[1] = omega - roll angular velocity
	 *      p[2] = roll radius
	 *      p[3] = roll gap
	 *      p[4] = width of feature
	 *      p[5] = time offset (PI/2*omega)
	 *      p[6] = radius of corner 1
	 *      p[7] = radius of corner 2
	 *      p[8] = vertical displacement
	 *      p[9] = added twist
	 *	xcen = x-coord of radius center
	 *	ycen = y-coord of radius center
	 *	alpha = subtended half angle
	 */
	else 
	if (model_id == 18)
	{
	  dbl xpt1, ypt1, theta1, rad1, xcen1 , ycen1, alpha1;
	  dbl xpt2, ypt2, theta2, rad2, xcen2 , ycen2, alpha2;
	  dbl theta1m, theta2m, th1, th2, th2t;
	  double omega, roll_rad, gap, width, t_offset, time1, yc, xc;
	  double y_displ, twist, angle, l_arm;
	  omega = p[1];
	  roll_rad = p[2];
	  gap = p[3];
	  width = p[4];
	  t_offset = p[5];
	  rad1=p[6];
	  rad2=p[7];
	  y_displ=p[8];
	  twist=p[9];
	  l_arm=p[10];

	  time1 = time + t_offset;
	  angle = omega*time1;
	  yc = gap + 2.*roll_rad*(1.-sin(angle));
/*  account for displacement	*/
	  yc += y_displ*sin(angle);
	  xc = y_displ*cos(angle);
/*  account for twist	*/
	  yc += l_arm*(sin(angle)-sin(angle+twist));
	  xc += l_arm*(cos(angle)-cos(angle+twist));
	  angle += twist;
	  xpt1= xc + 0.5*width*sin(angle);
	  ypt1= yc - 0.5*width*cos(angle);
	  theta1= angle;
	  xpt2= xc - 0.5*width*sin(angle);
	  ypt2= yc + 0.5*width*cos(angle);
	  theta2= angle;
	  
	/*  slope of middle line		*/

	  theta1m = atan2(ypt2-ypt1,xpt2-xpt1);
	  theta1m = theta1m > theta1 ? theta1m : theta1m + 2*M_PIE;  
	  theta2m = atan2(ypt1-ypt2,xpt1-xpt2);
	  alpha1 = 0.5*(theta1m-theta1);
	  alpha2 = 0.5*(theta2-theta2m);

	  xcen1 = xpt1 + (rad1/sin(alpha1))*cos(theta1+alpha1);
	  ycen1 = ypt1 + (rad1/sin(alpha1))*sin(theta1+alpha1);
	  xcen2 = xpt2 + (rad2/sin(alpha2))*cos(theta2m+alpha2);
	  ycen2 = ypt2 + (rad2/sin(alpha2))*sin(theta2m+alpha2);
	  
	  /**   compute angle of point on curve from arc center **/
	  
	  th1 = atan2(x2-ycen1,x1-xcen1);
	  th2 = atan2(x2-ycen2,x1-xcen2);
	  th2t = th2 > 0.0 ? th2 : th2 + 2*M_PIE;  

	  /**  use different f depending on theta  **/
	  
	  if( (theta1-0.5*M_PIE) <= th1 && th1 <= theta1)
	    {
	      f = (x2-ypt1)*cos(theta1) - (x1-xpt1)*sin(theta1);
	    }
	  else if ( theta2 <= th2t && (th2t - 0.5*M_PIE) <= theta2)
	    {
	      f = (x2-ypt2)*cos(theta2) - (x1-xpt2)*sin(theta2);
	    }
	  else if ( theta2m <= (th1+0.5*M_PIE) && th1 <= (theta1-0.5*M_PIE))
	    {
	      f = SQUARE(x1-xcen1)+SQUARE(x2-ycen1)-SQUARE(rad1);
	    }
	  else if ( (theta2m-0.5*M_PIE) <= th2 && (th1+0.5*M_PIE) <= theta2m)
	    {
	      f = (x2-ypt1)*cos(theta1m) - (x1-xpt1)*sin(theta1m);
	    }
	  else
	    {
	      f = SQUARE(x1-xcen2)+SQUARE(x2-ycen2)-SQUARE(rad2);
	    }
      
	}   /*  end of micro-printing (twist) */
	/*   micro-flexo-printing - axi-symmetric - one radius
	 *      p[1] = omega - roll angular velocity
	 *      p[2] = roll radius
	 *      p[3] = roll gap
	 *      p[4] = width of feature
	 *      p[5] = time offset (PI/2*omega)
	 *      p[6] = radius of corner 1
	 *      p[7] = radius of corner 2
	 *	xcen = x-coord of radius center
	 *	ycen = y-coord of radius center
	 *	alpha = subtended half angle
	 */
	else 
	if (model_id == 19)
	{
	  dbl xpt, ypt, theta1, theta2, rad, xcen , ycen, alpha;
	  dbl theta, tdc_delay;
	  double omega, roll_rad, gap, width, t_offset, time1, yc;
	  omega = p[1];
	  roll_rad = p[2];
	  gap = p[3];
	  width = p[4];
	  t_offset = p[5];
	  rad=p[6];
	  theta2=p[7];
	  tdc_delay=p[8];

	  time1 = time + t_offset;
	  if (time <= 0.0)	
		{ yc = gap + 2.*roll_rad*(1.-sin(omega*time1));}
	  else if(time <= tdc_delay)
		{ yc = gap ;}
	  else 
		{ time1 -= tdc_delay;
		 yc = gap + 2.*roll_rad*(1.-sin(omega*time1));
		}
	  xpt=yc;
	  ypt=0.5*width;
	  theta1= -0.5*M_PIE;
	  
	  alpha = 0.5*(theta2-theta1);
	  xcen = xpt + (rad/sin(alpha))*cos(theta1+alpha);
	  ycen = ypt + (rad/sin(alpha))*sin(theta1+alpha);
	  
	  /**   compute angle of point on curve from arc center **/
	  
	  theta = atan2(x2-ycen,x1-xcen);
	  theta = theta > theta2-1.5*M_PIE ? theta : theta + 2*M_PIE;  
	  
	  /**  use different f depending on theta  **/
	  
	  if( (theta1-0.5*M_PIE) <= theta && theta <= (theta1+alpha))
	    {
	      f = (x2-ypt)*cos(theta1) - (x1-xpt)*sin(theta1);
	    }
	  else if ( (theta1+alpha) <= theta && (theta - 0.5*M_PIE) <= theta2)
	    {
	      f = (x2-ypt)*cos(theta2) - (x1-xpt)*sin(theta2);
	    }
	  else
	    {
	      f = SQUARE(x1-xcen)+SQUARE(x2-ycen)-SQUARE(rad);
	    }
      
	}   /*  end of micro-printing (axi-symmetric) */
	/*   micro-flexo-printing bottom roll surface
	 *      p[1] = omega - roll angular velocity
	 *      p[2] = roll radius
	 *      p[3] = time offset (PI/2*omega)
	 */
	else 
	if (model_id == 20)
	{
	  double omega, roll_rad, t_offset, coat_th, time1;
	  omega = p[1];
	  roll_rad = p[2];
	  t_offset = p[3];
	  coat_th = p[4];

	  time1 = time + t_offset;
	  f = x1*(x1-2.*roll_rad*cos(omega*time1)) +
		x2*(x2+2.*roll_rad*sin(omega*time1))
		- coat_th*(coat_th+2.*roll_rad);
	}
	/*   micro-flexo-printing  normal to bottom roll surface
	 *      p[1] = omega - roll angular velocity
	 *      p[2] = roll radius
	 *      p[3] = time offset (PI/2*omega)
	 *      p[4] = x- offset 
	 */
	else 
	if (model_id == 21)
	{
	  double omega, roll_rad, t_offset, time1, x_offset;
	  omega = p[1];
	  roll_rad = p[2];
	  t_offset = p[3];
	  x_offset = p[4];

	  time1 = time + t_offset;
	  f = (x1-roll_rad*cos(omega*time1))
			*sin(omega*time1+x_offset/roll_rad)
		+ (x2+roll_rad*sin(omega*time1))
			*cos(omega*time1+x_offset/roll_rad);
	}
	/*   two radiused lip corners plus micro-flexo-printing
	 *     with accounting for displacement and twist
	 *	general wrt rolls, speeds, angles
	 *      p[1] = omega_lower - roll angular velocity
	 *      p[2] = omega_upper - roll angular velocity
	 *      p[3] = lower roll radius
	 *      p[4] = upper roll radius
	 *      p[5] = roll gap
	 *      p[6] = width of feature
	 *      p[7] = time offset (PI/2*omega)
	 *      p[8] = radius of corner 1
	 *      p[9] = radius of corner 2
	 *      p[10] = angle of sides
	 *      p[11] = vertical displacement
	 *      p[12] = added twist
	 *      p[13] = length of lever arm
	 *	xcen = x-coord of radius center
	 *	ycen = y-coord of radius center
	 *	alpha = subtended half angle
	 */
	else 
	if (model_id == 22)
	{
	  dbl xpt1, ypt1, theta1, rad1, xcen1 , ycen1, alpha1;
	  dbl xpt2, ypt2, theta2, rad2, xcen2 , ycen2, alpha2;
	  dbl theta1m, theta2m, th1, th2, th2t;
	  double omega_l, omega_u, roll_rad_l, roll_rad_u, gap, width;
	  double t_offset, time1, yc, xc, alpha;
	  double y_displ, twist, angle_l, angle_u, l_arm;
	  omega_l = p[1];
	  omega_u = p[2];
	  roll_rad_l = p[3];
	  roll_rad_u = p[4];
	  gap = p[5];
	  width = p[6];
	  t_offset = p[7];
	  rad1=p[8];
	  rad2=p[9];
	  alpha=p[10];
	  y_displ=p[11];
	  twist=p[12];
	  l_arm=p[13];

	  time1 = time + t_offset;
	  angle_l = omega_l*time1;
	  angle_u = omega_u*time1 + 0.5*M_PIE*(1.-omega_u/omega_l);
	  yc = gap + roll_rad_l*(1.-sin(angle_l)) + roll_rad_u*(1.-sin(angle_u));
	  xc = -roll_rad_u*cos(angle_u) + roll_rad_l*cos(angle_l);
/*  account for displacement	*/
	  yc += y_displ*sin(angle_u);
	  xc += y_displ*cos(angle_u);
/*  account for twist	*/
	  yc += l_arm*(sin(angle_u)-sin(angle_u+twist));
	  xc += l_arm*(cos(angle_u)-cos(angle_u+twist));
	  angle_u += twist;
	  xpt1= xc + 0.5*width*sin(angle_u);
	  ypt1= yc - 0.5*width*cos(angle_u);
	  theta1= angle_u - alpha;
	  xpt2= xc - 0.5*width*sin(angle_u);
	  ypt2= yc + 0.5*width*cos(angle_u);
	  theta2= angle_u + alpha;
	  
	/*  slope of middle line		*/

	  theta1m = atan2(ypt2-ypt1,xpt2-xpt1);
	  theta1m = theta1m > theta1 ? theta1m : theta1m + 2*M_PIE;  
	  theta2m = atan2(ypt1-ypt2,xpt1-xpt2);
	  alpha1 = 0.5*(theta1m-theta1);
	  alpha2 = 0.5*(theta2-theta2m);

	  xcen1 = xpt1 + (rad1/sin(alpha1))*cos(theta1+alpha1);
	  ycen1 = ypt1 + (rad1/sin(alpha1))*sin(theta1+alpha1);
	  xcen2 = xpt2 + (rad2/sin(alpha2))*cos(theta2m+alpha2);
	  ycen2 = ypt2 + (rad2/sin(alpha2))*sin(theta2m+alpha2);
	  
	  /**   compute angle of point on curve from arc center **/
	  
	  th1 = atan2(x2-ycen1,x1-xcen1);
	  th2 = atan2(x2-ycen2,x1-xcen2);
	  th2t = th2 > 0.0 ? th2 : th2 + 2*M_PIE;  

	  /**  use different f depending on theta  **/
	  
	  if( (theta1-0.5*M_PIE) <= th1 && th1 <= theta1)
	    {
	      f = (x2-ypt1)*cos(theta1) - (x1-xpt1)*sin(theta1);
	    }
	  else if ( theta2 <= th2t && (th2t - 0.5*M_PIE) <= theta2)
	    {
	      f = (x2-ypt2)*cos(theta2) - (x1-xpt2)*sin(theta2);
	    }
	  else if ( theta2m <= (th1+0.5*M_PIE) && th1 <= (theta1-0.5*M_PIE))
	    {
	      f = SQUARE(x1-xcen1)+SQUARE(x2-ycen1)-SQUARE(rad1);
	    }
	  else if ( (theta2m-0.5*M_PIE) <= th2 && (th1+0.5*M_PIE) <= theta2m)
	    {
	      f = (x2-ypt1)*cos(theta1m) - (x1-xpt1)*sin(theta1m);
	    }
	  else
	    {
	      f = SQUARE(x1-xcen2)+SQUARE(x2-ycen2)-SQUARE(rad2);
	    }
	}   /*  end of micro-printing (general) */
	/*   fluid bearing die lip by Slah Jendoubi
         *      BC = GEOM SS 111 14 {x6} {y6} {x5} {y5} {x29} {y29} {lip_rad_u} 
         *      BC = GEOM SS 111 14 {xx1} {yy1} {xx2} {yy2} {xcen} {ycen} {lip_rad_u} 
	 *      p[1] = x6 - xcoord of (the point on the flat surface)
	 *      p[2] = y6 - ycoord of (the point on the flat surface)
	 *      p[3] = x5 - xcoord of (intersection of flat and radiused surfaces)
	 *      p[4] = y5 - ycoord of (intersection of flat and radiused surfaces)
	 *	p[5] = xcen = x-coord of radius center
	 *	p[6] = ycen = y-coord of radius center
	 *      p[7] = rad = radius of corner
	 *	theta*= angle
	 *	alpha  = angle
         *
         *
         *                                           * 1
         *                                           |
         *                                           |
         *                                           |
         *                                           |
         *                                           |straight line
         *                                           |
         *              *(xc,yc)                     |
         *               \                           |
         *                 \       *(xcnew,ycnew)    |
         *                   \                       |
         *                     \                     |
         *                       \rad                |
         *                         \                 |
         *                           \               * 2
         *                             \           /
         *                               \       /
         *                                 \  /
         *                                 /circular Arc
         *                              /
         *                           /
         *                       /
         *                   / 
         *                 * 3
         *
         *
         *
	 */
	else
	if (model_id == 25)
	{
	  dbl xx1, yy1, xx2, yy2, xx3, yy3, xcen , ycen, rad, alpha, theta, theta1;
	  dbl xcen_new,ycen_new;
	  xx1=p[1];
	  yy1=p[2];
	  xx2=p[3];
	  yy2=p[4];
	  xx3=p[5];
	  yy3=p[6];
 	  xcen=p[7];
 	  ycen=p[8];
 	  rad=p[9];
	  
	  /**   compute angle of point on curve from (xcen_new,ycen_new) **/
	  
	  xcen_new = (xx1+xx3)/2;
	  ycen_new = (yy1+yy3)/2.;
	  theta  = atan2(x2-ycen_new,x1-xcen_new);
	  theta1 = atan2(yy2-ycen_new,xx2-xcen_new);
	  alpha  = atan2(yy1-yy2,xx1-xx2);

          /**  use different f depending on theta  **/

          if( theta <= theta1)
            {
              f = SQUARE(x1-xcen)+SQUARE(x2-ycen)-SQUARE(rad);
            }
          else
            {
              f = (x2-yy2)*cos(alpha ) - (x1-xx2)*sin(alpha );
            }
              /*printf("theta %lf theta1 %lf DT %lf x %lf y %lf f %lf \n\n", theta,theta1,theta-theta1,x1,x2,f);*/

        }   /*  end of fluid bearing die */
	/*    syringe with curved lip  by Slah Jendoubi
	 *      p[1] = x1 - xcoord of point 1
	 *      p[2] = y1 - ycoord of point 1
	 *      p[3] = x2 - xcoord of point 2
	 *      p[4] = y2 - ycoord of point 2
	 *      p[5] = x3 - xcoord of point 3
	 *      p[6] = y3 - ycoord of point 3
	 *      p[7] = x4 - xcoord of point 4
	 *      p[8] = y4 - ycoord of point 4
	 *      rad  = radius of lip 
	 *	theta = angle (reference angle)
         *
         *
         *
         *                        y=c
         *               3------------------------4
         *               /
         *             /
         *           /
         *          /
         *         /
         *        /circular Arc
         *        |
         *        |         *(xc,yc)
         *        |
         *        \
         *         \
         *          \ 
         *           \ 
         *             \ 
         *               \ 
         *               2------------------------1
         *                        y=c
         *
         *
         */
        else
	if (model_id == 26)
        {
          dbl xx2, yy2, xx3, yy3, xcen , ycen, rad, theta;
          xx2=p[3];
          yy2=p[4];
          xx3=p[5];
          yy3=p[6];

          /**   compute angle of point on curve from (xcen,ycen) **/

          xcen = (xx2+xx3)/2;
          ycen = (yy2+yy3)/2.;
          rad = sqrt(SQUARE(yy2-yy3)+SQUARE(xx2-xx3))/2;
          theta = atan2(x2-ycen,x1-xcen);

          /**  use different f depending on theta  **/

          if(-M_PIE/2 <= theta && theta <= 0)
            {
              f = (x2-yy2);
            }
          else if(-M_PIE <= theta && theta <= -M_PIE/2)
            {
              f = SQUARE(x1-xcen)+SQUARE(x2-ycen)-SQUARE(rad);
            }
          else if(M_PIE/2 <= theta && theta <= M_PIE)
            {
              f = SQUARE(x1-xcen)+SQUARE(x2-ycen)-SQUARE(rad);
            }
          else if(0 <= theta && theta <= M_PIE/2)
            {
              f = x2-yy3;
            }

        }   /*  end of syringe with curved lip  */
	/*   multiple radiused lip corners
	 *      p[1] = number of corners
	 *      p[2] = theta0
	 *      p[3+3i] = xptl
	 *      p[4+3i] = ypt1
	 *      p[5+3i] = radius1
	 *      p[3+3n] = thetan
	 *	xcen = x-coord of radius center
	 *	ycen = y-coord of radius center
	 *	alpha = subtended half angle
	 */
	else 
	if (model_id == 27)
	{
	  dbl xpt[10],ypt[10],theta[10],rad[10],xcen[10],ycen[10],alpha[10];
	  dbl dist[10], rad_sign[10],th1, th1t;
          double theta_tmp=0, dist_min, theta1, theta2;
          int i,i_min,sign_bit,n_corner=0;
	  n_corner = ((int)p[1]);
	  sign_bit = ((int)p[2]);
          if(n_corner > 10)
          	{EH(-1,"too many corners in multiple radius lip\n");}
          i_min = 1;
          for(i=0 ; i<n_corner ; i++)
          	{
          	rad_sign[i] = 1.;
          	if(sign_bit & i_min) {  rad_sign[i] = -1.;}
          	i_min *= 2;
          	}
          theta[0] = p[3];
          for(i=0 ; i<n_corner ; i++)
          	{
          	xpt[i] = p[4+3*i];
          	ypt[i] = p[5+3*i];
          	rad[i] = p[6+3*i];
          	}
          theta[n_corner] = p[4+3*n_corner];
	  
	/*  slope of middle line		*/

          for(i=1 ; i<n_corner ; i++)
          	{
	  	theta[i] = atan2(ypt[i]-ypt[i-1],xpt[i]-xpt[i-1]);
          	}
          for(i=0 ; i<n_corner ; i++)
          	{
	  	theta[i] = theta[i] >= 0 ? theta[i] : theta[i] + 2*M_PIE;  
	  	theta[i] = theta[i] <= 2*M_PIE ? theta[i] : theta[i] - 2*M_PIE;  
          	}
          for(i=0 ; i<n_corner ; i++)
          	{
          	theta_tmp = theta[i]+M_PIE;
	  	theta_tmp = theta_tmp <= 2*M_PIE ? theta_tmp : theta_tmp-2*M_PIE;  
	  	alpha[i] = 0.5*(theta_tmp-theta[i+1]);
	  	alpha[i] = alpha[i] <= 0.5*M_PIE ? alpha[i] : M_PIE-alpha[i];  
	  	xcen[i] = xpt[i]+(rad_sign[i]*rad[i]/sin(alpha[i]))*cos(theta[i+1]+alpha[i]);
	  	ycen[i] = ypt[i]+(rad_sign[i]*rad[i]/sin(alpha[i]))*sin(theta[i+1]+alpha[i]);
          	}

	  /**   compute angle of point on curve from arc center **/
	  
          dist_min = 10000000.;   i_min = 0;
          for(i=0 ; i<n_corner ; i++)
          	{
          	dist[i] = sqrt(SQUARE(x2-ypt[i]) + SQUARE(x1-xpt[i]));
          	if(dist[i] < dist_min)
          		{dist_min = dist[i];  i_min = i;}
          	}
	  th1 = atan2(x2-ycen[i_min],x1-xcen[i_min]);
	  th1t = th1 > 0.0 ? th1 : th1 + 2*M_PIE;  

          if(rad_sign[i_min] > 0)
          	{ theta1 = theta[i_min+1];  theta2 = theta[i_min]+M_PIE;}
          else
          	{ theta1 = theta[i_min]+M_PIE;  theta2 = theta[i_min+1];}

	  /**  use different f depending on theta  **/
	  
if(x1>0.00001)fprintf(stderr,"user %d %g %g %g %g %g %g\n",i_min,x1,x2,xpt[i_min],ypt[i_min],th1,th1t);
	  if( (theta1 <= M_PIE && (theta1-0.5*M_PIE) <= th1 && th1 <= theta1)
	  || (theta1 > M_PIE && (theta1-0.5*M_PIE) <= th1t && th1t <= theta1))
	    {
	      f = (x2-ypt[i_min])*cos(theta1) -(x1-xpt[i_min])*sin(theta1);
if(x1>0.00001)fprintf(stderr,"line1  %g %g %g %g\n",f,theta1,th1,th1t);
	    }
	  else 
	  if( (theta2 <= 0.5*M_PIE && (theta2+0.5*M_PIE) >= th1 && th1 >= theta2)
	  || (theta2 > 0.5*M_PIE && (theta2+0.5*M_PIE) >= th1t && th1t >= theta2))
	    {
	      f = (x2-ypt[i_min])*cos(theta2) - (x1-xpt[i_min])*sin(theta2);
if(x1>0.00001)fprintf(stderr,"line2  %g %g %g %g\n",f,theta2,th1,th1t);
	    }
	  else
	    {
	      f = SQUARE(x1-xcen[i_min])+SQUARE(x2-ycen[i_min])-SQUARE(rad[i_min]);
if(x1>0.00001)fprintf(stderr,"radius  %g %g %g %g\n",f,rad[i_min],xcen[i_min],ycen[i_min]);
	    }
      
	}   /*  end of multiple radiused lip corners */
           /*  falling circle  */
	else 
	if (model_id == 30)
	{
	  dbl xcen1 , ycen1, xvelo, yvelo, radius;
          dbl xcenter, ycenter;
	  xcen1=p[1];
	  ycen1=p[2];
	  radius=p[3];
	  xvelo=p[4];
	  yvelo=p[5];
	  xcenter = xcen1 + xvelo*time;
	  ycenter = ycen1 + yvelo*time;
	      f = SQUARE(x1-xcenter)+SQUARE(x2-ycenter)-SQUARE(radius);
	}
           /*  moving circle  */
	else 
	if (model_id == 31)
	{
	  dbl xcen1 , ycen1, deltax, deltay, radius;
          dbl xcenter, ycenter;
	  xcen1=p[1];
	  ycen1=p[2];
	  radius=p[3];
	  deltax=p[4];
	  deltay=p[5];
	  xcenter = xcen1 + deltax;
	  ycenter = ycen1 + deltay;
	      f = SQUARE(x1-xcenter)+SQUARE(x2-ycenter)-SQUARE(radius);
	}
           /*  radiused upstream bead  */
	else 
	if (model_id == 32)
	{
	  dbl radius,dca,sca, lip_angle, alpha1, alpha2;
          dbl xcenter, ycenter, pos_dcl[3]={0,0,0},pos_scl[3]={0,0,0};
          int nset_dcl, nset_scl, nsp, k, j, i, dir;
	  nset_dcl = ((int)p[1]);
	  nset_scl = ((int)p[2]);
	  dca=M_PIE*p[3]/180.0;
	  sca=M_PIE*p[4]/180.0;
	  lip_angle=M_PIE*p[5]/180.0;
          alpha1 = 0.5*M_PIE - sca + lip_angle;
          alpha2 = dca - 0.5*M_PIE;
         
          /*  Get DCL and SCL coordinates  */
          nsp = match_nsid(nset_dcl);
          k = Proc_NS_List[Proc_NS_Pointers[nsp]];
          for (j = 0; j < Proc_NS_Count[nsp]; j++)
             {
               k = Proc_NS_List[Proc_NS_Pointers[nsp]+j];
               for (dir = 0; dir < pd->Num_Dim ; dir++)
                   {
                    i = Index_Solution (k, MESH_DISPLACEMENT1+dir, 0, 0, -1);
                    EH(i, "Could not resolve index_solution.");
/*                    pos_dcl[dir] = Coor[dir][k] + x[i+dir];*/
                   }
             }
          nsp = match_nsid(nset_scl);
          k = Proc_NS_List[Proc_NS_Pointers[nsp]];
          for (j = 0; j < Proc_NS_Count[nsp]; j++)
             {
               k = Proc_NS_List[Proc_NS_Pointers[nsp]+j];
               for (dir = 0; dir < pd->Num_Dim ; dir++)
                   {
                    i = Index_Solution (k, MESH_DISPLACEMENT1+dir, 0, 0, -1);
                    EH(i, "Could not resolve index_solution.");
/*                    pos_scl[dir] = Coor[dir][k] + x[i+dir];*/
                   }
             }

          /* compute upstream bead radius  */
          radius = ((pos_dcl[0]-pos_scl[0])*sin(-lip_angle)
                          +(pos_dcl[1]-pos_scl[1])*cos(-lip_angle))/
                   (sin(alpha1-lip_angle)-sin(alpha2-lip_angle));
	  xcenter = pos_dcl[0] + radius*cos(alpha2);
	  ycenter = pos_dcl[1] + radius*sin(alpha2);
	  f = SQUARE(x1-xcenter)+SQUARE(x2-ycenter)-SQUARE(radius);
	}
  else
    {
      EH(-1,"invalid model_id in user_bc");
    }

        return(f);

}
/*****************************************************************************/
dbl dfncd1(x1, x2,  x3, p, time)
     const dbl x1;
     const dbl x2;
     const dbl x3;
     const dbl p[];
     const dbl time;
{
/* for PI use M_PIE Constant from std.h include file. */
	int model_id;
	dbl dfdx1;

	model_id = ((int)p[0]);


	if(model_id == 1)
	{
	dbl rad,ampl,th0,wid,web_sp,arc_position,delta_rad;
        rad=p[1];
        ampl=p[2];
        web_sp=p[5];
        th0=p[3]+(web_sp/rad)*time;
        wid=p[4];
	arc_position = rad*(atan2(x2-p[7],x1-p[6])-th0)/wid;
	delta_rad = ampl*exp(-arc_position*arc_position);
        dfdx1=2.*(x1-p[6])-2.*(rad+delta_rad)*delta_rad*(-2.)
              *arc_position*(rad/wid)*
			(-(x2-p[7])/(SQUARE(x1-p[6])+SQUARE(x2-p[7])));
	}  /*  end of Gaussian bump if block */
	else if(model_id == 2)
	{
	dbl xp, slope, ysurfds=0, dsdx;
        dbl cube_dep, cube_lng, flat, radius, xveloc, zangle, xoffset;
        dbl xcenter, ycenter, roll_rad, arc_pos, ysurf=0, flat_bot, xref=0;

        double x_dist[4], theta;
        cube_dep = p[1];
        cube_lng = p[2];
        flat = p[3];
        flat_bot = p[4];
        radius = p[5];
        xveloc = p[6];
        zangle = p[7];
        xoffset = p[8];
        roll_rad = p[9];
        xcenter = p[10];
        ycenter = p[11];

        theta = atan2(cube_lng*0.5,cube_dep);
        x_dist[0] = radius/tan(0.5*theta+0.25*M_PIE);
        x_dist[1] = x_dist[0]*sin(theta);
        flat = MAX(flat,2*x_dist[0]);
        flat_bot = MAX(flat_bot,2*x_dist[0]);
        slope = 2.*cube_dep/cube_lng;
        arc_pos = roll_rad * (0.5*M_PIE - acos((x1-xcenter)
                /sqrt(SQUARE(x1-xcenter)+SQUARE(x2-ycenter))));

	dsdx=roll_rad*((x2-ycenter)/(SQUARE(x1-xcenter)+SQUARE(x2-ycenter)));

        xp = fmod(arc_pos-xveloc*time-zangle*x3-xoffset,cube_lng+flat+flat_bot);
        if(xp < 0.0){xp += cube_lng+flat+flat_bot;}

        if(xp >= (-0.5*flat+x_dist[0]) && xp <= (0.5*flat-x_dist[0]))
                { ysurfds = 0.0; }

        else if(xp >= (0.5*flat-x_dist[0]) && xp <= (0.5*flat+x_dist[1]))
                {
                xref = 0.5*flat-x_dist[0];
                ysurfds = (xref-xp)/sqrt(SQUARE(radius)-SQUARE(xp-xref));
                }

        else if(xp >= (0.5*flat+x_dist[1]) && xp <= (0.5*flat+0.5*cube_lng-x_dist[1]))
                {
                ysurfds =  -slope;
                }

        else if(xp >= (0.5*flat+0.5*cube_lng-x_dist[1]) &&
                xp <= (0.5*flat+0.5*cube_lng+x_dist[0]))
                {
                xref = 0.5*flat + 0.5*cube_lng + x_dist[0];
                ysurfds = (xp-xref)/sqrt(SQUARE(radius)-SQUARE(xp-xref));
                }

        else if(xp >= (0.5*flat+0.5*cube_lng+x_dist[0]) &&
                xp <= (0.5*flat+0.5*cube_lng+flat_bot-x_dist[0]))
                {
                ysurfds = 0.0;
                }

        else if(xp >= (0.5*flat+0.5*cube_lng+flat_bot-x_dist[0]) &&
                xp <= (0.5*flat+0.5*cube_lng+flat_bot+x_dist[1]))
                {
                xref = 0.5*flat + 0.5*cube_lng + flat_bot - x_dist[0];
                ysurfds = (xp-xref)/sqrt(SQUARE(radius)-SQUARE(xp-xref));
                }

        else if(xp >= (0.5*flat+0.5*cube_lng+flat_bot+x_dist[1]) &&
                xp <= (0.5*flat+cube_lng+flat_bot-x_dist[1]))
                {
                ysurfds = slope;
                }
        else if(xp >= (0.5*flat+cube_lng+flat_bot-x_dist[1]) &&
                xp <= (0.5*flat+cube_lng+flat_bot+x_dist[0]))
                {
                xref = 0.5*flat + cube_lng + flat_bot + x_dist[0];
                ysurfds = (xref-xp)/sqrt(SQUARE(radius)-SQUARE(xp-xref));
                }
        else if(xp >= (0.5*flat+cube_lng+flat_bot+x_dist[0]) &&
                xp <= (flat+cube_lng+flat_bot-x_dist[0]))
                {
                ysurfds = 0.0;
                }
	else if(xp >= (flat+cube_lng+flat_bot-x_dist[0]) && 
		xp <= (flat+cube_lng+flat_bot+x_dist[1]))
        	{ 
		xref = flat+cube_lng+flat_bot-x_dist[0];
                ysurfds = (xref-xp)/sqrt(SQUARE(radius)-SQUARE(xp-xref));
		}
        else
        	{
      		EH(-1,"invalid gravure user_bc position");
        	}

	dfdx1 = 2.*(x1-xcenter) - 2.*(roll_rad+ysurf)*ysurfds*dsdx;
	}
	else if(model_id == 3)
	{
	dbl xpt, ypt, theta1, theta2, rad, xcen , ycen, alpha;
	dbl theta;
        xpt=p[1];
        ypt=p[2];
        theta1=p[3];
        theta2=p[4];
        rad=p[5];

	alpha = 0.5*(theta2-theta1);
	xcen = xpt + (rad/sin(alpha))*cos(theta1+alpha);
	ycen = ypt + (rad/sin(alpha))*sin(theta1+alpha);

	/**   compute angle of point on curve from arc center **/

	  theta = atan2(x2-ycen,x1-xcen);
	  theta = theta > theta2-1.5*M_PIE ? theta : theta + 2*M_PIE;  

	/**  use different f depending on theta  **/

	if( (theta1-0.5*M_PIE) <= theta && theta <= (theta1+alpha))
		{
	dfdx1 =  -sin(theta1);
		}
	else if ( (theta1+alpha) <= theta && (theta - 0.5*M_PIE) <= theta2)
		{
	dfdx1 = -sin(theta2);
		}
	else
		{
	dfdx1 = 2.*(x1-xcen);
		}


	}   /*  end of radiused lip corner */

	else
	if (model_id == 4) {

	  dbl w1, dnum, dden, dnumdx, ddendx; 
	  const dbl *pn, *pd;
	  int i, n, xyz;
	  
	  n = abs(((int)p[1]));
	  xyz = ((int)p[2]);
	  
	  
	  pn = &p[3];
	  pd = &p[n+3];
	  
	  
	  switch (xyz)
		  {
		  case 21:
		  case 31:
			  w1=x1; 
	  dnum = pn[n-1];
	  dden = pd[n-2];
	  dnumdx = (n-1.)*pn[n-1];
	  ddendx = (n-2.)*pd[n-2];
	  for (i=n-2;i>=0;i--)
	    { dnum = pn[i]+w1*dnum; }
	  for (i=n-3;i>=0;i--)
	    { dden = pd[i]+w1*dden;
	      dnumdx = (i+1.)*pn[i+1]+w1*dnumdx;}
	  for (i=n-4;i>=0;i--)
	    { ddendx = (i+1.)*pd[i+1]+w1*ddendx; }
	  
	  dfdx1 = (dden*dnumdx - dnum*ddendx)/(dden*dden);
			  break;
	  
		  case 12:
		  case 13:
			  dfdx1 = -1.; break;
		  case 32:
		  case 23:
			  dfdx1 = 0.; break;
		  default:
		EH(-1,"invalid xyz in Pade function user_bc");

		}

	  } /* END OF FITTED CURVE */

	/*
	   CELL 
	   USES HYPERBOLIC TAN FUNCTIONS
	   
	   f(x) = y-(-d*0.5*(tanh(m*(x-s))+1.0)+d*0.5*tanh(m*(x-t))-1.0)+d)

	   HAVE TO PROVIDE COEFFICIENTS

	   p[0] = model_id
	   p[1] = d 
	   p[2] = m
	   p[3] = s
	   p[4] = width
	   p[5] = speed (transient analysis)
	   p[6] = coordinate direction (1 or 2)

	    */
	else
	if (model_id == 5) {
	    
	  dbl d, m, s, t, 
	      w1, w2;
          double coord;
          int dir;
	  
	  d = p[1];
	  m = p[2];
	  s = p[3]+p[5]*time;
	  t = s+p[4];
	  dir = ((int)p[6]);
	  if(dir == 1)	
		{
		coord = x1;
	  	w1 = cosh(m*(coord-s))*cosh(m*(coord-s));
	  	w2 = cosh(m*(coord-t))*cosh(m*(coord-t));
	        dfdx1 = d*m*0.5*(1.0/w1-1.0/w2); 
		}
	  else
		{dfdx1 = 1.0;}
	  } /* END OF CELL */

	    /*
		tapered surface for 3D manifolds
		driven by a Pade curve

	    p[0] = model_id
	    p[1] = theta, angle in xy-plane (degrees) 
	    p[2] = N, pade order
	    p[3 ... 2N+1] = pade coefficients

	    */
	else
	if (model_id == 6) {
	    
	  dbl th; 
	  
	  th     = M_PIE * p[1]/180.;
	  
	  /*
	  EVALUATE FUNCTION
	  */
	  
	  dfdx1 = -sin(th);

	  } /* END OF tapered surface */

	else
	if (model_id == 7)	{

	  dbl theta1, theta2, rad, xcen , ycen, alpha;
	  dbl theta;
	  dbl x, dnum, dden, padex, padey;
	  const dbl *pnx, *pdx, *pny, *pdy;
	  int i,nx, ny;
	  theta1=p[1];
	  theta2=p[2];
	  rad=p[3];
	  nx      = abs(((int)p[4]));
	  pnx = &p[5];
	  pdx = &p[nx+5];
	  ny      = abs(((int)p[nx+6]));
	  pny = &p[nx+7];
	  pdy = &p[nx+7+ny];
	  x = x3;

	  /*
	  EVALUATE FUNCTION
	  */
	  
	  dnum = pnx[nx-1];
	  dden = pdx[nx-2];
	  for (i=nx-2;i>=0;i--)
	    { dnum = pnx[i]+x*dnum; }
	  for (i=nx-3;i>=0;i--)
	    { dden = pdx[i]+x*dden; }
	  
	  padex = dnum/dden;

	  dnum = pny[ny-1];
	  dden = pdy[ny-2];
	  for (i=ny-2;i>=0;i--)
	    { dnum = pny[i]+x*dnum; }
	  for (i=ny-3;i>=0;i--)
	    { dden = pdy[i]+x*dden; }
	  
	  padey = dnum/dden;

	alpha = 0.5*(theta2-theta1);
	xcen = padex + (rad/sin(alpha))*cos(theta1+alpha);
	ycen = padey + (rad/sin(alpha))*sin(theta1+alpha);

	/**   compute angle of point on curve from arc center **/

	  theta = atan2(x2-ycen,x1-xcen);
	  theta = theta > theta2-1.5*M_PIE ? theta : theta + 2*M_PIE;  

	/**  use different f depending on theta  **/

	if( (theta1-0.5*M_PIE) <= theta && theta <= (theta1+alpha))
		{
	dfdx1 =  -sin(theta1);
		}
	else if ( (theta1+alpha) <= theta && (theta - 0.5*M_PIE) <= theta2)
		{
	dfdx1 = -sin(theta2);
		}
	else
		{
	dfdx1 = 2.*(x1-xcen);
		}
	  }

	    /*
		SINGLE V CELL

	    p[0] = model_id = 8
	    p[1] = x01
	    p[2] = y01 
	    p[3] = H
	    p[4] = L
	    p[5] = A
	    p[6] = R
	    p[7] = G

	    */
	else
	if (model_id == 8)	{

	  dbl 
	    x01, x02, x03, x04, x05, x06, x07,
	    y01, y02, y03, y04, y05,
	    RcA, RsA, tA, yR, yC, yL, 
	    H, A, R;

	x01 = p[1];
	y01 = p[2];
	H   = p[3];
	A   = p[5];
	R   = p[6];

	RcA = R*cos(A);
	RsA = R*sin(A);
	tA  = tan(A);

	y02 = y01+R-RsA;
	x02 = x01-RcA;

	y03 = y02;
	x03 = x01+RcA;
	y04 = y01+H-R+RsA;
	x04 = x02-(y04-y02)*tA;
	y05 = y04;
	x05 = x03+(y05-y03)*tA;
	x06 = x04-RcA;
	x07 = x05+RcA;

	yC = y01+R;
	yR = y05-RsA;
	yL = y04-RsA;

	  /*
	  EVALUATE FUNCTION
	  */
	  
	if (x2 <= x06)
	  {
	    dfdx1 = 1.0;
	  }
	else
	if ((x2 <= x04) && (x2 > x06))
	  {
	    dfdx1 = 2.0*(x1-yL);
	  }
	else
	if ((x2 <= x02) && (x2 > x04))
	  {
	    dfdx1 = (x04-x02);
	  }
	else
	if ((x2 <= x03) && (x2 > x02))
	  {
	    dfdx1 = 2.0*(x1-yC);
	  }
	else
	if ((x2 <= x05) && (x2 > x03))
	  {
	    dfdx1 = (x05-x03);
	  }
	else
	if ((x2 <= x07) && (x2 > x05))
	  {
	    dfdx1 = 2.0*(x1-yR);
	  }
	else
	if (x2 > x07)
	  {
	    dfdx1 = 1.0;
	  } }
	else
        if (model_id == 9)
        {
        dbl x0,ampl,speed,width,xpt;
          x0=p[1];
          ampl=p[3];
          speed=p[4];
          width=p[5];
          xpt = x0+speed*time;
  
         dfdx1= - ampl*(-2.*(x1-xpt)/SQUARE(width))*exp(-SQUARE((x1-xpt)/width));
  
        }   /*  end of Doug's Gaussian bump if block */
	else if(model_id == 10)
	{
	dbl xpt, ypt, theta1, theta2, rad, xcen , ycen, alpha;
	dbl theta;
        xpt=p[1];
        ypt=p[2];
        theta1=p[3];
        theta2=p[4];
        rad=p[5];

	alpha = 0.5*(theta2-theta1);
	xcen = xpt + (rad/sin(alpha))*cos(theta1+alpha);
	ycen = ypt + (rad/sin(alpha))*sin(theta1+alpha);

	/**   compute angle of point on curve from arc center **/

	  theta = atan2(x2-ycen,x1-xcen);
	  theta = theta > theta2-1.5*M_PIE ? theta : theta + 2*M_PIE;  

	/**  use different f depending on theta  **/

	if( (theta1-0.5*M_PIE) <= theta && theta <= (theta1+alpha))
		{
	dfdx1 =  cos(theta1);
		}
	else if ( (theta1+alpha) <= theta && (theta - 0.5*M_PIE) <= theta2)
		{
	dfdx1 = cos(theta2);
		}
	else
		{
	dfdx1 = -sin(theta);
		}

	}   /*  end of normal line to radiused lip corner */

	else 
	if (model_id == 11)
	{
	  dbl xpt1, ypt1, theta1, rad1, xcen1 , ycen1, alpha1;
	  dbl xpt2, ypt2, theta2, rad2, xcen2 , ycen2, alpha2;
	  dbl theta1m, theta2m, th1, th2, th2t;
	  xpt1=p[1];
	  ypt1=p[2];
	  theta1=p[3];
	  rad1=p[4];
	  xpt2=p[5];
	  ypt2=p[6];
	  theta2=p[7];
	  rad2=p[8];
	  
	/*  slope of middle line		*/

	  theta1m = atan2(ypt2-ypt1,xpt2-xpt1);
	  theta1m = theta1m > theta1 ? theta1m : theta1m + 2*M_PIE;  
	  theta2m = atan2(ypt1-ypt2,xpt1-xpt2);
	  alpha1 = 0.5*(theta1m-theta1);
	  alpha2 = 0.5*(theta2-theta2m);

	  xcen1 = xpt1 + (rad1/sin(alpha1))*cos(theta1+alpha1);
	  ycen1 = ypt1 + (rad1/sin(alpha1))*sin(theta1+alpha1);
	  xcen2 = xpt2 + (rad2/sin(alpha2))*cos(theta2m+alpha2);
	  ycen2 = ypt2 + (rad2/sin(alpha2))*sin(theta2m+alpha2);
	  
	  /**   compute angle of point on curve from arc center **/
	  
	  th1 = atan2(x2-ycen1,x1-xcen1);
	  th2 = atan2(x2-ycen2,x1-xcen2);
	  th2t = th2 > 0.0 ? th2 : th2 + 2*M_PIE;  

	  /**  use different f depending on theta  **/
	  
	  if( (theta1-0.5*M_PIE) <= th1 && th1 <= theta1)
	    {
	      dfdx1 =  -sin(theta1);
	    }
	  else if ( theta2 <= th2t && (th2t - 0.5*M_PIE) <= theta2)
	    {
	      dfdx1 = -sin(theta2);
	    }
	  else if ( theta2m <= (th1+0.5*M_PIE) && th1 <= (theta1-0.5*M_PIE))
	    {
	      dfdx1 = 2.*(x1-xcen1);
	    }
	  else if ( (theta2m-0.5*M_PIE) <= th2 && (th1+0.5*M_PIE) <= theta2m)
	    {
	      dfdx1 = -sin(theta1m);
	    }
	  else
	    {
	      dfdx1 = 2.*(x1-xcen2);
	    }
      
	}   /*  end of 2 radiused lip corners */
	else 
	if (model_id == 12)
	{
	  dbl xcen1;
	  xcen1=p[1];
		dfdx1 = 2.*(x1-xcen1);
	}
	else 
	if (model_id == 13)
	{
	  dbl xpt, ypt, theta1, theta2, rad, xcen , ycen, alpha;
	  dbl theta, xcirc, ycirc, dierad;
	  int iside;
	  xpt=p[1];
	  ypt=p[2];
	  theta1=p[3];
	  theta2=p[4];
	  rad=p[5];
	  dierad=p[6];
	  iside = ((int)p[7]);
	  
	  /**  find center of die face  **/

	if( iside == 1)
	  {
	  xcirc = xpt + dierad*sin(theta1);
	  ycirc = ypt - dierad*cos(theta1);
	  }	else	{
	  xcirc = xpt - dierad*sin(theta2);
	  ycirc = ypt + dierad*cos(theta2);
	  }

	  /**  find center of corner radius  **/

	if( iside == 1)
	  	{
		alpha = theta2 - asin((rad-dierad*cos(theta2-theta1))/(rad+dierad));
	  	}	else	{
		alpha = theta1 + asin((rad-dierad*cos(theta2-theta1))/(rad+dierad));
	  	}
		xcen = xcirc + (rad+dierad)*cos(alpha);
		ycen = ycirc + (rad+dierad)*sin(alpha);
	  
	  /**   compute angle of point on curve from arc center **/
	  
	  theta = atan2(x2-ycen,x1-xcen);
	  theta = theta > alpha-M_PIE ? theta : theta + 2*M_PIE;  
	  
	  /**  use different f depending on theta  **/
	  
	if( iside == 1)
	  {
	  if( (alpha-M_PIE) <= theta && theta <= theta1)
	    {
	      dfdx1 = 2.*(x1-xcirc);
	    }
	  else if ( theta2 <= theta && (theta - 0.5*M_PIE) <= theta2)
	    {
	      dfdx1 =  -sin(theta2);
	    }
	  else
	    {
	      dfdx1 = 2.*(x1-xcen);
	    }
	  }	else	
	  {
	  if( theta2 <= theta && theta <= (alpha+M_PIE))
	    {
	      dfdx1 = 2.*(x1-xcirc);
	    }
	  else if ( (theta1-0.5*M_PIE) <= theta && theta <= theta1)
	    {
	      dfdx1 = -sin(theta1);
	    }
	  else
	    {
	      dfdx1 = 2.*(x1-xcen);
	    }
	  }
      
	}   /*  end of radiused lip corner */
	/*
	   Needle collar based on a tanh modified plane equation:

		Ax + By +Cz + D + E(tanh(beta*(x-x0))+1)/2 = 0
	   USES HYPERBOLIC TAN FUNCTIONS
	    */
	else
	if (model_id == 14) {
	    
	  double A, E, beta, x0;
	  
	  A = p[1];
	  E = p[5];
	  beta = p[6];
	  x0 = p[7];

	  dfdx1 = A + E*0.5*beta/SQUARE(cosh(beta*(x1-x0)));
	  }
	/* micro-flexo-printing bottom surface - omega, roll_rad, gap */
	else
	if (model_id == 15) {
	    
	  double omega, t_offset, time1;
	  omega = p[1];
	  t_offset = p[4];

	  time1 = time + t_offset;

	  dfdx1 = sin(2.*omega*time1);
	  }
	/* micro-flexo-printing vertical surface - omega, roll_rad, gap, width */
	else
	if (model_id == 16) {
	    
	  double omega, t_offset, time1;
	  omega = p[1];
	  t_offset = p[5];

	  time1 = time + t_offset;

	  dfdx1 = - cos(2.*omega*time1);
	  }
	else 
	if (model_id == 17)
	{
	  dbl xpt1, ypt1, theta1, rad1, xcen1 , ycen1, alpha1;
	  dbl xpt2, ypt2, theta2, rad2, xcen2 , ycen2, alpha2;
	  dbl theta1m, theta2m, th1, th2, th2t;
	  double omega, roll_rad, gap, width, t_offset, time1, yc;
	  omega = p[1];
	  roll_rad = p[2];
	  gap = p[3];
	  width = p[4];
	  t_offset = p[5];

	  time1 = time + t_offset;
	  yc = gap + 2.*roll_rad*(1.-sin(omega*time1));
	  xpt1= -0.5*width*cos(2.*omega*time1);
	  ypt1= yc - 0.5*width*sin(2.*omega*time1);
	  theta1= omega*time1;
	  rad1=p[6];
	  xpt2= 0.5*width*cos(2.*omega*time1);
	  ypt2= yc + 0.5*width*sin(2.*omega*time1);
	  theta2= omega*time1;
	  rad2=p[7];
	  
	/*  slope of middle line		*/

	  theta1m = atan2(ypt2-ypt1,xpt2-xpt1);
	  theta1m = theta1m > theta1 ? theta1m : theta1m + 2*M_PIE;  
	  theta2m = atan2(ypt1-ypt2,xpt1-xpt2);
	  alpha1 = 0.5*(theta1m-theta1);
	  alpha2 = 0.5*(theta2-theta2m);

	  xcen1 = xpt1 + (rad1/sin(alpha1))*cos(theta1+alpha1);
	  ycen1 = ypt1 + (rad1/sin(alpha1))*sin(theta1+alpha1);
	  xcen2 = xpt2 + (rad2/sin(alpha2))*cos(theta2m+alpha2);
	  ycen2 = ypt2 + (rad2/sin(alpha2))*sin(theta2m+alpha2);
	  
	  /**   compute angle of point on curve from arc center **/
	  
	  th1 = atan2(x2-ycen1,x1-xcen1);
	  th2 = atan2(x2-ycen2,x1-xcen2);
	  th2t = th2 > 0.0 ? th2 : th2 + 2*M_PIE;  

	  /**  use different f depending on theta  **/
	  
	  if( (theta1-0.5*M_PIE) <= th1 && th1 <= theta1)
	    {
	      dfdx1 =  -sin(theta1);
	    }
	  else if ( theta2 <= th2t && (th2t - 0.5*M_PIE) <= theta2)
	    {
	      dfdx1 = -sin(theta2);
	    }
	  else if ( theta2m <= (th1+0.5*M_PIE) && th1 <= (theta1-0.5*M_PIE))
	    {
	      dfdx1 = 2.*(x1-xcen1);
	    }
	  else if ( (theta2m-0.5*M_PIE) <= th2 && (th1+0.5*M_PIE) <= theta2m)
	    {
	      dfdx1 = -sin(theta1m);
	    }
	  else
	    {
	      dfdx1 = 2.*(x1-xcen2);
	    }
      
	}   /*  end of microprinting (original) */
	else 
	if (model_id == 18)
	{
	  dbl xpt1, ypt1, theta1, rad1, xcen1 , ycen1, alpha1;
	  dbl xpt2, ypt2, theta2, rad2, xcen2 , ycen2, alpha2;
	  dbl theta1m, theta2m, th1, th2, th2t;
	  double omega, roll_rad, gap, width, t_offset, time1, yc, xc;
	  double y_displ, twist, angle, l_arm;
	  omega = p[1];
	  roll_rad = p[2];
	  gap = p[3];
	  width = p[4];
	  t_offset = p[5];
	  rad1=p[6];
	  rad2=p[7];
	  y_displ=p[8];
	  twist=p[9];
	  l_arm=p[10];

	  time1 = time + t_offset;
	  angle = omega*time1;
	  yc = gap + 2.*roll_rad*(1.-sin(angle));
/*  account for displacement	*/
	  yc += y_displ*sin(angle);
	  xc = y_displ*cos(angle);
/*  account for twist	*/
	  yc += l_arm*(sin(angle)-sin(angle+twist));
	  xc += l_arm*(cos(angle)-cos(angle+twist));
	  angle += twist;
	  xpt1= xc + 0.5*width*sin(angle);
	  ypt1= yc - 0.5*width*cos(angle);
	  theta1= angle;
	  xpt2= xc - 0.5*width*sin(angle);
	  ypt2= yc + 0.5*width*cos(angle);
	  theta2= angle;
	  
	/*  slope of middle line		*/

	  theta1m = atan2(ypt2-ypt1,xpt2-xpt1);
	  theta1m = theta1m > theta1 ? theta1m : theta1m + 2*M_PIE;  
	  theta2m = atan2(ypt1-ypt2,xpt1-xpt2);
	  alpha1 = 0.5*(theta1m-theta1);
	  alpha2 = 0.5*(theta2-theta2m);

	  xcen1 = xpt1 + (rad1/sin(alpha1))*cos(theta1+alpha1);
	  ycen1 = ypt1 + (rad1/sin(alpha1))*sin(theta1+alpha1);
	  xcen2 = xpt2 + (rad2/sin(alpha2))*cos(theta2m+alpha2);
	  ycen2 = ypt2 + (rad2/sin(alpha2))*sin(theta2m+alpha2);
	  
	  /**   compute angle of point on curve from arc center **/
	  
	  th1 = atan2(x2-ycen1,x1-xcen1);
	  th2 = atan2(x2-ycen2,x1-xcen2);
	  th2t = th2 > 0.0 ? th2 : th2 + 2*M_PIE;  

	  /**  use different f depending on theta  **/
	  
	  if( (theta1-0.5*M_PIE) <= th1 && th1 <= theta1)
	    {
	      dfdx1 = - sin(theta1);
	    }
	  else if ( theta2 <= th2t && (th2t - 0.5*M_PIE) <= theta2)
	    {
	      dfdx1 = - sin(theta2);
	    }
	  else if ( theta2m <= (th1+0.5*M_PIE) && th1 <= (theta1-0.5*M_PIE))
	    {
	      dfdx1 = 2.*(x1-xcen1);
	    }
	  else if ( (theta2m-0.5*M_PIE) <= th2 && (th1+0.5*M_PIE) <= theta2m)
	    {
	      dfdx1 = - sin(theta1m);
	    }
	  else
	    {
	      dfdx1 = 2.*(x1-xcen2);
	    }
      
	}   /*  end of micro-printing (twist) */
	else 
	if (model_id == 19)
	{
	  dbl xpt, ypt, theta1, theta2, rad, xcen , ycen, alpha;
	  dbl theta, tdc_delay;
	  double omega, roll_rad, gap, width, t_offset, time1, yc;
	  omega = p[1];
	  roll_rad = p[2];
	  gap = p[3];
	  width = p[4];
	  t_offset = p[5];
	  rad=p[6];
	  theta2=p[7];
	  tdc_delay=p[8];

	  time1 = time + t_offset;
	  if (time <= 0.0)	
		{ yc = gap + 2.*roll_rad*(1.-sin(omega*time1));}
	  else if(time <= tdc_delay)
		{ yc = gap ;}
	  else 
		{ time1 -= tdc_delay;
		 yc = gap + 2.*roll_rad*(1.-sin(omega*time1));
		}
	  xpt=yc;
	  ypt=0.5*width;
	  theta1= -0.5*M_PIE;
	  
	  alpha = 0.5*(theta2-theta1);
	  xcen = xpt + (rad/sin(alpha))*cos(theta1+alpha);
	  ycen = ypt + (rad/sin(alpha))*sin(theta1+alpha);
	  
	  /**   compute angle of point on curve from arc center **/
	  
	  theta = atan2(x2-ycen,x1-xcen);
	  theta = theta > theta2-1.5*M_PIE ? theta : theta + 2*M_PIE;  
	  
	  /**  use different f depending on theta  **/
	  
	  if( (theta1-0.5*M_PIE) <= theta && theta <= (theta1+alpha))
	    {
	      dfdx1 =  - sin(theta1);
	    }
	  else if ( (theta1+alpha) <= theta && (theta - 0.5*M_PIE) <= theta2)
	    {
	      dfdx1 = - sin(theta2);
	    }
	  else
	    {
	      dfdx1 = 2.*(x1-xcen);
	    }
      
	}   /*  end of micro-printing (axi-symmetric) */
	else 
	if (model_id == 20)
	{
	  double omega, roll_rad, t_offset, time1;
	  omega = p[1];
	  roll_rad = p[2];
	  t_offset = p[3];

	  time1 = time + t_offset;
	  dfdx1 = (x1-2.*roll_rad*cos(omega*time1)) + x1;
	}
	else 
	if (model_id == 21)
	{
	  double omega, roll_rad, t_offset, time1, x_offset;
	  omega = p[1];
	  roll_rad = p[2];
	  t_offset = p[3];
	  x_offset = p[4];

	  time1 = time + t_offset;
	  dfdx1 = sin(omega*time1+x_offset/roll_rad);
	}
	else 
	if (model_id == 22)
	{
	  dbl xpt1, ypt1, theta1, rad1, xcen1 , ycen1, alpha1;
	  dbl xpt2, ypt2, theta2, rad2, xcen2 , ycen2, alpha2;
	  dbl theta1m, theta2m, th1, th2, th2t;
	  double omega_l, omega_u, roll_rad_l, roll_rad_u, gap, width;
	  double t_offset, time1, yc, xc, alpha;
	  double y_displ, twist, angle_l, angle_u, l_arm;
	  omega_l = p[1];
	  omega_u = p[2];
	  roll_rad_l = p[3];
	  roll_rad_u = p[4];
	  gap = p[5];
	  width = p[6];
	  t_offset = p[7];
	  rad1=p[8];
	  rad2=p[9];
	  alpha=p[10];
	  y_displ=p[11];
	  twist=p[12];
	  l_arm=p[13];

	  time1 = time + t_offset;
	  angle_l = omega_l*time1;
	  angle_u = omega_u*time1 + 0.5*M_PIE*(1.-omega_u/omega_l);
	  yc = gap + roll_rad_l*(1.-sin(angle_l)) + roll_rad_u*(1.-sin(angle_u));
	  xc = -roll_rad_u*cos(angle_u) + roll_rad_l*cos(angle_l);
/*  account for displacement	*/
	  yc += y_displ*sin(angle_u);
	  xc += y_displ*cos(angle_u);
/*  account for twist	*/
	  yc += l_arm*(sin(angle_u)-sin(angle_u+twist));
	  xc += l_arm*(cos(angle_u)-cos(angle_u+twist));
	  angle_u += twist;
	  xpt1= xc + 0.5*width*sin(angle_u);
	  ypt1= yc - 0.5*width*cos(angle_u);
	  theta1= angle_u - alpha;
	  xpt2= xc - 0.5*width*sin(angle_u);
	  ypt2= yc + 0.5*width*cos(angle_u);
	  theta2= angle_u + alpha;
	  
	/*  slope of middle line		*/

	  theta1m = atan2(ypt2-ypt1,xpt2-xpt1);
	  theta1m = theta1m > theta1 ? theta1m : theta1m + 2*M_PIE;  
	  theta2m = atan2(ypt1-ypt2,xpt1-xpt2);
	  alpha1 = 0.5*(theta1m-theta1);
	  alpha2 = 0.5*(theta2-theta2m);

	  xcen1 = xpt1 + (rad1/sin(alpha1))*cos(theta1+alpha1);
	  ycen1 = ypt1 + (rad1/sin(alpha1))*sin(theta1+alpha1);
	  xcen2 = xpt2 + (rad2/sin(alpha2))*cos(theta2m+alpha2);
	  ycen2 = ypt2 + (rad2/sin(alpha2))*sin(theta2m+alpha2);
	  
	  /**   compute angle of point on curve from arc center **/
	  
	  th1 = atan2(x2-ycen1,x1-xcen1);
	  th2 = atan2(x2-ycen2,x1-xcen2);
	  th2t = th2 > 0.0 ? th2 : th2 + 2*M_PIE;  

	  /**  use different f depending on theta  **/
	  
	  if( (theta1-0.5*M_PIE) <= th1 && th1 <= theta1)
	    {
	      dfdx1 = - sin(theta1);
	    }
	  else if ( theta2 <= th2t && (th2t - 0.5*M_PIE) <= theta2)
	    {
	      dfdx1 = - sin(theta2);
	    }
	  else if ( theta2m <= (th1+0.5*M_PIE) && th1 <= (theta1-0.5*M_PIE))
	    {
	      dfdx1 = 2.*(x1-xcen1);
	    }
	  else if ( (theta2m-0.5*M_PIE) <= th2 && (th1+0.5*M_PIE) <= theta2m)
	    {
	      dfdx1 = - sin(theta1m);
	    }
	  else
	    {
	      dfdx1 = 2.*(x1-xcen2);
	    }
      
	}   /*  end of micro-printing (general) */
        else if(model_id == 25)
        {
          dbl xx1, yy1, xx2, yy2, xx3, yy3, xcen , alpha, theta, theta1;
          dbl xcen_new,ycen_new;
          xx1=p[1];
          yy1=p[2];
          xx2=p[3];
          yy2=p[4];
          xx3=p[5];
          yy3=p[6];
          xcen=p[7];

          /**   compute angle of point on curve from (xcen_new,ycen_new) **/

          xcen_new = (xx1+xx3)/2;
          ycen_new = (yy1+yy3)/2.;
          theta  = atan2(x2-ycen_new,x1-xcen_new);
          theta1 = atan2(yy2-ycen_new,xx2-xcen_new);
          alpha  = atan2(yy1-yy2,xx1-xx2);

          /**  use different f depending on theta  **/

          if( theta <= theta1)
            {
              dfdx1 = 2.*(x1-xcen);
            }
          else
            {
              dfdx1 = -sin(alpha );
            }

        }   /*  end of fluid bearing die */
        else
        if (model_id == 26)
        {
          dbl xx2, yy2, xx3, yy3, xcen , ycen, theta;
          xx2=p[3];
          yy2=p[4];
          xx3=p[5];
          yy3=p[6];

          /**   compute angle of point on curve from (xcen,ycen) **/

          xcen = (xx2+xx3)/2;
          ycen = (yy2+yy3)/2.;
          theta = atan2(x2-ycen,x1-xcen);

          /**  use different f depending on theta  **/

          if(-M_PIE/2 <= theta && theta <= 0)
            {
              dfdx1 = 0;
            }
          else if(-M_PIE <= theta && theta <= -M_PIE/2)
            {
              dfdx1 = 2*(x1-xcen);
            }
          else if(M_PIE/2 <= theta && theta <= M_PIE)
            {
              dfdx1 = 2*(x1-xcen);
            }
          else if(0 <= theta && theta <= M_PIE/2)
            {
              dfdx1 = 0;
            }

        }   /*  end of syringe with curved lip  */
	/*   multiple radiused lip corners
	 *      p[1] = number of corners
	 *      p[2] = theta0
	 *      p[3+3i] = xptl
	 *      p[4+3i] = ypt1
	 *      p[5+3i] = radius1
	 *      p[3+3n] = thetan
	 *	xcen = x-coord of radius center
	 *	ycen = y-coord of radius center
	 *	alpha = subtended half angle
	 */
	else 
	if (model_id == 27)
	{
	  dbl xpt[10],ypt[10],theta[10],rad[10],xcen[10],ycen[10],alpha[10];
	  dbl dist[10], rad_sign[10], th1, th1t;
          double theta_tmp=0, dist_min, theta1, theta2;
          int i,i_min,sign_bit,n_corner=0;
	  n_corner = ((int)p[1]);
	  sign_bit = ((int)p[2]);
          if(n_corner > 10)
          	{EH(-1,"too many corners in multiple radius lip\n");}
          i_min = 1;
          for(i=0 ; i<n_corner ; i++)
          	{
          	rad_sign[i] = 1.;
          	if(sign_bit & i_min) {  rad_sign[i] = -1.;}
          	i_min *= 2;
          	}
          theta[0] = p[3];
          for(i=0 ; i<n_corner ; i++)
          	{
          	xpt[i] = p[4+3*i];
          	ypt[i] = p[5+3*i];
          	rad[i] = p[6+3*i];
          	}
          theta[n_corner] = p[4+3*n_corner];
	  
	/*  slope of middle line		*/

          for(i=1 ; i<n_corner ; i++)
          	{
	  	theta[i] = atan2(ypt[i]-ypt[i-1],xpt[i]-xpt[i-1]);
          	}
          for(i=0 ; i<n_corner ; i++)
          	{
	  	theta[i] = theta[i] >= 0 ? theta[i] : theta[i] + 2*M_PIE;  
	  	theta[i] = theta[i] <= 2*M_PIE ? theta[i] : theta[i] - 2*M_PIE;  
          	}
          for(i=0 ; i<n_corner ; i++)
          	{
          	theta_tmp = theta[i]+M_PIE;
	  	theta_tmp = theta_tmp <= 2*M_PIE ? theta_tmp : theta_tmp-2*M_PIE;  
	  	alpha[i] = 0.5*(theta_tmp-theta[i+1]);
	  	alpha[i] = alpha[i] <= 0.5*M_PIE ? alpha[i] : M_PIE-alpha[i];  
	  	xcen[i] = xpt[i]+(rad_sign[i]*rad[i]/sin(alpha[i]))*cos(theta[i+1]+alpha[i]);
	  	ycen[i] = ypt[i]+(rad_sign[i]*rad[i]/sin(alpha[i]))*sin(theta[i+1]+alpha[i]);
          	}

	  /**   compute angle of point on curve from arc center **/
	  
          dist_min = 10000000.;   i_min = 0;
          for(i=0 ; i<n_corner ; i++)
          	{
          	dist[i] = sqrt(SQUARE(x2-ypt[i]) + SQUARE(x1-xpt[i]));
          	if(dist[i] < dist_min)
          		{dist_min = dist[i];  i_min = i;}
          	}
	  th1 = atan2(x2-ycen[i_min],x1-xcen[i_min]);
	  th1t = th1 > 0.0 ? th1 : th1 + 2*M_PIE;  

          if(rad_sign[i_min] > 0)
          	{ theta1 = theta[i_min+1];  theta2 = theta[i_min]+M_PIE;}
          else
          	{ theta1 = theta[i_min]+M_PIE;  theta2 = theta[i_min+1];}

	  /**  use different f depending on theta  **/
	  
	  if( (theta1 <= M_PIE && (theta1-0.5*M_PIE) <= th1 && th1 <= theta1)
	  || (theta1 > M_PIE && (theta1-0.5*M_PIE) <= th1t && th1t <= theta1))
	    {
	      dfdx1 = -sin(theta1);
	    }
	  else 
	  if( (theta2 <= 0.5*M_PIE && (theta2+0.5*M_PIE) >= th1 && th1 >= theta2)
	  || (theta2 > 0.5*M_PIE && (theta2+0.5*M_PIE) >= th1t && th1t >= theta2))
	    {
	      dfdx1 = -sin(theta2);
	    }
	  else
	    {
	      dfdx1 = 2.*(x1-xcen[i_min]);
if(x1>0.00001)fprintf(stderr,"radius  %g %g %g %g\n",dfdx1,rad[i_min],xcen[i_min],ycen[i_min]);
	    }
      
	}   /*  end of multiple radiused lip corners */
           /*  falling circle  */
	else 
	if (model_id == 30)
	{
	  dbl xcen1 , xvelo;
          dbl xcenter;
	  xcen1=p[1];
	  xvelo=p[4];
	  xcenter = xcen1 + xvelo*time;
	      dfdx1 = 2.*(x1-xcenter);
	}
	else 
	if (model_id == 31)
	{
	  dbl xcen1 , deltax;
          dbl xcenter;
	  xcen1=p[1];
	  deltax=p[4];
	  xcenter = xcen1 + deltax;
	      dfdx1 = 2.*(x1-xcenter);
	}
           /*  radiused upstream bead  */
	else 
	if (model_id == 32)
	{
	  dbl radius,dca,sca, lip_angle, alpha1, alpha2;
          dbl xcenter, pos_dcl[3]={0,0,0},pos_scl[3]={0,0,0};
          int nset_dcl, nset_scl, nsp, k, j, i, dir;
	  nset_dcl = ((int)p[1]);
	  nset_scl = ((int)p[2]);
	  dca=M_PIE*p[3]/180.0;
	  sca=M_PIE*p[4]/180.0;
	  lip_angle=M_PIE*p[5]/180.0;
          alpha1 = 0.5*M_PIE - sca + lip_angle;
          alpha2 = dca - 0.5*M_PIE;
         
          /*  Get DCL and SCL coordinates  */
          nsp = match_nsid(nset_dcl);
          k = Proc_NS_List[Proc_NS_Pointers[nsp]];
          for (j = 0; j < Proc_NS_Count[nsp]; j++)
             {
               k = Proc_NS_List[Proc_NS_Pointers[nsp]+j];
               for (dir = 0; dir < pd->Num_Dim ; dir++)
                   {
                    i = Index_Solution (k, MESH_DISPLACEMENT1+dir, 0, 0, -1);
                    EH(i, "Could not resolve index_solution.");
                  /*  pos_dcl[dir] = Coor[dir][k] + x[i+dir];  */
                   }
             }
          nsp = match_nsid(nset_scl);
          k = Proc_NS_List[Proc_NS_Pointers[nsp]];
          for (j = 0; j < Proc_NS_Count[nsp]; j++)
             {
               k = Proc_NS_List[Proc_NS_Pointers[nsp]+j];
               for (dir = 0; dir < pd->Num_Dim ; dir++)
                   {
                    i = Index_Solution (k, MESH_DISPLACEMENT1+dir, 0, 0, -1);
                    EH(i, "Could not resolve index_solution.");
/*                    pos_scl[dir] = Coor[dir][k] + x[i+dir];*/
                   }
             }

          /* compute upstream bead radius  */
          radius = ((pos_dcl[0]-pos_scl[0])*sin(-lip_angle)
                          +(pos_dcl[1]-pos_scl[1])*cos(-lip_angle))/
                   (sin(alpha1-lip_angle)-sin(alpha2-lip_angle));

	  xcenter = pos_dcl[0] + radius*cos(alpha2);
	  dfdx1 = 2.0*(x1-xcenter);
	}
        else
	{
      	EH(-1,"invalid model_id in user_bc");
	}
        return(dfdx1);
}
/*****************************************************************************/
dbl dfncd2(x1, x2,  x3, p, time)
     const dbl x1;
     const dbl x2;
     const dbl x3;
     const dbl p[];
     const dbl time;
{
/* for PI use M_PIE Constant from std.h include file. */
	int model_id;
	dbl dfdx2;

	model_id = ((int)p[0]);

if(af->Assemble_LSA_Mass_Matrix)  return 0 ;

	if(model_id == 1)
	{
	dbl rad,ampl,th0,wid,web_sp, arc_position, delta_rad;
        rad=p[1];
        ampl=p[2];
        web_sp=p[5];
        th0=p[3]+(web_sp/rad)*time;
        wid=p[4];
	arc_position = rad*(atan2(x2-p[7],x1-p[6])-th0)/wid;
	delta_rad = ampl*exp(-arc_position*arc_position);

	dfdx2=2.*(x2-p[7])-2.*(rad+delta_rad)*delta_rad*(-2.)
	  *arc_position*(rad/wid)*((x1-p[6])/(SQUARE(x1-p[6])+SQUARE(x2-p[7])));
	}  /*  end of Gaussian bump if block */
	else 
	if(model_id == 2)
	{
        dbl xp, slope, ysurfds=0, dsdy;
        dbl cube_dep, cube_lng, flat, radius, xveloc, zangle, xoffset;
        dbl xcenter, ycenter, roll_rad, arc_pos, ysurf=0, flat_bot, xref=0;

        double x_dist[4], theta;
        cube_dep = p[1];
        cube_lng = p[2];
        flat = p[3];
        flat_bot = p[4];
        radius = p[5];
        xveloc = p[6];
        zangle = p[7];
        xoffset = p[8];
        roll_rad = p[9];
        xcenter = p[10];
        ycenter = p[11];

        theta = atan2(cube_lng*0.5,cube_dep);
        x_dist[0] = radius/tan(0.5*theta+0.25*M_PIE);
        x_dist[1] = x_dist[0]*sin(theta);
        flat = MAX(flat,2*x_dist[0]);
        flat_bot = MAX(flat_bot,2*x_dist[0]);
        slope = 2.*cube_dep/cube_lng;
        arc_pos = roll_rad * (0.5*M_PIE - acos((x1-xcenter)
                /sqrt(SQUARE(x1-xcenter)+SQUARE(x2-ycenter))));

	dsdy=xcenter-x1;

        xp = fmod(arc_pos-xveloc*time-zangle*x3-xoffset,cube_lng+flat+flat_bot);
        if(xp < 0.0){xp += cube_lng+flat+flat_bot;}

        if(xp >= (-0.5*flat+x_dist[0]) && xp <= (0.5*flat-x_dist[0]))
                { ysurfds = 0.0; }

        else if(xp >= (0.5*flat-x_dist[0]) && xp <= (0.5*flat+x_dist[1]))
                {
                xref = 0.5*flat-x_dist[0];
                ysurfds = (xref-xp)/sqrt(SQUARE(radius)-SQUARE(xp-xref));
                }

        else if(xp >= (0.5*flat+x_dist[1]) && xp <= (0.5*flat+0.5*cube_lng-x_dist[1]))
                {
                ysurfds =  -slope;
                }

        else if(xp >= (0.5*flat+0.5*cube_lng-x_dist[1]) &&
                xp <= (0.5*flat+0.5*cube_lng+x_dist[0]))
                {
                xref = 0.5*flat + 0.5*cube_lng + x_dist[0];
                ysurfds = (xp-xref)/sqrt(SQUARE(radius)-SQUARE(xp-xref));
                }

        else if(xp >= (0.5*flat+0.5*cube_lng+x_dist[0]) &&
                xp <= (0.5*flat+0.5*cube_lng+flat_bot-x_dist[0]))
                {
                ysurfds = 0.0;
                }

        else if(xp >= (0.5*flat+0.5*cube_lng+flat_bot-x_dist[0]) &&
                xp <= (0.5*flat+0.5*cube_lng+flat_bot+x_dist[1]))
                {
                xref = 0.5*flat + 0.5*cube_lng + flat_bot - x_dist[0];
                ysurfds = (xp-xref)/sqrt(SQUARE(radius)-SQUARE(xp-xref));
                }

        else if(xp >= (0.5*flat+0.5*cube_lng+flat_bot+x_dist[1]) &&
                xp <= (0.5*flat+cube_lng+flat_bot-x_dist[1]))
                {
                ysurfds = slope;
                }
        else if(xp >= (0.5*flat+cube_lng+flat_bot-x_dist[1]) &&
                xp <= (0.5*flat+cube_lng+flat_bot+x_dist[0]))
                {
                xref = 0.5*flat + cube_lng + flat_bot + x_dist[0];
                ysurfds = (xref-xp)/sqrt(SQUARE(radius)-SQUARE(xp-xref));
                }
        else if(xp >= (0.5*flat+cube_lng+flat_bot+x_dist[0]) &&
                xp <= (flat+cube_lng+flat_bot-x_dist[0]))
                {
                ysurfds = 0.0;
                }
	else if(xp >= (flat+cube_lng+flat_bot-x_dist[0]) && 
		xp <= (flat+cube_lng+flat_bot+x_dist[1]))
        	{ 
		xref = flat+cube_lng+flat_bot-x_dist[0];
                ysurfds = (xref-xp)/sqrt(SQUARE(radius)-SQUARE(xp-xref));
		}
        else
        	{
      		EH(-1,"invalid gravure user_bc position");
        	}

	dfdx2 = 2.*(x2-ycenter) - 2.*(roll_rad+ysurf)*ysurfds*dsdy;

	}  /*   end of gravure if block */
	else if(model_id == 3)
	{
	dbl xpt, ypt, theta1, theta2, rad, xcen , ycen, alpha;
	dbl theta;
	xpt=p[1];
	ypt=p[2];
	theta1=p[3];
	theta2=p[4];
	rad=p[5];
	
	alpha = 0.5*(theta2-theta1);
	xcen = xpt + (rad/sin(alpha))*cos(theta1+alpha);
	ycen = ypt + (rad/sin(alpha))*sin(theta1+alpha);
	
	/**   compute angle of point on curve from arc center **/
	
	  theta = atan2(x2-ycen,x1-xcen);
	  theta = theta > theta2-1.5*M_PIE ? theta : theta + 2*M_PIE;  
	
	/**  use different f depending on theta  **/
	
	if( (theta1-0.5*M_PIE) <= theta && theta <= (theta1+alpha))
	  {
	    dfdx2 = cos(theta1);
	  }
	else if ( (theta1+alpha) <= theta && (theta - 0.5*M_PIE) <= theta2)
	  {
	    dfdx2 = cos(theta2);
	  }
	else
	  {
	    dfdx2 = 2.*(x2-ycen);
	  }
      
      }   /*  end of radiused lip corner */

	else
	if (model_id == 4) {

	  dbl w1, dnum, dden, dnumdx, ddendx; 
	  const dbl *pn, *pd;
	  int i, n, xyz;
	  
	  n = abs(((int)p[1]));
	  xyz = ((int)p[2]);
	  
	  
	  pn = &p[3];
	  pd = &p[n+3];
	  
	  
	  switch (xyz)
		  {
		  case 12:
		  case 32:
			  w1=x2; 
	  dnum = pn[n-1];
	  dden = pd[n-2];
	  dnumdx = (n-1.)*pn[n-1];
	  ddendx = (n-2.)*pd[n-2];
	  for (i=n-2;i>=0;i--)
	    { dnum = pn[i]+w1*dnum; }
	  for (i=n-3;i>=0;i--)
	    { dden = pd[i]+w1*dden;
	      dnumdx = (i+1.)*pn[i+1]+w1*dnumdx;}
	  for (i=n-4;i>=0;i--)
	    { ddendx = (i+1.)*pd[i+1]+w1*ddendx; }
	  
	  dfdx2 = (dden*dnumdx - dnum*ddendx)/(dden*dden);
			  break;
	  
		  case 21:
		  case 23:
			  dfdx2 = -1.; break;
		  case 31:
		  case 13:
			  dfdx2 = 0.; break;
		  default:
		EH(-1,"invalid xyz in Pade function user_bc");

		}
	   } /* END OF FITTED CURVE */

	/*
	   CELL 
	   USES HYPERBOLIC TAN FUNCTIONS
	   
	   f(x) = y-(-d*0.5*(tanh(m*(x-s))+1.0)+d*0.5*tanh(m*(x-t))-1.0)+d)

	   HAVE TO PROVIDE COEFFICIENTS

	   p[0] = model_id
	   p[1] = d 
	   p[2] = m
	   p[3] = s
	   p[4] = width
	   p[5] = speed (transient analysis)
	   p[6] = coordinate direction (1 or 2)

	    */
	else
	if (model_id == 5) {
	    
	  double d, m, s, t;
          double coord, w1, w2;
          int dir;
	  
	  d = p[1];
	  m = p[2];
	  s = p[3]+p[5]*time;
	  t = s+p[4];
	  dir = ((int)p[6]);
	  if(dir == 1)	
		{dfdx2 = 1.0;}
	  else
		{
		coord = x2;
	  	w1 = cosh(m*(coord-s))*cosh(m*(coord-s));
	  	w2 = cosh(m*(coord-t))*cosh(m*(coord-t));
	        dfdx2 = d*m*0.5*(1.0/w1-1.0/w2); 
		}
	  } /* END OF CELL */
	    /*
		tapered surface for 3D manifolds
		driven by a Pade curve

	    p[0] = model_id
	    p[1] = theta, angle in xy-plane (degrees) 
	    p[2] = N, pade order
	    p[3 ... 2N+1] = pade coefficients

	    */
	else
	if (model_id == 6) {
	    
	  dbl th; 
	  
	  th     = M_PIE * p[1]/180.;
	  /*
	  EVALUATE FUNCTION
	  */
	  
	  dfdx2 = cos(th);

	  } /* END OF tapered surface */

	else
	if (model_id == 7)	{

	  dbl theta1, theta2, rad, xcen , ycen, alpha;
	  dbl theta;
	  dbl x, dnum, dden, padex, padey;
	  const dbl *pnx, *pdx, *pny, *pdy;
	  int i,nx, ny;
	  theta1=p[1];
	  theta2=p[2];
	  rad=p[3];
	  nx      = abs(((int)p[4]));
	  pnx = &p[5];
	  pdx = &p[nx+5];
	  ny      = abs(((int)p[nx+6]));
	  pny = &p[nx+7];
	  pdy = &p[nx+7+ny];
	  x = x3;

	  /*
	  EVALUATE FUNCTION
	  */
	  
	  dnum = pnx[nx-1];
	  dden = pdx[nx-2];
	  for (i=nx-2;i>=0;i--)
	    { dnum = pnx[i]+x*dnum; }
	  for (i=nx-3;i>=0;i--)
	    { dden = pdx[i]+x*dden; }
	  
	  padex = dnum/dden;

	  dnum = pny[ny-1];
	  dden = pdy[ny-2];
	  for (i=ny-2;i>=0;i--)
	    { dnum = pny[i]+x*dnum; }
	  for (i=ny-3;i>=0;i--)
	    { dden = pdy[i]+x*dden; }
	  
	  padey = dnum/dden;

	alpha = 0.5*(theta2-theta1);
	xcen = padex + (rad/sin(alpha))*cos(theta1+alpha);
	ycen = padey + (rad/sin(alpha))*sin(theta1+alpha);

	/**   compute angle of point on curve from arc center **/

	  theta = atan2(x2-ycen,x1-xcen);
	  theta = theta > theta2-1.5*M_PIE ? theta : theta + 2*M_PIE;  

	/**  use different f depending on theta  **/

	if( (theta1-0.5*M_PIE) <= theta && theta <= (theta1+alpha))
		{
	dfdx2 =  cos(theta1);
		}
	else if ( (theta1+alpha) <= theta && (theta - 0.5*M_PIE) <= theta2)
		{
	dfdx2 = cos(theta2);
		}
	else
		{
	dfdx2 = 2.*(x2-ycen);
		}
	  }


	    /*
		SINGLE V CELL

	    p[0] = model_id = 8
	    p[1] = x01
	    p[2] = y01 
	    p[3] = H
	    p[4] = L
	    p[5] = A
	    p[6] = R
	    p[7] = G

	    */
	else
	if (model_id == 8)	{

	  dbl 
	    x01, x02, x03, x04, x05, x06, x07, y01, y02, y03, y04, y05, 
	    RcA, RsA, tA, xR, xC, xL, H, A, R;

	x01 = p[1];
	y01 = p[2];
	H   = p[3];
	A   = p[5];
	R   = p[6];

	RcA = R*cos(A);
	RsA = R*sin(A);
	tA  = tan(A);

	y02 = y01+R-RsA;
	x02 = x01-RcA;

	y03 = y02;
	x03 = x01+RcA;
	y04 = y01+H-R+RsA;
	x04 = x02-(y04-y02)*tA;
	y05 = y04;
	x05 = x03+(y05-y03)*tA;
	x06 = x04-RcA;
	x07 = x05+RcA;

	xC = x01;
	xR = x05+RcA;
	xL = x04-RcA;

	  /*
	  EVALUATE FUNCTION
	  */
	  
	if (x2 <= x06)
	  {
	    dfdx2 = 0.0;
	  }
	else
	if ((x2 <= x04) && (x2 > x06))
	  {
	    dfdx2 = 2.0*(x2-xL);
	  }
	else
	if ((x2 <= x02) && (x2 > x04))
	  {
	    dfdx2 = -(y04-y02);
	  }
	else
	if ((x2 <= x03) && (x2 > x02))
	  {
	    dfdx2 = 2.0*(x2-xC);
	  }
	else
	if ((x2 <= x05) && (x2 > x03))
	  {
	    dfdx2 = -(y05-y03);
	  }
	else
	if ((x2 <= x07) && (x2 > x05))
	  {
	    dfdx2 = 2.0*(x2-xR);
	  }
	else
	if (x2 > x07)
	  {
	    dfdx2 = 0.0;
	  } }
	else
        if (model_id == 9)
        {
         dfdx2=1.;
  
        }   /*  end of Doug's Gaussian bump if block */
	else if(model_id == 10)
	{
	dbl xpt, ypt, theta1, theta2, rad, xcen , ycen, alpha;
	dbl theta;
	xpt=p[1];
	ypt=p[2];
	theta1=p[3];
	theta2=p[4];
	rad=p[5];
	
	alpha = 0.5*(theta2-theta1);
	xcen = xpt + (rad/sin(alpha))*cos(theta1+alpha);
	ycen = ypt + (rad/sin(alpha))*sin(theta1+alpha);
	
	/**   compute angle of point on curve from arc center **/
	
	  theta = atan2(x2-ycen,x1-xcen);
	  theta = theta > theta2-1.5*M_PIE ? theta : theta + 2*M_PIE;  
	
	/**  use different f depending on theta  **/
	
	if( (theta1-0.5*M_PIE) <= theta && theta <= (theta1+alpha))
	  {
	    dfdx2 = sin(theta1);
	  }
	else if ( (theta1+alpha) <= theta && (theta - 0.5*M_PIE) <= theta2)
	  {
	    dfdx2 = sin(theta2);
	  }
	else
	  {
	    dfdx2 = cos(theta);
	  }
      
      }   /*  end of radiused lip corner */

	else 
	if (model_id == 11)
	{
	  dbl xpt1, ypt1, theta1, rad1, xcen1 , ycen1, alpha1;
	  dbl xpt2, ypt2, theta2, rad2, xcen2 , ycen2, alpha2;
	  dbl theta1m, theta2m, th1, th2, th2t;
	  xpt1=p[1];
	  ypt1=p[2];
	  theta1=p[3];
	  rad1=p[4];
	  xpt2=p[5];
	  ypt2=p[6];
	  theta2=p[7];
	  rad2=p[8];
	  
	/*  slope of middle line		*/

	  theta1m = atan2(ypt2-ypt1,xpt2-xpt1);
	  theta1m = theta1m > theta1 ? theta1m : theta1m + 2*M_PIE;  
	  theta2m = atan2(ypt1-ypt2,xpt1-xpt2);
	  alpha1 = 0.5*(theta1m-theta1);
	  alpha2 = 0.5*(theta2-theta2m);

	  xcen1 = xpt1 + (rad1/sin(alpha1))*cos(theta1+alpha1);
	  ycen1 = ypt1 + (rad1/sin(alpha1))*sin(theta1+alpha1);
	  xcen2 = xpt2 + (rad2/sin(alpha2))*cos(theta2m+alpha2);
	  ycen2 = ypt2 + (rad2/sin(alpha2))*sin(theta2m+alpha2);
	  
	  /**   compute angle of point on curve from arc center **/
	  
	  th1 = atan2(x2-ycen1,x1-xcen1);
	  th2 = atan2(x2-ycen2,x1-xcen2);
	  th2t = th2 > 0.0 ? th2 : th2 + 2*M_PIE;  

	  /**  use different f depending on theta  **/
	  
	  if( (theta1-0.5*M_PIE) <= th1 && th1 <= theta1)
	    {
	      dfdx2 = cos(theta1);
	    }
	  else if ( theta2 <= th2t && (th2t - 0.5*M_PIE) <= theta2)
	    {
	      dfdx2 = cos(theta2);
	    }
	  else if ( theta2m <= (th1+0.5*M_PIE) && th1 <= (theta1-0.5*M_PIE))
	    {
	      dfdx2 = 2.*(x2-ycen1);
	    }
	  else if ( (theta2m-0.5*M_PIE) <= th2 && (th1+0.5*M_PIE) <= theta2m)
	    {
	      dfdx2 = cos(theta1m);
	    }
	  else
	    {
	      dfdx2 = 2.*(x2-ycen2);
	    }
      
	}   /*  end of 2 radiused lip corners */
	else 
	if (model_id == 12)
	{
	  dbl ycen1;
	  ycen1=p[2];
		dfdx2 = 2.*(x2-ycen1);
	}
	else 
	if (model_id == 13)
	{
	  dbl xpt, ypt, theta1, theta2, rad, xcen , ycen, alpha;
	  dbl theta, xcirc, ycirc, dierad;
	  int iside;
	  xpt=p[1];
	  ypt=p[2];
	  theta1=p[3];
	  theta2=p[4];
	  rad=p[5];
	  dierad=p[6];
	  iside = ((int)p[7]);
	  
	  /**  find center of die face  **/

	if( iside == 1)
	  {
	  xcirc = xpt + dierad*sin(theta1);
	  ycirc = ypt - dierad*cos(theta1);
	  }	else	{
	  xcirc = xpt - dierad*sin(theta2);
	  ycirc = ypt + dierad*cos(theta2);
	  }

	  /**  find center of corner radius  **/

	if( iside == 1)
	  	{
		alpha = theta2 - asin((rad-dierad*cos(theta2-theta1))/(rad+dierad));
	  	}	else	{
		alpha = theta1 + asin((rad-dierad*cos(theta2-theta1))/(rad+dierad));
	  	}
		xcen = xcirc + (rad+dierad)*cos(alpha);
		ycen = ycirc + (rad+dierad)*sin(alpha);
	  
	  /**   compute angle of point on curve from arc center **/
	  
	  theta = atan2(x2-ycen,x1-xcen);
	  theta = theta > alpha-M_PIE ? theta : theta + 2*M_PIE;  
	  
	  /**  use different f depending on theta  **/
	  
	if( iside == 1)
	  {
	  if( (alpha-M_PIE) <= theta && theta <= theta1)
	    {
	      dfdx2 = 2.*(x2-ycirc);
	    }
	  else if ( theta2 <= theta && (theta - 0.5*M_PIE) <= theta2)
	    {
	      dfdx2 = cos(theta2);
	    }
	  else
	    {
	      dfdx2 = 2.*(x2-ycen);
	    }
	  }	else	
	  {
	  if( theta2 <= theta && theta <= (alpha+M_PIE))
	    {
	      dfdx2 = 2.*(x2-ycirc);
	    }
	  else if ( (theta1-0.5*M_PIE) <= theta && theta <= theta1)
	    {
	      dfdx2 = cos(theta1);
	    }
	  else
	    {
	      dfdx2 = 2.*(x2-ycen);
	    }
	  }
      
	}   /*  end of radiused lip corner */
	/*
	   Needle collar based on a tanh modified plane equation:

		Ax + By +Cz + D + E(tanh(beta*(x-x0))+1)/2 = 0
	   USES HYPERBOLIC TAN FUNCTIONS
	    */
	else
	if (model_id == 14) {
	    
	  dfdx2 = p[2];
	  }
	/* micro-flexo-printing bottom surface - omega, roll_rad, gap */
	else
	if (model_id == 15) {
	    
	  double omega, t_offset, time1;
	  omega = p[1];
	  t_offset = p[4];

	  time1 = time + t_offset;
	  dfdx2 =  - cos(2.*omega*time1);
	  }
	/* micro-flexo-printing vertical surface - omega, roll_rad, gap, width */
	else
	if (model_id == 16) {
	    
	  double omega, t_offset, time1;
	  omega = p[1];
 	  t_offset = p[5];

	  time1 = time + t_offset;
	  dfdx2 = - sin(2.*omega*time1);
	  }
	else 
	if (model_id == 17)
	{
	  dbl xpt1, ypt1, theta1, rad1, xcen1 , ycen1, alpha1;
	  dbl xpt2, ypt2, theta2, rad2, xcen2 , ycen2, alpha2;
	  dbl theta1m, theta2m, th1, th2, th2t;
	  double omega, roll_rad, gap, width, t_offset, time1, yc;
	  omega = p[1];
	  roll_rad = p[2];
	  gap = p[3];
	  width = p[4];
	  t_offset = p[5];

	  time1 = time + t_offset;
	  yc = gap + 2.*roll_rad*(1.-sin(omega*time1));
	  xpt1= -0.5*width*cos(2.*omega*time1);
	  ypt1= yc - 0.5*width*sin(2.*omega*time1);
	  theta1= omega*time1;
	  rad1=p[6];
	  xpt2= 0.5*width*cos(2.*omega*time1);
	  ypt2= yc + 0.5*width*sin(2.*omega*time1);
	  theta2= omega*time1;
	  rad2=p[7];
	  
	/*  slope of middle line		*/

	  theta1m = atan2(ypt2-ypt1,xpt2-xpt1);
	  theta1m = theta1m > theta1 ? theta1m : theta1m + 2*M_PIE;  
	  theta2m = atan2(ypt1-ypt2,xpt1-xpt2);
	  alpha1 = 0.5*(theta1m-theta1);
	  alpha2 = 0.5*(theta2-theta2m);

	  xcen1 = xpt1 + (rad1/sin(alpha1))*cos(theta1+alpha1);
	  ycen1 = ypt1 + (rad1/sin(alpha1))*sin(theta1+alpha1);
	  xcen2 = xpt2 + (rad2/sin(alpha2))*cos(theta2m+alpha2);
	  ycen2 = ypt2 + (rad2/sin(alpha2))*sin(theta2m+alpha2);
	  
	  /**   compute angle of point on curve from arc center **/
	  
	  th1 = atan2(x2-ycen1,x1-xcen1);
	  th2 = atan2(x2-ycen2,x1-xcen2);
	  th2t = th2 > 0.0 ? th2 : th2 + 2*M_PIE;  

	  /**  use different f depending on theta  **/
	  
	  if( (theta1-0.5*M_PIE) <= th1 && th1 <= theta1)
	    {
	      dfdx2 = cos(theta1);
	    }
	  else if ( theta2 <= th2t && (th2t - 0.5*M_PIE) <= theta2)
	    {
	      dfdx2 = cos(theta2);
	    }
	  else if ( theta2m <= (th1+0.5*M_PIE) && th1 <= (theta1-0.5*M_PIE))
	    {
	      dfdx2 = 2.*(x2-ycen1);
	    }
	  else if ( (theta2m-0.5*M_PIE) <= th2 && (th1+0.5*M_PIE) <= theta2m)
	    {
	      dfdx2 = cos(theta1m);
	    }
	  else
	    {
	      dfdx2 = 2.*(x2-ycen2);
	    }
      
	}   /*  end of microprinting (original) */
	else 
	if (model_id == 18)
	{
	  dbl xpt1, ypt1, theta1, rad1, xcen1 , ycen1, alpha1;
	  dbl xpt2, ypt2, theta2, rad2, xcen2 , ycen2, alpha2;
	  dbl theta1m, theta2m, th1, th2, th2t;
	  double omega, roll_rad, gap, width, t_offset, time1, yc, xc;
	  double y_displ, twist, angle, l_arm;
	  omega = p[1];
	  roll_rad = p[2];
	  gap = p[3];
	  width = p[4];
	  t_offset = p[5];
	  rad1=p[6];
	  rad2=p[7];
	  y_displ=p[8];
	  twist=p[9];
	  l_arm=p[10];

	  time1 = time + t_offset;
	  angle = omega*time1;
	  yc = gap + 2.*roll_rad*(1.-sin(angle));
/*  account for displacement	*/
	  yc += y_displ*sin(angle);
	  xc = y_displ*cos(angle);
/*  account for twist	*/
	  yc += l_arm*(sin(angle)-sin(angle+twist));
	  xc += l_arm*(cos(angle)-cos(angle+twist));
	  angle += twist;
	  xpt1= xc + 0.5*width*sin(angle);
	  ypt1= yc - 0.5*width*cos(angle);
	  theta1= angle;
	  xpt2= xc - 0.5*width*sin(angle);
	  ypt2= yc + 0.5*width*cos(angle);
	  theta2= angle;
	  
	/*  slope of middle line		*/

	  theta1m = atan2(ypt2-ypt1,xpt2-xpt1);
	  theta1m = theta1m > theta1 ? theta1m : theta1m + 2*M_PIE;  
	  theta2m = atan2(ypt1-ypt2,xpt1-xpt2);
	  alpha1 = 0.5*(theta1m-theta1);
	  alpha2 = 0.5*(theta2-theta2m);

	  xcen1 = xpt1 + (rad1/sin(alpha1))*cos(theta1+alpha1);
	  ycen1 = ypt1 + (rad1/sin(alpha1))*sin(theta1+alpha1);
	  xcen2 = xpt2 + (rad2/sin(alpha2))*cos(theta2m+alpha2);
	  ycen2 = ypt2 + (rad2/sin(alpha2))*sin(theta2m+alpha2);
	  
	  /**   compute angle of point on curve from arc center **/
	  
	  th1 = atan2(x2-ycen1,x1-xcen1);
	  th2 = atan2(x2-ycen2,x1-xcen2);
	  th2t = th2 > 0.0 ? th2 : th2 + 2*M_PIE;  

	  /**  use different f depending on theta  **/
	  
	  if( (theta1-0.5*M_PIE) <= th1 && th1 <= theta1)
	    {
	      dfdx2 = cos(theta1);
	    }
	  else if ( theta2 <= th2t && (th2t - 0.5*M_PIE) <= theta2)
	    {
	      dfdx2 = cos(theta2);
	    }
	  else if ( theta2m <= (th1+0.5*M_PIE) && th1 <= (theta1-0.5*M_PIE))
	    {
	      dfdx2 = 2.*(x2-ycen1);
	    }
	  else if ( (theta2m-0.5*M_PIE) <= th2 && (th1+0.5*M_PIE) <= theta2m)
	    {
	      dfdx2 = cos(theta1m);
	    }
	  else
	    {
	      dfdx2 = 2.*(x2-ycen2);
	    }
      
	}   /*  end of micro-printing (twist) */
	else 
	if (model_id == 19)
	{
	  dbl xpt, ypt, theta1, theta2, rad, xcen , ycen, alpha;
	  dbl theta, tdc_delay;
	  double omega, roll_rad, gap, width, t_offset, time1, yc;
	  omega = p[1];
	  roll_rad = p[2];
	  gap = p[3];
	  width = p[4];
	  t_offset = p[5];
	  rad=p[6];
	  theta2=p[7];
	  tdc_delay=p[8];

	  time1 = time + t_offset;
	  if (time <= 0.0)	
		{ yc = gap + 2.*roll_rad*(1.-sin(omega*time1));}
	  else if(time <= tdc_delay)
		{ yc = gap ;}
	  else 
		{ time1 -= tdc_delay;
		 yc = gap + 2.*roll_rad*(1.-sin(omega*time1));
		}
	  xpt=yc;
	  ypt=0.5*width;
	  theta1= -0.5*M_PIE;
	  
	  alpha = 0.5*(theta2-theta1);
	  xcen = xpt + (rad/sin(alpha))*cos(theta1+alpha);
	  ycen = ypt + (rad/sin(alpha))*sin(theta1+alpha);
	  
	  /**   compute angle of point on curve from arc center **/
	  
	  theta = atan2(x2-ycen,x1-xcen);
	  theta = theta > theta2-1.5*M_PIE ? theta : theta + 2*M_PIE;  
	  
	  /**  use different f depending on theta  **/
	  
	  if( (theta1-0.5*M_PIE) <= theta && theta <= (theta1+alpha))
	    {
	      dfdx2 = cos(theta1);
	    }
	  else if ( (theta1+alpha) <= theta && (theta - 0.5*M_PIE) <= theta2)
	    {
	      dfdx2 = cos(theta2);
	    }
	  else
	    {
	      dfdx2 = 2.*(x2-ycen);
	    }
      
	}   /*  end of micro-printing (axi-symmetric) */
	else 
	if (model_id == 20)
	{
	  double omega, roll_rad, t_offset, time1;
	  omega = p[1];
	  roll_rad = p[2];
	  t_offset = p[3];

	  time1 = time + t_offset;
	  dfdx2 = x2 + (x2+2.*roll_rad*sin(omega*time1));
	}
	else 
	if (model_id == 21)
	{
	  double omega, roll_rad, t_offset, time1, x_offset;
	  omega = p[1];
	  roll_rad = p[2];
	  t_offset = p[3];
	  x_offset = p[4];

	  time1 = time + t_offset;
	  dfdx2 = cos(omega*time1+x_offset/roll_rad);
	}
	else 
	if (model_id == 22)
	{
	  dbl xpt1, ypt1, theta1, rad1, xcen1 , ycen1, alpha1;
	  dbl xpt2, ypt2, theta2, rad2, xcen2 , ycen2, alpha2;
	  dbl theta1m, theta2m, th1, th2, th2t;
	  double omega_l, omega_u, roll_rad_l, roll_rad_u, gap, width;
	  double t_offset, time1, yc, xc, alpha;
	  double y_displ, twist, angle_l, angle_u, l_arm;
	  omega_l = p[1];
	  omega_u = p[2];
	  roll_rad_l = p[3];
	  roll_rad_u = p[4];
	  gap = p[5];
	  width = p[6];
	  t_offset = p[7];
	  rad1=p[8];
	  rad2=p[9];
	  alpha=p[10];
	  y_displ=p[11];
	  twist=p[12];
	  l_arm=p[13];

	  time1 = time + t_offset;
	  angle_l = omega_l*time1;
	  angle_u = omega_u*time1 + 0.5*M_PIE*(1.-omega_u/omega_l);
	  yc = gap + roll_rad_l*(1.-sin(angle_l)) + roll_rad_u*(1.-sin(angle_u));
	  xc = -roll_rad_u*cos(angle_u) + roll_rad_l*cos(angle_l);
/*  account for displacement	*/
	  yc += y_displ*sin(angle_u);
	  xc += y_displ*cos(angle_u);
/*  account for twist	*/
	  yc += l_arm*(sin(angle_u)-sin(angle_u+twist));
	  xc += l_arm*(cos(angle_u)-cos(angle_u+twist));
	  angle_u += twist;
	  xpt1= xc + 0.5*width*sin(angle_u);
	  ypt1= yc - 0.5*width*cos(angle_u);
	  theta1= angle_u - alpha;
	  xpt2= xc - 0.5*width*sin(angle_u);
	  ypt2= yc + 0.5*width*cos(angle_u);
	  theta2= angle_u + alpha;
	  
	/*  slope of middle line		*/

	  theta1m = atan2(ypt2-ypt1,xpt2-xpt1);
	  theta1m = theta1m > theta1 ? theta1m : theta1m + 2*M_PIE;  
	  theta2m = atan2(ypt1-ypt2,xpt1-xpt2);
	  alpha1 = 0.5*(theta1m-theta1);
	  alpha2 = 0.5*(theta2-theta2m);

	  xcen1 = xpt1 + (rad1/sin(alpha1))*cos(theta1+alpha1);
	  ycen1 = ypt1 + (rad1/sin(alpha1))*sin(theta1+alpha1);
	  xcen2 = xpt2 + (rad2/sin(alpha2))*cos(theta2m+alpha2);
	  ycen2 = ypt2 + (rad2/sin(alpha2))*sin(theta2m+alpha2);
	  
	  /**   compute angle of point on curve from arc center **/
	  
	  th1 = atan2(x2-ycen1,x1-xcen1);
	  th2 = atan2(x2-ycen2,x1-xcen2);
	  th2t = th2 > 0.0 ? th2 : th2 + 2*M_PIE;  

	  /**  use different f depending on theta  **/
	  
	  if( (theta1-0.5*M_PIE) <= th1 && th1 <= theta1)
	    {
	      dfdx2 = cos(theta1);
	    }
	  else if ( theta2 <= th2t && (th2t - 0.5*M_PIE) <= theta2)
	    {
	      dfdx2 = cos(theta2);
	    }
	  else if ( theta2m <= (th1+0.5*M_PIE) && th1 <= (theta1-0.5*M_PIE))
	    {
	      dfdx2 = 2.*(x2-ycen1);
	    }
	  else if ( (theta2m-0.5*M_PIE) <= th2 && (th1+0.5*M_PIE) <= theta2m)
	    {
	      dfdx2 = cos(theta1m);
	    }
	  else
	    {
	      dfdx2 = 2.*(x2-ycen2);
	    }
      
	}   /*  end of micro-printing (twist) */
        else if(model_id == 25)
        {
          dbl xx1, yy1, xx2, yy2, xx3, yy3, alpha, theta, theta1;
          dbl xcen_new,ycen_new;
          xx1=p[1];
          yy1=p[2];
          xx2=p[3];
          yy2=p[4];
          xx3=p[5];
          yy3=p[6];

          /**   compute angle of point on curve from (xcen_new,ycen_new) **/

          xcen_new = (xx1+xx3)/2;
          ycen_new = (yy1+yy3)/2.;
          theta  = atan2(x2-ycen_new,x1-xcen_new);
          theta1 = atan2(yy2-ycen_new,xx2-xcen_new);
          alpha  = atan2(yy1-yy2,xx1-xx2);

          /**  use different f depending on theta  **/

          if( theta <= theta1)
            {
              dfdx2 = 2.*(x2-ycen_new);
            }
          else
            {
              dfdx2 = cos(alpha );
            }

        }   /*  end of fluid bearing die */

        else
        if (model_id == 26)
        {
          dbl xx2, yy2, xx3, yy3, xcen , ycen, theta;
          xx2=p[3];
          yy2=p[4];
          xx3=p[5];
          yy3=p[6];

          /**   compute angle of point on curve from (xcen,ycen) **/

          xcen = (xx2+xx3)/2;
          ycen = (yy2+yy3)/2.;
          theta = atan2(x2-ycen,x1-xcen);

          /**  use different f depending on theta  **/

          if(-M_PIE/2 <= theta && theta <= 0)
            {
              dfdx2 = 1;
            }
          else if(-M_PIE <= theta && theta <= -M_PIE/2)
            {
              dfdx2 = 2*(x2-ycen);
            }
          else if(M_PIE/2 <= theta && theta <= M_PIE)
            {
              dfdx2 = 2*(x2-ycen);
            }
          else if(0 <= theta && theta <= M_PIE/2)
            {
              dfdx2 = 1;
            }

        }   /*  end of syringe with curved lip  */
	/*   multiple radiused lip corners
	 *      p[1] = number of corners
	 *      p[2] = theta0
	 *      p[3+3i] = xptl
	 *      p[4+3i] = ypt1
	 *      p[5+3i] = radius1
	 *      p[3+3n] = thetan
	 *	xcen = x-coord of radius center
	 *	ycen = y-coord of radius center
	 *	alpha = subtended half angle
	 */
	else 
	if (model_id == 27)
	{
	  dbl xpt[10],ypt[10],theta[10],rad[10],xcen[10],ycen[10],alpha[10];
	  dbl dist[10], rad_sign[10], th1, th1t;
          double theta_tmp=0, dist_min, theta1, theta2;
          int i,i_min,sign_bit,n_corner=0;
	  n_corner = ((int)p[1]);
	  sign_bit = ((int)p[2]);
          if(n_corner > 10)
          	{EH(-1,"too many corners in multiple radius lip\n");}
          i_min = 1;
          for(i=0 ; i<n_corner ; i++)
          	{
          	rad_sign[i] = 1.;
          	if(sign_bit & i_min) {  rad_sign[i] = -1.;}
          	i_min *= 2;
          	}
          theta[0] = p[3];
          for(i=0 ; i<n_corner ; i++)
          	{
          	xpt[i] = p[4+3*i];
          	ypt[i] = p[5+3*i];
          	rad[i] = p[6+3*i];
          	}
          theta[n_corner] = p[4+3*n_corner];
	  
	/*  slope of middle line		*/

          for(i=1 ; i<n_corner ; i++)
          	{
	  	theta[i] = atan2(ypt[i]-ypt[i-1],xpt[i]-xpt[i-1]);
          	}
          for(i=0 ; i<n_corner ; i++)
          	{
	  	theta[i] = theta[i] >= 0 ? theta[i] : theta[i] + 2*M_PIE;  
	  	theta[i] = theta[i] <= 2*M_PIE ? theta[i] : theta[i] - 2*M_PIE;  
          	}
          for(i=0 ; i<n_corner ; i++)
          	{
          	theta_tmp = theta[i]+M_PIE;
	  	theta_tmp = theta_tmp <= 2*M_PIE ? theta_tmp : theta_tmp-2*M_PIE;  
	  	alpha[i] = 0.5*(theta_tmp-theta[i+1]);
	  	alpha[i] = alpha[i] <= 0.5*M_PIE ? alpha[i] : M_PIE-alpha[i];  
	  	xcen[i] = xpt[i]+(rad_sign[i]*rad[i]/sin(alpha[i]))*cos(theta[i+1]+alpha[i]);
	  	ycen[i] = ypt[i]+(rad_sign[i]*rad[i]/sin(alpha[i]))*sin(theta[i+1]+alpha[i]);
          	}

	  /**   compute angle of point on curve from arc center **/
	  
          dist_min = 10000000.;   i_min = 0;
          for(i=0 ; i<n_corner ; i++)
          	{
          	dist[i] = sqrt(SQUARE(x2-ypt[i]) + SQUARE(x1-xpt[i]));
          	if(dist[i] < dist_min)
          		{dist_min = dist[i];  i_min = i;}
          	}
	  th1 = atan2(x2-ycen[i_min],x1-xcen[i_min]);
	  th1t = th1 > 0.0 ? th1 : th1 + 2*M_PIE;  

          if(rad_sign[i_min] > 0)
          	{ theta1 = theta[i_min+1];  theta2 = theta[i_min]+M_PIE;}
          else
          	{ theta1 = theta[i_min]+M_PIE;  theta2 = theta[i_min+1];}

	  /**  use different f depending on theta  **/
	  
	  if( (theta1 <= M_PIE && (theta1-0.5*M_PIE) <= th1 && th1 <= theta1)
	  || (theta1 > M_PIE && (theta1-0.5*M_PIE) <= th1t && th1t <= theta1))
	    {
	      dfdx2 = cos(theta1);
	    }
	  else 
	  if( (theta2 <= 0.5*M_PIE && (theta2+0.5*M_PIE) >= th1 && th1 >= theta2)
	  || (theta2 > 0.5*M_PIE && (theta2+0.5*M_PIE) >= th1t && th1t >= theta2))
	    {
	      dfdx2 = cos(theta2);
	    }
	  else
	    {
	      dfdx2 = 2.0*(x2-ycen[i_min]);
if(x1>0.00001)fprintf(stderr,"radius  %g %g %g %g\n",dfdx2,rad[i_min],xcen[i_min],ycen[i_min]);
	    }
      
	}   /*  end of multiple radiused lip corners */
           /*  falling circle  */
	else 
	if (model_id == 30)
	{
	  dbl ycen1, yvelo;
          dbl ycenter;
	  ycen1=p[2];
	  yvelo=p[5];
	  ycenter = ycen1 + yvelo*time;
	      dfdx2 = 2.*(x2-ycenter);
	}
	else 
	if (model_id == 31)
	{
	  dbl ycen1, deltay;
          dbl ycenter;
	  ycen1=p[2];
	  deltay=p[5];
	  ycenter = ycen1 + deltay;
	      dfdx2 = 2.*(x2-ycenter);
	}
           /*  radiused upstream bead  */
	else 
	if (model_id == 32)
	{
	  dbl radius,dca,sca, lip_angle, alpha1, alpha2;
          dbl ycenter, pos_dcl[3]={0,0,0},pos_scl[3]={0,0,0};
          int nset_dcl, nset_scl, nsp, k, j, i, dir;
	  nset_dcl = ((int)p[1]);
	  nset_scl = ((int)p[2]);
	  dca=M_PIE*p[3]/180.0;
	  sca=M_PIE*p[4]/180.0;
	  lip_angle=M_PIE*p[5]/180.0;
          alpha1 = 0.5*M_PIE - sca + lip_angle;
          alpha2 = dca - 0.5*M_PIE;
         
          /*  Get DCL and SCL coordinates  */
          nsp = match_nsid(nset_dcl);
          k = Proc_NS_List[Proc_NS_Pointers[nsp]];
          for (j = 0; j < Proc_NS_Count[nsp]; j++)
             {
               k = Proc_NS_List[Proc_NS_Pointers[nsp]+j];
               for (dir = 0; dir < pd->Num_Dim ; dir++)
                   {
                    i = Index_Solution (k, MESH_DISPLACEMENT1+dir, 0, 0, -1);
                    EH(i, "Could not resolve index_solution.");
                    /*pos_dcl[dir] = Coor[dir][k] + x[i+dir];  */
                   }
             }
          nsp = match_nsid(nset_scl);
          k = Proc_NS_List[Proc_NS_Pointers[nsp]];
          for (j = 0; j < Proc_NS_Count[nsp]; j++)
             {
               k = Proc_NS_List[Proc_NS_Pointers[nsp]+j];
               for (dir = 0; dir < pd->Num_Dim ; dir++)
                   {
                    i = Index_Solution (k, MESH_DISPLACEMENT1+dir, 0, 0, -1);
                    EH(i, "Could not resolve index_solution.");
/*                    pos_scl[dir] = Coor[dir][k] + x[i+dir];  */
                   }
             }

          /* compute upstream bead radius  */
          radius = ((pos_dcl[0]-pos_scl[0])*sin(-lip_angle)
                          +(pos_dcl[1]-pos_scl[1])*cos(-lip_angle))/
                   (sin(alpha1-lip_angle)-sin(alpha2-lip_angle));

	  ycenter = pos_dcl[1] + radius*sin(alpha2);
	  dfdx2 = 2.0*(x2-ycenter);
	}
        else
	  {
	    EH(-1,"invalid model_id in user_bc");
	  }
        return(dfdx2);
	
}
/*****************************************************************************/
dbl dfncd3(x1, x2,  x3, p, time)
     const dbl x1;
     const dbl x2;
     const dbl x3;
     const dbl p[];
     const dbl time;
{
/* for PI use M_PIE Constant from std.h include file. */
	int model_id;
	dbl dfdx3=0;

if(af->Assemble_LSA_Mass_Matrix)  return 0 ;

	model_id = ((int)p[0]);

	if (model_id == 2)
	{
        dbl xp, slope, ysurfds=0, dsdz;
        dbl cube_dep, cube_lng, flat, radius, xveloc, zangle, xoffset;
        dbl xcenter, ycenter, roll_rad, arc_pos, ysurf=0, flat_bot, xref=0;

        double x_dist[4], theta;
        cube_dep = p[1];
        cube_lng = p[2];
        flat = p[3];
        flat_bot = p[4];
        radius = p[5];
        xveloc = p[6];
        zangle = p[7];
        xoffset = p[8];
        roll_rad = p[9];
        xcenter = p[10];
        ycenter = p[11];

        theta = atan2(cube_lng*0.5,cube_dep);
        x_dist[0] = radius/tan(0.5*theta+0.25*M_PIE);
        x_dist[1] = x_dist[0]*sin(theta);
        flat = MAX(flat,2*x_dist[0]);
        flat_bot = MAX(flat_bot,2*x_dist[0]);
        slope = 2.*cube_dep/cube_lng;
        arc_pos = roll_rad * (0.5*M_PIE - acos((x1-xcenter)
                /sqrt(SQUARE(x1-xcenter)+SQUARE(x2-ycenter))));

	dsdz=-zangle;

        xp = fmod(arc_pos-xveloc*time-zangle*x3-xoffset,cube_lng+flat+flat_bot);
        if(xp < 0.0){xp += cube_lng+flat+flat_bot;}

        if(xp >= (-0.5*flat+x_dist[0]) && xp <= (0.5*flat-x_dist[0]))
                { ysurfds = 0.0; }

        else if(xp >= (0.5*flat-x_dist[0]) && xp <= (0.5*flat+x_dist[1]))
                {
                xref = 0.5*flat-x_dist[0];
                ysurfds = (xref-xp)/sqrt(SQUARE(radius)-SQUARE(xp-xref));
                }

        else if(xp >= (0.5*flat+x_dist[1]) && xp <= (0.5*flat+0.5*cube_lng-x_dist[1]))
                {
                ysurfds =  -slope;
                }

        else if(xp >= (0.5*flat+0.5*cube_lng-x_dist[1]) &&
                xp <= (0.5*flat+0.5*cube_lng+x_dist[0]))
                {
                xref = 0.5*flat + 0.5*cube_lng + x_dist[0];
                ysurfds = (xp-xref)/sqrt(SQUARE(radius)-SQUARE(xp-xref));
                }

        else if(xp >= (0.5*flat+0.5*cube_lng+x_dist[0]) &&
                xp <= (0.5*flat+0.5*cube_lng+flat_bot-x_dist[0]))
                {
                ysurfds = 0.0;
                }

        else if(xp >= (0.5*flat+0.5*cube_lng+flat_bot-x_dist[0]) &&
                xp <= (0.5*flat+0.5*cube_lng+flat_bot+x_dist[1]))
                {
                xref = 0.5*flat + 0.5*cube_lng + flat_bot - x_dist[0];
                ysurfds = (xp-xref)/sqrt(SQUARE(radius)-SQUARE(xp-xref));
                }

        else if(xp >= (0.5*flat+0.5*cube_lng+flat_bot+x_dist[1]) &&
                xp <= (0.5*flat+cube_lng+flat_bot-x_dist[1]))
                {
                ysurfds = slope;
                }
        else if(xp >= (0.5*flat+cube_lng+flat_bot-x_dist[1]) &&
                xp <= (0.5*flat+cube_lng+flat_bot+x_dist[0]))
                {
                xref = 0.5*flat + cube_lng + flat_bot + x_dist[0];
                ysurfds = (xref-xp)/sqrt(SQUARE(radius)-SQUARE(xp-xref));
                }
        else if(xp >= (0.5*flat+cube_lng+flat_bot+x_dist[0]) &&
                xp <= (flat+cube_lng+flat_bot-x_dist[0]))
                {
                ysurfds = 0.0;
                }
	else if(xp >= (flat+cube_lng+flat_bot-x_dist[0]) && 
		xp <= (flat+cube_lng+flat_bot+x_dist[1]))
        	{ 
		xref = flat+cube_lng+flat_bot-x_dist[0];
                ysurfds = (xref-xp)/sqrt(SQUARE(radius)-SQUARE(xp-xref));
		}
        else
        	{
      		EH(-1,"invalid gravure user_bc position");
        	}
	dfdx3 = -2.*(roll_rad+ysurf)*ysurfds*dsdz;

	}  /*  end of gravure if block  */
	else
	if (model_id == 4)	{

	  dbl w1, dnum, dden, dnumdx, ddendx; 
	  const dbl *pn, *pd;
	  int i, n, xyz;
	  
	  n = abs(((int)p[1]));
	  xyz = ((int)p[2]);
	  
	  
	  pn = &p[3];
	  pd = &p[n+3];
	  
	  
	  switch (xyz)
		  {
		  case 13:
		  case 23:
			  w1=x3; 
	  dnum = pn[n-1];
	  dden = pd[n-2];
	  dnumdx = (n-1.)*pn[n-1];
	  ddendx = (n-2.)*pd[n-2];
	  for (i=n-2;i>=0;i--)
	    { dnum = pn[i]+w1*dnum; }
	  for (i=n-3;i>=0;i--)
	    { dden = pd[i]+w1*dden;
	      dnumdx = (i+1.)*pn[i+1]+w1*dnumdx;}
	  for (i=n-4;i>=0;i--)
	    { ddendx = (i+1.)*pd[i+1]+w1*ddendx; }
	  
	  dfdx3 = (dden*dnumdx - dnum*ddendx)/(dden*dden);
			  break;
	  
		  case 31:
		  case 32:
			  dfdx3 = -1.; break;
		  case 21:
		  case 12:
			  dfdx3 = 0.; break;
		  default:
		EH(-1,"invalid xyz in Pade function user_bc");

		}
	   }
	    /*
		tapered surface for 3D manifolds
		driven by a Pade curve

	    p[0] = model_id
	    p[1] = theta, angle in xy-plane (degrees) 
	    p[2] = N, pade order
	    p[3 ... 2N+1] = pade coefficients

	    */
	else
	if (model_id == 6) {
	    
	  dbl x, dnum, dden, padedx, th; 
	  dbl dnumdx, ddendx;
	  const dbl *pn, *pd;
	  int i, n;
	  
	  th     = M_PIE * p[1]/180.;
	  n      = abs(((int)p[2]));
	  pn = &p[3];
	  pd = &p[n+3];
	  x = x3;
	  
	  /*
	  EVALUATE FUNCTION
	  */
	  
	  dnum = pn[n-1];
	  dden = pd[n-2];
	  dnumdx = (n-1.)*pn[n-1];
	  ddendx = (n-2.)*pd[n-2];
	  for (i=n-2;i>=0;i--)
	    { dnum = pn[i]+x*dnum; }
	  for (i=n-3;i>=0;i--)
	    { dden = pd[i]+x*dden;
	      dnumdx = (i+1.)*pn[i+1]+x*dnumdx;}
	  for (i=n-4;i>=0;i--)
	    { ddendx = (i+1.)*pd[i+1]+x*ddendx; }
	  
	  padedx = (dden*dnumdx - dnum*ddendx)/(dden*dden);
	  
	  dfdx3 = sin(th)*padedx;

	  } /* END OF tapered surface */

	else
	if (model_id == 7)	{

	  dbl theta1, theta2, rad, xcen , ycen, alpha;
	  dbl theta;
	  dbl x, dnum, dden, padex, padey;
	  dbl dnumdx, ddendx, padexdz, padeydz;
	  const dbl *pnx, *pdx, *pny, *pdy;
	  int i,nx, ny;
	  theta1=p[1];
	  theta2=p[2];
	  rad=p[3];
	  nx      = abs(((int)p[4]));
	  pnx = &p[5];
	  pdx = &p[nx+5];
	  ny      = abs(((int)p[nx+6]));
	  pny = &p[nx+7];
	  pdy = &p[nx+7+ny];
	  x = x3;

	  /*
	  EVALUATE FUNCTION
	  */
	  
	  dnum = pnx[nx-1];
	  dden = pdx[nx-2];
	  dnumdx = (nx-1.)*pnx[nx-1];
	  ddendx = (nx-2.)*pdx[nx-2];
	  for (i=nx-2;i>=0;i--)
	    { dnum = pnx[i]+x*dnum; }
	  for (i=nx-3;i>=0;i--)
	    { dden = pdx[i]+x*dden;
	      dnumdx = (i+1.)*pnx[i+1]+x*dnumdx;}
	  for (i=nx-4;i>=0;i--)
	    { ddendx = (i+1.)*pdx[i+1]+x*ddendx; }
	  
	  padex = dnum/dden;
	  padexdz = (dden*dnumdx - dnum*ddendx)/(dden*dden);

	  dnum = pny[ny-1];
	  dden = pdy[ny-2];
	  dnumdx = (ny-1.)*pny[ny-1];
	  ddendx = (ny-2.)*pdx[ny-2];
	  for (i=ny-2;i>=0;i--)
	    { dnum = pny[i]+x*dnum; }
	  for (i=ny-3;i>=0;i--)
	    { dden = pdy[i]+x*dden;
	      dnumdx = (i+1.)*pny[i+1]+x*dnumdx;}
	  for (i=ny-4;i>=0;i--)
	    { ddendx = (i+1.)*pdy[i+1]+x*ddendx; }
	  
	  padey = dnum/dden;
	  padeydz = (dden*dnumdx - dnum*ddendx)/(dden*dden);

	  alpha = 0.5*(theta2-theta1);
	  xcen = padex + (rad/sin(alpha))*cos(theta1+alpha);
	  ycen = padey + (rad/sin(alpha))*sin(theta1+alpha);

	  /**   compute angle of point on curve from arc center **/
	  
	  theta = atan2(x2-ycen,x1-xcen);
	  theta = theta > theta2-1.5*M_PIE ? theta : theta + 2*M_PIE;  
	  
	  /**  use different f depending on theta  **/
	  
	  if( (theta1-0.5*M_PIE) <= theta && theta <= (theta1+alpha))
	    {
	      dfdx3 = -padeydz*cos(theta1) + padexdz*sin(theta1);
	    }
	  else if ( (theta1+alpha) <= theta && (theta - 0.5*M_PIE) <= theta2)
	    {
	      dfdx3 = -padeydz*cos(theta2) + padexdz*sin(theta2);
	    }
	  else
	    {
	      dfdx3 = 2.*(x1-padex)*(-padexdz)+2.*(x2-padey)*(-padeydz);
	    }
      
	  }

	/*
	   Needle collar based on a tanh modified plane equation:

		Ax + By +Cz + D + E(tanh(beta*(x-x0))+1)/2 = 0
	   USES HYPERBOLIC TAN FUNCTIONS
	    */
	else
	if (model_id == 14) {
	  dfdx3 = p[3];
	  }
        else
	{
	  dfdx3 = 0.;
	}
        return(dfdx3);
}

/****************************************************************************/

void 
quser_surf (func, d_func, p, time)
     double func[DIM];
     double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE];
     double p[];  /* parameters to parameterize heat transfer model*/
     dbl time;
/******************************************************************************
*
*  Function which calculates the surface integral for user-defined heat 
*  transfer model.
*
******************************************************************************/
{
  
/* Local variables */
  
  int j_id;
  int var;
  double phi_j;
  double heat_tran_coeff;
  double epsilon, sigma, T_inf=0;
 
#if 0
  double ht_neg, ht_pos;
#endif
  double x_transition, transition_width;
  int coord_dir;
  
  double T_first, T_second, T_profile, qflow, dqflow_dT, dqflow_dX;

  /* Comment this out FIRST!!!!! */
/*  EH(-1,"No Q_USER model implemented");  */
/**************************** EXECUTION BEGINS *******************************/
  

  /* Add your function and sensitivities here */

#if 1
  heat_tran_coeff = p[0];
  T_first = p[1];
  T_second = p[2];
  x_transition = p[3];
  transition_width = p[4];
  coord_dir = ((int) p[5]);
  epsilon = p[6];
  sigma = p[7];
  dqflow_dX = 0;
  if( fv->x[coord_dir] <= (x_transition-transition_width) )
 	{
 	qflow = heat_tran_coeff*(fv->T - T_first);
 	}
 	else if (fv->x[coord_dir] >= (x_transition + transition_width) )
 	{
 	qflow = heat_tran_coeff*(fv->T - T_second);
 	}
 	else 
 	{
	T_profile = 0.5*(T_first + T_second) + 0.5*(T_first-T_second)*
 		sin(M_PIE*(fv->x[coord_dir]-x_transition)/(2*transition_width));
 	qflow = heat_tran_coeff * (fv->T - T_profile);
        dqflow_dX = 0.5*(T_first-T_second)*cos(M_PIE*(fv->x[coord_dir]-x_transition)/(2*transition_width))
		*M_PIE/(2*transition_width);
 	}
  dqflow_dT = heat_tran_coeff;
#else
  ht_neg = p[0];
  T_inf = p[1];
  epsilon = p[2];
  sigma = p[3];
  ht_pos = p[4];
  x_transition = p[5];
  transition_width = p[6];
  coord_dir = ((int) p[7]);
 
  if( fv->x[coord_dir] <= (x_transition-transition_width) )
 	{
 	heat_tran_coeff = ht_neg;
 	}
 	else if (fv->x[coord_dir] >= (x_transition + transition_width) )
 	{
 	heat_tran_coeff = ht_pos;
 	}
 	else 
 	{
 	heat_tran_coeff = 0.5*(ht_pos+ht_neg)  + 0.5*(ht_pos-ht_neg)*
 		sin(M_PIE*(fv->x[coord_dir]-x_transition)/(2*transition_width));
 	}
#endif
 
  if(af->Assemble_LSA_Mass_Matrix)
    return;
 
  if (af->Assemble_Jacobian) {

 /* sum the contributions to the global stiffness matrix */
 
    var=TEMPERATURE;
     if (pd->v[var])
        {
         for (j_id = 0; j_id < ei->dof[var]; j_id++) {
           phi_j = bf[var]->phi[j_id];
#if 1
           d_func[0][var][j_id] -= dqflow_dT * phi_j;
#else
           d_func[0][var][j_id] -= heat_tran_coeff * phi_j;
#endif
           d_func[0][var][j_id] -= 4.*epsilon*sigma * pow(fv->T,3.0) * phi_j;
           }
         }
     var = MESH_DISPLACEMENT1;
     if (pd->v[var])
        {
         for (j_id = 0; j_id < ei->dof[var]; j_id++) {
            phi_j = bf[var]->phi[j_id];
            d_func[0][var+coord_dir][j_id] -= dqflow_dX*phi_j;
           }
        }
  }
 
/* Calculate the residual contribution                                       */
   
#if 1
  func[0] = -qflow + 
            epsilon*sigma * (pow(T_inf,4.0) - pow(fv->T,4.0));
#else
  func[0] = heat_tran_coeff * (T_inf - fv->T) + 
            epsilon*sigma * (pow(T_inf,4.0) - pow(fv->T,4.0));
#endif
   
  return;
} /* END of routine quser_surf                                              */
/****************************************************************************/

/******************************************************************************
 *
 * tuser() - compute surface integral for user-defined heat transfer
 *
 * Function which calculates the surface integral for user-defined heat 
 * transfer model.
 *
 ****************************************************************************/

void 
tuser(double *func,
      double d_func[],		/* defined [MAX_VARIABLE_TYPES + MAX_CONC] */
      const double u_bc[],	/* to parameterize temperature eqn model*/
      const double time)
{
  /* 
  int var;
  double time_hr;
  */ 
  /* Comment this out FIRST!!!!! */
   EH(-1,"No T_USER  model implemented"); 

/**************************** EXECUTION BEGINS *******************************/

/* Add your function and sensitivities here */

/* Sample tuser function */
/****************************************************************************/
/*  time_hr = time/3600.;  */  /* convert time from seconds to hours */
/* */
/*  if(time_hr <= 1.)  */
/*    { */
/*      *func = fv->T - 344.;  */  /* all temperatures are in Kelvin */
/*    } */
/*  else if((time_hr > 1.) && (time_hr <= 1.2)) */
/*    { */
/*      *func = fv->T - (574 - 230*time_hr); */
/*    } */
/*  else if((time_hr > 1.2) && (time_hr <= 2.0)) */
/*    { */
/*      *func = fv->T - 298.; */
/*    } */
/*  else if((time_hr > 2.0) && (time_hr <= 11.)) */
/*    { */
/*      *func = fv->T - (7.77778*time_hr + 282.4444); */
/*    } */
/*  else if(time_hr > 11.) */
/*    { */
/*      *func = fv->T - 368.; */
/*    } */
/* */
/*  var = TEMPERATURE; */
/*  if (pd->v[var]) */
/*    { */
/*      d_func[var] = 1.; */
/*    } */

  return;
} /* END of routine tuser                                                   */
/****************************************************************************/

void 
yuser_surf(double *func,
	   double d_func[DIM][MAX_VARIABLE_TYPES+MAX_CONC][MDE],
	   const int species,
	   const double u_bc[],   /* to parameterize temperature eqn model */
	   const double time)

/******************************************************************************
 *
 * Function which calculates the surface integral for user-defined heat 
 * transfer model.
 *
 *****************************************************************************/
{
  /* 
  int var, j;
  double radius, phiw, phim, xi, alpha;
  */
/* Comment this out FIRST!!!!! */
   EH(-1,"No Y_USER model implemented"); 

/* Add your function and sensitivities here */

/* Sample yuser function */
/****************************************************************************/
/*  radius = u_bc[0]; */
/*  phim = u_bc[1]; */
/*  phiw = u_bc[2]; */
/*  xi = fv->x[1]/radius; */
/*  alpha = ( phim - phiw )/phiw; */
/* */
/*  *func = fv->c[species] - phim/( alpha*xi + 1 ); */
/* */
/*  var = MASS_FRACTION; */
/* */
/*  if (pd->v[var]) */
/*    { */
/*      for( j=0 ; j<ei->dof[var]; j++) */
/*	{ */
/*	  d_func[0][MAX_VARIABLE_TYPES + species][j] = bf[var]->phi[j]; */
/* 	}  */
/*     }  */

  return;
} /* END of routine yuser_surf                                              */
/****************************************************************************/

void 
uuser_surf (func, d_func, u_bc, time)

     double func[DIM];
     double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE];
     double u_bc[];  /* parameters to parameterize heat transfer model*/
     const dbl time;
/******************************************************************************
*
*  Function which calculates the surface integral for user-defined velocity
*
******************************************************************************/
{
/* Local variables 
  int j, j_id;
  int var;
*/
/* Comment this out FIRST!!!!! */
   EH(-1,"No U_USER model implemented"); 
  
  
/* 
 if (time <= u_bc[0])
    {
      func[0] = fv->v[0] - u_bc[2];
      var = VELOCITY1; 
      if (pd->v[var])
        { 
          for( j=0 ; j<ei->dof[var]; j++)
            { 
              d_func[0][VELOCITY1][j] = bf[var]->phi[j];
            }  
        } 
    }
  else if (time > u_bc[0] && time <= u_bc[1])
    {
      func[0] = 0.;
    }
  else
    {
      func[0] = fv->v[0] - u_bc[2];
      var = VELOCITY1; 
      if (pd->v[var])
        { 
          for( j=0 ; j<ei->dof[var]; j++)
            { 
              d_func[0][VELOCITY1][j] = bf[var]->phi[j];
            }  
        }
    }
*/

/* Add your function and sensitivities here */

  
  return;
} /* END of routine uuser_surf                                              */

/****************************************************************************/

void
uuser_colloc_surf ( double *func,
                    double d_func[],
                    const double u_bc[],  /* parameters to parameterize velocity model*/
                    const int id,         /* node ID of the collocated surface */
                    const double time )
/******************************************************************************
*
*  Function which calculates the boundary collocation for user-defined velocity
*
******************************************************************************/
{

/* Comment this out FIRST!!!!! */
   EH(-1,"No U_USER_COLLOC model implemented");

} /* END of routine uuser_colloc_surf                                       */


/****************************************************************************/

void 
vuser_surf (func, d_func, u_bc, time)

     double func[DIM];
     double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE];
     double u_bc[];  /* parameters to parameterize heat transfer model*/
     const dbl time;
/******************************************************************************
*
*  Function which calculates the surface integral for user-defined velocity
*
******************************************************************************/
{
/* Local variables */
  
/* Comment this out FIRST!!!!! */
   EH(-1,"No V_USER  model implemented"); 
  
/* Add your function and sensitivities here */


  return;
} /* END of routine vuser_surf                                              */
/****************************************************************************/

/****************************************************************************/

void
vuser_colloc_surf ( double *func,
                    double d_func[],
                    const double u_bc[],  /* parameters to parameterize velocity model*/
                    const int id,         /* node ID of the collocated surface */
                    const double time )
/******************************************************************************
*
*  Function which calculates the boundary collocation for user-defined velocity
*
******************************************************************************/
{

/* Comment this out FIRST!!!!! */
   EH(-1,"No V_USER_COLLOC model implemented");

} /* END of routine vuser_colloc_surf                                       */

/****************************************************************************/

void 
wuser_surf (func, d_func, u_bc, time)

     double func[DIM];
     double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE];
     double u_bc[];  /* parameters to parameterize heat transfer model*/
     const dbl time;
/******************************************************************************
*
*  Function which calculates the surface integral for user-defined velocity
*
******************************************************************************/
{
/* Local variables */
/**************************** EXECUTION BEGINS *******************************/
  
/* Comment this out FIRST!!!!! */
   EH(-1,"No W_USER  model implemented"); 
  
/* Add your function and sensitivities here */


  return;
} /* END of routine wuser_surf                                               */
/*****************************************************************************/

void
wuser_colloc_surf ( double *func,
                    double d_func[],
                    const double u_bc[],  /* parameters to parameterize velocity model*/
                    const int id,         /* node ID of the collocated surface */
                    const double time )
/******************************************************************************
*
*  Function which calculates the boundary collocation for user-defined velocity
*
******************************************************************************/
{

/* Comment this out FIRST!!!!! */
   EH(-1,"No W_USER_COLLOC model implemented");

} /* END of routine vuser_colloc_surf                                       */

/****************************************************************************/

void 
dx_user_surf (double *func, 
	      double d_func[],
	      const double u_bc[],
	      const double time)

/******************************************************************************
*
*  Function which calculates the surface integral for user-defined x -displacement
*
******************************************************************************/
{
  /* 
  int j_id;
  int var;
  double theta;
  */ 
  /* Comment this out FIRST!!!!! */
  EH(-1,"No DX_USER model implemented"); 

/* Nice exmample for solid body rotation */
/*
  if(fv->x0[1] >= 0.) 
    {
      theta = acos(fv->x0[0]/u_bc[0]);
    }
  else if (fv->x0[1] < 0. && fv->x0[0] < 0. ) 
    {
      theta = -asin(fv->x0[1]/u_bc[0]) + M_PIE;
    }
  else if (fv->x0[1] < 0. && fv->x0[0] >= 0.) 
    {
      theta = asin(fv->x0[1]/u_bc[0]) + 2.*M_PIE;
    }

  *func = fv->d[0]-(u_bc[0]*cos(u_bc[1]*time + theta) - fv->x0[0]);
  d_func[MESH_DISPLACEMENT1] = 1.;
 
*/

/* Add your function and sensitivities here */

  
  return;
} /* END of routine dx_user_surf                                              */
/****************************************************************************/
/****************************************************************************/

void 
dy_user_surf (double *func, 
	      double d_func[],
	      const double u_bc[],
	      const double time)
/******************************************************************************
*
*  Function which calculates the surface integral for user-defined y -displacement
*
******************************************************************************/
{
  /* 
  double theta;
  */ 
/* Comment this out FIRST!!!!! */
  EH(-1,"No DY_USER model implemented"); 

/* Add your function and sensitivities here */

/*
 if(fv->x0[1] >= 0.) theta = acos(fv->x0[0]/u_bc[0]);
  else if (fv->x0[1] < 0. && fv->x0[0] < 0. ) theta = -asin(fv->x0[1]/u_bc[0]) + M_PIE;
  else if (fv->x0[1] < 0. && fv->x0[0] >= 0.) theta = asin(fv->x0[1]/u_bc[0]) + 2.*M_PIE;

  *func = fv->d[1]-(u_bc[0]*sin(u_bc[1]*time + theta) - fv->x0[1]);
 d_func[MESH_DISPLACEMENT2] = 1.;
*/
  
  return;
} /* END of routine dy_user_surf                                              */
/****************************************************************************/
/****************************************************************************/

void 
dz_user_surf (double *func, 
	      double d_func[],
	      const double u_bc[],
	      const double time)
/******************************************************************************
*
*  Function which calculates the surface integral for user-defined z -displacement
*
******************************************************************************/
{
  /*
  */
  
/* Comment this out FIRST!!!!! */
   EH(-1,"No DZ_USER model implemented"); 

/* Add your function and sensitivities here */
  
  
  return;
} /* END of routine dz_user_surf                                              */
/****************************************************************************/

/****************************************************************************/

void 
p_liq_user_surf (double *func, 
	      double d_func[],
	      const double u_bc[],
	      const double time)
/******************************************************************************
*
*  Function which sets liquid phase pressure over a sideset
*
******************************************************************************/
{
  /*
  int j_id;
  int var;
  double phi_j;
  */
  
/* Comment this out FIRST!!!!! */
   EH(-1,"No P_LIQ_USER model implemented"); 

/*  Example.....
  if(time < u_bc[2])
    {
      *func = fv->p_liq  - u_bc[0] + (u_bc[0] - u_bc[1])*time/u_bc[2];
     d_func[POR_LIQ_PRES] = 1.;
    }
  else
    {
      *func = fv->p_liq  - u_bc[1];
      d_func[POR_LIQ_PRES] = 1.;
    }
*/

/* Add your function and sensitivities here */

  
  return;
} /* END of routine p_liq_user_surf */

/****************************************************************************/

void 
shell_p_open_user_surf (double *func, 
	      double d_func[],
	      const double u_bc[],
	      const double time)
/******************************************************************************
*
*  Function which sets liquid phase open shell pressure over a sideset
*
******************************************************************************/
{
  /*
  int j_id;
  int var;
  double phi_j;
  */
  
  static int first_time = FALSE;
  double time_off = 0.;

/* Comment this out FIRST!!!!! */
  //EH(-1,"No SHELL_P_OPEN_USER model implemented"); 
   
  if(!first_time)
    {
      if(Porous_liq_inventory <= u_bc[1])
	{
	  *func = fv->sh_p_open - u_bc[0];
	  d_func[SHELL_PRESS_OPEN] = 1.;
	}
      else
	{
	  *func = fv->sh_p_open - u_bc[0];
	  d_func[SHELL_PRESS_OPEN] = 1.;
	  first_time = 1;
	  time_off = time;
	  EH(-1,"exit here we have reached the finite load");
	}
    }
  else
    {
      if(time <= (time_off + u_bc[2]))
	{
	  *func = (fv->sh_p_open - u_bc[0])* (1. - (time-time_off)/u_bc[2]);
	  d_func[SHELL_PRESS_OPEN] = (1.- (time-time_off)/u_bc[2]);
	}
      else
	{
	  *func = 0.;
	  d_func[SHELL_PRESS_OPEN] = 0.;
	}
    }
   

/*  Example.....
  if(time < u_bc[2])
    {
      *func = fv->sh_p_open  - u_bc[0] + (u_bc[0] - u_bc[1])*time/u_bc[2];
     d_func[SHELL_PRESS_OPEN] = 1.;
    }
  else
    {
      *func = fv->sh_p_open  - u_bc[1];
      d_func[SHELL_PRESS_OPEN] = 1.;
    }
*/

/* Add your function and sensitivities here */

  
  return;
} /* END of routine shell_p_open_user_surf                                                                                        */
/****************************************************************************/
void 
mass_flux_user_surf(dbl mass_flux[MAX_CONC],
		    dbl d_mass_flux[MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC],
		    const int wspec,
		    const double p[], /* Vector of constants from input card */
		    const double time)	/* time             */
/******************************************************************************
*
*  Function which calculates the surface integral for user-defined 
*  mass flux.
*
******************************************************************************/
     
{
  /* Local variables you might want 
  int j_id;
  int var;
  double phi_j;
  double Y_w; 
*/ 
/* local concentration of current species */
/* Comment this out FIRST!!!!! */
   EH(-1,"No YFLUX_USER model implemented"); 

/* Add your function and sensitivities here */

  /* Example code snipet */
  /* for a simple linear flux */
  /*     Y_w = fv->c[wspec]; */
  /*     mass_flux[wspec] = p[0]*Y_w; */
  /*     if (af->Assemble_Jacobian )  */
  /*     { */
  /*       var=MASS_FRACTION; */
  /* 	d_mass_flux[wspec][MAX_VARIABLE_TYPES + wspec] = p[0] ; */
  /*     } */
    
  


  return;
} /* END of routine mass_flux_user_surf                                      */
/*****************************************************************************/

void
fn_dot_T_user (func, d_func, u_bc, time)
     double func[DIM];
     double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE];
     const double u_bc[];
     const dbl time;
/******************************************************************************
*
*  Function which calculates the pressure variation on a boundary 
*    for inlet conditions with variable pressure
*
******************************************************************************/
     
{
  
/* Local variables */
  
  double coat_th, gap, web_sp,visc, t_offset, pgrad;
  double press;
  int a;

/* Comment this out FIRST!!!!! 
   EH(-1,"No PRESSURE_USER model implemented"); 
*/
/* Add your function and sensitivities here */


  /*
   * Example:
   *
   *  DeformingMesh = pd->e[R_MESH1];     Catch bad references to moving 
   *				          mesh which isn't.
   *  eqn = VELOCITY1;
   *
  double pb;		
  double pa;		
  double wavelength;	
   * pb = u_bc[0];
   * pa = u_bc[1];
   * wavelength = u_bc[2];
   */
	coat_th = u_bc[0];
	gap = u_bc[1];
	web_sp = u_bc[2];
	visc = u_bc[3];
	t_offset = u_bc[4];
	pgrad = (6*visc*web_sp/SQUARE(gap))*(1.-2.*coat_th/gap);
	if(time < t_offset || fv->F < 0)
	   {press = 0.0;}
	else
	   {press = pgrad*web_sp*(time-t_offset);}
   for (a=0; a<ei->ielem_dim; a++)
   	{  
   	  func[a] -= press * fv->snormal[a]; 
   	}

	return;
  /*
   * Example:
   *
   *  DeformingMesh = pd->e[R_MESH1];     Catch bad references to moving 
   *				          mesh which isn't.
   *  eqn = VELOCITY1;
   *
   * pb = u_bc[0];
   * pa = u_bc[1];
   * wavelength = u_bc[2];
   *
   *  if (af->Assemble_Jacobian)
   *    {
   * 
   *  Evaluate sensitivity to displacements d()/dx 
   *
   *      for (jvar=0; jvar<ei->ielem_dim; jvar++)
   *	{
   *	  var = MESH_DISPLACEMENT1 + jvar;
   *	  if (pd->v[var]) 
   *	    {
   *	      for ( j=0; j<ei->dof[var]; j++)
   *		{
   *		  for (a=0; a<ei->ielem_dim; a++)
   *		    {
   *		      
   * d_press = -pa * 2 * M_PIE / wavelength 
   *               * sin(fv->x[0] * 2 * M_PIE / wavelength);
   * d_func[a][var][j] -= press * 
   *                       fv->dsnormal_dx[a][jvar][j] +
   * 	              d_press * 
   *                       fv->snormal[a] 
   *                      * bf[var]->phi[j] * delta(jvar,0);
   * pressure only a function of x 
   *
   *		    }
   *		}
   *	    }
   *	}
   *    }
   *
   * Equate boundary stress in the fluid mechanics momentum equations to 
   * external pressure
   *
   *      press = pb + pa * cos(fv->x[0] * 2 * M_PIE / wavelength);
   * pressure only a function of x 
   *
   *      for (a=0; a<ei->ielem_dim; a++)
   *	{  
   *
   *	  func[a] -= press * fv->snormal[a]; 
   *
   *	}
   */

} /* END of routine fn_dot_T_user                                            */
/*****************************************************************************/

void flow_n_dot_T_user (func, d_func, u_BC, time)
     double func[DIM];
     double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE];
     const double u_BC[];			/* Parameters from input deck */
     const dbl time;
/******************************************************************************
*
*  Function which calculates the pressure variation on a boundary 
*    for inlet conditions with variable pressure
*
******************************************************************************/
     
{
/* Local variables */

/* Comment this out FIRST!!!!! */
   EH(-1,"No FLOW_PRESSURE_USER model implemented"); 

/* Add your function and sensitivities here */


  return;		/* Here's a good default behavior! */

  /*
   * Example:
   *
   *  if (af->Assemble_Jacobian)
   *    {
   *      for (jvar=0; jvar<ei->ielem_dim; jvar++)
   *	{
   *	  var = MESH_DISPLACEMENT1 + jvar;
   *	  if (pd->v[var]) 
   *	    {
   *	      for ( j=0; j<ei->dof[var]; j++)
   *		{
   *		  for (a=0; a<ei->ielem_dim; a++)
   *		    {
   *		      
   * d_press = -u_BC[0] * 2 * M_PIE / u_BC[2] 
   *               * sin(fv->x[0] * 2 * M_PIE /u_BC[2] );
   * d_func[a][var][j] -= d_press * 
   *                       fv->snormal[a] 
   *                      * bf[var]->phi[j] * delta(jvar,0);
   *
   *		    }
   *		}
   *	    }
   *	}
   *    }
   *      press = u_BC[1] + u_BC[0] * cos(fv->x[0] * 2 * M_PIE / u_BC[2]);
   *
   *      for (a=0; a<ei->ielem_dim; a++)
   *	{  
   *	  *func -= press; 
   *	}
   */

} /* END of routine flow_n_dot_T_user                                        */
/*****************************************************************************/

double 
var_CA_user(double Ca_local, 
	    int num, 
	    const double *a, 
	    double *d_cos_CA_Ca_local)
{
  double  cos_CA;
  double  static_CA;
  double  cT;

  static_CA = a[0]*M_PIE/180.0;
  cT = a[1];

  cos_CA = cos( static_CA) - cT * Ca_local;

  *d_cos_CA_Ca_local = cT;


  /*
  double  lamda;
  double  B,C ;

  static_CA = a[0]*M_PIE/180.0;

  lamda = a[1];
  C     = a[2];

  B  = C - ( C + 1.0)* cos( static_CA) ;
  B /= C + ( C + 1.0)* cos( static_CA) ;

  cos_CA =  B*exp( lamda*Ca_local ) / ( B*exp( lamda*Ca_local) + 1.0 );

  cos_CA -= C*exp( -lamda*Ca_local) / ( B*exp( -lamda*Ca_local) + 1.0 );

  *d_cos_CA_Ca_local = B*lamda*exp( lamda*Ca_local ) / ( B*exp( lamda*Ca_local) + 1.0 ) /
                                                       ( B*exp( lamda*Ca_local) + 1.0 );
  
  *d_cos_CA_Ca_local += C*lamda*exp( -lamda*Ca_local ) / ( C*exp( -lamda*Ca_local) + 1.0 ) /
                                                         ( C*exp( -lamda*Ca_local) + 1.0 );
							 */
  return ( cos_CA );

} /* End of routine var_CA_user                                              */
/*****************************************************************************/

int
user_gibbs_criterion(const double fsnormal[MAX_PDIM], /* Vector of free surface
						       * normal components */
		     const double ssnormal[MAX_PDIM], /* Vector of solid 
						       * surface normal 
						       * components */
		     const int imodel, /* Flag which tracks which model used */
		     int *ipin  , /* Flag which tracks whether pinned or not */
		     const double p[]) /* User defined parameter list, or model
					* spec. list */

/******************************************************************************
*
*  Function which evaluates the Gibb's inequality criterion for pinning
*  and releasing a contact line, viz., 
*    (theta - theta_s)*(dist_from_salient_point) = 0 with
*        theta_s-theta >= 0, hs >= 0
*
*    If contact line should be or stay released (return 1)
*    If contact line should be fixed (return 0)
*
*            Author: P. R. Schunk    (12/16/97)
******************************************************************************/
     
{
  
/* Local variables */
  
  int a;
  dbl actual_angle, dot_prod, pos;
  dbl contact_angle;
  dbl circ_center_x, circ_center_y, circ_center_z, r_circ;
  dbl sign_orig ;
  
/***************************** EXECUTION BEGINS *******************************/

  
  if(imodel == CIRCLE)
    {
      /****Unpack things *******************/
     contact_angle = p[0];          /*  Static or dynamic contact angle             */
     circ_center_x = p[4];
     circ_center_y = p[5];
     circ_center_z = p[6];
     r_circ        = fabs(p[7]);

  /***********************************/
  /* Compute distance from sharpe edge */



  pos = fabs(r_circ) - sqrt(SQUARE(fv->x[0] - circ_center_x) +
			 SQUARE(fv->x[1] - circ_center_y) +
			 SQUARE(fv->x[2] - circ_center_z));
  
  /* 2D only for now */


if (pd->Num_Dim < 3) EH(-1,"USE CA_OR_FIX instead of CA_EDGE_OR_FIX");

/* N.B. In the distance function is encoded the original
   position of the meniscus in the sign of r_circ.  If r_circ,
   which is the 7th parameter on the card, is negative, then
   the meniscus contact line lies outside the circle.  Otherwise
   it is assumed to lie within. 
   */

sign_orig = sign_of(p[7]);

if((sign_of(pos) == sign_of(sign_orig)) && 
   (fabs(pos) > 1.e-6)) 
  {
    *ipin = 0;
    return(1);
    
  }

else 
  {
    if(!*ipin) 
       {
	 *ipin = 1;
	 return(0);
       }
    else
      {
        *ipin = 1;
      }

  /* if dist is basically zero, or the line has gone past the feature,
     the contact line should be fixed from the previous iteration or should
     be fixed now.  Evaluate Gibbs Criterion */

    dot_prod = 0.;
    for( a=0; a<pd->Num_Dim; a++)
      {
	dot_prod += fsnormal[a]*ssnormal[a];
      }
    actual_angle = 180.*acos(dot_prod)/M_PIE;
  /* evaluate gibbs criterion here */
  if (actual_angle >= (contact_angle+1.e-3))
    {
      printf("Unpinning node actual angle = %f at pos %f %f %f\n",
	     actual_angle,fv->x[0], fv->x[1], fv->x[2]);

      *ipin = 0;
      return(1);
    }
  else
    {
      *ipin = 1;
      return(0);
    }
  }
    }
  else if (imodel == B_USER)
    {
      EH(-1,"No user CA_EDGE_OR_FIX model implemented");
    }
  else
    {
      EH(-1,"Don't Recognize CA_EDGE_OR_FIX model");
    }
 return(1);
  
} /* End of routine user_gibbs_criterion                                    */ 
/****************************************************************************/

/******************************************************************************
*
*  Function which calculates the surface integral for user-defined force.
*
******************************************************************************/

void 
force_user_surf(double func[DIM],
		double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
		double p[],
		dbl time)
{

/* Local variables */



/****************** EXECUTION BEGINS added by slah jendoubi on 12/28/2006*****/
  double efield[DIM], efield_sqr, es[DIM][DIM], local_q, perm;
  int a, b;
  perm = p[1];
  func[0] = 0; func[1] = 0.; func[2]=0;

                /*
                 * computing the electric stress tensor upfront.
                 */ 
                        /*memset( es, 0, sizeof(dbl)*DIM*DIM);
                      if(pd->e[R_MESH1] || pd->e[R_MESH2])
                         {
                         }*/
                          efield_sqr = 0.0;
                          for ( a=0; a<VIM; a++)
                            {
                                  efield[a] = - fv->grad_V[a];
                                  efield_sqr += efield[a]*efield[a];
                            }

                          for ( a=0; a<VIM; a++)
                            {
                              for ( b=0; b<VIM; b++)
                                {
                                  es[a][b] = efield[a]*efield[b] - 0.5*efield_sqr*delta(a,b);
                                }
                            }

                          local_q=0;

                          for ( a=0; a<VIM; a++)
                            {
                             local_q += (-perm*es[0][a]*fv->snormal[a]);
                            }
                      func[0] = local_q;
                          local_q=0;
                          for ( a=0; a<VIM; a++)
                            {
                                   local_q += (-perm*es[1][a]*fv->snormal[a]);
                            }
                      func[1] = local_q;


/**********************Example Here for electrostatic/spring force **********/
/*  if(time < p[0])
        func[0] =  +p[1]/(p[4] - fv->x[0]) - p[2]*(fv->x[0] - p[3]);
  else if (time > p[0])
        func[0] =  -p[2]*(fv->x[0] - p[3]);
  else
        EH(-1," ran out of bounds in time baby");


  if (af->Assemble_Jacobian)
    {
      double phi_j;
      int var;
      if(time < p[0])
        {
          var = MESH_DISPLACEMENT1;
          if (pd->v[var])
            {
              for ( j=0; j<ei->dof[var]; j++)
                {
                  phi_j = bf[var]->phi[j];
                  d_func[0][var][j] +=  -p[1]*(-phi_j)/SQUARE(p[4] - fv->x[0])
                                        -p[2]*phi_j;
                }
            }
        }
      else if (time > p[0])
        {
          var = MESH_DISPLACEMENT1;
          if (pd->v[var])
            {
              for ( j=0; j<ei->dof[var]; j++)
                {
                  phi_j = bf[var]->phi[j];
                  d_func[0][var][j] +=  -p[2]*phi_j;
                }
            }

        }
      else
        EH(-1," ran out of bounds in time baby");
    }
*/

  return;
} /* END of routine force_user_surf                                          */

/*****************************************************************************/
void
volt_user_surf (double func[DIM],
                double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                double p[],
                const dbl time)

/******************************************************************************
*
*  Function which calculates the surface integral for user-defined
   voltage or potential bc model.

   Ken S. Chen (2/2006)

*
******************************************************************************/
{
/* Local variables */
  const double R = 8.314;         /* Universal gas constant in units of J/mole K */
  const double F = 96487.0;       /* Faraday's constant in units of C/equiv. */
  double RTF;                     /* R*T/F */
  int j, j_id, w1, var;
  double phi_j;
  double PHI;
  double i, i0, ai0a, Ha, cH2, cH2ref, aa, ac, T;
  double cratio;

/***************************** EXECUTION BEGINS *******************************/

/*  electrolyte potential from the linearized kinectic model for the HOR */
    i = p[0];
    ai0a = p[1];
    Ha = p[2];
    i0 = ai0a*Ha;
    cH2ref = p[3];
    aa = p[4];
    ac = p[5];
    T = p[6];

    RTF = R*T/F;
    cH2 = fv->c[0];
    if(cH2 < 0.0) cH2 = 0.0;
    if(cH2 == 0) cH2 = cH2ref;
    cratio = pow(cH2ref/cH2, 0.5);
    PHI = fv->V;

    *func = PHI + (i/i0)*cratio*(2.0/(aa+ac))*RTF;

  /* J_s_c --- sensitivity wrt species concentrations */
  var=MASS_FRACTION;
  for (j_id = 0; j_id < ei->dof[var]; j_id++)
    {
      phi_j = bf[var]->phi[j_id];

      for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++ )
        {
          d_func[0][MAX_VARIABLE_TYPES + w1][j_id] = 0.0;
        }

          d_func[0][MAX_VARIABLE_TYPES][j_id] =
                  -(i/i0)*(cratio/cH2)*RTF/(aa+ac)*phi_j;
     }

  /* J_s_T --- sensitivity wrt electrolyte solution temperature */
  var=TEMPERATURE;
  if (pd->v[var])
    {
      for (j = 0; j < ei->dof[var]; j++)
        {
          j_id = j;
          phi_j = bf[var]->phi[j_id];

          d_func[0][var][j_id] = 0.0;
        }

     }

  /* J_s_V --- sensitivity wrt electrolyte potential */
  var=VOLTAGE;
  if (pd->v[var])
    {
      for (j = 0; j < ei->dof[var]; j++)
        {
          j_id = j;
          phi_j = bf[var]->phi[j_id];

          d_func[0][var][j_id] =  phi_j;

        }
    }


  return;
} /* END of routine volt_user_surf */

/*****************************************************************************/
void 
current_user_surf (double func[DIM],
                   double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                   double p[],
                   const dbl time)

/******************************************************************************
*
*  Function which calculates the surface integral for user-defined current 
*    density model. 
*
******************************************************************************/
{
/* Local variables */
  

double wire_voltage, constant_A, Volt_s, n_power, constant_B;
double R, angle_sq;
double origin[3],dir_angle[3], costheta;
double axis_pt[3],rad_dir[3], t;

double dt_dx[3], dax_dx[3][3], dR_dx[3], drad_dx[3][3]; /* mesh displ. deriv  */

int var, j_id, w1, b, q;
double phi_j, tmp;

  
/***************************** EXECUTION BEGINS *******************************/

/***********  Electrostatic Pinning Wire ******/
wire_voltage = p[0];
constant_A = p[7];
Volt_s = p[8];
n_power = p[9];
constant_B = p[10];

/**  Use ROTATIONAL_3D snippet for pinning wire geometry	**/
/*  origin and direction of rotation axis       */
origin[0] = p[1];  origin[1] = p[2]; origin[2] = p[3];
dir_angle[0] = p[4];  dir_angle[1] = p[5]; dir_angle[2] = p[6];

/*  find intersection of axis with normal plane - i.e., locate point on
        axis that intersects plane normal to axis that contains local point. */

angle_sq = SQUARE(dir_angle[0]) + SQUARE(dir_angle[1]) + SQUARE(dir_angle[2]);
t = (dir_angle[0]*(fv->x[0]-origin[0]) + dir_angle[1]*(fv->x[1]-origin[1])
        + dir_angle[2]*(fv->x[2]-origin[2]))/angle_sq;
axis_pt[0] = origin[0]+dir_angle[0]*t;
axis_pt[1] = origin[1]+dir_angle[1]*t;
axis_pt[2] = origin[2]+dir_angle[2]*t;

/*  compute radius and radial direction */

R = sqrt( SQUARE(fv->x[0]-axis_pt[0]) + SQUARE(fv->x[1]-axis_pt[1]) +
                SQUARE(fv->x[2]-axis_pt[2]) );
rad_dir[0] = (fv->x[0]-axis_pt[0])/R;
rad_dir[1] = (fv->x[1]-axis_pt[1])/R;
rad_dir[2] = (fv->x[2]-axis_pt[2])/R;
costheta = -(fv->snormal[0]*rad_dir[0] + fv->snormal[1]*rad_dir[1] 
		+ fv->snormal[2]*rad_dir[2]);


    *func = (constant_A*(wire_voltage - fv->V - Volt_s) + constant_B)
			*pow(costheta,n_power);

  /* J_s_c --- sensitivity wrt species concentrations */
  var=MASS_FRACTION;
  if (pd->v[var])
  {
  for (j_id = 0; j_id < ei->dof[var]; j_id++)
    {
      phi_j = bf[var]->phi[j_id];
      for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++ )
        {
          d_func[0][MAX_VARIABLE_TYPES + w1][j_id] = 0.0;
        }
     }
  }

  /* J_s_V --- sensitivity wrt electrolyte potential */
  var=VOLTAGE;
  if (pd->v[var])
    {
      for (j_id = 0; j_id < ei->dof[var]; j_id++)
        {
          phi_j = bf[var]->phi[j_id];

          d_func[0][var][j_id] =  (constant_A*( -phi_j))*pow(costheta,n_power);

        }
    }

  if (pd->v[MESH_DISPLACEMENT1] )
    {
	for(b=0 ; b<3 ; b++)	{ dt_dx[b] = dir_angle[b]/angle_sq;}
	for(q=0 ; q<3 ; q++)	{
		for(b=0 ; b<3 ; b++)	
			{ dax_dx[q][b] = dir_angle[q]*dt_dx[b];}
		}
	memset( dR_dx, 0, sizeof(double)*3);
	for(q=0 ; q<3 ; q++)	{
		for(b=0 ; b<3 ; b++)	
			{ dR_dx[q] += (fv->x[b]-axis_pt[b])*(delta(b,q)-dax_dx[b][q]);}
		}
	for(q=0 ; q<3 ; q++)	{
		for(b=0 ; b<3 ; b++)	
			{ drad_dx[q][b] = (R*(delta(b,q)-dax_dx[q][b]) - 
				(fv->x[q] - axis_pt[q])*dR_dx[b])/SQUARE(R);}
		}
    	tmp = constant_A*(wire_voltage - fv->V - Volt_s) + constant_B;
       for ( b=0; b<VIM; b++)
         {
           var = MESH_DISPLACEMENT1+b;
                for (j_id=0; j_id<ei->dof[var]; j_id++)
                    {
          		phi_j = bf[var]->phi[j_id];
                         for (q=0; q<VIM; q++)
                            {
          d_func[0][var][j_id] += tmp*n_power*pow(costheta,n_power-1.)*
		(-(fv->snormal[q]*drad_dx[q][b]*phi_j + 
			fv->dsnormal_dx[q][b][j_id]*rad_dir[q]));
                             }
                    }
         }
    }

 
  return;
} /* END of routine current_user_surf */

/*****************************************************************************/
/* END of file user_bc.c */
/*****************************************************************************/
