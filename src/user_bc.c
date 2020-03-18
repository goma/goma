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

#include <stdio.h>
#include <math.h>

/* GOMA include files */

#include "std.h"
#include "rf_fem_const.h"
#include "rf_io.h"
#include "el_elm.h"
#include "rf_bc_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_eh.h"

/*
 * Prototype declarations of functions defined in this file.
 */

#include "user_bc.h"

/*
 * Prototype declarations of functions not defined in this file.
 */

/*
 * Function definitions.
 */

/*****************************************************************************/
/*       Functions for solid boundary description                            */
/*****************************************************************************/

dbl velo_vary_fnc(const int velo_condition,
                  const dbl x1,
                  const dbl x2,
                  const dbl x3,
                  const dbl p[],
                  const dbl time)
{
  /*  dbl v_max, gap, u, midpt, channel; */
  double f = 0.0;
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
      /* if(time<900.)
	{
	  a1 = 0.002;
	} 
      else if(time>=900.&&time<1800.)
	{
	  a1 = 0.02;
	} 
      else if(time>=1800.&&time<2700.)
	{
	  a1 = 0.2;
	}
      else if(time>=2700.&&time<3600.)
	{
	  a1 = 0.02;
	}
      else if(time>=3600.&&time<4500.)
	{
	  a1 = 0.002;
	}
      else if(time>=4500.)
	{
	  a1 = 0.;
	}
      */
    double y = x2, z = x3;

    //origin of circle
    double z0 = (2.162810-1.21031)*0.5 + 1.21031;
    double y0 = 0.0;

    double R = 0.469015; // Radius of tube

    double v_max = -4.52418;

    double coeff = v_max*(1/(R*R));

    double r = sqrt((y-y0)*(y-y0) + (z-z0)*(z-z0));

    f = coeff * (R*R - r*r);

      /*  f = -a2*x2/radius; */
    }
  else if ( velo_condition == VVARY_BC )
    {
      
      f =0;
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
dbl dvelo_vary_fnc_d1(const int velo_condition,
                      const dbl x1,
                      const dbl x2,
                      const dbl x3,
                      const dbl p[],
                      const dbl time)
{
  dbl f = 0.0;
  dbl a2;
  /* dbl radius, drdx1; */


  /* Cylindrical swirling flow 
   * Radius: a1
   * Velocity: a2
   */

  a2 = p[1];

  /* radius = sqrt( x1*x1 + x2*x2 );
  drdx1  = x1/radius; */

  if( velo_condition == VVARY_BC)
    {
     /*  f = a2*(1.0 - x1*drdx1/radius)/radius; */
      f =2.0;
    }
  else if ( velo_condition == UVARY_BC )
    {
      f = 0.;
   /*   f = -a2*x2*drdx1/radius/radius; */
    }
  else if ( velo_condition == WVARY_BC )
    {
      a2 = p[0];
      f = a2*sin(x2);
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
dbl dvelo_vary_fnc_d2(const int velo_condition,
                      const dbl x1,
                      const dbl x2,
                      const dbl x3,
                      const dbl p[],
                      const dbl time)
{
/* for PI use M_PIE Constant from std.h include file. */

  dbl f=0.0;
  dbl a2;
  /* dbl a1, radius, drdx2; */


  /* Cylindrical swirling flow 
   * Radius: a1
   * Velocity: a2
   */

  /* radius = sqrt( x1*x1 + x2*x2 );
     drdx2  = x2/radius; */
  
  a2 = p[1];
  
  if( velo_condition == VVARY_BC)
    {
/*       f = -a2*x1*drdx2/radius/radius; */
      f= 0.;
    }
  else if ( velo_condition == UVARY_BC )
    {
      /* if(time<900.)
	{
	  a1 = 0.002;
	} 
      else if(time>=900.&&time<1800.)
	{
	  a1 = 0.02;
	} 
      else if(time>=1800.&&time<2700.)
	{
	  a1 = 0.2;
	}
      else if(time>=2700.&&time<3600.)
	{
	  a1 = 0.02;
	}
      else if(time>=3600.&&time<4500.)
	{
	  a1 = 0.002;
	}
      else if(time>=4500.)
	{
	  a1 = 0.;
	}
      */
      f = -2.0;

      /*   f = ( 1. - x2*drdx2/radius)/radius; */
    }
  else if ( velo_condition == WVARY_BC ) 
    {
      a2 = p[0];
      f = a2*x1*cos(x2);
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
dbl dvelo_vary_fnc_d3(const int velo_condition,
                      const dbl x1,
                      const dbl x2,
                      const dbl x3,
                      const dbl p[],
                      const dbl time)
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
dbl fnc( const dbl x1, const dbl x2,  const dbl x3, const dbl p[], const dbl time)
{
/* for PI use M_PIE Constant from std.h include file. */

  /* cylinder with axis parallel to z */
  // dbl cx = p[0], cy=p[1], r=p[2];
  dbl f=0;


  // f = (x1-cx)*(x1-cx) + ( x2 - cy)*(x2-cy) - r*r; 
  // f = p[0] + p[1] *(x1 - p[2])/p[3] + p[4]*sin(p[5]*x1) - x2; 


  return(f);		/* Here's a good default behavior! */

}
/*****************************************************************************/
dbl dfncd1(const dbl x1,
           const dbl x2,
           const dbl x3,
           const dbl p[],
           const dbl time)
{
/* for PI use M_PIE Constant from std.h include file. */
  dbl f=0;



  // f = 2.0*(x1 - p[0] ); /* cylinder with axis pararell to z */
  //f=p[1]/p[3] + p[4]*p[5]*sin(p[5]*x1);

  return(f);		/* Here's a good default behavior! */

  /* Example code fragments:
   *
   *  dfdx1 = 1.0;
   *  return(dfdx1);
   *  f = -1;   2d fiber  
   *  f = 2.*x1;   circle 
   *  return f; 
   */

}
/*****************************************************************************/
dbl dfncd2(const dbl x1,
           const dbl x2,
           const dbl x3,
           const dbl p[],
           const dbl time)
{
/* for PI use M_PIE Constant from std.h include file. */
  dbl f=0;			/* dfdx2; */

  /* f = p[1];           time translating plane */

  //f = 2.0*(x2 - p[1]);  /* cylinder with axis pararell to z */
  
  return(f);		/* Here's a good default behavior! */

  /* Example code fragments:
   *
   *  dfdx2 = 2*x2;
   *  return(dfdx2);
   *
   *  f = 1./13 -sin(M_PIE*x2)/24. -M_PIE*x2*cos(M_PIE*x2)/24;   2d fiber 
   *  f =  2.*x2;   circle 
   *  return f; 
   */

}
/*****************************************************************************/
dbl dfncd3(const dbl x1, const dbl x2,  const dbl x3, const dbl p[], const dbl time)
{
/* for PI use M_PIE Constant from std.h include file. */
  dbl f=0;

  /*  f = p[2];        time translating plane */

  f = 0.0;  /* expanding sphere */

  return(f);		/* Here's a good default behavior! */

  /* Example code fragments:
   *
   * 
   *  f = -(1-x3*x3/200)/50.;
   *  f = 0.;
   *  return f;
   */

}

/****************************************************************************/

void 
quser_surf (double func[DIM],
            double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
            double p[],  /* parameters to parameterize heat transfer model*/
            const dbl time)
/******************************************************************************
*
*  Function which calculates the surface integral for user-defined heat 
*  transfer model.
*
******************************************************************************/
{
  /*
  int j_id;
  int var;
  double phi_j;
  double heat_xfer_coeff, d_heat_xfer_coeff_dT;
  */  

  /* Comment this out FIRST!!!!! */
  EH(-1,"No Q_USER model implemented");


  /* Add your function and sensitivities here */


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
/*  if (pd->v[pg->imtrx][var]) */
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
/*  if (pd->v[pg->imtrx][var]) */
/*    { */
/*      for( j=0 ; j<ei[pg->imtrx]->dof[var]; j++) */
/*	{ */
/*	  d_func[0][MAX_VARIABLE_TYPES + species][j] = bf[var]->phi[j]; */
/* 	}  */
/*     }  */

  return;
} /* END of routine yuser_surf                                              */
/****************************************************************************/

void 
uuser_surf (double func[DIM],
            double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
            double u_bc[],  /* parameters to parameterize heat transfer model*/
            const dbl time)
/******************************************************************************
*
*  Function which calculates the surface integral for user-defined velocity
*
******************************************************************************/
{
  /* 
  int j, j_id;
  int var;
  double phi_j;
  */ 
  
/* Comment this out FIRST!!!!! */
   EH(-1,"No U_USER model implemented"); 
  
  
/* 
 if (time <= u_bc[0])
    {
      func[0] = fv->v[0] - u_bc[2];
      var = VELOCITY1; 
      if (pd->v[pg->imtrx][var])
        { 
          for( j=0 ; j<ei[pg->imtrx]->dof[var]; j++)
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
      if (pd->v[pg->imtrx][var])
        { 
          for( j=0 ; j<ei[pg->imtrx]->dof[var]; j++)
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
vuser_surf (double func[DIM],
            double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
            double u_bc[],  /* parameters to parameterize heat transfer model*/
            const dbl time)
/******************************************************************************
*
*  Function which calculates the surface integral for user-defined velocity
*
******************************************************************************/
{
  /*
  int j_id;
  int var;
  double phi_j;
  */
  
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
wuser_surf (
     double func[DIM],
     double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
     double u_bc[],  /* parameters to parameterize heat transfer model*/
     const dbl time)
/******************************************************************************
*
*  Function which calculates the surface integral for user-defined velocity
*
******************************************************************************/
{
  /*
  int j_id;
  int var;
  double phi_j;
  */
  
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
  double phi_j;
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
  int j_id;
  int var;
  double phi_j;
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
  int j_id;
  int var;
  double phi_j;
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
    
  /* Another example:  time depending sink concentration */
  /* Comment this out FIRST!!!!! */
  /* EH(-1,"No user YFLUX_USER  model implemented"); */
  /*  p[1] = end of ramp time (beginning starts at 0)
      p[2] = starting humidity
      p[3] - ending humidity */
  /*   if(time <= p[1]) */
  /*     { */
  /*       y_inf = p[2] - (time)*(p[2]-p[3])/p[1]; */
  /*     } */
  /*   else */
  /*     { */
  /*       y_inf=p[3]; */
  /*     } */
  /*        Y_w = fv->c[wspec]; */
  /*        mass_flux[wspec]=p[0]*(Y_w - y_inf); */
  /*        if (af->Assemble_Jacobian ) */
  /*        { */
  /*          var=MASS_FRACTION; */
  /*          d_mass_flux[wspec][MAX_VARIABLE_TYPES + wspec] = p[0] ; */
  /*        } */

  return;
} /* END of routine mass_flux_user_surf                                      */
/*****************************************************************************/

void
fn_dot_T_user (double func[DIM],
               double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
               const double u_bc[],
               const dbl time)
/******************************************************************************
*
*  Function which calculates the pressure variation on a boundary 
*    for inlet conditions with variable pressure
*
******************************************************************************/
     
{
  /*  int j, j_id, i, id, var, a, eqn, I, ldof, w; */
  /*  int p, q, jvar;                   Degree of freedom counter            */
  /*  int v;                            variable counter                     */
  /*  int DeformingMesh;		Logical. */
  /*  double press, d_press;*/
  /*  double pb;			baseline of applied pressure */
  /*  double pa;			amplitude of pressure variation */
  /*  double wavelength;		wavelength of pressure  */

/* Comment this out FIRST!!!!! */
   EH(-1,"No PRESSURE_USER model implemented"); 

/* Add your function and sensitivities here */


  return;		/* Here's a good default behavior! */

  /*
   * Example:
   *
   *  DeformingMesh = pd->e[pg->imtrx][R_MESH1];     Catch bad references to moving 
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
   *      for (jvar=0; jvar<ei[pg->imtrx]->ielem_dim; jvar++)
   *	{
   *	  var = MESH_DISPLACEMENT1 + jvar;
   *	  if (pd->v[pg->imtrx][var]) 
   *	    {
   *	      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
   *		{
   *		  for (a=0; a<ei[pg->imtrx]->ielem_dim; a++)
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
   *      for (a=0; a<ei[pg->imtrx]->ielem_dim; a++)
   *	{  
   *
   *	  func[a] -= press * fv->snormal[a]; 
   *
   *	}
   */

} /* END of routine fn_dot_T_user                                            */
/*****************************************************************************/

void flow_n_dot_T_user (double func[DIM],
                        double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                        const double u_BC[],			/* Parameters from input deck */
                        const dbl time)
/******************************************************************************
*
*  Function which calculates the pressure variation on a boundary 
*    for inlet conditions with variable pressure
*
******************************************************************************/
     
{
  /*
  int j, j_id, i, id, var, a, eqn, I, ldof, w;
  int p, q, jvar;                      
  int v;
  int DeformingMesh;
  double press, d_press;
  */

/* Comment this out FIRST!!!!! */
   EH(-1,"No FLOW_PRESSURE_USER model implemented"); 

/* Add your function and sensitivities here */


  return;		/* Here's a good default behavior! */

  /*
   * Example:
   *
   *  if (af->Assemble_Jacobian)
   *    {
   *      for (jvar=0; jvar<ei[pg->imtrx]->ielem_dim; jvar++)
   *	{
   *	  var = MESH_DISPLACEMENT1 + jvar;
   *	  if (pd->v[pg->imtrx][var]) 
   *	    {
   *	      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
   *		{
   *		  for (a=0; a<ei[pg->imtrx]->ielem_dim; a++)
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
   *      for (a=0; a<ei[pg->imtrx]->ielem_dim; a++)
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
  int a;
  dbl actual_angle, dot_prod, pos;
  dbl contact_angle;
  dbl circ_center_x, circ_center_y, circ_center_z, r_circ;
  dbl sign_orig ;
  
/***************************** EXECUTION BEGINS ******************************/

  
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
  /*
  int j, j_id;
  int var;
  double phi_j;
  double heat_xfer_coeff, d_heat_xfer_coeff_dT;
  */
  
/* Comment this out FIRST!!!!! */
 EH(-1,"No FORCE_USER model implemented. Check-routine. You may be commented out!");

/**************************** EXECUTION BEGINS *******************************/
  func[0] = 0.;
/* Add your function and sensitivities here */

/*  if(pd->Num_Dim == 3)
  {
     func[1] = p[0]*sin(2.*M_PI*p[1]*time - M_PI/2.)/2. + p[0]/2.;
     func[2] = p[0]*sin(2.*M_PI*p[1]*time - M_PI/2.)/2./700. + p[0]/2./700.;
     func[0]=0;
     } 
*/
/*
  if(pd->Num_Dim == 2)
  {
     func[0] = p[0]*sin(2.*M_PI*p[1]*time - M_PI/2.)/2. + p[0]/2.;
     func[1] = func[2] = 0.;
  }
*/
/*  if (time < p[1])
    {
      func[1] = p[0];
    }
  else if (time > p[1] && time < p[2])
    {
      func[1] = 0.;
    }
  else
    {
      func[1]=p[0];
    }
*/

/**********************Example Here for electrostatic/spring force **********/
/*  if(time < p[0])
	func[0] =  +p[1]/(p[4] - fv->x[0]) - p[2]*(fv->x[0] - p[3]);
  else if (time > p[0])
	func[0] =  -p[2]*(fv->x[0] - p[3]);
  else
	EH(-1," ran out of bounds in time baby");

  if (af->Assemble_Jacobian) 
    {
      if(time < p[0])
	{
	  var = MESH_DISPLACEMENT1;
	  if (pd->v[pg->imtrx][var])
	    {
	      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
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
	  if (pd->v[pg->imtrx][var])
	    {
	      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
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
  for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++)
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
  if (pd->v[pg->imtrx][var])
    {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
        {
          j_id = j;
          phi_j = bf[var]->phi[j_id];

          d_func[0][var][j_id] = 0.0;
        }

     }

  /* J_s_V --- sensitivity wrt electrolyte potential */
  var=VOLTAGE;
  if (pd->v[pg->imtrx][var])
    {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
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
  
  
/***************************** EXECUTION BEGINS *******************************/

/***********Very Simple Example for VAR Model with J(r) on Melt Pool ******/
  /* *func=p[0]*exp(-3.0*pow(fv->x[1]/p[1],2));  */

 
  return;
} /* END of routine current_user_surf */

/*****************************************************************************/
/* END of file user_bc.c */
/*****************************************************************************/
