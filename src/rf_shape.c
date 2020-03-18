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
 *$Id: rf_shape.c,v 5.4 2010-04-06 15:32:46 hkmoffa Exp $
 */

#ifdef USE_RCSID
static char rcsid[] = "$Id: rf_shape.c,v 5.4 2010-04-06 15:32:46 hkmoffa Exp $";
#endif

#include <stdlib.h>
#include <stdio.h>

#include "el_elm.h"
#include "mm_eh.h"
#include "rf_shape.h"

#define GOMA_RF_SHAPE_C


/* 
   shape() -- basic Lagrangian basis functions for finite elements

   Routine to calculate the value of the shape function psi, and
   it's partial derivatives, dpsi_s, dpsi_t and dpsi_u if applicable
   at the point (s,t,u) on the master element [-1,1]x[-1,1]x[-1,1].

   John N. Shadid  1421
   Date:		3/12/92
   Revised:        11/21/92 (SAH)
*/

double
shape (const double s,		/* quadrature point coordinates */
       const double t,
       const double u,
       const int Ielem_type,	/* element type */
       const int Iquant,	/* desired quantity (phi, phi_s, etc. */
       const int Inode)		/* current element node */
{
  double value=0;
  double temp;

  switch( Ielem_type ){	     /* select element */

  case LINEAR_TRI:         /* triangle shape functions */
  case BILINEAR_TRISHELL:

    switch( Iquant ){  /* select quantity */
    case PSI:          /* shape function */
      switch( Inode ){ /* select specific shape function */
      case 0:
	value = s;
	break;
      case 1:
	value = t;
	break;
      case 2:
	value = 1.0 - s - t;
	break;
      }
      break;

    case DPSI_S:       /* partial of shape fn w.r.t. s */
      switch( Inode ){ /* select specific shape function */
      case 0:
	value = 1.0;
	break;
      case 1:
	value = 0.0;
	break;
      case 2:
	value = -1.0;
	break;
      }
      break;

    case DPSI_T:       /* partial of shape fn w.r.t. t */
      switch( Inode ){ /* select specific shape function */
      case 0:
	value = 0.0;
	break;
      case 1:
	value = 1.0;
	break;
      case 2:
	value = -1.0;
	break;
      }
      break;

    default:
      fprintf(stderr, "Bad LINEAR_TRI case: %d!\n", Iquant);
      EH( -1, "Bad selection of phi,dphids, etc.");
      break;
    }
    break;

  case QUAD_TRI:         /* triangle shape functions */
  case QUAD6_TRI:
    temp = 1. - s - t;
    switch( Iquant ){  /* select quantity */
    case PSI:          /* shape function */
      switch( Inode ){ /* select specific shape function */
      case 0:
	value = s*(2.*s - 1.);
	break;
      case 1:
	value = t*(2.*t - 1.);
	break;
      case 2:
	value = temp*(2.*temp - 1.);
	break;
      case 3:
	value = 4.*s*t;
	break;
      case 4:
	value = 4.*t*temp;
	break;
      case 5:
	value = 4.*s*temp;
	break;
      }
      break;

    case DPSI_S:       /* partial of shape fn w.r.t. s */
      switch( Inode ){ /* select specific shape function */
      case 0:
	value = 4.*s - 1.;
	break;
      case 1:
	value = 0.0;
	break;
      case 2:
	value = 1. - 4.*temp;
	break;
      case 3:
	value = 4.*t;
	break;
      case 4:
	value = -4.*t;
	break;
      case 5:
	value = 4.*(temp-s);
	break;
      }
      break;

    case DPSI_T:       /* partial of shape fn w.r.t. t */
      switch( Inode ){ /* select specific shape function */
      case 0:
	value = 0.;
	break;
      case 1:
	value = 4.*t - 1.;
	break;
      case 2:
	value = 1. - 4.*temp;
	break;
      case 3:
	value = 4.*s;
	break;
      case 4:
	value = 4.*(temp-t);
	break;
      case 5:
	value = -4.*s;
	break;
      }
      break;

    default:
      fprintf(stderr, "Bad LINEAR_TRI case: %d!\n", Iquant);
      EH( -1, "Bad selection of phi,dphids, etc.");
      break;
    }
    break;
    
  case BILINEAR_QUAD:            /* bilinear shape functions */
  case C_BILINEAR_QUAD:
  case BILINEAR_SHELL:

    if (Inode >= 4){
      (void) printf("This element has 4 nodes!  Exiting...");
      exit(-1);
    }

    switch( Iquant ){  /* select quantity */
    case PSI:          /* shape function */
      switch( Inode ){ /* select specific shape function */
      case 0:
	value = 0.25*(1.0 - s)*(1.0 - t);
	break;
      case 1:
	value = 0.25*(1.0 + s)*(1.0 - t);
	break;
      case 2:
	value = 0.25*(1.0 + s)*(1.0 + t);
	break;
      case 3:
	value = 0.25*(1.0 - s)*(1.0 + t);
	break;
      }
      break;

    case DPSI_S:       /* partial of shape fn w.r.t. s */
      switch( Inode ){ /* select specific shape function */
      case 0:
	value = -0.25*(1.0 - t);
	break;
      case 1:
	value =  0.25*(1.0 - t);
	break;
      case 2:
	value =  0.25*(1.0 + t);
	break;
      case 3:
	value = -0.25*(1.0 + t);
	break;
      }
      break;

    case DPSI_T:       /* partial of shape fn w.r.t. t */
      switch( Inode ){ /* select specific shape function */
      case 0:
	value = -0.25*(1.0 - s);
	break;
      case 1:
	value = -0.25*(1.0 + s);
	break;
      case 2:
	value =  0.25*(1.0 + s);
	break;
      case 3:
	value =  0.25*(1.0 - s);
	break;
      }
      break;
      
    default:
      fprintf(stderr, "Bad BILINEAR_QUAD case: %d!\n", Iquant);
      fprintf(stderr, "\ts = %f\n", s);
      fprintf(stderr, "\tt = %f\n", t);
      fprintf(stderr, "\tu = %f\n", u);
      fprintf(stderr, "\tIelem_type = %d\n", Ielem_type);
      fprintf(stderr, "\tIquant     = %d\n", Iquant);
      fprintf(stderr, "\tInode      = %d\n", Inode);
      EH( -1, "Bad selection of phi,dphids, etc.");
      break;
    }
    break;

  case S_BIQUAD_QUAD: /* 8 node serendipty shape functions */
    
    if (Inode >= 8){
      (void) printf("This element has 8 nodes!  Exiting...");
      exit(-1);
    }

    switch( Iquant ){  /* select quantity */
    case PSI:          /* shape function */
      switch( Inode ){ /* select specific shape function */
      case 0:
	value = .25*(1.0 - s)*(1.0 - t)*(-s - t -1.0);
	break;
      case 1:
	value = .25*(1.0 + s)*(1.0 - t)*( s - t -1.0);
	break;
      case 2:
	value = .25*(1.0 + s)*(1.0 + t)*( s + t -1.0);
	break;
      case 3:
	value = .25*(1.0 - s)*(1.0 + t)*(-s + t -1.0);
	break;
      case 4:
	value = .5*(1.0 - s*s)*(1.0 - t);
	break;
      case 5:
	value = .5*(1.0 + s)*(1.0 - t*t);
	break;
      case 6:
	value = .5*(1.0 - s*s)*(1.0 + t);
	break;
      case 7:
	value = .5*(1.0 - s)*(1.0 - t*t);
	break;
      }
      break;

    case DPSI_S:       /* partial of shape fn w.r.t. s */
      switch( Inode ){ /* select specific shape function */
      case 0:
	value =-.25*(1.0 - t)*(-2.0*s - t);
	break;
      case 1:
	value = .25*(1.0 - t)*( 2.0*s - t);
	break;
      case 2:
	value = .25*(1.0 + t)*( 2.0*s + t);
	break;
      case 3:
	value =-.25*(1.0 + t)*(-2.0*s + t);
	break;
      case 4:
	value =-s*(1.0 - t);
	break;
      case 5:
	value = .5*(1.0 - t*t);
	break;
      case 6:
	value =-s*(1.0 + t);
	break;
      case 7:
	value =-.5*(1.0 - t*t);
	break;
      }
      break;

    case DPSI_T:       /* partial of shape fn w.r.t. t */
      switch( Inode ){ /* select specific shape function */
      case 0:
	value =-.25*(1.0 - s)*(-s - 2.0*t);
	break;
      case 1:
	value =-.25*(1.0 + s)*( s - 2.0*t); 
	break;
      case 2:
	value = .25*(1.0 + s)*( s + 2.0*t);
	break;
      case 3:
	value = .25*(1.0 - s)*(-s + 2.0*t);
	break;
      case 4:
	value =-.5*(1.0 - s*s);
	break;
      case 5:
	value =-(1.0 + s)*t;
	break;
      case 6:
	value = .5*(1.0 - s*s);
	break;
      case 7:
	value =-(1.0 - s)*t;
	break;
      }
      break;

    default:
      fprintf(stderr, "Bad S_BIQUAD_QUAD case: %d!\n", Iquant);
      EH( -1, "Bad selection of phi,dphids, etc.");
      break;
    }
    break;

  case BIQUAD_QUAD:    /* 9 node biquadratic shape fucntions */
  case BIQUAD_SHELL:
  case BIQUAD_QUAD_LS: /* biquadratic quadrilateral for level set */

    if (Inode >= 9){
      (void) fprintf(stderr, "This element has 9 nodes!  Exiting...");
      exit(-1);
    }

    switch( Iquant ){  /* select quantity */
    case PSI:          /* shape function */
      switch( Inode ){ /* select specific shape function */
      case 0:
	value = .25*(1.0-s)*(1.0 - t)*(-s - t -1.0)+.25*(1.0-s*s)*(1.0- t*t);
	break;
      case 1:
	value = .25*(1.0+s)*(1.0 - t)*( s - t -1.0)+.25*(1.0-s*s)*(1.0- t*t);
	break;
      case 2:
	value = .25*(1.0+s)*(1.0 + t)*( s + t -1.0)+.25*(1.0-s*s)*(1.0- t*t);
	break;
      case 3:
	value = .25*(1.0-s)*(1.0 + t)*(-s + t -1.0)+.25*(1.0-s*s)*(1.0- t*t);
	break;
      case 4:
	value = .5*(1.0 - s*s)*(1.0 - t) - .5*(1.0-s*s)*(1.0- t*t);
	break;
      case 5:
	value = .5*(1.0 + s)*(1.0 - t*t) - .5*(1.0-s*s)*(1.0- t*t);
	break;
      case 6:
	value = .5*(1.0 - s*s)*(1.0 + t) - .5*(1.0-s*s)*(1.0- t*t);
	break;
      case 7:
	value = .5*(1.0 - s)*(1.0 - t*t) - .5*(1.0-s*s)*(1.0- t*t);
	break;
      case 8:
	value = (1.0-s*s)*(1.0- t*t);
	break;
      }
      break;

    case DPSI_S:       /* partial of shape fn w.r.t. s */
      switch( Inode ){ /* select specific shape function */
      case 0:
	value =-.25*(1.0 - t)*(-2.0*s - t) -0.5*s*(1- t*t);
	break;
      case 1:
	value = .25*(1.0 - t)*( 2.0*s - t) -0.5*s*(1- t*t);
	break;
      case 2:
	value = .25*(1.0 + t)*( 2.0*s + t) -0.5*s*(1- t*t);
	break;
      case 3:
	value =-.25*(1.0 + t)*(-2.0*s + t) -0.5*s*(1- t*t);
	break;
      case 4:
	value =-s*(1.0 - t) + s*(1- t*t);
	break;
      case 5:
	value = .5*(1.0 - t*t) + s*(1- t*t);
	break;
      case 6:
	value =-s*(1.0 + t) + s*(1- t*t);
	break;
      case 7:
	value =-.5*(1.0 - t*t) + s*(1- t*t);
	break;
      case 8:
	value = -2.0*s*(1.0 - t*t);
	break;
      }
      break;

    case DPSI_T:       /* partial of shape fn w.r.t. t */
      switch( Inode ){ /* select specific shape function */
      case 0:
	value =-.25*(1.0 - s)*(-s - 2.0*t) -0.5*t*(1.0 - s*s);
	break;
      case 1:
	value =-.25*(1.0 + s)*( s - 2.0*t) -0.5*t*(1.0 - s*s);
	break;
      case 2:
	value = .25*(1.0 + s)*( s + 2.0*t) -0.5*t*(1.0 - s*s);
	break;
      case 3:
	value = .25*(1.0 - s)*(-s + 2.0*t) -0.5*t*(1.0 - s*s);
	break;
      case 4:
	value =-.5*(1.0 - s*s) + t*(1.0 - s*s);
	break;
      case 5:
	value =-(1.0 + s)*t  + t*(1.0 - s*s);
	break;
      case 6:
	value = .5*(1.0 - s*s)  + t*(1.0 - s*s);
	break;
      case 7:
	value =-(1.0 - s)*t  + t*(1.0 - s*s);
	break;
      case 8:
	value =-2.0*t*(1.0 - s*s);
	break;
      }
      break;

    default:
      fprintf(stderr, "Bad BIQUAD_QUAD case: %d!\n", Iquant);
      EH( -1, "Bad selection of phi,dphids, etc.");
      break;
    }
    break;

  case LINEAR_TET:
    switch (Iquant) {     /* select quantity */
    case PSI:           /* shape function */
      switch( Inode ) { /* select specific shape function */
      case 0:
	value = 1.0 - s - t - u; break;
      case 1:
	value = s;               break;
      case 2:
	value = t;               break;
      case 3:
	value = u;               break;
      }
      break;

    case DPSI_S:        /* partial of shape fn w.r.t. s */
      switch( Inode ) { /* select specific shape function */
      case 0:
	value = -1.0; break;
      case 1:
	value =  1.0; break;
      case 2:
      case 3:
	value =  0.0; break;
      }
      break;

    case DPSI_T:        /* partial of shape fn w.r.t. t */
      switch( Inode ) { /* select specific shape function */
      case 0:
	value = -1.0; break;
      case 1:
	value =  0.0; break;
      case 2:
	value =  1.0; break;
      case 3:
	value =  0.0; break;
      }
      break;

    case DPSI_U:        /* partial of shape fn w.r.t. u */
      switch( Inode ) { /* select specific shape function */
      case 0:
	value = -1.0; break;
      case 1:
      case 2:
	value =  0.0; break;
      case 3:
	value =  1.0; break;
      }
      break;

    default:
      fprintf(stderr, "ERROR: incorrect quantity %d.\n", Iquant);
      exit(-1);
      break;

    }
    
    break;  /*Case LINEAR_TET */


  
  case TRILINEAR_HEX:   /* trilinear shape functions */
  case C_TRILINEAR_HEX:

    if (Inode >= 8)
      {
	(void) printf("This element has 8 nodes!  Exiting...");
	exit(-1);
      }

    switch( Iquant ){  /* select quantity */
    case PSI:          /* shape function */
      switch( Inode ){ /* select specific shape function */
      case 0:
	value = .125*(1.0 - s)*(1.0 - t)*(1.0 - u);
	break;
      case 1:
	value = .125*(1.0 + s)*(1.0 - t)*(1.0 - u);
	break;
      case 2:
	value = .125*(1.0 + s)*(1.0 + t)*(1.0 - u);
	break;
      case 3:
	value = .125*(1.0 - s)*(1.0 + t)*(1.0 - u);
	break;
      case 4:
	value = .125*(1.0 - s)*(1.0 - t)*(1.0 + u);
	break;
      case 5:
	value = .125*(1.0 + s)*(1.0 - t)*(1.0 + u);
	break;
      case 6:
	value = .125*(1.0 + s)*(1.0 + t)*(1.0 + u);
	break;
      case 7:
	value = .125*(1.0 - s)*(1.0 + t)*(1.0 + u);
	break;
      }
      break;

    case DPSI_S:       /* partial of shape fn w.r.t. s */
      switch( Inode ){ /* select specific shape function */
      case 0:
	value = -0.125*(1.0 - t)*(1.0 - u);
	break;
      case 1:
	value =  0.125*(1.0 - t)*(1.0 - u);
	break;
      case 2:
	value =  0.125*(1.0 + t)*(1.0 - u);
	break;
      case 3:
	value = -0.125*(1.0 + t)*(1.0 - u);
	break;
      case 4:
	value = -0.125*(1.0 - t)*(1.0 + u);
	break;
      case 5:
	value =  0.125*(1.0 - t)*(1.0 + u);
	break;
      case 6:
	value =  0.125*(1.0 + t)*(1.0 + u);
	break;
      case 7:
	value = -0.125*(1.0 + t)*(1.0 + u);
	break;
      }
      break;

    case DPSI_T:       /* partial of shape fn w.r.t. t */
      switch( Inode ){ /* select specific shape function */
      case 0:
	value = -0.125*(1.0 - s)*(1.0 - u);
	break;
      case 1:
	value = -0.125*(1.0 + s)*(1.0 - u);
	break;
      case 2:
	value =  0.125*(1.0 + s)*(1.0 - u);
	break;
      case 3:
	value =  0.125*(1.0 - s)*(1.0 - u);
	break;
      case 4:
	value =  -0.125*(1.0 - s)*(1.0 + u);
	break;
      case 5:
	value =  -0.125*(1.0 + s)*(1.0 + u);
	break;
      case 6:
	value =   0.125*(1.0 + s)*(1.0 + u);
	break;
      case 7:
	value =   0.125*(1.0 - s)*(1.0 + u);
	break;
      }
      break;

    case DPSI_U:       /* partial of shape fn w.r.t. t */
      switch( Inode ){ /* select specific shape function */
      case 0:
	value = -0.125*(1.0 - s)*(1.0 - t);
	break;
      case 1:
	value = -0.125*(1.0 + s)*(1.0 - t);
	break;
      case 2:
	value = -0.125*(1.0 + s)*(1.0 + t);
	break;
      case 3:
	value = -0.125*(1.0 - s)*(1.0 + t);
	break;
      case 4:
	value =  0.125*(1.0 - s)*(1.0 - t);
	break;
      case 5:
	value =  0.125*(1.0 + s)*(1.0 - t);
	break;
      case 6:
	value =   0.125*(1.0 + s)*(1.0 + t);
	break;
      case 7:
	value =   0.125*(1.0 - s)*(1.0 + t);
	break;
      }
      break;

    default:
      fprintf(stderr, "Bad TRILINEAR_HEX case: %d!\n", Iquant);
      EH( -1, "Bad selection of phi,dphids, etc.");
      break;
    }
    break;

  case S_TRIQUAD_HEX: /* 20 node serendipty shape functions */

    if (Inode >= 20){
      (void) printf("This element has 20 nodes!  Exiting...");
      exit(-1);
    }

    switch( Iquant ){  /* select quantity */
    case PSI:          /* shape function */
      switch( Inode ){ /* select specific shape function */
      case 0:
	value = .125*(1.0 - s)*(1.0 - t)*(1.0 - u)*(-s - t - u -2.0);
	break;
      case 1:
	value = .125*(1.0 + s)*(1.0 - t)*(1.0 - u)*(s - t - u -2.0);
	break;
      case 2:
	value = .125*(1.0 + s)*(1.0 + t)*(1.0 - u)*(s + t - u -2.0);
	break;
      case 3:
	value = .125*(1.0 - s)*(1.0 + t)*(1.0 - u)*(-s + t - u -2.0);
	break;
      case 4:
	value = .125*(1.0 - s)*(1.0 - t)*(1.0 + u)*(-s - t + u -2.0);
	break;
      case 5:
	value = .125*(1.0 + s)*(1.0 - t)*(1.0 + u)*(s - t + u -2.0);
	break;
      case 6:
	value = .125*(1.0 + s)*(1.0 + t)*(1.0 + u)*(s + t + u -2.0);
	break;
      case 7:
	value = .125*(1.0 - s)*(1.0 + t)*(1.0 + u)*(-s + t + u -2.0);
	break;
      case 8:
	value = .25*(1.0 - s*s)*(1.0 - t)*(1 - u);
	break;
      case 9:
	value = .25*(1.0 + s)*(1.0 - t*t)*(1 - u);
	break;
      case 10:
	value = .25*(1.0 - s*s)*(1.0 + t)*(1 - u);
	break;
      case 11:
	value = .25*(1.0 - s)*(1.0 - t*t)*(1 - u);
	break;
      case 12:
	value = .25*(1.0 - s)*(1.0 - t)*(1 - u*u);
	break;
      case 13:
	value = .25*(1.0 + s)*(1.0 - t)*(1 - u*u);
	break;
      case 14:
	value = .25*(1.0 + s)*(1.0 + t)*(1 - u*u);
	break;
      case 15:
	value = .25*(1.0 - s)*(1.0 + t)*(1 - u*u);
	break;
      case 16:
	value = .25*(1.0 - s*s)*(1.0 - t)*(1 + u);
	break;
      case 17:
	value = .25*(1.0 + s)*(1.0 - t*t)*(1 + u);
	break;
      case 18:
	value = .25*(1.0 - s*s)*(1.0 + t)*(1 + u);
	break;
      case 19:
	value = .25*(1.0 - s)*(1.0 - t*t)*(1 + u);
	break;
      }
      break;

    case DPSI_S:       /* partial of shape fn w.r.t. s */
      switch( Inode ){ /* select specific shape function */
      case 0:
	value = .125*(1.0 - t)*(1.0 - u)*(2.0*s + t + u + 1.0);
	break;
      case 1:
	value = .125*(1.0 - t)*(1.0 - u)*(2.0*s - t - u - 1.0);
	break;
      case 2:
	value = .125*(1.0 + t)*(1.0 - u)*(2.0*s + t - u - 1.0);
	break;
      case 3:
	value = .125*(1.0 + t)*(1.0 - u)*(2.0*s - t + u + 1.0);
	break;
      case 4:
	value = .125*(1.0 - t)*(1.0 + u)*(2.0*s + t - u + 1.0);
	break;
      case 5:
	value = .125*(1.0 - t)*(1.0 + u)*(2.0*s - t + u -1.0);
	break;
      case 6:
	value = .125*(1.0 + t)*(1.0 + u)*(2.0*s + t + u -1.0);
	break;
      case 7:
	value = .125*(1.0 + t)*(1.0 + u)*(2.0*s - t - u + 1.0);
	break;
      case 8:
	value = -.5*s*(1.0 - t)*(1 - u);
	break;
      case 9:
	value =  .25*(1.0 - t*t)*(1 - u);
	break;
      case 10:
	value = -.5*s*(1.0 + t)*(1 - u);
	break;
      case 11:
	value = -.25*(1.0 - t*t)*(1 - u);
	break;
      case 12:
	value = -.25*(1.0 - t)*(1 - u*u);
	break;
      case 13:
	value =  .25*(1.0 - t)*(1 - u*u);
	break;
      case 14:
	value =  .25*(1.0 + t)*(1 - u*u);
	break;
      case 15:
	value = -.25*(1.0 + t)*(1 - u*u);
	break;
      case 16:
	value = -.5*s*(1.0 - t)*(1 + u);
	break;
      case 17:
	value =  .25*(1.0 - t*t)*(1 + u);
	break;
      case 18:
	value = -.5*s*(1.0 + t)*(1 + u);
	break;
      case 19:
	value = -.25*(1.0 - t*t)*(1 + u);
	break;
      }
      break;

    case DPSI_T:       /* partial of shape fn w.r.t. t */
      switch( Inode ){ /* select specific shape function */
      case 0:
	value = .125*(1.0 - s)*(1.0 - u)*(s + 2.0*t + u +1.0);
	break;
      case 1:
	value = .125*(1.0 + s)*(1.0 - u)*(-s + 2.0*t + u +1.0);
	break;
      case 2:
	value = .125*(1.0 + s)*(1.0 - u)*(s + 2.0*t - u -1.0);
	break;
      case 3:
	value = .125*(1.0 - s)*(1.0 - u)*(-s + 2.0*t - u -1.0);
	break;
      case 4:
	value = .125*(1.0 - s)*(1.0 + u)*(s + 2.0*t - u +1.0);
	break;
      case 5:
	value = .125*(1.0 + s)*(1.0 + u)*(-s +2.0*t - u +1.0);
	break;
      case 6:
	value = .125*(1.0 + s)*(1.0 + u)*(s + 2.0*t + u -1.0);
	break;
      case 7:
	value = .125*(1.0 - s)*(1.0 + u)*(-s + 2.0*t + u -1.0);
	break;
      case 8:
	value = -.25*(1.0 - s*s)*(1 - u);
	break;
      case 9:
	value = -.5*t*(1.0 + s)*(1 - u);
	break;
      case 10:
	value = .25*(1.0 - s*s)*(1 - u);
	break;
      case 11:
	value = -.5*t*(1.0 - s)*(1 - u);
	break;
      case 12:
	value = -.25*(1.0 - s)*(1 - u*u);
	break;
      case 13:
	value = -.25*(1.0 + s)*(1 - u*u);
	break;
      case 14:
	value = .25*(1.0 + s)*(1 - u*u);
	break;
      case 15:
	value = .25*(1.0 - s)*(1 - u*u);
	break;
      case 16:
	value = -.25*(1.0 - s*s)*(1 + u);
	break;
      case 17:
	value = -.5*t*(1.0 + s)*(1 + u);
	break;
      case 18:
	value = .25*(1.0 - s*s)*(1 + u);
	break;
      case 19:
	value = -.5*t*(1.0 - s)*(1 + u);
	break;
      }
      break;

    case DPSI_U:       /* partial of shape fn w.r.t. t */
      switch( Inode ){ /* select specific shape function */
      case 0:
	value = .125*(1.0 - s)*(1.0 - t)*(s + t + 2.0*u + 1.0);
	break;
      case 1:
	value = .125*(1.0 + s)*(1.0 - t)*(-s + t + 2.0*u + 1.0);
	break;
      case 2:
	value = .125*(1.0 + s)*(1.0 + t)*(-s - t + 2.0*u + 1.0);
	break;
      case 3:
	value = .125*(1.0 - s)*(1.0 + t)*(s - t + 2.0*u + 1.0);
	break;
      case 4:
	value = .125*(1.0 - s)*(1.0 - t)*(-s - t + 2.0*u - 1.0);
	break;
      case 5:
	value = .125*(1.0 + s)*(1.0 - t)*(s - t + 2.0*u - 1.0);
	break;
      case 6:
	value = .125*(1.0 + s)*(1.0 + t)*(s + t + 2.0*u - 1.0);
	break;
      case 7:
	value = .125*(1.0 - s)*(1.0 + t)*(-s + t + 2.0*u - 1.0);
	break;
      case 8:
	value = -.25*(1.0 - s*s)*(1.0 - t);
	break;
      case 9:
	value = -.25*(1.0 + s)*(1.0 - t*t);
	break;
      case 10:
	value = -.25*(1.0 - s*s)*(1.0 + t);
	break;
      case 11:
	value = -.25*(1.0 - s)*(1.0 - t*t);
	break;
      case 12:
	value = -.5*u*(1.0 - s)*(1.0 - t);
	break;
      case 13:
	value = -.5*u*(1.0 + s)*(1.0 - t);
	break;
      case 14:
	value = -.5*u*(1.0 + s)*(1.0 + t);
	break;
      case 15:
	value = -.5*u*(1.0 - s)*(1.0 + t);
	break;
      case 16:
	value = .25*(1.0 - s*s)*(1.0 - t);
	break;
      case 17:
	value = .25*(1.0 + s)*(1.0 - t*t);
	break;
      case 18:
	value = .25*(1.0 - s*s)*(1.0 + t);
	break;
      case 19:
	value = .25*(1.0 - s)*(1.0 - t*t);
	break;
      }
      break;

    default:
      fprintf(stderr, "Bad S_TRIQUAD_HEX case: %d!\n", Iquant);
      EH( -1, "Bad selection of phi,dphids, etc.");
      break;
    }
    break;
    
  case TRIQUAD_HEX: /* 27 node shape functions */

    if (Inode >= 27){
      (void) printf("This element has 27 nodes!  Exiting...");
      exit(-1);
    }
    
    switch( Iquant ){  /* select quantity */
    case PSI:          /* shape function */
      switch( Inode ){ /* select specific shape function */
      case 0:
	value = -0.125*s*t*u*(1.0 - s)*(1.0 - t)*(1.0 - u);
	break;
      case 1:
	value =  0.125*s*t*u*(1.0 + s)*(1.0 - t)*(1.0 - u);
	break;
      case 2:
	value = -0.125*s*t*u*(1.0 + s)*(1.0 + t)*(1.0 - u);
	break;
      case 3:
	value =  0.125*s*t*u*(1.0 - s)*(1.0 + t)*(1.0 - u);
	break;
      case  4:
	value =  0.125*s*t*u*(1.0 - s)*(1.0 - t)*(1.0 + u);
	break;
      case  5:
	value = -0.125*s*t*u*(1.0 + s)*(1.0 - t)*(1.0 + u);
	break;
      case  6:
	value =  0.125*s*t*u*(1.0 + s)*(1.0 + t)*(1.0 + u);
	break;
      case  7:
	value = -0.125*s*t*u*(1.0 - s)*(1.0 + t)*(1.0 + u);
	break;
      case 8:
	value =  0.25*t*u*(1.0 - s*s)*(1.0 - t)*(1.0 - u);
	break;
      case  9:
	value = -0.25*s*u*(1.0 + s)*(1.0 - t*t)*(1.0 - u);
	break;
      case 10:
	value = -0.25*t*u*(1.0 - s*s)*(1.0 + t)*(1.0 - u);
	break;
      case 11:
	value =  0.25*s*u*(1.0 - s)*(1.0 - t*t)*(1.0 - u);
	break;
      case 12:
	value =  0.25*s*t*(1.0 - s)*(1.0 - t)*(1.0 - u*u);
	break;
      case 13:
	value = -0.25*s*t*(1.0 + s)*(1.0 - t)*(1.0 - u*u);
	break;
      case 14:
	value =  0.25*s*t*(1.0 + s)*(1.0 + t)*(1.0 - u*u);
	break;
      case 15:
	value = -0.25*s*t*(1.0 - s)*(1.0 + t)*(1.0 - u*u);
	break;
      case 16:
	value = -0.25*t*u*(1.0 - s*s)*(1.0 - t)*(1.0 + u);
	break;
      case 17:
	value =  0.25*s*u*(1.0 + s)*(1.0 - t*t)*(1.0 + u);
	break;
      case 18:
	value =  0.25*t*u*(1.0 - s*s)*(1.0 + t)*(1.0 + u);
	break;
      case 19:
	value = -0.25*s*u*(1.0 - s)*(1.0 - t*t)*(1.0 + u);
	break;
      case 20:
	value = (1.0 - s*s)*(1.0 - t*t)*(1.0 - u*u);
	break;
      case 21:
	value = -0.5*u*(1.0 - s*s)*(1.0 - t*t)*(1.0 - u);
	break;
      case 22:
	value =  0.5*u*(1.0 - s*s)*(1.0 - t*t)*(1.0 + u);
	break;
      case 23:
	value = -0.5*s*(1.0 - s)*(1.0 - t*t)*(1.0 - u*u);
	break;
      case 24:
	value =  0.5*s*(1.0 + s)*(1.0 - t*t)*(1.0 - u*u);
	break;
      case 25:
	value = -0.5*t*(1.0 - s*s)*(1.0 - t)*(1.0 - u*u);
	break;
      case 26:
	value =  0.5*t*(1.0 - s*s)*(1.0 + t)*(1.0 - u*u);
	break;
      }
      break;
      
    case DPSI_S:       /* partial of shape fn w.r.t. s */
      switch( Inode ){ /* select specific shape function */
      case 0:
	value = -0.125*t*u*(1.0 - 2.0*s)*(1.0 - t)*(1.0 - u);
	break;
      case 1:
	value =  0.125*t*u*(1.0 + 2.0*s)*(1.0 - t)*(1.0 - u);
	break;
      case 2:
	value = -0.125*t*u*(1.0 + 2.0*s)*(1.0 + t)*(1.0 - u);
	break;
      case 3:
	value =  0.125*t*u*(1.0 - 2.0*s)*(1.0 + t)*(1.0 - u);
	break;
      case  4:
	value =  0.125*t*u*(1.0 - 2.0*s)*(1.0 - t)*(1.0 + u);
	break;
      case  5:
	value = -0.125*t*u*(1.0 + 2.0*s)*(1.0 - t)*(1.0 + u);
	break;
      case  6:
	value =  0.125*t*u*(1.0 + 2.0*s)*(1.0 + t)*(1.0 + u);
	break;
      case  7:
	value = -0.125*t*u*(1.0 - 2.0*s)*(1.0 + t)*(1.0 + u);
	break;
      case 8:
	value = -0.5*s*t*u*(1.0 - t)*(1.0 - u);
	break;
      case 9:
	value = -0.25*u*(1.0 + 2.0*s)*(1.0 - t*t)*(1.0 - u);
	break;
      case 10:
	value =  0.5*s*t*u*(1.0 + t)*(1.0 - u);
	break;
      case 11:
	value =  0.25*u*(1.0 - 2.0*s)*(1.0 - t*t)*(1.0 - u);
	break;
      case 12:
	value =  0.25*t*(1.0 - 2.0*s)*(1.0 - t)*(1.0 - u*u);
	break;
      case 13:
	value = -0.25*t*(1.0 + 2.0*s)*(1.0 - t)*(1.0 - u*u);
	break;
      case 14:
	value =  0.25*t*(1.0 + 2.0*s)*(1.0 + t)*(1.0 - u*u);
	break;
      case 15:
	value = -0.25*t*(1.0 - 2.0*s)*(1.0 + t)*(1.0 - u*u);
	break;
      case 16:
	value = 0.5*s*t*u*(1.0 - t)*(1.0 + u);
	break;
      case 17:
	value =  0.25*u*(1.0 + 2.0*s)*(1.0 - t*t)*(1.0 + u);
	break;
      case 18:
	value =  -0.5*s*t*u*(1.0 + t)*(1.0 + u);
	break;
      case 19:
	value = -0.25*u*(1.0 - 2.0*s)*(1.0 - t*t)*(1.0 + u);
	break;
      case 20:
	value = -2.0*s*(1.0 - t*t)*(1.0 - u*u);
	break;
      case 21:
	value =  s*u*(1.0 - t*t)*(1.0 - u);
	break;
      case 22:
	value =  -s*u*(1.0 - t*t)*(1.0 + u);
	break;
      case 23:
	value = -0.5*(1.0 - 2.0*s)*(1.0 - t*t)*(1.0 - u*u);
	break;
      case 24:
	value =  0.5*(1.0 + 2.0*s)*(1.0 - t*t)*(1.0 - u*u);
	break;
      case 25:
	value = s*t*(1.0 - t)*(1.0 - u*u);
	break;
      case 26:
	value =  -s*t*(1.0 + t)*(1.0 - u*u);
	break;
      }
      break;

    case DPSI_T:          /* partial of shape fn w.r.t. t */
      switch( Inode ){ /* select specific shape function */
      case 0:
	value = -0.125*s*u*(1.0 - s)*(1.0 - 2.0*t)*(1.0 - u);
	break;
      case 1:
	value =  0.125*s*u*(1.0 + s)*(1.0 - 2.0*t)*(1.0 - u);
	break;
      case 2:
	value = -0.125*s*u*(1.0 + s)*(1.0 + 2.0*t)*(1.0 - u);
	break;
      case 3:
	value =  0.125*s*u*(1.0 - s)*(1.0 + 2.0*t)*(1.0 - u);
	break;
      case 4:
	value =  0.125*s*u*(1.0 - s)*(1.0 - 2.0*t)*(1.0 + u);
	break;
      case 5:
	value = -0.125*s*u*(1.0 + s)*(1.0 - 2.0*t)*(1.0 + u);
	break;
      case  6:
	value =  0.125*s*u*(1.0 + s)*(1.0 + 2.0*t)*(1.0 + u);
	break;
      case  7:
	value = -0.125*s*u*(1.0 - s)*(1.0 + 2.0*t)*(1.0 + u);
	break;
      case 8:
	value =  0.25*u*(1.0 - s*s)*(1.0 - 2.0*t)*(1.0 - u);
	break;
      case  9:
	value = 0.5*s*t*u*(1.0 + s)*(1.0 - u);
	break;
      case 10:
	value = -0.25*u*(1.0 - s*s)*(1.0 + 2.0*t)*(1.0 - u);
	break;
      case 11:
	value =  -0.5*s*t*u*(1.0 - s)*(1.0 - u);
	break;
      case 12:
	value =  0.25*s*(1.0 - s)*(1.0 - 2.0*t)*(1.0 - u*u);
	break;
      case 13:
	value = -0.25*s*(1.0 + s)*(1.0 - 2.0*t)*(1.0 - u*u);
	break;
      case 14:
	value =  0.25*s*(1.0 + s)*(1.0 + 2.0*t)*(1.0 - u*u);
	break;
      case 15:
	value = -0.25*s*(1.0 - s)*(1.0 + 2.0*t)*(1.0 - u*u);
	break;
      case 16:
	value = -0.25*u*(1.0 - s*s)*(1.0 - 2.0*t)*(1.0 + u);
	break;
      case 17:
	value = -0.5*s*t*u*(1.0 + s)*(1.0 + u);
	break;
      case 18:
	value =  0.25*u*(1.0 - s*s)*(1.0 + 2.0*t)*(1.0 + u);
	break;
      case 19:
	value = 0.5*s*t*u*(1.0 - s)*(1.0 + u);
	break;
      case 20:
	value = -2.0*t*(1.0 - s*s)*(1.0 - u*u);
	break;
      case 21:
	value =  t*u*(1.0 - s*s)*(1.0 - u);
	break;
      case 22:
	value = -t*u*(1.0 - s*s)*(1.0 + u);
	break;
      case 23:
	value =  s*t*(1.0 - s)*(1.0 - u*u);
	break;
      case 24:
	value = -s*t*(1.0 + s)*(1.0 - u*u);
	break;
      case 25:
	value = -0.5*(1.0 - s*s)*(1.0 - 2.0*t)*(1.0 - u*u);
	break;
      case 26:
	value =  0.5*(1.0 - s*s)*(1.0 + 2.0*t)*(1.0 - u*u);
	break;
      }
      break;
      
    case DPSI_U:       /* partial of shape fn w.r.t. u */
      switch( Inode ){ /* select specific shape function */
      case 0:
	value = -0.125*s*t*(1.0 - s)*(1.0 - t)*(1.0 - 2.0*u);
	break;
      case 1:
	value =  0.125*s*t*(1.0 + s)*(1.0 - t)*(1.0 - 2.0*u);
	break;
      case 2:
	value = -0.125*s*t*(1.0 + s)*(1.0 + t)*(1.0 - 2.0*u);
	break;
      case 3:
	value =  0.125*s*t*(1.0 - s)*(1.0 + t)*(1.0 - 2.0*u);
	break;
      case 4:
	value =  0.125*s*t*(1.0 - s)*(1.0 - t)*(1.0 + 2.0*u);
	break;
      case 5:
	value = -0.125*s*t*(1.0 + s)*(1.0 - t)*(1.0 + 2.0*u);
	break;
      case  6:
	value =  0.125*s*t*(1.0 + s)*(1.0 + t)*(1.0 + 2.0*u);
	break;
      case 7:
	value = -0.125*s*t*(1.0 - s)*(1.0 + t)*(1.0 + 2.0*u);
	break;
      case 8:
	value =  0.25*t*(1.0 - s*s)*(1.0 - t)*(1.0 - 2.0*u);
	break;
      case  9:
	value = -0.25*s*(1.0 + s)*(1.0 - t*t)*(1.0 - 2.0*u);
	break;
      case 10:
	value = -0.25*t*(1.0 - s*s)*(1.0 + t)*(1.0 - 2.0*u);
	break;
      case 11:
	value =  0.25*s*(1.0 - s)*(1.0 - t*t)*(1.0 - 2.0*u);
	break;
      case 12:
	value = -0.5*s*t*u*(1.0 - s)*(1.0 - t);
	break;
      case 13:
	value =  0.5*s*t*u*(1.0 + s)*(1.0 - t);
	break;
      case 14:
	value = -0.5*s*t*u*(1.0 + s)*(1.0 + t);
	break;
      case 15:
	value =  0.5*s*t*u*(1.0 - s)*(1.0 + t);
	break;
      case 16:
	value = -0.25*t*(1.0 - s*s)*(1.0 - t)*(1.0 + 2.0*u);
	break;
      case 17:
	value =  0.25*s*(1.0 + s)*(1.0 - t*t)*(1.0 + 2.0*u);
	break;
      case 18:
	value =  0.25*t*(1.0 - s*s)*(1.0 + t)*(1.0 + 2.0*u);
	break;
      case 19:
	value = -0.25*s*(1.0 - s)*(1.0 - t*t)*(1.0 + 2.0*u);
	break;
      case 20:
	value = -2.0*u*(1.0 - s*s)*(1.0 - t*t);
	break;
      case 21:
	value = -0.5*(1.0 - s*s)*(1.0 - t*t)*(1.0 - 2.0*u);
	break;
      case 22:
	value =  0.5*(1.0 - s*s)*(1.0 - t*t)*(1.0 + 2.0*u);
	break;
      case 23:
	value =  s*u*(1.0 - s)*(1.0 - t*t);
	break;
      case 24:
	value = -s*u*(1.0 + s)*(1.0 - t*t);
	break;
      case 25:
	value =  t*u*(1.0 - s*s)*(1.0 - t);
	break;
      case 26:
	value = -t*u*(1.0 - s*s)*(1.0 + t);
	break;
      }
      break;
      
    default:
      fprintf(stderr, "Bad TRIQUAD_HEX case: %d!\n", Iquant);
      EH( -1, "Bad selection of phi,dphids, etc.");
      break;
    }

    break;

    
  case LINEAR_BAR:            /* linear shape functions */
    if (Inode >= 2){
      (void) printf("This element has 2 nodes!  Exiting...");
      exit(-1);
    }

    switch( Iquant ){  /* select quantity */
    case PSI:          /* shape function */
      switch( Inode ){ /* select specific shape function */
      case 0:
	value = 0.5 * (1.0 - s);
	break;
      case 1:
	value = 0.5 * (1.0 + s);
	break;
      }
      break;

    case DPSI_S:       /* partial of shape fn w.r.t. s */
      switch( Inode ){ /* select specific shape function */
      case 0:
	value = -0.5;
	break;
      case 1:
	value =  0.5;
	break;
      }
      break;
      
    default:
      fprintf(stderr, "Bad LINEAR_BAR case: %d!\n", Iquant);
      fprintf(stderr, "\ts = %f\n", s);
      fprintf(stderr, "\tt = %f\n", t);
      fprintf(stderr, "\tu = %f\n", u);
      fprintf(stderr, "\tIelem_type = %d\n", Ielem_type);
      fprintf(stderr, "\tIquant     = %d\n", Iquant);
      fprintf(stderr, "\tInode      = %d\n", Inode);
      EH( -1, "Bad selection of phi,dphids, etc.");
      break;
    }
    break;

  case QUAD_BAR:    /* 3 node quadratic shape fucntions */
    if (Inode >= 3){
      (void) fprintf(stderr, "This element has 9 nodes!  Exiting...");
      exit(-1);
    }

    switch( Iquant ){  /* select quantity */
    case PSI:          /* shape function */
      switch( Inode ){ /* select specific shape function */
      case 0:
	value = .50*(1.0 - s) *(-s) + (1.0- s*s);
	break;
      case 1:
	value = .50*(1.0 + s) * ( s) + (1.0-s*s);
	break;
      case 2:
	value =  (1.0-s*s);
	break;
      }
      break;

    case DPSI_S:       /* partial of shape fn w.r.t. s */
      switch( Inode ){ /* select specific shape function */
      case 0:
	value = s - 0.5 - 2.0 * s;
	break;
      case 1:
	value = s + 0.5 - 2.0 * s;
	break;
     
      case 2:
	value = -2.0*s;
	break;
      }
      break;

 
    default:
      fprintf(stderr, "Bad QUAD_BAR case: %d!\n", Iquant);
      EH( -1, "Bad selection of phi,dphids, etc.");
      break;
    }
    break;

  default:
    (void) fprintf(stderr, 
           "Element type not defined or not yet implemented.\n");
    exit(-1);
    break;
  }
  
  return value; 
  
} /* END of routine shape */
/*****************************************************************************/
/*  END of file rf_shape.c  */
/*****************************************************************************/

