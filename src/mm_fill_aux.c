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
 
#include "mm_fill_aux.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "std.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_mp.h"
#include "rf_io_const.h"
#include "rf_io.h"
#include "el_elm.h"
#include "el_geom.h"
#include "rf_vars_const.h"
#include "mm_mp_const.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_mp_structs.h"
#include "mm_mp.h"
#include "exo_struct.h"
#include "dpi.h"
#include "el_elm_info.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_util.h"
#include "mm_post_def.h"
#include "mm_unknown_map.h"
#include "mm_viscosity.h"
#include "mpi.h"
#include "rf_node_const.h"
#include "mm_qtensor_model.h"

extern	dbl *p0;		/* Defined in mm_as_alloc.c */

#include "mm_eh.h"

#define GOMA_MM_FILL_AUX_C
/* 
 *  flag to use default normal and determinant in CARTESIAN
 *  coordinates - otherwise multiply by coordinate scale factors 
 */
#define  DEFAULT_CARTESIAN  0

/*
 * load_coordinate_scales -- at given local coordinates q in an orthogonal
 *			     curvilinear coordinate system, fill up the
 *			     scale factors h and the derivatives of the
 *			     scale factors with respect to each coordinate
 *			     for a couple of commonly-recognized coordinate
 *			     systems.  Also, fill up the derivatives of the 
 *			     unit vectors with respect to each coordinate.
 *
 * Reference: pp 583-584 of Bird, Armstrong, Hassager,
 *	      Dynamics of Polymeric Liquids, volume 1
 *
 *	Scale factors: h[i] = h
 *			       i
 *
 *	and a differential volume element is
 *			dV = h  h  h  dq  dq  dq 
 *			      1  2  3   1   2   3
 *
 *	Quantities loaded here...
 *
 *	h[i]		--	scale factors
 *
 *	h3		--	differential volume element is product
 *				= h h h
 *				   1 2 3
 *
 *	hq[i][j]	--	d h_i
 *				-----
 *				d q_j
 *				
 *      dh3dq[i]	--	d h3
 *				----
 *				d q_i
 *				
 *				    2
 *	hqq[i][j][k]	--	   d h_i
 *				-----------
 *				d q_j d q_k
 *
 *
 *	grad_e[a][p][q] --	( e_p e_q ) : grad(e_a)
 *
 *
 *	d_grad_e_dq[a][p][q][b]	--	d [ ( e_p e_q ) : grad(e_a) ]
 *					-----------------------------
 *						d q_b
 *
 *
 * Note that the scale factors are generally functions of position,
 * i.e.,
 *
 *	h_i = h_i(q_1, q_2, q_3).
 *
 * input:  f -- pointer to an fv structure to be loaded up...
 *	   f->x is already filled with the physical coordinates (the "q")
 *
 * output: parts of fv structure are filled in properly with scale factors
 *	   and their derivatives at this point in space.
 *
 * Revised: Tue Feb 14 12:39 MST 1995 pasacki@sandia.gov
 *
 * Revised: 9/2/99 MMH Added e_cross_e and curl_e for later computing
 *          vorticity.
 * @param c   Coordinate system
 * @param f   Field_Variable pointer.
 */

int load_coordinate_scales(const int c, struct Field_Variables *f)
{
  int a,b;
  int i,j;
  int p,q;
  int status;
  size_t siz;
#ifdef DEBUG_TIME
  int k, l;
#endif

  dbl r, theta, xi, eta, R;
  dbl ceta, seta, sxi, cxi;
  static int i_warning = 0;
  static int i_warning1 = 0;


  status = 0;

  /*
   * Initialize with default Cartesian coordinate system...
   *
   *	q = [ x, y, z ]
   */
#ifdef DEBUG_TIME
  for (i = 0; i < VIM; i++) {
    f->h[i] = 1.;
    f->dh3dq[i] = 0.;
    for (j = 0; j < VIM; j++) {
      f->hq[i][j] = 0.;
      f->curl_e[i][j] = 0.0;
      for ( k=0; k<VIM; k++) {
        f->hqq[i][j][k]    = 0.;
        f->grad_e[i][j][k] = 0.;
	for ( l=0; l<VIM; l++) {
	  f->d_grad_e_dq[i][j][k][l]  = 0.;
	}
      }
    }
  }  
#else
  
  f->h[0] = 1.;
  f->h[1] = 1.;
  f->h[2] = 1.;
  
  siz = sizeof(double)*DIM;
  memset(f->dh3dq, 0, siz);
  
  siz *= DIM;
  memset(f->hq, 0, siz);
  memset(f->curl_e, 0, siz);

  siz *= DIM;
  memset(f->hqq, 0, siz);
  memset(f->grad_e, 0, siz);
  memset(f->d_curl_e_dq, 0, siz);

#endif

  f->h3 = 1.;

  /*
   * For Cartesian coordinates we can return now. For non-Cartesian
   * coordinate systems, we simply replace the parts that need
   * replacing. Several examples are done here: adding more coordinate
   * systems is as simple as writing an order for the coordinates
   * and providing scale factors and their derivatives.
   *
   * We can also skip out early for projected Cartesian
   * coordinates... All of the scale factors here are exactly as they
   * would be in regular Cartesian. -MMH
   */

  if (c == CARTESIAN ||
      c == PROJECTED_CARTESIAN ||
      c == CARTESIAN_2pt5D)
    {
      return (status);
    }

  /*
   * Cylindrical coordinates (this ordering facilitates the solution
   *			      of axisymmetric problems where only
   *			      the first two coordinates are used.)
   *
   *	q = [ z, r, theta]
   */

  if ( c == CYLINDRICAL || c == SWIRLING )
    {
      r = f->x[1];
      if (( r == 0.) && (i_warning == 0) &&(Debug_Flag !=0))
	{
	  i_warning = 1;
	  fprintf(stderr, "Bad cylindrical coordinate -- @ r=0!\n");
	}
      f->h[2]     = r;
      f->hq[2][1] = 1.0;

      /*
       * kludge to fix up r=0 boundary. Good compilers won't
       * let you divide by zero
       */
      if (r < 1.0E-30)
	{
   	  f->h[2] = 1.0;
	  f->hq[2][1] = 0.;
	} 

    }

  /*
   * Polar coordinates (re-ordered cylindrical, the idea being that
   *			only the first two coordinates are used...)
   *
   *	q = [ r, theta, z]
   *
   */

  if ( c == POLAR )
    {
      r = f->x[0];
      if ( r == 0. && i_warning1 == 0)
	{
	  i_warning1 = 1;
	  fprintf(stderr, "Bad polar coordinate -- @ r=0!\n");
	}
      f->h[1]     = r;
      f->hq[1][0] = 1.;
    }      
    
  /*
   * Spherical coordinates...
   *
   *	q = [ r, theta, phi ]
   */

  if ( c == SPHERICAL )
    {
      r     = f->x[0];
      theta = f->x[1];
      if ( r == 0. )
	{
	  EH(-1, "Bad spherical coordinate -- @ r=0!");
	}
      if ( theta == 0. || theta == M_PIE )
	{
	  EH(-1, "Bad spherical coordinate -- @ theta=0,pi!");
	}
      f->h[1]     = r;
      f->h[2]     = r*sin(theta);

      f->hq[1][0] = 1.;
      f->hq[2][0] = sin(theta);
      f->hq[2][1] = r*cos(theta);

      f->hqq[2][0][1] = cos(theta);
      f->hqq[2][1][0] = cos(theta);

      f->hqq[2][1][1] = -r*sin(theta);
    }

  /*
   * Elliptic cylinder coordinates, (see Happel&Brenner, p.495)
   *
   * q = [ xi, eta, z ]
   */

  if (c == ELLIPTIC_CYLINDER)
    {
      xi   = f->x[0];
      eta  = f->x[1];
      sxi  = sinh(xi);
      cxi  = cosh(xi);
      seta = sin(eta);
      ceta = cos(eta);
      R    = sqrt( SQUARE(sxi) + SQUARE(seta) );

      f->h[0] = R;
      f->h[1] = R;
      
      f->hq[0][0] = sxi*cxi/R;
      f->hq[0][1] = seta*ceta/R;

      f->hq[1][0] = sxi*cxi/R;
      f->hq[1][1] = seta*ceta/R;

      f->hqq[0][0][0] = - sxi*sxi*cxi*cxi/(R*R*R)+(cxi*cxi+sxi*sxi)/R;
      f->hqq[0][0][1] = - sxi*seta*cxi*ceta/(R*R*R);
      f->hqq[0][1][0] = - sxi*seta*cxi*ceta/(R*R*R);
      f->hqq[0][1][1] = - seta*seta*ceta*ceta/(R*R*R)+(ceta*ceta-seta*seta)/R;

      f->hqq[1][0][0] = - sxi*sxi*cxi*cxi/(R*R*R)+(cxi*cxi+sxi*sxi)/R;
      f->hqq[1][0][1] = - sxi*seta*cxi*ceta/(R*R*R);
      f->hqq[1][1][0] = - sxi*seta*cxi*ceta/(R*R*R);
      f->hqq[1][1][1] = - seta*seta*ceta*ceta/(R*R*R)+(ceta*ceta-seta*seta)/R;
    }

  /*
   * For generic orthogonal curvilinear coordinate system we
   * have the differential volume element that will be needed for
   * evaluating integrals...
   */

  f->h3 = 1.;
  for (i = 0; i < VIM; i++)
    {
      f->h3 *= f->h[i];
    }

  /*
   * For mesh derivatives, we'll need any spatial dependence of the
   * volume element...
   */
  
  for (i = 0; i < VIM; i++)
    {
      f->dh3dq[i] = 0.;
      for (j = 0; j < VIM; j++)
	{
	  f->dh3dq[i] += (f->hq[j][i])/(f->h[j]);
	}
      f->dh3dq[i] *= f->h3;
    }

  /*
   * If all the h, dh/dq are known, then the gradient of the unit
   * vectors can be computed...
   */

#ifndef DO_NO_UNROLL
  for (a = 0; a < VIM; a++)
    {
      for (p = 0; p < VIM; p++)
	{
	  f->grad_e[a][p][p]  = f->hq[p][a] / ( f->h[p] * f->h[a] );
	}
    }

  for (a = 0; a < VIM; a++)
    {
      for (q = 0; q < VIM; q++)
	{
	  f->grad_e[a][a][q] -= f->hq[a][q] / ( f->h[a] * f->h[q] );
	}
    }
#else
  f->grad_e[0][0][0]  = f->hq[0][0] / ( f->h[0] * f->h[0] );
  f->grad_e[1][1][1]  = f->hq[1][1] / ( f->h[1] * f->h[1] );
  f->grad_e[0][1][1]  = f->hq[1][0] / ( f->h[1] * f->h[0] );
  f->grad_e[1][0][0]  = f->hq[0][1] / ( f->h[0] * f->h[1] );


  f->grad_e[0][0][0] -= f->hq[0][0] / ( f->h[0] * f->h[0] );
  f->grad_e[1][1][1] -= f->hq[1][1] / ( f->h[1] * f->h[1] );
  f->grad_e[0][0][1] -= f->hq[0][1] / ( f->h[0] * f->h[1] );
  f->grad_e[1][1][0] -= f->hq[1][0] / ( f->h[1] * f->h[0] );
  
  if (VIM == 3) 
    {
      f->grad_e[2][2][2]  = f->hq[2][2] / ( f->h[2] * f->h[2] );
      f->grad_e[2][0][0]  = f->hq[0][2] / ( f->h[0] * f->h[2] );
      f->grad_e[2][1][1]  = f->hq[1][2] / ( f->h[1] * f->h[2] );
      f->grad_e[1][2][2]  = f->hq[2][1] / ( f->h[2] * f->h[1] );
      f->grad_e[0][2][2]  = f->hq[2][0] / ( f->h[2] * f->h[0] );
	  
      f->grad_e[2][2][2] -= f->hq[2][2] / ( f->h[2] * f->h[2] );
      f->grad_e[2][2][0] -= f->hq[2][0] / ( f->h[2] * f->h[0] );
      f->grad_e[2][2][1] -= f->hq[2][1] / ( f->h[2] * f->h[1] );
      f->grad_e[1][1][2] -= f->hq[1][2] / ( f->h[1] * f->h[2] );
      f->grad_e[0][0][2] -= f->hq[0][2] / ( f->h[0] * f->h[2] );
    }
#endif
		 
  /* MMH: Compute the curl of the basis vectors.
   */ 
  /* MMH: Always compute all three components of the vorticity vector.
   * If we are really in 2D, then only the third component in the
   * output is of interest.
   */

  if (CURL_V != -1)
    {

      /* Note d_curl_e_dq isn't computed in the DO_NOT_UNROLL block
       */
#ifdef DO_NOT_UNROLL
      for (a = 0; a < DIM; a++)
	for (p = 0; p < DIM; p++) /* VIM */
	  for (q = 0; q < DIM; q++) /* VIM */
	    f->curl_e[a][q] -= permute(a,p,q) / f->h[a] / f->h[p] * f->hq[a][p];
#else
      if (c == CYLINDRICAL || c == SWIRLING)
	{      
	  f->curl_e[2][0] = 1.0 / f->h[2];

	  f->d_curl_e_dq[2][0][1] = - f->hq[2][1] / (f->h[2] * f->h[2]);
	}
      else
	{
	  /* f->curl_e[a][q] -= (permute(a,0,q) / f->h[a] / f->h[0] * f->hq[a][0] +
	                         permute(a,1,q) / f->h[a] / f->h[1] * f->hq[a][1] +
	                         permute(a,2,q) / f->h[a] / f->h[2] * f->hq[a][2];   */

	  f->curl_e[0][1] = -(permute(0,2,1) / f->h[0] / f->h[2] * f->hq[0][2]);
					   
	  f->curl_e[1][0] = -(permute(1,2,0) / f->h[1] / f->h[2] * f->hq[1][2]);

	  f->curl_e[0][2] = -(permute(0,1,2) / f->h[0] / f->h[1] * f->hq[0][1]);

	  f->curl_e[1][2] = -(permute(1,0,2) / f->h[1] / f->h[0] * f->hq[1][0]);

	  f->curl_e[2][0] = -(permute(2,1,0) / f->h[2] / f->h[1] * f->hq[2][1]);

	  f->curl_e[2][1] = -(permute(2,0,1) / f->h[2] / f->h[0] * f->hq[2][0]);

	  for (b = 0; b < DIM; b++) 
	    {
	      f->d_curl_e_dq[0][1][b] = ( - f->curl_e[0][1] / f->h[0] * f->hq[0][b]
					  - f->curl_e[0][1] / f->h[2] * f->hq[2][b]
					  - permute(0,2,1)  / f->h[0] / f->h[2] * f->hqq[0][2][b]);


	      f->d_curl_e_dq[1][0][b] = ( - f->curl_e[1][0] / f->h[1] * f->hq[1][b]
					  - f->curl_e[1][0] / f->h[2] * f->hq[2][b]
					  - permute(1,2,0)  / f->h[1] / f->h[2] * f->hqq[1][2][b]);

	      f->d_curl_e_dq[0][2][b] = ( - f->curl_e[0][2] / f->h[0] * f->hq[0][b]
					  - f->curl_e[0][2] / f->h[1] * f->hq[1][b]
					  - permute(0,1,2)  / f->h[0] / f->h[1] * f->hqq[0][1][b]);

	      f->d_curl_e_dq[1][2][b] = ( - f->curl_e[1][2] / f->h[1] * f->hq[1][b]
					  - f->curl_e[1][2] / f->h[0] * f->hq[0][b]
					  - permute(1,0,2)  / f->h[1] / f->h[0] * f->hqq[1][0][b]);

	      f->d_curl_e_dq[2][0][b] = ( - f->curl_e[2][0] / f->h[2] * f->hq[2][b]
					  - f->curl_e[2][0] / f->h[1] * f->hq[1][b]
					  - permute(2,1,0)  / f->h[2] / f->h[1] * f->hqq[2][1][b]);

	      f->d_curl_e_dq[2][1][b] = ( - f->curl_e[2][1] / f->h[2] * f->hq[2][b]
					  - f->curl_e[2][1] / f->h[0] * f->hq[0][b]
					  - permute(2,0,1)  / f->h[2] / f->h[0] * f->hqq[2][0][b]);
	    }
	}
#endif
    }

  /*
   * Finally, we need to know the physical space gradient of the gradient
   * of the unit vectors...
   */
  
  /*   for ( a=0; a<VIM; a++)
   *     {
   *       for ( p=0; p<VIM; p++)
   *         {
   *           for (q=0; q<VIM; q++)
   *             {
   *               for ( b=0; b<VIM; b++)
   *                 {
   *                    f->d_grad_e_dq[a][p][q][b] =
   *                               ( delta(p,q) / ( SQUARE(f->h[p]) * SQUARE(f->h[a]) ) )
   *                                 *
   *                                  ( - (f->h[a])*(f->hq[p][b])*(f->hq[p][a])
   *                                    - (f->h[p])*(f->hq[a][b])*(f->hq[p][a])
   *                                    + (f->h[p])*(f->h[a])*(f->hqq[p][a][b]) )
   *                               -
   *                               ( delta(a,p) / ( SQUARE(f->h[p]) * SQUARE(f->h[q]) ) )
   *                                 *
   *                                 ( - (f->h[q])*(f->hq[p][b])*(f->hq[a][q])
   *                                   - (f->h[p])*(f->hq[q][b])*(f->hq[a][q])
   *                                   + (f->h[p])*(f->h[q])*(f->hqq[a][q][b]) );
   *                }
   *            }
   *        }
   *    }   
   */

	
  if( pd->v[pg->imtrx][POLYMER_STRESS11] || ei[pg->imtrx]->deforming_mesh )   /* these are the only two cases where d_grad_e_dq is used */
    {
      siz = sizeof(double) * DIM * DIM * DIM * DIM;
      memset(f->d_grad_e_dq, 0, siz);

      for ( a=0; a<VIM; a++)
	{
#ifndef DO_NO_UNROLL
	  for ( p=0; p<VIM; p++)
	    {
	      for ( b=0; b<VIM; b++)
		{
		  f->d_grad_e_dq[a][p][p][b] = ( 1. / ( (f->h[p])*(f->h[p]) *(f->h[a])*(f->h[a]) ) )
		    * ( - (f->h[a])*(f->hq[p][b])*(f->hq[p][a])
			- (f->h[p])*(f->hq[a][b])*(f->hq[p][a])
			+ (f->h[p])*(f->h[a])*(f->hqq[p][a][b]) );
		}
	    }
#else

	  f->d_grad_e_dq[a][0][0][0] = ( 1. / ( (f->h[0])*(f->h[0]) *(f->h[a])*(f->h[a]) ) )
	    * ( - (f->h[a])*(f->hq[0][0])*(f->hq[0][a])
		- (f->h[0])*(f->hq[a][0])*(f->hq[0][a])
		+ (f->h[0])*(f->h[a])*(f->hqq[0][a][0]) );

	  f->d_grad_e_dq[a][1][1][1] = ( 1. / ( (f->h[1])*(f->h[1]) *(f->h[a])*(f->h[a]) ) )
	    * ( - (f->h[a])*(f->hq[1][1])*(f->hq[1][a])
		- (f->h[1])*(f->hq[a][1])*(f->hq[1][a])
		+ (f->h[1])*(f->h[a])*(f->hqq[1][a][1]) );

	  f->d_grad_e_dq[a][0][0][1] = ( 1. / ( (f->h[0])*(f->h[0]) *(f->h[a])*(f->h[a]) ) )
	    * ( - (f->h[a])*(f->hq[0][1])*(f->hq[0][a])
		- (f->h[0])*(f->hq[a][1])*(f->hq[0][a])
		+ (f->h[0])*(f->h[a])*(f->hqq[0][a][1]) );
			   
	  f->d_grad_e_dq[a][1][1][0] = ( 1. / ( (f->h[1])*(f->h[1]) *(f->h[a])*(f->h[a]) ) )
	    * ( - (f->h[a])*(f->hq[1][0])*(f->hq[1][a])
		- (f->h[1])*(f->hq[a][0])*(f->hq[1][a])
		+ (f->h[1])*(f->h[a])*(f->hqq[1][a][0]) );
		
	  if( VIM == 3) 
	    {
	      f->d_grad_e_dq[a][2][2][2] = ( 1. / ( (f->h[2])*(f->h[2]) *(f->h[a])*(f->h[a]) ) )
		* ( - (f->h[a])*(f->hq[2][2])*(f->hq[2][a])
		    - (f->h[2])*(f->hq[a][2])*(f->hq[2][a])
		    + (f->h[2])*(f->h[a])*(f->hqq[2][a][2]) );

	      f->d_grad_e_dq[a][2][2][0] = ( 1. / ( (f->h[2])*(f->h[2]) *(f->h[a])*(f->h[a]) ) )
		* ( - (f->h[a])*(f->hq[2][0])*(f->hq[2][a])
		    - (f->h[2])*(f->hq[a][0])*(f->hq[2][a])
		    + (f->h[2])*(f->h[a])*(f->hqq[2][a][0]) );
			   
	      f->d_grad_e_dq[a][2][2][1] = ( 1. / ( (f->h[2])*(f->h[2]) *(f->h[a])*(f->h[a]) ) )
		* ( - (f->h[a])*(f->hq[2][1])*(f->hq[2][a])
		    - (f->h[2])*(f->hq[a][1])*(f->hq[2][a])
		    + (f->h[2])*(f->h[a])*(f->hqq[2][a][1]) );

	      f->d_grad_e_dq[a][1][1][2] = ( 1. / ( (f->h[1])*(f->h[1]) *(f->h[a])*(f->h[a]) ) )
		* ( - (f->h[a])*(f->hq[1][2])*(f->hq[1][a])
		    - (f->h[1])*(f->hq[a][2])*(f->hq[1][a])
		    + (f->h[1])*(f->h[a])*(f->hqq[1][a][2]) );

	      f->d_grad_e_dq[a][0][0][2] = ( 1. / ( (f->h[0])*(f->h[0]) *(f->h[a])*(f->h[a]) ) )
		* ( - (f->h[a])*(f->hq[0][2])*(f->hq[0][a])
		    - (f->h[0])*(f->hq[a][2])*(f->hq[0][a])
		    + (f->h[0])*(f->h[a])*(f->hqq[0][a][2]) );
	    }
#endif
	}
  
      for ( a=0; a<VIM; a++)
	{
#ifndef DO_NO_UNROLL
	  for (q=0; q<VIM; q++)
	    {
	      for ( b=0; b<VIM; b++)
		{
		  f->d_grad_e_dq[a][a][q][b] -= ( 1. / ( (f->h[a])*(f->h[a]) *(f->h[q]) * (f->h[q]) ) ) 
		    * ( - (f->h[q])*(f->hq[a][b])*(f->hq[a][q])
			- (f->h[a])*(f->hq[q][b])*(f->hq[a][q])
			+ (f->h[a])*(f->h[q])*(f->hqq[a][q][b]) );
		}
	    }
#else

	  f->d_grad_e_dq[a][a][0][0] -= ( 1. / ( (f->h[a])*(f->h[a]) *(f->h[0]) * (f->h[0]) ) ) 
	    * ( - (f->h[0])*(f->hq[a][0])*(f->hq[a][0])
		- (f->h[a])*(f->hq[0][0])*(f->hq[a][0])
		+ (f->h[a])*(f->h[0])*(f->hqq[a][0][0]) );
			   
	  f->d_grad_e_dq[a][a][1][1] -= ( 1. / ( (f->h[a])*(f->h[a]) *(f->h[1]) * (f->h[1]) ) ) 
	    * ( - (f->h[1])*(f->hq[a][1])*(f->hq[a][1])
		- (f->h[a])*(f->hq[1][1])*(f->hq[a][1])
		+ (f->h[a])*(f->h[1])*(f->hqq[a][1][1]) );
			   
	  f->d_grad_e_dq[a][a][0][1] -= ( 1. / ( (f->h[a])*(f->h[a]) *(f->h[0]) * (f->h[0]) ) ) 
	    * ( - (f->h[0])*(f->hq[a][1])*(f->hq[a][0])
		- (f->h[a])*(f->hq[0][1])*(f->hq[a][0])
		+ (f->h[a])*(f->h[0])*(f->hqq[a][0][1]) );
		
	  f->d_grad_e_dq[a][a][1][0] -= ( 1. / ( (f->h[a])*(f->h[a]) *(f->h[1]) * (f->h[1]) ) ) 
	    * ( - (f->h[1])*(f->hq[a][0])*(f->hq[a][1])
		- (f->h[a])*(f->hq[1][0])*(f->hq[a][1])
		+ (f->h[a])*(f->h[1])*(f->hqq[a][1][0]) );
	
	  if (VIM == 3)
	    {

	      f->d_grad_e_dq[a][a][2][2] -= ( 1. / ( (f->h[a])*(f->h[a]) *(f->h[2]) * (f->h[2]) ) ) 
		* ( - (f->h[2])*(f->hq[a][2])*(f->hq[a][2])
		    - (f->h[a])*(f->hq[2][2])*(f->hq[a][2])
		    + (f->h[a])*(f->h[2])*(f->hqq[a][2][2]) );

	      f->d_grad_e_dq[a][a][2][0] -= ( 1. / ( (f->h[a])*(f->h[a]) *(f->h[2]) * (f->h[2]) ) ) 
		* ( - (f->h[2])*(f->hq[a][0])*(f->hq[a][2])
		    - (f->h[a])*(f->hq[2][0])*(f->hq[a][2])
		    + (f->h[a])*(f->h[2])*(f->hqq[a][2][0]) );

	      f->d_grad_e_dq[a][a][2][1] -= ( 1. / ( (f->h[a])*(f->h[a]) *(f->h[2]) * (f->h[2]) ) ) 
		* ( - (f->h[2])*(f->hq[a][1])*(f->hq[a][2])
		    - (f->h[a])*(f->hq[2][1])*(f->hq[a][2])
		    + (f->h[a])*(f->h[2])*(f->hqq[a][2][1]) );

	      f->d_grad_e_dq[a][a][1][2] -= ( 1. / ( (f->h[a])*(f->h[a]) *(f->h[1]) * (f->h[1]) ) ) 
		* ( - (f->h[1])*(f->hq[a][2])*(f->hq[a][1])
		    - (f->h[a])*(f->hq[1][2])*(f->hq[a][1])
		    + (f->h[a])*(f->h[1])*(f->hqq[a][1][2]) );
	
	      f->d_grad_e_dq[a][a][0][2] -= ( 1. / ( (f->h[a])*(f->h[a]) *(f->h[0]) * (f->h[0]) ) ) 
		* ( - (f->h[0])*(f->hq[a][2])*(f->hq[a][0])
		    - (f->h[a])*(f->hq[0][2])*(f->hq[a][0])
		    + (f->h[a])*(f->h[0])*(f->hqq[a][0][2]) );
	
	    }
#endif

	}
    }
  
  return(status);
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

double 
global_velocity_norm(const dbl x[], /* solution vector */
		     const Exo_DB *exo,	/* ptr to FE db */
		     const Dpi *dpi) /* ptr to dist. info */
    
     /*******************************************************************
      *
      * global_velocity_norm():
      *
      *  Obtain some average velocity for the problem
      *  Crudely speaking, e.g., in 3D, provide the square of the
      *  L2 norm of the speed:
      *
      *	<U> = 1/(3*Num_Nodes) * Sum over nodes(Vx^2 + Vy^2 + Vz^2)
      *
      *  that can be used in a PSPG formulation.
      *  Note: which material to gather the velocity is not specified.
      *        Therefore, we take the  "first" velocity specififed at
      *        a node, irrespective of what material it belongs to.
      *        In the future, it may be better to gather an average
      *        velocity for each material in the problem, separately,
      *        for application to PSPG formulation.
      *
      *  Modified to communicate across distributed processors using MPI.
      *
      *  Created: ~1995 rrrao
      *  Revised: 1997/08/25 08:39 MDT pasacki@sandia.gov
      ********************************************************************/
{
  int norm_unknowns = 0; /* norm_unknowns is the actual number of
			    variables used in the norm calculation  */
#ifdef PARALLEL  
  int global_norm_unknowns;	/* sum over all processors */
  double global_U_norm;		/* average over all processors */
#endif
  int eqn, ie, i;
  double U_norm = 0.0;
  NODE_INFO_STRUCT *node;

  for (eqn = R_MOMENTUM1; eqn <= R_MOMENTUM3; eqn++) {
    if (Num_Var_In_Type[pg->imtrx][eqn]) {
      for (i = 0;
	   i < (dpi->num_internal_nodes + dpi->num_boundary_nodes);
	   i++) {
	node = Nodes[i];
	if (get_nv_ndofs(node->Nodal_Vars_Info[pg->imtrx], eqn)) {
	  ie = Index_Solution(i, eqn, 0, 0, -2, pg->imtrx);
	  if (ie != -1) {
	    U_norm += x[ie] * x[ie];
	    norm_unknowns += 1;
	  }
	}
      }
    }  
  }
#ifdef PARALLEL
  /*
   * Sum up the U_norm contribution from all processors,
   * likewise sum up the actual number of nodal contributions
   * from each processor.
   */
  MPI_Allreduce(&U_norm, &global_U_norm, 1, MPI_DOUBLE, MPI_SUM, 
		MPI_COMM_WORLD);
  MPI_Allreduce(&norm_unknowns, &global_norm_unknowns, 1,
		MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  U_norm        = global_U_norm;
  norm_unknowns = global_norm_unknowns;
#endif  
  U_norm /= (double) norm_unknowns;
  return(U_norm);
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

dbl element_viscosity(void)

     /*****************************************************************
      *
      *  element_viscosity()
      *
      * For now, this routine only returns the zero
      * shear-rate viscosity. It will need to be
      * enhanced for shear-thinning, pressure dependent
      * fluids etc. That will have complicated derivatives.
      *****************************************************************/
{
  int i, a, b, p, q, mode;
  dbl mu_avg;
  dbl gamma[DIM][DIM];  /* shrearrate tensor based on velocity */

  /*
   * Variables for vicosity and derivative
   */
  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct;  /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;

  /* polymer viscosity and derivatives */

  dbl mup;
  VISCOSITY_DEPENDENCE_STRUCT d_mup_struct;  /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mup = &d_mup_struct;

  /* load up shearrate tensor based on nonzero initialization */
  for ( a=0; a<VIM; a++)
    {
      for ( b=0; b<VIM; b++)
	{
	  gamma[a][b] = 1.e-5;
	}
    }
  
  /* load something into fill variable */
  fv->F = .5;
  
  if( ls != NULL  ) fv->F = 0.0;
  
  fv->T = mp->reference[TEMPERATURE];
  for (i = 0; i < pd->Num_Species_Eqn; i++) 
      fv->c[i]=mp->reference_concn[i];
  for (i = 0; i < ei[pg->imtrx]->dof[MASS_FRACTION]; i++)
    bf[MASS_FRACTION]->phi[i]=0.0;  

  /* initialize grad_phi_e for shear-thinning models */  
  for ( i=0; i<ei[pg->imtrx]->dof[VELOCITY1]; i++)
    {
      for ( p=0; p<VIM; p++)
	{
	  for ( q=0; q<VIM; q++)
	    {
	      for ( a=0; a<VIM; a++)
		{
		  bf[VELOCITY1]->grad_phi_e[i][a] [p][q] = 0.;
		}
	    }
	}
    }
  
  mu_avg = viscosity(gn, gamma, d_mu);

  /* get polymer viscosity */
  if ( pd->v[pg->imtrx][POLYMER_STRESS11] )
    {
      for ( mode=0; mode<vn->modes; mode++)
	{
	  mup = viscosity(ve[mode]->gn, gamma, d_mup);
	  mu_avg += mup;
	}
    }
  return(mu_avg);
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

void 
element_velocity(dbl v_avg[DIM], dbl dv_dnode[DIM][MDE],
		 const Exo_DB *exo)

    /*
     * this routine is used to calculate an average element velocity
     * for SUPG.  Current apps are for the energy equation. 
     */
{
  int i, p, I, dofs, centroid_node, dim;
  dbl ddofs;
  
  dim =  pd->Num_Dim;
  
  /* parameter variables are initialized in matrix_fill */
  
  if (pd->i[pg->imtrx][VELOCITY1]==I_Q1)
    {
      if (cr->MeshMotion == ARBITRARY) {
	  for (p = 0; p < dim; p++)
	    {

	      dofs     = ei[pg->imtrx]->dof[VELOCITY1];
	      ddofs = dofs;  

	      for (i = 0; i < dofs; i++)
		{
		  v_avg[p] += *esp->v[p][i]/ddofs;
		  dv_dnode[p][i] = 1.0/ddofs;
		}
	    }
	}
      else
	{
	  /* use the velocity of the stress-free-state */
	  dofs     = ei[pg->imtrx]->dof[MESH_DISPLACEMENT1];
	  ddofs = dofs;  
	  for (i=0; i<dofs; i++) 
	    {
	      I = exo->node_list[ei[pg->imtrx]->iconnect_ptr + i];

	      if ((cr->MeshMotion == LAGRANGIAN ||
		   cr->MeshMotion == DYNAMIC_LAGRANGIAN))
		{
		  if (elc->v_mesh_sfs_model == ROTATIONAL ||
			elc->v_mesh_sfs_model == ROTATIONAL_3D)
		    {
		      (void) V_mesh_sfs_model(elc->u_v_mesh_sfs, 
						elc->v_mesh_sfs, 
						elc->v_mesh_sfs_model, I);
		    }
		} 
	      else if (cr->MeshMotion == TOTAL_ALE)
		{
		  if (elc_rs->v_mesh_sfs_model == ROTATIONAL ||
			elc_rs->v_mesh_sfs_model == ROTATIONAL_3D)
		  {
		    (void) V_mesh_sfs_model(elc_rs->u_v_mesh_sfs,
					    elc_rs->v_mesh_sfs, 
					    elc_rs->v_mesh_sfs_model, I);
		  }
		}
	      for( p=0; p<dim; p++)		
		{  
		  v_avg[p] += elc->v_mesh_sfs[p]/ddofs;
		  dv_dnode[p][i] = 0.0;
		}
	    }

	}
	  
    }
  else if (pd->i[pg->imtrx][VELOCITY1]==I_Q2)
    {    
      if ( cr->MeshMotion == ARBITRARY)
	{
	  dofs     = ei[pg->imtrx]->dof[VELOCITY1];
	  centroid_node = dofs-1;

	  for( p=0; p<dim; p++)
	    {
	      v_avg[p] =  *esp->v[p][centroid_node];
	      dv_dnode[p][centroid_node] = 1.0;
	    }
	}
      else
	{
	  dofs     = ei[pg->imtrx]->dof[MESH_DISPLACEMENT1];
	  centroid_node = dofs-1;

	  I = exo->node_list[ei[pg->imtrx]->iconnect_ptr + centroid_node];
	  if ((cr->MeshMotion == LAGRANGIAN ||
	       cr->MeshMotion == DYNAMIC_LAGRANGIAN))
	    {
	      if (elc->v_mesh_sfs_model == ROTATIONAL ||
			elc->v_mesh_sfs_model == ROTATIONAL_3D)
	      {
		(void) V_mesh_sfs_model(elc->u_v_mesh_sfs, 
					elc->v_mesh_sfs, 
					elc->v_mesh_sfs_model, I);
	      }
	    } 
	  else if (cr->MeshMotion == TOTAL_ALE)
	    {
	      if (elc_rs->v_mesh_sfs_model == ROTATIONAL ||
			elc_rs->v_mesh_sfs_model == ROTATIONAL_3D)
	      {
		(void) V_mesh_sfs_model(elc_rs->u_v_mesh_sfs,
	 				    elc_rs->v_mesh_sfs, 
					    elc_rs->v_mesh_sfs_model, I);
	      }
	    }
	  for( p=0; p<dim; p++)
	    {
	      v_avg[p] =  elc->v_mesh_sfs[p];
	      dv_dnode[p][centroid_node] = 0.0;
	    }
	}
    }
} /* end of v_elem_avg */
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

void 
h_elem_siz(dbl hsquared[DIM], dbl hh[DIM][DIM],
	   dbl dhh_dxnode[DIM][MDE], const int DeformingMesh)

    /*********************************************************************
     * h_elem_siz():
     *
     * This routine calculates the approximate size of the element
     * in the local coordinate directions for use in SUPG stabilization
     * of advection dominated problems. It is
     * also used to determine an overall length scale for an element
     * for PSPG schemes. 
     *
     * The heuristics of this function are only
     * designed for quad/hex element types. If other elements
     * are used, a new algorithm must be developed.
     *
     * Output
     * ------------
     *   hsquared[p] = This is the square of the width of the element
     *                 along a local element  coordinate direction, 
     *                 measured from the center of the 
     *                 appropriate faces of the element.
     *   hh[p][i]    = This is the difference in the ith real coordinate
     *                 values of the face centroids along the pth 
     *                 local element direction.
     *                 hh[p] can also be considered to be a vector
     *                 oriented in the pth local coodinate direction.
     *                 When divided by the square root of hsquared[p],
     *                 it provides a unit vector oriented in the pth
     *                 local element coordinate direction.
     *   dhh_dxnode[p][j] = Derivative of the value hh[p][i] due to the
     *                 change in the mesh position of the jth basis
     *                 function along the ith real coordinate direction.
     *
     *********************************************************************/
{
  int i, j, p, idof, dim, index, var;
  dbl xnode[DIM][MDE];
  dbl p1[DIM], p2[DIM], p3[DIM], p4[DIM], p5[DIM], p6[DIM];
  dim =  pd->Num_Dim;
  int elem_type = ei[pg->imtrx]->ielem_type;
  int elem_shape = type2shape(elem_type);
  int DeformingMeshShell = 0;

  /* initialize xnode */
  for (i = 0; i < DIM; i++) {
    for (j = 0; j < MDE; j++) {
      xnode[i][j] = 0;
    }
  }

  /*
   * Find the local coodinates of the nodes in the current element
   * -> Can this be pushed to a more generic place ?
   */

  if (mp->FSIModel == FSI_MESH_CONTINUUM ||
      mp->FSIModel == FSI_REALSOLID_CONTINUUM ||
      mp->FSIModel == FSI_SHELL_ONLY_MESH) DeformingMeshShell = 1;

  j = Proc_Connect_Ptr[ei[pg->imtrx]->ielem];
  for (p = 0; p < dim; p++) {
    if (DeformingMesh || DeformingMeshShell) {
      var = MESH_DISPLACEMENT1 + p;
      for (i = 0; i < ei[pd->mi[var]]->num_local_nodes; i++) {
	idof = ei[pd->mi[var]]->ln_to_dof[var][i];
	index = Proc_Elem_Connect[j + i];
	if(idof != -1) xnode[p][i]= Coor[p][index] + *esp->d[p][idof];
      }
    } else {
      for (i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
	index = Proc_Elem_Connect[j  + i];
	xnode[p][i]= Coor[p][index];
      }
    }
  }
  
  if (dim == 2 && (elem_shape == TRIANGLE)) {
    WH(-1,"Beware that SUPG for trishells is held constant so you need a uniform, non-stretchnig grid");

    hsquared[0]=hsquared[1]=hsquared[2] = pow((xnode[0][0] - xnode[0][1]), 2) + pow((xnode[1][0] - xnode[1][1]),2);


    if (af->Assemble_Jacobian) {
      memset((void *)dhh_dxnode, 0, dim*MDE*sizeof(double));
      if (DeformingMesh || DeformingMeshShell) {
        dhh_dxnode[0][0] = -0.0;
        dhh_dxnode[0][1] =  0.0;
        dhh_dxnode[0][2] =  0.0;

        dhh_dxnode[1][0] =  0.0;
        dhh_dxnode[1][1] =  0.0;
        dhh_dxnode[1][2] = -0.0;

        dhh_dxnode[2][0] =  0.0;
        dhh_dxnode[2][1] =  0.0;
        dhh_dxnode[2][2] = -0.0;

      }
    }
  }
  else if (dim == 2) {
    /*
     * Calculate the midpoint positions on each of the faces
     */
    for (p = 0; p < dim; p++) {
      p1[p] = 0.5 * (xnode[p][0] + xnode[p][3]);
      p2[p] = 0.5 * (xnode[p][1] + xnode[p][2]);
      p3[p] = 0.5 * (xnode[p][0] + xnode[p][1]);
      p4[p] = 0.5 * (xnode[p][2] + xnode[p][3]);
    }
      
    hh[0][0] = (p2[0] - p1[0]);
    hh[0][1] = (p2[1] - p1[1]);
    hsquared[0] = hh[0][0] * hh[0][0] + hh[0][1] * hh[0][1];
 
    hh[1][0] = (p3[0] - p4[0]);
    hh[1][1] = (p3[1] - p4[1]);
    hsquared[1] = hh[1][0] * hh[1][0] + hh[1][1] * hh[1][1];
    if (af->Assemble_Jacobian) {
      memset((void *)dhh_dxnode, 0, dim*MDE*sizeof(double));
      if (DeformingMesh) {
	dhh_dxnode[0][0] = -0.5;
	dhh_dxnode[0][1] =  0.5;
	dhh_dxnode[0][2] =  0.5;
	dhh_dxnode[0][3] = -0.5;
	  
	dhh_dxnode[1][0] =  0.5;
	dhh_dxnode[1][1] =  0.5;
	dhh_dxnode[1][2] = -0.5;
	dhh_dxnode[1][3] = -0.5;
      }
    }
  } else if (dim == 3 && (elem_shape == SHELL)) {

     /*Special Case */
    /*
     * Calculate the midpoint positions on each of the faces
     */
    for (p = 0; p < dim; p++) {
      p1[p] = 0.5 * (xnode[p][0] + xnode[p][3]);
      p2[p] = 0.5 * (xnode[p][1] + xnode[p][2]);
      p3[p] = 0.5 * (xnode[p][0] + xnode[p][1]);
      p4[p] = 0.5 * (xnode[p][2] + xnode[p][3]);
    }
      
    hh[0][0] = (p2[0] - p1[0]);
    hh[0][1] = (p2[1] - p1[1]);
    hh[0][2] = (p2[2] - p1[2]);
    hsquared[0] = hh[0][0] * hh[0][0] + hh[0][1] * hh[0][1] + hh[0][2] * hh[0][2];
 
    hh[1][0] = (p3[0] - p4[0]);
    hh[1][1] = (p3[1] - p4[1]);
    hh[1][2] = (p3[2] - p4[2]);
    hsquared[1] = hh[1][0] * hh[1][0] + hh[1][1] * hh[1][1] + hh[1][2] * hh[1][2];


    if (af->Assemble_Jacobian) {
      memset((void *)dhh_dxnode, 0, dim*MDE*sizeof(double));
      if (DeformingMesh || DeformingMeshShell) {
	dhh_dxnode[0][0] = -0.5;
	dhh_dxnode[0][1] =  0.5;
	dhh_dxnode[0][2] =  0.5;
	dhh_dxnode[0][3] = -0.5;
	  
	dhh_dxnode[1][0] =  0.5;
	dhh_dxnode[1][1] =  0.5;
	dhh_dxnode[1][2] = -0.5;
	dhh_dxnode[1][3] = -0.5;

      }
    }
  } else if (dim == 3 && (elem_shape == TRISHELL || elem_shape == TETRAHEDRON)) {

    /*
     * Use a constant for now between local nodes 1 and 2 
     */
    WH(-1,"\nBeware that SUPG for trishells/tetrahedrons is held constant so...\nYou need a uniform, non-stretching grid\n");

    hsquared[0]=hsquared[1]=hsquared[2] = pow((xnode[0][0] - xnode[0][1]), 2) + pow((xnode[1][0] - xnode[1][1]), 2) + pow((xnode[2][0] - xnode[2][1]), 2);


    if (af->Assemble_Jacobian) {
      memset((void *)dhh_dxnode, 0, dim*MDE*sizeof(double));
      if (DeformingMesh || DeformingMeshShell) {
	dhh_dxnode[0][0] = -0.0;
	dhh_dxnode[0][1] =  0.0;
	dhh_dxnode[0][2] =  0.0;
	  
	dhh_dxnode[1][0] =  0.0;
	dhh_dxnode[1][1] =  0.0;
	dhh_dxnode[1][2] = -0.0;

	dhh_dxnode[2][0] =  0.0;
	dhh_dxnode[2][1] =  0.0;
	dhh_dxnode[2][2] = -0.0;

      }
    }
  } else if (dim == 3 && elem_shape != TETRAHEDRON) {

    for (p = 0; p < dim; p++) {
      p1[p] = 0.25 * (xnode[p][0] + xnode[p][1] + xnode[p][2] + xnode[p][3]);
      p2[p] = 0.25 * (xnode[p][1] + xnode[p][2] + xnode[p][5] + xnode[p][6]);
      p3[p] = 0.25 * (xnode[p][2] + xnode[p][3] + xnode[p][6] + xnode[p][7]);
      p4[p] = 0.25 * (xnode[p][0] + xnode[p][1] + xnode[p][4] + xnode[p][5]);
      p5[p] = 0.25 * (xnode[p][0] + xnode[p][3] + xnode[p][4] + xnode[p][7]);
      p6[p] = 0.25 * (xnode[p][4] + xnode[p][5] + xnode[p][6] + xnode[p][7]);
    }
    for (p = 0; p < dim; p++) {
      hh[0][p] = (p2[p] - p5[p]);
      hh[1][p] = (p3[p] - p4[p]);	  
      hh[2][p] = (p1[p] - p6[p]);
    }
    hsquared[0] = hh[0][0]*hh[0][0] + hh[0][1]*hh[0][1] + hh[0][2]*hh[0][2];
    hsquared[1] = hh[1][0]*hh[1][0] + hh[1][1]*hh[1][1] + hh[1][2]*hh[1][2];
    hsquared[2] = hh[2][0]*hh[2][0] + hh[2][1]*hh[2][1] + hh[2][2]*hh[2][2];
    if (af->Assemble_Jacobian) {      
      memset((void *)dhh_dxnode, 0, (dim*MDE*sizeof(double)));
      if (DeformingMesh) {
	dhh_dxnode[0][0] = -0.25;
	dhh_dxnode[0][1] =  0.25;
	dhh_dxnode[0][2] =  0.25;
	dhh_dxnode[0][3] = -0.25;
	dhh_dxnode[0][4] = -0.25;
	dhh_dxnode[0][5] =  0.25;
	dhh_dxnode[0][6] =  0.25;
	dhh_dxnode[0][7] = -0.25;
	  
	dhh_dxnode[1][0] = -0.25;
	dhh_dxnode[1][1] = -0.25;
	dhh_dxnode[1][2] =  0.25;
	dhh_dxnode[1][3] =  0.25;
	dhh_dxnode[1][4] = -0.25;
	dhh_dxnode[1][5] = -0.25;
	dhh_dxnode[1][6] =  0.25;
	dhh_dxnode[1][7] =  0.25;

	dhh_dxnode[2][0] =  0.25;
	dhh_dxnode[2][1] =  0.25;
	dhh_dxnode[2][2] =  0.25;
	dhh_dxnode[2][3] =  0.25;
	dhh_dxnode[2][4] = -0.25;
	dhh_dxnode[2][5] = -0.25;
	dhh_dxnode[2][6] = -0.25;
	dhh_dxnode[2][7] = -0.25;
      }
    }
  }
  else
    {
      EH(-1,"SUPG not allowed for tetrahedral elements yet, or whatever weird element you have");
    }
}
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

void 
const_h_elem_siz(dbl h[DIM], dbl hh[DIM][DIM], 
		 dbl dhh_dxnode[DIM][MDE])
 
     /*
      * this routine returns a constant element size of 1 and the
      * PSPG scaling term in the input deck should be used to indicate the 
      * average element length. Later, it would be nice to add a routine
      * that calculates the average mesh element size at each iteration.
      */
{
  int p, q, dim = pd->Num_Dim;
  for (p = 0; p < dim; p++) {
    h[p] = 1.;
    for(q = 0; q < dim; q++) {
      hh[p][q] = 1.; 
    }
    for (q = 0; q < MDE; q++) {
      dhh_dxnode[p][q] = 0.;
    }
  }
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

dbl 
global_h_elem_siz(dbl x[], dbl x_old[], dbl xdot[], dbl resid_vector[],
		  Exo_DB *exo, Dpi *dpi)

/* global_h_elem_siz -- calculates  global average element length
 *
 * For the SUPG terms in the EVSS formulation and for the PSPG terms in the
 * GLS stabilized formulation it is necessary to come up with some kind of
 * representative element size. An average over the elements is performed.
 *
 * The distributed processing version does not calculate the precise
 * global average. Rather, it computes contributions from every private element
 * and (1/2) contributions from every shared element. Since some shared 
 * elements may be shared by more than two processors, those elements will 
 * contribute to the average more than your typical element. The hope is that 
 * this will not upset your stabilized scheme all that much.
 *
 * Epilogue: Sam reports that it DOES matter, that observable differences in
 *           residuals, particularly for the continuity equation, result in
 *           parallel simulations that are evidently attributable to using
 *           different values of tau_pspg that, in turn, result from this
 *           h_elem_size being different.
 *
 *           Fortunately, with the need of the discontinuous Galerkin methods
 *           for a unique element decomposition and assigment, the distributed
 *           version of GOMA can make use of that assigment so that the
 *           average element size in parallel should be identical with that
 *           in serial.
 *
 * Written: rrrao
 *
 * Revised: 1997/08/19 09:26 MDT pasacki@sandia.gov
 *
 * Revised: 1999/10/01 07:48 MDT pasacki@sandia.gov
 */
{
  int e, dim, p;
  dbl h, h_elem, hsquared[DIM], hhv[DIM][DIM], dhv_dxnode[DIM][MDE];
  dbl weight;                   /* 1 usually, except for multiprocessing */
  dbl smele_mun;                /* (1/num_elems) */
#ifdef PARALLEL
  dbl h_global_really_i_mean_it; /* You have any doubt, now? */
#endif
  
  dim =  pd->Num_Dim;
  /*
   * Default values for weighting and number of global elements for serial
   * processing.
   */
  smele_mun = 1.0 / (double) dpi->num_elems_global;

  /*
   * Look through each element block, then at each element in each block.
   *
   * If the element is owned by this processor, then get it's contribution
   * to h_elem_size.
   *
   */
  h = 0.0;
  for (e = 0; e < exo->num_elems; e++) {
    (void) load_elem_dofptr(e, exo , x, x_old, xdot, xdot, 1);
    if(ei[pg->imtrx]->ielem_dim != 1) h_elem_siz(hsquared, hhv, dhv_dxnode, 0);
    h_elem = 0.;
    for (p = 0; p < ei[pg->imtrx]->ielem_dim; p++) {
      h_elem += hsquared[p];
    }
    h_elem = sqrt(h_elem/dim); 
    weight = 1.;
    if (Num_Proc > 1) {
      if ( dpi->elem_owner[e] != ProcID) {
	weight = 0.0;
      }
    }
    h += weight * h_elem * smele_mun;
  }
#ifdef PARALLEL
  MPI_Allreduce(&h, &h_global_really_i_mean_it, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  h = h_global_really_i_mean_it;
#endif
  return(h);
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
void

surface_determinant_and_normal(
 const int ielem,		/* current element number               */
 const int iconnect_ptr,	/* Pointer to beginning of connectivity
				 * list for current element             */
 const int nodes_per_elem,	/* number of nodes in the element       */
 const int ielem_surf_dim,	/* physical dimension of the element
				 * surface (0, 1, 2)                    */
 const int id_side,		/* ID of element side (exo/patran convention)  */
 const int num_nodes_on_side,	/* number of nodes on side of element   */
 const int local_elem_node_id [] ) /* vector of local element node numbers
				    * on the side of the element           */

    /************************************************************************
     *
     * surface_determinant_and_normal()
     *
     *      Function which calculates the surface determinant and surface
     *      normal at a local surface quadrature point. The function also
     *      calculates the sensitivities of those quantities to the
     *      mesh positions, if the mesh positions are part of the solution
     *      vector.
     *
     *
     * Insulate for more robust behavior for fixed grid problems - 950306pas
     *
     * Incorporate coordinate scale factors into the surface determinant
     * and combine with normal calculation
     *               Richard Cairncross 9/24/96
     *
     * Add tetrahedral elements and fix hex elements. - SAR 2011-11-03
     *
     *  Returns:
     * ------------
     *   fv->sdet = surface determinant at the quadrature point
     *   fv->snormal[] = surface normal at the quadrature point
     *   fv->dsurfdet_dx[][] = sensitivity of fv->sdet wrt mesh displacements.
     *   fv-.dsnormal_dx[][] = sensitivity of fv->snormal[]
     *                         wrt mesh displacements
     ***********************************************************************/
{
  /* TAB certifies that this function conforms to the exo/patran side
     numbering convention 11/10/98. */

  int 		i, id, inode, a, b, p, q;
  int		ShapeVar, ldof;
  int		DeformingMesh;
  double        r_det, det_h01, r_det_h01, d_det_h01_x;
  double        phi_i;
  int siz;
  double        T[DIM-1][DIM], t[DIM-1][DIM];  /* t = J . T */
  double        dt_x[DIM-1][DIM][DIM][MDE]; /* d(t) / d(x_j) */
  struct Basis_Functions *map_bf;
  double        signID;
  double tmp;
  int dim =  pd->Num_Dim;

  DeformingMesh = pd->gv[R_MESH1];
  ShapeVar = pd->ShapeVar;

  siz = MAX_PDIM*MDE*sizeof(double);
  memset(fv->dsurfdet_dx,0,siz);
  siz = MAX_PDIM*MDE*MAX_PDIM*sizeof(double);
  memset(fv->dsnormal_dx,0,siz);

  /*
  si = in_list(pd->IntegrationMap, 0, Num_Interpolations, Unique_Interpolations);
  map_bf = bfd[si];
  */

  map_bf = bf[ShapeVar];

  if (ielem_surf_dim == 0)
    {
      /*
       *  This code should be considered to be untested. It worked for two cases.
       */
      // ok we fill up snormal
      shell_determinant_and_normal(ei[pg->imtrx]->ielem, ei[pg->imtrx]->iconnect_ptr, ei[pg->imtrx]->num_local_nodes,
				   ei[pg->imtrx]->ielem_dim, 1);
      if (id_side == 1) {
	      signID = -1.0;
      } else if (id_side == 2) {
        signID = 1.0;
      } else {
	      signID = 0.0;
      }
      /*
       * We turn the normal into a tangent and then assign it back, multiplying by the signID
       *       tangent[0] = -normal[1]
       *       tangent[1] =  normal[0]
       *
       *  In 2D n x t = k, these vectors satisfy the rh rule.
       */

      tmp =   fv->snormal[0];
      fv->snormal[0] = - signID *  fv->snormal[1]; 
      fv->snormal[1] =   signID *  tmp; 

      for (i = 0; i < ei[pg->imtrx]->num_local_nodes; i++)
	{
	  inode = Proc_Elem_Connect[iconnect_ptr + i];
	  ldof = ei[pg->imtrx]->ln_to_dof[ShapeVar][i];
	  if (ldof >= 0)
	    {
	      if (Dolphin[pg->imtrx][inode][MESH_DISPLACEMENT1] > 0 )
              {
                for (a = 0; a < dim; a++)
                {
                  tmp =  fv->dsnormal_dx[0][a][ldof];
                  fv->dsnormal_dx[0][a][ldof] = - signID * fv->dsnormal_dx[1][a][ldof];
                  fv->dsnormal_dx[1][a][ldof] =   signID * tmp;
                }
              }
            }
        }
      fv->sdet = 1.0;
      memset(fv->dsurfdet_dx, 0.0, sizeof(double)*DIM*MDE);
      return;
    }

  /* define space of surface */  
  switch (ielem_surf_dim) {
  case 1:
    switch (ei[pg->imtrx]->ielem_shape) {
    case TRIANGLE:
    case TRISHELL:
      if ( id_side == 1 )
        {
          /*
          T[0][0] = -Ref; T[0][1] = Ref;
          */
          T[0][0] = -1.; T[0][1] = 1.;
        }
      else if ( id_side == 2 )
        {
          T[0][0] = 0.; T[0][1] = -1.;
        }
      else if ( id_side == 3 )
        {
          T[0][0] = 1.; T[0][1] = 0.;
        }
      else
        {
          EH(-1, "Incorrect side for TRIANGLE");
        }
      break;
    case QUADRILATERAL:
    case SHELL:
      if ( id_side == 1 )
        {
          T[0][0] = 1.; T[0][1] = 0.;
        }
      else if ( id_side == 2 )
        {
          T[0][0] = 0.; T[0][1] = 1.;
        }
      else if ( id_side == 3 )
        {
          T[0][0] = -1.; T[0][1] = 0.;
        }
      else if ( id_side == 4 )
        {
          T[0][0] = 0.; T[0][1] = -1.;
        }
      else
        {
          EH(-1, "Incorrect side for QUADRILATERAL");
        }
      break;
    default:
      EH(-1, "Invalid shape");
    }
    break;

  case 2:
    switch (ei[pg->imtrx]->ielem_shape) {
    case HEXAHEDRON:
      /*
       * In trying to figure out how these T matrices work, I think that
       * this calculation is incorrect for sides 4 and 5.  My modifications
       * are below.  Making this change did not affect the area calculation,
       * normal calculation, and sdet calculation in my test problem, 
       * and the test suite passed, so I'll keep it for now.
       * Scott A Roberts, 2011-11-03
       */
      if ( id_side == 1 )
        {
          T[0][0] =  1.; T[0][1] =  0.; T[0][2] =  0.;
          T[1][0] =  0.; T[1][1] =  0.; T[1][2] =  1.;
        }
      else if (id_side == 2)
        {
          T[0][0] =  0.; T[0][1] =  1.; T[0][2] =  0.;
          T[1][0] =  0.; T[1][1] =  0.; T[1][2] =  1.;
        }
      else if (id_side == 3)
        {
          T[0][0] = -1.; T[0][1] =  0.; T[0][2] =  0.;
          T[1][0] =  0.; T[1][1] =  0.; T[1][2] =  1.;
        }
      else if (id_side == 4)
        {
	  //T[0][0] =  0.; T[0][1] = -1.; T[0][2] =  0.;
	  //T[1][0] =  0.; T[1][1] =  0.; T[1][2] =  1.;
          T[0][0] =  0.; T[0][1] =  0.; T[0][2] =  1.;
          T[1][0] =  0.; T[1][1] =  1.; T[1][2] =  0.;
        }
      else if (id_side == 5)
        {
	  //T[0][0] =  1.; T[0][1] =  0.; T[0][2] =  0.;
	  //T[1][0] =  0.; T[1][1] = -1.; T[1][2] =  0.;
          T[0][0] =  0.; T[0][1] =  1.; T[0][2] =  0.;
          T[1][0] =  1.; T[1][1] =  0.; T[1][2] =  0.;
        }
      else if (id_side == 6)
        {
          T[0][0] =  1.; T[0][1] =  0.; T[0][2] =  0.;
          T[1][0] =  0.; T[1][1] =  1.; T[1][2] =  0.;
        }
      else
        {
          EH(-1, "Incorrect side for HEXAHEDRAL");
        }
      break;

    case TETRAHEDRON:
      if ( id_side == 1 )
        {
          T[0][0] =  1.; T[0][1] =  0.; T[0][2] =  0.;
          T[1][0] = -1.; T[1][1] =  0.; T[1][2] =  1.;
        }
      else if (id_side == 2)
        {
          T[0][0] = -1.; T[0][1] =  1.; T[0][2] =  0.;
          T[1][0] =  0.; T[1][1] = -1.; T[1][2] =  1.;
        }
      else if (id_side == 3)
        {
          T[0][0] =  0.; T[0][1] =  0.; T[0][2] =  1.;
          T[1][0] =  0.; T[1][1] =  1.; T[1][2] = -1.;
        }
      else if (id_side == 4)
        {
          T[0][0] =  0.; T[0][1] =  1.; T[0][2] =  0.;
          T[1][0] =  1.; T[1][1] = -1.; T[1][2] =  0.;
        }
      else
        {
          EH(-1, "Incorrect side for TETRAHEDRON");
        }
      break;

    default:
      EH(-1, "Invalid shape");
    }
    break;
  }

  /* transform from element to physical coords */
  for (p = 0; p < ielem_surf_dim; p++)
    {
      for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++)
        {
          t[p][a] = 0.;
          for (b = 0; b < ei[pg->imtrx]->ielem_dim; b++)
            {
              t[p][a] += map_bf->J[b][a] * T[p][b] * fv->h[b];
            }
        }
    }

  if (af->Assemble_Jacobian && DeformingMesh)
    {
      for (p = 0; p < ielem_surf_dim; p++)
        {
          for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++)
            {
              for (i = 0; i < num_nodes_on_side; i++)
                {
                  id   = local_elem_node_id[i];
                  inode = Proc_Elem_Connect[iconnect_ptr + id];
                  ldof = ei[pg->imtrx]->ln_to_dof[ShapeVar][id];
                  if (Dolphin[pg->imtrx][inode][MESH_DISPLACEMENT1] > 0)
                    {
                      phi_i = map_bf->phi[ldof];

                      for (q = 0; q < ei[pg->imtrx]->ielem_dim; q++)
                        {
                          dt_x[p][a][q][ldof] = 0.0;
                          for (b = 0; b < ei[pg->imtrx]->ielem_dim; b++)
                            {
                              dt_x[p][a][q][ldof] += T[p][b] * (map_bf->dJ[b][a][q][ldof] * fv->h[b] +
                                                                map_bf->J[b][a] * phi_i * fv->hq[b][q]);
                            }
                        }
                    }
                }
            }
        }
    }

  if (ielem_surf_dim == 1)
    {
      /* calculate surface determinant using the coordinate scale factors
       * for orthogonal curvilinear coordinates */
      det_h01 = sqrt(t[0][0]*t[0][0] + t[0][1]*t[0][1]);

      /*
       * If the shells are not aligned in the x-y plane, det_h01 is zero
       * When we set up the surface normal/determinant,  we assume the 
       * shells are aligned with the x-y plane.  This avoids dividing by 
       * zero and giving junk for the surface normal/determinant
       * DSH 03/24/2016
       */
      if(det_h01==0)
	{
	  EH(-1, "The shell elements need to be aligned in the X-Y plane for this problem to work");
	}

      r_det_h01 = 1. / det_h01;
      
      fv->sdet = fv->h[2] * det_h01;
      r_det = 1. / fv->sdet;

      /* calculate surface normal using the coordinate scale factors
       * for orthogonal curvilinear coordinates */
      fv->snormal[0] =  t[0][1]*r_det_h01;
      fv->snormal[1] = -t[0][0]*r_det_h01;
      fv->snormal[2] =  0.;

      /* Calculate sensitivity w.r.t. mesh, if applicable */
      if (af->Assemble_Jacobian && DeformingMesh)
        {
          for (i=0; i<num_nodes_on_side; i++)
            {
              id   = (int) local_elem_node_id[i];
              inode = Proc_Elem_Connect[iconnect_ptr + id];
              ldof = ei[pg->imtrx]->ln_to_dof[ShapeVar][id];
              if (Dolphin[pg->imtrx][inode][MESH_DISPLACEMENT1] > 0 )
                {
                  phi_i = map_bf->phi[ldof];
 
                  for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++)
                    {
                      d_det_h01_x = 0.;
                      for ( b = 0; b < ei[pg->imtrx]->ielem_dim; b++)
                        {
                          d_det_h01_x += r_det_h01 * t[0][b] * dt_x[0][b][a][ldof];
                        }

                      /* calculate sensitivity of surface determinant using the coordinate scale factors
                       * for orthogonal curvilinear coordinates */ 
                      fv->dsurfdet_dx[a][ldof] = fv->hq[2][a] * phi_i * det_h01 +
                                                 fv->h[2] * d_det_h01_x;

                      /* calculate sensitivity of surface normal using the coordinate scale factors
                       * for orthogonal curvilinear coordinates */
                      fv->dsnormal_dx[0][a][ldof] = r_det_h01 * ( dt_x[0][1][a][ldof] -
                                                                  r_det_h01 * t[0][1] * d_det_h01_x );
                      fv->dsnormal_dx[1][a][ldof] = r_det_h01 * (-dt_x[0][0][a][ldof] +
                                                                  r_det_h01 * t[0][0] * d_det_h01_x );
                    }   
                }
            }
        }
    }
  else if ( ielem_surf_dim == 2 )
    {
      double nx = t[0][1] * t[1][2] - t[0][2] * t[1][1];
      double ny = t[0][2] * t[1][0] - t[0][0] * t[1][2];
      double nz = t[0][0] * t[1][1] - t[0][1] * t[1][0];
      double d_nx_x, d_ny_x, d_nz_x;
      
      /* calculate surface determinant using the coordinate scale factors
       * for orthogonal curvilinear coordinates */
      fv->sdet = sqrt( nx * nx + ny * ny + nz * nz );
      r_det = 1. / fv->sdet;

      /* calculate surface normal using the coordinate scale factors
       * for orthogonal curvilinear coordinates */
      fv->snormal[0] = r_det * nx;
      fv->snormal[1] = r_det * ny;
      fv->snormal[2] = r_det * nz;

      /* Calculate sensitivity w.r.t. mesh, if applicable */
      if (af->Assemble_Jacobian && DeformingMesh)
        {
          for (i=0; i<num_nodes_on_side; i++)
            {
              id   = (int) local_elem_node_id[i];
              inode = Proc_Elem_Connect[iconnect_ptr + id];
              ldof = ei[pd->mi[MESH_DISPLACEMENT1]]->ln_to_dof[ShapeVar][id];
              if (Dolphin[pg->imtrx][inode][MESH_DISPLACEMENT1] > 0 )
                {
                  phi_i = map_bf->phi[ldof];

                  for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++)
                    {
                      d_nx_x = t[0][1] * dt_x[1][2][a][ldof] + dt_x[0][1][a][ldof] * t[1][2] -
                               t[0][2] * dt_x[1][1][a][ldof] - dt_x[0][2][a][ldof] * t[1][1];
                      d_ny_x = t[0][2] * dt_x[1][0][a][ldof] + dt_x[0][2][a][ldof] * t[1][0] -
                               t[0][0] * dt_x[1][2][a][ldof] - dt_x[0][0][a][ldof] * t[1][2];
                      d_nz_x = t[0][0] * dt_x[1][1][a][ldof] + dt_x[0][0][a][ldof] * t[1][1] -
                               t[0][1] * dt_x[1][0][a][ldof] - dt_x[0][1][a][ldof] * t[1][0];
                                      
                      /* calculate sensitivity of surface determinant using the coordinate scale factors
                       * for orthogonal curvilinear coordinates */
                      fv->dsurfdet_dx[a][ldof] = r_det * ( nx * d_nx_x + ny * d_ny_x + nz * d_nz_x );

                      /* calculate sensitivity of surface normal using the coordinate scale factors
                       * for orthogonal curvilinear coordinates */
                      fv->dsnormal_dx[0][a][ldof] =  r_det * ( d_nx_x -
                                                               r_det * nx * fv->dsurfdet_dx[a][ldof] );
                      fv->dsnormal_dx[1][a][ldof] =  r_det * ( d_ny_x -
                                                               r_det * ny * fv->dsurfdet_dx[a][ldof] );
                      fv->dsnormal_dx[2][a][ldof] =  r_det * ( d_nz_x -
                                                               r_det * nz * fv->dsurfdet_dx[a][ldof] );
                    }
                }
            }
        }
    } 
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void
edge_determinant_and_vectors(
 const int ielem,		/* current element number                    */
 const int iconnect_ptr,	/* Pointer into the beginning of the
				 * connectivity list for the current
				 * element			             */
 const int nodes_per_elem,	/* number of nodes in the element            */
 const int ielem_surf_dim,	/* the physical dimension of the
				 * surface of the element (0, 1, 2)          */
 const int id_side,		/* ID of the side of the element             */
 const int num_nodes_on_side,	/* number of nodes on the side of the
				 * element			             */
 const int local_elem_node_id[], /* vector of the local element node 
				  * numbers on the side of the element       */
 const int id_edge,		/* ID of the edge of the element             */
 const int num_nodes_on_edge,	/* number of nodes on the side of the
				 * element			             */
 const int edge_elem_node_id[],	/* vector of the local element node 
				 * numbers on the side of the element        */
 const int param_dir)		/* local coordinate which parameterizes edge */

/* 
 *        Function which calculates the inverse of the Jacobian matrix and
 *	the determinant of the Jacobian.
 *
 *       Author:          Harry Moffat (1421)
 *
 *	Insulate for more robust behavior for fixed grid problems - 950306pas
 *
 *      Incorporate coordinate scale factors into the surface determinant 
 *         and combine with normal calculation
 *               Richard Cairncross 9/24/96
 */

{
/* TAB certifies that this function conforms to the exo/patran side numbering convention 11/10/98. */

  int 		i, id, inode, i_basis, p, dim;
  int		ShapeVar, ldof;
  int		DeformingMesh;
  double        dxdalpha, dydalpha, dzdalpha;
  double        det;
  double        phi_i, sign;
  int siz;

  /* first step - make sure surface determinant and normal are calculated 
  *
  *  The surface determinant and normal are calculated with respect to the primary surface
  *  for this edge (i.e. the first surface listed by the user in defining the edge)
  */
 surface_determinant_and_normal(ielem, iconnect_ptr, nodes_per_elem, ielem_surf_dim,  
	       id_side, num_nodes_on_side, local_elem_node_id);

 dim = ielem_surf_dim + 1;
  DeformingMesh = pd->e[pg->imtrx][R_MESH1];
  ShapeVar = pd->ShapeVar;

  /* initialize variables */
  siz = MAX_PDIM*MDE*sizeof(double);
  memset( fv->dedgedet_dx,0,siz);
  siz = MAX_PDIM*MDE*MAX_PDIM*sizeof(double);
  memset( fv->dstangent_dx[0],0,siz);
  memset( fv->dstangent_dx[1],0,siz);
  det = fv->edge_det = 0.0;
  fv->stangent[0][0] = 0.;
  fv->stangent[0][1] = 0.;
  fv->stangent[0][2] = 0.;
  fv->stangent[1][0] = 0.;
  fv->stangent[1][1] = 0.;
  fv->stangent[1][2] = 0.;

  switch (ielem_surf_dim) {

  case 0:
    EH(-1, "Edges undefined in 1D");
    return;

    /* 2D Problem . . . */
  case 1:
    if (pd->CoordinateSystem == CARTESIAN && DEFAULT_CARTESIAN) {
      det = fv->edge_det = 1.;
    } else {
      det = fv->edge_det = fv->h[2];
      if (af->Assemble_Jacobian && DeformingMesh) {
	for (i=0; i<num_nodes_on_edge; i++)
	  {
	    id   = (int) edge_elem_node_id[i];
	    inode = Proc_Elem_Connect[iconnect_ptr + id];
	    ldof = ei[pg->imtrx]->ln_to_dof[ShapeVar][id];
	    if (Dolphin[pg->imtrx][inode][MESH_DISPLACEMENT1] > 0 )
	      {		  
		phi_i = bf[ShapeVar]->phi[ldof];
		    fv->dedgedet_dx[0][ldof] = fv->hq[2][0] * phi_i;
		    fv->dedgedet_dx[1][ldof] = fv->hq[2][1] * phi_i;
	      }
	  }
      }
    }

    /* find sign of edge tangent so binormal will be outward pointing */
    sign = 0.;
    switch (id_side) {
    case(1):
      switch (id_edge) {
      case(2):
	sign = 1.;
	break;
      case(1):
	sign = -1.;
	break;
      default:
	EH(-1,"Side not connected to edge");
      }
      break;
    case(2):
      switch (id_edge) {
      case(3):
	sign = 1.;
	break;
      case(2):
	sign = -1.;
	break;
      default:
	EH(-1,"Side not connected to edge");
      }
      break;
    case(3):
      switch (id_edge) {
      case(4):
	sign = 1.;
	break;
      case(3):
	sign = -1.;
	break;
      default:
	EH(-1,"Side not connected to edge");
      }
      break;
    case(4):
      switch (id_edge) {
      case(1):
	sign = 1.;
	break;
      case(4):
	sign = -1.;
      break;
      default:
	EH(-1,"Side not connected to edge");
      }
      break;
      
    default:
      EH(-1,"Edge not found");
    }
    
    fv->stangent[1][0]= 0.;
    fv->stangent[1][1]= 0.;
    fv->stangent[1][2]= sign;
    /* find binormal */
    fv->stangent[0][0]= fv->snormal[1] * sign;
    fv->stangent[0][1]= - fv->snormal[0] * sign;
    fv->stangent[0][2]= 0.;

    if (af->Assemble_Jacobian && DeformingMesh) {
      for (i=0; i<num_nodes_on_edge; i++)
	{
	  id   = (int) edge_elem_node_id[i];
	  inode = Proc_Elem_Connect[iconnect_ptr + id];
	  ldof = ei[pg->imtrx]->ln_to_dof[ShapeVar][id];
	  for (p=0; p<dim; p++) {
	      if (Dolphin[pg->imtrx][inode][MESH_DISPLACEMENT1 + p] > 0 )
		{	
		  fv->dstangent_dx[0][0][p][ldof] = fv->dsnormal_dx[1][p][ldof] * sign;
		  fv->dstangent_dx[0][1][p][ldof] = - fv->dsnormal_dx[0][p][ldof] * sign;
		}
	  }
	}
    }	  

    return;

    /* 3D Problem . . .  */
  case 2:
    dxdalpha = 0.;
    dydalpha = 0.;
    dzdalpha = 0.;

    /* find parametric basis along edge and sign of direction */
    i_basis = param_dir;

    for (i = 0; i < num_nodes_on_edge; i++)
      {
	id    = (int) edge_elem_node_id[i];
	inode = Proc_Elem_Connect[iconnect_ptr + id];
	ldof  = ei[pg->imtrx]->ln_to_dof[ShapeVar][id];
	if ( DeformingMesh )
	  {
	    if (Dolphin[pg->imtrx][inode][MESH_DISPLACEMENT1] > 0)
	      {
		dxdalpha += bf[ShapeVar]->dphidxi[ldof][i_basis]*(Coor[0][inode] + *esp->d[0][ldof]);
		dydalpha += bf[ShapeVar]->dphidxi[ldof][i_basis]*(Coor[1][inode] + *esp->d[1][ldof]);
		dzdalpha += bf[ShapeVar]->dphidxi[ldof][i_basis]*(Coor[2][inode] + *esp->d[2][ldof]);
	      }
	  } /*end of Baby_Dolphin[pg->imtrx] */
	else
	  {
	    dxdalpha += bf[ShapeVar]->dphidxi[ldof][i_basis]*Coor[0][inode];
	    dydalpha += bf[ShapeVar]->dphidxi[ldof][i_basis]*Coor[1][inode];
	    dzdalpha += bf[ShapeVar]->dphidxi[ldof][i_basis]*Coor[2][inode];
	  }
      }


    if (pd->CoordinateSystem == CARTESIAN && DEFAULT_CARTESIAN) {
    /* calculate surface determinant in cartesian coords */
      det = fv->edge_det = sqrt(dxdalpha*dxdalpha + dydalpha*dydalpha + dzdalpha*dzdalpha);
    } else {
      /* calculate surface determinant using the coordinate scale factors
       * for orthogonal curvilinear coordinates */
      det = fv->edge_det = sqrt(fv->h[0] * fv->h[0] * dxdalpha*dxdalpha + 
				fv->h[1] * fv->h[1] * dydalpha*dydalpha + 
				fv->h[2] * fv->h[2] * dzdalpha*dzdalpha);
    }
  
    /* find sign for tangent direction */
  /* determine sign on parametric curve to make Normal, line tangent, and 
   * binormal (outward pointing normal to egde) a right-handed basis */
  sign = 0.;
  switch (id_edge) {
  case(1):
    switch (id_side) {
    case(4):
      sign = 1.;
      break;
    case(5):
      sign = -1.;
      break;
    default:
      EH(-1,"Side not connected to edge");
    }
    break;

  case(2):
    switch (id_side) {
    case(5):
      sign = 1.;
      break;
    case(2):
      sign = -1.;
      break;
    default:
      EH(-1,"Side not connected to edge");
    }
    break;

  case(3):
    switch (id_side) {
    case(5):
      sign = 1.;
      break;
    case(1):
      sign = -1.;
      break;
    default:
      EH(-1,"Side not connected to edge");
    }
    break;

  case(4):
    switch (id_side) {
    case(3):
      sign = 1.;
      break;
    case(5):
      sign = -1.;
      break;
    default:
      EH(-1,"Side not connected to edge");
    }
    break;

  case(5):
    switch (id_side) {
    case(6):
      sign = 1.;
      break;
    case(4):
      sign = -1.;
      break;
    default:
      EH(-1,"Side not connected to edge");
    }
    break;

  case(6):
    switch (id_side) {
    case(2):
      sign = 1.;
      break;
    case(6):
      sign = -1.;
      break;
    default:
      EH(-1,"Side not connected to edge");
    }
    break;

  case(7):
    switch (id_side) {
    case(1):
      sign = 1.;
      break;
    case(6):
      sign = -1.;
      break;
    default:
      EH(-1,"Side not connected to edge");
    }
    break;

  case(8):
    switch (id_side) {
    case(6):
      sign = 1.;
      break;
    case(3):
      sign = -1.;
      break;
    default:
      EH(-1,"Side not connected to edge");
    }
    break;

  case(9):
    switch (id_side) {
    case(1):
      sign = 1.;
      break;
    case(4):
      sign = -1.;
      break;
    default:
      EH(-1,"Side not connected to edge");
    }
    break;

  case(10):
    switch (id_side) {
    case(4):
      sign = 1.;
      break;
    case(3):
      sign = -1.;
      break;
    default:
      EH(-1,"Side not connected to edge");
    }
    break;

  case(11):
    switch (id_side) {
    case(2):
      sign = 1.;
      break;
    case(1):
      sign = -1.;
      break;
    default:
      EH(-1,"Side not connected to edge");
    }
    break;

  case(12):
    switch (id_side) {
    case(3):
      sign = 1.;
      break;
    case(2):
      sign = -1.;
      break;
    default:
      EH(-1,"Side not connected to edge");
    }
    break;

    default:
      EH(-1,"Edge not found");
  }

    if (pd->CoordinateSystem == CARTESIAN && DEFAULT_CARTESIAN) {
      /* calculate tangent to surface */
      fv->stangent[1][0]=sign*dxdalpha/det; 
      fv->stangent[1][1]=sign*dydalpha/det; 
      fv->stangent[1][2]=sign*dzdalpha/det; 
    } else {
      /* calculate edge tangent using the coordinate scale factors
       * for orthogonal curvilinear coordinates 
       * NOTE that this isn't always the same as the above normal*/
      fv->stangent[1][0]=sign*dxdalpha*fv->h[0]/det;
      fv->stangent[1][1]=sign*dydalpha*fv->h[1]/det;
      fv->stangent[1][2]=sign*dzdalpha*fv->h[2]/det;
    }

    /* Calculate sensitivity w.r.t. mesh, if applicable */
    if (af->Assemble_Jacobian && DeformingMesh) {
      for (i=0; i<num_nodes_on_edge; i++)
	{
	  id   = (int) edge_elem_node_id[i];
	  inode = Proc_Elem_Connect[iconnect_ptr + id];
	  ldof = ei[pg->imtrx]->ln_to_dof[ShapeVar][id];
	      if (Dolphin[pg->imtrx][inode][MESH_DISPLACEMENT1] > 0 )
		{		  
		  phi_i = bf[ShapeVar]->phi[ldof];

		  if (pd->CoordinateSystem == CARTESIAN && DEFAULT_CARTESIAN) {
		    fv->dedgedet_dx[0][ldof] = (1./det)*(bf[ShapeVar]->dphidxi[ldof][i_basis]*dxdalpha);
		    fv->dedgedet_dx[1][ldof] = (1./det)*(bf[ShapeVar]->dphidxi[ldof][i_basis]*dydalpha);
		    fv->dedgedet_dx[2][ldof] = (1./det)*(bf[ShapeVar]->dphidxi[ldof][i_basis]*dzdalpha);
		  } else {
		    /* calculate sensitivity of surface determinant using the coordinate scale factors
		     * for orthogonal curvilinear coordinates */
		    fv->dedgedet_dx[0][ldof] = (1./det)*
		      ( fv->h[0] * fv->h[0] * bf[ShapeVar]->dphidxi[ldof][i_basis]*dxdalpha + 
			( fv->h[0] * fv->hq[0][0] * dxdalpha*dxdalpha +
			  fv->h[1] * fv->hq[1][0] * dydalpha*dydalpha + 
			  fv->h[2] * fv->hq[2][0] * dzdalpha*dzdalpha) * phi_i);
		    fv->dedgedet_dx[1][ldof] = (1./det)*
		      ( fv->h[1] * fv->h[1] * bf[ShapeVar]->dphidxi[ldof][i_basis]*dydalpha + 
			( fv->h[0] * fv->hq[0][1] * dxdalpha*dxdalpha +
			  fv->h[1] * fv->hq[1][1] * dydalpha*dydalpha + 
			  fv->h[2] * fv->hq[2][1] * dzdalpha*dzdalpha) * phi_i);
		    fv->dedgedet_dx[2][ldof] = (1./det)*
		      ( fv->h[2] * fv->h[2] * bf[ShapeVar]->dphidxi[ldof][i_basis]*dzdalpha + 
			( fv->h[0] * fv->hq[0][2] * dxdalpha*dxdalpha +
			  fv->h[1] * fv->hq[1][2] * dydalpha*dydalpha + 
			  fv->h[2] * fv->hq[2][2] * dzdalpha*dzdalpha) * phi_i);

		  }

		  if (pd->CoordinateSystem == CARTESIAN && DEFAULT_CARTESIAN) {
		      		  
		    fv->dstangent_dx[1][0][0][ldof]= sign 
		      * ( bf[ShapeVar]->dphidxi[ldof][i_basis]/det -
			  dxdalpha*fv->dedgedet_dx[0][ldof]/(det*det));
		      		  
		    fv->dstangent_dx[1][0][1][ldof]= sign 
		      * ( 0. -
			  dxdalpha*fv->dedgedet_dx[1][ldof]/(det*det));
		      		  
		    fv->dstangent_dx[1][0][2][ldof]= sign 
		      * ( 0. -
			  dxdalpha*fv->dedgedet_dx[2][ldof]/(det*det));
		      		  
		    fv->dstangent_dx[1][1][0][ldof]= sign 
		      * ( 0. -
			  dydalpha*fv->dedgedet_dx[0][ldof]/(det*det));
		      		  
		    fv->dstangent_dx[1][1][1][ldof]= sign 
		      * ( bf[ShapeVar]->dphidxi[ldof][i_basis]/det -
			  dydalpha*fv->dedgedet_dx[1][ldof]/(det*det));
		      		  
		    fv->dstangent_dx[1][1][2][ldof]= sign 
		      * ( 0. -
			  dydalpha*fv->dedgedet_dx[2][ldof]/(det*det));
		      		  
		    fv->dstangent_dx[1][2][0][ldof]= sign 
		      * ( 0. -
			  dzdalpha*fv->dedgedet_dx[0][ldof]/(det*det));
		      		  
		    fv->dstangent_dx[1][2][1][ldof]= sign 
		      * ( 0. -
			  dzdalpha*fv->dedgedet_dx[1][ldof]/(det*det));
		      		  
		    fv->dstangent_dx[1][2][2][ldof]= sign 
		      * ( bf[ShapeVar]->dphidxi[ldof][i_basis]/det -
			  dzdalpha*fv->dedgedet_dx[2][ldof]/(det*det));
		    
		  } else { 
		    /* calculate sensitivity of edge tangent using the coordinate scale factors
		     * for orthogonal curvilinear coordinates 
		     * NOTE that this isn't always the same as the above normal*/
		      		  
		    fv->dstangent_dx[1][0][0][ldof]= sign 
		      * ( bf[ShapeVar]->dphidxi[ldof][i_basis]*fv->h[0]/det -
			  dxdalpha*fv->dedgedet_dx[0][ldof]*fv->h[0]/(det*det) +
			  dxdalpha*fv->hq[0][0]*phi_i/det);
		      		  
		    fv->dstangent_dx[1][0][1][ldof]= sign 
		      * ( 0. -
			  dxdalpha*fv->dedgedet_dx[1][ldof]*fv->h[0]/(det*det) +
			  dxdalpha*fv->hq[0][1]*phi_i/det);
		      		  
		    fv->dstangent_dx[1][0][2][ldof]= sign 
		      * ( 0. -
			  dxdalpha*fv->dedgedet_dx[2][ldof]*fv->h[0]/(det*det) +
			  dxdalpha*fv->hq[0][2]*phi_i/det);
		      		  
		    fv->dstangent_dx[1][1][0][ldof]= sign 
		      * ( 0. -
			  dydalpha*fv->dedgedet_dx[0][ldof]*fv->h[1]/(det*det) +
			  dydalpha*fv->hq[1][0]*phi_i/det);
		      		  
		    fv->dstangent_dx[1][1][1][ldof]= sign 
		      * ( bf[ShapeVar]->dphidxi[ldof][i_basis]*fv->h[1]/det -
			  dydalpha*fv->dedgedet_dx[1][ldof]*fv->h[1]/(det*det) +
			  dydalpha*fv->hq[1][1]*phi_i/det);
		      		  
		    fv->dstangent_dx[1][1][2][ldof]= sign 
		      * ( 0. -
			  dydalpha*fv->dedgedet_dx[2][ldof]*fv->h[1]/(det*det) +
			  dydalpha*fv->hq[1][2]*phi_i/det);
		      		  
		    fv->dstangent_dx[1][2][0][ldof]= sign 
		      * ( 0. -
			  dzdalpha*fv->dedgedet_dx[0][ldof]*fv->h[2]/(det*det) +
			  dzdalpha*fv->hq[2][0]*phi_i/det);
		      		  
		    fv->dstangent_dx[1][2][1][ldof]= sign 
		      * ( 0. -
			  dzdalpha*fv->dedgedet_dx[1][ldof]*fv->h[2]/(det*det) +
			  dzdalpha*fv->hq[2][1]*phi_i/det);
		      		  
		    fv->dstangent_dx[1][2][2][ldof]= sign 
		      * ( bf[ShapeVar]->dphidxi[ldof][i_basis]*fv->h[2]/det -
			  dzdalpha*fv->dedgedet_dx[2][ldof]*fv->h[2]/(det*det) +
			  dzdalpha*fv->hq[2][2]*phi_i/det);
		  }
		}/* end of Baby_Dolphin[pg->imtrx] */
	}
    } /* end of Jacobian calculations */

    /* find binormal */
    fv->stangent[0][0]= fv->snormal[1] * fv->stangent[1][2] 
                      - fv->snormal[2] * fv->stangent[1][1];
    fv->stangent[0][1]= fv->snormal[2] * fv->stangent[1][0] 
                      - fv->snormal[0] * fv->stangent[1][2];
    fv->stangent[0][2]= fv->snormal[0] * fv->stangent[1][1] 
                      - fv->snormal[1] * fv->stangent[1][0];

    /* calculate sensitivities */
    if (af->Assemble_Jacobian && DeformingMesh) {
      for (i=0; i<num_nodes_on_edge; i++)
	{
	  id   = (int) edge_elem_node_id[i];
	  inode = Proc_Elem_Connect[iconnect_ptr + id];
	  ldof = ei[pg->imtrx]->ln_to_dof[ShapeVar][id];
	  for (p=0; p<dim; p++) {
	  if (Dolphin[pg->imtrx][inode][MESH_DISPLACEMENT1 + p] > 0 )
	    {	
	      fv->dstangent_dx[0][0][p][ldof] = 
		fv->dsnormal_dx[1][p][ldof] * fv->stangent[1][2] 
		- fv->dsnormal_dx[2][p][ldof] * fv->stangent[1][1] 
		+ fv->snormal[1] * fv->dstangent_dx[1][2][p][ldof] 
		- fv->snormal[2] * fv->dstangent_dx[1][1][p][ldof];
	      fv->dstangent_dx[0][1][p][ldof] = 
		fv->dsnormal_dx[2][p][ldof] * fv->stangent[1][0] 
		- fv->dsnormal_dx[0][p][ldof] * fv->stangent[1][2]
		+ fv->snormal[2] * fv->dstangent_dx[1][0][p][ldof]
		- fv->snormal[0] * fv->dstangent_dx[1][2][p][ldof];
	      fv->dstangent_dx[0][2][p][ldof] = 
		fv->dsnormal_dx[0][p][ldof] * fv->stangent[1][1] 
		- fv->dsnormal_dx[1][p][ldof] * fv->stangent[1][0]
		+ fv->snormal[0] * fv->dstangent_dx[1][1][p][ldof] 
		- fv->snormal[1] * fv->dstangent_dx[1][0][p][ldof];
	    }
	  }
	}
    }	  
    return;
    
  default:
    EH(-1, "Bad dimension to calc_edge_det.");
    break;
  }

} /* edge_determinant_and_vectors */


/*
 * Routine determines normal to contact line defined by tangent vector
 * and substrate normal
 */

void
calc_CL_normal ( double snormal[DIM],
		 double dsnormal_dx[DIM][DIM][MDE],
		 double fsnormal[DIM],
		 double dfsnormal_dx[DIM][DIM][MDE],
		 double tangent[DIM],
		 double dtangent_dx[DIM][DIM][MDE],		 
		 int elem,
		 int edge_elem_node_id[],
		 int elem_dim,
		 int num_edge_nodes,
		 double clnormal[DIM],
		 double dclnormal_dx[DIM][DIM][MDE],
		 const Exo_DB *exo)
{
  double dot;
  double sign=1.0;
  int dim = elem_dim;
  char err_msg[MAX_CHAR_IN_INPUT];
  if( dim < 3 ) EH(-1,"VAR_CA_EDGE invalid for 2D simulations ");

  clnormal[0]= snormal[1] * tangent[2] - snormal[2] * tangent[1];
  clnormal[1]= snormal[2] * tangent[0] - snormal[0] * tangent[2];
  clnormal[2]= snormal[0] * tangent[1] - snormal[1] * tangent[0];

  memset((void *)dclnormal_dx, 0, sizeof(double)*DIM*DIM*MDE);

  /*
   * Check sense of clnormal by forming dot product with fsnormal
   */
   
  if ( (dot = clnormal[0]*fsnormal[0] +
	      clnormal[1]*fsnormal[1] +
	      clnormal[2]*fsnormal[2] ) < 0. )
    {
      /* reverse clnormal */
      sign = -1.0;
      clnormal[0] *= sign;
      clnormal[1] *= sign;
      clnormal[2] *= sign;
    }
  else if ( fabs( dot ) < 1.e-15 )
    {
      sprintf(err_msg,"Free surface normal and contact line normal orthogonal at elem %d \n",elem);
      EH(-1,err_msg);
    }
  
  
  /* calculate sensitivities */
  if (af->Assemble_Jacobian && pd->e[pg->imtrx][R_MESH1] ) 
    {
      int i, id, Inode,ldof,p;
      int iconnect_ptr = exo->elem_ptr[elem];

      for (i=0; i<num_edge_nodes; i++)
	{
	  id   =  edge_elem_node_id[i];
	  Inode = exo->node_list[iconnect_ptr + id];
	  ldof = ei[pg->imtrx]->ln_to_dof[pd->ShapeVar][id];
	    
	  for (p=0; p<dim; p++) 
	    {
	      if (Dolphin[pg->imtrx][Inode][MESH_DISPLACEMENT1 + p] > 0 )
		{	
		  dclnormal_dx[0][p][ldof] = sign * 
		    (      dsnormal_dx[1][p][ldof] * tangent[2] 
			 - dsnormal_dx[2][p][ldof] * tangent[1] 
			 + snormal[1] * dtangent_dx[2][p][ldof] 
			 - snormal[2] * dtangent_dx[1][p][ldof] );

		  dclnormal_dx[1][p][ldof] = sign * 
		    (      dsnormal_dx[2][p][ldof] * tangent[0] 
			 - dsnormal_dx[0][p][ldof] * tangent[2] 
			 + snormal[2] * dtangent_dx[0][p][ldof] 
			 - snormal[0] * dtangent_dx[2][p][ldof] );
		    

		  dclnormal_dx[2][p][ldof] = sign * 
		      (    dsnormal_dx[0][p][ldof] * tangent[1] 
			 - dsnormal_dx[1][p][ldof] * tangent[0] 
			 + snormal[0] * dtangent_dx[1][p][ldof] 
			 - snormal[1] * dtangent_dx[0][p][ldof] );
		}
	    }
	}
    }
}


/*
 * Routine to compute the necessary terms for SU/PG
 */
void
get_supg_stuff(dbl *supg_term,
	       dbl vcent[DIM],
	       dbl d_vcent_du[DIM][MDE][DIM],
	       dbl d_supg_term_du[MDE][DIM],
	       dbl d_supg_term_dx[MDE][DIM],
	       const int DeformingMesh)
     
    /*********************************************************************
     * get_supg_stuff()
     *
     * This routine calculates terms (and derivatives of said terms)
     * for the discontinuous contribution to the SU/PG FEM function.
     *
     * The heuristics of this function are only designed for quad/hex
     * element types. If other elements are used, a new algorithm must
     * be developed.
     *
     * Output
     * ------
     *   supg_term   = This is the discontinuous contribution to the 
     *                 SU/PG weight function with the exception of the 
     *                 u.grad(w) term:
     *
     *                 w + p = w + SUM { supg_term * u_j * (dw/dx_j) }
     *                              j
     *
     *   vcent[i]    = The i-th component of the centroid velocity.
     *
     *   d_vcent_du[i][j][k] = The derivative of vcent[i] w.r.t. the k-th 
     *                 component of the velocity at node j.
     *
     *   d_supg_term_du[i][j] = The derivative of supg_term w.r.t. the
     *                 j-th component of the velocity at node i.
     *
     *   d_supg_term_dx[i][j] = The derivative of supg_term w.r.t. the
     *                 j-th component of the position vector of node i.
     *
     *
     * References
     * ----------
     * (1) A. N. Brooks and T. J. R. Hughes, "Streamline upwint/Petrov-
     *     Galerkin formulations for convection dominated flows with 
     *     particular emphasis on the incompressible Navier-Stokes equations.", 
     *     Comp. Meth. Appl. Mech., v32, 199--259 (1982).
     *
     *********************************************************************/
{
  int i, j, k, p, q, idof, dim, index, var, velodim, dofs = 0;
  dbl x_coeff[DIM][MDE],v_coeff[MDE];
  dbl xnode[DIM][MDE],d_vect[DIM][DIM],v[DIM][MDE];
  dbl ex,emx,foo,vmag;
  dbl hlen[DIM],alpha[DIM],cotha[DIM],d_alpha_dx[DIM][MDE][DIM];
  dbl d_alpha_du[DIM][MDE][DIM];
  dbl small_number;
  dbl vmag2;
  
  dim =  pd->Num_Dim;

  /**********************************************************************
   *                                                                    *
   *        EXTRACT THE LOCAL COORDINATE AND VELCOITY COMPONENTS        *
   *                                                                    *
   **********************************************************************/
  /*
   * Find the local coodinates of the nodes in the current element
   * -> Can this be pushed to a more generic place ?
   */
  j = Proc_Connect_Ptr[ei[pg->imtrx]->ielem];
  for (p = 0; p < dim; p++)
    {
      if (DeformingMesh)
	{
	  var = MESH_DISPLACEMENT1 + p;
	  for (i = 0; i < ei[pg->imtrx]->num_local_nodes; i++)
	    {
	      idof = ei[pg->imtrx]->ln_to_dof[var][i];
	      index = Proc_Elem_Connect[j + i];
	      xnode[p][i]= Coor[p][index] + *esp->d[p][idof];
	    }
	}
      else
	{
	  for (i = 0; i < ei[pg->imtrx]->num_local_nodes; i++)
	    {
	      index = Proc_Elem_Connect[j  + i];
	      xnode[p][i]= Coor[p][index];
	    }
	}
    }
  
  /*
   * Velocity code pulled from load_fv() and then severly flogged by PKN:
   */
  
  velodim = dim;		/* Later, this might include v_theta... */
  
  /*
   * MMH: Change here in case the 3rd velocity component is non-zero
   * (non LSA situation).
   */
  if(pd->CoordinateSystem == SWIRLING ||
     pd->CoordinateSystem == PROJECTED_CARTESIAN ||
     pd->CoordinateSystem == CARTESIAN_2pt5D)
    velodim = dim + 1; /* Later is Now!  Woo!!! */
  
  for ( p=0; p<velodim; p++)
    {
      j = VELOCITY1 + p;
      if ( pd->v[pg->imtrx][j] )
	{
	  dofs = ei[pg->imtrx]->dof[j]; 
	  for(i=0; i<dofs; i++) v[p][i] = *esp->v[p][i];
	}
      else
	{
	  /*
	   * If there isn't a velocity unknown, we shouldn't 
	   * be using SUPG!
	   */
	  printf("\nget_supg_stuff(): Hmmm, are we really supposed to be here?\n");
	  for(i=0; i<MDE; i++) v[p][i] = 0.;
	}
    }
  
  /**********************************************************************
   *                                                                    *
   *       BUILD ARRAYS OF COEFFICIENTS -- MAKES JACOBIANS EASIER       *
   *                                                                    *
   **********************************************************************/

  /*
   * Set the velocity coefficients.
   * For now, use the average velocity of all nodes.
   */
  memset(v_coeff, 0, sizeof(double)*MDE);
  if ( pd->v[pg->imtrx][VELOCITY1] )
    {
      foo = 1.0 / ei[pg->imtrx]->dof[VELOCITY1];
      for ( i=0; i < dofs; i++ )
	{
	  v_coeff[i] = foo;
	}
    }
  
  /*
   * Set the coordinate-direction-vector coefficients
   */
  memset(x_coeff, 0, sizeof(double)*MDE*DIM);
  if ( dim == 2 )
    {
      /* d[0] = 0.5*(x_1+x_2) - 0.5*(x_0+x_3) */
      x_coeff[0][0] = -0.5;
      x_coeff[0][1] = +0.5;
      x_coeff[0][2] = +0.5;
      x_coeff[0][3] = -0.5;
      /* d[1] = 0.5*(x_0+x_1) - 0.5*(x_2+x_3) */
      x_coeff[1][0] = +0.5;
      x_coeff[1][1] = +0.5;
      x_coeff[1][2] = -0.5;
      x_coeff[1][3] = -0.5;
    }
  else if ( dim == 3 )
    {
      /* d[0] = 0.25*(x_1+x_2+x_5+x_6) - 0.25*(x_0+x_3+x_4+x_7) */
      x_coeff[0][0] = -0.25;
      x_coeff[0][1] = +0.25;
      x_coeff[0][2] = +0.25;
      x_coeff[0][3] = -0.25;
      x_coeff[0][4] = -0.25;
      x_coeff[0][5] = +0.25;
      x_coeff[0][6] = +0.25;
      x_coeff[0][7] = -0.25;
      /* d[1] = 0.25*(x_2+x_3+x_6+x_7) - 0.25*(x_0+x_1+x_4+x_5) */
      x_coeff[1][0] = -0.25;
      x_coeff[1][1] = -0.25;
      x_coeff[1][2] = +0.25;
      x_coeff[1][3] = +0.25;
      x_coeff[1][4] = -0.25;
      x_coeff[1][5] = -0.25;
      x_coeff[1][6] = +0.25;
      x_coeff[1][7] = +0.25;
      /* d[2] = 0.25*(x_4+x_5+x_6+x_7) - 0.25*(x_0+x_1+x_2+x_3) */
      x_coeff[2][0] = -0.25;
      x_coeff[2][1] = -0.25;
      x_coeff[2][2] = -0.25;
      x_coeff[2][3] = -0.25;
      x_coeff[2][4] = +0.25;
      x_coeff[2][5] = +0.25;
      x_coeff[2][6] = +0.25;
      x_coeff[2][7] = +0.25;
    }
  else
    {
      EH(-99,"get_supg_stuff() is confused.\n");
    }
   
  /**********************************************************************
   *                                                                    *
   *     BUILD THE NECESSARY TERMS FOR THE SU/PG WEIGHT FUNCTIONS       *
   *                                                                    *
   **********************************************************************/

  /*
   * Form position vectors for the mid points of edges (for 2D)
   * or faces (for 3D).  Only valid for Quads or Hexes.
   *
   */

  /* Centroid velocity */
  memset(vcent,0,sizeof(double)*DIM);
  vmag = 0.;
  for ( p=0; p < dim; p++ )
    {
      j = VELOCITY1 + p;
      if ( pd->v[pg->imtrx][j] )
	{
	  for(i=0; i<ei[pg->imtrx]->dof[j]; i++)
	    {
	      vcent[p] = v_coeff[i] * v[p][i];
	    }
	}
      vmag += vcent[p]*vcent[p];
    }
  vmag = sqrt(vmag);

  /* Un-normalized directional vectors */
  /* p=0 -> s-coordinate; p=1 -> t-coordinate; etc */
  /* i=0 -> 1st vector component; s=1 -> 2nd vector component; etc. */
  memset(d_vect,0,sizeof(double)*DIM*DIM);
  for ( p=0; p < dim; p++ )
    { 
    for( i=0; i < dim; i++ )
      {
	for ( j=0; j < ei[pg->imtrx]->num_local_nodes; j++ )
	  {
	    d_vect[p][i] += x_coeff[p][j] * xnode[i][j];
	  }
      }
    }

  /* direction lengths */
  memset(hlen,0,sizeof(double)*DIM);
  for ( p=0; p < dim; p++ )
    {
      for ( i=0; i < dim; i++ )
	{
	  hlen[p] += d_vect[p][i]*d_vect[p][i];
	}
      hlen[p] = sqrt(hlen[p]);
    }

  /*
   * The alpha's:
   *
   * Take "k" == 1 (see Ref. 1)
   *
   * If alpha[p] is really small, cotha[p] -> oo.
   * But, the final supg_term has alpha*coth(alpha) -> 0
   * so just set cotha[p]=0 if alpha is "small enough".
   *
   * I'm not sure if this is the best choice for small_number.
   */
  small_number = 1.0e-6;

  memset(alpha, 0, sizeof(double)*DIM);
  memset(cotha, 0, sizeof(double)*DIM);
  for ( p=0; p < dim; p++ )
    {
      for ( i=0; i < dim; i++ )
	{
	  alpha[p] += d_vect[p][i]*vcent[i];
	}
      alpha[p] *= 0.5;
      if ( fabs( alpha[p] ) > small_number )
	{
	  ex	   = exp(alpha[p]);
	  emx	   = exp(-alpha[p]);
	  cotha[p] = (ex+emx) / (ex-emx); /* coth(alpha[p]) */
	}
    }

  /* SU/PG diffusivity */
  *supg_term = 0.;
  vmag2     = vmag*vmag;
  if ( vmag2 < small_number )
    {
      *supg_term = 0.;
    }
  else
    {
      for ( p=0; p < dim; p++ )
	{
	  *supg_term += alpha[p]*cotha[p] - 1.0;
	}
      *supg_term *= 1.0 / vmag2;
    }
      
  /**********************************************************************
   *                                                                    *
   *         BUILD SOME ARRAYS FOR EASY JACOBIAN CALCULATIONS           *
   *                                                                    *
   **********************************************************************/

  /*
   * d (alpha_i) / d x_(k)q
   */
  /* i: local coordinate direction: i=0: alpha_s; i=1: alpha_t; i=2: alpha_s */
  memset(d_alpha_dx, 0, sizeof(double)*DIM*MDE*DIM);
  if ( DeformingMesh )
    {
      for ( i=0; i < dim; i++ )
	{                  
	  /* k: x_k is the position vector of the kth node */
	  for ( k=0; k < ei[pg->imtrx]->num_local_nodes; k++ )
	    {  
	      /* q: x_(k)q is the q-th component of x_k */
	      for ( q=0; q < dim; q++ )
		{
		  d_alpha_dx[i][k][q] = 0.5 * vcent[q] * x_coeff[i][k];
		}
	    }
	}
    }

  /*
   * d (alpha_i) / d v_(k)q
   */
  /* i: local coordinate direction: i=0: alpha_s; i=1: alpha_t; i=2: alpha_s */
  memset(d_alpha_du, 0, sizeof(double)*DIM*MDE*DIM);
  for ( i=0; i < dim; i++ )
    {                  
      /* k: v_k is the velocity vector at the kth node */
      for ( k=0; k < ei[pg->imtrx]->num_local_nodes; k++ )
	{  
	  /* q: v_(k)q is the q-th component of v_k */
	  for ( q=0; q < dim; q++ )
	    {
	      d_alpha_du[i][k][q] = 0.5 * d_vect[i][q] * v_coeff[k];
	    }
	}
    }
  
  /*
   * d vmag / d v_(k)q = 2 * u_q * v_coeff[k]
   */

  /*
   * d u_i / d v_(k)q = b_k delta_iq
   */
  memset(d_vcent_du,0,sizeof(double)*DIM*MDE*DIM);
  for ( i=0; i < dim; i++ )
    {
      for ( k=0; k < ei[pg->imtrx]->num_local_nodes; k++ )
	{
	  d_vcent_du[i][k][i] = v_coeff[k];
	}
    }

  /*
   * d_supg_term_du[k][q]
   */
  memset(d_supg_term_du,0,sizeof(double)*MDE*DIM);
  if ( vmag2 > small_number )
    {
      foo = 1.0 / vmag2;
      for ( k=0; k < ei[pg->imtrx]->num_local_nodes; k++ )
	{
	  for ( q=0; q < dim; q++ )
	    {
	      for ( i=0; i < dim; i++ )
		{
		  d_supg_term_du[k][q] += foo * (cotha[i] + alpha[i]*(1-cotha[i]*cotha[i]))
		    * d_alpha_du[i][k][q];
		  d_supg_term_du[k][q] += foo*foo * 2.0 * v_coeff[k] * vcent[q]
		    * (alpha[i]*cotha[i] - 1.0);
		}
	    }
	}
    }
  
  /*
   * d_supg_term_dx[k][q]
   */
  memset(d_supg_term_dx,0,sizeof(double)*MDE*DIM);
  if ( DeformingMesh && vmag2 > small_number )
    {
      foo = 1.0 / (vmag * vmag);
      for ( k=0; k < ei[pg->imtrx]->num_local_nodes; k++ )
	{
	  for ( q=0; q < dim; q++ )
	    {
	      for ( i=0; i < dim; i++ )
		{
		  d_supg_term_dx[k][q] += (cotha[i] + alpha[i]*(1. - cotha[i]*cotha[i]) )
		    * d_alpha_dx[i][k][q];
		}
	      d_supg_term_dx[k][q] *= foo;
	    }
	}
    }
  
  return;
}

/*****************************************************************************/
/*			END of mm_fill_aux.c				     */
/*****************************************************************************/
