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
 
#ifdef USE_RCSID
static char rcsid[] = "$Id: user_post.c,v 5.1 2007-09-18 18:53:49 prschun Exp $";
#endif

/* Standard include files */

#include <stdio.h>

/* GOMA include files */

#include "std.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "user_post.h"

#define GOMA_USER_POST_C

/*********** R O U T I N E S   I N   T H I S   F I L E ************************
*
*       NAME            TYPE            CALLED_BY
*    ------------             ---------               --------------
*
*    usr_post           ()   int      calc_standard_fields
*    usr_ptracking      ()   int      post_process_nodal
******************************************************************************/
/*
* NB: Dependencies on Gradients are not a problem here, because post processing
*     functions don't need sensitivities
*/  

/*********** R E C I P E   F O R   U S A G E **********************************
* The User post processing routine is responsible for:
*
*     (1)Calculating the value of that post processing variable at the current gauss 
*        integration point.    
*
*******************************************************************************/
/*
 * int user_post (param)
 *
 * ------------------------------------------------------------------------------
 * This routine is responsible for calculating a user defined variable
 * at the current gauss point:
 *     intput:    param - array of constants input on the property card.  
 *
 *     output:  post_value   => value of post-processing variable at 
 *                                                  current gauss point
 */

double 
user_post(dbl *param)		/* ptr to the user-defined parameter list */
{
  /* Local Variables */
  double post_value;

  /* int a;
     int i; */

  /* dbl X[DIM], T, C[MAX_CONC], V[DIM], P; */ /* Convenient local variables */


  static int warning = 0;

  /* Begin Execution */
 /**********************************************************/

 /* Comment out our remove this line if using this routine */
  if (warning == 0)
    {
      fprintf(stderr,"\n\n#############\n"
	    "# WARNING!! #  No user_defined post processing model implemented"
	      "\n#############\n");
      warning = 1;
    }

 /************Initialize everything for saftey**************/
  post_value = 0;                           /*Do not touch */
 /**********************************************************/

 /***********Load up convenient local variables*************/
 /*NB This ought to be done once for all fields at gauss pt*/

  /* T = fv->T; */                                    /*Do not touch */
  /*  P = fv->P;  */                                   /*Do not touch */
  /*  for(a=0; a<DIM; a++) X[a] = fv->x[a]; */		 /*Do not touch */
  /*  for(a=0; a<DIM; a++) V[a] = fv->v[a]; */		 /*Do not touch */
  /*  for(i=0; i<pd->Num_Species; i++) C[i] = fv->c[i]; */ /*Do not touch */

/* All gradients of field variables are accesable through fv->?_grad[][]
 * and material properties are accesable through mp->?
 *
 * parameters accesible through param[i], i=  0,...nparam
 */

/* Now, calculate post-process variable and put it in post_value */
/* E.G. if you want the product of temperature and concentration */
/*    if (pd->v[VELOCITY1]) */
/*        post_value = V[0] - 1.; */
/*        post_value = V[0] - log(X[1] / 10)/ log(0.9); */

/* E.G. pore radius  */
/*  if (mp->PorousMediaType == POROUS_UNSATURATED) */
/*       post_value = pmv->r_pore; */

  return post_value;
} /* End of routine usr_post                                                 */
/*****************************************************************************/
  
/*
 *   routine usr_ptracking
 *      for computing quantities along a particle trace
 *
 *      the user has access to variables loaded into the fv structure:
 *
 *              fv->
 *                      T -- temperature        (scalar)
 *                      v -- velocity           (vector)
 *                      d -- mesh displacement  (vector)
 *                      c -- concentration      (multiple scalars)
 *                      P -- pressure           (scalar)
 *                      S -- polymer stress     (tensor)
 *                      G -- velocity gradient  (tensor)
 *                     pv -- particle velocity  (vector)
 *                     pG -- particle velocity gradient (tensor)
 *
 *              fv->
 *                      grad_T -- temperature gradient
 *                      grad_P -- pressure gradient
 *                      grad_c -- species concentration gradient
 *                      grad_F -- fill gradient
 *                      grad_V -- voltage potential gradient
 *                      div_v  -- divergence of velocity
 *                      grad_v -- velocity gradient tensor
 *                      curl_v -- curl of velocity, a.k.a. vorticity
 *                      div_d  -- divergence of displacement ( dilatation )
 *                      grad_d -- gradient of displacement ( "strain" )
 *                      grad_S -- gradient of the polymer stress tensor
 *                      grad_G -- gradient of the velocity gradient tensor
 *                     grad_pv -- gradient of particle velocity
 *
 *      BEWARE OF THE ANTI-BSL CONVENTION:
 *              e.g. fv->grad_v[i][j] = dv[j]/dx[i]
 *
 */

int usr_ptracking (FILE *  jfp,         /*  filename for output */
                   const int part_id,           /*  particle id - starts at 1 */
                   const double part_x[],       /*  current particle coords   */
                   const double part_v[],       /*  current particle velocity */
                   const double part_xold[],    /*  past coords         */
                   const double part_vold[],    /*  past velocity       */
                   const int heading,           /*  flag for writing headings */
                   const double time_value,     /*  porticle time       */
                   const double time_step)      /*  time step           */

{
        if(heading)
        {
        /*  write file headings and initial values */

        fprintf(jfp," Particle Path %d\n",part_id);
        fprintf(jfp," time   x   y    vx   vy  dv11  dv12\n");
        fprintf(jfp," %g %g %g %g %g %g %g\n", time_value,
                part_x[0],part_x[1],part_v[0],part_v[1],
                fv->grad_v[0][0],fv->grad_v[0][1]);
        }
        else
        {
        /*   write data to file */

        fprintf(jfp," %g %g %g %g %g %g %g\n", time_value,
                part_x[0],part_x[1],part_v[0],part_v[1],
                fv->grad_v[0][0],fv->grad_v[0][1]);
        }

return (1);
} /* End of routine usr_ptracking                                            */
/*****************************************************************************/
/* END of file user_post.c */
/*****************************************************************************/
