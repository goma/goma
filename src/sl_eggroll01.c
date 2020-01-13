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
 * $Id: sl_eggroll01.c,v 5.2 2007-09-18 18:53:47 prschun Exp $
 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "md_timer.h"
#include "mm_eh.h"
#include "sl_auxutil.h"
#include "sl_eggroll.h"
#include "sl_eggroll_def.h"
#include "std.h"

/* Solves for leading eigenvalues
 *
 * J z = t M z
 *
 * Friendly warning: do not edit this unless you know what you are doing!!
 *
 * Originally written by Ian Gates.
 */
void
gevp_solver_rc(int nj, 
	       int mm,
	       int maxit,
	       dbl tol,
	       int filter,
	       int *ev_n,
	       dbl *ev_r,
	       dbl *ev_i,
	       dbl *ev_e,
	       dbl *ev_x,
	       int *lead,
	       int ev_jac,
	       int nev_want,
	       int *nev_found,
	       dbl **gamma,
	       dbl **evect,
	       int recycle,
	       dbl ivector,
	       int raw_residual,
	       int *rcflag,
	       int *action,
	       dbl *dwork,
	       dbl *v1,
	       dbl *v2)
{
  static int i, j, k, ip, ir, mmm, n_shift, nev_want_to_get, 
    transformation;
  static dbl xm, w1, w2, 
    *r_shift, *i_shift, 
    *rev, *iev, *res, 
    *q1, *q2, *q3, *q4, *qres;

  /* Reverse communication jumps
   */
  switch(*rcflag)
    {
    case  0: break;
    case  1: goto l01;
    case  2: goto l02;
    case 20:
    case 21:
    case 22:
    case 23:
    case 24:
    case 25:
    case 26: 
    case 30:
    case 31:
    case 32:
    case 33:
    case 34:
    case 35:
    case 36: goto l20;
    case 41: goto l41;
    case 42: goto l42;
    case 43: goto l43;
    case 44: goto l44;
    case 45: goto l45;
    case 46: goto l46;
    case 51: goto l51;
    case 52: goto l52;
    case 61: goto l61;
    case 62: goto l62;
    case 63: goto l63;
    case 64: goto l64;
    case 65: goto l65;
    case 66: goto l66;
    default:
      EH(-1, "Uh-oh!  I shouldn't be here!");
      break;
    }

  /* Initialize 
   */
  mmm = mm+5;
  transformation = 1;
  nev_want_to_get = nev_want;
  n_shift = *ev_n;
  switch(ev_jac)
    {
    case 0: /* Linear stability analysis: Jx = qMx */
      transformation = 1;
      break;
    case 1: /* Jacobian eigenproblem: Jx = qIx */
      transformation = 3;
      break;
    case 2: /* Linear stability analysis: J^Tx = qM^Tx */
      transformation = 1;
      break;
    case 3: /* Jacobian eigenproblem: J^Tx = qIx */
      transformation = 3;
      break;
    default:
      EH(-1, "Uh-oh!  I shouldn't be here!");
      break;
    }

  /* Allocate work vectors
   */
  r_shift = Dvector_birth(mmm);
  i_shift = Dvector_birth(mmm);
  rev     = Dvector_birth(mmm);
  iev     = Dvector_birth(mmm);
  res     = Dvector_birth(mmm);
  q1      = Dvector_birth(nj+5);
  q2      = Dvector_birth(nj+5);
  q3      = Dvector_birth(nj+5);
  q4      = Dvector_birth(nj+5);
  qres    = Dvector_birth(nj+5);

  /* Set shifts (real for now)
   */
  vcopy(n_shift, &r_shift[0], 1.0, &ev_r[0]);
  vcopy(n_shift, &i_shift[0], 0.0, &ev_i[0]);

  *lead = 0;

  /* Estimate of leading eigenvalue */
  /*     
  (*lead) = 0;
  for (j=0;j<n_shift;j++) {
    if (r_shift[j] > xm) {
      (*lead) = j;
      xm = r_shift[j]; }}
  if (xm == -IDG_BIG)
    { xm = -IDG_STEP; }	
  r_shift[0] = xm+0.0001;
  i_shift[0] = 0.0;
  r_shift[1] = xm-0.0001;
  i_shift[1] = 0.0; 
  */

  if ((n_shift < 1) || (n_shift > 4))
    n_shift = 4;

  /* Set initial vector
   * v0 = M [ past leading eigenvector + random number:[0,1) ]
   */
  srand48(ut());
  /* Uncomment this if you want to see the effect of not really having
   * a random initial vector, from run to run. */
  /*
  srand48((long)0);
  */
  for(i = 0; i < nj; i++)
    evect[2][i] = drand48();
  xm = 1.0/nnorm(nj, &evect[2][0]);
  vsproduct(nj, &evect[2][0], xm);
  v2sum(nj, &evect[1][0], ivector, &evect[2][0], (1.0-ivector), &v1[0]);

  /* RC flag set to 1: get v2 = J*v1
   */
  *rcflag = 1;
  *action = 1;
  vcopy(nj, &v1[0], 1.0, &evect[1][0]);
  return;
l01:
  vcopy(nj, &evect[0][0], 0.0, &v2[0]);

  /* RC flag set to 2: get v2 = M*v1
   */
  *rcflag = 2;
  *action = 2;
  vcopy(nj, &v1[0], 1.0, &evect[1][0]);
  return;
l02:
  v1add(nj, &evect[0][0], 1.0, &v2[0]);
  
  /* Arnoldi eigenvalue extractor
   */
  ir = 0;
  do {
    ir++;
    *rcflag = 10;
    if (recycle == 1)
      {
	switch (ir)
	  {
	  case 1: nev_want_to_get = 5;
	    break;
	  case 2: nev_want_to_get = nev_want;
	    break;
	  default:
	    EH(-1, "Uh-oh!  I shouldn't be here!");
	    break;
	  }
      }
    else 
      nev_want_to_get = nev_want;
    if (ir > 1) 
      recycle = 0;
l20:
    
    /* Estimate eigenvalues and eigenvectors
     */
    gevp_arnoldi_rc(nj, mm, transformation, 
		    maxit, filter, nev_want_to_get, 
		    n_shift, &r_shift[0], &i_shift[0], 
		    tol, nev_found, &rev[0], &iev[0], &res[0], 
		    evect, gamma, 
		    rcflag, action, &dwork[0], &v1[0], &v2[0]);

    /* RC flag set to 20: get inv(J-sM) if using LU
     * RC flag set to 21: get v2 = inv(J-sM)*M*v1
     * RC flag set to 22: get v2 = inv(J-sM)*M*v1
     * RC flag set to 23: get v2 = inv(J-sM)*M*v1
     * RC flag set to 30: get v2 = J*v1
     * RC flag set to 31: get v2 = M*v1
     * RC flag set to 32: get v2 = J*v1
     * RC flag set to 33: get v2 = M*v1
     * RC flag set to -1: arnoldi done
     */	
    if (*rcflag > 0)
      return;

    /* No converged eigenvalues and eigenvectors
     */
    if ((*nev_found) == 0)
      {
	/* Shifts
	 */
	vcopy(n_shift, &r_shift[0], 1.0, &rev[0]);
	vcopy(n_shift, &i_shift[0], 0.0, &iev[0]);
	
	/* Initial vector
	 */
	for(k = 0; k < nj; k++)
	  evect[mm][k] = drand48();
	for(k = 0; k < mm; k++)
	  v1add(nj, &evect[mm][0], res[k], &evect[k][0]);
	vcopy(nj, &evect[0][0], 1.0, &evect[mm][0]);
      }

    /* Recycle to improve eigenvalue and eigenvector estimates
     */
    else if(recycle == 1)
      {
	/* Set eigenspectrum
	 */
	*ev_n = (*nev_found);
	vcopy(*nev_found, &ev_r[0], 1.0, &rev[0]);
	vcopy(*nev_found, &ev_i[0], 1.0, &iev[0]);

	/* Raw residual
	 */
	for(i = 0; i < (*ev_n); i++)
	  {
	    vzero(nj, &qres[0]);
	    
	    /* Real case
	     */
	    if(ev_i[i] == 0.0)
	      {
		/* RC flag set to 41: get v2 = J*v1
		 */
		*rcflag = 41;
		*action = 1;
		vcopy(nj, &v1[0], 1.0, &evect[i][0]);
		return;
l41:
		vcopy(nj, &res[0], 1.0, &v2[0]);
		
		/* RC flag set to 42: get v2 = M*v1
		 */	
		*rcflag = 42;
		*action = 2;
		vcopy(nj, &v1[0], 1.0, &evect[i][0]);
		return;
l42:
		v1add(nj, &res[0], -ev_r[i], &v2[0]);
		ev_e[i] = nnorm(nj, &qres[0]);
	      }

	    /* Complex case
	     */
	    else
	      {
		ip = i+1;

		/* RC flag set to 43: get v2 = J*v1
		 */
		*rcflag = 43;
		*action = 1;
		vcopy(nj, &v1[0], 1.0, &evect[i][0]);
		return;
l43:
		vcopy(nj, &q1[0], 1.0, &v2[0]);

		/* RC flag set to 44: get v2 = M*v1
		 */	
		*rcflag = 44;
		*action = 2;
		vcopy(nj, &v1[0], 1.0, &evect[i][0]);
		return;
l44:
		vcopy(nj, &q2[0], 1.0, &v2[0]);

		/* RC flag set to 45: get v2 = J*v1
		 */
		*rcflag = 45;
		*action = 1;
		vcopy(nj, &v1[0], 1.0, &evect[ip][0]);
		return;
l45:
		vcopy(nj, &q3[0], 1.0, &v2[0]);

		/* RC flag set to 46: get v2 = M*v1
		 */	
		*rcflag = 46;
		*action = 2;
		vcopy(nj, &v1[0], 1.0, &evect[ip][0]);
		return;
l46:
		vcopy(nj, &q4[0], 1.0, &v2[0]);
		v3sum(nj, &qres[0], 1.0, &q1[0], -ev_r[i], &q2[0],  ev_i[i], &q4[0]);
		ev_e[i] = nnorm(nj, &qres[0]);
		v3sum(nj, &qres[0], 1.0, &q3[0], -ev_r[i], &q4[0], -ev_i[i], &q2[0]);
		ev_e[i] += nnorm(nj, &qres[0]);
		ev_e[ip] = ev_e[i]; 
		i++;
	      }
	  }
	
	/* Re-order eigenvalues and eigenvectors
	 */
	gevp_order(nj, 
		   *ev_n, &ev_r[0], &ev_i[0], &ev_e[0], &ev_x[0], 
		   evect, gamma);

	/* Leading eigenvalue
	 */
	(*lead) = 0;
	xm = ev_r[0];
	for(j = 1; j < (*nev_found); j++)
	  if((ev_r[i] > xm) && (fabs(ev_e[i]) < tol))
	    {
	      (*lead) = j;
	      xm = rev[j];
	    }
	(*nev_found) = 0;

	/* Shifts
	 */
	vcopy(n_shift, &r_shift[0], 1.0, &rev[0]);
	vcopy(n_shift, &i_shift[0], 0.0, &iev[0]);
	r_shift[0] = xm+IDG_STEP;
	i_shift[0] = 0.0;
	r_shift[1] = xm-IDG_STEP;
	i_shift[1] = 0.0;
	r_shift[2] = xm-1.0;
	i_shift[2] = 0.0;

	/* Initial vector
	 */
	w1 = 0.001;
	w2 = 1.0 - w1;
	vzero(nj, &evect[mm][0]);

	for(k = 0; k < *nev_found; k++)
	  v1add(nj, &evect[mm][0], res[k], &evect[k][0]);

	if((xm = nnorm(nj, &evect[mm][0])) == 0.0)
	  xm = 1.0;

	vcopy(nj, &evect[1][0], w1/xm, &evect[mm][0]);
	for(k = 0; k < nj; k++)
	  evect[mm][k] = drand48();

	v1add(nj, &evect[1][0], w2, &evect[mm][0]);

	/* RC flag set to 51: get v2 = J*v1
	 */
	*rcflag = 51;
	*action = 1;
	vcopy(nj, &v1[0], 1.0, &evect[1][0]);
	return;
l51:
	vcopy(nj, &evect[0][0], 1.0, &v2[0]);

	/* RC flag set to 52: get v2 = M*v1
	 */
	*rcflag = 52;
	*action = 2;
	vcopy(nj, &v1[0], 1.0, &evect[1][0]);
	return;
l52:
	v1add(nj, &evect[0][0], 1.0, &v2[0]);
      }
  } while ((recycle == 1) && (ir < maxit));

  /* Exit if no eigenvalues found
   */
  if (*nev_found == 0)
    {
      puts("\n\n  Error:\n");
      puts("    Convergence not achieved in gevp_arnoldi().");
      puts("    Program terminated.\n");
      exit(0);
    }
  if (*nev_found >  nev_want)
    *nev_found = nev_want;

  /* Set new eigenspectrum and backup old (real) eigenvalues 
   * (for future shifts)
   */
  *ev_n = (*nev_found);
  vcopy(*nev_found, &ev_r[0], 1.0, &rev[0]);
  vcopy(*nev_found, &ev_i[0], 1.0, &iev[0]);
  vcopy(*nev_found, &ev_e[0], 1.0, &res[0]);

  /* Raw residual
   */
  if(raw_residual == 1)
    {
      for(i = 0; i < (*ev_n); i++)
	{
	  vzero(nj, &qres[0]);
	  /* Real case
	   */
	  if(ev_i[i] == 0.0)
	    {
	    /* RC flag set to 61: get v2 = J*v1
	     */
	      *rcflag = 61;
	      *action = 1;
	      vcopy(nj, &v1[0], 1.0, &evect[i][0]);
	      return;
l61:
	      vcopy(nj, &qres[0], 1.0, &v2[0]);

	      /* RC flag set to 62: get v2 = M*v1
	       */	
	      *rcflag = 62;
	      *action = 2;
	      vcopy(nj, &v1[0], 1.0, &evect[i][0]);
	      return;
l62:
	      v1add(nj, &qres[0], -ev_r[i], &v2[0]);
	      ev_e[i] = nnorm(nj, &qres[0]);
	    }

	  /* Complex case
	   */
	  else
	    {
	      ip = i+1;

	      /* RC flag set to 63: get v2 = J*v1
	       */
	      *rcflag = 63;
	      *action = 1;
	      vcopy(nj, &v1[0], 1.0, &evect[i][0]);
	      return;
l63:
	      vcopy(nj, &q1[0], 1.0, &v2[0]);

	      /* RC flag set to 64: get v2 = M*v1
	       */	
	      *rcflag = 64;
	      *action = 2;
	      vcopy(nj, &v1[0], 1.0, &evect[i][0]);
	      return;
l64:
	      vcopy(nj, &q2[0], 1.0, &v2[0]);

	      /* RC flag set to 65: get v2 = J*v1
	       */
	      *rcflag = 65;
	      *action = 1;
	      vcopy(nj, &v1[0], 1.0, &evect[ip][0]);
	      return;
l65:
	      vcopy(nj, &q3[0], 1.0, &v2[0]);

	      /* RC flag set to 66: get v2 = M*v1
	       */	
	      *rcflag = 66;
	      *action = 2;
	      vcopy(nj, &v1[0], 1.0, &evect[ip][0]);
	      return;
l66:
	      vcopy(nj, &q4[0], 1.0, &v2[0]);

	      v3sum(nj, &qres[0], 1.0, &q1[0], -ev_r[i], &q2[0],  ev_i[i], &q4[0]);
	      ev_e[i] = nnorm(nj, &qres[0]);
	      v3sum(nj, &qres[0], 1.0, &q3[0], -ev_r[i], &q4[0], -ev_i[i], &q2[0]);
	      ev_e[i] += nnorm(nj, &qres[0]);
	      ev_e[ip] = ev_e[i]; 
	      i++;
	    } /* complex case */
	} /* for(i = 0; i < (*ev_n); i++) */
    } /* if(raw_residual == 1) */

  /* Re-order eigenvalues and eigenvectors - increasing modulus
   */
  gevp_order(nj, 
	     *ev_n, &ev_r[0], &ev_i[0], &ev_e[0], &ev_x[0], 
	     evect, gamma);
  
  /* Filter out infinite eigenvalues
   * ev_x[] < 0 means the eigenvalue is "infinite"
   */
  for(i = 0; i < (*ev_n); i++)
    {
      res[i] = 1.0;
      ev_x[i] = sqrt(SQUARE(ev_r[i])+SQUARE(ev_i[i]));
    }
  xm = ev_x[0];
  if((w1 = log10(fabs(xm))) < 1.0)
    w1 = 1.0;
  for(i = 0; i < (*ev_n); i++)
    {
      if ((w2 = log10(fabs(ev_x[i]))) < 1.0)
	w2 = 1.0;
      if (w2/w1 > 3.0)
	res[i] = -1.0;
    }
  for(i = 0; i < (*ev_n); i++)
    ev_x[i] = res[i];

  /* Find leading eigenvalue and eigenvector id
   */
  (*lead) = 0;
  xm = ev_r[0];
  for(i = 1; i < (*ev_n); i++)
    if((ev_r[i] > xm) && (fabs(ev_e[i]) < tol) && (ev_x[i] > 0.0))
      {
	(*lead) = i;
	xm = ev_r[i];
      }

  /* De-allocate memory
   */
  Dvector_death(&r_shift[0], mmm);
  Dvector_death(&i_shift[0], mmm);
  Dvector_death(&rev[0], mmm);
  Dvector_death(&iev[0], mmm);
  Dvector_death(&res[0], mmm);
  Dvector_death(&q1[0], nj+5);
  Dvector_death(&q2[0], nj+5);
  Dvector_death(&q3[0], nj+5);
  Dvector_death(&q4[0], nj+5);
  Dvector_death(&qres[0], nj+5);

  /* End 
   */
  *rcflag = 0;
  return;
}
/* END of routine gevp_solver_rc */

/******************************************************************************/
/* END of file sl_eggroll01.c */
/******************************************************************************/

