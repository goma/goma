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
 * $Id: sl_eggroll02.c,v 5.2 2007-09-18 18:53:47 prschun Exp $
 */

/*
 * $Log: not supported by cvs2svn $
 * Revision 5.1  2007/04/04 21:15:46  tabaer
 *
 * Changes norm to nnorm calls to avoid conflict with chaco.
 *
 * Revision 5.0  2006/07/06 16:18:57  edwilke
 *
 * New Goma version: 'Farewell CRMPR, hello $3 a gallon gas!'
 *
 * Revision 4.3  2004/01/12 18:51:24  tabaer
 * More minor changes for Darwin compatibility.
 *
 * Revision 4.2  2003/11/25 23:16:02  dalabre
 * The copyright statement has been updated for Goma and the version ratcheted
 * to 4.4.0.
 *
 * The makefile Goma.mk has been updated for Linux and Sun so that these
 * versions can be built easily (until the configure-make production is
 * complete); the Linux version is default.
 *
 * Added a prototype for function assemble_interface_extension_velocity_sic
 * in mm_fill_terms.h.
 *
 * Revision 4.1  2003/09/23 18:19:34  drnoble
 * Another fine test of GOMA's configuration management.  Here sl_umfutil.c
 * and sl_umfutil.h are removed.  (sl_umfutil.c wasn't being compiled
 * previously but sl_umfutil.h was being used for prototypes for code
 * in sl_auxutil.c.  Got it?) sl_umf.h is created for prototypes of code
 * in sl_umf.c.
 *
 * Revision 4.0  2001/12/21 06:01:54  dalabre
 * Up Goma source code repository to V4.0. This identification 'coincides'
 * with our documentation upgrade. It will be tagged "Tora_Bora" in
 * recognition of the deep, unknown recesses that still remain in Goma.
 *
 * Revision 3.7  2000/05/18 05:41:17  dalabre
 * Multi-Platform changes, fixes and corrections (Part 1).
 * --------------------------------------------------------------------------
 *
 * Revision 3.6  2000/01/14 17:57:12  mmhopki
 * Ooops!  Don't need to exit() after EH(-1, ).  rf_eigensolver.h is chucked.
 *
 * Revision 3.5  2000/01/14 17:49:39  mmhopki
 * Mongo update:
 * - Included missing copyright notices, log, and id strings.
 * - Stylized the sl_eggroll* files to look more like C.
 * - Cleaned up some header file redundancy.
 * - Put in default: cases missing in switch statements.
 * - Removed unused or redundant code.
 * - Double -> dbl conversions.
 * - De-uppercased (?) some non-#defined variables (LINEAR_STABILITY,
 *   FILTER, VISC_SENS).
 * - Various other style updates to increase conformity to the prevailing
 *   Goma style.
 *
 */


#include <stdio.h>
#include <math.h>

#include "mm_eh.h"
#include "sl_auxutil.h"
#include "sl_eggroll.h"
#include "sl_eggroll_def.h"
#include "std.h"

/* Iterative arnoldi eigenvalue extraction that uses shift-and-invert
 * strategy
 *
 * Friendly warning: do not edit this unless you know what you are doing!!
 *
 * Originally written by Ian Gates.
 */
void
gevp_arnoldi_rc(int nj,
		int mm,
		int transformation,
		int maxit,
		int filter,
		int nev_want_to_get,
		int n_shift,
		dbl *r_shift,
		dbl *i_shift,
		dbl tol,
		int *nev_found,
		dbl *rev,
		dbl *iev,
		dbl *res,
		dbl **evect,
		dbl **schur,
		int *rcflag,
		int *action,
		dbl *dwork,
		dbl *v1,
		dbl *v2)
{
  static int i, j, k, l, ic, ip, jp, j0, j1, mmm, wts, 
    pass_number, new_shift, n_sigma, 
    nev_converged, nev_converged_last_pass, 
    nev_to_get_left,
    *order;
  static dbl tm, 
    w1, w2, w4, w5, w6, wa, wb, 
    xm, beta, hnorm, r_sigma, i_sigma,
    r_ev_a, i_ev_a, r_ev_b, i_ev_b,
    *q1, *q2, 
    **hh, **rr;
  static EV ev;

  /* Reverse communication jumps
   */
  switch(*rcflag)
    {
    case 10: break;
    case 20: goto l20;
    case 21: goto l21;
    case 22: goto l22;
    case 23: goto l23;
    case 30: goto l30;
    case 31: goto l31;
    case 32: goto l32;
    case 33: goto l33;
    default: 
      EH(-1, "Uh-oh!  I shouldn't be here!");
      break;
    }

  /* Initialize
   */
  wts           = 1;
  pass_number   = 0;
  new_shift     = 1;
  nev_converged = nev_converged_last_pass = 0;
  if(wts != 0) 
    puts("\n Arnoldi Eigenvalue Extractor\n");

  /* Allocate work vectors
   */
  q1    = Dvector_birth(nj+5);
  q2    = Dvector_birth(nj+5);
  if(wts != 0) 
    puts(" Allocated work vectors");

  /* Allocate Hessenberg and Schur matrices
   */
  mmm   = mm+5;
  hh    = Dmatrix_birth(mmm, mmm);
  rr    = Dmatrix_birth(mmm, mmm);
  order = Ivector_birth(mmm);
  ev.r  = Dvector_birth(mmm);
  ev.i  = Dvector_birth(mmm);
  ev.e  = Dvector_birth(mmm);
  ev.u  = Dmatrix_birth(mmm, mmm);
  if(wts != 0) 
    puts(" Allocated work martices");

  /* Initial shift
   */
  n_sigma = 0;
  r_sigma = r_shift[n_sigma];
  i_sigma = i_shift[n_sigma];

  /* Initial vector
   */
  for(i = 0; i < mm; i++)
    vzero(nj, &schur[i][0]);
  vcopy(nj, &schur[0][0], 1.0, &evect[0][0]);

  do {
    pass_number++;
    if(wts != 0)
      switch(new_shift)
	{
	case 0: printf(" Pass %2d\n", pass_number);
	  break;
	case 1: printf(" Pass %2d  New Shift = % 10.6e %+10.6e i\n", 
		       pass_number, r_sigma, i_sigma);
	break;
	default:
	  EH(-1, "Uh-oh!  I shouldn't be here!");
	  break;
	}

    /* Formation and factorization of shifted matrix
     */
    if (new_shift == 1)
      {
	/* RC flag set to 20: get inv(J-sM)
	 */
	*rcflag = 20;
	*action = 3;
	dwork[0] = r_sigma;
	dwork[1] = i_sigma;
	return;
      }
l20:
#ifdef DEBUG_gevp_arnoldi
    puts(" LU DECOMPOSITION");
#endif

    /* Filtering
     */
    if(pass_number == 1)
      for(i = 0; i < filter; i++)
	{
	  /* RC flag set to 21: get inv(J-sM)
	   */
	  *rcflag = 21;
	  *action = 4;
	  dwork[0] = r_sigma;
	  dwork[1] = i_sigma;
	  vcopy(nj, &v1[0], 1.0, &schur[0][0]);
	  return;
l21:
	  vcopy(nj, &schur[0][0], 1.0, &v2[0]);
	}
#ifdef DEBUG_gevp_arnoldi
    puts(" FILTER");
#endif

    /* Deflation
     */
    k = 0;
    if(new_shift == 1)
      {
	j0 = 0;
	for(i = 0; i < mmm; i++)
	  vzero(mmm, &rr[i][0]); 
      }
    else
      if(nev_converged_last_pass != 0) 
	j0 = nev_converged-nev_converged_last_pass;
      else 
	k = 1;

    if(k == 0)
      for(j = j0; j < nev_converged; j++)
	{
	  /* RC flag set to 22: get inv(J-sM)
	   */
	  *rcflag = 22;
	  *action = 4;
	  dwork[0] = r_sigma;
	  dwork[1] = i_sigma;
	  vcopy(nj, &v1[0], 1.0, &schur[j][0]);
	  return;
l22:
	  if (ev.i[j] > 0.0)
	    j1 = j+2;
	  else 
	    j1 = j+1;
	  for(i = 0; i < j1; i++)
	    rr[i][j] = dot_product(nj, &v2[0], &schur[i][0]);
	}
    for(i = 0; i < mmm; i++)
      vzero(mmm, &hh[i][0]);
    for(i = 0; i < mm; i++)
      vcopy(nev_converged, &hh[i][0], 1.0, &rr[i][0]);
#ifdef DEBUG_gevp_arnoldi
    puts(" DEFLATION");
#endif
    
    /* Arnoldi iteration
     */
    nev_to_get_left = nev_want_to_get-nev_converged;
    
    /* Normalize initial vector
     */
    xm = nnorm(nj, &schur[nev_converged][0]);
#ifdef DEBUG_gevp_arnoldi
    printf("Norm of Schur vector is %g.\n", xm);
#endif
    vcopy(nj, &q1[0], 1.0, &schur[nev_converged][0]);
    xm = nnorm(nj, &q1[0]);
    vcopy(nj, &schur[nev_converged][0], 1.0/xm, &q1[0]);
    
    /* Form KSS
     */
    xm = 0.0;
    for(i = nev_converged; i < mm; i++)
      {
	ip = i+1;

	/* Shift-invert matrix-vector product
	 * RC flag set to 23: get inv(J-sM)
	 */
	*rcflag = 23;
	*action = 4;
	dwork[0] = r_sigma;
	dwork[1] = i_sigma;
	vcopy(nj, &v1[0], 1.0, &schur[i][0]);
	return;
l23:
	vcopy(nj, &schur[ip][0], 1.0, &v2[0]);

	/* Orthogonalize schur[ip][] against all previous vectors schur[i][], i<ip
	 */
	for(j = 0; j < ip; j++)
	  hh[j][i] = 0.0;
#ifdef DEBUG_gevp_arnoldi
	puts(" MVP");
#endif
	j = 0;
	do {
	  j++;
	  hnorm = 0.0;
	  for(k = 0; k < ip; k++)
	    {
	      tm = dot_product(nj, &schur[k][0], &schur[ip][0]);
	      hnorm += SQUARE(tm); 
	      hh[k][i] += tm;
	      v1add(nj, &schur[ip][0], -tm, &schur[k][0]);
	    }
	  tm = dot_product(nj, &schur[ip][0], &schur[ip][0]);
	} while((j < 2) && (hnorm > 10.0*tm));
	tm = sqrt(tm);
	hh[ip][i] = tm;
	if(tm == 0.0)
	  tm = 1.0;
	tm = 1.0/tm;
	vsproduct(nj, &schur[ip][0], tm);
	for(j = 0; j < ip; j++)
	  xm += SQUARE(hh[j][i]);
      }
    beta = hh[mm][mm-1];
    xm = sqrt(xm);
#ifdef DEBUG_gevp_arnoldi
    puts(" KSS");
#endif

    /* QR method for eigenvalue extraction
     */
    eigenvv(mm, hh, &ev);
#ifdef DEBUG_gevp_arnoldi
    puts(" QR");
#endif

    /* Residuals
     */
    k = mm-1;
    for(j = 0; j < mm; j++)
      {
	/* Real case
	 */
	if(ev.i[j] == 0.0)
	  {
	    w1 = 0.0;
	    for(i = 0; i < mm; i++)
	      w1 += SQUARE(ev.u[i][j]);
	    w1 = sqrt(w1);
	    ev.e[j] = fabs(beta*ev.u[k][j]/w1);
	  }

	/* Complex case
	 */
	else
	  {
	    jp = j+1;
	    w1 = w2 = 0.0;
	    for(i = 0; i < mm; i++)
	      {
		w1 += SQUARE(ev.u[i][j]);
		w2 += SQUARE(ev.u[i][jp]);
	      }
	    w1 = sqrt(w1);
	    w2 = sqrt(w2);
	    ev.e[j]  = fabs(beta)*sqrt(SQUARE(ev.u[k][j])+SQUARE(ev.u[k][jp]))/(w1+w2);
	    ev.e[jp] = ev.e[j];
	    j++;
	  }
      } /* for(j = 0; j < mm; j++) */
#ifdef DEBUG_gevp_arnoldi
    puts(" RES");
#endif

    /* Ordering sequences - eliminate previously converged eigenpairs
     */
    vzero(nev_converged, &ev.r[0]);
    vzero(nev_converged, &ev.i[0]);
    vinit(nev_converged, &ev.e[0], IDG_BIG);

    /* Reorder eigenvalues - decreasing moduli
     */
    _heapsort(mm, &ev.r[0], &ev.i[0], &ev.e[0], &order[0], 0);

    /* Reorder eigenvectors
     */
    for(i = 0; i < mm; i++)
      {
	k = order[i];
	for(j = 0; j < mm; j++)
	  hh[j][i] = ev.u[j][k];
      }
    for(i = 0; i < mm; i++)
      vcopy(mm, &ev.u[i][0], 1.0, &hh[i][0]);

    /* Make sure eliminated all complex pairs from permuted eigenvalue lists
     */
    if(ev.i[nev_to_get_left] > 0.0)
      nev_to_get_left++;

    /* Reorder eigenvalues - decreasing residual
     */
    _heapsort(mm, &ev.r[0], &ev.i[0], &ev.e[0], &order[0], 1);

    /* Reorder eigenvectors
     */
    for(i = 0; i < mm; i++)
      {
	k = order[i];
	for(j = 0; j < mm; j++)
	  hh[j][i] = ev.u[j][k];
      }
    for(i = 0; i < mm; i++)
      vcopy(mm, &ev.u[i][0], 1.0, &hh[i][0]);

    /* Determine converged eigenvalues
     */
    nev_converged_last_pass = 0;
    for(j = 0; j < nev_to_get_left; j++)
      {
	tm = sqrt(SQUARE(ev.r[j])+SQUARE(ev.i[j]));
	ev.e[j] *= 2.0/tm;
	if((fabs(ev.e[j]) < tol) && (nev_converged_last_pass == j))
	  nev_converged_last_pass++;
      }

    /* Determine eigenvectors 
     */
    for(j = 0; j < nj; j++)
      {
	for(i = 0; i < nev_to_get_left; i++)
	  {
	    xm = 0.0;
	    k = i;
	    for(l = 0; l < mm; l++) 
	      xm += schur[l][j]*ev.u[l][k];
	    q1[i] = xm;
	  }
	for(i = 0; i < nev_to_get_left; i++)
	  evect[nev_converged+i][j] = schur[nev_converged+i][j] = q1[i];
      }
#ifdef DEBUG_gevp_arnoldi
    puts(" EVECTORS1");
#endif
    
    /* Rayleigh quotient
     */
    for(k = 0; k <= nev_converged_last_pass; k++)
      {
	i = k+nev_converged;
	ip = i+1;

	/* RC flag set to 30: get v2 = J*v1
	 */
	*rcflag = 30;
	*action = 1;
	vcopy(nj, &v1[0], 1.0, &schur[i][0]);
	return;
l30:
	vcopy(nj, &q1[0], 1.0, &v2[0]);
	w1 = dot_product(nj, &schur[i][0], &schur[i][0]);
	w2 = dot_product(nj, &q1[0], &schur[i][0]);
	xm = dot_product(nj, &q1[0], &q1[0]);

	/* Real eigenvalue
	 */
	if(ev.i[k] == 0.0)
	  {
	    r_ev_a = w2/w1;
	    i_ev_a = 0.0;
	  }

	/* Complex eigenvalue
	 */
	else
	  {
	    /* RC flag set to 31: get v2 = J*v1
	     */
	    *rcflag = 31;
	    *action = 1;
	    vcopy(nj, &v1[0], 1.0, &schur[ip][0]);
	    return;
l31:
	    vcopy(nj, &q2[0], 1.0, &v2[0]);
	    w1 += dot_product(nj, &schur[ip][0], &schur[ip][0]);
	    w4  = dot_product(nj, &q1[0], &schur[ip][0]);
	    w5  = dot_product(nj, &q2[0],  &schur[i][0]);
	    w6  = dot_product(nj, &q2[0], &schur[ip][0]);
	    xm += dot_product(nj, &q2[0], &q2[0]);
	    r_ev_a = (w2+w6)/w1;
	    i_ev_a = fabs((w5-w4)/w1);
	  }
	wa = sqrt(xm/w1);
	wb = sqrt(SQUARE(r_ev_a)+SQUARE(i_ev_a));
	if(wa < IDG_EPS * wb) 
	  wa = 0.0;

	/* RC flag set to 32: get v2 = M*v1
	 */
	*rcflag = 32;
	*action = 2;
	vcopy(nj, &v1[0], 1.0, &schur[i][0]);
	return;
l32:
	vcopy(nj, &q1[0], 1.0, &v2[0]);
	w1 = dot_product(nj, &schur[i][0], &schur[i][0]);
	w2 = dot_product(nj, &q1[0], &schur[i][0]);
	xm = dot_product(nj, &q1[0], &q1[0]);

	/* Real eigenvalue
	 */
	if(ev.i[k] == 0.0)
	  {
	    r_ev_b = w2/w1;
	    i_ev_b = 0.0;
	  }

	/* Complex eigenvalue
	 */
	else
	  {
	    /* RC flag set to 33: get v2 = M*v1
	     */
	    *rcflag = 33;
	    *action = 2;
	    vcopy(nj, &v1[0], 1.0, &schur[ip][0]);
	    return;
l33:
	    vcopy(nj, &q2[0], 1.0, &v2[0]);
	    w1 += dot_product(nj, &schur[ip][0], &schur[ip][0]);
	    w4  = dot_product(nj, &q1[0], &schur[ip][0]);
	    w5  = dot_product(nj, &q2[0],  &schur[i][0]);
	    w6  = dot_product(nj, &q2[0], &schur[ip][0]);
	    xm += dot_product(nj, &q2[0], &q2[0]);
	    r_ev_b = (w2+w6)/w1;
	    i_ev_b = fabs((w5-w4)/w1);
	  }
	wa = sqrt(xm/w1);
	wb = sqrt(SQUARE(r_ev_b)+SQUARE(i_ev_b));
	if (wa < IDG_EPS * wb) 
	  wa = 0.0;
	tm = SQUARE(r_ev_b)+SQUARE(i_ev_b);
	rev[i] = (r_ev_a*r_ev_b+i_ev_a*i_ev_b)/tm;
	iev[i] = (i_ev_a*r_ev_b-r_ev_a*i_ev_b)/tm;
	res[i] = ev.e[k];
	if(ev.i[k] != 0.0)
	  {
	    rev[ip] =  rev[i];
	    iev[ip] = -iev[i];
	    res[ip] =  res[i];
	    k++;
	  }
      } /* for(k = 0; k <= nev_converged_last_pass; k++) */
#ifdef DEBUG_gevp_arnoldi
    puts(" RAYLEIGH QUOTIENT");
#endif

    /* Orthogonalize newly converged eigenvectors against previously
     * converged eigenvectors
     */
    for(i = 0; i < nev_converged_last_pass; i++)
      {
	ip = i+nev_converged;
	j = 0;
	do {
	  j++;
	  hnorm = 0.0;
	  for(k = 0; k < ip; k++)
	    {
	      tm = dot_product(nj, &schur[k][0], &schur[ip][0]);
	      hnorm += SQUARE(tm);
	      v1add(nj, &schur[ip][0], -tm, &schur[k][0]);
	    }
	  tm = dot_product(nj, &schur[ip][0], &schur[ip][0]);
	} while((j < 2) && (hnorm > 10.0*tm));
	tm = sqrt(tm);
	if (tm == 0.0)
	  tm = 1.0;
	tm = 1.0/tm;
	vsproduct(nj, &schur[ip][0], tm);
      }
#ifdef DEBUG_gevp_arnoldi
    puts(" EVECTORS2");
#endif

    /* Increment counter
     */
    nev_converged += nev_converged_last_pass;

    /* New shift
     */
    new_shift = 0;
    xm = rev[nev_converged-1];

    /* Set shift to close to leading eigenvalue
     */
    i = 0;
    for(j = 0; j < nev_converged; j++)
      if((rev[j] > xm) && (fabs(res[j]) < tol))
	{
	  i = j;
	  xm = rev[j];
	}
    if((nev_converged_last_pass == 0) && (pass_number == 1))
      {
	new_shift = 1;
	n_sigma++;
	r_sigma = r_shift[n_sigma];
	i_sigma = i_shift[n_sigma];
      }
    else if((nev_converged_last_pass == 0) && (n_sigma < n_shift))
      {
	new_shift = 1;
	r_sigma = xm+IDG_EPS;
	i_sigma = 0.0;
      }
    else
      new_shift = 0;

    /* Print latest pass results
     */
    if(wts > 1)
      {
	if(nev_converged_last_pass > 0) 
	  printf(" ARN: %d new eigenvalues converged.\n", 
		 nev_converged_last_pass);
	for(j = 0; j < nev_converged; j++)
	  printf(" EV[%3d]= % 10.6e  %+10.6e i\n", j, rev[j], iev[j]);
	puts(" *********************************************************************");
      }

    /* Exit condition
     */
    ic = 0;
    if((nev_converged < nev_want_to_get) && (pass_number < maxit))
      ic = 1;
  } while (ic == 1);
  (*nev_found) = nev_converged;
  
  /* De-allocate memory
   */
  Dvector_death(&q1[0], nj+5);
  Dvector_death(&q2[0], nj+5);
  Dmatrix_death(rr, mmm, mmm);
  Dmatrix_death(hh, mmm, mmm);
  Ivector_death(&order[0], mmm);
  Dvector_death(&ev.r[0], mmm);
  Dvector_death(&ev.i[0], mmm);
  Dvector_death(&ev.e[0], mmm);
  Dmatrix_death(ev.u, mmm, mmm);
  if (wts != 0)
    puts(" De-allocated work storage");

  /* End
   */
  if (wts != 0)
    puts(" ARN-IT done. ");
  *rcflag = -1;
  return;
}
/* END of routine gevp_arnoldi_rc */

/******************************************************************************/
/* END of file sl_eggroll02.c */
/******************************************************************************/




