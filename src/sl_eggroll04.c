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
 * $Id: sl_eggroll04.c,v 5.1 2007-09-18 18:53:47 prschun Exp $
 */

/*
 * $Log: not supported by cvs2svn $
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
 * Revision 4.0  2001/12/21 06:01:55  dalabre
 * Up Goma source code repository to V4.0. This identification 'coincides'
 * with our documentation upgrade. It will be tagged "Tora_Bora" in
 * recognition of the deep, unknown recesses that still remain in Goma.
 *
 * Revision 3.5  2000/05/18 05:41:17  dalabre
 * Multi-Platform changes, fixes and corrections (Part 1).
 * --------------------------------------------------------------------------
 *
 * Revision 3.4  2000/01/14 17:57:13  mmhopki
 * Ooops!  Don't need to exit() after EH(GOMA_ERROR, ).  rf_eigensolver.h is chucked.
 *
 * Revision 3.3  2000/01/14 17:49:40  mmhopki
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

/* gets eigenvalues for matrix a[][] real and imaginary eigenvalues in
 * e.r[],e.i[] vectors.
 *
 * Friendly warning: do not edit this unless you know what you are
 * doing!!
 *
 * C version by Ian Gates, adapted from NR [Press et al 1986].
 */
void
eigenvv(int n,
	dbl **a,
	EV *e)
{
  int i, j, nn, 
    *ich;
  dbl xm, 
    *scl;

  /* Allocate
   */
  ich = Ivector_birth(n+5);
  scl = Dvector_birth(n+5);

  /* Initialize
   */
  nn = n+1;
  for(i = 0; i < nn; i++)
    {
      ich[i] = 0;
      scl[i] = 1.0;
    }
  for(i = n; i > 0; i--)
    for(j = n; j > 0; j--)
      a[i][j] = a[i-1][j-1];
  balanc(n, a, &scl[0]);
  elmhes(n, a, &ich[0]);
  eltran(n, a, &ich[0], e);
  eighqr(n, a, e);
  balbak(n, &scl[0], e);
  for(i = 0; i < n; i++)
    {
      e->r[i] = e->r[i+1];
      e->i[i] = e->i[i+1];
    }
  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      e->u[i][j] = e->u[i+1][j+1];

  /* Normalize (real) eigenvectors
   */
  for(i = 0; i < n; i++)
    if(e->i[i] == 0.0)
      {
	xm = 0.0;
	for(j = 0; j < n; j++)
	  xm += SQUARE(e->u[j][i]);
	xm = sqrt(xm);
	for(j = 0; j < n; j++)
	  e->u[j][i] /= xm;
      }

  /* De-allocate
   */
  Ivector_death(&ich[0], n+5);
  Dvector_death(&scl[0], n+5);
}

void
balanc(int n,
       dbl **a,
       dbl *scl)
{
  int i, j, nn, last;
  dbl radix, s, r, g, f, c, sqrdx;

  /* Initialize
   */
  nn = n+1;
  radix = 16.0;
  sqrdx = SQUARE(radix);
  for (i=0;i<nn;i++)
    scl[i] = 1.0;

  /* Balance
   */
  do {
    last = 0;
    for(i = 1; i <= n; i++)
      {
	c = r = 0.0;
	for(j = 1; j <= n; j++)
	  if(j != i)
	    {
	      c += fabs(a[j][i]);
	      r += fabs(a[i][j]);
	    }
	if((c == 0.0) || (r == 0.0))
	  continue;
	g = r/radix;
	f = 1.0;
	s = c+r;
	while(c < g)
	  {
	    f *= radix;
	    c *= sqrdx;
	  }
	g = r*radix;
	while(c > g)
	  {
	    f /= radix;
	    c /= sqrdx;
	  }
	if((c+r)/f > 0.95*s)
	  continue;
	last = 1;
	g = 1.0/f;
	scl[i] *= f;
	for(j = 1; j <= n; j++)
	  a[i][j] *= g;
	for(j = 1; j <= n; j++)
	  a[j][i] *= f;
      } /* for(i = 1; i <= n; i++) */
  } while(last == 1);
}

void
elmhes(int n,
       dbl **a,
       int *ich)
{
  int i, j, m, mm;
  dbl x, y;

  if(n <= 2)
    return;
  for(m = 2; m < n; m++)
    {
      mm = m-1;
      x = 0.0;
      i = m;
      for(j = m; j <= n; j++)
	if(fabs(a[j][mm]) > fabs(x))
	  {
	    x = a[j][mm];
	    i = j;
	  }
      ich[m] = i;
      if(i != m)
	{
	  for(j = mm; j <= n; j++)
	    {
	      y = a[i][j];
	      a[i][j] = a[m][j];
	      a[m][j] = y;
	    }
	  for(j = 1; j <= n; j++)
	    {
	      y = a[j][i];
	      a[j][i] = a[j][m];
	      a[j][m] = y;
	    }
	}
      if(x == 0.0)
	continue;
      for(i = m+1; i <= n; i++)
	{
	  y = a[i][mm];
	  if(y == 0.0)
	    continue;
	  y /= x;
	  a[i][mm] = y;
	  for(j = m; j <= n; j++)
	    a[i][j] -= y*a[m][j];
	  for(j = 1; j <= n; j++)
	    a[j][m] += y*a[j][i];
	}
    } /* for(m = 2; m < n; m++) */
}

void
eltran(int n,
       dbl **a,
       int *ich,
       EV *e)
{
  int i, j, k, l, m, nn;

  /* Initialize
   */
  nn = n+1;
  for(i = 0; i < nn; i++)
    {
      for(j = 0; j < nn; j++)
	e->u[i][j] = 0.0;
      e->u[i][i] = 1.0;
    }
  for(j = n-1; j > 1; j--)
    {
      k = j+1;
      for(l = k; l <= n; l++) 
	e->u[l][j] = a[l][j-1];
      i = ich[j];
      if(i != j)
	{
	  for(m = j; m <= n; m++)
	    {
	      e->u[j][m] = e->u[i][m];
	      e->u[i][m] = 0.0;
	    }
	  e->u[i][j] = 1.0;
	}
    }
}

void
balbak(int n,
       dbl *scl,
       EV *e)
{
  int i, j;
  dbl s;
  for(i = 1; i <= n; i++)
    {
      s = scl[i];
      for(j = 1; j <= n; j++)
	e->u[i][j] *= s;
    }
}	

void
eighqr(int n,
       dbl **a,
       EV *e)
{
  int i, j, k, l, m,
    ni, nl = 0, nn, nm = 0, mp, its = 0, mmin;
  dbl anorm, p, q, r, s, t, w, x, y, z, t1, t2, ra, sa, vr, vi, w1, w2;

  /* Initialize
   */
  ni = 0;
  nn = n;
  p = q = r = s = t = w = x = y = z = 0.0;

  /* Matrix norm
   */
  anorm = fabs(a[1][1]);
  for(i = 2; i <= n; i++)
    for(j = i-1; j <= n; j++)
      anorm += fabs(a[i][j]);

  /* Find eigenvalues
   */
  while(nn > 0)
    {
      if(ni == 0)
	{
	  its = 0;
	  nm = nn-1;
	  nl = nm-1;
	}
      for(l = nn; l > 0; l--)
	{
	  if(l == 1)
	    break;
	  s = fabs(a[l-1][l-1]) + fabs(a[l][l]);
	  if(s == 0.0)
	    s = anorm;
	  t1 = s;
	  t2 = t1+fabs(a[l][l-1]);
	  if(t2 == t1)
	    break;
	}
      x = a[nn][nn];
      y = a[nm][nm];
      w = a[nn][nm]*a[nm][nn];
      ni = 0;
      if(l == nn)
	{
	  e->r[nn] = a[nn][nn] = x+t;
	  e->i[nn] = 0.0;
	  nn = nm;
	}
      else if (l == nm)
	{
	  p = 0.5*(y-x);
	  q = SQUARE(p)+w;
	  z = sqrt(fabs(q));
	  a[nn][nn] = x+t;
	  a[nm][nm] = y+t;
	  x += t;
	  if(q >= 0.0)
	    {
	      z = p + SGN(p) * fabs(z);
	      e->r[nn] = e->r[nm] = x+z;
	      if(z != 0.0)
		e->r[nn] = x-w/z;
	      e->i[nn] = e->i[nm] = 0.0;
	      x = a[nn][nm];
	      s = fabs(x)+fabs(z);
	      p = x/s;
	      q = z/s;
	      r = sqrt(SQUARE(p)+SQUARE(q));
	      p /= r;
	      q /= r;
	      for(j = nm; j <= n; j++)
		{
		  z = a[nm][j];
		  a[nm][j] = q*z+p*a[nn][j];
		  a[nn][j] = q*a[nn][j]-p*z;
		}
	      for(j = 1; j <= nn; j++)
		{
		  z = a[j][nm];
		  a[j][nm] = q*z+p*a[j][nn];
		  a[j][nn] = q*a[j][nn]-p*z;
		}
	      for(j = 1; j <= n; j++)
		{
		  z = e->u[j][nm];
		  e->u[j][nm] = q*z+p*e->u[j][nn];
		  e->u[j][nn] = q*e->u[j][nn]-p*z;
		}
	    } /* if (q>= 0.0) */
	  else
	    {
	      e->r[nn] = e->r[nm] = x+p;
	      e->i[nn] = -z;
	      e->i[nm] =  z;
	    }
	  nn = nl;
	} /* else if (l == nm) */
      else
	{
	  ni = 1;
	  if(its == 30)
	    {
	      puts(" E: hqr()-> Too many iterations.\n");
	      return;
	    }
	  if((its == 10) || (its == 20))
	    {
	      t += x;
	      for(i = 1; i <= nn; i++)
		a[i][i] -= x;
	      s = fabs(a[nn][nm])+fabs(a[nm][nl]);
	      x = 0.75*s;
	      y = x;
	      w = -0.4375*SQUARE(s);
	    }
	  its++;
	  for(m = nn-2; m > l-1; m--)
	    {
	      mp = m+1;
	      z = a[m][m];
	      r = x-z;
	      s = y-z;
	      p = (r*s-w)/a[mp][m]+a[m][mp];
	      q = a[mp][mp]-z-r-s;
	      r = a[m+2][mp];
	      s = fabs(p)+fabs(q)+fabs(r);
	      p /= s;
	      q /= s;
	      r /= s;
	      if (m == l)
		break;
	      t1 = fabs(p)*(fabs(a[m-1][m-1])+fabs(z)+fabs(a[mp][mp]));
	      t2 = t1+fabs(a[m][m-1])*(fabs(q)+fabs(r));
	      if (t2 == t1)
		break;
	    }
	  for(i = m+2; i <= nn; i++)
	    {
	      a[i][i-2] = 0.0;
	      if(i != (m+2))
		a[i][i-3] = 0.0;
	    }
	  for(k = m; k <= nm; k++)
	    {
	      if(k != m)
		{
		  p = a[k][k-1];
		  q = a[k+1][k-1];
		  r = 0.0;
		  if(k != nm)
		    r = a[k+2][k-1];
		  x = fabs(p)+fabs(q)+fabs(r);
		  if(x == 0.0)
		    continue;
		  p /= x;
		  q /= x;
		  r /= x;
		}
	      s = SGN(p) * sqrt(SQUARE(p)+SQUARE(q)+SQUARE(r));
	      if(k != m)
		a[k][k-1] = -s*x;
	      else if((k == m) && (l != m))
		a[k][k-1] = -a[k][k-1];
	      p += s;
	      x = p/s;
	      y = q/s;
	      z = r/s;
	      q /= p;
	      r /= p;
	      for(j = k; j <= n; j++)
		{
		  p = a[k][j]+q*a[k+1][j];
		  if (k != nm)
		    {
		      p += r*a[k+2][j];
		      a[k+2][j] -= p*z;
		    }
		  a[k+1][j] -= p*y;
		  a[k][j] -= p*x;
		}
	      mmin = MIN(nn, k+3);
	      for(j = 1; j <= mmin; j++)
		{
		  p = x*a[j][k]+y*a[j][k+1];
		  if (k != nm)
		    {
		      p += z*a[j][k+2];
		      a[j][k+2] -= p*r;
		    }
		  a[j][k+1] -= p*q;
		  a[j][k] -= p;
		}
	      for(j = 1; j <=n; j++)
		{
		  p = x*e->u[j][k]+y*e->u[j][k+1];
		  if (k != nm)
		    {
		      p += z*e->u[j][k+2];
		      e->u[j][k+2] -= p*r;
		    }
		  e->u[j][k+1] -= p*q;
		  e->u[j][k] -= p;
		}
	    } /* for(k = m; k <= nm; k++) */
	} /* else if (l != nn and l != nm) */
    } /* while(nn > 0) */

  /* Find eigenvectors
   */
  for(nn = n; nn > 0; nn--)
    {
      nm = nn-1;
      p = e->r[nn];
      q = e->i[nn];
      if(q < 0.0)
	{
	  m = nm;
	  if(fabs(a[nn][nm]) > fabs(a[nm][nn]))
	    {
	      a[nm][nm] = q/a[nn][nm];
	      a[nm][nn] = -(a[nn][nn]-p)/a[nn][nm];
	    }
	  else
	    {
	      cdiv(0.0, -a[nm][nn], a[nm][nm]-p, q, &w1, &w2);
	      a[nm][nm] = w1;
	      a[nm][nn] = w2;
	    }
	  a[nn][nm] = 0.0;
	  a[nn][nn] = 1.0;
	  nl = nm-1;
	  if(nl == 0)
	    continue;
	  for(i = nl; i > 0; i--)
	    {
	      w = a[i][i]-p;
	      ra = sa = 0.0;
	      for(j = m; j <= nn; j++)
		{
		  ra += a[i][j]*a[j][nm];
		  sa += a[i][j]*a[j][nn];
		}
	      if(e->i[i] < 0.0)
		{
		  z = w;
		  r = ra;
		  s = sa;
		  continue;
		}
	      m = i;
	      if(e->i[i] != 0.0)
		{
		  x = a[i][i+1];
		  y = a[i+1][i];
		  vr = SQUARE(e->r[i]-p)+SQUARE(e->i[i])-SQUARE(q);
		  vi = (e->r[i]-p)*2.0*q;
		  if((vr == 0.0) && (vi == 0.0))
		    {
		      t1 = anorm*(fabs(w)+fabs(q)+fabs(x)+fabs(y)+fabs(z));
		      vr = t1;
		      do {
			vr *= 0.01;
			t2 = t1+vr;
		      } while(t2 > t1);
		    }
		  cdiv(x*r-z*ra+q*sa, x*s-z*sa-q*ra, vr, vi, &w1, &w2);
		  a[i][nm] = w1;
		  a[i][nn] = w2;
		  if(fabs(x) > (fabs(z)+fabs(q)))
		    {
		      a[i+1][nm] = (-ra-w*a[i][nm]+q*a[i][nn])/x;
		      a[i+1][nn] = (-sa-w*a[i][nn]-q*a[i][nm])/x;
		    }
		  else
		    {
		      cdiv(-r-y*a[i][nm], -s-y*a[i][nn], z, q, &w1, &w2);
		      a[i+1][nm] = w1;
		      a[i+1][nn] = w2;
		    }
		} /* if(e->i[i] != 0.0) */
	      else
		{
		  cdiv(-ra, -sa, w, q, &w1, &w2);
		  a[i][nm] = w1;
		  a[i][nn] = w2;
		}
	      if(fabs(a[i][nm]) > fabs(a[i][nn]))
		t = fabs(a[i][nm]);
	      else
		t = fabs(a[i][nn]);
	      if(t == 0.0)
		continue;
	      t1 = t;
	      t2 = t1+1.0/t1;
	      if(t2 > t1)
		continue;
	      for(j = i; j <= nn; j++)
		{
		  a[j][nm] /= t;
		  a[j][nn] /= t;
		}
	    } /* for(i = nl; i > 0; i--) */
	} /* if(q < 0.0) */
      else if(q == 0.0)
	{
	  m = nn;
	  a[nn][nn] = 1.0;
	  if(nm == 0)
	    continue;
	  for(i = nm; i > 0; i--)
	    {
	      w = a[i][i]-p;
	      r = 0.0;
	      for(j = m; j <= nn; j++)
		r += a[i][j]*a[j][nn];
	      if(e->i[i] < 0.0)
		{
		  z = w;
		  s = r;
		  continue;
		}
	      m = i;
	      if(e->i[i] != 0.0)
		{
		  x = a[i][i+1];
		  y = a[i+1][i];
		  q = SQUARE(e->r[i]-p)+SQUARE(e->i[i]);
		  t = (x*s-z*r)/q;
		  a[i][nn] = t;
		  if(fabs(x) > fabs(z))
		    a[i+1][nn] = (-r-w*t)/x;
		  else
		    a[i+1][nn] = (-s-y*t)/z;
		}
	      else
		{
		  t = w;
		  if(t == 0.0)
		    {
		      t = t1 = anorm;
		      do {
			t *= 0.01;
			t2 = anorm+t;
		      } while(t2 > t1);
		    }
		  a[i][nn] = -r/t;
		}
	      t = fabs(a[i][nn]);
	      if(t == 0.0)
		continue;
	      t1 = t;
	      t2 = t1+1.0/t1;
	      if(t2 > t1)
		continue;
	      for(j = i; j <= nn; j++)
		a[j][nn] /= t;
	    } /* for(i = nm; i > 0; i--) */
	} /* else if(q == 0.0) */
    } /* for(nn = n; nn > 0; nn--) */
  for(j = n; j > 0; j--)
    {
      m = MIN(j, n);
      for(i = 1; i <= n; i++)
	{
	  z = 0.0;
	  for(k = 1; k <= m; k++)
	    z += e->u[i][k]*a[k][j];
	  e->u[i][j] = z;
	}
    }
}

/* Heapsort 
 *	
 * ic	Sort based on		
 * 0	decreasing modulus
 * 1	increasing residual
 * 2	decreasing real part
 * 3	decreasing imag part
 *
 * Originally written by Ian Gates, adapted from NR [Press et al
 * 1986], Nov 1994
 * 
 * I mangled the name for compatibility with Darwin OS, tab jan 2004
 */
void
_heapsort(int n,
	  dbl *ev_r,
	  dbl *ev_i,
	  dbl *ev_e,
	  int *order,
	  int  ic)
{
  int i, j, k, ir, lm, nn, 
    *iord;
  dbl xm, 
    *x;

  /* Initialize
   */
  nn = n+5;

  /* Allocate
   */

  iord = Ivector_birth(nn);
  x    = Dvector_birth(nn);

  /* Load
   */

  for(i = 0, k = 1; i < n; i++, k++)
    {
      order[k] = iord[k] = k;
      switch(ic)
	{
	case 0:
	  x[k] = sqrt(SQUARE(ev_r[i])+SQUARE(ev_i[i]));
	  break;
	case 1:
	  x[k] = ev_e[i];
	  break;
	case 2:
	  x[k] = ev_r[i];
	  break;
	case 3:
	  x[k] = ev_i[i];
	  break;
	default:
	  EH(GOMA_ERROR, "Uh-oh!  I shouldn't be here!");
	  break;
	}
    }
  k = n/2+1;
  ir = n;
  while(1) {
    if(k > 1)
      {
	k--;
	lm = order[k];
	xm = x[k];
      }
    else
      {
	lm = order[ir];
	xm = x[ir];
	x[ir] = x[1];
	order[ir] = order[1];
	ir--;
	if(ir == 1)
	  {
	    order[1] = lm;
	    x[1] = xm;
	    break;
	  }
      }
    i = k;
    j = k+k;
    while(j <= ir) {
      if(j < ir)
	if(x[j] < x[j+1])
	  j++;
      if(xm < x[j])
	{
	  order[i] = order[j];
	  x[i] = x[j];
	  i = j;
	  j = j+j;
	}
      else
	j = ir+1;
    }
    order[i] = lm;
    x[i] = xm;
  }

  /* Order */
  for(i = 0; i < n; i++)
    order[i] = iord[i] = order[i+1]-1;
  if(ic != 1)
    {
      for(i = 0; i < n; i++)
	x[n-1-i] = (dbl)(order[i]);
      for(i = 0; i < n; i++)
	order[i] = (int)(x[i]);
    }

  /* Real part eigenvalue */
  for(i = 0; i < n; i++)
    x[i] = ev_r[order[i]];
  for(i = 0; i < n; i++)
    ev_r[i] = x[i];
  /* Imag part eigenvalue */
  for(i = 0; i < n; i++)
    x[i] = ev_i[order[i]];
  for(i = 0; i < n; i++)
    ev_i[i] = x[i];

  /* Residual */
  for(i = 0; i < n; i++)
    x[i] = ev_e[order[i]];
  for(i = 0; i < n; i++)
    ev_e[i] = x[i];

  /* De-allocate
   */
  Ivector_death(&iord[0], nn);
  Dvector_death(&x[0],    nn);
}

/* Jeapsort 
 *
 * ic	sort based on		
 * 0	increasing modulus
 * 1	decreasing residual
 * 2	increasing real part
 * 3	increasing imag part
 *
 * Originally written by Ian Gates, adapted from NR [Press et al,
 * 1986], Nov 1994.
 */
void
jeapsort(int n,
	 dbl *ev_r,
	 dbl *ev_i,
	 dbl *ev_e,
	 int *order,
	 int ic)
{
  int i, j, k, ir, lm, nn, 
    *iord;
  dbl xm, 
    *x;

  /* Initialize
   */
  nn = n+5;

  /* Allocate
   */
  iord = Ivector_birth(nn);
  x    = Dvector_birth(nn);

  /* Load
   */
  for(i = 0, k = 1; i < n; i++, k++)
    {
      order[k] = iord[k] = k;
      switch(ic)
	{
	case 0:
	  x[k] = sqrt(SQUARE(ev_r[i])+SQUARE(ev_i[i]));
	  break;
	case 1:
	  x[k] = ev_e[i];
	  break;
	case 2:
	  x[k] = ev_r[i];
	  break;
	case 3:
	  x[k] = ev_i[i];
	  break;
	default:
	  EH(GOMA_ERROR, "Uh-oh!  I shouldn't be here!");
	  break;
	}
    }
  k = n/2+1;
  ir = n;
  while(1) {
    if(k > 1)
      {
	k--;
	lm = order[k];
	xm = x[k];
      }
    else
      {
	lm = order[ir];
	xm = x[ir];
	x[ir] = x[1];
	order[ir] = order[1];
	ir--;
	if(ir == 1)
	  {
	    order[1] = lm;
	    x[1] = xm;
	    break;
	  }
      }
    i = k;
    j = k+k;
    while(j <= ir) {
      if(j < ir)
	if(x[j] < x[j+1])
	  j++;
      if(xm < x[j])
	{
	  order[i] = order[j];
	  x[i] = x[j];
	  i = j;
	  j = j+j;
	}
      else
	j = ir+1;
    }
    order[i] = lm;
    x[i] = xm;
  }

  /* Order */
  for(i = 0; i < n; i++)
    order[i] = iord[i] = order[i+1]-1;
  if(ic == 1)
    {
      for(i = 0; i < n; i++)
	x[n-1-i] = (dbl)(order[i]);
      for(i = 0; i < n; i++)
	order[i] = (int)(x[i]);
    }

  /* Real part eigenvalue */
  for(i = 0; i < n; i++)
    x[i] = ev_r[order[i]];
  for(i = 0; i < n; i++)
    ev_r[i] = x[i];
  
  /* Imag part eigenvalue */
  for(i = 0; i < n; i++)
    x[i] = ev_i[order[i]];
  for(i = 0; i < n; i++)
    ev_i[i] = x[i];

  /* Residual */
  for(i = 0; i < n; i++)
    x[i] = ev_e[order[i]];
  for(i = 0; i < n; i++)
    ev_e[i] = x[i];

  /* De-allocate
   */
  Ivector_death(&iord[0], nn);
  Dvector_death(&x[0], nn);
}
