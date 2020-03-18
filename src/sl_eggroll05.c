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
 * $Id: sl_eggroll05.c,v 5.1 2007-09-18 18:53:47 prschun Exp $
 */

#include "sl_auxutil.h"
#include "sl_eggroll.h"
#include "std.h"
#ifdef USE_RCSID
static const char rcs_id[] = "$Id: sl_eggroll05.c,v 5.1 2007-09-18 18:53:47 prschun Exp $";
#endif



/* Re-order eigenvalues and eigenvectors by real part
 *
 * Friendly warning: do not edit this unless you know what you are
 * doing!!
 *
 * Originally written by Ian Gates.
 *
 */
void
gevp_order(int nj,
	   int ev_n,
	   dbl *ev_r,
	   dbl *ev_i,
	   dbl *ev_e,
	   dbl *ev_x,
	   dbl **evect,
	   dbl **schur)
{
  int i, k, nn, 
    *order;
  dbl **wk;

  /* Initialize
   */
  nn = ev_n+5;

  /* Allocate
   */
  order = Ivector_birth(nn);
  for(i = 0; i < nn; i++) 
    order[i] = i;
  wk = Dmatrix_birth(ev_n, nj);
  
  /* Reorder eigenvalues - increasing modulus
   */
  jeapsort(ev_n, &ev_r[0], &ev_i[0], &ev_e[0], &order[0], 0);

  /* Reorder eigenvectors
   */

  for(i = 0; i < ev_n; i++)
    {
      k = order[i];
      vcopy(nj, &wk[i][0], 1.0, &evect[k][0]);
    }
  for(i = 0; i < ev_n; i++)
    vcopy(nj, &evect[i][0], 1.0, &wk[i][0]);
  for(i = 0; i < ev_n; i++)
    {
      k = order[i];
      vcopy(nj, &wk[i][0], 1.0, &schur[k][0]);
    }
  for(i = 0; i < ev_n; i++)
    vcopy(nj, &schur[i][0], 1.0, &wk[i][0]);

  /* De-allocate temp space
   */
  Ivector_death(&order[0], nn);
  Dmatrix_death(wk, ev_n, nj);
}



