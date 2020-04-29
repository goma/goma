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
 * $Id: sl_eggroll03.c,v 5.1 2007-09-18 18:53:47 prschun Exp $
 */

/*
 * $Log: not supported by cvs2svn $
 * Revision 5.0  2006/07/06 16:18:57  edwilke
 *
 * New Goma version: 'Farewell CRMPR, hello $3 a gallon gas!'
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
 * Revision 3.6  2001/02/08 15:36:59  mmhopki
 * - Main checkin for 3d stability of a 2d flow capability.  (I'll be
 * sending out an e-mail momentarily.)
 *
 *  - Also fixes UMF
 *
 * Revision 3.5  2000/05/18 05:41:17  dalabre
 * Multi-Platform changes, fixes and corrections (Part 1).
 * --------------------------------------------------------------------------
 *
 * Revision 3.4  2000/01/14 20:44:07  mmhopki
 * #include <stdio.h> required under Linux for FILE in a prototype in
 * another header file somewhere else.
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

#include "sl_auxutil.h"
#include "sl_eggroll.h"
#include "sl_umf.h"
#include "std.h"



/* Matrix-vector product for generalized eigenvalue problem
 *
 * MSR format only
 *
 * Friendly warning: do not edit this unless you know what you are
 * doing!!
 *
 * Originally written by Ian Gates
 */
void
gevp_transformation(int UMF_system_id,
		    int first,
		    int fflag,
                    int format,
                    int transformation, 
	            int nj, 
	            int nnz,
                    int *ija,
                    dbl *jac,
                    dbl *mas,
                    dbl *mat, 
			 /*	int soln_tech,  */
                    dbl *w,
                    dbl *v,
                    dbl r_sigma,
                    dbl i_sigma)
{
  dbl *z;

  /* Allocate work vectors
   */

  z = Dvector_birth(nj+5);

  /* z = M v 
   */
  MV_MSR(&nj, &ija[0], &mas[0], &v[0], &z[0]);

  /* Real shift matrix-vector product
   */
  UMF_system_id = SL_UMF(UMF_system_id,
			 &first, 
			 &fflag, 
			 &format, 
			 &nj, 
			 &nnz, 
			 ija, 
			 ija, 
			 mat, 
			 &z[0], 
			 &w[0]);

  /* De-allocate work storage
   */
  Dvector_death(&z[0], nj+5);
}








