/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2022 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/

/*
 * $Id: sl_eggrollutil.c,v 5.1 2007-09-18 18:53:47 prschun Exp $
 */

#include "sl_eggroll.h"
#include "std.h"

/* Utility routines for sl_eggroll*.c.
 *
 * Originall written by Ian Gates.
 *
 * Hacked by MMH to conform with Goma style.
 */

/* Complex division
 */
void cdiv(dbl ar, dbl ai, dbl br, dbl bi, dbl *cr, dbl *ci) {
  dbl s;
  s = SQUARE(br) + SQUARE(bi);
  (*cr) = (ar * br + ai * bi) / s;
  (*ci) = (ai * br - ar * bi) / s;
}
