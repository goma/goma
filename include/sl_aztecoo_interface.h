/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) Goma Developers, National Technology & Engineering        *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/

/*
 * sl_aztecoo_interface.h
 *
 *  Created on: Oct 27, 2014
 *      Author: wortiz
 */
#ifdef GOMA_ENABLE_AZTEC

#ifndef INCLUDE_SL_AZTECOO_INTERFACE_H_
#define INCLUDE_SL_AZTECOO_INTERFACE_H_

#ifdef __cplusplus
extern "C" {
#endif

void aztecoo_solve_epetra(struct GomaLinearSolverData *ams, double *x_, double *b_);

#ifdef __cplusplus
} /* End extern "C" */
#endif
#endif /* INCLUDE_SL_AZTECOO_INTERFACE_H_ */

#endif