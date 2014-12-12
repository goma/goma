/*
 * sl_aztecoo_interface.h
 *
 *  Created on: Oct 27, 2014
 *      Author: wortiz
 */

#ifndef INCLUDE_SL_AZTECOO_INTERFACE_H_
#define INCLUDE_SL_AZTECOO_INTERFACE_H_

#ifdef __cplusplus
  extern "C" {
#endif

void
aztecoo_solve_epetra(struct Aztec_Linear_Solver_System *ams, double *x_, double *b_);

#ifdef __cplusplus
  } /* End extern "C" */
#endif
#endif /* INCLUDE_SL_AZTECOO_INTERFACE_H_ */
