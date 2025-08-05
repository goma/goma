/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2025 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/

#ifndef GOMA_BC_EQN_INDEX_H
#define GOMA_BC_EQN_INDEX_H

#include "rf_vars_const.h"

int bc_eqn_index(int id,          /* local node number                 */
                 int I,           /* processor node number             */
                 int bc_input_id, /* boundary condition number         */
                 int curr_matID,  /* Current material ID */
                 int kdir,        /* coordinate index for mesh or velocity
                                   * - normally this is zero unless a
                                   *   bc is a vector condition                 */
                 int *eqn,        /* eqn to which this condition is applied     */
                 int *matID_retn, /* material ID to apply this eqn on           */
                 VARIABLE_DESCRIPTION_STRUCT **vd_retn);

int bc_eqn_index_stress(int id,          /* local node number                 */
                        int I,           /* processor node number             */
                        int bc_input_id, /* boundary condition number         */
                        int curr_matID,  /* Current material ID */
                        int kdir,        /* coordinate index for stress components */
                        int mode,        /* Stress mode number */
                        int *eqn,        /* eqn to which this condition is applied     */
                        int *matID_retn, /* material ID to apply this eqn on           */
                        VARIABLE_DESCRIPTION_STRUCT **vd_retn);

#endif /* GOMA_BC_EQN_INDEX_H */