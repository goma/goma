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

#include "bc/eqn_index.h"
#include "bc/rotate_coordinates.h"
#include "mm_as.h"
#include "mm_bc.h"
#include "rf_bc.h"
#include "rf_node_const.h"

int bc_eqn_index(int id,          /* local node number                 */
                 int I,           /* processor node number             */
                 int bc_input_id, /* boundary condition number         */
                 int curr_matID,  /* Current material ID */
                 int kdir,        /* coordinate index for mesh or velocity
                                   * - normally this is zero unless a
                                   *   bc is a vector condition                 */
                 int *eqn,        /* eqn to which this condition is applied     */
                 int *matID_retn, /* material ID to apply this eqn on           */
                 VARIABLE_DESCRIPTION_STRUCT **vd_retn)

/*************************************************************************
 *
 * bc_eqn_index():
 *
 * Check to see if this BC on this node is applicable
 *     (i.e. no other overriding Dirichlet conditions),
 * Find the global unknown number for applying this condition.
 *
 * NOTE: This must be called after the element basis functions have
 *       been set up.
 *
 * Input
 * -------
 *
 * id = local node number
 * I  = processor node number
 * bc_input_id = boundary condition number
 * curr_matID = Current material index of the element that we
 *              calling from
 * BC_desc  = description structure for current boundary condition
 *
 *
 * Output
 * -------
 *  *eqn = Variable type of the equation to which this boundary condition
 *        is applied.
 *  *matID_retn = Material number of the equation to which this
 *                boundary condition is applied.
 *
 * Return
 * --------
 *   This function returns the index in the solution vector of the
 *   equation to which this boundary condition will be applied.
 *   If the boundary condition is not to be applied, this function
 *   return the value of -1.
 *
 *  *matID_retn = returns the material id corresponding to the unknown
 *                on which this boundary condition is applied. Even if
 *                the underyling variable type has a value of -1, this
 *                variable will not be -1.
 *************************************************************************/
{
  int ieqn, i_calc, jeqn, bc_node, offset_j;
  int matID, index_eqn, num, node_offset;
  int i, index;
  NODAL_VARS_STRUCT *nv;
  struct Boundary_Condition *bc = BC_Types + bc_input_id;
  struct BC_descriptions *bc_desc = bc->desc;
  NODE_INFO_STRUCT *node = Nodes[I];
  VARIABLE_DESCRIPTION_STRUCT *vd, *vd2 = NULL;
  nv = node->Nodal_Vars_Info[pg->imtrx];

  /*
   *  Find equation number and species number from the BC description
   *  structure
   */
  ieqn = bc_desc->equation;
  matID = curr_matID;

  /* need to change equation index for mesh or velocity */
  if (kdir != 0) {
    if (ieqn == R_MESH1)
      ieqn += kdir;
    else if (ieqn == R_MOMENTUM1)
      ieqn += kdir;
    else if (ieqn == R_SOLID1)
      ieqn += kdir;
    else if (ieqn == R_LAGR_MULT1)
      ieqn += kdir;
    else if (ieqn == R_EM_H1_REAL)
      ieqn += kdir; // AMC: These vector bcs are not rotated..
    else if (ieqn == R_EM_H1_IMAG)
      ieqn += kdir;
    else if (ieqn == R_EM_E1_REAL)
      ieqn += kdir;
    else if (ieqn == R_EM_E1_IMAG)
      ieqn += kdir;
  }

  /*
   *  Find out the offset under non-conflicting conditions
   *  If a boundary condition isn't to be applied under non-
   *  conflicting conditions, it won't be applied under
   *  conflicting conditions. So exit here if bc is not
   *  to be applied.
   */
  node_offset = find_bc_unk_offset(bc, curr_matID, I, kdir, &matID, &vd);
  if (node_offset < 0) {
    return -1;
  }

  /*
   * check the processor node number, I, against BC_dup_nodes[]
   * to see if this node is at an intersection of side-sets
   *
   * HKM -> the search below seems like a time-waste. In lieu of
   *        that, we could set up an added  bit field in
   *        Node_Info struct or elsewhere where we could store the
   *        results of this search (or just the need to do the
   *        search).
   */
  if (BC_dup_ptr != NULL) {
    bc_node = in_list(I, 0, BC_dup_ptr[pg->imtrx] + 1, BC_dup_nodes[pg->imtrx]);
    if (bc_node != -1) {
      /*
       * current node is in the BC_duplication list.
       * The current boundary condition may have been deleted.
       * If the current boundary condition is a rotated boundary
       * condition, it may have been moved to another coordinate
       * direction.
       */
      i_calc = -1;
      if (ieqn <= LAST_REAL_EQ) {
        i_calc = search_bc_dup_list(bc_input_id, BC_dup_list[pg->imtrx][bc_node][node_offset]);

      } else if ((ieqn >= R_MESH_NORMAL) && (ieqn <= R_MESH_TANG2)) {
        if (goma_automatic_rotations.automatic_rotations) {
          if (ieqn == R_MESH_NORMAL) {
            int best_dir = -1;
            double dot_max = 0;
            for (int dir = 0; dir < 3; dir++) {
              double dot = 0;
              for (int j = 0; j < 3; j++) {
                dot +=
                    goma_automatic_rotations.rotation_nodes[I].rotated_coord[dir]->normal->data[j] *
                    fv->snormal[j];
              }
              if (fabs(dot) > dot_max) {
                best_dir = dir;
                dot_max = fabs(dot);
              }
            }
            offset_j = get_nodal_unknown_offset(nv, R_MESH1 + best_dir, matID, 0, &vd2);
            if (offset_j >= 0) {
              node_offset = offset_j;
              ieqn = R_MESH1 + best_dir;
              vd = vd2;
              i_calc = 0;
            }
          } else {
            GOMA_EH(GOMA_ERROR, "Automatic rotations only for normal conditions");
          }
        } else {
          /*
           * If it's a rotated boundary condition, the boundary
           * condition could have been moved to another coordinate
           * within check_for_bc_conflicts2D(). Therefore, we need
           * to check all coordinate directions for the presence of the
           * current boundary condition.
           */
          for (jeqn = R_MESH1; jeqn <= R_MESH3; jeqn++) {
            offset_j = get_nodal_unknown_offset(nv, jeqn, matID, 0, &vd2);
            if (offset_j >= 0) {
              i_calc = search_bc_dup_list(bc_input_id, BC_dup_list[pg->imtrx][bc_node][offset_j]);
              if (i_calc != -1) {
                node_offset = offset_j;
                ieqn = jeqn;
                vd = vd2;
                break;
              }
            }
          }
        }
      } else if ((ieqn >= R_MOM_NORMAL) && (ieqn <= R_MOM_TANG2)) {
        if (goma_automatic_rotations.automatic_rotations) {
          if (ieqn == R_MOM_NORMAL) {
            int best_dir = -1;
            double dot_max = 0;
            for (int dir = 0; dir < 3; dir++) {
              double dot = 0;
              for (int j = 0; j < 3; j++) {
                dot +=
                    goma_automatic_rotations.rotation_nodes[I].rotated_coord[dir]->normal->data[j] *
                    fv->snormal[j];
              }
              if (fabs(dot) > dot_max) {
                best_dir = dir;
                dot_max = fabs(dot);
              }
            }
            offset_j = get_nodal_unknown_offset(nv, R_MOMENTUM1 + best_dir, matID, 0, &vd2);
            if (offset_j >= 0) {
              node_offset = offset_j;
              ieqn = R_MOMENTUM1 + best_dir;
              vd = vd2;
              i_calc = 0;
            }
          } else {
            GOMA_EH(GOMA_ERROR, "Automatic rotations only for normal conditions");
          }
        } else {
          for (jeqn = R_MOMENTUM1; jeqn <= R_MOMENTUM3; jeqn++) {
            offset_j = get_nodal_unknown_offset(nv, jeqn, matID, 0, &vd2);
            if (offset_j >= 0) {
              i_calc = search_bc_dup_list(bc_input_id, BC_dup_list[pg->imtrx][bc_node][offset_j]);
              if (i_calc != -1) {
                node_offset = offset_j;
                ieqn = jeqn;
                vd = vd2;
                break;
              }
            }
          }
        }
      } else if ((ieqn >= R_SOLID_NORMAL) && (ieqn <= R_SOLID_TANG2)) {
        for (jeqn = R_SOLID1; jeqn <= R_SOLID3; jeqn++) {
          offset_j = get_nodal_unknown_offset(nv, jeqn, matID, 0, &vd2);
          if (offset_j >= 0) {
            i_calc = search_bc_dup_list(bc_input_id, BC_dup_list[pg->imtrx][bc_node][offset_j]);
            if (i_calc != -1) {
              node_offset = offset_j;
              ieqn = jeqn;
              vd = vd2;
              break;
            }
          }
        }
      } else {
        GOMA_EH(GOMA_ERROR, "Equation not in list ");
      }

      /*
       * The current boundary condition has been deleted from application
       * at the current node due to a conflict with another boundary
       * condition. Return a -1.
       */
      if (i_calc == -1) {
        return -1;
      }
    }
  }
  /*
   * Not at an intersection:
   *    For boundary conditions which refer to rotated
   *    coordinates, assign an equation number to them
   *    base upon a (normal, tang1, tang2) = ( x, y, z)
   *    mapping.
   */
  if (goma_automatic_rotations.automatic_rotations) {
    if (ieqn == R_MESH_NORMAL || ieqn == R_MOM_NORMAL) {
      int best_dir = -1;
      double dot_max = 0;
      for (int dir = 0; dir < 3; dir++) {
        double dot = 0;
        for (int j = 0; j < 3; j++) {
          dot += goma_automatic_rotations.rotation_nodes[I].rotated_coord[dir]->normal->data[j] *
                 fv->snormal[j];
        }
        if (fabs(dot) > dot_max) {
          best_dir = dir;
          dot_max = fabs(dot);
        }
      }
      if (ieqn == R_MESH_NORMAL) {
        offset_j = get_nodal_unknown_offset(nv, R_MESH1 + best_dir, matID, 0, &vd2);
        if (offset_j >= 0) {
          node_offset = offset_j;
          ieqn = R_MESH1 + best_dir;
          vd = vd2;
          i_calc = 0;
        }
      } else if (ieqn == R_MOM_NORMAL) {
        offset_j = get_nodal_unknown_offset(nv, R_MESH1 + best_dir, matID, 0, &vd2);
        if (offset_j >= 0) {
          node_offset = offset_j;
          ieqn = R_MOMENTUM1 + best_dir;
          vd = vd2;
          i_calc = 0;
        }
      } else {
        GOMA_EH(GOMA_ERROR, "Not sure how we got here");
      }
    }
  }

  if ((ieqn >= R_MESH_NORMAL) && (ieqn <= R_MESH_TANG2)) {
    ieqn = ieqn - R_MESH_NORMAL + R_MESH1;
  }
  if ((ieqn >= R_SOLID_NORMAL) && (ieqn <= R_SOLID_TANG2)) {
    ieqn = ieqn - R_SOLID_NORMAL + R_SOLID1;
  }
  if ((ieqn >= R_MOM_NORMAL) && (ieqn <= R_MOM_TANG2)) {
    ieqn = ieqn - R_MOM_NORMAL + R_MOMENTUM1;
  }

  /*
   * If the regular unknown equation is replaced by a Dirichlet
   * condition at this node, return -1
   */
  if (node->DBC[pg->imtrx]) {
    if (((int)node->DBC[pg->imtrx][node_offset]) != -1) {
      return -1;
    }
  }

  /*
   * Set up the return variables
   */
  *matID_retn = matID;
  *eqn = ieqn;
  index_eqn = node->First_Unknown[pg->imtrx] + node_offset;

  /*
   * HKM -> we can put this section into
   *        a subroutine
   */
  if (vd_retn) {
    *vd_retn = NULL;
    num = nv->Num_Var_Desc_Per_Type[ieqn];
    for (i = 0; i < num; i++) {
      index = nv->Var_Type_Index[ieqn][i];
      vd = nv->Var_Desc_List[index];
      if (matID == vd->MatID || (!(*vd_retn) && vd->MatID == -1)) {
        *vd_retn = vd;
      }
    }
  }

  return index_eqn;
} /* END of routine bc_eqn_index                                             */

int bc_eqn_index_stress(int id,          /* local node number                 */
                        int I,           /* processor node number             */
                        int bc_input_id, /* boundary condition number         */
                        int curr_matID,  /* Current material ID */
                        int kdir,        /* coordinate index for stress components */
                        int mode,        /* Stress mode number */
                        int *eqn,        /* eqn to which this condition is applied     */
                        int *matID_retn, /* material ID to apply this eqn on           */
                        VARIABLE_DESCRIPTION_STRUCT **vd_retn)

/*************************************************************************
 *
 * bc_eqn_index_stress():
 *
 * Check to see if this BC on this node is applicable
 *     (i.e. no other overriding Dirichlet conditions),
 * Find the global unknown number for applying this condition. This is clone
 * of bc_eqn_index repurposed to handle stress equations
 *
 * NOTE: This must be called after the element basis functions have
 *       been set up.
 *
 * Input
 * -------
 *
 * id = local node number
 * I  = processor node number
 * bc_input_id = boundary condition number
 * curr_matID = Current material index of the element that we
 *              calling from
 * BC_desc  = description structure for current boundary condition
 *
 *
 * Output
 * -------
 *  *eqn = Variable type of the equation to which this boundary condition
 *        is applied.
 *  *matID_retn = Material number of the equation to which this
 *                boundary condition is applied.
 *
 * Return
 * --------
 *   This function returns the index in the solution vector of the
 *   equation to which this boundary condition will be applied.
 *   If the boundary condition is not to be applied, this function
 *   return the value of -1.
 *
 *  *matID_retn = returns the material id corresponding to the unknown
 *                on which this boundary condition is applied. Even if
 *                the underyling variable type has a value of -1, this
 *                variable will not be -1.
 *************************************************************************/
{
  int ieqn, i_calc, bc_node;
  int matID, index_eqn, num, node_offset;
  int i, index;
  int ndofs_ieqn;
  NODAL_VARS_STRUCT *nv;
  struct Boundary_Condition *bc = BC_Types + bc_input_id;
  struct BC_descriptions *bc_desc = bc->desc;
  NODE_INFO_STRUCT *node = Nodes[I];
  VARIABLE_DESCRIPTION_STRUCT *vd;
  nv = node->Nodal_Vars_Info[pg->imtrx];

  /*
   *  Find equation number from the BC description
   *  structure
   */
  ieqn = bc_desc->equation;
  /* Sanity check */
  if (ieqn != R_STRESS11)
    GOMA_EH(GOMA_ERROR, "You can't be here");

  matID = curr_matID;

  switch (mode) {

  case 0:
    ieqn = R_STRESS11 + kdir;
    break;

  case 1:
    ieqn = R_STRESS11_1 + kdir;
    break;

  case 2:
    ieqn = R_STRESS11_2 + kdir;
    break;

  case 3:
    ieqn = R_STRESS11_3 + kdir;
    break;

  case 4:
    ieqn = R_STRESS11_4 + kdir;
    break;

  case 5:
    ieqn = R_STRESS11_5 + kdir;
    break;

  case 6:
    ieqn = R_STRESS11_6 + kdir;
    break;

  case 7:
    ieqn = R_STRESS11_7 + kdir;
    break;

  default:
    GOMA_EH(GOMA_ERROR, "Maximum 8 modes allowed here!");
    break;
  }

  /*
   * If the equation has zero degrees of freedom at this node
   * return. Otherwise, get the number of variable description
   * structures with this variable type. Don't count more
   * than one species unknown per material at this point.
   */
  if ((ndofs_ieqn = get_nv_ndofs_modMF(nv, ieqn)) <= 0) {
    return -1;
  }

  /*
   *  Find out the offset under non-conflicting conditions
   *  If a boundary condition isn't to be applied under non-
   *  conflicting conditions, it won't be applied under
   *  conflicting conditions. So exit here if bc is not
   *  to be applied.
   */
  node_offset = get_nodal_unknown_offset(nv, ieqn, matID, 0, &vd);
  if (node_offset < 0) {
    return -1;
  }

  /*
   * check the processor node number, I, against BC_dup_nodes[]
   * to see if this node is at an intersection of side-sets
   *
   * HKM -> the search below seems like a time-waste. In lieu of
   *        that, we could set up an added  bit field in
   *        Node_Info struct or elsewhere where we could store the
   *        results of this search (or just the need to do the
   *        search).
   */
  bc_node = in_list(I, 0, BC_dup_ptr[pg->imtrx] + 1, BC_dup_nodes[pg->imtrx]);
  if (bc_node != -1) {
    /*
     * current node is in the BC_duplication list.
     * The current boundary condition may have been deleted.
     * If the current boundary condition is a rotated boundary
     * condition, it may have been moved to another coordinate
     * direction.
     */
    i_calc = -1;
    if (ieqn <= LAST_REAL_EQ) {
      i_calc = search_bc_dup_list(bc_input_id, BC_dup_list[pg->imtrx][bc_node][node_offset]);
    } else {
      GOMA_EH(GOMA_ERROR, "Equation not in list ");
    }

    /*
     * The current boundary condition has been deleted from application
     * at the current node due to a conflict with another boundary
     * condition. Return a -1.
     */
    if (i_calc == -1) {
      return -1;
    }
  }

  /*
   * If the regular unknown equation is replaced by a Dirichlet
   * condition at this node, return -1
   */
  if (node->DBC[pg->imtrx]) {
    if (((int)node->DBC[pg->imtrx][node_offset]) != -1) {
      return -1;
    }
  }

  /*
   * Set up the return variables
   */
  *matID_retn = matID;
  *eqn = ieqn;
  index_eqn = node->First_Unknown[pg->imtrx] + node_offset;

  /*
   * HKM -> we can put this section into
   *        a subroutine
   */
  if (vd_retn) {
    *vd_retn = NULL;
    num = nv->Num_Var_Desc_Per_Type[ieqn];
    for (i = 0; i < num; i++) {
      index = nv->Var_Type_Index[ieqn][i];
      vd = nv->Var_Desc_List[index];
      if (matID == vd->MatID || (!(*vd_retn) && vd->MatID == -1)) {
        *vd_retn = vd;
      }
    }
  }

  return index_eqn;
} /* END of routine bc_eqn_index_stress                                      */