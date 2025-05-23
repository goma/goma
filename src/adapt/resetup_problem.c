
#include "adapt/resetup_problem.h"

#include <mm_bc.h>
#include <rf_bc.h>
#include <rf_pre_proc.h>
#include <rf_solve.h>
#include <string.h>

#include "dp_map_comm_vec.h"
#include "dp_types.h"
#include "dpi.h"
#include "linalg/sparse_matrix.h"
#include "mm_as.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "mm_unknown_map.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_io.h"
#include "rf_masks.h"
#include "rf_mp.h"
#include "rf_node_const.h"
#include "rf_solver.h"
#include "rf_util.h"
#include "rf_vars_const.h"
#include "sl_util_structs.h"

int resetup_problem(Exo_DB *exo, /* ptr to the finite element mesh database */
                    Dpi *dpi)    /* distributed processing information */

/********************************************************************
 *
 * setup_problem():
 *
 *      Setup_problem() determines the degrees of freedom at each
 * node and formulates the solution vector. It then determines the
 * communications pattern for exchanging that solution vector between
 * processors.
 *     Lastly, it sets up structures that help to carry out the
 * boundary condition integrals on side sets.
 *
 * NOTE:
 *   This function was formed by taking common parts out of
 **********************************************************************/
{

  pre_process(exo);
  /*
   *  Initialize nodal based structures pertaining to properties
   *  and boundary conditions
   */
  init_nodes(exo, dpi);
  /*
   * Enumerate my own degrees of freedom on this processor.
   */
  setup_local_nodal_vars(exo, dpi);

  /*
   * setup communications patterns between ghost and owned
   * nodes.
   */
  setup_nodal_comm_map(exo, dpi, cx);

  /*
   * Find the global maximum number of unknowns located at any one
   * node on any processor
   */
  MaxVarPerNode = find_MaxUnknownNode();

  /*
   * Exchange my idea of what materials I have at each node with
   * my surrounding processors. Make sure we are all in sync
   */
  setup_external_nodal_matrls(exo, dpi, cx[0]);

  /*
   * Exchange my idea of what degrees of freedom I have with my
   * surrounding processors. Make sure we are all in sync.
   * Owned nodes tell ghost nodes what variables are active
   * at that node.
   */
  setup_external_nodal_vars(exo, dpi, cx);

  /*
   * Finish setting the unknown map on this processor
   */
  set_unknown_map(exo, dpi);

  /*
   * Now determine the communications pattern that is necessary to
   * exchange the solution vector in as efficient a manner as
   * possible
   */
  // log_msg("setup_dof_comm_map...");
  setup_dof_comm_map(exo, dpi, cx);

  /*
   * Output some statistics concerning the communications pattern
   */
  if (Num_Proc > 1)
    output_comm_stats(dpi, cx);

  /*
   * I extracted this from setup_fill_comm_map because some of the
   * renormalization routines make use of them.
   */
  num_fill_unknowns = count_vardofs(FILL, dpi->num_universe_nodes);
  internal_fill_unknowns = count_vardofs(FILL, dpi->num_internal_nodes);
  owned_fill_unknowns = count_vardofs(FILL, (dpi->num_internal_nodes + dpi->num_boundary_nodes));
  boundary_fill_unknowns = owned_fill_unknowns - internal_fill_unknowns;
  external_fill_unknowns = num_fill_unknowns - owned_fill_unknowns;

  /*
   *  Possibly increase the number of variable descriptions to include
   *  those that are not part of the solution vector
   */
  vdesc_augment();

  /*
   *  Setup Boundary condition inter-connectivity
   *  -> Stefan flow bc's need to know what bc's contain the reaction
   *     info.
   *  -> YFLUX bc's need to be connected -> future implementation
   */
  set_up_BC_connectivity();

  /*
   *  Set up the structures necessary to carry out
   *  surface integrals
   */
  //  log_msg("set_up_Surf_BC...");
  set_up_Surf_BC(First_Elem_Side_BC_Array, exo, dpi);

  /*
   *  Set up the Edge boundary condition structures
   */
  //  log_msg("set_up_Edge_BC...");
  set_up_Edge_BC(First_Elem_Edge_BC_Array, exo, dpi);

  /*
   * Set up "boundary" conditions on level set surfaces
   */
  //  set_up_Embedded_BC();

  /* Special 1D
   * Set up surface integral boundary conditions
   * that apply at single nodes.
   */

  //  setup_Point_BC(First_Elem_Side_BC_Array, exo, dpi);

  /*
   *  Print out the edge boudary condition structures
   *  if necessary
   */
  //  if (Debug_Flag) {
  //    for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
  //      print_setup_Surf_BC(First_Elem_Side_BC_Array[pg->imtrx]);
  //    }
  //  }

  /*
   *  Malloc structures of size Num_Var_Info_Records
   */
  //  for (int imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
  //    ei[pg->imtrx]->VDindex_to_Lvdesc = alloc_int_1(Num_Var_Info_Records, -1);
  //  }
  //  for (int i = 0; i < MAX_ELEMENT_INDICES_RELATED; i++) {
  //    eiRelated[i]->VDindex_to_Lvdesc = alloc_int_1(Num_Var_Info_Records, -1);
  //  }

  /*
   * Setup the storage for temporary quantities of interest that
   * are storred on a "per processor element" basis. These include
   * temporary storage of volumetric quadrature information
   */
  //  setup_element_storage();

  /*
   * Setup some structures for solving problems with shell elements.
   */
  //  init_shell_element_blocks(exo);

  /* Communicate non-shared but needed BC information */
  //  exchange_bc_info();

  return 0;
}

int resetup_matrix(struct GomaLinearSolverData **ams, Exo_DB *exo, Dpi *dpi) {
  if ((strcmp(Matrix_Format, "tpetra") == 0) || (strcmp(Matrix_Format, "epetra") == 0)) {
    for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
      GomaSparseMatrix goma_matrix = ams[pg->imtrx]->GomaMatrixData;
      GomaSparseMatrix_Destroy(&goma_matrix);
      GomaSparseMatrix_CreateFromFormat(&goma_matrix, Matrix_Format);
      ams[pg->imtrx]->GomaMatrixData = goma_matrix;
      int local_nodes = exo->num_nodes;
      GomaSparseMatrix_SetProblemGraph(goma_matrix, num_internal_dofs[pg->imtrx],
                                       num_boundary_dofs[pg->imtrx], num_external_dofs[pg->imtrx],
                                       local_nodes, Nodes, MaxVarPerNode, Matilda, Inter_Mask, exo,
                                       dpi, cx[pg->imtrx], pg->imtrx, Debug_Flag, ams[JAC]);
    }
    pg->imtrx = 0;
    ams[pg->imtrx]->solveSetup = 0;
  } else {
    GOMA_EH(-1, "Unsupported matrix storage format use epetra");
  }
  return 0;
}

// vim: expandtab sw=2 ts=8
