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
 *$Id: mm_shell_util.c,v 5.17 2010-07-21 16:39:27 hkmoffa Exp $
 */

#ifdef USE_RCSID
static char rcsid[] = "$Id: mm_shell_util.c,v 5.17 2010-07-21 16:39:27 hkmoffa Exp $";
#endif

/* Standard include files */
#define GOMA_MM_SHELL_UTIL_C
#include "mm_shell_util.h"
#include "az_aztec.h"
#include "density.h"
#include "el_elm.h"
#include "el_geom.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_shell.h"
#include "mm_fill_terms.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "mm_std_models_shell.h"
#include "mm_viscosity.h"
#include "rd_mesh.h"
#include "rf_allo.h"
#include "rf_bc.h"
#include "rf_bc_const.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_fill_const.h"
#include "rf_io.h"
#include "rf_io_const.h"
#include "rf_masks.h"
#include "rf_mp.h"
#include "rf_solver.h"
#include "rf_solver_const.h"
#include "rf_vars_const.h"
#include "shell_tfmp_util.h"
#include "sl_util.h"
#include "std.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <limits.h>

#define DEBUG_SHELL 0
#define PRINTPROC   0
#define DSPRINTF                          \
  if (DEBUG_SHELL && ProcID == PRINTPROC) \
  printf
#define D2SPRINTF                             \
  if (DEBUG_SHELL > 1 && ProcID == PRINTPROC) \
  printf

/*********** R O U T I N E S   I N   T H I S   F I L E ***********************
 *
 *						-All routines in this
 *				NOTE:		 file are called only by
 *						 mm_fill.c: matrix_fill
 *						 except possibly for static
 *						 functions.
 *
 *       NAME			TYPE			CALLED BY
 *  -----------------------------------------------------------------
 *
 *  init_shell_element_blocks()  void
 *  is_shell_element_type()      int
 *  is_shell_element()           int
 *  is_shell_block()             int
 *  find_stu_from_xyz()          int
 *  solve_2x2()                  int
 *  solve_3x3()                  int
 *  bulk_side_id_and_stu()       int
 *  load_neighbor_var_data()     int
 *  find_stu_on_shell()          int
 *  find_stu_on_bulk()           int
 *  shell_normal_div_s()         int
 *  shell_determinant_and_normal() int
 *  shell_tangents               void
 *  shell_tangents_isoparametric void
 *  shell_tangents_seeded        void
 *  calculate_lub_q_v()          void
 *  calculate_lub_q_v_old()      void
 *  shell_saturation_pressure_curve() dbl
 *
 */

/******************************************************************************
 *
 * init_shell_element_blocks()
 *
 * Synopsis
 * ========
 * This function first checks for the exists of blocks of shell elements.
 * If such blocks are found, this function creates and populates an array
 * of Shell_Block structs.  This function also creates the element arrays
 * mum_elem_friends[elem] and elem_friends[elem][friend_elem].
 *
 * Input
 * =====
 * exo = Pointer to the exodus database struct
 *
 * Output
 * ======
 * None. Creates and populates global 'shell_blocks' if necessary.
 *
 * Returns
 * =======
 * void
 *
 * Revision History
 * ================
 * 7 May 2002 - Patrick Notz - Creation
 *
 ******************************************************************************/
int **elem_friends = NULL;
int *num_elem_friends = NULL;
int num_shell_blocks = 0;

void init_shell_element_blocks(const Exo_DB *exo) {
  int i, elem, n, nbr;
  int bindex, be, bn, bnn, bnpe, bnel, bid, bnoff;
  int sindex, se, sn, snn, snpe, snel, sid, snoff, si;
  int sn_found, all_sn_found;
  int snodes[9];
  int *nbr_elem_ids;
  int *bconn;
  int *friends_count;
  int num_found;
  int bID, mnParent, numMnFound;
  int candidate, id;
  int *mnList = alloc_int_1(MAX_NUMBER_MATLS, -1);
  struct Shell_Block *sb;
  const int num_elems = exo->num_elems;
  /* Global Variables: num_shell_blocks, shell_blocks, bulk_shell_nbrs */

  DSPRINTF("Initializing shell-element handling: ");

  /* Loop over the element blocks and look for shell blocks */
  num_shell_blocks = 0;
  for (i = 0; i < exo->num_elem_blocks; ++i) {
    num_shell_blocks += is_shell_element_type(exo->eb_elem_itype[i]);
  }
  DSPRINTF("%d shell blocks found.\n", num_shell_blocks);

  int *eb_ProcessingOrder = alloc_int_1(exo->num_elem_blocks, -1);
  int *eb_Used = alloc_int_1(exo->num_elem_blocks, 0);

  /*
   *  In order to get a consistent element ordering of friends for mp problems,
   *  we order the element block search algorithm based on the element block
   *  ids instead of the order of the element blocks in the exo object. This has been
   *  shown to fix the surface normal calculation in shells. The reason for this
   *  is that only the first element friend of a shell element is used. That had
   *  better be the correct (and consistent) side of the shell element.
   *
   *  A more consistent treatment of this issue would be to choose the shell side
   *  which is to become the first element friend of shell elements, and move that
   *  choice to the input file. However, this has not been done yet.
   */
  for (i = 0; i < exo->num_elem_blocks; ++i) {
    si = INT_MAX;
    candidate = -1;
    for (bindex = 0; bindex < exo->num_elem_blocks; ++bindex) {
      id = exo->eb_id[bindex];
      if (id < si) {
        if (eb_Used[bindex] == 0) {
          si = id;
          candidate = bindex;
        }
      }
    }
    eb_ProcessingOrder[i] = candidate;
    eb_Used[candidate] = 1;
  }

  /* No shells? Do nothing */
  if (num_shell_blocks == 0) {
    safe_free(mnList);
    safe_free(eb_Used);
    safe_free(eb_ProcessingOrder);
    return;
  }

  /* Allocate some space for the list of Shell_Block structures. */
  /* Also do some basic initialization. */
  shell_blocks = (struct Shell_Block **)smalloc(num_shell_blocks * sizeof(struct Shell_Block *));
  for (bindex = 0, si = -1; bindex < exo->num_elem_blocks; ++bindex) {
    if (is_shell_element_type(exo->eb_elem_itype[bindex])) {
      ++si;
      shell_blocks[si] = (struct Shell_Block *)smalloc(sizeof(struct Shell_Block));
      shell_blocks[si]->elemblock_index = bindex;
      shell_blocks[si]->elemblock_id = exo->eb_id[bindex];
      shell_blocks[si]->num_nbr_blocks = 0;
      shell_blocks[si]->mn = Matilda[bindex];
    }
  }

  /*
   * Find all of the neighboring blocks for each shell block.
   *
   * NB: A shell face is a curve (3D) or a vertex (2D) so you can just
   * compare the faces of neigboring elements.  For example, in 3D, a
   * shell's neighbor has a face which is essentially the whole shell
   * element.  Thus, we're looking for a match between N+1 and N dimensional
   * mesh entities.
   *
   * For each shell element Es with nodes Ns find the list of
   * elements Eb which contain ALL of Ns.
   */
  for (si = 0; si < num_shell_blocks; ++si) {
    sb = shell_blocks[si]; /* A pointer to the current shell block */
    sindex = sb->elemblock_index;
    sid = sb->elemblock_id;
    snel = exo->eb_num_elems[sindex];
    snpe = exo->eb_num_nodes_per_elem[sindex];
    if (snpe > 9)
      GOMA_EH(GOMA_ERROR, "Blecch! Too many nodes in this shell element. Muy malo.\n");

    DSPRINTF("Processing shell block: si=%d, sindex=%d, sid=%d (snel=%d)\n", si, sindex, sid, snel);

    /* Check this shell block against every other block */
    for (candidate = 0; candidate < exo->num_elem_blocks; ++candidate) {
      bindex = eb_ProcessingOrder[candidate];
      /* Don't check against self. */
      if (bindex == sindex)
        continue;

      nbr_elem_ids = (int *)smalloc(snel * sizeof(int));
      bid = exo->eb_id[bindex];
      bnpe = exo->eb_num_nodes_per_elem[bindex];
      bnel = exo->eb_num_elems[bindex];
      bconn = exo->eb_conn[bindex];

      DSPRINTF("Checking against block bid=%d\n", bid);

      /* Try to find a neighbor for each shell element */
      /* NB: we should find 0 or snel matches */
      num_found = -1;
      for (se = 0; se < snel; ++se) {
        DSPRINTF("Searching for neighbors of the shell element %d\n", exo->eb_ptr[sindex] + se + 1);

        /* Form the list of nodes to search for */
        snoff = se * snpe;
        DSPRINTF("\tSearching for nodes: ");
        for (sn = 0; sn < snpe; ++sn) {
          snodes[sn] = exo->eb_conn[sindex][snoff + sn];
          DSPRINTF("%d ", snodes[sn] + 1);
        }
        DSPRINTF("\n");

        /* OK, check against each element in the "other" block */
        for (be = 0; be < bnel; ++be) {
          D2SPRINTF("\tChecking against element %d\n", exo->eb_ptr[bindex] + be + 1);
          /* see if this element has all of the shell nodes */
          all_sn_found = 1;
          for (sn = 0; sn < snpe; ++sn) {
            snn = snodes[sn];
            sn_found = 0;
            bnoff = be * bnpe;
            for (bn = 0; bn < bnpe; ++bn) {
              bnn = bconn[bnoff + bn];
              D2SPRINTF("\t\tsnn=%d, bnn=%d\n", snn + 1, bnn + 1);
              sn_found = bnn == snn;
              if (sn_found)
                break;
            }
            if (sn_found)
              continue;

            /* We'll only get this far if a node wasn't found. */
            all_sn_found = 0;
            break;
          }
          if (all_sn_found) {
            /* Hey, be is a neighbor of se! */
            /* save this element's info */
            elem = exo->eb_ptr[bindex] + be;
            nbr_elem_ids[se] = elem;
            ++num_found;
#if DEBUG_SHELL
            DSPRINTF("\tElement %d is a neighbor, nodes: ", exo->eb_ptr[bindex] + be + 1);
            for (bn = 0; bn < bnpe; ++bn) {
              bnn = bconn[bnoff + bn];
              DSPRINTF("%d", bnn + 1);
              if (bn < bnpe - 1) {
                DSPRINTF(", ");
              }
            }
            DSPRINTF("\n");
#endif
            /*
             * For now, to make the bookkeeping easier,
             * I won't allow for looped meshes (one
             * shell element has two neighbors from the
             * same element block) so I'll break out...
             */
            break;
          }
        } /* be loop: for each element in the other block */

        if (num_found == -1) {
          /*
           * OK, we couldn't find a match for the first 'se' so this
           * must not be a neighboring block.
           */
          break;
        }
        /* Sanity checking */
        /* EDW: This was causing trouble for problems with partially wetted blocks.
           else if(num_found != se)
           {
           fprintf(stderr,"\n\n");
           fprintf(stderr,"What the shell is going on???\n");
           fprintf(stderr,"Shell Element ID  = %d\n",se);
           fprintf(stderr,"Shell Block ID    = %d\n",sid);
           fprintf(stderr,"Shell Block Index = %d\n",sindex);
           fprintf(stderr,"Other Block ID    = %d\n",bid);
           fprintf(stderr,"Other Block Index = %d\n",bindex);
           GOMA_EH(GOMA_ERROR,"Ack. It appears as though this shell's neighbor moved out.\n");
           }
        */

      } /* se loop: for each element in the shell block */

      /* At this point num_found is either 0 or snel */
      if (num_found == (snel - 1)) {
        n = sb->num_nbr_blocks;
        if (n >= MAX_SHELL_NBRS)
          GOMA_EH(GOMA_ERROR, "Ack. MAX_SHELL_NBRS isn't big enough.\n");
        sb->nbr_elem_ids[n] = nbr_elem_ids;
        sb->num_nbr_blocks += 1;
#if DEBUG_SHELL
        DSPRINTF("shell element: neighbor\n");
        for (se = 0; se < snel; ++se) {
          DSPRINTF("%d: %d\n", exo->eb_ptr[sb->elemblock_index] + se + 1, nbr_elem_ids[se] + 1);
        }
#endif
      } else {
        safe_free((void *)nbr_elem_ids);
      }

    } /* bindex loop: other element blocks */

  } /* si loop: shell blocks */

  /*
   * OK, we now have, for each element in each block of shells, a list of
   * that elements neighbors.  These elements will need to share
   * degrees of freedom (DoFs) even though each elements' DoFs may not
   * exist in the others.  We'll call such elements "friends" since we'll
   * have to bloat their lec->, esp->, fv-> etc. structures to make a
   * meta-element so we can handle sensitivies properly.
   *
   * So, next, for each shell element we'll make a list of that element's friends.
   */
  DSPRINTF("\n ***** CONSTRUCTING THE ELEMENT FRIENDS LIST ***** \n");

  /* Allocate some space for the elem_friends arrays. */
  num_elem_friends = alloc_int_1(num_elems, 0);
  friends_count = alloc_int_1(num_elems, 0);
  elem_friends = (int **)alloc_ptr_1(num_elems);

  /* Populate the num_elem_friends[elem] array. */
  DSPRINTF(" -> Generting the num_elem_friends[elem] array...\n");
  for (si = 0; si < num_shell_blocks; ++si) {
    sb = shell_blocks[si]; /* A pointer to the current shell block */
    sindex = sb->elemblock_index;
    sid = sb->elemblock_id;
    snel = exo->eb_num_elems[sindex];

    DSPRINTF("     + Processing shell block %d\n", sid);

    for (bindex = 0; bindex < sb->num_nbr_blocks; ++bindex) {

      DSPRINTF("     + Processing (block) neighbor %d\n", bindex);
      for (se = 0; se < snel; ++se) {
        elem = exo->eb_ptr[sindex] + se;
        nbr = sb->nbr_elem_ids[bindex][se];
        DSPRINTF("     + elem=%d, nbr=%d\n", elem + 1, nbr + 1);

        /* This element's neighbor counts as a friend */
        num_elem_friends[elem]++;

        /* Likewise, this element is a friend of its neighbor. */
        /* However, if the neighbor is a shell, we must not double count. */
        if (!is_shell_element(nbr, exo))
          num_elem_friends[nbr]++;
      }
    }
  }

  /* Now we know how many neighbors each element has so we can alloc */
  for (elem = 0; elem < num_elems; ++elem) {
    /* Only worry about elems with neighbors */
    if (num_elem_friends[elem] == 0)
      continue;

    n = num_elem_friends[elem];

    /* Alloc the space */
    elem_friends[elem] = (int *)smalloc(n * sizeof(int));

    /* Initialize the space */
    for (i = 0; i < n; ++i)
      elem_friends[elem][i] = -1;
  }

  /* Populate the elem_friends[elem][friend] array. */
  DSPRINTF(" -> Generting the elem_frieds[elem][nbr] array...\n");
  for (si = 0; si < num_shell_blocks; ++si) {
    sb = shell_blocks[si]; /* A pointer to the current shell block */
    sindex = sb->elemblock_index;
    sid = sb->elemblock_id;
    snel = exo->eb_num_elems[sindex];

    DSPRINTF("     + Processing shell block %d\n", sid);

    for (bindex = 0; bindex < sb->num_nbr_blocks; ++bindex) {
      DSPRINTF("     + Processing (block) neighbor %d\n", bindex);
      for (se = 0; se < snel; ++se) {
        elem = exo->eb_ptr[sindex] + se;
        nbr = sb->nbr_elem_ids[bindex][se];
        DSPRINTF("     + elem=%d, nbr=%d", elem + 1, nbr + 1);

        /* This shell's neighbor */
        n = friends_count[elem];
        if (n >= num_elem_friends[elem])
          GOMA_EH(GOMA_ERROR, "Ack. I have more neighbors than I thought!\n");
        elem_friends[elem][n] = nbr;
        friends_count[elem]++;
        DSPRINTF(" (elem updated)");

        /* Reverse */
        if (!is_shell_element(nbr, exo)) {
          DSPRINTF(" (nbr updated too)");
          n = friends_count[nbr];
          if (n >= num_elem_friends[nbr])
            GOMA_EH(GOMA_ERROR, "Ack. I have more neighbors than I thought!\n");
          elem_friends[nbr][n] = elem;
          friends_count[nbr]++;
        }
        DSPRINTF("\n");
      } /* se */
    }   /* bindex */
  }     /* si */

  safe_free((void *)friends_count);

#if DEBUG_SHELL
  DSPRINTF("\n ***** THE ELEMENT-FRIENDS RELATIONSHIPS ***** \n");
  for (elem = 0; elem < exo->num_elems; ++elem) {
    DSPRINTF("\tElement %5d has %d friends: ", elem + 1, num_elem_friends[elem]);
    for (n = 0; n < num_elem_friends[elem]; ++n)
      DSPRINTF("%5d ", elem_friends[elem][n] + 1);
    DSPRINTF("\n");
  }

#endif

  /*
   *  SECTION TO HANDLE Shell Element Block on Parent Element Block Logic
   */
  for (si = 0; si < num_shell_blocks; ++si) {
    sb = shell_blocks[si]; /* A pointer to the current shell block */
    sindex = sb->elemblock_index;
    sid = sb->elemblock_id;
    snel = exo->eb_num_elems[sindex];
    numMnFound = 0;

    DSPRINTF("     + Processing shell block %d\n", sid);
    for (bindex = 0; bindex < sb->num_nbr_blocks; ++bindex) {
      for (se = 0; se < snel; ++se) {
        elem = exo->eb_ptr[sindex] + se;
        nbr = sb->nbr_elem_ids[bindex][se];
        bID = find_elemblock_index(nbr, exo);
        mnParent = Matilda[bID];
        int notFound = 1;
        for (i = 0; i < numMnFound; i++) {
          if (mnParent == mnList[i]) {
            notFound = 0;
          }
        }
        if (notFound) {
          mnList[numMnFound] = mnParent;
          numMnFound++;
        }
      }
    }
    /*
     * We now have a list of parent materials adjacent to
     * this element block.
     */
    pd = pd_glob[sb->mn];
    int *pdv = pd->v[pg->imtrx];
    for (i = 0; i < numMnFound; i++) {
      mnParent = mnList[i];
      PROBLEM_DESCRIPTION_STRUCT *PDParent = pd_glob[mnParent];
      int *pdvParent = PDParent->v[pg->imtrx];

      /*
       *   We will now specify that if mesh equations aren't turned on
       *   in the shell, but they are turned on in the parent element,
       *   then the variables are flagged as V_SPECIFIED in the
       *   shell equations.
       *
       *   Note, the intension in some of the documentation is to do
       *   this for all variables. However, let's start with the significant
       *   case of mesh displacements and see where we are.
       *
       *   We should also note that pd->e[pg->imtrx][MESH_DISPLACEMENT1] is not
       *   turned on in this treatment.
       */
      if (pdv[MESH_DISPLACEMENT1] == 0) {
        if (pdvParent[MESH_DISPLACEMENT1] && V_SOLNVECTOR) {
          pdv[MESH_DISPLACEMENT1] = V_SPECIFIED;
          pd->i[pg->imtrx][MESH_DISPLACEMENT1] = PDParent->i[pg->imtrx][MESH_DISPLACEMENT1];
          /*
           * Ok, we turned on the interpolants for the mesh equation
           * We will need to adjust the ShapeVar and the InterpolantType
           * in order to set them back to the mesh equation
           */
          determine_ShapeVar(pd);
          determine_ProjectionVar(pd);
        }
      }
      if (pdv[MESH_DISPLACEMENT2] == 0) {
        if (pdvParent[MESH_DISPLACEMENT2] && V_SOLNVECTOR) {
          pdv[MESH_DISPLACEMENT2] = V_SPECIFIED;
          pd->i[pg->imtrx][MESH_DISPLACEMENT2] = PDParent->i[pg->imtrx][MESH_DISPLACEMENT2];
        }
      }
      if (pdv[MESH_DISPLACEMENT3] == 0) {
        if (pdvParent[MESH_DISPLACEMENT3] && V_SOLNVECTOR) {
          pdv[MESH_DISPLACEMENT3] = V_SPECIFIED;
          pd->i[pg->imtrx][MESH_DISPLACEMENT3] = PDParent->i[pg->imtrx][MESH_DISPLACEMENT3];
        }
      }
    }
  }

  safe_free(mnList);
  safe_free(eb_Used);
  safe_free(eb_ProcessingOrder);

  return;

} /* init_shell_element_blocks() */

/******************************************************************************
 * Check if a given element type is a type of shell element
 ******************************************************************************/
int is_shell_element_type(const int elem_itype) {
  /* These two element types apply to 2D problems and now 3D problems */
  return (elem_itype == LINEAR_BAR || elem_itype == QUAD_BAR || elem_itype == BILINEAR_SHELL ||
          elem_itype == BIQUAD_SHELL || elem_itype == BILINEAR_TRISHELL || elem_itype == P0_SHELL ||
          elem_itype == P1_SHELL);
}

/******************************************************************************
 * Check if a given element is a shell element
 ******************************************************************************/
int is_shell_element(const int elem, const Exo_DB *exo) {
  const int elem_itype = exo->eb_elem_itype[exo->elem_eb[elem]];
  return is_shell_element_type(elem_itype);
}

/******************************************************************************
 * Check if a given element block is a block of shell elements
 ******************************************************************************/
int is_shell_block(const int block_id, const Exo_DB *exo) {
  int eb_index;

  /* Find the block index */
  for (eb_index = 0; eb_index < exo->num_elem_blocks; eb_index++)
    if (exo->eb_id[eb_index] == block_id)
      break;

  /* Sanity check */
  if (eb_index == exo->num_elem_blocks) {
    fprintf(stderr, "While looking for block_id %d\n", block_id);
    GOMA_EH(GOMA_ERROR, "Invalid block id");
  }

  return is_shell_element_type(exo->eb_elem_itype[eb_index]);
}

/******************************************************************************
 *
 * find_stu_from_xyz()
 *
 * Synopsis
 * ========
 * Given coordinates (x,y,z) find the element coordinates (s,t,u) within a
 * given element.  This is solved using Newton's method as follows:
 *
 * xi    = [s,t,u], Elemental coordinates
 * xg    = [x,y,z], Global coordinates
 * X_j   = [x_j, y_j, z_j], Interpolation coefficients for the mesh at node j
 * phi_j = j'th basis function
 *
 * xg = SUM_j X_j phi_j
 *
 * Find xi using Newton's method:
 *
 * R = SUM_j{ X_j phi_j } - xg
 *
 * J = [ @R/@s  @R/@t  @R/@u ] = [ SUM_j{ X_j @phi_j/@s } ... ]
 *
 * J d = -R
 *
 * etc...
 *
 * Input
 * =====
 * elem = Element number to search in.
 * xg[] = [x,y,z], Global coordinates of the point in question.
 * xi[] = [s,t,u], Initial guess for [s,t,u].
 *
 * Output
 * ======
 * xi[] = [s,t,u], Element coordinates corresponding to [x,y,z].
 *
 * Returns
 * =======
 * -1 = No converged solution
 *  0 = Success, (s,t,u) c [-1,1]
 *  1 = Converged, but at least one of (s,t,u) are outide of [-1,1]
 *
 * Revision History
 * ================
 * 8 May 2002 - Patrick Notz - Creation
 *
 ******************************************************************************/
int find_stu_from_xyz(const int elem, const double xg[DIM], double xi[DIM], const Exo_DB *exo) {
  double x[DIM];
  double dxstu[DIM][DIM];
  double dxi[DIM];
  double R[DIM];
  double J[DIM][DIM];
  double dnorm, rnorm;
  int i, j, k, node, index, dofs, iter, var, dof_list_index;
  int sane;
  int dim;
  const double nm_tol = 1.0e-6;
  const int max_iter = 5;
  const int el1 = ei[pg->imtrx]->ielem;

  dim = ei[pg->imtrx]->ielem_dim;

  DSPRINTF("\t        dim = %d\n", dim);
  DSPRINTF("\t        DIM = %d\n", DIM);
  DSPRINTF("\tpd->Num_Dim = %d\n", pd->Num_Dim);

  /* Make the element being searched in the current element */
  load_ei(elem, exo, 0, pg->imtrx);

  /* Top of search iteration loop */
  for (iter = 0; iter < max_iter; iter++) {
    for (i = 0; i < DIM; i++) {
      dxi[i] = 0.0;
      x[i] = 0.0;
      R[i] = 0.0;
      for (j = 0; j < DIM; j++) {
        dxstu[i][j] = 0.0;
        J[i][j] = 0.0;
      }
    }

    DSPRINTF("\tIter %d of %d max iterations:\n", iter, max_iter);
    /* Compute the basis function at xi[] */

    DSPRINTF("\tload basis functions\n");
    if (xi == NULL)
      printf(" LOADING BF WITH NULL xi!\n");
    if (bfd == NULL)
      printf(" LOADING BF WITH NULL bfd!\n");
    load_basis_functions(xi, bfd);

    DSPRINTF("\tget derivs for j\n");
    DSPRINTF("\t\txi  = (%g, %g, %g)\n", xi[0], xi[1], xi[2]);
    /*
     * At this xi[], find x[] and it's derivs w.r.t. xi[] using
     * the finite element basis functions
     */
    for (j = 0; j < pd->Num_Dim; j++) {
      var = MESH_DISPLACEMENT1 + j;
      if (pd->v[pg->imtrx][var]) {
        dof_list_index = R_MESH1;
      } else {
        var = pd->ShapeVar;
        dof_list_index = var;
      }
      dofs = ei[pg->imtrx]->dof[var];
      ;
      for (i = 0; i < dofs; i++) {
        node = ei[pg->imtrx]->dof_list[dof_list_index][i];
        index = Proc_Elem_Connect[Proc_Connect_Ptr[elem] + node];
        x[j] += Coor[j][index] * bf[var]->phi[i];
        for (k = 0; k < pd->Num_Dim; k++) {
          dxstu[j][k] += Coor[j][index] * bf[var]->dphidxi[i][k];
        }
      }
    }

    DSPRINTF("\t\tassemble r and j\n");
    /* Assemble R & J */
    for (j = 0; j < pd->Num_Dim; j++) {
      R[j] = x[j] - xg[j];
      for (k = 0; k < pd->Num_Dim; k++)
        J[j][k] = dxstu[j][k];
    }

    DSPRINTF("R = \n\t[");
    for (j = 0; j < pd->Num_Dim; j++)
      DSPRINTF(" %g ", R[j]);
    DSPRINTF("]\n");

    DSPRINTF("J = \n");
    for (j = 0; j < pd->Num_Dim; j++) {
      DSPRINTF("\t[");
      for (k = 0; k < pd->Num_Dim; k++)
        DSPRINTF(" %g ", J[j][k]);
      DSPRINTF("]\n");
    }

    DSPRINTF("\t\tsolve for dxi\n");
    /* Solve for dxi */
    if (pd->Num_Dim == 2)
      solve_2x2(J, R, dxi);
    else
      solve_3x3(J, R, dxi);
    DSPRINTF("\t\tdxi = (%g, %g, %g)\n", dxi[0], dxi[1], dxi[2]);

    DSPRINTF("\t\tupdate xi\n");
    /* Update xi */
    for (i = 0; i < pd->Num_Dim; i++) {
      xi[i] -= dxi[i];
    }

    /* Compute norms */
    rnorm = 0.0;
    dnorm = 0.0;
    for (i = 0; i < pd->Num_Dim; i++) {
      dnorm += dxi[i] * dxi[i];
      rnorm += R[i] * R[i];
    }
    dnorm = sqrt(dnorm);
    DSPRINTF("\t\tcheck convergence\n");
    /* Check for convergence */
    if (dnorm < nm_tol && rnorm < nm_tol)
      break;
  }
  /* Bottom of search iteration loop */

  if (iter == max_iter) {
    /* didn't converge */
    GOMA_EH(GOMA_ERROR, "Didn't converge on a solution.");
    return -1;
  }

  /* Check for a sensible solution */
  sane = 1;
  for (i = 0; i < dim; i++) {
    if (xi[i] < -1.0 || xi[i] > 1.0)
      sane = 0;
  }

  /* Restore original current element */
  load_ei(el1, exo, 0, pg->imtrx);

  if (sane)
    return 0;
  else
    return 1;
}

/*
 * solve_2x2 - analytical solution of a system of two equations
 * NOTE: The matrix/vectors have dimension DIM
 */
int solve_2x2(double a[DIM][DIM], double b[DIM], double x[DIM]) {
  double det, deti;

  det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
  if (fabs(det) < 1.0e-12)
    return -1;

  deti = 1.0 / det;

  x[0] = deti * (a[1][1] * b[0] - a[0][1] * b[1]);
  x[1] = deti * (-a[1][0] * b[0] + a[0][0] * b[1]);

  return 0;
}

/*
 * solve_3x3 - analytical solution of a system of three equations
 * NOTE: The matrix/vectors have dimension DIM
 */
int solve_3x3(double a[DIM][DIM], double b[DIM], double x[DIM]) {
  double det, deti;

  det = +a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
        a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
        a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);

  if (fabs(det) < 1.0e-12)
    return -1;

  deti = 1.0 / det;

  /*
   * The solution is the transpose of the conjugate matrix
   * divided by the determinant.
   */

  x[0] = +(a[1][1] * a[2][2] - a[1][2] * a[2][1]) * b[0] -
         (a[0][1] * a[2][2] - a[0][2] * a[2][1]) * b[1] +
         (a[0][1] * a[1][2] - a[0][2] * a[1][1]) * b[2];

  x[1] = -(a[1][0] * a[2][2] - a[1][2] * a[2][0]) * b[0] +
         (a[0][0] * a[2][2] - a[0][2] * a[2][0]) * b[1] -
         (a[0][0] * a[1][2] - a[0][2] * a[1][0]) * b[2];

  x[2] = +(a[1][0] * a[2][1] - a[1][1] * a[2][0]) * b[0] -
         (a[0][0] * a[2][1] - a[0][1] * a[2][0]) * b[1] +
         (a[0][0] * a[1][1] - a[0][1] * a[1][0]) * b[2];

  x[0] *= deti;
  x[1] *= deti;
  x[2] *= deti;

  return 0;
}

int bulk_side_id_and_stu(const int bulk_elem,
                         const int shell_elem,
                         const double xi[DIM],
                         double xi2[DIM],
                         const Exo_DB *exo)
/*
 * Determines which bulk element nodes are shared with the shell element,
 * then deduces the corresponding bulk side ID and converts the shell
 * stu coordinates to those for the same point in the bulk element.
 *
 * Output
 * -------------
 *      xi2[DIM]  stu coordinates in the underlying bulk element
 *      return    side id of the bulk element
 *
 * PRS: as of 11-30-2011 I certify to the best of my ability that this routine
 * is now correct.     Famous last words.  See me if you want to have my model to
 * check it yourself.
 */
{
  int id = -1;
  int i, n, n1, n2, nno;
  double sign, Hsign = 0.0, Lsign = 0.0, r;
  int bulk_dim, shell_dim, nbulk = 0, nshell = 0, nmatch;
  int iH = -1, iL = -1, s_first = 0, bulk_n1, bulk_n2;
  int bulk_type = Elem_Type(exo, bulk_elem);
  int bulk_eptr = Proc_Connect_Ptr[bulk_elem];
  int shell_type = Elem_Type(exo, shell_elem);
  int shell_eptr = Proc_Connect_Ptr[shell_elem];

  /* Initialize */
  nno = 0;

  /* Get element dimensions */
  bulk_dim = elem_info(NDIM, bulk_type);
  shell_dim = elem_info(NDIM, shell_type);

  /* Limited to 2D bulk and 1D shell or 3D bulk and 1D shell for now */
  if (bulk_dim == 3 && shell_dim == 2) {
    switch (shell_type) {
    case BILINEAR_SHELL:
      nbulk = 8;
      nshell = 4;
      break;
    case BILINEAR_TRISHELL:
      nbulk = 4;
      nshell = 3;
      break;
    case BIQUAD_SHELL:
      nbulk = 27;
      nshell = 9;
      break;
    default:
      GOMA_EH(GOMA_ERROR, "Shell element type not supported.");
      break;
    }
  } else if (bulk_dim == 2 && shell_dim == 1) {
    switch (shell_type) {
    case LINEAR_BAR:
      nbulk = 4;
      nshell = 2;
      break;
    case QUAD_BAR:
      nbulk = 9;
      nshell = 3;
      break;
    default:
      GOMA_EH(GOMA_ERROR, "Shell element type not supported.");
      break;
    }
  } else
    GOMA_EH(GOMA_ERROR, "Only 2D bulk / 1D shell and 3D bulk / 2D shell are allowed!");

  int bulk_nodes[nbulk], shell_nodes[nshell], matches[nbulk];
  int bulk_node_order[nshell];

  /* Create node lists */
  for (i = 0; i < nbulk; i++) {
    matches[i] = -1;
    bulk_nodes[i] = Proc_Elem_Connect[bulk_eptr + i];
  }
  for (i = 0; i < nshell; i++) {
    shell_nodes[i] = Proc_Elem_Connect[shell_eptr + i];
  }

  /* Find indices of bulk nodes shared with the shell element */
  /* n1 and n2 keep track of the first and last shared nodes */
  nmatch = 0;
  n1 = 9999;
  n2 = -9999;
  for (i = 0; i < nshell; i++) {
    n = in_list(shell_nodes[i], 0, nbulk, &bulk_nodes[0]);
    if (n >= 0) {
      nmatch++;
      matches[i] = n;
      if (n > n2)
        n2 = n;
      if (n < n1)
        n1 = n;
    }
  }

  /* For tet / tri elements, let's also find which global node */
  /* does not lie on the shell (helps identify the bulk side). */
  if (shell_type == BILINEAR_TRISHELL) {
    for (i = 0; i < nbulk; i++) {
      n = in_list(bulk_nodes[i], 0, nshell, &shell_nodes[0]);
      if (n == -1) {
        nno = i;
        break;
      }
    }
  }

  /* Consistency checks */
  if ((nmatch != nshell) || (n1 >= n2)) {
    GOMA_EH(GOMA_ERROR, "Failed bulk/shell node match consistency check!");
  }

  /*
   * First, handle the 2D/1D case (this is the easier one by far).
   *
   * Fall through the local node indices for each possible bulk side:
   * Get id_side from the first and last bulk local node indices matched,
   * Use local node ordering to determine if the relevant xi coordinate
   * runs in the same or opposite direction in the two elements,
   * then use the sign to convert the running xi coordinate to the bulk
   * while assigning the fixed coordinate to +1 or -1, depending on id_side.
   * This is some "bulk"y logic, but it will get the job done!
   */
  if (bulk_dim == 2 && shell_dim == 1) {
    if (shell_type == LINEAR_BAR) {
      if (n1 == 0 && n2 == 1) {
        id = 1;
        if (matches[0] == 1 && matches[1] == 0) {
          sign = -1.0;
        } else {
          sign = 1.0;
        }
        xi2[0] = sign * xi[0];
        xi2[1] = -1.0;
      } else if (n1 == 1 && n2 == 2) {
        id = 2;
        if (matches[0] == 2 && matches[1] == 1) {
          sign = -1.0;
        } else {
          sign = 1.0;
        }
        xi2[0] = 1.0;
        xi2[1] = sign * xi[0];
      } else if (n1 == 2 && n2 == 3) {
        id = 3;
        if (matches[0] == 2 && matches[1] == 3) {
          sign = -1.0;
        } else {
          sign = 1.0;
        }
        xi2[0] = sign * xi[0];
        xi2[1] = 1.0;
      } else if (n1 == 0 && n2 == 3) {
        id = 4;
        if (matches[0] == 3 && matches[1] == 0) {
          sign = -1.0;
        } else {
          sign = 1.0;
        }
        xi2[0] = -1.0;
        xi2[1] = sign * xi[0];
      }
    } else if (shell_type == QUAD_BAR) {
      if (n1 == 0 && n2 == 4) {
        id = 1;
        if (matches[0] == 1 && matches[1] == 0) {
          sign = -1.0;
        } else {
          sign = 1.0;
        }
        xi2[0] = sign * xi[0];
        xi2[1] = -1.0;
      } else if (n1 == 1 && n2 == 5) {
        id = 2;
        if (matches[0] == 2 && matches[1] == 1) {
          sign = -1.0;
        } else {
          sign = 1.0;
        }
        xi2[0] = 1.0;
        xi2[1] = sign * xi[0];
      } else if (n1 == 2 && n2 == 6) {
        id = 3;
        if (matches[0] == 2 && matches[1] == 3) {
          sign = -1.0;
        } else {
          sign = 1.0;
        }
        xi2[0] = sign * xi[0];
        xi2[1] = 1.0;
      } else if (n1 == 0 && n2 == 7) {
        id = 4;
        if (matches[0] == 3 && matches[1] == 0) {
          sign = -1.0;
        } else {
          sign = 1.0;
        }
        xi2[0] = -1.0;
        xi2[1] = sign * xi[0];
      }
    }
  }

  /*
   * Handle the 3D/2D case here as follows:
   *
   * a) Determine the side ID from the high and low local bulk
   *    coordinates (n1, n2) within a case block.
   * b) Determine the indices of the two running xi coordinates (iH, iL),
   *    then set the constant cordinate to +1 or -1 according to side ID.
   * c) Set the base orientation for shared nodes with respect to the
   *    bulk element by placing the local node ID's in order.
   * d) The signs (+1.0, -1.0) for coordinate conversion are determined
   *    by which local bulk node matches the first local shell node.
   * e) The bulk and shell xi coordinates are then matched by determining
   *    the direction from the first local shell element to the second.
   * f) The two running bulk xi coordinates are then filled in using the
   *    results of steps d and e.
   *
   * For this case, there are 48 possible relative orientations of the
   * two elements - this is an attempt to handle all of them as concisely
   * as possible.
   */
  else if ((bulk_dim == 3) && (shell_dim == 2) &&
           ((shell_type == BILINEAR_SHELL) || (shell_type == BIQUAD_SHELL))) {

    /* Quickly determine the side ID from n1 and n2 */

    if (shell_type == BILINEAR_SHELL) {
      if (n1 == 0 && n2 == 5) {
        id = 1;
      } else if (n1 == 1 && n2 == 6) {
        id = 2;
      } else if (n1 == 2 && n2 == 7) {
        id = 3;
      } else if (n1 == 0 && n2 == 7) {
        id = 4;
      } else if (n1 == 0 && n2 == 3) {
        id = 5;
      } else if (n1 == 4 && n2 == 7) {
        id = 6;
      } else {
        GOMA_EH(GOMA_ERROR, "Couldn't determine 3D element side ID for HEX8!");
      }
    }

    else if (shell_type == BIQUAD_SHELL) {
      if (n1 == 0 && n2 == 25) {
        id = 1;
      } else if (n1 == 1 && n2 == 24) {
        id = 2;
      } else if (n1 == 2 && n2 == 26) {
        id = 3;
      } else if (n1 == 0 && n2 == 23) {
        id = 4;
      } else if (n1 == 0 && n2 == 21) {
        id = 5;
      } else if (n1 == 4 && n2 == 22) {
        id = 6;
      } else {
        GOMA_EH(GOMA_ERROR, "Couldn't determine 3D element side ID for HEX27!");
      }
    }

    /* Quickly determine the side ID from n1 and n2 */

    /* Set parameters according to side ID */
    switch (id) {
    case 1:
      xi2[1] = -1.0;
      iH = 0;
      iL = 2;
      bulk_node_order[0] = 0;
      bulk_node_order[1] = 1;
      bulk_node_order[2] = 4;
      bulk_node_order[3] = 5;
      break;
    case 2:
      xi2[0] = 1.0;
      iH = 1;
      iL = 2;
      bulk_node_order[0] = 1;
      bulk_node_order[1] = 2;
      bulk_node_order[2] = 5;
      bulk_node_order[3] = 6;
      break;
    case 3:
      xi2[1] = 1.0;
      iH = 0;
      iL = 2;
      bulk_node_order[0] = 3;
      bulk_node_order[1] = 2;
      bulk_node_order[2] = 7;
      bulk_node_order[3] = 6;
      break;
    case 4:
      xi2[0] = -1.0;
      iH = 1;
      iL = 2;
      bulk_node_order[0] = 0;
      bulk_node_order[1] = 3;
      bulk_node_order[2] = 4;
      bulk_node_order[3] = 7;
      break;
    case 5:
      xi2[2] = -1.0;
      iH = 0;
      iL = 1;
      bulk_node_order[0] = 0;
      bulk_node_order[1] = 1;
      bulk_node_order[2] = 3;
      bulk_node_order[3] = 2;
      break;
    case 6:
      xi2[2] = 1.0;
      iH = 0;
      iL = 1;
      bulk_node_order[0] = 4;
      bulk_node_order[1] = 5;
      bulk_node_order[2] = 7;
      bulk_node_order[3] = 6;
      break;
    }

    /* Find ordered bulk nodes which match the first two local shell nodes */
    bulk_n1 = in_list(matches[0], 0, 4, &bulk_node_order[0]);
    GOMA_EH(bulk_n1, "Node matching inconsistency (shell->bulk)!");
    bulk_n2 = in_list(matches[1], 0, 4, &bulk_node_order[0]);
    GOMA_EH(bulk_n2, "Node matching inconsistency (shell->bulk)!");

    /* Assign conversion signs and coordinate pairing */
    switch (bulk_n1) {
    case 0:
      Hsign = 1.0;
      Lsign = 1.0;
      if (bulk_n2 == 1) {
        s_first = TRUE;
      } else if (bulk_n2 == 2) {
        s_first = FALSE;
      } else {
        GOMA_EH(GOMA_ERROR, "SNAFU in bulk_side_id_and_stu!");
      }
      break;
    case 1:
      Hsign = -1.0;
      Lsign = 1.0;
      if (bulk_n2 == 0) {
        s_first = TRUE;
      } else if (bulk_n2 == 3) {
        s_first = FALSE;
      } else {
        GOMA_EH(GOMA_ERROR, "SNAFU in bulk_side_id_and_stu!");
      }
      break;
    case 2:
      Hsign = 1.0;
      Lsign = -1.0;
      if (bulk_n2 == 3) {
        s_first = TRUE;
      } else if (bulk_n2 == 0) {
        s_first = FALSE;
      } else {
        GOMA_EH(GOMA_ERROR, "SNAFU in bulk_side_id_and_stu!");
      }
      break;
    case 3:
      Hsign = -1.0;
      Lsign = -1.0;
      if (bulk_n2 == 2) {
        s_first = TRUE;
      } else if (bulk_n2 == 1) {
        s_first = FALSE;
      } else {
        GOMA_EH(GOMA_ERROR, "SNAFU in bulk_side_id_and_stu!");
      }
      break;
    }

    /* Now set the two running coordinates of xi2 for bulk element */
    if (s_first) {
      xi2[iH] = Hsign * xi[0];
      xi2[iL] = Lsign * xi[1];
    } else {
      xi2[iH] = Hsign * xi[1];
      xi2[iL] = Lsign * xi[0];
    }
  }

  /* This is the case for triangular shell / tetrahedral bulk elements */
  else if (shell_type == BILINEAR_TRISHELL) {

    /* Find bulk side ID */
    switch (nno) {
    case 0:
      id = 2;
      break;
    case 1:
      id = 3;
      break;
    case 2:
      id = 1;
      break;
    case 3:
      id = 4;
      break;
    }

    /* Find bulk node index for shell nodes 1 and 2 */
    n1 = in_list(shell_nodes[1], 0, nbulk, &bulk_nodes[0]);
    n2 = in_list(shell_nodes[2], 0, nbulk, &bulk_nodes[0]);

    /* Calculate the bulk coordinate that is 0.0 (if any) from excluded node */
    if (nno > 0)
      xi2[nno - 1] = 0.0;

    /* Calculate bulk coordinates coorresponding to xi1[0] and xi1[1] */
    if (n1 > 0)
      xi2[n1 - 1] = xi[0];
    if (n2 > 0)
      xi2[n2 - 1] = xi[1];

    /* Calculate which bulk coordinate must be r, if any */
    r = 1.0 - xi[0] - xi[1];
    for (i = 0; i < nshell; i++) {
      if (((n1 - 1) != i) && ((n2 - 1) != i) && ((nno - 1) != i))
        xi2[i] = r;
    }
  }

  /* If ID not found, error (-1) is returned */
  return id;
}

/****************************************************************************************
 * load_neighbor_var_data:
 *
 * Synopsis
 * ========
 * This routine will fill up values in fv and bf structures from the parent element for child
 *elements (or vica-versa).
 *
 * This is called for the local element, el1. The local
 * element may be either the bulk element or the shell element. The Remote element
 * is then the opposite.
 *
 * Ouptut
 * =========
 *  ndofs[]   =  Number of degrees of freedom for each variable in el2
 *  dof_map[i] = Mapping between the local node id and the local node id
 *               of the remote (e.g., parent) element. For example if there are 9 local nodes
 *               in the parent, and 3 in the shell, then this array will be a vector of
 *               length 3, with the values being the three nodes in the parent that
 *               correspond to the local node ids of the shell local nodes.
 *                 eg:   dof_map[0] = 3
 *                       dof_map[1] = 1
 *                       dof_map[2] = 7
 *               The mapping occurs between the shape variables for the shell and the
 *               parent elements.
 *  ndofptr[var,j] Map from (iVartype, j) to the local element degree of freedom in the remote
 *               element.
 *
 *
 * When el1 is the bulk, ddim = 1.
 * When el1 is the shell, ddim = -1.
 *
 ****************************************************************************************/
int load_neighbor_var_data(int el1,    // Element number of the local element
                           int el2,    // element number of the remote element
                           int *ndofs, // Number of degrees of freedom for each variable in el2
                           int *dof_map,
                           int ndofptr[MAX_VARIABLE_TYPES][MDE], // map from (iVartype, j) to local
                                                                 // element degree of freedom
                           const int id,
                           double xi[DIM],
                           const Exo_DB *exo) {
  double xi2[DIM];
  int i, j, node;
  int dim1, dim2, ddim;
  int nodes_per_side;
  int id_side;
  int local_elem_node_id[MAX_NODES_PER_SIDE];
  int svar, svar1, svar2, bulk_gnn[MDE], shell_gnn[MDE];

  /* Load coordinates from current element and initialize xi2 */
  for (i = 0; i < DIM; i++) {
    xi2[i] = 0.0;
  }

  /* Determine element dimensions */
  dim1 = elem_info(NDIM, Elem_Type(exo, el1)); /* Local element */
  dim2 = elem_info(NDIM, Elem_Type(exo, el2)); /* Remote element */
  ddim = dim1 - dim2;

  /*
   * Get stu coordinates for this Gauss point from neighbor element
   * Two cases to handle here:
   *
   * If bulk element is local and shell element is remote, i.e.
   * assembling a bulk element BC using shell variables,
   * then elem_side_bc provides the side id, which can be used
   * to convert the bulk stu coordinates to the shell stu
   * coordinates without iteration. This should work for
   * either 2D or 3D problems.
   */
  if (ddim == 1 && id != -1) {
    find_stu_on_shell(el1, id, el2, dim1, xi[0], xi[1], xi[2], xi2, exo);
    id_side = id;
  }

  /* Otherwise, find remote stu by iteration*/
  else {
    id_side = bulk_side_id_and_stu(el2, el1, xi, xi2, exo);
    if (id_side == -1) {
      printf("Error returned from bulk_side_id_and_stu()\n");
      GOMA_EH(GOMA_ERROR, "load_neighbor_var_data ");
    }
  }

  /* This call will populate the neighbor element fv/bf structures */
  setup_shop_at_point(el2, xi2, exo);

  /* svar2 is the shape variable for the second (i.e., parent) element */
  svar2 = pd->ShapeVar;

  /* Save active DOF counts for variables in the neighbor element */
  for (i = 0; i < MAX_VARIABLE_TYPES; i++) {
    ndofs[i] = ei[pg->imtrx]->dof[i];
    for (j = 0; j < ndofs[i]; j++) {
      ndofptr[i][j] = ei[pg->imtrx]->ieqn_ledof[ei[pg->imtrx]->lvdof_to_ledof[i][j]];
    }
  }

  if (ddim == -1) {
    /* compute surface normal and derivs on bulk elem */
    get_side_info(ei[pg->imtrx]->ielem_type, id_side, &nodes_per_side, local_elem_node_id);
    surface_determinant_and_normal(el2, ei[pg->imtrx]->iconnect_ptr, ei[pg->imtrx]->num_local_nodes,
                                   ei[pg->imtrx]->ielem_dim - 1, id_side, nodes_per_side,
                                   local_elem_node_id);
    if (ei[pg->imtrx]->ielem_dim != 3) {
      calc_surf_tangent(el2, ei[pg->imtrx]->iconnect_ptr, ei[pg->imtrx]->num_local_nodes,
                        ei[pg->imtrx]->ielem_dim - 1, nodes_per_side, local_elem_node_id);
    }
    /* get info for node mapping if needed */
    if (dof_map != NULL) {
      svar = pd->ShapeVar;
      for (i = 0; i < MDE; i++)
        dof_map[i] = -1;
      for (i = 0; i < ei[pg->imtrx]->dof[svar]; i++) {
        bulk_gnn[i] = ei[pg->imtrx]->gnn_list[svar][i];
      }
    }
  }
  if (ddim == 1) {
    /* get info for node mapping if needed */
    if (dof_map != NULL) {
      svar = pd->ShapeVar;
      for (i = 0; i < MDE; i++)
        dof_map[i] = -1;
      for (i = 0; i < ei[pg->imtrx]->dof[svar]; i++) {
        shell_gnn[i] = ei[pg->imtrx]->gnn_list[svar][i];
      }
    }
  }

  /* This call will repopulate the local element fv/bf structures */
  if (ddim == -1) {
    InShellElementWithParentElementCoverage = 1;
    (void)copy_int_1(ShellElementParentElementCoverageForVariable, pd->v, MAX_VARIABLE_TYPES);
  }
  /* Go back to the current element (i.e., the shell element) */
  setup_shop_at_point(el1, xi, exo);

  if (ddim == -1) {
    InShellElementWithParentElementCoverage = 0;
    zero_int_1(ShellElementParentElementCoverageForVariable, MAX_VARIABLE_TYPES);
  }

  svar1 = pd->ShapeVar;
  /* Construct node map if needed */
  if (ddim == -1 && dof_map != NULL) {
    for (i = 0; i < ei[pg->imtrx]->dof[svar1]; i++) {
      node = ei[pg->imtrx]->gnn_list[svar1][i];
      for (j = 0; j < ndofs[svar2]; j++) {
        if (node == bulk_gnn[j])
          dof_map[i] = j;
      }
    }
  }
  if (ddim == 1 && dof_map != NULL) {
    for (i = 0; i < ei[pg->imtrx]->dof[svar1]; i++) {
      node = ei[pg->imtrx]->gnn_list[svar1][i];
      for (j = 0; j < ndofs[svar2]; j++) {
        if (node == shell_gnn[j])
          dof_map[j] = i;
      }
    }
  }
  /* Done */
  return 0;
}

int find_stu_on_shell(const int bulk_elem,
                      const int id_side,
                      const int shell_elem,
                      const int bulk_dim,
                      const double s,
                      const double t,
                      const double u,
                      double xi2[DIM],
                      const Exo_DB *exo)
/*
 * Extracts shell element stu from bulk element stu
 * based on Exo/Patran side ID (no iteration necessary).
 */
{
  int err = 0;
  int bulk_eptr = Proc_Connect_Ptr[bulk_elem];
  int shell_eptr = Proc_Connect_Ptr[shell_elem];
  int bulk_type = Elem_Type(exo, bulk_elem);
  int shell_type = Elem_Type(exo, shell_elem);

  /* 2D bulk element, 1D shell element */
  if (bulk_dim == 2) {
    if (bulk_type == LINEAR_TRI)
      GOMA_EH(GOMA_ERROR, "find_stu_on_shell not fixed for 2d bulk triangles");

    xi2[1] = xi2[2] = 0.0;
    switch (id_side) {
    case 1:
      if (Proc_Elem_Connect[bulk_eptr + 0] == Proc_Elem_Connect[shell_eptr])
        xi2[0] = s;
      else if (Proc_Elem_Connect[bulk_eptr + 1] == Proc_Elem_Connect[shell_eptr])
        xi2[0] = -s;
      else
        GOMA_EH(GOMA_ERROR, "Error in find_stu_on_shell");
      break;
    case 2:
      if (Proc_Elem_Connect[bulk_eptr + 1] == Proc_Elem_Connect[shell_eptr])
        xi2[0] = t;
      else if (Proc_Elem_Connect[bulk_eptr + 2] == Proc_Elem_Connect[shell_eptr])
        xi2[0] = -t;
      else
        GOMA_EH(GOMA_ERROR, "Error in find_stu_on_shell");
      break;
    case 3:
      if (Proc_Elem_Connect[bulk_eptr + 2] == Proc_Elem_Connect[shell_eptr])
        xi2[0] = -s;
      else if (Proc_Elem_Connect[bulk_eptr + 3] == Proc_Elem_Connect[shell_eptr])
        xi2[0] = s;
      else
        GOMA_EH(GOMA_ERROR, "Error in find_stu_on_shell");
      break;
    case 4:
      if (Proc_Elem_Connect[bulk_eptr + 3] == Proc_Elem_Connect[shell_eptr])
        xi2[0] = -t;
      else if (Proc_Elem_Connect[bulk_eptr + 0] == Proc_Elem_Connect[shell_eptr])
        xi2[0] = t;
      else
        GOMA_EH(GOMA_ERROR, "Error in find_stu_on_shell");
      break;
    default:
      err = -1;
      break;
    }
  }

  /* 3D bulk element, 2D shell element */
  else if (bulk_dim == 3 && shell_type == BILINEAR_SHELL) {
    /* GOMA_EH(GOMA_ERROR,"DRN believes that this is too simplified, more tests are needed as done
     * above for 2D"); */
    xi2[2] = 0.0;
    switch (id_side) {
    case 1:
      xi2[0] = s;
      xi2[1] = u;
      break;
    case 2:
      xi2[0] = t;
      xi2[1] = u;
      break;
    case 3:
      xi2[0] = -s;
      xi2[1] = u;
      break;
    case 4:
      xi2[0] = -t;
      xi2[1] = u;
      break;
    case 5:
      xi2[0] = s;
      xi2[1] = -t;
      break;
    case 6:
      xi2[0] = s;
      xi2[1] = t;
      break;
    default:
      err = -1;
      break;
    }
  } else if (bulk_dim == 3 && shell_type == BILINEAR_TRISHELL) {
    /* GOMA_EH(GOMA_ERROR,"DRN believes that this is too simplified, more tests are needed as done
     * above for 2D"); */
    xi2[2] = 0.0;
    switch (id_side) {
    case 1:
      xi2[0] = s;
      xi2[1] = u;
      break;
    case 2:
      xi2[0] = t;
      xi2[1] = u;
      break;
    case 3:
      xi2[0] = u;
      xi2[1] = t;
      break;
    case 4:
      xi2[0] = t;
      xi2[1] = s;
      break;
    default:
      err = -1;
      break;
    }
  }

  /* Only these two cases are handled here */
  else {
    err = -1;
  }

  return err;
} /* End of find_stu_on_shell() */

int find_stu_on_bulk(
    const int id_side, const int bulk_dim, const double s, const double t, double xi2[DIM])
/*
 * Extracts bulk element stu from shell element stu
 * based on Exo/Patran side ID (no iteration necessary).
 */
{
  int err = 0;

  /* 2D bulk element, 1D shell element */
  if (bulk_dim == 2) {
    xi2[2] = 0.0;
    switch (id_side) {
    case 1:
      xi2[0] = s;
      xi2[1] = -1.0;
      break;
    case 2:
      xi2[0] = 1.0;
      xi2[1] = s;
      break;
    case 3:
      xi2[0] = -s;
      xi2[1] = 1.0;
      break;
    case 4:
      xi2[0] = -1.0;
      xi2[1] = -s;
      break;
    default:
      err = -1;
      break;
    }
  }

  /* 3D bulk element, 2D shell element */
  else if (bulk_dim == 3) {
    switch (id_side) {
    case 1:
      xi2[0] = s;
      xi2[1] = 1.0;
      xi2[2] = t;
      break;
    case 2:
      xi2[0] = 1.0;
      xi2[1] = s;
      xi2[2] = t;
      break;
    case 3:
      xi2[0] = -s;
      xi2[1] = 1.0;
      xi2[2] = t;
      break;
    case 4:
      xi2[0] = -1.0;
      xi2[1] = -s;
      xi2[2] = t;
      break;
    case 5:
      xi2[0] = s;
      xi2[1] = -t;
      xi2[2] = -1.0;
      break;
    case 6:
      xi2[0] = s;
      xi2[1] = t;
      xi2[2] = 1.0;
      break;
    default:
      err = -1;
      break;
    }
  }

  /* Only these two cases are handled here */
  else {
    err = -1;
  }

  return err;
} /* End of find_stu_on_bulk() */
/****************************************************************************/
int shell_normal_div_s(dbl *div_s_nv, dbl d_div_s_nv_dnv[DIM][MDE], dbl d_div_s_nv_dmesh[DIM][MDE])

/************************************************************************
*
* shell_normal_div_s
*
*
*
*  Returns:
* ------------
*   div_s_nv,
*   dbl d_div_s_nv_dnv[DIM][MDE],
*   dbl d_div_s_nv_dmesh[DIM][MDE])

***********************************************************************/
{
  int eqn = R_SHELL_NORMAL1;
  int i, j, k, p, q, r, dofs, pdim, node, index;
  int b;
  double phi_j;
  // double sJdet;
  double grad_nv[DIM][DIM], d_grad_nv_dnv[DIM][DIM][DIM][MDE];
  double d_grad_nv_dmesh[DIM][DIM][DIM][MDE];
  double gradphi[MDE][DIM], sJ[DIM], sB[DIM], local_x[DIM][MDE];
  double dsB[DIM][DIM][MDE], d_gradphi_dmesh[DIM][DIM][DIM][MDE];
  double f;
  BASIS_FUNCTIONS_STRUCT *bfe = bf[eqn];

  if (!pd->e[pg->imtrx][eqn])
    GOMA_EH(GOMA_ERROR, "No normal equations on shell!");
  dofs = ei[pg->imtrx]->dof[eqn];
  pdim = pd->Num_Dim;

  /* Algorithm not tested for 3D yet! */
  if (pdim == 3)
    GOMA_EH(GOMA_ERROR, "Cannot handle shell elements in 3D yet!");

  memset(&(sJ[0]), 0, DIM * sizeof(double));
  memset(&(sB[0]), 0, DIM * sizeof(double));
  memset(&(dsB[0][0][0]), 0, DIM * DIM * MDE * sizeof(double));
  memset(&(gradphi[0][0]), 0, MDE * DIM * sizeof(double));
  memset(&(d_gradphi_dmesh[0][0][0][0]), 0, DIM * DIM * DIM * MDE * sizeof(double));
  memset(&(grad_nv[0][0]), 0, DIM * DIM * sizeof(double));
  memset(&(d_grad_nv_dnv[0][0][0][0]), 0, DIM * DIM * DIM * MDE * sizeof(double));
  memset(&(d_grad_nv_dmesh[0][0][0][0]), 0, DIM * DIM * DIM * MDE * sizeof(double));

  /*
   * Since beer_belly() and load_bf_grad() do not supply all components
   * of the physical coordinate derivatives of basis functions and
   * gradients for shell elements, calculate them locally for now.
   */
  // double sJsum = 0.0;
  for (j = 0; j < pdim; j++) {
    for (k = 0; k < dofs; k++) {
      node = ei[pg->imtrx]->dof_list[eqn][k];
      index = Proc_Elem_Connect[Proc_Connect_Ptr[ei[pg->imtrx]->ielem] + node];
      local_x[j][k] = Coor[j][index] + *esp->d[j][k];
      sJ[j] += local_x[j][k] * bf[eqn]->dphidxi[k][0];
    }
    // sJsum += sJ[j] * sJ[j];
  }
  //  sJdet = sqrt(sJsum);

  /* Inverse Jacobian */
  for (j = 0; j < pdim; j++) {
    // sB[j] = sJ[j] / sJdet;
    //  This is B[j][0] in beer_beely
    sB[j] = 1.0 / sJ[j];
  }

  /* Local gradphi */
  for (p = 0; p < pdim; p++) {
    for (j = 0; j < dofs; j++) {
      gradphi[j][p] = bf[eqn]->dphidxi[j][0] * sB[p] / fv->h[p];
    }
  }

  /* Inverse Jacobian mesh derivatives */
  for (i = 0; i < pdim; i++) {
    for (j = 0; j < pdim; j++) {
      r = ((i == j) ? 1 : -1);
      for (k = 0; k < dofs; k++) {
        // dsB[i][j][k] = r * sJ[1-i] * sJ[1-j] * bf[eqn]->dphidxi[k][0]
        //	/ (sJdet * sJdet * sJdet);
        //
        //  This is dB[j][0][j][k] from beer_belly
        dsB[i][j][k] = -bf[eqn]->dphidxi[k][0] / (sJ[j] * sJ[j]);
      }
    }
  }

  /* Local gradphi mesh derivatives */
  /*
  for (p = 0; p < pdim; p++)
    {
      for (i = 0; i < dofs; i++)
        {
          for (j = 0 ; j < pdim; j++)
            {
              for (k = 0; k < dofs; k++)
                {
                  d_gradphi_dmesh[i][p][j][k] = dsB[p][j][k] * bf[eqn]->dphidxi[i][0] / fv->h[p];
                }
            }
        }
    }
  */
  double g[DIM], g2[DIM];
  for (p = 0; p < VIM; p++) {
    g[p] = 1.0 / fv->h[p];
    g2[p] = g[p] * g[p];
  }

  for (p = 0; p < pdim; p++) {
    for (i = 0; i < dofs; i++) {
      for (b = 0; b < pdim; b++) {
        f = -g2[p] * (fv->hq[p][b]);
        phi_j = bf[eqn]->phi[j];
        for (j = 0; j < dofs; j++) {
          phi_j = bf[eqn]->phi[j];

          d_gradphi_dmesh[i][p][b][j] =
              (f * phi_j * sB[p] + g[p] * bfe->dB[p][0][b][j]) * bf[eqn]->dphidxi[i][0];
        }
      }
    }
  }

  /* Local normal vector gradient, mesh and self sensitivity */
  for (p = 0; p < VIM; p++) {
    for (q = 0; q < VIM; q++) {
      grad_nv[p][q] = 0.0;

      for (b = 0; b < pdim; b++) {
        for (j = 0; j < dofs; j++) {
          grad_nv[p][q] += *(esp->n[b][j]) * bf[eqn]->grad_phi_e[j][b][p][q];
          // grad_nv[p][q] += *(esp->n[q][i]) * gradphi[i][p];

          // d_grad_nv_dnv[p][q][r][i] += (double)delta(q,r) * gradphi[i][p];
          d_grad_nv_dnv[p][q][b][j] += bf[eqn]->grad_phi_e[j][b][p][q];
          for (r = 0; r < pdim; r++) {
            for (i = 0; i < dofs; i++) {
              //  d_grad_nv_dmesh[p][q][b][k] += *esp->n[b][i] * d_gradphi_dmesh[i][p][b][k];
              d_grad_nv_dmesh[p][q][b][j] +=
                  *esp->n[r][i] * bf[eqn]->d_grad_phi_e_dmesh[i][r][p][q][b][j];
            }
          }
        }
      }
    }
  }
#ifdef DEBUG_HKM
  /*
    for (p = 0; p < VIM; p++)
    {
    for (q = 0; q < VIM; q++)
    {
    printf("shell_normal_div_s: grad_nv[%d][%d] = %g\n", p, q, grad_nv[p][q]);
    }
    }
  */
#endif
  *div_s_nv = 0.0;
  memset(&(d_div_s_nv_dnv[0][0]), 0, DIM * MDE * sizeof(double));
  for (p = 0; p < VIM; p++) {
    for (q = 0; q < VIM; q++) {
      *div_s_nv += ((double)delta(p, q) - fv->n[p] * fv->n[q]) * grad_nv[p][q];

      for (j = 0; j < dofs; j++) {
        phi_j = bf[eqn]->phi[j];
        for (r = 0; r < pdim; r++) {
          d_div_s_nv_dmesh[r][j] +=
              ((double)delta(p, q) - fv->n[q] * fv->n[p]) * d_grad_nv_dmesh[q][p][r][j];
          d_div_s_nv_dnv[r][j] +=
              ((double)delta(p, q) - fv->n[q] * fv->n[p]) * d_grad_nv_dnv[q][p][r][j];
        }
      }
    }
  }

  for (p = 0; p < VIM; p++) {
    for (q = 0; q < VIM; q++) {
      for (j = 0; j < dofs; j++) {
        phi_j = bf[eqn]->phi[j];
        d_div_s_nv_dnv[p][j] -= fv->n[q] * phi_j * grad_nv[q][p];
      }
    }
  }

  for (p = 0; p < VIM; p++) {
    for (q = 0; q < VIM; q++) {
      for (j = 0; j < dofs; j++) {
        phi_j = bf[eqn]->phi[j];
        d_div_s_nv_dnv[q][j] -= fv->n[p] * phi_j * grad_nv[q][p];
      }
    }
  }

#ifdef DEBUG_HKM
  // printf("shell_normal_div_s: *div_s_nv = %g\n", *div_s_nv);
#endif
  return 0;
}

/* End of file mm_shell_util.c */

/****************************************************************************/

void shell_determinant_and_normal(
    const int ielem,          /* current element number               */
    const int iconnect_ptr,   /* Pointer to beginning of connectivity
                               * list for current element             */
    const int nodes_per_elem, /* number of nodes in the element       */
    const int ielem_surf_dim, /* physical dimension of the element
                               * surface (0, 1, 2)                    */
    const int id_side)        /* shell side (bottom or top) (exo/patran convention)  */

/************************************************************************
 *
 * shell_determinant_and_normal()
 *
 *      Function which calculates the shell surface determinant and shell surface
 *      normal at a local surface quadrature point. The function also
 *      calculates the sensitivities of those quantities to the
 *      mesh positions, if the mesh positions are part of the solution
 *      vector.
 *
 *      Adapted from surface_determinant_and_normal by PR Schunk (9/2/2009)
 *      Looks like id_side is always 1, so pretty much not used - RBS
 *
 *  Returns:
 * ------------
 *   fv->sdet = surface determinant at the quadrature point
 *   fv->snormal[] = surface normal at the quadrature point
 *   fv->dsurfdet_dx[][] = sensitivity of fv->sdet wrt mesh displacements.
 *   fv->dsnormal_dx[][] = sensitivity of fv->snormal[]
 *                         wrt mesh displacements
 ***********************************************************************/
{
  int i, inode, a, b, p, q;
  int ShapeVar, ldof;
  int DeformingMesh;
  double r_det, det_h01, r_det_h01, d_det_h01_x;
  double phi_i;
  int siz;
  double T[DIM - 1][DIM], t[DIM - 1][DIM]; /* t = J . T */
  double dt_x[DIM - 1][DIM][DIM][MDE];     /* d(t) / d(x_j) */
  struct Basis_Functions *map_bf;

  DeformingMesh = pd->e[pg->imtrx][R_MESH1];
  DeformingMesh = (upd->ep[pg->imtrx][R_MESH1] != -1);
  ShapeVar = pd->ShapeVar;

  siz = MAX_PDIM * MDE * sizeof(double);
  memset(fv->dsurfdet_dx, 0, siz);
  siz = MAX_PDIM * MDE * MAX_PDIM * sizeof(double);
  memset(fv->dsnormal_dx, 0, siz);

  map_bf = bf[ShapeVar];

  /* Here ielem_surf_dim is the shell dimension, which is either 1 or 2 */
  /* It is set as ei[pg->imtrx]->ielem_dim (not pd->Num_Dim which is usual 1 dimension more) */

  if (ielem_surf_dim == 0) /* get out quickly */
  {
    fv->sdet = 1.0;
    fv->snormal[0] = 1.0;
    return;
  }

  /* define space of surface */
  siz = (DIM - 1) * DIM * sizeof(double);
  memset(T, 0, siz);
  memset(t, 0, siz);
  /*  since T & t are zeroed, only need to set nonzero elements */
  switch (ielem_surf_dim) {
  case 1:
    switch (ei[pg->imtrx]->ielem_shape) {
    case LINE_SEGMENT:
      T[0][0] = 1.;
      break;
    default:
      GOMA_EH(GOMA_ERROR, "Invalid shape");
      break;
    }
    break;
  case 2:
    switch (ei[pg->imtrx]->ielem_shape) {
    case SHELL:
    case TRISHELL:
      /* if (id_side == 5 )
        {
          T[0][0] =  1.; T[0][1] =  0.; T[0][2] =  0.;
          T[1][0] =  0.; T[1][1] = -1.; T[1][2] =  0.;
        }
      else if ( id_side == 6 )
        {
          T[0][0] =  1.; T[0][1] =  0.; T[0][2] =  0.;
          T[1][0] =  0.; T[1][1] =  1.; T[1][2] =  0.;
        }
      else
        {
          GOMA_EH(GOMA_ERROR, "Incorrect side for HEXAHEDRAL");
          } */

      /*Not sure how we want to handle this, viz. compatible with the
        bulk element surf or just define the isoparametric integration.  For
        now we will hardwire to "patran id_side 6" in which Zeta is fixed */

      /* Here we want to load up a 3D Jacobian matrix which we can use the T matrix to
       * for gradient conversion.  beer_belly loads up 2D jacobian for the shell-only case, so
       * we cannot use that.
       */

      T[0][0] = 1.;
      T[1][1] = 1.;
      break;

    default:
      GOMA_EH(GOMA_ERROR, "Invalid shape. Lubrication capability requires shell elements");
      break;
    }
    break;
  }

  /* Use T matrix which 2x3 and J which is 3x3 to transform grad operator */
  /* This is a load_bf_derivs equivalent for shells */

  /* transform from element to physical coords */
  /* NOTE: This does not work for 2D bars because map_bf->J
     is hardwired in beerbelly to already compute this */

  for (p = 0; p < ielem_surf_dim; p++) {
    for (a = 0; a < pd->Num_Dim; a++) {
      t[p][a] = 0.;
      for (b = 0; b < pd->Num_Dim; b++) {
        t[p][a] += map_bf->J[b][a] * T[p][b] * fv->h[b];
      }
    }
  }

  if (af->Assemble_Jacobian && DeformingMesh) {
    for (p = 0; p < ielem_surf_dim; p++) {
      for (a = 0; a < pd->Num_Dim; a++) {
        for (i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
          inode = Proc_Elem_Connect[iconnect_ptr + i];
          ldof = ei[pg->imtrx]->ln_to_dof[ShapeVar][i];
          if (Dolphin[pg->imtrx][inode][MESH_DISPLACEMENT1] > 0) {
            phi_i = map_bf->phi[ldof];

            for (q = 0; q < pd->Num_Dim; q++) {
              dt_x[p][a][q][ldof] = 0.;
              for (b = 0; b < ei[pg->imtrx]->ielem_dim; b++) {
                dt_x[p][a][q][ldof] += T[p][b] * (map_bf->dJ[b][a][q][ldof] * fv->h[b] +
                                                  map_bf->J[b][a] * phi_i * fv->hq[b][q]);
              }
            }
          }
        }
      }
    }
  }

  if (ielem_surf_dim == 1) {

    /* N.B. NEXT Stage of Change: Correct Mapbf->J for 1D case in beer_belly */
    // GOMA_EH(GOMA_ERROR, "you should not be here until mapbf->J in beerbelly is corrected for 1D
    // case");
    /* calculate surface determinant using the coordinate scale factors
     * for orthogonal curvilinear coordinates */
    det_h01 = sqrt(t[0][0] * t[0][0] + t[0][1] * t[0][1]);
    r_det_h01 = 1. / det_h01;

    fv->sdet = fv->h[2] * det_h01;
    r_det = 1. / fv->sdet;

    /* calculate surface normal using the coordinate scale factors
     * for orthogonal curvilinear coordinates */
    fv->snormal[0] = t[0][1] * r_det_h01;
    fv->snormal[1] = -t[0][0] * r_det_h01;
    fv->snormal[2] = 0.;

    /* Calculate sensitivity w.r.t. mesh, if applicable */
    if (af->Assemble_Jacobian && DeformingMesh) {
      for (i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
        inode = Proc_Elem_Connect[iconnect_ptr + i];
        ldof = ei[pg->imtrx]->ln_to_dof[ShapeVar][i];
        if (Dolphin[pg->imtrx][inode][MESH_DISPLACEMENT1] > 0) {
          phi_i = map_bf->phi[ldof];

          for (a = 0; a < pd->Num_Dim; a++) {
            d_det_h01_x = 0.;
            for (b = 0; b < pd->Num_Dim; b++) {
              d_det_h01_x += r_det_h01 * t[0][b] * dt_x[0][b][a][ldof];
            }

            /* calculate sensitivity of surface determinant using the coordinate scale factors
             * for orthogonal curvilinear coordinates */
            fv->dsurfdet_dx[a][ldof] = fv->hq[2][a] * phi_i * det_h01 + fv->h[2] * d_det_h01_x;

            /* calculate sensitivity of surface normal using the coordinate scale factors
             * for orthogonal curvilinear coordinates */
            fv->dsnormal_dx[0][a][ldof] =
                r_det_h01 * (dt_x[0][1][a][ldof] - r_det_h01 * t[0][1] * d_det_h01_x);
            fv->dsnormal_dx[1][a][ldof] =
                r_det_h01 * (-dt_x[0][0][a][ldof] + r_det_h01 * t[0][0] * d_det_h01_x);
          }
        }
      }
    }

    if (mp->ehl_integration_kind == SIK_S) {
      double det_J;
      double d_det_J_dmeshkj[DIM][MDE];
      memset(d_det_J_dmeshkj, 0.0, sizeof(double) * DIM * MDE);
      detJ_2d_bar(&det_J, d_det_J_dmeshkj);
      fv->sdet = det_J;
      for (int k = 0; k < DIM; k++) {
        for (int j = 0; j < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; j++) {
          fv->dsurfdet_dx[k][j] = d_det_J_dmeshkj[k][j];
        }
      }
    }
  } else if (ielem_surf_dim == 2) {
    double nx = t[0][1] * t[1][2] - t[0][2] * t[1][1];
    double ny = t[0][2] * t[1][0] - t[0][0] * t[1][2];
    double nz = t[0][0] * t[1][1] - t[0][1] * t[1][0];
    double d_nx_x, d_ny_x, d_nz_x;

    /* calculate surface determinant using the coordinate scale factors
     * for orthogonal curvilinear coordinates */
    fv->sdet = sqrt(nx * nx + ny * ny + nz * nz);
    r_det = 1. / fv->sdet;

    /* calculate surface normal using the coordinate scale factors
     * for orthogonal curvilinear coordinates */
    fv->snormal[0] = r_det * nx;
    fv->snormal[1] = r_det * ny;
    fv->snormal[2] = r_det * nz;

    /* Calculate sensitivity w.r.t. mesh, if applicable */
    if (af->Assemble_Jacobian && DeformingMesh) {
      for (i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
        inode = Proc_Elem_Connect[iconnect_ptr + i];
        ldof = ei[pg->imtrx]->ln_to_dof[ShapeVar][i];
        if (Dolphin[pg->imtrx][inode][MESH_DISPLACEMENT1] > 0) {
          phi_i = map_bf->phi[ldof];

          for (a = 0; a < pd->Num_Dim; a++) {
            d_nx_x = t[0][1] * dt_x[1][2][a][ldof] + dt_x[0][1][a][ldof] * t[1][2] -
                     t[0][2] * dt_x[1][1][a][ldof] - dt_x[0][2][a][ldof] * t[1][1];
            d_ny_x = t[0][2] * dt_x[1][0][a][ldof] + dt_x[0][2][a][ldof] * t[1][0] -
                     t[0][0] * dt_x[1][2][a][ldof] - dt_x[0][0][a][ldof] * t[1][2];
            d_nz_x = t[0][0] * dt_x[1][1][a][ldof] + dt_x[0][0][a][ldof] * t[1][1] -
                     t[0][1] * dt_x[1][0][a][ldof] - dt_x[0][1][a][ldof] * t[1][0];

            /* calculate sensitivity of surface determinant using the coordinate scale factors
             * for orthogonal curvilinear coordinates */
            fv->dsurfdet_dx[a][ldof] = r_det * (nx * d_nx_x + ny * d_ny_x + nz * d_nz_x);

            /* calculate sensitivity of surface normal using the coordinate scale factors
             * for orthogonal curvilinear coordinates */
            fv->dsnormal_dx[0][a][ldof] = r_det * (d_nx_x - r_det * nx * fv->dsurfdet_dx[a][ldof]);
            fv->dsnormal_dx[1][a][ldof] = r_det * (d_ny_x - r_det * ny * fv->dsurfdet_dx[a][ldof]);
            fv->dsnormal_dx[2][a][ldof] = r_det * (d_nz_x - r_det * nz * fv->dsurfdet_dx[a][ldof]);
          }
        }
      }
    }
  }
}

/****************************************************************************/

void shell_tangents(dbl t0[DIM],
                    dbl t1[DIM],
                    dbl dt0_dx[DIM][DIM][MDE],
                    dbl dt1_dx[DIM][DIM][MDE],
                    dbl dt0_dnormal[DIM][DIM][MDE],
                    dbl dt1_dnormal[DIM][DIM][MDE])
/******************************************************************************
 *
 * shell_tangents()
 *
 *    Wrapper routine to select which function is used to calculate shell
 *    tangents and curvatures and their sensitivities.
 *
 *    t0 = dx/ds normalized
 *    t1 = dx/dt normalized
 *
 * Andrew Cochrane <acochrane@gmail.com>
 *
 ******************************************************************************/
{
  switch (mp->shell_tangent_model) {
  case ISOPARAMETRIC:
    shell_tangents_isoparametric(t0, t1, dt0_dx, dt1_dx);
    break;
  case SEEDED:
    shell_tangents_seeded(t0, t1, dt0_dnormal, dt1_dnormal);
    break;
  default:
    shell_tangents_isoparametric(t0, t1, dt0_dx, dt1_dx);
    break;
  }
  return;
}

void shell_tangents_isoparametric(dbl t0[DIM],
                                  dbl t1[DIM],
                                  dbl dt0_dx[DIM][DIM][MDE],
                                  dbl dt1_dx[DIM][DIM][MDE])
/******************************************************************************
 *
 * shell_tangents_isoparametric()
 *
 *    Routine to calculate shell tangents and curvatures with their sensitivities
 *    Here, we use isoparametric coordinates s and t to form tangent vectors
 *
 *    t0 = dx/ds normalized
 *    t1 = dx/dt normalized
 *
 * Kristianto Tjiptowidjojo tjiptowi@unm.edu
 *
 ******************************************************************************/
{
  /*
   * Integers and indices
   */

  int var, dim, a, b;
  int j;

  /*
   * Local quantities
   */

  dbl t0_norm, t1_norm;
  dbl d_t0_norm_dx[DIM][MDE];
  dbl d_t1_norm_dx[DIM][MDE];

  dim = pd->Num_Dim;

  /******* TANGENT **********/
  for (a = 0; a < dim; a++) {
    t0[a] = bf[MESH_DISPLACEMENT1]->J[0][a];
    t1[a] = bf[MESH_DISPLACEMENT1]->J[1][a];
  }

  t0_norm = sqrt(t0[0] * t0[0] + t0[1] * t0[1] + t0[2] * t0[2]);
  t1_norm = sqrt(t1[0] * t1[0] + t1[1] * t1[1] + t1[2] * t1[2]);

  for (a = 0; a < dim; a++) {
    t0[a] = t0[a] / t0_norm;
    t1[a] = t1[a] / t1_norm;
  }

  if ((dt0_dx != NULL) && (dt1_dx != NULL)) {
    var = MESH_DISPLACEMENT1;
    if (pd->v[pg->imtrx][var]) {
      memset(d_t0_norm_dx, 0.0, sizeof(double) * DIM * MDE);
      memset(d_t1_norm_dx, 0.0, sizeof(double) * DIM * MDE);

      for (b = 0; b < dim; b++) {
        var = MESH_DISPLACEMENT1 + b;

        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_t0_norm_dx[b][j] = (1.0 / t0_norm) * bf[var]->J[0][b] * bf[var]->dphidxi[j][0];
          d_t1_norm_dx[b][j] = (1.0 / t1_norm) * bf[var]->J[1][b] * bf[var]->dphidxi[j][1];
        }
      }

      memset(dt0_dx, 0.0, sizeof(double) * DIM * DIM * MDE);
      memset(dt1_dx, 0.0, sizeof(double) * DIM * DIM * MDE);

      for (b = 0; b < dim; b++) {
        var = MESH_DISPLACEMENT1 + b;

        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (a = 0; a < dim; a++) {
            dt0_dx[a][b][j] = (-1.0 / (t0_norm * t0_norm)) * d_t0_norm_dx[b][j] * bf[var]->J[0][a] +
                              (1.0 / t0_norm) * bf[var]->dphidxi[j][0] * delta(a, b);

            dt1_dx[a][b][j] = (-1.0 / (t1_norm * t1_norm)) * d_t1_norm_dx[b][j] * bf[var]->J[1][a] +
                              (1.0 / t1_norm) * bf[var]->dphidxi[j][1] * delta(a, b);
          }
        }
      }
    }
  }
} /* End of shell_tangents */

/****************************************************************************/

void shell_tangents_seeded(dbl t0[DIM],
                           dbl t1[DIM],
                           dbl dt0_dnormal[DIM][DIM][MDE],
                           dbl dt1_dnormal[DIM][DIM][MDE])
/******************************************************************************
 *
 * shell_tangents_seeded()
 *
 *    Routine to calculate shell tangents and their sensitivities
 *
 *    Here we use a seed vector S to generate first tangent vector: t0 = (I - nn) . s
 *
 *    Second tangent vector t1 = - N x t0
 *
 * Kristianto Tjiptowidjojo tjiptowi@unm.edu
 *
 ******************************************************************************/
{
  /*
   * Integers and indices
   */

  int var, dim, a, b;
  int j;

  /*
   * Local quantities
   */

  dbl seed[DIM];
  dbl n_dot_seed;
  dbl d_n_dot_seed_dnormal[DIM][MDE];

  dbl t0_unscaled[DIM];
  dbl d_t0_unscaled_dnormal[DIM][DIM][MDE];
  dbl t0_norm;
  dbl d_t0_norm_dnormal[DIM][MDE];

  dim = pd->Num_Dim;

  /********* SEED *************/

  seed[0] = mp->shell_tangent_seed_vec_const[0];
  seed[1] = mp->shell_tangent_seed_vec_const[1];
  seed[2] = mp->shell_tangent_seed_vec_const[2];

  /******** NORMAL . SEED ******/
  n_dot_seed = 0.0;
  for (a = 0; a < dim; a++) {
    n_dot_seed += fv->n[a] * seed[a];
  }

  memset(d_n_dot_seed_dnormal, 0.0, sizeof(double) * DIM * MDE);
  for (b = 0; b < dim; b++) {
    var = SHELL_NORMAL1 + b;
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      d_n_dot_seed_dnormal[b][j] = bf[var]->phi[j] * seed[b];
    }
  }

  /******* TANGENT 1 *********/
  memset(t0_unscaled, 0.0, sizeof(double) * DIM);
  for (a = 0; a < dim; a++) {
    t0_unscaled[a] = seed[a] - fv->n[a] * n_dot_seed;
  }

  memset(d_t0_unscaled_dnormal, 0.0, sizeof(double) * DIM * DIM * MDE);
  for (b = 0; b < dim; b++) {
    var = SHELL_NORMAL1 + b;
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      d_t0_unscaled_dnormal[0][b][j] =
          -delta(0, b) * bf[var]->phi[j] * n_dot_seed - fv->n[0] * d_n_dot_seed_dnormal[b][j];
      d_t0_unscaled_dnormal[1][b][j] =
          -delta(1, b) * bf[var]->phi[j] * n_dot_seed - fv->n[1] * d_n_dot_seed_dnormal[b][j];
      d_t0_unscaled_dnormal[2][b][j] =
          -delta(2, b) * bf[var]->phi[j] * n_dot_seed - fv->n[2] * d_n_dot_seed_dnormal[b][j];
    }
  }

  t0_norm = 0.0;
  for (a = 0; a < dim; a++) {
    t0_norm += t0_unscaled[a] * t0_unscaled[a];
  }
  t0_norm = sqrt(t0_norm);

  memset(d_t0_norm_dnormal, 0.0, sizeof(double) * DIM * MDE);
  for (b = 0; b < dim; b++) {
    var = SHELL_NORMAL1 + b;
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      d_t0_norm_dnormal[b][j] = (1.0 / t0_norm) * (t0_unscaled[0] * d_t0_unscaled_dnormal[0][b][j] +
                                                   t0_unscaled[1] * d_t0_unscaled_dnormal[1][b][j] +
                                                   t0_unscaled[2] * d_t0_unscaled_dnormal[2][b][j]);
    }
  }

  memset(t0, 0.0, sizeof(double) * DIM);
  for (a = 0; a < dim; a++) {
    t0[a] = t0_unscaled[a] / t0_norm;
  }

  if (dt0_dnormal != NULL) {
    memset(dt0_dnormal, 0.0, sizeof(double) * DIM * DIM * MDE);
    for (a = 0; a < dim; a++) {
      for (b = 0; b < dim; b++) {
        var = SHELL_NORMAL1 + b;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          dt0_dnormal[a][b][j] = (1.0 / t0_norm) * d_t0_unscaled_dnormal[a][b][j] -
                                 t0_unscaled[a] / (t0_norm * t0_norm) * d_t0_norm_dnormal[b][j];
        }
      }
    }
  }
  /******* TANGENT 2  = -NORMAL X TANGENT 1 *********/

  memset(t1, 0.0, sizeof(double) * DIM);
  t1[0] = -(fv->n[1] * t0[2] - fv->n[2] * t0[1]);
  t1[1] = -(fv->n[2] * t0[0] - fv->n[0] * t0[2]);
  t1[2] = -(fv->n[0] * t0[1] - fv->n[1] * t0[0]);

  if (dt1_dnormal != NULL) {
    memset(dt1_dnormal, 0.0, sizeof(double) * DIM * DIM * MDE);
    for (j = 0; j < ei[pg->imtrx]->dof[SHELL_NORMAL1]; j++) {
      dt1_dnormal[0][0][j] = -(fv->n[1] * dt0_dnormal[2][0][j] - fv->n[2] * dt0_dnormal[1][0][j]);
      dt1_dnormal[0][1][j] = -(bf[SHELL_NORMAL2]->phi[j] * t0[2] + fv->n[1] * dt0_dnormal[2][1][j] -
                               fv->n[2] * dt0_dnormal[1][1][j]);
      dt1_dnormal[0][2][j] = -(fv->n[1] * dt0_dnormal[2][2][j] - bf[SHELL_NORMAL3]->phi[j] * t0[1] -
                               fv->n[2] * dt0_dnormal[1][2][j]);

      dt1_dnormal[1][0][j] = -(fv->n[2] * dt0_dnormal[0][0][j] - bf[SHELL_NORMAL1]->phi[j] * t0[2] -
                               fv->n[0] * dt0_dnormal[2][0][j]);
      dt1_dnormal[1][1][j] = -(fv->n[2] * dt0_dnormal[0][1][j] - fv->n[0] * dt0_dnormal[2][1][j]);
      dt1_dnormal[1][2][j] = -(bf[SHELL_NORMAL3]->phi[j] * t0[0] + fv->n[2] * dt0_dnormal[0][2][j] -
                               fv->n[0] * dt0_dnormal[2][2][j]);

      dt1_dnormal[2][0][j] = -(bf[SHELL_NORMAL1]->phi[j] * t0[1] + fv->n[0] * dt0_dnormal[1][0][j] -
                               fv->n[1] * dt0_dnormal[0][0][j]);
      dt1_dnormal[2][1][j] = -(fv->n[0] * dt0_dnormal[1][1][j] - bf[SHELL_NORMAL2]->phi[j] * t0[0] -
                               fv->n[1] * dt0_dnormal[0][1][j]);
      dt1_dnormal[2][2][j] = -(fv->n[0] * dt0_dnormal[1][2][j] - fv->n[1] * dt0_dnormal[0][2][j]);
    }
  }
} /* End of shell_tangents_seeded */

/****************************************************************************/

void shell_stress_tensor(dbl TT[DIM][DIM],
                         dbl dTT_dx[DIM][DIM][DIM][MDE],
                         dbl dTT_dnormal[DIM][DIM][DIM][MDE])
/******************************************************************************
 *
 * shell_stress_tensor()
 *
 *    Routine to calculate shell stress tensor TT and its sensitivities
 *
 * Kristianto Tjiptowidjojo tjiptowi@unm.edu
 *
 ******************************************************************************/
{
  int var, dim, a, b, p;
  int j;

  dbl K, nu;

  dbl t0[DIM];
  dbl t1[DIM];
  dbl dt0_dx[DIM][DIM][MDE];
  dbl dt1_dx[DIM][DIM][MDE];
  dbl dt0_dnormal[DIM][DIM][MDE];
  dbl dt1_dnormal[DIM][DIM][MDE];

  dbl t0_dot_grad_d[DIM];
  dbl t1_dot_grad_d[DIM];
  dbl d_t0_dot_grad_d_dx[DIM][DIM][MDE];
  dbl d_t1_dot_grad_d_dx[DIM][DIM][MDE];
  dbl d_t0_dot_grad_d_dnormal[DIM][DIM][MDE];
  dbl d_t1_dot_grad_d_dnormal[DIM][DIM][MDE];

  dbl du_ds;
  dbl d_du_ds_dnormal[DIM][MDE];
  dbl d_du_ds_dx[DIM][MDE];

  dbl du_dt;
  dbl d_du_dt_dx[DIM][MDE];
  dbl d_du_dt_dnormal[DIM][MDE];

  dbl dv_ds;
  dbl d_dv_ds_dx[DIM][MDE];
  dbl d_dv_ds_dnormal[DIM][MDE];

  dbl dv_dt;
  dbl d_dv_dt_dnormal[DIM][MDE];
  dbl d_dv_dt_dx[DIM][MDE];

  /* Unpack variables from structures for local convenience... */
  dim = pd->Num_Dim;

  K = elc->exten_stiffness;
  nu = elc->poisson;

  /******* TANGENTS AND CURVATURES **********/

  memset(dt0_dx, 0.0, sizeof(double) * DIM * DIM * MDE);
  memset(dt1_dx, 0.0, sizeof(double) * DIM * DIM * MDE);
  memset(dt0_dnormal, 0.0, sizeof(double) * DIM * DIM * MDE);
  memset(dt1_dnormal, 0.0, sizeof(double) * DIM * DIM * MDE);

  shell_tangents(t0, t1, dt0_dx, dt1_dx, dt0_dnormal, dt1_dnormal);

  /******* SHELL DISPLACEMENTS AND THEIR SENSITIVITIES **********/

  memset(t0_dot_grad_d, 0.0, sizeof(double) * DIM);
  memset(t1_dot_grad_d, 0.0, sizeof(double) * DIM);
  for (a = 0; a < dim; a++) {
    for (b = 0; b < dim; b++) {
      t0_dot_grad_d[a] += t0[b] * fv->grad_d[b][a];
      t1_dot_grad_d[a] += t1[b] * fv->grad_d[b][a];
    }
  }

  du_ds = 0.0;
  du_dt = 0.0;
  dv_ds = 0.0;
  dv_dt = 0.0;
  for (a = 0; a < dim; a++) {
    du_ds += t0[a] * t0_dot_grad_d[a];

    du_dt += t0[a] * t1_dot_grad_d[a];

    dv_ds += t1[a] * t0_dot_grad_d[a];

    dv_dt += t1[a] * t1_dot_grad_d[a];
  }

  memset(d_t0_dot_grad_d_dx, 0.0, sizeof(double) * DIM * DIM * MDE);
  memset(d_t1_dot_grad_d_dx, 0.0, sizeof(double) * DIM * DIM * MDE);
  for (a = 0; a < dim; a++) {
    for (b = 0; b < dim; b++) {
      for (p = 0; p < dim; p++) {
        var = MESH_DISPLACEMENT1 + p;

        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_t0_dot_grad_d_dx[a][p][j] +=
              t0[b] * fv->d_grad_d_dmesh[b][a][p][j] + dt0_dx[b][p][j] * fv->grad_d[b][a];

          d_t1_dot_grad_d_dx[a][p][j] +=
              t1[b] * fv->d_grad_d_dmesh[b][a][p][j] + dt1_dx[b][p][j] * fv->grad_d[b][a];
        }
      }
    }
  }

  memset(d_t0_dot_grad_d_dnormal, 0.0, sizeof(double) * DIM * DIM * MDE);
  memset(d_t1_dot_grad_d_dnormal, 0.0, sizeof(double) * DIM * DIM * MDE);
  for (a = 0; a < dim; a++) {
    for (b = 0; b < dim; b++) {
      for (p = 0; p < dim; p++) {
        var = SHELL_NORMAL1 + p;

        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_t0_dot_grad_d_dnormal[a][p][j] += dt0_dnormal[b][p][j] * fv->grad_d[b][a];

          d_t1_dot_grad_d_dnormal[a][p][j] += dt1_dnormal[b][p][j] * fv->grad_d[b][a];
        }
      }
    }
  }

  memset(d_du_ds_dnormal, 0.0, sizeof(double) * DIM * MDE);
  memset(d_dv_dt_dnormal, 0.0, sizeof(double) * DIM * MDE);
  memset(d_du_dt_dnormal, 0.0, sizeof(double) * DIM * MDE);
  memset(d_dv_ds_dnormal, 0.0, sizeof(double) * DIM * MDE);
  for (b = 0; b < dim; b++) {
    var = SHELL_NORMAL1 + b;

    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

      for (a = 0; a < dim; a++) {
        d_du_ds_dnormal[b][j] +=
            dt0_dnormal[a][b][j] * t0_dot_grad_d[a] + t0[a] * d_t0_dot_grad_d_dnormal[a][b][j];

        d_du_dt_dnormal[b][j] +=
            dt0_dnormal[a][b][j] * t1_dot_grad_d[a] + t0[a] * d_t1_dot_grad_d_dnormal[a][b][j];

        d_dv_ds_dnormal[b][j] +=
            dt1_dnormal[a][b][j] * t0_dot_grad_d[a] + t1[a] * d_t0_dot_grad_d_dnormal[a][b][j];

        d_dv_dt_dnormal[b][j] +=
            dt1_dnormal[a][b][j] * t1_dot_grad_d[a] + t1[a] * d_t1_dot_grad_d_dnormal[a][b][j];
      }
    }
  }

  memset(d_du_ds_dx, 0.0, sizeof(double) * DIM * MDE);
  memset(d_du_dt_dx, 0.0, sizeof(double) * DIM * MDE);
  memset(d_dv_ds_dx, 0.0, sizeof(double) * DIM * MDE);
  memset(d_dv_dt_dx, 0.0, sizeof(double) * DIM * MDE);
  for (b = 0; b < dim; b++) {
    var = MESH_DISPLACEMENT1 + b;

    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      for (a = 0; a < dim; a++) {
        d_du_ds_dx[b][j] +=
            dt0_dx[a][b][j] * t0_dot_grad_d[a] + t0[a] * d_t0_dot_grad_d_dx[a][b][j];

        d_du_dt_dx[b][j] +=
            dt0_dx[a][b][j] * t1_dot_grad_d[a] + t0[a] * d_t1_dot_grad_d_dx[a][b][j];

        d_dv_ds_dx[b][j] +=
            dt1_dx[a][b][j] * t0_dot_grad_d[a] + t1[a] * d_t0_dot_grad_d_dx[a][b][j];

        d_dv_dt_dx[b][j] +=
            dt1_dx[a][b][j] * t1_dot_grad_d[a] + t1[a] * d_t1_dot_grad_d_dx[a][b][j];
      }
    }
  }

  /******* SHELL STRESS TENSOR AND THEIR SENSITIVITIES **********/

  TT[0][0] = K * (du_ds + nu * dv_dt);
  TT[0][1] = 0.5 * K * (1.0 - nu) * (du_dt + dv_ds);
  TT[1][0] = 0.5 * K * (1.0 - nu) * (du_dt + dv_ds);
  TT[1][1] = K * (nu * du_ds + dv_dt);

  for (b = 0; b < dim; b++) {
    var = MESH_DISPLACEMENT1 + b;

    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      dTT_dx[0][0][b][j] = K * (d_du_ds_dx[b][j] + nu * d_dv_dt_dx[b][j]);

      dTT_dx[0][1][b][j] = 0.5 * K * (1.0 - nu) * (d_du_dt_dx[b][j] + d_dv_ds_dx[b][j]);

      dTT_dx[1][0][b][j] = 0.5 * K * (1.0 - nu) * (d_du_dt_dx[b][j] + d_dv_ds_dx[b][j]);

      dTT_dx[1][1][b][j] = K * (nu * d_du_ds_dx[b][j] + d_dv_dt_dx[b][j]);
    }
  }

  for (b = 0; b < dim; b++) {
    var = SHELL_NORMAL1 + b;

    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      dTT_dnormal[0][0][b][j] = K * (d_du_ds_dnormal[b][j] + nu * d_dv_dt_dnormal[b][j]);

      dTT_dnormal[0][1][b][j] =
          0.5 * K * (1.0 - nu) * (d_du_dt_dnormal[b][j] + d_dv_ds_dnormal[b][j]);
      dTT_dnormal[1][0][b][j] =
          0.5 * K * (1.0 - nu) * (d_du_dt_dnormal[b][j] + d_dv_ds_dnormal[b][j]);

      dTT_dnormal[1][1][b][j] = K * (nu * d_du_ds_dnormal[b][j] + d_dv_dt_dnormal[b][j]);
    }
  }
} /* End of shell_stress_tensor */

/****************************************************************************/

void shell_moment_tensor(dbl M[DIM][DIM],
                         dbl dM_dx[DIM][DIM][DIM][MDE],
                         dbl dM_dnormal[DIM][DIM][DIM][MDE],
                         dbl dM_dcurv0[DIM][DIM][MDE],
                         dbl dM_dcurv1[DIM][DIM][MDE])
/******************************************************************************
 *
 * shell_moment_tensor()
 *
 *    Routine to calculate shell moment tensor M and its sensitivities
 *
 * Kristianto Tjiptowidjojo tjiptowi@unm.edu
 *
 ******************************************************************************/
{
  int var, dim, a, b, p;
  int j;

  dbl D, nu;
  dbl K0, K1;

  dbl t0[DIM];
  dbl t1[DIM];
  dbl dt0_dx[DIM][DIM][MDE];
  dbl dt1_dx[DIM][DIM][MDE];
  dbl dt0_dnormal[DIM][DIM][MDE];
  dbl dt1_dnormal[DIM][DIM][MDE];

  dbl t0_dot_grad_d[DIM];
  dbl t1_dot_grad_d[DIM];
  dbl d_t0_dot_grad_d_dx[DIM][DIM][MDE];
  dbl d_t1_dot_grad_d_dx[DIM][DIM][MDE];
  dbl d_t0_dot_grad_d_dnormal[DIM][DIM][MDE];
  dbl d_t1_dot_grad_d_dnormal[DIM][DIM][MDE];

  dbl du_ds;
  dbl d_du_ds_dnormal[DIM][MDE];
  dbl d_du_ds_dx[DIM][MDE];

  dbl du_dt;
  dbl d_du_dt_dx[DIM][MDE];
  dbl d_du_dt_dnormal[DIM][MDE];

  dbl dv_ds;
  dbl d_dv_ds_dx[DIM][MDE];
  dbl d_dv_ds_dnormal[DIM][MDE];

  dbl dv_dt;
  dbl d_dv_dt_dnormal[DIM][MDE];
  dbl d_dv_dt_dx[DIM][MDE];

  dbl u, v;

  dbl d_u_dx[DIM][MDE];
  dbl d_v_dx[DIM][MDE];
  dbl d_u_dnormal[DIM][MDE];
  dbl d_v_dnormal[DIM][MDE];

  dbl dK0_ds;
  dbl dK1_ds;
  dbl dK0_dt;
  dbl dK1_dt;

  dbl d_dK0_ds_dx[DIM][MDE];
  dbl d_dK1_ds_dx[DIM][MDE];
  dbl d_dK0_dt_dx[DIM][MDE];
  dbl d_dK1_dt_dx[DIM][MDE];

  dbl d_dK0_ds_dnormal[DIM][MDE];
  dbl d_dK1_ds_dnormal[DIM][MDE];
  dbl d_dK0_dt_dnormal[DIM][MDE];
  dbl d_dK1_dt_dnormal[DIM][MDE];

  dbl d_dK0_ds_dcurv0[MDE];
  dbl d_dK1_ds_dcurv1[MDE];
  dbl d_dK0_dt_dcurv0[MDE];
  dbl d_dK1_dt_dcurv1[MDE];

  dbl phi_j;

  /* Unpack variables from structures for local convenience... */
  dim = pd->Num_Dim;

  D = elc->bend_stiffness;
  nu = elc->poisson;

  /******* TANGENTS **********/

  memset(dt0_dx, 0.0, sizeof(double) * DIM * DIM * MDE);
  memset(dt1_dx, 0.0, sizeof(double) * DIM * DIM * MDE);
  memset(dt0_dnormal, 0.0, sizeof(double) * DIM * DIM * MDE);
  memset(dt1_dnormal, 0.0, sizeof(double) * DIM * DIM * MDE);

  shell_tangents(t0, t1, dt0_dx, dt1_dx, dt0_dnormal, dt1_dnormal);

  /******* SHELL DISPLACEMENTS AND THEIR SENSITIVITIES **********/

  memset(t0_dot_grad_d, 0.0, sizeof(double) * DIM);
  memset(t1_dot_grad_d, 0.0, sizeof(double) * DIM);
  for (a = 0; a < dim; a++) {
    for (b = 0; b < dim; b++) {
      t0_dot_grad_d[a] += t0[b] * fv->grad_d[b][a];
      t1_dot_grad_d[a] += t1[b] * fv->grad_d[b][a];
    }
  }

  du_ds = 0.0;
  du_dt = 0.0;
  dv_ds = 0.0;
  dv_dt = 0.0;
  for (a = 0; a < dim; a++) {
    du_ds += t0[a] * t0_dot_grad_d[a];

    du_dt += t0[a] * t1_dot_grad_d[a];

    dv_ds += t1[a] * t0_dot_grad_d[a];

    dv_dt += t1[a] * t1_dot_grad_d[a];
  }

  memset(d_t0_dot_grad_d_dx, 0.0, sizeof(double) * DIM * DIM * MDE);
  memset(d_t1_dot_grad_d_dx, 0.0, sizeof(double) * DIM * DIM * MDE);
  for (a = 0; a < dim; a++) {
    for (b = 0; b < dim; b++) {
      for (p = 0; p < dim; p++) {
        var = MESH_DISPLACEMENT1 + p;

        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_t0_dot_grad_d_dx[a][p][j] +=
              t0[b] * fv->d_grad_d_dmesh[b][a][p][j] + dt0_dx[b][p][j] * fv->grad_d[b][a];

          d_t1_dot_grad_d_dx[a][p][j] +=
              t1[b] * fv->d_grad_d_dmesh[b][a][p][j] + dt1_dx[b][p][j] * fv->grad_d[b][a];
        }
      }
    }
  }

  memset(d_t0_dot_grad_d_dnormal, 0.0, sizeof(double) * DIM * DIM * MDE);
  memset(d_t1_dot_grad_d_dnormal, 0.0, sizeof(double) * DIM * DIM * MDE);
  for (a = 0; a < dim; a++) {
    for (b = 0; b < dim; b++) {
      for (p = 0; p < dim; p++) {
        var = SHELL_NORMAL1 + p;

        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_t0_dot_grad_d_dnormal[a][p][j] += dt0_dnormal[b][p][j] * fv->grad_d[b][a];

          d_t1_dot_grad_d_dnormal[a][p][j] += dt1_dnormal[b][p][j] * fv->grad_d[b][a];
        }
      }
    }
  }

  memset(d_du_ds_dnormal, 0.0, sizeof(double) * DIM * MDE);
  memset(d_dv_dt_dnormal, 0.0, sizeof(double) * DIM * MDE);
  memset(d_du_dt_dnormal, 0.0, sizeof(double) * DIM * MDE);
  memset(d_dv_ds_dnormal, 0.0, sizeof(double) * DIM * MDE);
  for (b = 0; b < dim; b++) {
    var = SHELL_NORMAL1 + b;

    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      phi_j = bf[var]->phi[j];

      for (a = 0; a < dim; a++) {
        d_du_ds_dnormal[b][j] +=
            dt0_dnormal[a][b][j] * t0_dot_grad_d[a] + t0[a] * d_t0_dot_grad_d_dnormal[a][b][j];

        d_du_dt_dnormal[b][j] +=
            dt0_dnormal[a][b][j] * t1_dot_grad_d[a] + t0[a] * d_t1_dot_grad_d_dnormal[a][b][j];

        d_dv_ds_dnormal[b][j] +=
            dt1_dnormal[a][b][j] * t0_dot_grad_d[a] + t1[a] * d_t0_dot_grad_d_dnormal[a][b][j];

        d_dv_dt_dnormal[b][j] +=
            dt1_dnormal[a][b][j] * t1_dot_grad_d[a] + t1[a] * d_t1_dot_grad_d_dnormal[a][b][j];
      }
    }
  }

  memset(d_du_ds_dx, 0.0, sizeof(double) * DIM * MDE);
  memset(d_du_dt_dx, 0.0, sizeof(double) * DIM * MDE);
  memset(d_dv_ds_dx, 0.0, sizeof(double) * DIM * MDE);
  memset(d_dv_dt_dx, 0.0, sizeof(double) * DIM * MDE);
  for (b = 0; b < dim; b++) {
    var = MESH_DISPLACEMENT1 + b;

    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      phi_j = bf[var]->phi[j];

      for (a = 0; a < dim; a++) {
        d_du_ds_dx[b][j] +=
            dt0_dx[a][b][j] * t0_dot_grad_d[a] + t0[a] * d_t0_dot_grad_d_dx[a][b][j];

        d_du_dt_dx[b][j] +=
            dt0_dx[a][b][j] * t1_dot_grad_d[a] + t0[a] * d_t1_dot_grad_d_dx[a][b][j];

        d_dv_ds_dx[b][j] +=
            dt1_dx[a][b][j] * t0_dot_grad_d[a] + t1[a] * d_t0_dot_grad_d_dx[a][b][j];

        d_dv_dt_dx[b][j] +=
            dt1_dx[a][b][j] * t1_dot_grad_d[a] + t1[a] * d_t1_dot_grad_d_dx[a][b][j];
      }
    }
  }

  u = 0;
  v = 0;

  for (a = 0; a < dim; a++) {
    u += t0[a] * fv->d[a];
    v += t1[a] * fv->d[a];
  }

  memset(d_u_dx, 0.0, sizeof(double) * DIM * MDE);
  memset(d_v_dx, 0.0, sizeof(double) * DIM * MDE);

  for (b = 0; b < dim; b++) {
    var = MESH_DISPLACEMENT1 + b;
    for (a = 0; a < dim; a++) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        phi_j = bf[var]->phi[j];

        d_u_dx[b][j] += dt0_dx[a][b][j] * fv->d[a] + t0[a] * delta(a, b) * phi_j;
        d_v_dx[b][j] += dt1_dx[a][b][j] * fv->d[a] + t1[a] * delta(a, b) * phi_j;
      }
    }
  }

  memset(d_u_dnormal, 0.0, sizeof(double) * DIM * MDE);
  memset(d_v_dnormal, 0.0, sizeof(double) * DIM * MDE);

  for (b = 0; b < dim; b++) {
    var = SHELL_NORMAL1 + b;
    for (a = 0; a < dim; a++) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_u_dnormal[b][j] += dt0_dnormal[a][b][j] * fv->d[a];
        d_v_dnormal[b][j] += dt1_dnormal[a][b][j] * fv->d[a];
      }
    }
  }

  /******* SHELL CURVATURES AND THEIR SENSITIVITIES **********/

  K0 = fv->sh_K;
  K1 = fv->sh_K2;

  dK0_ds = 0.0;
  dK1_ds = 0.0;
  dK0_dt = 0.0;
  dK1_dt = 0.0;

  for (a = 0; a < dim; a++) {
    dK0_ds += t0[a] * fv->grad_sh_K[a];
    dK1_ds += t0[a] * fv->grad_sh_K2[a];
    dK0_dt += t1[a] * fv->grad_sh_K[a];
    dK1_dt += t1[a] * fv->grad_sh_K2[a];
  }

  memset(d_dK0_ds_dx, 0.0, sizeof(double) * DIM * MDE);
  memset(d_dK1_ds_dx, 0.0, sizeof(double) * DIM * MDE);
  memset(d_dK0_dt_dx, 0.0, sizeof(double) * DIM * MDE);
  memset(d_dK1_dt_dx, 0.0, sizeof(double) * DIM * MDE);

  for (b = 0; b < dim; b++) {
    var = MESH_DISPLACEMENT1 + b;
    for (a = 0; a < dim; a++) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_dK0_ds_dx[b][j] +=
            dt0_dx[a][b][j] * fv->grad_sh_K[a] + t0[a] * fv->d_grad_sh_K_dmesh[a][b][j];
        d_dK1_ds_dx[b][j] +=
            dt0_dx[a][b][j] * fv->grad_sh_K2[a] + t0[a] * fv->d_grad_sh_K2_dmesh[a][b][j];
        d_dK0_dt_dx[b][j] +=
            dt1_dx[a][b][j] * fv->grad_sh_K[a] + t1[a] * fv->d_grad_sh_K_dmesh[a][b][j];
        d_dK1_dt_dx[b][j] +=
            dt1_dx[a][b][j] * fv->grad_sh_K2[a] + t1[a] * fv->d_grad_sh_K2_dmesh[a][b][j];
      }
    }
  }

  memset(d_dK0_ds_dnormal, 0.0, sizeof(double) * DIM * MDE);
  memset(d_dK1_ds_dnormal, 0.0, sizeof(double) * DIM * MDE);
  memset(d_dK0_dt_dnormal, 0.0, sizeof(double) * DIM * MDE);
  memset(d_dK1_dt_dnormal, 0.0, sizeof(double) * DIM * MDE);

  for (b = 0; b < dim; b++) {
    var = SHELL_NORMAL1 + b;
    for (a = 0; a < dim; a++) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_dK0_ds_dnormal[b][j] += dt0_dnormal[a][b][j] * fv->grad_sh_K[a];
        d_dK1_ds_dnormal[b][j] += dt0_dnormal[a][b][j] * fv->grad_sh_K2[a];
        d_dK0_dt_dnormal[b][j] += dt1_dnormal[a][b][j] * fv->grad_sh_K[a];
        d_dK1_dt_dnormal[b][j] += dt1_dnormal[a][b][j] * fv->grad_sh_K2[a];
      }
    }
  }

  memset(d_dK0_ds_dcurv0, 0.0, sizeof(double) * MDE);
  memset(d_dK1_ds_dcurv1, 0.0, sizeof(double) * MDE);
  memset(d_dK0_dt_dcurv0, 0.0, sizeof(double) * MDE);
  memset(d_dK1_dt_dcurv1, 0.0, sizeof(double) * MDE);

  var = SHELL_CURVATURE;
  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
    for (a = 0; a < dim; a++) {
      d_dK0_ds_dcurv0[j] += t0[a] * bf[var]->grad_phi[j][a];
      d_dK0_dt_dcurv0[j] += t1[a] * bf[var]->grad_phi[j][a];
    }
  }

  var = SHELL_CURVATURE2;
  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
    for (a = 0; a < dim; a++) {
      d_dK1_ds_dcurv1[j] += t0[a] * bf[var]->grad_phi[j][a];
      d_dK1_dt_dcurv1[j] += t1[a] * bf[var]->grad_phi[j][a];
    }
  }

  /******* SHELL MOMENT TENSOR AND THEIR SENSITIVITIES **********/

  int my_smt_kind = mp->shell_moment_tensor_model;

  switch (my_smt_kind) {
  case SMT_SIMPLE:
    M[0][0] = -D * (K0 + nu * K1);
    M[0][1] = -0.5 * D * (1.0 - nu) * (K0 + K1);
    M[1][0] = -0.5 * D * (1.0 - nu) * (K0 + K1);
    M[1][1] = -D * (nu * K0 + K1);
    break;
  case SMT_EXPANDED:
    M[0][0] = D * (K0 * du_ds + dK0_ds * u + nu * (K1 * dv_dt + dK1_dt * v));
    M[0][1] = 0.5 * D * (1.0 - nu) * (K0 * du_dt + dK0_dt * u + K1 * dv_ds + dK1_ds * v);
    M[1][0] = 0.5 * D * (1.0 - nu) * (K0 * du_dt + dK0_dt * u + K1 * dv_ds + dK1_ds * v);
    M[1][1] = D * (nu * (K0 * du_ds + dK0_ds * u) + K1 * dv_dt + dK1_dt * v);
    break;
  }

  for (b = 0; b < dim; b++) {
    var = MESH_DISPLACEMENT1 + b;

    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      switch (my_smt_kind) {
      case SMT_SIMPLE:
        dM_dx[0][0][b][j] = 0.0;
        dM_dx[0][1][b][j] = 0.0;
        dM_dx[1][0][b][j] = 0.0;
        dM_dx[1][1][b][j] = 0.0;
        break;
      case SMT_EXPANDED:
        dM_dx[0][0][b][j] =
            D * (K0 * d_du_ds_dx[b][j] + d_dK0_ds_dx[b][j] * u + dK0_ds * d_u_dx[b][j] +
                 (K1 * d_dv_dt_dx[b][j] + d_dK1_dt_dx[b][j] * v + dK1_dt * d_v_dx[b][j]) * nu);

        dM_dx[0][1][b][j] = 0.5 * D * (1.0 - nu) *
                            (K0 * d_du_dt_dx[b][j] + d_dK0_dt_dx[b][j] * u + dK0_dt * d_u_dx[b][j] +
                             K1 * d_dv_ds_dx[b][j] + d_dK1_ds_dx[b][j] * v + dK1_ds * d_v_dx[b][j]);

        dM_dx[1][0][b][j] = 0.5 * D * (1.0 - nu) *
                            (K0 * d_du_dt_dx[b][j] + d_dK0_dt_dx[b][j] * u + dK0_dt * d_u_dx[b][j] +
                             K1 * d_dv_ds_dx[b][j] + d_dK1_ds_dx[b][j] * v + dK1_ds * d_v_dx[b][j]);

        dM_dx[1][1][b][j] =
            D * ((K0 * d_du_ds_dx[b][j] + d_dK0_ds_dx[b][j] * u + dK0_ds * d_u_dx[b][j]) * nu +
                 K1 * d_dv_dt_dx[b][j] + d_dK1_dt_dx[b][j] * v + dK1_dt * d_v_dx[b][j]);

        break;
      }
    }
  }

  for (b = 0; b < dim; b++) {
    var = SHELL_NORMAL1 + b;
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      switch (my_smt_kind) {
      case SMT_SIMPLE:
        dM_dnormal[0][0][b][j] = 0.0;

        dM_dnormal[0][1][b][j] = 0.0;

        dM_dnormal[1][0][b][j] = 0.0;

        dM_dnormal[1][1][b][j] = 0.0;

        break;
      case SMT_EXPANDED:
        dM_dnormal[0][0][b][j] = D * (K0 * d_du_ds_dnormal[b][j] + d_dK0_ds_dnormal[b][j] * u +
                                      dK0_ds * d_u_dnormal[b][j] +
                                      (K1 * d_dv_dt_dnormal[b][j] + d_dK1_dt_dnormal[b][j] * v +
                                       dK1_dt * d_v_dnormal[b][j]) *
                                          nu);

        dM_dnormal[0][1][b][j] =
            0.5 * D * (1.0 - nu) *
            (K0 * d_du_dt_dnormal[b][j] + d_dK0_dt_dnormal[b][j] * u + dK0_dt * d_u_dnormal[b][j] +
             K1 * d_dv_ds_dnormal[b][j] + d_dK1_ds_dnormal[b][j] * v + dK1_ds * d_v_dnormal[b][j]);

        dM_dnormal[1][0][b][j] =
            0.5 * D * (1.0 - nu) *
            (K0 * d_du_dt_dnormal[b][j] + d_dK0_dt_dnormal[b][j] * u + dK0_dt * d_u_dnormal[b][j] +
             K1 * d_dv_ds_dnormal[b][j] + d_dK1_ds_dnormal[b][j] * v + dK1_ds * d_v_dnormal[b][j]);

        dM_dnormal[1][1][b][j] = D * ((K0 * d_du_ds_dnormal[b][j] + d_dK0_ds_dnormal[b][j] * u +
                                       dK0_ds * d_u_dnormal[b][j]) *
                                          nu +
                                      K1 * d_dv_dt_dnormal[b][j] + d_dK1_dt_dnormal[b][j] * v +
                                      dK1_dt * d_v_dnormal[b][j]);
        break;
      }
    }
  }

  var = SHELL_CURVATURE;
  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
    phi_j = bf[var]->phi[j];
    switch (my_smt_kind) {
    case SMT_SIMPLE:
      dM_dcurv0[0][0][j] = -D * (phi_j);
      dM_dcurv0[0][1][j] = -0.5 * D * (1.0 - nu) * (phi_j);
      dM_dcurv0[1][0][j] = -0.5 * D * (1.0 - nu) * (phi_j);
      dM_dcurv0[1][1][j] = -D * (phi_j)*nu;
      break;
    case SMT_EXPANDED:
      dM_dcurv0[0][0][j] = D * (phi_j * du_ds + d_dK0_ds_dcurv0[j] * u);
      dM_dcurv0[0][1][j] = 0.5 * D * (1.0 - nu) * (phi_j * du_dt + d_dK0_dt_dcurv0[j] * u);
      dM_dcurv0[1][0][j] = 0.5 * D * (1.0 - nu) * (phi_j * du_dt + d_dK0_dt_dcurv0[j] * u);
      dM_dcurv0[1][1][j] = D * (phi_j * du_ds + d_dK0_ds_dcurv0[j] * u) * nu;
      break;
    }
  }

  var = SHELL_CURVATURE2;
  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
    phi_j = bf[var]->phi[j];
    switch (my_smt_kind) {
    case SMT_SIMPLE:
      dM_dcurv1[0][0][j] = -D * (phi_j)*nu;
      dM_dcurv1[0][1][j] = -0.5 * D * (1.0 - nu) * (phi_j);
      dM_dcurv1[1][0][j] = -0.5 * D * (1.0 - nu) * (phi_j);
      dM_dcurv1[1][1][j] = -D * (phi_j);
      break;
    case SMT_EXPANDED:
      dM_dcurv1[0][0][j] = D * (phi_j * dv_dt + d_dK1_dt_dcurv1[j] * v) * nu;
      dM_dcurv1[0][1][j] = 0.5 * D * (1.0 - nu) * (phi_j * dv_ds + d_dK1_ds_dcurv1[j] * v);
      dM_dcurv1[1][0][j] = 0.5 * D * (1.0 - nu) * (phi_j * dv_ds + d_dK1_ds_dcurv1[j] * v);
      dM_dcurv1[1][1][j] = D * (phi_j * dv_dt + d_dK1_dt_dcurv1[j] * v);
      break;
    }
  }

} /* End of shell_moment_tensor */

/****************************************************************************/
void lubrication_shell_initialize(int *n_dof,        // Degrees of freedom
                                  int *dof_map,      // Map of DOFs
                                  int id_side,       // Side ID
                                  double xi[DIM],    // Local STU coordinates
                                  const Exo_DB *exo, // Exodus database
                                  int use_def        // Use deformed normal anyway
                                  )
/******************************************************************************
 *
 *
 *
 *
 *    Routine to set up all of the necessary shell normals and heights for
 *    lubrication shells.  Ideally, anything that needs to be done in multiple
 *    locations will be put here.  This will include loading of the heights
 *    and velocities, proper calculations of FSI and any mesh derivatives.
 *
 * Scott A Roberts (1514) sarober@sandia.gov
 *
 ******************************************************************************/
{
  int a, b, i, j, k, p;
  int el1 = -1, el2 = -1, nf;
  int FSIModel = -1;
  int n_dofptr[MAX_VARIABLE_TYPES][MDE];
  dbl wt = fv->wt;

  int ShapeVar = pd->ShapeVar;
  int mdof = ei[pg->imtrx]->dof[ShapeVar];
  int edim = ei[pg->imtrx]->ielem_dim;
  int pdim = pd->Num_Dim;
  int node, index;
  double Jlocal[DIM][DIM];
  double T[DIM - 1][DIM], t[DIM - 1][DIM]; /* t = J . T */
  int siz;
  double nx, ny, nz;
  double r_det;

  /*** LOOK FOR NEIGHBOR ELEMENTS AND SELECT CORRECT MODEL ********************/

  /* Find friends */
  el1 = ei[pg->imtrx]->ielem;
  nf = num_elem_friends[el1];

  /* Deal with number of friends */
  switch (nf) {
  case 0:
    if ((mp->FSIModel != FSI_SHELL_ONLY) && (mp->FSIModel != FSI_SHELL_ONLY_MESH) &&
        (mp->FSIModel != FSI_SHELL_ONLY_UNDEF))
      GOMA_EH(GOMA_ERROR, "ERROR:  What happened to my little friend?");
    FSIModel = mp->FSIModel;
    break;
  case 1:
    el2 = elem_friends[el1][0];
    if (mp->FSIModel == 0) {
      a = find_elemblock_index(el2, exo);
      b = Matilda[a];
      FSIModel = mp_glob[b]->FSIModel;
    } else {
      FSIModel = mp->FSIModel;
    }
    break;
  default:
    GOMA_EH(GOMA_ERROR, "ERROR: Not set up for more than one element friend!");
    break;
  }

  /* Make sure we have a model */
  if (FSIModel == 0)
    GOMA_EH(GOMA_ERROR, "ERROR: Could not find FSI Model");

  /* Use deformed normal for boundary */
  if ((use_def == 1) && (FSIModel == FSI_MESH_UNDEF))
    FSIModel = FSI_MESH_CONTINUUM;

  /*** CALCULATE SHELL NORMALS AND POPULATE FV AND BF STRUCTURES **************/
  switch (FSIModel) {

    /*** SHELL ONLY ***/
  case FSI_SHELL_ONLY:

    shell_determinant_and_normal(ei[pg->imtrx]->ielem, ei[pg->imtrx]->iconnect_ptr,
                                 ei[pg->imtrx]->num_local_nodes, ei[pg->imtrx]->ielem_dim, 1);
    n_dof[MESH_DISPLACEMENT1] = 0;
    n_dof[MESH_DISPLACEMENT2] = 0;
    n_dof[MESH_DISPLACEMENT3] = 0;

    /* calc_surf_tangent (ei[pg->imtrx]->ielem, ei[pg->imtrx]->iconnect_ptr,
       ei[pg->imtrx]->num_local_nodes, ei[pg->imtrx]->ielem_dim, ei[pg->imtrx]->num_local_nodes,
       &temp);
    */
    break;

  case FSI_SHELL_ONLY_MESH:

    if (!pd->e[pg->imtrx][R_MESH1])
      GOMA_EH(GOMA_ERROR, "ERROR:  FSI_SHELL_ONLY_MESH requires mesh equation turned on!");

    shell_determinant_and_normal(ei[pg->imtrx]->ielem, ei[pg->imtrx]->iconnect_ptr,
                                 ei[pg->imtrx]->num_local_nodes, ei[pg->imtrx]->ielem_dim, 1);

    /* Populate ndof array */
    n_dof[MESH_DISPLACEMENT1] = ei[pg->imtrx]->dof[MESH_DISPLACEMENT1];
    n_dof[MESH_DISPLACEMENT2] = ei[pg->imtrx]->dof[MESH_DISPLACEMENT2];
    n_dof[MESH_DISPLACEMENT3] = ei[pg->imtrx]->dof[MESH_DISPLACEMENT3];
    if (pd->e[pg->imtrx][R_SHELL_NORMAL1]) {
      n_dof[SHELL_NORMAL1] = ei[pg->imtrx]->dof[SHELL_NORMAL1];
      n_dof[SHELL_NORMAL2] = ei[pg->imtrx]->dof[SHELL_NORMAL2];
      n_dof[SHELL_NORMAL3] = ei[pg->imtrx]->dof[SHELL_NORMAL3];
    }

    /* Populate a trivial dof_map array */
    for (i = 0; i < ei[pg->imtrx]->dof[pd->ShapeVar]; i++) {
      dof_map[i] = i;
    }
    break;

  case FSI_SHELL_ONLY_UNDEF:

    if (!pd->e[pg->imtrx][R_MESH1])
      GOMA_EH(GOMA_ERROR, "ERROR:  FSI_SHELL_ONLY_UNDEF requires mesh equation turned on!");

    /* Populate surface determinants (detJ) and its sensitivities first */

    shell_determinant_and_normal(ei[pg->imtrx]->ielem, ei[pg->imtrx]->iconnect_ptr,
                                 ei[pg->imtrx]->num_local_nodes, ei[pg->imtrx]->ielem_dim, 1);

    /* Then calculate normal using the original configuration */

    /* Calculate Jacobian of transformation */
    for (i = 0; i < edim; i++) {
      for (j = 0; j < pdim; j++) {
        Jlocal[i][j] = 0.0;
        for (k = 0; k < mdof; k++) {
          node = ei[pg->imtrx]->dof_list[ShapeVar][k];
          index = Proc_Elem_Connect[Proc_Connect_Ptr[ei[pg->imtrx]->ielem] + node];
          Jlocal[i][j] += Coor[j][index] * bf[ShapeVar]->dphidxi[k][i];
        }
        Jlocal[2][j] = (j + 1) * 1.0;
      }
    }

#if 0
    /* Big T calculation*/
    T[0][0] = 1.;
    T[0][1] = 0.;
    T[0][2] = 0.;
    T[1][0] = 0.;
    T[1][1] = 1.;
    T[1][2] = 0.;

#else
    siz = (DIM - 1) * DIM * sizeof(double);
    memset(T, 0, siz);
    memset(t, 0, siz);
    /*  since T & t are zeroed, only need to set nonzero elements */
    /* revert back by setting id_side = 6  */
    switch (ei[pg->imtrx]->ielem_shape) {
    case SHELL:
    case TRISHELL:
      switch (id_side) {
      case 1:
        T[0][0] = 1.;
        T[1][2] = 1.;
        break;
      case 2:
        T[0][1] = 1.;
        T[1][2] = 1.;
        break;
      case 3:
        T[0][0] = -1.;
        T[1][2] = 1.;
        break;
      case 4:
        T[0][1] = -1.;
        T[1][2] = 1.;
        break;
      case 5:
        T[0][0] = 1.;
        T[1][1] = -1.;
        break;
      case 6:
      case -1:
        T[0][0] = 1.;
        T[1][1] = 1.;
        break;
      default:
        GOMA_EH(GOMA_ERROR, "Incorrect side for SHELL");
        break;
      }
      break;
    }
#endif
    /* Little t calculation */
    for (p = 0; p < 2; p++) {
      for (a = 0; a < pd->Num_Dim; a++) {
        t[p][a] = 0.;
        for (b = 0; b < pd->Num_Dim; b++) {
          t[p][a] += Jlocal[b][a] * T[p][b] * fv->h[b];
        }
      }
    }

    /* N calculation */
    nx = t[0][1] * t[1][2] - t[0][2] * t[1][1];
    ny = t[0][2] * t[1][0] - t[0][0] * t[1][2];
    nz = t[0][0] * t[1][1] - t[0][1] * t[1][0];

    /* Dets */
    r_det = 1. / sqrt(nx * nx + ny * ny + nz * nz);

    /* Overwrite normals */
    fv->snormal[0] = r_det * nx;
    fv->snormal[1] = r_det * ny;
    fv->snormal[2] = r_det * nz;

    /* Populate ndof array */
    n_dof[MESH_DISPLACEMENT1] = ei[pg->imtrx]->dof[MESH_DISPLACEMENT1];
    n_dof[MESH_DISPLACEMENT2] = ei[pg->imtrx]->dof[MESH_DISPLACEMENT2];
    n_dof[MESH_DISPLACEMENT3] = ei[pg->imtrx]->dof[MESH_DISPLACEMENT3];

    /* Populate a trivial dof_map array */
    for (i = 0; i < ei[pg->imtrx]->dof[pd->ShapeVar]; i++) {
      dof_map[i] = i;
    }

    /* Zero out sensitivities */
    int jk;
    for (a = 0; a < pdim; a++) {
      for (b = 0; b < pdim; b++) {
        for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
          jk = dof_map[k];
          fv->dsnormal_dx[a][b][jk] = 0.0;
        }
      }
    }

    break;

    /*** CONTINUUM DEFORMED MESH ***/
  case FSI_MESH_CONTINUUM:
  case FSI_MESH_ONEWAY:
  case FSI_REALSOLID_CONTINUUM:

    load_neighbor_var_data(el1, el2, n_dof, dof_map, n_dofptr, id_side, xi, exo);

    break;
    /*** CONTINUUM DEFORMED MESH ***/
  case FSI_MESH_UNDEF:

    // Load DOF count and repopulate fv structure from neighbor element
    load_neighbor_var_data(el1, el2, n_dof, dof_map, n_dofptr, id_side, xi, exo);

    // Calculate Jacobian of transformation
    for (i = 0; i < edim; i++) {
      for (j = 0; j < pdim; j++) {
        Jlocal[i][j] = 0.0;
        for (k = 0; k < mdof; k++) {
          node = ei[pg->imtrx]->dof_list[ShapeVar][k];
          index = Proc_Elem_Connect[Proc_Connect_Ptr[ei[pg->imtrx]->ielem] + node];
          Jlocal[i][j] += Coor[j][index] * bf[ShapeVar]->dphidxi[k][i];
        }
        Jlocal[2][j] = (j + 1) * 1.0;
      }
    }

    // Bit T calculation
    T[0][0] = 1.;
    T[0][1] = 0.;
    T[0][2] = 0.;
    T[1][0] = 0.;
    T[1][1] = 1.;
    T[1][2] = 0.;

    // Little t calculation
    for (p = 0; p < 2; p++) {
      for (a = 0; a < pd->Num_Dim; a++) {
        t[p][a] = 0.;
        for (b = 0; b < pd->Num_Dim; b++) {
          t[p][a] += Jlocal[b][a] * T[p][b] * fv->h[b];
        }
      }
    }

    // N calculation
    nx = t[0][1] * t[1][2] - t[0][2] * t[1][1];
    ny = t[0][2] * t[1][0] - t[0][0] * t[1][2];
    nz = t[0][0] * t[1][1] - t[0][1] * t[1][0];

    // Dets
    fv->sdet = sqrt(nx * nx + ny * ny + nz * nz);
    r_det = 1. / fv->sdet;

    // Calculate normals
    fv->snormal[0] = r_det * nx;
    fv->snormal[1] = r_det * ny;
    fv->snormal[2] = r_det * nz;

    // Zero out sensitivities
    int ldof;
    for (i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
      ldof = ei[pg->imtrx]->ln_to_dof[ShapeVar][2];
      for (a = 0; a < pdim; a++) {
        fv->dsurfdet_dx[a][ldof] = 0.0;
        for (j = 0; j < 3; j++) {
          fv->dsnormal_dx[j][a][ldof] = 0.0;
        }
      }
    }

    break;
    /*** UNIMPLEMENTED METHODS ***/
  case FSI_MESH_BOTH:
  case FSI_MESH_SHELL:

    GOMA_EH(GOMA_ERROR, "ERROR:  lubrication_shell_initialize() - Error in FSI Model");
    break;
  }

  /*** Reset weight ***/
  fv->wt = wt;
  return;
} /* End of lubrication_shell_initialize */

void Inn(double v[DIM], // Input vector
         double w[DIM]  // Output rotated vector
         )
/******************************************************************************
 *
 * Inn()
 *
 * Function to rotate into shell cordinates by the following transformation:
 *       w = (I-nn)*v
 *
 * Scott A Roberts (1514) sarober@sandia.gov
 *
 ******************************************************************************/
{
  int i, j;
  for (i = 0; i < DIM; i++) {
    w[i] = 0.0;
    for (j = 0; j < DIM; j++) {
      w[i] += (v[j] * delta(i, j) - v[j] * fv->snormal[i] * fv->snormal[j]);
    }
  }
  return;
} /* End of Inn */

void ShellRotate(double v[DIM],                // Input vector
                 double dv_dmx[DIM][DIM][MDE], // Input vector mesh sensitivity
                 double w[DIM],                // Rotated output vector
                 double dw_dmx[DIM][DIM][MDE], // Rotated output vector mesh sensitivity
                 int ndof                      // Number of DOFs for mesh equations
                 )
/******************************************************************************
 *
 * ShellRotate()
 *
 * Function to rotate into shell cordinates by the following transformation:
 *       w = (I-nn)*v
 * Also calculates the mesh sensitivities for the vector.
 *
 * Scott A Roberts (1514) sarober@sandia.gov
 *
 ******************************************************************************/
{

  /* Define variables */
  int i, j, k, l;

  /* Initialize vectors */
  memset(w, 0.0, sizeof(double) * DIM);
  memset(dw_dmx, 0.0, sizeof(double) * DIM * DIM * MDE);
  int ldim = pd->Num_Dim;
  /* Rotate vector and calculate mesh sensitivity */
  for (i = 0; i < ldim; i++) {
    for (j = 0; j < ldim; j++) {
      w[i] += v[j] * delta(i, j);
      w[i] -= v[j] * fv->snormal[i] * fv->snormal[j];
      for (k = 0; k < ldim; k++) {
        for (l = 0; l < ndof; l++) {
          dw_dmx[i][k][l] += dv_dmx[j][k][l] * delta(i, j);
          dw_dmx[i][k][l] -= dv_dmx[j][k][l] * fv->snormal[i] * fv->snormal[j];
          dw_dmx[i][k][l] -= v[j] * fv->dsnormal_dx[i][k][l] * fv->snormal[j];
          dw_dmx[i][k][l] -= v[j] * fv->snormal[i] * fv->dsnormal_dx[j][k][l];
        }
      }
    }
  }

} /* End of ShellRotate() */

void calculate_lub_q_v(const int EQN, double time, double dt, double xi[DIM], const Exo_DB *exo)
/******************************************************************************
 *
 * calculate_lub_q_v()
 *
 * Function to calculate flow rate per unit width (q) and average velocity (v)
 * in lubrication flow
 *
 *
 * Kris Tjiptowidjojo tjiptowi@unm.edu
 *
 * EDITED:
 * 2010-11-30: Scott Roberts - sarober@sandia.gov
 *    Re-wrote the lubrication section to include more key physics.
 *    Added calculation of full Jacobian entries.
 * REVAMPED to allow for shear-thinning and other viscosity models
 *    Robert Secor - rbs@hirdeal.com - October 4, 2023
 ******************************************************************************/
{
  int i, j, k, jk, w;
  dbl H;
  dbl veloL[DIM], veloU[DIM];
  dbl mu, dmu_dc = 0., dmu_dT = 0., srate = 0.;
  dbl *dmu_df;
  dbl rho;
  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct; /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;
  int VAR;
  int err;

  /* Problem dimensions */
  int dim = pd->Num_Dim;

  /* Calculate flow rate and average velocity with their sensitivities
   * depending on the lubrication model employed
   */

  /* Confined lubrication flow - Newtonian */
  /* The next else block is for film flow) */
  VAR = FILL;
  dmu_df = d_mu->F;
  if ((EQN == R_LUBP) || (EQN == R_LUBP_2)) {

    /* Set proper fill variable first.   If in lub_p layer, then use FILL,
     * but if in LUBP_2 layer use PHASE1.  This will have to be made more general
     * if you wanted to do both LS phases and R_phaseN fields in same layer.
     * We will leave that to the next sucker to develop
     */

    if (EQN == R_LUBP_2) {
      VAR = PHASE1;
      dmu_df = d_mu->pf[0];
    }

    /***** INITIALIZE LUBRICATION COMPONENTS AND LOAD IN VARIABLES *****/

    /* Setup lubrication shell constructs */
    dbl wt_old = fv->wt;
    int dof_map[MDE];
    int *n_dof = (int *)array_alloc(1, MAX_VARIABLE_TYPES, sizeof(int));
    lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

    /* Load viscosity and density */
    rho = density(d_rho, time);
    if (gn->ConstitutiveEquation == POWER_LAW || gn->ConstitutiveEquation == BINGHAM ||
        gn->ConstitutiveEquation == HERSCHEL_BULKLEY || gn->ConstitutiveEquation == CARREAU ||
        gn->ConstitutiveEquation == CARREAU_WLF || gn->ConstitutiveEquation == BINGHAM_WLF) {
      mu = gn->mu0;
    } else {
      mu = viscosity(gn, NULL, d_mu);
      dmu_dc = mp->d_viscosity[SHELL_PARTC];
      if (gn->ConstitutiveEquation == THERMAL || gn->ConstitutiveEquation == TABLE) {
        dmu_dT = mp->d_viscosity[SHELL_TEMPERATURE];
      }
    }

    /* Extract wall velocities */
    velocity_function_model(veloU, veloL, time, dt);

    /* Extract wall heights */
    dbl H_U, dH_U_dtime, H_L, dH_L_dtime;
    dbl dH_U_dX[DIM], dH_L_dX[DIM], dH_U_dp, dH_U_ddh, dH_dF[MDE];
    H = height_function_model(&H_U, &dH_U_dtime, &H_L, &dH_L_dtime, dH_U_dX, dH_L_dX, &dH_U_dp,
                              &dH_U_ddh, dH_dF, time, dt);

    /***** DEFORM HEIGHT AND CALCULATE SENSITIVITIES *****/

    /* Define variables */
    dbl D_H_DX[DIM][MDE], D_H_DP, D_H_ddh;
    dbl D_H_DRS[DIM][MDE];
    dbl D_H_DNORMAL[DIM][MDE];
    memset(D_H_DX, 0.0, sizeof(double) * DIM * MDE);
    memset(D_H_DRS, 0.0, sizeof(double) * DIM * MDE);
    memset(D_H_DNORMAL, 0.0, sizeof(double) * DIM * MDE);

    /* Deform height */
    switch (mp->FSIModel) {
    case FSI_MESH_CONTINUUM:
    case FSI_MESH_UNDEF:
    case FSI_SHELL_ONLY_UNDEF:
      for (i = 0; i < dim; i++) {
        H -= fv->snormal[i] * fv->d[i];
      }
      break;
    case FSI_SHELL_ONLY_MESH:
      if (pd->e[pg->imtrx][R_SHELL_NORMAL1] && pd->e[pg->imtrx][R_SHELL_NORMAL2] &&
          pd->e[pg->imtrx][R_SHELL_NORMAL3]) {
        for (i = 0; i < dim; i++) {
          H -= fv->n[i] * fv->d[i];
        }
      } else {
        for (i = 0; i < dim; i++) {
          H -= fv->snormal[i] * fv->d[i];
        }
      }
      break;
    case FSI_REALSOLID_CONTINUUM:
      for (i = 0; i < dim; i++) {
        H -= fv->snormal[i] * fv->d_rs[i];
      }
      break;
    }

    /* Calculate height sensitivity to mesh */
    switch (mp->FSIModel) {
    case FSI_MESH_CONTINUUM:
    case FSI_MESH_UNDEF:
    case FSI_SHELL_ONLY_UNDEF:
      for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
          for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
            jk = dof_map[k];
            D_H_DX[j][jk] +=
                delta(i, j) * (dH_U_dX[i] - dH_L_dX[i]) * bf[MESH_DISPLACEMENT1]->phi[k];
            D_H_DX[j][jk] -= fv->dsnormal_dx[i][j][jk] * fv->d[i];
            D_H_DX[j][jk] -= fv->snormal[i] * delta(i, j) * bf[MESH_DISPLACEMENT1]->phi[k];
          }
        }
      }
      break;
    case FSI_SHELL_ONLY_MESH:
      if ((pd->e[pg->imtrx][R_SHELL_NORMAL1]) && (pd->e[pg->imtrx][R_SHELL_NORMAL2]) &&
          (pd->e[pg->imtrx][R_SHELL_NORMAL3])) {
        for (i = 0; i < dim; i++) {
          for (j = 0; j < dim; j++) {
            for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
              jk = dof_map[k];
              D_H_DX[j][jk] +=
                  delta(i, j) * (dH_U_dX[i] - dH_L_dX[i]) * bf[MESH_DISPLACEMENT1]->phi[k];
              D_H_DX[j][jk] -= fv->n[i] * delta(i, j) * bf[MESH_DISPLACEMENT1]->phi[k];
            }
          }
        }
      } else {
        for (i = 0; i < dim; i++) {
          for (j = 0; j < dim; j++) {
            for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
              jk = dof_map[k];
              D_H_DX[j][jk] +=
                  delta(i, j) * (dH_U_dX[i] - dH_L_dX[i]) * bf[MESH_DISPLACEMENT1]->phi[k];
              D_H_DX[j][jk] -= fv->dsnormal_dx[i][j][jk] * fv->d[i];
              D_H_DX[j][jk] -= fv->snormal[i] * delta(i, j) * bf[MESH_DISPLACEMENT1]->phi[k];
            }
          }
        }
      }
      break;
    case FSI_REALSOLID_CONTINUUM:
      for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
          for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
            jk = dof_map[k];
            D_H_DX[j][jk] +=
                delta(i, j) * (dH_U_dX[i] - dH_L_dX[i]) * bf[MESH_DISPLACEMENT1]->phi[k];
            D_H_DX[j][jk] -= fv->dsnormal_dx[i][j][jk] * fv->d_rs[i];
          }
          for (k = 0; k < ei[pg->imtrx]->dof[SOLID_DISPLACEMENT1]; k++) {
            jk = dof_map[k];
            D_H_DRS[j][jk] -= fv->snormal[i] * delta(i, j) * bf[SOLID_DISPLACEMENT1]->phi[jk];
          }
        }
      }
      break;
    }

    /* Calculate height sensitivity to shell normal */
    switch (mp->FSIModel) {

    case FSI_SHELL_ONLY_MESH:
      if ((pd->e[pg->imtrx][R_SHELL_NORMAL1]) && (pd->e[pg->imtrx][R_SHELL_NORMAL2]) &&
          (pd->e[pg->imtrx][R_SHELL_NORMAL3])) {
        for (i = 0; i < dim; i++) {
          for (j = 0; j < dim; j++) {
            for (k = 0; k < ei[pg->imtrx]->dof[SHELL_NORMAL1]; k++) {
              D_H_DNORMAL[j][k] -= delta(i, j) * bf[SHELL_NORMAL1]->phi[k] * fv->d[i];
            }
          }
        }
      }
      break;
    }

    /* Calculate height sensitivity to pressure */
    D_H_DP = dH_U_dp;
    D_H_ddh = dH_U_ddh;

    /* Calculate height sensitivity to melting */

    /*
     * Mesh deformation does not yet affect the slopes, nor is that
     * Jacobian in there.  Maybe it should be.
     */

    /***** CALCULATE PRESSURE GRADIENT AND SENSITIVITIES *****/

    /* Define variables */
    dbl GRADP[DIM];
    dbl D_GRADP_DX[DIM][DIM][MDE];
    dbl d_grad_lubp_dmesh[DIM][DIM][MDE];
    dbl d_grad_lubp_2_dmesh[DIM][DIM][MDE];

    if (EQN == R_LUBP) {
      memset(d_grad_lubp_dmesh, 0.0, sizeof(double) * DIM * DIM * MDE);
    } else {
      memset(d_grad_lubp_2_dmesh, 0.0, sizeof(double) * DIM * DIM * MDE);
    }
    memset(GRADP, 0.0, sizeof(double) * DIM);
    memset(D_GRADP_DX, 0.0, sizeof(double) * DIM * DIM * MDE);

    if (EQN == R_LUBP) {
      /* Rotate and calculate mesh sensitivity */
      for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
          for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
            jk = dof_map[k];
            d_grad_lubp_dmesh[i][j][jk] =
                fv->d_grad_lubp_dmesh[i][j][k]; // PRS: NEED TO DO SOMETHING HERE
          }
        }
      }
      ShellRotate(fv->grad_lubp, d_grad_lubp_dmesh, GRADP, D_GRADP_DX, n_dof[MESH_DISPLACEMENT1]);
    } else {
      /* Rotate and calculate mesh sensitivity */
      for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
          for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
            jk = dof_map[k];
            d_grad_lubp_2_dmesh[i][j][jk] = fv->d_grad_lubp_2_dmesh[i][j][k];
          }
        }
      }
      ShellRotate(fv->grad_lubp_2, d_grad_lubp_2_dmesh, GRADP, D_GRADP_DX,
                  n_dof[MESH_DISPLACEMENT1]);
    }

    /***** CALCULATE HEAVISIDE GRADIENT AND SENSITIVITIES *****/

    /* Define variables */
    dbl GRADH[DIM];
    dbl D_GRADH_DF[DIM][MDE], D_GRADH_DX[DIM][DIM][MDE];
    memset(GRADH, 0.0, sizeof(double) * DIM);
    memset(D_GRADH_DF, 0.0, sizeof(double) * DIM * MDE);
    memset(D_GRADH_DX, 0.0, sizeof(double) * DIM * DIM * MDE);

    /* Rotate and calculate mesh sensitivity */
    dbl d_grad_Hside_dmx[DIM][DIM][MDE];
    memset(d_grad_Hside_dmx, 0.0, sizeof(double) * DIM * DIM * MDE);
    if (pd->v[pg->imtrx][VAR]) {
      if (mp->Lub_Curv_NormalModel) {
        load_lsi(ls->Length_Scale);
        if (!mp->Lub_Curv_Modulation || lsi->near) {
          load_lsi_derivs();
          for (i = 0; i < dim; i++) {
            for (j = 0; j < dim; j++) {
              for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
                jk = dof_map[k];
                d_grad_Hside_dmx[i][j][jk] = lsi->d_gradHn_dmesh[i][j][k];
              }
            }
          }
          ShellRotate(lsi->gradHn, d_grad_Hside_dmx, GRADH, D_GRADH_DX, n_dof[MESH_DISPLACEMENT1]);

          /* Calculate F sensitivity */
          for (i = 0; i < dim; i++) {
            for (j = 0; j < dim; j++) {
              for (k = 0; k < ei[pg->imtrx]->dof[VAR]; k++) {
                D_GRADH_DF[i][k] += lsi->d_gradHn_dF[j][k] * delta(i, j);
                D_GRADH_DF[i][k] -= lsi->d_gradHn_dF[j][k] * fv->snormal[i] * fv->snormal[j];
              }
            }
          }
        }
      } else {
        double deltan[DIM];
        load_lsi(ls->Length_Scale);
        if (!mp->Lub_Curv_Modulation || lsi->near) {
          load_lsi_derivs();
          for (i = 0; i < dim; i++) {
            deltan[i] = lsi->delta * lsi->normal[i];
            for (j = 0; j < dim; j++) {
              for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
                jk = dof_map[k];
                d_grad_Hside_dmx[i][j][jk] = lsi->d_delta_dmesh[j][k] * lsi->normal[i] +
                                             lsi->delta * lsi->d_normal_dmesh[i][j][k];
              }
            }
          }
          ShellRotate(deltan, d_grad_Hside_dmx, GRADH, D_GRADH_DX, n_dof[MESH_DISPLACEMENT1]);

          /* Calculate F sensitivity */
          for (i = 0; i < dim; i++) {
            for (j = 0; j < dim; j++) {
              for (k = 0; k < ei[pg->imtrx]->dof[VAR]; k++) {
                D_GRADH_DF[i][k] +=
                    (lsi->delta * lsi->d_normal_dF[i][k] + lsi->d_delta_dF[k] * lsi->normal[i]) *
                    (delta(i, j) - fv->snormal[i] * fv->snormal[j]);
              }
            }
          }
        }
      }
    }

    /***** CALCULATE CURVATURE AND SENSITIVITIES *****/

    /* Define variables */
    double CURV = 0.0;
    double D_CURV_DH = 0.0;
    double D_CURV_DF[MDE], D_CURV_DX[DIM][MDE], D_CURV_DNORMAL[DIM][MDE];
    /*  setting gap for capillary conditions to be the non-wall effect value */
    double H_cap, dHc_U_dX[DIM], dHc_L_dX[DIM], D_Hc_DX[DIM][MDE];
    /* Will need to reset dH_U_dX[] too at some point */
    memset(D_CURV_DF, 0.0, sizeof(double) * MDE);
    memset(D_CURV_DX, 0.0, sizeof(double) * DIM * MDE);
    memset(D_CURV_DNORMAL, 0.0, sizeof(double) * DIM * MDE);

    if (mp->HeightUFunctionModel == WALL_DISTMOD || mp->HeightUFunctionModel == WALL_DISTURB) {
      H_cap = MAX(H_U - H_L, DBL_SEMI_SMALL);
      memset(dHc_U_dX, 0.0, sizeof(double) * DIM);
      memset(dHc_L_dX, 0.0, sizeof(double) * DIM);
      memset(D_Hc_DX, 0.0, sizeof(double) * DIM * MDE);
    } else {
      H_cap = H;
      for (i = 0; i < dim; i++) {
        dHc_U_dX[i] = dH_U_dX[i];
        dHc_L_dX[i] = dH_L_dX[i];
        for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
          jk = dof_map[k];
          D_Hc_DX[i][jk] = D_H_DX[i][jk];
        }
      }
    }
    /* Curvature - analytic in the "z" direction  */
    dbl dcaU, dcaL, slopeU, slopeL;
    dcaU = dcaL = slopeU = slopeL = 0;
    if (pd->v[pg->imtrx][VAR]) {
      dcaU = mp->dcaU * M_PIE / 180.0;
      dcaL = mp->dcaL * M_PIE / 180.0;
      slopeU = slopeL = 0.;
      for (i = 0; i < dim; i++) {
        slopeU += dHc_U_dX[i] * lsi->normal[i];
        slopeL += dHc_L_dX[i] * lsi->normal[i];
      }
      /*  Positive sign for convex meniscus, negative for concave meniscus,
            this sign convention is opposite of generally accepted one for curvature
          i.e., 2H = grad-dot-normal_vector vs. 2H = -grad-dot-normal_vector       */
      CURV = -(cos(dcaU + atan(slopeU)) + cos(dcaL + atan(-slopeL))) / H_cap;

      /* Curvature - numerical in planview direction */
      if (pd->e[pg->imtrx][SHELL_LUB_CURV]) {
        CURV += fv->sh_l_curv;
      }
      if (pd->e[pg->imtrx][SHELL_LUB_CURV_2]) {
        CURV += fv->sh_l_curv_2;
      }
      if (pd->e[pg->imtrx][NORMAL1]) {
        CURV += fv->div_n;
      }

      /* Sensitivity to height */
      if (!pd->e[pg->imtrx][SHELL_LUB_CURV]) {
        D_CURV_DH = (cos(dcaU + atan(slopeU)) + cos(dcaL + atan(-slopeL))) / (H_cap * H_cap);

        /* Sensitivity to level set F */
        for (i = 0; i < ei[pg->imtrx]->dof[VAR]; i++) {
          for (j = 0; j < DIM; j++) {
            D_CURV_DF[i] += sin(dcaU + atan(slopeU)) / (H_cap * (1 + slopeU * slopeU)) *
                            dHc_U_dX[j] * lsi->d_normal_dF[j][i];
            D_CURV_DF[i] += sin(dcaL + atan(-slopeL)) / (H_cap * (1 + slopeL * slopeL)) *
                            dHc_L_dX[j] * lsi->d_normal_dF[j][i];
          }
          D_CURV_DF[i] += D_CURV_DH * dH_dF[i];
        }

        /* Sensitivity to mesh */
        for (i = 0; i < dim; i++) {
          for (j = 0; j < n_dof[MESH_DISPLACEMENT1]; j++) {
            D_CURV_DX[i][j] += D_CURV_DH * D_Hc_DX[i][j];
          }
        }

        /* Sensitivity to shell normal */
        for (i = 0; i < dim; i++) {
          for (j = 0; j < ei[pg->imtrx]->dof[SHELL_NORMAL1]; j++) {
            D_CURV_DNORMAL[i][j] += D_CURV_DH * D_H_DNORMAL[i][j];
          }
        }
      }
    }

    /***** CALCULATE GRAVITY AND LORENTZ (OTHER lubmomsource) / BODY FORCE AND SENSITIVITIES *****/

    /* Define variables */
    dbl bodf[DIM], GRAV[DIM];

    dbl D_GRAV_DF[DIM][MDE], D_GRAV_DX[DIM][DIM][MDE];
    memset(bodf, 0.0, sizeof(double) * DIM);
    memset(GRAV, 0.0, sizeof(double) * DIM);
    memset(D_GRAV_DF, 0.0, sizeof(double) * DIM * MDE);
    memset(D_GRAV_DX, 0.0, sizeof(double) * DIM * DIM * MDE);

    /* Calculate and rotate body force, calculate mesh derivatives */
    dbl Bouss[DIM];
    memset(Bouss, 0.0, sizeof(double) * DIM);
    MOMENTUM_SOURCE_DEPENDENCE_STRUCT dBouss_struct; /* Body force dependence */
    MOMENTUM_SOURCE_DEPENDENCE_STRUCT *dBouss = &dBouss_struct;

    /* Calculate and rotate body force, calculate mesh derivatives */
    dbl d_bodf_dmx[DIM][DIM][MDE];
    memset(d_bodf_dmx, 0.0, sizeof(double) * DIM * DIM * MDE);
    for (i = 0; i < dim; i++) {
      bodf[i] = mp->momentum_source[i] * rho;
    }

    ShellRotate(bodf, d_bodf_dmx, GRAV, D_GRAV_DX, n_dof[MESH_DISPLACEMENT1]);

    /* Sensitivity to level set F, then rotate */
    dbl d_bodf_df[DIM][MDE];
    for (i = 0; i < dim; i++) {
      for (j = 0; j < ei[pg->imtrx]->dof[VAR]; j++) {
        d_bodf_df[i][j] = mp->momentum_source[i] * d_rho->F[j];
      }
    }
    for (k = 0; k < ei[pg->imtrx]->dof[VAR]; k++) {
      for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
          D_GRAV_DF[i][k] += d_bodf_df[j][k] * delta(i, j);
          D_GRAV_DF[i][k] -= d_bodf_df[j][k] * fv->snormal[i] * fv->snormal[j];
        }
      }
    }

    /* Sensitivity to species if buoyancy force matters */
    if (mp->MomentumSourceModel == BOUSSINESQ) {
      err = bouss_momentum_source(Bouss, dBouss, 0, TRUE);
      GOMA_EH(err, "Problems in bouss_momentum_source");
    }

    /********** PREPARE VISCOSITY DERIVATIVES **********/
    dbl D_MU_DX[DIM][MDE];
    memset(D_MU_DX, 0.0, sizeof(double) * DIM * MDE);
    for (i = 0; i < dim; i++) {
      for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
        jk = dof_map[k];
        D_MU_DX[i][jk] = d_mu->X[i][k];
      }
    }
    if (pd->v[pg->imtrx][SHELL_PARTC]) {
      dmu_dc = mp->d_viscosity[SHELL_PARTC];
    }

    /********** CALCULATE FLOW RATE AND AVERAGE VELOCITY **********/

    /* Set some coefficients */
    dbl k_turb = 12., d_k_turb_dmu = 0., d_k_turb_dH = 0.;
    dbl vsqr, q_mag = 0., v_mag = 0., tau_w, vis_w = 1., pre_delP = 0., vpre_delP;
    dbl dq_gradp = 1., dv_gradp = 1.;
    dbl dq_dH = 0., dv_dH = 0., H_inv = 1. / H;
    dbl dqmag_dF[MDE], factor, ratio = 0., q_mag2;
    dbl q[DIM], ev[DIM], pgrad, pg_cmp[DIM], dev_dpg[DIM][DIM];
    dbl v_avg[DIM], dq_dT = 0.;
    double DQ_DH[DIM];
    double D_Q_DF[DIM][MDE], D_V_DF[DIM][MDE], DGRADP_DF[DIM][MDE], DGRADP_DK = 0.;
    double DGRADP_DX[DIM][DIM][MDE], DGRADP_DNORMAL[DIM][DIM][MDE];
    double D_Q_DGRADP[DIM][DIM], D_V_DGRADP[DIM][DIM];
    double D_Q_DT[DIM], D_V_DT[DIM];
    dbl D_Q_DP2[DIM];
    dbl D_V_DP2[DIM];
    int nonmoving_model =
        (gn->ConstitutiveEquation == BINGHAM || gn->ConstitutiveEquation == BINGHAM_WLF ||
         gn->ConstitutiveEquation == CARREAU || gn->ConstitutiveEquation == CARREAU_WLF);
    int movwall_model =
        (gn->ConstitutiveEquation == BINGHAM || gn->ConstitutiveEquation == BINGHAM_WLF ||
         gn->ConstitutiveEquation == CARREAU || gn->ConstitutiveEquation == CARREAU_WLF ||
         gn->ConstitutiveEquation == POWER_LAW || gn->ConstitutiveEquation == HERSCHEL_BULKLEY);
    int movingwall = FALSE;

    memset(D_Q_DP2, 0.0, sizeof(double) * DIM);
    memset(D_V_DP2, 0.0, sizeof(double) * DIM);
    memset(D_Q_DGRADP, 0.0, sizeof(double) * DIM * DIM);
    memset(D_V_DGRADP, 0.0, sizeof(double) * DIM * DIM);
    memset(D_Q_DF, 0.0, sizeof(double) * DIM * MDE);
    memset(D_V_DF, 0.0, sizeof(double) * DIM * MDE);
    memset(DGRADP_DF, 0.0, sizeof(double) * DIM * MDE);
    memset(DGRADP_DX, 0.0, sizeof(double) * DIM * DIM * MDE);
    memset(DGRADP_DNORMAL, 0.0, sizeof(double) * DIM * DIM * MDE);

    /* Calculate flow rate and velocity */
    memset(ev, 0.0, sizeof(double) * DIM);
    memset(dev_dpg, 0.0, sizeof(double) * DIM * DIM);
    memset(dqmag_dF, 0.0, sizeof(double) * MDE);
    for (i = 0; i < dim; i++) {
      pg_cmp[i] = GRADP[i] - GRAV[i] - Bouss[i];
      if (pd->v[pg->imtrx][VAR] && !isnan(CURV)) {
        pg_cmp[i] += GRADH[i] * CURV * mp->surface_tension;
        for (j = 0; j < ei[pg->imtrx]->dof[VAR]; j++) {
          DGRADP_DF[i][j] +=
              mp->surface_tension * (GRADH[i] * D_CURV_DF[j] + CURV * D_GRADH_DF[i][j]) -
              D_GRAV_DF[i][j];
        }
        DGRADP_DK += GRADH[i] * mp->surface_tension;
      }
      for (j = 0; j < dim; j++) {
        for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
          DGRADP_DX[i][j][k] += D_GRADP_DX[i][j][k] - D_GRAV_DX[i][j][k];
          if (pd->v[pg->imtrx][VAR] && !isnan(CURV)) {
            DGRADP_DX[i][j][k] +=
                mp->surface_tension * (D_GRADH_DX[i][j][k] * CURV + GRADH[i] * D_CURV_DX[j][k]);
          }
        }
        for (k = 0; k < ei[pg->imtrx]->dof[SHELL_NORMAL1]; k++) {
          if (pd->v[pg->imtrx][VAR] && !isnan(CURV)) {
            DGRADP_DNORMAL[i][j][k] += mp->surface_tension * GRADH[i] * D_CURV_DNORMAL[j][k];
          }
        }
      }
    }
    pgrad = 0.;
    vsqr = 0.;
    for (i = 0; i < dim; i++) {
      pgrad += SQUARE(pg_cmp[i]);
      vsqr += SQUARE(veloL[i]);
      vsqr += SQUARE(veloU[i]);
    }
    movingwall = DOUBLE_NONZERO(vsqr);
    pgrad = sqrt(pgrad);
    if (pgrad > DBL_SEMI_SMALL) {
      for (i = 0; i < dim; i++) {
        ev[i] = pg_cmp[i] / pgrad;
      }
    } else {
      ev[0] = 1.;
    }
    dev_dpg[0][0] = ev[1] * ev[1] + ev[2] * ev[2];
    dev_dpg[1][1] = ev[0] * ev[0] + ev[2] * ev[2];
    dev_dpg[2][2] = ev[1] * ev[1] + ev[0] * ev[0];
    dev_dpg[0][1] = dev_dpg[1][0] = -ev[0] * ev[1];
    dev_dpg[0][2] = dev_dpg[2][0] = -ev[0] * ev[2];
    dev_dpg[1][2] = dev_dpg[2][1] = -ev[1] * ev[2];

    tau_w = 0.5 * H * pgrad;
    if (!movingwall) {
      /*  First non-Newtonian models with analytical viscosity integration */
      if (gn->ConstitutiveEquation == POWER_LAW) {
        double nexp = gn->nexp;
        k_turb = 4. * (2. + 1. / nexp);
        tau_w = MAX(tau_w, DBL_SEMI_SMALL);
        q_mag = -2. * SQUARE(H) / k_turb * pow(tau_w / mu, 1. / nexp);
        pre_delP = -CUBE(H) / (k_turb * mu) * pow(tau_w / mu, 1. / nexp - 1.);
        dq_gradp = pre_delP / nexp;
        dq_dH = (2. + 1. / nexp) * H_inv * q_mag;
        srate = pow(fabs(tau_w) / mu, 1. / nexp);
        vis_w = tau_w / srate;
      } else if (gn->ConstitutiveEquation == HERSCHEL_BULKLEY) {
        double nexp = gn->nexp, yield = gn->tau_y, eps_rate = 0.000001;
        if (tau_w - yield > 0.) {
          double f, f_c, f_term, f_termd;
          f = yield / tau_w;
          f_c = 1. - f;
          f_term = 1. - SQUARE(f) - SQUARE(f_c) / (2. * nexp + 1.) - 2. * f * f_c / (2. + nexp);
          f_termd = SQUARE(f) + f * f_c / (2 * nexp + 1.) + f * (1. - 2 * f) / (2 + nexp);
          srate = pow((tau_w - yield) / mu, 1. / nexp);
          vis_w = tau_w / srate;
          q_mag = -0.25 * SQUARE(H) * srate * f_term;
          pre_delP = 0.5 * H * q_mag / tau_w;
          dq_gradp = -0.125 * CUBE(H) * f_term / (nexp * f_c * vis_w);
          dq_gradp += -0.25 * CUBE(H) * f_termd / vis_w;
          dq_dH = (2 + 1. / nexp / f_c) * q_mag / H - 0.5 * H * srate * f_termd;
        } else {
          srate = 0.;
          vis_w = yield / eps_rate;
          q_mag = 0.;
          pre_delP = 0.;
          dq_gradp = -0.25 * CUBE(H) / vis_w * SQUARE(nexp + 1.) / nexp / (nexp + 2.);
          dq_dH = 0.;
        }
        /*  Next  non-Newtonian models with numerical viscosity integration */
      } else if (nonmoving_model) {
        if (isnan(tau_w))
          GOMA_WH(GOMA_ERROR, "Trouble, tau_w is nan...\n");
        err = lub_viscosity_integrate(tau_w, H, &q_mag, &dq_gradp, &dq_dH, &srate, &pre_delP,
                                      &vis_w, &dq_dT);
        if (isnan(srate))
          DPRINTF(stderr, "lub_srate isnan %g %g %g %g\n", tau_w, q_mag, srate, vis_w);
        if (err < 0) {
          GOMA_WH(GOMA_ERROR, "Some trouble with Numerical Lubrication...\n");
        }
      } else { /*  Newtonian type models - nonmoving wall part  */
        k_turb = 12.;
        dq_gradp = pre_delP = -CUBE(H) / (k_turb * mu);
        q_mag = pre_delP * pgrad;
        dq_dH = -3. * SQUARE(H) / (k_turb * mu) * pgrad;
        dq_dH += q_mag * (-d_k_turb_dH / k_turb);
        srate = fabs(tau_w / mu);
        vis_w = mu;
        for (j = 0; j < ei[pg->imtrx]->dof[VAR]; j++) {
          dqmag_dF[j] += q_mag * (-d_k_turb_dmu * dmu_df[j] / k_turb - dmu_df[j] / mu);
        }
        /*  Only temperature-dependent Newtonian type model is THERMAL  */
        if (gn->ConstitutiveEquation == THERMAL || gn->ConstitutiveEquation == TABLE) {
          dq_dT = q_mag * (-dmu_dT / mu);
        }
      }
      /* modulate q (stationary wall part) if level-set interface present */
      if (pd->v[pg->imtrx][VAR] && mp->mp2nd->ViscosityModel == RATIO) {
        ratio = 1. / mp->mp2nd->viscosity; /* Assuming model = RATIO for now */
        q_mag2 = q_mag * ratio;
        q_mag = ls_modulate_property(q_mag, q_mag2, ls->Length_Scale,
                                     (double)mp->mp2nd->viscositymask[0],
                                     (double)mp->mp2nd->viscositymask[1], dqmag_dF, &factor);
        factor *= (1. - ratio);
        factor += ratio;
        dq_gradp *= factor;
        dq_dH *= factor;
        dq_dT *= factor;
        pre_delP *= factor;
        vis_w /= factor;
      }
      memset(q, 0.0, sizeof(double) * DIM);
      for (i = 0; i < dim; i++) {
        q[i] += q_mag * ev[i];
      }
      v_mag = q_mag / H;
      dv_gradp = dq_gradp / H;
      dv_dH = dq_dH / H - q_mag / SQUARE(H);
      vpre_delP = pre_delP / H;
      /* Convert to more general nomenclature  */
      for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
          D_Q_DGRADP[i][j] = dq_gradp * ev[i] * ev[j] + pre_delP * dev_dpg[i][j];
        }
        DQ_DH[i] = dq_dH * ev[i];
        DQ_DH[i] += 0.5 * (veloL[i] + veloU[i]);
      }
      /* moving wall parts  */
    } else {
      if (movwall_model) { /*  non-Newtonian models with a moving wall */
        double wstrs, relax = 0.5;
        int guess = 0;
        err = lub2D_flow2D(gn, pg_cmp, ev, dev_dpg, q, D_Q_DGRADP, DQ_DH, H, &srate, veloL, veloU,
                           guess, &wstrs, relax);
        if (isnan(srate))
          DPRINTF(stderr, "lub_srate isnan %g %g %g %g\n", wstrs, q[0], q[1], srate);
        if (err < 0) {
          GOMA_WH(GOMA_ERROR, "Some trouble with Numerical Lubrication...\n");
        }
      } else { /*  moving wall part of Newtonian type models */
        k_turb = 12.;
        dq_gradp = pre_delP = -CUBE(H) / (k_turb * mu);
        q_mag = pre_delP * pgrad;
        dq_dH = -3. * SQUARE(H) / (k_turb * mu) * pgrad;
        srate = fabs(tau_w / mu);
        vis_w = mu;
        for (j = 0; j < ei[pg->imtrx]->dof[VAR]; j++) {
          dqmag_dF[j] += q_mag * (-d_k_turb_dmu * dmu_df[j] / k_turb - dmu_df[j] / mu);
        }
        /*  Only temperature-dependent Newtonian type model is THERMAL  */
        if (gn->ConstitutiveEquation == THERMAL || gn->ConstitutiveEquation == TABLE) {
          dq_dT = q_mag * (-dmu_dT / mu);
        }
        memset(q, 0.0, sizeof(double) * DIM);
        for (i = 0; i < dim; i++) {
          q[i] += q_mag * ev[i];
          q[i] += 0.5 * H * (veloL[i] + veloU[i]);
        }
        /* Convert to more general nomenclature  */
        for (i = 0; i < dim; i++) {
          for (j = 0; j < dim; j++) {
            D_Q_DGRADP[i][j] = dq_gradp * ev[i] * ev[j] + pre_delP * dev_dpg[i][j];
          }
          DQ_DH[i] = dq_dH * ev[i];
          DQ_DH[i] += 0.5 * (veloL[i] + veloU[i]);
        }
      } /*  End of Viscosity Models **/

      /* modulate q (moving wall part) if level-set interface present */
      if (pd->v[pg->imtrx][VAR] && mp->mp2nd->ViscosityModel == RATIO) {
        ratio = 1. / mp->mp2nd->viscosity; /* Assuming model = RATIO for now */
        q_mag2 = q_mag * ratio;
        q_mag = ls_modulate_property(q_mag, q_mag2, ls->Length_Scale,
                                     (double)mp->mp2nd->viscositymask[0],
                                     (double)mp->mp2nd->viscositymask[1], dqmag_dF, &factor);
        factor *= (1. - ratio);
        factor += ratio;
        /* Possibly lots different here since dq_gradp not really used for moving wall */
        dq_gradp *= factor;
        dq_dH *= factor;
        pre_delP *= factor;
        vis_w /= factor;
      }

      v_mag = q_mag / H;
      dv_gradp = dq_gradp / H;
      dv_dH = dq_dH / H - q_mag / SQUARE(H);
      vpre_delP = pre_delP / H;
      /** Will need conversion to general here too  */
    } /* End of moving wall part  */

    memset(v_avg, 0.0, sizeof(double) * DIM);
    for (i = 0; i < dim; i++) {
      v_avg[i] += q[i] / H;
    }

    /* Sensitivity w.r.t. height */
    dbl D_Q_DH[DIM] = {0.0};
    dbl D_V_DH[DIM] = {0.0};
    if (movingwall) {
      if (movwall_model) {
        for (i = 0; i < dim; i++) {
          D_Q_DH[i] += DQ_DH[i];
        }
        for (i = 0; i < dim; i++) {
          D_V_DH[i] += DQ_DH[i] / H - q[i] / SQUARE(H);
        }
      } else {
        for (i = 0; i < dim; i++) {
          D_Q_DH[i] += dq_dH * ev[i];
          D_Q_DH[i] += q_mag * (-d_k_turb_dH / k_turb) * ev[i];
          if (gn->ConstitutiveEquation == NEWTONIAN || 1)
            D_Q_DH[i] += 0.5 * (veloL[i] + veloU[i]);
        }
        for (i = 0; i < dim; i++) {
          D_V_DH[i] += dv_dH * ev[i];
          D_V_DH[i] += v_mag * (-d_k_turb_dH / k_turb) * ev[i];
        }
      }
    }

    /* Sensitivity w.r.t. pressure */

    if (movingwall && movwall_model) {
      for (i = 0; i < dim; i++) {
        D_Q_DP2[i] += D_Q_DH[i] * D_H_DP;
        D_V_DP2[i] = (D_Q_DH[i] - q[i] / H) / H * D_H_DP;
        for (j = 0; j < dim; j++) {
          D_V_DGRADP[i][j] = D_Q_DGRADP[i][j] / H;
        }
      }
    } else {
      for (i = 0; i < dim; i++) {
        D_Q_DP2[i] += D_Q_DH[i] * D_H_DP;
        D_V_DP2[i] += D_V_DH[i] * D_H_DP;
      }
      for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
          D_V_DGRADP[i][j] = dv_gradp * ev[i] * ev[j] + vpre_delP * dev_dpg[i][j];
        }
      }
    }

    /* Sensitivity w.r.t. level set */
    memset(D_Q_DF, 0.0, sizeof(double) * DIM * MDE);
    memset(D_V_DF, 0.0, sizeof(double) * DIM * MDE);
    if (pd->v[pg->imtrx][VAR]) {
      for (i = 0; i < dim; i++) {
        for (j = 0; j < ei[pg->imtrx]->dof[VAR]; j++) {
          D_Q_DF[i][j] += dqmag_dF[j] * ev[i];
          D_Q_DF[i][j] += D_Q_DH[i] * dH_dF[j];
          for (k = 0; k < dim; k++) {
            D_Q_DF[i][j] += D_Q_DGRADP[i][k] * DGRADP_DF[k][j];
          }
        }
      }
      for (i = 0; i < dim; i++) {
        for (j = 0; j < ei[pg->imtrx]->dof[VAR]; j++) {
          D_V_DF[i][j] += dqmag_dF[j] / H * ev[i];
          D_V_DF[i][j] += D_V_DH[i] * dH_dF[j];
          for (k = 0; k < dim; k++) {
            D_V_DF[i][j] += D_V_DGRADP[i][k] * DGRADP_DF[k][j];
          }
        }
      }
    }

    /* Sensitivity w.r.t. curvature */
    dbl D_Q_DK[DIM], D_V_DK[DIM];

    memset(D_Q_DK, 0.0, sizeof(double) * DIM);
    memset(D_V_DK, 0.0, sizeof(double) * DIM);
    for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
        D_Q_DK[i] += D_Q_DGRADP[i][j] * GRADH[j] * mp->surface_tension;
        D_V_DK[i] += D_V_DGRADP[i][j] * GRADH[j] * mp->surface_tension;
      }
    }

    /* Sensitivity w.r.t. mesh and/or real-solid */
    dbl D_Q_DX[DIM][DIM][MDE], D_V_DX[DIM][DIM][MDE];
    dbl D_Q_DRS[DIM][DIM][MDE], D_V_DRS[DIM][DIM][MDE];
    memset(D_Q_DX, 0.0, sizeof(double) * DIM * DIM * MDE);
    memset(D_V_DX, 0.0, sizeof(double) * DIM * DIM * MDE);
    memset(D_Q_DRS, 0.0, sizeof(double) * DIM * DIM * MDE);
    memset(D_V_DRS, 0.0, sizeof(double) * DIM * DIM * MDE);
    switch (mp->FSIModel) {
    case FSI_MESH_CONTINUUM:
    case FSI_MESH_UNDEF:
    case FSI_SHELL_ONLY_MESH:
    case FSI_SHELL_ONLY_UNDEF:
      for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
          for (k = 0; k < n_dof[MESH_DISPLACEMENT1]; k++) {
            D_Q_DX[i][j][k] += D_Q_DH[i] * D_H_DX[j][k];
            D_Q_DX[i][j][k] += q_mag * (-d_k_turb_dmu * D_MU_DX[j][k] / k_turb) * ev[i];
            D_Q_DX[i][j][k] += q_mag * (-D_MU_DX[j][k] / mu) * ev[i];
            for (jk = 0; jk < dim; jk++) {
              D_Q_DX[i][j][k] += D_Q_DGRADP[i][jk] * D_GRADP_DX[jk][j][k];
            }
          }
        }
      }
      for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
          for (k = 0; k < n_dof[MESH_DISPLACEMENT1]; k++) {
            D_V_DX[i][j][k] += D_V_DH[i] * D_H_DX[j][k];
            D_V_DX[i][j][k] += v_mag * (-d_k_turb_dmu * D_MU_DX[j][k] / k_turb) * ev[i];
            D_V_DX[i][j][k] += v_mag * (-D_MU_DX[j][k] / mu) * ev[i];
            for (jk = 0; jk < dim; jk++) {
              D_V_DX[i][j][k] += D_V_DGRADP[i][jk] * D_GRADP_DX[jk][j][k];
            }
          }
        }
      }
      break;

    case FSI_REALSOLID_CONTINUUM:
      for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
          for (k = 0; k < n_dof[MESH_DISPLACEMENT1]; k++) {
            D_Q_DX[i][j][k] += D_Q_DH[i] * D_H_DX[j][k];
            D_Q_DX[i][j][k] += q_mag * (-d_k_turb_dmu * D_MU_DX[j][k] / k_turb) * ev[i];
            D_Q_DX[i][j][k] += q_mag * (-D_MU_DX[j][k] / mu) * ev[i];
            for (jk = 0; jk < dim; jk++) {
              D_Q_DX[i][j][k] += D_Q_DGRADP[i][jk] * D_GRADP_DX[jk][j][k];
            }
            for (k = 0; k < n_dof[SOLID_DISPLACEMENT1]; k++) {
              D_Q_DRS[i][j][k] += D_Q_DH[i] * D_H_DRS[j][k];
            }
          }
        }
      }
      for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
          for (k = 0; k < n_dof[MESH_DISPLACEMENT1]; k++) {
            D_V_DX[i][j][k] += D_V_DH[i] * D_H_DX[j][k];
            D_V_DX[i][j][k] += v_mag * (-d_k_turb_dmu * D_MU_DX[j][k] / k_turb) * ev[i];
            D_V_DX[i][j][k] += v_mag * (-D_MU_DX[j][k] / mu) * ev[i];
            for (jk = 0; jk < dim; jk++) {
              D_V_DX[i][j][k] += D_V_DGRADP[i][jk] * D_GRADP_DX[jk][j][k];
            }
          }
          for (k = 0; k < n_dof[MESH_DISPLACEMENT1]; k++) {
            D_V_DRS[i][j][k] += D_V_DH[i] * D_H_DRS[j][k];
          }
        }
      }
      break;
    }

    /* Sensitivity w.r.t. shell normal */
    dbl D_Q_DNORMAL[DIM][DIM][MDE], D_V_DNORMAL[DIM][DIM][MDE];
    memset(D_Q_DNORMAL, 0.0, sizeof(double) * DIM * DIM * MDE);
    memset(D_V_DNORMAL, 0.0, sizeof(double) * DIM * DIM * MDE);
    if ((pd->v[pg->imtrx][SHELL_NORMAL1]) && (pd->v[pg->imtrx][SHELL_NORMAL1]) &&
        (pd->v[pg->imtrx][SHELL_NORMAL1])) {
      for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
          for (k = 0; k < ei[pg->imtrx]->dof[SHELL_NORMAL1]; k++) {
            D_Q_DNORMAL[i][j][k] += D_Q_DH[i] * D_H_DNORMAL[j][k];
            if (pd->v[pg->imtrx][VAR]) {
              for (jk = 0; jk < dim; jk++) {
                D_Q_DNORMAL[i][j][k] += D_Q_DGRADP[i][jk] * DGRADP_DNORMAL[jk][j][k];
              }
            }
          }
        }
      }
      for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
          for (k = 0; k < ei[pg->imtrx]->dof[SHELL_NORMAL1]; k++) {
            D_V_DNORMAL[i][j][k] += D_V_DH[i] * D_H_DNORMAL[j][k];
            if (pd->v[pg->imtrx][VAR]) {
              for (jk = 0; jk < dim; jk++) {
                D_V_DNORMAL[i][j][k] += D_V_DGRADP[i][jk] * DGRADP_DNORMAL[jk][j][k];
              }
            }
          }
        }
      }
    }

    /* Sensitivity w.r.t. dH */
    dbl D_Q_DdH[DIM], D_V_DdH[DIM];
    if (pd->v[pg->imtrx][SHELL_DELTAH] && (mp->HeightUFunctionModel == CONSTANT_SPEED_DEFORM ||
                                           mp->HeightUFunctionModel == CONSTANT_SPEED_MELT ||
                                           mp->HeightUFunctionModel == ROLL_ON_MELT ||
                                           mp->HeightUFunctionModel == FLAT_GRAD_FLAT_MELT ||
                                           mp->HeightUFunctionModel == CIRCLE_MELT)) {
      for (i = 0; i < dim; i++) {
        D_Q_DdH[i] = D_Q_DH[i] * D_H_ddh;
        D_V_DdH[i] = D_V_DH[i] * D_H_ddh;
      }
    }

    /* Sensitivity w.r.t. sh_pc , i.e., scalar viscosity dependencies */
    dbl D_Q_DC[DIM][MDE], D_V_DC[DIM][MDE];
    memset(D_Q_DC, 0.0, sizeof(double) * DIM * MDE);
    memset(D_V_DC, 0.0, sizeof(double) * DIM * MDE);
    if (pd->v[pg->imtrx][SHELL_PARTC]) {
      for (i = 0; i < dim; i++) {
        for (j = 0; j < ei[pg->imtrx]->dof[SHELL_PARTC]; j++) {
          D_Q_DC[i][j] += q_mag * (-dmu_dc / mu) * ev[i];
        }
      }
      for (i = 0; i < DIM; i++) {
        for (j = 0; j < ei[pg->imtrx]->dof[SHELL_PARTC]; j++) {
          D_V_DC[i][j] += v_mag * (-dmu_dc / mu) * ev[i];
        }
      }
    }

    /* Sensitivity w.r.t. sh_t , i.e., scalar viscosity dependencies */
    if (pd->v[pg->imtrx][SHELL_TEMPERATURE]) {
      for (i = 0; i < dim; i++) {
        D_Q_DT[i] = dq_dT * ev[i];
        D_V_DT[i] = H_inv * dq_dT * ev[i];
      }
    }

    /* Sensitivity w.r.t. species for buoyancy force */
    dbl D_Q_D_CONC[DIM][DIM][MAX_CONC][MDE], D_V_D_CONC[DIM][DIM][MAX_CONC][MDE];
    memset(D_Q_D_CONC, 0.0, sizeof(double) * DIM * DIM * MAX_CONC * MDE);
    memset(D_V_D_CONC, 0.0, sizeof(double) * DIM * DIM * MAX_CONC * MDE);
    if (pd->v[pg->imtrx][MASS_FRACTION]) {
      for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
          for (w = 0; w < pd->Num_Species_Eqn; w++) {
            for (k = 0; k < ei[pg->imtrx]->dof[MASS_FRACTION]; k++) {
              D_Q_D_CONC[i][j][w][k] += D_Q_DGRADP[i][j] * (-dBouss->C[j][w][k]);
            }
          }
        }
      }
      for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
          for (w = 0; w < pd->Num_Species_Eqn; w++) {
            for (k = 0; k < ei[pg->imtrx]->dof[MASS_FRACTION]; k++) {
              D_V_D_CONC[i][j][w][k] += D_V_DGRADP[i][j] * (-dBouss->C[j][w][k]);
            }
          }
        }
      }
    }

    /******* STORE THE INFORMATION TO LUBRICATION AUXILIARIES STRUCTURE ***********/

    memset(LubAux->dq_dx, 0.0, sizeof(double) * DIM * DIM * MDE);
    memset(LubAux->dq_dnormal, 0.0, sizeof(double) * DIM * DIM * MDE);

    LubAux->H = H;
    LubAux->dH_dp = D_H_DP;
    LubAux->dH_ddh = D_H_ddh;
    LubAux->gradP_mag = 0;
    LubAux->srate = srate;
    LubAux->mu_star = vis_w;
    for (i = 0; i < dim; i++) {
      LubAux->q[i] = q[i];
      LubAux->v_avg[i] = v_avg[i];
      LubAux->gradP_mag += SQUARE(GRADP[i] - GRAV[i] - Bouss[i]);

      LubAux->dq_dp2[i] = D_Q_DP2[i];
      LubAux->dv_avg_dp2[i] = D_V_DP2[i];
      LubAux->dq_dk[i] = D_Q_DK[i];
      LubAux->dv_avg_dk[i] = D_V_DK[i];
      LubAux->dq_ddh[i] = D_Q_DdH[i];
      LubAux->dv_avg_ddh[i] = D_V_DdH[i];
      LubAux->dq_dT[i] = D_Q_DT[i];
      LubAux->dv_avg_dT[i] = D_V_DT[i];
      for (j = 0; j < dim; j++) {
        LubAux->dq_dgradp[i][j] = D_Q_DGRADP[i][j];
        LubAux->dv_dgradp[i][j] = D_V_DGRADP[i][j];
      }
      for (j = 0; j < ei[pg->imtrx]->dof[VAR]; j++) {
        LubAux->dq_df[i][j] = D_Q_DF[i][j];
        LubAux->dv_avg_df[i][j] = D_V_DF[i][j];
      }

      for (j = 0; j < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; j++) {
        jk = dof_map[j];
        for (k = 0; k < dim; k++) {
          LubAux->dq_dx[i][k][j] = D_Q_DX[i][k][jk];
          LubAux->dv_avg_dx[i][k][j] = D_V_DX[i][k][jk];
        }
      }
      if ((pd->v[pg->imtrx][SHELL_NORMAL1]) && (pd->v[pg->imtrx][SHELL_NORMAL2]) &&
          (pd->v[pg->imtrx][SHELL_NORMAL3])) {
        for (j = 0; j < ei[pg->imtrx]->dof[SHELL_NORMAL1]; j++) {
          for (k = 0; k < dim; k++) {
            LubAux->dq_dnormal[i][k][j] = D_Q_DNORMAL[i][k][j];
            LubAux->dv_avg_dnormal[i][k][j] = D_V_DNORMAL[i][k][j];
          }
        }
      }
      for (j = 0; j < ei[pg->imtrx]->dof[SOLID_DISPLACEMENT1]; j++) {
        jk = dof_map[j];
        for (k = 0; k < dim; k++) {
          LubAux->dq_drs[i][k][j] = D_Q_DRS[i][k][jk];
          LubAux->dv_avg_drs[i][k][j] = D_V_DRS[i][k][jk];
        }
      }
      for (j = 0; j < ei[pg->imtrx]->dof[SHELL_PARTC]; j++) {
        LubAux->dq_dc[i][j] = D_Q_DC[i][j];
        LubAux->dv_avg_dc[i][j] = D_V_DC[i][j];
      }
      if ((pd->v[pg->imtrx][MASS_FRACTION])) {
        for (j = 0; j < dim; j++) {
          for (w = 0; w < pd->Num_Species_Eqn; w++) {
            for (k = 0; k < ei[pg->imtrx]->dof[MASS_FRACTION]; k++) {
              LubAux->dq_dconc[i][j][w][k] = D_Q_D_CONC[i][j][w][k];
              LubAux->dv_avg_dconc[i][j][w][k] = D_V_D_CONC[i][j][w][k];
            }
          }
        }
      }
    }
    LubAux->gradP_mag = sqrt(LubAux->gradP_mag);

    for (j = 0; j < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; j++) {
      jk = dof_map[j];
      for (k = 0; k < dim; k++) {
        LubAux->dH_dmesh[k][j] = D_H_DX[k][jk];
        LubAux->dH_drealsolid[k][j] = D_H_DRS[k][jk];
      }
    }

    // Cleanup
    fv->wt = wt_old;
    safe_free((void *)n_dof);

  } else if (EQN == R_SHELL_FILMP) {

    /******* PRECALCULATE ALL NECESSARY COMPONENTS ***********/

    /* Setup lubrication shell constructs */
    dbl wt_old = fv->wt;
    int dof_map[MDE];
    int *n_dof = (int *)array_alloc(1, MAX_VARIABLE_TYPES, sizeof(int));
    lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

    /* Load viscosity */
    if (gn->ConstitutiveEquation == POWER_LAW || gn->ConstitutiveEquation == BINGHAM ||
        gn->ConstitutiveEquation == HERSCHEL_BULKLEY || gn->ConstitutiveEquation == CARREAU ||
        gn->ConstitutiveEquation == CARREAU_WLF || gn->ConstitutiveEquation == BINGHAM_WLF) {
      mu = gn->mu0;
    } else {
      mu = viscosity(gn, NULL, d_mu);
      dmu_dc = mp->d_viscosity[SHELL_PARTC];
      if (gn->ConstitutiveEquation == THERMAL || gn->ConstitutiveEquation == TABLE) {
        dmu_dT = mp->d_viscosity[SHELL_TEMPERATURE];
      }
    }

    /* Extract bottom wall velocity */
    velocity_function_model(veloU, veloL, time, dt);

    /* Get slip coefficient */
    double beta_slip;
    beta_slip = mp->SlipCoeff;

    /***** CALCULATE HEIGHT, SLOPES, AND SENSITIVITIES *****/

    /* Get lower height from height function model */
    dbl H_U, dH_U_dtime, H_L, dH_L_dtime;
    dbl dH_U_dX[DIM], dH_L_dX[DIM], dH_U_dp, dH_U_ddh, dH_dF[MDE];
    H = height_function_model(&H_U, &dH_U_dtime, &H_L, &dH_L_dtime, dH_U_dX, dH_L_dX, &dH_U_dp,
                              &dH_U_ddh, dH_dF, time, dt);

    /* Get the net film thickness */
    dbl H = fv->sh_fh - H_L;

    /* Define variables */
    dbl D_H_DX[DIM][MDE], D_H_DSH_FH[MDE];
    dbl D_H_DRS[DIM][MDE];
    dbl D_H_DNORMAL[DIM][MDE];
    memset(D_H_DX, 0.0, sizeof(double) * DIM * MDE);
    memset(D_H_DRS, 0.0, sizeof(double) * DIM * MDE);
    memset(D_H_DNORMAL, 0.0, sizeof(double) * DIM * MDE);
    memset(D_H_DSH_FH, 0.0, sizeof(double) * MDE);

    /* Deform height */
    switch (mp->FSIModel) {
    case FSI_MESH_CONTINUUM:
    case FSI_MESH_UNDEF:
    case FSI_SHELL_ONLY_UNDEF:
      for (i = 0; i < dim; i++) {
        H -= fv->snormal[i] * fv->d[i];
      }
      break;
    case FSI_SHELL_ONLY_MESH:
      if (pd->e[pg->imtrx][R_SHELL_NORMAL1] && pd->e[pg->imtrx][R_SHELL_NORMAL2] &&
          pd->e[pg->imtrx][R_SHELL_NORMAL3]) {
        for (i = 0; i < dim; i++) {
          H -= fv->n[i] * fv->d[i];
        }
      } else {
        for (i = 0; i < dim; i++) {
          H -= fv->snormal[i] * fv->d[i];
        }
      }
      break;
    case FSI_REALSOLID_CONTINUUM:
      for (i = 0; i < dim; i++) {
        H -= fv->snormal[i] * fv->d_rs[i];
      }
      break;
    }

    /* Calculate height sensitivity to mesh */
    switch (mp->FSIModel) {
    case FSI_MESH_CONTINUUM:
    case FSI_MESH_UNDEF:
    case FSI_SHELL_ONLY_UNDEF:
      for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
          for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
            jk = dof_map[k];
            D_H_DX[j][jk] += delta(i, j) * (-dH_L_dX[i]) * bf[MESH_DISPLACEMENT1]->phi[k];
            D_H_DX[j][jk] -= fv->dsnormal_dx[i][j][jk] * fv->d[i];
            D_H_DX[j][jk] -= fv->snormal[i] * delta(i, j) * bf[MESH_DISPLACEMENT1]->phi[k];
          }
        }
      }
      break;
    case FSI_SHELL_ONLY_MESH:
      if ((pd->e[pg->imtrx][R_SHELL_NORMAL1]) && (pd->e[pg->imtrx][R_SHELL_NORMAL2]) &&
          (pd->e[pg->imtrx][R_SHELL_NORMAL3])) {
        for (i = 0; i < dim; i++) {
          for (j = 0; j < dim; j++) {
            for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
              jk = dof_map[k];
              D_H_DX[j][jk] += delta(i, j) * (-dH_L_dX[i]) * bf[MESH_DISPLACEMENT1]->phi[k];
              D_H_DX[j][jk] -= fv->n[i] * delta(i, j) * bf[MESH_DISPLACEMENT1]->phi[k];
            }
          }
        }
      } else {
        for (i = 0; i < dim; i++) {
          for (j = 0; j < dim; j++) {
            for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
              jk = dof_map[k];
              D_H_DX[j][jk] += delta(i, j) * (-dH_L_dX[i]) * bf[MESH_DISPLACEMENT1]->phi[k];
              D_H_DX[j][jk] -= fv->dsnormal_dx[i][j][jk] * fv->d[i];
              D_H_DX[j][jk] -= fv->snormal[i] * delta(i, j) * bf[MESH_DISPLACEMENT1]->phi[k];
            }
          }
        }
      }
      break;
    }

    /* Calculate height sensitivity to shell normal */
    switch (mp->FSIModel) {
    case FSI_SHELL_ONLY_MESH:
      if ((pd->e[pg->imtrx][R_SHELL_NORMAL1]) && (pd->e[pg->imtrx][R_SHELL_NORMAL2]) &&
          (pd->e[pg->imtrx][R_SHELL_NORMAL3])) {
        for (i = 0; i < dim; i++) {
          for (j = 0; j < dim; j++) {
            for (k = 0; k < ei[pg->imtrx]->dof[SHELL_NORMAL1]; k++) {
              D_H_DNORMAL[j][k] -= delta(i, j) * bf[SHELL_NORMAL1]->phi[k] * fv->d[i];
            }
          }
        }
      }
      break;
    }

    /* Calculate height sensitivity to shell height */

    for (j = 0; j < ei[pg->imtrx]->dof[SHELL_FILMH]; j++) {
      D_H_DSH_FH[j] = 1.0;
    }

    /* Calculate slopes and sensitivities */

    dbl grad_H[DIM];
    dbl GRADH[DIM];
    dbl D_GRADH_DSH_FH[DIM][MDE], D_GRADH_DX[DIM][DIM][MDE];
    dbl d_grad_H_dmesh[DIM][DIM][MDE];

    memset(grad_H, 0.0, sizeof(double) * DIM);
    for (i = 0; i < dim; i++) {
      grad_H[i] = fv->grad_sh_fh - dH_L_dX; /* Net film thickness gradient */
    }

    memset(d_grad_H_dmesh, 0.0, sizeof(double) * DIM * DIM * MDE);
    memset(GRADH, 0.0, sizeof(double) * DIM);
    memset(D_GRADH_DSH_FH, 0.0, sizeof(double) * DIM * MDE);
    memset(D_GRADH_DX, 0.0, sizeof(double) * DIM * DIM * MDE);

    /* Rotate and calculate mesh sensitivity */
    for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
        for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
          jk = dof_map[k];
          d_grad_H_dmesh[i][j][jk] =
              fv->d_grad_sh_fh_dmesh[i][j][k]; // Ignore dH_L_dX dependency on mesh for now
        }
      }
    }
    ShellRotate(grad_H, d_grad_H_dmesh, GRADH, D_GRADH_DX, n_dof[MESH_DISPLACEMENT1]);

    /* Calculate sensitivity w.r.t. shell height */
    for (j = 0; j < dim; j++) {
      for (k = 0; k < ei[pg->imtrx]->dof[SHELL_FILMH]; k++) {
        D_GRADH_DSH_FH[j][k] = 1.0;
      }
    }

    /***** CALCULATE GRAVITY AND LORENTZ (OTHER lubmomsource) / BODY FORCE AND SENSITIVITIES
     * *****/

    dbl GRAV[DIM];
    memset(GRAV, 0.0, sizeof(double) * DIM);

    /* Calculate and rotate body force, calculate mesh derivatives */
    for (i = 0; i < DIM; i++) {
      GRAV[i] = mp->momentum_source[i];
    }

    /***** CALCULATE PRESSURE GRADIENT AND SENSITIVITIES *****/

    /* Define variables */
    dbl GRADP[DIM];
    dbl d_grad_P_dmesh[DIM][DIM][MDE];
    dbl D_GRADP_DX[DIM][DIM][MDE];

    memset(GRADP, 0.0, sizeof(double) * DIM);

    memset(d_grad_P_dmesh, 0.0, sizeof(double) * DIM * DIM * MDE);

    /* Rotate and calculate mesh sensitivity */
    for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
        for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
          jk = dof_map[k];
          d_grad_P_dmesh[i][j][jk] = fv->d_grad_sh_fp_dmesh[i][j][k];
        }
      }
    }

    ShellRotate(fv->grad_sh_fp, d_grad_P_dmesh, GRADP, D_GRADP_DX, n_dof[MESH_DISPLACEMENT1]);

    /***** CALCULATE DISJOINING PRESSURE GRADIENT AND SENSITIVITIES *****/

    /* Define variables */

    dbl GRAD_DISJ_PRESS[DIM];
    dbl D_GRAD_DISJ_PRESS_DH1[DIM][MDE], D_GRAD_DISJ_PRESS_DH2[DIM][MDE];
    dbl D_GRAD_DISJ_PRESS_DH[DIM][MDE];
    memset(GRAD_DISJ_PRESS, 0.0, sizeof(double) * DIM);
    memset(D_GRAD_DISJ_PRESS_DH1, 0.0, sizeof(double) * DIM * MDE);
    memset(D_GRAD_DISJ_PRESS_DH2, 0.0, sizeof(double) * DIM * MDE);

    /* Evaluate disjoining pressure and its sensitivities */
    disjoining_pressure_model(fv->sh_fh, fv->grad_sh_fh, n_dof, dof_map, GRAD_DISJ_PRESS,
                              D_GRAD_DISJ_PRESS_DH1, D_GRAD_DISJ_PRESS_DH2, D_GRAD_DISJ_PRESS_DH);
    /******* CALCULATE FLOW RATE AND AVERAGE VELOCITY ***********/

    double q[DIM], v_avg[DIM], pg_cmp[DIM], pgrad, ev[DIM], dev_dpg[DIM][DIM];
    double DGRADP_DX[DIM][DIM][MDE];
    double H_inv = 1. / H;
    memset(q, 0.0, sizeof(double) * DIM);
    memset(v_avg, 0.0, sizeof(double) * DIM);
    dbl q_mag = 0., pre_delP = 0., dq_gradp = 1., vpre_delP = 1.;
    dbl k_turb = 3., tau_w, vis_w, dq_gradpt, dv_gradp;
    dbl vsqr, dq_dH = 0., dv_dH = 0., dq_dT = 0.;

    for (i = 0; i < dim; i++) {
      pg_cmp[i] = GRADP[i] - GRAV[i] - GRAD_DISJ_PRESS[i];
      for (j = 0; j < dim; j++) {
        for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
          DGRADP_DX[i][j][k] += D_GRADP_DX[i][j][k];
        }
      }
    }
    pgrad = 0.;
    vsqr = 0.;
    for (i = 0; i < dim; i++) {
      pgrad += SQUARE(pg_cmp[i]);
      vsqr += SQUARE(veloL[i]);
    }
    pgrad = sqrt(pgrad);
    if (DOUBLE_NONZERO(pgrad)) {
      for (i = 0; i < dim; i++) {
        ev[i] = pg_cmp[i] / pgrad;
      }
    } else {
      ev[0] = 1.;
    }
    dev_dpg[0][0] = ev[1] * ev[1] + ev[2] * ev[2];
    dev_dpg[1][1] = ev[0] * ev[0] + ev[2] * ev[2];
    dev_dpg[2][2] = ev[1] * ev[1] + ev[0] * ev[0];
    dev_dpg[0][1] = dev_dpg[1][0] = -ev[0] * ev[1];
    dev_dpg[0][2] = dev_dpg[2][0] = -ev[0] * ev[2];
    dev_dpg[1][2] = dev_dpg[2][1] = -ev[1] * ev[2];
    /* Evaluate flow rate and average velocity */
    tau_w = H * pgrad;
    if (gn->ConstitutiveEquation == POWER_LAW) {
      double nexp = gn->nexp;
      k_turb = 2. + 1. / nexp;
      if (DOUBLE_NONZERO(pgrad)) {
        q_mag = -SQUARE(H) / k_turb * pow(tau_w / mu, 1. / nexp);
        pre_delP = -CUBE(H) / (k_turb * mu) * pow(tau_w / mu, 1. / nexp - 1.);
        dq_gradp = pre_delP / nexp;
        dq_dH = (2. + 1. / nexp) * H_inv * q_mag;
        srate = pow(fabs(tau_w) / mu, 1. / nexp);
        vis_w = tau_w / srate;
      } else {
        q_mag = dq_dH = 0.;
        dq_gradp = pre_delP = -CUBE(H) / (k_turb * mu);
        srate = 0.;
        vis_w = mu;
      }
    } else if (gn->ConstitutiveEquation == BINGHAM || gn->ConstitutiveEquation == BINGHAM_WLF ||
               gn->ConstitutiveEquation == CARREAU || gn->ConstitutiveEquation == CARREAU_WLF) {
      if (isnan(tau_w))
        GOMA_WH(GOMA_ERROR, "Trouble, tau_w is nan...\n");
      err = lub_viscosity_integrate(tau_w, H, &q_mag, &dq_gradp, &dq_dH, &srate, &pre_delP, &vis_w,
                                    &dq_dT);
      if (isnan(srate))
        fprintf(stderr, "srate isnan %g %g %g %g\n", tau_w, q_mag, srate, vis_w);
      if (err < 0) {
        GOMA_WH(GOMA_ERROR, "Some trouble with Numerical Lubrication...\n");
      }

      /** Make corrections for film flow from confined calculations **/
      q_mag *= 2.;
      pre_delP *= 4.;
      dq_gradp *= 4.;
      dq_dH *= 2.;
      q_mag -= beta_slip * SQUARE(H) * pgrad;
      pre_delP -= beta_slip * SQUARE(H);
      dq_dH -= beta_slip * 2. * H * pgrad;
    } else if (gn->ConstitutiveEquation == HERSCHEL_BULKLEY) {
      double nexp = gn->nexp, yield = gn->tau_y;
      k_turb = (2 * nexp + 1.);
      if (tau_w > yield) {
        pre_delP = -CUBE(H) / k_turb * pow((tau_w - yield) / mu, 1. / nexp) *
                   (nexp / tau_w - nexp / (1. + nexp) * yield / SQUARE(tau_w) *
                                       (1. + nexp * yield / SQUARE(tau_w)));
        q_mag = pre_delP * pgrad;
        dq_gradp =
            -CUBE(H) / (k_turb * mu) * pow((tau_w - yield) / mu, 1. / nexp - 1.) *
            (1. + (nexp - 1.) / (nexp + 1.) * yield / tau_w *
                      (1. + 2. * nexp * yield / tau_w * (1. + nexp / (1. - nexp) * yield / tau_w)));
        dq_dH = 2 * q_mag * H_inv + dq_gradp * H_inv * pgrad;
        srate = pow((fabs(tau_w) - yield) / mu, 1. / nexp);
      } else {
        q_mag = 0.;
        dq_gradp = pre_delP = -CUBE(H) / (k_turb * mu) - beta_slip * SQUARE(H);
      }
    } else { /*  Newtonian type models */
      k_turb = 3.;
      dq_gradp = pre_delP = -CUBE(H) / (k_turb * mu) - beta_slip * SQUARE(H);
      q_mag = pre_delP * pgrad;
      dq_dH = (-3. * SQUARE(H) / (k_turb * mu) - 2. * H * beta_slip) * pgrad;
      srate = fabs(tau_w / mu);
      vis_w = mu;
      /*  Only temperature-dependent Newtonian type model is THERMAL  */
      if (gn->ConstitutiveEquation == THERMAL || gn->ConstitutiveEquation == TABLE) {
        dq_dT = q_mag * (-dmu_dT / mu);
      }
    } /*  End of Viscosity Models **/

    dq_gradpt = dq_gradp - beta_slip * SQUARE(H);
    dv_gradp = dq_gradp * H_inv - beta_slip * H;
    dv_dH = dq_dH * H_inv - q_mag / SQUARE(H);
    vpre_delP = pre_delP * H_inv;

    for (i = 0; i < dim; i++) {
      q[i] += q_mag * ev[i];
    }
    if (DOUBLE_NONZERO(vsqr) &&
        (gn->ConstitutiveEquation == BINGHAM || gn->ConstitutiveEquation == BINGHAM_WLF ||
         gn->ConstitutiveEquation == CARREAU || gn->ConstitutiveEquation == CARREAU_WLF ||
         gn->ConstitutiveEquation == POWER_LAW || gn->ConstitutiveEquation == HERSCHEL_BULKLEY)) {
      GOMA_EH(GOMA_ERROR, "Shear-thining moving wall model not finished yet.\n");
    } else {
      for (i = 0; i < dim; i++) {
        q[i] += H * veloL[i];
      }
    }

    memset(v_avg, 0.0, sizeof(double) * DIM);
    for (i = 0; i < dim; i++) {
      v_avg[i] += q[i] * H_inv;
    }

    /******* CALCULATE FLOW RATE SENSITIVITIES ***********/

    /*Evaluate flowrate sensitivity w.r.t. pressure */
    dbl D_Q_DP[DIM][MDE], D_V_DP[DIM][MDE];
    dbl D_Q_DGRADP[DIM][DIM];
    dbl D_V_DGRADP[DIM][DIM];
    dbl phi_j, grad_phi_j[DIM], grad_II_phi_j[DIM], d_grad_II_phi_j_dmesh[DIM][DIM][MDE];
    memset(D_Q_DP, 0.0, sizeof(double) * DIM * MDE);
    memset(D_V_DP, 0.0, sizeof(double) * DIM * MDE);
    memset(D_Q_DGRADP, 0.0, sizeof(double) * DIM * DIM);
    memset(D_V_DGRADP, 0.0, sizeof(double) * DIM * DIM);
    for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
        D_Q_DGRADP[i][j] += dq_gradpt * ev[i] * ev[j] + pre_delP * dev_dpg[i][j];
        D_V_DGRADP[i][j] += dv_gradp * ev[i] * ev[j] + vpre_delP * dev_dpg[i][j];
      }
    }
    for (i = 0; i < dim; i++) {
      for (j = 0; j < ei[pg->imtrx]->dof[SHELL_FILMP]; j++) {
        ShellBF(SHELL_FILMH, j, &phi_j, grad_phi_j, grad_II_phi_j, d_grad_II_phi_j_dmesh,
                n_dof[MESH_DISPLACEMENT1], dof_map);
        for (k = 0; k < dim; k++) {
          D_Q_DP[i][j] += D_Q_DGRADP[i][k] * grad_II_phi_j[k];
          D_V_DP[i][j] += D_V_DGRADP[i][k] * grad_II_phi_j[k];
        }
      }
    }
    /*Evaluate flowrate sensitivity w.r.t. height */
    dbl D_Q_DH1[DIM][MDE], D_V_DH1[DIM][MDE];
    dbl D_Q_DH2[DIM][MDE], D_V_DH2[DIM][MDE];
    dbl D_Q_DH[DIM][MDE], D_V_DH[DIM][MDE];
    memset(D_Q_DH1, 0.0, sizeof(double) * DIM * MDE);
    memset(D_Q_DH2, 0.0, sizeof(double) * DIM * MDE);
    memset(D_Q_DH, 0.0, sizeof(double) * DIM * MDE);
    memset(D_V_DH1, 0.0, sizeof(double) * DIM * MDE);
    memset(D_V_DH2, 0.0, sizeof(double) * DIM * MDE);
    memset(D_V_DH, 0.0, sizeof(double) * DIM * MDE);
    for (i = 0; i < dim; i++) {
      for (k = 0; k < ei[pg->imtrx]->dof[SHELL_FILMH]; k++) {
        for (j = 0; j < dim; j++) {
          D_Q_DH1[i][k] += D_Q_DGRADP[i][j] * (-D_GRAD_DISJ_PRESS_DH1[j][k]);
          D_Q_DH2[i][k] += D_Q_DGRADP[i][j] * (-D_GRAD_DISJ_PRESS_DH2[i][k]);
          D_Q_DH[i][k] += D_Q_DGRADP[i][j] * (-D_GRAD_DISJ_PRESS_DH[i][k]);
          D_V_DH1[i][k] += D_V_DGRADP[i][j] * (-D_GRAD_DISJ_PRESS_DH1[j][k]);
          D_V_DH2[i][k] += D_V_DGRADP[i][j] * (-D_GRAD_DISJ_PRESS_DH2[i][k]);
          D_V_DH[i][k] += D_V_DGRADP[i][j] * (-D_GRAD_DISJ_PRESS_DH[i][k]);
        }
        D_Q_DH2[i][k] += dq_dH * ev[i];
        D_Q_DH2[i][k] += veloL[i];
        D_V_DH2[i][j] += dv_dH * ev[i];

        ShellBF(SHELL_FILMH, k, &phi_j, grad_phi_j, grad_II_phi_j, d_grad_II_phi_j_dmesh,
                n_dof[MESH_DISPLACEMENT1], dof_map);

        D_Q_DH[i][k] += dq_dH * phi_j * ev[i];
        D_V_DH[i][j] += dv_dH * phi_j * ev[i];
        D_Q_DH[i][k] += veloL[i] * phi_j;
      }
    }

    /*Evaluate flowrate sensitivity w.r.t. particles volume fraction, if applicable */
    dbl D_Q_DC[DIM][MDE], D_V_DC[DIM][MDE];
    memset(D_V_DC, 0.0, sizeof(double) * DIM * MDE);
    memset(D_Q_DC, 0.0, sizeof(double) * DIM * MDE);
    if (pd->v[pg->imtrx][SHELL_PARTC]) {
      for (i = 0; i < dim; i++) {
        for (j = 0; j < ei[pg->imtrx]->dof[SHELL_PARTC]; j++) {
          D_Q_DC[i][j] += q_mag * (-dmu_dc / mu) * ev[i];
          D_V_DC[i][j] += vpre_delP * pgrad * (-dmu_dc / mu) * ev[i];
        }
      }
    }

    /* Sensitivity w.r.t. sh_t , i.e., scalar viscosity dependencies */
    dbl D_Q_DT[DIM], D_V_DT[DIM];
    if (pd->v[pg->imtrx][SHELL_TEMPERATURE]) {
      for (i = 0; i < dim; i++) {
        D_Q_DT[i] = dq_dT * ev[i];
        D_V_DT[i] = H_inv * dq_dT * ev[i];
      }
    }
    /*Evaluate flowrate sensitivity w.r.t. mesh and/or real-solid */
    dbl D_Q_DX[DIM][DIM][MDE];
    dbl D_Q_DRS[DIM][DIM][MDE];
    dbl D_V_DX[DIM][DIM][MDE];
    memset(D_V_DX, 0.0, sizeof(double) * DIM * DIM * MDE);
    memset(D_Q_DX, 0.0, sizeof(double) * DIM * DIM * MDE);
    memset(D_Q_DRS, 0.0, sizeof(double) * DIM * DIM * MDE);

    switch (mp->FSIModel) {
    case FSI_MESH_CONTINUUM:
    case FSI_MESH_UNDEF:
    case FSI_SHELL_ONLY_MESH:
    case FSI_SHELL_ONLY_UNDEF:
      for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
          for (k = 0; k < n_dof[MESH_DISPLACEMENT1]; k++) {
            D_Q_DX[i][j][k] += dq_dH * ev[i] * D_H_DX[j][k];
            D_V_DX[i][j][k] += dv_dH * ev[i] * D_H_DX[j][k];
            for (jk = 0; jk < dim; jk++) {
              D_Q_DX[i][j][k] += D_Q_DGRADP[i][jk] * D_GRADP_DX[jk][j][k];
              D_V_DX[i][j][k] += D_V_DGRADP[i][jk] * D_GRADP_DX[jk][j][k];
            }
            // Ignore dependency of GRAD_DISJ_PRESS w.r.t. mesh and height for now

            D_Q_DX[i][j][k] += veloL[i] * D_H_DX[j][k];
          }
        }
      }
      break;
    }

    /******* STORE THE INFORMATION TO LUBRICATION AUXILIARIES STRUCTURE ***********/

    LubAux->H = H;
    LubAux->gradP_mag = 0;
    LubAux->srate = srate;
    LubAux->mu_star = vis_w;
    for (i = 0; i < dim; i++) {
      LubAux->q[i] = q[i];
      LubAux->v_avg[i] = v_avg[i];
      LubAux->gradP_mag += SQUARE(pg_cmp[i]);
      LubAux->dq_dT[i] = D_Q_DT[i];
      LubAux->dv_avg_dT[i] = D_V_DT[i];

      for (j = 0; j < ei[pg->imtrx]->dof[SHELL_FILMH]; j++) {
        LubAux->dq_dh1[i][j] = D_Q_DH1[i][j];
        LubAux->dq_dh2[i][j] = D_Q_DH2[i][j];
        LubAux->dq_dh[i][j] = D_Q_DH[i][j];

        LubAux->dv_avg_dh1[i][j] = D_V_DH1[i][j];
        LubAux->dv_avg_dh2[i][j] = D_V_DH2[i][j];
        LubAux->dv_avg_dh[i][j] = D_V_DH[i][j];
      }

      for (j = 0; j < ei[pg->imtrx]->dof[SHELL_FILMP]; j++) {
        LubAux->dq_dp[i][j] = D_Q_DP[i][j];
      }
      for (j = 0; j < dim; j++) {
        LubAux->dq_dgradp[i][j] = D_Q_DGRADP[i][j];
      }

      if (pd->v[pg->imtrx][SHELL_PARTC]) {
        for (j = 0; j < ei[pg->imtrx]->dof[SHELL_PARTC]; j++) {
          LubAux->dq_dc[i][j] = D_Q_DC[i][j];
          LubAux->dv_avg_dc[i][j] = D_V_DC[i][j];
        }
      }
      for (j = 0; j < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; j++) {
        jk = dof_map[j];
        for (k = 0; k < dim; k++) {
          LubAux->dq_dx[i][k][j] = D_Q_DX[i][k][jk];
          LubAux->dv_avg_dx[i][k][j] = D_V_DX[i][k][jk];
        }
      }
    }
    LubAux->gradP_mag = sqrt(LubAux->gradP_mag);

    // Cleanup
    fv->wt = wt_old;
    safe_free((void *)n_dof);
  }

  return;

} /* End of calculate_lub_q_v */

void calculate_lub_q_v_old(
    const int EQN, double time_old, double dt_old, double xi[DIM], const Exo_DB *exo)
/******************************************************************************
 *
 * calculate_lub_q_v_old()
 *
 * Function to calculate flow rate per unit width (q) and average velocity (v)
 * in lubrication flow at old time step - for Taylor Galerkin stabilization
 * purpose.
 *
 *
 * Kris Tjiptowidjojo tjiptowi@unm.edu
 *
 * EDITED:
 * 2010-11-30: Scott Roberts - sarober@sandia.gov
 *    Re-wrote the lubrication section to include more key physics.
 *    Added calculation of full Jacobian entries.
 *
 ******************************************************************************/
{
  int i;
  dbl ev[DIM], pgrad, pg_cmp[DIM];
  dbl q_old[DIM];
  dbl v_avg_old[DIM];
  dbl H_old;
  dbl veloL_old[DIM], veloU_old[DIM];
  dbl mu_old;
  dbl rho_old;
  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct; /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;
  int VAR;

  /* Calculate flow rate and average velocity with their sensitivities
   * depending on the lubrication model employed
   */

  /* Confined lubrication flow - Newtonian */
  if ((EQN == R_LUBP) || (EQN == R_LUBP_2)) {

    /* Set proper fill variable first.   If in lub_p layer, then use FILL,
     * but if in LUBP_2 layer use PHASE1.  This will have to be made more general
     * if you wanted to do both LS phases and R_phaseN fields in same layer.
     * We will leave that to the next sucker to develop
     */

    VAR = FILL;
    if (EQN == R_LUBP_2) {
      VAR = PHASE1;
    }

    /***** INITIALIZE LUBRICATION COMPONENTS AND LOAD IN VARIABLES *****/

    /* Problem dimensions */
    int dim = pd->Num_Dim;

    /* Setup lubrication shell constructs */
    dbl wt_old = fv->wt;
    int dof_map[MDE];
    int *n_dof = (int *)array_alloc(1, MAX_VARIABLE_TYPES, sizeof(int));
    lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

    /* Load viscosity and density */
    rho_old = density(d_rho, time_old);
    if (gn->ConstitutiveEquation == POWER_LAW || gn->ConstitutiveEquation == BINGHAM ||
        gn->ConstitutiveEquation == HERSCHEL_BULKLEY || gn->ConstitutiveEquation == CARREAU ||
        gn->ConstitutiveEquation == CARREAU_WLF || gn->ConstitutiveEquation == BINGHAM_WLF) {
      mu_old = gn->mu0;
    } else {
      mu_old = viscosity(gn, NULL, d_mu);
    }

    /* Extract wall velocities */
    velocity_function_model(veloU_old, veloL_old, time_old, dt_old);

    /* Extract wall heights */
    dbl H_U, dH_U_dtime, H_L, dH_L_dtime;
    dbl dH_U_dX[DIM], dH_L_dX[DIM], dH_U_dp, dH_U_ddh, dH_dF[MDE];
    H_old = height_function_model(&H_U, &dH_U_dtime, &H_L, &dH_L_dtime, dH_U_dX, dH_L_dX, &dH_U_dp,
                                  &dH_U_ddh, dH_dF, time_old, dt_old);

    /***** DEFORM HEIGHT AND CALCULATE K_TURB VALUE*****/

    /* Deform height */
    switch (mp->FSIModel) {
    case FSI_MESH_CONTINUUM:
    case FSI_MESH_UNDEF:
    case FSI_SHELL_ONLY_UNDEF:
      for (i = 0; i < dim; i++) {
        H_old -= fv->snormal[i] * fv_old->d[i];
      }
      break;
    case FSI_SHELL_ONLY_MESH:
      if (pd->e[pg->imtrx][R_SHELL_NORMAL1] && pd->e[pg->imtrx][R_SHELL_NORMAL2] &&
          pd->e[pg->imtrx][R_SHELL_NORMAL3]) {
        for (i = 0; i < dim; i++) {
          H_old -= fv->n[i] * fv_old->d[i];
        }
      } else {
        for (i = 0; i < dim; i++) {
          H_old -= fv->snormal[i] * fv_old->d[i];
        }
      }
      break;
    case FSI_REALSOLID_CONTINUUM:
      for (i = 0; i < dim; i++) {
        H_old -= fv->snormal[i] * fv_old->d_rs[i];
      }
      break;
    }

    /***** CALCULATE PRESSURE GRADIENT *****/

    /* Define variables */
    dbl GRADP[DIM];
    memset(GRADP, 0.0, sizeof(double) * DIM);

    Inn(fv_old->grad_lubp, GRADP);

    /***** CALCULATE GRAVITY AND LORENTZ (OTHER lubmomsource) / BODY FORCE AND SENSITIVITIES
     * *****/

    dbl GRAV[DIM];
    memset(GRAV, 0.0, sizeof(double) * DIM);

    /* Calculate and rotate body force, calculate mesh derivatives */
    for (i = 0; i < DIM; i++) {
      GRAV[i] = rho_old * mp->momentum_source[i];
    }

    /***** CALCULATE HEAVISIDE GRADIENT*****/

    dbl GRADH[DIM];
    memset(GRADH, 0.0, sizeof(double) * DIM);

    if (pd->v[pg->imtrx][VAR]) {
      for (i = 0; i < DIM; i++) {
        GRADH[i] = lsi->gradHn_old[i];
      }
    }

    /***** CALCULATE CURVATURE  *****/

    dbl CURV = 0.0;
    double H_cap, dHc_U_dX[DIM], dHc_L_dX[DIM];

    /* Curvature - analytic in the "z" direction  */
    dbl dcaU, dcaL, slopeU, slopeL;
    dcaU = dcaL = slopeU = slopeL = 0.0;
    if (mp->HeightUFunctionModel == WALL_DISTMOD || mp->HeightUFunctionModel == WALL_DISTURB) {
      H_cap = MAX(H_U - H_L, DBL_SEMI_SMALL);
      memset(dHc_U_dX, 0.0, sizeof(double) * DIM);
      memset(dHc_L_dX, 0.0, sizeof(double) * DIM);
    } else {
      H_cap = H_old;
      for (i = 0; i < dim; i++) {
        dHc_U_dX[i] = dH_U_dX[i];
        dHc_L_dX[i] = dH_L_dX[i];
      }
    }
    if (pd->v[pg->imtrx][VAR]) {
      dcaU = mp->dcaU * M_PIE / 180.0;
      dcaL = mp->dcaL * M_PIE / 180.0;
      slopeU = slopeL = 0.;
      for (i = 0; i < dim; i++) {
        slopeU += dHc_U_dX[i] * lsi->normal[i];
        slopeL += dHc_L_dX[i] * lsi->normal[i];
      }
      CURV = -(cos(dcaU + atan(slopeU)) + cos(dcaL + atan(-slopeL))) / H_cap;
    }

    /* Curvature - numerical in planview direction */
    if (pd->e[pg->imtrx][SHELL_LUB_CURV]) {
      CURV += fv_old->sh_l_curv;
    }
    if (pd->e[pg->imtrx][SHELL_LUB_CURV_2]) {
      CURV += fv_old->sh_l_curv_2;
    }
    if (pd->e[pg->imtrx][NORMAL1]) {
      CURV += fv_old->div_n;
    }

    /******* CALCULATE FLOW RATE AND AVERAGE VELOCITY ***********/

    /* Set some coefficients */
    dbl k_turb = 12.;
    dbl vsqr, q_mag = 0., tau_w, q_mag2;
    dbl dqmag_dF[MDE], factor, ratio;
    int nonmoving_model =
        (gn->ConstitutiveEquation == BINGHAM || gn->ConstitutiveEquation == BINGHAM_WLF ||
         gn->ConstitutiveEquation == CARREAU || gn->ConstitutiveEquation == CARREAU_WLF);
    int movwall_model =
        (gn->ConstitutiveEquation == BINGHAM || gn->ConstitutiveEquation == BINGHAM_WLF ||
         gn->ConstitutiveEquation == CARREAU || gn->ConstitutiveEquation == CARREAU_WLF ||
         gn->ConstitutiveEquation == POWER_LAW || gn->ConstitutiveEquation == HERSCHEL_BULKLEY);
    int movingwall = FALSE;

    /* Calculate flow rate and velocity */
    memset(ev, 0.0, sizeof(double) * DIM);
    for (i = 0; i < dim; i++) {
      pg_cmp[i] = GRADP[i] - GRAV[i];
      if (pd->v[pg->imtrx][VAR])
        pg_cmp[i] += GRADH[i] * CURV * mp->surface_tension;
    }
    pgrad = 0.;
    vsqr = 0.;
    for (i = 0; i < dim; i++) {
      pgrad += SQUARE(pg_cmp[i]);
      vsqr += SQUARE(veloL_old[i]);
      vsqr += SQUARE(veloU_old[i]);
    }
    pgrad = sqrt(pgrad);
    if (DOUBLE_NONZERO(pgrad)) {
      for (i = 0; i < dim; i++) {
        ev[i] = pg_cmp[i] / pgrad;
      }
    } else {
      ev[0] = 1.;
    }

    if (!movingwall) {
      /*  First non-Newtonian models with analytical viscosity integration */
      tau_w = 0.5 * H_old * pgrad;
      if (gn->ConstitutiveEquation == POWER_LAW) {
        double nexp = gn->nexp;
        k_turb = 4. * (2. + 1. / nexp);
        if (DOUBLE_NONZERO(pgrad)) {
          q_mag = -2. * SQUARE(H_old) / k_turb * pow(tau_w / mu_old, 1. / nexp);
        } else {
          q_mag = 0.;
        }
      } else if (gn->ConstitutiveEquation == HERSCHEL_BULKLEY) {
        double nexp = gn->nexp, yield = gn->tau_y;
        k_turb = 4. * (2 * nexp + 1.);
        if (tau_w > yield) {
          q_mag = -CUBE(H_old) / k_turb * pow((tau_w - yield) / mu_old, 1. / nexp) *
                  (nexp / tau_w - nexp / (1. + nexp) * yield / SQUARE(tau_w) *
                                      (1. + nexp * yield / SQUARE(tau_w)));
          q_mag *= pgrad;
        } else {
          q_mag = 0.;
        }
      } else if (nonmoving_model) {
        int err;
        if (isnan(tau_w))
          GOMA_WH(GOMA_ERROR, "Trouble, tau_w is nan...\n");
        err = lub_viscosity_integrate(tau_w, H_old, &q_mag, NULL, NULL, NULL, NULL, NULL, NULL);
        if (err < 0) {
          GOMA_WH(GOMA_ERROR, "Some trouble with Numerical Lubrication...\n");
        }
        /* moving wall parts  */
      } else { /*  Newtonian type models */
        k_turb = 12.;
        q_mag = -CUBE(H_old) / (k_turb * mu_old);
        q_mag *= pgrad;
      } /*  End of Viscosity Models **/
    } else {
      if (movwall_model) { /*  non-Newtonian models with a moving wall */
        double relax = 0.5;
        int guess = 0, err;
        err = lub2D_flow2D(gn, pg_cmp, ev, NULL, q_old, NULL, NULL, H_old, NULL, veloL_old,
                           veloU_old, guess, NULL, relax);
        if (err < 0) {
          GOMA_WH(GOMA_ERROR, "Some trouble with Numerical Lubrication...\n");
        }
      } else { /*  moving wall part of Newtonian type models */
        k_turb = 12.;
        q_mag = -CUBE(H_old) / (k_turb * mu_old) * pgrad;
      } /*  End of Viscosity Models **/
    }

    if (pd->v[pg->imtrx][VAR] && mp->mp2nd->ViscosityModel == RATIO) {
      ratio = 1. / mp->mp2nd->viscosity; /* Assuming model = RATIO for now */
      q_mag2 = q_mag * ratio;
      q_mag =
          ls_modulate_property(q_mag, q_mag2, ls->Length_Scale, (double)mp->mp2nd->viscositymask[0],
                               (double)mp->mp2nd->viscositymask[1], dqmag_dF, &factor);
      factor *= (1. - ratio);
      factor += ratio;
    }

    memset(q_old, 0.0, sizeof(double) * DIM);
    for (i = 0; i < dim; i++) {
      q_old[i] += q_mag * ev[i];
    }
    if (DOUBLE_NONZERO(vsqr) && movwall_model) {
      GOMA_EH(GOMA_ERROR, "Shear-thining moving wall model not finished yet.\n");
    } else {
      for (i = 0; i < dim; i++) {
        q_old[i] += 0.5 * H_old * (veloL_old[i] + veloU_old[i]);
      }
    }
    memset(v_avg_old, 0.0, sizeof(double) * DIM);
    for (i = 0; i < dim; i++) {
      v_avg_old[i] += q_old[i] / H_old;
    }

    /******* STORE THE INFORMATION TO LUBRICATION AUXILIARIES STRUCTURE ***********/

    for (i = 0; i < DIM; i++) {
      LubAux_old->q[i] = q_old[i];
      LubAux_old->v_avg[i] = v_avg_old[i];
    }

    // Cleanup
    fv->wt = wt_old;
    safe_free((void *)n_dof);

  }

  /* Film flow - Newtonian */
  else if (pd->e[pg->imtrx][R_SHELL_FILMP]) {

    /******* PRECALCULATE ALL NECESSARY COMPONENTS ***********/
    /* Problem dimensions */
    int dim = pd->Num_Dim;

    /* Setup lubrication shell constructs */
    dbl wt_old = fv->wt;
    int dof_map[MDE];
    int *n_dof = (int *)array_alloc(1, MAX_VARIABLE_TYPES, sizeof(int));
    lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

    /* Load viscosity */
    viscosity(gn, NULL, d_mu);
    mu_old = mp_old->viscosity;

    /* Get slip coefficient */
    double beta_slip;
    beta_slip = mp->SlipCoeff;

    /* Extract bottom wall velocity */
    velocity_function_model(veloU_old, veloL_old, time_old, dt_old);

    /***** CALCULATE HEIGHT, SLOPES, AND SENSITIVITIES *****/

    dbl GRADH[DIM];
    memset(GRADH, 0.0, sizeof(double) * DIM);

    /* Extract film thickness */
    H_old = fv_old->sh_fh;

    /* Perfrom I - nn gradient */
    Inn(fv_old->grad_sh_fh, GRADH);

    /***** CALCULATE GRAVITY AND LORENTZ (OTHER lubmomsource) / BODY FORCE AND SENSITIVITIES
     * *****/

    dbl GRAV[DIM];
    memset(GRAV, 0.0, sizeof(double) * DIM);

    /* Calculate and rotate body force, calculate mesh derivatives */
    for (i = 0; i < DIM; i++) {
      GRAV[i] = mp->momentum_source[i];
    }

    /***** CALCULATE PRESSURE GRADIENT AND SENSITIVITIES *****/

    /* Define variables */
    dbl GRADP[DIM];
    memset(GRADP, 0.0, sizeof(double) * DIM);

    /* Perfrom I - nn gradient */
    Inn(fv_old->grad_sh_fp, GRADP);

    /***** CALCULATE DISJOINING PRESSURE GRADIENT AND SENSITIVITIES *****/

    /* Define variables */

    dbl GRAD_DISJ_PRESS[DIM];
    dbl D_GRAD_DISJ_PRESS_DH1[DIM][MDE], D_GRAD_DISJ_PRESS_DH2[DIM][MDE];
    dbl D_GRAD_DISJ_PRESS_DH[DIM][MDE];
    memset(GRAD_DISJ_PRESS, 0.0, sizeof(double) * DIM);
    memset(D_GRAD_DISJ_PRESS_DH1, 0.0, sizeof(double) * DIM * MDE);
    memset(D_GRAD_DISJ_PRESS_DH2, 0.0, sizeof(double) * DIM * MDE);
    memset(D_GRAD_DISJ_PRESS_DH, 0.0, sizeof(double) * DIM * MDE);

    /* Evaluate disjoining pressure and its sensitivities */
    disjoining_pressure_model(fv_old->sh_fh, fv_old->grad_sh_fh, n_dof, dof_map, GRAD_DISJ_PRESS,
                              D_GRAD_DISJ_PRESS_DH1, D_GRAD_DISJ_PRESS_DH2, D_GRAD_DISJ_PRESS_DH);

    /******* CALCULATE FLOW RATE AND AVERAGE VELOCITY ***********/

    double q_old[DIM], v_avg_old[DIM], pg_cmp[DIM], pgrad, ev[DIM];
    memset(q_old, 0.0, sizeof(double) * DIM);
    memset(v_avg_old, 0.0, sizeof(double) * DIM);
    dbl q_mag = 0., k_turb = 3., tau_w, vis_w, srate, vsqr;

    for (i = 0; i < dim; i++) {
      pg_cmp[i] = GRADP[i] - GRAV[i] - GRAD_DISJ_PRESS[i];
    }
    pgrad = 0.;
    vsqr = 0.;
    for (i = 0; i < dim; i++) {
      pgrad += SQUARE(pg_cmp[i]);
      vsqr += SQUARE(veloL_old[i]);
    }
    pgrad = sqrt(pgrad);
    if (DOUBLE_NONZERO(pgrad)) {
      for (i = 0; i < dim; i++) {
        ev[i] = pg_cmp[i] / pgrad;
      }
    } else {
      ev[0] = 1.;
    }
    /* Evaluate flow rate and average velocity */
    tau_w = H_old * pgrad;
    if (gn->ConstitutiveEquation == POWER_LAW) {
      double nexp = gn->nexp;
      k_turb = 2. + 1. / nexp;
      if (DOUBLE_NONZERO(pgrad)) {
        q_mag = -SQUARE(H_old) / k_turb * pow(tau_w / mu_old, 1. / nexp);
        srate = pow(fabs(tau_w) / mu_old, 1. / nexp);
        vis_w = tau_w / srate;
      } else {
        q_mag = 0.;
        srate = 0.;
        vis_w = mu_old;
      }
    } else if (gn->ConstitutiveEquation == BINGHAM || gn->ConstitutiveEquation == BINGHAM_WLF ||
               gn->ConstitutiveEquation == CARREAU || gn->ConstitutiveEquation == CARREAU_WLF) {
      int err;
      if (isnan(tau_w))
        GOMA_WH(GOMA_ERROR, "Trouble, tau_w is nan...\n");
      err = lub_viscosity_integrate(tau_w, H_old, &q_mag, NULL, NULL, &srate, NULL, &vis_w, NULL);
      if (isnan(srate))
        fprintf(stderr, "srate isnan %g %g %g %g\n", tau_w, q_mag, srate, vis_w);
      if (err < 0) {
        GOMA_WH(GOMA_ERROR, "Some trouble with Numerical Lubrication...\n");
      }

      /** Make corrections for film flow from confined calculations **/
      q_mag *= 2.;
      q_mag -= beta_slip * SQUARE(H_old) * pgrad;
    } else if (gn->ConstitutiveEquation == HERSCHEL_BULKLEY) {
      double nexp = gn->nexp, yield = gn->tau_y;
      k_turb = (2 * nexp + 1.);
      if (tau_w > yield) {
        q_mag = -pgrad * CUBE(H_old) / k_turb * pow((tau_w - yield) / mu_old, 1. / nexp) *
                (nexp / tau_w -
                 nexp / (1. + nexp) * yield / SQUARE(tau_w) * (1. + nexp * yield / SQUARE(tau_w)));
        srate = pow((fabs(tau_w) - yield) / mu_old, 1. / nexp);
      } else {
        q_mag = 0.;
      }
    } else { /*  Newtonian type models */
      k_turb = 3.;
      q_mag = -pgrad * (CUBE(H_old) / (k_turb * mu_old) + beta_slip * SQUARE(H_old));
      srate = fabs(tau_w / mu_old);
      vis_w = mu_old;
    } /*  End of Viscosity Models **/

    for (i = 0; i < dim; i++) {
      q_old[i] += q_mag * ev[i];
    }
    if (DOUBLE_NONZERO(vsqr) &&
        (gn->ConstitutiveEquation == BINGHAM || gn->ConstitutiveEquation == BINGHAM_WLF ||
         gn->ConstitutiveEquation == CARREAU || gn->ConstitutiveEquation == CARREAU_WLF ||
         gn->ConstitutiveEquation == POWER_LAW || gn->ConstitutiveEquation == HERSCHEL_BULKLEY)) {
      GOMA_EH(GOMA_ERROR, "Shear-thining moving wall model not finished yet.\n");
    } else {
      for (i = 0; i < dim; i++) {
        q_old[i] += H_old * veloL_old[i];
      }
    }

    memset(v_avg_old, 0.0, sizeof(double) * DIM);
    for (i = 0; i < dim; i++) {
      v_avg_old[i] += q_old[i] / H_old;
    }

    /******* STORE THE INFORMATION TO LUBRICATION AUXILIARIES STRUCTURE ***********/

    for (i = 0; i < DIM; i++) {
      LubAux_old->q[i] = q_old[i];
      LubAux_old->v_avg[i] = v_avg_old[i];
    }

    // Cleanup
    fv->wt = wt_old;
    safe_free((void *)n_dof);
  }

  return;

} /* End of calculate_lub_q_v_old */

double shell_saturation_pressure_curve(double Pliq, double *dSdP, double *dSdP_P)
/******************************************************************************
 *
 * shell_saturation_pressure_curve()
 *
 * Calculates the pore saturation for the shell_sat_open porous shell model
 * from the pore pressure and other necessary parameters.  We will try
 * to load most of these parameters internally.  Assumes a normal distribution
 * centered around Rmean with (Rmean-Rmin) = 2*StdDev.
 *
 * Scott A Roberts (1514) sarober@sandia.gov
 * January 14, 2011
 *
 * REVISED:
 * 2011-01-20:  Added a proper derivation using a normal distribution
 *
 ******************************************************************************/
{

  /* Hard-code physical parameters. */
  /* To scale these for with a 0-1 external field with knuckles and pillows:
   * sigma = fv->external_field[efv->ev_SAT_index] * (mp->u_saturation[0] - mp->u_saturation[4]) -
   * mp->u_saturation[4])
   */
  /* This needs to be pushed out to the input deck. */
  dbl sigma = mp->u_saturation[0]; // Surface tension
  dbl theta = mp->u_saturation[1]; // Contact angle
  dbl Patm = mp->PorousShellPatm;  // Atmospheric pressure
  //  dbl R = porous_shell_closed_radius_model();                                    // Pillar
  //  radius
  dbl Rmin = mp->u_saturation[2]; // Minimum pore radius
  dbl Rmax = mp->u_saturation[3]; // Maximum pore radius
  //  dbl phi = mp->porosity;                                  // Porosity (open phase)

  /* Define variables */
  dbl S;

  /* Calculate various pressures */
  dbl Pcap = Patm - Pliq;
  dbl dPdP = -1.0;

#define _USE_MATH_DEFINES

  /* Calculate capillary radius */
  dbl Rc = 2 * sigma * cos(theta / 180.0 * PI) / (Pcap + DBL_SEMI_SMALL);
  dbl Rc_P = -2 * sigma * cos(theta / 180.0 * PI) / pow(Pcap + DBL_SEMI_SMALL, 2) * dPdP;
  dbl Rc_PP = 4 * sigma * cos(theta / 180.0 * PI) / pow(Pcap + DBL_SEMI_SMALL, 3) * pow(dPdP, 2);

  /* Define limiting values */
  dbl Rmax3 = pow(Rmax, 3);
  dbl Rmin3 = pow(Rmin, 3);
  dbl SatCut = 1 - 0.05;
  // dbl RCut = pow( (Rmax3-Rmin3)*(Rmin3*(SatCut-1)-Rmax3*SatCut)/(Rmin3-Rmax3) , 1.0/3.0);
  dbl RCut = Rmin + (Rmax - Rmin) * SatCut;
  dbl PCut = -2 * sigma * cos(theta / 180.0 * PI) / RCut;

  /* Define modulation function */
  // dbl ModK    = 0.1;
  dbl ModK = 1;
  dbl ModF = (1 + tanh(ModK * (Pliq - PCut))) / 2.0;
  dbl cosh_modk_p = 0;
  if (fabs(ModK * (Pliq - PCut)) < 700) {
    cosh_modk_p = cosh(ModK * (Pliq - PCut));
  } else {
    cosh_modk_p = cosh(SGN(ModK * (Pliq - PCut)) * 700);
  }
  dbl ModF_P = ModK * pow(cosh_modk_p, -2) / 2.0;
  dbl ModF_PP = -pow(ModK, 2) * pow(cosh_modk_p, -2) * tanh(ModK * (Pliq - PCut));

  /* Define basic saturation function */
  // dbl BS    = (pow(Rc,3)-Rmin3)/(Rmax3-Rmin3);
  // dbl BS_P  = 3*pow(Rc,2)*Rc_P/(Rmax3-Rmin3);
  // dbl BS_PP = ( 6*Rc*pow(Rc_P,2) + 3*pow(Rc,2)*Rc_PP ) /(Rmax3-Rmin3);
  dbl BS = (Rc - Rmin) / (Rmax - Rmin);
  dbl BS_P = Rc_P / (Rmax - Rmin);
  dbl BS_PP = Rc_PP / (Rmax3 - Rmin3);

  /* Modulate saturation function */
  S = (1 - ModF) * BS + ModF;
  *dSdP = (1 - ModF) * BS_P - ModF_P * BS + ModF_P;
  *dSdP_P = (1 - ModF) * BS_PP - 2 * ModF_P * BS_P - ModF_PP * BS + ModF_PP;

  /* Calculate saturation and derivatives */

  /* Uniform distribution (from SAND1996-2149) */
  /* if ( Rc < Rmin ) { */
  /*   S       = (pow(Rminlim,3)-pow(Rmin,3))/(pow(Rmax,3)-pow(Rmin,3)); */
  /*   *dSdP   = 3*pow(Rminlim,2)*Rc_P/(pow(Rmax,3)-pow(Rmin,3)); */
  /*   *dSdP_P = ( 6*Rminlim*pow(Rc_P,2) + 3*pow(Rminlim,2)*Rc_PP ) /(pow(Rmax,3)-pow(Rmin,3)); */
  /* } else if ( Rc > Rmax ) { */
  /*   S       = (pow(Rmaxlim,3)-pow(Rmin,3))/(pow(Rmax,3)-pow(Rmin,3)); */
  /*   *dSdP   = 3*pow(Rmaxlim,2)*Rc_P/(pow(Rmax,3)-pow(Rmin,3)); */
  /*   *dSdP_P = ( 6*Rmaxlim*pow(Rc_P,2) + 3*pow(Rmaxlim,2)*Rc_PP ) /(pow(Rmax,3)-pow(Rmin,3)); */
  /* } else { */
  /*   S       = (pow(Rc,3)-pow(Rmin,3))/(pow(Rmax,3)-pow(Rmin,3)); */
  /*   *dSdP   = 3*pow(Rc,2)*Rc_P/(pow(Rmax,3)-pow(Rmin,3)); */
  /*   *dSdP_P = ( 6*Rc*pow(Rc_P,2) + 3*pow(Rc,2)*Rc_PP ) /(pow(Rmax,3)-pow(Rmin,3)); */
  /* } */

  /* Calculate extrema radii */
  /* dbl Rmin = -2*R + sqrt(2*PI)*pow(R,2)/(pow(3.0,0.25)*sqrt(pow(R,2)*(1-phi))); */
  /* dbl Rmax = 2*R*(sqrt(3.0)-1) + sqrt(3.0)*Rmin; */
  /* dbl Rmean = (Rmin + Rmax) / 2.0; */
  /* dbl Rstdv = (Rmean - Rmin) / 2.0; */

  /* Calculate saturation and derivatives */

  /* Weighted exponential distribution (from SAND1996-2149) */
  /* S = 4.0/3.0*pow(Rc/Rmean,3) + 2*pow(Rc/Rmean,2) + 2*Rc/Rmean + 1; */
  /* S = 1 - S*exp(-2*Rc/Rmean); */
  /* *dSdP = 8*exp(-2*Rc/Rmean)*pow(Rc,3)*Rc_P/(3*pow(Rmean,4)); */
  /* *dSdP_P  = (3*Rmean-2*Rc)*pow(Rc_P,2) + Rmean*Rc*Rc_PP; */
  /* *dSdP_P *= 8*exp(-2*Rc/Rmean)*pow(Rc,2)/(3*pow(Rmean,5)); */

  /* Self-derived Gaussian distribution */
  /* S  = erf(Rmean/(sqrt(2.0)*Rstdv)) - erf((Rmean-Rc)/(sqrt(2.0)*Rstdv)); */
  /* S *= exp(pow(Rmean,2)/(2.0*pow(Rstdv,2)))*sqrt(2*PI)*(pow(Rmean,2)+pow(Rstdv,2)); */
  /* S += 2.0*Rmean*Rstdv; */
  /* S *= exp(pow(Rmean-Rc,2)/(2.0*pow(Rstdv,2))); */
  /* S -= 2.0*exp(pow(Rmean,2)/(2.0*pow(Rstdv,2)))*Rstdv*(Rmean+Rc); */
  /* S *= exp(-(2.0*pow(Rmean,2)-2*Rmean*Rc+pow(Rc,2))/(2.0*pow(Rstdv,2))); */
  /* S /= 2.0*sqrt(2.0*PI)*( */
  /* 			 exp(-pow(Rmean,2)/(2.0*pow(Rstdv,2)))*Rmean*Rstdv/sqrt(2.0*PI) */
  /* 			 + (pow(Rmean,2)+pow(Rstdv,2))*(1+erf(Rmean/(sqrt(2.0)*Rstdv)))/2.0 */
  /* 			 ); */
  /* *dSdP  = sqrt(2.0)*exp((2*Rmean-Rc)*Rc/(2.0*pow(Rstdv,2)))*pow(Rc,2)*Rc_P; */
  /* *dSdP /= Rstdv*( */
  /* 		  sqrt(2.0)*Rmean*Rstdv +  */
  /* 		  exp(pow(Rmean,2)/(2.0*pow(Rstdv,2)))*sqrt(PI)*(pow(Rmean,2)+pow(Rstdv,2))*(1+erf(Rmean/(sqrt(2.0)*Rstdv)))
   */
  /* 		  ); */
  /* *dSdP_P  = (2.0*pow(Rstdv,2) + (Rmean-Rc)*Rc)*pow(Rc_P,2) + pow(Rstdv,2)*Rc*Rc_PP; */
  /* *dSdP_P *= sqrt(2.0)*exp((2*Rmean-Rc)*Rc/pow(2.0*Rstdv,2))*Rc; */
  /* *dSdP_P /= pow(Rstdv,3)*( */
  /* 			   sqrt(2.0)*Rmean*Rstdv +  */
  /* 			   exp(pow(Rmean,2)/(2.0*pow(Rstdv,2)))*sqrt(PI)*(pow(Rmean,2)+pow(Rstdv,2))*(1+erf(Rmean/(sqrt(2.0)*Rstdv)))
   */
  /* 			   ); */

  /* Finalize */
  return (S);

} /* End of shell_saturation_pressure_curve() */

void ShellBF(int ev, // Equation or variable to fetch basis functions
             int ii, // Integer for which DOF to fetch
             double *phi_i,
             double grad_phi_i[DIM],
             double gradII_phi_i[DIM],
             double d_gradII_phi_i_dx[DIM][DIM][MDE],
             int ndof,
             int dof_map[MDE])
/******************************************************************************
 *
 * ShellBF()
 *
 * Calculate all basis functions and mesh sensitivities for a given equation
 * and DOF.  Returns both the regular gradients and II gradients.
 *
 * Scott A Roberts (1514) sarober@sandia.gov
 *
 ******************************************************************************/
{
  /* Variable definitions */
  int i, j, k, jk;
  double d_grad_phi_i_dx[DIM][DIM][MDE];

  /* Basis function */
  *phi_i = bf[ev]->phi[ii];

  // check if bar2 problem on s_integration
  if (pd->Num_Dim == 2 && mp->ehl_integration_kind == SIK_S) {
    ShellBF_2d_bar(ev, ii, gradII_phi_i, d_gradII_phi_i_dx);
    return;
  }

  /* Grad of basis function */
  for (i = 0; i < pd->Num_Dim; i++) {
    grad_phi_i[i] = bf[ev]->grad_phi[ii][i];
  }

  /* Mesh derivatives of basis function grad */
  memset(d_grad_phi_i_dx, 0.0, sizeof(double) * DIM * DIM * MDE);
  for (i = 0; i < pd->Num_Dim; i++) {
    for (j = 0; j < pd->Num_Dim; j++) {
      for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
        jk = dof_map[k];
        d_grad_phi_i_dx[i][j][jk] = bf[ev]->d_grad_phi_dmesh[ii][i][j][k];
      }
    }
  }

  /* Rotate basis function and calculate mesh */
  ShellRotate(grad_phi_i, d_grad_phi_i_dx, gradII_phi_i, d_gradII_phi_i_dx, ndof);

  return;
} /*** END OF ShellBF ***/

void calculate_lub_q_v_nonnewtonian(double time, double dt)
/******************************************************************************
 *
 * calculate_lub_q_v_nonnewtonian()
 *
 * Function to calculate flow rate per unit width (q) and average velocity (v)
 * in lubrication flow of non Newtonian liquid
 *
 *
 * Kris Tjiptowidjojo tjiptowi@unm.edu
 *
 ******************************************************************************/
{
  int i, j;
  dbl q[DIM];
  dbl v_avg[DIM];
  dbl grad_P[DIM], grad_II_P[DIM];
  dbl veloL[DIM], veloU[DIM];
  dbl mu0, nexp;
  dbl shear_top, shear_bot, cross_shear, gradP_mag;
  dbl H, H_U, dH_U_dtime, H_L, dH_L_dtime, dH_U_ddH;
  dbl dH_U_dX[DIM], dH_L_dX[DIM], dH_U_dp, dH_dF[MDE];
  dbl gradP_tangent[DIM], gradP_normal[DIM], gradP_normal_init[DIM];

  dbl epsilon = 1.0e-5;

  /******* PRECALCULATE ALL NECESSARY COMPONENTS ***********/

  /* Extract shear rates */
  shear_top = fv->sh_shear_top;
  shear_bot = fv->sh_shear_bot;
  cross_shear = fv->sh_cross_shear;

  /* Extract pressure gradient magnitude */
  for (i = 0; i < DIM; i++) {
    grad_P[i] = fv->grad_lubp[i];
  }

  /* Perfrom I - nn gradient */
  Inn(grad_P, grad_II_P);

  /* Extract top and bottom walls velocities */
  velocity_function_model(veloU, veloL, time, dt);

  /*Find tangent vector of pressure gradient */

  gradP_mag = 0.0;
  for (i = 0; i < DIM; i++) {
    gradP_mag += grad_II_P[i] * grad_II_P[i];
  }

  gradP_mag = sqrt(gradP_mag);

  if (gradP_mag > epsilon) {
    for (i = 0; i < DIM; i++) {
      gradP_tangent[i] = grad_II_P[i] / (gradP_mag);
    }
  }

  else {
    for (i = 0; i < DIM; i++) {
      gradP_tangent[i] = 0.0;
    }
  }

  /* Find tangent vector perpendicular to pressure gradient */

  gradP_normal_init[0] = -gradP_tangent[1];
  gradP_normal_init[1] = gradP_tangent[0];
  gradP_normal_init[2] = gradP_tangent[2];

  for (i = 0; i < DIM; i++) {
    gradP_normal[i] = 0.0;
    for (j = 0; j < DIM; j++) {
      gradP_normal[i] += (gradP_normal_init[j] * delta(i, j) -
                          gradP_normal_init[j] * fv->snormal[i] * fv->snormal[j] -
                          gradP_normal_init[j] * gradP_tangent[i] * gradP_tangent[j]);
    }
  }

  /* Extract wall heights */
  H = height_function_model(&H_U, &dH_U_dtime, &H_L, &dH_L_dtime, dH_U_dX, dH_L_dX, &dH_U_dp,
                            &dH_U_ddH, dH_dF, time, dt);

  /******* CALCULATE FLOW RATE AND AVERAGE VELOCITY ***********/

  /* Evaluate flow rate and average velocity based on viscosity model */

  if (gn->ConstitutiveEquation == POWER_LAW) {
    mu0 = gn->mu0;
    nexp = gn->nexp;

    for (i = 0; i < DIM; i++) {
      q[i] = 0.0;

      q[i] +=
          (mu0 / (gradP_mag + epsilon)) * (veloU[i] * pow(fabs(shear_top), nexp - 1.) * shear_top -
                                           veloL[i] * pow(fabs(shear_bot), nexp - 1.) * shear_bot);
      q[i] += -(mu0 * mu0) / pow(gradP_mag + epsilon, 2) * gradP_tangent[i] * nexp /
              (2. * nexp + 1.) *
              (pow(fabs(shear_top), 2. * nexp - 2.) * pow(shear_top, 3) -
               pow(fabs(shear_bot), 2. * nexp - 2.) * pow(shear_bot, 3));
      q[i] += (mu0 * cross_shear) / pow(gradP_mag + epsilon, 2) * gradP_normal[i] * nexp /
              (nexp + 1.) *
              (pow(fabs(shear_top), nexp - 1.) * shear_top * shear_top -
               pow(fabs(shear_bot), nexp - 1.) * shear_bot * shear_bot);
    }

    for (i = 0; i < DIM; i++) {
      v_avg[i] = 0.0;
      v_avg[i] +=
          (mu0 / (gradP_mag * H)) * (veloU[i] * pow(fabs(shear_top), nexp - 1.) * shear_top -
                                     veloL[i] * pow(fabs(shear_bot), nexp - 1.) * shear_bot);
      v_avg[i] += -(mu0 * mu0) / (2. * gradP_mag * gradP_mag * H) * gradP_tangent[i] * (2. * nexp) /
                  (2. * nexp + 1.) *
                  (pow(fabs(shear_top), 2. * nexp - 2.) * pow(shear_top, 3) -
                   pow(fabs(shear_bot), 2. * nexp - 2.) * pow(shear_bot, 3));
      v_avg[i] += (mu0 * cross_shear) / (gradP_mag * gradP_mag * H) * gradP_normal[i] * nexp /
                  (nexp + 1.) *
                  (pow(fabs(shear_top), nexp - 1.) * shear_top * shear_top -
                   pow(fabs(shear_bot), nexp - 1.) * shear_bot * shear_bot);
    }
  }

  else {
    GOMA_EH(GOMA_ERROR, "Not a supported constitutive equation ");
  }

  for (i = 0; i < DIM; i++) {
    LubAux->q[i] = q[i];
    LubAux->v_avg[i] = v_avg[i];
    LubAux->gradP_mag = gradP_mag;
    LubAux->gradP_tangent[i] = gradP_tangent[i];
    LubAux->gradP_normal[i] = gradP_normal[i];
  }

  return;

} /* End of calculate_lub_q_v_nonnewtonian */

void calculate_lub_q_v_nonnewtonian_sens(double mu,
                                         double H,
                                         double veloU[DIM],
                                         double veloL[DIM],
                                         double grad_II_P[DIM],
                                         double phi_j,
                                         double grad_phi_j[DIM],
                                         double shear_top_plus,
                                         double shear_top_minus,
                                         double shear_bot_plus,
                                         double shear_bot_minus,
                                         double eps_top,
                                         double eps_bot)
/******************************************************************************
 *
 * calculate_lub_q_v_nonnewtonian_sens()
 *
 * Function to calculate flow rate per unit width (q) and average velocity (v)
 * of non Newtonian liquid's sensitivities w.r.t. pretty much anything you can think of
 *
 *
 * Kris Tjiptowidjojo tjiptowi@unm.edu
 *
 ******************************************************************************/
{
  int i, j;
  dbl shear_top, shear_bot, cross_shear;
  dbl gradP_mag, gradP_tangent[DIM];
  dbl gradP_normal[DIM];
  dbl dq_dP[DIM];
  dbl dq_dshear_top[DIM], dq_dshear_bot[DIM], dq_dcross_shear[DIM];
  dbl dgradP_mag_dP;
  dbl dgradP_tangent_dP[DIM], dgradP_normal_dP[DIM], dgradP_normal_init_dP[DIM];

  dbl q_plus[DIM], q_minus[DIM];

  dbl epsilon = 1.0e-5;

  shear_top = fv->sh_shear_top;     /* Top wall shear rate */
  shear_bot = fv->sh_shear_bot;     /* Bottom wall shear rate */
  cross_shear = fv->sh_cross_shear; /* Cross stream shear stress */

  for (i = 0; i < DIM; i++) {
    gradP_tangent[i] = LubAux->gradP_tangent[i];
    gradP_normal[i] = LubAux->gradP_normal[i];
  }

  gradP_mag = LubAux->gradP_mag;

  /******* CALCULATE PRESSURE GRADIENT AND RELATED VECTORS SENSITIVITIES ***********/

  /*Evaluate pressure gradient magnitude sensitivity w.r.t. pressure */

  dgradP_mag_dP = 0.0;
  for (i = 0; i < DIM; i++) {
    dgradP_mag_dP += grad_II_P[i] * grad_phi_j[i];
  }

  if (gradP_mag > epsilon) {
    dgradP_mag_dP *= 1. / gradP_mag;
  } else {
    dgradP_mag_dP = 0.0;
  }

  /*Evaluate pressure gradient tangent vector sensitivity w.r.t. pressure */

  if (gradP_mag > epsilon) {
    for (i = 0; i < DIM; i++) {
      dgradP_tangent_dP[i] =
          grad_phi_j[i] / (gradP_mag)-grad_II_P[i] * dgradP_mag_dP / pow(gradP_mag, 2);
    }
  } else {
    for (i = 0; i < DIM; i++) {
      dgradP_tangent_dP[i] = 0.0;
    }
  }

  /*Evaluate pressure gradient perpendicular vector sensitivity w.r.t. pressure */

  dgradP_normal_init_dP[0] = -dgradP_tangent_dP[1];
  dgradP_normal_init_dP[1] = dgradP_tangent_dP[0];
  dgradP_normal_init_dP[2] = dgradP_tangent_dP[2];

  if (gradP_mag > epsilon) {
    for (i = 0; i < DIM; i++) {
      dgradP_normal_dP[i] =
          dgradP_normal_init_dP[i] -
          gradP_tangent[i] *
              (2. * grad_II_P[2] * grad_phi_j[2] / pow(gradP_mag, 2) -
               2. * grad_II_P[2] * grad_II_P[2] * dgradP_mag_dP / pow(gradP_mag, 3)) -
          (grad_II_P[2] * grad_II_P[2]) * dgradP_tangent_dP[i];
    }
  }

  else {
    for (i = 0; i < DIM; i++) {
      dgradP_normal_dP[i] = 0.0;
    }
  }

  /******* CALCULATE FLOW RATE SENSITIVITIES ***********/

  /* Power Law */
  if (gn->ConstitutiveEquation == POWER_LAW) {
    dbl mu0 = gn->mu0;
    dbl nexp = gn->nexp;

    /*Evaluate flowrate sensitivity w.r.t. top wall shear rate */

    for (i = 0; i < DIM; i++) {
      q_plus[i] = 0.0;

      q_plus[i] += (mu0 / (gradP_mag + epsilon)) *
                   (veloU[i] * pow(fabs(shear_top_plus), nexp - 1.) * shear_top_plus -
                    veloL[i] * pow(fabs(shear_bot), nexp - 1.) * shear_bot);
      q_plus[i] += -(mu0 * mu0) / pow(gradP_mag + epsilon, 2) * gradP_tangent[i] * nexp /
                   (2. * nexp + 1.) *
                   (pow(fabs(shear_top_plus), 2. * nexp - 2.) * pow(shear_top_plus, 3) -
                    pow(fabs(shear_bot), 2. * nexp - 2.) * pow(shear_bot, 3));
      q_plus[i] += (mu0 * cross_shear) / pow(gradP_mag + epsilon, 2) * gradP_normal[i] * nexp /
                   (nexp + 1.) *
                   (pow(fabs(shear_top_plus), nexp - 1.) * shear_top_plus * shear_top_plus -
                    pow(fabs(shear_bot), nexp - 1.) * shear_bot * shear_bot);
    }

    for (i = 0; i < DIM; i++) {
      q_minus[i] = 0.0;

      q_minus[i] += (mu0 / (gradP_mag + epsilon)) *
                    (veloU[i] * pow(fabs(shear_top_minus), nexp - 1.) * shear_top_minus -
                     veloL[i] * pow(fabs(shear_bot), nexp - 1.) * shear_bot);
      q_minus[i] += -(mu0 * mu0) / pow(gradP_mag + epsilon, 2) * gradP_tangent[i] * nexp /
                    (2. * nexp + 1.) *
                    (pow(fabs(shear_top_minus), 2. * nexp - 2.) * pow(shear_top_minus, 3) -
                     pow(fabs(shear_bot), 2. * nexp - 2.) * pow(shear_bot, 3));
      q_minus[i] += (mu0 * cross_shear) / pow(gradP_mag + epsilon, 2) * gradP_normal[i] * nexp /
                    (nexp + 1.) *
                    (pow(fabs(shear_top_minus), nexp - 1.) * shear_top_minus * shear_top_minus -
                     pow(fabs(shear_bot), nexp - 1.) * shear_bot * shear_bot);
    }

    for (i = 0; i < DIM; i++) {
      dq_dshear_top[i] = (q_plus[i] - q_minus[i]) / eps_top;
    }

    /*Evaluate flowrate sensitivity w.r.t. bottom wall shear rate */

    for (i = 0; i < DIM; i++) {
      q_plus[i] = 0.0;

      q_plus[i] += (mu0 / (gradP_mag + epsilon)) *
                   (veloU[i] * pow(fabs(shear_top), nexp - 1.) * shear_top -
                    veloL[i] * pow(fabs(shear_bot_plus), nexp - 1.) * shear_bot_plus);
      q_plus[i] += -(mu0 * mu0) / pow(gradP_mag + epsilon, 2) * gradP_tangent[i] * nexp /
                   (2. * nexp + 1.) *
                   (pow(fabs(shear_top), 2. * nexp - 2.) * pow(shear_top, 3) -
                    pow(fabs(shear_bot_plus), 2. * nexp - 2.) * pow(shear_bot_plus, 3));
      q_plus[i] += (mu0 * cross_shear) / pow(gradP_mag + epsilon, 2) * gradP_normal[i] * nexp /
                   (nexp + 1.) *
                   (pow(fabs(shear_top), nexp - 1.) * shear_top * shear_top -
                    pow(fabs(shear_bot_plus), nexp - 1.) * shear_bot_plus * shear_bot_plus);
    }

    for (i = 0; i < DIM; i++) {
      q_minus[i] = 0.0;

      q_minus[i] += (mu0 / (gradP_mag + epsilon)) *
                    (veloU[i] * pow(fabs(shear_top), nexp - 1.) * shear_top -
                     veloL[i] * pow(fabs(shear_bot_minus), nexp - 1.) * shear_bot_minus);
      q_minus[i] += -(mu0 * mu0) / pow(gradP_mag + epsilon, 2) * gradP_tangent[i] * nexp /
                    (2. * nexp + 1.) *
                    (pow(fabs(shear_top), 2. * nexp - 2.) * pow(shear_top, 3) -
                     pow(fabs(shear_bot_minus), 2. * nexp - 2.) * pow(shear_bot_minus, 3));
      q_minus[i] += (mu0 * cross_shear) / pow(gradP_mag + epsilon, 2) * gradP_normal[i] * nexp /
                    (nexp + 1.) *
                    (pow(fabs(shear_top), nexp - 1.) * shear_top * shear_top -
                     pow(fabs(shear_bot_minus), nexp - 1.) * shear_bot_minus * shear_bot_minus);
    }

    for (i = 0; i < DIM; i++) {
      dq_dshear_bot[i] = (q_plus[i] - q_minus[i]) / eps_bot;
    }

    /*Evaluate flowrate sensitivity w.r.t. cross stream shear stress */

    for (i = 0; i < DIM; i++) {

      dq_dcross_shear[i] = (mu0 * phi_j) / pow(gradP_mag + epsilon, 2) * gradP_normal[i] * nexp /
                           (nexp + 1.) *
                           (pow(fabs(shear_top), nexp - 1.) * shear_top * shear_top -
                            pow(fabs(shear_bot), nexp - 1.) * shear_bot * shear_bot);
    }

    /*Evaluate flowrate sensitivity w.r.t. pressure */

    for (i = 0; i < DIM; i++) {

      dq_dP[i] = -mu0 / pow(gradP_mag + epsilon, 2) *
                 (veloU[i] * pow(fabs(shear_top), nexp - 1.) * shear_top -
                  veloL[i] * pow(fabs(shear_bot), nexp - 1.) * shear_bot) *
                 dgradP_mag_dP;

      dq_dP[i] += (2. * mu0 * mu0) / pow(gradP_mag + epsilon, 3) * nexp / (2. * nexp + 1.) *
                  (pow(fabs(shear_top), 2. * nexp - 2.) * pow(shear_top, 3) -
                   pow(fabs(shear_bot), 2. * nexp - 2.) * pow(shear_bot, 3)) *
                  dgradP_mag_dP * gradP_tangent[i];

      dq_dP[i] += -(mu0 * mu0) / pow(gradP_mag + epsilon, 2) * nexp / (2. * nexp + 1.) *
                  (pow(fabs(shear_top), 2. * nexp - 2.) * pow(shear_top, 3) -
                   pow(fabs(shear_bot), 2. * nexp - 2.) * pow(shear_bot, 3)) *
                  dgradP_tangent_dP[i];

      dq_dP[i] += -(2. * cross_shear * mu0) / pow(gradP_mag + epsilon, 3) * nexp / (nexp + 1.) *
                  (pow(fabs(shear_top), nexp - 1.) * shear_top * shear_top -
                   pow(fabs(shear_bot), nexp - 1.) * shear_bot * shear_bot) *
                  dgradP_mag_dP * gradP_normal[i];

      dq_dP[i] += (cross_shear * mu0) / pow(gradP_mag + epsilon, 2) * nexp / (nexp + 1.) *
                  (pow(fabs(shear_top), nexp - 1.) * shear_top * shear_top -
                   pow(fabs(shear_bot), nexp - 1.) * shear_bot * shear_bot) *
                  dgradP_normal_dP[i];
    }
  }

  else {
    GOMA_EH(GOMA_ERROR, "Not a supported constitutive equation ");
  }

  LubAux->dgradP_mag_dP = dgradP_mag_dP;

  for (i = 0; i < DIM; i++) {
    LubAux->dgradP_tangent_dP[i] = dgradP_tangent_dP[i];
    LubAux->dgradP_normal_dP[i] = dgradP_normal_dP[i];

    for (j = 0; j < ei[pg->imtrx]->dof[SHELL_SHEAR_TOP]; j++) {
      LubAux->dq_dshear_top[i][j] = dq_dshear_top[i];
    }
    for (j = 0; j < ei[pg->imtrx]->dof[SHELL_SHEAR_BOT]; j++) {
      LubAux->dq_dshear_bot[i][j] = dq_dshear_bot[i];
    }
    for (j = 0; j < ei[pg->imtrx]->dof[SHELL_CROSS_SHEAR]; j++) {
      LubAux->dq_dcross_shear[i][j] = dq_dcross_shear[i];
    }
    LubAux->dq_dp2[i] = dq_dP[i];
  }

  return;

} /* End of calculate_lub_q_v_nonnewtonian_sens */

int lub_viscosity_integrate(const double strs,
                            const double H,
                            double *flow_mag,
                            double *dq_gradp,
                            double *dq_dh,
                            double *srate,
                            double *pre_P,
                            double *mu_star,
                            double *dq_dT)
/******************************************************************************
 *
 * lub_viscosity_integrate()
 *
 * Function to calculate viscosity integrals for lubrication flow
 *
 * Robert Secor
 *
 ******************************************************************************/
{
  double shr, shrw, vis_w = 1., visd = 0.;
  double nexp = gn->nexp, lam, aexp = gn->aexp, muinf;
  double yield = gn->tau_y, F, mu0, P_eps = gn->epsilon;
  double eps, res, TOL_CEIL = 1.e-6, res_tol, soln_tol;
  int iter, ITERMAX = 50, jdi, JDI_MAX = 25, ierr = 0, a_visc_type;
  double xint = 0., xint_a = 0., visc_a = 0., xintold = 0., temp, at = 1.;
  double xintdT = 0.;

  res_tol = MIN(TOL_CEIL, Epsilon[pg->imtrx][0]);
  soln_tol = MIN(TOL_CEIL, Epsilon[pg->imtrx][2]);

  if (mp->LubIntegrationModel == LUB_VISCINT_ANALYTICAL) {
    a_visc_type = TRUE;
  } else {
    a_visc_type = FALSE;
  }

  if (pd->gv[SHELL_TEMPERATURE]) {
    temp = fv->sh_t;
  } else {
    temp = upd->Process_Temperature;
  }

  if (DOUBLE_NONZERO(temp) && DOUBLE_NONZERO(mp->reference[TEMPERATURE])) {
    if (gn->ConstitutiveEquation == BINGHAM) {
      at = exp(gn->atexp * (1. / temp - 1. / mp->reference[TEMPERATURE]));
    } else if (gn->ConstitutiveEquation == CARREAU_WLF || gn->ConstitutiveEquation == BINGHAM_WLF) {
      at = exp(gn->atexp * (mp->reference[TEMPERATURE] - temp) /
               (gn->wlfc2 + temp - mp->reference[TEMPERATURE]));
    }
  }
  mu0 = at * gn->mu0;
  muinf = at * gn->muinf;
  lam = at * gn->lam;
  F = at * gn->fexp;

  /**  Take care of de-generate case first **/
  if (DOUBLE_ZERO(strs)) {
    switch (gn->ConstitutiveEquation) {
    case CARREAU:
    case CARREAU_WLF:
      vis_w = mu0;
      break;
    case BINGHAM:
    case BINGHAM_WLF:
      vis_w = mu0 + yield * F;
      break;
    default:
      GOMA_EH(GOMA_ERROR, "Missing Lub Viscosity model!");
    }
    *flow_mag = 0.;
    if (dq_gradp != NULL)
      *dq_gradp = -CUBE(H) / (12. * vis_w);
    if (pre_P != NULL)
      *pre_P = -CUBE(H) / (12. * vis_w);
    if (dq_dh != NULL)
      *dq_dh = 0.;
    if (srate != NULL)
      *srate = 0.;
    if (mu_star != NULL)
      *mu_star = vis_w;
    if (dq_dT != NULL)
      *dq_dT = 0.;
    return (ierr);
  }

  /** First iterate to find shearrate that corresponds to stress */
  shr = strs / mu0;
  for (iter = 0; iter < ITERMAX; iter++) {
    double xfact, tmp, tp1, tp2, tpe, tpe_d, xj, delta, P_sig;
    shrw = fabs(shr);
    tp1 = lam * shrw;
    xfact = 1. + pow(tp1, aexp);
    tmp = pow(xfact, (1. - nexp) / aexp);
    switch (gn->ConstitutiveEquation) {
    case CARREAU:
    case CARREAU_WLF:
      vis_w = muinf + (mu0 - muinf) / tmp;
      visd = (mu0 - muinf) * (nexp - 1.) * pow(tp1, aexp) / (tmp * xfact);
      break;
    case BINGHAM:
    case BINGHAM_WLF:
      tp2 = F * shrw;
      P_sig = pow(1. + tp2, P_eps);
      tpe = (1. - exp(-tp2)) / shrw * P_sig;
      tpe_d = (exp(-tp2) * (1. + tp2) - 1.) * P_sig / shrw + P_eps * tp2 * tpe / (1. + tp2);
      vis_w = muinf + (mu0 - muinf + yield * tpe) / tmp;
      visd = (mu0 - muinf + yield * tpe) * (nexp - 1.) * pow(tp1, aexp) / (tmp * xfact);
      visd += yield * tpe_d / tmp;
      break;
    default:
      GOMA_EH(GOMA_ERROR, "Missing Lub Viscosity model!");
    }

    res = vis_w * shr - strs;
    xj = vis_w + visd;
    if ((1. + visd / vis_w) < 0.01) {
      GOMA_WH(GOMA_ERROR, "Trouble... Viscosity function nearly singular!");
    }
    delta = -res / xj;
    if (iter < ITERMAX / 4) {
      shr += delta;
    } else if (iter < ITERMAX / 2) {
      shr += 0.2 * delta;
    } else if (iter < 3 * ITERMAX / 4) {
      shr += 0.5 * delta;
    } else {
      shr += delta;
    }
    eps = fabs(delta) / (1. + fabs(shr));
    if (fabs(res) < res_tol && eps < soln_tol)
      break;
  }
  if (eps > soln_tol) {
    if (eps < sqrt(soln_tol)) {
      ierr = -2;
      GOMA_WH(GOMA_ERROR, "Viscosity iteration sort of converged!");
    } else {
      int iconv = 1, ifail = 0, ifin = 0, icnt, ICNTMAX = 20;
      double th = 0., step = 0.1, shr0, shrd = 0.0, Fc, Yc;
      /** initiate continuation path instead **/
      shr0 = strs / mu0;
      for (icnt = 0; icnt < ICNTMAX; icnt++) {
        if ((th + step) >= 1.) {
          step = 1. - th;
          ifin = 1;
        }
        th += step;
        Fc = th * F;
        Yc = th * yield;
        shr = shr0 + step * shrd;
        for (iter = 0; iter < ITERMAX; iter++) {
          double xfact, tmp, tp1, tp2, tpe, tpe_d, xj, delta, P_sig;
          shrw = fabs(shr);
          tp1 = lam * shrw;
          xfact = 1. + pow(tp1, aexp);
          tmp = pow(xfact, (1. - nexp) / aexp);
          switch (gn->ConstitutiveEquation) {
          case CARREAU:
          case CARREAU_WLF:
            vis_w = muinf + (mu0 - muinf) / tmp;
            visd = (mu0 - muinf) * (nexp - 1.) * pow(tp1, aexp) / (tmp * xfact);
            break;
          case BINGHAM:
          case BINGHAM_WLF:
            tp2 = Fc * shrw;
            P_sig = pow(1. + tp2, P_eps);
            tpe = (1. - exp(-tp2)) / shrw * P_sig;
            tpe_d = (exp(-tp2) * (1. + tp2) - 1.) * P_sig / shrw + P_eps * tp2 * tpe / (1. + tp2);
            vis_w = muinf + (mu0 - muinf + Yc * tpe) / tmp;
            visd = (mu0 - muinf + Yc * tpe) * (nexp - 1.) * pow(tp1, aexp) / (tmp * xfact);
            visd += Yc * tpe_d / tmp;
            break;
          default:
            GOMA_EH(GOMA_ERROR, "Missing Lub Viscosity model!");
          }
          res = vis_w * shr - strs;
          xj = vis_w + visd;
          delta = -res / xj;
          if (iter < ITERMAX / 4) {
            shr += delta;
          } else if (iter < ITERMAX / 2) {
            shr += 0.2 * delta;
          } else if (iter < 3 * ITERMAX / 4) {
            shr += 0.5 * delta;
          } else {
            shr += delta;
          }
          eps = fabs(delta) / (1. + fabs(shr));
          if (fabs(res) < res_tol && eps < soln_tol) {
            iconv = 1;
            break;
          }
        } /* end of iteration loop*/
        if (iconv == 0) {
          ifin = 0;
          th -= step;
          step *= 0.1;
          ifail++;
          if (ifail >= 6) {
            GOMA_EH(GOMA_ERROR, "Max continuation failures!");
          }
        } else {
          ifail = 0;
          if (step != 0.0) {
            shrd = (shr - shr0) / step;
          } else {
            shrd = 0.0;
          }
          shr0 = shr;
          step *= pow(10., ((pow(2., 4 - iter) - 1.) / 4.));
        }
        if (ifin == 1)
          break;
      } /*  End of continuation Loop */
      if (!ifin || th < 1.) {
        ierr = -1;
        GOMA_EH(GOMA_ERROR, "Viscosity iteration not converged!");
      }
    }
  }

  /** Second compute viscosity integral (stationary walls) */
  shrw = fabs(shr);
  /** test analytical viscosity integration types **/
  if (a_visc_type) {
    switch (gn->ConstitutiveEquation) {
    case BINGHAM:
    case BINGHAM_WLF: {
      double shrF = 1. / F, shrY = pow(yield * pow(lam, 1. - nexp) / (mu0 - muinf), 1. / nexp);
      double tmp_pl3 = 1. / CUBE(lam * shrw);
      double tmp_cy1 = pow(1. + CUBE(lam * shrw), (2. + nexp) / 3.);
      double tmp_cyY = pow(1. + CUBE(lam * shrY), (2. * nexp) / 3.);
      double tmp_cy2 = pow(1. + CUBE(lam * shrw), (1. + 2. * nexp) / 3.);
      double tmp_cyY2 = pow(1. + CUBE(lam * shrY), (1. + 2. * nexp) / 3.);
      double shr1 = MIN(shrw, shrF), shr2 = MIN(shrw, shrY);

      xint_a = SQUARE(F * yield + mu0) * CUBE(shr1 / shrw) / 3.;
      if (shrw > shrF) {
        xint_a += (SQUARE(yield) * (shr2 - shrF) + mu0 * yield * (SQUARE(shr2) - SQUARE(shrF)) +
                   SQUARE(mu0) / 3. * (CUBE(shr2) - CUBE(shrF))) /
                  CUBE(shrw);
      }
      if (shrw > shrY) {
        xint_a += SQUARE(muinf) / 3. * (1. - CUBE(shrY / shrw)) +
                  2. * muinf * tmp_pl3 * (mu0 - muinf) / (2. + nexp) * (tmp_cy1 - tmp_cyY) +
                  tmp_pl3 * SQUARE(mu0 - muinf) / (1. + 2. * nexp) * (tmp_cy2 - tmp_cyY2);
      }
    } break;
    case CARREAU:
    case CARREAU_WLF: {
      double tmp_pl3 = 1. / CUBE(lam * shrw);
      double tmp_cy1 = pow(1. + CUBE(lam * shrw), (2. + nexp) / 3.) - 1.;
      double tmp_cy2 = pow(1. + CUBE(lam * shrw), (2. * nexp + 1.) / 3.) - 1.;
      xint_a = SQUARE(muinf) / 3. + 2. * muinf * tmp_pl3 * (mu0 - muinf) / (2. + nexp) * tmp_cy1 +
               tmp_pl3 * SQUARE(mu0 - muinf) / (1. + 2. * nexp) * tmp_cy2;
    } break;
    }
  }
  for (jdi = 0; jdi < JDI_MAX; jdi++) {
    double cee, x0, delx, vis = 1., jdiv, xfact, tmp, tpe, tp2, P_sig;
    int idiv, l;
    double tmp_pl;
    jdiv = pow(2., jdi);
    delx = 1. / jdiv;
    x0 = 0.0;
    xint = 0.;
    for (idiv = 0; idiv < jdiv; idiv++) {
      for (l = 0; l < mp->LubInt_NGP; l++) {
        cee = x0 + mp->Lub_gpts[l] * delx;
        xfact = 1. + pow(lam * cee * shrw, aexp);
        tmp = 1. / pow(xfact, (1. - nexp) / aexp);
        switch (gn->ConstitutiveEquation) {
        case CARREAU:
        case CARREAU_WLF:
          vis = muinf + (mu0 - muinf) * tmp;
          break;
        case BINGHAM:
        case BINGHAM_WLF:
          tp2 = F * shrw;
          P_sig = pow(1. + tp2, P_eps);
          tpe = (1. - exp(-tp2)) / shrw * P_sig;
          vis = muinf + (mu0 - muinf + yield * tpe) * tmp;
          break;
        default:
          GOMA_EH(GOMA_ERROR, "Missing Lub Viscosity model!");
        }
        if (a_visc_type) {
          switch (gn->ConstitutiveEquation) {
          case BINGHAM:
          case BINGHAM_WLF: {
            double shrF = 1. / F;
            double shrY = pow(yield * pow(lam, 1. - nexp) / (mu0 - muinf), 1. / nexp);
            if (cee * shrw < shrF) {
              visc_a = F * yield + mu0;
            } else if (cee * shrw < shrY) {
              visc_a = mu0 + yield / (cee * shrw);
            } else {
              tmp_pl = 1. / pow(1. + CUBE(cee * lam * shrw), (1. - nexp) / 3.);
              visc_a = muinf + (mu0 - muinf) * tmp_pl;
            }
          } break;
          case CARREAU:
          case CARREAU_WLF:
            tmp_pl = 1. / pow(1. + CUBE(cee * lam * shrw), (1. - nexp) / 3.);
            visc_a = muinf + (mu0 - muinf) * tmp_pl;
            break;
          }
        }
        xint += (vis * vis - visc_a * visc_a) * SQUARE(cee) * delx * mp->Lub_wts[l];
      }
      x0 += delx;
    }
    xint += xint_a;
    xint /= SQUARE(vis_w);
    eps = fabs(xint - xintold);
    xintold = xint;
    if (eps < soln_tol)
      break;
  }
  if (eps > soln_tol) {
    ierr = -1;
    GOMA_EH(GOMA_ERROR, "Viscosity Integral not converged!");
  }

  /**  Compute temperature sensitivity integral **/
  if (dq_dT != NULL && pd->gv[SHELL_TEMPERATURE]) {
    double dlnat_dT = 0., dvis_dT = 0., dlnvis_w_dT = 0., dlnvis = 0.;
    double xfact, tmp, tp2, tpe, P_sig;
    if (DOUBLE_NONZERO(temp) && DOUBLE_NONZERO(mp->reference[TEMPERATURE])) {
      if (gn->ConstitutiveEquation == BINGHAM) {
        dlnat_dT = -gn->atexp / SQUARE(temp);
      } else if (gn->ConstitutiveEquation == CARREAU_WLF ||
                 gn->ConstitutiveEquation == BINGHAM_WLF) {
        dlnat_dT = -gn->atexp * gn->wlfc2 / SQUARE(gn->wlfc2 + temp - mp->reference[TEMPERATURE]);
      }
    }
    xintold = 0.;
    switch (gn->ConstitutiveEquation) {
    case CARREAU_WLF:
      dlnvis_w_dT = dlnat_dT * (vis_w + aexp * pow(lam * shrw, aexp));
      break;
    case BINGHAM:
    case BINGHAM_WLF:
      tp2 = F * shrw;
      P_sig = pow(1. + tp2, P_eps);
      tpe = (1. - exp(-tp2)) / shrw * P_sig;
      dlnvis_w_dT = dlnat_dT * (vis_w + (mu0 - muinf + yield * tpe) * aexp * pow(lam * shrw, aexp));
      break;
    default:
      GOMA_EH(GOMA_ERROR, "Missing Lub Viscosity model!");
    }
    for (jdi = 0; jdi < JDI_MAX; jdi++) {
      double cee, x0, delx, vis = 1., jdiv;
      int idiv, l;
      jdiv = pow(2., jdi);
      delx = 1. / jdiv;
      x0 = 0.0;
      xintdT = 0.;
      for (idiv = 0; idiv < jdiv; idiv++) {
        for (l = 0; l < mp->LubInt_NGP; l++) {
          cee = x0 + mp->Lub_gpts[l] * delx;
          xfact = 1. + pow(lam * cee * shrw, aexp);
          tmp = 1. / pow(xfact, (1. - nexp) / aexp);
          switch (gn->ConstitutiveEquation) {
          case CARREAU_WLF:
            vis = muinf + (mu0 - muinf) * tmp;
            dvis_dT = dlnat_dT * vis * (1. + nexp * pow(lam * cee * shrw, aexp)) / xfact;
            dlnvis = (mu0 - muinf) * (nexp - 1.) * tmp / xfact * pow(lam * cee * shrw, aexp);
            break;
          case BINGHAM:
          case BINGHAM_WLF:
            tp2 = F * cee * shrw;
            P_sig = pow(1. + tp2, P_eps);
            tpe = (1. - exp(-tp2)) / (cee * shrw) * P_sig;
            vis = muinf + (mu0 - muinf + yield * tpe) * tmp;
            break;
          default:
            GOMA_EH(GOMA_ERROR, "Missing Lub Viscosity model!");
          }
          xintdT += 2 * dvis_dT * (vis - dlnvis) * SQUARE(cee) * delx * mp->Lub_wts[l];
        }
        x0 += delx;
      }
      xintdT /= SQUARE(vis_w);
      eps = fabs(xintdT - xintold);
      xintold = xintdT;
      if (eps < soln_tol)
        break;
    }
    if (eps > soln_tol) {
      ierr = -1;
      GOMA_EH(GOMA_ERROR, "Temperature Sensitivity Integral not converged!");
    }
    xintdT += -2. / SQUARE(vis_w) * dlnvis_w_dT * xint;
    *dq_dT = 0.25 * SQUARE(H) * shrw * (xintdT * dlnvis_w_dT * (1. - xint));
  }

  /**  Compute flow magnitude  **/
  *flow_mag = -0.25 * SQUARE(H) * shr * (1.0 - xint);
  if (dq_gradp != NULL) {
    if (Include_Visc_Sens) {
      *dq_gradp = -0.25 * CUBE(H) * xint / vis_w;
    } else {
      *dq_gradp = -0.25 * CUBE(H) * xint / vis_w * (1. + visd / vis_w);
    }
  }
  if (pre_P != NULL)
    *pre_P = -0.125 * CUBE(H) * (1. - xint) / vis_w;
  if (dq_dh != NULL)
    *dq_dh = -0.5 * H * shr;
  if (srate != NULL)
    *srate = shr;
  if (mu_star != NULL)
    *mu_star = vis_w;
  return (ierr);

} /* End of lub_viscosity_integrate */
/**************************
     Lub viscosity function - hopefully faster than viscosity
 **************************/
double
lub_viscos_fcn(const struct Generalized_Newtonian *gn_local, const double gammadot, double *visd) {
  double visw = 0., temp, at = 1.;
  double nexp = gn_local->nexp, lam = gn_local->lam, aexp = gn_local->aexp, muinf = gn_local->muinf;
  double yield = gn_local->tau_y, F = gn_local->fexp, mu0 = gn_local->mu0,
         P_eps = gn_local->epsilon;
  double tp1 = 0., xfact = 0., tmp = 0., tp2 = 0., P_sig = 0., tpe = 0., tpe_d = 0.;

  if (pd->gv[SHELL_TEMPERATURE]) {
    temp = fv->sh_t;
  } else {
    temp = upd->Process_Temperature;
  }

  if (DOUBLE_NONZERO(temp) && DOUBLE_NONZERO(mp->reference[TEMPERATURE])) {
    if (gn_local->ConstitutiveEquation == BINGHAM) {
      at = exp(gn_local->atexp * (1. / temp - 1. / mp->reference[TEMPERATURE]));
    } else if (gn_local->ConstitutiveEquation == CARREAU_WLF ||
               gn_local->ConstitutiveEquation == BINGHAM_WLF) {
      at = exp(gn_local->atexp * (mp->reference[TEMPERATURE] - temp) /
               (gn_local->wlfc2 + temp - mp->reference[TEMPERATURE]));
    }
  }
  muinf = at * gn_local->muinf;
  mu0 *= at;
  lam = at * gn_local->lam;
  F = at * gn_local->fexp;

  switch (gn_local->ConstitutiveEquation) {
  case POWER_LAW:
    visw = mu0 * pow(gammadot, 1. / nexp);
    *visd = visw / (nexp * SQUARE(gammadot));
    /*if(*visd != NULL)  *visd = visw / (nexp * SQUARE(gammadot));*/
    break;
  case CARREAU:
  case CARREAU_WLF:
    tp1 = lam * gammadot;
    xfact = 1. + pow(tp1, aexp);
    tmp = pow(xfact, (1. - nexp) / aexp);
    visw = muinf + (mu0 - muinf) / tmp;
    *visd = (mu0 - muinf) * (nexp - 1.) * pow(tp1, aexp) / (tmp * xfact);
    break;
  case BINGHAM:
  case BINGHAM_WLF:
    tp1 = lam * gammadot;
    xfact = 1. + pow(tp1, aexp);
    tmp = pow(xfact, (1. - nexp) / aexp);
    tp2 = F * gammadot;
    P_sig = pow(1. + tp2, P_eps);
    tpe = (1. - exp(-tp2)) / gammadot * P_sig;
    visw = muinf + (mu0 - muinf + yield * tpe) / tmp;
    if (1) {
      *visd = (mu0 - muinf + yield * tpe) * (nexp - 1.) * pow(tp1, aexp) / (tmp * xfact);
      *visd += yield * tpe_d / tmp;
    }
    break;
  default:
    GOMA_EH(GOMA_ERROR, "Missing Lub Viscosity model!");
  }
  return (visw);
}
/**************************
     flow subroutine - for knife gaps including web velocity
 **************************/
int lub2D_flow2D(struct Generalized_Newtonian *gn_loc,
                 const double gradp[DIM],
                 const double ev[DIM],
                 const double dev_dpg[DIM][DIM],
                 double qfl[DIM],
                 double dqfl_dp[DIM][DIM],
                 double dqfl_dh[DIM],
                 const double h,
                 double *wrate,
                 const double veloL[DIM],
                 const double veloU[DIM],
                 int iguess,
                 double *wstrs,
                 double relax_rate) {

  double pgrad, visw, vis0, vis1, str_scale, pg_inv;
  double shr, dvis0, dvis1, dvis0dK, dvis1dK;
  double eps, res, TOL_CEIL = 1.e-6, res_tol, soln_tol;
  int i, j, iter, ITERMAX = 50, ierr, ICNTMAX = 50;
  double xint1, xint2, xint1dK, xint2dK, xint3;
  double r[3], xj[3][3], u[3];
  int colinear = 0;
  struct Generalized_Newtonian *dum_ptr, *gn_tmp;
  double epsxint = 0.001, epsvisc = 1.e-6;
  static double rate0, rate1, Kconst;
  double fside, fv[DIM], evelo;
  double rate00, rate10, Kconst0, rate0d, rate1d, Kconstd, step, th, th1;
  int icnt, ifail, iconv, ifin, sz;
  double xjnum[3][3], xjs[3][3];
  double rate0f, rate1f, shear, res0, resid, relax;
  double shear0, shear1, wrate1;
  double term1, term2, term3, term4;
  double dfv_dpg[DIM][DIM], dqfl_depg[DIM][DIM];

  res_tol = MIN(TOL_CEIL, Epsilon[pg->imtrx][0]);
  soln_tol = MIN(TOL_CEIL, Epsilon[pg->imtrx][2]);

  /*        check if the flow is only drag flow  */

  shr = 0.;
  pgrad = 0.;
  for (i = 0; i < DIM; i++) {
    shr += SQUARE(veloU[i] - veloL[i]);
    pgrad += SQUARE(gradp[i]);
  }
  shr = sqrt(shr) / h;
  pgrad = sqrt(pgrad);
  str_scale = pgrad * h;
  visw = lub_viscos_fcn(gn_loc, shr, NULL);
  vis0 = lub_viscos_fcn(gn_loc, fabs(rate0), NULL);
  vis1 = lub_viscos_fcn(gn_loc, fabs(rate1), NULL);

  if (fabs(vis1 - vis0) / visw < res_tol) {
    double dqdp;
    dqdp = -CUBE(h) / (12. * visw);
    for (i = 0; i < DIM; i++) {
      qfl[i] = 0.5 * h * (veloU[i] + veloL[i]) + dqdp * gradp[i];
    }
    *wrate = MAX(rate1, rate0);
    *wstrs = visw * (*wrate);
    return (1);
  }

  fv[0] = -ev[1];
  fv[1] = ev[0];
  fv[2] = -ev[2];
  dfv_dpg[0][0] = -dev_dpg[1][0];
  dfv_dpg[0][1] = -dev_dpg[1][1];
  dfv_dpg[0][2] = -dev_dpg[1][2];
  dfv_dpg[1][0] = dev_dpg[0][0];
  dfv_dpg[1][1] = dev_dpg[0][1];
  dfv_dpg[1][2] = dev_dpg[0][2];
  dfv_dpg[2][0] = -dev_dpg[2][0];
  dfv_dpg[2][1] = -dev_dpg[2][1];
  dfv_dpg[2][2] = -dev_dpg[2][2];
  fside = 0.;
  for (i = 0; i < DIM; i++) {
    fside += fv[i] * (veloU[i] - veloL[i]);
  }
  Kconst = visw / h * fside;

  /*        Check for colinearity  */
  colinear = fabs(fside) <= res_tol;

  sz = sizeof(struct Generalized_Newtonian *);
  gn_tmp = (struct Generalized_Newtonian *)array_alloc(1, 1, sz);
  if (gn_tmp == NULL) {
    GOMA_EH(GOMA_ERROR, "Problem getting memory for Generalized_Newtonian material");
  }
  memset(gn_tmp, 0, sz);

  dum_ptr = gn_loc;
  gn_tmp = gn_loc;
  gn_loc = dum_ptr;

  /*   solve for the shear rate at the top and bottom surface
           and crosswise shear stress value  */

  /*        initial guesses based solutions for Newtonian liquids  */

  evelo = 0.;
  for (i = 0; i < DIM; i++) {
    evelo += ev[i] * (veloU[i] - veloL[i]);
  }
  rate00 = evelo / h - h * pgrad / (2. * visw);
  rate10 = evelo / h + h * pgrad / (2. * visw);
  Kconst0 = Kconst;
  if (iguess) {
    rate0 = rate00;
    rate1 = rate10;
    Kconst = Kconst0;
  }
  rate0d = rate0 - rate00;
  rate1d = rate1 - rate10;
  Kconstd = Kconst - Kconst0;
  /********************************************
         perform continuation from newtonian to viscosity parameters
   ********************************************/
  iconv = 1;
  step = 1.0;
  ifail = 0;
  th = 0.;
  ifin = 0;
  th1 = 1.;
  for (icnt = 0; icnt < ICNTMAX; icnt++) {
    double deter, detinv;
    if (th + step >= th1) {
      step = th1 - th;
      ifin = 1;
    }
    th += step;
    gn_tmp->nexp = 1.0 + th * (gn_loc->nexp - 1.0);
    gn_tmp->aexp = 2.0 + th * (gn_loc->aexp - 2.0);
    gn_tmp->tau_y = th * gn_loc->tau_y;
    gn_tmp->atexp = th * gn_loc->atexp;
    rate0 = rate00 + step * rate0d;
    rate1 = rate10 + step * rate1d;
    Kconst = Kconst0 + step * Kconstd;
    for (iter = 0; iter < ITERMAX; iter++) {
      if (colinear) {
        Kconst = 0.;
        vis0 = lub_viscos_fcn(gn_tmp, fabs(rate0), &dvis0);
        dvis0dK = 0.;
        vis1 = lub_viscos_fcn(gn_tmp, fabs(rate1), &dvis1);
        dvis1dK = 0.;
      } else {
        rate0f = lub2D_crsrate(gn_tmp, rate0, Kconst, epsvisc);
        shear = sqrt(SQUARE(rate0) + SQUARE(rate0f));
        vis0 = lub_viscos_fcn(gn_tmp, shear, &dvis0);
        dvis0dK = Kconst * vis0 * dvis0 / (CUBE(vis0) + SQUARE(Kconst) * dvis0);
        rate1f = lub2D_crsrate(gn_tmp, rate1, Kconst, epsvisc);
        shear = sqrt(SQUARE(rate1) + SQUARE(rate1f));
        vis1 = lub_viscos_fcn(gn_tmp, shear, &dvis1);
        dvis1dK = Kconst * vis1 * dvis1 / (CUBE(vis1) + SQUARE(Kconst) * dvis1);
      }

      /*       evaluate first viscosity integral  */

      ierr = lub2D_viscint_2D(gn_tmp, rate0, rate1, Kconst, &xint1, &xint2, &xint1dK, &xint2dK,
                              epsxint);
      if (!ierr)
        GOMA_EH(GOMA_ERROR, "2D viscosity integral failed.");

      /*        evaluate residuals, jacobian entries  */

      r[0] = h * pgrad - vis1 * rate1 + vis0 * rate0;
      r[1] = pgrad * (ev[0] * (veloU[0] - veloL[0]) + ev[1] * (veloU[1] - veloL[1])) -
             (vis1 * SQUARE(rate1) - vis0 * SQUARE(rate0) - xint1);
      r[2] = 0.;

      /*      jacobian entries   */

      if (Include_Visc_Sens || iter > ((int)(0.5 / relax_rate))) {
        xj[0][0] = vis0 + SQUARE(rate0) * dvis0;
        xj[0][1] = -vis1 - SQUARE(rate1) * dvis1;
        xj[1][0] = rate0 * (vis0 + SQUARE(rate0) * dvis0);
        xj[1][1] = -rate1 * (vis1 + SQUARE(rate1) * dvis1);
        xj[0][2] = rate0 * dvis0dK - rate1 * dvis1dK;
        xj[1][2] = SQUARE(rate0) * dvis0dK - SQUARE(rate1) * dvis1dK + xint1dK;
      } else {
        xj[0][0] = vis0;
        xj[0][1] = -vis1;
        xj[1][0] = vis0 * rate0;
        xj[1][1] = -vis1 * rate1;
        xj[0][2] = 0.;
        xj[1][2] = 0.;
      }
      if (!colinear) {
        if (Include_Visc_Sens || iter > ((int)(0.5 / relax_rate))) {
          xj[2][0] = Kconst * (1. + SQUARE(rate0) * dvis0 / vis0);
          xj[2][1] = -Kconst * (1. + SQUARE(rate1) * dvis1 / vis1);
        } else {
          xj[2][0] = Kconst;
          xj[2][1] = -Kconst;
        }
        xj[2][2] = -(rate1 * (1. + log(vis1)) - rate0 * (1. + log(vis0)) - xint2) -
                   Kconst * (rate1 / vis1 * dvis1dK - rate0 / vis0 * dvis0dK - xint2dK);
        r[2] = pgrad * (fv[0] * (veloU[0] - veloL[0]) + fv[1] * (veloU[1] - veloL[1])) -
               Kconst * (rate1 * (1. + log(vis1)) - rate0 * (1. + log(vis0)) - xint2);
      } else {
        xj[2][2] = 0.;
      }

      /*       Check numerical Jacobian  */

      if (0) {
        int i1, j1;
        ierr = lub2D_numjac(h, pgrad, rate0, rate1, Kconst, epsxint, ev, fv, veloU, veloL, gn_tmp,
                            xjnum);
        for (i1 = 0; i1 < 3; i1++) {
          for (j1 = 0; j1 < 3; j1++) {
            xjs[i1][j1] = xj[i1][j1];
            xj[i1][j1] = xjnum[i1][j1];
          }
        }
      }

      /*        solve 3x3 linear system  */

      if (colinear) {
        deter = xj[0][0] * xj[1][1] - xj[1][0] * xj[0][1];
        detinv = 1. / deter;
        u[0] = detinv * (-r[0] * xj[1][1] + r[1] * xj[0][1]);
        u[1] = detinv * (r[0] * xj[1][0] - r[1] * xj[0][0]);
        u[2] = 0.;
      } else {
        deter = xj[0][0] * (xj[1][1] * xj[2][2] - xj[2][1] * xj[1][2]) -
                xj[1][0] * (xj[0][1] * xj[2][2] - xj[2][1] * xj[0][2]) +
                xj[2][0] * (xj[1][2] * xj[0][1] - xj[1][1] * xj[0][2]);
        detinv = 1. / deter;
        u[0] = detinv * ((-r[0]) * (xj[2][2] * xj[1][1] - xj[2][1] * xj[1][2]) +
                         (-r[1]) * (xj[2][1] * xj[0][2] - xj[2][2] * xj[0][1]) +
                         (-r[2]) * (xj[1][2] * xj[0][1] - xj[1][1] * xj[0][2]));
        u[1] = detinv * (-r[0] * (xj[2][0] * xj[1][2] - xj[2][2] * xj[1][0]) -
                         r[1] * (xj[2][2] * xj[0][0] - xj[2][0] * xj[0][2]) -
                         r[2] * (xj[1][0] * xj[0][2] - xj[1][2] * xj[0][0]));
        u[2] = detinv * (-r[0] * (xj[2][1] * xj[1][0] - xj[2][0] * xj[1][1]) -
                         r[1] * (xj[2][0] * xj[0][1] - xj[2][1] * xj[0][0]) -
                         r[2] * (xj[1][1] * xj[0][0] - xj[1][0] * xj[0][1]));
      }
      res = SQUARE(u[0]) + SQUARE(u[1]);
      resid = sqrt(SQUARE(r[0]) + SQUARE(r[1]) + SQUARE(r[2]));
      res0 = SQUARE(rate0) + SQUARE(rate1);
      eps = sqrt(res / (1. + res0)) + fabs(u[2] / str_scale);
      relax = MIN(1., (double)(iter)*relax_rate);
      rate0 += relax * u[0];
      rate1 += relax * u[1];
      Kconst += relax * u[2];
      if (iter >= ITERMAX / 4 && icnt >= 3) {
        DPRINTF(stderr, "flow2 icnt %d %d %g %g %g\n", icnt, iter, eps, resid, relax);
        if (!colinear) {
          if (0) {
            int i1, j1;
            for (i1 = 0; i1 < 3; i1++) {
              for (j1 = 0; j1 < 3; j1++) {
                DPRINTF(stderr, "numjac %d %d %g %g %g\n", i1, j1, xjs[i1][j1], xjnum[i1][j1],
                        (xjs[i1][j1] - xjnum[i1][j1]) / xjnum[i1][j1]);
              }
            }
          }
        }
      }
      if (eps <= soln_tol && resid < res_tol) {
        iconv = 1;
        break;
      }
      if (res > 1.0E+20) {
        iconv = 0;
        break;
      }
    }
    if (eps > soln_tol || resid > res_tol)
      iconv = 0;

    /*     if not converged : cut step, go back  */

    if (!iconv) {
      ifin = 0;
      th -= step;
      step *= 0.1;
      ifail += 1;
      if (ifail >= 6) {
        GOMA_EH(GOMA_ERROR, "flow2D has 6 consecutive failures");
      }
      break;
    } else {
      ifail = 0;
      if (DOUBLE_NONZERO(step)) {
        rate0d = (rate0 - rate00) / step;
        rate1d = (rate1 - rate10) / step;
        Kconstd = (Kconst - Kconst0) / step;
      } else {
        rate0d = 0.;
        rate1d = 0.;
        Kconstd = 0.;
      }
      rate00 = rate0;
      rate10 = rate1;
      Kconst0 = Kconst;
      if (iter < ITERMAX / 4)
        step *= 1.8;
      if (iter >= ITERMAX / 2)
        step *= 0.5;
    }
    if (ifin)
      break;
  }
  DPRINTF(stderr, "flow2 continuation loop failed %g pct\n", 100 * th);

  /*    solve for flowrate, knowing both wall shear rates and
             the value of the crosswise shear stress  */

  xint3 = lub2D_flowint_2D(gn_loc, rate0, rate1, Kconst, epsxint);

  pg_inv = 1. / pgrad;
  term1 = veloU[0] * vis1 * rate1 - veloL[0] * vis0 * rate0;
  term2 = SQUARE(vis1) * CUBE(rate1) - SQUARE(vis0) * CUBE(rate0) - xint3;
  term3 = Kconst * (vis1 * SQUARE(rate1) - vis0 * SQUARE(rate0) - xint1);
  term4 = veloU[1] * vis1 * rate1 - veloL[1] * vis0 * rate0;
  qfl[0] = term1 * pg_inv - SQUARE(pg_inv) * (0.5 * ev[0] * term2 + fv[0] * term3);
  qfl[1] = term4 * pg_inv - SQUARE(pg_inv) * (0.5 * ev[1] * term2 + fv[1] * term3);

  /*        calculate maximum wall shear rate   */

  shear0 = lub2D_crsrate(gn_loc, rate0, Kconst, epsvisc);
  shear1 = lub2D_crsrate(gn_loc, rate1, Kconst, epsvisc);
  *wrate = sqrt(SQUARE(shear0) + SQUARE(rate0));
  wrate1 = sqrt(SQUARE(shear1) + SQUARE(rate1));
  *wrate = MAX(wrate1, *wrate);
  *wstrs = MAX(wrate1 * vis1, (*wrate) * vis0);
  iguess = 0;

  /*  Jacobian elements  */
  if (dqfl_dp != NULL) {
    double dq_dgamma[DIM][DIM], dq_gradp[DIM], dgamma_gradp[DIM];
    double rate0dh, rate1dh;
    memset(dq_dgamma, 0.0, sizeof(double) * DIM * DIM);

    /*        compute derivatives wrt rate0,rate1
              - removing termvx from both numerator and denominator  */

    for (i = 0; i < DIM; i++) {
      dq_dgamma[i][0] = -veloL[i] * pg_inv +
                        (ev[i] * vis0 * SQUARE(rate0) + fv[i] * Kconst * rate0) * SQUARE(pg_inv);
      dq_dgamma[i][1] = veloU[i] * pg_inv -
                        (ev[i] * vis1 * SQUARE(rate1) + fv[i] * Kconst * rate1) * SQUARE(pg_inv);
      /*         derivatives wrt pgrad  */
      dq_gradp[i] = -term1 * SQUARE(pg_inv) + (ev[i] * term2 + 2. * fv[i] * term3) * CUBE(pg_inv);
      /*        derivatives wrt px,py through ex,ey,fx,fy   */
      for (j = 0; j < DIM; j++) {
        dqfl_depg[i][j] = -(0.5 * dev_dpg[i][j] * term2 + dfv_dpg[i][j] * term3) * CUBE(pg_inv);
      }
    }

    /*        derivatives of rates wrt pgrad   */

    dgamma_gradp[0] = (evelo - h * rate1) / (rate1 - rate0);
    dgamma_gradp[1] = (evelo - h * rate0) / (rate1 - rate0);
    dgamma_gradp[2] = 0.;

    for (i = 0; i < DIM; i++) {
      for (j = 0; j < DIM; j++) {
        dqfl_dp[i][j] = (dq_dgamma[i][j] * dgamma_gradp[j] + dq_gradp[i]) * ev[j] + dqfl_depg[i][j];
      }
    }
    /*   qxpx = (qxd0*pd0 + qxd1*pd1 + qxdp)*ex + qxdpx
       qxpy = (qxd0*pd0 + qxd1*pd1 + qxdp)*ey + qxdpy
       qypx = (qyd0*pd0 + qyd1*pd1 + qydp)*ex + qydpx
       qypy = (qyd0*pd0 + qyd1*pd1 + qydp)*ey + qydpy  */
    iguess = 0;
    rate1dh = pgrad * rate0 / (rate1 - rate0);
    rate0dh = pgrad * rate1 / (rate1 - rate0);
    dqfl_dh[0] = ev[0] * (dq_dgamma[0][0] * rate0dh + dq_dgamma[0][1] * rate1dh);
    dqfl_dh[1] = ev[1] * (dq_dgamma[1][0] * rate0dh + dq_dgamma[1][0] * rate1dh);
    dqfl_dh[2] = 0.;
  }

  return (1);
}
/*************************************
         iteration routine for finding cross shear rate
 ************************************/
double lub2D_crsrate(const struct Generalized_Newtonian *gn_local,
                     const double rate,
                     const double strs,
                     const double tolerance) {
  double shr, shrw, vis_w = 1., visd = 0.;
  double eps, res, xj, delta;
  int iter, ITERMAX = 50;

  /** First iterate to find shearrate that corresponds to stress */
  shr = strs / gn_local->mu0;
  for (iter = 0; iter < ITERMAX; iter++) {
    shrw = sqrt(SQUARE(rate) + SQUARE(shr));
    vis_w = lub_viscos_fcn(gn_local, shrw, &visd);

    res = strs - vis_w * shr;
    xj = -vis_w - SQUARE(shr) * visd;
    if ((1. + visd / vis_w) < 0.01) {
      GOMA_WH(GOMA_ERROR, "Trouble... Viscosity function nearly singular!");
    }
    delta = -res / xj;
    if (iter < ITERMAX / 4) {
      shr += delta;
    } else if (iter < ITERMAX / 2) {
      shr += 0.2 * delta;
    } else if (iter < 3 * ITERMAX / 4) {
      shr += 0.5 * delta;
    } else {
      shr += delta;
    }
    eps = fabs(delta) / (1. + fabs(shr));
    if (fabs(res) < tolerance && eps < tolerance) {
      break;
    }
  }
  if (eps > tolerance) {
    GOMA_EH(GOMA_ERROR, "2D Cross-rate iteration not converged!");
  }
  return (shr);
}
/***************************
     flow2D numerical jacobian for 3x3
 ****************************/
int lub2D_numjac(const double h,
                 const double pgrad,
                 const double rate0,
                 const double rate1,
                 const double Kconst,
                 const double epstol,
                 const double ev[DIM],
                 const double fv[DIM],
                 const double veloU[DIM],
                 const double veloL[DIM],
                 struct Generalized_Newtonian *gn_loc,
                 double xj[3][3]) {
  double u[3], u0[3], r_plus[3], r_neg[3], denom;
  double STEP = 0.000001;
  int j, i, ierr;
  double xint1, xint2;
  double rate0f, shear, vis0, rate1f, vis1, dvis0, dvis1;

  u0[0] = u[0] = rate0;
  u0[1] = u[1] = rate1;
  u0[2] = u[2] = Kconst;

  for (j = 0; j < 3; j++) {
    if (DOUBLE_NONZERO(u0[j])) {
      u[j] = u0[j] * (1. + STEP);
    } else {
      u[j] = u0[j] + STEP;
    }

    rate0f = lub2D_crsrate(gn_loc, u[0], u[2], epstol);
    shear = sqrt(SQUARE(u[0]) + SQUARE(rate0f));
    vis0 = lub_viscos_fcn(gn_loc, shear, &dvis0);
    rate1f = lub2D_crsrate(gn_loc, u[1], u[2], epstol);
    shear = sqrt(SQUARE(u[1]) + SQUARE(rate1f));
    vis1 = lub_viscos_fcn(gn_loc, shear, &dvis1);

    ierr = lub2D_viscint_2D(gn_loc, u[0], u[1], u[2], &xint1, &xint2, NULL, NULL, epstol);
    if (!ierr)
      GOMA_EH(GOMA_ERROR, "2D viscosity integral failed.");
    r_plus[0] = h * pgrad - vis1 * u[1] + vis0 * u[0];
    r_plus[1] = pgrad * (ev[0] * (veloU[0] - veloL[0]) + ev[1] * (veloU[1] - veloL[1])) -
                (vis1 * SQUARE(u[1]) - vis0 * SQUARE(u[0]) - xint1);
    r_plus[2] = pgrad * (fv[0] * (veloU[0] - veloL[0]) + fv[1] * (veloU[1] - veloL[1])) -
                u[2] * (u[1] * (1. + log(vis1)) - u[0] * (1. + log(vis0)) - xint2);

    if (DOUBLE_NONZERO(u0[j])) {
      u[j] = u0[j] * (1. - STEP);
    } else {
      u[j] = u0[j] - STEP;
    }

    rate0f = lub2D_crsrate(gn_loc, u[0], u[2], epstol);
    shear = sqrt(SQUARE(u[0]) + SQUARE(rate0f));
    vis0 = lub_viscos_fcn(gn_loc, shear, &dvis0);
    rate1f = lub2D_crsrate(gn_loc, u[1], u[2], epstol);
    shear = sqrt(SQUARE(u[1]) + SQUARE(rate1f));
    vis1 = lub_viscos_fcn(gn_loc, shear, &dvis1);
    ierr = lub2D_viscint_2D(gn_loc, u[0], u[1], u[2], &xint1, &xint2, NULL, NULL, epstol);
    if (!ierr)
      GOMA_EH(GOMA_ERROR, "2D viscosity integral failed.");
    r_neg[0] = h * pgrad - vis1 * u[1] + vis0 * u[0];
    r_neg[1] = pgrad * (ev[0] * (veloU[0] - veloL[0]) + ev[1] * (veloU[1] - veloL[1])) -
               (vis1 * SQUARE(u[1]) - vis0 * SQUARE(u[0]) - xint1);
    r_neg[2] = pgrad * (fv[0] * (veloU[0] - veloL[0]) + fv[1] * (veloU[1] - veloL[1])) -
               u[2] * (u[1] * (1. + log(vis1)) - u[0] * (1. + log(vis0)) - xint2);

    if (DOUBLE_NONZERO(u0[j])) {
      denom = 2. * u0[j] * STEP;
    } else {
      denom = 2. * STEP;
    }
    for (i = 0; i < 3; i++) {
      xj[i][j] = (r_plus[i] - r_neg[i]) / denom;
    }
    u[j] = u0[j];
  }
  return (1);
}
/**************************
    flow2D viscosity integrals
***************************/
int lub2D_viscint_2D(const struct Generalized_Newtonian *gn_loc,
                     const double rate0,
                     const double rate1,
                     const double Kconst,
                     double *xint1,
                     double *xint2,
                     double *xint1dK,
                     double *xint2dK,
                     const double tolerance) {
  double shr;
  double eps, res, res0;
  int jdi, JDIPOWER_MAX = 16, l;
  double drate, xintold, x0, idiv, jdiv, delx, cee;

  /*       evaluate first viscosity integral  */

  if (rate0 * rate1 >= 0.) {
    double rate, shear, dvisdK, vis, dvis;

    /*       rates of the same sign - do as one integral */

    drate = rate1 - rate0;
    xintold = 0.;
    for (jdi = 0; jdi < JDIPOWER_MAX; jdi++) {
      jdiv = pow(2., jdi);
      delx = 1. / jdiv;
      x0 = 0.0;
      *xint1 = *xint1dK = 0.;
      for (idiv = 0; idiv < jdiv; idiv++) {
        for (l = 0; l < mp->LubInt_NGP; l++) {
          cee = x0 + mp->Lub_gpts[l] * delx;
          rate = drate * cee + rate0;
          shear = lub2D_crsrate(gn_loc, rate, Kconst, tolerance);
          shr = sqrt(SQUARE(rate) + SQUARE(shear));
          vis = lub_viscos_fcn(gn_loc, shr, &dvis);
          *xint1 += vis * rate * drate * delx * mp->Lub_wts[l];
          dvisdK = Kconst * vis * dvis / (CUBE(vis) + SQUARE(Kconst) * dvis);
          *xint1dK += dvisdK * rate * drate * delx * mp->Lub_wts[l];
        }
        x0 += delx;
      }
      res = fabs(*xint1 - xintold);
      res0 = fabs(*xint1);
      eps = res / (1. + res0);
      if (eps < tolerance)
        break;
      xintold = *xint1;
    }
    if (eps > tolerance) {
      GOMA_EH(GOMA_ERROR, "xint1 same sign integration diverged!");
    }
  } else {
    double ratea, rateb, sheara, shearb, dvisdKa, dvisdKb;
    double visa, dvisa, visb, dvisb;

    /*       rates of opposite sign - do as two integrals  */

    xintold = 0.;
    for (jdi = 0; jdi < JDIPOWER_MAX; jdi++) {
      jdiv = pow(2., jdi);
      delx = 1. / jdiv;
      x0 = 0.;
      *xint1 = *xint1dK = 0.;
      for (idiv = 0; idiv < jdiv; idiv++) {
        for (l = 0; l < mp->LubInt_NGP; l++) {
          cee = x0 + mp->Lub_gpts[l] * delx;
          ratea = cee * rate1;
          sheara = lub2D_crsrate(gn_loc, ratea, Kconst, tolerance);
          shr = sqrt(SQUARE(ratea) + SQUARE(sheara));
          visa = lub_viscos_fcn(gn_loc, shr, &dvisa);
          rateb = cee * rate0;
          shearb = lub2D_crsrate(gn_loc, rateb, Kconst, tolerance);
          shr = sqrt(SQUARE(rateb) + SQUARE(shearb));
          visb = lub_viscos_fcn(gn_loc, shr, &dvisb);
          *xint1 += (visa * ratea * rate1 - visb * rateb * rate0) * delx * mp->Lub_wts[l];
          dvisdKa = Kconst * visa * dvisa / (CUBE(visa) + SQUARE(Kconst) * dvisa);
          dvisdKb = Kconst * visb * dvisb / (CUBE(visb) + SQUARE(Kconst) * dvisb);
          *xint1dK += (dvisdKa * ratea * rate1 - dvisdKb * rateb * rate0) * delx * mp->Lub_wts[l];
        }
        x0 += delx;
      }
      res = fabs(*xint1 - xintold);
      res0 = fabs(*xint1);
      eps = res / (1. + res0);
      if (eps < tolerance)
        break;
      xintold = *xint1;
    }
    if (eps > tolerance)
      GOMA_EH(GOMA_ERROR, "xint1 opposite sign integration diverged");
  }

  /*       evaluate second viscosity integrals  */

  if (rate0 * rate1 >= 0.) {
    double rate, shear, dvisdK, vis, dvis;

    /*       rates of the same sign - do as one integral */

    drate = rate1 - rate0;
    xintold = 0.;
    for (jdi = 0; jdi < JDIPOWER_MAX; jdi++) {
      jdiv = pow(2., jdi);
      delx = 1. / jdiv;
      x0 = 0.;
      *xint2 = *xint2dK = 0.;
      for (idiv = 0; idiv < jdiv; idiv++) {
        for (l = 0; l < mp->LubInt_NGP; l++) {
          cee = x0 + mp->Lub_gpts[l] * delx;
          rate = drate * cee + rate0;
          shear = lub2D_crsrate(gn_loc, rate, Kconst, tolerance);
          shr = sqrt(SQUARE(rate) + SQUARE(shear));
          vis = lub_viscos_fcn(gn_loc, shr, &dvis);
          *xint2 += log(vis) * drate * delx * mp->Lub_wts[l];
          dvisdK = Kconst * vis * dvis / (CUBE(vis) + SQUARE(Kconst) * dvis);
          *xint2dK += dvisdK / vis * drate * delx * mp->Lub_wts[l];
        }
        x0 += delx;
      }
      res = fabs(*xint2 - xintold);
      res0 = fabs(*xint2);
      eps = res / (1. + res0);
      if (eps < tolerance)
        break;
      xintold = *xint2;
    }
    if (eps > tolerance)
      GOMA_EH(GOMA_ERROR, "xint2 same sign integration diverged");
  } else {
    double ratea, rateb, sheara, shearb, dvisdKa, dvisdKb;
    double visa, dvisa, visb, dvisb;

    /*       rates of opposite sign - do as two integrals */

    xintold = 0.;
    for (jdi = 0; jdi < JDIPOWER_MAX; jdi++) {
      jdiv = pow(2., jdi);
      delx = 1. / jdiv;
      x0 = 0.0;
      *xint2 = *xint2dK = 0.;
      for (idiv = 0; idiv < jdiv; idiv++) {
        for (l = 0; l < mp->LubInt_NGP; l++) {
          cee = x0 + mp->Lub_gpts[l] * delx;
          ratea = cee * rate1;
          sheara = lub2D_crsrate(gn_loc, ratea, Kconst, tolerance);
          shr = sqrt(SQUARE(ratea) + SQUARE(sheara));
          visa = lub_viscos_fcn(gn_loc, shr, &dvisa);
          rateb = cee * rate0;
          shearb = lub2D_crsrate(gn_loc, rateb, Kconst, tolerance);
          shr = sqrt(SQUARE(rateb) + SQUARE(shearb));
          visb = lub_viscos_fcn(gn_loc, shr, &dvisb);
          *xint2 += (log(visa) * rate1 - log(visb) * rate0) * delx * mp->Lub_wts[l];
          dvisdKa = Kconst * visa * dvisa / (CUBE(visa) + SQUARE(Kconst) * dvisa);
          dvisdKb = Kconst * visb * dvisb / (CUBE(visb) + SQUARE(Kconst) * dvisb);
          *xint2dK += (dvisdKa / visa * rate1 - dvisdKb / visb * rate0) * delx * mp->Lub_wts[l];
        }
        x0 += delx;
      }
      res = fabs(*xint2 - xintold);
      res0 = fabs(*xint2);
      eps = res / (1. + res0);
      if (eps < tolerance)
        break;
      xintold = *xint2;
    }
    if (eps > tolerance)
      GOMA_EH(GOMA_ERROR, "xint2 opposite sign integration diverged");
  }
  return (1);
}
/***************************
     flow2D flow integral
****************************/
double lub2D_flowint_2D(const struct Generalized_Newtonian *gn_loc,
                        const double rate0,
                        const double rate1,
                        const double Kconst,
                        const double tolerance) {
  double drate, xintold, xint3, x0, delx, cee, res, res0;
  int l, jdi, jdiv, idiv, JDIPOWER_MAX = 16;

  /*       solve for flowrate, knowing both wall shear rates and
          the value of the crosswise shear stress   */

  if (rate0 * rate1 >= 0.) {
    double rate, shear, shr, vis;

    /*       rates of the same sign - do as one integral  */

    drate = rate1 - rate0;
    xintold = 0.;
    for (jdi = 0; jdi < JDIPOWER_MAX; jdi++) {
      jdiv = pow(2, jdi);
      delx = 1. / ((double)jdiv);
      xint3 = 0.;
      x0 = 0.;
      for (idiv = 0; idiv < jdiv; idiv++) {
        for (l = 0; l < mp->LubInt_NGP; l++) {
          cee = x0 + mp->Lub_gpts[l] * delx;
          rate = drate * cee + rate0;

          /*       determine cross rate which corresponds to this shearx  */

          shear = lub2D_crsrate(gn_loc, rate, Kconst, tolerance);
          shr = sqrt(SQUARE(rate) + SQUARE(shear));
          vis = lub_viscos_fcn(gn_loc, shr, NULL);
          xint3 += SQUARE(vis * rate) * drate * delx * mp->Lub_wts[l];
        }
      }
      res = fabs(xint3 - xintold);
      res0 = fabs(xint3);
      eps = res / (1. + res0);
      if (eps < tolerance)
        break;
      xintold = xint3;
    }
    if (eps > tolerance)
      GOMA_EH(GOMA_ERROR, "xint3 same sign integration diverged");
  } else {
    double ratea, rateb, sheara, shearb, shr, visa, visb;

    /*       rates of opposite sign - do as two integrals  */

    xintold = 0.;
    for (jdi = 0; jdi < JDIPOWER_MAX; jdi++) {
      jdiv = pow(2, jdi);
      delx = 1. / ((double)jdiv);
      xint3 = 0.;
      x0 = 0.;
      for (idiv = 0; idiv < jdiv; idiv++) {
        for (l = 0; l < mp->LubInt_NGP; l++) {
          cee = x0 + mp->Lub_gpts[l] * delx;
          ratea = cee * rate1;
          rateb = cee * rate0;

          /*       determine cross rate which corresponds to this shearx  */

          sheara = lub2D_crsrate(gn_loc, ratea, Kconst, tolerance);
          shearb = lub2D_crsrate(gn_loc, rateb, Kconst, tolerance);
          shr = sqrt(SQUARE(ratea) + SQUARE(sheara));
          visa = lub_viscos_fcn(gn_loc, shr, NULL);
          shr = sqrt(SQUARE(rateb) + SQUARE(shearb));
          visb = lub_viscos_fcn(gn_loc, shr, NULL);
          xint3 +=
              (SQUARE(visa * ratea) * rate1 - SQUARE(visb * rateb) * rate0) * delx * mp->Lub_wts[l];
        }
      }
      res = fabs(xint3 - xintold);
      res0 = fabs(xint3);
      eps = res / (1. + res0);
      if (eps < tolerance)
        break;
      xintold = xint3;
    }
    if (eps > tolerance)
      GOMA_EH(GOMA_ERROR, "xint3 opposite sign integration diverged");
  }
  return (xint3);
}
/*===========================================================================*/
