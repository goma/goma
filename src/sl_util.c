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

/* sl_util.c -- hi-level, one-shot-forall initialization; other Aztec stuff
 *
 *
 * Usage:
 *		sl_init(flag, a, ia, a2, ia2)
 *
 *	Various linear solver packages require initialization that describes
 * the size of the arrays and other initialization of pointers, etc. Examples
 * include the sparse 1.3 package of Kundert as well as the Aztec 2.1 package
 * from Hutchinson, Tuminaro & Shadid at SNL.
 *
 * This routine is called from rf_solve before entering the time stepping
 * loop and, hence, should provide a more unambiguous initialization for
 * time dependent problems that does not rely so much on static variables
 * to check for the "first call".
 *
 * This routine assumes the matrices will be described by either the
 * Modified Sparse Row format (MSR) or the Variable Block Row (VBR) format
 * with 0-based array indeces and names for the lowest numbered unknowns
 * and equations.
 *
 * The basic idea is that during the time stepping loop or during the
 * Newton iteration, a linear system of the form ... ?
 *
 *
 */

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "dpi.h"
#include "exo_struct.h"
#include "rf_allo.h"
#include "rf_fem.h"
#include "rf_io.h" /* for Debug_Flag */
#include "rf_mp.h"
#include "rf_solver.h"
#include "rf_solver_const.h"
#include "std.h"
#include "stdbool.h"

/*
 * This function has been updated to reflect the deprecation of Aztec
 * version 2.0 and earlier.
 */

#ifndef AZ_MAX_MEMORY_SIZE
#define AZ_MAX_MEMORY_SIZE 536870912 /* 512 MB (2^29) (*/
#endif

#ifdef PARALLEL
#ifndef MPI
#define MPI /* otherwise az_aztec.h trounces MPI_Request */
#endif
#endif

#include "sl_util.h" /* Variables of interest */

#include "az_aztec.h"
#include "dp_types.h"
#include "mm_as.h"
#include "mm_as_structs.h"
#include "mm_eh.h" /* Error handler. */
#include "sl_util_structs.h"

static int Num_Calls = 0;

/*
 * Initialize one or more Aztec matrix systems.
 *
 *	n = 1     Initialize the full Jacobian system.
 *
 *	n = 2	  Initialize the explict fill equation system.
 *
 *
 * Therefore,
 *
 *	n = 3	  Initialize both systems (1+2).
 *
 * Assumptions:   Aztec v. 2.1 or later is compiled in (now required
 *		  regardless of which linear solver is selected).
 *
 * Many parameters for Aztec can be set using the same boiler plate. The
 * few things that differ from one matrix system to another are marked
 * using the comment "Whoa, mule!" comment.
 *
 * If you want to add a 3rd, 4th, etc. linear system, it should be
 * straightforward. Find the next unused power of 2 as an associated place
 * in the n integer for your system. Also, it's up to you in the rf_solve.c
 * file to actually allocate space for val a.k.a. a and for bindx a.k.a. ija.
 *
 * The static variable First_Call is altered so the first call is treated
 * differently from others.
 */

void sl_init(unsigned int option_mask, /* option flag */
             struct GomaLinearSolverData *ams[],
             Exo_DB *exo,
             Dpi *dpi,
             Comm_Ex cx[]) {
  /* LOCAL VARIABLES */
  int i;
  int err;
  int length;
  int p; /* processor index */
  int Do_Jacobian;
  struct GomaLinearSolverData *A;
  int imtrx = pg->imtrx;

  Do_Jacobian = (option_mask & 1);

  if (Num_Calls >= upd->Total_Num_Matrices) {
    GOMA_EH(GOMA_ERROR, "Calls should match the number of matrices");
  } else {
    /* One of these blocks for each matrix system to be solved. */
    if (Do_Jacobian) {

      A = ams[imtrx]; /* Whoa, mule! */

      /*
       * Manual version of AZ_processor_info(). Should work for both
       * serial and parallel...
       */

      A->proc_config[AZ_node] = ProcID;
      A->proc_config[AZ_N_procs] = Num_Proc;
      A->proc_config[AZ_dim] = 0;

      /*
       * Manual version of AZ_read_update().
       */

      /*
       * Number of equations in the matrix that this processor owns.
       */

      A->N_update = num_internal_dofs[imtrx] + num_boundary_dofs[imtrx]; /* Whoa, mule! */
      A->update = NULL;

      AZ_defaults(A->options, A->params);

      /*
       * Any goma specific preferences for solver options.
       */

      set_aztec_options_params(A->options, A->params);

      /*
       * These arrays are normally set by AZ_transform(). For two difft
       * matrix systems of different sizes, however, it's not obvious that
       * routine will do the right thing. So, instead, just fill these
       * arrays manually for each matrix system, assuming serial solution.
       */

      if (Num_Proc == 1) {
        A->external = NULL; /* Serial ? Nothing to be sent. */
        A->update_index = NULL;
        /*	      A->update_index = A->update; */

        A->extern_index = NULL; /* A->external; */

        /* Need one for self, evidently. */
        A->data_org = (int *)array_alloc(1, AZ_COMM_SIZE + 0, sizeof(int));

        if (strcmp(Matrix_Format, "vbr") == 0) {
          A->data_org[AZ_matrix_type] = AZ_VBR_MATRIX;
          A->mat_type = AZ_VBR_MATRIX;
          A->data_org[AZ_N_int_blk] = dpi->num_internal_nodes;
          A->data_org[AZ_N_bord_blk] = dpi->num_boundary_nodes;
          A->data_org[AZ_N_ext_blk] = dpi->num_external_nodes;
        } else {
          A->data_org[AZ_matrix_type] = AZ_MSR_MATRIX;
          A->mat_type = AZ_MSR_MATRIX;
          A->data_org[AZ_N_int_blk] = A->N_update;
          A->data_org[AZ_N_bord_blk] = 0;
          A->data_org[AZ_N_ext_blk] = 0;
        }

        A->data_org[AZ_N_internal] = A->N_update;
        A->data_org[AZ_N_border] = 0;
        A->data_org[AZ_N_external] = 0;
        A->data_org[AZ_N_neigh] = 0;
        A->data_org[AZ_total_send] = 0;
        A->data_org[AZ_name] = 1 + imtrx; /* Whoa, mule! */
        A->data_org[AZ_neighbors] = 0;
        A->data_org[AZ_rec_length] = 0;
        A->data_org[AZ_send_length] = 0;
        /*	  A->data_org[AZ_send_list]   = 0; */

      } else /* parallel case... */
      {

        /*
         * None of these are needed if advanced planning orders the
         * local degrees of freedom such that:
         *	(1) the internal dofs occur first
         *	(2) the boundary dofs occur next (that other processors
         *	    might need, but this processor owns and updates.
         *	(3) the external dofs occur last, and that they
         *            are clumped together so that contiguous chunks
         *            are sent to neighboring processors.
         */

        /*
         * GOMA numbers dofs in contiguous chunks: internal, boundary,
         * external.  The
         * internal dofs possess monotone global dof names automatically
         * due to the way that dofs are numbered according to nodes
         * and the way that internal nodes for a proc are collected
         * in order from the global set.
         *
         * External dofs that are owned by other processors are
         * numbered contiguously. Aztec likes to have
         * each neighbor processor's dofs listed in contiguous chunks.
         *
         */

        length = AZ_COMM_SIZE + ptr_dof_send[imtrx][dpi->num_neighbors];

        A->data_org = (int *)smalloc(length * sizeof(int));

        A->data_org[AZ_N_internal] = num_internal_dofs[imtrx];
        A->data_org[AZ_N_border] = num_boundary_dofs[imtrx];
        A->data_org[AZ_N_external] = num_external_dofs[imtrx];

        if (strcmp(Matrix_Format, "vbr") == 0) {
          A->data_org[AZ_matrix_type] = AZ_VBR_MATRIX;
          A->mat_type = AZ_VBR_MATRIX;
          A->data_org[AZ_N_int_blk] = dpi->num_internal_nodes;
          A->data_org[AZ_N_bord_blk] = dpi->num_boundary_nodes;
          A->data_org[AZ_N_ext_blk] = dpi->num_external_nodes;
        } else {
          A->data_org[AZ_matrix_type] = AZ_MSR_MATRIX;
          A->mat_type = AZ_MSR_MATRIX;
          /*
           * For MSR, these are the required placebos to VBR info.
           */
          A->data_org[AZ_N_int_blk] = A->data_org[AZ_N_internal];
          A->data_org[AZ_N_bord_blk] = A->data_org[AZ_N_border];
          A->data_org[AZ_N_ext_blk] = A->data_org[AZ_N_external];
        }

        A->data_org[AZ_N_neigh] = dpi->num_neighbors;
        A->data_org[AZ_total_send] = ptr_dof_send[imtrx][dpi->num_neighbors];
        A->data_org[AZ_name] = 1 + imtrx; /* Whoa, mule! */

        /*
         * The names of my neighbor processors and how much to send
         * and receive for each one. You can tell that Aztec has
         * hardwired in a maximum number of neighbors here.
         */

        for (p = 0; p < dpi->num_neighbors; p++) {
          A->data_org[AZ_neighbors + p] = cx[p].neighbor_name;
          A->data_org[AZ_rec_length + p] = cx[p].num_dofs_recv;
          A->data_org[AZ_send_length + p] = cx[p].num_dofs_send;
        }

        /*
         * Names of dofs to send to those processors. Assume these
         * are indeces into the x[] vector with local meaning.
         */

        for (i = 0; i < ptr_dof_send[imtrx][dpi->num_neighbors]; i++) {
          A->data_org[AZ_send_list + i] = list_dof_send[imtrx][i];
        }
      }

      if (Debug_Flag > 0) {
        err = AZ_check_input(A->data_org, A->options, A->params, A->proc_config);

        if (err != 0)
          AZ_print_error(err);

        if (strcmp(Matrix_Format, "vbr") == 0) {
          AZ_check_vbr(A->N_update, A->N_external, AZ_GLOBAL, A->bindx, A->bpntr, A->cpntr,
                       A->rpntr, A->proc_config);
        } else {
          AZ_check_msr(A->bindx, A->N_update, A->N_external, AZ_GLOBAL, A->proc_config);
        }
      }
#ifdef MATRIX_DUMP
      A->Number_Jac_Dump = Number_Jac_Dump;
#endif /* MATRIX_DUMP */
    }
  }
  Num_Calls++;
} /* END of routine sl_init */
/******************************************************************************/

/*
 * sl_free() -- free up memory allocated to control Aztec soln process
 *
 *
 * Description:
 *		The pointers to Aztec_Linear_Solver_System structures are
 *		used to access arrays that were allocated during sl_init().
 *		These arrays are freed. The option flag is used to determine
 *		which linear systems need to have their arrays freed up.
 *
 *
 * Created: 1997/03/10 07:32 MST pasacki@sandia.gov
 *
 * Revised: 1997/08/23 16:30 MDT pasacki@sandia.gov
 */

void sl_free(unsigned int option_mask, struct GomaLinearSolverData *ams[]) {
  if (option_mask & 1) {
    free_ams(ams[JAC]);
  }

  return;
}
/******************************************************************************/
/* free_ams() -- free up dynamically allocated memory for an Aztec matrix system
 *
 * So far, the linear systems have shared enough similarities that this one
 * routine may be used for both, even for more than these two if you allocate
 * similar pieces.
 *
 * Created: 1997/08/23 16:42 MDT pasacki@sandia.gov
 */

void free_ams(struct GomaLinearSolverData *a) {
  safer_free((void **)&(a->data_org));
  safer_free((void **)&(a->val));
  safer_free((void **)&(a->indx));
  safer_free((void **)&(a->val_old));
  safer_free((void **)&(a->bindx));
  safer_free((void **)&(a->bpntr));
  safer_free((void **)&(a->cpntr));
  safer_free((void **)&(a->rpntr));
  return;
}

/* set_aztec_options_params -- set options prior to AZ_solve
 *
 * input:  int options[] -- should already have space enough to hold options.
 *			    valid options are described in the Aztec User's
 *			    Guide, SAND 95-1559.
 *
 *	   dbl params[]  -- should already have space enough to hold parameters.
 *			    valid options are described in the Aztec User's
 *			    Guide, SAND 95-1559.
 *
 * output: int options[] -- gets filled with options based on strings that
 *			    were collected from the input deck
 *
 *	   dbl params[]  -- gets filled with values based on strings from
 *			    the input deck.
 *
 * return: void
 *
 * Created: 1996/06/06 10:06 MDT pasacki@sandia.gov
 *
 * Revised: 1997/01/30 12:39 MST pasacki@sandia.gov
 * Revised: 1999/07/21 09:19 MDT pasacki@sandia.gov
 */

void set_aztec_options_params(int options[], double params[]) {
  /* LOCAL VARIABLES */
  int i;
  int iread;
  double f;
  /*
   * Set options & parameters for AZTEC linear solver package based upon
   * directions given in input deck and saved in various variables.
   *
   * For further information see "Aztec User's Guide", pp. 2-6
   * (SAND95-1559) Hutchinson, Shadid, Tuminaro.
   */

  /*
   * Solver
   */

  if (strcmp(Matrix_Solver, "cg") == 0) {
    options[AZ_solver] = AZ_cg;
  } else if (strcmp(Matrix_Solver, "gmres") == 0) {
    options[AZ_solver] = AZ_gmres;
  } else if (strcmp(Matrix_Solver, "cgs") == 0) {
    options[AZ_solver] = AZ_cgs;
  } else if (strcmp(Matrix_Solver, "tfqmr") == 0) {
    options[AZ_solver] = AZ_tfqmr;
  } else if (strcmp(Matrix_Solver, "bicgstab") == 0) {
    options[AZ_solver] = AZ_bicgstab;
  } else if (strcmp(Matrix_Solver, "y12m") == 0) {
    options[AZ_solver] = AZ_lu;
  } else if (strcmp(Matrix_Solver, "ma28") == 0) {
    Linear_Solver = MA28;
    options[AZ_solver] = -1;
  } else if (strcmp(Matrix_Solver, "lu") == 0) {
    Linear_Solver = SPARSE13a;
    options[AZ_solver] = -1;
  }
  /*IGBRK*/
  else if (strcmp(Matrix_Solver, "umf") == 0) {
    Linear_Solver = UMFPACK2;
    options[AZ_solver] = -1;
  }
  /*IGBRK*/
  else if (strcmp(Matrix_Solver, "umff") == 0) {
    Linear_Solver = UMFPACK2F;
    options[AZ_solver] = -1;
  } else if (strcmp(Matrix_Solver, "front") == 0) {
    Linear_Solver = FRONT;
    options[AZ_solver] = -1;
  } else if (strcmp(Matrix_Solver, "amesos") == 0) {
    Linear_Solver = AMESOS;
    options[AZ_solver] = -1;
  } else if (strcmp(Matrix_Solver, "stratimikos") == 0) {
    Linear_Solver = STRATIMIKOS;
    options[AZ_solver] = -1;
  } else if (strcmp(Matrix_Solver, "petsc") == 0) {
    Linear_Solver = PETSC_SOLVER;
    options[AZ_solver] = -1;
  } else if (strcmp(Matrix_Solver, "petsc_complex") == 0) {
    Linear_Solver = PETSC_COMPLEX_SOLVER;
    options[AZ_solver] = -1;
  } else if (strcmp(Matrix_Solver, "aztecoo") == 0) {
    Linear_Solver = AZTECOO;
    if (strcmp(AztecOO_Solver, "cg") == 0) {
      options[AZ_solver] = AZ_cg;
    } else if (strcmp(AztecOO_Solver, "gmres") == 0) {
      options[AZ_solver] = AZ_gmres;
    } else if (strcmp(AztecOO_Solver, "cgs") == 0) {
      options[AZ_solver] = AZ_cgs;
    } else if (strcmp(AztecOO_Solver, "tfqmr") == 0) {
      options[AZ_solver] = AZ_tfqmr;
    } else if (strcmp(AztecOO_Solver, "bicgstab") == 0) {
      options[AZ_solver] = AZ_bicgstab;
    } else if (strcmp(AztecOO_Solver, "y12m") == 0) {
      options[AZ_solver] = AZ_lu;
    }
  } else {
    fprintf(stderr, "Valid solvers: cg, gmres, cgs, tfqmr, bicgstab, y12m, lu, umf, umff\n");
    GOMA_EH(GOMA_ERROR, "Unknown linear solver specification.");
  }

  /*
   * Matrix scaling
   */

  if (strcmp(Matrix_Scaling, "none") == 0 || strcmp(Matrix_Scaling, "") == 0) {
    options[AZ_scaling] = AZ_none;
  } else if (strcmp(Matrix_Scaling, "Jacobi") == 0) {
    options[AZ_scaling] = AZ_Jacobi;
  } else if (strcmp(Matrix_Scaling, "BJacobi") == 0) {
    options[AZ_scaling] = AZ_BJacobi;
  } else if (strcmp(Matrix_Scaling, "row_sum") == 0) {
    options[AZ_scaling] = AZ_row_sum;
  } else if (strcmp(Matrix_Scaling, "sym_diag") == 0) {
    options[AZ_scaling] = AZ_sym_diag;
  } else if (strcmp(Matrix_Scaling, "sym_row_sum") == 0) {
    options[AZ_scaling] = AZ_sym_row_sum;
  } else {
    fprintf(stderr, "Valid scalings are: none, Jacobi, BJacobi, row_sum, sym_diag, sym_row_sum\n");
    GOMA_EH(GOMA_ERROR, "Unknown matrix scaling specification.");
  }

  /*
   * Preconditioner
   */
  options[AZ_precond] = AZ_none;

  if (strcmp(Matrix_Preconditioner, "none") == 0 || strcmp(Matrix_Preconditioner, "") == 0) {
    options[AZ_precond] = AZ_none;
  } else if (strcmp(Matrix_Preconditioner, "Jacobi") == 0) {
    options[AZ_precond] = AZ_Jacobi;
  } else if (strcmp(Matrix_Preconditioner, "Neumann") == 0) {
    options[AZ_precond] = AZ_Neumann;
  } else if (strcmp(Matrix_Preconditioner, "ls") == 0) {
    options[AZ_precond] = AZ_ls;
  } else if (strcmp(Matrix_Preconditioner, "sym_GS") == 0) {
    options[AZ_precond] = AZ_sym_GS;
  }

  /*
   * In Aztec 1.x these were preconditioners. In Aztec 2.x they are
   * options to AZ_subdomain_solve...
   */

  else if (strcmp(Matrix_Preconditioner, "dom_decomp") == 0) {
    options[AZ_precond] = AZ_dom_decomp;
  }

  /*
   * Strive for some backwards compatibility...
   */

  else if (strcmp(Matrix_Preconditioner, "lu") == 0) {
    options[AZ_precond] = AZ_dom_decomp;
    if (Debug_Flag > 1) {
      DPRINTF(stderr, "Deprecated Aztec 1.x usage of Preconditioner.  Prefer:\n");
      DPRINTF(stderr, "Preconditioner                  = dom_decomp\n");
      DPRINTF(stderr, "Matrix subdomain solver          = %s\n", Matrix_Preconditioner);
    }
    /*
     * Save for subsequent parsing...
     */

    strcpy(Matrix_Subdomain_Solver, Matrix_Preconditioner);
  }

  else if (strcmp(Matrix_Preconditioner, "ilu") == 0) {
    options[AZ_precond] = AZ_dom_decomp;
    if (Debug_Flag > 1) {
      DPRINTF(stderr, "Deprecated Aztec 1.x usage of Preconditioner.  Prefer:\n");
      DPRINTF(stderr, "Preconditioner                  = dom_decomp\n");
      DPRINTF(stderr, "Matrix subdomain solver          = %s\n", Matrix_Preconditioner);
    }
    /*
     * Save for subsequent parsing...
     */

    strcpy(Matrix_Subdomain_Solver, Matrix_Preconditioner);
  }

  else if (strcmp(Matrix_Preconditioner, "ilut") == 0) {
    options[AZ_precond] = AZ_dom_decomp;
    if (Debug_Flag > 1) {
      DPRINTF(stderr, "Deprecated Aztec 1.x usage of Preconditioner.  Prefer:\n");
      DPRINTF(stderr, "Preconditioner                  = dom_decomp\n");
      DPRINTF(stderr, "Matrix subdomain solver          = %s\n", Matrix_Preconditioner);
    }
    /*
     * Save for subsequent parsing...
     */

    strcpy(Matrix_Subdomain_Solver, Matrix_Preconditioner);
  }

  else if (strcmp(Matrix_Preconditioner, "bilu") == 0) {
    if (strcmp(Matrix_Format, "vbr"))
      GOMA_EH(-1, "BILU requires VBR matrix format!");
    options[AZ_precond] = AZ_dom_decomp;
    if (Debug_Flag > 1) {
      DPRINTF(stderr, "Deprecated Aztec 1.x usage of Preconditioner.  Prefer:\n");
      DPRINTF(stderr, "Preconditioner                  = dom_decomp\n");
      DPRINTF(stderr, "Matrix subdomain solver          = %s\n", Matrix_Preconditioner);
    }

    /*
     * Save for subsequent parsing...
     */

    strcpy(Matrix_Subdomain_Solver, Matrix_Preconditioner);
  }

  else if (strcmp(Matrix_Preconditioner, "rilu") == 0) {
    options[AZ_precond] = AZ_dom_decomp;
    if (Debug_Flag > 1) {
      DPRINTF(stderr, "Deprecated Aztec 1.x usage of Preconditioner.  Prefer:\n");
      DPRINTF(stderr, "Preconditioner                  = dom_decomp\n");
      DPRINTF(stderr, "Matrix subdomain solver          = %s\n", Matrix_Preconditioner);
    }

    /*
     * Save for subsequent parsing...
     */

    strcpy(Matrix_Subdomain_Solver, Matrix_Preconditioner);
  }

  else if (strcmp(Matrix_Preconditioner, "icc") == 0) {
    options[AZ_precond] = AZ_dom_decomp;
    if (Debug_Flag > 1) {
      DPRINTF(stderr, "Deprecated Aztec 1.x usage of Preconditioner.  Prefer:\n");
      DPRINTF(stderr, "Preconditioner                  = dom_decomp\n");
      DPRINTF(stderr, "Matrix subdomain solver          = %s\n", Matrix_Preconditioner);
    }

    /*
     * Save for subsequent parsing...
     */

    strcpy(Matrix_Subdomain_Solver, Matrix_Preconditioner);
  }

  else {
    fprintf(
        stderr,
        "\nValid Preconditioners:  none, Jacobi, Neumann, ls, sym_GS, dom_decomp (Aztec 2.1)\n");
    GOMA_EH(GOMA_ERROR, "Unknown Aztec preconditioner specification.");
  }

  if (options[AZ_precond] == AZ_dom_decomp) {
    if (strcmp(Matrix_Subdomain_Solver, "lu") == 0) {
      options[AZ_subdomain_solve] = AZ_lu;
    } else if (strcmp(Matrix_Subdomain_Solver, "ilut") == 0) {
      options[AZ_subdomain_solve] = AZ_ilut;
    } else if (strcmp(Matrix_Subdomain_Solver, "ilu") == 0) {
      options[AZ_subdomain_solve] = AZ_ilu;
    } else if (strcmp(Matrix_Subdomain_Solver, "rilu") == 0) {
      options[AZ_subdomain_solve] = AZ_rilu;
    } else if (strcmp(Matrix_Subdomain_Solver, "bilu") == 0) {
      if (strcmp(Matrix_Format, "vbr"))
        GOMA_EH(-1, "BILU requires VBR matrix format!");
/* NOTE: Trilinos does not support the AZ_bilu_ifp option at this time! */
#if 0
          options[AZ_subdomain_solve] = AZ_bilu_ifp;
#endif
      options[AZ_subdomain_solve] = AZ_bilu;
    } else if (strcmp(Matrix_Subdomain_Solver, "icc") == 0) {
      options[AZ_subdomain_solve] = AZ_icc;
    } else {
      fprintf(stderr, "Matrix subdomain solver = %s ?!?\n", Matrix_Subdomain_Solver);
      GOMA_EH(
          GOMA_ERROR,
          "Unknown specification. Valid Matrix subdomain solver = lu, ilu, ilut, rilu, bilu, icc");
    }
  }

  /*
   * Matrix residual norm expression
   */

  if (strcmp(Matrix_Residual_Norm_Type, "r0") == 0 || strcmp(Matrix_Residual_Norm_Type, "") == 0) {
    options[AZ_conv] = AZ_r0;
  } else if (strcmp(Matrix_Residual_Norm_Type, "rhs") == 0) {
    options[AZ_conv] = AZ_rhs;
  } else if (strcmp(Matrix_Residual_Norm_Type, "Anorm") == 0) {
    options[AZ_conv] = AZ_Anorm;
  } else if (strcmp(Matrix_Residual_Norm_Type, "noscaled") == 0) {
    options[AZ_conv] = AZ_noscaled;
  } else if (strcmp(Matrix_Residual_Norm_Type, "sol") == 0) {
    options[AZ_conv] = AZ_sol;
  } else if (strcmp(Matrix_Residual_Norm_Type, "weighted") == 0) {
    options[AZ_conv] = AZ_weighted;
    GOMA_EH(GOMA_ERROR, "Weighted norm unavailable; allocate bigger params[].");
  } else {
    fprintf(stderr, "Valid norms are: r0, rhs, Anorm, sol, (weighted)\n");
    GOMA_EH(GOMA_ERROR, "Unknown matrix residual norm type.");
  }

  /*
   * Matrix output type.
   */

  if (strcmp(Matrix_Output_Type, "all") == 0) {
    options[AZ_output] = AZ_all;
  } else if (strcmp(Matrix_Output_Type, "none") == 0) {
    options[AZ_output] = AZ_none;
  } else if (strcmp(Matrix_Output_Type, "last") == 0) {
    options[AZ_output] = AZ_last;
  } else if (strcmp(Matrix_Output_Type, "warnings") == 0) {
    options[AZ_output] = AZ_warnings;
  } else if (strcmp(Matrix_Output_Type, "") == 0) {
    options[AZ_output] = 1;
  } else if (sscanf(Matrix_Output_Type, "%d", &i) > 0) {
    if (i < 0) {
      GOMA_EH(GOMA_ERROR, "Matrix output type must specify a *positive* integer.");
    }
    options[AZ_output] = i;
  } else {
    fprintf(stderr, "Valid outputs are: all, none, warnings(2.x) last, <integer>\n");
    GOMA_EH(GOMA_ERROR, "Unknown matrix output type.");
  }

  /*
   * Factorization information from previous calls.
   *
   * This option will have a nontrivial setup on 2nd and subsequent
   * calls to AZ_solve().
   */

  options[AZ_pre_calc] = AZ_calc;

  /*
   * Matrix factorization reuse.
   */

  if (strcmp(Matrix_Factorization_Reuse, "calc") == 0 ||
      strcmp(Matrix_Factorization_Reuse, "") == 0) {
    options[AZ_pre_calc] = AZ_calc;
  } else if (strcmp(Matrix_Factorization_Reuse, "recalc") == 0) {
    options[AZ_pre_calc] = AZ_recalc;
  } else if (strcmp(Matrix_Factorization_Reuse, "reuse") == 0) {
    options[AZ_pre_calc] = AZ_reuse;
  } else {
    fprintf(stderr, "Valid factorization reuses are: calc, recalc, reuse\n");
    GOMA_EH(GOMA_ERROR, "Unknown factorization reuse.");
  }

  /*
   * Level of graph fill-in (Aztec 2.x) for incomplete factorizations...
   */

  iread = sscanf(Matrix_Graph_Fillin, "%d", &i);
  if (iread == 1) {
    if (i < 0) {
      GOMA_EH(GOMA_ERROR, "Matrix graph fillin should be positive.");
    }
    options[AZ_graph_fill] = i;
  } else {
    fprintf(stderr, "\nYou said \"%s\"\n", Matrix_Graph_Fillin);
    fprintf(stderr, "Found %d args\n", iread);
    fprintf(stderr, "Valid argument is a single nonnegative integer.\n");
    GOMA_EH(GOMA_ERROR, "Undecipherable Matrix graph fillin argument.");
  }

  /*
   * Maximum number of iterations.
   */

  iread = sscanf(Matrix_Maximum_Iterations, "%d", &i);
  if (iread == 1) {
    /*
     * This ought to be like "Maxits"
     */
    options[AZ_max_iter] = i;
  } else {
    fprintf(stderr, "\nYou said \"%s\"\n", Matrix_Maximum_Iterations);
    fprintf(stderr, "Found %d args\n", iread);
    fprintf(stderr, "Valid argument is a single positive integer.\n");
    GOMA_EH(GOMA_ERROR, "Undecipherable maximum linear solve iteration argument.");
  }

  /*
   * Polynomial order when using polynomial preconditioning.
   */

  if (sscanf(Matrix_Polynomial_Order, "%d", &i) == 1) {
    options[AZ_poly_ord] = i;
  }

  if (strcmp(Matrix_Factor_Overlap, "none") == 0 || strcmp(Matrix_Factor_Overlap, "") == 0) {
    options[AZ_overlap] = AZ_none;
  } else if (strcmp(Matrix_Factor_Overlap, "diag") == 0) {
    options[AZ_overlap] = AZ_diag;
  } else if (sscanf(Matrix_Factor_Overlap, "%d", &i) == 1) {
    if (i < 0) {
      GOMA_EH(GOMA_ERROR, "Matrix factorization overlap must be nonnegative.");
    }
    options[AZ_overlap] = i;
  } else {
    fprintf(stderr, "Valid Aztec 2.1 factor overlaps: none(1), diag, k(2).\n");
    GOMA_EH(GOMA_ERROR, "Unknown factor overlap specification.");
  }

  /*
   * Matrix overlap type (Aztec 2)
   */

  if (strcmp(Matrix_Overlap_Type, "standard") == 0) {
    options[AZ_type_overlap] = AZ_standard;
  } else if (strcmp(Matrix_Overlap_Type, "symmetric") == 0) {
    options[AZ_type_overlap] = AZ_symmetric;
  } else {
    fprintf(stderr, "Valid Matrix overlap types: standard, symmetric.\n");
    GOMA_EH(GOMA_ERROR, "Unknown Matrix overlap type specification.");
  }

  /*
   * Size of Krylov subspace.
   */

  if (sscanf(Matrix_Krylov_Subspace, "%d", &i) == 1) {
    options[AZ_kspace] = i;
  }

  /*
   * Orthogonalization scheme for GMRES.
   */

  if (strcmp(Matrix_Orthogonalization, "classic") == 0 ||
      strcmp(Matrix_Orthogonalization, "classical") == 0) {
    options[AZ_orthog] = AZ_classic;
  } else if (strcmp(Matrix_Orthogonalization, "modified") == 0) {
    options[AZ_orthog] = AZ_modified;
  } else {
    fprintf(stderr, "\nValid Aztec orthogonalizations: classic, modified.\n");
    GOMA_EH(GOMA_ERROR, "Unknown orthogonalization specification.");
  }

  /*
   * Auxiliary vector (reqd for some iterative schemes.)
   */

  if (strcmp(Matrix_Auxiliary_Vector, "resid") == 0) {
    options[AZ_aux_vec] = AZ_resid;
  } else if (strcmp(Matrix_Auxiliary_Vector, "rand") == 0) {
    options[AZ_aux_vec] = AZ_rand;
  } else {
    fprintf(stderr, "Valid auxiliary vectors: resid, rand.\n");
    GOMA_EH(GOMA_ERROR, "Unknown auxiliary vector specification.");
  }

  /*
   * Matrix factorization information save....(Aztec 2)
   */

  if (sscanf(Matrix_Factorization_Save, "%d", &i) == 1) {
    if (i < 0 || i > 1) {
      GOMA_EH(GOMA_ERROR, "Matrix factorization save must be 0 or 1.");
    }
    options[AZ_keep_info] = i;
  } else {
    fprintf(stderr, "Valid Matrix factorization save: 0, 1\n");
    GOMA_EH(GOMA_ERROR, "Unknown Matrix factorization save specification.");
  }

  /*
   * Matrix reorder -- none is the default, but Aztec 2.0 offers
   * an RCM that is frequently helpful.
   */

  if (strcmp(Matrix_Reorder, "none") == 0 || strcmp(Matrix_Reorder, "") == 0) {
    options[AZ_reorder] = 0;
  } else if (strcmp(Matrix_Reorder, "rcm") == 0) {
    options[AZ_reorder] = 1;
  } else {
    fprintf(stderr, "Valid Matrix reorder options: none, rcm.\n");
    GOMA_EH(GOMA_ERROR, "Unknown Matrix reorder specification.");
  }

  if (sscanf(Matrix_Convergence_Tolerance, "%lf", &f) == 1) {
    params[AZ_tol] = f;
  } else {
    fprintf(stderr, "Valid Residual Ratio Tolerance is 1 float.\n");
    GOMA_EH(GOMA_ERROR, "Unknown Residual Ratio Tolerance.");
  }

  if (sscanf(Matrix_Drop_Tolerance, "%lf", &f) == 1) {
    params[AZ_drop] = f;
  } else {
    fprintf(stderr, "Valid matrix drop tolerance is 1 float.\n");
    GOMA_EH(GOMA_ERROR, "Unknown matrix drop tolerance.");
  }

  if (sscanf(Matrix_ILUT_Fill_Factor, "%lf", &f) == 1) {
    params[AZ_ilut_fill] = f;
  } else {
    fprintf(stderr, "Valid Matrix ILUT fill factor is 1 float.\n");
    GOMA_EH(GOMA_ERROR, "Unknown Matrix ILUT fill factor specification.");
  }

  if (sscanf(Matrix_RILU_Relax_Factor, "%lf", &f) == 1) {
    params[AZ_omega] = f;
  } else {
    fprintf(stderr, "Valid Matrix RILU relax factor is 1 float.\n");
    GOMA_EH(GOMA_ERROR, "Unknown Matrix RILU relax factor specification.");
  }

#ifdef TRILINOS
  if (sscanf(Matrix_BILU_Threshold, "%lf", &f) == 1) {
    params[AZ_rthresh] = f;
    params[AZ_athresh] = f;
  } else {
    fprintf(stderr, "Valid Matrix BILU Threshold is 1 float.\n");
    GOMA_EH(GOMA_ERROR, "Unknown Matrix BILU Threshold specification.");
  }

  if (sscanf(Matrix_Relative_Threshold, "%lf", &f) == 1) {
    params[AZ_rthresh] = f;
  } else {
    fprintf(stderr, "Valid Matrix Relative Threshold is 1 float.\n");
    GOMA_EH(GOMA_ERROR, "Unknown Matrix Relative Threshold specification.");
  }

  if (sscanf(Matrix_Absolute_Threshold, "%lf", &f) == 1) {
    params[AZ_athresh] = f;
  } else {
    fprintf(stderr, "Valid Matrix Absolute Threshold is 1 float.\n");
    GOMA_EH(GOMA_ERROR, "Unknown Matrix Absolute Threshold specification.");
  }
#endif

  /* MMH
   * This was never set so I did it.
   *
   * I think Aztec sets this internally based on options[AZ_output].
   * Perhaps not sufficiently customizable? -PAS
   */
  options[AZ_print_freq] = 0;

  if (strcmp(Matrix_Subdomain_Solver, "bilu") == 0) {
    if (strcmp(Matrix_Format, "vbr") != 0) {
      GOMA_EH(-1, "bilu Matrix subdomain solver requires the vbr Matrix storage format.");
    }
#ifdef TRILINOS
    if (params[AZ_rthresh] < 0) {
      GOMA_EH(-1, "Matrix BILU Threshold must be a positive floating point number.");
    }
#endif
  }

  return;
} /* END of routine set_aztec_options_params */
/******************************************************************************/

/* dump_aztec_status -- print status after AZ_solve()
 *
 * input:       dbl status[] -- presumed to be full of interesting results
 *				subsequent to some call to AZ_solve()
 *
 * output:	[none]	     -- <stdout> gets informative interpretation of
 *				this status vector based on Aztec User's Guide.
 * return:	void
 *
 * Created:	1996/06/06 10:11 MDT pasacki@sandia.gov
 *
 * Revised:
 */

void dump_aztec_status(double status[]) {
  fprintf(stdout, "\n");
  fprintf(stdout, "Az_solve() summary\n");
  fprintf(stdout, "------------------\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "\tNumber of iterations taken by iterative method = %.0f\n", status[AZ_its]);
  fprintf(stdout, "\tReason for AZ_solve() termination:\n");

  switch ((int)status[AZ_why]) {
  case AZ_normal:
    fprintf(stdout, "\t\tConvergence criterion satisfied.\n");
    break;

  case AZ_param:
    fprintf(stdout, "\t\tUser requested option unavailable.\n");
    break;

  case AZ_breakdown:
    fprintf(stdout, "\t\tNumerical breakdown occurred.\n");
    break;

  case AZ_loss:
    fprintf(stdout, "\t\tLoss of precision occurred.\n");
    break;

  case AZ_maxits:
    fprintf(stdout, "\t\tMaximum iterations occurred w/o convergence.\n");
    break;

  default:
    fprintf(stdout, "\t\tI have not the slightest idea.\n");
    break;
  }

  fprintf(stdout, "\tTrue residual norm                             = %.3e\n", status[AZ_r]);

  fprintf(stdout, "\tTrue residual ratio                            = %.3e\n", status[AZ_scaled_r]);

  fprintf(stdout, "\tNorm final residual (or estimate)              = %.3e\n", status[AZ_rec_r]);

  return;
} /* END of routine dump_aztec_status */

/* hide_external() -- move external rows from MSR sparse matrix to hiding place
 *
 * The intent is to be able to use a bloated ija, a matrix to accept all of the
 * manipulations that goma performs and then, once the matrix solver is ready
 * to be called, to erase that part of the matrix corresponding to rows that
 * are owned by other processors.
 *
 * It is possible to use a smaller matrix to begin, but then there are many
 * routes where we need to check if the particular operation is attempting
 * to fix a matrix entry for an external dof. The approach herein requires
 * substantially less programming effort with some sacrifice in increased
 * memory usage for bloated (ija,a).
 *
 * Notes: The names of the rows and columns are assumed to begin at zero and
 *        end at one less than the number of rows and columns.
 *
 *        The save areas should be allocated large enough to hold the extra
 *        equations, if you want to save them. If the a_save area is NULL, then
 *        the entries for those equations are lost forever. If the ija_save area
 *        is NULL, then those column pointers are forever lost.
 *
 *        For goma applications, we'll discard the a matrix entries, but save
 *        the column pointers so that assembly during subsequent iterations
 *        proceeds smoothly.
 *
 *
 *
 * Created: 1997/10/13 13:29 MDT pasacki@sandia.gov
 *
 * Revised: 1997/10/30 13:08 MST pasacki@sandia.gov
 */

void hide_external(int n,     /* order of the original system      (in) */
                   int m,     /* order of the truncated system     (in) */
                   int *ija,  /* original column pointers          (in) */
                   int *ijas, /* save ija area */
                   double *a) /* original nonzero matrix values    (in) */
{
  int i;                /* temporary index */
  int d;                /* difference */
  int s;                /* start of saving of external off diag names */
  bool go_fast = FALSE; /* start slow, finish fast */

  /*
   * Standard checks for sanity.
   */

  if (n < 0) {
    GOMA_EH(GOMA_ERROR, "Original system is too small.");
  }

  if (m < 0) {
    GOMA_EH(GOMA_ERROR, "Truncated system is too small.");
  }

  if (m > n) {
    GOMA_EH(GOMA_ERROR, "Attempt to truncate to a LARGER system?!?");
  }

  if (ija == NULL || a == NULL) {
    GOMA_EH(GOMA_ERROR, "Either ija[] or a[] is empty.");
  }

  if (ija[0] != (n + 1)) {
    GOMA_EH(GOMA_ERROR, "Order n is inconsistent with ija[0] notion.");
  }

  if (m == n)
    return;

  /*
   * Save all of the original diag pointers. This includes the initial set of
   * pointers for rows 0,...,m-1 that will be modified to reflect the
   * collapse into a smaller matrix, as well as pointers to external
   * row's offdiagonal column names that will disappear entirely.
   */

  if (!go_fast) {
    for (i = 0; i < n + 1; i++) {
      ijas[i] = ija[i];
    }
  } else {
    memmove(ijas, ija, (n + 1) * sizeof(int));
  }

  /*
   * Save column names for off diagonal entries for external rows that
   * are being deleted. (Necessary?)
   */

  s = n + 1 - ija[m];

  if (!go_fast) {
    for (i = ija[m]; i < ija[n]; i++) {
      ijas[s + i] = ija[i];
    }
  } else {
    memmove(&ijas[n + 1], &ija[ija[m]], (ija[n] - ija[m]) * sizeof(int));
  }

  /*
   * Now collapse ija[] and a[] on top of the old diagonal terms. Don't
   * even bother with the high end values.
   */

  d = n - m;

  if (!go_fast) {
    for (i = n + 1; i < ija[m]; i++) {
      ija[i - d] = ija[i];
      a[i - d] = a[i];
    }
  } else {
    memmove(&ija[m + 1], &ija[n + 1], (ija[m] - (n + 1)) * sizeof(int));
    memmove(&a[m + 1], &a[n + 1], (ija[m] - (n + 1)) * sizeof(double));
  }

  /*
   * Now adjust the pointer values for rows 0,1,2,...,m-1 to reflect the
   * collapsed structure that starts at ija[m+1] instead of ija[n+1]
   */

  for (i = 0; i < m + 1; i++) {
    ija[i] -= d;
  }

  return;
}

/* show_external() -- restore external rows from hiding place to matrix
 *
 * The intent is to be able to use a bloated ija, a matrix to accept all of the
 * manipulations that goma performs and then, once the matrix solver is ready
 * to be called, to erase that part of the matrix corresponding to rows that
 * are owned by other processors.
 *
 * It is possible to use a smaller matrix to begin, but then there are many
 * routes where we need to check if the particular operation is attempting
 * to fix a matrix entry for an external dof. The approach herein requires
 * substantially less programming effort with some sacrifice in increased
 * memory usage for bloated (ija,a).
 *
 * Notes: The names of the rows and columns are assumed to begin at zero and
 *        end at one less than the number of rows and columns.
 *
 *        The save areas should be allocated large enough to hold the extra
 *        equations, if you want to save them. If the a_save area is NULL, then
 *        the entries for those equations are lost forever. If the ija_save area
 *        is NULL, then those column pointers are forever lost.
 *
 *        For goma applications, we'll discard the a matrix entries, but save
 *        the column pointers so that assembly during subsequent iterations
 *        proceeds smoothly.
 *
 *
 *
 * Created: 1997/10/13 13:29 MDT pasacki@sandia.gov
 *
 * Revised: 1997/11/03 09:15 MST pasacki@sandia.gov
 */

void show_external(int n,     /* order of the original system      (in) */
                   int m,     /* order of the truncated system     (in) */
                   int *ija,  /* original column pointers          (in) */
                   int *ijas, /* save ija area */
                   double *a) /* original nonzero matrix values    (in) */
{
  int i; /* temporary index */
  int d; /* difference */
  int s; /* start of saving of external off diag names */
  /*  bool go_fast=FALSE;		 start slow, finish fast */

  /*
   * Standard checks for sanity.
   */

  if (n < 0) {
    GOMA_EH(GOMA_ERROR, "Original system is too small.");
  }

  if (m < 0) {
    GOMA_EH(GOMA_ERROR, "Truncated system is too small.");
  }

  if (m > n) {
    GOMA_EH(GOMA_ERROR, "Attempt to truncate to a LARGER system?!?");
  }

  if (ija == NULL || a == NULL) {
    GOMA_EH(GOMA_ERROR, "Either ija[] or a[] is empty.");
  }

  if (ija[0] != (m + 1)) {
    GOMA_EH(GOMA_ERROR, "Order m is inconsistent with ija[0] notion.");
  }

  if (m == n)
    return;

  /*
   * 0. The diagonal pointers for rows 0,..,m-1 need to be increment since
   *    the column names will now begin at ija[n+1] instead of at ija[m+1].
   *
   *    This is unnecessary since we have the full old state of ija stored
   *    in ijas[].
   */

  d = n - m;

  /*
  for ( i=0; i<m+1; i++)
    {
      ija[i] += d;
    }
  */

  /*
   * 1. Make room for the old diagonal pointers for rows m,..,n-1 that
   *    were collapsed out. Move column names back up higher into ija[].
   */

  for (i = ija[m]; i > m; i--) {
    ija[i + d] = ija[i];
  }

  /*
   * 2. Restore the diagonal pointers for all rows.
   */

  for (i = 0; i < n + 1; i++) {
    ija[i] = ijas[i];
  }

  /*
   * 3. Restore the column names for high numbered rows. This may not be
   *    strictly necessary if the high numbered column names were unaffected
   *    by prior operations on ija[].
   *
   *    If a is still important for, say, quasi-Newton methods, then move
   *    back the entries to their appropriate places.
   */

  s = n + 1 - ija[m];

  for (i = ija[m]; i < ija[n]; i++) {
    ija[i] = ijas[s + i];
  }

  return;
}

void aztec_stringer(int why, double lits, char *stringer)
/*
 * Generates 3-character status string for result of AZ_solve call
 * If linear solve converged, this is the number of Aztec iterations taken.
 * Otherwise, it is a letter code which shows the reason for
 * non-convergence.
 */
{
  int nits = (int)lits;

  /*
   * This chunk of code was moved here from the various places
   * where AZ_solve is called.
   */
  strcpy(stringer, "   ");
  switch (why) {
  case AZ_normal:
    if (lits < 1000) {
      sprintf(stringer, "%3d", nits);
    } else if (lits < 10000) {
      sprintf(stringer, "%2dh", (int)(lits / 1e2));
    } else if (lits < 100000) {
      sprintf(stringer, "%2dk", (int)(lits / 1e3));
    } else {
      sprintf(stringer, "%3.0e", lits);
    }
    break;
  case AZ_param:
    strcpy(stringer, "bad");
    break;
  case AZ_breakdown:
    strcpy(stringer, "brk");
    break;
  case AZ_loss:
    strcpy(stringer, "los");
    break;
  case AZ_maxits:
    strcpy(stringer, "max");
    break;
  case AZ_ill_cond:
    strcpy(stringer, "ill");
    break;
  default:
    strcpy(stringer, "???");
    break;
  }

  return;
}

/*
 * Solves A x = b by Gaussian elimination with partial (row)
 * pivoting and row sum scaling.  A is dense.
 *
 * A is expected to be an M x row_size matrix, where M is at least
 * rank and no more than row_size.  rank specifies the size of the
 * problem you're solving, while row_size designates the actual
 * dimension of the storage being passed in A (so only the upper rank
 * x rank submatrix of A is used if rank < row_size).  b is the right
 * hand side, and x is returned as the solution.  b and x must be at
 * least length rank.  So, if solving A x = b, where x and b are of
 * length 2, but A happened to be allocated as a 3x3, then call
 * solve_NxN_system(&A[0][0], b, x, 2, 3).  Only the 2x2 part of A
 * would be involved.  On exit A and b are modified and x has the
 * solution.  It would not be hard to modify this to allow for
 * multiple right hand sides, but I didn't need it ain't here!
 *
 * Author: Matt Hopkins, 10/30/2002.
 */
void solve_NxN_system(dbl *A, dbl *b, dbl *x, const int rank, const int row_size) {
  int p[MAX_NXN_RANK]; /* indirection/permuation vector for pivoting. */
  int i, j, k;
  int pivot_index, index, index2;
  dbl pivot, factor, scale;

  if (rank > MAX_NXN_RANK) {
    fprintf(stderr, "You requested a dense matrix of size %d x %d.\n", rank, rank);
    GOMA_EH(GOMA_ERROR, "Increase MAX_NXN_RANK in sl_util.h.");
  }

  for (i = 0; i < rank; i++)
    p[i] = i;

  /* First, row sum scale all rows. */
  for (i = 0; i < rank; i++) {
    scale = 0.0;
    index = i * row_size;
    for (j = 0; j < rank; j++)
      scale += fabs(A[index++]);
    scale += fabs(b[i]);
    if (scale < 1.0e-14)
      GOMA_EH(GOMA_ERROR, "row sum < 1.0e-14.");
    index = i * row_size;
    for (j = 0; j < rank; j++)
      A[index++] /= scale;
    b[i] /= scale;
  }

  /* Now forward eliminate. */
  for (i = 0; i < rank - 1; i++) {
    /* Find pivot row, assuming it might be just #i */
    pivot_index = i;
    index = p[i] * row_size + i;
    pivot = fabs(A[index]);
    for (j = i + 1; j < rank; j++) {
      /* Looking at row #j */
      index = p[j] * row_size + i;
      if (fabs(A[index]) > pivot) {
        pivot_index = j;
        pivot = fabs(A[index]);
      }
    }

    /* Change row indirection if our pivot row isn't row #i. */
    if (pivot_index != i) /* strange, this was p[i] before -- I think that was wrong! */
    {
      j = p[i];
      p[i] = p[pivot_index];
      p[pivot_index] = j;
    }

    /* Perform elimination step #i */
    for (j = i + 1; j < rank; j++) {
      /* Get row multiplication factor */
      index = p[i] * row_size + i;
      index2 = p[j] * row_size + i;
      factor = A[index2] / A[index];
      A[index2] = 0.0;
      /* Carry through with the row multiplication. */
      for (k = i + 1; k < rank; k++)
        A[++index2] -= factor * A[++index];
      b[p[j]] -= factor * b[p[i]];
    }
  }

  /* And finally, back substitute. */
  for (i = rank - 1; i >= 0; i--) {
    /* Point to the last column in row p[i]. */
    index = p[i] * row_size + rank - 1;
    for (j = rank - 1; j > i; j--)
      b[p[i]] -= A[index--] * x[j];
    if (fabs(A[index]) < 1.0e-14)
      GOMA_EH(GOMA_ERROR, "A[index] < 1.0e-14.");
    x[i] = b[p[i]] / A[index];
  }
}

/* Alternate amesos_solve_msr() for builds without Trilinos & Amesos */

#if defined(GOMA_ENABLE_AMESOS) && defined(TRILINOS)
/* Use the function in sl_amesos_interface.C; do nothing here! */

#else
/* Just give an Amesos error message and abort! */
void amesos_solve_msr(
    char *choice, struct Aztec_Linear_Solver_System *ams, double *x_, double *b_, int flag) {
  fprintf(
      stderr,
      "Error: Need to compile with GOMA_ENABLE_AMESOS flag before using AMESOS solver packages.");
  fprintf(stderr, "   Also make sure appropriate libraries are linked for solver packages.");
  exit(-1);
}
#endif

/******************************************************************************/
/* END of file sl_util.c */
/******************************************************************************/
