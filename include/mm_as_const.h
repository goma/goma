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
 * mm_as_const.h:
 *
 *   This routine contains  constants used in the specification of
 *   properties in structures in mm_as_structs.h.
 *   They help in the assembly of equations and Jacobian terms
 */

/*
 *$Id: mm_as_const.h,v 5.2 2009-07-30 17:12:58 hkmoffa Exp $
 */

/*
 * T_{TERM} FLAGS
 * ----------------------------------------------------------------------
 *
 * These are likely to be the terms in the equations that we assemble.
 * This is our convention for using these flags:
 *
 *		[1] Say we are interested in the energy equation. If
 *                  we are interested in doing anything with it, then
 *                  if ( pd->e[pg->imtrx][R_ENERGY] )
 *                     {
 *                        ...stuff if energy eqn is active...
 *                     }
 *
 *		[2] If something is to be done for the energy equation,
 *                  then we can "and" this variable to see if the term
 *                  of interest is to be assembled.
 *
 *			if ( pd->e[pg->imtrx][R_ENERGY] & T_MASS )
 * 			   {
 *                            ... set up capacitance term of
 *                                energy equation...
 *			   }
 *
 *		[3] Term multipliers can be obtained via...
 *
 * 			f_em = pd->etm[pg->imtrx][R_ENERGY][LOG2_MASS]
 *
 *
 *
 * Thus, here are flags for turning on terms in an equation:
 */

#ifndef GOMA_MM_AS_CONST_H
#define GOMA_MM_AS_CONST_H

#define LOG2_SOMETHING    (0) /* Yes, bother to assemble. */
#define LOG2_MASS         (1) /* Mass matrix: d()/dt - v_mesh.grad() */
#define LOG2_ADVECTION    (2) /* Advection:   v_matl.grad()  */
#define LOG2_BOUNDARY     (3) /* Boundary:    n.grad()  */
#define LOG2_DIFFUSION    (4) /* Diffusion:   grad(phi_i).grad() */
#define LOG2_SOURCE       (5) /* Source:      phi_i r() */
#define LOG2_DIVERGENCE   (6) /* Divergence:  div() (continuity) */
#define LOG2_DIAG_FIX     (7) /* Strong diagonal (Dirichlet) */
#define LOG2_OFF_DIAG_FIX (8) /* Strong offdiagonal (Dirichlet) (DC) */
#define LOG2_POROUS_BRINK (9) /* Porous (Brinkman) term in N-S eqns */

#define T_NOTHING      (0)
#define T_SOMETHING    (1 << LOG2_SOMETHING)
#define T_MASS         (1 << LOG2_MASS)
#define T_ADVECTION    (1 << LOG2_ADVECTION)
#define T_BOUNDARY     (1 << LOG2_BOUNDARY)
#define T_DIFFUSION    (1 << LOG2_DIFFUSION)
#define T_SOURCE       (1 << LOG2_SOURCE)
#define T_DIVERGENCE   (1 << LOG2_DIVERGENCE)
#define T_DIAG_FIX     (1 << LOG2_DIAG_FIX)
#define T_OFF_DIAG_FIX (1 << LOG2_OFF_DIAG_FIX)
#define T_POROUS_BRINK (1 << LOG2_POROUS_BRINK)

#define T_ANYTHING (0x11111111) /*  */

/*
 *  V_TERM FLAGS
 * ----------------------------------------------------------------------
 *
 * These constants will be used in the v[] field in the
 * Problem_Description structure
 *
 * Sometimes we will want to solve weakly-coupled problems. Let's describe
 * the activity of a fixed-value variable, which is still spatially
 * varying, with these flags:
 * We will also use this variable to specify which variables are specific
 * to the current material. These variables may or may not be
 * discontinuous at material interfaces (it's a prerequisite though).
 */

#define V_NOTHING 0x00000000 /* Does not appear at all. */
#define V_SOLNVECTOR                               \
  0x00000001 /* This variable is solved for in the \
                solution vector */
#define V_FIXED                                       \
  0x00000002 /* Variable is fixed at constant value.  \
                It is assumed to be continuous across \
                material boundaries. */
#define V_SPECIFIED                                    \
  0x00000004 /* Variable is not part of the solution   \
                variable, but it does vary across the  \
                domain. Value is calculated via interp \
                from nodal values */
#define V_MATSPECIFIC                                    \
  0x00000008 /* This variable is unique to this material \
              * -> It will not be contiguous across      \
              * material boundaries */

/*
 * I_FLAG CONSTANTS
 * ------------------------------------------------------------------------
 *
 * These constants are used in the i[] field in the Problem_Description
 * structure
 *
 * Interpolation functions for variables...
 */

#define I_NOTHING 0  /* Nothing at all. No variable, even.*/
#define I_G0      1  /* Globally zero. */
#define I_G1      2  /* Globally constant. */
#define I_P0      3  /* Piecewise constant on ea element. */
#define I_P1      4  /* Piecewise linear on ea element. */
#define I_Q1      5  /* Lagrangian linear. */
#define I_Q2      6  /* Lagrangian quadratic. */
#define I_H3      7  /* Hermite cubic polynomials. */
#define I_S2      8  /* Serendipity quadratics (2D, 3D). */
#define I_B3      9  /* Cubic splines. */
#define I_Q3      10 /* Lagrangian cubic. */
#define I_Q4      11 /* Lagrangian quartic. */
#define I_SP                                    \
  12 /* Subparametric - linear on the interior, \
                quadratic on  surfaces */
#define I_Q2_D                            \
  13 /* Lagrangian quadratic with special \
                surface dofs  */
#define I_Q1_D                                          \
  14                  /* Lagrangian linear with special \
                                 surface dofs  */
#define I_PQ1      15 /* Bilinear discontinuous */
#define I_PQ2      16 /* Biquadratic discontinous */
#define I_Q2_LSA   17 /* Special Q2 interpolation for 3D of 2D LSA */
#define I_Q2_D_LSA 18 /* Special Q2_D interpolation for 3D of 2D LSA */
#define I_P0_G     19 /* Piecewise constant with ghost values. */
#define I_P1_G     20 /* Piecewise linear with ghost values. */
#define I_Q1_G     21 /* Lagrangian linear with ghost values. */
#define I_Q2_G     22 /* Lagrangian quadratic with ghost values. */
#define I_P0_XV    23 /* Piecewise constant with enrichment for jump in value. */
#define I_P1_XV    24 /* Piecewise linear with enrichment for jump in value. */
#define I_Q1_XV    25 /* Lagrangian linear with enrichment for jump in value. */
#define I_Q2_XV    26 /* Lagrangian quadratic with enrichment for jump in value. */
#define I_P1_XG    27 /* Piecewise linear with enrichment for jump in gradient. */
#define I_Q1_XG    28 /* Lagrangian linear with enrichment for jump in gradient. */
#define I_Q2_XG    29 /* Lagrangian quadratic with enrichment for jump in gradient. */
#define I_P0_GP    30 /* Piecewise constant confined to positive side of interface. */
#define I_P1_GP    31 /* Piecewise linear confined to positive side of interfaces. */
#define I_Q1_GP    32 /* Lagrangian linear confined to positive side of interface. */
#define I_Q2_GP    33 /* Lagrangian quadratic confined to positive side of interface. */
#define I_P0_GN    34 /* Piecewise constant confined to negative side of interface. */
#define I_P1_GN    35 /* Piecewise linear confined to negative side of interfaces. */
#define I_Q1_GN    36 /* Lagrangian linear confined to negative side of interface. */
#define I_Q2_GN    37 /* Lagrangian quadratic confined to negative side of interface. */
#define I_Q1_HG    38 /* Lagrangian linear with discontinuous enrichment for jump in gradient. */
#define I_Q1_HV    39 /* Lagrangian linear with discontinuous enrichment for jump in value. */
#define I_Q1_HVG \
  40 /* Lagrangian linear with discontinuous enrichment for jump in value and gradient. */
#define I_Q2_HG 41 /* Lagrangian quadratic with discontinuous enrichment for jump in gradient. */
#define I_Q2_HV 42 /* Lagrangian quadratic with discontinuous enrichment for jump in value. */
#define I_Q2_HVG \
  43 /* Lagrangian quadratic with discontinuous enrichment for jump in value and gradient. */
#define I_TABLE          44 /* Table Interpolation	*/
#define I_N1             45 // Nedelec First Kind
#define MAX_INTERP_TYPES 46

/*
 * Time dependence flags...
 */

#define DT_NOTHING 0x00000000 /* Steady state, pure and simple. */
#define DT_FE      0x00000001 /* Forward Euler. */
#define DT_BE      0x00000002 /* Backward Euler. */
#define DT_CN      0x00000004 /* Crank Nicholson. */
#define DT_AB      0x00000008 /* Adams Bashforth. */
#define DT_AM      0x00000010 /* Adams Moulton. */

/*
 * Constitutive relation flags...
 */

/*
 * Heat flux...
 */
#define CR_HF_FOURIER_0 0x00000001 /* q=k.grad(T); k=constant */
#define CR_HF_FOURIER_1 0x00000002 /* q=k(T).grad(T) */
#define CR_HF_USER      0x00000004 /*  q=user_fcn(gradT)       */

/*
 * Mesh stress/strain relation (displacement flux)...
 */
#define CR_DF_HOOKS_0     0x00000001 /* Linear elastic solid */
#define CR_DF_NONLINEAR_0 0x00000002 /* Non-Linear elastic solid */

/*
 * Momentum flux...
 */
#define CR_MF_NEWTON_0 0x00000001 /* Newton's Law of viscosity */
#define CR_MF_NEWTON_1 0x00000002 /* mu=mu(T) */
#define CR_MF_NEWTON_2 0x00000004 /* mu=mu(I2) */

/*
 * Mass flux...
 */

#define CR_SF_FICKS_0 0x00000001 /* Fick's Law of binary diffusion */

/* n_w = D.grad(c_w) */

/*
 * Element topology
 */
#define BODY   0
#define FACE   1
#define CURVE  2
#define VERTEX 3

/*****************************************************************************/
#endif
