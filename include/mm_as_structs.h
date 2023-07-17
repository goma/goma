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
 *$Id: mm_as_structs.h,v 5.35 2010-04-07 22:27:00 prschun Exp $
 */

/*
 * These structures are meant to hold various useful quantities that are
 * needed when assembling the residual equations and Jacobian entries for
 * all the different permutations of different problems that we might be
 * interested in solving. The idea is to be able to import this information
 * down to a low level of subroutine where it will be useful.
 *
 * The structures are revealed here, and declared "extern" in most places
 * via the include file mm_as.h. Then, in mm_as_alloc() they are
 * actually defined
 */

#ifndef GOMA_MM_AS_STRUCTS_H
#define GOMA_MM_AS_STRUCTS_H

#include "el_elm.h"
#include "mm_elem_block_structs.h"
#include "mm_mp_const.h"
#include "rf_bc_const.h"
#include "rf_io_const.h"
#include "rf_vars_const.h"
#include "sl_util_structs.h"
#include "std.h"

#ifndef MNROT
#define MNROT                                       \
  30 /* maximum number of rotation vector           \
      * sensitivities.  This is basically the       \
      * largest number of nodes whose displacements \
      * can affect a vector in the surface.         \
      */
#endif

/* define VIM according to problem type */
/* #define VIM	3  */ /* Vector dimensions so that loops over */
                      /* vector products do this many components */
                      /* eg. v.gradT loops over 3 components for */
                      /* all problems. This is convenient for */
                      /* axisymmetric problems and works OK for */
                      /* Cartesian problems too if third */
                      /* component is set =0.*/

#ifndef GOMA_CK_NAME_DEF
#define GOMA_CK_NAME_DEF
typedef char CK_NAME[64]; /* Typedefs for common names used for naming
                             domains and species */
typedef char CK_NAME_STR[64];
#endif

/*____________________________________________________________________________*/

/*
 *   Element_Indices
 *
 *    These element indeces are specific to each element in mesh. They
 *    provide mapping information for degrees of freedoms of
 *    variable types defined on nodes owned by the element.
 *
 *
 *    A couple of key concepts concerning two types of element degrees of
 *    freedom:
 *
 *      lvdof or ldof -> For each equation type, this is an index into the
 *                       degree of freedom present at one of the nodes in
 *                       the element. The index spans all nodes in the local
 *                       element, and there
 *                       may be more than one dof of the same equation type
 *                       at a single node.
 *                       Note, the degree of freedom may or
 *                       may not be "active" in the element. Active means
 *                       that it participates in the element interpolation
 *                       of that variable type within that element.
 *                       The maximum size of an lvdof index can not exceed
 *                       the static size, MDE, which is compiled into
 *                       the code.
 *
 *
 *      ledof -> This is a new concept. This number is a unique number for
 *               each degree of freedom in an element, irrespective of the
 *               variable type.
 *               Right now, each variable type will have consequative
 *               ledof numbers. Therefore, one can obtain ledof by looping
 *               over variable types, and then over lvdof's, while keeping
 *               an incremental counter.
 *               The maximum size of ledof can not exceed MDE * MAX_PROB_VAR.
 *               However, in the future, we will attempt to calculate
 *               an MDE equivalent dynamically.
 *
 *  Implementation Plan
 * -------------------------
 *  I've been thinking that the trying to specify degrees of freedom by their
 *  variable type ID alone is the wrong way to go. Instead, one should convert
 *  or replace most of the structures below to use uniqueness of the variable
 *  type structures as the base indexing methodology. This means in effect
 *  that the an index will be created that is unique for each variable type
 *  and matID (and species number) present in the element or specified on the
 *  sides of the element.
 *
 *   Local Variable Description Number:
 *         This is the unique index for a variable description in an element.
 *   Local Variable Description Dof Index:
 *         For each  Local Variable Description Number, this index rasters over
 *         all of the dof's in the element that belong to the pertaining
 *         variable description.
 *
 */

/* Define values for determining side of LS */

#define LS_OFF      0
#define LS_POSITIVE 1
#define LS_NEGATIVE 2

/*
 * Define for right now the maximum number of local variable type descriptions
 * to be equal to the maximum number of different variable types in the problem
 */
#define MAX_LOCAL_VAR_DESC (MAX_PROB_VAR + MAX_CONC)

#define LEC_R_INDEX(peqn_macro, index_macro) ((lec->max_dof * (peqn_macro)) + index_macro)

#define LEC_J_INDEX(peqn_macro, pvar_macro, index_i, index_j)                  \
  (((MAX_PROB_VAR + MAX_CONC) * lec->max_dof * lec->max_dof) * (peqn_macro)) + \
      ((lec->max_dof * lec->max_dof) * (pvar_macro)) + (lec->max_dof * (index_i)) + index_j

#define LEC_J_STRESS_INDEX(peqn, pvar, index_i, index_j)          \
  ((lec->max_dof * MAX_LOCAL_VAR_DESC * lec->max_dof) * (peqn)) + \
      ((MAX_LOCAL_VAR_DESC * lec->max_dof) * (pvar)) + (lec->max_dof * (index_i)) + index_j

#ifndef MAX_PHASE_FUNC
#define MAX_PHASE_FUNC 5
#endif

struct Element_Indices {
  int iconnect_ptr;                /* To find way back to global scheme . This is an
                                    * index into the element connectivity list for
                                    * the current element:
                                    *
                                    *   e.g., Proc_Elem_Connect[iconnect_ptr + i]
                                    */
  int ielem;                       /*
                                    *  The id of the current element
                                    */
  int ielem_dim;                   /* is this a 2d or 3d problem */
  int ielem_type;                  /* what kind of element */
  int ielem_shape;                 /* element shape, fundamental topology (new) */
  int elem_blk_index;              /* element block index of block containing
                                    * current element */
  int elem_blk_id;                 /* Element block ID of block containing
                                    * current element */
  ELEM_BLK_STRUCT *current_EB_ptr; /* Pointer to the current element block */
  int mn;                          /* material number corresponding to the
                                    * current element */
  int num_local_nodes;             /* how many local nodes in element */
  int num_sides;                   /* how many sides does the current element
                                      have (at least for this processor) */
  int deforming_mesh;              /* Can nodes move?  NOTE: For shell
                                      elements, this may be TRUE even though
                                      mesh equations are not defined on the
                                      shell block! */

  /*
   * A note about the order of variable descriptions:
   *
   *       They are ordered wrt Var Type first
   *          matID next -> current matID is always first however
   *                        and after that matID order is unspecified.
   *             species number or idof next. (these are grouped together).
   */

  int Num_Lvdesc;               /* Number of local variable descriptions in the
                                 * element, whether they are interpolated in
                                 * the element or just exist on the sides
                                 * of the elements.
                                 */
  int *Num_Lvdesc_Per_Var_Type; /* Number of local variable descriptions per
                                 * variable types in the current element
                                 *  MallocSize = MAX_VARIABLE_TYPES + MAX_CONC
                                 *  actualSize = MAX_VARIABLE_TYPES +
                                 * Max_Num_Species
                                 */
  int *Lvdesc_First_Var_Type;   /* Mapping of the variable type to the
                                 * first local variable description number
                                 * pertaining to that variable type.
                                 * The first variable type
                                 * will in all possible cases be the variable
                                 * description that is an active volumetric
                                 * interpolation.
                                 *
                                 *  MallocSize = MAX_VARIABLE_TYPES + MAX_CONC
                                 *  actualSize = MAX_VARIABLE_TYPES +
                                 * Max_Num_Species
                                 */
  int **Lvdesc_to_ledof;        /* Mapping between the (local variable description
                                 * number, local variable description
                                 * dof index) pair to the local element dof index,
                                 * ledof.
                                 *   MallocSize = [MAX_LOCAL_VAR_DESC][MDE]
                                 *   ActualSize = [Num_lvdesc][Lvdesc_Numdof[i_vdesc]]
                                 */
  int **Lvdesc_to_lvdof;        /* Mapping between the (local variable description
                                 * number, local variable description
                                 * dof index) pair to the local variable dof number,
                                 * lvdof, of the same variable type. For MF's with a
                                 * subvariable types greater than one, the mapping
                                 * is made to the subvariable type=0 lvdof.
                                 *   MallocSize = [MAX_LOCAL_VAR_DESC][MDE]
                                 *   ActualSize = [Num_lvdesc][Lvdesc_Numdof[i_vdesc]]
                                 */
  int *Lvdesc_to_Var_Type;      /* Mapping between the local variable description
                                 * number to the variable type.
                                 *   MallocSize = [MAX_LOCAL_VAR_DESC]
                                 *   ActualSize = [Num_Lvdesc]
                                 */
  int *Lvdesc_to_MatID;         /* Mapping between the local variable description
                                 * number to the matID.
                                 *   MallocSize = [MAX_LOCAL_VAR_DESC]
                                 *   ActualSize = [Num_Lvdesc]
                                 */
  int *Lvdesc_Numdof;           /* Number of degrees of freedom for each variable
                                 * description number in the current element.
                                 *   MallocSize = [MAX_LOCAL_VAR_DESC]
                                 *   ActualSize  = [Num_Lvdesc]
                                 */
  int **Lvdesc_to_Gnn;          /* Mapping between the (local variable description
                                 * number, local variable description
                                 * dof index) pair to the processor's node
                                 * number.
                                 *   MallocSize = [MAX_LOCAL_VAR_DESC][MDE]
                                 *   ActualSize = [Num_Lvdesc][Lvdesc_Numdof[i_vdesc]]
                                 */
  int **Lvdesc_to_Gun;          /* Mapping between the (local variable description
                                 * number, local variable description
                                 * dof index) pair to the processor's global unknown
                                 * number.
                                 *   MallocSize = [MAX_LOCAL_VAR_DESC][MDE]
                                 *   ActualSize = [Num_Lvdesc][Lvdesc_Numdof[i_vdesc]]
                                 */
  int *Lvdesc_to_MFSubvar;      /*  For variable types with subvar indeces
                                 *  ie Mass Fraction ,
                                 *  we need another variable to describe the
                                 *  subvariable index (i.e., the species
                                 *  number) corresponding to the current
                                 *  variable description.
                                 *   MallocSize = [MAX_LOCAL_VAR_DESC]
                                 *   ActualSize  = [Num_Lvdesc]
                                 */
  int **Lvdesc_to_Lnn;          /* Mapping between the (local variable description
                                 * number, local variable description
                                 * dof index) pair to the local node number
                                 * within the element.
                                 *   MallocSize = [MAX_LOCAL_VAR_DESC][MDE]
                                 *   ActualSize = [Num_Lvdesc][Lvdesc_NumDof[i_vdesc]]
                                 */
  int **Lvdesc_Lnn_to_Offset;   /* Mapping between the (local variable description
                                 * number, local node number) pair to the offset
                                 * from the beginning of the solution vector at the
                                 * current node.
                                 *   MallocSize = [MAX_LOCAL_VAR_DESC][MDE]
                                 *   ActualSize = [Num_Lvdesc][num_local_nodes]
                                 */
  int **Lvdesc_Lnn_to_lvdof;    /* Mapping between the (local variable description
                                 * number, local node number) pair to the local
                                 * variable dof index.  For MF's with a
                                 * subvariable types greater than one, the mapping
                                 * is made to the subvariable type=0 lvdof.
                                 *   MallocSize = [MAX_LOCAL_VAR_DESC][MDE]
                                 *   ActualSize = [Num_Lvdesc][num_local_nodes]
                                 */
  int **lvdof_to_lvdesc;        /* Mapping between the local variable index and
                                 * a variable_type pair to the corresponding
                                 * local variable description index.
                                 *   MallocSize
                                 * = [MAX_VARIABLE_TYPES+MAX_CONC][MDE]
                                 *   ActualSize
                                 *     = [MAX_VARIABLE_TYPES+Max_Species]
                                 *       [ei->dof[v]]
                                 */
  int **lvdof_to_lvdesc_dof;    /*  Mapping between the local variable index and
                                 *  a variable_type pair to the corresponding
                                 *  dof number of the local variable description
                                 *  index.
                                 *   MallocSize
                                 * = [MAX_VARIABLE_TYPES+MAX_CONC][MDE]
                                 *   ActualSize
                                 *     = [MAX_VARIABLE_TYPES+Max_Species]
                                 *       [ei->dof[v]]
                                 */
  int **Lvdesc_Lnn_Numdof;      /*  Number of degrees of freedom for the variable
                                 *  description number at the local node number
                                 *   MallocSize = [MAX_LOCAL_VAR_DESC][MDE]
                                 *   ActualSize = [Num_Lvdesc][num_local_nodes]
                                 */
  VARIABLE_DESCRIPTION_STRUCT **Lvdesc_vd_ptr;
  /*  Pointer to the variable description
   *  structure corresponding to the lvdesc index.
   *   MallocSize = [MAX_LOCAL_VAR_DESC]
   *   ActualSize  = [Num_Lvdesc]
   */
  int *VDindex_to_Lvdesc; /*  Index pointing from the variable description
                           *  index to the current local variable description
                           *  index for the local variable
                           */

  int *owningElementForColVar; /* Indicator of who owns the column variable
                                * type. We separate and distinguish all columns
                                * in the local jacobian, according to whether
                                * they refer to equations that belong to the
                                * parent, the child, or the owning element.
                                *   [MAX_VARIABLE_TYPES+MAX_CONC]
                                *
                                *     ei->owningElementForColVar[varType] = elem
                                *
                                * where value is one of these:
                                *
                                *    -1     not this element
                                *    elem   This element. This means that the
                                * variable has an equation on this element, and
                                *           that there is an interpolation of
                                * the variable using this element's basis
                                * functions. elemP  Parent or child of this
                                * element. This means that the variable has an
                                * equation on this element's Parent, and that
                                * there is an interpolation of the variable
                                * using this element Parent's basis functions.
                                */

  struct Element_Indices *owningElement_ei_ptr[MAX_VARIABLE_TYPES];
  /* Pointer to the Element_Indices struct that contains
   * the info for the element. The element number is the
   * listed in owningElementForColVar[].
   */

  /*******************************************************/
  /* Slated to be AXED: */

  int *dof;                 /* Number of dof's in the element's nodes for each
                             * variable type (both, active and inactive)
                             *   Length = MAX_VARIABLE_TYPES */
  int *dof_ext;             /* Number of dof's in the elem for each external
                             * fixed variable
                             *   Length = MAX_EXTERNAL_TYPES */
  int **Baby_Dolphin;       /* Offset for the local variable in the
                             * solution vector in cases where there are
                             * discontinuous variables and multiple dofs
                             * of the same variable type at the same node.
                             * Note: doesn't take into account multiple
                             *       subvariable types for mass fractions.
                             *   Size = [MAX_VARIABLE_TYPES][MDE] */
  int **dof_list;           /* list of local nodes for ea dof in the elem
                             *      Size [MAX_VARIABLE_TYPES][MDE] */
  int **gnn_list;           /* list of global node numbers for ea dof
                             *      Size [MAX_VARIABLE_TYPES][MDE] */
  int **gun_list;           /* list of global unknown numbers for ea dof
                             *      Size [MAX_VARIABLE_TYPES][MDE] */
  int **ln_to_dof;          /* Mapping between the (variable type, local
                             * element node) pair to the local variable
                             * degree of freedom in the element pertaining
                             * to the active interpolation of that
                             * variable within the element
                             *   MallocSize
                             *     = [MAX_VARIABLE_TYPES][MDE]
                             *   ActualSize
                             *     = [MAX_VARIABLE_TYPES][num_local_nodes] */
  int **ln_to_first_dof;    /* list of local dof numbers for the first
                             * dof at ea local node
                             *      Size [MAX_VARIABLE_TYPES][MDE] */
  int *active_interp_ledof; /* This variable is equal to one if this
                             * ledof is part of an active interpolation
                             * of the variable within the element.
                             * It is equal to zero if it is not,
                             * (if for example, it is located on the
                             *  side of an element, and is part of another
                             *  material domain).
                             *      Size [MDE * MAX_PROB_VAR] */
  int **lvdof_to_row_lvdof; /* Identification of the row that an lvdof
                             * volumetric contribution gets assigned to.
                             * Usually, it is a one-to-one mapping. However,
                             * some boundary conditions require volumetric
                             * contributions from one material be assigned
                             * to unknowns from another material. In this
                             * case and this case alone, this field
                             * diverts from a one-to-one mapping. Note,
                             * this field need only contain valid information
                             * for active degrees of freedom in the
                             * element.
                             * (This field replaces the old field,
                             *  first_active_dof).
                             * size [MAX_VARIABLE_TYPES][MDE] */
  int **lvdof_to_ledof;     /*  Mapping between the local variable dof
                             *  number
                             *  the local element dof number, ledof
                             *      Size [MAX_VARIABLE_TYPES][MDE] */

  /*******************************************************/

  int *owned_ledof;                               /*  This is a boolean variable indicating
                                                   *  whether
                                                   *  this local element dof is owned by the
                                                   *   processor or not.
                                                   *      Size [MDE * MAX_PROB_VAR] */
  int *ieqn_ledof;                                /*  Processor equation number for the local
                                                   *  element  dof. Note, a dof will have
                                                   *  an equation number irrespective of whether
                                                   *  or not it is owned by the processor.
                                                   *      Size [MDE * MAX_PROB_VAR] */
  int *matID_ledof;                               /*  Material ID corresponding to the local
                                                   *  element dof. Here, the generic material
                                                   *  ID, -1, will not be used. Instead, the
                                                   *  element material ID will be used in its
                                                   *  place. Thus, the only time a ledof will
                                                   *  not have a matID equal to the element
                                                   *  matID is if it is inactive in the element
                                                   *  and corresponds to an active degree of
                                                   *  freedom
                                                   *  in another element that has a different
                                                   *  material ID.
                                                   *      Size [MDE * MAX_PROB_VAR] */
  int **MFsubvar_Offset;                          /*  For variable types with subvar indeces,
                                                   *  ie Mass Fraction,
                                                   *  we need another variable to describe the
                                                   *  offset of each of the subvariables wrt
                                                   *  the base subvar=0 index for the current
                                                   *  local variable dof corresponding to the
                                                   *  mass fraction variable type.
                                                   *      Size [MAX_CONC] [MDE]
                                                   */
  int linkedEIelems[MAX_ELEMENT_INDICES_RELATED]; /* List of elements ids that
                                                   * are slaves to this element
                                                   */
};

/*____________________________________________________________________________*/

struct Element_Variable_Pointers {
  /*
   * These are an abbreviated version of the Element_Stiffness_Pointers
   * containing only the nodal variable.
   *
   * this structure was designed to hold the old nodal variable and dot
   * variables in a convenient way
      *
   * Pointers to *unknowns*....(i.e, nodal point values for variables) that
   * are important in this element.  -----
   *
   * For variables defined at Gauss points, see the field variable sections...
   *

   * Idea --
   *		*T_dot[i] = x_dot[something];
   *
   *		   i      == local degree of freedom counter for temperature
   *		something == equivalent position in global unknown vector
   *
   * Usage:	*esp_dot->T[i] points to right place in x_dot[] to get the value
   you *		need.
   */

  dbl *T[MDE];                                  /* temperature */
  dbl *v[DIM][MDE];                             /* velocity */
  dbl *v_star[DIM][MDE];                        /* velocity* segregated */
  dbl *d[DIM][MDE];                             /* mesh displacement */
  dbl *d_rs[DIM][MDE];                          /* real solid displacement */
  dbl *c[MAX_CONC][MDE];                        /* concentration     */
  dbl *external_field[MAX_EXTERNAL_FIELD][MDE]; /* external field variables at nodes */
  dbl *initial_displacements[2 * DIM][MDE];     /* initial xyz displacement fields
                                                   for real solid and pseudo solid
                                                   for annealing */
  dbl *P[MDE];                                  /* pressure */
  dbl *P_star[MDE];                             /* pressure */
  dbl *S[MAX_MODES][DIM][DIM][MDE];             /* polymeric stress tensor, for each mode */
  dbl *G[DIM][DIM][MDE];                        /* velocity gradient tensor */
  dbl *F[MDE];                                  /* Fill */
  dbl *V[MDE];                                  /* Potential; added by KSC: 2/3/99 */
  dbl *qs[MDE];                                 /* Surface charge density */
  dbl *Enorm[MDE];                              /* |E| for dielectrophoresis. */
  dbl *pv[DIM][MDE];                            /* Particle velocity */

  dbl *p_liq[MDE];    /* *p_liq[MDE], liquid-phase pressure in porous media */
  dbl *p_gas[MDE];    /* *p_gas[MDE], liquid-phase pressure in porous media */
  dbl *porosity[MDE]; /* *porosity[MDE],liquid-phase pressure porous media */

  dbl *vd[DIM][MDE]; /* Vorticity principle flow direction (NOT
                      *  vorticity vector) */
  dbl *vlambda[MDE]; /* Eigenvalue associated with vorticity direction */
  dbl *nn[MDE];      /* This is the bond evolution */

  dbl *lm[DIM][MDE]; /* Lagrange Multiplier field */

  dbl *ext_v[MDE]; /* Extension velocity, normal direction */

  dbl *E_field[DIM][MDE]; /* Electric field */

  dbl *H[DIM]; /* Level Set Curvature */

  dbl *n[DIM][MDE]; /* level set normal OR shell normal */

  dbl *sh_K[MDE];            /* Shell curvature */
  dbl *sh_K2[MDE];           /* Shell second curvature */
  dbl *sh_tens[MDE];         /* Shell tension */
  dbl *sh_x[MDE];            /* Shell y coordinate */
  dbl *sh_y[MDE];            /* Shell x coordinate */
  dbl *sh_u[MDE];            /* Shell user */
  dbl *sh_ang[DIM - 1][MDE]; /* Shell orientation angles */
  dbl *div_s_v[MDE];         /* sundry pieces (next 4) for surface rheological
                                constitutive eqn */
  dbl *curv[MDE];
  dbl *grad_v_dot_n[DIM][MDE];   /* grad_s (v_) dot n  Vector variable  */
  dbl *n_dot_curl_s_v[MDE];      /* n dot (curl_s v) Scalar variable used in shell
                                    equations - curl_s is surface curl */
  dbl *pF[MAX_PHASE_FUNC][MDE];  /* phase function */
  dbl *sh_J[MDE];                /* sh_J[MDE], Shell surface diffusion flux */
  dbl *sh_Kd[MDE];               /* sh_Kd[MDE], Shell surface curvature */
  dbl *apr[MDE];                 /* acoustic pressure real part */
  dbl *api[MDE];                 /* acoustic pressure imag part */
  dbl *epr[MDE];                 /* em lagr pressure real part */
  dbl *epi[MDE];                 /* em lagr pressure imag part */
  dbl *ars[MDE];                 /* acoustic reynolds stress */
  dbl *sink_mass[MDE];           /* porous sink mass*/
  dbl *sh_bv[MDE];               /* acoustic boundary velocity */
  dbl *sh_p[MDE];                /* shell lub pressure */
  dbl *lubp[MDE];                /* lub pressure */
  dbl *lubp_2[MDE];              /* second lub pressure */
  dbl *sh_fp[MDE];               /* lub pressure in thin film */
  dbl *sh_fh[MDE];               /* film thickness */
  dbl *sh_pc[MDE];               /* particles concentration */
  dbl *sh_sat_closed[MDE];       /*  Porous shell saturation - Closed cells - SAR */
  dbl *sh_p_open[MDE];           /*  Porous shell pressure - Open closed - SAR */
  dbl *sh_p_open_2[MDE];         /*  Second Porous shell pressure - Open closed -PRS */
  dbl *sh_t[MDE];                /*  Shell temperature -- PRS*/
  dbl *sh_dh[MDE];               /*  Shell delta gap   -- PRS*/
  dbl *sh_l_curv[MDE];           /* Lubrication shell curvature - SAR */
  dbl *sh_l_curv_2[MDE];         /* Lubrication shell curvature 2 - PRS */
  dbl *sh_sat_gasn[MDE];         /*  Porous shell saturation - Gas compression - SAR */
  dbl *sh_shear_top[MDE];        /* Top wall shear rate */
  dbl *sh_shear_bot[MDE];        /* Bottom wall shear rate */
  dbl *sh_cross_shear[MDE];      /* Cross stream shear stress */
  dbl *max_strain[MDE];          /* Maximum Von Mises strain */
  dbl *cur_strain[MDE];          /* Von Mises strain */
  dbl *poynt[DIM][MDE];          /* Poynting Vector for light intensity */
  dbl *tfmp_pres[MDE];           /* thin-film multi-phase lubrication pressure */
  dbl *tfmp_sat[MDE];            /* thin-film multi-phase saturation */
  dbl *moment[MAX_MOMENTS][MDE]; /* moments */
  dbl *rho[MDE];
  dbl *restime[MDE];    /* residence time function field */
  dbl *em_er[DIM][MDE]; /* EMwave Electric Field real part */
  dbl *em_ei[DIM][MDE]; /* EMwave Electric Field imag part */
  dbl *em_hr[DIM][MDE]; /* EMwave Magnetic Field real part */
  dbl *em_hi[DIM][MDE]; /* EMwave Magnetic Field imag part */
  dbl *sh_sat_1[MDE];   /* Porous shell saturation layer 1 */
  dbl *sh_sat_2[MDE];   /* Porous shell saturation layer 2 */
  dbl *sh_sat_3[MDE];   /* Porous shell saturation layer 3 */

  dbl *eddy_nu[MDE]; /* Eddy viscosity for turbulent flow */
};

/*___________________________________________________________________________*/

struct Element_Stiffness_Pointers {
  /*
   * These are basically BULK degrees of freedom here. If you have need to
   * account for boundary sensitivities, then you will need to think
   * through the best procedure for doing so...
   *
   * Pointers to *unknowns*....(i.e, nodal point values for variables) that
   * are important in this element.  -----
   *
   * For variables defined at Gauss points, see the field variable sections...
   *
   * Idea --
   *		*T[i] = x[something];
   *
   *		   i      == local degree of freedom counter for temperature
   *		something == equivalent position in global unknown vector
   *
   * Usage:	*esp->T[i] points to right place in x[] to get the value you
   *		need.
   *
   * space for these arrays is allocated in mm_as_alloc.c if the variable
   * are defined
   */

  dbl **T;       /* *T[MDE], temperature */
  dbl ***v;      /* *v[DIM][MDE], velocity */
  dbl ***v_star; /* *v_star[DIM][MDE], velocity* segregated */
  dbl ***d;      /* *d[DIM][MDE], mesh displacement */
  dbl ***d_rs;   /* *d_rs[DIM][MDE], real solid displacement */
  dbl ***c;      /* *c[MAX_CONC][MDE], concentration  */
  dbl **P;       /* *P[MDE], pressure */
  dbl **P_star;
  dbl *****S; /* *S[MAX_MODES][DIM][DIM][MDE], polymeric
                 stress tensor */
  dbl ****G;  /* *G[DIM][DIM][MDE], velocity gradient
                 tensor */
  dbl **V;    /* *V[MDE], voltage potential */
  dbl **qs;   /* *qs[MDE], surface charge density */
  dbl **F;    /* *F[MDE], fill */
  dbl **SH;   /* *SH[MDE], shear rate from second
                   invariant */

  dbl **Enorm; /* Enorm[MDE], |E| */
  dbl **H;     /* H[MDE] curvature of level set function */

  dbl ***pv; /* *pv[DIM][MDE], particle velocity */

  dbl **p_liq;    /* *p_liq[MDE], liquid-phase pressure in porous media */
  dbl **p_gas;    /* *p_gas[MDE], liquid-phase pressure in porous media */
  dbl **porosity; /* *porosity[MDE],liquid-phase pressure porous media */
  dbl ***vd;      /* *vd[DIM][MDE}, vorticity prinicple flow direction. */
  dbl **vlambda;  /* *vlambda[MDE], eigenvalue associated with vd. */
  dbl **nn;       /* *nn[MDE], bond evolution */

  dbl ***lm; /* lm_gama[DIM][MDE] is the lagrange multiplier vector field */

  dbl **ext_v; /* Extension velocity */

  dbl ***E_field; /* Electric field */

  dbl ***n; /* n[DIM][MDE],  level set normal OR shell normal */

  dbl **sh_K;    /* sh_K[MDE],   Shell curvature */
  dbl **sh_K2;   /* sh_K2[MDE],  Shell second curvature */
  dbl **sh_tens; /* sh_tens[MDE], Shell tensions */
  dbl **sh_x;    /* sh_x[MDE],   Shell x coordinate */
  dbl **sh_y;    /* sh_y[MDE], Shell y coordinate */
  dbl **sh_u;    /* sh_u[MDE], Shell user */
  dbl ***sh_ang; /* sh_ang[DIM-1][MDE], Shell orientation angles */
  dbl **div_s_v; /* sundry pieces (next 4) for surface rheological constitutive
                    eqn */
  dbl **curv;
  dbl ***grad_v_dot_n;  /* grad_s_v_dot_n[DIM][MDE] , grad_s_v_dot_n[] defined at
                           a surface */
  dbl **n_dot_curl_s_v; /* n dot (curl_s v) Scalar variable used in shell
                           equations - curl_s is surface curl */
  dbl ***pF;            /* *pF[MAX_PHASE_FUNC][MDE], phase function */
  dbl **sh_J;           /* sh_J[MDE], Shell surface diffusion flux */
  dbl **sh_Kd;          /* sh_Kd[MDE], Shell surface curvature */
  dbl **apr;            /* *apr[MDE], acoustic pressure */
  dbl **api;            /* *api[MDE], acoustic pressure */
  dbl **ars;            /* *ars[MDE], acoustic reynolds stress */
  dbl **epr;            /* *epr[MDE], em pressure */
  dbl **epi;            /* *epi[MDE], em pressure */
  dbl **sink_mass;      /* Porous sink mass */
  dbl **sh_bv;          /* sh_bv[MDE], acoustic bdy velocity */
  dbl **sh_p;           /* sh_p[MDE], lub pressure */
  dbl **lubp;           /* lubp[MDE], lub pressure */
  dbl **lubp_2;         /* lubp_2[MDE], second lub pressure */
  dbl **sh_fp;          /* sh_fp[MDE], lub pressure in the thin film */
  dbl **sh_fh;          /* sh_fh[MDE], film thickness */
  dbl **sh_pc;          /* sh_pc[MDE], particles concentration */
  dbl **sh_sat_closed;  /* sh_sat_closed[MDE], porous shell saturation - closed
                           cells - SAR */
  dbl **sh_p_open;      /* sh_p_open[MDE], porous shell pressure - open cells - SAR */
  dbl **sh_p_open_2;    /* sh_p_open_2[MDE], porous shell pressure - open cells -
                           PRS */
  dbl **sh_t;           /* sh_t[MDE], Shell temperature -- PRS */
  dbl **sh_dh;          /* sh_dh[MDE], Shell delta_h -- PRS */
  dbl **sh_l_curv;      /* sh_l_curv[MDE], Lubrication shell curvature - SAR */
  dbl **sh_l_curv_2;    /* sh_l_curv_2[MDE], Lubrication_2 shell curvature - PRS */
  dbl **sh_sat_gasn;    /* sh_sat_gasn[MDE], porous shell saturation - gas
                           compression - SAR */
  dbl **sh_shear_top;   /* sh_shear_top[MDE], top wall shear rate */
  dbl **sh_shear_bot;   /* sh_shear_bot[MDE], bottom wall shear rate */
  dbl **sh_cross_shear; /* sh_cross_shear[MDE], cross stream shear stress */
  dbl **max_strain;     /* max_strain[MDE], maximum Von Mises strain */
  dbl **cur_strain;     /* cur_strain[MDE], Von Mises strain */
  dbl ***poynt;         /* *v[DIM][MDE], velocity */
  dbl **tfmp_pres;      /*  thin-film multi-phase lubrication pressure */
  dbl **tfmp_sat;       /* thin-film multi-phase saturation */
  dbl ***moment;        /* *moment[MAX_MOMENTS][MDE], moments */
  dbl **rho;
  dbl **restime;  /* Residence Time Function Field */
  dbl ***em_er;   /* *em_xx[DIM][MDE], em_wave*/
  dbl ***em_ei;   /* *em_xx[DIM][MDE], em_wave*/
  dbl ***em_hr;   /* *em_xx[DIM][MDE], em_wave*/
  dbl ***em_hi;   /* *em_xx[DIM][MDE], em_wave*/
  dbl **sh_sat_1; /* Porous shell saturation layer 1 */
  dbl **sh_sat_2; /* Porous shell saturation layer 2 */
  dbl **sh_sat_3; /* Porous shell saturation layer 3 */

  dbl **eddy_nu; /* Eddy viscosity for turbulent flow */

  /*
   * These are for debugging purposes...
   */

  dbl *a0; /* point to beginning of global Jacobian */
};

/*___________________________________________________________________________*/

/*
 * Information for calculation of element quality (distortion) metrics
 * and stop/remesh criterion based on quality.
 */

struct Element_Quality_Metrics {
  int do_jac;     /* Jacobian metric flag   (if requested) */
  int do_vol;     /* Volume change metric flag (if requested) */
  int do_ang;     /* Angle    metric flag   (if requested) */
  int do_tri;     /* Triangle metric flag   (if requested) */
  double wt_jac;  /* Jacobian metric weight (if requested) */
  double wt_vol;  /* Volume change metric weight (if requested) */
  double wt_ang;  /* Angle    metric weight (if requested) */
  double wt_tri;  /* Triangle metric weight (if requested) */
  double eq_jac;  /* Jacobian metric value		 */
  double eq_vol;  /* Volume change metric value		 */
  double eq_ang;  /* Angle    metric value 		 */
  double eq_tri;  /* Triangle metric value		 */
  double eq_avg;  /* Weighted metric average		 */
  double eq_low;  /* Limiting (lowest) quality metric	 */
  double eq_tol;  /* Stop criterion for quality		 */
  double vol_sum; /* Volume change global sum		 */
  double vol_low; /* Volume change global minimum		 */
  int vol_count;  /* Volume change Gauss point counter	 */
  int tol_type;   /* Tolerance type indicator		 */
};

/*___________________________________________________________________________*/

/*
 * Where the "esp" were just pointers into x[], resid_vector[], and a[], here
 * we load up actual contributions.
 */

struct Local_Element_Contributions {
  int max_dof;
  dbl *R;
  dbl *J;
  /* For face m and  mode k we have for mode imode
     d(tau_12_i)/d(tau_12_j) =
       J_stress_neighbor[m][i][POLYMER_STRESS11_k][j]
  */
  dbl *J_stress_neighbor;

  /*
   * NOTE: concentration entries in local element arrays are stored at
   *       the end of
   *       the equation and variable lists (i.e. for species w use l
   *       ec->R[MAX_VARIABLE_TYPES + w][i]
   *       to get entry in local residual array
   */

  /*
   * Unused for now; this could contain handy local copies of the
   * global unknowns...
   *
   *  dbl x[MAX_VARIABLE_TYPES] [MDE];
   */
};

/*___________________________________________________________________________*/

/*
 * These might be useful for different calls to the "do-all" assembly routines
 * to just do residual calculations (eg, for numerical Jacobian checking).
 * Later, you could add parts to do new columns on the RHS or for the Jacobian
 * matrix.
 */

struct Action_Flags {
  int Assemble_Residual;
  int Assemble_Jacobian;
  int Assemble_LSA_Jacobian_Matrix; /* Whether or not to compute the
                                     * Jacobian (J) matrix for the
                                     * generalized eigenvalue problem for
                                     * linear stability analysis, J x =
                                     * \lambda B x
                                     */
  int Assemble_LSA_Mass_Matrix;     /* Whether or not to compute the
                                     * "mass" (B) matrix for the
                                     * generalized eigenvalue problem for
                                     * linear stability analysis,
                                     * J x = \lambda B x
                                     */
  int Sat_hyst_reevaluate;          /* This placeholder is used to initiate
                                     * a re-evaluation of the hysteresis
                                     * saturation curve parameters based on
                                     * some chosen criteria.   The idea is
                                     * to use this to control the evaluation,
                                     * which can either be by Newton iteration
                                     * or by time step. See rf_solve.c for it's
                                     * initialization and load_saturation for it's
                                     * use.
                                     */
  /*
   * Unused for now.
   *  int Assemble_2nd_RHS;
   */
};

typedef struct turbulent_information {
  double *wall_distances;
  int *side_set_ids;
  int *node_set_ids;
  int num_side_sets;
  int num_node_sets;
  int use_internal_wall_distance;
} turbulent_information;

/*
 * This contains information that is uniformaly relevant
 * to all portions of the problem without regard to
 * block id or material number
 *
 */

struct Uniform_Problem_Description {
  int Total_Num_Matrices; /* Total number of problem graphs to be solved */

  int Total_Num_EQ[MAX_NUM_MATRICES];  /* This is used in conjunction with the
                                        * ep[] array.  The number of nonzero
                                        * entries in this array will  equal the
                                        * number of non neg 1 entries in the
                                        * ep[][]  array.
                                        */
  int Total_Num_Var[MAX_NUM_MATRICES]; /* This is used in conjunction with the
                                        * vp[] array. The number of nonzero
                                        * entries in this array will equal the
                                        * number of non neg 1 entries in the
                                        * vp[][] array.
                                        */
  int CoordinateSystem;
  int vp[MAX_NUM_MATRICES]
        [MAX_VARIABLE_TYPES + MAX_CONC]; /* Mapping from the actual variable type index
                                          * to a uniform problem variable
                                          * index valid for all materials. If a variable type is
                                          * active anywhere in the domain, then its corresponding
                                          * entry in this array will be nonzero and contain a unique
                                          * index.
                                          */
  int ep[MAX_NUM_MATRICES]
        [MAX_EQNS + MAX_CONC]; /* Mapping from the actual equation variable type index
                                * to a uniform problem equation
                                * index valid for all materials. If a variable type is
                                * active anywhere in the domain, then its corresponding
                                * entry in this array will be nonzero and contain a unique
                                * index.
                                */
  int Max_Num_Species;         /* The maximum number of species in any one
                                  volumetric materials domain in the problem */
  int Max_Num_Species_Eqn;     /* The maximum number of species equations
                                  in any one  volumetric materials domain
                                  in the problem */
  int Tot_Num_VolSpecies;      /* Total number of different volumetric species in
                                  all of the domains */
  int Num_Mat;                 /* Total number of materials, eventually this will be
                                  distinct from the the number of element blocks
                                  It will be less than or equal to the number of
                                  element blocks. */
  int Species_Var_Type;        /* Default type of the species variable, i.e., mass
                                  fraction,    mole fraction,  concentration, capillary
                                  pressure,    etc employed for the independent variable.
                                  This may    be overwritten by the value in the materials
                                  structure    for the present material. The acceptable
                                  values are listed    in rf_fem_const.h. This variable
                                  influences all aspects    of the species conservation
                                  equation, as well as    the names that are put into the
                                  output files. It    determines the units for the species
                                  equation, for    example. */
  double Pressure_Datum;       /* Set the pressure Datum for use in thermodynamic
                                  equations of state calculations.
                                  This is an additive constant that get added onto the
                                  pressure field before calculation of all
                                  thermodynamic   equations of state calculations. It is
                                  a constant   over the entire domain. Therefore, it is
                                  not included   in any one materials file. The default
                                  units for   the quantity are cgs units, and the
                                  default value   for the quantity is 1 atmosphere
                                  = 1.01325E6   gm cm-1 sec-2  (dyne cm-2).
                                  (conversion factor is to an exact standard atm) */
  int Max_Num_Porous_Eqn;      /* max number of porous media Equations */
  dbl Acoustic_Frequency;      /* Frequency for Acoustic Harmonic Eqns */
  dbl EM_Frequency;            /* Frequency for Time-Harmonic Maxwell Eqns */
  dbl Free_Space_Permittivity; /* Free space permittivity for Time-Harmonic Maxwell Eqns */
  dbl Free_Space_Permeability; /* Free space permeability for Time-Harmonic Maxwell Eqns */
  dbl Light_Cosmu;             /* Inclination of Incident Light */
  dbl Process_Temperature;     /* Temperature for thermal property data */
                               /*   for isothermal problems */
  int XFEM;                    /* Flag indicating that XFEM is in use */
  int SegregatedSolve;         /* Flag indicating segregated solve should be used */
  int SegregatedSubcycles;
  int PSPG_advection_correction;
  int matrix_index[MAX_VARIABLE_TYPES];
  int petsc_solve_post_proc;
  void *petsc_post_proc_data;
  int devss_traceless_gradient;
  turbulent_information *turbulent_info;
};
typedef struct Uniform_Problem_Description UPD_STRUCT;
/*____________________________________________________________________________*/

/*
 * Problem_Graph Structure:
 *
 * Problem graph related structure containing information needed for segregated
 * solver
 *
 */
struct Problem_Graph {
  int imtrx; /* Current active matrix index */

  /* Temporarily make some things global */
  struct Matrix_Data *matrices;
  struct Matrix_Data *sub_step_solutions;
  int subcycle_fraction[MAX_NUM_MATRICES];
  double delta_t_fraction[MAX_NUM_MATRICES];

  int time_step_control_disabled[MAX_NUM_MATRICES];
  int matrix_subcycle_count[MAX_NUM_MATRICES];
  double sub_delta_t[MAX_NUM_MATRICES];
  double sub_delta_t_old[MAX_NUM_MATRICES];
  double sub_delta_t_older[MAX_NUM_MATRICES];
};
typedef struct Problem_Graph PROBLEM_GRAPH_STRUCT;
/*____________________________________________________________________________*/

/*
 * Problem_Description Structure:
 *
 * Values of equation and variable activity for the problem in the current
 * element block.
 */

struct Problem_Description {
  int Num_Matrices; /* Number of matrices in each element block */
  int Matrix_Activity[MAX_NUM_MATRICES];
  /* Matrix activity field in each element block
   *  0  -> Matrix is off
   *  1  -> Matrix is on
   */
  int Num_EQ[MAX_NUM_MATRICES];                /* number of active equations */
  int e[MAX_NUM_MATRICES][MAX_EQNS];           /* This is a vector containing  the
                                                * active equation terms for each equation.
                                                * within  the current element block.
                                                * Each bit in the integer refers to a
                                                * different term that is either on or off
                                                * in corresponding equation - see mm_as_const.h.
                                                * The index is over the equation number referenced
                                                * rf_fem_const.h   */
  int v[MAX_NUM_MATRICES][MAX_VARIABLE_TYPES]; /* Variable activity bit field
                                                * Bit - Purpose
                                                *   0 -> Variable isn't active nor is its value
                                                * even defined in the problem  1 -> This variable
                                                * occurs in the solution vector.  It is solved for.
                                                *   2 -> Variable is a constant in this material.
                                                *   4 -> Variable is not part of the solution
                                                * variable  for this element for this material,  but
                                                * it does vary across the domain.  Value is
                                                * calculated via interp from nodal values  8 ->
                                                * This variable is unique to this material.  It
                                                * will not be contiguous across material
                                                * boundaries.  At interfaces between materials,  the
                                                * value of the variable will have a discontinuity
                                                *        across the material interface.
                                                * -- see mm_as_const.h for more info
                                                */
  int mi[MAX_VARIABLE_TYPES];                  /* Matrix index for a given variable
                                                * -1 -> Not present in any matrix
                                                * >= 0 -> index into global matrix number
                                                */
  int gv[MAX_VARIABLE_TYPES];                  /* If this variable is on in any matrix (for field
                                                  variable access) 0 -> not in any matrix 1 -> in
                                                  a matrix
                                                */

  int w[MAX_NUM_MATRICES][MAX_EQNS];                   /* Weight function for equations */
  int i[MAX_NUM_MATRICES][MAX_VARIABLE_TYPES];         /* Interpolation type for each unknown
                                                        * in the current element block */
  int m[MAX_NUM_MATRICES][MAX_EQNS];                   /* Mapping from input file to real names. */
  dbl etm[MAX_NUM_MATRICES][MAX_EQNS][MAX_TERM_TYPES]; /* equation term multipliers */
  int CoordinateSystem;                                /* Cartesian, cylindrical, etc. */
  int MeshMotion;                                      /* Arbitrary or lagrangian or total ALE*/
  int MeshInertia;                                     /* addional inertia due to convection
                                                        * in the stress free state */
  int RealSolidFluxModel;                              /* linear or nonlinear */
  int MassFluxModel;                                   /* Fickian, Stefan-Maxwell, etc. */
  int MomentumFluxModel;                               /* Newtonian, Carreau, Powerlaw, etc. */
  int PorousFluxModel;                                 /* Fickian-Darcy, Darcy */
  int Num_Dim;                                         /* Number of spatial dimensions (2 or 3) */
  int TimeIntegration;                                 /* Steady, transient scheme */
  int Continuation;                                    /* First, second order scheme */
  int AugmentingConditions;                            /* Augmenting conditions */
  int IntegrationMap;                                  /* Iso- or sub-parametric mapping */
  int ShapeVar;                                        /* Variable whose basis functions will be
                                                        * used for manipulation of the geometry
                                                        * of the elements within this element
                                                        * block. */
  int ProjectionVar;                                   /* Variable whose basis functions will be
                                                        * used for projection of the field values
                                                        * onto the nodal values for elements within
                                                        * this element block */
  char MaterialName[MAX_MATLNAME];                     /* Names of Materials*/
  int Num_Species;                                     /* Number of species in present material */
  int Num_Species_Eqn;                                 /* Number of species Equations solved for in
                                                        * the present material - usually one less
                                                        * than the total number of species */
  int Species_Var_Type;                                /* Overrides of the default species var type
                                                        * for a particular material. CAUTION, may
                                                        * cause units problems at internal interfaces */
  int Num_Rxn;                                         /* Number of chemical reactions in present
                                                        * material */
  int VolumeIntegral;                                  /* Augmenting volume integral */
  int LSVelocityIntegral;   /* Augmenting sevel set velocity integral flag */
  int Num_Porous_Eqn;       /* number of porous media Equations */
  int Num_Porous_Shell_Eqn; /* number of porous media shell Equations */
  int Do_Surf_Geometry;     /* Problem needs a bundle of surface geometry defined on
                               it */
};
typedef struct Problem_Description PROBLEM_DESCRIPTION_STRUCT;

/*
 *  Define some common manipulations using the Problem_Description
 *  structure
 *
 */
#define VARIABLE_IN_THE_EB_SOLN_VECTOR(PDS, imtrx, var) ((PDS)->v[imtrx][(var)] & 1)

/*
 * External_Field_Variables Modified and overloaded this structure for use
 * in JAS coupling (1/20/2003) and external pixel fields to be mapped to
 *  mesh (11/29/2010)
 */

struct External_Field_Variables {
  int Num_external_field;       /* number of external fields read in and fixed*/
  int Num_external_pixel_field; /*number of external pixel fields in to be
                                  mapped*/
  int ev;                       /* external (fixed) field variable activity*/
  int ev_porous_decouple;       /* external field displacements from decoupled
                                   poroelastic flow activity */
  int ev_dpordt_index;          /* Index for the d porosity/dt external field */

  int ev_etch_area;  /* external area fraction in which etching reaction takes place */
  int ev_etch_depth; /* external depth field due to etching reaction */

  char name[MAX_EXTERNAL_FIELD][32];
  /* names of external field variables*/
  char file_nm[MAX_EXTERNAL_FIELD][MAX_FNL];
  /* names of exodus or pixel files with variables */
  int i[MAX_EXTERNAL_FIELD];           /* Interpolation of variables */
  int ipix[MAX_EXTERNAL_FIELD];        /* 0 for exoII file and 1 for pix file */
  int ipix_matid[MAX_EXTERNAL_FIELD];  /*Mat ID for pixel field to be mapped */
  dbl empty_value[MAX_EXTERNAL_FIELD]; /*Set field to this value if it's outside
                                          the voxel field DSB 7/30/13*/
  dbl *xyz_data[DIM];                  /*Array holding coordinates of each pixel */
  dbl *f_data;                         /*Array holding value of each pixel */
  dbl *ext_fld_ndl_val[MAX_EXTERNAL_FIELD];
  /* Array holding actual field nodal values */
  dbl *ext_fld_ndl_val_old[MAX_EXTERNAL_FIELD];
  dbl *ext_fld_ndl_val_older[MAX_EXTERNAL_FIELD];
  int TALE; /* boolean for whether TALE is active or not */
  dbl *init_displacement_ndl_val[2 * DIM];
  /* Array holding initial displacments for */
  /* mesh and real solid upon startup */
  /* These are needed for mesh annealing with */
  /* TALE, and are only allocated for TALE problems */
  char field_type[MAX_EXTERNAL_FIELD][15];
  /* type of external field to read; steady or transient */
  /* SMD 1/24/11 */
};

/*
 * Library_IO: Information passed into Goma from an external host code
 *		(applicable only when compiled with LIBRARY_MODE flag)
 */
struct Library_IO {
  double *xnv_in;     /* Imported nodal variables */
  double *xev_in;     /* Imported element variables */
  double *xsoln;      /* Exported solution variables */
  double *xpost;      /* Exported post-processing variables */
  int animas_step;    /* Call number to Goma */
  int print_flag;     /* How often to write ExodusII solution in Goma */
  int solve_steps;    /* How many Goma time steps to take per call */
  int goma_first;     /* Indicates if Goma is called before other code */
  double t_start;     /* Start time passed in from driver */
  double t_end;       /* End time passed in from driver */
  double last_step;   /* Size of last Animas time step */
  double decelerator; /* Decleration factor applied to Goma time step */
};

struct Transient_Information {
  /*
   * Contains constants used for transient analysis
   */
  int step;
  int MaxSteadyStateSteps;
  int MaxTimeSteps;
  int Fill_Weight_Fcn;                  /* Weight function to use on the transient fill equation
                                         */
  int Fill_Equation;                    /* Equation for fill-level set */
  dbl Delta_t0;                         /* initial time step */
  dbl Delta_t_min;                      /* minimum time step size */
  dbl Delta_t_max;                      /* maximum time step size */
  dbl time_step_decelerator;            /* factor used to make time step smaller when a
                                           time step fails to converge */
  dbl resolved_delta_t_min;             /* if dt < resolved_delta_t_min, accept any
                                           converged soln  regardless of time step error */
  dbl TimeMax;                          /* time at which to end integration */
  dbl theta;                            /* time step parameter: theta = 0. => Backward Euler
                                                                theta = 1. => Forward Euler
                                                                theta = .5 => Crack-Nicholson  */
  dbl eps;                              /* time step error  */
  int use_var_norm[MAX_VARIABLE_TYPES]; /* Booleans used for time step
                                           truncation error control */
  int fix_freq;
  int print_freq;
  int march_to_steady_state;     /* boolean if problem should be marched to steady
                                    state */
  double steady_state_tolerance; /* Tolerance for march to steady state */
  double print_delt;
  double print_delt2_time;
  double print_delt2;
  double init_time;
  int const_dt_after_failure;
  int Restart_Time_Integ_After_Renorm;

  /* Quantities of displacement acceleration.   This is added
   * here rather than in the normal way xdot and xdot_old are handled due
   * to the labyrinth of routines that require the argument list modification
   * just to get the quantity down to load_elem_dof_ptr.    What ought to
   * be done is to augment this struct with all of these quantities and then
   * just pass the pointer down once and for all.   viz. , x_old, x_older,
   * x_oldest, xdot, xdot_old, etc.    For now just make global
   */
  double delta_t;
  double *xdbl_dot;
  double *xdbl_dot_old;
  int solid_inertia;
  double newmark_beta;
  double newmark_gamma;
  double Courant_Limit;
  double time_value;
  double delta_t_avg;

  double delta_t_old;
  double time_value_old;

  int ale_adapt;
  int ale_adapt_freq;
  double ale_adapt_iso_size;
};

struct Eigensolver_Info
/*
 * Contains inputs to be used by either eggroll or ARPACK.
 */
{
  int Eigen_Algorithm;
  int Eigen_NEV_WANT;
  int Eigen_Maximum_Iterations;
  int Eigen_Maximum_Outer_Iterations;
  int Eigen_Filter;
  int Eigen_Krylov_Subspace;
  int Eigen_Recycle;
  int Eigen_Record_Modes;
  int Eigen_Matrix_Output;
  int Eigen_Solve_Freq;
  int Eigen_Write_Freq;
  dbl Eigen_Tolerance;
  dbl Eigen_IV_Wt;
  dbl Eigen_Shifts[4];
  dbl Eigen_Cayley_Sigma;
  dbl Eigen_Cayley_Mu;
  dbl Eigen_SI_Tol_Param;
  dbl Eigen_Relative_Tol;
  dbl Eigen_Linear_Tol;
  char Eigen_Output_File[85];
};

struct Continuation_Information {
  /*
   * Contains constants used for continuation analysis
   */
  int MaxPathSteps;
  int PathIntr;
  dbl Delta_s0;
  dbl Delta_s_min;
  dbl Delta_s_max;
  dbl PathMax;
  dbl alpha, beta, gamma, delta, theta;
  dbl eps;
  int use_var_norm[MAX_VARIABLE_TYPES];
  int print_freq;
  int fix_freq;
  double print_delt;
  double print_delt2_path;
  double print_delt2;
  double radius;
  double BegParameterValue;
  double EndParameterValue;
  double InitDir;
  /*  */
  int upType;
  /*  */
  int upBCID;
  int upDFID;
  int upDHID;
  /*  */
  int upMTID;
  int upMPID;
  int upMDID;
  int upMFID;
  /*  */
  int sensvec_id;
  /*  */
  double tmp1;
  double tmp2;
  double tmp3;
};

struct Loca_Input {
  /*
   * Contains inputs for LOCA.
   */

  int Cont_Alg;                   /* Specific LOCA algorithm - see ac_con_const.h */
  int Cont_Order;                 /* Continuation order: presently 0, 1, or 2     */
  double StepAggr;                /* Parameter for increasing step size           */
  double perturb;                 /* Perturbation size for bordering algorithms   */
  int debug;                      /* LOCA print level:  higher = more output      */
  double DpDs2;                   /* Desired solution contribution to arc length  */
  double DpDsHi;                  /* High value of dp_ds at which to rescale      */
  double Texp;                    /* Exponent used to calculate tangent factor    */
  double MaxTS;                   /* Maximum step change in tangent factor        */
  int TPupType;                   /* Turning point parameter type (BC or MT)      */
  int TPupBCID;                   /* ID tag of BC type turning point parameter    */
  int TPupDFID;                   /* Float ID of BC type turning point parameter  */
  int TPupMTID;                   /* Matl ID of MT type turning point parameter   */
  int TPupMPID;                   /* Property ID of MT type turning point parameter */
  int TPupMDID;                   /* Subindex ID of MT type turning point parameter */
  double TPGuess;                 /* Initial guess of parameter value at turning point */
  double TPFinal;                 /* Final TP parameter value			   */
  int NVRestart;                  /* Restart flag: read previous null vector if true */
  char NV_exoII_infile[MAX_FNL];  /* Exodus file name for null vector for starting*/
                                  /* TP or pitchfork tracking algorithm	   */
  char NV_imag_infile[MAX_FNL];   /* Exodus file name for null vector (imag. part) */
                                  /* for starting Hopf tracking algorithm         */
  int NVSave;                     /* Flag to save current TP/PF null vector       */
  char NV_exoII_outfile[MAX_FNL]; /* Exodus file name for saving final null vector*/
                                  /* from TP tracking algorithm		   */
  char NV_imag_outfile[MAX_FNL];  /* Exodus file name for saving imaginary      */
                                  /* part of null vector from Hopf algorithm      */
  int NV_time_index;              /* Time index to read Null vector from above file */
  float **PF_Nod_Vals;            /* Temporary array for storing nodal values of null
                                     vector read in from PF_exoII_file for
                                     pitchfork tracking runs */
  double *X_pitchfork;            /* Null vector for pitchfork tracking runs,
                                     dimensioned exactly like the solution vector x */
  float **HP_Nod_Vals;            /* Temporary array for storing nodal values of
                                     complex part of eigenvector read in from
                                     PF_exoII_file for Hopf tracking runs */
  double *X_hopf;                 /* Eigenvector for Hopf tracking runs,
                                     dimensioned exactly like the solution vector x */
  double omega;                   /* Imaginary part of Eigenvalue for Hopf tracking
                                     problems */
  int Mass_Derivatives;           /* Flag which determines whether to calculate
                                     Mass Matrix derivatives for Hopf tracking
                                     problems */
};

//! Structure containing parameter information for a single augmenting
//! condition
struct AC_Information {
  /*
   * Contains constants used for augmenting conditions
   */
  int nAC;
  /*  */
  int iread;
  dbl theta;
  dbl eps;
  /* Type of the augmented condition */
  int Type;

  /*
   *  Identification of the unknown in the augmentation condition
   *   - index of the BC where the unknown exists. Ordering is
   *     dependent on the ordering in the input deck
   */
  int BCID;
  /*
   *  Identification of the unknown in the augmentation condition
   *   - index of the float on the BC card that is the unknown
   */
  int DFID;

  // Node Set ID designating a position
  int DHID;

  /*
   * AC_VOLUME:
   *  volume constraint integers - sets the type of volume constraint
   *    1 volume
   *    2 mass
   *    3 mass of a particular species
   *
   *  AC_POSITION -> coordinate direction to be used
   */
  int VOLID;

  /*
   *  AC_VOLUME:
   * Integer parameter that identifies the species number to be
   * used when evaluating a vc_type 3 constraint equation.
   *  AC_POSITION:
   *  Integer parameter identifying the form of the residual equation.
   */
  int COMPID;

  /* Flux constraint integers */
  int SSID;
  int SSID2;
  int VAR;
  /* Level Set Velocity integers */
  int LSPHASE;
  int DIR;
  /*
   * Element block index of the material to be used in volume constraint
   */
  int MTID;
  int MPID;
  int MDID;
  int MFID;
  /*  float list */
  int len_AC;
  double *DataFlt;
  /*
   *  Current value of the unknown associated with this agumented condition
   */
  double tmp1;

  /*
   * Current valu eof the time derivative of the unknown associated with this
   * augmented condition
   */
  double tmp2;

  /* Old value of the unknown associated with this augmented condition
   *
   */
  double tmp3;

  /* volume constraint variables */
  double evol;

  /*
   * Constant to be used in the residual expression
   *
   * AC_POSITION : Absolute value of the real position
   */
  double CONSTV;

  double LewisNum;

  double *d_evol_dx; /*pointer to the derivative array */
  /* Level Set Velocity variables */
  double lsvel;
  double lsvol;
  double *d_lsvel_dx; /*point to derivative array wrt. F, x & v */
  double *d_lsvol_dx; /*point to derivative array wrt. F, x & v */
  /* These are for overlap BC's applied as AC's */
  int fluid_eb;
  int solid_eb;
  int lm_eb;
  int lm_elem;
  int lm_side;
  int lm_dim;
  double lm_value;
  double lm_resid;

  /*   file name and parameter name for aprepro parameters  */
  char Params_File[128];
  char AP_param[64];
  char *Aprepro_lib_string;
  int Aprepro_lib_string_len;
};

struct Continuation_Conditions {
  /*
   * Contains information for multiple conditions which depend on
   * a continuation parameter. (See former hunting conditions)
   */
  int nCC;
  /*  */
  dbl ratio;
  dbl old_value;
  dbl value;
  /*  */
  int Type;
  int fn_flag;
  /*  */
  int BCID;
  int DFID;
  /*
   *
   */
  int MTID;
  int MPID;
  int MDID;
  /*  */
  double Beg_CC_Value;
  double End_CC_Value;
  double coeff_0;
  double coeff_1;
  double coeff_2;
  /*  */
  int sensvec_id;
};

struct User_Continuation_Info {
  /*
   * Contains ID information only for multiple user-defined
   * continuation condition functions.
   */
  int nUC;
  int Type;
  int BCID;
  int DFID;
  int MTID;
  int MPID;
  int MDID;
  dbl old_value;
  dbl value;
};

struct HC_Information {
  /*
   * Contains constants used for hunting conditions
   */
  int nHC;
  /*  */
  dbl theta;
  dbl eps;
  /*  */
  int Type;
  int ramp;
  /*  */
  int BCID;
  int DFID;
  int DHID;

  /*
   * Material identification number (element block id in the mesh)
   */
  int MTID;

  int MPID;
  int MDID;
  int MFID;
  /*  */
  double BegParameterValue;
  double EndParameterValue;
  double Delta_s0;
  double Delta_s_min;
  double Delta_s_max;
  /*  */
  int sensvec_id;
  /*  */
  double tmp1;
  double tmp2;
  double tmp3;
};

/*____________________________________________________________________________*/

/*
 *  Basis_Functions Structure :
 *
 * Values of
 * 	(i) basis functions,
 *	(ii) local spatial derivatives of bfs.
 *	(iii) global physical space derivatives of bf's.
 *	(iv) mesh derivatives of (iii)
 */

struct Basis_Functions {
  int ielem_type;             /* old SHM identifier of elements... */
  int interpolation;          /* eg., I_Q1, ... */
  int element_shape;          /* eg., QUADRILATERAL, ...*/
  int Max_Dofs_Interpolation; /* How many degrees of freedom are involved
                               * in the interpolation of this element?
                               * For variable numbers of dofs, such as
                               * the SP interpolation, assume the maximum
                               * value. However, don't include dofs for
                               * variables that are not part of the
                               * interpolation for this element. Thus,
                               * I_Q1_D has 4 dofs */
  int *Var_Type_MatID;        /* Var_Type_MatID[mn] is the representative
                               * variable type that is interpolated
                               * using the current basis function in
                               * material index, mn. Note, we need a
                               * material index, because this value
                               * can vary between different materials */
  /*
   * load_basis_functions() fills in this stuff...
   */
  dbl phi[MDE];          /* phi_i */
  dbl dphidxi[MDE][DIM]; /* d(phi_i)/d(xi_j) */

  // Nedelec / vector Basis
  dbl phi_e[MDE][DIM];     /* vector phi_i e_k */
  dbl ref_phi_e[MDE][DIM]; /* vector phi_i e_k */
  dbl curl_e[MDE][DIM];
  dbl curl_phi[MDE][DIM];

  /*
   * beer_belly() fills in these elemental Jacobian things...
   */
  dbl J[DIM][DIM];
  /*
   *  determinant of the jacobian of the matrix transformation
   *  of the ShapeVar shape function.
   */
  int shape_dof;
  dbl detJ;
  dbl B[DIM][DIM]; /* inverse Jacobian */
  dbl d_det_J_dm[DIM][MDE];
  dbl dJ[DIM][DIM][DIM][MDE]; /* d( J[i][j] ) / d (d_k,l) */
  dbl dB[DIM][DIM][DIM][MDE];

  /*
   * These two things are the same in Cartesian coordinates, but not
   * in nonCartesian coordinate systems with nontrivial scale factors
   * and spatially-varying unit vectors...
   *
   * Strictly, e_a . grad(phi_i) =    1    d ( phi[i] )
   *				   ------  ------------
   *				    h[a]   d ( x_a )
   * where:
   *		h[a] == scale factors
   *		x_a  == physical coordinates (eg., z,r,theta)
   *
   *
   * Thus, there are two transformations...
   *
   *	  d phi[i]             d phi[i]		  1    d phi[i]
   *      --------    ---->    --------   ----> -----  --------
   *      d xi[j]	       d x[j]		 h[j]  d x[j]
   *
   *		    elemental		  scale
   *		    Jacobian		  factors
   */

  dbl d_phi[MDE][DIM];    /* d_phi[i][a]    = d(phi_i)/d(q_a) */
  dbl grad_phi[MDE][DIM]; /* grad_phi[i][a] = e_a . grad(phi_i) */

  dbl grad_phi_e[MDE][DIM][DIM][DIM]; /* grad_phi_e[i][a][p][q] */
                                      /* = (e_p e_q): grad(phi_i e_a) */

  /*
   *  curl_phi_e[i][a][p] = e_p dot curl(phi_i e_a)
   */
  dbl curl_phi_e[MDE][DIM][DIM];

  /*
   * d_d_phi_dmesh[i][a] [b][j] = d ( d_phi[i][a] )
   *				     --------------------
   *				     d ( d_b,j )
   */
  dbl d_d_phi_dmesh[MDE][DIM][DIM][MDE];

  /*
   * d_grad_phi_dmesh[i][a] [b][j] = d ( grad_phi[i][a] )
   *				     --------------------
   *				     d ( d_b,j )
   */

  dbl d_grad_phi_dmesh[MDE][DIM][DIM][MDE];

  /*
   * d_grad_phi_e_dmesh[i][a] [p][q] [b][j] = d ( grad_phi_e[i][a][p][q] )
   *					      ----------------------------
   *					      d ( d_b,j )
   *
   */
  dbl d_grad_phi_e_dmesh[MDE][DIM][DIM][DIM][DIM][MDE];
};
typedef struct Basis_Functions BASIS_FUNCTIONS_STRUCT;

/*____________________________________________________________________________*/

/*
 * These are field variables at the Gauss points of interest. They get loaded
 * up prior to each volume integration loop in each element. They might also
 * be loaded up on surface integration loops.
 *
 */

struct Field_Variables {
  dbl wt; /* Gauss weight. */

  dbl x[DIM];  /* Position in physical space. */
  dbl x0[DIM]; /* Initial Position in physical space. */

  /*
   * Add some useful quantities for curvilinear orthogonal coordinate
   * systems...note the difference between raw derivatives and the gradient
   * operator...(see mm_fill_aux.c for explanations of each of these variables)
   */
  dbl h[DIM];                          /* Scale factors. */
  dbl hq[DIM][DIM];                    /* Derivatives of scale factors. */
  dbl hqq[DIM][DIM][DIM];              /* 2nd derivatives of scale factors. */
  dbl curl_e[DIM][DIM];                /* Curl of unit vectors. */
  dbl d_curl_e_dq[DIM][DIM][DIM];      /* Derivative of Curl of unit vectors wrt q_b. */
  dbl grad_e[DIM][DIM][DIM];           /* Gradient of unit vectors. */
  dbl d_grad_e_dq[DIM][DIM][DIM][DIM]; /* 2nd derivatives of unit vectors. */
                                       /* Note this is not grad(grad(e_a)). */
  dbl h3;                              /* Volume element factor. */
  dbl dh3dq[DIM];                      /* Derivative of volume element factor */
  /* wrt each coordinate in this system.*/

  dbl dh3dmesh[DIM][MDE]; /* Derivative of volume element factor */
  /* wrt mesh displacement "b" with dof "j" */

  dbl T;               /* Temperature. */
  dbl v[DIM];          /* Velocity. */
  dbl v_star[DIM];     /* AUX Velocity, segregated */
  dbl pv[DIM];         /* Particle velocity. */
  dbl d[DIM];          /* Mesh displacement. */
  dbl x_first[DIM];    /* Initial mesh displacement on startup */
  dbl x_rs_first[DIM]; /* Initial solid displacement on startup */
  dbl d_rs[DIM];       /* real solid displacement. */
  dbl d_rs_first[DIM]; /* Initial solid displacement on startup */
  dbl c[MAX_CONC];     /* Concentration(s). */
  dbl P;               /* Pressure. */
  dbl P_star;
  dbl S[MAX_MODES][DIM][DIM]; /* Polymer Stress, for each mode */
  dbl G[DIM][DIM];            /* Velocity Gradient */
  dbl F;                      /* Fill */
  dbl V;                      /* Voltage */
  dbl qs;                     /* Surface charge density (shell element) */
  dbl SH;                     /* Shear rate from second invariant of rate-of-strain */
  dbl H;                      /* curvature of level set function */
  dbl n[DIM];                 /* LS function normal OR shell normal */
  dbl Enorm;                  /* potential field norm. */
  dbl p_liq;                  /* liquid-phase pressure, porous media variables(s). */
  dbl p_gas;                  /* gas-phase pressure, porous media variables(s). */
  dbl porosity;               /* porosity, porous media variables(s). */

  dbl vd[DIM]; /* Vorticity principle flow direction. */
  dbl vlambda; /* Eigenvalue associated with dv. */
  dbl nn;      /* This is the bond evolution*/

  dbl ext_v; /* Extension velocity */

  dbl E_field[DIM]; /* Electric field */

  dbl lm[DIM]; /* Lagrange Multiplier vector variable */

  dbl sh_K;                 /* Shell region curvature */
  dbl sh_K2;                /* Shell region second curvature */
  dbl sh_tens;              /* Shell region tension */
  dbl sh_x;                 /* Shell region x coordinate */
  dbl sh_y;                 /* Shell region y coordinate */
  dbl sh_u;                 /* Shell user */
  dbl sh_ang[DIM - 1];      /* Shell orientation angles */
  dbl div_s_v;              /* The scalar field evaluated on a shell element is div_s of v
                               or (( I - n n) dot del) dot v      */
  dbl curv;                 /* The scalar field evaluated on a shell element is the curvature */
  dbl grad_v_dot_n[DIM];    /* This vector field is the del_s v dotted into the
                               surface normal */
  dbl n_dot_curl_s_v;       /* n dot (curl_s v) Scalar variable used in shell
                               equations - curl_s is surface curl */
  dbl pF[MAX_PHASE_FUNC];   /* phase function */
  dbl sh_J;                 /* Shell surface diffusion flux */
  dbl sh_Kd;                /* Shell surface curvature */
  dbl apr, api, ars, sh_bv; /* Acoustic pressure */
  dbl epr, epi;             /* LAGR MULT EM continuity */
  dbl sink_mass;            /* porous sink mass */

  dbl external_field[MAX_EXTERNAL_FIELD];      /* External field to be read and held
                                                  const*/
  dbl grad_ext_field[MAX_EXTERNAL_FIELD][DIM]; /* Gradient of external field...just becuase */
  dbl initial_displacements[2 * DIM];          /* Initial displacements to be read and
                                                  held const */

  dbl sh_p;           /* shell lub approx. */
  dbl lubp;           /* lub approx. */
  dbl lubp_2;         /* lub_2 approx. */
  dbl sh_fp;          /* lub pressure approx in thin film */
  dbl sh_fh;          /* film thickness approx */
  dbl sh_pc;          /* particles concentration */
  dbl sh_sat_closed;  /* closed shell saturation - SAR */
  dbl sh_p_open;      /* open shell pressure - SAR */
  dbl sh_p_open_2;    /* open shell pressure 2 - PRS*/
  dbl sh_t;           /* shell temperature - PRS */
  dbl sh_dh;          /* shell delta h     - PRS */
  dbl sh_l_curv;      /* Lubrication shell curvature - SAR */
  dbl sh_l_curv_2;    /* Lubrication 2 shell curvature - SAR */
  dbl sh_sat_gasn;    /* shell saturation, gas compression - SAR */
  dbl sh_shear_top;   /* Top wall shear rate */
  dbl sh_shear_bot;   /* Bottom wall shear rate */
  dbl sh_cross_shear; /* Cross stream shear stress */
  dbl max_strain;     /* Maximum Von Mises strain */
  dbl cur_strain;     /* Von Mises strain */
  dbl poynt[DIM];     /* Poynting Vector */
  dbl tfmp_pres;      /* thin-film multi-phase lubrication pressure */
  dbl tfmp_sat;       /* thin-film multi-phase saturation */
  dbl restime;        /* residence time function field */
  dbl moment[MAX_MOMENTS];
  dbl rho;

  dbl em_er[DIM]; /* EM Electric Field Vector (real)*/
  dbl em_ei[DIM]; /* EM Electric Field Vector (imag)*/
  dbl em_hr[DIM]; /* EM Magnetic Field Vector (real)*/
  dbl em_hi[DIM]; /* EM Magnetic Field Vector (imag)*/
  dbl sh_sat_1;   /* Porous shell saturation layer 1 */
  dbl sh_sat_2;   /* Porous shell saturation layer 2 */
  dbl sh_sat_3;   /* Porous shell saturation layer 3 */

  dbl eddy_nu;       /* Eddy viscosity for turbulent flow */
  dbl wall_distance; /* Distance to nearest wall */

  /*
   * Grads of scalars...
   */

  dbl grad_T[DIM];                   /* Gradient of temperature. */
  dbl grad_P[DIM];                   /* Gradient of pressure. */
  dbl grad_P_star[DIM];              /* Gradient of pressure. */
  dbl grad_c[MAX_CONC][DIM];         /* Gradient of concentration(s). */
  dbl grad_moment[MAX_MOMENTS][DIM]; /* Gradient of moments */
  dbl grad_rho[DIM];                 /* Gradient of density. */
  dbl grad_F[DIM];                   /* Gradient of fill. */
  dbl grad_H[DIM];                   /* Gradient of curvature. */
  dbl grad_V[DIM];                   /* Gradient of voltage potential. */
  dbl grad_qs[DIM];                  /* Gradient of surface charge density. */
  dbl grad_SH[DIM];                  /* Gradient of shear rate from second invariant of
                                        rate-of-strain  */
  dbl grad_Enorm[DIM];               /* Gradient of the potential field norm. */
  dbl grad_p_liq[DIM];               /* Gradient of porous liq-phase pressure variable. */
  dbl grad_p_gas[DIM];               /* Gradient of porous gas-phase pressure variable. */
  dbl grad_porosity[DIM];            /* Gradient of porous  porosity variable. */
  dbl grad_nn[DIM];                  /* Gradient of bond evolution. */
  dbl grad_ext_v[DIM];               /* Extension velocity */
  dbl grad_sh_K[DIM];                /* Gradient of shell curvature */
  dbl grad_sh_K2[DIM];               /* Gradient of shell second curvature */
  dbl grad_sh_tens[DIM];             /* Gradient of shell tension */
  dbl grad_pF[MAX_PHASE_FUNC][DIM];  /* Gradient of phase function */
  dbl grad_sh_J[DIM];                /* Gradient of shell surface diffusion flux */
  dbl grad_apr[DIM], grad_api[DIM], grad_ars[DIM]; /* Gradient of Acoustic pressure */
  dbl grad_sh_bv[DIM];                             /* Gradient of shell boundary velocity	*/
  dbl grad_sh_p[DIM];                              /* Gradient of shell lub pressure       */
  dbl grad_lubp[DIM];                              /* Gradient of lub pressure       */
  dbl grad_lubp_2[DIM];                            /* Gradient of second lub pressure       */
  dbl grad_sh_fp[DIM];                             /* Gradient of lub pressure in the thin film */
  dbl grad_sh_fh[DIM];                             /* Gradient of film thickness */
  dbl grad_sh_pc[DIM];                             /* Gradient of particles concentration */
  dbl grad_sh_t[DIM];                              /* Gradient of shell temperature */
  dbl grad_sh_l_curv[DIM];                         /* Gradient of shell curvature */
  dbl grad_sh_l_curv_2[DIM];                       /* Gradient of shell curvature_2 */
  dbl grad_sh_p_open[DIM];                         /* Gradient of open porous shell pressure */
  dbl grad_sh_p_open_2[DIM];                       /* Gradient of open porous shell pressure */
  dbl grad_tfmp_pres[DIM]; /* Gradient of the thin-film multi-phase lubrication pressure */
  dbl grad_tfmp_sat[DIM];  /* Gradient of the thin-film multi-phase lubrication saturation */
  dbl grad_restime[DIM];   /* Gradient of the residence time function */
  dbl grad_sh_sat_1[DIM];  /* Gradient of porous shell saturation layer 1 */
  dbl grad_sh_sat_2[DIM];  /* Gradient of porous shell saturation layer 2 */
  dbl grad_sh_sat_3[DIM];  /* Gradient of porous shell saturation layer 3 */

  dbl grad_eddy_nu[DIM]; /* Gradient of Eddy viscosity */

  /*
   * Grads of vectors...
   */

  dbl div_v;                 /* Divergence of velocity. */
  dbl grad_v[DIM][DIM];      /* Gradient of velocity.  d (v_i) / d (x_j) */
  dbl div_v_star;            /* Divergence of velocity*. */
  dbl grad_v_star[DIM][DIM]; /* Velocity* segregated */
  dbl curl_v[DIM];           /* Curl of velocity, aka vorticity. */

  dbl div_pv;            /* Divergence of particle velocity. */
  dbl grad_pv[DIM][DIM]; /* Gradient of particle velocity. */

  dbl div_d;                /* Divergence of mesh displacement. */
  dbl div_d_dot;            /* Divergence of mesh velocity     */
  dbl grad_d[DIM][DIM];     /* Gradient of mesh displacement. */
  dbl grad_d_dot[DIM][DIM]; /* Gradient tensor of mesh velocity */

  dbl div_d_rs;            /* Divergence of solid displacement. */
  dbl grad_d_rs[DIM][DIM]; /* Gradient of solid displacement. */

  dbl grad_vd[DIM][DIM]; /* Gradient of vorticity principle flow direction. */
  dbl div_vd;            /* Divergence of vorticity direction. */

  dbl grad_E_field[DIM][DIM]; /* Electric field */

  dbl grad_n[DIM][DIM];  /* Normal to level set function OR shell normal */
  dbl d_n_dxi[DIM][DIM]; /* Derivative of normal w.r.t. isoparametric coordinates */

  dbl div_n;                         /* Divergence of LS normal field */
  dbl div_s_n;                       /* Surface divergence of LS normal field */
  dbl surfCurvatureDyadic[DIM][DIM]; /* Surface Curvature dyadic = b = - (I - n
                                        n ) grad(n) */
  dbl grad_poynt[DIM][DIM];          /* Gradient of Poynting.  d (P_i) / d (x_j) */
  dbl grad_em_er[DIM][DIM];          /* Gradient of EM Efield (real) */
  dbl grad_em_ei[DIM][DIM];          /* Gradient of EM Efield (imag) */
  dbl grad_em_hr[DIM][DIM];          /* Gradient of EM Hfield (real) */
  dbl grad_em_hi[DIM][DIM];          /* Gradient of EM Hfield (imag) */
  dbl curl_em_er[DIM];               /* Curl of EM Efield (real) */
  dbl curl_em_ei[DIM];               /* Curl of EM Efield (imag) */

  /* these gradients of tensors are complete for Cartesian coordinates,
   * and currently work for axisymmetic coordinates, in context,
   * but must be augmented for other coordinate systems ... we really need a
   * grad_phi_e_e!
   */

  dbl grad_S[MAX_MODES][DIM][DIM][DIM]; /* Gradient of polymer stress tensor( or most of it!) */
  dbl div_S[MAX_MODES][DIM];            /* Divergence of polymer stress tensor */
  dbl grad_G[DIM][DIM][DIM];            /* Gradient of velocity tensor ( or most of it!) */
  dbl grad_Gt[DIM][DIM][DIM];           /* Gradient of the transpose of the velocity tensor */
  dbl div_G[DIM];                       /* Divergence of velocity gradient tensor */
  dbl div_Gt[DIM]; /* Divergence of the transpose of velocity gradient tensor */

  dbl grad_n_dot_curl_s_v[DIM];            /* This is the normal gradient of a scalar field
                                              defined on a shell. The scalar field is n dot
                                              curl_s v or n dotted into the (I-nn) Del
                                              cross v see - apply_surface_viscosity()   */
  dbl grad_div_s_v[DIM];                   /* This is the normal gradient of a scalar field
                                              defined on a shell.        The scalar field is div_s of v
                                              or (( I - n n) dot del) dot v
                                              see - apply_surface_viscosity()   */
  dbl grad_curv[DIM];                      /* This is the normal gradient of a scalar field defined
                                             on a shell. The scalar field is curv  or - 1/2 Del_s dot
                                             v see - apply_surface_viscosity()   */
  dbl serialgrad_grad_s_v_dot_n[DIM][DIM]; /* This is the normal gradient of a
                                     scalar field defined on a shell. The scalar
                                     field is grad_s_v_dot_n[b]
                                     */

  dbl density;                     /* total density of material at gauss point */
  dbl d_density_dc[MAX_CONC][MDE]; /* Derivative of density wrt species unknown vector */
  dbl d_density_dmesh[DIM][MDE];   /* Derivative of density wrt mesh position
                                      unknown vector */
  dbl d_density_dT[MDE];           /* Derivative of density wrt temperature vector */
  dbl d_density_dP[MDE];           /* Derivative of density wrt pressure unknown vectors*/

  /*
   * Mesh derivatives of field variable gradients...
   * These require corresponding mesh derivatives of basis
   * functions to be evaluated.
   */

  dbl d_grad_T_dmesh[DIM][DIM][MDE];
  dbl d_grad_P_dmesh[DIM][DIM][MDE];
  dbl d_grad_nn_dmesh[DIM][DIM][MDE];

  dbl d_grad_V_dmesh[DIM][DIM][MDE];
  dbl d_grad_qs_dmesh[DIM][DIM][MDE];
  dbl d_grad_F_dmesh[DIM][DIM][MDE];

  dbl d_grad_moment_dmesh[MAX_MOMENTS][DIM][DIM][MDE];

  dbl d_grad_SH_dmesh[DIM][DIM][MDE];

  dbl d_grad_c_dmesh[DIM][MAX_CONC][DIM][MDE];

  dbl d_grad_ext_v_dmesh[DIM][DIM][MDE];
  dbl d_grad_E_field_dmesh[DIM][DIM][DIM][MDE];
  dbl d_grad_poynt_dmesh[DIM][DIM][DIM][MDE];
  dbl d_grad_em_er_dmesh[DIM][DIM][DIM][MDE];
  dbl d_grad_em_ei_dmesh[DIM][DIM][DIM][MDE];
  dbl d_grad_em_hr_dmesh[DIM][DIM][DIM][MDE];
  dbl d_grad_em_hi_dmesh[DIM][DIM][DIM][MDE];

  dbl d_grad_v_dmesh[DIM][DIM][DIM][MDE];
  dbl d_div_v_dmesh[DIM][MDE];

  dbl d_grad_n_dmesh[DIM][DIM][DIM][MDE];
  dbl d_div_n_dmesh[DIM][MDE];

  dbl d_grad_sh_K_dmesh[DIM][DIM][MDE];
  dbl d_grad_sh_K2_dmesh[DIM][DIM][MDE];
  dbl d_grad_sh_tens_dmesh[DIM][DIM][MDE];
  dbl d_grad_sh_J_dmesh[DIM][DIM][MDE];

  dbl d_grad_pv_dmesh[DIM][DIM][DIM][MDE];

  dbl d_grad_d_dmesh[DIM][DIM][DIM][MDE];
  dbl d_div_d_dmesh[DIM][MDE];
  dbl d_grad_d_dot_dmesh[DIM][DIM][DIM][MDE];

  dbl d_grad_d_rs_dmesh[DIM][DIM][DIM][MDE];
  dbl d_grad_d_rs_dd_rs[DIM][DIM][DIM][MDE];
  dbl d_div_d_rs_dmesh[DIM][MDE];

  dbl d_grad_S_dmesh[MAX_MODES][DIM][DIM][DIM][DIM][MDE];
  dbl d_div_S_dmesh[MAX_MODES][DIM][DIM][MDE];

  dbl d_grad_G_dmesh[DIM][DIM][DIM][DIM][MDE];
  dbl d_div_G_dmesh[DIM][DIM][MDE];
  dbl d_div_Gt_dmesh[DIM][DIM][MDE];

  dbl d_grad_p_liq_dmesh[DIM][DIM][MDE];
  dbl d_grad_p_gas_dmesh[DIM][DIM][MDE];
  dbl d_grad_porosity_dmesh[DIM][DIM][MDE];

  dbl d_grad_vd_dmesh[DIM][DIM][DIM][MDE];
  dbl d_div_vd_dmesh[DIM][MDE];

  dbl d_grad_n_dot_curl_s_v_dmesh[DIM][DIM]
                                 [MDE];    /* This is the mesh derivatives for grad(n_dot_curl_s_v)
                                              n_dot_curl_s_v is a shell variable.
                                              Therefore, this is only calculated on shell elements.
                                              The gradient is a full gradient, and the mesh unknowns
                                              refer to local unknowns on the shell element.  */
  dbl d_grad_div_s_v_dmesh[DIM][DIM][MDE]; /* This is the mesh derivatives for
                                              grad(div_s_v) div_s_v is a shell variable.
                                              Therefore, this is only calculated on shell
                                              elements. The gradient is a full gradient,
                                              and the mesh unknowns refer to local
                                              unknowns on the shell element.  */
  dbl d_grad_curv_dmesh[DIM][DIM][MDE];    /* This is the mesh derivatives for grad(curv)
                                              curv is a shell variable.
                                              Therefore, this is only calculated on shell
                                              elements. The gradient is a full gradient, and
                                              the mesh unknowns refer to local unknowns on
                                              the shell element.  */
  dbl d_serialgrad_grad_s_v_dot_n_dmesh[DIM][DIM][DIM]
                                       [MDE]; /* This is the mesh derivatives for
                                       grad(grad_s_v_dot_n) grad_s_v_dot_n is a shell vector
                                       variable. Therefore, this is only calculated on shell
                                       elements. The gradient is a full gradient, and the mesh
                                       unknowns refer to local unknowns on the shell element.  */

  dbl d_grad_apr_dmesh[DIM][DIM][MDE];
  dbl d_grad_api_dmesh[DIM][DIM][MDE];
  dbl d_grad_ars_dmesh[DIM][DIM][MDE];
  dbl d_grad_sh_bv_dmesh[DIM][DIM][MDE];
  dbl d_grad_sh_p_dmesh[DIM][DIM][MDE];
  dbl d_grad_lubp_dmesh[DIM][DIM][MDE];
  dbl d_grad_lubp_2_dmesh[DIM][DIM][MDE];
  dbl d_grad_sh_fp_dmesh[DIM][DIM][MDE];
  dbl d_grad_sh_fh_dmesh[DIM][DIM][MDE];
  dbl d_grad_sh_pc_dmesh[DIM][DIM][MDE];
  dbl d_grad_sh_t_dmesh[DIM][DIM][MDE];
  dbl d_grad_sh_l_curv_dmesh[DIM][DIM][MDE];
  dbl d_grad_sh_l_curv_2_dmesh[DIM][DIM][MDE];
  dbl d_grad_sh_p_open_dmesh[DIM][DIM][MDE];
  dbl d_grad_sh_p_open_2_dmesh[DIM][DIM][MDE];
  dbl d_max_strain_dmesh[DIM][MDE];
  dbl d_cur_strain_dmesh[DIM][MDE];
  dbl d_grad_eddy_nu_dmesh[DIM][DIM][MDE];
  dbl d_grad_restime_dmesh[DIM][DIM][MDE];
  /*
   * Values at surfaces for integrated boundary conditions
   */
  /* surface normal with coordinate scale vectors */
  double snormal[DIM]; /* Vector holding surface normal components     */
  double dsnormal_dx[DIM][DIM][MDE];

  /* vector for surface tangents - well defined in 2-D (only
   *   stangent[0] is calculated)
   * in 3D I try not to use these unless I'm on an edge,
   * in which case stangent[0] is the binormal (perpendicular to both snormal
   *     and the edge tangent stangent[1]) and stangent[1] is the parametric
   *   edge tangent along the edge  */
  double stangent[2][DIM];
  double dstangent_dx[2][DIM][DIM][MDE];

  /* surface determinant with coordinate scale vectors */
  dbl sdet;
  double dsurfdet_dx[DIM][MDE];

  /* edge determinant with coordinate scale vectors */
  dbl edge_det;
  double dedgedet_dx[DIM][MDE];

  /*   double cart_dsurfdet_dx[DIM][MDE]; */

  /* Deformation gradients and strain tensors that are used in
   * solid mechanics
   */
  dbl volume_change;
  dbl d_volume_change_dx[DIM][MDE];
  dbl d_volume_change_drs[DIM][MDE];
  dbl d_volume_change_dp[MDE];
  dbl volume_strain;
  dbl d_volume_strain_dx[DIM][MDE];
  dbl d_volume_strain_drs[DIM][MDE];
  dbl d_volume_strain_dp[MDE];
  dbl strain[DIM][DIM];
  dbl d_strain_dx[DIM][DIM][DIM][MDE];
  dbl d_strain_drs[DIM][DIM][DIM][MDE];
  dbl d_strain_dp[DIM][DIM][MDE];

  dbl deform_grad[DIM][DIM];
  dbl d_deform_grad_dx[DIM][DIM][DIM][MDE];

  dbl deform_grad_rs[DIM][DIM];
  dbl d_deform_grad_rs_dx[DIM][DIM][DIM][MDE];
  dbl d_deform_grad_rs_drs[DIM][DIM][DIM][MDE];

  /* The Stefan-Maxwell fluxes and the inverse of the Stefan-Maxwell
     flux-equation coefficient matrix
     as used in the thermal-battery model; KSC: 10/22/98 */
  dbl SM_flux[DIM * MAX_CONC];                       /* the Stefan_Maxwell flux vector */
  dbl SM_matrix_inv[DIM * MAX_CONC][DIM * MAX_CONC]; /* inverse of S-M flux-equation coff. matrix */

  dbl d_grad_tfmp_pres_dmesh[DIM][DIM][MDE];
  dbl giant_C_matrix[MAX_CONC][MDE][DIM * MAX_CONC][DIM * MAX_CONC]; /* matrix
          used to compute Jacobians in mm_fill_potential.c -- RSL 3/31/00 */
  dbl d_grad_tfmp_sat_dmesh[DIM][DIM][MDE];
};

/*
 * These are old and dot field variables at the Gauss points of interest.
 * Not all the information is needed for the old and dot terms ...
 */

struct Diet_Field_Variables {
  dbl x[DIM]; /* Position in physical space. */
  dbl T;      /* Temperature. */
  dbl v[DIM]; /* Velocity. */
  dbl v_star[DIM];
  dbl pv[DIM];     /* Particle velocity. */
  dbl d[DIM];      /* Mesh displacement. */
  dbl d_rs[DIM];   /* SOLID displacement. */
  dbl c[MAX_CONC]; /* Concentration(s). */
  dbl P;           /* Pressure. */
  dbl P_star;
  dbl F;                                  /* Fill. */
  dbl V;                                  /* Potential; added by KSC: 2/4/99 */
  dbl qs;                                 /* Surface charge density (shell element) */
  dbl Enorm;                              /* Norm of potential field. */
  dbl H;                                  /* Curvature of Level Set function */
  dbl n[DIM];                             /* normal vector to level set field OR shell normal */
  dbl S[MAX_MODES][DIM][DIM];             /* Polymer Stress, for each modes */
  dbl G[DIM][DIM];                        /* Velocity Gradient */
  dbl nn;                                 /* This is the bond evolution */
  dbl p_liq;                              /* porous media liq-pressure variable. */
  dbl p_gas;                              /* porous media gas-pressure variable. */
  dbl porosity;                           /* porous media porosity variable */
  dbl external_field[MAX_EXTERNAL_FIELD]; /* External field to be read and held
                                             const*/
  dbl ext_v;                              /* Extension velocity */
  dbl lm[DIM];
  dbl sh_K;            /*shell element curvature */
  dbl sh_K2;           /*shell element second curvature */
  dbl sh_tens;         /*shell element tension */
  dbl sh_x;            /*shell element x coordinate */
  dbl sh_y;            /*shell element y coordinate */
  dbl sh_u;            /* Shell user */
  dbl sh_ang[DIM - 1]; /* Shell orientation angles */
  dbl div_s_v;         /* sundry pieces (next 4) for surface rheological constitutive
                          eqn */
  dbl curv;
  dbl grad_v_dot_n[DIM];    /* grad_s_v_dot_n[DIM] */
  dbl n_dot_curl_s_v;       /* n dot (curl_s v) Scalar variable used in shell
                               equations - curl_s is surface curl */
  dbl pF[MAX_PHASE_FUNC];   /* Phase function */
  dbl sh_J;                 /* shell surface diffusion flux */
  dbl sh_Kd;                /* shell surface curvature */
  dbl apr, api, ars, sh_bv; /* Acoustic pressure */
  dbl epr, epi;
  dbl sink_mass;             /* porous sink mass */
  dbl sh_p;                  /* lub approx. */
  dbl lubp;                  /* lub approx. */
  dbl lubp_2;                /* lub 2 approx. */
  dbl grad_lubp[DIM];        /* lub pressure gradient approx */
  dbl grad_lubp_2[DIM];      /* lub pressure gradient approx */
  dbl sh_fp;                 /* lub pressure approx in the thin film */
  dbl grad_sh_fp[DIM];       /* lub pressure gradient approx in the thin film */
  dbl sh_fh;                 /* film thickness approx */
  dbl grad_sh_fh[DIM];       /* film thickness gradient approx */
  dbl sh_pc;                 /* particles concentration */
  dbl sh_sat_closed;         /* porous shell saturation - closed cells - SAR */
  dbl sh_p_open;             /* porous shell pressure - open cells - SAR */
  dbl sh_p_open_2;           /* porous shell pressure - open cells - SAR */
  dbl grad_sh_p_open[DIM];   /* gradient in porous shell pressure */
  dbl grad_sh_p_open_2[DIM]; /* gradient in porous shell pressure 2 */
  dbl sh_t;                  /* shell temperature */
  dbl sh_dh;                 /* shell delta h */
  dbl sh_l_curv;             /* Lubrication shell curvature - SAR */
  dbl sh_l_curv_2;           /* Lubrication shell curvature 2 - PRS */
  dbl sh_sat_gasn;           /* porous shell saturation - gas compression - SAR */
  dbl sh_shear_top;          /* Top wall shear rate */
  dbl sh_shear_bot;          /* Bottom wall shear rate */
  dbl sh_cross_shear;        /* Cross stream shear stress */
  dbl max_strain;            /* Maximum Von Mises strain */
  dbl cur_strain;            /* Von Mises strain */
  dbl poynt[DIM];            /* Poynting Vector */
  dbl tfmp_pres;             /* thin-film multi-phase lubrication pressure */
  dbl tfmp_sat;              /* thin-film multi-phase saturation */
  dbl moment[MAX_MOMENTS];
  dbl rho;
  dbl restime;    /* residence time field */
  dbl em_er[DIM]; /* EM wave Fields */
  dbl em_ei[DIM]; /* EM wave Fields */
  dbl em_hr[DIM]; /* EM wave Fields */
  dbl em_hi[DIM]; /* EM wave Fields */
  dbl sh_sat_1;   /* Porous shell saturation layer 1 */
  dbl sh_sat_2;   /* Porous shell saturation layer 2 */
  dbl sh_sat_3;   /* Porous shell saturation layer 3 */

  dbl eddy_nu; /* Eddy viscosity for turbulent flow */

  dbl grad_em_er[DIM][DIM]; /* EM wave Fields */
  dbl grad_em_ei[DIM][DIM]; /* EM wave Fields */
  dbl grad_em_hr[DIM][DIM]; /* EM wave Fields */
  dbl grad_em_hi[DIM][DIM]; /* EM wave Fields */
  /*
   * Gradients... concentration is the only one we use in the
   * old form for VOF/Taylor-Galerkin stuff
   */
  dbl grad_c[MAX_CONC][DIM];        /* Gradient of concentration(s). */
  dbl grad_F[DIM];                  /* Gradient of Fill variable. */
  dbl grad_pF[MAX_PHASE_FUNC][DIM]; /* Gradient of phase function */
  dbl div_v;                        // Divergence of velocity
  dbl grad_v[DIM][DIM];             // Gradient of velocity
  dbl grad_P[DIM];                  // Gradient of pressure
  dbl grad_P_star[DIM];
  dbl grad_p_liq[DIM];    /* Gradient of porous liq-phase pressure variable. */
  dbl grad_p_gas[DIM];    /* Gradient of porous gas-phase pressure variable. */
  dbl grad_porosity[DIM]; /* Gradient of porous  porosity variable. */

  dbl grad_T[DIM];         /* Gradient of porous  temperature variable. */
  dbl grad_d[DIM][DIM];    /* Gradient of mesh displacement. */
  dbl grad_d_rs[DIM][DIM]; /* Gradient of solid displacement. */

  dbl grad_tfmp_pres[DIM]; /* Gradient of the thin-film multi-phase lubrication pressure */
  dbl grad_tfmp_sat[DIM];  /* Gradient of the thin-film multi-phase lubrication saturation */

  dbl grad_n[DIM][DIM]; /* Normal to level set function OR shell normal */
  dbl div_n;            /* Divergence of LS normal field */

  /* Material tensors used at old time values */
  dbl strain[DIM][DIM]; /* Strain tensor */
  dbl volume_change;    /* Volume change */
  dbl volume_strain;
  dbl deform_grad[DIM][DIM];
  dbl d_deform_grad_dx[DIM][DIM][DIM][MDE];
  dbl d_volume_change_dx[DIM][MDE];
  dbl d_volume_strain_dx[DIM][MDE];
  dbl d_strain_dx[DIM][DIM][DIM][MDE];
  dbl d_grad_d_dmesh[DIM][DIM][DIM][MDE];
  dbl d_grad_d_dot_dmesh[DIM][DIM][DIM][MDE];
  dbl grad_restime[DIM]; /* Gradient of the Residence time field */

  dbl grad_moment[MAX_MOMENTS][DIM];
};

struct Rotation_Vectors {
  double vector[DIM];                  /* THREE vectors used in rotation */
  double d_vector_dx[DIM][DIM][MNROT]; /* sensitivity w.r.t. global displacements */
  int d_vector_J[MNROT];               /* global node numbers of the displacement
                                        * sensitivities     can be different for each rotation
                                        * vector */
  int d_vector_n;                      /* number of global node numbers in sensitivity */
  int ok; /* flag indicating that a rotation vector has been calculated */
};
typedef struct Rotation_Vectors ROTATION_VECTORS_STRUCT;

/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/
/*
 * Intermediary Variables used in calculating flow through a porous media.
 * These values are calculated at the current quadrature point.
 *
 * Our current variable types are:
 *                   POR_LIQ_PRES
 *                   POR_GAS_PRES
 *                   POR_POROSITY
 *                   POR_TEMP
 *                   POR_SATURATION
 *
 * (Don't know of a problem that actually uses POR_SATURATION as a variable
 *  type, yet).
 */

struct Porous_Media_Variables {
  double cap_pres;

  /* Recently added for two-phase flow in porous media, nonisothermal */

  double enthalpy[3];
  double d_enthalpy[3][MAX_VARIABLE_TYPES + MAX_CONC];
  double d_d_enthalpy[3][MAX_VARIABLE_TYPES + MAX_CONC][MAX_VARIABLE_TYPES + MAX_CONC];

  double rhog, d_rhog[MAX_PMV], d_drhog[MAX_PMV][MAX_PMV];
  double d_Ywg[MAX_PMV], d_dYwg[MAX_PMV][MAX_PMV];
  double d_Yag[MAX_PMV], d_dYag[MAX_PMV][MAX_PMV];

  /* Before we launch into the difficult task of sorting through all of the
   * relevant variables required for the porous media calculations, we must
   * understand some fundamentals of phases and components, at least for all of
   * this to make sense. First of all, our primary variables p_liq, p_gas, and
   * porosity really represent the 3 phases we are restricting our calculations
   * to, viz. liquid, gas, and solid. We are fully well aware of Gibb's phase
   * rule, and realize that in multicomponent liquids you can have more than one
   * liquid phase, but WE are not allowing for such things right now.  The Darcy
   * equations we solve for each of these phases is really an overall mass
   * balance for that phase.  Yet another way to think of them are as component
   * balances for N-1 component of the phase, which we will refer to as the
   * "solvent".
   *
   * Let's take an example. For a solid/water/air system, where the solid is
   * deformable, we must track 3 phases.   The primitive variables we have
   * chosen for each of these phases are p_liq, p_gas, and porosity,
   * respectively.  Trust us, there are consitutive equations relating these
   * quantities to volume fraction of each phase in the mixture. So in the
   * liquid phase we have Water, which is single component, but in the gas phase
   * we have air and water vapor (2 components) and in the solid phase we assume
   * just one insoluble component.  To account for all concentrations of all
   * components at all places and at all times we need 1 equation for the solid
   * phase, two for the gas phase, and one for the liquid phase. If we further
   * assume phase equilibrium between water liquid and water vapor throughout,
   * we have one additional equation to the overall balances for p_liq, p_gas,
   * and porosity, hence 4 total.  In GOMA, this this equilibriumis assumed,
   * viz. we do not allow super-heated steam near liquid regions.
   *
   * Note for that lucky developer who gets to extend this for multicomponent
   * liquids, and hence more than 2 components in the gas: for each additional
   * liquid-phase volatile component you will need to add an additional species
   * transport, convective diffusion equation much like R_SPECIES, where the
   * velocity field will come from the overall Darcy-law velocity.   The place
   * holders for these probably would be in the PMV structures, and not the
   * species structures, due to the fact that they will have different
   * multipliers on them and different provisions for a discontinuous phase.
   *
   * A final note: much of what was here was geared towards this example, with
   * phases-and components interchanged.  The vestiges of that code still exist
   * with this MAX_PMV stuff below, as we are looping over the
   * phases/components.
   */
  /*
   * gas_density_solvents[i] - This is the local density of the "solvent"
   * component, i, in the GAS phase. The index is over the solvent for each
   * phase represented by the porous media variables index. The units for this
   * term are gm cm-3.
   *
   *        i = i_gas ->  gas_density_solvents[i_gas] = rho_g * Y_air
   *        i = i_pore->  gas_density_solvents[i_pore]=0 since solid matrix
   *                                             is insoluble in gas
   *        i = i_liq ->  gas_density_solvents[i_liq] = rho_g * Y_water
   *                               (where Y_water is calculated from an
   *                                equilibrium expression)
   *                     (rho_g is density of gas phase and Y_* are
   *                      mass fractions)
   *
   */
  double gas_density_solvents[MAX_PMV];
  double d_gas_density_solvents[MAX_PMV][MAX_VARIABLE_TYPES + MAX_CONC];

  /*
   * liq_Xvol_solvents[i] - This is the volume fraction of the solvent
   component, i,
   *                  in the liquid phase. The index is over the solvent for
   *                  each phase represented by the porous media variables
   index.
   *
   *        i = i_gas ->  liq_xvol_solvents[i_gas] =0 since air is insoluble in
   water
   *        i = i_pore->  liq_Xvol_solvents[i_pore]=0 since solid matrix
   *                                             is insoluble in water
   *        i = i_liq ->  liq_Xvol_solvents[i_liq]= 1.0 since solvent species
   makes up
   *                             the whole solvent phase until multicomponent
   *                             capability is installed
   *
   * -> Thus liq_Xvol_solvents[] is largely a placeholder until more complexity
   is    *    added.
   *
   *
   *    For the energy equation, this variable holds the liquid enthalpy to
   maintain the
   *    similarities in the equation structure, thus allowing looping with
   noncontributing
   *    terms set to zero.



   *
   */
  double liq_Xvol_solvents[MAX_PMV];
  double d_liq_Xvol_solvents[MAX_PMV][MAX_VARIABLE_TYPES + MAX_CONC];

  /*
   * bulk_density[i] - This is the local density of the solvent component, i,
   *                in all phases. The index is over the solvent for
   *                each phase represented by the porous media variables index.
   *                The units for this term are gm cm-3.
   *
   *        i = i_gas ->  bulk_density[i_gas]
   *        i = i_pore->  bulk_density[i_pore]
   *        i = i_liq ->  bulk_density[i_liq]
   *
   */
  double bulk_density[MAX_PMV];
  double d_bulk_density[MAX_PMV][MAX_VARIABLE_TYPES + MAX_CONC];
  /*
   * gas_darcy_velocity[a] - This is the Darcy velocity of the gas
   *             phase due to Darcy's flow possibly modified by
   *             the gravitational body force term.
   *             units - cm/sec
   */
  double gas_darcy_velocity[DIM];
  double d_gas_darcy_velocity[DIM][MAX_PMV][MDE];
  /*
   * liq_darcy_velocity[a] - This is the Darcy velocity of the liquid
   *             phase due to Darcy's flow possibly modif
   *             the gravitational body force term.
   *             units - cm/sec
   */
  double liq_darcy_velocity[DIM];
  double d_liq_darcy_velocity[DIM][MAX_PMV][MDE];
  double d_liq_darcy_velocity_dSM[DIM][MDE];

  /* variables for special pore models */
  double r_pore;
  double d_r_pore[MAX_VARIABLE_TYPES + MAX_CONC];
  double d_d_r_pore[MAX_VARIABLE_TYPES + MAX_CONC][MAX_VARIABLE_TYPES + MAX_CONC];
  double r_cap;
  double d_r_cap[MAX_VARIABLE_TYPES + MAX_CONC];
  double d_d_r_cap[MAX_VARIABLE_TYPES + MAX_CONC][MAX_VARIABLE_TYPES + MAX_CONC];

  /* need second derivatives (with cross-terms) of some of these 'derived'
   * quantities so that we can get analytical Jacobians of gradients of these
   * quantities
   */
  double d_d_gas_vol_frac[MAX_PMV][MAX_VARIABLE_TYPES + MAX_CONC][MAX_VARIABLE_TYPES + MAX_CONC];
  double d_d_gas_density_solvents[MAX_PMV][MAX_VARIABLE_TYPES + MAX_CONC]
                                 [MAX_VARIABLE_TYPES + MAX_CONC];

  double d_d_liq_Xvol_solvents[MAX_PMV][MAX_VARIABLE_TYPES + MAX_CONC]
                              [MAX_VARIABLE_TYPES + MAX_CONC];

  double d_d_bulk_density[MAX_PMV][MAX_VARIABLE_TYPES + MAX_CONC][MAX_VARIABLE_TYPES + MAX_CONC];
  double d_d_bulk_density_dt[MAX_PMV][MAX_VARIABLE_TYPES + MAX_CONC][MAX_VARIABLE_TYPES + MAX_CONC];

  double rel_mass_flux[MAX_PMV][DIM];
  double d_rel_mass_flux_dpmv[MAX_PMV][DIM][MAX_PMV][MDE];
  double d_rel_mass_flux_dmesh[MAX_PMV][DIM][DIM][MDE];
  double d_rel_mass_flux_dT[MAX_PMV][DIM][MDE];
  double d_rel_mass_flux_dSM[MAX_PMV][DIM][MDE];

  double U_supg[DIM];
  double U_supg_hnorm[DIM];
  double d_U_supg_hnorm_dmde[DIM][MDE];
  double h_veloc_lcd[DIM];
  double zeta;
  double k_art_diff;
  double U_supg_squared;
};

/************************************************************************/
/************************************************************************/
/************************************************************************/
/*
 *  Structure used in calculating the individual terms in the
 *  porous media equations, when Mass Lumping is employed.
 *
 *  The inner loop is always over the degree of freedom.
 */
struct Porous_Media_Variables_ML {

  double Bulk_Density[MDE][MAX_PMV];
  double d_Bulk_Density[MDE][MAX_PMV][MAX_VARIABLE_TYPES];
  double Bulk_Density_old[MDE][MAX_PMV];
  double d_Bulk_Density_old[MDE][MAX_PMV][MAX_VARIABLE_TYPES];

  double Inventory_Solvent[MDE][MAX_PMV];
  double Inventory_Solvent_old[MDE][MAX_PMV];
  double Inventory_Solvent_dot[MDE][MAX_PMV];
  double Inventory_Solvent_dot_old[MDE][MAX_PMV];
  double d_Inventory_Solvent_dot_dpmv[MDE][MAX_PMV][MAX_PMV];
};
typedef struct Porous_Media_Variables_ML PMV_ML_STRUCT;

/************************************************************************/
/************************************************************************/
/************************************************************************/
/*
 * Structure used in calculating the individual terms in the
 * porous media equations.
 *
 * Big NOTE: Really these are the quantities
 * which make up the individual terms of Darcy's equations.   Really,
 * these equations are being used to solve for p_liq, p_gas, and porosity,
 * but the terms are really pieces of a component mass balance.
 * WE MUST DEFINE WHAT A COMPONENT/SOLVENT IS FOR THIS TO ALL MAKE SENSE!
 *
 *
 * We hereafter define a "SOLVENT" as the primary component in each "PHASE"
 * for which the Darcy equation is accounting for.  In our simple case of
 * Water/Air/Solid, the "Solvent" for the liquid phase is water, the "solvent"
 * for the gas phase is "Air", and the solvent for the solid phase is, of
 * course, the solid.  Our component balances are for water, air, and solid in
 * this partially saturated case, NOT for liquid/gas/solid.   It just so happens
 * that in single species component liquid and gas, a "phase" and a "component"
 * are one in the same, hence the confusion in the equations.
 *
 * When we say "Inventory_solvent[i_pl]", we mean the concentration of the
 * primary component of the original Liquid phase in both GAS and LIQUID at that
 * point (assuming it is insoluble in the solid.  In the water/air system, this
 * means water.   This concept is a must for multicomponent liquids and gases,
 * for which we will be augmenting the system with additional species equations.
 *       PRS (5/9/01)
 */
struct Porous_Media_Terms {
  dbl Inventory_solvent[MAX_PMV];         /* Gas+liquid inventory of "solvent" for each
                                             phase */
  dbl Inventory_solvent_old[MAX_PMV];     /* value at last time step of quantity
                                             above */
  dbl Inventory_solvent_dot[MAX_PMV];     /* quantity above wrt time */
  dbl Inventory_solvent_dot_old[MAX_PMV]; /* value at last time step of quantity
                                             above */
  dbl Inventory_solvent_dot_dc[MAX_PMV][MAX_CONC][MDE];
  dbl d_Inventory_solvent_dot_dpmv[MAX_PMV][MAX_PMV][MDE];
  /* Porous Media Unknown Vector time
     derivative wrt porous media vars. */
  dbl d_Inventory_sol_dpmv[MAX_PMV][MAX_PMV][MDE];
  /*sensitivity of liq solvent inventory wrt porosity */
  dbl d_PM_dot_dP[MAX_PMV][MDE]; /* Porous Media Unknown Vector at last time step. */
  dbl grad_PM[MAX_PMV][DIM];     /* Porous Media Unknown Vector gradient. */
  dbl grad_PM_old[MAX_PMV][DIM]; /* Porous Media Unknown Vector gradient at last
                                    time step*/

  /* Here, rather than unrolling these, if the first dimension is 0, then we
   * mean the flux of liquid solvent, 1, then we mean gas solvent, or 2 we mean
   * solid solvent See comment above for this rather bizarre definition of
   * solvent.
   */
  dbl diff_flux[MAX_PMV][DIM]; /* Diffusive flux vector -> gm cm-2 sec-1 */
  dbl d_diff_flux_dc[MAX_PMV][DIM][MAX_CONC][MDE];
  dbl d_diff_flux_dpmv[MAX_PMV][DIM][MAX_PMV][MDE];
  dbl d_diff_flux_dmesh[MAX_PMV][DIM][DIM][MDE];
  dbl d_diff_flux_dv[MAX_PMV][DIM][DIM][MDE];
  dbl d_diff_flux_dT[MAX_PMV][DIM][MDE];
  dbl d_diff_flux_dSM[MAX_PMV][DIM][MDE];

  dbl taylor_flux[MAX_PMV][DIM]; /* Taylor Galerkin Diff.-like flux */
  dbl d_taylor_flux_dc[MAX_PMV][DIM][MAX_CONC][MDE];
  dbl d_taylor_flux_dpmv[MAX_PMV][DIM][MAX_PMV][MDE];
  dbl d_taylor_flux_dmesh[MAX_PMV][DIM][DIM][MDE];
  dbl d_taylor_flux_dv[MAX_PMV][DIM][DIM][MDE];
  dbl d_taylor_flux_dT[MAX_PMV][DIM][MDE];

  dbl taylor_flux_wt[MDE]; /* Taylor Galerkin wt fnc. */
  dbl d_taylor_flux_wt_dmesh[MDE][DIM][MDE];
  dbl d_taylor_flux_wt_dv[MDE][DIM][MDE];
  dbl d_taylor_flux_wt_dT[MDE][MDE];

  dbl conv_flux[MAX_PMV][DIM]; /* convection flux vector. */
  dbl d_conv_flux_dc[MAX_PMV][DIM][MAX_CONC][MDE];
  dbl d_conv_flux_dpmv[MAX_PMV][DIM][MAX_PMV][MDE];
  dbl d_conv_flux_dmesh[MAX_PMV][DIM][DIM][MDE];
  dbl d_conv_flux_dv[MAX_PMV][DIM][DIM][MDE];
  dbl d_conv_flux_dT[MAX_PMV][DIM][MDE];
  dbl d_conv_flux_dSM[MAX_PMV][DIM][MDE];

  dbl MassSource[MAX_PMV]; /* source . */
  dbl d_MassSource_dc[MAX_PMV][MAX_CONC][MDE];
  dbl d_MassSource_dpmv[MAX_PMV][MAX_PMV][MDE];
  dbl d_MassSource_dmesh[MAX_PMV][DIM][MDE];
  dbl d_MassSource_dv[MAX_PMV][DIM][MDE];
  dbl d_MassSource_dpv[MAX_PMV][DIM][MDE];
  dbl d_MassSource_dT[MAX_PMV][MDE];
  dbl d_MassSource_dsh[MAX_PMV][MDE];
  dbl d_MassSource_dV[MAX_PMV][MDE];
  dbl d_MassSource_dSM[MAX_PMV][MDE]; /*sink mass sensitivity */

  dbl pi_supg[MDE];
  dbl d_pi_supg_dpmv[MDE][MAX_PMV][MDE];
  dbl conv_flux_supg[MAX_PMV];
  dbl d_conv_flux_supg_dpmv[MAX_PMV][MAX_PMV][MDE];
};

/************************************************************************/
/************************************************************************/
/************************************************************************/
/*
 *  Structure used to store additional nodal variables
 *  used in hysteretic porous media equations
 *
 *
 */
struct Porous_Media_Variables_Hysteresis {

  /*  Each member has dimension of [MAX_POR_SHELL][NUM_NODES]  */

  int *curve_type[MAX_POR_SHELL]; /*Draining = 1, imbibition = 0*/
  int *curve_type_old[MAX_POR_SHELL];
  int *curve_switch[MAX_POR_SHELL]; /*Trigger to switch the curves */
  int *num_switch[MAX_POR_SHELL];   /*Number of curve switch on each node*/
  double *sat_switch[MAX_POR_SHELL] /*Saturation value at switch point*/;
  double *cap_pres_switch[MAX_POR_SHELL]; /*Capillary pressure value at switch point*/
  double *sat_min_imbibe[MAX_POR_SHELL];  /* Minimum saturation value of imbibition curve*/
  double *sat_max_drain[MAX_POR_SHELL];   /* Maximum saturation value of draining curve*/
};

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/*
 *  Common parameters associated with stabilization schemes for
 *  SUPG and PSPG
 *
 */
struct Stabilization_Params {

  double hsquared[DIM]; /* square of the element size variable */
  double hhv[DIM][DIM]; /* Vector of directional element sizes
                         * first coordinate is the local element
                         * coordinate number, while the second is
                         * the cartesian coordinate number.  Thus,
                         * hhv[p] is the midpoint to midpoint vector
                         * of the element in the p'th local
                         * element coordinate direction */
  double dhv_dxnode[DIM][MDE];
  double h_veloc_elem;
  double Grid_Peclet_Number[MAX_VARIABLE_TYPES + MAX_CONC];
};
typedef struct Stabilization_Params STABILIZATION_PARAMS_STRUCT;

/*____________________________________________________________________________*/

/*
 * These are a big deal now.
 */
struct Constitutive_Relations {
  int HeatFluxModel;
  int MeshFluxModel;
  int RealSolidFluxModel;
  int MeshMotion;
  int MassFluxModel;
  int MomentumFluxModel;
  int PorousFluxModel;
};

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/*
 * Structure used in calculating the individual terms in the species
 * conservation equations so that it is general for both continuous and porous
 * media
 */
struct Species_Conservation_Terms {
  dbl Y[MAX_CONC];         /* Species Unknown Vector */
  dbl Y_old[MAX_CONC];     /* Species Unknown Vector at last time step. */
  dbl Y_dot_old[MAX_CONC]; /* old Species Unknown Vector derivative wrt
                            * time.                                     */
  dbl Y_dot[MAX_CONC];     /* Species Unknown Vector derivative wrt time*/

  dbl d_Y_dot_dc[MAX_CONC][MAX_CONC][MDE];
  dbl d_Y_dot_dpmv[MAX_CONC][MAX_PMV][MDE];
  dbl d_Y_dot_dP[MAX_CONC][MDE];

  dbl grad_Y[MAX_CONC][DIM];     /* Species Unknown Vector gradient.          */
  dbl grad_Y_old[MAX_CONC][DIM]; /* Species Unknown Vector gradient at last   *
                                  * time step.                                */
  dbl diff_flux[MAX_CONC][DIM];  /* Diffusion flux vector.                    */
  dbl d_diff_flux_dc[MAX_CONC][DIM][MAX_CONC][MDE];
  dbl d_diff_flux_dmesh[MAX_CONC][DIM][DIM][MDE];
  dbl d_diff_flux_dv[MAX_CONC][DIM][DIM][MDE];
  dbl d_diff_flux_dvd[MAX_CONC][DIM][DIM][MDE];
  dbl d_diff_flux_dT[MAX_CONC][DIM][MDE];
  dbl d_diff_flux_dP[MAX_CONC][DIM][MDE];
  dbl d_diff_flux_dV[MAX_CONC][DIM][MDE]; /* derivative of total flux wrt potential */
  dbl d_diff_flux_dSH[MAX_CONC][DIM][MDE];
  dbl d_diff_flux_dG[MAX_CONC][DIM][DIM][DIM][MDE];
  dbl d_diff_flux_dpmv[MAX_CONC][DIM][MAX_PMV][MDE];

  dbl taylor_flux[MAX_CONC][DIM]; /* Taylor Galerkin Diff.-like flux */
  dbl d_taylor_flux_dc[MAX_CONC][DIM][MAX_CONC][MDE];
  dbl d_taylor_flux_dmesh[MAX_CONC][DIM][DIM][MDE];
  dbl d_taylor_flux_dv[MAX_CONC][DIM][DIM][MDE];
  dbl d_taylor_flux_dT[MAX_CONC][DIM][MDE];
  dbl d_taylor_flux_dP[MAX_CONC][DIM][MDE];
  dbl d_taylor_flux_dpmv[MAX_CONC][DIM][MAX_PMV][MDE];

  dbl taylor_flux_wt[MDE]; /* Taylor Galerkin wt fnc. */
  dbl d_taylor_flux_wt_dmesh[MDE][DIM][MDE];
  dbl d_taylor_flux_wt_dv[MDE][DIM][MDE];
  dbl d_taylor_flux_wt_dT[MDE][MDE];

  dbl conv_flux[MAX_CONC][DIM]; /* convection flux vector. */
  dbl d_conv_flux_dc[MAX_CONC][DIM][MAX_CONC][MDE];
  dbl d_conv_flux_dmesh[MAX_CONC][DIM][DIM][MDE];
  dbl d_conv_flux_dv[MAX_CONC][DIM][DIM][MDE];
  dbl d_conv_flux_dT[MAX_CONC][DIM][MDE];
  dbl d_conv_flux_dP[MAX_CONC][DIM][MDE];
  dbl d_conv_flux_dpmv[MAX_CONC][DIM][MAX_PMV][MDE];

  dbl MassSource[MAX_CONC]; /* source . */
  dbl d_MassSource_dc[MAX_CONC][MAX_CONC][MDE];
  dbl d_MassSource_dmesh[MAX_CONC][DIM][MDE];
  dbl d_MassSource_dv[MAX_CONC][DIM][MDE];
  dbl d_MassSource_dpv[MAX_CONC][DIM][MDE];
  dbl d_MassSource_dT[MAX_CONC][MDE];
  dbl d_MassSource_dP[MAX_CONC][MDE];
  dbl d_MassSource_dsh[MAX_CONC][MDE];
  dbl d_MassSource_dI[MAX_CONC][MDE];
  dbl d_MassSource_dV[MAX_CONC][MDE];
  dbl d_MassSource_dpmv[MAX_CONC][MAX_PMV][MDE];
  dbl d_MassSource_dF[MAX_CONC][MDE];
  dbl d_MassSource_drst[MAX_CONC][MDE];
};
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
//! A structure filled with common arguments.
/*!
 * This structure is populated at the top of solve_nonlinear_problem()
 */
struct Matrix_Fill_Arguments {
  struct GomaLinearSolverData *ams;
  double *x;
  double *resid;
  double *x_old;
  double *x_older;
  double *xdot;
  double *xdot_old;
  double *x_update;
  double *delta_t;
  double *theta_;
  struct elem_side_bc_struct **first_elem_side_bc;
  double *time;
  Exo_DB *exo;
  Dpi *dpi;
  int *num_total_nodes;
  double *h_elem_avg;
  double *U_norm;
  double *estifm;
};
typedef struct Matrix_Fill_Arguments MF_Args;

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/*
 *  Note: we have already used the assumption that the species dependencies
 *        are full and contiguous in the structure below. In other words
 *        they are always grouped together and there are mp->Num_Species
 *        number of them.
 */
struct propertyJac {
  int NUM_TERMS_MALLOC; /* Number of terms, used in the malloc of this
                         * structure
                         */
  int Num_Terms;        /* Number of terms in the Jacobian dependence
                         * vector below.
                         */
  int Species_Type;     /* Value of the species var type for
                         * the source term and dependent variable assumed
                         * in the Fields below.
                         */
  VARIABLE_DESCRIPTION_STRUCT **Var_List;
  /* Vector of pointers to variable description
   * structures for each of the Jacobian terms
   * (malloced length of NUM_TERMS_MALLOC);
   */
  int *idof;     /* Vector of idof's. This is usually equal to
                  * zero, except for Variable descriptions
                  * that pertain to more than one dof.
                  * Also doubles for kspec for vintage MASS_FRACTION
                  * variable types
                  * (malloced length of NUM_TERMS_MALLOC);
                  */
  int *Var_Type; /* Var type of of the vector.
                  * (malloced length of NUM_TERMS_MALLOC);
                  */
  int *MatID;    /* MatID of the variable
                  *(malloced length of NUM_TERMS_MALLOC);
                  */
  double Property_Value;
  /* Current value of the property whose
   * deriviative is being evaluated.
   */
  double *Var_Value; /* Current value of the variable corresponding
                      * to the entry. Species var type agrees with
                      * the Species_Type field
                      * (malloced length of NUM_TERMS_MALLOC);
                      */
  double *JacVector; /* Jacobian vector of the source term.
                      * (malloced length of NUM_TERMS_MALLOC);
                      */
};
typedef struct propertyJac PROPERTYJAC_STRUCT;

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/*
 * Prototpyes for functions in mm_propertyJac.c
 */

extern void propertyJac_realloc(PROPERTYJAC_STRUCT **, int);
extern void propertyJac_free(PROPERTYJAC_STRUCT *jac);
extern void propertyJac_destroy(PROPERTYJAC_STRUCT **);
extern void propertyJac_addEnd(PROPERTYJAC_STRUCT *, int, int, int, double, double);
extern void propertyJac_searchadd(PROPERTYJAC_STRUCT *, int, int, int, double, double);
extern int propertyJac_find_species_unk(PROPERTYJAC_STRUCT *);
extern void propertyJac_add1SpEnd(PROPERTYJAC_STRUCT *, MATRL_PROP_STRUCT *, int, double, double);
/******************************************************************************/
/******************************************************************************/
/*
 * This structure packages data regarding the level set initialization method
 * as well as info about geometric forms and their parameters used in
 * initializing the the level set.
 */

struct LS_Surf_Closest_Point {
  double distance;
  double confidence;
  int inflection;
  int elem;
  int elem_side;
  double x[3];
  double xi[3];
};

struct LS_Surf {
  int type;
  void *data;
  struct LS_Surf_Closest_Point *closest_point;
  struct LS_Surf_List *subsurf_list;
  struct LS_Surf *next;
};

struct LS_Surf_List {
  int size;
  struct LS_Surf *start;
  struct LS_Surf *current;
  struct LS_Surf *end;
};

struct LS_Surf_Point_Data {
  double x[3];
  int elem;
  double xi[3];
  int inflection;
};

struct LS_Embedded_BC {
  int bc_input_id;
  struct LS_Embedded_BC *next;
};

struct LS_Surf_Facet_Data {
  int num_points;
  int elem;
  int elem_side;
};

struct LS_Surf_SS_Data {
  int ss_id;
};

struct LS_Surf_NS_Data {
  int ns_id;
  int PosEB_id;
};

struct LS_Surf_Iso_Data {
  int isovar;
  double isoval;
};

struct LS_Surf_Plane_Data {
  double n[3];
  double d;
};

struct LS_Surf_Sphere_Data {
  double center[3];
  double r;
};

struct LS_Surf_Arc_Data {
  double center[3];
  double r;
  double n[3];
  double d;
  double sign;
};

struct LS_Surf_User_Data {
  int Int_Data[5];
  double Real_Data[10];
};

/*
 * Data structure that contains good and useful information about
 * the level set tracking parameters
 */

struct Level_Set_Data {
  int var;
  int Use_Level_Set;
  int Evolution;
  int Contact_Inflection;
  int Isosurface_Subsurf_Type;
  int Init_Method;
  int Num_Var_Init;
  double Length_Scale;
  int adapt;
  double adapt_inner_size;
  double adapt_outer_size;
  double adapt_width;
  int adapt_freq;
  double Control_Width;
  double Renorm_Tolerance;
  int Renorm_Method;
  int Search_Option;
  int Grid_Search_Depth;
  int Integration_Depth;
  int Interface_Output;
  char *output_file;
  int Renorm_Freq;
  int Renorm_Countdown;
  int Force_Initial_Renorm;
  double Mass_Value;
  int Mass_Sign;
  double Contact_Tolerance;
  int Fluid_Solid;
  int Fluid_Sign;
  int Solid_Sign;
  char *sm_object_name;
  char *sm_object_type;
  struct LS_Embedded_BC *embedded_bc;
  struct LS_Surf_List *init_surf_list;
  struct LS_Surf_List *last_surf_list;
  int Elem_Sign;
  int elem_overlap_state;
  int SubElemIntegration;
  int AdaptIntegration;
  int Adaptive_Order;
  int CrossMeshQuadPoints;
  int on_sharp_surf;
  int Extension_Velocity;
  int CalcSurfDependencies;
  int MatrixNum; /* segregated problems */
  int SubcyclesAfterRenorm;
  double Neg_Vol;
  double Pos_Vol;
  double Surface_Area;
  double Initial_LS_Displacement;
  int Ignore_F_deps;
  struct LS_Surf_List *Ext_Vel_Surf_List;
  int Periodic_Planes;
  double Periodic_Plane_Loc[6];
  int Ghost_Integ_Active;
  int Ghost_Integ;
  int PSPP_filter;
  int Sat_Hyst_Renorm_Lockout;
  int ghost_stress;
  int Toure_Penalty;
  int Huygens_Freeze_Nodes;
  int Enable_Div_Term;
  int Semi_Implicit_Integration;
  int YZbeta;
  dbl YZbeta_scale;
};

/*
 * This data structure holds useful level-set related functions and
 * their derivatives with respect to FILL and MESH_DISPLACEMENT{1,2,3}.
 */
struct Level_Set_Interface {
  /* Flag indicating if we're in the interfacial region. */
  int near;

  /* Half the interfacial thickness (alpha = 0.5 * width) */
  double alpha;

  /* Heaviside function; smooth form of: H=0 for F<0, H=1 for F>0 */
  double H;
  double dH;
  double d_H_dF[MDE];
  double d_H_dmesh[DIM][MDE];

  /* Heaviside function as above, but evaluated using FEM basis functions */
  double Hn;
  double Hn_old;
  double gradHn[DIM];
  double gradHn_old[DIM];
  double d_Hn_dF[MDE];
  double d_gradHn_dF[DIM][MDE];
  double d_Hn_dmesh[DIM][MDE];
  double d_gradHn_dmesh[DIM][DIM][MDE];

  /*
   * Delta function: smooth form of: delta(F) = 1 for F=0, =0 for F != 0
   * N.B. This delta has a correction for cases where F is not a pure
   * distance function.
   */
  double delta;
  double d_delta_dF[MDE];
  double d_delta_dmesh[DIM][MDE];
  double delta_max;

  /*
   * Normal vector: typically normal = grad(F); here we use normal =
   * grad(F) / |grad(F)| to be safe.
   */
  double normal[DIM];
  double d_normal_dF[DIM][MDE];
  double d_normal_dmesh[DIM][DIM][MDE];

  /* Magnitude of grad_F. */
  double gfmag;
  double d_gfmag_dF[MDE];
  double d_gfmag_dmesh[DIM][MDE];

  /* Magnitude of grad_F inverse. */
  double gfmaginv;
  double d_gfmaginv_dF[MDE];
  double d_gfmaginv_dmesh[DIM][MDE];
};

struct Search_Grid_Structure {
  struct Element_Indices *ei;  /* parent F E element ei struct */
  int dim;                     /* dimension of search_grid */
  int level;                   /* Level of division of this grid */
  struct Shape_Fcn_Tree *tree; /* Shape function tree that correspondes to this grid */
  int num_verts;               /* number of vertices in this grid */
  double LS_value[8];          /* Values of level set function at search grid vertices */
  struct Search_Grid_Structure **neighbors; /* neighboring search grids of this
                                               grid - Currently unused.  */
  int num_subgrids;
  struct Search_Grid_Structure **subgrids; /* children of the search grid */
};

typedef struct Search_Grid_Structure SGRID;

struct Shape_Fcn_Tree {
  int dim;                          /* dimension of tree struct */
  int level;                        /* Level of division of this tree */
  int num_verts;                    /* number of vertices in this tree */
  double (*xi)[DIM];                /* s,u,t coordinates of grid vertices */
  struct Basis_Functions *bf;       /* point to bfd master basis function list
                                       according to interpolation */
  double num_fcns;                  /* number of shape functions */
  double (*phi)[MDE];               /* value of shape all shape functions at grid nodes */
  int num_gpts;                     /* number of integration points on this tree */
  double (*s)[DIM];                 /* s,u,t coordinates of gauss pts on this tree */
  double *wt;                       /* weight on each integration point */
  int num_subtrees;                 /* number of children */
  struct Shape_Fcn_Tree **subtrees; /* children of the tree, scions if you will */
};

struct Shape_Fcn_Tree_Int {
  int active;
  int ip_total;     /* number of integration points on this tree on this elem */
  double (*s)[DIM]; /* s,u,t coordinates of gauss pts on this tree on this elem */
  double *wt;       /* weight on each integration point on this tree on this elem */
  int *ip_sign;     /* ls->Elem_Sign for integration point */
};

typedef struct Shape_Fcn_Tree NTREE;
typedef struct Shape_Fcn_Tree_Int NTREE_INT;

struct Phase_Function_Jacobian_Info {
  int length;
  double *d_pf_lm; /* Sensitivity of phase function residuals wrt to lagrange
                      multiplier */
  double *d_lm_pf; /* Sensitibity of constraint wrt to phase function unknowns */
  double d_lm_lm;  /* Sensitivity of the constraint wrt to the lagrange mulplier */
};

typedef struct Phase_Function_Jacobian_Info PF_JAC_INFO;

struct Phase_Function_Data {
  int num_phase_funcs;
  int Use_Phase_Field;
  struct Level_Set_Data **ls; /* individual LS data structures for each phase function */
  int Use_Constraint;
  int Constraint_Method;
  double Constraint_Integral;
  double shift[MAX_PHASE_FUNC];
  PF_JAC_INFO *jac_info;
};

/* The maximum number of neighboring element blocks for each shell block */
#define MAX_SHELL_NBRS 6

/* A data structure to help work with blocks of shell elements. */
struct Shell_Block {
  int elemblock_index;               /* The index in the list of element blocks. */
  int elemblock_id;                  /* The number known to the user. */
  int num_nbr_blocks;                /* Number of neighboring element blocks. */
  int *nbr_elem_ids[MAX_SHELL_NBRS]; /* Neighbor element numbers */
  int mn;                            /* This is the material number corresponding to
                                        this shell element block */
};

/* struct for d_Pi */
struct stress_dependence {
  double v[DIM][DIM][DIM][MDE];
  double vd[DIM][DIM][DIM][MDE];
  double X[DIM][DIM][DIM][MDE];
  double C[DIM][DIM][MAX_CONC][MDE];
  double T[DIM][DIM][MDE];
  double nn[DIM][DIM][MDE];
  double F[DIM][DIM][MDE];
  double P[DIM][DIM][MDE];
  double g[DIM][DIM][DIM][DIM][MDE];
  double S[DIM][DIM][MAX_MODES][DIM][DIM][MDE];
  double pf[DIM][DIM][MAX_PHASE_FUNC][MDE];
  double degrade[DIM][DIM][MDE];
  double eddy_nu[DIM][DIM][MDE];
};
typedef struct stress_dependence STRESS_DEPENDENCE_STRUCT;

/* struct for d_q */
struct heat_flux_dependence {
  double T[DIM][MDE];
  double C[DIM][MAX_CONC][MDE];
  double X[DIM][DIM][MDE];
  double F[DIM][MDE];
  double moment[MAX_MOMENTS][DIM][MDE];
};
typedef struct heat_flux_dependence HEAT_FLUX_DEPENDENCE_STRUCT;

/* struct for d_q */
struct acoustic_flux_dependence {
  double P[DIM][MDE];
  double T[DIM][MDE];
  double C[DIM][MAX_CONC][MDE];
  double X[DIM][DIM][MDE];
  double F[DIM][MDE];
};
typedef struct acoustic_flux_dependence ACOUSTIC_FLUX_DEPENDENCE_STRUCT;

/* struct for df */
struct momentum_source_dependence {
  double T[DIM][MDE];                  /* temperature dependence. */
  double X[DIM][DIM][MDE];             /* spatial dependence. */
  double C[DIM][MAX_CONC][MDE];        /* conc dependence. */
  double v[DIM][DIM][MDE];             /* velocity dependence. */
  double F[DIM][MDE];                  /* FILL dependence. */
  double E[DIM][DIM][MDE];             /* electric field dependence */
  double pf[DIM][MAX_PHASE_FUNC][MDE]; /* phase function dependence */
  double ars[DIM][MDE];                /* acoustic energy density. */
};
typedef struct momentum_source_dependence MOMENTUM_SOURCE_DEPENDENCE_STRUCT;

/* struct for d_mu */
struct viscosity_dependence {
  double v[DIM][MDE];             /* velocity dependence. */
  double X[DIM][MDE];             /* mesh dependence. */
  double T[MDE];                  /* temperature dependence. */
  double C[MAX_CONC][MDE];        /* conc dependence. */
  double P[MDE];                  /* pressure dependence. */
  double F[MDE];                  /* FILL dependence. */
  double nn[MDE];                 /* bond concentration dependence */
  double gd;                      /* strain rate dependence */
  double pf[MAX_PHASE_FUNC][MDE]; /* phase function */
  double degrade[MDE];            /* amount of degradation */
  double eddy_nu[MDE];            /* Turbulent viscosity */
};
typedef struct viscosity_dependence VISCOSITY_DEPENDENCE_STRUCT;

/* struct for d_saramito */
struct saramito_coefficient_dependence {
  double s[DIM][DIM]; /* stress dependence. */
  double tau_y;       /* yield stress dependence. */
};
typedef struct saramito_coefficient_dependence SARAMITO_DEPENDENCE_STRUCT;

/* struct for d_dilMu */
struct dilViscosity_dependence {
  double v[DIM][MDE];             /* velocity dependence. */
  double X[DIM][MDE];             /* mesh dependence. */
  double T[MDE];                  /* temperature dependence. */
  double C[MAX_CONC][MDE];        /* conc dependence. */
  double P[MDE];                  /* pressure dependence. */
  double F[MDE];                  /* FILL dependence. */
  double nn[MDE];                 /* bond concentration dependence */
  double gd;                      /* strain rate dependence */
  double pf[MAX_PHASE_FUNC][MDE]; /* phase function */
  double degrade[MDE];            /* amount of degradation */
};
typedef struct dilViscosity_dependence DILVISCOSITY_DEPENDENCE_STRUCT;

/* struct for d_vconv */
struct convection_velocity_dependence {
  double X[DIM][DIM][MDE];      /* mesh dependence. */
  double rs[DIM][DIM][MDE];     /* real solid motion dependence. */
  double v[DIM][DIM][MDE];      /* velocity dependence. */
  double T[DIM][MDE];           /* temperature dependence. */
  double C[DIM][MAX_CONC][MDE]; /* conc dependence. */
};
typedef struct convection_velocity_dependence CONVECTION_VELOCITY_DEPENDENCE_STRUCT;

/* struct for d_vnorm */
struct normal_velocity_dependence {
  double v[DIM][MDE];      /* velocity dependence. */
  double T[MDE];           /* temperature dependence. */
  double C[MAX_CONC][MDE]; /* conc dependence. */
  double V[MDE];           /* voltage dependence. */
  double F[MDE];           /* FILL dependence. */
  double X[DIM][MDE];      /* mesh dependence. */
};
typedef struct normal_velocity_dependence NORMAL_VELOCITY_DEPENDENCE_STRUCT;

/* struct for d_enorm */
struct normal_energy_dependence {
  double v[DIM][MDE];      /* velocity dependence. */
  double T[MDE];           /* temperature dependence. */
  double C[MAX_CONC][MDE]; /* conc dependence. */
  double V[MDE];           /* voltage dependence. */
  double F[MDE];           /* FILL dependence. */
};
typedef struct normal_energy_dependence NORMAL_ENERGY_DEPENDENCE_STRUCT;

/* struct for d_pspg */
struct pspg_dependence {
  double v[DIM][DIM][MDE]; /* velocity dependence. */
  double T[DIM][MDE];      /* temperature dependence. */
  double eddy_nu[DIM][MDE];
  double P[DIM][MDE];           /* pressure dependence. */
  double C[DIM][MAX_CONC][MDE]; /* conc dependence. */
  double X[DIM][DIM][MDE];      /* mesh dependence. */
  double g[DIM][DIM][DIM][MDE];
  double S[DIM][MAX_MODES][DIM][DIM][MDE]; /* stress mode dependence. */
};
typedef struct pspg_dependence PSPG_DEPENDENCE_STRUCT;

/* struct for d_cont_gls */
struct cont_gls_dependence {
  double v[DIM][MDE]; /* velocity dependence */
  double X[DIM][MDE]; /* mesh dependence */
};
typedef struct cont_gls_dependence CONT_GLS_DEPENDENCE_STRUCT;

struct Petrov_Galerkin_Data {
  double h[DIM];
  double h_elem_avg;
  double hsquared[DIM];
  double hh[DIM][DIM];
  double dh_dxnode[DIM][MDE];
  double U_norm;
  double mu_avg;
  double hhv[DIM][DIM];
  double dhv_dxnode[DIM][MDE];
  double rho_avg;
  double v_avg[DIM];
  double dv_dnode[DIM][MDE];
};

typedef struct Petrov_Galerkin_Data PG_DATA;

/*
 * Auxiliaries Variables used in calculating flow rate and average velocity
 * with their sensitivities in lubrication flow.
 * These values are calculated at the current quadrature point.
 */
struct Lubrication_Auxiliaries {
  double q[DIM];             /* Volumetric flow rate per unit width */
  double v_avg[DIM];         /* Average velocity, i.e. q divided by height */
  double gradP_mag;          /* Magnitude of pressure gradient */
  double gradP_tangent[DIM]; /* Tangent vector of the pressure gradient */
  double gradP_normal[DIM];  /* Unit vector perpendicular to the pressure */
  double H;                  /* Lubrication Gap Height */
  double srate;              /* Lubrication Characteristic Shear Rate */
  double mu_star;            /* Lubrication Characteristic Viscosity */

  double dgradP_mag_dP;          /* Pressure gradient magnitude sensitivities w.r.t.
                                    pressure */
  double dgradP_tangent_dP[DIM]; /* Pressure gradient tangent sensitivities
                                    w.r.t. pressure */
  double dgradP_normal_dP[DIM];  /* Pressure gradient normal sensitivities w.r.t.
                                    pressure */

  double dq_dh[DIM][MDE];              /* Flow rate sensitivities w.r.t. height */
  double dq_dh1[DIM][MDE];             /* Flow rate sensitivities w.r.t. height */
  double dq_dh2[DIM][MDE];             /* Flow rate sensitivities w.r.t. height */
  double dq_dp[DIM][MDE];              /* Flow rate sensitivities w.r.t. lubrication pressure */
  double dq_dp1[DIM][MDE];             /* Flow rate sensitivities w.r.t. lubrication pressure */
  double dq_dp2[DIM][MDE];             /* Flow rate sensitivities w.r.t. lubrication pressure */
  double dq_df[DIM][MDE];              /* Flow rate sensitivities w.r.t. level set */
  double dq_dk[DIM][MDE];              /* Flow rate sensitivities w.r.t. curvature */
  double dq_dx[DIM][DIM][MDE];         /* Flow rate sensitivities w.r.t. mesh deformation */
  double dq_dnormal[DIM][DIM][MDE];    /* Flow rate sensitivities w.r.t. shell normal */
  double dq_drs[DIM][DIM][MDE];        /* Flow rate sensitivities w.r.t. real solid
                                          deformation */
  double dq_ddh[DIM][MDE];             /* Flow rate sensitivities w.r.t. heat transport */
  double dq_dc[DIM][MDE];              /* Flow rate sensitivities w.r.t. particles volume
                                          fraction */
  double dq_dconc[DIM][MAX_CONC][MDE]; /* Flow rate sensitivities w.r.t. species concentration */
  double dq_dshear_top[DIM][MDE];      /* Flow rate sensitivities w.r.t. top wall
                                          shear rate */
  double dq_dshear_bot[DIM][MDE];      /* Flow rate sensitivities w.r.t. bottom wall
                                          shear rate */
  double dq_dcross_shear[DIM][MDE];    /* Flow rate sensitivities w.r.t. cross
                                          stream shear stress */
  double dq_dgradp[DIM][DIM][MDE];     /* Flow rate sensitivities w.r.t. pressure gradient */

  double dv_avg_dh[DIM][MDE];  /* Average velocity sensitivities w.r.t. height */
  double dv_avg_dh1[DIM][MDE]; /* Average velocity sensitivities w.r.t. height */
  double dv_avg_dh2[DIM][MDE]; /* Average velocity sensitivities w.r.t. height */
  double dv_avg_dp[DIM][MDE];  /* Average velocity sensitivities w.r.t. lubrication pressure */
  double dv_avg_dp1[DIM][MDE]; /* Average velocity sensitivities w.r.t.
                                  lubrication pressure */
  double dv_avg_dp2[DIM][MDE]; /* Average velocity sensitivities w.r.t.
                                  lubrication pressure */
  double dv_avg_dnormal[DIM][DIM][MDE]; /* Average velocity sensitivities w.r.t. mesh deformation */
  double dv_avg_df[DIM][MDE];           /* Average velocity sensitivities w.r.t. level set */
  double dv_avg_dk[DIM][MDE];           /* Average veloctiy sensitivities w.r.t. curvature */
  double dv_avg_dx[DIM][DIM][MDE];      /* Average velocity sensitivities w.r.t. mesh
                                           deformation */
  double dv_avg_drs[DIM][DIM][MDE];     /* Average velocity sensitivities w.r.t.
                                           real solid deformation*/
  double dv_avg_ddh[DIM][MDE];          /* Average velocity sensitivities w.r.t. heat
                                           transport */
  double dv_avg_dc[DIM][MDE];           /* Average velocity sensitivities w.r.t. particles
                                           volume fraction */
  double dv_avg_dconc[DIM][MAX_CONC]
                     [MDE]; /* Average velocity sensitivities w.r.t. species concentration */
  double dv_avg_dshear_top[DIM][MDE];   /* Average velocity sensitivities w.r.t.
                                           top wall shear rate */
  double dv_avg_dshear_bot[DIM][MDE];   /* Average velocity sensitivities w.r.t.
                                           bottom wall shear rate */
  double dv_avg_dcross_shear[DIM][MDE]; /* Average velocity sensitivities w.r.t.
                                           cross stream shear stress */
  double dv_dgradp[DIM][DIM][MDE]; /* Average velocity sensitivities w.r.t. pressure gradient */
  double dH_dmesh[DIM][MDE];       /* lubrication gap sensitivities w.r.t. mesh */
  double dH_drealsolid[DIM][MDE];  /* lubrication gap sensitivities w.r.t. real
                                      solid */
};

typedef struct Lubrication_Auxiliaries LUBRICATION_AUXILIARIES_STRUCT;

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
#endif /* GOMA_MM_AS_STRUCTS_H */
