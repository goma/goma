/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2014 Sandia Corporation.                                  *
*                                                                         *
* Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,  *
* the U.S. Government retains certain rights in this software.            *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
\************************************************************************/

/* 
 * rf_vars_const.h:
 *
 *     This include file contains information pertinent to the
 *     specification of variables and variable types at a node
 *
 *
 *    $Id: rf_vars_const.h,v 5.1 2007-09-18 18:53:47 prschun Exp $
 */
#ifndef GOMA_RF_VARS_CONST_H_
#define GOMA_RF_VARS_CONST_H_

#include "std.h"
#include "rf_fem_const.h"

/*
 * MAX_MATERIALS_AT_A_POINT:
 *
 * Maximum number of materials in which discontinuous interpolation
 * of variables are used that can meet at a single point. This
 * number is the maximum number of duplicate unknowns due to
 * the discontinuous variable treatment that can exist at a
 * single node.
 */
#ifndef MAX_MATERIALS_AT_A_POINT
#define MAX_MATERIALS_AT_A_POINT 4
#endif

/*
 *  NDOF_MAX = Maximum value of Ndof allowed in the structure below
 */
#ifndef MAX_NDOF
#define MAX_NDOF 4
#endif

/*
 * Variable_Description structure
 *
 *  This contains all of the information to uniquely identify a
 *  variable. Note there is only a small number of these structures
 *  malloced per problem.
 */

struct Variable_Description {
    int List_Index;          /*
			      * This is an index key. Its unique across
			      * all variable_descriptions on a single
			      * processor.
			      */
    int Variable_Type;       /* Global variable type, e.g., velocity*/
    int Base_Variable_Type;  /* variable type without the subindex
			      *	expansion (ie. SPECIES_UNK_23 corresponds
			      * to base type of SPECIES_UNK_0) */
    int Ndof;                /* Number of degrees of freedom for that
			      * variable type
                              * Delineated cases where this value is nonunity:
			      * 1) Discontinuous pressure interpolations
			      *    at centroid node
			      * (can't think of any others)
			      */
    short int Subvar_Index;  /*
			      * subvariable index (if applicable for the
			      * variable)
			      *
			      * NOTE -> In process of being phased out
			      *    -> only applicable for MASS_FRACTION
			      *       variable type.
			      *    Ndof = 1 by definition here too.
			      */
    int MatID;               /*
			      * Material index for the variable described
			      * by this structures
			      *
			      * -> A value of -1 indicates that this is a
			      *    generic variable, applicable for
			      *    multiple materials, and contiguous
			      *    across the interface between these
			      *    materials
			      */
    char *Var_Name[MAX_NDOF];/* Name of the degree of freedom
			      */
    /*
     * The next three variables refer to properties need by the nonlinear
     * solution algorithm to perform global bounds checking and damping
     * algorithms on a per variable_type basis.
     */
    double Upper_Bound[MAX_NDOF]; /* Upper bound for this variable */
    double Lower_Bound[MAX_NDOF]; /* Lower bound for this variable
				   * - violations of bounds causes damping
				   */
    double Delta_Bound[MAX_NDOF]; /* If a variable makes a jump exceeding
				   * this relative multiplicative factor, then
				   * we can assume that the previous Jacobian
				   * was unrepresentative of the new state
				   * of the solution vector. In this case,
				   * we would damp the result normally.
				   */
    double Common_Value[MAX_NDOF];/* This is used to help specify numerical
				   * derivatives where applicable
				   * Numerical derivatives are greatly
				   * aided by the right order of magnitude
				   * in the value of delta.
				   */
    double *double_work;          /* future expansion */
    double *integer_work;         /* future expansion */
};
typedef struct Variable_Description VARIABLE_DESCRIPTION_STRUCT;

/*
 * Declarations for global definitions in rf_vars_defn.h
 */
extern VARIABLE_DESCRIPTION_STRUCT **Var_Info;
extern int Num_Var_Info_Records;

/*
 * Nodal_Vars Structure;
 *
 *     This structure holds the information about what variables are
 *     defined at a particular node.
 *
 *     There will only be a small number of these Nodal_Vars structures
 *     defined in a problem. i.e., there will not be a number equal to the
 *     number of nodes defined in the processor.
 */

struct Nodal_Vars {
    int      list_index;    /* Unique id for this nodal var type
			     * on this processor. You can look up this
			     * entry in Proc_Nodal_Vars_List[] using
			     * this index.
			     */
    int *Var_Type_Index[V_LAST];
                            /* This is an array of integers, representing
                             * the indecese into Node_Vars_List that
			     * contain variable descriptions with
			     * variable types equal to the array index.
                             */
    short int Num_Var_Desc_Per_Type[V_LAST];
                            /* Number of Variable Description structures
                             * for each variable type
                             */
    int Num_Var_Desc;       /* Total number of variable description
			     * structures defined at this node
			     */
    int Num_Unknowns;       /* Total number of unknowns defined at this
			     * node.
			     */
    VARIABLE_DESCRIPTION_STRUCT **Var_Desc_List;
                            /* List of Pointers to Variable Description
                             * structures defined at this node. This is the
			     * prefered way to raster through variables defined
			     * at a node.
			     * Length = Num_Var_Desc
                             */
    int *Nodal_Offset;      /* Vector of Offsets for each variable type
			     * from the list of variables in the solution
			     * vector corresponding to this node.
			     */
};
typedef struct Nodal_Vars NODAL_VARS_STRUCT;

/*
 * Proc_Nodal_Vars_list:
 *   Processor list of unique Nodal_Vars structures
 *   Proc_Nodal_Vars_List[i] points to the ith unique Nodal_Vars_List
 */

extern NODAL_VARS_STRUCT **Nodal_Vars_List;
extern int                 Nodal_Vars_List_Length;

/*
 * vd_packed and nv_packed structures:
 *       These structures are used for communicating variable description
 *       information from owned nodes to ghost nodes in an mp problem. They
 *       are as short as permissible, for this reason.
 */
struct vd_packed {
    short int Variable_Type;
    short int Ndof;
    short int Subvar_Index;
    short int MatID;
};

struct nv_packed {
    short int num;
    struct vd_packed VDP_i;
};

/*
 * prototypes for functions in the rf_vars.c file
 */
extern int find_base_variable_type(int);
extern int variable_description_comparison(VARIABLE_DESCRIPTION_STRUCT *,
				           VARIABLE_DESCRIPTION_STRUCT *);
extern int variable_description_init(VARIABLE_DESCRIPTION_STRUCT **);
extern int variable_description_default(VARIABLE_DESCRIPTION_STRUCT *);
extern int variable_description_alloc(VARIABLE_DESCRIPTION_STRUCT *,
                                      const int, const int,  const int,
				      const int);
extern void variable_description_free(VARIABLE_DESCRIPTION_STRUCT *);
extern void variable_description_destroy(VARIABLE_DESCRIPTION_STRUCT **);

extern VARIABLE_DESCRIPTION_STRUCT *
variable_description_create(const int, const int, const int, const int);

extern VARIABLE_DESCRIPTION_STRUCT *
find_or_create_vd(const int, const int, const int, const int, const int);

extern VARIABLE_DESCRIPTION_STRUCT *
Variable_Description_Match(VARIABLE_DESCRIPTION_STRUCT *);

extern VARIABLE_DESCRIPTION_STRUCT *get_vd_ptr(int, int, int);
extern void vdesc_augment(void);

extern void assign_species_desc_suffix(const int, char *);
extern void assign_species_name(const int n,
				MATRL_PROP_STRUCT *mp_local,
				char *, char *, int);
extern void assign_var_name(const int, const int, MATRL_PROP_STRUCT *,
			    char *, char *, int);
#endif
