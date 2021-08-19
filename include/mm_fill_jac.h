/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2021 Sandia Corporation.                                  *
*                                                                         *
* Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,  *
* the U.S. Government retains certain rights in this software.            *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
\************************************************************************/
 
/*
 *$Id: mm_fill_jac.h,v 5.1 2007-09-18 18:53:43 prschun Exp $
 */

#ifndef GOMA_MM_FILL_JAC_H_
#define GOMA_MM_FILL_JAC_H_

/*
 *  Jacobian_Var_Desc:
 *
 *     This structure is used to store the dependence of a single function
 *  value on the independent variables in the problem. It is meant to be
 *  memory and speed efficient.
 *
 */

struct Jacobian_Var_Desc {
    double Func_Value;   /* This is the value of the entry, whose
			  * Jacobian values are storred in the
			  * rest of the structure.
			  * (optional, but can be useful, if this 
			  * structure is storred using QP_STORAGE
			  * structs and used later. Or, it can
			  * be used for internal consistency checking.
			  */
    /*
     * First way of storring jacobian entries
     *   -> we don't expand the dependence of the dependent variable
     *      on its underlying basis function interpolation. Instead
     *      we assume a "standard" basis function interpolation
     *      form here.
     */
    int NUM_LVDESC_MALLOC;
                         /* This is the number of var types used in
			  * the malloc.
			  */
    int Num_lvdesc;      /* Number of terms in the Jacobian term vector
			  * below corresponding to different local
			  * variable descriptions.
			  *   Note, we do not include the
			  * dependencies over basis functions at this level.
			  * These are added in at the time that the contents
			  * of this structure are added into the local
			  * element stiffness matrix.
			  */
    int *Lvdesc_Index;   /* Vector of indexes into the local variable
			  * descriptions for each entry below.
			  */
    double *JacCol;      /* Jacobian entries for each local variable
			  * description index.
			  */
    /*
     * Second way of storring jacobian entries
     * -> Note we need this because the dependence of surface normals
     *    on mesh positions are based on nodal coordinates that MUST
     *    be expressed at the lvdof level. All other information
     *    is best handled at the local variable description
     *    level described above.
     */
    int NUM_LVDOF_MALLOC;/* This is the number of direct lvdof entries
			  * used in the malloc
			  */
    int Num_lvdof;       /*  Number of local variable degress of freedom 
			  *  dependencies. Note, the mesh dependence variables
			  *  can not be handled via a num_lvdesc approach, 
			  *  because the basis functions themselves depend
			  *  upon the mesh variables. Thus, we need this added
			  *  functionality for them. And, the only variable
			  *  types that need to be handled this way are the 
			  *  mesh position variable types.
			  */
    int *Lvdof_var_type; /* Variable type for the direct lvdof entry
			  */
    int *Lvdof_lvdof;    /* Local variable degree of freedom index
			  * for the direct lvdof entry
			  */
    double *Jac_lvdof;   /* Entry for the direct lvdof entry
			  */
};
typedef struct Jacobian_Var_Desc JACOBIAN_VAR_DESC_STRUCT;

/*
 * Prototype definitions for functions in mm_fill_jac.c
 */
extern void jacobianVD_realloc(JACOBIAN_VAR_DESC_STRUCT **, int, int);
extern void jacobianVD_free(JACOBIAN_VAR_DESC_STRUCT *);
extern void jacobianVD_destroy(JACOBIAN_VAR_DESC_STRUCT **);
extern void jacobianVD_zero(JACOBIAN_VAR_DESC_STRUCT *);
extern void jacobianVD_addNewEntry(JACOBIAN_VAR_DESC_STRUCT *, int,
				   double);
extern void jacobianLVDOF_addNewEntry(JACOBIAN_VAR_DESC_STRUCT *,
				      int, int, double);
extern void jacobianVD_addEntry(JACOBIAN_VAR_DESC_STRUCT *, int,
				double);
#endif
