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
 *$Id: mm_interface.h,v 5.2 2007-09-18 18:53:45 prschun Exp $
 */

#ifndef GOMA_MM_INTERFACE_H_
#define GOMA_MM_INTERFACE_H_

#include "el_elm.h"
#include "mm_fill_jac.h"
#include "rf_bc_const.h"
#include "rf_vars_const.h"
#include "std.h"
/*
 *  Interface_Source Structure
 * ------------------------------
 *
 *        This structure is used to communicate with routines that
 *  calculate and track complicated interfacial source terms. The
 *  variables that each source term needs are model specific.
 *  However, several principles are generic. 
 *        The source term need only be calculated once for each 
 *  quadrature point. The results of the calculation may be used
 *  on one side of the interface after it has been calculated and
 *  used on the other side of the interface.
 *        The contents of the state variables are carried through 
 *  from one nonlinear iteration to the next and from one time
 *  step to the next. Thus, they may contain state information for
 *  the interface such as the total excess mass and excess species
 *  arial concentrations at the interface. They can also evolve
 *  with time.
 *        Each Interface_Source structure should only be set up
 *  once for each quadrature point on a surface. At that setup time,
 *  the Var_List vector is set up and thereafter not changed or
 *  freed.
 */

struct Interface_Source {
    int Num_Terms;       /* Number of terms in the source term vector
			  * below.
			  */
    int IS_Type;         /* Type of interface storage - uniquely 
			  * identifies what is storred here
			  */
    int SpeciesVT;       /* Value of the species var type for
			  * the source term and dependent variable assumed
			  * in the Fields below.
			  */ 
    VARIABLE_DESCRIPTION_STRUCT **Var_List;
			 /* Vector of pointers to variable description
			  * structures for each of the source terms
			  * (malloced length of Num_Terms);
			  */ 
    int *idof;           /* Vector of idof's. This is usually equal to
			  * zero, except for Variable descriptions
			  * that pertain to more than one dof. 
			  * Also doubles for kspec for vintage MASS_FRACTION
			  * variable types
			  * (malloced length of Num_Terms);
			  */
    double *Var_Value;   /* Current value of the variable corresponding
			  * to the entry. Species var type agrees with
			  * the Species_Type field
			  */
    double *SourceTerm;  /* Vector of the current value of the source term
			  * The vector is contiguous wrt variable types 
			  * described above
			  * (malloced length of Num_Terms);
			  */
    double **JacMatrix;  /* Jacobian matrix of the source term. It is
			  * assumed to be a dense interaction.
			  * (malloced length of [Num_Terms][Num_Terms] with
			  *  standard row-is-the-inner-loop numbering)
			  *  Jac_ij = i + Num_Terms*j = JacMatrix[i][j]
			  *  the first index is the source term, while the
			  *  second index is the dependent variable
			  */
    void *StateInterface; /* Pointer to memory associated with describing
			   * the time dependent state of the interface
			   * (pointer type and length is specific to the
			   *  model type);
			   */
    void *StateInterfaceOld; /* Pointer to memory associated with describing
			      * the time dependent state of the interface
			      * at the previous time step.		
			      * (pointer type and length is specific to the
			      *  model type);
			      */
    int *Processed;        /* This is true if the interfacial source term
			   * has been calculated for the current nonlinear
			   * calculation at the current time step
			   * It is false if not.
			   */
    int Do_Jac;           /* If true, the Jacobian terms are malloced
			   * zeroed, and calculated along with the 
			   * source term
			   */
    int BC_ID;		  /* sideset ID that identifies the interface*/
};
typedef struct Interface_Source INTERFACE_SOURCE_STRUCT;

struct Interface_Storage {

    struct Interface_Source  *Chemkin_Source;
    struct Interface_Source  *Misc_Source;


};

/*
 * Prototypes for functions in mm_interface.c
 */
extern INTERFACE_SOURCE_STRUCT *interface_source_alloc(int, int, int);
extern void interface_source_free(INTERFACE_SOURCE_STRUCT *);
extern void interface_source_destroy(INTERFACE_SOURCE_STRUCT **);
extern void interface_source_zero(INTERFACE_SOURCE_STRUCT *is);
extern void is_change1_speciesVT(INTERFACE_SOURCE_STRUCT *, const int, 
				 const int, MATRL_PROP_STRUCT *, 
				 const int, const double, const int);
extern void is_change1_lastspecies(INTERFACE_SOURCE_STRUCT *, const int, 
				  const int, MATRL_PROP_STRUCT *, const int);
extern int match_interface_source_string(char *);

/*
 * Prototypes for functions in mm_fill_interface.c
 */
extern INTERFACE_SOURCE_STRUCT *
is_masstemp_create(MATRL_PROP_STRUCT *, MATRL_PROP_STRUCT *, int);

extern double 
raoults_law_prxn(JACOBIAN_VAR_DESC_STRUCT *, BOUNDARY_CONDITION_STRUCT *,
		 int, ELEM_SIDE_BC_STRUCT *, double [MAX_PDIM],
		 const double, const double, const double, const int);
extern void
source_vle_prxn(INTERFACE_SOURCE_STRUCT *, BOUNDARY_CONDITION_STRUCT *,
		MATRL_PROP_STRUCT *, MATRL_PROP_STRUCT *,  int, const int);

extern void 
source_is_equil_prxn(INTERFACE_SOURCE_STRUCT *, BOUNDARY_CONDITION_STRUCT *, 
		     MATRL_PROP_STRUCT *, MATRL_PROP_STRUCT *,  int, const int);
extern double
is_equil_prxn(JACOBIAN_VAR_DESC_STRUCT *,
	      BOUNDARY_CONDITION_STRUCT *, int, ELEM_SIDE_BC_STRUCT *,
	      double [MAX_PDIM], const double, const double, const double,
              const int);

extern void 
is_masstemp_fillin(INTERFACE_SOURCE_STRUCT *, MATRL_PROP_STRUCT *,
		   MATRL_PROP_STRUCT *, int, double, const int);
#endif
