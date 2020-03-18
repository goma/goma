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
 *$Id: rf_node_const.h,v 5.2 2007-09-18 18:53:47 prschun Exp $
 */

#ifndef GOMA_NODE_CONST_H_
#define GOMA_NODE_CONST_H_

#include "exo_struct.h"
#include "mm_as_structs.h"
#include "rf_vars_const.h"

/*****************************************************************************/
/*                       STRUCTURE DEFINITIONS                               */
/*****************************************************************************/

/*
 * BCCondMask: masks for specifying certain boundary condition flags on a node.
 *
 * This is a bit-field structure defined for every BC.  If the bit is set, this
 * implies that the equation corresponding to the unknown at that node is
 * flagged.  The actual meaning of the flag depends upon the allocated arrays
 * of this structure.  Current identities for the bit fields are:
 *
 *       DBCU1     : U velocity
 *       DBCU2     : V velocity
 *       DBCU3     : W velocity
 *       DBCT      : Energy equation
 *       DBCY      : All gas-species equations
 *       DBCY_SPEC : Some, but not all, gas-species equations.
 *       DBCS      : All surface site fraction equations
 *       DBCP      : Pressure equation
 *       DBCK      : k-transport equation
 *       DBCE      : e-transport equation
 *       DBCTV     : TV-transport equation
 *       DBCES     : Electrostatic equation
 *       DBCMS     : Magnetostatic equation
 */

union BCCondMask {
  unsigned int BCCondInt;
  struct {
    unsigned int DBCU1     : 1;
    unsigned int DBCU2     : 1;
    unsigned int DBCU3     : 1;
    unsigned int DBCUN     : 1;
    unsigned int DBCUT1    : 1;
    unsigned int DBCUT2    : 1;
    unsigned int DBCU1_STEF_VEL : 1;   /* Set to 1 if the U1, U2, or U3  */
    unsigned int DBCU2_STEF_VEL : 1;   /* Dirichlet BCs enforced at node */
    unsigned int DBCU3_STEF_VEL : 1;   /* are surface_chemkin_bc.        */
    unsigned int DBCT      : 1;
    unsigned int DBCY      : 1;
    unsigned int DBCY_SPEC : 1;        /* Some, but not all, species are
                                          flagged */
    unsigned int DBCS      : 1;
    unsigned int DBCP      : 1;
    unsigned int DBCK      : 1;        /* 1-bit for k-equation EBC's     */
    unsigned int DBCE      : 1;        /* 1-bit for e-equation EBC's     */
    unsigned int DBCTV     : 1;        /* 1-bit for TV-equation EBC's    */

    unsigned int DBCES     : 1;        /* 1-bit for electrostatic BC's   */
    unsigned int DBCMS     : 1;        /* 1-bit for magnetostatic BC's   */
  } BCCondBits;
};

/****   Data type used for boundary condition masks.  ****/

typedef unsigned short MASK_TYPE;

/****   Structure for node-type information   ****/
struct Node_Type {                /* Flags indicating the type of the node.*/
  unsigned int Internal :1;       /* TRUE if node is an INTERNAL node.     */
  unsigned int Border   :1;       /* TRUE if node is a BORDER node.        */
  unsigned int External :1;       /* TRUE if node is an EXTERNAL node.     */
  unsigned int Owned    :1;       /* TRUE if node is OWNED (i.e., INTERNAL */
};                                 /* or BORDER.)                           */


struct UMI_list {
  int *List;
  int Length;
};
typedef struct UMI_list UMI_LIST_STRUCT;

/*
 *  This structure is precalculated from the boundary conditions
 *  and applied to the node level. contributions from
 *  var_type, Var_Type, in material MatID_Volume are placed in
 *  the row corresponding to var_type, Var_Type, in material,
 *  MatID_Row_Redirect. Currently DVI_VDDSIG boundary conditions
 *  do this
 */
struct Volume_Contrib_Row_Redirection {
    int Var_Type;
    int MatID_Volume;
    int MatID_Row_Redirect;
};
typedef struct Volume_Contrib_Row_Redirection VCRR_STRUCT;

/*
 *   Structure for holding the results of precalculated
 *   indexing issues and temporary calculations
 *   that occur during the residual and Jacobian fills
 */
struct Nodal_Resid_Wksp {
    VCRR_STRUCT **Vcrr;
    void *SurfaceCalcStorage;
    /* void *VolumeCalcStorage; */
};
typedef struct Nodal_Resid_Wksp NODAL_RESID_WKSP_STRUCT;

/*
 *       Structure for node-based information
 */

struct Node_Info {
    int Proc_Node_Num;              /* Node's (local) number on processor for 
                                       a given mesh. This field
                                       indicates the order in which nodes are
                                       processed.  With the dynamic data 
                                       structure, this field may change
                                       after the mesh is modified.
                                       Array index is mesh->Mesh_ID.         */
    int *First_Unknown;             /* Location in soln vector of node's
                                       first unknown
				       HKM -> Note there is no mesh
				       dependence for this number
				       in the goma version                   */
    struct Node_Type Type;          /* Flags indicating the type of the node.*/


    UMI_LIST_STRUCT Mat_List;       /* This is a list of material indecises
				     * present at this node */
    /* UMI_LIST_STRUCT Element_List; */ 

    NODAL_VARS_STRUCT **Nodal_Vars_Info;
                                    /* Pointer to the Nodal_Vars struct for
                                     * this node.  The Nodal_Vars struct has
				     * info on which variables are active at
				     * this node and how many sub-variables 
				     */
    int Proc;                        /* Processor that owns the node.  We'll
				     * try to find a way around this field,
				     * as it is needed only for External
				     * nodes. 
				     */
    int Global_Node_Num;            /* Node's Global node number.
				     *      (was GNodes[])                   */
    NODAL_RESID_WKSP_STRUCT *Resid_Wksp;
                                    /* Pointer to hold the temporary results
				     * computed in the residual and Jacobian
				     * fill operations.                     
				     */

    /* Bit patterns denoting boundary condition information */
    unsigned int DBCA     : 1;
    unsigned int DBCA_FIX : 1;
    unsigned int DBST     : 1;
    unsigned int DBSTS    : 1;
    unsigned int DBSH_SLOPE_X : 1;
    unsigned int DBSH_SLOPE_Y : 1;
  /* bit mask for SHEET_ENDSLOPE special condition */
	unsigned int DBSES : 1;

    /* bit masks for distinguishing  mesh movement conditions */
    unsigned int DBCDD1 : 1;
    unsigned int DBCDD2 : 1;
    unsigned int DBCDD3 : 1;
    /*
     * EDGE is a bit mask that denotes the fact that the current node
     * is part of side set. Note, it doesn't matter whether or
     * not that side set is used in an actual boundary condition
     * or not
     */
    unsigned int EDGE : 1;
    /*
     * DISC_BNDRY is a bit mask to denote whether this node is part
     * of a side set where discontinuous boundary conditions of some
     * sort are being applied.
     */
    unsigned int DISC_BNDRY : 1;
    /*
     * DBC[i]:
     *      If a Dirichlet condition for the ith unknown in the solution
     *  vector at this node is assigned, the DBC[i] will contain the
     *  index of the boundary condition pertaining to this Dirichlet
     *  condition. Note, multiple species unknowns do count in the numbering
     *  of the unknowns at the interface.
     *  Length = Number of unknowns at this node
     */
    short int **DBC;
};
typedef struct Node_Info NODE_INFO_STRUCT;

/* Declarations for global unknowns defined in rf_node.h */

/*
 *
 *  Global Pointer to the Vector of Pointers to Node Structures
 *  for all of the nodes on the current processor
 *  Index of the array is the processor (local) node number.
 *
 *  Nodes[node_num] points to a Node_Info structure, with information
 *  unique to each node.
 *
 */
extern NODE_INFO_STRUCT **Nodes;

/*
 * Prototypes for functions defined in rf_node.c
 */
extern int bin_search_max(const int [], const int, const int);
extern void add_to_umi_int_list(UMI_LIST_STRUCT *, int);
extern void init_nodes(Exo_DB *, Dpi *);
extern void free_nodes(void);
extern void make_node_to_elem_matrl_map(Exo_DB *);
extern int first_matID_at_node(const int);
extern int index_umi_list(const UMI_LIST_STRUCT *, const int);
extern int node_matrl_index(const int, const int);
extern void free_umi_list(UMI_LIST_STRUCT *);
extern void node_info_tmp_free (void);
extern NODAL_RESID_WKSP_STRUCT *nodal_resid_wksp_alloc (void);
extern void nodal_resid_wksp_destroy(NODAL_RESID_WKSP_STRUCT **);
extern void nodal_resid_wksp_free(NODAL_RESID_WKSP_STRUCT *);
extern void vcrr_add(VCRR_STRUCT ***, const int, const int,
		     const int);
extern int vcrr_lookup(VCRR_STRUCT **, const int, const int);
extern void nullify_dirichlet_bcs(void);

/*
 * Prototypes for the rf_node_vars.c file
 */
extern int dof_lnode_var_type(const int, const int, const int,
			      const int, PROBLEM_DESCRIPTION_STRUCT *, const int);
extern int num_varType_at_node(const int, const int);
extern int get_nv_ndofs(NODAL_VARS_STRUCT *, const int);
extern int get_nv_ndofs_modMF(NODAL_VARS_STRUCT *, const int);
extern int get_nv_index_varType(NODAL_VARS_STRUCT *, const int,
		                const int, const int);

extern int get_nv_offset_idof(NODAL_VARS_STRUCT *, const int, 
			      const int, int, VARIABLE_DESCRIPTION_STRUCT **);

extern int get_nodal_unknown_offset(NODAL_VARS_STRUCT *, const int,
				    const int, const int,
				    VARIABLE_DESCRIPTION_STRUCT **);
extern void get_nv_vd_from_offset(NODAL_VARS_STRUCT *, int,
	 	                  VARIABLE_DESCRIPTION_STRUCT **, int *);
extern int nodal_vars_comparison(NODAL_VARS_STRUCT *,
				 NODAL_VARS_STRUCT *);
extern int MatchAdd_Node_Vars_List(NODAL_VARS_STRUCT *,
			           NODAL_VARS_STRUCT **);
extern void add_var_to_nv_struct(VARIABLE_DESCRIPTION_STRUCT *,
		                 NODAL_VARS_STRUCT *);
extern NODAL_VARS_STRUCT *copy_and_rearrange_nv(NODAL_VARS_STRUCT *);
extern void nodal_vars_free(NODAL_VARS_STRUCT *);
extern void nodal_vars_destroy(NODAL_VARS_STRUCT **);
extern int find_var_type_index_in_nv(const int, NODAL_VARS_STRUCT *);
int find_vd_index_in_nv(const int, const int, NODAL_VARS_STRUCT *);
extern void pack_nv (struct nv_packed *, NODAL_VARS_STRUCT *);
extern NODAL_VARS_STRUCT *unpack_nv (struct nv_packed *, const int);
#endif
