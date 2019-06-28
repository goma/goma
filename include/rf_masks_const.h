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
 
#ifndef GOMA_RF_MASKS_CONST_H_
#define GOMA_RF_MASKS_CONST_H_


/*
 * Maximum number of materials in which discontinuous interpolation
 * of variables are used that can meet at a single point. This number is
 * the maximum number of duplicate unknowns due to the discontinuous variable
 * treatment that can exist at a single node.
 */
#define MAX_MATERIALS_AT_A_POINT 4

/*
 * Node_Var_Struct: structure holding specific variable information with the
 * following members: variable type (defined constant), vector size (e.g., 3D
 * fluid-flow => 3 velocity components), ordering at node and ghost variable
 * flag (on REALM boundaries for discontinous variables).
 */

struct Node_Var_Struct {
  int           Variable_Type;     /* Global variable type, e.g., velocity   */
  unsigned short Num_Subvar;        /* Number of sub variables, e.g., 3-D
                                      velocity => 3 subvariables             */
    /* HKM -> Goma has a specific variable type for each velocity. */
    /*        Thus the only variable that will have                */
    /*        Num_Subvar > 1 is the mass fraction vector           */

  unsigned short Num_Disc_Var_Copies;
  short int     Position;          /* Position of the first subvariable in the
                                      list of valid variables at the node    */
  int           Mat_ID[MAX_MATERIALS_AT_A_POINT];

};

/*
 * Var_Struct: structure indicating nodal variable existence, ordering and
 * location information, etc.  If the value of a particular variable is NULL,
 * this implies that the variables associated with that value are not valid.
 */

struct Var_Struct {

  /* NOTE: All these pointers are initialized to NULL => if you add new ones,
     you must initialize them also.  Note also that the order and number of
     these variables MUST match those defined in rf_fem_const.h */

  union {
    struct {
      NODE_VAR_STRUCT *Velocity;           /* Fluid Velocity */
      NODE_VAR_STRUCT *Pressure;           /* Pressure */
      NODE_VAR_STRUCT *Temperature;        /* Temperature */
      NODE_VAR_STRUCT *Mass_Fraction;      /* Species Mass Fraction */
      NODE_VAR_STRUCT *Surface;
      NODE_VAR_STRUCT *k;                  /* Turb. Kinetic Energy */
      NODE_VAR_STRUCT *e;                  /* Turb. Dissipation Rate */
      NODE_VAR_STRUCT *TV;                 /* Turb. Kinematic Viscosity */
      NODE_VAR_STRUCT *Mesh_Displacement;
      NODE_VAR_STRUCT *Site_Fraction;
      NODE_VAR_STRUCT *Site_Density;
      NODE_VAR_STRUCT *Bulk_Mole_Fraction;
      NODE_VAR_STRUCT *ES_Potential;       /* Electrostatic */
      NODE_VAR_STRUCT *MS_Potential;       /* Magnetostatic */
      NODE_VAR_STRUCT *E_Field;
      NODE_VAR_STRUCT *B_Field;
      NODE_VAR_STRUCT *Number_Density;     /* Plasma species number density */
    } Node_Vars;

    NODE_VAR_STRUCT *Node_Vars_Array[MAX_VARIABLE_TYPES];
                                           /* Array of pointers to
                                              Node_Var_Structs */

  } Var_Union;

  short Var_Offsets[MAX_VARIABLE_TYPES]; /* Integer "indexable" offsets for
                                            these nodes.  Note this information
                                            is essentially duplicate of each
                                            variable structure's position
                                            argument but is integer indexable
                                            allowing it to be used in loop
                                            constructs. */

  unsigned short          Num_Vars;      /* Number of distinct variables at
                                            this location.  DOES NOT ACCOUNT
                                            FOR VECTOR ENTRIES (e.g. Velocity
                                            counts as 1) */
  unsigned short          Num_Unknowns;  /* Number of unknowns at this location
                                            INCLUDING VECTOR ENTRIES
                                            (e.g. Velocity in 3D counts as 3)
                                            */
};




/*
 * Externs
 */

extern struct Var_Struct *Var_Info;           /* rf_mask.h */
extern unsigned short   Num_Var_Tags;       /* rf_mask.h */

#endif
