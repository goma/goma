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
 *$Id: bc_dirich.c,v 5.3 2010-07-21 16:39:26 hkmoffa Exp $
 */

/* Standard include files */

#include <stdio.h>

/* GOMA include files */

#include "std.h"
#include "stdlib.h"
#include "rf_fem_const.h"
#include "rf_fem.h"

#include "rf_io_const.h"
#include "rf_io_structs.h"
#include "rf_io.h"
#include "rf_mp.h"
#include "rf_solver.h"
#include "el_elm.h"
#include "el_geom.h"

#include "rf_masks.h"
#include "rf_bc_const.h"
#include "rf_solver_const.h"
#include "rf_fill_const.h"
#include "rf_vars_const.h"
#include "mm_mp_const.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"

#include "mm_eh.h"

#define _BC_DIRICH_C
#include "goma.h"

/************************************************************************************/
/************************************************************************************/
/************************************************************************************/

int
put_dirichlet_in_matrix(double x[], const int num_total_nodes)

    /**************************************************************************
     *
     * put_dirichlet_in_matrix():
     *
     * Function which fills the stiffness matrices and the 
     * right-hand side for simple dirichlet conditions. Specifically, if it
     * finds that a dirichlet condition must be applied for a variable, it will
     * set the row of the local Jacobian stiffness matrix to zero with a one
     * the diagonal, and then set the row of the local residual vector to the
     * residual of the Dirichlet condition.
     *
     *  Return
     * --------
     *  always returns 0
     *************************************************************************/
{
  int i, ieqn, ibc;		/* local node number in current element */
  int eqn, var;                 /* squished equation and variable number */
  int var_type;
  int I;			/* processor node number */
  int ldof_eqn;                 /* conversion from node to local dof */
  int offset, lvdesc, idof, matID, j, ledof, matID_dof, found;
  double V_set;
  NODE_INFO_STRUCT *node;
  NODAL_VARS_STRUCT *nv;
  VARIABLE_DESCRIPTION_STRUCT *vd;

  /*
   * Return for numerical jacobian options.
   */
  if (Debug_Flag < 0) return 0;

  for (i = 0; i < ei->num_local_nodes; i++) {
    I = Proc_Elem_Connect[ei->iconnect_ptr + i];
    node = Nodes[I];
    if (node->DBC && I < num_total_nodes) {
      nv = node->Nodal_Vars_Info;
      offset = 0;
      for (lvdesc = 0; lvdesc < nv->Num_Var_Desc; lvdesc++) {
	vd = nv->Var_Desc_List[lvdesc];
	for (idof = 0; idof < vd->Ndof; idof++) { 
	  if (node->DBC[offset] != -1) {
	    ibc = node->DBC[offset];
	    var_type = vd->Variable_Type;
	    matID = vd->MatID;
	    /*
	     * We are tacitly assuming here that all 
	     * lvdofs for the same variable type and subvariable
	     * type at the same node are grouped contiguously
	     * in the local element stiffness matrix lvdof
	     * ordering.
	     */
	    ldof_eqn = ei->ln_to_first_dof[var_type][i];
	    for (j = 0, found = FALSE; j < Dolphin[I][var_type]; j++) {
	      ledof = ei->lvdof_to_ledof[var_type][ldof_eqn];
	      matID_dof = ei->matID_ledof[ledof];
	      if (matID == matID_dof || matID == -1) {
		found = TRUE;
		break;
	      }
	      ldof_eqn++;
	    }
	    if (!found) {
	      EH(-1,"ERROR");
	    }
	    if (var_type == MASS_FRACTION) {
	      eqn = MAX_PROB_VAR + vd->Subvar_Index;
	      var = eqn;
	    } else {
	      eqn = upd->ep[var_type];
	      var = upd->vp[var_type];
	    }

	    if (pd->e[var_type])  /*Big test here.  We are no longer applying 
				   dirichlets from materials that the variables 
				   isn't defined */
	      {
	       	// Add for F_DIODE_BC if((var_type != FILL) )
		// Add for F_DIOD_BC  {

                if(BC_Types[ibc].BC_Name != DX_NOTHING_BC &&
                   BC_Types[ibc].BC_Name != DY_NOTHING_BC &&
                   BC_Types[ibc].BC_Name != DZ_NOTHING_BC )
                   {
		    zero_lec_row(lec->J, eqn, ldof_eqn);
		    if (!(af->Assemble_LSA_Mass_Matrix)) {
                      lec->J[LEC_J_INDEX(eqn,var,ldof_eqn,ldof_eqn)] = DIRICHLET_PENALTY;
		    }
		    if (BC_Types[ibc].BC_relax == -1.0) {
                      lec->R[LEC_R_INDEX(eqn,ldof_eqn)] = 0.0;
		    } else {
		      ieqn  = node->First_Unknown + offset;
		      V_set = BC_Types[ibc].BC_Data_Float[0];
                      lec->R[LEC_R_INDEX(eqn,ldof_eqn)] = DIRICHLET_PENALTY * (x[ieqn] - V_set);
		    }
                   }    /* end of DX_NOTHING test  */
		    // Add for F_DIODE_BC   }
	      }
		
	  } 
	  offset++;
	}
      }
    }
  }
  return(0);
} /* END of put_dirichlet_in_matrix                                           */
/******************************************************************************/
/******************************************************************************/
