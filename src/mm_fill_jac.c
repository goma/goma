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
 *$Id: mm_fill_jac.c,v 5.2 2008-10-02 15:36:28 hkmoffa Exp $
 */

#include "mm_fill_jac.h"

#include <stdio.h>
#include <string.h>

#include "std.h"
#include "rf_allo.h"

/************************************************************************/
/************************************************************************/
/************************************************************************/

void
jacobianVD_realloc(JACOBIAN_VAR_DESC_STRUCT **jacVD_hdl, int num_lvdesc,
		   int num_lvdof)

   /********************************************************************
    *
    * interface_source_alloc():
    *
    *    Allocate and zero a jacobian variable description structure
    *    We malloc the underlying structures at this time as well.
    *
    * Input
    * ------
    *  num_lvdesc -> number of local variable description entries 
    *                (i.e., the number of dependencies on variables
    *                 in the structure). Note, we do not include the
    *                 dependencies over basis functions at this level.
    *                 These are added in at the time that the contents
    *                 of this structure are added into the local
    *                 element stiffness matrix.
    *  num_lvdof -> Number of local variable degress of freedom 
    *               dependencies. Note, the mesh dependence variables
    *               can not be handled via a num_lvdesc approach, 
    *               because the basis functions themselves depend
    *               upon the mesh variables. Thus, we need this added
    *               functionality for them. And, the only variable
    *               types that need to be handled this way are the 
    *               mesh position variable types.
    *
    * Return:
    *   This function returns the pointer to the malloced and initialized
    *   structure.
    *********************************************************************/
{
  JACOBIAN_VAR_DESC_STRUCT *jacVD;
  if (*jacVD_hdl == NULL) {
    *jacVD_hdl = alloc_struct_1(JACOBIAN_VAR_DESC_STRUCT, 1);
    /*
     * Change the initial function value to some odd number so that
     * zero isn't matched
     */
    (*jacVD_hdl)->Func_Value = -245.23;
  }
  jacVD = *jacVD_hdl;
  if (jacVD->NUM_LVDESC_MALLOC < num_lvdesc) {
    realloc_int_1(&(jacVD->Lvdesc_Index), num_lvdesc,
		  jacVD->NUM_LVDESC_MALLOC);
    realloc_dbl_1(&(jacVD->JacCol), num_lvdesc,
		  jacVD->NUM_LVDESC_MALLOC);
    jacVD->NUM_LVDESC_MALLOC = num_lvdesc;
  }
  if (jacVD->NUM_LVDOF_MALLOC < num_lvdof) {
    realloc_int_1(&(jacVD->Lvdof_var_type), num_lvdof,
		  jacVD->NUM_LVDOF_MALLOC);
    realloc_int_1(&(jacVD->Lvdof_lvdof), num_lvdof,
		  jacVD->NUM_LVDOF_MALLOC);
    realloc_dbl_1(&(jacVD->Jac_lvdof), num_lvdof,
		  jacVD->NUM_LVDOF_MALLOC);
   jacVD->NUM_LVDOF_MALLOC = num_lvdof;
  }
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void
jacobianVD_free(JACOBIAN_VAR_DESC_STRUCT *jacVD)

    /********************************************************************
     *
     * jacobianVD_free():
     *
     *    This routine frees all of the underlying memory associated
     *  with a Jacobian VD struct.
     ********************************************************************/
{
  safer_free((void **) &(jacVD->Lvdesc_Index));
  safer_free((void **) &(jacVD->JacCol));
  safer_free((void **) &(jacVD->Lvdof_var_type));
  safer_free((void **) &(jacVD->Lvdof_lvdof));
  safer_free((void **) &(jacVD->Jac_lvdof));
  jacVD->Num_lvdesc = 0;
  jacVD->NUM_LVDESC_MALLOC = 0;
  jacVD->Num_lvdof = 0;
  jacVD->NUM_LVDOF_MALLOC = 0;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void
jacobianVD_destroy(JACOBIAN_VAR_DESC_STRUCT **jacVD_hdl)

    /********************************************************************
     *
     *jacobianVD_destroy():
     *
     *    This routine frees all of the underlying memory associated
     *  with Jacobian Variable Description structure. It then goes on to
     *  free the structure itself.
     *********************************************************************/
{
  jacobianVD_free(*jacVD_hdl);
  safer_free((void **) jacVD_hdl);
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void
jacobianVD_zero(JACOBIAN_VAR_DESC_STRUCT *jacVD)

   /********************************************************************
    *
    * jacobianVD_zero():
    *
    *    Zero the Jacobian and function value in the structure. 
    * Other information in the structure isn't touched.
    *********************************************************************/
{
  jacVD->Func_Value = -234.45;
  (void) memset((void *) jacVD->JacCol, 0,
		sizeof(double) * jacVD->Num_lvdesc);
  (void) memset((void *) jacVD->Jac_lvdof, 0,
		sizeof(double) * jacVD->Num_lvdof);
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void
jacobianVD_addNewEntry(JACOBIAN_VAR_DESC_STRUCT *jacVD, int index_lvdesc,
		       double entry)

   /********************************************************************
    *
    * jacobianVD_addNewentry():
    *
    *   Add an entry to the Jacobian_Var_Desc structure. Don't bother to
    *   check whether the entry existed previously.
    *********************************************************************/
{
  int num = jacVD->Num_lvdesc;
  if (num >= jacVD->NUM_LVDESC_MALLOC) {
    jacobianVD_realloc(&jacVD, num + 1, 0);
  }
  jacVD->Lvdesc_Index[num] = index_lvdesc;
  jacVD->JacCol[num] = entry;
  jacVD->Num_lvdesc++;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void
jacobianLVDOF_addNewEntry(JACOBIAN_VAR_DESC_STRUCT *jacVD,
			  int var_type, int lvdof, double entry)

   /********************************************************************
    *
    * jacobianLVDOF_addNewEntry():
    *
    *   Add an entry to the Jacobian_Var_Desc structure. Don't bother to
    *   check whether the entry existed previously.
    *   Add the entry as a direct variable type , lvdof entry.
    *********************************************************************/
{
  int num = jacVD->Num_lvdof;
  if (num >= jacVD->NUM_LVDOF_MALLOC) {
    jacobianVD_realloc(&jacVD, 0, num + 1);
  }
  jacVD->Lvdof_var_type[num] = var_type;
  jacVD->Lvdof_lvdof[num] = lvdof;
  jacVD->Jac_lvdof[num] = entry;
  jacVD->Num_lvdof++;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void
jacobianVD_addEntry(JACOBIAN_VAR_DESC_STRUCT *jacVD, int index_lvdesc,
		    double entry)

   /********************************************************************
    *
    * jacobianVD_addentry():
    *
    *   Add a contribution to the Jacobian_Var_Desc structure.
    *   Check whether the entry already exists. If it does, add the
    *   entry to the existing entry and return.
    *********************************************************************/
{
  int i, num, *iptr;
  if (DOUBLE_ZERO(entry)) return;
  iptr = jacVD->Lvdesc_Index;
  num = jacVD->Num_lvdesc;
  for (i = 0; i < num; i++, iptr++) {
    if (index_lvdesc == *iptr) {
       jacVD->JacCol[i] += entry;
       return;
    }
  }
  if (num >= jacVD->NUM_LVDESC_MALLOC) {
    jacobianVD_realloc(&jacVD, num + 1, 0);
  }
  jacVD->Lvdesc_Index[num] = index_lvdesc;
  jacVD->JacCol[num] = entry;
  jacVD->Num_lvdesc++;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/
