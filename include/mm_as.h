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
 *$Id: mm_as.h,v 5.3 2009-07-30 19:52:08 hkmoffa Exp $
 */

/* 	For definitions of these structures, see the entries in
 * 	files "mm_as_const.h" and "mm_as_structs.h"
 */

#ifndef _MM_AS_H
#define _MM_AS_H

extern struct Element_Indices			*ei;
extern struct Element_Indices		       **eiRelated;
extern struct Element_Stiffness_Pointers	*esp;
extern struct Element_Quality_Metrics           *eqm;
extern struct Element_Variable_Pointers	        *esp_old, *esp_dot, *esp_dbl_dot, *evp ;
extern UPD_STRUCT                               *upd;

/*
 *   pd_glob
 *       pd_glob[mn] is a pointer to the Problem_Description structure for the material
 *       mn. mn is the material number index, extending from 0 to the number of materials-1
 *       defined in the problem. mn is distinct from the material index, because it is
 *       required to be contiguous. 
 *
 *    pd is the pointer into pd_glob[] pertinent to the current material. It is redefined
 *       where needed.
 */
extern PROBLEM_DESCRIPTION_STRUCT              **pd_glob, *pd;

extern PROBLEM_GRAPH_STRUCT                     *pg;

extern struct Action_Flags			*af;
extern BASIS_FUNCTIONS_STRUCT			**bf;
extern BASIS_FUNCTIONS_STRUCT			**bfd;
extern BASIS_FUNCTIONS_STRUCT			**bfi;
extern BASIS_FUNCTIONS_STRUCT			**bfex;
extern struct Shell_Block                       **shell_blocks;
extern struct Field_Variables			*fv, *fv_sens;
extern struct Diet_Field_Variables		*fv_dot_dot, *fv_dot_dot_old;
extern struct Diet_Field_Variables		*fv_old, *fv_dot, *fv_dot_old;
extern struct External_Field_Variables          *efv;
extern struct Porous_Media_Variables            *pmv, *pmv_old;
extern PMV_ML_STRUCT                            *pmv_ml;
extern struct Constitutive_Relations		**cr_glob, *cr;
extern struct Local_Element_Contributions	*lec;
extern struct Transient_Information             *tran;
extern struct Library_IO                        *libio;
extern struct Eigensolver_Info                  *eigen;
extern struct Continuation_Information          *cont;
extern struct Loca_Input                        *loca_in;
extern struct AC_Information                    *augc;
extern struct Continuation_Conditions           *cpcc, *tpcc;
extern struct User_Continuation_Info            *cpuc, *tpuc;
extern struct HC_Information                    *hunt;
extern struct Rotation_Vectors                  ****rotation;
extern struct Level_Set_Data                    *ls ;
extern struct Level_Set_Interface               *lsi;
extern STABILIZATION_PARAMS_STRUCT              *Stab;
extern struct Phase_Function_Data               *pfd;

extern struct Lubrication_Auxiliaries           *LubAux;
extern struct Lubrication_Auxiliaries           *LubAux_old;

#endif /* _MM_AS_H */
