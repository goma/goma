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

#ifndef GOMA_RF_ELEMENT_STORAGE_CONST_H
#define GOMA_RF_ELEMENT_STORAGE_CONST_H

/*
 * prototypes for functions in rf_element_storage.c
 */
extern void setup_element_storage(void);
extern void init_element_storage(ELEM_BLK_STRUCT *, int);
extern void set_init_Element_Storage(ELEM_BLK_STRUCT *, int);
extern void free_element_blocks(Exo_DB *exo);
extern void free_element_storage(Exo_DB *exo);
extern void free_elemStorage(ELEM_BLK_STRUCT *);
extern double get_nodalSat_tnm1_FromES(int);
extern double get_Sat_tnm1_FromES(int);
extern void put_nodalSat_tn_IntoES(int, double);
extern void put_Sat_tn_IntoES(int, double);

#endif
