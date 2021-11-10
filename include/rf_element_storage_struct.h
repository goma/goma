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

#ifndef GOMA_RF_ELEMENT_STORAGE_STRUCT_H
#define GOMA_RF_ELEMENT_STORAGE_STRUCT_H

struct Element_Storage {

  /*
   * int Glob_elem_number; Such a beast is available in the dpi struct.
   *                       It would be logical to include it here.
   *                       Note, the processor element number is
   *                       implicitly available via the index that
   *                       is used to access this structure, so we don't
   *                       need to duplicate that here.
   *                       However, we don't at the moment need this
   *                       field, so I will leave it commented out.
   */
  double *Sat_QP_tn; /*
                      * Storage of the saturation at the quadrature
                      * points. This is a vector of doubles.
                      * (current time)
                      */
  double *p_cap_QP;  /*
                      * Storage of the capillary pressure
                      * at the quadrature
                      * points. This is a vector of doubles.
                      * (previous time step)
                      */

  /*PRS: These should be ints, but would necessitate rewriting all of
   *HKM's allocation in rf_element_storage.c.  Pardon the slop.
   */
  double *sat_curve_type;     /* Curve type at current gauss point and time step.
                               * =1. (true) means we are on wetting curve
                               * =0. (false) means we are on the drying curve.
                               */
  double *sat_curve_type_old; /* Switching status at old time at gauss point.
                               * =1. (true) means we are on wetting curve
                               * =0. (false) means we are on the drying curve.
                               */

  double *solidified; /* 0 means material at gaus point has not solidified
                       * 1 means solidified
                       */
};
typedef struct Element_Storage ELEMENT_STORAGE_STRUCT;

#endif
