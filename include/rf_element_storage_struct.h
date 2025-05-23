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
                               * =1. (true) means we are on draining/drying curve
                               * =0. (false) means we are on the wetting/imbibition curve.
                               */
  double *sat_curve_type_old; /* Switching status at old time at gauss point.
                               * =1. (true) means we are on draining/drying curve
                               * =0. (false) means we are on the wetting/imbibition curve.
                               */

  double *solidified; /* 0 means material at gaus point has not solidified
                       * 1 means solidified
                       */

  double *sat_node;   /*
                       * Storage of the saturation at the node
                       * points. This is a vector of doubles.
                       * (current time)
                       */
  double *p_cap_node; /*
                       * Storage of the capillary pressure
                       * at the node
                       * points. This is a vector of doubles.
                       * (previous time step)
                       */

  int *num_switch; /* number of curve switching at current Gauss point and time step */
  int *switch_now; /* Toggle on whether we are switching curve at current Gauss step and time step
                    * 0 means we are not switching at this time step
                    * 1 means we are switching
                    */

  double *sat_min_wet; /*
                        * Storage of the minimum saturation at the quadrature
                        * points when on wetting curve . This is a vector of doubles.
                        * (current time)
                        */

  double *sat_max_dry; /*
                        * Storage of the maximum saturation at the quadrature
                        * points when on drying curve . This is a vector of doubles.
                        * (current time)
                        */
};
typedef struct Element_Storage ELEMENT_STORAGE_STRUCT;

#endif
