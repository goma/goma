#include "el_elm.h"
#include "rf_fem_const.h"

#ifndef SHELL_TFMP_STRUCT_H
#define SHELL_TFMP_STRUCT_H

struct Gap_Descriptor {
  // h is the thickness of the gap
  double h;
  double dh_dmesh[DIM][MDE];
  double dh_dnormal[DIM][MDE];
  double dh_dtime;
  double d2h_dtime_dmesh[DIM][MDE];
  double d2h_dtime_dnormal[DIM][MDE];
  double gradII_h[DIM];
  double d_gradIIh_dmesh[DIM][DIM][MDE];
  double d_gradIIh_dnormal[DIM][DIM][MDE];
  // extra info needed for computation of the above values
  // time bits
  double time;
  double tt;
  double delta_t;
  // positional bits
  // needed when integration is wrt s where ds = sqrt(dx^2 +dy^2)
  //         AND h0 model is ROLLER
  double gradII_x[DIM];
  double d_gradII_x_dmesh[DIM][MDE];
  int *n_dof;
  int *dof_map;
};
typedef struct Gap_Descriptor GAP_STRUCT;

#endif // SHELL_TFMP_STRUCT_H
