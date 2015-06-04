#ifndef _MM_SEGREGATED_STRUCTS_H
#define _MM_SEGREGATED_STRUCTS_H

#include "el_elm.h"
#include "mm_mp_const.h"
#include "rf_bc_const.h"
#include "rf_vars_const.h"
#include "mm_elem_block_structs.h"
#include "mm_segregated_structs.h"

struct SplitB_Coupled_Field_Variables {
  dbl grad_v_old[MDE][DIM];
  dbl grad_v_star[MDE][DIM];
  dbl grad_P_star[MDE];
  dbl div_v_star;
  dbl P_star;
  dbl P_old;
  dbl v_star[DIM];
  dbl v_old[DIM];
};

struct SplitB_Element_Stiffness_Pointers {
  dbl v_star[DIM][MDE];                    /* v_star[DIM][MDE], velocity* segregated */
  dbl v_old[DIM][MDE];                      /* v[DIM][MDE], velocity */
  dbl P_old[MDE];                       /* P[MDE], pressure */
  dbl P_star[MDE];                       /* P_star[MDE], pressure */
};

#endif /* _MM_SEGREGATED_STRUCTS_H */
