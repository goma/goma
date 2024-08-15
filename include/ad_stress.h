#include "ad_turbulence.h"

#ifdef __cplusplus
extern "C" {
#endif

#include "std.h"

int ad_assemble_momentum(dbl time,       /* current time */
                         dbl tt,         /* parameter to vary time integration from
                                            explicit (tt = 1) to implicit (tt = 0) */
                         dbl dt,         /* current time step size */
                         dbl h_elem_avg, /* average global element size for PSPG*/
                         const PG_DATA *pg_data,
                         double xi[DIM], /* Local stu coordinates */
                         const Exo_DB *exo);

int ad_assemble_continuity(dbl time_value, /* current time */
                           dbl tt,         /* parameter to vary time integration from
                                              explicit (tt = 1) to implicit (tt = 0)    */
                           dbl dt,         /* current time step size                    */
                           const PG_DATA *pg_data);
int ad_assemble_stress_sqrt_conf(dbl tt, /* parameter to vary time integration from
                                          * explicit (tt = 1) to implicit (tt = 0) */
                                 dbl dt, /* current time step size */
                                 PG_DATA *pg_data);

dbl ad_viscosity_wrap(struct Generalized_Newtonian *gn_local);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
void ad_ve_polymer_stress(ADType gamma[DIM][DIM], ADType stress[DIM][DIM]);
ADType ad_viscosity(struct Generalized_Newtonian *gn_local, ADType gamma_dot[DIM][DIM]);
ADType ad_vec_dot(const int n1, ADType *v1, ADType *v2);
int ad_tensor_dot(ADType t1[DIM][DIM],
                  ADType t2[DIM][DIM],
                  ADType t1_dot_t2[DIM][DIM],
                  const int dim);
#endif