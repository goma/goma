#ifndef GOMA_MM_FILL_TURBULENT_H
#define GOMA_MM_FILL_TURBULENT_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MM_FILL_TURBULENT_C
#define EXTERN
#else
#define EXTERN extern
#endif

#include "exo_struct.h"
#include "mm_as_structs.h"
#include "std.h"

EXTERN int assemble_spalart_allmaras(dbl time_value, /* current time */
                                     dbl tt,         /* parameter to vary time integration from
                                                      explicit (tt = 1) to implicit (tt = 0)    */
                                     dbl dt,         /* current time step size                    */
                                     const PG_DATA *pg_data);

#endif
