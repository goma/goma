#ifndef GOMA_MM_FILL_SPLIT_H
#define GOMA_MM_FILL_SPLIT_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MM_FILL_SPLIT_C
#define EXTERN
#else
#define EXTERN extern
#endif

#include "exo_struct.h"
#include "mm_as_structs.h"
#include "std.h"

EXTERN int assemble_ustar(dbl time_value, /* current time */
                          dbl tt,         /* parameter to vary time integration from
                                             explicit (tt = 1) to implicit (tt = 0)    */
                          dbl dt,         /* current time step size                    */
                          const PG_DATA *pg_data);

EXTERN int assemble_pstar(dbl time_value, /* current time */
                          dbl tt,         /* parameter to vary time integration from
                                             explicit (tt = 1) to implicit (tt = 0)    */
                          dbl dt,         /* current time step size                    */
                          const PG_DATA *pg_data);

EXTERN int assemble_continuity_segregated(dbl time_value, /* current time */
                                          dbl tt,         /* parameter to vary time integration from
                                                             explicit (tt = 1) to implicit (tt = 0)    */
                                          dbl dt, /* current time step size                    */
                                          const PG_DATA *pg_data);

EXTERN int assemble_momentum_segregated(dbl time, /* current time */
                                        dbl tt,   /* parameter to vary time integration from
                                                     explicit (tt = 1) to implicit (tt = 0) */
                                        dbl dt,   /* current time step size */
                                        const PG_DATA *pg_data);

#endif
