#ifndef GOMA_AD_POROUS_H
#define GOMA_AD_POROUS_H

#ifdef GOMA_ENABLE_SACADO

#ifdef __cplusplus
extern "C" {
#endif

#include "mm_as_structs.h"
#include "mm_fill_stabilization.h"
#include "mm_mp_structs.h"
#include "std.h"

int ad_assemble_porous_shell_saturation(dbl tt,           // Time integration form
                                     dbl dt,           // Time step size
                                     dbl xi[DIM],      // Current coordinates
                                     const Exo_DB *exo // ExoII handle
);

#ifdef __cplusplus
}
#endif

#endif

#endif // GOMA_AD_POROUS_H