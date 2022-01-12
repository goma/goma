#ifndef GOMA_DP_GHOST_H
#define GOMA_DP_GHOST_H

#ifdef __cplusplus
extern "C" {
#endif

#include "exo_struct.h"
#include "dpi.h"
#include "mm_eh.h"

goma_error generate_ghost_elems(Exo_DB *exo, Dpi *dpi);

#ifdef __cplusplus
};
#endif

#endif // GOMA_DP_GHOST_H
