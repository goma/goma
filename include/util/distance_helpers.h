#ifndef GOMA_DISTANCE_HELPERS_H
#define GOMA_DISTANCE_HELPERS_H

#ifdef __cplusplus
extern "C" {
#endif
#include "dpi.h"
#include "exo_struct.h"
#include "mm_eh.h"

// Find the current distance for all nodes given the nodesets and sidesets to compare against
// returns distance array of size exo->num_nodes
goma_error find_current_distances(Exo_DB *exo,
                                  Dpi *dpi,
                                  double *solution_vector,
                                  bool apply_displacements,
                                  int num_ns,
                                  int *ns_ids,
                                  int num_ss,
                                  int *ss_ids,
                                  double *distances);

#ifdef __cplusplus
}
#endif

#endif // GOMA_DISTANCE_HELPERS_H