#ifndef GOMA_BASE_MESH_H
#define GOMA_BASE_MESH_H

#include "dpi.h"
#include "exo_struct.h"
#include "mm_eh.h"

goma_error setup_base_mesh(Dpi *dpi, Exo_DB *exo, int num_proc);
goma_error free_base_mesh(Exo_DB *exo);

#endif // GOMA_BASE_MESH_H
