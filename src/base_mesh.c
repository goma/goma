#include "base_mesh.h"
#include "stdlib.h"
#include "string.h"

#define EXO_TO_BASE(member) base->member = exo->member;

#define MALLOC_COPY(dst, src, size) \
  do {                              \
    size_t sz = size;               \
    dst = malloc(sz);               \
    memcpy(dst, src, sz);           \
  } while (0)

static goma_error setup_base_mesh_serial(Dpi *dpi, Exo_DB *exo);
static goma_error setup_base_mesh_parallel(Dpi *dpi, Exo_DB *exo);

static goma_error free_base_mesh_serial(Exo_DB *exo);
static goma_error free_base_mesh_parallel(Exo_DB *exo);

goma_error setup_base_mesh(Dpi *dpi, Exo_DB *exo, int num_proc) {
  goma_error error;
  if (num_proc == 1) {
    error = setup_base_mesh_serial(dpi, exo);
  } else {
    error = setup_base_mesh_parallel(dpi, exo);
  }
  return error;
}

goma_error free_base_mesh(Exo_DB *exo) {
  goma_error error;
  if (exo->base_mesh_is_serial) {
    error = free_base_mesh_serial(exo);
  } else {
    error = free_base_mesh_parallel(exo);
  }
  return error;
}

static goma_error setup_base_mesh_serial(Dpi *dpi, Exo_DB *exo) {
  // serial mesh is same as Exo_DB so just point to Exo_DB
  struct Exodus_Base *base = malloc(sizeof(struct Exodus_Base));
  exo->base_mesh = base;
  EXO_TO_BASE(num_dim);
  EXO_TO_BASE(num_elems);
  EXO_TO_BASE(num_elem_blocks);
  EXO_TO_BASE(num_node_sets);
  EXO_TO_BASE(num_side_sets);
  EXO_TO_BASE(num_nodes);
  EXO_TO_BASE(x_coord);
  EXO_TO_BASE(y_coord);
  EXO_TO_BASE(z_coord);
  EXO_TO_BASE(eb_id);
  EXO_TO_BASE(eb_elem_type);
  EXO_TO_BASE(eb_num_elems);
  EXO_TO_BASE(eb_num_nodes_per_elem);
  EXO_TO_BASE(eb_conn);
  EXO_TO_BASE(eb_num_attr);

  EXO_TO_BASE(ns_node_len);
  EXO_TO_BASE(ns_id);
  EXO_TO_BASE(ns_num_nodes);
  EXO_TO_BASE(ns_num_distfacts);
  EXO_TO_BASE(ns_node_index);
  EXO_TO_BASE(ns_distfact_index);
  EXO_TO_BASE(ns_node_list);
  EXO_TO_BASE(ns_distfact_list);

  EXO_TO_BASE(ss_elem_len);
  EXO_TO_BASE(ss_node_len);

  EXO_TO_BASE(ss_id);
  EXO_TO_BASE(ss_num_sides);
  EXO_TO_BASE(ss_num_distfacts);
  EXO_TO_BASE(ss_elem_index);
  EXO_TO_BASE(ss_distfact_index);
  EXO_TO_BASE(ss_elem_list);
  EXO_TO_BASE(ss_side_list);
  EXO_TO_BASE(ss_distfact_list);

  EXO_TO_BASE(ns_num_props);
  EXO_TO_BASE(ss_num_props);
  EXO_TO_BASE(eb_num_props);

  EXO_TO_BASE(eb_prop_name);
  EXO_TO_BASE(ns_prop_name);
  EXO_TO_BASE(ss_prop_name);

  EXO_TO_BASE(eb_prop);
  EXO_TO_BASE(ns_prop);
  EXO_TO_BASE(ss_prop);

  EXO_TO_BASE(elem_var_tab);
  exo->base_mesh_is_serial = true;
  return GOMA_SUCCESS;
}

static goma_error setup_base_mesh_parallel(Dpi *dpi, Exo_DB *exo) {
  // we are before rd_dpi so we should not have any ghosted elements
  struct Exodus_Base *base = malloc(sizeof(struct Exodus_Base));
  exo->base_mesh = base;
  EXO_TO_BASE(num_dim);
  EXO_TO_BASE(num_elems);
  EXO_TO_BASE(num_elem_blocks);
  EXO_TO_BASE(num_node_sets);
  EXO_TO_BASE(num_side_sets);
  EXO_TO_BASE(num_nodes);

  MALLOC_COPY(base->x_coord, exo->x_coord, sizeof(dbl) * exo->num_nodes);
  if (base->num_dim > 1) {
    MALLOC_COPY(base->y_coord, exo->y_coord, sizeof(dbl) * exo->num_nodes);
  }
  if (base->num_dim > 2) {
    MALLOC_COPY(base->z_coord, exo->z_coord, sizeof(dbl) * exo->num_nodes);
  }

  base->eb_elem_type = malloc(sizeof(char *) * exo->num_elem_blocks);
  base->eb_conn = malloc(sizeof(char *) * exo->num_elem_blocks);
  for (int i = 0; i < exo->num_elem_blocks; i++) {
    base->eb_elem_type[i] = (char *)malloc(MAX_STR_LENGTH * sizeof(char));
    memcpy(base->eb_elem_type[i], exo->eb_elem_type[i], MAX_STR_LENGTH * sizeof(char));
    if (exo->eb_num_elems[i] > 0) {
      base->eb_conn[i] = malloc(sizeof(int) * exo->eb_num_elems[i] * exo->eb_num_nodes_per_elem[i]);
      memcpy(base->eb_conn[i], exo->eb_conn[i],
             sizeof(int) * exo->eb_num_elems[i] * exo->eb_num_nodes_per_elem[i]);
    } else {
      base->eb_conn[i] = NULL;
    }
  }
  MALLOC_COPY(base->eb_id, exo->eb_id, sizeof(int) * exo->num_elem_blocks);
  MALLOC_COPY(base->eb_num_elems, exo->eb_num_elems, sizeof(int) * exo->num_elem_blocks);
  MALLOC_COPY(base->eb_num_nodes_per_elem, exo->eb_num_nodes_per_elem,
              sizeof(int) * exo->num_elem_blocks);
  MALLOC_COPY(base->eb_num_attr, exo->eb_num_attr, sizeof(int) * exo->num_elem_blocks);

  EXO_TO_BASE(ns_node_len);
  MALLOC_COPY(base->ns_id, exo->ns_id, sizeof(int) * exo->num_node_sets);
  MALLOC_COPY(base->ns_num_nodes, exo->ns_num_nodes, sizeof(int) * exo->num_node_sets);
  MALLOC_COPY(base->ns_num_distfacts, exo->ns_num_distfacts, sizeof(int) * exo->num_node_sets);
  MALLOC_COPY(base->ns_node_index, exo->ns_node_index, sizeof(int) * exo->num_node_sets);
  MALLOC_COPY(base->ns_distfact_index, exo->ns_distfact_index, sizeof(int) * exo->num_node_sets);
  MALLOC_COPY(base->ns_node_list, exo->ns_node_list, sizeof(int) * exo->ns_node_len);
  MALLOC_COPY(base->ns_distfact_list, exo->ns_distfact_list, sizeof(dbl) * exo->ns_distfact_len);

  EXO_TO_BASE(ss_elem_len);
  MALLOC_COPY(base->ss_id, exo->ss_id, sizeof(int) * exo->num_side_sets);
  MALLOC_COPY(base->ss_num_sides, exo->ss_num_sides, sizeof(int) * exo->num_side_sets);
  MALLOC_COPY(base->ss_num_distfacts, exo->ss_num_distfacts, sizeof(int) * exo->num_side_sets);
  MALLOC_COPY(base->ss_elem_index, exo->ss_elem_index, sizeof(int) * exo->num_side_sets);
  MALLOC_COPY(base->ss_distfact_index, exo->ss_distfact_index, sizeof(int) * exo->num_side_sets);
  MALLOC_COPY(base->ss_side_list, exo->ss_side_list, sizeof(int) * exo->ss_elem_len);
  MALLOC_COPY(base->ss_elem_list, exo->ss_elem_list, sizeof(int) * exo->ss_elem_len);
  MALLOC_COPY(base->ss_distfact_list, exo->ss_distfact_list, sizeof(dbl) * exo->ss_distfact_len);

  EXO_TO_BASE(ns_num_props);
  EXO_TO_BASE(ss_num_props);
  EXO_TO_BASE(eb_num_props);

  EXO_TO_BASE(eb_prop_name);
  EXO_TO_BASE(ns_prop_name);
  EXO_TO_BASE(ss_prop_name);

  EXO_TO_BASE(eb_prop);
  EXO_TO_BASE(ns_prop);
  EXO_TO_BASE(ss_prop);

  exo->base_mesh_is_serial = false;
  return GOMA_SUCCESS;
}

static goma_error free_base_mesh_serial(Exo_DB *exo) {
  free(exo->base_mesh->node_map);
  free(exo->base_mesh->elem_map);
  free(exo->base_mesh);
  return GOMA_SUCCESS;
}
static goma_error free_base_mesh_parallel(Exo_DB *exo) {
  struct Exodus_Base *base = exo->base_mesh;
  free(base->x_coord);
  if (base->num_dim > 1) {
    free(base->y_coord);
  }
  if (base->num_dim > 2) {
    free(base->z_coord);
  }

  for (int i = 0; i < exo->num_elem_blocks; i++) {
    free(base->eb_elem_type[i]);
    free(base->eb_conn[i]);
  }
  free(base->eb_elem_type);
  free(base->eb_conn);

  free(base->eb_id);
  free(base->eb_num_elems);
  free(base->eb_num_nodes_per_elem);
  free(base->eb_num_attr);

  free(base->ns_id);
  free(base->ns_num_nodes);
  free(base->ns_num_distfacts);
  free(base->ns_node_index);
  free(base->ns_distfact_index);
  free(base->ns_node_list);
  free(base->ns_distfact_list);

  free(base->ss_id);
  free(base->ss_num_sides);
  free(base->ss_num_distfacts);
  free(base->ss_elem_index);
  free(base->ss_distfact_index);
  free(base->ss_side_list);
  free(base->ss_elem_list);
  free(base->ss_distfact_list);

  free(base->node_map);
  free(base->elem_map);
  free(exo->base_mesh);

  return GOMA_SUCCESS;
}
