
#include <exodusII.h>
#include <ext/alloc_traits.h>
#include <mpi.h>
#include <stdlib.h>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cstdio>
#include <cstring>
#include <array>
#include <iterator>
#include <memory>
#include <utility>

#include "dp_ghost.h"
#include "el_elm.h"
#include "rf_mp.h"
#include "std.h"

extern "C" {
#define DISABLE_CPP
#include "el_elm_info.h"
#include "mm_elem_block_structs.h"
#include "rd_mesh.h"
#include "mm_mp.h"

struct Material_Properties;

#undef DISABLE_CPP
}

struct shared_elem {
  int block;
  int type;
  int id;
  int n_nodes;
};

struct shared_node {
  double coord[3];
  int id;
};

inline bool operator<(const shared_node &lhs, const shared_node &rhs) { return lhs.id < rhs.id; }
inline bool operator==(const shared_node &lhs, const shared_node &rhs) { return lhs.id == rhs.id; }
inline bool operator<(const shared_elem &lhs, const shared_elem &rhs) { return lhs.id < rhs.id; }
inline bool operator==(const shared_elem &lhs, const shared_elem &rhs) { return lhs.id == rhs.id; }

template <> struct std::hash<shared_node> {
  std::size_t operator()(shared_node const &node) const noexcept { return node.id; }
};

template <> struct std::hash<shared_elem> {
  std::size_t operator()(shared_elem const &elem) const noexcept { return elem.id; }
};

goma_error generate_ghost_elems(Exo_DB *exo, Dpi *dpi) {
  MPI_Datatype shared_elem_mpi;
  MPI_Type_contiguous(sizeof(shared_elem), MPI_BYTE, &shared_elem_mpi);
  MPI_Type_commit(&shared_elem_mpi);

  MPI_Datatype shared_node_mpi;
  MPI_Type_contiguous(sizeof(shared_node), MPI_BYTE, &shared_node_mpi);
  MPI_Type_commit(&shared_node_mpi);

  // setup num_nodes_per_elem
  std::vector<int> eb_num_nodes_per_elem(exo->num_elem_blocks);
  std::fill(eb_num_nodes_per_elem.begin(), eb_num_nodes_per_elem.end(), 0);
  MPI_Allreduce(exo->eb_num_nodes_per_elem, eb_num_nodes_per_elem.data(), exo->num_elem_blocks,
                MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  for (unsigned int i = 0; i < exo->num_elem_blocks; i++) {
    exo->eb_num_nodes_per_elem[i] = eb_num_nodes_per_elem[i];
  }

  // setup eb_num_attr
  std::vector<int> eb_num_attr(exo->num_elem_blocks);
  MPI_Allreduce(exo->eb_num_attr, eb_num_attr.data(), exo->num_elem_blocks, MPI_INT, MPI_MAX,
                MPI_COMM_WORLD);
  for (unsigned int i = 0; i < exo->num_elem_blocks; i++) {
    exo->eb_num_attr[i] = eb_num_attr[i];
  }

  // setup eb_elem_type
  std::vector<std::array<char, MAX_STR_LENGTH>> local_elem_type(exo->num_elem_blocks);
  std::vector<std::array<char, MAX_STR_LENGTH>> global_elem_type(Num_Proc * exo->num_elem_blocks);
  for (int i = 0; i < exo->num_elem_blocks; i++) {
    memcpy(local_elem_type[i].data(), exo->eb_elem_type[i], MAX_STR_LENGTH);
  }

  MPI_Allgather(local_elem_type.data(), MAX_STR_LENGTH * exo->num_elem_blocks, MPI_BYTE,
                global_elem_type.data(), MAX_STR_LENGTH * exo->num_elem_blocks, MPI_BYTE,
                MPI_COMM_WORLD);

  for (int i = 0; i < exo->num_elem_blocks; i++) {
    if (strncmp((const char *)local_elem_type[i].data(), "NULL", MAX_STR_LENGTH) == 0) {
      for (int j = 0; j < Num_Proc; j++) {
        if (strncmp((const char *)global_elem_type[j * exo->num_elem_blocks + i].data(), "NULL",
                    MAX_STR_LENGTH) != 0) {
          memcpy(exo->eb_elem_type[i], global_elem_type[j * exo->num_elem_blocks + i].data(),
                 MAX_STR_LENGTH);
        }
      }
    }
  }

  // collect nodes/elems this processor needs to send to lower processors
  std::vector<std::vector<shared_node>> external_nodes;

  std::vector<std::vector<shared_elem>> shared_elems_neighbor(dpi->num_neighbors);
  std::vector<std::vector<shared_node>> shared_nodes_neighbor(dpi->num_neighbors);
  std::vector<std::unordered_set<int>> saved_nodes(dpi->num_neighbors);
  std::vector<std::vector<int>> connectivity(dpi->num_neighbors);

  int offset = 0;
  int elem_offset = 0;
  for (int i = 0; i < exo->num_elem_blocks; i++) {
    if (i > 0) {
      elem_offset += exo->eb_num_elems[i - 1];
    }
    for (int j = 0; j < exo->eb_num_elems[i]; j++) {
      int owner = dpi->elem_owner[elem_offset + j];
      if (owner != ProcID) {
        int n_index = -1;
        for (int k = 0; k < dpi->num_neighbors; k++) {
          if (dpi->neighbor[k] == owner) {
            n_index = k;
            break;
          }
        }
        int nnode_per_elem = exo->eb_num_nodes_per_elem[i];
        GOMA_ASSERT(n_index != -1);
        int type =
            get_type(exo->eb_elem_type[i], exo->eb_num_nodes_per_elem[i], exo->eb_num_attr[i]);
        shared_elem elem{dpi->global_elem_block_ids[i], type, dpi->elem_index_global[offset + j],
                         exo->eb_num_nodes_per_elem[i]};
        shared_elems_neighbor[n_index].emplace_back(elem);
        for (int k = 0; k < nnode_per_elem; k++) {
          int local_node = exo->eb_conn[i][j * nnode_per_elem + k];
          connectivity[n_index].push_back(dpi->node_index_global[local_node]);
          if (saved_nodes[n_index].find(local_node) == saved_nodes[n_index].end()) {
            saved_nodes[n_index].insert(local_node);
            double coord[3] = {0.0, 0.0, 0.0};
            coord[0] = exo->x_coord[local_node];
            if (exo->num_dim > 1) {
              coord[1] = exo->y_coord[local_node];
            }
            if (exo->num_dim > 2) {
              coord[2] = exo->z_coord[local_node];
            }

            struct shared_node node;
            node.coord[0] = coord[0];
            node.coord[1] = coord[1];
            node.coord[2] = coord[2];
            node.id = dpi->node_index_global[local_node];
            shared_nodes_neighbor[n_index].emplace_back(node);
          }
        }
      }
    }
  }

  std::vector<std::array<int, 3>> recv_sizes(dpi->num_neighbors);
  std::vector<std::array<int, 3>> send_sizes;
  std::vector<MPI_Request> requests(dpi->num_neighbors * 2);

  for (int i = 0; i < dpi->num_neighbors; i++) {
    std::array<int, 3> sizes{static_cast<int>(shared_elems_neighbor[i].size()),
                             static_cast<int>(connectivity[i].size()),
                             static_cast<int>(shared_nodes_neighbor[i].size())};
    send_sizes.emplace_back(sizes);

    MPI_Irecv(recv_sizes[i].data(), 3, MPI_INT, dpi->neighbor[i], 2500, MPI_COMM_WORLD,
              &requests[i]);
    MPI_Isend(send_sizes[i].data(), 3, MPI_INT, dpi->neighbor[i], 2500, MPI_COMM_WORLD,
              &requests[dpi->num_neighbors + i]);
  }

  MPI_Waitall(dpi->num_neighbors, requests.data(), MPI_STATUSES_IGNORE);

  std::vector<std::vector<shared_elem>> recv_elems_neighbor(dpi->num_neighbors);
  std::vector<std::vector<shared_node>> recv_shared_node(dpi->num_neighbors);
  std::vector<std::vector<int>> recv_connectivity(dpi->num_neighbors);
  std::vector<MPI_Request> many_requests(3 * dpi->num_neighbors * 2);
  int req_index = 0;
  for (int i = 0; i < dpi->num_neighbors; i++) {
    if (recv_sizes[i][0] > 0) {
      recv_elems_neighbor[i].resize(recv_sizes[i][0]);
      MPI_Irecv(recv_elems_neighbor[i].data(), recv_sizes[i][0] * sizeof(shared_elem), MPI_BYTE,
                dpi->neighbor[i], 2500, MPI_COMM_WORLD, &many_requests[req_index++]);
    }
    if (send_sizes[i][0] > 0) {
      GOMA_ASSERT(send_sizes[i][0] == static_cast<int>(shared_elems_neighbor[i].size()));
      MPI_Isend(shared_elems_neighbor[i].data(), send_sizes[i][0] * sizeof(shared_elem), MPI_BYTE,
                dpi->neighbor[i], 2500, MPI_COMM_WORLD, &many_requests[req_index++]);
    }

    if (recv_sizes[i][1] > 0) {
      recv_connectivity[i].resize(recv_sizes[i][1]);
      MPI_Irecv(recv_connectivity[i].data(), recv_sizes[i][1], MPI_INT, dpi->neighbor[i], 2501,
                MPI_COMM_WORLD, &many_requests[req_index++]);
    }
    if (send_sizes[i][1] > 0) {
      GOMA_ASSERT(send_sizes[i][1] == static_cast<int>(connectivity[i].size()));
      MPI_Isend(connectivity[i].data(), send_sizes[i][1], MPI_INT, dpi->neighbor[i], 2501,
                MPI_COMM_WORLD, &many_requests[req_index++]);
    }
    if (recv_sizes[i][2] > 0) {
      recv_shared_node[i].resize(recv_sizes[i][2]);
      MPI_Irecv(recv_shared_node[i].data(), recv_sizes[i][2] * sizeof(shared_node), MPI_BYTE,
                dpi->neighbor[i], 2502, MPI_COMM_WORLD, &many_requests[req_index++]);
    }
    if (send_sizes[i][2] > 0) {
      GOMA_ASSERT(send_sizes[i][2] == static_cast<int>(shared_nodes_neighbor[i].size()));
      MPI_Isend(shared_nodes_neighbor[i].data(), send_sizes[i][2] * sizeof(shared_node), MPI_BYTE,
                dpi->neighbor[i], 2502, MPI_COMM_WORLD, &many_requests[req_index++]);
    }
  }

  MPI_Waitall(req_index, many_requests.data(), MPI_STATUSES_IGNORE);

  // These elems should all be unique, there should probably be some shared nodes however
  std::vector<int> global_nodes;
  std::unordered_set<int> global_nodes_set;
  for (int i = 0; i < exo->num_nodes; i++) {
    global_nodes.push_back(dpi->node_index_global[i]);
    global_nodes_set.insert(dpi->node_index_global[i]);
  }

  std::unordered_set<shared_node> new_nodes;
  std::unordered_set<int> my_boundary_nodes;
  std::unordered_set<int> new_nodes_ids;
  std::unordered_map<int, int> node_owner_new;
  for (int i = 0; i < dpi->num_neighbors; i++) {
    for (auto &ref : recv_shared_node[i]) {
      if (global_nodes_set.find(ref.id) == global_nodes_set.end()) {
        if (new_nodes_ids.find(ref.id) == new_nodes_ids.end()) {
          node_owner_new.insert(std::make_pair(ref.id, dpi->neighbor[i]));
          global_nodes.push_back(ref.id);
          new_nodes.insert(ref);
          new_nodes_ids.insert(ref.id);
        }
      } else {
        if (dpi->neighbor[i] > ProcID) {
          my_boundary_nodes.insert(ref.id);
        }
      }
    }
  }
  GOMA_ASSERT(new_nodes.size() == new_nodes_ids.size());

  // mapping of local to global
  std::unordered_map<int, int> global_to_local;
  std::unordered_map<int, int> local_to_global;
  for (int i = 0; i < static_cast<int>(global_nodes.size()); i++) {
    global_to_local.insert(std::make_pair(global_nodes[i], i));
    local_to_global.insert(std::make_pair(i, global_nodes[i]));
  }

  int num_nodes_new = exo->num_nodes;
  int num_old_nodes = exo->num_nodes;
  for (auto &node : new_nodes) {
    num_nodes_new++;
  }
  GOMA_ASSERT(num_nodes_new == static_cast<int>(global_nodes.size()));

  std::unordered_map<int, int> elem_local_to_global;

  // assume we only have one block to start
  int local_elem_offset = 0;
  for (int block = 0; block < dpi->num_elem_blocks_global; block++) {
    std::vector<shared_elem> elems_new;
    std::vector<int> block_connectivity;
    for (int i = 0; i < dpi->num_neighbors; i++) {
      int conn_offset = 0;
      for (auto &ref : recv_elems_neighbor[i]) {
        if (ref.block == dpi->global_elem_block_ids[block]) {
          elems_new.push_back(ref);
          for (int j = 0; j < ref.n_nodes; j++) {
            block_connectivity.push_back(global_to_local.at(recv_connectivity[i][conn_offset + j]));
          }
        }
        conn_offset += ref.n_nodes;
      }
    }

    for (int i = 0; i < exo->eb_num_elems[block]; i++) {
      elem_local_to_global.insert(std::make_pair(local_elem_offset, dpi->elem_index_global[i]));
      local_elem_offset++;
    }
    if (elems_new.size() == 0) {
      continue;
    }

    int n_nodes = exo->eb_num_nodes_per_elem[block];
    int *int_ptr =
        (int *)realloc(exo->eb_conn[block], sizeof(int) * (exo->eb_num_elems[block] * n_nodes +
                                                           elems_new.size() * n_nodes));
    GOMA_ASSERT(int_ptr != NULL);
    exo->eb_conn[block] = int_ptr;

    // setup new connectivity
    int block_offset = 0;
    for (int i = 0; i < elems_new.size(); i++) {
      auto &elem = elems_new[i];
      elem_local_to_global.insert(std::make_pair(local_elem_offset, elem.id));
      GOMA_ASSERT(elem.n_nodes == exo->eb_num_nodes_per_elem[block]);
      for (int j = 0; j < n_nodes; j++) {
        exo->eb_conn[block][block_offset * n_nodes + j] = block_connectivity[i * n_nodes + j];
      }
      block_offset++;
      local_elem_offset++;
    }

    exo->eb_num_elems[block] += elems_new.size();
  }

  int total_elems = 0;
  for (int block = 0; block < dpi->num_elem_blocks_global; block++) {
    total_elems += exo->eb_num_elems[block];
  }

  exo->num_elems = total_elems;
  int *int_ptr = (int *)realloc(dpi->elem_index_global, sizeof(int) * total_elems);
  GOMA_ASSERT(int_ptr != NULL);
  dpi->elem_index_global = int_ptr;
  for (int i = 0; i < exo->num_elems; i++) {
    dpi->elem_index_global[i] = elem_local_to_global.at(i);
  }
  int_ptr = (int *)realloc(exo->elem_map, sizeof(int) * total_elems);
  GOMA_ASSERT(int_ptr != NULL);
  exo->elem_map = int_ptr;
  memcpy(exo->elem_map, dpi->elem_index_global, sizeof(int) * total_elems);

  exo->num_nodes = num_nodes_new;
  int_ptr = (int *)realloc(dpi->node_index_global, sizeof(int) * num_nodes_new);
  GOMA_ASSERT(int_ptr != NULL);
  dpi->node_index_global = int_ptr;
  for (int i = 0; i < exo->num_nodes; i++) {
    dpi->node_index_global[i] = local_to_global.at(i);
  }
  int_ptr = (int *)realloc(exo->node_map, sizeof(int) * exo->num_nodes);
  GOMA_ASSERT(int_ptr != NULL);
  exo->node_map = int_ptr;
  memcpy(exo->node_map, dpi->node_index_global, sizeof(int) * exo->num_nodes);

  double *double_ptr = (double *)realloc(exo->x_coord, sizeof(double) * exo->num_nodes);
  GOMA_ASSERT(double_ptr != NULL);
  exo->x_coord = double_ptr;

  if (exo->num_dim > 1) {
    double_ptr = (double *)realloc(exo->y_coord, sizeof(double) * exo->num_nodes);
    GOMA_ASSERT(double_ptr != NULL);
    exo->y_coord = double_ptr;
  }

  if (exo->num_dim > 2) {
    double_ptr = (double *)realloc(exo->z_coord, sizeof(double) * exo->num_nodes);
    GOMA_ASSERT(double_ptr != NULL);
    exo->z_coord = double_ptr;
  }

  for (auto &node : new_nodes) {
    int local_id = global_to_local.at(node.id);
    exo->x_coord[local_id] = node.coord[0];
    if (exo->num_dim > 1) {
      exo->y_coord[local_id] = node.coord[1];
    }
    if (exo->num_dim > 2) {
      exo->z_coord[local_id] = node.coord[2];
    }
  }

  int_ptr = (int *)realloc(exo->elem_eb, sizeof(int) * exo->num_elems);
  GOMA_ASSERT(int_ptr != NULL);
  exo->elem_eb = int_ptr;

  offset = 0;
  for (int i = 0; i < exo->num_elem_blocks; i++) {
    for (int j = 0; j < exo->eb_num_elems[i]; j++) {
      exo->elem_eb[offset++] = i;
    }
    exo->eb_ptr[i + 1] = exo->eb_ptr[i] + exo->eb_num_elems[i];
  }

  int_ptr = (int *)realloc(dpi->node_owner, sizeof(int) * exo->num_nodes);
  GOMA_ASSERT(int_ptr != NULL);
  dpi->node_owner = int_ptr;

  for (int i = num_old_nodes; i < exo->num_nodes; i++) {
    int global_id = local_to_global[i];
    dpi->node_owner[i] = node_owner_new.at(global_id);
  }

  std::unordered_set<int> new_external_nodes;
  // all new nodes
  for (int i = num_old_nodes; i < exo->num_nodes; i++) {
    new_external_nodes.insert(i);
  }
  // all nodes with owners < ProcID
  int changed_to_external = 0;
  for (int i = 0; i < num_old_nodes; i++) {
    if (dpi->node_owner[i] < ProcID) {
      new_external_nodes.insert(i);
      changed_to_external++;
    }
  }

  dpi->num_external_nodes = exo->num_nodes - num_old_nodes + changed_to_external;

  std::unordered_set<int> new_boundary_nodes;
  for (int i = 0; i < dpi->num_neighbors; i++) {
    for (int j = 0; j < static_cast<int>(shared_nodes_neighbor[i].size()); j++) {
      int local_id = global_to_local.at(shared_nodes_neighbor[i][j].id);
      if (dpi->node_owner[local_id] >= ProcID) {
        new_boundary_nodes.insert(shared_nodes_neighbor[i][j].id);
      }
    }
  }

  for (auto &bn : my_boundary_nodes) {
    int local_id = global_to_local.at(bn);
    if (dpi->node_owner[local_id] >= ProcID) {
      new_boundary_nodes.insert(bn);
    }
  }
  for (auto &id : new_external_nodes) {
    if (new_boundary_nodes.find(dpi->node_index_global[id]) != new_boundary_nodes.end()) {
      printf("found external in boundary %d\n", id);
    }
  }

  // Reorder nodes to internal, boundary, external
  std::unordered_map<int, int> old_to_new_map;
  std::unordered_set<int> added_bn;
  int boundary_offset = num_old_nodes - changed_to_external - new_boundary_nodes.size();
  int external_offset = num_old_nodes - changed_to_external;
  int internal_offset = 0;
  for (int i = 0; i < num_old_nodes; i++) {
    if (new_external_nodes.find(i) != new_external_nodes.end()) {
      old_to_new_map.insert(std::make_pair(i, external_offset));
      external_offset++;
    } else if (new_boundary_nodes.find(dpi->node_index_global[i]) != new_boundary_nodes.end()) {
      old_to_new_map.insert(std::make_pair(i, boundary_offset));
      added_bn.insert(dpi->node_index_global[i]);
      boundary_offset++;
    } else {
      old_to_new_map.insert(std::make_pair(i, internal_offset));
      internal_offset++;
    }
  }
  std::unordered_set<int> result;
  std::set_difference(new_boundary_nodes.begin(), new_boundary_nodes.end(), added_bn.begin(),
                      added_bn.end(), std::inserter(result, result.end()));

  for (auto &item : result) {
    printf("missing %d\n", item);
  }

  GOMA_ASSERT(external_offset == num_old_nodes);
  GOMA_ASSERT(boundary_offset == num_old_nodes - changed_to_external);
  GOMA_ASSERT(internal_offset ==
              (num_old_nodes - changed_to_external - static_cast<int>(new_boundary_nodes.size())));

  for (int i = num_old_nodes; i < exo->num_nodes; i++) {
    old_to_new_map.insert(std::make_pair(i, i));
  }

  int *tmp_node_index_global = (int *)malloc(sizeof(int) * exo->num_nodes);
  memcpy(tmp_node_index_global, dpi->node_index_global, sizeof(int) * exo->num_nodes);
  int *tmp_node_owner = (int *)malloc(sizeof(int) * exo->num_nodes);
  memcpy(tmp_node_owner, dpi->node_owner, sizeof(int) * exo->num_nodes);
  double *tmp_x = (double *)malloc(sizeof(double) * exo->num_nodes);
  memcpy(tmp_x, exo->x_coord, sizeof(double) * exo->num_nodes);
  double *tmp_y = NULL;
  if (exo->num_dim > 1) {
    tmp_y = (double *)malloc(sizeof(double) * exo->num_nodes);
    memcpy(tmp_y, exo->y_coord, sizeof(double) * exo->num_nodes);
  }
  double *tmp_z = NULL;
  if (exo->num_dim > 2) {
    tmp_z = (double *)malloc(sizeof(double) * exo->num_nodes);
    memcpy(tmp_z, exo->z_coord, sizeof(double) * exo->num_nodes);
  }

  for (int i = 0; i < exo->num_elem_blocks; i++) {
    int *tmp_eb_conn =
        (int *)malloc(sizeof(int) * exo->eb_num_elems[i] * exo->eb_num_nodes_per_elem[i]);
    memcpy(tmp_eb_conn, exo->eb_conn[i],
           sizeof(int) * exo->eb_num_elems[i] * exo->eb_num_nodes_per_elem[i]);
    for (int j = 0; j < exo->eb_num_elems[i] * exo->eb_num_nodes_per_elem[i]; j++) {
      exo->eb_conn[i][j] = old_to_new_map.at(tmp_eb_conn[j]);
    }
  }

  for (int i = 0; i < exo->num_nodes; i++) {
    int new_id = old_to_new_map.at(i);
    dpi->node_index_global[new_id] = tmp_node_index_global[i];
    dpi->node_owner[new_id] = tmp_node_owner[i];

    exo->x_coord[new_id] = tmp_x[i];
    if (exo->num_dim > 1) {
      exo->y_coord[new_id] = tmp_y[i];
    }
    if (exo->num_dim > 2) {
      exo->z_coord[new_id] = tmp_z[i];
    }
  }

  dpi->num_internal_nodes = internal_offset;
  dpi->num_boundary_nodes = exo->num_nodes - internal_offset - new_external_nodes.size();
  dpi->num_external_nodes = new_external_nodes.size();

  for (int i = 0; i < exo->num_node_sets; i++) {
    for (int j = 0; j < exo->ns_num_nodes[i]; j++) {
      int idx = exo->ns_node_index[i] + j;
      int new_id = old_to_new_map.at(exo->ns_node_list[idx]);
      exo->ns_node_list[idx] = new_id;
    }
  }

  //  for (int i = 0; i < exo->ss_node_len; i++) {
  //    exo->ss_node_dist[i] = old_to_new_map.at(exo->ss_node_list[i]);
  //  }

  int_ptr = (int *)realloc(dpi->elem_owner, exo->num_elems * sizeof(int));
  GOMA_ASSERT(int_ptr != NULL);
  dpi->elem_owner = int_ptr;
  offset = 0;
  for (int i = 0; i < exo->num_elem_blocks; i++) {
    for (int j = 0; j < exo->eb_num_elems[i]; j++) {
      int min_proc = ProcID;
      int nnode_per_elem = exo->eb_num_nodes_per_elem[i];
      for (int k = 0; k < nnode_per_elem; k++) {
        int local_node = exo->eb_conn[i][j * nnode_per_elem + k];
        int proc = dpi->node_owner[local_node];
        if (proc < min_proc) {
          min_proc = proc;
        }
      }
      dpi->elem_owner[offset + j] = min_proc;
    }
    offset += exo->eb_num_elems[i];
  }

  dpi->num_owned_nodes = dpi->num_internal_nodes + dpi->num_boundary_nodes;
  dpi->num_universe_nodes =
      dpi->num_internal_nodes + dpi->num_boundary_nodes + dpi->num_external_nodes;
  free(tmp_node_index_global);
  free(tmp_node_owner);
  free(tmp_x);
  free(tmp_y);
  free(tmp_z);

  /*
         *  Fill in the information in the current
         *  element block structure
   */
  for (int i = 0; i < exo->num_elem_blocks; i++) {
    ELEM_BLK_STRUCT *eb_ptr = NULL;
    eb_ptr = Element_Blocks + i;
    eb_ptr->Elem_Blk_Num = i;
    eb_ptr->Elem_Blk_Id = exo->eb_id[i];
    if (exo->eb_num_elems[i] > 0) {
      eb_ptr->Elem_Type =
          get_type(exo->eb_elem_type[i], exo->eb_num_nodes_per_elem[i], exo->eb_num_attr[i]);
    } else {
      eb_ptr->Elem_Type = NULL_ELEM_TYPE;
    }
    eb_ptr->Num_Nodes_Per_Elem = exo->eb_num_nodes_per_elem[i];
    eb_ptr->Num_Attr_Per_Elem = exo->eb_num_attr[i];
    int mindex = map_mat_index(eb_ptr->Elem_Blk_Id);
    if (mindex < 0) {
      eb_ptr->MatlProp_ptr = NULL;
    } else {
      eb_ptr->MatlProp_ptr = mp_glob[mindex];
    }
    eb_ptr->ElemStorage = NULL;
    eb_ptr->Num_Elems_In_Block = exo->eb_num_elems[i];
    if (exo->eb_num_elems[i] > 0) {
      eb_ptr->IP_total = elem_info(NQUAD, eb_ptr->Elem_Type);
    } else {
      eb_ptr->IP_total = 0;
    }
  }

  return GOMA_SUCCESS;
}
