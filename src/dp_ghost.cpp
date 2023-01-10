/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2022 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/
#include <algorithm>
#include <array>
#include <cstdio>
#include <cstring>
#include <exodusII.h>
#include <ext/alloc_traits.h>
#include <iterator>
#include <memory>
#include <mpi.h>
#include <stdlib.h>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// not needed except to avoid including as a C file
#include "sl_epetra_interface.h"

#include "dp_ghost.h"

extern "C" {
#include "el_elm.h"
#include "el_elm_info.h"
#include "exo_conn.h"
#include "mm_elem_block_structs.h"
#include "mm_fill_fill.h"
#include "mm_mp.h"
#include "rd_mesh.h"
#include "rf_mp.h"
#include "std.h"
struct Material_Properties;
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
  int owner;
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
  for (int i = 0; i < exo->num_elem_blocks; i++) {
    exo->eb_num_nodes_per_elem[i] = eb_num_nodes_per_elem[i];
  }

  // setup eb_num_attr
  std::vector<int> eb_num_attr(exo->num_elem_blocks);
  MPI_Allreduce(exo->eb_num_attr, eb_num_attr.data(), exo->num_elem_blocks, MPI_INT, MPI_MAX,
                MPI_COMM_WORLD);
  for (int i = 0; i < exo->num_elem_blocks; i++) {
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
      int nnode_per_elem = exo->eb_num_nodes_per_elem[i];
      // check if neighbors need an elem
      std::unordered_set<int> neighbor_index;
      for (int k = 0; k < nnode_per_elem; k++) {
        int local_node = exo->eb_conn[i][j * nnode_per_elem + k];
        int owner = dpi->node_owner[local_node];
        if (dpi->node_owner[local_node] != ProcID) {
          // neighbor needs element find which negibhor
          int n_index = -1;
          for (int k = 0; k < dpi->num_neighbors; k++) {
            if (dpi->neighbor[k] == owner) {
              neighbor_index.insert(k);
              n_index = k;
            }
          }
          GOMA_ASSERT_ALWAYS(n_index != -1);
        }
      }

      for (auto n_index : neighbor_index) {
        int type =
            get_type(exo->eb_elem_type[i], exo->eb_num_nodes_per_elem[i], exo->eb_num_attr[i]);
        shared_elem elem{dpi->global_elem_block_ids[i], type,
                         dpi->elem_index_global[elem_offset + j], exo->eb_num_nodes_per_elem[i]};
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
            node.owner = dpi->node_owner[local_node];
            shared_nodes_neighbor[n_index].emplace_back(node);
          }
        }
      }
    }
  }

  std::vector<std::array<int, 3>> recv_sizes(dpi->num_neighbors);
  std::vector<std::array<int, 3>> send_sizes(dpi->num_neighbors);
  std::vector<MPI_Request> requests(dpi->num_neighbors * 2);

  for (int i = 0; i < dpi->num_neighbors; i++) {
    std::array<int, 3> sizes{static_cast<int>(shared_elems_neighbor[i].size()),
                             static_cast<int>(connectivity[i].size()),
                             static_cast<int>(shared_nodes_neighbor[i].size())};
    for (unsigned int j = 0; j < sizes.size(); j++) {
      send_sizes[i][j] = sizes[j];
    }

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
          node_owner_new.insert(std::make_pair(ref.id, ref.owner));
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
  num_nodes_new += new_nodes.size();
  GOMA_ASSERT(num_nodes_new == static_cast<int>(global_nodes.size()));

  std::unordered_map<int, int> elem_local_to_global;
  std::unordered_map<int, int> old_elem_to_new;

  // assume we only have one block to start
  int old_elem_offset = 0;
  int local_elem_offset = 0;
  for (int block = 0; block < dpi->num_elem_blocks_global; block++) {
    int block_offset = 0;
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
      elem_local_to_global.insert(
          std::make_pair(local_elem_offset, dpi->elem_index_global[old_elem_offset]));
      old_elem_to_new[old_elem_offset] = local_elem_offset;
      old_elem_offset++;
      local_elem_offset++;
      block_offset++;
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
    for (unsigned int i = 0; i < elems_new.size(); i++) {
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
  std::unordered_map<int, int> elem_global_to_local;
  dpi->elem_index_global = int_ptr;
  for (int i = 0; i < exo->num_elems; i++) {
    dpi->elem_index_global[i] = elem_local_to_global.at(i);
    elem_global_to_local[dpi->elem_index_global[i]] = i;
  }

  exo->num_nodes = num_nodes_new;
  int_ptr = (int *)realloc(dpi->node_index_global, sizeof(int) * num_nodes_new);
  GOMA_ASSERT(int_ptr != NULL);
  dpi->node_index_global = int_ptr;
  for (int i = 0; i < exo->num_nodes; i++) {
    dpi->node_index_global[i] = local_to_global.at(i);
  }

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

  // our neighbor list might have changed so now we have to find out which processors need our nodes
  std::unordered_set<int> neighbors_set;
  for (int i = num_old_nodes; i < exo->num_nodes; i++) {
    int global_id = local_to_global[i];
    dpi->node_owner[i] = node_owner_new.at(global_id);
  }
  for (int i = 0; i < exo->num_nodes; i++) {
    if (dpi->node_owner[i] != ProcID) {
      neighbors_set.insert(dpi->node_owner[i]);
    }
  }

  std::vector<int> neighbor_list(neighbors_set.begin(), neighbors_set.end());

  int my_neighbors = neighbor_list.size();
  int max_neighbors;

  MPI_Allreduce(&my_neighbors, &max_neighbors, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  std::vector<int> all_neighbors(max_neighbors * Num_Proc);
  neighbor_list.resize(max_neighbors);
  for (int i = my_neighbors; i < max_neighbors; i++) {
    neighbor_list[i] = -1;
  }

  MPI_Allgather(neighbor_list.data(), max_neighbors, MPI_INT, all_neighbors.data(), max_neighbors,
                MPI_INT, MPI_COMM_WORLD);

  // loop over the neighbors and decide if we have another neighbor
  for (int i = 0; i < Num_Proc; i++) {
    if (i == ProcID)
      continue;
    for (int j = max_neighbors * i; j < (max_neighbors * i + max_neighbors); j++) {
      if (all_neighbors[j] == ProcID) {
        neighbors_set.insert(i);
      }
    }
  }

  neighbor_list.resize(neighbors_set.size());
  neighbor_list.assign(neighbors_set.begin(), neighbors_set.end());
  std::sort(neighbor_list.begin(), neighbor_list.end());

  dpi->num_neighbors = neighbor_list.size();
  int_ptr = (int *)realloc(dpi->neighbor, sizeof(int) * dpi->num_neighbors);
  GOMA_ASSERT(int_ptr != NULL);
  dpi->neighbor = int_ptr;
  for (int i = 0; i < dpi->num_neighbors; i++) {
    dpi->neighbor[i] = neighbor_list[i];
  }

  requests.resize(2 * dpi->num_neighbors);
  std::vector<MPI_Request> requests_sendrecv(2 * dpi->num_neighbors);
  std::vector<int> num_send_nodes(dpi->num_neighbors);
  std::vector<int> num_recv_nodes(dpi->num_neighbors);

  std::vector<std::vector<int>> global_send_nodes(dpi->num_neighbors);
  std::vector<std::vector<int>> global_recv_nodes(dpi->num_neighbors);
  for (int i = 0; i < dpi->num_neighbors; i++) {
    MPI_Irecv(&num_send_nodes[i], 1, MPI_INT, dpi->neighbor[i], 206, MPI_COMM_WORLD,
              &requests[dpi->num_neighbors + i]);
    // printf("Proc %d recv from Proc %d tag %d", dpi->neighbor[i], ProcID, 206);
  }
  for (int i = 0; i < dpi->num_neighbors; i++) {
    num_recv_nodes[i] = 0;
    for (int j = 0; j < exo->num_nodes; j++) {
      if (dpi->node_owner[j] == dpi->neighbor[i]) {
        num_recv_nodes[i] += 1;
      }
    }
    MPI_Isend(&num_recv_nodes[i], 1, MPI_INT, dpi->neighbor[i], 206, MPI_COMM_WORLD, &requests[i]);
    // printf("Proc %d sending to Proc %d tag %d", ProcID, dpi->neighbor[i], 206);
  }

  MPI_Waitall(dpi->num_neighbors * 2, requests.data(), MPI_STATUSES_IGNORE);

  for (int i = 0; i < dpi->num_neighbors; i++) {
    global_send_nodes[i].resize(num_send_nodes[i]);
    MPI_Irecv(global_send_nodes[i].data(), num_send_nodes[i], MPI_INT, dpi->neighbor[i], 207,
              MPI_COMM_WORLD, &requests[dpi->num_neighbors + i]);
    // printf("Proc %d recv from Proc %d tag %d", dpi->neighbor[i], ProcID, 207);
  }

  for (int i = 0; i < dpi->num_neighbors; i++) {
    global_recv_nodes[i].resize(num_recv_nodes[i]);
    int recv_index = 0;
    for (int j = 0; j < exo->num_nodes; j++) {
      if (dpi->node_owner[j] == dpi->neighbor[i]) {
        global_recv_nodes[i][recv_index++] = dpi->node_index_global[j];
      }
    }
    MPI_Isend(global_recv_nodes[i].data(), num_recv_nodes[i], MPI_INT, dpi->neighbor[i], 207,
              MPI_COMM_WORLD, &requests[i]);
    // printf("Proc %d sending to Proc %d tag %d", ProcID, dpi->neighbor[i], 207);
  }

  MPI_Waitall(dpi->num_neighbors * 2, requests.data(), MPI_STATUSES_IGNORE);

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
    for (int j = 0; j < static_cast<int>(global_send_nodes[i].size()); j++) {
      int local_id = global_to_local.at(global_send_nodes[i][j]);
      if (dpi->node_owner[local_id] >= ProcID) {
        new_boundary_nodes.insert(dpi->node_index_global[local_id]);
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
    free(tmp_eb_conn);
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

  // exchange node sets
  // find boundary nodes that are part of a nodeset
  std::vector<std::vector<int>> ns_boundary_nodes(exo->num_node_sets);
  for (int i = 0; i < exo->num_node_sets; i++) {
    for (int j = 0; j < exo->ns_num_nodes[i]; j++) {
      int offset = j + exo->ns_node_index[i];
      int ln = exo->ns_node_list[offset];
      if (new_boundary_nodes.find(dpi->node_index_global[ln]) != new_boundary_nodes.end()) {
        ns_boundary_nodes[i].push_back(dpi->node_index_global[ln]);
      }
    }
  }

  std::vector<std::vector<int>> ns_recv_nodes(dpi->num_neighbors);

  std::vector<MPI_Request> ns_requests(2 * dpi->num_neighbors * exo->num_node_sets);

  int req_offset = 0;
  // exchange with neighbors
  for (int i = 0; i < dpi->num_neighbors; i++) {
    ns_recv_nodes[i].resize(exo->num_node_sets);
    MPI_Irecv(ns_recv_nodes[i].data(), exo->num_node_sets, MPI_INT, dpi->neighbor[i], 208,
              MPI_COMM_WORLD, &ns_requests[req_offset++]);
  }
  std::vector<int> num_send(exo->num_node_sets);
  for (int j = 0; j < exo->num_node_sets; j++) {
    num_send[j] = ns_boundary_nodes[j].size();
  }
  for (int i = 0; i < dpi->num_neighbors; i++) {
    MPI_Isend(&num_send[0], exo->num_node_sets, MPI_INT, dpi->neighbor[i], 208, MPI_COMM_WORLD,
              &ns_requests[req_offset++]);
  }
  MPI_Waitall(req_offset, ns_requests.data(), MPI_STATUSES_IGNORE);

  std::vector<std::vector<std::vector<int>>> ns_neighbor_external_nodes(dpi->num_neighbors);
  req_offset = 0;
  for (int i = 0; i < dpi->num_neighbors; i++) {
    ns_neighbor_external_nodes[i].resize(exo->num_node_sets);
    for (int ns = 0; ns < exo->num_node_sets; ns++) {
      ns_neighbor_external_nodes[i][ns].resize(ns_recv_nodes[i][ns]);
      MPI_Irecv(ns_neighbor_external_nodes[i][ns].data(), ns_recv_nodes[i][ns], MPI_INT,
                dpi->neighbor[i], 208, MPI_COMM_WORLD, &ns_requests[req_offset++]);
      MPI_Isend(ns_boundary_nodes[ns].data(), num_send[ns], MPI_INT, dpi->neighbor[i], 208,
                MPI_COMM_WORLD, &ns_requests[req_offset++]);
    }
  }
  MPI_Waitall(req_offset, ns_requests.data(), MPI_STATUSES_IGNORE);

  std::unordered_map<int, int> new_global_to_local;
  for (int i = 0; i < exo->num_nodes; i++) {
    new_global_to_local.insert(std::make_pair(dpi->node_index_global[i], i));
  }
  // create nodeset node lists
  std::vector<std::vector<int>> ns_nodes(exo->num_node_sets);
  int total_ns_nodes = 0;
  for (int i = 0; i < exo->num_node_sets; i++) {
    for (int j = 0; j < exo->ns_num_nodes[i]; j++) {
      ns_nodes[i].push_back(exo->ns_node_list[exo->ns_node_index[i] + j]);
    }
    // add external if we have that node
    for (int j = 0; j < dpi->num_neighbors; j++) {
      for (auto &gnn : ns_neighbor_external_nodes[j][i]) {
        if (new_global_to_local.find(gnn) != new_global_to_local.end()) {
          ns_nodes[i].push_back(new_global_to_local[gnn]);
        }
      }
    }
    std::sort(ns_nodes[i].begin(), ns_nodes[i].end());
    total_ns_nodes += ns_nodes[i].size();
  }

  exo->ns_node_len = total_ns_nodes;
  exo->ns_distfact_len = total_ns_nodes;

  int_ptr = (int *)realloc(exo->ns_node_list, sizeof(int) * exo->ns_node_len);
  GOMA_ASSERT(int_ptr != NULL);
  exo->ns_node_list = int_ptr;

  int index = 0;
  for (int i = 0; i < exo->num_node_sets; i++) {
    exo->ns_node_index[i] = index;
    exo->ns_distfact_index[i] = index;
    exo->ns_num_nodes[i] = ns_nodes[i].size();
    exo->ns_num_distfacts[i] = ns_nodes[i].size();

    for (unsigned int j = 0; j < ns_nodes[i].size(); j++) {
      exo->ns_node_list[index] = ns_nodes[i][j];
      index++;
    }
  }

  dbl *dbl_ptr = (dbl *)realloc(exo->ns_distfact_list, sizeof(dbl) * exo->ns_node_len);
  GOMA_ASSERT(dbl_ptr != NULL);
  exo->ns_distfact_list = dbl_ptr;

  for (int i = 0; i < exo->ns_node_len; i++) {
    exo->ns_distfact_list[i] = 0.0;
  }

  // fix ss lists
  for (int ins = 0; ins < exo->num_side_sets; ins++) {
    for (int side_index = 0; side_index < exo->ss_num_sides[ins]; side_index++) {
      for (int lni = exo->ss_node_side_index[ins][side_index];
           lni < exo->ss_node_side_index[ins][side_index + 1]; lni++) {
        int inode = exo->ss_node_list[ins][lni];
        exo->ss_node_list[ins][lni] = old_to_new_map[inode];
      }
    }
  }

  for (int ins = 0; ins < exo->num_side_sets; ins++) {
    for (int j = exo->ss_elem_index[ins]; j < exo->ss_elem_index[ins] + exo->ss_num_sides[ins];
         j++) {
      if (old_elem_to_new.find(exo->ss_elem_list[j]) != old_elem_to_new.end()) {
        exo->ss_elem_list[j] = old_elem_to_new[exo->ss_elem_list[j]];
      }
    }
  }

  // exchange_ss
  std::vector<std::vector<int>> ss_elems(exo->num_side_sets);
  std::vector<std::vector<int>> num_recv_sides(dpi->num_neighbors);
  std::vector<int> num_send_sides(exo->num_side_sets);
  std::vector<std::vector<int>> ss_sides(exo->num_side_sets);
  for (int ins = 0; ins < exo->num_side_sets; ins++) {
    for (int j = exo->ss_elem_index[ins]; j < exo->ss_elem_index[ins] + exo->ss_num_sides[ins];
         j++) {
      ss_elems[ins].push_back(dpi->elem_index_global[exo->ss_elem_list[j]]);
      ss_sides[ins].push_back(exo->ss_side_list[j]);
    }
    num_send_sides[ins] = ss_sides[ins].size();
  }

  std::vector<MPI_Request> ss_requests(2 * dpi->num_neighbors * exo->num_side_sets);

  req_offset = 0;
  // exchange with neighbors
  for (int i = 0; i < dpi->num_neighbors; i++) {
    num_recv_sides[i].resize(exo->num_side_sets);
    MPI_Irecv(num_recv_sides[i].data(), exo->num_side_sets, MPI_INT, dpi->neighbor[i], 210,
              MPI_COMM_WORLD, &ss_requests[req_offset++]);
  }
  for (int i = 0; i < dpi->num_neighbors; i++) {
    MPI_Isend(&num_send_sides[0], exo->num_side_sets, MPI_INT, dpi->neighbor[i], 210,
              MPI_COMM_WORLD, &ss_requests[req_offset++]);
  }
  MPI_Waitall(req_offset, ss_requests.data(), MPI_STATUSES_IGNORE);

  std::vector<std::vector<std::vector<int>>> ss_neighbor_sides(dpi->num_neighbors);
  std::vector<std::vector<std::vector<int>>> ss_neighbor_elem(dpi->num_neighbors);
  req_offset = 0;
  for (int i = 0; i < dpi->num_neighbors; i++) {
    ss_neighbor_sides[i].resize(exo->num_side_sets);
    for (int ss = 0; ss < exo->num_side_sets; ss++) {
      ss_neighbor_sides[i][ss].resize(num_recv_sides[i][ss]);
      MPI_Irecv(ss_neighbor_sides[i][ss].data(), num_recv_sides[i][ss], MPI_INT, dpi->neighbor[i],
                211, MPI_COMM_WORLD, &ss_requests[req_offset++]);
      MPI_Isend(ss_sides[ss].data(), num_send_sides[ss], MPI_INT, dpi->neighbor[i], 211,
                MPI_COMM_WORLD, &ss_requests[req_offset++]);
    }
  }
  MPI_Waitall(req_offset, ss_requests.data(), MPI_STATUSES_IGNORE);

  req_offset = 0;
  for (int i = 0; i < dpi->num_neighbors; i++) {
    ss_neighbor_elem[i].resize(exo->num_side_sets);
    for (int ss = 0; ss < exo->num_side_sets; ss++) {
      ss_neighbor_elem[i][ss].resize(num_recv_sides[i][ss]);
      MPI_Irecv(ss_neighbor_elem[i][ss].data(), num_recv_sides[i][ss], MPI_INT, dpi->neighbor[i],
                211, MPI_COMM_WORLD, &ss_requests[req_offset++]);
      MPI_Isend(ss_elems[ss].data(), num_send_sides[ss], MPI_INT, dpi->neighbor[i], 211,
                MPI_COMM_WORLD, &ss_requests[req_offset++]);
    }
  }
  MPI_Waitall(req_offset, ss_requests.data(), MPI_STATUSES_IGNORE);

  // find neighbor ss_elems that we have and add to our ss
  std::vector<std::vector<int>> my_ss_elems(exo->num_side_sets);
  std::vector<std::vector<int>> my_ss_sides(exo->num_side_sets);
  int total_ss_sides = 0;
  for (int ins = 0; ins < exo->num_side_sets; ins++) {
    for (int j = exo->ss_elem_index[ins]; j < exo->ss_elem_index[ins] + exo->ss_num_sides[ins];
         j++) {
      my_ss_elems[ins].push_back(exo->ss_elem_list[j]);
      my_ss_sides[ins].push_back(exo->ss_side_list[j]);
      total_ss_sides++;
    }

    for (int i = 0; i < dpi->num_neighbors; i++) {
      for (unsigned int j = 0; j < ss_neighbor_elem[i][ins].size(); j++) {
        if (elem_global_to_local.find(ss_neighbor_elem[i][ins][j]) != elem_global_to_local.end()) {
          my_ss_elems[ins].push_back(elem_global_to_local[ss_neighbor_elem[i][ins][j]]);
          my_ss_sides[ins].push_back(ss_neighbor_sides[i][ins][j]);
          total_ss_sides++;
        }
      }
    }
  }

  // resetup elem_list, side_list, indices
  int_ptr = (int *)realloc(exo->ss_elem_list, sizeof(int) * total_ss_sides);
  GOMA_ASSERT(int_ptr != NULL);
  exo->ss_elem_list = int_ptr;

  int_ptr = (int *)realloc(exo->ss_side_list, sizeof(int) * total_ss_sides);
  GOMA_ASSERT(int_ptr != NULL);
  exo->ss_side_list = int_ptr;

  exo->ss_elem_len = total_ss_sides;
  int ss_offset = 0;
  for (int ins = 0; ins < exo->num_side_sets; ins++) {
    exo->ss_elem_index[ins] = ss_offset;
    for (unsigned int j = 0; j < my_ss_elems[ins].size(); j++) {
      exo->ss_elem_list[ss_offset] = my_ss_elems[ins][j];
      exo->ss_side_list[ss_offset] = my_ss_sides[ins][j];
      ss_offset++;
    }
    exo->ss_num_sides[ins] = ss_offset - exo->ss_elem_index[ins];
  }

  // find nodes for sidesets
  std::vector<std::vector<int>> ss_nodes(exo->num_side_sets);
  for (int ins = 0; ins < exo->num_side_sets; ins++) {
    free(exo->ss_node_side_index[ins]);
    free(exo->ss_node_cnt_list[ins]);
    exo->ss_node_side_index[ins] = (int *)malloc((exo->ss_num_sides[ins] + 1) * sizeof(int));
    if (exo->ss_num_sides[ins] > 0) {
      exo->ss_node_cnt_list[ins] = (int *)malloc((exo->ss_num_sides[ins]) * sizeof(int));
    }
    for (int j = exo->ss_elem_index[ins]; j < exo->ss_elem_index[ins] + exo->ss_num_sides[ins];
         j++) {
      int elem = exo->ss_elem_list[j];
      int side = exo->ss_side_list[j];
      int block = exo->elem_eb[elem];

      int elem_type = get_type(exo->eb_elem_type[block], exo->eb_num_nodes_per_elem[block],
                               exo->eb_num_attr[block]);
      int local_side_node_list[MAX_NODES_PER_SIDE];

      /* find SIDE info for primary side */
      int num_nodes_on_side = 0;
      get_side_info(elem_type, side, &num_nodes_on_side, local_side_node_list);

      // conn location;
      int elem_offset = 0;
      for (int i = 0; i < block; i++) {
        elem_offset += exo->eb_num_elems[i];
      }

      for (int i = 0; i < num_nodes_on_side; i++) {
        int node = exo->eb_conn[block][(elem - elem_offset) * exo->eb_num_nodes_per_elem[block] +
                                       local_side_node_list[i]];
        ss_nodes[ins].push_back(node);
      }
      exo->ss_node_cnt_list[ins][j - exo->ss_elem_index[ins]] = num_nodes_on_side;
    }

    exo->ss_node_side_index[ins][0] = 0;
    if (exo->ss_num_sides[ins] > 0) {
      for (int j = 0; j < exo->ss_num_sides[ins]; j++) {
        exo->ss_node_side_index[ins][j + 1] =
            (exo->ss_node_side_index[ins][j] + exo->ss_node_cnt_list[ins][j]);
      }
    }
  }

  int total_ss_nodes = 0;
  for (int ins = 0; ins < exo->num_side_sets; ins++) {
    if (ss_nodes[ins].size() > 0) {
      int_ptr = (int *)realloc(exo->ss_node_list[ins], sizeof(int) * ss_nodes[ins].size());
      GOMA_ASSERT(int_ptr != NULL);
      exo->ss_node_list[ins] = int_ptr;
      for (unsigned int j = 0; j < ss_nodes[ins].size(); j++) {
        exo->ss_node_list[ins][j] = ss_nodes[ins][j];
      }
    }
    total_ss_nodes += ss_nodes[ins].size();
  }

  // setup distfacts that we use but expect to be 0
  dbl_ptr = (dbl *)realloc(exo->ss_distfact_list, sizeof(dbl) * total_ss_nodes);
  GOMA_ASSERT(dbl_ptr != NULL);
  for (int j = 0; j < total_ss_nodes; j++) {
    dbl_ptr[j] = 0.0;
  }
  exo->ss_distfact_list = dbl_ptr;
  ss_offset = 0;
  for (int ins = 0; ins < exo->num_side_sets; ins++) {
    exo->ss_distfact_index[ins] = ss_offset;
    exo->ss_num_distfacts[ins] = ss_nodes[ins].size();
    ss_offset += ss_nodes[ins].size();
  }

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

goma_error setup_ghost_to_base(Exo_DB *exo, Dpi *dpi) {
  exo->ghost_node_to_base = (int *)malloc(sizeof(int) * exo->num_nodes);
  std::unordered_map<int, int> old_map;
  for (int i = 0; i < exo->base_mesh->num_nodes; i++) {
    old_map[exo->base_mesh->node_map[i]] = i;
  }

  for (int i = 0; i < exo->num_nodes; i++) {
    if (old_map.find(dpi->node_index_global[i]) != old_map.end()) {
      exo->ghost_node_to_base[i] = old_map.at(dpi->node_index_global[i]);
    } else {
      exo->ghost_node_to_base[i] = -1;
    }
  }

  exo->eb_ghost_elem_to_base = (int **)malloc(sizeof(int *) * exo->num_elem_blocks);
  int offset_old = 0;
  int offset_new = 0;
  for (int i = 0; i < exo->num_elem_blocks; i++) {
    std::unordered_map<int, int> old_elem_map;
    for (int j = 0; j < exo->base_mesh->eb_num_elems[i]; j++) {
      old_elem_map[exo->base_mesh->elem_map[offset_old + j]] = j;
    }
    offset_old += exo->base_mesh->eb_num_elems[i];

    if (exo->eb_num_elems[i] > 0) {
      exo->eb_ghost_elem_to_base[i] = (int *)malloc(sizeof(int) * exo->eb_num_elems[i]);
    } else {
      exo->eb_ghost_elem_to_base[i] = NULL;
    }

    for (int j = 0; j < exo->eb_num_elems[i]; j++) {
      if (old_elem_map.find(dpi->elem_index_global[offset_new + j]) != old_elem_map.end()) {
        exo->eb_ghost_elem_to_base[i][j] = old_elem_map.at(dpi->elem_index_global[offset_new + j]);
      } else {
        exo->eb_ghost_elem_to_base[i][j] = -1;
      }
    }
    offset_new += exo->eb_num_elems[i];
  }
  return GOMA_SUCCESS;
}
