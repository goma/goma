#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <memory>
#include <nanoflann.hpp>
#include <unordered_set>
#include <vector>

#include "util/distance_helpers.h"
extern "C" {
#define DISABLE_CPP
#include "mm_as.h"
#include "mm_eh.h"
#include "mm_unknown_map.h"
#include "rf_fem_const.h"
#undef DISABLE_CPP
}

using coordinates_type = std::vector<std::array<double, 3>>;

/** A simple vector-of-vectors adaptor for nanoflann, without duplicating the
 * storage. The i'th vector represents a point in the state space.
 *
 *  \tparam DIM If set to >0, it specifies a compile-time fixed dimensionality
 *      for the points in the data set, allowing more compiler optimizations.
 *  \tparam num_t The type of the point coordinates (typ. double or float).
 *  \tparam Distance The distance metric to use: nanoflann::metric_L1,
 *          nanoflann::metric_L2, nanoflann::metric_L2_Simple, etc.
 *  \tparam IndexType The type for indices in the KD-tree index
 *         (typically, size_t of int)
 */
template <class VectorOfArraysType,
          typename NumberType = double,
          int KDTreeDim = -1,
          class Distance = nanoflann::metric_L2,
          typename IndexType = size_t>
struct KDTreeVectorOfArraysAdaptor {
  using SelfType = KDTreeVectorOfArraysAdaptor<VectorOfArraysType, NumberType, KDTreeDim, Distance>;
  using MetricType = typename Distance::template traits<NumberType, SelfType>::distance_t;
  using KDTreeIndexType =
      nanoflann::KDTreeSingleIndexAdaptor<MetricType, SelfType, KDTreeDim, IndexType>;

  /** The kd-tree index for the user to call its methods as usual with any
   * other FLANN index */
  std::unique_ptr<KDTreeIndexType> index;

  /// Constructor: takes a const ref to the vector of vectors object with the
  /// data points
  KDTreeVectorOfArraysAdaptor(const size_t,
                              const VectorOfArraysType &mat,
                              const size_t dim,
                              const int leaf_max_size = 10)
      : vector_of_arrays_(mat) {
    assert(mat.size() != 0 && mat[0].size() != 0);
    const size_t dims = dim;
    if (KDTreeDim > 0 && static_cast<int>(dims) != KDTreeDim)
      throw std::runtime_error("Data set dimensionality does not match the 'kd_tree_dim' template "
                               "argument");
    index = std::unique_ptr<KDTreeIndexType>(
        new KDTreeIndexType(static_cast<int>(dims), *this /* adaptor */,
                            nanoflann::KDTreeSingleIndexAdaptorParams(leaf_max_size)));
  }

  const VectorOfArraysType &vector_of_arrays_;

  /** Query for the \a num_closest closest points to a given point
   *  (entered as query_point[0:dim-1]).
   *  Note that this is a short-cut method for index->findNeighbors().
   *  The user can also call index->... methods as desired.
   *
   * \note nChecks_IGNORED is ignored but kept for compatibility with
   * the original FLANN interface.
   */
  inline void query(const NumberType *query_point,
                    const size_t num_closest,
                    IndexType *out_indices,
                    NumberType *out_distances_sq,
                    const int nChecks_IGNORED = 10) const {
    nanoflann::KNNResultSet<NumberType, IndexType> resultSet(num_closest);
    resultSet.init(out_indices, out_distances_sq);
    index->findNeighbors(resultSet, query_point, nanoflann::SearchParams());
  }

  /** @name Interface expected by KDTreeSingleIndexAdaptor
   * @{ */

  const SelfType &derived() const { return *this; }
  SelfType &derived() { return *this; }

  // Must return the number of data points
  inline size_t kdtree_get_point_count() const { return vector_of_arrays_.size(); }

  // Returns the dim'th component of the idx'th point in the class:
  inline NumberType kdtree_get_pt(const size_t idx, const size_t dim) const {
    return vector_of_arrays_[idx][dim];
  }

  // return false to use default
  template <class BBOX> bool kdtree_get_bbox(BBOX & /*bb*/) const { return false; }

  /** @} */
};

template <int dim>
void find_distances(coordinates_type side_coordinates,
                    coordinates_type distance_coordinates,
                    std::vector<double> &distances) {
  KDTreeVectorOfArraysAdaptor<coordinates_type, double, dim> kd_tree(side_coordinates.size(),
                                                                     side_coordinates, dim);

  for (size_t i = 0; i < distance_coordinates.size(); i++) {
    std::vector<size_t> ret_indexes(dim);
    std::vector<double> out_dists_sqr(dim);
    kd_tree.query(&distance_coordinates[i][0], 1, &ret_indexes[0], &out_dists_sqr[0]);
    distances[i] = sqrt(out_dists_sqr[0]);
  }
}

extern "C" goma_error find_current_distances(Exo_DB *exo,
                                             Dpi *dpi,
                                             double *solution_vector,
                                             bool apply_displacements,
                                             int num_ns,
                                             int *ns_ids,
                                             int num_ss,
                                             int *ss_ids,
                                             double *distances) {

  const int dim = exo->num_dim;

  if (dim != 2 && dim != 3) {
    return GOMA_ERROR;
  }

  // get the coordinates of all the nodes
  coordinates_type coordinates(exo->num_nodes);
  for (int i = 0; i < exo->num_nodes; ++i) {
    coordinates[i][0] = exo->x_coord[i];
    coordinates[i][1] = exo->y_coord[i];
    if (dim == 3) {
      coordinates[i][2] = exo->z_coord[i];
    } else {
      coordinates[i][2] = 0.0;
    }

    // apply displacements to the nodes
    if (apply_displacements) {
      int index =
          Index_Solution(i, MESH_DISPLACEMENT1, 0, 0, -1, upd->matrix_index[MESH_DISPLACEMENT1]);
      GOMA_ASSERT_ALWAYS(index != -1);
      coordinates[i][0] += solution_vector[index];
      index =
          Index_Solution(i, MESH_DISPLACEMENT2, 0, 0, -1, upd->matrix_index[MESH_DISPLACEMENT2]);
      GOMA_ASSERT_ALWAYS(index != -1);
      coordinates[i][1] += solution_vector[index];
      if (dim == 3) {
        index =
            Index_Solution(i, MESH_DISPLACEMENT3, 0, 0, -1, upd->matrix_index[MESH_DISPLACEMENT3]);
        GOMA_ASSERT_ALWAYS(index != -1);
        coordinates[i][2] += solution_vector[index];
      }
    }
  }

  std::unordered_set<int> set_nodes;

  // Get nodes from sidesets
  for (int in_ss_index = 0; in_ss_index < num_ss; in_ss_index++) {
    int ss_id = ss_ids[in_ss_index];
    int ss_index = -1;
    for (int i = 0; i < exo->num_side_sets; i++) {
      if (exo->ss_id[i] == ss_id) {
        ss_index = i;
        break;
      }
    }
    if (ss_index == -1) {
      GOMA_EH(GOMA_ERROR, "Side set %d not found", ss_id);
    }
    for (int side = 0; side < exo->ss_num_sides[ss_index]; side++) {
      for (int lni = exo->ss_node_side_index[ss_index][side];
           lni < exo->ss_node_side_index[ss_index][side + 1]; lni++) {
        int inode = exo->ss_node_list[ss_index][lni];
        if (dpi->num_proc == 0 || dpi->node_owner[inode] == dpi->rank) {
          set_nodes.insert(inode);
        }
      }
    }
  }

  // Get nodes from nodesets
  for (int in_ns_index = 0; in_ns_index < num_ns; in_ns_index++) {
    int ns_id = ns_ids[in_ns_index];
    int ns_index = -1;
    for (int i = 0; i < exo->num_node_sets; i++) {
      if (exo->ns_id[i] == ns_id) {
        ns_index = i;
        break;
      }
    }
    if (ns_index == -1) {
      GOMA_EH(GOMA_ERROR, "Node set %d not found", ns_id);
    }
    for (int lni = 0; lni < exo->ns_num_nodes[ns_index]; lni++) {
      int inode = exo->ns_node_list[exo->ns_node_index[ns_index] + lni];
      if (dpi->num_proc == 0 || dpi->node_owner[inode] == dpi->rank) {
        set_nodes.insert(inode);
      }
    }
  }

  std::vector<int> local_side_nodes(set_nodes.begin(), set_nodes.end());

  coordinates_type side_coordinates;
  if (dpi->num_proc > 0) {
    std::vector<size_t> proc_node_counts(dpi->num_proc);
    size_t my_count = local_side_nodes.size();
    size_t max_count = 0;
    std::vector<size_t> global_counts(dpi->num_proc);
    MPI_Allreduce(&my_count, &max_count, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allgather(&my_count, 1, MPI_UINT64_T, &proc_node_counts[0], 1, MPI_UINT64_T,
                  MPI_COMM_WORLD);
    size_t total_count = 0;
    for (size_t i = 0; i < proc_node_counts.size(); i++) {
      total_count += proc_node_counts[i];
    }
    std::vector<double> global_coordinates(dpi->num_proc * max_count * dim);
    std::vector<double> my_side_coordinates(max_count * dim);
    std::fill(my_side_coordinates.begin(), my_side_coordinates.end(), NAN);
    for (size_t i = 0; i < local_side_nodes.size(); i++) {
      for (int j = 0; j < dim; j++) {
        my_side_coordinates[i * dim + j] = coordinates[local_side_nodes[i]][j];
      }
    }
    MPI_Allgather(my_side_coordinates.data(), max_count * dim, MPI_DOUBLE,
                  global_coordinates.data(), max_count * dim, MPI_DOUBLE, MPI_COMM_WORLD);

    side_coordinates.resize(total_count);
    size_t offset_local = 0;
    for (int i = 0; i < dpi->num_proc; i++) {
      size_t offset = i * max_count * dim;
      for (size_t j = 0; j < proc_node_counts[i]; j++) {
        for (int k = 0; k < dim; k++) {
          side_coordinates[offset_local][k] = global_coordinates[offset + j * dim + k];
        }
        offset_local++;
      }
    }
  } else {
    side_coordinates.resize(local_side_nodes.size());
    for (size_t i = 0; i < local_side_nodes.size(); i++) {
      for (int j = 0; j < dim; j++) {
        side_coordinates[i][j] = coordinates[local_side_nodes[i]][j];
      }
    }
  }

  std::vector<double> distance_vector(exo->num_nodes);
  if (dim == 3) {
    find_distances<3>(side_coordinates, coordinates, distance_vector);
  } else if (dim == 2) {
    find_distances<2>(side_coordinates, coordinates, distance_vector);
  } else {
    fprintf(stderr, "%s ERROR: Unsupported dimensionality: %d\n", __FUNCTION__, dim);
    return GOMA_ERROR;
  }

  for (size_t i = 0; i < distance_vector.size(); i++) {
    distances[i] = distance_vector[i];
  }

  return GOMA_SUCCESS;
}