#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include "std.h"
#include "rf_io_const.h"
#include "rf_io_structs.h"
#include "rf_io.h"
#include "rf_solver.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_mp.h"
#include "rf_allo.h"
#include "rf_bc_const.h"
#include "mm_more_utils.h"

#include "rf_solver_const.h"
#include "sl_matrix_util.h"

#include "mm_qp_storage.h"

#include "sl_util_structs.h"

#include "sl_amesos_interface.h"

#include "sl_aztecoo_interface.h"

#include "sl_stratimikos_interface.h"
#define GOMA_DG_UTILS_C

#include "dg_utils.h"

#include "goma.h"

int has_discontinuous_interp(PROBLEM_DESCRIPTION_STRUCT *pd, int var, int imtrx)
{
  if (pd->mi[var] == imtrx)
    {
      switch (pd->i[imtrx][var]) {
      case I_P0:
      case I_P1:
      case I_PQ1:
      case I_PQ2:
        return TRUE;
      default:
        return FALSE;
      }
    }
  return FALSE;
}

int parallel_discontinuous_galerkin_enabled(PROBLEM_DESCRIPTION_STRUCT **pd_ptrs, UPD_STRUCT *upd)
{
  int enabled = 0;
  for (int imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    for (int mat = 0; mat < upd->Num_Mat; mat++) {
      if (has_discontinuous_interp(pd_ptrs[mat], MASS_FRACTION, imtrx)) {
        enabled = 1;
      } else if (has_discontinuous_interp(pd_ptrs[mat], POLYMER_STRESS11, imtrx)) {
        enabled = 1;
      }
    }
  }

  int global_enabled = 0;
  MPI_Allreduce(&enabled, &global_enabled, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  return global_enabled;
}

void count_neighbor_needed_elems(Exo_DB *exo, Dpi *dpi, int *neighbor_elem_count)
{
  for (int elem = exo->eb_ptr[0];  elem < exo->eb_ptr[exo->num_elem_blocks]; elem++) {
    for (int neighbor_index = exo->elem_elem_pntr[elem]; neighbor_index < exo->elem_elem_pntr[elem+1]; neighbor_index++) {
      if (exo->elem_elem_list[neighbor_index] == -1 &&
          dpi->elem_elem_list_global[neighbor_index] != -1) {
        // globally we have a neighbor but locally we don't have the info
        int proc = dpi->elem_elem_proc_global[neighbor_index];
        for (int neighbor = 0; neighbor < dpi->num_neighbors; neighbor++) {
          if (proc == dpi->neighbor[neighbor]) {
            neighbor_elem_count[neighbor] += 1;
          }
        }
      }
    }
  }
}

int proc_to_neighbor(int proc, Dpi *dpi) {
  for (int i = 0; i < dpi->num_neighbors; i++) {
    if (dpi->neighbor[i] == proc) {
      return i;
    }
  }

  return -1;
}

int all_elements_of_same_type(Dpi *dpi)
{
  int elem_type =  get_type(dpi->eb_elem_type_global[0],
      dpi->eb_num_nodes_per_elem_global[0],
      dpi->eb_num_attr_global[0]);

  for (int eb = 0; eb < dpi->num_elem_blocks_global; eb++) {
    int other_type = get_type(dpi->eb_elem_type_global[eb],
                              dpi->eb_num_nodes_per_elem_global[eb],
                              dpi->eb_num_attr_global[eb]);
    if (other_type != elem_type) {
      return FALSE;
    }
  }

  return TRUE;

}

// returns true if interp is unset as well
int parallel_all_interp_same(PROBLEM_DESCRIPTION_STRUCT **pd_ptrs,
                             UPD_STRUCT *upd,
                             int var)
{
  int root_interp_type = -1;
  int interp_not_same = 0;
  if (ProcID == 0) {
    for (int imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
      for (int mat = 0; mat < upd->Num_Mat; mat++) {
        if (root_interp_type == -1 && pd_ptrs[mat]->i[imtrx][var]) {
          root_interp_type = pd_ptrs[mat]->i[imtrx][var];
        } else if (pd_ptrs[mat]->i[imtrx][var] && pd_ptrs[mat]->i[imtrx][var] != root_interp_type) {
          root_interp_type = -2;
        }
      }
    }
  }

  MPI_Bcast(&root_interp_type, 1, MPI_INT, 0, MPI_COMM_WORLD);

  for (int imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    for (int mat = 0; mat < upd->Num_Mat; mat++) {
      if (pd_ptrs[mat]->i[imtrx][var] && pd_ptrs[mat]->i[imtrx][var] != root_interp_type) {
        interp_not_same = 1;
      }
    }
  }

  int global_max = 0;

  if (root_interp_type == 0) {
    return TRUE;
  }

  MPI_Allreduce(&interp_not_same, &global_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  return !global_max;
}

void dg_get_neighbor_needed_elems(Exo_DB *exo, Dpi *dpi, int *neighbor_send_elem_count, int *neighbor_recv_elem_count)
{
  count_neighbor_needed_elems(exo, dpi, neighbor_recv_elem_count);

  MPI_Request *recv_requests = calloc( dpi->num_neighbors, sizeof(MPI_Request));
  MPI_Request *send_requests = calloc( dpi->num_neighbors, sizeof(MPI_Request));
  MPI_Status *status = calloc(dpi->num_neighbors, sizeof(MPI_Status));

  for (int neighbor = 0; neighbor < dpi->num_neighbors; neighbor++)
  {
    MPI_Irecv(&neighbor_send_elem_count[neighbor], 1, MPI_INT, dpi->neighbor[neighbor],
              877877, MPI_COMM_WORLD, &recv_requests[neighbor]);
  }

  for (int neighbor = 0; neighbor < dpi->num_neighbors; neighbor++) {
    MPI_Isend(&neighbor_recv_elem_count[neighbor], 1, MPI_INT, dpi->neighbor[neighbor],
              877877, MPI_COMM_WORLD, &send_requests[neighbor]);
  }

  for (int neighbor = 0; neighbor < dpi->num_neighbors; neighbor++) {
    MPI_Wait(&send_requests[neighbor], &status[neighbor]);
  }

  for (int neighbor = 0; neighbor < dpi->num_neighbors; neighbor++) {
    MPI_Wait(&recv_requests[neighbor], &status[neighbor]);
  }

  free(recv_requests);
  free(send_requests);
  free(status);
}

void dg_get_global_send_elems(Dpi *dpi, int* recv_elem_index, int* global_recv_elems,
                              int *neighbor_send_elem_count, int *neighbor_recv_elem_count,
                              int* send_elem_index, int* global_send_elems)
{
  MPI_Request *recv_requests = calloc( dpi->num_neighbors, sizeof(MPI_Request));
  MPI_Request *send_requests = calloc( dpi->num_neighbors, sizeof(MPI_Request));
  MPI_Status *status = calloc(dpi->num_neighbors, sizeof(MPI_Status));

  for (int neighbor = 0; neighbor < dpi->num_neighbors; neighbor++) {
    MPI_Irecv(&(global_send_elems[send_elem_index[neighbor]]),
        neighbor_send_elem_count[neighbor], MPI_INT, dpi->neighbor[neighbor],
        877877, MPI_COMM_WORLD, &recv_requests[neighbor]);
  }

  for (int neighbor = 0; neighbor < dpi->num_neighbors; neighbor++) {
    MPI_Isend(&(global_recv_elems[recv_elem_index[neighbor]]),
        neighbor_recv_elem_count[neighbor], MPI_INT, dpi->neighbor[neighbor],
        877877, MPI_COMM_WORLD, &send_requests[neighbor]);
  }

  for (int neighbor = 0; neighbor < dpi->num_neighbors; neighbor++) {
    MPI_Wait(&send_requests[neighbor], &status[neighbor]);
  }

  for (int neighbor = 0; neighbor < dpi->num_neighbors; neighbor++) {
    MPI_Wait(&recv_requests[neighbor], &status[neighbor]);
  }

  free(recv_requests);
  free(send_requests);
  free(status);
}

int problem_has_stress_equations(PROBLEM_DESCRIPTION_STRUCT **pd_ptrs,
                                 UPD_STRUCT *upd)
{
  int stress_enabled = 0;
  for (int mat = 0; mat < upd->Num_Mat; mat++) {
    if (pd_ptrs[mat]->mi[POLYMER_STRESS11] != -1) {
      stress_enabled = 1;
      break;
    }
  }

  return stress_enabled;
}

int problem_has_constant_number_of_modes(struct Viscoelastic_Nonmodal **vn_ptrs,
                                         PROBLEM_DESCRIPTION_STRUCT **pd_ptrs,
                                         UPD_STRUCT *upd)
{
  if (!problem_has_stress_equations(pd_ptrs, upd)) {
    return TRUE;
  }

  int modes = vn_glob[0]->modes;
  for (int mat = 0; mat < upd->Num_Mat; mat++) {
    if (vn_glob[mat]->modes != modes) {
      return FALSE;
    }
  }

  return TRUE;
}

int setup_dg_neighbor_data(PROBLEM_DESCRIPTION_STRUCT **pd_ptrs,
                           UPD_STRUCT *upd,
                           Exo_DB *exo,
                           Dpi *dpi,
                           dg_neighbor_type *dg_neighbor_data)
{

  if (!all_elements_of_same_type(dpi)) {
    EH(-1, "Currently discontinuous galerkin requires all elements are of the same type");
    return -1;
  }

  if (!parallel_all_interp_same(pd_ptrs, upd, MASS_FRACTION) ||
      !parallel_all_interp_same(pd_ptrs, upd, POLYMER_STRESS11)) {
    EH(-1, "Currently parallel discontinuous galerkin assumes that interp for MASS_FRACTION/POLYMER_STRESS are the same on all materials");
  }

  int *neighbor_recv_elem_count = calloc(dpi->num_neighbors, sizeof(int));
  int *neighbor_send_elem_count = calloc(dpi->num_neighbors, sizeof(int));
  dg_get_neighbor_needed_elems(exo, dpi, neighbor_send_elem_count, neighbor_recv_elem_count);

  int * send_elem_index = calloc(dpi->num_neighbors + 1, sizeof(int));
  int * recv_elem_index = calloc(dpi->num_neighbors + 1, sizeof(int));

  recv_elem_index[0] = 0;
  send_elem_index[0] = 0;
  for (int neighbor = 0; neighbor < dpi->num_neighbors; neighbor++) {
    send_elem_index[neighbor + 1] = send_elem_index[neighbor];
    recv_elem_index[neighbor + 1] = recv_elem_index[neighbor];

    send_elem_index[neighbor + 1] += neighbor_send_elem_count[neighbor];
    recv_elem_index[neighbor + 1] += neighbor_recv_elem_count[neighbor];
  }

  int * recv_elem_count = calloc(dpi->num_neighbors, sizeof(int));
  int * global_recv_elems = calloc(recv_elem_index[dpi->num_neighbors], sizeof(int));

  for (int elem = exo->eb_ptr[0];  elem < exo->eb_ptr[exo->num_elem_blocks]; elem++) {
    for (int neighbor_index = exo->elem_elem_pntr[elem]; neighbor_index < exo->elem_elem_pntr[elem+1]; neighbor_index++) {
      if (exo->elem_elem_list[neighbor_index] == -1 &&
          dpi->elem_elem_list_global[neighbor_index] != -1) {
        // globally we have a neighbor but locally we don't have the info
        int proc = dpi->elem_elem_proc_global[neighbor_index];
        int local_neighbor = proc_to_neighbor(proc, dpi);

        // Only care if the proc is in the neighbor list as that is the only way
        // nodal dof's can be affected
        if (local_neighbor != -1) {

          global_recv_elems[recv_elem_index[local_neighbor] + recv_elem_count[local_neighbor]] = dpi->elem_elem_list_global[neighbor_index];

          recv_elem_count[local_neighbor] += 1;
        }
      }
    }
  }

  for (int local_neighbor = 0; local_neighbor < dpi->num_neighbors; local_neighbor++) {
    if (recv_elem_count[local_neighbor] >
        (recv_elem_index[local_neighbor+1] -recv_elem_index[local_neighbor])) {
      EH(-1, "Wrong recv count");
    }
  }

  free(recv_elem_count);

  int * global_send_elems = calloc(send_elem_index[dpi->num_neighbors], sizeof(int));


  dg_get_global_send_elems(dpi, recv_elem_index, global_recv_elems, neighbor_send_elem_count,
                           neighbor_recv_elem_count, send_elem_index, global_send_elems);

  int *local_send_elems = calloc(send_elem_index[dpi->num_neighbors], sizeof(int));
  memset(local_send_elems, -1, send_elem_index[dpi->num_neighbors] * sizeof(int));

  for (int elem = exo->eb_ptr[0];  elem < exo->eb_ptr[exo->num_elem_blocks]; elem++) {
    for (int idx = 0; idx < send_elem_index[dpi->num_neighbors]; idx++) {
      if (dpi->elem_index_global[elem] == global_send_elems[idx]) {
        local_send_elems[idx] = elem;
      }
    }
  }

  free(global_send_elems);

  // at this point we can assume that interp types and element types are the same globally
  int elem_type = get_type(dpi->eb_elem_type_global[0],
      dpi->eb_num_nodes_per_elem_global[0],
      dpi->eb_num_attr_global[0]);

  int var = MASS_FRACTION;
  int mass_fraction_num_dof = 0;
  // Checks were performed above to make sure using 0 indices are consistent here
  if (pd_glob[0]->mi[var] != -1) {
    for (int local_node = 0; local_node < exo->eb_num_nodes_per_elem[0]; local_node++) {
      for (int w = 0; w < pd->Num_Species; w++) {
        mass_fraction_num_dof += dof_lnode_interp_type(local_node,
                                        elem_type,
                                        pd_glob[0]->i[pd_ptrs[0]->mi[var]][var], 0);
      }
    }
  }

  var = POLYMER_STRESS11;
  int stress_num_dof = 0;

  if (!problem_has_constant_number_of_modes(vn_glob, pd_ptrs, upd)) {
    EH(-1, "Stress expected to have same number of modes across all materials");
    return -1;
  }

  // Checks were performed above to make sure using 0 indices are consistent here
  if (pd_glob[0]->mi[var] != -1) {
    for (int local_node = 0; local_node < exo->eb_num_nodes_per_elem[0]; local_node++) {
      for (int mode = 0; mode < vn_glob[0]->modes; mode++) {
        for (int v = POLYMER_STRESS11; v < POLYMER_STRESS33; v++) {
          if (pd_glob[0]->mi[v] != -1) {
            stress_num_dof += dof_lnode_interp_type(local_node,
                                                    elem_type,
                                                    pd_glob[0]->i[pd_ptrs[0]->mi[var]][var], 0);
          }
        }
      }
    }
  }

  dg_neighbor_data->neighbor_send_elem_index = send_elem_index;
  dg_neighbor_data->neighbor_recv_elem_index = recv_elem_index;
  dg_neighbor_data->local_elems_to_send = local_send_elems;
  dg_neighbor_data->mass_fraction_num_dof = mass_fraction_num_dof;
  dg_neighbor_data->stress_num_dof = stress_num_dof;
  dg_neighbor_data->global_recv_elem_list = global_recv_elems;

  dg_neighbor_data->nodal_recv_mass_fraction_values = NULL;
  dg_neighbor_data->nodal_recv_stress_values = NULL;
  dg_neighbor_data->nodal_send_mass_fraction_values = NULL;
  dg_neighbor_data->nodal_send_stress_values = NULL;

  if (mass_fraction_num_dof > 0) {
    int num_recv_elems = recv_elem_index[dpi->num_neighbors];
    int num_send_elems = send_elem_index[dpi->num_neighbors];
    dg_neighbor_data->nodal_recv_mass_fraction_values = calloc(mass_fraction_num_dof * num_recv_elems,
                                                               sizeof(double));
    dg_neighbor_data->nodal_send_mass_fraction_values = calloc(mass_fraction_num_dof * num_send_elems,
                                                               sizeof(double));
  }

  if (stress_num_dof > 0) {
    int num_recv_elems = recv_elem_index[dpi->num_neighbors];
    int num_send_elems = send_elem_index[dpi->num_neighbors];
    dg_neighbor_data->nodal_recv_stress_values = calloc(stress_num_dof * num_recv_elems,
                                                               sizeof(double));
    dg_neighbor_data->nodal_send_stress_values = calloc(stress_num_dof * num_send_elems,
                                                               sizeof(double));
  }

  dg_alloc_and_fill_coord(exo, dpi, dg_neighbor_data);


  free(neighbor_recv_elem_count);
  free(neighbor_send_elem_count);

  return 0;
}

void dg_alloc_and_fill_coord(Exo_DB *exo, Dpi *dpi, dg_neighbor_type *neighbor_data)
{
  int num_send_elem = neighbor_data->neighbor_send_elem_index[dpi->num_neighbors];
  int num_recv_elem = neighbor_data->neighbor_recv_elem_index[dpi->num_neighbors];

  int node_per_elem = dpi->eb_num_nodes_per_elem_global[0];

  int num_send_nodes = num_send_elem * node_per_elem;
  int num_recv_nodes = num_recv_elem * node_per_elem;


  int dim = exo->num_dim;
  double **send_card = calloc(dim, sizeof(double *));

  neighbor_data->nodal_recv_coord = calloc(dim, sizeof(double *));
  for (int card = 0; card < dim; card++) {
    send_card[card] = calloc(num_send_nodes, sizeof(double));
    neighbor_data->nodal_recv_coord[card] = calloc(num_recv_nodes, sizeof(double));
  }

  for (int card = 0; card < dim; card++) {
    for (int eidx = 0; eidx < num_send_elem; eidx++) {
      int local_elem = neighbor_data->local_elems_to_send[eidx];
      for (int local_node = 0; local_node < node_per_elem; local_node++) {
        int offset = eidx * node_per_elem + local_node;
        int node = exo->elem_node_list[exo->elem_node_pntr[local_elem] + local_node];
        switch (card) {
        case 0:
          send_card[card][offset] = exo->x_coord[node];
          break;
        case 1:
          send_card[card][offset] = exo->y_coord[node];
          break;
        case 2:
          send_card[card][offset] = exo->z_coord[node];
          break;
        default:
          EH(-1, "Unknown number of dimensions encountered");
          break;
        }
      }
    }
  }

  MPI_Request *recv_requests = calloc(dim * dpi->num_neighbors, sizeof(MPI_Request));
  MPI_Request *send_requests = calloc(dim * dpi->num_neighbors, sizeof(MPI_Request));
  MPI_Status *status = calloc(dim * dpi->num_neighbors, sizeof(MPI_Status));

  for (int neighbor = 0; neighbor < dpi->num_neighbors; neighbor++) {
    int offset = node_per_elem * neighbor_data->neighbor_recv_elem_index[neighbor];
    int count = node_per_elem * neighbor_data->neighbor_recv_elem_index[neighbor+1] - offset;

    for (int card = 0; card < dim; card++) {
      MPI_Irecv(&(neighbor_data->nodal_recv_coord[card][offset]),
                count, MPI_DOUBLE, dpi->neighbor[neighbor],
                800+card, MPI_COMM_WORLD,
                &(recv_requests[dpi->num_neighbors*card + neighbor]));
    }
  }


  for (int neighbor = 0; neighbor < dpi->num_neighbors; neighbor++) {

    int offset = node_per_elem * neighbor_data->neighbor_send_elem_index[neighbor];
    int count = node_per_elem * neighbor_data->neighbor_send_elem_index[neighbor+1] - offset;

    for (int card = 0; card < dim; card++) {
      MPI_Isend(&(send_card[card][offset]),
                count, MPI_DOUBLE, dpi->neighbor[neighbor],
                800+card, MPI_COMM_WORLD,
                &(send_requests[dpi->num_neighbors*card + neighbor]));
    }


  }

  MPI_Waitall(dpi->num_neighbors*dim, send_requests, status);
  MPI_Waitall(dpi->num_neighbors*dim, recv_requests, status);

  free(recv_requests);
  free(send_requests);
  free(status);
  for (int card = 0; card < dim; card++) {
    free(send_card[card]);
  }
  free(send_card);

}

void dg_fill_send_data(Exo_DB *exo,
                       Dpi *dpi,
                       dg_neighbor_type *neighbor_data,
                       UPD_STRUCT *upd,
                       double *x,
                       int imtrx,
                       NODE_INFO_STRUCT **Nodes,
                       struct Viscoelastic_Nonmodal **vn_ptrs)

{

  int var;

  var = MASS_FRACTION;
  if (upd->ep[imtrx][var] >= 0) {
    for (int index = 0; index < neighbor_data->neighbor_send_elem_index[dpi->num_neighbors]; index++) {

      int elem = neighbor_data->local_elems_to_send[index];
      int offset = index * neighbor_data->mass_fraction_num_dof;
      for (int w = 0; w < pd->Num_Species; w++) {
        for (int local_node = 0; local_node < dpi->eb_num_nodes_per_elem_global[0]; local_node++) {
          int node_number = exo->elem_node_list[exo->elem_node_pntr[elem] + local_node];
          NODE_INFO_STRUCT *node = Nodes[node_number];
          int nvdof = get_nv_ndofs_modMF(node->Nodal_Vars_Info[imtrx], var);
          for (int j = 0; j < nvdof; j++)
          {
            int ie = Index_Solution(node_number, var, w, 0, -2, imtrx) + j;
            EH(ie, "Could not find vbl in sparse matrix.");
            neighbor_data->nodal_send_mass_fraction_values[offset] = x[ie];
            offset++;
          }
        }
      }
    }
  }

  var = POLYMER_STRESS11;
  if (upd->ep[imtrx][var] >= 0) {
    for (int index = 0; index < neighbor_data->neighbor_send_elem_index[dpi->num_neighbors]; index++) {
      EH(-1, "Stress not implemented");
      int elem = neighbor_data->local_elems_to_send[index];
      int offset = elem * neighbor_data->mass_fraction_num_dof;
      for (int mode = 0; mode < vn_ptrs[0]->modes; mode++) {
        for (int local_node = 0; local_node < dpi->eb_num_nodes_per_elem_global[0]; local_node++) {
          int node_number = exo->elem_node_list[exo->elem_node_pntr[elem] + local_node];
          NODE_INFO_STRUCT *node = Nodes[node_number];
          int nvdof = get_nv_ndofs_modMF(node->Nodal_Vars_Info[imtrx], var);
          for (int j = 0; j < nvdof; j++)
          {
            int ie = Index_Solution(node_number, var, 0, 0, -2, imtrx) + j;
            EH(ie, "Could not find vbl in sparse matrix.");
            neighbor_data->nodal_send_mass_fraction_values[offset] = x[ie];
            offset++;
          }
        }
      }
    }
  }
}

void dg_communicate_neighbor_data(Exo_DB *exo,
                                  Dpi *dpi,
                                  dg_neighbor_type *neighbor_data,
                                  UPD_STRUCT *upd,
                                  double *x,
                                  int imtrx,
                                  NODE_INFO_STRUCT **Nodes,
                                  struct Viscoelastic_Nonmodal **vn_ptrs)
{

  dg_fill_send_data(exo, dpi, neighbor_data, upd, x, imtrx, Nodes, vn_ptrs);


  MPI_Request *recv_requests = calloc( 2*dpi->num_neighbors, sizeof(MPI_Request));
  MPI_Request *send_requests = calloc( 2*dpi->num_neighbors, sizeof(MPI_Request));
  MPI_Status *status = calloc(2*dpi->num_neighbors, sizeof(MPI_Status));

  for (int neighbor = 0; neighbor < dpi->num_neighbors; neighbor++) {
    int recv_elem_offset = neighbor_data->neighbor_recv_elem_index[neighbor];
    int count =  neighbor_data->neighbor_recv_elem_index[neighbor+1] - recv_elem_offset;
    if ((upd->ep[imtrx][MASS_FRACTION] >= 0) && (neighbor_data->mass_fraction_num_dof > 0)) {
      int mass_frac_count = count * neighbor_data->mass_fraction_num_dof;
      int mass_frac_offset = recv_elem_offset * neighbor_data->mass_fraction_num_dof;
      MPI_Irecv(&(neighbor_data->nodal_recv_mass_fraction_values[mass_frac_offset]),
                mass_frac_count, MPI_DOUBLE, dpi->neighbor[neighbor],
                888, MPI_COMM_WORLD, &recv_requests[neighbor]);
    }


    if ((upd->ep[imtrx][POLYMER_STRESS11] >= 0) && (neighbor_data->stress_num_dof > 0)) {
      printf("Sanity check\n");
      int stress_count = count * neighbor_data->stress_num_dof;
      int stress_offset = recv_elem_offset * neighbor_data->stress_num_dof;
      MPI_Irecv(&(neighbor_data->nodal_recv_stress_values[stress_offset]),
                stress_count, MPI_DOUBLE, dpi->neighbor[neighbor],
                877, MPI_COMM_WORLD, &recv_requests[dpi->num_neighbors + neighbor]);
    }
  }


  for (int neighbor = 0; neighbor < dpi->num_neighbors; neighbor++) {
    int send_elem_offset = neighbor_data->neighbor_send_elem_index[neighbor];
    int count = neighbor_data->neighbor_send_elem_index[neighbor + 1] - send_elem_offset;
    if ((upd->ep[imtrx][MASS_FRACTION] >= 0) && (neighbor_data->mass_fraction_num_dof > 0)) {
      int mass_frac_count = count * neighbor_data->mass_fraction_num_dof;
      int mass_frac_offset = send_elem_offset * neighbor_data->mass_fraction_num_dof;
      MPI_Isend(&(neighbor_data->nodal_send_mass_fraction_values[mass_frac_offset]),
                mass_frac_count,
                MPI_DOUBLE, dpi->neighbor[neighbor], 888,
                MPI_COMM_WORLD, &send_requests[neighbor]);
    }

    if (neighbor_data->stress_num_dof > 0) {
      int stress_count = count * neighbor_data->stress_num_dof;
      int stress_offset = send_elem_offset * neighbor_data->stress_num_dof;
      MPI_Isend(&(neighbor_data->nodal_send_stress_values[stress_offset]),
                stress_count,
                MPI_DOUBLE, dpi->neighbor[neighbor], 877,
                MPI_COMM_WORLD, &send_requests[dpi->num_neighbors + neighbor]);
    }
  }

  for (int neighbor = 0; neighbor < dpi->num_neighbors; neighbor++) {
    if ((upd->ep[imtrx][MASS_FRACTION] >= 0) && (neighbor_data->mass_fraction_num_dof > 0)) {
       MPI_Wait(&send_requests[neighbor], &status[neighbor]);
    }

    if (neighbor_data->stress_num_dof > 0) {
       MPI_Wait(&send_requests[dpi->num_neighbors + neighbor], &status[dpi->num_neighbors + neighbor]);
    }

  }

  for (int neighbor = 0; neighbor < dpi->num_neighbors; neighbor++) {
    if ((upd->ep[imtrx][MASS_FRACTION] >= 0) && (neighbor_data->mass_fraction_num_dof > 0)) {
       MPI_Wait(&recv_requests[neighbor], &status[neighbor]);
    }

    if ((upd->ep[imtrx][POLYMER_STRESS11] >= 0) && (neighbor_data->stress_num_dof > 0)) {
       MPI_Wait(&recv_requests[dpi->num_neighbors + neighbor], &status[dpi->num_neighbors + neighbor]);
    }
  }

  free(recv_requests);
  free(send_requests);
  free(status);

}
