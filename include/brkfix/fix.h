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
#ifndef FIX_H
#define FIX_H

struct fix_data {
  int ns_node_len_global;
  int ns_distfact_len_global;
  int ss_elem_len_global;
  int ss_distfact_len_global;

  char **eb_elem_type_global;
  int *eb_index_global;

  // int *eb_num_elems_global;
  // int *ns_num_nodes_global;
  // int *ss_num_sides_global;
  // int *ss_num_distfacts_global;
  // int *ns_num_distfacts_global;
  int *eb_num_nodes_per_elem;
  int *ns_distfact_list_index_global;
  int *ss_distfact_list_index_global;
  int *ns_node_index_global;
  int *ns_distfact_index_global;
  int *ss_elem_index_global;
  int *ss_distfact_index_global;
  int *ns_node_list_index_global;
  int *ss_elem_list_index_global;
  int num_elem_blocks;

  int *ns_node_list;
  double *ns_distfact_list;
  int *ss_elem_list;
  int *ss_side_list;
  double *ss_distfact_list;
};

void fix_output(void);
int fix_exo_file(int num_procs, const char *exo_mono_name);

#endif /* FIX_H */
