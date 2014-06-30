/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2014 Sandia Corporation.                                  *
*                                                                         *
* Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,  *
* the U.S. Government retains certain rights in this software.            *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
\************************************************************************/

/*
 * When building brk/fix as a library to be called from Goma, it is necessary to mangle a number of the brk/fix function 
 * names so there will not be link conflicts with identically named Goma routines.
 */

#ifdef BUILD_LIB
#define main                      _brk_
#define rd_exo                    _rd_exo_
#define copy_exo                  _copy_exo_
#define free_exo                  _free_exo_
#define init_exo_struct           _init_exo_struct_
#define free_exo_ev               _free_exo_ev_
#define free_exo_gv               _free_exo_gv_
#define free_exo_nv               _free_exo_nv_
#define alloc_exo_ev              _alloc_exo_ev_
#define alloc_exo_gv              _alloc_exo_gv_
#define alloc_init_exo_nv_indeces _alloc_init_exo_nv_indeces_
#define alloc_exo_nv              _alloc_exo_nv_
#define init_dpi_version          _init_dpi_version_
#define init_dpi_struct           _init_dpi_struct_
#define exo_dpi_clone             _exo_dpi_clone_
#define rd_dpi                    _rd_dpi_
#define getdid                    _getdid_
#define getvid                    _getvid_
#define getdim                    _getdim_
#define uni_dpi                   _uni_dpi_
#define free_dpi                  _free_dpi_
#define free_dpi_uni              _free_dpi_uni_
#define get_variable              _get_variable_
#define wr_dpi                    _wr_dpi_
#define wr_mesh_exo               _wr_mesh_exo_
#define wr_resetup_exo	          _wr_resetup_exo_
#define wr_result_exo             _wr_result_exo_
#define count_node_node_interactions	_count_node_node_interactions
#define in_list                   _in_list_
#define findex_mono               _findex_mono_
#define fence_post                _fence_post_
#define gcf                       _gcf_
#define get_node_index            _get_node_index_
#define get_internal_boundary_index _get_internal_boundary_index_
#define proc_sort                 _proc_sort_
#define proc_ident                _proc_ident_
#define isort                     _isort_
#define is_shell_element          _is_shell_element_
#define get_filename_num_procs    _get_filename_num_procs_
#define get_min_val_index         _get_min_val_index_
#define get_max_val_index         _get_max_val_index_
#define safe_malloc               _safe_malloc_
#define safe_free                 _safe_free_
#define zero_base                 _zero_base_
#define one_base                  _one_base_
#define shape2sides               _shape2sides_
#define get_element_shape         _get_element_shape_
#define build_elem_node           _build_elem_node_
#define build_node_node           _build_node_node_
#define build_node_elem           _build_node_elem_
#define build_elem_elem           _build_elem_elem_
#define build_side_node_list      _build_side_node_list_
#define sides2nodes               _sides2nodes_
#define int_intersect             _int_intersect_
#define demo_node_elem_conn       _demo_node_elem_conn_
#define assign_elem_ownership     _assign_elem_ownership_
#define build_elem_elem_xtra      _build_elem_elem_xtra_
#define build_elem_elem_dpi       _build_elem_elem_dpi_
#define pre_process               _pre_process_

#ifdef _BRK_C
int 
_brk_ (int , char **, char **);
#endif

#endif
