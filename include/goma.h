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
 
#ifndef _GOMA_H
#define _GOMA_H

/*
 * These include files should be just constant defs and new datatypes.
 */

#include "std.h"

#include "el_elm.h"
#include "el_geom.h"

#include "rf_allo.h"

#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_io_const.h"
#include "rf_io_structs.h"
#include "rf_io.h"
#include "rf_solver.h"

#include "rf_mp.h"
 
#include "rf_masks.h"
#include "rf_bc_const.h"
#include "rf_bc.h"
#include "rf_element_storage_struct.h"
#include "mm_elem_block_structs.h"
#include "rf_solver_const.h"
#include "rf_fill_const.h"
#include "rf_vars_const.h"
#include "mm_mp_const.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"

#include "mm_mp.h"
#include "mm_mp_structs.h"

#include "mm_species.h"

#include "mm_fill_jac.h"
#include "mm_interface.h"

#include "mm_post_def.h"
 
#include "mm_eh.h"

#include "exo_struct.h"
#include "dpi.h"
#include "dp_types.h"

/*
 * The files below are filled with mostly function prototype declarations.
 */

#include "bc_colloc.h"
#include "bc_curve.h"
#include "bc_dirich.h"
#include "bc_integ.h"
#include "bc_rotate.h"
#include "bc_special.h"
#include "bc_contact.h"
#include "bc_surfacedomain.h"
#include "dp_comm.h"
#include "dp_map_comm_vec.h"
#include "dp_utils.h"
#include "dp_vif.h"
#include "el_elm_info.h"
#include "el_quality.h"
#include "exo_conn.h"
#include "md_timer.h"
#include "mm_as_alloc.h"
#include "mm_augc_util.h"
#include "mm_bc.h"
#include "mm_chemkin.h"
#include "mm_fill.h"
#include "mm_fill_util.h"
#include "mm_fill_aux.h"
#include "mm_fill_fill.h"
#include "mm_fill_ls.h"
#include "mm_fill_porous.h"
#include "mm_fill_potential.h"
#include "mm_fill_pthings.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_rs.h"
#include "mm_fill_shell.h"
#include "mm_fill_solid.h"
#include "mm_fill_species.h"
#include "mm_fill_stress.h"
#include "mm_fill_terms.h"
#include "mm_flux.h"
#include "mm_input.h"
#include "mm_more_utils.h"
#include "mm_numjac.h"
#include "mm_ns_bc.h"
#include "mm_post_proc.h"
#include "mm_prob_def.h"
#include "mm_shell_util.h"
#include "mm_sol_nonlinear.h"
#include "mm_std_models.h"
#include "mm_qtensor_model.h"
#include "mm_unknown_map.h"
#include "mm_viscosity.h"
#include "mm_dil_viscosity.h"
#include "rd_dpi.h"
#include "rd_exo.h"
#include "rd_mesh.h"
#include "mm_shell_bc.h"

#include "ac_conti.h"
#include "ac_hunt.h"
#include "ac_update_parameter.h"
#include "ac_stability.h"
#include "ac_stability_util.h"
#include "ac_particles.h"

#include "loca_const.h"
#include "rf_element_storage_const.h" 


#include "rf_node_const.h"
#include "rf_pre_proc.h"
#include "rf_shape.h"
#include "rf_solve.h"
#include "rf_util.h"
#include "sl_aux.h"
#include "sl_auxutil.h"
#include "sl_eggroll.h"
#include "sl_lu.h"
#include "sl_matrix_util.h"
#include "sl_umf.h"
#include "sl_util.h"
#include "user_ac.h"
#include "user_bc.h"
#include "user_mp.h"
#include "user_mp_gen.h"
#include "user_post.h"
#include "user_pre.h"
#include "wr_dpi.h"
#include "wr_exo.h"
#include "wr_side_data.h"
#include "wr_soln.h"

#endif
