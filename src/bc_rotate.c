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
 *$Id: bc_rotate.c,v 5.4 2008-10-07 14:33:41 hkmoffa Exp $
 */

/* Standard include files */
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "std.h"
#include "rf_allo.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_io_const.h"
#include "rf_io.h"
#include "el_elm.h"
#include "el_geom.h"
#include "exo_struct.h"
#include "rf_bc_const.h"
#include "rf_bc.h"
#include "rf_solver.h"
#include "rf_vars_const.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_eh.h"
#include "rotate_util.h"
#include "ac_stability.h"
#include "ac_stability_util.h"
#include "bc_rotate.h"
#include "el_elm_info.h"
#include "gds/gds_vector.h"
#include "mm_as_alloc.h"
#include "mm_bc.h"
#include "mm_elem_block_structs.h"
#include "mm_fill_aux.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_terms.h"
#include "mm_fill_util.h"
#include "mm_post_proc.h"
#include "mm_unknown_map.h"
#include "rd_mesh.h"
#include "rf_node_const.h"
#include "sl_util_structs.h"

/*
 *  Variable Definitions
 */
struct Rotation_Vectors ****rotation = NULL;
int dup_blks_list[MAX_MAT_PER_SS+1];

#define GOMA_BC_ROTATE_C
#include "sl_epetra_interface.h"


/*********** R O U T I N E S   I N   T H I S   F I L E *************************
 *
 *       NAME            TYPE            CALLED_BY
 *    ------------             ---------               --------------
 *
 *  apply_rotated_bc  ()         int         matrix_fill
 *  rotate_res_jac_mesh  ()      void        rf_fill.c:  matrix_fill: apply_rotated_bc
 *  rotate_res_jac_mom   ()      void        rf_fill.c:  matrix_fill: apply_rotated_bc
 *
 *******************************************************************************/
/*ARGSUSED*/
int
apply_rotated_bc (
     double resid_vector[],   /* Residual vector for the current processor         */
     struct elem_side_bc_struct *first_elem_side_BC_array[],
                              /* An array of pointers to the first surface integral defined
                                 for each element. Its length is equal to the total number
                                 of elements defined on the current processor.     */
     const int ielem,               /* element number                                    */
     const int ielem_type,          /* element type                                      */
     const int num_local_nodes,
     const int ielem_dim,
     const int iconnect_ptr,
     const int num_total_nodes,
     const Exo_DB *exo)
    
    /****************************************************************************
     *  
     *  apply_rotated_bc
     *
     *      This routine handles rotation of mesh and momentum equations
     *      for 2D problems.
     ****************************************************************************/
{
#ifdef DEBUG
  static char *yo = "apply_rotated_bc";
  int k;
#endif
  char err_msg[MAX_CHAR_IN_INPUT]; 
  int i, I, ibc, j, id, blk_index, ss_index; /* counters */
  int err;                    /* status variable for functions */
  int status = 0;
  int bc_input_id;

  struct elem_side_bc_struct *elem_side_bc;
  int rotate_this_side;
  int rot_mesh[MAX_NODES_PER_SIDE];
  int rot_mom[MAX_NODES_PER_SIDE];
  int mesh_already_rotated[MDE];
  int mom_already_rotated[MDE];
  double sign_rot_mesh[MAX_NODES_PER_SIDE];
  double sign_rot_mom[MAX_NODES_PER_SIDE];
  int ss_rot, irot, num_nodes_on_side;
  double s, t;			/* Gaussian-quadrature point locations */
  double xi[DIM];		/* Local element coordinates of Gauss point. */

  /*****************************************************************************/

  /* set up arrays to avoid double-rotation of equations */
  for (id = 0; id < num_local_nodes; id++) {
    mesh_already_rotated[id] = 0;
    mom_already_rotated[id]  = 0;
  }
  
  /*
   *  begining of do while construct which loops over the sides of this
   *  element that have boundary conditions 
   *  - determine which side to rotate wrt 
   *    look at each side of the element and determine if rotation is
   *    needed w.r.t. this side 
   */
  elem_side_bc = first_elem_side_BC_array[ielem];
  do {
    num_nodes_on_side = (int) elem_side_bc->num_nodes_on_side;
    rotate_this_side = 0;
    for (i = 0; i < num_nodes_on_side; i++) {
      rot_mesh[i] = 0;
      rot_mom[i]  = 0;
      sign_rot_mesh[i] = 1.0;
      sign_rot_mom[i]  = 1.0;
    }
    for (i = 0; i < num_nodes_on_side; i++) {
      id = (int) elem_side_bc->local_elem_node_id[i];
      I = Proc_Elem_Connect[iconnect_ptr + id]; /* Find global node */
      /* check to see if this global node is in the momentum rotation list and
       * make sure this nodal equation hasn't been rotated yet */
      if (((irot = in_list(I, 0, num_mom_rotate[pg->imtrx], mom_rotate_node[pg->imtrx])) != -1) && 
	  mom_already_rotated[id] == 0) {
	/* determine if SS for rotation is on this side */
	ss_rot = mom_rotate_ss[pg->imtrx][irot];
	for (ibc = 0; 
	     (bc_input_id = (int) elem_side_bc->BC_input_id[ibc]) != -1; ibc++) {
	    if (BC_Types[bc_input_id].BC_ID == ss_rot) {
	      rot_mom[i] = 1;
	      rotate_this_side = 1;
	      mom_already_rotated[id]  = 1;
	      /*
	       *  Determine if the tangent should be negated for a consistent
	       * force balance 
               */
	      if ((ss_index = 
		   in_list(ss_rot, 0, Proc_Num_Side_Sets, ss_to_blks[0])) == -1) {
		EH(-1,"Cannot match side set id with that in ss_to_blks array");
	      } else {
		for (j = 0; j < MAX_MAT_PER_SS; j++) {
		  dup_blks_list[j] = ss_to_blks[j+1][ss_index];
		}
	      }
	      if ((blk_index = in_list(Current_EB_ptr->Elem_Blk_Id, 0,
				       MAX_MAT_PER_SS+1,
				       dup_blks_list)) == -1) {
		EH(-1,"Cannot match current element block with that in ss_to_blks array");
	      }
	      blk_index++;
	      if (SS_Internal_Boundary[ss_index] != -1 && ((blk_index % 2) == 1)) {
		sign_rot_mom[i] = -1.;
	      }
	    }
	  }
      }
      
      if (((irot = in_list(I, 0, num_mesh_rotate[pg->imtrx], mesh_rotate_node[pg->imtrx])) != -1) && 
	  mesh_already_rotated[id] == 0) {
	/* determine if SS for rotation is on this side */
	ss_rot = mesh_rotate_ss[pg->imtrx][irot];
        /* check to see if this global node is in the mesh rotation list and
         * make sure this nodal equation hasn't been rotated yet */
	for (ibc = 0; (bc_input_id = (int) elem_side_bc->BC_input_id[ibc]) != -1; ibc++)
	  {
	    if (BC_Types[bc_input_id].BC_ID == ss_rot) {
	      rot_mesh[i] = 1;
	      rotate_this_side = 1;
	      mesh_already_rotated[id]  = 1;
	      /*  ---- Determine if the tangent should be negated for a consistent force balance */
	      if ((ss_index = 
		   in_list(ss_rot, 0, exo->num_side_sets, ss_to_blks[0])) == -1) {
		sprintf(err_msg, "Could not find SS %d in ss_to_blks",
			ss_rot);
		EH(-1, err_msg);
	      } else {
		for (j=0; j<MAX_MAT_PER_SS; j++) {
		  dup_blks_list[j] = ss_to_blks[j+1][ss_index];
		}
	      }
	      if ((blk_index = in_list(Current_EB_ptr->Elem_Blk_Id, 0,
				       MAX_MAT_PER_SS+1, dup_blks_list)) == -1) {
#ifdef DEBUG
		fprintf(stderr, 
			"%s: Found ssindex=%d, SSID=%d, dup_blks_list = ",
			yo, ss_index, exo->ss_id[ss_index]);
		for ( k=0; k<MAX_MAT_PER_SS; k++)
		{
		  fprintf(stderr, "%d ", dup_blks_list[k]);
		}
		fprintf(stderr, "\n");
#endif
		sprintf(err_msg, "Could not find EB %d in dup_blks_list",
			Current_EB_ptr->Elem_Blk_Id);
		EH(-1, err_msg);
	      }
	      blk_index++;
	      if( SS_Internal_Boundary[ss_index] != -1 && ((blk_index % 2) == 1) ) {
		sign_rot_mesh[i] = -1.;
	      }
	      /* Special cludge for Bridgman case (the Pi case)  */
	      /*if (blk_index == 2 &&  BC_Types[bc_input_id].BC_Name == PLANE_BC
               *    && BC_Types[bc_input_id].BC_ID == 1013) {
	       *   sign_rot_mesh[i] = -1.;
	       *}
               */
	    }
	  }
      }
    }

    if (rotate_this_side) {
      /*  Evaluate n-t for rotation. Chose the centroid of the element surface */
      s = 0.0 ; t = 0.0 ;
      xi[0] = 0.0; xi[1] = 0.0; xi[2] = 0.0; /* need to update for 3D ??*/
      
      /* find the quadrature point locations for current ip */
      find_surf_center_st (ielem_type, elem_side_bc->id_side, pd->Num_Dim, xi, &s, &t);
      err = load_basis_functions(xi, bfd);
      EH( err, "problem from load_basis_functions");
      
      err = beer_belly();
      EH( err, "beer_belly");
	  
/*	  err = load_coordinate_scales(pd->CoordinateSystem, fv);
	  EH(err, "load_coordinate_scales(fv)");
*/      
      /* calculate the surface determinant of the surface jacobian */
      surface_determinant_and_normal(ielem, iconnect_ptr, num_local_nodes, 
				     ielem_dim - 1, 
				     (int) elem_side_bc->id_side,
				     (int) elem_side_bc->num_nodes_on_side,
				     (elem_side_bc->local_elem_node_id) );
      
      /* calculate the components of the surface normal and mesh displacement
       * derivatives*/
      
      if (ielem_dim !=3) {
	calc_surf_tangent(ielem, iconnect_ptr, num_local_nodes, ielem_dim-1,
			  (int) elem_side_bc->num_nodes_on_side,
			  (elem_side_bc->local_elem_node_id));      
      } else {
	EH(-1,"Illegal dimension in old rotation scheme");
      }
	
      do_LSA_mods(LSA_SURFACE);

      /* 
       * Loop over the number of local element nodes on this side 
       * and rotate the equations
       */
      for (i = 0; i < (int) elem_side_bc->num_nodes_on_side; i++) {
	
	/* Find the local element node number for the current node */
	id = (int) elem_side_bc->local_elem_node_id[i];
	
	/* Find the local node number given the local element node number,  'i'     */
	I = Proc_Elem_Connect[iconnect_ptr + id];
	
	if (rot_mesh[i] &&
	    ei[pg->imtrx]->ln_to_dof[MESH_DISPLACEMENT1][id] != -1) {
	  rotate_res_jac_mesh(ei[pg->imtrx]->ln_to_dof[MESH_DISPLACEMENT1][id],
			      I, iconnect_ptr,
			      num_local_nodes, ielem_dim-1,
			      fv->snormal, fv->dsnormal_dx, 
			      fv->stangent, fv->dstangent_dx, sign_rot_mesh[i]);
	}
	if (rot_mom[i] && ei[pg->imtrx]->ln_to_dof[VELOCITY1][id] != -1 && 
	    pd->i[pg->imtrx][VELOCITY1] != I_Q2_D &&
            pd->i[pg->imtrx][VELOCITY1] != I_Q1_D)	  { 
		rotate_res_jac_mom(ei[pg->imtrx]->ln_to_dof[VELOCITY1][id], 
			     I, iconnect_ptr, num_local_nodes, ielem_dim-1,
			     fv->snormal, fv->dsnormal_dx, fv->stangent,
			     fv->dstangent_dx, sign_rot_mom[i]); 
	  }
	else if (rot_mom[i] && ei[pg->imtrx]->ln_to_dof[VELOCITY1][id] != -1)
	  {
	    if (((Current_EB_ptr->Elem_Blk_Id + 1) % 2) == 0)
	      {
		rotate_res_jac_mom(ei[pg->imtrx]->ln_to_dof[VELOCITY1][id] - 1, 
				   I, iconnect_ptr,
				   num_local_nodes, ielem_dim-1,
				   fv->snormal, fv->dsnormal_dx, fv->stangent, fv->dstangent_dx, 
				   sign_rot_mom[i]); 
	      }
	    else
	      {
		rotate_res_jac_mom(ei[pg->imtrx]->ln_to_dof[VELOCITY1][id] - 1, 
				   I, iconnect_ptr,
				   num_local_nodes, ielem_dim-1,
				   fv->snormal, fv->dsnormal_dx, fv->stangent, fv->dstangent_dx, 
				   sign_rot_mom[i]);
	      } 
	  }
      } /* end of loop over boundary nodes */
    }
  } while ( (elem_side_bc = elem_side_bc->next_side_bc) != NULL );
  /* End of do  while () construct				      */
  return(status);

} /* END of routine apply_rotated_bc */
/*********************************************************************************/
/*********************************************************************************/
/*********************************************************************************/
/*ARGSUSED*/
void
rotate_res_jac_mesh (
    int irow_index,                     /* Elemental stiffness matrix row index */
    int I,                              /* Global node number                   */
    int iconnect_ptr,  		        /* Pointer to beginning of connectivity list
				           for the current element	         */
    int nodes_per_elem,                 /* number of nodes in the element       */
    int ielem_surf_dim,                 /* physical dim of elem surface - 0,1,2 */
    double snormal[MAX_PDIM],           /* vector components of surface normal  */
    double dsnormal_dx[MAX_PDIM][MAX_PDIM][MDE],
                                        /* Array of derivatives ([i][j][k]) of
                                           surface normal: component-i
                                           wrt displacement-j at node-k         */
    double stangent[2][MAX_PDIM],       /* vector components of 2 mutually orthogonal
					   surface tangent vectors            */
    double dstangent_dx[2][MAX_PDIM][MAX_PDIM][MDE],
                                        /* Array of derivatives ([i][j][k]) of surface
                                           tangent vector 1 or 2: component-i
                                           wrt displacement-j at node-k         */
    double sign )                       /* sign of tangent vector               */

/* 
 * Function which corrects the global residual vector "resid_vect" and
 * the global Jacobian vector "a" so that the vector mesh equations for surface nodes 
 * are projected into a normal and tangential coordinate system.
 *
 *       Author:          P. R. Schunk (1511)
 *       Revisor:         R. R. Rao  20 July 1999 (lec revolution)
 */

{
  int       var,  peq, pvar;
  int 	    id, j, n, Id, ldof, ShapeVar;
  unsigned int siz;
  int       kdir, ldir, w, jvar;
  double    svector[MAX_PDIM][MAX_PDIM];
              /* set of vectors corresponding to normal (1st vector) and tangents
                 at the surface - for convenience in the following loops */
  double    dsvector_dx[MAX_PDIM][MAX_PDIM][MAX_PDIM][MDE];
              /* sensitivity of surface vectors with respect to nodal positions */

  double    rotated_resid[MDE];
  double    rotated_jacobian_vector[MAX_PDIM][MAX_PDIM][MDE];
  double    rotated_jacobian_scalar[MAX_PDIM][MDE];
  double    rotated_jacobian_conc[MAX_PDIM][MAX_CONC][MDE];

  /* Load up values of surface vectors */
  /* first surface vector is the normal 
   *    - i.e. first mesh equation is the normal component */
  /* second (and third) surface vector(s) is(are) the tangent(s)
   *    - i.e. second and third mesh equations are the tangential components */

  ShapeVar = pd->ShapeVar;
  
  switch (ielem_surf_dim+1) 
    {
      
    case 1:
      svector[0][0] = snormal[0];
      for ( id=0; id<ei[pg->imtrx]->dof[ShapeVar]; id++ ) 
        {
	  Id = Proc_Elem_Connect[iconnect_ptr + id];
	  ldof = ei[pg->imtrx]->ln_to_dof[ShapeVar][id];
	  if (Dolphin[pg->imtrx][I][MESH_DISPLACEMENT1] > 0 && Dolphin[pg->imtrx][Id][MESH_DISPLACEMENT1] > 0 )
	    { 
	      dsvector_dx[0][0][0][ldof] = dsnormal_dx[0][0][ldof];
	    }
	}
      break;
      
    case 2:
      for (kdir = 0; kdir < ielem_surf_dim+1; kdir++)
	{
	  svector[0][kdir] = snormal[kdir];
	  svector[1][kdir] = sign * stangent[0][kdir];
	  for (ldir = 0; ldir < ielem_surf_dim+1; ldir++)
	    {
	      for ( id=0; id<ei[pg->imtrx]->dof[ShapeVar]; id++ ) 
		{
		  Id = Proc_Elem_Connect[iconnect_ptr + id];
		  ldof = ei[pg->imtrx]->ln_to_dof[ShapeVar][id];
		  if (Dolphin[pg->imtrx][I][MESH_DISPLACEMENT1] > 0  && Dolphin[pg->imtrx][Id][MESH_DISPLACEMENT1] > 0 )
		    { 
		      dsvector_dx[0][kdir][ldir][ldof] = dsnormal_dx[kdir][ldir][ldof];
		      dsvector_dx[1][kdir][ldir][ldof] = sign * dstangent_dx[0][kdir][ldir][ldof];
		    }
		}
	    }
	}
      break;
      
    case 3:
      for (kdir = 0; kdir < ielem_surf_dim+1; kdir++)
	{
	  svector[0][kdir] = snormal[kdir];
	  svector[1][kdir] = sign * stangent[0][kdir];
	  svector[2][kdir] = sign * stangent[1][kdir];
	  for (ldir = 0; ldir < ielem_surf_dim+1; ldir++)
	    {
	      for ( id=0; id<ei[pg->imtrx]->dof[ShapeVar]; id++ ) 
		{
		  Id = Proc_Elem_Connect[iconnect_ptr + id];
		  ldof = ei[pg->imtrx]->ln_to_dof[ShapeVar][id];
		  if (Dolphin[pg->imtrx][I][MESH_DISPLACEMENT1] > 0  && Dolphin[pg->imtrx][Id][MESH_DISPLACEMENT1] > 0 )
		    { 
		      
		      dsvector_dx[0][kdir][ldir][ldof] = dsnormal_dx[kdir][ldir][ldof];
		      dsvector_dx[1][kdir][ldir][ldof] = sign * dstangent_dx[0][kdir][ldir][ldof];
		      dsvector_dx[2][kdir][ldir][ldof] = sign * dstangent_dx[1][kdir][ldir][ldof];
		    }
		}
	    }
	}
      break;
      
    default:
      exit(-1);
      break;
    }

  /* Correct residual equation first at local node "irow_index" or global node "I" */
  /*                Rx -> Rn    and Ry -> Rt                                       */
  /*       i.e.,    Rn = nx*Rx + ny*Ry + nz*Rz                                     */
  /*                Rt1 = t1x*Rx + t1y*Ry + t1z*Rz                                 */
  /*                Rt2 = t2x*Rx + t2y*Ry + t2z*Rz                                 */
  
  /*Project residual into n-t space */
  for (kdir = 0; kdir < ielem_surf_dim+1; kdir++)
    {
      rotated_resid[kdir] = 0.;
      for (ldir = 0; ldir < ielem_surf_dim+1; ldir++)
	{
	  rotated_resid[kdir] +=
	    svector[kdir][ldir]*lec->R[upd->ep[pg->imtrx][R_MESH1 + ldir]][irow_index] ;
	}
    }
  
  
  /*   Now correct Jacobian                                               */
  
  if (af->Assemble_Jacobian) 
    {
      /* First, we have no global Jacobian, 
       * so zero it, but use it for a place holder to 
       * rotate the local jacobians 
       */
      
      siz = sizeof(double)* DIM*DIM*MDE;
      memset( rotated_jacobian_vector,0,siz);    

      /* Now correct for rotation */
      for (kdir = 0; kdir < ielem_surf_dim+1; kdir++)
	{
	  /* loop over sensitivities */
	  for ( j=0; j<ielem_surf_dim+1; j++)
	    {
	      
	      var=R_MESH1+j;
	      pvar = upd->vp[pg->imtrx][var];

	      for ( n=0; n< ei[pg->imtrx]->dof[var]; n++ ) 
		{
		  
		  for (ldir = 0; ldir < ielem_surf_dim+1; ldir++)
		    {
		      rotated_jacobian_vector[kdir][j][n]+=
			svector[kdir][ldir] * lec->J[upd->ep[pg->imtrx][R_MESH1 + ldir]][pvar][irow_index][n];
		    }
		}
	      
	       for ( id=0; id<ei[pg->imtrx]->dof[ShapeVar]; id++ ) 
		{
		  Id = Proc_Elem_Connect[iconnect_ptr + id];
		  ldof = ei[pg->imtrx]->ln_to_dof[ShapeVar][id];
		  if (Dolphin[pg->imtrx][I][MESH_DISPLACEMENT1] > 0  && Dolphin[pg->imtrx][Id][MESH_DISPLACEMENT1] > 0 )
		    { 
		      
		      for (ldir = 0; ldir < ielem_surf_dim+1; ldir++)
			{
			  rotated_jacobian_vector[kdir][j][ldof] +=  
			    dsvector_dx[kdir][ldir][j][ldof]
			       *lec->R[upd->ep[pg->imtrx][R_MESH1+ldir]][irow_index] ;
			}
		    } /* end of Baby_Dolphin[pg->imtrx] */
		}	  
	    } /* end of loop over sensitivities */
	}

      /* reinject back into lec-J for global assembly */
      for (kdir = 0; kdir < ielem_surf_dim+1; kdir++)
	{
	  for ( j=0; j<ielem_surf_dim+1; j++)
	    {
	      
	      var=R_MESH1+j;
	      pvar = upd->vp[pg->imtrx][var];
	      
	      for ( n=0; n< ei[pg->imtrx]->dof[var]; n++ ) 
		{
		  
		  lec->J[upd->ep[pg->imtrx][R_MESH1 + kdir]][pvar][irow_index][n]
		    = rotated_jacobian_vector[kdir][j][n];
		}
	    }
	} 


      /* mesh wrt. pressure */
      var = PRESSURE;
      if (pd->v[pg->imtrx][var]){
	pvar = upd->vp[pg->imtrx][var];
	for ( n=0; n<ei[pg->imtrx]->dof[var]; n++) {
	  
	  for (kdir = 0; kdir < ielem_surf_dim+1; kdir++)
	    {
	      rotated_jacobian_scalar[kdir][n] = 0.;
	    }
	  
	  for (ldir = 0; ldir < ielem_surf_dim+1; ldir++)
	    {
	      peq = upd->ep[pg->imtrx][R_MESH1+ldir];
	      for (kdir = 0; kdir < ielem_surf_dim+1; kdir++)
		{
		  rotated_jacobian_scalar[kdir][n] += 
		    svector[kdir][ldir]*lec->J[peq][pvar][irow_index][n];
		}
	    }
	  
	  /*reinject back into lec-J for global assembly */
	  for (ldir = 0; ldir < ielem_surf_dim+1; ldir++)
	    {
	      peq = upd->ep[pg->imtrx][R_MESH1+ldir];
	      lec->J[peq][pvar][irow_index][n] = rotated_jacobian_scalar[ldir][n];
	    }
	} /* end of loop over nodes */
      } /* end of if variable */
      
      /* mesh wrt. temperature */
      var = TEMPERATURE;
      if (pd->v[pg->imtrx][var]){
	pvar = upd->vp[pg->imtrx][var];
	for ( n=0; n<ei[pg->imtrx]->dof[var]; n++) {
	  for (kdir = 0; kdir < ielem_surf_dim+1; kdir++)
	    {
	      rotated_jacobian_scalar[kdir][n] = 0.;
	    }
	  
	  for (ldir = 0; ldir < ielem_surf_dim+1; ldir++)
	    {
	      peq = upd->ep[pg->imtrx][R_MESH1+ldir];
	      for (kdir = 0; kdir < ielem_surf_dim+1; kdir++)
		{
		  rotated_jacobian_scalar[kdir][n] += 
		    svector[kdir][ldir]*lec->J[peq][pvar][irow_index][n];
		}
	    }

	  /*reinject back into lec-J for global assembly */
	  for (ldir = 0; ldir < ielem_surf_dim+1; ldir++)
	    {
	      peq = upd->ep[pg->imtrx][R_MESH1+ldir];
	      lec->J[peq][pvar][irow_index][n] = rotated_jacobian_scalar[ldir][n];
	    }
	} /* end of loop over nodes */
      } /* end of if variable */
      
      
      /* mesh wrt. velocity */
      for (jvar = 0; jvar < ielem_surf_dim+1; jvar++)
	{
	  var = VELOCITY1 + jvar;
	  if (pd->v[pg->imtrx][var]){
	    pvar = upd->vp[pg->imtrx][var];
	    for ( n=0; n<ei[pg->imtrx]->dof[var]; n++) {
	      
	      for (kdir = 0; kdir < ielem_surf_dim+1; kdir++)
		{
		  rotated_jacobian_vector[kdir][jvar][n] = 0.;
		}
	      
	      for (ldir = 0; ldir < ielem_surf_dim+1; ldir++)
		{
		  peq = upd->ep[pg->imtrx][R_MESH1+ldir];
		  for (kdir = 0; kdir < ielem_surf_dim+1; kdir++)
		    {
		      rotated_jacobian_vector[kdir][jvar][n] += 
			svector[kdir][ldir]*lec->J[peq][pvar][irow_index][n];
		    }
		}
	      /*reinject back into lec-J for global assembly */
	      for (ldir = 0; ldir < ielem_surf_dim+1; ldir++)
		{
		  peq = upd->ep[pg->imtrx][R_MESH1+ldir];
		  lec->J[peq][pvar][irow_index][n] = rotated_jacobian_vector[ldir][jvar][n];
		}
	    } /* end of loop over nodes */
	  } /* end of if variable */

	} /* end of loop over jvar direction */



      
      /* mesh wrt. species concentration */
      var = MASS_FRACTION;
      if (pd->v[pg->imtrx][var]){
	for (w=0; w<pd->Num_Species_Eqn; w++)
	  {
	    for ( n=0; n<ei[pg->imtrx]->dof[var]; n++) {
	      
	      for (kdir = 0; kdir < ielem_surf_dim+1; kdir++)
		{
		  rotated_jacobian_conc[kdir][w][n] = 0.;
		}
	      
	      for (ldir = 0; ldir < ielem_surf_dim+1; ldir++)
		{
		  peq = upd->ep[pg->imtrx][R_MESH1+ldir];
		  for (kdir = 0; kdir < ielem_surf_dim+1; kdir++)
		    {
		      
		      rotated_jacobian_conc[kdir][w][n] += 
			svector[kdir][ldir]*lec->J[peq][MAX_PROB_VAR + w][irow_index][n];
		    }
		}
	      
	      /*reinject back into lec-J for global assembly */
	      for (ldir = 0; ldir < ielem_surf_dim+1; ldir++)
		{
		  peq = upd->ep[pg->imtrx][R_MESH1+ldir];
		  lec->J[peq][MAX_PROB_VAR + w][irow_index][n] = rotated_jacobian_conc[ldir][w][n];
		}
	    } /* end of loop over nodes */
	  }   /* end of loop over concentration */
      } /* end of if variable */
      
    } /* end of if Newton */

  /* Put rotated residual back into lec for scattering into global matrix */
  for (kdir = 0; kdir < ielem_surf_dim+1; kdir++)
    {
      lec->R[ upd->ep[pg->imtrx][R_MESH1+kdir] ][irow_index] = rotated_resid[kdir];
    } 
  
} /* END of rotate_res_jac_mesh */
/********************************************************************************/
/********************************************************************************/
/********************************************************************************/
/*ARGSUSED*/
void 
rotate_mesh_eqn (
    int id,                             /* Elemental stiffness matrix row index */
    int I,                              /* Global node number                   */
    int iconnect_ptr,  		        /* Pointer to beginning of connectivity
					   list for the current element	        */
    const int dim,		        /* physical dim of problem              */
    struct Aztec_Linear_Solver_System *ams)

    /****************************************************************************
     *
     * rotate_mesh_eqn():
     *
     * Function which transforms the local element mesh equations and corresponding
     * local element jacobian from unrotated form to rotated form for a 
     * single local node number i, corresponding to processor node number I.
     *
     * The following transformations are performed:
     *
     *   Rn =  nx*Rx +  ny*Ry +  nz*Rz                                     
     *  Rt1 = t1x*Rx + t1y*Ry + t1z*Rz                
     *  Rt2 = t2x*Rx + t2y*Ry + t2z*Rz
     *
     *  Rn will occupy the MESH1 position, Rt1 will occupy the MESH2 position
     *  and Rt2 will occupy the MESH3 position. If the 'ok' field in the
     *  Rotation_Vectors structure is not set, then the corresponding entry
     *  in the mesh local element residual vector and corresponding local
     *  element jacobian terms will not be touched.
     *
     *  Note, direct injection into the global matrix is carried out for the
     *  case of a dependence of a rotation vector on the position of a node
     *  where the node is not part of the current element. This is unavoidable,
     *  because the local stiffness matrix does not have a column entry for a
     *  node which is not part of the current element.
     *
     *  NOTE: The algorithm in this routine needs to be enhanced to rotate
     *        the jacobian terms for all variable types. Currently, it just
     *        does a couple of variable types.
     *
     *       Author:          P. R. Schunk (1511)
     *       Revisor:         R. A. Cairncross 5 September 1996
     *       Revisor:         R. R. Rao  20 July 1999 (lec revolution)
     ****************************************************************************/
{
  int index, index_eqn, index_var, ndof;
  const int eq = VECT_EQ_MESH;
  int var, eqn, pvar, peqn, j, n, b, J, v, pv, kv;
  int kdir, ldir, w;
  double *a = ams->val;
  int *ija  = ams->bindx; 
  int *rpntr, *bpntr, *bindx, *indx;
  int row_dofs, blk_row, K, blk_col;
  int peqn_mesh[MAX_PDIM]; /* Stores the UPD Index for the mesh variables */
  double rotated_resid[MAX_PDIM];
  double rotated_jacobian_scalar[MAX_PDIM][MDE];
  double  rotated_jacobian_vector[MAX_PDIM][MAX_PDIM][MDE];
  NODAL_VARS_STRUCT *nv = Nodes[I]->Nodal_Vars_Info[pg->imtrx];
  ROTATION_VECTORS_STRUCT *rot, *rot_n = rotation[I][eq][0];
  ROTATION_VECTORS_STRUCT *rot_t1 = rotation[I][eq][1];
  ROTATION_VECTORS_STRUCT *rot_t2 = 0;
  int doSpecialInjection = FALSE;

  if (dim > 2) {
    rot_t2 = rotation[I][eq][2];
  }

  /*
   * calculate some temporary variables and check for 
   * cases for which this subroutine isn't sufficiently
   * general.
   */
  for (ldir = 0; ldir < dim; ldir++) {
    ndof = get_nv_ndofs(nv, R_MESH1 + ldir);
    if (ndof != 1) {
      sprintf(Err_Msg,
	      "rotate_mesh_eqn: Limited to 1 mesh dof: found %d dofs\n",
	      ndof);
      EH(-1, Err_Msg);
    }
    peqn_mesh[ldir] = upd->ep[pg->imtrx][R_MESH1 + ldir];
  }

  /* Now add on projection into n-t space */

  for (kdir = 0; kdir < dim; kdir++) {
    rot = rotation[I][eq][kdir];
    if (rot->ok) {
      rotated_resid[kdir] = 0.0;     
      for (ldir = 0; ldir < dim; ldir++) {
	rotated_resid[kdir] += (rot->vector[ldir] *
				lec->R[peqn_mesh[ldir]][id]);
      }
    }
  }
 
  /* 
   *   Now correct the Jacobian 
   */
  if (af->Assemble_Jacobian) {
      
    /*
     * Handle direct insertion of rotation vector sensitivities into 
     * global sparse matrices first
     */
      
    /*
     * Loop over the rows that are to be rotated
     */
    for (kdir = 0; kdir < dim; kdir++) {
      eqn = R_MESH1 + kdir;
      rot = rotation[I][eq][kdir];
      if (rot->ok) {
	/*
	 * Loop over the mesh columns in the Jacobian
	 */
	for (b = 0; b < dim; b++) {	      
	  var = R_MESH1 + b;
	  pvar = upd->vp[pg->imtrx][var];

	  if (strcmp(Matrix_Format, "msr") == 0) {
	    /*
	     * Find the global equation number
	     */
	    if ((index_eqn = Index_Solution(I, eqn, 0, 0, -2, pg->imtrx)) == -1) { 
	      EH(-1, "Cant find eqn index");
	    }
	    /*
	     * Loop over the nodes that determine the value of the
	     * current rotation vector, J is the global node number
	     */
	    for (j = 0; j < rot->d_vector_n; j++ )  {
	      J = rot->d_vector_J[j];
	      if (Dolphin[pg->imtrx][I][MESH_DISPLACEMENT1] > 0  &&
		  Dolphin[pg->imtrx][J][MESH_DISPLACEMENT1] > 0 ) { 
		/* 
		 * find entry in global matrix - note that the sensitivities of the 
		 * rotation vector may be in a different element than the current 
		 * element 
		 */
		if ((index_var = Index_Solution(J, var, 0, 0, -2, pg->imtrx)) != -1) { 
		  index = (index_eqn == index_var) ? index_eqn :
		      in_list(index_var, ija[index_eqn], ija[index_eqn+1], ija);
			  
		  for (ldir = 0; ldir < dim; ldir++) {
		    a[index] += 
			rot->d_vector_dx[ldir][b][j] * lec->R[peqn_mesh[ldir]][id];
		  }
		}
	      }
	    }
	  } else if (strcmp(Matrix_Format, "vbr") == 0) {
	    rpntr = ams->rpntr;
	    bpntr = ams->bpntr;
	    bindx = ams->bindx;
	    indx  = ams->indx;
	    a     = ams->val;
	    nv  = Nodes[I]->Nodal_Vars_Info[pg->imtrx];
	    blk_row = get_nv_offset_idof(nv, eqn, 0, 0, NULL);
	    row_dofs = rpntr[I+1] - rpntr[I];
	    for (j = 0; j < rot->d_vector_n; j++) {
	      J = rot->d_vector_J[j];
	      if (Dolphin[pg->imtrx][I][MESH_DISPLACEMENT1] > 0 &&
		  Dolphin[pg->imtrx][J][MESH_DISPLACEMENT1] > 0) { 
		K = in_list(J, bpntr[I], bpntr[I+1], bindx);
#ifdef DEBUG_HKM
		if (K == -1) {
		  fprintf(stderr,"Couldn't find proc blk_col %d in blk_row %d\n",
			  J, I);
		  EH(-1, "rotate_mesh_eqn ERROR: vbr matrix structure error");
		}
#endif
		nv = Nodes[J]->Nodal_Vars_Info[pg->imtrx];
		blk_col = get_nv_offset_idof(nv, var, 0, 0, NULL);
		/*
		 * Calculate the index into the matrix:
		 *     here we assume only 1 mesh dof per node
		 */
		index = indx[K] + row_dofs * blk_col + blk_row; 
		for (ldir = 0; ldir < dim; ldir++) {
		  a[index] += (rot->d_vector_dx[ldir][b][j] *
			       lec->R[peqn_mesh[ldir]][id]);
		}
	      }
	    }  
	  } else if (strcmp(Matrix_Format, "epetra") == 0) {
            /*
             * Find the global equation number
             */
            if ((index_eqn = Index_Solution(I, eqn, 0, 0, -2, pg->imtrx)) == -1) {
              EH(-1, "Cant find eqn index");
            }
            /*
             * Loop over the nodes that determine the value of the
             * current rotation vector, J is the global node number
             */
            for (j = 0; j < rot->d_vector_n; j++ )  {
              double sum_val;
              int global_row;
              int global_col;

              J = rot->d_vector_J[j];
              if (Dolphin[pg->imtrx][I][MESH_DISPLACEMENT1] > 0  &&
                  Dolphin[pg->imtrx][J][MESH_DISPLACEMENT1] > 0 ) {
                /*
                 * find entry in global matrix - note that the sensitivities of the
                 * rotation vector may be in a different element than the current
                 * element
                 */
                if ((index_var = Index_Solution(J, var, 0, 0, -2, pg->imtrx)) != -1) {
                  sum_val = 0;
                  for (ldir = 0; ldir < dim; ldir++) {
                    sum_val +=
                        rot->d_vector_dx[ldir][b][j] * lec->R[peqn_mesh[ldir]][id];
                  }
                  global_row = ams->GlobalIDs[index_eqn];
                  global_col = ams->GlobalIDs[index_var];
                  EpetraSumIntoGlobalRowMatrix(ams->RowMatrix, global_row, 1, &sum_val, &global_col);
                }
              }
            }
          } else {
	    /* FRONTAL SOLVER SECTION:
	     *
	     * This loop adds only a portion of the rotation vector sensitivities 
	     * but does not have any direct injection 
	     */
	    doSpecialInjection = TRUE;
	    for (w = 0; w < ei[pg->imtrx]->dof[var]; w++) {
	      rotated_jacobian_vector[kdir][b][w] = 0.0;
	      J = ei[pg->imtrx]->gnn_list[var][w];
	      j = in_list(J, 0, rot->d_vector_n, rot->d_vector_J);
	      if (j != -1) {  
		for (ldir = 0; ldir < dim; ldir++) {
		  rotated_jacobian_vector[kdir][b][w] += 
		      (rot->d_vector_dx[ldir][b][j] *
		       lec->R[peqn_mesh[ldir]][id]);
		}
	      }
	    }
	  }
	} /* end of loop over sensitivities */
      }
    }

    /*
     * Let's rotate columns of the local stiffness matrix
     * one column at a time
     */
    for (v = V_FIRST; v < V_LAST; v++) {
      pv = upd->vp[pg->imtrx][v];
      if (pv != -1) {
	if (v == MASS_FRACTION) {
	  for (kv = 0; kv < upd->Max_Num_Species_Eqn; kv++) {	  
	    pv = MAX_PROB_VAR + kv;
	    switch (dim) {
	    case 2:
		for (n = 0; n < ei[pg->imtrx]->dof[v]; n++) {
		  rotated_jacobian_scalar[0][n] = 
		      rot_n->vector[0] * lec->J[peqn_mesh[0]][pv][id][n] +
		      rot_n->vector[1] * lec->J[peqn_mesh[1]][pv][id][n];
		  rotated_jacobian_scalar[1][n] = 
		      rot_t1->vector[0] * lec->J[peqn_mesh[0]][pv][id][n] +
		      rot_t1->vector[1] * lec->J[peqn_mesh[1]][pv][id][n];
		}	    
		break;
	    case 3:
		for (n = 0; n < ei[pg->imtrx]->dof[v]; n++) {
		  rotated_jacobian_scalar[0][n] = 
		      rot_n->vector[0] * lec->J[peqn_mesh[0]][pv][id][n] +
		      rot_n->vector[1] * lec->J[peqn_mesh[1]][pv][id][n] +
		      rot_n->vector[2] * lec->J[peqn_mesh[2]][pv][id][n];
		  rotated_jacobian_scalar[1][n] =
		      rot_t1->vector[0] * lec->J[peqn_mesh[0]][pv][id][n] +
		      rot_t1->vector[1] * lec->J[peqn_mesh[1]][pv][id][n] +
		      rot_t1->vector[2] * lec->J[peqn_mesh[2]][pv][id][n];
		  rotated_jacobian_scalar[2][n] =
		      rot_t2->vector[0] * lec->J[peqn_mesh[0]][pv][id][n] +
		      rot_t2->vector[1] * lec->J[peqn_mesh[1]][pv][id][n] +
		      rot_t2->vector[2] * lec->J[peqn_mesh[2]][pv][id][n];
		}
		break;
	    }
	    /*
	     * Now reinject the column into the matrix
	     */
	    for (kdir = 0; kdir < dim; kdir++) {
	      if (rotation[I][eq][kdir]->ok) {
		peqn = peqn_mesh[kdir];
		for (n = 0; n < ei[pg->imtrx]->dof[v]; n++) {	      
		  lec->J[peqn][pv][id][n] = rotated_jacobian_scalar[kdir][n];
		}
	      }
	    }
	  } /* end of loop over Max_Num_Species */
	} else if ( v >= MESH_DISPLACEMENT1 && v <= MESH_DISPLACEMENT3 ) {
	  switch (dim) {
	  case 2:
	      for (n = 0; n < ei[pg->imtrx]->dof[v]; n++) {
		rotated_jacobian_scalar[0][n] = 
		    rot_n->vector[0] * lec->J[peqn_mesh[0]][pv][id][n] +
		    rot_n->vector[1] * lec->J[peqn_mesh[1]][pv][id][n];
		rotated_jacobian_scalar[1][n] = 
		    rot_t1->vector[0] * lec->J[peqn_mesh[0]][pv][id][n] +
		    rot_t1->vector[1] * lec->J[peqn_mesh[1]][pv][id][n];
	      }	    
	      break;
	  case 3:
	      for (n = 0; n < ei[pg->imtrx]->dof[v]; n++) {
		rotated_jacobian_scalar[0][n] = 
		    rot_n->vector[0] * lec->J[peqn_mesh[0]][pv][id][n] +
		    rot_n->vector[1] * lec->J[peqn_mesh[1]][pv][id][n] +
		    rot_n->vector[2] * lec->J[peqn_mesh[2]][pv][id][n];
		rotated_jacobian_scalar[1][n] =
		    rot_t1->vector[0] * lec->J[peqn_mesh[0]][pv][id][n] +
		    rot_t1->vector[1] * lec->J[peqn_mesh[1]][pv][id][n] +
		    rot_t1->vector[2] * lec->J[peqn_mesh[2]][pv][id][n];
		rotated_jacobian_scalar[2][n] =
		    rot_t2->vector[0] * lec->J[peqn_mesh[0]][pv][id][n] +
		    rot_t2->vector[1] * lec->J[peqn_mesh[1]][pv][id][n] +
		    rot_t2->vector[2] * lec->J[peqn_mesh[2]][pv][id][n];
	      }
	      break;
	  }
	  /*
	   * Now reinject the column into the matrix
	   */
	  for (kdir = 0; kdir < dim; kdir++) {
	    if (rotation[I][eq][kdir]->ok) {
	      peqn = peqn_mesh[kdir];
	      for (n = 0; n < ei[pg->imtrx]->dof[v]; n++) {	      
		lec->J[peqn][pv][id][n] = rotated_jacobian_scalar[kdir][n];
	      }
	    }
	  }
	}
      }
    }

    /* 
     * reinject back into lec-J the special part of the
     * matrix that the frontal solver can handle for global assembly
     */
    if (doSpecialInjection) {
      for (kdir = 0; kdir < dim; kdir++) {
	peqn = peqn_mesh[kdir];
	rot = rotation[I][eq][kdir];
	if (rot->ok) {
	  for (j = 0; j < dim; j++) {
	    var = R_MESH1 + j;
	    pvar = upd->vp[pg->imtrx][var];
	    for (n = 0; n < ei[pg->imtrx]->dof[var]; n++ ) {
	      lec->J[peqn][pvar][id][n] += rotated_jacobian_vector[kdir][j][n];
	    }
	  }
	}
      }
    }
  }

  /*
   *  Put rotated residual back into lec for scattering into global matrix
   */
  for (kdir = 0; kdir < dim; kdir++) {
    rot = rotation[I][eq][kdir];
    if (rot->ok) {
      lec->R[peqn_mesh[kdir]][id] = rotated_resid[kdir];
    }
  }    
} /* END of rotate_mesh_eqn */
/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
/*ARGSUSED*/
void 
rotate_momentum_eqn (
    int id,                             /* Elemental stiffness matrix row index */
    int I,                              /* Global node number                   */
    int iconnect_ptr,  		        /* Pointer to beginning of connectivity list
				           for the current element	         */
    int dim,                            /* physical dim of problem              */
    struct Aztec_Linear_Solver_System *ams )  

/* 
 * Function which corrects the global residual vector "resid_vect" and
 * the global Jacobian vector "a" so that the vector mesh equations for surface nodes 
 * are projected into a normal and tangential coordinate system.
 *
 *       Author:          P. R. Schunk (1511)
 *       Revisor:         R. A. Cairncross 5 September 1996
 *       Revisor:         R. R. Rao  20 July 1999 (lec revolution)
 */

{
/* LOCAL VARIABLES */
  int       eq, var, eqn, pvar, peq;
  int 	    j, n, b, index, J;
  int       kdir, ldir, w, jvar;
  unsigned int       siz;
  int      *ija = ams->bindx;

  double   *a = ams->val;
  double    rotated_resid[MDE];
  double    rotated_jacobian_vector[MAX_PDIM][MAX_PDIM][MDE];
  double    rotated_jacobian_scalar[MAX_PDIM][MDE];
  double    rotated_jacobian_conc[MAX_PDIM][MAX_CONC][MDE];

/************************ EXECUTION BEGINS **********************************/
  eq = VECT_EQ_MOM;

/* Correct residual equation first at local node "id" or global node "I" */
/*                Rx -> Rn    and Ry -> Rt                                   */
/*       i.e.,    Rn = nx*Rx + ny*Ry + nz*Rz                                 */
/*                Rt1 = t1x*Rx + t1y*Ry + t1z*Rz                             */
/*                Rt2 = t2x*Rx + t2y*Ry + t2z*Rz                             */


  /* Now add on projection into n-t space */
  for (kdir = 0; kdir < dim; kdir++)
    {
      if (rotation[I][eq][kdir]->ok) {
	rotated_resid[kdir] = 0.;     
	for (ldir = 0; ldir < dim; ldir++)
	  {
	    peq = upd->ep[pg->imtrx][R_MOMENTUM1+ldir];
	    rotated_resid[kdir] += 
	      rotation[I][eq][kdir]->vector[ldir]*lec->R[peq][id];
	  }
      }
    } /* end of loop over direction */


  /*                                                                      */
  /*   Now correct Jacobian                                               */
   if (af->Assemble_Jacobian) 
    {
      
      siz = sizeof(double)* DIM*DIM*MDE;
      memset( rotated_jacobian_vector,0,siz);    
      
      
      /* Now correct for rotation */
      
      for (kdir = 0; kdir < dim; kdir++) {
	if (rotation[I][eq][kdir]->ok) {
	  eqn=R_MOMENTUM1+kdir;
	  /* loop over sensitivities */
	  for ( b=0; b<dim; b++)
	    {
	      
	      var=MESH_DISPLACEMENT1+b;
	      pvar = upd->vp[pg->imtrx][var];
	      
	      for ( j=0; j< ei[pg->imtrx]->dof[var]; j++ ) 
		{
		  
		  for (ldir = 0; ldir < dim; ldir++)
		    {
		      rotated_jacobian_vector[kdir][b][j] += 
			  rotation[I][eq][kdir]->vector[ldir] * 
			  lec->J[upd->ep[pg->imtrx][R_MOMENTUM1 + ldir]][pvar][id][j];
		    }
		}

	      /* direct insertion of rotation vector sensitivities into global sparse matrices */

	      if (strcmp(Matrix_Format, "msr") == 0 ) {
		for (j = 0; j < rotation[I][eq][kdir]->d_vector_n; j++ )  {
		  int ktype, ndof, index_eqn, index_var;
		  J = rotation[I][eq][kdir]->d_vector_J[j];
		  if (Dolphin[pg->imtrx][I][R_MOMENTUM1] > 0  && Dolphin[pg->imtrx][J][MESH_DISPLACEMENT1] > 0) { 
		    /* find entry in global matrix - note that the sensitivities of the 
		     * rotation vector may be in a different element than the current 
		     * element */
		    ktype = 0;
		    ndof = 0;
		    if ((index_eqn =  Index_Solution(I, eqn, ktype, ndof, -2, pg->imtrx)) == -1 ) { 
		      EH(-1, "Cant find eqn index");
		    }
		    if ((index_var =  Index_Solution(J, var, ktype, ndof, -2, pg->imtrx)) == -1 ) { 
		      EH(-1, "Cant find var index");
		    }
		      
		    index = (index_eqn == index_var) ? index_eqn :
			in_list(index_var, ija[index_eqn], ija[index_eqn+1], ija);
			  
		    /* !!!!!!!!!!!!!!!!! Warning !!!!!!!!!!!!!!!!! */
		    /* I need to fix this direct injection into a, but it isn't clear how. */
		    for (ldir = 0; ldir < dim; ldir++) {
		      a[index] +=  
			  rotation[I][eq][kdir]->d_vector_dx[ldir][b][j]
			  * lec->R[upd->ep[pg->imtrx][R_MESH1+ldir]][id];
		    }
		  } /* end of Baby_Dolphin[pg->imtrx] */
		}
	      }
	      else if ( strcmp(Matrix_Format, "vbr") == 0 )
		{
		  int *rpntr = ams->rpntr;
		  int *bpntr = ams->bpntr;
		  int *bindx = ams->bindx;
		  int *indx  = ams->indx;
		  double *a  = ams->val;

		  int row_dofs, blk_row;

		  int K, blk_col;
		  NODAL_VARS_STRUCT *nv;
		  nv = Nodes[I]->Nodal_Vars_Info[pg->imtrx];
		  blk_row = get_nv_offset_idof(nv, eqn, 0, 0, NULL);
		  row_dofs = rpntr[I+1] - rpntr[I];

		  for ( j=0; j<rotation[I][eq][kdir]->d_vector_n; j++ ) 
		    {
		      J = rotation[I][eq][kdir]->d_vector_J[j];

		      if (Dolphin[pg->imtrx][I][R_MOMENTUM1] > 0  && Dolphin[pg->imtrx][J][MESH_DISPLACEMENT1] > 0 )
			{ 

			  K = in_list( J, bpntr[I], bpntr[I+1], bindx );   /* K is the block index */
			  nv = Nodes[J]->Nodal_Vars_Info[pg->imtrx];
			  blk_col = get_nv_offset_idof(nv, var, 0, 0, NULL);

			  index = indx[K] + row_dofs* blk_col + blk_row; /* here we assume only 1 mesh dof per node */

			  for (ldir = 0; ldir < dim; ldir++)
			    {
			      a[index] +=  
				rotation[I][eq][kdir]->d_vector_dx[ldir][b][j]
				* lec->R[upd->ep[pg->imtrx][R_MESH1+ldir]][id];
			    }
			}
		    }
			  
		}
	      else if (strcmp(Matrix_Format, "epetra") == 0) {
	        // Direct translation from MSR
                for (j = 0; j < rotation[I][eq][kdir]->d_vector_n; j++) {
                  double sum_val;
                  int global_row;
                  int global_col;
                  int ktype, ndof, index_eqn, index_var;
                  J = rotation[I][eq][kdir]->d_vector_J[j];
                  if (Dolphin[pg->imtrx][I][R_MOMENTUM1] > 0
                      && Dolphin[pg->imtrx][J][MESH_DISPLACEMENT1] > 0) {
                    /* find entry in global matrix - note that the sensitivities of the
                     * rotation vector may be in a different element than the current
                     * element */
                    ktype = 0;
                    ndof = 0;
                    if ((index_eqn = Index_Solution(I, eqn, ktype, ndof, -2, pg->imtrx))
                        == -1) {
                      EH(-1, "Cant find eqn index");
                    }
                    if ((index_var = Index_Solution(J, var, ktype, ndof, -2, pg->imtrx))
                        == -1) {
                      EH(-1, "Cant find var index");
                    }

                    /* !!!!!!!!!!!!!!!!! Warning !!!!!!!!!!!!!!!!! */
                    /* I need to fix this direct injection into a, but it isn't clear how. */
                    sum_val = 0;

                    for (ldir = 0; ldir < dim; ldir++) {
                      sum_val += rotation[I][eq][kdir]->d_vector_dx[ldir][b][j]
                          * lec->R[upd->ep[pg->imtrx][R_MESH1 + ldir]][id];
                    }
                    global_row = ams->GlobalIDs[index_eqn];
                    global_col = ams->GlobalIDs[index_var];
                    EpetraSumIntoGlobalRowMatrix(ams->RowMatrix, global_row, 1, &sum_val, &global_col);
                  } /* end of Baby_dolphin */
                }
	      }
	      else
		{
		  /* 
		   * This loop addes only a portion of the rotation vector sensitivities 
		   * but does not have any direct injection 
		   */
		  
		  for( w=0; w< ei[pg->imtrx]->dof[var];w++)
		    {
		      J = ei[pg->imtrx]->gnn_list[var][w];
		      
		      if( ( j = in_list(  J,0, rotation[I][eq][kdir]->d_vector_n, rotation[I][eq][kdir]->d_vector_J ) ) != -1 )
			{
			  
			  for (ldir = 0; ldir < dim; ldir++)
			    {
			      rotated_jacobian_vector[kdir][b][w] += 
				rotation[I][eq][kdir]->d_vector_dx[ldir][b][j] * lec->R[upd->ep[pg->imtrx][R_MESH1+ldir]][id];
			    }
			}
		    }
		}
	    } /* end of loop over sensitivities */
	}
      }

      /* reinject back into lec-J for global assembly */
      for (kdir = 0; kdir < dim; kdir++)
	{
	  if (rotation[I][eq][kdir]->ok) {
	    /* loop over sensitivities */
	    for ( j=0; j<dim; j++)
	      {
		var= MESH_DISPLACEMENT1+j;
		pvar = upd->vp[pg->imtrx][var];
		for ( n=0; n< ei[pg->imtrx]->dof[var]; n++ ) 
		  {
		    lec->J[upd->ep[pg->imtrx][R_MOMENTUM1 + kdir]][pvar][id][n]
		      = rotated_jacobian_vector[kdir][j][n];
		  }
	      }
	  } /* end of loop over sensitivities */
	}

      
      /* momentum wrt. pressure */
      var = PRESSURE;
      if (pd->v[pg->imtrx][var]){
	pvar = upd->vp[pg->imtrx][var];
	for ( n=0; n<ei[pg->imtrx]->dof[var]; n++) {

	  rotated_jacobian_scalar[0][n] = 0.;
	  rotated_jacobian_scalar[1][n] = 0.;
	  rotated_jacobian_scalar[2][n] = 0.;
	  
	  for (kdir = 0; kdir < dim; kdir++) {
	    if (rotation[I][eq][kdir]->ok) {
	      eqn=R_MOMENTUM1+kdir;
	      for (ldir = 0; ldir < dim; ldir++)
		{
		  rotated_jacobian_scalar[kdir][n] += 
		    rotation[I][eq][kdir]->vector[ldir] * lec->J[upd->ep[pg->imtrx][R_MOMENTUM1+ldir]][pvar][id][n];
		}
	      }
	    }

	} /* end of loop over nodes */
	
	/* reinject d_mesh/d_pressure back into lec-J for global assembly */
	for (kdir = 0; kdir < dim; kdir++)
	  {
	    if (rotation[I][eq][kdir]->ok) {
	      /* loop over sensitivities */
	      for ( n=0; n< ei[pg->imtrx]->dof[var]; n++ ) 
		{
		  
		  lec->J[upd->ep[pg->imtrx][R_MOMENTUM1+ kdir]][pvar][id][n]
		    = rotated_jacobian_scalar[kdir][n];
		}
	    }
	  }

      } /* end of if variable */
      
      /* momentum wrt. temperature */
      var = TEMPERATURE;
      if (pd->v[pg->imtrx][var]){
	pvar = upd->vp[pg->imtrx][var];
	for ( n=0; n<ei[pg->imtrx]->dof[var]; n++) {
	  	  
	  rotated_jacobian_scalar[0][n] = 0.;
	  rotated_jacobian_scalar[1][n] = 0.;
	  rotated_jacobian_scalar[2][n] = 0.;
	    
	  for (kdir = 0; kdir < dim; kdir++)
	    {
	      if (rotation[I][eq][kdir]->ok) {
		eqn=R_MOMENTUM1+kdir;
		for (ldir = 0; ldir < dim; ldir++)
		  {
		    rotated_jacobian_scalar[kdir][n] += 
		      rotation[I][eq][kdir]->vector[ldir]*lec->J[upd->ep[pg->imtrx][R_MOMENTUM1+ldir]][pvar][id][n];
		    /*                                should this be a ldir or a kdir? | */
		  }
		}
	    }

	} /* end of loop over nodes */

	/* reinject d_mesh/d_temperature back into lec-J for global assembly */
	for (kdir = 0; kdir < dim; kdir++)
	  {
	    if (rotation[I][eq][kdir]->ok) {
	      /* loop over sensitivities */
	      for ( n=0; n< ei[pg->imtrx]->dof[var]; n++ ) 
		{
		  
		  lec->J[upd->ep[pg->imtrx][R_MOMENTUM1+ kdir]][pvar][id][n]
		    = rotated_jacobian_scalar[kdir][n];
		}
	    }
	  }
	


      } /* end of if variable */
#ifdef COUPLED_FILL
      /* momentum wrt. temperature */
      var = LS;
      if (pd->v[pg->imtrx][var]){
	pvar = upd->vp[pg->imtrx][var];
	for ( n=0; n<ei[pg->imtrx]->dof[var]; n++) {
	  
	  rotated_jacobian_scalar[0][n] = 0.;
	  rotated_jacobian_scalar[1][n] = 0.;
	  rotated_jacobian_scalar[2][n] = 0.;
	    
	  
	  for (kdir = 0; kdir < dim; kdir++)
	    {
	      if (rotation[I][eq][kdir]->ok) {
		eqn=R_MOMENTUM1+kdir;
		for (ldir = 0; ldir < dim; ldir++)
		  {
		    rotated_jacobian_scalar[kdir][n] += 
		      rotation[I][eq][kdir]->vector[ldir]*lec->J[upd->ep[pg->imtrx][R_MOMENTUM1+ldir]][pvar][id][n];
		    /*                                should this be a ldir or a kdir? | */
		  }
		}
	    }

	} /* end of loop over nodes */

	/* reinject d_mesh/d_temperature back into lec-J for global assembly */
	for (kdir = 0; kdir < dim; kdir++)
	  {
	    if (rotation[I][eq][kdir]->ok) {
	      /* loop over sensitivities */
	      for ( n=0; n< ei[pg->imtrx]->dof[var]; n++ ) 
		{
		  
		  lec->J[upd->ep[pg->imtrx][R_MOMENTUM1+ kdir]][pvar][id][n]
		    = rotated_jacobian_scalar[kdir][n];
		}
	    }
	  }
	


      } /* end of if variable */
#endif      
      
      /* momentum wrt. velocity */
      for (jvar = 0; jvar < dim; jvar++)
	{
	  var = VELOCITY1 + jvar;
	  if (pd->v[pg->imtrx][var]){
	    pvar = upd->vp[pg->imtrx][var];
	    for ( n=0; n<ei[pg->imtrx]->dof[var]; n++) {
	      
	      for (kdir = 0; kdir < dim; kdir++)
		{
		  rotated_jacobian_vector[kdir][jvar][n] = 0.;
		}
	      
	      for (kdir = 0; kdir < dim; kdir++)
		{
		  if (rotation[I][eq][kdir]->ok) {
		    eqn=R_MOMENTUM1+kdir;
		    for (ldir = 0; ldir < dim; ldir++)
		      {
			rotated_jacobian_vector[kdir][jvar][n] += 
			  rotation[I][eq][kdir]->vector[ldir]*lec->J[upd->ep[pg->imtrx][R_MOMENTUM1+ldir]][pvar][id][n];
		      }
		    }
		}

	    } /* end of loop over nodes */
	  } /* end of if variable */
	} /* end of loop over jvar direction */

      /* reinject back into lec-J for global assembly */
      for (kdir = 0; kdir < dim; kdir++)
	{
	  if (rotation[I][eq][kdir]->ok) {
	    /* loop over sensitivities */
	    for ( j=0; j<dim; j++)
	      {
		var= R_MOMENTUM1+j;
		pvar = upd->vp[pg->imtrx][var];
		for ( n=0; n< ei[pg->imtrx]->dof[var]; n++ ) 
		  {
		    lec->J[upd->ep[pg->imtrx][R_MOMENTUM1 + kdir]][pvar][id][n]
		      = rotated_jacobian_vector[kdir][j][n];
		  }
	      }
	  } /* end of loop over sensitivities */
	}

      
      /* momentum wrt. species concentration */
      var = MASS_FRACTION;
      if (pd->v[pg->imtrx][var]){
	for (w=0; w<pd->Num_Species_Eqn; w++)
	  {
	    for ( n=0; n<ei[pg->imtrx]->dof[var]; n++) {
	      
	      for (kdir = 0; kdir < dim; kdir++)
		{
		  rotated_jacobian_conc[kdir][w][n] = 0.;
		}
	      
	      for (kdir = 0; kdir < dim; kdir++)
		{
		  if (rotation[I][eq][kdir]->ok) {
		    eqn=R_MOMENTUM1+kdir;
		    for (ldir = 0; ldir < dim; ldir++)
		      {
			rotated_jacobian_conc[kdir][w][n] +=
			  rotation[I][eq][kdir]->vector[ldir]*
			  lec->J[upd->ep[pg->imtrx][R_MOMENTUM1+ldir]][MAX_PROB_VAR + w][id][n];
		      }
		  }
		}
	      
	    } /* end of loop over nodes */
	  } /* end of loop over species */
	
	/* reinject d_mesh/d_mass back into lec-J for global assembly */
	for (w=0; w<pd->Num_Species_Eqn; w++)
	  {
	    for (kdir = 0; kdir < dim; kdir++)
	      {
		if (rotation[I][eq][kdir]->ok) {
		  /* loop over sensitivities */
		  
		  pvar = MAX_PROB_VAR + w;
		  
		  for ( n=0; n< ei[pg->imtrx]->dof[var]; n++ ) 
		    {
		      
		      lec->J[upd->ep[pg->imtrx][R_MOMENTUM1 + kdir]][pvar][id][n]
			= rotated_jacobian_scalar[kdir][n];
		    }
		}
	      }
	  }

      } /* end of if variable */
      
    } /* end of if Newton */

  /* Put rotated residual back into lec for scattering into global matrix */
  /* Note this is the last thing we do so our chain rule still works for
     mesh derivatives wrt rotation vector */
  for (kdir = 0; kdir < dim; kdir++)
    {
      if (rotation[I][eq][kdir]->ok) {
	peq = upd->ep[pg->imtrx][R_MOMENTUM1+kdir];
	lec->R[peq][id] = rotated_resid[kdir];
      }
    } 
  
} /* END of rotate_momentum_eqn */
/*************************************************************************** */
/*************************************************************************** */
/*ARGSUSED*/
void
rotate_res_jac_mom (
    int irow_index,                     /* Elemental stiffness matrix row index */
    int I,                              /* Global node number                   */
    int iconnect_ptr,  		        /* Pointer to beginning of connectivity list
                                            for the current element	         */
    int nodes_per_elem,                 /* number of nodes in the element       */
    int ielem_surf_dim,                 /* physical dim of elem surface - 0,1,2 */
    double snormal[MAX_PDIM],           /* vector components of surface normal  */
    double dsnormal_dx[MAX_PDIM][MAX_PDIM][MDE],
                                         /* Array of derivatives ([i][j][k]) of
                                            surface normal: component-i
                                            wrt displacement-j at node-k         */
    double stangent[2][MAX_PDIM],        /* vector components of 2 mutually orthogonal
                                              surface tangent vectors            */
    double dstangent_dx[2][MAX_PDIM][MAX_PDIM][MDE],
                                         /* Array of derivatives ([i][j][k]) of surface
                                            tangent vector 1 or 2: component-i
                                            wrt displacement-j at node-k         */
    double sign )                        /* sign of tangent vector               */

     
/* 
 * Function which corrects the global residual vector "resid_vect" and
 * the global Jacobian vector "a" so that the vector momentum equations for surface nodes 
 * are projected into a normal and tangential coordinate system.
 *
 *       Author:          P. R. Schunk (1511)
 *       Revisor:         R. R. Rao  20 July 1999 (lec revolution)
 */

{
/* LOCAL VARIABLES */
  int       var;
  int       id, j, n, Id, ShapeVar, ldof;
  int       peq, pvar;
  int       kdir, ldir, w, jvar;
  unsigned int       siz;
  double    svector[MAX_PDIM][MAX_PDIM];
  /* set of vectors corresponding to normal (1st vector) and tangents
     at the surface - for convenience in the following loops */
  double    dsvector_dx[MAX_PDIM][MAX_PDIM][MAX_PDIM][MDE];
  /* sensitivity of surface vectors with respect to nodal positions */

  double    rotated_resid[MDE];
  double    rotated_jacobian_vector[MAX_PDIM][MAX_PDIM][MDE];
  double    rotated_jacobian_scalar[MAX_PDIM][MDE];
  double    rotated_jacobian_conc[MAX_PDIM][MAX_CONC][MDE];
  
  /***************************** execution begins **********************************/
  
  ShapeVar = pd->ShapeVar;
  
  
  /* Correct residual equation first at local node "irow_index" or global node "I" */
  /*                Rx -> Rn    and Ry -> Rt                                       */
  /*       i.e.,    Rn = nx*Rx + ny*Ry + nz*Rz                                     */
  /*                Rt1 = t1x*Rx + t1y*Ry + t1z*Rz                                 */
  /*                Rt2 = t2x*Rx + t2y*Ry + t2z*Rz                                 */
  
  /* Load up values of surface vectors */
  /* first surface vector is the normal 
     - i.e. first mesh equation is the normal component */
  /* second (and third) surface vector(s) is(are) the tangent(s)
     - i.e. second and third mesh equations are the tangential components */
  
  switch (ielem_surf_dim+1) 
    {
      
    case 1:
      svector[0][0] = snormal[0];
      for ( id=0; id<ei[pg->imtrx]->dof[ShapeVar]; id++ )
        {
          Id = Proc_Elem_Connect[iconnect_ptr + id];
          ldof = ei[pg->imtrx]->ln_to_dof[ShapeVar][id];
          if (Dolphin[pg->imtrx][I][MESH_DISPLACEMENT1] > 0 && Dolphin[pg->imtrx][Id][MESH_DISPLACEMENT1] > 0 )
            {
              dsvector_dx[0][0][0][ldof] = dsnormal_dx[0][0][ldof];
            } /*end of Baby_Dolphin[pg->imtrx] */
        }
      break;
      
    case 2:
      for (kdir = 0; kdir < ielem_surf_dim+1; kdir++)
	{
          svector[0][kdir] = snormal[kdir];
          svector[1][kdir] = sign * stangent[0][kdir];
          for (ldir = 0; ldir < ielem_surf_dim+1; ldir++)
            {
              for ( id=0; id<ei[pg->imtrx]->dof[ShapeVar]; id++ )
                {
                  Id = Proc_Elem_Connect[iconnect_ptr + id];
                  ldof = ei[pg->imtrx]->ln_to_dof[ShapeVar][id];
                  if (Dolphin[pg->imtrx][I][MESH_DISPLACEMENT1] > 0 && Dolphin[pg->imtrx][Id][MESH_DISPLACEMENT1] > 0 )
                    {
                      dsvector_dx[0][kdir][ldir][ldof] = dsnormal_dx[kdir][ldir][ldof];
                      dsvector_dx[1][kdir][ldir][ldof] = sign * dstangent_dx[0][kdir][ldir][ldof];
                    } /*end of Baby_Dolphin[pg->imtrx] */
                }
            }
        }
      break;
      
    case 3:
      for (kdir = 0; kdir < ielem_surf_dim+1; kdir++)
	{
          svector[0][kdir] = snormal[kdir];
          svector[1][kdir] = sign * stangent[0][kdir];
          svector[2][kdir] = sign * stangent[1][kdir];
          for (ldir = 0; ldir < ielem_surf_dim+1; ldir++)
            {
              for ( id=0; id<ei[pg->imtrx]->dof[ShapeVar]; id++ )
                {
                  Id = Proc_Elem_Connect[iconnect_ptr + id];
                  ldof = ei[pg->imtrx]->ln_to_dof[ShapeVar][id];
                  if (Dolphin[pg->imtrx][I][MESH_DISPLACEMENT1] > 0 && Dolphin[pg->imtrx][Id][MESH_DISPLACEMENT1] > 0 )
                    {

                      dsvector_dx[0][kdir][ldir][ldof] = dsnormal_dx[kdir][ldir][ldof];
                      dsvector_dx[1][kdir][ldir][ldof] = sign * dstangent_dx[0][kdir][ldir][ldof];
                      dsvector_dx[2][kdir][ldir][ldof] = sign * dstangent_dx[1][kdir][ldir][ldof];
                    } /*end of Baby_Dolphin[pg->imtrx] */
                }
            }
        }
      break;
      
    default:
      exit(-1);
      break;
    }
  
  /* Correct residual equation first at local node "irow_index" or global node "I" */
  /*                Rx -> Rn    and Ry -> Rt                                       */
  /*       i.e.,    Rn = nx*Rx + ny*Ry + nz*Rz                                     */
  /*                Rt1 = t1x*Rx + t1y*Ry + t1z*Rz                                 */
  /*                Rt2 = t2x*Rx + t2y*Ry + t2z*Rz                                 */
  
  /* Now add on projection into n-t space */
  
  for (kdir = 0; kdir < ielem_surf_dim+1; kdir++)
    {
      rotated_resid[kdir] = 0.;
      for (ldir = 0; ldir < ielem_surf_dim+1; ldir++)
	{
	  rotated_resid[kdir] +=
	    svector[kdir][ldir]*lec->R[upd->ep[pg->imtrx][R_MOMENTUM1 + ldir]][irow_index];
	}
      
    } /* end of loop over direction */
  /*                                                                      */
  /*   Now correct Jacobian                                               */
  
  if (af->Assemble_Jacobian) 
    {
      siz = sizeof(double)* DIM*DIM*MDE;
      memset( rotated_jacobian_vector,0,siz);    
      
      /* Now correct for rotation */
      for (kdir = 0; kdir < ielem_surf_dim+1; kdir++)
	{
	  /* loop over sensitivities */
	  for ( j=0; j<ielem_surf_dim+1; j++)
	    {
	      
	      var=MESH_DISPLACEMENT1+j;
	      if (pd->v[pg->imtrx][var])
		{
		  pvar = upd->vp[pg->imtrx][var];
		  for ( n=0; n< ei[pg->imtrx]->dof[var]; n++ ) 
		    {
		      
		      for (ldir = 0; ldir < ielem_surf_dim+1; ldir++)
			{
			  
			  rotated_jacobian_vector[kdir][j][n]+=
			    svector[kdir][ldir]*lec->J[upd->ep[pg->imtrx][R_MOMENTUM1 + ldir]][pvar][irow_index][n];
			}
		    }
		  
		  for ( id=0; id<ei[pg->imtrx]->dof[ShapeVar]; id++ ) 
		    {
		      Id = Proc_Elem_Connect[iconnect_ptr + id];
		      ldof = ei[pg->imtrx]->ln_to_dof[ShapeVar][id];
		      if (Dolphin[pg->imtrx][I][MESH_DISPLACEMENT1] > 0 && Dolphin[pg->imtrx][Id][MESH_DISPLACEMENT1] > 0  )
			{ 
			  
			  for (ldir = 0; ldir < ielem_surf_dim+1; ldir++)
			    {
			      rotated_jacobian_vector[kdir][j][ldof] +=  
				dsvector_dx[kdir][ldir][j][ldof]
				*lec->R[upd->ep[pg->imtrx][R_MOMENTUM1+ldir]][irow_index];
			    }
			} /*end of Baby_Dolphin[pg->imtrx] */
		    }	  
		}
	    } /* end of loop over sensitivities */
	}
	
      /* reinject back into lec-J for global assembly */
      for (kdir = 0; kdir < ielem_surf_dim+1; kdir++)
	{
	  /* loop over sensitivities */
	  for ( j=0; j<ielem_surf_dim+1; j++)
	    {
	      
	      var=R_MESH1+j;
	      pvar = upd->vp[pg->imtrx][var];
	      
	      for ( n=0; n< ei[pg->imtrx]->dof[var]; n++ ) 
		{
		  lec->J[upd->ep[pg->imtrx][R_MOMENTUM1 + kdir]][pvar][irow_index][n]
		    = rotated_jacobian_vector[kdir][j][n];
		}
	      
	    } /* end of loop over sensitivities */
	}
      
      /* momentum wrt. pressure */
      var = PRESSURE;
      if (pd->v[pg->imtrx][var]){
	pvar = upd->vp[pg->imtrx][var];
	for ( n=0; n<ei[pg->imtrx]->dof[var]; n++) {
	  
	  for (kdir = 0; kdir < ielem_surf_dim+1; kdir++)
	    {
	      rotated_jacobian_scalar[kdir][n] = 0.;
	    }
	  
	  for (ldir = 0; ldir < ielem_surf_dim+1; ldir++)
	    {
	      peq = upd->ep[pg->imtrx][R_MOMENTUM1+ldir];
	      for (kdir = 0; kdir < ielem_surf_dim+1; kdir++)
		{
		  rotated_jacobian_scalar[kdir][n] += 
		    svector[kdir][ldir]*lec->J[peq][pvar][irow_index][n];
		}
	    }
	  /*reinject back into lec-J  for global assembly */
	  for (ldir = 0; ldir < ielem_surf_dim+1; ldir++)
	    {
	      peq = upd->ep[pg->imtrx][R_MOMENTUM1+ldir];
	      lec->J[peq][pvar][irow_index][n] = rotated_jacobian_scalar[ldir][n];
	    }
	} /* end of loop over nodes */
      } /* end of if variable */
      
      /* momentum wrt. temperature */
      var = TEMPERATURE;
      if (pd->v[pg->imtrx][var]){
	pvar = upd->vp[pg->imtrx][var];
	for ( n=0; n<ei[pg->imtrx]->dof[var]; n++) {
	  for (kdir = 0; kdir < ielem_surf_dim+1; kdir++)
	    {
	      rotated_jacobian_scalar[kdir][n] = 0.;
	    }
	  
	  for (ldir = 0; ldir < ielem_surf_dim+1; ldir++)
	    {
	      peq = upd->ep[pg->imtrx][R_MOMENTUM1+ldir];
	      for (kdir = 0; kdir < ielem_surf_dim+1; kdir++)
		{
		  rotated_jacobian_scalar[kdir][n] += 
		    svector[kdir][ldir]*lec->J[peq][pvar][irow_index][n];
		}
	    }
	  /*reinject back into lec-J for frontal solver */
	  for (ldir = 0; ldir < ielem_surf_dim+1; ldir++)
	    {
	      peq = upd->ep[pg->imtrx][R_MOMENTUM1+ldir];
	      lec->J[peq][pvar][irow_index][n] = rotated_jacobian_scalar[ldir][n];
	    }
	} /* end of loop over nodes */
      } /* end of if variable */
      
#ifdef COUPLED_FILL
      /* momentum wrt. temperature */
      var = FILL;
      if (pd->v[pg->imtrx][var]){
	pvar = upd->vp[pg->imtrx][var];
	for ( n=0; n<ei[pg->imtrx]->dof[var]; n++) {
	  for (kdir = 0; kdir < ielem_surf_dim+1; kdir++)
	    {
	      rotated_jacobian_scalar[kdir][n] = 0.;
	    }
	  
	  for (ldir = 0; ldir < ielem_surf_dim+1; ldir++)
	    {
	      peq = upd->ep[pg->imtrx][R_MOMENTUM1+ldir];
	      for (kdir = 0; kdir < ielem_surf_dim+1; kdir++)
		{
		  rotated_jacobian_scalar[kdir][n] += 
		    svector[kdir][ldir]*lec->J[peq][pvar][irow_index][n];
		}
	    }
	  /*reinject back into lec-J for frontal solver */
	  for (ldir = 0; ldir < ielem_surf_dim+1; ldir++)
	    {
	      peq = upd->ep[pg->imtrx][R_MOMENTUM1+ldir];
	      lec->J[peq][pvar][irow_index][n] = rotated_jacobian_scalar[ldir][n];
	    }
	} /* end of loop over nodes */
      } /* end of if variable */
#endif /* COUPLED_FILL */      
      
      /* momentum wrt. velocity */
      for (jvar = 0; jvar < ielem_surf_dim+1; jvar++)
	{
	  var = VELOCITY1 + jvar;
	  if (pd->v[pg->imtrx][var]){
	    pvar = upd->vp[pg->imtrx][var];
	    for ( n=0; n<ei[pg->imtrx]->dof[var]; n++) {
	      
	      for (kdir = 0; kdir < ielem_surf_dim+1; kdir++)
		{
		  rotated_jacobian_vector[kdir][jvar][n] = 0.;
		}
	      
	      for (ldir = 0; ldir < ielem_surf_dim+1; ldir++)
		{
		  peq = upd->ep[pg->imtrx][R_MOMENTUM1+ldir];
		  for (kdir = 0; kdir < ielem_surf_dim+1; kdir++)
		    {
		      rotated_jacobian_vector[kdir][jvar][n] +=
			svector[kdir][ldir]*lec->J[peq][pvar][irow_index][n];
		    }
		}
	      /*reinject back into lec-J for global assembly */
	      for (ldir = 0; ldir < ielem_surf_dim+1; ldir++)
		{
		  peq = upd->ep[pg->imtrx][R_MOMENTUM1+ldir];
		  lec->J[peq][pvar][irow_index][n] = rotated_jacobian_vector[ldir][jvar][n];
		}
	      
	    } /* end of loop over nodes */
	  } /* end of if variable */
	} /* end of loop over jvar direction */
      
      /* momentum wrt. species concentration */
      var = MASS_FRACTION;
      if (pd->v[pg->imtrx][var]){
	for (w=0; w<pd->Num_Species_Eqn; w++)
	  {
	    for ( n=0; n<ei[pg->imtrx]->dof[var]; n++) {
	      
	      for (kdir = 0; kdir < ielem_surf_dim+1; kdir++)
		{
		  rotated_jacobian_conc[kdir][w][n] = 0.;
		}
	      
	      for (ldir = 0; ldir < ielem_surf_dim+1; ldir++)
		{
		  peq = upd->ep[pg->imtrx][R_MOMENTUM1+ldir];
		  for (kdir = 0; kdir < ielem_surf_dim+1; kdir++)
		    {
		      
		      rotated_jacobian_conc[kdir][w][n] += 
			svector[kdir][ldir]*lec->J[peq][MAX_PROB_VAR + w][irow_index][n];
		    }
		}
	      /*reinject back into lec-J for frontal solver */
	      for (ldir = 0; ldir < ielem_surf_dim+1; ldir++)
		{
		  peq = upd->ep[pg->imtrx][R_MOMENTUM1+ldir];
		  lec->J[peq][MAX_PROB_VAR + w][irow_index][n] = rotated_jacobian_conc[ldir][w][n];
		}
	    } /* end of loop over nodes */
	  }
      } /* end of if variable */
      
    } /* end of if Newton */
  
  /* Put rotated residual back into lec for scattering into global matrix */
  for (kdir = 0; kdir < ielem_surf_dim+1; kdir++)
    {
      lec->R[ upd->ep[pg->imtrx][R_MOMENTUM1+kdir] ][irow_index] = rotated_resid[kdir];
    } 
  
  
} /* END of rotate_res_jac_mom */
/*****************************************************************************/

static int rotation_allocated = FALSE;


void
calculate_all_rotation_vectors (Exo_DB *exo,		/* the mesh */
				double x[] )		/* Soln vector */

/*
 * Routine for determining which equations require which vectors for rotation of the 
 * bulk equations (prior to applying strong boundary conditions).  This is done ahead 
 * of time to avoid some annoying problems that arise from unstructured grids and 
 * multiple materials
 *
 * Written by Richard Cairncross, 3 September 1996
 */

{
/* TAB certifies that this function conforms to the exo/patran side numbering convention 11/10/98. */

  int i, j, k, I, J, p, q, b, num_total_nodes, irc, eq=-1;

  /* 
   * The next sets of declarations are for the various types of vectors which 
   * can be used for rotations - these declarations are just place-holders which will
   * be "pointed" to the correct place before the values are calculated
   */
  struct Rotation_Vectors *normal, *normal2, *normal3, *tangent1, *tangent2;
  struct Rotation_Vectors *line_tangent, *binormal;
  /* these dummy arrays allow calculation of some of the vectors that are needed for
   * the real rotation vectors */
  struct Rotation_Vectors dum_vect[DIM], vector[DIM];

  int ielem, iconnect_ptr, ielem_type; 
  int num_local_nodes, ielem_dim, dim, err, id_side, num_nodes_on_side;
  int num_ROT_nodes, *elem_node_id=NULL;
  int id, j_id, num_rots_finished, j_new;
  double xi[DIM];               /* Local element coordinates of Gauss point. */
  int num_nodes_on_edge=-1;
  int local_edge_node_list[MAX_NODES_PER_SIDE];
  int edge_elem_node_id[MAX_NODES_PER_SIDE];
  int local_side_node_list[MAX_NODES_PER_SIDE];
  int side_elem_node_id[MAX_NODES_PER_SIDE];
  int param_dir, id_edge=-1, v_id, reseed, kdir;
  double seed[DIM];

  /***************************************************************************/
  /* BEGIN EXECUTION */
  if (Debug_Flag > 0) DPRINTF(stderr, "Starting to calculate rotation vectors\n");

  /* initialize the rotation_vector array of rotation structures */
  num_total_nodes = Num_Internal_Nodes + Num_Border_Nodes + Num_External_Nodes;
  dim = pd_glob[0]->Num_Dim;

  if ( ! rotation_allocated )
    {
      rotation = (struct Rotation_Vectors ****) 
	smalloc(num_total_nodes*sizeof(struct Rotation_Vectors ***));

      for ( i=0; i<num_total_nodes; i++)
	{
	  rotation[i] = (struct Rotation_Vectors ***)
	    smalloc(NUM_VECTOR_EQUATIONS*sizeof(struct Rotation_Vectors **));
	  for ( j=0; j<NUM_VECTOR_EQUATIONS; j++)
	    {
	      rotation[i][j] = NULL;

	      if ( ROT_list[i] != NULL && ROT_list[i][j] != -1)
		{
		  rotation[i][j] = (struct Rotation_Vectors **) 
		    smalloc(dim*sizeof(struct Rotation_Vectors*));
		  for ( k=0; k<dim; k++)
		    {
		      rotation[i][j][k] = (struct Rotation_Vectors*)
			smalloc(sizeof(struct Rotation_Vectors));
		    }
		}
	    }
	}
      rotation_allocated = TRUE;
    }

  /* INITIALIZE */
  /* loop over nodes */
  for(I=0; I<num_total_nodes; I++) {
    for (j=0; j<NUM_VECTOR_EQUATIONS; j++) {
      if (rotation[I][j] != NULL){
	/* INITIALIZE */
	for (p=0; p<dim; p++) {
	  rotation[I][j][p]->ok = 0;
	  rotation[I][j][p]->d_vector_n = 0;
	  for (j_id=0; j_id<MNROT; j_id++) {
	    rotation[I][j][p]->d_vector_J[j_id] = -1;
	  }
	  for (q=0; q<dim; q++) {
	    rotation[I][j][p]->vector[q] = 0.;
	    for (b=0; b<dim; b++) {
	      for (j_id=0; j_id<MNROT; j_id++) {
		rotation[I][j][p]->d_vector_dx[q][b][j_id] = 0.;
	      }
	    }
	  }
	}
      }
    }
  } /* end of INITIALIZATION */

  /* LOOP over rotation conditions in input deck */
  for (irc = 0; irc < Num_ROT; irc++) {
    if (ROT_Types[irc].eq_type == R_MESH1) eq = VECT_EQ_MESH;
    if (ROT_Types[irc].eq_type == R_MOMENTUM1) eq = VECT_EQ_MOM;

    /* quick escape if no rotation needed */
    if (ROT_Types[irc].ROTATE) {

      if (ROT_Types[irc].elems == NULL) {
	/* Load list of elements to which the ROTATION condition applies*/
	get_ss_element_list(ROT_Types[irc].ss_id[0], ROT_Types[irc].ss_id[1], 
			    ROT_Types[irc].ss_id[2], &ROT_Types[irc].num_elem, 
			    &ROT_Types[irc].elems, &ROT_Types[irc].ss_ptr, 
			    eq, irc, exo);
      }

      for(i=0; i< ROT_Types[irc].num_elem; i++) {         
	ielem = ROT_Types[irc].elems[i];

        if (x == x_static) /* be the least disruptive possible */
          {
            err = load_elem_dofptr(ielem, exo, x_static, x_old_static,
                                   xdot_static, xdot_old_static, 0);
          }
        else
          {
            err = load_elem_dofptr(ielem, exo, x, x, x, x, 0);
          }
	err = bf_mp_init(pd);

	iconnect_ptr = ei[pg->imtrx]->iconnect_ptr;
	ielem_type   = ei[pg->imtrx]->ielem_type;
	num_local_nodes = ei[pg->imtrx]->num_local_nodes;
	ielem_dim       = ei[pg->imtrx]->ielem_dim;
	dim             = ielem_dim;


	/* find SIDE info for primary side */
	num_nodes_on_side = count_nodes_on_SS(ROT_Types[irc].ss_id[0], -1, -1, 
					      ei[pg->imtrx]->iconnect_ptr, ielem, 
					      num_local_nodes,
					      local_side_node_list, 
					      side_elem_node_id);

	id_side = find_id_side (ei[pg->imtrx]->ielem, num_nodes_on_side, 
				local_side_node_list, side_elem_node_id, exo);
	/* Calculates the ID side correctly for tets */
	int current_id; 

	current_id = in_list(ROT_Types[irc].ss_id[0], 0, exo->num_side_sets, exo->ss_id);
	if ( ei[pg->imtrx]->ielem_type == LINEAR_TET ) id_side = find_id_side_SS(ei[pg->imtrx]->ielem, current_id, exo);

	/* find EDGE info for primary edge */
	if (ROT_Types[irc].type == CURVE || ROT_Types[irc].type == VERTEX) {
	  num_nodes_on_edge = count_nodes_on_SS(ROT_Types[irc].ss_id[0], ROT_Types[irc].ss_id[1],
						-1, 
						ei[pg->imtrx]->iconnect_ptr, ielem, num_local_nodes,
						local_edge_node_list, edge_elem_node_id);

	  if ( num_nodes_on_edge > 0 )
	    {
	      param_dir = -1;
	      id_edge = find_id_edge(ielem, num_nodes_on_edge, 
				     local_edge_node_list, edge_elem_node_id, 
				     &param_dir, exo);

	      if ( ei[pg->imtrx]->ielem_type == LINEAR_TET ) id_edge = find_id_edge_TET(ielem, num_nodes_on_edge, 
				     local_edge_node_list, edge_elem_node_id, 
				     &param_dir, exo);
	    }
	}

	/* 
	 * LOOP over NODES to which this condition applies
	 */

	num_ROT_nodes = 0;

	
	if (ROT_Types[irc].type == FACE) {
	  num_ROT_nodes = num_nodes_on_side;
	  elem_node_id = &(side_elem_node_id[0]);
	} else if (ROT_Types[irc].type == CURVE) {
	  num_ROT_nodes = num_nodes_on_edge;
	  elem_node_id = &(edge_elem_node_id[0]);
	} else if (ROT_Types[irc].type == VERTEX &&
		   num_nodes_on_edge > 0 ) {
	  num_ROT_nodes = 1;
	  /* determine id of vertex node */
	  v_id = 0;
	  while(Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + v_id] != ROT_Types[irc].node) v_id++;
	  elem_node_id = &(v_id);
	  
	}

	  /* use nodal points only!! */
	  for ( k=0; k<num_ROT_nodes; k++) {
	    /* Find the local element node number for the current node */
	    id = elem_node_id[k];
	    I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + id];

	    if (ROT_list[I] != NULL && ROT_list[I][eq] == irc) {

	      /* set up pointers to spaces in the rotation structure array */
	      num_rots_finished = set_pointers_to_vectors(&normal, &tangent1, &tangent2, 
				      &normal2, &normal3, &line_tangent, &binormal,
				      dum_vect, vector, I, eq, dim, irc);

	      /* make sure we still need to calculate rotation vectors */
	      if (num_rots_finished < dim) {
		find_nodal_stu(id, ielem_type, &xi[0], &xi[1], &xi[2]);
		err = load_basis_functions( xi, bfd);
		EH( err, "problem from load_basis_functions");
		err = beer_belly();
		EH( err, "beer_belly");
		err = load_fv();
		EH( err, "load_fv");
		err = load_bf_grad();
		EH( err, "load_bf_grad");
		err = load_bf_mesh_derivs(); 
		EH( err, "load_bf_mesh_derivs");

		/* put NORMAL vector into array */
		if (normal != NULL) {
		  /* calculate the determinant of the surface jacobian  and the normal to 
		   * the surface all at one time */
		    surface_determinant_and_normal (ei[pg->imtrx]->ielem, iconnect_ptr, num_local_nodes, 
						    ielem_dim - 1, id_side, num_nodes_on_side,
						    side_elem_node_id);

		  for (p=0; p<dim; p++) {
		    normal->vector[p] = fv->snormal[p];
		    normal->ok = 1;
		    if (af->Assemble_Jacobian && pd->v[pg->imtrx][R_MESH1]) {
		      normal->d_vector_n = num_nodes_on_side;
		      for ( j=0; j<num_nodes_on_side; j++) {
			/* Find the local element node number for the current node */
			j_id = side_elem_node_id[j];
			J = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + j_id];
			normal->d_vector_J[j] = J;
			for (q=0; q<dim; q++) {
			  normal->d_vector_dx[p][q][j] = fv->dsnormal_dx[p][q][j_id];
			}
		      }
		    }
		  }
		} else EH(-1, "can't calculate normal vector");

		do_LSA_mods(LSA_SURFACE);

		reseed = 0;
		if (ROT_Types[irc].method == ROT_BASIS_RESEED) {
		  /* need some special logic to use previously calculated
		   * tangent-along-basis to reseed the calculation of
		   * the new tangent - this allows averaging of tangents
		   * in adjacent elements when their basis vectors do 
		   * not line up */

		  /* determine which direction holds tangent1 */
		  kdir = -1;
		  for (p=0; p<dim; p++) if (ROT_Types[irc].BC_Type[p] == ROT_T1) kdir = p;
		  EH(kdir,"cant find tangent1");

		  /* check to see if tangent1 has already been calculated */
		  if (rotation[I][eq][kdir]->ok) {
		    reseed = 1;
		    /* if it has, load up seed vector */
		    for (q=0; q<dim; q++) {
		      seed[q] = rotation[I][eq][kdir]->vector[q];
		    }
		  }
		}

		/* calculate SEEDED tangent, if requested */
		if ( tangent1 != NULL || tangent2 != NULL ) {
		  if (ROT_Types[irc].method == ROT_SEED ) {
		    calc_tangent_from_seed (tangent1, normal, ROT_Types[irc].seed, dim);
		  } else if (ROT_Types[irc].method == ROT_BASIS_RESEED && reseed) {
		    calc_tangent_from_seed (tangent1, normal, seed, dim);
		  }
		} /* end of SEEDED tangent */

		/* calculate UNSEEDED tangent, if required */
		if (tangent1 != NULL || tangent2 != NULL) { 
		  if (ROT_Types[irc].method == ROT_BASIS || 
		     ROT_Types[irc].method == ROT_BASIS_ONCE || 
		     (ROT_Types[irc].method == ROT_BASIS_RESEED && !reseed)) {

		    calc_tangent_along_basis (tangent1, normal, dim, 
					      id_side, num_nodes_on_side, side_elem_node_id);
		  }
		} /* end of UNSEEDED tangent */

		/* calculate SECOND tangnt = N X T1, if required */
		if (tangent2 != NULL) {

		  cross_vectors (tangent2, tangent1, normal, dim);

		} /* end of SECOND tangent */

		/* Find LINE_TANGENT along EDGE */
		if (line_tangent != NULL || binormal != NULL) {
		  
		  if ( ei[pg->imtrx]->ielem_type == LINEAR_TET )
		    {
		     calc_unseeded_edge_tangents_TET (line_tangent, 
					       ei[pg->imtrx]->iconnect_ptr, dim,  
					       id_side, 
					       id_edge, num_nodes_on_edge, 
					       edge_elem_node_id,
					       param_dir);
		    }
		  else
		    {
		      calc_unseeded_edge_tangents (line_tangent, 
					       ei[pg->imtrx]->iconnect_ptr, dim,  
					       id_side, 
					       id_edge, num_nodes_on_edge, 
					       edge_elem_node_id,
					       param_dir);
		    }

		}

		/* Find BINORMAL along EDGE */
		if (binormal != NULL) {
		  cross_vectors (binormal, normal, line_tangent, dim);
		}

		/* ADD the results into rotation vectors */
		for (p=0; p<dim; p++) {
		  if (vector[p].ok ) {
		    if (ROT_Types[irc].method == ROT_BASIS_ONCE &&  rotation[I][eq][p]->ok) {
		      /* Don't add in the new vector */
		    } else {
		      rotation[I][eq][p]->ok++;
		      for (q=0; q<dim; q++) {
			rotation[I][eq][p]->vector[q] += vector[p].vector[q];
			for (b=0; b<dim; b++) {
			  for (j_id=0; j_id<vector[p].d_vector_n; j_id++) {
			    /* determine if global sensitivity is in list already
			     *  or add it to list */
			    if ((j_new = in_list(vector[p].d_vector_J[j_id], 0, 
						 rotation[I][eq][p]->d_vector_n, 
						 rotation[I][eq][p]->d_vector_J))
				== -1) {
			      if (rotation[I][eq][p]->d_vector_n >= MNROT) {
				/* illegal number of sensitivities */
				EH(-1, "Increase the number of available rotation sensitivities (MNROT)");
			      }
			      j_new = rotation[I][eq][p]->d_vector_n++;
			      rotation[I][eq][p]->d_vector_J[j_new] = 
				vector[p].d_vector_J[j_id];
			    }
			  
			    rotation[I][eq][p]->d_vector_dx[q][b][j_new] +=
			      vector[p].d_vector_dx[q][b][j_id];
			  }
			}
		      }
		    }
		  }
		}

	      }
	    } /* end of nodal rotation check */
	  } /* end of loop over local nodes, k */
	} /* end of loop over side elements, i */

 
    } /* end of quick escape */
  } /* end of loop over rotation conditions, irc */

  /* ADJUST magnitude for multiple vectors */
  /* loop over nodes */
  for(I=0; I<num_total_nodes; I++) {
    for (j=0; j<NUM_VECTOR_EQUATIONS; j++) {
      if (rotation[I][j] != NULL){
	/* ADJUST */
	for (p=0; p<dim; p++) {
	  if (rotation[I][j][p]->ok > 1) {
	    for (q=0; q<dim; q++) {
	      rotation[I][j][p]->vector[q] /= (double) rotation[I][j][p]->ok;
	      for (b=0; b<dim; b++) {
		for (j_id=0; j_id<rotation[I][j][p]->d_vector_n; j_id++) {
		  rotation[I][j][p]->d_vector_dx[q][b][j_id] /= (double) rotation[I][j][p]->ok;
		  
		}
	      }
	    }
	    rotation[I][j][p]->ok = 1;
	  }
	}
      }
    }
  } /* end of ADJUST */

  if (Debug_Flag > 0) DPRINTF(stderr, "Done calculating rotation vectors\n");

  return;
} /* END of calculate_all_rotation_vectors  */




void
calculate_2D_rotation_vectors (Exo_DB *exo,		/* the mesh */
			       double x[] )		/* Soln vector */
					
     /*
      *   Routine for determining the rotation vectors to be applied at nodes for 2D problems
      *   Normally, for 2D problems rotation vectors can be determined "on the fly" because of the simplicity and localization
      *   of the interactions.
      *   However, for meshes that involve triangular elements, it is less easy to handle the rotation vectors on an element by element
      *   approach since there are likely to be element node residuals that must be rotated, but do not have a corresponding element face 
      *   on the boundary from which to determine the appropriate rotation vector
      *   Thus, it is necessary to determine a set of globally-known rotation vectors for each node at which a rotated condition is 
      *   applied.  This routine adapts the rotation structure developed for three dimensional problems to two dimensional problems.
      */ 
{
  int i, j, jvar, j_id,q,b;
  int inode, ssid, iss;
  int elem, id_side,lnn;
  int num_total_nodes;
  int dim;
  int mesh_count=0, mom_count=0;
  static int **local_ROT_list ;
  int *gnn_side_list, local_node_side_list[MAX_NODES_PER_SIDE];
  int num_nodes_on_side, num_local_nodes, iconnect_ptr;
	
  num_total_nodes = Num_Internal_Nodes + Num_Border_Nodes + Num_External_Nodes;
  dim = pd_glob[0]->Num_Dim;


  if ( ! rotation_allocated )
    {
      rotation = (struct Rotation_Vectors ****) 
	smalloc(num_total_nodes*sizeof(struct Rotation_Vectors ***));
      local_ROT_list = (int **) alloc_ptr_1(num_total_nodes);
		  
      for ( i=0; i<num_total_nodes; i++)
	{
	  rotation[i] = (struct Rotation_Vectors ***)
	    smalloc(NUM_VECTOR_EQUATIONS*sizeof(struct Rotation_Vectors **));
			  
	  local_ROT_list[i] = alloc_int_1(NUM_VECTOR_EQUATIONS,-1);
			  
	  for ( j=0; j<NUM_VECTOR_EQUATIONS; j++)
	    {
	      rotation[i][j] = NULL;
	    }
	}
      rotation_allocated = TRUE;
    }
	  
  while ( mesh_count < num_mesh_rotate[pg->imtrx] || mom_count < num_mom_rotate[pg->imtrx] )
    {
  
      if( mom_count < num_mom_rotate[pg->imtrx] )
	{
	  inode = mom_rotate_node[pg->imtrx][mom_count];
	  local_ROT_list[inode][0] = mom_rotate_ss[pg->imtrx][mom_count];
	}
		
      if( mesh_count < num_mesh_rotate[pg->imtrx] )
	{
	  inode = mesh_rotate_node[pg->imtrx][mesh_count];
	  local_ROT_list[inode][1] = mesh_rotate_ss[pg->imtrx][mesh_count];
	}
		
      mesh_count++; mom_count++;
    }
	
  for( inode=0; inode<num_total_nodes; inode++)
    {
      for( jvar=0; jvar<NUM_VECTOR_EQUATIONS; jvar++)
	{
	  if( rotation[inode][jvar] != NULL )
	    {
	      for (i=0; i<dim; i++) {
		rotation[inode][jvar][i]->ok = 0;
		rotation[inode][jvar][i]->d_vector_n = 0;
		for (j_id=0; j_id<MNROT; j_id++) {
		  rotation[inode][jvar][i]->d_vector_J[j_id] = -1;
		}
		for (q=0; q<dim; q++) {
		  rotation[inode][jvar][i]->vector[q] = 0.;
		  for (b=0; b<dim; b++) {
		    for (j_id=0; j_id<MNROT; j_id++) {
		      rotation[inode][jvar][i]->d_vector_dx[q][b][j_id] = 0.;
		    }
		  }
		}
	      }						
				
	    }
			
	  if( ( ssid = local_ROT_list[inode][jvar] ) != -1 )
	    {
	      if (  ( iss = in_list( ssid, 0, exo->num_side_sets, exo->ss_id ) ) == -1 ) EH( -1,"Error trap out");
				
	      for( i=0; i< exo->ss_num_sides[iss]; i++)
		{					
					
		  num_nodes_on_side = exo->ss_node_side_index[iss][i+1] - exo->ss_node_side_index[iss][i];
		  gnn_side_list = &( exo->ss_node_list[iss][exo->ss_node_side_index[iss][i]] );
					
		  if( (lnn = in_list( inode, 0, num_nodes_on_side, gnn_side_list )  ) != -1 )
		    {
		      double xi[3] = {0.,0.,0.};
		      int err;
		      int add_vectors = TRUE;
						
		      /* We found a side with this inode */
							
		      elem = exo->ss_elem_list[exo->ss_elem_index[iss]+i];
						
		      if (x == x_static) 
                        err = load_elem_dofptr(elem, exo, x_static, x_old_static, xdot_static, xdot_old_static, 0);
		      else
						
                        err = load_elem_dofptr(elem, exo, x, x, x, x, 0);
					
		      if( SS_Internal_Boundary[iss] != -1 )
			{
			  int curr_eb = Current_EB_ptr->Elem_Blk_Id;
			  int index = 0;
							
			  while( ss_to_blks[index+1][iss] != curr_eb && index <(MAX_MAT_PER_SS+1)) index++;
							
			  if( index == (MAX_MAT_PER_SS+1) ) EH(-1," WTH.  \n");
							
			  if ( (index % 2) != 0 ) add_vectors = FALSE;
						
			}


		      if( add_vectors ) 
			{
			  id_side = exo->ss_side_list[exo->ss_elem_index[iss]+i];
							
			  id_side = find_id_side( elem, num_nodes_on_side, gnn_side_list, local_node_side_list, exo) ;

			  if ( ei[pg->imtrx]->ielem_type == LINEAR_TET ) id_side = find_id_side_SS(ei[pg->imtrx]->ielem, iss, exo);
			  /*find_surf_center_st (Elem_Type(exo, elem) , id_side, pd->Num_Dim, xi, &s, &t);	*/						
							
			  num_local_nodes = elem_info(NNODES, Elem_Type(exo, elem) );
							
			  iconnect_ptr = exo->elem_node_pntr[elem];
							
			  /*lnn = local_node_side_list[lnn];*/
							
			  lnn = in_list( inode, 0, num_local_nodes, &Proc_Elem_Connect[iconnect_ptr]);
							
			  find_nodal_stu( lnn, Elem_Type(exo, elem), xi, xi+1, xi+2 );
							
			  err = load_basis_functions(xi, bfd);
			  EH(err, "problem from load_basis_functions");
							
			  err = beer_belly();
			  EH(err, "beer_belly");
							
			  err = load_coordinate_scales(pd->CoordinateSystem, fv);
			  EH(err, "load_coordinate_scales(fv)");
							
			  surface_determinant_and_normal( elem, iconnect_ptr, num_local_nodes, ei[pg->imtrx]->ielem_dim - 1, id_side, 
							  num_nodes_on_side, local_node_side_list );
							
			  calc_surf_tangent( elem, iconnect_ptr, num_local_nodes, ei[pg->imtrx]->ielem_dim-1,
					     num_nodes_on_side,local_node_side_list);
						
			  if( rotation[inode][jvar] == NULL ) 
			    {
			      rotation[inode][jvar] = (struct Rotation_Vectors **) smalloc(pd->Num_Dim*sizeof(struct Rotation_Vectors*));
			      rotation[inode][jvar][0] = NULL;
			      rotation[inode][jvar][1] = NULL;
			    }
							
			  append_vectors ( inode, num_nodes_on_side, iconnect_ptr, local_node_side_list, rotation[inode][jvar] );
			}
		      /*if( inode == 4 )
			{
			int vvec=1, comp=0;
			printf("I: %d, jvar: %d,  N(%7.4g, %7,4g), TN(%7.4g, %7,4g) \n", inode+1, jvar,  \
			rotation[inode][jvar][0]->vector[0],  rotation[inode][jvar][0]->vector[1], rotation[inode][jvar][1]->vector[0],  rotation[inode][jvar][1]->vector[1]);
			printf("\t %d   %7.4g : %d   %7.4g : %d   %7.4g \n", rotation[inode][jvar][vvec]->d_vector_J[0]+1, rotation[inode][jvar][vvec]->d_vector_dx[0][comp][0], \
			rotation[inode][jvar][vvec]->d_vector_J[1]+1, rotation[inode][jvar][vvec]->d_vector_dx[0][comp][1], \
			rotation[inode][jvar][vvec]->d_vector_J[2]+1, rotation[inode][jvar][vvec]->d_vector_dx[0][comp][2] );
			}
			if( inode == 1)										   
			printf("I: %d, elem: %d, side: %d,  N(%7.4g, %7,4g), TN(%7.4g, %7,4g) \n", inode, elem, id_side, fv->snormal[0], fv->snormal[1], fv->stangent[0][0], fv->stangent[0][1] );
		      */
		    }
		}
	    }
			
	  if( rotation[inode][jvar] != NULL )
	    {


	      /*printf("I: %d, jvar: %d,  N(%7.4g, %7,4g), TN(%7.4g, %7,4g) \n", inode+1, jvar,  \
		rotation[inode][jvar][0]->vector[0],  rotation[inode][jvar][0]->vector[1], rotation[inode][jvar][1]->vector[0],  rotation[inode][jvar][1]->vector[1]);

	      */
	      /*printf("\t %d   %7.4g : %d   %7.4g : %d   %7.4g \n", rotation[inode][jvar][1]->d_vector_J[0]+1, rotation[inode][jvar][1]->d_vector_dx[0][0][0], \
		rotation[inode][jvar][1]->d_vector_J[1]+1, rotation[inode][jvar][1]->d_vector_dx[0][0][1], \
		rotation[inode][jvar][1]->d_vector_J[2]+1, rotation[inode][jvar][1]->d_vector_dx[0][0][2] );
		}
		printf("I: %d, jvar: %d, d_vector_n[0]: %d, d_vector_n[1]: %d  \n", inode+1, jvar,  rotation[inode][jvar][0]->d_vector_n, rotation[inode][jvar][1]->d_vector_n  ); 
		for(i=0; i< rotation[inode][jvar][0]->d_vector_n; i++)  printf("\t %d ", rotation[inode][jvar][0]->d_vector_J[i]+1, rotation[inode][jvar][0]->d_vector_dx[ ) ; printf("\n");
		for(i=0; i< rotation[inode][jvar][1]->d_vector_n; i++)  printf("\t %d ", rotation[inode][jvar][1]->d_vector_J[i]+1 ) ; printf("\n");*/

	    }
	}
    }
  return;
}

									
void append_vectors( int inode, 
					 int num_nodes_on_side,
					 int iconnect_ptr,
					 int *local_node_side_list,
					 ROTATION_VECTORS_STRUCT **vectors )
					 

{
	int i,j,k, vect,id, J, j_id;
	int dim = pd->Num_Dim;
	
	for( vect=0; vect< dim; vect++) {
		
		if ( vectors[vect] == NULL ) 
		{
			vectors[vect] = (struct Rotation_Vectors*)smalloc(sizeof(struct Rotation_Vectors));
			
			for(i=0;i<dim;i++) {
				vectors[vect]->vector[i] = 0;
				memset(vectors[vect]->d_vector_dx[i],0,DIM*MNROT*sizeof(double));
			}
			for(i=0; i<MNROT;i++) vectors[vect]->d_vector_J[i] = -1;
			vectors[vect]->d_vector_n=0;
			vectors[vect]->ok = 0;
		}
		
		for( i=0; i<dim;i++) {
			
			if( vect == 0 )
				vectors[vect]->vector[i] += fv->snormal[i];
			else
		       vectors[vect]->vector[i] += fv->stangent[0][i];	
		}
		vectors[vect]->ok++;
		
		if(af->Assemble_Jacobian && pd->e[pg->imtrx][R_MESH1] )
		{			
			for ( j=0; j<num_nodes_on_side; j++) {
				
				id = local_node_side_list[j];
				J= Proc_Elem_Connect[iconnect_ptr + id];
				
				j_id = in_list( J, 0, vectors[vect]->d_vector_n,   vectors[vect]->d_vector_J );
				
				if( j_id == -1 )
				{
					j_id = vectors[vect]->d_vector_n;
				
					vectors[vect]->d_vector_J[j_id] = J;
					
					vectors[vect]->d_vector_n++;
				}
				
					for( i=0; i<dim; i++) {
						for( k=0; k<dim; k++) {
						
						if( vect == 0 )
							vectors[vect]->d_vector_dx[i][k][j_id] +=  fv->dsnormal_dx[i][k][id];
						else
						vectors[vect]->d_vector_dx[i][k][j_id] += fv->dstangent_dx[0][i][k][id];

						
					}
				}
			}
		}
	}
	return;
}
	


/* get_ss_element_list -- load list of elements along 1,2,3 sidesets...
 *
 * Usage:
 *
 *	Given one meaningful SSID and two zeros for the second & third args,
 * this returns the list of elements on the sideset.
 *
 *      Given two meaningful SSID's and a zero for SSID3, it attempts to
 * return the intersection set. In 3D, for reasonable meshes, a list of the
 * elements along an edge is return. 
 *
 *	Given three meaningful SSID's, it attempts to locate the vertex element
 * that is the intersection of all three sidesets.
 *
 * Rewritten to use information that has already been read in from the
 * EXODUS II database, new input argument "exo" points to the database.
 *
 * Based on code from Randy Schunk's mm_flux
 *
 * Written by Richard Cairncross 4 September 1996
 *
 * Revised: 1997/08/19 08:55 MDT pasacki@sandia.gov
 *
 * Revised: 1999/01/28 14:57 MST pasacki@sandia.gov
 */

void
get_ss_element_list(const int ss_id1, /*  side set ID                   (in) */
		    const int ss_id2, /*  side set ID                   (in) */
		    const int ss_id3, /*  side set ID                   (in) */
		    int *num_elem, /*  length of element list          (out) */
		    int **elems, /*  pointer to list of elements       (out) */
		    int *ss_ptr, /*  ss index number of SS ID 1        (out) */
		    const int eq, /*  number of vector equation type    (in) */
		    const int irc, /* rotation condition index          (in) */
		    const Exo_DB *exo) /* the mesh                      (in) */

/*
 * Function to load a list of elements on a side-set or side-set intersections
 * into an array.
 */

{
  int common_elem_count;
  int i;
  int k;
  int inode;

  int *elem_list;
  int *edge_elem_list;

  int nodes_on_edge;

  int ss_1_index=-1;
  int ss_2_index=-1;
  int ss_3_index=-1;

  int ss_id1_found;
  int ss_id2_found;
  int ss_id3_found;

  int is_surface=FALSE;
  int is_edge=FALSE;
  int is_vertex=FALSE;

  int irc_node;			/* derived from ROT_list[candidate_node] */
  int *node_ssid_list;		/* derived from ROT_list[candidate_node] */

  int e;
  int elem;

  int index_node_in_ssid2;
  int index_node_in_ssid3;

  int len;
  int len2;
  int len3;

  int *list2;
  int *list3;

  int offset;

  int *vertex_elem_list;

  static const char yo[] = "get_ss_element_list"; /* this routine name */

  /*
   * Real oddball stuff that should not happen.
   */

  if ( ss_id1 < 1 )
    {
      log_err("Improper input SS ID 1 = %d", ss_id1);
    }

  if ( ss_id2 == 0 && ss_id3 != 0 )
    {
      log_err("Bad usage SS ID 2 = %d, SS ID 3 = %d", ss_id2, ss_id3);
    }

  if ( ss_id1 < 0 || ss_id2 < 0 || ss_id3 < 0 )
    {
      log_err("An SS ID is negative in %d, %d, %d", ss_id1, ss_id2, ss_id3);
    }

  is_surface = ( ss_id1 > 0 && ss_id2 == 0 && ss_id3 == 0 );
  is_edge    = ( ss_id1 > 0 && ss_id2 >  0 && ss_id3 == 0 );
  is_vertex  = ( ss_id1 > 0 && ss_id2 >  0 && ss_id3 >  0 );

  /*
   * Search to see if these side set ID's exist ON THIS PROCESSOR...
   */

  ss_1_index = in_list(ss_id1, 0, exo->num_side_sets, exo->ss_id);
  ss_2_index = in_list(ss_id2, 0, exo->num_side_sets, exo->ss_id);
  ss_3_index = in_list(ss_id3, 0, exo->num_side_sets, exo->ss_id);

  ss_id1_found = ( ss_1_index > -1 );
  ss_id2_found = ( ss_2_index > -1 );
  ss_id3_found = ( ss_3_index > -1 );

  /*
   * Set some default return values to answer some questions...
   *
   *	(1) There aren't any elements like you want.
   *    (2) They have no names.
   *    (3) They don't refer to any recognized side set index.
   */

  *num_elem = 0;
  *elems    = NULL;
  *ss_ptr   = -1;

  if ( is_surface )
    {
      if ( ! ss_id1_found )
	{
#ifdef PARALLEL
	  log_dbg("Did not find surface for SS ID %d for ROT %d\n",
		  ss_id1, irc);
	  return;
#endif
#ifndef PARALLEL
	  log_err("Expected to find SS ID %d for ROT %d", ss_id1, irc);
#endif
	}
    }

  /*
   * In parallel, this is possible, not having the entire mesh at hand...
   */

  if ( ! ss_id1_found )
    {
#ifndef PARALLEL      
      log_err("Not a surface, but did not find SS ID %d, ROT %d", ss_id1, irc);
#endif
#ifdef PARALLEL
      return;
#endif
    }

  len       = exo->ss_num_sides[ss_1_index];
  elem_list = (int *) smalloc(len*sizeof(int));

  offset    = exo->ss_elem_index[ss_1_index];
  for ( i=0; i<len; i++)
    {
      elem_list[i] = exo->ss_elem_list[offset+i];
    }

  if ( is_surface ) 
    {
      *ss_ptr   = ss_1_index;
      *num_elem = len;
      *elems    = elem_list;
      return;
    }

  /*
   * Edge. First, do some quick easy tests.
   */

  if ( is_edge && ( ! ss_id2_found || ! ss_id1_found ) )
    {
      *ss_ptr   = -1;
      *num_elem = 0;
      *elems    = NULL;
      safe_free(elem_list);
#ifdef PARALLEL
      log_dbg("Did not find edge for SS ID %d %d for ROT [%d]\n",
	      ss_id1, ss_id2, irc);
#endif
#ifndef PARALLEL
      //log_err("ROT [%d] wanted an edge, but either SS ID %d or %d gone!",
      //	      irc, ss_id1, ss_id2);
      log_dbg("Did not find edge for SS ID %d %d for ROT [%d]\n",
	      ss_id1, ss_id2, irc);
#endif
      return;	
    }

  /*
   * Edge, Second, examine ssid1 and ssid2 in detail for common nodes.
   *
   * In parallel it is possible for an edge not to be found even
   * if it has parts of both SS_ID 1 and SS_ID 2. Rigorously, sweep
   * through every element of SS_ID 1 and check for any nodes that
   * belong to SS_ID 2...
   */

  if ( is_edge && ( ss_id1_found && ss_id2_found ) )
    {

      /*
       * Create enough space to hold the edge element list...
       */

      len = MIN(exo->ss_num_sides[ss_1_index], exo->ss_num_sides[ss_2_index]);

      edge_elem_list = (int *) smalloc(len * sizeof(int));

      for ( i=0; i<len; i++)
	{
	  edge_elem_list[i] = -1;
	}

      common_elem_count = 0;

      /*
       * Set these up now for convenience in later searches...
       */

      len2  = exo->ss_node_side_index[ss_2_index]
	                             [exo->ss_num_sides[ss_2_index]];

      list2 = &(exo->ss_node_list[ss_2_index][0]);

      for ( e=0; e<exo->ss_num_sides[ss_1_index]; e++)
	{
	  elem          = exo->ss_elem_list[exo->ss_elem_index[ss_1_index]+e];

	  nodes_on_edge = 0;

	  for ( k=exo->ss_node_side_index[ss_1_index][e];
		k<exo->ss_node_side_index[ss_1_index][e+1]; k++ )
	    {
	      inode = exo->ss_node_list[ss_1_index][k];

	      /*
	       * This is not sufficient. Really, we need to check
	       * the nodes to see:
	       *
	       *        1) Do they have some ROT's associated with them?
	       *           (We do this below alread.)
	       *
	       *        2) Does the list of nodes for ssid2 ON THIS PROCESSOR
	       *           have this particular node?
	       */

	      index_node_in_ssid2 = in_list(inode, 0, len2, list2);

	      /*
	       * This means we have a potential edge. (It might be
	       * a pathological case with an element corner just touching
	       * the edge, but check later for nodes_on_edge>1 to eliminate
	       * that case.)
	       *
	       * However, this edge might not have any ROT condition 
	       * associated with it. In particular, any candidate nodes
	       * ought to have *this* ROT condition marked on it.
	       */

	      if ( index_node_in_ssid2 > -1 )
		{
		  if ( ROT_list[inode] != NULL )
		    {
		      /*
		       * If the node we found lists this irc, then
		       * the circle of trust is complete.
		       */
		      if ( ROT_list[inode][eq] == irc )
			{
			  nodes_on_edge++;
			}
		      /*
		       * But a VERTEX condition also qualifies and we
		       * have to check it because it stands first in line.
		       */

		      irc_node       = ROT_list[inode][eq];
		      node_ssid_list = ROT_Types[irc_node].ss_id;

		      if ( ROT_Types[irc_node].type == VERTEX &&
			   in_list(ss_id1, 0, 3, node_ssid_list) != -1 &&
			   in_list(ss_id2, 0, 3, node_ssid_list) != -1 )
			{
			  nodes_on_edge++; 
			}
		    }
		}
	    }

	  /* 
	   * Only count elements that have at least 2 nodes on the
	   * edge in 3D...
	   */

	  if ( nodes_on_edge > 1 )
	    {
	      edge_elem_list[common_elem_count] = elem;
	      common_elem_count++;
	    }
	}


      if ( common_elem_count > 0 )
	{
	  *num_elem = common_elem_count;
	  edge_elem_list = (int *) realloc(edge_elem_list, 
					   common_elem_count * sizeof(int));
	  *elems = edge_elem_list;
	  safe_free(elem_list);
	  return;
	}
      else
	{
#ifdef PARALLEL
	  log_dbg("ROT [%d] on SS ID %d & %d not found.", irc, ss_id1, ss_id2);
	  safe_free(elem_list);
	  return;
#endif
#ifndef PARALLEL
	  /*
	   * You had all your chances.
	   */
	  log_dbg("ROT [%d] on SS ID %d & %d not found.", irc, ss_id1, ss_id2);
	  safe_free(elem_list);
	  return;
	  //  log_err("ROT [%d] on SS ID %d & %d not found.", irc, ss_id1, ss_id2);
#endif
	}
    }

  /*
   * Vertex. First, do quick easy checks.
   */

  if ( is_vertex && ( ! ss_id2_found || ! ss_id1_found || ! ss_id3_found ) )
    {
      *ss_ptr   = -1;
      *num_elem = 0;
      *elems    = NULL;
#ifdef PARALLEL
      log_dbg("Did not find vertex for SS ID %d %d %d for ROT [%d]\n",
	      ss_id1, ss_id2, ss_id3, irc);
      safe_free(elem_list);
      return;	
#endif
#ifndef PARALLEL
      log_err("ROT [%d] wanted an edge, but one of SS ID %d %d %d gone!",
	      irc, ss_id1, ss_id2, ss_id3);
#endif
    }
    
  /*
   * Vertex. Look in detail.
   */

  if ( is_vertex )
    {

      len = MIN(exo->ss_num_sides[ss_1_index], exo->ss_num_sides[ss_2_index]);
      len = MIN(len, exo->ss_num_sides[ss_3_index]);

      vertex_elem_list = (int *) smalloc(len * sizeof(int));

      for ( i=0; i<len; i++)
	{
	  vertex_elem_list[i] = -1;
	}

      common_elem_count = 0;

      /*
       * Set these up now for convenience in later cross searches...
       */

      len2  = exo->ss_node_side_index[ss_2_index]
	                             [exo->ss_num_sides[ss_2_index]];

      list2 = &(exo->ss_node_list[ss_2_index][0]);

      len3  = exo->ss_node_side_index[ss_3_index]
	                             [exo->ss_num_sides[ss_3_index]];

      list3 = &(exo->ss_node_list[ss_3_index][0]);

      for ( e=0; e<exo->ss_num_sides[ss_1_index]; e++)
	{
	  elem          = exo->ss_elem_list[exo->ss_elem_index[ss_1_index]+e];

	  for ( k=exo->ss_node_side_index[ss_1_index][e];
		k<exo->ss_node_side_index[ss_1_index][e+1]; k++ )
	    {
	      inode = exo->ss_node_list[ss_1_index][k];

	      index_node_in_ssid2 = in_list(inode, 0, len2, list2);

	      index_node_in_ssid3 = in_list(inode, 0, len3, list3);

	      /*
	       * That node was found in both SS ID 2 and SS ID 3
	       */

	      if ( index_node_in_ssid2 > -1 &&
		   index_node_in_ssid3 > -1 )
		{
		  /*
		   * Maybe this geometric vertex doesn't have any conditions
		   * applied to it. Check and see...
		   */

		  if ( ROT_list[inode] != NULL )
		    {
		      /*
		       * Maybe there are multiple verteces and they even
		       * have ROT conditions, but not all of them have the
		       * same ROT conditions!
		       */

		      if ( ROT_list[inode][eq] == irc )
			{
			  /*
			   * I'm not sure why the following check is needed.
			   */

			  if ( inode == ROT_Types[irc].node )
			    {

			      /*
			       * But I do know that we only want to add
			       * such an element once. No telling whether
			       * some fool made a side set that wraps around
			       * the block several times...
			       */

			      if ( -1 == in_list(elem, 0, common_elem_count, 
						 vertex_elem_list) )
				{
				  vertex_elem_list[common_elem_count] = elem;
				  common_elem_count++;
				}
			    }
			}
		    }
		}
	    }
	}

      /*
       * What have we got after looking through all the side of SS ID 1?
       */

      if ( common_elem_count > 0 )
	{
	  *num_elem        = common_elem_count;
	  vertex_elem_list = (int *) realloc(vertex_elem_list, 
					     common_elem_count*sizeof(int));
	  *elems           = vertex_elem_list;
	  safe_free(elem_list);
	  return;
	}
      else
	{
#ifdef PARALLEL
	  log_dbg("ROT [%d] on SS ID %d & %d & %d not found.", irc, ss_id1, 
		  ss_id2, ss_id3);
	  safe_free(elem_list);
	  return;
#endif
#ifndef PARALLEL
	  /*
	   * You had all your chances.
	   */
	  log_err("ROT [%d] on SS ID %d & %d & %d not found.", irc, ss_id1, 
		  ss_id2, ss_id3);
#endif
	}
    }
      
  log_msg("A bad fall through.");
  log_err("Blame? ss_id1=%d, ss_id2=%d, ss_id3=%d, ROT [%d] on eqn %d",
	  ss_id1, ss_id2, ss_id3, irc, eq);

  return;
} /* END of get_ss_element_list */

/************************************************************************/
/************************************************************************/
/************************************************************************/




int
set_pointers_to_vectors (
    struct Rotation_Vectors **normal, 	/* defined in mm_as_structs.h */   
    struct Rotation_Vectors **tangent1, 	      
    struct Rotation_Vectors **tangent2, 	      
    struct Rotation_Vectors **normal2, 	      
    struct Rotation_Vectors **normal3, 	      
    struct Rotation_Vectors **line_tangent, 	      
    struct Rotation_Vectors **binormal, 	      
    struct Rotation_Vectors *dum_vect, 	      
    struct Rotation_Vectors *vector, 	      
    int I,                                  /* global node number                */
    int eq,                                 /* number of vector equation type    */
    int dim,                                /* dimensions                        */
    int irc )                               /* counter for rotation Specs        */

{
/* LOCAL VARIABLES */
  int j_id, p, q, b, num_rots_finished;

  /* Specify which types of vectors end up in which slots for 
   * rotation vectors and set up pointers so this becomes transparent */
  if (ROT_Types[irc].type == FACE) {
    /* defaults for SURFACE */
    *normal = dum_vect;
    *tangent1 = dum_vect + 1;
    *tangent2 = dum_vect + 2;
    *normal2 = NULL; *normal3 = NULL; *line_tangent = NULL; *binormal = NULL;

  } else if (ROT_Types[irc].type == CURVE) {
    /* defaults for EDGES */
    *normal = dum_vect;
    *line_tangent  = dum_vect + 1;
    *binormal  = dum_vect + 2;
    *normal2 = NULL; *normal3 = NULL; *tangent1 = NULL; *tangent2 = NULL;

  } else if (ROT_Types[irc].type == VERTEX) {
    /* defaults for VERTEX */
    *normal = dum_vect;
    *line_tangent  = dum_vect + 1;
    *binormal  = dum_vect + 2;
    *normal2 = NULL; *normal3 = NULL; *tangent1 = NULL; *tangent2 = NULL;

  } else EH(-1,"Illegal topology type");
  
  for (p=0; p<dim; p++) {
    dum_vect[p].ok = 0;
    dum_vect[p].d_vector_n = 0;
    vector[p].ok = 0;
    vector[p].d_vector_n = 0;
    for (j_id=0; j_id<MNROT; j_id++) {
      dum_vect[p].d_vector_J[j_id] = -1;
      vector[p].d_vector_J[j_id] = -1;
    }
    for (q=0; q<dim; q++) {
      dum_vect[p].vector[q] = 0.;
      vector[p].vector[q] = 0.;
      for (b=0; b<dim; b++) {
	for (j_id=0; j_id<MNROT; j_id++) {
	  dum_vect[p].d_vector_dx[q][b][j_id] = 0.;
	  vector[p].d_vector_dx[q][b][j_id] = 0.;
	}
      }
    }
  }
  
  num_rots_finished = 0;
  for (p=0; p<dim; p++) {
    if (ROT_Types[irc].BC_Type[p] > 0) {
/*     if (ROT_Types[irc].BC_Type[p] > 0 || rotation[I][eq][p]->ok) { */
      /* No rotation at this node */
      num_rots_finished++;
    } else {
      switch(ROT_Types[irc].BC_Type[p]) {
      case ROT_NONE:
	/* do nothing */
	num_rots_finished++;
	break;
	
      case ROT_N:
	*normal = vector + p;
	break;
	
      case ROT_N2:
	*normal2 = vector + p;
	EH(-1,"Second Normal Rotation not enabled yet");
	break;
	
      case ROT_N3:
	*normal3 = vector + p;
	EH(-1,"Third Normal Rotation not enabled yet");
	break;
	
      case ROT_T:
	*line_tangent = vector + p;
	if (ROT_Types[irc].type == FACE) 
	       EH(-1,"Need to specify T1 or T2 on Surfaces");
	break;
	
      case ROT_T1:
	*tangent1 = vector + p;
	break;
	
      case ROT_T2:
	*tangent2 = vector + p;
	break;
	
      case ROT_B:
	*binormal = vector + p;
	if (ROT_Types[irc].type == FACE) 
	       EH(-1,"Need to specify T1 or T2 on Surfaces");
	break;
	
      case ROT_S:
	rotation[I][eq][p]->vector[0] = ROT_Types[irc].seed[0];
	rotation[I][eq][p]->vector[1] = ROT_Types[irc].seed[1];
	rotation[I][eq][p]->vector[2] = ROT_Types[irc].seed[2];
	rotation[I][eq][p]->ok = 1;
	num_rots_finished++;
	break;
	
      case ROT_X:
	rotation[I][eq][p]->vector[0] = 1.;
	rotation[I][eq][p]->vector[1] = 0.;
	rotation[I][eq][p]->vector[2] = 0.;
	rotation[I][eq][p]->ok = 1;
	num_rots_finished++;
	break;
	
      case ROT_Y:
	rotation[I][eq][p]->vector[0] = 0.;
	rotation[I][eq][p]->vector[1] = 1.;
	rotation[I][eq][p]->vector[2] = 0.;
	rotation[I][eq][p]->ok = 1;
	num_rots_finished++;
	break;
	
      case ROT_Z:
	rotation[I][eq][p]->vector[0] = 0.;
	rotation[I][eq][p]->vector[1] = 0.;
	rotation[I][eq][p]->vector[2] = 1.;
	rotation[I][eq][p]->ok = 1;
	num_rots_finished++;
	break;
	
      default:
	EH(-1,"illegal rotation vector");
      } /* end of switch */
    }
  }	    
  
  return(num_rots_finished);
} /* END of set_pointers_to_vectors */


void
rotate_eqns_at_node_2D( int iconn,
			int dim,
			int num_local_nodes,
			struct Aztec_Linear_Solver_System *ams )  
{
  int i, id_mesh, id_mom,  I;
	
  for(i=0; i<num_local_nodes; i++)
    {
      I = Proc_Elem_Connect[iconn+i];
		
      if(rotation[I][VECT_EQ_MOM] != NULL && pd->e[pg->imtrx][R_MOMENTUM1] ) {
	id_mom = ei[pg->imtrx]->ln_to_first_dof[VELOCITY1][i];
			
	if(id_mom > -1 ) rotate_momentum_eqn(id_mom,  I, iconn, dim ,ams);
      }

      if(rotation[I][VECT_EQ_MESH] != NULL && pd->e[pg->imtrx][R_MESH1] ) {
	id_mesh = ei[pg->imtrx]->ln_to_first_dof[MESH_DISPLACEMENT1][i];
			
	if(id_mesh > -1 ) rotate_mesh_eqn(id_mesh,  I, iconn, dim ,ams);
      }
    }
  return;
}
						
/******************************************************************************************/

void rotate_momentum_auto(int id,  /* Elemental stiffness matrix row index */
                          int I,   /* Global node number                   */
                          int dim, /* physical dim of problem              */
                          struct Aztec_Linear_Solver_System *ams)

/*
 * Function which corrects the global residual vector "resid_vect" and
 * the global Jacobian vector "a" so that the vector momentum equations for surface nodes
 * are projected into a normal and tangential coordinate system.
 */

{
  /* LOCAL VARIABLES */
  int eq, pvar, peq;
  int n;
  int kdir, ldir;
  double rotated_resid[MDE];
  double rotated_jacobian_scalar[MAX_PDIM][MDE];

  /************************ EXECUTION BEGINS **********************************/
  eq = VECT_EQ_MOM;

  if (goma_automatic_rotations.rotation_nodes[I].n_normals == 0) {
    return; // not a rotated node
  }

  /* Correct residual equation first at local node "id" or global node "I" */
  /*                Rx -> Rn    and Ry -> Rt                                   */
  /*       i.e.,    Rn = nx*Rx + ny*Ry + nz*Rz                                 */
  /*                Rt1 = t1x*Rx + t1y*Ry + t1z*Rz                             */
  /*                Rt2 = t2x*Rx + t2y*Ry + t2z*Rz                             */

  double rc[DIM][DIM];
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      rc[i][j] = gds_vector_get(goma_automatic_rotations.rotation_nodes[I].rotated_coord[i], j);
    }
  }

  /* Now add on projection into n-t space */
  for (kdir = 0; kdir < dim; kdir++) {
    rotated_resid[kdir] = 0.;
    for (ldir = 0; ldir < dim; ldir++) {
      peq = upd->ep[pg->imtrx][R_MOMENTUM1 + ldir];
      rotated_resid[kdir] += rc[kdir][ldir] * lec->R[peq][id];
    }
  } /* end of loop over direction */

  for (kdir = 0; kdir < dim; kdir++) {
    peq = upd->ep[pg->imtrx][R_MOMENTUM1 + kdir];
    lec->R[peq][id] = rotated_resid[kdir];
  }

  /*                                                                      */
  /*   Now correct Jacobian                                               */
  if (af->Assemble_Jacobian) {
    for (int var = V_FIRST; var < V_LAST; var++) {
      if (pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];
        for (n = 0; n < ei[pg->imtrx]->dof[var]; n++) {

          rotated_jacobian_scalar[0][n] = 0.;
          rotated_jacobian_scalar[1][n] = 0.;
          rotated_jacobian_scalar[2][n] = 0.;

          for (kdir = 0; kdir < dim; kdir++) {
            for (ldir = 0; ldir < dim; ldir++) {
              rotated_jacobian_scalar[kdir][n] +=
                  rc[kdir][ldir] * lec->J[upd->ep[pg->imtrx][R_MOMENTUM1 + ldir]][pvar][id][n];
            }
          }

        } /* end of loop over nodes */

        /* reinject d_mesh/d_pressure back into lec-J for global assembly */
        for (kdir = 0; kdir < dim; kdir++) {
          /* loop over sensitivities */
          for (n = 0; n < ei[pg->imtrx]->dof[var]; n++) {

            lec->J[upd->ep[pg->imtrx][R_MOMENTUM1 + kdir]][pvar][id][n] =
                rotated_jacobian_scalar[kdir][n];
          }
        }

      } /* end of if variable */
    }
  }

  /* Put rotated residual back into lec for scattering into global matrix */
  /* Note this is the last thing we do so our chain rule still works for
     mesh derivatives wrt rotation vector */

} /* END of rotate_momentum_eqn */
