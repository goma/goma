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
 
/* wr_side_data -- print user-specified information at ea time step to file
 * 
 * Typically, this is for pulling nodal point data out into a convenient
 * "x,y" data file format.
 *
 * Created: 1998.12.12 16:39 MST prschun
 *
 */


/* Needed to declare POSIX function drand48 */
#define _XOPEN_SOURCE

#include <stdio.h>
#include <math.h>
#include <strings.h>

#include "std.h"
#include "rf_allo.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_io_const.h"
#include "rf_mp.h"
#include "el_geom.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_eh.h"
#include "mm_mp_structs.h"
#include "mm_post_def.h"
#include "mm_std_models_shell.h"
#include "dp_utils.h"
#include "el_elm.h"
#include "el_elm_info.h"
#include "exo_struct.h"
#include "md_timer.h"
#include "mm_fill_aux.h"
#include "mm_fill_fill.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_terms.h"
#include "mm_fill_util.h"
#include "mm_mp.h"
#include "mm_shell_util.h"
#include "mm_unknown_map.h"
#include "mpi.h"
#include "rd_mesh.h"
#include "wr_side_data.h"
 
extern double time_goma_started; /* def'd and set in main.c */

char err_msg[MAX_CHAR_IN_INPUT];

FILE	*uf;			/* file pointer for this user data */

int
ns_data_print(pp_Data * p, 
	      double x[], 
	      const Exo_DB * exo, 
	      const double time_value,
	      const double time_step_size)
{
  const int quantity       = p->data_type;
  int mat_num        = p->mat_num;
  const int elemBlock_id   = p->elem_blk_id;
  const int node_set_id    = p->ns_id;
  const int species_id     = p->species_number;
  const char * filenm      = p->data_filenm;
  const char * qtity_str   = p->data_type_name;
  const char * format_flag = p->format_flag;
  int * first_time         = &(p->first_time);

  static int err=0;
  int num_nodes_on_side;
  int ebIndex_first = -1;
  int local_side[2];
  int side_nodes[3];		/* Assume quad has no more than 3 per side. */
  int elem_list[4], elem_ct=0, face, ielem, node2;
  int local_node[4];
  int node = -1;
  int idx, idy, idz, id_var;
  int iprint;
  int nsp;			/* node set pointer for this node set */
  dbl x_pos, y_pos, z_pos;
  int j, wspec;
  int doPressure = 0;

#ifdef PARALLEL
  double some_time=0.0;
#endif
  double abscissa=0;
  double ordinate=0;
  double n1[3], n2[3];
  double xi[3];

  /*
   * Find an element block that has the desired material id.
   */
  if (elemBlock_id != -1) {
    for (j = 0; j < exo->num_elem_blocks; j++) {
      if (elemBlock_id == exo->eb_id[j]) {
	ebIndex_first = j;
	break;
      }
    }
    if (ebIndex_first == -1) {
      sprintf(err_msg, "Can't find an element block with the elem Block id %d\n", elemBlock_id);
    if (Num_Proc == 1) {
      EH(GOMA_ERROR, err_msg);
    }
    }
    mat_num = Matilda[ebIndex_first];
    p->mat_num = mat_num;
    pd = pd_glob[mat_num];
  } else {
    mat_num = -1;
    p->mat_num = -1;
    pd = pd_glob[0];
  }

  nsp = match_nsid(node_set_id);  

  if( nsp != -1 )
    {
      node = Proc_NS_List[Proc_NS_Pointers[nsp]];
    }
  else
    {
      sprintf(err_msg, "Node set ID %d not found.", node_set_id);
      if( Num_Proc == 1 ) EH(GOMA_ERROR,err_msg);
    }

  /* first right time stamp or run stamp to separate the sets */

  print_sync_start(FALSE);

  if (*first_time)
    {
      if ( format_flag[0] != '\0' ) 
	{
	  if (ProcID == 0)
	    {
	      uf = fopen(filenm,"a");
	      if (uf != NULL)
		{
		  fprintf(uf,"# %s %s @ nsid %d node (%d) \n", 
			  format_flag, qtity_str, node_set_id, node );
		  *first_time = FALSE;
		  fclose(uf);
		}
	    }
	}
    }

  if (format_flag[0] == '\0')
    {
      if (ProcID == 0)
	{
	  if ((uf = fopen(filenm,"a")) != NULL)
	    {
	      fprintf(uf,"Time/iteration = %e \n", time_value);
	      fprintf(uf,"  %s   Node_Set  %d Species %d\n", qtity_str,node_set_id,species_id);
	      fflush(uf);
	      fclose(uf);
	    }
	}
    }

  if (nsp != -1 ) {

    for (j = 0; j < Proc_NS_Count[nsp]; j++) {
      node = Proc_NS_List[Proc_NS_Pointers[nsp]+j];
      if (node < num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx] ) {
        idx = Index_Solution(node, MESH_DISPLACEMENT1, 0, 0, -1, pg->imtrx);
        if (idx == -1) {
          x_pos = Coor[0][node];
          WH(idx, "Mesh variable not found.  May get undeformed coords.");
        } else {
          x_pos = Coor[0][node] + x[idx];
        }
        idy = Index_Solution(node, MESH_DISPLACEMENT2, 0, 0, -1, pg->imtrx);
        if (idy == -1) {
          y_pos = Coor[1][node];
        } else {
          y_pos = Coor[1][node] + x[idy];
        }
        z_pos = 0.;
        if(pd->Num_Dim == 3) {
          idz = Index_Solution(node, MESH_DISPLACEMENT3, 0, 0, -1, pg->imtrx);
          if (idz == -1) {
	    z_pos = Coor[2][node];
          }  else{
	    z_pos = Coor[2][node] + x[idz];
          }
        }
        if (quantity == MASS_FRACTION) {
          id_var = Index_Solution(node, quantity, species_id, 0, mat_num, pg->imtrx);
        } else if (quantity < 0) {
          id_var = -1;
        } else {
          id_var = Index_Solution(node, quantity, 0, 0, mat_num, pg->imtrx);
        }

	/*
	 * In the easy case, the variable can be found somewhere in the
	 * big vector of unknowns. But sometimes we want a derived quantity
	 * that requires a little more work than an array reference.
	 *
	 * For now, save the good result if we have it.
	 */

	if ( id_var != -1 )
	  {
	    ordinate = x[id_var];
	    iprint = 1;
	  }
	else
	  {
	    /*
	     *  If we have an element based interpolation, let's calculate the interpolated value
	     */
	    if (quantity == PRESSURE) {
	      if ((pd->i[pg->imtrx][PRESSURE] == I_P1) || ( (pd->i[pg->imtrx][PRESSURE] > I_PQ1) && (pd->i[pg->imtrx][PRESSURE] < I_Q2_HVG) )) {
		doPressure = 1;
	      }
	    }
	    iprint = 0;
	  }

	/*
	 * If the quantity is "theta", an interior angle that only
	 * makes sense at a point, in 2D, we'll need to compute it.
	 */

	if ( strncasecmp(qtity_str, "theta", 5 ) == 0 || doPressure)
	  {
	    /*
	     * Look for the two sides connected to this node...?
	     *
	     * Premise:
	     *	1. The node appears in only one element(removed-RBS,6/14/06)
	     *          2. Exactly two sides emanate from the node.
	     *          3. Quadrilateral.
	     *
	     * Apologies to people who wish to relax premise 1. I know
	     * there are some obtuse angles out there that benefit from
	     * having more than one element at a vertex. With care, this
	     * procedure could be extended to cover that case as well.
	     */
	    
	    if ( ! exo->node_elem_conn_exists )
	      {
		EH(GOMA_ERROR, "Cannot compute angle without node_elem_conn.");
	      }
	    
	    elem_list[0] = exo->node_elem_list[exo->node_elem_pntr[node]];

	    /*
	     * Find out where this node appears in the elements local
	     * node ordering scheme...
	     */

	    local_node[0] = in_list(node, exo->elem_node_pntr[elem_list[0]], 
				    exo->elem_node_pntr[elem_list[0]+1],
				    exo->elem_node_list);

	    EH(local_node[0], "Can not find node in elem node connectivity!?! ");
	    local_node[0] -= exo->elem_node_pntr[elem_list[0]];
	    /* check for neighbors*/
	    if( mat_num == find_mat_number(elem_list[0], exo))
	      {elem_ct = 1;}
	    else
	      {WH(-1,"block id doesn't match first element");}
	    for (face=0 ; face<ei[pg->imtrx]->num_sides ; face++)
	      {
		ielem = exo->elem_elem_list[exo->elem_elem_pntr[elem_list[0]]+face];
		if (ielem != -1)
		  {
		    node2 = in_list(node, exo->elem_node_pntr[ielem], 
				    exo->elem_node_pntr[ielem+1],
				    exo->elem_node_list);
		    if (node2 != -1 && (mat_num == find_mat_number(ielem, exo)))
		      {
			elem_list[elem_ct] = ielem;
			local_node[elem_ct] = node2;
			local_node[elem_ct] -= exo->elem_node_pntr[ielem];
			elem_ct++;
		      }
		  }
	      }

	    /*
	     * Note that indeces are zero based!
	     */

	    ordinate = 0.0;
	    for (ielem = 0 ; ielem < elem_ct ; ielem++)
	      {
		if ( local_node[ielem] < 0 || local_node[ielem] > 3 ) 
		  {
		    if (strncasecmp(qtity_str, "theta", 5 ) == 0) {
		      EH(GOMA_ERROR, "Node out of bounds.");
		    }
		  }

		/*
		 * Now, determine the local name of the sides adjacent to this
		 * node...this works for the exo patran convention for quads...
		 *
		 * Again, local_node and local_side are zero based...
		 */

		local_side[0] = (local_node[ielem]+3)%4;
		local_side[1] = local_node[ielem];

		/*
		 * With the side names, we can find the normal vector.
		 * Again, assume the sides live on the same element.
		 */
		load_ei(elem_list[ielem], exo, 0, pg->imtrx);

		/*
		 * We abuse the argument list under the conditions that
		 * we're going to do read-only operations and that
		 * we're not interested in old time steps, time derivatives
		 * etc.
		 */
		if (x == x_static) /* be the least disruptive possible */
		  {
		    err = load_elem_dofptr(elem_list[ielem], exo, x_static, x_old_static,
                                           xdot_static, xdot_old_static, 1);
		  }
		else
		  {
                    err = load_elem_dofptr(elem_list[ielem], exo, x, x, x, x, 1);
		  }

		/*
		 * What are the local coordinates of the nodes in a quadrilateral?
		 */

		find_nodal_stu(local_node[ielem], ei[pg->imtrx]->ielem_type, xi, xi+1, xi+2);

		err = load_basis_functions(xi, bfd);

		EH( err, "problem from load_basis_functions");
	    
		err = beer_belly();
		EH( err, "beer_belly");
	    
		err = load_fv();
		EH( err, "load_fv");
	    
		err = load_bf_grad();
		EH( err, "load_bf_grad");

		err = load_bf_mesh_derivs(); 
		EH(err, "load_bf_mesh_derivs");
		
		if (doPressure) {
		  ordinate = fv->P;
		  iprint = 1;
		} else {

		/* First, one side... */

		get_side_info(ei[pg->imtrx]->ielem_type, local_side[0]+1, &num_nodes_on_side, 
			      side_nodes);

		surface_determinant_and_normal(elem_list[ielem], 
					       exo->elem_node_pntr[elem_list[ielem]],
					       ei[pg->imtrx]->num_local_nodes,  
					       ei[pg->imtrx]->ielem_dim-1, 
					       local_side[0]+1, 
					       num_nodes_on_side,
					       side_nodes);

		n1[0] = fv->snormal[0];
		n1[1] = fv->snormal[1];

		/* Second, the adjacent side of the quad... */

		get_side_info(ei[pg->imtrx]->ielem_type, local_side[1]+1, &num_nodes_on_side, 
			      side_nodes);

		surface_determinant_and_normal(elem_list[ielem], 
					       exo->elem_node_pntr[elem_list[ielem]],
					       ei[pg->imtrx]->num_local_nodes,  
					       ei[pg->imtrx]->ielem_dim-1, 
					       local_side[1]+1, 
					       num_nodes_on_side,
					       side_nodes);

		n2[0] = fv->snormal[0];
		n2[1] = fv->snormal[1];

		/* cos (theta) = n1.n2 / ||n1|| ||n2|| */

		ordinate += 180. - (180./M_PI)*acos((n1[0]*n2[0] + n1[1]*n2[1])/
						    (sqrt(n1[0]*n1[0]+n1[1]*n1[1])*
						     sqrt(n2[0]*n2[0]+n2[1]*n2[1])));
		}
		iprint = 1;
	      }	/*ielem loop	*/
	  }
	else if ( strncasecmp(qtity_str, "timestepsize", 12 ) == 0 )
	  {
	    ordinate = time_step_size;
	    iprint = 1;
	  }
	else if ( strncasecmp(qtity_str, "cputime", 7 ) == 0 )
	  {
	    ordinate = ut();
	    iprint = 1;
	  }
	else if ( strncasecmp(qtity_str, "wallclocktime", 13 ) == 0 )
	  {
	    /* Get these from extern via main...*/
#ifdef PARALLEL
	    some_time = MPI_Wtime();
	    ordinate = some_time - time_goma_started;
#endif
#ifndef PARALLEL
            time_t now=0;
	    (void)time(&now);
	    ordinate = (double)(now) - time_goma_started;
#endif
	    iprint = 1;
	  }
	else if ( strncasecmp(qtity_str, "speed", 5 ) == 0 )
	  {
            id_var = Index_Solution(node, VELOCITY1, 0, 0, mat_num, pg->imtrx);
	    ordinate = SQUARE(x[id_var]);
            id_var = Index_Solution(node, VELOCITY2, 0, 0, mat_num, pg->imtrx);
	    ordinate += SQUARE(x[id_var]);
            id_var = Index_Solution(node, VELOCITY3, 0, 0, mat_num, pg->imtrx);
	    ordinate += SQUARE(x[id_var]);
	    ordinate = sqrt(ordinate);
	    iprint = 1;
	  }
        else if ( strncasecmp(qtity_str, "ac_pres", 7 ) == 0 )
          {
            id_var = Index_Solution(node, ACOUS_PREAL, 0, 0, mat_num, pg->imtrx);
            ordinate = SQUARE(x[id_var]);
            id_var = Index_Solution(node, ACOUS_PIMAG, 0, 0, mat_num, pg->imtrx);
            ordinate += SQUARE(x[id_var]);
            ordinate = sqrt(ordinate);
            iprint = 1;
          }
        else if ( strncasecmp(qtity_str, "light_comp", 10 ) == 0 )
          {
            id_var = Index_Solution(node, LIGHT_INTP, 0, 0, mat_num, pg->imtrx);
            ordinate = x[id_var];
            id_var = Index_Solution(node, LIGHT_INTM, 0, 0, mat_num, pg->imtrx);
            ordinate += x[id_var];
            iprint = 1;
          }
        else if ( strncasecmp(qtity_str, "lub_height", 10 ) == 0 )
          {
             dbl H, H_U, dH_U_dtime, H_L, dH_L_dtime;
             dbl dH_U_dX[DIM],dH_L_dX[DIM], dH_U_dp, dH_U_ddh;
             int *n_dof = NULL;
             int dof_map[MDE];
	    if ( ! exo->node_elem_conn_exists )
	      {
		EH(GOMA_ERROR, "Cannot compute angle without node_elem_conn.");
	      }
	    
	    elem_list[0] = exo->node_elem_list[exo->node_elem_pntr[node]];

	    /*
	     * Find out where this node appears in the elements local
	     * node ordering scheme...
	     */

	    local_node[0] = in_list(node, exo->elem_node_pntr[elem_list[0]], 
				    exo->elem_node_pntr[elem_list[0]+1],
				    exo->elem_node_list);

	    EH(local_node[0], "Can not find node in elem node connectivity!?! ");
	    local_node[0] -= exo->elem_node_pntr[elem_list[0]];

		load_ei(elem_list[0], exo, 0, pg->imtrx);

		find_nodal_stu(local_node[0], ei[pg->imtrx]->ielem_type, xi, xi+1, xi+2);
		err = load_basis_functions(xi, bfd);

		EH( err, "problem from load_basis_functions");
	    
		err = beer_belly();
		EH( err, "beer_belly");
	    
		err = load_fv();
		EH( err, "load_fv");
	    
            n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
            lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);
                                                                 

             /* Lubrication height from model */
                H = height_function_model(&H_U, &dH_U_dtime, &H_L, &dH_L_dtime, 
                     dH_U_dX, dH_L_dX, &dH_U_dp, &dH_U_ddh, time_value, time_step_size); 
            ordinate = H;
            iprint = 1;
          }
        else if ( strncasecmp(qtity_str, "external_field", 14 ) == 0 )
          {
            id_var = Index_Solution(node, MASS_FRACTION, species_id, 0, mat_num, pg->imtrx);
            ordinate = efv->ext_fld_ndl_val[species_id][node];
            iprint = 1;
          }
        else if ( strncasecmp(qtity_str, "untracked", 11 ) == 0 )
          {
            double density_tot=0.;
            switch(mp_glob[mat_num]->Species_Var_Type)   {
               case SPECIES_CONCENTRATION:
                    density_tot = calc_density(mp_glob[mat_num], FALSE, NULL, 0.0);
                    ordinate = density_tot;
	            for(wspec = 0 ; wspec < pd->Num_Species_Eqn ; wspec++)
		        {
            	          id_var = Index_Solution(node, MASS_FRACTION, wspec, 0, mat_num, pg->imtrx);
            	          ordinate -= x[id_var]*mp_glob[mat_num]->molecular_weight[wspec];
		        }
            	     ordinate /= mp_glob[mat_num]->molecular_weight[pd->Num_Species_Eqn];
               break;
               case SPECIES_DENSITY:
                    density_tot = calc_density(mp_glob[mat_num], FALSE, NULL, 0.0);
                    ordinate = density_tot;
	            for(wspec = 0 ; wspec < pd->Num_Species_Eqn ; wspec++)
		        {
            	          id_var = Index_Solution(node, MASS_FRACTION, wspec, 0, mat_num, pg->imtrx);
            	          ordinate -= x[id_var];
		        }
               break;
               case SPECIES_MASS_FRACTION:
               case SPECIES_UNDEFINED_FORM:
                    ordinate = 1.0;
	            for(wspec = 0 ; wspec < pd->Num_Species_Eqn ; wspec++)
		        {
            	          id_var = Index_Solution(node, MASS_FRACTION, wspec, 0, mat_num, pg->imtrx);
            	          ordinate -= x[id_var];
		        }
               break;
               default:
                    WH(-1,"Undefined Species Type in untracked species\n");
               }
            iprint = 1;
          }
	else
	  {
	    WH(id_var,
	       "Requested print variable is not defined at all nodes. May get 0.");
	    if(id_var == -1) iprint = 0;
	  }

        if ((uf=fopen(filenm,"a")) != NULL)
	  {
	    if ( format_flag[0] == '\0' )
	      {
		if (iprint)
		  {
		    fprintf(uf,"  %e %e %e %e \n", x_pos, y_pos, z_pos, ordinate);
		  }
	      }
	    else
	      {
		if ( strncasecmp(format_flag, "t", 1) == 0 )
		  {
		    abscissa = time_value;
		  }
		else if (  strncasecmp(format_flag, "x", 1) == 0 )
		  {
		    abscissa = x_pos;		      
		  }
		else if (  strncasecmp(format_flag, "y", 1) == 0 )
		  {
		    abscissa = y_pos;		      
		  }
		else if (  strncasecmp(format_flag, "z", 1) == 0 )
		  {
		    abscissa = z_pos;		      
		  }
		else
		  {
		    abscissa = 0;
		  }
		if (iprint)
		  {
		    fprintf(uf, "%.16g\t%.16g\n", abscissa, ordinate);
		  }
	      }
	    fclose(uf);
	  }
      }
    }
  }
  print_sync_end(FALSE);

  return(1);
} /* END of routine ns_data_print */
/******************************************************************************/

int
ns_data_sens_print(const struct Post_Processing_Data_Sens *p,
		   const double x[], /* solution vector */
		   double **x_sens,  /* solution sensitivity vector */
		   const double time_value) /* current time */
{
  const int node_set_id = p->ns_id;
  const int quantity    = p->data_type;
  const int mat_id      = p->mat_id;
  const int species_id  = p->species_number;
  const int sens_type   = p->sens_type;
  const int sens_id     = p->sens_id;
  const int sens_flt    = p->sens_flt;
  const int sens_flt2   = p->sens_flt2;
  const char *filenm    = p->data_filenm;
  const char *qtity_str = p->data_type_name;
  const int sens_ct     = p->vector_id;

  int node;
  int idx, idy, idz, id_var;
  int nsp;                      /* node set pointer for this node set */
  dbl x_pos, y_pos, z_pos;
  int j;

  nsp            = match_nsid(node_set_id);
  if( nsp != -1 ) {
     node           = Proc_NS_List[Proc_NS_Pointers[nsp]];
  }
  else
  {
    sprintf(err_msg, "Node set ID %d not found.", node_set_id);
    if( Num_Proc == 1 ) EH(GOMA_ERROR,err_msg);
  }

  /* first right time stamp or run stamp to separate the sets */

  print_sync_start(TRUE);

  if (ProcID == 0 && (uf=fopen(filenm,"a")) != NULL)
    {
      fprintf(uf,"Time/iteration = %e \n", time_value);
      fprintf(uf,"  %s   Node_Set  %d \n", qtity_str,node_set_id);
	if(sens_type == 1)
	{
      fprintf(uf,"Sensitivity_type BC  ID  %d  Float  %d\n",sens_id,sens_flt);
	}
	else
	if(sens_type == 2)
	{
      fprintf(uf,"Sensitivity_type MT  NO  %d  Prop.  %d\n",sens_id,sens_flt);
	}
	else
	if(sens_type == 3)
	{
      fprintf(uf,"Sensitivity_type AC  ID  %d  Float  %d\n",sens_id,sens_flt);
	}
	else
	if(sens_type == 4)
	{
      fprintf(uf,"Sensitivity_type UM  NO  %d  Prop.  %d %d\n",sens_id,sens_flt,sens_flt2);
	}
	else
	if(sens_type == 5)
	{
      fprintf(uf,"Sensitivity_type UF  ID  %d  Float  %d\n",sens_id,sens_flt);
	}
	else
	if(sens_type == 6)
	{
      fprintf(uf,"Sensitivity_type AN  ID  %d  Float  %d\n",sens_id,sens_flt);
	}
      fflush(uf);
      fclose(uf);
    }

  if( nsp != -1 ) {

    for (j = 0; j < Proc_NS_Count[nsp]; j++) {
      node = Proc_NS_List[Proc_NS_Pointers[nsp]+j];
      if( node < num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx] ) {
        idx = Index_Solution (node, MESH_DISPLACEMENT1, 0, 0, -1, pg->imtrx);
        if (idx == -1) {
	  x_pos = Coor[0][node];
          WH(idx, "Mesh variable not found.  May get undeformed coords.");
        } else {
	  x_pos = Coor[0][node] + x[idx];
        }
        idy = Index_Solution (node, MESH_DISPLACEMENT2, 0, 0, -1, pg->imtrx);
        if (idy == -1) {
	  y_pos = Coor[1][node];
        } else {
	  y_pos = Coor[1][node] + x[idy];
        }
        z_pos = 0.;
        if (pd->Num_Dim == 3) {
	  idz = Index_Solution(node, MESH_DISPLACEMENT3, 0, 0, -1, pg->imtrx);
	  if (idz == -1) {
	    z_pos = Coor[2][node];
	  } else {
	    z_pos = Coor[2][node] + x[idz];
	  }
        }
        if(quantity == MASS_FRACTION) {
	  id_var = Index_Solution(node, quantity, species_id, 0, mat_id, pg->imtrx);
        } else {
	  id_var = Index_Solution(node, quantity, 0, 0, mat_id, pg->imtrx);
        }
        WH(id_var,
	   "Requested print variable is not defined at all nodes. May get 0.");

        if ((uf=fopen(filenm,"a")) != NULL) {
	  if (id_var != -1) {
	    fprintf(uf, "  %e %e %e %e \n",
		    x_pos, y_pos, z_pos, x_sens[sens_ct][id_var]);
	  }
	  fclose(uf);
        }
      }
    }
  }

  print_sync_end(TRUE);

  return(1);
} /* END of routine ns_data_sens_print */
/*****************************************************************************/

/*
 * match_nsid() -- find index of nodesets possessing given NS_ID
 */

int
match_nsid ( int id )

{
  int i;
  
  for ( i=0; i<Proc_Num_Node_Sets; i++)
    {
      if (Proc_NS_Ids[i] == id )
	{
	  return(i);
	}
    }
  return(-1);
} /* END of routine match_nsid */
/******************************************************************************/

/*
 * psid2nn() -- find global node number given a point node set id, 
 */

int
psid2nn ( int psid )

{
  int nsp;			/* node set pointer for this pointset */

  nsp     = match_nsid(psid);	/* node set pointer this node set */

  if( nsp == -1 )
    {
      if( Num_Proc == 1 )
       {
        EH(nsp, "Failed to match nodeset.");
       }
      else
        return(-1);
    }

  return( Proc_NS_List[ Proc_NS_Pointers[nsp] ] ); /* global node num */

}/* END of routine psid2nn */
/******************************************************************************/


/*
 * nsid2nn() -- find global node numbers given a generalized node set id, 
 */

int
nsid2nn ( int psid )

{
  int nsp;			/* node set pointer for this pointset */

  nsp     = match_nsid(psid);	/* node set pointer this node set */

  if( nsp == -1 )
    {
      if( Num_Proc == 1 )
       {
        EH(nsp, "Failed to match nodeset.");
       }
      else
        return(-1);
    }

  return( Proc_NS_List[ Proc_NS_Pointers[nsp] ] ); /* global node num */

}/* END of routine nsid2nn */
/******************************************************************************/
/* END of file wr_side_data */
/******************************************************************************/
