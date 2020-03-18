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
 

/* Standard include files */
 
#include <string.h>
#include <math.h>
 
/* GOMA include files */

#include "std.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "el_elm.h"
#include "el_geom.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_eh.h"
#include "el_elm_info.h"
#include "exo_struct.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_util.h"
#include "mm_post_proc.h"
#include "mm_unknown_map.h"
#include "rd_mesh.h"

#define GOMA_MM_FLUX_C
#define GOMA_MM_POST_PROC_UTIL_C



/* find_id_elem() -- Find global element id given the global coordinates 
 *                   of a point.   
 *
 * Author:          P. R. Schunk
 * Date:            8 Nov. 1999
 *
 * (o) Routine to search and find the global element id containing the
 *     coordinates of the point that is input 
 *
 */

int
find_id_elem(const double x,	  /* x-coordinate */
	     const double y,      /* y-coordinate */
	     const double z,      /* z-coordinate */
             double  xv[],
             const Exo_DB * exo,  /* ptr to basic exodus ii mesh information */
	     const int e_start,   /* start of search */  
	     const int e_end)     /* end of search */  
{
#ifdef DEBUG
  static char *yo = "find_id_elem";
#endif

  int element_no = -1;		/*local element number */
  int i;			/*local element counters */
  dbl sum_xcoord, sum_ycoord, sum_zcoord; 
  dbl x_center, y_center, z_center; 
  int id, I;           /*global node number */
  int i_elem;
  int mode = 0;        /* temporary.  You may want to pass in someday */
  dbl separation, closest_distance;
  int iconn_ptr, ebn, mn, ielem_type, num_local_nodes;


if(!mode)  
  {
    element_no = -1; 
    closest_distance = 1.e30;
    
    /*    for (i_elem = 0; i_elem < exo->num_elems; i_elem++ ) */
    /*PRS: Need to fix this up to generalize this for blocks 
      vs. global 0 to numelem*/
    for (i_elem = e_start; i_elem < e_end; i_elem++ )
      {
	/*Compute coordinates of "element center of mass" */
	/* TOO EXPENSIVE TO load_ei
	err = load_ei( i_elem, exo);
	EH( err,"load_ei" );
        */
        ebn = find_elemblock_index(i_elem, exo);
        mn = Matilda[ebn];
        iconn_ptr = exo->elem_ptr[i_elem];
        ielem_type = Elem_Type(exo, i_elem);
        num_local_nodes = elem_info(NNODES, ielem_type);
	
	sum_xcoord = sum_ycoord = sum_zcoord = 0.;
	for(id=0; id<num_local_nodes; id++)
	  {
	    I = exo->node_list[iconn_ptr + id];
            i = Index_Solution(I, MESH_DISPLACEMENT1, 0, 0, -1, pg->imtrx);
	    if (i == -1 && pd_glob[mn]->IntegrationMap != SUBPARAMETRIC)
		{
            	    sum_xcoord += Coor[0][I];
       		    sum_ycoord += Coor[1][I];
	            if(exo->num_dim > 2) sum_zcoord += Coor[2][I];
		}
	    if (i != -1 )
		{
            	    sum_xcoord += Coor[0][I] + xv[i];
                    sum_ycoord += Coor[1][I] + xv[i+1];
                    if(exo->num_dim > 2) sum_zcoord += Coor[2][I] + xv[i+2];
		}
	  }
	x_center = sum_xcoord/num_local_nodes;
	y_center = sum_ycoord/num_local_nodes;
	z_center = sum_zcoord/num_local_nodes;

	separation = sqrt((x-x_center)*(x-x_center)+
			  (y-y_center)*(y-y_center)+
			  (z-z_center)*(z-z_center));
	if(separation <= closest_distance) 
	  {
	    element_no = i_elem;
	    closest_distance = separation;
	  }
      }
  }
else if (mode) /*do a local search using exo->node_elem */
  {
    EH(-1, " Local node-elem list search not available for find_id_elem");
  }
 
else
  {
    EH(-1,"problem in find_elem_id: need current element number");
  }
  return(element_no);		/* failsafe default? */
}


/*
 *  invert_isoparametric_map() - Used to search element number and compute local
 *                               isoparametric coordinates given the global 
 *                               coordinates as input 
 *
 *    Author:          P. R. Schunk
 *    Date:            8 Nov. 1999
 */


int
invert_isoparametric_map(  int *current_ielem,  /* initial element of search */
                           const double coordinate[],
                           double xi[],          /* s-coordinate output */
                           const Exo_DB * exo,
                           double xv[], 	 /*  solution vector  */
			   int *velo_interp)	/* velocity basis fcns  */
{
  /* Local variables */
  int err;
  int a, b, c; 			/*  loop counters */
  int dim;			/*  problem dimension */
  int converged, inewton;  /*  convergence flag, iteration counter */
  double R[MAX_PDIM];		/* residual vector for invert */
  double update[MAX_PDIM] = {0};	/* update vector for xi */
  double norm;			/* convergence norm */
  double max_xi, tmp;		/*  element switching vars, flags  */
  int i_max_xi, face=0, newface=0, old_ielem ;
  
  double xi_tmp[MAX_PDIM];	/*  extra vector for isopar. coords */
  double Jinv[MAX_PDIM][MAX_PDIM];    	/*  store old element Jac inverse */
  double rot[MAX_PDIM][MAX_PDIM];	/*  isopar. coord rotation tensor */

  
  dim    = ei[pg->imtrx]->ielem_dim;

  if (xv == x_static) /* be the least disruptive possible */
    {
      err = load_elem_dofptr(*current_ielem, exo, x_static, x_old_static,
                             xdot_static, xdot_old_static, 0);
    }
  else
    {
      err = load_elem_dofptr(*current_ielem, exo, xv, xv, xv, xv, 0);
    }

  /* Initialize */

  converged = 0;
  inewton = 0;
 
      

  while ( ( !converged ) && ( inewton < 50 ) )
    {
      err = load_basis_functions( xi, bfd );
      EH( err, "problem from load_basis_functions");

      /*
       * This has elemental Jacobian transformation and some 
       * basic mesh derivatives...
       *
       * The physical space derivatives from beer_belly are still
       * "raw", in the sense that they do not yet include the 
       * scale factors necessary to make them into *gradients*.
       * That is done in load_fv.
       */
      
      err = beer_belly();
      EH( err, "beer_belly");
      if( neg_elem_volume ) return -1;

      /* form residual equations */

      for (a = 0 ; a < dim; a++)
	{
          R[a] = coordinate[a] - fv->x[a];
	}

      /* Solve matrix system with inverted Jacobian from beer_belly */

        for(a = 0 ; a < dim ; a++){update[a]=0.;}

      for (a = 0 ; a < dim; a++)
	{
	  for (b = 0 ; b < dim; b++)
	    {
	      
              update[a] += bfd[*velo_interp]->B[b][a]*R[b];
	    }
	}
      
      for (a = 0 ; a < dim; a++) xi[a] += update[a];

      /*Element switch if local coordinates xi fall outside
	current element */

        max_xi = 0;
        i_max_xi = 0;
        for(a = 0 ; a < dim ; a++)
                {
                if(fabs(xi[a]) > max_xi){max_xi = fabs(xi[a]);i_max_xi=a;}
                }


        if(max_xi > 1.)

                /*  we have left the current element  */
        {
                if(xi[i_max_xi] < -1.)
                {
                        switch (i_max_xi)
                        {
                        case 0:
                                face = 3; break;
                        case 1:
                                face = 0; break;
                        case 2:
                                face = 4; break;
                        }
                }
                else
                {
                        switch (i_max_xi)
                        {
                        case 0:
                                face = 1; break;
                        case 1:
                                face = 2; break;
                        case 2:
                                face = 5; break;
                        }
                }

/*  before we switch elements, store Jac inverse */
	/*  find centroid of leaving face */
	
	for(b=0;b<MAX_PDIM;b++){xi_tmp[b]=0;}
	switch (face)
		{
		case 0:
			xi_tmp[1]=-1.; break;
		case 1:
			xi_tmp[0]=1.; break;
		case 2:
			xi_tmp[1]=1.; break;
		case 3:
			xi_tmp[0]=-1.; break;
		case 4:
			xi_tmp[2]=-1.; break;
		case 5:
			xi_tmp[2]=1.; break;
		}
      err = load_basis_functions( xi_tmp, bfd );
      EH( err, "problem from load_basis_functions");

      err = beer_belly();
      EH( err, "beer_belly");
      
	for(a=0;a<dim;a++)
		{
		for(b=0;b<dim;b++){Jinv[a][b]= bfd[*velo_interp]->B[a][b];}
		}

/*   switch to the next element    */

       a =  exo->elem_elem_pntr[*current_ielem] + face;

/*   test to see if there is an adjacent element  */

        if(exo->elem_elem_list[a] != -1)
           {
                old_ielem = *current_ielem;
                *current_ielem = exo->elem_elem_list[a];
        for (a = 0 ; a < 2*dim ; a++)
                {
                b = exo->elem_elem_pntr[*current_ielem] + a;
                if(exo->elem_elem_list[b] == old_ielem)newface = a;
                }

        load_ei(*current_ielem, exo, 0, pg->imtrx);

        if (xv == x_static) /* be the least disruptive possible */
          {
            err = load_elem_dofptr(*current_ielem, exo, x_static, x_old_static,
                                   xdot_static, xdot_old_static, 0);
          }
        else
          {
            err = load_elem_dofptr(*current_ielem, exo, xv, xv, xv, xv, 0);
          }


/**  shouldn't have to change elem dimensions
  dim    = ei[pg->imtrx]->ielem_dim;

  for(b=0; b< Num_Basis_Functions; b++)
    {
      vi = VELOCITY1;
      if(pd_glob[ei[pg->imtrx]->mn]->i[vi] == bfd[b]->interpolation);
         {
           velo_interp = b;
         }
    }
 **/
    
/*   compute J for new element   */


	for(b=0;b<MAX_PDIM;b++){xi_tmp[b]=0;}
	switch (newface)
		{
		case 0:
			xi_tmp[1]=-1.; break;
		case 1:
			xi_tmp[0]=1.; break;
		case 2:
			xi_tmp[1]=1.; break;
		case 3:
			xi_tmp[0]=-1.; break;
		case 4:
			xi_tmp[2]=-1.; break;
		case 5:
			xi_tmp[2]=1.; break;
		}

      err = load_basis_functions( xi_tmp, bfd );
      EH( err, "problem from load_basis_functions");

      err = beer_belly();
      EH( err, "beer_belly");

/*   transform xi to new element orientation plus offset  */

	memset( rot, 0, sizeof(double)*MAX_PDIM*MAX_PDIM);

	for(a=0;a<dim;a++)
	    {
	     for(b=0;b<dim;b++)
		{
		 for(c=0;c<dim;c++)
		     {
			rot[a][c] += bfd[*velo_interp]->J[a][b]
				*Jinv[b][c];
		     }
		}
	     }

/*	scale the rotation tensor and zero elements - we should have 
	one nonzero element per row which is either 1 or -1         */

        for(b = 0 ; b < dim ; b++)
             {
        	max_xi = -1.;
        	i_max_xi = -1;
	        for(a = 0 ; a < dim ; a++)
                {
		tmp=fabs(rot[b][a]);
                if(tmp > max_xi){max_xi = tmp;i_max_xi=a;}
                }
	        for(a = 0 ; a < dim ; a++) 
		    { 
			if(a == i_max_xi){rot[b][a]= rot[b][a]/max_xi;}
			else {rot[b][a]=0.;}
		    }
	     }  
	for(a=0;a<dim;a++)
	    {
	     xi_tmp[a]=2.*xi_tmp[a];
		 for(c=0;c<dim;c++)
		     {
			xi_tmp[a] += rot[a][c]*xi[c];
		     }
	     }
	for(a=0;a<dim;a++){ xi[a]=xi_tmp[a];}

      
           }   /*  end of check left domain */
        }      /*  end of switched elements */

      norm = 0.;
      for(a=0; a<dim; a++)norm+=update[a]*update[a];
      converged = (norm < 1.e-10); 
      inewton++;
    }
      
        return(converged);

  /* invert element mapping and evaluate at the input physical space coordinates */
}
