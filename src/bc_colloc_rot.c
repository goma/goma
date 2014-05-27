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
 
#include <stdio.h>
 
/* SALSA include files */
 
#include "rf_salsa.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_io_const.h"
#include "rf_io_structs.h"
#include "rf_io.h"
#include "rf_mp.h"
#include "el_elm.h"
#include "el_geom.h"
 
#include "rf_masks.h"
#include "rf_bc_const.h"
#include "rf_solver_const.h"
#include "rf_fill_const.h"
 
#include "rf_chemkin_const.h"
#include "rf_chemkin.h"
 
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
 
#include "mm_eh.h"


int
apply_point_colloc_rotated_bc (ija, a, x, resid_vector, x_old, xdot, delta_t, theta,
         first_elem_side_BC_array, fill_flag,
      ielem, ip_total, ielem_type, num_local_nodes, ielem_dim,
      iconnect_ptr, elem_side_bc, num_total_nodes, local_node_list, local_node_list_ptr)
 
     int    ija[];            /* Vector of integer pointers into the vector, a*/
     double a[];              /* Vector of non-zero entries in the coefficient
                                 matrix                               */
     double x[];              /* Solution vector for the current processor    */
     double x_old[];          /* Solution vector at last time step
                                 (old analog of x)                            */
     double xdot[];           /* xdot of current solution                     */
     double delta_t;          /* current time step size                       */
     double theta;            /* parameter to vary time integration from
                                 explicit (theta = 1) to implicit (theta = 0) */
     double resid_vector[];   /* Residual vector for the current processor    */
     struct elem_side_bc_struct *first_elem_side_BC_array[];
                  /* This is an array of pointers to the first
                 surface integral defined for each element.
                 It has a length equal to the total number
                 of elements defined on the current processor */
     enum e_fill_func fill_flag;
                  /* Flag which indicates the functionality within
                                 the routine :
                RESIDUAL : Calculate the residual, only
                NEWTON   : Calculate residual and matrix*/
  int ielem;          /* element number */
  int ip_total;       /* total number of gauss points */
  int ielem_type;          /* element type */
  int num_local_nodes;
  int ielem_dim;
  int iconnect_ptr;
  struct elem_side_bc_struct *elem_side_bc ;
                 /* Pointer to an element side boundary condition
                structure                     */
  int num_total_nodes;
  int local_node_list[MAX_NP_ELEM];
                             /* List to keep track of nodes at which residual and
                                Jacobian corrections have been made           */
  int *local_node_list_ptr;   /* Pointer to position in local_node_list array  */

 
/*******************************************************************************
  Function which fills the FEM stiffness matrices and the right-hand side
  with contributions from boundary conditions
 
  Author:          Harry Moffat, Scott Hutchinson (1421) and others
                   transferred to separate files by Rich Cairncross
  Date:            12 December 1994
*******************************************************************************/
 
{
/* Local variables */
int ip, w, i, I, ibc, k, l, j, id, ij, index, icount,idj; /* counters */
int ieqn_penalt;
int rotate_mesh;
int rotate_mom;
int err;         /* status variable for functions */
int status = 0;
int ndof=0;           /* Hardwire, for now. */
int bc_input_id;
int ifixed, ii, decrement;
 
double s, t, u;             /* Gaussian-quadrature point locations          */
 
double phi[MAX_NP_ELEM];    /* Value of each of the basis functions at the
                 current quadrature point in local coordinate
                 system                                       */
double grad_phi_L[MAX_NP_ELEM][MAX_PDIM];
                              /* Gradient of the basis function wrt the local
                 coordinates                                  */
double grad_phi[MAX_NP_ELEM][MAX_PDIM];
                              /* Gradient of the basis function wrt the
                 physical coordinates                         */
double grad_phi_old[MAX_NP_ELEM][MAX_PDIM];
                              /* Gradient of the basis function wrt the
                 physical coordinates at last time step       */
double det;                 /* Value of the determinant of the Jacobian
                 transformation                               */
double xi[DIM];       /* Local element coordinates of Gauss point. */
  double x_dot[MAX_PDIM];
  double snormal[MAX_PDIM];   /* Vector holding surface normal components     */
  double dsnormal_dx[MAX_PDIM][MAX_PDIM][MAX_NP_ELEM];
                              /* Array holding derivatives ([i][j][k])
                                 of surface normal
                                 component i wrt displacement j at node k     */
  double stangent[2][MAX_PDIM];/* vector components of 2 mutually
                                     orthogonal surface tangent vectors*/
  double dstangent_dx[2][MAX_PDIM][MAX_PDIM][MAX_NP_ELEM];
                              /* Array holding derivatives ([i][j][k])
                                 of surface tangent vector 1 or 2
                                 component i wrt displacement j at node k     */
  double dsurfdet_dx[MAX_PDIM][MAX_NP_ELEM];
                              /* Array holding derivatives ([i][j])
                                 of surface jacobian determinant
                                 wrt displacement i at node j                 */
  double cart_dsurfdet_dx[MAX_PDIM][MAX_NP_ELEM];
                              /* Array holding derivatives ([i][j])
                                 of surface jacobian determinant
                                 wrt displacement i at node j
  /****************************************************************************/
  /*  These next six entries the same as the previous six, except these       */
  /*  act as dummies to hold the values at the surface centroid.  These c_    */
  /*  variables are used in the residual rotation algorithm in Loops 9' & 9'' */
  /****************************************************************************/
  double c_det;                 /* Value of the determinant of the Jacobian
                 transformation                               */
  double c_snormal[MAX_PDIM];   /* Vector holding surface normal components     */
  double c_dsnormal_dx[MAX_PDIM][MAX_PDIM][MAX_NP_ELEM];
                              /* Array holding derivatives ([i][j][k])
                                 of surface normal
                                 component i wrt displacement j at node k     */
  double c_stangent[2][MAX_PDIM];/* vector components of 2 mutually
                                     orthogonal surface tangent vectors*/
  double c_dstangent_dx[2][MAX_PDIM][MAX_PDIM][MAX_NP_ELEM];
                              /* Array holding derivatives ([i][j][k])
                                 of surface tangent vector 1 or 2
                                 component i wrt displacement j at node k     */
  double c_dsurfdet_dx[MAX_PDIM][MAX_NP_ELEM];
                              /* Array holding derivatives ([i][j])
                                 of surface jacobian determinant
                                 wrt displacement i at node j                 */
  /****************************************************************************/
  double wt;                  /* Quadrature weights
                 units - ergs/(sec*cm*K) = g*cm/(sec^3*K)     */
  double temp;                /* Temporary variable                           */
  double rcoord;

 /* EXTERNAL VARIABLES */
 
  extern struct Boundary_Condition *BC_Types;  /* defined in rf_bc.h          */
                  /* Vector of boundary condition structures input
                 from the input file                  */

  extern int jelly_belly();
  extern int load_elem_var();
  extern int load_elem_grd();


/******************************************************************************/
/*                              BLOCK 9'                                      */
/*     START OF SURFACE LOOPS THAT DON'T REQUIRE INTEGRATION (STRONG SENSE)   */
/*         LOOP OVER SURFACE SIDES WITH A DISTINGUISHING COND. ON MESH        */
/* THIS SECTION ROTATES EQUATIONS ON THOSE SURFACES TO NORMAL-TANGENTIAL FORM */
/******************************************************************************/

         /*  For now chose the centroid of the element surface for the evaluation of */
         /*  the surface normal.                                                     */
         s = 0.0 ; t = 0.0 ;
/*         printf("\nnum_local_nodes = %d\n",num_local_nodes);
         printf("num_local_nodes = %d\n",ielem_dim - 1);
         printf("num_local_nodes = %d\n",ielem_type);
         printf("num_local_nodes = %f\n",phi[0] );
         printf("num_local_nodes = %f\n",grad_phi_L[0][0]);
         printf("num_local_nodes = %f\n",s);
         printf("num_local_nodes = %f\n",t);
         printf("num_local_nodes = %d\n", (int) elem_side_bc->id_side);
         printf("num_local_nodes = %d\n", (int) elem_side_bc->num_nodes_on_side);
         printf("num_local_nodes = %d\n", (elem_side_bc->local_elem_node_id[0]));*/
         fill_surf_shape (num_local_nodes, ielem_dim - 1, ielem_type, phi, 
	              grad_phi_L, s, t, 
		      (int) elem_side_bc->id_side,
		      (int) elem_side_bc->num_nodes_on_side,
		      (elem_side_bc->local_elem_node_id));
         /* calculate the determinant of the surface jacobian */
         det = calc_surf_det (ielem, iconnect_ptr, num_local_nodes, 
		          ielem_dim - 1, grad_phi_L, 
		          (int) elem_side_bc->id_side,
		          (int) elem_side_bc->num_nodes_on_side,
		          (elem_side_bc->local_elem_node_id) );

        /* calculate the components of the surface normal and mesh displacement
                                     derivatives*/
         calc_surf_normal (ielem, iconnect_ptr, num_local_nodes, 
		      ielem_dim - 1, grad_phi_L,
		      (int) elem_side_bc->id_side,
		      (int) elem_side_bc->num_nodes_on_side,
		      (elem_side_bc->local_elem_node_id), 
                      snormal, det, phi, dsnormal_dx, dsurfdet_dx);
         calc_surf_tangent (ielem, iconnect_ptr, num_local_nodes, ielem_dim-1,  
	             grad_phi_L, 
		     (int) elem_side_bc->id_side,
		     (int) elem_side_bc->num_nodes_on_side,
		     (elem_side_bc->local_elem_node_id), 
                     stangent, phi, dstangent_dx, dsurfdet_dx,
                     snormal, dsnormal_dx);
        /* Adjust tangent directions for correct balancing of forces.  This is*/
        /* to assure that the tangents computed at a common material boundary */
        /* between two elements are in the same direction.                    */
         if( (elem_block_id[ielem] % 2) == 1) {
/*           for(i=0; i < 2; i++) { looks fishy - should have a dimension here */
           for(i=0; i < ielem_dim - 1; i++)                             {
              for(k = 0; k< ielem_dim; k++)                 {
                 stangent[i][k] = -stangent[i][k]; 
                 for(l=0; l< ielem_dim; l++)                {
                    for ( j = 0 ; j < num_local_nodes; j++) {
                      dstangent_dx[i][k][l][j] = -dstangent_dx[i][k][l][j];
		     }
	         }
	       }
	    }
	 }

          

        for( ibc=0; (bc_input_id = (int) elem_side_bc-> BC_input_id[ibc]) != -1; ibc++) {

        /* Here we will rotate all x-, y- and z-mesh position residual equations    */
        /* to a coordinate system whose basis includes the normal vector and        */
        /* mutually orthogonal tangent vectors to the current surface on the current*/
        /* element. This eases the application of the boundary conditions below     */
        /*                                                                          */
        /* Loop over the number of local element nodes defined in the
				 current element                                   */
        for (i = 0; i < (int) elem_side_bc->num_nodes_on_side; i++) {
           
        /* Find the local element node number for the current node */
	  id = (int) elem_side_bc->local_elem_node_id[i];



        /* Find the local node number given the local element node number,  'i'     */
	    I = Proc_Elem_Connect[iconnect_ptr + id];

        /* Test to see if this node is at the intersection with another side set.   */
        /* Or with a node set containing a DX or DY "fixed" bc.                     */
        /* If so, assign appropriate equation with ieqn_penalt parameter            */

           ieqn_penalt = 0; rotate_mesh = 1;
           ifixed = Variable_Mask[I].DBCD1;
           ii = in_list(I, 0, i_dup_ptr+1, ss_node_dup_list); 
           if (ii != -1)  { 
             if(dup_list_ss[0][ii] == BC_Types[bc_input_id].BC_ID) {
               ieqn_penalt = 0; rotate_mesh = 0;  
               ifixed = Variable_Mask[I].DBCD1;
             } else if (dup_list_ss[1][ii] ==  BC_Types[bc_input_id].BC_ID) {
               ieqn_penalt = 1; rotate_mesh = 0;
               ifixed = Variable_Mask[I].DBCD2;
             } else if (dup_list_ss[2][ii] == BC_Types[bc_input_id].BC_ID && ielem_dim == 3) {
               ieqn_penalt =  2; rotate_mesh = 0;
               ifixed = Variable_Mask[I].DBCD3;
             } else {
             printf("Problem in mm_fill.c: cannot match BC_ID and current node\n");
	   }
           } 


	  if (BC_Types[bc_input_id].BC_Name == DISTNG_BC && !ifixed) 

               {
                /* Check to see if equations at this node have already been corrected      */


                distng                (fill_flag, id, I, iconnect_ptr,
                                       ielem_dim, ija, a , resid_vector,
                                       x, MESH_DISPLACEMENT1 + ieqn_penalt,
                                       (int) elem_side_bc->num_nodes_on_side,
                                       &(elem_side_bc->local_elem_node_id[0]),
                                      (double) BC_Types[bc_input_id].BC_Data_Float[0]);

		} 
  
	  if (BC_Types[bc_input_id].BC_Name == PLANE_BC && !ifixed) 
               {
                /* Check to see if equations at this node have already been corrected      */

                plane_n                (fill_flag, id, I, iconnect_ptr,
                                       ielem_dim, ija, a , resid_vector,
                                       x, MESH_DISPLACEMENT1+ieqn_penalt,
                                       (int) elem_side_bc->num_nodes_on_side,
                                       &(elem_side_bc->local_elem_node_id[0]),
                                      (double) BC_Types[bc_input_id].BC_Data_Float[0],
                                      (double) BC_Types[bc_input_id].BC_Data_Float[1],
                                      (double) BC_Types[bc_input_id].BC_Data_Float[2],
                                      (double) BC_Types[bc_input_id].BC_Data_Float[3]);

		}  
	  if (BC_Types[bc_input_id].BC_Name == SPLINE_BC && !ifixed) 
	        {
                /* Check to see if equations at this node have already been corrected      */


                spline_n               (fill_flag, id, I, iconnect_ptr,
                                       ielem_dim, ija, a , resid_vector,
                                       x, MESH_DISPLACEMENT1+ieqn_penalt,
                                       (int) elem_side_bc->num_nodes_on_side,
                                       &(elem_side_bc->local_elem_node_id[0]),
                                      (double) BC_Types[bc_input_id].BC_Data_Float[0],
                                      (double) BC_Types[bc_input_id].BC_Data_Float[1],
                                      (double) BC_Types[bc_input_id].BC_Data_Float[2],
                                      (double) BC_Types[bc_input_id].BC_Data_Float[3]);
                } 
         }  /* end for (i=0; i< num_nodes_on_side; i++) */
      } /*(end for ibc) */
/******************************************************************************/

  return(status);
} /* end of routine apply_integrated_bc */
