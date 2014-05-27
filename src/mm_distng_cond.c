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
 *$Id: mm_distng_cond.c,v 1.6 1995-03-10 23:13:00 prschun Exp $
 */

/*
 * Revision history:
 *
 * $Log: not supported by cvs2svn $
 * Revision 1.5  1995/03/07  21:46:24  pasacki
 * (o) created LOG2_term indeces to replace all those ilog2() calls
 * (o) extended most vector dimension loops to VIM(=3) as aid to
 *     collecting reqd terms in nonCartesian coordinates
 * (o) added Harwell MA28 linear solver option (cf Makefile and sl_ma28.c)
 * (o) checked port back to SunOS 4.1.3
 * (o) eliminated various unused routines (eg, heat_conduct_with_mesh)
 * (o) generic orthogonal curvilinear coordinates (note scale factors in
 *     fv structure, bulk integrals, but *not* in boundary conditions)
 * (o) better insulation for fixed grid problems (fixed many UMR's)
 * (o) made a default "point to a zero" for unused esp entries
 * (o) used PROTO to help make code palatable to strict K&R compilers
 * (o) added more test problems
 *
 * Revision 1.4  1994/12/22  22:03:52  racairn
 * This is a 'purified' version of GOMA with boundary conditions separated
 * into four separate functions (the first step in cleaning up this section).
 * Also, this version includes lagrangian mesh motion using linear or
 * nonlinear constitutive equations.
 *
 * Revision 1.3  1994/11/28  20:30:58  racairn
 * Fixed some minor Bugs in the KIN_LEAK condition and implemented the
 * Write Intermediate Results card to print data after each Newton Iteration
 *
 * Revision 1.2  1994/11/03  18:20:26  pasacki
 * AIX 3.2 port; works for simple problems, not all problems.
 *
 * Revision 1.1.1.1  1994/05/24  19:56:30  pasacki
 * initial checkin
 *
 */

static char rcsid[] = "$Id: mm_distng_cond.c,v 1.6 1995-03-10 23:13:00 prschun Exp $";

/* Standard include files */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* SALSA include files */

#include "rf_salsa.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_masks.h"
#include "rf_io_const.h"
#include "rf_io_structs.h"
#include "rf_io.h"
#include "rf_mp.h"
#include "el_elm.h"
#include "el_geom.h"
#include "rf_bc_const.h"
#include "rf_solver_const.h"
#include "rf_fill_const.h"

/*********** R O U T I N E S   I N   T H I S   F I L E *************************
*
*						-All routines in this
*				NOTE:		 file are called only by
*						 rf_fill.c: matrix_fill
*						 except possibly for static
*						 functions.
*
*       NAME			TYPE			CALLED BY
*  -----------------------------------------------------------------
*  Tmelting_pt_distng			void	
*  Tmelting_pt_distng_res		void 		
*
*
*******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void 
Tmelting_pt_distng (irow_index, I, num_local_nodes, iconnect_ptr, ielem_dim,
	      ija, a, x, meshvar)

     int irow_index;          /* Elemental stiffness matrix row index         */
     int I;                   /* Global node number                           */
     int num_local_nodes;     /* number of nodes in the current element       */
     int iconnect_ptr;        /* pointer for this element into the connectivity
				 vector                                       */
     int ielem_dim;           /* the physical dimension of the element,
				 ie., 1, 2, 3                                 */
     int    ija[];            /* column pointer array                         */
     double a[];              /* nonzero array                                */
     double x[];              /* Vector containing the current solution       */
     int  meshvar;            /* integer 1,2, or 3, denoting which displacement
                                 equation is the current one                  */

/*******************************************************************************

  Function which calculates the contribution of the melting point condition
  to the overall stiffness matrix.  Only called for those elements with a 
  side flagged as having a distinguishing condition applied to it.

  Author:          Randy Schunk (1511)
  Date:            8 March 1992
  Revised: 

  ----------------------------------------------------------------------------

  Functions called:

  in_list -- Function which returns location of an integer value in a 
             subsection of a vector
  index_solution --
  Proc_Elem_Connect --

  ----------------------------------------------------------------------------

*******************************************************************************/
     
{
  
/* Local variables */
  
  int j, index_a;
  
  int J;                      /* Global interaction node number               */
  int i_eqn;                  /* Row index into the global stiffness matrix   */
  int j_eqn;                  /* Column index into the global stiffness matrix*/
  
/* Function and externals definitions */
  
  extern int index_solution ();     /* mm_unknown_map.c			      */
  extern int *Proc_Elem_Connect;    /* el_geom.h 			      */
  
/***************************** EXECUTION BEGINS *******************************/
  
  /* find out where in 'a' this contribution goes.  This is simply a map to the
     MSR sparse storage format */

  if ( (i_eqn =  index_solution (I, meshvar, 0, 0)) > -1 ) {

      if ( (j_eqn = index_solution (I, TEMPERATURE, 0, 0)) > -1) {
	index_a = in_list(j_eqn, ija[i_eqn], ija[i_eqn+1], ija);
	if (index_a == -1) {
	  (void) fprintf (stderr, "Tmelting_pt_distng: ERROR - bad index \n");
	  (void) fprintf (stderr, 
			  "Error! j_eqn = %d in processor %d not found in  ",
		          j_eqn, Proc);
	  (void) fprintf (stderr, "list %d thru %d\n", ija[i_eqn], 
		          ija[i_eqn+1]);
	  (void) fprintf (stderr, "\tGlobal node number, I, = %d\n", I); 
	  (void) fprintf (stderr, 
	    		"\tEquation number for the current row, i_eqn, = %d\n", 
			  i_eqn);
	  exit(-1);
	}
      }
      else {
	(void) fprintf (stderr, "Tmelting_pt_distng: ERROR - bad index \n");
	(void) fprintf (stderr, "\tUnable to find proper equation numbers.\n");
	exit(-1);
      }
  } /* END if (i_eqn =  index_solution(I, meshvar, 0, 0)) > -1 )             */
  else {
    (void) fprintf (stderr, "Tmelting_pt_distng: ERROR - bad index \n");
    (void) fprintf (stderr, "\tUnable to find proper equation numbers.\n");
    exit(-1);
  } /* END else (i_eqn =  index_solution(I, TEMPERATURE, 0, 0)) > -1 )           */

  
  /* sum the contributions to the global stiffness matrix */


        a[index_a] += BIG_PENALTY;

  
} /* END of routine Tmelting_pt_distng                                              */

/******************************************************************************//******************************************************************************//******************************************************************************/

void 
Tmelting_pt_distng_res (irow_index, I, num_local_nodes, iconnect_ptr, ielem_dim,
	      x, resid_vector, meshvar)

     int irow_index;          /* Local element node number                    */
     int I;                   /* Local node number                            */
     int num_local_nodes;     /* Number of nodes in the current element       */
     int iconnect_ptr;        /* Pointer for this element into the connectivity
				 vector                                       */
     int ielem_dim;           /* the physical dimension of the element,
				 ie., 1, 2, 3                                 */
     double x[];              /* Vector containing the current solution       */
     double resid_vector[];     /* Residual vector                              */
     int  meshvar;            /* integer 1,2, or 3, denoting which displacement
                                 equation is the current one                  */

/*******************************************************************************

  Function which calculates the contribution of Fickian energy diffusion to
  the residual.


  Functions called:

  in_list -- Function which returns location of an integer value in a 
             subsection of a vector
  index_solution --
  Proc_Elem_Connect --

  ----------------------------------------------------------------------------

*******************************************************************************/
     
{
  
/* Local variables */
  
  int j;
  int J;                      /* Global interaction node number               */
  int i_eqn;                  /* Row index into the global stiffness matrix   */
  int j_eqn;                  /* Column index into the global stiffness matrix*/
  double grad_T[MAX_PDIM];    /* Value of the gradient of T at the local 
				 gauss point				      */
  double TMLTNG_PT=0.30      ;     /* Value of the melting point at which the 
                                 free boundary is located   */

/* Function and externals definitions */
  
  extern int index_solution ();     /* mm_unknown_map.c			      */
  extern int *Proc_Elem_Connect;    /* el_geom.h 			      */
  
/***************************** EXECUTION BEGINS *******************************/
  
  /* find out where in 'a' this contribution goes.  This is simply a map to the
     MSR sparse storage format */

  if ( (i_eqn =  index_solution (I, meshvar, 0, 0)) == -1 ) {
    (void) fprintf (stderr, "Tmelting_pt_distng_res: ERROR - bad index \n");
    (void) fprintf (stderr, "\tUnable to find proper equation number.\n");
    exit(-1);
  }

      if ( (j = index_solution (I, TEMPERATURE, 0, 0))==-1){
        (void) fprintf (stderr, "Tmelting_pt_distng_res: ERROR - bad index \n");
        (void) fprintf (stderr, "\tUnable to find proper equation number.\n");
        exit(-1);
      } else {
        resid_vector[i_eqn] +=  BIG_PENALTY*(x[j]-TMLTNG_PT);
      }
  
} /* END of routine Tmelting_pt_cond_res                               */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
static int
set_temperature_mltng_pt (inode, boundary_condition,TMLTNG_PT )

int 			   inode;
struct Boundary_Condition *boundary_condition;
double TMLTNG_PT      ;     /* Value of the melting point at which the 
                                 free boundary is located   */

/*
*    set_termperature_mltng_pt: 
*  This function sets the value of the melting point used in the
*  distinguishing condition.  This function returns a 0 if the 
*  distinguishing boundary condition is 
*  set successfully.  
*
*/

{


/* Local Variables */
  int  			ieqn, error_cond = 0;

  /* Set melting point from value input on data card */
  TMLTNG_PT = boundary_condition->BC_Data_Float[0];
  printf("SETTING MELTING POINT TO %e\n",TMLTNG_PT);

  return (error_cond);

}
void 
Tmelting_pt_distng_ss (fill_flag, irow_index, I, 
              iconnect_ptr, ielem_dim,
	      ija, a, resid_vector, x, meshvar, 
              num_nodes_on_side, local_elem_node_id, t_melting_pt)

     enum e_fill_func fill_flag;
     int irow_index;          /* Elemental stiffness matrix row index         */
     int I;                   /* Global node number                           */
     int iconnect_ptr;        /* pointer for this element into the connectivity
				 vector                                       */
     int ielem_dim;           /* the physical dimension of the element,
				 ie., 1, 2, 3                                 */
     int    ija[];            /* column pointer array                         */
     double a[];              /* nonzero array                                */
     double resid_vector[];   /* Residual vector			      */
     double x[];              /* Vector containing the current solution       */
     int  meshvar;            /* integer 1,2, or 3, denoting which displacement
                                 equation is the current one                  */
     int   num_nodes_on_side; /* Number of nodes on the side of the element   */
     int local_elem_node_id[MAX_NODES_PER_SIDE];
			      /* Vector of the local elements node numbers
			 	 located on the side of the element 	      */
     double t_melting_pt;     /* melting point temperature to be followed     */

/*******************************************************************************
  Function is the same as Tmelting_pt_distng except used for side sets, 
  not node sets.

  Function which calculates the contribution of the melting point condition
  to the overall stiffness matrix.  Only called for those elements with a 
  side flagged as having a distinguishing condition applied to it.

  Author:          Randy Schunk (1511)
  Date:            8 March 1992
  Revised: 

  ----------------------------------------------------------------------------

  Functions called:

  in_list -- Function which returns location of an integer value in a 
             subsection of a vector
  index_solution --
  Proc_Elem_Connect --

  ----------------------------------------------------------------------------

*******************************************************************************/
     
{
  
/* Local variables */
  
  int j, j_id, index_a;
  
  int J;                      /* Global interaction node number               */
  int i_eqn;                  /* Row index into the global stiffness matrix   */
  int j_eqn;                  /* Column index into the global stiffness matrix*/
  
/* Function and externals definitions */
  
  extern int index_solution ();     /* mm_unknown_map.c			      */
  extern int *Proc_Elem_Connect;    /* el_geom.h 			      */
  
/***************************** EXECUTION BEGINS *******************************/
  
  /* find out where in 'a' this contribution goes.  This is simply a map to the
     MSR sparse storage format */

if (fill_flag == NEWTON) {
  if ( (i_eqn =  index_solution (I, meshvar, 0, 0)) > -1 ) {


      if ( (j_eqn = index_solution (I, TEMPERATURE, 0, 0)) > -1) {
	index_a = in_list(j_eqn, ija[i_eqn], ija[i_eqn+1], ija);
	if (index_a == -1) {
	  (void) fprintf (stderr, "Tmelting_pt_distng_ss: ERROR - bad index \n");
	  (void) fprintf (stderr, 
			  "Error! j_eqn = %d in processor %d not found in  ",
		          j_eqn, Proc);
	  (void) fprintf (stderr, "list %d thru %d\n", ija[i_eqn], 
		          ija[i_eqn+1]);
	  (void) fprintf (stderr, "\tGlobal node number, I, = %d\n", I); 
	  (void) fprintf (stderr, 
	    		"\tEquation number for the current row, i_eqn, = %d\n", 
			  i_eqn);
	  exit(-1);
	}
      }
      else {
	(void) fprintf (stderr, "Tmelting_pt_distng_ss: ERROR - bad index \n");
	(void) fprintf (stderr, "\tUnable to find proper equation numbers.\n");
	exit(-1);
      }
  } /* END if (i_eqn =  index_solution(I, meshvar, 0, 0)) > -1 )             */
  else {
    (void) fprintf (stderr, "Tmelting_pt_distng_ss: ERROR - bad index \n");
    (void) fprintf (stderr, "\tUnable to find proper equation numbers.\n");
    exit(-1);
  } /* END else (i_eqn =  index_solution(I, TEMPERATURE, 0, 0)) > -1 )           */

  
  /* sum the contributions to the global stiffness matrix */

    a[index_a] += BIG_PENALTY;
 
} else {
  i_eqn =  index_solution (I, meshvar, 0, 0);
}

/* Calculate the residual contribution					     */
  if ( (j = index_solution (I, TEMPERATURE, 0, 0))==-1){
    (void) fprintf (stderr, "Tmelting_pt_distng_ss: ERROR - bad index \n");
    (void) fprintf (stderr, "\tUnable to find proper equation number.\n");
    exit(-1);
  } else {
    resid_vector[i_eqn] +=  BIG_PENALTY*(x[j]-t_melting_pt);
  }
  
  
} /* END of routine Tmelting_pt_distng_ss                                     */
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
distng (fill_flag, irow_index, I, 
              iconnect_ptr, ielem_dim,
	      ija, a, resid_vector, x, meshvar, 
              num_nodes_on_side, local_elem_node_id, t_melting_pt)

     enum e_fill_func fill_flag;
     int irow_index;          /* Elemental stiffness matrix row index         */
     int I;                   /* Global node number                           */
     int iconnect_ptr;        /* pointer for this element into the connectivity
				 vector                                       */
     int ielem_dim;           /* the physical dimension of the element,
				 ie., 1, 2, 3                                 */
     int    ija[];            /* column pointer array                         */
     double a[];              /* nonzero array                                */
     double resid_vector[];   /* Residual vector			      */
     double x[];              /* Vector containing the current solution       */
     int  meshvar;            /* integer 1,2, or 3, denoting which displacement
                                 equation is the current one                  */
     int   num_nodes_on_side; /* Number of nodes on the side of the element   */
     int local_elem_node_id[MAX_NODES_PER_SIDE];
			      /* Vector of the local elements node numbers
			 	 located on the side of the element 	      */
     double t_melting_pt;     /* melting point temperature to be followed     */

/*******************************************************************************
  Function is the same as accept used for side sets, 
  not node sets.

  Function which calculates the contribution of the melting point condition
  to the overall stiffness matrix.  Only called for those elements with a 
  side flagged as having a distinguishing condition applied to it.

  Author:          Randy Schunk (1511)
  Date:            8 March 1992
  Revised: 

  ----------------------------------------------------------------------------

  Functions called:

  in_list -- Function which returns location of an integer value in a 
             subsection of a vector
  index_solution --
  Proc_Elem_Connect --

  ----------------------------------------------------------------------------

*******************************************************************************/
     
{
  
/* Local variables */
  
  int j, j_id, index_a;
  
  int J;                      /* Global interaction node number               */
  int i_eqn;                  /* Row index into the global stiffness matrix   */
  int j_eqn;                  /* Column index into the global stiffness matrix*/
  
/* Function and externals definitions */
  
  extern int index_solution ();     /* mm_unknown_map.c			      */
  extern int *Proc_Elem_Connect;    /* el_geom.h 			      */
  
/***************************** EXECUTION BEGINS *******************************/
  
  /* find out where in 'a' this contribution goes.  This is simply a map to the
     MSR sparse storage format */

if (fill_flag == NEWTON) {
  if ( (i_eqn =  index_solution (I, meshvar, 0, 0)) > -1 ) {


      if ( (j_eqn = index_solution (I, TEMPERATURE, 0, 0)) > -1) {
	index_a = in_list(j_eqn, ija[i_eqn], ija[i_eqn+1], ija);
	if (index_a == -1) {
	  (void) fprintf (stderr, "distng: ERROR - bad index \n");
	  (void) fprintf (stderr, 
			  "Error! j_eqn = %d in processor %d not found in  ",
		          j_eqn, Proc);
	  (void) fprintf (stderr, "list %d thru %d\n", ija[i_eqn], 
		          ija[i_eqn+1]);
	  (void) fprintf (stderr, "\tGlobal node number, I, = %d\n", I); 
	  (void) fprintf (stderr, 
	    		"\tEquation number for the current row, i_eqn, = %d\n", 
			  i_eqn);
	  exit(-1);
	}
      }
      else {
	(void) fprintf (stderr, "distng: ERROR - bad index \n");
	(void) fprintf (stderr, "\tUnable to find proper equation numbers.\n");
	exit(-1);
      }
  } /* END if (i_eqn =  index_solution(I, meshvar, 0, 0)) > -1 )             */
  else {
    (void) fprintf (stderr, "distng: ERROR - bad index \n");
    (void) fprintf (stderr, "\tUnable to find proper equation numbers.\n");
    exit(-1);
  } /* END else (i_eqn =  index_solution(I, TEMPERATURE, 0, 0)) > -1 )           */

  
  /* sum the contributions to the global stiffness matrix */
  if ( (j = index_solution (I, TEMPERATURE, 0, 0))==-1){
    (void) fprintf (stderr, "distng: ERROR - bad index \n");
    (void) fprintf (stderr, "\tUnable to find proper equation number.\n");
    exit(-1);
  } else {
    a[index_a] += BIG_PENALTY;
  } 
} else {
  i_eqn =  index_solution (I, meshvar, 0, 0);
}

/* Calculate the residual contribution					     */
  if ( (j = index_solution (I, TEMPERATURE, 0, 0))==-1){
    (void) fprintf (stderr, "distng: ERROR - bad index \n");
    (void) fprintf (stderr, "\tUnable to find proper equation number.\n");
    exit(-1);
  } else {
    resid_vector[i_eqn] +=  BIG_PENALTY*(x[j]-t_melting_pt);
  }
  
  
} /* END of routine distng                                     */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void 
plane_n (fill_flag, irow_index, I, 
              iconnect_ptr, ielem_dim,
	      ija, a, resid_vector, x, meshvar, 
              num_nodes_on_side, local_elem_node_id, a1, a2, a3, a4)

     enum e_fill_func fill_flag;
     int irow_index;          /* Elemental stiffness matrix row index         */
     int I;                   /* Global node number                           */
     int iconnect_ptr;        /* pointer for this element into the connectivity
				 vector                                       */
     int ielem_dim;           /* the physical dimension of the element,
				 ie., 1, 2, 3                                 */
     int    ija[];            /* column pointer array                         */
     double a[];              /* nonzero array                                */
     double resid_vector[];   /* Residual vector			      */
     double x[];              /* Vector containing the current solution       */
     int  meshvar;            /* integer 1,2, or 3, denoting which displacement
                                 equation is the current one                  */
     int   num_nodes_on_side; /* Number of nodes on the side of the element   */
     int local_elem_node_id[MAX_NODES_PER_SIDE];
			      /* Vector of the local elements node numbers
			 	 located on the side of the element 	      */
     double a1, a2, a3, a4  ;         /*  function parameters from data card  */
/*******************************************************************************
  Function is the same as Tmelting_pt_distng except used for side sets, 
  not node sets.

  Function which calculates the contribution of the melting point condition
  to the overall stiffness matrix.  Only called for those elements with a 
  side flagged as having a distinguishing condition applied to it.

  Author:          Randy Schunk (1511)
  Date:            8 March 1992
  Revised: 

  ----------------------------------------------------------------------------

  Functions called:

  in_list -- Function which returns location of an integer value in a 
             subsection of a vector
  index_solution --
  Proc_Elem_Connect --
  plane_fnc --     Function for plane  boundary description
  dplane_fnc_d1 -- Derivative of plane_fnc wrt first coordinate
  dplane_fnc_d2 -- Derivative of plane_fnc wrt second coordinate
  dplane_fnc_d3 -- Derivative of plane_fnc wrt third coordinate

  ----------------------------------------------------------------------------

*******************************************************************************/
{
  
/* Local variables  */
  
  int j, j1, j2, j3, j_id, index_a1, index_a2, index_a3 ;
  
  int J;                      /* Global interaction node number               */
  int i_eqn;                  /* Row index into the global stiffness matrix   */
  int j_eqn;                  /* Column index into the global stiffness matrix*/
  double PI=3.1415927;
  
/* Function and externals definitions */
  
  extern int index_solution ();     /* mm_unknown_map.c			      */
  extern int *Proc_Elem_Connect;    /* el_geom.h 			      */

  /*
   * Insulate prototypes so ANSI and K&R comilers can both handle this...
   */
  dbl plane_fnc		PROTO((
			       dbl c1, 
			       dbl c2, 
			       dbl c3,
			       dbl x1, 
			       dbl x2, 
			       dbl x3, 
			       dbl a1, 
			       dbl a2, 
			       dbl a3, 
			       dbl a4
));    
  
  dbl dplane_fnc_d1	PROTO((
			       dbl c1,
			       dbl c2,
			       dbl c3,
			       dbl x1,
			       dbl x2,
			       dbl x3, 
			       dbl a1,
			       dbl a2,
			       dbl a3,
			       dbl a4
));    

  dbl dplane_fnc_d2	PROTO((
			       dbl c1,
			       dbl c2,
			       dbl c3,
			       dbl x1,
			       dbl x2,
			       dbl x3,
			       dbl a1,
			       dbl a2,
			       dbl a3,
			       dbl a4
));     

  dbl dplane_fnc_d3	PROTO((
			       dbl c1,
			       dbl c2,
			       dbl c3,
			       dbl x1,
			       dbl x2,
			       dbl x3,
			       dbl a1,
			       dbl a2,
			       dbl a3,
			       dbl a4
));     

  
/***************************** EXECUTION BEGINS *******************************/
  
  /* find out where in 'a' this contribution goes.  This is simply a map to the
     MSR sparse storage format */

if (fill_flag == NEWTON) {
  if ( (i_eqn =  index_solution (I, meshvar, 0, 0)) > -1 ) {


      if ( (j_eqn = index_solution (I, MESH_DISPLACEMENT1, 0, 0)) > -1) {
	index_a1 = (i_eqn == j_eqn) ? i_eqn :
	  in_list(j_eqn, ija[i_eqn], ija[i_eqn+1], ija);
	if (index_a1 == -1) {
	  (void) fprintf (stderr, "plane_n: ERROR - bad index \n");
	  (void) fprintf (stderr, 
			  "Error! j_eqn = %d in processor %d not found in  ",
		          j_eqn, Proc);
	  (void) fprintf (stderr, "list %d thru %d\n", ija[i_eqn], 
		          ija[i_eqn+1]);
	  (void) fprintf (stderr, "\tGlobal node number, I, = %d\n", I); 
	  (void) fprintf (stderr, 
	    		"\tEquation number for the current row, i_eqn, = %d\n", 
			  i_eqn);
	  exit(-1);
	}
      }
      else {
	(void) fprintf (stderr, "plane_n: ERROR - bad index \n");
	(void) fprintf (stderr, "\tUnable to find proper equation numbers.\n");
	exit(-1);
      }

      if ( (j_eqn = index_solution (I, MESH_DISPLACEMENT2, 0, 0)) > -1) {
	index_a2 = (i_eqn == j_eqn) ? i_eqn :
	  in_list(j_eqn, ija[i_eqn], ija[i_eqn+1], ija);
	if (index_a2 == -1) {
	  (void) fprintf (stderr, "plane_n: ERROR - bad index \n");
	  (void) fprintf (stderr, 
			  "Error! j_eqn = %d in processor %d not found in  ",
		          j_eqn, Proc);
	  (void) fprintf (stderr, "list %d thru %d\n", ija[i_eqn], 
		          ija[i_eqn+1]);
	  (void) fprintf (stderr, "\tGlobal node number, I, = %d\n", I); 
	  (void) fprintf (stderr, 
	    		"\tEquation number for the current row, i_eqn, = %d\n", 
			  i_eqn);
	  exit(-1);
	}
      }
      else {
	(void) fprintf (stderr, "plane_n: ERROR - bad index \n");
	(void) fprintf (stderr, "\tUnable to find proper equation numbers.\n");
	exit(-1);
      }

      if (ielem_dim == 3) {
        if ( (j_eqn = index_solution (I, MESH_DISPLACEMENT3, 0, 0)) > -1) {
	  index_a3 = (i_eqn == j_eqn) ? i_eqn :
	    in_list(j_eqn, ija[i_eqn], ija[i_eqn+1], ija);
	  if (index_a3 == -1) {
	    (void) fprintf (stderr, "plane_n: ERROR - bad index \n");
	    (void) fprintf (stderr, 
			  "Error! j_eqn = %d in processor %d not found in  ",
		          j_eqn, Proc);
	    (void) fprintf (stderr, "list %d thru %d\n", ija[i_eqn], 
		          ija[i_eqn+1]);
	    (void) fprintf (stderr, "\tGlobal node number, I, = %d\n", I); 
	    (void) fprintf (stderr, 
	    		"\tEquation number for the current row, i_eqn, = %d\n", 
			  i_eqn);
	    exit(-1);
	  }
        } 
        else {
	  (void) fprintf (stderr, "plane_n: ERROR - bad index \n");
	  (void) fprintf (stderr, "\tUnable to find proper equation numbers.\n");
	  exit(-1);
        }
     }

  } /* END if (i_eqn =  index_solution(I, MESH_DISPLACEMENT2, 0, 0)) > -1 )             */
  else {
    (void) fprintf (stderr, "plane_n: ERROR - bad index \n");
    (void) fprintf (stderr, "\tUnable to find proper equation numbers.\n");
    exit(-1);
  } /* END else (i_eqn =  index_solution(I, meshvar, 0, 0)) > -1 )           */

  
  /* sum the contributions to the global stiffness matrix */
    j1  = index_solution (I, MESH_DISPLACEMENT1, 0, 0);
    j2  = index_solution (I, MESH_DISPLACEMENT2, 0, 0);

    a[index_a1] += BIG_PENALTY*dplane_fnc_d1(Coor[0][I],Coor[1][I],0.0,
                                     x[j1],x[j2],0.0,a1,a2,a3,a4);

    
    a[index_a2] += BIG_PENALTY*dplane_fnc_d2(Coor[0][I],Coor[1][I],0.0,
                                     x[j1],x[j2],0.0,a1,a2,a3,a4);

    if(ielem_dim == 3) {
     j3  = index_solution (I, MESH_DISPLACEMENT3, 0, 0); 
     a[index_a3] += BIG_PENALTY*dplane_fnc_d3(Coor[0][I],Coor[1][I],Coor[2][I],
                                     x[j1],x[j2],x[j3],a1,a2,a3,a4);
   }
 
} else {
  i_eqn =  index_solution (I, meshvar, 0, 0);
}

/* Calculate the residual contribution					     */
 if (ielem_dim == 2) {
  if ( (j1 = index_solution (I, MESH_DISPLACEMENT1, 0, 0))==-1
      || (j2 = index_solution(I,MESH_DISPLACEMENT2,0, 0)) == -1){
    (void) fprintf (stderr, "plane_n: ERROR - bad index \n");
    (void) fprintf (stderr, "\tUnable to find proper equation number.\n");
    exit(-1);
  } else {
    
    resid_vector[i_eqn] += BIG_PENALTY* plane_fnc(Coor[0][I],Coor[1][I],0.0,
                                                  x[j1],x[j2],0.0,a1,a2,a3,a4);
  }
}

 if (ielem_dim == 3) {
  if ( (j1 = index_solution (I, MESH_DISPLACEMENT1, 0, 0))==-1
      || (j2 = index_solution(I,MESH_DISPLACEMENT2,0, 0)) == -1
      || (j3 = index_solution(I,MESH_DISPLACEMENT3,0, 0)) == -1){
    (void) fprintf (stderr, "plane_n: ERROR - bad index \n");
    (void) fprintf (stderr, "\tUnable to find proper equation number.\n");
    exit(-1);
  } else {
    

    resid_vector[i_eqn] += BIG_PENALTY* plane_fnc(Coor[0][I],Coor[1][I],Coor[2][I],
                                                  x[j1], x[j2], x[j3], a1,a2,a3,a4);

  
  }
}

  
}
/* END of routine plane_n */

/******************************************************************************/

void 
spline (fill_flag, irow_index, I, 
              iconnect_ptr, ielem_dim,
	      ija, a, resid_vector, x, meshvar, 
              num_nodes_on_side, local_elem_node_id, a1, a2, a3, a4)

     enum e_fill_func fill_flag;
     int irow_index;          /* Elemental stiffness matrix row index         */
     int I;                   /* Global node number                           */
     int iconnect_ptr;        /* pointer for this element into the connectivity
				 vector                                       */
     int ielem_dim;           /* the physical dimension of the element,
				 ie., 1, 2, 3                                 */
     int    ija[];            /* column pointer array                         */
     double a[];              /* nonzero array                                */
     double resid_vector[];   /* Residual vector			      */
     double x[];              /* Vector containing the current solution       */
     int  meshvar;            /* integer 1,2, or 3, denoting which displacement
                                 equation is the current one                  */
     int   num_nodes_on_side; /* Number of nodes on the side of the element   */
     int local_elem_node_id[MAX_NODES_PER_SIDE];
			      /* Vector of the local elements node numbers
			 	 located on the side of the element 	      */
     double a1, a2, a3, a4  ;         /*  function parameters from data card  */
/*******************************************************************************
  Function is the same as Tmelting_pt_distng except used for side sets, 
  not node sets.

  Function which calculates the contribution of the melting point condition
  to the overall stiffness matrix.  Only called for those elements with a 
  side flagged as having a distinguishing condition applied to it.

  Author:          Randy Schunk (1511)
  Date:            8 March 1992
  Revised: 

  ----------------------------------------------------------------------------

  Functions called:

  in_list -- Function which returns location of an integer value in a 
             subsection of a vector
  index_solution --
  Proc_Elem_Connect --
  fnc --           Function for general solid boundary description 
                   (user defined and found in mm_user.c)
  dfncd1 --        Derivative of fnc wrt first coordinate
  dfncd2 --        Derivative of fnc wrt second coordinate
  dfncd3 --        Derivative of fnc wrt third coordinate

  ----------------------------------------------------------------------------

*******************************************************************************/
{
  
/* Local variables  */
  
  int j, j1, j2, j3, j_id, index_a1, index_a2, index_a3 ;
  
  int J;                      /* Global interaction node number               */
  int i_eqn;                  /* Row index into the global stiffness matrix   */
  int j_eqn;                  /* Column index into the global stiffness matrix*/
  double PI=3.1415927;
  
/* Function and externals definitions */
  
  extern int index_solution ();     /* mm_unknown_map.c			      */
  extern int *Proc_Elem_Connect;    /* el_geom.h 			      */

  dbl fnc	PROTO((
		       dbl c1,
		       dbl c2,
		       dbl c3,
		       dbl x1,
		       dbl x2,
		       dbl x3, 
		       dbl a1,
		       dbl a2,
		       dbl a3,
		       dbl a4
));    

  dbl dfncd1	PROTO((
		       dbl c1,
		       dbl c2,
		       dbl c3,
		       dbl x1,
		       dbl x2,
		       dbl x3, 
		       dbl a1,
		       dbl a2,
		       dbl a3,
		       dbl a4
));    

  dbl dfncd2	PROTO((
		       dbl c1,
		       dbl c2,
		       dbl c3,
		       dbl x1,
		       dbl x2,
		       dbl x3,
		       dbl a1,
		       dbl a2,
		       dbl a3,
		       dbl a4
));

  dbl dfncd3	PROTO((
		       dbl c1,
		       dbl c2,
		       dbl c3,
		       dbl x1,
		       dbl x2,
		       dbl x3,
		       dbl a1,
		       dbl a2,
		       dbl a3,
		       dbl a4
));     

  
/***************************** EXECUTION BEGINS *******************************/
  
  /* find out where in 'a' this contribution goes.  This is simply a map to the
     MSR sparse storage format */

if (fill_flag == NEWTON) {
  if ( (i_eqn =  index_solution (I, meshvar, 0, 0)) > -1 ) {


      if ( (j_eqn = index_solution (I, MESH_DISPLACEMENT1, 0, 0)) > -1) {
	index_a1 = (i_eqn == j_eqn) ? i_eqn :
	  in_list(j_eqn, ija[i_eqn], ija[i_eqn+1], ija);
	if (index_a1 == -1) {
	  (void) fprintf (stderr, "spline: ERROR - bad index \n");
	  (void) fprintf (stderr, 
			  "Error! j_eqn = %d in processor %d not found in  ",
		          j_eqn, Proc);
	  (void) fprintf (stderr, "list %d thru %d\n", ija[i_eqn], 
		          ija[i_eqn+1]);
	  (void) fprintf (stderr, "\tGlobal node number, I, = %d\n", I); 
	  (void) fprintf (stderr, 
	    		"\tEquation number for the current row, i_eqn, = %d\n", 
			  i_eqn);
	  exit(-1);
	}
      }
      else {
	(void) fprintf (stderr, "spline: ERROR - bad index \n");
	(void) fprintf (stderr, "\tUnable to find proper equation numbers.\n");
	exit(-1);
      }

      if ( (j_eqn = index_solution (I, MESH_DISPLACEMENT2, 0, 0)) > -1) {
	index_a2 = (i_eqn == j_eqn) ? i_eqn :
	  in_list(j_eqn, ija[i_eqn], ija[i_eqn+1], ija);
	if (index_a2 == -1) {
	  (void) fprintf (stderr, "spline: ERROR - bad index \n");
	  (void) fprintf (stderr, 
			  "Error! j_eqn = %d in processor %d not found in  ",
		          j_eqn, Proc);
	  (void) fprintf (stderr, "list %d thru %d\n", ija[i_eqn], 
		          ija[i_eqn+1]);
	  (void) fprintf (stderr, "\tGlobal node number, I, = %d\n", I); 
	  (void) fprintf (stderr, 
	    		"\tEquation number for the current row, i_eqn, = %d\n", 
			  i_eqn);
	  exit(-1);
	}
      }
      else {
	(void) fprintf (stderr, "spline: ERROR - bad index \n");
	(void) fprintf (stderr, "\tUnable to find proper equation numbers.\n");
	exit(-1);
      }

      if (ielem_dim == 3) {
        if ( (j_eqn = index_solution (I, MESH_DISPLACEMENT3, 0, 0)) > -1) {
	  index_a3 = (i_eqn == j_eqn) ? i_eqn :
	    in_list(j_eqn, ija[i_eqn], ija[i_eqn+1], ija);
	  if (index_a3 == -1) {
	    (void) fprintf (stderr, "spline: ERROR - bad index \n");
	    (void) fprintf (stderr, 
			  "Error! j_eqn = %d in processor %d not found in  ",
		          j_eqn, Proc);
	    (void) fprintf (stderr, "list %d thru %d\n", ija[i_eqn], 
		          ija[i_eqn+1]);
	    (void) fprintf (stderr, "\tGlobal node number, I, = %d\n", I); 
	    (void) fprintf (stderr, 
	    		"\tEquation number for the current row, i_eqn, = %d\n", 
			  i_eqn);
	    exit(-1);
	  }
        } 
        else {
	  (void) fprintf (stderr, "spline: ERROR - bad index \n");
	  (void) fprintf (stderr, "\tUnable to find proper equation numbers.\n");
	  exit(-1);
        }
     }

  } /* END if (i_eqn =  index_solution(I, MESH_DISPLACEMENT2, 0, 0)) > -1 )             */
  else {
    (void) fprintf (stderr, "spline: ERROR - bad index \n");
    (void) fprintf (stderr, "\tUnable to find proper equation numbers.\n");
    exit(-1);
  } /* END else (i_eqn =  index_solution(I, meshvar, 0, 0)) > -1 )           */

  
  /* sum the contributions to the global stiffness matrix */
    j1  = index_solution (I, MESH_DISPLACEMENT1, 0, 0);
    j2  = index_solution (I, MESH_DISPLACEMENT2, 0, 0);

    a[index_a1] += BIG_PENALTY*dfncd1(Coor[0][I],Coor[1][I],0.0,
                                     x[j1],x[j2],0.0,a1,a2,a3,a4);

    
    a[index_a2] += BIG_PENALTY*dfncd2(Coor[0][I],Coor[1][I],0.0,
                                     x[j1],x[j2],0.0,a1,a2,a3,a4);

    if(ielem_dim == 3) {
     j3  = index_solution (I, MESH_DISPLACEMENT3, 0, 0); 
     a[index_a3] += BIG_PENALTY*dfncd3(Coor[0][I],Coor[1][I],Coor[2][I],
                                     x[j1],x[j2],x[j3],a1,a2,a3,a4);
   }
 
} else {
  i_eqn =  index_solution (I, meshvar, 0, 0);
}

/* Calculate the residual contribution					     */
 if (ielem_dim == 2) {
  if ( (j1 = index_solution (I, MESH_DISPLACEMENT1, 0, 0))==-1
      || (j2 = index_solution(I,MESH_DISPLACEMENT2,0, 0)) == -1){
    (void) fprintf (stderr, "spline: ERROR - bad index \n");
    (void) fprintf (stderr, "\tUnable to find proper equation number.\n");
    exit(-1);
  } else {
    
    resid_vector[i_eqn] += BIG_PENALTY* fnc(Coor[0][I],Coor[1][I],0.0,
                                            x[j1],x[j2],0.0,a1,a2,a3,a4);
  }
}

 if (ielem_dim == 3) {
  if ( (j1 = index_solution (I, MESH_DISPLACEMENT1, 0, 0))==-1
      || (j2 = index_solution(I,MESH_DISPLACEMENT2,0, 0)) == -1
      || (j3 = index_solution(I,MESH_DISPLACEMENT3,0, 0)) == -1){
    (void) fprintf (stderr, "spline: ERROR - bad index \n");
    (void) fprintf (stderr, "\tUnable to find proper equation number.\n");
    exit(-1);
  } else {
    

    resid_vector[i_eqn] += BIG_PENALTY* fnc(Coor[0][I],Coor[1][I],Coor[2][I],
                                            x[j1], x[j2], x[j3], a1,a2,a3,a4);

  
  }
}

  
} /* END of routine spline                                                   */
/*****************************************************************************/
/******************************************************************************/
void 
spline_n (fill_flag, irow_index, I, 
              iconnect_ptr, ielem_dim,
	      ija, a, resid_vector, x, meshvar, 
              num_nodes_on_side, local_elem_node_id, a1, a2, a3, a4)

     enum e_fill_func fill_flag;
     int irow_index;          /* Elemental stiffness matrix row index         */
     int I;                   /* Global node number                           */
     int iconnect_ptr;        /* pointer for this element into the connectivity
				 vector                                       */
     int ielem_dim;           /* the physical dimension of the element,
				 ie., 1, 2, 3                                 */
     int    ija[];            /* column pointer array                         */
     double a[];              /* nonzero array                                */
     double resid_vector[];   /* Residual vector			      */
     double x[];              /* Vector containing the current solution       */
     int  meshvar;            /* integer 1,2, or 3, denoting which displacement
                                 equation is the current one                  */
     int   num_nodes_on_side; /* Number of nodes on the side of the element   */
     int local_elem_node_id[MAX_NODES_PER_SIDE];
			      /* Vector of the local elements node numbers
			 	 located on the side of the element 	      */
     double a1, a2, a3, a4  ;         /*  function parameters from data card  */
/*******************************************************************************
  Function is the same as Tmelting_pt_distng except used for side sets, 
  not node sets.

  Function which calculates the contribution of the melting point condition
  to the overall stiffness matrix.  Only called for those elements with a 
  side flagged as having a distinguishing condition applied to it.

  Author:          Randy Schunk (1511)
  Date:            8 March 1992
  Revised: 

  ----------------------------------------------------------------------------

  Functions called:

  in_list -- Function which returns location of an integer value in a 
             subsection of a vector
  index_solution --
  Proc_Elem_Connect --
  fnc --           Function for general solid boundary description 
                   (user defined and found in mm_user.c)
  dfncd1 --        Derivative of fnc wrt first coordinate
  dfncd2 --        Derivative of fnc wrt second coordinate
  dfncd3 --        Derivative of fnc wrt third coordinate

  ----------------------------------------------------------------------------

*******************************************************************************/
{
  
/* Local variables  */
  
  int j, j1, j2, j3, j_id, index_a1, index_a2, index_a3 ;
  
  int J;                      /* Global interaction node number               */
  int i_eqn;                  /* Row index into the global stiffness matrix   */
  int j_eqn;                  /* Column index into the global stiffness matrix*/
  double PI=3.1415927;
  
/* Function and externals definitions */
  
  extern int index_solution ();     /* mm_unknown_map.c			      */
  extern int *Proc_Elem_Connect;    /* el_geom.h 			      */

  dbl fnc	PROTO((
		       dbl c1,
		       dbl c2,
		       dbl c3,
		       dbl x1,
		       dbl x2,
		       dbl x3, 
		       dbl a1,
		       dbl a2,
		       dbl a3,
		       dbl a4
));    

  dbl dfncd1	PROTO((
		       dbl c1,
		       dbl c2,
		       dbl c3,
		       dbl x1,
		       dbl x2,
		       dbl x3, 
		       dbl a1,
		       dbl a2,
		       dbl a3,
		       dbl a4
));    

  dbl dfncd2	PROTO((
		       dbl c1,
		       dbl c2,
		       dbl c3,
		       dbl x1,
		       dbl x2,
		       dbl x3,
		       dbl a1,
		       dbl a2,
		       dbl a3,
		       dbl a4
));

  dbl dfncd3	PROTO((
		       dbl c1,
		       dbl c2,
		       dbl c3,
		       dbl x1,
		       dbl x2,
		       dbl x3,
		       dbl a1,
		       dbl a2,
		       dbl a3,
		       dbl a4
));     

/***************************** EXECUTION BEGINS *******************************/
  
  /* find out where in 'a' this contribution goes.  This is simply a map to the
     MSR sparse storage format */

if (fill_flag == NEWTON) {
  if ( (i_eqn =  index_solution (I, meshvar, 0, 0)) > -1 ) {


      if ( (j_eqn = index_solution (I, MESH_DISPLACEMENT1, 0, 0)) > -1) {
	index_a1 = (i_eqn == j_eqn) ? i_eqn :
	  in_list(j_eqn, ija[i_eqn], ija[i_eqn+1], ija);
	if (index_a1 == -1) {
	  (void) fprintf (stderr, "spline_n: ERROR - bad index \n");
	  (void) fprintf (stderr, 
			  "Error! j_eqn = %d in processor %d not found in  ",
		          j_eqn, Proc);
	  (void) fprintf (stderr, "list %d thru %d\n", ija[i_eqn], 
		          ija[i_eqn+1]);
	  (void) fprintf (stderr, "\tGlobal node number, I, = %d\n", I); 
	  (void) fprintf (stderr, 
	    		"\tEquation number for the current row, i_eqn, = %d\n", 
			  i_eqn);
	  exit(-1);
	}
      }
      else {
	(void) fprintf (stderr, "spline_n: ERROR - bad index \n");
	(void) fprintf (stderr, "\tUnable to find proper equation numbers.\n");
	exit(-1);
      }

      if ( (j_eqn = index_solution (I, MESH_DISPLACEMENT2, 0, 0)) > -1) {
	index_a2 = (i_eqn == j_eqn) ? i_eqn :
	  in_list(j_eqn, ija[i_eqn], ija[i_eqn+1], ija);
	if (index_a2 == -1) {
	  (void) fprintf (stderr, "spline_n: ERROR - bad index \n");
	  (void) fprintf (stderr, 
			  "Error! j_eqn = %d in processor %d not found in  ",
		          j_eqn, Proc);
	  (void) fprintf (stderr, "list %d thru %d\n", ija[i_eqn], 
		          ija[i_eqn+1]);
	  (void) fprintf (stderr, "\tGlobal node number, I, = %d\n", I); 
	  (void) fprintf (stderr, 
	    		"\tEquation number for the current row, i_eqn, = %d\n", 
			  i_eqn);
	  exit(-1);
	}
      }
      else {
	(void) fprintf (stderr, "spline_n: ERROR - bad index \n");
	(void) fprintf (stderr, "\tUnable to find proper equation numbers.\n");
	exit(-1);
      }

      if (ielem_dim == 3) {
        if ( (j_eqn = index_solution (I, MESH_DISPLACEMENT3, 0, 0)) > -1) {
	  index_a3 = (i_eqn == j_eqn) ? i_eqn :
	    in_list(j_eqn, ija[i_eqn], ija[i_eqn+1], ija);
	  if (index_a3 == -1) {
	    (void) fprintf (stderr, "spline_n: ERROR - bad index \n");
	    (void) fprintf (stderr, 
			  "Error! j_eqn = %d in processor %d not found in  ",
		          j_eqn, Proc);
	    (void) fprintf (stderr, "list %d thru %d\n", ija[i_eqn], 
		          ija[i_eqn+1]);
	    (void) fprintf (stderr, "\tGlobal node number, I, = %d\n", I); 
	    (void) fprintf (stderr, 
	    		"\tEquation number for the current row, i_eqn, = %d\n", 
			  i_eqn);
	    exit(-1);
	  }
        } 
        else {
	  (void) fprintf (stderr, "spline_n: ERROR - bad index \n");
	  (void) fprintf (stderr, "\tUnable to find proper equation numbers.\n");
	  exit(-1);
        }
     }

  } /* END if (i_eqn =  index_solution(I, MESH_DISPLACEMENT2, 0, 0)) > -1 )             */
  else {
    (void) fprintf (stderr, "spline_n: ERROR - bad index \n");
    (void) fprintf (stderr, "\tUnable to find proper equation numbers.\n");
    exit(-1);
  } /* END else (i_eqn =  index_solution(I, meshvar, 0, 0)) > -1 )           */

  
  /* sum the contributions to the global stiffness matrix */
    j1  = index_solution (I, MESH_DISPLACEMENT1, 0, 0);
    j2  = index_solution (I, MESH_DISPLACEMENT2, 0, 0);

    a[index_a1] += BIG_PENALTY*dfncd1(Coor[0][I],Coor[1][I],0.0,
                                     x[j1],x[j2],0.0,a1,a2,a3,a4);

    
    a[index_a2] += BIG_PENALTY*dfncd2(Coor[0][I],Coor[1][I],0.0,
                                     x[j1],x[j2],0.0,a1,a2,a3,a4);

    if(ielem_dim == 3) {
     j3  = index_solution (I, MESH_DISPLACEMENT3, 0, 0); 
     a[index_a3] += BIG_PENALTY*dfncd3(Coor[0][I],Coor[1][I],Coor[2][I],
                                     x[j1],x[j2],x[j3],a1,a2,a3,a4);
   }
 
} else {
  i_eqn =  index_solution (I, meshvar, 0, 0);
}

/* Calculate the residual contribution					     */
 if (ielem_dim == 2) {
  if ( (j1 = index_solution (I, MESH_DISPLACEMENT1, 0, 0))==-1
      || (j2 = index_solution(I,MESH_DISPLACEMENT2,0, 0)) == -1){
    (void) fprintf (stderr, "spline_n: ERROR - bad index \n");
    (void) fprintf (stderr, "\tUnable to find proper equation number.\n");
    exit(-1);
  } else {
    
    resid_vector[i_eqn] += BIG_PENALTY* fnc(Coor[0][I],Coor[1][I],0.0,
                                            x[j1],x[j2],0.0,a1,a2,a3,a4);
  }
}

 if (ielem_dim == 3) {
  if ( (j1 = index_solution (I, MESH_DISPLACEMENT1, 0, 0))==-1
      || (j2 = index_solution(I,MESH_DISPLACEMENT2,0, 0)) == -1
      || (j3 = index_solution(I,MESH_DISPLACEMENT3,0, 0)) == -1){
    (void) fprintf (stderr, "spline_n: ERROR - bad index \n");
    (void) fprintf (stderr, "\tUnable to find proper equation number.\n");
    exit(-1);
  } else {
    

    resid_vector[i_eqn] += BIG_PENALTY* fnc(Coor[0][I],Coor[1][I],Coor[2][I],
                                            x[j1], x[j2], x[j3], a1,a2,a3,a4);

  
  }
}

  
} /* END of routine spline_n                                                   */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void 
plane (fill_flag, irow_index, I, 
              iconnect_ptr, ielem_dim,
	      ija, a, resid_vector, x, meshvar, 
              num_nodes_on_side, local_elem_node_id, a1, a2, a3, a4)

     enum e_fill_func fill_flag;
     int irow_index;          /* Elemental stiffness matrix row index         */
     int I;                   /* Global node number                           */
     int iconnect_ptr;        /* pointer for this element into the connectivity
				 vector                                       */
     int ielem_dim;           /* the physical dimension of the element,
				 ie., 1, 2, 3                                 */
     int    ija[];            /* column pointer array                         */
     double a[];              /* nonzero array                                */
     double resid_vector[];   /* Residual vector			      */
     double x[];              /* Vector containing the current solution       */
     int  meshvar;            /* integer 1,2, or 3, denoting which displacement
                                 equation is the current one                  */
     int   num_nodes_on_side; /* Number of nodes on the side of the element   */
     int local_elem_node_id[MAX_NODES_PER_SIDE];
			      /* Vector of the local elements node numbers
			 	 located on the side of the element 	      */
     double a1, a2, a3, a4  ;         /*  function parameters from data card  */
/*******************************************************************************
  Function is the same as Tmelting_pt_distng except used for side sets, 
  not node sets.

  Function which calculates the contribution of the melting point condition
  to the overall stiffness matrix.  Only called for those elements with a 
  side flagged as having a distinguishing condition applied to it.

  Author:          Randy Schunk (1511)
  Date:            8 March 1992
  Revised: 

  ----------------------------------------------------------------------------

  Functions called:

  in_list -- Function which returns location of an integer value in a 
             subsection of a vector
  index_solution --
  Proc_Elem_Connect --
  plane_fnc --     Function for plane  boundary description
  dplane_fnc_d1 -- Derivative of plane_fnc wrt first coordinate
  dplane_fnc_d2 -- Derivative of plane_fnc wrt second coordinate
  dplane_fnc_d3 -- Derivative of plane_fnc wrt third coordinate

  ----------------------------------------------------------------------------

*******************************************************************************/
{
  
/* Local variables  */
  
  int j, j1, j2, j3, j_id, index_a1, index_a2, index_a3 ;
  
  int J;                      /* Global interaction node number               */
  int i_eqn;                  /* Row index into the global stiffness matrix   */
  int j_eqn;                  /* Column index into the global stiffness matrix*/
  double PI=3.1415927;
  
/* Function and externals definitions */
  
  extern int index_solution ();     /* mm_unknown_map.c			      */
  extern int *Proc_Elem_Connect;    /* el_geom.h 			      */

  dbl plane_fnc		PROTO((
			       dbl c1,
			       dbl c2,
			       dbl c3,
			       dbl x1,
			       dbl x2,
			       dbl x3, 
			       dbl a1,
			       dbl a2,
			       dbl a3,
			       dbl a4
));    

  dbl dplane_fnc_d1	PROTO((
			       dbl c1,
			       dbl c2,
			       dbl c3,
			       dbl x1,
			       dbl x2,
			       dbl x3, 
			       dbl a1,
			       dbl a2,
			       dbl a3,
			       dbl a4
));    

  dbl dplane_fnc_d2	PROTO((
			       dbl c1,
			       dbl c2,
			       dbl c3,
			       dbl x1,
			       dbl x2,
			       dbl x3,
			       dbl a1,
			       dbl a2,
			       dbl a3,
			       dbl a4
));     

  dbl dplane_fnc_d3	PROTO((
			       dbl c1,
			       dbl c2,
			       dbl c3,
			       dbl x1,
			       dbl x2,
			       dbl x3,
			       dbl a1,
			       dbl a2,
			       dbl a3,
			       dbl a4
));     

  
/***************************** EXECUTION BEGINS *******************************/
  
  /* find out where in 'a' this contribution goes.  This is simply a map to the
     MSR sparse storage format */

if (fill_flag == NEWTON) {
  if ( (i_eqn =  index_solution (I, meshvar, 0, 0)) > -1 ) {


      if ( (j_eqn = index_solution (I, MESH_DISPLACEMENT1, 0, 0)) > -1) {
	index_a1 = (i_eqn == j_eqn) ? i_eqn :
	  in_list(j_eqn, ija[i_eqn], ija[i_eqn+1], ija);
	if (index_a1 == -1) {
	  (void) fprintf (stderr, "plane: ERROR - bad index \n");
	  (void) fprintf (stderr, 
			  "Error! j_eqn = %d in processor %d not found in  ",
		          j_eqn, Proc);
	  (void) fprintf (stderr, "list %d thru %d\n", ija[i_eqn], 
		          ija[i_eqn+1]);
	  (void) fprintf (stderr, "\tGlobal node number, I, = %d\n", I); 
	  (void) fprintf (stderr, 
	    		"\tEquation number for the current row, i_eqn, = %d\n", 
			  i_eqn);
	  exit(-1);
	}
      }
      else {
	(void) fprintf (stderr, "plane: ERROR - bad index \n");
	(void) fprintf (stderr, "\tUnable to find proper equation numbers.\n");
	exit(-1);
      }

      if ( (j_eqn = index_solution (I, MESH_DISPLACEMENT2, 0, 0)) > -1) {
	index_a2 = (i_eqn == j_eqn) ? i_eqn :
	  in_list(j_eqn, ija[i_eqn], ija[i_eqn+1], ija);
	if (index_a2 == -1) {
	  (void) fprintf (stderr, "plane: ERROR - bad index \n");
	  (void) fprintf (stderr, 
			  "Error! j_eqn = %d in processor %d not found in  ",
		          j_eqn, Proc);
	  (void) fprintf (stderr, "list %d thru %d\n", ija[i_eqn], 
		          ija[i_eqn+1]);
	  (void) fprintf (stderr, "\tGlobal node number, I, = %d\n", I); 
	  (void) fprintf (stderr, 
	    		"\tEquation number for the current row, i_eqn, = %d\n", 
			  i_eqn);
	  exit(-1);
	}
      }
      else {
	(void) fprintf (stderr, "plane: ERROR - bad index \n");
	(void) fprintf (stderr, "\tUnable to find proper equation numbers.\n");
	exit(-1);
      }

      if (ielem_dim == 3) {
        if ( (j_eqn = index_solution (I, MESH_DISPLACEMENT3, 0, 0)) > -1) {
	  index_a3 = (i_eqn == j_eqn) ? i_eqn :
	    in_list(j_eqn, ija[i_eqn], ija[i_eqn+1], ija);
	  if (index_a3 == -1) {
	    (void) fprintf (stderr, "plane: ERROR - bad index \n");
	    (void) fprintf (stderr, 
			  "Error! j_eqn = %d in processor %d not found in  ",
		          j_eqn, Proc);
	    (void) fprintf (stderr, "list %d thru %d\n", ija[i_eqn], 
		          ija[i_eqn+1]);
	    (void) fprintf (stderr, "\tGlobal node number, I, = %d\n", I); 
	    (void) fprintf (stderr, 
	    		"\tEquation number for the current row, i_eqn, = %d\n", 
			  i_eqn);
	    exit(-1);
	  }
        } 
        else {
	  (void) fprintf (stderr, "plane: ERROR - bad index \n");
	  (void) fprintf (stderr, "\tUnable to find proper equation numbers.\n");
	  exit(-1);
        }
     }

  } /* END if (i_eqn =  index_solution(I, MESH_DISPLACEMENT2, 0, 0)) > -1 )             */
  else {
    (void) fprintf (stderr, "plane: ERROR - bad index \n");
    (void) fprintf (stderr, "\tUnable to find proper equation numbers.\n");
    exit(-1);
  } /* END else (i_eqn =  index_solution(I, meshvar, 0, 0)) > -1 )           */

  
  /* sum the contributions to the global stiffness matrix */
    j1  = index_solution (I, MESH_DISPLACEMENT1, 0, 0);
    j2  = index_solution (I, MESH_DISPLACEMENT2, 0, 0);

    a[index_a1] += BIG_PENALTY*dplane_fnc_d1(Coor[0][I],Coor[1][I],0.0,
                                     x[j1],x[j2],0.0,a1,a2,a3,a4);

    
    a[index_a2] += BIG_PENALTY*dplane_fnc_d2(Coor[0][I],Coor[1][I],0.0,
                                     x[j1],x[j2],0.0,a1,a2,a3,a4);

    if(ielem_dim == 3) {
     j3  = index_solution (I, MESH_DISPLACEMENT3, 0, 0); 
     a[index_a3] += BIG_PENALTY*dplane_fnc_d3(Coor[0][I],Coor[1][I],Coor[2][I],
                                     x[j1],x[j2],x[j3],a1,a2,a3,a4);
   }
 
} else {
  i_eqn =  index_solution (I, meshvar, 0, 0);
}

/* Calculate the residual contribution					     */
 if (ielem_dim == 2) {
  if ( (j1 = index_solution (I, MESH_DISPLACEMENT1, 0, 0))==-1
      || (j2 = index_solution(I,MESH_DISPLACEMENT2,0, 0)) == -1){
    (void) fprintf (stderr, "plane: ERROR - bad index \n");
    (void) fprintf (stderr, "\tUnable to find proper equation number.\n");
    exit(-1);
  } else {
    
    resid_vector[i_eqn] += BIG_PENALTY* plane_fnc(Coor[0][I],Coor[1][I],0.0,
                                                  x[j1],x[j2],0.0,a1,a2,a3,a4);
  }
}

 if (ielem_dim == 3) {
  if ( (j1 = index_solution (I, MESH_DISPLACEMENT1, 0, 0))==-1
      || (j2 = index_solution(I,MESH_DISPLACEMENT2,0, 0)) == -1
      || (j3 = index_solution(I,MESH_DISPLACEMENT3,0, 0)) == -1){
    (void) fprintf (stderr, "plane: ERROR - bad index \n");
    (void) fprintf (stderr, "\tUnable to find proper equation number.\n");
    exit(-1);
  } else {
    

    resid_vector[i_eqn] += BIG_PENALTY* plane_fnc(Coor[0][I],Coor[1][I],Coor[2][I],
                                                  x[j1], x[j2], x[j3], a1,a2,a3,a4);

  
  }
}

  
} /* END of routine plane                                                   */
/******************************************************************************/
/*       Functions for solid boundary description                            */
/******************************************************************************/
