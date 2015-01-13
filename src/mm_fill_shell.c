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

/* Residual and sensitivity assembly routines for equations on shell elements */ 

/*
 *$Id: mm_fill_shell.c,v 5.62 2010-07-30 21:14:52 prschun Exp $
 */

#ifdef USE_RCSID
static char rcsid[] =
"$Id: mm_fill_shell.c,v 5.62 2010-07-30 21:14:52 prschun Exp $";
#endif

/* Standard include files */
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "std.h"
#include "az_aztec.h"
#include "rf_allo.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_masks.h"
#include "rf_io_const.h"
#include "rf_io.h"
#include "rf_mp.h"
#include "el_elm.h"
#include "el_geom.h"
#include "rf_bc.h"
#include "rf_bc_const.h"
#include "rf_solver_const.h"
#include "rf_solver.h"
#include "rf_fill_const.h"
#include "rf_vars_const.h"
#include "mm_mp_const.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_eh.h"
#include "mm_mp.h"
#include "mm_mp_structs.h"
#include "sl_util.h"
#include "mm_fill_shell.h"
#include "mm_std_models_shell.h"
#include "mm_std_models.h"

#define _MM_FILL_SHELL_C
#include "goma.h"

/*
 * Here is a RECIPE for adding new equation/variable sets to Goma.
 * There are a large number of files involved, so they are listed here.
 * While going through this list, always observe how the entries are
 * made in each place for existing equations and variables, and follow
 * the existing formats (FEF) as much as possible. This can be used as a
 * general procedure or guidelines, but there are a few additional
 * steps for shell elements which are indicated below.
 *
 * NOTE: If the new variable is a vector or tensor, then entries
 *       must be made for each individual component - follow the
 *       examples of VELOCITY or VELOCITY_GRADIENT in this case.
 *
 *  STEP 1:  File rf_fem_const.h
 *    Increase constant V_LAST value by one.
 *    Add a "#define" statement for the new variable to give it a global
 *    name. Assign the integer value one higher than the last (highest)
 *    one previously assigned at the end of the variable list.
 *    Add a corresponding "#define" statement for the new equation at the
 *    end of the equation/residual list. Use the name R_<var name> and
 *    the SAME integer value as the new variable.
 *
 *  STEP 2:  File mm_names.h
 *    Add an entry for the new equation to the EQ_Name structure initializer
 *    list, in the position just after the last previously-defined equation
 *    and just before the species unknowns (R_Y0 to R_Y29). Increase the
 *    (commented) indices below this point by one.
 *    Here is where the short equation label (1 to 3
 *    characters) is assigned, which will be used in the input card. Be
 *    sure to use the same global name as assigned in rf_fem_const.h.
 *    Add a similar entry for the new variable to the Var_Name structure
 *    initializer list between the previously-defined variables and the
 *    "extra" variables (starting with MESH_POSITION1). Update the
 *    (commented) indices of the entries below.
 *    Also add entries to the Exo_Var_Names[] and Var_Units[] structure
 *    initializer lists.
 *
 *  STEP 3:  File mm_as_structs.h
 *    Add entries for the new variable to the following structure
 *    definitions:
 *	Element_Stiffness_Pointers
 *	Element_Variable_Pointers
 *      Field_Variables
 *      Diet_Field_Variables
 *    Include any gradient, divergence, and/or mesh derivative terms which
 *    may be needed for assembly routines.
 *
 *  STEP 4:  File mm_as_alloc.c
 *    Add code to allocate space in Element_Stiffness_Pointers for the new
 *    variable when it is active.
 *
 *  STEP 5:  File mm_fill_ptrs.c
 *    In function load_elem_dofptr(), add a new call to function
 *    load_varType_Interpolation_ptrs() for the new variable.
 *
 *  STEP 6:  File bc_colloc.c
 *    In function load_variable(), add a case for the new variable to
 *    the switch(jvar) block.
 *
 *  STEP 7:  File mm_input.c
 *    rd_eqn_specs, part 1: Add a set_eqn() call for the new equation.
 *    rd_eqn_specs, part 3: Add a set_var() call for the new variable.
 *    -- Before going on --
 *    Decide which terms your equation will use (mass, advection,
 *    diffusion, source, boundary, porous, divergence). The number of
 *    these which will be used is the number of equation term multiplier
 *    (ETM) terms which must be parsed in the input file for this
 *    equation's card.
 *    rd_eqn_specs, part 5: Add a case statement in the appropriate place
 *    for the number of ETM terms to be parsed.
 *
 *  STEP 8:  File mm_input_util.c
 *    In function variable_string_to_int(), add a line for the new variable.
 *
 *  STEP 9:  File mm_prob_def.c
 *    In function setup_pd(), add lines of code to set member e of the
 *    Problem Description (pd->) structure. Include one entry for each
 *    ETM term which is used. This assignment is tricky since bitwise
 *    operators are used, so refer carefully to existing entries. 
 *
 *  STEP 10:  File mm_fill_terms.c
 *    In function load_fv(), add a block of code to load the new variable
 *    into the fv structure. See the comments there about zeroing out
 *    previous values.
 *    In function load_fv_grads(), add code to load the gradient terms for
 *    the new variable, if applicable and added to structures in Step 3.
 *    In function load_fv_mesh_derivs(), add code to load the mesh
 *    derivative terms for the new variable, if applicable and added to
 *    structures in Step 3.
 *    
 *
 *  STEP 11:  File mm_flux.c
 *    In function load_fv_sens(), add code to load the new variable into
 *    the fv_sens structure.
 *
 *  STEP 12:  File mm_unknown_map.c
 *    In function set_interaction_masks(), add a case for the new equation
 *    to the switch(e) block. Include entries for the new variable and
 *    any other variables which may appear in the new equation.
 *    Mesh displacements should almost always be included.
 *    If any other equations or BC's will be modified to use the new
 *    variable, add entries for this variable to the cases for each
 *    affected equation, including those which may use a modified BC.
 *
 *  STEP 13:  File mm_fill_shell.c or mm_fill.c
 *    Put the assembly function for the new equation here, to load
 *    lec->R and lec->J as appropriate. Use standard methods from other
 *    equation assembly functions. See notes below for shell equations.
 *    Don't forget to add a call to this assemble function in mm_fill.c. 
 *
 *  STEP 14:  File mm_fill_shell.h or mm_fill.h or mm_fill_rs.h, etc. 
 *    Provide a prototype for the new equation assembly function.
 *
 *  STEP 15:  File ac_stability_util.c
 *    In function modify_fv_mesh_derivs_for_LSA_3D_of_2D(), add the
 *    appropriate lines of code for the new variable, perhaps by
 *    copying the case for a variable of the same dimension.
 *
 *  STEP 16:  Files as needed
 *    If the new variable is to be used in any other equation or
 *    boundary condition, modify the appropriate assembly functions.
 *    See notes below for shell elements.
 *
 *  STEP 17:  File rf_util.c
 *    This is where variable types are grouped for calculation of
 *    variable norms for time integration. Once this section is brought
 *    up to date, some instructions will be provided here.
 *    Until then, nothing is necessary here.
 *
 *  STEP 18:  Use and enjoy your new equation and variable!
 *
 * SPECIAL NOTES FOR SHELL ELEMENTS:
 *  When creating an equation/variable pair to be applied on a block of
 *  shell elements, and the equation will use bulk variables,
 *  place the new assembly function in this file below. Include the
 *  following in the function argument list:
 *      const double wt,        (Gauss point weight)
 *      double xi[DIM],         (stu coordinates for local [bulk] element)
 *	const Exo_DB *exo       (Exodus database ptr)
 *
 *  There are two general procedures for writing shell equation assembly
 *  functions. The first applies to cases where the equation depends on
 *  variables in only one bulk phase. For this case, the function should be
 *  similar to those for other equations, but must also include the following:
 *      Create a local array xi2[DIM].
 *      Create and allocate local array int *n_dof, which will need
 *	  MAX_VARIABLE_TYPES entries.
 *      Get remote (shell) element number from global elem_friends array.
 *      Call function load_neighbor_var_data(), which will require the
 *        above arguments to be passed in.
 *      Assemble the Residual and Jacobian terms AFTER the above steps.
 *      Free the n_dof array before returning.
 *
 *  For the assembly, the field variables on both element blocks are
 *  available in the usual fv-type structures. The shell elements will
 *  have their own basis functions, which can be accessed via bf[var]
 *  as usual. The n_dof array will contain the dof counts for the
 *  remote (shell) element variables, so when looping over these
 *  variable dofs, use n_dof[var] rather than ei->dof[var].
 *
 *  The second procedure is used when the shell equation depends on variables
 *  on both sides, such as a "jump" condition. This will require a special
 *  three-pass assembly: The contributions for each side are assembled
 *  separately in the first two passes and stored in temporary res[] and
 *  jac[][][] arrays while the "shop" is setup on the respective bulk
 *  elements on each side. Then, a third pass on the shell element is done
 *  in which the shell variable contributions are added in and terms such
 *  as wt, h3, and det_J are included. Then, the terms are transferred into
 *  the lec arrays. This procedure does not use load_neighbor_var_data(),
 *  but requires some additional initializations. See the function
 *  assemble surface_charge() below as an example.
 *
 *  If a shell equation uses the GRADIENT of a bulk variable, a special
 *  procedure must be used to enable the cross-block sensitivities, as
 *  there will be nonzero contributions from DOFs of the variable which
 *  do not lie on the shell interface. In this case, it will be necessary
 *  to split up the equation into parts which depend on bulk variable
 *  gradients, and handle those separately. A special type of boundary
 *  condition (WEAK_SHELL_GRAD) has been implemented which assembles
 *  these parts of shell equations as boundary conditions and are
 *  invoked from the bulk side, so that the sensitivity terms can be
 *  transported to the Jacobian through lec->J. See the function 
 *  surface_electric_field_bc() below as an example. These functions
 *  are called from bc_special.c; see the comments there.
 *
 * SPECIAL NOTES FOR MATERIAL PROPERTIES:
 * 
 *  If it is necessary to add new material cards, be sure to add
 *  appropriate lines to the file dp_vif.c.  Follow existing entries
 *  to understand what to add.
 * 
 */

int InShellElementWithParentElementCoverage = 0;
int ShellElementParentElementCoverageForVariable[MAX_VARIABLE_TYPES] = {MAX_VARIABLE_TYPES * 0};


/*****************************************************************************
 *									     *
 *               PLACE EQUATIONS FOR SHELL ELEMENTS HERE.                    *
 *									     *
 *****************************************************************************/

/******************************************************************************
 * assemble_shell_a - Assembles the residual and Jacobian equations for a
 *                      simple (meaningless) equation for testing shell
 *                      elements.
 *
 * Currently, the following equation is solved.  This is likely to change as
 * different things are tested.  "Use the source, Luke".
 *
 *    A = 2V
 *
 * where A is the shell variable and V is the voltage potential.
 *
 * Input
 * =====
 * time_value = The current time.
 * theta      = The implicit-explicit time stepping parameter.
 * delta_t    = The current step size.
 * 
 * Output
 * ======
 * (none)
 * 
 * Returns
 * ======
 * 0  = Success
 * -1 = Failure
 *
 * Revision History
 * ================
 * 7 May 2002 - Patrick Notz - Creation.
 *
 ******************************************************************************/
/*
 * This test function was set up by PKN. It is being retained
 * for use as an example for creating new shell equations.
 */
/*Sa
  int
  assemble_shell_a(double time_value,  Sa*/ /* Time */
/*Sa		 double theta,       Sa*/ /* Time stepping parameter */
/*Sa		 double delta_t,     Sa*/ /* Time step size */
/*Sa             double xi[DIM],     Sa*/ /* Local stu coordinates */
/*Sa                 const Exo_DB *exo)
  {
  int i, j, peqn, var, pvar;
  double phi_i, phi_j;
  double res, jac;
  int eqn = R_SHELL_A; Sa*/

/* Declare some neighbor structures */
/*Sa  int el1 = ei->ielem;
  int el2, nf, err;
  int *n_dof = NULL;  Sa*/

/* Unpack variables from structures for local convenience. */
/*Sa  double wt    = fv->wt;
  double h3    = fv->h3;
  double det_J = bf[eqn]->detJ; Sa*/
  
/* See if there is a friend for this element */
/*Sa  nf = num_elem_friends[el1];
  if (nf == 0) return 0; Sa*/

/*
 * Get neighbor element number.
 * NOTE: For now, only one friend can be handled. If there is
 *       more than one, the first one will be used. This 
 *       will be fixed at a later date.
 */
/*Sa  el2 = elem_friends[el1][0];
  if (nf != 1) WH(-1, "WARNING: Not set up for more than one element friend!");  Sa*/

/*
 * DOF counts will be needed to construct sensitivities of
 * local shell equations to neighbor bulk element variables.
 * Allocate this array to save them in.
 */
/*Sa  n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int)); Sa*/

/*
 * This call will collect all field variable and basis function data
 * for the neighbor element "el2" and populate the fv structure,
 * then repopulate the local element structure. Thus, as long as
 * no variable is defined on both blocks, all data from both
 * elements will be available in the fv structure.
 * For example, pressure in the neighbor element is accessed with fv->P.
 * This is done to simplify the job of writing shell equations.
 */
/*Sa  err = load_neighbor_var_data(el1, el2, n_dof, dof_map, -1, xi, exo); Sa*/

/* Assemble the residual equation */
/* Here, variables from both elements are used. */ 
/*Sa  if ( af->Assemble_Residual )
  {
  peqn = upd->ep[eqn];
  var  = SHELL_A;

  for(i = 0; i < ei->dof[eqn]; i++)
  {
  phi_i = bf[eqn]->phi[i];

  res = (fv->A - 2.0 * fv->V) * phi_i * wt * det_J * h3;
  lec->R[peqn][i] += res;
  }
  }  Sa*/

/* Assemble the sensitivity equations */
/*Sa  if ( af->Assemble_Jacobian )
  {
  peqn = upd->ep[eqn];
  for(i = 0; i < ei->dof[eqn]; i++)
  {
	  
  phi_i = bf[eqn]->phi[i];  Sa*/

/* d_(SHELL_A)_d_A:  LOCAL sensitivity */
/*Sa	  var = SHELL_A;
  if ( pd->v[var] )
  {
  pvar = upd->vp[var];   Sa*/

/* Here, the unknown is defined on the local element.
 * So get the DOF count from ei->dof. */
/*Sa	      for(j = 0; j < ei->dof[var]; j++) 
  {
  phi_j = bf[var]->phi[j];
  jac = phi_j * phi_i * wt * det_J * h3;
  lec->J[peqn][pvar][i][j] += jac;
  }
  }  Sa*/

/* d_(SHELL_A)_d_V:  REMOTE senstivity */
/* Here, the unknown is defined on the neighbor element. */
/*Sa	  var = VOLTAGE;
  if ( n_dof[var] > 0)  Sa*/ /* NOTE: Cannot use pd->v here! */
/*Sa	    {
  pvar = upd->vp[var];   Sa*/

/* Here, the unknown is defined on the remote element.
 * So get the DOF count from n_dof. */
/*Sa	      for(j = 0; j < n_dof[var]; j++) 
  {
  phi_j = bf[var]->phi[j];
  jac = -2.0 * phi_j * phi_i * wt * det_J * h3;
  lec->J[peqn][pvar][i][j] += jac;
  }
  }
  }
  } Sa*/

/* Clean up */
/*Sa  safe_free((void *) n_dof);

return(0);
} Sa*/
/* End of assemble_shell_a() */


/******************************************************************************
 * assemble_surface charge - Assembles the residual and Jacobian equations
 *                             for surface charge distribution on a
 *                             shell element block.
 *                      elements.
 *
 * Currently, the following equation is invoked for testing purposes.
 *               __
 *    d(qs)/dt + \/.j + [n.J] = 0 , [n.J] is constructed in surface_electric_field_bc
 *
 * where n is the surface notmal and D is the electric displacement vector
 * which is permittivity * grad V and qs is the surface charge unknown.
 *
 *
 * Input
 * =====
 * time_value = The current time.
 * theta      = The implicit-explicit time stepping parameter.
 * delta_t    = The current step size.
 * 
 * Output
 * ======
 * (none)
 * 
 * Returns
 * ======
 * 0  = Success
 * -1 = Failure
 *
 * Revision History
 * ================
 * 29 May 2003 - Edward Wilkes - Creation.
 *
 ******************************************************************************/

int
assemble_surface_charge(double time_value,  /* Time */
		        double theta,       /* Time stepping parameter */
		        double delta_t,     /* Time step size */
                        const double wt,    /* Gauss point weight */
                        double xi[DIM],     /* Local stu coordinates */
                        const Exo_DB *exo,
			const int eqn)
{
  int i, j, p, peqn, var, pvar, node, index;
  double phi_i, phi_j, phi_t;
  double d_phi_dxi[MDE], d_qs_ds, d_phi_ds[MDE],dx_dxi[DIM], xs[DIM];
  double surf_diff = 0.0;
  double mass = 0.0;
  double diffusion = 0.0;
  double res[MDE], jac[MDE][MAX_PROB_VAR][MDE];
  int *n_dof = NULL;
  /*  int eqn = R_SURF_CHARGE;*/

  /* Declare some neighbor structures */
  int el0 = ei->ielem;
  int nf;

  /* These are needed to get convection velocity and mesh derivatives */

  double dh3J_dm[DIM][MDE];

  /*
   * These will point to the pd structures for the shell block and
   * one or two neighbor bulk blocks, to simplify keeping track.
   */
  PROBLEM_DESCRIPTION_STRUCT *pd0;

  /* Unpack variables from structures for local convenience. */
  double h3;
  double det_J;
  peqn = upd->ep[eqn];
  
  /* See if there are friends for this element (maximum 2) */
  nf = num_elem_friends[el0];
  if (nf == 0) return 0;

  /* Initialize temporary arrays */
  memset (res, 0, sizeof(double)*MDE);
  memset (jac, 0, sizeof(double)*MDE*MAX_PROB_VAR*MDE);
  memset (dh3J_dm, 0, sizeof(double)*DIM*MDE);

  pd0 = pd;

  /*
   * DOF counts will be needed to construct sensitivities of
   * local shell equations to neighbor bulk element variables.
   * Allocate this array to save them in.
   */
  n_dof = (int *)array_alloc(1, MAX_VARIABLE_TYPES, sizeof(int));

  /* Load coordinates from current element and initialize xi2 */
  for (i=0; i<DIM; i++)
    {
      xs[i] = 0.0;
    }

  /* Get neighbor element number(s). */
  if (nf > 2) EH(-1, "Not set up for more than two element friends!");

  /* Now return to original element (el0) on shell block to finish assembly */
  setup_shop_at_point(el0, xi, exo);

  /* compute shell basis functions and coordinates  */
  for (i = 0; i < ei->dof[eqn]; i++)
    {
      phi_i = bf[eqn]->phi[i];
      d_phi_dxi[i] = bf[eqn]->dphidxi[i][0];
    }
  dx_dxi[0] = dx_dxi[1] = 0.;
  if (ei->deforming_mesh )
    {  
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  node = ei->dof_list[R_MESH1][i];
	  index = Proc_Elem_Connect[ei->iconnect_ptr +node];
          //HKM dx_dxi[i] is exactly the same thing as J[0][i], which is calculated in beer_belly()
          //HKM   thinking that this is all a duplicate of beer_belly functionality
	  dx_dxi[0] +=    (Coor[0][index] + *esp->d[0][i]) * d_phi_dxi[i];
	  dx_dxi[1] +=    (Coor[1][index] + *esp->d[1][i]) * d_phi_dxi[i];
	  xs[0] +=    (Coor[0][index] + *esp->d[0][i]) * bf[eqn]->phi[i];
	  xs[1] +=    (Coor[1][index] + *esp->d[1][i]) * bf[eqn]->phi[i];
	}
    }
  else
    {
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  node = ei->dof_list[eqn][i];
	  index = Proc_Elem_Connect[ei->iconnect_ptr +node];
	  dx_dxi[0] +=    Coor[0][index] * d_phi_dxi[i];
	  dx_dxi[1] +=    Coor[1][index] * d_phi_dxi[i];
	  xs[0] +=    Coor[0][index] * bf[eqn]->phi[i];
	  xs[1] +=    Coor[1][index] * bf[eqn]->phi[i];
	}
    }
 
  //HKM  Thinking this is exactly the same calculation as in beer_belly 
  det_J = sqrt(dx_dxi[0]*dx_dxi[0] + dx_dxi[1]*dx_dxi[1]);
  h3 = 1.;
  if(upd->CoordinateSystem == CYLINDRICAL ||
     upd->CoordinateSystem == SWIRLING)	h3 = xs[1];
  if (ei->deforming_mesh )
    {  
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  dh3J_dm[0][i] += h3*dx_dxi[0]*d_phi_dxi[i]/det_J;
	  dh3J_dm[1][i] += h3*dx_dxi[1]*d_phi_dxi[i]/det_J;
	}
      if(upd->CoordinateSystem == CYLINDRICAL ||
	 upd->CoordinateSystem == SWIRLING)
	{  
	  for (i = 0; i < ei->dof[eqn]; i++)
	    {
	      dh3J_dm[1][i] += det_J*bf[eqn]->phi[i];
	    }
	}  
    }  
 
  for (i = 0; i < ei->dof[eqn]; i++)
    { d_phi_ds[i] = d_phi_dxi[i]/det_J; }
 
/*  gather any needed material properties, etc.   */

  switch(eqn)
    {
    case R_SURF_CHARGE:
      surf_diff = mp->elect_surf_diffusivity;
      break;
    case R_SHELL_USER:
      break;
    case R_SHELL_BDYVELO:
      break;
    case R_SHELL_LUBP:
      break;
    default:
      EH(-1,"shell equation not present");
      break;
    }
  switch(eqn)
    {
    case R_SURF_CHARGE:
  
      d_qs_ds = 0;
      for (i = 0; i < ei->dof[eqn]; i++)
	{ d_qs_ds += *esp->qs[i] * d_phi_ds[i]; }
 
      /* Include residual contributions from shell variables */
      if (af->Assemble_Residual)
	{
	  for (i = 0; i < ei->dof[eqn]; i++)
	    {
	      phi_i = bf[eqn]->phi[i];

	      /* Mass term */
	      if (pd0->e[eqn] & T_MASS)
		{
 		  mass = phi_i * fv_dot->qs;
		  mass *= pd0->etm[eqn][(LOG2_MASS)];
 		  diffusion = d_phi_ds[i] * surf_diff * d_qs_ds;
  		  diffusion *= pd0->etm[eqn][(LOG2_DIFFUSION)];
		}
	      res[i] += (mass+diffusion)*wt*h3*det_J;
	      lec->R[peqn][i] += res[i];
	    }
	}

      /* Include Jacobian contributions from shell variables */
      if (af->Assemble_Jacobian)
	{
	  for (i = 0; i < ei->dof[eqn]; i++)
	    {
	      phi_i = bf[eqn]->phi[i];

	      /* J_qs_qs:  Shell sensitivity */
	      var = SURF_CHARGE;
	      pvar = upd->vp[var];
	      for (j = 0; j < ei->dof[var]; j++)
		{
		  phi_j = bf[var]->phi[j];
		  phi_t = phi_j;

		  /* Mass term only */
		  mass = 0.0;
		  if (pd->TimeIntegration != STEADY)
		    {
		      if (pd0->e[eqn] & T_MASS)
			{
			  phi_t *= ( (1.0 + 2.0 * theta) / delta_t);
			  mass += phi_i*phi_t;
			  mass *= pd0->etm[eqn][(LOG2_MASS)];
			}
		    }
		  diffusion = 0.0;
 		  if (pd0->e[eqn] & T_MASS)
 		    {
 		      diffusion += d_phi_ds[i] * surf_diff * d_phi_ds[j];
  		      diffusion *= pd0->etm[eqn][(LOG2_DIFFUSION)];
 		    }
		  jac[i][pvar][j] += (mass + diffusion)*wt*h3*det_J;
		  lec->J[peqn][pvar][i][j] += jac[i][pvar][j];
		}

	      /* J_qs_x:  Mesh displacement sensitivity term completion & transfer */
	      for (p = 0; p < pd->Num_Dim; p++)
		{
		  var = MESH_DISPLACEMENT1 + p;;
		  if (ei->deforming_mesh )
		    {
		      pvar = upd->vp[var];
		      for (j = 0; j < ei->dof[var]; j++)
			{
			  double diff_dx = 0.0; 
			  double det_J_dmesh;
			  /* Mass term */
			  if (pd0->e[eqn] & T_MASS)
			    {
			      mass = phi_i*fv_dot->qs;
			      mass *= pd0->etm[eqn][(LOG2_MASS)];
			      det_J_dmesh = dx_dxi[p]*d_phi_dxi[j]/det_J;     
			      diffusion = d_phi_ds[i]*surf_diff*d_qs_ds;
			      diffusion *= pd0->etm[eqn][(LOG2_DIFFUSION)];
			      diff_dx = surf_diff*d_phi_ds[i]*d_qs_ds*
                                        (-det_J_dmesh/det_J);
			      diff_dx *= pd0->etm[eqn][(LOG2_DIFFUSION)];
			    }
			  jac[i][pvar][j] += ((mass + diffusion) * wt * dh3J_dm[p][j]
					      + diff_dx * wt * h3 * det_J);
			  lec->J[peqn][pvar][i][j] += jac[i][pvar][j];
			}
		    }
		}

	    }
	}
      break;
    case R_SHELL_BDYVELO:
      d_qs_ds = 0;
      for (i = 0; i < ei->dof[eqn]; i++)
	{ d_qs_ds += *esp->sh_bv[i] * d_phi_ds[i]; }
      if (af->Assemble_Residual)
	{
	  for (i = 0; i < ei->dof[eqn]; i++)
	    {
	      phi_i = bf[eqn]->phi[i];

	      /* Mass term */
              if (pd0->e[eqn] & T_MASS)
                {
                  mass = phi_i*(fv->sh_bv);
                  mass *= pd0->etm[eqn][(LOG2_MASS)];
                  diffusion = 0.0;
                  diffusion *= pd0->etm[eqn][(LOG2_DIFFUSION)];
                }
	      res[i] += (mass+diffusion)*wt*h3*det_J;
	      lec->R[peqn][i] += res[i];
	    }
	}
      /* Include Jacobian contributions from shell variables */
      if (af->Assemble_Jacobian)
	{
	  for (i = 0; i < ei->dof[eqn]; i++)
	    {
	      phi_i = bf[eqn]->phi[i];

	      /* J_qs_qs:  Shell sensitivity */
	      var = SHELL_BDYVELO;
	      pvar = upd->vp[var];
	      for (j = 0; j < ei->dof[var]; j++)
		{
		  phi_j = bf[var]->phi[j];
		  phi_t = phi_j;

		  /* Mass term only */
		  mass = 0.0;
                  if (pd0->e[eqn] & T_MASS)
                    {
                      mass += phi_i*phi_t;
                      mass *= pd0->etm[eqn][(LOG2_MASS)];
                    }
		  diffusion = 0.0;
                  if (pd0->e[eqn] & T_MASS)
                    {
                      diffusion *= pd0->etm[eqn][(LOG2_DIFFUSION)];
                    }
		  jac[i][pvar][j] += (mass + diffusion)*wt*h3*det_J;
		  lec->J[peqn][pvar][i][j] += jac[i][pvar][j];
		}
	      /* J_qs_x:  Mesh displacement sensitivity term completion & transfer */
	      for (p = 0; p < pd->Num_Dim; p++)
		{
		  var = MESH_DISPLACEMENT1 + p;;
		  if (ei->deforming_mesh )
		    {
		      pvar = upd->vp[var];
		      for (j = 0; j < ei->dof[var]; j++)
			{
			  // double diff_dx, det_J_dmesh;
			  /* Mass term */
			  if (pd0->e[eqn] & T_MASS)
			    {
			      mass = phi_i*(fv->sh_bv);
			      mass *= pd0->etm[eqn][(LOG2_MASS)];
			      diffusion = 0.0;
			      diffusion *= pd0->etm[eqn][(LOG2_DIFFUSION)];
			    }
			  jac[i][pvar][j] += (mass+diffusion)*wt*dh3J_dm[p][j];
			  lec->J[peqn][pvar][i][j] += jac[i][pvar][j];
			}
		    }
		}
	    }
	}
      break;
    case R_SHELL_USER:

#include "user_shell.h"
 
      break;
    }
 
  /* Clean up */
  safe_free((void *) n_dof);

  return(0);
} /* assemble_surface_charge() */

/******************************************************************************
 * assemble_shell_structure - Assembles the residual and Jacobian equations
 *                             for inextensible shell structure equations.
 *                            THIS IS A 2D IMPLEMENTATION ONLY! This means that
 *                             h3 scale factor is 1.0 and we don't bother with
 *                             mesh sensitivities of h3. 
 *
 * Currently, the following equations are implemented.
 *
 *    d_2(K)/ds_2 + KT +nn:T =0
 *    dT/ds + KdK/ds + nt:T = 0
 *
 * where K is the surface curvature and T is the tension.  s is the 
 * arclength along the surface.
 *
 *
 * Input
 * =====
 * time_value = The current time.
 * theta      = The implicit-explicit time stepping parameter.
 * delta_t    = The current step size.
 * 
 * Output
 * ======
 * (none)
 * 
 * Returns
 * ======
 * 0  = Success
 * -1 = Failure
 *
 * Revision History
 * ================
 * 22 October 2003 - P. R. Schunk- Creation.
 *
 ******************************************************************************/
int
assemble_shell_structure(double time_value,  /* Time */
			 double theta,       /* Time stepping parameter */
			 double delta_t,     /* Time step size */
			 const double wt,    /* Gauss point weight */
			 double xi[DIM],     /* Local stu coordinates */
			 const Exo_DB *exo)   
{
  int i, j, peqn, var, pvar;
  double phi_i, phi_j;

  /* note the following definitions restrict this to 1D bar elements */
  double d_phi_dxi[MDE], d_sh_K_dxi, d_sh_tens_dxi, d_sh_x_dxi, d_sh_y_dxi; 
  double diffusion;
  double res[MDE], jac[MDE][MAX_PROB_VAR][MDE];
  int *n_dof = NULL;
  int eqn;
  double d_det_J_dmeshbj, det_J_sh;

  PROBLEM_DESCRIPTION_STRUCT *pd0;
  int node, index;

  /* After reviewing the methodology presented in the original shell equations
   * addressing surface charge, and considering that the structural shell equations
   * are more complicated in the sense that they deal with spatial derivatives in 
   * the surface, we will take the following approach to assemble
   *
   * 1) first set up necessary surface arclength derivative quantities for curvature
   *    variable and tension variable.   To do this you will need to compute the derivative
   *    wrt the isoparametric coordinate that aligns with the mesh edge, ascertained from the
   *    bulk element friend on one side of the interface.   
   * 2) Unlike the surface charge type equations above, there is no need to jump to the bulk
   *    and evaluate the normal and tangential stresses, as these will be applied as rotated weak
   *    boundary conditions using a clever manipulation of the liquid momentum residual using
   *    rotations. So you can stay at home in el0 for this routine.
   * 3) Evaluate residual pieces as necessary and add up.   
   */


  /* Declare some neighbor structures */
  int el0 = ei->ielem;
  /* int el1; */
  int nf;

  /* These are needed to get convection velocity and mesh derivatives */

  /* Even though this routine assemble 2 shell equations, we will assume for
   * now that the basis functions are the same 
   */

  /* Unpack variables from structures for local convenience. */
  double h3    = fv->h3;

  /* Initialize d_phi_dxi */
  for (i = 0; i < MDE; i++) {
    d_phi_dxi[i] = 0;
  }

  eqn = R_SHELL_CURVATURE;

  pd0 = pd;  /*set prob description to current shell material */
  
  /* See if there are friends for this element (maximum 2) */
  nf = num_elem_friends[el0];

  /* Initialize temporary arrays */
  memset (res, 0, sizeof(double)*MDE);
  memset (jac, 0, sizeof(double)*MDE*MAX_PROB_VAR*MDE);
  /*
   * Array dh3J_dm is a multiplier to convert the (h3 * det_J)
   * terms of the residuals to their mesh derivatives.
   * This will require mesh equations to be defined on the bulk block!
   */

  /*
   * DOF counts will be needed to construct sensitivities of
   * local shell equations to neighbor bulk element variables.
   * Allocate this array to save them in.
   */
  n_dof = (int *)array_alloc(1, MAX_VARIABLE_TYPES, sizeof(int));

  /* Get neighbor element number(s). */
  if(&elem_friends[el0][0] != NULL) /*Unwetted shells don't have friends*/
    {
      if (nf > 2) EH(-1, "Not set up for more than two element friends!");
    }

  /* 
   * Now that the preliminaries are done, let us compute the necessary building
   * blocks for the structural shells, viz. d(T_sh)/d_xi, d(K_sh)/d_xi, d_phi_d_xi etc.
   */

  for (i = 0; i < ei->dof[eqn]; i++)
    {
      phi_i = bf[eqn]->phi[i];
      d_phi_dxi[i] = bf[eqn]->dphidxi[i][0];
    }

  d_sh_K_dxi = d_sh_tens_dxi = 0.0;
  for (i = 0; i < ei->dof[eqn]; i++)
    {
      d_sh_K_dxi +=    *esp->sh_K[i] * d_phi_dxi[i];
      d_sh_tens_dxi += *esp->sh_tens[i] * d_phi_dxi[i];
    }


  d_sh_x_dxi = d_sh_y_dxi = 0.;
  for (i = 0; i < ei->dof[eqn]; i++)
    {
      node = ei->dof_list[R_MESH1][i];
      index = Proc_Elem_Connect[ei->iconnect_ptr +node];

      d_sh_x_dxi +=    (Coor[0][index] + *esp->d[0][i]) * d_phi_dxi[i];
      d_sh_y_dxi +=    (Coor[1][index] + *esp->d[1][i]) * d_phi_dxi[i];
    }

  det_J_sh = sqrt(d_sh_x_dxi*d_sh_x_dxi + d_sh_y_dxi*d_sh_y_dxi);



  /* First process the side belonging to element el1. */



  /* Assemble the residual equation (Side 1) */
  if (af->Assemble_Residual )
    {
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  phi_i = bf[eqn]->phi[i];

	  /* I only use "diffusion" term to keep things as similar
	   * as possible to the other bulk equation assemble.  Ain't
	   * no real diffusion here  */

	  /*First assemble Normal component */
	  /* PRS: note you need to add the mat properties to these */

	  peqn = upd->ep[R_SHELL_CURVATURE];

	  diffusion  = 0.0;
	  if (pd0->e[eqn])
	    {
	      diffusion = -elc->bend_stiffness * d_phi_dxi[i] * d_sh_K_dxi/det_J_sh 
		- fv->sh_K * fv->sh_tens * phi_i * det_J_sh;
	      diffusion *= pd0->etm[eqn][(LOG2_DIFFUSION)];
	    }

	  /* Sum the terms into res[] */
	  lec->R[peqn][i] += diffusion * wt * h3;

	  /*Now the tangential shell component of momentum */
	  eqn = R_SHELL_TENSION;
	  peqn = upd->ep[eqn];

	  diffusion  = 0.0;
	  if (pd0->e[eqn])
	    {
	      diffusion = phi_i * d_sh_tens_dxi +
		elc->bend_stiffness * phi_i * fv->sh_K * d_sh_K_dxi;
		
	      diffusion *= pd0->etm[eqn][(LOG2_DIFFUSION)];
	    }

	  /* Sum the terms into res[] */
	  lec->R[peqn][i] += diffusion * wt * h3;

	}
    }

  /* Assemble the sensitivity equations (Side 1) */
  if (af->Assemble_Jacobian)
    {
      eqn = R_SHELL_CURVATURE; /* first the normal component */
      peqn = upd->ep[eqn];
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  phi_i = bf[eqn]->phi[i];

	  /* J_sh_K_sh_K: */
          var  = SHELL_CURVATURE;
	  if ( pd->v[var] )
	    {
	      pvar = upd->vp[var];
 
	      /* diffusion term only */
	      diffusion = 0.0;
	      if (pd->e[eqn])
		{
	          for (j = 0; j < ei->dof[var]; j++) 
		    {
		      phi_j = bf[var]->phi[j];

		      diffusion = -elc->bend_stiffness * d_phi_dxi[i] * d_phi_dxi[j]/det_J_sh 
			-phi_j * fv->sh_tens * phi_i * det_J_sh;

		      diffusion *= pd0->etm[eqn][(LOG2_DIFFUSION)];
                      lec->J[peqn][pvar][i][j] += diffusion * wt * h3;
                    }
                }
            }

	  /* J_sh_K_sh_tens: */
          var  = SHELL_TENSION;
	  if ( pd->v[var] )
	    {
	      pvar = upd->vp[var];

	      /* diffusion term only */
	      diffusion = 0.0;
	      if (pd->e[eqn])
		{
	          for (j = 0; j < ei->dof[var]; j++) 
		    {
		      phi_j = bf[var]->phi[j];
		      diffusion = -fv->sh_K * phi_j * phi_i * det_J_sh;

		      diffusion *= pd0->etm[eqn][(LOG2_DIFFUSION)];
                      lec->J[peqn][pvar][i][j] += diffusion * wt * h3;
                    }
                }
            }

	 

	  /* J_k_sh_x:  Side 1 sensitivity */
	  var  = MESH_DISPLACEMENT1;
	  if ( pd0->v[var] )
	    {
	      pvar = upd->vp[var];
	      n_dof[pvar] = ei->dof[var];

	      for(j = 0; j < ei->dof[var]; j++) 
		{
		  d_det_J_dmeshbj = 0.5*(2.*d_sh_x_dxi*d_phi_dxi[j])/det_J_sh;
		  /* Diffusion term */
		  diffusion = 0.0;
		  if (pd0->e[eqn])
		    {
		      diffusion = elc->bend_stiffness * d_phi_dxi[i] * d_sh_K_dxi * d_det_J_dmeshbj / det_J_sh / det_J_sh
			-fv->sh_K * fv->sh_tens * phi_i * d_det_J_dmeshbj;
			 
		      diffusion *= pd0->etm[eqn][(LOG2_DIFFUSION)];
		    }

		  lec->J[peqn][pvar][i][j] += (diffusion) * wt * h3;
		}
	    }

	  /* J_k_sh_y:  Side 1 sensitivity */
	  var  = MESH_DISPLACEMENT2;
	  if ( pd0->v[var] )
	    {
	      pvar = upd->vp[var];
	      n_dof[pvar] = ei->dof[var];

	      for(j = 0; j < ei->dof[var]; j++) 
		{
		  d_det_J_dmeshbj = 0.5*(2.*d_sh_y_dxi*d_phi_dxi[j])/det_J_sh;
		  /* Diffusion term */
		  diffusion = 0.0;
		  if (pd0->e[eqn])
		    {
		      diffusion = elc->bend_stiffness * d_phi_dxi[i] * d_sh_K_dxi * d_det_J_dmeshbj / det_J_sh / det_J_sh
			-fv->sh_K * fv->sh_tens * phi_i * d_det_J_dmeshbj;
			 
		      diffusion *= pd0->etm[eqn][(LOG2_DIFFUSION)];
		    }

		  lec->J[peqn][pvar][i][j] += (diffusion) * wt * h3;
		}
	    }
	}
      
      eqn = R_SHELL_TENSION; /* now the tangential component */
      peqn = upd->ep[eqn];
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  phi_i = bf[eqn]->phi[i];

	  /* J_sh_tens_sh_K: */
          var  = SHELL_CURVATURE;
	  if ( pd->v[var] )
	    {
	      pvar = upd->vp[var];
              n_dof[pvar] = ei->dof[var];

	      /* diffusion term only */
	      diffusion = 0.0;
	      if (pd->e[eqn])
		{
	          for (j = 0; j < ei->dof[var]; j++) 
		    {
		      phi_j = bf[var]->phi[j];

		      diffusion = elc->bend_stiffness * phi_i *
			(phi_j * d_sh_K_dxi +
			 fv->sh_K * d_phi_dxi[j]);

		      diffusion *= pd0->etm[eqn][(LOG2_DIFFUSION)];
                      lec->J[peqn][pvar][i][j] += diffusion * wt * h3;
                    }
                }
            }

	  /* J_sh_tens_sh_tens: */
          var  = SHELL_TENSION;
	  if ( pd->v[var] )
	    {
	      pvar = upd->vp[var];
              n_dof[pvar] = ei->dof[var];

	      /* diffusion term only */
	      diffusion = 0.0;
	      if (pd->e[eqn])
		{
	          for (j = 0; j < ei->dof[var]; j++) 
		    {
		      phi_j = bf[var]->phi[j];
		      diffusion = phi_i * d_phi_dxi[j];

		      diffusion *= pd0->etm[eqn][(LOG2_DIFFUSION)];
                      lec->J[peqn][pvar][i][j] += diffusion * wt * h3;
                    }
                }
            }

	 

	  /* J_k_tens_x:  Side 1 sensitivity */
          /* hMM, THESE ARE CONVENIENTLY ZERO */
	}
    }


  /* Clean up */
  safe_free((void *) n_dof);

  return(0);
} /* assemble_shell_structure() */

/******************************************************************************
 * assemble_shell_tension - Assembles the residual and Jacobian equations
 *                          for inextensible shell structure equation
 *                          SHELL_TENSION when used without SHELL_CURVATURE.
 *                          THIS IS A 2D IMPLEMENTATION ONLY! This means that
 *                          h3 scale factor is 1.0 and we don't bother with
 *                          mesh sensitivities of h3. 
 *
 * Currently, the following equation is implemented:
 *
 *    dT/ds + 0 + nt:T = 0
 *
 * where T is the shell tension and s is the arclength along the surface.
 * This corresponds to the second equation in assemble_shell_structure()
 * with the curvature term omitted.
 *
 *
 * Input
 * =====
 * time_value = The current time.
 * theta      = The implicit-explicit time stepping parameter.
 * delta_t    = The current step size.
 * 
 * Output
 * ======
 * (none)
 * 
 * Returns
 * ======
 * 0  = Success
 * -1 = Failure
 *
 * Revision History
 * ================
 * 22 October 2003 - P. R. Schunk- Creation.
 *
 ******************************************************************************/
int
assemble_shell_tension(double time_value,  /* Time */
		       double theta,       /* Time stepping parameter */
		       double delta_t,     /* Time step size */
                       const double wt,    /* Gauss point weight */
                       double xi[DIM],     /* Local stu coordinates */
                       const Exo_DB *exo)   
{
  int i, j, peqn, var, pvar;
  double phi_i;

  /* note the following definitions restrict this to 1D bar elements */
  double d_phi_dxi[MDE], d_sh_tens_dxi, d_sh_x_dxi, d_sh_y_dxi; 
  double diffusion;
  double res[MDE], jac[MDE][MAX_PROB_VAR][MDE];
  int  dof_map[MDE];
  int *n_dof = NULL;
  int eqn;

  int node, index;
  // double Pi[DIM][DIM];
  //STRESS_DEPENDENCE_STRUCT d_Pi;
  double dTL_dX[DIM][MDE];
  double dTL_dv[DIM][MDE];
  double dTL_dP[MDE];


  /* After reviewing the methodology presented in the original shell equations
   * addressing surface charge, and considering that the structural shell equations
   * are more complicated in the sense that they deal with spatial derivatives in 
   * the surface, we will take the following approach to assemble
   *
   * 1) first set up necessary surface arclength derivative quantities for curvature
   *    variable and tension variable.   To do this you will need to compute the derivative
   *    wrt the isoparametric coordinate that aligns with the mesh edge, ascertained from the
   *    bulk element friend on one side of the interface.   
   * 2) Unlike the surface charge type equations above, there is no need to jump to the bulk
   *    and evaluate the normal and tangential stresses, as these will be applied as rotated weak
   *    boundary conditions using a clever manipulation of the liquid momentum residual using
   *    rotations. So you can stay at home in el0 for this routine.
   * 3) Evaluate residual pieces as necessary and add up.   
   */


  /* Declare some neighbor structures */
#if 0
  int el0 = ei->ielem;
  int nf;
  /* See if there are friends for this element (maximum 2) */
  nf = num_elem_friends[el0];
#endif

  /* These are needed to get convection velocity and mesh derivatives */

  /* Even though this routine assemble 2 shell equations, we will assume for
   * now that the basis functions are the same 
   */

  /* Unpack variables from structures for local convenience. */
  double h3    = fv->h3;

  eqn = R_SHELL_TENSION;

  /* Initialize temporary arrays */
  memset(res, 0, sizeof(double)*MDE);
  memset(dof_map, 0, sizeof(int)*MDE);
  memset(jac, 0, sizeof(double)*MDE*MAX_PROB_VAR*MDE);

  /*
   * Array dh3J_dm is a multiplier to convert the (h3 * det_J)
   * terms of the residuals to their mesh derivatives.
   * This will require mesh equations to be defined on the bulk block!
   */

  /*
   * DOF counts will be needed to construct sensitivities of
   * local shell equations to neighbor bulk element variables.
   * Allocate this array to save them in.
   */
  n_dof = (int *)array_alloc(1, MAX_VARIABLE_TYPES, sizeof(int));

#if 0 /* if setting 1 uncomment el1 declaration */
  /* Get neighbor element number(s). */
  if(&elem_friends[el0][0] != NULL) /*Unwetted shells don't have friends*/
    {
      el1 = elem_friends[el0][0];
      if (nf > 2) EH(-1, "Not set up for more than two element friends!");
    }
	
  if (nf == 1 ) /* Load up data from bulk neighbor element in order to compute tangent loading */
    {
      err = load_neighbor_var_data(el0, el1, n_dof, dof_map,
                                   n_dofptr, id_side, xi, exo);
	
      fluid_stress( Pi, &d_Pi );
    }
#endif

  /* 
   * Now that the preliminaries are done, let us compute the necessary building
   * blocks for the structural shells, viz. d(T_sh)/d_xi, d(K_sh)/d_xi, d_phi_d_xi etc.
   */

  for (i = 0; i < ei->dof[eqn]; i++)
    {
      d_phi_dxi[i] = bf[eqn]->dphidxi[i][0];
    }
	
  d_sh_tens_dxi = 0.0;
  for (i = 0; i < ei->dof[eqn]; i++)
    {
      d_sh_tens_dxi += *esp->sh_tens[i] * d_phi_dxi[i];
    }


  
  d_sh_x_dxi = d_sh_y_dxi = 0.;
  for (i = 0; i < ei->dof[eqn]; i++)
    {
      node = ei->dof_list[eqn][i];
      index = Proc_Elem_Connect[ei->iconnect_ptr +node];
	  
      if( n_dof[R_MESH1] > 0  ) {
	d_sh_x_dxi +=    (Coor[0][index] + *esp->d[0][ dof_map[i] ]) * d_phi_dxi[i];
	d_sh_y_dxi +=    (Coor[1][index] + *esp->d[1][ dof_map[i] ]) * d_phi_dxi[i];
      }
      else {
	d_sh_x_dxi +=    (Coor[0][index] ) * d_phi_dxi[i];
	d_sh_y_dxi +=    (Coor[1][index] ) * d_phi_dxi[i];
      }
    }



  /* Compute tangent loading and appropriate sensitivities */

  memset( dTL_dv, 0, sizeof(double)*DIM*MDE);
  memset( dTL_dX, 0, sizeof(double)*DIM*MDE);
  memset( dTL_dP, 0, sizeof(double)*MDE);




  /* First process the side belonging to element el1. */

  peqn = upd->ep[eqn];


  /* Assemble the residual equation (Side 1) */
  if ( af->Assemble_Residual )
    {
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  phi_i = bf[eqn]->phi[i];


	  /* The tangential shell component of momentum */

	  diffusion  = 0.0;
	  if (pd->e[eqn])
	    {
	      diffusion = phi_i * d_sh_tens_dxi;
		
	      diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
	    }

	  /* Sum the terms into res[] */
	  lec->R[peqn][i] += diffusion * wt * h3;

	}
    }

  /* Assemble the sensitivity equations (Side 1) */
  if ( af->Assemble_Jacobian )
    {
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  phi_i = bf[eqn]->phi[i];
		  
	  /* J_sh_tens_sh_tens: */
          var  = SHELL_TENSION;
	  if ( pd->v[var] )
	    {
	      pvar = upd->vp[var];
			  
	      /* diffusion term only */
	      diffusion = 0.0;
			  
	      for (j = 0; j < ei->dof[var]; j++) 
		{

		  diffusion = phi_i * d_phi_dxi[j];
				  
		  diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
		  lec->J[peqn][pvar][i][j] += diffusion * wt * h3;
		}
			  
	    }
		  


	}
    }


  /* Clean up */
  safe_free((void *) n_dof);

  return(0);
} /* assemble_shell_tension() */

/******************************************************************************
 * assemble_shell_coordinates - Assembles the residual and Jacobian equations
 *                             for inextensible shell coordinate equations
 *                            THIS IS A 2D IMPLEMENTATION ONLY!
 *                             h3 scale factor is 1.0 and we don't bother with
 *                             mesh sensitivities of h3.
 *
 * Input
 * =====
 * time_value = The current time.
 * theta      = The implicit-explicit time stepping parameter.
 * delta_t    = The current step size.
 * 
 * Output
 * ======
 * (none)
 * 
 * Returns
 * ======
 * 0  = Success
 * -1 = Failure
 *
 * Revision History
 * ================
 * 10 November 2003 - P. R. Schunk- Creation.
 *
 ******************************************************************************/
int
assemble_shell_coordinates(double time_value,  /* Time */
			   double theta,       /* Time stepping parameter */
			   double delta_t,     /* Time step size */
			   const double wt,    /* Gauss point weight */
			   double xi[DIM],     /* Local stu coordinates */
			   const Exo_DB *exo)
{
  int i, j, peqn, var, pvar;
  double phi_i, phi_j;

  /* note the following definitions restrict this to 1D bar elements */
  double d_phi_dxi[MDE], d_sh_x_dxi, d_sh_y_dxi; 
  double diffusion;
  int eqn;
  double d_det_J_dmeshbj, det_J_sh;
 
  PROBLEM_DESCRIPTION_STRUCT *pd0;
  int node, index;

  /* Declare some neighbor structures needed to look across from the shells to see
   * if bulk elements are attached, and whether there are mesh equations. If this is
   * true, we do not apply the arclength equation*/
  int nf;
  int el0 = ei->ielem;
  double h3    = fv->h3;

  /* Some prelims to see if we are wetted or unwetted */

  /* See if there is a friend for this element */
  /* In this case we assume if there is a friend, there are mesh equations
   * in that friend's prob description, and hence we've already a mesh2 equation
   * and don't need an arclength equation, so don't apply below! For structural
   * shells a mesh equation in the bulk is a sure bet!
   */
  nf = num_elem_friends[el0]; 



  /* 
   * Please see comments for assemble_shell_structure for editorial comments.  
   * This routine assembles the equations that define structural shell 
   * coordinates given the curvature as a function of arclength, viz.
   *
   * d_2(x)/d_s_2 + K*dy/ds = 0   ;   d_2_(y)/d_s_2 - K*dx/ds = 0
   *
   * Note that relating these coordinates to the actual displacements is 
   * accomplished throught the boundary conditions applied to mesh1 and mesh2 
   * equations, viz.
   *  dx = x - X and dy = y - Y   
   * where X, Y are the mesh coordinates and dx and dy are the 
   * displacment varialbes. 
   */


  eqn = R_MESH1;

  pd0 = pd;  /*set prob description to current shell material */
  
  for (i = 0; i < ei->dof[eqn]; i++)
    {
      d_phi_dxi[i] = bf[eqn]->dphidxi[i][0];
    }

  d_sh_x_dxi = d_sh_y_dxi = 0.;
  for (i = 0; i < ei->dof[eqn]; i++)
    {
      node = ei->dof_list[R_MESH1][i];
      index = Proc_Elem_Connect[ei->iconnect_ptr +node];

      d_sh_x_dxi +=    (Coor[0][index] + *esp->d[0][i]) * d_phi_dxi[i];
      d_sh_y_dxi +=    (Coor[1][index] + *esp->d[1][i]) * d_phi_dxi[i];
    }

  det_J_sh = sqrt(d_sh_x_dxi*d_sh_x_dxi + d_sh_y_dxi*d_sh_y_dxi);


  /* First process the side belonging to element el1. */



  /* Assemble the residual equation (Side 1) */
  if ( af->Assemble_Residual )
    {
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  phi_i = bf[eqn]->phi[i];

	  /* I only use "diffusion" term to keep things as similar
	   * as possible to the other bulk equation assemble.  Ain't
	   * no real diffusion here  */

	  /*First assemble Normal component */
	  /* PRS: note you need to add the mat properties to these */

	  peqn = upd->ep[R_MESH1];

	  diffusion  = 0.0;
	  if (pd0->e[eqn] && nf == 0)
	    {
	      /* diffusion =  (d_sh_x_dxi * d_phi_dxi[i])/det_J_sh
		 -phi_i * fv->sh_K * d_sh_y_dxi; */

	      /* Const Node Space version */
	      diffusion = -det_J_sh * d_phi_dxi[i]; 

	      diffusion *= pd0->etm[eqn][(LOG2_DIFFUSION)];
	    }

	  /* Sum the terms into res[] */
	  lec->R[peqn][i] += diffusion * wt * h3;

	  /*Now the tangential shell component of momentum */
	  eqn = R_MESH2;
	  peqn = upd->ep[eqn];

	  diffusion  = 0.0;
	  if (pd0->e[eqn])
	    {
	      diffusion = -(d_sh_y_dxi * d_phi_dxi[i])/det_J_sh
		- phi_i * fv->sh_K * d_sh_x_dxi;
	      
	      diffusion *= BIG_PENALTY*pd0->etm[eqn][(LOG2_DIFFUSION)];
	    }

	  /* Sum the terms into res[] */
	  lec->R[peqn][i] += diffusion * wt * h3;
	    
	}
    }

  /* Assemble the sensitivity equations (Side 1) */
  if (af->Assemble_Jacobian)
    {
      eqn = R_MESH1; /* first the X equation */
      peqn = upd->ep[eqn];
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  phi_i = bf[eqn]->phi[i];

	  /* J_sh_x_sh_K: */
          var  = SHELL_CURVATURE;
	  if (pd->v[var] && nf == 0)
	    {
	      pvar = upd->vp[var];

	      /* diffusion term only */
	      diffusion = 0.0;
	      if (pd->e[eqn])
		{
	          for (j = 0; j < ei->dof[var]; j++) 
		    {
		      phi_j = bf[var]->phi[j];

		      diffusion = -phi_i * phi_j * d_sh_y_dxi; 
		      diffusion = 0.;  /*const node space version */

		      diffusion *= pd0->etm[eqn][(LOG2_DIFFUSION)];
                      lec->J[peqn][pvar][i][j] += diffusion * wt * h3;
                    }
                }
            }

	  /* J_sh_x_mesh_x:  Side 1 sensitivity */
	  var  = MESH_DISPLACEMENT1;
	  if (pd0->v[var] && nf == 0)
	    {
	      pvar = upd->vp[var];
 
	      for (j = 0; j < ei->dof[var]; j++) 
		{
		  d_det_J_dmeshbj = 0.5*(2.*d_sh_x_dxi*d_phi_dxi[j])/det_J_sh;
		  /* if(d_det_J_dmeshbj <= 1.e-10) d_det_J_dmeshbj = 1.e-20; */
		  /* Diffusion term */
		  diffusion = 0.0;
		  if (pd0->e[eqn])
		    {
		      /*diffusion = -(d_sh_x_dxi * d_phi_dxi[i]) * d_det_J_dmeshbj / det_J_sh / det_J_sh
			+ (d_phi_dxi[j] * d_phi_dxi[i])/det_J_sh;*/
		      diffusion = -d_det_J_dmeshbj * d_phi_dxi[i];
			 
		      diffusion *= pd0->etm[eqn][(LOG2_DIFFUSION)];
		    }

		  lec->J[peqn][pvar][i][j] += (diffusion) * wt * h3;
		}
	    }
	  var  = MESH_DISPLACEMENT2;
	  if (pd0->v[var] && nf == 0)
	    {
	      pvar = upd->vp[var];
 
	      for (j = 0; j < ei->dof[var]; j++) 
		{
		  d_det_J_dmeshbj = 0.5*(2.*d_sh_y_dxi*d_phi_dxi[j])/det_J_sh;
		  /* if(d_det_J_dmeshbj <= 1.e-10) d_det_J_dmeshbj = 1.e-20;*/
		  /* Diffusion term */
		  diffusion = 0.0;
		  if (pd0->e[eqn])
		    {
		      /*diffusion = -(d_sh_x_dxi * d_phi_dxi[i]) * d_det_J_dmeshbj / det_J_sh / det_J_sh
			-phi_i * fv->sh_K * d_phi_dxi[j];*/
		      diffusion = -d_det_J_dmeshbj * d_phi_dxi[i];
			 
		      diffusion *= pd0->etm[eqn][(LOG2_DIFFUSION)];
		    }

		  lec->J[peqn][pvar][i][j] += (diffusion) * wt * h3;
		}
	    }
  	}
      
      eqn = R_MESH2; /* now the tangential component */
      peqn = upd->ep[eqn];

      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  phi_i = bf[eqn]->phi[i];

	  /* J_sh_y_sh_K: */
          var  = SHELL_CURVATURE;
	  if (pd->v[var])
	    {
	      pvar = upd->vp[var];

	      /* diffusion term only */
	      diffusion = 0.0;
	      if (pd->e[eqn])
		{
	          for (j = 0; j < ei->dof[var]; j++) 
		    {
		      phi_j = bf[var]->phi[j];

		      diffusion = -phi_i * phi_j * d_sh_x_dxi;

		      diffusion *= BIG_PENALTY*pd0->etm[eqn][(LOG2_DIFFUSION)];
                      lec->J[peqn][pvar][i][j] += diffusion * wt * h3;
                    }
                }
            }

	  /* J_sh_y_mesh_x:  Side 1 sensitivity */
 	  var  = MESH_DISPLACEMENT1;
	  if (pd0->v[var])
	    {
	      pvar = upd->vp[var];
 
	      for (j = 0; j < ei->dof[var]; j++) 
		{
		  d_det_J_dmeshbj = 0.5*(2.*d_sh_x_dxi*d_phi_dxi[j])/det_J_sh;
		  /*if(d_det_J_dmeshbj <= 1.e-10) d_det_J_dmeshbj = 1.e-20;*/

		  /* Diffusion term */
		  diffusion = 0.0;
		  if (pd0->e[eqn])
		    {
		      diffusion = (d_sh_y_dxi * d_phi_dxi[i]) * d_det_J_dmeshbj / det_J_sh / det_J_sh
			-(phi_i * fv->sh_K * d_phi_dxi[j]);
			 
		      diffusion *= BIG_PENALTY* pd0->etm[eqn][(LOG2_DIFFUSION)];
		    }

		  lec->J[peqn][pvar][i][j] += (diffusion) * wt * h3;
		}
	    }
	  var  = MESH_DISPLACEMENT2;
	  if (pd0->v[var])
	    {
	      pvar = upd->vp[var];
 
	      for (j = 0; j < ei->dof[var]; j++) 
		{
		  d_det_J_dmeshbj = 0.5*(2.*d_sh_y_dxi*d_phi_dxi[j])/det_J_sh;
		  /*if(d_det_J_dmeshbj <= 1.e-10) d_det_J_dmeshbj = 1.e-20; */
		  /* Diffusion term */
		  diffusion = 0.0;
		  if (pd0->e[eqn])
		    {
		      diffusion = (d_sh_y_dxi * d_phi_dxi[i]) * d_det_J_dmeshbj / det_J_sh / det_J_sh
			-(d_phi_dxi[j] * d_phi_dxi[i])/det_J_sh; 
			 
		      diffusion *= BIG_PENALTY*pd0->etm[eqn][(LOG2_DIFFUSION)];
		    }

		  lec->J[peqn][pvar][i][j] += (diffusion) * wt * h3;
		}
	    }
	}
    }


  /* Clean up */

  return(0);
} /* assemble_shell_coordinates() */

/******************************************************************************
 * assemble_shell_diffusion - Assembles the residual and Jacobian terms
 *                            for the inextensible shell diffusion equation:
 *                            SHELL_DIFF_FLUX
 *                            THIS IS A 2D IMPLEMENTATION ONLY! This means that
 *                            h3 scale factor is 1.0 and we don't bother with
 *                            mesh sensitivities of h3 (for now).
 *
 * Currently, the following equation is implemented:
 *
 *    Js + Ds * Gs * Va * K * (dJs/ds) = 0
 *
 * where:
 *      Js = surface diffusive flux unknown
 *      K  = surface curvature unknown
 *      Ds = surface diffusion coefficient (temperature-dependent)
 *      Gs = surface energy
 *      Av = atomic volume
 *       s = arclength along the surface.
 *
 *
 * Input
 * =====
 * time_value = The current time.
 * theta      = The implicit-explicit time stepping parameter.
 * delta_t    = The current step size.
 *
 * Output
 * ======
 * (none)
 *
 * Returns
 * ======
 * 0  = Success
 * -1 = Failure
 *
 * Revision History
 * ================
 * 28 June 2005 - E. D. Wilkes- Creation.
 *
 ******************************************************************************/
int
assemble_shell_diffusion(double time_value,  /* Time */
                         double theta,       /* Time stepping parameter */
                         double delta_t,     /* Time step size */
                         const double wt,    /* Gauss point weight */
                         double xi[DIM],     /* Local stu coordinates */
                         const Exo_DB *exo)
{
  int err, i, j, p;
  int peqn = -1;
  int var, pvar;
  double phi_i, phi_j;

  /* note the following definitions restrict this to 1D bar elements */
  double d_phi_dxi[MDE];
  double diffusion;
  double Ds, Gs, Va;
  double res[MDE], jac[MDE][MAX_PROB_VAR][MDE];
  int *n_dof = NULL;
  int n_dofptr[MAX_VARIABLE_TYPES][MDE];
  int eqn;
  double  det_J_sh;

  PROBLEM_DESCRIPTION_STRUCT *pd0;

  /* After reviewing the methodology presented in the original shell equations
   * addressing surface charge, and considering that the structural shell equations
   * are more complicated in the sense that they deal with spatial derivatives in
   * the surface, we will take the following approach to assemble
   *
   * 1) first set up necessary surface arclength derivative quantities for curvature
   *    variable and tension variable.   To do this you will need to compute the derivative
   *    wrt the isoparametric coordinate that aligns with the mesh edge, ascertained from the
   *    bulk element friend on one side of the interface.
   * 2) Unlike the surface charge type equations above, there is no need to jump to the bulk
   *    and evaluate the normal and tangential stresses, as these will be applied
   as rotated weak
   *    boundary conditions using a clever manipulation of the liquid momentum residual using
   *    rotations. So you can stay at home in el0 for this routine.
   * 3) Evaluate residual pieces as necessary and add up.
   */

  /* Declare some neighbor structures */
  int el0 = ei->ielem;
  int el1 = -1;
  int nf;

  /* These are needed to get convection velocity and mesh derivatives */

  /* Even though this routine assemble 2 shell equations, we will assume for
   * now that the basis functions are the same
   */

  /* Unpack variables from structures for local convenience. */
  double h3 = fv->h3;

  eqn = R_SHELL_DIFF_FLUX;

  pd0 = pd;  /*set prob description to current shell material */

  /* See if there are friends for this element (maximum 2) */
  nf = num_elem_friends[el0];

  /* Initialize temporary arrays */
  memset(res, 0, sizeof(double)*MDE);
  memset(jac, 0, sizeof(double)*MDE*MAX_PROB_VAR*MDE);
  /*
   * Array dh3J_dm is a multiplier to convert the (h3 * det_J)
   * terms of the residuals to their mesh derivatives.
   * This will require mesh equations to be defined on the bulk block!
   */

  /*
   * DOF counts will be needed to construct sensitivities of
   * local shell equations to neighbor bulk element variables.
   * Allocate this array to save them in.
   */
  n_dof = (int *)array_alloc(1, MAX_VARIABLE_TYPES, sizeof(int));

  /* Get neighbor element number(s). */
  if (&elem_friends[el0][0] != NULL) /*Unwetted shells don't have friends*/
    {
      el1 = elem_friends[el0][0];
      if (nf > 2) EH(-1, "Not set up for more than two element friends!");
    }

  /* Load bulk element data */
  err = load_neighbor_var_data(el0, el1, n_dof, NULL, n_dofptr, -1, xi, exo);
  EH(err, "Problem loading bulk element data in assemble_shell_diffusion!");

  /*
   * Now that the preliminaries are done, let us compute the necessary building
   * blocks for the shells, viz. d(sh_J)/d_xi, d(sh_Kd)/d_xi, d_phi_d_xi etc.
   */

  for (i = 0; i < ei->dof[eqn]; i++)
    {
      phi_i = bf[eqn]->phi[i];
      d_phi_dxi[i] = bf[eqn]->dphidxi[i][0];
    }

  det_J_sh = fv->sdet;

  /* Load surface material properties */
  /* err = load_shell_diffusivity(&Ds);     [EDW: Not ready yet!]  */
  Ds = mp->cur_diffusivity[0];  /* EDW: using curvature diffusivity for now! */
  Gs = mp->surface_tension;
  Va = mp->density;
  Ds = 1.0;

  /* Assemble the residual equations  */
  if (af->Assemble_Residual)
    {
      for (i = 0; i < ei->dof[eqn]; i++)
        {
          phi_i = bf[eqn]->phi[i];

          /* I only use "diffusion" term to keep things as similar
           * as possible to the other bulk equation assemble.  Ain't
           * no real diffusion here  */

          peqn = upd->ep[eqn];

          diffusion  = 0.0;
          if (pd0->e[eqn])
            {
              diffusion = fv->sh_J * phi_i * det_J_sh
		+ Ds * Gs * Va * fv->sh_Kd * d_phi_dxi[i];
              diffusion *= pd0->etm[eqn][(LOG2_DIFFUSION)];
            }

          /* Sum the terms into res[] */
          lec->R[peqn][i] += diffusion * wt * h3;
        }
    }

  /* Assemble the sensitivity equations (Side 1) */
  if (af->Assemble_Jacobian)
    {
      for (i = 0; i < ei->dof[eqn]; i++)
        {
          phi_i = bf[eqn]->phi[i];

          /* J_sh_J_sh_Kd: */
          var  = SHELL_DIFF_CURVATURE;
          if (pd->v[var])
            {
              pvar = upd->vp[var];

              /* diffusion term only */
              diffusion = 0.0;
              for (j = 0; j < ei->dof[var]; j++)
                {
                  phi_j = bf[var]->phi[j];

                  diffusion = Ds * Gs * Va * phi_j * d_phi_dxi[i];
                  diffusion *= pd0->etm[eqn][(LOG2_DIFFUSION)];
                  lec->J[peqn][pvar][i][j] += diffusion * wt * h3;
                }
            }

          /* J_sh_J_sh_J: */
          var  = SHELL_DIFF_FLUX;
          if (pd->v[var])
            {
              pvar = upd->vp[var];

              /* diffusion term only */
              diffusion = 0.0;
              for (j = 0; j < ei->dof[var]; j++)
                {
                  phi_j = bf[var]->phi[j];

                  diffusion = phi_j * phi_i * det_J_sh;
                  diffusion *= pd0->etm[eqn][(LOG2_DIFFUSION)];
                  lec->J[peqn][pvar][i][j] += diffusion * wt * h3;
                }
            }

          /* J_k_sh_x:  Side 1 sensitivity */
          for (p = 0; p < pd->Num_Dim; p++)
            {
              var  = MESH_DISPLACEMENT1 + p;
              pvar = upd->vp[var];
              if (pvar > -1)
                {

                  for(j = 0; j < n_dof[var]; j++)
                    {
                      /* Diffusion term */
                      diffusion = fv->sh_J * phi_i * fv->dsurfdet_dx[p][j];
                      diffusion *= pd0->etm[eqn][(LOG2_DIFFUSION)];
                      lec->J[peqn][pvar][i][j] += diffusion * wt * h3;
                    }
                }
            }
        }
    }

  /* Clean up */
  safe_free((void *) n_dof);

  return(0);
} /* assemble_shell_diffusion() */

/******************************************************************************
 * assemble_shell_geometry -  Assembles the residual and Jacobian equations
 *                            for inextensible shell geometry equations:
 *                            SHELL_NORMAL1(2) and SHELL_DIFF_CURVATURE
 *                            THIS IS A 2D IMPLEMENTATION ONLY! This means that
 *                            h3 scale factor is 1.0 and we don't bother with
 *                            mesh sensitivities of h3 (for now).
 *                   HKM -> Ok big error, as it isn't good with cylindrical geometry
 *
 * Currently, the following equations are implemented.
 *
 *    N - fv->snormal = 0 (2 components)
 *    K - div_s_N = 0
 *
 * where:
 *       N = surface normal vector unknown
 *       K = surface curvature unknown
 *       div_s_N = (I - NN) : grad_N (surface divergence defn)
 *
 * Note: Normal vector equations are used to avoid second derivative
 *       terms in the curvature equation.
 *
 * Input
 * =====
 * time_value = The current time.
 * theta      = The implicit-explicit time stepping parameter.
 * delta_t    = The current step size.
 *
 * Output
 * ======
 * (none)
 *
 * Returns
 * ======
 * 0  = Success
 * -1 = Failure
 *
 * Revision History
 * ================
 * 08 June 2005 - E. D. Wilkes- Creation.
 *
 ******************************************************************************/
int
assemble_shell_geometry(double time_value,  /* Time */
                        double theta,       /* Time stepping parameter */
                        double delta_t,     /* Time step size */
                        const double wt,    /* Gauss point weight */
                        double xi[DIM],     /* Local stu coordinates */
                        const Exo_DB *exo)
{
  int err, i, j, jk, p, b, peqn, var, pvar;
  double phi_i, phi_j;

  /* note the following definitions restrict this to 1D bar elements */
  double diffusion;
  double div_s_nv, d_div_s_nv_dnv[DIM][MDE], d_div_s_nv_dmesh[DIM][MDE];
  double res[MDE], jac[MDE][MAX_PROB_VAR][MDE];
  int *n_dof = NULL;
  int n_dofptr[MAX_VARIABLE_TYPES][MDE];
  int eqn;

  /*
   * Find the variable whose basis functions are used to do element integrations of shape
   * functions for the shell equation
   */
  int shapeVar = pd->ShapeVar;
  double det_J = bf[shapeVar]->detJ;

  double *p_div_s_nv;
  PROBLEM_DESCRIPTION_STRUCT *pd0;
  BASIS_FUNCTIONS_STRUCT *bfn;
  int *dof_map;

  /* Declare some neighbor structures */
  int el0 = ei->ielem;
  int el1 = -1;
  int nf;
#ifdef DEBUG_HKM
  struct Element_Indices *ei_ptr;
#endif
  /*
   * Even though this routine assemble 3 shell equations, we will assume for
   * now that the basis functions are the same
   */

  /* Unpack variables from structures for local convenience. */
  double h3    = fv->h3;

  pd0 = pd;  /*set prob description to current shell material */
  p_div_s_nv = &div_s_nv;  /* Surface divergence pointer */

  /* See if there are friends for this element (maximum 2) */
  nf = num_elem_friends[el0];

  /* Initialize temporary arrays */
  dof_map = (int *)array_alloc(1, MDE, sizeof(int));
  memset (res, 0, sizeof(double)*MDE);
  memset (jac, 0, sizeof(double)*MDE*MAX_PROB_VAR*MDE);
  memset (d_div_s_nv_dmesh, 0, sizeof(double)*DIM*MDE);
  memset (d_div_s_nv_dnv, 0, sizeof(double)*DIM*MDE);

  /*
   * Array dh3J_dm is a multiplier to convert the (h3 * det_J)
   * terms of the residuals to their mesh derivatives.
   * This will require mesh equations to be defined on the bulk block!
   */

  /*
   * DOF counts will be needed to construct sensitivities of
   * local shell equations to neighbor bulk element variables.
   * Allocate this array to save them in.
   */
  n_dof = (int *)array_alloc(1, MAX_VARIABLE_TYPES, sizeof(int));

  /* Get neighbor element number(s). */
  if(&elem_friends[el0][0] != NULL) /*Unwetted shells don't have friends*/
    {
      el1 = elem_friends[el0][0];
      if (nf > 2) EH(-1, "Not set up for more than two element friends!");
    }

  /* Load bulk element data */
  err = load_neighbor_var_data(el0, el1, n_dof, dof_map, n_dofptr, -1, xi, exo);
  EH(err, "Problem loading bulk element data in assemble_shell_diffusion!");

  /* Load normal vector surface divergence and sensitivities */
  err = shell_normal_div_s(p_div_s_nv, d_div_s_nv_dnv, d_div_s_nv_dmesh);
  eqn = R_SHELL_DIFF_CURVATURE;
  bfn = bf[eqn];

  /*
   * Now that the preliminaries are done, let us compute the necessary building
   * blocks for the shells, viz.
   */
  for (i = 0; i < ei->dof[eqn]; i++)
    {
      phi_i = bf[eqn]->phi[i];
    }

  // Curvature H = -0.5 (del_s dot n) -> so we add in a -0.5 fac
  double fac = -0.5;

  /* Assemble the residual equations  */
  if (af->Assemble_Residual)
    {
      /* Assemble curvature equation first */
      eqn = R_SHELL_DIFF_CURVATURE;
      if (pd->e[eqn])
        {
          bfn = bf[eqn];
          peqn = upd->ep[R_SHELL_DIFF_CURVATURE];

          diffusion  = 0.0;
          for (i = 0; i < ei->dof[eqn]; i++)
            {
              phi_i = bfn->phi[i];
              diffusion = (fv->sh_Kd - fac * div_s_nv) * phi_i * det_J;
              diffusion *= pd0->etm[eqn][(LOG2_DIFFUSION)];

              /* Sum the terms into res[] */
              lec->R[peqn][i] += diffusion * wt * h3;
            }
        }
#ifdef DEBUG_HKM
      /*
       *   Here we assert a relationship that is important. the surface determinant, which is
       *   calculated from the surface integral of the bulk parent element, is equal to the
       *   determinant of the shell jacobian * h3 of the shell element. 
       *   
       */
      if (fabs(det_J * h3 - fv->sdet) > 1.0E-9) {
	printf("we have a problem: assertion on sdet\n");
        exit(-1);
      }
#endif
      /* Now the shell normal vector (2 components) */
      for (p = 0; p < pd->Num_Dim; p++)
        {
          eqn = R_SHELL_NORMAL1 + p;
          if (pd0->e[eqn])
            {
              bfn = bf[eqn];
              peqn = upd->ep[eqn];
              for (i = 0; i < ei->dof[eqn]; i++)
                {
                  phi_i = bfn->phi[i];
                  diffusion = (fv->n[p] - fv->snormal[p]) * phi_i * det_J;
                  diffusion *= pd0->etm[eqn][(LOG2_DIFFUSION)];

                  /* Sum the terms into res[] */
                  lec->R[peqn][i] += diffusion * wt * h3;
                }
            }
        }
    }

  /* Assemble the sensitivity equations (Side 1) */
  if (af->Assemble_Jacobian)
    {
      eqn = R_SHELL_DIFF_CURVATURE;
      if (pd->e[eqn])
        {
          peqn = upd->ep[eqn];
          bfn = bf[eqn];

          for (i = 0; i < ei->dof[eqn]; i++)
            {
              phi_i = bfn->phi[i];

              /* J_sh_Kd_sh_Kd: */
              var  = SHELL_DIFF_CURVATURE;
              if (pd->v[var])
                {
                  pvar = upd->vp[var];

                  /* diffusion term only */
                  diffusion = 0.0;
                  for (j = 0; j < ei->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];
                      diffusion = phi_i * phi_j * det_J;
                      diffusion *= pd0->etm[eqn][(LOG2_DIFFUSION)];
                      lec->J[peqn][pvar][i][j] += diffusion * wt * h3;
                    }
                }

              for (p = 0; p < pd->Num_Dim; p++)
                {
                  /* J_sh_Kd_n: */
                  var = SHELL_NORMAL1 + p;
                  if (pd->v[var])
                    {
                      pvar = upd->vp[var];
                      for (j = 0; j < ei->dof[var]; j++)
                        {
                          diffusion = -d_div_s_nv_dnv[p][j] * phi_i * det_J;
                          diffusion *= pd0->etm[eqn][(LOG2_DIFFUSION)];
                          lec->J[peqn][pvar][i][j] += diffusion * wt * h3;
                        }
                    }

                  /* J_sh_Kd_d: */
		  /*
		   *  HKM Before this can be updated to the new method, I have to 
		   *      update d_div_s_nv_dmesh[][]. This calculation is done
		   *      in the routine shell_normal_div_s. The whole routine
		   *      needs to be changed from shell mesh basis functions to
		   *      bulk basis functions. Thus, this step has to be done
		   *      in a fairly systematic fashion. 1/5/2009
		   */
		  var  = MESH_DISPLACEMENT1 + p;
                  pvar = upd->vp[var];
                  if (pvar > -1)
                    {
                      for (j = 0; j < ei->dof[var]; j++)
                        {
                          jk = dof_map[j];
                          phi_j = bf[var]->phi[j];
                          diffusion = ( (fv->sh_Kd - fac * div_s_nv) * phi_i
                                        * fv->dsurfdet_dx[p][jk]
                                        - fac * d_div_s_nv_dmesh[p][j] * phi_i * det_J);
                          diffusion *= pd0->etm[eqn][(LOG2_DIFFUSION)];

                          lec->J[peqn][pvar][i][j] += diffusion * wt * h3;
                        }
                    }
                }
            }
        }

      for (p = 0; p < pd->Num_Dim; p++)
        {
          eqn = R_SHELL_NORMAL1 + p;
          if (pd->e[eqn])
            {
              peqn = upd->ep[eqn];
              bfn = bf[eqn];

              for (i = 0; i < ei->dof[eqn]; i++)
                {
                  phi_i = bfn->phi[i];

                  /* J_n_n: */
                  var = eqn;
                  pvar = upd->vp[var];

                  for (j = 0; j < ei->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];

                      diffusion = phi_i * phi_j * det_J;
                      diffusion *= pd0->etm[eqn][(LOG2_DIFFUSION)];
                      lec->J[peqn][pvar][i][j] += diffusion * wt * h3;
                    }

                  /* J_n_d */
                  for (b = 0; b < pd->Num_Dim; b++)
                    {
                      var = MESH_DISPLACEMENT1 + b;
#ifdef DEBUG_HKM
		      ei_ptr = ((ei->owningElement_ei_ptr[var]) ? (ei->owningElement_ei_ptr[var]) : ei);
		      if (n_dof[var] != ei_ptr->dof[var]) {
			printf("found a logic error\n");
			exit(-1);
		      }
#endif
		      /*
		       * Here, the unknown is defined on the remote element.
		       * So get the DOF count from n_dof.
		       */
		      if (n_dof[var] > 0) 
			{
			  pvar = upd->vp[var];
			  
			  for (j = 0; j < n_dof[var]; j++) 
			    {			   
			      lec->J[peqn][pvar][i][j] += -fv->dsnormal_dx[p][b][j] * 
				phi_i * pd0->etm[eqn][(LOG2_DIFFUSION)] * wt * det_J * h3;
			      
			      lec->J[peqn][pvar][i][j] += (fv->n[p] - fv->snormal[p]) * 
				phi_i * pd0->etm[eqn][(LOG2_DIFFUSION)] * wt * fv->dsurfdet_dx[b][j];
			    }
			}
                    }
                }
            }
        }
    }

  /* Clean up */
  safe_free((void *) n_dof);

  return(0);
} /* assemble_shell_geometry() */

/*
 * Here's a RECIPE for adding new boundary conditions so you don't have any
 * excuses not to add new ones.  The changes should be made in at least
 * four files (rf_bc_const.h, mm_names.h, mm_input.c, and bc_[method].c)
 * for some boundary conditions you may want to make additions elsewhere also.
 * One example of extra additions is in el_exoII_io.c where the ss_dup_list
 * is created - you may want to adapt the logic for rotation of new conditions.
 * (note that these lines are repeated at each place where you need to 
 * make changes):
 *  Boundary Condition  Step 1: add Macro Substitution for your new boundary
 *                              condition in rf_bc_const.h - this is how you
 *      rf_bc_const.h           will refer to the new boundary condition 
 *                              throughout the code.  Make sure the integer
 *                              you choose is unique from all the other BC
 *                              types.
 *  Boundary Condition  Step 2: add a description of your boundary condition
 *                              to the BC_Desc structure array in mm_names.h.
 *      mm_names.h              This structure includes information about the 
 *                              type of boundary condition, which equation it
 *                              applies to, what variables it is sensitive to,
 *                              whether to rotate mesh or momentum equations, 
 *                              etc.  It is very important that you fill out
 *                              this structure carefully, otherwise the code
 *                              won't know what to do.
 *  Boundary Condition  Step 3: add your BC case to the correct format listing 
 *                              for reading the necessary arguments from the
 *      mm_input_bc.c           input file in mm_input_bc.c.
 *
 *  Boundary Condition  Step 4: Add a function call (and a function) in the 
 *                              correct routine for evaluating your boundary
 *      bc_colloc.c             condition.  This will probably in bc_colloc.c
 *      bc_integ.c              for collocated conditions or bc_integ.c for 
 *                              strong or weak integrated conditions.
 *  Boundary Conditions Step 5: Add type into BC conflict lists.
 *  Boundary Condition  Step 6: use and enjoy your new boundary condition
 *
 * SPECIAL NOTES FOR SHELL ELEMENTS:
 *  When creating a boundary condition to be applied on a bulk equation
 *  which will use variables defined on an adjacent shell element block,
 *  place the new assembly function in this file below. Include the
 *  following in the function argument list:
 *	const int id_side,	(Get from elem_side_bc array)
 *      double xi[DIM],         (stu coordinates for local [bulk] element)
 *	const Exo_DB *exo       (Exodus database ptr)
 *
 *  The function should be similar to those for other integrated BC's,
 *  but must also include the following:
 *      Create a local array xi2[DIM].
 *      Create and allocate local array int *n_dof, which will need
 *	  MAX_VARIABLE_TYPES entries.
 *      Get remote (shell) element number from global elem_friends array.
 *      Call function load_neighbor_var_data(), which will require the
 *        above arguments to be passed in.
 *      Assemble the Residual and Jacobian terms AFTER the above steps.
 *      Free the n_dof array before returning.
 *
 *  For the assembly, the field variables on both element blocks are
 *  available in the usual fv-type structures. The shell elements will
 *  have their own basis functions, which can be accessed via bf[var]
 *  as usual. The n_dof array will contain the dof counts for the
 *  remote (shell) element variables, so when looping over these
 *  variable dofs, use n_dof[var] rather than ei->dof[var].
 */

/*****************************************************************************
 *									     *
 *        PLACE BULK BC ASSEMBLY FUNCTIONS USING SHELL VARIABLES HERE.       *
 *									     *
 *****************************************************************************/

/*
 * shell_surface_charge_bc(): Integrated BC on bulk potential equation.
 */
void
shell_surface_charge_bc(double func[DIM],
                        double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                        const double x_dot[MAX_PDIM],
                        const double tt,
                        const double dt,
                        const int id_side,
                        const double wt,
                        double xi[DIM],
                        const Exo_DB *exo,
 			const int sic_flag)
{
  //int eqn = R_POTENTIAL;
  int j, p, var, pvar;
  int el1 = ei->ielem;
  int el2, nf;
  int *n_dof = NULL;
  int n_dofptr[MAX_VARIABLE_TYPES][MDE], dof_map[MDE];
  double *n_esp;

  int a, b;
  dbl perm;                     /* electrical permittivity */
  dbl d_p_dT[MDE];              /* Temperature derivative of permittivity. */
  dbl d_p_dV[MDE];              /* Potential derivative of permittivity. */
  dbl d_p_dC[MAX_CONC][MDE];    /* Concentration derivative of permittivity. */
  dbl d_p_dX[DIM][MDE];         /* Spatial derivatives of permittivity. */
  dbl efield[MAX_PDIM];      /* electric field */
  double grad_phi_j, shell_qs;
 
  /* See if there is a friend for this element */
  nf = num_elem_friends[el1];
  if (nf == 0) return;

  /*
   * Get neighbor element number.
   * NOTE: For now, only one friend can be handled. If there is
   *       more than one, the first one will be used. This 
   *       will be fixed at a later date.
   */
  el2 = elem_friends[el1][0];
  if (nf != 1) WH(-1, "WARNING: Not set up for more than one element friend!");

  /*
   * DOF counts will be needed to construct sensitivities of
   * local shell equations to neighbor bulk element variables.
   * Allocate this array to save them in.
   */
  n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));

  /*
   * This call will collect all field variable and basis function data
   * for the neighbor element "el2" and populate the fv structure,
   * then repopulate the local element structure. Thus, as long as
   * no variable is defined on both blocks, all data from both
   * elements will be available in the fv structure.
   * For example, pressure in the neighbor element is accessed with fv->P.
   * This is done to simplify the job of writing shell equations.
   */
  load_neighbor_var_data(el1, el2, n_dof, dof_map, n_dofptr, id_side, xi, exo);

  /* Assemble the residual equation */
  /* Here, variables from the remote elements are used. */ 
 
  /* Get electrical conductivity */
  perm = mp->permittivity;
  if (mp->VoltageFormulation != V_PERMITTIVITY  && !sic_flag)
    {
      EH(-1,"trouble - SHELL_SURFACE_CHARGE with conductivity formulation");
    }
 
  memset( d_p_dV,0, sizeof(double)*MDE);
  memset( d_p_dT,0, sizeof(double)*MDE);
  memset( d_p_dX,0, sizeof(double)*MDE*DIM);
  memset( d_p_dC,0, sizeof(double)*MAX_CONC*MDE);
 
  for ( a=0; a<pd->Num_Dim; a++)
    { efield[a] = - fv->grad_V[a]; }

  shell_qs = 0;
  for (j=0 ; j<n_dof[SURF_CHARGE] ; j++)
    {
      n_esp = x_static + n_dofptr[SURF_CHARGE][j];
      shell_qs += *n_esp * bf[SURF_CHARGE]->phi[j];  
    }
 
  if ( af->Assemble_Residual )
    {
      func[0] = 0.5*shell_qs;
      if( sic_flag)
	{
	  for ( a=0 ; a<pd->Num_Dim ; a++)    {
	    func[0] += perm * fv->snormal[a] * efield[a] ;
	  }
	}
    }

  /* Assemble Jacobian sensitivities */
  if ( af->Assemble_Jacobian )
    {

      /* Surface charge: REMOTE sensitivity */
      var = SURF_CHARGE;
      pvar = upd->vp[var];
      if (pvar != -1)
	{
	  for (j=0; j<n_dof[var]; j++)
	    {
	      /*  Find the right variable in the bulk context */
	      for(el1=0 ; el1<n_dof[var] ; el1++)   {
		if(dof_map[j] == ei->dof_list[var][el1])el2=el1;
	      }
	      d_func[0][var][el2] += 0.5*bf[var]->phi[j];
	    }
	}
      if( sic_flag)
	{
	  var = VOLTAGE;
	  if (pd->v[var])
	    {
	      for(j=0; j<ei->dof[var]; j++)
		{
		  for ( p=0; p<pd->Num_Dim; p++)
		    { grad_phi_j = bf[var]->grad_phi[j][p]; 
		    d_func[0][var][j] += fv->snormal[p]*
		      ( -perm * grad_phi_j + d_p_dV[j] * efield[p]);
		    }
		}
	    }
     
	  var = TEMPERATURE;
	  if (pd->v[var])
	    {
	      for (j=0; j<ei->dof[var]; j++)
		{
		  for ( p=0; p<pd->Num_Dim; p++)
		    { grad_phi_j = bf[var]->grad_phi[j][p]; 
		    d_func[0][var][j] += fv->snormal[p]*
		      ( d_p_dT[j] * efield[p]);
		    }
		}
	    }
 
	  if (pd->v[MESH_DISPLACEMENT1])
	    {
	      for ( b=0; b<pd->Num_Dim; b++)
		{
		  var = MESH_DISPLACEMENT1+b;
		  for (j=0; j<ei->dof[var]; j++)
		    {
		      for (p=0; p<pd->Num_Dim; p++)
			{
			  d_func[0][var][j] += fv->snormal[p]*
			    ( -perm*fv->d_grad_V_dmesh[p][b][j]
			      +d_p_dX[b][j]*efield[p])
			    + fv->dsnormal_dx[p][b][j]*
			    (perm*efield[p]);
			}
		    }
		}
	    }
 
	}
    } /* End of if(af->Assemble_Jacobian) */

  /* Done */
  safe_free ((void *) n_dof);
  return;
}  /* shell_surface_charge_bc() */

void
surface_electric_field_bc(double R[MAX_PROB_VAR+MAX_CONC][MAX_NODES_PER_SIDE],
                          double J[MAX_PROB_VAR+MAX_CONC][MAX_PROB_VAR+MAX_CONC][MAX_NODES_PER_SIDE][MDE],
                          const int bm,
                          const int *n_dof,
                          const double wt)
{
  int i, j, p, q, dofs;
  int eqn, peqn, var, pvar;
  double boundary, em, det_J, phi_i, k;

  /* Check for active SURFACE_CHARGE equation */
  eqn = R_SURF_CHARGE;

  /* Only BOUNDARY terms are assembled here */
  if ( !(pd->e[eqn] & T_BOUNDARY) ) return;

  /* Unpack variables from structures for local convenience. */
  det_J = fv->sdet;
  peqn  = upd->ep[eqn];
  em    = pd->etm[eqn][(LOG2_BOUNDARY)];

  /* Get electrical conductivity */
  k = mp_glob[bm]->electrical_conductivity;

  /* Residuals */
  if (af->Assemble_Residual)
    {
      for (i = 0; i < ei->dof[eqn]; i++)
        {
          phi_i = bf[eqn]->phi[i];

          /* Boundary (jump) term */
          boundary = 0.0;
          for (p = 0; p < pd->Num_Dim; p++)
            {
              boundary += k * fv->snormal[p] * fv->grad_V[p];
            }
          boundary *= (phi_i * wt * det_J * em);
          R[peqn][i] += boundary;
        }
    }

  /* Jacobian */
  if (af->Assemble_Jacobian)
    {
      for (i = 0; i < ei->dof[eqn]; i++)
        {
          phi_i = bf[eqn]->phi[i];

          /* J_qs_V:  Remote voltage sensitivity */
          var = VOLTAGE;
          if (pd_glob[bm]->v[var])
            {
              pvar = upd->vp[var];
              dofs = n_dof[var];

              for (j = 0; j < dofs; j++)
                {
                  boundary = 0.0;
                  for (p = 0; p < pd->Num_Dim; p++)
                    {
                      boundary += k * fv->snormal[p] * bf[var]->grad_phi[j][p];
                    }
                  boundary *= (phi_i * wt * det_J * em);
                  J[peqn][pvar][i][j] += boundary;
                }
            }

          /* J_qs_x:  Remote mesh sensitivity */
          /* NOTE: It may be necessary to restrict this to one
           *       side of a shell if there are two sides! */
          for (q = 0; q < pd->Num_Dim; q++)
            {
              var = MESH_DISPLACEMENT1 + q;
              if (pd_glob[bm]->v[var])
                {
                  pvar = upd->vp[var];
                  dofs = n_dof[var];

                  for (j = 0; j < dofs; j++)
                    {
                      boundary = 0.0;
                      for (p = 0; p < pd->Num_Dim; p++)
                        {
                          boundary += det_J * ( k *
						(fv->dsnormal_dx[p][q][j] * fv->grad_V[p]
						 + fv->snormal[p] * fv->d_grad_V_dmesh[p][q][j]))
			    + k*fv->snormal[p]*fv->grad_V[p]*fv->dsurfdet_dx[q][j]; 
                        }
                      boundary *= (phi_i * wt * em);
                      J[peqn][pvar][i][j] += boundary;
                    }
                }
            }
        }
    }

  return;
}  /* End of surface_electric_field_bc() */

/********************************************************************************/
void
surface_acoustic_velocity_bc(double R[MAX_PROB_VAR+MAX_CONC][MAX_NODES_PER_SIDE],
			     double J[MAX_PROB_VAR+MAX_CONC][MAX_PROB_VAR+MAX_CONC][MAX_NODES_PER_SIDE][MDE],
			     const int bm,
			     const int *n_dof,
			     const double wt,
			     const double time_value)
{
  int i, j, p, q, dofs, w;
  int eqn, peqn, var, pvar;
  double boundary, em, det_J, phi_i;
  dbl R_imped;                               /* Acoustic impedance */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_R_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_R = &d_R_struct;
                                                                                  
  dbl wnum, kR_inv;                               /* Acoustic wavenumber */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_wnum_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_wnum = &d_wnum_struct;


  /* Check for active shell equation */
  eqn = R_SHELL_BDYVELO;

  /* Only BOUNDARY terms are assembled here */
  if ( !(pd->e[eqn] & T_BOUNDARY) ) return;

  /* Unpack variables from structures for local convenience. */
  det_J = fv->sdet;
  peqn  = upd->ep[eqn];
  em    = pd->etm[eqn][(LOG2_BOUNDARY)];

  /* Get acoustic properties */
  R_imped = acoustic_impedance( d_R, time_value );
  wnum = wave_number( d_wnum, time_value );
  kR_inv = 1./(wnum*R_imped);

  /* Residuals */
  if (af->Assemble_Residual)
    {
      for (i = 0; i < ei->dof[eqn]; i++)
        {
          phi_i = bf[eqn]->phi[i];

          /* Boundary (jump) term */
          boundary = 0.0;
          for (p = 0; p < pd->Num_Dim; p++)
            {
	      boundary -= (1.-SQUARE(fv->snormal[p]))*
		(SQUARE(fv->grad_apr[p]) + SQUARE(fv->grad_api[p]));
            }
	  boundary *= SQUARE(kR_inv);
          boundary *= (phi_i * wt * det_J * em);
          R[peqn][i] += boundary;
        }
    }

  /* Jacobian */
  if (af->Assemble_Jacobian)
    {
      for (i = 0; i < ei->dof[eqn]; i++)
        {
          phi_i = bf[eqn]->phi[i];

          /* real acoustic pressure sensitivity */
          var = ACOUS_PREAL;
          if (pd_glob[bm]->v[var])
            {
              pvar = upd->vp[var];
              dofs = n_dof[var];

              for (j = 0; j < dofs; j++)
                {
                  boundary = 0.0;
                  for (p = 0; p < pd->Num_Dim; p++)
                    {
                      boundary -= (1.-SQUARE(fv->snormal[p]))
			*2.*fv->grad_apr[p]*bf[var]->grad_phi[j][p];
                    }
	  	  boundary *= SQUARE(kR_inv);
                  boundary *= (phi_i * wt * det_J * em);
                  J[peqn][pvar][i][j] += boundary;
                }
            }

          /* imaginary acoustic pressure sensitivity */
          var = ACOUS_PIMAG;
          if (pd_glob[bm]->v[var])
            {
              pvar = upd->vp[var];
              dofs = n_dof[var];

              for (j = 0; j < dofs; j++)
                {
                  boundary = 0.0;
                  for (p = 0; p < pd->Num_Dim; p++)
                    {
                      boundary -= (1.-SQUARE(fv->snormal[p]))
			*2.*fv->grad_api[p]*bf[var]->grad_phi[j][p];
                    }
	  	  boundary *= SQUARE(kR_inv);
                  boundary *= (phi_i * wt * det_J * em);
                  J[peqn][pvar][i][j] += boundary;
                }
            }

          /* J_qs_x:  Remote mesh sensitivity */
          /* NOTE: It may be necessary to restrict this to one
           *       side of a shell if there are two sides! */
          for (q = 0; q < pd->Num_Dim; q++)
            {
              var = MESH_DISPLACEMENT1 + q;
              if (pd_glob[bm]->v[var])
                {
                  pvar = upd->vp[var];
                  dofs = n_dof[var];

                  for (j = 0; j < dofs; j++)
                    {
                      boundary = 0.0;
                      for (p = 0; p < pd->Num_Dim; p++)
                        {
			  boundary -= det_J * ((1.-SQUARE(fv->snormal[p]))*
					       2.*(fv->grad_apr[p]*fv->d_grad_apr_dmesh[p][q][j]
						   + fv->grad_api[p]*fv->d_grad_api_dmesh[p][q][j])
					       + (SQUARE(fv->grad_api[p])+SQUARE(fv->grad_apr[p]))*
					       (-2.*fv->snormal[p]*fv->dsnormal_dx[p][q][j]))
			    + fv->dsurfdet_dx[q][j]*(1.-SQUARE(fv->snormal[p]))
			    *(SQUARE(fv->grad_apr[p])+SQUARE(fv->grad_api[p]));
                        }
	  	      boundary *= SQUARE(kR_inv);
          	      boundary *= (phi_i * wt * em);
                      J[peqn][pvar][i][j] += boundary;
                    }
                }
            }
          /* temperature sensitivity */
          var = TEMPERATURE;
          if (pd_glob[bm]->v[var])
            {
              pvar = upd->vp[var];
              dofs = n_dof[var];

              for (j = 0; j < dofs; j++)
                {
                  boundary = 0.0;
                  for (p = 0; p < pd->Num_Dim; p++)
                    {
                      boundary -= (1.-SQUARE(fv->snormal[p]))
			*(SQUARE(fv->grad_apr[p])+SQUARE(fv->grad_api[p]));
                    }
	  	  boundary *= -kR_inv*SQUARE(kR_inv)*
		    (wnum*d_R->T[j]+R_imped*d_wnum->T[j]);
                  boundary *= (phi_i * wt * det_J * em);
                  J[peqn][pvar][i][j] += boundary;
                }
            }
          /* species sensitivity */
          var = MASS_FRACTION;
          if (pd_glob[bm]->v[var])
            {
              dofs = n_dof[var];
	      for ( w=0; w<pd->Num_Species_Eqn; w++)
		{
		  for (j = 0; j < dofs; j++)
		    {
		      boundary = 0.0;
		      for (p = 0; p < pd->Num_Dim; p++)
			{
			  boundary -= (1.-SQUARE(fv->snormal[p]))
			    *(SQUARE(fv->grad_apr[p])+SQUARE(fv->grad_api[p]));
			}
		      boundary *= -kR_inv*SQUARE(kR_inv)*
			(wnum*d_R->C[w][j]+R_imped*d_wnum->C[w][j]);
		      boundary *= (phi_i * wt * det_J * em);
		      J[peqn][MAX_PROB_VAR + w][i][j] += boundary;
		    }
		}
            }
        }
    }

  return;
}  /* End of surface_acoustic_velocity_bc() */
/********************************************************************************/
/*****************************************************************************/

void 
apply_surface_viscosity(double cfunc[MDE][DIM],
			double d_cfunc[MDE][DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
			double surf_shear_visc, /* Surface shear viscosity */
			double surf_dilat_visc, /* Surface dilitational visc */
                        const double time_start,      /* start to turn on surface viscosity */
                        const double time_full,       /* time that effect is fully on */
                        const double time_value,      /* current time */
			struct elem_side_bc_struct *elem_side_bc,
			const double wt, 
			double xi[DIM],
			const Exo_DB *exo,
			const int iconnect_ptr)
  
/******************************************************************************
 *
 *  Function which calculates the contributions of surface viscosity to a capillary
 *  free surface. This function is called by the CAPILLARY_SHEAR_VISC card
 *
 *  It is called from within a bulk element, when that element is evaluating surface
 *  boundary conditions.
 *
 *   cfunc[ldof][a]
 *
 ******************************************************************************/
{
  int i, j = 0, jj, k, id, var, a, b, eqn, I, ldof, q;
  int jvar;			/* Degree of freedom counter                 */
  int el1 = ei->ielem;
  int el2, nf;
  int *n_dof = NULL;
  int dof_map[MDE];
  int n_dofptr[MAX_VARIABLE_TYPES][MDE];
  double fac; 
  int dim = pd->Num_Dim;
  /*
   *  Debugging the signs on each term  -> note theory indicates that all of these
   *  signs should be positive.
   */
  // 
  //  Believe the sign is positive based on an analytical model involving an expanding bubble.
  //  The idea is that the viscosity term should impede the expansion of the bubble.
  double sgn1 = 1;
  //  Heuristically, I found that the sign on #2 should be negative. I found this out by ramping
  //  up the dilational viscosity coefficient for a quiestent bubble which was in equilibrium.
  //  For one value of bulk viscosity (with a positive sgn2) the flow started in one direction
  //  along the interface. For another value of the bulk viscosity the flow started going in the
  //  other direction. The flow only increased in strength. Therefore, we are creating energy with
  //  a positive sgn. Again, I don't see any theoretical basis for this switch. There may be an
  //  internal error within Goma having to do with the rh rule, since this is a tangential term.
  double sgn2 = -1;
  // Sign not explored
  double sgn3 = 1;
  // Setting this to one for the time being
  // Heuristically, I tried multiple things here and none worked
  double sgn4 = 1;

  // Setting this to neg 1 for the time being. Both positive and negative worked here.
  double sgn5 = -1;

  double d_grad_n_dn[DIM][DIM][DIM][MDE];
  double d_surfCurvatureDyadic_dn[DIM][DIM][DIM][MDE];
  double d_surfCurvatureDyadic_dmesh[DIM][DIM][DIM][MDE];
  double d_grad_s_v_dv[DIM][DIM][DIM][MDE];

  double d_grad_s_v_dmesh[DIM][DIM][DIM][MDE];
  int dofs;
  BASIS_FUNCTIONS_STRUCT *bfn;
  int r, p;

  double phi_j, dotdotTmp;
#ifdef DEBUG_HKM
  double normFixed[3];
  normFixed[0] = 1.0;
  normFixed[1] = 0.0;
  normFixed[2] = 0.0;
#endif

  /*
   * Adjust values for a time ramp
   */
  if (time_full > 0.0) {
    if (time_value < time_full) {
      if (time_value < time_start) {
         surf_shear_visc = 0.0;
         surf_dilat_visc = 0.0;
      } else {
        fac = (time_value - time_start) / (time_full - time_start);
        surf_shear_visc *= fac;
        surf_dilat_visc *= fac;
      }
    }
  }
  // Calculate the number of velocity dimensions
  int wim   = dim;
  if (pd->CoordinateSystem == SWIRLING || pd->CoordinateSystem == PROJECTED_CARTESIAN) 
    {
      wim = 3;
    }
  /***************************** EXECUTION BEGINS ******************************/

  /* See if there is a friend for this element */
  nf = num_elem_friends[el1];
  if (nf == 0) return;

  /*
   * Get neighbor element number.
   * NOTE: For now, only one friend can be handled. If there is
   *       more than one, the first one will be used. This 
   *       will be fixed at a later date.
   */
  el2 = elem_friends[el1][0];
  if (nf != 1) WH(-1, "WARNING: Not set up for more than one element friend!");

  /*
   * DOF counts will be needed to construct sensitivities of
   * local shell equations to neighbor bulk element variables.
   * Allocate this array to save them in.
   */
  n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));

  /*
   * This call will collect all field variable and basis function data
   * for the neighbor element "el2" and populate the fv structure,
   * then repopulate the local element structure. Thus, as long as
   * no variable is defined on both blocks, all data from both
   * elements will be available in the fv structure.
   * For example, pressure in the neighbor element is accessed with fv->P.
   * This is done to simplify the job of writing shell equations.
   */
  load_neighbor_var_data(el1, el2, n_dof, dof_map, n_dofptr, 
			 elem_side_bc->id_side, xi, exo);

  /* Assemble the residual equation */
  /* Here, variables from the remote elements are used. */ 



  /* Form grad_s (grad_s dot v) */
  double grad_s_div_s_v[DIM];
  for (a = 0; a < VIM; a++)
    {
      grad_s_div_s_v[a] = fv->grad_div_s_v[a];
      for (b = 0; b < VIM; b++)
	{
	  grad_s_div_s_v[a] += - fv->snormal[a] * fv->snormal[b] * fv->grad_div_s_v[b];
	}
    }
  

  // Formulate grad_s_v[a][b]         -- Gradient of velocity.  d (v_i) / d (x_j)
  //     grad_s_v[][]   = (I - n n ) grad_v[][] - latter operation is a matrix matrix multiply
  double grad_s_v[DIM][DIM];
  for (a = 0; a < VIM; a++)
    {
      for (b = 0; b < VIM; b++)
	{
	  grad_s_v[a][b] = fv->grad_v[a][b];
	  for (q = 0; q < VIM; q++)
	    {
	      grad_s_v[a][b] += - fv->snormal[a] * fv->snormal[q] * fv->grad_v[q][b];
	    }	    
	}
    }

  // Formulate Is = I - n n
  double Is[DIM][DIM];
  for (a = 0; a < VIM; a++)
    {
      for (b = 0; b < VIM; b++)
	{
	  Is[a][b] = delta(a,b) -  fv->n[a] * fv->n[b];
	}
    }

  // Formulate d_grad_n_dn

  memset(d_grad_n_dn, 0, sizeof(double)*DIM*DIM*DIM*MDE);
  memset(d_surfCurvatureDyadic_dn, 0, sizeof(double)*DIM*DIM*DIM*MDE);
  memset(d_surfCurvatureDyadic_dmesh, 0, sizeof(double)*DIM*DIM*DIM*MDE);
  memset(d_grad_s_v_dv, 0, sizeof(double)*DIM*DIM*DIM*MDE);
  memset(d_grad_s_v_dmesh, 0, sizeof(double)*DIM*DIM*DIM*MDE);


  eqn = SHELL_NORMAL1;
  // Check to see if dofs has the right num
  dofs = ei->dof[eqn];
  bfn =  bf[eqn];
  for (i = 0; i < VIM; i++)
    {
      for (j = 0; j < VIM; j++)
	{
	  for (r = 0; r < dim; r++)
	    {
	      for (k = 0; k < dofs; k++)
		{
		  d_grad_n_dn[i][j] [r][k] = bfn->grad_phi_e[k][r][i][j];
		}
	    }
	}
    }

  for (p = 0; p < VIM; p++)
    {
      for (q = 0; q < VIM; q++)
	{
	  for (r = 0; r < dim; r++)
	    {
	      for (j = 0; j < dofs; j++)
		{
		  d_surfCurvatureDyadic_dn[p][q] [r][j] = - d_grad_n_dn[p][q] [r][j];
		}
	    }
	  for ( a = 0; a < VIM; a++) 
	    {
	      for (j = 0; j < dofs; j++)
		{
		  phi_j = bfn->phi[j];
		  for (r = 0; r < dim; r++)
		    {
		      d_surfCurvatureDyadic_dn[p][q] [r][j] += fv->n[p] * fv->n[a] * d_grad_n_dn[a][q] [r][j];
		    }
		  d_surfCurvatureDyadic_dn[p][q] [p][j] += phi_j * fv->n[a] * fv->grad_n[a][q];
		  d_surfCurvatureDyadic_dn[p][q] [a][j] += fv->n[p] * phi_j * fv->grad_n[a][q];
		}

	    }
	}
    }

  var = MESH_DISPLACEMENT1;
  for (p = 0; p < VIM; p++)
    {
      for (q = 0; q < VIM; q++)
	{
	  for (r = 0; r < dim; r++)
	    {
	      for (j = 0; j < ei->dof[var]; j++)
		{
		  int jShell = -1;
		  int jjj;
		  int ndofMeshShell = n_dof[MESH_DISPLACEMENT1];
		  for (jjj = 0; jjj < ndofMeshShell; jjj++) {
		    if (dof_map[jjj] == j) {
		      jShell = jjj;
		    }
		  }
		  if (jShell >= 0) {
		    d_surfCurvatureDyadic_dmesh[p][q] [r][j] = - fv->d_grad_n_dmesh[p][q] [r][jShell];
		  } else {
		    d_surfCurvatureDyadic_dmesh[p][q] [r][j] = 0.0;
		  }
		}
	    }
	  for (a = 0; a < VIM; a++) 
	    {
	      for (j = 0; j < ei->dof[var]; j++)
		{
		  int jShell = -1;
		  int jjj;
		  int ndofMeshShell = n_dof[MESH_DISPLACEMENT1];
		  for (jjj = 0; jjj < ndofMeshShell; jjj++) {
		    if (dof_map[jjj] == j) {
		      jShell = jjj;
		    }
		  }
		  if (jShell >= 0) {
		    for (r = 0; r < dim; r++)
		      {
			d_surfCurvatureDyadic_dmesh[p][q] [r][j] += fv->n[p] * fv->n[a] * fv->d_grad_n_dmesh[a][q] [r][jShell];
		      }
		  }
		}
	    }
	}
    }


 
  bfn = bf[VELOCITY1];
  for (a = 0; a < VIM; a++)
    {
      for (b = 0; b < VIM; b++)
	{ 
	  for (r = 0; r < wim; r++)
	    {
	      var = VELOCITY1 + r;
	      for (j = 0; j < ei->dof[var]; j++)
		{
		  //  grad_v += (*esp->v[r][j]) * bfn->grad_phi_e[j][r] [p][q];
		  d_grad_s_v_dv[a][b] [r][j] =  bfn->grad_phi_e[j][r] [a][b];
		  for (q = 0; q < VIM; q++)
		    {
		      // grad_s_v[a][b] += - fv->snormal[a] * fv->snormal[q] * fv->grad_v[q][b];
		      d_grad_s_v_dv[a][b] [r][j] +=  - fv->snormal[a] * fv->snormal[q] * bfn->grad_phi_e[j][r] [q][b];
		    }	    
		}
	    }
	}
    }

  bfn = bf[MESH_DISPLACEMENT1];
  for (a = 0; a < VIM; a++)
    {
      for (b = 0; b < VIM; b++)
	{ 
	  for (r = 0; r < dim; r++)
	    {
	      for (j = 0; j < ei->dof[var]; j++)
		{
		  d_grad_s_v_dmesh[a][b] [r][j] = fv->d_grad_v_dmesh[a][b][r][j];
		  for (q = 0; q < VIM; q++)
		    {
		      d_grad_s_v_dmesh[a][b] [r][j] += - fv->snormal[a] * fv->snormal[q] * fv->d_grad_v_dmesh[q][b][r][j];
		      d_grad_s_v_dmesh[a][b] [r][j] += - fv->dsnormal_dx[a][r][j] * fv->snormal[q] * fv->grad_v[q][b];
		      d_grad_s_v_dmesh[a][b] [r][j] += - fv->snormal[a] * fv->dsnormal_dx[q][r][j] * fv->grad_v[q][b];
 		    }	    
		}
	    }
	}
    }


  /*  compute n x (I-nn).grad(curl_s_v_dot_n)  */

  eqn = VELOCITY1;
  for (i = 0; i < (int) elem_side_bc->num_nodes_on_side; i++)  
    {
      // Find the local node number of the node on the correct side
      id = (int) elem_side_bc->local_elem_node_id[i];
      // Find the global node number of the node
      I = Proc_Elem_Connect[iconnect_ptr + id];
      // Find the local degree of freedom of the first velocity unknown on that node
      ldof  = ei->ln_to_dof[eqn][id];
      if (Dolphin[I][VELOCITY1] > 0) 
	{
	  for (a = 0; a < VIM; a++)
	    {  
	      for (j = 0; j < ei->ielem_dim; j++)
		{
		  for (k = 0; k < VIM; k++)
		    {
		      // This piece has been fully tested
		      /* = mu_s * [n x I . grad(fv->n_dot_curl_v] */
		      cfunc[ldof][a] += sgn3 * surf_shear_visc * fv->snormal[j] * permute(j,k,a) * fv->grad_n_dot_curl_s_v[k];
		    }
		}
	    }
	}
    }
   
  /*  compute (kappa_s + mu_s) grad_s(grad_s dot v) = (kappa_s + mu_s) grad_s(div_s v)  */
  eqn = VELOCITY1;
  for (i = 0; i < (int) elem_side_bc->num_nodes_on_side; i++)  
    {
      // Find the local node number of the node on the correct side
      id = (int) elem_side_bc->local_elem_node_id[i];
      // Find the global node number of the node
      I = Proc_Elem_Connect[iconnect_ptr + id];
      // Find the local degree of freedom of the first velocity unknown on that node
      ldof  = ei->ln_to_dof[eqn][id];
      if (Dolphin[I][VELOCITY1] > 0) 
	{
	  for (a = 0; a < VIM; a++)
	    {
	      cfunc[ldof][a] += sgn2 * (surf_dilat_visc + surf_shear_visc) * grad_s_div_s_v[a];
	    }
	}
    }

  /*  compute (kappa_s + mu_s) 2 H n (div_s v)   */

  eqn = VELOCITY1;
  for (i = 0; i < (int) elem_side_bc->num_nodes_on_side; i++)  
    {
      // Find the local node number of the node on the correct side
      id = (int) elem_side_bc->local_elem_node_id[i];
      // Find the global node number of the node
      I = Proc_Elem_Connect[iconnect_ptr + id];
      // Find the local degree of freedom of the first velocity unknown on that node
      ldof  = ei->ln_to_dof[eqn][id];
      if (Dolphin[I][VELOCITY1] > 0) 
	{
	  for (a = 0; a < VIM; a++)
	    {
	      cfunc[ldof][a] += sgn1 * (surf_dilat_visc + surf_shear_visc) * 2 * fv->curv * fv->snormal[a] * fv->div_s_v;
	    }
	}
    }


  /*  compute  mu_s ( 2 n ((b[][] - 2 H Is[][]) DotDot (grad_s_v[][]) ))  */

  eqn = VELOCITY1;
  for (i = 0; i < (int) elem_side_bc->num_nodes_on_side; i++)  
    {
      // Find the local node number of the node on the correct side
      id = (int) elem_side_bc->local_elem_node_id[i];
      // Find the global node number of the node
      I = Proc_Elem_Connect[iconnect_ptr + id];
      // Find the local degree of freedom of the first velocity unknown on that node
      ldof  = ei->ln_to_dof[eqn][id];
      if (Dolphin[I][VELOCITY1] > 0) 
	{
	  for (a = 0; a < VIM; a++)
	    {
	      dotdotTmp = 0.0;
	      for (p = 0; p < VIM; p++) 
		{
		  for (q = 0; q < VIM; q++) 
		    {
		      dotdotTmp += (fv->surfCurvatureDyadic[p][q] - 2 * fv->curv * Is[p][q]) * grad_s_v[q][p];
		    }
		}
	      cfunc[ldof][a] += sgn4 * surf_shear_visc * 2.0 * dotdotTmp * fv->n[a];
	    }
	}
    }

  /*  compute  mu_s ( - 2 ((b[][] - 2 H Is[][]) Dot (grad_s_v[][]) Dot n))  */

  eqn = VELOCITY1;
  for (i = 0; i < (int) elem_side_bc->num_nodes_on_side; i++)  
    {
      // Find the local node number of the node on the correct side
      id = (int) elem_side_bc->local_elem_node_id[i];
      // Find the global node number of the node
      I = Proc_Elem_Connect[iconnect_ptr + id];
      // Find the local degree of freedom of the first velocity unknown on that node
      ldof  = ei->ln_to_dof[eqn][id];
      if (Dolphin[I][VELOCITY1] > 0) 
	{
	  for (a = 0; a < VIM; a++)
	    {
	      for (q = 0; q < VIM; q++) 
		{
		  cfunc[ldof][a] -= sgn5 * (surf_shear_visc) * 2.000 * (fv->surfCurvatureDyadic[a][q] - 2.00 * fv->curv * Is[a][q]) * fv->grad_v_dot_n[q];
		}
	    }
	}
    }




  if (af->Assemble_Jacobian)
    { 
      eqn = VELOCITY1;
      for (i = 0; i < (int) elem_side_bc->num_nodes_on_side; i++)  {
	id = (int) elem_side_bc->local_elem_node_id[i];
	I = Proc_Elem_Connect[iconnect_ptr + id];
	ldof  = ei->ln_to_dof[eqn][id];
	if (Dolphin[I][VELOCITY1] > 0 ) {
	  
	  /* 
	   *  Evaluate sensitivity to displacements d()/dx 
	   */
	  for (jvar = 0; jvar < ei->ielem_dim; jvar++)
	    {
	      var = MESH_DISPLACEMENT1 + jvar;
	      if (pd->v[var]) 
		{
		  /*
		   *  Here j refers to the bulk element dofs', which go from 0 to 8.
		   */
		  for (j = 0; j < ei->dof[var]; j++)
		    {
		      /*
		       *  Find the inverse of the dof_map that comes out of load_neighbor_var_data()
		       */
		      int jShell = -1;
		      int jjj;
		      int ndofMeshShell = n_dof[MESH_DISPLACEMENT1];
		      for (jjj = 0; jjj < ndofMeshShell; jjj++) {
			if (dof_map[jjj] == j) {
			  jShell = jjj;
			}
		      }
		      for (a = 0; a < VIM; a++)
			{  
			  for (jj = 0; jj < ei->ielem_dim; jj++)
			    {
			      for (k = 0; k < VIM; k++)
				{
				  //   HKM 3/8/10 -> This has been checked out to be correct
				  /*   residual_piece  = mu_s * [n x  grad(fv->n_dot_curl_v[p]] */
				  /*   cfunc[ldof][a] += surf_shear_visc * fv->snormal[j] * permute(j,k,a) * fv->grad_n_dot_curl_s_v[k]; */
				  if (jShell >= 0)
				    {
				      d_cfunc[ldof][a][var][j] += sgn3 * surf_shear_visc * 
					(fv->snormal[jj] * permute(jj,k,a) * fv->d_grad_n_dot_curl_s_v_dmesh[k][jvar][jShell]);
				    }
				  d_cfunc[ldof][a][var][j] += sgn3 * surf_shear_visc *
				    (fv->dsnormal_dx[jj][jvar][j] * permute(jj,k,a) * fv->grad_n_dot_curl_s_v[k]);				  
				}
			    }
			}
		    }
		}
	    }

	  /*
	   *  We don't check for pd->var[] here because the shell variable will not be active on the
	   *  bulk element. However, ei->dof[] will be nonzero. 
	   *  We must check that we are using the correct indexing into this variable, as it will only exist
	   *  on the shell
	   *   HKM 3/8/2010 -> This interaction has been shown to be correct via the numerical jacobian checker
	   *           There is a lot going on here. j is the local dof number in the shell. bf[var] is the basis function
	   *           for the shell. 
	   *            j turns out to be the ldof number in the bulk element for the shell unknowns. It turns out that it gets
	   *            put into the Jacobian in the right spot.
	   *            I wonder if this is fortuitous? I could probably check dof_map[] versus local element node numbers for ei->dof[var]
	   *            and check this out.
	   */
	  var = N_DOT_CURL_V;
	  for (j = 0; j < ei->dof[var]; j++)
	    {
	      for (a = 0; a < VIM; a++)
		{  
		  for (jj = 0; jj < VIM; jj++)
		    {
		      for (k = 0; k < VIM; k++)
			{
			  //    double tmpI =  surf_shear_visc * (fv->snormal[jj] * permute(jj,k,a) * bf[var]->grad_phi[j][k]); 
			  //    if (tmpI != 0.0) {
			  //      printf("[ldof=%d][a=%d][var=%d][j=%d] = %g\n", ldof, a, var, j, tmpI);
			  //    }
			  //   HKM 3/8/10 -> This has been checked out to be correct
			  /*   residual_piece  = mu_s * [n x  grad(fv->n_dot_curl_v[p]] */
			  /*   cfunc[ldof][a] += surf_shear_visc * fv->snormal[j] * permute(j,k,a) * fv->grad_n_dot_curl_s_v[k]; */
			  d_cfunc[ldof][a][var][j] += sgn3 * surf_shear_visc * (fv->snormal[jj] * permute(jj,k,a) * bf[var]->grad_phi[j][k]); 
			}
		    }
		}
	    }
	}
      }

    for (i = 0; i < (int) elem_side_bc->num_nodes_on_side; i++)  {
        id = (int) elem_side_bc->local_elem_node_id[i];
	I = Proc_Elem_Connect[iconnect_ptr + id];
	ldof  = ei->ln_to_dof[eqn][id];
	if (Dolphin[I][VELOCITY1] > 0 ) {
	  
	  /* 
	   *  Evaluate sensitivity to displacements d()/dx 
	   */
	  for (jvar = 0; jvar < ei->ielem_dim; jvar++)
	    {
	      var = MESH_DISPLACEMENT1 + jvar;
	      if (pd->v[var]) 
		{
		  /*
		   *  Here j refers to the bulk element dofs', which go from 0 to 8.
		   */
		  for (j = 0; j < ei->dof[var]; j++)
		    {
		      /*
		       *  Find the inverse of the dof_map that comes out of load_neighbor_var_data()
		       */
		      int jShell = -1;
		      int jjj;
		      int ndofMeshShell = n_dof[MESH_DISPLACEMENT1];
		      for (jjj = 0; jjj < ndofMeshShell; jjj++) {
			if (dof_map[jjj] == j) {
			  jShell = jjj;
			}
		      }
		      for (a = 0; a < VIM; a++)
			{			 
			  for (b = 0; b < VIM; b++)
			    {
			      //   HKM 3/8/10 -> This has been checked out to be correct
			      /*   residual_piece  = (surf_dilat_visc + surf_shear_visc) * grad_s_div_s_v[a]; */
			     
			      if (jShell >= 0)
				{
				  d_cfunc[ldof][a][var][j] += sgn2 * (surf_dilat_visc + surf_shear_visc)  * 
				    (delta(a,b) * fv->d_grad_div_s_v_dmesh[b][jvar][jShell] - fv->snormal[a] * fv->snormal[b] * fv->d_grad_div_s_v_dmesh[b][jvar][jShell]);
				}
			      /* - fv->snormal[a] * fv->snormal[b] * fv->grad_div_s_v[b] */;
			      d_cfunc[ldof][a][var][j] += - sgn2 * (surf_dilat_visc + surf_shear_visc) * 
				((fv->dsnormal_dx[a][jvar][j] * fv->snormal[b] * fv->grad_div_s_v[b]) +	(fv->snormal[a] * fv->dsnormal_dx[b][jvar][j] * fv->grad_div_s_v[b]));
			    }
			  
			}
		    }
		}
	    }

	  /*  (surf_dilat_visc + surf_shear_visc) * grad_s_div_s_v[a]; 
	   *        Gradient of this term wrt the SHELL_SURF_DIV_V variables which are defined on the shell element      
	   */
	  var = SHELL_SURF_DIV_V;
	  for (j = 0; j < ei->dof[var]; j++)
	    {
	      for (a = 0; a < VIM; a++)
		{  
		  d_cfunc[ldof][a][var][j] += sgn2 * (surf_dilat_visc + surf_shear_visc) * (bf[var]->grad_phi[j][a]); 

		  for (b = 0; b < VIM; b++)
		    {
		      d_cfunc[ldof][a][var][j] += - sgn2 * (surf_dilat_visc + surf_shear_visc) * fv->snormal[a] * fv->snormal[b] *( bf[var]->grad_phi[j][b]); 
		    }
		}

	    }  
	}
      }


    /* (surf_dilat_visc + surf_shear_visc) * 2 * fv->curv * fv->snormal[a] * fv->div_s_v; */
    for (i = 0; i < (int) elem_side_bc->num_nodes_on_side; i++)  {
	id = (int) elem_side_bc->local_elem_node_id[i];
	I = Proc_Elem_Connect[iconnect_ptr + id];
	ldof  = ei->ln_to_dof[eqn][id];
	if (Dolphin[I][VELOCITY1] > 0 ) {
	  
	  /* 
	   *  Evaluate sensitivity to displacements d()/dx 
	   *     ei->ielem_dim -> check this is equal to 2!
	   */
	  for (jvar = 0; jvar < ei->ielem_dim; jvar++)
	    {
	      var = MESH_DISPLACEMENT1 + jvar;
	      if (pd->v[var]) 
		{
		  /*
		   *  Here j refers to the bulk element dofs', which go from 0 to 8.
		   */
		  for (j = 0; j < ei->dof[var]; j++)
		    {
		      for (a = 0; a < VIM; a++)
			{			 
			  d_cfunc[ldof][a][var][j] += sgn1 * (surf_dilat_visc + surf_shear_visc) * 2 * fv->curv * fv->div_s_v
			    * ( fv->dsnormal_dx[a][jvar][j] );
			}
		    }
		}
	    }

	  /*  (surf_dilat_visc + surf_shear_visc) * 2 * fv->curv * fv->snormal[a] * fv->div_s_v; 
	   *        Gradient of this term wrt the SHELL_SURF_DIV_V variables which are defined on the shell element      
	   */
	  var = SHELL_SURF_DIV_V;
	  for (j = 0; j < ei->dof[var]; j++)
	    {
	      for (a = 0; a < VIM; a++)
		{  
		  d_cfunc[ldof][a][var][j] += sgn1 * (surf_dilat_visc + surf_shear_visc) * 2 * fv->curv * fv->snormal[a] * (bf[var]->phi[j]); 
		}

	    } 

	  /*  (surf_dilat_visc + surf_shear_visc) * grad_s_div_s_v[a]; 
	   *        Gradient of this term wrt the SHELL_SURF_DIV_V variables which are defined on the shell element      
	   */
	  var = SHELL_SURF_CURV;
	  for (j = 0; j < ei->dof[var]; j++)
	    {
	      for (a = 0; a < VIM; a++)
		{  
		  d_cfunc[ldof][a][var][j] += sgn1 * (surf_dilat_visc + surf_shear_visc) * 2 * fv->div_s_v * fv->snormal[a] * (bf[var]->phi[j]); 
		}

	    }
	}
      }


    // Jacobian terms for the last surface rheology expression
    //                 cfunc[ldof][a] += (surf_shear_visc) * 2.0 * dotdotTmp * fv->n[a];
    for (i = 0; i < (int) elem_side_bc->num_nodes_on_side; i++)
      {
	// Find the local node number of the node on the correct side
	id = (int) elem_side_bc->local_elem_node_id[i];
	// Find the global node number of the node
	I = Proc_Elem_Connect[iconnect_ptr + id];
	// Find the local degree of freedom of the first velocity unknown on that node
	ldof  = ei->ln_to_dof[eqn][id];
	for (b = 0; b < pd->Num_Dim; b++)
	  {
	    var = VELOCITY1 + b;
	    for (j = 0; j < ei->dof[var]; j++)
	      if (Dolphin[I][VELOCITY1] > 0) 
		{ 
		  phi_j = bf[var]->phi[j];
		  for (a = 0; a < VIM; a++)
		    {
		      for (p = 0; p < VIM; p++) 
			{
			  for (q = 0; q < VIM; q++) 
			    {
			      // HKM  - looks good
			      d_cfunc[ldof][a][var][j] += sgn4 * (surf_shear_visc * 2.0) * 
				(fv->surfCurvatureDyadic[p][q] - 2 * fv->curv * Is[p][q]) * d_grad_s_v_dv[q][p] [b][j] * fv->n[a];
			    }
			}    
		    }
		}
	  }

	for (b = 0; b < pd->Num_Dim; b++)
	  {
	    var = SHELL_NORMAL1 + b;
	    for (j = 0; j < ei->dof[var]; j++)
	      {
		phi_j = bf[var]->phi[j];
		for (a = 0; a < VIM; a++)
		  {
		    dotdotTmp = 0.0;
		    for (p = 0; p < VIM; p++) 
		      {
			for (q = 0; q < VIM; q++) 
			  {
			    dotdotTmp += (fv->surfCurvatureDyadic[p][q] - 2 * fv->curv * Is[p][q]) * grad_s_v[q][p];
			    d_cfunc[ldof][a][var][j] += sgn4 * (surf_shear_visc) * 2.0 * (d_surfCurvatureDyadic_dn[p][q] [b][j] * grad_s_v[q][p]) *  fv->n[a];
			    d_cfunc[ldof][a][var][j] += sgn4 * (surf_shear_visc) * 2.0 * (- 2 * fv->curv * (- fv->n[p] * phi_j * delta(q,b)) * grad_s_v[q][p]) * fv->n[a];
			    d_cfunc[ldof][a][var][j] += sgn4 * (surf_shear_visc) * 2.0 * (- 2 * fv->curv * (- fv->n[q] * phi_j * delta(p,b)) * grad_s_v[q][p]) * fv->n[a];
			  }
		      }    
		    
		    d_cfunc[ldof][a][var][j] += sgn4 * (surf_shear_visc) * 2.0 * dotdotTmp * phi_j * delta(a,b);
		  }
		
	      }
	  }

	var = SHELL_SURF_CURV;
	for (j = 0; j < ei->dof[var]; j++)
	  {
	    phi_j = bf[var]->phi[j];
	    for (a = 0; a < VIM; a++)
	      {  
		for (p = 0; p < VIM; p++) 
		  {
		    for (q = 0; q < VIM; q++) 
		      {
			// curvature dependence to first term
			d_cfunc[ldof][a][var][j] += sgn4 * (surf_shear_visc) * 2.0 * fv->n[a] * (-2.0 * phi_j * Is[p][q] * grad_s_v[q][p]);
		      }
		  }
	      }
	  }

	for (b = 0; b < dim; b++)
	  {
	    var = MESH_DISPLACEMENT1 + b;
	    for (j = 0; j < ei->dof[var]; j++)
	      if (Dolphin[I][MESH_DISPLACEMENT1] > 0) 
		{  
	
		  phi_j = bf[var]->phi[j];
		  for (a = 0; a < VIM; a++)
		    { 
		      for (p = 0; p < VIM; p++) 
			{
			  for (q = 0; q < VIM; q++) 
			    {   
			      // dotdotTmp += (fv->surfCurvatureDyadic[p][q] - 2 * fv->curv * Is[p][q]) * grad_s_v[q][p];
			      d_cfunc[ldof][a][var][j] += sgn4 * (surf_shear_visc * 2.0) * (d_surfCurvatureDyadic_dmesh[p][q][b][j] * grad_s_v[q][p]) * fv->n[a];
			      d_cfunc[ldof][a][var][j] += sgn4 * (surf_shear_visc * 2.0) * ((fv->surfCurvatureDyadic[p][q] - 2 * fv->curv * Is[p][q]) * d_grad_s_v_dmesh[q][p][b][j]) * fv->n[a];
			    }
			}		      
		    }
		}
	  }



	for (b = 0; b < pd->Num_Dim; b++)
	  {
	    var = GRAD_S_V_DOT_N1 + b;
	    if (ei->dof[var] == 0) {
	      printf("shouldn't behere \n");
	      exit(-1);
	    }
	    for (j = 0; j < ei->dof[var]; j++)
	      {
		phi_j = bf[var]->phi[j];
		for (a = 0; a < VIM; a++)
		  {  
		    for (q = 0; q < VIM; q++) 
		      {
			// cfunc[ldof][a] -= sgn5 * (surf_shear_visc) * 2.000 * (fv->surfCurvatureDyadic[a][q] - 2.00 * fv->curv * Is[a][q]) * fv->grad_v_dot_n[q];
			d_cfunc[ldof][a][var][j] -= sgn5 * (surf_shear_visc) * 2.000 * (fv->surfCurvatureDyadic[a][q] - 2.00 * fv->curv * Is[a][q]) * delta(q,b) * phi_j;
		      }
		  }
	      }
	  }

	var = SHELL_SURF_CURV;
	for (j = 0; j < ei->dof[var]; j++)
	  {
	    phi_j = bf[var]->phi[j];
	    for (a = 0; a < VIM; a++)
	      {  
		for (q = 0; q < VIM; q++) 
		  {
		    // cfunc[ldof][a] -= sgn5 * (surf_shear_visc) * 2.000 * (fv->surfCurvatureDyadic[a][q] - 2.00 * fv->curv * Is[a][q]) * fv->grad_v_dot_n[q];
		    d_cfunc[ldof][a][var][j] -= sgn5 * (surf_shear_visc) * 2.000 * (fv->surfCurvatureDyadic[a][q] - 2.00 * phi_j * Is[a][q]) * fv->grad_v_dot_n[q];
		  }
	      }
	  }
	

	for (b = 0; b < pd->Num_Dim; b++)
	  {
	    var = SHELL_NORMAL1 + b;
	    for (j = 0; j < ei->dof[var]; j++)
	      {
		phi_j = bf[var]->phi[j];
		for (a = 0; a < VIM; a++)
		  {
		    for (q = 0; q < VIM; q++) 
		      {			 
			d_cfunc[ldof][a][var][j] -= sgn5 * (surf_shear_visc) * 2.000 * (d_surfCurvatureDyadic_dn[a][q] [b][j])              * fv->grad_v_dot_n[q];
			d_cfunc[ldof][a][var][j] -= sgn5 * (surf_shear_visc) * 2.000 * (- 2 * fv->curv * (- fv->n[a] * phi_j * delta(q,b))) * fv->grad_v_dot_n[q];
			d_cfunc[ldof][a][var][j] -= sgn5 * (surf_shear_visc) * 2.000 * (- 2 * fv->curv * (- fv->n[q] * phi_j * delta(a,b))) * fv->grad_v_dot_n[q];
		      }
		  }
	      }
	  }

	for (b = 0; b < dim; b++)
	  {
	    var = MESH_DISPLACEMENT1 + b;
	    for (j = 0; j < ei->dof[var]; j++)
	      if (Dolphin[I][MESH_DISPLACEMENT1] > 0) 
		{ 
		  phi_j = bf[var]->phi[j];
		  for (a = 0; a < VIM; a++)
		    { 
		      for (q = 0; q < VIM; q++) 
			{
			  d_cfunc[ldof][a][var][j] -= sgn5 * (surf_shear_visc) * 2.000 * (d_surfCurvatureDyadic_dmesh[a][q][b][j]) * fv->grad_v_dot_n[q];
			}
		    }
		}
	  }

      }
    

    }
  
} /* END of routine apply_surface_viscosity */
/*****************************************************************************/

/*
 * This routine is used for the boundary condition SHELL_FLUID_STRESS that results in
 * a normal and tangential fluid stress loading on the 2D structural shell equations. 
 * It is written in the spirit of put_fluid_stress_in_solid invoked by the FLUID_SOLID
 * BC. The one exception is that this routine actually performs the rotation of the
 * residual rather than relying on the rotation mechanism in GOMA.  NOTE:  ONLY WRITTEN
 * FOR 2D SIMULATION!!!
 */

void 
put_fluid_stress_on_shell(int id, /* local element node number for the 
				   * current node whose residual contribution
				   * is being sought                        */
			  int id_shell, /* local shell element node number corresponding to id */
			  int I, /* Global node number                      */
			  int ielem_dim, /* physical dimension of the elem  */
			  double resid_vector[], /* Residual vector         */
			  int local_node_list_fs[], /* MDE list to keep track
						     * of nodes at which 
						     * solid contributions 
						     * have been transfered
						     * to liquid (fluid-solid
						     * boundaries)          */
			  double scale) /* Scale factor, nondimension       */
{
  int j_id, dim, var, pvar, q;
  int id_dofmom1, id_dofmom2,  id_dofshell, offset;
  int peqn_mom1, peqn_mom2, peqn_shell;
  int ieqn_mom1, ieqn_mom2, ieqn_shell;
  int curvature_fixed, tension_fixed;
  NODE_INFO_STRUCT *node = Nodes[I];
  NODAL_VARS_STRUCT *nv = node->Nodal_Vars_Info;

  dim = pd->Num_Dim;

  if (node->DBC) {
     
    offset = get_nodal_unknown_offset(nv, R_SHELL_CURVATURE, -2, 0, NULL);
    curvature_fixed = (offset >= 0) && (node->DBC[offset] != -1);
    offset = get_nodal_unknown_offset(nv, R_SHELL_TENSION, -2, 0, NULL);
    tension_fixed   = (offset >= 0) && (node->DBC[offset] != -1);	
  } else {
    curvature_fixed = tension_fixed = 0;
  }


  /*
   * if you are in the shell phase, return without doing anything
   * In the shell phase, there are no fluid momentum equations.
   * This probably is not needed since it is impossible with CUBIT to
   * put attach a sideset to a shell element, and hence all sidesets
   * are one-sided. 
   */
  if (!pd->e[R_MOMENTUM1]) return;

  /* Hmm, unlike in the FLUID_SOLID case, we've no shell equations in the 
   * bulk so we cannot use ei struct like this 
   * id_dofmom = ei->ln_to_dof[R_MOMENTUM1][id];
   * id_dofshell = ei->ln_to_dof[R_SHELL_CURVATURE][i];
   *
   * So try this because it is only 2D (assumption that there is a shell
   * variable at every surface node
   */
  id_dofmom1 = ei->ln_to_dof[R_MOMENTUM1][id];
  id_dofshell = id_shell;  /* This is the imperfect assumption */


  /*
   * if this nodal contribution has already been added to fluid momentum
   * equation (i.e. we are at a corner on the second side) return
   * without doing anything
   */
  if (local_node_list_fs[id] == -1) {
    local_node_list_fs[id] = 1;
  } else {
    return;
  }

  /*
   * check to make sure that velocity
   * dofs and shell dofs exist at this node
   */
  if (Dolphin[I][R_MOMENTUM1] <= 0)     return;
  if (Dolphin[I][R_SHELL_CURVATURE] <= 0) return;
  if (Dolphin[I][R_SHELL_TENSION] <= 0) return;

  /*
   * Add local contribution to fluid
   * momentum into local contribution for shell momentum.
   */
  if (af->Assemble_Residual) {
    /* don't do this if the curvature is fixed for this node */
    if (!curvature_fixed ) {
      ieqn_shell = R_SHELL_CURVATURE;
      ieqn_mom1   = R_MOMENTUM1;
      ieqn_mom2   = R_MOMENTUM2;
      id_dofmom1 = ei->ln_to_dof[ieqn_mom1][id];
      id_dofmom2 = ei->ln_to_dof[ieqn_mom2][id];
      id_dofshell = id_shell; 
      lec->R[upd->ep[ieqn_shell]][id_dofshell] +=
	scale * (fv->snormal[0] * lec->R[upd->ep[ieqn_mom1]][id_dofmom1] + 
		 fv->snormal[1] * lec->R[upd->ep[ieqn_mom2]][id_dofmom2]);
    }
    if (!tension_fixed ) {
      ieqn_shell = R_SHELL_TENSION;
      ieqn_mom1   = R_MOMENTUM1;
      ieqn_mom2   = R_MOMENTUM2;
      id_dofmom1 = ei->ln_to_dof[ieqn_mom1][id];
      id_dofmom2 = ei->ln_to_dof[ieqn_mom2][id];
      id_dofshell = id_shell; 
      lec->R[upd->ep[ieqn_shell]][id_dofshell] +=
	scale * (fv->stangent[0][0] * lec->R[upd->ep[ieqn_mom1]][id_dofmom1] + 
		 fv->stangent[0][1] * lec->R[upd->ep[ieqn_mom2]][id_dofmom2]);
    }
  }
    
  /*
   * loop over directions and add local contribution to solid
   * momentum into local contribution for liquid momentum
   */
  if (af->Assemble_Jacobian)
    {
      if (! curvature_fixed)
	{
	  ieqn_shell = R_SHELL_CURVATURE;
	  peqn_shell = upd->ep[ieqn_shell];
	  ieqn_mom1   = R_MOMENTUM1;
	  ieqn_mom2   = R_MOMENTUM2;
	  id_dofmom1 = ei->ln_to_dof[ieqn_mom1][id];
	  peqn_mom1 = upd->ep[ieqn_mom1];
	  id_dofmom2 = ei->ln_to_dof[ieqn_mom2][id];
	  peqn_mom2 = upd->ep[ieqn_mom2];
	  id_dofshell = id_shell; 


	  /* Add contributions due to all nodal sensitivities in shell element */
	    
	  /*
	   * local J_m_d -> J_d_d
	   */
	  for ( q=0; q<dim; q++)
	    {
	      var = MESH_DISPLACEMENT1+q;
	      if ( pd->v[var] )
		{
		  pvar = upd->vp[var];
		  for ( j_id=0; j_id<ei->dof[var]; j_id++)
		    {
		      lec->J[peqn_shell][pvar][id_dofshell][j_id] +=
			scale*(fv->snormal[0] * lec->J[peqn_mom1][pvar][id_dofmom1][j_id] +
			       fv->snormal[1] * lec->J[peqn_mom2][pvar][id_dofmom2][j_id] +
			       fv->dsnormal_dx[0][q][j_id] * lec->R[upd->ep[ieqn_mom1]][id_dofmom1] +
			       fv->dsnormal_dx[1][q][j_id] * lec->R[upd->ep[ieqn_mom2]][id_dofmom2]);
		    }
		}
	    }
	    
	  /*
	   * local J_m_P -> J_d_P
	   */
	  var = PRESSURE;
	  if ( pd->v[var] )
	    {
	      pvar = upd->vp[var];
	      for ( j_id=0; j_id<ei->dof[var]; j_id++)
		{			
		  lec->J[peqn_shell][pvar][id_dofshell][j_id] +=
		    scale*(fv->snormal[0] * lec->J[peqn_mom1][pvar][id_dofmom1][j_id] +
			   fv->snormal[1] * lec->J[peqn_mom2][pvar][id_dofmom2][j_id]);
		}
	    }
	    
	    
	  /*
	   * local J_m_v -> J_d_v  
	   */
	  for ( q=0; q<dim; q++)
	    {
	      var = VELOCITY1+q;
	      if ( pd->v[var] )
		{
		  pvar = upd->vp[var];
		  for ( j_id=0; j_id<ei->dof[var]; j_id++)
		    {
		      lec->J[peqn_shell][pvar][id_dofshell][j_id] +=
			scale*(fv->snormal[0] * lec->J[peqn_mom1][pvar][id_dofmom1][j_id] +
			       fv->snormal[1] * lec->J[peqn_mom2][pvar][id_dofmom2][j_id]);
		    }
		}
	    }
	    
	}

      if (! tension_fixed)
	{

	  ieqn_shell = R_SHELL_TENSION;
	  peqn_shell = upd->ep[ieqn_shell];
	  ieqn_mom1   = R_MOMENTUM1;
	  ieqn_mom2   = R_MOMENTUM2;
	  id_dofmom1 = ei->ln_to_dof[ieqn_mom1][id];
	  peqn_mom1 = upd->ep[ieqn_mom1];
	  id_dofmom2 = ei->ln_to_dof[ieqn_mom2][id];
	  peqn_mom2 = upd->ep[ieqn_mom2];
	  id_dofshell = id_shell; 


	  /* Add contributions due to all nodal sensitivities in shell element */
	    
	  /*
	   * local J_m_d -> J_d_d
	   */
	  for ( q=0; q<dim; q++)
	    {
	      var = MESH_DISPLACEMENT1+q;
	      if ( pd->v[var] )
		{
		  pvar = upd->vp[var];
		  for ( j_id=0; j_id<ei->dof[var]; j_id++)
		    {
		      lec->J[peqn_shell][pvar][id_dofshell][j_id] +=
			scale*(fv->stangent[0][0] * lec->J[peqn_mom1][pvar][id_dofmom1][j_id] +
			       fv->stangent[0][1] * lec->J[peqn_mom2][pvar][id_dofmom2][j_id] +
			       fv->dstangent_dx[0][0][q][j_id] * lec->R[upd->ep[ieqn_mom1]][id_dofmom1] +
			       fv->dstangent_dx[0][1][q][j_id] * lec->R[upd->ep[ieqn_mom2]][id_dofmom2]);
		    }
		}
	    }
	    
	  /*
	   * local J_m_P -> J_d_P
	   */
	  var = PRESSURE;
	  if ( pd->v[var] )
	    {
	      pvar = upd->vp[var];
	      for ( j_id=0; j_id<ei->dof[var]; j_id++)
		{			
		  lec->J[peqn_shell][pvar][id_dofshell][j_id] +=
		    scale*(fv->stangent[0][0] * lec->J[peqn_mom1][pvar][id_dofmom1][j_id] +
			   fv->stangent[0][1] * lec->J[peqn_mom2][pvar][id_dofmom2][j_id]);
		}
	    }
	    
	    
	  /*
	   * local J_m_v -> J_d_v  
	   */
	  for ( q=0; q<dim; q++)
	    {
	      var = VELOCITY1+q;
	      if ( pd->v[var] )
		{
		  pvar = upd->vp[var];
		  for ( j_id=0; j_id<ei->dof[var]; j_id++)
		    {
		      lec->J[peqn_shell][pvar][id_dofshell][j_id] +=
			scale*(fv->stangent[0][0] * lec->J[peqn_mom1][pvar][id_dofmom1][j_id] +
			       fv->stangent[0][1] * lec->J[peqn_mom2][pvar][id_dofmom2][j_id]);
		    }
		}
	    }
	    
	}
    } /* end of Jacobian entries */
    
} /* end of routine put_fluid_stress_on_shell */




/*
 * This routine is a specialization of the put_fluid_stress_on_shell.  It is used to 
 * apply the fluid shear stress term onto the shell tension equation.  It is called
 * when the boundary condition SH_SHEAR_STRESS is applied on a boundary and should only 
 * for problems involving just the shell tension unknown.  Problems that involve both 
 * shell structural variables (CURVATURE and TENSION ) should use SH_FLUID_STRESS which applies
 * both normal and tangential stress components.
 * It is written in the spirit of put_fluid_stress_in_solid invoked by the FLUID_SOLID
 * BC. The one exception is that this routine actually performs the rotation of the
 * residual rather than relying on the rotation mechanism in GOMA.  NOTE:  ONLY WRITTEN
 * FOR 2D SIMULATION!!!
 *
 * Appended:  This currently isn't called from anywhere.  The function of this routine was given 
 *     to the boundary condtion SH_FLUID_SHEAR.  This code can be revived quickly enough by making the appropriate
 *     hooks in the code.
 */

void 
put_shear_stress_on_shell(int id, /* local element node number for the 
				   * current node whose residual contribution
				   * is being sought                        */
			  int id_shell, /* local shell element node number corresponding to id */
			  int I, /* Global node number                      */
			  int ielem_dim, /* physical dimension of the elem  */
			  int local_node_list_fs[], /* MDE list to keep track
						     * of nodes at which 
						     * solid contributions 
						     * have been transfered
						     * to liquid (fluid-solid
						     * boundaries)          */
			  double scale) /* Scale factor, nondimension       */
{
  int j_id, dim, var, pvar,  q;
  int id_dofmom1, id_dofmom2, id_dofshell, offset;
  int peqn_mom1, peqn_mom2, peqn_shell;
  int ieqn_mom1, ieqn_mom2, ieqn_shell;
  int tension_fixed;
  NODE_INFO_STRUCT *node = Nodes[I];
  NODAL_VARS_STRUCT *nv = node->Nodal_Vars_Info;

  dim = pd->Num_Dim;

  if (node->DBC) {
     
    offset = get_nodal_unknown_offset(nv, R_SHELL_TENSION, -2, 0, NULL);
    tension_fixed   = (offset >= 0) && (node->DBC[offset] != -1);	
  } else {
    tension_fixed = 0;
  }


  /*
   * if you are in the shell phase, return without doing anything
   * In the shell phase, there are no fluid momentum equations.
   * This probably is not needed since it is impossible with CUBIT to
   * put attach a sideset to a shell element, and hence all sidesets
   * are one-sided. 
   */
  if (!pd->e[R_MOMENTUM1]) return;

  /* Hmm, unlike in the FLUID_SOLID case, we've no shell equations in the 
   * bulk so we cannot use ei struct like this 
   * id_dofmom = ei->ln_to_dof[R_MOMENTUM1][id];
   * id_dofshell = ei->ln_to_dof[R_SHELL_CURVATURE][i];
   *
   * So try this because it is only 2D (assumption that there is a shell
   * variable at every surface node
   */
  id_dofmom1 = ei->ln_to_dof[R_MOMENTUM1][id];
  id_dofshell = id_shell;  /* This is the imperfect assumption */


  /*
   * if this nodal contribution has already been added to fluid momentum
   * equation (i.e. we are at a corner on the second side) return
   * without doing anything
   */
  if (local_node_list_fs[id] == -1) {
    local_node_list_fs[id] = 1;
  } else {
    return;
  }

  /*
   * check to make sure that velocity
   * dofs and shell dofs exist at this node
   */
  if (Dolphin[I][R_MOMENTUM1] <= 0)     return;
  if (Dolphin[I][R_SHELL_TENSION] <= 0) return;

  /*
   * Add local contribution to fluid
   * momentum into local contribution for shell momentum.
   */
  if (af->Assemble_Residual) {
    /* don't do this if the curvature is fixed for this node */

    if (!tension_fixed ) {
      ieqn_shell = R_SHELL_TENSION;
      ieqn_mom1   = R_MOMENTUM1;
      ieqn_mom2   = R_MOMENTUM2;
      id_dofmom1 = ei->ln_to_dof[ieqn_mom1][id];
      id_dofmom2 = ei->ln_to_dof[ieqn_mom2][id];
      id_shell = ei->ln_to_first_dof[ieqn_shell][id];
      /*id_dofshell = id_shell; */

      lec->R[upd->ep[ieqn_shell]][id_dofshell] +=	  scale * (fv->stangent[0][0] * lec->R[upd->ep[ieqn_mom1]][id_dofmom1] + 
								   fv->stangent[0][1] * lec->R[upd->ep[ieqn_mom2]][id_dofmom2]);
		   
    }
  }
    
  /*
   * loop over directions and add local contribution to solid
   * momentum into local contribution for liquid momentum
   */
  if (af->Assemble_Jacobian)
    {


      if (! tension_fixed)
	{
		
	  ieqn_shell = R_SHELL_TENSION;
	  peqn_shell = upd->ep[ieqn_shell];
	  ieqn_mom1   = R_MOMENTUM1;
	  ieqn_mom2   = R_MOMENTUM2;
	  id_dofmom1 = ei->ln_to_dof[ieqn_mom1][id];
	  peqn_mom1 = upd->ep[ieqn_mom1];
	  id_dofmom2 = ei->ln_to_dof[ieqn_mom2][id];
	  peqn_mom2 = upd->ep[ieqn_mom2];
	  id_shell = ei->ln_to_first_dof[ieqn_shell][id];
	  /* id_dofshell = id_shell; */
		
		
	  /* Add contributions due to all nodal sensitivities in shell element */
	    
	  /*
	   * local J_m_d -> J_d_d
	   */
	  for ( q=0; q<dim; q++)
	    {
	      var = MESH_DISPLACEMENT1+q;
	      if ( pd->v[var] )
		{
		  pvar = upd->vp[var];
		  for ( j_id=0; j_id<ei->dof[var]; j_id++)
		    {
		      lec->J[peqn_shell][pvar][id_dofshell][j_id] +=
			scale*(fv->stangent[0][0] * lec->J[peqn_mom1][pvar][id_dofmom1][j_id] +
			       fv->stangent[0][1] * lec->J[peqn_mom2][pvar][id_dofmom2][j_id] +
			       fv->dstangent_dx[0][0][q][j_id] * lec->R[upd->ep[ieqn_mom1]][id_dofmom1] +
			       fv->dstangent_dx[0][1][q][j_id] * lec->R[upd->ep[ieqn_mom2]][id_dofmom2]);
		    }
		}
	    }
	    
	  /*
	   * local J_m_P -> J_d_P
	   */
	  var = PRESSURE;
	  if ( pd->v[var] )
	    {
	      pvar = upd->vp[var];
	      for ( j_id=0; j_id<ei->dof[var]; j_id++)
		{	
		  /*			lec->J[peqn_shell][pvar][id_dofshell][j_id] += scale* lec->J[peqn_mom1][pvar][id_dofmom1][j_id]	;	*/
		  lec->J[peqn_shell][pvar][id_dofshell][j_id] +=
		    scale*(fv->stangent[0][0] * lec->J[peqn_mom1][pvar][id_dofmom1][j_id] +
			   fv->stangent[0][1] * lec->J[peqn_mom2][pvar][id_dofmom2][j_id]);
		}
	    }
	    
	    
	  /*
	   * local J_m_v -> J_d_v  
	   */
	  for ( q=0; q<dim; q++)
	    {
	      var = VELOCITY1+q;
	      if ( pd->v[var] )
		{
		  pvar = upd->vp[var];
		  for ( j_id=0; j_id<ei->dof[var]; j_id++)
		    {
		      /*	lec->J[peqn_shell][pvar][id_dofshell][j_id] += scale* lec->J[peqn_mom1][pvar][id_dofmom1][j_id]; */
		      lec->J[peqn_shell][pvar][id_dofshell][j_id] +=
			scale*(fv->stangent[0][0] * lec->J[peqn_mom1][pvar][id_dofmom1][j_id] +
			       fv->stangent[0][1] * lec->J[peqn_mom2][pvar][id_dofmom2][j_id]);
		    }
		}
	    }
	    
	}
    } /* end of Jacobian entries */
    
} /* end of routine put_shear_stress_on_shell */




/********************************************************************************/

int 
assemble_shell_angle(double time_value,   /* Time */
		     double theta,        /* Time stepping parameter */
		     double delta_t,      /* Time step size */
		     double xi[DIM],      /* Local stu coordinates */
		     const Exo_DB *exo)
{
  int i, j, peqn, var, pvar;
  int p, b;
  double phi_i, phi_j;
  int dim = pd->Num_Dim;
  int eqn = R_SHELL_ANGLE1; 

  /* Declare some neighbor structures */
  int el1 = ei->ielem;
  int el2, nf, err=0;
  int n_dof[MAX_VARIABLE_TYPES];  
  int n_dofptr[MAX_VARIABLE_TYPES][MDE];

  /* Unpack variables from structures for local convenience. */
  double wt    = fv->wt;
  double h3    = fv->h3;
  double det_J = bf[eqn]->detJ; 
  double sn_ang[DIM-1];
  double d_sn_ang_dsn[DIM-1][DIM];
  
  /* See if there is a friend for this element */
  nf = num_elem_friends[el1];
  if (nf == 0) return(err); 

  /*
   * Get neighbor element number.
   * NOTE: For now, only one friend can be handled. If there is
   *       more than one, the first one will be used. This 
   *       will be fixed at a later date.
   */
  el2 = elem_friends[el1][0];
  if (nf != 1) WH(-1, "WARNING: Not set up for more than one element friend!");  

  /*
   * DOF counts will be needed to construct sensitivities of
   * local shell equations to neighbor bulk element variables.
   * This info goes in n_dof
   */

  /*
   * This call will collect all field variable and basis function data
   * for the neighbor element "el2" and populate the fv structure,
   * then repopulate the local element structure. Thus, as long as
   * no variable is defined on both blocks, all data from both
   * elements will be available in the fv structure.
   * For example, pressure in the neighbor element is accessed with fv->P.
   * This is done to simplify the job of writing shell equations.
   */
  err = load_neighbor_var_data(el1, el2, n_dof, NULL, n_dofptr, -1, xi, exo); 
  
  if ( dim == 2 )
    {
      sn_ang[0] = atan2( fv->snormal[1], fv->snormal[0] );
      d_sn_ang_dsn[0][0] = -fv->snormal[1];
      d_sn_ang_dsn[0][1] = fv->snormal[0];
      if ( sn_ang[0] - fv->sh_ang[0] > M_PIE ) sn_ang[0] -= 2. * M_PIE;
      else if ( sn_ang[0] - fv->sh_ang[0] < -M_PIE ) sn_ang[0] += 2. * M_PIE; 
    }
  else if ( dim == 3 )
    {
      double xylen2 = fv->snormal[0] * fv->snormal[0] +
	fv->snormal[1] * fv->snormal[1];
      sn_ang[0] = atan2( fv->snormal[1], fv->snormal[0] );
      sn_ang[1] = acos( fv->snormal[2] );
      if ( xylen2 > 1.e-8 )
        {
	  d_sn_ang_dsn[0][0] = -fv->snormal[1] / xylen2;
          d_sn_ang_dsn[0][1] =  fv->snormal[0] / xylen2;
	  d_sn_ang_dsn[0][2] = 0.;
	  d_sn_ang_dsn[1][0] = 0.;
          d_sn_ang_dsn[1][1] = 0.;
	  d_sn_ang_dsn[1][2] = 1. / sqrt( xylen2 );
	}
      else
        {
	  d_sn_ang_dsn[0][0] = 0.;
          d_sn_ang_dsn[0][1] = 0.;
	  d_sn_ang_dsn[0][2] = 0.;
	  d_sn_ang_dsn[1][0] = 0.;
          d_sn_ang_dsn[1][1] = 0.;
	  d_sn_ang_dsn[1][2] = 0.;
	}
      if ( sn_ang[0] - fv->sh_ang[0] > M_PIE ) sn_ang[0] -= 2. * M_PIE;
      else if ( sn_ang[0] - fv->sh_ang[0] < -M_PIE ) sn_ang[0] += 2. * M_PIE;
      if ( sn_ang[1] - fv->sh_ang[1] > M_PIE ) sn_ang[1] -= 2. * M_PIE;
      else if ( sn_ang[1] - fv->sh_ang[1] < -M_PIE ) sn_ang[1] += 2. * M_PIE;
    }

  /* Assemble the residual equation */
  /* Here, variables from both elements are used. */ 
  if ( af->Assemble_Residual )
    {
      for ( p = 0; p < dim-1; p++ )
        {
	  eqn = R_SHELL_ANGLE1 + p;
	  peqn = upd->ep[eqn];

          for(i = 0; i < ei->dof[eqn]; i++)
	    {
	      phi_i = bf[eqn]->phi[i];
	      
	      lec->R[peqn][i] += ( fv->sh_ang[p] - sn_ang[p] ) * phi_i * wt * det_J * h3;
	    }
	}
    }  

  /* Assemble the sensitivity equations */
  if ( af->Assemble_Jacobian )
    {
      for ( p = 0; p < dim-1; p++ )
        {
	  eqn = R_SHELL_ANGLE1 + p;
	  peqn = upd->ep[eqn];

          for(i = 0; i < ei->dof[eqn]; i++)
	    {
	      phi_i = bf[eqn]->phi[i];
	      
	      for ( b = 0; b < dim-1; b++ )
	        {
		  var = eqn;
		  if ( pd->v[var] )
		    {
	              pvar = upd->vp[var];   

	              for(j = 0; j < ei->dof[var]; j++) 
		        {
			  phi_j = bf[var]->phi[j];
			  
			  lec->J[peqn][pvar][i][j] += phi_j * phi_i * wt * det_J * h3;
		        }
		    }
		}
		
	      /* d_(SHELL_ANGLE)_d_mesh:  REMOTE senstivity */
              /* Here, the unknown is defined on the neighbor element. */
	      for ( b = 0; b < dim; b++ )
	        {
		  var = MESH_DISPLACEMENT1 + b;
		  if ( n_dof[var] > 0)   /* NOTE: Cannot use pd->v here! */
		    {
	              pvar = upd->vp[var];   

                      /* Here, the unknown is defined on the remote element.
                       * So get the DOF count from n_dof. */
	              for(j = 0; j < n_dof[var]; j++) 
		        {
			  lec->J[peqn][pvar][i][j] += d_sn_ang_dsn[p][b] * phi_i *
			    wt * det_J * h3;
		        }
		    }
		}
	    }
	} 
    }
  return(err);
} /* end of assemble_shell_angle */

/********************************************************************************************************************************************/
/* End of file mm_fill_shell.c */
int 
assemble_shell_surface_rheo_pieces(double time_value,   /* Time */
		                   double theta,        /* Time stepping parameter */
		                   double delta_t,      /* Time step size */
		                   double xi[DIM],      /* Local stu coordinates */
		                   const Exo_DB *exo)
{
  int i, j, k, l, m, nn, r, jk;
  int peqn = -1;
  int dofs;
  //! Number of velocities dimensions
  int wim;
  int var, pvar;
  int p, b, a;
  double phi_i, phi_j;
  int dim = pd->Num_Dim;
  int eqn = R_N_DOT_CURL_V; 
  int dof_map[MDE];
#ifdef NEWMETHOD_SHELL_SURF_DIV_V
  double ndotv, tmp_i;
#endif
  double tmp;
  double surfaceDiffusionCoeff = mp->SurfaceDiffusionCoeffProjectionEqn;
  double surfaceDiffusionCoeff2 = mp->SurfaceDiffusionCoeffProjectionEqn;
  double surfaceDiffusionCoeff3 = mp->SurfaceDiffusionCoeffProjectionEqn;
  //double surfaceDiffusionCoeff3 = 0.0;

  /* Declare some neighbor structures */
  int el1 = ei->ielem;
  int el2, nf, err=0;
  /*
   * Number of degrees of freedom in the underlying regular element for each equation type
   */
  int n_dof[MAX_VARIABLE_TYPES];  
  int n_dofptr[MAX_VARIABLE_TYPES][MDE];

  /* Unpack variables from structures and define variables for local convenience. */
  double wt    = fv->wt;
  double h3    = fv->h3;
  /*
   * Find the variable whose basis functions are used to do element integrations of shape
   * functions for the shell equation
   */
  int shapeVar = pd->ShapeVar;
  double det_J = bf[shapeVar]->detJ;
 
  //! Value of the surface divergence of the velocity, evaluated at the gauss point
  double div_s_V = 0.0;

  double grad_v_dot_n[DIM];
  double grad_s_v_dot_n[DIM];
  double n_dot_curl_v = 0.0;
  double nn_dot_curl_v = 0.0;
  double n_dot_curl_s_v = 0.0;
  double curv = 0.0;
  double d_div_s_V_dx[DIM][MDE], d_div_s_V_dv[DIM][MDE];
  double d_grad_v_dot_n_dx[DIM][DIM][MDE], d_grad_v_dot_n_dv[DIM][DIM][MDE];
  double d_grad_s_v_dot_n_dx[DIM][DIM][MDE], d_grad_s_v_dot_n_dv[DIM][DIM][MDE];
  double d_n_dot_curl_s_v_dx[DIM][MDE], d_n_dot_curl_s_v_dv[DIM][MDE];
  double d_nn_dot_curl_v_dx[DIM][MDE], d_nn_dot_curl_v_dv[DIM][MDE];
  double d_n_dot_curl_v_dx[DIM][MDE], d_n_dot_curl_v_dv[DIM][MDE];
  double d_div_s_n_dx[DIM][MDE];
  double d_div_s_n_dn[DIM][MDE];
  double d_grad_n_dn[DIM][DIM][DIM][MDE];
  BASIS_FUNCTIONS_STRUCT *bfn;
  BASIS_FUNCTIONS_STRUCT *bfv;

  if (CURL_V == -1) {
    EH(-1,"ERROR: inconsistency: need to set Vorticity Vector = yes in Post processing section");
  }

  // Calculate the number of velocity dimensions
  wim   = dim;
  if (pd->CoordinateSystem == SWIRLING || pd->CoordinateSystem == PROJECTED_CARTESIAN) 
    {
      wim = 3;
    }
  
  /* See if there is a friend for this element */
  nf = num_elem_friends[el1];
  if (nf == 0) return(err); 

  /*
   * Get neighbor element number.
   * NOTE: For now, only one friend can be handled. If there is
   *       more than one, the first one will be used. This 
   *       will be fixed at a later date.
   */
  el2 = elem_friends[el1][0];
  if (nf != 1) WH(-1, "WARNING: Not set up for more than one element friend!");  
  memset(d_grad_n_dn, 0, sizeof(double)*DIM*DIM*DIM*MDE);
  if (pd->v[SHELL_NORMAL1]) 
    {
      eqn = SHELL_NORMAL1;
      dofs = ei->dof[eqn];
      bfn =  bf[eqn];
      for (i = 0; i < VIM; i++)
	{
	  for (j = 0; j < VIM; j++)
	    {
	      for (r = 0; r < dim; r++)
		{
		  for (k = 0; k < dofs; k++)
		    {
		      d_grad_n_dn[i][j] [r][k] = bfn->grad_phi_e[k][r][i][j];
		    }
		}
	    }
	}  
    }

  /*
   * DOF counts will be needed to construct sensitivities of
   * local shell equations to neighbor bulk element variables.
   * This info goes in n_dof
   */

  /*
   * This call will collect all field variable and basis function data
   * for the neighbor element "el2" and populate the fv structure,
   * then repopulate the local element structure. Thus, as long as
   * no variable is defined on both blocks, all data from both
   * elements will be available in the fv structure.
   * For example, pressure in the neighbor element is accessed with fv->P.
   * This is done to simplify the job of writing shell equations.
   */
  err = load_neighbor_var_data(el1, el2, n_dof, dof_map, n_dofptr, -1, xi, exo); 
  

  if (dim == 3)
    {
      printf("WARNING: assemble_shell_surface_rheo_pieces() Shell capability in 3D is untested\n");
    }

  /* Form components with surface gradient operator */

  div_s_V = 0.0;
  curv = n_dot_curl_v = 0.0;
  memset(d_div_s_V_dx, 0, sizeof(double)*DIM*MDE);
  memset(d_div_s_V_dv, 0, sizeof(double)*DIM*MDE);
  memset(d_div_s_n_dx, 0, sizeof(double)*DIM*MDE);
  memset(d_div_s_n_dn, 0, sizeof(double)*DIM*MDE);


  // HKM -> Changed below from dim to VIM
  bfv = bf[VELOCITY1];

  for (i = 0; i < VIM; i++)
    {
      for (j = 0; j < VIM; j++)
	{
	  // Form grad_s dot v:
	  //            note Grad_s = (I - n n) dot grad
	  //                   fv->grad_v[j][i]  -> d (v[i]) / d (x[j])
	  div_s_V += (delta(i,j) * fv->grad_v[j][i] - fv->snormal[i] * fv->snormal[j] * fv->grad_v[j][i]);
	      
	  // Form the mesh dependencies for div_s_V
	  for (b = 0; b < dim; b++)
	    {
	      for (k = 0; k < n_dof[MESH_DISPLACEMENT1]; k++)
		{
		  d_div_s_V_dx[b][k] += (delta(i,j) * fv->d_grad_v_dmesh[i][j][b][k]
					 - fv->dsnormal_dx[i][b][k] * fv->snormal[j] * fv->grad_v[j][i]
					 - fv->snormal[i] * fv->dsnormal_dx[j][b][k] * fv->grad_v[j][i]
					 - fv->snormal[i] * fv->snormal[j] * fv->d_grad_v_dmesh[j][i][b][k]);
		}
	    }

	  // Form the velocity dependencies 
	  for (b = 0; b < wim; b++)
	    {
	      for (k = 0; k < n_dof[VELOCITY1]; k++)
		{
		  d_div_s_V_dv[b][k] += (bfv->grad_phi_e[k][b][i][j] * delta(i,j) -
					 fv->snormal[i] * fv->snormal[j] * bfv->grad_phi_e[k][b][i][j]);
		}
	    }
   
	  /*
	   *   This equation is correct according to the numerical debugger. It needs some explanation. 
	   *   We have formulated fv->d_grad_n_dmesh using the local dofs for the shell equations.
	   *   However, the mesh equations are being accumulated using the column dofs corresponding
	   *   to the bulk mesh unknowns. This is very convoluted. Therefore, we nned to use the
	   *   dof_map[] variable to adjust the column variables from the shell positions into the
	   *   bulk positions.
	   */
	  for (b = 0; b < dim; b++)
	    {
	      for (k = 0; k < ei->dof[MESH_DISPLACEMENT1]; k++)
		{
		  jk = dof_map[k];
		  d_div_s_n_dx[b][jk] += (delta(i,j) * fv->d_grad_n_dmesh[i][j][b][k]
					  - fv->n[i] * fv->n[j] * fv->d_grad_n_dmesh[j][i][b][k]);
		}		  
	    }
	}
    }
  // HKM -> SHELL_NORMAL jacobian terms now work.
  eqn = R_SHELL_NORMAL1;
 
  for (i = 0; i < VIM; i++)
    {
      for (j = 0; j < VIM; j++)
	{
	  for (b = 0; b < dim; b++)
	    {
	      for (k = 0; k < ei->dof[eqn]; k++)
		{		     
		  d_div_s_n_dn[b][k] += (d_grad_n_dn[i][j] [b][k] * delta(i,j) -
					 fv->n[i] * fv->n[j] * d_grad_n_dn[j][i] [b][k] );
		}
	    }
	}
    }
   
  for (i = 0; i < VIM; i++)
    {
      for (j = 0; j < VIM; j++)
	{    
	  for (k = 0; k < ei->dof[eqn]; k++)
	    {  
	      phi_j = bf[eqn]->phi[k];		 
	      d_div_s_n_dn[i][k] +=  - phi_j * fv->n[j] * fv->grad_n[j][i];
	    }
	}
    }

  for (i = 0; i < VIM; i++)
    {
      for (j = 0; j < VIM; j++)
	{
	  for (k = 0; k < ei->dof[eqn]; k++)
	    {
	      phi_j = bf[eqn]->phi[k];
	      d_div_s_n_dn[j][k] +=  - fv->n[i] * phi_j * fv->grad_n[j][i] ;
	    }
	}
    }

  /*
   *  HKM -> experimental section to debug
   */
  // shell_normal_div_s(&div_s_nv_alt, d_div_s_nv_dnv_alt, d_div_s_nv_dmesh_alt);
  
  curv = - 0.5 * (fv->div_s_n);

  /* n . (curl_s_v) */
  /*  = n . ((I - nn).curl_v = n . (curl_v - n(n.grad)Xv) */
  memset(d_n_dot_curl_s_v_dx, 0.0, sizeof(double)*DIM*MDE);
  memset(d_n_dot_curl_v_dx,   0.0, sizeof(double)*DIM*MDE);
  memset(d_nn_dot_curl_v_dx,  0.0, sizeof(double)*DIM*MDE);
  memset(d_n_dot_curl_s_v_dv, 0.0, sizeof(double)*DIM*MDE);
  memset(d_n_dot_curl_v_dv,   0.0, sizeof(double)*DIM*MDE);
  memset(d_nn_dot_curl_v_dv,  0.0, sizeof(double)*DIM*MDE);


  // dim is ok here because snormal will always be zero in the third dimension for cylindrical coordinate
  n_dot_curl_v   = 0.0;
  for (i = 0; i < dim; i++)
    {
      n_dot_curl_v += (fv->curl_v[i] * fv->snormal[i]); 
    }

  bfv = bf[VELOCITY1];
  nn_dot_curl_v  = 0.0;
  for (i = 0; i < VIM; i++)
    {
      for (j = 0; j < dim; j++)
	{
	  for (k = 0; k < VIM; k++)
	    {
	      for (nn = 0; nn < dim; nn++)
		{
		  nn_dot_curl_v += fv->snormal[j] * fv->snormal[nn] * permute(j, k, i) * fv->snormal[i] * fv->grad_v[nn][k];

		  for (l = 0; l < dim; l++)
		    {
		      for (m = 0; m < n_dof[MESH_DISPLACEMENT1]; m++)
			{
			  d_nn_dot_curl_v_dx[l][m] += permute(j,k,i) *
			    (fv->dsnormal_dx[j][l][m]*fv->snormal[nn]* fv->snormal[i] * fv->grad_v[nn][k]
			     + fv->snormal[j] * fv->dsnormal_dx[nn][l][m] * fv->snormal[i] * fv->grad_v[nn][k]
			     + fv->snormal[j]*fv->snormal[nn]*fv->dsnormal_dx[i][l][m] * fv->grad_v[nn][k]
			     + fv->snormal[j]*fv->snormal[nn]* fv->snormal[i] * fv->d_grad_v_dmesh[nn][k][l][m]);
			}
		    }
		  for (l = 0; l < wim; l++)
		    { 
		      for (m = 0; m < n_dof[VELOCITY1]; m++)
			{
			  d_nn_dot_curl_v_dv[l][m] += permute(j,k,i) *
			    (fv->snormal[j]*fv->snormal[nn]* fv->snormal[i] * bfv->grad_phi_e[m][l][nn][k]);
			}
		    }
		}
	    }
	}
    }

  n_dot_curl_s_v = n_dot_curl_v - nn_dot_curl_v;

  /* Jacobian pieces for later */
  for (i = 0; i < dim; i++)
    {
      for (b = 0; b < dim; b++)
	{
	  for (j = 0; j < n_dof[MESH_DISPLACEMENT1]; j++)
	    {
	      d_n_dot_curl_v_dx[b][j] += ( fv->curl_v[i] * fv->dsnormal_dx[i][b][j]);  

	      for (k = 0; k < VIM; k++)
		{
		  for (nn = 0; nn < VIM; nn++)
		    {		
		      d_n_dot_curl_v_dx[b][j] +=  fv->snormal[i] * 
			permute(i,nn,k) * fv->d_grad_v_dmesh[nn][k][b][j];
		    }
		}
	    }
	}
      for (b = 0; b < wim; b++)
	{
	  for (j = 0; j < n_dof[VELOCITY1]; j++)
	    {
	      d_n_dot_curl_v_dv[b][j] += (bfv->curl_phi_e[j][b][i] * fv->snormal[i]);
	    }
	}
    }

#ifdef DEBUG_HKM
  double tester[3];
  tester[0] = 0.0;
  tester[1] = 0.0;
  tester[2] = 0.0;
  for (i = 0; i < VIM; i++)
    {
      for (k = 0; k < VIM; k++)
	{
	  for (nn = 0; nn < VIM; nn++)
	    {
	      tester[i] += permute(i,nn,k) *  fv->grad_v[nn][k];
	    }
	}
      if (fabs(tester[i] - fv->curl_v[i]) > 1.0E-9) {
	printf("WARNING %d: assertion on curl_v is wrong, permute grad [%d] = %g, curl_v[%d] = %g",
               ProcID, i, tester[i], i, fv->curl_v[i]);
      }
    }
#endif

    
  for (b = 0; b < dim; b++)
    {
      for (j = 0; j < n_dof[MESH_DISPLACEMENT1]; j++)
	{
	  d_n_dot_curl_s_v_dx[b][j] = d_n_dot_curl_v_dx[b][j] - d_nn_dot_curl_v_dx[b][j];
	}
    }
  for (b = 0; b < wim; b++)
    {
      for (j = 0; j < n_dof[VELOCITY1]; j++)
	{
	  d_n_dot_curl_s_v_dv[b][j] = d_n_dot_curl_v_dv[b][j] - d_nn_dot_curl_v_dv[b][j];
	}
    }
      

  memset(grad_v_dot_n, 0, sizeof(double)*DIM);
  memset(grad_s_v_dot_n, 0, sizeof(double)*DIM);
  memset(d_grad_v_dot_n_dx, 0, sizeof(double)*DIM*MDE*DIM);
  memset(d_grad_v_dot_n_dv, 0, sizeof(double)*DIM*MDE*DIM);
  memset(d_grad_s_v_dot_n_dv, 0, sizeof(double)*DIM*MDE*DIM);
  memset(d_grad_s_v_dot_n_dx, 0, sizeof(double)*DIM*MDE*DIM);
  for (i = 0; i < VIM; i++)
    {
      for (j = 0; j < dim; j++)
	{
	  grad_v_dot_n[i] += (fv->grad_v[i][j] * fv->snormal[j]);
	  for (b = 0; b < dim; b++)
	    {
	      for (k = 0; k < n_dof[MESH_DISPLACEMENT1]; k++)
		{
		  d_grad_v_dot_n_dx[i][b][k] += fv->d_grad_v_dmesh[i][j][b][k] * fv->snormal[j] + fv->grad_v[i][j] * fv->dsnormal_dx[j][b][k];
		}
	    }
	  for (b = 0; b < wim; b++)
	    {
	      var = VELOCITY1 + b;
	      for (k = 0; k < n_dof[var]; k++)
		{
		  d_grad_v_dot_n_dv[i][b][k] += bf[var]->grad_phi_e[k][b][i][j] * fv->snormal[j];
		}
	    }

	}
    }


  for (i = 0; i < VIM; i++)
    {
      grad_s_v_dot_n[i] = grad_v_dot_n[i];
      for (j = 0; j < VIM; j++)
	{
	  grad_s_v_dot_n[i] -= fv->snormal[i] * fv->snormal[j] * grad_v_dot_n[j];
	  
	  for (b = 0; b < wim; b++)
	    {
	      var = VELOCITY1 + b;
	      for (k = 0; k < n_dof[var]; k++)
		{
		  d_grad_s_v_dot_n_dv[i][b][k] += delta(i,j) * d_grad_v_dot_n_dv[i][b][k] - fv->snormal[i] * fv->snormal[j] * d_grad_v_dot_n_dv[j][b][k];
		}
	    }

	  for (b = 0; b < dim; b++)
	    {
	      for (k = 0; k < n_dof[MESH_DISPLACEMENT1]; k++)
		{
		  d_grad_s_v_dot_n_dx[i][b][k] += delta(i,j) * d_grad_v_dot_n_dx[i][b][k] - fv->snormal[i] * fv->snormal[j] * d_grad_v_dot_n_dx[j][b][k];
		  d_grad_s_v_dot_n_dx[i][b][k] -= fv->dsnormal_dx[i][b][k] * fv->snormal[j] * grad_v_dot_n[j];
		  d_grad_s_v_dot_n_dx[i][b][k] -= fv->snormal[i] * fv->dsnormal_dx[j][b][k] * grad_v_dot_n[j];
		}
	    }
	  
	  
	}
    }




  /* Assemble the residual equation */
  /* Here, variables from both elements are used. */ 
  if (af->Assemble_Residual)
    {
      /*
       *  Surface Divergence of the velocity
       *       div_s_V  -> Previously calculated at the top of this function
       *       fv_div_s_v -> This is the variable whose gp value has been evalulated in ?routine?
       */
      eqn = R_SHELL_SURF_DIV_V;
      peqn = upd->ep[eqn];
      bfn = bf[eqn];
#ifdef NEWMETHOD_SHELL_SURF_DIV_V
#ifdef MASSLUMP_PROJECTION_EQUATIONS
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  phi_i = bfn->phi[i];
	  lec->R[peqn][i] += ( *(esp->div_s_v[i])) * phi_i * wt * det_J * h3;
	  tmp_i = 0.0;
	  ndotv = 0.0;
	  for (b = 0; b < dim; b++) 
	    {
	      ndotv += fv->n[b] * fv->v[b];
	      tmp_i += fv->v[b] * bfn->grad_phi[i][b];
	      for (a = 0; a < dim; a++) 
		{
		  tmp_i -= fv->n[b] * fv->n[a] * fv->v[a] * bfn->grad_phi[i][b];
		}
	    }
	  lec->R[peqn][i] += (2.0 * fv->curv * ndotv) * phi_i * wt * det_J * h3;	 
	  lec->R[peqn][i] += tmp_i * wt * det_J * h3;
	}
#else
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  phi_i = bf[eqn]->phi[i];
	  lec->R[peqn][i] += (fv->div_s_v ) * phi_i * wt * det_J * h3;
	  tmp_i = 0.0;
	  ndotv = 0.0;
	  for (b = 0; b < dim; b++) 
	    {
	      ndotv += fv->n[b] * fv->v[b];
	      tmp_i += fv->v[b] * bfn->grad_phi[i][b];
	      for (a = 0; a < dim; a++) 
		{
		  tmp_i -= fv->n[b] * fv->n[a] * fv->v[a] * bfn->grad_phi[i][b];
		}
	    }
	  lec->R[peqn][i] += (2.0 * fv->curv * ndotv) * phi_i * wt * det_J * h3;
	  lec->R[peqn][i] += tmp_i * wt * det_J * h3;
	}

#endif
#else
#ifdef MASSLUMP_PROJECTION_EQUATIONS
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  phi_i = bf[eqn]->phi[i];
	  // div_s_v_ml = *(esp->div_s_v[i]);
	  lec->R[peqn][i] += ( *(esp->div_s_v[i]) - div_s_V) * phi_i * wt * det_J * h3;
	}
#else
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  phi_i = bf[eqn]->phi[i];
	  lec->R[peqn][i] += (fv->div_s_v - div_s_V) * phi_i * wt * det_J * h3;
	}
#endif
#endif
      /*
       *  Add surface diffusion term
       */
      for (i = 0; i < ei->dof[eqn]; i++)
	{ 
	  tmp = 0.0;
	  for (b = 0; b < dim; b++)
	    {
	      tmp += fv->n[b] *  bfn->grad_phi[i][b];
	    } 
	  for (b = 0; b < dim; b++) 
	    {
	      lec->R[peqn][i] -= surfaceDiffusionCoeff * (fv->grad_div_s_v[b] * (bfn->grad_phi[i][b] - tmp * fv->n[b])) * wt * det_J * h3;
	    }
	}


      eqn = R_SHELL_SURF_CURV;
      peqn = upd->ep[eqn];
      bfn = bf[eqn];
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  phi_i = bf[eqn]->phi[i];
       	  lec->R[peqn][i] += (fv->curv - curv) * phi_i * wt * det_J * h3;

	}
      /*
       *  Add surface diffusion term
       */
      for (i = 0; i < ei->dof[eqn]; i++)
	{ 
	  tmp = 0.0;
	  for (b = 0; b < dim; b++)
	    {
	      tmp += fv->n[b] * bfn->grad_phi[i][b];
	    } 
	  for (b = 0; b < dim; b++) 
	    {
	      lec->R[peqn][i] -= surfaceDiffusionCoeff3 * (fv->grad_curv[b] * (bfn->grad_phi[i][b] - tmp * fv->n[b]))  * wt * det_J * h3;
	    }
	}


#ifdef DEBUG_HKM
      /*
       *   Here we assert a relationship that is important. the surface determinant, which is
       *   calculated from the surface integral of the bulk parent element, is equal to the
       *   determinant of the shell jacobian * h3 of the shell element. 
       *   
       */
      if (fabs(det_J * h3 - fv->sdet) > 1.0E-9) {
	printf("we have a problem: assertion on sdet\n");
	exit(-1);
      }
      // if ( fabs(n_dot_curl_s_v) > 1.0E-15) {
      //	printf("we should not be here for non-swirling flow\n");
      // }
#endif
      eqn = R_N_DOT_CURL_V;
      peqn = upd->ep[eqn];

      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  phi_i = bf[eqn]->phi[i];	      
	
	  lec->R[peqn][i] += (fv->n_dot_curl_s_v - n_dot_curl_s_v) * phi_i * wt * det_J * h3; 
	}


      /*
       *  GRAD_S_V_DOT_N1 terms
       */
      for (p = 0; p < dim; p++)
	{
	  eqn = R_GRAD_S_V_DOT_N1 + p;
	  peqn = upd->ep[eqn];
	  bfn = bf[eqn];
	  for (i = 0; i < ei->dof[eqn]; i++)
	    {
	      phi_i = bf[eqn]->phi[i];
	      lec->R[peqn][i] += (fv->grad_v_dot_n[p] - grad_s_v_dot_n[p]) * phi_i * wt * det_J * h3;
	      
	      /*
	       *  Add surface diffusion term
	       */

	      tmp = 0.0;
	      for (b = 0; b < dim; b++)
		{
		  tmp += fv->n[b] *  bfn->grad_phi[i][b];
		} 
	      for (b = 0; b < VIM; b++) 
		{
		  lec->R[peqn][i] -= surfaceDiffusionCoeff2 * (fv->serialgrad_grad_s_v_dot_n[p][b] * (bfn->grad_phi[i][b] - tmp * fv->n[b])) * wt * det_J * h3;
		}	      
	    }
	}
 
    }
  /* Assemble the sensitivity equations */
  if (af->Assemble_Jacobian)
    {
      for (p = 0; p < dim; p++ )
	{
	  eqn = R_GRAD_S_V_DOT_N1 + p;
	  peqn = upd->ep[eqn];
	  bfn = bf[eqn];
	  // lec->R[peqn][i] += ( fv->grad_v_dot_n[p] - grad_s_v_dot_n[p] ) * phi_i * wt * det_J * h3;
	  for (i = 0; i < ei->dof[eqn]; i++)
	    {
	      phi_i = bf[eqn]->phi[i];
	      
	      var = eqn;
	      if (pd->v[var])
		{
		  pvar = upd->vp[var];   

		  for (j = 0; j < ei->dof[var]; j++) 
		    {
		      phi_j = bf[var]->phi[j];			  
		      lec->J[peqn][pvar][i][j] += phi_j * phi_i * wt * det_J * h3;
		    }
		}
		
	      for (b = 0; b < dim; b++ )
		{
		  var = MESH_DISPLACEMENT1 + b;
		  if (n_dof[var] > 0)   /* NOTE: Cannot use pd->v here! */
		    {
		      pvar = upd->vp[var];   

		      /* Here, the unknown is defined on the remote element.
		       * So get the DOF count from n_dof. */
		      for (j = 0; j < n_dof[var]; j++) 
			{
			  lec->J[peqn][pvar][i][j] -= d_grad_s_v_dot_n_dx[p][b][j] * phi_i * wt * det_J * h3;
			}
		      for (j = 0; j < n_dof[var]; j++) 
			{
			  lec->J[peqn][pvar][i][j] += (fv->grad_v_dot_n[p] - grad_s_v_dot_n[p]) * phi_i * wt * (fv->dsurfdet_dx[b][j]);
			}
		    }
		}
	      for (b = 0; b < VIM; b++ )
		{
		  var = VELOCITY1 + b;
		  if (n_dof[var] > 0)   /* NOTE: Cannot use pd->v here! */
		    {
		      pvar = upd->vp[var];   

		      /* Here, the unknown is defined on the remote element.
		       * So get the DOF count from n_dof. */
		      for (j = 0; j < n_dof[var]; j++) 
			{
			  lec->J[peqn][pvar][i][j] -= d_grad_s_v_dot_n_dv[p][b][j] * phi_i * wt * det_J * h3;
			}
		    }
		}
	    
	
	      /*
	       * Add the surface diffusion term jacobians
	       */
	  
	      tmp = 0.0;
	      for (a = 0; a < dim; a++)
		{
		  tmp += fv->n[a] * bfn->grad_phi[i][a];
		}
	
	      var = eqn;
	      if (pd->v[var])
		{
		  pvar = upd->vp[var];
		  for (b = 0; b < VIM; b++)
		    {
		      for (j = 0; j < ei->dof[var]; j++)
			{			  
			  //      lec->R[peqn][i] -= surfaceDiffusionCoeff2 * (fv->serialgrad_grad_s_v_dot_n[p][b] * (bfn->grad_phi[i][b] - tmp * fv->n[b])) * wt * det_J * h3;
			  lec->J[peqn][pvar][i][j] -= surfaceDiffusionCoeff2 * bfn->grad_phi[j][b] * (bfn->grad_phi[i][b] - tmp * fv->n[b])  * wt * det_J * h3;
			}
		    }
		}

	      for (b = 0; b < dim; b++)
		{
		  var = SHELL_NORMAL1 + b;
		  if (pd->v[var])
		    {
		      pvar = upd->vp[var];
		      for (j = 0; j < ei->dof[var]; j++)
			{
			  phi_j = bf[var]->phi[j];
			  lec->J[peqn][pvar][i][j] -= surfaceDiffusionCoeff2 * fv->serialgrad_grad_s_v_dot_n[p][b] * (- tmp * phi_j)  * wt * det_J * h3;
		      
			  for (a = 0; a < dim; a++)
			    {
			      lec->J[peqn][pvar][i][j] -= surfaceDiffusionCoeff2 * fv->serialgrad_grad_s_v_dot_n[p][a] * (- phi_j * bfn->grad_phi[i][b] * fv->n[a])  * wt * det_J * h3;
			    } 
			}
		    }
		}

	      for (r = 0; r < dim; r++ )
		{
		  var = MESH_DISPLACEMENT1 + r;
		  if (n_dof[var] > 0)
		    {
		      pvar = upd->vp[var];   
	
		      for (j = 0; j < n_dof[var]; j++) 
			{
			  for (b = 0; b < dim; b++)
			    {
			      lec->J[peqn][pvar][i][j] -= surfaceDiffusionCoeff2 * fv->serialgrad_grad_s_v_dot_n[p][b] * (bfn->grad_phi[i][b] - tmp * fv->n[b]) * wt * (fv->dsurfdet_dx[r][j]);	
			    }
			}
	
		      for (j = 0; j < ei->dof[var]; j++)
			{
			  jk = dof_map[j];
			  for (b = 0; b < dim; b++) 
			    {
			      lec->J[peqn][pvar][i][jk] -= surfaceDiffusionCoeff2 * fv->serialgrad_grad_s_v_dot_n[p][b] * (bfn->d_grad_phi_dmesh[i][b][r][j]) * wt *  det_J * h3;
			      lec->J[peqn][pvar][i][jk] -= surfaceDiffusionCoeff2 * fv->d_serialgrad_grad_s_v_dot_n_dmesh[p][b][r][j] * (bfn->grad_phi[i][b] - tmp * fv->n[b]) * wt * det_J * h3;
			      for (a = 0; a < dim; a++) 
				{
				  lec->J[peqn][pvar][i][jk] -= surfaceDiffusionCoeff2 * fv->serialgrad_grad_s_v_dot_n[p][b] * (- fv->n[a] * bfn->d_grad_phi_dmesh[i][a][r][j] * fv->n[b]) * wt * det_J * h3;
				}
			    }
			}
		    }
		}
	    }
	}

     
      eqn = R_SHELL_SURF_DIV_V;
      peqn = upd->ep[eqn];
      bfn = bf[eqn];
      //      lec->R[peqn][i] += (fv->div_s_v - div_s_V) * phi_i * wt * det_J * h3;
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  phi_i = bfn->phi[i];
#ifdef NEWMETHOD_SHELL_SURF_DIV_V
	  tmp_i = 0.0;
	  ndotv = 0.0;
	  for (b = 0; b < dim; b++) 
	    {
	      ndotv += fv->n[b] * fv->v[b];
	      tmp_i += fv->v[b] * bfn->grad_phi[i][b];
	      for (a = 0; a < dim; a++) 
		{
		  tmp_i -= fv->n[b] * fv->n[a] * fv->v[a] * bfn->grad_phi[i][b];
		}
	    }
#endif
	      
	  var = eqn;
	  if (pd->v[var])
	    {
	      pvar = upd->vp[var];
#ifdef MASSLUMP_PROJECTION_EQUATIONS
	      lec->J[peqn][pvar][i][i] += phi_i * wt * det_J * h3;
#else
	      for (j = 0; j < ei->dof[var]; j++) 
		{
		  phi_j = bf[var]->phi[j];			  
		  lec->J[peqn][pvar][i][j] += phi_j * phi_i * wt * det_J * h3;
		}
#endif
	    }
#ifdef NEWMETHOD_SHELL_SURF_DIV_V
	  //     lec->R[peqn][i] += (2.0 * fv->curv * fv->n[b] * fv->v[b]) * phi_i * wt * det_J * h3;
	  for (b = 0; b < dim; b++)
	    {
	      var = VELOCITY1 + b;
	      if (n_dof[var] > 0)   /* NOTE: Cannot use pd->v here! */
		{
		  pvar = upd->vp[var];
		  /* Here, the unknown is defined on the remote element.
		   * So get the DOF count from n_dof. */
		  for (j = 0; j < n_dof[var]; j++) 
		    {
		      phi_j = bfv->phi[j];
		      lec->J[peqn][pvar][i][j] += (2.0 * fv->curv * fv->n[b] * phi_j) * phi_i * wt * det_J * h3;
		      lec->J[peqn][pvar][i][j] += phi_j * bfn->grad_phi[i][b] * wt * det_J * h3;
		      for (a = 0; a < dim; a++) 
			{ 
			  lec->J[peqn][pvar][i][j] -= fv->n[a] * fv->n[b] * phi_j * bfn->grad_phi[i][a] * wt * det_J * h3;
			}
		    }
		}
	    }

	  for (b = 0; b < dim; b++)
	    {
	      var = SHELL_NORMAL1 + b;
	      if (pd->v[var])
		{
		  pvar = upd->vp[var];
		  for (j = 0; j < ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];
		      lec->J[peqn][pvar][i][j] += (2.0 * fv->curv * phi_j * fv->v[b]) * phi_i * wt * det_J * h3;
		      for (a = 0; a < dim; a++) 
			{
			  lec->J[peqn][pvar][i][j] -= phi_j * fv->n[a] *(fv->v[a] * bfn->grad_phi[i][b] + fv->v[b] * bfn->grad_phi[i][a]) * wt * det_J * h3;
			}
		    }
		}
	    }

	  var = SHELL_SURF_CURV;
	  if (pd->v[var]) 
	    {
	      for (b = 0; b < dim; b++)
		{
		  pvar = upd->vp[var];
		  for (j = 0; j < ei->dof[var]; j++) 
		    {
		      phi_j = bf[var]->phi[j];
		      lec->J[peqn][pvar][i][j] += (2.0 * phi_j * fv->n[b] * fv->v[b]) * phi_i * wt * det_J * h3;
		    }
		}
	    }

	  ndotv = 0.0;
	  for (b = 0; b < dim; b++ )
	    {
	      ndotv += fv->n[b] * fv->v[b];
	    }
	  for (p = 0; p < dim; p++ )
	    {
	      var = MESH_DISPLACEMENT1 + p;
	      if (n_dof[var] > 0)   /* NOTE: Cannot use pd->v here! */
		{
		  pvar = upd->vp[var];   
		  /* Here, the unknown is defined on the remote element.
		   * So get the DOF count from n_dof. */
		  for (j = 0; j < n_dof[var]; j++) 
		    {
#ifdef MASSLUMP_PROJECTION_EQUATIONS
		      tmp = *(esp->div_s_v[i]);
#else
		      tmp = fv->div_s_v;
#endif
		      tmp += (2.0 * fv->curv * ndotv);
		      lec->J[peqn][pvar][i][j] += (tmp_i + tmp * phi_i) * wt * (fv->dsurfdet_dx[p][j]);
		    }
		  /*
		   *  This has been checked out in the debugger and it passes. What's going on?
		   *  The dgrad_phi_dmesh is the derivative of the gradients of the shell basis function for the shell variable with
		   *  respect to the mesh variables. Now the mesh variables are delineated by the bulk element. However, all derivatives
		   *  wrt unknowns of the mesh variables that are not on the shell must be zero (yes this is true). 
		   *  What's actually in dgrad_phi_dmesh is the variation of the gradient of the shell variables with respect to the
		   *  displacement unknowns for nodes which are on the shell (in the shell local node number convention). We then
		   *  use dof_map[] to change the shell local node number convention to the bulk element local node convention.
		   *  this then gets put in the correct place in load_lec().  
		   *  We use ei->dof[var] here because it correctly identifies the number of mesh unknowns that are on the shell (3)
		   *  even though mesh unknowns are only defined on the bulk elements.
		   */
		  for (j = 0; j < ei->dof[var]; j++)
		    {
		      jk = dof_map[j];
		      for (b = 0; b < dim; b++) 
			{
			  lec->J[peqn][pvar][i][jk] += fv->v[b] * bfn->d_grad_phi_dmesh[i][b][p][j] * wt * det_J * h3;
			  for (a = 0; a < dim; a++) 
			    {
			      lec->J[peqn][pvar][i][jk] -= fv->n[b] * fv->n[a] * fv->v[a] * bfn->d_grad_phi_dmesh[i][b][p][j] * wt * det_J * h3;
			    }
			}
		    }
		}
	    }

#else // NEW
	  for (b = 0; b < dim; b++ )
	    {
	      var = MESH_DISPLACEMENT1 + b;
	      if (n_dof[var] > 0)   /* NOTE: Cannot use pd->v here! */
		{
		  pvar = upd->vp[var];   

		  /* Here, the unknown is defined on the remote element.
		   * So get the DOF count from n_dof. */
		  for (j = 0; j < n_dof[var]; j++) 
		    {
		      lec->J[peqn][pvar][i][j] -= d_div_s_V_dx[b][j] * phi_i * wt * det_J * h3;
		    }
#ifdef MASSLUMP_PROJECTION_EQUATIONS
		  for (j = 0; j < n_dof[var]; j++) 
		    {
		      lec->J[peqn][pvar][i][j] += (*(esp->div_s_v[i]) - div_s_V) * phi_i * wt * (fv->dsurfdet_dx[b][j]);
		    }
#else
		  for (j = 0; j < n_dof[var]; j++) 
		    {
		      lec->J[peqn][pvar][i][j] += (fv->div_s_v - div_s_V) * phi_i * wt * (fv->dsurfdet_dx[b][j]);
		    }
#endif
		}
	    }
	  for (b = 0; b < VIM; b++)
	    {
	      var = VELOCITY1 + b;
	      if (n_dof[var] > 0)   /* NOTE: Cannot use pd->v here! */
		{
		  pvar = upd->vp[var];   

		  /* Here, the unknown is defined on the remote element.
		   * So get the DOF count from n_dof. */
		  for (j = 0; j < n_dof[var]; j++) 
		    {
		      lec->J[peqn][pvar][i][j] -= d_div_s_V_dv[b][j] * phi_i * wt * det_J * h3;
		    }
		}
	    }
#endif


	  /*
	   *  Add surface diffusion term
	   */
	
	  tmp = 0.0;
	  for (a = 0; a < dim; a++)
	    {
	      tmp += fv->n[a] * bfn->grad_phi[i][a];
	    } 

       	  //   lec->R[peqn][i] -= surfaceDiffusionCoeff * (fv->grad_div_s_v[b] * (bfn->grad_phi[i][b] - tmp * fv->n[b]))  * wt * det_J * h3;
	
	  var = eqn;
	  if (pd->v[var])
	    {
	      pvar = upd->vp[var];
	      for (b = 0; b < dim; b++)
		{
		  for (j = 0; j < ei->dof[var]; j++) 
		    {
		      phi_j = bf[var]->phi[j];	
		      lec->J[peqn][pvar][i][j] -= surfaceDiffusionCoeff * bfn->grad_phi[j][b] * (bfn->grad_phi[i][b] - tmp * fv->n[b])  * wt * det_J * h3;
		    }
		}
	    }

	  for (b = 0; b < dim; b++)
	    {
	      var = SHELL_NORMAL1 + b;
	      if (pd->v[var])
		{
		  pvar = upd->vp[var];
		  for (j = 0; j < ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];
		      lec->J[peqn][pvar][i][j] -= surfaceDiffusionCoeff * fv->grad_div_s_v[b] * (- tmp * phi_j)  * wt * det_J * h3;

		      for (a = 0; a < dim; a++)
			{
			  lec->J[peqn][pvar][i][j] -= surfaceDiffusionCoeff * fv->grad_div_s_v[a] * (- phi_j * bfn->grad_phi[i][b] * fv->n[a])  * wt * det_J * h3;
			} 
		    }
		}
	    }

	  tmp = 0.0;
	  for (b = 0; b < dim; b++)
	    {
	      tmp += fv->n[b] *  bfn->grad_phi[i][b];
	    } 
	  //   lec->R[peqn][i] -= surfaceDiffusionCoeff * (fv->grad_div_s_v[b] * (bfn->grad_phi[i][b] - tmp * fv->n[b]))  * wt * det_J * h3;
	  for (p = 0; p < dim; p++ )
	    {
	      var = MESH_DISPLACEMENT1 + p;
	      if (n_dof[var] > 0)   /* NOTE: Cannot use pd->v here! */
		{
		  pvar = upd->vp[var];   
		  /*
		   *  Here, the unknown is defined on the remote element. So get the DOF count from n_dof.
		   */
		  for (j = 0; j < n_dof[var]; j++) 
		    {
		      for (b = 0; b < dim; b++)
			{
			  lec->J[peqn][pvar][i][j] -= surfaceDiffusionCoeff * (fv->grad_div_s_v[b] * (bfn->grad_phi[i][b] - tmp * fv->n[b]))   * wt * (fv->dsurfdet_dx[p][j]);
			}
		    }
		  /*
		   *  This has been checked out in the debugger and it passes. What's going on?
		   *  The dgrad_phi_dmesh is the derivative of the gradients of the shell basis function for the shell variable with
		   *  respect to the mesh variables. Now the mesh variables are delineated by the bulk element. However, all derivatives
		   *  wrt unknowns of the mesh variables that are not on the shell must be zero (yes this is true). 
		   *  What's actually in dgrad_phi_dmesh is the variation of the gradient of the shell variables with respect to the
		   *  displacement unknowns for nodes which are on the shell (in the shell local node number convention). We then
		   *  use dof_map[] to change the shell local node number convention to the bulk element local node convention.
		   *  this then gets put in the correct place in load_lec().  
		   *  We use ei->dof[var] here because it correctly identifies the number of mesh unknowns that are on the shell (3)
		   *  even though mesh unknowns are only defined on the bulk elements.
		   */
		  for (j = 0; j < ei->dof[var]; j++)
		    {
		      jk = dof_map[j];
		      for (b = 0; b < dim; b++) 
			{
			  lec->J[peqn][pvar][i][jk] -= surfaceDiffusionCoeff * fv->grad_div_s_v[b] * (bfn->d_grad_phi_dmesh[i][b][p][j]) * wt *  det_J * h3 ;
			  lec->J[peqn][pvar][i][jk] -= surfaceDiffusionCoeff * fv->d_grad_div_s_v_dmesh[b][p][j] * (bfn->grad_phi[i][b] - tmp * fv->n[b]) * wt *  det_J * h3;
		
			  for (a = 0; a < dim; a++) 
			    {
			      lec->J[peqn][pvar][i][jk] -= surfaceDiffusionCoeff * fv->grad_div_s_v[b] * (- fv->n[a] * bfn->d_grad_phi_dmesh[i][a][p][j] * fv->n[b]) * wt *  det_J * h3;
			    }
			}
		    }
		}
	    }


	}
      //-----------------------------------------------------------------------------------------------------------------------------------------------------------------
      eqn = R_SHELL_SURF_CURV;
      peqn = upd->ep[eqn];
      bfn = bf[eqn];
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  phi_i = bf[eqn]->phi[i];
	  //  lec->R[peqn][i] += (fv->curv - curv) * phi_i * wt * det_J * h3;
	  var = eqn;
	  if (pd->v[var])
	    {
	      pvar = upd->vp[var];   
	      for (j = 0; j < ei->dof[var]; j++) 
		{
		  phi_j = bf[var]->phi[j];
		  lec->J[peqn][pvar][i][j] += phi_j * phi_i * wt * det_J * h3;
		}
	    }

	  for (b = 0; b < dim; b++ )
	    {
	      var = SHELL_NORMAL1 + b;
	      if (n_dof[var] > 0)  
		{
		  pvar = upd->vp[var];   
		  for (j = 0; j < ei->dof[var]; j++) 
		    {
		      lec->J[peqn][pvar][i][j] -=  (-0.5) * d_div_s_n_dn[b][j] * phi_i * wt * det_J * h3;
		    }
		}
	    }

	  for (b = 0; b < dim; b++ )
	    {
	      var = MESH_DISPLACEMENT1 + b;
	      if (n_dof[var] > 0)   /* NOTE: Cannot use pd->v here! */
		{
		  pvar = upd->vp[var];   

		  /* Here, the unknown is defined on the remote element.
		   * So get the DOF count from n_dof. */

		  for (j = 0; j < n_dof[var]; j++) 
		    {
		      //      if (fabs( d_div_s_n_dx[b][j]) != 0.0) {
		      //	printf("we are here\n");
		      //    }
		      lec->J[peqn][pvar][i][j] -=  (-0.5) * d_div_s_n_dx[b][j] * phi_i * wt * det_J * h3;
		    }

		  for (j = 0; j < n_dof[var]; j++)
		    {
		      lec->J[peqn][pvar][i][j] += (fv->curv - curv) * phi_i * wt * (fv->dsurfdet_dx[b][j]);
		    }
		}
	    }

	  /*
	   *  Add surface diffusion Jacobian terms
	   */
	  tmp = 0.0;
	  for (a = 0; a < dim; a++)
	    {
	      tmp += fv->n[a] * bfn->grad_phi[i][a];
	    } 

       	  //   lec->R[peqn][i] -= surfaceDiffusionCoeff * (fv->grad_curv[b] * (bfn->grad_phi[i][b] - tmp * fv->n[b]))  * wt * det_J * h3;
	
	  var = eqn;
	  if (pd->v[var])
	    {
	      pvar = upd->vp[var];
	      for (b = 0; b < dim; b++)
		{
		  for (j = 0; j < ei->dof[var]; j++) 
		    {
		      lec->J[peqn][pvar][i][j] -= surfaceDiffusionCoeff3 * bfn->grad_phi[j][b] * (bfn->grad_phi[i][b] - tmp * fv->n[b])  * wt * det_J * h3;
		    }
		}
	    }

	  for (b = 0; b < dim; b++)
	    {
	      var = SHELL_NORMAL1 + b;
	      if (pd->v[var])
		{
		  pvar = upd->vp[var];
		  for (j = 0; j < ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];
		      lec->J[peqn][pvar][i][j] -= surfaceDiffusionCoeff3 * fv->grad_curv[b] * (- tmp * phi_j)  * wt * det_J * h3;
		      for (a = 0; a < dim; a++)
			{
			  lec->J[peqn][pvar][i][j] -= surfaceDiffusionCoeff3 * fv->grad_curv[a] * (- phi_j * bfn->grad_phi[i][b] * fv->n[a]) * wt * det_J * h3;
			}
		    }
		}
	    }

	  tmp = 0.0;
	  for (b = 0; b < dim; b++)
	    {
	      tmp += fv->n[b] *  bfn->grad_phi[i][b];
	    } 
	  //   lec->R[peqn][i] -= surfaceDiffusionCoeff * (fv->grad_curv[b] * (bfn->grad_phi[i][b] - tmp * fv->n[b]))  * wt * det_J * h3;
	  for (p = 0; p < dim; p++ )
	    {
	      var = MESH_DISPLACEMENT1 + p;
	      if (n_dof[var] > 0)   /* NOTE: Cannot use pd->v here! */
		{
		  pvar = upd->vp[var];   
		  /*
		   *  Here, the unknown is defined on the remote element. So get the DOF count from n_dof.
		   */
		  for (j = 0; j < n_dof[var]; j++) 
		    {
		      for (b = 0; b < dim; b++)
			{
			  lec->J[peqn][pvar][i][j] -= surfaceDiffusionCoeff3 * (fv->grad_curv[b] * (bfn->grad_phi[i][b] - tmp * fv->n[b])) * wt * (fv->dsurfdet_dx[p][j]);
			}
		    }
		  /*
		   *  This has been checked out in the debugger and it passes. What's going on?
		   *  The dgrad_phi_dmesh is the derivative of the gradients of the shell basis function for the shell variable with
		   *  respect to the mesh variables. Now the mesh variables are delineated by the bulk element. However, all derivatives
		   *  wrt unknowns of the mesh variables that are not on the shell must be zero (yes this is true). 
		   *  What's actually in dgrad_phi_dmesh is the variation of the gradient of the shell variables with respect to the
		   *  displacement unknowns for nodes which are on the shell (in the shell local node number convention). We then
		   *  use dof_map[] to change the shell local node number convention to the bulk element local node convention.
		   *  this then gets put in the correct place in load_lec().  
		   *  We use ei->dof[var] here because it correctly identifies the number of mesh unknowns that are on the shell (3)
		   *  even though mesh unknowns are only defined on the bulk elements.
		   */
		  for (j = 0; j < ei->dof[var]; j++)
		    {
		      jk = dof_map[j];
		      for (b = 0; b < dim; b++) 
			{
			  lec->J[peqn][pvar][i][jk] -= surfaceDiffusionCoeff3 * fv->grad_curv[b] * (bfn->d_grad_phi_dmesh[i][b][p][j]) * wt *  det_J * h3 ;
			  lec->J[peqn][pvar][i][jk] -= surfaceDiffusionCoeff3 * fv->d_grad_curv_dmesh[b][p][j] * (bfn->grad_phi[i][b] - tmp * fv->n[b]) * wt *  det_J * h3;
			  for (a = 0; a < dim; a++) 
			    {
			      lec->J[peqn][pvar][i][jk] -= surfaceDiffusionCoeff3 * fv->grad_curv[b] * (- fv->n[a] * bfn->d_grad_phi_dmesh[i][a][p][j] * fv->n[b]) * wt * det_J * h3;
			    }
			}
		    }
		}
	    }
	}

      //-----------------------------------------------------------------------------------------------------------------------------------------------------------------

      // HKM This has all been checked out and it,s correct
      eqn = R_N_DOT_CURL_V;
      peqn = upd->ep[eqn];
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  phi_i = bf[eqn]->phi[i];
	  /* Diagonal contribution */
	  var = eqn;
	  if (pd->v[var])
	    {
	      pvar = upd->vp[var];  
	      for (j = 0; j < ei->dof[var]; j++) 
		{
		  phi_j = bf[var]->phi[j];  
		  lec->J[peqn][pvar][i][j] += phi_j * phi_i * wt * det_J * h3;
		}
	    }
	  
	  for (b = 0; b < dim; b++)
	    {
	      var = MESH_DISPLACEMENT1 + b;

	      if (n_dof[var] > 0)   /* NOTE: Cannot use pd->v here! */
		{
		  pvar = upd->vp[var];   

		  /* Here, the unknown is defined on the remote element.
		   * So get the DOF count from n_dof. */

		  for (j = 0; j < n_dof[var]; j++) 
		    {
		      lec->J[peqn][pvar][i][j] -= d_n_dot_curl_s_v_dx[b][j] * phi_i *
			wt * det_J * h3;
		    }
		  /*
		   * HKM there is an obvious error in the jacobian here. It's missing dependence
		   * on det_J and h3. 
		   *  Add in the dependence of det_J * h3.
		   *
		   *  det_J is the determinant of the shell equation basis functions wrt to the
		   *  mesh equation unknowns, which are defined on the parent element. This is 
		   *  pretty complicated as it involves two different sets of basis functions.
		   *  What this involves is an understanding of the mapping of the local shell
		   *  basis functions where this determinant is specified into the local node
		   *  numbers of the parent basis functions.
		   *
		   *  However, this dependence has already been calculated in the routine
		   *   surface_determinant_and_normal() and storred in fv->dsurfdet_dx[][]
		   */
		  for (j = 0; j < n_dof[var]; j++) 
		    {
		      lec->J[peqn][pvar][i][j] += (fv->n_dot_curl_s_v - n_dot_curl_s_v) * phi_i *
			wt * (fv->dsurfdet_dx[b][j]);
		    }

		}
	    }

	  for (b = 0; b < VIM; b++)
	    {
	      var = VELOCITY1 + b;

	      if (n_dof[var] > 0)   /* NOTE: Cannot use pd->v here! */
		{
		  pvar = upd->vp[var];

		  /* Here, the unknown is defined on the remote element.
		   * So get the DOF count from n_dof. */
		  /*
		   * HKM -> checked this jacobian contribution. It is now spot on
		   */
		  for (j = 0; j < n_dof[var]; j++) 
		    {
		      lec->J[peqn][pvar][i][j] -= d_n_dot_curl_s_v_dv[b][j] * phi_i * wt * det_J * h3;
		    }
		}
	    }
	}
    }
	
  return(err); 
} /* end of assemble_rheo_pieces */

  /*
    shell_diff_kinematic_bc(): Integrated kinematic BC on bulk mesh equations
    with diffusion on a shell surface.
  */
void
shell_diff_kinematic_bc(double func[DIM],
			double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
			const double x_dot[MAX_PDIM],
			const double tt,
			const double dt,
			const int id_side,
			const double wt,
			double xi[DIM],
			const Exo_DB *exo,
			const int sic_flag)
{
  int eqn = R_MESH1;
  int  j, p, q, var, pvar;
  double mass = 0.0;
  double diffusion = 0.0;
  int el1 = ei->ielem;
  int el2, nf;
  int *n_dof = NULL;
  int n_dofptr[MAX_VARIABLE_TYPES][MDE];
  //BASIS_FUNCTIONS_STRUCT *bfsh = bf[SHELL_DIFF_FLUX];
  PROBLEM_DESCRIPTION_STRUCT *pd0;
                                                                                
  /* Unpack variables from structures for local convenience. */
  // double h3    = fv->h3;
  //double det_J = bf[eqn]->detJ;
  double phi_t = (1.0 + 2.0 * tt) / dt;
  pd0 = pd;
                                                                                
  /* See if there is a friend for this element */
  nf = num_elem_friends[el1];
  if (nf == 0) return;
                                                                                
  /*
   * Get neighbor element number.
   * NOTE: For now, only one friend can be handled. If there is
   *       more than one, the first one will be used. This
   *       will be fixed at a later date.
   */
  el2 = elem_friends[el1][0];
  if (nf != 1) WH(-1, "WARNING: Not set up for more than one element friend!");
                                                                                
  /*
   * DOF counts will be needed to construct sensitivities of
   * local shell equations to neighbor bulk element variables.
   * Allocate this array to save them in.
   */
  n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));

  /*
   * This call will collect all field variable and basis function data
   * for the neighbor element "el2" and populate the fv structure,
   * then repopulate the local element structure. Thus, as long as
   * no variable is defined on both blocks, all data from both
   * elements will be available in the fv structure.
   * For example, pressure in the neighbor element is accessed with fv->P.
   * This is done to simplify the job of writing shell equations.
   */
  load_neighbor_var_data(el1, el2, n_dof, NULL, n_dofptr, id_side, xi, exo);
                                                                                
  /* Assemble the residual equation */
  /* Here, variables from the remote elements are used. */
  if (af->Assemble_Residual)
    {
                                                                                
      /* Mass term */
      if (pd0->e[eqn] & T_MASS)
	{
	  mass = 0.0;
	  for (p = 0; p < pd->Num_Dim; p++)
	    {
	      mass += fv->n[p] * fv_dot->x[p];
	    }
	  mass *= pd0->etm[eqn][(LOG2_MASS)];
	}

      /* Diffusion term */
      if (pd0->e[eqn] & T_DIFFUSION)
	{
	  diffusion = -fv->grad_sh_J[0] / fv->sdet; 
	  diffusion *= pd0->etm[eqn][(LOG2_DIFFUSION)]; 
	}
      func[0] += (mass + diffusion);
    }
                                                                                
  /* Assemble the Jacobian terms */
  /* Here, variables from the remote elements are used. */
  if (af->Assemble_Jacobian)
    {
                                                                                
      /* Mass term */

      /* J_kbc_n */
      /*
	mass = 0.0;
	for (p = 0; p < pd->Num_Dim; p++)
	{
	var = SHELL_NORMAL1 + p;
	pvar = upd->vp[var];
	if (pvar != -1)
	{
	for (j = 0; j < n_dof[var]; j++)
	{
	mass = bf[var]->phi[j] * fv_dot->x[p];
	mass *= pd->etm[eqn][(LOG2_MASS)];
	d_func[0][var][j] += mass;
	}
	}
	}
      */

      /* J_kbc_d */
      for (p = 0; p < pd->Num_Dim; p++)
	{
	  var = MESH_DISPLACEMENT1 + p;
	  pvar = upd->vp[var];
	  if (pvar != -1)
	    {
	      for (j = 0; j < ei->dof[var]; j++)
		{
		  mass = 0.0;
		  diffusion = 0.0;
		  for (q = 0; q < pd->Num_Dim; q++)
		    {
		      mass  = ( fv->dsnormal_dx[q][p][j] * fv_dot->x[q] );
		    }
		  mass += (fv->n[p] * phi_t * bf[var]->phi[j]);
		  mass *= pd->etm[eqn][(LOG2_MASS)];

		  diffusion = -fv->d_grad_sh_J_dmesh[0][p][j] / fv->sdet
		    + fv->grad_sh_J[0] * fv->dsurfdet_dx[p][j] / (fv->sdet * fv->sdet);
		  /*
		    diffusion = fv->grad_sh_J[0] * fv->dsurfdet_dx[p][j]
		    / (fv->sdet * fv->sdet);
		  */
		  diffusion *= pd0->etm[eqn][(LOG2_DIFFUSION)];

		  d_func[0][var][j] += (mass + diffusion);
		}
	    }
	}

      /* J_kbc_sh_J */
      diffusion = 0.0;
      var = SHELL_DIFF_FLUX;
      pvar = upd->vp[var];
      if (pvar != -1)
	{
	  for (j = 0; j < n_dof[var]; j++)
	    {
	      diffusion = -bf[var]->dphidxi[j][0] / fv->sdet;
	      diffusion *= pd0->etm[eqn][(LOG2_DIFFUSION)];

	      d_func[0][var][j] += diffusion;
	    }
	}
    }

  safe_free((void *) n_dof);
  return;
}  /* shell_diff_kinematic_bc */

/*
 * shell_surface_charge_bc(): Integrated BC on bulk potential equation.
 */
void
rep_force_shell_n_dot_f_bc(double func[DIM],
			   double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
			   const double x_dot[MAX_PDIM],
			   const double theta,
			   const double delta_t,
			   const int ip,
			   const int ip_total,
			   const int id_side,
			   const double wt,
			   double xi[DIM],
			   const Exo_DB *exo,
			   const int sic_flag)
{
  int  j, var, pvar, jvar, dim, j_id;
  double phi_j;
  int el1 = ei->ielem;
  int el2, nf, err;
  int *n_dof = NULL;
  int n_dofptr[MAX_VARIABLE_TYPES][MDE], dof_map[MDE];
  double *n_esp;

  int a, b, c;
  double shell_p, shell_pgrad, shear_stress;

  double TT[MAX_PDIM][MAX_PDIM];   /**  solid stresses  **/
  double dTT_drs[DIM][DIM][DIM][MDE];
  double dTT_dx[MAX_PDIM][MAX_PDIM][MAX_PDIM][MDE];
  double dTT_dp[MAX_PDIM][MAX_PDIM][MDE];
  double dTT_dc[MAX_PDIM][MAX_PDIM][MAX_CONC][MDE];
  double dTT_dp_liq[DIM][DIM][MDE];/* Sensitivity of stress tensor...
				      to nodal porous liquid pressure*/
  double dTT_dp_gas[DIM][DIM][MDE];/* Sensitivity of stress tensor...
				      to nodal porous gas pressure*/
  double dTT_dporosity[DIM][DIM][MDE];/* Sensitivity of stress tensor...
					 to nodal porosity*/
  double dTT_dsink_mass[DIM][DIM][MDE];/* Sensitivity of stress tensor...
					  to sink_mass*/
  double dTT_dT[MAX_PDIM][MAX_PDIM][MDE];
  double dTT_dmax_strain[MAX_PDIM][MAX_PDIM][MDE];
  double dTT_dcur_strain[MAX_PDIM][MAX_PDIM][MDE];
  double elast_modulus;

  double vconv[MAX_PDIM], vconv_old[MAX_PDIM];/*Calculated convection velocity*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;


  /*Vweb = mp->u_shell_user_par[3];
    roll_rad = mp->u_shell_user_par[4];
    x0 = mp->u_shell_user_par[5];
    gap_nom = mp->u_shell_user_par[6];
    gap = gap_nom+roll_rad - sqrt(SQUARE(roll_rad)-SQUARE(fv->x[0]-x0)) - fv->x[1];
    dgap_dx[0] = (fv->x[0]-x0)/sqrt(SQUARE(roll_rad)-SQUARE(fv->x[0]-x0));
    dgap_dx[1] = -1.;
    Vwebx = Vweb*sqrt(SQUARE(roll_rad)-SQUARE(fv->x[0]-x0))/roll_rad;
    dVwebx_dx[0] = (Vweb/roll_rad)*(x0-fv->x[0])/sqrt(SQUARE(roll_rad)-SQUARE(fv->x[0]-x0)
    );
    dVwebx_dx[1] = 0.0; 
  */


  if(af->Assemble_LSA_Mass_Matrix)
    return;

 
  /* See if there is a friend for this element */
  nf = num_elem_friends[el1];
  if (nf == 0) return;

  /*
   * Get neighbor element number.
   */
  el2 = elem_friends[el1][0];
  if (nf != 1) WH(-1, "WARNING: Not set up for more than one element friend!");

  /*
   * DOF counts will be needed to construct sensitivities of
   * local shell equations to neighbor bulk element variables.
   * Allocate this array to save them in.
   */
  n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));

  /*
   * This call will collect all field variable and basis function data
   * for the neighbor element "el2" and populate the fv structure,
   * then repopulate the local element structure. Thus, as long as
   * no variable is defined on both blocks, all data from both
   * elements will be available in the fv structure.
   * For example, pressure in the neighbor element is accessed with fv->P.
   * This is done to simplify the job of writing shell equations.
   */
  err = load_neighbor_var_data(el1, el2, n_dof, dof_map, n_dofptr, id_side, xi, exo);

  /* Assemble the residual equation */
  /*  initialize variables */ 
  dim = pd->Num_Dim;

  /* initialize some arrays */
  memset( TT, 0, sizeof(double)*DIM*DIM);
  if (af->Assemble_Jacobian) {
    memset( dTT_dx,         0, sizeof(double)*DIM*DIM*DIM*MDE);
    memset( dTT_dp,         0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_drs,        0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dc,         0, sizeof(double)*DIM*DIM*MAX_CONC*MDE);
    memset( dTT_dp_liq,     0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dp_gas,     0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dporosity,  0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dsink_mass, 0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dT,         0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dmax_strain,0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dcur_strain,0, sizeof(double)*DIM*DIM*MDE);
  }


  if(pd->e[R_MESH1] && cr->MeshMotion != ARBITRARY)
    {
      err = belly_flop(elc->lame_mu);
      EH(err, "error in belly flop");

      err = mesh_stress_tensor(TT, dTT_dx, dTT_dp, dTT_dc, dTT_dp_liq,
			       dTT_dp_gas, dTT_dporosity, dTT_dsink_mass, dTT_dT, dTT_dmax_strain, dTT_dcur_strain,
			       elc->lame_mu, elc->lame_lambda, delta_t, ei->ielem,
			       ip, ip_total);
      if ( cr->MeshMotion == LAGRANGIAN ||
	   cr->MeshMotion == DYNAMIC_LAGRANGIAN)
	{
	  err = get_convection_velocity(vconv, vconv_old, d_vconv, delta_t, theta);
	}
      else /* No inertia in an Arbitrary Mesh */
	{
	  memset( vconv, 0, sizeof(double)*MAX_PDIM );
	  if (pd->v[MESH_DISPLACEMENT1])
	    memset(d_vconv->X, 0, DIM*DIM*MDE*sizeof(dbl));
	  if (pd->v[VELOCITY1] || pd->v[POR_LIQ_PRES])
	    memset(d_vconv->v, 0, DIM*DIM*MDE*sizeof(dbl));
	  if (pd->v[MASS_FRACTION] || pd->v[POR_LIQ_PRES])
	    memset(d_vconv->C, 0, DIM*MAX_CONC*MDE*sizeof(dbl));
	  if (pd->v[TEMPERATURE])
	    memset(d_vconv->T, 0, DIM*MDE*sizeof(dbl) );
	}
      /* For LINEAR ELASTICITY */
      if (cr->MeshFluxModel == LINEAR)
	{
          if (dim == 2)
	    {
	      TT[2][2] = 1.;
	      TT[1][2] = 0.;
	      TT[0][2] = 0.;
	    }
	}

      /*  For Hookian Elasticity and shrinkage */
      else
	{
          if (dim == 2)
	    {
	      elast_modulus = elc->lame_mu;
	      if (cr->MeshMotion == ARBITRARY)
		{
		  TT[2][2] = (1. - fv->volume_change) * elast_modulus;
		}
	      else
		{
		  if (cr->MeshFluxModel == NONLINEAR ||
		      cr->MeshFluxModel == HOOKEAN_PSTRAIN ||
		      cr->MeshFluxModel == INCOMP_PSTRAIN )
		    TT[2][2] = (1. - pow(fv->volume_change,2./3.))
		      * elast_modulus - fv->P;
		  /*              if (cr->MeshFluxModel == INCOMP_PSTRESS) */
		  else  TT[2][2] = 0.;
		}
	      TT[1][2] = 0.;
	      TT[0][2] = 0.;
	    }
	}
    } /* end of STRESS_TENSOR */

  /* calculate real-solid stress here !!*/
  if(pd->e[R_SOLID1] && cr->MeshMotion != ARBITRARY)
    {
      err = belly_flop_rs(elc_rs->lame_mu);
      EH(err, "error in belly flop");

      err = solid_stress_tensor(TT, dTT_dx, dTT_drs, dTT_dp,
				dTT_dc,  dTT_dp_liq, dTT_dp_gas, dTT_dporosity,
				dTT_dT, dTT_dmax_strain, elc_rs->lame_mu, elc_rs->lame_lambda);
      if (  cr->MeshMotion == LAGRANGIAN ||
	    cr->MeshMotion == DYNAMIC_LAGRANGIAN ||
	    cr->MeshMotion == TOTAL_ALE)
	{
	  err = get_convection_velocity_rs(vconv, vconv_old, d_vconv, delta_t, theta);
	  if ( pd->TimeIntegration != STEADY &&
	       pd->etm[R_SOLID1][(LOG2_ADVECTION)] )
	    {
	      for ( a=0; a<dim; a++)
		{
		  vconv[a] = fv_dot->d_rs[a];
		  for(b=0; b < VIM; b++)
		    {
		      vconv[a] -= fv_dot->d_rs[b]*fv->grad_d_rs[b][a];
		    }
		}
	      for (a = 0; a < VIM; a++)
		{
		  for (b = 0; b < VIM; b++)
		    {
		      var = SOLID_DISPLACEMENT1 + b;
		      for (j = 0; j <ei->dof[var]; j++)
			{
			  phi_j = bf[var]->phi[j];
			  d_vconv->rs[a][b][j] += phi_j * (1.+2.*theta)*delta_t*delta(a,b) ;
			  d_vconv->rs[a][b][j] -= (1.+2.*theta)*phi_j*delta_t*fv->grad_d_rs[b][a] ;
			  for(c=0; c < VIM; c++)

			    {
			      /*grad_phi_e[i][j][k][l] = e_k e_l : grad(phi_i e_j )*/
			      d_vconv->rs[a][b][j] -= fv_dot->d_rs[c]* bf[var]->grad_phi_e[j][b][c][a];
			    }
			}
		    }
		}
	    }
	}
      else /* No inertia in an Arbitrary Mesh */
	{
	  for ( a=0; a<dim; a++)
	    {
	      vconv[a] = 0.;
	    }
	}

      if (dim == 2)
	{
	  elast_modulus = elc_rs->lame_mu;
	  if (cr->RealSolidFluxModel == NONLINEAR ||
	      cr->RealSolidFluxModel == HOOKEAN_PSTRAIN ||
	      cr->RealSolidFluxModel == INCOMP_PSTRAIN )
	    TT[2][2] = (1.-pow(fv->volume_change,2./3.))*elast_modulus-fv->P;
	  /*              if (cr->MeshFluxModel == INCOMP_PSTRESS) */
	  else  TT[2][2] = 0.;
	  TT[1][2] = 0.;
	  TT[0][2] = 0.;
	}
    } /* end of REAL_STRESS_TENSOR */


      /* Here, variables from the remote elements are used. */ 
 
  shell_p = 0;
  shell_pgrad = 0;
  for (j=0 ; j<n_dof[SHELL_LUBP] ; j++)
    {
      n_esp = x_static + n_dofptr[SHELL_LUBP][j];
      shell_p += *n_esp * bf[SHELL_LUBP]->phi[j];  
      shell_pgrad += *n_esp * bf[SHELL_LUBP]->dphidxi[j][0];  
    }
  shell_pgrad /= (-fv->stangent[0][0]*fv->sdet);
  shear_stress = 0.0;
#if 0
  fprintf(stderr,"rep_force %g %g %g %g %g\n",fv->x[0],shell_pgrad,shell_p,vconv[0],vconv[1]);
#endif
  if ( af->Assemble_Residual )
    {
      for ( a=0 ; a<pd->Num_Dim ; a++)
	{
	  func[a] += -shell_p * fv->snormal[a];
	  func[a] += -shear_stress*fv->stangent[0][a];
	  if (ei->ielem_dim == 3)
	    {
	      func[a] += 0.0*fv->stangent[1][a];
	    }
	}
      if( sic_flag)
	{
	  for ( a=0 ; a<pd->Num_Dim ; a++)    
	    {
	      for ( b=0; b<dim; b++)
		{
		  func[a] += fv->snormal[a]*TT[a][b]*fv->snormal[b];
		  func[a] += fv->snormal[a]*TT[a][b]*fv->stangent[0][b];
		  if (ei->ielem_dim == 3)
		    {
		      func[a] += fv->snormal[a]*TT[a][b]*fv->stangent[1][b];
		    }
		}
	    }
	}
    }

  /* Assemble Jacobian sensitivities */
  if ( af->Assemble_Jacobian )
    {

      /* Surface charge: REMOTE sensitivity */
      var = SHELL_LUBP;
      pvar = upd->vp[var];
      if (pvar != -1)
	{
	  for (j=0; j<n_dof[var]; j++)
	    {
	      /*  Find the right variable in the bulk context */
	      for(el1=0 ; el1<n_dof[var] ; el1++)   {
		if(dof_map[j] == ei->dof_list[var][el1])el2=el1;
	      }
	      for (a=0; a<ei->ielem_dim; a++)
		{
                  d_func[a][var][el2] -= bf[var]->phi[j] * fv->snormal[a];
                  d_func[a][var][el2] += 0.0 * fv->stangent[0][a];
		}
	    }
	}
      /* 
       *  Evaluate sensitivity to displacements d()/dx 
       */
      for (jvar=0; jvar<dim; jvar++)
        {
          var = MESH_DISPLACEMENT1 + jvar;
          if (pd->v[var])
            {
              for ( j=0; j<ei->dof[var]; j++)
                {
		  phi_j = bf[var]->phi[j];
                  for (a=0; a<ei->ielem_dim; a++)
                    {

                      d_func[a][var][j] += -shell_p*fv->dsnormal_dx[a][jvar][j];
		      /*		      dshear_dx[jvar] = phi_j*dgap_dx[jvar]*(mp->viscosity*
					      (Vwebx-vconv[0])*(-1./SQUARE(gap)) - shell_pgrad/2.);
					      dshear_dx[jvar] += mp->viscosity*(dVwebx_dx[jvar]*phi_j - 
					      d_vconv->X[0][jvar][j])/gap;
					      dshear_dx[jvar] += mp->viscosity*gap/2.*
					      (shell_pgrad*
					      (-fv->stangent[0][0]*fv->dsurfdet_dx[jvar][j]
					      -fv->sdet*fv->dstangent_dx[0][0][jvar][j])
					      /(-fv->stangent[0][0]*fv->sdet));
					      d_func[a][var][j] += 
					      -shear_stress*fv->dstangent_dx[0][a][jvar][j] 
					      - fv->stangent[0][a]*dshear_dx[jvar];
		      */
		      if (ei->ielem_dim == 3)
        		{
			  d_func[a][var][j] += 0.0*fv->dstangent_dx[1][a][jvar][j];
        		}

                    }
                }
            }
        }

      if( sic_flag)
	{
	  /*   	     var = TEMPERATURE;
		     if (pd->v[var] )
		     {
		     for (j=0; j<ei->dof[var]; j++)
         	     {
		     for ( p=0; p<pd->Num_Dim; p++)
		     { grad_phi_j = bf[var]->grad_phi[j][p]; 
		     d_func[0][var][j] += fv->snormal[p]*
		     ( d_p_dT[j] * efield[p]);
		     }
		     }
		     }
	  */
 
	  for(jvar=0; jvar<dim; jvar++)
	    {
	      var = MESH_DISPLACEMENT1 + jvar;
	      if (pd->v[var])
		{
		  for (j_id = 0; j_id < ei->dof[var]; j_id++)
		    {
		      for(a=0; a<dim; a++)
                        {
			  for(b=0; b<dim; b++)
                            {
			      d_func[a][var][j_id] +=
				fv->dsnormal_dx[a][jvar][j_id]*TT[a][b]*fv->snormal[b] +
				fv->snormal[a]*dTT_dx[a][b][jvar][j_id]*fv->snormal[b] +
				fv->snormal[a]*TT[a][b]*fv->dsnormal_dx[b][jvar][j_id];
			      d_func[a][var][j_id] +=
				fv->dsnormal_dx[a][jvar][j_id]*TT[a][b]*fv->stangent[0][b] +
				fv->snormal[a]*dTT_dx[a][b][jvar][j_id]*fv->stangent[0][b] +
				fv->snormal[a]*TT[a][b]*fv->dstangent_dx[0][b][jvar][j_id];
			      if (ei->ielem_dim == 3)
				{
				  d_func[a][var][j_id] +=
				    fv->dsnormal_dx[a][jvar][j_id]*TT[a][b]*fv->stangent[1][b] +
				    fv->snormal[a]*dTT_dx[a][b][jvar][j_id]*fv->stangent[1][b] +
				    fv->snormal[a]*TT[a][b]*fv->dstangent_dx[1][b][jvar][j_id];
				}
                            }
                        }
		    }
		}
	    }
	}  /*sic_flag */
    } /* End of if(af->Assemble_Jacobian) */

  /* Done */
  safe_free ((void *) n_dof);
  return;
}  /* rep_force_shell_n_dot_f */
void
surface_user_shell_bc(double R[MAX_PROB_VAR+MAX_CONC][MAX_NODES_PER_SIDE],
		      double J[MAX_PROB_VAR+MAX_CONC][MAX_PROB_VAR+MAX_CONC][MAX_NODES_PER_SIDE][MDE],
		      const int bm,
		      const int *n_dof,
		      const double wt,
		      const double theta,
		      const double delta_t,
		      const double *coord)
{
  int i, j, q, dofs, err;
  int eqn, peqn, var, pvar;
  double boundary, em, det_J, phi_i, phi_j;
  double vconv[MAX_PDIM], vconv_old[MAX_PDIM];/*Calculated convection velocity*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;

  double gap, pgrad_lub, dpgrad_lub_dx[DIM];
  double Vweb, roll_rad, x0, gap_nom;
  double Vwebx, dgap_dx[DIM], dVwebx_dx[DIM];
  double ups_xloc, dns_xloc, ups_width, dns_width, heavi, F_dns, F_ups;

  Vweb = mp->u_shell_user_par[3];
  roll_rad = mp->u_shell_user_par[4];
  x0 = mp->u_shell_user_par[5];
  gap_nom = mp->u_shell_user_par[6];
  ups_xloc = mp->u_shell_user_par[7];
  dns_xloc = mp->u_shell_user_par[8];
  ups_width = mp->u_shell_user_par[9];
  dns_width = mp->u_shell_user_par[10];
  gap = gap_nom+roll_rad - sqrt(SQUARE(roll_rad)-SQUARE(coord[0]-x0)) - coord[1];
  dgap_dx[0] = (coord[0]-x0)/sqrt(SQUARE(roll_rad)-SQUARE(coord[0]-x0));
  dgap_dx[1] = -1.;
  Vwebx = Vweb*sqrt(SQUARE(roll_rad)-SQUARE(coord[0]-x0))/roll_rad;
  dVwebx_dx[0] = (Vweb/roll_rad)*(x0-coord[0])/sqrt(SQUARE(roll_rad)-SQUARE(coord[0]-x0)
						    );
  dVwebx_dx[1] = 0.0; 
  F_dns = dns_xloc - coord[0];
  if (fabs(F_dns) > 0.5*dns_width)
    { heavi = ( F_dns < 0 ) ? 0.0 :1.0;}
  else
    {heavi = 0.5*(1.+2.*F_dns/dns_width+sin(M_PIE*2.*F_dns/dns_width)/M_PIE);}
  F_ups = coord[0] - ups_xloc;
  if (F_ups < -0.5*ups_width)
    { heavi = 0.0;}
  else if (F_ups < 0.5*ups_width)
    {heavi = 0.5*(1.+2.*F_ups/ups_width+sin(M_PIE*2.*F_ups/ups_width)/M_PIE);}

  /* Check for active SURFACE_CHARGE equation */
  eqn = R_SHELL_USER;

  /* Only BOUNDARY terms are assembled here */
  if ( !(pd->e[eqn] & T_BOUNDARY) ) return;

  /* Unpack variables from structures for local convenience. */
  det_J = fv->sdet;
  peqn  = upd->ep[eqn];
  em    = pd->etm[eqn][(LOG2_BOUNDARY)];

  if (pd_glob[bm]->e[R_MESH1] && cr->MeshMotion != ARBITRARY)
    {
      err = belly_flop(elc->lame_mu);
      EH(err, "error in belly flop");

      if ( cr->MeshMotion == LAGRANGIAN ||
	   cr->MeshMotion == DYNAMIC_LAGRANGIAN)
	{
	  err = get_convection_velocity(vconv, vconv_old, d_vconv, delta_t, theta);
	}
      else /* No inertia in an Arbitrary Mesh */
	{
	  memset( vconv, 0, sizeof(double)*MAX_PDIM );
	  if (pd_glob[bm]->v[MESH_DISPLACEMENT1])
	    memset(d_vconv->X, 0, DIM*DIM*MDE*sizeof(dbl));
	  if (pd_glob[bm]->v[VELOCITY1] || pd_glob[bm]->v[POR_LIQ_PRES])
	    memset(d_vconv->v, 0, DIM*DIM*MDE*sizeof(dbl));
	  if (pd_glob[bm]->v[MASS_FRACTION] || pd_glob[bm]->v[POR_LIQ_PRES])
	    memset(d_vconv->C, 0, DIM*MAX_CONC*MDE*sizeof(dbl));
	  if (pd_glob[bm]->v[TEMPERATURE])
	    memset(d_vconv->T, 0, DIM*MDE*sizeof(dbl) );
	}
    } /* end of STRESS_TENSOR */


  pgrad_lub = 6.*mp->viscosity/SQUARE(gap)*(vconv[0]-Vwebx);  
#if 0
  fprintf(stderr,"grad %g %g %g\n",pgrad_lub,vconv[0],Vwebx);
  fprintf(stderr,"F %g %g %g %g\n",coord[0],F_dns,F_ups, heavi);
#endif
  /* Residuals */
  if (af->Assemble_Residual)
    {
      for (i = 0; i < ei->dof[eqn]; i++)
        {
          phi_i = bf[eqn]->dphidxi[i][0];  

          /* Boundary (jump) term */
	  boundary = heavi*pgrad_lub*fv->stangent[0][0]*det_J;
          boundary *= (phi_i * wt * em);
          R[peqn][i] += boundary;
        }
    }

  /* Jacobian */
  if (af->Assemble_Jacobian)
    {
      for (i = 0; i < ei->dof[eqn]; i++)
        {
          phi_i = bf[eqn]->dphidxi[i][0];  

          /* J_qs_x:  Remote mesh sensitivity */
          /* NOTE: It may be necessary to restrict this to one
           *       side of a shell if there are two sides! */
          for (q = 0; q < pd->Num_Dim; q++)
            {
              var = MESH_DISPLACEMENT1 + q;
              if (pd_glob[bm]->v[var])
                {
                  pvar = upd->vp[var];
                  dofs = n_dof[var];

                  for (j = 0; j < dofs; j++)
                    {
		      phi_j = bf[eqn]->phi[j];
		      dpgrad_lub_dx[q] = -12.*mp->viscosity*dgap_dx[q]
			*(vconv[0]-Vwebx)/pow(gap,3)*phi_j;
		      dpgrad_lub_dx[q] += 6.*mp->viscosity/SQUARE(gap)
			*(d_vconv->X[0][q][j] - dVwebx_dx[q]*phi_j);
		      boundary = det_J*fv->stangent[0][0]*dpgrad_lub_dx[q]
			+ pgrad_lub*(fv->stangent[0][0]*
				     fv->dsurfdet_dx[q][j] + 
				     det_J*fv->dstangent_dx[0][0][q][j]);
		      boundary *= heavi;
                      boundary *= (phi_i * wt * em);
                      J[peqn][pvar][i][j] += boundary;
                    }
                }
            }
        }
    }

  return;
}  /* End of surface_user_shell_bc() */

void
surface_lubrication_shell_bc(double R[MAX_PROB_VAR+MAX_CONC][MAX_NODES_PER_SIDE],
			     double J[MAX_PROB_VAR+MAX_CONC][MAX_PROB_VAR+MAX_CONC][MAX_NODES_PER_SIDE][MDE],
			     const int bm,
			     const int *n_dof,
			     const double wt,
			     const double theta,
			     const double delta_t,
			     const double *coord,
			     int dof_map[MDE],
			     int n_dofptr[MAX_VARIABLE_TYPES][MDE])
{
  int i, j, q, dofs, err;
  int eqn, peqn, var, pvar;
  double boundary, em, det_J, phi_i, phi_j;
  double vconv[MAX_PDIM], vconv_old[MAX_PDIM];/*Calculated convection velocity*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;

  double *n_esp;
  /*  variables for lubrication approximation	*/
  double gap, shell_p, shell_pgrad;
  int model_id;
  double flow_inlet, Vweb, roll_rad, x0, gap_nom, flow_target;
  double Vwebx;
  double ups_xloc, dns_xloc, ups_width, dns_width, heavi, F_dns, F_ups;

  double qlub[DIM],dq_gradP[DIM][DIM],dq_dX[DIM][DIM],grad_P[DIM];
  double dq_dVb[DIM][DIM];
#ifdef SECOR_HEAT_FLUX
  double dq_dVt[DIM][DIM], Vt[DIM], Vb[DIM], dgap_dx[DIM];
#endif
  model_id = (int) mp->u_shell_user_par[0];
  if( model_id == 1)
    {
      flow_inlet = mp->u_shell_user_par[3];
      Vweb = mp->u_shell_user_par[4];
      roll_rad = mp->u_shell_user_par[5];
      x0 = mp->u_shell_user_par[6];
      gap_nom = mp->u_shell_user_par[7];
      ups_xloc = mp->u_shell_user_par[8];
      dns_xloc = mp->u_shell_user_par[9];
      ups_width = mp->u_shell_user_par[10];
      dns_width = mp->u_shell_user_par[11];

      gap = gap_nom+roll_rad - sqrt(SQUARE(roll_rad)-SQUARE(coord[0]-x0)) - coord[1];
#ifdef SECOR_HEAT_FLUX
      dgap_dx[0] = (coord[0]-x0)/sqrt(SQUARE(roll_rad)-SQUARE(coord[0]-x0));
      dgap_dx[1] = -1.;
#endif
      Vwebx = Vweb*sqrt(SQUARE(roll_rad)-SQUARE(coord[0]-x0))/roll_rad;
      F_dns = dns_xloc - coord[0];
      if (fabs(F_dns) > 0.5*dns_width)
	{ heavi = ( F_dns < 0 ) ? 0.0 :1.0;}
      else
	{heavi = 0.5*(1.+2.*F_dns/dns_width+sin(M_PIE*2.*F_dns/dns_width)/M_PIE);}
      F_ups = coord[0] - ups_xloc;
      if (F_ups < -0.5*ups_width)
	{ heavi = 0.0;}
      else if (F_ups < 0.5*ups_width)
	{heavi = 0.5*(1.+2.*F_ups/ups_width+sin(M_PIE*2.*F_ups/ups_width)/M_PIE);}
    }
  else
    { EH(-1,"invalid lubrication model_id\n");}

  /* Check for active SURFACE_CHARGE equation */
  eqn = R_SHELL_LUBP;

  /* Only BOUNDARY terms are assembled here */
  if ( !(pd->e[eqn] & T_ADVECTION) ) return;

  /* Unpack variables from structures for local convenience. */
  det_J = fv->sdet;
  peqn  = upd->ep[eqn];
  em    = pd->etm[eqn][(LOG2_ADVECTION)];

  if(pd_glob[bm]->e[R_MESH1] && cr->MeshMotion != ARBITRARY)
    {
      err = belly_flop(elc->lame_mu);
      EH(err, "error in belly flop");

      if ( cr->MeshMotion == LAGRANGIAN ||
	   cr->MeshMotion == DYNAMIC_LAGRANGIAN)
	{
	  err = get_convection_velocity(vconv, vconv_old, d_vconv, delta_t, theta);
	}
      else /* No inertia in an Arbitrary Mesh */
	{
	  memset( vconv, 0, sizeof(double)*MAX_PDIM );
	  if (pd_glob[bm]->v[MESH_DISPLACEMENT1])
	    memset(d_vconv->X, 0, DIM*DIM*MDE*sizeof(dbl));
	  if (pd_glob[bm]->v[VELOCITY1] || pd_glob[bm]->v[POR_LIQ_PRES])
	    memset(d_vconv->v, 0, DIM*DIM*MDE*sizeof(dbl));
	  if (pd_glob[bm]->v[MASS_FRACTION] || pd_glob[bm]->v[POR_LIQ_PRES])
	    memset(d_vconv->C, 0, DIM*MAX_CONC*MDE*sizeof(dbl));
	  if (pd_glob[bm]->v[TEMPERATURE])
	    memset(d_vconv->T, 0, DIM*MDE*sizeof(dbl) );
	}
    } /* end of STRESS_TENSOR */

      /* Here, variables from the remote elements are used. */ 
  /* See if there is a friend for this element */
 
  shell_p = 0;
  shell_pgrad = 0;
  for (j=0 ; j<n_dof[SHELL_LUBP] ; j++)
    {
      n_esp = x_static + n_dofptr[SHELL_LUBP][j];
      shell_p += *n_esp * bf[SHELL_LUBP]->phi[j];  
      shell_pgrad += *n_esp * bf[SHELL_LUBP]->dphidxi[j][0];  
    }
  shell_pgrad /= (fv->stangent[0][0]*fv->sdet);

  grad_P[0] = shell_pgrad;
  grad_P[1] = 0.;


#if defined SECOR_HEAT_FLUX
  Vt[0] = Vwebx;  Vt[1] = 0.;
  Vb[0] = vconv[0]; Vb[1] = 0.;
  usr_heat_flux(grad_P, qlub, dq_gradP, dq_dX, 0.0, 
			gap, dgap_dx, Vb, Vt, dq_dVb, dq_dVt);
#else
  usr_heat_flux(grad_P, qlub, dq_gradP, dq_dX, 0.0);
  printf("untested\n");
  exit(-1);
#endif
#if 0
  fprintf(stderr,"lub_shell %g %g %g %g %g\n",coord[0],shell_pgrad,shell_p,vconv[0],vconv[1]);
  fprintf(stderr,"grad %g %g %g \n",shell_pgrad, fv->stangent[0][0],fv->sdet);
#endif

  flow_target = heavi*flow_inlet + (1.-heavi)*gap*(Vwebx+vconv[0])*0.5;
#if 0
  fprintf(stderr,"grad %g %g %g\n",pgrad_lub,vconv[0],Vwebx);
  fprintf(stderr,"F %g %g %g %g\n",coord[0],F_dns,F_ups, heavi);
#endif
  /* Residuals */
  if (af->Assemble_Residual)
    {
      for (i = 0; i < ei->dof[eqn]; i++)
        {
	  /*          phi_i = bf[eqn]->phi[i];    */
          phi_i = bf[eqn]->dphidxi[i][0];  

          /* Boundary (jump) term */
          boundary = (qlub[0]-flow_target)*det_J;
          boundary *= (phi_i * wt * em);
          R[peqn][i] += boundary;
	  /*fprintf(stderr,"resid %d %g %g \n",i, boundary, R[peqn][i]); */
        }
    }

  /* Jacobian */
  if (af->Assemble_Jacobian)
    {
      for (i = 0; i < ei->dof[eqn]; i++)
        {
	  /*          phi_i = bf[eqn]->phi[i];    */
          phi_i = bf[eqn]->dphidxi[i][0];  

          /* Surface charge: REMOTE sensitivity */
          var = SHELL_LUBP;
          pvar = upd->vp[var];
          if (pvar != -1)
            {
              for (j=0; j<n_dof[var]; j++)
                {
		  /*  Find the right variable in the bulk context */
		  phi_j = bf[var]->dphidxi[j][0];  
		  boundary = dq_gradP[0][0]*phi_j/(fv->stangent[0][0]);
		  boundary *= (phi_i * wt * em);
		  J[peqn][pvar][i][j] += boundary;
                }
            }
          /* J_qs_x:  Remote mesh sensitivity */
          /* NOTE: It may be necessary to restrict this to one
           *       side of a shell if there are two sides! */
          for (q = 0; q < pd->Num_Dim; q++)
            {
              var = MESH_DISPLACEMENT1 + q;
              if (pd_glob[bm]->v[var])
                {
                  pvar = upd->vp[var];
                  dofs = n_dof[var];

                  for (j = 0; j < dofs; j++)
                    {
		      phi_j = bf[eqn]->phi[j];
		      boundary = det_J*(dq_dX[0][q] + dq_dVb[0][0]*d_vconv->X[0][q][j]
					+ dq_dVb[0][1]*d_vconv->X[1][q][j])*phi_j;
		      boundary += (qlub[0]-flow_target)*fv->dsurfdet_dx[q][j];
		      boundary *= (phi_i * wt * em);
                      J[peqn][pvar][i][j] += boundary;
                    }
                }
            }
        }
    }

  return;
}  /* End of surface_lubrication_shell_bc() */


/*****************************************************************************/
/***assemble_lubrication******************************************************/
/*  _______________________________________________________________________  */

/* assemble_lubrication -- assemble terms (Residual & Jacobian) for scalar Reynolds
 *                         Lubrication equation
 *
 * in:
 * 	ei -- pointer to Element Indices	structure
 *	pd -- pointer to Problem Description	structure
 *	af -- pointer to Action Flag		structure
 *	bf -- pointer to Basis Function		structure
 *	fv -- pointer to Field Variable		structure
 *  fv_old -- pointer to old Diet Field Variable	structure
 *  fv_dot -- pointer to dot Diet Field Variable	structure
 *	cr -- pointer to Constitutive Relation	structure
 *	md -- pointer to Mesh Derivative	structure
 *	me -- pointer to Material Entity	structure
 *
 * out:
 *	a   -- gets loaded up with proper contribution 
 *	lec -- gets loaded up with local contributions to resid, Jacobian
 * 	r   -- residual RHS vector
 *
 * Created:	Friday March 13 07:48:01 MST 2009 prschun@sandia.gov
 *
 */
/*ARGSUSED*/
int 
assemble_lubrication(const int EQN,     /* equation type: either R_LUBP or R_LUBP2 */
		     double time,	/* present time value */ 
		     double tt,   	/* parameter to vary time integration from 
			   	         * explicit (tt = 1) to implicit (tt = 0)    */
		     double dt,         /* current time step size */
		     double xi[DIM],    /* Local stu coordinates */
		     const Exo_DB *exo)  
{
  int eqn, var, peqn, pvar, p, a, b, k, jk;
  int i = -1, j, status; //, err;
  int *n_dof = NULL;
  int dof_map[MDE];
  
  // dbl toggle_dh_dependence = 0.;

  dbl H, dH_dtime; 
  dbl H_U, dH_U_dtime, H_L, dH_L_dtime;
  dbl dH_U_dX[DIM],dH_L_dX[DIM], dH_dtime_dmesh[DIM][MDE];
  dbl dH_dtime_drealsolid[DIM][MDE];
  dbl dH_U_dp, dH_U_ddh;
  dbl veloU[DIM], veloL[DIM];
  dbl diffusion,  source;

  /*
   * Basis functions and derivatives
   */
  dbl phi_i, grad_phi_i[DIM], grad_II_phi_i[DIM], d_grad_II_phi_i_dmesh[DIM][DIM][MDE];
  dbl phi_j, grad_phi_j[DIM], grad_II_phi_j[DIM], d_grad_II_phi_j_dmesh[DIM][DIM][MDE];

  /*
   * Bail out fast if there's nothing to do...
   */
  status = 0;
  // eqn   = R_LUBP;  //PRS: NEED TO DO SOMETHING HERE
  eqn   = EQN;  
  if (! pd->e[eqn]) return(status);

  /*
   * Load Gauss point weights before looking for friends
   */
  dbl dim = pd->Num_Dim;
  dbl wt = fv->wt;
  dbl h3 = fv->h3;

  /*
   * Prepare geometry
   */
  n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
  lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);
                                                                 
  /* Load proper FEM weights */
  dbl det_J = fv->sdet;            

  /* Load up source models -- momentum*/
  // err = load_lubrication_momentum_source(time, dt);

  /* Load up source models -- mass */
  // No calls yet as only constat models exist. See mm_input_mp.c 

  /* Time settings */
  if(pd->TimeIntegration != TRANSIENT) {
    tt = -0.5; dt = 1.0;
  }

  /*** CALCULATE FLOW RATE FROM FUNCTION **************************************/
  calculate_lub_q_v(EQN,  time, dt, xi, exo); //PRS: NEED TO DO SOMETHING HERE


  /*** CALCULATE PHYSICAL PROPERTIES AND SENSITIVITIES ************************/

  /* Lubrication height from model */
  H = height_function_model(&H_U, &dH_U_dtime, &H_L, &dH_L_dtime, dH_U_dX, dH_L_dX, &dH_U_dp, &dH_U_ddh, time, dt); 
  dH_dtime = dH_U_dtime - dH_L_dtime;
  /*
  if (pd->v[SHELL_DELTAH] && 
      (mp->HeightUFunctionModel == CONSTANT_SPEED_DEFORM ||
       mp->HeightUFunctionModel == CONSTANT_SPEED_MELT ||
       mp->HeightUFunctionModel == FLAT_GRAD_FLAT_MELT ||
       mp->HeightUFunctionModel == CIRCLE_MELT )) toggle_dh_dependence = 1.;
  */

  /* Deform lubrication height for FSI interaction */
  switch ( mp->FSIModel ) {
  case FSI_MESH_CONTINUUM:
  case FSI_MESH_UNDEF:
    for ( i = 0; i < dim; i++) {
      H -= fv->snormal[i] * fv->d[i];
      if ( pd->TimeIntegration == TRANSIENT ) {
	dH_dtime -= fv->snormal[i] * fv_dot->d[i];
      }
    }
    break;

  case FSI_REALSOLID_CONTINUUM:
    for ( i = 0; i < dim; i++) {
      H -= fv->snormal[i] * fv->d_rs[i];
      if ( pd->TimeIntegration == TRANSIENT ) {
	dH_dtime -= fv->snormal[i] * fv_dot->d_rs[i];
      }
    }  
    break;
  }

  /* Check for nehative lubrication height, if so, get out */
  if(H <= 0.0)
   {
    neg_lub_height = TRUE;

#ifdef PARALLEL
    fprintf(stderr,"\nP_%d: Lubrication height =  %e\n",ProcID,H);
#else
    fprintf(stderr,"\n Lubrication height =  %e\n",H);
#endif

    status = 2;
    return(status);
   }


  /* Lubrication wall velocity from model */
  velocity_function_model(veloU, veloL, time, dt);

  /* Lubrication height - mesh sensitivity */
  memset(dH_dtime_dmesh, 0.0, sizeof(double)*DIM*MDE);
  memset(dH_dtime_drealsolid, 0.0, sizeof(double)*DIM*MDE);
  switch ( mp->FSIModel ) {
  case FSI_MESH_CONTINUUM:
  case FSI_MESH_UNDEF:
    for ( i = 0; i < VIM; i++) {
      for ( b = 0; b < dim; b++) {
	for ( k = 0; k < ei->dof[MESH_DISPLACEMENT1]; k++) {
	  jk = dof_map[k];
	  dH_dtime_dmesh[b][k] -= fv->dsnormal_dx[i][b][jk] * fv_dot->d[i];
	  dH_dtime_dmesh[b][k] -= fv->snormal[i] * delta(i,b) * bf[MESH_DISPLACEMENT1]->phi[k] * (1+2*tt)/dt;
	}
      }
    }
    break;
  case FSI_REALSOLID_CONTINUUM:
    for ( i = 0; i < VIM; i++) {
      for ( b = 0; b < dim; b++) {
	for ( k = 0; k < ei->dof[MESH_DISPLACEMENT1]; k++) {
	  jk = dof_map[k];
	  dH_dtime_dmesh[b][k] -= fv->dsnormal_dx[i][b][jk] * fv_dot->d_rs[i];
	}
	for ( k = 0; k < ei->dof[SOLID_DISPLACEMENT1]; k++) {
	  jk = dof_map[k];
	  dH_dtime_drealsolid[b][k] -= fv->snormal[i] * delta(i,b) * bf[SOLID_DISPLACEMENT1]->phi[jk] * (1+2*tt)/dt;

	}
      }
    }
    break;
  }

  /*** RESIDUAL ASSEMBLY ******************************************************/
  if (af->Assemble_Residual) {
    // eqn = R_LUBP; //PRS: NEED TO DO SOMETHING HERE
    eqn = EQN;
    peqn = upd->ep[eqn];
    
    /*** Loop over DOFs (i) ***/
    for ( i=0; i<ei->dof[eqn]; i++) {
      
      /* Prepare basis funcitons */
      ShellBF( eqn, i, &phi_i, grad_phi_i, grad_II_phi_i, d_grad_II_phi_i_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
      
      /* Assemble diffusion term */
      diffusion = 0.0;
      if (pd->e[eqn] & T_DIFFUSION) {
	for (p = 0; p < dim; p++) {
	  diffusion += LubAux->q[p] * grad_II_phi_i[p];
	}
	diffusion *= det_J * wt * h3 * pd->etm[eqn][(LOG2_DIFFUSION)];
      } 
      
      /* Assemble source term */
      source = 0.0;
      if ( pd->e[eqn] & T_SOURCE ) {
	source = (mp->lubsource);
	source += -dH_dtime;     
	source += veloU[0]*dH_U_dX[0] + veloU[1]*dH_U_dX[1];
	source -= veloL[0]*dH_L_dX[0] + veloL[1]*dH_L_dX[1];
	source *= phi_i;
      }
      source *= det_J * wt * h3 * pd->etm[eqn][(LOG2_SOURCE)];
      
      lec->R[peqn][i] += diffusion + source;  
    } /* end of loop over i */
  } /* end of Assemble_Residual */
  

  /*** JACOBIAN ASSEMBLY ******************************************************/

  if (af->Assemble_Jacobian) {
    // eqn   = R_LUBP; //PRS: NEED TO DO SOMETHING HERE
    eqn   = EQN;
    peqn = upd->ep[eqn];

    /*** Loop over DOFs (i) ***/
    for ( i=0; i<ei->dof[eqn]; i++) {

      /* Prepare basis functions (i) */
      ShellBF( eqn, i, &phi_i, grad_phi_i, grad_II_phi_i, d_grad_II_phi_i_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );

      /* 
       * J_lubp_p or J_lubp2_p2  --the diagonal piece. 
       */
      if(EQN==R_LUBP) 
	{
	  var = LUBP;
	}
      else if (EQN==R_LUBP_2) 
	{
	  var = LUBP_2;
	}
      else EH(-1,"Mucho problema: Shouldn't be here.");

      if ( pd->v[var] ) {
	pvar = upd->vp[var];
	
	/*** Loop over DOFs (j) ***/
	for ( j=0; j<ei->dof[var]; j++) {

	  /* Load basis functions (j) */
	  ShellBF( var, j, &phi_j, grad_phi_j, grad_II_phi_j, d_grad_II_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
	  
	  /* Add diffusion term */
	  diffusion = 0.0;
	  if ( pd->e[eqn] && T_DIFFUSION ) {
	    for ( a = 0; a < dim; a++) {
	      diffusion += LubAux->dq_dp1[a][j] * grad_II_phi_j[a] * grad_II_phi_i[a];
	      diffusion += LubAux->dq_dp2[a][j] * phi_j * grad_II_phi_i[a];
	    }
	  }
	  diffusion *= det_J * wt * h3 * pd->etm[eqn][(LOG2_DIFFUSION)];

	  lec->J[peqn][pvar][i][j] += diffusion;
	} // End of loop over j
      } // End of J_lubp_p


      /* 
       * J_lubp_curv
       */
      var = SHELL_LUB_CURV;
      if ( pd->v[var] ) {
	pvar = upd->vp[var];
	
	/*** Loop over DOFs (j) ***/
	for ( j=0; j<ei->dof[var]; j++) {

	  /* Load basis functions (j) */
	  ShellBF( var, j, &phi_j, grad_phi_j, grad_II_phi_j, d_grad_II_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );

	  /* Add diffusion term */
	  diffusion = 0.0;
	  if ( pd->e[eqn] && T_DIFFUSION ) {
	    for ( b = 0; b < dim; b++ ) {
	      diffusion += LubAux->dq_dk[b][j] * grad_II_phi_i[b] * phi_j;
	    }
	  }
	  diffusion *= det_J * wt * h3 * pd->etm[eqn][(LOG2_DIFFUSION)];

	  lec->J[peqn][pvar][i][j] += diffusion;
	} // End of loop over j
      } // End of J_lubp_curv

     /* 
       * J_lubp_curv_2
       */
      var = SHELL_LUB_CURV_2;
      if ( pd->v[var] ) {
	pvar = upd->vp[var];
	
	/*** Loop over DOFs (j) ***/
	for ( j=0; j<ei->dof[var]; j++) {

	  /* Load basis functions (j) */
	  ShellBF( var, j, &phi_j, grad_phi_j, grad_II_phi_j, d_grad_II_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );

	  /* Add diffusion term */
	  diffusion = 0.0;
	  if ( pd->e[eqn] && T_DIFFUSION ) {
	    for ( b = 0; b < dim; b++ ) {
	      diffusion += LubAux->dq_dk[b][j] * grad_II_phi_i[b] * phi_j;
	    }
	  }
	  diffusion *= det_J * wt * h3 * pd->etm[eqn][(LOG2_DIFFUSION)];

	  lec->J[peqn][pvar][i][j] += diffusion;
	} // End of loop over j
      } // End of J_lubp_curv_2


      /* 
       * J_lubp_LS or J_lubp_phase1  depending on lubp or lubp2
       */
      var = LS;
      if(EQN == R_LUBP_2)  var = PHASE1;

      if ( pd->v[var] ) {
	pvar = upd->vp[var];
	
	/*** Loop over DOFs (j) ***/
	for ( j=0; j<ei->dof[var]; j++) {

	  /* Load basis functions (j) */
	  ShellBF( var, j, &phi_j, grad_phi_j, grad_II_phi_j, d_grad_II_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );

	  /* Add diffusion term */
	  diffusion = 0.0;
	  if ( pd->e[eqn] && T_DIFFUSION ) {
	    for ( b = 0; b < dim; b++ ) {
	      diffusion += LubAux->dq_df[b][j] * grad_II_phi_i[b];
	    }
	  }
	  diffusion *= det_J * wt * h3 * pd->etm[eqn][(LOG2_DIFFUSION)];

	  lec->J[peqn][pvar][i][j] += diffusion;
	} // End of loop over j
      } // End of J_lubp_LS




      /* 
       * J_lubp_DMX
       */
      var = MESH_DISPLACEMENT1;
      if (pd->v[var] && (
			  mp->FSIModel == FSI_MESH_CONTINUUM ||
			  mp->FSIModel == FSI_REALSOLID_CONTINUUM ||
			  mp->FSIModel == FSI_MESH_UNDEF
			  )) {
	pvar = upd->vp[var];

	/*** Loop over dimensions of mesh displacement ***/
	for ( b = 0; b < dim; b++) {
	  var = MESH_DISPLACEMENT1 + b;
	  pvar = upd->vp[var];

	  /*** Loop over DOFs (j) ***/
	  for ( j=0; j < ei->dof[var]; j++) {
	    jk = dof_map[j];
	    
	    /* Load basis functions (j) */
	    ShellBF( eqn, j, &phi_j, grad_phi_j, grad_II_phi_j, d_grad_II_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
	    
	    /* Add diffusion term */
	    diffusion = 0.0;
	    if ( pd->e[eqn] & T_DIFFUSION ) {
	      for ( p = 0; p < dim; p++) {
		diffusion += det_J * LubAux->dq_dx[p][b][j] * grad_II_phi_i[p];
		diffusion += det_J * LubAux->q[p] * d_grad_II_phi_i_dmesh[p][b][jk];
		diffusion += fv->dsurfdet_dx[b][jk] * LubAux->q[p] * grad_II_phi_i[p];
	      }
	    }
	    diffusion *= wt * h3 * pd->etm[eqn][(LOG2_DIFFUSION)];

	    /* Add source term */
	    source = 0.0;
	    if ( pd->e[eqn] & T_SOURCE ) {
	      source += -dH_dtime_dmesh[b][j] * det_J;
	      source += (mp->lubsource - dH_dtime) * fv->dsurfdet_dx[b][jk];
	      source += (veloU[0]*dH_U_dX[0] + veloU[1]*dH_U_dX[1])*fv->dsurfdet_dx[b][jk];
	      source -= (veloL[0]*dH_L_dX[0] + veloL[1]*dH_L_dX[1])*fv->dsurfdet_dx[b][jk];
	      source *= phi_i;
	    }
	    source *= wt * h3 * pd->etm[eqn][(LOG2_SOURCE)];
	    
	    lec->J[peqn][pvar][i][jk] += diffusion + source;
	  } // End of loop over j
	} // End of loop over b
      } // End of J_lubp_mesh

      /* 
       * J_lubp_DRS
       */
      var = SOLID_DISPLACEMENT1;
      if ( upd->vp[var] >= 0 && (
			  mp->FSIModel == FSI_REALSOLID_CONTINUUM )) {
	pvar = upd->vp[var];

	/*** Loop over dimensions of mesh displacement ***/
	for ( b = 0; b < dim; b++) {
	  var = SOLID_DISPLACEMENT1 + b;
	  pvar = upd->vp[var];

	  /*** Loop over DOFs (j) ***/
	  for ( j=0; j < ei->dof[var]; j++) {
	    jk = dof_map[j];
	    
	    /* Load basis functions (j) */
	    ShellBF( eqn, j, &phi_j, grad_phi_j, grad_II_phi_j, d_grad_II_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
	    
	    /* Add diffusion term */
	    diffusion = 0.0;
	    if ( pd->e[eqn] & T_DIFFUSION ) {
	      for ( p = 0; p < dim; p++) {
		diffusion += det_J * LubAux->dq_drs[p][b][j] * grad_II_phi_i[p];
	      }
	    }
	    diffusion *= wt * h3 * pd->etm[eqn][(LOG2_DIFFUSION)];

	    /* Add source term */
	    source = 0.0;
	    if ( pd->e[eqn] & T_SOURCE ) {
	      source += -dH_dtime_drealsolid[b][j] * det_J;
	      source *= phi_i;
	    }
	    source *= wt * h3 * pd->etm[eqn][(LOG2_SOURCE)];
	    
	    lec->J[peqn][pvar][i][jk] += diffusion + source;
	  } // End of loop over j
	} // End of loop over b
      } // End of J_lubp_drs


      /* 
       * J_lubp_D_sh_dh
       */
      var = SHELL_DELTAH;
      if ( pd->v[var]) {
	pvar = upd->vp[var];
	
	/*** Loop over DOFs (j) ***/
	for ( j=0; j < ei->dof[var]; j++) {
	  
	  /* Load basis functions (j) */
	  ShellBF( var, j, &phi_j, grad_phi_j, grad_II_phi_j, d_grad_II_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
	  
	  /* Add diffusion term */
	  diffusion = 0.0;
	  if ( pd->e[eqn] & T_DIFFUSION ) {
	    for ( p = 0; p < dim; p++) {
	      diffusion += det_J * LubAux->dq_ddh[p][j] * phi_j * grad_II_phi_i[p];
	    }
	  }
	  diffusion *= wt * h3 * pd->etm[eqn][(LOG2_DIFFUSION)];
	  
	  /* Add source term */
	  source = 0.0;
	  if ( pd->e[eqn] & T_SOURCE ) {
	    // dh_time no longer has dependence here, as of 4/11/2011. Talk to PRS. 
	    // If you wanted to add some volume expansion, however, there would be
	    // a boost here. 
	    // source += -0.*toggle_dh_dependence*(1 + 2. * tt)*phi_j/dt;
	     source *= phi_i;
	  }
	  source *= det_J * wt * h3 * pd->etm[eqn][(LOG2_SOURCE)];
	  
	  lec->J[peqn][pvar][i][j] += diffusion + source;
	} // End of loop over j
      } // End of J_lubp_dDeltah

      /* 
       * J_lubp_D_sh_pc
       */

      var = SHELL_PARTC;

      if (pd->v[var])
	{
	 pvar = upd->vp[var];
	 for ( j=0; j<ei->dof[var]; j++)
            {
	     phi_j = bf[var]->phi[j];

             diffusion = 0.;
	     if (pd->e[eqn] & T_DIFFUSION)
	       {
		for (p = 0; p < VIM; p++) 
	           {
	            diffusion += LubAux->dq_dc[p][j] * phi_j * grad_II_phi_i[p];
                   }
		    
		 diffusion *= det_J * wt;
		 diffusion *= h3;
		 diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
               }
	     lec->J[peqn][pvar][i][j] += diffusion;
            }// End of loop over j
        }// End of J_lubp_dsh_pc


      
    } /* end of loop over i */
  } /* end of Assemble_Jacobian */
  
    /* clean-up */
  fv->wt = wt;  /* load_neighbor_var_data screws this up */
  safe_free((void *) n_dof);
  return(status);
} /* end of assemble_lubrication */

/*****************************************************************************/
/***assemble_shell_energy******************************************************/
/*  _______________________________________________________________________  */

/* assemble_shell_energy -- assemble terms (Residual & Jacobian) for scalar shell
 *                          energy equation, assuming constant temperature across shell. 
 *
 * in:
 * 	ei -- pointer to Element Indices	structure
 *	pd -- pointer to Problem Description	structure
 *	af -- pointer to Action Flag		structure
 *	bf -- pointer to Basis Function		structure
 *	fv -- pointer to Field Variable		structure
 *  fv_old -- pointer to old Diet Field Variable	structure
 *  fv_dot -- pointer to dot Diet Field Variable	structure
 *	cr -- pointer to Constitutive Relation	structure
 *	md -- pointer to Mesh Derivative	structure
 *	me -- pointer to Material Entity	structure
 *
 * out:
 *	a   -- gets loaded up with proper contribution 
 *	lec -- gets loaded up with local contributions to resid, Jacobian
 * 	r   -- residual RHS vector
 *
 * Created:	Friday August 13, 2010 prschun@sandia.gov
 *
 */
/*ARGSUSED*/
int 
assemble_shell_energy(double time,	/* present time value */ 
		      double tt,   	/* parameter to vary time integration from 
			   	         * explicit (tt = 1) to implicit (tt = 0)    */
		      double dt,         /* current time step size */
		      double xi[DIM],    /* Local stu coordinates */
		      const PG_DATA *pg_data, /*Upwinding stuff */
		      const Exo_DB *exo)  
{
  int eqn, var, peqn, pvar, dim, p, a, b, k, jk;
  int *n_dof = NULL;
  int dof_map[MDE];
  int i = -1, ii;
  int j, jj, status;
  
  dbl curv = 0, H, dH_dtime;    /* Temperature derivative of viscosity */
  dbl H_U, dH_U_dtime, H_L, dH_L_dtime;
  dbl dH_U_dX[DIM],dH_L_dX[DIM], dH_dmesh[DIM][MDE], dH_dtime_dmesh[DIM][MDE];
  dbl dH_drealsolid[DIM][MDE], dH_dtime_drealsolid[DIM][MDE];
  dbl dH_U_dp, dH_U_ddh;

  dbl q_tot = 0.;
  dbl dqdh[DIM];

  dbl grad_P[DIM];  /* Lubrication pressure gradient. */
  dbl grad_T[DIM];  /* Shell temperature gradient. */
  dbl q[DIM], dqdp[DIM], dqdp_1[DIM], dq_lub_dH[DIM], dqdmu[DIM], dqdF[DIM][MDE], dq_dmesh[DIM][DIM][MDE];
  dbl dq_drs[DIM][DIM][MDE];

  dbl veloU[DIM], veloL[DIM];
  dbl  mu, rho, t_cond, Cp;  
  dbl k_eff, k_turb, d_k_turb_dmu, d_k_eff_dh;  /* sundry factors*/
  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct;  /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;

  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct; /* Thermal conductivity dependence. */
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  HEAT_CAPACITY_DEPENDENCE_STRUCT d_Cp_struct; /* Heat capacity dependence. */
  HEAT_CAPACITY_DEPENDENCE_STRUCT *d_Cp = &d_Cp_struct;

 /* density derivatives */
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;  /* density dependence */
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;

  const double *hsquared = pg_data->hsquared ;
  const double *vcent = pg_data->v_avg; /* Average element velocity, which is the
	  centroid velocity for Q2 and the average
	  of the vertices for Q1. It comes from
	  the routine "element_velocity." */


  dbl mass, advection, advection_b, diffusion,  source ;

  /*
   * Galerkin weighting functions for i-th shell residuals 
   * and some of their derivatives...
   */
  dbl phi_i, phi_j;
  dbl grad_II_phi_i[DIM];
  dbl grad_phi_j[DIM], grad_II_phi_j[DIM];
  dbl d_grad_T_dmesh[DIM][DIM][MDE];
  dbl d_grad_P_dmesh[DIM][DIM][MDE];
  dbl d_grad_II_phi_i_dmesh[DIM][DIM][MDE];

  dbl gradII_Hside[DIM], gradII_Hside_F[DIM][MDE];

  /*
   * Petrov-Galerkin weighting functions for i-th residuals
   * and some of their derivatives...
   */

  dbl wt_func;

/* SUPG variables */
  dbl h_elem=0, h_elem_inv=0, h_elem_deriv=0., h_elem_inv_deriv=0.;
  dbl supg = 0, d_wt_func;


  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl h3;			/* Volume element (scale factors). */
  dbl det_J;
  dbl wt;

  /*   static char yo[] = "assemble_shell_energy";*/
  status = 0;

  /*
   * Bail out fast if there's nothing to do...
   */
  eqn   = R_SHELL_ENERGY;
  if (! pd->e[eqn]) return(status);

  /*
   * Load Gauss point weights before looking for friends
   */
  dim = pd->Num_Dim;
  wt = fv->wt;		 	/* Gauss point weight. */
  h3 = fv->h3;			/* Differential volume element. */

  /* 
   * More upwinding stuff
   */

  if( mp->Ewt_funcModel == GALERKIN)
    {
      supg = 0.;
    }
  else if( mp->Ewt_funcModel == SUPG)
    {
      if( !pd->e[R_MOMENTUM1]) 
	EH(-1, " must have momentum equation velocity field for shell_energy upwinding");
      supg = mp->Ewt_func;
    }

  /*
   * Prepare geometry
   */
  n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
  /* NOT SURE THIS IS NEEDED like in assemble lubrication */
  lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);
                                                                 
  /*
   * Load proper FEM weights
   */ 
  det_J = fv->sdet;            


  /* 
   * Load material property constants 
   */ 
  mu = viscosity(gn, NULL, d_mu); 

  rho = density(d_rho, time);

  t_cond = conductivity( d_k, time );

  Cp = heat_capacity( d_Cp, time );


  /*
   * Time settings
   */
  if(pd->TimeIntegration != TRANSIENT) {
    tt = -0.5; dt = 1.0;
  }

  /*** CALCULATE PHYSICAL PROPERTIES AND SENSITIVITIES ************************/

  /* Lubrication height from model */
  H = height_function_model(&H_U, &dH_U_dtime, &H_L, &dH_L_dtime, dH_U_dX, dH_L_dX, &dH_U_dp, &dH_U_ddh, time, dt); 
  dH_dtime = dH_U_dtime - dH_L_dtime;

  /* Deform lubrication height for FSI interaction */
  switch ( mp->FSIModel ) {
  case FSI_MESH_CONTINUUM:
  case FSI_MESH_UNDEF:
    for ( i = 0; i < dim; i++) {
      H -= fv->snormal[i] * fv->d[i];
      if ( pd->TimeIntegration == TRANSIENT ) {
	dH_dtime -= fv->snormal[i] * fv_dot->d[i];
      }
    }
    break;
  case FSI_REALSOLID_CONTINUUM:
    for ( i = 0; i < dim; i++) {
      H -= fv->snormal[i] * fv->d_rs[i];
      if ( pd->TimeIntegration == TRANSIENT ) {
	dH_dtime -= fv->snormal[i] * fv_dot->d_rs[i];
      }
    }
    break;
  }

 /* Lubrication wall velocity from model */
  velocity_function_model(veloU, veloL, time, dt);

  k_eff = t_cond;  /* default laminar case */
  k_turb = 12.;
  d_k_turb_dmu = 0.;
  d_k_eff_dh = 0.;


  /* Deform wall slope for FSI interaction */

  /* Lubrication height - mesh sensitivity */
  memset(dH_dmesh, 0.0, sizeof(double)*DIM*MDE);
  memset(dH_dtime_dmesh, 0.0, sizeof(double)*DIM*MDE);
  switch ( mp->FSIModel ) {
  case FSI_MESH_CONTINUUM:
  case FSI_MESH_UNDEF:
    for ( i = 0; i < VIM; i++) {
      for ( b = 0; b < dim; b++) {
	for ( k = 0; k < ei->dof[MESH_DISPLACEMENT1]; k++) {
	  jk = dof_map[k];
	  dH_dmesh[b][jk] -= fv->dsnormal_dx[i][b][jk] * fv->d[i]
	    + delta(i,b) * fv->snormal[i] * bf[MESH_DISPLACEMENT1]->phi[k];
	  dH_dtime_dmesh[b][jk] -= fv->dsnormal_dx[i][b][k] * fv_dot->d[i]
	    + delta(i,b) * fv->snormal[i] * bf[MESH_DISPLACEMENT1]->phi[k] * (1+2*tt)/dt;
	}
      }
    }
    break;
  case FSI_REALSOLID_CONTINUUM:
    memset(dH_drealsolid, 0.0, sizeof(double)*DIM*MDE);
    memset(dH_dtime_drealsolid, 0.0, sizeof(double)*DIM*MDE);
    for ( i = 0; i < VIM; i++) {
      for ( b = 0; b < dim; b++) {
	for ( k = 0; k < ei->dof[MESH_DISPLACEMENT1]; k++) {
	  jk = dof_map[k];
	  dH_dmesh[b][k] -= fv->dsnormal_dx[i][b][jk] * fv->d_rs[i];
	  dH_dtime_dmesh[b][k] -= fv->dsnormal_dx[i][b][jk] * fv_dot->d_rs[i];
	  }
	for ( k = 0; k < ei->dof[SOLID_DISPLACEMENT1]; k++) {
	  jk = dof_map[k];
	  dH_drealsolid[b][k] -= delta(i,b) * fv->snormal[i] * bf[SOLID_DISPLACEMENT1]->phi[jk];
	  dH_dtime_drealsolid[b][k] -= delta(i,b) * fv->snormal[i] * bf[SOLID_DISPLACEMENT1]->phi[jk] * (1+2*tt)/dt;
	  }
      }
    }
    break; 
  }
  

  double if_liquid = 1.0; /*bolean to determine if in liquid phase (F<0) or not */
  if(pd->v[FILL])
    {
      load_lsi (ls->Length_Scale);
      load_lsi_derivs();
      if_liquid = (1.-lsi->H);
    }

  /* Temperature gradient */
  for ( i = 0; i < VIM; i++) {
    grad_T[i] = 0;
    for ( j = 0; j < VIM; j++) {
      grad_T[i] += ( fv->grad_sh_t[j] * delta(i,j) -
		     fv->grad_sh_t[j] * ( fv->snormal[i] * fv->snormal[j] ) );
    }
  }

  /* Temperature gradient - mesh sensitivity */
  memset(d_grad_T_dmesh, 0.0, sizeof(double)*DIM*DIM*MDE);
  memset(d_grad_P_dmesh, 0.0, sizeof(double)*DIM*DIM*MDE);
  switch ( mp->FSIModel ) {
  case FSI_MESH_CONTINUUM:
  case FSI_REALSOLID_CONTINUUM:
    for ( i = 0; i < VIM; i++) {
      for ( j = 0; j < VIM; j++) {
	for ( b = 0; b < dim; b++) {
	  for ( k = 0; k < ei->dof[MESH_DISPLACEMENT1]; k++) {
	    jk = dof_map[k];
	    d_grad_T_dmesh[i][b][jk] += ( fv->d_grad_sh_t_dmesh[j][b][k] * delta(i,j) -
					 fv->d_grad_sh_t_dmesh[j][b][k] * ( fv->snormal[i] * fv->snormal[j] ) )
	      - fv->dsnormal_dx[i][b][jk] * fv->snormal[j] * fv->grad_sh_t[j]
	      - fv->snormal[i] * fv->dsnormal_dx[j][b][jk] * fv->grad_sh_t[j];
	     if (pd->e[R_LUBP])
	       {
		 d_grad_P_dmesh[i][b][jk] += ( fv->d_grad_lubp_dmesh[j][b][k] * delta(i,j) -
					      fv->d_grad_lubp_dmesh[j][b][k] * ( fv->snormal[i] * fv->snormal[j] ) )
		   - fv->dsnormal_dx[i][b][jk] * fv->snormal[j] * fv->grad_lubp[j]
		   - fv->snormal[i] * fv->dsnormal_dx[j][b][jk] * fv->grad_lubp[j];
	       }
	  }
	}
      }
    }
    break;
  }
  
  /* velocity field from flow rate (normal q divided by H */

  if (pd->e[R_LUBP])
    {
      for ( i = 0; i < VIM; i++) 
	{
	  grad_P[i] = 0;
	  for ( j = 0; j < VIM; j++) 
	    {
	      grad_P[i] += ( fv->grad_lubp[j] * delta(i,j) -
			     fv->grad_lubp[j] * ( fv->snormal[i] * fv->snormal[j] ) );
	    }
	}
      for ( i = 0; i < dim; i++) 
	{
	  q[i]  = 0.0;
	  q[i] -= pow(H,3)*grad_P[i]/(k_turb*mu);
	  q[i] += H*(veloL[i]+veloU[i])/2;
	  if (pd->v[FILL]) q[i] -= pow(H,3)*gradII_Hside[i]*curv*mp->surface_tension/(k_turb*mu);
	}

      if(supg!=0.)
	{
	  h_elem = 0.;
	  for ( p=0; p<dim; p++)
	    {
	      if(hsquared[p] != 0.) h_elem += vcent[p] * vcent[p] / hsquared[p];
	    }
	  h_elem = sqrt(h_elem)/2.;
	  if(h_elem == 0.)
	    {
	      h_elem_inv=0.;
	    }
	  else
	    {
	      h_elem_inv=1./h_elem;
	    }

	}

      /* Flow rate - mesh sensitivity */
      memset(dq_dmesh, 0.0, sizeof(double)*DIM*DIM*MDE);
      memset(dq_drs, 0.0, sizeof(double)*DIM*DIM*MDE);
      switch ( mp->FSIModel ) {
      case FSI_MESH_CONTINUUM:
      case FSI_MESH_UNDEF:
	for ( i = 0; i < VIM; i++) {
	  for ( b = 0; b < dim; b++) {
	    for ( k = 0; k < n_dof[MESH_DISPLACEMENT1]; k++) {
	      dq_dmesh[i][b][k] += -3.0 * pow(H,2) * dH_dmesh[b][k] * grad_P[i] / (k_turb*mu)
		- pow(H,3) * d_grad_P_dmesh[i][b][k] / (k_turb*mu)
		+ dH_dmesh[b][k] * (veloL[i]+veloU[i]) / 2.0;
	    }
	  }
	}
	break;

      case FSI_REALSOLID_CONTINUUM:
	for ( i = 0; i < VIM; i++) {
	  for ( b = 0; b < dim; b++) {
	    for ( k = 0; k < n_dof[MESH_DISPLACEMENT1]; k++) {
	      dq_dmesh[i][b][k] += -3.0 * pow(H,2) * dH_dmesh[b][k] * grad_P[i] / (k_turb*mu)
		- pow(H,3) * d_grad_P_dmesh[i][b][k] / (k_turb*mu)
		+ dH_dmesh[b][k] * (veloL[i]+veloU[i]) / 2.0;

	      dq_drs[i][b][k] += -3.0 * pow(H,2) * dH_drealsolid[b][k] * grad_P[i] / (k_turb*mu)
				+ dH_drealsolid[b][k] * (veloL[i]+veloU[i]) / 2.0;
	    }
	  }
	}
	break;
      }

      /* Flow rate - parameter sensitivities */
      memset(dqdF, 0.0, sizeof(double)*DIM*MDE);
      memset(dq_lub_dH, 0.0, sizeof(double)*DIM);
      memset(dqdp_1, 0.0, sizeof(double)*DIM);
      for ( i = 0; i < dim; i++) {
	dqdp[i] = -pow(H,3)/(k_turb*mu);
	/*this next piece is to pick up dependence on H on lubp from the deform spring equation */
	dqdp_1[i]=-3.*pow(H,2)*dH_U_dp*grad_P[i]/(k_turb*mu) + dH_U_dp*(veloL[i]+veloU[i])/2.;
	if (pd->v[FILL]) dqdp_1[i] -= 3.*pow(H,2)*dH_U_dp*gradII_Hside[i]*curv*mp->surface_tension/(k_turb*mu);

	dqdmu[i] = pow(H,3)*grad_P[i]*(d_k_turb_dmu*mu + k_turb)/(k_turb*k_turb*mu*mu); 
	if (pd->v[FILL]) dqdmu[i] += pow(H,3)*gradII_Hside[i]*curv*mp->surface_tension*(d_k_turb_dmu*mu + k_turb)/(k_turb*k_turb*mu*mu);
	if (pd->v[SHELL_DELTAH] && 
	    (mp->HeightUFunctionModel == CONSTANT_SPEED_DEFORM ||
	     mp->HeightUFunctionModel == CONSTANT_SPEED_MELT ||
	     mp->HeightUFunctionModel == FLAT_GRAD_FLAT_MELT ||
	     mp->HeightUFunctionModel == CIRCLE_MELT)) dq_lub_dH[i] += -3*pow(H,2)*grad_P[i]/(k_turb*mu) + (veloL[i]+veloU[i])/2.; 
	if (pd->v[FILL])
	  {
	    if (pd->v[SHELL_DELTAH]) dqdh[i] -= 3.*pow(H,2)*gradII_Hside[i]*curv*mp->surface_tension/(k_turb*mu);
	    for ( j = 0; j < ei->dof[eqn]; j++ ) {
	      dqdF[i][j]  = -pow(H,3)*gradII_Hside_F[i][j]*curv*mp->surface_tension/(k_turb*mu);
	      
	      dqdF[i][j] +=  dqdmu[i]*d_mu->F[j];

	    }
	  }
      }
    }
   else 
     {
       memset(dq_dmesh, 0.0, sizeof(double)*DIM*DIM*MDE);
       memset(dq_drs, 0.0, sizeof(double)*DIM*DIM*MDE);
       memset(q, 0.0, sizeof(double)*DIM);
     }

     
 /* Finall, Load up all heat source models */

  /* Note all of this depend on H_lub.  You need to add those sensitivities */ 
  /* No source terms available right now. Feel free to add your own.  */
  q_tot = 0.;



  /*** RESIDUAL ASSEMBLY ******************************************************/
  if (af->Assemble_Residual) {
    eqn = R_SHELL_ENERGY;
    peqn = upd->ep[eqn];
    
    /*** Loop over DOFs (i) ***/
    for ( i=0; i<ei->dof[eqn]; i++) {
      
      /* Prepare basis funcitons */
      phi_i = bf[eqn]->phi[i];
      for ( ii = 0; ii < VIM; ii++) {
	grad_II_phi_i[ii] = 0.0;
	for ( jj = 0; jj < VIM; jj++) {
	  grad_II_phi_i[ii] += ( bf[eqn]->grad_phi[i][jj] * delta(ii,jj) -
				 bf[eqn]->grad_phi[i][jj] * ( fv->snormal[ii] * fv->snormal[jj] ) );
	}
      }

      /* Mass term */
      mass = 0.;
      if ( pd->TimeIntegration != STEADY )
	{
	  if ( pd->e[eqn] & T_MASS )
	    {
	      mass  = fv_dot->sh_t;
	      mass *= - H * phi_i * rho * Cp * det_J * wt;
	      mass *= h3;
	      mass *= pd->etm[eqn][(LOG2_MASS)];
	    }
	}

      advection = 0.;
      wt_func = bf[eqn]->phi[i];
      
/* add Petrov-Galerkin terms as necessary */
	  if(supg!=0.)
	    {
	      for(p=0; p<dim; p++)
		{
		  wt_func += supg * h_elem_inv * (q[p]/H) * grad_II_phi_i[p] ;
		}
	    }
	  if ( pd->e[eqn] & T_ADVECTION )
	    {

	      for ( p=0; p<VIM; p++)
		{
		  advection += q[p]     * grad_T[p];
		}

	      advection *= - rho * Cp * det_J * wt_func * wt;
	      advection *= h3;
	      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
	    }
      
	  /* Assemble diffusion term */
	  diffusion = 0.0;
	  if (pd->e[eqn] & T_DIFFUSION) {
	
	    for (p = 0; p < dim; p++) {
	      diffusion += H * k_eff * grad_T[p] * grad_II_phi_i[p];
	    }
	
	    diffusion *= det_J * wt * h3 * pd->etm[eqn][(LOG2_DIFFUSION)];
	  } 
      
	  /* Assemble source term */
	  source = 0.0;
	  if ( pd->e[eqn] & T_SOURCE ) {
	    source = q_tot;
	    source *= phi_i;
	  }
	  source *= if_liquid * det_J * wt * h3 * pd->etm[eqn][(LOG2_SOURCE)];
      
	  lec->R[peqn][i] += mass + advection + diffusion + source;  
    } /* end of loop over i */
  } /* end of Assemble_Residual */
  

  /*** JACOBIAN ASSEMBLY ******************************************************/

  if (af->Assemble_Jacobian) {
    eqn   = R_SHELL_ENERGY;
    peqn = upd->ep[eqn];

    /*** Loop over DOFs (i) ***/
    for ( i=0; i<ei->dof[eqn]; i++) {

      /* Prepare basis functions (i) */
      phi_i = bf[eqn]->phi[i];
      for ( ii = 0; ii < VIM; ii++) {
	grad_II_phi_i[ii] = 0.0;
	for ( jj = 0; jj < VIM; jj++) {
	  grad_II_phi_i[ii] += ( bf[eqn]->grad_phi[i][jj] * delta(ii,jj) -
				 bf[eqn]->grad_phi[i][jj] * ( fv->snormal[ii] * fv->snormal[jj] ) );
	}
      }

      /* Basis functions - mesh sensitivities */
      memset(d_grad_II_phi_i_dmesh, 0, sizeof(double)*DIM*DIM*MDE);
      switch ( mp->FSIModel ) {
      case FSI_MESH_CONTINUUM:
      case FSI_REALSOLID_CONTINUUM:
	for ( ii = 0; ii < VIM; ii++) {
	  for ( jj = 0; jj < VIM; jj++) {
	    for ( b = 0; b < dim; b++) {
	      for ( k = 0; k < ei->dof[MESH_DISPLACEMENT1]; k++) {
		jk = dof_map[k];
		d_grad_II_phi_i_dmesh[ii][b][jk] +=  bf[eqn]->d_grad_phi_dmesh[i][jj][b][k] * delta(ii,jj) 
		  - bf[eqn]->d_grad_phi_dmesh[i][jj][b][k] * ( fv->snormal[ii] * fv->snormal[jj] )
		  - fv->dsnormal_dx[ii][b][jk] * fv->snormal[jj] * bf[eqn]->grad_phi[i][jj]
		  - fv->snormal[ii] * fv->dsnormal_dx[jj][b][jk] * bf[eqn]->grad_phi[i][jj];
	      }
	    }
	  }
	}
	break;
      }
      wt_func = bf[eqn]->phi[i];
      /* add Petrov-Galerkin terms as necessary */
      if(supg!=0.)
	{
	  for(p=0; p<dim; p++)
	    {
	      wt_func += supg * h_elem_inv * (q[p]/H) * grad_II_phi_i[p];
	    }
	}

      /* 
       * J_shell_energy_d_sh_t 
       */
      var = SHELL_TEMPERATURE;
      if ( pd->v[var] ) {
	pvar = upd->vp[var];
	
	/*** Loop over DOFs (j) ***/
	for ( j=0; j<ei->dof[var]; j++) {

	  /* Load basis functions (j) */
	  phi_j = bf[eqn]->phi[j];
	  for ( ii = 0; ii < VIM; ii++) {
	    grad_phi_j[ii] = bf[eqn]->grad_phi[j][ii];
	    grad_II_phi_j[ii] = 0.0;
	    for ( jj = 0; jj < VIM; jj++) {
	      grad_II_phi_j[ii] += ( bf[eqn]->grad_phi[j][jj] * delta(ii,jj) -
				     bf[eqn]->grad_phi[j][jj] * ( fv->snormal[ii] * fv->snormal[jj] ) );
	    }
	  }
	  
	  /* Add mass term */
	  mass = 0.0;
	  if( pd->TimeIntegration != STEADY)
	    {
	      if ( pd->e[eqn] & T_MASS )
		{
		  mass =  (1 + 2. * tt) * phi_j / dt;
		  mass *= - H*phi_i * rho * Cp * det_J * wt;
		  mass *= h3 * pd->etm[eqn][(LOG2_MASS)];
		}
	    }

	  /* Add advection term */
	  advection = 0.;
	  if ( pd->e[eqn] & T_ADVECTION )
	    {

	      for ( a=0; a<VIM; a++)
		{
		  advection += q[a]     * grad_II_phi_j[a];
		}

	      advection *= - rho * Cp * det_J * wt * wt_func;
	      advection *= h3;
	      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
	    }
	
	  /* Add diffusion term */
	  diffusion = 0.0;
	  if ( pd->e[eqn] && T_DIFFUSION ) {
	    for ( a = 0; a < dim; a++) {
	      diffusion += H * k_eff *  grad_II_phi_i[a] * grad_II_phi_j[a];
	    }
	  }
	  diffusion *= det_J * wt * h3 * pd->etm[eqn][(LOG2_DIFFUSION)];

	  /* Add source term */
	  source = 0;
	    if ( pd->e[eqn] & T_SOURCE ) {
	      /* PRS Note these are initialized to zero in mm_input_mp.c */
	      source = 0.;  /* add your dq_dt sensitivities here */
	      source *= phi_i;
	    }
	  source *= if_liquid * det_J * wt * h3 * pd->etm[eqn][(LOG2_SOURCE)];

	  lec->J[peqn][pvar][i][j] += mass + advection + diffusion + source ;
	} // End of loop over j
      } // End of J_shell_energy_d_shell_temperature

      /* 
       * J_shell_energy_LS
       */
      var = LS;
      if ( pd->v[var] ) {
	pvar = upd->vp[var];
	
	/*** Loop over DOFs (j) ***/
	for ( j=0; j<ei->dof[var]; j++) {

	  /* Load basis functions (j) */
	  phi_j = bf[eqn]->phi[j];
	  for ( ii = 0; ii < VIM; ii++) {
	    grad_phi_j[ii] = bf[eqn]->grad_phi[j][ii];
	    grad_II_phi_j[ii] = 0.0;
	    for ( jj = 0; jj < VIM; jj++) {
	      grad_II_phi_j[ii] +=  ( grad_phi_j[jj] * delta(ii,jj) -
				      grad_phi_j[jj] * ( fv->snormal[ii] * fv->snormal[jj] ) );
	    }
	  }

	  WH(-1," Haven't added LS sensitivities to shell energy equation yet");

	} // End of loop over j
      } // End of J_lubp_LS



      /* 
       * J_shell_energy_dmx
       */
      var = MESH_DISPLACEMENT1;
      if ( pd->v[var] && (
			  mp->FSIModel == FSI_MESH_CONTINUUM ||
			  mp->FSIModel == FSI_REALSOLID_CONTINUUM ||
			  mp->FSIModel == FSI_MESH_UNDEF
			  )) {
	pvar = upd->vp[var];

	/*** Loop over dimensions of mesh displacement ***/
	for ( b = 0; b < dim; b++) {
	  var = MESH_DISPLACEMENT1 + b;
	  pvar = upd->vp[var];

	  /*** Loop over DOFs (j) ***/
	  /***PRS/SAR should be ei->dof */
	  /***And all d_mesh pieces like dsurfdet_dx etc. should all be jk=dof_map[j] ***/
	  for ( j=0; j < ei->dof[var]; j++) {
	    jk = dof_map[j];
	    
	    /* Load basis functions (j) */
	    phi_j = bf[eqn]->phi[j];
	    for ( ii = 0; ii < VIM; ii++) {
	      grad_phi_j[ii] = bf[eqn]->grad_phi[j][ii];
	      grad_II_phi_j[ii] = 0.0;
	      for ( jj = 0; jj < VIM; jj++) {
		grad_II_phi_j[ii] +=  ( grad_phi_j[jj] * delta(ii,jj) -
					grad_phi_j[jj] * ( fv->snormal[ii] * fv->snormal[jj] ) );
	      }
	    }

	    if(supg!=0.)
	      {
		h_elem_deriv = 0.;
		h_elem_inv_deriv=0.;
		for( p=0; p<dim; p++ )
		  {
		    if(pg_data->hhv[p][b] != 0.)
		      {
			h_elem_deriv -=
			  vcent[p] * vcent[p] *
			  pg_data->dhv_dxnode[p][j] * pg_data->hhv[p][b] * h_elem_inv/4./hsquared[p]/hsquared[p];
		      }
		  }
		if (h_elem != 0.) h_elem_inv_deriv = -h_elem_deriv/h_elem/h_elem;
		// h_elem_inv_deriv = 0.; /* PRS: NOT SURE WHY THIS IS NOT RIGHT, SO SET TO ZERO */

	      }

	    /* Add mass term */
		mass = 0.;
	    if ( pd->TimeIntegration != STEADY )
	      {
		if ( pd->e[eqn] & T_MASS )
		  {
		    mass  = fv_dot->sh_t;
		    mass *= (- dH_dmesh[b][j]* phi_i * rho * Cp * det_J               * wt
			     -  H            * phi_i * rho * Cp * fv->dsurfdet_dx[b][jk] * wt);
		    mass *= h3;
		    mass *= pd->etm[eqn][(LOG2_MASS)];
		  }
	      }

	    /* Add advection term */
	    advection = 0.;
	    if ( pd->e[eqn] & T_ADVECTION )
	      {

		for ( p=0; p<VIM; p++)
		  {
		    advection += det_J * q[p]                  * d_grad_T_dmesh[p][b][jk];
		    advection += det_J * dq_dmesh[p][b][jk]     * grad_T[p];
		    advection += fv->dsurfdet_dx[b][jk] *q[p] * grad_T[p];
		  }

		advection *= -  rho * Cp * wt * wt_func;
		advection *= h3;
		
		if(supg != 0.)
		  {
		    d_wt_func = 0.;
		    for( p=0; p<dim; p++ )
		      {
			d_wt_func += supg
			  * (h_elem_inv *(q[p]/H)* d_grad_II_phi_i_dmesh[p] [b][jk]
			     +  h_elem_inv_deriv * (q[p]/H) *grad_II_phi_i[p] 
			     -  h_elem_inv * (q[p]/H/H) * grad_II_phi_i[p] * dH_dmesh[b][j]);

		      }
		    for (p=0; p< dim; p++)
		      {
			advection += (-rho * Cp * d_wt_func * h3 * det_J * wt)*q[p] * grad_T[p];
		      }

		  }
		advection *= pd->etm[eqn][(LOG2_ADVECTION)];

	      }
	    
	    /* Add diffusion term */
	    diffusion = 0.0;
	    if ( pd->e[eqn] & T_DIFFUSION ) {
	      for ( p = 0; p < dim; p++) {
		diffusion += H * det_J * d_grad_T_dmesh[p][b][jk] * grad_II_phi_i[p];
		diffusion += H * det_J * grad_T[p] * d_grad_II_phi_i_dmesh[p][b][jk];
		diffusion += H * fv->dsurfdet_dx[b][jk] * grad_T[p] * grad_II_phi_i[p];
		diffusion += dH_dmesh[b][j] * det_J * grad_T[p] * grad_II_phi_i[p];
	      }
	    }
	    diffusion *= k_eff * wt * h3 * pd->etm[eqn][(LOG2_DIFFUSION)];

	    /* Add source term */
	    source = 0.0;
	    if ( pd->e[eqn] & T_SOURCE ) {
	      source += q_tot * fv->dsurfdet_dx[b][jk];
	      source += 0.;   /* Add on your d_q_tot_dx sensitivities here */
	      source *= phi_i;
	    }
	    source *= if_liquid * wt * h3 * pd->etm[eqn][(LOG2_SOURCE)];
	    
	    lec->J[peqn][pvar][i][jk] += mass+ advection + diffusion + source;
	  } // End of loop over j
	} // End of loop over b
      } // End of J_sh_energy_mesh

     /* 
       * J_shell_energy_drs
       */
      var = SOLID_DISPLACEMENT1;
      if ( upd->vp[var] >= 0 && (
			  mp->FSIModel == FSI_REALSOLID_CONTINUUM )) {
	pvar = upd->vp[var];

	/*** Loop over dimensions of mesh displacement ***/
	for ( b = 0; b < dim; b++) {
	  var = SOLID_DISPLACEMENT1 + b;
	  pvar = upd->vp[var];

	  /*** Loop over DOFs (j) ***/
	  /***PRS/SAR should be ei->dof */
	  /***And all d_mesh pieces like dsurfdet_dx etc. should all be jk=dof_map[j] ***/
	  for ( j=0; j < ei->dof[var]; j++) {
	    jk = dof_map[j];
	    
	    /* Load basis functions (j) */
	    phi_j = bf[eqn]->phi[j];
	    for ( ii = 0; ii < VIM; ii++) {
	      grad_phi_j[ii] = bf[eqn]->grad_phi[j][ii];
	      grad_II_phi_j[ii] = 0.0;
	      for ( jj = 0; jj < VIM; jj++) {
		grad_II_phi_j[ii] +=  ( grad_phi_j[jj] * delta(ii,jj) -
					grad_phi_j[jj] * ( fv->snormal[ii] * fv->snormal[jj] ) );
	      }
	    }

	    if(supg!=0.)
	      {
		h_elem_deriv = 0.;
		h_elem_inv_deriv=0.;
		for( p=0; p<dim; p++ )
		  {
		    if(pg_data->hhv[p][b] != 0.)
		      {
			h_elem_deriv -=
			  vcent[p] * vcent[p] *
			  pg_data->dhv_dxnode[p][j] * pg_data->hhv[p][b] * h_elem_inv/4./hsquared[p]/hsquared[p];
		      }
		  }
		if (h_elem != 0.) h_elem_inv_deriv = -h_elem_deriv/h_elem/h_elem;
		// h_elem_inv_deriv = 0.; /* PRS: NOT SURE WHY THIS IS NOT RIGHT, SO SET TO ZERO */

	      }

	    /* Add mass term */
		mass = 0.;
	    if ( pd->TimeIntegration != STEADY )
	      {
		if ( pd->e[eqn] & T_MASS )
		  {
		    mass  = fv_dot->sh_t;
		    mass *= (- dH_drealsolid[b][j]* phi_i * rho * Cp * det_J * wt);
		    mass *= h3;
		    mass *= pd->etm[eqn][(LOG2_MASS)];
		  }
	      }

	    /* Add advection term */
	    advection = 0.;
	    if ( pd->e[eqn] & T_ADVECTION )
	      {

		for ( p=0; p<VIM; p++)
		  {
		    advection += det_J * dq_drs[p][b][jk]     * grad_T[p];
		  }

		advection *= -  rho * Cp * wt * wt_func;
		advection *= h3;
		
		if(supg != 0.)
		  {
		    d_wt_func = 0.;
		    for( p=0; p<dim; p++ )
		      {
			d_wt_func -= supg
			  * (h_elem_inv * (q[p]/H/H) * grad_II_phi_i[p] * dH_drealsolid[b][j]);

		      }
		    for (p=0; p< dim; p++)
		      {
			advection += (-rho * Cp * d_wt_func * h3 * det_J * wt)*q[p] * grad_T[p];
		      }

		  }
		advection *= pd->etm[eqn][(LOG2_ADVECTION)];

	      }
	    
	    /* Add diffusion term */
	    diffusion = 0.0;
	    if ( pd->e[eqn] & T_DIFFUSION ) {
	      for ( p = 0; p < dim; p++) {
		diffusion += dH_drealsolid[b][j] * det_J * grad_T[p] * grad_II_phi_i[p];
	      }
	    }
	    diffusion *= k_eff * wt * h3 * pd->etm[eqn][(LOG2_DIFFUSION)];

	    /* Add source term */
	    source = 0.0;
	    if ( pd->e[eqn] & T_SOURCE ) {
	      source += 0.;  /* Add your d_q_tot_dH terms here */
	      source *= phi_i;
	    }
	    source *= if_liquid * wt * h3 * pd->etm[eqn][(LOG2_SOURCE)];
	    
	    lec->J[peqn][pvar][i][jk] += mass+ advection + diffusion + source;
	  } // End of loop over j
	} // End of loop over b
      } // End of J_sh_energy_rs

     /* 
       * J_shell_energy_dvelocity
       */
      var = VELOCITY1;
      if ( pd->v[var] ) {
	pvar = upd->vp[var];

	/*** Loop over dimensions of velocity ***/
	for ( b = 0; b < dim; b++) {
	  var = VELOCITY1 + b;
	  pvar = upd->vp[var];

	  /*** Loop over DOFs (j) ***/
	  for ( j=0; j < ei->dof[var]; j++) {
	    
	    /* Load basis functions (j) */
	    phi_j = bf[eqn]->phi[j];
	    for ( ii = 0; ii < VIM; ii++) {
	      grad_phi_j[ii] = bf[eqn]->grad_phi[j][ii];
	      grad_II_phi_j[ii] = 0.0;
	      for ( jj = 0; jj < VIM; jj++) {
		grad_II_phi_j[ii] +=  ( grad_phi_j[jj] * delta(ii,jj) -
					grad_phi_j[jj] * ( fv->snormal[ii] * fv->snormal[jj] ) );
	      }
	    }

	    if(supg!=0.)
	      {
		h_elem_deriv = 0.;
		if(hsquared[b] != 0.)
		  {
		    h_elem_deriv =
		      vcent[b] * pg_data->dv_dnode[b][j] * h_elem_inv/4./hsquared[b];
		  }
		if (h_elem != 0.) h_elem_inv_deriv = -h_elem_deriv/h_elem/h_elem;

	      }



	    /* Add advection term */
	    advection = 0.;
	    if ( pd->e[eqn] & T_ADVECTION )
	      {

		if(supg != 0.)
		  {
		    d_wt_func = 0.;
		    for( p=0; p<dim; p++ )
		      {
			d_wt_func += supg* h_elem_inv_deriv * (q[p]/H) *grad_II_phi_i[p];
		      }
			     
		    for( p=0; p<dim; p++)
		      {
			advection += (-rho * Cp * d_wt_func * h3 * det_J * wt)*q[p] * grad_T[p];
		      }
		  }
		advection *= pd->etm[eqn][(LOG2_ADVECTION)];

	      }	    
	    lec->J[peqn][pvar][i][j] += advection;
	  } // End of loop over j
	} // End of loop over b
      } // End of J_sh_energy_dvelocity
      
      /* 
       * J_shell_energy_d_dh
       */
      var = SHELL_DELTAH;
      if ( pd->v[var]) {
	pvar = upd->vp[var];

	/*** Loop over dimensions of mesh displacement ***/
	  var = SHELL_DELTAH;
	  pvar = upd->vp[var];

	  /*** Loop over DOFs (j) ***/
	  for ( j=0; j < ei->dof[var]; j++) {
	    
	    /* Load basis functions (j) */
	    phi_j = bf[eqn]->phi[j];
	    for ( ii = 0; ii < VIM; ii++) {
	      grad_phi_j[ii] = bf[eqn]->grad_phi[j][ii];
	      grad_II_phi_j[ii] = 0.0;
	      for ( jj = 0; jj < VIM; jj++) {
		grad_II_phi_j[ii] +=  ( grad_phi_j[jj] * delta(ii,jj) -
					grad_phi_j[jj] * ( fv->snormal[ii] * fv->snormal[jj] ) );
	      }
	    }

	    /* Add mass term */
	    mass = 0.;
	    if ( pd->TimeIntegration != STEADY )
	      {
		if ( pd->e[eqn] & T_MASS )
		  {
		    mass  = fv_dot->sh_t;
		    mass *= (-phi_j* phi_i * rho * Cp * det_J * wt);
		    mass *= h3;
		    mass *= pd->etm[eqn][(LOG2_MASS)];
		  }
	      }

	    /* Add advection term */
	    advection = 0.;
	    if ( pd->e[eqn] & T_ADVECTION )
	      {
		d_wt_func = 0.;
		if(supg != 0.)
		  {
		    for( p=0; p<dim; p++ )
		      {
			d_wt_func -= supg * h_elem_inv * (q[p]/H/H) * grad_II_phi_i[p] *phi_j;
			d_wt_func += supg * h_elem_inv * (dq_lub_dH[p]/H) * grad_II_phi_i[p] * phi_j;

		      }
		  }

		for ( p=0; p<VIM; p++)
		  {
		    advection += dq_lub_dH[p] * phi_j * grad_T[p] * wt_func;
		    advection += q[p] * grad_T[p] * d_wt_func;
		  }

		advection *= -  rho * Cp * det_J * wt ;
		advection *= h3;
		advection *= pd->etm[eqn][(LOG2_ADVECTION)];
	      }
	    
	    /* Add diffusion term */
	    diffusion = 0.0;
	    if ( pd->e[eqn] & T_DIFFUSION ) {
	      for ( p = 0; p < dim; p++) {

		diffusion += k_eff * phi_j * det_J * grad_T[p] * grad_II_phi_i[p] * dH_U_ddh;
		diffusion += d_k_eff_dh * dH_U_ddh * phi_j * H * det_J* grad_T[p] * grad_II_phi_i[p];
	      }
	    }
	    diffusion *= wt * h3 * pd->etm[eqn][(LOG2_DIFFUSION)];

	    /* Add source term */
	    source = 0.0;
	    if ( pd->e[eqn] & T_SOURCE ) {
	      source += 0.;  /* Add on your dq_tot_dh sensitivities here */
	      source *= phi_i;
	    }
	    source *= if_liquid * det_J * wt * h3 * pd->etm[eqn][(LOG2_SOURCE)];
	    
	    lec->J[peqn][pvar][i][j] += mass+ advection + diffusion + source;
	  } // End of loop over j
      } // End of J_sh_energy_d_dh

    /* 
     * J_shell_energy_d_lubp 
     */
      var = LUBP;
      if ( pd->v[var] ) {
	pvar = upd->vp[var];
	
	/*** Loop over DOFs (j) ***/
	for ( j=0; j<ei->dof[var]; j++) {

	  /* Load basis functions (j) */
	  phi_j = bf[eqn]->phi[j];

	  /* Mass term */
	  mass = 0.;	
	  if ( pd->TimeIntegration != STEADY )
	    {
	      if ( pd->e[eqn] & T_MASS )
		{
		  mass  = fv_dot->sh_t;
		  mass *= - dH_U_dp * phi_j * phi_i * rho * Cp * det_J * wt;
		  mass *= h3;
		  mass *= pd->etm[eqn][(LOG2_MASS)];
		}
	    }

	  /* add advection term */
	  advection = 0.;
	  if ( pd->e[eqn] & T_ADVECTION )
	    {

	      for ( a=0; a<VIM; a++)
		{
		  advection += dqdp[a] * grad_II_phi_j[a]     * grad_T[a];
		  advection += dqdp_1[a] * phi_j     * grad_T[a];
		}

	      advection *= - rho * Cp * det_J * wt * wt_func;
	      advection *= h3;
	      advection *= pd->etm[eqn][(LOG2_ADVECTION)];

	      advection_b = 0.;
	      if(supg != 0.)
		{
		  for ( a=0; a<VIM; a++)
		    {
		      for(p=0; p<VIM; p++)
			{
			  advection_b += q[a] * grad_T[a] * (supg*h_elem_inv*(dqdp[p]*grad_II_phi_j[p]/H)*grad_II_phi_i[p]);
			  advection_b += q[a] * grad_T[a] * (supg*h_elem_inv*(dqdp_1[p]*phi_j/H)*grad_II_phi_i[p]);
			  advection_b -= q[a] * grad_T[a] * (supg*h_elem_inv*(q[p]/H/H * dH_U_dp * phi_j)*grad_II_phi_i[p]);
			}
		    }
		  advection_b *= -rho * Cp * det_J * wt;
		   
		  
		}
	      advection += advection_b;
	    }

	  /* Assemble diffusion term */
	  diffusion = 0.0;
	  if (pd->e[eqn] & T_DIFFUSION) {
	
	    for (p = 0; p < dim; p++) {
	      diffusion += dH_U_dp * phi_j * k_eff * grad_T[p] * grad_II_phi_i[p];
	    }
	
	    diffusion *= det_J * wt * h3 * pd->etm[eqn][(LOG2_DIFFUSION)];
	  } 
	  source = 0.0;
	  if ( pd->e[eqn] & T_SOURCE ) {
	    source += 0.; /* Add on your dq_tot_d_lubp terms here */
	    source *= phi_i;
	  }
	  source *= if_liquid * wt * det_J* h3 * pd->etm[eqn][(LOG2_SOURCE)];
	  
	  lec->J[peqn][pvar][i][j] +=  mass + advection + diffusion + source ;
	}
      }
    } /* end of loop over i */
  } /* end of Assemble_Jacobian */
  
    /* clean-up */
  fv->wt = wt;  /* load_neighbor_var_data screws this up */
  safe_free((void *) n_dof);
  return(status);
} /* end of assemble_shell_energy */

/*****************************************************************************/
/***assemble_film******************************************************/
/*  _______________________________________________________________________  */

/* assemble_film -- assemble terms (Residual & Jacobian) for scalar film profile
 *                  equation
 *
 * in:
 * 	ei -- pointer to Element Indices	structure
 *	pd -- pointer to Problem Description	structure
 *	af -- pointer to Action Flag		structure
 *	bf -- pointer to Basis Function		structure
 *	fv -- pointer to Field Variable		structure
 *  fv_old -- pointer to old Diet Field Variable	structure
 *  fv_dot -- pointer to dot Diet Field Variable	structure
 *	cr -- pointer to Constitutive Relation	structure
 *	md -- pointer to Mesh Derivative	structure
 *	me -- pointer to Material Entity	structure
 *
 * out:
 *	a   -- gets loaded up with proper contribution 
 *	lec -- gets loaded up with local contributions to resid, Jacobian
 * 	r   -- residual RHS vector
 *
 * Created:	Wednesday January 6 2010 tjiptowi@unm.edu
 *
 */
/*ARGSUSED*/
int 
assemble_film( double time,	/* present time value */ 
	       double tt,	/* parameter to vary time integration from 
				   explicit (tt = 1) to implicit (tt = 0)    */
	       double dt,	/* current time step size */  
               double xi[DIM],  /* Local stu coordinates */
               const Exo_DB *exo)

{
  int eqn; 
  int var, peqn, pvar, dim, p;
  int i = -1, ii;
  int j, status;
  int *n_dof = NULL;
  int dof_map[MDE];

  dbl grad_II_P[DIM];  /* Lubrication pressure gradient. */
  dbl grad_II_H[DIM];  /* Film thickness gradient */
  dbl wt;

  dbl P, H, C, H_dot;
  dbl sigma; 
  dbl EvapRate, dEvapRate_dC, dEvapRate_dH;
  dbl veloU[DIM], veloL[DIM];
  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct;  /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;

  dbl mass, diffusion, source;



  /*
   * Galerkin weighting functions for i-th shell residuals 
   * and some of their derivatives...
   */
  dbl phi_i;
  dbl grad_phi_i[DIM];
  dbl grad_II_phi_i[DIM];

  /*
   * Interpolation functions for variables and some of their derivatives.
   */
  dbl phi_j;
  dbl grad_phi_j[DIM];
  dbl grad_II_phi_j[DIM];

  dbl h3;			/* Volume element (scale factors). */
  dbl det_J;

  /*   static char yo[] = "assemble_film";*/

  status = 0;


  /*
   * Prepare geometry
   */
  n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
  lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);


  /* Unpack variables from structures for local convenience... */
  dim = pd->Num_Dim;
  P = fv->sh_fp;  /* Lubrication pressure at the thin film */
  H = fv->sh_fh;  /* Film thickness */

  if (pd->v[SHELL_PARTC])
   {
     C = fv->sh_pc;
   } 
  else
   {
     C = 0.0;
   }


  wt = fv->wt;       /* Gauss weight */
  h3 = fv->h3;	     /* Differential volume element, = 1 when CARTESIAN. */
  

  det_J = fv->sdet;

  /*
   * Material property constants, etc. Any variations at this Gauss point 
   * are calculated with user-defined subroutine options.
   * 
   */

  if (pd->v[SHELL_PARTC])
   {
     viscosity(gn, NULL, d_mu);
   }

   sigma = mp->surface_tension;

    
  velocity_function_model(veloU, veloL, time, dt);

  EvapRate = film_evaporation_model(C, &dEvapRate_dC, H, &dEvapRate_dH);
 

   if (pd->TimeIntegration != STEADY) {
    H_dot = fv_dot->sh_fh;
  } else {
    H_dot = 0.0;
  }
  

/* Calculate lubrication pressure and thickness gradients */

  Inn(fv->grad_sh_fh, grad_II_H);
  Inn(fv->grad_sh_fp, grad_II_P);


/* Calculate flow rate and sensitivities */

  calculate_lub_q_v (R_LUBP, time, dt, xi, exo );
  


  /*
   * Residuals_________________________________________________________________
   */

  if (af->Assemble_Residual)
    {

/* ************** ASSEMBLE RESIDUAL OF LUBRICATION PRESSURE ********* */

      eqn = R_SHELL_FILMP;
      peqn = upd->ep[eqn];

      for ( i=0; i<ei->dof[eqn]; i++)
	{
	  phi_i = bf[eqn]->phi[i];
	  for (p = 0; p < dim; p++) 
	     {
	      grad_phi_i[p] = bf[eqn]->grad_phi[i][p];
              grad_II_phi_i[p] = 0.0;
	     }

          Inn(grad_phi_i, grad_II_phi_i);



/* Assemble mass term */

	  mass = 0.0;
	  if ( pd->TimeIntegration != STEADY )
            {
              if ( pd->e[eqn] & T_MASS )
                {
                  mass  = phi_i * H_dot;
                  mass *= det_J * wt;
                  mass *= h3;
                  mass *= pd->etm[eqn][(LOG2_MASS)];
                }
            }


/* Assemble diffusion term */

	  diffusion = 0.0;
	  if (pd->e[eqn] & T_DIFFUSION) 
	  {
	    for (p = 0; p < dim; p++) 
	      {
		diffusion += - LubAux->q[p]*grad_II_phi_i[p];
	      }

	    diffusion *= det_J * wt * h3 * pd->etm[eqn][(LOG2_DIFFUSION)];
	  }

/* Assemble source term */

	  source = 0.0;
	  if (pd->e[eqn] & T_SOURCE)
	    {
	      source  = phi_i * EvapRate;     
	      source *= det_J * wt * h3 * pd->etm[eqn][(LOG2_SOURCE)];
	    }

/* Combine them all */

	  lec->R[peqn][i] += mass + diffusion + source;  
	  
	} /* end of loop over i */


/* ************* ASEMBLE RESIDUAL OF FILM THICKNESS ************ */

      eqn = R_SHELL_FILMH;
      peqn = upd->ep[eqn];

      for ( i=0; i<ei->dof[eqn]; i++)
	{
	  phi_i = bf[eqn]->phi[i];

	  for (p = 0; p < dim; p++) 
	     {
	      grad_phi_i[p] = bf[eqn]->grad_phi[i][p];
              grad_II_phi_i[p] = 0.0;
	     }


          Inn(grad_phi_i, grad_II_phi_i);



/* Assemble diffusion term */

	  diffusion = 0.0;
	  if (pd->e[eqn] & T_DIFFUSION) 
	  {
	    for (p = 0; p < dim; p++) 
	      {
		diffusion += - sigma*grad_II_phi_i[p]*grad_II_H[p];
	      }

	    diffusion *= det_J * wt * h3 * pd->etm[eqn][(LOG2_DIFFUSION)];
	  }

/* Assemble source term */

	  source = 0.0;
	  if (pd->e[eqn] & T_SOURCE)
	    {
	      source  = P;
              source *= phi_i;
	      source *= det_J * wt * h3 * pd->etm[eqn][(LOG2_SOURCE)];
	    }

/* Combine them all */

	  lec->R[peqn][i] +=  diffusion + source;  
	  
	} /* end of loop over i */


    } /* end of Assemble_Residual */



  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian)
    {

/* ************* ASEMBLE JACOBIAN OF LUBRICATION PRESSURE ************ */

      eqn   = R_SHELL_FILMP;
      peqn = upd->ep[eqn];

      for ( i=0; i<ei->dof[eqn]; i++)
	{

	  phi_i = bf[eqn]->phi[i];

	  /*
	   * Set up some preliminaries that are needed for the (a,i)
	   * equation for bunches of (b,j) column variables...
	   */

	  for ( p=0; p<dim; p++)
	    {
	      grad_phi_i[p] = bf[eqn]->grad_phi[i][p];
              grad_II_phi_i[p] = 0.0;
	    }

          Inn(grad_phi_i, grad_II_phi_i);


	  /* SENSITIVITY W.R.T. LUBRICATION PRESSURE */

	  var = SHELL_FILMP;

	  if ( pd->v[var] )
	    {
	      pvar = upd->vp[var];
	      for ( j=0; j<ei->dof[var]; j++)
		{
	           phi_j = bf[var]->phi[j];
                  
		  for ( p=0; p<dim; p++)
		    {
		      grad_phi_j[p] = bf[var]->grad_phi[j][p];
                      grad_II_phi_j[p] = 0.0;
		    }

                  Inn(grad_phi_j, grad_II_phi_j);


		  diffusion = 0.0;
		  if (pd->e[eqn] & T_DIFFUSION)
		    {
		      for (ii = 0; ii < VIM; ii++) 
			{
                          diffusion += - LubAux->dq_dp1[ii][j] * grad_II_phi_i[ii] * grad_II_phi_j[ii]; 
			}

		      diffusion *= det_J * wt;
		      diffusion *= h3;
		      diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
		    }


		  lec->J[peqn][pvar][i][j] += diffusion;
		}
	    }

	  /* SENSITIVITY W.R.T. FILM THICKNESS */


	  var = SHELL_FILMH;

	  if (pd->v[var])
	    {
	      pvar = upd->vp[var];
	      for ( j=0; j<ei->dof[var]; j++)
		{

	          phi_j = bf[var]->phi[j];

		  for ( p=0; p<dim; p++)
		    {
		      grad_phi_j[p] = bf[var]->grad_phi[j][p];
                      grad_II_phi_j[p] = 0.0;
		    }

                  Inn(grad_phi_j, grad_II_phi_j);

                  mass = 0.;
                  if ( pd->TimeIntegration != STEADY )
                    {
                      if ( pd->e[eqn] & T_MASS )
                        {
                          mass += (1. + 2. * tt) * phi_j / dt;
                          mass *= phi_i * det_J * wt;
                          mass *= h3;
                          mass *= pd->etm[eqn][(LOG2_MASS)];
                        }
                    }

		  diffusion = 0.;
		  if (pd->e[eqn] & T_DIFFUSION)
		    {

		     for (ii = 0; ii < VIM; ii++) 
			{
	                  diffusion += - LubAux->dq_dh1[ii][j] * grad_II_phi_j[ii] * grad_II_phi_i[ii];
	                  diffusion += - LubAux->dq_dh2[ii][j] * phi_j * grad_II_phi_i[ii];
                        }
		    
		      diffusion *= det_J * wt;
		      diffusion *= h3;
		      diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];

		    }


                  source = 0.;
                  if ( pd->e[eqn] & T_SOURCE )
                    {
                       source += dEvapRate_dH * phi_j;
                       source *= phi_i * det_J * wt;
                       source *= h3;
                       source *= pd->etm[eqn][(LOG2_SOURCE)];
                    }


		  lec->J[peqn][pvar][i][j] += mass + diffusion + source;
		}

            }

	  /* SENSITIVITY W.R.T. PARTICLES CONCENTRATION */


	  var = SHELL_PARTC;

	  if (pd->v[var])
	    {
	      pvar = upd->vp[var];
	      for ( j=0; j<ei->dof[var]; j++)
		{
	          phi_j = bf[var]->phi[j];

		  diffusion = 0.;
		  if (pd->e[eqn] & T_DIFFUSION)
		    {

		     for (ii = 0; ii < VIM; ii++) 
			{
	                  diffusion += - LubAux->dq_dc[ii][j] * phi_j * grad_II_phi_i[ii];
                        }
		    
		      diffusion *= det_J * wt;
		      diffusion *= h3;
		      diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];

		    }

                  source = 0.;
                  if ( pd->e[eqn] & T_SOURCE )
                    {
                       source += dEvapRate_dC * phi_j;
                       source *= phi_i * det_J * wt;
                       source *= h3;
                       source *= pd->etm[eqn][(LOG2_SOURCE)];
                    }

		  lec->J[peqn][pvar][i][j] += diffusion + source;   

		}

	    }	 


	  var = R_MESH1;

	  if ( pd->v[var] )
	    {
	      EH(-1,"See comment in code where this line is printed. You need to add mesh sensitivities to sh_fp");
	      /*
	       * Mesh sensitivities will be to grad_phi and det_j.   these are available now for shells. Start
	       * looking in load_bf_grad and shell_determinant_and_normal routines.
	       */
	    }

        }

/* ************* ASEMBLE JACOBIAN OF FILM THICKNESS ************ */

      eqn   = R_SHELL_FILMH;
      peqn = upd->ep[eqn];

      for ( i=0; i<ei->dof[eqn]; i++)
	{

	  phi_i = bf[eqn]->phi[i];

	  /*
	   * Set up some preliminaries that are needed for the (a,i)
	   * equation for bunches of (b,j) column variables...
	   */

	  for ( p=0; p<dim; p++)
	    {
	      grad_phi_i[p] = bf[eqn]->grad_phi[i][p];
              grad_II_phi_i[p] = 0.0;
	    }

          Inn(grad_phi_i, grad_II_phi_i);



	  /* SENSITIVITY W.R.T. LUBRICATION PRESSURE */

	  var = SHELL_FILMP;

	  if ( pd->v[var] )
	    {
	      pvar = upd->vp[var];
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  phi_j = bf[var]->phi[j];

		  source = 0.0;
		  if (pd->e[eqn] & T_SOURCE)
		    {
                      source += phi_i * phi_j; 
		      source *= det_J * wt;
		      source *= h3;
		      source *= pd->etm[eqn][(LOG2_DIFFUSION)];
		    }


		  lec->J[peqn][pvar][i][j] += source;
		}
	    }

	  /* SENSITIVITY W.R.T. FILM THICKNESS */


	  var = SHELL_FILMH;

	  if (pd->v[var])
	    {
	      pvar = upd->vp[var];
	      for ( j=0; j<ei->dof[var]; j++)
		{

	          phi_j = bf[var]->phi[j];

		  for ( p=0; p<dim; p++)
		    {
		      grad_phi_j[p] = bf[var]->grad_phi[j][p];
                      grad_II_phi_j[p] = 0.0;
		    }

                  Inn(grad_phi_j, grad_II_phi_j);


		  diffusion = 0.;
		  if (pd->e[eqn] & T_DIFFUSION)
		    {
		      for (ii = 0; ii < VIM; ii++) 
                        {
                             diffusion += - sigma * grad_II_phi_i[ii] * grad_II_phi_j[ii];
			}

		      diffusion *= det_J * wt;
		      diffusion *= h3;
		      diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];

		    }

		  lec->J[peqn][pvar][i][j] += diffusion;
		}
	    }

	      
	  /* SENSITIVITY W.R.T. PARTICLES CONCENTRATION */
          /* To be pursued when coupling with surface tension is finished */	 


	  var = R_MESH1;

	  if ( pd->v[var] )
	    {
	      EH(-1,"See comment in code where this line is printed. You need to add mesh sensitivities to sh_fh");
	      /*
	       * Mesh sensitivities will be to grad_phi and det_j.   these are available now for shells. Start
	       * looking in load_bf_grad and shell_determinant_and_normal routines.
	       */
	    }

	}


    } /* end of Assemble_Jacobian */

  safe_free((void *) n_dof);

  return(status);
} /* end of assemble_film */




/*****************************************************************************/
/***assemble_film_particles******************************************************/
/*  _______________________________________________________________________  */

/* assemble_film_particles -- assemble terms (Residual & Jacobian) for scalar z-averaged
 *                            particles conservation  equation
 *
 * in:
 * 	ei -- pointer to Element Indices	structure
 *	pd -- pointer to Problem Description	structure
 *	af -- pointer to Action Flag		structure
 *	bf -- pointer to Basis Function		structure
 *	fv -- pointer to Field Variable		structure
 *  fv_old -- pointer to old Diet Field Variable	structure
 *  fv_dot -- pointer to dot Diet Field Variable	structure
 *	cr -- pointer to Constitutive Relation	structure
 *	md -- pointer to Mesh Derivative	structure
 *	me -- pointer to Material Entity	structure
 *
 * out:
 *	a   -- gets loaded up with proper contribution 
 *	lec -- gets loaded up with local contributions to resid, Jacobian
 * 	r   -- residual RHS vector
 *
 * Created:	Wednesday January 6 2010 tjiptowi@unm.edu
 *
 */
/*ARGSUSED*/
int 
assemble_film_particles(  double time,	/* present time value */ 
	                  double tt,	/* parameter to vary time integration from 
				           explicit (tt = 1) to implicit (tt = 0)    */
	                  double dt,	/* current time step size */
                          double xi[DIM],  /* Local stu coordinates */
                          const PG_DATA *pg_data, /*Upwinding stuff */
                          const Exo_DB *exo  )

{
  int eqn; 
  int var, peqn, pvar, dim, p;
  int i = -1, ii;
  int j, status;
  int *n_dof = NULL;
  int dof_map[MDE];

  dbl grad_II_C[DIM];  /* Particles concentration gradient */
  dbl wt;

  dbl H = 0, C, C_dot; 
  dbl H_U, dH_U_dtime, H_L, dH_L_dtime, dH_U_dp, dH_U_ddh;
  dbl dH_U_dX[DIM],dH_L_dX[DIM];
  dbl q_old[DIM], q[DIM], v[DIM];
  dbl dq_dp1[DIM][MDE], dq_dp2[DIM][MDE], dq_dh1[DIM][MDE], dq_dh2[DIM][MDE], dq_dc[DIM][MDE];
  dbl dv_dp1[DIM][MDE], dv_dp2[DIM][MDE], dv_dh1[DIM][MDE], dv_dh2[DIM][MDE], dv_dc[DIM][MDE];
  dbl mu, dmu_dc;
  dbl EvapRate, dEvapRate_dC, dEvapRate_dH;
  dbl diff_coeff, ddiff_dmu, ddiff_dc;
  dbl veloU[DIM], veloL[DIM];
  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct;  /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;
  dbl mass, advection, diffusion, source;

  const double *hsquared = pg_data->hsquared ;
  const double *vcent = pg_data->v_avg; /* Average element velocity, which is the
                                           centroid velocity for Q2 and the average
                                           of the vertices for Q1. It comes from
                                           the routine "element_velocity." */

  /*
   * Petrov-Galerkin weighting functions for i-th residuals
   * and some of their derivatives...
   */

  dbl wt_func;

/* SUPG variables */
  dbl h_elem=0, h_elem_inv=0;
  dbl supg = 0.;

  /*
   * Galerkin weighting functions for i-th shell residuals 
   * and some of their derivatives...
   */
  dbl phi_i;
  dbl grad_phi_i[DIM];
  dbl grad_II_phi_i[DIM];

  /*
   * Interpolation functions for variables and some of their derivatives.
   */
  dbl phi_j;
  dbl grad_phi_j[DIM];
  dbl grad_II_phi_j[DIM];

  dbl h3;			/* Volume element (scale factors). */
  dbl det_J;

  /*   static char yo[] = "assemble_film_particles";*/

  status = 0;

   if( mp->Spwt_funcModel == GALERKIN)
    {
      supg = 0.;
    }
  else if( mp->Spwt_funcModel == SUPG)
    {
      if( !pd->e[R_MOMENTUM1])
        EH(-1, " must have momentum equation velocity field for shell_particle_films upwinding");
      supg = mp->Spwt_func;
    }

  /*
   * Prepare geometry
   */
  n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
  lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);



  /* Unpack variables from structures for local convenience... */
  dim = pd->Num_Dim;

  if ( pd->v[SHELL_FILMH] )
    { 
     H = fv->sh_fh;  /* Film thickness */
    }
  else if ( pd->v[LUBP] )
    {
     H = height_function_model(&H_U, &dH_U_dtime, &H_L, &dH_L_dtime, dH_U_dX, dH_L_dX, 
                               &dH_U_dp, &dH_U_ddh, time, dt);
    }
  else
    {
     EH(-1,"In assemble_film_particles: Can't find an appropriate lubrication height");
    }

  C = fv->sh_pc;  /* Particles concentration */

  wt = fv->wt;       /* Gauss weight */
  h3 = fv->h3;	     /* Differential volume element, = 1 when CARTESIAN. */
  

  shell_determinant_and_normal(ei->ielem, ei->iconnect_ptr, ei->num_local_nodes, 
			       ei->ielem_dim, 1);
  det_J = fv->sdet;

  /*
   * Material property constants, etc. Any variations at this Gauss point 
   * are calculated with user-defined subroutine options.
   * 
   */

   mu = viscosity(gn, NULL, d_mu);
   diff_coeff = diffusion_coefficient_model(mu, &ddiff_dmu);

   dmu_dc = mp->d_viscosity[SHELL_PARTC];
   ddiff_dc = ddiff_dmu * dmu_dc;

    
  velocity_function_model(veloU, veloL, time, dt);
  
  if ( pd->v[SHELL_FILMH] )
    { 
     EvapRate = film_evaporation_model(C, &dEvapRate_dC, H, &dEvapRate_dH);
    }
  else
    {
     EvapRate = 0.0;
     dEvapRate_dC = 0.0;
     dEvapRate_dH = 0.0;
    }

   if (pd->TimeIntegration != STEADY) {
    C_dot = fv_dot->sh_pc;
  } else {
    C_dot = 0.0;
  }
  

/* Calculate particles concentration  gradients */

   Inn(fv->grad_sh_pc, grad_II_C);


/* Calculate flow rate and sensitivities */

   calculate_lub_q_v (R_LUBP, time, dt, xi, exo );
   calculate_lub_q_v_old ( tran->time_value_old, tran->delta_t_old, xi, exo );

   for (p=0; p<dim; p++)
      {
       q[p] = LubAux->q[p];
       v[p] = LubAux->v_avg[p];
       q_old[p] = LubAux_old->q[p];
       for (j=0; j<ei->dof[SHELL_PARTC]; j++)
          {
           dq_dp1[p][j] = LubAux->dq_dp1[p][j];
           dq_dp2[p][j] = LubAux->dq_dp2[p][j];
           dq_dh1[p][j] = LubAux->dq_dh1[p][j];
           dq_dh2[p][j] = LubAux->dq_dh2[p][j];
           dq_dc[p][j]  = LubAux->dq_dc[p][j]; 
           dv_dp1[p][j] = LubAux->dv_avg_dp1[p][j];
           dv_dp2[p][j] = LubAux->dv_avg_dp2[p][j];
           dv_dh1[p][j] = LubAux->dv_avg_dh1[p][j];
           dv_dh2[p][j] = LubAux->dv_avg_dh2[p][j];
           dv_dc[p][j]  = LubAux->dv_avg_dc[p][j];
          }
      }


/* Calculate SUPG related variables */

   if (supg!=0.)
     {
      h_elem = 0.;
      for ( p=0; p<dim; p++)
         {
          if(hsquared[p] != 0.) h_elem += vcent[p] * vcent[p] / hsquared[p];
         }
      h_elem = sqrt(h_elem)/2.;
      if (h_elem == 0.)
        {
          h_elem_inv=0.;
        }
      else
        {
         h_elem_inv=1./h_elem;
        }

     }
  

  /*
   * Residuals_________________________________________________________________
   */

  if (af->Assemble_Residual)
    {


/* ************** ASSEMBLE RESIDUAL OF PARTICLES CONCENTRATION ********* */

      eqn = R_SHELL_PARTC;
      peqn = upd->ep[eqn];

      for ( i=0; i<ei->dof[eqn]; i++)
	{
	  phi_i = bf[eqn]->phi[i];
	  for (p = 0; p < dim; p++) 
	     {
	      grad_phi_i[p] = bf[eqn]->grad_phi[i][p];
              grad_II_phi_i[p] = 0.0;
	     }

          Inn(grad_phi_i, grad_II_phi_i);


/* Assemble mass term */

	  mass = 0.0;
	  if ( pd->TimeIntegration != STEADY )
            {
              if ( pd->e[eqn] & T_MASS )
                {
                  mass  = phi_i * C_dot * H;
                  mass *= det_J * wt;
                  mass *= h3;
                  mass *= pd->etm[eqn][(LOG2_MASS)];
                }
            }


/* Assemble advection term */

          wt_func = bf[eqn]->phi[i];

#if 0

          /* add Petrov-Galerkin terms as necessary */
          if (supg!=0.)
            {
              for(p=0; p<dim; p++)
                {
                  wt_func += supg * h_elem_inv * LubAux->v_avg[p] * grad_II_phi_i[p] ;
                }
            }

	  advection = 0.0;
	  if (pd->e[eqn] & T_ADVECTION) 
	  {
	    for (p = 0; p < dim; p++) 
	      {
		advection +=  wt_func * LubAux->q[p] * grad_II_C[p];
	      }

#endif

	  advection = 0.0;
	  if (pd->e[eqn] & T_ADVECTION) 
	  {
	    for (p = 0; p < dim; p++) 
	      {
		advection +=   (1.5 * q[p] - 0.5 * q_old[p]) * grad_II_C[p] * wt_func
	                     + q[p] * grad_II_C[p] * v[p] * grad_II_phi_i[p] * 0.5 * dt;
              }
	    advection *= det_J * wt * h3 * pd->etm[eqn][(LOG2_ADVECTION)];
	  }


/* Assemble diffusion term */

	  diffusion = 0.0;
	  if (pd->e[eqn] & T_DIFFUSION) 
	  {
	    for (p = 0; p < dim; p++) 
	      {
                diffusion +=   diff_coeff * H * grad_II_C[p] * grad_II_phi_i[p];
	      }

	    diffusion *= det_J * wt * h3 * pd->etm[eqn][(LOG2_DIFFUSION)];
	  }

/* Assemble source term */

	  source = 0.0;
	  if (pd->e[eqn] & T_SOURCE) 
	  {
            source +=  - phi_i * C * EvapRate;
	    source *=    det_J * wt * h3 * pd->etm[eqn][(LOG2_SOURCE)];
	  }


/* Combine them all */

	  lec->R[peqn][i] += mass + advection + diffusion + source;  
	  
	} /* end of loop over i */

    } /* End of assembly residuals */


  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian)
    {


/* ************* ASSEMBLE JACOBIAN OF PARTICLES CONCENTRATION ************ */

      eqn   = R_SHELL_PARTC;
      peqn = upd->ep[eqn];

      for ( i=0; i<ei->dof[eqn]; i++)
	{

	  phi_i = bf[eqn]->phi[i];

	  /*
	   * Set up some preliminaries that are needed for the (a,i)
	   * equation for bunches of (b,j) column variables...
	   */

	  for ( p=0; p<dim; p++)
	    {
	      grad_phi_i[p] = bf[eqn]->grad_phi[i][p];
              grad_II_phi_i[p] = 0.0;
            }

           Inn(grad_phi_i, grad_II_phi_i);

           wt_func = bf[eqn]->phi[i];

          /* add Petrov-Galerkin terms as necessary */
          if (supg!=0.)
            {
              for(p=0; p<dim; p++)
                {
                  wt_func += supg * h_elem_inv * LubAux->v_avg[p] * grad_II_phi_i[p] ;
                }
            }


	  /* SENSITIVITY W.R.T. LUBRICATION PRESSURE */

	  var = LUBP;

	  if ( pd->v[var] )
	    {
	      pvar = upd->vp[var];
	      for ( j=0; j<ei->dof[var]; j++)
		{
	           phi_j = bf[var]->phi[j];
                  
		  for ( p=0; p<dim; p++)
		    {
		      grad_phi_j[p] = bf[var]->grad_phi[j][p];
                      grad_II_phi_j[p] = 0.0;
		    }

                  Inn(grad_phi_j, grad_II_phi_j);

		  advection = 0.0;

		  if (pd->e[eqn] & T_ADVECTION)
		    {

		      for (ii = 0; ii < VIM; ii++) 
			{
                          advection +=  1.5 * dq_dp1[ii][j] * grad_II_phi_j[ii] * grad_II_C[ii]* wt_func
                                      + 1.5 * dq_dp2[ii][j] * phi_j * grad_II_C[ii]* wt_func
                                      + dq_dp1[ii][j] * grad_II_phi_j[ii] * grad_II_C[ii] * v[ii] * grad_II_phi_i[ii] * 0.5 * dt
                                      + dq_dp2[ii][j] * phi_j * grad_II_C[ii] * v[ii] * grad_II_phi_i[ii] * 0.5 * dt
                                      + q[ii] * grad_II_C[ii] * dv_dp1[ii][j] * grad_II_phi_j[ii] * grad_II_phi_i[ii] * 0.5 * dt
                                      + q[ii] * grad_II_C[ii] * dv_dp2[ii][j] * phi_j * grad_II_phi_i[ii] * 0.5 * dt; 
			}

		      advection *= det_J * wt;
		      advection *= h3;
		      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
		    }

		  lec->J[peqn][pvar][i][j] += advection;
		}
	    }


	  /* SENSITIVITY W.R.T. FILM PRESSURE */

	  var = SHELL_FILMP;

	  if ( pd->v[var] )
	    {
	      pvar = upd->vp[var];
	      for ( j=0; j<ei->dof[var]; j++)
		{
	           phi_j = bf[var]->phi[j];
                  
		  for ( p=0; p<dim; p++)
		    {
		      grad_phi_j[p] = bf[var]->grad_phi[j][p];
                      grad_II_phi_j[p] = 0.0;
		    }

                  Inn(grad_phi_j, grad_II_phi_j);

		  advection = 0.0;
#if 0
		  if (pd->e[eqn] & T_ADVECTION)
		    {

		      for (ii = 0; ii < VIM; ii++) 
			{
                          advection += LubAux->dq_dp1[ii][j] * grad_II_phi_j[ii] * grad_II_C[ii] * wt_func; 
			}

                      if (supg != 0.)
                        {
                         for ( ii = 0; ii < VIM; ii++)
                            {
                             for (p = 0; p < VIM; p++)
                                {
                                 advection += LubAux->q[ii] * grad_II_C[ii] * 
                                              (supg * h_elem_inv * LubAux->dv_avg_dp1[p][j] * grad_II_phi_j[p] * grad_II_phi_i[p]);
                                }
                            }                   
                        } 

		      advection *= det_J * wt;
		      advection *= h3;
		      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
		    }
#endif
		  if (pd->e[eqn] & T_ADVECTION)
		    {

		      for (ii = 0; ii < VIM; ii++) 
			{
                          advection +=  1.5 * dq_dp1[ii][j] * grad_II_phi_j[ii] * grad_II_C[ii]* wt_func
                                      + dq_dp1[ii][j] * grad_II_phi_j[ii] * grad_II_C[ii] * v[ii] * grad_II_phi_i[ii] * 0.5 * dt
                                      + q[ii] * grad_II_C[ii] * dv_dp1[ii][j] * grad_II_phi_j[ii] * grad_II_phi_i[ii] * 0.5 * dt; 
			}

		      advection *= det_J * wt;
		      advection *= h3;
		      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
		    }


		  lec->J[peqn][pvar][i][j] += advection;
		}
	    }

	  /* SENSITIVITY W.R.T. FILM THICKNESS */


	  var = SHELL_FILMH;

	  if (pd->v[var])
	    {
	      pvar = upd->vp[var];
	      for ( j=0; j<ei->dof[var]; j++)
		{

	          phi_j = bf[var]->phi[j];

		  for ( p=0; p<dim; p++)
		    {
		      grad_phi_j[p] = bf[var]->grad_phi[j][p];
                      grad_II_phi_j[p] = 0.0;
		    }

                  Inn(grad_phi_j, grad_II_phi_j);


                  mass = 0.;
                  if ( pd->TimeIntegration != STEADY )
                    {
                      if ( pd->e[eqn] & T_MASS )
                        {
                          mass += phi_j * C_dot;
                          mass *= phi_i * det_J * wt;
                          mass *= h3;
                          mass *= pd->etm[eqn][(LOG2_MASS)];
                        }
                    }

		  advection = 0.;
#if 0
		  if (pd->e[eqn] & T_ADVECTION)
		    {
                     for (ii = 0; ii < VIM; ii++)
                        {
	                  advection +=  LubAux->dq_dh1[ii][j] * grad_II_phi_j[ii] * grad_II_C[ii] * wt_func;
	                  advection +=  LubAux->dq_dh2[ii][j] * phi_j * grad_II_C[ii] * wt_func;
                        }

                      if (supg != 0.)
                        {
                         for ( ii = 0; ii < VIM; ii++)
                            {
                             for (p = 0; p < VIM; p++)
                                {
                                 advection += LubAux->q[ii] * grad_II_C[ii] * 
                                              (supg * h_elem_inv * 
                                               (  LubAux->dv_avg_dh1[p][j] * grad_II_phi_j[p] * grad_II_phi_i[p]
                                                + LubAux->dv_avg_dh2[p][j] * phi_j * grad_II_phi_i[p] ) );
                                }
                            }                   
                        } 
		    
		      advection *= det_J * wt;
		      advection *= h3;
		      advection *= pd->etm[eqn][(LOG2_ADVECTION)];

		    }
#endif

		  if (pd->e[eqn] & T_ADVECTION)
		    {
		      for (ii = 0; ii < VIM; ii++) 
			{
                          advection +=  1.5 * dq_dh1[ii][j] * grad_II_phi_j[ii] * grad_II_C[ii] * wt_func
                                      + 1.5 * dq_dh2[ii][j] * phi_j * grad_II_C[ii] * wt_func
                                      + dq_dh1[ii][j] * grad_II_phi_j[ii] * grad_II_C[ii] * v[ii] * grad_II_phi_i[ii] * 0.5 * dt 
                                      + dq_dh2[ii][j] * phi_j * grad_II_C[ii] * v[ii] * grad_II_phi_i[ii] * 0.5 * dt
                                      + q[ii] * grad_II_C[ii] * dv_dh1[ii][j] * grad_II_phi_j[ii] * grad_II_phi_i[ii] * 0.5 * dt
                                      + q[ii] * grad_II_C[ii] * dv_dh2[ii][j] * phi_j * grad_II_phi_i[ii] * 0.5 * dt;
			}

		      advection *= det_J * wt;
		      advection *= h3;
		      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
                    } 


		  diffusion = 0.;
		  if (pd->e[eqn] & T_DIFFUSION)
		    {

		     for (ii = 0; ii < VIM; ii++) 
			{
                          diffusion +=   grad_II_phi_i[ii] * phi_j * diff_coeff * grad_II_C[ii];
                        }
		    
		      diffusion *= det_J * wt;
		      diffusion *= h3;
		      diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];

		    }

		  source = 0.;
		  if (pd->e[eqn] & T_SOURCE)
		    {
                      source += - C * dEvapRate_dH * phi_j * phi_i;
		    
		      source *= det_J * wt;
		      source *= h3;
		      source *= pd->etm[eqn][(LOG2_SOURCE)];

		    }

		  lec->J[peqn][pvar][i][j] += mass + advection + diffusion + source;
		}
	    }

	      
	  /* SENSITIVITY W.R.T. PARTICLES CONCENTRATION */
	 

	  var = SHELL_PARTC;

	  if (pd->v[var])
	    {
	      pvar = upd->vp[var];
	      for ( j=0; j<ei->dof[var]; j++)
		{

	          phi_j = bf[var]->phi[j];


		  for ( p=0; p<dim; p++)
		    {
		      grad_phi_j[p] = bf[var]->grad_phi[j][p];
                      grad_II_phi_j[p] = 0.0;
		    }

                  Inn(grad_phi_j, grad_II_phi_j);

                  mass = 0.;
                  if ( pd->TimeIntegration != STEADY )
                    {
                      if ( pd->e[eqn] & T_MASS )
                        {
                          mass += H * (1. + 2. * tt) * phi_j / dt;
                          mass *= phi_i * det_J * wt;
                          mass *= h3;
                          mass *= pd->etm[eqn][(LOG2_MASS)];
                        }
                    }

		  advection = 0.;
		  if (pd->e[eqn] & T_ADVECTION)
		    {
#if 0
		     for (ii = 0; ii < VIM; ii++) 
			{
	                  advection += grad_II_phi_j[ii] * LubAux->q[ii] * wt_func;
                          advection += grad_II_C[ii] * LubAux->dq_dc[ii][j] * phi_j * wt_func;
                        }
#endif		    
		     for (ii = 0; ii < VIM; ii++) 
			{
	                  advection +=  1.5 * dq_dc[ii][j] * grad_II_C[ii] * wt_func
                                      + (1.5 * q[ii] - 0.5 * q_old[ii]) * grad_II_phi_j[ii] * wt_func
                                      + dq_dc[ii][j] * grad_II_C[ii] * v[ii] * grad_II_phi_i[ii] * 0.5 * dt
                                      + q[ii] * grad_II_phi_j[ii] * v[ii] * grad_II_phi_i[ii] * 0.5 * dt
                                      + q[ii] * grad_II_C[ii] * dv_dc[ii][j] * grad_II_phi_i[ii] * 0.5 * dt;
                        }
		      advection *= det_J * wt;
		      advection *= h3;
		      advection *= pd->etm[eqn][(LOG2_ADVECTION)];

		    }


		  diffusion = 0.;
		  if (pd->e[eqn] & T_DIFFUSION)
		    {

		     for (ii = 0; ii < VIM; ii++) 
			{
                          diffusion +=  grad_II_phi_i[ii] * H * ddiff_dc * phi_j * grad_II_C[ii]; 
                          diffusion +=  grad_II_phi_i[ii] * H * diff_coeff * grad_II_phi_j[ii];
                        }
		    
		      diffusion *= det_J * wt;
		      diffusion *= h3;
		      diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];

		    }

		  source = 0.;
		  if (pd->e[eqn] & T_SOURCE)
		    {
	              source += - phi_j * EvapRate * phi_i;
                      source += - C * dEvapRate_dC * phi_j * phi_i;
		    
		      source *= det_J * wt;
		      source *= h3;
		      source *= pd->etm[eqn][(LOG2_SOURCE)];

		    }

		  lec->J[peqn][pvar][i][j] += mass + advection + diffusion + source;
		}
	    }


          /* SENSITIVITY W.R.T. VELOCITY */


          var = VELOCITY1;
#if 0
          if (pd->v[var])
            {

              /*** Loop over dimensions of velocity ***/
              for ( b = 0; b < dim; b++) 
                 {
                  var = VELOCITY1 + b;
                  pvar = upd->vp[var];

                  /*** Loop over DOFs (j) ***/
                  for ( j=0; j < ei->dof[var]; j++) 
                     {
                      phi_j = bf[var]->phi[j];


                      if (supg!=0.)
                        {
                         h_elem_deriv = 0.;

                         if (hsquared[b] != 0.)
                           {
                            h_elem_deriv = vcent[b] * pg_data->dv_dnode[b][j] * h_elem_inv/ (4. * hsquared[b]);
                           }
                         if (h_elem != 0.) 
                           {
                            h_elem_inv_deriv = -h_elem_deriv/(h_elem * h_elem);
                           }
                        }

                      advection = 0.;
                      if ( pd->e[eqn] & T_ADVECTION )
                        {

                         if (supg != 0.)
                           {
                            d_wt_func = 0.;
                            for ( p=0; p<dim; p++ )
                               {
                                d_wt_func += supg* h_elem_inv_deriv * LubAux->v_avg[p] * grad_II_phi_i[p];
                               }
                            for ( ii = 0; ii < dim; ii++ )
                               {
                                 advection += LubAux->q[ii] * grad_II_C[ii] * d_wt_func;
                               }
                           }
                         advection *= h3 * det_J * wt;
                         advection *= pd->etm[eqn][(LOG2_ADVECTION)];

                        }
                      lec->J[peqn][pvar][i][j] += advection;            
                     }
                 }
            }
#endif
	  var = R_MESH1;

	  if ( pd->v[var] )
	    {
	      EH(-1,"See comment in code where this line is printed. You need to add mesh sensitivities to sh_pc");
	      /*
	       * Mesh sensitivities will be to grad_phi and det_j.   these are available now for shells. Start
	       * looking in load_bf_grad and shell_determinant_and_normal routines.
	       */
	    }

        }


    } /* end of Assemble_Jacobian */

   safe_free((void *) n_dof);

  return(status);
} /* end of assemble_film_particles */


/******************************************************************************/


/*
 * dPdz_function()
 * Actually calculates the pressure gradient and sensitivities based on values
 * passed to it.  The logic for determining regularization, in/out of drop, etc.
 * is not included here, just the basic pressure gradient physics.
 * Scott A Roberts
 */
void dPdz_function(
		   dbl tt,
		   dbl dt,
		   dbl S,
		   dbl Plub,
		   dbl nbar,
		   dbl H,
		   dbl P0,
		   dbl sigma,
		   dbl theta,
		   dbl R,
		   dbl dS,
		   dbl *dPdz,
		   dbl *dPdz_S,
		   dbl *dPdz_Plub,
		   dbl *dPdz_nbar,
		   dbl *dPdz_S_Plub,
		   dbl *dPdz_S_nbar
		   )
{ 

  /* Parameters */
  dbl Pgas, Pgas_S, Pgas_nbar, Pgas_S_nbar;
  dbl Pcap, Pcap_S;

  /* Load properties */
  dbl Patm = mp->PorousShellPatm;
  dbl Pref = mp->PorousShellPref;

  /* Calculate wetting velocity */
  dbl V   = H * fv_dot->sh_sat_closed;
  dbl V_S = H * (1+2*tt)/ dt;

  /* Calculate dynamic contact angle */
  dbl cosT, cosT_V, a1, a2;
  dynamic_contact_angle_model(&cosT, &a1, V, &cosT_V, &a2);

  /* Calculate capillary pressure */
  Pcap        = 2*sigma*cosT/R;
  Pcap_S      = 2*sigma*cosT_V*V_S/R;

  /* Calculate minimum lubrication pressure */
  Pgas        = (P0)*nbar/(1-dS)-(Pref-Patm);
  
  /* Calculate gas pressure */
  Pgas        = P0*nbar/(1-S) + (Pref-Patm);
  Pgas_S      = P0*nbar/pow(1-S,2);
  Pgas_nbar   = P0/(1-S);
  Pgas_S_nbar = P0/pow(1-S,2);
  //Pgas        = P0*nbar/(1-S/100) + (Pref-Patm);
  //Pgas_S      = P0*nbar/pow(1-S/100,2)/100;
  //Pgas_nbar   = P0/(1-S/100);
  //Pgas_S_nbar = P0/pow(1-S/100,2);

  /* Calculate pressure drop */
  *dPdz        = ( Pgas - Plub - Pcap ) / (S*H);
  *dPdz_S      = ( Pgas_S - Pcap_S ) / (S*H) - ( Pgas - Plub - Pcap ) / (S*S*H);
  *dPdz_Plub   = -1 / (S*H);
  *dPdz_nbar   = Pgas_nbar / (S*H);
  *dPdz_S_Plub = 1 / (S*S*H);
  *dPdz_S_nbar = -Pgas_nbar / (S*S*H) + Pgas_S_nbar / (S*H);
  
  return;
}

/*
 * dPdz_calc()
 * This function calculates the pressure gradient for a shell of structured
 * closed pores.  It internally loads all of the necessary values and only
 * returns the pressure gradients and sensitivities.  It contains logic
 * to properly regularize and to deal with level set drop fields.
 *
 * Recent additions allow
 * evaluation at either the nodes or at the gauss points.  To load at the gauss
 * point, pass a negative value for ni.  To load at the nodes, just pass 
 * the node index.
 */
void dPdz_calc(
	       dbl tt,
	       dbl dt,
	       dbl *dPdz,
	       dbl *dPdz_S,
	       dbl *dPdz_Plub,
	       dbl *dPdz_nbar,
	       dbl *dPdz_F,
	       dbl *dPdz_curv
	       )
{
  int i;

  /* Load static material properties */
  dbl sigma = mp->surface_tension;                     // Surface tension
  dbl p0 = mp->PorousShellClosedP0;                    // Initial gas pressure
  dbl theta = mp->dcaU*M_PIE/180.0;                    // Contact angle

  /* Load external field parameters */
  dbl phi = porous_shell_closed_porosity_model();      // Porosity
  dbl r = porous_shell_closed_radius_model();          // Pore radius
  dbl H = porous_shell_closed_height_model();          // Pore height
  
  /* Load field variables */
  dbl S = fv->sh_sat_closed;                           // Saturation
  dbl plub = fv->lubp;                                 // Lubrication pressure
  dbl nbar = 1;                                        // Gas volume multiplier
  if ( pd->e[R_SHELL_SAT_GASN] ) nbar = fv->sh_sat_gasn;

  /* Lubrication height from model */
  dbl Hlub, H_U, dH_U_dtime, H_L, dH_L_dtime, dH_U_dp, dH_U_ddh;
  dbl dH_U_dX[DIM], dH_L_dX[DIM];
  Hlub = height_function_model(&H_U, &dH_U_dtime, &H_L, &dH_L_dtime, dH_U_dX, dH_L_dX, &dH_U_dp, &dH_U_ddh, 0.0, 0.0); 

  /* Level set curvature */
  int b;
  dbl dcaU, dcaL, slopeU, slopeL, curv = 0.0;
  dbl slopeU_F[MDE] = {0.0};
  dbl slopeL_F[MDE] = {0.0};
  dbl curv_F[MDE] = {0.0};
  dbl oldlsls = 0.0;
  if ( pd->v[FILL] ) {
    oldlsls = ls->Length_Scale;
    load_lsi( ls->Length_Scale );
    load_lsi_derivs();
    dcaU = mp->dcaU*M_PIE/180.0;
    dcaL = mp->dcaL*M_PIE/180.0;
    slopeU = slopeL = 0.;
    for ( b = 0; b < DIM; b++) {
      slopeU += dH_U_dX[b]*lsi->normal[b];
      slopeL += dH_L_dX[b]*lsi->normal[b];
      for ( i = 0; i < ei->dof[FILL]; i++) {
	slopeU_F[i] += dH_U_dX[b] * lsi->d_normal_dF[b][i];
	slopeL_F[i] += dH_L_dX[b] * lsi->d_normal_dF[b][i];
      }
    }
    curv= (cos(M_PIE-dcaU-atan(slopeU)) + cos(M_PIE-dcaL-atan(-slopeL)))/Hlub;
    for ( i = 0; i < ei->dof[FILL]; i++) {
      curv_F[i] -= sin(-M_PIE+dcaU+atan(slopeU))*slopeU_F[i]/(Hlub*(1+pow(slopeU,2)));
      curv_F[i] += sin(-M_PIE+dcaL-atan(slopeL))*slopeL_F[i]/(Hlub*(1+pow(slopeL,2)));
    }
    if ( pd->e[SHELL_LUB_CURV] ) {
      curv += fv->sh_l_curv;
    }
  }

  /* Correct pressure field for surface tension */
  dbl plub_F[MDE] = {0.0};
  if (pd->v[FILL]) {
    plub = plub + sigma*curv*lsi->Hn;
    for ( i = 0; i < ei->dof[FILL]; i++) {
      //plub_F[i] = sigma*curv*lsi->d_Hn_dF[i]*bf[FILL]->phi[i] + sigma*curv_F[i]*lsi->Hn;
      plub_F[i] = sigma*curv*lsi->d_Hn_dF[i] + sigma*curv_F[i]*lsi->Hn;
    }
  }

  /* Determine the mode of calculation */
  int pores_on = phi > 1.0e-4;                         // Are there pores above me?

  /* Calculate boundary values */
  dbl dS = 1e-2;
  dbl Smin = dS;
  dbl Smax = 1-dS;

  /* Calculate pressure gradient */
  dbl dPdz_S_Plub, dPdz_S_nbar;
  if ( S < -Smin ) {
    dPdz_function(tt,dt,Smin,plub,nbar,H,p0,sigma,theta,r,dS,dPdz,dPdz_S,dPdz_Plub,dPdz_nbar,&dPdz_S_Plub,&dPdz_S_nbar);
    *dPdz      = (*dPdz      + *dPdz_S*(S-Smin));
    *dPdz_S    = (*dPdz_S);
    *dPdz_Plub = (*dPdz_Plub + dPdz_S_Plub*(S-Smin));
    *dPdz_nbar = (*dPdz_nbar + dPdz_S_nbar*(S-Smin));
  } else if ( !pores_on ) {
    *dPdz      = 0.0;
    *dPdz_S    = 0.0;
    *dPdz_Plub = 0.0;
    *dPdz_nbar = 0.0;
  } else if ( S < Smin ) {
    dPdz_function(tt,dt,Smin,plub,nbar,H,p0,sigma,theta,r,dS,dPdz,dPdz_S,dPdz_Plub,dPdz_nbar,&dPdz_S_Plub,&dPdz_S_nbar);
    *dPdz      = (*dPdz      + *dPdz_S*(S-Smin));
    *dPdz_S    = (*dPdz_S);
    *dPdz_Plub = (*dPdz_Plub + dPdz_S_Plub*(S-Smin));
    *dPdz_nbar = (*dPdz_nbar + dPdz_S_nbar*(S-Smin));
  } else if ( S > Smax ) {
    dPdz_function(tt,dt,Smax,plub,nbar,H,p0,sigma,theta,r,dS,dPdz,dPdz_S,dPdz_Plub,dPdz_nbar,&dPdz_S_Plub,&dPdz_S_nbar);
    *dPdz      = (*dPdz      + *dPdz_S*(S-Smax));
    *dPdz_S    = (*dPdz_S);
    *dPdz_Plub = (*dPdz_Plub + dPdz_S_Plub*(S-Smax));
    *dPdz_nbar = (*dPdz_nbar + dPdz_S_nbar*(S-Smax));
  } else {
    dPdz_function(tt,dt,S,plub,nbar,H,p0,sigma,theta,r,dS,dPdz,dPdz_S,dPdz_Plub,dPdz_nbar,&dPdz_S_Plub,&dPdz_S_nbar);
  }

  /* Calculate sensitivity for curvature */
  if ( pd->v[FILL] ) *dPdz_curv = *dPdz_Plub * sigma*lsi->Hn;

  /* Calculate weighting for level set fields */
  dbl HsideM = 1.0, d_HsideM_dF[MDE] = {0.0};
  if ( pd->v[FILL] ) {
    HsideM = 1 - 2*lsi->Hn;
    for ( i = 0; i < ei->dof[R_FILL]; i++) {
      d_HsideM_dF[i] = -2*lsi->d_Hn_dF[i];
    }
  }

  /* Bang with heaviside function to not use it outside the drop */
  dbl bang;
  if ( pd->v[FILL] ) {
    bang = HsideM;
    if ( bang < 0 ) {
      *dPdz      = 0.0;
      *dPdz_S    = 0.0;
      *dPdz_Plub = 0.0;
      *dPdz_nbar = 0.0;
      *dPdz_curv = 0.0;
      for (i=0;i<ei->dof[R_FILL];i++) dPdz_F[i] = 0.0;
    } else {
      *dPdz      *= bang;
      *dPdz_S    *= bang;
      *dPdz_Plub *= bang;
      *dPdz_nbar *= bang;
      *dPdz_curv *= bang;
      for ( i = 0; i < ei->dof[R_FILL]; i++) {
  	dPdz_F[i] = *dPdz_Plub*plub_F[i]*bang + *dPdz*d_HsideM_dF[i];
      }
    }
  }

  /* Reset level set field */
  if ( pd->v[FILL] ) load_lsi( oldlsls );  

  return;
}
	       
/*****************************************************************************/
/***assemble_porous_shell_closed**********************************************/
/*  _______________________________________________________________________  */

/* assemble_porous_shell_closed -- Assemble terms (Residual & Jacobian) for a
 *                                 structured porous shell equation for closed
 *                                 pores
 *
 * Created:	Thursday February 25 04:30:00 MST 2010 sarober@sandia.gov 
 * Revised:     Friday March 26, 2010 sarober@sandia.gov
 *                                    Renamed for closed pores
 *
 */
/*ARGSUSED*/ 
int
assemble_porous_shell_closed(
		      dbl tt,                  // Time integration form
		      dbl dt,                  // Time step size
		      dbl xi[DIM],
		      const Exo_DB *exo
) {
  
  /* --- Initialization -----------------------------------------------------*/    

  // Variable definitions
  int eqn, peqn, var, pvar;                    // Internal integers
  int i, j, jj, b;                             // Counter variables
  dbl phi_i, phi_j;                            // Basis functions
  dbl mass, divergence, source;                // Residual terms

  // Initialize output status function
  int status = 0;

  // Bail out of function if there is nothing to do
  eqn = R_SHELL_SAT_CLOSED;
  if (!pd->e[eqn]) return(status);

  /* Setup lubrication */
  int *n_dof = NULL;
  int dof_map[MDE];
  dbl wt_old = fv->wt;
  n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
  lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

  // Unpack FEM variables from global structures
  dbl wt = fv->wt;                             // Gauss point weight
  dbl h3 = fv->h3;                             // Differential volume element
  dbl det_J = fv->sdet;                        // Jacobian of transformation

  // Load material properties
  dbl mu = mp->viscosity;                      // Constant viscosity

  // Load external file parameters
  dbl phi = porous_shell_closed_porosity_model();  // Porosity
  dbl r = porous_shell_closed_radius_model();      // Pore radius
  dbl H = porous_shell_closed_height_model();      // Pore height

  // Prepare switches
  int mass_on = pd->e[eqn] & T_MASS;
  int divergence_on = pd->e[eqn] & T_DIVERGENCE;
  int transient_on = pd->TimeIntegration != STEADY;

  // FEM weights
  dbl dA = det_J * wt * h3;
  dbl etm_mass = pd->etm[eqn][(LOG2_MASS)];
  dbl etm_divergence = pd->etm[eqn][(LOG2_DIVERGENCE)];

  /* --- Calculate equation components ---------------------------------------*/    

  // Calculate permeability
  dbl k = phi*r*r/8;

  // Now, calculate velocities at the GP, for use as the lubrication source
  dbl dPdz, dPdz_S, dPdz_plub, dPdz_nbar, dPdz_F[MDE], dPdz_curv;
  dPdz_calc( tt, dt, &dPdz, &dPdz_S, &dPdz_plub, &dPdz_nbar, dPdz_F, &dPdz_curv);

  // Calculate velocities
  dbl vz_F[MDE];
  dbl vz      = -k/mu*dPdz;
  dbl vz_S    = -k/mu*dPdz_S;
  dbl vz_plub = -k/mu*dPdz_plub;
  dbl vz_nbar = -k/mu*dPdz_nbar;
  for ( i = 0; i < ei->dof[R_FILL]; i++) vz_F[i] = -k/mu*dPdz_F[i];
  dbl vz_curv = -k/mu*dPdz_curv;

  /* --- Assemble residuals --------------------------------------------------*/

  int mytest[MDE];
  for ( i = 0; i < ei->dof[eqn]; i++) {
    if ( pd->e[R_FILL] ) {
      mytest[i] = fv->F < 0;
    } else {
      mytest[i] = 1;
    }
  }

  // Assemble residual contribution to this equation
  eqn = R_SHELL_SAT_CLOSED;
  if (af->Assemble_Residual) {
    peqn = upd->ep[eqn];
    
    // Loop over DOF (i)
    for ( i = 0; i < ei->dof[eqn]; i++) {
      
      // Load basis functions
      phi_i = bf[eqn]->phi[i];
      
      // Assemble mass term
      mass = 0.0;
      if ( mass_on && transient_on ) {
	//mass += fv_dot->sh_sat_closed * phi_i;
	mass += *esp_dot->sh_sat_closed[i] * phi_i;
      }
      mass *= dA * etm_mass;

      // Assemble divergence term
      divergence = 0.0;
      if (divergence_on && mytest[i]) {
	divergence += -vz/(phi*H) * phi_i;
      }
      divergence *= dA * etm_divergence;
      
      // Assemble full residual
      lec->R[peqn][i] += mass + divergence;
      
    }  // End of loop over DOF (i)

  } // End of residual assembly of R_SHELL_SAT_CLOSED
  

  // Assemble residual contribution to the lubrication equation
  eqn = R_LUBP;
  if ( af->Assemble_Residual & pd->e[eqn] & transient_on ) {
    peqn = upd->ep[eqn];
    
    // Loop over DOF (i)
    for ( i = 0; i < ei->dof[eqn]; i++ ) { 
      
      // Load basis functions
      phi_i = bf[eqn]->phi[i];
      
      // Assemble source term
      source = 0.0;
      if ( mytest[i] ) {
	source -= vz * phi_i;
      }
      source *= dA * pd->etm[eqn][(LOG2_SOURCE)];
      
      // Add to residual vector
      lec->R[peqn][i] += source;
      
    } // End lf loop over DOF (i)
    
  } // End of residual assembly of R_LUBP 
  
  
  /* --- Assemble Jacobian --------------------------------------------------*/
  eqn = R_SHELL_SAT_CLOSED;
  if (af->Assemble_Jacobian) {
    peqn = upd->ep[eqn];
    
    // Loop over DOF (i)
    for ( i = 0; i < ei->dof[eqn]; i++) {
      
      // Load basis functions
      phi_i = bf[eqn]->phi[i];
      
      // Assemble sensitivities for SHELL_SAT_CLOSED
      var = SHELL_SAT_CLOSED;
      if (pd->v[var]) {
  	pvar = upd->vp[var];

  	// Loop over DOF (j)
  	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  phi_j = bf[eqn]->phi[j];
	  
	  // Assemble mass terms
	  mass = 0.0;
	  if ( mass_on && transient_on ) { 
	    //mass += phi_i * phi_j * (1+2*tt)/dt;
	    if ( i == j ) mass += phi_i * (1+2*tt)/dt;
	  }
	  mass *= dA * etm_mass;
	  
	  // Assemble divergence term
	  divergence = 0.0;
	  if (divergence_on && mytest[i]) {
	    divergence += -vz_S/(phi*H) * phi_i * phi_j;
	  }
	  divergence *= dA * etm_divergence;
	  
	  // Assemble full Jacobian
	  lec->J[peqn][pvar][i][j] += mass + divergence;
	  
	} // End of loop over DOF (j)
	
      } // End of SH_SAT_CLOSED sensitivities
      
      // Assemble sensitivities for LUBP
      var = LUBP;
      if (pd->v[var]) {
  	pvar = upd->vp[var];

  	// Loop over DOF (j)
  	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  phi_j = bf[eqn]->phi[j];
	  
	  // Assemble mass terms
	  mass = 0.0;
	  
	  // Assemble divergence term
	  divergence = 0.0;
	  if (divergence_on && mytest[i]) {
	    divergence += -vz_plub/(phi*H) * phi_i * phi_j;
	  }
	  divergence *= dA * etm_divergence;
	  
	  // Assemble full Jacobian
	  lec->J[peqn][pvar][i][j] += mass + divergence;
	  
	} // End of loop over DOF (j)
	
      } // End of LUBP sensitivities
      
      // Assemble sensitivities for F
      var = FILL;
      if (pd->v[var]) {
  	pvar = upd->vp[var];

  	// Loop over DOF (j)
  	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  phi_j = bf[var]->phi[j];
	  
	  // Assemble mass terms
	  mass = 0.0;
	  
	  // Assemble divergence term
	  divergence = 0.0;
	  if (divergence_on && mytest[i]) {
	    divergence += -vz_F[j]/(phi*H) * phi_i;
	  }
	  divergence *= dA * etm_divergence;
	  
	  // Assemble full Jacobian
	  lec->J[peqn][pvar][i][j] += mass + divergence;
	  
	} // End of loop over DOF (j)
	
      } // End of F sensitivities
      
      // Assemble sensitivities for SHELL_LUB_CURV
      var = SHELL_LUB_CURV;
      if (pd->v[var]) {
  	pvar = upd->vp[var];

  	// Loop over DOF (j)
  	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  phi_j = bf[eqn]->phi[j];
	  
	  // Assemble mass terms
	  mass = 0.0;
	  
	  // Assemble divergence term
	  divergence = 0.0;
	  if (divergence_on && mytest[i]) {
	    divergence += -vz_curv/(phi*H) * phi_i * phi_j;
	  }
	  divergence *= dA * etm_divergence;
	  
	  // Assemble full Jacobian
	  lec->J[peqn][pvar][i][j] += mass + divergence;
	  
	} // End of loop over DOF (j)
	
      } // End of SHELL_LUB_CURV sensitivities
      
      // Assemble sensitivities for SHELL_SAT_GASN
      var = SHELL_SAT_GASN;
      if (pd->v[var]) {
  	pvar = upd->vp[var];

  	// Loop over DOF (j)
  	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  phi_j = bf[eqn]->phi[j];
	  
	  // Assemble mass terms
	  mass = 0.0;
	  
	  // Assemble divergence term
	  divergence = 0.0;
	  if (divergence_on && mytest[i]) {
	    divergence += -vz_nbar/(phi*H) * phi_i * phi_j;
	  }
	  divergence *= dA * etm_divergence;
	  
	  // Assemble full Jacobian
	  lec->J[peqn][pvar][i][j] += mass + divergence;
	  
	} // End of loop over DOF (j)
	
      } // End of SHELL_SAT_GASN sensitivities

      // Assemble sensitivities for DMX
      var = MESH_DISPLACEMENT1;
      if (pd->v[var] && (
			 mp->FSIModel == FSI_MESH_CONTINUUM ||
			 mp->FSIModel == FSI_MESH_UNDEF
			 )) {
	pvar = upd->vp[var];

	// Loop over dimensions of mesh displacement
	for ( b = 0; b < DIM; b++) {
	  var = MESH_DISPLACEMENT1 + b;
	  pvar = upd->vp[var];
	  
	  // Loop over DOF (j)
	  for ( j = 0; j < ei->dof[var]; j++) {
	    jj = dof_map[j];
	    
	    // Assemble mass term
	    mass = 0.0;
	    if ( mass_on && transient_on ) {   
	      mass += fv_dot->sh_sat_closed * phi_i;
	    }
	    mass *= pd->etm[eqn][(LOG2_MASS)];

	    // Assemble divergence term
	    divergence = 0.0;
	    if (divergence_on && mytest[i]) {
	      divergence += -vz/(phi*H) * phi_i;
	    }
	    divergence *= pd->etm[eqn][(LOG2_DIVERGENCE)];
	    
	    // Assemble Jacobian terms
	    lec->J[peqn][pvar][i][jj] += (mass + divergence) * wt * h3 * fv->dsurfdet_dx[b][jj];
	    
	  } // End of loop over DOF (j)

	} // End of loop over mesh dimensions

      } // End of sensitivity for DMX
                 
    } // End of loop over DOF (i)

  } // End of Jacobian assembly of R_SHELL_SAT_CLOSED

  // Assemble Jacobian entries for R_LUBP
  eqn = R_LUBP;
  if ( af->Assemble_Jacobian & pd->e[eqn] & transient_on ) {
    peqn = upd->ep[eqn];

    // Loop over DOF (i)
    for ( i = 0; i < ei->dof[eqn]; i++) {

      // Load basis functions
      phi_i = bf[eqn]->phi[i];

      // Assemble sensitivities for SHELL_SAT_CLOSED
      var = SHELL_SAT_CLOSED;
      if (pd->v[var]) {
	pvar = upd->vp[var];

	// Loop over DOF (j)
	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  phi_j = bf[var]->phi[j];

	  // Assemble source term
	  source  = 0.0;
	  if ( mytest[i] ) {
	    source -= vz_S * phi_i * phi_j;
	  }
	  source *= dA * pd->etm[eqn][(LOG2_SOURCE)];

	  // Assemble Jacobian terms
	  lec->J[peqn][pvar][i][j] += source;

	} // End of loop over DOF (j)

      } // End of sensitivity for SHELL_SAT_CLOSED

      // Assemble sensitivities for LUBP
      var = LUBP;
      if (pd->v[var]) {
	pvar = upd->vp[var];

	// Loop over DOF (j)
	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  phi_j = bf[eqn]->phi[j];

	  // Assemble source term
	  source  = 0.0;
	  if ( mytest[i] ) {
	    source -= vz_plub * phi_i * phi_j;
	  }
	  source *= dA * pd->etm[eqn][(LOG2_SOURCE)];

	  // Assemble Jacobian terms
	  lec->J[peqn][pvar][i][j] += source;

	} // End of loop over DOF (j)

      } // End of sensitivity for LUBP

      // Assemble sensitivities for F
      var = FILL;
      if (pd->v[var]) {
	pvar = upd->vp[var];

	// Loop over DOF (j)
	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  phi_j = bf[eqn]->phi[j];

	  // Assemble source term
	  source  = 0.0;
	  if ( mytest[i] ) {
	    source -= vz_F[j] * phi_i;
	  }
	  source *= dA * pd->etm[eqn][(LOG2_SOURCE)];

	  // Assemble Jacobian terms
	  lec->J[peqn][pvar][i][j] += source;

	} // End of loop over DOF (j)

      } // End of sensitivity for F

      // Assemble sensitivities for SHELL_LUB_CURV
      var = SHELL_LUB_CURV;
      if (pd->v[var]) {
	pvar = upd->vp[var];

	// Loop over DOF (j)
	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  phi_j = bf[eqn]->phi[j];

	  // Assemble source term
	  source  = 0.0;
	  if ( mytest[i] ) {
	    source -= vz_curv * phi_i * phi_j;
	  }
	  source *= dA * pd->etm[eqn][(LOG2_SOURCE)];

	  // Assemble Jacobian terms
	  lec->J[peqn][pvar][i][j] += source;

	} // End of loop over DOF (j)

      } // End of sensitivity for SHELL_LUB_CURV

      // Assemble sensitivities for SHELL_SAT_GASN
      var = SHELL_SAT_GASN;
      if (pd->v[var]) {
	pvar = upd->vp[var];

	// Loop over DOF (j)
	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  phi_j = bf[eqn]->phi[j];

	  // Assemble source term
	  source  = 0.0;
	  if ( mytest[i] ) {
	    source -= vz_nbar * phi_i * phi_j;
	  }
	  source *= dA * pd->etm[eqn][(LOG2_SOURCE)];

	  // Assemble Jacobian terms
	  lec->J[peqn][pvar][i][j] += source;

	} // End of loop over DOF (j)

      } // End of sensitivity for SHELL_SAT_GASN

      // Assemble sensitivities for DMX
      var = MESH_DISPLACEMENT1;
      if (pd->v[var] && (
			 mp->FSIModel == FSI_MESH_CONTINUUM ||
			 mp->FSIModel == FSI_MESH_UNDEF
			 )) {
	pvar = upd->vp[var];

	// Loop over dimensions of mesh displacement
	for ( b = 0; b < DIM; b++) {
	  var = MESH_DISPLACEMENT1 + b;
	  pvar = upd->vp[var];
	  
	  // Loop over DOF (j)
	  for ( j = 0; j < ei->dof[var]; j++) {
	    jj = dof_map[j];
	    
	    // Load basis functions
	    phi_j = bf[eqn]->phi[j];
	    
	    // Assemble source term
	    source  = 0.0;
	    if ( mytest[i] ) {
	      source -= vz * wt * h3 * fv->dsurfdet_dx[b][jj] * phi_i;
	    }
	    source *= pd->etm[eqn][(LOG2_SOURCE)];
	    
	    // Assemble Jacobian terms
	    lec->J[peqn][pvar][i][jj] += source;
	    
	  } // End of loop over DOF (j)

	} // End of loop over mesh dimensions

      } // End of sensitivity for DMX

    } // End of loop over DOF (i)
    
  } // End of Jacobian assembly of R_LUBP 
  
  // Cleanup 
  fv->wt = wt_old;
  safe_free((void *) n_dof);

  // Finish function
  return(status);

} // End of assemble_porous_shell_closed



	       
/*****************************************************************************/
/***assemble_porous_shell_gasn************************************************/
/*  _______________________________________________________________________  */

/* assemble_porous_shell_gasn -- Assemble terms for a gas compression
 *                               equation for the closed porous shells
 *
 * Created:	December 15, 2010 sarober@sandia.gov 
 *
 */
/*ARGSUSED*/ 
int
assemble_porous_shell_gasn(
			   dbl tt,                  // Time integration form
			   dbl dt,                  // Time step size
			   dbl xi[DIM],             // Coordinates
			   const Exo_DB *exo        // ExodusII database
) {

  /* --- Initialization ------------------------------------------------------*/ 

  // Variable definitions
  int eqn, peqn, var, pvar;                       // Internal integers
  int i, j, b, jk;                             // Counters
  dbl mass, source;                               // Residual / Jacobian terms
  int status = 0;                                 // Function status
  dbl phi_i;                                      // Basis functions (i)
  dbl phi_j;                                      // Basis functions (j)

  // Make sure that you're supposed to be here, bail out otherwise
  eqn = R_SHELL_SAT_GASN;
  if ( !pd->e[eqn] || !pd->e[R_SHELL_SAT_CLOSED] )
    EH(-1, "Woah, you should not be in this function, assemble_porous_shell_gasn().");

  /* Setup lubrication */
  int *n_dof = NULL;
  int dof_map[MDE];
  dbl wt_old = fv->wt;
  n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
  lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

  // Unpack FEM variables from global structures
  dbl wt = fv->wt;
  dbl h3 = fv->h3;
  dbl det_J = fv->sdet;

  // Prepare switches
  int mass_on = pd->e[eqn] & T_MASS;
  int source_on = pd->e[eqn] & T_SOURCE;

  // Load variables
  dbl GASN = fv->sh_sat_gasn;
  dbl SAT  = fv->sh_sat_closed; 

  // Load material parameters
  dbl H    = porous_shell_closed_height_model();
  dbl P0   = mp->PorousShellClosedP0;
  dbl Patm = mp->PorousShellPatm;
  dbl kH   = mp->PorousShellHenry;
  dbl RT   = mp->PorousShellRT;
  dbl D    = mp->PorousShellDiffusivity;

  // Prepare heaviside function for multiphase flow
  dbl HsideM = 1.0, d_HsideM_dF[DIM] = {0.0};
  if ( pd->v[FILL] ) {
    load_lsi( ls->Length_Scale );
    HsideM = 1 - lsi->Hn;
    for ( i = 0; i < ei->dof[R_FILL]; i++) d_HsideM_dF[i] = -lsi->d_Hn_dF[i];
  }


  /* --- Calculate primary components ----------------------------------------*/ 

  // Calculate equation constants
  dbl Z   = RT*D/(H*H*kH);
  dbl Pbr = Patm/P0;

  // Calculate flux
  dbl flux, flux_N, flux_S, flux_F[DIM] = {0.0};
  dbl si, fi, fi_S, fi_N, fi_SN;
  dbl dS = 1e-2;
  if ( (SAT<dS) || (SAT>(1-dS)) ) {
    if ( SAT < dS   ) si = dS;
    if ( SAT > 1-dS ) si = 1-dS;
    fi     = -Z / si * ( GASN/(1-si) - Pbr );
    fi_S   =  Z / pow(si,2) * ( GASN/(1-si) - Pbr ) - Z / si * GASN/pow(1-si,2);
    fi_SN  =  Z / pow(si,2) * ( 1/(1-si) ) - Z / si * 1/pow(1-si,2);
    fi_N   = -Z / si / (1-si);
    flux   = fi   + (SAT-si)*fi_S;
    flux_S =                 fi_S;
    flux_N = fi_N + (SAT-si)*fi_SN;
  } else {
    flux   = -Z / SAT * ( GASN/(1-SAT) - Pbr );
    flux_N = -Z / SAT / (1-SAT);
    flux_S =  Z / pow(SAT,2) * ( GASN/(1-SAT) - Pbr ) - Z / SAT * GASN/pow(1-SAT,2);
  }

  // Hit with heaviside function for multiphase flow
  if ( pd->v[FILL] ) {
    flux   *= HsideM;
    flux_N *= HsideM;
    flux_S *= HsideM;
    for ( i = 0; i < ei->dof[R_FILL]; i++) flux_F[i] = flux * d_HsideM_dF[i];
  }

  /* --- Assemble residuals --------------------------------------------------*/

  // Assemble residual contribution to this equation
  eqn = R_SHELL_SAT_GASN;
  if (af->Assemble_Residual) {
    peqn = upd->ep[eqn];
    
    // Loop over DOF (i)
    for ( i = 0; i < ei->dof[eqn]; i++) {         
      
      // Load basis functions
      phi_i = bf[eqn]->phi[i];

      // Assemble mass terms
      mass = 0.0;
      if ( mass_on ) {
	mass += fv_dot->sh_sat_gasn * phi_i;
      }
      mass *= wt * h3 * det_J * pd->etm[eqn][(LOG2_MASS)];

      // Assemble source terms
      source = 0.0;
      if ( source_on ) {
	source -= flux * phi_i;
      }
      source *= wt * h3 * det_J * pd->etm[eqn][(LOG2_SOURCE)];
      
      // Assemble full residual
      lec->R[peqn][i] += mass + source;
      
    }  // End of loop over DOF (i)

  } // End of residual assembly of R_SHELL_SAT_GASN



  /* --- Assemble Jacobians --------------------------------------------------*/

  // Assemble jacobian contribution to this equation
  eqn = R_SHELL_SAT_GASN;
  if (af->Assemble_Jacobian) {
    peqn = upd->ep[eqn];
    
    // Loop over DOF (i)
    for ( i = 0; i < ei->dof[eqn]; i++) {         
      
      // Load basis functions
      phi_i = bf[eqn]->phi[i];
      
      // Assemble sensitivities for SHELL_SAT_GASN
      var = SHELL_SAT_GASN;
      if (pd->v[var]) {
  	pvar = upd->vp[var];

  	// Loop over DOF (j)
  	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  phi_j = bf[eqn]->phi[j];

	  // Assemble mass terms
	  mass = 0.0;
	  if ( mass_on ) {
	    mass += (1.0 + 2.0*tt)/dt * phi_i * phi_j;
	  }
	  mass *= wt * h3 * det_J * pd->etm[eqn][(LOG2_MASS)];
	  
	  // Assemble source terms
	  source = 0.0;
	  if ( source_on ) {
	    source -= flux_N * phi_i * phi_j;
	  }
	  source *= wt * h3 * det_J * pd->etm[eqn][(LOG2_SOURCE)];
	  
	  // Assemble full residual
	  lec->J[peqn][pvar][i][j] += mass + source;

	} // End of loop over DOF (j)
	
      } // End of SH_SAT_GASN sensitivities
      
      // Assemble sensitivities for SHELL_SAT_CLOSED
      var = SHELL_SAT_CLOSED;
      if (pd->v[var]) {
  	pvar = upd->vp[var];

  	// Loop over DOF (j)
  	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  phi_j = bf[eqn]->phi[j];

	  // Assemble mass terms
	  mass = 0.0;
	  
	  // Assemble source terms
	  source = 0.0;
	  if ( source_on ) {
	    source -= flux_S * phi_i * phi_j;
	  }
	  source *= wt * h3 * det_J * pd->etm[eqn][(LOG2_SOURCE)];
	  
	  // Assemble full residual
	  lec->J[peqn][pvar][i][j] += mass + source;

	} // End of loop over DOF (j)
	
      } // End of SH_SAT_CLOSED sensitivities
      
      // Assemble sensitivities for FILL
      var = FILL;
      if (pd->v[var]) {
  	pvar = upd->vp[var];

  	// Loop over DOF (j)
  	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  phi_j = bf[eqn]->phi[j];

	  // Assemble mass terms
	  mass = 0.0;
	  
	  // Assemble source terms
	  source = 0.0;
	  if ( source_on ) {
	    source -= flux_F[j] * phi_i * phi_j;
	  }
	  source *= wt * h3 * det_J * pd->etm[eqn][(LOG2_SOURCE)];
	  
	  // Assemble full residual
	  lec->J[peqn][pvar][i][j] += mass + source;

	} // End of loop over DOF (j)
	
      } // End of FILL sensitivities
      
      // Assemble sensitivities for MESH_DISPLACEMENT
      var = MESH_DISPLACEMENT1;
      if (pd->v[var]) {
  	pvar = upd->vp[var];

	// Loop over dimensions of mesh displacement 
	for ( b = 0; b < DIM; b++) {
	  var = MESH_DISPLACEMENT1 + b;
	  pvar = upd->vp[var];

	  // Loop over DOF (j)
	  for ( j = 0; j < ei->dof[var]; j++) {
	    jk = dof_map[j];
	    
	    // Load basis functions
	    phi_j = bf[eqn]->phi[j];
	    
	    // Assemble mass terms
	    mass = 0.0;
	    if ( mass_on ) {
	      mass += fv_dot->sh_sat_gasn * phi_i;
	    }
	    mass *= wt * h3 * fv->dsurfdet_dx[b][jk] * pd->etm[eqn][(LOG2_MASS)];
	    
	    // Assemble source terms
	    source = 0.0;
	    if ( source_on ) {
	      source -= flux * phi_i;
	    }
	    source *= wt * h3 * fv->dsurfdet_dx[b][jk] * pd->etm[eqn][(LOG2_SOURCE)];
	    
	    // Assemble full residual
	    lec->J[peqn][pvar][i][jk] += mass + source;
	    
	  } // End of loop over DOF (j)

	} // End of loop over mesh variables
	  
      } // End of MESH_DISPLACEMENT sensitivities

    }  // End of loop over DOF (i)

  } // End of jacobian assembly of R_SHELL_SAT_GASN equation
   

  // Cleanup 
  fv->wt = wt_old;
  safe_free((void *) n_dof);

  // Finish function
  return(status);

} // End of assemble_porous_shell_gasn



/*****************************************************************************/
/***assemble_porous_shell_open************************************************/
/*  _______________________________________________________________________  */

/* assemble_porous_shell_open -- Assemble terms (Residual & Jacobian) for a
 *                               structured porous shell equation for open
 *                               pores
 *
 * Created:     Friday March 26, 2010 sarober@sandia.gov
 *
 */
/*ARGSUSED*/ 
int
assemble_porous_shell_open(
			   dbl tt,                     // Time integration form
			   dbl dt,                     // Time step size
			   dbl xi[DIM],                // Current coordinates
			   const Exo_DB *exo           // ExoII handle
) {

  /* --- Initialization -----------------------------------------------------*/    

  // Variable definitions
  int eqn, peqn, var, pvar;                       // Equation / variables
  int i, j, a, b;                                 // Counter variables
  dbl phi_i, grad_phi_i[DIM], gradII_phi_i[DIM];  // Basis funcitons (i)
  dbl d_gradII_phi_i_dmesh[DIM][DIM][MDE];
  dbl phi_j, grad_phi_j[DIM], gradII_phi_j[DIM];  // Basis funcitons (j)
  dbl d_gradII_phi_j_dmesh[DIM][DIM][MDE];
  dbl mass, diff, sour;                           // Residual terms

  // Initialize output status function
  int status = 0;

  // Bail out of function if there is nothing to do
  eqn = R_SHELL_SAT_OPEN;

  if (!pd->e[eqn]) return(status);

  // Setup lubrication
  int *n_dof = NULL;
  int dof_map[MDE];
  dbl wt_old = fv->wt;
  n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
  lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

  // Unpack FEM variables from global structures
  dbl wt = fv->wt;                                // Gauss point weight
  dbl h3 = fv->h3;                                // Differential volume element
  dbl det_J = fv->sdet;                           // Jacobian of transformation
  dbl dA = det_J * wt * h3;
  dbl etm_mass = pd->etm[eqn][(LOG2_MASS)];
  dbl etm_diff = pd->etm[eqn][(LOG2_DIFFUSION)];
  dbl etm_sour = pd->etm[eqn][(LOG2_SOURCE)];

  // Load liquid material properties
  dbl mu = mp->viscosity;                         // Viscosity

  // Load porous medium parameters
  dbl phi = mp->porosity;                         // Porosity
  dbl H =  porous_shell_closed_height_model();    // Pore height (vertical)

  // Load field variables - PRS NOTE: NEED cross BC for integrating the two (set-up-shop)
  //  dbl P = fv->sh_p_open;                          // Porous pressure

  struct Level_Set_Data *ls_old;
    
  /* --- Calculate equation components ---------------------------------------*/    


  if (mp->SaturationModel != SHELL_TANH && 
      mp->SaturationModel != TANH &&
      mp->SaturationModel != TANH_EXTERNAL &&
      mp->SaturationModel != TANH_HYST) 
    {
      EH(-1,"Pacito problema: Only shell_tanh and tanh and tanh_external and tanh_hyst  models available for shell open pore. ");
      // PRS: just need to expand the nodal call on the mass term 
    }

  dbl S, dSdP;
  S = mp->saturation;
  dSdP = mp->d_saturation[SHELL_PRESS_OPEN];

  // Load permeability tensor
  dbl kappa = 0.0;

  // May want to add this as a kappa model. 
  // for ( i = 0; i < DIM; i++) kappa -= fv->snormal[i]*mp->perm_tensor[2][i];

  if (mp->PorousShellCrossKappaModel == CONSTANT)
    {
      kappa                 = mp->PorousShellCrossKappa;
    }
  else if (mp->PorousShellCrossKappaModel == EXTERNAL_FIELD)
    {
      EH(mp->Xperm_external_field_index, "Cross Permeability external field not found!");
      kappa = mp->PorousShellCrossKappa = 
	mp->u_PorousShellCrossKappa_function_constants[0]*fv->external_field[mp->Xperm_external_field_index];
      if (pd->TimeIntegration == TRANSIENT)
	{
	  mp_old->PorousShellCrossKappa =mp->u_PorousShellCrossKappa_function_constants[0]*fv->external_field[mp->Xperm_external_field_index];
	}
    }
  else
    {
      EH(-1,"Unrecognized Cross Kappa model");
    }

  // Load heaviside for level set weighting
  dbl Hside = 1.0, d_Hside_dF[DIM] = {0.0};
  if ( pd->v[FILL] ) {
    load_lsi( ls->Length_Scale );
    Hside = 1 - lsi->Hn;
    for ( i = 0; i < ei->dof[FILL]; i++) d_Hside_dF[i] = -lsi->d_Hn_dF[i];
  }

  // Rotate gradient
  dbl gradIIp[DIM] = {0.0};
  for ( i = 0; i < DIM; i++) {
    for ( j = 0; j < DIM; j++) {
      gradIIp[i] += fv->grad_sh_p_open[j] * delta(i,j);
      gradIIp[i] -= fv->grad_sh_p_open[j] * fv->snormal[i] * fv->snormal[j];
    }
  }

  dbl E_MASS[MDE], E_MASS_P[MDE];
  dbl Pnode, dSdPnode, dSdP_Pnode;
  dbl d_cap_pres[2], cap_pres;
  d_cap_pres[0] = d_cap_pres[1] = 0.;
  dbl Patm = mp->PorousShellPatm;

  for ( i = 0; i < ei->dof[eqn]; i++) {
    Pnode = *esp->sh_p_open[i];
    cap_pres = Patm - Pnode;
    //Snode = shell_saturation_pressure_curve(Pnode, &dSdPnode, &dSdP_Pnode);
   
    /* CHECK FOR REMOVAL */
    load_saturation(phi, cap_pres, d_cap_pres);
    dSdPnode = mp->d_saturation[SHELL_PRESS_OPEN];
    dSdP_Pnode = mp->d_d_saturation[SHELL_PRESS_OPEN][SHELL_PRESS_OPEN];

    E_MASS[i]    = H * phi * dSdPnode   * *esp_dot->sh_p_open[i];
    E_MASS_P[i]  = H * phi * dSdPnode   * (1.0+2.0*tt)/dt;
    E_MASS_P[i] += H * phi * dSdP_Pnode * *esp_dot->sh_p_open[i];
  }

  // Calculate DIFFUSION terms
  dbl E_DIFF[DIM] = {0.0};
  dbl E_DIFF_P[DIM][DIM] = {{0.0}};
  for ( a = 0; a < DIM; a++) {
    for ( b = 0; b < DIM; b++) {
      E_DIFF[a]      += -H / mu * mp->perm_tensor[a][b] * fv->grad_sh_p_open[b];
      E_DIFF_P[a][b] += -H / mu * mp->perm_tensor[a][b];
    }
  }

  // Calculate SOURCE term
  dbl E_SOUR, E_SOUR_P, E_SOUR_PLUB;
  dbl E_SOUR_F[MDE] = {0.0};
  dbl E_SOUR_2, E_SOUR_P_2, E_SOUR_2_PLUB_2;
  dbl E_SOUR_2_PF[MDE] = {0.0};
  dbl Pmin, Peff, Peff2;

  Pmin = mp->PorousShellInitPorePres; // PRS Note: ahh, this is the key lynch pin in the re-emergence
                                      // criteria.  This should not be const.  We have tofigure this out.
  Peff  = fv->lubp*Hside + Pmin*(1-Hside);

  E_SOUR      = kappa / mu * (fv->sh_p_open - Peff) / (2*S*H);
  E_SOUR_P    = kappa / mu / (2*S*H);
  E_SOUR_P   -= kappa / mu * (fv->sh_p_open - Peff) / (2*pow(S,2)*H) * dSdP;
  E_SOUR_PLUB =-kappa / mu / (2*S*H) * Hside;
  if ( pd->e[R_FILL] ) {
    for ( i = 0; i < ei->dof[FILL]; i++) {
      E_SOUR_F[i] = kappa / mu * (-1) / (2*S*H) * (fv->lubp - Pmin) * d_Hside_dF[i];
    }
  }

  // HACK to keep liquid from being sucked back to the lubrication layer
  if ( E_SOUR > 0.0 ) {
    // E_SOUR = 0.0;
    //  E_SOUR_P = 0.0;
    // E_SOUR_PLUB = 0.0;
    // for ( i = 0; i < ei->dof[FILL]; i++) E_SOUR_F[i] = 0.0;
  }

  // OK, if we have multilayer, we need to get source terms from the second-story lubrication field, whose 
  // footprint is dictated by a phase field.   We could use setup_shop but the element friend thing is not 
  // function with this really pathological shell-on-shell stack.   So we will have to do things by hand.

  dbl Hside_2 = 1.0, d_Hside_2_dpF[DIM] = {0.0}; dbl lubp_2 = 0.0; dbl pF=0.0;
  E_SOUR_2 = 0.; E_SOUR_P_2 = 0.; E_SOUR_2_PLUB_2 = 0.;
  if(upd->ep[R_LUBP_2] >=0 )
    {
      //Compute Heaviside function for phase field
      ls_old = ls;
      if(pfd != NULL) ls = pfd->ls[0]; 
      if ( upd->vp[PHASE1] >= 0) {
	load_lsi_shell_second( ls->Length_Scale );
	Hside_2 = 1 - lsi->Hn;
	for ( i = 0; i < ei->dof[FILL]; i++) d_Hside_2_dpF[i] = -lsi->d_Hn_dF[i];
      }

      //Compute lubp_2 field and phase field variables by hand here since we are not in the right element.
      lubp_2 = 0.0; pF=0.0;
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  lubp_2 +=    *esp->lubp_2[i] * bf[R_LUBP_2]->phi[i];
	  pF +=    *esp->pF[0][i] * bf[R_PHASE1]->phi[i];
	} 


      Peff2=  lubp_2*Hside_2 + Pmin*(1-Hside_2);
      E_SOUR_2   = kappa / mu * (fv->sh_p_open - Peff2) / (2*S*H);
      E_SOUR_P_2    = kappa / mu / (2*S*H);
      E_SOUR_P_2   -= kappa / mu * (fv->sh_p_open - Peff) / (2*pow(S,2)*H) * dSdP;
      E_SOUR_2_PLUB_2 = -kappa / mu / (2*S*H) * Hside_2;
	if (upd->ep[R_PHASE1] >= 0) {
	  for ( i = 0; i < ei->dof[FILL]; i++) {
	    E_SOUR_2_PF[i] = kappa / mu * (-1) / (2*S*H) * (lubp_2 - Pmin) * d_Hside_2_dpF[i];
	  }
	}

      //Now convert back to level set field   
      ls = ls_old;
      if ( upd->vp[FILL] >= 0) {
	load_lsi( ls->Length_Scale );
      }

    }

  // Assemble test for LS weight
  int mytest[MDE];
  for ( i = 0; i < ei->dof[eqn]; i++) {
    if ( pd->e[R_FILL] ) {
      mytest[i] = fv->F < 0;
      //mytest[i] = Hside > 0.95;
      //mytest[i] = 1;
    } else {
      mytest[i] = 1;
    }
  }
  // Assemble test for PHASE1 weight
  int mytest_2[MDE];
  for ( i = 0; i < ei->dof[eqn]; i++) {
    if ( upd->ep[R_PHASE1] >= 0) {
      mytest_2[i] = pF < 0;
      //mytest[i] = Hside > 0.95;
      //mytest[i] = 1;
    } else {
      mytest_2[i] = 1;
    }
  }

  /* --- Assemble residuals --------------------------------------------------*/

  // Assemble residual contribution to this equation
  eqn = R_SHELL_SAT_OPEN;
  if (af->Assemble_Residual) {
    peqn = upd->ep[eqn];
    
    // Loop over DOF (i)
    for ( i = 0; i < ei->dof[eqn]; i++) {         
      
      // Load basis functions
      ShellBF( eqn, i, &phi_i, grad_phi_i, gradII_phi_i, d_gradII_phi_i_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
      
      // Assemble mass term
      mass = 0.0;
      if ( T_MASS ) {
	mass += E_MASS[i] * phi_i;
      }
      mass *= dA * etm_mass;
      
      // Assemble diffusion term
      diff = 0.0;
      if ( T_DIFFUSION ) {
	for ( a = 0; a < DIM; a++) {
	  diff -= E_DIFF[a] * gradII_phi_i[a];
	}
      }
      diff *= dA * etm_diff;
      
      // Assemble source term
      sour = 0.0;
      if ( T_SOURCE ) {
	sour += (E_SOUR*mytest[i] + E_SOUR_2*mytest_2[i]) * phi_i;
      }
      sour *= dA * etm_sour;

      // Assemble full residual
      lec->R[peqn][i] += mass + diff + sour;
      
    }  // End of loop over DOF (i)

  } // End of residual assembly of R_SHELL_SAT_OPEN

  // Assemble residual contribution to the lubrication equation
  eqn = R_LUBP;
  if (af->Assemble_Residual & pd->e[eqn]) {
    peqn = upd->ep[eqn];
    
    // Loop over DOF (i)
    for ( i = 0; i < ei->dof[eqn]; i++) {         
      
      // Load basis functions
      ShellBF( eqn, i, &phi_i, grad_phi_i, gradII_phi_i, d_gradII_phi_i_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
            
      // Assemble source term
      sour = 0.0;
      if ( T_SOURCE ) {
	sour += E_SOUR * phi_i;
      }
      sour *= dA * etm_sour * mytest[i];

      // Assemble full residual
      lec->R[peqn][i] += sour;
      
    }  // End of loop over DOF (i)

  } // End of residual assembly of R_LUBP

  // Assemble residual contribution to the second lubrication equation
  eqn = R_LUBP_2;
  if ((af->Assemble_Residual & upd->ep[eqn]) >= 0) {
    peqn = upd->ep[eqn];
    
    // Loop over DOF (i)
    for ( i = 0; i < ei->dof[eqn]; i++) {         
      
      // Load basis functions
      ShellBF( eqn, i, &phi_i, grad_phi_i, gradII_phi_i, d_gradII_phi_i_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
            
      // Assemble source term
      sour = 0.0;
      if ( T_SOURCE ) {
	sour += E_SOUR_2 * phi_i;
      }
      sour *= dA * etm_sour * mytest_2[i];

      // Assemble full residual
      lec->R[peqn][i] += sour;
      
    }  // End of loop over DOF (i)

  } // End of residual assembly of R_LUBP_2


  /* --- Assemble Jacobian --------------------------------------------------*/
  eqn = R_SHELL_SAT_OPEN;
  if (af->Assemble_Jacobian) {
    peqn = upd->ep[eqn];
    
    // Loop over DOF (i)
    for ( i = 0; i < ei->dof[eqn]; i++) {
      
      // Load basis functions
      ShellBF( eqn, i, &phi_i, grad_phi_i, gradII_phi_i, d_gradII_phi_i_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
      
      // Assemble sensitivities for SHELL_PRESS_OPEN
      var = SHELL_PRESS_OPEN;
      if (pd->v[var]) {
  	pvar = upd->vp[var];

  	// Loop over DOF (j)
  	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  ShellBF( var, j, &phi_j, grad_phi_j, gradII_phi_j, d_gradII_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
      
	  // Assemble mass term
	  mass = 0.0;
	  if ( T_MASS ) {
	    //mass += E_MASS_P * phi_i * phi_j;
	    if ( i == j ) mass += E_MASS_P[i] * phi_i;
	  }
	  mass *= dA * etm_mass;
	  
	  // Assemble diffusion term
	  diff = 0.0;
	  if ( T_DIFFUSION ) {
	    for ( a = 0; a < DIM; a++) {
	      for ( b = 0; b < DIM; b++) {
		diff -= E_DIFF_P[a][b] * gradII_phi_i[a] * gradII_phi_j[b];
	      }
	    }
	  }
	  diff *= dA * etm_diff;
	  
	  // Assemble source term
	  sour = 0.0;
	  if ( T_SOURCE ) {   
	    sour += (E_SOUR_P*mytest[i] + E_SOUR_P_2*mytest_2[i]) * phi_i * phi_j;
	  }
	  sour *= dA * etm_sour;
	  
	  // Assemble full Jacobian
	  lec->J[peqn][pvar][i][j] += mass + diff + sour;
	  
	} // End of loop over DOF (j)
	
      } // End of SHELL_PRESS_OPEN sensitivities
      
      // Assemble sensitivities for LUBP
      var = LUBP;
      if (pd->v[var]) {
  	pvar = upd->vp[var];

  	// Loop over DOF (j)
  	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  ShellBF( var, j, &phi_j, grad_phi_j, gradII_phi_j, d_gradII_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
	  
	  // Assemble source term
	  sour = 0.0;
	  if ( T_SOURCE ) {   
	    sour += E_SOUR_PLUB * phi_i * phi_j;
	  }
	  sour *= dA * etm_sour * mytest[i];
	  
	  // Assemble full Jacobian
	  lec->J[peqn][pvar][i][j] += sour;
	  
	} // End of loop over DOF (j)
	
      } // End of LUBP sensitivities

      // Assemble sensitivities for LUBP
      var = LUBP_2;
      if (upd->vp[var] >= 0) {
  	pvar = upd->vp[var];

  	// Loop over DOF (j)
  	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  ShellBF( var, j, &phi_j, grad_phi_j, gradII_phi_j, d_gradII_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
	  
	  // Assemble source term
	  sour = 0.0;
	  if ( T_SOURCE ) {   
	    sour += E_SOUR_2_PLUB_2 * phi_i * phi_j;
	  }
	  sour *= dA * etm_sour * mytest_2[i];
	  
	  // Assemble full Jacobian
	  lec->J[peqn][pvar][i][j] += sour;
	  
	} // End of loop over DOF (j)
	
      } // End of LUBP sensitivities
      
      // Assemble sensitivities for FILL
      var = FILL;
      if (pd->v[var]) {
  	pvar = upd->vp[var];

  	// Loop over DOF (j)
  	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  ShellBF( var, j, &phi_j, grad_phi_j, gradII_phi_j, d_gradII_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
	  
	  // Assemble source term
	  sour = 0.0;
	  if ( T_SOURCE ) {   
	    sour += E_SOUR_F[j] * phi_i;
	  }
	  sour *= dA * etm_sour * mytest[i];
	  
	  // Assemble full Jacobian
	  lec->J[peqn][pvar][i][j] += sour;
	  
	} // End of loop over DOF (j)
	
      } // End of FILL sensitivities

 // Assemble sensitivities for PHASE1
      var = PHASE1;
      if (upd->vp[var] >=0) {
  	pvar = upd->vp[var];

  	// Loop over DOF (j)
  	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  ShellBF( var, j, &phi_j, grad_phi_j, gradII_phi_j, d_gradII_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
	  
	  // Assemble source term
	  sour = 0.0;
	  if ( T_SOURCE ) {   
	    sour += E_SOUR_2_PF[j] * phi_i;
	  }
	  sour *= dA * etm_sour * mytest_2[i];
	  
	  // Assemble full Jacobian
	  lec->J[peqn][pvar][i][j] += sour;
	  
	} // End of loop over DOF (j)
	
      } // End of PHASE1 sensitivities
                 
    } // End of loop over DOF (i)

  } // End of Jacobian assembly of R_SHELL_SAT_OPEN

  eqn = R_LUBP;
  if (af->Assemble_Jacobian & pd->e[eqn]) {
    peqn = upd->ep[eqn];
    
    // Loop over DOF (i)
    for ( i = 0; i < ei->dof[eqn]; i++) {
      
      // Load basis functions
      ShellBF( eqn, i, &phi_i, grad_phi_i, gradII_phi_i, d_gradII_phi_i_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
      
      // Assemble sensitivities for SHELL_PRESS_OPEN
      var = SHELL_PRESS_OPEN;
      if (pd->v[var]) {
  	pvar = upd->vp[var];

  	// Loop over DOF (j)
  	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  ShellBF( var, j, &phi_j, grad_phi_j, gradII_phi_j, d_gradII_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
 	  
	  // Assemble source term
	  sour = 0.0;
	  if ( T_SOURCE ) {   
	    sour += E_SOUR_P * phi_i * phi_j;
	  }
	  sour *= dA * etm_sour * mytest[i];
	  
	  // Assemble full Jacobian
	  lec->J[peqn][pvar][i][j] += sour;
	  
	} // End of loop over DOF (j)
	
      } // End of SHELL_PRESS_OPEN sensitivities
      
      // Assemble sensitivities for LUBP
      var = LUBP;
      if (pd->v[var]) {
  	pvar = upd->vp[var];

  	// Loop over DOF (j)
  	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  ShellBF( var, j, &phi_j, grad_phi_j, gradII_phi_j, d_gradII_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
	  
	  // Assemble source term
	  sour = 0.0;
	  if ( T_SOURCE ) {   
	    sour += E_SOUR_PLUB * phi_i * phi_j;
	  }
	  sour *= dA * etm_sour * mytest[i];
	  
	  // Assemble full Jacobian
	  lec->J[peqn][pvar][i][j] += sour;
	  
	} // End of loop over DOF (j)
	
      } // End of LUBP sensitivities
      
      // Assemble sensitivities for FILL
      var = FILL;
      if (pd->v[var]) {
  	pvar = upd->vp[var];

  	// Loop over DOF (j)
  	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  ShellBF( var, j, &phi_j, grad_phi_j, gradII_phi_j, d_gradII_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
	  
	  // Assemble source term
	  sour = 0.0;
	  if ( T_SOURCE ) {   
	    sour += E_SOUR_F[j] * phi_i;
	  }
	  sour *= dA * etm_sour * mytest[i];
	  
	  // Assemble full Jacobian
	  lec->J[peqn][pvar][i][j] += sour;
	  
	} // End of loop over DOF (j)
	
      } // End of FILL sensitivities
                 
    } // End of loop over DOF (i)

  } // End of Jacobian assembly of R_LUBP

  eqn = R_LUBP_2;
  if ((af->Assemble_Jacobian & upd->ep[eqn]) >= 0) {
    peqn = upd->ep[eqn];
    
    // Loop over DOF (i)
    for ( i = 0; i < ei->dof[eqn]; i++) {
      
      // Load basis functions
      ShellBF( eqn, i, &phi_i, grad_phi_i, gradII_phi_i, d_gradII_phi_i_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
      
      // Assemble sensitivities for SHELL_PRESS_OPEN
      var = SHELL_PRESS_OPEN;
      if (pd->v[var]) {
  	pvar = upd->vp[var];

  	// Loop over DOF (j)
  	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  ShellBF( var, j, &phi_j, grad_phi_j, gradII_phi_j, d_gradII_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
 	  
	  // Assemble source term
	  sour = 0.0;
	  if ( T_SOURCE ) {   
	    sour += E_SOUR_P_2 * phi_i * phi_j;
	  }
	  sour *= dA * etm_sour * mytest_2[i];
	  
	  // Assemble full Jacobian
	  lec->J[peqn][pvar][i][j] += sour;
	  
	} // End of loop over DOF (j)
	
      } // End of SHELL_PRESS_OPEN sensitivities
      
      // Assemble sensitivities for LUBP
      var = LUBP_2;
      if (upd->vp[var] >=0) {
  	pvar = upd->vp[var];

  	// Loop over DOF (j)
  	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  ShellBF( var, j, &phi_j, grad_phi_j, gradII_phi_j, d_gradII_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
	  
	  // Assemble source term
	  sour = 0.0;
	  if ( T_SOURCE ) {   
	    sour += E_SOUR_2_PLUB_2 * phi_i * phi_j;
	  }
	  sour *= dA * etm_sour * mytest_2[i];
	  
	  // Assemble full Jacobian
	  lec->J[peqn][pvar][i][j] += sour;
	  
	} // End of loop over DOF (j)
	
      } // End of LUBP sensitivities
      
      // Assemble sensitivities for FILL
      var = PHASE1;
      if (upd->vp[var] >= 0) {
  	pvar = upd->vp[var];

  	// Loop over DOF (j)
  	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  ShellBF( var, j, &phi_j, grad_phi_j, gradII_phi_j, d_gradII_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
	  
	  // Assemble source term
	  sour = 0.0;
	  if ( T_SOURCE ) {   
	    sour += E_SOUR_2_PF[j] * phi_i;
	  }
	  sour *= dA * etm_sour * mytest_2[i];
	  
	  // Assemble full Jacobian
	  lec->J[peqn][pvar][i][j] += sour;
	  
	} // End of loop over DOF (j)
	
      } // End of FILL sensitivities
                 
    } // End of loop over DOF (i)

  } // End of Jacobian assembly of R_LUBP_2

  // Finalize and exit
  fv->wt = wt_old;
  safe_free((void *) n_dof);
  return(status);

} // End of assemble_porous_shell_open

/*****************************************************************************/
/*****************************************************************************/
/***assemble_porous_shell_open_2************************************************/
/*  _______________________________________________________________________  */

/* assemble_porous_shell_open_2 -- Assemble terms (Residual & Jacobian) for a
 *                               structured porous shell equation for open
 *                               pores
 *
 * Created:     Friday August 3, 2012 prschun
 *
 */
/*ARGSUSED*/ 
int
assemble_porous_shell_open_2(
			   dbl tt,                     // Time integration form
			   dbl dt,                     // Time step size
			   dbl xi[DIM],                // Current coordinates
			   const Exo_DB *exo           // ExoII handle
			   ) 
  /* PRS Note:  this open porous shell routine only interacts with R_LUBP_2 and
   * phase-fields.   That is, unlike its brother assemble_porous_shell_open, which interacts
   * with R_LUBP and eventually R_LUBP_2 from the other side of the layered stack, this 
   * routine is for a shell that has no other lubrication layer above it, for now.  (8/16/2012)
   */
{

  /* --- Initialization -----------------------------------------------------*/    

  // Variable definitions
  int eqn, peqn, var, pvar;                       // Equation / variables
  int i, j, a, b;                                 // Counter variables
  dbl phi_i, grad_phi_i[DIM], gradII_phi_i[DIM];  // Basis funcitons (i)
  dbl d_gradII_phi_i_dmesh[DIM][DIM][MDE];
  dbl phi_j, grad_phi_j[DIM], gradII_phi_j[DIM];  // Basis funcitons (j)
  dbl d_gradII_phi_j_dmesh[DIM][DIM][MDE];
  dbl mass, diff, sour;                           // Residual terms

  // Initialize output status function
  int status = 0;

  // Bail out of function if there is nothing to do
  eqn = R_SHELL_SAT_OPEN_2;

  if (!pd->e[eqn]) return(status);

  // For that matter, also bail out if there is no R_LUBP_2 equation here, as then it is 
  // a moot point even being in here;

  if (!pd->e[R_LUBP_2]) return(status);

  // Setup lubrication
  int *n_dof = NULL;
  int dof_map[MDE];
  dbl wt_old = fv->wt;
  n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
  lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

  // Unpack FEM variables from global structures
  dbl wt = fv->wt;                                // Gauss point weight
  dbl h3 = fv->h3;                                // Differential volume element
  dbl det_J = fv->sdet;                           // Jacobian of transformation
  dbl dA = det_J * wt * h3;
  dbl etm_mass = pd->etm[eqn][(LOG2_MASS)];
  dbl etm_diff = pd->etm[eqn][(LOG2_DIFFUSION)];
  dbl etm_sour = pd->etm[eqn][(LOG2_SOURCE)];

  // Load liquid material properties
  dbl mu;                         // Viscosity

  // Load porous medium parameters
  dbl phi = mp->porosity;                         // Porosity
  dbl H =  porous_shell_closed_height_model();    // Pore height (vertical)

  // Load field variables - PRS NOTE: NEED cross BC for integrating the two (set-up-shop)

  struct Level_Set_Data *ls_old;

  /* Load up properties */
  mu = mp->viscosity;
  /* OK this is a cludge for the second ply.  the viscosity function in the gap needs to be modulated,but here we just need 
   * the liquid imbibing phase viscosity.  PRS (9/21/2012)
   */
  if (mp->ViscosityModel == CONST_PHASE_FUNCTION) mu = mp->u_viscosity[2];
    
  /* --- Calculate equation components ---------------------------------------*/    


  if (mp->SaturationModel != SHELL_TANH && 
      mp->SaturationModel != TANH &&
      mp->SaturationModel != TANH_EXTERNAL &&
      mp->SaturationModel != TANH_HYST) 
    {
      EH(-1,"Pacito problema: Only shell_tanh, tanh, tanh_external, and tanh_hyst model available for shell open pore. Not much work to remedy this, though");
      // PRS: just need to expand the nodal call on the mass term 
    }

  dbl S, dSdP;
  S = mp->saturation;
  dSdP = mp->d_saturation[SHELL_PRESS_OPEN_2];

  // Load permeability tensor
  dbl kappa = 0.0;

  if (mp->PorousShellCrossKappaModel == CONSTANT)
    {
      kappa                 = mp->PorousShellCrossKappa;
    }
  else if (mp->PorousShellCrossKappaModel == EXTERNAL_FIELD)
    {
      EH(mp->Xperm_external_field_index, "Cross Permeability external field not found!");
      kappa = mp->PorousShellCrossKappa = fv->external_field[mp->Xperm_external_field_index];
      if (pd->TimeIntegration == TRANSIENT)
	{
	  mp_old->PorousShellCrossKappa = fv->external_field[mp->Xperm_external_field_index];
	}
    }
  else
    {
      EH(-1,"Unrecognized Cross Kappa model");
    }

  // Load heaviside for phase-field  weighting
  dbl Hside = 1.0, d_Hside_dF[DIM] = {0.0};
  if ( pd->v[PHASE1] ) {
    ls_old = ls;
    if(pfd != NULL) ls = pfd->ls[0]; 
    load_lsi( ls->Length_Scale );
    Hside = 1 - lsi->Hn;
    for ( i = 0; i < ei->dof[PHASE1]; i++) d_Hside_dF[i] = -lsi->d_Hn_dF[i];
    ls = ls_old;
  }

  // Rotate gradient
  dbl gradIIp[DIM] = {0.0};
  for ( i = 0; i < DIM; i++) {
    for ( j = 0; j < DIM; j++) {
      gradIIp[i] += fv->grad_sh_p_open_2[j] * delta(i,j);
      gradIIp[i] -= fv->grad_sh_p_open_2[j] * fv->snormal[i] * fv->snormal[j];
    }
  }

  dbl E_MASS[MDE], E_MASS_P[MDE];
  dbl Pnode, dSdPnode, dSdP_Pnode;
  dbl d_cap_pres[2], cap_pres;
  d_cap_pres[0] = d_cap_pres[1] = 0.;
  dbl Patm = mp->PorousShellPatm;

  for ( i = 0; i < ei->dof[eqn]; i++) {
    Pnode = *esp->sh_p_open_2[i];
    cap_pres = Patm - Pnode;
    //Snode = shell_saturation_pressure_curve(Pnode, &dSdPnode, &dSdP_Pnode);
   
    /* CHECK FOR REMOVAL */
    load_saturation(phi, cap_pres, d_cap_pres);
    dSdPnode = mp->d_saturation[SHELL_PRESS_OPEN_2];
    dSdP_Pnode = mp->d_d_saturation[SHELL_PRESS_OPEN_2][SHELL_PRESS_OPEN_2];

    E_MASS[i]    = H * phi * dSdPnode   * *esp_dot->sh_p_open_2[i];
    E_MASS_P[i]  = H * phi * dSdPnode   * (1.0+2.0*tt)/dt;
    E_MASS_P[i] += H * phi * dSdP_Pnode * *esp_dot->sh_p_open_2[i];
  }

  // Calculate DIFFUSION terms
  dbl E_DIFF[DIM] = {0.0};
  dbl E_DIFF_P[DIM][DIM] = {{0.0}};
  for ( a = 0; a < DIM; a++) {
    for ( b = 0; b < DIM; b++) {
      E_DIFF[a]      += -H / mu * mp->perm_tensor[a][b] * fv->grad_sh_p_open_2[b];
      E_DIFF_P[a][b] += -H / mu * mp->perm_tensor[a][b];
    }
  }

  // Calculate SOURCE term
  dbl E_SOUR, E_SOUR_P, E_SOUR_PLUB;
  dbl E_SOUR_F[MDE] = {0.0};
  dbl Pmin, Peff;

  Pmin = mp->PorousShellInitPorePres; // PRS Note: ahh, this is the key lynch pin in the re-emergence
                                      // criteria.  This should not be const.  We have tofigure this out.
  Peff  = fv->lubp_2*Hside + Pmin*(1-Hside);
  E_SOUR      = kappa / mu * (fv->sh_p_open_2 - Peff) / (2*S*H);
  E_SOUR_P    = kappa / mu / (2*S*H);
  E_SOUR_P   -= kappa / mu * (fv->sh_p_open_2 - Peff) / (2*pow(S,2)*H) * dSdP;
  E_SOUR_PLUB =-kappa / mu / (2*S*H) * Hside;
  if ( pd->e[R_PHASE1] ) {
    for ( i = 0; i < ei->dof[R_PHASE1]; i++) {
      E_SOUR_F[i] = kappa / mu * (-1) / (2*S*H) * (fv->lubp_2 - Pmin) * d_Hside_dF[i];
    }
  }

  // HACK to keep liquid from being sucked back to the lubrication layer
  if ( E_SOUR > 0.0 ) {
    // E_SOUR = 0.0;
    //  E_SOUR_P = 0.0;
    // E_SOUR_PLUB = 0.0;
    // for ( i = 0; i < ei->dof[PHASE1]; i++) E_SOUR_F[i] = 0.0;
  }

  // Assemble test for LS weight
  int mytest[MDE];
  for ( i = 0; i < ei->dof[eqn]; i++) {
    if ( pd->e[R_PHASE1] ) {
      mytest[i] = fv->pF[0] < 0;
      //mytest[i] = Hside > 0.95;
      //mytest[i] = 1;
    } else {
      mytest[i] = 1;
    }
  }

  /* --- Assemble residuals --------------------------------------------------*/

  // Assemble residual contribution to this equation
  eqn = R_SHELL_SAT_OPEN_2;
  if (af->Assemble_Residual) {
    peqn = upd->ep[eqn];
    
    // Loop over DOF (i)
    for ( i = 0; i < ei->dof[eqn]; i++) {         
      
      // Load basis functions
      ShellBF( eqn, i, &phi_i, grad_phi_i, gradII_phi_i, d_gradII_phi_i_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
      
      // Assemble mass term
      mass = 0.0;
      if ( T_MASS ) {
	mass += E_MASS[i] * phi_i;
      }
      mass *= dA * etm_mass;
      
      // Assemble diffusion term
      diff = 0.0;
      if ( T_DIFFUSION ) {
	for ( a = 0; a < DIM; a++) {
	  diff -= E_DIFF[a] * gradII_phi_i[a];
	}
      }
      diff *= dA * etm_diff;
      
      // Assemble source term
      sour = 0.0;
      if ( T_SOURCE ) {
	sour += E_SOUR * phi_i;
      }
      sour *= dA * etm_sour * mytest[i];

      // Assemble full residual
      lec->R[peqn][i] += mass + diff + sour;
      
    }  // End of loop over DOF (i)

  } // End of residual assembly of R_SHELL_SAT_OPEN

  // Assemble residual contribution to the lubrication equation
  eqn = R_LUBP_2;
  if (af->Assemble_Residual & pd->e[eqn]) {
    peqn = upd->ep[eqn];
    
    // Loop over DOF (i)
    for ( i = 0; i < ei->dof[eqn]; i++) {         
      
      // Load basis functions
      ShellBF( eqn, i, &phi_i, grad_phi_i, gradII_phi_i, d_gradII_phi_i_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
            
      // Assemble source term
      sour = 0.0;
      if ( T_SOURCE ) {
	sour += E_SOUR * phi_i;
      }
      sour *= dA * etm_sour * mytest[i];

      // Assemble full residual
      lec->R[peqn][i] += sour;
      
    }  // End of loop over DOF (i)

  } // End of residual assembly of R_LUBP_2

  /* PRS Note you will have another source here for the other layer  */

  /* --- Assemble Jacobian --------------------------------------------------*/
  eqn = R_SHELL_SAT_OPEN_2;
  if (af->Assemble_Jacobian) {
    peqn = upd->ep[eqn];
    
    // Loop over DOF (i)
    for ( i = 0; i < ei->dof[eqn]; i++) {
      
      // Load basis functions
      ShellBF( eqn, i, &phi_i, grad_phi_i, gradII_phi_i, d_gradII_phi_i_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
      
      // Assemble sensitivities for SHELL_PRESS_OPEN_2
      var = SHELL_PRESS_OPEN_2;
      if (pd->v[var]) {
  	pvar = upd->vp[var];

  	// Loop over DOF (j)
  	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  ShellBF( var, j, &phi_j, grad_phi_j, gradII_phi_j, d_gradII_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
      
	  // Assemble mass term
	  mass = 0.0;
	  if ( T_MASS ) {
	    //mass += E_MASS_P * phi_i * phi_j;
	    if ( i == j ) mass += E_MASS_P[i] * phi_i;
	  }
	  mass *= dA * etm_mass;
	  
	  // Assemble diffusion term
	  diff = 0.0;
	  if ( T_DIFFUSION ) {
	    for ( a = 0; a < DIM; a++) {
	      for ( b = 0; b < DIM; b++) {
		diff -= E_DIFF_P[a][b] * gradII_phi_i[a] * gradII_phi_j[b];
	      }
	    }
	  }
	  diff *= dA * etm_diff;
	  
	  // Assemble source term
	  sour = 0.0;
	  if ( T_SOURCE ) {   
	    sour += E_SOUR_P * phi_i * phi_j;
	  }
	  sour *= dA * etm_sour * mytest[i];
	  
	  // Assemble full Jacobian
	  lec->J[peqn][pvar][i][j] += mass + diff + sour;
	  
	} // End of loop over DOF (j)
	
      } // End of SHELL_PRESS_OPEN_2 sensitivities
      
      // Assemble sensitivities for LUBP_2
      var = LUBP_2;
      if (pd->v[var]) {
  	pvar = upd->vp[var];

  	// Loop over DOF (j)
  	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  ShellBF( var, j, &phi_j, grad_phi_j, gradII_phi_j, d_gradII_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
	  
	  // Assemble source term
	  sour = 0.0;
	  if ( T_SOURCE ) {   
	    sour += E_SOUR_PLUB * phi_i * phi_j;
	  }
	  sour *= dA * etm_sour * mytest[i];
	  
	  // Assemble full Jacobian
	  lec->J[peqn][pvar][i][j] += sour;
	  
	} // End of loop over DOF (j)
	
      } // End of LUBP_2 sensitivities
      
      // Assemble sensitivities for PHASE1
      var = PHASE1;
      if (pd->v[var]) {
  	pvar = upd->vp[var];

  	// Loop over DOF (j)
  	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  ShellBF( var, j, &phi_j, grad_phi_j, gradII_phi_j, d_gradII_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
	  
	  // Assemble source term
	  sour = 0.0;
	  if ( T_SOURCE ) {   
	    sour += E_SOUR_F[j] * phi_i;
	  }
	  sour *= dA * etm_sour * mytest[i];
	  
	  // Assemble full Jacobian
	  lec->J[peqn][pvar][i][j] += sour;
	  
	} // End of loop over DOF (j)
	
      } // End of phase1 sensitivities
                 
    } // End of loop over DOF (i)

  } // End of Jacobian assembly of R_SHELL_SAT_OPEN_2

  eqn = R_LUBP_2;
  if (af->Assemble_Jacobian & pd->e[eqn]) {
    peqn = upd->ep[eqn];
    
    // Loop over DOF (i)
    for ( i = 0; i < ei->dof[eqn]; i++) {
      
      // Load basis functions
      ShellBF( eqn, i, &phi_i, grad_phi_i, gradII_phi_i, d_gradII_phi_i_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
      
      // Assemble sensitivities for SHELL_PRESS_OPEN
      var = SHELL_PRESS_OPEN_2;
      if (pd->v[var]) {
  	pvar = upd->vp[var];

  	// Loop over DOF (j)
  	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  ShellBF( var, j, &phi_j, grad_phi_j, gradII_phi_j, d_gradII_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
 	  
	  // Assemble source term
	  sour = 0.0;
	  if ( T_SOURCE ) {   
	    sour += E_SOUR_P * phi_i * phi_j;
	  }
	  sour *= dA * etm_sour * mytest[i];
	  
	  // Assemble full Jacobian
	  lec->J[peqn][pvar][i][j] += sour;
	  
	} // End of loop over DOF (j)
	
      } // End of SHELL_PRESS_OPEN sensitivities
      
      // Assemble sensitivities for LUBP_2
      var = LUBP_2;
      if (pd->v[var]) {
  	pvar = upd->vp[var];

  	// Loop over DOF (j)
  	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  ShellBF( var, j, &phi_j, grad_phi_j, gradII_phi_j, d_gradII_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
	  
	  // Assemble source term
	  sour = 0.0;
	  if ( T_SOURCE ) {   
	    sour += E_SOUR_PLUB * phi_i * phi_j;
	  }
	  sour *= dA * etm_sour * mytest[i];
	  
	  // Assemble full Jacobian
	  lec->J[peqn][pvar][i][j] += sour;
	  
	} // End of loop over DOF (j)
	
      } // End of LUBP sensitivities
      
      // Assemble sensitivities for PHASE1
      var = PHASE1;
      if (pd->v[var]) {
  	pvar = upd->vp[var];

  	// Loop over DOF (j)
  	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  ShellBF( var, j, &phi_j, grad_phi_j, gradII_phi_j, d_gradII_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
	  
	  // Assemble source term
	  sour = 0.0;
	  if ( T_SOURCE ) {   
	    sour += E_SOUR_F[j] * phi_i;
	  }
	  sour *= dA * etm_sour * mytest[i];
	  
	  // Assemble full Jacobian
	  lec->J[peqn][pvar][i][j] += sour;
	  
	} // End of loop over DOF (j)
	
      } // End of PHASE1 sensitivities
                 
    } // End of loop over DOF (i)

  } // End of Jacobian assembly of R_LUBP

  // Finalize and exit
  fv->wt = wt_old;
  safe_free((void *) n_dof);
  return(status);

} // End of assemble_porous_shell_open_2


/* End of file mm_fill_shell.c */

  /*
    shell_lubr_solid_struct_bc(): Balance of Shell-lubrication forces with
    abutting continuum solid, like for sliding-structure problems
  */
void
shell_lubr_solid_struct_bc(double func[DIM],
			   double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
			   const double x_dot[MAX_PDIM],
			   const double tt,
			   const double dt,
			   const int id_side,
			   const double wt,
			   double xi[DIM],
			   const Exo_DB *exo,
			   const double scale)
{
  int  j, jvar, p, var;
  int *n_dof = NULL;
                                                                                
  /* Unpack variables from structures for local convenience. */

  /*
   * Prepare geometry
   */
  int dof_map[MDE];
  n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
  lubrication_shell_initialize(n_dof, dof_map, id_side, xi, exo, 1);
                                                                                
  /* Assemble the residual equation */
  /* Here, variables from the remote elements are used. */
  if (af->Assemble_Residual)
    {
      for(p=0; p<pd->Num_Dim; p++)
	{
	  func[p] -= scale*fv->snormal[p]*fv->lubp;
	}                    
    }
  
  /* Assemble the Jacobian terms */
  /* Here, variables from the remote elements are used. */
  if (af->Assemble_Jacobian)
    {

      /* dfunc/dx */
      for(jvar=0; jvar<pd->Num_Dim; jvar++)
	{
	  var = MESH_DISPLACEMENT1 + jvar;
	  if (pd->v[var])
	    {
	      for (j = 0; j < ei->dof[var]; j++)
		{
		  for(p=0; p<pd->Num_Dim; p++)
		    {
		      d_func[p][var][j] -=
			scale*fv->dsnormal_dx[p][jvar][j]*fv->lubp;
		    }
		}
	    }
	}

      /*dfunc/dp */
      /* you need to really check this out. We are only on the solid-continuum
       * element and you are taking sensitivities wrt lubp, on the friend
       * shell element. pd[var] is negative here, so we pull a fast one with 
       * upd.  
       */
      var = LUBP;
      if (upd->vp[var] != -1)
	{
	  for (j = 0; j < ei->dof[var]; j++)
	    {
	      for(p=0; p<pd->Num_Dim; p++)
		{
		  d_func[p][var][j] -= scale*fv->snormal[p]*bf[var]->phi[j];
		}
	    }
	}

    }

  safe_free((void *) n_dof);
  return;
}  /* shell_lubr_solid_struct_bc */


/*****************************************************************************/
/***assemble_shell_deltah******************************************************/
/*  _______________________________________________________________________  */

/* assemble_shell_deltah -- assemble terms (Residual & Jacobian) for scalar shell
 *                          gap evolution equation.  This routine can be a user-routine
 *                          for any shell-gap evolution which depends on local shell variables.
 *                          Current version here written for melting problems using the 
 *                          shell temperature. 
 *
 * in:
 * 	ei -- pointer to Element Indices	structure
 *	pd -- pointer to Problem Description	structure
 *	af -- pointer to Action Flag		structure
 *	bf -- pointer to Basis Function		structure
 *	fv -- pointer to Field Variable		structure
 *  fv_old -- pointer to old Diet Field Variable	structure
 *  fv_dot -- pointer to dot Diet Field Variable	structure
 *	cr -- pointer to Constitutive Relation	structure
 *	md -- pointer to Mesh Derivative	structure
 *	me -- pointer to Material Entity	structure
 *
 * out:
 *	a   -- gets loaded up with proper contribution 
 *	lec -- gets loaded up with local contributions to resid, Jacobian
 * 	r   -- residual RHS vector
 *
 * Created:	Friday August 15, 2010 prschun@sandia.gov
 *
 */
/*ARGSUSED*/
int 
assemble_shell_deltah(double time,	/* present time value */ 
		     double tt,   	/* parameter to vary time integration from 
			   	         * explicit (tt = 1) to implicit (tt = 0)    */
		     double dt,         /* current time step size */
		     double xi[DIM],    /* Local stu coordinates */
		     const Exo_DB *exo)  
{
  EH(-1, "assemble_shell_deltah disabled. Contact prschun@sandia.gov");
  return(-1);
} /* end of assemble_shell_deltah */





/*****************************************************************************/
/***assemble_lubrication_curvature********************************************/
/*  _______________________________________________________________________  */

/* assemble_lubrication_Curvature -- Calculates the curvature of the level 
 *                                   set field for a lubricaiton field.  
 *
 * Created:     October 19, 2010 sarober@sandia.gov
 *
 */
int
assemble_lubrication_curvature(
			       double time,              /* present time value */ 
			       double tt,   	         /* parameter to vary time integration  */
			       double dt,                /* current time step size */
			       const PG_DATA *pg_data,   /* Element scales */
			       double xi[DIM],           /* Local stu coordinates */
			       const Exo_DB *exo)  {     /* Exodus database */

  /* --- Initialization -----------------------------------------------------*/    
  
  /* Variable definitions */
  int eqn = R_SHELL_LUB_CURV;                    
  int peqn, var, pvar;
  int status = 0;                                
  int i, j, k, a, jj, b;
  dbl phi_i, grad_phi_i[DIM], grad_II_phi_i[DIM], d_grad_II_phi_i_dmesh[DIM][DIM][MDE];
  dbl phi_j, grad_phi_j[DIM], grad_II_phi_j[DIM], d_grad_II_phi_j_dmesh[DIM][DIM][MDE];
  dbl mass, diff, div;

  /* Bail out fast if there's nothing to do */
  if ( !pd->e[eqn]  ) return(status);
  if ( !pd->e[FILL] ) EH(-1, "Must activate level set equation to calculate curvature.");

  /* Prepare shell geometry */
  dbl wt_old = fv->wt;
  int *n_dof = NULL;
  int dof_map[MDE];
  n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
  lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

  /* Load up FEM weightings */
  dbl wt = fv->wt;                                // Gauss point weight
  dbl h3 = fv->h3;                                // Differential volume element
  dbl det_J = fv->sdet;                           // Jacobian of transformation

  /* --- Calculate problem parameters ---------------------------------------*/

  /* Load level set fields */
  load_lsi( ls->Length_Scale );
  load_lsi_derivs();

  /* Rotate grad(F) and grad(kappa) to shell coordinates */
  dbl gradII_F[DIM], gradII_kappa[DIM];
  dbl d_grad_F_dmesh[DIM][DIM][MDE], d_gradII_F_dmesh[DIM][DIM][MDE];
  dbl d_grad_kappa_dmesh[DIM][DIM][MDE], d_gradII_kappa_dmesh[DIM][DIM][MDE];
  memset(d_grad_F_dmesh, 0.0, sizeof(double)*DIM*DIM*MDE);
  memset(d_grad_kappa_dmesh, 0.0, sizeof(double)*DIM*DIM*MDE);
  for ( a = 0; a < DIM; a++) {
    for ( b = 0; b < DIM; b++) {
      for ( i = 0; i < ei->dof[MESH_DISPLACEMENT1]; i++) {
	j = dof_map[i];
	d_grad_F_dmesh[a][b][j]     = fv->d_grad_F_dmesh[a][b][i];
	d_grad_kappa_dmesh[a][b][j] = fv->d_grad_sh_l_curv_dmesh[a][b][i];
      }
    }
  }
  ShellRotate(fv->grad_F, d_grad_F_dmesh, gradII_F, d_gradII_F_dmesh, n_dof[MESH_DISPLACEMENT1]);
  ShellRotate(fv->grad_sh_l_curv, d_grad_kappa_dmesh, gradII_kappa, d_gradII_kappa_dmesh, n_dof[MESH_DISPLACEMENT1]);

  /* Calculate rotated level set normal */
  dbl LSnormal[DIM], LSnormal_mag = 0, LSnormal_maginv;
  for ( i = 0; i < DIM; i++ ) {
    LSnormal_mag += gradII_F[i] * gradII_F[i];
  }
  LSnormal_mag = sqrt(LSnormal_mag);
  LSnormal_maginv = ( LSnormal_mag == 0.0 ) ? 1.0 : 1.0 / LSnormal_mag;
  for ( i = 0; i < DIM; i++ ) {
    LSnormal[i] = gradII_F[i] * LSnormal_maginv;
  }

  /* Calculate sensitivity of level set normal to mesh */
  dbl d_LSnormal_dmesh[DIM][DIM][MDE];
  dbl d_LSnormal_mag_dmesh[DIM][MDE];
  memset(d_LSnormal_dmesh,     0.0, sizeof(double)*DIM*DIM*MDE);
  memset(d_LSnormal_mag_dmesh, 0.0, sizeof(double)*DIM*MDE);
  for ( i = 0; i < DIM; i++ ) {
    for ( j = 0; j < DIM; j++) {
      for ( k = 0; k < n_dof[MESH_DISPLACEMENT1]; k++) {
	d_LSnormal_mag_dmesh[j][k] += gradII_F[i] * d_gradII_F_dmesh[i][j][k] * LSnormal_maginv;
      }
    }
  }
  for ( i = 0; i < DIM; i++) {
    for ( j = 0; j < DIM; j++) {
      for ( k = 0; k < n_dof[MESH_DISPLACEMENT1]; k++) {
	d_LSnormal_dmesh[i][j][k] += d_gradII_F_dmesh[i][j][k] * LSnormal_maginv;
	d_LSnormal_dmesh[i][j][k] -= gradII_F[i] * LSnormal_maginv * LSnormal_maginv * d_LSnormal_mag_dmesh[j][k];
      }
    }
  }

  /* Calculate sensitivity of level set normal to F */
  dbl d_LSnormal_mag_dF[MDE] ;
  dbl d_LSnormal_dF[DIM][MDE];
  memset(d_LSnormal_mag_dF, 0.0, sizeof(double)*MDE);
  memset(d_LSnormal_dF, 0.0, sizeof(double)*MDE*DIM);
  var = FILL;
  for ( i = 0; i < ei->dof[FILL]; i++ ) {

    /* Basis functions */
    ShellBF( eqn, i, &phi_i, grad_phi_i, grad_II_phi_i, d_grad_II_phi_i_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );

    /* Derivative of magnitude */
    for ( a = 0; a < DIM; a++ ) {
      d_LSnormal_mag_dF[i] += gradII_F[a] * grad_II_phi_i[a] * LSnormal_maginv;
    }

    /* Derivative of normal vector */
    for ( a = 0; a < DIM; a++ ) {
      d_LSnormal_dF[a][i]  = grad_II_phi_i[a] * LSnormal_maginv;
      d_LSnormal_dF[a][i] -= gradII_F[a] * d_LSnormal_mag_dF[i] * pow(LSnormal_maginv,2);
    }

  }
    
  /* Prepare weighting for artificial diffusion term */
  const dbl *hsquared = pg_data->hsquared;  

  /* --- Residual assembly --------------------------------------------------*/
  if (af->Assemble_Residual) {
    eqn = R_SHELL_LUB_CURV;
    peqn = upd->ep[eqn];

    /* Loop over DOFs (i) */
    for ( i=0; i<ei->dof[eqn]; i++) {
      
      /* Prepare basis funcitons (i) */
      ShellBF( eqn, i, &phi_i, grad_phi_i, grad_II_phi_i, d_grad_II_phi_i_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );

      /* Assemble mass term */
      mass = 0.0;
      if ( T_MASS ) {
	mass += *esp->sh_l_curv[i] * phi_i;
      }
      mass *= det_J * wt * h3 * pd->etm[eqn][(LOG2_MASS)];

      /* Assemble diffusion terms */
      diff = 0.0;
      if ( T_DIFFUSION ) {
	for ( a = 0; a < VIM; a++ ) {
	  diff += *hsquared * gradII_kappa[a] * grad_II_phi_i[a];
	}
      }
      diff *= det_J * wt * h3 * pd->etm[eqn][(LOG2_DIFFUSION)];

      /* Assemble divergence terms */
      div = 0.0;
      if ( T_DIVERGENCE ) {
	for ( a = 0; a < VIM; a++ ) {
	  div += LSnormal[a] * grad_II_phi_i[a];
	}
      }
      div *= det_J * wt * h3 * pd->etm[eqn][(LOG2_DIVERGENCE)];

      /* Assemble residual */
      lec->R[peqn][i] += mass + diff + div;

    } // End of loop over DOFs (i)
  } // End of residual assembly


  /* --- Jacobian assembly --------------------------------------------------*/
  if (af->Assemble_Jacobian) {
    eqn = R_SHELL_LUB_CURV;
    peqn = upd->ep[eqn];

    /* Loop over DOFs (i) */
    for ( i=0; i<ei->dof[eqn]; i++) {
      
      /* Prepare basis funcitons (i) */
      ShellBF( eqn, i, &phi_i, grad_phi_i, grad_II_phi_i, d_grad_II_phi_i_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );


      /*** SHELL_LUB_CURV ***/
      var = SHELL_LUB_CURV;
      if ( pd->v[var] ) {
	pvar = upd->vp[var];

	/* Loop over DOFs (j) */
	for ( j = 0; j < ei->dof[var]; j++ ) {
      
	  /* Prepare basis funcitons (j) */
	  ShellBF( var, j, &phi_j, grad_phi_j, grad_II_phi_j, d_grad_II_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
	  	  
	  /* Assemble mass term */
	  mass = 0.0;
	  if ( T_MASS ) {
	    //mass += phi_i * phi_j;
	    if ( i == j ) mass += phi_i;
	  }
	  mass *= det_J * wt * h3 * pd->etm[eqn][(LOG2_MASS)];
      
	  /* Assemble diffusion terms */
	  diff = 0.0;
	  if ( T_DIFFUSION ) {
	    for ( a = 0; a < VIM; a++ ) {
	      diff += *hsquared * grad_II_phi_i[a] * grad_II_phi_j[a];
	    }
	  }
	  diff *= det_J * wt * h3 * pd->etm[eqn][(LOG2_DIFFUSION)];
	  
	  /* Assemble jacobian */
	  lec->J[peqn][pvar][i][j] += mass + diff;

	} // End of loop over DOFs (j)
      } // End of SH_LUB_CURV assembly


      /*** MESH_DISPLACEMENT1 ***/
      var = MESH_DISPLACEMENT1;
      if ( pd->v[var]  && (
			   mp->FSIModel == FSI_MESH_CONTINUUM ||
			   mp->FSIModel == FSI_MESH_UNDEF
			   )) {
	pvar = upd->vp[var];
	
	/*** Loop over dimensions of mesh displacement ***/
	for ( b = 0; b < DIM; b++) {
	  var = MESH_DISPLACEMENT1 + b;
	  pvar = upd->vp[var];
	  
	  /* Loop over DOFs (j) */
	  for ( j = 0; j < ei->dof[var]; j++ ) {
	    jj = dof_map[j];
	    
	    /* Prepare basis funcitons (j) */
	    ShellBF( var, j, &phi_j, grad_phi_j, grad_II_phi_j, d_grad_II_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
	    
	    /* Assemble mass term */
	    mass = 0.0;
	    if ( T_MASS ) {
	      mass += *esp->sh_l_curv[i] * phi_i;
	    }
	    mass *= fv->dsurfdet_dx[b][jj] * wt * h3 * pd->etm[eqn][(LOG2_MASS)];

	    /* Assemble diffusion terms */
	    diff = 0.0;
	    if ( T_DIFFUSION ) {
	      for ( a = 0; a < VIM; a++ ) {
		diff += *hsquared * d_gradII_kappa_dmesh[a][b][jj] * grad_II_phi_i[a] * det_J;
		diff += *hsquared * gradII_kappa[a] * d_grad_II_phi_i_dmesh[a][b][jj] * det_J;
		diff += *hsquared * gradII_kappa[a] * grad_II_phi_i[a] * fv->dsurfdet_dx[b][jj];
	      }
	    }
	    diff *= wt * h3 * pd->etm[eqn][(LOG2_DIFFUSION)];
	    
	    /* Assemble divergence terms */
	    div = 0.0;
	    if ( T_DIVERGENCE ) {
	      for ( a = 0; a < VIM; a++ ) {
		div += d_LSnormal_dmesh[a][b][jj] * grad_II_phi_i[a] * det_J;
		div += LSnormal[a] * d_grad_II_phi_i_dmesh[a][b][jj] * det_J;
		div += LSnormal[a] * grad_II_phi_i[a] * fv->dsurfdet_dx[b][jj];
	      }
	    }
	    div *= wt * h3 * pd->etm[eqn][(LOG2_DIVERGENCE)];
	    
	    /* Assemble jacobian */
	    lec->J[peqn][pvar][i][jj] += mass + diff + div;
	    
	  } // End of loop over DOFs (j)
	} // End of loop over mesh dimensions
      } // End of DMX assembly
      
      
      /*** FILL ***/
      var = FILL;
      if ( pd->v[var] ) {
	pvar = upd->vp[var];
	
	/* Loop over DOFs (j) */
	for ( j = 0; j < ei->dof[var]; j++ ) {
	  
	  /* Assemble divergence terms */
	  div = 0.0;
	  if ( T_DIVERGENCE ) {
	    for ( a = 0; a < VIM; a++ ) {
	      div += d_LSnormal_dF[a][j] * grad_II_phi_i[a];
	    }
	  }
	  div *= det_J * wt * h3 * pd->etm[eqn][(LOG2_DIVERGENCE)];
	  
	  /* Assemble jacobian */
	  lec->J[peqn][pvar][i][j] += div;

	} // End of loop over DOFs (j)
      } // End of FILL assembly

    } // End of loop over DOFs (i)
  } // End of jacobian assembly
  
  /* Clean up */
  fv->wt = wt_old;
  safe_free((void *) n_dof);

  return(status);
} // End of assemble_lubrication_curvature

/*****************************************************************************/
/***assemble_lubrication_curvature********************************************/
/*  _______________________________________________________________________  */

/* assemble_lubrication_Curvature_2 -- Calculates the curvature of the phase 
 *                                   field for a lubricaiton_2 field.  
 *
 * Created:     August 30, 2012 prschun@sandia.gov
 *
 */
int
assemble_lubrication_curvature_2(
				 double time,              /* present time value */ 
				 double tt,   	         /* parameter to vary time integration  */
				 double dt,                /* current time step size */
				 const PG_DATA *pg_data,   /* Element scales */
				 double xi[DIM],           /* Local stu coordinates */
				 const Exo_DB *exo)  {     /* Exodus database */

  /* --- Initialization -----------------------------------------------------*/    
  
  /* Variable definitions */
  int eqn = R_SHELL_LUB_CURV_2;                    
  int peqn, var, pvar;
  int status = 0;                                
  int i, j, k, a, jj, b;
  dbl phi_i, grad_phi_i[DIM], grad_II_phi_i[DIM], d_grad_II_phi_i_dmesh[DIM][DIM][MDE];
  dbl phi_j, grad_phi_j[DIM], grad_II_phi_j[DIM], d_grad_II_phi_j_dmesh[DIM][DIM][MDE];
  dbl mass, diff, div;

  /* Bail out fast if there's nothing to do */
  if ( !pd->e[eqn]  ) return(status);
  if ( !pd->e[PHASE1] ) EH(-1, "Must activate phase1 equation to calculate curvature.");

  /* Prepare shell geometry */
  dbl wt_old = fv->wt;
  int *n_dof = NULL;
  int dof_map[MDE];
  n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
  lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

  /* Load up FEM weightings */
  dbl wt = fv->wt;                                // Gauss point weight
  dbl h3 = fv->h3;                                // Differential volume element
  dbl det_J = fv->sdet;                           // Jacobian of transformation

 struct Level_Set_Data *ls_old;

  /* --- Calculate problem parameters ---------------------------------------*/

  /* Load level set fields */

  ls_old = ls;
  if(pfd != NULL) ls = pfd->ls[0];   
  load_lsi( ls->Length_Scale );
  load_lsi_derivs();
  ls = ls_old;

  /* Rotate grad(F) and grad(kappa) to shell coordinates */
  if(upd->ep[R_MESH1] > -1) EH(-1," Must add mesh dependence to phase field");
  dbl gradII_F[DIM], gradII_kappa[DIM];
  dbl d_grad_F_dmesh[DIM][DIM][MDE], d_gradII_F_dmesh[DIM][DIM][MDE];
  dbl d_grad_kappa_dmesh[DIM][DIM][MDE], d_gradII_kappa_dmesh[DIM][DIM][MDE];
  memset(d_grad_F_dmesh, 0.0, sizeof(double)*DIM*DIM*MDE);
  memset(d_grad_kappa_dmesh, 0.0, sizeof(double)*DIM*DIM*MDE);
  //for ( a = 0; a < DIM; a++) {
  //  for ( b = 0; b < DIM; b++) {
  //    for ( i = 0; i < ei->dof[MESH_DISPLACEMENT1]; i++) {
  //	j = dof_map[i];
  //	d_grad_F_dmesh[a][b][j]     = fv->d_grad_F_dmesh[a][b][i];
  //	d_grad_kappa_dmesh[a][b][j] = fv->d_grad_sh_l_curv_dmesh[a][b][i];
  //    }
  //  }
  //}
  ShellRotate(fv->grad_pF[0], d_grad_F_dmesh, gradII_F, d_gradII_F_dmesh, n_dof[MESH_DISPLACEMENT1]);
  ShellRotate(fv->grad_sh_l_curv_2, d_grad_kappa_dmesh, gradII_kappa, d_gradII_kappa_dmesh, n_dof[MESH_DISPLACEMENT1]);

  /* Calculate rotated level set normal */
  dbl LSnormal[DIM], LSnormal_mag = 0, LSnormal_maginv;
  for ( i = 0; i < DIM; i++ ) {
    LSnormal_mag += gradII_F[i] * gradII_F[i];
  }
  LSnormal_mag = sqrt(LSnormal_mag);
  LSnormal_maginv = ( LSnormal_mag == 0.0 ) ? 1.0 : 1.0 / LSnormal_mag;
  for ( i = 0; i < DIM; i++ ) {
    LSnormal[i] = gradII_F[i] * LSnormal_maginv;
  }

  /* Calculate sensitivity of level set normal to mesh */
  /* Just multiplying zeros right now until you implement a d_grad_pf_dmesh field PRS-9/7/2012 */
  dbl d_LSnormal_dmesh[DIM][DIM][MDE];
  dbl d_LSnormal_mag_dmesh[DIM][MDE];
  memset(d_LSnormal_dmesh,     0.0, sizeof(double)*DIM*DIM*MDE);
  memset(d_LSnormal_mag_dmesh, 0.0, sizeof(double)*DIM*MDE);
  for ( i = 0; i < DIM; i++ ) {
    for ( j = 0; j < DIM; j++) {
      for ( k = 0; k < n_dof[MESH_DISPLACEMENT1]; k++) {
	d_LSnormal_mag_dmesh[j][k] += gradII_F[i] * d_gradII_F_dmesh[i][j][k] * LSnormal_maginv;
      }
    }
  }
  for ( i = 0; i < DIM; i++) {
    for ( j = 0; j < DIM; j++) {
      for ( k = 0; k < n_dof[MESH_DISPLACEMENT1]; k++) {
	d_LSnormal_dmesh[i][j][k] += d_gradII_F_dmesh[i][j][k] * LSnormal_maginv;
	d_LSnormal_dmesh[i][j][k] -= gradII_F[i] * LSnormal_maginv * LSnormal_maginv * d_LSnormal_mag_dmesh[j][k];
      }
    }
  }

  /* Calculate sensitivity of phase field normal to pf */
  dbl d_LSnormal_mag_dF[MDE] = {0.0};
  dbl d_LSnormal_dF[DIM][MDE] = {{0.0}};
  var = PHASE1;
  for ( i = 0; i < ei->dof[PHASE1]; i++ ) {

    /* Basis functions */
    ShellBF( eqn, i, &phi_i, grad_phi_i, grad_II_phi_i, d_grad_II_phi_i_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );

    /* Derivative of magnitude */
    for ( a = 0; a < DIM; a++ ) {
      d_LSnormal_mag_dF[i] += gradII_F[a] * grad_II_phi_i[a] * LSnormal_maginv;
    }

    /* Derivative of normal vector */
    for ( a = 0; a < DIM; a++ ) {
      d_LSnormal_dF[a][i]  = grad_II_phi_i[a] * LSnormal_maginv;
      d_LSnormal_dF[a][i] -= gradII_F[a] * d_LSnormal_mag_dF[i] * pow(LSnormal_maginv,2);
    }

  }
    
  /* Prepare weighting for artificial diffusion term */
  const dbl *hsquared = pg_data->hsquared;  

  /* --- Residual assembly --------------------------------------------------*/
  if (af->Assemble_Residual) {
    eqn = R_SHELL_LUB_CURV_2;
    peqn = upd->ep[eqn];

    /* Loop over DOFs (i) */
    for ( i=0; i<ei->dof[eqn]; i++) {
      
      /* Prepare basis funcitons (i) */
      ShellBF( eqn, i, &phi_i, grad_phi_i, grad_II_phi_i, d_grad_II_phi_i_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );

      /* Assemble mass term */
      mass = 0.0;
      if ( T_MASS ) {
	mass += *esp->sh_l_curv_2[i] * phi_i;
      }
      mass *= det_J * wt * h3 * pd->etm[eqn][(LOG2_MASS)];

      /* Assemble diffusion terms */
      diff = 0.0;
      if ( T_DIFFUSION ) {
	for ( a = 0; a < VIM; a++ ) {
	  diff += *hsquared * gradII_kappa[a] * grad_II_phi_i[a];
	}
      }
      diff *= det_J * wt * h3 * pd->etm[eqn][(LOG2_DIFFUSION)];

      /* Assemble divergence terms */
      div = 0.0;
      if ( T_DIVERGENCE ) {
	for ( a = 0; a < VIM; a++ ) {
	  div += LSnormal[a] * grad_II_phi_i[a];
	}
      }
      div *= det_J * wt * h3 * pd->etm[eqn][(LOG2_DIVERGENCE)];

      /* Assemble residual */
      lec->R[peqn][i] += mass + diff + div;

    } // End of loop over DOFs (i)
  } // End of residual assembly


  /* --- Jacobian assembly --------------------------------------------------*/
  if (af->Assemble_Jacobian) {
    eqn = R_SHELL_LUB_CURV_2;
    peqn = upd->ep[eqn];

    /* Loop over DOFs (i) */
    for ( i=0; i<ei->dof[eqn]; i++) {
      
      /* Prepare basis funcitons (i) */
      ShellBF( eqn, i, &phi_i, grad_phi_i, grad_II_phi_i, d_grad_II_phi_i_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );


      /*** SHELL_LUB_CURV ***/
      var = SHELL_LUB_CURV_2;
      if ( pd->v[var] ) {
	pvar = upd->vp[var];

	/* Loop over DOFs (j) */
	for ( j = 0; j < ei->dof[var]; j++ ) {
      
	  /* Prepare basis funcitons (j) */
	  ShellBF( var, j, &phi_j, grad_phi_j, grad_II_phi_j, d_grad_II_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
	  	  
	  /* Assemble mass term */
	  mass = 0.0;
	  if ( T_MASS ) {
	    //mass += phi_i * phi_j;
	    if ( i == j ) mass += phi_i;
	  }
	  mass *= det_J * wt * h3 * pd->etm[eqn][(LOG2_MASS)];
      
	  /* Assemble diffusion terms */
	  diff = 0.0;
	  if ( T_DIFFUSION ) {
	    for ( a = 0; a < VIM; a++ ) {
	      diff += *hsquared * grad_II_phi_i[a] * grad_II_phi_j[a];
	    }
	  }
	  diff *= det_J * wt * h3 * pd->etm[eqn][(LOG2_DIFFUSION)];
	  
	  /* Assemble jacobian */
	  lec->J[peqn][pvar][i][j] += mass + diff;

	} // End of loop over DOFs (j)
      } // End of SH_LUB_CURV assembly


      /*** MESH_DISPLACEMENT1 ***/
      var = MESH_DISPLACEMENT1;
      if ( pd->v[var]  && (
			   mp->FSIModel == FSI_MESH_CONTINUUM ||
			   mp->FSIModel == FSI_MESH_UNDEF
			   )) {
	pvar = upd->vp[var];
	
	/*** Loop over dimensions of mesh displacement ***/
	for ( b = 0; b < DIM; b++) {
	  var = MESH_DISPLACEMENT1 + b;
	  pvar = upd->vp[var];
	  
	  /* Loop over DOFs (j) */
	  for ( j = 0; j < ei->dof[var]; j++ ) {
	    jj = dof_map[j];
	    
	    /* Prepare basis funcitons (j) */
	    ShellBF( var, j, &phi_j, grad_phi_j, grad_II_phi_j, d_grad_II_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map );
	    
	    /* Assemble mass term */
	    mass = 0.0;
	    if ( T_MASS ) {
	      mass += *esp->sh_l_curv_2[i] * phi_i;
	    }
	    mass *= fv->dsurfdet_dx[b][jj] * wt * h3 * pd->etm[eqn][(LOG2_MASS)];

	    /* Assemble diffusion terms */
	    diff = 0.0;
	    if ( T_DIFFUSION ) {
	      for ( a = 0; a < VIM; a++ ) {
		diff += *hsquared * d_gradII_kappa_dmesh[a][b][jj] * grad_II_phi_i[a] * det_J;
		diff += *hsquared * gradII_kappa[a] * d_grad_II_phi_i_dmesh[a][b][jj] * det_J;
		diff += *hsquared * gradII_kappa[a] * grad_II_phi_i[a] * fv->dsurfdet_dx[b][jj];
	      }
	    }
	    diff *= wt * h3 * pd->etm[eqn][(LOG2_DIFFUSION)];
	    
	    /* Assemble divergence terms */
	    div = 0.0;
	    if ( T_DIVERGENCE ) {
	      for ( a = 0; a < VIM; a++ ) {
		div += d_LSnormal_dmesh[a][b][jj] * grad_II_phi_i[a] * det_J;
		div += LSnormal[a] * d_grad_II_phi_i_dmesh[a][b][jj] * det_J;
		div += LSnormal[a] * grad_II_phi_i[a] * fv->dsurfdet_dx[b][jj];
	      }
	    }
	    div *= wt * h3 * pd->etm[eqn][(LOG2_DIVERGENCE)];
	    
	    /* Assemble jacobian */
	    lec->J[peqn][pvar][i][jj] += mass + diff + div;
	    
	  } // End of loop over DOFs (j)
	} // End of loop over mesh dimensions
      } // End of DMX assembly
      
      
      /*** PHASE1 ***/
      var = PHASE1;
      if ( pd->v[var] ) {
	pvar = upd->vp[var];
	
	/* Loop over DOFs (j) */
	for ( j = 0; j < ei->dof[var]; j++ ) {
	  
	  /* Assemble divergence terms */
	  div = 0.0;
	  if ( T_DIVERGENCE ) {
	    for ( a = 0; a < VIM; a++ ) {
	      div += d_LSnormal_dF[a][j] * grad_II_phi_i[a];
	    }
	  }
	  div *= det_J * wt * h3 * pd->etm[eqn][(LOG2_DIVERGENCE)];
	  
	  /* Assemble jacobian */
	  lec->J[peqn][pvar][i][j] += div;

	} // End of loop over DOFs (j)
      } // End of FILL assembly

    } // End of loop over DOFs (i)
  } // End of jacobian assembly
  
  /* Clean up */
  fv->wt = wt_old;
  safe_free((void *) n_dof);

  return(status);
} // End of assemble_lubrication_curvature

int
assemble_lubrication_power_law( double time,    /* present time value */
                                double tt,      /* parameter to vary time integration from
                                                   explicit (tt = 1) to implicit (tt = 0)    */
                                double dt,      /* current time step size */  
                                double xi[DIM], /* Local stu coordinates */
                                const Exo_DB *exo) 

{
  int eqn;
  int var, peqn, pvar, dim;
  int i, ii;
  int j, jj, status;

  dbl wt;

  dbl H, H_U, dH_U_dtime, H_L, dH_L_dtime, dH_dtime, dH_U_ddh;
  dbl dH_U_dX[DIM],dH_L_dX[DIM], dH_U_dp;
  dbl shear_top, shear_bot, cross_shear, gradP_mag;
  dbl *grad_P, *grad_II_P;
  grad_P = (double *) malloc(DIM * sizeof(double));
  grad_II_P = (double *) malloc(DIM * sizeof(double));
  dbl mu, mu0, nexp;
  dbl veloU[DIM], veloL[DIM];

  dbl advection, diffusion, source;

  /*
   * Galerkin weighting functions for i-th shell residuals
   * and some of their derivatives...
   */
  dbl phi_i;
  dbl grad_phi_i[DIM];
  dbl *grad_II_phi_i;
  grad_II_phi_i = (double *) malloc(DIM * sizeof(double));

  /*
   * Interpolation functions for variables and some of their derivatives.
   */
  dbl phi_j;
  dbl grad_phi_j[DIM];
  dbl *grad_II_phi_j;
  grad_II_phi_j = (double *) malloc(DIM * sizeof(double));


  dbl h3;                       /* Volume element (scale factors). */
  dbl det_J;

  dbl shear_top_plus, shear_top_minus, shear_bot_plus, shear_bot_minus;
  dbl source_plus, source_minus;
  dbl eps_top, eps_bot;
  int iii;

  /*   static char yo[] = "assemble_lubrication_power_law";*/

  status = 0;

  /* Unpack variables from structures for local convenience... */
  dim = pd->Num_Dim;
  shear_top = fv->sh_shear_top;  /* Top wall shear rate */
  shear_bot = fv->sh_shear_bot;  /* Bottom wall shear rate */
  cross_shear = fv->sh_cross_shear; /* Cross stream shear stress */
  grad_P = fv->grad_lubp; /* Pressure gradient */


  eps_top = 2.0e-6 * fabs(shear_top);
  eps_bot = 2.0e-6 * fabs(shear_bot);

  wt = fv->wt;       /* Gauss weight */
  h3 = fv->h3;       /* Differential volume element, = 1 when CARTESIAN. */

  /* Tangent vectors to the lubrication plane */
  dbl tangent1_init[DIM], tangent2_init[DIM];
  dbl *tangent1, *tangent2;
  tangent1 = (double *) malloc(DIM * sizeof(double));
  tangent2 = (double *) malloc(DIM * sizeof(double));
  dbl tangent1_mag, tangent2_mag;

  shell_determinant_and_normal(ei->ielem, ei->iconnect_ptr, ei->num_local_nodes,
                               ei->ielem_dim, 1);
  det_J = fv->sdet;

  Inn(grad_P, grad_II_P);

  gradP_mag = 0.0;

  for (ii = 0; ii < DIM; ii++)
     {
      gradP_mag += grad_II_P[ii] * grad_II_P[ii];
     }

  gradP_mag = sqrt(gradP_mag);


  /*
   * Material property constants, etc. Any variations at this Gauss point
   * are calculated with user-defined subroutine options.
   *
   */

  H = height_function_model(&H_U, &dH_U_dtime, &H_L, &dH_L_dtime, dH_U_dX, dH_L_dX, &dH_U_dp, &dH_U_ddh, time, dt);
  dH_dtime = dH_U_dtime - dH_L_dtime;

   mu0 = gn->mu0;
   mu = mu0;
   nexp = gn->nexp;


  velocity_function_model(veloU, veloL, time, dt);


  /* Calculate flow rate */

  calculate_lub_q_v_nonnewtonian(time, dt);


 /*Evaluate tangent vectors of the lubrication planes */

 tangent1_init[0] = 1./sqrt(3.);
 tangent1_init[1] = 1./sqrt(3.);
 tangent1_init[2] = 1./sqrt(3.);

 tangent2_init[0] = -1./sqrt(3.);
 tangent2_init[1] = 1./sqrt(3.);

 Inn(tangent1_init, tangent1);
 Inn(tangent2_init, tangent2);

 tangent1_mag = 0.0;
 tangent2_mag = 0.0;
 for (i = 0; i < DIM; i++)
    {
      tangent1_mag += tangent1[i] * tangent1[i];
      tangent2_mag += tangent2[i] * tangent2[i];
    }

 tangent1_mag = sqrt(tangent1_mag);
 tangent2_mag = sqrt(tangent2_mag);

 for (i = 0; i < DIM; i++)
    {
      tangent1[i] *= 1.0/tangent1_mag;
      tangent2[i] *= 1.0/tangent2_mag;
    }

   /*
   * Residuals_________________________________________________________________
   */

   if (af->Assemble_Residual)
    {


/* ************** ASSEMBLE RESIDUAL OF TOP WALL SHEAR RATE ********* */

      eqn = R_SHELL_SHEAR_TOP;
      peqn = upd->ep[eqn];

      for ( i=0; i<ei->dof[eqn]; i++)
        {
          phi_i = bf[eqn]->phi[i];


/* Assemble advection term */

          advection = 0.0;
          if (pd->e[eqn] & T_ADVECTION)
            {
             for (ii = 0; ii < DIM; ii++)
                {
                 advection += phi_i * tangent1[ii] * (veloU[ii] - veloL[ii]) * (gradP_mag);
                }
             advection *= det_J * wt * h3 * pd->etm[eqn][(LOG2_ADVECTION)];
            }

/* Assemble source term */

          source = 0.0;
          if (pd->e[eqn] & T_SOURCE)
            {

             for (ii = 0; ii < DIM; ii++)
                {

                 source += - mu0 * tangent1[ii] * LubAux->gradP_tangent[ii] * nexp/(nexp + 1.) *
                           (  pow( fabs(shear_top), nexp - 1.) * shear_top * shear_top
                            - pow( fabs(shear_bot), nexp - 1.) * shear_bot * shear_bot );
                 source += - cross_shear * tangent1[ii] * LubAux->gradP_normal[ii] * nexp *
                           (shear_top - shear_bot);
                }

             source *= phi_i;
             source *= det_J * wt * h3 * pd->etm[eqn][(LOG2_SOURCE)];
            }

/* Combine them all */

          lec->R[peqn][i] += advection + source;

        } /* end of loop over i */

/* ************** ASSEMBLE RESIDUAL OF BOTTOM WALL SHEAR RATE ********* */

      eqn = R_SHELL_SHEAR_BOT;
      peqn = upd->ep[eqn];

      for ( i=0; i<ei->dof[eqn]; i++)
        {
          phi_i = bf[eqn]->phi[i];


/* Assemble advection term */


          advection = 0.0;
          if (pd->e[eqn] & T_ADVECTION)
            {
             for ( ii = 0; ii < DIM; ii++)
                {
                 advection += phi_i * tangent2[ii] * (veloU[ii] - veloL[ii]) * (gradP_mag);
                }
             advection *= det_J * wt * h3 * pd->etm[eqn][(LOG2_ADVECTION)];

            }


/* Assemble source term */


          source = 0.0;
          if (pd->e[eqn] & T_SOURCE)
            {

             for (ii = 0; ii < DIM; ii++)
                {

                 source += - mu0 * tangent2[ii] * LubAux->gradP_tangent[ii] * nexp/(nexp + 1.) *
                           (  pow( fabs(shear_top), nexp - 1.) * shear_top * shear_top
                            - pow( fabs(shear_bot), nexp - 1.) * shear_bot * shear_bot );
                 source += - cross_shear * tangent2[ii] * LubAux->gradP_normal[ii] * nexp *
                           (shear_top - shear_bot);
                }

             source *= phi_i;
             source *= det_J * wt * h3 * pd->etm[eqn][(LOG2_SOURCE)];
            }

/* Combine them all */

          lec->R[peqn][i] += advection + source;

        } /* end of loop over i */

/* ************** ASSEMBLE RESIDUAL OF CROSS STREAM SHEAR STRESS ********* */

      eqn = R_SHELL_CROSS_SHEAR;
      peqn = upd->ep[eqn];

      for ( i=0; i<ei->dof[eqn]; i++)
        {
          phi_i = bf[eqn]->phi[i];


/* Assemble advection term */

          advection = 0.0;
          if (pd->e[eqn] & T_ADVECTION)
            {

             advection += phi_i * H * (gradP_mag);
             advection *= det_J * wt * h3 * pd->etm[eqn][(LOG2_ADVECTION)];

            }


/* Assemble source term */

          source = 0.0;
          if (pd->e[eqn] & T_SOURCE)
            {

             source += - phi_i * mu0 *
                       (  pow( fabs(shear_top), nexp - 1.) * shear_top
                        - pow( fabs(shear_bot), nexp - 1.) * shear_bot );
             source *= det_J * wt * h3 * pd->etm[eqn][(LOG2_SOURCE)];

            }


/* Combine them all */

          lec->R[peqn][i] += advection + source;

        } /* end of loop over i */

/* ************** ASSEMBLE RESIDUAL OF PRESSURE ********* */

      eqn = R_LUBP;
      peqn = upd->ep[eqn];

      for ( i=0; i<ei->dof[eqn]; i++)
        {
          phi_i = bf[eqn]->phi[i];

          for (ii = 0; ii < dim; ii++)
             {
              grad_phi_i[ii] = bf[eqn]->grad_phi[i][ii];
              grad_II_phi_i[ii] = 0.0;
             }

          /* Perform I-nn grad operation */
          Inn(grad_phi_i, grad_II_phi_i);


/* Assemble diffusion term */

          diffusion = 0.0;
          if (pd->e[eqn] & T_DIFFUSION)
          {
            for (ii = 0; ii < dim; ii++)
              {
                diffusion += - LubAux->q[ii] * grad_II_phi_i[ii];
              }

            diffusion *= det_J * wt * h3 * pd->etm[eqn][(LOG2_DIFFUSION)];
          }


/* Assemble source term */

          source = 0.0;
          if (pd->e[eqn] & T_SOURCE)
          {

            source += dH_dtime;

            source *= phi_i;
            source *= det_J * wt * h3 * pd->etm[eqn][(LOG2_SOURCE)];
          }


/* Combine them all */

          lec->R[peqn][i] += diffusion + source;

        } /* end of loop over i */

    } /* end of assemble_residual */



  if (af->Assemble_Jacobian)
    {


/* ************* ASSEMBLE JACOBIAN OF TOP WALL SHEAR RATE ************ */

      eqn   = R_SHELL_SHEAR_TOP;
      peqn = upd->ep[eqn];

      for ( i=0; i<ei->dof[eqn]; i++)
        {

          phi_i = bf[eqn]->phi[i];


          /* SENSITIVITY W.R.T. TOP WALL SHEAR RATE */

          var = SHELL_SHEAR_TOP;

          if ( pd->v[var] )
            {
              pvar = upd->vp[var];
              for ( j=0; j<ei->dof[var]; j++)
                {
                  phi_j = bf[var]->phi[j];

                  shear_top_plus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_top_plus += (*(esp->sh_shear_top[iii]) + 0.5 * eps_top) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_top_plus += *(esp->sh_shear_top[iii]) * bf[var]->phi[iii];
                        }
                     }

                  shear_top_minus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_top_minus += (*(esp->sh_shear_top[iii]) - 0.5 * eps_top) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_top_minus += *(esp->sh_shear_top[iii]) * bf[var]->phi[iii];
                        }
                     }

                  source_plus = 0.0;
                  if (pd->e[eqn] & T_SOURCE)
                    {
                     for (ii = 0; ii < DIM; ii++)
                        {

                         source_plus += - mu0 * tangent1[ii] * LubAux->gradP_tangent[ii] * nexp/(nexp + 1.) *
                                        (  pow( fabs(shear_top_plus), nexp - 1.) * shear_top_plus * shear_top_plus
                                        - pow( fabs(shear_bot), nexp - 1.) * shear_bot * shear_bot );
                         source_plus += - cross_shear * tangent1[ii] * LubAux->gradP_normal[ii] * nexp *
                                        (shear_top_plus - shear_bot);
                        }
                      source_plus *= phi_i* det_J * h3 * wt;
                      source_plus *= pd->etm[eqn][(LOG2_SOURCE)];
                    }


                  source_minus = 0.0;
                  if (pd->e[eqn] & T_SOURCE)
                    {
                     for (ii = 0; ii < DIM; ii++)
                        {

                         source_minus += - mu0 * tangent1[ii] * LubAux->gradP_tangent[ii] * nexp/(nexp + 1.) *
                                         (  pow( fabs(shear_top_minus), nexp - 1.) * shear_top_minus * shear_top_minus
                                         - pow( fabs(shear_bot), nexp - 1.) * shear_bot * shear_bot );
                         source_minus += - cross_shear * tangent1[ii] * LubAux->gradP_normal[ii] * nexp *
                                         (shear_top_minus - shear_bot);
                        }
                      source_minus *= phi_i* det_J * h3 * wt;
                      source_minus *= pd->etm[eqn][(LOG2_SOURCE)];
                    }

                  lec->J[peqn][pvar][i][j] += (source_plus - source_minus)/(eps_top);

                }
            }

          /* SENSITIVITY W.R.T. BOTTOM WALL SHEAR RATE */


          var = SHELL_SHEAR_BOT;

          if ( pd->v[var] )
            {
              pvar = upd->vp[var];
              for ( j=0; j<ei->dof[var]; j++)
                {
                  phi_j = bf[var]->phi[j];


                  shear_bot_plus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_bot_plus += (*(esp->sh_shear_bot[iii]) + 0.5 * eps_bot) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_bot_plus += *(esp->sh_shear_bot[iii]) * bf[var]->phi[iii];
                        }
                     }

                  shear_bot_minus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_bot_minus += (*(esp->sh_shear_bot[iii]) - 0.5 * eps_bot) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_bot_minus += *(esp->sh_shear_bot[iii]) * bf[var]->phi[iii];
                        }
                     }

                  source_plus = 0.0;
                  if (pd->e[eqn] & T_SOURCE)
                    {
                     for (ii = 0; ii < DIM; ii++)
                        {

                         source_plus += - mu0 * tangent1[ii] * LubAux->gradP_tangent[ii] * nexp/(nexp + 1.) *
                                        (  pow( fabs(shear_top), nexp - 1.) * shear_top * shear_top
                                         - pow( fabs(shear_bot_plus), nexp - 1.) * shear_bot_plus * shear_bot_plus );
                         source_plus += - cross_shear * tangent1[ii] * LubAux->gradP_normal[ii] * nexp *
                                        (shear_top - shear_bot_plus);
                        }
                      source_plus *= phi_i* det_J * h3 * wt;
                      source_plus *= pd->etm[eqn][(LOG2_SOURCE)];
                    }

                  source_minus = 0.0;
                  if (pd->e[eqn] & T_SOURCE)
                    {
                     for (ii = 0; ii < DIM; ii++)
                        {

                         source_minus += - mu0 * tangent1[ii] * LubAux->gradP_tangent[ii] * nexp/(nexp + 1.) *
                                         (  pow( fabs(shear_top), nexp - 1.) * shear_top * shear_top
                                          - pow( fabs(shear_bot_minus), nexp - 1.) * shear_bot_minus * shear_bot_minus );
                         source_minus += - cross_shear * tangent1[ii] * LubAux->gradP_normal[ii] * nexp *
                                         (shear_top - shear_bot_minus);
                        }
                      source_minus *= phi_i* det_J * h3 * wt;
                      source_minus *= pd->etm[eqn][(LOG2_SOURCE)];
                    }

                  lec->J[peqn][pvar][i][j] += (source_plus - source_minus)/(eps_bot);

                }
            }

          /* SENSITIVITY W.R.T. CROSS STREAM SHEAR STRESS */


          var = SHELL_CROSS_SHEAR;

          if (pd->v[var])
            {
              pvar = upd->vp[var];
              for ( j=0; j<ei->dof[var]; j++)
                {

                  phi_j = bf[var]->phi[j];

                  source = 0.;
                  if (pd->e[eqn] & T_SOURCE)
                    {

                      for (ii = 0; ii < DIM; ii++)
                         {
                          source +=  - phi_j * tangent1[ii] * LubAux->gradP_normal[ii] * nexp *
                                       (shear_top - shear_bot);
                         }

                      source *= phi_i * det_J * h3 * wt;
                      source *= pd->etm[eqn][(LOG2_SOURCE)];

                    }
                  lec->J[peqn][pvar][i][j] += source;
                }
            }

          /* SENSITIVITY W.R.T. PRESSURE */


          var = LUBP;

          if (pd->v[var])
            {
              pvar = upd->vp[var];
              for ( j=0; j<ei->dof[var]; j++)
                {

                  phi_j = bf[var]->phi[j];


                  for (jj = 0; jj < dim; jj++)
                     {
                      grad_phi_j[jj] = bf[var]->grad_phi[j][jj];
                      grad_II_phi_j[jj] = 0.0;
                     }

                  /* Perform I-nn grad operation */
                  Inn(grad_phi_j, grad_II_phi_j);

                  shear_top_plus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_top_plus += (*(esp->sh_shear_top[iii]) + 0.5 * eps_top) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_top_plus += *(esp->sh_shear_top[iii]) * bf[var]->phi[iii];
                        }
                     }

                  shear_top_minus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_top_minus += (*(esp->sh_shear_top[iii]) - 0.5 * eps_top) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_top_minus += *(esp->sh_shear_top[iii]) * bf[var]->phi[iii];
                        }
                     }
                  shear_bot_plus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_bot_plus += (*(esp->sh_shear_bot[iii]) + 0.5 * eps_bot) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_bot_plus += *(esp->sh_shear_bot[iii]) * bf[var]->phi[iii];
                        }
                     }

                  shear_bot_minus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_bot_minus += (*(esp->sh_shear_bot[iii]) - 0.5 * eps_bot) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_bot_minus += *(esp->sh_shear_bot[iii]) * bf[var]->phi[iii];
                        }
                     }
                  calculate_lub_q_v_nonnewtonian_sens ( mu,H, veloU, veloL, grad_II_P, phi_j, grad_II_phi_j,
                                                        shear_top_plus, shear_top_minus, shear_bot_plus, shear_bot_minus,
                                                        eps_top, eps_bot );


                  advection = 0.;
                  if ( pd->e[eqn] & T_ADVECTION )
                    {
                     for (ii = 0; ii < DIM; ii++)
                        {
                         advection += tangent1[ii] * (veloU[ii] - veloL[ii]) * LubAux->dgradP_mag_dP;
                        }
                     advection *= phi_i * det_J * h3 * wt;
                     advection *= pd->etm[eqn][(LOG2_ADVECTION)];
                    }


                  source = 0.;
                  if ( pd->e[eqn] & T_SOURCE )
                    {
                     for (ii = 0; ii < DIM; ii++)
                        {

                         source += - mu0 * tangent1[ii] * LubAux->dgradP_tangent_dP[ii] * nexp/(nexp + 1.) *
                                   (  pow( fabs(shear_top), nexp - 1.) * shear_top * shear_top
                                    - pow( fabs(shear_bot), nexp - 1.) * shear_bot * shear_bot );
                         source += - cross_shear * tangent1[ii] * LubAux->dgradP_normal_dP[ii] * nexp *
                                     (shear_top - shear_bot);

                        }

                     source *= phi_i * det_J * h3 * wt;
                     source *= pd->etm[eqn][(LOG2_SOURCE)];
                    }
                  lec->J[peqn][pvar][i][j] += advection + source;
                }
            }
        }                        

/* ************* ASSEMBLE JACOBIAN OF BOTTOM WALL SHEAR RATE ************ */

      eqn   = R_SHELL_SHEAR_BOT;
      peqn = upd->ep[eqn];

      for ( i=0; i<ei->dof[eqn]; i++)
        {

          phi_i = bf[eqn]->phi[i];

          /* SENSITIVITY W.R.T. TOP WALL SHEAR RATE */

          var = SHELL_SHEAR_TOP;

          if ( pd->v[var] )
            {
              pvar = upd->vp[var];
              for ( j=0; j<ei->dof[var]; j++)
                {
                  phi_j = bf[var]->phi[j];


                  shear_top_plus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_top_plus += (*(esp->sh_shear_top[iii]) + 0.5 * eps_top) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_top_plus += *(esp->sh_shear_top[iii]) * bf[var]->phi[iii];
                        }
                     }

                  shear_top_minus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_top_minus += (*(esp->sh_shear_top[iii]) - 0.5 * eps_top) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_top_minus += *(esp->sh_shear_top[iii]) * bf[var]->phi[iii];
                        }
                     }

                  source_plus = 0.0;
                  if (pd->e[eqn] & T_SOURCE)
                    {
                     for (ii = 0; ii < DIM; ii++)
                        {

                         source_plus += - mu0 * tangent2[ii] * LubAux->gradP_tangent[ii] * nexp/(nexp + 1.) *
                                        (  pow( fabs(shear_top_plus), nexp - 1.) * shear_top_plus * shear_top_plus
                                        - pow( fabs(shear_bot), nexp - 1.) * shear_bot * shear_bot );
                         source_plus += - cross_shear * tangent2[ii] * LubAux->gradP_normal[ii] * nexp *
                                        (shear_top_plus - shear_bot);
                        }
                      source_plus *= phi_i* det_J * h3 * wt;
                      source_plus *= pd->etm[eqn][(LOG2_SOURCE)];
                    }

                  source_minus = 0.0;
                  if (pd->e[eqn] & T_SOURCE)
                    {
                     for (ii = 0; ii < DIM; ii++)
                        {

                         source_minus += - mu0 * tangent2[ii] * LubAux->gradP_tangent[ii] * nexp/(nexp + 1.) *
                                         (  pow( fabs(shear_top_minus), nexp - 1.) * shear_top_minus * shear_top_minus
                                         - pow( fabs(shear_bot), nexp - 1.) * shear_bot * shear_bot );
                         source_minus += - cross_shear * tangent2[ii] * LubAux->gradP_normal[ii] * nexp *
                                         (shear_top_minus - shear_bot);
                        }
                      source_minus *= phi_i* det_J * h3 * wt;
                      source_minus *= pd->etm[eqn][(LOG2_SOURCE)];
                    }

                  lec->J[peqn][pvar][i][j] += (source_plus - source_minus)/(eps_top);

                }
            }

          /* SENSITIVITY W.R.T. BOTTOM WALL SHEAR RATE */


          var = SHELL_SHEAR_BOT;

          if (pd->v[var])
            {
              pvar = upd->vp[var];
              for ( j=0; j<ei->dof[var]; j++)
                {

                  phi_j = bf[var]->phi[j];

                  shear_bot_plus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_bot_plus += (*(esp->sh_shear_bot[iii]) + 0.5 * eps_bot) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_bot_plus += *(esp->sh_shear_bot[iii]) * bf[var]->phi[iii];
                        }
                     }

                  shear_bot_minus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_bot_minus += (*(esp->sh_shear_bot[iii]) - 0.5 * eps_bot) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_bot_minus += *(esp->sh_shear_bot[iii]) * bf[var]->phi[iii];
                        }
                     }

                  source_plus = 0.0;
                  if (pd->e[eqn] & T_SOURCE)
                    {
                     for (ii = 0; ii < DIM; ii++)
                        {

                         source_plus += - mu0 * tangent2[ii] * LubAux->gradP_tangent[ii] * nexp/(nexp + 1.) *
                                        (  pow( fabs(shear_top), nexp - 1.) * shear_top * shear_top
                                         - pow( fabs(shear_bot_plus), nexp - 1.) * shear_bot_plus * shear_bot_plus );
                         source_plus += - cross_shear * tangent2[ii] * LubAux->gradP_normal[ii] * nexp *
                                        (shear_top - shear_bot_plus);
                        }
                      source_plus *= phi_i* det_J * h3 * wt;
                      source_plus *= pd->etm[eqn][(LOG2_SOURCE)];
                    }


                  source_minus = 0.0;
                  if (pd->e[eqn] & T_SOURCE)
                    {
                     for (ii = 0; ii < DIM; ii++)
                        {

                         source_minus += - mu0 * tangent2[ii] * LubAux->gradP_tangent[ii] * nexp/(nexp + 1.) *
                                         (  pow( fabs(shear_top), nexp - 1.) * shear_top * shear_top
                                          - pow( fabs(shear_bot_minus), nexp - 1.) * shear_bot_minus * shear_bot_minus );
                         source_minus += - cross_shear * tangent2[ii] * LubAux->gradP_normal[ii] * nexp *
                                         (shear_top - shear_bot_minus);
                        }
                      source_minus *= phi_i* det_J * h3 * wt;
                      source_minus *= pd->etm[eqn][(LOG2_SOURCE)];
                    }

                  lec->J[peqn][pvar][i][j] += (source_plus - source_minus)/(eps_bot);


                }
            }

          /* SENSITIVITY W.R.T. CROSS STREAM SHEAR STRESS */


          var = SHELL_CROSS_SHEAR;

          if (pd->v[var])
            {
              pvar = upd->vp[var];
              for ( j=0; j<ei->dof[var]; j++)
                {

                  phi_j = bf[var]->phi[j];

                  source = 0.;
                  if (pd->e[eqn] & T_SOURCE)
                    {

                      for (ii = 0; ii < DIM; ii++)
                         {
                          source +=  - phi_j * tangent2[ii] * LubAux->gradP_normal[ii] * nexp *
                                       (shear_top - shear_bot);
                         }
                      source *= phi_i * det_J * h3 * wt;
                      source *= pd->etm[eqn][(LOG2_SOURCE)];

                    }

                  lec->J[peqn][pvar][i][j] += source;
                }
            }

          /* SENSITIVITY W.R.T. PRESSURE */


          var = LUBP;

          if (pd->v[var])
            {
              pvar = upd->vp[var];
              for ( j=0; j<ei->dof[var]; j++)
                {

                  phi_j = bf[var]->phi[j];


                  for (jj = 0; jj < dim; jj++)
                     {
                      grad_phi_j[jj] = bf[var]->grad_phi[j][jj];
                      grad_II_phi_j[jj] = 0.0;
                     }

                  /* Perform I-nn grad operation */
                  Inn(grad_phi_j, grad_II_phi_j);

                  shear_top_plus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_top_plus += (*(esp->sh_shear_top[iii]) + 0.5 * eps_top) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_top_plus += *(esp->sh_shear_top[iii]) * bf[var]->phi[iii];
                        }
                     }

                  shear_top_minus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_top_minus += (*(esp->sh_shear_top[iii]) - 0.5 * eps_top) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_top_minus += *(esp->sh_shear_top[iii]) * bf[var]->phi[iii];
                        }
                     }

                  shear_bot_plus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_bot_plus += (*(esp->sh_shear_bot[iii]) + 0.5 * eps_bot) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_bot_plus += *(esp->sh_shear_bot[iii]) * bf[var]->phi[iii];
                        }
                     }

                  shear_bot_minus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_bot_minus += (*(esp->sh_shear_bot[iii]) - 0.5 * eps_bot) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_bot_minus += *(esp->sh_shear_bot[iii]) * bf[var]->phi[iii];
                        }
                     }


                  calculate_lub_q_v_nonnewtonian_sens ( mu,H, veloU, veloL, grad_II_P, phi_j, grad_II_phi_j,
                                                        shear_top_plus, shear_top_minus, shear_bot_plus, shear_bot_minus,
                                                        eps_top, eps_bot );

                  advection = 0.;
                  if ( pd->e[eqn] & T_ADVECTION )
                    {
                     for (ii = 0; ii < DIM; ii++)
                        {
                         advection += tangent2[ii] * (veloU[ii] - veloL[ii]) * LubAux->dgradP_mag_dP;
                        }
                     advection *= phi_i * det_J * h3 * wt;
                     advection *= pd->etm[eqn][(LOG2_ADVECTION)];
                    }

                  source = 0.;
                  if ( pd->e[eqn] & T_SOURCE )
                    {
                     for (ii = 0; ii < DIM; ii++)
                        {

                         source += - mu0 * tangent2[ii] * LubAux->dgradP_tangent_dP[ii] * nexp/(nexp + 1.) *
                                   (  pow( fabs(shear_top), nexp - 1.) * shear_top * shear_top
                                    - pow( fabs(shear_bot), nexp - 1.) * shear_bot * shear_bot );
                         source += - cross_shear * tangent2[ii] * LubAux->dgradP_normal_dP[ii] * nexp *
                                     (shear_top - shear_bot);

                        }

                     source *= phi_i * det_J * h3 * wt;
                     source *= pd->etm[eqn][(LOG2_SOURCE)];
                    }

                  lec->J[peqn][pvar][i][j] += advection + source;
                }
            }
        }


/* ************* ASSEMBLE JACOBIAN OF CROSS STREAM SHEAR STRESS ************ */

      eqn   = R_SHELL_CROSS_SHEAR;
      peqn = upd->ep[eqn];

      for ( i=0; i<ei->dof[eqn]; i++)
        {

          phi_i = bf[eqn]->phi[i];


          /* SENSITIVITY W.R.T. TOP WALL SHEAR RATE */

          var = SHELL_SHEAR_TOP;

          if ( pd->v[var] )
            {
              pvar = upd->vp[var];
              for ( j=0; j<ei->dof[var]; j++)
                {
                   phi_j = bf[var]->phi[j];


                  shear_top_plus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_top_plus += (*(esp->sh_shear_top[iii]) + 0.5 * eps_top) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_top_plus += *(esp->sh_shear_top[iii]) * bf[var]->phi[iii];
                        }
                     }

                  shear_top_minus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_top_minus += (*(esp->sh_shear_top[iii]) - 0.5 * eps_top) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_top_minus += *(esp->sh_shear_top[iii]) * bf[var]->phi[iii];
                        }
                     }

                  source_plus = 0.0;
                  if (pd->e[eqn] & T_SOURCE)
                    {
                     source_plus += - phi_i * mu0 *
                                    (  pow( fabs(shear_top_plus), nexp - 1.) * shear_top_plus
                                     - pow( fabs(shear_bot), nexp - 1.) * shear_bot );

                     source_plus *= det_J * h3 * wt;
                     source_plus *= pd->etm[eqn][(LOG2_SOURCE)];
                    }


                  source_minus = 0.0;
                  if (pd->e[eqn] & T_SOURCE)
                    {
                     source_minus += - phi_i * mu0 *
                                     (  pow( fabs(shear_top_minus), nexp - 1.) * shear_top_minus
                                      - pow( fabs(shear_bot), nexp - 1.) * shear_bot );

                     source_minus *= det_J * h3 * wt;
                     source_minus *= pd->etm[eqn][(LOG2_SOURCE)];
                    }

                  lec->J[peqn][pvar][i][j] += (source_plus - source_minus)/eps_top;

                }
            }

          /* SENSITIVITY W.R.T. BOTTOM WALL SHEAR RATE */


          var = SHELL_SHEAR_BOT;

          if (pd->v[var])
            {
              pvar = upd->vp[var];
              for ( j=0; j<ei->dof[var]; j++)
                {

                  phi_j = bf[var]->phi[j];


                  shear_bot_plus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_bot_plus += (*(esp->sh_shear_bot[iii]) + 0.5 * eps_bot) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_bot_plus += *(esp->sh_shear_bot[iii]) * bf[var]->phi[iii];
                        }
                     }

                  shear_bot_minus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_bot_minus += (*(esp->sh_shear_bot[iii]) - 0.5 * eps_bot) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_bot_minus += *(esp->sh_shear_bot[iii]) * bf[var]->phi[iii];
                        }
                     }

                  source_plus = 0.0;
                  if (pd->e[eqn] & T_SOURCE)
                    {
                     source_plus += - phi_i * mu0 *
                                    (  pow( fabs(shear_top), nexp - 1.) * shear_top
                                     - pow( fabs(shear_bot_plus), nexp - 1.) * shear_bot_plus );

                     source_plus *= det_J * h3 * wt;
                     source_plus *= pd->etm[eqn][(LOG2_SOURCE)];
                    }


                  source_minus = 0.0;
                  if (pd->e[eqn] & T_SOURCE)
                    {
                     source_minus += - phi_i * mu0 *
                                     (  pow( fabs(shear_top), nexp - 1.) * shear_top
                                      - pow( fabs(shear_bot_minus), nexp - 1.) * shear_bot_minus );

                     source_minus *= det_J * h3 * wt;
                     source_minus *= pd->etm[eqn][(LOG2_SOURCE)];
                    }

                  lec->J[peqn][pvar][i][j] += (source_plus - source_minus)/eps_bot;

                }
            }

          /* SENSITIVITY W.R.T. PRESSURE */


          var = LUBP;

          if (pd->v[var])
            {
              pvar = upd->vp[var];
              for ( j=0; j<ei->dof[var]; j++)
                {

                  phi_j = bf[var]->phi[j];


                  for (jj = 0; jj < dim; jj++)
                     {
                      grad_phi_j[jj] = bf[var]->grad_phi[j][jj];
                      grad_II_phi_j[jj] = 0.0;
                     }

                  /* Perform I-nn grad operation */
                  Inn(grad_phi_j, grad_II_phi_j);


                  shear_top_plus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_top_plus += (*(esp->sh_shear_top[iii]) + 0.5 * eps_top) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_top_plus += *(esp->sh_shear_top[iii]) * bf[var]->phi[iii];
                        }
                     }

                  shear_top_minus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_top_minus += (*(esp->sh_shear_top[iii]) - 0.5 * eps_top) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_top_minus += *(esp->sh_shear_top[iii]) * bf[var]->phi[iii];
                        }
                     }

                  shear_bot_plus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_bot_plus += (*(esp->sh_shear_bot[iii]) + 0.5 * eps_bot) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_bot_plus += *(esp->sh_shear_bot[iii]) * bf[var]->phi[iii];
                        }
                     }

                  shear_bot_minus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_bot_minus += (*(esp->sh_shear_bot[iii]) - 0.5 * eps_bot) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_bot_minus += *(esp->sh_shear_bot[iii]) * bf[var]->phi[iii];
                        }
                     }


                  calculate_lub_q_v_nonnewtonian_sens ( mu,H, veloU, veloL, grad_II_P, phi_j, grad_II_phi_j,
                                                        shear_top_plus, shear_top_minus, shear_bot_plus, shear_bot_minus,
                                                        eps_top, eps_bot );



                  advection = 0.0;
                  if (pd->e[eqn] & T_ADVECTION)
                    {
                     advection += phi_i * H * LubAux->dgradP_mag_dP;

                     advection *= det_J * h3 * wt;
                     advection *= pd->etm[eqn][(LOG2_ADVECTION)];
                    }


                  lec->J[peqn][pvar][i][j] += advection;
                }
            }
        }

/* ************* ASSEMBLE JACOBIAN OF PRESSURE ************ */

      eqn   = R_LUBP;
      peqn = upd->ep[eqn];

      for ( i=0; i<ei->dof[eqn]; i++)
        {

          phi_i = bf[eqn]->phi[i];

          for (ii = 0; ii < dim; ii++)
             {
              grad_phi_i[ii] = bf[eqn]->grad_phi[i][ii];
              grad_II_phi_i[ii] = 0.0;
             }

          /* Perform I-nn grad operation */
          Inn(grad_phi_i, grad_II_phi_i);



          /* SENSITIVITY W.R.T. TOP WALL SHEAR RATE */


          var = SHELL_SHEAR_TOP;

          if (pd->v[var])
            {
              pvar = upd->vp[var];
              for ( j=0; j<ei->dof[var]; j++)
                {

                  phi_j = bf[var]->phi[j];


                  for (jj = 0; jj < dim; jj++)
                     {
                      grad_phi_j[jj] = bf[var]->grad_phi[j][jj];
                      grad_II_phi_j[jj] = 0.0;
                     }

                  /* Perform I-nn grad operation */
                  Inn(grad_phi_j, grad_II_phi_j);

                  shear_top_plus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_top_plus += (*(esp->sh_shear_top[iii]) + 0.5 * eps_top) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_top_plus += *(esp->sh_shear_top[iii]) * bf[var]->phi[iii];
                        }
                     }

                  shear_top_minus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_top_minus += (*(esp->sh_shear_top[iii]) - 0.5 * eps_top) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_top_minus += *(esp->sh_shear_top[iii]) * bf[var]->phi[iii];
                        }
                     }

                  shear_bot_plus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_bot_plus += (*(esp->sh_shear_bot[iii]) + 0.5 * eps_bot) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_bot_plus += *(esp->sh_shear_bot[iii]) * bf[var]->phi[iii];
                        }
                     }

                  shear_bot_minus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_bot_minus += (*(esp->sh_shear_bot[iii]) - 0.5 * eps_bot) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_bot_minus += *(esp->sh_shear_bot[iii]) * bf[var]->phi[iii];
                        }
                     }


                  calculate_lub_q_v_nonnewtonian_sens ( mu,H, veloU, veloL, grad_II_P, phi_j, grad_II_phi_j,
                                                        shear_top_plus, shear_top_minus, shear_bot_plus, shear_bot_minus,
                                                        eps_top, eps_bot );

                  diffusion = 0.0;
                  if (pd->e[eqn] & T_DIFFUSION)
                    {

                      for (ii = 0; ii < dim; ii++)
                         {
                          diffusion += - grad_II_phi_i[ii] * LubAux->dq_dshear_top[ii][j];
                         }

                      diffusion *= det_J * h3 * wt;
                      diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
                    }


                  lec->J[peqn][pvar][i][j] += diffusion;
                }
            }

          /* SENSITIVITY W.R.T. BOTTOM WALL SHEAR RATE */


          var = SHELL_SHEAR_BOT;

          if (pd->v[var])
            {
              pvar = upd->vp[var];
              for ( j=0; j<ei->dof[var]; j++)
                {

                  phi_j = bf[var]->phi[j];


                  for (jj = 0; jj < dim; jj++)
                     {
                      grad_phi_j[jj] = bf[var]->grad_phi[j][jj];
                      grad_II_phi_j[jj] = 0.0;
                     }

                  /* Perform I-nn grad operation */
                  Inn(grad_phi_j, grad_II_phi_j);

                  shear_top_plus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_top_plus += (*(esp->sh_shear_top[iii]) + 0.5 * eps_top) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_top_plus += *(esp->sh_shear_top[iii]) * bf[var]->phi[iii];
                        }
                     }

                  shear_top_minus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_top_minus += (*(esp->sh_shear_top[iii]) - 0.5 * eps_top) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_top_minus += *(esp->sh_shear_top[iii]) * bf[var]->phi[iii];
                        }
                     }

                  shear_bot_plus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_bot_plus += (*(esp->sh_shear_bot[iii]) + 0.5 * eps_bot) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_bot_plus += *(esp->sh_shear_bot[iii]) * bf[var]->phi[iii];
                        }
                     }

                  shear_bot_minus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_bot_minus += (*(esp->sh_shear_bot[iii]) - 0.5 * eps_bot) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_bot_minus += *(esp->sh_shear_bot[iii]) * bf[var]->phi[iii];
                        }
                     }


                  calculate_lub_q_v_nonnewtonian_sens ( mu,H, veloU, veloL, grad_II_P, phi_j, grad_II_phi_j,
                                                        shear_top_plus, shear_top_minus, shear_bot_plus, shear_bot_minus,
                                                        eps_top, eps_bot );

                  diffusion = 0.0;
                  if (pd->e[eqn] & T_DIFFUSION)
                    {

                      for (ii = 0; ii < dim; ii++)
                         {
                          diffusion += - grad_II_phi_i[ii] * LubAux->dq_dshear_bot[ii][j];
                         }

                      diffusion *= det_J * h3 * wt;
                      diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
                    }


                  lec->J[peqn][pvar][i][j] += diffusion;
                }
            }

          /* SENSITIVITY W.R.T. CROSS STREAM SHEAR STRESS */


          var = SHELL_CROSS_SHEAR;

          if (pd->v[var])
            {
              pvar = upd->vp[var];
              for ( j=0; j<ei->dof[var]; j++)
                {

                  phi_j = bf[var]->phi[j];


                  for (jj = 0; jj < dim; jj++)
                     {
                      grad_phi_j[jj] = bf[var]->grad_phi[j][jj];
                      grad_II_phi_j[jj] = 0.0;
                     }

                  /* Perform I-nn grad operation */
                  Inn(grad_phi_j, grad_II_phi_j);

                  shear_top_plus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_top_plus += (*(esp->sh_shear_top[iii]) + 0.5 * eps_top) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_top_plus += *(esp->sh_shear_top[iii]) * bf[var]->phi[iii];
                        }
                     }

                  shear_top_minus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_top_minus += (*(esp->sh_shear_top[iii]) - 0.5 * eps_top) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_top_minus += *(esp->sh_shear_top[iii]) * bf[var]->phi[iii];
                        }
                     }

                  shear_bot_plus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_bot_plus += (*(esp->sh_shear_bot[iii]) + 0.5 * eps_bot) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_bot_plus += *(esp->sh_shear_bot[iii]) * bf[var]->phi[iii];
                        }
                     }

                  shear_bot_minus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_bot_minus += (*(esp->sh_shear_bot[iii]) - 0.5 * eps_bot) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_bot_minus += *(esp->sh_shear_bot[iii]) * bf[var]->phi[iii];
                        }
                     }


                  calculate_lub_q_v_nonnewtonian_sens ( mu,H, veloU, veloL, grad_II_P, phi_j, grad_II_phi_j,
                                                        shear_top_plus, shear_top_minus, shear_bot_plus, shear_bot_minus,
                                                        eps_top, eps_bot );

                  diffusion = 0.0;
                  if (pd->e[eqn] & T_DIFFUSION)
                    {

                      for (ii = 0; ii < dim; ii++)
                         {
                          diffusion += - grad_II_phi_i[ii] * LubAux->dq_dcross_shear[ii][j];
                         }

                      diffusion *= det_J * h3 * wt;
                      diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
                    }


                  lec->J[peqn][pvar][i][j] += diffusion;
                }
            }

          /* SENSITIVITY W.R.T. PRESSURE */


          var = LUBP;

          if (pd->v[var])
            {
              pvar = upd->vp[var];
              for ( j=0; j<ei->dof[var]; j++)
                {

                  phi_j = bf[var]->phi[j];


                  for (jj = 0; jj < dim; jj++)
                     {
                      grad_phi_j[jj] = bf[var]->grad_phi[j][jj];
                      grad_II_phi_j[jj] = 0.0;
                     }

                  /* Perform I-nn grad operation */
                  Inn(grad_phi_j, grad_II_phi_j);


                  shear_top_plus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_top_plus += (*(esp->sh_shear_top[iii]) + 0.5 * eps_top) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_top_plus += *(esp->sh_shear_top[iii]) * bf[var]->phi[iii];
                        }
                     }

                  shear_top_minus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_top_minus += (*(esp->sh_shear_top[iii]) - 0.5 * eps_top) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_top_minus += *(esp->sh_shear_top[iii]) * bf[var]->phi[iii];
                        }
                     }

                  shear_bot_plus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_bot_plus += (*(esp->sh_shear_bot[iii]) + 0.5 * eps_bot) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_bot_plus += *(esp->sh_shear_bot[iii]) * bf[var]->phi[iii];
                        }
                     }

                  shear_bot_minus = 0.0;
                  for ( iii = 0; iii < ei->dof[var]; iii++)
                     {
                      if (iii == j)
                        {
                          shear_bot_minus += (*(esp->sh_shear_bot[iii]) - 0.5 * eps_bot) * bf[var]->phi[iii];
                        }
                      else
                        {
                          shear_bot_minus += *(esp->sh_shear_bot[iii]) * bf[var]->phi[iii];
                        }
                     }


                  calculate_lub_q_v_nonnewtonian_sens ( mu,H, veloU, veloL, grad_II_P, phi_j, grad_II_phi_j,
                                                        shear_top_plus, shear_top_minus, shear_bot_plus, shear_bot_minus,
                                                        eps_top, eps_bot );

                  diffusion = 0.0;
                  if (pd->e[eqn] & T_DIFFUSION)
                    {

                      for (ii = 0; ii < dim; ii++)
                         {
                          diffusion += - grad_II_phi_i[ii] * LubAux->dq_dp1[ii][j];
                         }

                      diffusion *= det_J * h3 * wt;
                      diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
                    }


                  lec->J[peqn][pvar][i][j] += diffusion;
                }
            }
        }



    } /* end of Assemble_Jacobian */




  return(status);
} /* end of assemble_lubrication_power_law */

/* End of mm_fill_shell.c */
