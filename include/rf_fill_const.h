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

/*
 *$Id: rf_fill_const.h,v 5.1 2007-09-18 18:53:46 prschun Exp $
 */

/*
 * Revision history:
 *
 * $Log: not supported by cvs2svn $
 * Revision 5.0  2006/07/06 16:18:54  edwilke
 *
 * New Goma version: 'Farewell CRMPR, hello $3 a gallon gas!'
 *
 * Revision 4.1  2003/11/25 23:16:00  dalabre
 * The copyright statement has been updated for Goma and the version ratcheted
 * to 4.4.0.
 *
 * The makefile Goma.mk has been updated for Linux and Sun so that these
 * versions can be built easily (until the configure-make production is
 * complete); the Linux version is default.
 *
 * Added a prototype for function assemble_interface_extension_velocity_sic
 * in mm_fill_terms.h.
 *
 * Revision 4.0  2001/12/21 06:01:41  dalabre
 * Up Goma source code repository to V4.0. This identification 'coincides'
 * with our documentation upgrade. It will be tagged "Tora_Bora" in
 * recognition of the deep, unknown recesses that still remain in Goma.
 *
 * Revision 3.4  1999/12/09 15:27:09  pasacki
 * Function prototype cleanup, removal of unused variables, initialize some.
 *
 * Revision 3.3  1999/11/16 15:50:23  hkmoffa
 * Chemkin-Goma merger checkin.
 * This passed the FULL test suite last night, except for the single problem
 * that Duane pointed out had NaN's on. It also passed my test suite which
 * actually checks the goma answer against a blessed solution file.
 * That test suite has two mp problems in it, as well as 2 chemkin
 * problems.
 *
 * Revision 1.1.1.3  1999/11/03 22:44:28  hkmoffa
 * Import of goma src from engsci
 *
 * Revision 3.2  1999/10/20 15:09:40  pasacki
 * Centralized prototypes for rd_{mesh,exo,dpi}.c
 *
 * Revision 3.1  1999/08/10 20:49:44  tabaer
 * Installation of VBR sparse matrix format
 *
 * Revision 3.0  1999/01/07 06:14:33  dalabre
 * The IMPEACHMENT Version.
 *
 * Revision 2.9  1998/11/17 11:47:42  pasacki
 * More code clean-up. Note new include files!
 *
 * Revision 2.8  1998/03/17 18:12:56  rrrao
 * forgot another! Sorry Phil.
 *
 * Revision 2.7  1998/03/05 21:57:29  pasacki
 * o Fixed up bugs in distributed processing modules.
 * o Minor cosmetic fix of ex_opts in mm_post_proc.c eliminates warnings of
 *   nonexistent element order map.
 * o New logic in mm_bc.c diverts duplicate BC listing to a file if
 *   threshhold criteria are exceeded (too many nodesets, sidesets or BCs).
 *
 * Revision 2.6  1998/01/15 00:03:47  prschun
 * Some more odds-and-ends.  Bug fixes to modified newton.  OSF1...
 *
 * Revision 2.5  1997/09/19 05:40:15  pasacki
 * Prelude to PUMBAA. Significant EXODUS II and other changes.
 *
 * Revision 2.4  1997/09/02 19:03:40  dalabre
 * This is the second round of PROTOtype updates to Goma, primarily in files
 * el_*.*, md_*.c and rf_*.*; other changes are companion changes needed for
 * compatibility.
 *
 * Revision 2.3  1997/07/17 22:54:59  rrrao
 * Fill eqn for q1 elements
 * Suspension model changes
 * Probably more...
 *
 * Revision 2.2  1997/05/14 20:33:06  rrrao
 * added variable interpolation for fill equation
 *
 * Revision 2.1  1996/12/13 23:11:12  rrrao
 * added second sparse matrix for vof stuff.
 *
 * Revision 2.0  1996/09/30 22:03:07  prschun
 * FY97-REVISION 2.0
 *
 * Revision 1.11  1996/09/27 22:35:35  prschun
 * Hold on, I must commit all of this files to get rid of salsa.h
 *
 * Revision 1.10  1995/08/01  19:01:26  racairn
 * Improvements in mesh motion scheme - nonlinear arbitrary mesh motion
 * separates shear (mu) and dilational (lambda) components; now mu and lambda
 * correspond to the shear and bulk moduli.  Also sub-parametric mapping
 * appears to work.  New porous media stuff. . .
 *
 * Revision 1.9  1995/06/05  22:02:39  rrrao
 * minor-time bug fixes to velo_tangent ...
 *
 * Revision 1.8  1995/06/02  16:11:50  racairn
 * Further changes in implementation of boundary conditions - the main goal
 * here was to standardize the bookkeeping for the bc's and make it easier
 * to add bc's or update the methods used for the old bc's
 *
 * Revision 1.7  1995/04/13  23:07:28  prschun
 * Main addition here is the USER option for the property model (previously
 * all taken as CONSTANT) has been activated for several properties.
 * The general dependency of mat props on all dependent variables necessitated
 * some significant changes to the core volume assemble_ routines.
 * Other misc changes, bug fixes etc. include:
 * 	-Storage format for d_mat_prop arrays of the mp struct.
 * 	   (removed MAX_CONC subscript)
 * 	-Added a BOUSS model option on the Navier-Stokes Source card
 * 	  for the Boussinesq model. (Thermal and solutal buoyancy)
 * 	-Added new .c file called mm_std_models.c  which currently holds
 * 	  just the boussinesq model.  All other standard models should be
 * 	  placed there.
 * 	-Fixed bug in parse_parameter routine.
 * 	-Added some more error checking in mm_bc.c
 *
 * Revision 1.6  1995/03/07  21:46:40  pasacki
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
 * Revision 1.5  1995/02/22  22:22:07  racairn
 * Added pressure datum, convected lagrangian, and maybe some other stuff too?
 *
 * Revision 1.4  1995/01/26  22:49:37  racairn
 * I was starting this update because I thought we needed to fix a Mysterious
 * problem in the Dip coater (with unmovable free surface) apparently you shouldn't
 * run these types of problems (transient movable free surface) on the IBM.
 *
 * Revision 1.3  1994/12/22  22:04:10  racairn
 * This is a 'purified' version of GOMA with boundary conditions separated
 * into four separate functions (the first step in cleaning up this section).
 * Also, this version includes lagrangian mesh motion using linear or
 * nonlinear constitutive equations.
 *
 * Revision 1.2  1994/11/22  16:54:49  racairn
 * This version includes equations for mass transport of multiple chemical
 * species and dirichlet and flux boundary conditions. RAC 11/22/94
 * Also some new post processing option
 *
 * Revision 1.1.1.1  1994/05/24  19:56:33  pasacki
 * initial checkin
 *
 */

#ifndef GOMA_RF_FILL_CONST_H
#define GOMA_RF_FILL_CONST_H

#ifndef EXTERN
#define EXTERN extern
#endif

double calc_surf_det(int,                /* ielem  */
                     int,                /* iconnect_ptr */
                     int,                /* nodes_per_elem  */
                     int,                /* ielem_surf_dim  */
                     int,                /* id_side  */
                     int,                /* num_nodes_on_side */
                     int[]);             /* local_elem_node_id [] */

void fill_surf_shape(int,                /* num_local_nodes, number of nodes in this element   */
                     int,                /* ielem_surf_dim, physical dimension of element surface */
                     int,                /* ielem_type, type-code for this element              */
                     double[],           /* phi[], vector of shape functions for this element   */
                     double[][MAX_PDIM], /* grad_phi_L[MDE][MAX_PDIM],
                                          Gradients of shape functions for this element.
                                          Note, this function refers to coordinates
                                          alpha and beta (if 3D), i.e., the
                                          coordinates of the surface parameterization
                                          and NOT the local element coordinates.     */
                     double,             /* alpha   current integration point            */
                     double,             /* beta    coordinates of the surface
                                                      alpha is used for 2D and 3D
                                                      beta  is used only for 3D	      */
                     int,                /* id_side,  identity of element side           */
                     int,                /* num_nodes_on_side, number of nodes on element side */
                     int[]);             /* local_elem_node_id[], vector containing id's
                                               of local element nodes on the current side */

#if 0
void
alloc_sparse_arrays
   ( int **,           /* **ija, column pointer array                  */
            double **,        /* **a, nonzero array                           */
            double **,        /* **a_old, nonzero array                           */
            int,              /* Fill, flag to allocate space for either the
				    the whole banana or just the fill eqn     */
            int [],		/* node_to_fill[]                             */
	    Exo_DB *,		/* exo - ptr to whole mesh data */

	    Dpi *);		/* dpi - distr proc info */
#endif
#if 0
void
set_diag_to_one
   ( int,              /* I, Global node number                        */
            int,              /* var_type, identity of equation on the global node  */
            int,              /* var_type_subindex, subindex of variable type var_type */
            struct Aztec_Linear_Solver_System * );     /* a[],  nonzero array                          */
#endif

/*************** PROTOTYPE FUNCTIONS FOR  rf_fill_terms.c ********************/

double Ysurf_quad(double[],                     /* x[]  */
                  double[],                     /* phi[MDE]  */
                  struct elem_side_bc_struct *, /* *elem_side_bc  */
                  int);                         /* w  */

#if 0                                           /* old */
void 
fvelo_tangential_bc
   ( double *,         /* *func,   */
            double [][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                              /*  d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE] */
            double [],        /* x[]  */
            int,              /* dcl_node  */
            double,           /* vtangent  */
            double,           /* beta  */
            double,           /* alpha  */
            double *,         /* *xsurf  */
            double *,         /* *x_dot  */
            double,           /* tt  */
            double  );       /* dt  */
#endif

#endif
