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
 

#ifdef USE_RCSID
static const char rcsid[] = "$Id: user_continuation.c,v 5.3 2008-03-22 00:55:51 hkmoffa Exp $";
#endif

/* GOMA include files */

#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_eh.h"
#include "ac_update_parameter.h"
#include "dp_types.h"
#include "dpi.h"
#include "exo_struct.h"

#define GOMA_USER_AC_C
#include "user_continuation.h"

#define BC 1
#define MT 2
#define AC 3
#define UM 4

void
update_user_parameter(double lambda, double *x, double *xdot, double *x_AC,
                      Comm_Ex *cx, Exo_DB *exo, Dpi *dpi)
/************************************************************************
 *                                                                      *
 *    Created by Ed Wilkes 11/28/2001					*
 *                                                                      *
 *    This function allows the user to specify a set of arbitrary	*
 *    continuation conditions (quantities which must be updated at	*
 *    each step of a continuation run).	 These quantities can be	*
 *    boundary condition floats(BC), material properties(MT), or	*
 *    parameters in user-defined property models(UM), as in single	*
 *    parameter continuation.						*
 *    The main requirement is that the update values must be		*
 *    functions of the continuation parameter lambda, which is		*
 *    defined in the input deck. Each time a continuation algorithm	*
 *    requests a parameter update, the current value of lambda is	*
 *    passed to here so that the update values can be evaluated.	*
 *    For each BC, MT, AC, or UM update required, the user must provide	*
 *    1) The information required to identify the quantity:		*
 *            Type: BC, MT, AC, or UM.					*
 *		BC: provide BCID, DFID					*
 *		MT: provide MTID, MPID, MDID (MDID not used)		*
 *		AC: provide BCID, DFID					*
 *		UM: provide MTID, MPID, MDID				*
 *	 These are the same ID's which would be required in the		*
 *	 Continuation or Hunting sections of the input deck.		*
 *	 Note that in order to use a UM condition, a material		*
 *	 property model must be provided in file "user_mp.c"		*
 *	 and the property model must be specified as "USER" in		*
 *	 the *.mat file, followed by a user parameter list.		*
 *    2) The necessary line(s) of C code to calculate each update	*
 *	 value as a function of lambda [i.e. value = f(lambda); ]	*
 *    3) A standard call to function do_user_update. This in turn	*
 *	 calls the appropriate update_{BC|MT|UM}_parameter function	*
 *	 with all of the arguments needed for that function.		*
 *    Each evaluation and update will be performed sequentially		*
 *    in the order given - no loop is used.				*
 *                                                                      *
 *    This function applies only to the first continuation parameter    *
 *    When LOCA bifurcation tracking algorithms are used, the user      *
 *    may supply update functions for the second (TP) parameter in      *
 *    the second function "update_user_TP_parameter" just below this    *
 *    one - be careful to place each update in the correct function.    *
 *                                                                      *
 *    A (commented out) set of example entries for each continuation	*
 *    type is provided below.						*
 *                                                                      *
 *    NOTES:                                                            *
 *    1) This function is invoked when the Continuation Type		*
 *	 specified in the input deck is "UF" - in this case,		*
 *	 hunting or continuation condition cards are not needed		*
 *	 because the information is provided here. Also, the ID		*
 *	 type cards from the input deck are not used here.		*
 *    2) If all update quantities are linear functions of the		*
 *	 continuation parameter lambda (i.e. value = a + b * lambda),	*
 *	 then a set of continuation condition (CC or TC) cards may be	*
 *	 provided in the input deck in lieu of using this function.	*
 *    3) To save time, intermediate quantities which are used in more	*
 *	 than one update value function, such as lambda^2, can be	*
 *	 pre-calculated and stored in temporary variables for use by	*
 *	 these functions.  Just be sure to declare any such variables.	*
 *	 The same may be done for any physical constants (i.e. pi, g).	*
 *    4) The solution vector *x, time derivatives *xdot, and augmenting *
 *       condition extra unknowns *x_AC are also passed into this       *
 *       routine and can be used in calculating the update values!      *
 *                                                                      *
 *    EXAMPLE:                                                          *
 *    Here, the angle theta is a continuation parameter for a		*
 *    non-isothermal rotating fluid disk. Its density is a function	*
 *    of theta, but its surface tension is temperature-dependent,	*
 *    and a model of type sigma = a - b * T must be specified. The	*
 *    parameter b in that model is also dependent on theta. The third	*
 *    boundary condition card takes sin(theta) as its second float.	*
 *    All three of these quantities must be updated at each		*
 *    continuation step. The user surface tension model has been	*
 *    entered in user_mp.c and specified in the .mat file. The		*
 *    Continuation Type = UF card is in the input deck, and the cards	*
 *    for value and step size ranges for theta are there.		*
 *                                                                      *
 *    This function would then be filled in as follows:			*
 *    (Note that since the entire example is commented out,		*
 *     the comment delimiters are changed to "**")			*

{
  static int first_cp = 1;
  int first_tp = -1;
  int n = 0;
  int Type, BCID, DFID, MTID, MPID, MDID;
  double value;
** Declare any additional variables here **
  double stheta, ctheta, rho_0, g;

** Evaluate any intermediate quantities and/or assign constants here **
  stheta = sin(lambda);
  ctheta = cos(lambda);
  rho_0 = 0.0125;
  g = 980.6;

** Enter each continuation condition in this space **

** Density (MT): **

** ID's **
  Type = MT;
  MTID = 1;
  MPID = 1700;
  MDID = 0;

** Value **
  value = 1.0 + rho_0 * ctheta;

** Update call **
  n = do_user_update(n, first_cp, first_tp, Type, BCID, DFID,
                     MTID, MPID, MDID, value, cx, exo, dpi);


** Surface tension (UM): **

** ID's **
  Type = UM;
  MTID = 1;
  MPID = 1400;
  MDID = 1;

** Value **
  value = 0.1 * lambda;

** Update call **
  n = do_user_update(n, first_cp, first_tp, Type, BCID, DFID,
                     MTID, MPID, MDID, value, cx, exo, dpi);


** Force (BC): **

** ID's **
  Type = BC;
  BCID = 2;
  DFID = 1;

** Value **
  value = g * ctheta;

** Update call - copy from example **
  n = do_user_update(n, first_cp, first_tp, Type, BCID, DFID,
                     MTID, MPID, MDID, value, cx, exo, dpi);


** Done **
  first_cp = 0;
  return;
}

 END OF EXAMPLE */

/************************************************************************/

/* Actual function begins here */
{
 // static int first_cp = 1;
 // int first_tp = -1;
 // int n = 0;
 // int Type, BCID, DFID, MTID, MPID, MDID;
 // double value;
/* Declare any additional variables here */


/* If using this function, comment out this line. */
  EH(-1, "No user continuation conditions entered!");


/* Evaluate any intermediate quantities and/or assign constants here */


/* Enter each continuation condition in sequence in this space */


/* ID's */


/* Value */


/* Update call - copy from example */


/* Done */
 // first_cp = 0;
  return;

} /* END of routine update_user_parameter */
/****************************************************************************/

void
update_user_TP_parameter(double lambda, double *x, double *xdot,
                         double *x_AC, Comm_Ex *cx, Exo_DB *exo, Dpi *dpi)
/*
 * This function handles updates to the second (TP) parameter
 * when LOCA bifurcation tracking algorithms (turning point,
 * pitchfork, or Hopf are used. It works the same as the above
 * function, but here lambda is the TP parameter.
 *
 * The examples in the above function also apply identically to this
 * one, and for brevity are not repeated. The standard call to
 * do_user_update is also the same as above. The main difference is
 * the variables first_cp and first_tp are set up to indicate which
 * parameter is being updated and whether this is the first update
 * call for that parameter. The user needn't worry about this detail.
 */

{
/* Actual function begins here (refer to previous example) */
  //static int first_tp = 1;
  // int first_cp = -1;
  // int n = 0;
  // int Type, BCID, DFID, MTID, MPID, MDID;
  // double value;
/* Declare any additional variables here */


/* If using this function, comment out this line. */
  EH(-1, "No user TP continuation conditions entered!");


/* Evaluate any intermediate quantities and/or assign constants here */


/* Enter each continuation condition in sequence in this space */


/* ID's */


/* Value */


/* Update call - copy from example */


/* Done */
  //first_tp = 0;
  return;

} /* END of routine update_user_TP_parameter */
/****************************************************************************/


int
do_user_update(int n, int first_cp, int first_tp,
               int Type, int BCID, int DFID, int MTID, int MPID, int MDID,
               double value, Comm_Ex *cx, Exo_DB *exo, Dpi *dpi)
/*
 * This function calls the appropriate parameter update function
 * for each user-defined continuation condition type and saves
 * the calculated values. On the first call from each condition,
 * the ID's are also saved. The incremented function counter is returned.
 */
{
  int param_type = -1;    /* 1 = first (CP), 2 = second (TP) */

/* Determine which parameter is being updated */
  if (first_cp >= 0 && first_tp == -1)
    {
      param_type = 1;
    }
  else if (first_tp >= 0 && first_cp == -1)
    {
      param_type = 2;
    }
  else
    {
      EH(-1, "User parameter type settings inconsistent!");
    }

/* Store new and old values */
  if (param_type == 1)
    {
      cpuc[n].old_value = cpuc[n].value;
      cpuc[n].value = value;
    }
  else if (param_type == 2)
    {
      tpuc[n].old_value = tpuc[n].value;
      tpuc[n].value = value;
    }

/* Call the appropriate update function */
  switch (Type) {
    case BC:
      update_BC_parameter(value, BCID, DFID, cx, exo, dpi);

/* Store ID info on first call */
      if (first_cp == 1)
        {
          cpuc[n].Type = 1;
          cpuc[n].BCID = BCID;
          cpuc[n].DFID = DFID;
        }
      else if (first_tp == 1)
        {
          tpuc[n].Type = 1;
          tpuc[n].BCID = BCID;
          tpuc[n].DFID = DFID;
        }
      break;

    case MT:
      update_MT_parameter(value, MTID-1, MPID, MDID, cx, exo, dpi);

/* Store ID info on first call */
      if (first_cp == 1)
        {
          cpuc[n].Type = 2;
          cpuc[n].MTID = MTID-1;
          cpuc[n].MPID = MPID;
        }
      else if (first_tp == 1)
        {
          tpuc[n].Type = 2;
          tpuc[n].MTID = MTID-1;
          tpuc[n].MPID = MPID;
        }
      break;

    case AC:
      update_AC_parameter(value, BCID, DFID, cx, exo, dpi);

/* Store ID info on first call */
      if (first_cp == 1)
        {
          cpuc[n].Type = 3;
          cpuc[n].BCID = BCID;
          cpuc[n].DFID = DFID;
        }
      else if (first_tp == 1)
        {
          tpuc[n].Type = 3;
          tpuc[n].BCID = BCID;
          tpuc[n].DFID = DFID;
        }
      break;

    case UM:
      update_UM_parameter(value, MTID-1, MPID, MDID, cx, exo, dpi);

/* Store ID info on first call */
      if (first_cp == 1)
        {
          cpuc[n].Type = 4;
          cpuc[n].MTID = MTID-1;
          cpuc[n].MPID = MPID;
          cpuc[n].MDID = MDID;
        }
      else if (first_tp == 1)
        {
          tpuc[n].Type = 4;
          tpuc[n].MTID = MTID-1;
          tpuc[n].MPID = MPID;
          tpuc[n].MDID = MDID;
        }
      break;
    }

/* Update counter and return */
  n++;
  return n;
}
/****************************************************************************/
/* End of file user_continuation.c */
