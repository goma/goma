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
 *  PLACID:
 *    This program solves for the surface production terms from a
 *  a chemkin reaction mechanism
 */

/*
 *$Id: mm_placid.c,v 5.1 2007-09-18 18:53:45 prschun Exp $
 */

/* Putting in a single #include here to get rid of that damn "empty
 * translation unit" warning on compilation... */

/*
 *  The whole program is put in an ifdef block around the usage of 
 *  chemkin. It doesn't make sense to compile this program if chemkin
 *  is not used.
 *
 */
#ifdef USE_CHEMKIN

/* Standard include files */

#include <stdlib.h>
#include <math.h>

/*
 * Include the GOMA standard header information
 */
#include "std.h"
/* Chemkin include files */
#include "cpc_defs.h"
#include "ck_chemkin_const.h"

/* GOMA include files  */

#include "rf_allo.h"
#include "rf_fem_const.h"
#include "rf_vars_const.h"
#include "mm_mp_const.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_mp_structs.h"
#include "mm_mp.h"
#include "dp_types.h"
#include "rf_util.h"
#include "mm_chemkin.h"

/******************************************************************************
*       STATIC ROUTINES DEFINED IN THIS FILE
******************************************************************************/

static void   calc_activity(double [], double []);
static double calc_damping(double *x, double *dx, int dim, int *);
static double calc_update_norm(double [], double dx[], int, double, double);
static double calc_t(double [], double [], int *, double);

#ifdef DEBUG_PLACID
static void print_stuff(int, double, int, double *[], double [], CK_NAME []);
static void print_stuff2(int, double, int, int, double, double, int, double,
                         double, double [], double [], double [], double [],
                         int, BOOLEAN, CK_NAME []);
static void print_header(int, int, int, double, int, double, double, double [],
                         double, double, double[], double [], double [], double
                         [], double [], CK_NAME []);
#endif

/******************************************************************************
*                    LAPACK PROTOTYPES
******************************************************************************/

extern FSUB_TYPE dgetrf_(int *, int *, double *, int *, int [], int *);
extern FSUB_TYPE dgetrs_(char *, int *, int *, double *, int *, int [],
                         double [], int *, int *, unsigned int);

/******************************************************************************
*   PROTOTYPES and PREPROC DIRECTIVES FOR MISC. ROUTINES NEEDED BY placid
******************************************************************************/


#ifndef MAX
#  define MAX(x,y) (( (x) > (y) ) ? (x) : (y))     /* max function */
#endif

#ifndef DAMPING
#  define DAMPING TRUE
#endif

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int 
placid_sd(int ifunc, int bulkFunc, double time_scale, double XGas[],
	  double TKelvin, double Ptherm, double soln_init[], int dim,
	  double reltol, double abstol, double sdot[], CK_NAME sname[])

    /*************************************************************************
     *
     * placid_sd():
     *
     *  Surface domain interface to placid
     *
     *
     *************************************************************************/
{
  int retn;

  retn = placid(ifunc, bulkFunc, time_scale,  XGas, TKelvin, Ptherm,
		soln_init, dim, reltol, abstol, sdot, sname);
  return retn;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int 
placid(int ifunc, int bulkFunc, double time_scale, double XGas[],
       double TKelvin, double Ptherm, double soln_init[], int dim,
       double reltol, double abstol, double sdot[], CK_NAME sname[])

/*****************************************************************************
 *
 * placid():
 *
 * The following calculation is a dimension dim=KkSurf+KkBulk Newton's 
 * method to get the surface fractions of the surface and bulk species 
 * by requiring that the surface species production rate = 0 and that 
 * the bulk fractions are  proportional to their production rates.  
 * The routine then returns surf_frac[]
 * and the corresponding sdot[] vector.
 *
 * Authors:    Andrew Salinger Dept. 1421  MP Comp. Sci.
 *             Harry Moffat    Dept. 1126, Chemical Processing Sciences
 *             Sandia National Labs
 *             Albuquerque, NM 87185
 *
 * Parameters:
 * ---------------------------------------------------------------------------
 *
 * ifunc: Functionality that the program is suppose to perform:
 *
 * 1: SFLUX_INITIALIZE   = This assumes that the initial guess supplied to the
 *                         routine is far from the correct one. Substantial
 *                         work plus transient time-stepping is to be expected
 *                         to find a solution.
 * 2:  SFLUX_RESIDUAL    = Need to solve the surface problem in order to
 *                         calculate the surface fluxes of gas-phase species.
 *                         (Can expect a moderate change in the solution
 *                          vector -> try to solve the system by direct methods
 *                          with no damping first -> then, try time-stepping
 *                          if the first method fails)
 *                         A "time_scale" supplied here is used in the
 *                         algorithm to determine when to shut off
 *                         time-stepping.
 * 3:  SFLUX_JACOBIAN    = Calculation of the surface problem is due to the
 *                         need for a numerical jacobian for the gas-problem.
 *                         The solution is expected to be very close to the
 *                         initial guess, and accuracy is needed.
 * 4:  SFLUX_TRANSIENT   = The transient calculation is performed here for an
 *                         amount of time specified by "time_scale".  It is
 *                         not garraunted to be time-accurate - just stable
 *                         and fairly fast. The solution after del_t time is
 *                         returned, whether it's converged to a steady
 *                         state or not.
 *
 *bulkFunc:     Functionality expected from the bulk phase.  This changes the
 *              equations that will be used to solve for the bulk mole
 *              fractions.
 * 1:  BULK_DEPOSITION   = deposition of a bulk phase is to be expected.
 *                         Bulk mole fractions are determined from ratios of
 *                         growth rates of bulk species.
 * 2:  BULK_ETCH         = Etching of a bulk phase is to be expected.
 *                         Bulk mole fractions are assumed constant, and given
 *                         by the initial conditions.
 *
 *time_scale:   Time over which to integrate the surface equations.
 *
 *XGas[]:       Value of the gas-phase mole fractions (KkGas long).
 *
 *TKelvin:      Temperature of the surface in Kelvin
 *
 *PTherm:       Pressure of the gas phase (cgs units).
 *
 *soln_init[]:  On input, it is the vector of initial guess for the surface
 *              problem. It has length, dim.
 *              On output, it is the solution of the surface problem.
 *              The first KkSur entries are the surface site fractions, while
 *              the next KkBulk entries are the bulk mole fractions.
 *
 *dim:          The size of the surface problem = KkSur + KkBulk.
 *
 *reltol:       Relative tolerance being used in the non-linear solver for
 *              the gas phase.  This program will try to solve the surface
 *              equations to a relative accuracy of 0.1*reltol.
 *abstol:       Absolute tolerance being used in the non-linear solver for
 *              the gas phase.  This program will try to solve the surface
 *              equations to a relative accuracy of 0.1*abstol.
 *
 *sdot[]:       This variable has length KkTot.
 *              on return, it returns the surface production rates (mole/cm**2)
 *              for all variables due to the surface reactions.
 *
 *ioflag:         0 - No output is produced by routine.
 *                1 - Small output is produced
 *                2 - Large output is produced
 *                3 - Jacobian is also dumped out.
 *                If the macro DEBUG_PLACID isn't predefined, no output is
 *                generated from this routine.  This feature cuts down on the
 *                total size of the program.
 *sname[]       Species names (vector of length KkTot) - each entry has
 *              type CK_NAME, which now, is a character array of length 16.
 *
 *Return Value:
 *---------------------------------------------------------------------------
 *
 * Extern variables assumed to exist as part of the CHEMKIN_GLOBS_STRUCT:
 *---------------------------------------------------------------------------

  int KkTot, KkGas, KkSurf, KkBulk:
                 Total number of species, gas species, surf species and bulk
                 species.
  int NnPhase:   Total number of phases
  int NfSurf, NlSurf: ID of first and last surface phase.
  int NfBulk, NlBulk: ID of first and last bulk phase.
  int KFirst[], KLast[], KkPhas[]: Vector of the first, the last and the total
                                   number of species in each phase.
  double  SDen0[]:  Intial densities of the surface sites.
  int     Kcov[]:   Number of surface sites each surface species covers
  int    *Iskwrk_p: Pointers to the Surface Chemkin work space
  double *Rskwrk_p:

 Important internal variables:
 ----------------------------------------------------------------------------

  (1) Reaction and Thermo Stuff

  act[]:        Array of activities (mole fractions for gas species, site
                fractions for surface species, and activities of bulk species)
                that function skrat_ uses to calculate production rates.
                The surface and bulk components of act are currently the same
                as surf_frac, but a hook is in place for non-unity activity
                coeffs.

  GRbulk:       Growth rate of current phase (mole/(cm**2*sec))

  dGRdXj[]      Derivative of GRbulk wrt solution unknowns
                Vector of length dim.

  sitdot[]:     Molar production rate of phase, as defined in Surface Chemkin
                input file.  The values are not currently used, but contain
                valuable information such as deposition rate.

  DsdotDX[][]:  Partial derivative of production rates wrt mole fractions.
                KkTot by KkTot in length.

  Xsoln[]:      Vector of length KkGas+dim containing the current solution,
                mole fraction units.  The first KkGas are gas phase mole
                fractions.  The next KkSurf are surface site fractions. The
                next KkBulk equations are bulk mole fractions.

  sden_total:   Total surface site density. Used as the time constant for the
                relaxation of the bulk mole fractions during a time transient.

  surf_frac[]:  Vector of length dim containing the current solution. surf_frac
                points into the Xsoln[] vector.


  (2) Newton's Method and Linear Solve Stuff

  damp:         The damping factor calculated by calc_damping to be multiplied
                by the update vector (resid returned by linear solve).

  info:         Error flag returned by LAPACK, currently ignored.

  ipiv[]:       Integer array of pivot rows passed between LAPACK routines.

  iter:         Iteration counter for local Newton's method.

  iter_max:     Maximum number of iterations in Newton's method.

  EXTRA_ACCURACY:A constant that is the ratio of the required update norm in
                this Newton iteration compared to that in the nonlinear solver.
                A value of 0.1 is used so surface species are safely
                overconverged.

  Jac[][]:      The "dim" by "dim" computed Jacobian matrix for the
                local Newton's method.

  resid[]:      The residual vector of length "dim" that, that has the value
                of "sdot" for surface species.  The residuals for the bulk
                species are a function of the sdots for all species in the bulk
                phase. The last residual of each phase enforces {Sum(fractions)
                = 1}. After linear solve (dgetrf_ & dgetrs_), resid holds the
                update vector.

  resid_norm:   Weighted L2 norm of the residual.  Currently, this is only
                used for IO purposes.  It doesn't control convergence.
                Therefore, it is turned off when DEBUG_PLACID isn't defined.

  update_norm:  Scalar holding the update norm to check convergence of the
                Newton iteration.

  inv_t:        1.0/Delta_t - the current time step.  Currently, this is
                determined from sdot[], so that the time step doesn't ever
                change the value of a variable by more than 100%.

  t_real:       Total value of the time integrated.


  Functions called:
  ----------------------------------------------------------------------------

  calc_activity   -- A hook to calculate activities from fractions.

  calc_damping    -- Routine to calculate damping factors so that the
                     fractions stay bounded from 0 to 1.

  skdsdx_    -- Surface Chemkin routine that returns the derivative of the
                sdots with respect to activities.  It also returns sdot and
                sitdot.

  skrat_     -- Surface Chemkin routine that returns the molar production rates
                of each species (gas, surphase, and bulk), sdot,  and the
                production rates of each phase, sitdot.

  dgetrf_    -- First half of LAPACK direct solve of a full Matrix

  dgetrs_    -- Second half of LAPACK direct solve of a full matrix. Returns
                solution vector in the right-hand-side vector, resid.

  ----------------------------------------------------------------------------

  ****************************************************************************/
{

#define EXTRA_ACCURACY  0.1

  static char cflag = 'N';
  int        kspec, kstart, kend, istart, iend, ispecial,
    irow, jcol, krow, kcol, phase, info;
  int        label_t=-1, label_d; /* Species IDs for time and damping control*/
  int        iter=0, iter_max=1000, nrhs=1;
  double     damp=1.0, sden_tot, GRbulk, RU, tmp, resid_norm,
    sum, inv_t, t_real = 0.0, update_norm = 1.0E6;
  double  **DsdotDX, **Jac, *Xsoln, *act, *resid, *soln_old, *dGRdXj, *sitdot,
    *surf_frac;
  int       *ipiv;
  BOOLEAN    do_time = FALSE, not_converged = TRUE;

  /*-------------------------------------------------------------------------*/

#ifdef DEBUG_PLACID
  double         t1;
  extern double  second(void);
  int ioflag = 2;
#endif
  /*-------------------------------------------------------------------------*/

#ifdef DEBUG_PLACID
  if (ioflag)  t1 = second();
#endif

  /*
   *       Set the initial value of the do_time parameter
   */

  if (ifunc == SFLUX_INITIALIZE || ifunc == SFLUX_TRANSIENT) do_time = TRUE;

  /*
   *  Allocate temporary storage -
   *      DsdotDX  = partial derivatives wrt unknowns  (KkTot, KkTot)
   *      Jac      = Jacobian of the surface problem (dim*dim)
   *      Xsoln    = Storage for the solution vector (KkGas+dim)
   *      act      = Vector of activities (KkTot)
   *      resid    = Vector of residuals (dim)
   *      soln_old = storage for previous solution vector (dim)
   *      dGRdXj   = storage for partial deriv of GR wrt solnvect (dim)
   *      sitdot   = storage for rates of change of site densitites (dim)
   *      ipiv     = pivot storage for direct solver (dim)
   */

  Jac      = alloc_dbl_2(dim, dim, DBL_NOINIT);
  Xsoln    = alloc_dbl_1(4*dim + Ck.KkGas + Ck.KkTot + Ck.NnPhase, 
                         DBL_NOINIT);
  act      = Xsoln       + (Ck.KkGas + dim);
  resid    = act         +  Ck.KkTot;
  soln_old = resid       +  dim;
  dGRdXj   = soln_old    +  dim;
  sitdot   = dGRdXj      +  dim;
  ipiv     = alloc_int_1(dim, INT_NOINIT);
  DsdotDX  = alloc_dbl_2(Ck.KkTot, Ck.KkTot, DBL_NOINIT);
  if (DsdotDX ==  NULL) {
    fprintf(stderr, "placid ERROR: out of memory");
    exit (-1);
  }

  /*
   *    Align surf_frac vector in Xsoln
   */

  surf_frac = Xsoln + Ck.KkGas;

  /*
   *    Get the gas constant, RU, in cgs units.
   */

  (void) ckrp_(Ck.Ickwrk_p, Ck.Rckwrk_p, &RU, &tmp, &sum);

  /*
   *  Calculate the sum of the surface site densities.  This will be used
   *  as the Lewis number for the bulk mole fractions.
   */

  sden_tot = 0.0;
  for (phase = Ck.NfSurf; phase <= Ck.NlSurf; phase++) {
    sden_tot += Ck.SDen0[phase];
  }

  /*
   * Activities of gas phase are constant, because gas phase is fixed.
   * Therefore, just have to store them in first positions of act[]
   * and forget them.
   */

  for (kspec = 0; kspec < Ck.KkGas; kspec++) {
    Xsoln[kspec] = MAX(XGas[kspec],0.0);
    act[kspec] = Xsoln[kspec];
  }

  /*
   *    Store the initial guess for the surface problem in the soln vector
   */

  for (irow = 0; irow < dim; irow++)  surf_frac[irow] = soln_init[irow];

#ifdef DEBUG_PLACID
  if (ioflag) {
    print_header(ioflag, ifunc, bulkFunc, time_scale, DAMPING,
		 reltol, abstol, XGas, TKelvin, Ptherm, soln_init, 
		 sdot, sitdot, Xsoln, act, sname);
  }
#endif

  /* ------------------------------------------------------------------
   *                         Start of Newton's method
   * ------------------------------------------------------------------
   */

  while (not_converged && iter < iter_max) {
    iter++;

    /*
     *    Store previous iteration's solution in the old solution vector
     */

    for (irow = 0; irow < dim; irow++)  soln_old[irow]  = surf_frac[irow];

    /*
     *  Calculate surface and bulk activities as a function of fractions
     *  -- just a hook for now
     *
     */
    calc_activity(act, Xsoln);

    /*
     *  Surface Chemkin calculates the derivatives of sdot wrt activities
     *  It also returns the corresponding values of sdot[] and sitdot[]
     */
    (void) skdsdx_(&Ptherm, &TKelvin, Xsoln, act, Ck.SDen0, Ck.Iskwrk_p,
                   Ck.Rskwrk_p, *DsdotDX, &Ck.KkTot, sdot, sitdot);

    /*
     *    Calculate the value of the time step
     *       - heuristics to stop large oscillations in deltaT
     */

    if (do_time) {
      tmp = calc_t(sdot, surf_frac, &label_t, time_scale);
      if (iter < 10)
        inv_t = tmp;
      else if (tmp > 2.0*inv_t)
        inv_t =  2.0*inv_t;
      else
        inv_t = tmp;

      /*
       *   Check end condition
       */

      if (ifunc == SFLUX_TRANSIENT) {
        tmp = t_real + 1.0/inv_t;
        if (tmp > time_scale) inv_t = 1.0/(time_scale - t_real);
      }
    }
    else
      inv_t = 0.0;

    /*
     *     Zero the Jacobian
     */
    init_vec_value(Jac[0], 0.0, dim*dim);

    /*
     *        Loop Over the Row Index Now
     *        Forming the residual and jacobian entry at the same time
     *       - Do this by phase, so each phase can be handled individually
     *         Skip the gas phase!
     */

    for (phase = 1; phase < Ck.NnPhase; phase ++) {

      /*
       *  Calculate some useful offsets for the phase
       */

      kstart = Ck.KFirst[phase];
      kend   = Ck.KLast[phase];
      istart = kstart        - Ck.KkGas;
      iend   = kend          - Ck.KkGas;

      /*
       *    Find the species with the maximum mole fraction
       *      - make that the special species
       */

      ispecial = istart;
      for (irow = istart; irow <= iend; irow++) {
        if (surf_frac[irow] > surf_frac[ispecial]) ispecial = irow;
      }

      /*
       *    Preprocessing for bulk phases with growth
       *    - Calculate the growth rate (mole/cm**2*sec)
       *    - Note: species with negative growth rates don't contribute
       *            negatively to the growth of that phase.
       *    - Calculate the derivative of the growth rate wrt mole fractions
       */

      if ((phase >= Ck.NfBulk) && (Ck.KkPhas[phase] > 1)) {
        GRbulk = 0.0;
        init_vec_value(dGRdXj, 0.0, dim);
        for (krow = kstart; krow <= kend; krow++) {
          if (sdot[krow] > 0.0) {
            GRbulk +=  sdot[krow];
            for (jcol = 0; jcol < dim; jcol++) {
              kcol = jcol + Ck.KkGas;
              dGRdXj[jcol] += DsdotDX[kcol][krow];
            }
          }
        }
      }

      /*
       *    Loop over the rows in the phase
       */

      for (irow = istart; irow <= iend; irow++) {

        /*
         *    Find the row number in the residual and jacobian
         */

        krow = irow + Ck.KkGas;

        /*
         *    Check to see if this is the special species in the phase
         *    - if it is, have to impose sum constraint
         */

        if (irow == ispecial) {
          sum = -1.0;
          for (jcol = istart; jcol <= iend; jcol++) {
            sum += surf_frac[jcol];
            Jac[jcol][irow] = 1.0;
          }
          resid[irow] = sum;
        }

        /*
         *    Next handle Surface Species Rows
         */

        else if (phase <= Ck.NlSurf) {

          resid[irow] =  -sdot[krow]*Ck.Kcov[krow];
          for (jcol = 0; jcol < dim; jcol++) {
            kcol = jcol + Ck.KkGas;
            Jac[jcol][irow] =  -DsdotDX[kcol][krow]*Ck.Kcov[krow];
          }
          if (do_time) {
            resid[irow]     +=
              inv_t * Ck.SDen0[phase] * (surf_frac[irow] - soln_old[irow]);
            Jac[irow][irow] += inv_t * Ck.SDen0[phase];
          }
        }

        /*
         *    Handle rows corresponding to bulk mole fractions
         */

        else if (phase <= Ck.NlBulk) {

          if (GRbulk == 0.0 || bulkFunc == BULK_ETCH) {
            resid[irow]     = surf_frac[irow] - soln_old[irow];
            Jac[irow][irow] = 1.0;
          }
          else {

            /*
             *        For phases where deposition occurs,
             *        we don't sdot's < 0 to contribution,
             *        because they would lead to negative mole fractions.
             */

            if (sdot[krow] >= 0.0) {
              resid[irow] = surf_frac[irow] * GRbulk - sdot[krow];
              for (jcol = 0; jcol < dim; jcol++) {
                kcol = jcol + Ck.KkGas;
                Jac[jcol][irow] = surf_frac[irow] * dGRdXj[jcol]
                  - DsdotDX[kcol][krow] ;
              }
            }
            else {
              resid[irow] = surf_frac[irow] * GRbulk;
              for (jcol = 0; jcol < dim; jcol++) {
                Jac[jcol][irow] = surf_frac[irow] * dGRdXj[jcol];
              }
            }
            Jac[irow][irow] += GRbulk;

            /*
             *     Add in the time dependence
             *     - The time constant is fairly arbitrary, here
             *       sden_tot at least has the correct units.
             */

            if (do_time) {
              resid[irow]     +=
                inv_t*sden_tot*(surf_frac[irow]-soln_old[irow]);
              Jac[irow][irow] += inv_t*sden_tot;
            }
          }

        } /* End of rows corresponding to bulk phases */
      }
    }

    /*
     *    Find the weighted norm of the residual
     */

    resid_norm = calc_update_norm(surf_frac, resid, dim, 1.0, 1.0e-6);

#ifdef DEBUG_PLACID
    /*
     *    Print out the residual and jacobian
     */

    print_stuff(ioflag, resid_norm, dim, Jac, resid, sname);
#endif

    /*
     *  Solve Linear system (with LAPACK).  The solution is in resid[]
     */

/* 
    if (resid_norm > 1.0e-7) {
*/
      (void) dgetrf_(&dim, &dim, *Jac, &dim, ipiv, &info);

    if (info==0) {
      (void) dgetrs_(&cflag, &dim, &nrhs, *Jac, &dim, ipiv, resid, &dim,
                     &info, 1);
    }
    /*
     *    Force convergence if residual is small to avoid
     *    "nan" results from the linear solve.
     */
    else {
      printf("Zero pivot, assuming converged: %g (%d)\n",resid_norm, info);
      for (jcol = 0; jcol < dim; jcol++) resid[jcol] = 0.0;
      if (do_time) t_real += time_scale;
#ifdef DEBUG_PLACID
      printf("\nResidual is small, forcing convergence!\n");
#endif
    }

    /*
     *    Calculate the Damping factor needed to keep all unknowns
     *    between 0 and 1, and not allow too large a change (factor of 2)
     *    in any unknown.
     */

#   if (DAMPING == TRUE)
    damp = calc_damping(surf_frac, resid, dim, &label_d);
#   endif

    /*
     *    Calculate the weighted norm of the update vector
     */

    update_norm = calc_update_norm(surf_frac, resid, dim, reltol, abstol);

    /*
     *    Update the solution vector and real time
     */

    for (irow = 0; irow < dim; irow++) surf_frac[irow] -= damp*resid[irow];

    if (do_time) t_real += damp/inv_t;

#ifdef DEBUG_PLACID
    if (ioflag) {
      print_stuff2(ioflag, damp, label_d, label_t,  inv_t, t_real, iter,
                   update_norm, resid_norm, sdot, surf_frac, resid, XGas,
                   dim, do_time, sname);
    }
#endif

    if (ifunc == SFLUX_TRANSIENT)
      not_converged = (t_real < time_scale);
    else {
      if (do_time) {
        if (t_real > time_scale) do_time = FALSE;
      }
      else {
        not_converged = (update_norm > EXTRA_ACCURACY);
      }
    }
  }  /* End of Newton's Method while statement */

  /*
   *  End Newton's method.  If not converged, print error message and
   *  recalculate sdot's at equal site fractions.
   */

#ifdef DEBUG_PLACID
  if (not_converged) {
    printf("#$#$#$# Error in placid $#$#$#$ \n");
    printf("Newton iter on surface species did not converge, "
           "update_norm = %e \n", update_norm);
    printf("Continuing anyway\n");
  }
  if (ioflag) printf("\nEnd of solve, time used: %e\n", second()-t1);
#endif

  /*
   *        Decide on what to return in the solution vector
   *        - right now, will always return the last solution
   *          no matter how bad
   */

  /*
   *        Calculate a consistent sdot for the solution vector to be returned
   */

  calc_activity(act, Xsoln);
  (void) skrat_(&Ptherm, &TKelvin, act, Ck.SDen0, Ck.Iskwrk_p, 
		Ck.Rskwrk_p, sdot, sitdot);

  /*
   *        Transfer the solution to be returned to the return vector
   */

  for (irow = 0; irow<dim; irow++) soln_init[irow] = surf_frac[irow];

  /*
   *        Free memory allocated within the routine.
   */

  safe_free((void **) &DsdotDX);
  safe_free((void **) &Jac);
  safe_free((void **) &Xsoln);
  safe_free((void **) &ipiv);

  /*
   *        Return with the appropriate flag
   */

  if (update_norm > EXTRA_ACCURACY)
    return -1;
  else
    return 1;

} /* placid */

#undef DAMPING

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

static void 
calc_activity(double act[], double Xsoln[])

/*
 * This function calculates the surface and bulk activities as a function of
 * the surface and bulk fractions.  Currently, this is just a hook.
 */

{
  int i;
  extern int KkGas, KkTot;

  for (i = KkGas; i < KkTot; i++) act[i] = Xsoln[i];
}

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

#define APPROACH 0.80

static double 
calc_damping(double x[], double dxneg[], int dim, int *label)

   /* This function calculates a damping factor for the Newton iteration update
    * vector, dxneg, to insure that all site and bulk fractions, x, remain
    * bounded between zero and one.
    *
    *      dxneg[] = negative of the update vector.
    *
    * The constant "APPROACH" sets the fraction of the distance to the boundary
    * that the step can take.  If the full step would not force any fraction
    * outside of 0-1, then Newton's method is allowed to operate normally.
    */
{
  int       i;
  double    damp = 1.0, xnew, xtop, xbot;
  static double damp_old = 1.0;

  *label = -1;

  for (i = 0; i < dim; i++) {

    /*
     * Calculate the new suggested new value of x[i]
     */

    xnew = x[i] - damp * dxneg[i];

    /*
     *  Calculate the allowed maximum and minimum values of x[i]
     *   - Only going to allow x[i] to converge to zero by a
     *     single order of magnitude at a time
     */

    xtop = 1.0 - 0.1*fabs(1.0-x[i]);
    xbot = fabs(x[i]*0.1) - 1.0e-16;
    if (xnew > xtop )  {
      damp = - APPROACH * (1.0 - x[i]) / dxneg[i];
      *label = i;
    }
    else if (xnew < xbot) {
      damp = APPROACH * x[i] / dxneg[i];
      *label = i;
    } else  if (xnew > 3.0*MAX(x[i], 1.0E-10)) {
      damp = - 2.0 * MAX(x[i], 1.0E-10) / dxneg[i];
      *label = i;
    }
  }

  if (damp < 1.0e-6) damp = 1.0e-6;
  /*
   * Only allow the damping parameter to increase by a factor of three each
   * iteration. Heuristic to avoid oscillations in the value of damp
   */

  if (damp > damp_old*3) {
    damp = damp_old*3;
    *label = -1;
  }

  /*
   *      Save old value of the damping parameter for use
   *      in subsequent calls.
   */

  damp_old = damp;
  return damp;

} /* calc_damping */
#undef APPROACH

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

static double 
calc_update_norm(double x[], double dx[], int dim, double reltol,
		 double abstol)

  /*
   *    This function calculates the update norm using the same reltol
   *    and abstol as the nonlinear solver.  The convergence of the
   *    surface species is forced to converge better than the nonlinear
   *    solver by the variable factor (factor < 1)
   */
{
  int         i;
  double      norm = 0.0, tmp;

  for (i = 0; i < dim; i++) {
    tmp   = dx[i] / (reltol * fabs(x[i]) + abstol);
    norm += tmp * tmp;
  }

  return (sqrt(norm/dim));

} /* calc_update_norm */

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

static double 
calc_t(double sdot[], double surf_frac[], int *label, double time_scale)

  /*
   *    This routine calculates a pretty conservative 1/del_t based
   *    on  MAX_i(sdot_i/(X_i*SDen0)).  This probably guarantees
   *    diagonal dominance.
   *
   *     Small surface fractions are allowed to intervene in the del_t
   *     determination, no matter how small.  This may be changed.
   *     Now minimum changed to 1.0e-12,
   *
   *     Maximum time step set to time_scale.
   */
{
  int      phase, i, k, istart, iend;
  double   inv_t = 0.0, tmp;

  for (phase = 1; phase <= Ck.NlSurf; phase ++) {
    istart = Ck.KFirst[phase]  - Ck.KkGas;
    iend   = Ck.KLast[phase]   - Ck.KkGas;
    for (i = istart; i <= iend; i++) {
      k = i + Ck.KkGas;
      if (surf_frac[i] <= 1.0E-12)
        tmp = 1.0E-12;
      else
        tmp = surf_frac[i];
      tmp = fabs(Ck.Kcov[k]*sdot[k]) /(tmp*Ck.SDen0[phase]);
      if (tmp > inv_t) {
        inv_t = tmp;
        *label = i;
      }
    }
  }

  return MAX(inv_t, 1.0 / time_scale);

} /* calc_t */

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

#ifdef DEBUG_PLACID

static void
 print_stuff(int ioflag, double norm, int dim, double *Jac[],
	     double resid[], CK_NAME sname[])
{
  int i, j;

  if (ioflag > 1) {
    printf("Printout of residual and jacobian\n");
    printf("\tResidual: norm = %10.4e\n", norm);
    for (i = 0; i < dim; i++) {
      printf("\t%d: %.16s: %e\n", i, sname[i+KkGas], resid[i]);
    }
  }
  if (ioflag >2) {
    printf("\tJacobian:\n");
    for (i = 0; i < dim; i++) {
      printf("Row %d:%.16s:\n", i, sname[i+KkGas]);
      for (j = 0; j < dim; j++) {
        printf("%10.4e ", Jac[j][i]);
      }
      printf("\n");
    }
  }

} /* print_stuff */

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

static void
print_header(int ioflag, int ifunc, int bulkFunc,
	     double time_scale, int damping, double reltol,
	     double abstol, double XGas[], double TKelvin,
	     double Ptherm, double soln_init[], double sdot[],
	     double sitdot[], double Xsoln[], double act[],
	     CK_NAME sname[])
{
  int i, k;
  extern int KkGas, KkBulk, KkSurf;

  if (ioflag) {

    if (ifunc == SFLUX_INITIALIZE) {
      printf("\nPLACID Called with Initialization turned on\n");
      printf("  Time scale input = %9.3e\n", time_scale);
    }
    else if (ifunc == SFLUX_RESIDUAL) {
      printf("\n PLACID Called to calculate steady state residual\n");
      printf( "         from a good initial guess\n");
    }
    else if (ifunc == SFLUX_JACOBIAN)  {
      printf("\n PLACID Called to calculate steady state jacobian\n");
      printf( "         from a good initial guess\n");
    }
    else if (ifunc == SFLUX_TRANSIENT) {
      printf("\n PLACID Called to integrate surface in time\n");
      printf( "         for a total of %9.3e sec\n", time_scale);
    }
    else {
      fprintf(stderr,"Unknown ifunc flag = %d\n", ifunc);
      exit (-1);
    }

    if (bulkFunc == BULK_DEPOSITION)
      printf("  Deposition is to be expected\n");
    else if (bulkFunc == BULK_ETCH)
      printf("  Etching is to be expected\n");
    else {
      fprintf(stderr,"Unknown bulkFunc flag = %d\n", bulkFunc);
      exit (-1);
    }

    if (damping)
      printf("  Damping is ON   \n");
    else
      printf("  Damping is OFF  \n");
    printf("  Reltol = %9.3e, Abstol = %9.3e\n", reltol, abstol);

  }

  /*
   *   Print out the initial guess
   */

  if (ioflag > 1) {
    calc_activity(act, Xsoln);
    skrat_(&Ptherm, &TKelvin, act, SDen0, Iskwrk_p, Rskwrk_p, sdot, sitdot);
    printf("\n============================= INITIAL GUESS "
           "======================\n");
    printf("     Temperature = %10.3e Kelvin\n", TKelvin);
    printf("     Pressure    = %10.3e ergs/cm2\n\n", Ptherm);
    printf("  Name            Prod_Rate   Fraction\n");
    printf("--------------------------------------------------------------"
           "---\n");
    for (k = 0; k < KkGas; k++) {
      printf("%.16s %10.3e %10.3e\n", sname[k], sdot[k], XGas[k]);
    }
    for (i = 0; i < KkSurf+KkBulk; i++) {
      k = i + KkGas;
      printf("%.16s %10.3e %10.3e\n", sname[k], sdot[k], soln_init[i]);
    }
    printf("---------------------------------------------------------------"
           "--\n");
  }

  if (ioflag == 1) {
    printf("\n\n Iter    Time       Del_t      Damp      DelX    Resid       "
           "Name-Time    Name-Damp\n");
  }

} /* print_header */

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

static void 
print_stuff2(int ioflag, double damp, int label_d, int label_t,
	     double inv_t, double t_real, int iter,
	     double update_norm, double resid_norm,
	     double sdot[], double surf_frac[],
	     double resid[], double XGas[], int dim,
	     BOOLEAN do_time, CK_NAME sname[])
{
  int i, k;
  if (ioflag == 1) {

    printf("%6d ", iter);
    if (do_time)
      printf("%9.4e %9.4e ", t_real, 1.0/inv_t);
    else
      for (i = 0; i < 20; i++) printf(" ");
    if (damp < 1.0)
      printf("%9.4e ", damp);
    else
      for (i = 0; i < 10; i++) printf(" ");
    printf("%9.4e %9.4e", update_norm, resid_norm);
    if (do_time)
      printf(" %.16s", sname[KkGas+label_t]);
    else
      for (i = 0; i < 16; i++) printf(" ");
    if (label_d >= 0) printf(" %.16s", sname[KkGas+label_d]);
    printf("\n");
  }
  else if (ioflag > 1) {
    printf("\n=============================Iteration %5d "
           "========================\n",
           iter);
    if (do_time) {
      printf(" Real Time = %10.4e sec\n", t_real);
      printf(" Delta t   = %10.4e sec\n", 1.0/inv_t);
    }
    if (damp < 1.0)
      printf(" Damping value =  %10.4e\n", damp);
    printf(" Weighted norm of update = %10.4e\n", update_norm);
    printf(" Weighted norm of residual = %10.4e\n\n", resid_norm);

    printf("  Name            Prod_Rate   Fraction    Fraction_Old      ");
    if (damp < 1.0) printf(" UnDamped_Fraction");
    printf("\n");
    printf("---------------------------------------------------------------"
           "--\n");
    for (k = 0; k < KkGas; k++) {
      printf("%.16s %10.3e %10.3e\n", sname[k], sdot[k], XGas[k]);
    }
    for (i = 0; i < dim; i++) {
      k = i + KkGas;
      printf("%.16s %10.3e %10.3e %10.3e ", sname[k],
             sdot[k], surf_frac[i], surf_frac[i]+damp*resid[i]);
      if (damp < 1.0) {
        printf("%10.4e ", surf_frac[i]+(damp-1.0)*resid[i]);
        if (label_d == i) printf("Damp ");
      }
      if (label_t == i) printf("Tctrl");
      printf("\n");
    }
    printf("---------------------------------------------------------------"
           "--\n");
  }

} /* print_stuff2 */

#endif  /*  DEBUG_PLACID */
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

#endif  /*  USE_CHEMKIN  */
