/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2025 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/

#include "bc/user_geom.h"
#include "mm_as.h"

/*****************************************************************************/
/*       Functions for user defined solid boundary description               */
/*****************************************************************************/

dbl fnc(const dbl x1 MAYBE_UNUSED,
        const dbl x2 MAYBE_UNUSED,
        const dbl x3 MAYBE_UNUSED,
        const dbl p[] MAYBE_UNUSED,
        const dbl time MAYBE_UNUSED) {

  /* for PI use M_PIE Constant from std.h include file. */
  /* cylinder with axis parallel to z */
  // dbl cx = p[0], cy=p[1], r=p[2];
  dbl f = 0;

  return (f); /* Here's a good default behavior! */
}

dbl dfncd1(const dbl x1 MAYBE_UNUSED,
           const dbl x2 MAYBE_UNUSED,
           const dbl x3 MAYBE_UNUSED,
           const dbl p[] MAYBE_UNUSED,
           const dbl time MAYBE_UNUSED) {
  /* for PI use M_PIE Constant from std.h include file. */
  dbl f = 0;

  return (f); /* Here's a good default behavior! */

  /* Example code fragments:
   *
   *  dfdx1 = 1.0;
   *  return(dfdx1);
   *  f = -1;   2d fiber
   *  f = 2.*x1;   circle
   *  return f;
   */
}

dbl dfncd2(const dbl x1 MAYBE_UNUSED,
           const dbl x2 MAYBE_UNUSED,
           const dbl x3 MAYBE_UNUSED,
           const dbl p[] MAYBE_UNUSED,
           const dbl time MAYBE_UNUSED) {
  /* for PI use M_PIE Constant from std.h include file. */
  dbl f = 0;  /* dfdx2; */
  return (f); /* Here's a good default behavior! */

  /* Example code fragments:
   *
   *  dfdx2 = 2*x2;
   *  return(dfdx2);
   *
   *  f = 1./13 -sin(M_PIE*x2)/24. -M_PIE*x2*cos(M_PIE*x2)/24;   2d fiber
   *  f =  2.*x2;   circle
   *  return f;
   */
}

dbl dfncd3(const dbl x1 MAYBE_UNUSED,
           const dbl x2 MAYBE_UNUSED,
           const dbl x3 MAYBE_UNUSED,
           const dbl p[] MAYBE_UNUSED,
           const dbl time MAYBE_UNUSED) {
  /* for PI use M_PIE Constant from std.h include file. */
  dbl f = 0;

  f = 0.0; /* expanding sphere */

  return (f); /* Here's a good default behavior! */

  /* Example code fragments:
   *
   *
   *  f = -(1-x3*x3/200)/50.;
   *  f = 0.;
   *  return f;
   */
}

void fspline(int ielem_dim,
             double *func,
             double d_func[], /* dimensioned [MAX_VARIABLE_TYPES + MAX_CONC] */
             double p[],      /* parameters to parameterize temperature eqn model*/
             double time)     /* time  at which bc's are evaluated     */
{
  if (af->Assemble_LSA_Mass_Matrix)
    return;

  d_func[MESH_DISPLACEMENT1] = dfncd1(fv->x[0], fv->x[1], fv->x[2], p, time);

  d_func[MESH_DISPLACEMENT2] = dfncd2(fv->x[0], fv->x[1], fv->x[2], p, time);

  if (ielem_dim == 3)
    d_func[MESH_DISPLACEMENT3] = dfncd3(fv->x[0], fv->x[1], fv->x[2], p, time);

  *func = fnc(fv->x[0], fv->x[1], fv->x[2], p, time);

} /* END of routine fspline                                                  */

void fspline_rs(int ielem_dim,
                double *func,
                double d_func[], /* dimensioned [MAX_VARIABLE_TYPES + MAX_CONC] */
                double p[],      /* parameters to parameterize temperature eqn model*/
                double time)     /* time  at which bc's are evaluated     */
{
  if (af->Assemble_LSA_Mass_Matrix)
    return;

  d_func[SOLID_DISPLACEMENT1] = dfncd1(fv->x[0], fv->x[1], fv->x[2], p, time);

  d_func[SOLID_DISPLACEMENT2] = dfncd2(fv->x[0], fv->x[1], fv->x[2], p, time);

  if (ielem_dim == 3)
    d_func[SOLID_DISPLACEMENT3] = dfncd3(fv->x[0], fv->x[1], fv->x[2], p, time);

  *func = fnc(fv->x[0], fv->x[1], fv->x[2], p, time);

} /* END of routine fspline_rs */