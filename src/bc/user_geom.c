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