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
 *$Id: rf_io_structs.h,v 5.1 2007-09-18 18:53:46 prschun Exp $
 */

#ifndef _RF_IO_STRUCTS_H
#define _RF_IO_STRUCTS_H

#include "std.h"		/* define proto if not already done so. */

#include "rf_io_const.h"	/* Just in case you haven't done so already. */

/*
 * Here's a handy dandy structure for describing the results that are in
 * an Exodus II database...
 */

struct Results_Description
{
  /*
   * Nodal variables (nv)...
   *       Note: if time_derivatives are included in the output file then
   *                  nnv = 2* TotalNVSolnOutput + TotalNVPostOutput
   *             else
   *                  nnv =  TotalNVSolnOutput + TotalNVPostOutput
   */
 int	nnv;			             /* Total number of nodal values
						in the ExodusII database */
 int    TotalNVSolnOutput;                   /* Total number of nv calced from
						extract_soln_vec() that aren't time
						derivatives */
 int    TotalNVPostOutput;                    /* Total number of nodal variables
						 calculated in post_process routine */
 int	nvtype[MAX_NNV];	             /* nv types */
 int	nvkind[MAX_NNV];	             /* nv kinds (i.e., species) */
 int    nvmatID[MAX_NNV];                    /* material index of the variable */
 char	nvname[MAX_NNV][MAX_VAR_NAME_LNGTH]; /* nv names */
 char	nvdesc[MAX_NNV][MAX_VAR_DESC_LNGTH]; /* nv long descriptions */
 char	nvunit[MAX_NNV][MAX_VAR_DESC_LNGTH]; /* nv physical units */
 int    nvderivative[MAX_NNV];               /* Boolean indicating whether the nodal
				      	        variable refers to a time deriviatve
					        or not */
 /*
  * Element variables (ev)...
  */
 int	nev;			             /* number of ev */
 int	evtype[MAX_NEV];	             /* ev types */
 char	evname[MAX_NEV][MAX_VAR_NAME_LNGTH]; /* ev names */
 char	evdesc[MAX_NEV][MAX_VAR_DESC_LNGTH]; /* nv long descriptions */
 char	evunit[MAX_NEV][MAX_VAR_DESC_LNGTH]; /* nv physical units */
 int    evderivative[MAX_NEV];               /* Boolean indicating whether the nodal
				      	        variable refers to a time deriviatve
					        or not */
 /*
  * Global variables (gv)...
  */
 int	ngv;			             /* number of gv */
 char	gvname[MAX_NGV][MAX_VAR_NAME_LNGTH]; /* gv names */
 /*
  * History variables (hv)...
  */
 int	nhv;				     /* number of hv */
 char	hvname[MAX_NHV][MAX_VAR_NAME_LNGTH]; /* names of hv */
};
typedef struct Results_Description RESULTS_DESCRIPTION_STRUCT;

struct  Command_line_command
{
  int    type;
  char   *string;
  int    i_val;
  double r_val;
};
char aprepro_command[1024];

#endif

