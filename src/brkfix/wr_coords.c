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

/* write_coords() -- write out node point coordinates to file
 *
 * Notes:
 *
 *	[1]	The file format is based on that required for input to
 *		Chaco 2.0 
 *
 * Created: 1997/05/08 08:19 MDT pasacki@sandia.gov
 *
 * Revised:
 */

#define _WR_COORDS_C

#include <stdio.h>
#include <stdlib.h>

#include "map_names.h"
#include "std.h"
#include "eh.h"
#include "exo_struct.h"
#include "wr_coords.h"

void
write_coords(char *fn,		/* filename for coordinates */
	     Exo_DB *x)		/* pointer to EXODUS II db structure */
{
  int i;
  int err;
  FILE *s;

  s = fopen(fn, "w");

  if ( s == NULL )
    {
      EH(-1, "Problem opening coordinate file.");
    }

  fprintf(s, "%% coordinates of verteces for inertial partition\n");
  fprintf(s, "%% \n");
  fprintf(s, "%% Number dimensions = %d\n", x->num_dim);
  fprintf(s, "%% Number nodes      = %d\n", x->num_nodes);
  fprintf(s, "%% \n");

#ifdef DEBUG  
  fprintf(stderr, "Coordinates:\n");

  for ( i=0; i<x->num_nodes; i++)
    {
      fprintf(stderr, "\t%g\t%g\n", x->x_coord[i], x->y_coord[i]);
    }
  
  fprintf(stderr, "This is a %d dimensional finite element mesh.\n", 
	  x->num_dim);
#endif

  switch ( x->num_dim )
    {
    case 1:
      for ( i=0; i<x->num_nodes; i++)
	{
	  fprintf(s, "%g\n", x->x_coord[i]);
	}
      break;
	  
    case 2:
      for ( i=0; i<x->num_nodes; i++)
	{
	  fprintf(s, "%g %g\n", x->x_coord[i], x->y_coord[i]);
	}
      break;
      
    case 3:
      for ( i=0; i<x->num_nodes; i++)
	{
	  fprintf(s, "%g %g %g\n", x->x_coord[i], x->y_coord[i], x->z_coord[i]);
	}
      break;
      
    default:
      EH(-1, "Bad num_dim in EXODUSII FE database.");
      break;
    }
  err = fclose(s);

  if ( err != 0 ) exit(2);
}

