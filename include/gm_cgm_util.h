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
 
#ifndef _GM_CGM_UTIL_H
#define _GM_CGM_UTIL_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _GM_CGM_UTIL_C
#define EXTERN /* nothing */
#endif

#ifndef _GM_CGM_UTIL_C
#define EXTERN extern /* nothing */
#endif

#include "gm_cgm_typedefs.h"

EXTERN void print_edge_info
(void);

EXTERN void cgm_initialize
(void);

EXTERN void create_cgm_geometry
(void);

EXTERN int read_ACIS_file
(char *fname,		/* name of file to read from. */
       BodyHandle *);		/* CGM body to put it in. */

#define MAX_CGM_INPUT_STRING_LENGTH 40000

extern char cgm_exported_filename[MAX_FNL];
extern char cgm_input_string[MAX_CGM_INPUT_STRING_LENGTH];
extern int cgm_input_string_length; /* Usually a lot smaller than MAX_CGM_INPUT_STRING_LENGTH */
#endif
