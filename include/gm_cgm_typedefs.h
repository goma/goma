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
 
/***************************************************************************
   Filename      : mm_cgm_typedefs.h
  
   Purpose       : This file contains the typedefs associated with the 
                 : CGM C Interface functions.
  
   Special Notes :
  
   Creator       : Robert A. Kerr
  
   Creation Date : 09/26/2001
  
   Current Owner : Robert A. Kerr
****************************************************************************/

#ifndef GOMA_MM_CGM_TYPEDEFS_H
#define GOMA_MM_CGM_TYPEDEFS_H

/******************** BEGIN STANDARD INCLUDES   ****************************/
/******************** END STANDARD INCLUDES     ****************************/

/******************** BEGIN CGM INCLUDES   ****************************/
#include "CubitDefines.h"
#include "GeometryDefines.h"
/******************** END CGM INCLUDES     ****************************/

/******************** BEGIN STRUCT DECLARATIONS ****************************/

/* Handles to various CGM Topology Entities.
 * These Handles are passed through the CGM C Interface. 
*/

/*
 * Topology Entities
 */
typedef struct BodyHandle_    BodyHandle;
typedef struct VolumeHandle_  VolumeHandle;
typedef struct FaceHandle_    FaceHandle;
typedef struct EdgeHandle_    EdgeHandle;
typedef struct VertexHandle_  VertexHandle;

/*
 * Geometry Entities
 */
typedef struct PlaneHandle_  PlaneHandle;

/******************** END STRUCT DECLARATIONS   ****************************/

/******************** BEGIN EXTERN FUNCTIONS    ****************************/
/******************** END EXTERN FUNCTIONS    ****************************/

#endif

