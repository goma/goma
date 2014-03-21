/************************************************************************ *
 * Copyright (c) 2013 Sandia Corporation. All rights reserved.            *
 *                                                                        *
 * Under the terms of Contract DE-AC04-94AL85000, there is a              *
 * non-exclusive license for use of this work by or on behalf of the      *
 * U.S. Government.                                                       *
 *                                                                        *
 *                                NOTICE:				  * 
 * For five (5) years from 10/03/2013, the United States Government is 	  *
 * granted for itself and others acting on its behalf a paid-up,    	  *
 * nonexclusive, irrevocable worldwide license in this data to reproduce, *
 * prepare derivative works, and perform publicly and display publicly,   *
 * by or on behalf of the Government. There is provision for the possible *
 * extension of the term of this license. Subsequent to that period or    *
 * any extension granted, the United States Government is granted for     *
 * itself and others acting on its behalf a paid-up, nonexclusive,        *
 * irrevocable worldwide license in this data to reproduce, prepare       *
 * derivative works, distribute copies to the public, perform publicly    *
 * and display publicly, and to permit others to do so. The specific      *
 * term of the license can be identified by inquiry made to               *
 * Sandia Corporation or DOE.                                             *
 *	                                                                  *
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT *
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES   *
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL RESPONSIBILITY  *
 * FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION,      *
 * APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE   *
 * WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.                             *
 *
 * Any licensee of this software has the obligation and responsibility to *
 * abide by the applicable export control laws, regulations, and general  *
 * prohibitions relating to the export of technical data. Failure to 	  *
 * obtain an export control license or other authority from the Government*
 *  may result in criminal liability under U.S. laws.  	    	          *
 *	 	                                                          *
 *                             (End of Notice)                            *
 *					         			*
\************************************************************************/
//
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Filename      : gm_cgm_c_interface.C
//
// Purpose       : This file contains the definitions of the functions
//                 to access the CGM functionality through C functions.
//
// Special Notes : These functions are C++ functions and must be compiled
//                 using the C++ compiler
//
// Creator       : Robert A. Kerr
//
// Creation Date : 09/26/2001
//
// Current Owner : Robert A. Kerr
//---------------------------------------------------------------------------

#ifdef USE_CGM

// If CGM is not to be used, do not produce any code for this file.

// ****************** BEGIN STANDARD INCLUDES   *****************************

#include <stdio.h>

// ****************** END STANDARD INCLUDES     *****************************

// ****************** BEGIN GOMA INCLUDES       *****************************

#include <gm_cgm_c_interface.h>

// ****************** END GOMA INCLUDES         *****************************

// ****************** BEGIN CGM INCLUDES      *****************************

#include "GeometryQueryTool.hpp"
#include "GeometryModifyTool.hpp"
#include "AcisQueryEngine.hpp"
#include "AcisModifyEngine.hpp"
#include "VirtualGeometryEngine.hpp"
#include "CubitObserver.hpp"
#include "CubitMessage.hpp"

#include <CastTo.hpp>
#include <Body.hpp>

#include <RefEdge.hpp>
#include <RefFace.hpp>
#include <RefGroup.hpp>
#include <RefVertex.hpp>
#include <RefVolume.hpp>
#include <RefEntityName.hpp>
#include <CompositeTool.hpp>
#include <DLIList.hpp>
#include <MergeTool.hpp>
#include <CubitAttrib.hpp>
#include <RefEntityFactory.hpp>


// ****************** END CGM INCLUDES        *****************************

// ****************** BEGIN ACIS INCLUDES      *****************************

#include <kernel/kernapi/api/api.hxx>
#include <kernel/kernapi/api/api.err>

// ****************** END ACIS INCLUDES      *****************************

// ****************** BEGIN FORWARD DECLARATIONS ****************************
// ****************** END FORWARD DECLARATIONS   *****************************

// ****************** BEGIN STRUCT DECLARATIONS *****************************

// Handles to various CGM Topology and Geometry Entities.
// These Handles are passed through the CGM C Interface.
// MJP Note:
// I think I will need to consolidate these into a smaller number of
// Handles - "base class" handles rather than "specific class" handles.
// The down side of this, of course, is the need for RTTI.

// Topology Entities

// Body Handle
struct BodyHandle_
{
  Body* bodyPtr_;
};

// Volume Handle
struct VolumeHandle_
{
  RefVolume* refVolumePtr_;
};

// Face Handle
struct FaceHandle_
{
  RefFace* refFacePtr_;
};

// Edge Handle
struct EdgeHandle_
{
  RefEdge* refEdgePtr_;
};

// Vertex Handle
struct VertexHandle_
{
  RefVertex* refVertexPtr_;
};

// Geometry Entities

// Plane Handle
// struct PlaneHandle_
// {
//   J_GMT_Plane* planePtr_;
// };

static DLIList<CubitString*> names_to_export;

static const CubitString cgm_version = "7.0.0";

// ****************** END STRUCT DECLARATIONS   *****************************

// ****************** BEGIN EXTERN FUNCTIONS    *****************************

//***********************BEGIN Setup Functions******************************//

extern "C" void gl_cleanup()
{}

int cpp_cgm_initialize()
{
  // This turns off *all* output.  I've turned it on because there
  // doesn't appear to be a more refined way to choose a diagnostic
  // level.
  int flag = 0;
  INFO_FLAG(flag);

  // Initialize the CGM

  CubitObserver::init_static_observers();

  // Initialize the GeometryTool
  
  GeometryQueryTool *gti = GeometryQueryTool::instance();
  if(gti == NULL)
    {
      PRINT_ERROR("GeometryQueryTool failed to instantiate!!!\n");
      return -1;
    }
  GeometryModifyTool *gmi = GeometryModifyTool::instance();
  if(gmi == NULL)
    {
      PRINT_ERROR("GeometryModifyTool failed to instantiate!!!\n");
      return -1;
    }

  if(AcisQueryEngine::instance() == NULL)
    {
      PRINT_ERROR("AcisQueryEngine failed to instantiate!!!\n");
      return -1;
    }
       
  if(AcisModifyEngine::instance() == NULL)
    {
      PRINT_ERROR("AcisModifyEngine failed to instantiate!!!\n");
      return -1;
    }

  if(VirtualGeometryEngine::instance() == NULL)
    {
      PRINT_ERROR("VirtualGeometryEngine failed to instantiate!!!\n");
      return -1;
    }

  // Matt sez: I'm leaving this in for Bob's own debugging purposes.
  // The user's version of this is via an input card "Exported file
  // name = <fname>".
  //_debugExport = 1;
  _debugExport = 0;

  CubitBoolean pass_flag = CUBIT_TRUE;
  CubitAttributeType att_type = CubitAttrib::attrib_type("name");
  //  CubitAttrib::auto_flag(CUBIT_TRUE);
  CubitAttrib::auto_actuate_flag(att_type, pass_flag);
  CubitAttrib::auto_update_flag(att_type, pass_flag);
 
  return(0);
}


//***********************END Setup Functions******************************//


//*************BEGIN Construction and Destruction Functions*****************//

// void cgm_body_delete(BodyHandle* bodyHandle)
// {
//   if ( NULL == bodyHandle )
//   {
//      return;
//   }

//   //
//   // Remove the CGM Body and its associated TopologyEntity's
//   //
//   if ( bodyHandle->bodyPtr_ != NULL )
//   {
//     //
//     // API_BEGIN
//     //
//     // bodyHandle->bodyPtr_->remove();
//     //
//     // API_END
//   }

//   //
//   // Now delete the BodyHandle
//   //
//   delete bodyHandle;
// }

int cgm_geometry_read_SAT( char const* ACISSATFileName,
                           BodyHandle** bodyHandle )
{
  // Debug messages
  //  PRINT_INFO("cgm_geometry_read_SAT called\n");
  
  //  printf( "   ACIS SAT File Name = %s\n", ACISSATFileName );
  
  // Read the Bodies from the SAT file
  FILE* file_ptr = fopen(ACISSATFileName, "r");
  if(!file_ptr)
    {
      PRINT_ERROR("Cannot open file %s\n",ACISSATFileName);
      return (-1);
    }
  
  int id;
  id = RefEntityFactory::instance()->num_ref_volumes();
  CubitStatus status =
    GeometryQueryTool::instance()->import_solid_model( file_ptr,
						       ACISSATFileName,
						       "ACIS_SAT");
  //#ifdef DEBUG
  int num_volumes = id;
  
  id = RefEntityFactory::instance()->num_ref_volumes();
  
  num_volumes -= id;
  num_volumes *= -1;
  //  PRINT_INFO("Number of volumes read in is %d\n",num_volumes);
  
  for(int jj = num_volumes; jj>0; jj--)
    {
      // Although I commented out the "diagnostic" messages, it might still be useful to try
      // to reference the just-read-in volumes. 
      RefVolume* ref_volume =  GeometryQueryTool::instance()->get_ref_volume(id + 1 - jj);
      
//       double volume = ref_volume->measure();
//       PRINT_INFO( "Volume of the input volume  %d (name = %s) = %f\n",
//       		  num_volumes - (jj - 1), ref_volume->entity_name().c_str(), volume);
      
//       CubitVector center = ref_volume->center_point();
//       PRINT_INFO( "Center of the input volume %d (name = %s) = [ %f %f %f ]\n",
//       		  num_volumes - (jj - 1), ref_volume->entity_name().c_str(), center.x(), center.y(), center.z());
    }
  //#endif
  
  // Return success
  return(0);
}

int cgm_plane_construct( double a, double b, double c, double d,
                         PlaneHandle** planeHdl )
{
  PRINT_INFO("cgm_plane_construct not done\n");
  return -1;
  
  //   // Construct a J_GMT_Plane object using the input coefficients
  //   J_GMT_Plane* planePtr = new J_GMT_Plane( a, b, c, d );
  
  //   // Insert it into the Plane Handle
  //   if ( NULL == planePtr )
  //   {
  //     return(-1);
  //   }
  
  //   else
  //   {
  //     // Construct a PlaneHandle
  //     *planeHdl = new PlaneHandle;
  
  //     (*planeHdl)->planePtr_ = planePtr;
  
  //     // Return success
  //     return(0);
  //   }
}

int cgm_vertex_construct( GeometryType refVertexType,
                          char const* name,
                          int dimension,
                          double* coordinates,
                          VertexHandle** vertexHdl )
{
  if(is_name_used(name, 1))
    return (-1);
  //  PRINT_INFO("CGM_vertex_construct called\n");
  
  // Extract the coordinates of the RefVertex
  CubitVector vertexPosition;
  vertexPosition.x( coordinates[0] );
  vertexPosition.y( coordinates[1] );
  
  // Build the appropriately dimensioned position
  if ( dimension == 2 )
    vertexPosition.z(0.0);
  else if ( dimension == 3 )
    vertexPosition.z( coordinates[2] );
  else
    {
      // Oops!
      return(-1);
    }

  RefVertex* refVertexPtr = NULL;

  //   // Construct a RefVertex
  refVertexPtr =
    GeometryModifyTool::instance()->make_RefVertex(UNDEFINED_POINT_TYPE,
						   vertexPosition);
  // Make sure we got one!
  if ( NULL == refVertexPtr )
    {
      PRINT_ERROR("Vertex creation failed!!!!!\n");
      return(-1);
    }

  // Name (label) the RefVertex
  CubitString vertexName = name;
  refVertexPtr->entity_name( vertexName );


  // Stuff it into a VertexHandle structure
  *vertexHdl = new VertexHandle;
  (*vertexHdl)->refVertexPtr_ = refVertexPtr;

  // Return success
  return(0);
}

int is_name_used(char const* name, int print_warning)
{
  RefEntity* temp_entity =
    RefEntityName::instance()->get_refentity( name );
  if(temp_entity != NULL)
    {
      if(print_warning)
	{
	  PRINT_ERROR("Name %s is already used\n",name);
	}
      return 1;
    }
  else
    return 0;
}
  
int cgm_straight_edge_construct( GeometryType refEdgeType,
                                 char const* name,
                                 char const* vertex1Name,
                                 char const* vertex2Name,
                                 EdgeHandle** edgeHdl )
{
  //  PRINT_INFO("CGM_straight_edge_construct called\n");
  // Extract the 2 end RefVertex'es
  if(is_name_used(name, 1))
    return (-1);
  
  RefEntity* temp_entity =
    RefEntityName::instance()->get_refentity( vertex1Name );
  RefVertex* vertex1Ptr = CAST_TO(temp_entity, RefVertex);
  if ( NULL == vertex1Ptr )
    {
      PRINT_ERROR("Couldn't get first RefVertex by name\n");
      return(-1);
    }
  temp_entity = 
    RefEntityName::instance()->get_refentity( vertex2Name );
  RefVertex* vertex2Ptr = CAST_TO(temp_entity, RefVertex);
  if ( NULL == vertex2Ptr )
    {
      PRINT_ERROR("Couldn't get second RefVertex by name\n");
      return(-1);
    }

  RefEdge* refEdgePtr = NULL;

  // Construct a straight RefEdge
  refEdgePtr =
    GeometryModifyTool::instance()->make_RefEdge( refEdgeType,
						  vertex1Ptr,
						  vertex2Ptr );

  // Make sure we got one!
  if ( NULL == refEdgePtr )
    {
      PRINT_ERROR("Straight edge creation failed!!!!!\n");
      return(-1);
  
    }

  // Name (label) the RefEdge
  CubitString refEdgeName = name;
  refEdgePtr->entity_name( refEdgeName );

  // Stuff it into a RefEdge Handle structure
  *edgeHdl = new EdgeHandle;
  (*edgeHdl)->refEdgePtr_ = refEdgePtr;

  //Debug function to export acis sat file of this edge
  if(get_debug_export())
    {
      DLIList<RefEntity*> entity_list;
      entity_list.append(refEdgePtr);
      const char* file_type = "ACIS_SAT";
      refEdgeName += ".sat";
      PRINT_INFO("Exporting %s\n",refEdgeName.c_str());
      int status = GeometryQueryTool::instance()->
	export_solid_model(entity_list,
			   file_type,
			   refEdgeName.c_str(),
			   cgm_version);
    }
  
  // Return success
  return(0);
}

int cgm_quadratic_edge_construct( GeometryType refEdgeType,
                                  char const* name,
                                  char const* vertex1Name,
                                  char const* vertex2Name,
                                  int dimension,
                                  double* coordinates,
                                  CubitSense edgeDirection,
                                  EdgeHandle** edgeHdl )
{
  //  PRINT_INFO("CGM_quadratic_edge_construct called\n");
  // Extract the 2 end RefVertex'es
  if(is_name_used(name, 1))
    return (-1);
  RefEntity* temp_entity =
    RefEntityName::instance()->get_refentity( vertex1Name );
  RefVertex* vertex1Ptr = CAST_TO(temp_entity, RefVertex);
  if ( NULL == vertex1Ptr )
    {
      PRINT_ERROR("Couldn't get first RefVertex by name\n");
      return(-1);
    }
  temp_entity = 
    RefEntityName::instance()->get_refentity( vertex2Name );
  RefVertex* vertex2Ptr = CAST_TO(temp_entity, RefVertex);
  if ( NULL == vertex2Ptr )
    {
      PRINT_ERROR("Couldn't get second RefVertex by name\n");
      return(-1);
    }

  // Extract the coordinates of the intermediate point
  CubitVector point;
  point.x( coordinates[0] );
  point.y( coordinates[1] );

  // Build the appropriately dimensioned point
  if ( dimension == 2 )
    {
      point.z(0.0);
    }

  else if ( dimension == 3 )
    {
      point.z( coordinates[2] );
    }

  else
    {
      // Oops!
      PRINT_ERROR("Unrecognized dimension %d in cgm_quadratic_edge_construct\n", dimension);
    
      return(-1);
    }

  RefEdge* refEdgePtr = NULL;

  // Construct a straight RefEdge
  refEdgePtr =
    GeometryModifyTool::instance()->make_RefEdge( refEdgeType,
						  vertex1Ptr,
						  vertex2Ptr,
						  &point,
						  edgeDirection );
  // Make sure we got one!
  if ( NULL == refEdgePtr )
    {
      PRINT_ERROR("Straight edge creation failed!!!!!\n");
      return(-1);
    }

  // Name (label) the RefEdge
  CubitString refEdgeName = name;
  refEdgePtr->entity_name( refEdgeName );

  // Stuff it into a RefEdge Handle structure
  *edgeHdl = new EdgeHandle;
  (*edgeHdl)->refEdgePtr_ = refEdgePtr;

  //Debug function to export acis sat file of this edge
  if(get_debug_export())
    {
      DLIList<RefEntity*> entity_list;
      entity_list.append(refEdgePtr);
      const char* file_type = "ACIS_SAT";
      refEdgeName += ".sat";
      PRINT_INFO("Exporting %s\n",refEdgeName.c_str());
      int status = GeometryQueryTool::instance()->
	export_solid_model(entity_list,
			   file_type,
			   refEdgeName.c_str(),
			   cgm_version);
    }
  
  // Return success
  return(0);
}

int cgm_composite_edge_construct( char const* name,
                                  int numEdges,
                                  char const** edgeNameList,
                                  EdgeHandle** edgeHdl )
{
  if(is_name_used(name, 1))
    return (-1);
  // Extract the constituent RefEdges
  DLIList<RefEdge*> refEdgeList;
  DLIList<RefVertex*> refVertList;
  
  RefEdge* refEdgePtr = NULL;

  RefEntityName* ren = RefEntityName::instance();
  
  for (int i = 0; i < numEdges; i++)
    {
      // Get the RefEdge given its name
      RefEntity* temp_entity =
	ren->get_refentity( edgeNameList[i] );

      refEdgePtr = CAST_TO(temp_entity, RefEdge);
    
      // Make sure we have a valid RefEdge
      if ( NULL == refEdgePtr )
	{
	  PRINT_ERROR("NULL refEdgePtr in CGM_composite_edge_construct!!!\n");
	  return(-1);
	}
      //else
      //PRINT_INFO("Found edge %d named %s for composite\n",
      //refEdgePtr->id(),refEdgePtr->entity_name().c_str() );
    
    
      refVertList.append(refEdgePtr->startRefVertex());
      refVertList.append(refEdgePtr->endRefVertex());
    
      //PRINT_INFO("Endpoints are:\n");
      //refEdgePtr->startRefVertex()->coordinates().print_me();
      //PRINT_INFO(" and\n");
      //refEdgePtr->endRefVertex()->coordinates().print_me();

      // Insert it into the list
      refEdgeList.append( refEdgePtr );
    }

  RefEdge* crePtr = NULL;

  //Must merge vertices first:
  //PRINT_INFO("Trying to merge vertices # ");
  //for(int kk = refVertList.size(); kk>0; kk--)
  //PRINT_INFO("%d ",refVertList.get_and_step()->id());
  
  MergeTool::instance()->merge_refvertices( refVertList, CUBIT_FALSE, CUBIT_FALSE);
  //PRINT_INFO("\nFinished merging vertices\n");
  
  // Now that we have the list of RefEdges, construct a CompositeRefEdge
  crePtr = CompositeTool::instance()->composite( refEdgeList );

  if ( NULL == crePtr )
    {
      PRINT_ERROR("NULL composite RefEdge in CGM_composite_edge_construct!!!\n");
      return (-1);
    }

  // Name (label) the CompositeRefEdge
  CubitString refEdgeName = name;
  crePtr->entity_name( refEdgeName );

  // Stuff it into a RefEdge Handle structure
  *edgeHdl = new EdgeHandle;
  (*edgeHdl)->refEdgePtr_ = crePtr;

  // Debug function to export acis sat file of this edge
  if(get_debug_export())
    {
      DLIList<RefEntity*> entity_list;
      entity_list.append(crePtr);
      const char* file_type = "ACIS_SAT";
      refEdgeName += ".sat";
      PRINT_INFO("Exporting composite edge %s\n",refEdgeName.c_str());
      int status = GeometryQueryTool::instance()->
	export_solid_model(entity_list,
			   file_type,
			   refEdgeName.c_str(),
			   cgm_version);
    }

  // Return success
  return(0);
}


int cgm_face_polygon_construct_by_coords( GeometryType refFaceType,
                                          char const * name,
                                          int num_verts,
                                          double * double_list,
                                          FaceHandle** faceHandle)
{
  if(is_name_used(name, 1))
    return (-1);
  //  PRINT_INFO("Into cgm_face_polygon_construct_by_coords\n");
  
  // Okay, we should have an array of num_coords doubles for our
  // coordinates.  First we must make vertices, then send it to the
  // other calls
  GeometryType vert_type = UNDEFINED_POINT_TYPE;
  VertexHandle** verthandles = new VertexHandle*[num_verts];
  int ii;
  CubitString vname;
  char** vertNameList = new char* [num_verts];
  
  for(ii = 0; ii<num_verts; ii++)
    {
      //double dime[2];
      //dime[0] = double_list[ii * 2];
      //dime[1] = double_list[ii * 2 + 1];
    
      vname = name;
      vname += "vert";
      vname += ii;
      int rv = cgm_vertex_construct(vert_type, vname.c_str(), 3,
				    &double_list[3*ii],
				    &verthandles[ii]);
      if(rv == 0)
	{
	  //        PRINT_INFO("created vertex %s with coordinates %f %f 0\n",
	  //                     vname.c_str(), dime[0], dime[1]);
	}
    
      vertNameList[ii] = new char [strlen(vname.c_str()) +1];
      strcpy(vertNameList[ii], vname.c_str());
    }
  int rv = cgm_face_polygon_construct_by_verts(refFaceType, name, num_verts,
                                               (char const **)vertNameList, faceHandle);
  delete [] verthandles;
  for(ii = 0; ii< num_verts; ii++)
    delete [] vertNameList[ii];
  delete [] vertNameList;
  
     
  if(rv == 0)
    {
      //    PRINT_INFO("Created surf %s\n",name);
      return 0;
    }
  else
    {
      PRINT_ERROR("Failed to create surf %s\n",name);
      return -1;
    }
}


int cgm_face_polygon_construct_by_verts( GeometryType refFaceType,
                                         char const* name,
                                         int num_verts,
                                         char const** vertNameList,
                                         FaceHandle** faceHandle)
{
  if(is_name_used(name, 1))
    return (-1);
  // PRINT_INFO("Into cgm_face_polygon_construct_by_verts\n");
  
  // Okay, we should have an array of vertices.  First we make the
  // edges, then call the other function
  EdgeHandle** edgehandles = new EdgeHandle*[num_verts];
  CubitString cname;
  char** edgeNameList = new char* [num_verts];
  int ii;
  
  for(ii = 0; ii<num_verts; ii++)
    {
      cname = name;
      cname += "curve";
      cname += ii;
      int rv = cgm_straight_edge_construct(STRAIGHT_CURVE_TYPE,
					   cname.c_str(),
					   vertNameList[ii],
					   vertNameList[(ii + 1) % num_verts],
					   &edgehandles[ii]);
      if(rv == 0)
	{
	  //  PRINT_INFO("created curve %s with verts named %s %s\n",
	  //       cname.c_str(),vertNameList[ii], vertNameList[(ii + 1) % num_verts]);
	}
    
      edgeNameList[ii] = new char [strlen(cname.c_str()) + 1];
    
      strcpy(edgeNameList[ii], cname.c_str());
  
    }
  int rv = cgm_face_construct_by_edges(refFaceType, name, num_verts,
                                       (char const **)edgeNameList, faceHandle);
  delete [] edgehandles;
  for(ii = 0; ii< num_verts; ii++)
    delete [] edgeNameList[ii];
  delete [] edgeNameList;
  
     
  if(rv == 0)
    {
      //    PRINT_INFO("Created surf %s\n",name);
      return 0;
    }
  else
    {
      PRINT_ERROR("Failed to create surf %s\n",name);
      return -1;
    }
}

int cgm_face_disk_construct( GeometryType refFaceType,
                             char const* name,
                             double* coords,
                             double radius,
                             FaceHandle** faceHandle)
{
  if(is_name_used(name, 1))
    return (-1);
  //      PRINT_INFO("Into cgm_face_disk_construct\n");
  //      PRINT_INFO("coords[0] = %f coords[1] = %f radius = %f\n",coords[0],coords[1],radius);
  
  // Okay, first we creat a vertex at the center
  CubitVector centerPos(coords[0], coords[1], 0);
  RefVertex* verts[3];
  verts[0] = GeometryModifyTool::instance()->make_RefVertex(UNDEFINED_POINT_TYPE,
                                                            centerPos);
  if(verts[0] != NULL)
    {
      //     PRINT_INFO("Created vertex with coords %f %f %f\n",centerPos.x(), centerPos.y(),
      //                centerPos.z());
    }
  
  
  //Now we create verts at an extreme in the x direction
  CubitVector max_x(coords[0] + radius, coords[1], 0);
  verts[1] = GeometryModifyTool::instance()->make_RefVertex(UNDEFINED_POINT_TYPE,
                                                            max_x);
  if(verts[1] != NULL)
    {
      //     PRINT_INFO("Created vertex with coords %f %f %f\n",verts[1]->coordinates().x(), verts[1]->coordinates().y(), verts[1]->coordinates().z());
    }
  
  //Now we create verts at an extreme in the y direction
  CubitVector min_x(coords[0], coords[1] + radius, 0);
  verts[2] = GeometryModifyTool::instance()->make_RefVertex(UNDEFINED_POINT_TYPE,
                                                            min_x);
  if(verts[2] != NULL)
    {
      //     PRINT_INFO("Created vertex with coords %f %f %f\n",min_x.x(), min_x.y(), min_x.z());
    }
 
  //Now we create curve arc center vert 1 2 3
  RefEdge* new_edge = GeometryModifyTool::instance()->
    create_arc_center_edge(verts[0], verts[1],
			   verts[2], radius, CUBIT_TRUE);
  //   PRINT_INFO("New edge %d has length %f\n",new_edge->id(), new_edge->measure());
  
  DLIList<RefEdge*> edges;
  edges.append(new_edge);
      
  //Now we create curve arc center vert 1 3 2
  //   new_edge = GeometryModifyTool::instance()->
  //      create_arc_center_edge(verts[0], verts[2],
  //                             verts[1], radius, CUBIT_FALSE );
  // //   PRINT_INFO("New edge %d has length %f\n",new_edge->id(), new_edge->measure());
  //   edges.append(new_edge);

  Body* bodyPtr = GeometryModifyTool::instance()->make_Body( refFaceType,
                                                             edges );
  if ( bodyPtr != NULL)
    {
      DLIList<RefFace*> face_list;
      RefFace *new_face;
      bodyPtr->ref_faces( face_list );
      new_face = face_list.get();
      if ( new_face != NULL )
	{
	  PRINT_INFO("\nCreation of %s (Surface %d) Successful.\n"
		     "This is sheet body %d (%s)\n",
		     new_face->entity_name().c_str(),
		     new_face->id(),
		     bodyPtr->id(),
		     bodyPtr->entity_name().c_str());
	}
      else
	PRINT_WARNING(">>>WARNING:Surface creation was unsuccessful.<<<\n");
    }
   
  //Now we regularize the body
  if(bodyPtr == NULL)
    {
      PRINT_ERROR("Failed to create DISK surface\n");
      return (-1);
    }
  
  //   Body* new_body = NULL;
  //   GeometryModifyTool::instance()->regularize_body( bodyPtr, new_body );

  //   if( new_body )
  //   {
  DLIList<RefFace*> owned_surfs;
  bodyPtr->ref_faces(owned_surfs);
  if(owned_surfs.size() >1)
    {
      PRINT_ERROR("More than one surface created!!!!!\n");
      return (-1);
    }
    
  RefFace* facePtr = owned_surfs.get();
    
  if ( NULL == facePtr )
    {
      PRINT_ERROR("Face creation failed!!!!!!!\n");
      return(-1);
    }
    
  CubitString refFaceName = name;
  facePtr->entity_name( refFaceName );

  *faceHandle = new FaceHandle;
  (*faceHandle)->refFacePtr_ = facePtr;
    
  //Debug function to export acis sat file of this edge
  if(get_debug_export())
    {
      DLIList<RefEntity*> entity_list;
      entity_list.append(bodyPtr);
      const char* file_type = "ACIS_SAT";
      refFaceName += ".sat";
      PRINT_INFO("Exporting %s\n",refFaceName.c_str());
      int status = GeometryQueryTool::instance()->
	export_solid_model(entity_list,
			   file_type,
			   refFaceName.c_str(),
			   cgm_version);
    }
    
  return (0);
    
  //   }
  //   else
  //   {
  //     PRINT_ERROR("Couldn't regularize new body\n");
  //     return(-1);
  //   }

}

int cgm_face_construct_by_edges( GeometryType refFaceType,
                                 char const* name,
                                 int numEdges,
                                 char const** edgeNameList,
                                 FaceHandle** faceHdl )
{
  
  if(is_name_used(name, 1))
    return (-1);
  DLIList<RefEdge*> refEdgeList;
  RefEdge* refEdgePtr = NULL;

  for (int i = 0; i < numEdges; i++)
    {
      // Get the RefEdge given its name
      RefEntity* temp_entity =
	RefEntityName::instance()->get_refentity( edgeNameList[i] );
      refEdgePtr = CAST_TO(temp_entity, RefEdge);

      // Make sure we have a valid RefEdge
      if ( NULL == refEdgePtr )
	{
	  return(-1);
	}

      // Insert it into the list
      refEdgeList.append( refEdgePtr );
    }

  RefFace* facePtr = NULL;

  //   API_BEGIN

  //   // Now that we have the list of RefEdges, construct a Face
  //   facePtr = J_GMT_API::instance()->
  //                           makeRefFace( refFaceType, refEdgeList );

  Body* bodyPtr = GeometryModifyTool::instance()->make_Body( refFaceType,
                                                             refEdgeList );

  DLIList<RefFace*> owned_surfs;
  bodyPtr->ref_faces(owned_surfs);
  if(owned_surfs.size() >1)
    {
      PRINT_ERROR("More than one surface created!!!!!\n");
      return (-1);
    }
  
  facePtr = owned_surfs.get();
    
  if ( NULL == facePtr )
    {
      PRINT_ERROR("Face creation failed!!!!!!!\n");
      return(-1);
    }

  CubitString refFaceName = name;
  facePtr->entity_name( refFaceName );

  //   API_END

  // Stuff it into a RefFace Handle structure
  *faceHdl = new FaceHandle;
  (*faceHdl)->refFacePtr_ = facePtr;

  //Debug function to export acis sat file of this edge
  if(get_debug_export())
    {
      DLIList<RefEntity*> entity_list;
      entity_list.append(bodyPtr);
      const char* file_type = "ACIS_SAT";
      refFaceName += ".sat";
      PRINT_INFO("Exporting %s\n",refFaceName.c_str());
      int status = GeometryQueryTool::instance()->
	export_solid_model(entity_list,
			   file_type,
			   refFaceName.c_str(),
			   cgm_version);
    }
  // Return success
  return(0);
}

int cgm_single_face_body_construct(
				   char const* name,
				   char const* singleFaceName,
				   BodyHandle** bodyHdl )
{
  PRINT_INFO("cgm_single_face_body_construct not done\n");
  return -1;
  
  //   // Extract the single constituent RefFace, given its name
  //   Model* modelPtr = Model::instance();
  //   J_GMT_RefFace* refFacePtr =
  //                      modelPtr->get_ref_face( singleFaceName );

  //   // Make sure we have a valid RefFace
  //   if ( NULL == refFacePtr )
  //   {
  //     return(-1);
  //   }

  //   J_GMT_Body* bodyPtr = NULL;

  //   API_BEGIN

  //   // Now that we have the RefFace, construct a RefVolume
  //   DLJ_GMT_RefFaceList refFaceList;
  //   refFaceList.append(refFacePtr);
  //   J_GMT_RefVolume* refVolumePtr = J_GMT_API::instance()->
  //                           makeRefVolume( NORMAL_LUMP_TYPE, refFaceList );
  //   // Make sure we have a valid RefVolume
  //   if ( NULL == refVolumePtr )
  //   {
  //     return(-1);
  //   }

  //   // Now construct the Body
  //   DLJ_GMT_RefVolumeList refVolumeList;
  //   refVolumeList.append(refVolumePtr);
  //   bodyPtr = J_GMT_API::instance()->makeBody( refVolumeList );

  //   // Make sure we have a valid Body
  //   if ( NULL == bodyPtr )
  //   {
  //     result = outcome(API_FAILED);
  //   }

  //   // Name (label) the Body
  //   CmtString refBodyName = name;
  //   bodyPtr->setLabel( refBodyName );

  //   API_END

  //   if ( NULL == bodyPtr )
  //   {
  //     return(-1);
  //   }

  //   // Stuff it into a Body Handle structure
  //   *bodyHdl = new BodyHandle;
  //   (*bodyHdl)->bodyPtr_ = bodyPtr;

  //   // Return success
  //   return(0);
}

int cgm_multi_face_body_construct( char const* name,
                                   int numFaces,
                                   char const** faceNameArray,
                                   BodyHandle** bodyHdl )
{
  PRINT_INFO("cgm_multi_face_body_construct not done\n");
  return -1;
  
  //   DLJ_GMT_RefVolumeList refVolumeList;
  //   J_GMT_Body* bodyPtr = NULL;

  //   API_BEGIN

  //   // Extract the RefFaces in the input array and construct a RefVolume for each
  //   // one
  //   Model* modelPtr = Model::instance();
  //   for ( int i = 0; i < numFaces; i++ )
  //   {
  //     // Extract a RefFace, given its name
  //     J_GMT_RefFace* refFacePtr =
  //                        modelPtr->get_ref_face( faceNameArray[i] );

  //     // Make sure we have a valid RefFace
  //     if ( NULL == refFacePtr )
  //     {
  //       return(-1);
  //     }

  //     // Construct a RefVolume using this RefFace
  //     DLJ_GMT_RefFaceList refFaceList;
  //     refFaceList.append(refFacePtr);
  //     J_GMT_RefVolume* refVolumePtr = J_GMT_API::instance()->
  //                             makeRefVolume( NORMAL_LUMP_TYPE, refFaceList );

  //     // Make sure we have a valid RefVolume
  //     if ( NULL == refVolumePtr )
  //     {
  //       return(-1);
  //     }

  //     // Insert the new RefVolume into the list
  //     refVolumeList.append(refVolumePtr);
  //   }

  //   // Now construct the Body from the above, newly created RefVolumes
  //   bodyPtr = J_GMT_API::instance()->makeBody( refVolumeList );

  //   // Make sure we have a valid Body
  //   if ( NULL == bodyPtr )
  //   {
  //     result = outcome(API_FAILED);
  //   }
  //   else
  //   {
  //     // Name (label) the Body
  //     CmtString refBodyName = name;
  //     bodyPtr->setLabel( refBodyName );
  //   }

  //   API_END

  //   if ( NULL == bodyPtr )
  //   {
  //     return(-1);
  //   }

  //   // Stuff it into a Body Handle structure
  //   *bodyHdl = new BodyHandle;
  //   (*bodyHdl)->bodyPtr_ = bodyPtr;

  //   // Return success
  //   return(0);
}

// //*************END Construction and Destruction Functions*****************//

// //**********************BEGIN Query Functions*****************************//

int cgm_get_vertex_by_name( char const* name, VertexHandle** vertexHdl )
{
  // Get the RefVertex from the geometry container class.
  RefEntity* temp_entity =
    RefEntityName::instance()->get_refentity( name );
  RefVertex* refVertexPtr = CAST_TO(temp_entity, RefVertex);
  
  if ( NULL == refVertexPtr )
    {
      return(-1);
    }

  // Stuff it into a VertexHandle structure
  *vertexHdl = new VertexHandle;
  (*vertexHdl)->refVertexPtr_ = refVertexPtr;

  // Return success
  return(0);
}

int cgm_get_edge_by_name( char const* name, EdgeHandle** edgeHdl )
{
  //  PRINT_INFO("CGM_get_edge_by_name called\n");
  
  // Get the RefEdge from the geometry container class.
  RefEntity* temp_entity =
    RefEntityName::instance()->get_refentity( name );
   
  RefEdge* refEdgePtr = CAST_TO(temp_entity, RefEdge);
  
  if ( NULL == refEdgePtr )
    {
      //    PRINT_ERROR("Couldn't get edge by name!!!\n");
      return(-1);
    }

  // Stuff it into a EdgeHandle structure
  *edgeHdl = new EdgeHandle;
  (*edgeHdl)->refEdgePtr_ = refEdgePtr;

  // Return success
  return(0);
}

int cgm_get_volume_by_name( char const* name, VolumeHandle** volumeHdl )
{
  // Get the volume from the geometry container class.
  RefEntity* temp_entity =
    RefEntityName::instance()->get_refentity( name );
  RefVolume* refVolumePtr = CAST_TO(temp_entity, RefVolume);
  if ( NULL == refVolumePtr )
    {
      return(-1);
    }

  // Stuff it into a VolumeHandle structure
  *volumeHdl = new VolumeHandle;
  (*volumeHdl)->refVolumePtr_ = refVolumePtr;

  // Return success
  return(0);
}

int cgm_get_face_by_name(char const * name,
                         FaceHandle** faceHdl)
{
  // Get the RefEdge from the geometry container class.
  RefEntity * temp_entity =
    RefEntityName::instance()->get_refentity(name);
   
  RefFace * refFacePtr = CAST_TO(temp_entity, RefFace);
  
  if(NULL == refFacePtr)
    return(-1);

  // Stuff it into a FaceHandle structure
  *faceHdl = new FaceHandle;
  (*faceHdl)->refFacePtr_ = refFacePtr;

  // Return success
  return(0);
}

int cgm_get_body_by_name( char const* name, BodyHandle** bodyHdl )
{
  return (-1);
  
  //     // Get the Body from the geometry container class.
  //   RefEntity* temp_entity =
  //      RefEntityName::instance()->get_refentity( name );
   
  //   RefEdge* refEdgePtr = CAST_TO(temp_entity, RefEdge);
  
  //   if ( NULL == refEdgePtr )
  //   {
  // //    PRINT_ERROR("Couldn't get edge by name!!!\n");
  //     return(-1);
  //   }

  //   // Stuff it into a EdgeHandle structure
  //   *edgeHdl = new EdgeHandle;
  //   (*edgeHdl)->refEdgePtr_ = refEdgePtr;

  //   // Return success
  //   return(0);

}

// //**************BEGIN Vertex Query Functions**************/

// char const* cgm_vertex_get_name( VertexHandle const* entityHdl )
// {
//   // Get the Entity
//   Result status;
//   J_GMT_RefVertex* refEntityPtr = entityHdl->refVertexPtr_;

//   // Return its name
//   return refEntityPtr->getLabel();
// }

int cgm_vertex_get_coordinates( VertexHandle const* vertexHdl,
                                int* num_coords,
                                double** coordinates )
{
  PRINT_ERROR("cgm_vertex_get_coordinates not done\n");
  return -1;
  
  //   // Get the RefVertex
  //   Result status;
  //   J_GMT_RefVertex* refVertexPtr = vertexHdl->refVertexPtr_;

  //   // Get the coordinates
  //   CubitVector coords;
  //   coords = refVertexPtr->coordinates();

  //   // Allocate space for the coordinates
  //   *num_coords = 3;
  //   *coordinates = new double[*num_coords];
  //   if ( NULL == *coordinates )
  //   {
  //     return(-1);
  //   }
  //   (*coordinates)[0] = coords.x();
  //   (*coordinates)[1] = coords.y();
  //   (*coordinates)[2] = coords.z();
}

// //**************END Vertex Query Functions**************/

// //**************BEGIN Edge Query Functions**************/

// char const* cgm_edge_get_name( EdgeHandle const* entityHdl )
// {
//   // Get the Entity
//   Result status;
//   J_GMT_RefEdge* refEntityPtr = entityHdl->refEdgePtr_;

//   // Return its name
//   return refEntityPtr->getLabel();
// }

int cgm_edge_get_closest_point ( EdgeHandle const* edgeHdl,
                                 double point_X,
                                 double point_Y,
                                 double point_Z,
                                 double* closest_point_X,
                                 double* closest_point_Y,
                                 double* closest_point_Z,
                                 double* tangent_vector,
                                 double* distance )
{
  // Get the Edge
  RefEdge* refEdgePtr = edgeHdl->refEdgePtr_;
  assert( refEdgePtr != NULL );

  // Put the input arguments into CGM objects
  CubitVector inputPoint( point_X, point_Y, point_Z );

  CubitVector closestPoint;
  CubitVector tangentVector;

  // Let the RefEdge do the hard work :-)
  CubitStatus result = refEdgePtr->closest_point( inputPoint,
                                                  closestPoint,
                                                  &tangentVector );
  if ( result == CUBIT_FAILURE )
    {
      return(-1);
    }

  // Compute the distance
  *distance = inputPoint.distance_between(closestPoint );

  // Set the closest point and tangent vector output arguments
  if ( tangent_vector != NULL )
    {
      tangent_vector[0] = tangentVector.x();
      tangent_vector[1] = tangentVector.y();
      tangent_vector[2] = tangentVector.z();
    }

  *closest_point_X = closestPoint.x();
  *closest_point_Y = closestPoint.y();
  *closest_point_Z = closestPoint.z();

  // Simple!
  return(0);
}

int cgm_edge_get_closest_point_trimmed ( EdgeHandle const* edgeHdl,
                                         double point_X,
                                         double point_Y,
                                         double point_Z,
                                         double* closest_point_X,
                                         double* closest_point_Y,
                                         double* closest_point_Z,
                                         double* tangent_vector,
                                         double* distance )
{
  
  // Get the Edge
  RefEdge* refEdgePtr = edgeHdl->refEdgePtr_;
  assert( refEdgePtr != NULL );

  // Put the input arguments into CGM objects
  CubitVector inputPoint( point_X, point_Y, point_Z );

  CubitVector closestPoint;
  CubitVector tangentVector;

  // Let the RefEdge do the hard work :-)
  refEdgePtr->closest_point_trimmed( inputPoint, closestPoint );

  // Compute the distance
  *distance = inputPoint.distance_between(closestPoint );

  // Get the tangent vector at this point
  if ( tangent_vector != NULL )
    {
      refEdgePtr->tangent( closestPoint, tangentVector );
      tangent_vector[0] = tangentVector.x();
      tangent_vector[1] = tangentVector.y();
      tangent_vector[2] = tangentVector.z();
    }

  // Set the closest point output arguments
  *closest_point_X = closestPoint.x();
  *closest_point_Y = closestPoint.y();
  *closest_point_Z = closestPoint.z();

  // Simple!
  return(0);
}

int cgm_edge_evaluate_trimmed ( EdgeHandle const* edgeHdl,
                                double inputCoordinate1,
                                unsigned int coord1Index,
                                double inputCoordinate2,
                                unsigned int coord2Index,
                                int* numCoord3Found,
                                double** coord3Array )
{
  PRINT_INFO("cgm_edge_evaluate_trimmed not done\n");
  return -1;
  
  //   // Get the Edge
  //   J_GMT_RefEdge* refEdgePtr = edgeHdl->refEdgePtr_;
  //   assert( refEdgePtr != NULL );

  //   // Make sure some of the other input parameters are reasonable :-)
  //   if ( coord1Index > 2 || coord2Index > 2 )
  //   {
  //     return(-1);
  //   }

  //   Result status;
  //   CmtSet<CubitVector> intersectionPointSet;

  //   status = J_GMT_API::instance()->
  //                    getEdgeLineIntersection( *refEdgePtr,
  //                                             inputCoordinate1,
  //                                             coord1Index,
  //                                             inputCoordinate2,
  //                                             coord2Index,
  //                                             intersectionPointSet );

  //   // Number of intersection points found
  //   *numCoord3Found = intersectionPointSet.size();

  //   // No intersection points were found (on the Edge) that match the input coordinates.
  //   if ( *numCoord3Found == 0 )
  //   {
  //     return(-1);
  //   }

  //   // Allocate memory for the array of intersection coordinates that is to be returned
  //   *coord3Array = new double[(unsigned int)(*numCoord3Found)];
  //   if ( NULL == *coord3Array )
  //   {
  //     assert( *coord3Array != NULL );
  //   }

  //   // Return the array of "3rd" coordinates of the intersection points that were found
  //   int index = 3 - (coord1Index + coord2Index);
  //   int counter = 0;
  //   CmtSetIterator<CubitVector> pointIterator(intersectionPointSet);
  //   for ( pointIterator.startIteration();
  //         pointIterator.anyMore();
  //         pointIterator.stepForward())
  //   {
  //     // Get the Point
  //     CubitVector intersectionPoint = pointIterator.getCurrentItem();

  //     // Extract the "3rd" coordinate and insert into the output array
  //     // MJP Note: There needs to be a member function of CubitVector that extracts
  //     //           one of its components given an index value.
  //     if ( index == 0 )
  //     {
  //       (*coord3Array)[counter] = intersectionPoint.x();
  //     }

  //     else if ( index == 1 )
  //     {
  //       (*coord3Array)[counter] = intersectionPoint.y();
  //     }

  //     else
  //     {
  //       (*coord3Array)[counter] = intersectionPoint.z();
  //     }

  //     counter++;
  //   }

  //   // Done!
  //   return(0);
}

int cgm_edge_get_vertices( EdgeHandle const* edgeHdl,
                           VertexHandle** startVertexHandle,
                           VertexHandle** endVertexHandle )
{
  PRINT_ERROR("cgm_edge_get_vertices not done\n");
  return -1;
  
  //   // Extract the RefEdge from the EdgeHandle
  //   J_GMT_RefEdge* refEdgePtr = edgeHdl->refEdgePtr_;

  //   // Get the Start and End Vertices of the RefEdge
  //   J_GMT_RefVertex* startRefVertexPtr =
  //                       refEdgePtr->startRefVertex();

  //   J_GMT_RefVertex* endRefVertexPtr =
  //                       refEdgePtr->endRefVertex();

  //   // Make sure we have them both
  //   if ( NULL == startRefVertexPtr || NULL == endRefVertexPtr )
  //   {
  //     return (-1);
  //   }

  //   // Construct the VertexHandles for the start and end RefVertex objects
  //   *startVertexHandle = new VertexHandle;
  //   (*startVertexHandle)->refVertexPtr_ = startRefVertexPtr;
  //   *endVertexHandle = new VertexHandle;
  //   (*endVertexHandle)->refVertexPtr_ = endRefVertexPtr;

  //   // Make sure the Handles were appropriately constructed
  //   if ( NULL != *startVertexHandle && NULL != *endVertexHandle )
  //   {
  //     return (0);
  //   }
  //   else
  //   {
  //     return (-1);
  //   }
}

// //**************END Edge Query Functions**************/

// //**************BEGIN Face Query Functions**************/

// char const* cgm_face_get_name( FaceHandle const* entityHdl )
// {
//   // Get the Entity
//   Result status;
//   J_GMT_RefFace* refEntityPtr = entityHdl->refFacePtr_;

//   // Return its name
//   return refEntityPtr->getLabel();
// }

// void cgm_face_get_closest_point_trimmed (
//                                       FaceHandle const* faceHdl,
//                                       double point_X,
//                                       double point_Y,
//                                       double point_Z,
//                                       double* closest_point_X,
//                                       double* closest_point_Y,
//                                       double* closest_point_Z )
// {
//   printf( "Function: cgm_face_get_closest_point_trimmed - not implemented.\n" );
//   return;
// }

int cgm_named_face_get_closest_boundary( char const* name, 
                                         double point_x,
                                         double point_y,
                                         double point_z,
                                         double* distance)
{
  FaceHandle* temp_face;
  if(cgm_get_face_by_name( name, &temp_face) == 0)
    return((cgm_face_get_closest_boundary(temp_face, point_x, point_y,
					  point_z, distance)));
  else
    return -1;
}

int cgm_face_get_closest_boundary( FaceHandle const* faceHdl,
                                   double point_x,
                                   double point_y,
                                   double point_z,
                                   double* distance)
{
  //Okay, first we need to get distance, then signedness (or vice
  //versa would be fine, but that's the way I want to deal with it

  //Okay, first we need the boundaries
  RefFace* ref_face = faceHdl->refFacePtr_;
  //PRINT_INFO("Got ref face %d from faceHandle\n",ref_face->id());
  DLIList<RefEdge*> edges;
  ref_face->ref_edges(edges);
  if(edges.size() < 1)
    {
      PRINT_ERROR("No refEdges in refFace %d\n", ref_face->id());
      return -1;
    }
  double closest_distance = CUBIT_DBL_MAX;
  CubitVector location(point_x, point_y, point_z);
  CubitVector new_location;
  
  for(int ii = edges.size(); ii>0; ii--)
    {
      RefEdge* this_edge = edges.get_and_step();
      if(this_edge->closest_point_trimmed(location, new_location) != CUBIT_SUCCESS)
	{
	  PRINT_ERROR("Error in closest_point_trimmed\n");
	  return -1;
	}
      closest_distance = CUBIT_MIN(closest_distance,
				   new_location.distance_between(location));
    }
  //Okay, we should have the closest_distance now, let's get signage
  CubitPointContainment is_in = ref_face->point_containment(location);
  switch(is_in)
    {
    case CUBIT_PNT_OUTSIDE:
      closest_distance *= -1.0;
      break;
    case CUBIT_PNT_BOUNDARY:
      closest_distance = 0;
      break;
    case CUBIT_PNT_INSIDE:
      //do nothing, since it is already positive
      break;
    case CUBIT_PNT_UNKNOWN:
      PRINT_ERROR("Unable to classify point for containment\n");
      return (-1);
    }
  
  *distance = closest_distance;
  return 1;
}
  
  
int cgm_named_face_get_closest_boundary_keep(char const * name, // FaceHandle const* faceHdl,
					     double src_point[3],
					     double * distance,
					     double closest[3])
{
  FaceHandle * temp_face;
  CubitVector location_keep;

  if(cgm_get_face_by_name( name, &temp_face))
    return -1;

  //Okay, first we need to get distance, then signedness (or vice
  //versa would be fine, but that's the way I want to deal with it

  //Okay, first we need the boundaries
  RefFace* ref_face = temp_face->refFacePtr_;
  //PRINT_INFO("Got ref face %d from faceHandle\n",ref_face->id());
  DLIList<RefEdge*> edges;
  ref_face->ref_edges(edges);
  if(edges.size() < 1)
    {
      PRINT_ERROR("No refEdges in refFace %d\n", ref_face->id());
      return -1;
    }
  double closest_distance = CUBIT_DBL_MAX;
  CubitVector location = src_point;
  CubitVector new_location;
  
  for(int ii = edges.size(); ii>0; ii--)
    {
      RefEdge* this_edge = edges.get_and_step();
      if(this_edge->closest_point_trimmed(location, new_location) != CUBIT_SUCCESS)
	{
	  PRINT_ERROR("Error in closest_point_trimmed\n");
	  return -1;
	}
      double this_distance = new_location.distance_between(location);
      if(this_distance < closest_distance)
	{
	  closest_distance = this_distance;
	  location_keep = new_location;
	}
    }
  //Okay, we should have the closest_distance now, let's get signage
  CubitPointContainment is_in = ref_face->point_containment(location);
  switch(is_in)
    {
    case CUBIT_PNT_OUTSIDE:
      closest_distance *= -1.0;
      break;
    case CUBIT_PNT_BOUNDARY:
      closest_distance = 0;
      break;
    case CUBIT_PNT_INSIDE:
      //do nothing, since it is already positive
      break;
    case CUBIT_PNT_UNKNOWN:
      PRINT_ERROR("Unable to classify point for containment\n");
      return (-1);
    }
  
  * distance = closest_distance;
  closest[0] = location_keep.x();
  closest[1] = location_keep.y();
  closest[2] = location_keep.z();
  return 1;
}
  
  
int cgm_named_volume_get_closest_boundary( char const* name, 
                                           double point_x,
                                           double point_y,
                                           double point_z,
                                           double* distance)
{
  VolumeHandle* temp_vol;
  if(cgm_get_volume_by_name( name, &temp_vol) == 0)
    return((cgm_volume_get_closest_boundary(temp_vol, point_x, point_y,
					    point_z, distance)));
  else
    return -1;
}
  
int cgm_volume_get_closest_boundary( VolumeHandle const* volHdl,
                                     double point_x,
                                     double point_y,
                                     double point_z,
                                     double* distance)
{
  //Okay, first we need to get distance, then signedness (or vice
  //versa would be fine, but that's the way I want to deal with it

  //Okay, first we need the boundaries
  RefVolume* ref_vol = volHdl->refVolumePtr_;
  //PRINT_INFO("Got ref volume %d from volumeHandle\n",ref_vol->id());
  DLIList<RefFace*> faces;
  ref_vol->ref_faces(faces);
  if(faces.size() < 1)
    {
      PRINT_ERROR("No refFaces in refVolume %d\n", ref_vol->id());
      return -1;
    }
  double closest_distance = CUBIT_DBL_MAX;
  CubitVector location(point_x, point_y, point_z);
  CubitVector new_location;
  
  for(int ii = faces.size(); ii>0; ii--)
    {
      RefFace* this_face = faces.get_and_step();
      this_face->find_closest_point_trimmed(location, new_location);
      closest_distance = CUBIT_MIN(closest_distance,
				   new_location.distance_between(location));
    }
  //Okay, we should have the closest_distance now, let's get signage
  CubitPointContainment is_in = ref_vol->point_containment(location);
  switch(is_in)
    {
    case CUBIT_PNT_OUTSIDE:
      closest_distance *= -1.0;
      break;
    case CUBIT_PNT_BOUNDARY:
      closest_distance = 0;
      break;
    case CUBIT_PNT_INSIDE:
      //do nothing, since it is already positive
      break;
    case CUBIT_PNT_UNKNOWN:
      PRINT_ERROR("Unable to classify point for containment\n");
      return (-1);
    }
  
  *distance = closest_distance;
  return 1;
}


int cgm_named_volume_get_closest_boundary_keep(const char * name,
					       double point[3],
					       double * distance,
					       double closest[3])
{
  VolumeHandle * temp_vol;
  CubitVector location_keep;

  if(cgm_get_volume_by_name( name, &temp_vol))
    return -1;
  
  //Okay, first we need to get distance, then signedness (or vice
  //versa would be fine, but that's the way I want to deal with it

  //Okay, first we need the boundaries
  RefVolume* ref_vol = temp_vol->refVolumePtr_;
  //PRINT_INFO("Got ref volume %d from volumeHandle\n",ref_vol->id());
  DLIList<RefFace*> faces;
  ref_vol->ref_faces(faces);
  if(faces.size() < 1)
    {
      PRINT_ERROR("No refFaces in refVolume %d\n", ref_vol->id());
      return -1;
    }
  double closest_distance = CUBIT_DBL_MAX;
  CubitVector location(point[0], point[1], point[2]);
  CubitVector new_location;
  
  for(int ii = faces.size(); ii>0; ii--)
    {
      RefFace* this_face = faces.get_and_step();
      this_face->find_closest_point_trimmed(location, new_location);
      double this_distance = new_location.distance_between(location);
      if(this_distance < closest_distance)
	{
	  closest_distance = this_distance;
	  location_keep = new_location;
	}
    }
  //Okay, we should have the closest_distance now, let's get signage
  CubitPointContainment is_in = ref_vol->point_containment(location);
  switch(is_in)
    {
    case CUBIT_PNT_OUTSIDE:
      closest_distance *= -1.0;
      break;
    case CUBIT_PNT_BOUNDARY:
      closest_distance = 0;
      break;
    case CUBIT_PNT_INSIDE:
      //do nothing, since it is already positive
      break;
    case CUBIT_PNT_UNKNOWN:
      PRINT_ERROR("Unable to classify point for containment\n");
      return (-1);
    }
  
  *distance = closest_distance;
  closest[0] = location_keep.x();
  closest[1] = location_keep.y();
  closest[2] = location_keep.z();
  
  return 1;
}


// void cgm_face_get_normal ( FaceHandle const* faceHdl,
//                            double point_X,
//                            double point_Y,
//                            double point_Z,
//                            double* closest_point_X,
//                            double* closest_point_Y,
//                            double* closest_point_Z,
//                            double* unit_normal_vector[] )
// {
//   printf( "Function: cgm_face_get_normal - not implemented.\n" );
//   return;
// }

// void cgm_face_get_principal_curvatures (
//                                   FaceHandle const* faceHdl,
//                                   double point_X,
//                                   double point_Y,
//                                   double point_Z,
//                                   double* closest_point_X,
//                                   double* closest_point_Y,
//                                   double* closest_point_Z,
//                                   double* curvature1,
//                                   double* curvature2 )
// {
//   printf( "Function: cgm_face_get_principal_curvatures - not implemented.\n" );
//   return;
// }

// //**************END Face Query Functions**************/

// //**************BEGIN Plane Query Functions**************/

// void cgm_plane_get_shortest_distance_to_point (
//                                       PlaneHandle const* planeHdl,
//                                       double point_X,
//                                       double point_Y,
//                                       double point_Z,
//                                       double* distance )
// {
//   // Construct a Vector that represents the input point
//   CubitVector point( point_X, point_Y, point_Z );

//   // Extract the J_GMT_Plane
//   J_GMT_Plane* planePtr = planeHdl->planePtr_;

//   // Get the distance from the point to the plane
//   *distance = planePtr->distance( point );

//   return;
// }

// void cgm_plane_get_coefficients( PlaneHandle const* planeHdl,
//                                  double* a,
//                                  double* b,
//                                  double* c,
//                                  double* d )
// {
//   // Extract the J_GMT_Plane
//   J_GMT_Plane* planePtr = planeHdl->planePtr_;

//   // Get the coefficients of the plane equation:
//   //    Ax + By + Cz + D = 0

//   // Get A, B and C
//   // Note: The normal is a normalized vector (unit vector)
//   CubitVector normal = planePtr->normal();
//   *a = normal.x();
//   *b = normal.y();
//   *c = normal.z();

//   // Get D
//   *d = planePtr->coefficient();

//   return;
// }

// //**************END Plane Query Functions**************/

// //**************BEGIN Body Query Functions**************/

// char const* cgm_body_get_name( BodyHandle const* entityHdl )
// {
//   // Get the Entity
//   Result status;
//   J_GMT_Body* refEntityPtr = entityHdl->bodyPtr_;

//   // Return its name
//   return refEntityPtr->getLabel();
// }

// int cgm_body_get_volumes( BodyHandle const* bodyHdl,
//                           int* num_volumes,
//                           VolumeHandle** volumeHandles )
// {
//   // Get the Volumes
//   Result status;
//   CmtSet<J_GMT_RefVolume*> refVolumeSet;
//   J_GMT_Body* bodyPtr = bodyHdl->bodyPtr_;
//   J_DBC_RefVolume dbcRefVolume;
//   status = dbcRefVolume.getRefVolumes(*bodyPtr, refVolumeSet);
//   if ( status.ok() == CMT_FAILURE )
//   {
//     return(-1);
//   }

//   // Find out how many RefVolumes there are and allocate space in the
//   // array of VolumeHandles
//   *num_volumes = refVolumeSet.size();
//   if ( *num_volumes == 0 )
//   {
//     return(0);
//   }

//   else
//   {
//     // Allocate space for the array of VolumeHandle pointers that
//     // will be constructed and returned
//     *volumeHandles = new VolumeHandle[*num_volumes];
//     if ( NULL == (*volumeHandles) )
//     {
//       return(-1);
//     }
//   }

//   // Now construct Handles for the returned RefVolumes and return them
//   CmtSetIterator<J_GMT_RefVolume*> refVolumeSetIter(refVolumeSet);

//   int counter = 0;
//   for (refVolumeSetIter.startIteration();
//        refVolumeSetIter.anyMore();
//        refVolumeSetIter.stepForward())
//   {
//     // Get the RefVolume
//     J_GMT_RefVolume* refVolumePtr = refVolumeSetIter.getCurrentItem();

//     // Construct a new VolumeHandle and populate it
//     VolumeHandle* volumeHdl = new VolumeHandle;
//     if ( NULL == volumeHdl )
//     {
//       return(-1);
//     }
//     volumeHdl->refVolumePtr_ = refVolumePtr;

//     // Add it to the VolumeHandle array
//     volumeHandles[counter] = volumeHdl;
//     counter++;
//   }

//   return(0);
// }

// int cgm_body_get_faces( BodyHandle const* bodyHdl,
//                         int* num_faces,
//                         FaceHandle** faceHandles )
// {
//   printf( "Function: cgm_body_get_faces - not implemented.\n" );
//   return(-1);
// }

// int cgm_body_get_edges( BodyHandle const* bodyHdl,
//                         int* num_edges,
//                         EdgeHandle** edgeHandles )
// {
//   printf( "Function: cgm_body_get_edges - not implemented.\n" );
//   return(-1);
// }

// int cgm_body_get_vertices( BodyHandle const* bodyHdl,
//                            int* num_vertices,
//                            VertexHandle** vertexHandles )
// {
//   printf( "Function: cgm_body_get_vertices - not implemented.\n" );
//   return(-1);
// }

// int cgm_body_classify_point ( BodyHandle const* bodyHdl,
//                               double point_X,
//                               double point_Y,
//                               double point_Z )
// {
//   // Get the Body
//   J_GMT_Body* bodyPtr = bodyHdl->bodyPtr_;

//   // Construct a Vector that represents the input point
//   CubitVector point( point_X, point_Y, point_Z );

//   // Classify the point wrt the Body
//   PointClassification classify = J_GMT_API::instance()->
//                                   classifyPoint( point, bodyPtr );

//   int result;
//   if ( classify == CMT_INSIDE )
//   {
//     result = -1;
//   }
//   else if ( classify == CMT_ON_BOUNDARY )
//   {
//     result = 0;
//   }
//   else if ( classify == CMT_OUTSIDE )
//   {
//     result = 1;
//   }
//   else
//   {
//     result = -2;
//   }

//   return result;
// }

// //**************END Body Query Functions**************/

// //**************BEGIN Point Query Functions**************/

// double cgm_point_get_closest_distance_to_body_boundary (
//                                      double point_X,
//                                      double point_Y,
//                                      double point_Z,
//                                      BodyHandle const* bodyHdl,
//                                      double* closest_point_X,
//                                      double* closest_point_Y,
//                                      double* closest_point_Z )
// {
//   // Get the Body
//   J_GMT_Body* bodyPtr = bodyHdl->bodyPtr_;

//   // Construct Vectors that represent the input points
//   CubitVector point( point_X, point_Y, point_Z );
//   CubitVector closestPoint( *closest_point_X, *closest_point_Y, *closest_point_Z );

//   // Get the closest distance (and classify the point)
//   double distance = 0.0;
//   J_GMT_BasicTopologyEntity* hitBTEPtr = NULL;
//   PointClassification classify = J_GMT_API::instance()->
//                                   classifyPoint( point,
//                                                  bodyPtr,
//                                                  distance,
//                                                  closestPoint,
//                                                  hitBTEPtr );

//   return distance;
// }

// //**************END Point Query Functions**************/

// //**********************END Query Functions*****************************//

// //**************BEGIN Body Boolean Functions**************/

int cgm_volume_subtract( char const* name,
                         VolumeHandle const* vol1Hdl,
                         VolumeHandle const* vol2Hdl,
                         VolumeHandle** newVolumeHdl )
{
  // Extract the input volumes
  RefVolume* vol1Ptr = vol1Hdl->refVolumePtr_;
  RefVolume* vol2Ptr = vol2Hdl->refVolumePtr_;
  //convert volumes to Bodies
  DLIList<Body*> bodies;
  DLIList<Body*> new_bodies;

  vol1Ptr->bodies(bodies);

  if ( bodies.size() > 1 )
    {
      PRINT_ERROR("RefVolume %d is owned by more than 1 Body.\n",
		  vol1Ptr->id());
      return (-1);
    }
  Body* bod1Ptr = bodies.get();
  bodies.clean_out();
  vol2Ptr->bodies(bodies);

  if ( bodies.size() > 1 )
    {
      PRINT_ERROR("RefVolume %d is owned by more than 1 Body.\n",
		  vol1Ptr->id());
      return (-1);
    }

  Body* newBodyPtr = NULL;

  // Perform the non-destructive boolean operation on the input Bodies.

  CubitStatus status =
    GeometryModifyTool::instance()->subtract( bod1Ptr,
					      bodies,
					      new_bodies,
					      true );

  if ( status == CUBIT_FAILURE )
    {
      return (-1);
    }

  else
    {
      if(new_bodies.size() > 1)
	{
	  PRINT_ERROR("Subtract produced more than one body\n");
	  return (-1);
	}
      DLIList<RefVolume*> ref_vols;
      new_bodies.get()->ref_volumes(ref_vols);
      if(ref_vols.size() >1)
	{
	  PRINT_ERROR("Subtract produced more than one volume\n");
	  return (-1);
	}
      RefVolume *new_vol = ref_vols.get();
    
      // Name (label) the volume
      CubitString refVolumeName = name;
      new_vol->entity_name( refVolumeName );

      // Stuff it into a Volume Handle structure
      *newVolumeHdl = new VolumeHandle;
      (*newVolumeHdl)->refVolumePtr_ = new_vol;
    }

  // Return success
  return(0);
}

int cgm_volume_unite( char const* name,
                      VolumeHandle const* vol1Hdl,
                      VolumeHandle const* vol2Hdl,
                      VolumeHandle** newVolumeHdl )
{
  // Extract the input volumes
  RefVolume* vol1Ptr = vol1Hdl->refVolumePtr_;
  RefVolume* vol2Ptr = vol2Hdl->refVolumePtr_;
  //convert volumes to Bodies
  DLIList<Body*> bodies;
  DLIList<Body*> new_bodies;

  vol1Ptr->bodies(bodies);

  if ( bodies.size() > 1 )
    {
      PRINT_ERROR("RefVolume %d is owned by more than 1 Body.\n",
		  vol1Ptr->id());
      return (-1);
    }
  Body* bod1Ptr = bodies.get();
  bodies.clean_out();
  vol2Ptr->bodies(bodies);

  if ( bodies.size() > 1 )
    {
      PRINT_ERROR("RefVolume %d is owned by more than 1 Body.\n",
		  vol1Ptr->id());
      return (-1);
    }
  Body* bod2Ptr = bodies.get();
   

  Body* newBodyPtr = NULL;

  // Perform the non-destructive boolean operation on the input Bodies.

  CubitStatus status =
    GeometryModifyTool::instance()->unite( bod1Ptr,
					   bod2Ptr,
					   newBodyPtr,
					   true );

  if ( status == CUBIT_FAILURE )
    {
      return (-1);
    }

  else
    {
      DLIList<RefVolume*> ref_vols;
      newBodyPtr->ref_volumes(ref_vols);
      if(ref_vols.size() >1)
	{
	  PRINT_ERROR("Unite produced more than one volume\n");
	  return (-1);
	}
      RefVolume *new_vol = ref_vols.get();
    
      // Name (label) the volume
      CubitString refVolumeName = name;
      new_vol->entity_name( refVolumeName );

      // Stuff it into a Volume Handle structure
      *newVolumeHdl = new VolumeHandle;
      (*newVolumeHdl)->refVolumePtr_ = new_vol;
    }

  // Return success
  return(0);
}

// //**************END Body Boolean Functions**************/
//********* Begin Helper Functions *********************

int get_debug_export()
{
  return _debugExport;
}

void set_debug_export(int new_state)
{
  _debugExport = new_state;
}

int add_me_to_export(const char* name_of_entity)
{
  if(name_of_entity[0] != '\0')
    {
      names_to_export.append(new CubitString(name_of_entity));
    }
  return names_to_export.size();
}

int all_done_export_now(const char * name_of_export_file)
{
  if(name_of_export_file[0] == '\0')
    {
      PRINT_ERROR("NULL filename in all_done_export_now\n");
      return -1;
    }
  DLIList<RefEntity*> entity_list;

  for(int ii = names_to_export.size(); ii>0; ii--)
    {
      entity_list.append(RefEntityName::instance()->
			 get_refentity(names_to_export.get_and_step()->c_str()));
    
    }
  const char* file_type = "ACIS_SAT";
  int status = GeometryQueryTool::instance()->
    export_solid_model(entity_list,
		       file_type,
		       name_of_export_file,
		       cgm_version);
  return status;
}

void reset_export_list()
{
  for(int ii = names_to_export.size(); ii>0; ii--)
    //delete names_to_export.pop();
    {
      names_to_export.last();
      delete names_to_export.remove();
    }
  names_to_export.reset();
}


// // ****************** END EXTERN FUNCTIONS      *****************************

#endif
