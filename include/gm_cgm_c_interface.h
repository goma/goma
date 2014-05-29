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
   Filename      : mm_cgm_c_interface.h

   Purpose       : This file contains the declarations of the functions
                   to access the CGM functionality from C functions.

   Special Notes :

   Creator       : Robert A. Kerr

   Creation Date : 09/26/2001

   Current Owner : Robert A. Kerr
****************************************************************************/

#ifndef _MM_CGM_C_INTERFACE_H
#define _MM_CGM_C_INTERFACE_H

#ifdef USE_CGM

/* If CGM is not to be used, do not process any part of this file. */

/******************** BEGIN STANDARD INCLUDES   ****************************/
/******************** END STANDARD INCLUDES     ****************************/

/******************** BEGIN GOMA INCLUDES   ****************************/

#include "gm_cgm_typedefs.h"

/******************** END GOMA INCLUDES     ****************************/

/******************** BEGIN STRUCT DECLARATIONS ****************************/
/******************** END STRUCT DECLARATIONS   ****************************/

/******************** BEGIN EXTERN FUNCTIONS    ****************************/

/*
 * Use the necessary keywords to make these functions callable and
 * implementable in C++
 */

#ifdef __cplusplus
extern "C"
{
#endif

/***********************BEGIN Setup Functions******************************/

extern int cpp_cgm_initialize();
/*
 *R int
 *R- Return a status: 0 => Success
 *R-                 <0 => Failure
 *- Initializes the CGM library. This includes setting up the ERD for
 *- solid model objects in the CGM.
 *- Errors occur when:
 *-   The ERD Scheme file was not found.
*/

/***********************END Setup Functions******************************/

/*************BEGIN Construction and Destruction Functions*****************/

/* extern void cgm_body_delete(BodyHandle* bodyHdl); */
/*
 *I bodyHdl
 *I- Handle to the ACIS BODY to be deleted
 * This function deletes the CGM Body associated with the input bodyHdl
 * along with all the TopologyEntity's associated with the BODY.
 * The input BodyHandle is also deleted and hence, should not be accessed
 * again.
*/

extern int cgm_geometry_read_SAT( char const* ACISSATFileName,
                                BodyHandle** bodyHdl );
/*
 *R int
 *R- Return a status: 0 => Success
 *R-                 <0 => Failure
 *I ACISSATFileName
 *I- Name of the ACIS SAT file to be read
 *O bodyHdl
 *O- Handle to the ACIS BODY in the SAT file
 *-Function to read a specified ACIS SAT file and return a pointer to
 *-a BodyHandle struct. The C code uses this handle when referring to
 *-this Body.
 *-Errors:
 *-  ACIS SAT file does not exist
 *-  The file does not contain a BODY
 *-  The file contains more than 1 BODY
 *-If an error occurs, The returned BodyHandle is NULL.
 *-Note: The BodyHandle object is constructed within the CGM.
*/

extern int cgm_plane_construct( double a, double b, double c, double d,
                                PlaneHandle** planeHdl );
/*
 *R int
 *R- Return a status: 0 => Success
 *R-                 <0 => Failure
 *I a, b, c, d
 *I- Coefficients of the Plane Equation. The first 3 represent the vector that
 *I- is normal to the plane and the fourth represents the shortest distance from
 *I- the origin to the plane.
 *O planeHdl
 *O- Handle to the newly constructed Plane object.
 *- This function creates a plane object given the 4 coefficients of the equation:
 *-     ax + by + cz + d = 0
*/

extern int cgm_vertex_construct( enum GeometryType geometryType,
                                char const* name,
                                int dimension,
                                double* coordinates,
                                VertexHandle** vertexHdl );
/*
 *R int
 *R- Return a status: 0 => Success
 *R-                 <0 => Failure
 *I GeometryType
 *I- enum describing the type of geometry being constructed.
 *I- Allowed geometry type enums for this function:
 *I-   UNDEFINED_POINT_TYPE
 *I name
 *I- The name or label of the Vertex object. The object can be
 *I- referenced by name after it has been created. A handle to
 *I- this object can also be requested by name.
 *I dimension
 *I- An integer indicating the dimensionality of the geometric space
 *I- (2 for 2D and 3 for 3D).
 *I coordinates
 *I- The array of coordinates
 *O vertexHdl
 *O- Handle to the newly constructed Vertex object.
 *- This function creates a new, named vertex object given the coordinates of
 *- the desired point.
*/

extern int cgm_straight_edge_construct( enum GeometryType geometryType,
                                        char const* name,
                                        char const* vertex1Name,
                                        char const* vertex2Name,
                                        EdgeHandle** edgeHdl );
/*
 *R int
 *R- Return a status: 0 => Success
 *R-                 <0 => Failure
 *I GeometryType
 *I- enum describing the type of geometry being constructed.
 *I- Allowed geometry type enums for this function:
 *I-   STRAIGHT_CURVE_TYPE
 *I name
 *I- The name or label of the Edge object. The object can be
 *I- referenced by name after it has been created. A handle to
 *I- this object can also be requested by name.
 *I vertex1Name, vertex2Name
 *I- The names (labels) of the start and end vertex'es of the edge.
 *O edgeHdl
 *O- Handle to the newly constructed Edge object.
 *- This function creates a new, named edge object given the end Vertex'es and
 *- other required data.
*/

extern int cgm_quadratic_edge_construct( enum GeometryType geometryType,
                                        char const* name,
                                        char const* vertex1Name,
                                        char const* vertex2Name,
                                        int dimension,
                                        double* coordinates,
                                        enum CubitSense edgeDirection,
                                        EdgeHandle** edgeHdl );
/*
 *R int
 *R- Return a status: 0 => Success
 *R-                 <0 => Failure
 *I GeometryType
 *I- enum describing the type of geometry being constructed.
 *I- Allowed geometry type enums for this function:
 *I-   ELLIPSE_CURVE_TYPE
 *I-   PARABOLA_CURVE_TYPE
 *I name
 *I- The name or label of the Edge object. The object can be
 *I- referenced by name after it has been created. A handle to
 *I- this object can also be requested by name.
 *I vertex1Name, vertex2Name
 *I- The names (labels) of the start and end vertex'es of the edge.
 *I dimension
 *I- An integer indicating the dimensionality of the geometric space
 *I- (2 for 2D and 3 for 3D).
 *I coordinates
 *I- The array of coordinates for the center of the Ellipse.
 *I edgeDirection
 *I- This is applicable only if the ellipse is in 2D space. In this case,
 *I- if edgeDirection is FORWARD, the edge goes from Vertex1 to Vertex2
 *I- such that the right hand rule gives the direction of the Z axis. If
 *I- edgeDirection is REVERSED, the edge goes from Vertex1 to Vertex2
 *I- such that the right hand rule gives the negative of the Z axis.
 *I- This argument is ignored if it is not an ELLIPSE_CURVE_TYPE.
 *O edgeHdl
 *O- Handle to the newly constructed Edge object.
 *- This function creates a new, named edge object given the end Vertex'es and
 *- other required data.
 *-
 *- ELLIPSE_CURVE_TYPE
 *- If the coordinates are the center of a circular arc that
 *- passes through the 2 end Vertex'es, then a circular arc is
 *- generated.
 *- If the coordinates are not the center of either a circle
 *- or an ellipse passing through the 2 end Vertex'es, then a
 *- logarithmic spiral is generated.
 *-
 *- PARABOLA_CURVE_TYPE
 *- In this case, the coordinates must be those of the "top"
 *- or "bottom" of the parabola. The 3 points must form an isosceles
 *- triangle. This definition limits the user to the generation of
 *- the (symmetric) tip of parabolic shapes only.
 *-
*/

extern int cgm_composite_edge_construct( char const* name,
                                        int numEdges,
                                        char const** edgeNameList,
                                        EdgeHandle** edgeHdl );
/*
 *R int
 *R- Return a status: 0 => Success
 *R-                 <0 => Failure
 *I name
 *I- The name or label of the CompositeEdge object. The object can be
 *I- referenced by name after it has been created. A handle to
 *I- this object can also be requested by name.
 *I numEdges
 *I- Number of Edges in edgeNameList
 *I edgeNameList
 *I- An array of names (labels) of the edges that make up the CompositeEdge
 *I- to be constructed. It is important that these be ordered from the
 *I- first Edge to the last Edge.
 *O edgeHdl
 *O- Handle to the newly constructed CompositeEdge object.
 *- This function takes an ordered list of previously constructed Edges and
 *- generates a new, named CompositeEdge.
 *- Note that any subset of the input Edges can be CompositeEdges themselves.
 */

extern int cgm_face_polygon_construct_by_coords( enum GeometryType refFaceType,
						 char const* name,
						 int num_coords,
						 double* double_list,
						 FaceHandle** faceHandle);
  /*
 *R int
 *R- Return a status: 0 => Success
 *R-                 <0 => Failure
 *I name
 *I- The name or label of the Face object. The object can be
 *I- referenced by name after it has been created. A handle to
 *I- this object can also be requested by name.
 *I num_coords
 *I-  int, how many coords to expect
 *I double_list
 *I-  doubles, representing the corners of face
 *O faceHandle
 *O- Handle to the newly constructed Face object.
 *- This function takes an ordered list of doubles and
 *- generates a new, named Face.
*/
extern int cgm_face_polygon_construct_by_verts( enum GeometryType refFaceType,
						char const* name,
						int num_verts,
						char const** vertNameList,
						FaceHandle** faceHandle);
/*
 *R int
 *R- Return a status: 0 => Success
 *R-                 <0 => Failure
 *I name
 *I- The name or label of the Face object. The object can be
 *I- referenced by name after it has been created. A handle to
 *I- this object can also be requested by name.
 *I num_verts
 *I-  int, how many verts to expect
 *I vertNamelist
 *I- list of names of verts, representing the corners of rectangular face
 *O faceHandle
 *O- Handle to the newly constructed RectangleFace object.
 *- This function takes an ordered list of doubles and
 *- generates a new, named RectangleFace.
*/
   
extern int cgm_face_disk_construct( enum GeometryType refFaceType,
				    char const* name,
				    double* coords,
				    double radius,
				    FaceHandle** faceHandle);
/*
 *R int
 *R- Return a status: 0 => Success
 *R-                 <0 => Failure
 *I name
 *I- The name or label of the Face object. The object can be
 *I- referenced by name after it has been created. A handle to
 *I- this object can also be requested by name.
 *I coords
 *I-  double, center of circle
 *I radius
 *I- Okay, if you don't know what the radius is, there's not a whole
 *I- lot I can do to make it clearer here 
 *O faceHandle
 *O- Handle to the newly constructed RectangleFace object.
 *- This function takes an ordered list of doubles and
 *- generates a new, named RectangleFace.
*/
  
extern int cgm_face_construct_by_edges( enum GeometryType geometryType,
					char const* name,
					int numEdges,
					char const** edgeNameList,
					FaceHandle** faceHdl );
  /*
 *R int
 *R- Return a status: 0 => Success
 *R-                 <0 => Failure
  *I GeometryType
 *I- enum describing the type of geometry being constructed.
 *I- Allowed geometry type enums for this function:
 *I-   PLANE_SURFACE_TYPE
 *I name
 *I- The name or label of the Face object. The object can be
 *I- referenced by name after it has been created. A handle to
 *I- this object can also be requested by name.
 *I numEdges
 *I- Number of Edges in edgeNameList
 *I edgeNameList
 *I- An array of names (labels) of the edges that make up the Face
 *I- to be constructed. It is important that these be ordered from the
 *I- first Edge to the last Edge around the perimeter of the Face.
 *O faceHdl
 *O- Handle to the newly constructed Face object.
 *- This function takes an ordered list of previously constructed Edges and
 *- generates a new, named Face.
 *- Note that at this time, none of the input Edges can be CompositeEdges.
 *- Note that the new Face will have a single (external) Loop.
*/

extern int cgm_single_face_body_construct( char const* name,
					   char const* singleFaceName,
					   BodyHandle** bodyHdl );
/*
 *R int
 *R- Return a status: 0 => Success
 *R-                 <0 => Failure
 *I name
 *I- The name or label of the Body object. The object can be
 *I- referenced by name after it has been created. A handle to
 *I- this object can also be requested by name.
 *I singleFaceName
 *I- The name or label of the single Face object that the Body is constructed
 *I- from. This function should be able to reference this Face object
 *I- by name or it will fail to construct the Body. A handle to the Body
 *I- object can be requested by name after it is constructed.
 *O bodyHdl
 *O- Handle to the newly constructed Body object.
 *- This function takes a single previously-constructed Face and
 *- generates a new, named Body.
*/

extern int cgm_multi_face_body_construct( char const* name,
					  int numFaces,
					  char const** faceNameArray,
					  BodyHandle** bodyHdl );
/*
 *R int
 *R- Return a status: 0 => Success
 *R-                 <0 => Failure
 *I name
 *I- The name or label of the Body object. The object can be
 *I- referenced by name after it has been created. A handle to
 *I- this object can also be requested by name.
 *I numFaces
 *I- Number of Faces in the input array, face_name_array.
 *I faceNameArray
 *I- The array of Face object names that the Body is constructed
 *I- from. Each Face becomes an independent Lump. The Body is
 *I- a simple composition of all the Lumps. The Faces must be
 *I- spatially disjoint and must not share any lower level
 *I- topological entities.
 *I- This function should be able to reference the Face objects
 *I- by name or it will fail to construct the Body. A handle to the Body
 *I- object can be requested by name after it is constructed.
 *O bodyHdl
 *O- Handle to the newly constructed Body object.
 *- This function takes a single previously-constructed Face and
 *- generates a new, named Body.
*/

/*************END Construction and Destruction Functions*****************/


/**********************BEGIN Query Functions*****************************/

/**************BEGIN General Query Functions**************/

/* extern int cgm_get_vertex_by_name( char const* name, */
/*                                    VertexHandle** vertexHdl ); */
/*
 *R int
 *R- Return a status: 0 => Success
 *R-                 <0 => Failure
 *I name
 *I- The name or label of the requested Vertex object.
 *O vertexHdl
 *O- Handle to the Vertex object.
 *- This function returns a handle to the Vertex object with the
 *- input name. The handle is NULL if such a Vertex has not been
 *- constructed.
 *- Note: Memory for the returned Handle(s) is allocated by this routine and
 *-       must be freed by the calling code.
*/

extern int cgm_get_edge_by_name( char const* name,
                                 EdgeHandle** edgeHdl );

/*
 *R int
 *R- Return a status: 0 => Success
 *R-                 <0 => Failure
 *I name
 *I- The name or label of the requested Edge object
 *O edgeHdl
 *O- Handle to the Edge object (primitive or composite)
 *- This function returns a handle to the Edge object with the
 *- input name. The handle is NULL if such an Edge has not been
 *- constructed.
 *- Note: Memory for the returned Handle(s) is allocated by this routine and
 *-       must be freed by the calling code.
*/

extern int cgm_get_face_by_name(char const* name,
				FaceHandle** faceHdl);
/*
 *R int
 *R- Return a status: 0 => Success
 *R-                 <0 => Failure
 *I name
 *I- The name or label of the requested Face object.
 *O faceHdl
 *O- Handle to the Face object.
 *- This function returns a handle to the Face object with the
 *- input name. The handle is NULL if such an Face has not been
 *- constructed.
 *- Note: Memory for the returned Handle(s) is allocated by this routine and
 *-       must be freed by the calling code.
*/

extern int cgm_get_volume_by_name( char const* name,
                                   VolumeHandle** volumeHdl );
/*
 *R int
 *R- Return a status: 0 => Success
 *R-                 <0 => Failure
 *I name
 *I- The name or label of the requested Volume object.
 *O volumeHdl
 *O- Handle to the Volume object.
 *- This function returns a handle to the Volume object with the
 *- input name. The handle is NULL if such an Volume has not been
 *- constructed.
 *- Note: Memory for the returned Handle(s) is allocated by this routine and
 *-       must be freed by the calling code.
*/

extern int cgm_get_body_by_name( char const* name,
                                 BodyHandle** bodyHdl );
/*
 *R int
 *R- Return a status: 0 => Success
 *R-                 <0 => Failure
 *I name
 *I- The name or label of the requested Body object.
 *O bodyHdl
 *O- Handle to the Body object.
 *- This function returns a handle to the Body object with the
 *- input name. The handle is NULL if such an Body has not been
 *- constructed.
 *- Note: Memory for the returned Handle(s) is allocated by this routine and
 *-       must be freed by the calling code.
*/

/**************END General Query Functions**************/


/**************BEGIN Vertex Query Functions**************/

/* extern char const* cgm_vertex_get_name( VertexHandle const* vertexHdl ); */
/*
 *R char*
 *R- Pointer to a string.
 *I vertexHdl
 *I- Handle to a Vertex
 *-Function to get the name of the input geometric entity. If the entity does
 *-not have a name, it returns a NULL pointer.
 *-Note: Memory for the name is allocated within the CGM and should be freed
 *-      by the calling code.
*/

extern int cgm_vertex_get_coordinates( VertexHandle const* vertexHdl,
                                       int* num_coords,
                                       double** coordinates );
/*
 *R int
 *R- Return a status: 0 => Success
 *R-                 <0 => Failure
 *I vertexHdl
 *I- Handle to a Vertex
 *O num_coords
 *O- Number of coordinates of the Vertex
 *O coordinates
 *O- Array of num_coords doubles
 *-Function to get the coordinates of a Vertex. It returns failure
 *-if memory for the coordinates cannot be allocated.
 *-Note: Memory for the coordinates are allocated within the CGM.
*/

/**************END Vertex Query Functions**************/


/**************BEGIN Edge Query Functions**************/

/* extern char const* cgm_edge_get_name( EdgeHandle const* edgeHdl ); */
/*
 *R char*
 *R- Pointer to a string.
 *I edgeHdl
 *I- Handle to an Edge
 *-Function to get the name of the input geometric entity. If the entity does
 *-not have a name, it returns a NULL pointer.
 *-Note: Memory for the name is allocated within the CGM and should be freed
 *-      by the calling code.
*/

extern int cgm_edge_get_closest_point ( EdgeHandle const* edgeHdl,
					double point_X,
					double point_Y,
					double point_Z,
					double* closest_point_X,
					double* closest_point_Y,
					double* closest_point_Z,
					double* tangent_vector,
					double* distance );
/*
 *R int
 *R- Return a status: 0 => Success
 *R-                 <0 => Failure
 *I edgeHdl
 *I- Handle to an Edge (primitive or composite)
 *I point_X, point_Y, point_Z
 *I- Input location vector
 *O closest_point_X, closest_point_Y, closest_point_Z
 *O- Point on the Edge that is closest to the input point
 *O tangent_vector
 *O- The tangent vector to the Edge at the closest point. Note that
 *O- this is in the positive direction of the Edge (going from its
 *O- start Vertex to its end Vertex.
 *O distance
 *O- Shortest distance from the input point to the non-trimmed Edge
 *- Function to get the shortest distance from the input point to the Edge,
 *- the closest point on the Edge to the input point, and the tangent to
 *- the Edge at that point.
 *- Note that this function deals with the untrimmed Edge. Think of a
 *- topological Edge spilling out beyond its end Vertices - the geometry
 *- of the topological Edge is trimmed by the Edge. Hence, the closest
 *- trimmed point is not always on the topologically bounded Edge.
*/

extern int cgm_edge_get_closest_point_trimmed ( EdgeHandle const* edgeHdl,
						double point_X,
						double point_Y,
						double point_Z,
						double* closest_point_X,
						double* closest_point_Y,
						double* closest_point_Z,
						double* tangent_vector,
						double* distance );
/*
 *R int
 *R- Return a status: 0 => Success
 *R-                 <0 => Failure
 *I edgeHdl
 *I- Handle to an Edge (primitive or composite)
 *I point_X, point_Y, point_Z
 *I- Input location vector
 *O closest_point_X, closest_point_Y, closest_point_Z
 *O- Point on the "trimmed" (topology) Edge that is closest to the input point
 *O tangent_vector
 *O- The tangent vector to the Edge at the closest point. Note that
 *O- this is in the positive direction of the Edge (going from its
 *O- start Vertex to its end Vertex.
 *O distance
 *O- Shortest distance from the input point to the trimmed Edge
 *- Function to get the shortest distance from the input point to the Edge,
 *- the closest point on the Edge to the input point, and the tangent to
 *- the Edge at that point.
 *- Note that this function deals with the trimmed or topological Edge
 *- that is bounded by its end Vertices. Hence, the closest trimmed point
 *- is always on the topologically bounded Edge.
*/

extern int cgm_edge_evaluate_trimmed ( EdgeHandle const* edgeHdl,
                                       double inputCoordinate1,
                                       unsigned int coord1Index,
                                       double inputCoordinate2,
                                       unsigned int coord2Index,
                                       int* numCoord3Found,
                                       double** coord3Array );
/*
 *R int
 *R- Return a status: 0 => Success
 *R-                 <0 => Failure
 *R- Failure status is returned if no points are found on the Edge having
 *R- the input coordinates.
 *I edgeHdl
 *I- Handle to an Edge (primitive or composite)
 *I inputCoordinate1, inputCoordinate2
 *I- Input coordinates (2 of the 3 coordinates for points on a general Edge -
 *I- an Edge that is embedded in 3-Space)
 *I coord1Index, coord2Index
 *I- Indices of the input coordinates (0, 1 or 2 - implying e.g., X, Y or Z)
 *O numCoord3Found
 *O- Number of points found along the Edge (any positive number including 0)
 *O coord3Array
 *O- Array of 3rd coordinates of the points found along the Edge
 *O- NOTE: This function allocates space for the array. The caller is responsible
 *O-       for releasing this memory.
 *- Function to find the points along the (trimmed) Edge all of which have the
 *- input 2 coordinates. The function returns the list of the 3rd coordinate
 *- of each of the points.
 *- Note that the function could return 0 points. In this case, a Failure status
 *- is also returned. Success is returned for all other cases.
 *- Note that this function deals with the trimmed or topological Edge
 *- that is bounded by its end Vertices.
*/

extern int cgm_edge_get_vertices( EdgeHandle const* edgeHdl,
                                  VertexHandle** startVertexHandle,
                                  VertexHandle** endVertexHandle );
/*
 *R int
 *R- Return a status: 0 => Success
 *R-                 <0 => Failure
 *I edgeHdl
 *I- Handle to an Edge
 *O startVertexHandle, endVertexHandle
 *O- Start and End VertexHandle objects for the Edge
 *-Function to get the start and end Vertex objects of the input Edge.
 *-This function works for both Primitive and Composite Edges.
 *-Errors occur only if there is a problem with the input EdgeHandle,
 *-or if there was a problem retrieving the start and end Vertices of the Edge.
 *- Note: Memory for the returned Handle(s) is allocated by this routine and
 *-       must be freed by the calling code.
*/

/**************END Edge Query Functions**************/


/**************BEGIN Face Query Functions**************/

/* extern char const* cgm_face_get_name( FaceHandle const* faceHdl ); */
/*
 *R char*
 *R- Pointer to a string.
 *I faceHdl
 *I- Handle to a Face
 *-Function to get the name of the input geometric entity. If the entity does
 *-not have a name, it returns a NULL pointer.
 *-Note: Memory for the name is allocated within the CGM and should be freed
 *-      by the calling code.
*/

/* extern void cgm_face_get_closest_point_trimmed ( */
/*                                       FaceHandle const* faceHdl, */
/*                                       double point_X, */
/*                                       double point_Y, */
/*                                       double point_Z, */
/*                                       double* closest_point_X, */
/*                                       double* closest_point_Y, */
/*                                       double* closest_point_Z ); */
/*
 *I faceHdl
 *I- Handle to a Face
 *I point_X, point_Y, point_Z
 *I- Input location vector
 *O closest_point_X, closest_point_Y, closest_point_Z
 *O- Location Vector of the closest point on the Face
 *- Function to get the Location Vector of the closest point on the trimmed
 *- (topological) input Face.
*/

   extern int cgm_named_face_get_closest_boundary( char const* name, 
                                                     double point_x,
                                                     double point_y,
                                                     double point_z,
                                                     double* distance);
   
   extern int cgm_face_get_closest_boundary( FaceHandle const* faceHdl,
                                             double point_x,
                                             double point_y,
                                             double point_z,
                                             double* distance);

  extern int cgm_named_face_get_closest_boundary_keep(char const * name,
						      double point[3],
						      double * distance,
						      double closest[3]);
/*
 *R int
 *R- -1 = failure, 1 = success

 *I faceHdl
 *I- Handle to a Face
 *I point_X, point_Y, point_Z
 *I- Input location vector
 *O double
 *O distance to boundary of face.
 *O Note sign convention-- Positive means interior, negative means exterior
 *- Function to get the signed distance of a point from the boundary of
 *- a RefFace
 */
   extern int cgm_named_volume_get_closest_boundary( char const* name, 
                                                     double point_x,
                                                     double point_y,
                                                     double point_z,
                                                     double* distance);
   
   extern int cgm_volume_get_closest_boundary( VolumeHandle const* volHdl,
                                             double point_x,
                                             double point_y,
                                             double point_z,
                                             double* distance);

  extern int cgm_named_volume_get_closest_boundary_keep(char const * name,
							double point[3],
							double * distance,
							double closest[3]);
/*
 *R int
 *R- -1 = failure, 1 = success

 *I volumeHdl
 *I- Handle to a Volume
 *I point_X, point_Y, point_Z
 *I- Input location vector
 *O double
 *O distance to boundary of volume.
 *O Note sign convention-- Positive means interior, negative means exterior
 *- Function to get the signed distance of a point from the boundary of
 *- a RefVolume
 */
   
/* extern void cgm_face_get_normal ( FaceHandle const* faceHdl, */
/*                                   double point_X, */
/*                                   double point_Y, */
/*                                   double point_Z, */
/*                                   double* closest_point_X, */
/*                                   double* closest_point_Y, */
/*                                   double* closest_point_Z, */
/*                                   double* unit_normal_vector[] ); */
/*
 *I faceHdl
 *I- Handle to a Face
 *I point_X, point_Y, point_Z
 *I- Input location vector
 *O closest_point_X, closest_point_Y, closest_point_Z
 *O- Location Vector of the closest point on the Face to the input point.
 *O unit_normal_vector
 *O- Unit normal vector to the Face at the closest point.
 *- Function to get the Location Vector of the closest point on the trimmed
 *- (topological) input Face and then compute the unit normal vector to the Face
 *- at that point.
*/

/* extern void cgm_face_get_principal_curvatures ( */
/*                                   FaceHandle const* faceHdl, */
/*                                   double point_X, */
/*                                   double point_Y, */
/*                                   double point_Z, */
/*                                   double* closest_point_X, */
/*                                   double* closest_point_Y, */
/*                                   double* closest_point_Z, */
/*                                   double* curvature1, */
/*                                   double* curvature2 ); */
/*
 *I faceHdl
 *I- Handle to a Face
 *I point_X, point_Y, point_Z
 *I- Input location vector
 *O closest_point_X, closest_point_Y, closest_point_Z
 *O- Location Vector of the closest point on the Face to the input point.
 *O curvature1, curvature2
 *O- Values of the principal curvatures of the Face at the closest point.
 *- Function to get the Location Vector of the closest point on the trimmed
 *- (topological) input Face and then compute the principal curvatures
 *- of the Face at the closest point.
*/

/**************END Face Query Functions**************/


/**************BEGIN Plane Query Functions**************/

/* extern void cgm_plane_get_coefficients( PlaneHandle const* planeHdl, */
/*                                       double* a, */
/*                                       double* b, */
/*                                       double* c, */
/*                                       double* d ); */

/*
 *I planeHdl
 *I- Handle to a Plane
 *O a, b, c, d
 *O- Plane cooefficients
 *- Function to get the coefficients of the plane equation of the input Plane:
 *-     aX + bY + cZ + D = 0
*/

/* extern void cgm_plane_get_shortest_distance_to_point ( */
/*                                       PlaneHandle const* planeHdl, */
/*                                       double point_X, */
/*                                       double point_Y, */
/*                                       double point_Z, */
/*                                       double* distance ); */
/*
 *I planeHdl
 *I- Handle to a Plane
 *I point_X, point_Y, point_Z
 *I- Input location vector
 *O distance
 *O- Shortest distance from the input point to the Plane
 *- Function to get the shortest distance from the input point to the Plane.
*/

/**************END Plane Query Functions**************/


/**************BEGIN Body Query Functions**************/

/* extern char const* cgm_body_get_name( BodyHandle const* bodyHdl ); */
/*
 *R char*
 *R- Pointer to a string.
 *I bodyHdl
 *I- Handle to a Body
 *-Function to get the name of the input geometric entity. If the entity does
 *-not have a name, it returns a NULL pointer.
 *-Note: Memory for the name is allocated within the CGM and should be freed
 *-      by the calling code.
*/

/* extern int cgm_body_get_volumes( BodyHandle const* bodyHdl, */
/*                                  int* num_volumes, */
/*                                  VolumeHandle** volumeHandles ); */
/*
 *R int
 *R- Return a status: 0 => Success
 *R-                 <0 => Failure
 *I bodyHdl
 *I- Handle to a Body
 *O num_volumes
 *O- Number of Volumes attached to the input Body
 *O volumeHandles
 *O- Array of num_volumes VolumeHandle objects
 *-Function to get the Volumes attached to the input Body.
 *-Errors occur only if there is a problem with the input BodyHandle,
 *-not if there are no Volumes attached to the Body.
 *- Note: Memory for the returned Handle(s) is allocated by this routine and
 *-       must be freed by the calling code.
*/

/* extern int cgm_body_get_faces( BodyHandle const* bodyHdl, */
/*                                int* num_faces, */
/*                                FaceHandle** faceHandles ); */
/*
 *R int
 *R- Return a status: 0 => Success
 *R-                 <0 => Failure
 *I bodyHdl
 *I- Handle to a Body
 *O num_faces
 *O- Number of Faces attached to the input Body
 *O faceHandles
 *O- Array of num_faces FaceHandle objects
 *-Function to get the Faces attached to the input Body.
 *-Errors occur only if there is a problem with the input BodyHandle,
 *-not if there are no Faces attached to the Body.
 *- Note: Memory for the returned Handle(s) is allocated by this routine and
 *-       must be freed by the calling code.
*/

/* extern int cgm_body_get_edges( BodyHandle const* bodyHdl, */
/*                                int* num_edges, */
/*                                EdgeHandle** edgeHandles ); */
/*
 *R int
 *R- Return a status: 0 => Success
 *R-                 <0 => Failure
 *I bodyHdl
 *I- Handle to a Body
 *O num_edges
 *O- Number of Edges attached to the input Body
 *O edgeHandles
 *O- Array of num_edges EdgeHandle objects
 *-Function to get the Edges attached to the input Body.
 *-Errors occur only if there is a problem with the input BodyHandle,
 *-not if there are no Edges attached to the Body.
 *- Note: Memory for the returned Handle(s) is allocated by this routine and
 *-       must be freed by the calling code.
*/

/* extern int cgm_body_get_vertices( BodyHandle const* bodyHdl, */
/*                                   int* num_vertices, */
/*                                   VertexHandle** vertexHandles ); */
/*
 *R int
 *R- Return a status: 0 => Success
 *R-                 <0 => Failure
 *I bodyHdl
 *I- Handle to a Body
 *O num_vertices
 *O- Number of Vertices attached to the input Body
 *O vertexHandles
 *O- Array of num_vertices VertexHandle objects
 *-Function to get the Vertices attached to the input Body.
 *-Errors occur only if there is a problem with the input BodyHandle,
 *-not if there are no Vertices attached to the Body.
 *- Note: Memory for the returned Handle(s) is allocated by this routine and
 *-       must be freed by the calling code.
*/

/* extern int cgm_body_classify_point ( BodyHandle const* bodyHdl, */
/*                                      double point_X, */
/*                                      double point_Y, */
/*                                      double point_Z ); */
/*
 *R int
 *R- Returned point classification result:
 *R-                        -1 => Inside
 *R-                         0 => On The Boundary
 *R-                        +1 => Outside
 *R-                        <any other value> => Unknown classification
 *I bodyHdl
 *I- Handle to a Body
 *I point_X, point_Y, point_Z
 *I- Input point
 *- Function to classify the input point wrt the input Body. The result is
 *- returned as an integer.
*/

/**************END Body Query Functions**************/

/**************BEGIN Point Query Functions**************/

/* extern double cgm_point_get_closest_distance_to_body_boundary ( */
/*                                      double point_X, */
/*                                      double point_Y, */
/*                                      double point_Z, */
/*                                      BodyHandle const* bodyHdl, */
/*                                      double* closest_point_X, */
/*                                      double* closest_point_Y, */
/*                                      double* closest_point_Z ); */
/*
 *R double
 *R- Closest distance from the input point to the boundary of the input Body:
 *R-                                       < 0.0 => Point is Inside
 *R-                                       = 0.0 => Point is On The Boundary
 *R-                                       > 0.0 => Point is Outside
 *I point_X, point_Y, point_Z
 *I- Input point
 *I bodyHdl
 *I- Handle to a Body
 *I closest_point_X, closest_point_Y, closest_point_Z
 *I- Point on the Body that is closest to the input point
 *- Function to find the closest distance from a point to the boundary of
 *- a Body. The function works for multi-Lump, and multi-Shell/Loop 3D
 *- and 2D Bodies - i.e., the Body can be composed of a disjoint set of
 *- Lumps and can have "holes" in it. The calculation of the closest
 *- distance takes the holes into account. Besides the distance itself,
 *- the closest point on the Body is also returned.
*/

/**************END Point Query Functions**************/

/**********************END Query Functions*****************************/

/**************BEGIN Body Boolean Functions**************/

extern int cgm_volume_subtract( char const* name,
                              VolumeHandle const* body1Hdl,
                              VolumeHandle const* body2Hdl,
                              VolumeHandle** newVolumeHdl ); 
/*
 *R int
 *R- Return a status: 0 => Success
 *R-                 <0 => Failure
 *I name
 *I- Name (label) of the new Volume
 *I body1Hdl
 *I- Handle to a Volume
 *I body2Hdl
 *I- Handle to a Volume
 *O newVolumeHdl
 *O- Handle to the resulting new Volume
 *-Function to subtract body1 from body2. If successful, the result is
 *-the newVolume. The input Bodies are not modified. The new Volume can be
 *- accessed using the input "name".
*/
   
extern int cgm_volume_unite( char const* name,
                              VolumeHandle const* body1Hdl,
                              VolumeHandle const* body2Hdl,
                              VolumeHandle** newVolumeHdl ); 
/*
 *R int
 *R- Return a status: 0 => Success
 *R-                 <0 => Failure
 *I name
 *I- Name (label) of the new Volume
 *I body1Hdl
 *I- Handle to a Volume
 *I body2Hdl
 *I- Handle to a Volume
 *O newVolumeHdl
 *O- Handle to the resulting new Volume
 *-Function to unite body1 and body2. If successful, the result is
 *-the newVolume. The input Bodies are not modified. The new Volume can be
 *- accessed using the input "name".
*/

/**************END Body Boolean Functions**************/

/************ Begin helper Functions************/
extern int get_debug_export();

extern void set_debug_export(int new_state);

extern int is_name_used(char const* name, int print_warning);
   
extern int add_me_to_export(const char * name_of_entity);
/*
 *R int
 *R- Return size of name list
 *I name_of_entity
 *I- Name of entity to add to list--can be an empty "\0" string
 *-Function to append the name of an entity to the current list of entities
 *-that need to be exported
 */

extern int all_done_export_now(const char * name_of_export_file);
/*
 *R int
 *R- Return value is 0 if success, -1 if failure
 *I name_of_export_file
 *I- Name of file for exporting, can't be empty or returns error.
 *-Function to actually export the list of entities previously added to list
 */
extern void reset_export_list();

     /*
      *Does just what its name says
      */
       
      
#ifndef __cplusplus
  extern
#else
  static
#endif
  int _debugExport;
   
#ifdef __cplusplus
}
#endif


/******************** END EXTERN FUNCTIONS      ****************************/

#endif

#endif
