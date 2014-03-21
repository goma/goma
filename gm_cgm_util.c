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
 
#define _GM_CGM_UTIL_C

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "std.h"
#include "rf_mp.h"
#include "rf_io.h"
#include "goma.h"

#ifdef USE_CGM

char cgm_exported_filename[MAX_FNL];
char cgm_input_string[MAX_CGM_INPUT_STRING_LENGTH];
int cgm_input_string_length;

/* ========================================================================= */
/* ========================================================================= */
/*
 * Helper function for printing edge info -- for the purpose of debugging the
 * latest additions to goma.
 */
static void print_edge_info_helper(char const* edgeName)
{
  EdgeHandle* edgeHandle = NULL;
  double points[][3] =
  {
     {0.0, 10.0, 0.0},
     {-10.0, 10.0, 0.0},
     {-10.0, -20.0, 0.0},
     {20.0, 0.0, 0.0},
     {20.0, 20.0, 0.0},
     {20.0, -20.0, 0.0},
  };
  int num_points = sizeof(points)/(sizeof(double)*3);
  int ipoint = 0;
  double closest_point[3];
  double closest_point_trimmed[3];
  double tangent[3];
  double distance;
  int status = 0;

  /*
   * Get a handle to the edge with the above name
   */
/*   printf("Getting edge named %s\n",edgeName); */
  
  status = cgm_get_edge_by_name(edgeName, &edgeHandle);
  if ( status < 0 )
    {
      fprintf(stdout, "There is no edge named '%s'\n", edgeName);
      return;
    }

  /*
   * Display the closest points and closest points trimmed to
   * the edge given the input points.
   */
  for ( ipoint = 0; ipoint < num_points; ipoint++ )
  {
    /*
     * Find the closest point
     */ 
    status = cgm_edge_get_closest_point(edgeHandle,
                                        points[ipoint][0],
                                        points[ipoint][1],
                                        points[ipoint][2],
                                        closest_point,
                                        closest_point+1,
                                        closest_point+2,
                                        tangent,
                                        &distance);
    if ( status < 0 )
      {
        fprintf(stdout, "Unable to get the closest point to an Edge\n");
        continue;
      }

    /*
     * Find the closest trimmed point
     */
    status = cgm_edge_get_closest_point_trimmed(edgeHandle,
                                        points[ipoint][0],
                                        points[ipoint][1],
                                        points[ipoint][2],
                                        closest_point_trimmed,
                                        closest_point_trimmed+1,
                                        closest_point_trimmed+2,
                                        tangent,
                                        &distance);
    if ( status < 0 )
      {
        fprintf(stdout, "Unable to get the closest point to an Edge\n");
        continue;
      }

    /*
     * Print info
     */
    fprintf(stdout, "\n");
    fprintf(stdout, "Input point: %f %f %f\n",
                    points[ipoint][0],
                    points[ipoint][1],
	            points[ipoint][2]);
    fprintf(stdout, "Closest point: %f %f %f\n",
                    closest_point[0],
                    closest_point[1],
                    closest_point[2]);
    fprintf(stdout, "Closest point trimmed: %f %f %f\n",
                    closest_point_trimmed[0],
                    closest_point_trimmed[1],
                    closest_point_trimmed[2]);
  }
/*  fprintf(stdout, "gm_cgm_util::print_edge_info not implemented yet\n"); */
 
}

/* ========================================================================= */
/* ========================================================================= */
/*
 * Print some edge info -- for the purpose of debugging the
 * latest additions to goma.
 */
void print_edge_info()
{
  int status = 0;

  /*
   * Print info about the edge named "straight"
   */
  fprintf(stdout, "\n");
  print_edge_info_helper("straight");

  /*
   * Print info about the edge named "ell"
   */
  fprintf(stdout, "\n");
  print_edge_info_helper("ell");

  /*
   * Print info about the edge named "channel"
   */
  fprintf(stdout, "\n");
  print_edge_info_helper("channel");

  /*
   * Print info about the edge named "curved"
   */
  fprintf(stdout, "\n");
  print_edge_info_helper("curved");

}

/* ========================================================================= */
/* ========================================================================= */
void cgm_initialize()
{
  static const char yo[] = "cgm_initialize";
  int error;

  /*
   * Testing initialization of the CGM subsystem
   */
  error = cpp_cgm_initialize();
  if ( error < 0 )
    {
      /* error initializing the CGM subsystem */
      EH( error, "Problem initializing the CGM subsystem.");
#ifdef PARALLEL
      MPI_Finalize();
#endif
    }
  else   /* successfully initialized the CGM subsystem */
    {
      log_msg("Successfully initialized the CGM subsystem.");
#ifdef PARALLEL
      DPRINTF(stderr, "Successfully initialized the CGM subsystem on P0.\n");
      DPRINTF(stderr, "  Need to uncomment code to see if other processors were successful, too.\n"); 
      /*
	fprintf(stderr, "P%d successfully initialized the CGM subsystem.\n", ProcID);
      */
#else
      fprintf(stderr, "Successfully initialized the CGM subsystem.\n");
#endif
    }
}

/* ========================================================================= */
/* ========================================================================= */
int read_ACIS_file(char *fname,
		   BodyHandle *bodyHdl)
{
  static const char yo[] = "read_ACIS_file";
  int error;

  if(strlen(fname))
    {
      log_msg("Reading geometry from input ACIS SAT file...");

      /* I find it odd that a **BodyHandle is what we pass into the
       * cgm_geometery_read_SAT routine...
       */
      error = cgm_geometry_read_SAT(fname, &bodyHdl ); 
      if ( error < 0 )
	{
#ifdef DEBUG
	  fprintf(stderr, "P_%d at barrier before setup_pd()\n", ProcID);
	  error = MPI_Barrier(MPI_COMM_WORLD);
#endif

	  /* error reading geometry file */
	  /* missing files on any processor are detected at a lower level
	     forcing a return to the higher level
	     read_geometry_SAT -->  main
	     shutdown now, as the geometry file wasn't found on a processor */
	  
	  EH( error, "Problem reading ACIS SAT geometry file.");
#ifdef PARALLEL
	  MPI_Finalize();
#endif
	  return(-1);
	}
      else   /* successfully read the geometry file */
	{
	  log_msg("Successfully read input ACIS SAT file.");
	  /*
#ifdef PARALLEL
	  fprintf(stdout, "P_%d: Successfully read input ACIS SAT file %s.\n", ProcID, fname);
#else
	  fprintf(stdout, "Successfully read input ACIS SAT file %s.\n", fname);
#endif
	  */
	  return(1);
	}
    }
  
    /*
     * Print out some information about the edges we constructed
     */
/*
  print_edge_info();
*/
    /* Testing the point classification code in the CGM */
  
    /* Get handles to the 2 bodies */
  return(1);
}

/* ========================================================================= */
/* ========================================================================= */

/* This routine goes through the big string that all processors own
 * and creates CGM geometry from it. */  
void create_cgm_geometry()
{
  static const char yo[] = "create_cgm_geometry()";
  enum CubitSense edge_sense;
  char err_msg[255];
  char **composited_item_names;
  char *stok, *s1, *s2, *geom_type;
  int i, j, k, num_vertices, num_edges, num_faces, num_bodies,
    num_composited_things;
  dbl coordinate[3];
  dbl *coords;
  
  VertexHandle *vertex_handle;
  EdgeHandle *edge_handle;
  FaceHandle *face_handle;
  VolumeHandle *volume_handle, *volume_handle1, *volume_handle2;
  
  BodyHandle *body_handle, *body_handle1, *body_handle2;
  
  stok = strtok(cgm_input_string, " ");	/* Either "ACIS" or number of vertices. */
  if(!stok)			/* nothing was there except whitespace. */
    return;
  /* See if we begin with an "ACIS file <fname>" specification. */
  if(!strcmp(stok, "ACIS"))
    {
      read_ACIS_file(strtok(NULL, " "), body_handle); /* .sat filename */
      stok = strtok(NULL, " ");	/* number of vertices */
    }

  num_vertices = atoi(stok);
  for(i = 0; i < num_vertices; i++)
    {
      stok = strtok(NULL, " ");	/* vertex name */
      for(j = 0; j < 3; j++)
	coordinate[j] = atof(strtok(NULL, " ")); /* coordinate */
      if(cgm_vertex_construct(UNDEFINED_POINT_TYPE,
			      stok,
			      3,
			      coordinate,
			      &vertex_handle) < 0)
	{
	  sprintf(err_msg, "Could not create VERTEX named \"%s\"\n", stok);
	  EH(-1, err_msg);
	}
			   
    }

  num_edges = atoi(strtok(NULL, " ")); /* number of edges */
  for(i = 0; i < num_edges; i++)
    {
      stok = strtok(NULL, " ");	/* edge name */
      geom_type = strtok(NULL, " "); /* edge type */
      if(strcmp(geom_type, "COMPOSITE"))
	{
	  s1 = strtok(NULL, " "); /* name of vertex 1 */
	  s2 = strtok(NULL, " "); /* name of vertex 2 */
	}
      if(!strcmp(geom_type, "STRAIGHT"))
	{
	  if(cgm_straight_edge_construct(STRAIGHT_CURVE_TYPE,
					 stok, s1, s2, &edge_handle) < 0) 
	    {
	      sprintf(err_msg, "Could not create STRAIGHT EDGE named \"%s\"\n", stok);
	      EH(-1, err_msg);
	    }
	}
      else if(!strcmp(geom_type, "ELLIPSE"))
	{
	  for(j = 0; j < 3; j++)
	    coordinate[j] = atof(strtok(NULL, " "));
	  edge_sense = strcmp(strtok(NULL, " "), "FORWARD") ?
	    CUBIT_REVERSED : CUBIT_FORWARD;
	  if(cgm_quadratic_edge_construct(ELLIPSE_CURVE_TYPE,
					  stok, s1, s2, 3,
					  coordinate, edge_sense,
					  &edge_handle) < 0)
	    {
	      sprintf(err_msg, "Could not create ELLIPSE EDGE named \"%s\"\n", stok);
	      EH(-1, err_msg);
	    }
	}
      else if(!strcmp(geom_type, "PARABOLA"))
	{
	  for(j = 0; j < 3; j++)
	    coordinate[j] = (double)atof(strtok(NULL, " "));
	  edge_sense = strcmp(strtok(NULL, " "), "FORWARD") ?
	    CUBIT_REVERSED : CUBIT_FORWARD;
	  if(cgm_quadratic_edge_construct(PARABOLA_CURVE_TYPE,
					  stok, s1, s2, 3,
					  coordinate, edge_sense,
					  &edge_handle) < 0)
	    {
	      sprintf(err_msg, "Could not create PARABOLA EDGE named \"%s\"\n", stok);
	      EH(-1, err_msg);
	    }
	}
      else if(!strcmp(geom_type, "COMPOSITE"))
	{
	  num_composited_things = atoi(strtok(NULL, " "));
	  composited_item_names = (char **)malloc(num_composited_things * sizeof(char *));
	  for(j = 0; j < num_composited_things; j++)
	    {
	      s1 = strtok(NULL, " ");
	      composited_item_names[j] = (char *)malloc((strlen(s1) + 1) * sizeof(char));
	      strcpy(composited_item_names[j], s1);
	    }
	  if(cgm_composite_edge_construct(stok, num_composited_things,
					  (char const  **)composited_item_names,
					  &edge_handle) < 0)
	    {
	      sprintf(err_msg, "Could not create COMPOSITE EDGE named \"%s\"\n", stok);
	      EH(-1, err_msg);
	    }
	  for(j = 0; j < num_composited_things; j++)
	    free(composited_item_names[j]);
	  free(composited_item_names);
	}
      else
	{
	  sprintf(err_msg, "Huh?  How'd you get here?  Unknown EDGE type: %s\n", geom_type);
	  EH(-1, err_msg);
	}
    }
  
  num_faces = atoi(strtok(NULL, " ")); /* number of faces */
  for(i = 0; i < num_faces; i++)
    {
      stok = strtok(NULL, " ");	/* face name */
      geom_type = strtok(NULL, " "); /* face type */
      if(!strcmp(geom_type, "PLANE"))
	{
	  num_composited_things = atoi(strtok(NULL, " "));
	  composited_item_names = (char **)malloc(num_composited_things * sizeof(char *));
	  for(j = 0; j < num_composited_things; j++)
	    {
	      s1 = strtok(NULL, " ");
	      composited_item_names[j] = (char *)malloc((strlen(s1) + 1) * sizeof(char));
	      strcpy(composited_item_names[j], s1);
	      /*fprintf(stdout, "Added edge '%s' to boundary list.\n",
		composited_item_names[j]);*/
	    }
	  if(cgm_face_construct_by_edges(PLANE_SURFACE_TYPE,
					 stok, num_composited_things,
					 (char const **)composited_item_names,
					 &face_handle) < 0)
	    {
	      sprintf(err_msg, "Could not create PLANE FACE named \"%s\"\n", stok);
	      EH(-1, err_msg);
	    }
	  for(j = 0; j < num_composited_things; j++)
	    free(composited_item_names[j]);
	  free(composited_item_names);
	}
      else if(!strcmp(geom_type, "POLY"))
	{
	  num_composited_things = atoi(strtok(NULL, " "));
	  coords = (double *)malloc(num_composited_things * 3 * sizeof(double));
	  for(j = 0; j < num_composited_things; j++)
	    for(k = 0; k < 3; k++)
	      coords[j*3+k] = atof(strtok(NULL, " "));
	  if(cgm_face_polygon_construct_by_coords(PLANE_SURFACE_TYPE,
						  stok, num_composited_things, coords,
						  &face_handle) < 0)
	    {
	      sprintf(err_msg, "Could not create POLY FACE named \"%s\"\n", stok);
	      EH(-1, err_msg);
	    }
	  free(coords);
	}
      else if(!strcmp(geom_type, "POLY_VERT"))
	{
	  num_composited_things = atoi(strtok(NULL, " "));

	  composited_item_names = (char **)malloc(num_composited_things * sizeof(char *));
	  for(j = 0; j < num_composited_things; j++)
	    {
	      s1 = strtok(NULL, " ");
	      composited_item_names[j] = (char *)malloc((strlen(s1) + 1) * sizeof(char));
	      strcpy(composited_item_names[j], s1);
	    }

	  if(cgm_face_polygon_construct_by_verts(PLANE_SURFACE_TYPE,
						 stok, num_composited_things,
						 (const char **)composited_item_names,
						 &face_handle) < 0)
	    {
	      sprintf(err_msg, "Could not create POLY_VERT FACE named \"%s\"\n", stok);
	      EH(-1, err_msg);
	    }

	  for(j = 0; j < num_composited_things; j++)
	    free(composited_item_names[j]);
	  free(composited_item_names);
	}
      else if(!strcmp(geom_type, "DISK"))
	{
	  for(j = 0; j < 3; j++)
	    {
	      coords[j] = atof(strtok(NULL, " "));
	    }
	  if(cgm_face_disk_construct(PLANE_SURFACE_TYPE,
				     stok, coords, coords[2],
				     &face_handle) < 0)
	    {
	      sprintf(err_msg, "Could not create DISK FACE named \"%s\"\n", stok);
	      EH(-1, err_msg);
	    }
	}
      else
      {
	  sprintf(err_msg, "Huh?  How'd you get here?  Unknown FACE type: %s\n", geom_type);
	  EH(-1, err_msg);
	}
    }

  num_bodies = atoi(strtok(NULL, " "));
  for(i = 0; i < num_bodies; i++)
    {
      stok = strtok(NULL, " "); /* body name */
      geom_type = strtok(NULL, " "); /* body type */
      if(!strcmp(geom_type, "SINGLE_FACE"))
	{
	  s1 = strtok(NULL, " "); /* face name */
	  if(cgm_single_face_body_construct(stok, s1, &body_handle) < 0)
	    {
	      sprintf(err_msg, "Could not create SINGLE_FACE BODY named \"%s\"\n", stok);
	      EH(-1, err_msg);
	    }
	}
      else if(!strcmp(geom_type, "MULTI_FACE"))
	{
	  num_composited_things = atoi(strtok(NULL, " "));
	  composited_item_names = (char **)malloc(num_composited_things * sizeof(char *));
	  for(j = 0; j < num_composited_things; j++)
	    {
	      s1 = strtok(NULL, " ");
	      composited_item_names[j] = (char *)malloc((strlen(s1) + 1) * sizeof(char));
	      strcpy(composited_item_names[j], s1);
	    }
	  if(cgm_multi_face_body_construct(stok, num_composited_things,
					   (char const **)composited_item_names,
					   &body_handle) < 0)
	    {
	      sprintf(err_msg, "Could not create MULTI_FACE BODY named \"%s\"\n", stok);
	      EH(-1, err_msg);
	    }
	  for(j = 0; j < num_composited_things; j++)
	    free(composited_item_names[j]);
	  free(composited_item_names);
	}
      else if(!strcmp(geom_type, "SUBTRACT"))
	{
	  s1 = strtok(NULL, " "); /* name of original body. */
	  s2 = strtok(NULL, " "); /* name of subtracted body. */
	  if(cgm_get_volume_by_name((char const *)s1, &volume_handle1))
	    {
	      sprintf(err_msg, "Could not find BODY named \"%s\"\n", s1);
	      EH(-1, err_msg);
	    }
	  if(cgm_get_volume_by_name((char const *)s2, &volume_handle2))
	    {
	      sprintf(err_msg, "Could not find BODY named \"%s\"\n", s2);
	      EH(-1, err_msg);
	    }
	  if(cgm_volume_subtract(stok, volume_handle1, volume_handle2,
			       &volume_handle) < 0)
	    {
	      sprintf(err_msg, "Could not create SUBTRACT BODY named \"%s\"\n", stok);
	      EH(-1, err_msg);
	    }
	}
      else
	{
	  sprintf(err_msg, "Huh?  How'd you get here?  Unknown BODY type: %s\n", geom_type);
	  EH(-1, err_msg);
	}
    }
  if(cgm_exported_filename[0] != 0)
    all_done_export_now(cgm_exported_filename);
}


#endif /* GM_CGM_UTIL_C */
