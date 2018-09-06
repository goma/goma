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
 
/* mm_input_particles.c -- read in parameters for particles.
 *
 * Author: Matthew M. Hopkins
 * Created: 11/01/2001
 */

/*
 * $Id: mm_input_particles.c,v 5.1 2007-09-18 18:53:45 prschun Exp $
 */

#define MM_INPUT_PARTICLES_C

/* Standard include files */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>

/* GOMA include files */
#include "goma.h"

#define SLEN 256
#define SLEN_2 512

void
rd_particle_specs(FILE *ifp, char *input)
{
  PBC_t *PBC;
  char err_msg[MAX_CHAR_ERR_MSG], s_tmp[SLEN], s_tmp_save[SLEN_2], *s_ptr1, *s_ptr2;
  int i, j, iread;

  iread = look_for_optional(ifp, "Particle Specifications", input, '=');
  if(iread != 1)
    {
      Particle_Dynamics = 0;
      Particle_Number_Sample_Types = 0;
      Particle_Number_Samples_Existing = NULL;
      Particle_Number_Samples = NULL;
      Particle_Number_Output_Variables = NULL;
      Particle_Output_Variables = NULL;
      Particle_Filename_Template = NULL;
      Particle_Number_PBCs = 0;
      PBCs = NULL;
      return;
    }

  printf("\n--- Particle specifications ---\n");
  Particle_Dynamics = 1;	/* turn on particle dynamics */

  /* What kind of particle model do you want?  This determines all
   * kinds of things, such as FEM<->particle coupling, particle
   * trajectory calculation, etc. */
  look_for(ifp, "Particle model", input, '=');
  if(!fgets(s_tmp, SLEN-1, ifp))
    EH(-1, "Error reading Particle model card.");
  strip(s_tmp);			/* s_tmp still has '\n' at end */
  if(s_tmp[strlen(s_tmp)-1] == '\n')
    s_tmp[strlen(s_tmp)-1] = 0;
  strcpy(s_tmp_save, s_tmp);
  s_ptr1 = strpbrk(s_tmp_save, " \t\n");
  if(s_ptr1 != NULL)
    s_ptr1[0] = 0;
  if(!strncmp("INERTIAL_TRACER_EXPLICIT", s_tmp, 24))
    Particle_Model = INERTIAL_TRACER_EXPLICIT;
  else if(!strncmp("INERTIAL_TRACER_IMPLICIT", s_tmp, 24))
    Particle_Model = INERTIAL_TRACER_IMPLICIT;
  else if(!strncmp("TRACER_EXPLICIT", s_tmp, 15))
    Particle_Model = TRACER_EXPLICIT;
  else if(!strncmp("TRACER_IMPLICIT", s_tmp, 15))
    Particle_Model = TRACER_IMPLICIT;
  else if(!strncmp("SWIMMER_EXPLICIT", s_tmp, 16))
    Particle_Model = SWIMMER_EXPLICIT;
  else if(!strncmp("SWIMMER_IMPLICIT", s_tmp, 16))
    Particle_Model = SWIMMER_IMPLICIT;
  else if(!strncmp("CHARGED_TRACER_EXPLICIT", s_tmp, 23))
    Particle_Model = CHARGED_TRACER_EXPLICIT;
  else if(!strncmp("CHARGED_TRACER_IMPLICIT", s_tmp, 23))
    Particle_Model = CHARGED_TRACER_IMPLICIT;
  else if(!strncmp("DIELECTROPHORETIC_TRACER_IMPLICIT", s_tmp, 33))
    Particle_Model = DIELECTROPHORETIC_TRACER_IMPLICIT;
  else
    {
      sprintf(s_tmp_save, "Unknown Particle model: %s\n", s_tmp);
      EH(-1, s_tmp_save);
    }
  printf("Setting Particle_Model = %s\n", s_tmp_save);

  switch(Particle_Model)
    {
    case INERTIAL_TRACER_IMPLICIT:
    case INERTIAL_TRACER_EXPLICIT:
      s_ptr1 = s_tmp + 24;
      if(!(s_ptr2 = strtok(s_ptr1, " \t\n")))
	{
	  sprintf(err_msg, "Error reading random walk coefficient for %s particle model.", s_tmp);
	  EH(-1, err_msg);
	}
      Particle_Model_Data[0] = atof(s_ptr2);
      printf("Setting random walk coefficient = %g\n", Particle_Model_Data[0]);
      
      if(!(s_ptr2 = strtok(NULL, " \t\n")))
	{
	  sprintf(err_msg, "Error reading gravity coefficient for %s particle_model.", s_tmp);
	  EH(-1, err_msg);
	}
      Particle_Model_Data[1] = atof(s_ptr2);
      printf("Setting gravity coefficient = %g\n", Particle_Model_Data[1]);
      break;
    case TRACER_EXPLICIT:
    case TRACER_IMPLICIT:
      s_ptr1 = s_tmp + 15;
      if(!(s_ptr2 = strtok(s_ptr1, " \t\n")))
	{
	  sprintf(err_msg, "Error reading random walk coefficient for %s particle model.", s_tmp);
	  EH(-1, err_msg);
	}
      Particle_Model_Data[0] = atof(s_ptr2);
      printf("Setting random walk coefficient = %g\n", Particle_Model_Data[0]);
      break;
    case SWIMMER_EXPLICIT:
    case SWIMMER_IMPLICIT:
      s_ptr1 = s_tmp + 16;
      if(!(s_ptr2 = strtok(s_ptr1, " \t\n")))
	{
	  sprintf(err_msg, "Error reading random swimming speed coefficient for %s particle model.", s_tmp);
	  EH(-1, err_msg);
	}
      Particle_Model_Data[0] = atof(s_ptr2);
      printf("Setting random swimming speed coefficient = %g\n", Particle_Model_Data[0]);
      if(!(s_ptr2 = strtok(NULL, " \t\n")))
	{
	  sprintf(err_msg, "Error reading random cell orientation coefficient for %s particle model.", s_tmp);
	  EH(-1, err_msg);
	}
      Particle_Model_Data[1] = atof(s_ptr2);
      printf("Setting random cell orientation coefficient = %g\n", Particle_Model_Data[1]);
      if(!(s_ptr2 = strtok(NULL, " \t\n")))
	{
	  sprintf(err_msg, "Error reading upswimming speed for %s particle model.", s_tmp);
	  EH(-1, err_msg);
	}
      Particle_Model_Data[2] = atof(s_ptr2);
      printf("Setting particle upswimming speed = %g\n", Particle_Model_Data[2]);
      break;
    case CHARGED_TRACER_EXPLICIT:
    case CHARGED_TRACER_IMPLICIT:
      s_ptr1 = s_tmp + 23;
      if(!(s_ptr2 = strtok(s_ptr1, " \t\n")))
	{
	  sprintf(err_msg, "Error reading random walk coefficient for %s particle model.", s_tmp);
	  EH(-1, err_msg);
	}
      Particle_Model_Data[0] = atof(s_ptr2);
      printf("Setting random walk coefficient = %g\n", Particle_Model_Data[0]);
      if(!(s_ptr2 = strtok(NULL, " \t\n")))
	{
	  sprintf(err_msg, "Error reading charge coefficient for %s particle model.", s_tmp);
	  EH(-1, err_msg);
	}
      Particle_Model_Data[1] = atof(s_ptr2);
      printf("Setting charge coefficient = %g\n", Particle_Model_Data[1]);
      break;
    case DIELECTROPHORETIC_TRACER_IMPLICIT:
      s_ptr1 = s_tmp + 33;
      if(!(s_ptr2 = strtok(s_ptr1, " \t\n")))
	{
	  sprintf(err_msg, "Error reading random walk coefficient for %s particle model.", s_tmp);
	  EH(-1, err_msg);
	}
      Particle_Model_Data[0] = atof(s_ptr2);
      printf("Setting random walk coefficient = %g\n", Particle_Model_Data[0]);
      if(!(s_ptr2 = strtok(NULL, " \t\n")))
	{
	  sprintf(err_msg, "Error reading particle permittivity for %s particle model.", s_tmp);
	  EH(-1, err_msg);
	}
      Particle_Model_Data[1] = atof(s_ptr2);
      printf("Setting particle permittivity = %g\n", Particle_Model_Data[1]);
      if(!(s_ptr2 = strtok(NULL, " \t\n")))
	{
	  sprintf(err_msg, "Error reading medium (fluid) permittivity for %s particle model.", s_tmp);
	  EH(-1, err_msg);
	}
      Particle_Model_Data[2] = atof(s_ptr2);
      printf("Setting medium (fluid) permittivity = %g\n", Particle_Model_Data[2]);
      if(!(s_ptr2 = strtok(NULL, " \t\n")))
	{
	  sprintf(err_msg, "Error reading particle conductivity for %s particle model.", s_tmp);
	  EH(-1, err_msg);
	}
      Particle_Model_Data[3] = atof(s_ptr2);
      printf("Setting particle conductivity = %g\n", Particle_Model_Data[3]);
      if(!(s_ptr2 = strtok(NULL, " \t\n")))
	{
	  sprintf(err_msg, "Error reading medium (fluid) conductivity for %s particle model.", s_tmp);
	  EH(-1, err_msg);
	}
      Particle_Model_Data[4] = atof(s_ptr2);
      printf("Setting medium (fluid) conductivity = %g\n", Particle_Model_Data[4]);
      if(!(s_ptr2 = strtok(NULL, " \t\n")))
	{
	  sprintf(err_msg, "Error reading AC angular frequency for %s particle model.", s_tmp);
	  EH(-1, err_msg);
	}
      Particle_Model_Data[5] = atof(s_ptr2);
      printf("Setting AC angular frequency = %g\n", Particle_Model_Data[5]);
      if(!(s_ptr2 = strtok(NULL, " \t\n")))
	{
	  sprintf(err_msg, "Error reading Volt conversion factor for %s particle model.", s_tmp);
	  EH(-1, err_msg);
	}
      Particle_Model_Data[6] = atof(s_ptr2);
      printf("Setting Volt conversion factor = %g\n", Particle_Model_Data[6]);
      break;
    default:
      EH(-1, "How'd you get here?");
    }

  /* Get maximum number of particle time steps if we're running after
   * a steady-state solution... */
  if(TimeIntegration == STEADY)
    {
      if(look_for_optional(ifp, "Time steps", input, '=') == 1)
	{
	  if(fscanf(ifp, "%d", &Particle_Max_Time_Steps) != 1)
	    EH(-1, "Error reading Time steps card for post-steady problem.");
	  printf("Setting Particle_Max_Time_Steps = %d\n", Particle_Max_Time_Steps);
	}
      else
	Particle_Max_Time_Steps = 0;
    }
  else
    if(look_for_optional(ifp, "Time steps", input, '=') == 1)
      printf("CAUTION: ignoring Time steps card for transient problem.");

  /* How many particles would you like?  For now, they are always
   * introduced unformly throughout the domain. */
  look_for(ifp, "Number of particles", input, '=');
  if(fscanf(ifp, "%d", &Particle_Number) != 1)
    EH(-1, "Error reading Number of particles card.");
  printf("Setting Particle_Number = %d\n", Particle_Number);

  if(look_for_optional(ifp, "Restart file", input, '=') == 1)
    {
      if(!fgets(s_tmp, MAX_PARTICLE_FILENAME_LENGTH, ifp))
	EH(-1, "Error reading Restart file card.");
      strip(s_tmp);
      strncpy(Particle_Restart_Filename, s_tmp, strlen(s_tmp) - 1);
      printf("Restart file = %s\n", Particle_Restart_Filename);
    }
  else
    sprintf(Particle_Restart_Filename, "<not active>");

  /* When to output particle information.  There are two ways to
   * specify this, and they are incompatible (i.e., only one can be
   * specified).  Either you can output particles every so many Goma
   * time steps, or you can have particles output every so often
   * time-wise (e.g., every 0.1 seconds), whether within a Goma time
   * step or containing many Goma time steps.*/ 
  if(look_for_optional(ifp, "Output stride", input, '=') == 1)
    {
      if(fscanf(ifp, "%d", &Particle_Output_Stride) != 1)
	EH(-1, "Error reading Output stride card.");
      printf("Setting Particle_Output_Stride = %d\n", Particle_Output_Stride);
    }
  else
    Particle_Output_Stride = -1;

  if(look_for_optional(ifp, "Output time step", input, '=') == 1)
    {
      if(fscanf(ifp, "%lf", &Particle_Output_Time_Step) != 1)
	EH(-1, "Error reading Output time step card.");
      printf("Setting Particle_Output_Time_Step = %g\n", Particle_Output_Time_Step);
      if(Particle_Output_Stride != -1)
	EH(-1, "You cannot specify both Output stride and Output time step.");
    }
  else
    Particle_Output_Time_Step = 0.0;

  if(Particle_Output_Stride == -1 &&
     Particle_Output_Time_Step == 0.0)
    EH(-1, "One of the Output stride or Output time step cards must be present.");

  if(look_for_optional(ifp, "Output format", input, '=') == 1)
    {
      if(!fgets(s_tmp, SLEN-1, ifp))
	EH(-1, "Error reading Output format card.");
      strip(s_tmp);			/* s_tmp still has '\n' at end -- don't care. */
      if(!strncmp("TECPLOT", s_tmp, 7))
	Particle_Output_Format = TECPLOT;
      else if(!strncmp("FLAT_TEXT", s_tmp, 9))
	Particle_Output_Format = FLAT_TEXT;
      else
	{
	  sprintf(s_tmp_save, "Unknown Output format: %s\n", s_tmp);
	  EH(-1, s_tmp_save);
	}
    }
  else
    Particle_Output_Format = FLAT_TEXT;

  if(look_for_optional(ifp, "Particle density", input, '=') == 1)
    {
      if(fscanf(ifp, "%lf", &Particle_Density) != 1)
	EH(-1, "Error reading Particle density card.");
    }
  else
    Particle_Density = 1.0;
  printf("Setting Particle_Density = %g\n", Particle_Density);

  if(look_for_optional(ifp, "Particle radius", input, '=') == 1)
    {
      if(fscanf(ifp, "%lf", &Particle_Radius) != 1)
	EH(-1, "Error reading Particle radius card.");
    }
  else
    Particle_Radius = 1.0;
  printf("Setting Particle_Radius = %g\n", Particle_Radius);
  if(Particle_Radius <= 0.0)
    EH(-1, "Particle radius must be > 0.0.  Non-inertial tracer particles will ignore this setting.");

  if(look_for_optional(ifp, "Particle ratio", input, '=') == 1)
    {
      if(fscanf(ifp, "%lf", &Particle_Ratio) != 1)
	EH(-1, "Error reading Particle ratio card.");
    }
  else
    Particle_Ratio = 1.0;
  printf("Setting Particle_Ratio = %g\n", Particle_Ratio);

  Particle_Show_Debug_Info = 0;
  if(look_for_optional(ifp, "Show particle debug info", input, '=') == 1)
    {
      if(fscanf(ifp, "%s", s_tmp) != 1)
	EH(-1, "Need yes or no for Show particle debug info card.");
      if(!strncmp(s_tmp, "YES", 3))
	{
	  Particle_Show_Debug_Info = 1;
	  printf("Showing particle debug info\n");
	}
    }
  if(!Particle_Show_Debug_Info)
    printf("Not showing particle debug info\n");

  /* Initialization for particle creation and move domains. */
  for(i = 0; i < MAX_DOMAIN_REAL_VALUES; i++)
    {
      Particle_Creation_Domain_Reals[i] = 1.0e+10;
      Particle_Move_Domain_Reals[i] = 1.0e+10;
    }
  Particle_Creation_Domain = UNRESTRICTED;
  Particle_Move_Domain = UNRESTRICTED;
  sprintf(Particle_Creation_Domain_Filename, "<not active>");
  sprintf(Particle_Creation_Domain_Name, "<not active>");
  sprintf(Particle_Move_Domain_Filename, "<not active>");
  sprintf(Particle_Move_Domain_Name, "<not active>");

  /* Set optional creation domain. */
  if(look_for_optional(ifp, "Particle creation domain", input, '=') == 1)
    {
      if(fscanf(ifp, "%s", s_tmp) != 1)
	EH(-1, "Error reading Particle creation domain card.");
      if(!strncmp(s_tmp, "BRICK", 5))
	{
	  Particle_Creation_Domain = BRICK;
	  if(fscanf(ifp, "%lf %lf %lf %lf %lf %lf\n", &Particle_Creation_Domain_Reals[0],
		     &Particle_Creation_Domain_Reals[1], &Particle_Creation_Domain_Reals[2],
		     &Particle_Creation_Domain_Reals[3], &Particle_Creation_Domain_Reals[4],
		     &Particle_Creation_Domain_Reals[5]) != 6)
	    EH(-1, "Error reading min/max values on Particle creation domain card.");
	  printf("Using BRICK creation domain with x in [%g,%g], y in [%g,%g], z in [%g,%g]\n",
		 Particle_Creation_Domain_Reals[0], Particle_Creation_Domain_Reals[1],
		 Particle_Creation_Domain_Reals[2], Particle_Creation_Domain_Reals[3],
		 Particle_Creation_Domain_Reals[4], Particle_Creation_Domain_Reals[5]);
	}
      else if(!strncmp(s_tmp, "ACIS", 4))
	{
#ifdef USE_CGM
	  Particle_Creation_Domain = ACIS_OBJECT;
	  if(fscanf(ifp, "%s", Particle_Creation_Domain_Filename) != 1)
	    EH(-1, "Error reading filename on Particle creation domain card.");
	  if(fscanf(ifp, "%s", Particle_Creation_Domain_Name) != 1)
	    EH(-1, "Error reading geometry name on Particle creation domain card.");
	  printf("Setting Particle_Creation_Domain = ACIS,\n");
	  printf("        Particle_Creation_Domain_Filename = %s,\n", Particle_Creation_Domain_Filename);
	  printf("        Particle_Creation_Domain_Name = %s\n", Particle_Creation_Domain_Name);
#else
	  EH(-1, "Need the CGM library to use ACIS for Particle creation domain card.");
#endif
	}
      else
	{
	  sprintf(err_msg, "Error reading Particle creation domain card, unknown domain type: %s.", s_tmp);
	  EH(-1, err_msg);
	}
    }
  else
    printf("Using unrestricted Particle_Creation_Domain\n");

  /* Set optional move domain. */
  if(look_for_optional(ifp, "Particle move domain", input, '=') == 1)
    {
      if(fscanf(ifp, "%s", s_tmp) != 1)
	EH(-1, "Error reading Particle move domain card.");
      if(!strncmp(s_tmp, "BRICK", 5))
	{
	  Particle_Move_Domain = BRICK;
	  if(fscanf(ifp, "%lf %lf %lf %lf %lf %lf\n", &Particle_Move_Domain_Reals[0],
		    &Particle_Move_Domain_Reals[1], &Particle_Move_Domain_Reals[2],
		    &Particle_Move_Domain_Reals[3], &Particle_Move_Domain_Reals[4],
		    &Particle_Move_Domain_Reals[5]) != 6)
	    EH(-1, "Error reading min/max values on Particle move domain card.");
	  printf("Using BRICK move domain with x in [%g,%g], y in [%g,%g], z in [%g,%g]\n",
		 Particle_Move_Domain_Reals[0], Particle_Move_Domain_Reals[1],
		 Particle_Move_Domain_Reals[2], Particle_Move_Domain_Reals[3],
		 Particle_Move_Domain_Reals[4], Particle_Move_Domain_Reals[5]);
	}
      else if(!strncmp(s_tmp, "ACIS", 4))
	{
#ifdef USE_CGM
	  Particle_Move_Domain = ACIS_OBJECT;
	  if(fscanf(ifp, "%s", Particle_Move_Domain_Filename) != 1)
	    EH(-1, "Error reading filename on Particle move domain card.");
	  if(fscanf(ifp, "%s", Particle_Move_Domain_Name) != 1)
	    EH(-1, "Error reading geometry name on Particle move domain card.");
	  printf("Setting Particle_Move_Domain = ACIS,\n");
	  printf("        Particle_Move_Domain_Filename = %s,\n", Particle_Move_Domain_Filename);
	  printf("        Particle_Move_Domain_Name = %s\n", Particle_Move_Domain_Name);
#else
	  EH(-1, "Need the CGM library to use ACIS for Particle creation domain card.");
#endif
	}
      else
	{
	  sprintf(err_msg, "Error reading Particle move domain card, unknown domain type: %s.", s_tmp);
	  EH(-1, err_msg);
	}
    }
  else
    printf("Using unrestricted Particle_Move_Domain\n");

  /* This is used for restarts, too, and ignores output file type.  If
   * we're a steady-state Goma calculation, then this applies to the
   * particle timesteps taken afterwards.  If this is a transient Goma
   * calculation, then this is a multiple of Goma steps.  In
   * particular, if particles are output on a time strided basis, that
   * will be ignored here. */
  if(look_for_optional(ifp, "Full output stride", input, '=') == 1)
    {
      if(fscanf(ifp, "%d", &Particle_Full_Output_Stride) != 1)
	EH(-1, "Problem reading Full output stride card.");
      if(!fgets(s_tmp, MAX_PARTICLE_FILENAME_LENGTH, ifp))
	EH(-1, "Problem reading Full output stride card.");
      strip(s_tmp);
      strncpy(Particle_Full_Output_Filename, s_tmp, strlen(s_tmp) - 1);
      printf("Setting full output every %d steps to file %s\n", Particle_Full_Output_Stride, Particle_Full_Output_Filename);
    }
  else
    {
      Particle_Full_Output_Stride = 0;
      sprintf(Particle_Full_Output_Filename, "<not active>");
    }

  /* General description for number of particles to track, and how
   * many of them.  Can specify multiple selections, and the particles
   * will be separate samples.  Subject to total number of particles,
   * of course... */
  look_for(ifp, "Number of sample types", input, '=');
  if(fscanf(ifp, "%d", &Particle_Number_Sample_Types) != 1)
    EH(-1, "Problem reading Number of sample types card.");
  printf("Setting Particle_Number_Sample_Types = %d\n", Particle_Number_Sample_Types);

  if(Particle_Number_Sample_Types)
    {
      Particle_Number_Samples_Existing = (int *)array_alloc(1, Particle_Number_Sample_Types, sizeof(int));
      Particle_Number_Samples = (int *)array_alloc(1, Particle_Number_Sample_Types, sizeof(int));
      Particle_Number_Output_Variables = (int *)array_alloc(1, Particle_Number_Sample_Types, sizeof(int));
      Particle_Output_Variables = (particle_variable_s **)
	array_alloc(1, Particle_Number_Sample_Types, sizeof(particle_variable_s *));
      Particle_Filename_Template = (particle_filename_s *)array_alloc(1, Particle_Number_Sample_Types, sizeof(particle_filename_s));
    }
  else
    {
      Particle_Number_Samples_Existing = NULL;
      Particle_Number_Samples = NULL;
      Particle_Number_Output_Variables = NULL;
      Particle_Output_Variables = NULL;
      Particle_Filename_Template = NULL;
    }

  for(i = 0; i < Particle_Number_Sample_Types; i++)
    {
      look_for(ifp, "Sample", input, '=');
      if(fscanf(ifp, "%d", &Particle_Number_Samples[i]) != 1)
	EH(-1, "Error reading number of particle samples.");
      if(Particle_Number_Samples[i] > 0)
	printf("Sample %d contains %d particles", i, Particle_Number_Samples[i]);
      else
	EH(-1, "Number of particle samples must be > 0.");
      
      if(!fgets(s_tmp, SLEN-1, ifp))
	EH(-1, "Error reading sample variables and filename.");
      strncpy(s_tmp_save, s_tmp, SLEN-1);
      Particle_Number_Output_Variables[i] = 0;

      s_ptr1 = strtok(s_tmp, " \t\n");
      while((s_ptr1 = strtok(NULL, " \t\n")))
	Particle_Number_Output_Variables[i]++;

      s_ptr1 = strtok(s_tmp_save, " \t\n");
      if(Particle_Number_Output_Variables[i] > 0)
	{
	  printf(" tracking variable(s)");
	  Particle_Output_Variables[i] = (particle_variable_s *)calloc(Particle_Number_Output_Variables[i], sizeof(particle_variable_s));
	  /*
	    Particle_Output_Variables[i] = (char **)array_alloc(1, Particle_Number_Output_Variables[i], sizeof(char *));
	  */
	  for(j = 0; j < Particle_Number_Output_Variables[i]; j++)
	    {
	      if(strlen(s_ptr1) + 1 > MAX_PARTICLE_OUTPUT_VARIABLE_LENGTH)
		EH(-1, "Increase MAX_PARTICLE_OUTPUT_VARIABLE_LENGTH");
	      strncpy(Particle_Output_Variables[i][j], s_ptr1, strlen(s_ptr1) + 1);
	      printf(" %s", Particle_Output_Variables[i][j]);
	      s_ptr1 = strtok(NULL, " \t\n");
	    }
	}
      
      if(strlen(s_ptr1) + 1 > MAX_PARTICLE_FILENAME_LENGTH)
	EH(-1, "Increase MAX_PARTICLE_FILENAME_LENGTH");
      strncpy(Particle_Filename_Template[i], s_ptr1, strlen(s_ptr1) + 1);
      Particle_Filename_Template[i][strlen(s_ptr1)] = 0;
      printf(" to file %s.\n", Particle_Filename_Template[i]);
    }

  iread = look_for_optional(ifp, "START OF PBC", input, '=');
  if(iread == 1)
    {
      Particle_Number_PBCs = count_list(ifp, "PBC", input, '=', "END OF PBC");
      printf("Found %d PBC's\n", Particle_Number_PBCs);
      PBCs = (PBC_t *)calloc(Particle_Number_PBCs, sizeof(PBC_t));
      for(i = 0; i < Particle_Number_PBCs; i++)
	{
	  PBC = &PBCs[i];
	  look_for(ifp, "PBC", input, '=');
	  if(fscanf(ifp, "%80s", input) != 1)
	    EH(-1, "Error reading PBC.");
	  stringup(input);

	  /* I know it is stupid to do this for an unnecssary token,
	   * but it makes the input look better... */
	  if(fscanf(ifp, "%s", s_tmp) != 1)
	    EH(-1, "Missing SS token on PBC card.");
	  if(strncmp(s_tmp, "SS", 2))
	    EH(-1, "Missing SS token on PBC card.");

	  if(fscanf(ifp, "%d", &(PBC->SS_id)) != 1)
	    EH(-1, "Error reading SS_id on PBC card.");

	  if(!strncmp(input, "OUTFLOW", 7))
	    {
	      PBC->type = PBC_OUTFLOW;
	      puts("Found an OUTFLOW PBC.");
	    }
	  else if(!strncmp(input, "SOURCE", 6))
	    {
	      PBC->type = PBC_SOURCE;
	      if(fscanf(ifp, "%lf", &(PBC->real_data[0])) != 1)
		EH(-1, "Error reading SS specification and float on PBC = SOURCE card.");
	      printf("Found a SOURCE PBC with source = %g.\n",
		     PBC->real_data[0]);
	    }
	  else if(!strncmp(input, "TARGET", 6))
	    {
	      PBC->type = PBC_TARGET;
	      if(fscanf(ifp, "%s", s_tmp) != 1)
		EH(-1, "Missing filename on PBC TARGET card.");
	      if(strlen(s_tmp) + 1 > MAX_PBC_STRING_DATA_LENGTH)
		EH(-1, "Increase MAX_PBC_STRING_DATA_LENGTH");
	      strncpy(PBC->string_data, s_tmp, strlen(s_tmp) + 1);
	      printf("Found a TARGET PBC with filename = '%s'.\n",
		     PBC->string_data);
	    }
	  else if(!strncmp(input, "FREESTREAM_SOURCE", 17))
	    {
	      PBC->type = PBC_FREESTREAM_SOURCE;
	      if(fscanf(ifp, "%lf", &(PBC->real_data[0])) != 1)
		EH(-1, "Error reading SS specification and float on PBC = FREESTREAM_SOURCE card.");
	      printf("Found a FREESTREAM_SOURCE PBC with source = %g.\n",
		     PBC->real_data[0]);
	    }
	  else if(!strncmp(input, "IMPERMEABLE", 11))
	    {
	      PBC->type = PBC_IMPERMEABLE;
	      if(fscanf(ifp, "%lf", &(PBC->real_data[0])) != 1)
		EH(-1, "Error reading distance factor on PBC = IMPERMEABLE card.");
	      if(!fgets(s_tmp, SLEN-1, ifp))
		EH(-1, "Error reading geometry entity name on PBC = IMPERMEABLE card.");
	      strip(s_tmp);
	      if(s_tmp[strlen(s_tmp)-1] == '\n')
		s_tmp[strlen(s_tmp)-1] = 0;
	      if(strlen(s_tmp) + 1 > MAX_PBC_STRING_DATA_LENGTH)
		EH(-1, "Increase MAX_PBC_STRING_DATA_LENGTH");
	      strncpy(PBC->string_data, s_tmp, strlen(s_tmp) + 1);
	      printf("Found an IMPERMEABLE PBC with distance factor = %g to entity '%s'.\n",
		     PBC->real_data[0], PBC->string_data);
	    }
	  else
	    EH(-1, "Error reading PBC type.");
	}
    }
  else
    {
      Particle_Number_PBCs = 0;
      PBCs = NULL;
    }
}
