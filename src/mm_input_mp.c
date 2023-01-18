/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2022 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/

/*
 *$Id: mm_input_mp.c,v 5.37 2010-07-30 21:14:52 prschun Exp $
 */

/*************** R O U T I N E S   I N   T H E   F I L E ***********************
 *
 *    NAME				TYPE		CALLED_BY
 *--------------------------------------------------------------------
 *
 *    look_for_optional_string          int
 *    rd_mp_specs	       		void		read_input_file
 ******************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h> /* for strcasecmp */

#include "el_elm.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "mm_fill_em.h"
#include "mm_input.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "mm_post_def.h"
#include "rf_allo.h"
#include "rf_bc.h"
#include "rf_bc_const.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_io_const.h"
#include "std.h"

#define GOMA_MM_INPUT_C

#define NO_USER      NULL
#define NO_INPUT     0
#define SCALAR_INPUT 1
#define VECTOR_INPUT 3

extern Spfrtn sr; /* External declaration for the sprintf return variable,sr
                   *  sr is defined in mm_input.c */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

static void fallback_chemkin_generic_prop(int *model_read,
                                          int species_num,
                                          int *propID,
                                          int canBeGenericallyHandled,
                                          MATRL_PROP_STRUCT *mat_ptr)

/****************************************************************************
 *
 * fallback_chemkin_generic_prop
 *
 *        This routine will check to see whether an unspecified entry in
 *  the Goma's material file can be filled in with the generic CHEMKIN_MODEL
 *  property. If it can be used, then the model_read parameter is changed
 *  to "satisfied", and the property ID is set to CHEMKIN_MODEL if it has not
 *  previously been set nonzero by the chemkin initialization routine.
 *  The property ID is left alone, if it has been previously set by the
 *  chemkin initialization routine.  CHEMKIN_MODEL is a generic chemkin
 *  constitutive model for the current property.
 *        The routine checks to see whether the default database is set to
 *  chemkin to make this determination.
 *        This routine also checks to see whether it is called for the
 *  last species in the mechanism. If the mechanism is nondilute so that the
 *  number of species is different than the number of species equations, it
 *  changes the value of the model_read so as not to generate an error
 *  call on return. The parameters for the "non-condensible", i.e., last
 *  species in the mechanism are input via special cases later on. Therefore,
 *  we don't want to create an error condition here.
 *
 *       This code snipet also handles the case of input for nondilute
 *  formulations, where the number of species is 1 larger than the number
 *  of species equations. In this case, sometimes the property input for the
 *  last species in the mechanism (the one whose continuity equation is not
 *  part of the solutin set) is not specified. For backwards compatibility,
 *  we must allow this to be OK. (((.
 *
 *  Input
 * --------
 *  *model_read    = Address of the model id read in from the input
 *                   deck. If none was found this will be equal to -1.
 *                   This routine only considers the case where no
 *                   model was found.
 *  *species_num   = pointer to the species index in question.
 *  *propID        = Address in the integer property arrray corresponding
 *                   to the Material property constituitive equation model
 *                   for the species above.
 *  canBeGenericallyHandled =
 *                   This value is either true or false. If true it means,
 *                   for this particular property, that chemkin can handle
 *                   the calculation even if the propID wasn't initially
 *                   set by the chemkin initialization process.
 *  mat_ptr        = Pointer to the material structure.
 *
 *  Output
 * ---------
 * *model_read     = If chemkin can handle the missing property, this value
 *                   is changed from -1 to 1.
 *                   (this is also set to 1 if we are addressing the last
 *                    species in a nondilute material, where the number
 *                    of species equations is one less than the total
 *                    number of equations)
 * *propID         = If chemkin can handle a missing property, and the
 *                   property model value has not been previously specified
 *                   the property model is set to CHEMKIN_MODEL.
 *                   (this is also set to 0 if we are addressing the last
 *                    species in a nondilute material, where the number
 *                    of species equations is one less than the total
 *                    number of equations)
 ****************************************************************************/
{
  if (*model_read == -1) {
    if (mat_ptr->DefaultDatabase == DB_CHEMKIN_MAT) {
      if ((*propID == UNINITIALIZED_MODEL) || (*propID == NO_MODEL)) {
        if (canBeGenericallyHandled) {
          *propID = CHEMKIN_MODEL;
          *model_read = 1;
        }
      } else {
        *model_read = 1;
      }
    } else {
      if (mat_ptr->Num_Species_Eqn != mat_ptr->Num_Species) {
        if (species_num == mat_ptr->Num_Species_Eqn) {
          *model_read = 1;
          *propID = UNINITIALIZED_MODEL;
        }
      }
    }
  }
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int look_for_optional_string(FILE *ifp,
                             const char *match_string,
                             char retn_string[],
                             int retn_string_len)

/***************************************************************************
 *
 * look_for_optional_string:
 *
 *	Scan the input file (reading in strings according to 'read_string(ifp,)'
 *	specifications) until the character pattern in 'string' is matched.
 *	If 'string' is not matched, return a value of -1.
 *	If 'string' is matched, return the number of characters read
 *  into the return string not including the trailing null character.
 *
 *  The string value after an equals sign is returned in retn_string upon
 *  successful search of string. On an unsuccessful return, the
 *  retn_string is set to the null string.
 *
 *  This search starts at the current position in the input
 *  file and searches to the end until it finds the input string.
 *  If this routine can't find the input string, it returns to the
 *  file position where the search was started. If it finds the
 *	string pattern, it leaves the file pointer positioned after the
 *  new line character at the end of matched line.
 *
 *  Input
 * --------
 *  ifp          = file pointer to file "input"
 *  match_string = contains string pattern to be matched.
 *  retn_string_len = max string length of retn_string permissible
 *
 *  Output
 * ---------
 *  retn_string  = On return it contains the string after the
 *                 equals sign on a line containing the match_string
 *                 If the match_string is not found, this string
 *                 is not touched.
 *
 *  Example:
 *
 *   Line:
 *     Default Model = Any which way but loose
 *
 *   Calling statement
 *
 *     look_for_optional_string(ifp,"Default Model", retn_string);
 *
 *   On return, rntn_string will be equal to "Any which way but loose" and
 *   return value will be equal to 23
 ***************************************************************************/
{
  int retn;
  char tmp_string[MAX_CHAR_IN_INPUT];
  if ((retn_string_len < 1) || (!match_string) || (!ifp)) {
    GOMA_EH(GOMA_ERROR, "look_for_optional_string interface ERROR");
  }
  retn = look_forward_optional(ifp, match_string, tmp_string, '=');
  if (retn == 1) {
    retn = read_string(ifp, tmp_string, '\n');
    strip(tmp_string);
    (void)strncpy(retn_string, tmp_string, retn_string_len);
    retn_string[retn_string_len - 1] = '\0';
    retn = strlen(retn_string);
  }
  return retn;
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void rd_mp_specs(FILE *imp, char input[], int mn, char *echo_file)

/***************************************************************************
 *
 * rd_mp_specs -- read material properties specifications
 *
 * Comments:
 *
 * Input
 * -----------
 *
 * imp = File pointer for the material file for this material.
 * input = buffer array used to read the input file.
 *         MAX_CHAR_IN_INPUT long (256 chars in length, currently)
 * mn = Index of the current material in the list of materials defined for
 *      this problem.
 **************************************************************************/
{ /*start*/
  char err_msg[MAX_CHAR_ERR_MSG];
  int i, j, var;
  int imtrx;
  int ipore;
  static const char yo[] = "rd_mp_specs";
  struct Elastic_Constitutive *dum_ptr;

  int ConstitutiveEquation;
  int LameLambdaModel;
  dbl a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11;
  dbl v0[DIM];
  int i0;
  /* dummy variable to hold modal data before it is put into the ve struct */
  dbl *modal_data;
  int DiffusionConstitutiveEquation = -1;
  int PorousDiffusionConstitutiveEquation = -1;
  int ElasticConstitutiveEquation = -1;
  int PlasticConstitutiveEquation = -1;
  int MomentumSourceModel = -1;
  int HeatSourceModel = -1;
  int SpeciesSourceModel = -1;
  int SpeciesTimeIntegration;
  int PorousTimeIntegration;
  int Num_liquid_phase_components = 1;        /*PRS: hardwired for Porous (see
                                                below) 052901 */
  int Num_insoluble_gas_phase_components = 1; /*PRS: hardwired for Porous (see
                                                below) 052901 */
  int iread;
  int num_const, species_no;
  int porous_no, porous_shell_no;
  int PorousShellOn;
  int mm;                          /* modal counter */
  int ii, jj, n_dij, k, n_species; /* KSC: 7/98 */
  dbl dij, E, T0;                  /* KSC: 7/98, 9/04 */
  int n_rxn;                       /* KSC/GHE: 10/98 */

  int have_mesh_eqn;
  int have_por_liq_pres;
  int have_por_gas_pres;
  int have_por_porosity;
  int have_por_sink_mass;
  int have_por_energy;
  int have_shell_sat_open;
  int have_shell_sat_open2;
  int have_shear_rate;
  /* int have_vort_dir; */

  /*
   *  Pointers to Material property structures. "matr_ptr" points to the
   *  current structure that we will be filling up in this routine.
   *  mp_glob is a global vector of pointers to material prop structures.
   */
  /*  extern MATRL_PROP_STRUCT **mp_glob; cf mm_mp.h */
  MATRL_PROP_STRUCT *mat_ptr = mp_glob[mn];
  PROBLEM_DESCRIPTION_STRUCT *pd_ptr = pd_glob[mn];

  /*  fpos_t file_position; position in file at start of search */

  char model_name[MAX_CHAR_IN_INPUT];
  char *s; /* used to tokenize optional input string. */
  char echo_string[MAX_CHAR_ECHO_INPUT] = "\0";
  char search_string[MAX_CHAR_IN_INPUT];
  char *es = echo_string;

  int model_read, retn;
  int matl_model = 0;
  int NO_SPECIES = -1; /* Signal to parsing routine to not
                        * expect a species num                    */
  int n_ij, read_bc_mp = -1;
  dbl chi_ij, mw, mv;

  /*
   *  Copy the materials name into the materials structure from the
   *  problem structure
   */
  (void)strncpy(mat_ptr->Material_Name, pd_glob[mn]->MaterialName, MAX_MATLNAME);

  /*
   *  Copy the number of species and species equations into the material
   *  property structure. The number of species in each material may vary
   *  in the future.
   */
  mat_ptr->Num_Species_Eqn = pd_glob[mn]->Num_Species_Eqn;
  mat_ptr->Num_Species = pd_glob[mn]->Num_Species;
  mat_ptr->Num_Porous_Eqn = pd_glob[mn]->Num_Porous_Eqn;
  mat_ptr->Num_Porous_Shell_Eqn = pd_glob[mn]->Num_Porous_Shell_Eqn;

  /*
   *  Intialize to good default behavior
   */

  for (i = 0; i < mat_ptr->Num_Species; i++) {
    mat_ptr->SpeciesTimeIntegration[i] = STANDARD;
    mat_ptr->AdvectiveScalingModel[i] = CONSTANT;
    mat_ptr->ExtrinsicIndependentSpeciesVar[i] = 0;
    mat_ptr->AdvectiveScaling[i] = 1.0;
    mat_ptr->FreeVolSolvent[i] = TRUE;
  }

  /*
   *  Database Location Section --
   *
   *       This section assigns the DefaultDatabase  int variable in the
   *       material  structure.
   */
  mat_ptr->DefaultDatabase = DB_GOMA_MAT;
  retn = look_for_optional_string(imp, "Default Database", input, MAX_CHAR_IN_INPUT);
  if (retn > 0) {
    if (!strcasecmp(input, "chemkin_mat"))
      mat_ptr->DefaultDatabase = DB_CHEMKIN_MAT;
    else if (!strcasecmp(input, "goma_mat"))
      mat_ptr->DefaultDatabase = DB_GOMA_MAT;
    else {
      sr = sprintf(err_msg, "Unknown value for Default Database: %s, use chemkin_mat or goma_mat\n",
                   input);
      GOMA_EH(GOMA_ERROR, err_msg);

      SPF(es, "%s = %s", "Default Database", input);
      ECHO(es, echo_file);
    }
#ifndef USE_CHEMKIN
    if (!strcasecmp(input, "chemkin_mat")) {
      fprintf(stderr, "ERROR! The default database has been specified as chemkin.\n");
      fprintf(stderr, "\tHowever, CHEMKIN has not been linked in!\n");
      fprintf(stderr, "\tGOMA must be recompiled with the USE_CHEMKIN definition!\n");
      GOMA_EH(GOMA_ERROR, "chemkin not linked in\n");
    }
#endif
  }

  /*
   *  If we have a Chemkin material, then lets
   *  1) read in the chemkin database,  if it hasn't already been read in
   *  2) Initialize this material structure with the chemkin properties
   *  3) check for inconsistencies
   *
   *  For Goma Default Databases, lets do initializations dependent upon
   *  the number of species in this material.
   */
  if (mat_ptr->DefaultDatabase == DB_CHEMKIN_MAT) {
#ifdef USE_CHEMKIN
    retn = chemkin_mat_prop_init(mat_ptr, mn, pd_ptr);
    if (retn < 0) {
      sr = sprintf(err_msg, "%s%d had an inconsistency, BAIL!/n",
                   "Chemkin Material property specification for material ", mn);
      GOMA_EH(GOMA_ERROR, err_msg);
    }
#endif
  } else {
    retn = goma_mat_prop_init(mat_ptr, mn, pd_ptr);
    if (retn < 0) {
      sprintf(err_msg, "%s for mat %d had an inconsistency, BAIL!/n", "goma_mat_prop_int", mn);
      GOMA_EH(GOMA_ERROR, err_msg);
    }
  }

  /* Assume no second level set phase */

  if (ls != NULL)
    mat_ptr->mp2nd =
        (SECOND_LS_PHASE_PROP_STRUCT *)alloc_void_struct_1(sizeof(SECOND_LS_PHASE_PROP_STRUCT), 1);
  else
    mat_ptr->mp2nd = NULL;

  /*
   * Density section
   *
   *
   *        Read the "Density" line in the materials property data file.
   *        this will  handle the generic cases of CONSTANT, USER, and USER_GEN
   */

  ECHO("\n---Density\n", echo_file);

  model_name[0] = '\0';
  model_read = look_for_mat_prop(imp, "Density", &(mat_ptr->DensityModel), &(mat_ptr->density),
                                 &(mat_ptr->u_density), &(mat_ptr->len_u_density), model_name,
                                 SCALAR_INPUT, &NO_SPECIES, es);

  /*
   *    Handle specific density models, finish reading the input line.
   */
  if (model_read == -1 && !strcmp(model_name, "FILL")) {
    mat_ptr->DensityModel = DENSITY_FILL;
    num_const = read_constants(imp, &(mat_ptr->u_density), 0);
    if (num_const < 2) {
      sprintf(err_msg, "Material %s - expected at least 2 constants for %s %s model.\n",
              pd_ptr->MaterialName, "Density", "FILL");
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    mat_ptr->len_u_density = num_const;

    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_density);
  } else if (model_read == -1 && !strcmp(model_name, "SUSPENSION")) {
    mat_ptr->DensityModel = DENSITY_SUSPENSION;
    num_const = read_constants(imp, &(mat_ptr->u_density), 0);
    if (num_const < 3) {
      sprintf(err_msg, "Material %s - expected at least 3 constants for %s %s model.\n",
              pd_glob[mn]->MaterialName, "Density", "SUSPENSION");
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    mat_ptr->len_u_density = num_const;

    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_density);
  } else if (model_read == -1 && !strcmp(model_name, "SOLVENT_POLYMER")) {
    mat_ptr->DensityModel = SOLVENT_POLYMER;
    num_const = read_constants(imp, &(mat_ptr->u_density), 0);
    if (num_const < 1) {
      sprintf(err_msg, "Material %s - expected at least 1 constants for %s %s model.\n",
              pd_glob[mn]->MaterialName, "Density", "SOLVENT_POLYMER");
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    mat_ptr->len_u_density = num_const;
    mat_ptr->specific_volume[pd_glob[mn]->Num_Species_Eqn] = mat_ptr->u_density[0];

    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_density);

  } else if (model_read == -1 && !strcmp(model_name, "REACTIVE_FOAM")) {
    mat_ptr->DensityModel = REACTIVE_FOAM;
    num_const = read_constants(imp, &(mat_ptr->u_density), 0);
    if (num_const < 1) {
      sprintf(err_msg, "Material %s - expected at least 1 constants for %s %s model.\n",
              pd_glob[mn]->MaterialName, "Density", "REACTIVE_FOAM");
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    mat_ptr->len_u_density = num_const;
    mat_ptr->specific_volume[pd_glob[mn]->Num_Species_Eqn] = mat_ptr->u_density[0];
    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_density);
  }
  /* MMH */
  else if (model_read == -1 && !strcmp(model_name, "SUSPENSION_PM")) {
    mat_ptr->DensityModel = DENSITY_SUSPENSION_PM;
    num_const = read_constants(imp, &(mat_ptr->u_density), 0);
    if (num_const < 3) {
      sprintf(err_msg, "Material %s - expected at least 3 constants for %s %s model.\n",
              pd_glob[mn]->MaterialName, "Density", "SUSPENSION_PM");
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    mat_ptr->len_u_density = num_const;
    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_density);
  } else if (model_read == -1 && !strcmp(model_name, "THERMAL_BATTERY")) {
    mat_ptr->DensityModel = DENSITY_THERMAL_BATTERY;
    num_const = read_constants(imp, &(mat_ptr->u_density), 0);
    if (num_const < 2) {
      sprintf(err_msg, "Material %s - expected at least 2 constants for %s %s model.\n",
              pd_glob[mn]->MaterialName, "Density", "THERMAL_BATTERY");
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    mat_ptr->len_u_density = num_const;
    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_density);
  } else if (model_read == -1 && !strcmp(model_name, "IDEAL_GAS")) {
    mat_ptr->DensityModel = DENSITY_IDEAL_GAS;
    num_const = read_constants(imp, &(mat_ptr->u_density), 0);
    if (num_const < 2) {
      sprintf(err_msg, "Material %s - expected at least 2 constants for %s %s model.\n",
              pd_glob[mn]->MaterialName, "Density", "IDEAL_GAS");
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    mat_ptr->len_u_density = num_const;
    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_density);
  } else if (model_read == -1 && !strcmp(model_name, "LEVEL_SET")) {
    mat_ptr->DensityModel = DENSITY_LEVEL_SET;
    num_const = read_constants(imp, &(mat_ptr->u_density), 0);
    if (num_const < 3) {
      sprintf(err_msg, "Material %s - expected at least 3 constants for %s %s model.\n",
              pd_glob[mn]->MaterialName, "Density", "LEVEL_SET");
      GOMA_EH(GOMA_ERROR, err_msg);
    }

    if (mat_ptr->u_density[2] == 0.0)
      mat_ptr->u_density[2] = ls->Length_Scale / 2.0;

    mat_ptr->len_u_density = num_const;
    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_density);
  } else if (model_read == -1 && !strcmp(model_name, "CONST_PHASE_FUNCTION")) {
    mat_ptr->DensityModel = DENSITY_CONST_PHASE_FUNC;
    num_const = read_constants(imp, &(mat_ptr->u_density), 0);

    if (num_const < pfd->num_phase_funcs + 2) {
      sprintf(err_msg, "Material %s - expect at least %d constants for %s %s model.\n",
              pd_glob[mn]->MaterialName, pfd->num_phase_funcs + 2, "Density",
              "CONST_PHASE_FUNCTION");
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    mat_ptr->len_u_density = num_const;
    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_density);
  } else if (model_read == -1 && !strcmp(model_name, "FOAM")) {
    mat_ptr->DensityModel = DENSITY_FOAM;
    num_const = read_constants(imp, &(mat_ptr->u_density), 0);
    if (num_const < 6) {
      sprintf(err_msg, "Material %s - expected at least 6 constants for %s %s model.\n",
              pd_glob[mn]->MaterialName, "Density", "FOAM");
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    mat_ptr->len_u_density = num_const;
    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_density);
  } else if (model_read == -1 && !strcmp(model_name, "FOAM_PMDI_10")) {
    mat_ptr->DensityModel = DENSITY_FOAM_PMDI_10;
    num_const = read_constants(imp, &(mat_ptr->u_density), 0);
    if (num_const < 3) {
      sprintf(err_msg, "Material %s - expected at least 3 constants for %s %s model.\n",
              pd_glob[mn]->MaterialName, "Density", "FOAM_PMDI_10");
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    mat_ptr->len_u_density = num_const;
    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_density);
  } else if (model_read == -1 && !strcmp(model_name, "MOMENT_BASED")) {
    mat_ptr->DensityModel = DENSITY_MOMENT_BASED;
    num_const = read_constants(imp, &(mat_ptr->u_density), 0);
    if (num_const < 2) {
      sprintf(err_msg, "Material %s - expected at least 2 constants for %s %s model.\n",
              pd_glob[mn]->MaterialName, "Density", "FOAM_PMDI_10");
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    mat_ptr->len_u_density = num_const;
    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_density);
  } else if (model_read == -1 && !strcmp(model_name, "FOAM_CONC")) {
    mat_ptr->DensityModel = DENSITY_FOAM_CONC;
    num_const = read_constants(imp, &(mat_ptr->u_density), 0);
    if (num_const < 8) {
      sprintf(err_msg, "Material %s - expected at least 8 constants for %s %s model.\n",
              pd_glob[mn]->MaterialName, "Density", "FOAM_CONC");
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    mat_ptr->len_u_density = num_const;
    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_density);
  } else if (model_read == -1 && !strcmp(model_name, "FOAM_TIME")) {
    mat_ptr->DensityModel = DENSITY_FOAM_TIME;
    num_const = read_constants(imp, &(mat_ptr->u_density), 0);
    if (num_const < 4) {
      sprintf(err_msg, "Material %s - expected at least 4 constants for %s %s model.\n",
              pd_glob[mn]->MaterialName, "Density", "FOAM_TIME");
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    mat_ptr->len_u_density = num_const;
    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_density);
  } else if (model_read == -1 && !strcmp(model_name, "FOAM_TIME_TEMP")) {
    mat_ptr->DensityModel = DENSITY_FOAM_TIME_TEMP;
    num_const = read_constants(imp, &(mat_ptr->u_density), 0);
    if (num_const < 5) {
      sprintf(err_msg, "Material %s - expected at least 5 constants for %s %s model.\n",
              pd_glob[mn]->MaterialName, "Density", "FOAM_TIME_TEMP");
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    mat_ptr->len_u_density = num_const;
    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_density);
  } else if (model_read == -1 && !strcmp(model_name, "CONSTANT_LAST_CONC")) {
    mat_ptr->DensityModel = DENSITY_CONSTANT_LAST_CONC;
    num_const = read_constants(imp, &(mat_ptr->u_density), 0);
    if (num_const < 1) {
      sprintf(err_msg, "Material %s - expected at least 1 constants for %s %s model.\n",
              pd_glob[mn]->MaterialName, "Density", "CONSTANT_LAST_CONC");
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    mat_ptr->len_u_density = num_const;
    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_density);
  } else if (model_read == -1 && !strcmp(model_name, "FOAM_PBE")) {
    mat_ptr->DensityModel = DENSITY_FOAM_PBE;
    num_const = read_constants(imp, &(mat_ptr->u_density), 0);
    if (num_const != 3) {
      sprintf(err_msg, "Material %s - expected 3 constants for %s %s model.\n",
              pd_glob[mn]->MaterialName, "Density", "FOAM_PBE");
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    mat_ptr->len_u_density = num_const;
    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_density);
  } else if (model_read == -1 && !strcmp(model_name, "FOAM_PBE_EQN")) {
    mat_ptr->DensityModel = DENSITY_FOAM_PBE_EQN;
    num_const = read_constants(imp, &(mat_ptr->u_density), 0);
    if (num_const != 3) {
      sprintf(err_msg, "Material %s - expected 3 constants for %s %s model.\n",
              pd_glob[mn]->MaterialName, "Density", "FOAM_PBE");
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    mat_ptr->len_u_density = num_const;
    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_density);
  } else {
    sprintf(err_msg, "Material %s - unrecognized model for %s \"%s\" ???\n",
            pd_glob[mn]->MaterialName, "Density", model_name);
    GOMA_EH(model_read, err_msg);
  }

  ECHO(es, echo_file);

  model_read = look_for_mat_prop(imp, "Second Level Set Density", &(i0), &(a0), NO_USER, NULL,
                                 model_name, SCALAR_INPUT, &NO_SPECIES, es);

  if (model_read != -1) {
    if (ls == NULL)
      GOMA_EH(GOMA_ERROR, "Second Level Set Density requires activation of Level Set Tracking.\n");

    mat_ptr->mp2nd->DensityModel = i0;
    mat_ptr->mp2nd->density = a0;

    stringup(model_name);

    if (!strcmp(model_name, "CONSTANT")) {
      if (fscanf(imp, "%s", input) != 1) {
        GOMA_EH(GOMA_ERROR, "Expecting trailing keyword for Second Level Set Density.\n");
      }

      stringup(input);

      if (strncmp(input, "POSITIVE", 3) == 0) {
        mat_ptr->mp2nd->densitymask[0] = 0;
        mat_ptr->mp2nd->densitymask[1] = 1;
      } else if (strncmp(input, "NEGATIVE", 3) == 0) {
        mat_ptr->mp2nd->densitymask[0] = 1;
        mat_ptr->mp2nd->densitymask[1] = 0;
      } else {
        GOMA_EH(GOMA_ERROR, "Keyword must be POSITIVE or NEGATIVE for Second Level Set Density.\n");
      }
      SPF(endofstring(es), " %s", input);
      if (pfd != NULL) {
        for (i = 0; i < pfd->num_phase_funcs; i++) {
          if (fscanf(imp, "%lf", &(mat_ptr->mp2nd->density_phase[i])) != 1) {
            GOMA_EH(GOMA_ERROR, "error reading phase density");
          }
          SPF(endofstring(es), " %g", mat_ptr->mp2nd->density_phase[i]);
        }
      }

    } else {
      GOMA_EH(GOMA_ERROR, "Second Level Set Density model can only be CONSTANT.\n");
    }
  }

  ECHO(es, echo_file);

  /*
   * Solid Constitutive Equation
   */
  ECHO("\n----Solid Constitutive Equation\n", echo_file);

  model_read = look_for_mat_prop(imp, "Solid Constitutive Equation", &(ElasticConstitutiveEquation),
                                 &a0, NO_USER, NULL, model_name, NO_INPUT, &NO_SPECIES, es);
  if (model_read == 1)
    GOMA_EH(GOMA_ERROR, "Can't have USER or CONSTANT for Solid Constitutive Equation");

  if (!strcmp(model_name, "LINEAR")) {
    ElasticConstitutiveEquation = LINEAR;
  } else if (!strcmp(model_name, "NONLINEAR") || !strcmp(model_name, "NONLINEAR_PLANE_STRAIN")) {
    ElasticConstitutiveEquation = NONLINEAR;
  } else if (!strcmp(model_name, "INCOMP_PSTRAIN")) {
    ElasticConstitutiveEquation = INCOMP_PSTRAIN;
  } else if (!strcmp(model_name, "INCOMP_PSTRESS")) {
    ElasticConstitutiveEquation = INCOMP_PSTRESS;
  } else if (!strcmp(model_name, "HOOKEAN_PSTRAIN")) {
    ElasticConstitutiveEquation = HOOKEAN_PSTRAIN;
  } else if (!strcmp(model_name, "HOOKEAN_PSTRESS")) {
    ElasticConstitutiveEquation = HOOKEAN_PSTRESS;
  } else if (!strcmp(model_name, "INCOMP_3D")) {
    ElasticConstitutiveEquation = INCOMP_3D;
  } else if (!strcmp(model_name, "KELVIN_VOIGT")) {
    ElasticConstitutiveEquation = KELVIN_VOIGT;
  } else /* default to nonlinear */
  {
    ElasticConstitutiveEquation = NONLINEAR;
    SPF(endofstring(es), "\t(%s)", "Solid Constitutive Equation defaulting to NONLINEAR");
  }
  /*This seems a little redundant! prs (2/24/95) */
  cr_glob[mn]->MeshFluxModel = ElasticConstitutiveEquation;

  ECHO(es, echo_file);

  model_read = look_for_mat_prop(imp, "Plasticity Equation", &(PlasticConstitutiveEquation), &a0,
                                 NO_USER, NULL, model_name, NO_INPUT, &NO_SPECIES, es);

  if (model_read == 1)
    GOMA_EH(GOMA_ERROR, "Can't have USER or CONSTANT for Plasticity Constitutive Equation");
  if (!strcmp(model_name, "EVP_HYPER")) {
    PlasticConstitutiveEquation = EVP_HYPER;
  } else {
    PlasticConstitutiveEquation = NO_MODEL;
  }
  evpl_glob[mn]->ConstitutiveEquation = PlasticConstitutiveEquation;

  ECHO(es, echo_file);

  /* An optional additional type of mesh motion which includes an inertial term in the
     momentum equation */

  if (pd_glob[mn]->MeshMotion == LAGRANGIAN || pd_glob[mn]->MeshMotion == DYNAMIC_LAGRANGIAN ||
      pd_glob[mn]->MeshMotion == TOTAL_ALE) {
    model_read = look_for_mat_prop(imp, "Convective Lagrangian Velocity",
                                   &(pd_glob[mn]->MeshInertia), elc_glob[mn]->v_mesh_sfs, NO_USER,
                                   NULL, model_name, VECTOR_INPUT, &NO_SPECIES, es);
    if (model_read == -1) {
      pd_glob[mn]->MeshInertia = 0;
      elc_glob[mn]->v_mesh_sfs_model = 0;
      if (!strcmp(model_name, "ROTATIONAL")) {
        pd_glob[mn]->MeshInertia = 1;
        elc_glob[mn]->v_mesh_sfs_model = ROTATIONAL;

        num_const = read_constants(imp, &(elc_glob[mn]->u_v_mesh_sfs), NO_SPECIES);

        if (num_const < 4) {
          sr = sprintf(err_msg, "Matl %s expected at least 4 constants for %s %s model.\n",
                       pd_glob[mn]->MaterialName, "Convective Lagrangian Velocity", "ROTATIONAL");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
        elc_glob[mn]->len_u_v_mesh_sfs = num_const;

        SPF_DBL_VEC(endofstring(es), num_const, elc_glob[mn]->u_v_mesh_sfs);
      }
      if (!strcmp(model_name, "ROTATIONAL_3D")) {
        pd_glob[mn]->MeshInertia = 1;
        elc_glob[mn]->v_mesh_sfs_model = ROTATIONAL_3D;

        num_const = read_constants(imp, &(elc_glob[mn]->u_v_mesh_sfs), NO_SPECIES);

        if (num_const < 7) {
          sr =
              sprintf(err_msg, "Matl %s expected at least 7 constants for %s %s model.\n",
                      pd_glob[mn]->MaterialName, "Convective Lagrangian Velocity", "ROTATIONAL_3D");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
        elc_glob[mn]->len_u_v_mesh_sfs = num_const;

        SPF_DBL_VEC(endofstring(es), num_const, elc_glob[mn]->u_v_mesh_sfs);
      }
      if (!strcmp(model_name, "OSC_LINEAR")) {
        pd_glob[mn]->MeshInertia = 1;
        elc_glob[mn]->v_mesh_sfs_model = OSC_LINEAR;

        num_const = read_constants(imp, &(elc_glob[mn]->u_v_mesh_sfs), NO_SPECIES);

        if (num_const < 5) {
          sr = sprintf(err_msg, "Matl %s expected at least 5 constants for %s %s model.\n",
                       pd_glob[mn]->MaterialName, "Convective Lagrangian Velocity", "OSC_LINEAR");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
        elc_glob[mn]->len_u_v_mesh_sfs = num_const;

        SPF_DBL_VEC(endofstring(es), num_const, elc_glob[mn]->u_v_mesh_sfs);
      }
      if (!strcmp(model_name, "NONE")) {
        pd_glob[mn]->MeshInertia = 0;
        elc_glob[mn]->v_mesh_sfs_model = 0;
        elc_glob[mn]->v_mesh_sfs[0] = 0.;
        elc_glob[mn]->v_mesh_sfs[1] = 0.;
        elc_glob[mn]->v_mesh_sfs[2] = 0.;
      }

    } else {
      elc_glob[mn]->v_mesh_sfs_model = pd_glob[mn]->MeshInertia;
      if (pd_glob[mn]->TimeIntegration == TRANSIENT) {
        GOMA_WH(-1, "Currently can't do transient mesh motion with inertia");
      }
    }

    ECHO(es, echo_file);

  } else {
    pd_glob[mn]->MeshInertia = 0;
  }

  /* read in constants for constitutive equation if they are input */

  model_read = look_for_mat_proptable(imp, "Lame MU", &(elc_glob[mn]->lame_mu_model),
                                      &(elc_glob[mn]->lame_mu), &(elc_glob[mn]->u_mu),
                                      &(elc_glob[mn]->len_u_mu), &(elc_glob[mn]->lame_mu_tableid),
                                      model_name, SCALAR_INPUT, &NO_SPECIES, es);
  if (model_read == -1) {
    if (!strcmp(model_name, "POWER_LAW")) {
      elc_glob[mn]->lame_mu_model = POWER_LAW;
      num_const = read_constants(imp, &(elc_glob[mn]->u_mu), NO_SPECIES);
      if (num_const < 3) {
        sr = sprintf(err_msg, "Matl %s expected at least 3 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Lame MU", "POWER_LAW");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      elc_glob[mn]->len_u_mu = num_const;

    } else if (!strcmp(model_name, "CONTACT_LINE")) {
      elc_glob[mn]->lame_mu_model = CONTACT_LINE;
      num_const = read_constants(imp, &(elc_glob[mn]->u_mu), NO_SPECIES);
      if (num_const < 4) {
        sr = sprintf(err_msg, "Matl %s expected at least 4 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Lame MU", "CONTACT_LINE");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      elc_glob[mn]->len_u_mu = num_const;
    } else if (!strcmp(model_name, "SHEAR_HARDEN")) {
      elc_glob[mn]->lame_mu_model = SHEAR_HARDEN;
      num_const = read_constants(imp, &(elc_glob[mn]->u_mu), NO_SPECIES);
      if (num_const < 2) {
        sr = sprintf(err_msg, "Matl %s expected at least 2 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Lame MU", "SHEAR_HARDEN");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      elc_glob[mn]->len_u_mu = num_const;
    } else if (!strcmp(model_name, "EXPONENTIAL")) {
      elc_glob[mn]->lame_mu_model = EXPONENTIAL;
      num_const = read_constants(imp, &(elc_glob[mn]->u_mu), NO_SPECIES);
      if (num_const < 3) {
        sr = sprintf(err_msg, "Matl %s expected at least 3 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Lame MU", "EXPONENTIAL");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      elc_glob[mn]->len_u_mu = num_const;
    } else if (!strcmp(model_name, "DENSE_POWER_LAW")) {
      elc_glob[mn]->lame_mu_model = DENSE_POWER_LAW;
      num_const = read_constants(imp, &(elc_glob[mn]->u_mu), NO_SPECIES);
      if (num_const < 2) {
        sr = sprintf(err_msg, "Matl %s expected at least 2 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Lame MU", "DENSE_POWER_LAW");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      elc_glob[mn]->len_u_mu = num_const;
    } else {
      sr = sprintf(err_msg, "Material %s - unrecognized model for %s \"%s\" ???\n",
                   pd_glob[mn]->MaterialName, "Lame MU", model_name);
      GOMA_EH(model_read, err_msg);
    }

    SPF_DBL_VEC(endofstring(es), num_const, elc_glob[mn]->u_mu);
  }

  ECHO(es, echo_file);

  model_read =
      look_for_mat_prop(imp, "Lame LAMBDA", &(elc_glob[mn]->lame_lambda_model),
                        &(elc_glob[mn]->lame_lambda), &(elc_glob[mn]->u_lambda),
                        &(elc_glob[mn]->len_u_lambda), model_name, SCALAR_INPUT, &NO_SPECIES, es);
  if (model_read == -1) {
    if (!strcmp(model_name, "POWER_LAW")) {
      elc_glob[mn]->lame_lambda_model = POWER_LAW;
      num_const = read_constants(imp, &(elc_glob[mn]->u_lambda), NO_SPECIES);
      if (num_const < 3) {
        sr = sprintf(err_msg, "Matl %s expected at least 3 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Lame LAMBDA", "POWER_LAW");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      elc_glob[mn]->len_u_lambda = num_const;
    } else if (!strcmp(model_name, "EXPONENTIAL")) {
      elc_glob[mn]->lame_lambda_model = EXPONENTIAL;
      num_const = read_constants(imp, &(elc_glob[mn]->u_lambda), NO_SPECIES);
      if (num_const < 3) {
        sr = sprintf(err_msg, "Matl %s expected at least 3 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Lame LAMBDA", "EXPONENTIAL");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      elc_glob[mn]->len_u_lambda = num_const;
    } else if (!strcmp(model_name, "POISSON_RATIO")) {
      elc_glob[mn]->lame_lambda_model = POISSON_RATIO;
      num_const = read_constants(imp, &(elc_glob[mn]->u_lambda), NO_SPECIES);
      if (num_const < 1) {
        sr = sprintf(err_msg, "Matl %s expected at least 1 constant for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Lame LAMBDA", "POISSON_RATIO");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      elc_glob[mn]->len_u_lambda = num_const;
    } else {
      sr = sprintf(err_msg, "Material %s - unrecognized model for %s \"%s\" ???\n",
                   pd_glob[mn]->MaterialName, "Lame LAMBDA", model_name);
      GOMA_EH(model_read, err_msg);
    }
    SPF_DBL_VEC(endofstring(es), num_const, elc_glob[mn]->u_lambda);
  }

  ECHO(es, echo_file);

  /* Temperature shift multiplier for Lame Lambda and Lame Mu */

  /*Initialize for good default behavior */
  elc_glob[mn]->lameTempShiftModel = CONSTANT;
  elc_glob[mn]->lame_TempShift = 1.;

  if (elc_glob[mn]->lame_mu_model == TABLE || elc_glob[mn]->lame_lambda_model == TABLE) {
    /* Only Constant or Table for now */
    model_read = look_for_mat_proptable(
        imp, "Lame Temperature Shift", &(elc_glob[mn]->lameTempShiftModel),
        &(elc_glob[mn]->lame_TempShift), &(elc_glob[mn]->u_lame_TempShift),
        &(elc_glob[mn]->len_u_lame_TempShift), &(elc_glob[mn]->lame_TempShift_tableid), model_name,
        SCALAR_INPUT, &NO_SPECIES, es);
    if (model_read == -1) {
      if (!strcmp(model_name, "POWER_LAW")) {
        elc_glob[mn]->lameTempShiftModel = POWER_LAW;
        num_const = read_constants(imp, &(elc_glob[mn]->u_lame_TempShift), NO_SPECIES);
        if (num_const < 3) {
          sr = sprintf(err_msg, "Matl %s expected at least 3 constants for %s %s model.\n",
                       pd_glob[mn]->MaterialName, "Lame Temperature Shift", "POWER_LAW");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
        elc_glob[mn]->len_u_lame_TempShift = num_const;
      } else {
        sr = sprintf(err_msg, "Material %s - unrecognized model for %s \"%s\" ???\n",
                     pd_glob[mn]->MaterialName, "LAME Temperature Shift", model_name);
        GOMA_EH(model_read, err_msg);
      }

      SPF_DBL_VEC(endofstring(es), num_const, elc_glob[mn]->u_lame_TempShift);
    }

    ECHO(es, echo_file);
  }

  model_read = look_for_mat_prop(
      imp, "Solid Viscosity", &(elc_glob[mn]->solid_viscosity_model),
      &(elc_glob[mn]->solid_viscosity), &(elc_glob[mn]->u_solid_viscosity),
      &(elc_glob[mn]->len_u_solid_viscosity), model_name, SCALAR_INPUT, &NO_SPECIES, es);

  if (model_read == -1) {
    elc_glob[mn]->solid_viscosity_model = CONSTANT;
    elc_glob[mn]->solid_viscosity = 0.0;

    SPF(es, "\t(%s = CONSTANT %.4g)", "Solid Viscosity", 0.0);
  }
  ECHO(es, echo_file);

  model_read = look_for_mat_prop(imp, "Solid Dilational Viscosity",
                                 &(elc_glob[mn]->solid_dil_viscosity_model),
                                 &(elc_glob[mn]->solid_dil_viscosity), NO_USER, NULL, model_name,
                                 SCALAR_INPUT, &NO_SPECIES, es);

  if (model_read == -1) {
    elc_glob[mn]->solid_dil_viscosity_model = CONSTANT;
    elc_glob[mn]->solid_dil_viscosity = 0.0;

    SPF(es, "\t(%s = CONSTANT %.4g)", "Solid Dilational Viscosity", 0.0);
  }
  ECHO(es, echo_file);

  /* Bending stiffness of structural shells */
  model_read = look_for_mat_prop(
      imp, "Shell bending stiffness", &(elc_glob[mn]->bend_stiffness_model),
      &(elc_glob[mn]->bend_stiffness), &(elc_glob[mn]->u_bend_stiffness),
      &(elc_glob[mn]->len_u_bend_stiffness), model_name, SCALAR_INPUT, &NO_SPECIES, es);

  /* Model must be CONSTANT for now! */
  if (model_read == 1 && strcmp(model_name, "CONSTANT") != 0) {
    sr = sprintf(err_msg, "Material %s - unrecognized model for %s \"%s\" ???\n",
                 pd_glob[mn]->MaterialName, "Shell bending stiffness", model_name);
    GOMA_EH(model_read, err_msg);
  } else if (model_read == -1) {
    /* Default to CONSTANT(1) */
    elc_glob[mn]->bend_stiffness_model = CONSTANT;
    elc_glob[mn]->bend_stiffness = 1.0;
  }

  ECHO(es, echo_file);

  /* Extensional stiffness of structural shells */
  model_read = look_for_mat_prop(
      imp, "Shell extensional stiffness", &(elc_glob[mn]->exten_stiffness_model),
      &(elc_glob[mn]->exten_stiffness), &(elc_glob[mn]->u_exten_stiffness),
      &(elc_glob[mn]->len_u_exten_stiffness), model_name, SCALAR_INPUT, &NO_SPECIES, es);

  /* Model must be CONSTANT for now! */
  if (model_read == 1 && strcmp(model_name, "CONSTANT") != 0) {
    sr = sprintf(err_msg, "Material %s - unrecognized model for %s \"%s\" ???\n",
                 pd_glob[mn]->MaterialName, "Extensional bending stiffness", model_name);
    GOMA_EH(model_read, err_msg);
  } else if (model_read == -1) {
    /* Default to CONSTANT(1) */
    elc_glob[mn]->exten_stiffness_model = CONSTANT;
    elc_glob[mn]->exten_stiffness = 1.0;
  }

  ECHO(es, echo_file);

  /* Poisson ratio of structural shells */
  model_read =
      look_for_mat_prop(imp, "Shell Poisson ratio", &(elc_glob[mn]->poisson_model),
                        &(elc_glob[mn]->poisson), &(elc_glob[mn]->u_poisson),
                        &(elc_glob[mn]->len_u_poisson), model_name, SCALAR_INPUT, &NO_SPECIES, es);

  /* Model must be CONSTANT for now! */
  if (model_read == 1 && strcmp(model_name, "CONSTANT") != 0) {
    sr = sprintf(err_msg, "Material %s - unrecognized model for %s \"%s\" ???\n",
                 pd_glob[mn]->MaterialName, "Shell Poisson ratio", model_name);
    GOMA_EH(model_read, err_msg);
  } else if (model_read == -1) {
    /* Default to CONSTANT(1) */
    elc_glob[mn]->poisson_model = CONSTANT;
    elc_glob[mn]->poisson = 0.5;
  }

  ECHO(es, echo_file);

  /* check for shell tangent calculator model */
  if (pd_glob[mn]->gv[R_MESH1] && pd_glob[mn]->gv[R_SHELL_CURVATURE]) {
    char input[MAX_CHAR_IN_INPUT] = "zilch\0";
    strcpy(search_string, "Shell Tangent Computation Method");
    model_read = look_for_optional(imp, search_string, input, '=');

    if (model_read == 1) {
      if (fscanf(imp, "%s", model_name) != 1) {
        sr = sprintf(err_msg, "Error reading model name string in material file, property %s",
                     search_string);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      SPF(es, "%s = %s", search_string, model_name);
      if (model_read == 1 && !strcmp(model_name, "ISOPARAMETRIC")) {
        model_read = 1;
        mat_ptr->shell_tangent_model = ISOPARAMETRIC;
        //	num_const = read_constants(imp, &(mat_ptr->shell_tangent_seed_vec),
        // NO_SPECIES);

        if (num_const != 0) {
          GOMA_EH(GOMA_ERROR,
                  "ISOPARAMETRIC Shell Tangent Computation Method takes no other input parameters");
        }
      } else if (model_read == 1 && !strcmp(model_name, "SEEDED")) {
        model_read = 1;
        mat_ptr->shell_tangent_model = SEEDED;
        num_const = read_constants(imp, &(mat_ptr->shell_tangent_seed_vec_const), NO_SPECIES);
        mat_ptr->len_shell_tangent_seed_vec_const = num_const;

        SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->shell_tangent_seed_vec_const);

        if (num_const != 3) {
          GOMA_EH(GOMA_ERROR, "SEEDED Shell Tangent Computation Model requires input of the seed "
                              "as three components of a unit vector.");
        }
      }

      else {
        // default is isoparametric
        mat_ptr->shell_tangent_model = ISOPARAMETRIC;
        SPF(es, "%s = %s", search_string, "ISOPARAMETRIC");
      }

      ECHO(es, echo_file);
    }
  }

  /* check for shell moment tensor calculator model */
  if (pd_glob[mn]->gv[R_MESH1] && pd_glob[mn]->gv[R_SHELL_CURVATURE]) {
    char input[MAX_CHAR_IN_INPUT] = "zilch\0";
    strcpy(search_string, "Shell Moment Tensor Model");
    model_read = look_for_optional(imp, search_string, input, '=');

    if (model_read == 1) {
      if (fscanf(imp, "%s", model_name) != 1) {
        sr = sprintf(err_msg, "Error reading model name string in material file, property %s",
                     search_string);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      SPF(es, "%s = %s", search_string, model_name);
      if (model_read == 1 && !strcmp(model_name, "EXPANDED")) {
        model_read = 1;
        mat_ptr->shell_moment_tensor_model = SMT_EXPANDED;
        //	num_const = read_constants(imp, &(mat_ptr->shell_tangent_seed_vec),
        // NO_SPECIES);

        if (num_const != 0) {
          GOMA_EH(GOMA_ERROR, "EXPANDED Shell Moment Tensor Model takes no other input parameters");
        }
      } else

          if (model_read == 1 && !strcmp(model_name, "SIMPLE")) {
        model_read = 1;
        mat_ptr->shell_moment_tensor_model = SMT_SIMPLE;
        //	num_const = read_constants(imp, &(mat_ptr->shell_tangent_seed_vec),
        // NO_SPECIES);

        if (num_const != 0) {
          GOMA_EH(GOMA_ERROR, "SIMPLE Shell Moment Tensor Model takes no other input parameters");
        }
      }

      else {
        // default is simple
        mat_ptr->shell_tangent_model = SMT_SIMPLE;
        SPF(es, "%s = %s", search_string, "SIMPLE");
      }

      ECHO(es, echo_file);
    }
  }

  model_read = look_for_mat_prop(imp, "Stress Free Solvent Vol Frac", &(LameLambdaModel),
                                 &(elc_glob[mn]->Strss_fr_sol_vol_frac), NO_USER, NULL, model_name,
                                 SCALAR_INPUT, &NO_SPECIES, es);
  ECHO(es, echo_file);

  model_read = look_for_mat_proptable(
      imp, "Solid Thermal Expansion", &(elc_glob[mn]->thermal_expansion_model),
      &(elc_glob[mn]->thermal_expansion), &(elc_glob[mn]->u_thermal_expansion),
      &(elc_glob[mn]->len_u_thermal_expansion), &(elc_glob[mn]->thermal_expansion_tableid),
      model_name, SCALAR_INPUT, &NO_SPECIES, es);

  if (model_read == -1) {

    if (!strcmp(model_name, "SHRINKAGE")) {
      elc_glob[mn]->thermal_expansion_model = SHRINKAGE;
      num_const = read_constants(imp, &(elc_glob[mn]->u_thermal_expansion), NO_SPECIES);
      if (num_const < 2) {
        sr = sprintf(err_msg, "Matl %s expected at least 2 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Thermal Expansion", "SHRINKAGE");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      elc_glob[mn]->len_u_thermal_expansion = num_const;
    } else if (!strcmp(model_name, "IDEAL_GAS")) {
      elc_glob[mn]->thermal_expansion_model = IDEAL_GAS;
      num_const = read_constants(imp, &(elc_glob[mn]->u_thermal_expansion), NO_SPECIES);
      if (num_const < 1) {
        sr = sprintf(err_msg, "Matl %s expected at least 1 constant for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Thermal Expansion", "IDEAL_GAS");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      elc_glob[mn]->len_u_thermal_expansion = num_const;
    } else if (!strcmp(model_name, "THERMAL")) {
      elc_glob[mn]->thermal_expansion_model = THERMAL_HEAT;
      num_const = read_constants(imp, &(elc_glob[mn]->u_thermal_expansion), NO_SPECIES);
      if (num_const < 5) {
        sprintf(err_msg, "Material %s - expected at least 5 constants for %s %s model.\n",
                pd_glob[mn]->MaterialName, search_string, "THERMAL");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      elc_glob[mn]->len_u_thermal_expansion = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, elc_glob[mn]->u_thermal_expansion);
    } else if (!strcmp(model_name, "ORTHOTROPIC")) {
      elc_glob[mn]->thermal_expansion_model = ORTHOTROPIC;
      num_const = read_constants(imp, &(elc_glob[mn]->u_thermal_expansion), NO_SPECIES);
      if (num_const < 6) {
        sprintf(err_msg, "Material %s - expected at least 6 constants for %s %s model.\n",
                pd_glob[mn]->MaterialName, search_string, "ORTHOTROPIC");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      elc_glob[mn]->len_u_thermal_expansion = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, elc_glob[mn]->u_thermal_expansion);
    }

    else {
      elc_glob[mn]->thermal_expansion_model = CONSTANT;
      elc_glob[mn]->thermal_expansion = 0.0;
    }

    SPF(es, "\t(%s = CONSTANT %.4g)", "Solid Thermal Expansion", 0.0);
  }

  ECHO(es, echo_file);

  model_read = look_for_mat_prop(imp, "Solid Reference Temperature",
                                 &(elc_glob[mn]->solid_reference_temp_model),
                                 &(elc_glob[mn]->solid_reference_temp), NO_USER, NULL, model_name,
                                 SCALAR_INPUT, &NO_SPECIES, es);

  if (model_read == -1) {
    elc_glob[mn]->solid_reference_temp_model = CONSTANT;
    elc_glob[mn]->solid_reference_temp = 0.0;

    SPF(es, "\t(%s = CONSTANT %.4g)", "Solid Reference Temperature", 0.0);
  }
  ECHO(es, echo_file);

  if (evpl_glob[mn]->ConstitutiveEquation == EVP_HYPER) {
    /*get plastic viscosity and yield stress */
    model_read = look_for_mat_prop(imp, "Plastic Viscosity", &(evpl_glob[mn]->plastic_mu_model),
                                   &(evpl_glob[mn]->plastic_mu), NO_USER, NULL, model_name,
                                   SCALAR_INPUT, &NO_SPECIES, es);

    if (model_read == -1) {
      if (!strcmp(model_name, "LINEAR")) {
        evpl_glob[mn]->plastic_mu_model = LINEAR;
        num_const = read_constants(imp, &(evpl_glob[mn]->u_plastic_mu), NO_SPECIES);
        if (num_const < 2) {
          sr = sprintf(err_msg, "Matl %s expected at least 2 constants for %s %s model.\n",
                       pd_glob[mn]->MaterialName, "Plastic Viscosity", "LINEAR");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
        evpl_glob[mn]->len_u_plastic_mu = num_const;

      } else {
        sr = sprintf(err_msg, "Material %s - unrecognized model for %s \"%s\" ???\n",
                     pd_glob[mn]->MaterialName, "Plastic Viscosity", model_name);
        GOMA_EH(model_read, err_msg);
      }
      SPF_DBL_VEC(endofstring(es), num_const, evpl_glob[mn]->u_plastic_mu);
    }

    ECHO(es, echo_file);

    model_read = look_for_mat_prop(imp, "EVP Yield Stress", &(evpl_glob[mn]->yield_model),
                                   &(evpl_glob[mn]->yield), NO_USER, NULL, model_name, SCALAR_INPUT,
                                   &NO_SPECIES, es);

    if (model_read == -1 && !strcmp(model_name, "LINEAR")) {
      evpl_glob[mn]->yield_model = LINEAR;
      num_const = read_constants(imp, &(evpl_glob[mn]->u_yield), NO_SPECIES);
      if (num_const < 2) {
        sr = sprintf(err_msg, "Matl %s expected at least 2 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "EVP Yield Stress", "LINEAR");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      evpl_glob[mn]->len_u_yield = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, evpl_glob[mn]->u_yield);
    } else {
      sr = sprintf(err_msg, "Material %s - unrecognized model for %s \"%s\" ???\n",
                   pd_glob[mn]->MaterialName, "EVP Yield Stress", model_name);
      GOMA_EH(model_read, err_msg);
    }
    ECHO(es, echo_file);
  }
  /*********************************************************************/
  /*Before we continue, test to see if this is a Total ALE material, in
   *which case the solid props just parsed will be loaded into the real solid
   *properties, and some additional properties of the pseudo-solid will now
   *be loaded up.  This next section is ignored if you are ARBITRARY or LAGRANGIAN
   *material type
   */
  /*********************************************************************/

  if (pd_glob[mn]->MeshMotion == TOTAL_ALE) {

    /* First transfer already parsed real-properties to the appropriate
     * structure elements. Do this by swapping pointers. You need to maintain
     * though some of the elastic (elc) and constitutive relations (cr)
     * constants in both.
     */

    cr_glob[mn]->RealSolidFluxModel = ElasticConstitutiveEquation;
    dum_ptr = elc_rs_glob[mn];
    elc_rs_glob[mn] = elc_glob[mn];
    elc_glob[mn] = dum_ptr;
    elc_glob[mn]->Strss_fr_sol_vol_frac = elc_rs_glob[mn]->Strss_fr_sol_vol_frac;
    elc_glob[mn]->thermal_expansion = 0.;
    elc_glob[mn]->thermal_expansion_model = CONSTANT;
    elc_glob[mn]->solid_reference_temp_model = CONSTANT;
    elc_glob[mn]->solid_reference_temp = 25.;
    elc_glob[mn]->solid_viscosity = 0.;
    elc_glob[mn]->solid_viscosity_model = CONSTANT;
    elc_glob[mn]->solid_dil_viscosity = 0.;
    elc_glob[mn]->solid_dil_viscosity_model = CONSTANT;

    model_read =
        look_for_mat_prop(imp, "Pseudo-Solid Constitutive Equation", &(ElasticConstitutiveEquation),
                          &a0, NO_USER, NULL, model_name, NO_INPUT, &NO_SPECIES, es);

    if (model_read == 1)
      GOMA_EH(GOMA_ERROR, "Can't have USER or CONSTANT for Pseudo-Solid Constitutive Equation");

    if (!strcmp(model_name, "LINEAR")) {
      ElasticConstitutiveEquation = LINEAR;
    } else if (!strcmp(model_name, "NONLINEAR") || !strcmp(model_name, "NONLINEAR_PLANE_STRAIN")) {
      ElasticConstitutiveEquation = NONLINEAR;
    } else if (!strcmp(model_name, "INCOMP_PSTRAIN")) {
      ElasticConstitutiveEquation = INCOMP_PSTRAIN;
    } else if (!strcmp(model_name, "INCOMP_PSTRESS")) {
      ElasticConstitutiveEquation = INCOMP_PSTRESS;
    } else if (!strcmp(model_name, "HOOKEAN_PSTRAIN")) {
      ElasticConstitutiveEquation = HOOKEAN_PSTRAIN;
    } else if (!strcmp(model_name, "HOOKEAN_PSTRESS")) {
      ElasticConstitutiveEquation = HOOKEAN_PSTRESS;
    } else if (!strcmp(model_name, "INCOMP_3D")) {
      ElasticConstitutiveEquation = INCOMP_3D;
    } else if (!strcmp(model_name, "KELVIN_VOIGT")) {
      ElasticConstitutiveEquation = KELVIN_VOIGT;
    } else if (!strcmp(model_name, "EVP_HYPER")) {
      ElasticConstitutiveEquation = EVP_HYPER;
      GOMA_EH(GOMA_ERROR, "Pseudo-solid constitutive equationis elasto-viscoplastic. Hmm. Won't "
                          "continue with my better judgement");
    } else /* default to nonlinear */
    {
      ElasticConstitutiveEquation = NONLINEAR;
      GOMA_WH(-1, "Pseudo-solid Elastic Constitutive set to default NONLINEAR");
    }

    ECHO(es, echo_file);

    /*This seems a little redundant! prs (2/24/95) */
    cr_glob[mn]->MeshFluxModel = ElasticConstitutiveEquation;

    /* An optional additional type of mesh motion which includes an intertial term in the
       momentum equation */

    /* read in constants for constitutive equation if they are input */

    model_read =
        look_for_mat_prop(imp, "Pseudo-Solid Lame MU", &(elc_glob[mn]->lame_mu_model),
                          &(elc_glob[mn]->lame_mu), &(elc_glob[mn]->u_mu),
                          &(elc_glob[mn]->len_u_mu), model_name, SCALAR_INPUT, &NO_SPECIES, es);
    if (model_read == -1) {
      if (!strcmp(model_name, "CONTACT_LINE")) {
        elc_glob[mn]->lame_mu_model = CONTACT_LINE;
        num_const = read_constants(imp, &(elc_glob[mn]->u_mu), NO_SPECIES);
        if (num_const < 4) {
          sr = sprintf(err_msg, "Matl %s expected at least 4 constants for %s %s model.\n",
                       pd_glob[mn]->MaterialName, "Lame MU", "CONTACT_LINE");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
        elc_glob[mn]->len_u_mu = num_const;
        SPF_DBL_VEC(endofstring(es), num_const, elc_glob[mn]->u_mu);
      } else {
        sr = sprintf(err_msg, "Material %s - unrecognized model for %s \"%s\" ???\n",
                     pd_glob[mn]->MaterialName, "Pseudo-Solid Lame MU", model_name);
        GOMA_EH(model_read, err_msg);
      }

      ECHO(es, echo_file);
    }

    model_read =
        look_for_mat_prop(imp, "Pseudo-Solid Lame LAMBDA", &(elc_glob[mn]->lame_lambda_model),
                          &(elc_glob[mn]->lame_lambda), &(elc_glob[mn]->u_lambda),
                          &(elc_glob[mn]->len_u_lambda), model_name, SCALAR_INPUT, &NO_SPECIES, es);
    if (model_read == -1) {
      if (!strcmp(model_name, "POISSON_RATIO")) {
        elc_glob[mn]->lame_lambda_model = POISSON_RATIO;
        num_const = read_constants(imp, &(elc_glob[mn]->u_lambda), NO_SPECIES);
        if (num_const < 1) {
          sr = sprintf(err_msg, "Matl %s expected at least 1 constant for %s %s model.\n",
                       pd_glob[mn]->MaterialName, "Lame LAMBDA", "POISSON_RATIO");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
        elc_glob[mn]->len_u_lambda = num_const;
        SPF_DBL_VEC(endofstring(es), num_const, elc_glob[mn]->u_mu);
      } else {
        sr = sprintf(err_msg, "Material %s - unrecognized model for %s \"%s\" ???\n",
                     pd_glob[mn]->MaterialName, "Lame LAMBDA", model_name);
        GOMA_EH(model_read, err_msg);
      }
    }

    ECHO(es, echo_file);
  } // End TOTAL ALE

  ECHO("\n---Fluid Constitutive Equation \n", echo_file);

  /*********************************************************************/
  /*
   * Read in generalized Newtonian constitutive equation type together
   * with the relevant constants depending on whether the the fluid is:
   *
   *       Newtonian: mu
   *       Power Law: mu0 nexp
   *       Carreau: mu0 nexp muinf lam aexp
   */
  model_read = look_for_mat_prop(imp, "Liquid Constitutive Equation", &(ConstitutiveEquation),
                                 &(a0), NO_USER, NULL, model_name, NO_INPUT, &NO_SPECIES, es);
  if (model_read == 1) {
    GOMA_EH(GOMA_ERROR,
            "Liquid Constitutive Equations options CONSTANT, USER, USER_GEN are not valid\n");
  }

  if (!strcmp(model_name, "NEWTONIAN")) {
    ConstitutiveEquation = NEWTONIAN;
  } else if (!strcmp(model_name, "POWER_LAW")) {
    ConstitutiveEquation = POWER_LAW;
  } else if (!strcmp(model_name, "POWERLAW_SUSPENSION")) {
    ConstitutiveEquation = POWERLAW_SUSPENSION;
  } else if (!strcmp(model_name, "CARREAU")) {
    ConstitutiveEquation = CARREAU;
  } else if (!strcmp(model_name, "CARREAU_SUSPENSION")) {
    ConstitutiveEquation = CARREAU_SUSPENSION;
  } else if (!strcmp(model_name, "SUSPENSION")) {
    ConstitutiveEquation = SUSPENSION;
  } else if (!strcmp(model_name, "EPOXY")) {
    ConstitutiveEquation = EPOXY;
  } else if (!strcmp(model_name, "SYLGARD")) {
    ConstitutiveEquation = SYLGARD;
  } else if (!strcmp(model_name, "FILLED_EPOXY")) {
    ConstitutiveEquation = FILLED_EPOXY;
  } else if (!strcmp(model_name, "FOAM_EPOXY")) {
    ConstitutiveEquation = FOAM_EPOXY;
  } else if (!strcmp(model_name, "THERMAL")) {
    ConstitutiveEquation = THERMAL;
  } else if (!strcmp(model_name, "CURE")) {
    ConstitutiveEquation = CURE;
  } else if (!strcmp(model_name, "BINGHAM")) {
    ConstitutiveEquation = BINGHAM;
  } else if (!strcmp(model_name, "BINGHAM_WLF")) {
    ConstitutiveEquation = BINGHAM_WLF;
  } else if (!strcmp(model_name, "BINGHAM_MIXED")) {
    ConstitutiveEquation = BINGHAM_MIXED;
  } else if (!strcmp(model_name, "CARREAU_WLF")) {
    ConstitutiveEquation = CARREAU_WLF;
  } else if (!strcmp(model_name, "HERSCHEL_BULKLEY")) {
    ConstitutiveEquation = HERSCHEL_BULKLEY;
  } else if (!strcmp(model_name, "BOND")) {
    ConstitutiveEquation = BOND;
  } else if (!strcmp(model_name, "BOND_SH")) {
    ConstitutiveEquation = BOND_SH;
  } else if (!strcmp(model_name, "CARREAU_WLF_CONC_PL")) {
    ConstitutiveEquation = CARREAU_WLF_CONC_PL;
  } else if (!strcmp(model_name, "CARREAU_WLF_CONC_EXP")) {
    ConstitutiveEquation = CARREAU_WLF_CONC_EXP;
  } else if (!strcmp(model_name, "FOAM_PMDI_10")) {
    ConstitutiveEquation = FOAM_PMDI_10;
  } else {
    GOMA_EH(GOMA_ERROR, "Unrecognizable Constitutive Equation");
  }

  ECHO(es, echo_file);

  pd_glob[mn]->MomentumFluxModel = ConstitutiveEquation;
  cr_glob[mn]->MomentumFluxModel = CR_MF_NEWTON_0;
  gn_glob[mn]->ConstitutiveEquation = ConstitutiveEquation;
  mp_glob[mn]->DilationalViscosityModel = DILVISCM_KAPPAWIPESMU;

  /* read in constants for constitutive equation if they are input */

  if (ConstitutiveEquation == NEWTONIAN) {
    model_read = look_for_mat_proptable(
        imp, "Viscosity", &(mp_glob[mn]->ViscosityModel), &(mp_glob[mn]->viscosity),
        &(mp_glob[mn]->u_viscosity), &(mp_glob[mn]->len_u_viscosity),
        &(mp_glob[mn]->viscosity_tableid), model_name, SCALAR_INPUT, &NO_SPECIES, es);

    if (model_read == -1 && !strcmp(model_name, "FILL")) {
      mat_ptr->ViscosityModel = FILL;

      num_const = read_constants(imp, &(mat_ptr->u_viscosity), 0);

      if (num_const < 2) {
        sr = sprintf(err_msg, "Matl %s expected at least 2 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Viscosity", "FILL");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_viscosity = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_viscosity);
    } else if (model_read == -1 && !strcmp(model_name, "LEVEL_SET")) {
      mat_ptr->ViscosityModel = LEVEL_SET;

      num_const = read_constants(imp, &(mat_ptr->u_viscosity), 0);

      if (num_const < 3) {
        sr = sprintf(err_msg, "Matl %s expected at least 3 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Viscosity", "LEVEL_SET");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      if (mat_ptr->u_viscosity[2] == 0.0)
        mat_ptr->u_viscosity[2] = ls->Length_Scale / 2.0;

      mat_ptr->len_u_viscosity = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_viscosity);
    } else if (model_read == -1 && !strcmp(model_name, "LEVEL_SET_Q")) {
      mat_ptr->ViscosityModel = LS_QUADRATIC;

      num_const = read_constants(imp, &(mat_ptr->u_viscosity), 0);

      if (num_const < 3) {
        sr = sprintf(err_msg, "Matl %s expected at least 3 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Viscosity", "LEVEL_SET");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      if (mat_ptr->u_viscosity[2] == 0.0)
        mat_ptr->u_viscosity[2] = ls->Length_Scale / 2.0;

      mat_ptr->len_u_viscosity = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_viscosity);
    } else if (model_read == -1 && !strcmp(model_name, "SUSPENSION_PM")) {
      mat_ptr->ViscosityModel = SUSPENSION_PM;

      num_const = read_constants(imp, &(mat_ptr->u_viscosity), 0);
      mat_ptr->viscosity = mat_ptr->u_viscosity[0];

      if (num_const < 1) {
        sr = sprintf(err_msg, "Matl %s expected at least 1 constant for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Viscosity", "SUSPENSION_PM");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_viscosity = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_viscosity);
    } else if (model_read == -1 && !strcmp(model_name, "CONST_PHASE_FUNCTION")) {
      mat_ptr->ViscosityModel = CONST_PHASE_FUNCTION;
      num_const = read_constants(imp, &(mat_ptr->u_viscosity), 0);

      if (num_const < pfd->num_phase_funcs + 2) {
        sr = sprintf(err_msg, "Matl %s expected at least %d constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, pfd->num_phase_funcs + 2, "Viscosity",
                     "CONST_PHASE_FUNCTION");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_viscosity = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_viscosity);
    } else {
      GOMA_EH(model_read, "Viscosity");
      gn_glob[mn]->mu0 = mat_ptr->viscosity;
    }
    gn_glob[mn]->mu0 = mat_ptr->viscosity;

    ECHO(es, echo_file);
  }

  if (ConstitutiveEquation == POWER_LAW || ConstitutiveEquation == POWERLAW_SUSPENSION ||
      ConstitutiveEquation == CARREAU || ConstitutiveEquation == CARREAU_SUSPENSION ||
      ConstitutiveEquation == BINGHAM || ConstitutiveEquation == BINGHAM_WLF ||
      ConstitutiveEquation == CARREAU_WLF || ConstitutiveEquation == SUSPENSION ||
      ConstitutiveEquation == EPOXY || ConstitutiveEquation == SYLGARD ||
      ConstitutiveEquation == FILLED_EPOXY || ConstitutiveEquation == THERMAL ||
      ConstitutiveEquation == CURE || ConstitutiveEquation == HERSCHEL_BULKLEY ||
      ConstitutiveEquation == CARREAU_WLF_CONC_PL || ConstitutiveEquation == CARREAU_WLF_CONC_EXP ||
      ConstitutiveEquation == BOND || ConstitutiveEquation == BOND_SH ||
      ConstitutiveEquation == FOAM_EPOXY || ConstitutiveEquation == FOAM_PMDI_10) {
    model_read =
        look_for_mat_prop(imp, "Low Rate Viscosity", &(gn_glob[mn]->mu0Model), &(gn_glob[mn]->mu0),
                          NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);

    if (model_read == -1 && !strcmp(model_name, "LEVEL_SET")) {
      gn_glob[mn]->mu0Model = LEVEL_SET;

      num_const = read_constants(imp, &(gn_glob[mn]->u_mu0), 0);

      if (num_const < 3) {
        sr = sprintf(err_msg, "Matl %s expected at least 3 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Low Rate Viscosity", "LEVEL_SET");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      gn_glob[mn]->len_u_mu0 = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, gn_glob[mn]->u_mu0);

      if (gn_glob[mn]->u_mu0[2] == 0.0)
        gn_glob[mn]->u_mu0[2] = ls->Length_Scale / 2.0;

    } else {
      GOMA_EH(model_read, "Low Rate Viscosity");
    }

    ECHO(es, echo_file);
  }

  if (ConstitutiveEquation == POWER_LAW || ConstitutiveEquation == POWERLAW_SUSPENSION ||
      ConstitutiveEquation == CARREAU || ConstitutiveEquation == CARREAU_SUSPENSION ||
      ConstitutiveEquation == BINGHAM || ConstitutiveEquation == BINGHAM_WLF ||
      ConstitutiveEquation == CARREAU_WLF || ConstitutiveEquation == SUSPENSION ||
      ConstitutiveEquation == FILLED_EPOXY || ConstitutiveEquation == HERSCHEL_BULKLEY ||
      ConstitutiveEquation == CARREAU_WLF_CONC_PL || ConstitutiveEquation == CARREAU_WLF_CONC_EXP) {
    model_read = look_for_mat_prop(imp, "Power Law Exponent", &(gn_glob[mn]->nexpModel),
                                   &(gn_glob[mn]->nexp), NO_USER, NULL, model_name, SCALAR_INPUT,
                                   &NO_SPECIES, es);

    if (model_read == -1 && !strcmp(model_name, "LEVEL_SET")) {
      gn_glob[mn]->nexpModel = LEVEL_SET;

      num_const = read_constants(imp, &(gn_glob[mn]->u_nexp), 0);

      if (num_const < 3) {
        sr = sprintf(err_msg, "Matl %s expected at least 3 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Power law exponent", "LEVEL_SET");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      gn_glob[mn]->len_u_nexp = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, gn_glob[mn]->u_nexp);

      if (gn_glob[mn]->u_nexp[2] == 0.0)
        gn_glob[mn]->u_nexp[2] = ls->Length_Scale / 2.0;

    } else {
      GOMA_EH(model_read, "Power Law Exponent");
    }
    ECHO(es, echo_file);
  }

  if (ConstitutiveEquation == CARREAU || ConstitutiveEquation == CARREAU_SUSPENSION ||
      ConstitutiveEquation == CARREAU_WLF || ConstitutiveEquation == BINGHAM ||
      ConstitutiveEquation == BINGHAM_WLF || ConstitutiveEquation == CARREAU_WLF_CONC_PL ||
      ConstitutiveEquation == CARREAU_WLF_CONC_EXP || ConstitutiveEquation == BOND_SH ||
      ConstitutiveEquation == BOND || ConstitutiveEquation == BINGHAM_MIXED) {
    model_read = look_for_mat_prop(imp, "High Rate Viscosity", &(gn_glob[mn]->muinfModel),
                                   &(gn_glob[mn]->muinf), NO_USER, NULL, model_name, SCALAR_INPUT,
                                   &NO_SPECIES, es);

    if (model_read == -1 && !strcmp(model_name, "LEVEL_SET")) {
      gn_glob[mn]->muinfModel = LEVEL_SET;

      num_const = read_constants(imp, &(gn_glob[mn]->u_muinf), 0);

      if (num_const < 3) {
        sr = sprintf(err_msg, "Matl %s expected at least 3 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Low Rate Viscosity", "LEVEL_SET");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      gn_glob[mn]->len_u_muinf = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, gn_glob[mn]->u_muinf);

      if (gn_glob[mn]->u_muinf[2] == 0.0)
        gn_glob[mn]->u_muinf[2] = ls->Length_Scale / 2.0;

    } else {
      GOMA_EH(model_read, "High Rate Viscosity");
    }
    ECHO(es, echo_file);

    if (ConstitutiveEquation != BOND && ConstitutiveEquation != BINGHAM_MIXED) {
      model_read =
          look_for_mat_prop(imp, "Time Constant", &(gn_glob[mn]->lamModel), &(gn_glob[mn]->lam),
                            NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);

      if (model_read == -1 && !strcmp(model_name, "LEVEL_SET")) {
        gn_glob[mn]->lamModel = LEVEL_SET;

        num_const = read_constants(imp, &(gn_glob[mn]->u_lam), 0);

        if (num_const < 3) {
          sr = sprintf(err_msg, "Matl %s expected at least 3 constants for %s %s model.\n",
                       pd_glob[mn]->MaterialName, "Low Rate Viscosity", "LEVEL_SET");
          GOMA_EH(GOMA_ERROR, err_msg);
        }

        gn_glob[mn]->len_u_lam = num_const;
        SPF_DBL_VEC(endofstring(es), num_const, gn_glob[mn]->u_lam);

        if (gn_glob[mn]->u_lam[2] == 0.0)
          gn_glob[mn]->u_lam[2] = ls->Length_Scale / 2.0;

      } else {
        GOMA_EH(model_read, "Time Constant");
      }

      ECHO(es, echo_file);
    }
  }

  if (ConstitutiveEquation == CARREAU || ConstitutiveEquation == CARREAU_SUSPENSION ||
      ConstitutiveEquation == CARREAU_WLF || ConstitutiveEquation == BINGHAM ||
      ConstitutiveEquation == BINGHAM_WLF || ConstitutiveEquation == CARREAU_WLF_CONC_PL ||
      ConstitutiveEquation == CARREAU_WLF_CONC_EXP || ConstitutiveEquation == BOND_SH ||
      ConstitutiveEquation == BOND) {
    model_read = look_for_mat_prop(imp, "Aexp", &(gn_glob[mn]->aexpModel), &(gn_glob[mn]->aexp),
                                   NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);

    if (model_read == -1 && !strcmp(model_name, "LEVEL_SET")) {
      gn_glob[mn]->aexpModel = LEVEL_SET;

      num_const = read_constants(imp, &(gn_glob[mn]->u_aexp), 0);

      if (num_const < 3) {
        sr = sprintf(err_msg, "Matl %s expected at least 3 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Low Rate Viscosity", "LEVEL_SET");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      gn_glob[mn]->len_u_aexp = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, gn_glob[mn]->u_aexp);

      if (gn_glob[mn]->u_aexp[2] == 0.0)
        gn_glob[mn]->u_aexp[2] = ls->Length_Scale / 2.0;

    } else {
      GOMA_EH(model_read, "Aexp");
    }

    ECHO(es, echo_file);
  }

  if (ConstitutiveEquation == BINGHAM || ConstitutiveEquation == BINGHAM_WLF ||
      ConstitutiveEquation == POWERLAW_SUSPENSION || ConstitutiveEquation == CARREAU_SUSPENSION ||
      ConstitutiveEquation == CARREAU_WLF || ConstitutiveEquation == EPOXY ||
      ConstitutiveEquation == SYLGARD || ConstitutiveEquation == FILLED_EPOXY ||
      ConstitutiveEquation == CARREAU_WLF_CONC_PL || ConstitutiveEquation == CARREAU_WLF_CONC_EXP ||
      ConstitutiveEquation == THERMAL || ConstitutiveEquation == BOND ||
      ConstitutiveEquation == FOAM_EPOXY || ConstitutiveEquation == FOAM_PMDI_10) {
    model_read = look_for_mat_prop(imp, "Thermal Exponent", &(gn_glob[mn]->atexpModel),
                                   &(gn_glob[mn]->atexp), NO_USER, NULL, model_name, SCALAR_INPUT,
                                   &NO_SPECIES, es);

    if (model_read == -1 && !strcmp(model_name, "LEVEL_SET")) {
      gn_glob[mn]->atexpModel = LEVEL_SET;

      num_const = read_constants(imp, &(gn_glob[mn]->u_atexp), 0);

      if (num_const < 3) {
        sr = sprintf(err_msg, "Matl %s expected at least 3 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Thermal Exponent", "LEVEL_SET");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      gn_glob[mn]->len_u_atexp = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, gn_glob[mn]->u_atexp);

      if (gn_glob[mn]->u_atexp[2] == 0.0)
        gn_glob[mn]->u_atexp[2] = ls->Length_Scale / 2.0;

    } else {
      GOMA_EH(model_read, "Thermal Exponent");
    }
    ECHO(es, echo_file);
  }

  if (ConstitutiveEquation == CARREAU_WLF || ConstitutiveEquation == BINGHAM_WLF ||
      ConstitutiveEquation == CARREAU_WLF_CONC_PL || ConstitutiveEquation == BOND ||
      ConstitutiveEquation == CARREAU_WLF_CONC_EXP) {
    model_read = look_for_mat_prop(imp, "Thermal WLF Constant2", &(gn_glob[mn]->wlfc2Model),
                                   &(gn_glob[mn]->wlfc2), NO_USER, NULL, model_name, SCALAR_INPUT,
                                   &NO_SPECIES, es);

    if (model_read == -1 && !strcmp(model_name, "LEVEL_SET")) {
      gn_glob[mn]->wlfc2Model = LEVEL_SET;

      num_const = read_constants(imp, &(gn_glob[mn]->u_wlfc2), 0);

      if (num_const < 3) {
        sr = sprintf(err_msg, "Matl %s expected at least 3 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Thermal Exponent", "LEVEL_SET");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      gn_glob[mn]->len_u_wlfc2 = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, gn_glob[mn]->u_wlfc2);

      if (gn_glob[mn]->u_wlfc2[2] == 0.0)
        gn_glob[mn]->u_wlfc2[2] = ls->Length_Scale / 2.0;

    } else {
      GOMA_EH(model_read, "Thermal WLF Constant2");
    }
    ECHO(es, echo_file);
  }

  if (ConstitutiveEquation == BINGHAM || ConstitutiveEquation == BINGHAM_WLF ||
      ConstitutiveEquation == BOND || ConstitutiveEquation == HERSCHEL_BULKLEY ||
      ConstitutiveEquation == BINGHAM_MIXED) {
    model_read =
        look_for_mat_prop(imp, "Yield Stress", &(gn_glob[mn]->tau_yModel), &(gn_glob[mn]->tau_y),
                          &(gn_glob[mn]->u_tau_y), &(gn_glob[mn]->len_u_tau_y), model_name,
                          SCALAR_INPUT, &NO_SPECIES, es);
    GOMA_EH(model_read, "Yield Stress");
    ECHO(es, echo_file);
  }

  if (ConstitutiveEquation == BINGHAM || ConstitutiveEquation == BOND ||
      ConstitutiveEquation == BINGHAM_WLF) {
    model_read =
        look_for_mat_prop(imp, "Yield Exponent", &(gn_glob[mn]->fexpModel), &(gn_glob[mn]->fexp),
                          NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);
    GOMA_EH(model_read, "Yield Exponent");
    ECHO(es, echo_file);
  }

  if (ConstitutiveEquation == BINGHAM_MIXED || ConstitutiveEquation == BINGHAM ||
      ConstitutiveEquation == BINGHAM_WLF) {
    model_read = look_for_mat_prop(imp, "Epsilon Regularization", &(gn_glob[mn]->epsilonModel),
                                   &(gn_glob[mn]->epsilon), NO_USER, NULL, model_name, SCALAR_INPUT,
                                   &NO_SPECIES, es);
    if ((ConstitutiveEquation == BINGHAM || ConstitutiveEquation == BINGHAM_WLF) &&
        model_read == -1) {
      gn_glob[mn]->epsilon = 0.0;
      gn_glob[mn]->epsilonModel = CONSTANT;
    } else {
      GOMA_EH(model_read, "Epsilon Regularization");
    }
    ECHO(es, echo_file);
  }
  /*
   * For now, apply thixotrophy to just the shear-thinning models although
   * it should be general for all
   */

  if (ConstitutiveEquation == POWER_LAW || ConstitutiveEquation == POWERLAW_SUSPENSION ||
      ConstitutiveEquation == CARREAU || ConstitutiveEquation == CARREAU_SUSPENSION ||
      ConstitutiveEquation == BINGHAM || ConstitutiveEquation == BINGHAM_WLF ||
      ConstitutiveEquation == CARREAU_WLF || ConstitutiveEquation == SUSPENSION ||
      ConstitutiveEquation == EPOXY || ConstitutiveEquation == SYLGARD ||
      ConstitutiveEquation == FILLED_EPOXY || ConstitutiveEquation == THERMAL ||
      ConstitutiveEquation == CURE || ConstitutiveEquation == HERSCHEL_BULKLEY ||
      ConstitutiveEquation == CARREAU_WLF_CONC_PL || ConstitutiveEquation == CARREAU_WLF_CONC_EXP ||
      ConstitutiveEquation == BOND || ConstitutiveEquation == BOND_SH ||
      ConstitutiveEquation == FOAM_EPOXY) {
    model_read = look_for_mat_prop(imp, "Thixotropic Factor", &(gn_glob[mn]->thixoModel),
                                   &(gn_glob[mn]->thixo_factor), NO_USER, NULL, model_name,
                                   SCALAR_INPUT, &NO_SPECIES, es);

    if (model_read == -1 && !strcmp(model_name, "LEVEL_SET")) {
      gn_glob[mn]->thixoModel = LEVEL_SET;

      num_const = read_constants(imp, &(gn_glob[mn]->u_thixo_factor), 0);

      if (num_const < 3) {
        sr = sprintf(err_msg, "Matl %s expected at least 3 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Thixotropic Factor", "LEVEL_SET");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      gn_glob[mn]->len_u_thixo = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, gn_glob[mn]->u_thixo_factor);

      if (gn_glob[mn]->u_thixo_factor[2] == 0.0)
        gn_glob[mn]->u_thixo_factor[2] = ls->Length_Scale / 2.0;

    } else {
      gn_glob[mn]->thixoModel = CONSTANT;
      fprintf(stderr, "MAT %d Thixotropic Factor = %g\n", mn, gn_glob[mn]->thixo_factor);
      GOMA_WH(model_read, "Defaulting Thixotropic Factor to Zero");
    }

    ECHO(es, echo_file);
  }

  if (ConstitutiveEquation == SUSPENSION || ConstitutiveEquation == POWERLAW_SUSPENSION ||
      ConstitutiveEquation == CARREAU_SUSPENSION || ConstitutiveEquation == FILLED_EPOXY) {
    model_read = look_for_mat_prop(imp, "Suspension Maximum Packing", &(gn_glob[mn]->maxpackModel),
                                   &(gn_glob[mn]->maxpack), NO_USER, NULL, model_name, SCALAR_INPUT,
                                   &NO_SPECIES, es);
    GOMA_EH(model_read, "Suspension Maximum Packing");
    ECHO(es, echo_file);

    iread = look_for_optional(imp, "Suspension Species Number", input, '=');
    if (fscanf(imp, "%d", &species_no) != 1) {
      GOMA_EH(GOMA_ERROR, "error reading Suspension Species Number");
    }
    gn_glob[mn]->sus_species_no = species_no;

    SPF(es, "%s %d", input, species_no);
    ECHO(es, echo_file);
  }

  if (ConstitutiveEquation == CARREAU_WLF_CONC_PL || ConstitutiveEquation == CARREAU_WLF_CONC_EXP) {
    model_read = look_for_mat_prop(imp, "Suspension Maximum Packing", &(gn_glob[mn]->maxpackModel),
                                   &(gn_glob[mn]->maxpack), NO_USER, NULL, model_name, SCALAR_INPUT,
                                   &NO_SPECIES, es);
    GOMA_EH(model_read, "Suspension Maximum Packing");
    ECHO(es, echo_file);

    model_read =
        look_for_mat_prop(imp, "Yield Exponent", &(gn_glob[mn]->fexpModel), &(gn_glob[mn]->fexp),
                          NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);
    GOMA_EH(model_read, "Yield Exponent");
    ECHO(es, echo_file);
  }

  if (ConstitutiveEquation == CURE || ConstitutiveEquation == EPOXY ||
      ConstitutiveEquation == FILLED_EPOXY || ConstitutiveEquation == FOAM_PMDI_10) {
    model_read = look_for_mat_prop(imp, "Cure Gel Point", &(gn_glob[mn]->gelpointModel),
                                   &(gn_glob[mn]->gelpoint), NO_USER, NULL, model_name,
                                   SCALAR_INPUT, &NO_SPECIES, es);
    GOMA_EH(model_read, "Cure Gel Point");
    ECHO(es, echo_file);
    model_read = look_for_mat_prop(imp, "Cure A Exponent", &(gn_glob[mn]->cureaexpModel),
                                   &(gn_glob[mn]->cureaexp), NO_USER, NULL, model_name,
                                   SCALAR_INPUT, &NO_SPECIES, es);
    GOMA_EH(model_read, "Cure  A Exponent");
    ECHO(es, echo_file);
    model_read = look_for_mat_prop(imp, "Cure B Exponent", &(gn_glob[mn]->curebexpModel),
                                   &(gn_glob[mn]->curebexp), NO_USER, NULL, model_name,
                                   SCALAR_INPUT, &NO_SPECIES, es);
    GOMA_EH(model_read, "Cure B Exponent");

    look_for(imp, "Cure Species Number", input, '=');
    if (fscanf(imp, "%d", &species_no) != 1) {
      GOMA_EH(GOMA_ERROR, "error reading Cure Species Number");
    }

    SPF(endofstring(es), " %d", species_no);
    ECHO(es, echo_file);

    gn_glob[mn]->cure_species_no = species_no;

    model_read = look_for_mat_prop(imp, "Unreacted Gel Temperature", &(gn_glob[mn]->tgel0Model),
                                   &(gn_glob[mn]->tgel0), NO_USER, NULL, model_name, SCALAR_INPUT,
                                   &NO_SPECIES, es);
    GOMA_EH(model_read, "Unreacted Gel Temperature");
    ECHO(es, echo_file);
  }

  if (ConstitutiveEquation == BOND_SH) {

    iread = look_for_optional(imp, "Suspension Species Number", input, '=');
    if (fscanf(imp, "%d", &species_no) != 1) {
      GOMA_EH(GOMA_ERROR, "error reading Suspension Species Number");
    }
    gn_glob[mn]->sus_species_no = species_no;

    SPF(endofstring(es), " %d", species_no);
    ECHO(es, echo_file);
  }
  if (ConstitutiveEquation == FOAM_EPOXY) {

    iread = look_for_optional(imp, "Suspension Species Number", input, '=');
    if (fscanf(imp, "%d", &species_no) != 1) {
      GOMA_EH(GOMA_ERROR, "error reading Suspension Species Number");
    }
    gn_glob[mn]->sus_species_no = species_no;

    SPF(endofstring(es), " %d", species_no);
    ECHO(es, echo_file);

    model_read = look_for_mat_prop(imp, "Cure Gel Point", &(gn_glob[mn]->gelpointModel),
                                   &(gn_glob[mn]->gelpoint), NO_USER, NULL, model_name,
                                   SCALAR_INPUT, &NO_SPECIES, es);
    GOMA_EH(model_read, "Cure Gel Point");
    ECHO(es, echo_file);

    look_for(imp, "Cure Species Number", input, '=');
    if (fscanf(imp, "%d", &species_no) != 1) {
      GOMA_EH(GOMA_ERROR, "error reading Cure Species Number");
    }
    SPF(es, "%s = %d", input, species_no);
    ECHO(es, echo_file);

    gn_glob[mn]->cure_species_no = species_no;
  }

  if (ConstitutiveEquation == SYLGARD) {
    model_read = look_for_mat_prop(imp, "Cure Gel Point", &(gn_glob[mn]->gelpointModel),
                                   &(gn_glob[mn]->gelpoint), NO_USER, NULL, model_name,
                                   SCALAR_INPUT, &NO_SPECIES, es);
    GOMA_EH(model_read, "Cure Gel Point");
    ECHO(es, echo_file);

    model_read = look_for_mat_prop(imp, "Cure A Exponent", &(gn_glob[mn]->cureaexpModel),
                                   &(gn_glob[mn]->cureaexp), NO_USER, NULL, model_name,
                                   SCALAR_INPUT, &NO_SPECIES, es);
    GOMA_EH(model_read, "Cure  A Exponent");
    ECHO(es, echo_file);

    look_for(imp, "Cure Species Number", input, '=');
    if (fscanf(imp, "%d", &species_no) != 1) {
      GOMA_EH(GOMA_ERROR, "error reading Cure Species Number");
    }
    SPF(es, "%s = %d", input, species_no);
    ECHO(es, echo_file);

    gn_glob[mn]->cure_species_no = species_no;
  }

  if (ConstitutiveEquation == BOND) {
    look_for(imp, "Bond Evolution Parameters", input, '=');
    if (fscanf(imp, "%le %le %le %le %le %le", &gn_glob[mn]->k1, &gn_glob[mn]->k2, &gn_glob[mn]->n0,
               &gn_glob[mn]->pexp, &gn_glob[mn]->qexp, &gn_glob[mn]->diff) != 6) {
      GOMA_EH(GOMA_ERROR, "error reading Bond Evolution Parameters");
    }
    SPF(es, "%s = ", "Bond Evolution Parameters");
    SPF(endofstring(es), " %.4g %.4g %.4g %.4g %.4g %.4g", gn_glob[mn]->k1, gn_glob[mn]->k2,
        gn_glob[mn]->n0, gn_glob[mn]->pexp, gn_glob[mn]->qexp, gn_glob[mn]->diff);
  }

  // Set the default
  mp_glob[mn]->DilationalViscosityModel = DILVISCM_KAPPAWIPESMU;
  model_read = look_for_mat_proptable(
      imp, "Dilational Viscosity", &(mp_glob[mn]->DilationalViscosityModel),
      &(mp_glob[mn]->dilationalViscosity), &(mp_glob[mn]->u_dilationalViscosity),
      &(mp_glob[mn]->len_u_dilationalViscosity), &(mp_glob[mn]->dilationalViscosity_tableid),
      model_name, SCALAR_INPUT, &NO_SPECIES, es);
  if (model_read == 1) {
    mp_glob[mn]->DilationalViscosityModel = DILVISCM_KAPPACONSTANT;
  }
  if (model_read == -1 && !strcmp(model_name, "DILVISCM_KAPPAWIPESMU")) {
    mp_glob[mn]->dilationalViscosity = 0.0;
  } else if (model_read == -1 && !strcmp(model_name, "DILVISCM_KAPPACONSTANT")) {
    mp_glob[mn]->DilationalViscosityModel = DILVISCM_KAPPACONSTANT;
    num_const = read_constants(imp, &(mat_ptr->u_dilationalViscosity), 0);
    mat_ptr->dilationalViscosity = mat_ptr->u_dilationalViscosity[0];

    if (num_const < 1) {
      sr = sprintf(err_msg, "Matl %s expected at least 1 constant for %s %s model.\n",
                   pd_glob[mn]->MaterialName, "Dilational Viscosity", "DILVISCM_KAPPACONSTANT");
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    mat_ptr->len_u_dilationalViscosity = num_const;
    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_dilationalViscosity);
  } else if (model_read == -1 && !strcmp(model_name, "DILVISCM_KAPPAFIXEDRATIO")) {
    mp_glob[mn]->DilationalViscosityModel = DILVISCM_KAPPAFIXEDRATIO;
    num_const = read_constants(imp, &(mat_ptr->u_dilationalViscosity), 0);
    mat_ptr->dilationalViscosityRatio = mat_ptr->u_dilationalViscosity[0];
    mat_ptr->dilationalViscosity = 0.0;
    if (num_const < 1) {
      sr = sprintf(err_msg, "Matl %s expected at least 1 constant for %s %s model.\n",
                   pd_glob[mn]->MaterialName, "Dilational Viscosity", "DILVISCM_KAPPAFIXEDRATIO");
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    mat_ptr->len_u_dilationalViscosity = num_const;
    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_dilationalViscosity);
  } else if (model_read == -1 && !strcmp(model_name, "DILVISCM_KAPPABUBBLES")) {
    mp_glob[mn]->DilationalViscosityModel = DILVISCM_KAPPABUBBLES;
    num_const = read_constants(imp, &(mat_ptr->u_dilationalViscosity), 0);
    mat_ptr->len_u_dilationalViscosity = num_const;
    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_dilationalViscosity);
  } else if (strcmp(model_name, " ")) {
    // We're here if we found the card, but couldn't read it
    GOMA_EH(model_read, "Dilational Viscosity");
  }

  model_read = look_for_mat_prop(imp, "Dilational Viscosity Multiplier", &(i0), &(a0), NO_USER,
                                 NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);

  mat_ptr->dilationalViscosityMultiplier = 1.0;

  if (model_read != -1) {
    mat_ptr->dilationalViscosityMultiplier = a0;

    stringup(model_name);

    if (strcmp(model_name, "CONSTANT")) {
      GOMA_EH(GOMA_ERROR, "Dilational Viscosity Multiplier can only be CONSTANT.\n");
    }
    ECHO(es, echo_file);
  }

  model_read = look_for_mat_prop(imp, "Second Level Set Viscosity", &(i0), &(a0), NO_USER, NULL,
                                 model_name, SCALAR_INPUT, &NO_SPECIES, es);

  if (model_read != -1) {

    if (ls == NULL)
      GOMA_EH(GOMA_ERROR,
              "Second Level Set Viscosity requires activation of Level Set Tracking.\n");

    mat_ptr->mp2nd->ViscosityModel = i0;
    mat_ptr->mp2nd->viscosity = a0;

    stringup(model_name);

    if (!strcmp(model_name, "CONSTANT") || !strcmp(model_name, "RATIO")) {
      if (fscanf(imp, "%s", input) != 1) {
        GOMA_EH(GOMA_ERROR, "Expecting trailing keyword for Second Level Set Viscosity.\n");
      }

      stringup(input);

      if (strncmp(input, "POSITIVE", 3) == 0) {
        mat_ptr->mp2nd->viscositymask[0] = 0;
        mat_ptr->mp2nd->viscositymask[1] = 1;
      } else if (strncmp(input, "NEGATIVE", 3) == 0) {
        mat_ptr->mp2nd->viscositymask[0] = 1;
        mat_ptr->mp2nd->viscositymask[1] = 0;
      } else {
        GOMA_EH(GOMA_ERROR,
                "Keyword must be POSITIVE or NEGATIVE for Second Level Set Viscosity.\n");
      }

      SPF(endofstring(es), " %s", input);
      if (pfd != NULL) {
        for (i = 0; i < pfd->num_phase_funcs; i++) {
          if (fscanf(imp, "%lf", &(mat_ptr->mp2nd->viscosity_phase[i])) != 1) {
            GOMA_EH(GOMA_ERROR, "error reading phase viscosity");
          }
          SPF(endofstring(es), " %g", mat_ptr->mp2nd->viscosity_phase[i]);
        }
      }
    } else {
      GOMA_EH(GOMA_ERROR, "Second Level Set Viscosity model can only be CONSTANT or RATIO.\n");
    }
    ECHO(es, echo_file);
  }
  /** Momentum Equation weight Function	**/

  strcpy(search_string, "Momentum Weight Function");
  model_read =
      look_for_mat_prop(imp, search_string, &(mat_ptr->Mwt_funcModel), &(mat_ptr->Mwt_func),
                        NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);
  if (strncmp(model_name, " ", 1) != 0) {
    if (!strcmp(model_name, "GALERKIN")) {
      mat_ptr->Mwt_funcModel = GALERKIN;
      mat_ptr->Mwt_func = 0.;
    } else if (!strcmp(model_name, "SUPG")) {
      int err;
      mat_ptr->Mwt_funcModel = SUPG;
      err = fscanf(imp, "%lg", &(mat_ptr->Mwt_func));
      if (err != 1) {
        GOMA_EH(GOMA_ERROR, "Expected to read one double for Momentum Weight Function SUPG");
      }
      SPF(endofstring(es), " %.4g", mat_ptr->Mwt_func);
    } else if (!strcmp(model_name, "SUPG_SHAKIB")) {
      int err;
      mat_ptr->Mwt_funcModel = SUPG_SHAKIB;
      err = fscanf(imp, "%lg", &(mat_ptr->Mwt_func));
      if (err != 1) {
        GOMA_EH(GOMA_ERROR, "Expected to read one double for Momentum Weight Function SUPG");
      }
      SPF(endofstring(es), " %.4g", mat_ptr->Mwt_func);
    } else {
      SPF(err_msg, "Syntax error or invalid model for %s\n", search_string);
      GOMA_EH(GOMA_ERROR, err_msg);
    }
  } else {
    mat_ptr->Mwt_funcModel = GALERKIN;
    mat_ptr->Mwt_func = 0.;
    SPF(es, "\t(%s = %s)", search_string, "GALERKIN");
  }

  ECHO(es, echo_file);

  /*
   *  Polymer Constitutive Equation
   *
   */
  model_read = look_for_mat_prop(imp, "Polymer Constitutive Equation", &(ConstitutiveEquation),
                                 &(a0), NO_USER, NULL, model_name, NO_INPUT, &NO_SPECIES, es);
  stringup(model_name);

  if (!strcmp(model_name, "GIESEKUS")) {
    vn_glob[mn]->ConstitutiveEquation = GIESEKUS;
  } else if (!strcmp(model_name, "WHITE_METZNER")) {
    vn_glob[mn]->ConstitutiveEquation = WHITE_METZNER;
  } else if (!strcmp(model_name, "OLDROYDB")) {
    vn_glob[mn]->ConstitutiveEquation = OLDROYDB;
  } else if (!strcmp(model_name, "PTT") || !strcmp(model_name, "PHAN THIEN-TANNER") ||
             !strcmp(model_name, "PHAN-THIEN TANNER")) {
    vn_glob[mn]->ConstitutiveEquation = PTT;
  } else if (!strcmp(model_name, "SARAMITO_OLDROYDB")) {
    vn_glob[mn]->ConstitutiveEquation = SARAMITO_OLDROYDB;
  } else if (!strcmp(model_name, "SARAMITO_GIESEKUS")) {
    vn_glob[mn]->ConstitutiveEquation = SARAMITO_GIESEKUS;
  } else if (!strcmp(model_name, "SARAMITO_PTT")) {
    vn_glob[mn]->ConstitutiveEquation = SARAMITO_PTT;
  } else if (!strcmp(model_name, "MODIFIED_JEFFREYS")) {
    vn_glob[mn]->ConstitutiveEquation = MODIFIED_JEFFREYS;
  } else if (!strcmp(model_name, "NOPOLYMER")) {
    vn_glob[mn]->ConstitutiveEquation = NOPOLYMER;
    /* set defaults if the next section is not entered */
    vn_glob[mn]->wt_funcModel = GALERKIN;
    vn_glob[mn]->wt_func = 0.;
    vn_glob[mn]->evssModel = EVSS_G;
    vn_glob[mn]->dg_J_model = FALSE;
  } else {
    vn_glob[mn]->ConstitutiveEquation = NOPOLYMER;
    /* set defaults if the next section is not entered */
    vn_glob[mn]->wt_funcModel = GALERKIN;
    vn_glob[mn]->wt_func = 0.;
    vn_glob[mn]->evssModel = EVSS_G;
    vn_glob[mn]->dg_J_model = FALSE;
    strcpy(es, "\t(Polymer Constitutive Equation = NOPOLYMER)");
  }
  ECHO(es, echo_file);

  /* this will be true if a VE constitutive equation is specified
   * if the NOPOLYMER option is used we can skip this section
   */
  if (vn_glob[mn]->ConstitutiveEquation) {

    if (vn_glob[mn]->ConstitutiveEquation == PTT) {

      strcpy(search_string, "PTT Form");

      model_read = look_for_mat_prop(imp, search_string, &(vn_glob[mn]->evssModel), &(a0), NO_USER,
                                     NULL, model_name, NO_INPUT, &NO_SPECIES, es);
      if (!strcmp(model_name, "LINEAR")) {
        vn_glob[mn]->ptt_type = PTT_LINEAR;
        SPF(es, "\t%s = %s", "PTT Form", model_name);
      } else if (!strcmp(model_name, "EXPONENTIAL")) {
        vn_glob[mn]->ptt_type = PTT_EXPONENTIAL;
        SPF(es, "\t%s = %s", "PTT Form", model_name);
      } else {
        GOMA_WH(GOMA_ERROR,
                "Unrecognized PTT Form: %s, expected LINEAR or EXPONENTIAL defaulting EXPONENTIAL",
                model_name);
        vn_glob[mn]->ptt_type = PTT_EXPONENTIAL;
        SPF(es, "\t(%s = %s)", "PTT Form", "EXPONENTIAL");
      }
      ECHO(es, echo_file);
    } else {
      vn_glob[mn]->ptt_type = PTT_EXPONENTIAL;
    }

    strcpy(search_string, "Polymer Stress Formulation");

    model_read = look_for_mat_prop(imp, search_string, &(vn_glob[mn]->evssModel), &(a0), NO_USER,
                                   NULL, model_name, NO_INPUT, &NO_SPECIES, es);
    if (!strcmp(model_name, "EVSS_G")) {
      if (vn_glob[mn]->ConstitutiveEquation == PTT ||
          vn_glob[mn]->ConstitutiveEquation == SARAMITO_PTT)
        GOMA_EH(GOMA_ERROR, "Error: EVSS_G stress formulation is not implemented in this case.");

      vn_glob[mn]->evssModel = EVSS_G;
    } else if (!strcmp(model_name, "EVSS_F")) {
      vn_glob[mn]->evssModel = EVSS_F;
    } else if (!strcmp(model_name, "EVSS_GRADV")) {
      vn_glob[mn]->evssModel = EVSS_GRADV;
    } else if (!strcmp(model_name, "EVSS_L")) {
      vn_glob[mn]->evssModel = EVSS_L;
    } else if (!strcmp(model_name, "LOG_CONF")) {
      vn_glob[mn]->evssModel = LOG_CONF;
    } else if (!strcmp(model_name, "SQRT_CONF")) {
      vn_glob[mn]->evssModel = SQRT_CONF;
    } else if (!strcmp(model_name, "CONF")) {
      vn_glob[mn]->evssModel = CONF;
    } else if (!strcmp(model_name, "LOG_CONF_TRANSIENT")) {
      vn_glob[mn]->evssModel = LOG_CONF_TRANSIENT;
    } else if (!strcmp(model_name, "LOG_CONF_TRANSIENT_GRADV")) {
      vn_glob[mn]->evssModel = LOG_CONF_TRANSIENT_GRADV;
    } else if (!strcmp(model_name, "LOG_CONF_GRADV")) {
      vn_glob[mn]->evssModel = LOG_CONF_GRADV;
    } else {
      if (vn_glob[mn]->ConstitutiveEquation == PTT ||
          vn_glob[mn]->ConstitutiveEquation == SARAMITO_PTT)
        GOMA_EH(GOMA_ERROR, "Error: EVSS_G stress formulation is not implemented in this case.");

      vn_glob[mn]->evssModel = EVSS_G; /* default to Rajagopalan's
                                          formulation */

      SPF(es, "\t(%s = %s)", "Polymer Stress Formulation", "EVSS_G");
    }

    ECHO(es, echo_file);

    strcpy(search_string, "Polymer Weight Function");

    model_read = look_for_mat_prop(imp, search_string, &(vn_glob[mn]->wt_funcModel), &(a0), NO_USER,
                                   NULL, model_name, NO_INPUT, &NO_SPECIES, es);

    if (!strcmp(model_name, "GALERKIN")) {
      vn_glob[mn]->wt_funcModel = GALERKIN;
    } else if (!strcmp(model_name, "SUPG")) {
      vn_glob[mn]->wt_funcModel = SUPG;
    } else {
      vn_glob[mn]->wt_funcModel = GALERKIN;

      SPF(es, "\t(%s = %s)", search_string, "GALERKIN");
    }

    ECHO(es, echo_file);

    if (vn_glob[mn]->wt_funcModel == SUPG) {

      strcpy(search_string, "Polymer Weighting");

      model_read =
          look_for_mat_prop(imp, search_string, &(ConstitutiveEquation), &(vn_glob[mn]->wt_func),
                            NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);

      GOMA_EH(model_read, "Polymer Weighting not set");
    } else {
      vn_glob[mn]->wt_func = 0.;
      SPF(es, "\t(%s = %s %g)", search_string, "CONSTANT", vn_glob[mn]->wt_func);
    }

    strcpy(search_string, "Polymer Shock Capturing");
    model_read =
        look_for_mat_prop(imp, search_string, &(mat_ptr->Ewt_funcModel), &(mat_ptr->Ewt_func),
                          NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);
    if (strncmp(model_name, " ", 1) != 0) {
      if (!strcmp(model_name, "NONE")) {
        vn_glob[mn]->shockcaptureModel = SC_NONE;
        vn_glob[mn]->shockcapture = 0.0;
      } else if (!strcmp(model_name, "DCDD")) {
        int err;
        vn_glob[mn]->shockcaptureModel = SC_DCDD;
        err = fscanf(imp, "%lg", &(vn_glob[mn]->shockcapture));
        if (err != 1) {
          GOMA_EH(GOMA_ERROR, "Expected to read one double for Polymer Shock Capturing = DCDD");
        }
        SPF(endofstring(es), " %.4g", vn_glob[mn]->shockcapture);
      } else if (!strcmp(model_name, "YZBETA")) {
        int err;
        vn_glob[mn]->shockcaptureModel = SC_YZBETA;
        err = fscanf(imp, "%lg", &(vn_glob[mn]->shockcapture));
        if (err != 1) {
          GOMA_EH(GOMA_ERROR, "Expected to read one double for Polymer Shock Capturing = YZBETA");
        }
        SPF(endofstring(es), " %.4g", vn_glob[mn]->shockcapture);
      } else {
        GOMA_EH(GOMA_ERROR, "Syntax error or invalid model for %s\n", search_string);
      }
    } else {
      vn_glob[mn]->shockcaptureModel = SC_NONE;
      vn_glob[mn]->shockcapture = 0.0;
      SPF(es, "\t(%s = %s)", search_string, "NONE");
    }

    ECHO(es, echo_file);

    strcpy(search_string, "Polymer Shift Function");

    if (look_forward_optional(imp, search_string, input, '=') == 1) {
      if (fscanf(imp, "%s", model_name) != 1) {
        GOMA_EH(GOMA_ERROR, "Need option for Polymer Shift Function ");
      }

      SPF(es, "%s = %s", search_string, model_name);

      if (!strcmp(model_name, "CONSTANT")) {
        vn_glob[mn]->shiftModel = CONSTANT;

        num_const = read_constants(imp, &(vn_glob[mn]->shift), NO_SPECIES);
        if (num_const < 1) {
          log_err("Matl %s expected at least 1 constants for %s model.\n",
                  pd_glob[mn]->MaterialName, "CONSTANT_WLF shift factor");
        }
        vn_glob[mn]->len_shift = num_const;

        SPF_DBL_VEC(endofstring(es), num_const, vn_glob[mn]->shift);

      } else if (!strcmp(model_name, "MODIFIED_WLF")) {
        vn_glob[mn]->shiftModel = MODIFIED_WLF;

        num_const = read_constants(imp, &(vn_glob[mn]->shift), NO_SPECIES);

        if (num_const < 2) {
          log_err("Matl %s expected at least 2 constants for %s model.\n",
                  pd_glob[mn]->MaterialName, "MODIFIED_WLF shift factor");
        }
        vn_glob[mn]->len_shift = num_const;

        SPF_DBL_VEC(endofstring(es), num_const, vn_glob[mn]->shift);
      } else {
        vn_glob[mn]->shiftModel = CONSTANT;
        vn_glob[mn]->shift = alloc_dbl_1(1, 1.0);
        vn_glob[mn]->len_shift = 1;
        SPF(es, "\t(%s = %s %.4g)", search_string, "CONSTANT", 1.0);
      }
    } else {
      vn_glob[mn]->shiftModel = CONSTANT;
      vn_glob[mn]->shift = alloc_dbl_1(1, 1.0);
      vn_glob[mn]->len_shift = 1;
      SPF(es, "\t(%s = %s %.4g)", search_string, "CONSTANT", 1.0);
    }

    ECHO(es, echo_file);

    strcpy(search_string, "Discontinuous Jacobian Formulation");

    if (look_forward_optional(imp, search_string, input, '=') == 1) {
      if (fscanf(imp, "%s", model_name) != 1) {
        GOMA_EH(GOMA_ERROR, "Need option for Discontinuous Jacobian Formulation ");
      }

      SPF(es, "%s = %s", search_string, model_name);

      if (!strcmp(model_name, "FULL")) {
        vn_glob[mn]->dg_J_model = FULL_DG;
      } else if (!strcmp(model_name, "EXPLICIT")) {
        vn_glob[mn]->dg_J_model = EXPLICIT_DG;

        num_const = read_constants(imp, &(vn_glob[mn]->dg_J_model_wt), NO_SPECIES);

        if (num_const < 1) {
          log_err("Matl %s expected at least 1 constants for %s model.\n",
                  pd_glob[mn]->MaterialName, "EXPLICIT_DG weighting factor");
        }
        vn_glob[mn]->len_dg_J_model_wt = num_const;

        SPF_DBL_VEC(endofstring(es), num_const, vn_glob[mn]->dg_J_model_wt);

      } else if (!strcmp(model_name, "SEGREGATED")) {
        vn_glob[mn]->dg_J_model = SEGREGATED;

        num_const = read_constants(imp, &(vn_glob[mn]->dg_J_model_wt), NO_SPECIES);
        if (num_const < 1) {
          log_err("Matl %s expected at least 1 constants for %s model.\n",
                  pd_glob[mn]->MaterialName, "SEGREGATED weighting factor");
        }
        vn_glob[mn]->len_dg_J_model_wt = num_const;
        SPF_DBL_VEC(endofstring(es), num_const, vn_glob[mn]->dg_J_model_wt);
      } else {
        vn_glob[mn]->dg_J_model = FALSE;
        SPF(es, "\t(%s = %s)", search_string, "FALSE");
      }
    } else {
      vn_glob[mn]->dg_J_model = FALSE;
      SPF(es, "\t(%s = %s)", search_string, "FALSE");
    }

    ECHO(es, echo_file);

    /* read in adaptive viscosity scaling: if it is zero
       we have the normal formulation without numerical artifacts */

    strcpy(search_string, "Adaptive Viscosity Scaling");

    model_read = look_for_mat_prop(imp, search_string, &(ConstitutiveEquation), &(vn_glob[mn]->eps),
                                   NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);
    if (model_read == -1) {
      vn_glob[mn]->eps = 0.;
      SPF(es, "\t(%s = %s %.4g)", search_string, "CONSTANT", vn_glob[mn]->eps);
    }

    ECHO(es, echo_file);

    /* allocate space */

    strcpy(search_string, "Polymer Viscosity");

    if (vn_glob[mn]->modes == 0)
      GOMA_EH(GOMA_ERROR, "Need to specify number of VE modes in input deck");

    modal_data = (dbl *)array_alloc(1, vn_glob[mn]->modes, sizeof(dbl));

    model_read =
        look_for_modal_prop(imp, search_string, vn_glob[mn]->modes, &matl_model, modal_data, es);
    if (model_read < 1) {
      if (model_read == -1)
        SPF(err_msg, "%s card is missing.", search_string);
      if (model_read == -2)
        SPF(err_msg, "Only CONSTANT, POWER LAW and HERSCHEL_BULKLEY %s mode model supported.",
            search_string);
      fprintf(stderr, "%s\n", err_msg);
      exit(-1);
    }

    // in case of non-constant polymer viscosity, parse polymer viscosity parameters
    // For now, these all assume a single node
    const bool mupIsConstant = matl_model == CONSTANT;

    int nExpModel = CONSTANT;
    int aExpModel = CONSTANT;
    int fExpModel = CONSTANT;
    int mu0Model = CONSTANT;
    int muInfModel = CONSTANT;
    int lamModel = CONSTANT;
    int tauyModel = CONSTANT;
    dbl nExpVal = 0;
    dbl aExpVal = 0;
    dbl fExpVal = 0;
    dbl mu0Val = 0;
    dbl muInfVal = 0;
    dbl lamVal = 0;
    dbl tauyVal = 0;

    if (!mupIsConstant) {
      model_read = look_for_mat_prop(imp, "Polymer Low Rate Viscosity", &(mu0Model), &(mu0Val),
                                     NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);
      printf("Polymer Low Rate Viscosity model %s\n", model_name);
      printf("Polymer Low Rate Viscosity value %E\n", mu0Val);
      ECHO(es, echo_file);

      model_read = look_for_mat_prop(imp, "Polymer Power Law Exponent", &(nExpModel), &(nExpVal),
                                     NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);
      printf("Polymer Power Law Exponent model %s\n", model_name);
      printf("Polymer Power Law Exponent value %E\n", nExpVal);
      ECHO(es, echo_file);

      model_read = look_for_mat_prop(imp, "Polymer High Rate Viscosity", &(muInfModel), &(muInfVal),
                                     NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);
      printf("Polymer High Rate Viscosity model %s\n", model_name);
      printf("Polymer High Rate Viscosity value %E\n", muInfVal);
      ECHO(es, echo_file);

      model_read = look_for_mat_prop(imp, "Polymer Viscosity Time Constant", &(lamModel), &(lamVal),
                                     NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);
      printf("Polymer Viscosity Time Constant model %s\n", model_name);
      printf("Polymer Viscosity Time Constant value %E\n", lamVal);
      ECHO(es, echo_file);

      model_read = look_for_mat_prop(imp, "Polymer Aexp", &(aExpModel), &(aExpVal), NO_USER, NULL,
                                     model_name, SCALAR_INPUT, &NO_SPECIES, es);
      printf("Polymer Aexp model %s\n", model_name);
      printf("Polymer Aexp value %E\n", aExpVal);
      ECHO(es, echo_file);

      model_read = look_for_mat_prop(imp, "Polymer Yield Stress", &(tauyModel), &(tauyVal), NO_USER,
                                     NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);
      printf("Polymer Yield Stress model %s\n", model_name);
      printf("Polymer Yield Stress value %E\n", tauyVal);
      ECHO(es, echo_file);

      model_read =
          look_for_mat_prop(imp, "Polymer Viscosity Yield Exponent", &(fExpModel), &(fExpVal),
                            NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);
      printf("Polymer Viscosity Yield Exponent model %s\n", model_name);
      printf("Polymer Viscosity Yield Exponent value %E\n", fExpVal);
      ECHO(es, echo_file);
    }

    ECHO(es, echo_file);

    for (mm = 0; mm < vn_glob[mn]->modes; mm++) {
      ve_glob[mn][mm]->gn->ConstitutiveEquation = matl_model;
      ve_glob[mn][mm]->gn->mu0 = (mupIsConstant ? modal_data[mm] : mu0Val);
      ve_glob[mn][mm]->gn->muinf = muInfVal;
      ve_glob[mn][mm]->gn->muinfModel = muInfModel;
      ve_glob[mn][mm]->gn->lam = lamVal;
      ve_glob[mn][mm]->gn->lamModel = lamModel;
      ve_glob[mn][mm]->gn->aexp = aExpVal;
      ve_glob[mn][mm]->gn->aexpModel = aExpModel;
      ve_glob[mn][mm]->gn->nexp = nExpVal;
      ve_glob[mn][mm]->gn->nexpModel = nExpModel;
      ve_glob[mn][mm]->gn->tau_yModel = tauyModel;
      ve_glob[mn][mm]->gn->tau_y = tauyVal;
      ve_glob[mn][mm]->gn->fexpModel = fExpModel;
      ve_glob[mn][mm]->gn->fexp = fExpVal;
    }

    strcpy(search_string, "Positive Level Set Polymer Viscosity");

    model_read =
        look_for_modal_prop(imp, search_string, vn_glob[mn]->modes, &matl_model, modal_data, es);

    if (model_read == 1) {

      if (ls == NULL)
        GOMA_EH(
            GOMA_ERROR,
            "Positive Level Set Polymer Viscosity requires activation of Level Set Tracking.\n");

      for (mm = 0; mm < vn_glob[mn]->modes; mm++) {
        ve_glob[mn][mm]->gn->pos_ls_mup = modal_data[mm];
        ve_glob[mn][mm]->gn->mu0Model = VE_LEVEL_SET;
        ve_glob[mn][mm]->gn->ConstitutiveEquation = VE_LEVEL_SET;
      }

      ECHO(es, echo_file);
    } else if (model_read == -2) {
      SPF(err_msg, "Only CONSTANT %s mode model supported.", search_string);
      fprintf(stderr, "%s\n", err_msg);
      exit(-1);
    }

    strcpy(search_string, "Polymer Time Constant");

    model_read =
        look_for_modal_prop(imp, search_string, vn_glob[mn]->modes, &matl_model, modal_data, es);
    if (model_read < 1) {
      if (model_read == -1)
        SPF(err_msg, "%s card is missing.", search_string);
      if (model_read == -2)
        SPF(err_msg, "Only CONSTANT %s mode model supported.", search_string);
      fprintf(stderr, "%s\n", err_msg);
      exit(-1);
    }

    ECHO(es, echo_file);

    for (mm = 0; mm < vn_glob[mn]->modes; mm++) {
      ve_glob[mn][mm]->time_const = modal_data[mm];
      ve_glob[mn][mm]->time_constModel = matl_model;
    }

    strcpy(search_string, "Positive Level Set Polymer Time Constant");

    model_read =
        look_for_modal_prop(imp, search_string, vn_glob[mn]->modes, &matl_model, modal_data, es);

    if (model_read == 1) {

      if (ls == NULL)
        GOMA_EH(GOMA_ERROR, "Positive Level Set Polymer Time Constant requires activation of Level "
                            "Set Tracking.\n");

      for (mm = 0; mm < vn_glob[mn]->modes; mm++) {
        ve_glob[mn][mm]->pos_ls.time_const = modal_data[mm];
        if (ve_glob[mn][mm]->time_constModel != CONSTANT) {
          fprintf(stderr, "%s\n",
                  "Only CONSTANT Polymer Time Constant model supported for viscoelastic level set");
          exit(-1);
        }
        ve_glob[mn][mm]->time_constModel = VE_LEVEL_SET;
      }

      ECHO(es, echo_file);
    } else if (model_read == -2) {
      SPF(err_msg, "Only CONSTANT %s mode model supported.", search_string);
      fprintf(stderr, "%s\n", err_msg);
      exit(-1);
    }

    if (vn_glob[mn]->ConstitutiveEquation == MODIFIED_JEFFREYS) {
      strcpy(search_string, "Jeffreys Viscosity");

      model_read =
          look_for_modal_prop(imp, search_string, vn_glob[mn]->modes, &matl_model, modal_data, es);

      if (model_read < 1) {
        if (model_read == -1)
          SPF(err_msg, "%s card is missing", search_string);
        if (model_read == -2)
          SPF(err_msg, "Only CONSTANT %s  mode model supported.", search_string);
        fprintf(stderr, "%s\n", err_msg);
        exit(-1);
      }

      ECHO(es, echo_file);

      for (mm = 0; mm < vn_glob[mn]->modes; mm++) {
        ve_glob[mn][mm]->muJeffreys = modal_data[mm];
        ve_glob[mn][mm]->muJeffreysModel = matl_model;
      }
    }

    if (vn_glob[mn]->ConstitutiveEquation == GIESEKUS ||
        vn_glob[mn]->ConstitutiveEquation == SARAMITO_GIESEKUS) {
      strcpy(search_string, "Mobility Parameter");

      model_read = look_for_modal_prop(imp, "Mobility Parameter", vn_glob[mn]->modes, &matl_model,
                                       modal_data, es);

      if (model_read < 1) {
        if (model_read == -1)
          SPF(err_msg, "%s is missing", search_string);
        if (model_read == -2)
          SPF(err_msg, "Only CONSTANT %s mode models supported.", search_string);
        fprintf(stderr, "%s\n", err_msg);
        exit(-1);
      }

      for (mm = 0; mm < vn_glob[mn]->modes; mm++) {
        ve_glob[mn][mm]->alpha = modal_data[mm];
        ve_glob[mn][mm]->alphaModel = matl_model;
      }

      ECHO(es, echo_file);

      strcpy(search_string, "Positive Level Set Mobility Parameter");

      model_read =
          look_for_modal_prop(imp, search_string, vn_glob[mn]->modes, &matl_model, modal_data, es);

      if (model_read == 1) {

        if (ls == NULL)
          GOMA_EH(
              GOMA_ERROR,
              "Positive Level Set Mobility Parameter requires activation of Level Set Tracking.\n");

        for (mm = 0; mm < vn_glob[mn]->modes; mm++) {
          ve_glob[mn][mm]->pos_ls.alpha = modal_data[mm];
          ve_glob[mn][mm]->alphaModel = VE_LEVEL_SET;
        }

        ECHO(es, echo_file);
      } else if (model_read == -2) {
        SPF(err_msg, "Only CONSTANT %s mode model supported.", search_string);
        fprintf(stderr, "%s\n", err_msg);
        exit(-1);
      }

    } else {
      for (mm = 0; mm < vn_glob[mn]->modes; mm++) {
        ve_glob[mn][mm]->alpha = 0.;
        ve_glob[mn][mm]->pos_ls.alpha = 0.;
        ve_glob[mn][mm]->alphaModel = CONSTANT;
      }
    }

    /*
     * If one of the Saramito model combinations is enabled, ensure that a yield stress card is
     * present
     */
    if (vn_glob[mn]->ConstitutiveEquation == SARAMITO_OLDROYDB ||
        vn_glob[mn]->ConstitutiveEquation == SARAMITO_PTT ||
        vn_glob[mn]->ConstitutiveEquation == SARAMITO_GIESEKUS) {
      /* Should yield stress be a modal property? Let's assume not for now */
      dbl tau_y_val;
      dbl fexp_val;
      dbl nexp_val;

      strcpy(search_string, "Polymer Yield Stress");
      model_read = look_for_mat_prop(imp, search_string, &(ConstitutiveEquation), &tau_y_val,
                                     NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);

      if (model_read < 1) {
        if (model_read == -1)
          SPF(err_msg, "%s card is missing.", search_string);
        if (model_read == -2)
          SPF(err_msg, "Only CONSTANT %s mode model supported.", search_string);
        fprintf(stderr, "%s\n", err_msg);
        exit(-1);
      }

      strcpy(search_string, "Yield Exponent");
      model_read = look_for_mat_prop(imp, search_string, &(ConstitutiveEquation), &fexp_val,
                                     NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);

      if (model_read < 1) {
        if (model_read == -1)
          SPF(err_msg, "%s card is missing.", search_string);
        if (model_read == -2)
          SPF(err_msg, "Only CONSTANT %s mode model supported.", search_string);
        fprintf(stderr, "%s\n", err_msg);
        exit(-1);
      }
      strcpy(search_string, "Saramito Power Law Exponent");
      model_read = look_for_mat_prop(imp, search_string, &(ConstitutiveEquation), &nexp_val,
                                     NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);

      if (model_read < 1) {
        if (model_read == -1) {
          nexp_val = 1;
        } else if (model_read == -2) {
          SPF(err_msg, "Only CONSTANT %s mode model supported.", search_string);
          GOMA_EH(GOMA_ERROR, err_msg);
        }
      }

      for (mm = 0; mm < vn_glob[mn]->modes; mm++) {
        ve_glob[mn][mm]->gn->tau_y = tau_y_val;
        ve_glob[mn][mm]->gn->fexp = fexp_val;
        ve_glob[mn][mm]->gn->nexp = nexp_val;

        // set polymer viscosity to (consistency index)^(1/nexp) when nexp != 1
        if (nexp_val != 1)
          ve_glob[mn][mm]->gn->mu0 = pow(ve_glob[mn][mm]->gn->mu0, 1. / nexp_val);
      }
      ECHO(es, echo_file);
    }

    if (vn_glob[mn]->ConstitutiveEquation == PTT ||
        vn_glob[mn]->ConstitutiveEquation == SARAMITO_PTT) {
      strcpy(search_string, "PTT Xi parameter");

      model_read =
          look_for_modal_prop(imp, search_string, vn_glob[mn]->modes, &matl_model, modal_data, es);

      if (model_read < 1) {
        if (model_read == -1)
          SPF(err_msg, "%s card is missing", search_string);
        if (model_read == -2)
          SPF(err_msg, "Only CONSTANT %s  mode model supported.", search_string);
        fprintf(stderr, "%s\n", err_msg);
        exit(-1);
      }

      for (mm = 0; mm < vn_glob[mn]->modes; mm++) {
        ve_glob[mn][mm]->xi = modal_data[mm];
        ve_glob[mn][mm]->xiModel = matl_model;
      }

      ECHO(es, echo_file);

      strcpy(search_string, "Positive Level Set PTT Xi parameter");

      model_read =
          look_for_modal_prop(imp, search_string, vn_glob[mn]->modes, &matl_model, modal_data, es);

      if (model_read == 1) {

        if (ls == NULL)
          GOMA_EH(
              GOMA_ERROR,
              "Positive Level Set PTT Xi parameter requires activation of Level Set Tracking.\n");

        for (mm = 0; mm < vn_glob[mn]->modes; mm++) {
          ve_glob[mn][mm]->pos_ls.xi = modal_data[mm];
          ve_glob[mn][mm]->xiModel = VE_LEVEL_SET;
        }

        ECHO(es, echo_file);
      } else if (model_read == -2) {
        SPF(err_msg, "Only CONSTANT %s mode model supported.", search_string);
        fprintf(stderr, "%s\n", err_msg);
        exit(-1);
      }

      strcpy(search_string, "PTT Epsilon parameter");

      model_read =
          look_for_modal_prop(imp, search_string, vn_glob[mn]->modes, &matl_model, modal_data, es);

      if (model_read < 1) {
        if (model_read == -1)
          SPF(err_msg, "%s card is missing", search_string);
        if (model_read == -2)
          SPF(err_msg, "Only CONSTANT %s  mode model supported.", search_string);
        fprintf(stderr, "%s\n", err_msg);
        exit(-1);
      }

      ECHO(es, echo_file);

      for (mm = 0; mm < vn_glob[mn]->modes; mm++) {
        ve_glob[mn][mm]->eps = modal_data[mm];
        ve_glob[mn][mm]->epsModel = matl_model;
      }

      strcpy(search_string, "Positive Level Set PTT Epsilon parameter");

      model_read =
          look_for_modal_prop(imp, search_string, vn_glob[mn]->modes, &matl_model, modal_data, es);

      if (model_read == 1) {

        if (ls == NULL)
          GOMA_EH(GOMA_ERROR, "Positive Level Set PTT Epsilon parameter requires activation of "
                              "Level Set Tracking.\n");

        for (mm = 0; mm < vn_glob[mn]->modes; mm++) {
          ve_glob[mn][mm]->pos_ls.eps = modal_data[mm];
          ve_glob[mn][mm]->epsModel = VE_LEVEL_SET;
        }

        ECHO(es, echo_file);
      } else if (model_read == -2) {
        SPF(err_msg, "Only CONSTANT %s mode model supported.", search_string);
        fprintf(stderr, "%s\n", err_msg);
        exit(-1);
      }
    } else {
      for (mm = 0; mm < vn_glob[mn]->modes; mm++) {
        ve_glob[mn][mm]->xi = 0.;
        ve_glob[mn][mm]->pos_ls.xi = 0.;
        ve_glob[mn][mm]->xiModel = CONSTANT;
      }

      for (mm = 0; mm < vn_glob[mn]->modes; mm++) {
        ve_glob[mn][mm]->eps = 0.;
        ve_glob[mn][mm]->pos_ls.eps = 0.;
        ve_glob[mn][mm]->epsModel = CONSTANT;
      }
    }

    free(modal_data);
  }

  /* surface tension */
  strcpy(search_string, "Surface Tension");

  model_read = look_for_mat_prop(imp, search_string, &(mat_ptr->SurfaceTensionModel),
                                 &(mat_ptr->surface_tension), &(mat_ptr->u_surface_tension),
                                 &(mat_ptr->len_u_surface_tension), model_name, SCALAR_INPUT,
                                 &NO_SPECIES, es);

  if (model_read == -1) {
    if (!strcmp(model_name, "DILATION")) {
      mat_ptr->SurfaceTensionModel = DILATION;

      num_const = read_constants(imp, &(mat_ptr->u_surface_tension), NO_SPECIES);

      if (num_const < 2) {
        sr = snprintf(err_msg, MAX_CHAR_ERR_MSG,
                      "Matl %s expected at least 2 constants for %s %s model.\n",
                      pd_glob[mn]->MaterialName, search_string, model_name);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_surface_tension = num_const;

      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_surface_tension);

    } else if (!strcmp(model_name, "GIBBS_ISOTHERM")) {
      mat_ptr->SurfaceTensionModel = GIBBS_ISOTHERM;

      num_const = read_constants(imp, &(mat_ptr->u_surface_tension), NO_SPECIES);

      if (num_const < 3) {
        sr = snprintf(err_msg, MAX_CHAR_ERR_MSG,
                      "Matl %s expected at least 3 constants for %s %s model.\n",
                      pd_glob[mn]->MaterialName, search_string, model_name);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_surface_tension = num_const;

      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_surface_tension);

    } else {
      mat_ptr->SurfaceTensionModel = CONSTANT;
      mat_ptr->surface_tension =
          1.; /* NOTE TO FUTURE DEVELOPERS:  FOR THE LOVE OF ALL THAT IS HOLY AND RIGHTGEOUS ON
               * GOD'S GREEN EARTH, PLEASE MAKE CERTAIN THAT THE DEFAULT VALUE FOR SURFACE TENSION
               * SET IN THIS ROUTINE IS 1.0 THE FUTURE OF ALL FREEDOM-LOVING PEOPLE DEPENDS UPON
               * THIS
               */

      SPF(es, "\t(%s = %s %.4g)", search_string, "CONSTANT", mat_ptr->surface_tension);
    }
  }

  ECHO(es, echo_file);

  /* Projection Equation Surface Diffusivity */

  strcpy(search_string, "Projection Equation Surface Diffusion Coefficient");

  model_read =
      look_for_mat_prop(imp, search_string, &i, &(mat_ptr->SurfaceDiffusionCoeffProjectionEqn),
                        NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);
  if (model_read == -1) {
    mat_ptr->SurfaceDiffusionCoeffProjectionEqn = 0.0;
  } else {
    if (strcmp(model_name, "CONSTANT")) {
      sr = sprintf(err_msg, "Matl %s: Surface diffusion model must be const.\n",
                   pd_glob[mn]->MaterialName);
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    num_const = 1;
    SPF(es, "\t(%s = %s %.4g)", search_string, "CONSTANT",
        mat_ptr->SurfaceDiffusionCoeffProjectionEqn);
  }

  /*
   * Thermal ProBperties
   *
   */

  ECHO("\n----Thermal Properties\n", echo_file);

  if (mn == 0)
    mat_ptr->thermal_cond_external_field = -1;

  strcpy(search_string, "Heat Flux Model");

  model_read = look_forward_optional(imp, search_string, input, '=');

  if (model_read == 1) {
    if (fscanf(imp, "%s", model_name) != 1) {
      GOMA_EH(GOMA_ERROR, "Need option for Heat Flux Model ");
    }

    SPF(es, "%s = %s", search_string, model_name);

    if (!strcmp(model_name, "USER")) {
      cr_glob[mn]->HeatFluxModel = CR_HF_USER;
    } else {
      cr_glob[mn]->HeatFluxModel = CR_HF_FOURIER_0;
    }
    ECHO(es, echo_file);
  }

  strcpy(search_string, "Conductivity");

  model_read = look_for_mat_proptable(
      imp, search_string, &(mat_ptr->ConductivityModel), &(mat_ptr->thermal_conductivity),
      &(mat_ptr->u_thermal_conductivity), &(mat_ptr->len_u_thermal_conductivity),
      &(mat_ptr->thermal_conductivity_tableid), model_name, SCALAR_INPUT, &NO_SPECIES, es);

  if (model_read == -1) {

    if (!strcmp(model_name, "CONST_LS") || !strcmp(model_name, "LEVEL_SET")) {
      mat_ptr->ConductivityModel = LEVEL_SET;
      num_const = read_constants(imp, &(mat_ptr->u_thermal_conductivity), 0);
      if (num_const < 3) {
        sprintf(err_msg, "Material %s - expected at least 3 constants for %s %s model.\n",
                pd_glob[mn]->MaterialName, search_string, "LEVEL_SET");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      if (mat_ptr->u_thermal_conductivity[2] == 0.0)
        mat_ptr->u_thermal_conductivity[2] = ls->Length_Scale / 2.0;

      mat_ptr->len_u_thermal_conductivity = num_const;

      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_thermal_conductivity);
    } else if (!strcmp(model_name, "THERMAL")) {
      mat_ptr->ConductivityModel = THERMAL_HEAT;
      num_const = read_constants(imp, &(mat_ptr->u_thermal_conductivity), 0);
      if (num_const < 5) {
        sprintf(err_msg, "Material %s - expected at least 5 constants for %s %s model.\n",
                pd_glob[mn]->MaterialName, search_string, "THERMAL");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      mat_ptr->len_u_thermal_conductivity = num_const;

      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_thermal_conductivity);
    }

    else if (!strcmp(model_name, "FOAM_PBE")) {
      mat_ptr->ConductivityModel = FOAM_PBE;
      num_const = read_constants(imp, &(mat_ptr->u_thermal_conductivity), 0);
      /* if (num_const < 5)  */
      /*   { */
      /*     sprintf(err_msg,  */
      /* 	      "Material %s - expected at least 5 constants for %s %s model.\n", */
      /* 	      pd_glob[mn]->MaterialName, search_string, "THERMAL"); */
      /*     GOMA_EH(GOMA_ERROR, err_msg); */
      /*   } */
      mat_ptr->len_u_thermal_conductivity = num_const;

      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_thermal_conductivity);

    } else if (!strcmp(model_name, "FOAM_PMDI_10")) {
      mat_ptr->ConductivityModel = FOAM_PMDI_10;
      num_const = read_constants(imp, &(mat_ptr->u_thermal_conductivity), 0);
      if (num_const < 2) {
        sprintf(err_msg, "Material %s - expected at least 2 constants for %s %s model.\n",
                pd_glob[mn]->MaterialName, search_string, "FOAM_PMDI_10");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      mat_ptr->len_u_thermal_conductivity = num_const;

      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_thermal_conductivity);
    } else if (!strcmp(model_name, "EXTERNAL_FIELD")) {
      if (fscanf(imp, "%s", input) != 1) {
        GOMA_EH(GOMA_ERROR,
                "Expecting trailing keyword for Thermal Conductivity EXTERNAL_FIELD model.\n");
      }
      ii = 0;
      for (j = 0; j < efv->Num_external_field; j++) {
        if (!strcmp(efv->name[j], input)) {
          ii = 1;
          mat_ptr->thermal_cond_external_field = j;
          break;
        }
      }
      if (ii == 0) {
        GOMA_EH(GOMA_ERROR, "Cannot match the name with that in the external field file");
      }
      mat_ptr->ConductivityModel = EXTERNAL_FIELD;
      /* pick up scale factor for property */
      num_const = read_constants(imp, &(mat_ptr->u_thermal_conductivity), NO_SPECIES);
      mat_ptr->len_u_thermal_conductivity = num_const;
      if (num_const < 1) {
        sr = sprintf(err_msg, "Matl %s expected at least 1 constant for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Thermal Conductivity", "EXTERNAL_FIELD");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
    } else {
      SPF(err_msg, "%s card read error.  Card missing or unknown model.", search_string);
      GOMA_EH(GOMA_ERROR, err_msg);
    }
  }

  ECHO(es, echo_file);

  strcpy(search_string, "Second Level Set Conductivity");

  model_read = look_for_mat_prop(imp, search_string, &(i0), &(a0), NO_USER, NULL, model_name,
                                 SCALAR_INPUT, &NO_SPECIES, es);

  if (model_read != -1) {

    if (ls == NULL)
      GOMA_EH(GOMA_ERROR,
              "Second Level Set Conductivity requires activation of Level Set Tracking.\n");

    mat_ptr->mp2nd->ThermalConductivityModel = i0;
    mat_ptr->mp2nd->thermalconductivity = a0;

    stringup(model_name);

    if (!strcmp(model_name, "CONSTANT")) {
      if (fscanf(imp, "%s", input) != 1) {
        GOMA_EH(GOMA_ERROR, "Expecting trailing keyword for Second Level Set Conductivity.\n");
      }

      stringup(input);

      if (strncmp(input, "POSITIVE", 3) == 0) {
        mat_ptr->mp2nd->thermalconductivitymask[0] = 0;
        mat_ptr->mp2nd->thermalconductivitymask[1] = 1;
      } else if (strncmp(input, "NEGATIVE", 3) == 0) {
        mat_ptr->mp2nd->thermalconductivitymask[0] = 1;
        mat_ptr->mp2nd->thermalconductivitymask[1] = 0;
      } else {
        GOMA_EH(GOMA_ERROR,
                "Keyword must be POSITIVE or NEGATIVE for Second Level Set Conductivity.\n");
      }
      SPF(endofstring(es), " %s", input);
      if (pfd != NULL) {
        for (i = 0; i < pfd->num_phase_funcs; i++) {
          if (fscanf(imp, "%lf", &(mat_ptr->mp2nd->thermalconductivity_phase[i])) != 1) {
            GOMA_EH(GOMA_ERROR, "error reading phase thermal conductivity");
          }
          SPF(endofstring(es), " %g", mat_ptr->mp2nd->thermalconductivity_phase[i]);
        }
      }
    } else {
      GOMA_EH(GOMA_ERROR, "Second Level Set Conductivity model can only be CONSTANT.\n");
    }
  } else {
    if (strlen(es) != 0) {
      SPF(err_msg, "Syntax error or unsupported model for %s ", search_string);
      GOMA_EH(GOMA_ERROR, err_msg);
    }
  }

  ECHO(es, echo_file);

  strcpy(search_string, "Heat Capacity");

  model_read = look_for_mat_proptable(
      imp, search_string, &(mat_ptr->HeatCapacityModel), &(mat_ptr->heat_capacity),
      &(mat_ptr->u_heat_capacity), &(mat_ptr->len_u_heat_capacity),
      &(mat_ptr->heat_capacity_tableid), model_name, SCALAR_INPUT, &NO_SPECIES, es);

  if (!strcmp(model_name, "ENTHALPY")) {
    mat_ptr->HeatCapacityModel = ENTHALPY;

    if (fscanf(imp, "%lf %lf", &(mat_ptr->heat_capacity), &(mat_ptr->latent_heat_fusion[0])) != 2) {
      GOMA_EH(GOMA_ERROR, "Expecting 2 floats for ENTHALPY model on Heat Capacity card");
    }

    SPF(endofstring(es), " %.4g %.4g", mat_ptr->heat_capacity, mat_ptr->latent_heat_fusion[0]);
  } else if (!strcmp(model_name, "CONST_LS") || !strcmp(model_name, "LEVEL_SET")) {
    mat_ptr->HeatCapacityModel = LEVEL_SET;
    num_const = read_constants(imp, &(mat_ptr->u_heat_capacity), 0);
    if (num_const < 3) {
      sprintf(err_msg, "Material %s - expected at least 3 constants for %s %s model.\n",
              pd_glob[mn]->MaterialName, search_string, "LEVEL_SET");
      GOMA_EH(GOMA_ERROR, err_msg);
    }

    if (mat_ptr->u_heat_capacity[2] == 0.0)
      mat_ptr->u_heat_capacity[2] = ls->Length_Scale / 2.0;

    mat_ptr->len_u_heat_capacity = num_const;

    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_heat_capacity);

  } else if (!strcmp(model_name, "THERMAL")) {
    mat_ptr->HeatCapacityModel = THERMAL_HEAT;
    num_const = read_constants(imp, &(mat_ptr->u_heat_capacity), 0);
    if (num_const < 5) {
      sprintf(err_msg, "Material %s - expected at least 5 constants for %s %s model.\n",
              pd_glob[mn]->MaterialName, search_string, "THERMAL");
      GOMA_EH(GOMA_ERROR, err_msg);
    }

    mat_ptr->len_u_heat_capacity = num_const;

    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_heat_capacity);

  } else if (!strcmp(model_name, "FOAM_PMDI_10")) {
    mat_ptr->HeatCapacityModel = FOAM_PMDI_10;
    num_const = read_constants(imp, &(mat_ptr->u_heat_capacity), 0);
    if (num_const < 2) {
      sprintf(err_msg, "Material %s - expected at least 2 constants for %s %s model.\n",
              pd_glob[mn]->MaterialName, search_string, "FOAM_PMDI_10");
      GOMA_EH(GOMA_ERROR, err_msg);
    }

    mat_ptr->len_u_heat_capacity = num_const;

    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_heat_capacity);

  }

  else {
    GOMA_EH(model_read, "Heat Capacity");
  }

  ECHO(es, echo_file);

  strcpy(search_string, "Second Level Set Heat Capacity");

  model_read = look_for_mat_prop(imp, search_string, &(i0), &(a0), NO_USER, NULL, model_name,
                                 SCALAR_INPUT, &NO_SPECIES, es);

  if (model_read != -1) {

    if (ls == NULL)
      GOMA_EH(GOMA_ERROR,
              "Second Level Set Heat Capacity requires activation of Level Set Tracking.\n");

    mat_ptr->mp2nd->HeatCapacityModel = i0;
    mat_ptr->mp2nd->heatcapacity = a0;

    stringup(model_name);

    if (!strcmp(model_name, "CONSTANT")) {
      if (fscanf(imp, "%s", input) != 1) {
        GOMA_EH(GOMA_ERROR, "Expecting trailing keyword for Second Level Set Heat Capacity.\n");
      }

      stringup(input);

      if (strncmp(input, "POSITIVE", 3) == 0) {
        mat_ptr->mp2nd->heatcapacitymask[0] = 0;
        mat_ptr->mp2nd->heatcapacitymask[1] = 1;
      } else if (strncmp(input, "NEGATIVE", 3) == 0) {
        mat_ptr->mp2nd->heatcapacitymask[0] = 1;
        mat_ptr->mp2nd->heatcapacitymask[1] = 0;
      } else {
        GOMA_EH(GOMA_ERROR,
                "Keyword must be POSITIVE or NEGATIVE for Second Level Set Heat Capacity.\n");
      }
      SPF(endofstring(es), " %s", input);
      if (pfd != NULL) {
        for (i = 0; i < pfd->num_phase_funcs; i++) {
          if (fscanf(imp, "%lf", &(mat_ptr->mp2nd->heatcapacity_phase[i])) != 1) {
            GOMA_EH(GOMA_ERROR, "error reading phase heat capacity");
          }
          SPF(endofstring(es), " %g", mat_ptr->mp2nd->heatcapacity_phase[i]);
        }
      }
    } else {
      GOMA_EH(GOMA_ERROR, "Second Level Set Heat Capacity model can only be CONSTANT.\n");
    }
  }

  ECHO(es, echo_file);

  strcpy(search_string, "Volume Expansion");
  model_read = look_for_mat_prop(imp, search_string, &(mat_ptr->VolumeExpansionModel),
                                 &(mat_ptr->Volume_Expansion), NO_USER, NULL, model_name,
                                 SCALAR_INPUT, &NO_SPECIES, es);

  if (model_read == -1) {
    if (strncmp(model_name, " ", 1) != 0) {
      SPF(err_msg, "Syntax error or invalid model for %s\n", search_string);
      GOMA_EH(GOMA_ERROR, err_msg);
    } else {
      SPF(es, "\t(%s = ) card is missing.  No default assigned.", search_string);
    }
  }
  ECHO(es, echo_file);

  strcpy(search_string, "Reference Temperature");
  model_read = look_for_mat_prop(imp, search_string, &(mat_ptr->ReferenceModel[TEMPERATURE]),
                                 &(mat_ptr->reference[TEMPERATURE]), NO_USER, NULL, model_name,
                                 SCALAR_INPUT, &NO_SPECIES, es);

  if (model_read == -1) {
    if (strncmp(model_name, " ", 1) != 0) {
      SPF(err_msg, "Syntax error or invalid model for %s\n", search_string);
      GOMA_EH(GOMA_ERROR, err_msg);
    } else {
      SPF(es, "\t(%s = ) card is missing.  No default assigned.", search_string);
    }
  }
  ECHO(es, echo_file);

  strcpy(search_string, "Liquidus Temperature");
  model_read = look_for_mat_prop(imp, search_string, &(mat_ptr->LiquidusModel),
                                 &(mat_ptr->melting_point_liquidus), NO_USER, NULL, model_name,
                                 SCALAR_INPUT, &NO_SPECIES, es);

  if (model_read == -1) {
    if (strncmp(model_name, " ", 1) != 0) {
      SPF(err_msg, "Syntax error or invalid model for %s\n", search_string);
      GOMA_EH(GOMA_ERROR, err_msg);
    } else {
      SPF(es, "\t(%s = ) card is missing.  No default assigned.", search_string);
    }
  }
  ECHO(es, echo_file);

  strcpy(search_string, "Solidus Temperature");
  model_read = look_for_mat_prop(imp, search_string, &(mat_ptr->SolidusModel),
                                 &(mat_ptr->melting_point_solidus), NO_USER, NULL, model_name,
                                 SCALAR_INPUT, &NO_SPECIES, es);

  if (model_read == -1) {
    if (strncmp(model_name, " ", 1) != 0) {
      SPF(err_msg, "Syntax error or invalid model for %s\n", search_string);
      GOMA_EH(GOMA_ERROR, err_msg);
    } else {
      SPF(es, "\t(%s = ) card is missing.  No default assigned.", search_string);
    }
  }
  ECHO(es, echo_file);

  strcpy(search_string, "Energy Weight Function");
  model_read =
      look_for_mat_prop(imp, search_string, &(mat_ptr->Ewt_funcModel), &(mat_ptr->Ewt_func),
                        NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);
  if (strncmp(model_name, " ", 1) != 0) {
    if (!strcmp(model_name, "GALERKIN")) {
      mat_ptr->Ewt_funcModel = GALERKIN;
      mat_ptr->Ewt_func = 0.;
    } else if (!strcmp(model_name, "SUPG")) {
      int err;
      mat_ptr->Ewt_funcModel = SUPG;
      err = fscanf(imp, "%lg", &(mat_ptr->Ewt_func));
      if (err != 1) {
        GOMA_EH(GOMA_ERROR, "Expected to read one double for Energy Weight Function SUPG");
      }
      SPF(endofstring(es), " %.4g", mat_ptr->Ewt_func);
    } else {
      SPF(err_msg, "Syntax error or invalid model for %s\n", search_string);
      GOMA_EH(GOMA_ERROR, err_msg);
    }
  } else {
    mat_ptr->Ewt_funcModel = GALERKIN;
    mat_ptr->Ewt_func = 0.;
    SPF(es, "\t(%s = %s)", search_string, "GALERKIN");
  }

  ECHO(es, echo_file);

  strcpy(search_string, "Energy Divergence Term");
  model_read =
      look_for_mat_prop(imp, search_string, &(mat_ptr->Ewt_funcModel), &(mat_ptr->Ewt_func),
                        NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);
  if (strncmp(model_name, " ", 1) != 0) {
    if (!strcmp(model_name, "yes") || !strcmp(model_name, "on")) {
      mat_ptr->Energy_Div_Term = 1;
    } else if (!strcmp(model_name, "no") || !strcmp(model_name, "off")) {
      mat_ptr->Energy_Div_Term = 0;
    } else {
      SPF(err_msg, "Syntax error or invalid model for %s\n", search_string);
      GOMA_EH(GOMA_ERROR, err_msg);
    }
  } else {
    mat_ptr->Energy_Div_Term = 0;
    SPF(es, "\t(%s = %s)", search_string, "off");
  }

  ECHO(es, echo_file);
  strcpy(search_string, "Residence Time Weight Function");
  model_read =
      look_for_mat_prop(imp, search_string, &(mat_ptr->Rst_funcModel), &(mat_ptr->Rst_func),
                        NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);
  if (model_read == -1) {
    if (!strcmp(model_name, "LINEAR_TIMETEMP")) {
      mat_ptr->Rst_funcModel = LINEAR_TIMETEMP;
    } else if (!strcmp(model_name, "EXPONENTIAL_TIMETEMP")) {
      mat_ptr->Rst_funcModel = EXPONENTIAL_TIMETEMP;
      int err = fscanf(imp, "%lg", &(mat_ptr->Rst_func));
      if (err != 1) {
        GOMA_EH(GOMA_ERROR, "Expected to read one double for Residence Time Weight Function");
      }
      SPF(endofstring(es), " %.4g", mat_ptr->Rst_func);
    } else if (!strcmp(model_name, "DROP_EVAP")) {
      mat_ptr->Rst_funcModel = DROP_EVAP;
    } else {
      mat_ptr->Rst_funcModel = CONSTANT;
      mat_ptr->Rst_func = 1.;
      mat_ptr->Rst_func_supg = 0.;
    }
    if (fscanf(imp, "%lg %lg %lg", &mat_ptr->Rst_func, &mat_ptr->Rst_diffusion,
               &mat_ptr->Rst_func_supg) != 3) {
      mat_ptr->Rst_func = 1.;
      mat_ptr->Rst_diffusion = 1. / LITTLE_PENALTY;
      mat_ptr->Rst_func_supg = 0.;
    }

    SPF(endofstring(es), " %.4g %.4g %.4g", mat_ptr->Rst_func, mat_ptr->Rst_diffusion,
        mat_ptr->Rst_func_supg);
  } else {
    mat_ptr->Rst_funcModel = CONSTANT;
    mat_ptr->Rst_func = 1.;
    mat_ptr->Rst_diffusion = 1. / LITTLE_PENALTY;
    mat_ptr->Rst_func_supg = 0.;
    SPF(es, "\t(%s = %s)", search_string, "CONSTANT");
  }
  ECHO(es, echo_file);

  /*
   * Electrical Properties
   *
   */

  ECHO("\n----Electrical Properties\n", echo_file);

  if (mn == 0)
    mat_ptr->elec_cond_external_field = -1;

  strcpy(search_string, "Electrical Conductivity");
  model_read = look_for_mat_prop(
      imp, search_string, &(mat_ptr->Elec_ConductivityModel), &(mat_ptr->electrical_conductivity),
      &(mat_ptr->u_electrical_conductivity), &(mat_ptr->len_u_electrical_conductivity), model_name,
      SCALAR_INPUT, &NO_SPECIES, es);

  if (model_read == -1) {
    if (!strcmp(model_name, "LEVEL_SET")) {
      mat_ptr->Elec_ConductivityModel = LEVEL_SET;

      num_const = read_constants(imp, &(mat_ptr->u_electrical_conductivity), 0);

      if (num_const < 3) {
        sr = sprintf(err_msg, "Matl %s expected at least 3 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Thermal Exponent", "LEVEL_SET");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      mat_ptr->len_u_electrical_conductivity = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_electrical_conductivity);

      if (mat_ptr->u_electrical_conductivity[2] == 0.0)
        mat_ptr->u_electrical_conductivity[2] = ls->Length_Scale / 2.0;

    } else if (strncmp(model_name, " ", 1) != 0) {
      if (!strcmp(model_name, "ELECTRODE_KINETICS")) {
        mat_ptr->Elec_ConductivityModel = ELECTRODE_KINETICS;
      } else if (!strcmp(model_name, "ELECTRONEUTRALITY_SM")) {
        mat_ptr->Elec_ConductivityModel = ELECTRONEUTRALITY_SM;
      } else if (!strcmp(model_name, "ELECTRONEUTRALITY_FICKIAN")) {
        mat_ptr->Elec_ConductivityModel = ELECTRONEUTRALITY_FICKIAN;
      } else if (!strcmp(model_name, "EXTERNAL_FIELD")) {
        if (fscanf(imp, "%s", input) != 1) {
          GOMA_EH(GOMA_ERROR,
                  "Expecting trailing keyword for Electrical Conductivity EXTERNAL_FIELD model.\n");
        }
        ii = 0;
        for (j = 0; j < efv->Num_external_field; j++) {
          if (!strcmp(efv->name[j], input)) {
            ii = 1;
            mat_ptr->elec_cond_external_field = j;
            break;
          }
        }
        if (ii == 0) {
          GOMA_EH(GOMA_ERROR, "Cannot match the name with that in the external field file");
        }
        mat_ptr->Elec_ConductivityModel = EXTERNAL_FIELD;

        /* pick up scale factor for property */
        num_const = read_constants(imp, &(mat_ptr->u_electrical_conductivity), NO_SPECIES);
        mat_ptr->len_u_electrical_conductivity = num_const;
        if (num_const < 1) {
          sr = sprintf(err_msg, "Matl %s expected at least 1 constant for %s %s model.\n",
                       pd_glob[mn]->MaterialName, "Electrical Conductivity", "EXTERNAL_FIELD");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
      }
    } else {
      /* no card defaulting case */
      mat_ptr->Elec_ConductivityModel = CONSTANT;
      mat_ptr->electrical_conductivity = 1.0;
      SPF(es, "\t(%s = %s %.4g)", search_string, "CONSTANT", mat_ptr->electrical_conductivity);
    }
  }
  ECHO(es, echo_file);

  /*
   * There is now a defined electrical permittivity, so it is no
   * longer necessary to use electrical conductivity in its place.
   */
  rewind(imp);
  strcpy(search_string, "Electrical Permittivity");
  model_read =
      look_for_mat_prop(imp, search_string, &(mat_ptr->PermittivityModel), &(mat_ptr->permittivity),
                        &(mat_ptr->u_permittivity), &(mat_ptr->len_u_permittivity), model_name,
                        SCALAR_INPUT, &NO_SPECIES, es);
  if (model_read == -1) {
    if (strcmp(model_name, "COMPLEX_CONSTANT") == 0) {
      mat_ptr->PermittivityModel = COMPLEX_CONSTANT;
      if (fscanf(imp, "%lf %lf", &(mat_ptr->permittivity), &(mat_ptr->permittivity_imag)) != 2) {
        GOMA_EH(GOMA_ERROR, "Expected 2 constants for %s = %s", search_string, model_name);
      }
    } else if (strcmp(model_name, "RADIAL_PML") == 0) {
      mat_ptr->PermittivityModel = RADIAL_PML;
      num_const = read_constants(imp, &(mat_ptr->u_permittivity), 0);
      if (num_const != 5) {
        GOMA_EH(GOMA_ERROR, "Expected 5 constants for %s = %s", search_string, model_name);
      }
      mat_ptr->len_u_permittivity = num_const;
    } else if (strcmp(model_name, "REFRACTIVE_INDEX") == 0) {
      mat_ptr->PermittivityModel = REFRACTIVE_INDEX;
      mat_ptr->len_u_permittivity = 0;
    } else if (strncmp(model_name, " ", 1) != 0) {

      SPF(err_msg, "Syntax error or invalid model for %s\n", search_string);
      GOMA_EH(GOMA_ERROR, err_msg);
    } else {
      mat_ptr->PermittivityModel = CONSTANT;
      mat_ptr->permittivity = 1.;
      SPF(es, "\t(%s = %s %.4g)", search_string, "CONSTANT", mat_ptr->permittivity);
    }
  }
  ECHO(es, echo_file);

  /*
   * There is now a defined magnetic permeability.
   */
  rewind(imp);
  strcpy(search_string, "Magnetic Permeability");
  model_read = look_for_mat_prop(
      imp, search_string, &(mat_ptr->MagneticPermeabilityModel), &(mat_ptr->magnetic_permeability),
      &(mat_ptr->u_magnetic_permeability), &(mat_ptr->len_u_magnetic_permeability), model_name,
      SCALAR_INPUT, &NO_SPECIES, es);
  if (model_read == -1) {
    if (strcmp(model_name, "COMPLEX_CONSTANT") == 0) {
      mat_ptr->PermeabilityModel = COMPLEX_CONSTANT;
      if (fscanf(imp, "%lf %lf", &(mat_ptr->permeability), &(mat_ptr->permeability_imag)) != 2) {
        GOMA_EH(GOMA_ERROR, "Expected 2 constants for %s = %s", search_string, model_name);
      }
    } else if (strcmp(model_name, "RADIAL_PML") == 0) {
      mat_ptr->PermeabilityModel = RADIAL_PML;
      num_const = read_constants(imp, &(mat_ptr->u_permeability), 0);
      if (num_const != 5) {
        GOMA_EH(GOMA_ERROR, "Expected 5 constants for %s = %s", search_string, model_name);
      }
      mat_ptr->len_u_permeability = num_const;
    } else if (strncmp(model_name, " ", 1) != 0) {
      SPF(err_msg, "Syntax error or invalid model for %s\n", search_string);
      GOMA_EH(GOMA_ERROR, err_msg);
    } else {
      mat_ptr->MagneticPermeabilityModel = CONSTANT;
      mat_ptr->magnetic_permeability = 1.;
      SPF(es, "\t(%s = %s %.4g)", search_string, "CONSTANT", mat_ptr->magnetic_permeability);
    }
  }
  ECHO(es, echo_file);

  strcpy(search_string, "Electrical Surface Diffusivity");
  model_read = look_for_mat_prop(
      imp, search_string, &(mat_ptr->Elect_Surf_DiffusivityModel),
      &(mat_ptr->elect_surf_diffusivity), &(mat_ptr->u_elect_surf_diffusivity),
      &(mat_ptr->len_u_elect_surf_diffusivity), model_name, SCALAR_INPUT, &NO_SPECIES, es);
  if (model_read == -1) {
    if (strncmp(model_name, " ", 1) != 0) {
      SPF(err_msg, "Syntax error or invalid model for %s\n", search_string);
      GOMA_EH(GOMA_ERROR, err_msg);
    } else {
      mat_ptr->Elect_Surf_DiffusivityModel = CONSTANT;
      mat_ptr->elect_surf_diffusivity = 0.;
      SPF(es, "\t(%s = %s %.4g)", search_string, "CONSTANT", mat_ptr->elect_surf_diffusivity);
    }
  }
  ECHO(es, echo_file);

  rewind(imp);
  strcpy(search_string, "Electromagnetic Incident Wave");
  model_read =
      look_for_mat_prop(imp, search_string, &(mat_ptr->IncidentWaveModel),
                        &(mat_ptr->incident_wave), &(mat_ptr->u_incident_wave),
                        &(mat_ptr->len_u_incident_wave), model_name, SCALAR_INPUT, &NO_SPECIES, es);
  if (model_read == -1) {
    if (strcmp(model_name, "PLANE_Z_WAVE") == 0) {
      mat_ptr->IncidentWaveModel = EM_INC_PLANE_Z_WAVE;
      num_const = read_constants(imp, &mat_ptr->u_incident_wave, 0);
      mat_ptr->len_u_incident_wave = num_const;
      if (num_const != 1) {
        GOMA_EH(GOMA_ERROR, "Expected 1 constants for %s = %s", search_string, model_name);
      }
    } else {
      mat_ptr->IncidentWaveModel = CONSTANT;
      mat_ptr->len_u_incident_wave = 0;
      SPF(es, "\t(%s = %s %.4g)", search_string, "NONE", mat_ptr->incident_wave);
    }
  }

  ECHO(es, echo_file);
  strcpy(search_string, "Shell User Parameter");
  model_read = look_for_mat_prop(imp, search_string, &(mat_ptr->Shell_User_ParModel),
                                 &(mat_ptr->shell_user_par), &(mat_ptr->u_shell_user_par),
                                 &(mat_ptr->len_u_shell_user_par), model_name, SCALAR_INPUT,
                                 &NO_SPECIES, es);
  if (model_read == -1) {
    if (strncmp(model_name, " ", 1) != 0) {
      SPF(err_msg, "Syntax error or invalid model for %s\n", search_string);
      GOMA_EH(GOMA_ERROR, err_msg);
    } else {
      mat_ptr->Shell_User_ParModel = CONSTANT;
      mat_ptr->shell_user_par = 0.;
      SPF(es, "\t(%s = %s %.4g)", search_string, "CONSTANT", mat_ptr->shell_user_par);
    }
  }
  ECHO(es, echo_file);

  /*
   * However, this card will allow a choice of the k value to be used in
   * assemble_potential(): electrical conductivity(default) or permittivity.
   */

  strcpy(search_string, "Voltage Formulation");
  model_read = look_for_optional(imp, search_string, input, '=');

  if (model_read == 1) {
    (void)read_string(imp, input, '\n');
    strip(input);
    stringup(input);

    if (strcmp(input, "PERMITTIVITY") == 0) {
      mat_ptr->VoltageFormulation = V_PERMITTIVITY;
      /* However, this will require a CONSTANT electrical conductivity */
      if (mat_ptr->Elec_ConductivityModel != CONSTANT) {
        GOMA_EH(
            GOMA_ERROR,
            "permittivity voltage formulation requires a CONSTANT electrical conductivity model!");
      }
    } else if (strcmp(input, "CONDUCTIVITY") == 0) {
      mat_ptr->VoltageFormulation = V_CONDUCTIVITY;
    } else {
      GOMA_EH(GOMA_ERROR, "Unknown option for voltage formulation!");
    }
    SPF(es, "%s = %s", search_string, input);
  } else {
    mat_ptr->VoltageFormulation = V_CONDUCTIVITY;
    SPF(es, "\t(%s = %s)", search_string, "CONDUCTIVITY");
  }

  ECHO(es, echo_file);

  /*
   * Acoustic Properties
   *
   */

  ECHO("\n----Acoustic Properties\n", echo_file);

  strcpy(search_string, "Acoustic Wave Number");
  model_read = look_for_mat_proptable(
      imp, search_string, &(mat_ptr->wave_numberModel), &(mat_ptr->wave_number),
      &(mat_ptr->u_wave_number), &(mat_ptr->len_u_wave_number), &(mat_ptr->wave_number_tableid),
      model_name, SCALAR_INPUT, &NO_SPECIES, es);
  if (model_read == -1) {
    if (!strcmp(model_name, "CONST_LS") || !strcmp(model_name, "LEVEL_SET")) {
      mat_ptr->wave_numberModel = LEVEL_SET;
      num_const = read_constants(imp, &(mat_ptr->u_wave_number), 0);
      if (num_const < 3) {
        sprintf(err_msg, "Material %s - expected at least 3 constants for %s %s model.\n",
                pd_glob[mn]->MaterialName, "Acoustic Wave Number", "LEVEL_SET");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      if (mat_ptr->u_wave_number[2] == 0.0)
        mat_ptr->u_wave_number[2] = ls->Length_Scale / 2.0;

      mat_ptr->len_u_wave_number = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_wave_number);

    } else if (strncmp(model_name, " ", 1) == 0) {
      mat_ptr->wave_numberModel = CONSTANT;
      mat_ptr->wave_number = 1.;
      SPF(es, "\t(%s = %s %.4g)", search_string, "CONSTANT", mat_ptr->wave_number);
    } else {
      SPF(err_msg, "Invalid model or syntax error for %s", search_string);
      GOMA_EH(GOMA_ERROR, err_msg);
    }
  }

  ECHO(es, echo_file);

  strcpy(search_string, "Second Level Set Acoustic Wave Number");
  model_read = look_for_mat_prop(imp, search_string, &(i0), &(a0), NO_USER, NULL, model_name,
                                 SCALAR_INPUT, &NO_SPECIES, es);
  if (model_read != -1) {
    if (ls == NULL)
      GOMA_EH(GOMA_ERROR,
              "Second Level Set Acoustic Wave Number requires activation of Level Set Tracking.\n");

    mat_ptr->mp2nd->wavenumberModel = i0;
    mat_ptr->mp2nd->wavenumber = a0;

    stringup(model_name);

    if (!strcmp(model_name, "CONSTANT")) {
      if (fscanf(imp, "%s", input) != 1) {
        GOMA_EH(GOMA_ERROR, "Expecting trailing keyword for Second Level Set Wave Number.\n");
      }

      stringup(input);

      if (strncmp(input, "POSITIVE", 3) == 0) {
        mat_ptr->mp2nd->wavenumbermask[0] = 0;
        mat_ptr->mp2nd->wavenumbermask[1] = 1;
      } else if (strncmp(input, "NEGATIVE", 3) == 0) {
        mat_ptr->mp2nd->wavenumbermask[0] = 1;
        mat_ptr->mp2nd->wavenumbermask[1] = 0;
      } else {
        GOMA_EH(GOMA_ERROR,
                "Keyword must be POSITIVE or NEGATIVE for Second Level Set Wave Number.\n");
      }
      SPF(endofstring(es), " %s", input);
      if (pfd != NULL) {
        for (i = 0; i < pfd->num_phase_funcs; i++) {
          if (fscanf(imp, "%lf", &(mat_ptr->mp2nd->wavenumber_phase[i])) != 1) {
            GOMA_EH(GOMA_ERROR, "error reading phase wave number");
          }
          SPF(endofstring(es), " %g", mat_ptr->mp2nd->wavenumber_phase[i]);
        }
      }
    } else {
      GOMA_EH(GOMA_ERROR, "Second Level Set Wave Number model can only be CONSTANT.\n");
    }
  } else if (strncmp(model_name, " ", 1) != 0) {
    SPF(err_msg, "Syntax error or invalid model for %s", search_string);
    GOMA_EH(GOMA_ERROR, err_msg);
  }

  ECHO(es, echo_file);

  strcpy(search_string, "Acoustic Impedance");

  model_read = look_for_mat_proptable(
      imp, search_string, &(mat_ptr->Acoustic_ImpedanceModel), &(mat_ptr->acoustic_impedance),
      &(mat_ptr->u_acoustic_impedance), &(mat_ptr->len_u_acoustic_impedance),
      &(mat_ptr->acoustic_impedance_tableid), model_name, SCALAR_INPUT, &NO_SPECIES, es);
  if (model_read == -1) {
    if (!strcmp(model_name, "CONST_LS") || !strcmp(model_name, "LEVEL_SET")) {
      mat_ptr->Acoustic_ImpedanceModel = LEVEL_SET;
      num_const = read_constants(imp, &(mat_ptr->u_acoustic_impedance), 0);
      if (num_const < 3) {
        sprintf(err_msg, "Material %s - expected at least 3 constants for %s %s model.\n",
                pd_glob[mn]->MaterialName, "Acoustic Impedance", "LEVEL_SET");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      if (mat_ptr->u_acoustic_impedance[2] == 0.0)
        mat_ptr->u_acoustic_impedance[2] = ls->Length_Scale / 2.0;

      mat_ptr->len_u_acoustic_impedance = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_acoustic_impedance);
    } else {
      mat_ptr->Acoustic_ImpedanceModel = CONSTANT;
      mat_ptr->acoustic_impedance = 1.;
      SPF(es, "\t(%s = %s %.4g)", search_string, "CONSTANT", mat_ptr->acoustic_impedance);
    }
  }

  ECHO(es, echo_file);

  strcpy(search_string, "Second Level Set Acoustic Impedance");

  model_read = look_for_mat_prop(imp, search_string, &(i0), &(a0), NO_USER, NULL, model_name,
                                 SCALAR_INPUT, &NO_SPECIES, es);

  if (model_read == 1) {

    if (ls == NULL)
      GOMA_EH(GOMA_ERROR,
              "Second Level Set Acoustic Impedance requires activation of Level Set Tracking.\n");

    mat_ptr->mp2nd->AcousticImpedanceModel = i0;
    mat_ptr->mp2nd->acousticimpedance = a0;

    stringup(model_name);

    if (!strcmp(model_name, "CONSTANT")) {
      if (fscanf(imp, "%s", input) != 1) {
        GOMA_EH(GOMA_ERROR,
                "Expecting trailing keyword for Second Level Set Acoustic Impedance.\n");
      }

      stringup(input);

      if (strncmp(input, "POSITIVE", 3) == 0) {
        mat_ptr->mp2nd->acousticimpedancemask[0] = 0;
        mat_ptr->mp2nd->acousticimpedancemask[1] = 1;
      } else if (strncmp(input, "NEGATIVE", 3) == 0) {
        mat_ptr->mp2nd->acousticimpedancemask[0] = 1;
        mat_ptr->mp2nd->acousticimpedancemask[1] = 0;
      } else {
        GOMA_EH(GOMA_ERROR,
                "Keyword must be POSITIVE or NEGATIVE for Second Level Set Acoustic Impedance.\n");
      }
      SPF(endofstring(es), " %s", input);
      if (pfd != NULL) {
        for (i = 0; i < pfd->num_phase_funcs; i++) {
          if (fscanf(imp, "%lf", &(mat_ptr->mp2nd->acousticimpedance_phase[i])) != 1) {
            GOMA_EH(GOMA_ERROR, "error reading phase acoustic impedance");
          }
          SPF(endofstring(es), " %g", mat_ptr->mp2nd->acousticimpedance_phase[i]);
        }
      }
    } else {
      GOMA_EH(GOMA_ERROR, "Second Level Set Acoustic Impedance model can only be CONSTANT.\n");
    }
  } else if (strncmp(model_name, " ", 1) != 0) {
    SPF(err_msg, "Syntax error or invalid model for %s", search_string);
    GOMA_EH(GOMA_ERROR, err_msg);
  }

  ECHO(es, echo_file);

  strcpy(search_string, "Acoustic Absorption");

  model_read = look_for_mat_proptable(
      imp, search_string, &(mat_ptr->Acoustic_AbsorptionModel), &(mat_ptr->acoustic_absorption),
      &(mat_ptr->u_acoustic_absorption), &(mat_ptr->len_u_acoustic_absorption),
      &(mat_ptr->acoustic_absorption_tableid), model_name, SCALAR_INPUT, &NO_SPECIES, es);
  if (model_read == -1) {
    if (!strcmp(model_name, "CONST_LS") || !strcmp(model_name, "LEVEL_SET")) {
      mat_ptr->Acoustic_AbsorptionModel = LEVEL_SET;
      num_const = read_constants(imp, &(mat_ptr->u_acoustic_absorption), 0);
      if (num_const < 3) {
        sprintf(err_msg, "Material %s - expected at least 3 constants for %s %s model.\n",
                pd_glob[mn]->MaterialName, "Acoustic Absorption", "LEVEL_SET");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      if (mat_ptr->u_acoustic_absorption[2] == 0.0)
        mat_ptr->u_acoustic_absorption[2] = ls->Length_Scale / 2.0;

      mat_ptr->len_u_acoustic_absorption = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_acoustic_absorption);
    } else {
      mat_ptr->Acoustic_AbsorptionModel = CONSTANT;
      mat_ptr->acoustic_absorption = 1.;
      SPF(es, "\t(%s = %s %.4g)", search_string, "CONSTANT", mat_ptr->acoustic_absorption);
    }
  }

  ECHO(es, echo_file);

  strcpy(search_string, "Second Level Set Acoustic Absorption");

  model_read = look_for_mat_prop(imp, search_string, &(i0), &(a0), NO_USER, NULL, model_name,
                                 SCALAR_INPUT, &NO_SPECIES, es);

  if (model_read == 1) {

    if (ls == NULL)
      GOMA_EH(GOMA_ERROR,
              "Second Level Set Acoustic Absorption requires activation of Level Set Tracking.\n");

    mat_ptr->mp2nd->AcousticAbsorptionModel = i0;
    mat_ptr->mp2nd->acousticabsorption = a0;

    stringup(model_name);

    if (!strcmp(model_name, "CONSTANT")) {
      if (fscanf(imp, "%s", input) != 1) {
        GOMA_EH(GOMA_ERROR,
                "Expecting trailing keyword for Second Level Set Acoustic Absorption.\n");
      }

      stringup(input);

      if (strncmp(input, "POSITIVE", 3) == 0) {
        mat_ptr->mp2nd->acousticabsorptionmask[0] = 0;
        mat_ptr->mp2nd->acousticabsorptionmask[1] = 1;
      } else if (strncmp(input, "NEGATIVE", 3) == 0) {
        mat_ptr->mp2nd->acousticabsorptionmask[0] = 1;
        mat_ptr->mp2nd->acousticabsorptionmask[1] = 0;
      } else {
        GOMA_EH(GOMA_ERROR,
                "Keyword must be POSITIVE or NEGATIVE for Second Level Set Acoustic Absorption.\n");
      }
      SPF(endofstring(es), " %s", input);
      if (pfd != NULL) {
        for (i = 0; i < pfd->num_phase_funcs; i++) {
          if (fscanf(imp, "%lf", &(mat_ptr->mp2nd->acousticabsorption_phase[i])) != 1) {
            GOMA_EH(GOMA_ERROR, "error reading phase acoustic absorption");
          }
          SPF(endofstring(es), " %g", mat_ptr->mp2nd->acousticabsorption_phase[i]);
        }
      }
    } else {
      GOMA_EH(GOMA_ERROR, "Second Level Set Acoustic Absorption model can only be CONSTANT.\n");
    }
  } else if (strncmp(model_name, " ", 1) != 0) {
    SPF(err_msg, "Syntax error or invalid model for %s", search_string);
    GOMA_EH(GOMA_ERROR, err_msg);
  }

  ECHO(es, echo_file);

  strcpy(search_string, "Acoustic Ksquared Sign");

  model_read = look_for_mat_prop(imp, search_string, &(mat_ptr->Ksquared_SignModel),
                                 &(mat_ptr->acoustic_ksquared_sign), NO_USER, NULL, model_name,
                                 SCALAR_INPUT, &NO_SPECIES, es);

  if (model_read == -1) {
    mat_ptr->Ksquared_SignModel = CONSTANT;
    mat_ptr->acoustic_ksquared_sign = 1.0;

    SPF(es, "\t(%s = CONSTANT %.4g)", "Acoustic Ksquared Sign", 1.0);
  }

  ECHO(es, echo_file);

  strcpy(search_string, "Refractive Index");

  model_read = look_for_mat_proptable(
      imp, search_string, &(mat_ptr->Refractive_IndexModel), &(mat_ptr->refractive_index),
      &(mat_ptr->u_refractive_index), &(mat_ptr->len_u_refractive_index),
      &(mat_ptr->refractive_index_tableid), model_name, SCALAR_INPUT, &NO_SPECIES, es);
  if (model_read == -1) {
    if (!strcmp(model_name, "CONST_LS") || !strcmp(model_name, "LEVEL_SET")) {
      mat_ptr->Refractive_IndexModel = LEVEL_SET;
      num_const = read_constants(imp, &(mat_ptr->u_refractive_index), 0);
      if (num_const < 3) {
        sprintf(err_msg, "Material %s - expected at least 3 constants for %s %s model.\n",
                pd_glob[mn]->MaterialName, "Refractive Index", "LEVEL_SET");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      if (mat_ptr->u_refractive_index[2] == 0.0)
        mat_ptr->u_refractive_index[2] = ls->Length_Scale / 2.0;

      mat_ptr->len_u_refractive_index = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_refractive_index);
    } else {
      mat_ptr->Refractive_IndexModel = CONSTANT;
      mat_ptr->refractive_index = 1.;
      SPF(es, "\t(%s = %s %.4g)", search_string, "CONSTANT", mat_ptr->refractive_index);
    }
  }

  ECHO(es, echo_file);

  strcpy(search_string, "Second Level Set Refractive Index");

  model_read = look_for_mat_prop(imp, search_string, &(i0), &(a0), NO_USER, NULL, model_name,
                                 SCALAR_INPUT, &NO_SPECIES, es);

  if (model_read == 1) {

    if (ls == NULL)
      GOMA_EH(GOMA_ERROR,
              "Second Level Set Refractive Index requires activation of Level Set Tracking.\n");

    mat_ptr->mp2nd->RefractiveIndexModel = i0;
    mat_ptr->mp2nd->refractiveindex = a0;

    stringup(model_name);

    if (!strcmp(model_name, "CONSTANT")) {
      if (fscanf(imp, "%s", input) != 1) {
        GOMA_EH(GOMA_ERROR, "Expecting trailing keyword for Second Level Set Refractive Index.\n");
      }

      stringup(input);

      if (strncmp(input, "POSITIVE", 3) == 0) {
        mat_ptr->mp2nd->refractiveindexmask[0] = 0;
        mat_ptr->mp2nd->refractiveindexmask[1] = 1;
      } else if (strncmp(input, "NEGATIVE", 3) == 0) {
        mat_ptr->mp2nd->refractiveindexmask[0] = 1;
        mat_ptr->mp2nd->refractiveindexmask[1] = 0;
      } else {
        GOMA_EH(GOMA_ERROR,
                "Keyword must be POSITIVE or NEGATIVE for Second Level Set Refractive Index.\n");
      }
      SPF(endofstring(es), " %s", input);
      if (pfd != NULL) {
        for (i = 0; i < pfd->num_phase_funcs; i++) {
          if (fscanf(imp, "%lf", &(mat_ptr->mp2nd->refractiveindex_phase[i])) != 1) {
            GOMA_EH(GOMA_ERROR, "error reading phase refractive index");
          }
          SPF(endofstring(es), " %g", mat_ptr->mp2nd->refractiveindex_phase[i]);
        }
      }
    } else {
      GOMA_EH(GOMA_ERROR, "Second Level Set Refractive Index model can only be CONSTANT.\n");
    }
  } else if (strncmp(model_name, " ", 1) != 0) {
    SPF(err_msg, "Syntax error or invalid model for %s", search_string);
    GOMA_EH(GOMA_ERROR, err_msg);
  }

  ECHO(es, echo_file);

  strcpy(search_string, "Light Absorption");

  model_read = look_for_mat_proptable(
      imp, search_string, &(mat_ptr->Light_AbsorptionModel), &(mat_ptr->light_absorption),
      &(mat_ptr->u_light_absorption), &(mat_ptr->len_u_light_absorption),
      &(mat_ptr->light_absorption_tableid), model_name, SCALAR_INPUT, &NO_SPECIES, es);
  if (model_read == -1) {
    if (!strcmp(model_name, "CONST_LS") || !strcmp(model_name, "LEVEL_SET")) {
      mat_ptr->Light_AbsorptionModel = LEVEL_SET;
      num_const = read_constants(imp, &(mat_ptr->u_light_absorption), 0);
      if (num_const < 3) {
        sprintf(err_msg, "Material %s - expected at least 3 constants for %s %s model.\n",
                pd_glob[mn]->MaterialName, "Light Absorption", "LEVEL_SET");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      if (mat_ptr->u_light_absorption[2] == 0.0)
        mat_ptr->u_light_absorption[2] = ls->Length_Scale / 2.0;

      mat_ptr->len_u_light_absorption = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_light_absorption);
    } else {
      mat_ptr->Light_AbsorptionModel = CONSTANT;
      mat_ptr->light_absorption = 1.;
      SPF(es, "\t(%s = %s %.4g)", search_string, "CONSTANT", mat_ptr->light_absorption);
    }
  }

  ECHO(es, echo_file);

  strcpy(search_string, "Second Level Set Light Absorption");

  model_read = look_for_mat_prop(imp, search_string, &(i0), &(a0), NO_USER, NULL, model_name,
                                 SCALAR_INPUT, &NO_SPECIES, es);

  if (model_read == 1) {

    if (ls == NULL)
      GOMA_EH(GOMA_ERROR,
              "Second Level Set Light Absorption requires activation of Level Set Tracking.\n");

    mat_ptr->mp2nd->LightAbsorptionModel = i0;
    mat_ptr->mp2nd->lightabsorption = a0;

    stringup(model_name);

    if (!strcmp(model_name, "CONSTANT")) {
      if (fscanf(imp, "%s", input) != 1) {
        GOMA_EH(GOMA_ERROR, "Expecting trailing keyword for Second Level Set Light Absorption.\n");
      }

      stringup(input);

      if (strncmp(input, "POSITIVE", 3) == 0) {
        mat_ptr->mp2nd->lightabsorptionmask[0] = 0;
        mat_ptr->mp2nd->lightabsorptionmask[1] = 1;
      } else if (strncmp(input, "NEGATIVE", 3) == 0) {
        mat_ptr->mp2nd->lightabsorptionmask[0] = 1;
        mat_ptr->mp2nd->lightabsorptionmask[1] = 0;
      } else {
        GOMA_EH(GOMA_ERROR,
                "Keyword must be POSITIVE or NEGATIVE for Second Level Set Light Absorption.\n");
      }
      SPF(endofstring(es), " %s", input);
      if (pfd != NULL) {
        for (i = 0; i < pfd->num_phase_funcs; i++) {
          if (fscanf(imp, "%lf", &(mat_ptr->mp2nd->lightabsorption_phase[i])) != 1) {
            GOMA_EH(GOMA_ERROR, "error reading phase light absorption");
          }
          SPF(endofstring(es), " %g", mat_ptr->mp2nd->lightabsorption_phase[i]);
        }
      }
    } else {
      GOMA_EH(GOMA_ERROR, "Second Level Set Light Absorption model can only be CONSTANT.\n");
    }
  } else if (strncmp(model_name, " ", 1) != 0) {
    SPF(err_msg, "Syntax error or invalid model for %s", search_string);
    GOMA_EH(GOMA_ERROR, err_msg);
  }

  ECHO(es, echo_file);

  strcpy(search_string, "Extinction Index");

  model_read = look_for_mat_proptable(
      imp, search_string, &(mat_ptr->Extinction_IndexModel), &(mat_ptr->extinction_index),
      &(mat_ptr->u_extinction_index), &(mat_ptr->len_u_extinction_index),
      &(mat_ptr->extinction_index_tableid), model_name, SCALAR_INPUT, &NO_SPECIES, es);
  if (model_read == -1) {
    if (!strcmp(model_name, "CONST_LS") || !strcmp(model_name, "LEVEL_SET")) {
      mat_ptr->Extinction_IndexModel = LEVEL_SET;
      num_const = read_constants(imp, &(mat_ptr->u_extinction_index), 0);
      if (num_const < 3) {
        sprintf(err_msg, "Material %s - expected at least 3 constants for %s %s model.\n",
                pd_glob[mn]->MaterialName, "Extinction Index", "LEVEL_SET");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      if (mat_ptr->u_extinction_index[2] == 0.0)
        mat_ptr->u_extinction_index[2] = ls->Length_Scale / 2.0;

      mat_ptr->len_u_extinction_index = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_extinction_index);
    } else {
      mat_ptr->Extinction_IndexModel = CONSTANT;
      mat_ptr->extinction_index = 0.;
      SPF(es, "\t(%s = %s %.4g)", search_string, "CONSTANT", mat_ptr->extinction_index);
    }
  }

  ECHO(es, echo_file);

  strcpy(search_string, "Second Level Set Extinction Index");

  model_read = look_for_mat_prop(imp, search_string, &(i0), &(a0), NO_USER, NULL, model_name,
                                 SCALAR_INPUT, &NO_SPECIES, es);

  if (model_read == 1) {

    if (ls == NULL)
      GOMA_EH(GOMA_ERROR,
              "Second Level Set Extinction Index requires activation of Level Set Tracking.\n");

    mat_ptr->mp2nd->ExtinctionIndexModel = i0;
    mat_ptr->mp2nd->extinctionindex = a0;

    stringup(model_name);

    if (!strcmp(model_name, "CONSTANT")) {
      if (fscanf(imp, "%s", input) != 1) {
        GOMA_EH(GOMA_ERROR, "Expecting trailing keyword for Second Level Set Extinction Index.\n");
      }

      stringup(input);

      if (strncmp(input, "POSITIVE", 3) == 0) {
        mat_ptr->mp2nd->extinctionindexmask[0] = 0;
        mat_ptr->mp2nd->extinctionindexmask[1] = 1;
      } else if (strncmp(input, "NEGATIVE", 3) == 0) {
        mat_ptr->mp2nd->extinctionindexmask[0] = 1;
        mat_ptr->mp2nd->extinctionindexmask[1] = 0;
      } else {
        GOMA_EH(GOMA_ERROR,
                "Keyword must be POSITIVE or NEGATIVE for Second Level Set Extinction Index.\n");
      }
      SPF(endofstring(es), " %s", input);
      if (pfd != NULL) {
        for (i = 0; i < pfd->num_phase_funcs; i++) {
          if (fscanf(imp, "%lf", &(mat_ptr->mp2nd->extinctionindex_phase[i])) != 1) {
            GOMA_EH(GOMA_ERROR, "error reading phase extinction index");
          }
          SPF(endofstring(es), " %g", mat_ptr->mp2nd->extinctionindex_phase[i]);
        }
      }
    } else {
      GOMA_EH(GOMA_ERROR, "Second Level Set Extinction Index model can only be CONSTANT.\n");
    }
  } else if (strncmp(model_name, " ", 1) != 0) {
    SPF(err_msg, "Syntax error or invalid model for %s", search_string);
    GOMA_EH(GOMA_ERROR, err_msg);
  }

  ECHO(es, echo_file);

  /* New porous section */
  /*To avoid the ordering
   * complaints, let's rewind the mat file and start looking from the
   * beginning.
   */

  ECHO("\n----Porous Media Properties\n", echo_file);

  rewind(imp);
  mat_ptr->porosity_external_field_index = -1;
  mat_ptr->perm_external_field_index = -1;
  mat_ptr->Xperm_external_field_index = -1;
  mat_ptr->rel_liq_perm_external_field_index = -1;
  efv->ev_dpordt_index = -1;
  mat_ptr->SAT_external_field_index = -1;
  mat_ptr->por_shell_closed_porosity_ext_field_index = -1;
  mat_ptr->por_shell_closed_height_ext_field_index = -1;
  mat_ptr->por_shell_closed_radius_ext_field_index = -1;

  /* This initialization is for porous shell equations */
  for (i = 0; i < MAX_POR_SHELL; i++) {
    mat_ptr->por_shell_porosity_ext_field_index[i] = -1;
    mat_ptr->por_shell_height_ext_field_index[i] = -1;
    mat_ptr->por_shell_permeability_ext_field_index[i] = -1;
    mat_ptr->por_shell_rel_perm_ext_field_index[i] = -1;
    mat_ptr->por_shell_cap_pres_ext_field_index[i] = -1;
    mat_ptr->por_shell_cap_pres_hyst_curve_type_ext_field_index[i] = -1;
    mat_ptr->por_shell_cap_pres_hyst_num_switch_ext_field_index[i] = -1;
  }

  /*
   * Porous Section 1:  Microstructure Properties.  These are inputs
   * that are really a function of solid skeleton phase, and also some very
   * specific Darcy-law constants (e.g. permeabilities) that are only accessible
   * experimentally, in some cases.
   *
   * Also input here is the compressibility of the solvent liquid phase, which is
   * taken to be constant regardless of composition.  This compressibility is required
   * for numerical stability on impregnation problems.
   *
   * Microstructure stuff ************************************
   */

  /* Check for existance of porous media variables in any matrix*/

  have_por_liq_pres = 0;
  have_por_gas_pres = 0;
  have_por_porosity = 0;
  have_por_sink_mass = 0;
  have_por_energy = 0;
  have_shell_sat_open = 0;
  have_shell_sat_open2 = 0;
  int have_shell_sat_n = 0;
  PorousShellOn = 0;

  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    if (pd_glob[mn]->e[imtrx][R_POR_LIQ_PRES]) {
      have_por_liq_pres = 1;
    }
    if (pd_glob[mn]->e[imtrx][R_POR_GAS_PRES]) {
      have_por_gas_pres = 1;
    }
    if (pd_glob[mn]->e[imtrx][R_POR_POROSITY]) {
      have_por_porosity = 1;
    }
    if (pd_glob[mn]->e[imtrx][R_POR_ENERGY]) {
      have_por_energy = 1;
    }
    if (pd_glob[mn]->e[imtrx][R_SHELL_SAT_OPEN]) {
      have_shell_sat_open = 1;
    }
    if (pd_glob[mn]->e[imtrx][R_SHELL_SAT_OPEN_2]) {
      have_shell_sat_open2 = 1;
    }
    if ((pd_glob[mn]->e[imtrx][R_SHELL_SAT_1]) || (pd_glob[mn]->e[imtrx][R_SHELL_SAT_2]) ||
        (pd_glob[mn]->e[imtrx][R_SHELL_SAT_3])) {
      PorousShellOn = 1;
      have_shell_sat_n = 1;
    }
  }

  strcpy(search_string, "Media Type");
  model_read = look_for_mat_prop(imp, search_string, &(mat_ptr->PorousMediaType), &(a0), NO_USER,
                                 NULL, model_name, NO_INPUT, &NO_SPECIES, es);
  if (model_read == -1 && (!strcmp(model_name, "CONTINUOUS") || !strcmp(model_name, "NONE"))) {
    mat_ptr->PorousMediaType = CONTINUOUS;
    /*consistency checks */
    if (have_por_liq_pres == 1) {
      SPF(err_msg, "CONTINUOUS model for %s cannot be used with porous equation\n", search_string);
      GOMA_EH(GOMA_ERROR, err_msg);
    }
  } else if (model_read == -1 && !strcmp(model_name, "POROUS_SATURATED")) {
    mat_ptr->PorousMediaType = POROUS_SATURATED;
    /*consistency checks */
    if (have_por_liq_pres == 0) {
      GOMA_EH(
          GOMA_ERROR,
          "You cannot run a porous media simulation without selecting the porous media equations");
    }
    mat_ptr->CapStress = NO_MODEL; /*so proper effective stress principle is chosen*/
  } else if (model_read == -1 && !strcmp(model_name, "POROUS_TWO_PHASE")) {
    mat_ptr->PorousMediaType = POROUS_TWO_PHASE;
    if (have_por_liq_pres == 0 || have_por_gas_pres == 0) {
      GOMA_EH(GOMA_ERROR, "You cannot run a two-phase porous media simulation without selecting "
                          "both the porous_liq and porous_gas equations");
    }
  } else if (model_read == -1 && (!strcmp(model_name, "POROUS_UNSATURATED") ||
                                  !strcmp(model_name, "POROUS_PART_SAT"))) {
    mat_ptr->PorousMediaType = POROUS_UNSATURATED;
    if (have_por_liq_pres == 0) {
      GOMA_EH(
          GOMA_ERROR,
          "You cannot run a porous media simulation without selecting the porous media equations");
    }
  } else if (model_read == -1 && (!strcmp(model_name, "POROUS_BRINKMAN"))) {
    mat_ptr->PorousMediaType = POROUS_BRINKMAN;
    mat_ptr->i_ys = 0;
  } else if (model_read == -1 && (!strcmp(model_name, "POROUS_SHELL_UNSATURATED"))) {
    mat_ptr->PorousMediaType = POROUS_SHELL_UNSATURATED;
    mat_ptr->i_ys = 0;
    if (!have_shell_sat_open && !have_shell_sat_open2 && !have_shell_sat_n) {
      GOMA_EH(GOMA_ERROR, "You cannot run a porous shell media simulation without selecting the "
                          "porous shell equation(s)");
    }
  } else if ((model_read == -1) && (strncmp(model_name, " ", 1) == 0)) {
    mat_ptr->PorousMediaType = CONTINUOUS;
    /*consistency checks */
    if (have_por_liq_pres == 1) {
      SPF(err_msg, "CONTINUOUS model for %s cannot be used with porous equation\n", search_string);
      GOMA_EH(GOMA_ERROR, err_msg);
    }

    SPF(es, "\t(%s = %s)", search_string, "CONTINUOUS");
  } else {
    SPF(err_msg, "Incorrect syntax or model choice for %s\n", search_string);
    GOMA_EH(GOMA_ERROR, err_msg);
  }

  ECHO(es, echo_file);

  if ((mat_ptr->PorousMediaType != CONTINUOUS) &&
      (!PorousShellOn)) /* This only applies to continuum porous media equations */
  {
    /* create a list of the porous media equations
       active in this material                     */
    i = 0;

    for (j = 0; j < MAX_POROUS_NUM; j++) {
      for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
        if (pd_glob[mn]->e[imtrx][R_POR_LIQ_PRES + j]) {
          mat_ptr->Porous_Eqn[i] = R_POR_LIQ_PRES + j;
          i++;
        }
      }
    }
    if (have_shell_sat_open == 1) {
      mat_ptr->Porous_Eqn[i] = R_SHELL_SAT_OPEN;
      i++;
    }

    if (have_shell_sat_open2 == 1) {
      mat_ptr->Porous_Eqn[i] = R_SHELL_SAT_OPEN_2;
      i++;
    }

    if (i != pd_glob[mn]->Num_Porous_Eqn) {
      GOMA_EH(GOMA_ERROR, "Possible duplicate porous media equations");
    }

    strcpy(search_string, "Porosity");
    model_read =
        look_for_mat_prop(imp, search_string, &(mat_ptr->PorosityModel), &(mat_ptr->porosity),
                          NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);

    if ((have_por_porosity == 1) && strcmp(model_name, "DEFORM")) {
      GOMA_EH(GOMA_ERROR, " you must have a DEFORM porosity model with the pore porosity equation");
    }

    if (model_read == -1 && !strcmp(model_name, "DEFORM")) {
      if (have_por_porosity ==
          0) /* OK OK i get it...no R_POR_POROSITY without DEFORM and vice versa */
      {
        GOMA_EH(GOMA_ERROR,
                "You cannot calculate porosity variation without the pore porosity equation");
      }

      mat_ptr->PorosityModel = DEFORM;
      num_const = read_constants(imp, &(mat_ptr->u_porosity), NO_SPECIES);
      if (num_const < 1) {
        sr = sprintf(err_msg, "Matl %s %s %s model needs initial undeformed porosity \n",
                     pd_glob[mn]->MaterialName, "Porosity", "DEFORM");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_porosity = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_porosity);

    } else if (model_read == -1 && !strcmp(model_name, "EXTERNAL_FIELD")) {
      if (fscanf(imp, "%s", input) != 1) {
        GOMA_EH(GOMA_ERROR, "Expecting trailing keyword for Porosity EXTERNAL_FIELD model.\n");
      }
      ii = 0;
      for (j = 0; j < efv->Num_external_field; j++) {
        if (!strcmp(efv->name[j], input)) {
          ii = 1;
          mat_ptr->porosity_external_field_index = j;
        }
      }
      if (ii == 0) {
        GOMA_EH(GOMA_ERROR, "Cannot match the name with that in the external field file");
      }
      mat_ptr->PorosityModel = EXTERNAL_FIELD;

      /* pick up scale factor for property */
      num_const = read_constants(imp, &(mat_ptr->u_porosity), NO_SPECIES);
      mat_ptr->len_u_porosity = num_const;
      if (num_const < 1) {
        sr = sprintf(err_msg, "Matl %s expected at least 1 constant for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Porosity", "EXTERNAL_FIELD");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
    } else {
      GOMA_EH(model_read, "Porosity: Is card missing?");
    }

    ECHO(es, echo_file);

    /* Print out Porous Media Equation numbers used by Code */

    SPF(es, "\n(Active Porous Equations: pl = %d, pg = %d, porosity = %d, pe = %d)\n",
        have_por_liq_pres, have_por_gas_pres, have_por_porosity, have_por_energy);

    ECHO(es, echo_file);

    /* Porous (rock) compressibility */

    strcpy(search_string, "Porous Compressibility");
    model_read = look_for_mat_prop(
        imp, search_string, &(mat_ptr->PorousCompressibilityModel),
        &(mat_ptr->porous_compressibility), &(mat_ptr->u_porous_compressibility),
        &(mat_ptr->len_u_porous_compressibility), model_name, SCALAR_INPUT, &NO_SPECIES, es);
    /* model CONST_INIT also reads in an initial porosity value */
    if (model_read == -1 && !strcmp(model_name, "CONST_INIT")) {
      mat_ptr->PorousCompressibilityModel = POROUS_CONST_INIT;
      num_const = read_constants(imp, &mat_ptr->u_porous_compressibility, 0);
      if (num_const < 2) {
        sprintf(err_msg, "Material %s - expected at least 2 constants for %s %s model.\n",
                pd_glob[mn]->MaterialName, "Porous Compressibility", "CONST_INIT");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->porous_compressibility = mat_ptr->u_porous_compressibility[0];
      mat_ptr->initial_porosity = mat_ptr->u_porous_compressibility[1];
      mat_ptr->len_u_porous_compressibility = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_porous_compressibility);
    } else {
#ifdef LIBRARY_MODE
      if (efv->ev_porous_decouple) {
        sprintf(err_msg, "Material %s - unrecognized model for %s \"%s\" ???\n",
                pd_glob[mn]->MaterialName, "Porous Compressibility", model_name);
        GOMA_EH(model_read, err_msg);
      }
#endif
    }

    ECHO(es, echo_file);

    strcpy(search_string, "Permeability");
    model_read = look_for_mat_prop(imp, search_string, &(mat_ptr->PermeabilityModel),
                                   &(mat_ptr->permeability), &(mat_ptr->u_permeability),
                                   &(mat_ptr->len_u_permeability), model_name, SCALAR_INPUT,
                                   &NO_SPECIES, es);
    if (model_read == -1 && !strcmp(model_name, "PSD_VOL")) {
      mat_ptr->PermeabilityModel = PSD_VOL;
      num_const = read_constants(imp, &(mat_ptr->u_permeability), NO_SPECIES);
      if (num_const < 4) {
        sr = sprintf(err_msg, "Matl %s expected at least 4 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Permeability", "PSD_VOL");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_permeability = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_permeability);
    } else if (model_read == -1 && !strcmp(model_name, "PSD_WEXP")) {
      mat_ptr->PermeabilityModel = PSD_WEXP;
      num_const = read_constants(imp, &(mat_ptr->u_permeability), NO_SPECIES);
      if (num_const < 4) {
        sr = sprintf(err_msg, "Matl %s expected at least 4 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Permeability", "PSD_WEXP");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_permeability = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_permeability);
    } else if (model_read == -1 && !strcmp(model_name, "PSD_SEXP")) {
      mat_ptr->PermeabilityModel = PSD_SEXP;
      num_const = read_constants(imp, &(mat_ptr->u_permeability), NO_SPECIES);
      if (num_const < 4) {
        sr = sprintf(err_msg, "Matl %s expected at least 4 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Permeability", "PSD_SEXP");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_permeability = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_permeability);
    } else if (model_read == -1 && !strcmp(model_name, "SOLIDIFICATION")) {
      mat_ptr->PermeabilityModel = SOLIDIFICATION;
      num_const = read_constants(imp, &(mat_ptr->u_permeability), NO_SPECIES);
      if (num_const < 1) {
        sr = sprintf(err_msg, "Matl %s expected at least 1 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Permeability", "SOLIDIFICATION");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      mat_ptr->len_u_permeability = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_permeability);
    } else if (model_read == -1 && !strcmp(model_name, "TENSOR")) {
      mat_ptr->PermeabilityModel = K_TENSOR;
      num_const = read_constants(imp, &(mat_ptr->u_permeability), NO_SPECIES);
      if (pd_glob[mn]->Num_Dim == 3)
        GOMA_EH(GOMA_ERROR, "Tensor Permeability only 2D for now");
      if (num_const < 4) {
        sr = sprintf(err_msg, "Matl %s expected at least 4 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Permeability", "TENSOR");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      mat_ptr->perm_tensor[0][0] = mat_ptr->u_permeability[0];
      mat_ptr->perm_tensor[1][1] = mat_ptr->u_permeability[1];
      mat_ptr->perm_tensor[1][0] = mat_ptr->u_permeability[2];
      mat_ptr->perm_tensor[0][1] = mat_ptr->u_permeability[3];
      mat_ptr->permeability = mat_ptr->u_permeability[0];

      mat_ptr->len_u_permeability = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_permeability);
    } else if (model_read == -1 && !strcmp(model_name, "KOZENY_CARMAN")) {
      mat_ptr->PermeabilityModel = KOZENY_CARMAN;
      num_const = read_constants(imp, &(mat_ptr->u_permeability), NO_SPECIES);
      if (num_const < 2) {
        GOMA_EH(GOMA_ERROR, "expected 2 constants for KOZENY_CARMAN perm model");
      }
      mat_ptr->len_u_permeability = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_permeability);
    } else if (model_read == -1 && !strcmp(model_name, "SINK_MASS_PERM")) {
      if (have_por_sink_mass == 0)
        GOMA_EH(GOMA_ERROR, "SINK_MASS_PERM model not permitted without sink eqn");

      mat_ptr->PermeabilityModel = SINK_MASS_PERM;
      num_const = read_constants(imp, &(mat_ptr->u_permeability), NO_SPECIES);
      if (num_const < 4) {
        GOMA_EH(GOMA_ERROR, "expected 4 constants for SINK_MASS_PERM perm model");
      }
      mat_ptr->len_u_permeability = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_permeability);
    } else if (model_read == -1 && !strcmp(model_name, "ORTHOTROPIC")) {
      mat_ptr->PermeabilityModel = ORTHOTROPIC;
      num_const = read_constants(imp, &(mat_ptr->u_permeability), NO_SPECIES);
      if (num_const < 12) {
        sr = sprintf(err_msg, "Matl %s expected at least 12 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Permeability", "ORTHOTROPIC");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      memset(mat_ptr->perm_tensor, 0, sizeof(double) * DIM * DIM); /*these are loaded up later */
      mat_ptr->permeability = mat_ptr->u_permeability[0];          /*just in case */

      mat_ptr->len_u_permeability = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_permeability);
    }

    else if (model_read == -1 && !strcmp(model_name, "KC_TENSOR")) {
      mat_ptr->PermeabilityModel = KC_TENSOR;
      num_const = read_constants(imp, &(mat_ptr->u_permeability), NO_SPECIES);
      if (num_const < 13) {
        sr = sprintf(err_msg, "Matl %s expected at least 13 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Permeability", "KC_TENSOR");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      memset(mat_ptr->perm_tensor, 0, sizeof(double) * DIM * DIM); /*these are loaded up later */

      mat_ptr->len_u_permeability = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_permeability);
    } else if (model_read == -1 && !strcmp(model_name, "SM_TENSOR")) {
      mat_ptr->PermeabilityModel = SM_TENSOR;
      num_const = read_constants(imp, &(mat_ptr->u_permeability), NO_SPECIES);
      if (num_const < 15) {
        sr = sprintf(err_msg, "Matl %s expected at least 15 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Permeability", "SM_TENSOR");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      memset(mat_ptr->perm_tensor, 0, sizeof(double) * DIM * DIM); /*these are loaded up later */

      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_permeability);
    } else if (model_read == -1 && !strcmp(model_name, "SHELL_CYLINDER_SQUARE")) {
      mat_ptr->PermeabilityModel = SHELL_CYLINDER_SQUARE;
      num_const = read_constants(imp, &(mat_ptr->u_permeability), NO_SPECIES);

      memset(mat_ptr->perm_tensor, 0, sizeof(double) * DIM * DIM); /*these are loaded up later */
      mat_ptr->len_u_permeability = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_permeability);
    }

    else if (model_read == -1 && !strcmp(model_name, "EXTERNAL_FIELD")) {

      if (fscanf(imp, "%s", input) != 1) {
        GOMA_EH(GOMA_ERROR, "Expecting trailing keyword for Permeability EXTERNAL_FIELD model.\n");
      }
      ii = 0;
      for (j = 0; j < efv->Num_external_field; j++) {
        if (!strcmp(efv->name[j], input)) {
          ii = 1;
          mat_ptr->perm_external_field_index = j;
        }
      }
      if (ii == 0) {
        GOMA_EH(GOMA_ERROR, "Must activate external fields to use this Permeability model and "
                            "requested name must be an external field name");
      }
      mat_ptr->PermeabilityModel = EXTERNAL_FIELD;

      /* pick up scale factor for property */
      num_const = read_constants(imp, &(mat_ptr->u_permeability), NO_SPECIES);
      mat_ptr->len_u_permeability = num_const;
      if (num_const < 1) {
        sr = sprintf(err_msg, "Matl %s expected at least 1 constant for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Rel Liq Permeability", "EXTERNAL_FIELD");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
    }

    else {
      GOMA_EH(model_read, "Permeability: Card missing or wrong model? ");
    }

    ECHO(es, echo_file);
    /*
     * Compressibility of liquid phase
     */
    mat_ptr->PorousLiqCompress = 0.;

    strcpy(search_string, "Liquid phase compressibility");
    model_read = look_for_mat_prop(imp, search_string, &(mat_ptr->PorousLiquidCompressModel),
                                   &(mat_ptr->PorousLiqCompress), NO_USER, NULL, model_name,
                                   SCALAR_INPUT, &NO_SPECIES, es);

    ECHO(es, echo_file);

    /*
     * Reference pressure for compressibility of liquid phase
     */
    strcpy(search_string, "Liquid phase reference pressure");
    mat_ptr->PorousLiqRefPress = 0.;
    model_read = look_for_mat_prop(imp, search_string, &(mat_ptr->PorousLiqRefPressModel),
                                   &(mat_ptr->PorousLiqRefPress), NO_USER, NULL, model_name,
                                   SCALAR_INPUT, &NO_SPECIES, es);
    ECHO(es, echo_file);
  } /*end of if( (mat_ptr->PorousMediaType != CONTINUOUS) && (!PorousShellOn) */

  /*
   * Read in the special Brinkman Equation Parameters:
   *    FlowingLiquid Viscosity
   *    Inertia Coefficent
   *    (porosity and permeability are obtained from above)
   */
  if (mat_ptr->PorousMediaType == POROUS_BRINKMAN) {
    strcpy(search_string, "FlowingLiquid Viscosity");
    model_read = look_for_mat_prop(
        imp, search_string, &(mat_ptr->FlowingLiquidViscosityModel),
        &(mat_ptr->FlowingLiquid_viscosity), &(mat_ptr->u_FlowingLiquid_viscosity),
        &(mat_ptr->len_u_FlowingLiquid_viscosity), model_name, SCALAR_INPUT, &NO_SPECIES, es);

    if (!strcmp(model_name, "MOLTEN_GLASS")) {
      mat_ptr->FlowingLiquidViscosityModel = MOLTEN_GLASS;
      num_const = read_constants(imp, &(mat_ptr->u_FlowingLiquid_viscosity), NO_SPECIES);
      if (num_const < 3) {
        sprintf(err_msg, "Matl %s (conc %d) needs at least 3 constants for %s %s model.\n",
                pd_glob[mn]->MaterialName, species_no, search_string, "MOLTEN GLASS");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_FlowingLiquid_viscosity = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_FlowingLiquid_viscosity);

    } else if (!strcmp(model_name, "EPOXY")) {
      mat_ptr->FlowingLiquidViscosityModel = EPOXY;
      num_const = read_constants(imp, &(mat_ptr->u_FlowingLiquid_viscosity), NO_SPECIES);
      if (num_const < 6) {
        sprintf(err_msg, "Matl %s (conc %d) needs at least 6 constants for %s %s model.\n",
                pd_glob[mn]->MaterialName, species_no, search_string, "EPOXY");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_FlowingLiquid_viscosity = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_FlowingLiquid_viscosity);

    } else {
      GOMA_EH(model_read, "FlowingLiquid Viscosity");
    }

    if ((mat_ptr->FlowingLiquid_viscosity != 0.) &&
        ((pd_glob[mn]->gv[R_LUBP]) || (pd_glob[mn]->gv[R_LUBP_2])))
      GOMA_WH(-1, "ON POROUS_BRINKMAN: You will get erroneous results if you are using the "
                  "brinkman formulation as an expedient for lubrication velocity calculation. "
                  "FlowingliquidViscosity should be zero in that case");

    ECHO(es, echo_file);

    model_read = look_for_mat_prop(imp, "Second Level Set FlowingLiquid Viscosity", &(i0), &(a0),
                                   NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);

    if (model_read != -1) {
      if (ls == NULL)
        GOMA_EH(GOMA_ERROR, "Second Level Set FlowingLiquid Viscosity requires activation of Level "
                            "Set Tracking.\n");

      mat_ptr->mp2nd->FlowingLiquidViscosityModel = i0;
      mat_ptr->mp2nd->FlowingLiquid_viscosity = a0;

      stringup(model_name);

      if (!strcmp(model_name, "CONSTANT")) {
        if (fscanf(imp, "%s", input) != 1) {
          GOMA_EH(GOMA_ERROR,
                  "Expecting trailing keyword for Second Level Set FlowingLiquid Viscosity.\n");
        }

        stringup(input);

        if (strncmp(input, "POSITIVE", 3) == 0) {
          mat_ptr->mp2nd->FlowingLiquid_viscositymask[0] = 0;
          mat_ptr->mp2nd->FlowingLiquid_viscositymask[1] = 1;
        } else if (strncmp(input, "NEGATIVE", 3) == 0) {
          mat_ptr->mp2nd->FlowingLiquid_viscositymask[0] = 1;
          mat_ptr->mp2nd->FlowingLiquid_viscositymask[1] = 0;
        } else {
          GOMA_EH(GOMA_ERROR, "Keyword must be POSITIVE or NEGATIVE for Second Level Set "
                              "FlowingLiquid Viscosity.\n");
        }
        SPF(endofstring(es), " %s", input);
        if (pfd != NULL) {
          for (i = 0; i < pfd->num_phase_funcs; i++) {
            if (fscanf(imp, "%lf", &(mat_ptr->mp2nd->FlowingLiquid_viscosity_phase[i])) != 1) {
              GOMA_EH(GOMA_ERROR, "error reading phase FlowingLiquid viscosity");
            }
            SPF(endofstring(es), " %g", mat_ptr->mp2nd->FlowingLiquid_viscosity_phase[i]);
          }
        }

      } else {
        GOMA_EH(GOMA_ERROR,
                "Second Level Set FlowingLiquid Viscosity model can only be CONSTANT.\n");
      }
    }

    ECHO(es, echo_file);

    strcpy(search_string, "Inertia Coefficient");
    model_read = look_for_mat_prop(imp, search_string, &(mat_ptr->InertiaCoefficientModel),
                                   &(mat_ptr->Inertia_coefficient), NO_USER, NULL, model_name,
                                   SCALAR_INPUT, &NO_SPECIES, es);
    ECHO(es, echo_file);
  } else {
    mat_ptr->FlowingLiquid_viscosity = 1.0;
    mat_ptr->Inertia_coefficient = 0.0;
    SPF(es, "\t(%s = %s %.4g)", "FlowingLiquid Viscosity", "CONSTANT",
        mat_ptr->FlowingLiquid_viscosity);
    ECHO(es, echo_file);
    SPF(es, "\t(%s = %s %.4g)", "Inertia Coefficient", "CONSTANT", mat_ptr->Inertia_coefficient);
    ECHO(es, echo_file);
  }

  if ((mat_ptr->PorousMediaType == POROUS_UNSATURATED ||
       mat_ptr->PorousMediaType == POROUS_TWO_PHASE ||
       mat_ptr->PorousMediaType == POROUS_SHELL_UNSATURATED) &&
      (!PorousShellOn)) {
    model_read = look_for_mat_prop(imp, "Capillary Network Stress", &(mat_ptr->CapStress), &(a0),
                                   NO_USER, NULL, model_name, NO_INPUT, &NO_SPECIES, es);
    if (model_read == -1 && !strcmp(model_name, "WETTING")) {
      mat_ptr->CapStress = WETTING;
    } else if (model_read == -1 && !strcmp(model_name, "PARTIALLY_WETTING")) {
      mat_ptr->CapStress = PARTIALLY_WETTING;
    } else if (model_read == -1 && !strcmp(model_name, "COMPRESSIBLE")) {
      mat_ptr->CapStress = COMPRESSIBLE;
    } else if (model_read == -1 && !strcmp(model_name, "NONE")) {
      mat_ptr->CapStress = NO_CAP_STRESS;
    } else
      GOMA_EH(GOMA_ERROR, "Error getting Capillary Network Stress");

    ECHO(es, echo_file);

    model_read = look_for_mat_prop(imp, "Rel Gas Permeability", &(mat_ptr->RelGasPermModel),
                                   &(mat_ptr->rel_gas_perm), NO_USER, NULL, model_name,
                                   SCALAR_INPUT, &NO_SPECIES, es);
    if (model_read == -1 && !strcmp(model_name, "SUM_TO_ONE")) {
      mat_ptr->RelGasPermModel = SUM_TO_ONE;
      num_const = read_constants(imp, &(mat_ptr->u_rel_gas_perm), NO_SPECIES);
      if (num_const < 1) {
        sr = sprintf(err_msg, "Matl %s expected at least 1 constant for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Rel Gas Permeability", "SUM_TO_ONE");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_rel_gas_perm = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_rel_gas_perm);
    } else {
      GOMA_EH(model_read, "Rel Gas Permeability");
    }

    ECHO(es, echo_file);

    model_read = look_for_mat_prop(imp, "Rel Liq Permeability", &(mat_ptr->RelLiqPermModel),
                                   &(mat_ptr->rel_liq_perm), NO_USER, NULL, model_name,
                                   SCALAR_INPUT, &NO_SPECIES, es);
    if (model_read == -1 && !strcmp(model_name, "VAN_GENUCHTEN")) {
      mat_ptr->RelLiqPermModel = VAN_GENUCHTEN;
      num_const = read_constants(imp, &(mat_ptr->u_rel_liq_perm), NO_SPECIES);
      if (num_const < 4) {
        sr = sprintf(err_msg, "Matl %s expected at least 4 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Rel Liq Permeability", "VAN_GENUCHTEN");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_rel_liq_perm = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_rel_liq_perm);
    } else if (model_read == -1 && !strcmp(model_name, "PSD_VOL")) {
      mat_ptr->RelLiqPermModel = PSD_VOL;
      num_const = read_constants(imp, &(mat_ptr->u_rel_liq_perm), NO_SPECIES);
      if (num_const < 1) {
        sr = sprintf(err_msg, "Matl %s expected at least 1 constant for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Rel Liq Permeability", "PSD_VOL");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_rel_liq_perm = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_rel_liq_perm);
    } else if (model_read == -1 && !strcmp(model_name, "PSD_WEXP")) {
      mat_ptr->RelLiqPermModel = PSD_WEXP;
      num_const = read_constants(imp, &(mat_ptr->u_rel_liq_perm), NO_SPECIES);
      if (num_const < 1) {
        sr = sprintf(err_msg, "Matl %s expected at least 1 constant for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Rel Liq Permeability", "PSD_WEXP");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_rel_liq_perm = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_rel_liq_perm);
    } else if (model_read == -1 && !strcmp(model_name, "PSD_SEXP")) {
      mat_ptr->RelLiqPermModel = PSD_SEXP;
      num_const = read_constants(imp, &(mat_ptr->u_rel_liq_perm), NO_SPECIES);
      if (num_const < 1) {
        sr = sprintf(err_msg, "Matl %s expected at least 1 constant for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Rel Liq Permeability", "PSD_SEXP");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_rel_liq_perm = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_rel_liq_perm);
    } else if (model_read == -1 && !strcmp(model_name, "EXTERNAL_FIELD")) {

      if (fscanf(imp, "%s", input) != 1) {
        GOMA_EH(GOMA_ERROR,
                "Expecting trailing keyword for Rel Liq Permeability EXTERNAL_FIELD model.\n");
      }
      ii = 0;
      for (j = 0; j < efv->Num_external_field; j++) {
        if (!strcmp(efv->name[j], input)) {
          ii = 1;
          mat_ptr->rel_liq_perm_external_field_index = j;
        }
      }
      if (ii == 0) {
        GOMA_EH(GOMA_ERROR, "Must activate external fields to use this Rel Liq Permeability model "
                            "and requested name must be an external field name");
      }
      mat_ptr->RelLiqPermModel = EXTERNAL_FIELD;

      /* pick up scale factor for property */
      num_const = read_constants(imp, &(mat_ptr->u_rel_liq_perm), NO_SPECIES);
      mat_ptr->len_u_rel_liq_perm = num_const;
      if (num_const < 1) {
        sr = sprintf(err_msg, "Matl %s expected at least 1 constant for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Rel Liq Permeability", "EXTERNAL_FIELD");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
    } else if (model_read == -1 && !strcmp(model_name, "VAN_GENUCHTEN_EXTERNAL")) {
      if (fscanf(imp, "%s", input) != 1) {
        GOMA_EH(GOMA_ERROR, "Expecting trailing keyword for VAN_GENUCHTEN_EXTERNAL_FIELD model.\n");
      }
      ii = 0;
      for (j = 0; j < efv->Num_external_field; j++) {
        if (!strcmp(efv->name[j], input)) {
          ii = 1;
          mat_ptr->rel_liq_perm_external_field_index = j;
        }
      }
      if (ii == 0) {
        GOMA_EH(GOMA_ERROR,
                "Must activate external fields to use this VAN_GENUCHTEN_EXTERNAL model");
      }
      mat_ptr->RelLiqPermModel = VAN_GENUCHTEN_EXTERNAL;
      num_const = read_constants(imp, &(mat_ptr->u_rel_liq_perm), NO_SPECIES);
      if (num_const < 5) {
        sr = sprintf(err_msg, "Matl %s expected at least 5 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Rel Liq Permeability", "VAN_GENUCHTEN_EXTERNAL");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_rel_liq_perm = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_rel_liq_perm);
    } else {
      GOMA_EH(model_read, "Rel Liq Permeability");
    }
    ECHO(es, echo_file);

    /* Read in Saturation */ model_read = look_for_mat_proptable(
        imp, "Saturation", &(mat_ptr->SaturationModel), &(mat_ptr->saturation),
        &(mat_ptr->u_saturation), &(mat_ptr->len_u_saturation), &(mat_ptr->saturation_tableid),
        model_name, SCALAR_INPUT, &NO_SPECIES, es);

    if (model_read == -1 && !strcmp(model_name, "VAN_GENUCHTEN")) {
      mat_ptr->SaturationModel = VAN_GENUCHTEN;
      num_const = read_constants(imp, &(mat_ptr->u_saturation), NO_SPECIES);
      if (num_const < 4) {
        sr = sprintf(err_msg, "Matl %s expected at least 4 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Saturation", "VAN_GENUCHTEN");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_saturation = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_saturation);
    } else if (model_read == -1 && !strcmp(model_name, "PSD_VOL")) {
      mat_ptr->SaturationModel = PSD_VOL;
      num_const = read_constants(imp, &(mat_ptr->u_saturation), NO_SPECIES);
      if (num_const < 4) {
        sr = sprintf(err_msg, "Matl %s expected at least 4 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Saturation", "PSD_VOL");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_saturation = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_saturation);
    } else if (model_read == -1 && !strcmp(model_name, "PSD_WEXP")) {
      mat_ptr->SaturationModel = PSD_WEXP;
      num_const = read_constants(imp, &(mat_ptr->u_saturation), NO_SPECIES);
      if (num_const < 4) {
        sr = sprintf(err_msg, "Matl %s expected at least 4 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Saturation", "PSD_WEXP");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_saturation = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_saturation);
    } else if (model_read == -1 && !strcmp(model_name, "PSD_SEXP")) {
      mat_ptr->SaturationModel = PSD_SEXP;
      num_const = read_constants(imp, &(mat_ptr->u_saturation), NO_SPECIES);
      if (num_const < 4) {
        sr = sprintf(err_msg, "Matl %s expected at least 4 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Saturation", "PSD_SEXP");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_saturation = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_saturation);
    } else if (model_read == -1 && !strcmp(model_name, "TANH")) {
      mat_ptr->SaturationModel = TANH;
      num_const = read_constants(imp, &(mat_ptr->u_saturation), NO_SPECIES);
      if (num_const < 4) {
        sr = sprintf(err_msg, "Matl %s expected at least 4 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Saturation", "TANH");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_saturation = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_saturation);
    } else if (model_read == -1 && !strcmp(model_name, "SHELL_TANH")) {
      mat_ptr->SaturationModel = SHELL_TANH;
      num_const = read_constants(imp, &(mat_ptr->u_saturation), NO_SPECIES);
      if (num_const < 4) {
        sr = sprintf(err_msg, "Matl %s expected at least 4 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Saturation", "SHELL_TANH");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_saturation = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_saturation);
    } else if (model_read == -1 && !strcmp(model_name, "TANH_EXTERNAL")) {
      if (fscanf(imp, "%s", input) != 1) {
        GOMA_EH(GOMA_ERROR, "Expecting trailing keyword for TANH EXTERNAL_FIELD model.\n");
      }
      ii = 0;
      for (j = 0; j < efv->Num_external_field; j++) {
        if (!strcmp(efv->name[j], input)) {
          ii = 1;
          mat_ptr->SAT_external_field_index = j;
        }
      }
      if (ii == 0) {
        GOMA_EH(GOMA_ERROR, "Must activate external fields to use this TANH_EXTERNAL model");
      }
      mat_ptr->SaturationModel = TANH_EXTERNAL;
      num_const = read_constants(imp, &(mat_ptr->u_saturation), NO_SPECIES);
      if (num_const < 8) {
        sr = sprintf(err_msg, "Matl %s expected at least 8 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Saturation", "TANH_EXTERNAL");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_saturation = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_saturation);
    } else if (model_read == -1 && !strcmp(model_name, "TANH_HYST")) {
      mat_ptr->SaturationModel = TANH_HYST;
      num_const = read_constants(imp, &(mat_ptr->u_saturation), NO_SPECIES);
      if (num_const < 10) {
        sr = sprintf(err_msg, "Matl %s expected at least 10 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Saturation", "TANH_HYST");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_saturation = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_saturation);

    } else {
      GOMA_EH(model_read, "Saturation");
    }
    ECHO(es, echo_file);

    if (have_por_energy == 1) {
      /*
       * Density of solid matrix
       */
      mat_ptr->matrix_density = 1.;
      model_read = look_for_mat_prop(imp, "Matrix Density", &(mat_ptr->PorousMatrixDensityModel),
                                     &(mat_ptr->matrix_density), NO_USER, NULL, model_name,
                                     SCALAR_INPUT, &NO_SPECIES, es);
      if (model_read == -1 && strcmp(model_name, "CONSTANT")) {
        GOMA_EH(GOMA_ERROR, " Only a CONSTANT Matrix Density is allowed.");
      } else {
        GOMA_EH(model_read, "Matrix Density card needed for R_POR_ENERGY ");
      }

      ECHO(es, echo_file);

      /*
       * Specific heat of solid matrix
       */
      mat_ptr->specific_heat = 1.;
      model_read = look_for_mat_prop(
          imp, "Matrix Specific Heat", &(mat_ptr->PorousSpecificHeatModel),
          &(mat_ptr->matrix_density), NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);

      if (model_read == -1 && strcmp(model_name, "CONSTANT")) {
        GOMA_EH(GOMA_ERROR, " Only a CONSTANT Matrix Specific Heat is allowed.");
      } else {
        GOMA_EH(model_read, "Matrix Specific Heat needed for R_POR_ENERGY");
      }
      ECHO(es, echo_file);
    }
  }
  /*
   * Porous Section 2: Special Numerical treatment inputs are made here
   * for each porous equation.  Recall that each equation is really a species
   * conservation law, and so this will really be over 1 to the number of
   * possible species, in whatever phase, someday.  The input here is specifically
   * for the liquid phase solvent component, the gas-phase solvent component (air),
   * and the solid-phase solvent compont (just the solid skeleton).   Eventually you
   * will want the liquid-phase components to be looped when we add multicomponent
   * porous liquid capability. Possibly also you will loop over base gas phase
   * components (which will usually be air, with all other components in the gas
   * coming from the liquid
   */

  /* TAB here. just like to comment in passing that this whole porous media input section is
   * an unpleasant mess and needs fixing.  November 2005 */

  if (mat_ptr->PorousMediaType != CONTINUOUS && mat_ptr->PorousMediaType != POROUS_BRINKMAN &&
      (!PorousShellOn)) {

    /*
     *   Porous Weight Function:
     *
     *    This is where you specify GALERKIN or SUPG weighting functions
     *    To follow the theory of Brooks and Hughes, the weighting function
     *    value should be set equal to 1.0
     */
    model_read = look_for_mat_prop(imp, "Porous Weight Function", &(mat_ptr->Porous_wt_funcModel),
                                   &(mat_ptr->Porous_wt_func), NO_USER, NULL, model_name,
                                   SCALAR_INPUT, &NO_SPECIES, es);

    if (!strcmp(model_name, "GALERKIN")) {
      mat_ptr->Porous_wt_funcModel = GALERKIN;
      mat_ptr->Porous_wt_func = 0.0;
    } else if (!strcmp(model_name, "SUPG")) {
      int err;
      mat_ptr->Porous_wt_funcModel = SUPG;
      err = fscanf(imp, "%lg", &(mat_ptr->Porous_wt_func));
      if (err != 1) {
        GOMA_EH(GOMA_ERROR, "Expected to read one double for Porous Weight Function SUPG");
      }
      SPF(endofstring(es), " %.4g", mat_ptr->Porous_wt_func);
    } else {
      mat_ptr->Porous_wt_funcModel = GALERKIN;
      mat_ptr->Porous_wt_func = 0.0;

      SPF(es, "\t(%s = %s)", "Porous Weight Function", "GALERKIN");
    }
    ECHO(es, echo_file);

    /*
     *   Porous Mass Lumping:
     *
     *    This is where you specify whether you want to use
     *    Mass Lumping of the time derivative or whether you
     *    want to use a consistent time derivative treatment.
     *
     *    Note, the combination of Mass Lumping and SUPG has
     *    been shown to exhibit monotone behavior. Also note that
     *    mass lumping is only useful for partially saturated flows.
     */

    mat_ptr->Porous_Mass_Lump = FALSE;
    model_read = look_for_mat_prop(imp, "Porous Mass Lumping", &(mat_ptr->Porous_wt_funcModel),
                                   &(mat_ptr->Porous_wt_func), NO_USER, NULL, model_name,
                                   SCALAR_INPUT, &NO_SPECIES, es);

    if (mat_ptr->PorousMediaType != POROUS_SATURATED &&
        mat_ptr->PorousMediaType != POROUS_SHELL_UNSATURATED) {
      if (!strcasecmp(model_name, "yes") || !strcasecmp(model_name, "true")) {
        mat_ptr->Porous_Mass_Lump = TRUE;
      } else if (!strcasecmp(model_name, "no") || !strcasecmp(model_name, "false")) {
        mat_ptr->Porous_Mass_Lump = FALSE;
      } else {
        GOMA_EH(GOMA_ERROR, "Porous Mass Lumping must be set to TRUE, YES, FALSE, or NO!");
        mat_ptr->Porous_Mass_Lump = FALSE;
      }
    }
    /* Jim Simmons: make sure this consistency check is added to new parse code */
    if (mat_ptr->Porous_Mass_Lump == TRUE && mat_ptr->PorousMediaType == POROUS_SATURATED) {
      GOMA_EH(GOMA_ERROR,
              "Mass Lumping is NOT implemented for POROUS_SATURATED problems. Turn to FALSE");
    }

    ECHO(es, echo_file);

    /*
     * set porous number equal to max number of porous media phases
     * it is changed to the porous phase number of input property by
     * look_for_mat_prop
     */
    porous_no = mat_ptr->Num_Porous_Eqn;

    /*
     *   Porous Time Integration
     */
    model_read = look_for_mat_prop(imp, "Porous Time Integration", &(PorousTimeIntegration), &(a0),
                                   NO_USER, NULL, model_name, NO_INPUT, &porous_no, es);
    if (model_read == -1 && !strcmp(model_name, "STANDARD")) {
      mat_ptr->PorousTimeIntegration[porous_no] = STANDARD;
    } else if (model_read == -1 && !strcmp(model_name, "TAYLOR_GALERKIN")) {
      mat_ptr->PorousTimeIntegration[porous_no] = TAYLOR_GALERKIN;
    } else if (model_read == -1 && !strcmp(model_name, "TAYLOR_GALERKIN_EXP")) {
      mat_ptr->PorousTimeIntegration[porous_no] = TAYLOR_GALERKIN_EXP;
    } else {
      mat_ptr->PorousTimeIntegration[porous_no] = STANDARD;
      SPF(es, "\t(%s = %s)", "Porous Time Integration", "STANDARD");
    }
    ECHO(es, echo_file);

    /*
     * Porous Section 3:  Species Dependent properties are finally input
     * here.  Remember that really each porous media equation is a species balance
     * across all phases. Hence you will notice we split out the phases where
     * appropriate in separate card inputs, and loop someday over the species in
     * each base phase.     This is kind of confusing but eventually after specifying
     * the phase diffusion transport constitutive equation types, we then loop
     * over the base solvents of each phase.  In the case right now, there is only
     * one liquid solvent, but it can exist in 2 phases, hence the Liquid and Gas
     * qualifiers.   In the Gas, there is only one base component and it is insoluble
     * in liquid (e.g. air).  The finaly section here deals with specific properties
     * of air.
     */

    if (pd_glob[mn]->Num_Porous_Eqn > 0) {

      /*
       * Hey! If it's in there, then it should be fine. To avoid the ordering
       * complaints, let's rewind the mat file and start looking from the
       * beginning.
       */
      rewind(imp);

      /* Important Note from PRS: When we start allowing multicomponent (nonideal)
       * diffusion in the liquid phase and in the gas phase, we will need to make this
       * two cards (one for liquid and one for gas).  I could think of wanting to do
       * Stephan-Maxwell in the gas, and fickian only in the liquid, with, of course,
       * the microstructure tortuosity corrections. For now leave this be.  05/29/01
       */

      model_read = look_for_mat_prop(imp, "Porous Diffusion Constitutive Equation",
                                     &(PorousDiffusionConstitutiveEquation), &(a0), NO_USER, NULL,
                                     model_name, NO_INPUT, &NO_SPECIES, es);
      if (!strcmp(model_name, "DARCY")) {
        GOMA_EH(GOMA_ERROR,
                "DARCY model for Porous Diffusion Constitutive Equation no longer allowed");
      } else if (!strcmp(model_name, "DARCY_FICKIAN")) {
        PorousDiffusionConstitutiveEquation = DARCY_FICKIAN;
      } else if (!strcmp(model_name, "POWERLAW_DARCY_FICKIAN")) {
        PorousDiffusionConstitutiveEquation = POWERLAW_DARCY_FICKIAN;
      } else {
        sprintf(err_msg, "Matl %s, \"%s\" = \"%s\"?\n(valid: %s)\n[%s]", pd_glob[mn]->MaterialName,
                "Porous Diffusion Constitutive Equation", model_name, " DARCY_FICKIAN only",
                "If still bad - check orderings in the mat file!");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      pd_glob[mn]->PorousFluxModel = PorousDiffusionConstitutiveEquation;
      cr_glob[mn]->PorousFluxModel = PorousDiffusionConstitutiveEquation;
      ECHO(es, echo_file);
    }

    /*
     ***************************************************************************
     *  LOOP OVER THE SPECIES THAT can be present in gas and liquid phases.
     * These come from the base liquid composition as we are assuming the
     *  base gas phase components,  e.g. air, are insoluble in the liquid phase.
     * As of 05/29/01 we only  allow single component
     * and so this is a loop of 1 to 1.   Change will be made here for multicomponent.
     **************************************************************************
     */

    ECHO("\n----Species Transport Properties \n", echo_file);

    Num_liquid_phase_components = 1; /* Put this in the porous structure and read in
                                      * when you expand to more than one component
                                      */
    for (j = 0; j < Num_liquid_phase_components; j++) {
      /* N.B. porous_no is the species number, really.  Should be zero in
       * all cases, until we go multicomponent.  Name should be changed.
       * For now this card will be the vapor diffusivity value of the primary
       * liquid solvent component.    It is also assumed to be binary
       * diffusion.  However, this should also have a Porous Liquid Diffusivity
       * counterpart for the same species to account for the diffusion model
       * in the liquid phase.
       * PRS 052901
       */

      if (mat_ptr->PorousMediaType != POROUS_SATURATED) {
        model_read = look_for_porous_prop(
            imp, "Porous Gas Diffusivity", mat_ptr, mat_ptr->PorousDiffusivityModel,
            mat_ptr->porous_diffusivity, mat_ptr->u_porous_diffusivity,
            mat_ptr->len_u_porous_diffusivity, model_name, 0, &porous_no, es);
        /*
         * Postprocess unique Porous Diffusivity models
         */
        if (model_read == -1) {
          GOMA_EH(
              GOMA_ERROR,
              "Porous Gas Diffusivity: Bad Card syntax or need another set of porous mat cards?");
        } else if (model_read == 0) {
          if (!strcmp(model_name, "POROUS")) {
            mat_ptr->PorousDiffusivityModel[porous_no] = POROUS;
            num_const = mat_ptr->len_u_porous_diffusivity[porous_no];
            if (num_const < 5) {
              sprintf(err_msg, "Matl %s (eqn %d) needs at least 5 constants for %s %s model.\n",
                      pd_glob[mn]->MaterialName, porous_no, "Porous Gas Diffusivity", "POROUS");
              GOMA_EH(GOMA_ERROR, err_msg);
            }
          }

        } /* End of if (model_read == 0) */

        ECHO(es, echo_file);

        /*
         *  Porous Latent Heat of Vaporization Section
         *
         * set porous number equal to max number of porous media phases
         * it is changed to the porous phase number of input properly by
         * look_for_mat_prop
         */
        porous_no = mat_ptr->Num_Porous_Eqn;
        model_read =
            look_for_mat_prop(imp, "Porous Latent Heat Vaporization",
                              mat_ptr->PorousLatentHeatVapModel, mat_ptr->porous_latent_heat_vap,
                              NO_USER, NULL, model_name, SCALAR_INPUT, &porous_no, es);
        GOMA_EH(model_read, "Porous Latent Heat Vaporization");
        ECHO(es, echo_file);

        /*
         *  Porous Latent Heat of Fusion Section
         *
         * set porous number equal to max number of porous media phases
         * it is changed to the porous phase number of input property by
         * look_for_mat_prop
         */
        porous_no = mat_ptr->Num_Porous_Eqn;
        model_read = look_for_mat_prop(imp, "Porous Latent Heat Fusion",
                                       mat_ptr->PorousLatentHeatFusionModel,
                                       mat_ptr->porous_latent_heat_fusion, NO_USER, NULL,
                                       model_name, SCALAR_INPUT, &porous_no, es);
        GOMA_EH(model_read, "Porous Latent Heat Fusion");
        ECHO(es, echo_file);

        /*
         *  Vapor Pressure Section
         *
         * set porous number equal to max number of porous media phases
         * it is changed to the porous phase number of input property by
         * look_for_mat_prop
         */
        porous_no = mat_ptr->Num_Porous_Eqn;
        model_read =
            look_for_mat_prop(imp, "Porous Vapor Pressure", mat_ptr->PorousVaporPressureModel,
                              mat_ptr->porous_vapor_pressure, NO_USER, NULL, model_name,
                              SCALAR_INPUT, &porous_no, es);

        if (model_read == -1 && (!strcmp(model_name, "KELVIN") || !strcmp(model_name, "FLAT"))) {
          if (!strcmp(model_name, "KELVIN"))
            mat_ptr->PorousVaporPressureModel[porous_no] = KELVIN;
          if (!strcmp(model_name, "FLAT"))
            mat_ptr->PorousVaporPressureModel[porous_no] = FLAT;

          if (mat_ptr->PorousMediaType == POROUS_TWO_PHASE) {
            num_const = read_constants(imp, mat_ptr->u_porous_vapor_pressure, porous_no);
            if (num_const < 5) {
              sprintf(err_msg, "Matl %s (%s, conc %d) needs 5 constants for %s %s model.\n",
                      pd_glob[mn]->MaterialName, "porous 2-phase", porous_no,
                      "Porous Vapor Pressure", "KELVIN or FLAT");
              GOMA_EH(GOMA_ERROR, err_msg);
            }
            mat_ptr->len_u_porous_vapor_pressure[porous_no] = num_const;
            SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_porous_vapor_pressure[porous_no]);
          } else if (mat_ptr->PorousMediaType == POROUS_UNSATURATED) {
            num_const = read_constants(imp, mat_ptr->u_porous_vapor_pressure, porous_no);
            if (num_const < 5) {
              sprintf(err_msg, "Matl %s (%s, conc %d) needs 5 constants for %s %s model.\n",
                      pd_glob[mn]->MaterialName, "porous unsaturated", porous_no,
                      "Porous Vapor Pressure", "KELVIN OR FLAT");
              GOMA_EH(GOMA_ERROR, err_msg);
            }
            mat_ptr->len_u_porous_vapor_pressure[porous_no] = num_const;
            SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_porous_vapor_pressure[porous_no]);
          } else {
            sprintf(err_msg, "%s model invalid in matl %s unless %s is \"%s\" or \"%s\"\n",
                    "KELVIN or FLAT", pd_glob[mn]->MaterialName, "Media Type", "POROUS_TWO_PHASE",
                    "POROUS_UNSATURATED");
            GOMA_EH(GOMA_ERROR, err_msg);
          }
        } else if (model_read == -1 && !strcmp(model_name, "IDEAL_GAS")) {
          GOMA_EH(GOMA_ERROR, "No longer set the IDEAL_GAS properties/model on the "
                              "vapor pressure card. Use Gas Constants");
        } else if (model_read == -1 && !strcmp(model_name, "NON_VOLATILE")) {
          mat_ptr->PorousVaporPressureModel[porous_no] = NON_VOLATILE;
        } else {
          GOMA_EH(model_read, "Porous Vapor Pressure");
        }
        ECHO(es, echo_file);
      } /*end if(!POROUS_SATURATED) */

      /*
       *  Porous Liquid Volume Expansion Section
       *
       * set porous number equal to max number of porous media phases
       * it is changed to the porous phase number of input property  by
       *  look_for_mat_prop
       *
       * PRS (052901): This is a phase property that is component dependent.  You need to
       * add the Porous Gas Volume Expansion card for the gas phase as well, when
       * you go multicomponent.
       */
      if (mat_ptr->PorousMediaType != POROUS_SHELL_UNSATURATED) {
        model_read = look_for_porous_prop(imp, "Porous Liquid Volume Expansion", mat_ptr,
                                          mat_ptr->PorVolExpModel, mat_ptr->porous_vol_expansion,
                                          NO_USER, NULL, model_name, SCALAR_INPUT, &porous_no, es);
        GOMA_EH(model_read, "Porous Liquid Volume Expansion");
        ECHO(es, echo_file);
      }

    } /* End for(j=0;j<Num_liquid_phase_components; j++) */

    /*Final porous input stage */

    Num_insoluble_gas_phase_components = 1; /* Put this in the porous structure and read in
                                             * when you expand to more than one component
                                             */
    /*These cards are needed if there in all gas/liquid flow cases,
     * viz. TWO_PHASE and UNSATURATED.  In the Unsaturated case this
     * card is used to set the external CONSTANT ambient pressure*/
    if (mat_ptr->PorousMediaType != CONTINUOUS && mat_ptr->PorousMediaType != POROUS_SATURATED) {
      for (j = 0; j < Num_insoluble_gas_phase_components; j++) {
        /*
         *  Insoluble Porous Gas phase component  Section
         *  This supplants old way of using the Vapor Pressure
         *  card for this. YOu may want to make this 3 separate
         *  cards for the ideal gas model (viz. MW, Universal Gas Const, T)
         *
         */
        model_read =
            look_for_mat_prop(imp, "Porous Gas Constants", &(mat_ptr->PorousGasConstantsModel),
                              &(mat_ptr->porous_gas_constants), NO_USER, NULL, model_name,
                              SCALAR_INPUT, &NO_SPECIES, es);

        if (model_read == -1 && !strcmp(model_name, "IDEAL_GAS")) {
          mat_ptr->PorousGasConstantsModel = IDEAL_GAS;

          num_const = read_constants(imp, &(mat_ptr->u_porous_gas_constants), NO_SPECIES);
          if (num_const < 4) {
            sr = sprintf(err_msg, "Matl %s (conc %d) needs 4 constants for %s %s model.\n",
                         pd_glob[mn]->MaterialName, 0, "Porous Gas Constants", "IDEAL_GAS");
            GOMA_EH(GOMA_ERROR, err_msg);
          }
          mat_ptr->len_u_porous_gas_constants = num_const;
          SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_porous_gas_constants);
        } else if (!strcmp(model_name, "CONSTANT")) {
          GOMA_EH(GOMA_ERROR, "Ironically we don't allow a CONSTANT model for Porous Gas "
                              "Constants.  Try IDEAL_GAS");
        } else {
          GOMA_EH(model_read, "Porous Gas Constants");
        }
        ECHO(es, echo_file);
      }
    }
  } /*if ( mat_ptr->PorousMediaType != CONTINUOUS && != POROUS_BRINKMAN)*/

  /*
   *  Porous Shell Section
   *
   */

  if ((mat_ptr->PorousMediaType == POROUS_SHELL_UNSATURATED) &&
      (PorousShellOn)) /* This only applies to porous shell  equations */
  {

    /* create a list of the porous media equations
       active in this material                     */
    i = 0;

    for (int imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
      if (pd_glob[mn]->e[imtrx][R_SHELL_SAT_1]) {
        mat_ptr->Porous_Shell_Eqn[i] = R_SHELL_SAT_1;
        i++;
      }

      if (pd_glob[mn]->e[imtrx][R_SHELL_SAT_2]) {
        mat_ptr->Porous_Shell_Eqn[i] = R_SHELL_SAT_2;
        i++;
      }

      if (pd_glob[mn]->e[imtrx][R_SHELL_SAT_3]) {
        mat_ptr->Porous_Shell_Eqn[i] = R_SHELL_SAT_3;
        i++;
      }
    }

    if (i != pd_glob[mn]->Num_Porous_Shell_Eqn) {
      GOMA_EH(GOMA_ERROR, "Number of porous shell equations do not add up");
    }

    /***************************************************************************
     *  LOOP OVER THE POROUS SHELL LAYERS.
     **************************************************************************/

    /* Initialize porous shell no with maximum number of porous shell layers read from pd */

    for (ipore = 0; ipore < pd_glob[mn]->Num_Porous_Shell_Eqn; ipore++) {

      /******************
       *                *
       * Read porosity  *
       *                *
       *******************/
      porous_shell_no = pd_glob[mn]->Num_Porous_Shell_Eqn;
      strcpy(search_string, "Porosity");
      model_read = look_for_mat_prop(imp, search_string, mat_ptr->PorousShellPorosityModel,
                                     mat_ptr->PorousShellPorosity, NO_USER, NULL, model_name,
                                     SCALAR_INPUT, &porous_shell_no, es);

      if (model_read == -1 && !strcmp(model_name, "EXTERNAL_FIELD")) {
        if (fscanf(imp, "%s", input) != 1) {
          GOMA_EH(GOMA_ERROR, "Expecting trailing keyword for Porosity EXTERNAL_FIELD model.\n");
        }
        ii = 0;
        for (j = 0; j < efv->Num_external_field; j++) {
          if (!strcmp(efv->name[j], input)) {
            ii = 1;
            mat_ptr->por_shell_porosity_ext_field_index[ipore] = j;
          }
        }
        if (ii == 0) {
          GOMA_EH(GOMA_ERROR, "Cannot match the name with that in the external field file");
        }
        mat_ptr->PorousShellPorosityModel[ipore] = EXTERNAL_FIELD;

        /* pick up scale factor for property */
        num_const = read_constants(imp, mat_ptr->u_PorousShellPorosity, porous_shell_no);
        mat_ptr->len_u_PorousShellPorosity[ipore] = num_const;
        if (num_const < 1) {
          sr = sprintf(err_msg, "Matl %s expected at least 1 constant for %s %s model.\n",
                       pd_glob[mn]->MaterialName, "Porosity", "EXTERNAL_FIELD");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
      } else {
        GOMA_EH(model_read, "Porosity: Is card missing?");
      }

      if (porous_shell_no != ipore) {
        GOMA_EH(GOMA_ERROR, "Incomplete number of Porosity card");
      }

      ECHO(es, echo_file);

      /***************
       *             *
       * Read height *
       *             *
       ****************/
      porous_shell_no = pd_glob[mn]->Num_Porous_Shell_Eqn;
      strcpy(search_string, "Porous Shell Height");
      model_read = look_for_mat_prop(imp, search_string, mat_ptr->PorousShellHeightModel,
                                     mat_ptr->PorousShellHeight, NO_USER, NULL, model_name,
                                     SCALAR_INPUT, &porous_shell_no, es);
      if (model_read == -1 && !strcmp(model_name, "EXTERNAL_FIELD")) {
        if (fscanf(imp, "%s", input) != 1) {
          GOMA_EH(GOMA_ERROR,
                  "Expecting trailing keyword for Porous Shell Height EXTERNAL_FIELD model.\n");
        }
        ii = 0;
        for (j = 0; j < efv->Num_external_field; j++) {
          if (!strcmp(efv->name[j], input)) {
            ii = 1;
            mat_ptr->por_shell_height_ext_field_index[ipore] = j;
          }
        }
        if (ii == 0) {
          GOMA_EH(GOMA_ERROR, "Cannot match the name with that in the external field file");
        }
        mat_ptr->PorousShellHeightModel[ipore] = EXTERNAL_FIELD;

        /* pick up scale factor for property */
        num_const = read_constants(imp, mat_ptr->u_PorousShellHeight, porous_shell_no);
        mat_ptr->len_u_PorousShellHeight[ipore] = num_const;
        if (num_const < 1) {
          sr = sprintf(err_msg, "Matl %s expected at least 1 constant for %s %s model.\n",
                       pd_glob[mn]->MaterialName, "Porous Shell Height", "EXTERNAL_FIELD");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
      } else {
        GOMA_EH(model_read, "Porous Shell Height: Is card missing?");
      }

      if (porous_shell_no != ipore) {
        GOMA_EH(GOMA_ERROR, "Incomplete number of Porous Shell Height card");
      }
      ECHO(es, echo_file);

      /*********************
       *                   *
       * Read permeability *
       *                   *
       *********************/
      porous_shell_no = pd_glob[mn]->Num_Porous_Shell_Eqn;
      strcpy(search_string, "Permeability");
      model_read = look_for_mat_prop(imp, search_string, mat_ptr->PorousShellPermeabilityModel,
                                     mat_ptr->PorousShellPermeability, NO_USER, NULL, model_name,
                                     SCALAR_INPUT, &porous_shell_no, es);
      if (model_read == -1 && !strcmp(model_name, "EXTERNAL_FIELD")) {
        if (fscanf(imp, "%s", input) != 1) {
          GOMA_EH(GOMA_ERROR,
                  "Expecting trailing keyword for Permeability EXTERNAL_FIELD model.\n");
        }
        ii = 0;
        for (j = 0; j < efv->Num_external_field; j++) {
          if (!strcmp(efv->name[j], input)) {
            ii = 1;
            mat_ptr->por_shell_permeability_ext_field_index[ipore] = j;
          }
        }
        if (ii == 0) {
          GOMA_EH(GOMA_ERROR, "Cannot match the name with that in the external field file");
        }
        mat_ptr->PorousShellPermeabilityModel[ipore] = EXTERNAL_FIELD;

        /* pick up scale factor for property */
        num_const = read_constants(imp, mat_ptr->u_PorousShellPermeability, porous_shell_no);
        mat_ptr->len_u_PorousShellPermeability[ipore] = num_const;
        if (num_const < 1) {
          sr = sprintf(err_msg, "Matl %s expected at least 1 constant for %s %s model.\n",
                       pd_glob[mn]->MaterialName, "Permeability", "EXTERNAL_FIELD");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
      } else if (model_read == -1 && !strcmp(model_name, "ORTHOTROPIC")) {
        mat_ptr->PorousShellPermeabilityModel[ipore] = ORTHOTROPIC;
        num_const = read_constants(imp, mat_ptr->u_PorousShellPermeability, porous_shell_no);
        if (num_const < 12) {
          sr = sprintf(err_msg, "Matl %s expected at least 12 constants for %s %s model.\n",
                       pd_glob[mn]->MaterialName, "Permeability", "ORTHOTROPIC");
          GOMA_EH(GOMA_ERROR, err_msg);
        }

        memset(mat_ptr->PorousShellPermTensor[ipore], 0,
               sizeof(double) * DIM * DIM); /*these are loaded up later */
        mat_ptr->PorousShellPermeability[ipore] =
            mat_ptr->u_PorousShellPermeability[ipore][0]; /*just in case */

        mat_ptr->len_u_PorousShellPermeability[ipore] = num_const;
      } else {
        GOMA_EH(model_read, "Permeability: Is card missing?");
      }

      if (porous_shell_no != ipore) {
        GOMA_EH(GOMA_ERROR, "Incomplete number of Permeability card");
      }
      ECHO(es, echo_file);

      /******************************
       *                            *
       * Read relative permeability *
       *                            *
       * ****************************/
      porous_shell_no = pd_glob[mn]->Num_Porous_Shell_Eqn;
      strcpy(search_string, "Rel Liq Permeability");
      model_read = look_for_mat_prop(imp, search_string, mat_ptr->PorousShellRelPermModel,
                                     mat_ptr->PorousShellRelPerm, NO_USER, NULL, model_name,
                                     SCALAR_INPUT, &porous_shell_no, es);
      if (model_read == -1 && !strcmp(model_name, "EXTERNAL_FIELD")) {
        if (fscanf(imp, "%s", input) != 1) {
          GOMA_EH(GOMA_ERROR,
                  "Expecting trailing keyword for Rel Liq Permeability EXTERNAL_FIELD model.\n");
        }
        ii = 0;
        for (j = 0; j < efv->Num_external_field; j++) {
          if (!strcmp(efv->name[j], input)) {
            ii = 1;
            mat_ptr->por_shell_rel_perm_ext_field_index[ipore] = j;
          }
        }
        if (ii == 0) {
          GOMA_EH(GOMA_ERROR, "Cannot match the name with that in the external field file");
        }
        mat_ptr->PorousShellRelPermModel[ipore] = EXTERNAL_FIELD;
        /* pick up scale factor for property */
        num_const = read_constants(imp, mat_ptr->u_PorousShellRelPerm, porous_shell_no);
        mat_ptr->len_u_PorousShellRelPerm[ipore] = num_const;
        if (num_const < 1) {
          sr = sprintf(err_msg, "Matl %s expected at least 1 constant for %s %s model.\n",
                       pd_glob[mn]->MaterialName, "Permeability", "EXTERNAL_FIELD");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
      } else if (model_read == -1 && !strcmp(model_name, "VAN_GENUCHTEN")) {
        mat_ptr->PorousShellRelPermModel[ipore] = VAN_GENUCHTEN;
        num_const = read_constants(imp, mat_ptr->u_PorousShellRelPerm, porous_shell_no);
        if (num_const < 4) {
          sr = sprintf(err_msg, "Matl %s expected at least 4 constants for %s %s model.\n",
                       pd_glob[mn]->MaterialName, "Rel Liq Permeability", "VAN_GENUCHTEN");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
        mat_ptr->len_u_PorousShellRelPerm[ipore] = num_const;
      } else if (model_read == -1 && !strcmp(model_name, "VAN_GENUCHTEN_EXTERNAL")) {
        if (fscanf(imp, "%s", input) != 1) {
          GOMA_EH(GOMA_ERROR, "Expecting trailing keyword for Rel Liq Permeability "
                              "VAN_GENUCHTEN_EXTERNAL model.\n");
        }
        ii = 0;
        for (j = 0; j < efv->Num_external_field; j++) {
          if (!strcmp(efv->name[j], input)) {
            ii = 1;
            mat_ptr->por_shell_rel_perm_ext_field_index[ipore] = j;
          }
        }
        if (ii == 0) {
          GOMA_EH(GOMA_ERROR, "Cannot match the name with that in the external field file");
        }
        mat_ptr->PorousShellRelPermModel[ipore] = VAN_GENUCHTEN_EXTERNAL;
        num_const = read_constants(imp, mat_ptr->u_PorousShellRelPerm, porous_shell_no);
        if (num_const < 5) {
          sr = sprintf(err_msg, "Matl %s expected at least 5 constants for %s %s model.\n",
                       pd_glob[mn]->MaterialName, "Rel Liq Permeability", "VAN_GENUCHTEN_EXTERNAL");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
        mat_ptr->len_u_PorousShellRelPerm[ipore] = num_const;
      } else {
        GOMA_EH(model_read, "Rel Liq Permeability: Is card missing?");
      }

      if (porous_shell_no != ipore) {
        GOMA_EH(GOMA_ERROR, "Incomplete number of Rel Liq Permeability card");
      }
      ECHO(es, echo_file);

      /***************************
       *                         *
       * Read cross permeability *
       *                         *
       ***************************/
      porous_shell_no = pd_glob[mn]->Num_Porous_Shell_Eqn;
      strcpy(search_string, "Porous Shell Cross Permeability");
      model_read = look_for_mat_prop(imp, search_string, mat_ptr->PorousShellCrossPermeabilityModel,
                                     mat_ptr->PorousShellCrossPermeability, NO_USER, NULL,
                                     model_name, SCALAR_INPUT, &porous_shell_no, es);
      if (model_read == -1 && !strcmp(model_name, "EXTERNAL_FIELD")) {
        if (fscanf(imp, "%s", input) != 1) {
          GOMA_EH(GOMA_ERROR, "Expecting trailing keyword for Porous Shell Cross Permeability "
                              "EXTERNAL_FIELD model.\n");
        }
        ii = 0;
        for (j = 0; j < efv->Num_external_field; j++) {
          if (!strcmp(efv->name[j], input)) {
            ii = 1;
            mat_ptr->por_shell_cross_permeability_ext_field_index[ipore] = j;
          }
        }
        if (ii == 0) {
          GOMA_EH(GOMA_ERROR, "Cannot match the name with that in the external field file");
        }
        mat_ptr->PorousShellCrossPermeabilityModel[ipore] = EXTERNAL_FIELD;

        /* pick up scale factor for property */
        num_const = read_constants(imp, mat_ptr->u_PorousShellCrossPermeability, porous_shell_no);
        mat_ptr->len_u_PorousShellCrossPermeability[ipore] = num_const;
        if (num_const < 1) {
          sr = sprintf(err_msg, "Matl %s expected at least 1 constant for %s %s model.\n",
                       pd_glob[mn]->MaterialName, "Porous Shell Cross Permeability",
                       "EXTERNAL_FIELD");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
      } else {
        GOMA_EH(model_read, "Porous Shell Cross Permeability: Is card missing?");
      }

      if (porous_shell_no != ipore) {
        GOMA_EH(GOMA_ERROR, "Incomplete number of Porous Shell Cross Permeability card");
      }
      ECHO(es, echo_file);

      /***************************
       *                         *
       * Read capillary pressure *
       *                         *
       ***************************/
      porous_shell_no = pd_glob[mn]->Num_Porous_Shell_Eqn;
      strcpy(search_string, "Capillary Pressure");
      model_read = look_for_mat_prop(imp, search_string, mat_ptr->PorousShellCapPresModel,
                                     mat_ptr->PorousShellCapPres, NO_USER, NULL, model_name,
                                     SCALAR_INPUT, &porous_shell_no, es);

      if (model_read == -1 && !strcmp(model_name, "ATANH")) {
        mat_ptr->PorousShellCapPresModel[ipore] = ATANH;
        num_const = read_constants(imp, mat_ptr->u_PorousShellCapPres, porous_shell_no);
        if (num_const < 4) {
          sr = sprintf(err_msg, "Matl %s expected at least 4 constants for %s %s model.\n",
                       pd_glob[mn]->MaterialName, "Capillary Pressure", "ATANH");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
        mat_ptr->len_u_PorousShellCapPres[ipore] = num_const;
      } else if (model_read == -1 && !strcmp(model_name, "SINH")) {
        mat_ptr->PorousShellCapPresModel[ipore] = SINH;
        num_const = read_constants(imp, mat_ptr->u_PorousShellCapPres, porous_shell_no);
        if (num_const < 4) {
          sr = sprintf(err_msg, "Matl %s expected at least 4 constants for %s %s model.\n",
                       pd_glob[mn]->MaterialName, "Capillary Pressure", "SINH");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
        mat_ptr->len_u_PorousShellCapPres[ipore] = num_const;
      } else if (model_read == -1 && !strcmp(model_name, "VAN_GENUCHTEN")) {
        mat_ptr->PorousShellCapPresModel[ipore] = VAN_GENUCHTEN;
        num_const = read_constants(imp, mat_ptr->u_PorousShellCapPres, porous_shell_no);
        if (num_const < 4) {
          sr = sprintf(err_msg, "Matl %s expected at least 4 constants for %s %s model.\n",
                       pd_glob[mn]->MaterialName, "Capillary Pressure", "VAN_GENUCHTEN");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
        mat_ptr->len_u_PorousShellCapPres[ipore] = num_const;
      } else if (model_read == -1 && !strcmp(model_name, "VAN_GENUCHTEN_EXTERNAL")) {
        if (fscanf(imp, "%s", input) != 1) {
          GOMA_EH(
              GOMA_ERROR,
              "Expecting trailing keyword for Capillary Pressure VAN_GENUCHTEN_EXTERNAL model.\n");
        }
        ii = 0;
        for (j = 0; j < efv->Num_external_field; j++) {
          if (!strcmp(efv->name[j], input)) {
            ii = 1;
            mat_ptr->por_shell_cap_pres_ext_field_index[ipore] = j;
          }
        }
        if (ii == 0) {
          GOMA_EH(GOMA_ERROR, "Cannot match the name with that in the external field file");
        }
        mat_ptr->PorousShellCapPresModel[ipore] = VAN_GENUCHTEN_EXTERNAL;
        num_const = read_constants(imp, mat_ptr->u_PorousShellCapPres, porous_shell_no);
        if (num_const < 8) {
          sr = sprintf(err_msg, "Matl %s expected at least 8 constants for %s %s model.\n",
                       pd_glob[mn]->MaterialName, "Capillary Pressure", "VAN_GENUCHTEN_EXTERNAL");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
        mat_ptr->len_u_PorousShellCapPres[ipore] = num_const;
      } else if (model_read == -1 && !strcmp(model_name, "VAN_GENUCHTEN_HYST")) {
        mat_ptr->PorousShellCapPresModel[ipore] = VAN_GENUCHTEN_HYST;
        num_const = read_constants(imp, mat_ptr->u_PorousShellCapPres, porous_shell_no);
        if (num_const < 10) {
          sr = sprintf(err_msg, "Matl %s expected at least 10 constants for %s %s model.\n",
                       pd_glob[mn]->MaterialName, "Capillary Pressure", "VAN_GENUCHTEN_HYST");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
        mat_ptr->len_u_PorousShellCapPres[ipore] = num_const;

        /* New change in 2021, use external field to store history of the curve.
         * For now, use fixed names of the field
         */
        if (efv->ev) {
          for (i = 0; i < efv->Num_external_field; i++) {

            /* Search for number of curve switches */
            if (!strcmp(efv->name[i], "NUM_SWITCH_1")) {
              mat_ptr->por_shell_cap_pres_hyst_num_switch_ext_field_index[0] = i;
            }

            if (!strcmp(efv->name[i], "NUM_SWITCH_2")) {
              mat_ptr->por_shell_cap_pres_hyst_num_switch_ext_field_index[1] = i;
            }

            if (!strcmp(efv->name[i], "NUM_SWITCH_3")) {
              mat_ptr->por_shell_cap_pres_hyst_num_switch_ext_field_index[2] = i;
            }

            /* Search for curve type - drainage or imbibition */
            if (!strcmp(efv->name[i], "CURVE_TYPE_1")) {
              mat_ptr->por_shell_cap_pres_hyst_curve_type_ext_field_index[0] = i;
            }
            if (!strcmp(efv->name[i], "CURVE_TYPE_2")) {
              mat_ptr->por_shell_cap_pres_hyst_curve_type_ext_field_index[1] = i;
            }
            if (!strcmp(efv->name[i], "CURVE_TYPE_3")) {
              mat_ptr->por_shell_cap_pres_hyst_curve_type_ext_field_index[2] = i;
            }
          }
        } else {
          GOMA_EH(GOMA_ERROR, " You need external fields to use VAN_GENUCHTEN_HYST");
        }

      } else if (model_read == -1 && !strcmp(model_name, "VAN_GENUCHTEN_HYST_EXT")) {
        if (fscanf(imp, "%s", input) != 1) {
          GOMA_EH(
              GOMA_ERROR,
              "Expecting trailing keyword for Capillary Pressure VAN_GENUCHTEN_HYST_EXT model.\n");
        }
        ii = 0;
        for (j = 0; j < efv->Num_external_field; j++) {
          if (!strcmp(efv->name[j], input)) {
            ii = 1;
            mat_ptr->por_shell_cap_pres_ext_field_index[ipore] = j;
          }
        }
        if (ii == 0) {
          GOMA_EH(GOMA_ERROR, "Cannot match the name with that in the external field file");
        }

        mat_ptr->PorousShellCapPresModel[ipore] = VAN_GENUCHTEN_HYST_EXT;
        num_const = read_constants(imp, mat_ptr->u_PorousShellCapPres, porous_shell_no);
        if (num_const < 18) {
          sr = sprintf(err_msg, "Matl %s expected at least 18 constants for %s %s model.\n",
                       pd_glob[mn]->MaterialName, "Capillary Pressure", "VAN_GENUCHTEN_HYST_EXT");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
        mat_ptr->len_u_PorousShellCapPres[ipore] = num_const;
      } else {
        GOMA_EH(model_read, "Capillary Pressure: Is card missing?");
      }
      if (porous_shell_no != ipore) {
        GOMA_EH(GOMA_ERROR, "Incomplete number of Capillary Pressure card");
      }
      ECHO(es, echo_file);

    } /* End of loop over porous shell layers*/

    /*
     *   Porous Mass Lumping:
     *
     *    This is where you specify whether you want to use
     *    Mass Lumping of the time derivative or whether you
     *    want to use a consistent time derivative treatment.
     *    By default, mass lumping is always TRUE
     *
     */

    model_read = look_for_mat_prop(imp, "Porous Mass Lumping", NULL, NULL, NO_USER, NULL,
                                   model_name, SCALAR_INPUT, &NO_SPECIES, es);

    if (!strcasecmp(model_name, "yes") || !strcasecmp(model_name, "true")) {
      mat_ptr->Porous_Mass_Lump = TRUE;
    } else if (!strcasecmp(model_name, "no") || !strcasecmp(model_name, "false")) {
      mat_ptr->Porous_Mass_Lump = FALSE;
    } else {
      mat_ptr->Porous_Mass_Lump = TRUE;
    }

  } /*end of if( (mat_ptr->PorousMediaType == POROUS_SHELL_UNSATURATED) &&
                 (!PorousShellOn) ) */

  /* Special constants for sink models formulation -- PRS 8/19/05 */
  if (have_por_sink_mass == 1) {
    model_read =
        look_for_mat_prop(imp, "Sink Adsorption Rate Data", &(mat_ptr->PorousSinkConstantsModel),
                          &(mat_ptr->porous_sink_constants), NO_USER, NULL, model_name,
                          SCALAR_INPUT, &NO_SPECIES, es);
    if (model_read == -1 && !strcmp(model_name, "LINEAR")) {
      mat_ptr->PorousSinkConstantsModel = LINEAR;

      num_const = read_constants(imp, &(mat_ptr->u_porous_sink_constants), NO_SPECIES);
      if (num_const < 8) {
        sr = sprintf(err_msg, "Matl %s needs 8 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Sink Adsorption Rate Data", "LINEAR");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_porous_sink_constants = num_const;
    } else if (model_read == -1 && !strcmp(model_name, "POWER_LAW")) {
      mat_ptr->PorousSinkConstantsModel = POWER_LAW;

      num_const = read_constants(imp, &(mat_ptr->u_porous_sink_constants), NO_SPECIES);
      if (num_const < 4) {
        sr = sprintf(err_msg, "Matl %s needs 4 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Sink Adsorption Rate Data", "POWER_LAW");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_porous_sink_constants = num_const;
    } else if (!strcmp(model_name, "CONSTANT")) {
      GOMA_EH(GOMA_ERROR, "Ironically we don't allow a CONSTANT model for Sink Adsorption Rate "
                          "Data.  Try LINEAR");
    } else {
      GOMA_EH(model_read, "Sink Adsorption Rate Data");
    }
    ECHO(es, echo_file);
  }

  /*
   * *****************Species Stuff*******************************
   */

  /* Check for Boundary Condition which require constants from the
   * materials file.
   * Molecular Weight, Specific and Molar Volumes are required when
   * VL_POLY_BC or YFLUX_EQUIL_BC is specified.
   * Instead of reading the constants from the BC card itself, they
   * are specified in the material's file.
   * read_bc_mp != -1 indicates that some mat prop is needed for BC card.
   * ACSun 8/21/98
   */

  read_bc_mp = -1;
  for (i = 0; i < Num_BC; i++) {
    if (BC_Types[i].BC_Name == VL_POLY_BC || BC_Types[i].BC_Name == YFLUX_EQUIL_BC ||
        BC_Types[i].BC_Name == VL_EQUIL_PRXN_BC) {
      read_bc_mp = i;
    }
  }

  /* Check for existance of shear rate and vorticity variable in any matrix*/

  have_shear_rate = 0;
  //*have_vort_dir = 0;*/

  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    if (pd_glob[mn]->e[imtrx][R_SHEAR_RATE]) {
      have_shear_rate = 1;
    }
    /*if (pd_glob[mn]->e[imtrx][R_VORT_DIR1])
      {
       have_vort_dir = 1;
      }*/
  }

  if (pd_glob[mn]->Num_Species_Eqn > 0) {

    /*
     * Hey! If it's in there, then it should be fine. To avoid the ordering
     * complaints, let's rewind the mat file and start looking from the
     * beginning.
     */
    rewind(imp);

    /*
     *  Optionally read in the number of species in the material, and
     *  check for compatibility with the problem statement
     */
    (void)look_for_optional_int(imp, "Number of Species", &model_read, pd_glob[mn]->Num_Species);
    if (model_read != pd_glob[mn]->Num_Species) {
      if (model_read == pd_glob[mn]->Num_Species + 1) {
        pd_glob[mn]->Num_Species = model_read;
        mp_glob[mn]->Num_Species = model_read;
        if (upd->Max_Num_Species < mp_glob[mn]->Num_Species) {
          upd->Max_Num_Species = mp_glob[mn]->Num_Species;
        }
        SPF(es, "WARNING: Number of species in phase %d increased by mat file read to %d\n", mn,
            model_read);
        ECHO(es, echo_file);
      } else {
        fprintf(stderr, "ERROR: number of species in goma input and mat file differ: %d %d\n",
                pd_glob[mn]->Num_Species, model_read);
        GOMA_EH(GOMA_ERROR, "read_input_mp error");
      }
    }
    SPF(es, "%s = %d", "Number of Species", model_read);
    ECHO(es, echo_file);

    /*
     *     Read in the Diffusion Constitutive Equation Model
     */

    strcpy(search_string, "Diffusion Constitutive Equation");

    model_read = look_for_mat_prop(imp, search_string, &(DiffusionConstitutiveEquation), &(a0),
                                   NO_USER, NULL, model_name, NO_INPUT, &NO_SPECIES, es);
    /* EDW: Insert alternate call here! */
    if (!strcmp(model_name, "FICKIAN")) {
      DiffusionConstitutiveEquation = FICKIAN;
    } else if (!strcmp(model_name, "FICKIAN_SHELL")) {
      DiffusionConstitutiveEquation = FICKIAN_SHELL;
    } else if (!strcmp(model_name, "GENERALIZED_FICKIAN")) {
      DiffusionConstitutiveEquation = GENERALIZED_FICKIAN;
    } else if (!strcmp(model_name, "STEFAN_MAXWELL")) {
      DiffusionConstitutiveEquation = STEFAN_MAXWELL;
    } else if (!strcmp(model_name, "STEFAN_MAXWELL_CHARGED")) {
      DiffusionConstitutiveEquation = STEFAN_MAXWELL_CHARGED;
    } else if (!strcmp(model_name, "STEFAN_MAXWELL_VOLUME")) /*  RSL 6/29/00  */
    {
      DiffusionConstitutiveEquation = STEFAN_MAXWELL_VOLUME;
    } else if (!strcmp(model_name, "FICKIAN_CHARGED")) {
      DiffusionConstitutiveEquation = FICKIAN_CHARGED;
    } else if (!strcmp(model_name, "FICKIAN_CHARGED_X")) /*  RSL 9/18/00  */
    {
      DiffusionConstitutiveEquation = FICKIAN_CHARGED_X;
    } else if (!strcmp(model_name, "DARCY")) {
      DiffusionConstitutiveEquation = DARCY;
    } else if (!strcmp(model_name, "DARCY_FICKIAN")) {
      DiffusionConstitutiveEquation = DARCY_FICKIAN;
    } else if (!strcmp(model_name, "HYDRODYNAMIC")) {
      DiffusionConstitutiveEquation = HYDRODYNAMIC;

      if (have_shear_rate == 0)
        GOMA_EH(GOMA_ERROR, "HYDRODYNAMIC mass flux requires shear_rate dof in EQ list.");
    } else if (!strcmp(model_name, "HYDRODYNAMIC_QTENSOR")) {
      DiffusionConstitutiveEquation = HYDRODYNAMIC_QTENSOR;

      if (have_shear_rate == 0)
        GOMA_EH(GOMA_ERROR, "HYDRODYNAMIC_QTENSOR mass flux requires shear_rate dof in EQ list.");
    } else if (!strcmp(model_name, "HYDRODYNAMIC_QTENSOR_OLD")) {
      DiffusionConstitutiveEquation = HYDRODYNAMIC_QTENSOR_OLD;

      if (have_shear_rate == 0)
        GOMA_EH(GOMA_ERROR, "HYDRODYNAMIC_QTENSOR mass flux requires shear_rate dof in EQ list.");
    } else if (!strcmp(model_name, "SUSPENSION_BALANCE")) {
      DiffusionConstitutiveEquation = DM_SUSPENSION_BALANCE;

      /*   if( have_vort_dir == 0 ) */
      /* 	    GOMA_EH(GOMA_ERROR, "SUSPENSION_BALANCE mass flux requires a vorticity vector
       * in EQ list."); */
    } else if (!strcmp(model_name, "NONE")) {
      DiffusionConstitutiveEquation = NON_DIFFUSING;
    } else {
      sprintf(err_msg, "Matl %s, \"%s\" = \"%s\"?\n(valid: %s)\n[%s]", pd_glob[mn]->MaterialName,
              "Diffusion Constitutive Equation", model_name,
              "FICKIAN GENERALIZED_FICKIAN STEFAN_MAXWELL STEFAN_MAXWELL_CHARGED FICKIAN_CHARGED "
              "FICKIAN_CHARGED_X STEFAN_MAXWELL_VOLUME DARCY DARCY_FICKIAN HYDRODYNAMIC "
              "HYDRODYNAMIC_QTENSOR NONE", /*  RSL 9/18/00  */
              "If still bad - check orderings in the mat file!");
      GOMA_EH(GOMA_ERROR, err_msg);
    }

    pd_glob[mn]->MassFluxModel = DiffusionConstitutiveEquation;
    cr_glob[mn]->MassFluxModel = DiffusionConstitutiveEquation;

    ECHO(es, echo_file);

    strcpy(search_string, "PBE Blowing Agent Type");
    model_read = look_for_mat_prop(imp, search_string, &(mat_ptr->PBE_BA_Type), &(a0), NO_USER,
                                   NULL, model_name, NO_INPUT, &NO_SPECIES, es);

    if (!strcmp(model_name, "N_PENTANE")) {
      mat_ptr->PBE_BA_Type = PBE_N_PENTANE;
    } else if (!strcmp(model_name, "R_11")) {
      mat_ptr->PBE_BA_Type = PBE_R_11;
    } else if (!strcmp(model_name, " ")) // default
    {
      mat_ptr->PBE_BA_Type = PBE_N_PENTANE;
    } else {
      GOMA_EH(GOMA_ERROR, "Unknown PBA Blowing Agent Type");
    }

    ECHO(es, echo_file);
  }

  /* Parameters for Ryan's Q tensor model */
  model_read = look_for_mat_prop(
      imp, "Qtensor extension pressure", &(mat_ptr->QtensorExtensionPModel),
      &(mat_ptr->Qtensor_Extension_P), NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);
  ECHO(es, echo_file);

  model_read =
      look_for_mat_prop(imp, "Qtensor Nct", &(mat_ptr->QtensorNctModel), &(mat_ptr->Qtensor_Nct),
                        NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);
  ECHO(es, echo_file);

  model_read = look_for_mat_prop(imp, "Species Weight Function", &(mat_ptr->Spwt_funcModel),
                                 &(mat_ptr->Spwt_func), NO_USER, NULL, model_name, SCALAR_INPUT,
                                 &NO_SPECIES, es);

  if (!strcmp(model_name, "GALERKIN")) {
    mat_ptr->Spwt_funcModel = GALERKIN;
    mat_ptr->Spwt_func = 0.;
  } else if (!strcmp(model_name, "SUPG")) {
    int err;
    mat_ptr->Spwt_funcModel = SUPG;
    err = fscanf(imp, "%lg", &(mat_ptr->Spwt_func));
    if (err != 1) {
      GOMA_EH(GOMA_ERROR, "Expected to read one double for Species Weight Function SUPG");
    }
  } else if (!strcmp(model_name, "SUPG_GP")) {
    int err;
    mat_ptr->Spwt_funcModel = SUPG_GP;
    err = fscanf(imp, "%lg", &(mat_ptr->Spwt_func));
    if (err != 1) {
      GOMA_EH(GOMA_ERROR, "Expected to read one double for Species Weight Function SUPG_GP");
    }
  } else if (!strcmp(model_name, "SUPG_SHAKIB")) {
    int err;
    mat_ptr->Spwt_funcModel = SUPG_SHAKIB;
    err = fscanf(imp, "%lg", &(mat_ptr->Spwt_func));
    if (err != 1) {
      GOMA_EH(GOMA_ERROR, "Expected to read one double for Species Weight Function SUPG_SHAKIB");
    }
  } else {
    mat_ptr->Spwt_funcModel = GALERKIN;
    mat_ptr->Spwt_func = 0.;
    SPF(es, "\t(%s = %s)", "Species Weight Function", "GALERKIN");
  }
  ECHO(es, echo_file);

  model_read = look_for_mat_prop(imp, "Species SSPG Function", &(mat_ptr->SpSSPG_funcModel),
                                 &(mat_ptr->SpSSPG_func), NO_USER, NULL, model_name, SCALAR_INPUT,
                                 &NO_SPECIES, es);
  ECHO(es, echo_file);

  model_read = look_for_mat_prop(imp, "Species YZbeta Function", &(mat_ptr->SpYZbeta_funcModel),
                                 &(mat_ptr->SpYZbeta_func), NO_USER, NULL, model_name, SCALAR_INPUT,
                                 &NO_SPECIES, es);
  if (!strcmp(model_name, "ONE")) {
    mat_ptr->SpYZbeta_funcModel = YZBETA_ONE;
    if (fscanf(imp, "%lg", &(mat_ptr->SpYZbeta_func)) != 1) {
      GOMA_EH(GOMA_ERROR, "Could not read Scale for Species YZbeta Function YZBETA_ONE");
    }
  } else if (!strcmp(model_name, "TWO")) {
    mat_ptr->SpYZbeta_funcModel = YZBETA_TWO;
    if (fscanf(imp, "%lg", &(mat_ptr->SpYZbeta_func)) != 1) {
      GOMA_EH(GOMA_ERROR, "Could not read Scale for Species YZbeta Function YZBETA_TWO");
    }
  } else if (!strcmp(model_name, "MIXED")) {
    mat_ptr->SpYZbeta_funcModel = YZBETA_MIXED;
    if (fscanf(imp, "%lg", &(mat_ptr->SpYZbeta_func)) != 1) {
      GOMA_EH(GOMA_ERROR, "Could not read Scale for Species YZbeta Function YZBETA_MIXED");
    }
  } else if (!strcmp(model_name, "CUSTOM")) {
    mat_ptr->SpYZbeta_funcModel = YZBETA_CUSTOM;
    if (fscanf(imp, "%lg %lg", &(mat_ptr->SpYZbeta_func), &(mat_ptr->SpYZbeta_value)) != 2) {
      GOMA_EH(GOMA_ERROR,
              "Could not read Scale and beta value for Species YZbeta Function YZBETA_CUSTOM");
    }
  } else {
    mat_ptr->SpYZbeta_funcModel = SC_NONE;
    mat_ptr->SpYZbeta_func = 0.;
    SPF(es, "\t(%s = %s)", "Species YZbeta Function", "NONE");
  }
  ECHO(es, echo_file);

  /*
   *   Special section (LONG) to read in parameters associated with Stefan-Maxwell diffusion
   *   of neutral and/or charged species in concentrated solutions (as in thermal batteries).
   *
   *   KSC: added (7/98), revised (9/98, 9/2000).
   *
   */
  if (DiffusionConstitutiveEquation == STEFAN_MAXWELL ||
      DiffusionConstitutiveEquation == STEFAN_MAXWELL_CHARGED ||
      DiffusionConstitutiveEquation == STEFAN_MAXWELL_VOLUME) {

    if (pd_glob[mn]->Num_Species_Eqn + 1 != pd_glob[mn]->Num_Species) {
      fprintf(stderr, "ERROR stefan Maxwell diffusion model chosen but number of species, %d\n",
              pd_glob[mn]->Num_Species);
      fprintf(stderr, "\t isn't one less than number of species equations, %d\n",
              pd_glob[mn]->Num_Species_Eqn);
      GOMA_EH(GOMA_ERROR, "Error in number of species and species eqn specs");
    }
    n_species = pd_glob[mn]->Num_Species;

    if (n_species < 2) {
      GOMA_EH(
          GOMA_ERROR,
          "Error: Stefan_Maxwell model should be used for modeling transport of 2 or more species");
    }

    iread = look_for_optional(imp, "Diffusivity", input, '=');
    if (fscanf(imp, "%s", model_name) != 1) {
      GOMA_EH(GOMA_ERROR, "Error: need to specify diffusivity name, e.g. CONSTANT");
    }

    SPF(es, "%s = %s", "Diffusivity", model_name);

    if (DiffusionConstitutiveEquation == STEFAN_MAXWELL ||
        DiffusionConstitutiveEquation == STEFAN_MAXWELL_CHARGED ||
        DiffusionConstitutiveEquation == STEFAN_MAXWELL_VOLUME) /* added by RSL, 9/14/00 */
    {
      ECHO(es, echo_file);
      if (!strcmp(model_name, "CONSTANT")) {
        for (j = 0; j < n_species; j++)
          mat_ptr->DiffusivityModel[j] = CONSTANT;
        n_dij = (n_species * n_species - n_species) / 2;

        for (i = 0; i < n_dij; i++) /* reading the Stefan-Maxwell diffusivities */
        {
          if (fscanf(imp, "%d %d %lf", &ii, &jj, &dij) != 3) {
            GOMA_EH(GOMA_ERROR,
                    "Error in reading Stefan_Maxwell diffusivities: need to input i, j, and Dij");
          }
          mat_ptr->diffusivity_Stefan_Maxwell[ii][jj] = dij;
          mat_ptr->diffusivity_Stefan_Maxwell[jj][ii] = dij;
          SPF(es, "\t\t %d %d %.4g", ii, jj, dij);
          ECHO(es, echo_file);
        }
      } else if (!strcmp(model_name, "ARRHENIUS")) {
        for (j = 0; j < n_species; j++)
          mat_ptr->DiffusivityModel[j] = CONSTANT;
        n_dij = (n_species * n_species - n_species) / 2;
        for (i = 0; i < n_dij; i++) /* reading the Stefan-Maxwell diffusivities */
        {
          if (fscanf(imp, "%d %d %lf %lf %lf", &ii, &jj, &dij, &E, &T0) != 5) {
            GOMA_EH(GOMA_ERROR, "Error: need to input ii, jj, dij, E, T0");
          }
          mat_ptr->diffusivity_Stefan_Maxwell[ii][jj] = dij;
          mat_ptr->diffusivity_Stefan_Maxwell[jj][ii] = dij;
          mat_ptr->u_diffusivity_Stefan_Maxwell[ii][jj][0] = dij;
          mat_ptr->u_diffusivity_Stefan_Maxwell[jj][ii][0] = dij;
          mat_ptr->u_diffusivity_Stefan_Maxwell[ii][jj][1] = E;
          mat_ptr->u_diffusivity_Stefan_Maxwell[jj][ii][1] = E;
          mat_ptr->u_diffusivity_Stefan_Maxwell[ii][jj][2] = T0;
          mat_ptr->u_diffusivity_Stefan_Maxwell[jj][ii][2] = T0;
          SPF(es, "\t\t %d %d %.4g %.4g %.4g", ii, jj, dij, E, T0);
          ECHO(es, echo_file);
        }
      } else {
        GOMA_EH(GOMA_ERROR, "This S-M diffusivity model has not been implemented yet!");
      }
    }

    /*
     *  Get the number of chemical reactions involved in the present
     *  material block; KSC/GHE: 10/98
     */
    iread = look_for_optional(imp, "Number of chemical reactions", input, '=');
    if (fscanf(imp, "%d", &n_rxn) != 1) {
      GOMA_EH(GOMA_ERROR, "Expected to read 1 int for \"Number of chemical reactions\"");
    }

    pd_glob[mn]->Num_Rxn = n_rxn;
    if (n_rxn > MAX_RXN) {
      sprintf(err_msg,
              "Specified number of chemical rxns %d > %d, compiled limit - boost MAX_RXN in "
              "rf_fem_const.h",
              ii, MAX_RXN);
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    SPF(es, "%s = %d", "Number of chemical reactions", n_rxn);
    ECHO(es, echo_file);

    /*
     *     if n_rxn > 0, then read Butler Volmer kinetics parameters;
     *        KSC: 10/98; revised: KSC 3/99
     *
     *
     */
    if (n_rxn > 0) {
      model_read = look_for_mat_prop(imp, "Reaction Rate", &(mat_ptr->ReactionRateModel),
                                     &(mat_ptr->reaction_rate), NO_USER, NULL, model_name,
                                     SCALAR_INPUT, &NO_SPECIES, es);

      if (model_read == -1 && !strcmp(model_name, "ELECTRODE_KINETICS")) {
        mat_ptr->ReactionRateModel = ELECTRODE_KINETICS;
        num_const = read_constants(imp, &(mat_ptr->u_reaction_rate), 0);
        if (num_const < 3) {
          sprintf(err_msg, "Material %s - expected at least 3 constants for %s %s model.\n",
                  pd_glob[mn]->MaterialName, "Reaction Rate", "ELECTRODE_KINETICS");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
        mat_ptr->len_u_reaction_rate = num_const;
        SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_reaction_rate);
      } else {
        GOMA_EH(GOMA_ERROR, "Error: reaction-rate models other than ELECTRODE_KINETICS awaits "
                            "future implementation");
      }

      model_read =
          look_for_mat_prop(imp, "Thermodynamic Potential", &(mat_ptr->ThermodynamicPotentialModel),
                            &(mat_ptr->thermodynamic_potential), NO_USER, NULL, model_name,
                            SCALAR_INPUT, &NO_SPECIES, es);

      if (model_read == -1 && !strcmp(model_name, "FeS2")) {
        mat_ptr->ThermodynamicPotentialModel = FeS2;
        num_const = read_constants(imp, &(mat_ptr->u_thermodynamic_potential), 0);
        if (num_const < 8) {
          sprintf(err_msg, "Material %s - expected at least 8 constants for %s %s model.\n",
                  pd_glob[mn]->MaterialName, "Thermodynamic Potential", "FeS2");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
        mat_ptr->len_u_thermodynamic_potential = num_const;
        SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_thermodynamic_potential);
      } else if (model_read == -1 && !strcmp(model_name, "LiSi")) {
        mat_ptr->ThermodynamicPotentialModel = LiSi;
        num_const = read_constants(imp, &(mat_ptr->u_thermodynamic_potential), 0);
        if (num_const < 7) {
          sprintf(err_msg, "Material %s - expected at least 7 constants for %s %s model.\n",
                  pd_glob[mn]->MaterialName, "Thermodynamic Potential", "LiSi");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
        mat_ptr->len_u_thermodynamic_potential = num_const;
        SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_thermodynamic_potential);
      } else if (model_read == -1 && !strcmp(model_name, "CONSTANT")) {
        mat_ptr->ThermodynamicPotentialModel = CONSTANT;
        num_const = read_constants(imp, &(mat_ptr->u_thermodynamic_potential), 0);
        if (num_const < 1) {
          sprintf(err_msg, "Material %s - expected at least 1 constant for %s %s model.\n",
                  pd_glob[mn]->MaterialName, "Thermodynamic Potential", "CONSTANT");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
        mat_ptr->len_u_thermodynamic_potential = num_const;
        SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_thermodynamic_potential);
      }

      ECHO(es, echo_file);

      model_read = look_for_mat_prop(imp, "Interfacial Area", &(mat_ptr->InterfacialAreaModel),
                                     &(mat_ptr->interfacial_area), NO_USER, NULL, model_name,
                                     SCALAR_INPUT, &NO_SPECIES, es);

      if (model_read == -1 && !strcmp(model_name, "THERMAL_BATTERY")) {
        mat_ptr->InterfacialAreaModel = THERMAL_BATTERY;
        num_const = read_constants(imp, &(mat_ptr->u_interfacial_area), 0);
        if (num_const < 9) {
          sprintf(err_msg, "Material %s - expected at least 9 constants for %s %s model.\n",
                  pd_glob[mn]->MaterialName, "Interfacial Area", "THERMAL_BATTERY");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
        mat_ptr->len_u_interfacial_area = num_const;
        SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_interfacial_area);
      } else if (model_read == -1 && !strcmp(model_name, "CONSTANT")) {
        mat_ptr->InterfacialAreaModel = CONSTANT;
        num_const = read_constants(imp, &(mat_ptr->u_interfacial_area), 0);
        if (num_const < 1) {
          sprintf(err_msg, "Material %s - expected at least 1 constant for %s %s model.\n",
                  pd_glob[mn]->MaterialName, "Interfacial Area", "CONSTANT");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
        mat_ptr->len_u_interfacial_area = num_const;
        SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_interfacial_area);
      }

      ECHO(es, echo_file);

    } /* end if (n_rxn >0) statement; KSC: 10/98, 2/02 */

    /*
     *    Read the solution temperature for systems
     *    involving charged species; KSC (9/98, 9/2000).
     *
     *     KEN ! YOUR A NICE GUY BUT I'M GOING TO CUT YOU !
     *                                                      - tab
     *
     */
    if (DiffusionConstitutiveEquation == STEFAN_MAXWELL_CHARGED ||
        DiffusionConstitutiveEquation == STEFAN_MAXWELL_VOLUME) {
      model_read =
          look_for_mat_prop(imp, "Solution Temperature", &(mat_ptr->SolutionTemperatureModel),
                            &(mat_ptr->solution_temperature), NO_USER, NULL, model_name,
                            SCALAR_INPUT, &NO_SPECIES, es);

      if (model_read == -1 && !strcmp(model_name, "THERMAL_BATTERY")) {
        mat_ptr->SolutionTemperatureModel = THERMAL_BATTERY;
        num_const = read_constants(imp, &(mat_ptr->u_solution_temperature), 0);
        if (num_const <
            6) { /* should 6 constants for the thermal_battery solution-temperature model */
          sprintf(err_msg, "Material %s - expected at least 6 constants for %s %s model.\n",
                  pd_glob[mn]->MaterialName, "Solution Temperature", "THERMAL_BATTERY");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
        mat_ptr->len_u_solution_temperature = num_const;
        SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_solution_temperature);
      }
    } else {
      /*
       * for a solution with neutral species, set
       * the solution temperature to 298 K
       * The temperature is irrelevant in solution with neutral species
       */
      mat_ptr->SolutionTemperatureModel = CONSTANT;
      mat_ptr->solution_temperature = 298.0;
      SPF(es, "\t(%s = %s %.4g)", "Solution Temperature", "CONSTANT",
          mat_ptr->solution_temperature);
    }

    ECHO(es, echo_file);

    model_read = look_for_mat_prop(imp, "Porosity", &(mat_ptr->PorosityModel), &(mat_ptr->porosity),
                                   NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);
    if (model_read == -1 && !strcmp(model_name, "THERMAL_BATTERY")) {
      mat_ptr->PorosityModel = THERMAL_BATTERY;
      num_const = read_constants(imp, &(mat_ptr->u_porosity), 0);
      if (num_const < 2) {
        sprintf(err_msg, "Material %s - expected at least 2 constants for %s %s model.\n",
                pd_glob[mn]->MaterialName, "Porosity", "THERMAL_BATTERY");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_porosity = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_porosity);
    }

    ECHO(es, echo_file);

  } /* end of mp-input reading for STEFAN_MAXWELL, STEFAN_MAXWELL_ChARGED
       KSC: (7/98, 9/2000)
       *
       *  AND NOT A MOMENT TOO SOON !
       */

  /*
   * Read diffusivity and other properties for non-Stefan-Maxwell models;
   * it is captured in the while-loop to avoid conflict.  Added ACSun 9/98
   */

  /*
   ********************************************************************************
   *  LOOP OVER THE SPECIES EQUATIONS
   ********************************************************************************
   */

  /*
  * HKM  - Changed the following parameter to Num_Species instead of Num_Species_Eqn
  *        For dilute problems the two are the same. For nondilute problems, the
  *        material properties for the last species in the mechanism should be
  *        read in as well, even though there isn't an explicit conservation
  *        equation solved for it.
  *

  */

  for (j = 0; j < mat_ptr->Num_Species; j++) {
    if (DiffusionConstitutiveEquation != STEFAN_MAXWELL &&
        DiffusionConstitutiveEquation != STEFAN_MAXWELL_CHARGED &&
        DiffusionConstitutiveEquation != STEFAN_MAXWELL_VOLUME) {
      model_read = look_for_species_proptable(
          imp, "Diffusivity", mat_ptr, mat_ptr->DiffusivityModel, mat_ptr->diffusivity,
          mat_ptr->u_diffusivity, mat_ptr->len_u_diffusivity, &(mat_ptr->diffusivity_tableid[j]),
          model_name, 0, &species_no, es);

      fallback_chemkin_generic_prop(&model_read, j, &(mat_ptr->DiffusivityModel[j]), FALSE,
                                    mat_ptr);

      ECHO(es, echo_file);
      /*
       * Postprocess unique Diffusivity models
       */
      if (model_read == -1) {
        GOMA_EH(GOMA_ERROR,
                "Diffusivity: Bad Card syntax or need another set of species mat cards?");
      }

      else if (model_read == 0) {
        if (!strcmp(model_name, "POROUS")) {
          mat_ptr->DiffusivityModel[species_no] = POROUS;
          num_const = mat_ptr->len_u_diffusivity[species_no];
          if (num_const < 5) {
            sr = sprintf(err_msg, "Matl %s (conc %d) needs at least 5 constants for %s %s model.\n",
                         pd_glob[mn]->MaterialName, species_no, "Diffusivity", "POROUS");
            GOMA_EH(GOMA_ERROR, err_msg);
          }
        }

        /*
         * HYDRO is based on diffusive-flux model with contributions from shear gradient
         * (gammadot), viscosity gradient, curvature, and hindered settling.  Fickian
         * diffusion is available for rough convergence spots.
         *
         *  Thus, if HYDRO is specified, then the following additional parameters are
         *  looked for and must be supplied:
         *     "Shear Rate Diffusivity"
         *     "Viscosity Diffusivity"
         *     "Fickian Diffusivity"
         *     "Gravity-based Diffusivity"
         *
         *  For the Q-tensor diffusivity model, the following need to be defined:
         *     "Shear Rate Diffusivity"
         *     "Viscosity Diffusivity"
         *     "Q Tensor Diffusivity"
         *
         *  The Q Tensor diffusivity card takes the form of a model
         *  name (which currently must be CONSTANT or NONE), species
         *  number and then 3 floats (representing the flux modifier
         *  in the flow, normal, and vorticity directions).  i.e.,
         *
         *      Q Tensor Diffusivity = CONSTANT 0 1.0 1.0 0.5
         */
        else if (!strcmp(model_name, "HYDRO")) {

          mat_ptr->DiffusivityModel[species_no] = HYDRO;

          species_no = mat_ptr->Num_Species; /* set species number equal to max number of species
                                            it is changed to species number of input property
                                            by look_for_mat_prop */

          model_read = look_for_mat_prop(imp, "Shear Rate Diffusivity", mat_ptr->GamDiffType,
                                         mat_ptr->gam_diffusivity, mat_ptr->u_gadiffusivity, NULL,
                                         model_name, SCALAR_INPUT, &species_no, es);

          if (model_read == -1 && !strcmp(model_name, "LINEAR")) {
            mat_ptr->GamDiffType[species_no] = LINEAR;
            num_const = read_constants(imp, mat_ptr->u_gadiffusivity, species_no);
            mat_ptr->len_u_gadiffusivity[species_no] = num_const;
            SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_gadiffusivity[species_no]);
          } else if (model_read == -1 && !strcmp(model_name, "CONST_LS")) {
            mat_ptr->GamDiffType[species_no] = LEVEL_SET;
            num_const = read_constants(imp, mat_ptr->u_gadiffusivity, species_no);

            if (num_const < 3) {
              sr =
                  sprintf(err_msg, "Matl %s %s needs 3 constants for %s  model: Kg and exponent.\n",
                          pd_glob[mn]->MaterialName, "Shear Rate Diffusivity", "CONST_LS");
              GOMA_EH(GOMA_ERROR, err_msg);
            }

            mat_ptr->len_u_gadiffusivity[species_no] = num_const;

            SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_gadiffusivity[species_no]);
          } else {
            sr = sprintf(err_msg, "Material %s - unrecognized model for %s \"%s\" ???\n",
                         pd_glob[mn]->MaterialName, "Shear Rate Diffusivity", model_name);
            GOMA_EH(model_read, err_msg);
          }

          ECHO(es, echo_file);

          species_no = mat_ptr->Num_Species; /* set species number equal to max number of species
                                          it is changed to species number of input property
                                          by look_for_mat_prop */

          model_read = look_for_mat_prop(imp, "Viscosity Diffusivity", mat_ptr->MuDiffType,
                                         mat_ptr->mu_diffusivity, mat_ptr->u_mdiffusivity, NULL,
                                         model_name, SCALAR_INPUT, &species_no, es);

          if (model_read == -1 && !strcmp(model_name, "LINEAR")) {
            mat_ptr->MuDiffType[species_no] = LINEAR;
            num_const = read_constants(imp, mat_ptr->u_mdiffusivity, species_no);
            mat_ptr->len_u_mdiffusivity[species_no] = num_const;
            SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_mdiffusivity[species_no]);
          } else if (model_read == -1 && !strcmp(model_name, "CONST_LS")) {
            mat_ptr->MuDiffType[species_no] = LEVEL_SET;
            num_const = read_constants(imp, mat_ptr->u_mdiffusivity, species_no);

            if (num_const < 3) {
              sr =
                  sprintf(err_msg, "Matl %s %s needs 3 constants for %s  model: Kg and exponent.\n",
                          pd_glob[mn]->MaterialName, "Viscosity Diffusivity", "CONST_LS");
              GOMA_EH(GOMA_ERROR, err_msg);
            }
            mat_ptr->len_u_mdiffusivity[species_no] = num_const;
            SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_mdiffusivity[species_no]);
          }

          else {
            sr = sprintf(err_msg, "Material %s - unrecognized model for %s \"%s\" ???\n",
                         pd_glob[mn]->MaterialName, "Viscosity Diffusivity", model_name);
            GOMA_EH(model_read, err_msg);
          }

          ECHO(es, echo_file);

          species_no = mat_ptr->Num_Species; /* set species number equal to max number of species
                                          it is changed to species number of input property
                                          by look_for_mat_prop */

          model_read = look_for_mat_prop(imp, "Fickian Diffusivity", mat_ptr->FickDiffType,
                                         mat_ptr->f_diffusivity, mat_ptr->u_fdiffusivity, NULL,
                                         model_name, SCALAR_INPUT, &species_no, es);

          if (model_read == -1 && !strcmp(model_name, "ANISOTROPIC")) {
            mat_ptr->FickDiffType[species_no] = ANISOTROPIC;
            num_const = read_constants(imp, mat_ptr->u_fdiffusivity, species_no);
            if (num_const < 3) {
              sr = sprintf(err_msg, "Matl %s %s needs 3 constants for %s  model.\n",
                           pd_glob[mn]->MaterialName, "Fickian Diffusivity", "HYDRO Diffusivity");
              GOMA_EH(GOMA_ERROR, err_msg);
            }
            mat_ptr->len_u_fdiffusivity[species_no] = num_const;
            SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_fdiffusivity[species_no]);
          } else if (model_read == -1 && !strcmp(model_name, "EXP_DECAY")) {
            mat_ptr->FickDiffType[species_no] = EXP_DECAY;
            num_const = read_constants(imp, mat_ptr->u_fdiffusivity, species_no);
            if (num_const < 2) {
              sr = sprintf(err_msg, "Matl %s %s needs 2 constants for %s  model.\n",
                           pd_glob[mn]->MaterialName, "EXP_DECAY Diffusivity", "HYDRO Diffusivity");
              GOMA_EH(GOMA_ERROR, err_msg);
            }
            mat_ptr->len_u_fdiffusivity[species_no] = num_const;
            SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_fdiffusivity[species_no]);
          } else if (model_read == -1 && !strcmp(model_name, "SHOCK")) {
            mat_ptr->FickDiffType[species_no] = SHOCK;
            num_const = read_constants(imp, mat_ptr->u_fdiffusivity, species_no);
            if (num_const < 1) {
              sr = sprintf(err_msg, "Matl %s %s needs 1 constants for %s  model.\n",
                           pd_glob[mn]->MaterialName, "SHOCK Diffusivity", "HYDRO Diffusivity");
              GOMA_EH(GOMA_ERROR, err_msg);
            }
            mat_ptr->len_u_fdiffusivity[species_no] = num_const;
            SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_fdiffusivity[species_no]);
          } else {
            sr = sprintf(err_msg, "Material %s - unrecognized model for %s \"%s\" ???\n",
                         pd_glob[mn]->MaterialName, "Fickian Diffusivity", model_name);
            GOMA_EH(GOMA_ERROR, err_msg);
          }

          ECHO(es, echo_file);

          species_no = mat_ptr->Num_Species; /* set species number equal to max number of species
                                          it is changed to species number of input property
                                          by look_for_mat_prop */

          model_read = look_for_mat_prop(imp, "Gravity-based Diffusivity", mat_ptr->GravDiffType,
                                         mat_ptr->g_diffusivity, mat_ptr->u_gdiffusivity, NULL,
                                         model_name, SCALAR_INPUT, &species_no, es);

          if (model_read == -1 && !strcmp(model_name, "BISECTION")) {
            mat_ptr->GravDiffType[species_no] = BISECTION;
            num_const = read_constants(imp, mat_ptr->u_gdiffusivity, species_no);
            if (num_const < 3) {
              sr = sprintf(err_msg,
                           "Matl %s %s needs 3 constants for %s  model: Kg, avg_conc, slope.\n",
                           pd_glob[mn]->MaterialName, "Gravity-based Diffusivity", "BISECTION");
              GOMA_EH(GOMA_ERROR, err_msg);
            }
            mat_ptr->len_u_gdiffusivity[species_no] = num_const;
            SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_gdiffusivity[species_no]);
          } else if (model_read == -1 && !strcmp(model_name, "RZBISECTION")) {
            mat_ptr->GravDiffType[species_no] = RZBISECTION;
            num_const = read_constants(imp, mat_ptr->u_gdiffusivity, species_no);
            if (num_const < 4) {
              sr = sprintf(
                  err_msg,
                  "Matl %s %s needs 4 constants for %s  model: Kg, exponent, avg_conc, slope .\n",
                  pd_glob[mn]->MaterialName, "Gravity-based Diffusivity", "RZBISECTION");
              GOMA_EH(GOMA_ERROR, err_msg);
            }
            mat_ptr->len_u_gdiffusivity[species_no] = num_const;
            SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_gdiffusivity[species_no]);
          } else if (model_read == -1 && !strcmp(model_name, "RICHARDSON_ZAKI")) {
            mat_ptr->GravDiffType[species_no] = RICHARDSON_ZAKI;
            num_const = read_constants(imp, mat_ptr->u_gdiffusivity, species_no);
            if (num_const < 2) {
              sr = sprintf(
                  err_msg, "Matl %s %s needs 2 constants for %s  model: Kg and exponent.\n",
                  pd_glob[mn]->MaterialName, "Gravity-based Diffusivity", "RICHARDSON_ZAKI");
              GOMA_EH(GOMA_ERROR, err_msg);
            }
            mat_ptr->len_u_gdiffusivity[species_no] = num_const;
            SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_gdiffusivity[species_no]);
          } else if (model_read == -1 && !strcmp(model_name, "CONST_LS")) {
            mat_ptr->GravDiffType[species_no] = RICHARDSON_ZAKI;
            num_const = read_constants(imp, mat_ptr->u_gdiffusivity, species_no);
            if (num_const < 3) {
              sr =
                  sprintf(err_msg, "Matl %s %s needs 3 constants for %s  model: Kg and exponent.\n",
                          pd_glob[mn]->MaterialName, "Gravity-based Diffusivity", "CONST_LS");
              GOMA_EH(GOMA_ERROR, err_msg);
            }
            mat_ptr->len_u_gdiffusivity[species_no] = num_const;
            SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_gdiffusivity[species_no]);
          } else {
            sr = sprintf(err_msg, "Material %s - unrecognized model for %s \"%s\" ???\n",
                         pd_glob[mn]->MaterialName, "Gravity-based Diffusivity", model_name);
            GOMA_EH(model_read, err_msg);
          }

          ECHO(es, echo_file);

          /* MMH: I made this optional so it would be
             backwards-compatible with the other HYDRO diffusivity
             model material files. */
          input[0] = '\0';
          look_for_optional_string(imp, "Q Tensor Diffusivity", input, MAX_CHAR_IN_INPUT);
          s = &(input[0]);
          if (s[0]) {
            strcpy(model_name, strtok(s, " \t"));
            if (!strncmp(model_name, "CONSTANT", 8)) {
              species_no = atoi(strtok(NULL, " \t"));

              mat_ptr->QTensorDiffType[species_no] = CONSTANT;

              for (i = 0; i < 3; i++)
                mat_ptr->q_diffusivity[species_no][i] = atof(strtok(NULL, " \t"));

              mat_ptr->len_u_qdiffusivity[species_no] = 0;

              /* 		 if(mat_ptr->q_diffusivity[species_no][0] != 1.0 || */
              /* 			    mat_ptr->q_diffusivity[species_no][1] != 1.0 || */
              /* 			    mat_ptr->q_diffusivity[species_no][2] != 0.5) */
              /* 			   GOMA_EH(GOMA_ERROR, "Sorry, the Q tensor components are
               * internally hard-coded to 1.0, 1.0, 0.5."); */

              SPF(es, "%s = %s %d", "Q Tensor Diffusivity", model_name, species_no);
              SPF_DBL_VEC(endofstring(es), 3, mat_ptr->q_diffusivity[species_no]);
            } else if (!strncmp(model_name, "NONE", 4)) {
              species_no = atoi(strtok(NULL, " \t"));
              mat_ptr->QTensorDiffType[species_no] = NO_MODEL;
              mat_ptr->len_u_qdiffusivity[species_no] = 0;
              SPF(es, "%s = %s %d", "Q Tensor Diffusivity", "NONE", species_no);
            } else {
              sr = sprintf(err_msg, "Material %s - unrecognized model for %s \"%s\" ???\n",
                           pd_glob[mn]->MaterialName, "Q Tensor Diffusivity", model_name);
              GOMA_EH(GOMA_ERROR, err_msg);
            }
            ECHO(es, echo_file);
          }

          if (model_read == -1) {

            /*
             * And just check to make sure this is active, too.
             */

            if (have_shear_rate == 0) {
              sr = sprintf(err_msg, "Matl %s (conc %d) %s %s model needs the \"%s\" eqn active.\n",
                           pd_glob[mn]->MaterialName, species_no, "Diffusivity", "HYDRO",
                           "shear_rate");
              GOMA_EH(GOMA_ERROR, err_msg);
            }
          }
        }
        /*
         * SUSPENSION is based on diffusive-flux model with contributions from the velocity
         *  gradient tensor, viscosity gradient, particle stress, and hindered settling.  Fickian
         * diffusion is available for rough convergence spots. There are currently no adjustable
         * parameters for this model, but the 2/9a*a, where a is the particle size, must be
         * entered for the settling term and the Fickian diffusivity is included for stability.
         *
         *  Thus, if SUSPENSION is specified, then the following additional parameters are
         *  looked for and must be supplied:
         *     "Fickian Diffusivity"
         *     "Gravity-based Diffusivity"
         *     "Q Tensor Diffusivity"
         */
        else if (!strcmp(model_name, "SUSPENSION")) {
          mat_ptr->DiffusivityModel[species_no] = SUSP_BAL;

          species_no = mat_ptr->Num_Species; /* set species number equal to max number of species
                                                it is changed to species number of input property
                                                by look_for_mat_prop */

          model_read = look_for_mat_prop(imp, "Fickian Diffusivity", mat_ptr->FickDiffType,
                                         mat_ptr->f_diffusivity, mat_ptr->u_fdiffusivity, NULL,
                                         model_name, SCALAR_INPUT, &species_no, es);

          if (model_read == -1 && !strcmp(model_name, "ANISOTROPIC")) {
            mat_ptr->FickDiffType[species_no] = ANISOTROPIC;
            num_const = read_constants(imp, mat_ptr->u_fdiffusivity, species_no);
            if (num_const < 3) {
              sr = sprintf(err_msg, "Matl %s %s needs 3 constants for %s  model.\n",
                           pd_glob[mn]->MaterialName, "Fickian Diffusivity",
                           "SUSPENSION Diffusivity");
              GOMA_EH(GOMA_ERROR, err_msg);
            }
            mat_ptr->len_u_fdiffusivity[species_no] = num_const;
            SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_fdiffusivity[species_no]);
          } else if (model_read == -1 && !strcmp(model_name, "EXP_DECAY")) {
            mat_ptr->FickDiffType[species_no] = EXP_DECAY;
            num_const = read_constants(imp, mat_ptr->u_fdiffusivity, species_no);
            if (num_const < 2) {
              sr = sprintf(err_msg, "Matl %s %s needs 2 constants for %s  model.\n",
                           pd_glob[mn]->MaterialName, "Fickian Diffusivity",
                           "SUSPENSION Diffusivity");
              GOMA_EH(GOMA_ERROR, err_msg);
            }
            mat_ptr->len_u_fdiffusivity[species_no] = num_const;
            SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_fdiffusivity[species_no]);
          } else {
            sr = sprintf(err_msg, "Material %s - unrecognized model for %s \"%s\" ???\n",
                         pd_glob[mn]->MaterialName, "Fickian Diffusivity", model_name);
            GOMA_EH(GOMA_ERROR, err_msg);
          }

          ECHO(es, echo_file);

          species_no = mat_ptr->Num_Species; /* set species number equal to max number of species
                                                it is changed to species number of input property
                                                by look_for_mat_prop */

          model_read = look_for_mat_prop(imp, "Q Tensor Diffusivity", mat_ptr->QTensorDiffType,
                                         mat_ptr->q_diffusivity[0], mat_ptr->u_qdiffusivity, NULL,
                                         model_name, SCALAR_INPUT, &species_no, es);

          if (model_read == -1 && !strcmp(model_name, "ANISOTROPIC")) {

            mat_ptr->QTensorDiffType[species_no] = ANISOTROPIC;
            num_const = read_constants(imp, mat_ptr->u_qdiffusivity, species_no);
            if (num_const < 3) {
              sr = sprintf(err_msg, "Matl %s %s needs 3 constants for %s  model.\n",
                           pd_glob[mn]->MaterialName, "Q Tensor Diffusivity",
                           "SUSPENSION Diffusivity");
              GOMA_EH(GOMA_ERROR, err_msg);
            }
            mat_ptr->len_u_qdiffusivity[species_no] = num_const;
            SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_qdiffusivity[species_no]);

          } else if (!strncmp(model_name, "NONE", 4)) {
            species_no = atoi(strtok(NULL, " \t"));
            mat_ptr->QTensorDiffType[species_no] = NO_MODEL;
            mat_ptr->len_u_qdiffusivity[species_no] = 0;
            SPF(es, "%s = %s %d", "Q Tensor Diffusivity", "NONE", species_no);
          } else {
            sr = sprintf(err_msg, "Material %s - unrecognized model for %s \"%s\" ???\n",
                         pd_glob[mn]->MaterialName, "Q Tensor Diffusivity", model_name);
            GOMA_EH(GOMA_ERROR, err_msg);
          }

          ECHO(es, echo_file);

          species_no = mat_ptr->Num_Species; /* set species number equal to max number of species
                                                it is changed to species number of input property
                                                by look_for_mat_prop */

          model_read = look_for_mat_prop(imp, "Gravity-based Diffusivity", mat_ptr->GravDiffType,
                                         mat_ptr->g_diffusivity, mat_ptr->u_gdiffusivity, NULL,
                                         model_name, SCALAR_INPUT, &species_no, es);

          if (model_read == -1 && !strcmp(model_name, "BISECTION")) {
            mat_ptr->GravDiffType[species_no] = BISECTION;
            num_const = read_constants(imp, mat_ptr->u_gdiffusivity, species_no);
            if (num_const < 3) {
              sr = sprintf(err_msg,
                           "Matl %s %s needs 3 constants for %s  model: Kg, avg_conc, slope.\n",
                           pd_glob[mn]->MaterialName, "Gravity-based Diffusivity", "BISECTION");
              GOMA_EH(GOMA_ERROR, err_msg);
            }
            mat_ptr->len_u_gdiffusivity[species_no] = num_const;
            SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_gdiffusivity[species_no]);
          } else if (model_read == -1 && !strcmp(model_name, "RZBISECTION")) {
            mat_ptr->GravDiffType[species_no] = RZBISECTION;
            num_const = read_constants(imp, mat_ptr->u_gdiffusivity, species_no);
            if (num_const < 4) {
              sr = sprintf(
                  err_msg,
                  "Matl %s %s needs 4 constants for %s  model: Kg, exponent, avg_conc, slope .\n",
                  pd_glob[mn]->MaterialName, "Gravity-based Diffusivity", "RZBISECTION");
              GOMA_EH(GOMA_ERROR, err_msg);
            }
            mat_ptr->len_u_gdiffusivity[species_no] = num_const;
            SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_gdiffusivity[species_no]);
          } else if (model_read == -1 && !strcmp(model_name, "RICHARDSON_ZAKI")) {
            mat_ptr->GravDiffType[species_no] = RICHARDSON_ZAKI;
            num_const = read_constants(imp, mat_ptr->u_gdiffusivity, species_no);
            if (num_const < 2) {
              sr = sprintf(
                  err_msg, "Matl %s %s needs 2 constants for %s  model: Kg and exponent.\n",
                  pd_glob[mn]->MaterialName, "Gravity-based Diffusivity", "RICHARDSON_ZAKI");
              GOMA_EH(GOMA_ERROR, err_msg);
            }
            mat_ptr->len_u_gdiffusivity[species_no] = num_const;
            SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_gdiffusivity[species_no]);
          } else {
            sr = sprintf(err_msg, "Material %s - unrecognized model for %s \"%s\" ???\n",
                         pd_glob[mn]->MaterialName, "Gravity-based Diffusivity", model_name);
            GOMA_EH(model_read, err_msg);
          }

          ECHO(es, echo_file);

          species_no = mat_ptr->Num_Species; /* set species number equal to max number of species
                                                it is changed to species number of input property
                                                by look_for_mat_prop */

          model_read = look_for_mat_prop(imp, "Suspension Balance Length Scales", mat_ptr->SBM_Type,
                                         mat_ptr->SBM_Lengths, mat_ptr->SBM_Lengths2, NULL,
                                         model_name, SCALAR_INPUT, &species_no, es);

          mat_ptr->SBM_Type[species_no] = CONSTANT;
          mat_ptr->SBM_Length_enabled = 1;
          if (model_read == -1 && !strcmp(model_name, "INPUT")) {
            num_const = read_constants(imp, mat_ptr->SBM_Lengths2, species_no);
            if (num_const < 3) {
              sr = sprintf(err_msg,
                           "Matl %s %s needs 3 constants: particle radius, characteristic length "
                           "scale, and max velocity.\n",
                           pd_glob[mn]->MaterialName, "Suspension Balance Length Scales");
              GOMA_EH(GOMA_ERROR, err_msg);
            }
            mat_ptr->len_SBM_Lengths2[species_no] = num_const;
            SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->SBM_Lengths2[species_no]);

            ECHO(es, echo_file);
          } else {
            mat_ptr->SBM_Length_enabled = 0;
          }
        } else if (!strcmp(model_name, "FREE_VOL")) {
          mat_ptr->DiffusivityModel[species_no] = FREE_VOL;
          num_const = mat_ptr->len_u_diffusivity[species_no];
          if (num_const < 12) {
            sr =
                sprintf(err_msg, "Matl %s (conc %d) needs at least 12 constants for %s %s model.\n",
                        pd_glob[mn]->MaterialName, species_no, "Diffusivity", "FREE_VOL");
            GOMA_EH(GOMA_ERROR, err_msg);
          } else if (num_const < 13) {
            sr = sprintf(
                err_msg,
                "Matl %s (conc %d) needs u_diffusivity[w][12]=0 for Duda-Vrentas %s %s model.\n",
                pd_glob[mn]->MaterialName, species_no, "Diffusivity", "FREE_VOL");
            GOMA_EH(GOMA_ERROR, err_msg);
          } else if (num_const < 15 && (int)mat_ptr->u_diffusivity[species_no][12] != 0) {
            sr = sprintf(
                err_msg,
                "Matl %s (conc %d) needs at least 16 constants for non-Duda-Vrentas %s %s model.\n",
                pd_glob[mn]->MaterialName, species_no, "Diffusivity", "FREE_VOL");
            GOMA_EH(GOMA_ERROR, err_msg);
          } else if (num_const < 17 && (int)mat_ptr->u_diffusivity[species_no][12] == 4) {
            sr = sprintf(
                err_msg,
                "Matl %s (conc %d) needs at least 18 constants for friction-based %s %s model.\n",
                pd_glob[mn]->MaterialName, species_no, "Diffusivity", "FREE_VOL");
            GOMA_EH(GOMA_ERROR, err_msg);
          }

          switch ((int)mat_ptr->u_diffusivity[species_no][12]) {
          case 0:
            sr = SPF(err_msg, "\t## Matl %s (conc %d) using Vrentas & Duda %s %s model.##",
                     pd_glob[mn]->MaterialName, species_no, "Diffusivity", "FREE_VOL");
            ECHO(err_msg, echo_file);
            break;
          case 1:
            sr = SPF(err_msg, "\t ### Matl %s (conc %d) using Zelinsky & Hanley %s %s model.###\n",
                     pd_glob[mn]->MaterialName, species_no, "Diffusivity", "FREE_VOL");
            ECHO(err_msg, echo_file);
            break;
          case 2:
            sr = SPF(err_msg, "\t ### Matl %s (conc %d) using self-diffusion %s %s model.###\n",
                     pd_glob[mn]->MaterialName, species_no, "Diffusivity", "FREE_VOL");
            ECHO(err_msg, echo_file);
            break;
          case 3:
            sr = SPF(err_msg, "\t ### Matl %s (conc %d) using Alsoy & Duda %s %s model.###\n",
                     pd_glob[mn]->MaterialName, species_no, "Diffusivity", "FREE_VOL");
            ECHO(err_msg, echo_file);
            break;
          case 4:
            sr = SPF(err_msg, "\t ### Matl %s (conc %d) using friction-based %s %s model.### \n",
                     pd_glob[mn]->MaterialName, species_no, "Diffusivity", "FREE_VOL");
            ECHO(err_msg, echo_file);
            break;
          default:
            sr = SPF(err_msg, "\t ### Undefined %s %s model for Matl %s (conc %d) model.### \n",
                     "FREE_VOL", "Diffusivity", pd_glob[mn]->MaterialName, species_no);
            ECHO(err_msg, echo_file);
            break;
          }

          num_const = mat_ptr->Num_Species;
          mat_ptr->FreeVolSolvent[species_no] = TRUE;
          model_read = look_for_mat_prop(imp, "Free Volume Solvent", NULL, NULL, NO_USER, NULL,
                                         model_name, SCALAR_INPUT, &num_const, es);
          if (model_read == -1) {
            if (!strcasecmp(model_name, "yes") || !strcasecmp(model_name, "true")) {
              mat_ptr->FreeVolSolvent[num_const] = TRUE;
            } else if (!strcasecmp(model_name, "no") || !strcasecmp(model_name, "false")) {
              mat_ptr->FreeVolSolvent[num_const] = FALSE;
            } else {
              GOMA_WH(-1, "Defaulting Free Volume Solvent to TRUE");
              mat_ptr->FreeVolSolvent[species_no] = TRUE;
            }
          }
          fprintf(stderr, "Free Volume Solvent %d = %d\n", species_no,
                  mat_ptr->FreeVolSolvent[species_no]);
          /*		 if( pd_glob[mn]->Num_Species_Eqn != 1 )
                             {
                               GOMA_EH(GOMA_ERROR, "Binary Free volume models are for 2 components,
             or one tracked species.");
                             }   */
        }

        /* multicomponent, generalized fickian based formulation. */

        else if (!strcmp(model_name, "GENERALIZED_FREE_VOL")) {
          mat_ptr->DiffusivityModel[species_no] = GENERALIZED_FREE_VOL;
          num_const = mat_ptr->len_u_diffusivity[species_no];
          if (num_const < 12) {
            sr = sprintf(
                err_msg, "Matl %s (conc %d) needs at least 12 constants for %s %s model.\n",
                pd_glob[mn]->MaterialName, species_no, "Diffusivity", "GENERALIZED_FREE_VOL");
            GOMA_EH(GOMA_ERROR, err_msg);
          }
          if (pd_glob[mn]->Num_Species_Eqn < 1) {
            GOMA_EH(GOMA_ERROR, "Generalized models are for 2 or more BULK components.");
          }
        }
        /* Set a constant binary diffusivity if no concentration dependency
         * is known for generalized_fickian formulation.  It is known that
         * this is a poor approximation for multicomponent case.*/

        else if (!strcmp(model_name, "GENERALIZED")) {
          mat_ptr->DiffusivityModel[species_no] = GENERALIZED;
          num_const = mat_ptr->len_u_diffusivity[species_no];
          if (num_const < pd_glob[mn]->Num_Species_Eqn) {
            sr = sprintf(err_msg,
                         "Matl %s (conc %d) needs one constant for each i-j pair %s %s model.\n",
                         pd_glob[mn]->MaterialName, species_no, "Diffusivity", "GENERALIZED");
            GOMA_EH(GOMA_ERROR, err_msg);
          }
          if (pd_glob[mn]->Num_Species_Eqn < 2) {
            GOMA_EH(GOMA_ERROR, "Generalized diffusivity model is for 2 or more BULK components.");
          }
        } else if (!strcmp(model_name, "CONST_LS")) {
          mat_ptr->DiffusivityModel[species_no] = LEVEL_SET;
          num_const = mat_ptr->len_u_diffusivity[species_no];
        } else if (!strcmp(model_name, "CHAPMAN_GAS")) {
          mat_ptr->DiffusivityModel[species_no] = CHAPMAN_GAS;
          num_const = mat_ptr->len_u_diffusivity[species_no];
        } else {
          sprintf(err_msg, "Diffusivity. Unrecognized Diffusivity Model: %s", model_name);
          GOMA_EH(GOMA_ERROR, err_msg);
        }

      } /* End of if (model_read == 0) */

    } /*end of if(DiffusionConstitutiveEquation != STEFAN_MAXWELL &&
        DiffusionConstitutiveEquation != STEFAN_MAXWELL_CHARGED &&
        DiffusionConstitutiveEquation != STEFAN_MAXWELL_VOLUME) */

    species_no = mat_ptr->Num_Species; /* set species number equal to max number of species
                                          it is changed to species number of input property
                                          by look_for_mat_prop */

    int SpeciesSecondLevelSetDiffusivity;
    model_read = look_for_mat_prop(imp, "Species Second Level Set Diffusivity",
                                   &(SpeciesSecondLevelSetDiffusivity),
                                   mat_ptr->SpeciesSecondLevelSetDiffusivity, NO_USER, NULL,
                                   model_name, SCALAR_INPUT, &species_no, es);

    ECHO(es, echo_file);

    species_no = mat_ptr->Num_Species; /* set species number equal to max number of species
                                          it is changed to species number of input property
                                          by look_for_mat_prop */

    int SpeciesOnlyDiffusion;
    model_read = look_for_mat_prop(imp, "Species Level Set Diffusion Only", &(SpeciesOnlyDiffusion),
                                   &(a0), NO_USER, NULL, model_name, NO_INPUT, &species_no, es);

    mat_ptr->SpeciesOnlyDiffusion[species_no] = DIFF_OFF;
    if (model_read == -1 && !strcmp(model_name, "POSITIVE")) {
      mat_ptr->SpeciesOnlyDiffusion[species_no] = DIFF_POSITIVE;
    } else if (model_read == -1 && !strcmp(model_name, "NEGATIVE")) {
      mat_ptr->SpeciesOnlyDiffusion[species_no] = DIFF_NEGATIVE;
    }
    ECHO(es, echo_file);

    species_no = mat_ptr->Num_Species; /* set species number equal to max number of species
                                          it is changed to species number of input property
                                          by look_for_mat_prop */

    model_read = look_for_mat_prop(imp, "Species Time Integration", &(SpeciesTimeIntegration),
                                   &(a0), NO_USER, NULL, model_name, NO_INPUT, &species_no, es);
    if (model_read == -1 && !strcmp(model_name, "STANDARD")) {
      mat_ptr->SpeciesTimeIntegration[species_no] = STANDARD;
    } else if (model_read == -1 && !strcmp(model_name, "TAYLOR_GALERKIN")) {
      mat_ptr->SpeciesTimeIntegration[species_no] = TAYLOR_GALERKIN;
    } else if (model_read == -1 && !strcmp(model_name, "TAYLOR_GALERKIN_EXP")) {
      mat_ptr->SpeciesTimeIntegration[species_no] = TAYLOR_GALERKIN_EXP;
    } else {
      mat_ptr->SpeciesTimeIntegration[species_no] = STANDARD;
      SPF(es, "\t(%s = %s)", "Species Time Integration", "STANDARD");
    }

    ECHO(es, echo_file);

    model_read = look_for_mat_prop(imp, "Species Enable Div Term", &(SpeciesTimeIntegration), &(a0),
                                   NO_USER, NULL, model_name, NO_INPUT, &species_no, es);
    if (model_read == -1 && !strcmp(model_name, "on")) {
      mat_ptr->ExtrinsicIndependentSpeciesVar[species_no] = 1;
    } else if (model_read == -1 && !strcmp(model_name, "yes")) {
      mat_ptr->ExtrinsicIndependentSpeciesVar[species_no] = 1;
    } else {
      mat_ptr->ExtrinsicIndependentSpeciesVar[species_no] = 0;
      SPF(es, "\t(%s = %s)", "Species Enable Div Term", "off");
    }

    ECHO(es, echo_file);

    /* initialize to good default behavior */
    mat_ptr->AdvectiveScalingModel[species_no] = CONSTANT;
    mat_ptr->AdvectiveScaling[species_no] = 1.0;

    species_no = mat_ptr->Num_Species; /* set species number equal to max number of species
                                          it is changed to species number of input property
                                          by look_for_mat_prop */
    model_read = look_for_mat_prop(imp, "Advective Scaling", mat_ptr->AdvectiveScalingModel,
                                   mat_ptr->AdvectiveScaling, NO_USER, NULL, model_name,
                                   SCALAR_INPUT, &species_no, es);
    ECHO(es, echo_file);

    /*
     *  Latent Heat of Vaporization Section
     *
     *    set species number equal to max number of species
     *    it is changed to species number of input property  by look_for_mat_prop
     */
    species_no = mat_ptr->Num_Species;
    model_read = look_for_mat_prop(imp, "Latent Heat Vaporization", mat_ptr->LatentHeatVapModel,
                                   mat_ptr->latent_heat_vap, NO_USER, NULL, model_name,
                                   SCALAR_INPUT, &species_no, es);
    fallback_chemkin_generic_prop(&model_read, j, &(mat_ptr->LatentHeatVapModel[j]), TRUE, mat_ptr);
    GOMA_EH(model_read, "Latent Heat Vaporization");

    ECHO(es, echo_file);

    /*
     *  Latent Heat of Fusion Section
     *
     * set species number equal to max number of species
     * it is changed to species number of input property  by look_for_mat_prop
     */
    species_no = mat_ptr->Num_Species;
    model_read = look_for_mat_prop(imp, "Latent Heat Fusion", mat_ptr->LatentHeatFusionModel,
                                   mat_ptr->latent_heat_fusion, NO_USER, NULL, model_name,
                                   SCALAR_INPUT, &species_no, es);
    fallback_chemkin_generic_prop(&model_read, j, &(mat_ptr->LatentHeatFusionModel[j]), TRUE,
                                  mat_ptr);
    GOMA_EH(model_read, "Latent Heat Fusion");

    ECHO(es, echo_file);

    /*
     *  Vapor Pressure Section
     *
     * set species number equal to max number of species
     * it is changed to species number of input property  by look_for_mat_prop
     */
    species_no = mat_ptr->Num_Species;
    model_read = look_for_mat_prop(imp, "Vapor Pressure", mat_ptr->VaporPressureModel,
                                   mat_ptr->vapor_pressure, NO_USER, NULL, model_name, SCALAR_INPUT,
                                   &species_no, es);
    fallback_chemkin_generic_prop(&model_read, j, &(mat_ptr->VaporPressureModel[j]), TRUE, mat_ptr);
    if (model_read == -1 && (!strcmp(model_name, "KELVIN") || !strcmp(model_name, "FLAT"))) {
      if (!strcmp(model_name, "KELVIN"))
        mat_ptr->VaporPressureModel[species_no] = KELVIN;
      if (!strcmp(model_name, "FLAT"))
        mat_ptr->VaporPressureModel[species_no] = FLAT;

      if (mat_ptr->PorousMediaType == POROUS_TWO_PHASE) {
        num_const = read_constants(imp, mat_ptr->u_vapor_pressure, species_no);
        if (num_const < 5) {
          sr = sprintf(err_msg, "Matl %s (%s, conc %d) needs 5 constants for %s %s model.\n",
                       pd_glob[mn]->MaterialName, "porous 2-phase", species_no, "Vapor Pressure",
                       "KELVIN or FLAT");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
        mat_ptr->len_u_vapor_pressure[species_no] = num_const;

        SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_vapor_pressure[species_no]);

      } else if (mat_ptr->PorousMediaType == POROUS_UNSATURATED) {
        num_const = read_constants(imp, mat_ptr->u_vapor_pressure, species_no);
        if (num_const < 7) {
          sr = sprintf(err_msg, "Matl %s (%s, conc %d) needs 7 constants for %s %s model.\n",
                       pd_glob[mn]->MaterialName, "porous unsaturated", species_no,
                       "Vapor Pressure", "KELVIN OR FLAT");
          GOMA_EH(GOMA_ERROR, err_msg);
        }
        mat_ptr->len_u_vapor_pressure[species_no] = num_const;
        SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_vapor_pressure[species_no]);
      } else {
        sr = sprintf(err_msg, "%s model invalid in matl %s unless %s is \"%s\" or \"%s\"\n",
                     "KELVIN or FLAT", pd_glob[mn]->MaterialName, "Media Type", "POROUS_TWO_PHASE",
                     "POROUS_UNSATURATED");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
    } else if (model_read == -1 && !strcmp(model_name, "IDEAL_GAS")) {
      mat_ptr->VaporPressureModel[species_no] = IDEAL_GAS;

      num_const = read_constants(imp, mat_ptr->u_vapor_pressure, species_no);
      if (num_const < 3) {
        sr = sprintf(err_msg, "Matl %s (conc %d) needs 3 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, species_no, "Vapor Pressure", "IDEAL_GAS");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_vapor_pressure[species_no] = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_vapor_pressure[species_no]);

    } else if (model_read == -1 && !strcmp(model_name, "NON_VOLATILE")) {
      mat_ptr->VaporPressureModel[species_no] = NON_VOLATILE;
    } else if (model_read == -1 && !strcmp(model_name, "ANTOINE")) {
      mat_ptr->VaporPressureModel[species_no] = ANTOINE;

      num_const = read_constants(imp, mat_ptr->u_vapor_pressure, species_no);
      if (num_const < 6) {
        sr = sprintf(err_msg, "Matl %s (conc %d) needs 6 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, species_no, "Vapor Pressure", "ANTOINE");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_vapor_pressure[species_no] = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_vapor_pressure[species_no]);

    } else if (model_read == -1 && !strcmp(model_name, "RIEDEL")) {
      mat_ptr->VaporPressureModel[species_no] = RIEDEL;

      num_const = read_constants(imp, mat_ptr->u_vapor_pressure, species_no);
      if (num_const < 8) {
        sr = sprintf(err_msg, "Matl %s (conc %d) needs 6 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, species_no, "Vapor Pressure", "RIEDEL");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_vapor_pressure[species_no] = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_vapor_pressure[species_no]);
    } else {
      GOMA_EH(model_read, "Vapor Pressure");
    }

    ECHO(es, echo_file);

    /*
     *  Species Volume Expansion Section
     *
     * set species number equal to max number of species
     * it is changed to species number of input property  by look_for_mat_prop
     */
    model_read = look_for_species_prop(imp, "Species Volume Expansion", mat_ptr,
                                       mat_ptr->SpecVolExpModel, mat_ptr->species_vol_expansion,
                                       NO_USER, NULL, model_name, SCALAR_INPUT, &species_no, es);
    fallback_chemkin_generic_prop(&model_read, j, &(mat_ptr->SpecVolExpModel[j]), TRUE, mat_ptr);
    GOMA_EH(model_read, "Species Volume Expansion");

    ECHO(es, echo_file);

    /*
     * Specification of the Standard State Chemical Potential of the species
     * in the current material
     *   (mu_o(T) = This is a function of temperature only)
     */
    model_read = look_for_species_prop(imp, "Standard State Chemical Potential", mat_ptr,
                                       mat_ptr->SSChemPotModel, mat_ptr->SSChemPotData, NO_USER,
                                       NULL, model_name, SCALAR_INPUT, &species_no, es);
    fallback_chemkin_generic_prop(&model_read, j, &(mat_ptr->SSChemPotModel[j]), TRUE, mat_ptr);
    if (model_read == -1) {
      SPF(es, "\t(Standard State Chemical Potential defaulted)");
    }

    ECHO(es, echo_file);

    /*
     * Specification of the Chemical Potential of the species in the
     * solution
     */
    model_read = look_for_species_prop(imp, "Pure Species Chemical Potential", mat_ptr,
                                       mat_ptr->PSChemPotModel, mat_ptr->PSChemPotData, NO_USER,
                                       NULL, model_name, SCALAR_INPUT, &species_no, es);
    fallback_chemkin_generic_prop(&model_read, j, &(mat_ptr->PSChemPotModel[j]), TRUE, mat_ptr);
    if (model_read == 0) {
      if (!strcasecmp(model_name, "PRESSURE_INDEPENDENT")) {
        mat_ptr->PSChemPotModel[species_no] = PSCHEMPOT_PRESSURE_INDEPENDENT;
      } else if (!strcasecmp(model_name, "PRESSURE_IDEALGAS")) {
        mat_ptr->PSChemPotModel[species_no] = PSCHEMPOT_IDEALGAS;
      } else {
        GOMA_EH(GOMA_ERROR, "PSChemPot");
      }
    }

    if (model_read == -1) {
      SPF(es, "\t(Pure Species Chemical Potential defaulted)");
    }

    ECHO(es, echo_file);

    /*
     * Specification of the Chemical Potential of the species in the
     * solution
     */
    model_read = look_for_species_prop(imp, "Chemical Potential", mat_ptr, mat_ptr->ChemPotModel,
                                       mat_ptr->ChemPotData, NO_USER, NULL, model_name,
                                       SCALAR_INPUT, &species_no, es);
    fallback_chemkin_generic_prop(&model_read, j, &(mat_ptr->ChemPotModel[j]), TRUE, mat_ptr);
    if (model_read == 0) {
      if (!strcasecmp(model_name, "IDEAL_SOLUTION")) {
        mat_ptr->PSChemPotModel[species_no] = CHEMPOT_IDEALSOLN;
      } else if (!strcasecmp(model_name, "STOICHIOMETRIC_PHASE")) {
        mat_ptr->PSChemPotModel[species_no] = CHEMPOT_STOICHPHASE;
      } else {
        GOMA_EH(GOMA_ERROR, "ChemPot");
      }
    }
    if (model_read == -1) {
      SPF(es, "\t(Chemical Potential defaulted)");
    }

    ECHO(es, echo_file);

    /*
     *  Reference Concentration Section
     *
     * set species number equal to max number of species
     * it is changed to species number of input property  by look_for_mat_prop
     */
    model_read = look_for_species_prop(imp, "Reference Concentration", mat_ptr,
                                       mat_ptr->RefConcnModel, mat_ptr->reference_concn,
                                       mat_ptr->u_reference_concn, mat_ptr->len_u_reference_concn,
                                       model_name, SCALAR_INPUT, &species_no, es);

    fallback_chemkin_generic_prop(&model_read, j, &(mat_ptr->RefConcnModel[j]), TRUE, mat_ptr);
    if (mat_ptr->len_u_reference_concn[species_no] > 0)
      mat_ptr->reference_concn[species_no] = mat_ptr->u_reference_concn[species_no][0];

    GOMA_EH(model_read, "Reference Concentration");

    ECHO(es, echo_file);

    /*
     *                  Molecular Weight and Molar Volumes section
     *
     *
     * Molecular Weight and Molar Volumes are required when VL_POLY_BC
     * is specified
     * Instead of reading the constants from the BC card itself, they
     * are specified in the material's file.  8/21/98. ACSun
     *
     * HKM -> There should be checks to see that that they have
     *        been specified for those problem types
     */

    model_read = look_for_species_prop(imp, "Molecular Weight", mat_ptr,
                                       mat_ptr->MolecularWeightModel, mat_ptr->molecular_weight,
                                       NO_USER, NULL, model_name, SCALAR_INPUT, &species_no, es);
    fallback_chemkin_generic_prop(&model_read, j, &(mat_ptr->MolecularWeightModel[j]), TRUE,
                                  mat_ptr);
    ECHO(es, echo_file);

    model_read = look_for_species_prop(imp, "Specific Volume", mat_ptr,
                                       mat_ptr->SpecificVolumeModel, mat_ptr->specific_volume,
                                       NO_USER, NULL, model_name, SCALAR_INPUT, &species_no, es);
    fallback_chemkin_generic_prop(&model_read, j, &(mat_ptr->SpecificVolumeModel[j]), TRUE,
                                  mat_ptr);
    ECHO(es, echo_file);

    model_read = look_for_species_prop(imp, "Molar Volume", mat_ptr, mat_ptr->MolarVolumeModel,
                                       mat_ptr->molar_volume, NO_USER, NULL, model_name,
                                       SCALAR_INPUT, &species_no, es);
    fallback_chemkin_generic_prop(&model_read, j, &(mat_ptr->MolarVolumeModel[j]), TRUE, mat_ptr);
    ECHO(es, echo_file);

    model_read = look_for_species_prop(imp, "Charge Number", mat_ptr, mat_ptr->ChargeNumberModel,
                                       mat_ptr->charge_number, NO_USER, NULL, model_name,
                                       SCALAR_INPUT, &species_no, es);
    fallback_chemkin_generic_prop(&model_read, j, &(mat_ptr->ChargeNumberModel[j]), TRUE, mat_ptr);
    ECHO(es, echo_file);

  } /* End for(j=0;j<mat_ptr->Num_Species; j++) */

  if (DiffusionConstitutiveEquation == FICKIAN_CHARGED ||
      DiffusionConstitutiveEquation == FICKIAN_CHARGED_X) /*  RSL 6/23/02  */
  {
    model_read = look_for_mat_prop(
        imp, "Solution Temperature", &(mat_ptr->SolutionTemperatureModel),
        &(mat_ptr->solution_temperature), NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);
    GOMA_EH(model_read,
            "Solution Temperature must be specified when FICKIAN_CHARGED flux model is used.");
    ECHO(es, echo_file);
  }

  /*
   * Forward one pass parsing is choking. Go back and look better this
   * time.
   */
  rewind(imp);

  if (read_bc_mp != -1) {
    iread = look_for_optional(imp, "Non-condensable Molecular Weight", input, '=');
    if (iread != -1) {
      if (fscanf(imp, "%s %d %lf", model_name, &ii, &mw) != 3) {
        GOMA_EH(GOMA_ERROR, "Error reading non-condensable MW: e.g. CONSTANT species_no  MW");
      } else {
        mat_ptr->molecular_weight[mat_ptr->Num_Species_Eqn] = mw;
        SPF(es, "%s = %s %d %.4g", "Non-condensable Molecular Weight", model_name, ii, mw);
      }
      ECHO(es, echo_file);
    }

    iread = look_for_optional(imp, "Non-volatile Molar Volume", input, '=');
    if (iread != -1) {
      if (fscanf(imp, "%s %d %lf", model_name, &ii, &mv) != 3) {
        GOMA_EH(GOMA_ERROR,
                "Error reading non-volatile Molar Volume: e.g. CONSTANT  species_id  MV");
      } else {
        mat_ptr->molar_volume[pd_glob[mn]->Num_Species_Eqn] = mv;
        SPF(es, "%s = %s %d %.4g", "Non-volatile Molar Volume", model_name, ii, mw);
      }
      ECHO(es, echo_file);
    }

    iread = look_for_optional(imp, "Non-volatile Specific Volume", input, '=');
    if (iread != -1) {
      if (fscanf(imp, "%s %d %lf", model_name, &ii, &mv) != 3) {
        GOMA_EH(GOMA_ERROR,
                "Error reading non-volatile Specific Volume: e.g. CONSTANT  species_id  MV");
      } else {
        mat_ptr->specific_volume[pd_glob[mn]->Num_Species_Eqn] = mv;
        SPF(es, "%s = %s %d %.4g", "Non-volatile Specific Volume", model_name, ii, mw);
      }
      ECHO(es, echo_file);
    }

    iread = look_for_optional(imp, "Flory-Huggins parameters", input, '=');
    if (iread != -1) {
      n_species = pd_glob[mn]->Num_Species_Eqn + 1;
      /*number of independent interaction parameters */
      n_ij = (n_species * n_species - n_species) / 2;

      if (fscanf(imp, "%s", model_name) != 1) {
        GOMA_EH(GOMA_ERROR, "Error reading F-H parameter model name: e.g. CONSTANT");
      } else {
        for (i = 0; i < n_ij; i++) /* reading the chi parameters */
        {
          if (fscanf(imp, "%d %d %lf", &ii, &jj, &chi_ij) != 3) {
            GOMA_EH(GOMA_ERROR, "Error:must have three entries, i, j, and chi(i,j)");
          }
          mat_ptr->flory_param[ii][jj] = chi_ij;
          mat_ptr->flory_param[jj][ii] =
              chi_ij * mat_ptr->molar_volume[jj] / mat_ptr->molar_volume[ii];
        }
        for (k = 0; k < n_species; k++) {
          mat_ptr->flory_param[k][k] = 0.;
        }
        SPF(es, "%s = %s %d %d %.4g", "Flory-Huggins parameters", model_name, ii, jj, chi_ij);
      }
      ECHO(es, echo_file);
    }
  }

  /*
   * Moment Properties
   */

  if (pd_glob[mn]->gv[MOMENT0]) {
    ECHO("\n----Moment Properties\n", echo_file);

    model_read = look_for_mat_prop(imp, "Moment Weight Function", &(mat_ptr->Momentwt_funcModel),
                                   &(mat_ptr->Momentwt_func), NO_USER, NULL, model_name,
                                   SCALAR_INPUT, &NO_SPECIES, es);

    if (!strcmp(model_name, "SUPG")) {
      mat_ptr->Momentwt_funcModel = SUPG;
      if (fscanf(imp, "%lg", &(mat_ptr->Momentwt_func)) != 1) {
        GOMA_EH(GOMA_ERROR, "Could not read SUPG value for Moment Weight Function");
      }
      SPF(endofstring(es), "SUPG %.4g", mat_ptr->Mwt_func);
    } else {
      mat_ptr->Momentwt_funcModel = GALERKIN;
      mat_ptr->Momentwt_func = 0.;
      SPF(es, "\t(%s = %s)", "Moment Weight Function", "GALERKIN");
    }
    ECHO(es, echo_file);

    model_read = look_for_mat_prop(imp, "Moment SSPG Function", &(mat_ptr->MomentSSPG_funcModel),
                                   &(mat_ptr->MomentSSPG_func), NO_USER, NULL, model_name,
                                   SCALAR_INPUT, &NO_SPECIES, es);
    ECHO(es, echo_file);

    model_read = look_for_mat_prop(imp, "Moment Diffusivity", &(mat_ptr->MomentDiffusivityModel),
                                   &(mat_ptr->MomentDiffusivity), NO_USER, NULL, model_name,
                                   SCALAR_INPUT, &NO_SPECIES, es);
    ECHO(es, echo_file);

    model_read = look_for_mat_prop(imp, "Moment Second Level Set Diffusivity",
                                   &(mat_ptr->MomentSecondLevelSetDiffusivityModel),
                                   &(mat_ptr->MomentSecondLevelSetDiffusivity), NO_USER, NULL,
                                   model_name, SCALAR_INPUT, &NO_SPECIES, es);

    ECHO(es, echo_file);

    model_read = look_for_mat_prop(imp, "Moment Level Set Diffusion Only",
                                   &(mat_ptr->MomentLevelSetDiffusionOnly), &(a0), NO_USER, NULL,
                                   model_name, SCALAR_INPUT, &NO_SPECIES, es);

    if (!strcmp(model_name, "POSITIVE")) {
      mat_ptr->MomentLevelSetDiffusionOnly = DIFF_POSITIVE;
      SPF(endofstring(es), "POSITIVE");
    } else if (!strcmp(model_name, "NEGATIVE")) {
      mat_ptr->MomentLevelSetDiffusionOnly = DIFF_NEGATIVE;
      SPF(endofstring(es), "POSITIVE");
    } else {
      mat_ptr->MomentLevelSetDiffusionOnly = DIFF_OFF;
      SPF(endofstring(es), "OFF");
    }

    ECHO(es, echo_file);

    model_read = look_for_mat_prop(imp, "Moment Shock Function", &(mat_ptr->MomentShock_funcModel),
                                   &(mat_ptr->MomentShock_func), NO_USER, NULL, model_name,
                                   SCALAR_INPUT, &NO_SPECIES, es);

    if (!strcmp(model_name, "YZBETA")) {
      mat_ptr->MomentShock_funcModel = YZBETA_MIXED;
      if (fscanf(imp, "%lg %lg %lg %lg %lg", &(mat_ptr->MomentShock_func), &a0, &a1, &a2, &a3) !=
          5) {
        GOMA_EH(GOMA_ERROR,
                "Could not read YZbeta value for Moment Shock Function, expected 5 values");
      }
      mat_ptr->MomentShock_Ref[0] = a0;
      mat_ptr->MomentShock_Ref[1] = a1;
      mat_ptr->MomentShock_Ref[2] = a2;
      mat_ptr->MomentShock_Ref[3] = a3;
      SPF(endofstring(es), "YZBETA %.4g", mat_ptr->Mwt_func);
    } else {
      mat_ptr->MomentShock_funcModel = SC_NONE;
      mat_ptr->MomentShock_func = 0.;
      SPF(es, "\t(%s = %s)", "Moment Shock Function", "NONE");
    }
    ECHO(es, echo_file);

    model_read = look_for_mat_prop(imp, "Moment Growth Kernel", &(mat_ptr->moment_growth_model),
                                   &(mat_ptr->moment_growth_scale), NULL, NULL, model_name,
                                   SCALAR_INPUT, &NO_SPECIES, es);
    if (!strcmp(model_name, "VISCOSITY_SCALED")) {
      model_read = 1;
      mat_ptr->moment_growth_model = VISCOSITY_SCALED_GROWTH_RATE;
      if (fscanf(imp, "%lf", &a0) != 1) {
        sr = sprintf(err_msg, "Matl %s needs 1 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Moment Growth Kernel", "VISCOSITY_SCALED");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->moment_growth_scale = a0;
      SPF_DBL_VEC(endofstring(es), 1, &(mat_ptr->moment_growth_scale));
    } else if (!strcmp(model_name, "VISCOSITY_PRESSURE_SCALED")) {
      model_read = 1;
      mat_ptr->moment_growth_model = VISCOSITY_PRESSURE_GROWTH_RATE;
      if (fscanf(imp, "%lf %lf", &a0, &a1) != 2) {
        sr =
            sprintf(err_msg, "Matl %s needs 2 constants for %s %s model.\n",
                    pd_glob[mn]->MaterialName, "Moment Growth Kernel", "VISCOSITY_PRESSURE_SCALED");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->moment_growth_scale = a0;
      mat_ptr->moment_growth_reference_pressure = a1;
      SPF_DBL_VEC(endofstring(es), 1, &(mat_ptr->moment_growth_scale));
    } else {
      if (model_read == -1) {
        GOMA_EH(model_read, "Moment Growth Kernel invalid");
      }
      GOMA_EH(model_read, "Moment Growth Kernel");
    }

    ECHO(es, echo_file);

    model_read =
        look_for_mat_prop(imp, "Moment Coalescence Kernel", &(mat_ptr->moment_coalescence_model),
                          &(mat_ptr->moment_coalescence_scale), NULL, NULL, model_name,
                          SCALAR_INPUT, &NO_SPECIES, es);
    if (!strcmp(model_name, "ADDITION")) {
      model_read = 1;
      mat_ptr->moment_coalescence_model = ADDITION_COALESCENCE;
      if (fscanf(imp, "%lf", &a0) != 1) {
        sr =
            sprintf(err_msg, "Matl %s needs 1 constants for %s %s model.\n",
                    pd_glob[mn]->MaterialName, "Moment Coalescence Kernel", "ADDITION_COALESCENCE");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->moment_coalescence_scale = a0;
      SPF_DBL_VEC(endofstring(es), 1, &(mat_ptr->moment_coalescence_scale));
    } else if (!strcmp(model_name, "VISCOSITY_SCALED_ADDITION")) {
      model_read = 1;
      mat_ptr->moment_coalescence_model = VISCOSITY_SCALED_COALESCENCE;
      if (fscanf(imp, "%lf", &a0) != 1) {
        sr = sprintf(err_msg, "Matl %s needs 1 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Moment Coalescence Kernel",
                     "VISCOSITY_SCALED_ADDITION");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->moment_coalescence_scale = a0;
      SPF_DBL_VEC(endofstring(es), 1, &(mat_ptr->moment_coalescence_scale));
    } else if (!strcmp(model_name, "VISCOSITY_ADDITION_BUBBLE_RATIO")) {
      model_read = 1;
      mat_ptr->moment_coalescence_model = VISCOSITY_BUBBLE_RATIO_COALESCENCE;
      if (fscanf(imp, "%lf", &a0) != 1) {
        sr = sprintf(err_msg, "Matl %s needs 1 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Moment Coalescence Kernel",
                     "VISCOSITY_ADDITION_BUBBLE_RATIO");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->moment_coalescence_scale = a0;
      SPF_DBL_VEC(endofstring(es), 1, &(mat_ptr->moment_coalescence_scale));
    } else {
      if (model_read == -1) {
        GOMA_EH(model_read, "Moment Coalescence Kernel invalid");
      }
      GOMA_EH(model_read, "Moment Coalescence Kernel");
    }

    ECHO(es, echo_file);
  }

  /*
   * Source Terms
   */

  ECHO("\n----Volumetric Source Terms\n", echo_file);

  model_read = look_for_mat_prop(imp, "Navier-Stokes Source", &(mat_ptr->MomentumSourceModel),
                                 mat_ptr->momentum_source, &(mat_ptr->u_momentum_source),
                                 &(mat_ptr->len_u_momentum_source), model_name, VECTOR_INPUT,
                                 &NO_SPECIES, es);
  if (!strcmp(model_name, "BOUSS")) {
    MomentumSourceModel = BOUSS;
    model_read = 1;
    mat_ptr->MomentumSourceModel = MomentumSourceModel;
    if (fscanf(imp, "%lf %lf %lf", &a0, &a1, &a2) != 3) {
      sr = sprintf(err_msg, "Matl %s needs 3 constants for %s %s model.\n",
                   pd_glob[mn]->MaterialName, "Navier-Stokes Source", "BOUSS");
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    mat_ptr->momentum_source[0] = a0;
    mat_ptr->momentum_source[1] = a1;
    mat_ptr->momentum_source[2] = a2;
    SPF_DBL_VEC(endofstring(es), 3, mat_ptr->momentum_source);
  } else if (!strcmp(model_name, "VARIABLE_DENSITY")) {
    MomentumSourceModel = VARIABLE_DENSITY;
    model_read = 1;
    mat_ptr->MomentumSourceModel = MomentumSourceModel;
    if (fscanf(imp, "%lf %lf %lf", &a0, &a1, &a2) != 3) {
      sr = sprintf(err_msg, "Matl %s needs 3 constants for %s %s model.\n",
                   pd_glob[mn]->MaterialName, "Navier-Stokes Source", "VARIABLE_DENSITY");
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    mat_ptr->momentum_source[0] = a0;
    mat_ptr->momentum_source[1] = a1;
    mat_ptr->momentum_source[2] = a2;
    SPF_DBL_VEC(endofstring(es), 3, mat_ptr->momentum_source);
  } else if (!strcmp(model_name, "VARIABLE_DENSITY_NO_GAS")) {
    MomentumSourceModel = VARIABLE_DENSITY_NO_GAS;
    model_read = 1;
    mat_ptr->MomentumSourceModel = MomentumSourceModel;
    if (fscanf(imp, "%lf %lf %lf", &a0, &a1, &a2) != 3) {
      sr = sprintf(err_msg, "Matl %s needs 3 constants for %s %s model.\n",
                   pd_glob[mn]->MaterialName, "Navier-Stokes Source", "VARIABLE_DENSITY");
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    mat_ptr->momentum_source[0] = a0;
    mat_ptr->momentum_source[1] = a1;
    mat_ptr->momentum_source[2] = a2;
    SPF_DBL_VEC(endofstring(es), 3, mat_ptr->momentum_source);
  } else if (!strcmp(model_name, "GRAV_VIBRATIONAL")) {
    MomentumSourceModel = GRAV_VIBRATIONAL;
    model_read = 1;
    mat_ptr->MomentumSourceModel = MomentumSourceModel;
    if (fscanf(imp, "%lf %lf %lf %lf %lf", &a0, &a1, &a2, &a3, &a4) != 5) {
      sr = sprintf(err_msg, "Matl %s needs 5 constants for %s %s model.\n",
                   pd_glob[mn]->MaterialName, "Navier-Stokes Source", "GRAV_VIBRATIONAL");
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    mat_ptr->u_momentum_source = (dbl *)array_alloc(1, 2, sizeof(dbl));
    mat_ptr->len_u_momentum_source = 2;

    mat_ptr->momentum_source[0] = a0;
    mat_ptr->momentum_source[1] = a1;
    mat_ptr->momentum_source[2] = a2;
    mat_ptr->u_momentum_source[0] = a3;
    mat_ptr->u_momentum_source[1] = a4;
    SPF_DBL_VEC(endofstring(es), 3, mat_ptr->momentum_source);
    SPF_DBL_VEC(endofstring(es), 2, mat_ptr->u_momentum_source);
  } else if (!strcmp(model_name, "FILL")) {
    MomentumSourceModel = FILL_SRC;
    model_read = 1;
    mat_ptr->MomentumSourceModel = MomentumSourceModel;
    if (fscanf(imp, "%lf %lf %lf", &a0, &a1, &a2) != 3) {
      sr = sprintf(err_msg, "Matl %s needs 3 constants for %s %s model.\n",
                   pd_glob[mn]->MaterialName, "Navier-Stokes Source", "FILL");
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    mat_ptr->momentum_source[0] = a0;
    mat_ptr->momentum_source[1] = a1;
    mat_ptr->momentum_source[2] = a2;
    SPF_DBL_VEC(endofstring(es), 3, mat_ptr->momentum_source);
  } else if (!strcmp(model_name, "LEVEL_SET") || !strcmp(model_name, "PHASE_FUNCTION")) {
    MomentumSourceModel = LEVEL_SET;
    model_read = 1;
    mat_ptr->MomentumSourceModel = MomentumSourceModel;
    if (fscanf(imp, "%lf %lf %lf", &a0, &a1, &a2) != 3) {
      sr = sprintf(err_msg, "Matl %s needs 3 constants for %s %s model.\n",
                   pd_glob[mn]->MaterialName, "Navier-Stokes Source", "LEVEL_SET");
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    mat_ptr->momentum_source[0] = a0;
    mat_ptr->momentum_source[1] = a1;
    mat_ptr->momentum_source[2] = a2;
    SPF_DBL_VEC(endofstring(es), 3, mat_ptr->momentum_source);
  } else if (!strcmp(model_name, "SUSPEND")) {
    MomentumSourceModel = SUSPEND;
    model_read = 1;
    mat_ptr->MomentumSourceModel = MomentumSourceModel;
    if (fscanf(imp, "%lf %lf %lf %lf", &a0, &a1, &a2, &a3) != 4) {
      sr = sprintf(err_msg, "Matl %s needs 4 constants for %s %s model.\n",
                   pd_glob[mn]->MaterialName, "Navier-Stokes Source", "SUSPEND");
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    mat_ptr->u_momentum_source = (dbl *)array_alloc(1, 1, sizeof(dbl));

    mat_ptr->len_u_momentum_source = 1;

    mat_ptr->momentum_source[0] = a0; /* x component of gravity */
    mat_ptr->momentum_source[1] = a1; /* y component of gravity */
    mat_ptr->momentum_source[2] = a2; /* z component of gravity */

    mat_ptr->u_momentum_source[0] = a3;

    /* fluid and solid densities are specified with the SUSPENSION
       model on the density card */
    if (mat_ptr->DensityModel != SUSPENSION) {
      sr = sprintf(err_msg, "For matl %s, %s = \"%s\" needs %s = \"%s\".\n",
                   pd_glob[mn]->MaterialName, "Navier-Stokes Source", "SUSPEND", "Density",
                   "SUSPENSION");
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    SPF_DBL_VEC(endofstring(es), 3, mat_ptr->momentum_source);
    SPF_DBL_VEC(endofstring(es), 1, mat_ptr->u_momentum_source);
  } else if (!strcmp(model_name, "SUSPENSION")) {
    MomentumSourceModel = SUSPENSION;
    model_read = 1;
    mat_ptr->MomentumSourceModel = MomentumSourceModel;
    if (fscanf(imp, "%lf %lf %lf %lf", &a0, &a1, &a2, &a3) != 4) {
      sr = sprintf(err_msg, "Matl %s needs 4 constants for %s %s model.\n",
                   pd_glob[mn]->MaterialName, "Navier-Stokes Source", "SUSPENSION");
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    mat_ptr->u_momentum_source = (dbl *)array_alloc(1, 1, sizeof(dbl));

    mat_ptr->len_u_momentum_source = 1;

    mat_ptr->momentum_source[0] = a0; /* x component of gravitation */
    mat_ptr->momentum_source[1] = a1; /* y component of gravitation */
    mat_ptr->momentum_source[2] = a2; /* z component of gravitation */

    mat_ptr->u_momentum_source[0] = a3;

    /* fluid and solid densities are specified with the SUSPENSION
       model on the density card */
    if (mat_ptr->DensityModel != SUSPENSION) {
      sr = sprintf(err_msg, "For matl %s, %s = \"%s\" needs %s = \"%s\".\n",
                   pd_glob[mn]->MaterialName, "Navier-Stokes Source", "SUSPEND", "Density",
                   "SUSPENSION");
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    SPF_DBL_VEC(endofstring(es), 3, mat_ptr->momentum_source);
    SPF_DBL_VEC(endofstring(es), 1, mat_ptr->u_momentum_source);
  }
  /* MMH */
  else if (!strcmp(model_name, "SUSPENSION_PM")) {
    MomentumSourceModel = SUSPENSION_PM;
    model_read = 1;
    mat_ptr->MomentumSourceModel = MomentumSourceModel;
    if (fscanf(imp, "%lf %lf %lf", &a0, &a1, &a2) != 3) {
      sr = sprintf(err_msg, "Matl %s needs 3 constants for %s %s model.\n",
                   pd_glob[mn]->MaterialName, "Navier-Stokes Source", "SUSPENSION_PM");
      GOMA_EH(GOMA_ERROR, err_msg);
    }

    mat_ptr->momentum_source[0] = a0; /* x component of gravitation */
    mat_ptr->momentum_source[1] = a1; /* y component of gravitation */
    mat_ptr->momentum_source[2] = a2; /* z component of gravitation */

    /* For SUSPENSION_PM momentum source term, the densities are actually
     * constant within each phase.  The DensityModel card is SUSPENSION_PM,
     * which is currently just like Density = SUSPENSION.
     */
    if (mat_ptr->DensityModel != SUSPENSION_PM) {
      sprintf(err_msg, "For matl %s, %s = \"%s\" needs %s = \"%s\".\n", pd_glob[mn]->MaterialName,
              "Navier-Stokes Source", "SUSPENSION_PM", "Density", "SUSPENSION_PM");
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    SPF_DBL_VEC(endofstring(es), 3, mat_ptr->momentum_source);
  } else if (!strcmp(model_name, "BOUSS_JXB")) {
    MomentumSourceModel = BOUSS_JXB;
    model_read = 1;
    mat_ptr->MomentumSourceModel = MomentumSourceModel;
    if (fscanf(imp, "%lf %lf %lf %lf", &a0, &a1, &a2, &a3) != 4) {
      sr = sprintf(err_msg, "Matl %s needs 4 constants for %s=\"%s\" model.\n",
                   pd_glob[mn]->MaterialName, "Navier-Stokes Source", "BOUSS_JXB");
      GOMA_EH(GOMA_ERROR, err_msg);
    }

    mat_ptr->u_momentum_source = (dbl *)array_alloc(1, 1, sizeof(dbl));
    mat_ptr->len_u_momentum_source = 1;

    mat_ptr->momentum_source[0] = a0;
    mat_ptr->momentum_source[1] = a1;
    mat_ptr->momentum_source[2] = a2;

    mat_ptr->u_momentum_source[0] = a3;

    SPF_DBL_VEC(endofstring(es), 3, mat_ptr->momentum_source);
    SPF_DBL_VEC(endofstring(es), 1, mat_ptr->u_momentum_source);
  } else if (!strcmp(model_name, "BOUSSINESQ")) {
    MomentumSourceModel = BOUSSINESQ;
    model_read = 1;
    mat_ptr->MomentumSourceModel = MomentumSourceModel;
    if (fscanf(imp, "%lf %lf %lf", &a0, &a1, &a2) != 3) {
      sr = sprintf(err_msg, "Matl %s needs 3 constants for %s=\"%s\" model.\n",
                   pd_glob[mn]->MaterialName, "Navier-Stokes Source", "BOUSSINESQ");
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    /*
     * Kludge: save these 3 body force vector components as if they
     * were a momentum_source, even though later on we will construct
     * the true momentum source by using the density, the coefficient
     * of volume expansion, etc.
     */
    mat_ptr->momentum_source[0] = a0;
    mat_ptr->momentum_source[1] = a1;
    mat_ptr->momentum_source[2] = a2;
    SPF_DBL_VEC(endofstring(es), 3, mat_ptr->momentum_source);
  } else if (!strcmp(model_name, "EHD_POLARIZATION")) {
    MomentumSourceModel = EHD_POLARIZATION;
    model_read = 1;
    mat_ptr->MomentumSourceModel = MomentumSourceModel;
    if (fscanf(imp, "%lf", &a0) != 1) {
      sr = sprintf(err_msg, "Matl %s needs 1 constants for %s=\"%s\" model.\n",
                   pd_glob[mn]->MaterialName, "Navier-Stokes Source", "EHD_POLARIZATION");
      GOMA_EH(GOMA_ERROR, err_msg);
    }

    mat_ptr->momentum_source[0] = a0;
    SPF_DBL_VEC(endofstring(es), 1, mat_ptr->momentum_source);
  } else if (!strcmp(model_name, "ACOUSTIC")) {
    MomentumSourceModel = ACOUSTIC;
    model_read = 1;
    mat_ptr->MomentumSourceModel = MomentumSourceModel;
    if (fscanf(imp, "%lf %lf %lf %lf", &a0, &a1, &a2, &a3) != 4) {
      sr = sprintf(err_msg, "Matl %s needs 4 constants for %s=\"%s\" model.\n",
                   pd_glob[mn]->MaterialName, "Navier-Stokes Source", "ACOUSTIC");
      GOMA_EH(GOMA_ERROR, err_msg);
    }

    mat_ptr->u_momentum_source = (dbl *)array_alloc(1, 1, sizeof(dbl));
    mat_ptr->len_u_momentum_source = 1;

    mat_ptr->momentum_source[0] = a0;
    mat_ptr->momentum_source[1] = a1;
    mat_ptr->momentum_source[2] = a2;

    mat_ptr->u_momentum_source[0] = a3;

    SPF_DBL_VEC(endofstring(es), 3, mat_ptr->momentum_source);
    SPF_DBL_VEC(endofstring(es), 1, mat_ptr->u_momentum_source);

  } else {
    if (model_read == -1) {
      GOMA_EH(model_read, "Navier-Stokes Source model invalid");
    }
    GOMA_EH(model_read, "Navier-Stokes Source");
  }

  ECHO(es, echo_file);

  model_read = look_for_mat_prop(imp, "Second Level Set Momentum Source", &(i0), v0, NO_USER, NULL,
                                 model_name, VECTOR_INPUT, &NO_SPECIES, es);

  if (model_read != -1) {
    if (ls == NULL)
      GOMA_EH(GOMA_ERROR,
              "Second Level Set Momentum Source requires activation of Level Set Tracking.\n");

    mat_ptr->mp2nd->MomentumSourceModel = i0;
    mat_ptr->mp2nd->momentumsource[0] = v0[0];
    mat_ptr->mp2nd->momentumsource[1] = v0[1];
    mat_ptr->mp2nd->momentumsource[2] = v0[2];

    stringup(model_name);

    if (!strcmp(model_name, "CONSTANT")) {
      if (fscanf(imp, "%s", input) != 1) {
        GOMA_EH(GOMA_ERROR, "Expecting trailing keyword for Second Level Set Momentum Source.\n");
      }

      stringup(input);

      if (strncmp(input, "POSITIVE", 3) == 0) {
        mat_ptr->mp2nd->momentumsourcemask[0] = 0;
        mat_ptr->mp2nd->momentumsourcemask[1] = 1;
      } else if (strncmp(input, "NEGATIVE", 3) == 0) {
        mat_ptr->mp2nd->momentumsourcemask[0] = 1;
        mat_ptr->mp2nd->momentumsourcemask[1] = 0;
      } else {
        GOMA_EH(GOMA_ERROR,
                "Keyword must be POSITIVE or NEGATIVE for Second Level Set Momentum Source.\n");
      }
      SPF(endofstring(es), " %s", input);
      if (pfd != NULL) {
        for (i = 0; i < pfd->num_phase_funcs; i++) {
          if (fscanf(imp, "%lf %lf %lf", &(mat_ptr->mp2nd->momentumsource_phase[i][0]),
                     &(mat_ptr->mp2nd->momentumsource_phase[i][1]),
                     &(mat_ptr->mp2nd->momentumsource_phase[i][2])) != 3) {
            GOMA_EH(GOMA_ERROR, "error reading phase momentum source");
          }
          SPF(endofstring(es), " %g %g %g ", mat_ptr->mp2nd->momentumsource_phase[i][0],
              mat_ptr->mp2nd->momentumsource_phase[i][1],
              mat_ptr->mp2nd->momentumsource_phase[i][2]);
        }
      }

    } else {
      GOMA_EH(GOMA_ERROR, "Second Level Set Momentum Source model can only be CONSTANT.\n");
    }
  }

  ECHO(es, echo_file);

  model_read =
      look_for_mat_prop(imp, "Solid Body Source", &(mat_ptr->MeshSourceModel), mat_ptr->mesh_source,
                        NO_USER, NULL, model_name, VECTOR_INPUT, &NO_SPECIES, es);

  mat_ptr->RealSolidSourceModel = CONSTANT; // See total_ale discussion below
  mat_ptr->real_solid_source[0] = mat_ptr->real_solid_source[1] = mat_ptr->real_solid_source[2] =
      0.;

  ECHO(es, echo_file);

  model_read =
      look_for_mat_prop(imp, "Mass Source", &(mat_ptr->MassSourceModel), &(mat_ptr->mass_source),
                        NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);
  GOMA_EH(model_read, "Mass Source");
  ECHO(es, echo_file);

  model_read =
      look_for_mat_prop(imp, "Heat Source", &(mat_ptr->HeatSourceModel), &(mat_ptr->heat_source),
                        &(mat_ptr->u_heat_source), &(mat_ptr->len_u_heat_source), model_name,
                        SCALAR_INPUT, &NO_SPECIES, es);
  if (!strcmp(model_name, "JOULE")) {
    HeatSourceModel = JOULE;
    model_read = 1;
    mat_ptr->HeatSourceModel = HeatSourceModel;
    num_const = read_constants(imp, &(mat_ptr->u_heat_source), NO_SPECIES);
    mat_ptr->len_u_heat_source = num_const;
    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_heat_source);
  } else if (!strcmp(model_name, "VISC_DISS")) {
    HeatSourceModel = VISC_DISS;
    model_read = 1;
    mat_ptr->HeatSourceModel = HeatSourceModel;
    num_const = read_constants(imp, &(mat_ptr->u_heat_source), NO_SPECIES);
    mat_ptr->len_u_heat_source = num_const;
    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_heat_source);
  } else if (!strcmp(model_name, "DROP_EVAP")) {
    HeatSourceModel = DROP_EVAP;
    model_read = 1;
    mat_ptr->HeatSourceModel = HeatSourceModel;
    num_const = read_constants(imp, &(mat_ptr->u_heat_source), NO_SPECIES);
    mat_ptr->len_u_heat_source = num_const;
    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_heat_source);
  } else if (!strcmp(model_name, "EPOXY")) {
    HeatSourceModel = EPOXY;
    model_read = 1;
    mat_ptr->HeatSourceModel = HeatSourceModel;
    if (fscanf(imp, "%lf", &a0) != 1) {
      sprintf(err_msg, "Matl %s needs 1 constant for %s %s model.\n", pd_glob[mn]->MaterialName,
              "Heat Source", "EPOXY");
      GOMA_EH(GOMA_ERROR, err_msg);
    }

    mat_ptr->heat_source = a0; /* this will be the heat of reaction for the epoxy */
    SPF_DBL_VEC(endofstring(es), 1, &(mat_ptr->heat_source));
  } else if (!strcmp(model_name, "BUTLER_VOLMER")) /* added by KSC: 04/21/2006 */
  {
    HeatSourceModel = BUTLER_VOLMER;
    model_read = 1;
    mat_ptr->HeatSourceModel = HeatSourceModel;
    num_const = read_constants(imp, &(mat_ptr->u_heat_source), NO_SPECIES);
    mat_ptr->len_u_heat_source = num_const;
    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_heat_source);
  } else if (!strcmp(model_name, "PHOTO_CURING")) /* added by RBS: 06/05/2013 */
  {
    HeatSourceModel = PHOTO_CURING;
    model_read = 1;
    mat_ptr->HeatSourceModel = HeatSourceModel;
    num_const = read_constants(imp, &(mat_ptr->u_heat_source), NO_SPECIES);
    mat_ptr->len_u_heat_source = num_const;
    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_heat_source);
  } else if (!strcmp(model_name, "ELECTRODE_KINETICS")) /* added by KSC/GHE on 10/21/98 */
  {
    HeatSourceModel = ELECTRODE_KINETICS;
    model_read = 1;
    mat_ptr->HeatSourceModel = HeatSourceModel;
  } else if (!strcmp(model_name, "FOAM")) {
    HeatSourceModel = HS_FOAM;
    model_read = 1;
    mat_ptr->HeatSourceModel = HeatSourceModel;

    num_const = read_constants(imp, &(mat_ptr->u_heat_source), NO_SPECIES);
    mat_ptr->len_u_heat_source = num_const;
    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_heat_source);
  } else if (!strcmp(model_name, "FOAM_PMDI_10")) {
    HeatSourceModel = HS_FOAM_PMDI_10;
    model_read = 1;
    mat_ptr->HeatSourceModel = HeatSourceModel;

    num_const = read_constants(imp, &(mat_ptr->u_heat_source), NO_SPECIES);
    mat_ptr->len_u_heat_source = num_const;
    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_heat_source);
  } else if (!strcmp(model_name, "VISC_ACOUSTIC")) {
    HeatSourceModel = VISC_ACOUSTIC;
    model_read = 1;
    mat_ptr->HeatSourceModel = HeatSourceModel;
    num_const = read_constants(imp, &(mat_ptr->u_heat_source), NO_SPECIES);
    mat_ptr->len_u_heat_source = num_const;
    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_heat_source);
  } else if (!strcmp(model_name, "FOAM_PBE")) {
    HeatSourceModel = HS_FOAM_PBE;
    model_read = 1;
    mat_ptr->HeatSourceModel = HeatSourceModel;

    num_const = read_constants(imp, &(mat_ptr->u_heat_source), NO_SPECIES);
    mat_ptr->len_u_heat_source = num_const;
    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_heat_source);
  } else if (!strcmp(model_name, "EM_DISS")) {
    HeatSourceModel = EM_DISS;
    model_read = 1;
    mat_ptr->HeatSourceModel = HeatSourceModel;
    num_const = read_constants(imp, &(mat_ptr->u_heat_source), NO_SPECIES);
    mat_ptr->len_u_heat_source = num_const;
    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_heat_source);
  } else if (!strcmp(model_name, "EM_VECTOR_DISS")) {
    HeatSourceModel = EM_VECTOR_DISS;
    model_read = 1;
    mat_ptr->HeatSourceModel = HeatSourceModel;
    num_const = read_constants(imp, &(mat_ptr->u_heat_source), NO_SPECIES);
    mat_ptr->len_u_heat_source = num_const;
    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_heat_source);
  } else {
    if (model_read == -1) {
      GOMA_EH(model_read, "Heat Source model invalid");
    }
    GOMA_EH(model_read, "Heat Source");
  }

  ECHO(es, echo_file);

  model_read = look_for_mat_prop(imp, "Second Level Set Heat Source", &(i0), v0, NO_USER, NULL,
                                 model_name, SCALAR_INPUT, &NO_SPECIES, es);

  if (model_read != -1) {
    if (ls == NULL)
      GOMA_EH(GOMA_ERROR,
              "Second Level Set Heat Source requires activation of Level Set Tracking.\n");

    mat_ptr->mp2nd->HeatSourceModel = i0;
    mat_ptr->mp2nd->heatsource = *v0;

    stringup(model_name);

    if (!strcmp(model_name, "CONSTANT")) {
      if (fscanf(imp, "%s", input) != 1) {
        GOMA_EH(GOMA_ERROR, "Expecting trailing keyword for Second Level Set Heat Source.\n");
      }

      stringup(input);

      if (strncmp(input, "POSITIVE", 3) == 0) {
        mat_ptr->mp2nd->heatsourcemask[0] = 0;
        mat_ptr->mp2nd->heatsourcemask[1] = 1;
      } else if (strncmp(input, "NEGATIVE", 3) == 0) {
        mat_ptr->mp2nd->heatsourcemask[0] = 1;
        mat_ptr->mp2nd->heatsourcemask[1] = 0;
      } else {
        GOMA_EH(GOMA_ERROR,
                "Keyword must be POSITIVE or NEGATIVE for Second Level Set Heat Source.\n");
      }
      SPF(endofstring(es), " %s", input);
      if (pfd != NULL) {
        for (i = 0; i < pfd->num_phase_funcs; i++) {
          if (fscanf(imp, "%lf", &(mat_ptr->mp2nd->heatsource_phase[i])) != 1) {
            GOMA_EH(GOMA_ERROR, "error reading phase heat source");
          }
          SPF(endofstring(es), " %g", mat_ptr->mp2nd->heatsource_phase[i]);
        }
      }

    } else {
      GOMA_EH(GOMA_ERROR, "Second Level Set Heat Source model can only be CONSTANT.\n");
    }
  }

  ECHO(es, echo_file);

  /* Initialize for good behavior */
  efv->ev_etch_area = -1;
  efv->ev_etch_depth = -1;
  for (i = 0; i < mat_ptr->Num_Species; i++) {
    /*
     *  set species number equal to max number of species
     *  it is changed to species number of input property
     *  by look_for_mat_prop
     */
    species_no = pd_glob[mn]->Num_Species;
    model_read =
        look_for_mat_prop(imp, "Species Source", mat_ptr->SpeciesSourceModel,
                          mat_ptr->species_source, mat_ptr->u_species_source,
                          mat_ptr->len_u_species_source, model_name, SCALAR_INPUT, &species_no, es);

    fallback_chemkin_generic_prop(&model_read, j, &(mat_ptr->SpeciesSourceModel[i]), FALSE,
                                  mat_ptr);
    if (!strcmp(model_name, "EPOXY")) {
      SpeciesSourceModel = EPOXY;
      model_read = 1;
      mat_ptr->SpeciesSourceModel[species_no] = SpeciesSourceModel;
      if (fscanf(imp, "%lf %lf %lf %lf %lf %lf", &a0, &a1, &a2, &a3, &a4, &a5) != 6) {
        sprintf(err_msg, "Matl %s needs 6 constants for %s %s model.\n", pd_glob[mn]->MaterialName,
                "Species Source", "EPOXY");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      mat_ptr->u_species_source[species_no] = (dbl *)array_alloc(1, 6, sizeof(dbl));

      mat_ptr->len_u_species_source[species_no] = 6;

      mat_ptr->u_species_source[species_no][0] = a0; /* prefactor for k1 */
      mat_ptr->u_species_source[species_no][1] = a1; /* exponent for k1 */
      mat_ptr->u_species_source[species_no][2] = a2; /* prefactor for k2 */
      mat_ptr->u_species_source[species_no][3] = a3; /* exponent for k2 */
      mat_ptr->u_species_source[species_no][4] = a4; /* m power law coefficient */
      mat_ptr->u_species_source[species_no][5] = a5; /* n power law coefficient */

      SPF_DBL_VEC(endofstring(es), 6, mat_ptr->u_species_source[species_no]);
    }

    else if (!strcmp(model_name, "EPOXY_DEA")) {
      SpeciesSourceModel = EPOXY_DEA;
      model_read = 1;
      mat_ptr->SpeciesSourceModel[species_no] = SpeciesSourceModel;
      if (fscanf(imp, "%lf %lf %lf %lf %lf", &a0, &a1, &a2, &a3, &a4) != 5) {
        sr = sprintf(err_msg, "Matl %s needs  5 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Species Source", "EPOXY_DEA");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      mat_ptr->u_species_source[species_no] = (dbl *)array_alloc(1, 5, sizeof(dbl));

      mat_ptr->len_u_species_source[species_no] = 5;

      mat_ptr->u_species_source[species_no][0] = a0; /* prefactor for k1 */
      mat_ptr->u_species_source[species_no][1] = a1; /* exponent for k1 */
      mat_ptr->u_species_source[species_no][2] = a2; /* prefactor for k2, low T */
      mat_ptr->u_species_source[species_no][3] = a3; /* exponent for k2, low T */
      mat_ptr->u_species_source[species_no][4] = a4; /* prefactor for k2, mid T */

      SPF_DBL_VEC(endofstring(es), 5, mat_ptr->u_species_source[species_no]);
    }

    else if (!strcmp(model_name, "SSM_BOND")) {
      SpeciesSourceModel = SSM_BOND;
      model_read = 1;
      mat_ptr->SpeciesSourceModel[species_no] = SpeciesSourceModel;
      if (fscanf(imp, "%lf %lf %lf %lf %lf", &a0, &a1, &a2, &a3, &a4) != 5) {
        sr = sprintf(err_msg, "Matl %s needs  5 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Species Source", "SSM_BOND");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      mat_ptr->u_species_source[species_no] = (dbl *)array_alloc(1, 5, sizeof(dbl));

      mat_ptr->len_u_species_source[species_no] = 5;

      mat_ptr->u_species_source[species_no][0] = a0; /* rate for breakup k2 */
      mat_ptr->u_species_source[species_no][1] = a1; /* rate for aggregation k1 */
      mat_ptr->u_species_source[species_no][2] = a2; /* n0 for breakup eqn */
      mat_ptr->u_species_source[species_no][3] = a3; /* exponent for breakup eqn */
      mat_ptr->u_species_source[species_no][4] = a4; /* exponent for aggregation eqn */

      SPF_DBL_VEC(endofstring(es), 5, mat_ptr->u_species_source[species_no]);
    } else if (!strcmp(model_name, "FOAM_EPOXY")) {
      SpeciesSourceModel = FOAM_EPOXY;
      model_read = 1;
      mat_ptr->SpeciesSourceModel[species_no] = SpeciesSourceModel;
      if (fscanf(imp, "%lf %lf %lf  %lf %lf", &a0, &a1, &a2, &a3, &a4) != 5) {
        sr = sprintf(err_msg, "Matl %s needs  5 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Species Source", "FOAM_EPOXY");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      // Set this species type to be an extrinsic variable
      mat_ptr->ExtrinsicIndependentSpeciesVar[species_no] = 1;
      mat_ptr->u_species_source[species_no] = (dbl *)array_alloc(1, 5, sizeof(dbl));

      mat_ptr->len_u_species_source[species_no] = 5;

      mat_ptr->u_species_source[species_no][0] = a0; /* vapor pressure coef aT, p_vap = aT - bT/T */
      mat_ptr->u_species_source[species_no][1] = a1; /* vapor pressure coef bT, p_vap = aT - bT/T */
      mat_ptr->u_species_source[species_no][2] = a2; /* characteristic velocity */
      mat_ptr->u_species_source[species_no][3] = a3; /* coefficient for condensation */
      mat_ptr->u_species_source[species_no][4] = a4; /* coefficient for evaporation  */

      SPF_DBL_VEC(endofstring(es), 5, mat_ptr->u_species_source[species_no]);
    }

    else if (!strcmp(model_name, "FOAM")) {
      if (MAX_CONC <= 8)
        GOMA_EH(GOMA_ERROR, "MAX_CONC must be greater than 8 for FOAM model");
      SpeciesSourceModel = FOAM;
      model_read = 1;
      mat_ptr->SpeciesSourceModel[species_no] = SpeciesSourceModel;
      if (fscanf(imp, "%lf %lf %lf %lf %lf %lf %lf %lf", &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7) !=
          8) {
        sr = sprintf(err_msg, "Matl %s needs 8 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Species Source", "FOAM");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      mat_ptr->u_species_source[species_no] = (dbl *)array_alloc(1, 8, sizeof(dbl));

      mat_ptr->len_u_species_source[species_no] = 8;

      mat_ptr->u_species_source[species_no][0] = a0; /* prefactor for k1 */
      mat_ptr->u_species_source[species_no][1] = a1; /* exponent for k1 */
      mat_ptr->u_species_source[species_no][2] = a2; /* sigma for k1 */
      mat_ptr->u_species_source[species_no][3] = a3; /* prefactor for k2 */
      mat_ptr->u_species_source[species_no][4] = a4; /* exponent for k2 */
      mat_ptr->u_species_source[species_no][5] = a5; /* sigma for k2 */
      mat_ptr->u_species_source[species_no][6] = a6; /* low temperature limit */
      mat_ptr->u_species_source[species_no][7] = a7; /* high temperature limit */

      SPF_DBL_VEC(endofstring(es), 8, mat_ptr->u_species_source[species_no]);
    } else if (!strcmp(model_name, "FOAM_PMDI_10_RXN")) {
      SpeciesSourceModel = FOAM_PMDI_10_RXN;
      model_read = 1;
      mat_ptr->SpeciesSourceModel[species_no] = SpeciesSourceModel;
      if (fscanf(imp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &a0, &a1, &a2, &a3, &a4,
                 &a5, &a6, &a7, &a8, &a9, &a10, &a11) != 12) {
        sr = sprintf(err_msg, "Matl %s needs 12 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Species Source", "FOAM_PMDI_10_RXN");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      mat_ptr->u_species_source[species_no] = (dbl *)array_alloc(1, 12, sizeof(dbl));

      mat_ptr->len_u_species_source[species_no] = 12;

      mat_ptr->u_species_source[species_no][0] = a0;   /* k0 */
      mat_ptr->u_species_source[species_no][1] = a1;   /* w */
      mat_ptr->u_species_source[species_no][2] = a2;   /* Beta */
      mat_ptr->u_species_source[species_no][3] = a3;   /* C_1 */
      mat_ptr->u_species_source[species_no][4] = a4;   /* C_2 */
      mat_ptr->u_species_source[species_no][5] = a5;   /* m */
      mat_ptr->u_species_source[species_no][6] = a6;   /* n */
      mat_ptr->u_species_source[species_no][7] = a7;   /* b */
      mat_ptr->u_species_source[species_no][8] = a8;   /* T_g0 */
      mat_ptr->u_species_source[species_no][9] = a9;   /* T_ginf */
      mat_ptr->u_species_source[species_no][10] = a10; /* A */
      mat_ptr->u_species_source[species_no][11] = a11; /* En/R */
      SPF_DBL_VEC(endofstring(es), 12, mat_ptr->u_species_source[species_no]);
    } else if (!strcmp(model_name, "FOAM_PMDI_10_RXN_DIVV")) {
      SpeciesSourceModel = FOAM_PMDI_10_RXN;
      model_read = 1;
      mat_ptr->SpeciesSourceModel[species_no] = SpeciesSourceModel;
      if (fscanf(imp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &a0, &a1, &a2, &a3, &a4,
                 &a5, &a6, &a7, &a8, &a9, &a10, &a11) != 12) {
        sr = sprintf(err_msg, "Matl %s needs 12 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Species Source", "FOAM_PMDI_10_RXN");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      mat_ptr->u_species_source[species_no] = (dbl *)array_alloc(1, 12, sizeof(dbl));

      mat_ptr->len_u_species_source[species_no] = 12;

      mat_ptr->u_species_source[species_no][0] = a0;   /* k0 */
      mat_ptr->u_species_source[species_no][1] = a1;   /* w */
      mat_ptr->u_species_source[species_no][2] = a2;   /* Beta */
      mat_ptr->u_species_source[species_no][3] = a3;   /* C_1 */
      mat_ptr->u_species_source[species_no][4] = a4;   /* C_2 */
      mat_ptr->u_species_source[species_no][5] = a5;   /* m */
      mat_ptr->u_species_source[species_no][6] = a6;   /* n */
      mat_ptr->u_species_source[species_no][7] = a7;   /* b */
      mat_ptr->u_species_source[species_no][8] = a8;   /* T_g0 */
      mat_ptr->u_species_source[species_no][9] = a9;   /* T_ginf */
      mat_ptr->u_species_source[species_no][10] = a10; /* A */
      mat_ptr->u_species_source[species_no][11] = a11; /* En/R */
      mat_ptr->ExtrinsicIndependentSpeciesVar[species_no] = 1;

      SPF_DBL_VEC(endofstring(es), 12, mat_ptr->u_species_source[species_no]);
    } else if (!strcmp(model_name, "FOAM_PMDI_10_H2O")) {
      SpeciesSourceModel = FOAM_PMDI_10_H2O;
      model_read = 1;
      mat_ptr->SpeciesSourceModel[species_no] = SpeciesSourceModel;
      mat_ptr->ExtrinsicIndependentSpeciesVar[species_no] = 0;
      if (fscanf(imp, "%lf %lf %lf %lf", &a0, &a1, &a2, &a3) != 4) {
        sr = sprintf(err_msg, "Matl %s needs 4 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Species Source", "FOAM_PMDI_10_H2O");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      mat_ptr->u_species_source[species_no] = (dbl *)array_alloc(1, 4, sizeof(dbl));

      mat_ptr->len_u_species_source[species_no] = 4;

      mat_ptr->u_species_source[species_no][0] = a0; /* n */
      mat_ptr->u_species_source[species_no][1] = a1; /* t_nuc */
      mat_ptr->u_species_source[species_no][2] = a2; /* A */
      mat_ptr->u_species_source[species_no][3] = a3; /* En/R */
      SPF_DBL_VEC(endofstring(es), 4, mat_ptr->u_species_source[species_no]);
    } else if (!strcmp(model_name, "FOAM_PMDI_10_CO2")) {
      SpeciesSourceModel = FOAM_PMDI_10_CO2;
      model_read = 1;
      mat_ptr->SpeciesSourceModel[species_no] = SpeciesSourceModel;
      mat_ptr->ExtrinsicIndependentSpeciesVar[species_no] = 0;
    } else if (!strcmp(model_name, "FOAM_PMDI_10_CO2_LIQ")) {
      SpeciesSourceModel = FOAM_PMDI_10_CO2_LIQ;
      model_read = 1;
      mat_ptr->SpeciesSourceModel[species_no] = SpeciesSourceModel;
      mat_ptr->ExtrinsicIndependentSpeciesVar[species_no] = 0;
    } else if (!strcmp(model_name, "FOAM_PMDI_10_CO2_GAS")) {
      SpeciesSourceModel = FOAM_PMDI_10_CO2_GAS;
      model_read = 1;
      mat_ptr->SpeciesSourceModel[species_no] = SpeciesSourceModel;
      mat_ptr->ExtrinsicIndependentSpeciesVar[species_no] = 0;
    } else if (!strcmp(model_name, "FOAM_PMDI_10_H2O_DIVV")) {
      SpeciesSourceModel = FOAM_PMDI_10_H2O;
      model_read = 1;
      mat_ptr->SpeciesSourceModel[species_no] = SpeciesSourceModel;
      mat_ptr->ExtrinsicIndependentSpeciesVar[species_no] = 1;
      if (fscanf(imp, "%lf %lf %lf %lf", &a0, &a1, &a2, &a3) != 4) {
        sr = sprintf(err_msg, "Matl %s needs 4 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Species Source", "FOAM_PMDI_10_H2O");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      mat_ptr->u_species_source[species_no] = (dbl *)array_alloc(1, 4, sizeof(dbl));

      mat_ptr->len_u_species_source[species_no] = 4;

      mat_ptr->u_species_source[species_no][0] = a0; /* n */
      mat_ptr->u_species_source[species_no][1] = a1; /* t_nuc */
      mat_ptr->u_species_source[species_no][2] = a2; /* A */
      mat_ptr->u_species_source[species_no][3] = a3; /* En/R */
      SPF_DBL_VEC(endofstring(es), 4, mat_ptr->u_species_source[species_no]);
    } else if (!strcmp(model_name, "FOAM_PMDI_10_CO2_DIVV")) {
      SpeciesSourceModel = FOAM_PMDI_10_CO2;
      model_read = 1;
      mat_ptr->SpeciesSourceModel[species_no] = SpeciesSourceModel;
      mat_ptr->ExtrinsicIndependentSpeciesVar[species_no] = 1;
    } else if (!strcmp(model_name, "FOAM_PMDI_10_CO2_LIQ_DIVV")) {
      SpeciesSourceModel = FOAM_PMDI_10_CO2_LIQ;
      model_read = 1;
      mat_ptr->SpeciesSourceModel[species_no] = SpeciesSourceModel;
      mat_ptr->ExtrinsicIndependentSpeciesVar[species_no] = 1;
    } else if (!strcmp(model_name, "FOAM_PMDI_10_CO2_GAS_DIVV")) {
      SpeciesSourceModel = FOAM_PMDI_10_CO2_GAS;
      model_read = 1;
      mat_ptr->SpeciesSourceModel[species_no] = SpeciesSourceModel;
      mat_ptr->ExtrinsicIndependentSpeciesVar[species_no] = 1;
    } else if (!strcmp(model_name, "BUTLER_VOLMER")) {
      if (MAX_CONC <= 4)
        GOMA_EH(GOMA_ERROR, "MAX_CONC must be greater than 4 for Butler_volumer");
      SpeciesSourceModel = BUTLER_VOLMER;
      model_read = 1;
      mat_ptr->SpeciesSourceModel[species_no] = SpeciesSourceModel;
      if (fscanf(imp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7,
                 &a8) != 9) {
        sr = sprintf(err_msg, "Matl %s needs 9 floats for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Species Source", "BUTLER_VOLMER");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->u_species_source[species_no] = (dbl *)array_alloc(1, 9, sizeof(dbl));
      mat_ptr->len_u_species_source[species_no] = 9;
      mat_ptr->u_species_source[species_no][0] = a0; /* stoichiometric coefficient, s */
      mat_ptr->u_species_source[species_no][1] =
          a1; /* product of interfacial area by exchange current density, ai0 */
      mat_ptr->u_species_source[species_no][2] = a2; /* reaction order, beta */
      mat_ptr->u_species_source[species_no][3] = a3; /* reference concentration, cref */
      mat_ptr->u_species_source[species_no][4] = a4; /* anodic transfer coeficient, aa */
      mat_ptr->u_species_source[species_no][5] = a5; /* cathodic transfer coefficient, ac */
      mat_ptr->u_species_source[species_no][6] = a6; /* temperature, T */
      mat_ptr->u_species_source[species_no][7] = a7; /* open-circuit potential, V */
      mat_ptr->u_species_source[species_no][8] = a8; /* number of electrons, n */

      SPF_DBL_VEC(endofstring(es), 9, mat_ptr->u_species_source[species_no]);
    } else if (!strcmp(model_name, "PHOTO_CURING")) {
      SpeciesSourceModel = PHOTO_CURING;
      model_read = 1;
      mat_ptr->SpeciesSourceModel[species_no] = SpeciesSourceModel;
      if (fscanf(imp, "%lf %lf %lf %lf ", &a0, &a1, &a2, &a3) != 4) {
        sr = sprintf(err_msg, "Matl %s needs 4 floats for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Species Source", "PHOTO_CURING");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      if (fscanf(imp, "%lf", &a4) != 1) {
        a4 = 1.0;
      }
      if (fscanf(imp, "%lf", &a5) != 1) {
        a5 = 0.0;
      }
      mat_ptr->u_species_source[species_no] = (dbl *)array_alloc(1, 6, sizeof(dbl));
      mat_ptr->len_u_species_source[species_no] = 4;
      mat_ptr->u_species_source[species_no][0] = a0; /* model bit O2:Radical*/
      mat_ptr->u_species_source[species_no][1] = a1; /* intensity_coeff */
      mat_ptr->u_species_source[species_no][2] = a2; /* functionality*/
      mat_ptr->u_species_source[species_no][3] = a3; /* Rate Arrhenius */
      mat_ptr->u_species_source[species_no][4] = a4; /* Monomer 2nd Order */
      mat_ptr->u_species_source[species_no][5] = a5; /* All 2nd Order Coeff */

      SPF_DBL_VEC(endofstring(es), 6, mat_ptr->u_species_source[species_no]);
    }

    else if (!strcmp(model_name, "ELECTROOSMOTIC")) {
      SpeciesSourceModel = ELECTROOSMOTIC;
      model_read = 1;
      mat_ptr->SpeciesSourceModel[species_no] = SpeciesSourceModel;
      if (fscanf(imp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &a0, &a1, &a2, &a3, &a4, &a5,
                 &a6, &a7, &a8, &a9, &a10) != 11) {
        sr = sprintf(err_msg, "Matl %s needs 11 floats for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Species Source", "ELECTROOSMOTIC");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->u_species_source[species_no] = (dbl *)array_alloc(1, 11, sizeof(dbl));
      mat_ptr->len_u_species_source[species_no] = 11;
      mat_ptr->u_species_source[species_no][0] = a0; /* index of species involved in rxn */
      mat_ptr->u_species_source[species_no][1] = a1; /* stoichiometric coefficient, s */
      mat_ptr->u_species_source[species_no][2] =
          a2; /* product of interfacial area by exchange current density, ai0 */
      mat_ptr->u_species_source[species_no][3] = a3;   /* reaction order, beta */
      mat_ptr->u_species_source[species_no][4] = a4;   /* reference concentration, cref */
      mat_ptr->u_species_source[species_no][5] = a5;   /* anodic transfer coeficient, aa */
      mat_ptr->u_species_source[species_no][6] = a6;   /* cathodic transfer coefficient, ac */
      mat_ptr->u_species_source[species_no][7] = a7;   /* temperature, T */
      mat_ptr->u_species_source[species_no][8] = a8;   /* open-circuit potential, V */
      mat_ptr->u_species_source[species_no][9] = a9;   /* number of electrons, n */
      mat_ptr->u_species_source[species_no][10] = a10; /* electro-osmotic drag coef., nd */

      SPF_DBL_VEC(endofstring(es), 11, mat_ptr->u_species_source[species_no]);
    } else if (!strcmp(model_name, "FOAM_PBE_WATER")) {
      SpeciesSourceModel = FOAM_PBE_WATER;
      model_read = 1;
      mat_ptr->SpeciesSourceModel[species_no] = SpeciesSourceModel;
      if (fscanf(imp, "%lf %lf %lf  %lf", &a0, &a1, &a2, &a3) != 4) {
        sr = sprintf(err_msg, "Matl %s needs  4 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Species Source", "FOAM_PBE_WATER");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      // Set this species type to be an extrinsic variable
      // mat_ptr->ExtrinsicIndependentSpeciesVar[species_no] = 1;
      mat_ptr->u_species_source[species_no] = (dbl *)array_alloc(1, 4, sizeof(dbl));

      mat_ptr->len_u_species_source[species_no] = 4;

      mat_ptr->u_species_source[species_no][0] = a0; /* C0_W */
      mat_ptr->u_species_source[species_no][1] = a1; /* A_W */
      mat_ptr->u_species_source[species_no][2] = a2; /* E_W */
      mat_ptr->u_species_source[species_no][3] = a3; /* Delta_H_W */

      SPF_DBL_VEC(endofstring(es), 4, mat_ptr->u_species_source[species_no]);
    } else if (!strcmp(model_name, "FOAM_PBE_OH")) {
      SpeciesSourceModel = FOAM_PBE_OH;
      model_read = 1;
      mat_ptr->SpeciesSourceModel[species_no] = SpeciesSourceModel;
      if (fscanf(imp, "%lf %lf %lf %lf %lf %lf", &a0, &a1, &a2, &a3, &a4, &a5) != 6) {
        sr = sprintf(err_msg, "Matl %s needs  6 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Species Source", "FOAM_PBE_OH");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      // Set this species type to be an extrinsic variable
      // mat_ptr->ExtrinsicIndependentSpeciesVar[species_no] = 1;
      mat_ptr->u_species_source[species_no] = (dbl *)array_alloc(1, 6, sizeof(dbl));

      mat_ptr->len_u_species_source[species_no] = 6;

      mat_ptr->u_species_source[species_no][0] = a0; /* C0_OH */
      mat_ptr->u_species_source[species_no][1] = a1; /* A_OH */
      mat_ptr->u_species_source[species_no][2] = a2; /* E_OH */
      mat_ptr->u_species_source[species_no][3] = a3; /* Delta_H_OH */
      mat_ptr->u_species_source[species_no][4] = a4; /* C0_NCO */
      mat_ptr->u_species_source[species_no][5] = a5; /* Gelling Point */
      SPF_DBL_VEC(endofstring(es), 5, mat_ptr->u_species_source[species_no]);
    } else if (!strcmp(model_name, "FOAM_PBE_CO2_L")) {
      SpeciesSourceModel = FOAM_PBE_CO2_L;
      model_read = 1;
      mat_ptr->SpeciesSourceModel[species_no] = SpeciesSourceModel;
      if (fscanf(imp, "%lf", &a0) != 1) {
        sr = sprintf(err_msg, "Matl %s needs  1 constant for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Species Source", "FOAM_PBE_CO2_L");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      // Set this species type to be an extrinsic variable
      // mat_ptr->ExtrinsicIndependentSpeciesVar[species_no] = 1;
      mat_ptr->u_species_source[species_no] = (dbl *)array_alloc(1, 1, sizeof(dbl));

      mat_ptr->len_u_species_source[species_no] = 1;

      mat_ptr->u_species_source[species_no][0] = a0; /* M_CO2 */

      SPF_DBL_VEC(endofstring(es), 1, mat_ptr->u_species_source[species_no]);
    } else if (!strcmp(model_name, "FOAM_PBE_CO2_G")) {
      model_read = 1;
      mat_ptr->SpeciesSourceModel[species_no] = FOAM_PBE_CO2_G;
    } else if (!strcmp(model_name, "FOAM_PBE_BA_L")) {
      SpeciesSourceModel = FOAM_PBE_BA_L;
      model_read = 1;
      mat_ptr->SpeciesSourceModel[species_no] = SpeciesSourceModel;
      if (fscanf(imp, "%lf %lf %lf %lf %lf", &a0, &a1, &a2, &a3, &a4) != 5) {
        sr = sprintf(err_msg, "Matl %s needs  5 constant for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Species Source", "FOAM_PBE_BA_L");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      // Set this species type to be an extrinsic variable
      // mat_ptr->ExtrinsicIndependentSpeciesVar[species_no] = 1;
      mat_ptr->u_species_source[species_no] = (dbl *)array_alloc(1, 5, sizeof(dbl));

      mat_ptr->len_u_species_source[species_no] = 5;

      mat_ptr->u_species_source[species_no][0] = a0; /* M_BA */
      mat_ptr->u_species_source[species_no][1] = a1; /* Lambda, latent heat */
      mat_ptr->u_species_source[species_no][2] = a2; /* G0 growth rate */
      mat_ptr->u_species_source[species_no][3] = a3; /* T0 reference temp */
      mat_ptr->u_species_source[species_no][4] = a4; /* M_NCO */
      SPF_DBL_VEC(endofstring(es), 1, mat_ptr->u_species_source[species_no]);
    } else if (!strcmp(model_name, "FOAM_PBE_BA_G")) {
      model_read = 1;
      mat_ptr->SpeciesSourceModel[species_no] = FOAM_PBE_BA_G;
    }

    else if (!strcmp(model_name, "ELECTRODE_KINETICS")) {
      model_read = 1;
      mat_ptr->SpeciesSourceModel[species_no] = ELECTRODE_KINETICS;
    }

    else if (!strcmp(model_name, "ION_REACTIONS")) /*  RSL 3/19/01  */
    {
      model_read = 1;
      mat_ptr->SpeciesSourceModel[species_no] = ION_REACTIONS;
    }

    else if (!strcmp(model_name, "DROP_EVAP")) {
      SpeciesSourceModel = DROP_EVAP;
      model_read = 1;
      if (read_bc_mp == -1)
        read_bc_mp = DROP_EVAP;
      mat_ptr->SpeciesSourceModel[species_no] = SpeciesSourceModel;
      if (fscanf(imp, "%lf %lf %lf %lf %lf", &a0, &a1, &a2, &a3, &a4) != 5) {
        sr = sprintf(err_msg, "Matl %s needs 5 floats for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Species Source", "DROP_EVAP");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->u_species_source[species_no] = (dbl *)array_alloc(1, 5, sizeof(dbl));
      mat_ptr->len_u_species_source[species_no] = 5;
      mat_ptr->u_species_source[species_no][0] = a0; /* liquid droplet concentration*/
      mat_ptr->u_species_source[species_no][1] = a1; /* droplet radius*/
      mat_ptr->u_species_source[species_no][2] = a2; /* droplet number concentration */
      mat_ptr->u_species_source[species_no][3] = a3; /* fraction solvent evaporation */
      mat_ptr->u_species_source[species_no][4] = a4; /* Ideal Gas constant */

      SPF_DBL_VEC(endofstring(es), 5, mat_ptr->u_species_source[species_no]);
    } else if (!strcmp(model_name, "ETCHING_KOH")) {
      model_read = 1;
      mat_ptr->SpeciesSourceModel[species_no] = ETCHING_KOH;
    }

    else if (!strcmp(model_name, "ETCHING_KOH_EXTERNAL")) {
      if (fscanf(imp, "%s", input) != 1) {
        GOMA_EH(GOMA_ERROR, "Expecting trailing keyword for ETCHING_KOH_EXTERNAL model.\n");
      }

      ii = 0;
      for (j = 0; j < efv->Num_external_field; j++) {
        if (!strcmp(efv->name[j], input)) {
          mat_ptr->species_source_external_field_index = j;
          efv->ev_etch_area = j;
        }
      }

      if (fscanf(imp, "%s", input) != 1) {
        GOMA_EH(GOMA_ERROR, "Expecting trailing keyword for ETCHING_KOH_EXTERNAL model.\n");
      }

      for (j = 0; j < efv->Num_external_field; j++) {
        if (!strcmp(efv->name[j], input)) {
          ii = 1;
          efv->ev_etch_depth = j;
        }
      }

      if (ii == 0) {
        GOMA_EH(GOMA_ERROR,
                "Must activate two external fields to use this ETCHING_KOH_EXTERNAL model");
      }

      model_read = 1;
      mat_ptr->SpeciesSourceModel[species_no] = ETCHING_KOH_EXT;

    } else if (model_read == -1) {
      GOMA_EH(model_read, "Species Source model invalid. May need more cards for other species");
    }
    ECHO(es, echo_file);

    if (ls != NULL) {
      species_no = pd_glob[mn]->Num_Species;
      model_read = look_for_mat_prop(
          imp, "Second Level Set Species Source", mat_ptr->mp2nd->SpeciesSourceModel,
          mat_ptr->mp2nd->speciessource, NO_USER, NULL, model_name, SCALAR_INPUT, &species_no, es);

      if (model_read != -1) {
        if (ls == NULL)
          GOMA_EH(GOMA_ERROR,
                  "Second Level Set Species Source requires activation of Level Set Tracking.\n");

        // mat_ptr->mp2nd->SpeciesSourceModel[species_no] = i0;
        // mat_ptr->mp2nd->speciessource[species_no] = v0[species_no];

        stringup(model_name);

        if (!strcmp(model_name, "CONSTANT")) {
          if (fscanf(imp, "%s", input) != 1) {
            GOMA_EH(GOMA_ERROR,
                    "Expecting trailing keyword for Second Level Set Species Source.\n");
          }

          stringup(input);

          if (strncmp(input, "POSITIVE", 3) == 0) {
            mat_ptr->mp2nd->speciessourcemask[0][species_no] = 0;
            mat_ptr->mp2nd->speciessourcemask[1][species_no] = 1;
          } else if (strncmp(input, "NEGATIVE", 3) == 0) {
            mat_ptr->mp2nd->speciessourcemask[0][species_no] = 1;
            mat_ptr->mp2nd->speciessourcemask[1][species_no] = 0;
          } else {
            GOMA_EH(GOMA_ERROR,
                    "Keyword must be POSITIVE or NEGATIVE for Second Level Set Heat Source.\n");
          }
          SPF(endofstring(es), " %s", input);
          if (pfd != NULL) {
            for (i = 0; i < pfd->num_phase_funcs; i++) {
              if (fscanf(imp, "%lf", &(mat_ptr->mp2nd->speciessource_phase[i][species_no])) != 1) {
                GOMA_EH(GOMA_ERROR, "error reading phase species source");
              }
              SPF(endofstring(es), " %g", mat_ptr->mp2nd->speciessource_phase[i][species_no]);
            }
          }
        } else {
          GOMA_EH(GOMA_ERROR, "Second Level Set Species Source model can only be CONSTANT.\n");
        }
      }

      ECHO(es, echo_file);

      species_no = pd_glob[mn]->Num_Species;
      model_read = look_for_mat_prop(imp, "Level Set Species Width",
                                     mat_ptr->mp2nd->use_species_source_width,
                                     mat_ptr->mp2nd->species_source_width, NO_USER, NULL,
                                     model_name, SCALAR_INPUT, &species_no, es);

      if (model_read != -1) {
        if (ls == NULL)
          GOMA_EH(GOMA_ERROR,
                  "Level Set Species Width requires activation of Level Set Tracking.\n");

        mat_ptr->mp2nd->use_species_source_width[species_no] = 1;
        // mat_ptr->mp2nd->species_source_width[species_no] = v0[species_no];

        stringup(model_name);

        if (strcmp(model_name, "CONSTANT")) {
          GOMA_EH(GOMA_ERROR, "Level Set Species Width can only be CONSTANT.\n");
        }
      }

      ECHO(es, echo_file);
    }
  }

  model_read = look_for_mat_prop(imp, "Current Source", &(mat_ptr->CurrentSourceModel),
                                 &(mat_ptr->current_source), NO_USER, NULL, model_name,
                                 SCALAR_INPUT, &NO_SPECIES, es);

  if (!strcmp(model_name, "ELECTRODE_KINETICS")) /* added by KSC: 10/21/98 */
  {
    model_read = 1;
    mat_ptr->CurrentSourceModel = ELECTRODE_KINETICS;
  } else if (!strcmp(model_name, "BUTLER_VOLMER")) /* added by KSC: 04-21-2006 */

  {
    model_read = 1;
    mat_ptr->CurrentSourceModel = BUTLER_VOLMER;
    num_const = read_constants(imp, &(mat_ptr->u_current_source), NO_SPECIES);
    mat_ptr->len_u_current_source = num_const;
    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_current_source);
  } else if (!strcmp(model_name, "FICKIAN_CHARGED")) /* added by KSC: 10/5/00 */
  {
    model_read = 1;
    mat_ptr->CurrentSourceModel = FICKIAN_CHARGED;
  } else if (!strcmp(model_name, "STEFAN_MAXWELL_CHARGED")) /* added by KSC: 10/5/00 */
  {
    model_read = 1;
    mat_ptr->CurrentSourceModel = STEFAN_MAXWELL_CHARGED;
  } else if (!strcmp(model_name, "NET_CHARGE")) /* added by KSC: 5/11/02 */
  {
    model_read = 1;
    mat_ptr->CurrentSourceModel = NET_CHARGE;
  } else if (!strcmp(model_name, "DEBYE_HUCKEL")) /* added by ACS: 7/11/03 */
  {
    model_read = 1;
    mat_ptr->CurrentSourceModel = DEBYE_HUCKEL;
    num_const = read_constants(imp, &(mat_ptr->u_current_source), NO_SPECIES);
    mat_ptr->len_u_heat_source = num_const;
    SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_current_source);
  } else if (model_read == -1) {
    GOMA_EH(model_read, "Current Source model invalid");
  }

  ECHO(es, echo_file);

  model_read =
      look_for_mat_prop(imp, "Moment Source", &(mat_ptr->MomentSourceModel),
                        &(mat_ptr->moment_source), &(mat_ptr->u_moment_source),
                        &(mat_ptr->len_u_moment_source), model_name, SCALAR_INPUT, &NO_SPECIES, es);

  if (model_read == -1) {
    if (!strcmp(model_name, "FOAM_PMDI_10")) {
      mat_ptr->MomentSourceModel = FOAM_PMDI_10;
      model_read = 1;
      num_const = read_constants(imp, &(mat_ptr->u_moment_source), NO_SPECIES);

      /* Requires growth rate and coalescence rate constants */
      if (num_const < 2) {
        sr = sprintf(err_msg, "Matl %s needs 2 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Moment Source", "FOAM_PMDI_10");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_moment_source = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_moment_source);
    }
    if (!strcmp(model_name, "CONSTANT_GROWTH")) {
      mat_ptr->MomentSourceModel = MOMENT_CONSTANT_GROWTH;
      model_read = 1;
      num_const = read_constants(imp, &(mat_ptr->u_moment_source), NO_SPECIES);

      /* Requires growth rate and coalescence rate constants */
      if (num_const < 3) {
        sr = sprintf(err_msg, "Matl %s needs 3 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Moment Source", "CONSTANT_GROWTH");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_moment_source = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_moment_source);
    } else if (!strcmp(model_name, "FOAM_PBE")) {
      mat_ptr->MomentSourceModel = FOAM_PBE;
      model_read = 1;
      num_const = read_constants(imp, &(mat_ptr->u_moment_source), NO_SPECIES);
      mat_ptr->len_u_moment_source = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_moment_source);
    }
    ECHO(es, echo_file);
  }

  /*
   * Before we go on let's check to see if a Source term needed extra material Properties
   *    i.e., Non-condensable for DROP_EVAP models
   */
  rewind(imp);

  if (read_bc_mp == DROP_EVAP) {
    iread = look_for_optional(imp, "Non-condensable Molecular Weight", input, '=');
    if (iread != -1) {
      if (fscanf(imp, "%s %d %lf", model_name, &ii, &mw) != 3) {
        GOMA_EH(GOMA_ERROR, "Error reading non-condensable MW: e.g. CONSTANT species_no  MW");
      } else {
        mat_ptr->molecular_weight[mat_ptr->Num_Species_Eqn] = mw;
        SPF(es, "%s = %s %d %.4g", "Non-condensable Molecular Weight", model_name, ii, mw);
      }
      ECHO(es, echo_file);
    }

    iread = look_for_optional(imp, "Non-volatile Molar Volume", input, '=');
    if (iread != -1) {
      if (fscanf(imp, "%s %d %lf", model_name, &ii, &mv) != 3) {
        GOMA_EH(GOMA_ERROR,
                "Error reading non-volatile Molar Volume: e.g. CONSTANT  species_id  MV");
      } else {
        mat_ptr->molar_volume[pd_glob[mn]->Num_Species_Eqn] = mv;
        SPF(es, "%s = %s %d %.4g", "Non-volatile Molar Volume", model_name, ii, mw);
      }
      ECHO(es, echo_file);
    }

    iread = look_for_optional(imp, "Non-volatile Specific Volume", input, '=');
    if (iread != -1) {
      if (fscanf(imp, "%s %d %lf", model_name, &ii, &mv) != 3) {
        GOMA_EH(GOMA_ERROR,
                "Error reading non-volatile Specific Volume: e.g. CONSTANT  species_id  MV");
      } else {
        mat_ptr->specific_volume[pd_glob[mn]->Num_Species_Eqn] = mv;
        SPF(es, "%s = %s %d %.4g", "Non-volatile Specific Volume", model_name, ii, mw);
      }
      ECHO(es, echo_file);
    }

    iread = look_for_optional(imp, "Flory-Huggins parameters", input, '=');
    if (iread != -1) {
      n_species = pd_glob[mn]->Num_Species_Eqn + 1;
      /*number of independent interaction parameters */
      n_ij = (n_species * n_species - n_species) / 2;

      if (fscanf(imp, "%s", model_name) != 1) {
        GOMA_EH(GOMA_ERROR, "Error reading F-H parameter model name: e.g. CONSTANT");
      } else {
        for (i = 0; i < n_ij; i++) /* reading the chi parameters */
        {
          if (fscanf(imp, "%d %d %lf", &ii, &jj, &chi_ij) != 3) {
            GOMA_EH(GOMA_ERROR, "Error:must have three entries, i, j, and chi(i,j)");
          }
          mat_ptr->flory_param[ii][jj] = chi_ij;
          mat_ptr->flory_param[jj][ii] =
              chi_ij * mat_ptr->molar_volume[jj] / mat_ptr->molar_volume[ii];
        }
        for (k = 0; k < n_species; k++) {
          mat_ptr->flory_param[k][k] = 0.;
        }
        SPF(es, "%s = %s %d %d %.4g", "Flory-Huggins parameters", model_name, ii, jj, chi_ij);
      }
      ECHO(es, echo_file);
    }
  }
  /* End of material property re-read  */

  ECHO("\n---Initialization\n", echo_file);

  Num_Var_Init_Mat[mn] = 0;
  while ((iread = look_forward_optional(imp, "Initialize", input, '=')) == 1) {
    Var_init_mat[mn][Num_Var_Init_Mat[mn]].len_u_pars = -1;
    /*
     *  Read the variable name to be fixed
     */
    if (fscanf(imp, "%80s", input) != 1) {
      sprintf(err_msg, "Error reading variable for initialization in material, %s",
              mat_ptr->Material_Name);
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    (void)strip(input);
    var = variable_string_to_int(input, "Variable for matrl initialization");
    if (var >= 0) {
      Var_init_mat[mn][Num_Var_Init_Mat[mn]].var = var;
    } else {
      sprintf(err_msg, "Invalid choice of initialization variable in material, %s",
              mat_ptr->Material_Name);
      GOMA_EH(GOMA_ERROR, err_msg);
    }

    if (fscanf(imp, "%d %lf", &Var_init_mat[mn][Num_Var_Init_Mat[mn]].ktype,
               &Var_init_mat[mn][Num_Var_Init_Mat[mn]].init_val) != 2)
      GOMA_EH(GOMA_ERROR, "Error reading initialization data");

    SPF(es, "%s = %s %d %.4g", "Initialize", input, Var_init_mat[mn][Num_Var_Init_Mat[mn]].ktype,
        Var_init_mat[mn][Num_Var_Init_Mat[mn]].init_val);

    if (fscanf(imp, "%d", &Var_init_mat[mn][Num_Var_Init_Mat[mn]].slave_block) != 1) {
      Var_init_mat[mn][Num_Var_Init_Mat[mn]].slave_block = 0;
    } else
      SPF(endofstring(es), " %d", Var_init_mat[mn][Num_Var_Init_Mat[mn]].slave_block);

    Num_Var_Init_Mat[mn]++;
    ECHO(es, echo_file);
  }

  while ((iread = look_forward_optional(imp, "User Initialize", input, '=')) == 1) {
    int curr_var = Num_Var_Init_Mat[mn];
    double tmp;

    Var_init_mat[mn][curr_var].len_u_pars = -1;
    /*
     *  Read the variable name to be fixed
     */
    if (fscanf(imp, "%80s", input) != 1) {
      sprintf(err_msg, "Error reading variable for user initialization in material, %s",
              mat_ptr->Material_Name);
      GOMA_EH(GOMA_ERROR, err_msg);
    }
    (void)strip(input);
    var = variable_string_to_int(input, "Variable for matrl initialization");
    if (var >= 0) {
      Var_init_mat[mn][curr_var].var = var;
    } else {
      sprintf(err_msg, "Invalid choice of user initialization variable in material, %s",
              mat_ptr->Material_Name);
      GOMA_EH(GOMA_ERROR, err_msg);
    }

    if (fscanf(imp, "%d %lf", &Var_init_mat[mn][Num_Var_Init_Mat[mn]].ktype,
               &Var_init_mat[mn][curr_var].init_val) != 2)
      GOMA_EH(GOMA_ERROR, "Error reading initialization data");

    SPF(es, "%s = %s %d %.4g", "User Initialize", input, Var_init_mat[mn][curr_var].ktype,
        Var_init_mat[mn][curr_var].init_val);

    if (fscanf(imp, "%d", &Var_init_mat[mn][curr_var].slave_block) != 1) {
      Var_init_mat[mn][curr_var].slave_block = 0;
    } else
      SPF(endofstring(es), " %d", Var_init_mat[mn][curr_var].slave_block);

    /* add float list */
    Var_init_mat[mn][curr_var].u_pars = alloc_dbl_1(MAX_NUMBER_PARAMS, 0.0);
    Var_init_mat[mn][curr_var].len_u_pars = 0;
    while (fscanf(imp, "%lf ", &tmp) == 1) {
      i = Var_init_mat[mn][curr_var].len_u_pars;
      Var_init_mat[mn][curr_var].u_pars[i] = tmp;
      Var_init_mat[mn][curr_var].len_u_pars++;
      SPF(endofstring(echo_string), " %.4g", tmp);
    }

    Num_Var_Init_Mat[mn]++;
    ECHO(es, echo_file);
  }

  ECHO("\n---Special Inputs\n", echo_file); /* added by PRS 3/17/2009 */

  /************ SHELL PROPERTIES SECTION ********************/

  /*Initialize for good behavior */
  mat_ptr->HeightUFunctionModel = CONSTANT;
  mat_ptr->heightU = 0.0;

  mat_ptr->VeloUFunctionModel = CONSTANT;
  mat_ptr->veloU[0] = 0.0;
  mat_ptr->veloU[1] = 0.0;
  mat_ptr->veloU[2] = 0.0;

  if (pd_glob[mn]->gv[R_LUBP] || pd_glob[mn]->gv[R_LUBP_2] || pd_glob[mn]->gv[R_TFMP_MASS] ||
      pd_glob[mn]->gv[R_TFMP_BOUND] ||
      ((pd_glob[mn]->gv[R_MASS]) && (pd_glob[mn]->MassFluxModel == FICKIAN_SHELL))) {
    model_read = look_for_mat_proptable(
        imp, "Upper Height Function Constants", &(mat_ptr->HeightUFunctionModel),
        &(mat_ptr->heightU), &(mat_ptr->u_heightU_function_constants),
        &(mat_ptr->len_u_heightU_function_constants),
        &(mat_ptr->heightU_function_constants_tableid), model_name, SCALAR_INPUT, &NO_SPECIES, es);

    mat_ptr->heightU_ext_field_index = -1; // Default to NO external field

    if (model_read == -1 && !strcmp(model_name, "CONSTANT_SPEED")) {
      model_read = 1;
      mat_ptr->HeightUFunctionModel = CONSTANT_SPEED;
      num_const = read_constants(imp, &(mat_ptr->u_heightU_function_constants), NO_SPECIES);
      if (num_const < 2) {
        sr = sprintf(err_msg, "Matl %s needs at least 2 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Upper Height Function", "CONSTANT_SPEED");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      if (num_const > 2) {
        /* We may have an external field "height" we will be adding to this model.  Check
         * for it now and flag its existence through the material properties structure
         */
        mat_ptr->heightU_ext_field_index = -1; // Default to NO external field
        if (efv->ev) {
          for (i = 0; i < efv->Num_external_field; i++) {
            if (!strcmp(efv->name[i], "HEIGHT")) {
              mat_ptr->heightU_ext_field_index = i;
            }
          }
        } else {
          GOMA_EH(GOMA_ERROR, " You have a third float on Upper Height Function Constants card, "
                              "but NO external field");
        }
      }

      mat_ptr->len_u_heightU_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_heightU_function_constants);

    } else if (model_read == -1 && !strcmp(model_name, "CONSTANT_SPEED_MELT")) {
      model_read = 1;
      mat_ptr->HeightUFunctionModel = CONSTANT_SPEED_MELT;
      num_const = read_constants(imp, &(mat_ptr->u_heightU_function_constants), NO_SPECIES);
      if (num_const < 2) {
        sr = sprintf(err_msg, "Matl %s needs at least 2 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Upper Height Function", "CONSTANT_SPEED");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_heightU_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_heightU_function_constants);

    } else if (model_read == -1 && !strcmp(model_name, "CONSTANT_SPEED_DEFORM")) {
      model_read = 1;
      mat_ptr->HeightUFunctionModel = CONSTANT_SPEED_DEFORM;
      num_const = read_constants(imp, &(mat_ptr->u_heightU_function_constants), NO_SPECIES);
      if (num_const < 5) {
        sr = sprintf(err_msg, "Matl %s needs at least 5 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Upper Height Function", "CONSTANT_SPEED_DEFORM");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_heightU_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_heightU_function_constants);

    } else if (model_read == -1 && !strcmp(model_name, "ROLL_ON")) {
      model_read = 1;
      mat_ptr->HeightUFunctionModel = ROLL_ON;
      num_const = read_constants(imp, &(mat_ptr->u_heightU_function_constants), NO_SPECIES);
      if (num_const < 5) {
        sr = sprintf(err_msg, "Matl %s needs 5 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Upper Height Function", "ROLL_ON");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      if (num_const > 5) {
        /* We may have an external field "height" we will be adding to this model.  Check
  mat_ptr->veloU	       * for it now and flag its existence through the material properties
  structure
         */
        mat_ptr->heightU_ext_field_index = -1; // Default to NO external field
        if (efv->ev) {
          for (i = 0; i < efv->Num_external_field; i++) {
            if (!strcmp(efv->name[i], "HEIGHT")) {
              mat_ptr->heightU_ext_field_index = i;
            }
          }
        } else {
          GOMA_EH(GOMA_ERROR, " You have a sixth float on Upper Height Function Constants card, "
                              "but NO external field");
        }
      }
      mat_ptr->len_u_heightU_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_heightU_function_constants);
    } else if (model_read == -1 && !strcmp(model_name, "ROLL_ON_MELT")) {
      model_read = 1;
      mat_ptr->HeightUFunctionModel = ROLL_ON_MELT;
      num_const = read_constants(imp, &(mat_ptr->u_heightU_function_constants), NO_SPECIES);
      if (num_const < 5) {
        sr = sprintf(err_msg, "Matl %s needs 5 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Upper Height Function", "ROLL_ON_MELT");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_heightU_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_heightU_function_constants);
    }

    else if (model_read == -1 && !strcmp(model_name, "ROLL")) {
      model_read = 1;
      mat_ptr->HeightUFunctionModel = ROLL;
      num_const = read_constants(imp, &(mat_ptr->u_heightU_function_constants), NO_SPECIES);
      if (num_const < 8) {
        sr = sprintf(err_msg, "Matl %s needs 8 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Upper Height Function", "ROLL");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_heightU_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_heightU_function_constants);
    }

    else if (model_read == -1 && !strcmp(model_name, "CAP_SQUEEZE")) {
      model_read = 1;
      mat_ptr->HeightUFunctionModel = CAP_SQUEEZE;
      num_const = read_constants(imp, &(mat_ptr->u_heightU_function_constants), NO_SPECIES);
      if (num_const < 5) {
        sr = sprintf(err_msg, "Matl %s needs 5 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Upper Height Function", "CAP_SQUEEZE");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_heightU_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_heightU_function_constants);
    }

    else if (model_read == -1 && !strcmp(model_name, "FLAT_GRAD_FLAT")) {
      model_read = 1;
      mat_ptr->HeightUFunctionModel = FLAT_GRAD_FLAT;
      num_const = read_constants(imp, &(mat_ptr->u_heightU_function_constants), NO_SPECIES);
      if (num_const < 5) {
        sr = sprintf(err_msg, "Matl %s needs 5 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Upper Height Function", "FLAT_GRAD_FLAT");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_heightU_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_heightU_function_constants);
    }

    else if (model_read == -1 && !strcmp(model_name, "FLAT_GRAD_FLAT_MELT")) {
      model_read = 1;
      mat_ptr->HeightUFunctionModel = FLAT_GRAD_FLAT_MELT;
      num_const = read_constants(imp, &(mat_ptr->u_heightU_function_constants), NO_SPECIES);
      if (num_const < 5) {
        sr = sprintf(err_msg, "Matl %s needs 5 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Upper Height Function", "FLAT_GRAD_FLAT_MELT");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_heightU_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_heightU_function_constants);
    }

    else if (model_read == -1 && !strcmp(model_name, "CIRCLE_MELT")) {
      model_read = 1;
      mat_ptr->HeightUFunctionModel = CIRCLE_MELT;
      num_const = read_constants(imp, &(mat_ptr->u_heightU_function_constants), NO_SPECIES);
      if (num_const < 3) {
        sr = sprintf(err_msg, "Matl %s needs 3 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Upper Height Function", "CIRCLE_MELT");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_heightU_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_heightU_function_constants);
    }

    else if (model_read == -1 && !strcmp(model_name, "POLY_TIME")) {
      model_read = 1;
      mat_ptr->HeightUFunctionModel = POLY_TIME;
      num_const = read_constants(imp, &(mat_ptr->u_heightU_function_constants), NO_SPECIES);
      if (num_const < 1) {
        sr = sprintf(err_msg, "Matl %s needs at least 1 constant for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Upper Height Function", "POLY_TIME");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_heightU_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_heightU_function_constants);
    }

    else if (model_read == -1 && !strcmp(model_name, "JOURNAL")) {
      model_read = 1;
      mat_ptr->HeightUFunctionModel = JOURNAL;
      num_const = read_constants(imp, &(mat_ptr->u_heightU_function_constants), NO_SPECIES);
      if (num_const < 2) {
        sr = sprintf(err_msg, "Matl %s needs at least 2 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Upper Height Function", "JOURNAL");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_heightU_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_heightU_function_constants);
    } else if (model_read == -1 && !strcmp(model_name, "EXTERNAL_FIELD")) {
      if (fscanf(imp, "%s", input) != 1) {
        GOMA_EH(GOMA_ERROR,
                "Expecting trailing keyword for Upper height function EXTERNAL_FIELD model.\n");
      }
      model_read = 1;

      ii = 0;
      for (j = 0; j < efv->Num_external_field; j++) {
        if (strcmp(efv->name[j], input) == 0) {
          ii = 1;
          if (mat_ptr->heightU_ext_field_index == -1)
            mat_ptr->heightU_ext_field_index = j;
        }
      }
      if (ii == 0) {
        GOMA_EH(GOMA_ERROR, "Must activate external fields to use this Upper height function "
                            "model.  Field name needed for the EXTERNAL_FIELD");
      }
      mat_ptr->HeightUFunctionModel = EXTERNAL_FIELD;
      /* pick up scale factor for property */
      num_const = read_constants(imp, &(mat_ptr->u_heightU_function_constants), NO_SPECIES);

      mat_ptr->len_u_heightU_function_constants = num_const;
      if (num_const < 1) {
        sr = sprintf(err_msg, "Matl %s expected at least 3 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Upper Height Function", "EXTERNAL_FIELD");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
    }
    /*
     *  TABLE model added to upper height function constant to apply height function model
     *  as computed from videos of experiments of drop merger.
     *
     *  Columns :   Time  Height  dHeight/dTime
     *
     *  Takes scaling factor before filename so units of height can be scaled to problem units
     *
     *  1/23/2017 - AMC
     */
    else if (model_read == -1 && !strcmp(model_name, "TABLE")) {
      if (fscanf(imp, "%s", input) != 1) {
        GOMA_EH(GOMA_ERROR, "Expecting trailing filename for Upper height function TABLE model.\n");
      }
      model_read = 1;

      ii = 0;
      for (j = 0; j < efv->Num_external_field; j++) {
        if (strcmp(efv->name[j], input) == 0) {
          ii = 1;
          if (mat_ptr->heightU_ext_field_index == -1)
            mat_ptr->heightU_ext_field_index = j;
        }
      }
      if (ii == 0) {
        GOMA_EH(GOMA_ERROR, "Must activate external fields to use this Upper height function "
                            "model.  Field name needed for the EXTERNAL_FIELD");
      }
      mat_ptr->HeightUFunctionModel = EXTERNAL_FIELD; // TODO: this is probably wrong
                                                      /* pick up scale factor for property */
      num_const = read_constants(imp, &(mat_ptr->u_heightU_function_constants), NO_SPECIES);

      mat_ptr->len_u_heightU_function_constants = num_const;
      if (num_const < 1) {
        sr = sprintf(err_msg, "Matl %s expected at least 1 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Upper Height Function", "TABLE");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
    } else if (model_read == -1 && !strcmp(model_name, "ROLLER")) {
      model_read = 1;
      mat_ptr->HeightUFunctionModel = ROLLER;
      num_const = read_constants(imp, &(mat_ptr->u_heightU_function_constants), NO_SPECIES);
      if (num_const > 3) {
        /* We may have an external field "height" we will be adding to this model.  Check
         * for it now and flag its existence through the material properties structure
         */
        mat_ptr->heightU_ext_field_index = -1; // Default to NO external field
        if (efv->ev) {
          for (i = 0; i < efv->Num_external_field; i++) {
            if (!strcmp(efv->name[i], "HEIGHT")) {
              mat_ptr->heightU_ext_field_index = i;
            }
          }
        } else {
          GOMA_EH(GOMA_ERROR, " You have a fourth float on Upper Height Function Constants card, "
                              "but NO external field named 'HEIGHT'!");
        }
      }

      if (num_const < 3 || num_const > 4) {
        sr = sprintf(err_msg, "Matl %s needs 3 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Upper Height Function", "ROLLER");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_heightU_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_heightU_function_constants);
    }

    else if (model_read == -1) {
      GOMA_EH(model_read, "Upper Height Function model invalid");
    }
    ECHO(es, echo_file);
  } /* End of LUBP, LUBP_2, TFMP_MASS, and TFMP_BOUND cards */

  if (pd_glob[mn]->gv[R_LUBP] || pd_glob[mn]->gv[R_LUBP_2] ||
      (pd_glob[mn]->gv[R_SHELL_FILMP] && pd_glob[mn]->gv[R_SHELL_FILMH]) ||
      pd_glob[mn]->gv[R_TFMP_MASS] || pd_glob[mn]->gv[R_TFMP_BOUND]) {

    model_read = look_for_mat_proptable(
        imp, "Lower Height Function Constants", &(mat_ptr->HeightLFunctionModel),
        &(mat_ptr->heightL), &(mat_ptr->u_heightL_function_constants),
        &(mat_ptr->len_u_heightL_function_constants),
        &(mat_ptr->heightL_function_constants_tableid), model_name, SCALAR_INPUT, &NO_SPECIES, es);

    mat_ptr->heightL_ext_field_index = -1; // Default to NO external field

    if (model_read == -1 && !strcmp(model_name, "CONSTANT_SPEED")) {
      model_read = 1;
      mat_ptr->HeightLFunctionModel = CONSTANT_SPEED;
      num_const = read_constants(imp, &(mat_ptr->u_heightL_function_constants), NO_SPECIES);
      if (num_const < 2) {
        sr = sprintf(err_msg, "Matl %s needs at least 2 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Lower Height Function", "CONSTANT_SPEED");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      if (num_const > 2) {
        /* We may have an external field "height" we will be adding to this model.  Check
         * for it now and flag its existence through the material properties structure
         */
        mat_ptr->heightL_ext_field_index = -1; // Default to NO external field
        if (efv->ev) {
          for (i = 0; i < efv->Num_external_field; i++) {
            if (!strcmp(efv->name[i], "HEIGHT_L")) {
              mat_ptr->heightL_ext_field_index = i;
            }
          }
        } else {
          GOMA_EH(GOMA_ERROR, " You have a third float on Upper Height Function Constants card, "
                              "but NO external field");
        }
      }
      mat_ptr->len_u_heightL_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_heightL_function_constants);
    }

    else if (model_read == -1 && !strcmp(model_name, "ROLL_ON")) {
      model_read = 1;
      mat_ptr->HeightLFunctionModel = ROLL_ON;
      num_const = read_constants(imp, &(mat_ptr->u_heightL_function_constants), NO_SPECIES);
      if (num_const < 5) {
        sr = sprintf(err_msg, "Matl %s needs 5 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Lower Height Function", "ROLL_ON");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_heightL_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_heightL_function_constants);
    }

    else if (model_read == -1 && !strcmp(model_name, "ROLL")) {
      model_read = 1;
      mat_ptr->HeightLFunctionModel = ROLL;
      num_const = read_constants(imp, &(mat_ptr->u_heightL_function_constants), NO_SPECIES);
      if (num_const < 8) {
        sr = sprintf(err_msg, "Matl %s needs 8 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Lower Height Function", "ROLL");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_heightL_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_heightL_function_constants);
    }

    else if (model_read == -1 && !strcmp(model_name, "EXTERNAL_FIELD")) {
      if (fscanf(imp, "%s", input) != 1) {
        GOMA_EH(GOMA_ERROR,
                "Expecting trailing keyword for Lower height function EXTERNAL_FIELD model.\n");
      }
      model_read = 1;

      ii = 0;
      for (j = 0; j < efv->Num_external_field; j++) {
        if (strcmp(efv->name[j], input) == 0) {
          ii = 1;
          if (mat_ptr->heightL_ext_field_index == -1)
            mat_ptr->heightL_ext_field_index = j;
        }
      }
      if (ii == 0) {
        GOMA_EH(GOMA_ERROR, "Must activate external fields to use this Lower height function "
                            "model.  Field name needed for the EXTERNAL_FIELD");
      }
      mat_ptr->HeightLFunctionModel = EXTERNAL_FIELD;
      /* pick up scale factor for property */
      num_const = read_constants(imp, &(mat_ptr->u_heightL_function_constants), NO_SPECIES);

      mat_ptr->len_u_heightL_function_constants = num_const;
      if (num_const < 3) {
        sr = sprintf(err_msg, "Matl %s expected at least 3 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Lower Height Function", "EXTERNAL_FIELD");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
    }

    else if (model_read == -1) {
      GOMA_EH(model_read, "Lower Height Function model invalid");
    }
    ECHO(es, echo_file);

  } /* End of LUBP, LUBP_2, SHELL_FILMP, SHELL_FILMH, TFMP_MASS, and TFMP_BOUND cards */

  if (pd_glob[mn]->gv[R_LUBP] || pd_glob[mn]->gv[R_LUBP_2] || pd_glob[mn]->gv[R_TFMP_MASS] ||
      pd_glob[mn]->gv[R_TFMP_BOUND]) {
    model_read =
        look_for_mat_prop(imp, "Upper Velocity Function Constants", &(mat_ptr->VeloUFunctionModel),
                          mat_ptr->veloU, NO_USER, NULL, model_name, VECTOR_INPUT, &NO_SPECIES, es);

    if (model_read == -1 && !strcmp(model_name, "ROLL")) {
      model_read = 1;
      mat_ptr->VeloUFunctionModel = ROLL;
      num_const = read_constants(imp, &(mat_ptr->u_veloU_function_constants), NO_SPECIES);
      if (num_const < 1) {
        sr = sprintf(err_msg, "Matl %s needs 1 constant for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Upper Velocity Function", "ROLL");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_veloU_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_veloU_function_constants);
    }

    else if (model_read == -1 && !strcmp(model_name, "TANGENTIAL_ROTATE")) {
      model_read = 1;
      mat_ptr->VeloUFunctionModel = TANGENTIAL_ROTATE;
      num_const = read_constants(imp, &(mat_ptr->u_veloU_function_constants), NO_SPECIES);
      if (num_const < 5) {
        sr = sprintf(err_msg, "Matl %s needs 5 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Upper Velocity Function", "TANGENTIAL_ROTATE");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_veloU_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_veloU_function_constants);
    }

    else if (model_read == -1 && !strcmp(model_name, "LINEAR_TIME")) {
      model_read = 1;
      mat_ptr->VeloUFunctionModel = LINEAR_TIME;
      num_const = read_constants(imp, &(mat_ptr->u_veloU_function_constants), NO_SPECIES);
      if (num_const < 6) {
        sr = sprintf(err_msg, "Matl %s needs 6 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Upper Velocity Function", "LINEAR_TIME");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_veloU_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_veloU_function_constants);
    }

    else if (model_read == -1) {
      GOMA_EH(model_read, "Upper Velocity Function model invalid");
    }
    ECHO(es, echo_file);

  } /* End of LUBP, LUBP_2, TFMP_MASS, and TFMP_BOUND cards */

  if (pd_glob[mn]->gv[R_LUBP] || pd_glob[mn]->gv[R_LUBP_2] ||
      (pd_glob[mn]->gv[R_SHELL_FILMP] && pd_glob[mn]->gv[R_SHELL_FILMH]) ||
      pd_glob[mn]->gv[R_TFMP_MASS] || pd_glob[mn]->gv[R_TFMP_BOUND]) {
    model_read = look_for_mat_prop(
        imp, "Lower Velocity Function Constants", &(mat_ptr->VeloLFunctionModel), mat_ptr->veloL,
        &(mat_ptr->u_veloL_function_constants), &(mat_ptr->len_u_veloL_function_constants),
        model_name, VECTOR_INPUT, &NO_SPECIES, es);

    if (model_read == -1 && !strcmp(model_name, "ROLL")) {
      model_read = 1;
      mat_ptr->VeloLFunctionModel = ROLL;
      num_const = read_constants(imp, &(mat_ptr->u_veloL_function_constants), NO_SPECIES);
      if (num_const < 1) {
        sr = sprintf(err_msg, "Matl %s needs 1 constant for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Lower Velocity Function", "ROLL");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_veloL_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_veloL_function_constants);
    }

    else if (model_read == -1 && !strcmp(model_name, "TANGENTIAL_ROTATE")) {
      model_read = 1;
      mat_ptr->VeloLFunctionModel = TANGENTIAL_ROTATE;
      num_const = read_constants(imp, &(mat_ptr->u_veloL_function_constants), NO_SPECIES);
      if (num_const < 5) {
        sr = sprintf(err_msg, "Matl %s needs 5 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "lower Velocity Function", "TANGENTIAL_ROTATE");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_veloL_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_veloL_function_constants);
    }

    else if (model_read == -1 && !strcmp(model_name, "LINEAR_TIME")) {
      model_read = 1;
      mat_ptr->VeloLFunctionModel = LINEAR_TIME;
      num_const = read_constants(imp, &(mat_ptr->u_veloL_function_constants), NO_SPECIES);
      if (num_const < 6) {
        sr = sprintf(err_msg, "Matl %s needs 6 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "lower Velocity Function", "LINEAR_TIME");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_veloL_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_veloL_function_constants);
    }

    else if (model_read == -1 && !strcmp(model_name, "SLIDER_POLY_TIME")) {
      model_read = 1;
      mat_ptr->VeloLFunctionModel = SLIDER_POLY_TIME;
      num_const = read_constants(imp, &(mat_ptr->u_veloL_function_constants), NO_SPECIES);
      if (num_const < 1) {
        sr = sprintf(err_msg, "Matl %s needs at least 1 constant for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Lower Velocity Function", "SIDER_POLY_TIME");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_veloL_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_veloL_function_constants);
    }

    else if (model_read == -1) {
      GOMA_EH(model_read, "Lower Velocity Function model invalid");
    }
    ECHO(es, echo_file);

  } /* End of LUBP, LUBP_2, TFMP_MASS, SHELL_FILMP, SHELL_FILMH, and TFMP_BOUND cards */

  if (pd_glob[mn]->gv[R_LUBP] || pd_glob[mn]->gv[R_LUBP_2] || pd_glob[mn]->gv[R_TFMP_MASS] ||
      pd_glob[mn]->gv[R_TFMP_BOUND]) {
    model_read = look_for_mat_prop(imp, "Upper Contact Angle", &(mat_ptr->DcaUFunctionModel),
                                   &(mat_ptr->dcaU), NO_USER, NULL, model_name, SCALAR_INPUT,
                                   &NO_SPECIES, es);

    if (model_read == -1 && !strcmp(model_name, "DYNAMIC")) {
      model_read = 1;
      mat_ptr->DcaUFunctionModel = DYNAMIC_CA;
      num_const = read_constants(imp, &(mat_ptr->u_dcaU_function_constants), NO_SPECIES);
      if (num_const < 4) {
        sr = sprintf(err_msg, "Matl %s needs 4 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Upper Contact Angle", "DYNAMIC");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->dcaU = mat_ptr->u_dcaU_function_constants[0];
      mat_ptr->len_u_dcaU_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_dcaU_function_constants);

    } else if (model_read == -1 && !strcmp(model_name, "DYNAMIC_LINEAR")) {
      model_read = 1;
      mat_ptr->DcaUFunctionModel = DYNAMIC_LINEAR_CA;
      num_const = read_constants(imp, &(mat_ptr->u_dcaU_function_constants), NO_SPECIES);
      if (num_const < 4) {
        sr = sprintf(err_msg, "Matl %s needs 4 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Upper Contact Angle", "DYNAMIC_LINEAR");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->dcaU = mat_ptr->u_dcaU_function_constants[0];
      mat_ptr->len_u_dcaU_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_dcaU_function_constants);

    } else if (model_read == -1) {
      GOMA_EH(model_read, "Upper Contact Angle model invalid");
    }
    ECHO(es, echo_file);

    model_read = look_for_mat_prop(imp, "Lower Contact Angle", &(mat_ptr->DcaLFunctionModel),
                                   &(mat_ptr->dcaL), NO_USER, NULL, model_name, SCALAR_INPUT,
                                   &NO_SPECIES, es);

    if (model_read == -1 && !strcmp(model_name, "DYNAMIC")) {
      model_read = 1;
      mat_ptr->DcaLFunctionModel = DYNAMIC_CA;
      num_const = read_constants(imp, &(mat_ptr->u_dcaL_function_constants), NO_SPECIES);
      if (num_const < 4) {
        sr = sprintf(err_msg, "Matl %s needs 4 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Lower Contact Angle", "DYNAMIC");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->dcaL = mat_ptr->u_dcaL_function_constants[0];
      mat_ptr->len_u_dcaL_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_dcaL_function_constants);

    } else if (model_read == -1 && !strcmp(model_name, "DYNAMIC_LINEAR")) {
      model_read = 1;
      mat_ptr->DcaLFunctionModel = DYNAMIC_LINEAR_CA;
      num_const = read_constants(imp, &(mat_ptr->u_dcaL_function_constants), NO_SPECIES);
      if (num_const < 4) {
        sr = sprintf(err_msg, "Matl %s needs 4 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Lower Contact Angle", "DYNAMIC_LINEAR");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->dcaL = mat_ptr->u_dcaL_function_constants[0];
      mat_ptr->len_u_dcaL_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_dcaL_function_constants);

    } else if (model_read == -1) {
      GOMA_EH(model_read, "Lower Contact Angle model invalid");
    }
    ECHO(es, echo_file);

  } /* End of LUBP, LUBP_2, TFMP_MASS, and TFMP_BOUND cards */

  if (pd_glob[mn]->gv[R_LUBP] || pd_glob[mn]->gv[R_LUBP_2] ||
      (pd_glob[mn]->gv[R_SHELL_FILMP] && pd_glob[mn]->gv[R_SHELL_FILMH]) ||
      pd_glob[mn]->gv[R_TFMP_MASS] || pd_glob[mn]->gv[R_TFMP_BOUND]) {
    /* Optional lubrication fluid source term */

    strcpy(search_string, "Lubrication Fluid Source");
    mat_ptr->LubSourceModel = CONSTANT;
    mat_ptr->lubsource = 0.;
    memset(mp_glob[mn]->d_lubsource, 0, sizeof(dbl) * (MAX_VARIABLE_TYPES + MAX_CONC));

    if (look_forward_optional(imp, search_string, input, '=') == 1) {
      if (fscanf(imp, "%s", model_name) != 1) {
        GOMA_EH(GOMA_ERROR, "Need model for lubrication source term ");
      }

      SPF(es, "%s = %s", search_string, model_name);

      if (!strcmp(model_name, "CONSTANT")) {
        mat_ptr->LubSourceModel = CONSTANT;
        mat_ptr->u_lubsource_function_constants = alloc_dbl_1(1, 0.0);
        mat_ptr->len_lubsource = 1;

        if (fscanf(imp, "%lf", &(mat_ptr->lubsource)) != 1) {
          GOMA_EH(GOMA_ERROR, "Lubrication fluid source constant model expects 1 flt");
        }

        SPF(endofstring(es), "%g", mat_ptr->lubsource);

      } else if (!strcmp(model_name, "MELT")) {
        mp_glob[mn]->LubSourceModel = MELT;
        num_const = read_constants(imp, &(mat_ptr->u_lubsource_function_constants), NO_SPECIES);

        if (num_const < 3) {
          sr = sprintf(err_msg, "Matl %s needs 3 constants for %s %s model.\n",
                       pd_glob[mn]->MaterialName, "Lubrication source", "MELT");
          GOMA_EH(GOMA_ERROR, err_msg);
        }

        if (((mat_ptr->veloU[0] + mat_ptr->veloU[1] + mat_ptr->veloU[2]) == 0. &&
             (mat_ptr->veloL[0] + mat_ptr->veloL[1] + mat_ptr->veloL[2]) == 0.) ||
            (mat_ptr->VeloLFunctionModel != CONSTANT)) {
          GOMA_WH(-1, "Lubrication source model MELT must have nonzero slider velocities");
        }
        mat_ptr->len_lubsource = num_const;
        SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_lubsource_function_constants);
      }

      else if (!strcmp(model_name, "CONTINUUM_FLUID")) {
        mp_glob[mn]->LubSourceModel = CONTINUUM_FLUID;
        num_const = read_constants(imp, &(mat_ptr->u_lubsource_function_constants), NO_SPECIES);

        if (num_const < 3) {
          sr = sprintf(err_msg, "Matl %s needs 3 constants for %s %s model.\n",
                       pd_glob[mn]->MaterialName, "Lubrication source", "CONTINUUM_FLUID");
          GOMA_EH(GOMA_ERROR, err_msg);
        }

        mat_ptr->len_lubsource = num_const;
        SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_lubsource_function_constants);
      }

      else {
        GOMA_EH(GOMA_ERROR, "Unrecognized lubrication fluid source model");
      }
    }

  } /* End of LUBP, LUBP_2, SHELL_FILMP, SHELL_FILMH, TFMP_MASS, and TFMP_BOUND cards */

  if (pd_glob[mn]->gv[R_LUBP] || pd_glob[mn]->gv[R_LUBP_2] || pd_glob[mn]->gv[R_TFMP_MASS] ||
      pd_glob[mn]->gv[R_TFMP_BOUND]) {

    /* Optional lubrication momentum source term */

    strcpy(search_string, "Lubrication Momentum Source");
    mat_ptr->lubmomsource[0] = 0.;
    mat_ptr->lubmomsource[1] = 0.;
    mat_ptr->lubmomsource[2] = 0.;
    memset(mp_glob[mn]->d_lubmomsource, 0, sizeof(dbl) * (MAX_VARIABLE_TYPES + MAX_CONC));

    if (look_forward_optional(imp, search_string, input, '=') == 1) {
      if (fscanf(imp, "%s", model_name) != 1) {
        GOMA_EH(GOMA_ERROR, "Need model for lubrication momentum source term ");
      }

      SPF(es, "%s = %s", search_string, model_name);

      if (!strcmp(model_name, "CONSTANT")) {
        mp_glob[mn]->LubMomSourceModel = CONSTANT;

        if (fscanf(imp, "%lf %lf %lf", &(mat_ptr->lubmomsource[0]), &(mat_ptr->lubmomsource[1]),
                   &(mat_ptr->lubmomsource[2])) != 3)

        {
          GOMA_EH(GOMA_ERROR, "Lubrication momentum source constant model expects 3 flts");
        }
        num_const = mp_glob[mn]->len_lubmomsource = 3;

        SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->lubmomsource);

      } else {
        GOMA_EH(GOMA_ERROR, "Unrecognized lubrication momentum source  model");
      }
    }

  } /* End of LUBP, LUBP_2, SHELL_FILMP, SHELL_FILMH, TFMP_MASS, and TFMP_BOUND card */

  /* Shell Energy Cards - heat sources, sinks, etc. */

  if (pd_glob[mn]->gv[R_SHELL_ENERGY] == 1) {
    /* no source terms available.  Feel free to add some!  */
  } /* End of shell_energy cards */

  if (pd_glob[mn]->gv[R_SHELL_FILMP] && pd_glob[mn]->gv[R_SHELL_FILMH]) {

    model_read = look_for_mat_prop(imp, "Film Evaporation Model", &(mat_ptr->FilmEvapModel),
                                   &(mat_ptr->FilmEvap), NO_USER, NULL, model_name, SCALAR_INPUT,
                                   &NO_SPECIES, es);

    if (model_read == -1 && !strcmp(model_name, "CONC_POWER")) {
      model_read = 1;
      mat_ptr->FilmEvapModel = CONC_POWER;
      num_const = read_constants(imp, &(mat_ptr->u_FilmEvap_function_constants), NO_SPECIES);
      if (num_const < 3) {
        sr = sprintf(err_msg, "Matl %s needs 3 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Film Evaporation Model", "CONC_POWER");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_FilmEvap_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_FilmEvap_function_constants);
    }

    else if (model_read == -1) {
      GOMA_EH(model_read, "Film Evaporation Model invalid");
    }

    ECHO(es, echo_file);

    model_read =
        look_for_mat_prop(imp, "Upper Velocity Function Constants", &(mat_ptr->VeloUFunctionModel),
                          mat_ptr->veloU, NO_USER, NULL, model_name, VECTOR_INPUT, &NO_SPECIES, es);

    /* Optional slip term */

    strcpy(search_string, "Slip Coefficient Model");
    mat_ptr->SlipCoeff = 0.;
    memset(mp_glob[mn]->d_SlipCoeff, 0, sizeof(dbl) * (MAX_VARIABLE_TYPES + MAX_CONC));

    if (look_forward_optional(imp, search_string, input, '=') == 1) {
      if (fscanf(imp, "%s", model_name) != 1) {
        GOMA_EH(GOMA_ERROR, "Need model for Slip coefficient ");
      }

      SPF(es, "%s = %s", search_string, model_name);

      if (!strcmp(model_name, "CONSTANT")) {
        mat_ptr->SlipCoeffModel = CONSTANT;
        mat_ptr->u_SlipCoeff_function_constants = alloc_dbl_1(1, 0.0);
        mat_ptr->len_u_SlipCoeff_function_constants = 1;

        if (fscanf(imp, "%lf", &(mat_ptr->SlipCoeff)) != 1)

        {
          GOMA_EH(GOMA_ERROR, "Slip coefficient constant model expects 1 flt");
        }

        SPF(endofstring(es), "%g", mat_ptr->SlipCoeff);
      } else {
        GOMA_EH(GOMA_ERROR, "Slip coefficient model invalid");
      }
    }

    /* Disjoining pressure term */

    model_read = look_for_mat_prop(imp, "Disjoining Pressure Model", &(mat_ptr->DisjPressModel),
                                   &(mat_ptr->DisjPress), NO_USER, NULL, model_name, SCALAR_INPUT,
                                   &NO_SPECIES, es);

    if (model_read == -1 && !strcmp(model_name, "TWO_TERM")) {
      model_read = 1;
      mat_ptr->DisjPressModel = TWO_TERM;
      num_const = read_constants(imp, &(mat_ptr->u_DisjPress_function_constants), NO_SPECIES);
      if (num_const < 5) {
        sr = sprintf(err_msg, "Matl %s needs 5 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Disjoining Pressure Model", "TWO_TERM");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_DisjPress_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_DisjPress_function_constants);
    }

    if (model_read == -1 && !strcmp(model_name, "TWO_TERM_EXT_CA")) {
      model_read = 1;
      mat_ptr->DisjPressModel = TWO_TERM_EXT_CA;
      num_const = read_constants(imp, &(mat_ptr->u_DisjPress_function_constants), NO_SPECIES);
      if (num_const < 4) {
        sr = sprintf(err_msg, "Matl %s needs 4 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Disjoining Pressure Model", "TWO_TERM_EXT_CA");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_DisjPress_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_DisjPress_function_constants);
    }

    if (model_read == -1 && !strcmp(model_name, "ONE_TERM")) {
      model_read = 1;
      mat_ptr->DisjPressModel = ONE_TERM;
      num_const = read_constants(imp, &(mat_ptr->u_DisjPress_function_constants), NO_SPECIES);
      if (num_const < 3) {
        sr = sprintf(err_msg, "Matl %s needs 3 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Disjoining Pressure Model", "ONE_TERM");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_DisjPress_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_DisjPress_function_constants);
    }

    else if (model_read == -1) {
      GOMA_EH(model_read, "Disjoining Pressure model invalid");
    }

    ECHO(es, echo_file);
  } /* End of SHELL_FILMP and SHELL_FILMH cards */

  if (pd_glob[mn]->gv[R_SHELL_PARTC] == 1) {

    model_read = look_for_mat_prop(imp, "Diffusion Coefficient Model", &(mat_ptr->DiffCoeffModel),
                                   &(mat_ptr->DiffCoeff), NO_USER, NULL, model_name, SCALAR_INPUT,
                                   &NO_SPECIES, es);

    if (model_read == -1 && !strcmp(model_name, "STOKES_EINSTEIN")) {
      model_read = 1;
      mat_ptr->DiffCoeffModel = STOKES_EINSTEIN;
      num_const = read_constants(imp, &(mat_ptr->u_DiffCoeff_function_constants), NO_SPECIES);
      if (num_const < 3) {
        sr = sprintf(err_msg, "Matl %s needs 3 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Diffusion Coefficient Model", "STOKES_EINSTEIN");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_DiffCoeff_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->u_DiffCoeff_function_constants);
    }

    else if (model_read == -1) {
      GOMA_EH(model_read, "Diffusion Coefficient Model invalid");
    }

    ECHO(es, echo_file);
  }
  /*
   * Material input for fluid-structural interaction types for
   * lubrication shell elements
   */

  /* Check for mesh equation in any material and any matrix */
  have_mesh_eqn = 0;
  for (i = 0; i < upd->Total_Num_Matrices; i++) {
    if (upd->ep[i][R_MESH1] > -1) {
      have_mesh_eqn = 1;
    }
  }

  // mat_ptr->FSIModel = FSI_SHELL_ONLY;
  /*added by DSB 7/13; some post proc routines for 2D/shell problems
                                        get confused if there are shell elements and no FSI model is
     specified. Need to check that this doesn't break anything.*/
  if (pd_glob[mn]->gv[R_LUBP] || pd_glob[mn]->gv[R_LUBP_2] || pd_glob[mn]->gv[R_SHELL_FILMP] ||
      pd_glob[mn]->gv[R_SHELL_SAT_OPEN] || pd_glob[mn]->gv[R_SHELL_SAT_OPEN_2] ||
      pd_glob[mn]->gv[R_SHELL_SAT_1] ||
      (pd_glob[mn]->gv[R_SHELL_NORMAL1] && pd_glob[mn]->gv[R_SHELL_NORMAL2] &&
       pd_glob[mn]->gv[R_SHELL_NORMAL3]) ||
      pd_glob[mn]->gv[R_TFMP_MASS] || pd_glob[mn]->gv[R_TFMP_BOUND] ||
      ((pd_glob[mn]->gv[R_MASS]) && (pd_glob[mn]->MassFluxModel == FICKIAN_SHELL))) {
    model_read = look_for_mat_prop(imp, "FSI Deformation Model", &(mat_ptr->FSIModel), &(a0),
                                   NO_USER, NULL, model_name, NO_INPUT, &NO_SPECIES, es);

    if (!strcmp(model_name, "FSI_MESH_BOTH")) {
      mat_ptr->FSIModel = FSI_MESH_BOTH;
      GOMA_EH(model_read, "This FSI Deformation Model is not currently implemented!");

    } else if (!strcmp(model_name, "FSI_MESH_CONTINUUM")) {
      // PRS: Made this a GOMA_WH instead of a GOMA_EH as if the mesh equations are in a higher
      // block #, this trips.
      if (have_mesh_eqn == 0)
        GOMA_WH(-1, " Must have mesh continuum equations on somewhere for FSI_MESH_CONTINUUM");
      mat_ptr->FSIModel = FSI_MESH_CONTINUUM;

    } else if (!strcmp(model_name, "FSI_MESH_SHELL")) {
      mat_ptr->FSIModel = FSI_MESH_SHELL;
      GOMA_EH(model_read, "This FSI Deformation Model is not currently implemented!");

    } else if (!strcmp(model_name, "FSI_SHELL_ONLY")) {
      mat_ptr->FSIModel = FSI_SHELL_ONLY;

    } else if (!strcmp(model_name, "FSI_MESH_UNDEF")) {
      mat_ptr->FSIModel = FSI_MESH_UNDEF;

    } else if (!strcmp(model_name, "FSI_MESH_ONEWAY")) {
      mat_ptr->FSIModel = FSI_MESH_ONEWAY;

    } else if (!strcmp(model_name, "FSI_REALSOLID_CONTINUUM")) {
      mat_ptr->FSIModel = FSI_REALSOLID_CONTINUUM;

    } else if (!strcmp(model_name, "FSI_SHELL_ONLY_MESH")) {
      mat_ptr->FSIModel = FSI_SHELL_ONLY_MESH;

    } else if (!strcmp(model_name, "FSI_SHELL_ONLY_UNDEF")) {
      mat_ptr->FSIModel = FSI_SHELL_ONLY_UNDEF;

    } else {
      GOMA_EH(model_read, "This FSI Deformation Model is not valid!");
    }

    ECHO(es, echo_file);
    model_read =
        look_for_mat_prop(imp, "Lubrication Integration Model", &(mat_ptr->LubIntegrationModel),
                          &(a0), NO_USER, NULL, model_name, NO_INPUT, &NO_SPECIES, es);
    if (model_read == -1) {
      if (!strcmp(model_name, "GAUSSIAN")) {
        mat_ptr->LubIntegrationModel = LUB_VISCINT_GAUSSIAN;
        model_read = 1;
      } else if (!strcmp(model_name, "ANALYTICAL")) {
        mat_ptr->LubIntegrationModel = LUB_VISCINT_ANALYTICAL;
        model_read = 1;
      } else if (!strcmp(model_name, "POWER_LAW")) {
        mat_ptr->LubIntegrationModel = LUB_VISCINT_POWERLAW;
        model_read = 1;
      }
    }
    if (model_read == 1) {
      if (fscanf(imp, "%d", &(mat_ptr->LubInt_NGP)) != 1) {
        sr = sprintf(err_msg, "Error reading Lubrication Gauss pts %s", model_name);
        GOMA_EH(GOMA_ERROR, err_msg);
      } else {
        SPF(endofstring(es), " %d", mat_ptr->LubInt_NGP);
      }
      if (mat_ptr->LubIntegrationModel == LUB_VISCINT_GAUSSIAN ||
          mat_ptr->LubIntegrationModel == LUB_VISCINT_ANALYTICAL) {
        /* These are Gaussian integration constants shifted to the interval [0,1]
           Source: Abramowitz & Stegun, Handbook of Mathematical Functions, 9th printing, 1972*/
        if (mat_ptr->LubInt_NGP == 2) {
          mat_ptr->Lub_wts[0] = 0.5;
          mat_ptr->Lub_wts[1] = 0.5;
          mat_ptr->Lub_gpts[0] = 0.5 * (1. - sqrt(1. / 3.));
          mat_ptr->Lub_gpts[1] = 0.5 * (1. + sqrt(1. / 3.));
        } else if (mat_ptr->LubInt_NGP == 3) {
          mat_ptr->Lub_wts[0] = 5. / 18.;
          mat_ptr->Lub_wts[1] = 4. / 9;
          mat_ptr->Lub_wts[2] = 5. / 18.;
          mat_ptr->Lub_gpts[0] = 0.5 * (1. - sqrt(0.6));
          mat_ptr->Lub_gpts[1] = 0.5;
          mat_ptr->Lub_gpts[2] = 0.5 * (1. + sqrt(0.6));
        } else if (mat_ptr->LubInt_NGP == 4) {
          mat_ptr->Lub_wts[0] = 0.5 * 0.347854845137454;
          mat_ptr->Lub_wts[1] = 0.5 * 0.652145154862546;
          mat_ptr->Lub_wts[2] = mat_ptr->Lub_wts[1];
          mat_ptr->Lub_wts[3] = mat_ptr->Lub_wts[0];
          mat_ptr->Lub_gpts[0] = 0.5 * (1. - 0.861136311594053);
          mat_ptr->Lub_gpts[1] = 0.5 * (1. - 0.339981043584856);
          mat_ptr->Lub_gpts[2] = 0.5 * (1. + 0.339981043584856);
          mat_ptr->Lub_gpts[3] = 0.5 * (1. + 0.861136311594053);
        } else if (mat_ptr->LubInt_NGP == 5) {
          mat_ptr->Lub_wts[0] = 0.5 * 0.236926885056189;
          mat_ptr->Lub_wts[1] = 0.5 * 0.478628670499366;
          mat_ptr->Lub_wts[2] = 0.5 * 0.568888888888889;
          mat_ptr->Lub_wts[3] = mat_ptr->Lub_wts[1];
          mat_ptr->Lub_wts[4] = mat_ptr->Lub_wts[0];
          mat_ptr->Lub_gpts[0] = 0.5 * (1. - 0.906179845938664);
          mat_ptr->Lub_gpts[1] = 0.5 * (1. - 0.538469310105683);
          mat_ptr->Lub_gpts[2] = 0.5;
          mat_ptr->Lub_gpts[3] = 0.5 * (1. + 0.538469310105683);
          mat_ptr->Lub_gpts[0] = 0.5 * (1. + 0.906179845938664);
        } else {
          GOMA_EH(GOMA_ERROR, "Those integration points not defined yet!");
        }
      } else if (mat_ptr->LubIntegrationModel == LUB_VISCINT_POWERLAW) {
        double n = 0.5, beta0, beta1, beta2;
        if (fscanf(imp, "%lf", &(mat_ptr->LubInt_PL)) != 1) {
          sr = sprintf(err_msg, "Error reading POWER_LAW index %s", model_name);
          GOMA_EH(GOMA_ERROR, err_msg);
        } else {
          n = mat_ptr->LubInt_PL;
        }
        SPF(endofstring(es), " %g", mat_ptr->LubInt_PL);
        beta0 = 1. / (2. * n + 1);
        beta1 = SQUARE(beta0) / (2. * n + 3.);
        beta2 = beta0 + 2. * SQUARE(2. * n + 3) / (n + 1.) * beta1 -
                2. * (2. * n * n + 7. * n + 7.) / (n + 1.) / SQUARE(2. * n + 1) +
                SQUARE((2. * n + 3.) * (n + 2.) / ((n + 1.) * (2. * n + 1.))) / (2. * n + 5.);
        if (mat_ptr->LubInt_NGP == 2 && DOUBLE_NONZERO(n)) {
          double disc = sqrt((1. + n) / (2 * n + 3));
          double fcn11, fcn12;
          mat_ptr->Lub_gpts[0] = (n + 1. - disc) / (n + 2.);
          mat_ptr->Lub_gpts[1] = (n + 1. + disc) / (n + 2.);
          fcn11 = SQUARE(1. - (2. * n + 2) / (2. * n + 1) * mat_ptr->Lub_gpts[0]);
          fcn12 = SQUARE(1. - (2. * n + 2) / (2. * n + 1) * mat_ptr->Lub_gpts[1]);
          mat_ptr->Lub_wts[0] = (beta1 - beta0 * fcn12) / (fcn11 - fcn12);
          mat_ptr->Lub_wts[1] = beta0 - mat_ptr->Lub_wts[0];
        } else if (mat_ptr->LubInt_NGP == 3) { /* Trignometric Cubic Solution  */
          double p, q, r, a, b, three_th, mag;
          double f11, f12, f13, f21, f22, f23;

          p = -3 * (2 * n + 3) / (2 * n + 6);
          q = -p * (2 * n + 2) / (2 * n + 5);
          r = -q * (2 * n + 1) / (2 * n + 4) / 3.;
          a = (3 * q - SQUARE(p)) / 3.;
          b = (2 * CUBE(p) - 9 * p * q + 27 * r) / 27.;
          mag = 2 * sqrt(-a / 3.);
          three_th = acos(3 * b / a / mag);

          mat_ptr->Lub_gpts[0] = mag * cos(three_th / 3.) - p / 3.;
          mat_ptr->Lub_gpts[1] = mag * cos(three_th / 3. + 2. / 3. * M_PIE) - p / 3.;
          mat_ptr->Lub_gpts[2] = mag * cos(three_th / 3. + 4. / 3. * M_PIE) - p / 3.;
          f11 = SQUARE(1. - (2. * n + 2) / (2. * n + 1) * mat_ptr->Lub_gpts[0]);
          f12 = SQUARE(1. - (2. * n + 2) / (2. * n + 1) * mat_ptr->Lub_gpts[1]);
          f13 = SQUARE(1. - (2. * n + 2) / (2. * n + 1) * mat_ptr->Lub_gpts[2]);
          f21 =
              SQUARE(1. - 2 * (2 * n + 3) / (2 * n + 1) * mat_ptr->Lub_gpts[0] +
                     (2 * n + 3) / (2 * n + 1) * (n + 2) / (n + 1) * SQUARE(mat_ptr->Lub_gpts[0]));
          f22 =
              SQUARE(1. - 2 * (2 * n + 3) / (2 * n + 1) * mat_ptr->Lub_gpts[1] +
                     (2 * n + 3) / (2 * n + 1) * (n + 2) / (n + 1) * SQUARE(mat_ptr->Lub_gpts[1]));
          f23 =
              SQUARE(1. - 2 * (2 * n + 3) / (2 * n + 1) * mat_ptr->Lub_gpts[2] +
                     (2 * n + 3) / (2 * n + 1) * (n + 2) / (n + 1) * SQUARE(mat_ptr->Lub_gpts[2]));
          mat_ptr->Lub_wts[2] =
              (beta2 - f21 * beta0 - (f22 - f21) / (f12 - f11) * (beta1 - f11 * beta0)) /
              (f23 - f21 - (f22 - f21) / (f12 - f11) * (f13 - f11));
          mat_ptr->Lub_wts[1] =
              (beta1 - f11 * beta0 - (f13 - f11) * mat_ptr->Lub_wts[2]) / (f12 - f11);
          mat_ptr->Lub_wts[0] = beta0 - mat_ptr->Lub_wts[2] - mat_ptr->Lub_wts[1];
        } else if (mat_ptr->LubInt_NGP == 4) {
          GOMA_EH(GOMA_ERROR, "4-point POWER_LAW integration points not defined yet!");
        } else {
          GOMA_EH(GOMA_ERROR, "More than 4 POWER_LAW integration points not allowed yet!");
        }
      } else {
        GOMA_EH(model_read, "This Integration Model is not valid!");
      }
      ECHO(es, echo_file);
    }
  }

  /*
   * Input conditions for the structured porous shell equations
   * Added by SAR 2010-03-01
   */

  if ((pd_glob[mn]->gv[R_SHELL_SAT_CLOSED] == 1) || (have_shell_sat_open == 1) ||
      (have_shell_sat_open2 == 1)) {

    // Structured shell porosity
    model_read = look_for_mat_prop(imp, "Porous Shell Closed Porosity",
                                   &(mat_ptr->PorousShellClosedPorosityModel),
                                   &(mat_ptr->PorousShellClosedPorosity), NO_USER, NULL, model_name,
                                   SCALAR_INPUT, &NO_SPECIES, es);

    if (model_read == -1 && !strcmp(model_name, "EXTERNAL_FIELD")) {
      if (fscanf(imp, "%s", input) != 1) {
        GOMA_EH(
            GOMA_ERROR,
            "Expecting trailing keyword for Porous Shell Closed Porosity EXTERNAL_FIELD model.\n");
      }
      model_read = 1;

      // mat_ptr->u_PorousShellClosedPorosity_function_constants = alloc_dbl_1(1,0.0);

      ii = 0;
      for (j = 0; j < efv->Num_external_field; j++) {
        if (strcmp(efv->name[j], input) == 0) {
          ii = 1;
          mat_ptr->por_shell_closed_porosity_ext_field_index = j;
        }
      }
      if (ii == 0) {
        GOMA_EH(GOMA_ERROR,
                "Must activate external fields to use this PorousShellClosedPorosity model.  Field "
                "name needed for the EXTERNAL_FIELD Porous Shell Closed Porosity model");
      }
      mat_ptr->PorousShellClosedPorosityModel = EXTERNAL_FIELD;
      /* pick up scale factor for property */
      num_const = read_constants(imp, &(mat_ptr->u_PorousShellClosedPorosity_function_constants),
                                 NO_SPECIES);

      mat_ptr->len_u_PorousShellClosedPorosity_function_constants = num_const;
      if (num_const < 1) {
        sr = sprintf(err_msg, "Matl %s expected at least 1 constant for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Porous Shell Closed Porosity", "EXTERNAL_FIELD");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

    } else if (model_read == -1 && !strcmp(model_name, "MULTI_MODE")) {
      model_read = 1;
      mat_ptr->PorousShellClosedPorosityModel = MULTI_MODE;
      num_const = read_constants(imp, &(mat_ptr->u_PorousShellClosedPorosity_function_constants),
                                 NO_SPECIES);
      if (num_const < 2) {
        sr = sprintf(err_msg, "Matl %s needs at least 2 constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Porous Shell Closed Porosity", "MULTI_MODE");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_PorousShellClosedPorosity_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const,
                  mat_ptr->u_PorousShellClosedPorosity_function_constants);

    } else if (model_read == -1) {
      GOMA_EH(model_read, "Porous Shell Closed Porosity Invalid!");
    }
    ECHO(es, echo_file);

    // Structured shell height
    // This should eventually be removed as PRS changed to "Porous Shell Height"

    GOMA_WH(-1, " Change Porous Shell Closed Height to Porous Shell Height ");

    model_read =
        look_for_mat_prop(imp, "Porous Shell Height", &(mat_ptr->PorousShellClosedHeightModel),
                          &(mat_ptr->PorousShellClosedHeight), NO_USER, NULL, model_name,
                          SCALAR_INPUT, &NO_SPECIES, es);
    if (model_read == -1 && !strcmp(model_name, "EXTERNAL_FIELD")) {
      if (fscanf(imp, "%s", input) != 1) {
        GOMA_EH(GOMA_ERROR,
                "Expecting trailing keyword for Porous Shell Height EXTERNAL_FIELD model.\n");
      }
      model_read = 1;

      ii = 0;
      for (j = 0; j < efv->Num_external_field; j++) {
        if (strcmp(efv->name[j], input) == 0) {
          ii = 1;
          mat_ptr->por_shell_closed_height_ext_field_index = j;
        }
      }
      if (ii == 0) {
        GOMA_EH(GOMA_ERROR, "Must activate external fields to use this PorousShellheight model.  "
                            "Field name needed for the EXTERNAL_FIELD");
      }
      mat_ptr->PorousShellClosedHeightModel = EXTERNAL_FIELD;
      /* pick up scale factor for property */
      num_const =
          read_constants(imp, &(mat_ptr->u_PorousShellClosedHeight_function_constants), NO_SPECIES);

      mat_ptr->len_u_PorousShellClosedHeight_function_constants = num_const;
      if (num_const < 1) {
        sr = sprintf(err_msg, "Matl %s expected at least 1 constant for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Porous Shell Height", "EXTERNAL_FIELD");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

    } else if (model_read == -1) {
      GOMA_EH(model_read, "Porous Shell Height Invalid!");
    }
    ECHO(es, echo_file);

    // Pore radius distribution
    // This should eventually be removed as PRS changed to "Porous Shell Radius"
    GOMA_WH(-1, " Change Porous Shell Closed Radius to Porous Shell Radius ");

    // Pore radius distribution
    model_read =
        look_for_mat_prop(imp, "Porous Shell Radius", &(mat_ptr->PorousShellClosedRadiusModel),
                          &(mat_ptr->PorousShellClosedRadius), NO_USER, NULL, model_name,
                          SCALAR_INPUT, &NO_SPECIES, es);

    if (model_read == -1 && !strcmp(model_name, "EXTERNAL_FIELD")) {
      if (fscanf(imp, "%s", input) != 1) {
        GOMA_EH(GOMA_ERROR,
                "Expecting trailing keyword for Porous Shell Radius EXTERNAL_FIELD model.\n");
      }
      model_read = 1;

      ii = 0;
      for (j = 0; j < efv->Num_external_field; j++) {
        if (strcmp(efv->name[j], input) == 0) {
          ii = 1;
          mat_ptr->por_shell_closed_radius_ext_field_index = j;
        }
      }
      if (ii == 0) {
        GOMA_EH(GOMA_ERROR, "Must activate external fields to use this PorousShellRadius model.  "
                            "Field name needed for the EXTERNAL_FIELD");
      }

      mat_ptr->PorousShellClosedRadiusModel = EXTERNAL_FIELD;

      /* pick up scale factor for property */
      num_const =
          read_constants(imp, &(mat_ptr->u_PorousShellClosedRadius_function_constants), NO_SPECIES);

      mat_ptr->len_u_PorousShellClosedRadius_function_constants = num_const;

      if (num_const < 1) {
        sr = sprintf(err_msg, "Matl %s expected at least 1 constant for %s %s model.\n",
                     pd_glob[mn]->MaterialName, "Porous Shell Radius", "EXTERNAL_FIELD");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      model_read = 1;

    } else if (model_read == -1 && !strcmp(model_name, "MULTI_MODE")) {
      model_read = 1;
      mat_ptr->PorousShellClosedRadiusModel = MULTI_MODE;
      num_const =
          read_constants(imp, &(mat_ptr->u_PorousShellClosedRadius_function_constants), NO_SPECIES);
      if (num_const < mat_ptr->len_u_PorousShellClosedPorosity_function_constants - 1) {
        sr = sprintf(err_msg, "Matl %s needs at least %i constants for %s %s model.\n",
                     pd_glob[mn]->MaterialName,
                     mat_ptr->len_u_PorousShellClosedPorosity_function_constants - 1,
                     "Porous Shell Closed Radius", "MULTI_MODE");
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      mat_ptr->len_u_PorousShellClosedRadius_function_constants = num_const;
      SPF_DBL_VEC(endofstring(es), num_const,
                  mat_ptr->u_PorousShellClosedRadius_function_constants);

    } else if (model_read == -1) {
      GOMA_EH(model_read, "Porous Shell Radius Invalid!");
    }
    ECHO(es, echo_file);

    // Initial gas pressure
    model_read = look_for_mat_prop(
        imp, "Porous Shell Closed Gas Pressure", &(mat_ptr->PorousShellClosedP0Model),
        &(mat_ptr->PorousShellClosedP0), NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);
    if (model_read == -1) {
      GOMA_EH(model_read, "Porous Shell Closed Gas Pressure Invalid!");
    }
    ECHO(es, echo_file);

    // Atmospheric pressure
    model_read = look_for_mat_prop(imp, "Porous Shell Atmospheric Pressure",
                                   &(mat_ptr->PorousShellPatmModel), &(mat_ptr->PorousShellPatm),
                                   NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);
    if (model_read == -1) {
      GOMA_EH(model_read, "Porous Shell Atmospheric Pressure Invalid!");
    }
    ECHO(es, echo_file);

    // Reference pressure
    model_read = look_for_mat_prop(imp, "Porous Shell Reference Pressure",
                                   &(mat_ptr->PorousShellPrefModel), &(mat_ptr->PorousShellPref),
                                   NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);
    if (model_read == -1) {
      GOMA_EH(model_read, "Porous Shell Reference Pressure Invalid!");
    }
    ECHO(es, echo_file);

    // Cross Permeability for random shell media (not in plane). We will want
    // eventually an external field for this for variability
    model_read =
        look_for_mat_prop(imp, "Porous Shell Cross Permeability",
                          &(mat_ptr->PorousShellCrossKappaModel), &(mat_ptr->PorousShellCrossKappa),
                          NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);

    if (model_read == -1 && !strcmp(model_name, "EXTERNAL_FIELD")) {
      if (fscanf(imp, "%s", input) != 1) {
        GOMA_EH(GOMA_ERROR, "Expecting trailing keyword for Porous Shell Cross Permeability "
                            "EXTERNAL_FIELD model.\n");
      }
      ii = 0;
      for (j = 0; j < efv->Num_external_field; j++) {
        if (!strcmp(efv->name[j], input)) {
          ii = 1;
          mat_ptr->Xperm_external_field_index = j;
        }
      }
      if (ii == 0) {
        GOMA_EH(GOMA_ERROR, "Cannot match the Porous Shell Cross Permeability name with that in "
                            "the external field file");
      }
      mat_ptr->PorousShellCrossKappaModel = EXTERNAL_FIELD;

      /* pick up scale factor for property */
      num_const =
          read_constants(imp, &(mat_ptr->u_PorousShellCrossKappa_function_constants), NO_SPECIES);
      mat_ptr->len_u_PorousShellCrossKappa_function_constants = num_const;

      if (num_const < 1) {
        sr =
            sprintf(err_msg, "Matl %s expected at least 1 constant for %s %s model.\n",
                    pd_glob[mn]->MaterialName, "Porous Shell Cross Permeability", "EXTERNAL_FIELD");
        GOMA_EH(GOMA_ERROR, err_msg);
      }

    }

    else if (model_read == -1) {
      GOMA_EH(model_read, "Porous Shell Cross Permeability Invalid!");
    }
    ECHO(es, echo_file);

    // Initial Pore pressure in open porous media. this might be the same as
    // the dry state and may need to be overrriden when it is a restart.
    model_read = look_for_mat_prop(imp, "Porous Shell Initial Pore Pressure",
                                   &(mat_ptr->PorousShellInitPorePresModel),
                                   &(mat_ptr->PorousShellInitPorePres), NO_USER, NULL, model_name,
                                   SCALAR_INPUT, &NO_SPECIES, es);
    if (model_read == -1) {
      GOMA_EH(model_read, "Porous Shell Initial Pore Pressure model invalid");
    }
    ECHO(es, echo_file);

  } // End of structured porous shell inputs

  /*
   * Input conditions for the porous shell gas diffusion model
   * Added by SAR 2010-12-20
   */

  if (pd_glob[mn]->gv[R_SHELL_SAT_GASN] == 1) {

    // Gas diffusivity
    model_read = look_for_mat_prop(imp, "Porous Shell Gas Diffusivity",
                                   &(mat_ptr->PorousShellDiffusivityModel),
                                   &(mat_ptr->PorousShellDiffusivity), NO_USER, NULL, model_name,
                                   SCALAR_INPUT, &NO_SPECIES, es);
    if (model_read == -1) {
      GOMA_EH(model_read, "Porous Shell Gas Diffusivity Invalid!");
    }
    ECHO(es, echo_file);

    // Gas Temperature constant (RT)
    model_read = look_for_mat_prop(imp, "Porous Shell Gas Temperature Constant",
                                   &(mat_ptr->PorousShellRTModel), &(mat_ptr->PorousShellRT),
                                   NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);
    if (model_read == -1) {
      GOMA_EH(model_read, "Porous Shell Gas Temperature Constant Invalid!");
    }
    ECHO(es, echo_file);

    // Henry's Law Constant
    model_read = look_for_mat_prop(imp, "Porous Shell Henrys Law Constant",
                                   &(mat_ptr->PorousShellHenryModel), &(mat_ptr->PorousShellHenry),
                                   NO_USER, NULL, model_name, SCALAR_INPUT, &NO_SPECIES, es);
    if (model_read == -1) {
      GOMA_EH(model_read, "Porous Shell Henrys Law Invalid!");
    }
    ECHO(es, echo_file);

  } // End of porous shell gas diffusion constants

  /*
   * Inputs specific to thin film multiphase flow density and viscosity calculations
    // So far, the density card only pertains to the gas model.
    // The incompressible liquid density does not matter.
   */
  if (pd_glob[mn]->gv[R_TFMP_MASS] || pd_glob[mn]->gv[R_TFMP_BOUND]) {
    char input[MAX_CHAR_IN_INPUT] = "zilch\0";
    strcpy(search_string, "Thin Film Multiphase Density");
    model_read = look_for_optional(imp, search_string, input, '=');

    if (model_read == 1) {
      if (fscanf(imp, "%s", model_name) != 1) {
        sr = sprintf(err_msg, "Error reading model name string in material file, property %s",
                     search_string);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      SPF(es, "%s = %s", search_string, model_name);
      if (model_read == 1 && !strcmp(model_name, "CONSTANT")) {
        model_read = 1;
        mat_ptr->tfmp_density_model = CONSTANT;
        num_const = read_constants(imp, &(mat_ptr->tfmp_density_const), NO_SPECIES);

        SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->tfmp_density_const);
        mat_ptr->len_tfmp_density_const = num_const;
        if (num_const == 1) {
          mat_ptr->len_tfmp_density_const = 4;
          safe_free(mat_ptr->tfmp_density_const);
          mat_ptr->tfmp_density_const = alloc_dbl_1(4, 0.0);
          // make sure reasonable values are here, to prevent divide by zero and provide a uniform
          // ambient pressure value.
          mat_ptr->tfmp_density_const[1] = 1.0;
          mat_ptr->tfmp_density_const[2] = 1.0;
          mat_ptr->tfmp_density_const[3] = 0.0;
        }
      }
      if (model_read == 1 && !strcmp(model_name, "IDEAL_GAS")) {
        model_read = 1;
        mat_ptr->tfmp_density_model = IDEAL_GAS;
        num_const = read_constants(imp, &(mat_ptr->tfmp_density_const), NO_SPECIES);

        SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->tfmp_density_const);
        mat_ptr->len_tfmp_density_const = num_const;
        if (num_const != 4) {
          GOMA_EH(GOMA_ERROR, "The IDEAL_GAS model requires 4 values: molecular weight of gas, "
                              "universal gas constant, temperature[const], ambient pressure.");
        }
      }
    } else {
      GOMA_EH(GOMA_ERROR, "You must use the \"Thin Film Multiphase Density\" card to specify the "
                          "gas density for the tfmp equations.");
    }
    ECHO(es, echo_file);
  }

  if (pd_glob[mn]->gv[R_TFMP_MASS] || pd_glob[mn]->gv[R_TFMP_BOUND]) {
    char input[MAX_CHAR_IN_INPUT] = "zilch\0";
    strcpy(search_string, "Thin Film Multiphase Viscosity");
    model_read = look_for_optional(imp, search_string, input, '=');

    if (model_read == 1) {
      if (fscanf(imp, "%s", model_name) != 1) {
        sr = sprintf(err_msg, "Error reading model name string in material file, property %s",
                     search_string);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      SPF(es, "%s = %s", search_string, model_name);
      if (!strcmp(model_name, "CONSTANT")) {
        model_read = 1;
        mat_ptr->tfmp_viscosity_model = CONSTANT;
        num_const = read_constants(imp, &(mat_ptr->tfmp_viscosity_const), NO_SPECIES);

        SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->tfmp_viscosity_const);
        mat_ptr->len_tfmp_viscosity_const = num_const;
        if (num_const != 2) {
          sr = sprintf(err_msg, "Wrong number of parameters on property, %s", search_string);
          GOMA_EH(GOMA_ERROR, err_msg);
        }
      }
    } else {
      GOMA_EH(GOMA_ERROR, "No default fluid viscosities, to specify them use the \"Thin Film "
                          "Multiphase Viscosity\" card.");
    }
  }
  if (pd_glob[mn]->gv[R_TFMP_MASS] || pd_glob[mn]->gv[R_TFMP_BOUND]) {
    char input[MAX_CHAR_IN_INPUT] = "zilch\0";
    strcpy(search_string, "Thin Film Multiphase Diffusivity Model");
    model_read = look_for_optional(imp, search_string, input, '=');

    if (model_read == 1) {
      GOMA_WH(-1, "\"Thin Film Multiphase Diffusivity Model\" adds a diffusion term to the liquid "
                  "volume balance equation for the express purpose of numerical stabilization, "
                  "proceed with caution.");
      if (fscanf(imp, "%s", model_name) != 1) {
        sr = sprintf(err_msg, "Error reading model name string in material file, property %s",
                     search_string);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      SPF(es, "%s = %s", search_string, model_name);
      if (!strcmp(model_name, "CONSTANT")) {
        model_read = 1;

        mat_ptr->tfmp_diff_model = CONSTANT;
        mat_ptr->tfmp_diff_const = alloc_dbl_1(1, 0.0);
        if (fscanf(imp, "%lg", mat_ptr->tfmp_diff_const) != 1) {
          sr = sprintf(err_msg, "Wrong number of constants in material file, property %s",
                       search_string);
          GOMA_EH(GOMA_ERROR, err_msg);
        }

        mat_ptr->len_tfmp_diff_const = 1;
        SPF_DBL_VEC(endofstring(es), 1, mat_ptr->tfmp_diff_const);
      } else if (!strcmp(model_name, "PIECEWISE")) {
        mat_ptr->tfmp_diff_model = PIECEWISE;

        mat_ptr->len_tfmp_diff_const = read_constants(imp, &(mat_ptr->tfmp_diff_const), NO_SPECIES);
        SPF_DBL_VEC(endofstring(es), mat_ptr->len_tfmp_diff_const, mat_ptr->tfmp_diff_const);
      }
      ECHO(es, echo_file);
    } else {
      GOMA_WH(-1, "If you're having trouble try adding some numerical diffusion with the \"Thin "
                  "Film Multiphase Diffusivity Model\" material property.");
    }
  }

  if (pd_glob[mn]->gv[R_TFMP_BOUND]) {
    char input[MAX_CHAR_IN_INPUT] = "zilch\0";
    strcpy(search_string, "Thin Film Multiphase Dissolution Model");
    model_read = look_for_optional(imp, search_string, input, '=');
    if (model_read == 1) {
      if (fscanf(imp, "%s", model_name) != 1) {
        sr = sprintf(err_msg, "Error reading model name string in material file, property %s",
                     search_string);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      SPF(es, "%s = %s", search_string, model_name);
      num_const = read_constants(imp, &(mat_ptr->tfmp_dissolution_const), NO_SPECIES);
      if (!strcmp(model_name, "SQUARE")) {
        mat_ptr->tfmp_dissolution_model = TFMP_SQUARE;
        if (num_const != 3) {
          sr = sprintf(err_msg, "Error property %s only supports 3 input values.", search_string);
          GOMA_EH(GOMA_ERROR, err_msg);
        }
        mat_ptr->len_tfmp_dissolution_const = num_const;
        SPF_DBL_VEC(endofstring(es), mat_ptr->len_tfmp_dissolution_const,
                    mat_ptr->tfmp_dissolution_const);
      }
      ECHO(es, echo_file);
    } else {
      mat_ptr->tfmp_dissolution_model = NO_MODEL;
      mat_ptr->len_tfmp_dissolution_const = 0;
      GOMA_WH(-1, "By default, dissolution is inactive. Use \"Thin Film Multiphase Dissolution "
                  "Model\" to activate it.");
    }
  }

  if (pd_glob[mn]->gv[R_TFMP_MASS] || pd_glob[mn]->gv[R_TFMP_BOUND]) {
    char input[MAX_CHAR_IN_INPUT] = "zilch\0";
    strcpy(search_string, "Thin Film Multiphase Relative Permeability Model");
    model_read = look_for_optional(imp, search_string, input, '=');
    if (model_read == 1) {
      if (fscanf(imp, "%s", model_name) != 1) {
        sr = sprintf(err_msg, "Error reading model name string in material file, property %s",
                     search_string);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      SPF(es, "%s = %s", search_string, model_name);
      if (!strcmp(model_name, "LEVER")) {
        mat_ptr->tfmp_rel_perm_model = LEVER;
        mat_ptr->len_tfmp_rel_perm_const = 0;
      } else if (!strcmp(model_name, "SATURATION")) {
        mat_ptr->tfmp_rel_perm_model = SATURATION;
        mat_ptr->len_tfmp_rel_perm_const = 0;
      } else if (!strcmp(model_name, "PIECEWISE")) {
        mat_ptr->tfmp_rel_perm_model = PIECEWISE;
        mat_ptr->len_tfmp_rel_perm_const =
            read_constants(imp, &mat_ptr->tfmp_rel_perm_const, NO_SPECIES);
        SPF_DBL_VEC(endofstring(es), mat_ptr->len_tfmp_rel_perm_const,
                    mat_ptr->tfmp_rel_perm_const);
      } else {
        sr = sprintf(err_msg, "Invalid model name string in material file, property %s",
                     search_string);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      ECHO(es, echo_file);
    } else {
      GOMA_EH(GOMA_ERROR, "There are no defaults for \"Thin Film Multiphase Relative Permeability "
                          "Model\". You must set them.");
    }
  }

  if (pd_glob[mn]->gv[R_TFMP_BOUND] || pd_glob[mn]->gv[R_TFMP_MASS]) {
    strcpy(search_string, "Thin Film Multiphase Mass Lumping");
    model_read = look_for_optional(imp, search_string, input, '=');
    if (model_read == 1) {
      if (fscanf(imp, "%s", model_name) != 1) {
        sr = sprintf(err_msg, "Error reading model name string in material file, property %s",
                     search_string);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      SPF(es, "%s = %s", search_string, model_name);
      if (!strcasecmp(model_name, "yes") || !strcasecmp(model_name, "true")) {
        mat_ptr->tfmp_mass_lump = TRUE;
        SPF(es, "%s = %s", search_string, "TRUE");
      } else if (!strcasecmp(model_name, "no") || !strcasecmp(model_name, "false")) {
        mat_ptr->tfmp_mass_lump = FALSE;
        SPF(es, "%s = %s", search_string, "FALSE");
      } else {
        GOMA_EH(GOMA_ERROR,
                "Thin Film Multiphase Mass Lumping must be set to TRUE, YES, FALSE, or NO");
      }
    } else {
      GOMA_WH(-1, "Mass lumping is on by default.");
      mat_ptr->tfmp_mass_lump = TRUE;
      SPF(es, "%s = %s", search_string, "TRUE");
    }
    ECHO(es, echo_file);
  }

  if (pd_glob[mn]->gv[R_TFMP_BOUND]) {
    char input[MAX_CHAR_IN_INPUT] = "zilch\0";
    strcpy(model_name, "\0");
    strcpy(search_string, "Thin Film Multiphase Clipping");
    model_read = look_for_optional(imp, search_string, input, '=');
    if (model_read == 1) {
      if (fscanf(imp, "%s", model_name) != 1) {
        sr = sprintf(err_msg, "Error reading model name string in material file, property %s",
                     search_string);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      SPF(es, "%s = %s", search_string, model_name);
      //      GOMA_EH(GOMA_ERROR, model_name);
      if (!strcasecmp(model_name, "yes") || !strcasecmp(model_name, "true")) {
        mat_ptr->tfmp_clipping = TRUE;
        if (fscanf(imp, "%lg", &(mat_ptr->tfmp_clip_strength)) != 1) {
          sr = sprintf(err_msg, "Wrong number of constants in material file, property %s",
                       search_string);
          GOMA_EH(GOMA_ERROR, err_msg);
        }
        SPF(endofstring(es), " %.4g", mat_ptr->tfmp_clip_strength);

      } else if (!strcasecmp(model_name, "no") || !strcasecmp(model_name, "false")) {
        mat_ptr->tfmp_clipping = FALSE;
        GOMA_WH(-1, "Spurious oscillations can be mitigated with the \"Thin Film Multiphase "
                    "Clipping\" material property.");
      } else {
        SPF(err_msg, "Syntax error or invalid model for %s\n", search_string);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

    } else {
      mat_ptr->tfmp_clipping = FALSE;
      SPF(es, "%s = %s", search_string, "FALSE");
      GOMA_WH(-1, "\"Thin Film Multiphase Clipping\" is off by default.");
    }
    ECHO(es, echo_file);
  }

  if (pd_glob[mn]->gv[R_TFMP_MASS] || pd_glob[mn]->gv[R_TFMP_BOUND]) {
    char input[MAX_CHAR_IN_INPUT] = "zilch\0";
    strcpy(search_string, "Thin Film Multiphase Drop Lattice");
    model_read = look_for_optional(imp, search_string, input, '=');
    if (model_read == 1) {
      if (fscanf(imp, "%s", model_name) != 1) {
        sr = sprintf(err_msg, "Error reading model name string in material file, property %s",
                     search_string);
        GOMA_EH(GOMA_ERROR, err_msg);
      }
      SPF(es, "%s = %s", search_string, model_name);
      if (!strcmp(model_name, "SQUARE")) {
        mat_ptr->tfmp_drop_lattice_model = TFMP_SQUARE;
        num_const = read_constants(imp, &(mat_ptr->tfmp_drop_lattice_const), NO_SPECIES);
        if (num_const != 2) {
          GOMA_EH(GOMA_ERROR, "Thin Film Multiphase Drop Lattice 'SQUARE' requires two and only "
                              "two input values, lambda and Vd");
        }
        mat_ptr->len_tfmp_drop_lattice_const = num_const;
        SPF_DBL_VEC(endofstring(es), num_const, mat_ptr->tfmp_drop_lattice_const);
      }
    } else {
      GOMA_WH(-1, "By default, \"Thin Film Multiphase Drop Lattice\" is a square lattice with "
                  "spacing of 160 microns and 6 pL drop volume, units in cgs.");
    }
    ECHO(es, echo_file);
  }

  /* check for roller-web gap thickness model */
  if (pd_glob[mn]->gv[R_TFMP_BOUND]) {
    char input[MAX_CHAR_IN_INPUT] = "zilch\0";
    strcpy(search_string, "Elastohydrodynamic Lubrication Gap Model");
    model_read = look_for_optional(imp, search_string, input, '=');

    if (model_read == 1) {
      if (fscanf(imp, "%s", model_name) != 1) {
        sr = sprintf(err_msg, "Error reading model name string in material file, property %s",
                     search_string);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      SPF(es, "%s = %s", search_string, model_name);
      if (model_read == 1 && !strcmp(model_name, "NDOTD")) {
        model_read = 1;
        mat_ptr->ehl_gap_model = GM_NDOTD;

      } else

          if (model_read == 1 && !strcmp(model_name, "RADIAL")) {
        model_read = 1;
        mat_ptr->ehl_gap_model = GM_RADIAL;

      }

      else {
        // default is simple
        mat_ptr->ehl_gap_model = GM_RADIAL;
        SPF(es, "%s = %s", search_string, "RADIAL");
      }

      ECHO(es, echo_file);
    }
  }

  /* check for roller-web normal calculation method (only needed if gap model is NDOTD) */
  if (pd_glob[mn]->gv[R_MESH1] && pd_glob[mn]->gv[R_TFMP_BOUND] &&
      mat_ptr->ehl_gap_model == GM_NDOTD) {
    char input[MAX_CHAR_IN_INPUT] = "zilch\0";
    strcpy(search_string, "Elastohydrodynamic Lubrication Normal Calculation Method");
    model_read = look_for_optional(imp, search_string, input, '=');

    if (model_read == 1) {
      if (fscanf(imp, "%s", model_name) != 1) {
        sr = sprintf(err_msg, "Error reading model name string in material file, property %s",
                     search_string);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      SPF(es, "%s = %s", search_string, model_name);
      if (model_read == 1 && !strcmp(model_name, "MAPPING")) {
        model_read = 1;
        mat_ptr->ehl_normal_method = NCM_MAPPING;

      } else

          if (model_read == 1 && !strcmp(model_name, "SIK_S_WEB")) {
        model_read = 1;
        mat_ptr->ehl_normal_method = NCM_PRIMITIVE_S_WEB;

      } else if (model_read == 1 && !strcmp(model_name, "SIK_S_ROLLER")) {
        model_read = 1;
        mat_ptr->ehl_normal_method = NCM_PRIMITIVE_S_ROLLER;

      } else if (model_read == 1 && !strcmp(model_name, "SIK_XY")) {
        model_read = 1;
        mat_ptr->ehl_normal_method = NCM_PRIMITIVE_XY;

      }

      else {
        // default is normal of roller
        mat_ptr->ehl_normal_method = NCM_PRIMITIVE_S_ROLLER;
        SPF(es, "%s = %s", search_string, "SIK_S_ROLLER");
      }

      ECHO(es, echo_file);
    }
  }

  /* check for 2d bar integration kind */
  if (pd_glob[mn]->gv[R_TFMP_BOUND]) {
    char input[MAX_CHAR_IN_INPUT] = "zilch\0";
    strcpy(search_string, "Elastohydrodynamic Lubrication Shell Integration Kind");
    model_read = look_for_optional(imp, search_string, input, '=');

    if (model_read == 1) {
      if (fscanf(imp, "%s", model_name) != 1) {
        sr = sprintf(err_msg, "Error reading model name string in material file, property %s",
                     search_string);
        GOMA_EH(GOMA_ERROR, err_msg);
      }

      SPF(es, "%s = %s", search_string, model_name);
      if (model_read == 1 && !strcmp(model_name, "S")) {
        model_read = 1;
        mat_ptr->ehl_integration_kind = SIK_S;

      } else

          if (model_read == 1 && !strcmp(model_name, "XY")) {
        model_read = 1;
        mat_ptr->ehl_integration_kind = SIK_XY;

      }

      else {
        // default is XY
        mat_ptr->ehl_integration_kind = SIK_XY;
        SPF(es, "%s = %s", search_string, "XY");
      }

      ECHO(es, echo_file);
    }
  }

  /*********************************************************************/

  /*********************************************************************/
  /*********************************************************************/
  /*********************************************************************/
  /*
   *  Check the input for consistency
   */
  NO_SPECIES = 0;
  for (i = 0; i < Num_Var_Init_Mat[mn]; i++) {
    if (Var_init_mat[mn][i].var == MASS_FRACTION)
      NO_SPECIES++;
  }
  if (NO_SPECIES > pd_ptr->Num_Species) {
    sprintf(Err_Msg,
            "Attempt to initialize %d species in matl %s with only"
            "%d species active!",
            NO_SPECIES, pd_glob[mn]->MaterialName, pd_ptr->Num_Species_Eqn);
    GOMA_EH(GOMA_ERROR, Err_Msg);
  }

  /*
   *  Check for consistency in the specifications of the Species Variables
   *  now that we have all of the information about the constituitive
   *  modeling for the material
   */

  if (mat_ptr->Species_Var_Type == SPECIES_UNDEFINED_FORM) {
    mat_ptr->Species_Var_Type = pd_ptr->Species_Var_Type;
  }

  /*
   *  Calculate the value of Dropped_last_sepecies_eqn from the values
   *  of Num_Species_Eqn and Num_Species. We will assume that if they
   *  are one apart, then the last species continuity equation has
   *  been dropped from the equation system
   */
  if (mat_ptr->Num_Species_Eqn == mat_ptr->Num_Species - 1) {
    mat_ptr->Dropped_Last_Species_Eqn = TRUE;
  } else {
    mat_ptr->Dropped_Last_Species_Eqn = FALSE;
  }
  /*
   * HKM -> Possible spot in the code to include a section on
   *        the consistency of mass fraction values, i.e., they
   *        should sum to one if Species_Var_Type is a certain
   *        type.
   *        Also, we should check the range of the input ktypes.
   */
}
/*  END of rd_mp_specs -- read material properties specifications  */
