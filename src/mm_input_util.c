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
 *$Id: mm_input_util.c,v 5.6 2010-04-07 22:27:00 prschun Exp $
 */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************
 *
 *     Functions in this file
 *  -----------------------------------
 *    NAME                      DESCRIPTION
 * --------------------------------------------------------------------------
 *  read_int                    Reads and interprets an int from a file pointer
 *  read_dbl                    Reads and interprets a dbl from a file pointer
 *  read_1_boolean              Reads and interprets a string representing
 *                              a boolean from a file pointer
 *  variable_string_to_int      Translates a string representation of a
 *                              variable into an integer number
 *  species_type_str_to_int     Converts species type string into an int
 *  species_type_int_to_str     Converts an int into a species type string
 *  interpret_int               Translates a string into an int.
 *  interpret_double            Translates a string into a double
 *  indentify_species_ID_string Translates a string into a species ID.
 *                              Both ints and strings are used to identify
 *                              the string.
 *  stokenize                   Tokenizes a string
 *  tokenize_by_whsp            Tokenizes a string using white space defn.
 *  fillTokStruct               Fills up a Token structure given a string
 *  in_char_list                Searches a list of strings for a match
 *  read_string                 Searches forward for a string in a
 *                              file descriptor
 *  strip                       strips leading and trailing whitespace
 *
 *****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mm_eh.h"
#include "mm_input.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "rf_allo.h"
#include "rf_bc_const.h"
#include "rf_fem_const.h"
#include "rf_io_const.h"
#include "std.h"

/* Read up to n doubles from a file up to a newline


   ifp is the file ptr where the current line should be looked at

   n is the number of doubles to look for

   array is a double array of at least length n

   return number of doubles found and read
*/
int look_for_n_doubles(FILE *ifp, int n, double *array) {
  char strbuf[MAX_INPUT_LINE_LENGTH];
  char *fgerror;
  int error;
  int i;
  int buf_offset;
  int read_count;
  fgerror = fgets(strbuf, MAX_INPUT_LINE_LENGTH, ifp);
  if (fgerror == NULL) {
    GOMA_EH(GOMA_ERROR, "Error reading line");
  }

  buf_offset = 0;
  for (i = 0; i < n; i++) {
    error = sscanf(strbuf + buf_offset, "%lf%n", &array[i], &read_count);
    buf_offset += read_count;
    if (error != 1)
      break;
  }

  return i;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int look_for_optional_int(FILE *ifp, const char *stringMatch, int *intValue, const int intDefault)

/*************************************************************************
 *
 * look_for_optional_int():
 *
 * Scan the input file (reading in strings according to
 * 'read_string(ifp,)' specifications)
 * until the character pattern in 'stringMatch' is matched.
 * It then reads the first int after the '=' sign and stores it in
 * intValue. The rest of the line is discarded.
 * If 'stringMatch' is not matched, return a value of -1, and assign
 * intValue to its default value if intDefault isn't set to INT_NODEFAULT.
 * If 'strinMatch' is matched, return 1.
 *
 * This search starts at the current position in the input
 * file and searches to the end until it finds the input string.
 * If this routine can't find the input string, it returns ifp to the
 * file position where the search was started. If it finds the
 * string pattern, it leaves the file pointer positioned after the
 * first token on the matching line after the "=" sign.
 *
 * Parameter list:
 *  ifp      == pointer to file "input"
 *  string   == contains string pattern to be matched.
 *  intValue == pointer to the location which will receive the
 *              integer value
 *  intDefault == Default value of the intValue if the stringMatch
 *                can't be found. if intDefault == INT_NODEFAULT
 *                then the default is applied and the incomming
 *                value of *intValue is preserved.
 *
 * Return Value:
 *    1 if string is matched
 *   -1 if no match
 *
 *************************************************************************/
{
  int retn;
  char input[MAX_CHAR_IN_INPUT];
  if (intValue == NULL) {
    GOMA_EH(GOMA_ERROR, "look_for_optional_int ERROR: Incoming value of intValue is NULL");
  }
  if (stringMatch == NULL) {
    GOMA_EH(GOMA_ERROR, "look_for_optional_int ERROR: Incoming value of stringMatch is NULL");
  }
  retn = look_forward_optional(ifp, stringMatch, input, '=');
  if (retn == 1) {
    *intValue = read_int(ifp, stringMatch);
  } else {
    if (intDefault != INT_NOINIT) {
      *intValue = intDefault;
    }
  }
  return retn;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int read_int(FILE *ifp, const char *intName)

/***********************************************************************
 *
 * read_int:
 *
 *      Reads an integer from the file corresponding to the
 * file pointer ifp. Leaves the file pointer positioned after the
 * integer read. If an error does occur during the fscanf operation,
 * an error string is created an the error hanlder EH, is called
 * with a fatal error condition indicated.
 *
 * Parameters:
 *  ifp          File Pointer
 *  intName      character string identifying the integer
 *               (only used for fatal error messages)
 *
 * Return:
 *   Function returns the value of the integer read
 **********************************************************************/
{
  int intTmp, len;
  char *errString;
  if (fscanf(ifp, "%d", &intTmp) != 1) {
    len = strlen(intName) + 50;
    errString = (char *)smalloc(len);
    sprintf(errString, "read_int: Expected to read an int for \"%s\"", intName);
    GOMA_EH(GOMA_ERROR, errString);
    safer_free((void **)&errString);
  }
  return intTmp;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double read_dbl(FILE *ifp, const char *dblName)

/***********************************************************************
 *
 * read_dbl:
 *
 *      Reads a double from the file corresponding to the
 * file pointer ifp. Leaves the file pointer positioned after the
 * read. If an error does occur during the fscanf operation,
 * an error string is created and an the error hanlder EH, is called
 * with a fatal error condition indicated.
 *
 * Parameters:
 *  ifp          File Pointer
 *  dblName      character string identifying the double
 *               (only used for fatal error messages)
 *
 * Return:
 *   Function returns the value of the double read
 **********************************************************************/
{
  double dblTmp;
  int len;
  char *errString;
  if (fscanf(ifp, "%le", &dblTmp) != 1) {
    len = strlen(dblName) + 50;
    errString = (char *)smalloc(len);
    (void)sprintf(errString, "read_dbl: Expected to read a double for \"%s\"", dblName);
    GOMA_EH(-1, errString);
    safer_free((void **)&errString);
  }
  return dblTmp;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int read_1_boolean(FILE *ifp, const char *boolName)

/***********************************************************************
 *
 * read_1_boolean:
 *
 *      Reads a boolean from the file corresponding to the
 * file pointer ifp. Leaves the file pointer positioned after the
 * end of the line on which the boolean resides.             <- NOTE ====
 * If an error does occur,
 * an error string is created and an the error hanlder, EH, is called
 * with a fatal error condition indicated.
 *
 * Only the first character is read. t, T, or 1 is interpreted as true.
 * f, F, or 0 is interpreted as false. Every other character produces
 * an error message
 *
 * Parameters:
 *  ifp          File Pointer
 *  boolName      character string identifying the bool
 *               (only used for fatal error messages)
 *
 * Return:
 *   Function returns either 1 for true, or 0 for false.
 **********************************************************************/
{
  int success, retn = 0;
  char input[MAX_CHAR_IN_INPUT];
  int len;
  char *errString;
  success = read_string(ifp, input, '\n');
  if (success > 0) {
    success = strip(input);
    if (success > 0) {
      if (input[0] == 't' || input[0] == 'T' || input[0] == '1') {
        retn = TRUE;
      } else if (input[0] == 'f' || input[0] == 'F' || input[0] == '0') {
        retn = FALSE;
      } else {
        success = -2;
      }
    }
  }
  if (success <= 0) {
    len = strlen(boolName) + 80;
    errString = (char *)smalloc(len);
    (void)sprintf(errString,
                  "read_1_bool: Expected to read a boolean for \"%s\" \n\t\tfound: \"%s\"s\n",
                  boolName, input);
    GOMA_EH(-1, errString);
    safer_free((void **)&errString);
  }
  return retn;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int variable_string_to_int(const char *input, const char *err_string)

/***************************************************************************
 *
 * variable_string_to_int:
 *
 *  Does the conversion between the string representation of a variable
 *  and the number representation. All of the valid strings and numbers
 *  in rf_fem_const.h are included here.
 *  Right now The Initialize keyword section is the only place where this
 *  function is utilizied. However, it should find uses in multiple
 *  places.
 *
 *  Parameters:
 * ---------------
 *  input      - Character string representing the variable (and various
 *               instantiations of the variable as well ala distinguishing
 *               between mass, mole, and volume fractions for instance).
 *  err_string - additional error string to be printed out in case
 *               no match is found (use NULL if no additional io is
 *               required
 *  Return
 * -------------
 *  zero and positive values -> successful match has been made
 *  -1: special match to R_ANYTHING
 *  -2: no match was found.
 *
 *  IMPLEMENTATION NOTE
 *      Broke up the if statement into sections because of internal
 *      compiler parsing problems on HPUX. The compiler didn't like
 *      all those else if statements. A better approach would probably
 *      be to put all of these variables in a static array, which has
 *      been initialized in the .h file. Then, this file would become
 *      very short.
 **************************************************************************/
{
  int var = -10;
  if (!strcmp(input, "VELOCITY1"))
    var = VELOCITY1;
  else if (!strcmp(input, "VELOCITY2"))
    var = VELOCITY2;
  else if (!strcmp(input, "VELOCITY3"))
    var = VELOCITY3;
  else if (!strcmp(input, "TEMPERATURE"))
    var = TEMPERATURE;
  else if (!strcmp(input, "MASS_FRACTION"))
    var = MASS_FRACTION;
  else if (!strcmp(input, "MESH_DISPLACEMENT1"))
    var = MESH_DISPLACEMENT1;
  else if (!strcmp(input, "MESH_DISPLACEMENT2"))
    var = MESH_DISPLACEMENT2;
  else if (!strcmp(input, "MESH_DISPLACEMENT3"))
    var = MESH_DISPLACEMENT3;
  else if (!strcmp(input, "SURFACE"))
    var = SURFACE;
  else if (!strcmp(input, "PRESSURE"))
    var = PRESSURE;

  else if (!strcmp(input, "POLYMER_STRESS11"))
    var = POLYMER_STRESS11;
  else if (!strcmp(input, "POLYMER_STRESS12"))
    var = POLYMER_STRESS12;
  else if (!strcmp(input, "POLYMER_STRESS13"))
    var = POLYMER_STRESS13;
  else if (!strcmp(input, "POLYMER_STRESS22"))
    var = POLYMER_STRESS22;
  else if (!strcmp(input, "POLYMER_STRESS23"))
    var = POLYMER_STRESS23;
  else if (!strcmp(input, "POLYMER_STRESS33"))
    var = POLYMER_STRESS33;

  else if (!strcmp(input, "SOLID_DISPLACEMENT1"))
    var = SOLID_DISPLACEMENT1;
  else if (!strcmp(input, "SOLID_DISPLACEMENT2"))
    var = SOLID_DISPLACEMENT2;
  else if (!strcmp(input, "SOLID_DISPLACEMENT3"))
    var = SOLID_DISPLACEMENT3;

  else if (!strcmp(input, "VELOCITY_GRADIENT11"))
    var = VELOCITY_GRADIENT11;
  else if (!strcmp(input, "VELOCITY_GRADIENT12"))
    var = VELOCITY_GRADIENT12;
  else if (!strcmp(input, "VELOCITY_GRADIENT13"))
    var = VELOCITY_GRADIENT13;
  else if (!strcmp(input, "VELOCITY_GRADIENT21"))
    var = VELOCITY_GRADIENT21;
  else if (!strcmp(input, "VELOCITY_GRADIENT22"))
    var = VELOCITY_GRADIENT22;
  else if (!strcmp(input, "VELOCITY_GRADIENT23"))
    var = VELOCITY_GRADIENT23;
  else if (!strcmp(input, "VELOCITY_GRADIENT31"))
    var = VELOCITY_GRADIENT31;
  else if (!strcmp(input, "VELOCITY_GRADIENT32"))
    var = VELOCITY_GRADIENT32;
  else if (!strcmp(input, "VELOCITY_GRADIENT33"))
    var = VELOCITY_GRADIENT33;

  else if (!strcmp(input, "VOLTAGE"))
    var = VOLTAGE;
  else if (!strcmp(input, "FILL"))
    var = FILL;
  else if (!strcmp(input, "LS"))
    var = FILL;
  else if (!strcmp(input, "SHEAR_RATE"))
    var = SHEAR_RATE;

  else if (!strcmp(input, "PVELOCITY1"))
    var = PVELOCITY1;
  else if (!strcmp(input, "PVELOCITY2"))
    var = PVELOCITY2;
  else if (!strcmp(input, "PVELOCITY3"))
    var = PVELOCITY3;

  /*
   * Kluge to break up large if block. Problems with HP compiler
   */
  if (var != -10)
    return var;
  /*
   * Multimode viscoelasticity
   */
  if (!strcmp(input, "POLYMER_STRESS11_1"))
    var = POLYMER_STRESS11_1;
  else if (!strcmp(input, "POLYMER_STRESS12_1"))
    var = POLYMER_STRESS12_1;
  else if (!strcmp(input, "POLYMER_STRESS22_1"))
    var = POLYMER_STRESS22_1;
  else if (!strcmp(input, "POLYMER_STRESS13_1"))
    var = POLYMER_STRESS13_1;
  else if (!strcmp(input, "POLYMER_STRESS23_1"))
    var = POLYMER_STRESS23_1;
  else if (!strcmp(input, "POLYMER_STRESS33_1"))
    var = POLYMER_STRESS33_1;

  else if (!strcmp(input, "POLYMER_STRESS11_2"))
    var = POLYMER_STRESS11_2;
  else if (!strcmp(input, "POLYMER_STRESS12_2"))
    var = POLYMER_STRESS12_2;
  else if (!strcmp(input, "POLYMER_STRESS22_2"))
    var = POLYMER_STRESS22_2;
  else if (!strcmp(input, "POLYMER_STRESS13_2"))
    var = POLYMER_STRESS13_2;
  else if (!strcmp(input, "POLYMER_STRESS23_2"))
    var = POLYMER_STRESS23_2;
  else if (!strcmp(input, "POLYMER_STRESS33_2"))
    var = POLYMER_STRESS33_2;

  else if (!strcmp(input, "POLYMER_STRESS11_3"))
    var = POLYMER_STRESS11_3;
  else if (!strcmp(input, "POLYMER_STRESS12_3"))
    var = POLYMER_STRESS12_3;
  else if (!strcmp(input, "POLYMER_STRESS22_3"))
    var = POLYMER_STRESS22_3;
  else if (!strcmp(input, "POLYMER_STRESS13_3"))
    var = POLYMER_STRESS13_3;
  else if (!strcmp(input, "POLYMER_STRESS23_3"))
    var = POLYMER_STRESS23_3;
  else if (!strcmp(input, "POLYMER_STRESS33_3"))
    var = POLYMER_STRESS33_3;

  else if (!strcmp(input, "POLYMER_STRESS11_4"))
    var = POLYMER_STRESS11_4;
  else if (!strcmp(input, "POLYMER_STRESS12_4"))
    var = POLYMER_STRESS12_4;
  else if (!strcmp(input, "POLYMER_STRESS22_4"))
    var = POLYMER_STRESS22_4;
  else if (!strcmp(input, "POLYMER_STRESS13_4"))
    var = POLYMER_STRESS13_4;
  else if (!strcmp(input, "POLYMER_STRESS23_4"))
    var = POLYMER_STRESS23_4;
  else if (!strcmp(input, "POLYMER_STRESS33_4"))
    var = POLYMER_STRESS33_4;

  else if (!strcmp(input, "POLYMER_STRESS11_5"))
    var = POLYMER_STRESS11_5;
  else if (!strcmp(input, "POLYMER_STRESS12_5"))
    var = POLYMER_STRESS12_5;
  else if (!strcmp(input, "POLYMER_STRESS22_5"))
    var = POLYMER_STRESS22_5;
  else if (!strcmp(input, "POLYMER_STRESS13_5"))
    var = POLYMER_STRESS13_5;
  else if (!strcmp(input, "POLYMER_STRESS23_5"))
    var = POLYMER_STRESS23_5;
  else if (!strcmp(input, "POLYMER_STRESS33_5"))
    var = POLYMER_STRESS33_5;

  else if (!strcmp(input, "POLYMER_STRESS11_6"))
    var = POLYMER_STRESS11_6;
  else if (!strcmp(input, "POLYMER_STRESS12_6"))
    var = POLYMER_STRESS12_6;
  else if (!strcmp(input, "POLYMER_STRESS22_6"))
    var = POLYMER_STRESS22_6;
  else if (!strcmp(input, "POLYMER_STRESS13_6"))
    var = POLYMER_STRESS13_6;
  else if (!strcmp(input, "POLYMER_STRESS23_6"))
    var = POLYMER_STRESS23_6;
  else if (!strcmp(input, "POLYMER_STRESS33_6"))
    var = POLYMER_STRESS33_6;

  else if (!strcmp(input, "POLYMER_STRESS11_7"))
    var = POLYMER_STRESS11_7;
  else if (!strcmp(input, "POLYMER_STRESS12_7"))
    var = POLYMER_STRESS12_7;
  else if (!strcmp(input, "POLYMER_STRESS22_7"))
    var = POLYMER_STRESS22_7;
  else if (!strcmp(input, "POLYMER_STRESS13_7"))
    var = POLYMER_STRESS13_7;
  else if (!strcmp(input, "POLYMER_STRESS23_7"))
    var = POLYMER_STRESS23_7;
  else if (!strcmp(input, "POLYMER_STRESS33_7"))
    var = POLYMER_STRESS33_7;

  /*
   * Kluge to break up large if block. Problems with HP compiler!
   */
  if (var != -10)
    return var;
  /*
   * Special Species initialization keywords
   */
  else if (!strcmp(input, "SPECIES_MASS_FRACTION"))
    var = SPECIES_MASS_FRACTION;
  else if (!strcmp(input, "SPECIES_MOLE_FRACTION"))
    var = SPECIES_MOLE_FRACTION;
  else if (!strcmp(input, "SPECIES_VOL_FRACTION"))
    var = SPECIES_VOL_FRACTION;
  else if (!strcmp(input, "SPECIES_DENSITY"))
    var = SPECIES_DENSITY;
  else if (!strcmp(input, "SPECIES_CONCENTRATION"))
    var = SPECIES_CONCENTRATION;
  else if (!strcmp(input, "SPECIES_CAP_PRESSURE"))
    var = SPECIES_CAP_PRESSURE;
  else if (!strcmp(input, "SPECIES_UNDEFINED_FORM"))
    var = SPECIES_UNDEFINED_FORM;

  /*
   * Porous Media initialization keywords
   */
  else if (!strcmp(input, "POR_LIQ_PRES"))
    var = POR_LIQ_PRES;
  else if (!strcmp(input, "POR_GAS_PRES"))
    var = POR_GAS_PRES;
  else if (!strcmp(input, "POR_POROSITY"))
    var = POR_POROSITY;
  else if (!strcmp(input, "POR_TEMP"))
    var = POR_TEMP;
  else if (!strcmp(input, "POR_SATURATION"))
    var = POR_SATURATION;

  else if (!strcmp(input, "VORT_DIR1"))
    var = VORT_DIR1;
  else if (!strcmp(input, "VORT_DIR2"))
    var = VORT_DIR2;
  else if (!strcmp(input, "VORT_DIR3"))
    var = VORT_DIR3;
  else if (!strcmp(input, "VORT_LAMBDA"))
    var = VORT_LAMBDA;

  else if (!strcmp(input, "CURVATURE"))
    var = CURVATURE;

  else if (!strcmp(input, "BOND_EVOLUTION"))
    var = BOND_EVOLUTION;

  else if (!strcmp(input, "SURF_CHARGE"))
    var = SURF_CHARGE;

  else if (!strcmp(input, "EXT_VELOCITY"))
    var = EXT_VELOCITY;

  else if (!strcmp(input, "EFIELD1"))
    var = EFIELD1;
  else if (!strcmp(input, "EFIELD2"))
    var = EFIELD2;
  else if (!strcmp(input, "EFIELD3"))
    var = EFIELD3;

  else if (!strcmp(input, "ENORM"))
    var = ENORM;

  else if (!strcmp(input, "NORMAL1"))
    var = NORMAL1;
  else if (!strcmp(input, "NORMAL2"))
    var = NORMAL2;
  else if (!strcmp(input, "NORMAL3"))
    var = NORMAL3;

  else if (!strcmp(input, "SHELL_CURVATURE"))
    var = SHELL_CURVATURE;
  else if (!strcmp(input, "SHELL_CURVATURE2"))
    var = SHELL_CURVATURE2;
  else if (!strcmp(input, "SHELL_TENSION"))
    var = SHELL_TENSION;
  else if (!strcmp(input, "SHELL_X"))
    var = SHELL_X;
  else if (!strcmp(input, "SHELL_Y"))
    var = SHELL_Y;
  else if (!strcmp(input, "SHELL_USER"))
    var = SHELL_USER;
  else if (!strcmp(input, "PHASE1"))
    var = PHASE1;
  else if (!strcmp(input, "PHASE2"))
    var = PHASE2;
  else if (!strcmp(input, "PHASE3"))
    var = PHASE3;
  else if (!strcmp(input, "PHASE4"))
    var = PHASE4;
  else if (!strcmp(input, "PHASE5"))
    var = PHASE5;
  else if (!strcmp(input, "SHELL_ANGLE1"))
    var = SHELL_ANGLE1;
  else if (!strcmp(input, "SHELL_ANGLE2"))
    var = SHELL_ANGLE2;

  else if (!strcmp(input, "SHELL_SURF_DIV_V"))
    var = SHELL_SURF_DIV_V;
  else if (!strcmp(input, "SHELL_SURF_CURV"))
    var = SHELL_SURF_CURV;
  else if (!strcmp(input, "N_DOT_CURL_V"))
    var = N_DOT_CURL_V;
  else if (!strcmp(input, "GRAD_S_V_DOT_N1"))
    var = GRAD_S_V_DOT_N1;
  else if (!strcmp(input, "GRAD_S_V_DOT_N2"))
    var = GRAD_S_V_DOT_N2;
  else if (!strcmp(input, "GRAD_S_V_DOT_N3"))
    var = GRAD_S_V_DOT_N3;
  else if (!strcmp(input, "ACOUS_PREAL"))
    var = ACOUS_PREAL;
  else if (!strcmp(input, "ACOUS_PIMAG"))
    var = ACOUS_PIMAG;
  else if (!strcmp(input, "ACOUS_REYN_STRESS"))
    var = ACOUS_REYN_STRESS;
  else if (!strcmp(input, "SHELL_BDYVELO"))
    var = SHELL_BDYVELO;
  else if (!strcmp(input, "SHELL_LUBP"))
    var = SHELL_LUBP;
  else if (!strcmp(input, "LUBP"))
    var = LUBP;
  else if (!strcmp(input, "LUBP_2"))
    var = LUBP_2;
  else if (!strcmp(input, "SHELL_SAT_CLOSED"))
    var = SHELL_SAT_CLOSED;
  else if (!strcmp(input, "SHELL_PRESS_OPEN"))
    var = SHELL_PRESS_OPEN;
  else if (!strcmp(input, "SHELL_PRESS_OPEN_2"))
    var = SHELL_PRESS_OPEN_2;
  else if (!strcmp(input, "POR_SINK_MASS"))
    var = POR_SINK_MASS;
  else if (!strcmp(input, "SHELL_FILMP"))
    var = SHELL_FILMP;
  else if (!strcmp(input, "SHELL_FILMH"))
    var = SHELL_FILMH;
  else if (!strcmp(input, "SHELL_PARTC"))
    var = SHELL_PARTC;
  else if (!strcmp(input, "SHELL_TEMPERATURE"))
    var = SHELL_TEMPERATURE;
  else if (!strcmp(input, "SHELL_DELTAH"))
    var = SHELL_DELTAH;
  else if (!strcmp(input, "SHELL_LUB_CURV"))
    var = SHELL_LUB_CURV;
  else if (!strcmp(input, "SHELL_LUB_CURV_2"))
    var = SHELL_LUB_CURV_2;
  else if (!strcmp(input, "SHELL_SAT_GASN"))
    var = SHELL_SAT_GASN;
  else if (!strcmp(input, "SHELL_SHEAR_TOP"))
    var = SHELL_SHEAR_TOP;
  else if (!strcmp(input, "SHELL_SHEAR_BOT"))
    var = SHELL_SHEAR_BOT;
  else if (!strcmp(input, "SHELL_CROSS_SHEAR"))
    var = SHELL_CROSS_SHEAR;
  else if (!strcmp(input, "MAX_STRAIN"))
    var = MAX_STRAIN;
  else if (!strcmp(input, "CUR_STRAIN"))
    var = CUR_STRAIN;
  else if (!strcmp(input, "LIGHT_INTP"))
    var = LIGHT_INTP;
  else if (!strcmp(input, "LIGHT_INTM"))
    var = LIGHT_INTM;
  else if (!strcmp(input, "LIGHT_INTD"))
    var = LIGHT_INTD;
  else if (!strcmp(input, "MOMENT0"))
    var = MOMENT0;
  else if (!strcmp(input, "MOMENT1"))
    var = MOMENT1;
  else if (!strcmp(input, "MOMENT2"))
    var = MOMENT2;
  else if (!strcmp(input, "MOMENT3"))
    var = MOMENT3;
  else if (!strcmp(input, "DENSITY_EQN"))
    var = DENSITY_EQN;

  else if (!strcmp(input, "TFMP_PRES"))
    var = TFMP_PRES;
  else if (!strcmp(input, "TFMP_SAT"))
    var = TFMP_SAT;
  else if (!strcmp(input, "RESTIME"))
    var = RESTIME;
  else if (!strcmp(input, "SHELL_SAT_1"))
    var = SHELL_SAT_1;
  else if (!strcmp(input, "SHELL_SAT_2"))
    var = SHELL_SAT_2;
  else if (!strcmp(input, "SHELL_SAT_3"))
    var = SHELL_SAT_3;
  else if (!strcmp(input, "EM_E1_REAL"))
    var = EM_E1_REAL;
  else if (!strcmp(input, "EM_E2_REAL"))
    var = EM_E2_REAL;
  else if (!strcmp(input, "EM_E3_REAL"))
    var = EM_E3_REAL;
  else if (!strcmp(input, "EM_E1_IMAG"))
    var = EM_E1_IMAG;
  else if (!strcmp(input, "EM_E2_IMAG"))
    var = EM_E2_IMAG;
  else if (!strcmp(input, "EM_E3_IMAG"))
    var = EM_E3_IMAG;
  else if (!strcmp(input, "EM_H1_REAL"))
    var = EM_H1_REAL;
  else if (!strcmp(input, "EM_H2_REAL"))
    var = EM_H2_REAL;
  else if (!strcmp(input, "EM_H3_REAL"))
    var = EM_H3_REAL;
  else if (!strcmp(input, "EM_H1_IMAG"))
    var = EM_H1_IMAG;
  else if (!strcmp(input, "EM_H2_IMAG"))
    var = EM_H2_IMAG;
  else if (!strcmp(input, "EM_H3_IMAG"))
    var = EM_H3_IMAG;
  else if (!strcmp(input, "EM_CONT_REAL"))
    var = EM_CONT_REAL;
  else if (!strcmp(input, "EM_CONT_IMAG"))
    var = EM_CONT_IMAG;
  else if (!strcmp(input, "USTAR"))
    var = USTAR;
  else if (!strcmp(input, "VSTAR"))
    var = VSTAR;
  else if (!strcmp(input, "WSTAR"))
    var = WSTAR;
  else if (!strcmp(input, "EDDY_NU"))
    var = EDDY_NU;
  else if (!strcmp(input, "TURB_K"))
    var = TURB_K;
  else if (!strcmp(input, "TURB_OMEGA"))
    var = TURB_OMEGA;

  /*
   * Kluge to break up large if block. Problems with HP compiler!
   */
  if (var != -10)
    return var;
  /*
   * Catch all for external field that is not part of the solution vector
   */
  if (!strcmp(input, "EXTERNAL"))
    var = EXTERNAL;

  else if (!strcmp(input, "MESH_POSITION1"))
    var = MESH_POSITION1;
  else if (!strcmp(input, "MESH_POSITION2"))
    var = MESH_POSITION2;
  else if (!strcmp(input, "MESH_POSITION3"))
    var = MESH_POSITION3;
  else if (!strcmp(input, "R_MOM_NORMAL"))
    var = R_MOM_NORMAL;
  else if (!strcmp(input, "R_MOM_TANG1"))
    var = R_MOM_TANG1;
  else if (!strcmp(input, "R_MOM_TANG2"))
    var = R_MOM_TANG2;
  else if (!strcmp(input, "R_MESH_NORMAL"))
    var = R_MESH_NORMAL;
  else if (!strcmp(input, "R_MESH_TANG1"))
    var = R_MESH_TANG1;
  else if (!strcmp(input, "R_MESH_TANG2"))
    var = R_MESH_TANG2;
  else if (!strcmp(input, "VEL_NORM"))
    var = VEL_NORM;
  else if (!strcmp(input, "VEL_NORMAL"))
    var = VEL_NORM;
  else if (!strcmp(input, "VEL_TANG1"))
    var = VEL_TANG1;
  else if (!strcmp(input, "VEL_TANG2"))
    var = VEL_TANG2;
  else if (!strcmp(input, "MESH_NORM"))
    var = MESH_NORM;
  else if (!strcmp(input, "MESH_NORMAL"))
    var = MESH_NORM;
  else if (!strcmp(input, "MESH_TANG1"))
    var = MESH_TANG1;
  else if (!strcmp(input, "MESH_TANG2"))
    var = MESH_TANG2;
  else if (!strcmp(input, "SOLID_POSITION1"))
    var = SOLID_POSITION1;
  else if (!strcmp(input, "SOLID_POSITION2"))
    var = SOLID_POSITION2;
  else if (!strcmp(input, "SOLID_POSITION3"))
    var = SOLID_POSITION3;
  else if (!strcmp(input, "R_SOLID_NORMAL"))
    var = R_SOLID_NORMAL;
  else if (!strcmp(input, "R_SOLID_TANG1"))
    var = R_SOLID_TANG1;
  else if (!strcmp(input, "R_SOLID_TANG2"))
    var = R_SOLID_TANG2;

  else if (!strcmp(input, "R_ANYTHING"))
    var = R_ANYTHING;

  else if (!strcmp(input, "D_VEL1_DT"))
    var = D_VEL1_DT;
  else if (!strcmp(input, "D_VEL2_DT"))
    var = D_VEL2_DT;
  else if (!strcmp(input, "D_VEL3_DT"))
    var = D_VEL3_DT;
  else if (!strcmp(input, "D_T_DT"))
    var = D_T_DT;
  else if (!strcmp(input, "D_Y_DT"))
    var = D_Y_DT;
  else if (!strcmp(input, "D_X1_DT"))
    var = D_X1_DT;
  else if (!strcmp(input, "D_X2_DT"))
    var = D_X2_DT;
  else if (!strcmp(input, "D_X3_DT"))
    var = D_X3_DT;
  else if (!strcmp(input, "D_S_DT"))
    var = D_S_DT;
  else if (!strcmp(input, "D_P_DT"))
    var = D_P_DT;
  else if (!strcmp(input, "SHELL_NORMAL3"))
    var = SHELL_NORMAL3;

  else {
    if (err_string != NULL) {
      fprintf(stderr, "variable_string_to_int error: %s\n", err_string);
    }
    GOMA_EH(GOMA_ERROR, "Invalid variable type on card expecting string names of variables");
    var = -2;
  }
  return var;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int species_type_str_to_int(char *input)

/*************************************************************************
 *
 * species_type_str_to_int
 *
 *  This is a small utility routime to do the string to integer conversion
 *  for species var type input lines. Fatal crashes are induced if
 *  string matches are not achieved.
 ************************************************************************/
{
  int var = 0;
  if (!strcmp(input, "SPECIES_MASS_FRACTION"))
    var = SPECIES_MASS_FRACTION;
  else if (!strcmp(input, "SPECIES_MOLE_FRACTION"))
    var = SPECIES_MOLE_FRACTION;
  else if (!strcmp(input, "SPECIES_VOL_FRACTION"))
    var = SPECIES_VOL_FRACTION;
  else if (!strcmp(input, "SPECIES_DENSITY"))
    var = SPECIES_DENSITY;
  else if (!strcmp(input, "SPECIES_CONCENTRATION"))
    var = SPECIES_CONCENTRATION;
  else if (!strcmp(input, "SPECIES_CAP_PRESSURE"))
    var = SPECIES_CAP_PRESSURE;
  else if (!strcmp(input, "SPECIES_UNDEFINED_FORM"))
    var = SPECIES_UNDEFINED_FORM;
  if (var == 0) {
    printf("Unknown species type string: %s\n", input);
    GOMA_EH(GOMA_ERROR, "Invalid species type input");
  }
  return var;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void species_type_int_to_str(char *str, const int var)

/*************************************************************************
 *
 *  species_type_int_to_str:
 *
 *  This is a small utility routine to do the int to string conversion
 *  for species var type.
 *
 *  Input
 * --------
 * var -> Species Type variable
 *
 *  Output
 * --------
 * str -> String containing the description of the species variable type
 ************************************************************************/
{
  if (var == SPECIES_MASS_FRACTION) {
    (void)strcpy(str, "SPECIES_MASS_FRACTION");
  } else if (var == SPECIES_MOLE_FRACTION) {
    (void)strcpy(str, "SPECIES_MOLE_FRACTION");
  } else if (var == SPECIES_VOL_FRACTION) {
    (void)strcpy(str, "SPECIES_VOL_FRACTION");
  } else if (var == SPECIES_DENSITY) {
    (void)strcpy(str, "SPECIES_DENSITY");
  } else if (var == SPECIES_CONCENTRATION) {
    (void)strcpy(str, "SPECIES_CONCENTRATION");
  } else if (var == SPECIES_CAP_PRESSURE) {
    (void)strcpy(str, "SPECIES_CAP_PRESSURE");
  } else if (var == SPECIES_UNDEFINED_FORM) {
    (void)strcpy(str, "SPECIES_UNDEFINED_FORM");
  } else if (var == 0) {
    (void)strcpy(str, "SPECIES_UNSPECIFIED");
  } else {
    printf("WARNING: unknown species type!: %d\n", var);
    (void)strcpy(str, "?????????????");
  }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int interpret_string(const char *string, char *retn_string)

/*************************************************************************
 *
 * interpret_string:
 *
 *     Attempts to translate a string into a single integer. If it
 *     succeeds then the routine returns TRUE. If it fails then the
 *     routine returns FALSE.
 *
 *  Input
 * ---------
 *   string = String to be converted
 *
 *  Output
 * ---------
 *  retn_string  = string on return
 *  return       = TRUE on successful conversion, and false on an
 *                 unsuccessful conversion.
 ************************************************************************/
{
  int retn;
  if ((retn = sscanf(string, "%s", retn_string)) != 1) {
    return (FALSE);
  }
  return TRUE;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int interpret_int(const char *string, int *retn_value)

/*************************************************************************
 *
 * interpret_int:
 *
 *     Attempts to translate a string into a single integer. If it
 *     succeeds then the routine returns TRUE. If it fails then the
 *     routine returns FALSE.
 *
 *  Input
 * ---------
 *   string = String to be converged
 *
 *  Output
 * ---------
 *  retn_value   = int value on return
 *  return       = TRUE on successful conversion, and false on an
 *                 unsuccessful conversion.
 ************************************************************************/
{
  int retn;
  if ((retn = sscanf(string, "%d", retn_value)) != 1) {
    *retn_value = retn;
    return (FALSE);
  }
  return TRUE;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int interpret_double(const char *string, double *retn_value)

/*************************************************************************
 *
 * interpret_double:
 *
 *     Attempts to translate a string into a single double. If it
 *     succeeds then the routine returns TRUE. If it fails then the
 *     routine returns FALSE.
 *
 *  Input
 * ---------
 *   string = String to be converged
 *
 *  Output
 * ---------
 *  retn_value   = double value on return
 *  return       = TRUE on successful conversion, and false on an
 *                 unsuccessful conversion.
 ************************************************************************/
{
  int retn;
  if ((retn = sscanf(string, "%le", retn_value)) != 1) {
    *retn_value = retn;
    return (FALSE);
  }
  return TRUE;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int indentify_species_ID_string(const char *input, const char *list, const int numList)

/*************************************************************************
 *
 *  indentify_species_ID_string:
 *
 *  This utility tries to interpret a character string to identify
 *  a species ID. First, it attempts to treat it as an integer.
 *  Then, it tries to match against a list of strings.
 *
 *  Note the species indecises start with a 0, as per C indexing
 *  conventions.
 *
 *  Input
 * --------
 *  input -> String to be interpreted as a species index.
 *
 *  Return
 * --------
 *  Successful matches return the species index.
 *  Unsuccessful matches return the value -1.
 ************************************************************************/
{
  int species_index = -1;
  int nTokes;
  TOKEN_STRUCT tok;
  nTokes = fillTokStruct(&tok, input);
  if (nTokes != 1) {
    fprintf(stderr, "indentify_species_ID_string ERROR: was expecting one token: %s\n", input);
    GOMA_EH(GOMA_ERROR, "interface error");
  }
  if (!interpret_int(tok.tok_ptr[0], &species_index)) {
    species_index = in_char_list(tok.tok_ptr[0], list, numList);
  }
  return species_index;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

int look_for_species_prop(FILE *imp,
                          const char *search_string,
                          MATRL_PROP_STRUCT *mat_ptr,
                          int *materialModel,
                          dbl material_property[],
                          dbl **User_constants,
                          int *User_count,
                          char *model_name,
                          const int num_values_expected,
                          int *species_index,
                          char *echo_string)

/*************************************************************************
 *
 * look_for_species_prop -- generic code to read a species property
 *
 *
 *       Reads the string:
 *
 *  "search_string" = "model_name"  species_index   value
 *  "search_string" = "model_name"  species_index   value1 value2 value3
 *
 *  Alternatively, species_index can be replaced with the current
 *  string definition of the species name.
 *
 *    It understands the following generic model types:
 *        CONSTANT
 *        USER
 *        USER_GEN
 *
 *  "search_string" =  CONSTANT  species_index   value
 *  "search_string" =  USER      species_index   value1 value2 value3 ...
 *  "search_string" =  USER_GEN  species_index   value1 value2 value3 ...
 *
 *  Note:
 *     This routine always finishes reading the search_string line.
 *  All parameters are storred in the User_constants array.
 *
 *  Input
 * -------
 *
 *  imp             - ptr to input stream (in)
 * search_string,   -  search string (in)
 * mat_ptr          - Pointer to the current material structure
 * num_values_expected - number of values expected. If this is positive
 *                   then checking of the number of doubles read
 *                   on the command line is carried out. If zero,
 *                   this argument is ignored.
 *
 *  Output
 * -------
 * materialModel    - material model found. If the material
 *                    model is unknown, then this value is returned
 *                    as zero.
 * material_property - This is a vector over the species in the
 *                  material. If the materialModel is a scalar
 *                  then the one scalar is storred in the
 *                  location, material_property[species_ID],
 *                  on return to the calling program. Otherwise
 *                  this vector is not used.
 * User_constants, - ptr to vector of double constants for input
 *                   into the model dependent material property.
 *                   It is also used for user defined models.
 * User_count[]     - Length of the vector of double constants.
 *                    This is a vector over species index.
 * model_name      - The name of the model read
 * species_index   - On output this is equal to the species ID
 *
 *  Return Value
 * --------------
 *    -1 : This number is returned if the search string, "search_string"
 *         is not found in the material property file.
 *         This number is also returned if the model name is not one of the
 *         special cases it knows how to parse.
 *     1 : Successful completion of the routine using a model which
 *         this program knows about
 *     0 : Successful completion of the routine using a model which
 *         this program doesn't know about. It may need further
 *         checking.
 ***********************************************************************/
{
  static char yo[] = "look_for_species_prop";
  char input[MAX_CHAR_IN_INPUT]; /* storage for input strings */
  int iread = -1;                /* status flag  */
  int species_ID, fill_user_constant, i, retn, num_args, numTok;
  char *species_str;
  char err_mesg[MAX_CHAR_ERR_MSG];
  double *arg_list = NULL;
  TOKEN_STRUCT tok;

  input[0] = '\0';

  /*
   *  Search forward in the file for the search_string.
   *  Stop after the equals sign, if the search string is found
   */
  retn = look_forward_optional(imp, search_string, input, '=');

  /*
   *  If our search failed, then we set the model_name to a blank, the species
   *  index to -1 and we return -1 as the outcome. The materialModel is not
   *  touched, because it may contain valid information
   */
  if (retn == -1) {
    (void)sprintf(model_name, " ");
    *species_index = -1;
    echo_string[0] = '\0';
    return -1;
  }

  /*
   * Read the rest of the line and tokenize it
   */
  retn = read_string(imp, input, '\n');
  numTok = fillTokStruct(&tok, input);
  if (numTok < 2) {
    (void)fprintf(stderr, "%s error: reading model name string, mat file \"%s\", property %s", yo,
                  mat_ptr->Material_Name, search_string);
    (void)sprintf(err_mesg, "Number of tokens isn't sufficient: %s\n", input);
    GOMA_EH(GOMA_ERROR, err_mesg);
  }

  /*
   * Interpret the Model name as the first token after the equals sign
   */
  (void)strcpy(model_name, tok.tok_ptr[0]);
  /*
   * Interpret the species name as the second token after the equals sign.
   * Then interpret it as a species index.
   */
  species_str = tok.tok_ptr[1];
  species_ID =
      indentify_species_ID_string(species_str, mat_ptr->Species_Names[0], mat_ptr->Num_Species);
  if (species_ID == -1) {
    fprintf(stderr, "%s error: reading model name string, mat file \"%s\", property %s\n", yo,
            mat_ptr->Material_Name, search_string);
    sprintf(err_mesg, "\tUnsuccessful match of species string: %s\n", species_str);
    GOMA_EH(GOMA_ERROR, err_mesg);
  }
  if (species_ID >= mat_ptr->Num_Species || species_ID < 0) {
    fprintf(stderr, "%s error: reading model name string, mat file \"%s\", property %s\n", yo,
            mat_ptr->Material_Name, search_string);
    sprintf(err_mesg, "\tspecies index is out of bounds: %d\n", species_ID);
    GOMA_EH(GOMA_ERROR, err_mesg);
  }

  /*
   *  Check on the number of arguments
   */
  num_args = numTok - 2;
  if (num_values_expected > 0 && strcmp(model_name, "USER")) {
    if (num_args != num_values_expected) {
      fprintf(stderr, "%s Warning: reading model name string, mat file \"%s\", property %s\n", yo,
              mat_ptr->Material_Name, search_string);

      fprintf(stderr, "\tNumber of expected args, %d, differs from actual, %d\n",
              num_values_expected, num_args);
    }
  }

  SPF(echo_string, "%s = %s %s", search_string, model_name, species_str);

  /*
   * Set branching flags depending upont the Model Name read on the line
   *
   */

  fill_user_constant = TRUE;
  iread = 1;
  if (!strcmp(model_name, "CONSTANT")) {
    if (num_values_expected == 1) {
      fill_user_constant = FALSE;
    }
    materialModel[species_ID] = CONSTANT;
  } else if (!strcmp(model_name, "USER")) {
    materialModel[species_ID] = USER;
  } else if (!strcmp(model_name, "USER_GEN")) {
    materialModel[species_ID] = USER_GEN;
  } else {
    materialModel[species_ID] = 0;
    iread = 0;
  }

  /*
   * Read the rest of the line into the buffer, input. Then, tokenize it.
   */

  arg_list = alloc_dbl_1(num_args, DBL_NOINIT);
  for (i = 0; i < num_args; i++) {
    retn = interpret_double(tok.tok_ptr[i + 2], arg_list + i);
    if (!retn) {
      fprintf(stderr, "%s error: reading model name string, mat file \"%s\", property %s\n", yo,
              mat_ptr->Material_Name, search_string);
      sprintf(err_mesg, "\tUnsuccessful interpretation of double: %s\n", tok.tok_ptr[i]);
      GOMA_EH(GOMA_ERROR, err_mesg);
    }
  }

  SPF_DBL_VEC(endofstring(echo_string), num_args, arg_list);

  /*
   *  Decide whether to fill in material_property[] location
   */
  if (!strcmp(model_name, "CONSTANT") || num_args <= 1) {
    material_property[species_ID] = arg_list[0];
    fill_user_constant = FALSE;
    if (!strcmp(model_name, "PHOTO_CURING"))
      materialModel[species_ID] = PHOTO_CURING;
    if (!strcmp(model_name, "DROP_EVAP"))
      materialModel[species_ID] = DROP_EVAP;
  }

  /*
   * Transfer the argument list of doubles to the return address if needed.
   * If not needed, we need to free the space
   */
  if (fill_user_constant) {
    if (User_constants == NULL) {
      if (num_args <= 1) {
        safer_free((void **)&(arg_list));
      } else {
        fprintf(stderr, "%s error: reading model name string, mat file \"%s\", property %s\n", yo,
                mat_ptr->Material_Name, search_string);
        sprintf(err_mesg, "/tSpace for pointer vector over species needs to be malloced first\n");
        GOMA_EH(GOMA_ERROR, err_mesg);
      }
    } else {
      safer_free((void **)&(User_constants[species_ID]));
      User_constants[species_ID] = (dbl *)array_alloc(1, num_args, sizeof(dbl));
      User_constants[species_ID] = arg_list;
      /*    for(i=0 ; i<num_args ; i++)	{
          User_constants[species_ID][i] = arg_list[i];
              }  */
      User_count[species_ID] = num_args;
    }
  } else {
    safer_free((void **)&(arg_list));
  }
  *species_index = species_ID;
  return iread;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

int look_for_porous_prop(FILE *imp,
                         const char *search_string,
                         MATRL_PROP_STRUCT *mat_ptr,
                         int *materialModel,
                         dbl material_property[],
                         dbl **User_constants,
                         int *User_count,
                         char *model_name,
                         const int num_values_expected,
                         int *porous_index,
                         char *echo_string)

/*************************************************************************
 *
 * look_for_porous_prop -- generic code to read a porous media
 *                         phase property
 *
 *
 *       Reads the string:
 *
 *  "search_string" = "model_name"  porous_index   value
 *  "search_string" = "model_name"  porous_index   value1 value2 value3
 *
 *  Alternatively, porous_index can be replaced with the current
 *  string definition of the porous phase name.
 *
 *    It understands the following generic model types:
 *        CONSTANT
 *        USER
 *        USER_GEN
 *
 *  "search_string" =  CONSTANT  porous_index   value
 *  "search_string" =  USER      porous_index   value1 value2 value3 ...
 *  "search_string" =  USER_GEN  porous_index   value1 value2 value3 ...
 *
 *  Note:
 *     This routine always finishes reading the search_string line.
 *  All parameters are stored in the User_constants array.
 *
 *  Input
 * -------
 *
 *  imp             - ptr to input stream (in)
 * search_string,   -  search string (in)
 * mat_ptr          - Pointer to the current material structure
 * num_values_expected - number of values expected. If this is positive
 *                   then checking of the number of doubles read
 *                   on the command line is carried out. If zero,
 *                   this argument is ignored.
 *
 *  Output
 * -------
 * materialModel    - material model found. If the material
 *                    model is unknown, then this value is returned
 *                    as zero.
 * material_property - This is a vector over the porous phases in the
 *                  material. If the material Model is a scalar
 *                  then the one scalar is storred in the
 *                  location, material_property[porous_phase_ID],
 *                  on return to the calling program. Otherwise
 *                  this vector is not used.
 * User_constants, - ptr to vector of double constants for input
 *                   into the model dependent material property.
 *                   It is also used for user defined models.
 * User_count[]     - Length of the vector of double constants.
 *                    This is a vector over porous phase index.
 * model_name      - The name of the model read
 * porous_index   - On output this is equal to the porous phase ID
 *
 *  Return Value
 * --------------
 *    -1 : This number is returned if the search string, "search_string"
 *         is not found in the material property file.
 *         This number is also returned if the model name is not one of the
 *         special cases it knows how to parse.
 *     1 : Successful completion of the routine using a model which
 *         this program knows about
 *     0 : Successful completion of the routine using a model which
 *         this program doesn't know about. It may need further
 *         checking.
 ***********************************************************************/
{
  static char yo[] = "look_for_porous_prop";
  char input[MAX_CHAR_IN_INPUT]; /* storage for input strings */
  int iread = -1;                /* status flag  */
  int porous_ID, fill_user_constant, i, retn, num_args, numTok;
  char *porous_str;
  char err_mesg[MAX_CHAR_ERR_MSG];
  double *arg_list = NULL;
  TOKEN_STRUCT tok;

  input[0] = '\0';
  model_name[0] = '\0';

  /*
   *  Search forward in the file for the search_string.
   *  Stop after the equals sign, if the search string is found
   */
  retn = look_forward_optional(imp, search_string, input, '=');

  /*
   *  If our search failed, then we set the model_name to a blank, the porous
   *  phase index to -1 and we return -1 as the outcome. The material Model
   *  is not touched, because it may contain valid information
   */
  if (retn == -1) {
    (void)sprintf(model_name, " ");
    *porous_index = -1;
    return -1;
  }

  /*
   * Read the rest of the line and tokenize it
   */

  retn = read_string(imp, input, '\n');
  numTok = fillTokStruct(&tok, input);
  if (numTok < 2) {
    (void)fprintf(stderr, "%s error: reading model name string, mat file \"%s\", property %s", yo,
                  mat_ptr->Material_Name, search_string);
    (void)sprintf(err_mesg, "Number of tokens isn't sufficient: %s\n", input);
    GOMA_EH(GOMA_ERROR, err_mesg);
  }

  /*
   * Interpret the Model name as the first token after the equals sign
   */
  (void)strcpy(model_name, tok.tok_ptr[0]);

  /*
   * Interpret the species name as the second token after the equals sign.
   * Then interpret it as a porous phase index.
   */
  porous_str = tok.tok_ptr[1];
  porous_ID =
      indentify_species_ID_string(porous_str, mat_ptr->Porous_Names[0], mat_ptr->Num_Porous_Eqn);
  if (porous_ID == -1) {
    fprintf(stderr, "%s error: reading model name string, mat file \"%s\", property %s\n", yo,
            mat_ptr->Material_Name, search_string);
    sprintf(err_mesg, "\tUnsuccessful match of porous string: %s\n", porous_str);
    GOMA_EH(GOMA_ERROR, err_mesg);
  }

  if (porous_ID >= mat_ptr->Num_Porous_Eqn || porous_ID < 0) {
    fprintf(stderr, "%s error: reading model name string, mat file \"%s\", property %s\n", yo,
            mat_ptr->Material_Name, search_string);
    sprintf(err_mesg, "\tporous phase index is out of bounds: %d\n", porous_ID);
    GOMA_EH(GOMA_ERROR, err_mesg);
  }

  /*
   *  Check on the number of arguments
   */
  num_args = numTok - 2;
  if (num_values_expected > 0) {
    if (num_args != num_values_expected) {
      fprintf(stderr, "%s Warning: reading model name string, mat file \"%s\", property %s\n", yo,
              mat_ptr->Material_Name, search_string);

      fprintf(stderr, "\tNumber of expected args, %d, differs from actual, %d\n",
              num_values_expected, num_args);
    }
  }

  SPF(echo_string, "%s = %s %s", search_string, model_name, porous_str);

  /*
   * Set branching flags depending upont the Model Name read on the line
   *
   */
  fill_user_constant = TRUE;
  iread = 1;
  if (!strcmp(model_name, "CONSTANT")) {
    if (num_values_expected == 1) {
      fill_user_constant = FALSE;
    }
    materialModel[porous_ID] = CONSTANT;
  } else if (!strcmp(model_name, "USER")) {
    materialModel[porous_ID] = USER;
  } else if (!strcmp(model_name, "USER_GEN")) {
    materialModel[porous_ID] = USER_GEN;
  } else {
    materialModel[porous_ID] = 0;
    iread = 0;
  }

  /*
   * Read the rest of the line into the buffer, input. Then, tokenize it.
   */

  arg_list = alloc_dbl_1(num_args, DBL_NOINIT);
  for (i = 0; i < num_args; i++) {
    retn = interpret_double(tok.tok_ptr[i + 2], arg_list + i);
    if (!retn) {
      fprintf(stderr, "%s error: reading model name string, mat file \"%s\", property %s\n", yo,
              mat_ptr->Material_Name, search_string);
      sprintf(err_mesg, "\tUnsuccessful interpretation of double: %s\n", tok.tok_ptr[i]);
      GOMA_EH(GOMA_ERROR, err_mesg);
    }
  }

  SPF_DBL_VEC(endofstring(echo_string), num_args, arg_list);

  /*
   *  Decide whether to fill in material_property[] location
   */
  if (!strcmp(model_name, "CONSTANT")) {
    material_property[porous_ID] = arg_list[0];
  }

  /*
   * Transfer the argument list of doubles to the return address if needed.
   * If not needed, we need to free the space
   */
  if (fill_user_constant) {
    if (User_constants == NULL) {
      fprintf(stderr, "%s error: reading model name string, mat file \"%s\", property %s\n", yo,
              mat_ptr->Material_Name, search_string);
      sprintf(err_mesg,
              "/tSpace for pointer vector over porous phases needs to be malloced first\n");
      GOMA_EH(GOMA_ERROR, err_mesg);
    }
    if (User_constants[porous_ID] != NULL) {
      safer_free((void **)&(User_constants[porous_ID]));
    }
    User_constants[porous_ID] = arg_list;
    User_count[porous_ID] = num_args;
  } else {
    safer_free((void **)&(arg_list));
  }

  *porous_index = porous_ID;
  return iread;
} /* end of function   look_for_porous_prop   */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int look_for_species_proptable(FILE *imp,
                               char *search_string,
                               MATRL_PROP_STRUCT *mat_ptr,
                               int *materialModel,
                               dbl material_property[],
                               dbl **User_constants,
                               int *User_count,
                               int *tableindex,
                               char *model_name,
                               const int num_values_expected,
                               int *species_index,
                               char *echo_string)

/*************************************************************************
 *
 * look_for_species_prop -- generic code to read a species property
 *
 *
 *       Reads the string:
 *
 *  "search_string" = "model_name"  species_index   value
 *  "search_string" = "model_name"  species_index   value1 value2 value3
 *
 *  Alternatively, species_index can be replaced with the current
 *  string definition of the species name.
 *
 *    It understands the following generic model types:
 *        CONSTANT
 *        USER
 *        USER_GEN
 *        TABLE
 *
 *  "search_string" =  CONSTANT  species_index   value
 *  "search_string" =  USER      species_index   value1 value2 value3 ...
 *  "search_string" =  USER_GEN  species_index   value1 value2 value3 ...
 *
 *  Note:
 *     This routine always finishes reading the search_string line.
 *  All parameters are storred in the User_constants array.
 *
 *  Input
 * -------
 *
 *  imp             - ptr to input stream (in)
 * search_string,   -  search string (in)
 * mat_ptr          - Pointer to the current material structure
 * num_values_expected - number of values expected. If this is positive
 *                   then checking of the number of doubles read
 *                   on the command line is carried out. If zero,
 *                   this argument is ignored.
 *
 *  Output
 * -------
 * materialModel    - material model found. If the material
 *                    model is unknown, then this value is returned
 *                    as zero.
 * material_property - This is a vector over the species in the
 *                  material. If the materialModel is a scalar
 *                  then the one scalar is storred in the
 *                  location, material_property[species_ID],
 *                  on return to the calling program. Otherwise
 *                  this vector is not used.
 * User_constants, - ptr to vector of double constants for input
 *                   into the model dependent material property.
 *                   It is also used for user defined models.
 * User_count[]     - Length of the vector of double constants.
 *                    This is a vector over species index.
 * model_name      - The name of the model read
 * species_index   - On output this is equal to the species ID
 *
 *  Return Value
 * --------------
 *    -1 : This number is returned if the search string, "search_string"
 *         is not found in the material property file.
 *         This number is also returned if the model name is not one of the
 *         special cases it knows how to parse.
 *     1 : Successful completion of the routine using a model which
 *         this program knows about
 *     0 : Successful completion of the routine using a model which
 *         this program doesn't know about. It may need further
 *         checking.
 ***********************************************************************/
{
  static char yo[] = "look_for_species_proptable";
  char input[MAX_CHAR_IN_INPUT]; /* storage for input strings */
  int iread = -1;                /* status flag  */
  int species_ID = 0, fill_user_constant, i, retn, num_args = 0, numTok;
  char *species_str;
  char err_mesg[MAX_CHAR_ERR_MSG];
  double *arg_list = NULL;
  TOKEN_STRUCT tok;
  fpos_t file_position;
  struct Data_Table *table_local;

  input[0] = '\0';

  /*
   *  Search forward in the file for the search_string.
   *  Stop after the equals sign, if the search string is found
   */
  retn = look_forward_optional(imp, search_string, input, '=');

  /*
   *  If our search failed, then we set the model_name to a blank, the species
   *  index to -1 and we return -1 as the outcome. The materialModel is not
   *  touched, because it may contain valid information
   */
  if (retn == -1) {
    (void)sprintf(model_name, " ");
    *species_index = -1;
    echo_string[0] = '\0';
    return -1;
  }

  /*
   * Read the rest of the line and tokenize it
   */
#ifndef tflop
  fgetpos(imp, &file_position);
#else
  file_position = ftell(imp);
#endif
  retn = read_string(imp, input, '\n');
  numTok = fillTokStruct(&tok, input);
  if (numTok < 2) {
    (void)fprintf(stderr, "%s error: reading model name string, material file \"%s\", property %s",
                  yo, mat_ptr->Material_Name, search_string);
    (void)sprintf(err_mesg, "Number of tokens isn't sufficient: %s\n", input);
    GOMA_EH(GOMA_ERROR, err_mesg);
  }
  /*
   * Interpret the Model name as the first token after the equals sign
   */
  (void)strcpy(model_name, tok.tok_ptr[0]);

  SPF(echo_string, "%s = %s", search_string, model_name);

  /* If model_name is TABLE */
  if (strcmp(model_name, "TABLE")) {

    /*
     * Interpret the species name as the second token after the equals sign.
     * Then interpret it as a species index.
     */
    species_str = tok.tok_ptr[1];
    species_ID =
        indentify_species_ID_string(species_str, mat_ptr->Species_Names[0], mat_ptr->Num_Species);
    if (species_ID == -1) {
      fprintf(stderr, "%s error: reading model name string, mat file \"%s\", property %s\n", yo,
              mat_ptr->Material_Name, search_string);
      sprintf(err_mesg, "\tUnsuccessful match of species string: %s\n", species_str);
      GOMA_EH(GOMA_ERROR, err_mesg);
    }
    if (species_ID >= mat_ptr->Num_Species || species_ID < 0) {
      fprintf(stderr, "%s error: reading model name string, mat file \"%s\", property %s\n", yo,
              mat_ptr->Material_Name, search_string);
      sprintf(err_mesg, "\tspecies index is out of bounds: %d\n", species_ID);
      GOMA_EH(GOMA_ERROR, err_mesg);
    }
    /*
     *  Check on the number of arguments
     */
    num_args = numTok - 2;
    if (num_values_expected > 0) {
      if (num_args != num_values_expected) {
        fprintf(stderr, "%s Warning: reading model name string, mat file \"%s\", property %s\n", yo,
                mat_ptr->Material_Name, search_string);

        fprintf(stderr, "\tNumber of expected args, %d, differs from actual, %d\n",
                num_values_expected, num_args);
      }
    }
    SPF(endofstring(echo_string), " %s", species_str);
  } /*  End TABLE if */

  /*
   * Set branching flags depending upon the Model Name read on the line
   *
   */
  fill_user_constant = TRUE;
  iread = 1;

  if (!strcmp(model_name, "CONSTANT")) {
    if (num_values_expected == 1) {
      fill_user_constant = FALSE;
    }
    materialModel[species_ID] = CONSTANT;
  } else if (!strcmp(model_name, "USER")) {
    materialModel[species_ID] = USER;
  } else if (!strcmp(model_name, "USER_GEN")) {
    materialModel[species_ID] = USER_GEN;
  } else if (!strcmp(model_name, "TABLE")) {
    /*
     * Fall through for all unimplemented cases
     */
    if (num_MP_Tables == MAX_MP_TABLES) {
      fprintf(stderr, " Maximum MP TABLEs Exceeded");
      exit(-1);
    }

    /*  Reset position to after = and reread model_name */
#ifndef tflop
    fsetpos(imp, &file_position);
#else
    fseek(imp, file_position, SEEK_SET);
#endif
    if (fscanf(imp, "%s", model_name) != 1) {
      fprintf(stderr,
              "ERROR ReReading model_name for TABLE model for property %s in mat \"%s\" requires "
              "space!",
              mat_ptr->Material_Name, search_string);
    }

    table_local = (struct Data_Table *)smalloc(sizeof(struct Data_Table));
    MP_Tables[num_MP_Tables] = setup_table_MP(imp, table_local, search_string);
    species_ID = MP_Tables[num_MP_Tables]->species_eq;
    *tableindex = num_MP_Tables++;
    materialModel[species_ID] = TABLE;

  } else {
    materialModel[species_ID] = 0;
    iread = 0;
  }

  /* Skip if TABLE model */
  if (strcmp(model_name, "TABLE") != 0) {

    /*
     * Read the rest of the line into the buffer, input. Then, tokenize it.
     */

    arg_list = alloc_dbl_1(num_args, DBL_NOINIT);
    for (i = 0; i < num_args; i++) {
      retn = interpret_double(tok.tok_ptr[i + 2], arg_list + i);
      if (!retn) {
        fprintf(stderr, "%s error: reading model name string, mat file \"%s\", property %s\n", yo,
                mat_ptr->Material_Name, search_string);
        sprintf(err_mesg, "\tUnsuccessful interpretation of double: %s\n", tok.tok_ptr[i]);
        GOMA_EH(GOMA_ERROR, err_mesg);
      }
    }

    SPF_DBL_VEC(endofstring(echo_string), num_args, arg_list);

    /*
     *  Decide whether to fill in material_property[] location
     */
    if (!strcmp(model_name, "CONSTANT")) {
      material_property[species_ID] = arg_list[0];
    }

    /*
     * Transfer the argument list of doubles to the return address if needed.
     * If not needed, we need to free the space
     */
    if (fill_user_constant) {
      if (User_constants == NULL) {
        fprintf(stderr, "%s error: reading model name string, mat file \"%s\", property %s\n", yo,
                mat_ptr->Material_Name, search_string);
        sprintf(err_mesg, "/tSpace for pointer vector over species needs to be malloced first\n");
        GOMA_EH(GOMA_ERROR, err_mesg);
      }
      if (User_constants[species_ID] != NULL) {
        safer_free((void **)&(User_constants[species_ID]));
      }
      User_constants[species_ID] = arg_list;
      User_count[species_ID] = num_args;
    } else {
      safer_free((void **)&(arg_list));
    }

  } /* END Skip for TABLE */
  *species_index = species_ID;

  return iread;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int stokenize(char *string, const char *delimiters, char **token_ptrs, const int max_tokens)

/************************************************************************
 *
 * stokenize:
 *
 *    This function will break up a string into its respective "tokens".
 *    It is a wrapper around the strtok() ISO C library routine.
 *
 *  input
 *  ----------
 *    string        - String to be tokenized.  Note, that the string is
 *                    changed by this procedure
 *    delimiters    - String containing a list of delimiters.
 *                      e.g., char *delimiters  = " \t\n";
 *    max_tokens    - Maximum number of tokens that can be found
 *
 *  output
 * -----------
 *    token_ptrs    - Vector of pointers to strings. They contain the
 *                    input string's tokens.
 *                     char *token_ptrs[max_tokens]
 *
 *  Return Value
 * ---------------
 *                  - Number of tokens actually found. If the string
 *                    consists totally of delimiters, 0 will be returned.
 *************************************************************************/
{
  register int i = 0;
  if ((token_ptrs[0] = strtok(string, delimiters)) == NULL) {
    return (0);
  } else {
    do {
      if ((++i) == max_tokens)
        break;
    } while ((token_ptrs[i] = strtok(NULL, delimiters)) != NULL);
  }
  return i;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int tokenize_by_whsp(char *string, char **token_ptrs, const int max_tokens)

/***************************************************************************
 *
 * tokenize_by_whsp:
 *
 *      This function will break up a string into its respective "tokens"
 *  using white space, the equals sign, and a comma as the definition
 *  of a delimiter.
 *
 *  input
 *  ----------
 *    string        - String to be tokenized.  Note, that the string is
 *                    changed by this procedure
 *    max_tokens    - Maximum number of tokens that can be found
 *
 *  output
 * -----------
 *    token_ptrs    - Vector of pointers to strings. They contain the input
 *                    string's tokens.
 *                     char *token_ptrs[max_tokens]
 *    retn_value    - Number of tokens actually found. If the string
 *                    consists totally of delimiters, 0 will be returned.
 ****************************************************************************/
{
  char *delimiters = " \t\f\r\v\n=,";
  return stokenize(string, delimiters, token_ptrs, max_tokens);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int fillTokStruct(TOKEN_STRUCT *tp, const char *string)

/************************************************************************
 *
 * fillTokStruct:
 *
 *    Tokenizes a string and stores the results in a Token structure.
 *    Leaves the original string unchanged.
 *
 *  Input
 * -------
 *   string = Input string to be tokenized
 *
 *  Output
 * --------
 *   tp    = On input, it must point to a previously allocated token
 *           structure. On output, it contains the results of the
 *           tokenization procedure.
 *
 *  Return
 * --------
 *        Returns the number of tokens found.
 *        -1 if there is an interface error.
 *************************************************************************/
{
  if (tp == NULL)
    return -1;
  if (string == NULL) {
    tp->orig_str[0] = tp->tok_str[0] = '\0';
    tp->ntokes = 0;
    tp->tok_ptr[0] = NULL;
  } else {
    (void)strcpy(tp->orig_str, string);
    (void)strcpy(tp->tok_str, string);
    tp->ntokes = tokenize_by_whsp(tp->tok_str, tp->tok_ptr, MAXTOKENS);
  }
  return tp->ntokes;
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int in_char_list(const char *target, const char *list, const int num_list)

/*************************************************************************
 *
 * in_char_list()
 *
 *   Finds the location of target in the list of strings, list. Returns
 *   the index number if it finds target in the list, and -1 if it does
 *   not.
 *
 *   String comparisons are case and white-space sensitive.
 *
 *  Input
 * -------
 *    target = string to be looked up
 *    list   = Vector of pointers to strings
 *    num_list = length of the string list
 *
 *  Return
 * --------
 *    Position of the string in the list. -1 if the string is not
 *    found or target is the null pointer.
 ************************************************************************/
{
  int i;
  char const *list_ptr = list;
  if (target == NULL)
    return -1;
  for (i = 0; i < num_list; i++) {
    if (list_ptr != NULL) {
      if (!strcmp(target, list_ptr))
        return i;
    }
    list_ptr++;
  }
  return -1;
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int read_line(FILE *ifp, char input[], const int print_flag)

/***************************************************************************
 *   read_line()
 *
 *      Reads a line of input.   The line is
 *      printed to standard output, if print_flag is true.
 *           The line is returned in the character
 *      string pointed to by input. Leading and trailing white spaces are
 *      stripped from the line.
 *          The number of characters, excluding the null character, is
 *      returned, except when read_string encounters an error condition
 *      (negative return values).  Then, the error condition is returned.
 ***************************************************************************/
{
  int retn_value;

  /*
   *   read the file up to the next new line, read_string will return
   *   a 0 or positive number for a success, and a negative value for failure
   */

  retn_value = read_string(ifp, input, '\n');

  /*
   *   Print out the line before stripping it of comments and white space
   */

  if (print_flag)
    (void)printf("%s\n", input);

  /*
   *   Strip the return line of comments and leading/trailing white space
   *   Use the function strip to return the number of characters remaining
   */

  if (retn_value >= 0)
    return (strip(input));

  /*
   *   If an error condition occurred in read_string, return the error
   *   condition value instead of the character count
   */

  (void)strip(input);
  return (retn_value);
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int read_string(FILE *ifp, char string[], const char ch)

/******************************************************************
 *
 *  This routine reads the standard input until encountering
 *  the end-of-file, a newline,  the character 'ch' or until
 *  MAX_CHAR_IN_INPUT characters are read. The inputted characters
 *  are read into 'string'.
 *   If an error occurs, -1 is returned and an error message
 *  is written to standard error.  Upon successful
 *  completion, the number of characters read plus 1 for the
 *  null character at the end of the string is returned.
 *  After the routine, the file pointer is positioned after
 *  the terminating character.
 *
 *    Author:			Ray S. Tuminaro Div 1422
 *    Date:			8/8/90
 *    revised:		10/2/90 John N. Shadid
 *
 *     Parameter list:
 *
 *    ifp    == pointer to file "input". On return, the file pointer
 *              is positioned after the last character read.
 *    string == On output 'string' contains the characters read
 *  	          from the input stream, terminated by a null
 *              character
 *    ch     == Termination character. That is, input function
 *      	  stops when 'ch' is read.
 *
 *    return == -1 if 1) EOF is encountered
 *                    2) MAX_CHAR_IN_INPUT characters are read before
 *                       an end of line or ch character is read.
 *           == strlen(string) + 1 if not an error condition
 *******************************************************************/
{
  int i = 0;
  int new_ch;
  while ((i < MAX_CHAR_IN_INPUT) && ((new_ch = getc(ifp)) != ch) && (new_ch != '\n') &&
         (new_ch != EOF)) {
    string[i++] = new_ch;
  }
  if (new_ch == EOF)
    return -1;
  if (i == MAX_CHAR_IN_INPUT) {
    fprintf(stderr, "read_string: %d chars w/ no \"%c\" found in\n%s", MAX_CHAR_IN_INPUT, ch,
            string);
    return (-1);
  }
  string[i] = '\0';
  return (i + 1);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int strip(char string[])

/*************************************************************************
 *
 * This routine strips off blanks and tabs (only leading and trailing
 * characters) in 'string'.
 *
 * Author:			Ray S. Tuminaro Div 1422
 * Date:			8/8/90
 *
 * Parameter list:
 *
 * string == On output 'string' contains the same characters as on
 *           input except the leading and trailing blanks and tabs
 *           have been removed.
 * return == New strlen of string[] upon return.
 ************************************************************************/
{
  int i, j;
  char ch;

  /* find first real character */
  i = 0;
  while (((ch = string[i]) != '\0') && ((ch == ' ') || (ch == '\t')))
    i++;

  /* move real part of string to the front */
  j = 0;
  while ((ch = string[j + i]) != '\0') {
    string[j] = string[j + i];
    j++;
  }
  string[j] = string[j + i];
  j--;

  /* remove trailing blanks */

  while ((j != -1) && (((ch = string[j]) == ' ') || (ch == '\t')))
    j--;
  string[j + 1] = '\0';
  return strlen(string);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
