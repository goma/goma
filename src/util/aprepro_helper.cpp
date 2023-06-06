#ifdef GOMA_ENABLE_APREPRO_LIB
#include <aprepro.h>
#include <fstream>
#include <sstream>
// Include MPI here for cases when openmpi is built with CXX support
// otherwise imports will be broken in the extern "C" section
#include <mpi.h>
extern "C" {
#define DISABLE_CPP
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "rf_bc.h"
#undef DISABLE_CPP
}
#include "util/aprepro_helper.h"

extern "C" goma_error aprepro_parse_file(char *infile, char *outfile) {
  SEAMS::Aprepro aprepro;

  auto result = aprepro.parse_file(infile);

  if (!result) {
    return GOMA_ERROR;
  }

  std::ofstream aprepro_out(outfile, std::ofstream::out);

  aprepro_out << aprepro.parsing_results().str();
  return GOMA_SUCCESS;
}

extern "C" goma_error aprepro_parse_goma_file(char *filename) {
  SEAMS::Aprepro aprepro;
  std::string stfilename = filename;

  auto result = aprepro.parse_file(stfilename);

  if (!result) {
    return GOMA_ERROR;
  }

  std::string outfile = "tmp." + stfilename;
  std::ofstream aprepro_out(outfile, std::ofstream::out);

  aprepro_out << aprepro.parsing_results().str();
  return GOMA_SUCCESS;
}

extern "C" goma_error aprepro_parse_goma_single_var(char *filename, char *var, double val) {
  SEAMS::Aprepro aprepro;
  std::string stfilename = filename;

  aprepro.add_variable(var, val, true);
  auto result = aprepro.parse_file(stfilename);

  if (!result) {
    return GOMA_ERROR;
  }

  std::string outfile = "tmp." + stfilename;
  std::ofstream aprepro_out(outfile, std::ofstream::out);

  aprepro_out << aprepro.parsing_results().str();
  return GOMA_SUCCESS;
}

extern "C" goma_error aprepro_parse_goma_augc(struct AC_Information *ac, double val) {
  SEAMS::Aprepro aprepro;
  std::string var = ac->AP_param;
  aprepro.add_variable(var, val, true);
  std::string filedata = ac->Aprepro_lib_string;
  auto result = aprepro.parse_string(filedata);

  if (!result) {
    return GOMA_ERROR;
  }

  auto aprepro_out = std::istringstream(aprepro.parsing_results().str());

  std::string line;
  while (std::getline(aprepro_out, line)) {
    if (line[0] == '$' || line[0] == '#') {
      continue;
    }
    std::stringstream ss(line);
    int bcid, bcidx;
    double computed_value;
    if (ss >> bcid >> bcidx >> computed_value) {
      switch (BC_Types[bcid].BC_Name) {
      case SPLINE_BC:
      case SPLINEX_BC:
      case SPLINEY_BC:
      case SPLINEZ_BC:
      case SPLINE_RS_BC:
      case SPLINEX_RS_BC:
      case SPLINEY_RS_BC:
      case SPLINEZ_RS_BC:
      case FILLET_BC:
      case DOUBLE_RAD_BC:
      case FEATURE_ROLLON_BC:
      case ROLL_FLUID_BC:
      case UVARY_BC:
      case VVARY_BC:
      case WVARY_BC:
      case U_PARABOLA_BC:
      case V_PARABOLA_BC:
      case W_PARABOLA_BC:
      case PRESSURE_USER_BC:
      case FLOW_PRESS_USER_BC:
      case T_USER_BC:
      case UUSER_BC:
      case VUSER_BC:
      case WUSER_BC:
      case QUSER_BC:
      case DX_USER_BC:
      case DY_USER_BC:
      case DZ_USER_BC:
      case P_LIQ_USER_BC:
      case SH_P_OPEN_USER_BC:
      case VAR_CA_USER_BC:
      case CA_EDGE_OR_FIX_BC:
      case YFLUX_USER_BC:
      case SURFACE_CHARGE_BC:
      case YUSER_BC:
      case FORCE_USER_BC:
      case FORCE_USER_RS_BC:
        BC_Types[bcid].u_BC[bcidx] = computed_value;
        break;
      default:
        BC_Types[bcid].BC_Data_Float[bcidx] = computed_value;
        break;
      }
    } else {
      return GOMA_ERROR;
    }
  }
  return GOMA_SUCCESS;
}
#endif
