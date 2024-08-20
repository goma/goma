#include <cstddef>
#ifdef GOMA_ENABLE_STRATIMIKOS

#ifdef GOMA_ENABLE_CPP_FILESYSTEM
#include <filesystem>
#endif
#include <iostream>
#include <string>
#include <utility>

#include "Epetra_DataAccess.h"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Teuchos_ENull.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterListExceptions.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_VerbosityLevel.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "Teuchos_YamlParameterListCoreHelpers.hpp"
#include "Teuchos_config.h"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_LinearOpBase_decl.hpp"
#include "Thyra_LinearOpWithSolveBase_decl.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase_decl.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_OperatorVectorTypes.hpp"
#include "Thyra_SolveSupportTypes.hpp"
#include "Thyra_VectorBase.hpp"

#include "Thyra_LinearOpTester.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_LinearOpWithSolveTester.hpp"

#ifdef GOMA_ENABLE_TPETRA
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_TpetraVector.hpp"
#include "linalg/sparse_matrix_tpetra.h"
#endif

#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_VectorOut.h"

#ifdef GOMA_ENABLE_TEKO
// Teko-Package includes
#include "Teko_StratimikosFactory.hpp"
#endif

#ifdef HAVE_MPI
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Vector.h"
#include "linalg/sparse_matrix.h"
#include "linalg/sparse_matrix_epetra.h"
#include "sl_stratimikos_interface.h"
#include "sl_util_structs.h"

struct Stratimikos_Solver_Data {
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<double>> solver;
  Teuchos::RCP<Teuchos::ParameterList> solverParams;
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double>> solverFactory;
  Teuchos::RCP<const Thyra::LinearOpBase<double>> A;

  Stratimikos_Solver_Data() {
    solver = Teuchos::null;
    solverParams = Teuchos::null;
    solverFactory = Teuchos::null;
  }
};

static std::string get_file_extension(const std::string &filename) {
#ifdef GOMA_ENABLE_CPP_FILESYSTEM
  std::filesystem::path path(stratimikos_file[imtrx]);
  return path.extension();
#else
  size_t pos = filename.find_last_of('.');
  if (pos != std::string::npos) {
    return filename.substr(pos + 1);
  }
  return "";
#endif
}

extern "C" void stratimikos_solver_destroy(struct GomaLinearSolverData *ams) {
  auto solver_data = static_cast<Stratimikos_Solver_Data *>(ams->SolverData);
  solver_data->solver = Teuchos::null;
  solver_data->solverParams = Teuchos::null;
  delete solver_data;
}

using Teuchos::RCP;

static void stratimikos_solve_setup(RCP<const Thyra::LinearOpBase<double>> A,
                                    Stratimikos_Solver_Data *solver_data,
                                    std::string stratimikos_file,
                                    bool echo_params) {

  if (solver_data->solver.is_null()) {
    // Set up solver (only once per matrix)

    RCP<Teuchos::ParameterList> solverParams = solver_data->solverParams;

    // Set up base builder
    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;

#ifdef GOMA_ENABLE_TEKO
    Teko::addTekoToStratimikosBuilder(linearSolverBuilder);
#endif

    linearSolverBuilder.setParameterList(solverParams);

    auto valid_params = linearSolverBuilder.getValidParameters();
    Teuchos::writeParameterListToYamlFile(*valid_params, "valid_params.yaml");

    // set up solver factory using base/params
    RCP<Thyra::LinearOpWithSolveFactoryBase<double>> solverFactory =
        linearSolverBuilder.createLinearSolveStrategy("");
    solver_data->solverFactory = solverFactory;

    if (echo_params) {
      std::string echo_file(stratimikos_file);
      echo_file = "echo_" + echo_file;
      linearSolverBuilder.writeParamsFile(*solverFactory, echo_file);
#ifdef GOMA_STRATIMIKOS_WRITE_VALID_PARAMS
      if (imtrx == 0) {
        auto valid_params = linearSolverBuilder.getValidParameters();
        Teuchos::writeParameterListToYamlFile(*valid_params, "stratimikos_valid_params.yaml");
      }
#endif
    }

    Teuchos::RCP<Teuchos::FancyOStream> outstream = Teuchos::VerboseObjectBase::getDefaultOStream();
    // set output stream
    solverFactory->setOStream(outstream);

    // set solver verbosity
    solverFactory->setDefaultVerbLevel(Teuchos::VERB_NONE);

    solver_data->solver = solverFactory->createOp();
    Thyra::initializeOp(*(solver_data->solverFactory), A, solver_data->solver.ptr());
  } else {
    Thyra::initializeAndReuseOp(*(solver_data->solverFactory), A, solver_data->solver.ptr());
  }
}

extern "C" {
#ifdef GOMA_ENABLE_TPETRA
int stratimikos_solve_tpetra(struct GomaLinearSolverData *ams,
                             double *x_,
                             double *b_,
                             int *iterations,
                             char stratimikos_file[MAX_NUM_MATRICES][MAX_CHAR_IN_INPUT],
                             int imtrx) {
  using Teuchos::RCP;
  auto matrix = static_cast<GomaSparseMatrix>(ams->GomaMatrixData);
  auto *tpetra_data = static_cast<TpetraSparseMatrix *>(matrix->data);
  bool success = true;
  bool verbose = true;
  static bool param_echo[MAX_NUM_MATRICES] = {false};

  if (ams->SolverData == NULL) {
    ams->SolverData = new Stratimikos_Solver_Data();
    ams->DestroySolverData = stratimikos_solver_destroy;
  }
  auto solver_data = static_cast<Stratimikos_Solver_Data *>(ams->SolverData);

  try {
    RCP<const Tpetra::FECrsMatrix<double, LO, GO>> tpetra_A = tpetra_data->matrix;
    if (!tpetra_data->matrix->isFillComplete()) {
      tpetra_data->matrix->endAssembly();
    }

    RCP<Tpetra::Vector<double, LO, GO>> tpetra_x =
        rcp(new Tpetra::Vector<double, LO, GO>(tpetra_A->getDomainMap()));
    RCP<Tpetra::Vector<double, LO, GO>> tpetra_b =
        rcp(new Tpetra::Vector<double, LO, GO>(tpetra_A->getRangeMap()));

    for (size_t i = 0; i < tpetra_x->getLocalLength(); i++) {
      tpetra_x->replaceGlobalValue(matrix->global_ids[i], x_[i]);
    }
    for (size_t i = 0; i < tpetra_b->getLocalLength(); i++) {
      tpetra_b->replaceGlobalValue(matrix->global_ids[i], b_[i]);
    }
#if 0
    Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<double, LO, GO>>::writeSparseFile("A.mm",
                                                                                     tpetra_A);
    Tpetra::MatrixMarket::Writer<Tpetra::Vector<double, LO, GO>>::writeDenseFile("x.mm", tpetra_x);
    Tpetra::MatrixMarket::Writer<Tpetra::Vector<double, LO, GO>>::writeDenseFile("b.mm", tpetra_b);
#endif

    solver_data->A = Thyra::createConstLinearOp(
        Teuchos::rcp_dynamic_cast<const Tpetra::Operator<double, LO, GO>>(tpetra_A));

    RCP<Thyra::VectorBase<double>> x = Thyra::createVector(tpetra_x);
    RCP<const Thyra::VectorBase<double>> b = Thyra::createVector(tpetra_b);

    Teuchos::RCP<Teuchos::FancyOStream> outstream = Teuchos::VerboseObjectBase::getDefaultOStream();

    // Get parameters from file
    if (solver_data->solverParams.is_null()) {
      if (get_file_extension(stratimikos_file[imtrx]) == ".yaml") {
        solver_data->solverParams = Teuchos::getParametersFromYamlFile(stratimikos_file[imtrx]);
      } else {
        solver_data->solverParams = Teuchos::getParametersFromXmlFile(stratimikos_file[imtrx]);
      }
    }

    stratimikos_solve_setup(solver_data->A, solver_data, stratimikos_file[imtrx],
                            param_echo[imtrx]);
    param_echo[imtrx] = true;

    Thyra::SolveStatus<double> status =
        Thyra::solve<double>(*(solver_data->solver), Thyra::NOTRANS, *b, x.ptr());

    *iterations = 1;
    if (!status.extraParameters.is_null()) {
      try {
        *iterations = status.extraParameters.get()->get<int>("Iteration Count");
      } catch (const Teuchos::Exceptions::InvalidParameter &excpt) {
      }
    }

    /* Convert solution vector */
    int NumMyRows = matrix->n_rows;

    auto x_data = tpetra_x->getData();
    for (int i = 0; i < NumMyRows; i++) {
      x_[i] = x_data[i];
    }
    x = Teuchos::null;
    Thyra::uninitializeOp(*(solver_data->solverFactory), solver_data->solver.ptr());
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success)
  tpetra_data->matrix->beginAssembly();

  if (success) {
    return 0;
  } else {
    return -1;
  }
}
#else  /* GOMA_ENABLE_TPETRA */
int stratimikos_solve_tpetra(struct GomaLinearSolverData *ams,
                             double *x_,
                             double *b_,
                             int *iterations,
                             char stratimikos_file[MAX_NUM_MATRICES][MAX_CHAR_IN_INPUT],
                             int imtrx) {
  GOMA_EH(GOMA_ERROR, "Not built with Tpetra Stratimikos support!");
  return -1;
}
#endif /* GOMA_ENABLE_TPETRA */

int stratimikos_solve(struct GomaLinearSolverData *ams,
                      double *x_,
                      double *b_,
                      int *iterations,
                      char stratimikos_file[MAX_NUM_MATRICES][MAX_CHAR_IN_INPUT],
                      int imtrx) {
  using Teuchos::RCP;
  bool success = true;
  bool verbose = true;
  static bool param_echo[MAX_NUM_MATRICES] = {false};

  GomaSparseMatrix matrix = (GomaSparseMatrix)ams->GomaMatrixData;
  EpetraSparseMatrix *epetra_matrix = static_cast<EpetraSparseMatrix *>(matrix->data);

  if (ams->SolverData == NULL) {
    ams->SolverData = new Stratimikos_Solver_Data();
    ams->DestroySolverData = stratimikos_solver_destroy;
  }
  auto solver_data = static_cast<Stratimikos_Solver_Data *>(ams->SolverData);
  try {
    Epetra_Map map = epetra_matrix->matrix->RowMatrixRowMap();

    // Assign A with false so it doesn't get garbage collected.
    RCP<Epetra_CrsMatrix> epetra_A = epetra_matrix->matrix;
    RCP<Epetra_Vector> epetra_x = Teuchos::rcp(new Epetra_Vector(Copy, map, x_));
    RCP<Epetra_Vector> epetra_b = Teuchos::rcp(new Epetra_Vector(Copy, map, b_));

    solver_data->A = Thyra::epetraLinearOp(epetra_A);
    RCP<Thyra::VectorBase<double>> x = Thyra::create_Vector(epetra_x, solver_data->A->domain());
    RCP<const Thyra::VectorBase<double>> b =
        Thyra::create_Vector(epetra_b, solver_data->A->range());

    Teuchos::RCP<Teuchos::FancyOStream> outstream = Teuchos::VerboseObjectBase::getDefaultOStream();

    // Get parameters from file
    if (solver_data->solverParams.is_null()) {
      if (get_file_extension(stratimikos_file[imtrx]) == ".yaml") {
        solver_data->solverParams = Teuchos::getParametersFromYamlFile(stratimikos_file[imtrx]);
      } else {
        solver_data->solverParams = Teuchos::getParametersFromXmlFile(stratimikos_file[imtrx]);
      }
    }

    stratimikos_solve_setup(solver_data->A, solver_data, stratimikos_file[imtrx],
                            param_echo[imtrx]);
    param_echo[imtrx] = true;

    Thyra::SolveStatus<double> status =
        Thyra::solve<double>(*(solver_data->solver), Thyra::NOTRANS, *b, x.ptr());

    x = Teuchos::null;

    *iterations = 1;
    if (!status.extraParameters.is_null()) {
      try {
        *iterations = status.extraParameters.get()->get<int>("Iteration Count");
      } catch (const Teuchos::Exceptions::InvalidParameter &excpt) {
      }
    }

    /* Convert solution vector */
    int NumMyRows = map.NumMyElements();

    Epetra_Vector *raw_x = epetra_x.get();
    for (int i = 0; i < NumMyRows; i++) {
      x_[i] = (*raw_x)[i];
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success)

  if (success) {
    return 0;
  } else {
    return -1;
  }
}

} /* End extern "C" */

#else /* GOMA_ENABLE_STRATIMIKOS */

#include "mpi.h"

#include "sl_stratimikos_interface.h"
extern "C" {
#include "mm_eh.h"
#include "std.h"
int stratimikos_solve_tpetra(struct GomaLinearSolverData *ams,
                             double *x_,
                             double *b_,
                             int *iterations,
                             char stratimikos_file[MAX_NUM_MATRICES][MAX_CHAR_IN_INPUT],
                             int imtrx) {
  GOMA_EH(GOMA_ERROR, "Not built with stratimikos support!");
  return -1;
}

int stratimikos_solve(struct GomaLinearSolverData *ams,
                      double *x_,
                      double *b_,
                      int *iterations,
                      char stratimikos_file[MAX_NUM_MATRICES][MAX_CHAR_IN_INPUT],
                      int imtrx) {
  GOMA_EH(GOMA_ERROR, "Not built with stratimikos support!");
  return -1;
}
}
#endif /* GOMA_ENABLE_STRATIMIKOS */
