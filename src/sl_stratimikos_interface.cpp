#ifdef GOMA_ENABLE_STRATIMIKOS

#include <filesystem>
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

#include "Teuchos_YamlParameterListCoreHelpers.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_LinearOpWithSolveTester.hpp"
#include "sl_epetra_interface.h"

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
#include "sl_stratimikos_interface.h"
#include "sl_util_structs.h"

#ifdef GOMA_ENABLE_TPETRA
extern "C" {
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
  static RCP<Teuchos::ParameterList> solverParams_static[MAX_NUM_MATRICES];
  static bool param_set[MAX_NUM_MATRICES] = {false};
  static bool param_echo[MAX_NUM_MATRICES] = {false};

  try {
    RCP<const Tpetra::FECrsMatrix<double, LO, GO>> tpetra_A = tpetra_data->matrix;
    if (!tpetra_data->matrix->isFillComplete()) {
      tpetra_data->matrix->endAssembly();
    }

    RCP<Tpetra::Vector<double, LO, GO>> tpetra_x =
        rcp(new Tpetra::Vector<double, LO, GO>(tpetra_A->getDomainMap()));
    RCP<Tpetra::Vector<double, LO, GO>> tpetra_b =
        rcp(new Tpetra::Vector<double, LO, GO>(tpetra_A->getRangeMap()));

    auto n_rows = tpetra_data->matrix->getLocalNumRows();
    for (size_t i = 0; i < n_rows; i++) {
      tpetra_x->replaceGlobalValue(matrix->global_ids[i], x_[i]);
      tpetra_b->replaceGlobalValue(matrix->global_ids[i], b_[i]);
    }
#if 0
    Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<double, LO, GO>>::writeSparseFile("A.mm",
                                                                                     tpetra_A);
    Tpetra::MatrixMarket::Writer<Tpetra::Vector<double, LO, GO>>::writeDenseFile("x.mm", tpetra_x);
    Tpetra::MatrixMarket::Writer<Tpetra::Vector<double, LO, GO>>::writeDenseFile("b.mm", tpetra_b);
#endif

    RCP<const Thyra::LinearOpBase<double>> A = Thyra::createConstLinearOp(
        Teuchos::rcp_dynamic_cast<const Tpetra::Operator<double, LO, GO>>(tpetra_A));

    RCP<Thyra::VectorBase<double>> x = Thyra::createVector(tpetra_x);
    RCP<const Thyra::VectorBase<double>> b = Thyra::createVector(tpetra_b);

    Teuchos::RCP<Teuchos::FancyOStream> outstream = Teuchos::VerboseObjectBase::getDefaultOStream();

    // Get parameters from file
    if (!param_set[imtrx]) {
      param_set[imtrx] = true;
      std::filesystem::path path(stratimikos_file[imtrx]);
      if (path.extension() == ".yaml") {
        solverParams_static[imtrx] = Teuchos::getParametersFromYamlFile(stratimikos_file[imtrx]);
      } else {
        solverParams_static[imtrx] = Teuchos::getParametersFromXmlFile(stratimikos_file[imtrx]);
      }
    }

    RCP<Teuchos::ParameterList> solverParams = solverParams_static[imtrx];

    // Set up base builder
    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;

#ifdef GOMA_ENABLE_TEKO
    Teko::addTekoToStratimikosBuilder(linearSolverBuilder);
#endif

    linearSolverBuilder.setParameterList(solverParams);

    // set up solver factory using base/params
    RCP<Thyra::LinearOpWithSolveFactoryBase<double>> solverFactory =
        linearSolverBuilder.createLinearSolveStrategy("");

    if (!param_echo[imtrx]) {
      std::string echo_file(stratimikos_file[imtrx]);
      echo_file = "echo_" + echo_file;
      linearSolverBuilder.writeParamsFile(*solverFactory, echo_file);
      param_echo[imtrx] = true;
      if (imtrx == 0) {
        auto valid_params = linearSolverBuilder.getValidParameters();
        Teuchos::writeParameterListToYamlFile(*valid_params, "stratimikos_valid_params.yaml");
      }
    }

    // set output stream
    solverFactory->setOStream(outstream);

    // set solver verbosity
    solverFactory->setDefaultVerbLevel(Teuchos::VERB_NONE);

    RCP<Thyra::LinearOpWithSolveBase<double>> solver = Thyra::linearOpWithSolve(*solverFactory, A);

    Thyra::SolveStatus<double> status = Thyra::solve<double>(*solver, Thyra::NOTRANS, *b, x.ptr());

    x = Teuchos::null;

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
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success)
  tpetra_data->matrix->beginAssembly();

  if (success) {
    return 0;
  } else {
    return -1;
  }
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
  static RCP<Teuchos::ParameterList> solverParams_static[MAX_NUM_MATRICES];
  static bool param_set[MAX_NUM_MATRICES] = {false};
  static bool param_echo[MAX_NUM_MATRICES] = {false};

  try {
    Epetra_Map map = ams->RowMatrix->RowMatrixRowMap();

    // Assign A with false so it doesn't get garbage collected.
    RCP<Epetra_RowMatrix> epetra_A = Teuchos::rcp(ams->RowMatrix, false);
    RCP<Epetra_Vector> epetra_x = Teuchos::rcp(new Epetra_Vector(Copy, map, x_));
    RCP<Epetra_Vector> epetra_b = Teuchos::rcp(new Epetra_Vector(Copy, map, b_));

    RCP<const Thyra::LinearOpBase<double>> A = Thyra::epetraLinearOp(epetra_A);
    RCP<Thyra::VectorBase<double>> x = Thyra::create_Vector(epetra_x, A->domain());
    RCP<const Thyra::VectorBase<double>> b = Thyra::create_Vector(epetra_b, A->range());

    Teuchos::RCP<Teuchos::FancyOStream> outstream = Teuchos::VerboseObjectBase::getDefaultOStream();

    // Get parameters from file
    if (!param_set[imtrx]) {
      param_set[imtrx] = true;
      std::filesystem::path path(stratimikos_file[imtrx]);
      if (path.extension() == ".yaml") {
        solverParams_static[imtrx] = Teuchos::getParametersFromYamlFile(stratimikos_file[imtrx]);
      } else {
        solverParams_static[imtrx] = Teuchos::getParametersFromXmlFile(stratimikos_file[imtrx]);
      }
    }

    RCP<Teuchos::ParameterList> solverParams = solverParams_static[imtrx];

    // Set up base builder
    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;

#ifdef GOMA_ENABLE_TEKO
    Teko::addTekoToStratimikosBuilder(linearSolverBuilder);
#endif

    linearSolverBuilder.setParameterList(solverParams);

    // set up solver factory using base/params
    RCP<Thyra::LinearOpWithSolveFactoryBase<double>> solverFactory =
        linearSolverBuilder.createLinearSolveStrategy("");

    if (!param_echo[imtrx]) {
      std::string echo_file(stratimikos_file[imtrx]);
      echo_file = "echo_" + echo_file;
      linearSolverBuilder.writeParamsFile(*solverFactory, echo_file);
      param_echo[imtrx] = true;
    }

    // set output stream
    solverFactory->setOStream(outstream);

    // set solver verbosity
    solverFactory->setDefaultVerbLevel(Teuchos::VERB_NONE);

    RCP<Thyra::LinearOpWithSolveBase<double>> solver = Thyra::linearOpWithSolve(*solverFactory, A);

    Thyra::SolveStatus<double> status = Thyra::solve<double>(*solver, Thyra::NOTRANS, *b, x.ptr());

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

extern "C" {
#include "mm_eh.h"
#include "std.h"

int stratimikos_solve(struct Aztec_Linear_Solver_System *ams,
                      double *x_,
                      double *b_,
                      int *iterations,
                      char *stratimikos_file) {
  GOMA_EH(GOMA_ERROR, "Not built with stratimikos support!");
  return -1;
}
}
#endif /* GOMA_ENABLE_STRATIMIKOS */
