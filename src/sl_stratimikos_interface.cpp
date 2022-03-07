#ifdef GOMA_ENABLE_STRATIMIKOS

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
#include "sl_epetra_interface.h"

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

extern "C" {

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
      solverParams_static[imtrx] = Teuchos::getParametersFromXmlFile(stratimikos_file[imtrx]);
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
