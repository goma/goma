#ifdef HAVE_STRATIMIKOS

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"

#include "sl_util_structs.h"
#include "sl_stratimikos_interface.h"

#include "Thyra_SolveSupportTypes.hpp"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_VectorOut.h"


extern "C" {

int stratimikos_solve(struct Aztec_Linear_Solver_System *ams, double *x_,
    double *b_, int *iterations, char *stratimikos_file)
{
  using Teuchos::RCP;
  bool success = true;
  bool verbose = true;
  try {
    Epetra_Map map = ams->RowMatrix->RowMatrixRowMap();

    // Assign A with false so it doesn't get garbage collected.
    RCP<Epetra_RowMatrix> epetra_A = Teuchos::rcp(ams->RowMatrix, false);
    RCP<Epetra_Vector> epetra_x = Teuchos::rcp(
        new Epetra_Vector(Copy, map, x_));
    RCP<Epetra_Vector> epetra_b = Teuchos::rcp(
        new Epetra_Vector(Copy, map, b_));

    RCP<const Thyra::LinearOpBase<double> > A = Thyra::epetraLinearOp(epetra_A);
    RCP<Thyra::VectorBase<double> > x = Thyra::create_Vector(epetra_x,
        A->domain());
    RCP<const Thyra::VectorBase<double> > b = Thyra::create_Vector(epetra_b,
        A->range());

    Teuchos::RCP<Teuchos::FancyOStream> outstream =
        Teuchos::VerboseObjectBase::getDefaultOStream();

    // Get parameters from file
    RCP<Teuchos::ParameterList> solverParams;
    solverParams = Teuchos::getParametersFromXmlFile(stratimikos_file);

    // Set up base builder
    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
    linearSolverBuilder.setParameterList(solverParams);

    // set up solver factory using base/params
    RCP<Thyra::LinearOpWithSolveFactoryBase<double> > solverFactory =
        linearSolverBuilder.createLinearSolveStrategy("");

    // set output stream
    solverFactory->setOStream(outstream);

    // set solver verbosity
    solverFactory->setDefaultVerbLevel(Teuchos::VERB_NONE);

    RCP<Thyra::LinearOpWithSolveBase<double> > solver =
        Thyra::linearOpWithSolve(*solverFactory, A);

    Thyra::SolveStatus<double> status = Thyra::solve<double>(*solver,
        Thyra::NOTRANS, *b, x.ptr());

    x = Teuchos::null;

    *iterations = 1;
    if (!status.extraParameters.is_null()) {
      try {
        *iterations = status.extraParameters.get()->get<int> ("Iteration Count");
      } catch (const Teuchos::Exceptions::InvalidParameter &excpt) {}
    }

    if (status.solveStatus != Thyra::SOLVE_STATUS_CONVERGED) {
      *iterations = -1;
    }

    /* Convert solution vector */
    int NumMyRows = map.NumMyElements();

    Epetra_Vector *raw_x = epetra_x.get();
    for (int i = 0; i < NumMyRows; i++) {
      x_[i] = (*raw_x)[i];
    }

  } TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success)

  if (success) {
    return 0;
  } else {
    return -1;
  }
}

} /* End extern "C" */

#else /* HAVE_STRATIMIKOS */

#include "mpi.h"
extern "C" {
#include "std.h"
#include "mm_eh.h"

int stratimikos_solve(struct Aztec_Linear_Solver_System *ams, double *x_,
    double *b_, int *iterations, char *stratimikos_file)
{
  EH(-1, "Not built with stratimikos support!");
  return -1;
}

}
#endif /* HAVE_STRATIMIKOS */
