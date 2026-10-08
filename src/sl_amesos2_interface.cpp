#ifdef GOMA_ENABLE_AMESOS2

#include <Amesos2.hpp>
#include <Amesos2_Version.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterListExceptions.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>
#include <Teuchos_YamlParameterListCoreHelpers.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_FECrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <filesystem>

#include "linalg/sparse_matrix.h"
#include "linalg/sparse_matrix_tpetra.h"
#include "sl_amesos2_interface.h"

using MAT = Tpetra::CrsMatrix<double, LO, GO>;
using VEC = Tpetra::Vector<double, LO, GO>;
using MV = Tpetra::MultiVector<double, LO, GO>;

struct Amesos2_Solver_Data {
  Teuchos::RCP<Amesos2::Solver<MAT, MV>> solver;

  Amesos2_Solver_Data() { solver = Teuchos::null; }
};

extern "C" void amesos2_solver_destroy(struct GomaLinearSolverData *ams) {
  auto solver_data = static_cast<Amesos2_Solver_Data *>(ams->SolverData);
  delete solver_data;
}

extern "C" {
int amesos2_solve(struct GomaLinearSolverData *ams,
                  double *x_,
                  double *b_,
                  char *amesos2_solver,
                  char *amesos2_file) {
  using Teuchos::RCP;
  auto matrix = static_cast<GomaSparseMatrix>(ams->GomaMatrixData);
  auto *tpetra_data = static_cast<TpetraSparseMatrix *>(matrix->data);
  bool success = true;
  bool verbose = true;

  if (ams->SolverData == NULL) {
    ams->SolverData = new Amesos2_Solver_Data();
    ams->DestroySolverData = amesos2_solver_destroy;
  }
  auto solver_data = static_cast<Amesos2_Solver_Data *>(ams->SolverData);

  try {
    if (!tpetra_data->matrix->isFillComplete()) {
      tpetra_data->matrix->endAssembly();
    }

    RCP<VEC> tpetra_x = rcp(new VEC(tpetra_data->matrix->getDomainMap()));
    RCP<VEC> tpetra_b = rcp(new VEC(tpetra_data->matrix->getRangeMap()));

    for (size_t i = 0; i < tpetra_x->getLocalLength(); i++) {
      tpetra_x->replaceGlobalValue(matrix->global_ids[i], x_[i]);
    }
    for (size_t i = 0; i < tpetra_b->getLocalLength(); i++) {
      tpetra_b->replaceGlobalValue(matrix->global_ids[i], b_[i]);
    }

    if (solver_data->solver.is_null()) {
      Teuchos::RCP<MAT> crs_matrix = Teuchos::rcp_dynamic_cast<MAT>(tpetra_data->matrix);
      solver_data->solver = Amesos2::create<MAT, MV>(amesos2_solver, crs_matrix);
      if (amesos2_file != NULL && strlen(amesos2_file) > 0) {
        std::filesystem::path path(amesos2_file);
        Teuchos::RCP<Teuchos::ParameterList> amesos2_params;
        if (path.extension() == ".yaml") {
          amesos2_params = Teuchos::getParametersFromYamlFile(amesos2_file);
        } else {
          amesos2_params = Teuchos::getParametersFromXmlFile(amesos2_file);
        }
        solver_data->solver->setParameters(amesos2_params);
      }

      solver_data->solver->symbolicFactorization();
    }
    solver_data->solver->numericFactorization();

    solver_data->solver->solve(tpetra_x.ptr(), tpetra_b.ptr());

    Teuchos::RCP<Teuchos::FancyOStream> outstream = Teuchos::VerboseObjectBase::getDefaultOStream();

    /* Convert solution vector */
    int NumMyRows = matrix->n_rows;

    auto x_data = tpetra_x->getData();
    for (int i = 0; i < NumMyRows; i++) {
      x_[i] = x_data[i];
    }
    tpetra_data->matrix->beginAssembly();
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success)
  if (success) {
    return 0;
  }

  return 1;
}

} /* End extern "C" */

#else /* GOMA_ENABLE_AMESOS2 */

#include "mpi.h"

extern "C" {
#include "mm_eh.h"
#include "std.h"

int amesos2_solve(struct GomaLinearSolverData *ams,
                  double *x_,
                  double *b_,
                  int *iterations,
                  char *stratimikos_file) {
  GOMA_EH(GOMA_ERROR, "Not built with Amesos2 support!");
  return -1;
}
}
#endif /* GOMA_ENABLE_AMESOS2 */
