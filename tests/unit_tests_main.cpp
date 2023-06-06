#define CATCH_CONFIG_RUNNER
#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <mpi.h>

extern "C" {
int ProcID = 0;
int parallel_err = 0;
int Num_Proc = 1;
}

int main(int argc, char *argv[]) {
  // global setup...
  MPI_Init(&argc, &argv);

  auto error = MPI_Comm_size(MPI_COMM_WORLD, &Num_Proc);
  if (error) {
    return error;
  }
  error = MPI_Comm_rank(MPI_COMM_WORLD, &ProcID);
  if (error) {
    return error;
  }

  if (Num_Proc > 1) {
    if (ProcID == 0) {
      std::cerr << "Unit tests are currently expected to be run serially, found " << Num_Proc
                << " MPI Processes\n";
    }
    return -1;
  }

  int result = Catch::Session().run(argc, argv);

  // global clean-up...
  MPI_Finalize();

  return result;
}
