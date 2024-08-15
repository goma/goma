from tpl_tools.packages import packages
import os


class Package(packages.AutotoolsPackage):
    def __init__(self):
        self.name = "openmpi"
        self.version = "4.1.6"
        self.sha256 = "f740994485516deb63b5311af122c265179f5328a0d857a567b85db00b11e415"
        self.filename = "openmpi-" + self.version + "tar.bz2"
        self.url = (
            "https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-"
            + self.version
            + ".tar.bz2"
        )
        self.executables = ["mpicc", "mpifort", "mpicxx"]

    def setDependencies(self, builder):
        return

    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        registry.register_executable(
            os.path.join(builder.install_dir(), "bin", "mpicc")
        )
        registry.register_executable(
            os.path.join(builder.install_dir(), "bin", "mpicxx")
        )
        registry.register_executable(
            os.path.join(builder.install_dir(), "bin", "mpifort")
        )
        registry.set_environment_variable("MPI_HOME", builder.install_dir())
        registry.set_environment_variable(
            "MPI_C_COMPILER", os.path.join(builder.install_dir(), "bin", "mpicc")
        )
        registry.set_environment_variable(
            "MPI_CXX_COMPILER", os.path.join(builder.install_dir(), "bin", "mpicxx")
        )
        registry.set_environment_variable(
            "MPI_FORTRAN_COMPILER",
            os.path.join(builder.install_dir(), "bin", "mpifort"),
        )
        registry.set_environment_variable(
            "MPI_Fortran_COMPILER",
            os.path.join(builder.install_dir(), "bin", "mpifort"),
        )
        registry.append_environment_variable("CMAKE_PREFIX_PATH", builder.install_dir())
        registry.append_environment_variable(
            "PATH", os.path.join(builder.install_dir(), "bin")
        )
