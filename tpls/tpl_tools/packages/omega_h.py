from tpl_tools.packages import packages
import os


class Package(packages.CMakePackage):
    def __init__(self):
        self.name = "omega-h"
        self.version = "11.2.0"
        self.sha256 = "741fc5c911682e2d46b6726ae48216e213ae5d5d1d1875d900ea26d636f8fb0b"
        self.filename = "omega-h-scorec-v" + self.version + ".tar.gz"
        self.url = (
            "https://github.com/SCOREC/omega_h/archive/refs/tags/scorec-v"
            + self.version
            + ".tar.gz"
        )
        self.executables = []
        self.libraries = ["omega_h"]
        self.includes = ["Omega_h_adapt.hpp", "Omega_h_mesh.hpp"]
        self.dependencies = ["openmpi", "cmake"]

    def set_environment(self, builder):
        builder.env = builder._registry.get_environment().copy()
        builder.env["CC"] = builder._registry.get_executable("mpicc")
        builder.env["CXX"] = builder._registry.get_executable("mpicxx")
        builder.env["FC"] = builder._registry.get_executable("mpifort")

    def configure_options(self, builder):
        if builder.build_shared:
            builder.add_option("-DBUILD_SHARED_LIBS:BOOL=ON")
        else:
            builder.add_option("-DBUILD_SHARED_LIBS:BOOL=OFF")

        CC = builder.env["CC"]
        CXX = builder.env["CXX"]
        FC = builder.env["FC"]
        builder.add_option("-DCMAKE_C_COMPILER=" + CC)
        builder.add_option("-DCMAKE_CXX_COMPILER=" + CXX)
        builder.add_option("-DCMAKE_Fortran_COMPILER=" + FC)
        # Goma adapts distributed meshes; a non-MPI Omega_h crashes in parallel
        builder.add_option("-DOmega_h_USE_MPI:BOOL=ON")
        builder.add_option("-DBUILD_TESTING:BOOL=OFF")

    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        registry.set_environment_variable("OMEGA_H_DIR", builder.install_dir())
        registry.prepend_environment_variable(
            "CMAKE_PREFIX_PATH", builder.install_dir()
        )
