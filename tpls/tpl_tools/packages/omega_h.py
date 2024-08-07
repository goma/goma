from tpl_tools.packages import packages
import os

class Package(packages.CMakePackage):
    def __init__(self):
        self.name = 'omega-h'
        self.version = "scorec-v10.8.4"
        self.sha256 = "a8c85314acdc30e4680d97191674d3ab6f15ee3a59cbfe169c3bf541aa9cf9d8"
        self.filename = "omega-h-" + self.version + ".tar.gz"
        self.url = (
            "https://github.com/SCOREC/omega_h/archive/refs/tags/"
            + self.version
            + ".tar.gz"
        )
        self.executables = []
        self.libraries = ["omega_h"]
        self.includes = ["Omega_h_adapt.hpp", "Omega_h_mesh.hpp"]

    def setDependencies(self, builder):
        builder.set_dependency('packages.openmpi')
        return

    def set_environment(self, builder):
        builder.env = builder._registry.get_environment().copy()
        builder.env["CC"] = builder._registry.get_executable("mpicc")
        builder.env["CXX"] = builder._registry.get_executable("mpicxx")
        builder.env["FC"] = builder._registry.get_executable("mpifort")

    def configure_options(self, builder):
        if builder.build_shared:
            builder.add_option("-D=BUILD_SHARED_LIBS:BOOL=ON")
        else:
            builder.add_option("-D=BUILD_SHARED_LIBS:BOOL=OFF")

        CC = builder.env["CC"]
        CXX = builder.env["CXX"]
        FC = builder.env["FC"]
        builder.add_option("-DCMAKE_C_COMPILER=" + CC)
        builder.add_option("-DCMAKE_CXX_COMPILER=" + CXX)
        builder.add_option("-DCMAKE_Fortran_COMPILER=" + FC)

    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        registry.set_environment_variable("OMEGA_H_DIR", builder.install_dir())
        registry.append_environment_variable("CMAKE_PREFIX_PATH", builder.install_dir())