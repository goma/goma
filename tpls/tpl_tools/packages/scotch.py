from tpl_tools.packages import packages
from tpl_tools import utils
import os


class Package(packages.CMakePackage):
    def __init__(self):
        self.name = "scotch"
        self.version = "7.0.4"
        self.sha256 = "8ef4719d6a3356e9c4ca7fefd7e2ac40deb69779a5c116f44da75d13b3d2c2c3"
        self.filename = "scotch-" + self.version + ".tar.gz"
        self.url = (
            "https://gitlab.inria.fr/scotch/scotch/-/archive/v"
            +self.version +"/scotch-v"
            + self.version
            + ".tar.gz"
        )
        self.executables = []
        self.libraries = ["scotch", "ptscotch"]
        self.includes = []

    def setDependencies(self, builder):
        builder.set_dependency("packages.openmpi")
        builder.set_dependency("packages.hdf5")
        builder.set_dependency("packages.pnetcdf")
        builder.set_dependency("packages.netcdf")
        builder.set_dependency("packages.fmt")
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
    
    def install(self, builder):
        cmake = builder._registry.get_executable("cmake")
        builder.run_command([cmake, "--install", "build_tpl"])
        mi = os.path.join(builder.install_dir(), "include/metis.h")
        pi = os.path.join(builder.install_dir(), "include/parmetis.h")
        for f in [mi,pi]:
            if os.path.exists(f):
                os.unlink(f)

    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        registry.set_environment_variable("SCOTCH_DIR", builder.install_dir())
        registry.append_environment_variable("CMAKE_PREFIX_PATH", builder.install_dir())
