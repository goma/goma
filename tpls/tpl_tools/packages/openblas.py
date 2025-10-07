from tpl_tools.packages import packages
from tpl_tools import utils
import os
import sys


class Package(packages.CMakePackage):
    def __init__(self):
        self.name = "openblas"
        self.version = "0.3.30"
        self.sha256 = "27342cff518646afb4c2b976d809102e368957974c250a25ccc965e53063c95d"
        self.filename = "OpenBLAS" + self.version + ".tar.gz"
        self.compression = "gz"
        self.url = (
            "https://github.com/OpenMathLib/OpenBLAS/releases/download/v"
            + self.version
            + "/OpenBLAS-"
            + self.version
            + ".tar.gz"
        )
        self.libraries = ["openblas"]
        self.dependencies = []

    def configure_options(self, builder):
        if builder.build_shared:
            builder.add_option("-DBUILD_SHARED_LIBS:BOOL=ON")
        else:
            builder.add_option("-DBUILD_SHARED_LIBS:BOOL=OFF")
            builder.add_option("-DBUILD_STATIC_LIBS:BOOL=ON")

        builder.add_option("-DDYNAMIC_ARCH:BOOL=ON")
        builder.add_option("-DUSE_THREAD:BOOL=OFF")
        builder.add_option("-DUSE_OPENMP:BOOL=OFF")
        if sys.platform == "darwin":
            builder.add_option("-DNO_SVE=1")
        else:
            builder.add_option("-DTARGET=GENERIC")
            builder.add_option("-DDYNAMIC_OLDER:BOOL=ON")

    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        ext = utils.get_library_extension(builder.build_shared)
        registry.set_environment_variable("OPENBLAS_DIR", builder.install_dir())

        # check if the library file exists in lib or lib64
        if os.path.exists(
            os.path.join(builder.install_dir(), "lib64/libopenblas" + ext)
        ):
            # if it exists in lib64, set the library path to lib64
            registry.set_environment_variable(
                "BLAS_LIBRARIES",
                os.path.join(builder.install_dir(), "lib64/libopenblas" + ext),
            )
            registry.set_environment_variable(
                "LAPACK_LIBRARIES",
                os.path.join(builder.install_dir(), "lib64/libopenblas" + ext),
            )
        else:
            # otherwise, set the library path to lib
            registry.set_environment_variable(
                "BLAS_LIBRARIES",
                os.path.join(builder.install_dir(), "lib/libopenblas" + ext),
            )
            registry.set_environment_variable(
                "LAPACK_LIBRARIES",
                os.path.join(builder.install_dir(), "lib/libopenblas" + ext),
            )
        registry.prepend_environment_variable(
            "CMAKE_PREFIX_PATH", builder.install_dir()
        )


# class Package(packages.GenericPackage):
#    def __init__(self):
#        self.name = "openblas"
#        self.version = "0.3.30"
#        self.sha256 = "27342cff518646afb4c2b976d809102e368957974c250a25ccc965e53063c95d"
#        self.filename = "OpenBLAS" + self.version + ".tar.gz"
#        self.compression = "gz"
#        self.url = (
#            "https://github.com/OpenMathLib/OpenBLAS/releases/download/v"
#            + self.version
#            + "/OpenBLAS-"
#            + self.version
#            + ".tar.gz"
#        )
#        self.libraries = ["openblas"]
#        self.dependencies = []
#
#    def build(self, builder):
#        self.build_command = ["make"]
#        if sys.platform == "darwin":
#            self.build_command.append("NO_SVE=1")
#        else:
#            self.build_command.append("TARGET=GENERIC")
#            self.build_command.append("DYNAMIC_OLDER=1")
#        self.build_command.append("USE_THREAD=0")
#        self.build_command.append("USE_OPENMP=0")
#        self.build_command.append("DYNAMIC_ARCH=1")
#        if builder.build_shared:
#            self.build_command.append("NO_SHARED=0")
#            self.build_command.append("shared")
#        else:
#            self.build_command.append("NO_SHARED=1")
#        self.build_command.append("libs")
#        self.build_command.append("netlib")
#        builder.run_command(self.build_command, jobs_flag="-j", parallel=True)
#
#    def install(self, builder):
#        command = self.build_command.copy()
#        command.append("PREFIX=" + builder.install_dir())
#        command.append("install")
#        builder.run_command(command)
#
#    def register(self, builder):
#        registry = builder._registry
#        registry.register_package(self.name, builder.install_dir())
#        ext = utils.get_library_extension(builder.build_shared)
#        registry.set_environment_variable("OPENBLAS_DIR", builder.install_dir())
#        registry.set_environment_variable(
#            "BLAS_LIBRARIES",
#            os.path.join(builder.install_dir(), "lib/libopenblas" + ext),
#        )
#        registry.set_environment_variable(
#            "LAPACK_LIBRARIES",
#            os.path.join(builder.install_dir(), "lib/libopenblas" + ext),
#        )
#        registry.prepend_environment_variable(
#            "CMAKE_PREFIX_PATH", builder.install_dir()
#        )
