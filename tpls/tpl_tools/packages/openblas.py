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