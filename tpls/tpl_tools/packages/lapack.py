from tpl_tools.packages import packages
from tpl_tools import utils
import os


class Package(packages.CMakePackage):
    def __init__(self):
        self.name = "lapack"
        self.version = "3.12.0"
        self.sha256 = "eac9570f8e0ad6f30ce4b963f4f033f0f643e7c3912fc9ee6cd99120675ad48b"
        self.filename = "lapack-" + self.version + ".tar.gz"
        self.url = (
            "https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v"
            + self.version
            + ".tar.gz"
        )
        self.libraries = ["lapack", "blas"]

    def configure_options(self, builder):
        if builder.build_shared:
            builder.add_option("-DBUILD_SHARED_LIBS:BOOL=ON")
        else:
            builder.add_option("-DBUILD_SHARED_LIBS:BOOL=OFF")


    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        registry.set_environment_variable("LAPACK_DIR", builder.install_dir())
        ext = ".a"
        if builder.build_shared:
            ext = ".so"
        registry.set_environment_variable(
            "BLAS_LIBRARIES",
            os.path.join(builder.install_dir(), "lib/libblas" + ext),
        )
        registry.set_environment_variable(
            "LAPACK_LIBRARIES",
            os.path.join(builder.install_dir(), "lib/liblapack" + ext),
        )
        registry.append_environment_variable("CMAKE_PREFIX_PATH", builder.install_dir())
