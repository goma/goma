from tpl_tools.packages import packages
from tpl_tools import utils
import os


class Package(packages.CMakePackage):
    def __init__(self):
        self.name = "lapack"
        self.version = "3.12.1"
        self.sha256 = "2ca6407a001a474d4d4d35f3a61550156050c48016d949f0da0529c0aa052422"
        self.filename = "lapack-" + self.version + ".tar.gz"
        self.url = (
            "https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v"
            + self.version
            + ".tar.gz"
        )
        self.libraries = ["lapack", "blas"]
        self.dependencies = ["cmake"]

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
        libdir = 'lib'
        if os.path.exists(os.path.join(builder.install_dir(), "lib64/libblas" + ext)):
            libdir = 'lib64'

        registry.set_environment_variable(
            "BLAS_LIBRARIES",
            os.path.join(builder.install_dir(), libdir+"/libblas" + ext),
        )
        registry.set_environment_variable(
            "LAPACK_LIBRARIES",
            os.path.join(builder.install_dir(), libdir+"/liblapack" + ext),
        )
        registry.prepend_environment_variable(
            "CMAKE_PREFIX_PATH", builder.install_dir()
        )
