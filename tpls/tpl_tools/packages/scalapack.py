from tpl_tools.packages import packages
from tpl_tools import utils


class Package(packages.CMakePackage):
    def __init__(self):
        self.name = "scalapack"
        self.version = "2.2.2"
        self.sha256 = "a2f0c9180a210bf7ffe126c9cb81099cf337da1a7120ddb4cbe4894eb7b7d022"
        self.filename = "scalapack-" + self.version + ".tar.gz"
        self.url = (
            "https://github.com/Reference-ScaLAPACK/scalapack/archive/refs/tags/v"
            + self.version
            + ".tar.gz"
        )
        self.libraries = ["scalapack"]
        self.dependencies = ["cmake", "openmpi", "lapack"]

    def set_environment(self, builder):
        builder.env = builder._registry.get_environment().copy()
        builder.env["CC"] = builder._registry.get_executable("mpicc")
        builder.env["CXX"] = builder._registry.get_executable("mpicxx")
        builder.env["FC"] = builder._registry.get_executable("mpifort")
        # check gcc version
        compiler, version = utils.check_gcc_clang_version(builder.env["CC"])
        if (
            compiler == "gcc"
            and version >= (14, 0, 0)
            or compiler == "clang"
            and version >= (16, 0, 0)
        ):
            builder.env["CFLAGS"] = "-Wno-implicit-function-declaration -O2"

    def configure_options(self, builder):
        if builder.build_shared:
            builder.add_option("-DBUILD_SHARED_LIBS:BOOL=ON")
        else:
            builder.add_option("-DBUILD_SHARED_LIBS:BOOL=OFF")
        builder.add_option("-DBLAS_LIBRARIES=" + builder.env["BLAS_LIBRARIES"])
        builder.add_option("-DLAPACK_LIBRARIES=" + builder.env["LAPACK_LIBRARIES"])
        builder.add_option("-DSCALAPACK_BUILD_TESTS=OFF")

    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        registry.set_environment_variable("SCALAPACK_DIR", builder.install_dir())
        registry.prepend_environment_variable(
            "CMAKE_PREFIX_PATH", builder.install_dir()
        )
