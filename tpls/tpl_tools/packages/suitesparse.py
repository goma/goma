from tpl_tools.packages import packages


class Package(packages.CMakePackage):
    def __init__(self):
        self.name = "suitesparse"
        self.version = "7.7.0"
        self.sha256 = "529b067f5d80981f45ddf6766627b8fc5af619822f068f342aab776e683df4f3"
        self.filename = "suitesparse-" + self.version + ".tar.gz"
        self.url = (
            "https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v"
            + self.version
            + ".tar.gz"
        )
        self.libraries = ["umfpack", "suitesparseconfig"]
        self.includes = ["suitesparse/umfpack.h"]

    def setDependencies(self, builder):
        builder.set_dependency("packages.openmpi")
        builder.set_dependency("packages.openblas")
        builder.set_dependency("packages.metis")
        builder.set_dependency("packages.parmetis")
        return

    def configure_options(self, builder):
        if builder.build_shared:
            builder.add_option("-D=BUILD_SHARED_LIBS:BOOL=ON")
        else:
            builder.add_option("-D=BUILD_SHARED_LIBS:BOOL=OFF")
        builder.add_option("-D=BLAS_LIBRARIES=" + builder.env["BLAS_LIBRARIES"])
        builder.add_option("-D=LAPACK_LIBRARIES=" + builder.env["LAPACK_LIBRARIES"])
        builder.add_option(
            "-D=SUITESPARSE_ENABLE_PROJECTS=suitesparse_config;amd;camd;ccolamd;colamd;cholmod;umfpack"
        )

    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        registry.set_environment_variable("SUITESPARSE_DIR", builder.install_dir())
        registry.append_environment_variable("CMAKE_PREFIX_PATH", builder.install_dir())
