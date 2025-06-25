from tpl_tools.packages import packages
from tpl_tools import utils


class Package(packages.CMakePackage):
    def __init__(self):
        self.name = "hypre"
        self.version = "2.33.0"
        self.sha256 = "0f9103c34bce7a5dcbdb79a502720fc8aab4db9fd0146e0791cde7ec878f27da"
        self.filename = "hypre-" + self.version + ".tar.gz"
        self.url = (
            "https://github.com/hypre-space/hypre/archive/refs/tags/v"
            + self.version
            + ".tar.gz"
        )
        self.libraries = ["HYPRE"]
        self.dependencies = ["cmake", "openmpi", "lapack"]

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
        builder.add_option("-DBLAS_LIBRARIES=" + builder.env["BLAS_LIBRARIES"])
        builder.add_option("-DLAPACK_LIBRARIES=" + builder.env["LAPACK_LIBRARIES"])

    def configure(self, builder):
        cmake = builder._registry.get_executable("cmake")
        configure_options = [cmake, "-B", "build_tpl", "src"]
        configure_options.extend(builder.get_options())
        configure_options.append("-DCMAKE_INSTALL_PREFIX=" + builder.install_dir())
        configure_options.append("-DCMAKE_BUILD_TYPE=Release")
        configure_options.append("-DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE")
        builder.run_command(configure_options)

    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        registry.set_environment_variable("HYPRE_DIR", builder.install_dir())
        registry.prepend_environment_variable(
            "CMAKE_PREFIX_PATH", builder.install_dir()
        )
