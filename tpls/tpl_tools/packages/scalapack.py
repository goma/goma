from tpl_tools.packages import packages


class Package(packages.CMakePackage):
    def __init__(self):
        self.name = "scalapack"
        self.version = "0234af94c6578c53ac4c19f2925eb6e5c4ad6f0f"
        self.sha256 = "6252706d4ed2e69aa9531730541eb0bec2b7332e9a08d39024aa2ec96b224a42"
        self.filename = "scalapack-" + self.version + ".tar.gz"
        self.url = (
            "https://github.com/Reference-ScaLAPACK/scalapack/archive/"
            + self.version
            + ".tar.gz"
        )
        self.libraries = ["scalapack"]

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
        builder.add_option("-D=SCALAPACK_BUILD_TESTS=OFF")

    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        registry.set_environment_variable("SCALAPACK_DIR", builder.install_dir())
        registry.append_environment_variable("CMAKE_PREFIX_PATH", builder.install_dir())
