from tpl_tools.packages import packages
import os


class Package(packages.CMakePackage):
    def __init__(self):
        self.name = "strumpack"
        self.version = "8.0.0"
        self.sha256 = "11cc8645d622a16510b39a20efc64f34862b41976152d17f9fbf3e91f899766c"
        self.filename = "strumpack-v" + self.version + ".tar.gz"
        self.url = (
            "https://github.com/pghysels/STRUMPACK/archive/v" + self.version + ".tar.gz"
        )
        self.libraries = ["strumpack"]
        self.includes = ["StrumpackConfig.h"]
        self.dependencies = ["cmake", "metis", "scotch", "parmetis", "openmpi"]

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

    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        registry.set_environment_variable("STRUMPACK_DIR", builder.install_dir())
        registry.prepend_environment_variable(
            "CMAKE_PREFIX_PATH", builder.install_dir()
        )
