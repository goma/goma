from tpl_tools.packages import packages
import os


class Package(packages.CMakePackage):
    def __init__(self):
        self.name = "hdf5"
        self.version = "1.14.6"
        version_tag = "1.14.6"
        self.sha256 = "e4defbac30f50d64e1556374aa49e574417c9e72c6b1de7a4ff88c4b1bea6e9b"
        self.filename = "hdf5-" + self.version + ".tar.bz2"
        self.url = (
            "https://github.com/HDFGroup/hdf5/releases/download/hdf5_"
            + self.version
            + "/hdf5-"
            + version_tag
            + ".tar.gz"
        )
        self.libraries = ["hdf5"]
        self.dependencies = ["openmpi"]

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
        builder.add_option("-DHDF5_ENABLE_PARALLEL:BOOL=ON")
        builder.add_option("-DDEFAULT_API_VERSION=V18")
        builder.add_option("-DHDF5_ENABLE_Z_LIB_SUPPORT:BOOL=ON")
        builder.add_option("-DHDF5_BUILD_FORTRAN:BOOL=OFF")
        builder.add_option("-DHDF5_BUILD_HL_LIB:BOOL=ON")
        builder.add_option("-DHDF5_DISABLE_COMPILER_WARNINGS:BOOL=ON")

        builder.add_option("-DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON")
        builder.add_option("-DBUILD_TESTING:BOOL=OFF")
        builder.add_option("-DDEFAULT_API_VERSION=V18")
        builder.add_option("-DHDF5_ENABLE_NONSTANDARD_FEATURE_FLOAT16=OFF")
        builder.add_option("-DHDF5_ENABLE_PARALLEL:BOOL=ON")
        builder.add_option("-DHDF5_ENABLE_Z_LIB_SUPPORT:BOOL=ON")
        builder.add_option("-DHDF5_BUILD_CPP_LIB:BOOL=OFF")
        builder.add_option("-DHDF5_BUILD_FORTRAN:BOOL=OFF")
        builder.add_option("-DHDF5_BUILD_HL_LIB:BOOL=ON")
        builder.add_option("-DHDF5_DISABLE_COMPILER_WARNINGS:BOOL=ON")

    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        registry.set_environment_variable("HDF5_DIR", builder.install_dir())
        registry.prepend_environment_variable(
            "CMAKE_PREFIX_PATH", builder.install_dir()
        )
        registry.prepend_environment_variable(
            "PATH", os.path.join(builder.install_dir(), "bin")
        )
