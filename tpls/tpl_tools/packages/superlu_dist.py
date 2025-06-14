from tpl_tools.packages import packages
from tpl_tools import utils


class Package(packages.CMakePackage):
    def __init__(self):
        self.name = "superlu_dist"
        self.version = "9.1.0"
        self.sha256 = "1cb2c6dc7e8231b2ec30c1266e55e440ffca9f55527771d8df28f900dd179f9d"
        self.filename = "superlu_dist-" + self.version + ".tar.gz"
        self.url = (
            "https://github.com/xiaoyeli/superlu_dist/archive/refs/tags/v"
            + self.version
            + ".tar.gz"
        )
        self.libraries = ["superlu_dist"]
        self.includes = ["superlu_dist_config.h"]
        self.dependencies = ["cmake", "metis", "parmetis", "openmpi"]

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
        builder.add_option("-DTPL_ENABLE_Netcdf:BOOL=ON")
        builder.add_option("-DTPL_ENABLE_MPI:BOOL=ON")
        builder.add_option("-DTPL_ENABLE_X11:BOOL=ON")

        CC = builder.env["CC"]
        CXX = builder.env["CXX"]
        FC = builder.env["FC"]
        builder.add_option("-DCMAKE_C_COMPILER=" + CC)
        builder.add_option("-DCMAKE_CXX_COMPILER=" + CXX)
        builder.add_option("-DCMAKE_Fortran_COMPILER=" + FC)
        builder.add_option("-DCMAKE_C_STANDARD=11")
        builder.add_option("-DTPL_BLAS_LIBRARIES=" + builder.env["BLAS_LIBRARIES"])
        ext = utils.get_library_extension(builder.build_shared)
        builder.add_option(
            "-DTPL_PARMETIS_LIBRARIES="
            + builder.env["PARMETIS_DIR"]
            + "/lib/libparmetis"
            + ext
            + ";"
            + builder.env["METIS_DIR"]
            + "/lib/libmetis"
            + ext
        )
        builder.add_option("-Denable_openmp:BOOL=FALSE")
        builder.add_option(
            "-DTPL_PARMETIS_INCLUDE_DIRS="
            + builder.env["PARMETIS_DIR"]
            + "/include;"
            + builder.env["METIS_DIR"]
            + "/include"
        )

    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        registry.set_environment_variable("SUPERLU_DIST_DIR", builder.install_dir())
        registry.prepend_environment_variable(
            "CMAKE_PREFIX_PATH", builder.install_dir()
        )
