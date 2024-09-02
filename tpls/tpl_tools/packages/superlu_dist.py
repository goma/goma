from tpl_tools.packages import packages


class Package(packages.CMakePackage):
    def __init__(self):
        self.name = "superlu_dist"
        self.version = "9.0.0"
        self.sha256 = "aa43d33d4b1b0f5f7b5ad7685e9a6bc25088832c6c74d2ab8f75a2c9f4e9e955"
        self.filename = "superlu_dist-" + self.version + ".tar.gz"
        self.url = (
            "https://github.com/xiaoyeli/superlu_dist/archive/refs/tags/v"
            + self.version
            + ".tar.gz"
        )
        self.libraries = ["superlu_dist"]
        self.includes = ["superlu_dist_config.h"]

    def setDependencies(self, builder):
        builder.set_dependency("packages.openmpi")
        builder.set_dependency("packages.openblas")
        builder.set_dependency("packages.metis")
        builder.set_dependency("packages.parmetis")
        return

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
        builder.add_option("-DPNetCDF_ROOT:PATH=" + builder.env["PNETCDF_DIR"])
        builder.add_option("-DNetCDF_ROOT:PATH=" + builder.env["NETCDF_DIR"])
        builder.add_option("-DnetCDF_ROOT:PATH=" + builder.env["NETCDF_DIR"])
        builder.add_option("-DHDF5_ROOT:PATH=" + builder.env["HDF5_DIR"])
        builder.add_option("-DHDF5_DIR:PATH=" + builder.env["HDF5_DIR"])
        builder.add_option("-DTPL_ENABLE_Matio=OFF")
        builder.add_option("-DSeacas_ENABLE_ALL_PACKAGES:BOOL=ON")
        builder.add_option("-DSeacas_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON")
        builder.add_option("-DSeacas_ENABLE_SECONDARY_TESTED_CODE:BOOL=ON")
        builder.add_option("-DTPL_BLAS_LIBRARIES=" + builder.env["BLAS_LIBRARIES"])
        ext = ".a"
        if builder.build_shared:
            ext = ".so"
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
        registry.append_environment_variable("CMAKE_PREFIX_PATH", builder.install_dir())
