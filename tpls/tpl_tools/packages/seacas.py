from tpl_tools.packages import packages
from tpl_tools import utils
import os


class Package(packages.CMakePackage):
    def __init__(self):
        self.name = "seacas"
        self.version = "v2025-06-07"
        self.sha256 = "2974705f2859e30bca48b619fda078bb771c0e94381af9e624749afb9fd72780"
        self.filename = "seacas-" + self.version + ".tar.gz"
        self.url = (
            "https://github.com/sandialabs/seacas/archive/" + self.version + ".tar.gz"
        )
        self.executables = ["algebra", "aprepro", "mapvar", "explore"]
        self.libraries = ["exodus", "aprepro_lib"]
        self.includes = ["exodusII.h", "aprepro.h"]
        self.dependencies = ["openmpi", "cmake", "fmt", "hdf5", "netcdf", "pnetcdf"]

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
            builder.add_option("-DSeacas_EXTRA_LINK_FLAGS=z;dl;m")
            builder.add_option("-DSEACASExodus_ENABLE_SHARED:BOOL=OFF")

        builder.add_option("-DTPL_ENABLE_Netcdf:BOOL=ON")
        builder.add_option("-DTPL_ENABLE_MPI:BOOL=ON")
        if utils.check_for_x11(builder._extract_dir, builder.env["CC"]):
            builder.add_option("-DTPL_ENABLE_X11:BOOL=ON")
        else:
            builder.add_option("-DTPL_ENABLE_X11:BOOL=OFF")

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

    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        registry.set_environment_variable("ACCESS", builder.install_dir())
        registry.set_environment_variable("SEACAS_DIR", builder.install_dir())
        registry.prepend_environment_variable(
            "CMAKE_PREFIX_PATH", builder.install_dir()
        )
        registry.prepend_environment_variable(
            "PATH", os.path.join(builder.install_dir(), "bin")
        )
        if builder.build_shared:
            registry.prepend_environment_variable(
                "PYTHONPATH", os.path.join(builder.install_dir(), "lib")
            )
