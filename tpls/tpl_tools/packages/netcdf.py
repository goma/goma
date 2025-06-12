from tpl_tools.packages import packages
import os


class Package(packages.CMakePackage):
    def __init__(self):
        self.name = "netcdf"
        self.version = "4.9.3"
        self.sha256 = "a474149844e6144566673facf097fea253dc843c37bc0a7d3de047dc8adda5dd"
        self.filename = "netcdf-c-" + self.version + ".tar.gz"
        self.url = (
            "https://downloads.unidata.ucar.edu/netcdf-c/"
            + self.version
            + "/netcdf-c-"
            + self.version
            + ".tar.gz"
        )
        self.executables = ["ncdump"]
        self.libraries = ["netcdf"]
        self.includes = ["netcdf.h"]
        self.dependencies = ["openmpi", "hdf5", "pnetcdf"]

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
        builder.add_option("-DNETCDF_ENABLE_DAP=OFF")
        builder.add_option("-DNETCDF_ENABLE_BYTERANGE:BOOL=OFF")
        builder.add_option("-DNETCDF_ENABLE_PNETCDF:BOOL=ON")
        builder.add_option("-DNETCDF_ENABLE_CDF5=ON")
        builder.add_option("-DNETCDF_ENABLE_MMAP:BOOL=ON")
        builder.add_option("-DHDF5_ROOT:PATH=" + builder.env["HDF5_DIR"])
        builder.add_option("-DHDF5_DIR:PATH=" + builder.env["HDF5_DIR"])
        builder.add_option("-DNETCDF_ENABLE_HDF5:BOOL=ON")
        builder.add_option("-DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON") 
        builder.add_option("-DBUILD_TESTING:BOOL=OFF") 
        builder.add_option("-DNETCDF_ENABLE_NCZARR_FILTERS:BOOL=OFF") 
        builder.add_option("-DNETCDF_ENABLE_NCZARR:BOOL=OFF")
        builder.add_option("-DNETCDF_ENABLE_V2_API:BOOL=OFF")
        builder.add_option("-DNETCDF_ENABLE_FILTER_TESTING:BOOL=OFF")
        builder.add_option("-DNETCDF_ENABLE_TESTS:BOOL=OFF")
        builder.add_option("-DNETCDF_ENABLE_QUANTIZE:BOOL=ON")

    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        registry.set_environment_variable("NETCDF_DIR", builder.install_dir())
        registry.prepend_environment_variable(
            "CMAKE_PREFIX_PATH", builder.install_dir()
        )
        registry.prepend_environment_variable(
            "PATH", os.path.join(builder.install_dir(), "bin")
        )
