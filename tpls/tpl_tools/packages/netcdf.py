from tpl_tools.packages import packages
import os


class Package(packages.AutotoolsPackage):
    def __init__(self):
        self.name = "netcdf"
        self.version = "4.9.0"
        self.sha256 = "4c956022b79c08e5e14eee8df51b13c28e6121c2b7e7faadc21b375949400b49"
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

    def setDependencies(self, builder):
        builder.set_dependency("packages.openmpi")
        builder.set_dependency("packages.hdf5")
        builder.set_dependency("packages.pnetcdf")
        return

    def set_environment(self, builder):
        builder.env = builder._registry.get_environment().copy()
        builder.env["CC"] = builder._registry.get_executable("mpicc")
        builder.env["CXX"] = builder._registry.get_executable("mpicxx")
        builder.env["FC"] = builder._registry.get_executable("mpifort")
        builder.env["CFLAGS"] = "-I" + builder.env["HDF5_DIR"] + "/include"
        builder.env["CFLAGS"] += " -I" + builder.env["PNETCDF_DIR"] + "/include"
        builder.env["CPPFLAGS"] = "-I" + builder.env["HDF5_DIR"] + "/include"
        builder.env["CPPFLAGS"] += " -I" + builder.env["PNETCDF_DIR"] + "/include"
        builder.env["LDFLAGS"] = "-L" + builder.env["HDF5_DIR"] + "/lib"
        builder.env["LDFLAGS"] += " -L" + builder.env["PNETCDF_DIR"] + "/lib"

    def configure_options(self, builder):
        if builder.build_shared:
            builder.add_option("--enable-shared")
        else:
            builder.add_option("--enable-shared=off")
        builder.add_option("--enable-pnetcdf")
        builder.add_option("--disable-dap")

    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        registry.set_environment_variable("NETCDF_DIR", builder.install_dir())
        registry.append_environment_variable("CMAKE_PREFIX_PATH", builder.install_dir())
        registry.append_environment_variable(
            "PATH", os.path.join(builder.install_dir(), "bin")
        )
