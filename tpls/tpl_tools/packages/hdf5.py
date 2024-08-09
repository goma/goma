from tpl_tools.packages import packages
import os


class Package(packages.AutotoolsPackage):
    def __init__(self):
        self.name = "hdf5"
        self.version = "1.12.2"
        self.sha256 = "1a88bbe36213a2cea0c8397201a459643e7155c9dc91e062675b3fb07ee38afe"
        self.filename = "hdf5-" + self.version + ".tar.bz2"
        self.url = (
            "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-"
            + self.version
            + "/src/hdf5-"
            + self.version
            + ".tar.bz2"
        )
        self.libraries = ["hdf5"]

    def setDependencies(self, builder):
        builder.set_dependency("packages.openmpi")

    def set_environment(self, builder):
        builder.env = builder._registry.get_environment().copy()
        builder.env["CC"] = builder._registry.get_executable("mpicc")

    def configure_options(self, builder):
        if builder.build_shared:
            builder.add_option("--enable-shared")
        else:
            builder.add_option("--enable-shared=off")
        builder.add_option("--enable-parallel")
        builder.add_option("--with-default-api-version=v18")

    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        registry.set_environment_variable("HDF5_DIR", builder.install_dir())
        registry.append_environment_variable("CMAKE_PREFIX_PATH", builder.install_dir())
        registry.append_environment_variable(
            "PATH", os.path.join(builder.install_dir(), "bin")
        )
