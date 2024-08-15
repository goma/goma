from tpl_tools.packages import packages
import os


class Package(packages.AutotoolsPackage):
    def __init__(self):
        self.name = "hdf5"
        self.version = "1.14.4.3"
        version_tag = "1.14.4-3"
        self.sha256 = "019ac451d9e1cf89c0482ba2a06f07a46166caf23f60fea5ef3c37724a318e03"
        self.filename = "hdf5-" + self.version + ".tar.bz2"
        self.url = (
            "https://github.com/HDFGroup/hdf5/releases/download/hdf5_"
            + self.version
            + "/hdf5-"
            + version_tag
            + ".tar.gz"
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
