from tpl_tools.packages import packages


class Package(packages.AutotoolsPackage):
    def __init__(self):
        self.name = "pnetcdf"
        self.version = "1.13.0"
        self.sha256 = "aba0f1c77a51990ba359d0f6388569ff77e530ee574e40592a1e206ed9b2c491"
        self.filename = "pnetcdf-" + self.version + ".tar.gz"
        self.url = (
            "https://parallel-netcdf.github.io/Release/pnetcdf-"
            + self.version
            + ".tar.gz"
        )
        self.libraries = ["pnetcdf"]

    def setDependencies(self, builder):
        builder.set_dependency("packages.openmpi")
        return

    def set_environment(self, builder):
        builder.env = builder._registry.get_environment().copy()
        builder.env["CC"] = builder._registry.get_executable("mpicc")

    def configure_options(self, builder):
        if builder.build_shared:
            builder.add_option("--enable-shared")
        builder.add_option("--disable-fortran")
        builder.add_option("--disable-cxx")

    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        registry.set_environment_variable("PNETCDF_DIR", builder.install_dir())
        registry.append_environment_variable("CMAKE_PREFIX_PATH", builder.install_dir())
