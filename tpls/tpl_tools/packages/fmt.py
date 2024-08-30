from tpl_tools.packages import packages


class Package(packages.CMakePackage):
    def __init__(self):
        self.name = "fmt"
        self.version = "10.2.1"
        self.sha256 = "1250e4cc58bf06ee631567523f48848dc4596133e163f02615c97f78bab6c811"
        self.filename = "fmt-" + self.version + ".tar.gz"
        self.url = (
            "https://github.com/fmtlib/fmt/archive/refs/tags/"
            + self.version
            + ".tar.gz"
        )
        self.includes = ["fmt/core.h", "fmt/format.h"]
        self.libraries = ["fmt"]

    def setDependencies(self, builder):
        builder.set_dependency("packages.openmpi")
        builder.set_dependency("packages.hdf5")
        builder.set_dependency("packages.pnetcdf")
        builder.set_dependency("packages.netcdf")
        builder.set_dependency("packages.fmt")
        return

    def configure_options(self, builder):
        CXX = builder.env["CXX"]
        if builder.build_shared:
            builder.add_option("-DBUILD_SHARED_LIBS:BOOL=ON")
        else:
            builder.add_option("-DBUILD_SHARED_LIBS:BOOL=OFF")
        builder.add_option("-DCMAKE_CXX_COMPILER=" + CXX)

    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        registry.set_environment_variable("FMT_DIR", builder.install_dir())
        registry.append_environment_variable("CMAKE_PREFIX_PATH", builder.install_dir())
