from tpl_tools.packages import packages
from tpl_tools import utils


class Package(packages.AutotoolsPackage):
    def __init__(self):
        self.name = "pnetcdf"
        self.version = "1.14.0"
        self.sha256 = "575f189fb01c53f93b3d6ae0e506f46e19694807c81af0b9548e947995acf704"
        self.filename = "pnetcdf-" + self.version + ".tar.gz"
        self.url = (
            "https://parallel-netcdf.github.io/Release/pnetcdf-"
            + self.version
            + ".tar.gz"
        )
        self.libraries = ["pnetcdf"]
        self.dependencies = ["openmpi"]

    def set_environment(self, builder):
        builder.env = builder._registry.get_environment().copy()
        builder.env["CC"] = builder._registry.get_executable("mpicc")

        # check gcc version
        compiler, version = utils.check_gcc_clang_version(builder.env["CC"])
        extra_cflags = ""
        if compiler == "gcc" and version >= (15, 0, 0):
            extra_cflags = " -std=c11 -Wno-incompatible-pointer-types -Wno-implicit-function-declaration -Wno-int-conversion"
        builder.env["CFLAGS"] = "-O2" + extra_cflags
        builder.env["SEQ_CFLAGS"] = "-O2" + extra_cflags

    def configure_options(self, builder):
        if builder.build_shared:
            builder.add_option("--enable-shared")
        builder.add_option("--disable-fortran")
        builder.add_option("--disable-cxx")

    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        registry.set_environment_variable("PNETCDF_DIR", builder.install_dir())
        registry.prepend_environment_variable(
            "CMAKE_PREFIX_PATH", builder.install_dir()
        )
