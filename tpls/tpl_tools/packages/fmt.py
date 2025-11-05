from tpl_tools.packages import packages


class Package(packages.CMakePackage):
    def __init__(self):
        self.name = "fmt"
        self.version = "11.2.0"
        self.sha256 = "bc23066d87ab3168f27cef3e97d545fa63314f5c79df5ea444d41d56f962c6af"
        self.filename = "fmt-" + self.version + ".tar.gz"
        self.url = (
            "https://github.com/fmtlib/fmt/archive/refs/tags/"
            + self.version
            + ".tar.gz"
        )
        self.includes = ["fmt/core.h", "fmt/format.h"]
        self.libraries = ["fmt"]
        self.dependencies = ["cmake"]

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
        registry.prepend_environment_variable(
            "CMAKE_PREFIX_PATH", builder.install_dir()
        )
