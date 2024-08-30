from tpl_tools.packages import packages
from tpl_tools import utils
import os


class Package(packages.CMakePackage):
    def __init__(self):
        self.name = "catch2"
        self.version = "3.7.0"
        self.sha256 = "5b10cd536fa3818112a82820ce0787bd9f2a906c618429e7c4dea639983c8e88"
        self.filename = "catch2-" + self.version + ".tar.gz"
        self.url = (
            "https://github.com/catchorg/Catch2/archive/refs/tags/v"
            + self.version
            + ".tar.gz"
        )
        self.includes = ["catch2/catch_session.hpp"]

    def configure_options(self, builder):
        if builder.build_shared:
            builder.add_option("-DBUILD_SHARED_LIBS:BOOL=ON")
        else:
            builder.add_option("-DBUILD_SHARED_LIBS:BOOL=OFF")

    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        registry.set_environment_variable("CATCH2_DIR", builder.install_dir())
        registry.append_environment_variable("CMAKE_PREFIX_PATH", builder.install_dir())
