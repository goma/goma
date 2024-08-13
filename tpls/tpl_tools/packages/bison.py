from tpl_tools.packages import packages
import os


class Package(packages.AutotoolsPackage):
    def __init__(self):
        self.name = "bison"
        self.version = "3.8.2"
        self.sha256 = "06c9e13bdf7eb24d4ceb6b59205a4f67c2c7e7213119644430fe82fbd14a0abb"
        self.filename = "bison-" + self.version + "tar.bz2"
        self.url = (
            "http://ftp.gnu.org/gnu/bison/bison-"
            + self.version
            + ".tar.gz"
        )
        self.executables = ["bison"]

    def setDependencies(self, builder):
        return

    def register(self, builder):
        registry = builder._registry
        registry.append_environment_variable("BISON_DIR", builder.install_dir())
        registry.append_environment_variable("CMAKE_PREFIX_PATH", builder.install_dir())
        registry.append_environment_variable(
            "PATH", os.path.join(builder.install_dir(), "bin")
        )
