from tpl_tools.packages import packages
import os


class Package(packages.AutotoolsPackage):
    def __init__(self):
        self.name = "flex"
        self.version = "2.6.4"
        self.sha256 = "e87aae032bf07c26f85ac0ed3250998c37621d95f8bd748b31f15b33c45ee995"
        self.filename = "flex-" + self.version + "tar.bz2"
        self.url = (
            "https://github.com/westes/flex/releases/download/v" + self.version + "/flex-"
            + self.version
            + ".tar.gz"
        )
        self.executables = ["flex"]

    def setDependencies(self, builder):
        return

    def register(self, builder):
        registry = builder._registry
        registry.append_environment_variable("FLEX_DIR", builder.install_dir())
        registry.append_environment_variable("CMAKE_PREFIX_PATH", builder.install_dir())
        registry.append_environment_variable(
            "PATH", os.path.join(builder.install_dir(), "bin")
        )
