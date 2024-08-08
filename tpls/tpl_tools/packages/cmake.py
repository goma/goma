from tpl_tools.packages import packages
import os


class Package(packages.AutotoolsPackage):
    def __init__(self):
        self.name = "cmake"
        self.version = "3.30.1"
        self.sha256 = "df9b3c53e3ce84c3c1b7c253e5ceff7d8d1f084ff0673d048f260e04ccb346e1"
        self.filename = "cmake-" + self.version + ".tar.gz"
        self.compression = "gz"
        self.url = (
            "https://github.com/Kitware/CMake/releases/download/v"
            + self.version
            + "/cmake-"
            + self.version
            + ".tar.gz"
        )
        self.executables = ["cmake"]

    def configure(self, builder):
        configure_options = ["./bootstrap"]
        configure_options.extend(builder.get_options())
        configure_options.append("--prefix=" + builder.install_dir())
        builder.run_command(configure_options, jobs_flag="--parallel=", parallel=True)

    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        registry.register_executable(
            os.path.join(builder.install_dir(), "bin", "cmake")
        )
        registry.append_environment_variable(
            "PATH", os.path.join(builder.install_dir(), "bin")
        )
