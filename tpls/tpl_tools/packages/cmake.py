from tpl_tools.packages import packages
import os


class Package(packages.AutotoolsPackage):
    def __init__(self):
        self.name = "cmake"
        self.version = "3.31.9"
        self.sha256 = "5d4fdec04247ca8a8e8f63692f0d0f1e9d6d082a2bdd008dff8ab3ba7215aa83"
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
        self.dependencies = []

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
        registry.prepend_environment_variable(
            "CMAKE_DIR", os.path.join(builder.install_dir())
        )
        registry.prepend_environment_variable(
            "PATH", os.path.join(builder.install_dir(), "bin")
        )
