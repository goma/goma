from tpl_tools.packages import packages


class Package(packages.CMakePackage):
    def __init__(self):
        self.name = "arpack-ng"
        self.version = "3.9.1"
        self.sha256 = "f6641deb07fa69165b7815de9008af3ea47eb39b2bb97521fbf74c97aba6e844"
        self.filename = "arpack-ng-" + self.version + ".tar.gz"
        self.url = (
            "https://github.com/opencollab/arpack-ng/archive/refs/tags/"
            + self.version
            + ".tar.gz"
        )
        self.libraries = ["arpack", "parpack"]

    def set_environment(self, builder):
        builder.env = builder._registry.get_environment().copy()
        builder.env["CC"] = builder._registry.get_executable("mpicc")
        builder.env["CXX"] = builder._registry.get_executable("mpicxx")
        builder.env["FC"] = builder._registry.get_executable("mpifort")

    def configure_options(self, builder):
        if builder.build_shared:
            builder.add_option("-DBUILD_SHARED_LIBS:BOOL=ON")
        else:
            builder.add_option("-DBUILD_SHARED_LIBS:BOOL=OFF")
        builder.add_option("-DMPI=ON")

    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        registry.set_environment_variable("ARPACK_DIR", builder.install_dir())
        registry.append_environment_variable("CMAKE_PREFIX_PATH", builder.install_dir())
