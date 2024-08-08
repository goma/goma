from tpl_tools.packages import packages


class Package(packages.CMakePackage):
    def __init__(self):
        self.name = "metis"
        self.version = "5.1.0-p12"
        self.sha256 = "2f16abe87394d4fd1430d66a8995076c352af40f7a4d1c5300b3b1cc9d545663"
        self.filename = "petsc-pkg-metis-" + self.version + ".tar.gz"
        self.url = (
            "https://bitbucket.org/petsc/pkg-metis/get/v" + self.version + ".tar.gz"
        )
        self.libraries = ["metis"]
        self.includes = ["metis.h"]

    def configure_options(self, builder):
        if builder.build_shared:
            builder.add_option("-D=BUILD_SHARED_LIBS:BOOL=ON")
            builder.add_option("-D=SHARED:BOOL=ON")
        else:
            builder.add_option("-D=BUILD_SHARED_LIBS:BOOL=OFF")
            builder.add_option("-D=SHARED:BOOL=OFF")
        builder.add_option("-D=GKLIB_PATH=./GKlib")

    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        registry.set_environment_variable("METIS_DIR", builder.install_dir())
        registry.append_environment_variable("CMAKE_PREFIX_PATH", builder.install_dir())
