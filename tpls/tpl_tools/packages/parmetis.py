from tpl_tools.packages import packages


class Package(packages.CMakePackage):
    def __init__(self):
        self.name = "parmetis"
        self.version = "4.0.3-p9"
        self.sha256 = "612717e85992c984f09b0f5670be421bbb90a4c04145ab5b9a3358b92765d891"
        self.filename = "petsc-pkg-parmetis-" + self.version + ".tar.gz"
        self.url = (
            "https://bitbucket.org/petsc/pkg-parmetis/get/v" + self.version + ".tar.gz"
        )
        self.libraries = ["parmetis"]
        self.includes = ["parmetis.h"]

    def set_environment(self, builder):
        builder.env = builder._registry.get_environment().copy()
        builder.env["CC"] = builder._registry.get_executable("mpicc")
        builder.env["CXX"] = builder._registry.get_executable("mpicxx")
        builder.env["FC"] = builder._registry.get_executable("mpifort")

    def configure_options(self, builder):
        if builder.build_shared:
            builder.add_option("-DBUILD_SHARED_LIBS:BOOL=ON")
            builder.add_option("-DSHARED:BOOL=ON")
        else:
            builder.add_option("-DBUILD_SHARED_LIBS:BOOL=OFF")
            builder.add_option("-DSHARED:BOOL=OFF")
        builder.add_option("-DGKLIB_PATH=./headers")
        builder.add_option(
            "-DMETIS_PATH=" + builder._registry.environment["METIS_DIR"]
        )

    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        registry.set_environment_variable("PARMETIS_DIR", builder.install_dir())
        registry.append_environment_variable("CMAKE_PREFIX_PATH", builder.install_dir())
