from tpl_tools.packages import packages
from tpl_tools import utils
import os


class Package(packages.CMakePackage):
    def __init__(self):
        self.name = "mumps"
        # Use a specific commit to avoid issues with NetCDF 4.9.3-rc1
        self.version = "5.7.3.1"
        self.sha256 = "54ac6470bf045c7ef5c9766f8a185078aff9953bfe3acbca23a8e3a2fab7d361"
        self.filename = "mumps-" + self.version + ".tar.gz"
        self.url = (
            "https://github.com/scivision/mumps/archive/refs/tags/v"
            + self.version
            + ".tar.gz"
        )
        self.libraries = [
            "cmumps",
            "dmumps",
            "zmumps",
            "mumps_common",
            "pord",
            "smumps",
        ]

    def setDependencies(self, builder):
        builder.set_dependency("packages.openmpi")
        return

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
        builder.add_option("-Dparmetis=yes")
        builder.add_option("-Dptscotch=yes")
        builder.add_option("-Dscotch=yes")
        # builder.add_option("-Dmetis=yes")
        builder.add_option("-DMETIS_ROOT=" + builder.env["METIS_DIR"])
        builder.add_option("-DBUILD_COMPLEX=ON")
        builder.add_option("-DBUILD_COMPLEX16=ON")
        builder.add_option(
            "-DMETIS_INCLUDE_DIR="
            + builder.env["METIS_DIR"]
            + "/include;"
            + builder.env["PARMETIS_DIR"]
            + "/include"
        )

    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        registry.set_environment_variable("MUMPS_DIR", builder.install_dir())
        registry.append_environment_variable("CMAKE_PREFIX_PATH", builder.install_dir())
