from tpl_tools.packages import packages
import os


class Package(packages.GenericPackage):
    def __init__(self):
        self.name = "openblas"
        self.version = "0.3.28"
        self.sha256 = "f1003466ad074e9b0c8d421a204121100b0751c96fc6fcf3d1456bd12f8a00a1"
        self.filename = "OpenBLAS" + self.version + ".tar.gz"
        self.compression = "gz"
        self.url = (
            "https://github.com/OpenMathLib/OpenBLAS/releases/download/v"
            + self.version
            + "/OpenBLAS-"
            + self.version
            + ".tar.gz"
        )
        self.libraries = ["openblas"]

    def build(self, builder):
        self.build_command = ["make"]
        self.build_command.append("USE_THREAD=0")
        self.build_command.append("USE_OPENMP=0")
        if builder.build_shared:
            self.build_command.append("NO_SHARED=0")
        else:
            self.build_command.append("NO_SHARED=1")
        builder.run_command(self.build_command, jobs_flag="-j", parallel=True)

    def install(self, builder):
        command = self.build_command.copy()
        command.append("PREFIX=" + builder.install_dir())
        command.append("install")
        builder.run_command(command)

    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        ext = ".a"
        if builder.build_shared:
            ext = ".so"
        registry.set_environment_variable("OPENBLAS_DIR", builder.install_dir())
        registry.set_environment_variable(
            "BLAS_LIBRARIES",
            os.path.join(builder.install_dir(), "lib/libopenblas" + ext),
        )
        registry.set_environment_variable(
            "LAPACK_LIBRARIES",
            os.path.join(builder.install_dir(), "lib/libopenblas" + ext),
        )
        registry.append_environment_variable("CMAKE_PREFIX_PATH", builder.install_dir())
