from tpl_tools.packages import packages
from tpl_tools import utils
import os


class Package(packages.AutotoolsPackage):
    def __init__(self):
        self.name = "petsc-complex"
        self.version = "3.23.0"
        self.sha256 = "aeebd7094f4d583fd04700e73779caa7d9a3d54742e95eff2c3dd87768a79063"
        self.filename = "petsc-" + self.version + ".tar.gz"
        self.url = (
            "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-"
            + self.version
            + ".tar.gz"
        )
        self.includes = ["petsc"]
        self.libraries = ["petsc", "HYPRE"]
        self.dependencies = ["openmpi", "scalapack", "metis", "mumps", "scotch"]

    def set_environment(self, builder):
        builder.env = builder._registry.get_environment().copy()
        builder.env["PETSC_ARCH"] = "arch-c-complex"
        builder.env["PETSC_DIR"] = os.path.join(
            builder._extract_dir, builder._extracted_folder
        )

    def petsc_options(self, builder):
        configure_options = ["./configure"]
        configure_options.extend(builder.get_options())
        configure_options.append("--prefix=" + builder.install_dir())
        configure_options.append("--with-mpi")
        configure_options.append(
            "--with-cc=" + builder._registry.get_executable("mpicc")
        )
        configure_options.append(
            "--with-fc=" + builder._registry.get_executable("mpifort")
        )
        configure_options.append(
            "--with-cxx=" + builder._registry.get_executable("mpicxx")
        )
        if builder.build_shared:
            configure_options.append("--with-shared-libraries=1")
        else:
            configure_options.append("--with-shared-libraries=0")
        configure_options.append("--with-debugging=0")
        configure_options.append("--download-hypre")
        configure_options.append("--with-scalapack=1")
        configure_options.append("--with-scalapack-dir=" + builder.env["SCALAPACK_DIR"])
        configure_options.append("--with-metis=1")
        configure_options.append("--with-metis-dir=" + builder.env["METIS_DIR"])
        if "PARMETIS_DIR" in builder.env:
            configure_options.append("--with-parmetis=1")
            configure_options.append(
                "--with-parmetis-dir=" + builder.env["PARMETIS_DIR"]
            )
        configure_options.append("--with-ptscotch=1")
        configure_options.append("--with-ptscotch-dir=" + builder.env["SCOTCH_DIR"])
        configure_options.append("--with-blas-lib=" + builder.env["BLAS_LIBRARIES"])
        configure_options.append("--with-lapack-lib=" + builder.env["LAPACK_LIBRARIES"])
        configure_options.append("--with-mumps=1")
        configure_options.append("--with-mumps-dir=" + builder.env["MUMPS_DIR"])
        configure_options.append("--with-scalar-type=complex")
        compiler, version = utils.check_gcc_clang_version(builder.env["CC"])
        if compiler == "gcc" and version >= (14, 0, 0):
            configure_options.append("COPTFLAGS=-O3 -Wno-incompatible-pointer-types")
        else:
            configure_options.append("COPTFLAGS=-O3")

        configure_options.append("CXXOPTFLAGS=-O3")
        configure_options.append("FOPTFLAGS=-O3")
        configure_options.append("PETSC_ARCH=arch-c-complex")
        return configure_options

    def configure(self, builder):
        configure_options = self.petsc_options(builder)
        builder.run_command(configure_options)

    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        registry.set_environment_variable("PETSC_COMPLEX_DIR", builder.install_dir())
