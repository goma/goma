from tpl_tools.packages import packages
import os


class Package(packages.AutotoolsPackage):
    def __init__(self):
        self.name = "petsc"
        self.version = "3.21.4"
        self.sha256 = "a9ae076d4617c7d84ce2bed37194022319c19f19b3930edf148b2bc8ecf2248d"
        self.filename = "petsc-" + self.version + ".tar.gz"
        self.url = (
            "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-"
            + self.version
            + ".tar.gz"
        )
        self.includes = ["petsc"]
        self.libraries = ["petsc", "HYPRE", "strumpack"]

    def set_environment(self, builder):
        builder.env = builder._registry.get_environment().copy()
        builder.env["PETSC_ARCH"] = "arch-c-opt"
        builder.env["PETSC_DIR"] = os.path.join(
            builder._extract_dir, builder._extracted_folder
        )

    def configure(self, builder):
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
        configure_options.append("--download-strumpack")
        configure_options.append("--with-scalapack=1")
        configure_options.append("--with-scalapack-dir=" + builder.env["SCALAPACK_DIR"])
        if "PARMETIS_DIR" in builder.env:
            configure_options.append("--with-superlu_dist=1")
            configure_options.append(
                "--with-superlu_dist-dir=" + builder.env["SUPERLU_DIST_DIR"]
            )
            configure_options.append("--with-parmetis=1")
            configure_options.append("--with-parmetis-dir=" + builder.env["PARMETIS_DIR"])
        configure_options.append("--with-metis=1")
        configure_options.append("--with-metis-dir=" + builder.env["METIS_DIR"])
        configure_options.append("--with-ptscotch=1")
        configure_options.append("--with-ptscotch-dir=" + builder.env["SCOTCH_DIR"])
        configure_options.append("--with-blas-lib=" + builder.env["BLAS_LIBRARIES"])
        configure_options.append("--with-lapack-lib=" + builder.env["LAPACK_LIBRARIES"])
        configure_options.append("--with-mumps=1")
        configure_options.append("--with-mumps-dir=" + builder.env["MUMPS_DIR"])
        configure_options.append("COPTFLAGS=-O3")
        configure_options.append("CXXOPTFLAGS=-O3")
        configure_options.append("FOPTFLAGS=-O3")
        builder.run_command(configure_options)

    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        registry.set_environment_variable("PETSC_DIR", builder.install_dir())
