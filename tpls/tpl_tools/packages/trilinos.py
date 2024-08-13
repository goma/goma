from tpl_tools.packages import packages


class Package(packages.CMakePackage):
    def __init__(self):
        self.name = "trilinos"
        self.version = "16.0.0"
        self.sha256 = "46bfc40419ed2aa2db38c144fb8e61d4aa8170eaa654a88d833ba6b92903f309"
        self.filename = "trilinos-" + self.version + ".tar.gz"
        self.url = (
            "https://github.com/trilinos/Trilinos/archive/refs/tags/trilinos-release-"
            + self.version.replace(".", "-")
            + ".tar.gz"
        )
        self.libraries = ["amesos2", "belos", "aztecoo", "amesos"]
        self.includes = [
            "Amesos2.hpp",
            "BelosSolverFactory.hpp",
            "Sacado.hpp",
            "AztecOO.h",
            "Amesos.h",
        ]

    def setDependencies(self, builder):
        return

    def set_environment(self, builder):
        builder.env = builder._registry.get_environment().copy()
        builder.env["CC"] = builder._registry.get_executable("mpicc")
        builder.env["CXX"] = builder._registry.get_executable("mpicxx")
        builder.env["FC"] = builder._registry.get_executable("mpifort")

    def configure_options(self, builder):
        if builder.build_shared:
            builder.add_option("-D=BUILD_SHARED_LIBS:BOOL=ON")
        else:
            builder.add_option("-D=BUILD_SHARED_LIBS:BOOL=OFF")
        CC = builder.env["CC"]
        CXX = builder.env["CXX"]
        FC = builder.env["FC"]
        builder.add_option("-DCMAKE_C_COMPILER=" + CC)
        builder.add_option("-DCMAKE_CXX_COMPILER=" + CXX)
        builder.add_option("-DCMAKE_Fortran_COMPILER=" + FC)
        builder.add_option("-DTrilinos_SHOW_DEPRECATED_WARNINGS=OFF")
        builder.add_option("-DCMAKE_VERBOSE_MAKEFILE:BOOL=TRUE")
        builder.add_option("-DTPL_ENABLE_Boost:BOOL=OFF")
        builder.add_option("-DTrilinos_SHOW_DEPRECATED_WARNINGS:BOOL=OFF")
        builder.add_option("-DTrilinos_ENABLE_SECONDARY_TESTED_CODE=ON")
        builder.add_option("-DTrilinos_ENABLE_Triutils:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_SEACAS:BOOL=OFF")
        builder.add_option("-DTrilinos_ENABLE_Epetra:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_Xpetra:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_Ifpack:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_Ifpack2:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_Teuchos:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_ML:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_MueLu:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_Stratimikos:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_Teko:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_Belos:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_Amesos2:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_Amesos:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_AztecOO:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_Sacado:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_EpetraExt:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_Thyra:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_ThyraTpetraAdapters:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_ThyraEpetraAdapters:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_Tpetra:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_Stratimikos:BOOL=ON")
        builder.add_option("-DTrilinos_ENABLE_TESTS:BOOL=OFF")
        builder.add_option("-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON")
        builder.add_option("-DTPL_ENABLE_MPI:BOOL=ON ")
        builder.add_option("-DMPI_BASE_DIR:PATH=" + builder.env["MPI_HOME"])
        builder.add_option("-DEpetraExt_BUILD_GRAPH_REORDERINGS:BOOL=ON")
        builder.add_option("-DTPL_ENABLE_LAPACK:BOOL=ON")
        builder.add_option("-DLAPACK_LIBRARY_DIRS=" + builder.env["OPENBLAS_DIR"])
        builder.add_option("-DLAPACK_LIBRARY_NAMES=openblas")
        builder.add_option("-DTPL_ENABLE_BLAS:BOOL=ON ")
        builder.add_option("-DHAVE_EPETRA_LAPACK_GSSVD3:BOOL=ON ")
        builder.add_option("-DBLAS_LIBRARY_DIRS=" + builder.env["OPENBLAS_DIR"])
        builder.add_option("-DBLAS_LIBRARY_NAMES=openblas")
        builder.add_option("-DTPL_ENABLE_UMFPACK:BOOL=ON ")
        builder.add_option(
            "-DUMFPACK_LIBRARY_NAMES:STRING=umfpack;amd;suitesparseconfig;cholmod;colamd;ccolamd;camd"
        )
        builder.add_option(
            "-DUMFPACK_LIBRARY_DIRS:PATH=" + builder.env["SUITESPARSE_DIR"] + "/lib"
        )
        builder.add_option(
            "-DUMFPACK_INCLUDE_DIRS:PATH="
            + builder.env["SUITESPARSE_DIR"]
            + "/include/suitesparse"
        )
        builder.add_option("-DTPL_ENABLE_AMD:BOOL=ON ")
        builder.add_option("-DAMD_LIBRARY_NAMES:STRING=amd;suitesparseconfig")
        builder.add_option(
            "-DAMD_LIBRARY_DIRS:PATH=" + builder.env["SUITESPARSE_DIR"] + "/lib"
        )
        builder.add_option(
            "-DAMD_INCLUDE_DIRS:PATH="
            + builder.env["SUITESPARSE_DIR"]
            + "/include/suitesparse"
        )
        if "SUPERLU_DIST_DIR" in builder.env:
            builder.add_option("-D Amesos_ENABLE_SuperLUDist:BOOL=ON ")
            builder.add_option("-DTPL_ENABLE_SuperLUDist:BOOL=ON ")
            builder.add_option("-DSuperLUDist_LIBRARY_NAMES:STRING=superlu_dist")
            builder.add_option(
                "-DSuperLUDist_LIBRARY_DIRS:PATH="
                + builder.env["SUPERLU_DIST_DIR"]
                + "/lib"
            )
            builder.add_option(
                "-DSuperLUDist_INCLUDE_DIRS:PATH="
                + builder.env["SUPERLU_DIST_DIR"]
                + "/include"
            )
        ext = ".a"
        if builder.build_shared:
            ext = ".so"
        if "PARMETIS_DIR" in builder.env:
            builder.add_option("-DTPL_ENABLE_ParMETIS:BOOL=ON ")
            builder.add_option("-D Amesos_ENABLE_ParMETIS:BOOL=ON ")
            builder.add_option(
                "-DParMETIS_LIBRARY_DIRS:PATH="
                + builder.env["PARMETIS_DIR"]
                + "/lib;"
                + builder.env["METIS_DIR"]
                + "/lib"
            )
            builder.add_option(
                "-DTPL_ParMETIS_INCLUDE_DIRS:PATH="
                + builder.env["PARMETIS_DIR"]
                + "/include;"
                + builder.env["METIS_DIR"]
                + "/include"
            )
            builder.add_option(
               "-DTPL_ParMETIS_LIBRARIES="
               + builder.env["PARMETIS_DIR"]
               + "/lib/libparmetis"
               + ext
               + ";"
               + builder.env["METIS_DIR"]
               + "/lib/libmetis"
               + ext
            )
        else:
            builder.add_option("-DTPL_ENABLE_ParMETIS:BOOL=OFF")
        builder.add_option("-DTPL_ENABLE_MUMPS:BOOL=ON ")
        builder.add_option(
            "-DMUMPS_LIBRARY_NAMES:STRING=dmumps;mumps_common;pord;scalapack;ptesmumps;ptscotch;ptscotcherr;scotch;scotcherr"
        )
        builder.add_option(
            "-DMUMPS_LIBRARY_DIRS:PATH="
            + builder.env["MUMPS_DIR"]
            + "/lib;"
            + builder.env["SCALAPACK_DIR"]
            + "/lib;"
            + builder.env["SCOTCH_DIR"]
            + "/lib"
        )
        builder.add_option(
            "-DMUMPS_INCLUDE_DIRS:PATH=" + builder.env["MUMPS_DIR"] + "/include"
        )
        builder.add_option("-DCMAKE_CXX_FLAGS:STRING=-DMUMPS_5_0")
        builder.add_option("-DAmesos_ENABLE_SCALAPACK:BOOL=ON ")
        builder.add_option(
            "-DSCALAPACK_LIBRARY_DIRS:FILEPATH=" + builder.env["SCALAPACK_DIR"] + "/lib"
        )
        builder.add_option("-DSCALAPACK_LIBRARY_NAMES:STRING=scalapack")
        builder.add_option("-D Amesos_ENABLE_LAPACK:BOOL=ON ")
        builder.add_option("-D Amesos_ENABLE_KLU:BOOL=ON ")
        builder.add_option("-D Amesos_ENABLE_UMFPACK:BOOL=ON ")
        builder.add_option("-D Amesos_ENABLE_MUMPS:BOOL=ON ")
        builder.add_option("-D Tpetra_INST_INT_INT:BOOL=ON ")

    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        registry.set_environment_variable("TRILINOS_DIR", builder.install_dir())
        registry.append_environment_variable("CMAKE_PREFIX_PATH", builder.install_dir())
