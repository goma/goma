from tpl_tools.packages import packages
from tpl_tools import utils
import os
import shutil
import sys


class Package(packages.GenericPackage):
    def __init__(self):
        self.name = "mumps"
        self.version = "5.8.0"
        self.sha256 = "d762eb8b1d9843a0993b8cfc137d043d04c7c51877ad37c94560433a474340a0"
        self.filename = "mumps-" + self.version + ".tar.gz"
        self.url = "https://mumps-solver.org/MUMPS_" + self.version + ".tar.gz"
        self.libraries = ["dmumps", "zmumps"]
        self.dependencies = [
            "lapack",
            "scalapack",
            "scotch",
            "metis",
            "parmetis",
            "openmpi",
        ]

    def set_environment(self, builder):
        builder.env = builder._registry.get_environment().copy()
        builder.env["CC"] = builder._registry.get_executable("mpicc")
        builder.env["CXX"] = builder._registry.get_executable("mpicxx")
        builder.env["FC"] = builder._registry.get_executable("mpifort")

    def configure_options(self, builder):
        with open(
            os.path.join(
                builder._extract_dir, builder._extracted_folder, "Makefile.inc"
            ),
            "w",
        ) as f:
            f.write("LPORDDIR = $(topdir)/PORD/lib/\n")
            f.write("IPORD    = -I$(topdir)/PORD/include/\n")
            f.write("LPORD    = -L$(LPORDDIR) -lpord$(PLAT)\n")
            f.write("LMETISDIR = " + builder.env["METIS_DIR"] + "/lib\n")
            if "PARMETIS_DIR" in builder.env:
                f.write("LPARMETISDIR = " + builder.env["PARMETIS_DIR"] + "/lib\n")
            f.write("IMETIS = -I" + builder.env["METIS_DIR"] + "/include")

            if "PARMETIS_DIR" in builder.env:
                f.write(" -I" + builder.env["PARMETIS_DIR"] + "/include \n")
            else:
                f.write("\n")
            f.write("LMETIS    = -L$(LPARMETISDIR) -lparmetis -L$(LMETISDIR) -lmetis\n")
            f.write("SCOTCHDIR = " + builder.env["SCOTCH_DIR"] + "\n")
            f.write("ISCOTCH = -I" + builder.env["SCOTCH_DIR"] + "/include\n")
            f.write(
                "LSCOTCH = -L"
                + builder.env["SCOTCH_DIR"]
                + "/lib -lptesmumps -lptscotch -lptscotcherr -lscotch -lscotcherr\n"
            )
            f.write("ORDERINGSF  = -Dpord -Dmetis -Dptscotch")
            if "PARMETIS_DIR" in builder.env:
                f.write(" -Dparmetis\n")
            else:
                f.write("\n")
            f.write("ORDERINGSC  = $(ORDERINGSF)\n")
            f.write("LORDERINGS = $(LMETIS) $(LPORD) $(LSCOTCH)\n")
            f.write("IORDERINGSF = $(ISCOTCH)\n")
            f.write("IORDERINGSC = $(IMETIS) $(IPORD) $(ISCOTCH)\n")
            ext = utils.get_library_extension(True)
            f.write("LIBEXT_SHARED  = " + ext + "\n")
            if sys.platform == "darwin":
                f.write("SONAME = -install_name\n")
            else:
                f.write("SONAME = -soname\n")
            f.write(
                "SHARED_OPT = -shared -Wl,-rpath,"
                + builder.install_dir()
                + "/lib -Wl,-rpath,"
                + builder.env["SCOTCH_DIR"]
                + "/lib"
            )
            if builder.env["PARMETIS_DIR"]:
                f.write(" -Wl,-rpath," + builder.env["PARMETIS_DIR"] + "/lib")
            f.write(" -Wl,-rpath," + builder.env["METIS_DIR"] + "/lib\n")
            f.write("FPIC_OPT = -fPIC\n")
            f.write("LIBEXT  = .a\n")
            f.write("OUTC    = -o\n")
            f.write("OUTF    = -o\n")
            f.write("RM      = rm -f\n")
            f.write("# CC : C compiler\n")
            f.write("CC      = " + builder.env["CC"] + "\n")
            f.write("# FC : Fortran 90 compiler\n")
            f.write("FC      = " + builder.env["FC"] + "\n")
            f.write("FL      = " + builder.env["FC"] + "\n")
            f.write("AR      = ar vr \n")
            f.write("RANLIB  = ranlib \n")
            f.write(
                "LAPACK = "
                + builder.env["LAPACK_LIBRARIES"]
                + " "
                + builder.env["BLAS_LIBRARIES"]
                + "\n"
            )
            f.write(
                "SCALAP  = -L"
                + builder.env["SCALAPACK_DIR"]
                + "/lib -lscalapack "
                + builder.env["LAPACK_LIBRARIES"]
                + " "
                + builder.env["BLAS_LIBRARIES"]
                + "\n"
            )
            f.write("BLAS = " + builder.env["BLAS_LIBRARIES"] + "\n")
            f.write("LIBOTHERS = -lpthread\n")
            f.write("CDEFS = -DAdd_\n")
            f.write("FPIC_OPT = -fPIC\n")
            f.write("OPTF    = -O3 -fallow-argument-mismatch\n")
            f.write("OPTC    = -O3 -fPIC -I.\n")
            f.write("OPTL    = -O3 -fallow-argument-mismatch\n")
            f.write("INCPAR  = \n")
            f.write("LIBPAR  = $(SCALAP) $(LAPACK) -lmpi\n")
            f.write("INCS = $(INCPAR)\n")
            f.write("LIBS = $(LIBPAR)\n")
            f.write("LIBSEQNEEDED =\n")

    def build(self, builder):
        if builder.build_shared:
            builder.run_command(["make", "allshared"], jobs_flag="-j", parallel=True)
        else:
            builder.run_command(["make", "all"], jobs_flag="-j", parallel=True)

    def install(self, builder):
        if not os.path.exists(builder.install_dir()):
            os.makedirs(builder.install_dir())
        build_dir = os.path.join(builder._extract_dir, builder._extracted_folder)
        shutil.copytree(
            os.path.join(build_dir, "include"),
            os.path.join(builder.install_dir(), "include"),
        )
        shutil.copytree(
            os.path.join(build_dir, "lib"), os.path.join(builder.install_dir(), "lib")
        )

    def register(self, builder):
        registry = builder._registry
        registry.register_package(self.name, builder.install_dir())
        registry.set_environment_variable("MUMPS_DIR", builder.install_dir())
        registry.prepend_environment_variable(
            "CMAKE_PREFIX_PATH", builder.install_dir()
        )
