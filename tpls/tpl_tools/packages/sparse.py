from tpl_tools.packages import packages
import os
import shutil


class Package(packages.GenericPackage):
    def __init__(self):
        self.name = "sparse"
        self.version = "1.4b"
        self.sha256 = "63e6646244fd8f4d89f7f70fbf4cfd46b7688d21b22840a0ce57d294a7496d28"
        self.filename = "sparse-" + self.version + ".tar.gz"
        self.url = (
            "http://downloads.sourceforge.net/project/sparse/sparse/sparse"
            + self.version
            + "/sparse"
            + self.version
            + ".tar.gz"
        )
        self.includes = ["spMatrix.h"]
        self.libraries = ["sparse"]

    def configure_options(self, builder):
        with open(
            os.path.join(
                builder._extract_dir, builder._extracted_folder, "src", "Makefile"
            ),
            "w",
        ) as f:
            if builder.build_shared:
                f.write("CFLAGS = -O2 -fPIC\n")
            else:
                f.write("CFLAGS = -O2\n")
            f.write("LINTFLAGS = -lc -lm\n")
            f.write("SHELL = /bin/sh\n")
            f.write("CC = " + builder.env["CC"] + "\n")
            f.write("\n")
            f.write("HFILES = spConfig.h spDefs.h spMatrix.h\n")
            f.write(
                "CFILES = spAllocate.c spBuild.c spFactor.c spOutput.c spSolve.c spUtils.c \\n"
            )
            f.write("         spFortran.c\n")
            f.write(
                "OFILES = spAllocate.o spBuild.o spFactor.o spOutput.o spSolve.o spUtils.o \\n"
            )
            f.write("         spFortran.o\n")
            if builder.build_shared:
                f.write("LIBRARY = ../lib/libsparse.so\n")
            else:
                f.write("LIBRARY = ../lib/libsparse.a\n")
            f.write("DESTINATION = ../bin/sparse\n")
            f.write("TESTC = spTest.c\n")
            f.write("TESTO = spTest.o\n")
            f.write("\n")
            f.write("SOURCE = $(HFILES) $(CFILES)\n")
            f.write("\n")
            f.write("$(DESTINATION)  : $(LIBRARY) $(TESTO)\n")
            f.write(
                "        $(CC) $(CFLAGS) -o $(DESTINATION) $(TESTO) $(LIBRARY) -lm\n"
            )
            f.write("\n")
            if builder.build_shared:
                f.write("$(LIBRARY)      : $(OFILES)\n")
                f.write("        $(CC) -shared $(LIBRARY) $?\n")
            else:
                f.write("$(LIBRARY)      : $(OFILES)\n")
                f.write("        ar r   $(LIBRARY) $?\n")
                f.write("        ranlib $(LIBRARY)\n")
            f.write("\n")

    def build(self, builder):
        builder.run_command(["make"])

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
        registry.set_environment_variable("SPARSE_DIR", builder.install_dir())
        registry.append_environment_variable("CMAKE_PREFIX_PATH", builder.install_dir())
