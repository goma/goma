#!/usr/bin/env python3
from tpl_tools.builder import Builder
from tpl_tools.registry import Registry
import importlib
import tpl_tools.utils as utils
import os
import argparse
import pathlib


# Currently we don't check dependencies so this needs to be in a specific order
packages = [
    "cmake",
    "openmpi",
    "hdf5",
    "pnetcdf",
    "netcdf",
    "omega_h",
    "fmt",
    "seacas",
    "bison",
    "flex",
    "openblas",
    "metis",
    "parmetis",
    "scotch",
    "arpack-ng",
    "scalapack",
    "mumps",
    "superlu_dist",
    "suitesparse",
    "trilinos",
    "petsc",
    "petsc_complex",
    "sparse",
]


if __name__ == "__main__":
    CC = os.environ.get("CC")
    CXX = os.environ.get("CXX")
    FC = os.environ.get("FC")
    tpl_registry = Registry()
    logger = utils.PrintLogger()
    parser = argparse.ArgumentParser(
        description="""Third party library installer for the finite element code Goma"""
    )
    parser.add_argument(
        "--cc",
        help="C compiler to use, default is CC environment variable or gcc",
        type=pathlib.Path,
    )
    parser.add_argument(
        "--cxx",
        help="C++ compiler to use, default is CXX environment variable or g++",
        type=pathlib.Path,
    )
    parser.add_argument(
        "--fc",
        help="Fortran compiler to use, defualt is FC environment variable or gfortran",
        type=pathlib.Path,
    )
    parser.add_argument(
        "--download-dir", help="Download location of tarballs", type=pathlib.Path
    )
    parser.add_argument(
        "--extract-dir", help="Extract and Build location", type=pathlib.Path
    )
    parser.add_argument(
        "--build-shared", help="Build shared libraries (Default)", action='store_true'
    )
    parser.add_argument(
        "--build-static", help="Build static libraries", dest='enable_shared', action='store_false'
    )
    parser.set_defaults(build_shared=True)
    parser.add_argument(
        "-j", "--jobs", help="Number of parallel jobs", type=int, default=1
    )
    parser.add_argument(
        "--enable-parmetis", help="Build ParMETIS library, (Default, check license requirements)", action='store_true'
    )
    parser.add_argument(
        "--disable-parmetis", help="Disable ParMETIS library", dest='enable_parmetis', action='store_false'
    )
    parser.set_defaults(enable_parmetis=True)
    parser.add_argument(
        "INSTALL_DIR", help="Install location of TPLs", type=pathlib.Path
    )


    for p in packages:
        pm = importlib.import_module("tpl_tools." + ".".join(["packages", p]))
        pc = pm.Package()
        parser.add_argument(
            "--" + pc.name.replace("-", "_") + "-dir",
            help="System location of package {}".format(pc.name),
            type=pathlib.Path,
        )
    args = parser.parse_args()
    if not args.enable_parmetis:
        print("ParMETIS has been disabled ")
        print("\tDisabling ParMETIS")
        print("\tDisabling SuperLU_DIST as Trilinos requires it be built with ParMETIS")
        packages.remove("parmetis")
        packages.remove("superlu_dist")



    install_dir = os.path.abspath(os.path.expanduser(args.INSTALL_DIR))
    download_dir = os.path.join(install_dir, "downloads")
    if args.download_dir:
        download_dir = os.path.abspath(os.path.expanduser(args.download_dir))
    extract_dir = os.path.join(install_dir, "sources")
    if args.extract_dir:
        extract_dir = os.path.abspath(os.path.expanduser(args.extract_dir))

    jobs = args.jobs



    if args.cc:
        if CC:
            logger.log(
                "CC {} is an environment variable, overriding with --cc={}".format(
                    CC, args.cc
                )
            )
        CC = str(args.cc)
    if args.cxx:
        if CXX:
            logger.log(
                "CXX {} is an environment variable, overriding with --cxx={}".format(
                    CXX, args.cxx
                )
            )
        CXX = str(args.cxx)
    if args.fc:
        if FC:
            logger.log(
                "FC {} is an environment variable, overriding with --fc={}".format(
                    FC, args.fc
                )
            )
        FC = str(args.fc)

    if CC:
        tpl_registry.set_environment_variable("CC", CC)
    else:
        print("C compiler not set, defaulting to gcc, set with --cc")
        tpl_registry.set_environment_variable("CC", "gcc")

    if CXX:
        tpl_registry.set_environment_variable("CXX", CXX)
    else:
        print("C++ compiler not set, defaulting to g++, set with --cxx")
        tpl_registry.set_environment_variable("CXX", "g++")

    if FC:
        tpl_registry.set_environment_variable("FC", FC)
        tpl_registry.set_environment_variable("F77", FC)
    else:
        print("Fortran compiler not set, defaulting to gfortan, set with --fc")
        tpl_registry.set_environment_variable("FC", "gfortran")
        tpl_registry.set_environment_variable("F77", "gfortran")

    for p in packages:
        pm = importlib.import_module("tpl_tools." + ".".join(["packages", p]))
        pc = pm.Package()
        if getattr(args, pc.name.replace("-", "_") + "_dir"):
            package_dir = getattr(args, pc.name + "_dir")

            build = Builder(
                pc,
                jobs,
                download_dir,
                extract_dir,
                package_dir,
                logger,
                tpl_registry,
                args.build_shared,
            )
            if build.check(True):
                build.logger.log("Package {} found at {}".format(pc.name, package_dir))
            else:
                break
            build.register()
        else:
            build = Builder(
                pc,
                jobs,
                download_dir,
                extract_dir,
                install_dir,
                logger,
                tpl_registry,
                args.build_shared,
            )

            if build.check():
                build.logger.log(
                    "Package {} already built, skipping".format(build._package.name)
                )
                build.register()
                continue
            build.download()
            build.extract()
            build.build()
            build.install()
            if not build.check():
                break
            build.register()
    tpl_registry.config.write_config(os.path.join(install_dir, "config.sh"))
    logger.log(
        "Bash config written to {}, source with bash".format(
            os.path.join(install_dir, "config.sh")
        )
    )

    tpl_registry.config.write_config(
        os.path.join(install_dir, "config.fish"), shell="fish"
    )
    logger.log(
        "Fish config written to {}, source with fish".format(
            os.path.join(install_dir, "config.fish")
        )
    )
