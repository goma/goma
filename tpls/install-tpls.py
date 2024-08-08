#!/usr/bin/env python3
from tpl_tools.builder import Builder
from tpl_tools.registry import Registry
import importlib
import tpl_tools.utils as utils
import os
import argparse
import pathlib

packages = [
    "cmake",
    "openmpi",
    "hdf5",
    "pnetcdf",
    "netcdf",
    "fmt",
    "seacas",
    "openblas",
    "metis",
    "parmetis",
    "arpack-ng",
    "scalapack",
    "mumps",
    "superlu_dist",
    "suitesparse",
    "trilinos",
    "petsc",
    "omega_h",
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
    parser.add_argument("--cc", help="C compiler to use", type=pathlib.Path)
    parser.add_argument("--cxx", help="C++ compiler to use", type=pathlib.Path)
    parser.add_argument("--fc", help="Fortran compiler to use", type=pathlib.Path)
    parser.add_argument(
        "--download-dir", help="Download location of tarballs", type=pathlib.Path
    )
    parser.add_argument(
        "--extract-dir", help="Extract and Build location", type=pathlib.Path
    )
    parser.add_argument(
        "--build-shared", help="Build shared libraries", type=bool, default=True
    )
    parser.add_argument(
        "-j", "--jobs", help="Number of parallel jobs", type=int, default=1
    )
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
                "CC {} is environment variable, overriding with --cc={}".format(
                    CC, args.cc
                )
            )
        CC = str(args.cc)
    if args.cxx:
        CXX = str(args.cxx)
        logger.log(
            "CXX {} is environment variable, overriding with --cxx={}".format(
                CXX, args.cxx
            )
        )
    if args.fc:
        FC = str(args.fc)
        logger.log(
            "FC {} is environment variable, overriding with --fc={}".format(FC, args.fc)
        )

    if CC:
        tpl_registry.set_environment_variable("CC", CC)
    if CXX:
        tpl_registry.set_environment_variable("CXX", CXX)
    if FC:
        tpl_registry.set_environment_variable("FC", FC)
        tpl_registry.set_environment_variable("F77", FC)

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
