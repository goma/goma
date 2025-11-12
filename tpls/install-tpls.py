#!/usr/bin/env python3
from tpl_tools.builder import Builder
from tpl_tools.registry import Registry
from tpl_tools.builder import mkdir_p
import importlib
import tpl_tools.utils as utils
import os
import sys
import argparse
import pathlib
import time


default_packages = [
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
    "lapack",
    "metis",
    "parmetis",
    "scotch",
    "arpack-ng",
    "scalapack",
    "mumps",
    "superlu_dist",
    "suitesparse",
    "trilinos",
    "hypre",
    "strumpack",
    "petsc",
    "petsc_complex",
    "sparse",
    "catch2",
]


if __name__ == "__main__":
    packages = default_packages.copy()
    CC = os.environ.get("CC")
    CXX = os.environ.get("CXX")
    FC = os.environ.get("FC")
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
        "--log-dir", help="Directory location for logs", type=pathlib.Path
    )
    parser.add_argument(
        "--build-shared",
        help="Build shared libraries (Default)",
        action="store_true",
        dest="enable_shared",
    )
    parser.add_argument(
        "--build-static",
        help="Build static libraries",
        dest="enable_shared",
        action="store_false",
    )
    parser.set_defaults(enable_shared=True)
    parser.add_argument(
        "-j", "--jobs", help="Number of parallel jobs", type=int, default=1
    )
    parser.add_argument(
        "--netlib-blas",
        help="Build using reference BLAS/LAPACK",
        action="store_true",
    )
    parser.add_argument(
        "--openblas",
        help="Build using OpenBLAS (Default)",
        dest="netlib_blas",
        action="store_true",
    )
    parser.set_defaults(netlib_blas=False)
    parser.add_argument(
        "--enable-parmetis",
        help="Build ParMETIS library, (Default, check license requirements)",
        action="store_true",
    )
    parser.add_argument(
        "--disable-parmetis",
        help="Disable ParMETIS library",
        dest="enable_parmetis",
        action="store_false",
    )
    parser.set_defaults(enable_parmetis=True)
    parser.add_argument(
        "INSTALL_DIR", help="Install location of TPLs", type=pathlib.Path
    )
    parser.add_argument(
        "--skip-ssl-verify",
        help="Disable SSL checks on download",
        dest="skip_ssl_verify",
        action="store_true",
    )
    parser.set_defaults(skip_ssl_verify=False)
    parser.add_argument(
        "--write-dynamic-library-path",
        help="Writes (DY)LD_LIBRARY_PATH to config, default is off",
        dest="write_dynamic_library_path",
        action="store_true",
    )
    parser.set_defaults(write_dynamic_library_path=False)
    parser.add_argument(
        "--install-complex-petsc",
        help="Install a complex version of petsc alongside the regular",
        dest="petsc_complex",
        action="store_true",
    )
    parser.set_defaults(petsc_complex=False)

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
    log_dir = os.path.join(install_dir, "logs")
    if args.log_dir:
        log_dir = os.path.abspath(os.path.expanduser(args.log_dir))
    mkdir_p(log_dir)
    logger = utils.PrintAndFileLogger(os.path.join(log_dir, "install-tpls.log"))
    if not args.petsc_complex:
        packages.remove("petsc_complex")

    if not args.enable_parmetis:
        logger.log("ParMETIS has been disabled ")
        logger.log("\tDisabling ParMETIS")
        logger.log(
            "\tDisabling SuperLU_DIST as Trilinos requires it be built with ParMETIS"
        )
        packages.remove("parmetis")
        packages.remove("superlu_dist")

    if args.netlib_blas:
        packages.remove("openblas")
    else:
        packages.remove("lapack")

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

    tpl_registry = Registry(logger)
    if CC:
        tpl_registry.set_environment_variable("CC", CC)
    else:
        logger.log("C compiler not set, defaulting to gcc, set with --cc")
        tpl_registry.set_environment_variable("CC", "gcc")

    if CXX:
        tpl_registry.set_environment_variable("CXX", CXX)
    else:
        logger.log("C++ compiler not set, defaulting to g++, set with --cxx")
        tpl_registry.set_environment_variable("CXX", "g++")

    if FC:
        tpl_registry.set_environment_variable("FC", FC)
        tpl_registry.set_environment_variable("F77", FC)
    else:
        logger.log("Fortran compiler not set, defaulting to gfortan, set with --fc")
        tpl_registry.set_environment_variable("FC", "gfortran")
        tpl_registry.set_environment_variable("F77", "gfortran")

    dependency_graph = {
        p: pkg.dependencies
        for p, pkg in [
            (
                p,
                importlib.import_module(
                    "tpl_tools." + ".".join(["packages", p])
                ).Package(),
            )
            for p in packages
        ]
    }

    if "openblas" in packages:
        for p in dependency_graph.keys():
            if "lapack" in dependency_graph[p]:
                dependency_graph[p].remove("lapack")
                dependency_graph[p].append("openblas")

    if not args.enable_parmetis:
        for p in dependency_graph.keys():
            for n in ["parmetis", "superlu_dist"]:
                if n in dependency_graph[p]:
                    dependency_graph[p].remove(n)

    logger.log("Dependency graph: {}".format(dependency_graph))

    install_order = utils.topological_sort(dependency_graph)[::-1]

    logger.log("Install order: {}".format(install_order))

    start_time = time.perf_counter()
    package_timings = {}

    for p in install_order:
        pm = importlib.import_module("tpl_tools." + ".".join(["packages", p]))
        pc = pm.Package()
        # default to 0 for skipped packages
        package_timings[p] = 0.0
        package_start_time = time.perf_counter()
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
                args.enable_shared,
                log_dir=log_dir,
                prebuilt=True,
                skip_ssl_verify=args.skip_ssl_verify,
            )
            if build.check(True):
                build.logger.log("Package {} found at {}".format(pc.name, package_dir))
            else:
                logger.log(
                    "Package {} not found check directory or let script build".format(
                        pc.name
                    ),
                    file=sys.stderr,
                )
                exit(1)
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
                args.enable_shared,
                log_dir=log_dir,
                skip_ssl_verify=args.skip_ssl_verify,
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
                logger.log(
                    "Package {} not built, contact Goma developers".format(pc.name),
                    file=sys.stderr,
                )
                exit(1)
            build.register()
        package_end_time = time.perf_counter()
        package_timings[p] = package_end_time - package_start_time
        logger.log(p, "took {:.2f} seconds".format(package_timings[p]))

    end_time = time.perf_counter()
    total_time = end_time - start_time

    logger.log("All packages built and installed to {}".format(install_dir))
    logger.log("Total time to build all packages: {:.2f} seconds".format(total_time))
    logger.log("Time to build individual packages:")
    for p in package_timings.keys():
        logger.log("\t{}: {:.2f} seconds".format(p, package_timings[p]))

    tpl_registry.config.write_config(
        os.path.join(install_dir, "config.sh"), "bash", args.write_dynamic_library_path
    )
    logger.log(
        "Bash config written to {}, source with bash".format(
            os.path.join(install_dir, "config.sh")
        )
    )

    tpl_registry.config.write_config(
        os.path.join(install_dir, "config.fish"),
        "fish",
        args.write_dynamic_library_path,
    )
    logger.log(
        "Fish config written to {}, source with fish".format(
            os.path.join(install_dir, "config.fish")
        )
    )
