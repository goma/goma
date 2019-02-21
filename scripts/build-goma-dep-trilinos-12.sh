#!/usr/bin/env bash
#
# This software is distributed under the GNU General Public License.
#
# Goma dependency build script
#
# 2014, July: Created by Cory Miner based off of notes by Scott Roberts
# 2014: Weston Ortiz modified to include blacs, scalapack and mumps support in trilinos
# 2015, February: Andrew Cochrane modified to handle recent changes to SEACAS-latest and hdf5, also added some logic to help with troubleshooting
# 2017, June: Andrew Cochrane added compatibility for Intel Parallel Studio.
# 2017, July: Eli Stone updated many libraries and enabled standalone openmpi - modularized the script a little bit..
# 2017, Nov: Andrew added "user" CC_NAME option

# This script tries to build all of the libraries needed for goma
# in the specified directories.
#
# The libraries built should be suitable for a minimal settings.mk

# Set script to exit on error
# set -e

#These can be set here by uncommenting, or in environment variables beforehand (such as by running the Easy Goma Builder script), or by running the script as normal

# Tells the TPL builder which C compiler to use, choices are [intel, gnu, user]
#CC_NAME="user"

# if CC_NAME="user" the following variables must be specified manually

SYSTEM_COMPILER_BASE=/usr
SYSTEM_CC=$SYSTEM_COMPILER_BASE/bin/gcc
SYSTEM_CXX=$SYSTEM_COMPILER_BASE/bin/g++
SYSTEM_FC=$SYSTEM_COMPILER_BASE/bin/gfortran
SYSTEM_CPP=$SYSTEM_COMPILER_BASE/bin/cpp
f77=$SYSTEM_COMPILER_BASE/bin/gfortran
f90=$SYSTEM_COMPILER_BASE/bin/gfortran
LINKER=$SYSTEM_CXX
# don't use a version of ar that was
# compiled without --enable-plugins (und{erd}ocumented feature)
# -- see https://gcc.gnu.org/ml/gcc-help/2012-03/msg00100.html
ARCHIVER=ar

# If using gnu compiler toolchain that isn't in standard location it is super
# important to put it's libraries on the LD_LIBARY_PATH
#export LD_LIBRARY_PATH="$SYSTEM_COMPILER_BASE/lib64:$LD_LIBRARY_PATH"

#echo $LD_LIBRARY_PATH
#continue_check()

# if mpi is built by the script these get reset later.
MPI_C_COMPILER="mpicc" 
MPI_CXX_COMPILER="mpiCC"
# mpif90 and mpif77 are depricated in favor of just mpifort.
MPI_F90_COMPILER="mpifort"
MPI_F77_COMPILER="mpifort"
MPI_RUNNER="mpirun"
FORTRAN_LIBS="-lgfortran"


# Tells the TPL builder where MPI is installed.
# Usually /usr/lib/openmpi/
#export MPI_BASE_DIR=""

# Tells the TPL builder where the MPI compilers are.
# Usually /usr/bin
#export MPI_BIN_DIR=""


# Tells the TPL builder which MPI you are using.
# Possible options are: "intel" "open". "open" is recommended
#export MPI_NAME="open"

# Tells the TPL builder which math libraries to use.
# Possible options are: "intel" "netlib blas" "atlas". "netlib blas" is recommended
#export MATH_LIBRARIES="netlib blas"

# Tells the TPL builder where the Atlas/Intel math libraries are if those are used
#export MATH_PATH=""

#Set defaults to what we know works if values are unset.

if [ "$MPI_NAME" == "intel" ] && [ ! "$MATH_LIBRARIES" == "intel" ]; then
    echo "Unable to build Goma deps because non-intel math libraries are used with intel MPI"
    echo "Please switch to openMPI or intel MKLs"
    exit 24
fi

# All user interaction functions like continue_check are here
source "${0%/*}/user-interaction.sh"

function usage() {
    echo "Usage: build-goma-dependencies [options] [library install location]"
    echo "       Options:"
    echo "               -j N   N : Number of make jobs to run"
    echo "               "
    echo "               -h         Display this message"
    echo "               "
    echo "               -m         Disable configure menu and use default settings"
    echo "               "
    echo "       Environment Variables:"
    echo "               CC_NAME=$CC_NAME"
    echo "               MPI_BASE_DIR=$MPI_BASE_DIR"
    echo "               MPI_NAME=$MPI_NAME"
    echo "               MATH_LIBRARIES=$MATH_LIBRARIES"
    echo "               MATH_PATH=$MATH_PATH          (This should be blank unless Atlas math libraries are used)"
    export PRINTED_USAGE="true"
}

MAKE_JOBS=1

while getopts "hmj:" opt; do
    case "${opt}" in
        h)
            usage
            exit 1
            ;;
        j)
            MAKE_JOBS=$OPTARG
            ;;
        m)
            PRINT_MENU="false"
            ;;
        *)
            usage
            exit 1
            ;;
    esac
done
shift $((OPTIND-1))
if [ -z "$1" ]
then
    echo "ERROR: Missing library install location"
    echo
    usage
    exit 1
fi
export GOMA_LIB=`readlink --canonicalize $1`

if ulimit -S -n 2048 &> /dev/null; then
    echo
else
    echo "Compiling Trilinos may require opening up to 2048 files at once."
    echo "But your security settings limit this script to:"
    echo "`ulimit -n` files"
    echo "Please contact your system administrator about fixing this."
    echo "Or continue at your own risk."
    continue_check
fi


if [ "$PRINT_MENU" == "false" ]; then
    if which mpiicc &> /dev/null; then
        MPI_NAME="intel"
        MPI_BASE_DIR="inteldir"
    else
        MPI_NAME="open"
        if [ -d "/usr/lib64/openmpi" ] ; then
            MPI_BASE_DIR="/usr/lib64/openmpi"
        elif [ -d "/usr/lib/openmpi" ] ; then
            MPI_BASE_DIR="/usr/lib/openmpi"
        else
            MPI_BASE_DIR="BUILD"
        fi
    fi
    if [ -z "${CC_NAME}" ]; then
        CC_NAME="gnu"
    fi
    if [ -z "${MATH_LIBRARIES}" ]; then
        MATH_LIBRARIES="netlib blas"
    fi
    USED_MAKE_JOBS="$MAKE_JOBS"
    confirm Goma dependency script
fi

if [ -z "${CC_NAME}" ]; then
    ASKED_USER_ANYTHING="true"
    compiler_test "Goma TPLs"
fi

if [ -z "${MPI_BASE_DIR}" ]; then
    ASKED_USER_ANYTHING="true"
    mpi_test "Goma TPLs"
fi
if [ -z "${MATH_LIBRARIES}" ]; then
    ASKED_USER_ANYTHING="true"
    mkl_test "Goma TPLs"
fi


if [ "$CC_NAME" == "intel" ] && [ "$MPI_NAME" == "open" ]; then
    MPI_BASE_DIR="BUILD"
fi

export MPI_RUNNER="mpirun"

mkdir -p $1

# check that we were able to make that directory
if [ ! -d $1 ]
    then
    echo "ERROR: could not create directory $1"
    exit 1
fi
cd $1
OWNER=$USER

# This is to ensure that goma builds without interruption with the easy build script.
if [ "$ASKED_USER_ANYTHING" == "true" ]; then
    USED_MAKE_JOBS="$MAKE_JOBS"
    confirm Goma dependency script
fi


ARCHIVE_NAMES=("arpack96.tar.gz" \
"patch.tar.gz" \
"hdf5-1.8.20.tar.gz" \
"netcdf-4.4.1.1.tar.gz" \
"parmetis-4.0.3.tar.gz" \
"sparse.tar.gz" \
"superlu_dist-5.1.3.tar.gz" \
"y12m.tar.gz" \
"Trilinos-trilinos-release-12-12-1.tar.gz" \
"MUMPS_5.1.1.tar.gz" \
"SuiteSparse-4.5.5.tar.gz" \
"matio-1.5.10.tar.gz")

#y12m archive is skipped because it stores the number of downloads in the header
#meaning each y12m tar has a unique MD5SUM.
ARCHIVE_MD5SUMS=("fffaa970198b285676f4156cebc8626e" \
"14830d758f195f272b8594a493501fa2" \
"7f2d3fd67106968eb45d133f5a22150f" \
"503a2d6b6035d116ed53b1d80c811bda" \
"f69c479586bf6bb7aff6a9bc0c739628" \
"1566d914d1035ac17b73fe9bc0eed02a" \
"fdee368cba0e95cb0143b6d47915e7a1" \
"SKIP" \
"ecd4606fa332212433c98bf950a69cc7" \
"f15c6b5dd8c71b1241004cd19818259d" \
"0a5b38af0016f009409a9606d2f1b555" \
"d3b6e9d24a04c56036ef57e8010c80f1")

ARCHIVE_URLS=("http://www.caam.rice.edu/software/ARPACK/SRC/arpack96.tar.gz" \
"http://www.caam.rice.edu/software/ARPACK/SRC/patch.tar.gz" \
"https://support.hdfgroup.org/ftp/HDF5/current18/src/hdf5-1.8.20.tar.gz" \
"ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.4.1.1.tar.gz" \
"http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz" \
"http://downloads.sourceforge.net/project/sparse/sparse/sparse1.4b/sparse1.4b.tar.gz" \
"http://codeload.github.com/xiaoyeli/superlu_dist/tar.gz/v5.1.3" \
"http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz\\&filename=y12m%2Fy12m.f" \
"https://github.com/trilinos/Trilinos/archive/trilinos-release-12-12-1.tar.gz" \
"http://mumps.enseeiht.fr/MUMPS_5.1.1.tar.gz" \
"http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-4.5.5.tar.gz" \
"http://superb-dca2.dl.sourceforge.net/project/matio/matio/1.5.10/matio-1.5.10.tar.gz")

# You can't call the ARPACK patch ARPACK or it will think it is already extracted
# When in reality it isn't
ARCHIVE_DIR_NAMES=("ARPACK" \
"FAKE_DIR_FOR_ARPACK_PATCH" \
"hdf5-1.8.20" \
"netcdf-4.4.1.1" \
"parmetis-4.0.3" \
"sparse" \
"superlu_dist-5.1.3" \
"y12m" \
"Trilinos-trilinos-release-12-12-1" \
"MUMPS_5.1.1" \
"SuiteSparse" \
"matio-1.5.10")

ARCHIVE_HOMEPAGES=("http://www.caam.rice.edu/software/ARPACK/" \
"https://www.hdfgroup.org/" \
"https://www.unidata.ucar.edu/software/netcdf/" \
"http://glaros.dtc.umn.edu/gkhome/metis/metis/overview" \
"http://sparse.sourceforge.net/" \
"https://github.com/xiaoyeli/superlu_dist" \
"http://doi.org/10.1007/3-540-10874-2" \
"https://trilinos.org/" \
"http://mumps.enseeiht.fr/" \
"http://faculty.cse.tamu.edu/davis/suitesparse.html" \
"https://sourceforge.net/projects/matio/")

ARCHIVE_REAL_NAMES=("ARPACK96" \
"HDF5" \
"NetCDF" \
"Metis" \
"Sparse" \
"SuperLU_DIST" \
"y12m" \
"Trilinos" \
"MUMPS" \
"SuiteSparse" \
"MATIO")

if [[ "$MPI_BASE_DIR" == "BUILD" ]]; then
    ARCHIVE_NAMES+=("openmpi-2.1.1.tar.gz")
    ARCHIVE_MD5SUMS+=("bee36445ad1ef452a378f446f054804d")
    ARCHIVE_URLS+=("https://www.open-mpi.org/software/ompi/v2.1/downloads/openmpi-2.1.1.tar.gz")
    ARCHIVE_DIR_NAMES+=("openmpi_2.1.1")
    ARCHIVE_HOMEPAGES+=("https://www.open-mpi.org/")
    ARCHIVE_REAL_NAMES+=("OpenMPI")
#    If the user says to use the intel library but sets the base dir to BUILD then there's a conflict here
    MPI_NAME="open"
    MPI_BASE_DIR="$GOMA_LIB/openmpi-2.1.1"
    MPI_IS_BUILT_FROM_SCRATCH="true"
    echo
    echo "OpenMPI 2.1.1 will be built."
    echo "You will have to add ${MPI_BASE_DIR}/lib to your LD_LIBRARY_PATH"
    echo "And ${MPI_BASE_DIR}/bin to your PATH environment variables"
    echo "Before running Goma"
    echo
    continue_check
fi

if command -v cmake; then
    cmake_vers=$(cmake --version |grep "version" | awk '{print $NF}')
    if [[ "$cmake_vers" = $(echo -e "$cmake_vers\n2.8.12.2\n" | sort -V |tail -n1) ]]; then
	build_cmake="false"
    else
	build_cmake="true"
    fi
else
    build_cmake="true"
fi

if [ "$build_cmake" == "false" ] ; then
    echo "Native cmake found newer than 2.8.12.2, skipping download"
else
    ARCHIVE_NAMES+=("cmake-3.9.5.tar.gz")
    ARCHIVE_MD5SUMS+=("0922130d0e0c142a88e58c6e4fef4d7d")
    ARCHIVE_URLS+=("https://cmake.org/files/v3.9/cmake-3.9.5.tar.gz")
    ARCHIVE_DIR_NAMES+=("cmake-3.9.5")
    ARCHIVE_HOMEPAGES+=("https://cmake.org/")
    ARCHIVE_REAL_NAMES+=("CMake")
    echo "Cmake not found, will build."
    echo
fi



function setMPIvars() {
    if [[ "$MPI_NAME" == "intel" ]]; then
        export MPI_BASE_DIR="${I_MPI_ROOT}/intel64"
        export PATH=${MPI_BASE_DIR}/bin:${PATH}
        export LD_LIBRARY_PATH=${MPI_BASE_DIR}/lib:${LD_LIBRARY_PATH}
        export MPI_FORTRAN_LIB="-lmpifort -lmpicxx"
        export I_MPI_SHM_LMT=shm
        # Variable used by BLACS to determine MPI used
        export MKL_BLACS_MPI="INTELMPI"
    elif [[ "$MPI_NAME" == "open" ]]; then
        export MPI_FORTRAN_LIB="-lmpi_mpifh"
        export PATH=${MPI_BASE_DIR}/bin:${PATH}
	echo "in setMPIvars() LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
        export LD_LIBRARY_PATH=${MPI_BASE_DIR}/lib:${LD_LIBRARY_PATH}
	echo "in setMPIvars() LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
        unset MKL_BLACS_MPI
    else
        echo "Unknown MPI: $MPI_NAME. You can continue, but what will happen?"
        MPI_NAME="open"
        continue_check
    fi
}

function setCompilerVars() {
    if [[ "$CC_NAME" == "intel" ]]; then
        #Special flag only needed by intel compiler because reasons? No really why does GCC not need -fopenmp?
        export COMPILER_FLAG_MPI="-qopenmp"
        #TODO: Test if this fails with openmpi (probably)
        if [ "$MPI_NAME" == "open" ]; then
            export MPI_C_COMPILER="mpicc"
            export MPI_CXX_COMPILER="mpic++"
            export MPI_F90_COMPILER="mpifort"
            export MPI_F77_COMPILER="mpifort"
        else
            export MPI_C_COMPILER="mpiicc"
            export MPI_CXX_COMPILER="mpiicpc"
            export MPI_F90_COMPILER="mpiifort"
            export MPI_F77_COMPILER="mpiifort"
        fi
        export SYSTEM_CC="icc"
        SYSTEM_CXX="icpc"
        export SYSTEM_FC="ifort"
        SYSTEM_CPP="mpiicc -E"
        export LINKER="xild"
        export ARCHIVER="xiar"
        export f77="ifort"
        export f90="ifort"
        EXTRA_CXX_FLAGS=""
        export FORTRAN_LIBS="-L$INTEL_PARALLEL_STUDIO_ROOT/lib/intel64 -lifcore"
    elif [[ "$CC_NAME" == "gnu" ]]; then
        export MPI_C_COMPILER="mpicc"
        export MPI_CXX_COMPILER="mpiCC"
        # mpif90 and mpif77 are depricated in favor of just mpifort.
        if [ "$MPI_NAME" == "intel" ]; then
            # Intel calls their mpi fortran compiler mpifc (for some reason
            export MPI_F90_COMPILER="mpifc"
            export MPI_F77_COMPILER="mpifc"
        else
            export MPI_F90_COMPILER="mpifort"
            export MPI_F77_COMPILER="mpifort"
        fi
        export MPI_RUNNER="mpirun"
        export FORTRAN_LIBS="-lgfortran"
        export SYSTEM_CC=gcc
        SYSTEM_CXX=g++
        export SYSTEM_FC=gfortran
        SYSTEM_CPP=cpp
        export f77=gfortran
        export f90=gfortran
        export LINKER="ld"
        export ARCHIVER="ar"
    elif [[ "$CC_NAME" == "user" ]]; then
	echo "MPI_C_COMPILER=$MPI_C_COMPILER"
        echo "MPI_CXX_COMPILER=$MPI_CXX_COMPILER"
        echo "MPI_F90_COMPILER=$MPI_F90_COMPILER"
        echo "MPI_F77_COMPILER=$MPI_F77_COMPILER"
        echo "MPI_RUNNER=$MPI_RUNNER"
        echo "FORTRAN_LIBS=$FORTRAN_LIBS"
        echo "SYSTEM_CC=$SYSTEM_CC"
        echo "SYSTEM_CXX=$SYSTEM_CXX"
        echo "SYSTEM_FC=$SYSTEM_FC"
        echo "SYSTEM_CPP=$SYSTEM_CPP"
        echo "f77=$f77"
        echo "f90=$f90"
        echo "LINKER=$LINKER"
        echo "ARCHIVER=$ARCHIVER"
    else
        echo "Unsupported compiler: $CC_NAME. You can continue but the build probably will fail."
        CC_NAME="gnu"
        continue_check

    fi
}

function setMathVars {
    if [[ "$MATH_LIBRARIES" == "intel" ]]; then
        USING_MKLS="ON"
        MKL_LIBRARY_NAME="mkl_sequential;mkl_core"
        MKL_LIBRARY_DIR="${MKLROOT}/lib/intel64"
        MKL_INCLUDE_DIR="${MKLROOT}/include"
        # Libraries included as part of Intel Math Kernel Library
        export BLAS_LIBRARY_DIR="${MKLROOT}/lib/intel64"
        export BLAS_LIBRARY_NAME_ARG="-lmkl_blas95_ilp64 -lmkl_sequential -lmkl_core"
        BLAS_LIBRARY_NAME="mkl_blas95_ilp64;mkl_sequential;mkl_core" #In a format Trilinos can read
        export SCALAPACK_LIBRARY_DIR="${BLAS_LIBRARY_DIR}"
        export SCALAPACK_LIBRARY_NAME="mkl_scalapack_ilp64"
        export SCALAPACK_LIBRARY_NAME_ARG="-L${BLAS_LIBRARY_DIR} -lmkl_scalapack_ilp64"
        SCALAPACK_INCLUDE_DIR="${MKLROOT}/include"
        export BLACS_LIBRARY_NAME="mkl_blacs_intelmpi_ilp64"
        export BLACS_LIBRARY_PATH="${BLAS_LIBRARY_DIR}"
        export BLACS_LIBRARY_NAME_ARG="-L${SCALAPACK_LIBRARY_DIR} -lmkl_blacs_intelmpi_ilp64"
        if [ "$CC_NAME" == "intel" ]; then
            BLAS_LIBRARY_NAME="$BLAS_LIBRARY_NAME;mkl_intel_ilp64"
            export BLAS_LIBRARY_NAME_ARG="$BLAS_LIBRARY_NAME_ARG -lmkl_intel_ilp64"
            export BLAS_FLAGS="-mkl"
        else
            echo
            echo "It is highly recommended you use the intel compiler with intel MKLs."
            echo
            echo "Proceed at your own risk."
            continue_check
            BLAS_LIBRARY_NAME="$BLAS_LIBRARY_NAME;mkl_gf_ilp64"
            export BLAS_LIBRARY_NAME_ARG="$BLAS_LIBRARY_NAME_ARG -lmkl_gf_ilp64"
            # Link for non-intel compilers with intel MKL
            export BLAS_FLAGS="-m64 -DMKL_ILP64 -I${MKLROOT}/include"
            if [ "$MPI_NAME" == "open" ]; then
                export NON_INTEL_BLAS_LINK="${MKLROOT}/lib/intel64/libmkl_blas95_ilp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_openmpi_ilp64.a -Wl,--end-group -ldl -lz"
                export NON_INTEL_BLAS_LIBRARY="$NON_INTEL_BLAS_LINK"
                # Be sure to BUILD ScaLAPACK if there's no intelMPI or EVERYTHING BREAKS?! Probably, but not 100% sure yet
                # Also built ScaLAPACK includes blacs so unset that too
                export SCALAPACK_LIBRARY_DIR="${GOMA_LIB}/scalapack-2.0.2"
                export SCALAPACK_LIBRARY_NAME="scalapack"
                export SCALAPACK_LIBRARY_NAME_ARG="-L${SCALAPACK_LIBRARY_DIR} -lscalapack"
                SCALAPACK_INCLUDE_DIR="${GOMA_LIB}/scalapack-2.0.2/SRC"
                ARCHIVE_NAMES+=("scalapack-2.0.2.tgz")
                ARCHIVE_MD5SUMS+=("2f75e600a2ba155ed9ce974a1c4b536f")
                ARCHIVE_URLS+=("http://www.netlib.org/scalapack/scalapack-2.0.2.tgz")
                ARCHIVE_DIR_NAMES+=("scalapack-2.0.2")
                ARCHIVE_HOMEPAGES+=("http://www.netlib.org/scalapack/")
                ARCHIVE_REAL_NAMES+=("ScaLAPACK")
                export BLACS_LIBRARY_NAME="scalapack"
                unset BLACS_LIBRARY_PATH
                unset BLACS_LIBRARY_NAME_ARG
            else
                export NON_INTEL_BLAS_LINK="${MKLROOT}/lib/intel64/libmkl_blas95_ilp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_ilp64.a ${MKLROOT}/lib/intel64/libmkl_scalapack_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a -Wl,--end-group -ldl -lz"
                export NON_INTEL_BLAS_LIBRARY="$NON_INTEL_BLAS_LINK"
            fi
        fi
        export LAPACK_LIBRARY_DIR="${BLAS_LIBRARY_DIR}"
        export LAPACK_LIBRARY_NAME_ARG="-lmkl_lapack95_ilp64"
        LAPACK_LIBRARY_NAME="mkl_lapack95_ilp64"
    elif [[ "$MATH_LIBRARIES" == "netlib blas" ]]; then
        USING_MKLS="OFF"
        #Add packages that otherwise come preinstalled in intel compiler.
        ARCHIVE_NAMES+=("blas-3.7.1.tgz" \
        "lapack-3.8.0.tar.gz" \
        "scalapack-2.0.2.tgz")
        ARCHIVE_MD5SUMS+=("cd132aea6f7055a49aa48ca0a61e7cd5" \
        "96591affdbf58c450d45c1daa540dbd2" \
        "2f75e600a2ba155ed9ce974a1c4b536f" )
        ARCHIVE_URLS+=("http://www.netlib.org/blas/blas-3.7.1.tgz" \
        "http://www.netlib.org/lapack/lapack-3.8.0.tar.gz" \
        "http://www.netlib.org/scalapack/scalapack-2.0.2.tgz" )
        ARCHIVE_DIR_NAMES+=("BLAS-3.7.1" \
        "lapack-3.8.0" \
        "scalapack-2.0.2" )
        ARCHIVE_HOMEPAGES+=("http://www.netlib.org/blas/" \
        "http://www.netlib.org/lapack/" \
        "http://www.netlib.org/scalapack/")
        ARCHIVE_REAL_NAMES+=("BLAS" \
        "LAPACK" \
        "ScaLAPACK")

        BLAS_LIBRARY_DIR="${GOMA_LIB}/BLAS-3.7.1"
        BLAS_LIBRARY_NAME_ARG="-L${BLAS_LIBRARY_DIR} -lblas"
        BLAS_LIBRARY_NAME="blas"
        BLACS_LIBRARY_NAME="scalapack"
        LAPACK_LIBRARY_DIR="${GOMA_LIB}/lapack-3.8.0"
        LAPACK_LIBRARY_NAME_ARG="-L${LAPACK_LIBRARY_DIR} -llapack"
        LAPACK_LIBRARY_NAME="liblapack.a"
        SCALAPACK_LIBRARY_DIR="${GOMA_LIB}/scalapack-2.0.2"
        SCALAPACK_LIBRARY_NAME="scalapack"
        SCALAPACK_LIBRARY_NAME_ARG="-L${SCALAPACK_LIBRARY_DIR} -lscalapack"
        SCALAPACK_INCLUDE_DIR="${GOMA_LIB}/scalapack-2.0.2/SRC"
        # Alternative to the intel "blas_flags" variable
        NON_INTEL_BLAS_LIBRARY="${BLAS_LIBRARY_DIR}/libblas.a"
        NON_INTEL_BLAS_LINK="-L${BLAS_LIBRARY_DIR} -l${BLAS_LIBRARY_NAME} ${FORTRAN_LIBS}"
        SUITESPARSE_NON_INTEL_LAPACK_LINK="${LAPACK_LIBRARY_NAME_ARG} ${FORTRAN_LIBS}"


    elif [[ "$MATH_LIBRARIES" == "atlas" ]]; then
        USING_MKLS="OFF"
        if [ -n "$MATH_PATH" ]; then
            ATLAS_PATH=${MATH_PATH}
        else
            if [ -e /usr/lib64/libblas.a ]; then
                ATLAS_PATH=/usr/lib64
            elif [ -e /usr/lib/libblas.a ]; then
                ATLAS_PATH=/usr/lib
            else
                echo "Unable to find atlas libraries. Please enter their path"
                read ATLAS_PATH
            fi
        fi
        # Please replace with your blas and lapack library names
        BLAS_LIBRARY_NAME="blas"
        LAPACK_LIBRARY_NAME="lapack"

        export BLAS_LIBRARY_DIR=${ATLAS_PATH}
        export BLAS_LIBRARY_NAME_ARG="-lblas"
        export LAPACK_LIBRARY_DIR=${ATLAS_PATH}
        export LAPACK_LIBRARY_NAME_ARG="-L${ATLAS_PATH} -llapack"
        export BLACS_LIBRARY_NAME="scalapack"
        export SCALAPACK_LIBRARY_NAME="scalapack"
        export SCALAPACK_LIBRARY_DIR="${GOMA_LIB}/scalapack-2.0.2"
        export SCALAPACK_LIBRARY_NAME_ARG="-L{SCALAPACK_LIBRARY_DIR} -lscalapack"
        SCALAPACK_INCLUDE_DIR="${GOMA_LIB}/scalapack-2.0.2/SRC"
        export NON_INTEL_BLAS_LIBRARY="${BLAS_LIBRARY_DIR}/lib${BLAS_LIBRARY_NAME}.a"
        export NON_INTEL_BLAS_LINK="-L${BLAS_LIBRARY_DIR} -l${BLAS_LIBRARY_NAME} ${FORTRAN_LIBS}"
        ARCHIVE_NAMES+=("scalapack-2.0.2.tgz")
        ARCHIVE_MD5SUMS+=("2f75e600a2ba155ed9ce974a1c4b536f" )
        ARCHIVE_URLS+=("http://www.netlib.org/scalapack/scalapack-2.0.2.tgz" )
        ARCHIVE_DIR_NAMES+=("scalapack-2.0.2" )
        ARCHIVE_HOMEPAGES+=("http://www.netlib.org/scalapack/")
        ARCHIVE_REAL_NAMES+=("ScaLAPACK")

        export SUITESPARSE_NON_INTEL_LAPACK_LINK="-L${LAPACK_LIBRARY_NAME_ARG} ${FORTRAN_LIBS}"
        echo
        echo "Atlas libraries are named differently across every distribution."
        echo "Please edit the build-goma-dep-trilinos-12.sh file and provide"
        echo "BLAS_LIBRARY_NAME and LAPACK_LIBRARY_NAME the appropriate values"
        echo
        echo "If you have already done this then please continue"
        continue_check
    else
        echo "Unsupported math libraries: $MATH_LIBRARIES. You can continue"
        echo "But netlib blas libraries will be used instead"
        MATH_LIBRARIES="netlib blas"
        continue_check
    fi
}
echo "You are going to build with CC_NAME=$CC_NAME"
setCompilerVars
continue_check
setMPIvars
setMathVars

#export MPI_BASE_DIR=$MPI_BASE_DIR
export GOMA_LIB
export OWNER
export MPI_BASE_DIR
export MAKE_JOBS



read -d '' MUMPS_PATCH << "EOF"
5a6,8
> GOMA_LIB=__GOMA_LIB__
> 
> 
37,38c40,41
< #LMETISDIR = /local/metis/
< #IMETIS    = # Should be provided if you use parmetis
---
> LMETISDIR = $(GOMA_LIB)/ParMetis-3.1.1
> IMETIS    = -I$(GOMA_LIB)/ParMetis-3.1.1 
47c50
< #LMETIS    = -L$(LMETISDIR) -lparmetis -lmetis
---
> LMETIS    = -L$(LMETISDIR) -lparmetis -lmetis
52c55
< ORDERINGSF  = -Dpord
---
> ORDERINGSF  = -Dpord -Dparmetis -Dmetis
68,70c71,73
< CC = gcc
< FC = gfortran
< FL = gfortran
---
> CC = $(MPI_C_COMPILER)
> FC = $(MPI_F77_COMPILER)
> FL = $(MPI_CXX_COMPILER)
72,74c75,76
< #RANLIB = ranlib
< RANLIB  = echo
< SCALAP  = /local/SCALAPACK/libscalapack.a /local/BLACS/LIB/blacs_MPI-LINUX-0.a   /local/BLACS/LIB/blacsF77init_MPI-LINUX-0.a  /local/BLACS/LIB/blacs_MPI-LINUX-0.a
---
> RANLIB = ranlib
> SCALAP  = $(SCALAPACK_LIBRARY_DIR)/lib$(SCALAPACK_LIBRARY_NAME).a 
76,78c78,80
< INCPAR = -I/usr/local/mpich/include
< # LIBPAR = $(SCALAP)  -L/usr/local/lib/ -llammpio -llamf77mpi -lmpi -llam -lutil -ldl -lpthread
< LIBPAR = $(SCALAP)  -L/usr/local/mpich/lib/ -lmpich 
---
> MPI_HOME=$(MPI_BASE_DIR)
> INCPAR = -I$(MPI_BASE_DIR)/include
> LIBPAR = $(SCALAP)  -L$(MPI_BASE_DIR) -lmpi -lutil -ldl -lpthread $(MPI_FORTRAN_LIB)
83c85
< LIBBLAS = -L/local/BLAS -lblas
---
> LIBBLAS = -L$(BLAS_LIBRARY_DIR) $(BLAS_LIBRARY_NAME_ARG)
EOF

read -d '' NETCDF_PATCH << "EOF"
229c229
< #define NC_MAX_DIMS	1024
---
> #define NC_MAX_DIMS	65536
231c231
< #define NC_MAX_VARS	8192
---
> #define NC_MAX_VARS	524288
233c233
< #define NC_MAX_VAR_DIMS	1024 /**< max per variable dimensions */
---
> #define NC_MAX_VAR_DIMS	8 /**< max per variable dimensions */
EOF


read -d '' BLAS_PATCH << "EOF"
18c18
< FORTRAN  = gfortran
---
> FORTRAN  = $(SYSTEM_FC)
22c22
< LOADER   = gfortran
---
> LOADER   = $(SYSTEM_FC)
EOF


#Etime isn't provided anymore so comment it out
read -d '' ARPACK_SECOND_PATCH << "EOF"
24c24
<       EXTERNAL           ETIME
---
> *     EXTERNAL           ETIME
EOF


read -d '' LAPACK_PATCH << "EOF"
11c11
< CC     = gcc
---
> CC     = $(SYSTEM_CC)
22c22
< FORTRAN = gfortran
---
> FORTRAN = $(SYSTEM_FC)
30c30
< LOADER   = gfortran
---
> LOADER   = $(SYSTEM_FC)
50c50
< TIMER = INT_ETIME
---
> #TIMER = INT_ETIME
54c54
< #TIMER = INT_CPU_TIME
---
> TIMER = INT_CPU_TIME
81c81
< BLASLIB      = ../../librefblas.a
---
> BLASLIB      = $(GOMA_LIB)/BLAS-3.7.1/blas_LINUX.a
EOF


read -d '' SUPERLU_PATCH << "EOF"
19c19
< PLAT      = _power5
---
> PLAT      =
24c24
< DSuperLUroot  = $(HOME)/Release_Codes/superlu_dist-2.3
---
> DSuperLUroot  = $(GOMA_LIB)/superlu_dist-2.3
28c28
< BLASLIB       = -lessl
---
> BLASLIB       =
31,32c31,32
< METISLIB        = -L/project/projectdirs/sparse/xiaoye/parmetis-3.1/64 -lmetis
< PARMETISLIB   = -L/project/projectdirs/sparse/xiaoye/parmetis-3.1/64 -lparmetis
---
> METISLIB        = $(GOMA_LIB)/ParMetis-3.1.1 -lmetis
> PARMETISLIB   = $(GOMA_LIB)/ParMetis-3.1.1 -lparmetis
42,43c42,43
< ARCHFLAGS     = -X64 cr
< RANLIB        = ranlib
---
> ARCHFLAGS     = -rc
> RANLIB        = echo
48c48
< CC            = mpcc_r
---
> CC            = $(MPI_C_COMPILER)
50c50
< CFLAGS          = -D_SP -qarch=pwr5 -qalias=allptrs -q64
---
> CFLAGS          = -D_SP 
55c55
< NOOPTS        = -q64
---
> NOOPTS        = 
60,62c60,61
< FORTRAN         = mpxlf90_r
< FFLAGS          = -WF,-Dsp -O3 -Q -qstrict -qfixed -qinit=f90ptr -qarch=pwr5\
<                   -q64 -qintsize=8
---
> FORTRAN         = $(MPI_F90_COMPILER)
> FFLAGS        = -O3 -Q 
65c64
< LOADER    = mpxlf90_r
---
> LOADER    = $(MPI_C_COMPILER)
68c67
< LOADOPTS  = -q64
---
> LOADOPTS  = 
74c73
< CDEFS        = -DNoChange
---
> CDEFS        = -DAdd_

EOF

read -r -d '' MUMPS_MAKE_PATCH<<"EOF"
51c51
< \t$(FL) -o $@ $(OPTL) $@.o $(LIBDMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)
---
> \t$(CC) -o $@ $(OPTL) $@.o $(LIBDMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)

EOF


#https://www.dealii.org/8.5.0/external-libs/trilinos.html See this for details on this patch.
#Allows building Trilinos with modern versions of Lapack.
read -r -d '' TRILINOS_EPETRA_LAPACK_PATCH << "EOF"
175c175
< #define DGGSVD_F77  F77_BLAS_MANGLE(dggsvd,DGGSVD)
---
> #define DGGSVD_F77  F77_BLAS_MANGLE(dggsvd3,DGGSVD3)
230c230
< #define SGGSVD_F77  F77_BLAS_MANGLE(sggsvd,SGGSVD)
---
> #define SGGSVD_F77  F77_BLAS_MANGLE(sggsvd3,SGGSVD3)
EOF


cd $GOMA_LIB

function mychecksum {
    local count=$2
    local archive=$1
    local MD5SAVED=${ARCHIVE_MD5SUMS[count]}
    if [ ! "$MD5SAVED" == "SKIP" ]; then
        local MD5ARCHIVE=($(md5sum $archive))
        if [ $MD5SAVED != $MD5ARCHIVE ]; then
            echo "Issue checksum with archive:"
            echo "$MD5SAVED"
            echo "$MD5ARCHIVE"
            echo $archive
            continue_check
        fi
    fi
}

echo
echo -e "\e[31m \e[1m"
echo
echo "The following libraries will be built for Goma"
count=0
for i in ${ARCHIVE_REAL_NAMES[@]}; do
    echo "$i"
    echo "${ARCHIVE_HOMEPAGES[count]}"
    echo
    count=$(( $count + 1 ))
done
echo
echo "You are encouraged to investigate the licenses of these software packages before continuing."
echo
echo -e "\e[0m"
continue_check

mkdir -p tars
cd tars
#downloads
count=0
for i in ${ARCHIVE_NAMES[@]}; do
    echo "Check for $i at ${ARCHIVE_URLS[count]}"
    if ! [ -f $i ]
    then
        wget ${ARCHIVE_URLS[count]} -O $i
        mychecksum $i $count
    else
        mychecksum $i $count
    fi
    cd ..

    if [ -d ${ARCHIVE_DIR_NAMES[count]} ]
    then
	echo "already extracted ${i}"
    else
	if ! tar tf $i &> /dev/null; then
	    tar -xf tars/$i
	fi
    fi
    count=$(( $count + 1 ))
    cd tars
done

#continue_check
export CXX=${SYSTEM_CXX}

if [ "$build_cmake" == "false" ] ; then
    echo "Native cmake found of sufficient version, skipping build"
else
    cd $GOMA_LIB/cmake-3.9.5
    if [ -f bin/cmake ]
    then
        echo "cmake is already built"
    else
        CC=$SYSTEM_CC CXX=$SYSTEM_CXX FC=$SYSTEM_FC ./bootstrap --prefix=$GOMA_LIB/cmake-3.9.5
        make -j$MAKE_JOBS
        make install
    fi
    export PATH=$GOMA_LIB/cmake-3.9.5/bin:$PATH
fi


if [[ "$MPI_IS_BUILT_FROM_SCRATCH" == "true" ]]; then
    #continue_check
    #make openmpi
    cd $GOMA_LIB/openmpi-2.1.1
    if [ -f bin/ompi_info ]
    then
        echo "openmpi is already built"
    else
        CC=${SYSTEM_CC} CXX=${SYSTEM_CXX} FC=${SYSTEM_FC} CPP=${SYSTEM_CPP} ./configure --prefix=$GOMA_LIB/openmpi-2.1.1
        make -j$MAKE_JOBS CC=${SYSTEM_CC} CXX=${SYSTEM_CXX} FC=${SYSTEM_FC} CPP=${SYSTEM_CPP}
        make install
    fi
    cd ..
else
    echo "Using custom mpi: $mpi"
fi

#continue_check

cd $GOMA_LIB

#hdf5
if [ -e hdf5-1.8.20/lib/libhdf5.a ]
then
    echo "hdf5 already built"
else
    if ! [ -e hdf5-1.8.20/.goma-extracted ]
    then
	mv hdf5-1.8.20 tmpdir
	mkdir hdf5-1.8.20
	mv tmpdir hdf5-1.8.20/src
	touch hdf5-1.8.20/.goma-extracted
    fi

    cd hdf5-1.8.20/src
    CC="$MPI_C_COMPILER" CPP="$MPI_C_COMPILER -E" AR=${ARCHIVER} ./configure --enable-shared=off --prefix=$GOMA_LIB/hdf5-1.8.20 --enable-parallel
    make -j$MAKE_JOBS
    make install
    cd ../..
fi
#continue_check
cd $GOMA_LIB

#matio
if [ -e matio-1.5.10/lib/libmatio.a ]
then
    echo "matio already built"
else
    if ! [ -e matio-1.5.10/.goma-extracted]
    then
	mv matio-1.5.10 tmpdir
	mkdir matio-1.5.10
	mv tmpdir matio-1.5.10/src
	touch matio-1.5.10/.goma-extracted
    fi
    cd matio-1.5.10/src
    CC=${MPI_C_COMPILER} LD=${MPI_CXX_COMPILER} AR=${ARCHIVER} LIBS="-ldl" ./configure --with-hdf5=${GOMA_LIB}/hdf5-1.8.20 --prefix=${GOMA_LIB}/matio-1.5.10 --enable-shared=off
    make -j$MAKE_JOBS
    make install
fi
cd $GOMA_LIB

#continue_check

#netcdf
if [ -e netcdf-4.4.1.1/lib/libnetcdf.a ]
then
    echo "netcdf already built"
else
    if ! [ -e netcdf-4.4.1.1/.goma-extracted ]
    then
	mv netcdf-4.4.1.1 tmpdir
	mkdir netcdf-4.4.1.1
	mv tmpdir netcdf-4.4.1.1/src
	touch netcdf-4.4.1.1/.goma-extracted
    fi
    cd $GOMA_LIB/netcdf-4.4.1.1/src
    cd include
    echo "$NETCDF_PATCH" > netcdf.patch
    patch -f --ignore-whitespace netcdf.h < netcdf.patch
    cd ..
    export CPPFLAGS=-I$GOMA_LIB/hdf5-1.8.20/include
#    export LDFLAGS=-L$GOMA_LIB/hdf5-1.8.20/lib
    echo $CPPFLAGS
    echo $LDFLAGS

    CC=${MPI_C_COMPILER} CFLAGS="-I${GOMA_LIB}/hdf5-1.8.20/include" \
      CPP="${MPI_C_COMPILER} -E" \
      CPPFLAGS="-I${GOMA_LIB}/hdf5-1.8.20/include" \
      LDFLAGS="-L${GOMA_LIB}/hdf5-1.8.20/lib" \
      ./configure \
      --prefix=$GOMA_LIB/netcdf-4.4.1.1 \
      --enable-shared=off \
      --disable-dap

    make -j$MAKE_JOBS
    make install
    cd ../..
fi

#continue_check

 #make BLAS
if [[ "$MATH_LIBRARIES" == "netlib blas" ]]; then
    cd $GOMA_LIB/BLAS-3.7.1
    if [ -f libblas.a ]
    then
        echo "BLAS already built"
    else
        echo $("$BLAS_PATCH") > makepatch.inc
        patch -f make.inc < makepatch.inc
        make -j$MAKE_JOBS cc=${MPI_C_COMPILER} ccflags=${COMPILER_FLAG_MPI} cxx=${MPI_CXX_COMPILER} cxxflags=${COMPILER_FLAG_MPI}
        cp blas_LINUX.a libblas.a
    fi
    # For some reason this junk apple file is included with BLAS 3.7.1 from netlib.
    rm -f "$GOMA_LIB/._BLAS-3.7.1"
    export LD_LIBRARY_PATH="${GOMA_LIB}/BLAS-3.7.1:$LD_LIBRARY_PATH"
fi
#continue_check

#parmetis patch
read -d '' PARMETIS_PATCH << "EOF"
55c55,58
<     CONFIG_FLAGS += -DCMAKE_C_COMPILER=$(cc)
---
>     CONFIG_FLAGS += -DCMAKE_C_COMPILER="$(cc)"
> endif
> ifneq ($(ccflags), not-set)
>     CONFIG_FLAGS += -DCMAKE_C_FLAGS="$(ccflags)"
58c61,64
<     CONFIG_FLAGS += -DCMAKE_CXX_COMPILER=$(cxx)
---
>     CONFIG_FLAGS += -DCMAKE_CXX_COMPILER="$(cxx)"
> endif
> ifneq ($(cxxflags), not-set)
>     CONFIG_FLAGS += -DCMAKE_CXX_FLAGS="$(cxxflags)"
EOF

#make parMetis
cd $GOMA_LIB
if [ -e parmetis-4.0.3/lib/libparmetis.a ]
then
    echo "ParMetis already Built"
else
    if [ ! -e parmetis-4.0.3/src/Makefile ]; then
        mv parmetis-4.0.3 tmpdir
        mkdir parmetis-4.0.3
        mv tmpdir parmetis-4.0.3/src
    fi
    cd parmetis-4.0.3/src
    echo "$PARMETIS_PATCH" > parmetis_make.patch
    patch -f Makefile < parmetis_make.patch
    echo "make config parmetis?"
    make config cc=${MPI_C_COMPILER} ccflags=${COMPILER_FLAG_MPI} cxx=${MPI_CXX_COMPILER} cxxflags=${COMPILER_FLAG_MPI} prefix=$GOMA_LIB/parmetis-4.0.3
    make
    make install
    cd ..
    if [ -d include ]; then
	cp src/metis/include/metis.h include
    fi
    if [ -d lib ]; then
	cp src/build/Linux-x86_64/libmetis/libmetis.a lib/
    fi
fi


#make ARPACK
cd $GOMA_LIB/ARPACK
if [ -e libarpack_x86_64.a ]
then
    echo "ARPACK already built"
else
    cat > ARmake.patch << EOF
28c28
< home = \$(HOME)/ARPACK
---
> home = $GOMA_LIB/ARPACK
35c35
< PLAT = SUN4
---
> PLAT = x86_64
104,105c104,105
< FC      = f77
< FFLAGS	= -O -cg89
---
> FC      = $MPI_F77_COMPILER
> FFLAGS	= -O $BLAS_FLAGS
115c115
< MAKE    = /bin/make
---
> MAKE    = make
EOF

    patch ARmake.inc < ARmake.patch
    # Documentation says this should always be needed but in reality only intel MKL requires it.
    if [[ "$MATH_LIBRARIES" == "intel" ]]; then
        echo $("$ARPACK_SECOND_PATCH") > UTIL/second.f.patch
        patch UTIL/second.f < UTIL/second.f.patch
    fi
    make all
    mkdir lib
    cp libarpack_x86_64.a lib/libarpack.a
fi

#make SuperLU
cd $GOMA_LIB/superlu_dist-5.1.3
if [ -e lib/libsuperludist.a ]
then
    echo "SuperLU_DIST already built"
else
    cat > make.inc << EOF
SuperLUroot     =  $GOMA_LIB/superlu_dist-5.1.3
DSUPERLULIB     = $GOMA_LIB/superlu_dist-5.1.3/lib/libsuperludist.a

# BLASDEF       = -DUSE_VENDOR_BLAS

LIBS            = $GOMA_LIB/superlu_dist-5.1.3/lib/libsuperludist.a $NON_INTEL_BLAS_LIBRARY $GOMA_LIB/parmetis-4.0.3/lib/libparmetis.a $GOMA_LIB/parmetis-4.0.3/lib/libmetis.a ${FORTRAN_LIBS}

#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         = /usr/bin/ar
ARCHFLAGS    = cr
RANLIB       = /usr/bin/ranlib

CC           = $MPI_C_COMPILER
CFLAGS       = -O3 -DNDEBUG -DUSE_VENDOR_BLAS -DDEBUGlevel=0 -DPRNTlevel=0 -std=c99 -g -I$GOMA_LIB/parmetis-4.0.3/include
CDEFS = -DAdd_
# CFLAGS       += -D
# CFLAGS       +=
NOOPTS       = -O0
FORTRAN      = $MPI_F90_COMPILER

LOADER       = $MPI_CXX_COMPILER $BLAS_FLAGS
LOADOPTS     = -Wl,-rpath,$GOMA_LIB/superlu_dist-5.1.3/lib
EOF
    mkdir -p lib
    make
fi

if [[ "$MATH_LIBRARIES" == "netlib blas" ]]; then
    #make lapack
    cd $LAPACK_LIBRARY_DIR
    if [ -e liblapack.a ]
    then
        echo "LAPACK already built"
    else
        mv make.inc.example make.inc
        echo $("$LAPACK_PATCH") > make.patch
        patch make.inc < make.patch
        make lapacklib -j$MAKE_JOBS cc=${MPI_C_COMPILER} ccflags=${COMPILER_FLAG_MPI} cxx=${MPI_CXX_COMPILER} cxxflags=${COMPILER_FLAG_MPI}
#        cp lapack_LINUX.a liblapack.a
    fi
    export LD_LIBRARY_PATH="${LAPACK_LIBRARY_DIR}:$LD_LIBRARY_PATH"
fi

#make sparse
cd $GOMA_LIB/sparse
if [ -e lib/libsparse.a ]
then
    echo "Sparse already built"
else
    cd src
    sed -i "/CC/c\CC=$MPI_C_COMPILER" Makefile
    make -j$MAKE_JOBS
    cd ../lib/
    cp sparse.a libsparse.a
    cd ..
fi

#make SuiteSparse
cd $GOMA_LIB/SuiteSparse
if [ -e UMFPACK/Lib/libumfpack.a ]
then
    echo "SuiteSparse is already built"
else
    cd SuiteSparse_config
    echo "compiler flag MPI is ${COMPILER_FLAG_MPI}"

    cat > SuiteSparse_config.patch << EOF
148c150
<     LAPACK ?= -llapack
---
>     LAPACK = $SUITESPARSE_NON_INTEL_LAPACK_LINK
164c166
<             BLAS = -lopenblas
---
>             BLAS = $NON_INTEL_BLAS_LINK
250c252
<     UMFPACK_CONFIG ?=
---
>     UMFPACK_CONFIG = -DNCHOLMOD
EOF


    patch SuiteSparse_config.mk < SuiteSparse_config.patch
    cd ..
    echo ${MPI_C_COMPILER}
#    continue_check
    if [ -z "${BLAS_FLAGS}" ]; then
        make static AUTOCC="no" CC="${MPI_C_COMPILER}" \
             CXX="${MPI_CXX_COMPILER}" \
             AR="${ARCHIVER}"
    else
        make static AUTOCC="no" CC="${MPI_C_COMPILER} ${COMPILER_FLAG_MPI}" \
             CXX="${MPI_CXX_COMPILER} ${COMPILER_FLAG_MPI}" \
             AR="${ARCHIVER}" BLAS_FLAGS="${BLAS_FLAGS}"
    fi
    cd ${GOMA_LIB}/SuiteSparse/UMFPACK/Include
    ln -s ../../SuiteSparse_config/SuiteSparse_config.h UFconfig.h
    ln -s ../../SuiteSparse_config/SuiteSparse_config.h SuiteSparse_config.h
fi
# Otherwise Goma dynamically links to UMFPACK when intel is sourced and disabled
export LD_LIBRARY_PATH="${GOMA_LIB}/SuiteSparse/UMFPACK/Lib:$LD_LIBRARY_PATH"

#make y12m
cd $GOMA_LIB/y12m
if [ -e liby12m.a ]
then
    echo "y12m already built"
else
# The y12m we use does not include a makefile, so this is a simple one for compiling a static library.

cat > makefile << EOF
LIB=y12m
FFLAGS=-O
FC=$MPI_F77_COMPILER
OBJ = \
	y12m.o  \
	y12mae.o    \
	y12maf.o    \
	y12mbe.o    \
	y12mbf.o    \
	y12mce.o    \
	y12mcf.o    \
	y12mde.o    \
	y12mdf.o    \
	y12mfe.o    \
	y12mge.o    \
	y12mhe.o    \
	y12cck.o
lib\$(LIB).a:	\$(OBJ)
	ar ru lib\$(LIB).a \$?
	ranlib lib\$(LIB).a
EOF

    make
fi


if [[ "$MATH_LIBRARIES" == "intel" ]] && [[ ! "$SCALAPACK_LIBRARY_NAME" = "scalapack" ]]; then
    echo "Not building scalapack because intel MKL used"
else
    # make scalapack
    cd $GOMA_LIB/scalapack-2.0.2
    if [ -f libscalapack.a ]; then
        echo "scalapack already built"
    else
        cp SLmake.inc.example SLmake.inc

cat > scalapack.patch << EOF
58,59c59,60
< BLASLIB       = -lblas
< LAPACKLIB     = -llapack
---
> BLASLIB       =  $NON_INTEL_BLAS_LINK
> LAPACKLIB     =  $LAPACK_LIBRARY_NAME_ARG
EOF

        patch SLmake.inc < scalapack.patch
        make CC="${MPI_C_COMPILER}" FC="${MPI_F90_COMPILER}" # scalapack only compiles with 1 make job
    fi
    export LD_LIBRARY_PATH="${GOMA_LIB}/scalapack-2.0.2:$LD_LIBRARY_PATH"
fi

# make mumps
cd $GOMA_LIB/MUMPS_5.1.1
if [ -e lib/libdmumps.a ]
then
    echo "MUMPS already built"
else
    cat > Makefile.inc <<EOF
# Begin orderings
#LSCOTCHDIR = /usr/lib
#ISCOTCH   = -I/usr/include/scotch # only needed for ptscotch

#LSCOTCH   = -L\$(LSCOTCHDIR) -lptesmumps -lptscotch -lptscotcherr
#LSCOTCH   = -L\$(LSCOTCHDIR) -lesmumps -lscotch -lscotcherr

LPORDDIR = \$(topdir)/PORD/lib/
IPORD    = -I\$(topdir)/PORD/include/
LPORD    = -L\$(LPORDDIR) -lpord

LMETISDIR = $GOMA_LIB/parmetis-4.0.3
IMETIS    = -I\$(LMETISDIR)/include
LMETIS    = -L\$(LMETISDIR)/lib -lparmetis -lmetis

ORDERINGSF = -Dparmetis -Dpord
ORDERINGSC  = \$(ORDERINGSF)
LORDERINGS = \$(LMETIS) \$(LPORD)
IORDERINGSF =
IORDERINGSC = \$(IMETIS) \$(IPORD)
# End orderings
################################################################################

PLAT    =
LIBEXT  = .a
OUTC    = -o
OUTF    = -o
RM = rm -f
CC = $MPI_C_COMPILER $BLAS_FLAGS $FORTRAN_LIBS
FC = $MPI_F90_COMPILER
FL = $MPI_F90_COMPILER $BLAS_FLAGS
AR = $ARCHIVER vr 
RANLIB = ranlib
SCALAP  = $SCALAPACK_LIBRARY_NAME_ARG $LAPACK_LIBRARY_NAME_ARG $BLACS_LIBRARY_NAME_ARG $NON_INTEL_BLAS_LINK


INCPAR = -I$MPI_BASE_DIR/include
LIBPAR = \$(SCALAP) -L$MPI_BASE_DIR/lib -lmpi $MPI_FORTRAN_LIB

#ifndef $OPENMPI_TOP
#    LIBPAR = $SCALAP
#endif

INCSEQ = -I\$(topdir)/libseq
LIBSEQ  =  -L\$(topdir)/libseq -lmpiseq

LIBBLAS = \$(BLAS_FLAGS) \$(NON_INTEL_BLAS_LINK)
LIBOTHERS = 

#Preprocessor defs for calling Fortran from C (-DAdd_ or -DAdd__ or -DUPPER)
CDEFS   = -DAdd_

#Begin Optimized options
OPTF    = -O  -DALLOW_NON_INIT
OPTL    = -O
OPTC    = -O
#End Optimized options
INCS = \$(INCPAR)
LIBS = \$(LIBPAR)
LIBSEQNEEDED =

EOF

    if [[ "$CC_NAME" == "intel" ]]; then
        #TODO: Find if CC_NAME or MATH_LIBRARIES affects this
        echo -e $("$MUMPS_MAKE_PATCH") > examples/Makefile.patch
        patch -f examples/Makefile < examples/Makefile.patch
        #continue_check
    fi
    make -j$MAKE_JOBS
fi

#make trilinos
rm -rf $GOMA_LIB/trilinos-12.12.1-Temp
mkdir $GOMA_LIB/trilinos-12.12.1-Temp
cd $GOMA_LIB/trilinos-12.12.1-Temp

rm -f CMakeCache.txt

MPI_LIBS="-lmpi ${MPI_FORTRAN_LIB}"

HDF5_LIBS="-L${GOMA_LIB}/hdf5-1.8.20/lib -lhdf5_hl -lhdf5 -lz -ldl"
# Install directory
TRILINOS_INSTALL=$GOMA_LIB/trilinos-12.12.1-Built

if [ -e $GOMA_LIB/trilinos-12.12.1-Built/bin/aprepro ]; then
    echo "Trilinos is already built!"
else
    cmake \
-D CMAKE_AR=/usr/bin/ar \
-D CMAKE_RANLIB=/usr/bin/ranlib \
-D CMAKE_BUILD_TYPE:STRING=RELEASE \
-D CMAKE_CXX_COMPILER:FILEPATH=${MPI_CXX_COMPILER} \
-D CMAKE_C_COMPILER:FILEPATH=${MPI_C_COMPILER} \
-D CMAKE_Fortran_COMPILER:FILEPATH=${MPI_F90_COMPILER} \
-D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
-D CMAKE_VERBOSE_CONFIGURE:BOOL=TRUE \
-D BUILD_SHARED_LIBS:BOOL=OFF \
-D Trilinos_ENABLE_Triutils:BOOL=ON \
-D Trilinos_ENABLE_SEACAS:BOOL=ON \
-D Trilinos_ENABLE_Amesos:BOOL=ON \
-D Trilinos_ENABLE_Epetra:BOOL=ON \
-D Trilinos_ENABLE_Xpetra:BOOL=OFF \
-D Trilinos_ENABLE_Ifpack:BOOL=ON \
-D Trilinos_ENABLE_Teuchos:BOOL=ON \
-D Trilinos_ENABLE_ML:BOOL=ON \
-D Trilinos_ENABLE_AztecOO:BOOL=ON \
-D Trilinos_ENABLE_Stratimikos:BOOL=ON \
-D Trilinos_ENABLE_Teko:BOOL=ON \
-D Trilinos_ENABLE_KokkosClassic:BOOL=OFF \
-D Trilinos_ENABLE_STK:BOOL=OFF \
-D Trilinos_ENABLE_Amesos2:BOOL=OFF \
-D Trilinos_ENABLE_Zoltan2:BOOL=OFF \
-D Trilinos_ENABLE_Belos:BOOL=ON \
-D Trilinos_ENABLE_Sacado:BOOL=ON \
-D Trilinos_ENABLE_EpetraExt:BOOL=ON \
-D Trilinos_ENABLE_Thyra:BOOL=ON \
-D Trilinos_ENABLE_ThyraTpetraAdapters:BOOL=ON \
-D Trilinos_ENABLE_Tpetra:BOOL=ON \
-D Trilinos_ENABLE_Stratimikos:BOOL=ON \
-D Trilinos_ENABLE_TESTS:BOOL=ON \
-D Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
-D Trilinos_ENABLE_SECONDARY_STABLE_CODE:BOOL=ON \
      -D Netcdf_LIBRARY_DIRS:PATH="$GOMA_LIB/netcdf-4.4.1.1/lib" \
      -D TPL_ENABLE_Netcdf:BOOL=ON \
      -D TPL_Netcdf_INCLUDE_DIRS:PATH="$GOMA_LIB/netcdf-4.4.1.1/include" \
      -D Matio_LIBRARY_DIRS:PATH=$GOMA_LIB/matio-1.5.10/lib \
      -D Matio_INCLUDE_DIRS:PATH=$GOMA_LIB/matio-1.5.10/include \
-D TPL_ENABLE_MPI:BOOL=ON \
  -D MPI_COMPILER:FILEPATH=$MPI_BASE_DIR/bin/mpicc \
  -D MPI_EXECUTABLE:FILEPATH=$MPI_BASE_DIR/bin/mpirun \
  -D MPI_BASE_DIR:PATH=$MPI_BASE_DIR \
-D EpetraExt_BUILD_GRAPH_REORDERINGS:BOOL=ON \
-D EpetraExt_BUILD_BDF:BOOL=ON \
-D TPL_ENABLE_Boost:BOOL=OFF \
-D TPL_ENABLE_LAPACK:BOOL=ON \
-D LAPACK_LIBRARY_DIRS="${LAPACK_LIBRARY_DIR}" \
-D LAPACK_LIBRARY_NAMES="${LAPACK_LIBRARY_NAME}" \
-D TPL_ENABLE_BLAS:BOOL=ON \
-DHAVE_EPETRA_LAPACK_GSSVD3:BOOL=ON \
-D BLAS_LIBRARY_DIRS="${BLAS_LIBRARY_DIR}" \
-D BLAS_LIBRARY_NAMES="${BLAS_LIBRARY_NAME}" \
-D CMAKE_INSTALL_PREFIX:PATH=$TRILINOS_INSTALL \
-D Trilinos_EXTRA_LINK_FLAGS:STRING="$HDF5_LIBS $MPI_LIBS $FORTRAN_LIBS -lrt" \
-D TPL_ENABLE_UMFPACK:BOOL=ON \
  -D UMFPACK_LIBRARY_NAMES:STRING="umfpack;amd;suitesparseconfig" \
  -D UMFPACK_LIBRARY_DIRS:PATH="$GOMA_LIB/SuiteSparse/UMFPACK/Lib;$GOMA_LIB/SuiteSparse/AMD/Lib;$GOMA_LIB/SuiteSparse/SuiteSparse_config" \
  -D UMFPACK_INCLUDE_DIRS:PATH="$GOMA_LIB/SuiteSparse/UMFPACK/Include;$GOMA_LIB/SuiteSparse/AMD/Include;$GOMA_LIB/SuiteSparse/SuiteSparse_config" \
-D TPL_ENABLE_AMD:BOOL=ON \
  -D AMD_LIBRARY_NAMES:STRING="amd;suitesparseconfig" \
  -D AMD_LIBRARY_DIRS:PATH="$GOMA_LIB/SuiteSparse/AMD/Lib;$GOMA_LIB/SuiteSparse/SuiteSparse_config" \
  -D AMD_INCLUDE_DIRS:PATH="$GOMA_LIB/SuiteSparse/AMD/Include;$GOMA_LIB/SuiteSparse/SuiteSparse_config" \
-D TPL_ENABLE_SuperLUDist:BOOL=ON \
  -D SuperLUDist_LIBRARY_NAMES:STRING="superludist" \
  -D SuperLUDist_LIBRARY_DIRS:PATH=$GOMA_LIB/superlu_dist-5.1.3/lib \
  -D SuperLUDist_INCLUDE_DIRS:PATH=$GOMA_LIB/superlu_dist-5.1.3/SRC \
-D TPL_ENABLE_ParMETIS:BOOL=ON \
  -D ParMETIS_LIBRARY_DIRS:PATH=$GOMA_LIB/parmetis-4.0.3/lib \
  -D ParMETIS_INCLUDE_DIRS:PATH=$GOMA_LIB/parmetis-4.0.3/include \
  -D TPL_ParMETIS_INCLUDE_DIRS:PATH=$GOMA_LIB/parmetis-4.0.3/include \
  -D TPL_ENABLE_y12m:BOOL=ON \
  -D y12m_LIBRARY_NAMES:STRING="y12m" \
  -D y12m_LIBRARY_DIRS:PATH=$GOMA_LIB/y12m \
-D TPL_ENABLE_MUMPS:BOOL=ON \
  -D MUMPS_LIBRARY_NAMES:STRING="dmumps;mumps_common;pord;$BLACS_LIBRARY_NAME" \
  -D MUMPS_LIBRARY_DIRS:PATH="$GOMA_LIB/MUMPS_5.1.1/lib;$GOMA_LIB/MUMPS_5.1.1/PORD/lib;$SCALAPACK_LIBRARY_DIR" \
  -D MUMPS_INCLUDE_DIRS:PATH="$GOMA_LIB/MUMPS_5.1.1/include;$GOMA_LIB/MUMPS_5.1.1/PORD/include" \
  -D CMAKE_CXX_FLAGS:STRING="-DMUMPS_5_0 $BLAS_FLAGS $COMPILER_FLAG_MPI" \
  -D Amesos_ENABLE_SCALAPACK:BOOL=ON \
  -D SCALAPACK_INCLUDE_DIRS:FILEPATH="${SCALAPACK_INCLUDE_DIR}" \
  -D SCALAPACK_LIBRARY_DIRS:FILEPATH="${SCALAPACK_LIBRARY_DIR}" \
  -D SCALAPACK_LIBRARY_NAMES:STRING="${SCALAPACK_LIBRARY_NAME}" \
-D TPL_ENABLE_MKL:BOOL=$USING_MKLS \
  -D MKL_LIBRARY_NAMES:STRING="${MKL_LIBRARY_NAME}" \
  -D MKL_LIBRARY_DIRS:PATH="${MKL_LIBRARY_DIR}" \
  -D MKL_INCLUDE_DIRS:PATH="${MKL_INCLUDE_DIR}" \
-D Amesos_ENABLE_SuperLUDist:BOOL=on \
-D Amesos_ENABLE_ParMETIS:BOOL=on \
-D Amesos_ENABLE_LAPACK:BOOL=ON \
-D Amesos_ENABLE_KLU:BOOL=ON \
-D Amesos_ENABLE_UMFPACK:BOOL=ON \
-D Amesos_ENABLE_MUMPS:BOOL=ON \
$EXTRA_ARGS \
$GOMA_LIB/Trilinos-trilinos-release-12-12-1

#    continue_check
    make -j$MAKE_JOBS
    make install
fi
