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

BUILD_LOG=$GOMA_LIB/goma_tpl_build.log
COMPILE_LOG=$GOMA_LIB/goma_tpl_compile.log

function log_echo() {
    builtin echo "$@" 2>&1 | tee -a $BUILD_LOG
}

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
        elif [ -d "/usr/lib/x86_64-linux-gnu/openmpi" ] ; then
            MPI_BASE_DIR="/usr/lib/x86_64-linux-gnu/openmpi"
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

echo "Start Goma Build" >> $BUILD_LOG
echo "Start Goma Compile" >> $COMPILE_LOG

HDF5_VERSION="1.10.6"
HDF5_MD5="03095102a6118c32a75a9b9b40be66f2"

NETCDF_VERSION="c-4.7.4"
NETCDF_MD5="3e0a97e6abb9a989f8a8a2e395473597"

TRILINOS_VERSION="12.18.1"
TRILINOS_VERSION_DASH="12-18-1"
TRILINOS_MD5="9c1d151169949bca6cf203831e4d6aee"

MUMPS_VERSION="5.3.1"
MUMPS_MD5="31b64a11c1df6a56b1750411efb20986"

OPENMPI_VERSION="4.0.3"
OPENMPI_MD5="851553085013939f24cdceb1af06b828"
OPENMPI_ARCHIVE_URL="https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-$OPENMPI_VERSION.tar.bz2"
OPENMPI_EXTRA_CONFIGURE_FLAGS=""

CMAKE_VERSION="3.17.1"
CMAKE_MD5="958959aa5e0338144eed7320e9b48561"

SUITESPARSE_VERSION="4.5.6"
SUITESPARSE_MD5="eeb87a842a9b3b0425cf08d97fb3c5ec"

MATIO_VERSION="1.5.17"
MATIO_MD5="170075cce5c144e19f610af9b64cb63b"

SCALAPACK_VERSION="2.1.0"
SCALAPACK_MD5="3b239ef80353c67354a0a62d54946fa8"

ARCHIVE_NAMES=("arpack96.tar.gz" \
"patch.tar.gz" \
"hdf5-${HDF5_VERSION}.tar.bz2" \
"netcdf-${NETCDF_VERSION}.tar.gz" \
"parmetis-4.0.3.tar.gz" \
"sparse.tar.gz" \
"superlu_dist-5.1.3.tar.gz" \
"y12m.tar.gz" \
"Trilinos-trilinos-release-$TRILINOS_VERSION_DASH.tar.gz" \
"MUMPS_$MUMPS_VERSION.tar.gz" \
"SuiteSparse-$SUITESPARSE_VERSION.tar.gz" \
"matio-$MATIO_VERSION.tar.gz")

#y12m archive is skipped because it stores the number of downloads in the header

#meaning each y12m tar has a unique MD5SUM.
ARCHIVE_MD5SUMS=("fffaa970198b285676f4156cebc8626e" \
"14830d758f195f272b8594a493501fa2" \
"${HDF5_MD5}" \
"${NETCDF_MD5}" \
"f69c479586bf6bb7aff6a9bc0c739628" \
"1566d914d1035ac17b73fe9bc0eed02a" \
"fdee368cba0e95cb0143b6d47915e7a1" \
"SKIP" \
$TRILINOS_MD5 \
$MUMPS_MD5 \
"$SUITESPARSE_MD5" \
"$MATIO_MD5")

ARCHIVE_URLS=("http://www.caam.rice.edu/software/ARPACK/SRC/arpack96.tar.gz" \
"http://www.caam.rice.edu/software/ARPACK/SRC/patch.tar.gz" \
"https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-${HDF5_VERSION}/src/hdf5-${HDF5_VERSION}.tar.bz2" \
"https://www.unidata.ucar.edu/downloads/netcdf/ftp/netcdf-${NETCDF_VERSION}.tar.gz" \
"http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz" \
"http://downloads.sourceforge.net/project/sparse/sparse/sparse1.4b/sparse1.4b.tar.gz" \
"http://codeload.github.com/xiaoyeli/superlu_dist/tar.gz/v5.1.3" \
"http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz\\&filename=y12m%2Fy12m.f" \
"https://github.com/trilinos/Trilinos/archive/trilinos-release-$TRILINOS_VERSION_DASH.tar.gz" \
"http://mumps.enseeiht.fr/MUMPS_$MUMPS_VERSION.tar.gz" \
"http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-$SUITESPARSE_VERSION.tar.gz" \
"https://github.com/tbeu/matio/releases/download/v$MATIO_VERSION/matio-$MATIO_VERSION.tar.gz")

# You can't call the ARPACK patch ARPACK or it will think it is already extracted
# When in reality it isn't
ARCHIVE_DIR_NAMES=("ARPACK" \
"FAKE_DIR_FOR_ARPACK_PATCH" \
"hdf5-${HDF5_VERSION}" \
"netcdf-${NETCDF_VERSION}" \
"parmetis-4.0.3" \
"sparse" \
"superlu_dist-5.1.3" \
"y12m" \
"Trilinos-trilinos-release-$TRILINOS_VERSION_DASH" \
"MUMPS_$MUMPS_VERSION" \
"SuiteSparse" \
"matio-$MATIO_VERSION")

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
    ARCHIVE_NAMES+=("openmpi-$OPENMPI_VERSION.tar.bz2")
    ARCHIVE_MD5SUMS+=("$OPENMPI_MD5")
    ARCHIVE_URLS+=("$OPENMPI_ARCHIVE_URL")
    ARCHIVE_DIR_NAMES+=("openmpi-$OPENMPI_VERSION")
    ARCHIVE_HOMEPAGES+=("https://www.open-mpi.org/")
    ARCHIVE_REAL_NAMES+=("OpenMPI")
#    If the user says to use the intel library but sets the base dir to BUILD then there's a conflict here
    MPI_NAME="open"
    MPI_BASE_DIR="$GOMA_LIB/openmpi-$OPENMPI_VERSION"
    MPI_IS_BUILT_FROM_SCRATCH="true"
    echo
    echo "OpenMPI $OPENMPI_VERSION will be built."
    echo "You will have to add ${MPI_BASE_DIR}/lib to your LD_LIBRARY_PATH"
    echo "And ${MPI_BASE_DIR}/bin to your PATH environment variables"
    echo "Before running Goma"
    echo
    continue_check
fi

if command -v cmake; then
    cmake_vers=$(cmake --version |grep "version" | awk '{print $NF}')
    if [[ "$cmake_vers" = $(echo -e "$cmake_vers\n3.10.0\n" | sort -V |tail -n1) ]]; then
	build_cmake="false"
    else
	build_cmake="true"
    fi
else
    build_cmake="true"
fi

if [ "$build_cmake" == "false" ] ; then
    log_echo "Native cmake found newer than 3.10.0, skipping download"
else
    ARCHIVE_NAMES+=("cmake-$CMAKE_VERSION.tar.gz")
    ARCHIVE_MD5SUMS+=("$CMAKE_MD5")
    ARCHIVE_URLS+=("https://github.com/Kitware/CMake/releases/download/v$CMAKE_VERSION/cmake-$CMAKE_VERSION.tar.gz")
    ARCHIVE_DIR_NAMES+=("cmake-$CMAKE_VERSION")
    ARCHIVE_HOMEPAGES+=("https://cmake.org/")
    ARCHIVE_REAL_NAMES+=("CMake")
    log_echo "Cmake not found, will build."
    log_echo
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
                export SCALAPACK_LIBRARY_DIR="${GOMA_LIB}/scalapack-$SCALAPACK_VERSION/lib"
                export SCALAPACK_LIBRARY_NAME="scalapack"
                export SCALAPACK_LIBRARY_NAME_ARG="-L${SCALAPACK_LIBRARY_DIR} -lscalapack"
                SCALAPACK_INCLUDE_DIR="${GOMA_LIB}/scalapack-$SCALAPACK_VERSION/include"
                ARCHIVE_NAMES+=("scalapack-$SCALAPACK_VERSION.tgz")
                ARCHIVE_MD5SUMS+=("$SCALAPACK_MD5")
                ARCHIVE_URLS+=("http://www.netlib.org/scalapack/scalapack-$SCALAPACK_VERSION.tgz")
                ARCHIVE_DIR_NAMES+=("scalapack-$SCALAPACK_VERSION")
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
        "scalapack-$SCALAPACK_VERSION.tgz")
        ARCHIVE_MD5SUMS+=("cd132aea6f7055a49aa48ca0a61e7cd5" \
        "96591affdbf58c450d45c1daa540dbd2" \
        "$SCALAPACK_MD5" )
        ARCHIVE_URLS+=("http://www.netlib.org/blas/blas-3.7.1.tgz" \
        "http://www.netlib.org/lapack/lapack-3.8.0.tar.gz" \
        "http://www.netlib.org/scalapack/scalapack-$SCALAPACK_VERSION.tgz" )
        ARCHIVE_DIR_NAMES+=("BLAS-3.7.1" \
        "lapack-3.8.0" \
        "scalapack-$SCALAPACK_VERSION" )
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
        SCALAPACK_LIBRARY_DIR="${GOMA_LIB}/scalapack-$SCALAPACK_VERSION/lib"
        SCALAPACK_LIBRARY_NAME="scalapack"
        SCALAPACK_LIBRARY_NAME_ARG="-L${SCALAPACK_LIBRARY_DIR} -lscalapack"
        SCALAPACK_INCLUDE_DIR="${GOMA_LIB}/scalapack-$SCALAPACK_VERSION/include"
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
        export SCALAPACK_LIBRARY_DIR="${GOMA_LIB}/scalapack-$SCALAPACK_VERSION/lib"
        export SCALAPACK_LIBRARY_NAME_ARG="-L{SCALAPACK_LIBRARY_DIR} -lscalapack"
        SCALAPACK_INCLUDE_DIR="${GOMA_LIB}/scalapack-$SCALAPACK_VERSION/include"
        export NON_INTEL_BLAS_LIBRARY="${BLAS_LIBRARY_DIR}/lib${BLAS_LIBRARY_NAME}.a"
        export NON_INTEL_BLAS_LINK="-L${BLAS_LIBRARY_DIR} -l${BLAS_LIBRARY_NAME} ${FORTRAN_LIBS}"
        ARCHIVE_NAMES+=("scalapack-$SCALAPACK_VERSION.tgz")
        ARCHIVE_MD5SUMS+=("$SCALAPACK_MD5" )
        ARCHIVE_URLS+=("http://www.netlib.org/scalapack/scalapack-$SCALAPACK_VERSION.tgz" )
        ARCHIVE_DIR_NAMES+=("scalapack-$SCALAPACK_VERSION" )
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


cd $GOMA_LIB


function mychecksum {
    local count=$2
    local archive=$1
    local MD5SAVED=${ARCHIVE_MD5SUMS[count]}
    if [ ! "$MD5SAVED" == "SKIP" ]; then
        local MD5ARCHIVE=($(md5sum $archive))
        if [ $MD5SAVED != $MD5ARCHIVE ]; then
            log_echo "Issue checksum with archive:"
            log_echo "$MD5SAVED"
            log_echo "$MD5ARCHIVE"
            log_echo $archive
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
    log_echo "Check for $i at ${ARCHIVE_URLS[count]}"
    if ! [ -f $i ]
    then
        wget "${ARCHIVE_URLS[count]}" -O $i
        mychecksum $i $count
    else
        mychecksum $i $count
    fi
    cd ..

    if [ -d ${ARCHIVE_DIR_NAMES[count]} ]
    then
	log_echo "already extracted ${i}"
    else
	if ! tar tf $i &> /dev/null; then
	    tar -xf tars/$i
	fi
    fi
    count=$(( $count + 1 ))
    cd tars
done

export CXX=${SYSTEM_CXX}

if [ "$build_cmake" == "false" ] ; then
    log_echo "Native cmake found of sufficient version, skipping build"
else
    cd $GOMA_LIB/cmake-$CMAKE_VERSION
    if [ -f bin/cmake ]
    then
        log_echo "cmake is already built"
    else
        CC=$SYSTEM_CC CXX=$SYSTEM_CXX FC=$SYSTEM_FC ./bootstrap --prefix=$GOMA_LIB/cmake-$CMAKE_VERSION 2>&1 | tee -a $COMPILE_LOG
        make -j$MAKE_JOBS 2>&1 | tee -a $COMPILE_LOG
        make install 2>&1 | tee -a $COMPILE_LOG
        if [ -f $GOMA_LIB/cmake-$CMAKE_VERSION/bin/cmake ]
        then
            log_echo "Built CMake version $CMAKE_VERSION"
        else
            log_echo "Failed to build CMake version $CMAKE_VERSION"
            exit 1
        fi
    fi
    export PATH=$GOMA_LIB/cmake-$CMAKE_VERSION/bin:$PATH
fi


if [[ "$MPI_IS_BUILT_FROM_SCRATCH" == "true" ]]; then
    #continue_check
    #make openmpi
    cd $GOMA_LIB/openmpi-$OPENMPI_VERSION
    if [ -f bin/ompi_info ]
    then
        log_echo "openmpi is already built"
    else
        CC=${SYSTEM_CC} CXX=${SYSTEM_CXX} FC=${SYSTEM_FC} CPP=${SYSTEM_CPP} ./configure $OPENMPI_EXTRA_CONFIGURE_FLAGS --prefix=$GOMA_LIB/openmpi-$OPENMPI_VERSION 2>&1 | tee -a $COMPILE_LOG
        make -j$MAKE_JOBS CC=${SYSTEM_CC} CXX=${SYSTEM_CXX} FC=${SYSTEM_FC} CPP=${SYSTEM_CPP} 2>&1 | tee -a $COMPILE_LOG 
        make install 2>&1 | tee -a $COMPILE_LOG
        if [ -f $GOMA_LIB/openmpi-$OPENMPI_VERSION/bin/ompi_info ]
        then
            log_echo "Built OpenMPI $OPENMPI_VERSION"
        else
            log_echo "Failed to build OpenMPI $OPENMPI_VERSION"
        fi
    fi
    cd ..
else
    log_echo "Using custom mpi: $mpi"
fi
cd $GOMA_LIB

#hdf5
if [ -e hdf5-${HDF5_VERSION}/lib/libhdf5.a ]
then
    log_echo "hdf5 already built"
else
    if ! [ -e hdf5-${HDF5_VERSION}/.goma-extracted ]
    then
	mv hdf5-${HDF5_VERSION} tmpdir
	mkdir hdf5-${HDF5_VERSION}
	mv tmpdir hdf5-${HDF5_VERSION}/src
	touch hdf5-${HDF5_VERSION}/.goma-extracted
    fi

    cd hdf5-${HDF5_VERSION}/src
    CC="$MPI_C_COMPILER" CPP="$MPI_C_COMPILER -E" AR=${ARCHIVER} ./configure --enable-shared=off --prefix=$GOMA_LIB/hdf5-${HDF5_VERSION} --enable-parallel 2>&1 | tee -a $COMPILE_LOG
    make -j$MAKE_JOBS 2>&1 | tee -a $COMPILE_LOG
    make install 2>&1 | tee -a $COMPILE_LOG
    if [ -e $GOMA_LIB/hdf5-${HDF5_VERSION}/lib/libhdf5.a ]
    then
        log_echo "Built HDF5 $HDF5_VERSION"
    else
        log_echo "Failed to build HDF5 $HDF5_VERSION"
        exit 1
    fi
    cd ../..
fi
cd $GOMA_LIB

#matio
if [ -e matio-$MATIO_VERSION/lib/libmatio.a ]
then
    log_echo "matio already built"
else
    if ! [ -e matio-$MATIO_VERSION/.goma-extracted ]
    then
	mv matio-$MATIO_VERSION tmpdir
	mkdir matio-$MATIO_VERSION
	mv tmpdir matio-$MATIO_VERSION/src
	touch matio-$MATIO_VERSION/.goma-extracted
    fi
    cd matio-$MATIO_VERSION/src
    CC=${MPI_C_COMPILER} LD=${MPI_CXX_COMPILER} AR=${ARCHIVER} LIBS="-ldl" ./configure --with-hdf5=${GOMA_LIB}/hdf5-${HDF5_VERSION} --prefix=${GOMA_LIB}/matio-$MATIO_VERSION --enable-shared=off 2>&1 | tee -a $COMPILE_LOG
    make -j$MAKE_JOBS 2>&1 | tee -a $COMPILE_LOG
    make install 2>&1 | tee -a $COMPILE_LOG
    if [ -e $GOMA_LIB/matio-$MATIO_VERSION/lib/libmatio.a ]
    then
        log_echo "Built MATIO $MATIO_VERSION"
    else
        log_echo "Failed to build MATIO $MATIO_VERSION"
        exit 1
    fi
fi
cd $GOMA_LIB

#netcdf
if [ -e netcdf-${NETCDF_VERSION}/lib/libnetcdf.a ]
then
    log_echo "netcdf already built"
else
    if ! [ -e netcdf-${NETCDF_VERSION}/.goma-extracted ]
    then
	mv netcdf-${NETCDF_VERSION} tmpdir
	mkdir netcdf-${NETCDF_VERSION}
	mv tmpdir netcdf-${NETCDF_VERSION}/src
	touch netcdf-${NETCDF_VERSION}/.goma-extracted
    fi
    cd $GOMA_LIB/netcdf-${NETCDF_VERSION}/src
    export CPPFLAGS=-I$GOMA_LIB/hdf5-${HDF5_VERSION}/include
    log_echo $CPPFLAGS
    log_echo $LDFLAGS

    CC=${MPI_C_COMPILER} CFLAGS="-I${GOMA_LIB}/hdf5-${HDF5_VERSION}/include" \
      CPP="${MPI_C_COMPILER} -E" \
      CPPFLAGS="-I${GOMA_LIB}/hdf5-${HDF5_VERSION}/include" \
      LDFLAGS="-L${GOMA_LIB}/hdf5-${HDF5_VERSION}/lib" \
      ./configure \
      --prefix=$GOMA_LIB/netcdf-${NETCDF_VERSION} \
      --enable-shared=off \
      --disable-dap 2>&1 | tee -a $COMPILE_LOG

    make -j$MAKE_JOBS 2>&1 | tee -a $COMPILE_LOG
    make install 2>&1 | tee -a $COMPILE_LOG
    if [ -e $GOMA_LIB/netcdf-${NETCDF_VERSION}/lib/libnetcdf.a ]
    then
        log_echo "Built NetCDF $NETCDF_VERSION"
    else
        log_echo "Failed to build NetCDF $NETCDF_VERSION"
        exit 1
    fi
    cd ../..
fi

 #make BLAS
if [[ "$MATH_LIBRARIES" == "netlib blas" ]]; then
    cd $GOMA_LIB/BLAS-3.7.1
    if [ -f libblas.a ]
    then
        log_echo "BLAS already built"
    else
        log_echo $("$BLAS_PATCH") > makepatch.inc
        patch -f make.inc < makepatch.inc
        make -j$MAKE_JOBS cc=${MPI_C_COMPILER} ccflags=${COMPILER_FLAG_MPI} cxx=${MPI_CXX_COMPILER} cxxflags=${COMPILER_FLAG_MPI} 2>&1 | tee -a $COMPILE_LOG
        cp blas_LINUX.a libblas.a
    fi
    if [ -f $GOMA_LIB/BLAS-3.7.1/libblas.a ]
    then
        log_echo "Built BLAS 3.7.1"
    else
        log_echo "Failed to build BLAS 3.7.1"
        exit 1
    fi
    # For some reason this junk apple file is included with BLAS 3.7.1 from netlib.
    rm -f "$GOMA_LIB/._BLAS-3.7.1"
    export LD_LIBRARY_PATH="${GOMA_LIB}/BLAS-3.7.1:$LD_LIBRARY_PATH"
fi

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
    log_echo "ParMetis already Built"
else
    if [ ! -e parmetis-4.0.3/src/Makefile ]; then
        mv parmetis-4.0.3 tmpdir
        mkdir parmetis-4.0.3
        mv tmpdir parmetis-4.0.3/src
    fi
    cd parmetis-4.0.3/src
    log_echo "$PARMETIS_PATCH" > parmetis_make.patch
    patch -f Makefile < parmetis_make.patch
    log_echo "make config parmetis?"
    make config cc=${MPI_C_COMPILER} ccflags=${COMPILER_FLAG_MPI} cxx=${MPI_CXX_COMPILER} cxxflags=${COMPILER_FLAG_MPI} prefix=$GOMA_LIB/parmetis-4.0.3 2>&1 | tee -a $COMPILE_LOG
    make 2>&1 | tee -a $COMPILE_LOG
    make install 2>&1 | tee -a $COMPILE_LOG
    cd ..
    if [ -d include ]; then
	cp src/metis/include/metis.h include
    fi
    if [ -d lib ]; then
	cp src/build/Linux-x86_64/libmetis/libmetis.a lib/
    fi
    if [ -e $GOMA_LIB/parmetis-4.0.3/lib/libparmetis.a ]
    then
        log_echo "Built ParMetis 4.0.3"
    else
        log_echo "Failed to build ParMetis 4.0.3"
        exit 1
    fi
fi


#make ARPACK
cd $GOMA_LIB/ARPACK
if [ -e lib/libarpack.a ]
then
    log_echo "ARPACK already built"
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
    cat > UTIL/second.f.patch << EOF
24c24
<       EXTERNAL           ETIME
---
> *     EXTERNAL           ETIME
EOF
    patch UTIL/second.f < UTIL/second.f.patch
    make all 2>&1 | tee -a $COMPILE_LOG
    mkdir lib
    cp libarpack_x86_64.a lib/libarpack.a
    if [ -e $GOMA_LIB/ARPACK/lib/libarpack.a ]
    then
        log_echo "Built ARPACK"
    else
        log_echo "Failed to build ARPACK"
        exit 1
    fi
fi

#make SuperLU
cd $GOMA_LIB/superlu_dist-5.1.3
if [ -e lib/libsuperludist.a ]
then
    log_echo "SuperLU_DIST already built"
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
    make 2>&1 | tee -a $COMPILE_LOG
    if [ -e $GOMA_LIB/superlu_dist-5.1.3/lib/libsuperludist.a ]
    then
        log_echo "Built SuperLU_DIST 5.1.3"
    else
        log_echo "Failed to build SuperLU_DIST 5.1.3"
        exit 1
    fi
fi

if [[ "$MATH_LIBRARIES" == "netlib blas" ]]; then
    #make lapack
    cd $LAPACK_LIBRARY_DIR
    if [ -e liblapack.a ]
    then
        log_echo "LAPACK already built"
    else
        mv make.inc.example make.inc
        log_echo $("$LAPACK_PATCH") > make.patch
        patch make.inc < make.patch
        make lapacklib -j$MAKE_JOBS cc=${MPI_C_COMPILER} ccflags=${COMPILER_FLAG_MPI} cxx=${MPI_CXX_COMPILER} cxxflags=${COMPILER_FLAG_MPI} 2>&1 | tee -a $COMPILE_LOG
#        cp lapack_LINUX.a liblapack.a
        if [ -e $LAPACK_LIBRARY_DIR/liblapack.a ]
        then
            log_echo "Built LAPACK"
        else
            log_echo "Failed to build LAPACK"
            exit 1
        fi
    fi
    export LD_LIBRARY_PATH="${LAPACK_LIBRARY_DIR}:$LD_LIBRARY_PATH"
fi

#make sparse
cd $GOMA_LIB/sparse
if [ -e lib/libsparse.a ]
then
    log_echo "Sparse already built"
else
    cd src
    sed -i "/CC/c\CC=$MPI_C_COMPILER" Makefile
    make -j$MAKE_JOBS 2>&1 | tee -a $COMPILE_LOG
    cd ../lib/
    cp sparse.a libsparse.a
    cd ..
    if [ -e $GOMA_LIB/sparse/lib/libsparse.a ]
    then
        log_echo "Built sparse"
    else
        log_echo "Failed to build sparse"
        exit 1
    fi
fi

#make SuiteSparse
cd $GOMA_LIB/SuiteSparse
if [ -e UMFPACK/Lib/libumfpack.a ]
then
    log_echo "SuiteSparse is already built"
else
    cd SuiteSparse_config
    log_echo "compiler flag MPI is ${COMPILER_FLAG_MPI}"

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
    log_echo ${MPI_C_COMPILER}
    if [ -z "${BLAS_FLAGS}" ]; then
        make static AUTOCC="no" CC="${MPI_C_COMPILER}" \
             CXX="${MPI_CXX_COMPILER}" \
             AR="${ARCHIVER}" 2>&1 | tee -a $COMPILE_LOG
    else
        make static AUTOCC="no" CC="${MPI_C_COMPILER} ${COMPILER_FLAG_MPI}" \
             CXX="${MPI_CXX_COMPILER} ${COMPILER_FLAG_MPI}" \
             AR="${ARCHIVER}" BLAS_FLAGS="${BLAS_FLAGS}" 2>&1 | tee -a $COMPILE_LOG
    fi
    cd ${GOMA_LIB}/SuiteSparse/UMFPACK/Include
    ln -s ../../SuiteSparse_config/SuiteSparse_config.h UFconfig.h
    ln -s ../../SuiteSparse_config/SuiteSparse_config.h SuiteSparse_config.h
    if [ -e $GOMA_LIB/SuiteSparse/UMFPACK/Lib/libumfpack.a ]
    then
        log_echo "Built SuiteSparse $SUITESPARSE_VERSION"
    else
        log_echo "Failed to build SuiteSparse $SUITESPARSE_VERSION"
        exit 1
    fi
fi
# Otherwise Goma dynamically links to UMFPACK when intel is sourced and disabled
export LD_LIBRARY_PATH="${GOMA_LIB}/SuiteSparse/UMFPACK/Lib:$LD_LIBRARY_PATH"

#make y12m
cd $GOMA_LIB/y12m
if [ -e liby12m.a ]
then
    log_echo "y12m already built"
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

    make 2>&1 | tee -a $COMPILE_LOG
    if [ -e $GOMA_LIB/y12m/liby12m.a ]
    then
        log_echo "Built y12m"
    else
        log_echo "Failed to build y12m"
        exit 1
    fi
fi


if [[ "$MATH_LIBRARIES" == "intel" ]] && [[ ! "$SCALAPACK_LIBRARY_NAME" = "scalapack" ]]; then
    log_echo "Not building scalapack because intel MKL used"
else
    # make scalapack
    if [ -f $GOMA_LIB/scalapack-$SCALAPACK_VERSION/lib/libscalapack.a ]; then
        log_echo "scalapack already built"
    else
        mv $GOMA_LIB/scalapack-$SCALAPACK_VERSION src-scalapack
        mkdir $GOMA_LIB/scalapack-$SCALAPACK_VERSION
        mv src-scalapack $GOMA_LIB/scalapack-$SCALAPACK_VERSION/src
        cd $GOMA_LIB/scalapack-$SCALAPACK_VERSION/src
        mkdir build
        cd build
        cmake .. -DCMAKE_INSTALL_PREFIX=$GOMA_LIB/scalapack-$SCALAPACK_VERSION -DBUILD_SHARED=OFF -DCMAKE_C_COMPILER=$MPI_C_COMPILER -DCMAKE_CXX_COMPILER=$MPI_CXX_COMPILER -DCMAKE_Fortran_Compiler=$MPI_F90_COMPILER 2>&1 | tee -a $COMPILE_LOG
        make -j$MAKE_JOBS 2>&1 | tee -a $COMPILE_LOG
        make install 2>&1 | tee -a $COMPILE_LOG
        mkdir -p $GOMA_LIB/scalapack-$SCALAPACK_VERSION/include
        cd $GOMA_LIB
        if [ -f $GOMA_LIB/scalapack-$SCALAPACK_VERSION/lib/libscalapack.a ]; then
            log_echo "Build scalapack $SCALAPACK_VERSION"
        else
            log_echo "Error building scalapack $SCALAPACK_VERSION"
            exit 1
        fi
    fi
    export LD_LIBRARY_PATH="${GOMA_LIB}/scalapack-$SCALAPACK_VERSION/lib:$LD_LIBRARY_PATH"
fi

# make mumps
cd $GOMA_LIB/MUMPS_$MUMPS_VERSION
if [ -e lib/libdmumps.a ]
then
    log_echo "MUMPS already built"
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
OPTF    = -O2  -DALLOW_NON_INIT
OPTL    = -O2
OPTC    = -O2
#End Optimized options
INCS = \$(INCPAR)
LIBS = \$(LIBPAR)
LIBSEQNEEDED =

EOF

    if [[ "$CC_NAME" == "intel" ]]; then
        #TODO: Find if CC_NAME or MATH_LIBRARIES affects this
        log_echo -e $("$MUMPS_MAKE_PATCH") > examples/Makefile.patch
        patch -f examples/Makefile < examples/Makefile.patch
    fi
    make -j$MAKE_JOBS 2>&1 | tee -a $COMPILE_LOG
    if [ -e $GOMA_LIB/MUMPS_$MUMPS_VERSION/lib/libdmumps.a ]
    then
        log_echo "Built MUMPS $MUMPS_VERSION"
    else
        log_echo "Failed to build MUMP $MUMPS_VERSION"
        exit 1
    fi
fi

#make trilinos
rm -rf $GOMA_LIB/trilinos-$TRILINOS_VERSION-Temp
mkdir $GOMA_LIB/trilinos-$TRILINOS_VERSION-Temp
cd $GOMA_LIB/trilinos-$TRILINOS_VERSION-Temp

rm -f CMakeCache.txt

MPI_LIBS="-lmpi ${MPI_FORTRAN_LIB}"

HDF5_LIBS="-L${GOMA_LIB}/hdf5-${HDF5_VERSION}/lib -lhdf5_hl -lhdf5 -lz -ldl"
# Install directory
TRILINOS_INSTALL=$GOMA_LIB/trilinos-$TRILINOS_VERSION

if [ -e $TRILINOS_INSTALL/bin/aprepro ]; then
    log_echo "Trilinos is already built!"
else
    cmake \
-D CMAKE_AR=/usr/bin/ar \
-D CMAKE_RANLIB=/usr/bin/ranlib \
-D CMAKE_BUILD_TYPE:STRING=RELEASE \
-D CMAKE_CXX_COMPILER:FILEPATH=${MPI_CXX_COMPILER} \
-D CMAKE_C_COMPILER:FILEPATH=${MPI_C_COMPILER} \
-D CMAKE_Fortran_COMPILER:FILEPATH=${MPI_F90_COMPILER} \
-D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
-D BUILD_SHARED_LIBS:BOOL=OFF \
-D TPL_ENABLE_Boost:BOOL=OFF \
-D Trilinos_ENABLE_ShyLU:BOOL=OFF \
-D Trilinos_ENABLE_ShyLU_NodeTacho:BOOL=OFF \
-D Trilinos_ENABLE_ShyLU_NodeBasker:BOOL=OFF \
-D Trilinos_ENABLE_Triutils:BOOL=ON \
-D Trilinos_ENABLE_SEACAS:BOOL=ON \
-D Trilinos_ENABLE_Amesos:BOOL=ON \
-D Trilinos_ENABLE_Epetra:BOOL=ON \
-D Trilinos_ENABLE_Xpetra:BOOL=ON \
-D Trilinos_ENABLE_Ifpack:BOOL=ON \
-D Trilinos_ENABLE_Ifpack2:BOOL=ON \
-D Trilinos_ENABLE_Teuchos:BOOL=ON \
-D Trilinos_ENABLE_ML:BOOL=ON \
-D Trilinos_ENABLE_MueLu:BOOL=ON \
-D Trilinos_ENABLE_AztecOO:BOOL=ON \
-D Trilinos_ENABLE_Stratimikos:BOOL=ON \
-D Trilinos_ENABLE_Teko:BOOL=ON \
-D Trilinos_ENABLE_Amesos2:BOOL=ON \
-D Trilinos_ENABLE_Belos:BOOL=ON \
-D Trilinos_ENABLE_Sacado:BOOL=ON \
-D Trilinos_ENABLE_EpetraExt:BOOL=ON \
-D Trilinos_ENABLE_Thyra:BOOL=ON \
-D Trilinos_ENABLE_ThyraTpetraAdapters:BOOL=ON \
-D Trilinos_ENABLE_Tpetra:BOOL=ON \
-D Trilinos_ENABLE_Stratimikos:BOOL=ON \
-D Trilinos_ENABLE_TESTS:BOOL=ON \
-D Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
      -D HDF5_LIBRARY_DIRS:PATH="$GOMA_LIB/hdf5-${HDF5_VERSION}/lib" \
      -D TPL_ENABLE_HDF5:BOOL=ON \
      -D TPL_HDF5_INCLUDE_DIRS:PATH="$GOMA_LIB/hdf5-${HDF5_VERSION}/include" \
      -D HDF5_LIBRARY_DIRS:PATH="$GOMA_LIB/hdf5-${HDF5_VERSION}/lib" \
      -D HDF5_LIBRARY_NAMES:STRING="hdf5_hl;hdf5;z;dl" \
      -D Netcdf_LIBRARY_DIRS:PATH="$GOMA_LIB/netcdf-${NETCDF_VERSION}/lib" \
      -D TPL_ENABLE_Netcdf:BOOL=ON \
      -D TPL_Netcdf_INCLUDE_DIRS:PATH="$GOMA_LIB/netcdf-${NETCDF_VERSION}/include" \
      -D Matio_LIBRARY_DIRS:PATH=$GOMA_LIB/matio-$MATIO_VERSION/lib \
      -D Matio_INCLUDE_DIRS:PATH=$GOMA_LIB/matio-$MATIO_VERSION/include \
-D TPL_ENABLE_MPI:BOOL=ON \
-D MPI_BASE_DIR:PATH=$MPI_BASE_DIR \
-D EpetraExt_BUILD_GRAPH_REORDERINGS:BOOL=ON \
-D TPL_ENABLE_LAPACK:BOOL=ON \
-D LAPACK_LIBRARY_DIRS="${LAPACK_LIBRARY_DIR}" \
-D LAPACK_LIBRARY_NAMES="${LAPACK_LIBRARY_NAME}" \
-D TPL_ENABLE_BLAS:BOOL=ON \
-DHAVE_EPETRA_LAPACK_GSSVD3:BOOL=ON \
-D BLAS_LIBRARY_DIRS="${BLAS_LIBRARY_DIR}" \
-D BLAS_LIBRARY_NAMES="${BLAS_LIBRARY_NAME}" \
-D CMAKE_INSTALL_PREFIX:PATH=$TRILINOS_INSTALL \
-D Trilinos_EXTRA_LINK_FLAGS:STRING="$HDF5_LIBS $MPI_LIBS $LAPACK_LIBRARY_NAME_ARG $BLAS_LIBRARY_NAME_ARG $FORTRAN_LIBS -lrt -lm" \
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
  -D MUMPS_LIBRARY_DIRS:PATH="$GOMA_LIB/MUMPS_$MUMPS_VERSION/lib;$GOMA_LIB/MUMPS_$MUMPS_VERSION/PORD/lib;$SCALAPACK_LIBRARY_DIR" \
  -D MUMPS_INCLUDE_DIRS:PATH="$GOMA_LIB/MUMPS_$MUMPS_VERSION/include;$GOMA_LIB/MUMPS_$MUMPS_VERSION/PORD/include" \
  -D CMAKE_CXX_FLAGS:STRING="-DMUMPS_5_0 $BLAS_FLAGS $COMPILER_FLAG_MPI" \
  -D Amesos_ENABLE_SCALAPACK:BOOL=ON \
  -D SCALAPACK_INCLUDE_DIRS:FILEPATH="${SCALAPACK_INCLUDE_DIR}" \
  -D SCALAPACK_LIBRARY_DIRS:FILEPATH="${SCALAPACK_LIBRARY_DIR}" \
  -D SCALAPACK_LIBRARY_NAMES:STRING="${SCALAPACK_LIBRARY_NAME}" \
-D TPL_ENABLE_MKL:BOOL=$USING_MKLS \
  -D MKL_LIBRARY_NAMES:STRING="${MKL_LIBRARY_NAME}" \
  -D MKL_LIBRARY_DIRS:PATH="${MKL_LIBRARY_DIR}" \
  -D MKL_INCLUDE_DIRS:PATH="${MKL_INCLUDE_DIR}" \
-D Amesos_ENABLE_SuperLUDist:BOOL=ON \
-D Amesos_ENABLE_ParMETIS:BOOL=ON \
-D Amesos_ENABLE_LAPACK:BOOL=ON \
-D Amesos_ENABLE_KLU:BOOL=ON \
-D Amesos_ENABLE_UMFPACK:BOOL=ON \
-D Amesos_ENABLE_MUMPS:BOOL=ON \
$EXTRA_ARGS \
$GOMA_LIB/Trilinos-trilinos-release-$TRILINOS_VERSION_DASH 2>&1 | tee -a $COMPILE_LOG

    make -j$MAKE_JOBS 2>&1 | tee -a $COMPILE_LOG
    make install 2>&1 | tee -a $COMPILE_LOG
    if [ -e $TRILINOS_INSTALL/bin/aprepro ]; then
        log_echo "Built Trilinos $TRILINOS_VERSION"
    else
        log_echo "Failed to build Trilinos $TRILINOS_VERSION"
        exit 1
    fi
fi
