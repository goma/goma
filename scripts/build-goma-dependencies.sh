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
MPI_CXX_COMPILER="mpicxx"
# mpif90 and mpif77 are depricated in favor of just mpif90.
MPI_F90_COMPILER="mpif90"
MPI_F77_COMPILER="mpif77"
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

if command -v readlink &> /dev/null; then
    export GOMA_LIB=`readlink --canonicalize-missing $1`
else
    GOMA_LIB=$1
    echo "WARNING: readlink not found make sure GOMA_LIB is a full path"
    echo "Should be /home/username/gomalib not ~/gomalib, ./gomalib, ../gomalib, etc."
    echo "current GOMA_LIB=$GOMA_LIB"
    continue_check
fi

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

HDF5_VERSION="1.12.2"
HDF5_MD5="4f4c87981ee5e3db12975f0904b634e2"

NETCDF_VERSION="c-4.8.1"
NETCDF_MD5="8da8665f1f185b85bc8127bc0bad0ee6"

TRILINOS_VERSION="13.2.0"
TRILINOS_VERSION_DASH="13-2-0"
TRILINOS_MD5="099680cd3660dba5ec447ddc50a8406c"

MUMPS_VERSION="5.5.0"
MUMPS_MD5="08eb0887bfd1e570ce7faecbd8bf97c0"

OPENMPI_VERSION="4.1.3"
OPENMPI_MD5="292fcbce24b1e6c18218e0fdcdc080e4"
OPENMPI_ARCHIVE_URL="https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-$OPENMPI_VERSION.tar.bz2"
OPENMPI_EXTRA_CONFIGURE_FLAGS=""

CMAKE_VERSION="3.21.1"
CMAKE_MD5="1d8d33628f1c56b0c3cda67abddbea91"

SUITESPARSE_VERSION="5.12.0"
SUITESPARSE_MD5="08292a05b16acf37767090875d210ced"

MATIO_VERSION="1.5.21"
MATIO_MD5="afeb5d21b234699fd5b9dc4564afe1ca"

SCALAPACK_VERSION="2.2.0"
SCALAPACK_MD5="2397d36790d1445383bc3cdb1e18ca5f"

LAPACK_VERSION="3.8.0"
LAPACK_MD5="96591affdbf58c450d45c1daa540dbd2"

PETSC_VERSION="3.17.1"
PETSC_MD5="7a39f099b99f2b03edb9e02876cebb6d"

OMEGA_H_VERSION="9.34.8"
OMEGA_H_MD5="daa3efaf5ea3aed32d2d60760d1a928e"

SUPERLU_DIST_VERSION="7.2.0"
SUPERLU_DIST_MD5="10d20b97012e5ae89a6cc69b768f61b7"

ARCHIVE_NAMES=("arpack96.tar.gz" \
"patch.tar.gz" \
"hdf5-${HDF5_VERSION}.tar.bz2" \
"netcdf-${NETCDF_VERSION}.tar.gz" \
"parmetis-4.0.3.tar.gz" \
"sparse.tar.gz" \
"superlu_dist-$SUPERLU_DIST_VERSION.tar.gz" \
"y12m.tar.gz" \
"Trilinos-trilinos-release-$TRILINOS_VERSION_DASH.tar.gz" \
"MUMPS_$MUMPS_VERSION.tar.gz" \
"SuiteSparse-$SUITESPARSE_VERSION.tar.gz" \
"matio-$MATIO_VERSION.tar.gz"  \
"petsc-$PETSC_VERSION.tar.gz" \
"omega_h-v$OMEGA_H_VERSION.tar.gz")

#y12m archive is skipped because it stores the number of downloads in the header

#meaning each y12m tar has a unique MD5SUM.
ARCHIVE_MD5SUMS=("fffaa970198b285676f4156cebc8626e" \
"14830d758f195f272b8594a493501fa2" \
"${HDF5_MD5}" \
"${NETCDF_MD5}" \
"f69c479586bf6bb7aff6a9bc0c739628" \
"1566d914d1035ac17b73fe9bc0eed02a" \
"$SUPERLU_DIST_MD5" \
"SKIP" \
$TRILINOS_MD5 \
$MUMPS_MD5 \
"$SUITESPARSE_MD5" \
"$MATIO_MD5" \
"$PETSC_MD5" \
"$OMEGA_H_MD5")

ARCHIVE_URLS=("http://www.caam.rice.edu/software/ARPACK/SRC/arpack96.tar.gz" \
"http://www.caam.rice.edu/software/ARPACK/SRC/patch.tar.gz" \
"https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-${HDF5_VERSION}/src/hdf5-${HDF5_VERSION}.tar.bz2" \
"https://downloads.unidata.ucar.edu/netcdf-c/4.8.1/netcdf-${NETCDF_VERSION}.tar.gz" \
"http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz" \
"http://downloads.sourceforge.net/project/sparse/sparse/sparse1.4b/sparse1.4b.tar.gz" \
"http://codeload.github.com/xiaoyeli/superlu_dist/tar.gz/v$SUPERLU_DIST_VERSION" \
"http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz\\&filename=y12m%2Fy12m.f" \
"https://github.com/trilinos/Trilinos/archive/trilinos-release-$TRILINOS_VERSION_DASH.tar.gz" \
"http://mumps.enseeiht.fr/MUMPS_$MUMPS_VERSION.tar.gz" \
"https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v$SUITESPARSE_VERSION.tar.gz" \
"https://github.com/tbeu/matio/releases/download/v$MATIO_VERSION/matio-$MATIO_VERSION.tar.gz" \
"https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-$PETSC_VERSION.tar.gz" \
"https://github.com/sandialabs/omega_h/archive/refs/tags/v$OMEGA_H_VERSION.tar.gz")

# You can't call the ARPACK patch ARPACK or it will think it is already extracted
# When in reality it isn't
ARCHIVE_DIR_NAMES=("ARPACK" \
"FAKE_DIR_FOR_ARPACK_PATCH" \
"hdf5-${HDF5_VERSION}" \
"netcdf-${NETCDF_VERSION}" \
"parmetis-4.0.3" \
"sparse" \
"superlu_dist-$SUPERLU_DIST_VERSION" \
"y12m" \
"Trilinos-trilinos-release-$TRILINOS_VERSION_DASH" \
"MUMPS_$MUMPS_VERSION" \
"SuiteSparse-$SUITESPARSE_VERSION" \
"matio-$MATIO_VERSION" \
"petsc-$PETSC_VERSION" \
"omega_h-$OMEGA_H_VERSION")

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
"https://sourceforge.net/projects/matio/" \
"https://petsc.org" \
"https://github.com/sandialabs/omega_h/")

ARCHIVE_REAL_NAMES=("ARPACK96" \
"HDF5" \
"NetCDF" \
"Metis" \
"Sparse" \
"SuperLU_DIST" \
"y12m" \
"Trilinos" \
"MUMPS" \
"SuiteSparse-$SUITESPARSE_VERSION" \
"MATIO" \
"PETSc" \
"Omega_h")

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
    if [[ "$cmake_vers" = $(echo -e "$cmake_vers\n3.17.1\n" | sort -V |tail -n1) ]]; then
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
    ARCHIVE_NAMES+=("cmake-$CMAKE_VERSION-Linux-x86_64.tar.gz")
    ARCHIVE_MD5SUMS+=("$CMAKE_MD5")
    ARCHIVE_URLS+=("https://github.com/Kitware/CMake/releases/download/v$CMAKE_VERSION/cmake-$CMAKE_VERSION-Linux-x86_64.tar.gz")
    ARCHIVE_DIR_NAMES+=("cmake-$CMAKE_VERSION-linux-x86_64")
    ARCHIVE_HOMEPAGES+=("https://cmake.org/")
    ARCHIVE_REAL_NAMES+=("CMake")
    log_echo "Cmake not found, will download."
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
    export GCC_EXTRA_FFLAGS=""
    if [[ "$CC_NAME" == "intel" ]]; then
        #Special flag only needed by intel compiler because reasons? No really why does GCC not need -fopenmp?
        export COMPILER_FLAG_MPI="-qopenmp"
        #TODO: Test if this fails with openmpi (probably)
        if [ "$MPI_NAME" == "open" ]; then
            export MPI_C_COMPILER="mpicc"
            export MPI_CXX_COMPILER="mpic++"
            export MPI_F90_COMPILER="mpif90"
            export MPI_F77_COMPILER="mpif77"
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
        GCC_VERSION=$(gcc --version | grep ^gcc | sed 's/^.* //g')
        if [[ "$GCC_VERSION" = $(echo -e "$GCC_VERSION\n10.0.0\n" | sort -V |tail -n1) ]]; then
            export GCC_EXTRA_FFLAGS="-fallow-argument-mismatch"
        fi
        export MPI_C_COMPILER="mpicc"
        export MPI_CXX_COMPILER="mpicxx"
        # mpif90 and mpif77 are depricated in favor of just mpifort.
        if [ "$MPI_NAME" == "intel" ]; then
            # Intel calls their mpi fortran compiler mpifc (for some reason
            export MPI_F90_COMPILER="mpifc"
            export MPI_F77_COMPILER="mpifc"
        else
            export MPI_F90_COMPILER="mpif90"
            export MPI_F77_COMPILER="mpif77"
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
        ARCHIVE_NAMES+=("lapack-$LAPACK_VERSION.tar.gz" \
        "scalapack-$SCALAPACK_VERSION.tgz")
        ARCHIVE_MD5SUMS+=("$LAPACK_MD5" \
        "$SCALAPACK_MD5" )
        ARCHIVE_URLS+=("https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v$LAPACK_VERSION.tar.gz" \
        "http://www.netlib.org/scalapack/scalapack-$SCALAPACK_VERSION.tgz" )
        ARCHIVE_DIR_NAMES+=("lapack-$LAPACK_VERSION" \
        "scalapack-$SCALAPACK_VERSION" )
        ARCHIVE_HOMEPAGES+=("http://www.netlib.org/lapack/" \
        "http://www.netlib.org/scalapack/")
        ARCHIVE_REAL_NAMES+=("LAPACK" \
        "ScaLAPACK")

        LAPACK_DIR="${GOMA_LIB}/lapack-$LAPACK_VERSION"
        LAPACK_LIBRARY_DIR="$LAPACK_DIR/lib"
        LAPACK_LIBRARY_NAME_ARG="-L${LAPACK_LIBRARY_DIR} -llapack"
        LAPACK_LIBRARY_NAME="lapack"
        BLAS_LIBRARY_DIR="${LAPACK_LIBRARY_DIR}"
        BLAS_LIBRARY_NAME_ARG="-L${BLAS_LIBRARY_DIR} -lblas"
        BLAS_LIBRARY_NAME="blas"
        BLACS_LIBRARY_NAME="scalapack"
        SCALAPACK_LIBRARY_DIR="${GOMA_LIB}/scalapack-$SCALAPACK_VERSION/lib"
        SCALAPACK_LIBRARY_NAME="scalapack"
        SCALAPACK_LIBRARY_NAME_ARG="-L${SCALAPACK_LIBRARY_DIR} -lscalapack"
        SCALAPACK_INCLUDE_DIR="${GOMA_LIB}/scalapack-$SCALAPACK_VERSION/include"
        # Alternative to the intel "blas_flags" variable
        NON_INTEL_BLAS_LIBRARY="${BLAS_LIBRARY_DIR}/libblas.a"
        NON_INTEL_LAPACK_LIBRARY="${LAPACK_LIBRARY_DIR}/liblapack.a"
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
        wget --no-check-certificate "${ARCHIVE_URLS[count]}" -O $i
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
    cd $GOMA_LIB/cmake-$CMAKE_VERSION-linux-x86_64
    if [ -f bin/cmake ]
    then
        log_echo "cmake is already built"
    else
        log_echo "Downloaded cmake does not include cmake in bin"
        exit 1
    fi
    export PATH=$GOMA_LIB/cmake-$CMAKE_VERSION-linux-x86_64/bin:$PATH
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
    CC="$MPI_C_COMPILER" CPP="$MPI_C_COMPILER -E" AR=${ARCHIVER} ./configure --enable-shared=off --prefix=$GOMA_LIB/hdf5-${HDF5_VERSION} --enable-parallel --with-default-api-version=v18 2>&1 | tee -a $COMPILE_LOG
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
> FFLAGS	= -O $BLAS_FLAGS $GCC_EXTRA_FFLAGS
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


if [[ "$MATH_LIBRARIES" == "netlib blas" ]]; then
    #make lapack/blas
    if [ -e $LAPACK_LIBRARY_DIR/liblapack.a ]
    then
        log_echo "LAPACK already built"
    else
        tempdir=$(mktemp -u -p $GOMA_LIB)
        mv $LAPACK_DIR $tempdir
	mkdir $LAPACK_DIR
	mv $tempdir $LAPACK_DIR/src
	mkdir $LAPACK_DIR/src/build
        cd $LAPACK_DIR/src/build
        CC=$MPI_C_COMPILER CXX=$MPI_CXX_COMPILER FC=$MPI_F90_COMPILER cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DCMAKE_INSTALL_PREFIX=$LAPACK_DIR $LAPACK_DIR/src 2>&1 | tee -a $COMPILE_LOG
	make -j$MAKE_JOBS 2>&1 | tee -a $COMPILE_LOG
	make install 2>&1 | tee -a $COMPILE_LOG
	cd $GOMA_LIB
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

#make SuperLU
cd $GOMA_LIB/superlu_dist-$SUPERLU_DIST_VERSION
if [ -e lib64/libsuperlu_dist.a ] || [ -e lib/libsuperlu_dist.a ]
then
    log_echo "SuperLU_DIST already built"
else
    mkdir -p build
    CC=$MPI_C_COMPILER CXX=$MPI_CXX_COMPILER FC=$MPI_F90_COMPILER cmake -B build -DTPL_BLAS_LIBRARIES="${NON_INTEL_BLAS_LINK}" \
      -DTPL_PARMETIS_LIBRARIES="$GOMA_LIB/parmetis-4.0.3/lib/libparmetis.a;$GOMA_LIB/parmetis-4.0.3/lib/libmetis.a" \
      -DCMAKE_INSTALL_PREFIX=$GOMA_LIB/superlu_dist-$SUPERLU_DIST_VERSION \
      -Denable_openmp:BOOL=FALSE \
      -DTPL_PARMETIS_INCLUDE_DIRS="$GOMA_LIB/parmetis-4.0.3/include" 2>&1 | tee -a $COMPILE_LOG
    make -C build 2>&1 | tee -a $COMPILE_LOG
    make -C build install 2>&1 | tee -a $COMPILE_LOG
    if [ -e $GOMA_LIB/superlu_dist-$SUPERLU_DIST_VERSION/lib/libsuperlu_dist.a ] || [ -e $GOMA_LIB/superlu_dist-$SUPERLU_DIST_VERSION/lib64/libsuperlu_dist.a ]
    then
        log_echo "Built SuperLU_DIST $SUPERLU_DIST_VERSION"
    else
        log_echo "Failed to build SuperLU_DIST $SUPERLU_DIST_VERSION"
        exit 1
    fi
fi
SUPERLU_LIBDIR="lib64"
if [ -e $GOMA_LIB/superlu_dist-$SUPERLU_DIST_VERSION/lib/libsuperlu_dist.a ]; then
  SUPERLU_LIBDIR="lib"
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
cd $GOMA_LIB/SuiteSparse-$SUITESPARSE_VERSION
if [ -e $GOMA_LIB/SuiteSparse-$SUITESPARSE_VERSION/lib/libumfpack.a ]
then
    log_echo "SuiteSparse is already built"
else
    log_echo "compiler flag MPI is ${COMPILER_FLAG_MPI}"

    log_echo ${MPI_C_COMPILER}
    if [ -z "${BLAS_FLAGS}" ]; then
	targets=(SuiteSparse_config AMD BTF CAMD CCOLAMD COLAMD CHOLMOD CXSparse LDL KLU UMFPACK RBio SPQR)
	make CUDA="no" AUTOCC="no" CC="${MPI_C_COMPILER}" \
             CXX="${MPI_CXX_COMPILER}" \
             BLAS="$NON_INTEL_BLAS_LINK" \
             LAPACK="$SUITESPARSE_NON_INTEL_LAPACK_LINK" \
	     MY_METIS_LIB="-L$GOMA_LIB/parmetis-4.0.3/lib -lmetis" \
	     MY_METIS_INC="$GOMA_LIB/parmetis-4.0.3/include" \
	     AR="${ARCHIVER}" -C SuiteSparse_config library config 2>&1 | tee -a $COMPILE_LOG
	for target in "${targets[@]}"; do
	    make CUDA="no" AUTOCC="no" CC="${MPI_C_COMPILER}" \
                 CXX="${MPI_CXX_COMPILER}" \
                 BLAS="$NON_INTEL_BLAS_LINK" \
                 LAPACK="$SUITESPARSE_NON_INTEL_LAPACK_LINK" \
	         MY_METIS_LIB="-L$GOMA_LIB/parmetis-4.0.3/lib -lmetis" \
	         MY_METIS_INC="$GOMA_LIB/parmetis-4.0.3/include" \
	         AR="${ARCHIVER}" -C $target library 2>&1 | tee -a $COMPILE_LOG
	    make CUDA="no" AUTOCC="no" CC="${MPI_C_COMPILER}" \
                 CXX="${MPI_CXX_COMPILER}" \
                 BLAS="$NON_INTEL_BLAS_LINK" \
                 LAPACK="$SUITESPARSE_NON_INTEL_LAPACK_LINK" \
	         MY_METIS_LIB="-L$GOMA_LIB/parmetis-4.0.3/lib -lmetis" \
	         MY_METIS_INC="$GOMA_LIB/parmetis-4.0.3/include" \
	         AR="${ARCHIVER}" -C $target install 2>&1 | tee -a $COMPILE_LOG
        done
    else
        make config CUDA="no" AUTOCC="no" CC="${MPI_C_COMPILER} ${COMPILER_FLAG_MPI}" \
             CXX="${MPI_CXX_COMPILER} ${COMPILER_FLAG_MPI}" \
             JOBS="$MAKE_JOBS" \
             AR="${ARCHIVER}" BLAS_FLAGS="${BLAS_FLAGS}" 2>&1 | tee -a $COMPILE_LOG
        make CUDA="no" AUTOCC="no" CC="${MPI_C_COMPILER} ${COMPILER_FLAG_MPI}" \
             CXX="${MPI_CXX_COMPILER} ${COMPILER_FLAG_MPI}" \
             JOBS="$MAKE_JOBS" \
             AR="${ARCHIVER}" BLAS_FLAGS="${BLAS_FLAGS}" 2>&1 | tee -a $COMPILE_LOG
    fi
    find . -name '*.a' -exec cp "{}" $GOMA_LIB/SuiteSparse-$SUITESPARSE_VERSION/lib \;
    if [ -e $GOMA_LIB/SuiteSparse-$SUITESPARSE_VERSION/lib/libumfpack.a ]
    then
        log_echo "Built SuiteSparse $SUITESPARSE_VERSION"
    else
        log_echo "Failed to build SuiteSparse $SUITESPARSE_VERSION"
        exit 1
    fi
fi
# Otherwise Goma dynamically links to UMFPACK when intel is sourced and disabled
export LD_LIBRARY_PATH="${GOMA_LIB}/SuiteSparse-$SUITESPARSE_VERSION/UMFPACK/Lib:$LD_LIBRARY_PATH"

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
        if [ ! -d " $GOMA_LIB/scalapack-$SCALAPACK_VERSION/src" ]; then
            mv $GOMA_LIB/scalapack-$SCALAPACK_VERSION src-scalapack
            mkdir $GOMA_LIB/scalapack-$SCALAPACK_VERSION
            mv src-scalapack $GOMA_LIB/scalapack-$SCALAPACK_VERSION/src
        fi
        cd $GOMA_LIB/scalapack-$SCALAPACK_VERSION/src
        # Use PETSc's parallel make changes
cat > Makefile.objs <<EOF
objslamov = \
SRC/clamov.o \
SRC/dlamov.o \
SRC/slamov.o \
SRC/zlamov.o \

objsblacs = \
BLACS/SRC/igesd2d_.o BLACS/SRC/sgesd2d_.o BLACS/SRC/dgesd2d_.o BLACS/SRC/cgesd2d_.o BLACS/SRC/zgesd2d_.o \
BLACS/SRC/itrsd2d_.o BLACS/SRC/strsd2d_.o BLACS/SRC/dtrsd2d_.o BLACS/SRC/ctrsd2d_.o BLACS/SRC/ztrsd2d_.o \
BLACS/SRC/igerv2d_.o BLACS/SRC/sgerv2d_.o BLACS/SRC/dgerv2d_.o BLACS/SRC/cgerv2d_.o BLACS/SRC/zgerv2d_.o \
BLACS/SRC/itrrv2d_.o BLACS/SRC/strrv2d_.o BLACS/SRC/dtrrv2d_.o BLACS/SRC/ctrrv2d_.o BLACS/SRC/ztrrv2d_.o \
BLACS/SRC/igebs2d_.o BLACS/SRC/sgebs2d_.o BLACS/SRC/dgebs2d_.o BLACS/SRC/cgebs2d_.o BLACS/SRC/zgebs2d_.o \
BLACS/SRC/igebr2d_.o BLACS/SRC/sgebr2d_.o BLACS/SRC/dgebr2d_.o BLACS/SRC/cgebr2d_.o BLACS/SRC/zgebr2d_.o \
BLACS/SRC/itrbs2d_.o BLACS/SRC/strbs2d_.o BLACS/SRC/dtrbs2d_.o BLACS/SRC/ctrbs2d_.o BLACS/SRC/ztrbs2d_.o \
BLACS/SRC/itrbr2d_.o BLACS/SRC/strbr2d_.o BLACS/SRC/dtrbr2d_.o BLACS/SRC/ctrbr2d_.o BLACS/SRC/ztrbr2d_.o \
BLACS/SRC/igsum2d_.o BLACS/SRC/sgsum2d_.o BLACS/SRC/dgsum2d_.o BLACS/SRC/cgsum2d_.o BLACS/SRC/zgsum2d_.o \
BLACS/SRC/igamx2d_.o BLACS/SRC/sgamx2d_.o BLACS/SRC/dgamx2d_.o BLACS/SRC/cgamx2d_.o BLACS/SRC/zgamx2d_.o \
BLACS/SRC/igamn2d_.o BLACS/SRC/sgamn2d_.o BLACS/SRC/dgamn2d_.o BLACS/SRC/cgamn2d_.o BLACS/SRC/zgamn2d_.o \
BLACS/SRC/blacs_setup_.o BLACS/SRC/blacs_set_.o BLACS/SRC/blacs_get_.o \
BLACS/SRC/blacs_abort_.o BLACS/SRC/blacs_exit_.o BLACS/SRC/blacs_pnum_.o BLACS/SRC/blacs_pcoord_.o \
BLACS/SRC/ksendid_.o BLACS/SRC/krecvid_.o BLACS/SRC/kbsid_.o BLACS/SRC/kbrid_.o \
BLACS/SRC/dcputime00_.o BLACS/SRC/dwalltime00_.o BLACS/SRC/blacs_pinfo_.o \
BLACS/SRC/blacs_init_.o BLACS/SRC/blacs_map_.o BLACS/SRC/blacs_free_.o BLACS/SRC/blacs_grid_.o BLACS/SRC/blacs_info_.o \
BLACS/SRC/blacs_barr_.o BLACS/SRC/sys2blacs_.o BLACS/SRC/blacs2sys_.o BLACS/SRC/free_handle_.o

objs = \
BLACS/SRC/BI_Arecv.o \
BLACS/SRC/BI_ArgCheck.o \
BLACS/SRC/BI_Asend.o \
BLACS/SRC/BI_BeComb.o \
BLACS/SRC/BI_BlacsAbort.o \
BLACS/SRC/BI_BlacsErr.o \
BLACS/SRC/BI_BlacsWarn.o \
BLACS/SRC/BI_BuffIsFree.o \
BLACS/SRC/BI_ContxtNum.o \
BLACS/SRC/BI_EmergencyBuff.o \
BLACS/SRC/BI_GetBuff.o \
BLACS/SRC/BI_GetMpiGeType.o \
BLACS/SRC/BI_GetMpiTrType.o \
BLACS/SRC/BI_GlobalVars.o \
BLACS/SRC/BI_HypBR.o \
BLACS/SRC/BI_HypBS.o \
BLACS/SRC/BI_IdringBR.o \
BLACS/SRC/BI_IdringBS.o \
BLACS/SRC/BI_MpathBR.o \
BLACS/SRC/BI_MpathBS.o \
BLACS/SRC/BI_MringComb.o \
BLACS/SRC/BI_Pack.o \
BLACS/SRC/BI_Rsend.o \
BLACS/SRC/BI_Srecv.o \
BLACS/SRC/BI_SringBR.o \
BLACS/SRC/BI_SringBS.o \
BLACS/SRC/BI_Ssend.o \
BLACS/SRC/BI_TransDist.o \
BLACS/SRC/BI_TransUserComm.o \
BLACS/SRC/BI_TreeBR.o \
BLACS/SRC/BI_TreeBS.o \
BLACS/SRC/BI_TreeComb.o \
BLACS/SRC/BI_Unpack.o \
BLACS/SRC/BI_UpdateBuffs.o \
BLACS/SRC/BI_cMPI_amn.o \
BLACS/SRC/BI_cMPI_amn2.o \
BLACS/SRC/BI_cMPI_amx.o \
BLACS/SRC/BI_cMPI_amx2.o \
BLACS/SRC/BI_cMPI_sum.o \
BLACS/SRC/BI_cvvamn.o \
BLACS/SRC/BI_cvvamn2.o \
BLACS/SRC/BI_cvvamx.o \
BLACS/SRC/BI_cvvamx2.o \
BLACS/SRC/BI_cvvsum.o \
BLACS/SRC/BI_dMPI_amn.o \
BLACS/SRC/BI_dMPI_amn2.o \
BLACS/SRC/BI_dMPI_amx.o \
BLACS/SRC/BI_dMPI_amx2.o \
BLACS/SRC/BI_dmvcopy.o \
BLACS/SRC/BI_dvmcopy.o \
BLACS/SRC/BI_dvvamn.o \
BLACS/SRC/BI_dvvamn2.o \
BLACS/SRC/BI_dvvamx.o \
BLACS/SRC/BI_dvvamx2.o \
BLACS/SRC/BI_dvvsum.o \
BLACS/SRC/BI_iMPI_amn.o \
BLACS/SRC/BI_iMPI_amn2.o \
BLACS/SRC/BI_iMPI_amx.o \
BLACS/SRC/BI_iMPI_amx2.o \
BLACS/SRC/BI_imvcopy.o \
BLACS/SRC/BI_ivmcopy.o \
BLACS/SRC/BI_ivvamn.o \
BLACS/SRC/BI_ivvamn2.o \
BLACS/SRC/BI_ivvamx.o \
BLACS/SRC/BI_ivvamx2.o \
BLACS/SRC/BI_ivvsum.o \
BLACS/SRC/BI_sMPI_amn.o \
BLACS/SRC/BI_sMPI_amn2.o \
BLACS/SRC/BI_sMPI_amx.o \
BLACS/SRC/BI_sMPI_amx2.o \
BLACS/SRC/BI_smvcopy.o \
BLACS/SRC/BI_svmcopy.o \
BLACS/SRC/BI_svvamn.o \
BLACS/SRC/BI_svvamn2.o \
BLACS/SRC/BI_svvamx.o \
BLACS/SRC/BI_svvamx2.o \
BLACS/SRC/BI_svvsum.o \
BLACS/SRC/BI_zMPI_amn.o \
BLACS/SRC/BI_zMPI_amn2.o \
BLACS/SRC/BI_zMPI_amx.o \
BLACS/SRC/BI_zMPI_amx2.o \
BLACS/SRC/BI_zMPI_sum.o \
BLACS/SRC/BI_zvvamn.o \
BLACS/SRC/BI_zvvamn2.o \
BLACS/SRC/BI_zvvamx.o \
BLACS/SRC/BI_zvvamx2.o \
BLACS/SRC/BI_zvvsum.o \
PBLAS/SRC/pcagemv_.o \
PBLAS/SRC/pcahemv_.o \
PBLAS/SRC/pcamax_.o \
PBLAS/SRC/pcatrmv_.o \
PBLAS/SRC/pcaxpy_.o \
PBLAS/SRC/pccopy_.o \
PBLAS/SRC/pcdotc_.o \
PBLAS/SRC/pcdotu_.o \
PBLAS/SRC/pcgeadd_.o \
PBLAS/SRC/pcgemm_.o \
PBLAS/SRC/pcgemv_.o \
PBLAS/SRC/pcgerc_.o \
PBLAS/SRC/pcgeru_.o \
PBLAS/SRC/pchemm_.o \
PBLAS/SRC/pchemv_.o \
PBLAS/SRC/pcher2_.o \
PBLAS/SRC/pcher2k_.o \
PBLAS/SRC/pcher_.o \
PBLAS/SRC/pcherk_.o \
PBLAS/SRC/pcscal_.o \
PBLAS/SRC/pcsscal_.o \
PBLAS/SRC/pcswap_.o \
PBLAS/SRC/pcsymm_.o \
PBLAS/SRC/pcsyr2k_.o \
PBLAS/SRC/pcsyrk_.o \
PBLAS/SRC/pctradd_.o \
PBLAS/SRC/pctranc_.o \
PBLAS/SRC/pctranu_.o \
PBLAS/SRC/pctrmm_.o \
PBLAS/SRC/pctrmv_.o \
PBLAS/SRC/pctrsm_.o \
PBLAS/SRC/pctrsv_.o \
PBLAS/SRC/pdagemv_.o \
PBLAS/SRC/pdamax_.o \
PBLAS/SRC/pdasum_.o \
PBLAS/SRC/pdasymv_.o \
PBLAS/SRC/pdatrmv_.o \
PBLAS/SRC/pdaxpy_.o \
PBLAS/SRC/pdcopy_.o \
PBLAS/SRC/pddot_.o \
PBLAS/SRC/pdgeadd_.o \
PBLAS/SRC/pdgemm_.o \
PBLAS/SRC/pdgemv_.o \
PBLAS/SRC/pdger_.o \
PBLAS/SRC/pdnrm2_.o \
PBLAS/SRC/pdscal_.o \
PBLAS/SRC/pdswap_.o \
PBLAS/SRC/pdsymm_.o \
PBLAS/SRC/pdsymv_.o \
PBLAS/SRC/pdsyr2_.o \
PBLAS/SRC/pdsyr2k_.o \
PBLAS/SRC/pdsyr_.o \
PBLAS/SRC/pdsyrk_.o \
PBLAS/SRC/pdtradd_.o \
PBLAS/SRC/pdtran_.o \
PBLAS/SRC/pdtrmm_.o \
PBLAS/SRC/pdtrmv_.o \
PBLAS/SRC/pdtrsm_.o \
PBLAS/SRC/pdtrsv_.o \
PBLAS/SRC/pdzasum_.o \
PBLAS/SRC/pdznrm2_.o \
PBLAS/SRC/picopy_.o \
PBLAS/SRC/pilaenv.o \
PBLAS/SRC/psagemv_.o \
PBLAS/SRC/psamax_.o \
PBLAS/SRC/psasum_.o \
PBLAS/SRC/psasymv_.o \
PBLAS/SRC/psatrmv_.o \
PBLAS/SRC/psaxpy_.o \
PBLAS/SRC/pscasum_.o \
PBLAS/SRC/pscnrm2_.o \
PBLAS/SRC/pscopy_.o \
PBLAS/SRC/psdot_.o \
PBLAS/SRC/psgeadd_.o \
PBLAS/SRC/psgemm_.o \
PBLAS/SRC/psgemv_.o \
PBLAS/SRC/psger_.o \
PBLAS/SRC/psnrm2_.o \
PBLAS/SRC/psscal_.o \
PBLAS/SRC/psswap_.o \
PBLAS/SRC/pssymm_.o \
PBLAS/SRC/pssymv_.o \
PBLAS/SRC/pssyr2_.o \
PBLAS/SRC/pssyr2k_.o \
PBLAS/SRC/pssyr_.o \
PBLAS/SRC/pssyrk_.o \
PBLAS/SRC/pstradd_.o \
PBLAS/SRC/pstran_.o \
PBLAS/SRC/pstrmm_.o \
PBLAS/SRC/pstrmv_.o \
PBLAS/SRC/pstrsm_.o \
PBLAS/SRC/pstrsv_.o \
PBLAS/SRC/pzagemv_.o \
PBLAS/SRC/pzahemv_.o \
PBLAS/SRC/pzamax_.o \
PBLAS/SRC/pzatrmv_.o \
PBLAS/SRC/pzaxpy_.o \
PBLAS/SRC/pzcopy_.o \
PBLAS/SRC/pzdotc_.o \
PBLAS/SRC/pzdotu_.o \
PBLAS/SRC/pzdscal_.o \
PBLAS/SRC/pzgeadd_.o \
PBLAS/SRC/pzgemm_.o \
PBLAS/SRC/pzgemv_.o \
PBLAS/SRC/pzgerc_.o \
PBLAS/SRC/pzgeru_.o \
PBLAS/SRC/pzhemm_.o \
PBLAS/SRC/pzhemv_.o \
PBLAS/SRC/pzher2_.o \
PBLAS/SRC/pzher2k_.o \
PBLAS/SRC/pzher_.o \
PBLAS/SRC/pzherk_.o \
PBLAS/SRC/pzscal_.o \
PBLAS/SRC/pzswap_.o \
PBLAS/SRC/pzsymm_.o \
PBLAS/SRC/pzsyr2k_.o \
PBLAS/SRC/pzsyrk_.o \
PBLAS/SRC/pztradd_.o \
PBLAS/SRC/pztranc_.o \
PBLAS/SRC/pztranu_.o \
PBLAS/SRC/pztrmm_.o \
PBLAS/SRC/pztrmv_.o \
PBLAS/SRC/pztrsm_.o \
PBLAS/SRC/pztrsv_.o \
REDIST/SRC/pcgemr.o \
REDIST/SRC/pcgemr2.o \
REDIST/SRC/pctrmr.o \
REDIST/SRC/pctrmr2.o \
REDIST/SRC/pdgemr.o \
REDIST/SRC/pdgemr2.o \
REDIST/SRC/pdtrmr.o \
REDIST/SRC/pdtrmr2.o \
REDIST/SRC/pgemraux.o \
REDIST/SRC/pigemr.o \
REDIST/SRC/pigemr2.o \
REDIST/SRC/pitrmr.o \
REDIST/SRC/pitrmr2.o \
REDIST/SRC/psgemr.o \
REDIST/SRC/psgemr2.o \
REDIST/SRC/pstrmr.o \
REDIST/SRC/pstrmr2.o \
REDIST/SRC/pzgemr.o \
REDIST/SRC/pzgemr2.o \
REDIST/SRC/pztrmr.o \
REDIST/SRC/pztrmr2.o \
SRC/bdlaapp.o \
SRC/bdlaexc.o \
SRC/bdtrexc.o \
SRC/bslaapp.o \
SRC/bslaexc.o \
SRC/bstrexc.o \
SRC/cdbtf2.o \
SRC/cdbtrf.o \
SRC/cdttrf.o \
SRC/cdttrsv.o \
SRC/clahqr2.o \
SRC/clamsh.o \
SRC/clanv2.o \
SRC/claref.o \
SRC/cpttrsv.o \
SRC/csteqr2.o \
SRC/ctrmvt.o \
SRC/ddbtf2.o \
SRC/ddbtrf.o \
SRC/ddttrf.o \
SRC/ddttrsv.o \
SRC/dlamsh.o \
SRC/dlapst.o \
SRC/dlaqr6.o \
SRC/dlar1va.o \
SRC/dlaref.o \
SRC/dlarrb2.o \
SRC/dlarrd2.o \
SRC/dlarre2.o \
SRC/dlarre2a.o \
SRC/dlarrf2.o \
SRC/dlarrv2.o \
SRC/dlasorte.o \
SRC/dlasrt2.o \
SRC/dpttrsv.o \
SRC/dstegr2.o \
SRC/dstegr2a.o \
SRC/dstegr2b.o \
SRC/dstein2.o \
SRC/dsteqr2.o \
SRC/dtrmvt.o \
SRC/getpbbuf.o \
SRC/pbchkvect.o \
SRC/pcdbsv.o \
SRC/pcdbtrf.o \
SRC/pcdbtrs.o \
SRC/pcdbtrsv.o \
SRC/pcdtsv.o \
SRC/pcdttrf.o \
SRC/pcdttrs.o \
SRC/pcdttrsv.o \
SRC/pcgbsv.o \
SRC/pcgbtrf.o \
SRC/pcgbtrs.o \
SRC/pcgebd2.o \
SRC/pcgebrd.o \
SRC/pcgecon.o \
SRC/pcgeequ.o \
SRC/pcgehd2.o \
SRC/pcgehrd.o \
SRC/pcgelq2.o \
SRC/pcgelqf.o \
SRC/pcgels.o \
SRC/pcgeql2.o \
SRC/pcgeqlf.o \
SRC/pcgeqpf.o \
SRC/pcgeqr2.o \
SRC/pcgeqrf.o \
SRC/pcgerfs.o \
SRC/pcgerq2.o \
SRC/pcgerqf.o \
SRC/pcgesv.o \
SRC/pcgesvd.o \
SRC/pcgesvx.o \
SRC/pcgetf2.o \
SRC/pcgetrf.o \
SRC/pcgetri.o \
SRC/pcgetrs.o \
SRC/pcggqrf.o \
SRC/pcggrqf.o \
SRC/pcheev.o \
SRC/pcheevd.o \
SRC/pcheevr.o \
SRC/pcheevx.o \
SRC/pchegs2.o \
SRC/pchegst.o \
SRC/pchegvx.o \
SRC/pchengst.o \
SRC/pchentrd.o \
SRC/pchetd2.o \
SRC/pchetrd.o \
SRC/pchettrd.o \
SRC/pclabrd.o \
SRC/pclacgv.o \
SRC/pclacon.o \
SRC/pclaconsb.o \
SRC/pclacp2.o \
SRC/pclacp3.o \
SRC/pclacpy.o \
SRC/pclaevswp.o \
SRC/pclahqr.o \
SRC/pclahrd.o \
SRC/pclamr1d.o \
SRC/pclange.o \
SRC/pclanhe.o \
SRC/pclanhs.o \
SRC/pclansy.o \
SRC/pclantr.o \
SRC/pclapiv.o \
SRC/pclapv2.o \
SRC/pclaqge.o \
SRC/pclaqsy.o \
SRC/pclarf.o \
SRC/pclarfb.o \
SRC/pclarfc.o \
SRC/pclarfg.o \
SRC/pclarft.o \
SRC/pclarz.o \
SRC/pclarzb.o \
SRC/pclarzc.o \
SRC/pclarzt.o \
SRC/pclascl.o \
SRC/pclase2.o \
SRC/pclaset.o \
SRC/pclasmsub.o \
SRC/pclassq.o \
SRC/pclaswp.o \
SRC/pclatra.o \
SRC/pclatrd.o \
SRC/pclatrs.o \
SRC/pclatrz.o \
SRC/pclattrs.o \
SRC/pclauu2.o \
SRC/pclauum.o \
SRC/pclawil.o \
SRC/pcmax1.o \
SRC/pcpbsv.o \
SRC/pcpbtrf.o \
SRC/pcpbtrs.o \
SRC/pcpbtrsv.o \
SRC/pcpocon.o \
SRC/pcpoequ.o \
SRC/pcporfs.o \
SRC/pcposv.o \
SRC/pcposvx.o \
SRC/pcpotf2.o \
SRC/pcpotrf.o \
SRC/pcpotri.o \
SRC/pcpotrs.o \
SRC/pcptsv.o \
SRC/pcpttrf.o \
SRC/pcpttrs.o \
SRC/pcpttrsv.o \
SRC/pcrot.o \
SRC/pcsrscl.o \
SRC/pcstein.o \
SRC/pctrcon.o \
SRC/pctrevc.o \
SRC/pctrrfs.o \
SRC/pctrti2.o \
SRC/pctrtri.o \
SRC/pctrtrs.o \
SRC/pctzrzf.o \
SRC/pcung2l.o \
SRC/pcung2r.o \
SRC/pcungl2.o \
SRC/pcunglq.o \
SRC/pcungql.o \
SRC/pcungqr.o \
SRC/pcungr2.o \
SRC/pcungrq.o \
SRC/pcunm2l.o \
SRC/pcunm2r.o \
SRC/pcunmbr.o \
SRC/pcunmhr.o \
SRC/pcunml2.o \
SRC/pcunmlq.o \
SRC/pcunmql.o \
SRC/pcunmqr.o \
SRC/pcunmr2.o \
SRC/pcunmr3.o \
SRC/pcunmrq.o \
SRC/pcunmrz.o \
SRC/pcunmtr.o \
SRC/pddbsv.o \
SRC/pddbtrf.o \
SRC/pddbtrs.o \
SRC/pddbtrsv.o \
SRC/pddtsv.o \
SRC/pddttrf.o \
SRC/pddttrs.o \
SRC/pddttrsv.o \
SRC/pdgbsv.o \
SRC/pdgbtrf.o \
SRC/pdgbtrs.o \
SRC/pdgebal.o \
SRC/pdgebd2.o \
SRC/pdgebrd.o \
SRC/pdgecon.o \
SRC/pdgeequ.o \
SRC/pdgehd2.o \
SRC/pdgehrd.o \
SRC/pdgelq2.o \
SRC/pdgelqf.o \
SRC/pdgels.o \
SRC/pdgeql2.o \
SRC/pdgeqlf.o \
SRC/pdgeqpf.o \
SRC/pdgeqr2.o \
SRC/pdgeqrf.o \
SRC/pdgerfs.o \
SRC/pdgerq2.o \
SRC/pdgerqf.o \
SRC/pdgesv.o \
SRC/pdgesvd.o \
SRC/pdgesvx.o \
SRC/pdgetf2.o \
SRC/pdgetrf.o \
SRC/pdgetri.o \
SRC/pdgetrs.o \
SRC/pdggqrf.o \
SRC/pdggrqf.o \
SRC/pdhseqr.o \
SRC/pdlabad.o \
SRC/pdlabrd.o \
SRC/pdlacon.o \
SRC/pdlaconsb.o \
SRC/pdlacp2.o \
SRC/pdlacp3.o \
SRC/pdlacpy.o \
SRC/pdlaed0.o \
SRC/pdlaed1.o \
SRC/pdlaed2.o \
SRC/pdlaed3.o \
SRC/pdlaedz.o \
SRC/pdlaevswp.o \
SRC/pdlahqr.o \
SRC/pdlahrd.o \
SRC/pdlaiect.o \
SRC/pdlamch.o \
SRC/pdlamr1d.o \
SRC/pdlamve.o \
SRC/pdlange.o \
SRC/pdlanhs.o \
SRC/pdlansy.o \
SRC/pdlantr.o \
SRC/pdlapiv.o \
SRC/pdlapv2.o \
SRC/pdlaqge.o \
SRC/pdlaqr0.o \
SRC/pdlaqr1.o \
SRC/pdlaqr2.o \
SRC/pdlaqr3.o \
SRC/pdlaqr4.o \
SRC/pdlaqr5.o \
SRC/pdlaqsy.o \
SRC/pdlared1d.o \
SRC/pdlared2d.o \
SRC/pdlarf.o \
SRC/pdlarfb.o \
SRC/pdlarfg.o \
SRC/pdlarft.o \
SRC/pdlarz.o \
SRC/pdlarzb.o \
SRC/pdlarzt.o \
SRC/pdlascl.o \
SRC/pdlase2.o \
SRC/pdlaset.o \
SRC/pdlasmsub.o \
SRC/pdlasrt.o \
SRC/pdlassq.o \
SRC/pdlaswp.o \
SRC/pdlatra.o \
SRC/pdlatrd.o \
SRC/pdlatrs.o \
SRC/pdlatrz.o \
SRC/pdlauu2.o \
SRC/pdlauum.o \
SRC/pdlawil.o \
SRC/pdorg2l.o \
SRC/pdorg2r.o \
SRC/pdorgl2.o \
SRC/pdorglq.o \
SRC/pdorgql.o \
SRC/pdorgqr.o \
SRC/pdorgr2.o \
SRC/pdorgrq.o \
SRC/pdorm2l.o \
SRC/pdorm2r.o \
SRC/pdormbr.o \
SRC/pdormhr.o \
SRC/pdorml2.o \
SRC/pdormlq.o \
SRC/pdormql.o \
SRC/pdormqr.o \
SRC/pdormr2.o \
SRC/pdormr3.o \
SRC/pdormrq.o \
SRC/pdormrz.o \
SRC/pdormtr.o \
SRC/pdpbsv.o \
SRC/pdpbtrf.o \
SRC/pdpbtrs.o \
SRC/pdpbtrsv.o \
SRC/pdpocon.o \
SRC/pdpoequ.o \
SRC/pdporfs.o \
SRC/pdposv.o \
SRC/pdposvx.o \
SRC/pdpotf2.o \
SRC/pdpotrf.o \
SRC/pdpotri.o \
SRC/pdpotrs.o \
SRC/pdptsv.o \
SRC/pdpttrf.o \
SRC/pdpttrs.o \
SRC/pdpttrsv.o \
SRC/pdrot.o \
SRC/pdrscl.o \
SRC/pdstebz.o \
SRC/pdstedc.o \
SRC/pdstein.o \
SRC/pdsyev.o \
SRC/pdsyevd.o \
SRC/pdsyevr.o \
SRC/pdsyevx.o \
SRC/pdsygs2.o \
SRC/pdsygst.o \
SRC/pdsygvx.o \
SRC/pdsyngst.o \
SRC/pdsyntrd.o \
SRC/pdsytd2.o \
SRC/pdsytrd.o \
SRC/pdsyttrd.o \
SRC/pdtrcon.o \
SRC/pdtrord.o \
SRC/pdtrrfs.o \
SRC/pdtrsen.o \
SRC/pdtrti2.o \
SRC/pdtrtri.o \
SRC/pdtrtrs.o \
SRC/pdtzrzf.o \
SRC/pdzsum1.o \
SRC/pilaenvx.o \
SRC/pilaver.o \
SRC/piparmq.o \
SRC/pjlaenv.o \
SRC/pmpcol.o \
SRC/pmpim2.o \
SRC/pscsum1.o \
SRC/psdbsv.o \
SRC/psdbtrf.o \
SRC/psdbtrs.o \
SRC/psdbtrsv.o \
SRC/psdtsv.o \
SRC/psdttrf.o \
SRC/psdttrs.o \
SRC/psdttrsv.o \
SRC/psgbsv.o \
SRC/psgbtrf.o \
SRC/psgbtrs.o \
SRC/psgebal.o \
SRC/psgebd2.o \
SRC/psgebrd.o \
SRC/psgecon.o \
SRC/psgeequ.o \
SRC/psgehd2.o \
SRC/psgehrd.o \
SRC/psgelq2.o \
SRC/psgelqf.o \
SRC/psgels.o \
SRC/psgeql2.o \
SRC/psgeqlf.o \
SRC/psgeqpf.o \
SRC/psgeqr2.o \
SRC/psgeqrf.o \
SRC/psgerfs.o \
SRC/psgerq2.o \
SRC/psgerqf.o \
SRC/psgesv.o \
SRC/psgesvd.o \
SRC/psgesvx.o \
SRC/psgetf2.o \
SRC/psgetrf.o \
SRC/psgetri.o \
SRC/psgetrs.o \
SRC/psggqrf.o \
SRC/psggrqf.o \
SRC/pshseqr.o \
SRC/pslabad.o \
SRC/pslabrd.o \
SRC/pslacon.o \
SRC/pslaconsb.o \
SRC/pslacp2.o \
SRC/pslacp3.o \
SRC/pslacpy.o \
SRC/pslaed0.o \
SRC/pslaed1.o \
SRC/pslaed2.o \
SRC/pslaed3.o \
SRC/pslaedz.o \
SRC/pslaevswp.o \
SRC/pslahqr.o \
SRC/pslahrd.o \
SRC/pslaiect.o \
SRC/pslamch.o \
SRC/pslamr1d.o \
SRC/pslamve.o \
SRC/pslange.o \
SRC/pslanhs.o \
SRC/pslansy.o \
SRC/pslantr.o \
SRC/pslapiv.o \
SRC/pslapv2.o \
SRC/pslaqge.o \
SRC/pslaqr0.o \
SRC/pslaqr1.o \
SRC/pslaqr2.o \
SRC/pslaqr3.o \
SRC/pslaqr4.o \
SRC/pslaqr5.o \
SRC/pslaqsy.o \
SRC/pslared1d.o \
SRC/pslared2d.o \
SRC/pslarf.o \
SRC/pslarfb.o \
SRC/pslarfg.o \
SRC/pslarft.o \
SRC/pslarz.o \
SRC/pslarzb.o \
SRC/pslarzt.o \
SRC/pslascl.o \
SRC/pslase2.o \
SRC/pslaset.o \
SRC/pslasmsub.o \
SRC/pslasrt.o \
SRC/pslassq.o \
SRC/pslaswp.o \
SRC/pslatra.o \
SRC/pslatrd.o \
SRC/pslatrs.o \
SRC/pslatrz.o \
SRC/pslauu2.o \
SRC/pslauum.o \
SRC/pslawil.o \
SRC/psorg2l.o \
SRC/psorg2r.o \
SRC/psorgl2.o \
SRC/psorglq.o \
SRC/psorgql.o \
SRC/psorgqr.o \
SRC/psorgr2.o \
SRC/psorgrq.o \
SRC/psorm2l.o \
SRC/psorm2r.o \
SRC/psormbr.o \
SRC/psormhr.o \
SRC/psorml2.o \
SRC/psormlq.o \
SRC/psormql.o \
SRC/psormqr.o \
SRC/psormr2.o \
SRC/psormr3.o \
SRC/psormrq.o \
SRC/psormrz.o \
SRC/psormtr.o \
SRC/pspbsv.o \
SRC/pspbtrf.o \
SRC/pspbtrs.o \
SRC/pspbtrsv.o \
SRC/pspocon.o \
SRC/pspoequ.o \
SRC/psporfs.o \
SRC/psposv.o \
SRC/psposvx.o \
SRC/pspotf2.o \
SRC/pspotrf.o \
SRC/pspotri.o \
SRC/pspotrs.o \
SRC/psptsv.o \
SRC/pspttrf.o \
SRC/pspttrs.o \
SRC/pspttrsv.o \
SRC/psrot.o \
SRC/psrscl.o \
SRC/psstebz.o \
SRC/psstedc.o \
SRC/psstein.o \
SRC/pssyev.o \
SRC/pssyevd.o \
SRC/pssyevr.o \
SRC/pssyevx.o \
SRC/pssygs2.o \
SRC/pssygst.o \
SRC/pssygvx.o \
SRC/pssyngst.o \
SRC/pssyntrd.o \
SRC/pssytd2.o \
SRC/pssytrd.o \
SRC/pssyttrd.o \
SRC/pstrcon.o \
SRC/pstrord.o \
SRC/pstrrfs.o \
SRC/pstrsen.o \
SRC/pstrti2.o \
SRC/pstrtri.o \
SRC/pstrtrs.o \
SRC/pstzrzf.o \
SRC/pzaxpy.o \
SRC/pzdbsv.o \
SRC/pzdbtrf.o \
SRC/pzdbtrs.o \
SRC/pzdbtrsv.o \
SRC/pzdotc.o \
SRC/pzdotu.o \
SRC/pzdrscl.o \
SRC/pzdtsv.o \
SRC/pzdttrf.o \
SRC/pzdttrs.o \
SRC/pzdttrsv.o \
SRC/pzgbsv.o \
SRC/pzgbtrf.o \
SRC/pzgbtrs.o \
SRC/pzgebd2.o \
SRC/pzgebrd.o \
SRC/pzgecon.o \
SRC/pzgeequ.o \
SRC/pzgehd2.o \
SRC/pzgehrd.o \
SRC/pzgelq2.o \
SRC/pzgelqf.o \
SRC/pzgels.o \
SRC/pzgeql2.o \
SRC/pzgeqlf.o \
SRC/pzgeqpf.o \
SRC/pzgeqr2.o \
SRC/pzgeqrf.o \
SRC/pzgerfs.o \
SRC/pzgerq2.o \
SRC/pzgerqf.o \
SRC/pzgesv.o \
SRC/pzgesvd.o \
SRC/pzgesvx.o \
SRC/pzgetf2.o \
SRC/pzgetrf.o \
SRC/pzgetri.o \
SRC/pzgetrs.o \
SRC/pzggqrf.o \
SRC/pzggrqf.o \
SRC/pzheev.o \
SRC/pzheevd.o \
SRC/pzheevr.o \
SRC/pzheevx.o \
SRC/pzhegs2.o \
SRC/pzhegst.o \
SRC/pzhegvx.o \
SRC/pzhengst.o \
SRC/pzhentrd.o \
SRC/pzhetd2.o \
SRC/pzhetrd.o \
SRC/pzhettrd.o \
SRC/pzlabrd.o \
SRC/pzlacgv.o \
SRC/pzlacon.o \
SRC/pzlaconsb.o \
SRC/pzlacp2.o \
SRC/pzlacp3.o \
SRC/pzlacpy.o \
SRC/pzlaevswp.o \
SRC/pzlahqr.o \
SRC/pzlahrd.o \
SRC/pzlamr1d.o \
SRC/pzlange.o \
SRC/pzlanhe.o \
SRC/pzlanhs.o \
SRC/pzlansy.o \
SRC/pzlantr.o \
SRC/pzlapiv.o \
SRC/pzlapv2.o \
SRC/pzlaqge.o \
SRC/pzlaqsy.o \
SRC/pzlarf.o \
SRC/pzlarfb.o \
SRC/pzlarfc.o \
SRC/pzlarfg.o \
SRC/pzlarft.o \
SRC/pzlarz.o \
SRC/pzlarzb.o \
SRC/pzlarzc.o \
SRC/pzlarzt.o \
SRC/pzlascl.o \
SRC/pzlase2.o \
SRC/pzlaset.o \
SRC/pzlasmsub.o \
SRC/pzlassq.o \
SRC/pzlaswp.o \
SRC/pzlatra.o \
SRC/pzlatrd.o \
SRC/pzlatrs.o \
SRC/pzlatrz.o \
SRC/pzlattrs.o \
SRC/pzlauu2.o \
SRC/pzlauum.o \
SRC/pzlawil.o \
SRC/pzmax1.o \
SRC/pzpbsv.o \
SRC/pzpbtrf.o \
SRC/pzpbtrs.o \
SRC/pzpbtrsv.o \
SRC/pzpocon.o \
SRC/pzpoequ.o \
SRC/pzporfs.o \
SRC/pzposv.o \
SRC/pzposvx.o \
SRC/pzpotf2.o \
SRC/pzpotrf.o \
SRC/pzpotri.o \
SRC/pzpotrs.o \
SRC/pzptsv.o \
SRC/pzpttrf.o \
SRC/pzpttrs.o \
SRC/pzpttrsv.o \
SRC/pzrot.o \
SRC/pzstein.o \
SRC/pztrcon.o \
SRC/pztrevc.o \
SRC/pztrrfs.o \
SRC/pztrti2.o \
SRC/pztrtri.o \
SRC/pztrtrs.o \
SRC/pztzrzf.o \
SRC/pzung2l.o \
SRC/pzung2r.o \
SRC/pzungl2.o \
SRC/pzunglq.o \
SRC/pzungql.o \
SRC/pzungqr.o \
SRC/pzungr2.o \
SRC/pzungrq.o \
SRC/pzunm2l.o \
SRC/pzunm2r.o \
SRC/pzunmbr.o \
SRC/pzunmhr.o \
SRC/pzunml2.o \
SRC/pzunmlq.o \
SRC/pzunmql.o \
SRC/pzunmqr.o \
SRC/pzunmr2.o \
SRC/pzunmr3.o \
SRC/pzunmrq.o \
SRC/pzunmrz.o \
SRC/pzunmtr.o \
SRC/sdbtf2.o \
SRC/sdbtrf.o \
SRC/sdttrf.o \
SRC/sdttrsv.o \
SRC/slamsh.o \
SRC/slapst.o \
SRC/slaqr6.o \
SRC/slar1va.o \
SRC/slaref.o \
SRC/slarrb2.o \
SRC/slarrd2.o \
SRC/slarre2.o \
SRC/slarre2a.o \
SRC/slarrf2.o \
SRC/slarrv2.o \
SRC/slasorte.o \
SRC/slasrt2.o \
SRC/spttrsv.o \
SRC/sstegr2.o \
SRC/sstegr2a.o \
SRC/sstegr2b.o \
SRC/sstein2.o \
SRC/ssteqr2.o \
SRC/strmvt.o \
SRC/zdbtf2.o \
SRC/zdbtrf.o \
SRC/zdttrf.o \
SRC/zdttrsv.o \
SRC/zlahqr2.o \
SRC/zlamsh.o \
SRC/zlanv2.o \
SRC/zlaref.o \
SRC/zpttrsv.o \
SRC/zsteqr2.o \
SRC/ztrmvt.o \
TOOLS/LAPACK/clagge.o \
TOOLS/LAPACK/claghe.o \
TOOLS/LAPACK/clagsy.o \
TOOLS/LAPACK/clarnd.o \
TOOLS/LAPACK/clarnv.o \
TOOLS/LAPACK/clarot.o \
TOOLS/LAPACK/clatm1.o \
TOOLS/LAPACK/clatms.o \
TOOLS/LAPACK/dlagge.o \
TOOLS/LAPACK/dlagsy.o \
TOOLS/LAPACK/dlaran.o \
TOOLS/LAPACK/dlarnd.o \
TOOLS/LAPACK/dlarot.o \
TOOLS/LAPACK/dlatm1.o \
TOOLS/LAPACK/dlatms.o \
TOOLS/LAPACK/icopy.o \
TOOLS/LAPACK/slagge.o \
TOOLS/LAPACK/slagsy.o \
TOOLS/LAPACK/slaran.o \
TOOLS/LAPACK/slarnd.o \
TOOLS/LAPACK/slarot.o \
TOOLS/LAPACK/slatm1.o \
TOOLS/LAPACK/slatms.o \
TOOLS/LAPACK/zlagge.o \
TOOLS/LAPACK/zlaghe.o \
TOOLS/LAPACK/zlagsy.o \
TOOLS/LAPACK/zlarnd.o \
TOOLS/LAPACK/zlarnv.o \
TOOLS/LAPACK/zlarot.o \
TOOLS/LAPACK/zlatm1.o \
TOOLS/LAPACK/zlatms.o \
TOOLS/SL_gridreshape.o \
TOOLS/SL_init.o \
TOOLS/ccdotc.o \
TOOLS/ccdotu.o \
TOOLS/chk1mat.o \
TOOLS/clatcpy.o \
TOOLS/cmatadd.o \
TOOLS/dddot.o \
TOOLS/desc_convert.o \
TOOLS/descinit.o \
TOOLS/descset.o \
TOOLS/dlatcpy.o \
TOOLS/dmatadd.o \
TOOLS/dsasum.o \
TOOLS/dscasum.o \
TOOLS/dscnrm2.o \
TOOLS/dsnrm2.o \
TOOLS/iceil.o \
TOOLS/ilacpy.o \
TOOLS/ilcm.o \
TOOLS/indxg2l.o \
TOOLS/indxg2p.o \
TOOLS/indxl2g.o \
TOOLS/infog1l.o \
TOOLS/infog2l.o \
TOOLS/npreroc.o \
TOOLS/numroc.o \
TOOLS/pcchekpad.o \
TOOLS/pccol2row.o \
TOOLS/pcelget.o \
TOOLS/pcelset.o \
TOOLS/pcelset2.o \
TOOLS/pcfillpad.o \
TOOLS/pchkxmat.o \
TOOLS/pclaprnt.o \
TOOLS/pclaread.o \
TOOLS/pclawrite.o \
TOOLS/pcmatadd.o \
TOOLS/pcrow2col.o \
TOOLS/pctreecomb.o \
TOOLS/pdchekpad.o \
TOOLS/pdcol2row.o \
TOOLS/pdelget.o \
TOOLS/pdelset.o \
TOOLS/pdelset2.o \
TOOLS/pdfillpad.o \
TOOLS/pdlaprnt.o \
TOOLS/pdlaread.o \
TOOLS/pdlawrite.o \
TOOLS/pdmatadd.o \
TOOLS/pdrow2col.o \
TOOLS/pdtreecomb.o \
TOOLS/pichekpad.o \
TOOLS/picol2row.o \
TOOLS/pielget.o \
TOOLS/pielset.o \
TOOLS/pielset2.o \
TOOLS/pifillpad.o \
TOOLS/pilaprnt.o \
TOOLS/pirow2col.o \
TOOLS/pitreecomb.o \
TOOLS/pschekpad.o \
TOOLS/pscol2row.o \
TOOLS/pselget.o \
TOOLS/pselset.o \
TOOLS/pselset2.o \
TOOLS/psfillpad.o \
TOOLS/pslaprnt.o \
TOOLS/pslaread.o \
TOOLS/pslawrite.o \
TOOLS/psmatadd.o \
TOOLS/psrow2col.o \
TOOLS/pstreecomb.o \
TOOLS/pzchekpad.o \
TOOLS/pzcol2row.o \
TOOLS/pzelget.o \
TOOLS/pzelset.o \
TOOLS/pzelset2.o \
TOOLS/pzfillpad.o \
TOOLS/pzlaprnt.o \
TOOLS/pzlaread.o \
TOOLS/pzlawrite.o \
TOOLS/pzmatadd.o \
TOOLS/pzrow2col.o \
TOOLS/pztreecomb.o \
TOOLS/reshape.o \
TOOLS/slatcpy.o \
TOOLS/sltimer.o \
TOOLS/smatadd.o \
TOOLS/ssdot.o \
TOOLS/zlatcpy.o \
TOOLS/zmatadd.o \
TOOLS/zzdotc.o \
TOOLS/zzdotu.o \
PBLAS/SRC/PBBLAS/pbcmatadd.o \
PBLAS/SRC/PBBLAS/pbctran.o \
PBLAS/SRC/PBBLAS/pbctrget.o \
PBLAS/SRC/PBBLAS/pbctrnv.o \
PBLAS/SRC/PBBLAS/pbctrsrt.o \
PBLAS/SRC/PBBLAS/pbctrst1.o \
PBLAS/SRC/PBBLAS/pbcvecadd.o \
PBLAS/SRC/PBBLAS/pbdmatadd.o \
PBLAS/SRC/PBBLAS/pbdtran.o \
PBLAS/SRC/PBBLAS/pbdtrget.o \
PBLAS/SRC/PBBLAS/pbdtrnv.o \
PBLAS/SRC/PBBLAS/pbdtrsrt.o \
PBLAS/SRC/PBBLAS/pbdtrst1.o \
PBLAS/SRC/PBBLAS/pbdvecadd.o \
PBLAS/SRC/PBBLAS/pbsmatadd.o \
PBLAS/SRC/PBBLAS/pbstran.o \
PBLAS/SRC/PBBLAS/pbstrget.o \
PBLAS/SRC/PBBLAS/pbstrnv.o \
PBLAS/SRC/PBBLAS/pbstrsrt.o \
PBLAS/SRC/PBBLAS/pbstrst1.o \
PBLAS/SRC/PBBLAS/pbsvecadd.o \
PBLAS/SRC/PBBLAS/pbzmatadd.o \
PBLAS/SRC/PBBLAS/pbztran.o \
PBLAS/SRC/PBBLAS/pbztrget.o \
PBLAS/SRC/PBBLAS/pbztrnv.o \
PBLAS/SRC/PBBLAS/pbztrsrt.o \
PBLAS/SRC/PBBLAS/pbztrst1.o \
PBLAS/SRC/PBBLAS/pbzvecadd.o \
PBLAS/SRC/PTOOLS/PB_CGatherV.o \
PBLAS/SRC/PTOOLS/PB_CInOutV.o \
PBLAS/SRC/PTOOLS/PB_CInOutV2.o \
PBLAS/SRC/PTOOLS/PB_CInV.o \
PBLAS/SRC/PTOOLS/PB_CInV2.o \
PBLAS/SRC/PTOOLS/PB_COutV.o \
PBLAS/SRC/PTOOLS/PB_CScatterV.o \
PBLAS/SRC/PTOOLS/PB_CVMcontig.o \
PBLAS/SRC/PTOOLS/PB_CVMinit.o \
PBLAS/SRC/PTOOLS/PB_CVMloc.o \
PBLAS/SRC/PTOOLS/PB_CVMnpq.o \
PBLAS/SRC/PTOOLS/PB_CVMpack.o \
PBLAS/SRC/PTOOLS/PB_CVMswp.o \
PBLAS/SRC/PTOOLS/PB_CVMupdate.o \
PBLAS/SRC/PTOOLS/PB_Cabort.o \
PBLAS/SRC/PTOOLS/PB_Cainfog2l.o \
PBLAS/SRC/PTOOLS/PB_CargFtoC.o \
PBLAS/SRC/PTOOLS/PB_Cbinfo.o \
PBLAS/SRC/PTOOLS/PB_Cchkmat.o \
PBLAS/SRC/PTOOLS/PB_Cchkvec.o \
PBLAS/SRC/PTOOLS/PB_Cconjg.o \
PBLAS/SRC/PTOOLS/PB_Cctypeset.o \
PBLAS/SRC/PTOOLS/PB_Cdescribe.o \
PBLAS/SRC/PTOOLS/PB_Cdescset.o \
PBLAS/SRC/PTOOLS/PB_Cdtypeset.o \
PBLAS/SRC/PTOOLS/PB_Cfirstnb.o \
PBLAS/SRC/PTOOLS/PB_Cg2lrem.o \
PBLAS/SRC/PTOOLS/PB_Cgcd.o \
PBLAS/SRC/PTOOLS/PB_Cgetbuf.o \
PBLAS/SRC/PTOOLS/PB_Cindxg2p.o \
PBLAS/SRC/PTOOLS/PB_Cinfog2l.o \
PBLAS/SRC/PTOOLS/PB_Citypeset.o \
PBLAS/SRC/PTOOLS/PB_Clastnb.o \
PBLAS/SRC/PTOOLS/PB_Clcm.o \
PBLAS/SRC/PTOOLS/PB_Cmalloc.o \
PBLAS/SRC/PTOOLS/PB_Cnnxtroc.o \
PBLAS/SRC/PTOOLS/PB_Cnpreroc.o \
PBLAS/SRC/PTOOLS/PB_Cnumroc.o \
PBLAS/SRC/PTOOLS/PB_Cpaxpby.o \
PBLAS/SRC/PTOOLS/PB_CpaxpbyDN.o \
PBLAS/SRC/PTOOLS/PB_CpaxpbyND.o \
PBLAS/SRC/PTOOLS/PB_CpaxpbyNN.o \
PBLAS/SRC/PTOOLS/PB_Cpdot11.o \
PBLAS/SRC/PTOOLS/PB_CpdotND.o \
PBLAS/SRC/PTOOLS/PB_CpdotNN.o \
PBLAS/SRC/PTOOLS/PB_Cpgeadd.o \
PBLAS/SRC/PTOOLS/PB_CpgemmAB.o \
PBLAS/SRC/PTOOLS/PB_CpgemmAC.o \
PBLAS/SRC/PTOOLS/PB_CpgemmBC.o \
PBLAS/SRC/PTOOLS/PB_Cplacnjg.o \
PBLAS/SRC/PTOOLS/PB_Cplapad.o \
PBLAS/SRC/PTOOLS/PB_Cplapd2.o \
PBLAS/SRC/PTOOLS/PB_Cplaprnt.o \
PBLAS/SRC/PTOOLS/PB_Cplasca2.o \
PBLAS/SRC/PTOOLS/PB_Cplascal.o \
PBLAS/SRC/PTOOLS/PB_CpswapND.o \
PBLAS/SRC/PTOOLS/PB_CpswapNN.o \
PBLAS/SRC/PTOOLS/PB_Cpsym.o \
PBLAS/SRC/PTOOLS/PB_CpsymmAB.o \
PBLAS/SRC/PTOOLS/PB_CpsymmBC.o \
PBLAS/SRC/PTOOLS/PB_Cpsyr.o \
PBLAS/SRC/PTOOLS/PB_Cpsyr2.o \
PBLAS/SRC/PTOOLS/PB_Cpsyr2kA.o \
PBLAS/SRC/PTOOLS/PB_Cpsyr2kAC.o \
PBLAS/SRC/PTOOLS/PB_CpsyrkA.o \
PBLAS/SRC/PTOOLS/PB_CpsyrkAC.o \
PBLAS/SRC/PTOOLS/PB_Cptradd.o \
PBLAS/SRC/PTOOLS/PB_Cptran.o \
PBLAS/SRC/PTOOLS/PB_Cptrm.o \
PBLAS/SRC/PTOOLS/PB_CptrmmAB.o \
PBLAS/SRC/PTOOLS/PB_CptrmmB.o \
PBLAS/SRC/PTOOLS/PB_Cptrsm.o \
PBLAS/SRC/PTOOLS/PB_CptrsmAB.o \
PBLAS/SRC/PTOOLS/PB_CptrsmAB0.o \
PBLAS/SRC/PTOOLS/PB_CptrsmAB1.o \
PBLAS/SRC/PTOOLS/PB_CptrsmB.o \
PBLAS/SRC/PTOOLS/PB_Cptrsv.o \
PBLAS/SRC/PTOOLS/PB_Cspan.o \
PBLAS/SRC/PTOOLS/PB_Cstypeset.o \
PBLAS/SRC/PTOOLS/PB_Ctop.o \
PBLAS/SRC/PTOOLS/PB_Ctzahemv.o \
PBLAS/SRC/PTOOLS/PB_Ctzasymv.o \
PBLAS/SRC/PTOOLS/PB_Ctzatrmv.o \
PBLAS/SRC/PTOOLS/PB_Ctzhemm.o \
PBLAS/SRC/PTOOLS/PB_Ctzhemv.o \
PBLAS/SRC/PTOOLS/PB_Ctzher.o \
PBLAS/SRC/PTOOLS/PB_Ctzher2.o \
PBLAS/SRC/PTOOLS/PB_Ctzher2k.o \
PBLAS/SRC/PTOOLS/PB_Ctzherk.o \
PBLAS/SRC/PTOOLS/PB_Ctzsymm.o \
PBLAS/SRC/PTOOLS/PB_Ctzsymv.o \
PBLAS/SRC/PTOOLS/PB_Ctzsyr.o \
PBLAS/SRC/PTOOLS/PB_Ctzsyr2.o \
PBLAS/SRC/PTOOLS/PB_Ctzsyr2k.o \
PBLAS/SRC/PTOOLS/PB_Ctzsyrk.o \
PBLAS/SRC/PTOOLS/PB_Ctztrmm.o \
PBLAS/SRC/PTOOLS/PB_Ctztrmv.o \
PBLAS/SRC/PTOOLS/PB_Cwarn.o \
PBLAS/SRC/PTOOLS/PB_Cztypeset.o \
PBLAS/SRC/PTOOLS/PB_freebuf_.o \
PBLAS/SRC/PTOOLS/PB_topget_.o \
PBLAS/SRC/PTOOLS/PB_topset_.o \
PBLAS/SRC/PTZBLAS/cagemv.o \
PBLAS/SRC/PTZBLAS/cahemv.o \
PBLAS/SRC/PTZBLAS/casymv.o \
PBLAS/SRC/PTZBLAS/catrmv.o \
PBLAS/SRC/PTZBLAS/ccshft.o \
PBLAS/SRC/PTZBLAS/chescal.o \
PBLAS/SRC/PTZBLAS/cmmadd.o \
PBLAS/SRC/PTZBLAS/cmmcadd.o \
PBLAS/SRC/PTZBLAS/cmmdda.o \
PBLAS/SRC/PTZBLAS/cmmddac.o \
PBLAS/SRC/PTZBLAS/cmmddact.o \
PBLAS/SRC/PTZBLAS/cmmddat.o \
PBLAS/SRC/PTZBLAS/cmmtadd.o \
PBLAS/SRC/PTZBLAS/cmmtcadd.o \
PBLAS/SRC/PTZBLAS/crshft.o \
PBLAS/SRC/PTZBLAS/cset.o \
PBLAS/SRC/PTZBLAS/csymv.o \
PBLAS/SRC/PTZBLAS/csyr.o \
PBLAS/SRC/PTZBLAS/csyr2.o \
PBLAS/SRC/PTZBLAS/ctzcnjg.o \
PBLAS/SRC/PTZBLAS/ctzpad.o \
PBLAS/SRC/PTZBLAS/ctzpadcpy.o \
PBLAS/SRC/PTZBLAS/ctzscal.o \
PBLAS/SRC/PTZBLAS/cvvdotc.o \
PBLAS/SRC/PTZBLAS/cvvdotu.o \
PBLAS/SRC/PTZBLAS/dagemv.o \
PBLAS/SRC/PTZBLAS/dascal.o \
PBLAS/SRC/PTZBLAS/dasqrtb.o \
PBLAS/SRC/PTZBLAS/dasymv.o \
PBLAS/SRC/PTZBLAS/datrmv.o \
PBLAS/SRC/PTZBLAS/dcshft.o \
PBLAS/SRC/PTZBLAS/dmmadd.o \
PBLAS/SRC/PTZBLAS/dmmcadd.o \
PBLAS/SRC/PTZBLAS/dmmdda.o \
PBLAS/SRC/PTZBLAS/dmmddac.o \
PBLAS/SRC/PTZBLAS/dmmddact.o \
PBLAS/SRC/PTZBLAS/dmmddat.o \
PBLAS/SRC/PTZBLAS/dmmtadd.o \
PBLAS/SRC/PTZBLAS/dmmtcadd.o \
PBLAS/SRC/PTZBLAS/drshft.o \
PBLAS/SRC/PTZBLAS/dset.o \
PBLAS/SRC/PTZBLAS/dtzpad.o \
PBLAS/SRC/PTZBLAS/dtzpadcpy.o \
PBLAS/SRC/PTZBLAS/dtzscal.o \
PBLAS/SRC/PTZBLAS/dvasum.o \
PBLAS/SRC/PTZBLAS/dvvdot.o \
PBLAS/SRC/PTZBLAS/dzvasum.o \
PBLAS/SRC/PTZBLAS/immadd.o \
PBLAS/SRC/PTZBLAS/immdda.o \
PBLAS/SRC/PTZBLAS/immddat.o \
PBLAS/SRC/PTZBLAS/immtadd.o \
PBLAS/SRC/PTZBLAS/pxerbla.o \
PBLAS/SRC/PTZBLAS/sagemv.o \
PBLAS/SRC/PTZBLAS/sascal.o \
PBLAS/SRC/PTZBLAS/sasqrtb.o \
PBLAS/SRC/PTZBLAS/sasymv.o \
PBLAS/SRC/PTZBLAS/satrmv.o \
PBLAS/SRC/PTZBLAS/scshft.o \
PBLAS/SRC/PTZBLAS/scvasum.o \
PBLAS/SRC/PTZBLAS/smmadd.o \
PBLAS/SRC/PTZBLAS/smmcadd.o \
PBLAS/SRC/PTZBLAS/smmdda.o \
PBLAS/SRC/PTZBLAS/smmddac.o \
PBLAS/SRC/PTZBLAS/smmddact.o \
PBLAS/SRC/PTZBLAS/smmddat.o \
PBLAS/SRC/PTZBLAS/smmtadd.o \
PBLAS/SRC/PTZBLAS/smmtcadd.o \
PBLAS/SRC/PTZBLAS/srshft.o \
PBLAS/SRC/PTZBLAS/sset.o \
PBLAS/SRC/PTZBLAS/stzpad.o \
PBLAS/SRC/PTZBLAS/stzpadcpy.o \
PBLAS/SRC/PTZBLAS/stzscal.o \
PBLAS/SRC/PTZBLAS/svasum.o \
PBLAS/SRC/PTZBLAS/svvdot.o \
PBLAS/SRC/PTZBLAS/zagemv.o \
PBLAS/SRC/PTZBLAS/zahemv.o \
PBLAS/SRC/PTZBLAS/zasymv.o \
PBLAS/SRC/PTZBLAS/zatrmv.o \
PBLAS/SRC/PTZBLAS/zcshft.o \
PBLAS/SRC/PTZBLAS/zhescal.o \
PBLAS/SRC/PTZBLAS/zmmadd.o \
PBLAS/SRC/PTZBLAS/zmmcadd.o \
PBLAS/SRC/PTZBLAS/zmmdda.o \
PBLAS/SRC/PTZBLAS/zmmddac.o \
PBLAS/SRC/PTZBLAS/zmmddact.o \
PBLAS/SRC/PTZBLAS/zmmddat.o \
PBLAS/SRC/PTZBLAS/zmmtadd.o \
PBLAS/SRC/PTZBLAS/zmmtcadd.o \
PBLAS/SRC/PTZBLAS/zrshft.o \
PBLAS/SRC/PTZBLAS/zset.o \
PBLAS/SRC/PTZBLAS/zsymv.o \
PBLAS/SRC/PTZBLAS/zsyr.o \
PBLAS/SRC/PTZBLAS/zsyr2.o \
PBLAS/SRC/PTZBLAS/ztzcnjg.o \
PBLAS/SRC/PTZBLAS/ztzpad.o \
PBLAS/SRC/PTZBLAS/ztzpadcpy.o \
PBLAS/SRC/PTZBLAS/ztzscal.o \
PBLAS/SRC/PTZBLAS/zvvdotc.o \
PBLAS/SRC/PTZBLAS/zvvdotu.o
EOF
cat > Makefile.parallel <<EOF
include SLmake.inc
include Makefile.objs

objsblacsC = \$(objsblacs:.o=.oo)

lib: \${objslamov} \${objsblacs} \${objs}
	\$(ARCH) \$(ARCHFLAGS) \$(SCALAPACKLIB) \${objslamov} \${objsblacs} \${objs} \${objsblacsC}
	\$(RANLIB) \$(SCALAPACKLIB)

\${objslamov}:
	\$(CC) -o \$*.o -c \$(CFLAGS) \$(CDEFS) \$*.c

\${objsblacs}:
	\$(CC) -o \$*.oo -c \$(CDEFS) \$(CCFLAGS) -DCallFromC \$*.c
	\$(CC) -o \$*.o -c \$(CDEFS) \$(CCFLAGS) \$*.c

.f.o :
	\$(FC) -o \$*.o -c \$(FCFLAGS) \$*.f

.c.o :
	\$(CC) -o \$*.o -c \$(CDEFS) \$(CCFLAGS) \$*.c

cleanlib :
	rm -f \${objslamov} \${objsblacs} \${objs} \${objsblacsC}
	rm -f \$(SCALAPACKLIB)
EOF
cat > SLmake.inc <<EOF
CDEFS         = -DAdd_

#
#  The fortran and C compilers, loaders, and their flags
#

FC            = $MPI_F90_COMPILER
CC            = $MPI_C_COMPILER
NOOPT         = -O0
FCFLAGS       = -O3 $GCC_EXTRA_FFLAGS
CCFLAGS       = -O3
FCLOADER      = \$(FC)
CCLOADER      = \$(CC)
FCLOADFLAGS   = \$(FCFLAGS)
CCLOADFLAGS   = \$(CCFLAGS)

#
#  The archiver and the flag(s) to use when building archive (library)
#  Also the ranlib routine.  If your system has no ranlib, set RANLIB = echo
#

ARCH          = ar
ARCHFLAGS     = cr
RANLIB        = ranlib

#
#  The name of the ScaLAPACK library to be created
#

SCALAPACKLIB  = libscalapack.a

#
#  BLAS, LAPACK (and possibly other) libraries needed for linking test programs
#

BLASLIB       = $NON_INTEL_BLAS_LINK
LAPACKLIB     = $NON_INTEL_LAPACK_LINK
LIBS          = \$(LAPACKLIB) \$(BLASLIB)
EOF
        make -f Makefile.parallel -j$MAKE_JOBS
        mkdir -p $GOMA_LIB/scalapack-$SCALAPACK_VERSION/lib
        cp -a libscalapack.a $GOMA_LIB/scalapack-$SCALAPACK_VERSION/lib
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
OPTF    = -O2  -DALLOW_NON_INIT $GCC_EXTRA_FFLAGS
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
    make -j$MAKE_JOBS prerequisites 2>&1 | tee -a $COMPILE_LOG
    make -j$MAKE_JOBS all 2>&1 | tee -a $COMPILE_LOG

    if [ -e $GOMA_LIB/MUMPS_$MUMPS_VERSION/lib/libdmumps.a ]
    then
        log_echo "Built MUMPS $MUMPS_VERSION"
    else
        log_echo "Failed to build MUMP $MUMPS_VERSION"
        exit 1
    fi
fi

if [ -e $GOMA_LIB/Omega_h/lib/libomega_h.a ]; then
    log_echo "Omega_h is already built!"
else
    cd $GOMA_LIB/omega_h-$OMEGA_H_VERSION
    CC=mpicc CXX=mpicxx cmake -Bbuild -DCMAKE_INSTALL_PREFIX=$GOMA_LIB/Omega_h -DBUILD_SHARED_LIBS=OFF -DOmega_h_USE_MPI=ON 2>&1 | tee -a $COMPILE_LOG
    cmake --build build -j$MAKE_JOBS 2>&1 | tee -a $COMPILE_LOG
    make -C build install 2>&1 | tee -a $COMPILE_LOG
    if [ -e $GOMA_LIB/Omega_h/lib/libomega_h.a ]; then
        log_echo "Omega_h built!"
    else
        log_echo "Failed to build Omega_H $OMEGA_H_VERSION"
        exit 1
    fi
    cd $GOMA_LIB
fi

#make trilinos



cd $GOMA_LIB

rm -rf $GOMA_LIB/trilinos-$TRILINOS_VERSION-Temp
mkdir $GOMA_LIB/trilinos-$TRILINOS_VERSION-Temp

rm -f CMakeCache.txt

MPI_LIBS="-lmpi ${MPI_FORTRAN_LIB}"

HDF5_LIBS="-L${GOMA_LIB}/hdf5-${HDF5_VERSION}/lib -lhdf5_hl -lhdf5 -lz -ldl"
# Install directory
TRILINOS_INSTALL=$GOMA_LIB/trilinos-$TRILINOS_VERSION

if [ -e $TRILINOS_INSTALL/bin/aprepro ]; then
    log_echo "Trilinos is already built!"
else
cd $GOMA_LIB/Trilinos-trilinos-release-$TRILINOS_VERSION_DASH/packages/aztecoo/src
cat << "EOF" > az_aztec_h.patch 
731c731
<   extern char *AZ_allocate(unsigned int iii);
---
>   extern char *AZ_allocate(size_t iii);
1024c1024
<   extern double *AZ_manage_memory(unsigned int size, int action, int type,
---
>   extern double *AZ_manage_memory(size_t size, int action, int type,
1198c1198
<   extern char *AZ_realloc(void *ptr, unsigned int size);
---
>   extern char *AZ_realloc(void *ptr, size_t size);
EOF
patch -f az_aztec.h < az_aztec_h.patch
rm az_aztec_h.patch

cat << "EOF" > az_util_c.patch 
843c843
<   (void) AZ_manage_memory((unsigned int) 0, AZ_CLEAR, label, (char *) NULL,
---
>   (void) AZ_manage_memory((size_t) 0, AZ_CLEAR, label, (char *) NULL,
851c851
< double *AZ_manage_memory(unsigned int input_size, int action, int type,
---
> double *AZ_manage_memory(size_t input_size, int action, int type,
936c936
<     int     size;
---
>     size_t  size;
941c941
<   long int size;
---
>   size_t                size, aligned_size;
945,946c945,946
<   long int aligned_str_mem, aligned_j, aligned_size;
< double *dtmp;
---
>   unsigned int          aligned_str_mem, aligned_j;
>   double *dtmp;
949c949
<   size = (long int) input_size;
---
>   size = input_size;
997,998c997
<       dtmp = (double *) AZ_allocate((unsigned int) (aligned_str_mem+aligned_j+
<                                                 aligned_size) );
---
>       dtmp = (double *) AZ_allocate(aligned_str_mem+aligned_j+aligned_size);
1183,1184c1182
<     dtmp    = (double *) AZ_realloc((char *) dtmp,(unsigned int)
<                                     aligned_str_mem+aligned_j+aligned_size);
---
>     dtmp    = (double *) AZ_realloc((char *) dtmp, aligned_str_mem+aligned_j+aligned_size);
1872c1870
<    int size;
---
>    size_t size;
1902c1900
< char *AZ_allocate(unsigned int isize) {
---
> char *AZ_allocate(size_t isize) {
1917,1918c1915,1917
<     int *size_ptr, i;
<     unsigned int size;
---
>     int i;
>     size_t size;
>     size_t *size_ptr; 
1952c1951
<     size_ptr = (int *) ptr;
---
>     size_ptr = (size_t *) ptr;
1992c1991
<    int *iptr, size, i;
---
>    int i;
1993a1993,1994
>    size_t size;
>    size_t* iptr;
2023c2024
<            iptr = (int *) ptr;
---
>            iptr = (size_t *) ptr;
2073c2074
< char *AZ_realloc(void *vptr, unsigned int new_size) {
---
> char *AZ_realloc(void *vptr, size_t new_size) {
2076c2077
<    int i, *iptr, size, *new_size_ptr;
---
>    int i;
2079d2079
<    int newmsize, smaller;
2080a2081,2082
>    size_t size, newmsize, smaller;
>    size_t *iptr, *new_size_ptr; 
2108c2110
<            iptr = (int *) ptr;
---
>            iptr = (size_t *) ptr;
2128c2130
<     new_size_ptr = (int *) new_ptr;
---
>     new_size_ptr = (size_t *) new_ptr;
2175c2177
< char *AZ_allocate(unsigned int size) {
---
> char *AZ_allocate(size_t size) {
2185c2187
< char *AZ_realloc(void *ptr, unsigned int size) {
---
> char *AZ_realloc(void *ptr, size_t size) {
EOF
patch -f az_util.c < az_util_c.patch
rm az_util_c.patch 
cd $GOMA_LIB/trilinos-$TRILINOS_VERSION-Temp
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
-D Trilinos_ENABLE_Triutils:BOOL=ON \
-D Trilinos_ENABLE_SEACAS:BOOL=ON \
-D Trilinos_GOMA_ENABLE_AMESOS:BOOL=ON \
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
-D Trilinos_GOMA_ENABLE_AMESOS2:BOOL=ON \
-D Trilinos_ENABLE_Belos:BOOL=ON \
-D Trilinos_ENABLE_EpetraExt:BOOL=ON \
-D Trilinos_ENABLE_Thyra:BOOL=ON \
-D Trilinos_ENABLE_ThyraTpetraAdapters:BOOL=ON \
-D Trilinos_ENABLE_Tpetra:BOOL=ON \
-D Trilinos_ENABLE_Stratimikos:BOOL=ON \
-D Trilinos_ENABLE_TESTS:BOOL=OFF \
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
  -D UMFPACK_LIBRARY_NAMES:STRING="umfpack;amd;suitesparseconfig;cholmod;colamd;ccolamd;camd" \
  -D UMFPACK_LIBRARY_DIRS:PATH="$GOMA_LIB/SuiteSparse-$SUITESPARSE_VERSION/UMFPACK/Lib;$GOMA_LIB/SuiteSparse-$SUITESPARSE_VERSION/CHOLMOD/Lib;$GOMA_LIB/SuiteSparse-$SUITESPARSE_VERSION/AMD/Lib;$GOMA_LIB/SuiteSparse-$SUITESPARSE_VERSION/SuiteSparse_config;$GOMA_LIB/SuiteSparse-$SUITESPARSE_VERSION/CAMD/Lib;$GOMA_LIB/SuiteSparse-$SUITESPARSE_VERSION/COLAMD/Lib;$GOMA_LIB/SuiteSparse-$SUITESPARSE_VERSION/CCOLAMD/Lib" \
  -D UMFPACK_INCLUDE_DIRS:PATH="$GOMA_LIB/SuiteSparse-$SUITESPARSE_VERSION/include" \
-D TPL_ENABLE_AMD:BOOL=ON \
  -D AMD_LIBRARY_NAMES:STRING="amd;suitesparseconfig" \
  -D AMD_LIBRARY_DIRS:PATH="$GOMA_LIB/SuiteSparse-$SUITESPARSE_VERSION/AMD/Lib;$GOMA_LIB/SuiteSparse-$SUITESPARSE_VERSION/SuiteSparse_config" \
  -D AMD_INCLUDE_DIRS:PATH="$GOMA_LIB/SuiteSparse-$SUITESPARSE_VERSION/AMD/Include;$GOMA_LIB/SuiteSparse-$SUITESPARSE_VERSION/SuiteSparse_config" \
-D TPL_ENABLE_SuperLUDist:BOOL=ON \
  -D SuperLUDist_LIBRARY_NAMES:STRING="superlu_dist" \
  -D SuperLUDist_LIBRARY_DIRS:PATH=$GOMA_LIB/superlu_dist-$SUPERLU_DIST_VERSION/$SUPERLU_LIBDIR \
  -D SuperLUDist_INCLUDE_DIRS:PATH=$GOMA_LIB/superlu_dist-$SUPERLU_DIST_VERSION/include \
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
  -D CMAKE_Fortran_FLAGS:STRING="$GCC_EXTRA_FFLAGS" \
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

cd $GOMA_LIB/petsc-$PETSC_VERSION
export PETSC_DIR=$GOMA_LIB/petsc-$PETSC_VERSION
export PETSC_ARCH=arch-linux-c-opt

if [ -e $PETSC_DIR/$PETSC_ARCH/lib/libpetsc.a ]; then
    log_echo "PETSc is already built!"
else
    ./configure --with-shared-libraries=0 --with-cc=$(which mpicc) --with-cxx=$(which mpicxx) --with-fc=$(which mpif90) --with-debugging=0 COPTFLAGS='-O3' CXXOPTFLAGS='-O3' FOPTFLAGS='-O3' --download-hypre --with-scalapack=1 --with-scalapack-dir=$(readlink --canonicalize-missing "${SCALAPACK_LIBRARY_DIR}/..") --with-superlu_dist=1 --with-superlu_dist-dir=$GOMA_LIB/superlu_dist-$SUPERLU_DIST_VERSION --with-metis=1 --with-metis-dir=$GOMA_LIB/parmetis-4.0.3 --with-parmetis=1 --with-parmetis-dir=$GOMA_LIB/parmetis-4.0.3 --with-blas-lib=${NON_INTEL_BLAS_LIBRARY} --with-lapack-lib=${NON_INTEL_LAPACK_LIBRARY} --with-mumps=1 --with-mumps-dir="$GOMA_LIB/MUMPS_$MUMPS_VERSION"  2>&1 | tee -a $COMPILE_LOG
    make -j$MAKE_JOBS all 2>&1 | tee -a $COMPILE_LOG
    make check 2>&1 | tee -a $COMPILE_LOG
    if [ -e $PETSC_DIR/$PETSC_ARCH/lib/libpetsc.a ]; then
        log_echo "PETSc built!"
    else
        log_echo "Failed to build PETSc $PETSC_VERSION"
        exit 1
    fi
fi


# Generate a config file for bash
cd $GOMA_LIB
cat > config.sh <<EOF
export CMAKE_PREFIX_PATH=$TRILINOS_INSTALL:${GOMA_LIB}/Omega_h:\$CMAKE_PREFIX_PATH
export PATH=$TRILINOS_INSTALL/bin:\$PATH
export PATH=${GOMA_LIB}/openmpi-$OPENMPI_VERSION/bin:\$PATH
export METISDIR=${GOMA_LIB}/parmetis-4.0.3
export UMFPACK_DIR=${GOMA_LIB}/SuiteSparse-$SUITESPARSE_VERSION
export ARPACKDIR=${GOMA_LIB}/ARPACK
export SPARSEDIR=${GOMA_LIB}/sparse
export PETSC_DIR=${PETSC_DIR}
export PETSC_ARCH=${PETSC_ARCH}
EOF


if [ "$build_cmake" == "true" ] ; then
    echo "export PATH=$GOMA_LIB/cmake-$CMAKE_VERSION-linux-x86_64/bin:\$PATH" >> config.sh
fi

log_echo
log_echo "An example bash configuration file has been written to $GOMA_LIB/config.sh"
log_echo
log_echo "Activate with $ source $GOMA_LIB/config.sh"
