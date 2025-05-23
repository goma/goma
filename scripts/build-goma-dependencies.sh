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
# if using GCC >= 10 uncomment:
# export GCC_EXTRA_FFLAGS="-fallow-argument-mismatch"

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


echo "There is a new TPL build script in goma/tpls/ called install-tpls.py"
echo "Please transition to using that script."
echo "This script will be removed in a future release."
continue_check


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

NETCDF_VERSION="c-4.9.0"
NETCDF_VERSION_ONLY="4.9.0"
NETCDF_MD5="26cfd1b2e32d511bd82477803976365d"

PNETCDF_VERSION="1.13.0"
PNETCDF_MD5="31b94d39462b1f1f2293f735c9819bf2"
PNETCDF_URL="https://parallel-netcdf.github.io/Release/pnetcdf-$PNETCDF_VERSION.tar.gz"

TRILINOS_VERSION="15.1.0"
TRILINOS_VERSION_DASH="15-1-0"
TRILINOS_MD5="79237697af4fc42eaaf70f23104a8e12"

MUMPS_VERSION="5.6.2"
MUMPS_URL="https://graal.ens-lyon.fr/MUMPS/MUMPS_$MUMPS_VERSION.tar.gz"
MUMPS_MD5="adc38b6cca34dfcc8f7036b7d39cf260"

OPENMPI_VERSION="4.1.6"
OPENMPI_MD5="c9b1c974cfc23c77c0fbdb965cd58a1c"
OPENMPI_ARCHIVE_URL="https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-$OPENMPI_VERSION.tar.bz2"
OPENMPI_EXTRA_CONFIGURE_FLAGS=""

CMAKE_VERSION="3.24.2"
CMAKE_MD5="f6f79ec66c91bbc075757a70205128ca"

#SUITESPARSE_VERSION="7.7.0"
#SUITESPARSE_MD5="e659373ed5e9b961d2fcb6d67d250783"
SUITESPARSE_VERSION="5.13.0"
SUITESPARSE_MD5="e9e7bc594b77ae4b58d943cdc286d679"

MATIO_VERSION="1.5.21"
MATIO_MD5="afeb5d21b234699fd5b9dc4564afe1ca"

SCALAPACK_VERSION="0234af94c6578c53ac4c19f2925eb6e5c4ad6f0f"
SCALAPACK_MD5="221dab094055b73fa7bed645c426548c"

LAPACK_VERSION="3.8.0"
LAPACK_MD5="96591affdbf58c450d45c1daa540dbd2"

PETSC_VERSION="3.20.2"
PETSC_MD5="1e170a5f096433ca21aa643c80c749eb"

OMEGA_H_VERSION="9.34.13"
OMEGA_H_MD5="b6eee23870bb21b58f9bbf0f640bf12a"

SUPERLU_DIST_VERSION="7.2.0"
SUPERLU_DIST_MD5="10d20b97012e5ae89a6cc69b768f61b7"

OPENBLAS_VERSION="0.3.21"
OPENBLAS_MD5="ffb6120e2309a2280471716301824805"
OPENBLAS_URL="https://github.com/xianyi/OpenBLAS/releases/download/v$OPENBLAS_VERSION/OpenBLAS-$OPENBLAS_VERSION.tar.gz"

ARPACK_NG_VERSION="3.8.0"
ARPACK_NG_MD5="bb4cf448f2480a0ffe5517d579f980c3"
ARPACK_NG_URL="https://github.com/opencollab/arpack-ng/archive/refs/tags/$ARPACK_NG_VERSION.tar.gz"

METIS_VERSION="5.1.0-p10"
METIS_MD5="82e3b09d577580747fbcada30d9b4228"
METIS_EXTRACT_NAME="petsc-pkg-metis-c8d2dc1e751e"

PARMETIS_VERSION="4.0.3-p8"
PARMETIS_MD5="5549747a1a653f1e1909d40fff52ddb9"
PARMETIS_EXTRACT_NAME="petsc-pkg-parmetis-5777d7ec2084"

ARCHIVE_NAMES=("arpack-ng-$ARPACK_NG_VERSION.tar.gz" \
"hdf5-${HDF5_VERSION}.tar.bz2" \
"pnetcdf-$PNETCDF_VERSION.tar.gz" \
"netcdf-${NETCDF_VERSION}.tar.gz" \
"metis-$METIS_VERSION.tar.gz" \
"parmetis-$PARMETIS_VERSION.tar.gz" \
"sparse.tar.gz" \
"superlu_dist-$SUPERLU_DIST_VERSION.tar.gz" \
"Trilinos-trilinos-release-$TRILINOS_VERSION_DASH.tar.gz" \
"MUMPS_$MUMPS_VERSION.tar.gz" \
"SuiteSparse-$SUITESPARSE_VERSION.tar.gz" \
"matio-$MATIO_VERSION.tar.gz"  \
"petsc-$PETSC_VERSION.tar.gz" \
"omega_h-v$OMEGA_H_VERSION.tar.gz")

#y12m archive is skipped because it stores the number of downloads in the header

#meaning each y12m tar has a unique MD5SUM.
ARCHIVE_MD5SUMS=("$ARPACK_NG_MD5" \
"${HDF5_MD5}" \
"$PNETCDF_MD5" \
"${NETCDF_MD5}" \
"$METIS_MD5" \
"$PARMETIS_MD5" \
"1566d914d1035ac17b73fe9bc0eed02a" \
"$SUPERLU_DIST_MD5" \
$TRILINOS_MD5 \
$MUMPS_MD5 \
"$SUITESPARSE_MD5" \
"$MATIO_MD5" \
"$PETSC_MD5" \
"$OMEGA_H_MD5")

ARCHIVE_URLS=("$ARPACK_NG_URL" \
"https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-${HDF5_VERSION}/src/hdf5-${HDF5_VERSION}.tar.bz2" \
"$PNETCDF_URL" \
"https://downloads.unidata.ucar.edu/netcdf-c/${NETCDF_VERSION_ONLY}/netcdf-${NETCDF_VERSION}.tar.gz" \
"https://bitbucket.org/petsc/pkg-metis/get/v$METIS_VERSION.tar.gz" \
"https://bitbucket.org/petsc/pkg-parmetis/get/v$PARMETIS_VERSION.tar.gz" \
"http://downloads.sourceforge.net/project/sparse/sparse/sparse1.4b/sparse1.4b.tar.gz" \
"http://codeload.github.com/xiaoyeli/superlu_dist/tar.gz/v$SUPERLU_DIST_VERSION" \
"https://github.com/trilinos/Trilinos/archive/trilinos-release-$TRILINOS_VERSION_DASH.tar.gz" \
"$MUMPS_URL" \
"https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v$SUITESPARSE_VERSION.tar.gz" \
"https://github.com/tbeu/matio/releases/download/v$MATIO_VERSION/matio-$MATIO_VERSION.tar.gz" \
"https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-$PETSC_VERSION.tar.gz" \
"https://github.com/sandialabs/omega_h/archive/refs/tags/v$OMEGA_H_VERSION.tar.gz")

# You can't call the ARPACK patch ARPACK or it will think it is already extracted
# When in reality it isn't
ARCHIVE_DIR_NAMES=("arpack-ng-${ARPACK_NG_VERSION}" \
"hdf5-${HDF5_VERSION}" \
"pnetcdf-$PNETCDF_VERSION" \
"netcdf-${NETCDF_VERSION}" \
"metis-$METIS_VERSION" \
"parmetis-$PARMETIS_VERSION" \
"sparse" \
"superlu_dist-$SUPERLU_DIST_VERSION" \
"Trilinos-trilinos-release-$TRILINOS_VERSION_DASH" \
"MUMPS_$MUMPS_VERSION" \
"SuiteSparse-$SUITESPARSE_VERSION" \
"matio-$MATIO_VERSION" \
"petsc-$PETSC_VERSION" \
"omega_h-$OMEGA_H_VERSION")

ARCHIVE_HOMEPAGES=("https://github.com/opencollab/arpack-ng/" \
"https://www.hdfgroup.org/" \
"https://parallel-netcdf.github.io/" \
"https://www.unidata.ucar.edu/software/netcdf/" \
"http://glaros.dtc.umn.edu/gkhome/metis/metis/overview" \
"http://glaros.dtc.umn.edu/gkhome/metis/metis/overview" \
"http://sparse.sourceforge.net/" \
"https://github.com/xiaoyeli/superlu_dist" \
"https://trilinos.org/" \
"http://mumps.enseeiht.fr/" \
"http://faculty.cse.tamu.edu/davis/suitesparse.html" \
"https://sourceforge.net/projects/matio/" \
"https://petsc.org" \
"https://github.com/sandialabs/omega_h/")

ARCHIVE_REAL_NAMES=("arpack-ng" \
"HDF5" \
"PNetCDF" \
"NetCDF" \
"METIS" \
"ParMETIS" \
"Sparse" \
"SuperLU_DIST" \
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
    if [[ "$cmake_vers" = $(echo -e "$cmake_vers\n3.23.0\n" | sort -V |tail -n1) ]]; then
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
        export GCC_EXTRA_FFLAGS=""
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
                ARCHIVE_NAMES+=("$SCALAPACK_VERSION.tar.gz")
                ARCHIVE_MD5SUMS+=("$SCALAPACK_MD5")
                ARCHIVE_URLS+=("https://github.com/Reference-ScaLAPACK/scalapack/archive/$SCALAPACK_VERSION.tar.gz")
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
        "$SCALAPACK_VERSION.tar.gz")
        ARCHIVE_MD5SUMS+=("$LAPACK_MD5" \
        "$SCALAPACK_MD5" )
        ARCHIVE_URLS+=("https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v$LAPACK_VERSION.tar.gz" \
        "https://github.com/Reference-ScaLAPACK/scalapack/archive/$SCALAPACK_VERSION.tar.gz" )
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
    elif [[ "$MATH_LIBRARIES" == "openblas" ]]; then
        USING_MKLS="OFF"
        #Add packages that otherwise come preinstalled in intel compiler.
        ARCHIVE_NAMES+=("OpenBLAS-$OPENBLAS_VERSION.tar.gz" \
        "$SCALAPACK_VERSION.tar.gz")
        ARCHIVE_MD5SUMS+=("$OPENBLAS_MD5" \
        "$SCALAPACK_MD5" )
        ARCHIVE_URLS+=("$OPENBLAS_URL" \
        "https://github.com/Reference-ScaLAPACK/scalapack/archive/$SCALAPACK_VERSION.tar.gz" )
        ARCHIVE_DIR_NAMES+=("OpenBLAS-$OPENBLAS_VERSION" \
        "scalapack-$SCALAPACK_VERSION" )
        ARCHIVE_HOMEPAGES+=("https://github.com/xianyi/OpenBLAS/" \
        "http://www.netlib.org/scalapack/")
        ARCHIVE_REAL_NAMES+=("OpenBLAS" \
        "ScaLAPACK")

        LAPACK_DIR="${GOMA_LIB}/OpenBLAS-$OPENBLAS_VERSION"
        LAPACK_LIBRARY_DIR="$LAPACK_DIR/lib"
        LAPACK_LIBRARY_NAME_ARG="-L${LAPACK_LIBRARY_DIR} -lopenblas"
        LAPACK_LIBRARY_NAME="openblas"
        BLAS_LIBRARY_DIR="${LAPACK_LIBRARY_DIR}"
        BLAS_LIBRARY_NAME_ARG="-L${BLAS_LIBRARY_DIR} -lopenblas"
        BLAS_LIBRARY_NAME="openblas"
        BLACS_LIBRARY_NAME="scalapack"
        SCALAPACK_LIBRARY_DIR="${GOMA_LIB}/scalapack-$SCALAPACK_VERSION/lib"
        SCALAPACK_LIBRARY_NAME="scalapack"
        SCALAPACK_LIBRARY_NAME_ARG="-L${SCALAPACK_LIBRARY_DIR} -lscalapack"
        SCALAPACK_INCLUDE_DIR="${GOMA_LIB}/scalapack-$SCALAPACK_VERSION/include"
        # Alternative to the intel "blas_flags" variable
        NON_INTEL_BLAS_LIBRARY="${BLAS_LIBRARY_DIR}/libopenblas.a"
        NON_INTEL_LAPACK_LIBRARY="${LAPACK_LIBRARY_DIR}/libopenblas.a"
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
        ARCHIVE_NAMES+=("$SCALAPACK_VERSION.tar.gz")
        ARCHIVE_MD5SUMS+=("$SCALAPACK_MD5" )
        ARCHIVE_URLS+=("https://github.com/Reference-ScaLAPACK/scalapack/archive/$SCALAPACK_VERSION.tar.gz" )
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

# PNetCDF
cd $GOMA_LIB
if [ -e pnetcdf-${PNETCDF_VERSION}/lib/libpnetcdf.a ]
then
    log_echo "PNetCDF already built."
else
    if ! [ -e pnetcdf-${NETCDF_VERSION}/.goma-extracted ]
    then
    mv pnetcdf-${PNETCDF_VERSION} tmpdir
    mkdir pnetcdf-${PNETCDF_VERSION}
    mv tmpdir pnetcdf-${PNETCDF_VERSION}/src
    touch pnetcdf-${PNETCDF_VERSION}/.goma-extracted
    fi
    cd pnetcdf-${PNETCDF_VERSION}/src
    CC=${MPI_C_COMPILER} CXX=${MPI_CXX_COMPILER} FC=${MPI_F90_COMPILER} ./configure --disable-fortran --enable-shared=off --disable-cxx --prefix=$GOMA_LIB/pnetcdf-${PNETCDF_VERSION} 2>&1 | tee -a $COMPILE_LOG
    make -j$MAKE_JOBS 2>&1 | tee -a $COMPILE_LOG
    make install 2>&1 | tee -a $COMPILE_LOG
    if [ -e $GOMA_LIB/pnetcdf-${PNETCDF_VERSION}/lib/libpnetcdf.a ]
    then
        log_echo "Built PNetCDF $PNETCDF_VERSION"
    else
        log_echo "Failed to build PNetCDF $PNETCDF_VERSION"
        exit 1
    fi
fi


# NetCDF
cd $GOMA_LIB
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
    export CPPFLAGS="-I$GOMA_LIB/hdf5-${HDF5_VERSION}/include -I$GOMA_LIB/pnetcdf-${PNETCDF_VERSION}/include"
    export LDFLAGS="-L${GOMA_LIB}/hdf5-${HDF5_VERSION}/lib -L${GOMA_LIB}/pnetcdf-${PNETCDF_VERSION}/lib" 
    log_echo $CPPFLAGS
    log_echo $LDFLAGS

    CC=${MPI_C_COMPILER} CFLAGS="-I${GOMA_LIB}/hdf5-${HDF5_VERSION}/include -I${GOMA_LIB}/pnetcdf-${PNETCDF_VERSION}/include" \
      CPP="${MPI_C_COMPILER} -E" \
      CPPFLAGS="-I${GOMA_LIB}/hdf5-${HDF5_VERSION}/include -I${GOMA_LIB}/pnetcdf-${PNETCDF_VERSION}/include" \
      LDFLAGS="-L${GOMA_LIB}/hdf5-${HDF5_VERSION}/lib -L${GOMA_LIB}/pnetcdf-${PNETCDF_VERSION}/lib" \
      ./configure \
      --prefix=$GOMA_LIB/netcdf-${NETCDF_VERSION} \
      --enable-shared=off \
      --enable-pnetcdf \
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
fi

#parmetis patch
cd $GOMA_LIB
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

#make Metis
cd $GOMA_LIB
if [ -e metis-$METIS_VERSION/lib/libmetis.a ]
then
    log_echo "Metis already Built"
else
    cd $METIS_EXTRACT_NAME
    mkdir build
    cd build
    CC=$MPI_C_COMPILER CXX=$MPI_CXX_COMPILER FC=$MPI_F90_COMPILER cmake .. -DGKLIB_PATH=../GKlib -DSHARED=0  -DCMAKE_INSTALL_PREFIX=$GOMA_LIB/metis-$METIS_VERSION | tee -a $COMPILE_LOG
    make 2>&1 | tee -a $COMPILE_LOG
    make install 2>&1 | tee -a $COMPILE_LOG
    cd ..
    if [ -e $GOMA_LIB/metis-$METIS_VERSION/lib/libmetis.a ]
    then
        log_echo "Built Metis $METIS_VERSION"
    else
        log_echo "Failed to build Metis $METIS_VERSION"
        exit 1
    fi
fi
#make parMetis
cd $GOMA_LIB
if [ -e parmetis-$PARMETIS_VERSION/lib/libparmetis.a ]
then
    log_echo "ParMetis already Built"
else
    cd $PARMETIS_EXTRACT_NAME
    mkdir build
    cd build
    CC=$MPI_C_COMPILER CXX=$MPI_CXX_COMPILER FC=$MPI_F90_COMPILER cmake .. -DGKLIB_PATH=../headers -DMETIS_PATH=$GOMA_LIB/metis-$METIS_VERSION -DSHARED=0 -DMPI_INCLUDE_PATH=$MPI_BASE_DIR/include -DCMAKE_INSTALL_PREFIX=$GOMA_LIB/parmetis-$PARMETIS_VERSION | tee -a $COMPILE_LOG
    make 2>&1 | tee -a $COMPILE_LOG
    make install 2>&1 | tee -a $COMPILE_LOG
    cd ..
    if [ -e $GOMA_LIB/parmetis-$PARMETIS_VERSION/lib/libparmetis.a ]
    then
        log_echo "Built ParMetis 4.0.3"
    else
        log_echo "Failed to build ParMetis 4.0.3"
        exit 1
    fi
fi


#make ARPACK
#cd $GOMA_LIB/ARPACK
#if [ -e lib/libarpack.a ]
#then
#    log_echo "ARPACK already built"
#else
#    cat > ARmake.patch << EOF
#28c28
#< home = \$(HOME)/ARPACK
#---
#> home = $GOMA_LIB/ARPACK
#35c35
#< PLAT = SUN4
#---
#> PLAT = x86_64
#104,105c104,105
#< FC      = f77
#< FFLAGS	= -O -cg89
#---
#> FC      = $MPI_F77_COMPILER
#> FFLAGS	= -O $BLAS_FLAGS $GCC_EXTRA_FFLAGS
#115c115
#< MAKE    = /bin/make
#---
#> MAKE    = make
#EOF
#
#    patch ARmake.inc < ARmake.patch
#    # Documentation says this should always be needed but in reality only intel MKL requires it.
#    cat > UTIL/second.f.patch << EOF
#24c24
#<       EXTERNAL           ETIME
#---
#> *     EXTERNAL           ETIME
#EOF
#    patch UTIL/second.f < UTIL/second.f.patch
#    make all 2>&1 | tee -a $COMPILE_LOG
#    mkdir lib
#    cp libarpack_x86_64.a lib/libarpack.a
#    if [ -e $GOMA_LIB/ARPACK/lib/libarpack.a ]
#    then
#        log_echo "Built ARPACK"
#    else
#        log_echo "Failed to build ARPACK"
#        exit 1
#    fi
#fi

cd $GOMA_LIB

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

if [[ "$MATH_LIBRARIES" == "openblas" ]]; then
    #make lapack/blas
    if [ -e $LAPACK_DIR/lib/libopenblas.a ]
    then
        log_echo "OpenBLAS already built"
    else
        tempdir=$(mktemp -u -p $GOMA_LIB)
        mv $LAPACK_DIR $tempdir
	mkdir $LAPACK_DIR
	mv $tempdir $LAPACK_DIR/src
        cd $LAPACK_DIR/src
        CC=$MPI_C_COMPILER CXX=$MPI_CXX_COMPILER FC=$MPI_F90_COMPILER make USE_THREAD=0 USE_OPENMP=0 NO_SHARED=1 2>&1 | tee -a $COMPILE_LOG
	make PREFIX=$LAPACK_DIR USE_THREAD=0 USE_OPENMP=0 NO_SHARED=1 install 2>&1 | tee -a $COMPILE_LOG
	cd $GOMA_LIB
        if [ -e $LAPACK_LIBRARY_DIR/libopenblas.a ]
        then
            log_echo "Built OpenBLAS"
        else
            log_echo "Failed to build OpenBLAS"
            exit 1
        fi
    fi
    export LD_LIBRARY_PATH="${LAPACK_LIBRARY_DIR}:$LD_LIBRARY_PATH"
fi



#make arpack
if [ -e arpack-ng-$ARPACK_NG_VERSION/lib/libarpack.a ]
then
    log_echo "Arpack-ng already built"
else
    tempdir=$(mktemp -u -p $GOMA_LIB)
    mv arpack-ng-$ARPACK_NG_VERSION $tempdir
    mkdir arpack-ng-$ARPACK_NG_VERSION
    mv $tempdir arpack-ng-$ARPACK_NG_VERSION/src
    cd arpack-ng-$ARPACK_NG_VERSION/src/
cat << "EOF" > arpack-cmake.patch
--- CMakeLists.txt      2020-12-07 03:40:45.000000000 -0700
+++ CMakeLists.txt-blasfix      2023-01-23 12:27:17.893954624 -0700
@@ -58,7 +58,7 @@
     foreach(l ${${list_name}})
         get_filename_component(lwe ${l} NAME_WE)
         add_executable(${lwe} ${arpackexample_DIR}/${l} ${examples_EXTRA_SRCS})
-       target_link_libraries(${lwe} arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
+       target_link_libraries(${lwe} arpack LAPACK::LAPACK BLAS::BLAS ${EXTRA_LDFLAGS})
         add_test(NAME "${lwe}_ex" COMMAND ${lwe} WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
     endforeach()
 endfunction(examples)
@@ -306,8 +306,8 @@
 # use -DBUILD_SHARED_LIBS=ON|OFF to control static/shared
 add_library(arpack ${arpackutil_STAT_SRCS} ${arpacksrc_STAT_SRCS} ${arpacksrc_ICB})

-target_link_libraries(arpack BLAS::BLAS)
 target_link_libraries(arpack LAPACK::LAPACK ${EXTRA_LDFLAGS})
+target_link_libraries(arpack BLAS::BLAS)
 set_target_properties(arpack PROPERTIES OUTPUT_NAME arpack${LIBSUFFIX})
 set_target_properties(arpack PROPERTIES VERSION 2.1.0)
 set_target_properties(arpack PROPERTIES SOVERSION 2)
@@ -529,7 +529,7 @@

 add_executable(dnsimp_test TESTS/dnsimp.f TESTS/mmio.f TESTS/debug.h)
 set_target_properties( dnsimp_test PROPERTIES OUTPUT_NAME  dnsimp )
-target_link_libraries(dnsimp_test arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
+target_link_libraries(dnsimp_test arpack LAPACK::LAPACK BLAS::BLAS ${EXTRA_LDFLAGS})
 add_custom_command(TARGET dnsimp_test POST_BUILD
   COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/TESTS/testA.mtx testA.mtx
 )
@@ -537,39 +537,39 @@

 if (NOT ICB)
     add_executable(bug_1315_single TESTS/bug_1315_single.c)
-    target_link_libraries(bug_1315_single arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
+    target_link_libraries(bug_1315_single arpack LAPACK::LAPACK BLAS::BLAS ${EXTRA_LDFLAGS})
     add_test(bug_1315_single_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bug_1315_single)

     add_executable(bug_1315_double TESTS/bug_1315_double.c)
-    target_link_libraries(bug_1315_double arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
+    target_link_libraries(bug_1315_double arpack LAPACK::LAPACK BLAS::BLAS ${EXTRA_LDFLAGS})
     add_test(bug_1315_double_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bug_1315_double)
 endif()

 add_executable(bug_1323 TESTS/bug_1323.f)
-target_link_libraries(bug_1323 arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
+target_link_libraries(bug_1323 arpack LAPACK::LAPACK BLAS::BLAS ${EXTRA_LDFLAGS})
 add_test(bug_1323_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bug_1323)

 add_executable(bug_58_double TESTS/bug_58_double.f)
-target_link_libraries(bug_58_double arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
+target_link_libraries(bug_58_double arpack LAPACK::LAPACK BLAS::BLAS ${EXTRA_LDFLAGS})
 add_test(bug_58_double_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bug_58_double)

 add_executable(bug_79_double_complex TESTS/bug_79_double_complex.f)
-target_link_libraries(bug_79_double_complex arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
+target_link_libraries(bug_79_double_complex arpack LAPACK::LAPACK BLAS::BLAS ${EXTRA_LDFLAGS})
 add_test(bug_79_double_complex_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bug_79_double_complex)

 add_executable(bug_142 TESTS/bug_142.f)
-target_link_libraries(bug_142 arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
+target_link_libraries(bug_142 arpack LAPACK::LAPACK BLAS::BLAS ${EXTRA_LDFLAGS})
 add_test(bug_142_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bug_142)

 add_executable(bug_142_gen TESTS/bug_142_gen.f)
-target_link_libraries(bug_142_gen arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
+target_link_libraries(bug_142_gen arpack LAPACK::LAPACK BLAS::BLAS ${EXTRA_LDFLAGS})
 add_test(bug_142_gen_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bug_142_gen)

 if(MPI)
   set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/PARPACK/TESTS/MPI)

   add_executable(issue46 PARPACK/TESTS/MPI/issue46.f)
-  target_link_libraries(issue46 parpack arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
+  target_link_libraries(issue46 parpack arpack LAPACK::LAPACK BLAS::BLAS ${EXTRA_LDFLAGS})
   add_test(issue46_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/issue46)
 endif()

@@ -578,12 +578,12 @@

   add_executable(icb_arpack_c TESTS/icb_arpack_c.c)
   target_include_directories(icb_arpack_c PUBLIC ${PROJECT_SOURCE_DIR}/ICB) # Get arpack.h
-  target_link_libraries(icb_arpack_c arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
+  target_link_libraries(icb_arpack_c arpack LAPACK::LAPACK BLAS::BLAS ${EXTRA_LDFLAGS})
   add_test(icb_arpack_c_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/icb_arpack_c)

   add_executable(icb_arpack_cpp TESTS/icb_arpack_cpp.cpp)
   target_include_directories(icb_arpack_cpp PUBLIC ${PROJECT_SOURCE_DIR}/ICB) # Get arpack.hpp
-  target_link_libraries(icb_arpack_cpp arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
+  target_link_libraries(icb_arpack_cpp arpack LAPACK::LAPACK BLAS::BLAS ${EXTRA_LDFLAGS})
   add_test(icb_arpack_cpp_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/icb_arpack_cpp)

   if (ICBEXMM)
@@ -591,7 +591,7 @@

     add_executable(arpackmm EXAMPLES/MATRIX_MARKET/arpackmm.cpp)
     target_include_directories(arpackmm PUBLIC ${PROJECT_SOURCE_DIR}/ICB ${EIGEN3_INCLUDE_DIR}) # Get arpack.h + eigen
-    target_link_libraries(arpackmm arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
+    target_link_libraries(arpackmm arpack LAPACK::LAPACK BLAS::BLAS ${EXTRA_LDFLAGS})
     configure_file(EXAMPLES/MATRIX_MARKET/As.mtx ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/As.mtx)
     configure_file(EXAMPLES/MATRIX_MARKET/An.mtx ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/An.mtx)
     configure_file(EXAMPLES/MATRIX_MARKET/Az.mtx ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Az.mtx)
@@ -606,7 +606,7 @@
     target_compile_definitions(pyarpack PRIVATE PY_MAJOR_VERSION="3")
     set(pyarpack_HDR ${PROJECT_SOURCE_DIR}/ICB ${PROJECT_SOURCE_DIR}/EXAMPLES/MATRIX_MARKET ${PROJECT_SOURCE_DIR}/EXAMPLES/PYARPACK)
     target_include_directories(pyarpack PUBLIC ${pyarpack_HDR} ${EIGEN3_INCLUDE_DIR} ${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS})
-    target_link_libraries(pyarpack BLAS::BLAS LAPACK::LAPACK ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
+    target_link_libraries(pyarpack LAPACK::LAPACK BLAS::BLAS ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
     install(TARGETS pyarpack
             ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/pyarpack
             LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/pyarpack)
@@ -662,12 +662,12 @@

     add_executable(icb_parpack_c PARPACK/TESTS/MPI/icb_parpack_c.c)
     target_include_directories(icb_parpack_c PUBLIC ${PROJECT_SOURCE_DIR}/ICB MPI::MPI_C) # Get parpack.h mpi.h
-    target_link_libraries(icb_parpack_c parpack arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS} MPI::MPI_C)
+    target_link_libraries(icb_parpack_c parpack arpack LAPACK::LAPACK BLAS::BLAS ${EXTRA_LDFLAGS} MPI::MPI_C)
     add_test(icb_parpack_c_tst mpirun -n 2 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/icb_parpack_c)

     add_executable(icb_parpack_cpp PARPACK/TESTS/MPI/icb_parpack_cpp.cpp)
     target_include_directories(icb_parpack_cpp PUBLIC ${PROJECT_SOURCE_DIR}/ICB MPI::MPI_CXX) # Get parpack.hpp mpi.h
-    target_link_libraries(icb_parpack_cpp parpack arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS} MPI::MPI_CXX)
+    target_link_libraries(icb_parpack_cpp parpack arpack LAPACK::LAPACK BLAS::BLAS ${EXTRA_LDFLAGS} MPI::MPI_CXX)
     add_test(icb_parpack_cpp_tst mpirun -n 2 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/icb_parpack_cpp)
   endif()
 endif()
EOF
    patch CMakeLists.txt arpack-cmake.patch
    mkdir build
    cd build
    CC=$MPI_C_COMPILER CXX=$MPI_CXX_COMPILER FC=$MPI_F90_COMPILER cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DBLAS_LIBRARIES=$NON_INTEL_BLAS_LIBRARY -DLAPACK_LIBRARIES="$NON_INTEL_LAPACK_LIBRARY" $GOMA_LIB/arpack-ng-$ARPACK_NG_VERSION/src -DCMAKE_INSTALL_PREFIX=$GOMA_LIB/arpack-ng-$ARPACK_NG_VERSION -DBUILD_SHARED_LIBS=OFF 2>&1 | tee -a $COMPILE_LOG
    make -j$MAKE_JOBS 2>&1 | tee -a $COMPILE_LOG
    make install 2>&1 | tee -a $COMPILE_LOG
    cd $GOMA_LIB
    if [ -e $GOMA_LIB/arpack-ng-$ARPACK_NG_VERSION/lib/libarpack.a ]
    then
        log_echo "Built Arpack-NG"
    else
        log_echo "Failed to build arpack-ng"
        exit 1
    fi
fi
 export LD_LIBRARY_PATH="arpack-ng-${ARPACK_NG_VERSION}:$LD_LIBRARY_PATH"

#make SuperLU
cd $GOMA_LIB/superlu_dist-$SUPERLU_DIST_VERSION
if [ -e lib64/libsuperlu_dist.a ] || [ -e lib/libsuperlu_dist.a ]
then
    log_echo "SuperLU_DIST already built"
else
    mkdir -p build
    CC=$MPI_C_COMPILER CXX=$MPI_CXX_COMPILER FC=$MPI_F90_COMPILER cmake -B build -DTPL_BLAS_LIBRARIES="${NON_INTEL_BLAS_LINK}" \
      -DTPL_PARMETIS_LIBRARIES="$GOMA_LIB/parmetis-$PARMETIS_VERSION/lib/libparmetis.a;$GOMA_LIB/metis-$METIS_VERSION/lib/libmetis.a" \
      -DCMAKE_INSTALL_PREFIX=$GOMA_LIB/superlu_dist-$SUPERLU_DIST_VERSION \
      -Denable_openmp:BOOL=FALSE \
      -DTPL_PARMETIS_INCLUDE_DIRS="$GOMA_LIB/parmetis-$PARMETIS_VERSION/include;$GOMA_LIB/metis-$METIS_VERSION/include" 2>&1 | tee -a $COMPILE_LOG
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
	     MY_METIS_LIB="-L$GOMA_LIB/metis-$METIS_VERSION/lib -lmetis" \
	     MY_METIS_INC="$GOMA_LIB/metis-$METIS_VERSION/include" \
	     AR="${ARCHIVER}" -C SuiteSparse_config library config 2>&1 | tee -a $COMPILE_LOG
	for target in "${targets[@]}"; do
	    make CUDA="no" AUTOCC="no" CC="${MPI_C_COMPILER}" \
                 CXX="${MPI_CXX_COMPILER}" \
                 BLAS="$NON_INTEL_BLAS_LINK" \
                 LAPACK="$SUITESPARSE_NON_INTEL_LAPACK_LINK" \
	         MY_METIS_LIB="-L$GOMA_LIB/metis-$METIS_VERSION/lib -lmetis" \
	         MY_METIS_INC="$GOMA_LIB/metis-$METIS_VERSION/include" \
	         AR="${ARCHIVER}" -C $target library 2>&1 | tee -a $COMPILE_LOG
	    make CUDA="no" AUTOCC="no" CC="${MPI_C_COMPILER}" \
                 CXX="${MPI_CXX_COMPILER}" \
                 BLAS="$NON_INTEL_BLAS_LINK" \
                 LAPACK="$SUITESPARSE_NON_INTEL_LAPACK_LINK" \
	         MY_METIS_LIB="-L$GOMA_LIB/metis-$METIS_VERSION/lib -lmetis" \
	         MY_METIS_INC="$GOMA_LIB/metis-$METIS_VERSION/include" \
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
        mkdir -p $GOMA_LIB/scalapack-$SCALAPACK_VERSION/lib
        mkdir -p $GOMA_LIB/scalapack-$SCALAPACK_VERSION/include
        CC=$MPI_C_COMPILER CXX=$MPI_CXX_COMPILER FC=$MPI_F90_COMPILER cmake -DCMAKE_BUILD_TYPE=Release \
         -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DBLAS_LIBRARIES=$NON_INTEL_BLAS_LIBRARY \
         -DLAPACK_LIBRARIES="$NON_INTEL_LAPACK_LIBRARY" -DCMAKE_INSTALL_PREFIX=$GOMA_LIB/scalapack-$SCALAPACK_VERSION \
         -B build 2>&1 | tee -a $COMPILE_LOG
        cmake --build build -j $MAKE_JOBS 2>&1 | tee -a $COMPILE_LOG
        cmake --install build 2>&1 | tee -a $COMPILE_LOG
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

LMETISDIR = $GOMA_LIB/parmetis-$PARMETIS_VERSION
IMETIS    = -I\$(LMETISDIR)/include -I$GOMA_LIB/metis-$METIS_VERSION/include
LMETIS    = -L\$(LMETISDIR)/lib -lparmetis -L$GOMA_LIB/metis-$METIS_VERSION/lib -lmetis

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
-D Trilinos_ENABLE_Epetra:BOOL=ON \
-D Trilinos_ENABLE_Xpetra:BOOL=ON \
-D Trilinos_ENABLE_Ifpack:BOOL=ON \
-D Trilinos_ENABLE_Ifpack2:BOOL=ON \
-D Trilinos_ENABLE_Teuchos:BOOL=ON \
-D Trilinos_ENABLE_ML:BOOL=ON \
-D Trilinos_ENABLE_MueLu:BOOL=ON \
-D Trilinos_ENABLE_AztecOO:BOOL=ON \
-D Trilinos_ENABLE_Belos:BOOL=ON \
-D Trilinos_ENABLE_Amesos2:BOOL=ON \
-D Trilinos_ENABLE_Sacado:BOOL=ON \
-D Trilinos_ENABLE_Tpetra:BOOL=ON \
-D Trilinos_ENABLE_Stratimikos:BOOL=ON \
-D Trilinos_ENABLE_Teko:BOOL=ON \
-D Trilinos_ENABLE_Belos:BOOL=ON \
-D Trilinos_ENABLE_Amesos2:BOOL=ON \
-D Trilinos_ENABLE_Sacado:BOOL=ON \
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
      -D Netcdf_LIBRARY_DIRS:PATH="$GOMA_LIB/netcdf-${NETCDF_VERSION}/lib;$GOMA_LIB/pnetcdf-${PNETCDF_VERSION}/lib" \
      -D TPL_ENABLE_Netcdf:BOOL=ON \
      -D TPL_Netcdf_INCLUDE_DIRS:PATH="$GOMA_LIB/netcdf-${NETCDF_VERSION}/include;$GOMA_LIB/pnetcdf-${PNETCDF_VERSION}/include" \
      -D Netcdf_LIBRARY_NAMES:STRING="netcdf;pnetcdf" \
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
  -D ParMETIS_LIBRARY_DIRS:PATH="$GOMA_LIB/parmetis-$PARMETIS_VERSION/lib;$GOMA_LIB/metis-$METIS_VERSION/lib" \
  -D TPL_ParMETIS_INCLUDE_DIRS:PATH="$GOMA_LIB/parmetis-$PARMETIS_VERSION/include;$GOMA_LIB/metis-$METIS_VERSION/include" \
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
-D Tpetra_INST_INT_INT:BOOL=ON \
$EXTRA_ARGS \
$GOMA_LIB/Trilinos-trilinos-release-$TRILINOS_VERSION_DASH 2>&1 | tee -a $COMPILE_LOG
# can set LONG LONG but only one is valid
# -D Tpetra_INST_INT_LONG_LONG:BOOL=ON \

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
    ./configure --with-shared-libraries=0 --with-cc=$(which mpicc) --with-cxx=$(which mpicxx) --with-fc=$(which mpif90) --with-debugging=0 COPTFLAGS='-O3' CXXOPTFLAGS='-O3' FOPTFLAGS='-O3' --download-strumpack --download-hypre --with-scalapack=1 --with-scalapack-dir=$(readlink --canonicalize-missing "${SCALAPACK_LIBRARY_DIR}/..") --with-superlu_dist=1 --with-superlu_dist-dir=$GOMA_LIB/superlu_dist-$SUPERLU_DIST_VERSION --with-metis=1 --with-metis-dir=$GOMA_LIB/metis-$METIS_VERSION --with-parmetis=1 --with-parmetis-dir=$GOMA_LIB/parmetis-$PARMETIS_VERSION --with-blas-lib=${NON_INTEL_BLAS_LIBRARY} --with-lapack-lib=${NON_INTEL_LAPACK_LIBRARY} --with-mumps=1 --with-mumps-dir="$GOMA_LIB/MUMPS_$MUMPS_VERSION"  2>&1 | tee -a $COMPILE_LOG
    make -j$MAKE_JOBS all 2>&1 | tee -a $COMPILE_LOG
    make check 2>&1 | tee -a $COMPILE_LOG
    if [ -e $PETSC_DIR/$PETSC_ARCH/lib/libpetsc.a ]; then
        log_echo "PETSc built!"
    else
        log_echo "Failed to build PETSc $PETSC_VERSION"
        exit 1
    fi
fi

if [ -n "$BUILD_PETSC_COMPLEX" ]; then

export PETSC_ARCH=arch-linux-c-debug-complex
if [ -e $PETSC_DIR/$PETSC_ARCH/lib/libpetsc.a ]; then
    log_echo "PETSc is already built!"
else
    ./configure --with-shared-libraries=0 --with-cc=$(which mpicc) --with-cxx=$(which mpicxx) --with-fc=$(which mpif90) --with-debugging=1 --download-hypre --with-scalapack=1 --with-scalapack-dir=$(readlink --canonicalize-missing "${SCALAPACK_LIBRARY_DIR}/..") --with-superlu_dist=1 --with-superlu_dist-dir=$GOMA_LIB/superlu_dist-$SUPERLU_DIST_VERSION --with-metis=1 --with-metis-dir=$GOMA_LIB/parmetis-4.0.3 --with-parmetis=1 --with-parmetis-dir=$GOMA_LIB/parmetis-4.0.3 --with-blas-lib=${NON_INTEL_BLAS_LIBRARY} --with-lapack-lib=${NON_INTEL_LAPACK_LIBRARY} --with-mumps=1 --with-mumps-dir="$GOMA_LIB/MUMPS_$MUMPS_VERSION" --with-scalar-type=complex 2>&1 | tee -a $COMPILE_LOG
    make -j$MAKE_JOBS all 2>&1 | tee -a $COMPILE_LOG
    make check 2>&1 | tee -a $COMPILE_LOG
    if [ -e $PETSC_DIR/$PETSC_ARCH/lib/libpetsc.a ]; then
        log_echo "PETSc built!"
    else
        log_echo "Failed to build PETSc $PETSC_VERSION"
        exit 1
    fi
fi

export PETSC_ARCH=arch-linux-c-opt-complex
if [ -e $PETSC_DIR/$PETSC_ARCH/lib/libpetsc.a ]; then
    log_echo "PETSc is already built!"
else
    ./configure --with-shared-libraries=0 --with-cc=$(which mpicc) --with-cxx=$(which mpicxx) --with-fc=$(which mpif90) --with-debugging=0 COPTFLAGS='-O3' CXXOPTFLAGS='-O3' FOPTFLAGS='-O3' --download-hypre --with-scalapack=1 --with-scalapack-dir=$(readlink --canonicalize-missing "${SCALAPACK_LIBRARY_DIR}/..") --with-superlu_dist=1 --with-superlu_dist-dir=$GOMA_LIB/superlu_dist-$SUPERLU_DIST_VERSION --with-metis=1 --with-metis-dir=$GOMA_LIB/parmetis-4.0.3 --with-parmetis=1 --with-parmetis-dir=$GOMA_LIB/parmetis-4.0.3 --with-blas-lib=${NON_INTEL_BLAS_LIBRARY} --with-lapack-lib=${NON_INTEL_LAPACK_LIBRARY} --with-mumps=1 --with-mumps-dir="$GOMA_LIB/MUMPS_$MUMPS_VERSION" --with-scalar-type=complex 2>&1 | tee -a $COMPILE_LOG
    make -j$MAKE_JOBS all 2>&1 | tee -a $COMPILE_LOG
    make check 2>&1 | tee -a $COMPILE_LOG
    if [ -e $PETSC_DIR/$PETSC_ARCH/lib/libpetsc.a ]; then
        log_echo "PETSc built!"
    else
        log_echo "Failed to build PETSc $PETSC_VERSION"
        exit 1
    fi
fi

fi

# Generate a config file for bash
cd $GOMA_LIB
export PETSC_ARCH=arch-linux-c-opt
cat > config.sh <<EOF
export CMAKE_PREFIX_PATH=$TRILINOS_INSTALL:${GOMA_LIB}/Omega_h:\$CMAKE_PREFIX_PATH
export PATH=$TRILINOS_INSTALL/bin:\$PATH
export PATH=${GOMA_LIB}/openmpi-$OPENMPI_VERSION/bin:\$PATH
export METISDIR=${GOMA_LIB}/metis-$METIS_VERSION
export UMFPACK_DIR=${GOMA_LIB}/SuiteSparse-$SUITESPARSE_VERSION
export ARPACKDIR=${GOMA_LIB}/arpack-ng-$ARPACK_NG_VERSION
export SPARSEDIR=${GOMA_LIB}/sparse
export PETSC_DIR=${PETSC_DIR}
export PETSC_ARCH=${PETSC_ARCH}
EOF

cat > config.fish <<EOF
set -x CMAKE_PREFIX_PATH $TRILINOS_INSTALL ${GOMA_LIB}/Omega_h \$CMAKE_PREFIX_PATH
set -x PATH $TRILINOS_INSTALL/bin \$PATH
set -x PATH ${GOMA_LIB}/openmpi-$OPENMPI_VERSION/bin \$PATH
set -x METISDIR ${GOMA_LIB}/metis-$METIS_VERSION
set -x UMFPACK_DIR ${GOMA_LIB}/SuiteSparse-$SUITESPARSE_VERSION
set -x ARPACKDIR ${GOMA_LIB}/arpack-ng-$ARPACK_NG_VERSION
set -x SPARSEDIR ${GOMA_LIB}/sparse
set -x PETSC_DIR ${PETSC_DIR}
set -x PETSC_ARCH ${PETSC_ARCH}
EOF

if [ -n "$BUILD_PETSC_COMPLEX" ]; then
cd $GOMA_LIB
export PETSC_ARCH=arch-linux-c-opt-complex
cat > config-complex.sh <<EOF
export CMAKE_PREFIX_PATH=$TRILINOS_INSTALL:${GOMA_LIB}/Omega_h:\$CMAKE_PREFIX_PATH
export PATH=$TRILINOS_INSTALL/bin:\$PATH
export PATH=${GOMA_LIB}/openmpi-$OPENMPI_VERSION/bin:\$PATH
export METISDIR=${GOMA_LIB}/metis-$METIS_VERSION
export UMFPACK_DIR=${GOMA_LIB}/SuiteSparse-$SUITESPARSE_VERSION
export ARPACKDIR=${GOMA_LIB}/arpack-ng-$ARPACK_NG_VERSION
export SPARSEDIR=${GOMA_LIB}/sparse
export PETSC_DIR=${PETSC_DIR}
export PETSC_ARCH=${PETSC_ARCH}
EOF

cat > config-complex.fish <<EOF
set -x CMAKE_PREFIX_PATH $TRILINOS_INSTALL ${GOMA_LIB}/Omega_h \$CMAKE_PREFIX_PATH
set -x PATH $TRILINOS_INSTALL/bin \$PATH
set -x PATH ${GOMA_LIB}/openmpi-$OPENMPI_VERSION/bin \$PATH
set -x METISDIR ${GOMA_LIB}/metis-$METIS_VERSION
set -x UMFPACK_DIR ${GOMA_LIB}/SuiteSparse-$SUITESPARSE_VERSION
set -x ARPACKDIR ${GOMA_LIB}/arpack-ng-$ARPACK_NG_VERSION
set -x SPARSEDIR ${GOMA_LIB}/sparse
set -x PETSC_DIR ${PETSC_DIR}
set -x PETSC_ARCH ${PETSC_ARCH}
EOF

cd $GOMA_LIB
export PETSC_ARCH=arch-linux-c-opt-debug-complex
cat > config-debug-complex.sh <<EOF
export CMAKE_PREFIX_PATH=$TRILINOS_INSTALL:${GOMA_LIB}/Omega_h:\$CMAKE_PREFIX_PATH
export PATH=$TRILINOS_INSTALL/bin:\$PATH
export PATH=${GOMA_LIB}/openmpi-$OPENMPI_VERSION/bin:\$PATH
export METISDIR=${GOMA_LIB}/metis-$METIS_VERSION
export UMFPACK_DIR=${GOMA_LIB}/SuiteSparse-$SUITESPARSE_VERSION
export ARPACKDIR=${GOMA_LIB}/arpack-ng-$ARPACK_NG_VERSION
export SPARSEDIR=${GOMA_LIB}/sparse
export PETSC_DIR=${PETSC_DIR}
export PETSC_ARCH=${PETSC_ARCH}
EOF

cat > config-debug-complex.fish <<EOF
set -x CMAKE_PREFIX_PATH $TRILINOS_INSTALL ${GOMA_LIB}/Omega_h \$CMAKE_PREFIX_PATH
set -x PATH $TRILINOS_INSTALL/bin \$PATH
set -x PATH ${GOMA_LIB}/openmpi-$OPENMPI_VERSION/bin \$PATH
set -x METISDIR ${GOMA_LIB}/metis-$METIS_VERSION
set -x UMFPACK_DIR ${GOMA_LIB}/SuiteSparse-$SUITESPARSE_VERSION
set -x ARPACKDIR ${GOMA_LIB}/arpack-ng-$ARPACK_NG_VERSION
set -x SPARSEDIR ${GOMA_LIB}/sparse
set -x PETSC_DIR ${PETSC_DIR}
set -x PETSC_ARCH ${PETSC_ARCH}
EOF
fi


if [ "$build_cmake" == "true" ] ; then
    echo "export PATH=$GOMA_LIB/cmake-$CMAKE_VERSION-linux-x86_64/bin:\$PATH" >> config.sh
    echo "set -x PATH $GOMA_LIB/cmake-$CMAKE_VERSION-linux-x86_64/bin \$PATH" >> config.fish
fi

log_echo
log_echo "An example bash configuration file has been written to $GOMA_LIB/config.sh"
log_echo
log_echo "Activate with $ source $GOMA_LIB/config.sh"
if [ -n "$BUILD_PETSC_COMPLEX" ]; then
log_echo "An example bash configuration file has been written to $GOMA_LIB/config-complex.sh"
log_echo
log_echo "Activate with $ source $GOMA_LIB/config-complex.sh"
fi
