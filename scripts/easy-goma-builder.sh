#!/usr/bin/env bash
#
# This software is distributed under the GNU General Public License.
#
# Goma dependency build script builder
#

# This script tries to build all of the libraries needed for goma
# in the specified directories.
#
# The libraries built should be suitable for a minimal settings.mk

# Set script to exit on error
# set -e

#Try to set soft ulimit to needed number. Need to open approximately 1200 files and counting to build Goma
if ulimit -S -n 2048 &> /dev/null; then
    echo
else
    echo "Compiling Goma may require opening up to 2048 files."
    echo "But your security settings prevent me from opening more than:"
    echo "`ulimit -n` files"
    echo "Please contact your system administrator about fixing this."
    exit 1
fi

MPI_NAME="open"
if [ -d "/usr/lib64/openmpi" ] ; then
    MPI_BASE_DIR="/usr/lib64/openmpi"
elif [ -d "/usr/lib/openmpi" ] ; then
    MPI_BASE_DIR="/usr/lib/openmpi"
else
    MPI_BASE_DIR="BUILD"
fi
if which mpiicc &> /dev/null; then
    MPI_BASE_DIR="BUILD"
fi

CC_NAME="gnu"
CC_CMD="gcc"
COMPILER_VERSION=$(gcc -dumpversion)
[ "$COMPILER_VERSION" = "`echo -e "4.8.1\n$COMPILER_VERSION" | sort -V | tail -n1`" ] && CXX11="true" || CXX11="false"
if [[ "$CXX11" == "true" ]]; then
   CXX11="true"
else
   CXX11="false"
   echo "Cannot use easy-builder with a compiler that does not support c++11"
   exit 1
fi
MATH_LIBRARIES="netlib blas"
SYS_LIB="-lm -lz"
User_Flags="-O1"


cd $(dirname `readlink --canonicalize $0`)
source user-interaction.sh
GOMA_LIB=gomaTPLs
if [ "$MPI_BASE_DIR" == "BUILD" ]; then
    QUESTION="The recommended goma TPL settings for your system are: \n\
Compiler: $CC_NAME $COMPILER_VERSION    `which $CC_CMD`\n\
MPI: ${MPI_NAME}MPI   (to be built) \n\
C compiler: Built mpicc \n\
C++ compiler: Built mpicxx \n\
F compiler: Built mpifort \n\
Math Libraries: $MATH_LIBRARIES \n\
Library Location: `readlink --canonicalize ${GOMA_LIB}` \n\
\n\
Build Goma with these settings?"
elif [ "$MPI_NAME" == "intel" ]; then
    QUESTION="The recommended goma TPL settings for your system are: \n\
Compiler: $CC_NAME $COMPILER_VERSION    `which $CC_CMD`\n\
MPI: ${MPI_NAME}MPI   $MPI_BASE_DIR \n\
C compiler: `which mpicc` \n\
C++ compiler: `which mpicxx` \n\
F compiler: `which mpifc` \n\
Math Libraries: $MATH_LIBRARIES 3.7.1    (to be built) \n\
Library Location: `readlink --canonicalize ${GOMA_LIB}` \n\
\n\
Build Goma with these settings?"
else
    QUESTION="The recommended goma TPL settings for your system are: \n\
Compiler: $CC_NAME $COMPILER_VERSION    `which $CC_CMD`\n\
MPI: ${MPI_NAME}MPI   $MPI_BASE_DIR \n\
C compiler: `which mpicc` \n\
C++ compiler: `which mpicxx` \n\
F compiler: `which mpifort` \n\
Math Libraries: $MATH_LIBRARIES 3.7.1    (to be built) \n\
Library Location: `readlink --canonicalize ${GOMA_LIB}` \n\
\n\
Build Goma with these settings?"
fi
recommended="true"

ok_check
if [ "$isok" == "false" ]; then
    compiler_test Goma
    mpi_test Goma
    mkl_test Goma
    echo "Enter the build path to use"
    read GOMA_LIB
    if [ -z "$GOMA_LIB" ]; then
        GOMA_LIB="gomaTPLs"
    fi
    echo "Enter the number of processors you want to use to compile the libraries"
    read USED_MAKE_JOBS
    if [ -z "$USED_MAKE_JOBS" ]; then
        USED_MAKE_JOBS="1"
    fi
    confirm Easy Goma Builder script
else
    echo "Enter the number of processors you want to use to compile the libraries"
    read USED_MAKE_JOBS
    if [ -z "$USED_MAKE_JOBS" ]; then
        USED_MAKE_JOBS="1"
    fi
    QUESTION="Start goma TPL build?"
    recommended="true"
    ok_check
    if [ "$isok" == "false" ]; then
        thank_user Easy Goma Builder script
        exit 0
    fi
fi

export MPI_BASE_DIR
export CC_NAME
export MPI_NAME
export MATH_LIBRARIES
export MATH_PATH
export INTEL_PARALLEL_STUDIO_ROOT

if [ "$CXX11" = "true" ]; then
    echo "Starting TPL build..."
    ./build-goma-dep-trilinos-12.sh -j ${USED_MAKE_JOBS} "${GOMA_LIB}" 2>&1 | tee easygomabuild.txt
    # If user failed a continue_check exit program entirely
    if [ ${PIPESTATUS[0]} -eq 24 ]; then
        thank_user Easy Goma Builder script
        exit 0
    fi
    # In case cmake was built in the build-goma-deps script add it to the path. This won't do anything if cmake is already installed because it appends it to the end.
    export PATH="$PATH:$GOMA_LIB/cmake-2.8.12.2/bin"
else
    echo "Your compiler does not support c++11 and thus is not supported by this script at this time"
    thank_user Easy Goma Builder script
    exit 1
fi



if [ "$MPI_NAME" == "open" ]; then
    if [ "$MPI_BASE_DIR" == "BUILD" ]; then
        export PATH="`readlink --canonicalize ${GOMA_LIB}`/openmpi-2.1.1/bin:$PATH"
        export LD_LIBRARY_PATH="`readlink --canonicalize ${GOMA_LIB}`/openmpi-2.1.1/lib:$LD_LIBRARY_PATH"
        BUILD_MPI="true"
        # This lets you at least compile openMPI with intel compiler, but this still breaks some unit tests
        if [ "$CC_NAME" == "intel" ]; then
            MPI_LIB="-L`readlink --canonicalize ${GOMA_LIB}`/openmpi-2.1.1/lib -lmpi -lmpi_mpifh -lifcore"
        else
            MPI_LIB="-L`readlink --canonicalize ${GOMA_LIB}`/openmpi-2.1.1/lib -lmpi -lmpi_mpifh"
        fi
    else
        MPI_LIB="-L${MPI_BASE_DIR}/lib -lmpi -lmpi_mpifh"
    fi
    C_COMP=`which mpicc`
    CXX_COMP=`which mpicxx`
    FORT_COMP=`which mpifort`
else
    MPI_LIB="-L${I_MPI_ROOT}/intel64/lib -lmpi -lmpi_ilp64 -lmpifort -lmpigf -lifcore"
    if [ "$CC_NAME" == "intel" ]; then
        C_COMP=`which mpiicc`
        CXX_COMP=`which mpiicpc`
        FORT_COMP=`which mpiifort`
    else
        C_COMP=`which mpicc`
        CXX_COMP=`which mpicxx`
        if [ "$MPI_NAME" == "open" ]; then
            FORT_COMP=`which mpifort`
        else
            FORT_COMP=`which mpifc`
        fi
    fi
fi

if [ "$MATH_LIBRARIES" == "intel" ] && [ "$CC_NAME" == "gnu" ]; then
    SYS_LIB="${MKLROOT}/lib/intel64/libmkl_blas95_ilp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread"
    User_Flags="-O1 -I${MKLROOT}/include/intel64/ilp64 -m64 -I${MKLROOT}/include"
fi

export GOMA_LIBS=`readlink --canonicalize ${GOMA_LIB}`
rm -rf ../build
mkdir -p ../build
cd ../build

echo "Easy-Goma-Builder v 1.0: BUILD VARIABLES ARE AS FOLLOWS -DCMAKE_Fortran_COMPILER=$FORT_COMP -DCMAKE_C_COMPILER=$C_COMP -DCMAKE_CXX_COMPILER=$CXX_COMP -Dgoma_LIBS=${GOMA_LIBS} -Dgoma_MPI_DIR=${MPI_BASE_DIR} -Dgoma_MPI_LIB=$MPI_LIB -Dgoma_User_Flags=$User_Flags -Dgoma_SYS_LIB=$SYS_LIB -Dgoma_Trilinos_DIR=$GOMA_LIBS/trilinos-12.10.1-Built" | tee -a ../scripts/easygomabuild.txt

cmake .. -DCMAKE_Fortran_COMPILER="$FORT_COMP" -DCMAKE_C_COMPILER="$C_COMP" -DCMAKE_CXX_COMPILER="$CXX_COMP" -Dgoma_LIBS=${GOMA_LIBS} -Dgoma_MPI_DIR=${MPI_BASE_DIR} -Dgoma_MPI_LIB="$MPI_LIB" -Dgoma_User_Flags="$User_Flags" -Dgoma_SYS_LIB="$SYS_LIB" -Dgoma_Trilinos_DIR="$GOMA_LIBS/trilinos-12.10.1-Built" | tee -a ../scripts/easygomabuild.txt

if make -j ${USED_MAKE_JOBS} 2>&1 && [ -f goma ] | tee -a ../scripts/easygomabuild.txt; then
    mkdir -p ../bin
    mv goma ../bin
    echo "Goma has been built successfully!"
    echo "Visit https://goma.github.io for useful documentation & tutorials"
    echo
    echo "Goma is located at"
    echo $(readlink --canonicalize "`dirname $0`/../bin/goma")
    echo
    echo "It's probably necessary that you run:"
    echo "export PATH=\$PATH:$GOMA_LIBS/trilinos-12.10.1-Built/bin"
    echo "before running Goma"
    echo
    if [ "$BUILD_MPI" == "true" ]; then
        # It's important that openMPI comes first in the path or it will be overridden by intel
        echo "You have built MPI instead of using a native version"
        echo "You will need to run"
        echo "export LD_LIBRARY_PATH=$GOMA_LIBS/openmpi-2.1.1/lib:\$LD_LIBRARY_PATH"
        echo "export PATH=$GOMA_LIBS/openmpi-2.1.1/bin:\$PATH"
        echo
        echo "before running Goma."
    fi
    if [ "$MPI_NAME" == "intel" ]; then
        # Without this some unit tests hang FOREVER?!
        echo "You have used intel MPI when building Goma"
        echo "When running parallel problems, it is important that you first run"
        echo "export I_MPI_SHM_LMT=shm"
        echo
        echo "Or goma could possibly get stuck in an endless loop."
    fi
else
    echo "Something has gone wrong. Please create a GitHub issue."
    echo "Attach the logfile scripts/easygomabuild.txt,"
    echo "and provide any information which might be helpful to debugging."
fi
thank_user Easy Goma Builder script
