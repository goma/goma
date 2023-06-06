function thank_user {
    echo
    echo "Thank you for using the $@"
    echo "Have a nice day."
}

function continue_check {
    echo "Enter \"c\" to continue (any other letter to exit):"
    read -N 1 user_choice
    if [ "$user_choice" != 'c' ]
        then
        exit 24
    fi
}

function ok_check {
    echo
    if type dialog &> /dev/null; then
        if [ "$recommended" == "true" ]; then
            dialog --title "Confirm"  --backtitle "Press escape to quit" --yesno "${QUESTION}"  15 60
        else
            dialog --defaultno --title "Confirm" --backtitle "Press escape to quit" --yesno "${QUESTION}" 15 60
        fi
        choice=$?
        case "$choice" in
            0 ) isok="true" ;;
            1 ) isok="false" ;;
            255 ) thank_user Goma script suite && exit 0 ;;
            * ) exit 0 ;;
        esac
    else
        if [ "$recommended" == "true" ]; then
            printf "$QUESTION (Y/n) "
        else
            printf "$QUESTION (y/N) "
        fi
        read choice
        case "$choice" in
            y|Y ) isok="true" ;;
            n|N ) isok="false" ;;
            * ) isok="$recommended" &&
                if [ "$recommended" == "true" ]; then echo "y"
                else echo "n"
                fi ;;
        esac
    fi
}


function compiler_test {
    #First test for compilers
    #QUESTION="Do you want to use the Intel compiler? (not recommended)"
    #recommended="false"
    #ok_check
    isok="false"
    if [ "$isok" == "true" ]; then
        if icc --version &> /dev/null ; then
            echo "ICC detected. Installing with Intel compiler"
            CXX11="true"
            CC_NAME="intel"
            CC_CMD="icc"
            FORTRAN_LIB="-lifcore"
            echo "Please enter where you installed the Intel compiler: (EG /opt/intel)"
            read INTEL_PARALLEL_STUDIO_ROOT
        else
            echo "Unable to find an Intel compiler."
            echo "If you want to use it, please source compilervars.sh first"
            echo "Please see https://software.intel.com/en-us/articles/setting-up-the-build-environment-for-using-intel-c-or-fortran-compilers for more information"
            exit 1
        fi
    else
        if COMPILER_VERSION=$(gcc -dumpversion) ; then
            [ "$COMPILER_VERSION" = "`echo -e "4.8.1\n$COMPILER_VERSION" | sort -V | tail -n1`" ] && CXX11="true" || CXX11="false"
            CC_NAME="gnu"
            FORTRAN_LIB="-lc -lgcc -lstdc++ -lgfortran"
            if [[ "$CXX11" == "true" ]]; then
                echo "C++ 11 support detected on gcc compiler. Will compile $1 with it."
                CXX11="true"
            else
                echo "C++ 11 support not detected on compiler. Cannot complete $1 build."
                exit 1
            fi
        else
            echo "GCC not installed or in your path. At this time only GCC and the Intel compiler work"
            exit 1
        fi
    fi
}

function mpi_test {
    #Test for MPI
    if which mpiicc &> /dev/null ; then
        QUESTION="Intel MPI detected and on current PATH. Install using Intel MPI? This also requires using the intel compiler, and will override your current compiler setting. (not recommended) "
        recommended="false"
        ok_check
        if [ "$isok" == "true" ]; then
            MPI_NAME="intel"
            MPI_BASE_DIR="intelMPI"
            CC_NAME="intel"
        else
            MPI_NAME="open"
            MPI_BASE_DIR="BUILD"
            echo "Okay, openMPI will be downloaded and built instead."
        fi
    elif mpicc --version &> /dev/null && [ "$CC_NAME" == "gnu" ] ; then
        echo "Native non-Intel MPI detected."
        QUESTION="Install $1 using native MPI? (recommended)"
        recommended="true"
        ok_check
        if [ "$isok" == "true" ]; then
            MPI_NAME="open"
            if [ -d "/usr/lib64/openmpi" ] ; then
                MPI_BASE_DIR="/usr/lib64/openmpi"
            elif [ -d "/usr/lib/openmpi" ] ; then
                MPI_BASE_DIR="/usr/lib/openmpi"
            elif [ -d "/usr/lib/x86_64-linux-gnu/openmpi" ] ; then
                MPI_BASE_DIR="/usr/lib/x86_64-linux-gnu/openmpi"
            else
                echo "You have openMPI natively installed but I was unable to find it"
                echo "Please enter your openMPI path"
                echo "Or leave blank to download and build openMPI"
                MPI_BASE_DIR=read -r
                if [ -z "$MPI_BASE_DIR" ]; then
                    echo "Okay, openMPI will be downloaded and built instead."
                    MPI_BASE_DIR="BUILD"
                fi
            fi
        else
            echo "Okay, openMPI will be downloaded and built instead."
            MPI_BASE_DIR="BUILD"
        fi
    elif [ "$CC_NAME" == "intel" ]; then
        MPI_BASE_DIR="BUILD"
        echo "Building MPI from scratch to use with the intel compiler."
    else
        echo "I couldn't find MPI installed."
        QUESTION="I couldn't find MPI installed.\n\
Want me to build openMPI from scratch? (recommended)"
        recommended="true"
        ok_check
        if [ "$isok" == "true" ]; then
            MPI_BASE_DIR="BUILD"
        else
            echo "$1 require(s) MPI to function properly and it wasn't found."
            exit 1
        fi
    fi
}


function mkl_test {
    QUESTION="Do you want to use OpenBLAS with $1? (recommended)"
    recommended="true"
    ok_check
    if [ "$isok" == "true" ]; then
        echo "Using OpenBLAS for BLAS and LAPACK"
        MATH_LIBRARIES="openblas"
        SYS_LIB="-lm -lz"
    else
        echo "Blas will be downloaded and built from Netlib"
        MATH_LIBRARIES="netlib blas"
        SYS_LIB="-lm -lz"
        User_Flags="-O1"
    fi
    
}

function confirm {

    if [ "$CC_NAME" == "gnu" ]; then
        QUESTION="Compiler: $CC_NAME $COMPILER_VERSION    `which $CC_CMD`\n\
"
    else
        QUESTION="Compiler: $CC_NAME    `which $CC_CMD`\n\
"
    fi
    if [ "$MPI_BASE_DIR" == "BUILD" ]; then
        QUESTION="${QUESTION}MPI: ${MPI_NAME}MPI (to be built) \n\
C compiler: Built mpicc \n\
C++ compiler: Built mpicxx \n\
F compiler: Built mpifort \n\
"
    elif [ "$CC_NAME" == "gnu" ] && [ "$MPI_NAME" == "intel" ]; then
        QUESTION="${QUESTION}MPI: ${MPI_NAME}MPI (to be built) \n\
C compiler: `which mpicc` \n\
C++ compiler: `which mpicxx` \n\
F compiler: `which mpifc` \n\
"
    elif [ "$CC_NAME" == "intel" ] && [ "$MPI_NAME" == "intel" ]; then
        QUESTION="${QUESTION}MPI: ${MPI_NAME}MPI (to be built) \n\
C compiler: `which mpiicc` \n\
C++ compiler: `which mpiicpc` \n\
F compiler: `which mpiifort` \n\
"
    else
        QUESTION="${QUESTION}MPI: ${MPI_NAME}MPI   $MPI_BASE_DIR \n\
C compiler: `which mpicc` \n\
C++ compiler: `which mpicxx` \n\
F compiler: `which mpifort` \n\
"
    fi
    if [ "$MATH_LIBRARIES" == "netlib blas" ]; then
        QUESTION="${QUESTION}Math Libraries: $MATH_LIBRARIES 3.7.1 (to be built) \n\
"
    elif [ "$MATH_LIBRARIES" == "atlas" ]; then
        QUESTION="${QUESTION}Math Libraries: $MATH_LIBRARIES     $MATH_PATH \n\
"
    else
        QUESTION="${QUESTION}Math Libraries: $MATH_LIBRARIES     $I_MKL_ROOT \n\
"
    fi

    QUESTION="${QUESTION}Library Location: `readlink --canonicalize ${GOMA_LIB}` \n\
Compile Processorss: $USED_MAKE_JOBS \n\
\n\
Proceed?"
    recommended=true
    ok_check
    if [ "$isok" == "false" ]; then
        thank_user $@
        exit 0
    fi
}
