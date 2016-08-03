#!/usr/bin/env bash
#
# This software is distributed under the GNU General Public License.
#
# Goma dependency build script
#
# 2014, July: Created by Cory Miner based off of notes by Scott Roberts
# 2014: Weston Ortiz modified to include blacs, scalapack and mumps support in trilinos
# 2015, February: Andrew Cochrane modified to handle recent changes to SEACAS-latest and hdf5, also added some logic to help with troubleshooting

# This script tries to build all of the libraries needed for goma
# in the specified directories.
#
# The libraries built should be suitable for a minimal settings.mk

# Set script to exit on error
# set -e

function usage() {
    echo "Usage: build-goma-dependencies [options] [library install location]"
    echo "       Options:"
    echo "               -jN  N : Number of make jobs to run"
}

MAKE_JOBS=1

while getopts ":j:h" opt; do
    case $opt in
        j)
            MAKE_JOBS=$OPTARG
            ;;
        h)
            usage
            ;;
    esac
done

shift $((OPTIND - 1))

if [ -z "$1" ]
    then
    echo "ERROR: Missing library install location"
    usage
    exit 1
fi

cd $1
GOMA_LIB=`pwd`
export GOMA_LIB

OWNER=$USER

function continue_check {
    echo "Enter \"c\" to continue (any other letter to exit):"

    read -N 1 user_choice

    if [ $user_choice != 'c' ]
        then
        exit
    fi
}

echo "The current library configuration"
echo "Libraries will be installed at $GOMA_LIB."
echo "This script will use $MAKE_JOBS make jobs."
continue_check

FORTRAN_LIBS=-lgfortran

OPENMPI_TOP=$GOMA_LIB/openmpi-1.6.4
export PATH=$OPENMPI_TOP/bin:$PATH
export LD_LIBRARY_PATH=$OPENMPI_TOP/lib:$LD_LIBRARY_PATH
FORTRAN_COMPILER=$OPENMPI_TOP/bin/mpif90

export MPI_BASE_DIR=$OPENMPI_TOP
export GOMA_LIB
export OWNER
export OPENMPI_TOP
export MAKE_JOBS

ARCHIVE_NAMES=("arpack96.tar.gz" \
"patch.tar.gz" \
"blas-3.6.0.tgz" \
"cmake-2.8.12.2.tar.gz" \
"hdf5-1.8.15.tar.gz" \
"lapack-3.2.1.tgz" \
"netcdf-4.3.3.1.tar.gz" \
"openmpi-1.6.4.tar.gz" \
"parmetis-4.0.3.tar.gz" \
"sparse.tar.gz" 
"superlu_dist_5.0.0.tar.gz" \
"y12m-1.0.tar.gz" \
"trilinos-12.6.3-Source.tar.bz2" \
"scalapack-2.0.2.tgz" \
"MUMPS_5.0.1.tar.gz" \
"SuiteSparse-4.4.4.tar.gz" \
"matio-1.5.2.tar.gz")

ARCHIVE_MD5SUMS=("fffaa970198b285676f4156cebc8626e" \
"14830d758f195f272b8594a493501fa2" \
"0af668589d693beb91ad49ddd742a002" \
"17c6513483d23590cbce6957ec6d1e66" \
"03cccb5b33dbe975fdcd8ae9dc021f24" \
"a3202a4f9e2f15ffd05d15dab4ac7857" \
"5c9dad3705a3408d27f696e5b31fb88c" \
"70aa9b6271d904c6b337ca326e6613d1" \
"f69c479586bf6bb7aff6a9bc0c739628" \
"1566d914d1035ac17b73fe9bc0eed02a" \
"2b53baf1b0ddbd9fcf724992577f0670" \
"eed01310baca61f22fb8a88a837d2ae3" \
"d94e31193559b334fd41d05eb22f9285" \
"2f75e600a2ba155ed9ce974a1c4b536f" \
"b477573fdcc87babe861f62316833db0" \
"e0af74476935c9ff6d971df8bb6b82fc" \
"85b007b99916c63791f28398f6a4c6f1")

ARCHIVE_URLS=("http://www.caam.rice.edu/software/ARPACK/SRC/arpack96.tar.gz" \
"http://www.caam.rice.edu/software/ARPACK/SRC/patch.tar.gz" \
"http://www.netlib.org/blas/blas-3.6.0.tgz" \
"http://www.cmake.org/files/v2.8/cmake-2.8.12.2.tar.gz" \
"http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.15/src/hdf5-1.8.15.tar.gz" \
"http://www.netlib.org/lapack/lapack-3.2.1.tgz" \
"ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.3.3.1.tar.gz" \
"http://www.open-mpi.org/software/ompi/v1.6/downloads/openmpi-1.6.4.tar.gz" \
"http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz" \
"http://sourceforge.net/projects/sparse/files/sparse/sparse1.4b/sparse1.4b.tar.gz/download" \
"http://crd-legacy.lbl.gov/~xiaoye/SuperLU/superlu_dist_5.0.0.tar.gz" \
"http://sisyphus.ru/cgi-bin/srpm.pl/Branch5/y12m/getsource/0" \
"http://trilinos.csbsju.edu/download/files/trilinos-12.6.3-Source.tar.bz2" \
"http://www.netlib.org/scalapack/scalapack-2.0.2.tgz" \
"http://graal.ens-lyon.fr/MUMPS/MUMPS_5.0.1.tar.gz" \
"http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-4.4.4.tar.gz" \
"http://downloads.sourceforge.net/project/matio/matio/1.5.2/matio-1.5.2.tar.gz")

ARCHIVE_DIR_NAMES=("ARPACK" \
"ARPACK" \
"BLAS-3.6.0" \
"cmake-2.8.12.2" \
"hdf5-1.8.15" \
"lapack-3.2.1" \
"netcdf-4.3.3.1" \
"openmpi_1.6.4" \
"parmetis-4.0.3" \
"sparse" \
"SuperLU_DIST_5.0.0" \
"y12m-1.0" \
"trilinos-12.6.3-Source" \
"scalapack-2.0.2" \
"MUMPS_5.0.1" \
"SuiteSparse" \
"matio-1.5.2")

SEACAS_LINUX_PATCH="14c14
< /* #define GCC4GFORTRAN 1 */
---
> #define GCC4GFORTRAN 1 
266c266
< #define X11Includes -I/usr/include/X11R6
---
> #define X11Includes -I/usr/include/X11
"

SITE_PATCH="4c4
< #define Owner gdsjaar
---
> #define Owner $OWNER
11c11
< #define AccessRoot /Users/gdsjaar/src/SEACAS-SF
---
> #define AccessRoot $GOMA_LIB/SEACAS-2013-12-03
47a48,50
> #define ExcludeAnalysis
> #define Parallel 0
> 
48a52
> #define HasMatlab No
"

read -d '' SCALAPACK_PATCH << "EOF"
21a22
> GOMA_LIB     = __GOMA_LIB__
58,59c59,60
< BLASLIB       = -lblas
< LAPACKLIB     = -llapack
---
> BLASLIB       = -L$(GOMA_LIB)/BLAS-3.6.0 -lblas
> LAPACKLIB     = -L$(GOMA_LIB)/lapack-3.2.1 -llapack
EOF
SCALAPACK_PATCH=${SCALAPACK_PATCH/__GOMA_LIB__/$GOMA_LIB}

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
> CC = mpicc
> FC = mpif77
> FL = mpif77
72,74c75,76
< #RANLIB = ranlib
< RANLIB  = echo
< SCALAP  = /local/SCALAPACK/libscalapack.a /local/BLACS/LIB/blacs_MPI-LINUX-0.a   /local/BLACS/LIB/blacsF77init_MPI-LINUX-0.a  /local/BLACS/LIB/blacs_MPI-LINUX-0.a
---
> RANLIB = ranlib
> SCALAP  = $(GOMA_LIB)/scalapack-2.0.2/libscalapack.a 
76,78c78,80
< INCPAR = -I/usr/local/mpich/include
< # LIBPAR = $(SCALAP)  -L/usr/local/lib/ -llammpio -llamf77mpi -lmpi -llam -lutil -ldl -lpthread
< LIBPAR = $(SCALAP)  -L/usr/local/mpich/lib/ -lmpich 
---
> MPI_HOME=$(GOMA_LIB)/openmpi-1.6.4
> INCPAR = -I$(MPI_HOME)
> LIBPAR = $(SCALAP)  -L$(MPI_HOME)/lib/ -lmpi -lmpi_f77 -lutil -ldl -lpthread
83c85
< LIBBLAS = -L/local/BLAS -lblas
---
> LIBBLAS = -L$(GOMA_LIB)/BLAS-3.6.0 -lblas
EOF
MUMPS_PATCH=${MUMPS_PATCH/__GOMA_LIB__/$GOMA_LIB}

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

read -d '' ARMAKE_PATCH << "EOF"
28c28
< home = $(HOME)/ARPACK
---
> home = __GOMA_LIB5__/ARPACK
35c35
< PLAT = SUN4
---
> PLAT = x86_64
104,105c104,105
< FC      = f77
< FFLAGS	= -O -cg89
---
> FC      = gfortran
> FFLAGS	= -O
115c115
< MAKE    = /bin/make
---
> MAKE    = /usr/bin/make
EOF
ARMAKE_PATCH=${ARMAKE_PATCH/__GOMA_LIB5__/$GOMA_LIB}

# read -d '' UFCONFIG_PATCH << "EOF"
# 61c61,64
# < ARCHIVE = $(AR) $(ARFLAGS)
# ---
# > AR = ar cr
# > ARFLAGS = cr
# > ARCHIVE = ar cr
# > 
# 111,113c114,115
# < BLAS = -lblas -lgfortran
# < LAPACK = -llapack
# < 
# ---
# > BLAS = __GOMA_LIB4__/BLAS/blas_LINUX.a -lgfortran
# > LAPACK = __GOMA_LIB4__/lapack-3.2.1/lapack_LINUX.a
# 173c175
# < UMFPACK_CONFIG =
# ---
# > UMFPACK_CONFIG = -DNCHOLMOD
# 214c216
# < CHOLMOD_CONFIG =
# ---
# > CHOLMOD_CONFIG = 

# EOF
# UFCONFIG_PATCH=${UFCONFIG_PATCH/__GOMA_LIB4__/$GOMA_LIB}

read -d '' SUITESPARSE_CONFIG_PATCH << "EOF"
100a101,102
> GOMA_LIB = __GOMA_LIB__
> 
123c125
<   LAPACK = -llapack
---
> #  LAPACK = -llapack
134c136
<   BLAS = -lopenblas
---
> #  BLAS = -lopenblas
149a152,154
> BLAS = -L$(GOMA_LIB)/BLAS-3.6.0 -lblas -lgfortran
> LAPACK = -L$(GOMA_LIB)/lapack-3.2.1 -llapack -lgfortran
> 
234c239
< UMFPACK_CONFIG =
---
> UMFPACK_CONFIG = -DNCHOLMOD

EOF
SUITESPARSE_CONFIG_PATCH=${SUITESPARSE_CONFIG_PATCH/__GOMA_LIB__/$GOMA_LIB}


read -d '' LAPACK_PATCH << "EOF"
54c54
< BLASLIB      = ../../blas$(PLAT).a
---
> BLASLIB      = __GOMA_LIB3__/BLAS-3.6.0/blas_LINUX.a

EOF
LAPACK_PATCH=${LAPACK_PATCH/__GOMA_LIB3__/$GOMA_LIB}


read -d '' SUPERLU_PATCH << "EOF"
19c19
< PLAT      = _power5
---
> PLAT      =
24c24
< DSuperLUroot  = $(HOME)/Release_Codes/SuperLU_DIST_2.3
---
> DSuperLUroot  = __GOMA_LIB2__/SuperLU_DIST_2.3
28c28
< BLASLIB       = -lessl
---
> BLASLIB       = __GOMA_LIB2__/BLAS-3.6.0 -lessl
31,32c31,32
< METISLIB        = -L/project/projectdirs/sparse/xiaoye/parmetis-3.1/64 -lmetis
< PARMETISLIB   = -L/project/projectdirs/sparse/xiaoye/parmetis-3.1/64 -lparmetis
---
> METISLIB        = __GOMA_LIB2__/ParMetis-3.1.1 -lmetis
> PARMETISLIB   = __GOMA_LIB2__/ParMetis-3.1.1 -lparmetis
42,43c42,43
< ARCHFLAGS     = -X64 cr
< RANLIB        = ranlib
---
> ARCHFLAGS     = -rc
> RANLIB        = echo
48c48
< CC            = mpcc_r
---
> CC            = mpicc
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
> FORTRAN         = mpif90
> FFLAGS        = -O3 -Q 
65c64
< LOADER    = mpxlf90_r
---
> LOADER    = mpicc
68c67
< LOADOPTS  = -q64
---
> LOADOPTS  = 
74c73
< CDEFS        = -DNoChange
---
> CDEFS        = -DAdd_

EOF
SUPERLU_PATCH=${SUPERLU_PATCH/__GOMA_LIB2__/$GOMA_LIB}

mkdir -p $GOMA_LIB
cd $GOMA_LIB

function mychecksum {
    local count=$2
    local archive=$1
    local MD5SAVED=${ARCHIVE_MD5SUMS[count]}
    local MD5ARCHIVE=($(md5sum $archive))
    if [ $MD5SAVED != $MD5ARCHIVE ]; then
        echo "Issue checksum with archive:"
	echo "$MD5SAVED"
	echo "$MD5ARCHIVE"
        echo $archive
        continue_check
    fi
}

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
#move

# #Seacas
# if [ -d SEACAS-2013-12-03 ]
# then
#     echo "SEACAS directory already renamed"
# else
# mv 2013-12-03 SEACAS-2013-12-03
# fi

# if [ -d SEACAS-2013-12-03/TPL/hdf5/hdf5-source ]
# then
#     echo "hdf5 directory already placed in SEACAS/TPL"
# else
# mv hdf5-1.8.12 hdf5-source
# mv hdf5-source SEACAS-2013-12-03/TPL/hdf5/
# fi

# if [ -d SEACAS-2013-12-03/TPL/netcdf/netcdf-4.3.2 ]
# then
#     echo "netcdf directory already placed in SEACAS/TPL"
# else
# mv netcdf-4.3.2 SEACAS-2013-12-03/TPL/netcdf/
# fi



#continue_check
#UMFPACK
# if [ -d UMFPACK-5.4/AMD ]
# then
#     echo "UMFPACK directory structure already created"
# else
# mkdir UMFPACK-5.4
# mv UMFPACK UMFPACK-5.4
# mv UFconfig UMFPACK-5.4
# mv AMD UMFPACK-5.4
# fi

#continue_check
# make
# make cmake
export CXX=g++
cd $GOMA_LIB/cmake-2.8.12.2
if [ -f bin/cmake ]
then
    echo "cmake is already built"
else
    ./bootstrap --prefix=$GOMA_LIB/cmake-2.8.12.2
    make -j$MAKE_JOBS
    make install
fi

export PATH=$GOMA_LIB/cmake-2.8.12.2/bin:$PATH


#continue_check
#make openmpi
cd $GOMA_LIB/openmpi-1.6.4
if [ -f bin/ompi_info ]
then
    echo "openmpi is already built"
else
    ./configure --prefix=$GOMA_LIB/openmpi-1.6.4
    make -j$MAKE_JOBS
    make install
fi
cd ..

#hdf5
if [ -d hdf5-1.8.15/lib ]
then
    echo "hdf5 already built"
else
    mv hdf5-1.8.15 tmpdir
    mkdir hdf5-1.8.15
    mv tmpdir hdf5-1.8.15/src
    cd hdf5-1.8.15/src
    CC=$OPENMPI_TOP/bin/mpicc ./configure --enable-parallel --enable-shared=off --prefix=$GOMA_LIB/hdf5-1.8.15
    make -j$MAKE_JOBS
    make install
    cd ../..
fi

#matio
if [ -d matio-1.5.2/lib ]
then
    echo "matio already built"
else
    mv matio-1.5.2 tmpdir
    mkdir matio-1.5.2
    mv tmpdir matio-1.5.2/src
    cd matio-1.5.2/src
    ./configure --with-hdf5=$GOMA_LIB/hdf5-1.8.15 --prefix=$GOMA_LIB/matio-1.5.2 --enable-shared=off
    make -j$MAKE_JOBS
    make install
    cd ../..
fi

#netcdf
if [ -d netcdf-4.3.3.1/lib ]
then
    echo "netcdf already built"
else
    mv netcdf-4.3.3.1 tmpdir
    mkdir netcdf-4.3.3.1
    mv tmpdir netcdf-4.3.3.1/src
    cd $GOMA_LIB/netcdf-4.3.3.1/src
    cd include
    echo "$NETCDF_PATCH" > netcdf.patch
    patch -f --ignore-whitespace netcdf.h < netcdf.patch
    cd ..
    export CPPFLAGS=-I$GOMA_LIB/hdf5-1.8.15/include
    export LDFLAGS=-L$GOMA_LIB/hdf5-1.8.15/lib 
    echo $CPPFLAGS
    echo $LDFLAGS
    CC=mpicc ./configure --prefix=$GOMA_LIB/netcdf-4.3.3.1 --enable-shared=off --disable-dap --enable-parallel-tests
    make -j$MAKE_JOBS
    make install
    cd ../..
fi

    # CXXFLAGS=-I$GOMA_LIB/hdf5-1.8.12/hdf5/include \
    # FFFLAGS=-I$GOMA_LIB/hdf5-1.8.12/hdf5/include \
    # FCFLAGS=-I$GOMA_LIB/hdf5-1.8.12/hdf5/include \
    # FC=$GOMA_LIB/openmpi-1.6.4/bin/mpif90 \
    # CXX=$GOMA_LIB/openmpi-1.6.4/bin/mpicxx \
    # CC=$GOMA_LIB/openmpi-1.6.4/bin/mpicc \


#continue_check
#make Seacas

# cd $GOMA_LIB/SEACAS-2013-12-03
# if [ -f bin/aprepro ]
# then
#     echo "SEACAS already built"
# else
#     cd TPL/netcdf/netcdf-4.3.2/include
#     echo "$NETCDF_PATCH" > netcdf.patch
#     patch -f --ignore-whitespace netcdf.h < netcdf.patch
#     cd ../../
#     echo "$IMAKE_PATCH" > Imake.patch
#     patch -f Imakefile < Imake.patch
#     cd ../../
#     export ACCESS=$GOMA_LIB/SEACAS-2013-12-03
#     cd ACCESS/itools/config/cf
#     echo "$SITE_PATCH" > site.patch
#     patch -f site.def < site.patch
#     echo "$SEACAS_LINUX_PATCH" > linux.patch
#     patch -f linux.cf < linux.patch
#     cd ../../../../
#     ACCESS/scripts/buildSEACAS -auto
# fi

#continue_check
 #make BLAS
cd $GOMA_LIB/BLAS-3.6.0
if [ -f libblas.a ]
then
    echo "BLAS already built"
else
    make -j$MAKE_JOBS
    cp blas_LINUX.a libblas.a
fi

#continue_check
#make parMetis
cd $GOMA_LIB
if [ -d parmetis-4.0.3/lib ]
then
    echo "ParMetis already Built"
else
    mv parmetis-4.0.3 tmpdir
    mkdir parmetis-4.0.3
    mv tmpdir parmetis-4.0.3/src
    cd parmetis-4.0.3/src
    make config cc=mpicc cxx=mpicxx prefix=$GOMA_LIB/parmetis-4.0.3
    make
    make install
    cd ..
    cp src/metis/include/metis.h include
    cp src/build/Linux-x86_64/libmetis/libmetis.a lib/
fi

#continue_check
#make ARPACK
cd $GOMA_LIB/ARPACK
if [ -f libarpack_x86_64.a ]
then
    echo "ARPACK already built"
else
    echo "$ARMAKE_PATCH" > ARmake.patch
    patch ARmake.inc < ARmake.patch
    make all
    mkdir lib
    cp libarpack_x86_64.a lib/libarpack.a
fi

#continue_check
#make SuperLU
cd $GOMA_LIB/SuperLU_DIST_5.0.0
if [ -f lib/libsuperludist.a ]
then
    echo "SuperLU_DIST already built"
else
    cat > make.inc << EOF
SuperLUroot     =  $GOMA_LIB/SuperLU_DIST_5.0.0
DSUPERLULIB     = $GOMA_LIB/SuperLU_DIST_5.0.0/lib/libsuperludist.a

# BLASDEF       = -DUSE_VENDOR_BLAS

LIBS            = $GOMA_LIB/SuperLU_DIST_5.0.0/lib/libsuperludist.a $GOMA_LIB/BLAS-3.6.0/blas_LINUX.a $GOMA_LIB/parmetis-4.0.3/lib/libparmetis.a $GOMA_LIB/parmetis-4.0.3/lib/libmetis.a -lgfortran

#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         = /usr/bin/ar
ARCHFLAGS    = cr
RANLIB       = /usr/bin/ranlib

CC           = mpicc
CFLAGS       = -O3 -DNDEBUG -DUSE_VENDOR_BLAS -DDEBUGlevel=0 -DPRNTlevel=0 -std=c99 -g -I$GOMA_LIB/parmetis-4.0.3/include
CDEFS = -DAdd_
# CFLAGS       += -D
# CFLAGS       +=
NOOPTS       = -O0
FORTRAN      = /usr/bin/gfortran

LOADER       = mpicc
LOADOPTS     = -Wl,-rpath,$GOMA_LIB/SuperLU_DIST_5.0.0/lib
EOF
    make
fi

#continue_check
#make lapack
cd $GOMA_LIB/lapack-3.2.1
if [ -f liblapack.a ]
then
    echo "LAPACK already built"
else
    mv make.inc.example make.inc
    echo "$LAPACK_PATCH" > make.patch
    patch make.inc < make.patch
    make -j$MAKE_JOBS
    cp lapack_LINUX.a liblapack.a
fi

#continue_check
#make sparse
cd $GOMA_LIB/sparse
if [ -f lib/libsparse.a ]
then 
    echo "Sparse already built"
else
    cd src
    make -j$MAKE_JOBS
    cd ../lib/
    cp sparse.a libsparse.a
    cd ..
fi

#continue_check
#make UMFPACK
# cd $GOMA_LIB/UMFPACK-5.4
# if [ -f Lib/libamd.a ]
# then
#     echo "UMFPACK already built"
# else
#     cd UFconfig
#     echo "$UFCONFIG_PATCH" > UFconfig.patch
#     patch UFconfig.mk < UFconfig.patch
#     cd ../UMFPACK
#     make -j$MAKE_JOBS
#     cd ../
#     mkdir -p Include
#     mkdir -p Lib
#     cp UMFPACK/Include/* Include
#     cp UMFPACK/Lib/* Lib
#     cp AMD/Include/amd.h Include
#     cp AMD/Include/amd_internal.h Include
#     cp AMD/Lib/libamd.a Lib
#     cp UFconfig/UFconfig.h Include
# fi


#make SuiteSparse
cd $GOMA_LIB/SuiteSparse
if [ -f UMFPACK/Lib/libumfpack.a ]
then
    echo "SuiteSparse is already built"
else
    cd SuiteSparse_config
    echo `pwd`
    echo "$SUITESPARSE_CONFIG_PATCH" > SuiteSparse_config.patch
    patch SuiteSparse_config.mk < SuiteSparse_config.patch
    cd ..
    make CC=gcc
    cd UMFPACK/Include
    ln -s ../../SuiteSparse_config/SuiteSparse_config.h UFconfig.h
    ln -s ../../SuiteSparse_config/SuiteSparse_config.h SuiteSparse_config.h

fi
cd ../../..


#continue_check
#make y12m
cd $GOMA_LIB/
if [ -d $GOMA_LIB/y12m ]
then 
    echo "y12m directory already renamed"
else
mv $GOMA_LIB/y12m-1.0 $GOMA_LIB/y12m
fi

cd $GOMA_LIB/y12m
if [ -f liby12m.a ]
then
    echo "y12m already built"
else
    make FC=gfortran
fi

#continue_check
# make scalapack
cd $GOMA_LIB/scalapack-2.0.2
if [ -f libscalapack.a ]
then
    echo "scalapack already built"
else
    cp SLmake.inc.example SLmake.inc
    echo "$SCALAPACK_PATCH" > scalapack.patch
    patch SLmake.inc < scalapack.patch
    make # scalapack only compiles with 1 make job
fi

#continue_check
# make mumps
cd $GOMA_LIB/MUMPS_5.0.1
if [ -f lib/libdmumps.a ]
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
RM = /bin/rm -f
CC = mpicc
FC = mpif90
FL = mpif90
AR = ar vr 
RANLIB = ranlib
SCALAP  = -L$GOMA_LIB/scalapack-2.0.2 -lscalapack

INCPAR = -I$GOMA_LIB/openmpi-1.6.4/include

LIBPAR = \$(SCALAP)  -L$GOMA_LIB/openmpi-1.6.4/lib -lmpi -lmpi_f77

INCSEQ = -I\$(topdir)/libseq
LIBSEQ  =  -L\$(topdir)/libseq -lmpiseq

LIBBLAS = -L$GOMA_LIB/BLAS-3.6.0 -lblas
LIBOTHERS = -lpthread

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
    make -j$MAKE_JOBS
fi

#continue_check
#make trilinos
rm -rf $GOMA_LIB/trilinos-12.6.3-Temp
mkdir $GOMA_LIB/trilinos-12.6.3-Temp
cd $GOMA_LIB/trilinos-12.6.3-Temp

rm -f CMakeCache.txt

FORTRAN_COMPILER=mpif90

FORTRAN_LIBS=-lgfortran

export MPI_BASE_DIR=$GOMA_LIB/openmpi-1.6.4

# Have to set this to get TRY_RUN(...) commands to work
export LD_LIBRARY_PATH=$MPI_BASE_DIR/lib:$LD_LIBRARY_PATH

MPI_LIBS="-LMPI_BASE_DIR/lib -lmpi_f90 -lmpi_f77 -lmpi"
SEACAS_LIBS="-L${GOMA_LIB}/hdf5-1.8.15/lib -lhdf5_hl -lhdf5 -lz -lm"
# Install directory
TRILINOS_INSTALL=$GOMA_LIB/trilinos-12.6.3-Built
#continue_check


cmake \
-D CMAKE_AR=/usr/bin/ar \
-D CMAKE_RANLIB=/usr/bin/ranlib \
-D CMAKE_BUILD_TYPE:STRING=RELEASE \
-D CMAKE_CXX_COMPILER:FILEPATH=$MPI_BASE_DIR/bin/mpiCC \
-D CMAKE_C_COMPILER:FILEPATH=$MPI_BASE_DIR/bin/mpicc \
-D CMAKE_Fortran_COMPILER:FILEPATH=$MPI_BASE_DIR/bin/mpif90 \
-D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
-D Trilinos_ENABLE_Triutils:BOOL=ON \
-D Trilinos_ENABLE_SEACAS:BOOL=ON \
-D Trilinos_ENABLE_Amesos:BOOL=ON \
-D Trilinos_ENABLE_Epetra:BOOL=ON \
-D Trilinos_ENABLE_Xpetra:BOOL=OFF \
-D Trilinos_ENABLE_Ifpack:BOOL=ON \
-D Trilinos_ENABLE_Teuchos:BOOL=ON \
-D Trilinos_ENABLE_ML:BOOL=ON \
-D Trilinos_ENABLE_AztecOO:BOOL=ON \
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
      -D Netcdf_LIBRARY_DIRS:PATH="$GOMA_LIB/netcdf-4.3.3.1/lib" \
      -D TPL_ENABLE_Netcdf:BOOL=ON \
      -D TPL_Netcdf_INCLUDE_DIRS:PATH="$GOMA_LIB/netcdf-4.3.3.1/include" \
      -D Matio_LIBRARY_DIRS:PATH=$GOMA_LIB/matio-1.5.2/lib \
      -D Matio_INCLUDE_DIRS:PATH=$GOMA_LIB/matio-1.5.2/include \
-D TPL_ENABLE_MPI:BOOL=ON \
  -D MPI_COMPILER:FILEPATH=$MPI_BASE_DIR/bin/mpiCC \
  -D MPI_EXECUTABLE:FILEPATH=$MPI_BASE_DIR/bin/mpirun \
  -D MPI_BASE_DIR:PATH=$MPI_BASE_DIR \
-D EpetraExt_BUILD_GRAPH_REORDERINGS:BOOL=ON \
-D EpetraExt_BUILD_BDF:BOOL=ON \
-D TPL_ENABLE_Boost:BOOL=OFF \
-D TPL_ENABLE_LAPACK:BOOL=ON \
-D LAPACK_LIBRARY_DIRS=$GOMA_LIB/lapack-3.2.1 \
-D LAPACK_LIBRARY_NAMES="lapack" \
-D TPL_ENABLE_BLAS:BOOL=ON \
-D BLAS_LIBRARY_DIRS=$GOMA_LIB/BLAS-3.6.0 \
-D BLAS_LIBRARY_NAMES="blas" \
-D CMAKE_INSTALL_PREFIX:PATH=$TRILINOS_INSTALL \
-D Trilinos_EXTRA_LINK_FLAGS:STRING="$FORTRAN_LIBS $SEACAS_LIBS $MPI_LIBS" \
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
  -D SuperLUDist_LIBRARY_DIRS:PATH=$GOMA_LIB/SuperLU_DIST_5.0.0/lib \
  -D SuperLUDist_INCLUDE_DIRS:PATH=$GOMA_LIB/SuperLU_DIST_5.0.0/SRC \
-D TPL_ENABLE_ParMETIS:BOOL=ON \
  -D ParMETIS_LIBRARY_DIRS:PATH=$GOMA_LIB/parmetis-4.0.3/lib \
  -D ParMETIS_INCLUDE_DIRS:PATH=$GOMA_LIB/parmetis-4.0.3/include \
  -D TPL_ParMETIS_INCLUDE_DIRS:PATH=$GOMA_LIB/parmetis-4.0.3/include \
  -D TPL_ENABLE_y12m:BOOL=ON \
  -D y12m_LIBRARY_NAMES:STRING="y12m" \
  -D y12m_LIBRARY_DIRS:PATH=$GOMA_LIB/y12m \
-D TPL_ENABLE_MUMPS:BOOL=ON \
  -D MUMPS_LIBRARY_NAMES:STRING="dmumps;mumps_common;pord" \
  -D MUMPS_LIBRARY_DIRS:PATH="$GOMA_LIB/MUMPS_5.0.1/lib;$GOMA_LIB/MUMPS_5.0.1/PORD/lib" \
  -D MUMPS_INCLUDE_DIRS:PATH="$GOMA_LIB/MUMPS_5.0.1/include;$GOMA_LIB/MUMPS_5.0.1/PORD/include" \
  -D CMAKE_CXX_FLAGS:STRING="-DMUMPS_5_0" \
  -D Amesos_ENABLE_SCALAPACK:BOOL=ON \
  -D SCALAPACK_INCLUDE_DIRS:FILEPATH="$GOMA_LIB/scalapack-2.0.2/SRC" \
  -D SCALAPACK_LIBRARY_DIRS:FILEPATH="$GOMA_LIB/scalapack-2.0.2" \
  -D SCALAPACK_LIBRARY_NAMES:STRING="scalapack" \
-D Amesos_ENABLE_SuperLUDist:BOOL=on \
-D Amesos_ENABLE_ParMETIS:BOOL=on \
-D Amesos_ENABLE_LAPACK:BOOL=ON \
-D Amesos_ENABLE_KLU:BOOL=ON \
-D Amesos_ENABLE_UMFPACK:BOOL=ON \
-D Amesos_ENABLE_MUMPS:BOOL=ON \
$EXTRA_ARGS \
$GOMA_LIB/trilinos-12.6.3-Source



#continue_check
make -j$MAKE_JOBS
#continue_check
make install

