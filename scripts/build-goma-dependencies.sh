#!/usr/bin/env bash
#
# Goma dependency build script
#
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

GOMA_LIB=$1

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
FORTRAN_COMPILER=$OPENMPI_TOP/bin/mpif90

export MPI_BASE_DIR=$OPENMPI_TOP
export GOMA_LIB
export OWNER
export OPENMPI_TOP
export MAKE_JOBS

ARCHIVE_NAMES=("AMD-2.2.4.tar.gz" \
"arpack96.tar.gz" \
"patch.tar.gz" \
"blas.tgz" \
"cmake-2.8.12.2.tar.gz" \
"hdf5-1.8.13.tar.gz" \
"lapack-3.2.1.tgz" \
"netcdf-4.3.2.tar.gz" \
"openmpi-1.6.4.tar.gz" \
"ParMetis-3.1.1.tar.gz" \
"Seacas2013-12-03.tar.gz" \
"sparse.tar.gz" 
"superlu_dist_2.3.tar.gz" \
"UMFPACK-5.4.0.tar.gz" \
"UFconfig-3.7.1.tar.gz" \
"y12m-1.0.tar.gz" \
"trilinos-11.8.1-Source.tar.gz" \
"scalapack-2.0.2.tgz" \
"MUMPS_4.10.0.tar.gz")

ARCHIVE_MD5SUMS=("31943ff81859eb42d695d8a34d5b4201" \
"fffaa970198b285676f4156cebc8626e" \
"14830d758f195f272b8594a493501fa2" \
"5e99e975f7a1e3ea6abcad7c6e7e42e6" \
"17c6513483d23590cbce6957ec6d1e66" \
"c03426e9e77d7766944654280b467289" \
"a3202a4f9e2f15ffd05d15dab4ac7857" \
"2fd2365e1fe9685368cd6ab0ada532a0" \
"70aa9b6271d904c6b337ca326e6613d1" \
"57318dbaddff2c3d1ef820cff0bf87b0" \
"d7051624547de02ace05b2988c9df865" \
"1566d914d1035ac17b73fe9bc0eed02a" \
"8cf8cb964bb25ff6981f994b7e794f29" \
"4455192ff4bf0fd75796d08c0ed614a6" \
"145b8126b7ee957d1ea6b401be0bf798" \
"eed01310baca61f22fb8a88a837d2ae3" \
"3c9465b6d63d824e9dc0365ca73c3370" \
"2f75e600a2ba155ed9ce974a1c4b536f" \
"959e9981b606cd574f713b8422ef0d9f")

ARCHIVE_URLS=("http://www.cise.ufl.edu/research/sparse/amd/AMD-2.2.4.tar.gz" \
"http://www.caam.rice.edu/software/ARPACK/SRC/arpack96.tar.gz" \
"http://www.caam.rice.edu/software/ARPACK/SRC/patch.tar.gz" \
"http://www.netlib.org/blas/blas.tgz" \
"http://www.cmake.org/files/v2.8/cmake-2.8.12.2.tar.gz" \
"http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.13.tar.gz" \
"http://www.netlib.org/lapack/lapack-3.2.1.tgz" \
"ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.3.2.tar.gz" \
"http://www.open-mpi.org/software/ompi/v1.6/downloads/openmpi-1.6.4.tar.gz" \
"http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/OLD/ParMetis-3.1.1.tar.gz" \
"http://sourceforge.net/projects/seacas/files/latest/download" \
"http://sourceforge.net/projects/sparse/files/sparse/sparse1.4b/sparse1.4b.tar.gz/download" \
"http://crd-legacy.lbl.gov/~xiaoye/SuperLU/superlu_dist_2.3.tar.gz" \
"http://www.cise.ufl.edu/research/sparse/umfpack/UMFPACK-5.4.0.tar.gz" \
"http://www.cise.ufl.edu/research/sparse/UFconfig/UFconfig-3.7.1.tar.gz" \
"http://sisyphus.ru/cgi-bin/srpm.pl/Branch5/y12m/getsource/0" \
"http://trilinos.sandia.gov/download/files/trilinos-11.8.1-Source.tar.gz" \
"http://www.netlib.org/scalapack/scalapack-2.0.2.tgz" \
"http://mumps.enseeiht.fr/MUMPS_4.10.0.tar.gz")

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
> BLASLIB       = -L$(GOMA_LIB)/BLAS -lblas
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
> SCALAP  = $(GOMA_LIB)/scalapack-2.0.2/libscalapack.a $(GOMA_LIB)/BLACS/LIB/libblacs.a $(GOMA_LIB)/BLACS/LIB/libblacsCinit.a $(GOMA_LIB)/BLACS/LIB/libblacsF77init.a
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
> LIBBLAS = -L$(GOMA_LIB)/BLAS -lblas
EOF
MUMPS_PATCH=${MUMPS_PATCH/__GOMA_LIB__/$GOMA_LIB}

read -d '' IMAKE_PATCH << "EOF"
13c13
< V_NUM = netcdf-4.2.1.1
---
> V_NUM = netcdf-4.3.2
EOF


read -d '' NETCDF_PATCH << "EOF"
228c228
< #define NC_MAX_DIMS   1024    
---
> #define NC_MAX_DIMS   65536   
230c230
< #define NC_MAX_VARS   8192    
---
> #define NC_MAX_VARS   524288  
232c232
< #define NC_MAX_VAR_DIMS   1024 /**< max per variable dimensions */
---
> #define NC_MAX_VAR_DIMS   8 /**< max per variable dimensions */
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
104c104
< FC      = f77
---
> FC      = gfortran
115c115
< MAKE    = /bin/make
---
> MAKE    = /usr/bin/make
#130c130
#< #RANLIB  = touch
#---
#>
 
EOF
ARMAKE_PATCH=${ARMAKE_PATCH/__GOMA_LIB5__/$GOMA_LIB}

read -d '' UFCONFIG_PATCH << "EOF"
61c61,64
< ARCHIVE = $(AR) $(ARFLAGS)
---
> AR = ar cr
> ARFLAGS = cr
> ARCHIVE = ar cr
> 
111,113c114,115
< BLAS = -lblas -lgfortran
< LAPACK = -llapack
< 
---
> BLAS = __GOMA_LIB4__/BLAS/blas_LINUX.a -lgfortran
> LAPACK = __GOMA_LIB4__/lapack-3.2.1/lapack_LINUX.a
173c175
< UMFPACK_CONFIG =
---
> UMFPACK_CONFIG = -DNCHOLMOD
214c216
< CHOLMOD_CONFIG =
---
> CHOLMOD_CONFIG = 

EOF
UFCONFIG_PATCH=${UFCONFIG_PATCH/__GOMA_LIB4__/$GOMA_LIB}


read -d '' LAPACK_PATCH << "EOF"
54c54
< BLASLIB      = ../../blas$(PLAT).a
---
> BLASLIB      = __GOMA_LIB3__/BLAS/blas_LINUX.a

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
> BLASLIB       = __GOMA_LIB2__/BLAS -lessl
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
        echo $archive
        continue_check
    fi
}

mkdir -p tars
cd tars
#downloads
count=0
for i in ${ARCHIVE_NAMES[@]}; do
    echo "Check for $i"

    if ! [ -f $i ]
    then
        wget ${ARCHIVE_URLS[count]} -O $i
        mychecksum $i $count
    else
        mychecksum $i $count
    fi
    cd ..
    count=$(( $count + 1 ))

    if ! tar tf $i &> /dev/null; then
      tar -xzf tars/$i
    fi
    cd tars
done

#move
cd ../

#Seacas
mv 2013-12-03 SEACAS-2013-12-03
mv hdf5-1.8.13 hdf5-source
mv hdf5-source SEACAS-2013-12-03/TPL/hdf5/
mv netcdf-4.3.2 SEACAS-2013-12-03/TPL/netcdf/

#UMFPACK
mkdir UMFPACK-5.4
mv UMFPACK UMFPACK-5.4
mv UFconfig UMFPACK-5.4
mv AMD UMFPACK-5.4

# make
# make cmake
cd cmake-2.8.12.2
./bootstrap --prefix=$GOMA_LIB/cmake-2.8.12.2
make -j$MAKE_JOBS
make install
cd ../

#make openmpi
cd openmpi-1.6.4
./configure --prefix=$GOMA_LIB/openmpi-1.6.4
make -j$MAKE_JOBS
make install
cd ../

#make Seacas
cd SEACAS-2013-12-03/TPL/netcdf/netcdf-4.3.2/include
echo "$NETCDF_PATCH" > netcdf.patch
patch netcdf.h < netcdf.patch
cd ../../
echo "$IMAKE_PATCH" > Imake.patch
patch Imakefile < Imake.patch
cd ../../
export ACCESS=$GOMA_LIB/SEACAS-2013-12-03
cd ACCESS/itools/config/cf
echo "$SITE_PATCH" > site.patch
patch site.def < site.patch
echo "$SEACAS_LINUX_PATCH" > linux.patch
patch linux.cf < linux.patch
cd ../../../../
ACCESS/scripts/buildSEACAS -auto
cd ../

 #make BLAS
cd BLAS
make -j$MAKE_JOBS
cp blas_LINUX.a libblas.a
cd ../

#make parMetis
cd ParMetis-3.1.1
make -j$MAKE_JOBS
cd ../  

#make ARPACK
cd ARPACK
echo "$ARMAKE_PATCH" > ARmake.patch
patch ARmake.inc < ARmake.patch
make all
mkdir lib
cp libarpack_x86_64.a lib/libarpack.a
cd ..

#make SuperLU
cd SuperLU_DIST_2.3
make PLAT= DSuperLUroot=$GOMA_LIB/SuperLU_DIST_2.3 BLASLIB="$GOMA_LIB/BLAS/blas_LINUX.a" METISLIB="-L$GOMA_LIB/ParMetis-3.1.1 -lmetis" PARMETISLIB="-L$GOMA_LIB/ParMetis-3.1.1 -lparmetis" ARCHFLAGS=-rc RANLIB=echo CC=mpicc CFLAGS=-D_SP NOOPTS= FORTRAN=gfortran FFLAGS="-O3 -Q" LOADER=mpicxx LOADOPTS=-lgfortran CDEFS="-DAdd_"
cd lib/
cp libsuperlu_dist_2.3.a libsuperludist.a
cd ../../

#make lapack
cd lapack-3.2.1
mv make.inc.example make.inc
echo "$LAPACK_PATCH" > make.patch
patch make.inc < make.patch
make -j$MAKE_JOBS
cp lapack_LINUX.a liblapack.a
cd ..

#make sparse
cd sparse/src
make -j$MAKE_JOBS
cd ../lib/
cp sparse.a libsparse.a
 cd ../../

#make UMFPACK
cd UMFPACK-5.4/UFconfig
echo "$UFCONFIG_PATCH" > UFconfig.patch
patch UFconfig.mk < UFconfig.patch
cd ../UMFPACK
make -j$MAKE_JOBS
cd ../
mkdir -p Include
mkdir -p Lib
cp UMFPACK/Include/* Include
cp UMFPACK/Lib/* Lib
cp AMD/Include/amd.h Include
cp AMD/Include/amd_internal.h Include
cp AMD/Lib/libamd.a Lib
cp UFconfig/UFconfig.h Include
cd ../

#make y12m
mv y12m-1.0 y12m
cd y12m
make FC=gfortran
cd ..

# make scalapack
cd scalapack-2.0.2
cp SLmake.inc.example SLmake.inc
echo "$SCALAPACK_PATCH" > scalapack.patch
patch SLmake.inc < scalapack.patch
make # scalapack only compiles with 1 make job
cd ..

# make mumps
cd MUMPS_4.10.0
cp Make.inc/Makefile.gfortran.PAR Makefile.inc
echo "$MUMPS_PATCH" > mumps.patch
patch Makefile.inc < mumps.patch
make -j$MAKE_JOBS
cd ..

#make trilinos
mkdir trilinos-11.8.1-Temp
cd trilinos-11.8.1-Temp

rm -f CMakeCache.txt

FORTRAN_COMPILER=mpif90

FORTRAN_LIBS=-lgfortran

export MPI_BASE_DIR=$GOMA_LIB/openmpi-1.6.4

# Have to set this to get TRY_RUN(...) commands to work
export LD_LIBRARY_PATH=$MPI_BASE_DIR/lib:$LD_LIBRARY_PATH

export PATH=$GOMA_LIB/cmake-2.8.12.2/bin:$PATH

MPI_LIBS="-LMPI_BASE_DIR/lib -lmpi_f90 -lmpi_f77 -lmpi"
# Install directory
TRILINOS_INSTALL=$GOMA_LIB/trilinos-11.8.1-Built

cmake \
-D CMAKE_AR=/usr/bin/ar \
-D CMAKE_RANLIB=/usr/bin/ranlib \
-D CMAKE_BUILD_TYPE:STRING=RELEASE \
-D CMAKE_CXX_COMPILER:FILEPATH=$MPI_BASE_DIR/bin/mpiCC \
-D CMAKE_C_COMPILER:FILEPATH=$MPI_BASE_DIR/bin/mpicc \
-D CMAKE_Fortran_COMPILER:FILEPATH=$MPI_BASE_DIR/bin/mpif90 \
-D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
-D Trilinos_ENABLE_Triutils:BOOL=ON \
-D Trilinos_ENABLE_Amesos:BOOL=ON \
-D Trilinos_ENABLE_Epetra:BOOL=ON \
-D Trilinos_ENABLE_Ifpack:BOOL=ON \
-D Trilinos_ENABLE_Teuchos:BOOL=ON \
-D Trilinos_ENABLE_AztecOO:BOOL=ON \
-D Trilinos_ENABLE_TESTS:BOOL=ON \
-D Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
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
-D BLAS_LIBRARY_DIRS=$GOMA_LIB/BLAS \
-D BLAS_LIBRARY_NAMES="blas" \
-D CMAKE_INSTALL_PREFIX:PATH=$TRILINOS_INSTALL \
-D Trilinos_EXTRA_LINK_FLAGS:STRING="$FORTRAN_LIBS;$MPI_LIBS" \
-D TPL_ENABLE_UMFPACK:BOOL=ON \
  -D UMFPACK_LIBRARY_NAMES:STRING="umfpack;amd" \
  -D UMFPACK_LIBRARY_DIRS:PATH=$GOMA_LIB/UMFPACK-5.4/Lib \
  -D UMFPACK_INCLUDE_DIRS:PATH=$GOMA_LIB/UMFPACK-5.4/Include \
-D TPL_ENABLE_AMD:BOOL=ON \
  -D AMD_LIBRARY_NAMES:STRING="amd" \
  -D AMD_LIBRARY_DIRS:PATH=$GOMA_LIB/UMFPACK-5.4/Lib \
  -D AMD_INCLUDE_DIRS:PATH=$GOMA_LIB/UMFPACK-5.4/Include \
-D TPL_ENABLE_SuperLUDist:BOOL=ON \
  -D SuperLUDist_LIBRARY_NAMES:STRING="superludist" \
  -D SuperLUDist_LIBRARY_DIRS:PATH=$GOMA_LIB/SuperLU_DIST_2.3/lib \
  -D SuperLUDist_INCLUDE_DIRS:PATH=$GOMA_LIB/SuperLU_DIST_2.3/SRC \
-D TPL_ENABLE_ParMETIS:BOOL=ON \
  -D ParMETIS_LIBRARY_DIRS:PATH=$GOMA_LIB/ParMetis-3.1.1\
  -D ParMETIS_INCLUDE_DIRS:PATH=$GOMA_LIB/ParMetis-3.1.1 \
  -D TPL_ParMETIS_INCLUDE_DIRS:PATH=$GOMA_LIB/ParMetis-3.1.1 \
  -D TPL_ENABLE_y12m:BOOL=ON \
  -D y12m_LIBRARY_NAMES:STRING="y12m" \
  -D y12m_LIBRARY_DIRS:PATH=$GOMA_LIB/y12m \
-D TPL_ENABLE_MUMPS:BOOL=ON \
  -D MUMPS_LIBRARY_NAMES:STRING="dmumps;mumps_common;pord" \
  -D MUMPS_LIBRARY_DIRS:PATH=$GOMA_LIB/MUMPS_4.10.0/lib \
  -D MUMPS_INCLUDE_DIRS:PATH=$GOMA_LIB/MUMPS_4.10.0/include \
  -D CMAKE_CXX_FLAGS:STRING="-DMUMPS_4_9" \
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
$GOMA_LIB/trilinos-11.8.1-Source

make -j$MAKE_JOBS

make install

