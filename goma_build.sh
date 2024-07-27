#!/usr/bin/env bash
#

###    cmake -B <build directory> <path to goma src>
###  # cmake -B <build directory> <path to goma src> -DMDE=10 -DMAX_PROB_VAR=23 ...more options
###    make -C <build directory>

		####then make clean -C build/votd
##Btype=Release
Btype=RelWithDebInfo
cflag=0
G_Defs="-DMDE=27 -DMAX_PROB_VAR=16 -DMAX_EXTERNAL_FIELD=4 -DMAX_CONC=6"  
for arg in $*
    do
	if [ "$arg" = debug ]
		then Btype=Debug
	    elif [ "$arg" = asan ]
		then Btype=Asan
	    elif [ "$arg" = clean ]
		then cflag=1
	    elif [ "$arg" = guts ]
		then
		G_Defs="-DMDE=27 -DMAX_PROB_VAR=16 -DMAX_EXTERNAL_FIELD=4 -DMAX_CONC=6"  
	    elif [ "$arg" = votd ]
		then 
		G_Defs="-DMDE=27 -DMAX_PROB_VAR=16 -DMAX_EXTERNAL_FIELD=6 -DMAX_CONC=8"  
	    elif [ "$arg" = 3D ]
		then
		G_Defs="-DMDE=10 -DMAX_PROB_VAR=10 -DMAX_EXTERNAL_FIELD=4 -DMAX_CONC=4"  
	    elif [ "$arg" = ve ]
		then
		G_Defs="-DMDE=27 -DMAX_PROB_VAR=44 -DMAX_EXTERNAL_FIELD=6 -DMAX_CONC=8"  
	    elif [ "$arg" = tst ]
		then 
		G_Defs="-DMDE=27 -DMAX_PROB_VAR=16 -DMAX_EXTERNAL_FIELD=6 -DMAX_CONC=8"  
            elif [ "$arg" = mls ]
                then
                G_Defs="-DMDE=12 -DMAX_PROB_VAR=16 -DMAX_EXTERNAL_FIELD=6 -DMAX_CONC=8"
	    else echo "No Match"
	fi

    done
if [ $cflag -eq 1 ]
    then 
    echo "cleaning $arg build"
    make clean -C build/$arg
    rm -f build/$arg/CMakeCache.txt
fi
echo "building $arg"
cmake -B build/$arg build/$arg  -DCMAKE_BUILD_TYPE=$Btype  $G_Defs 
make -C build/$arg  -j8
echo "Done with $arg build"
