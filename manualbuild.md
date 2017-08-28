#### Build the TPLs

Use these commands:

`cd goma`

`gomadir=$(pwd)`

Set a temporary variable that helps the user complete the install at their site.

`mkdir TPLs`

`cd scripts`

run if your compiler supports c++11:

`nice ./build-goma-dep-trilinos-12.sh -j2 ${gomadir}/TPLs`

or for older compilers (deprecated):

`nice ./build-goma-dep-trilinos-12-noc++11.sh -j2 ${gomadir}/TPLs`

This will take a long time (hours). `nice` lowers the process priority, effectively backgrounding the build script and allowing continued use of the computer. The -j2 flag indicates that two processes will be used to build the TPLs.

`cd ..`

Test that these libraries built successfully by checking for aprepro.

`${gomadir}/TPLs/trilinos-12.10.1/bin/aprepro -v`

If the libraries built successfully, this will print the version of the SEACAS Algebraic Preprocessor, if not address any errors before moving on.

### Configure and build Goma

Now follow these commands/instructions to build and install Goma:

* `mkdir build`
* `cd build`
* `cp ../cmake-config-example ../cmake-config`
* Copy the output of `echo ${gomadir}/TPLs` to the clipboard. Edit `../cmake-config`, find `GOMA_LIBS=` and paste the copied path after the `=` sign. This sets the value of the GOMA_LIBS variable to the path the TPLs were installed to.
* `../cmake-config`
* Either: `make -j 2` To make Goma.
* Or: `make gomad -j 2` To make the debug version
* `make install`

### Run Goma

Many uses of Goma require the helper program `aprepro` to be on the path, ALSO the OpenMPI libraries must be on the library path. Make these environment changes


