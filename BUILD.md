Overview
===============================

This document explains how to build the program, in case that the pre-build binaries do not work
or are not available for your system. We currently provide binaries for Linux and Mac systems,
which are available at the
[Releases Page](https://github.com/algomaus/QuartetScores/releases) of the
[GitHub Repository](https://github.com/algomaus/QuartetScores).

Thus, first try to use these. If this does not work, follow the instructions below.

Get the source code
===============================

First, we need to get the source code for the program.

You can download the source code via the green "Clone or download" button at the
[GitHub Repository](https://github.com/algomaus/QuartetScores), then clicking on "Download ZIP".
Extract the zip file to your desired destination.

Alternatively, if you use Git, you can clone the Git repository:

    cd path/to/destination
    git clone --recursive https://github.com/algomaus/QuartetScores.git

Then, follow the instructions for either Mac or Linux systems.

Mac OSX
===============================

These are instructions on how to build the program on Mac OSX systems.
It was tested on OSX 10.12.6 (Sierra).

Unfortunately, because Mac OSX neither supports OpenMP nor static linking,
this process is quite complex.

Prerequisites
-------------------------------

First, we need CMake:

    brew install cmake

Then, we install a LLVM/Clang version that supports OpenMP, because the default one doesn't:

    brew install llvm

Without OpenMP, the program will be considerably slower, so this is highly recommended.
The printed output of this step is needed in the following step.

Change CMakeLists.txt
-------------------------------

This step is necessary if the installation of LLVM above printed a message similar to
**"This formula is keg-only"**. This means that LLVM was only installed locally,
so that we have to tell CMake where to find it.

In the main directory of QuartetScores (the one you cloned or unzipped above),
add the following lines to `CMakeLists.txt` to the end of the section "Compiler and Linker Options",
just above the section "Dependencies":

    # Paths to LLVM directories. These are the ones printed when installing LLVM.
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L/usr/local/opt/llvm/lib -Wl,-rpath,/usr/local/opt/llvm/lib ")
    set(CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS}        -I/usr/local/opt/llvm/include")

    # Somehow, CMake does not find OpenMP, although it is there.
    # Thus, hard code the necesary flags. Ugly, but seems to work.
    add_definitions( "-DGENESIS_OPENMP" )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp=libomp" )

    # Lastly, add libffi. We might not need it, but it also doesn't hurt to add it.
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L/usr/local/opt/libffi/lib")

It is possible that the paths in the first two `set` instructions differ from what we used here.
Check the output of the LLVM install step for the correct paths to use.
It will (hopefully) tell you which "LDFLAGS" and which "CPPFLAGS" to use.
Use those in the first two `set` lines in the above snippet. Adjust them accordingly:

  * "LDFLAGS" (potentially twice, where the second occurrence only contains some part of the path):
    This has to go in `CMAKE_EXE_LINKER_FLAGS`, replacing the parts that start with `-L/usr/...`,
    until the end of the line.
  * "CPPFLAGS": This has to be put in `CMAKE_CXX_FLAGS`, replacing the part starting with `-I/usr/...`.

Unfortunately, this step is quite complex. This is due to OSX not supporting OpenMP out of the box.

Build the binary
-------------------------------

Now its time to compile the code. We run CMake by hand in order to tell it to use our new Clang,
and to turn off static linking (which is not supported on OSX):

    cd path/to/destination/QuartetScores
    mkdir build
    cd build
    cmake -DCMAKE_C_COMPILER=/usr/local/opt/llvm/bin/clang -DCMAKE_CXX_COMPILER=/usr/local/opt/llvm/bin/clang++ -DBUILD_STATIC=OFF ..
    make

Here again you might need to change the paths to Clang.

If you want to compile with GCC instead (which has to be installed first), use this fourth line instead:

    cmake -DCMAKE_C_COMPILER=gcc-7 -DCMAKE_CXX_COMPILER=g++-7 ..

This is it. The program should now be located in `bin`. To change to that directory, call:

    cd path/to/destination/QuartetScores/bin/

Then, the program can be called as described in the `README.md` document.

Change dynamic linking
-------------------------------

If you additionally want to be able to copy the program to other Mac computers, we need to make sure
that all its dependencies are also copied with the program.

To check, run (from the main directory):

    cd path/to/destination/QuartetScores/bin
    otool -L QuartetScores

Which should yield something like:

    QuartetScores:
    /usr/local/opt/llvm/lib/libomp.dylib (compatibility version 5.0.0, current version 5.0.0)
    /usr/local/opt/llvm/lib/libc++.1.dylib (compatibility version 1.0.0, current version 1.0.0)
    /usr/lib/libSystem.B.dylib (compatibility version 1.0.0, current version 1238.0.0)
    /usr/lib/libc++abi.dylib (compatibility version 1.0.0, current version 307.2.0)

The first two libraries are in our local LLVM/Clang directory (which we installed above),
which is not good for portability.
Thus, we copy them to the `bin` directory, and change all linking paths so that they are relative:

    cp /usr/local/opt/llvm/lib/libomp.dylib .
    cp /usr/local/opt/llvm/lib/libc++.1.dylib .

    install_name_tool -change /usr/local/opt/llvm/lib/libomp.dylib @executable_path/libomp.dylib QuartetScores
    install_name_tool -change /usr/local/opt/llvm/lib/libc++.1.dylib @executable_path/libc++.1.dylib QuartetScores

    chmod 664 libomp.dylib libc++.1.dylib
    install_name_tool -id "@loader_path/libomp.dylib" libomp.dylib
    install_name_tool -id "@loader_path/libc++.1.dylib" libc++.1.dylib
    chmod 444 libc++.1.dylib libomp.dylib

Now, we should see different paths:

    otool -L QuartetScores

yields:

    QuartetScores:
    @executable_path/libomp.dylib (compatibility version 5.0.0, current version 5.0.0)
    @executable_path/libc++.1.dylib (compatibility version 1.0.0, current version 1.0.0)
    /usr/lib/libSystem.B.dylib (compatibility version 1.0.0, current version 1238.0.0)
    /usr/lib/libc++abi.dylib (compatibility version 1.0.0, current version 307.2.0)

and:

    otool -L libomp.dylib
    otool -L libc++.1.dylib

gives:

    libomp.dylib:
    @loader_path/libomp.dylib (compatibility version 5.0.0, current version 5.0.0)
    /usr/lib/libSystem.B.dylib (compatibility version 1.0.0, current version 1238.50.2)
    libc++.1.dylib:
    @loader_path/libc++.1.dylib (compatibility version 1.0.0, current version 1.0.0)
    /usr/lib/libc++abi.dylib (compatibility version 1.0.0, current version 307.3.0)
    /usr/lib/libSystem.B.dylib (compatibility version 1.0.0, current version 1238.50.2)

That looks better. Only dependencies to system libraries, which should be there on other systems, too.

Thus, we are done. If you now want to copy the program to other Mac computers, make sure that you
also copy the two libraries `libomp` and `libc++`, which are now also located in the `bin` directory.

For reference, here are the sources we used to get this to work:

  * https://blogs.oracle.com/dipol/dynamic-libraries,-rpath,-and-mac-os
  * http://thecourtsofchaos.com/2013/09/16/how-to-copy-and-relink-binaries-on-osx/
  * https://www.mikeash.com/pyblog/friday-qa-2009-11-06-linking-and-install-names.html

Linux
===============================

These are instructions on how to build the program on Linux/Unix systems.
It was tested on Ubuntu 14.04.

Prerequisites
-------------------------------

In order to compile the program, you need:

 *  [Make](https://www.gnu.org/software/make/) and [CMake](https://cmake.org/) 2.8.7 or higher.
 *  A fairly up-to-date C++11 compiler, e.g., [GCC](https://gcc.gnu.org/) 4.9 or higher,
    or [clang++](http://clang.llvm.org/) 3.6 or higher.

On typical Linux distributions (e.g., Ubuntu), those programs can be installed via

    sudo apt-get install build-essential cmake

Build the binary
-------------------------------

Change to the directory that you got via unzipping or cloning the repository, and build the program:

    cd path/to/destination/QuartetScores/
    make

If needed, this will also download any dependencies.

The compiled program is located in `bin`. To change to that directory, call:

    cd path/to/destination/QuartetScores/bin/

Then, the program can be called as described in the `README.md` document.

Troubleshooting
-------------------------------

By default, on Linux, we compile a fully static program, which has the advantage that you are able
to move and copy the program to other locations.
However, Clang seems to have some trouble statically linking OpenMP into the program.
If you use Clang and get errors in the linking step, either try to use a different compiler,
like GCC, or try to deactivate static linking of system libraries.
To do so, instead of the build process described above, use this to compile the program:

    cd path/to/destination/QuartetScores/
    mkdir build
    cd build
    cmake -DBUILD_STATIC=OFF ..
    make
    cd ..

This will most likely prevent you from moving the program to other computers, but should work at
least on your machine.
