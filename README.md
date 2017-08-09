QuartetScores
=========================

Code for computing various support scores for internodes.

See http://www.biorxiv.org/content/early/2017/07/27/168526 for the preprint.

Download & Install
-------------------------

You can download the source code via the green "Clone or download" button on
[this GitHub page](https://github.com/algomaus/QuartetScores), then clicking on "Download ZIP".
Extract the zip file to your desired destination.

Alternatively, if you use Git, you can clone the Git repository:

    cd path/to/destination
    git clone --recursive https://github.com/algomaus/QuartetScores.git

Then, in order to compile the program, change to the directory that you got via unzipping
or cloning the repository, and build the program:

    cd path/to/destination/QuartetScores/
    make

If needed, this will also download any dependencies.

The compiled program is located in `bin`. To change to that directory, call:

    cd path/to/destination/QuartetScores/bin/

Then, the program can be called as described in the following section.

Usage & Command Line Options
-------------------------

The command line options of the program are:

    ./QuartetScores  [-s] [-v] [-t <number>] -o <file_path> -e <file_path> -r <file_path> [--version] [-h]

Where:

`-r <file_path>`,  `--ref <file_path>`: (required)  Path to the reference tree

`-e <file_path>`,  `--eval <file_path>`: (required)  Path to the evaluation trees

`-o <file_path>`,  `--output <file_path>`: (required)  Path to the output file

`-s`, `--savemem`: Consume less memory, but with the cost of increased runtime (~50% more)

`-v`,  `--verbose`: Verbose mode

`-t <number>`,  `--threads <number>`: Maximum number of threads to use

`--version`: Displays version information and exits.

`-h`,  `--help`: Displays usage information and exits.

Requirements & Troubleshooting
-------------------------

In order to compile the program, you need:

 *  [Make](https://www.gnu.org/software/make/) and [CMake](https://cmake.org/) 2.8.7 or higher.
 *  A fairly up-to-date C++11 compiler, e.g., [GCC](https://gcc.gnu.org/) 4.9 or higher,
    or [clang++](http://clang.llvm.org/) 3.6 or higher.

On typical Linux distributions (e.g., Ubuntu), those programs can be installed via

    sudo apt-get install build-essential cmake

By default, we compile a fully static program, which has the advantage that you are able to move and
copy the program to other locations. However, Clang seems to have some trouble statically linking
OpenMP into the program. If you use Clang and get errors in the linking step, either try to use a
different compiler, like GNU, or try to deactivate static linking of system libraries.
To do so, instead of the build process described above, use this to compile the program:

    cd path/to/destination/QuartetScores/
    mkdir build
    cd build
    cmake -DBUILD_STATIC=OFF ..
    make
    cd ..

This will most likely prevent you from moving the program to other computers, but should work at
least on your machine.
