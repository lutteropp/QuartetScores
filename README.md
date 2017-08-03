QuartetScores
=========================

Code for computing various support scores for internodes.

See http://www.biorxiv.org/content/early/2017/07/27/168526 for the preprint.

Download & Install
-------------------------

You can download the source code via the green "Clone or download" button on
[GitHub](https://github.com/algomaus/QuartetScores), then clicking on "Download ZIP".
Extract the zip file to your desired destination.

Alternatively, if you use Git, you can clone the Git repository:

    git clone --recursive https://github.com/algomaus/QuartetScores.git

Then, in order to compile the program, run

    make

in the main directory that you got via unzipping or cloning the repository.
If needed, this will also download any dependencies.

Requirements:

 *  [Make](https://www.gnu.org/software/make/) and [CMake](https://cmake.org/) 2.8.7 or higher.
 *  A fairly up-to-date C++11 compiler, e.g., [clang++](http://clang.llvm.org/) 3.6 or higher,
    or [GCC](https://gcc.gnu.org/) 4.9 or higher.

The compiled program is then located in `bin`.


Usage
-------------------------

The program can be called as follows:

`./QuartetScores  [-s] [-v] [-t <uint>] -o <string> -e <string> -r <string> [--version] [-h]`


Where:

`-r <string>`,  `--ref <string>`: (required)  Path to the reference tree

`-e <string>`,  `--eval <string>`: (required)  Path to the evaluation trees

`-o <string>`,  `--output <string>`: (required)  Path to the output file

`-s`, `--savemem`: Consume less memory, but with the cost of increased runtime

`-v`,  `--verbose`: Verbose mode

`-t <uint>`,  `--threads <uint>`: Maximum number of threads to use

`--version`: Displays version information and exits.

`-h`,  `--help`: Displays usage information and exits.
