QuartetScores
=========================

Code for computing various support scores for internodes.
See http://www.biorxiv.org/content/early/2017/07/27/168526 for the preprint.

Download & Install
-------------------------

You can download the program via the green "Clone or download" button on
[GitHub](https://github.com/algomaus/QuartetScores), then clicking on "Download ZIP".
Extract the zip file to your desired destination.

Alternatively, you can clone the Git repository:

    git clone --recursive https://github.com/algomaus/QuartetScores.git

Then, in order to compile the program, run

    make

in the main directory that you got via unzipping or cloning the repository.

Requirements:

 *  [Make](https://www.gnu.org/software/make/) and [CMake](https://cmake.org/) 2.8.7 or higher.
 *  A fairly up-to-date C++11 compiler, e.g., [clang++](http://clang.llvm.org/) 3.6 or higher,
    or [GCC](https://gcc.gnu.org/) 4.9 and higher.

The program is then located in `bin`.


Usage
-------------------------

`./QuartetScores  [-s] [-v] [-t <uint>] -o <string> -e <string> -r <string> [--] [--version] [-h]`


Where:

`-e <string>`,  `--eval <string>`: (required)  Path to the evaluation trees

`-r <string>`,  `--ref <string>`: (required)  Path to the reference tree

`-o <string>`,  `--output <string>`: (required)  Path to the output file

`-s`, `--savemem`: Consume less memory, but with the cost of increased runtime

`-v`,  `--verbose`: Verbose mode

`-t <uint>`,  `--threads <uint>`: Maximum number of threads to use

`--`,  `--ignore_rest`: Ignores the rest of the labeled arguments following this flag.

`--version`: Displays version information and exits.

`-h`,  `--help`: Displays usage information and exits.
