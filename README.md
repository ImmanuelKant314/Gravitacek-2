# About Gravitáček 2

Gravitáček 2 is a program for studying chaos in geodesic dynamics in different 
spacetimes. The software enables integration of geodesic motion and some 
techniques for studying chaos.

This software was written for my master's thesis and was heavily based on the 
software written for my bachelor's thesis. The whole project is then inspired by
program *Gravitáček* by Miroslav Žáček.

# Setting up

## Prerequisites

To get the full support for this code the following is assumed:
- usage of Linux operating system (for Windows you can use
[WSL](https://github.com/microsoft/WSL))
- [CMake](https://cmake.org/) for building the software
- [Google Test Framework](https://github.com/google/googletest) for testing the
software
- [GSL](https://www.gnu.org/software/gsl/) for matrix operations
- [doxygen](https://doxygen.nl/) for creating documentation
- [LaTeX](https://en.wikipedia.org/wiki/LaTeX) for creating pdf documentation
(packages [amsmath](https://ctan.org/pkg/amsmath),
[physics](https://ctan.org/pkg/physics) and
[tensor](https://ctan.org/pkg/tensor))

## Downloading, building and testing the project

First create a folder, for example `gravitacek2` and download/save the
code from this project inside the folder. You can do this as
```bash
git clone https://github.com/ImmanuelKant314/Gravitacek-2
```

Before compiling the project you have to install
[GSL](https://www.gnu.org/software/gsl/) and [CMake](https://cmake.org/). 
After installing those you can download [Google Test
Framework](https://github.com/google/googletest) into folder `external` as follows:
```bash
mkdir external                                                  # make directory
cd external                                                     # go to the directory
git clone https://github.com/google/googletest.git -b v1.15.2   # download code
cd ..                                                           # go back
```

To build a project, make new directory `build` and build the project there using
[CMake](https://cmake.org/).
```bash
mkdir build     # create directory
cd build        # go to the directory
cmake ..        # prepare build
cmake --build . # build
```

To check if the code actually runs you can do a test using [Google Test
Framework](https://github.com/google/googletest). To run the test use command
`ctest`, but preferably use `ctest --timeout 10` to make sure that if there is a
problem, then each test runs 10 seconds at most.

## Create documentation

Documentation can be generated using [doxygen](https://doxygen.nl/). First
you need to install *doxygen*, this can be done using
[tutorial](https://doxygen.nl/manual/install.html). After that go to folder
`docs` and build documentation using *doxygen* and Doxyfile.
```bash
cd docs
doxygen Doxyfile
```
To build the pdf version (which is encouraged) you can use
[LaTeX](https://en.wikipedia.org/wiki/LaTeX). This can be done as follows:
```bash
cd latex
make
```

# Using program Gravitáček 2

## Basic idea

Gravitáček is a very simple program in which you can do two things: the first
one is defining macros and the second one is running implemented functions. Both
concepts are discussed in this document. The most important thing to mention is
how to exist the program - just type `end` or `END` or `exit`.

## Macros

## Using functions
