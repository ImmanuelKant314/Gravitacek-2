# About Gravitáček 2

Gravitáček 2 is a program focused on studying geodesic dynamics around black holes 
with rings or discs. The program is basicaly integrator with some build-in 
equations of motion and methods to study chaos.

As number 2 suggest, this program is continuation of the program *Gravitáček* by Miroslav Žáček.
There were two reasons to write new software. First, old program is 
over 20 years old,
and therefore in some ways is little bit outdated. Secondly, some parts of 
numerical side were not writen optimally and it could be improved. The result
is this software.

# Setting up

## Downloading and building

First create folder, for example `gravitacek2`. To this folder download/save code 
from this project. To build a project, make new directory `build` and build 
the project there using [CMake](https://cmake.org/).
```bash
mkdir build
cd build
cmake ..
cmake --build .
```

## Testing

For testing [Google Test Framework](https://github.com/google/googletest) is used. 
To test code yourself first build directory `external` and to this directory 
download *Google testing framework*.
```bash
mkdir external
cd external
git clone https://github.com/google/googletest.git -b v1.15.2
```
Now rebuild project
```bash
cd ../build
cmake --build .
```
To run test use command `ctest`.

## Create documentation

Documentation can be made generated using [doxygen](https://doxygen.nl/). 
First you need to install *doxygen*, this can be done using [tutorial](https://doxygen.nl/manual/install.html). 
After that go to folder `docs` and build documentation using *doxygen* and Doxyfile.
```bash
cd docs
doxygen Doxyfile
```
If you want to build pdf version of documentation, go to folder `latex` and there
build pdf-document.
```bash
cd latex
make
```

# Using program Gravitáček 2

## Basic idea

Gravitáček is very simple program in which we cand do two things. The first one 
is defining macros. The second one is then running implemented functions. We will
discuss both concepts in this document. But the most important thing to mention,
if we want to exit the program, just type `end` or `END` or `exit`.

## Macros

## Using functions

