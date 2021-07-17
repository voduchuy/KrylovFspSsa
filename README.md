# KrylovFspSsa

Fortran implementation of the  Krylov-FSP-SSA algorithm. 

Reference:
R.B.Sidje, H.D.Vo. "Solving the chemical master equation by a fast adaptive finite state projection based on the stochastic simulation algorithm". _Mathematical Biosciences_, Vol. 269, November 2015.

## System requirements
This repository has been compiled successfully on ubuntu18.04 with the following libraries (installed with ```sudo apt-get install```):
1) libblas
2) liblapack
3) gfortran 

In addition, you need to build system CMake (version 3.17) that can be downloaded from https://cmake.org/download/.

If your computer uses a different operating system, I recommend trying first to compile this repository within a Docker container.

## Compile the programs 
You need the latest CMake (v 3.17) to properly generate the Makefile for this repository. Once you have installed CMake on your system, take the following steps to compile the modules and the example programs
1) Open a terminal; navigate to this repository folder; Then 
    mkdir build
2) cd build; 
    cmake .. 
3) If things go well, you will see the following outputs from cmake:
```
-- Configuring done 
-- Generating done 
-- Build files have been written to: <repository_folder>/build
```
4) Type make 
5) After compilation, all program files will be available in the directory ```build```.
## Example of solving a model from input file 
Copy the compiled program ```TestSolverFromFile``` from ```build``` into the subfolder ```models``` of this repository. Navigate to ```models``` and execute the program. The program will parse a reaction network model from the file ```toggle_model.input``` and solve the CME defined from that model.
