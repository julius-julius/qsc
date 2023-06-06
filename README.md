# qsc
A numerical non-perturbative quantum spectral curve solver and database of states.

## Installation and Compilation

The first step is to install the requisite C++ compilers and packages and compile the C++ code in the /core directory.

### Compilers and Packages

The `g++` compiler, and the `libcln` and `libcln-dev` packages are necessary to be installed. Other compliers may be used depending on taste and packages may be installed using a package manager. 
In order to use the packages on a Windows operating system, one needs to install a Windows Subsystem for Linux (WSL), though one must bare in mind to type `wsl` before he Shell command. Below we present the steps to install the above on a Debian linux system. The following commands should be run on the terminal 

`> sudo apt install g++`

`> sudo apt install libcln6 `

`> sudo apt install libcln-dev`

### Compilation

In the /core directory, there are four *.cpp* 
