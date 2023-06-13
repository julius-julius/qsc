# qsc
A numerical non-perturbative quantum spectral curve solver and database of states.

## Installation and Compilation

The first step is to install the requisite C++ compilers and packages and compile the C++ code in the /core directory.

### Compilers and Packages

The `g++` compiler, and the `libcln` and `libcln-dev` packages are necessary to be installed. Other compliers may be used depending on taste and packages may be installed using a package manager. 
In order to use the packages on a Windows operating system, one needs to install a Windows Subsystem for Linux (WSL), though one must bare in mind to type `wsl` before any Shell command. Below we present the steps to install the above on a Debian linux system. The following commands should be run on the terminal 

`> sudo apt install g++`

`> sudo apt install libcln6 `

`> sudo apt install libcln-dev`

### Compilation

In the /core directory, there are four *.cpp* codes. Each code pertains to a particular type of state in planar N = 4 supersymmetric Yang-Mills theory (SYM):

- For type I states, use *TypeI_core.cpp*
- For type II states, use *TypeII_core.cpp*
- For type III states, use *TypeIII_core.cpp*
- For type IV states, use *TypeIV_core.cpp*

To compile any of the above codes, use the following command where *source.cpp* is the C++ file, and *executable.out* is the executable output file

`> g++ source.cpp -lm -lcln -o executable.out`

In order for the executables to be compatible with the Mathematica notebooks in the /core directory, the names of the execulatbles in the four cases should be as follows.

- For type I states, use *TypeI_exec.out*
- For type II states, use *TypeII_exec.out*
- For type III states, use *TypeIII_exec.out*
- For type IV states, use *TypeIV_exec.out*

Once compiled (it should compile without errors, warnings are OK), ensure that the exeutable file and the Mathematica notebook are in the same directory. The C++ sourcecode needs to be compiled only once. 

## Execution and Automation

In order to ensure that everything goes smoothly, and the code was compiled properly, we recommend to test the executables by running the corresponding Mathematica notebooks in /cores. 

### Execute from Mathematica notebook

Ensure that the executable file is in the same directry as the Mathematica notebooks available in the /core directory. Depending on your application, you should use

- *TypeI_example.nb* for type I states
- *TypeII_example.nb* for type II states
- *TypeIII_example.nb* for type III states
- *TypeIV_example.nb* for type IV states

Run the corresponding Mathematica notebook. If it runs without errors, then this step is successful. The specifics of the inputs and outputs of the run may be read from the comments in the Mathematica notebook.

### Automation using Python
