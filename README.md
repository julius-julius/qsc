# qsc
A numerical non-perturbative quantum spectral curve solver and database of non-perturbative values of the scaling dimensions of states/operators of planar ${\cal N} = 4$ supersymmetric Yang-Mills theory (SYM), in a wide range of the 't Hooft coupling $g\equiv \sqrt{\lambda}/4\pi \equiv g\_{\tt YM}^2 N /4\pi$, where $g_{\tt YM}$ is the Yang-Mills coupling and $N-1$ is the rank of the gauge group ${\rm SU}(N)$.

## Structure of Repository

This repository has 6 directories: */core*, */auxiliary*, */run*, */data*, */seed* and */demo*. 

### */core*

This directory contains the C++ code, for the various types of local opearators in planar ${\cal N =}$ 4 SYM.

### */auxiliary*

This directory contains Mathematica and Python packages/modules/scripts which assist in automation and parameter management.

### */run*

This directory contains *.ipnb* files, which is the main interface through which one is able to run the numerics.

### */data*

This directory contains 3 sub-directories: */numerical*, */perturbative* and */output*

***/numerical***

This sub-directory contains precomputed scaling dimensions for all 219 states in ${\cal N =}$ 4 SYM with bare dimension $\Delta_0 \leq 6$, in a wide range of the 't Hooft coupling. This data is readily available for anyone to use in their research. The spectral data for a state with 

$\texttt{State ID}:    {}\_{\Delta\_0}[n\_{{\bf b}\_1}\\;n\_{{\bf b}\_2}\\;n\_{{\bf f}\_1}\\;n\_{{\bf f}\_2}\\;n\_{{\bf f}\_3}\\;n\_{{\bf f}\_4}\\;n\_{{\bf a}\_1}\\;n\_{{\bf a}\_2}]\_{\tt sol}$ 

is stored in the *.mx* file called

*numerical\_data\_Delta0*$\\,\Delta_0\\,$*\_b1*$\\,n\_{{\bf b}\_1\\,}$*\_b2*$\\,n\_{{\bf b}\_2}\\,$*\_f1*$\\,n\_{{\bf f}\_1}\\,$*\_f2*$\\,n\_{{\bf f}\_2}\\,$*\_f3*$\\,n\_{{\bf f}\_3}\\,$*\_f4*$\\,n\_{{\bf f}\_4\\,}$*\_a1*$\\,n\_{{\bf a}\_1}\\,$*\_a2*$\\,n\_{{\bf a}\_2}\\,$*\_sol*$\\,{\tt sol}\\,$*.mx*.

The file contains a $2\times 2$ array `dataGH`. Each element of of this array is a tuple $(g,\Delta)$ where $g$ is the value of the 't Hooft coupling, and $\Delta$ is the value of the dimension of the state at that value of $g$. The precision of $\Delta$ is **... (J: to be added)**

Depending on the type of the state/operator, we are able to provide data for different ranges in $g$.

***/perturbative***

This sub-directory contains precomputed perturbative data to initialise the numerical algorithm for all 219 states in ${\cal N =}$ 4 SYM with bare dimension $\Delta_0 \leq 6$. The perturvative data for a state with 

$\texttt{State ID}:    {}\_{\Delta\_0}[n\_{{\bf b}\_1}\\;n\_{{\bf b}\_2}\\;n\_{{\bf f}\_1}\\;n\_{{\bf f}\_2}\\;n\_{{\bf f}\_3}\\;n\_{{\bf f}\_4}\\;n\_{{\bf a}\_1}\\;n\_{{\bf a}\_2}]\_{\tt sol}$ 

is stored in the *.mx* file called

*perturbative\_data\_Delta0*$\\,\Delta_0\\,$*\_b1*$\\,n\_{{\bf b}\_1\\,}$*\_b2*$\\,n\_{{\bf b}\_2}\\,$*\_f1*$\\,n\_{{\bf f}\_1}\\,$*\_f2*$\\,n\_{{\bf f}\_2}\\,$*\_f3*$\\,n\_{{\bf f}\_3}\\,$*\_f4*$\\,n\_{{\bf f}\_4\\,}$*\_a1*$\\,n\_{{\bf a}\_1}\\,$*\_a2*$\\,n\_{{\bf a}\_2}\\,$*\_sol*$\\,{\tt sol}\\,$*.mx*.

The file contains the substitution rule `sbWeak` which substitutes the perturbative expansion of $\Delta$ as well as the expansion coefficeints $c_{a,n}$ of ${\bf P}\_a$ and $c^{a,n}$ of ${\bf P}^a$.

***/output***

This is the default directory where outputs of the numerical runs get stored. The default output format of a state with 

$\texttt{State ID}:    {}\_{\Delta\_0}[n\_{{\bf b}\_1}\\;n\_{{\bf b}\_2}\\;n\_{{\bf f}\_1}\\;n\_{{\bf f}\_2}\\;n\_{{\bf f}\_3}\\;n\_{{\bf f}\_4}\\;n\_{{\bf a}\_1}\\;n\_{{\bf a}\_2}]\_{\tt sol}$ 

is an *.mx* file called

*numerical\_data\_Delta0*$\\,\Delta_0\\,$*\_b1*$\\,n\_{{\bf b}\_1\\,}$*\_b2*$\\,n\_{{\bf b}\_2}\\,$*\_f1*$\\,n\_{{\bf f}\_1}\\,$*\_f2*$\\,n\_{{\bf f}\_2}\\,$*\_f3*$\\,n\_{{\bf f}\_3}\\,$*\_f4*$\\,n\_{{\bf f}\_4\\,}$*\_a1*$\\,n\_{{\bf a}\_1}\\,$*\_a2*$\\,n\_{{\bf a}\_2}\\,$*\_sol*$\\,{\tt sol}\\,$*.mx*.

### */seed*

This directory contains a Mathematica notebook **(J: to be added)**

### */demo*

This directory contains a Mathematica notebook **(J: to be added)**

## Installation and Compilation

The first step is to install the requisite C++ compilers and packages and compile the C++ code in the */core* directory.

### Compilers and Packages

The `g++` compiler, and the `libcln` and `libcln-dev` packages are necessary to be installed. Other compliers may be used depending on taste and packages may be installed using a package manager. 

**Debian Linux**

Below we present the steps to install the above on a Debian linux system. The following commands should be run on the terminal 

`> sudo apt install g++`

`> sudo apt install libcln6 `

`> sudo apt install libcln-dev`

**Mac**

Below we present the steps to install the above on a Mac system. First, you need to install a package manager such as Homebrew (https://brew.sh/). In the sequel we assume that Homebrew is installed. In the terminal, run

`> brew install g++`

`> brew install cln `

**Windows**

In order to use the packages on a Windows operating system, one needs to install a Windows Subsystem for Linux (WSL), though one must bare in mind to type `wsl` before any Shell command. 


### Compilation

In the */core* directory, there are four *.cpp* files. Each code pertains to a particular type of state in planar ${\cal N} = 4$ SYM:

- For type I states, use *TypeI_core.cpp*
- For type II states, use *TypeII_core.cpp*
- For type III states, use *TypeIII_core.cpp*
- For type IV states, use *TypeIV_core.cpp*

To compile any of the above, use the following command where *source.cpp* is the C++ file, and *executable.out* is the executable output file

`> g++ source.cpp -lm -lcln -o executable.out`

In order for the executables to be compatible with the Mathematica notebooks in the */core* directory, the names of the execulatbles in the four cases should be as follows.

- For type I states, use *TypeI_exec.out*
- For type II states, use *TypeII_exec.out*
- For type III states, use *TypeIII_exec.out*
- For type IV states, use *TypeIV_exec.out*

Once compiled (it should compile without errors, warnings are OK), ensure that the exeutable file and the Mathematica notebook are in the same directory. The C++ sourcecode needs to be compiled only once. 

## Execution and Automation

In order to ensure that everything goes smoothly, and the code was compiled properly, we recommend to test the executables by running the corresponding Mathematica notebooks in */core*. 

### Execute from Mathematica notebook

Ensure that the executable file is in the same directry as the Mathematica notebooks available in */core*. Depending on your application, you should use

- *TypeI_example.nb* for type I states
- *TypeII_example.nb* for type II states
- *TypeIII_example.nb* for type III states
- *TypeIV_example.nb* for type IV states

Run the corresponding Mathematica notebook. If it runs without errors, then this step is successful. The specifics of the inputs and outputs of the run may be read from the comments in the Mathematica notebook.

### Automation using Python

The directory */run* contains four *.ipynb* files, which run and automatically manage the hyperparameters of a particular state:

- *TypeI_run.ipynb* for type I states
- *TypeII_run.ipynb* for type II states
- *TypeIII_run.ipynb* for type III states
- *TypeIV_run.ipynb* for type IV states

Select the *.ipynb* file according to the type of state which you wish to run. 

**Pre-requisites**

Before running the *.ipynb* ensure that all the paths are correctly adjusted based on your specific project organisation. The default path settings will work for a project that has the same strucutre as this repository, and should work without change if you pull this repository to your local system. There are various places where the path needs to be adjusted. We list them below:

- In */run/TypeI_run.ipynb*, you need to specify the path to *TypeI_module.ipynb*. By default *TypeI_module.ipynb* is located in */auxiliary* (same for *.ipynb* files of other types of states).
- In */auxiliary/TypeI_module.ipynb*, you need to specify the path to *TypeI_run.wls*. By default *TypeI_script.wls* is located in */auxiliary* (same for *.ipynb* files of other types of states).
- In */auxiliary/TypeI_script.wls* you need to specify the path to *TypeI_exec.out*, *TypeI_package.wl* and the location of your output data. By default *TypeI_exec.out* is located in */core*, *TypeI_package.wl* is located in */auxiliary*, and the output location is */data/output* (same for *.ipynb* files of other types of states).
- In */auxiliary/TypeI_package.wl* you need to specify the location of the perturbative data. By default the perturbative data is located in */data/perturbative* (same for *.ipynb* files of other types of states).

To initialise a state with a given $\texttt{State ID}$ you must ensure that there is a *.mx* file with perturbative data that you will need to initialise the numerics. For the 219 states with bare dimenion $\Delta_0\leq 6$, this is already precomputed, and such a file is available in */data/perturbative*. For states with $\Delta_0>6$, you can generate this data by running **(J: to be added)**. By default, the perturbative data generated will be stored in */data/perturbative*.

**Running**

In order to run a state with a given $\texttt{State ID}$, open the *.ipynb* associted with the type of state as a Jupyter notebook (https://jupyter.org/). Then specify the $\texttt{State ID}$ of the state which you want to run, follow the comments in the noteboook and run. If all the paths are specified properly, and the perturbative data exists, then it should run smoothly, and automatically start to produce spectral data for that state.

**(Mention Anaconda, Homebrew and Mathematica installation)**
