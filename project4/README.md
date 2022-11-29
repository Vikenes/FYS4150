# FYS4150 Project 4: The 2D Ising Model

We present the work of Vetle A. Vikenes, Johan M. Kroken and Nanna Bryne in FYS4150.

The project report is found as `latex/project4.pdf`, or [here](https://github.com/Vikenes/FYS4150/blob/main/project4/latex/project4.pdf)

## The code

All code is in C++ or Python and is found in the `code` directory.

### Dependencies
In order to run the C++ code the following libraries are needed:
```
sstream
string
iomanip
vector
fstream
iostream
cmath
armadillo
typeinfo
random
algorithm

```

### How to run
It is important to move into the `code` directory for the following to work. One may write:
```
cd code
```

The code is divided into specific tasks:

#### Analytic comparison:
Compile with:
```
make compile_anal
```
Run with:
```
make run_anal Nlog0=arg1 Nlog1=arg2 T=arg3
```
with arguments: 

* arg1: Starting point of MC-runs, $N=10^\mathrm{arg1}$
* arg2: End point of MC-runs, $N=10^\mathrm{arg2}$
* arg3: Temperature of the simulation.

#### Equilibration
Compile with
```
make compile_equi
```
Run with:
```
make run_equi N=arg1, T=arg2, order=arg3
```
with arguments:

* arg1: Number of MC-runs, $N=\mathrm{arg1}$
* arg2: Temperature $T=\mathrm{arg2}$
* arg3: Initialization: arg3=0 (random) or arg3=1 (ordered): $\mathrm{order}=\mathrm{arg3}$

#### Probability distribution
Compile with
```
make compile_pdf
```
Run with:
```
make run_pdf NMC=arg1, Neq=arg2, T=arg3 order=arg4
```
with arguments:

* arg1: Number of MC-runs, $N_\mathrm{MC}=\mathrm{arg1}$
* arg2: Equilibration cycles, $N_\mathrm{eq}=\mathrm{arg2}$
* arg3: Temperature, $T=\mathrm{arg3}$
* arg4: Initialization: arg4=0 (random) or arg4=1 (ordered): $\mathrm{order}=\mathrm{arg4}$

#### Serial and parallel runs
These are compiled similarly:
```
make parallel_compile
make serial_compile
```
and run:
```
make parallel_run L=arg1 T0=arg2 T1=arg3 nT=arg4 NMC=arg5 Neq=arg6
make serial_run L=arg1 T0=arg2 T1=arg3 nT=arg4 NMC=arg5 Neq=arg6
```
with arguments:

* arg1: Lattice size: $L=\mathrm{arg1}$
* arg2: Lower temperature limit ($T\in[T_0,T_1]$): $T_0 = \mathrm{arg2}$
* arg3: Upper temperature limit: $T_1 = \mathrm{arg3}$
* arg4: Number of temperatures to be evaluated: $\mathrm{nT} = \mathrm{arg4}$
* arg5: Number of MC-runs, $N_\mathrm{MC}=\mathrm{arg5}$
* arg6: Equilibration cycles, $N_\mathrm{eq}=\mathrm{arg6}$

####  Timed runs
These are compiled similarly:
```
make time_parallel_compile
make time_serial_compile
```
and run:
```
make timing_parallel nL=arg1 T0=arg2 T1=arg3 nT=arg4 NMC=arg5 Neq=arg6
make timing_serial nL=arg1 T0=arg2 T1=arg3 nT=arg4 NMC=arg5 Neq=arg6
```
with arguments:

* arg1: Number of lattice sizes to evaluate: $n_L=\mathrm{arg1}$
* arg2: Lower temperature limit ($T\in[T_0,T_1]$): $T_0 = \mathrm{arg2}$
* arg3: Upper temperature limit: $T_1 = \mathrm{arg3}$
* arg4: Number of temperatures to be evaluated: $\mathrm{nT} = \mathrm{arg4}$
* arg5: Number of MC-runs, $N_\mathrm{MC}=\mathrm{arg5}$
* arg6: Equilibration cycles, $N_\mathrm{eq}=\mathrm{arg6}$

#### Plots

To generate plots, write:
```
make plots
```
This command runs the script `code/plot/analysis.py`


