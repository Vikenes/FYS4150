# FYS4150 Project 5: The Schrödinger equation

We present the work of Vetle A. Vikenes, Johan M. Kroken and Nanna Bryne in the final project of FYS4150.

The project report is found as `latex/project5.pdf`, or [here](https://github.com/Vikenes/FYS4150/blob/main/project5/latex/project5.pdf).

## **CODE**

All code is in C++ or Python and is found in the `code` directory. We have a `makefile` making it easy to reproduce our results.

#### **`C++`**
We have the following files:
* `main.cpp`
* `src/Box.cpp` (header `include/Box.hpp`)
* `src/Simulation.cpp` (header `include/Simulation.hpp`)

In addition, there is a source file `src/utils.cpp` with appropriate header file `include/utils.hpp`.

#### **`Python`**
We have the following scripts:
* `plot/analysis.py` 
* `plot/plot.py`

### How to run

It is important to move into the `code` directory for the following to work. One may write:
```
cd code
```
We have sat up five simulations that are ready to be run. The command
~~~
make all
~~~
defaults to running all of them, with no optimiser flag. To compile with `-O3`, type:
~~~
make all OPTIMISER=yes
~~~
To run one simulation at a time, pass the simulation code XX (see table in results-section) as argument as follows:
~~~
make all WHICH=XX
~~~
(`XX`$=$`ALL` is the same as not giving an argument)

>#### Example:
>Run the single-slit experiment with `-O3`-flag:
>~~~
>make all WHICH=SS OPTIMISER=yes
>~~~
>The resulting binary file (large!) will be saved as `../output/binfiles/SS_arma_cube.bin`.


To generate plots, write:
```
make plots
```
To generate animations, write:
```
make animations
```
The same syntax as above can be used to specify what figures or animations to produce.
>#### Examples:
>Create figures related to DS2:
>~~~
>make plots WHICH=DS2
>~~~
>Make animation of TS:
>~~~
>make animations WHICH=TS
>~~~

>#### Sidenote: 
>_The default setting makes the figures/movies appear simultaneously. To save as `.png`- and/or `.pdf`- and `.mp4`-files in the output folder, set the global variables_ `TEMP = True` _and/or_ `SAVE = True` _in_ `plot/plot.py`_. To show every plot in order, set_ `SHOW = True` _in the same file._


### Prerequisites

CHECK THESE!!

In order to run the C++ code the following libraries are needed:

* `sstream`
* `string`
* `iomanip`
* `vector`
* `fstream`
* `iostream`
* `cmath`
* `armadillo`
* `typeinfo`
* `random`
* `algorithm`
* `omp`

In order to run the python code the following libraries are needed:

* `numpy`
* `matplotlib.pyplot`
* `pandas`
* `os`
* `pyarma`


## **RESULTS**

The simulations are characterised by a code that will be the prefix in the names of the files resulting from it. The parameters that vary between them are presented below, however, this should make much more sense in context of the [method-section of our report](https://github.com/Vikenes/FYS4150/blob/main/project5/latex/project5.pdf), specifically tables I and II. In the following, XX is used to denote the code of some simulation.

|No.|Name               | Prefix    | Number of slits   | Simulation duration   | Vertical extent   |
|--:|:---               | :---      | :----:            | :----:                | :----:            |
|1. |No slits           | NS        | 0                 | 0.008                 | 0.05              |
|2. |Double-slit (1)    | DS1       | 2                 | 0.008                 | 0.10              |
|3. |Double-slit (2)    | DS2       | 2                 | 0.002                 | 0.20              |
|4. |Single-slit        | SS        | 1                 | 0.004                 | 0.20              |
|5. |Triple-slit        | TS        | 3                 | 0.004                 | 0.20              |

**Important**: Due to the large sizes of the binary files resulting from the simulations, none of them is present in the repository. However, by running the code as explained above, the raw data is saved as `output/binfiles/XX_arma_cube.bin`.

The figures we have produced are available in `output/plots/temp/` as `.png`-files and in `output/plots/pdf/` as `.pdf`-files.

The animations are found as `output/videos/XX_anim.mp4`. Appreciate the dynamic colour bar; the probability density value on this bar has not been altered with, with the limits deciding the colour scale change with each time step.