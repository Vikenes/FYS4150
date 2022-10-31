# FYS4150 Project 4: The 2D Ising Model

We present the work of Vetle A. Vikenes, Johan M. Kroken and Nanna Bryne in FYS4150.

The project report is found as `latex/project4.pdf`.

## The code

All code is in C++ or Python and is found in the `code` directory.

### How to run
It is important to move into the `code` directory for the following to work. One may write:
```
cd code
```

To run all code in `main.cpp`, write:
```
make all
```
Swap `all` with `all2` or `all3` to compile with `-O2` or `-O3` flag.


To genereate plots, write:
```
make plots
```
This command runs the script `code/plot/analysis.py` which uses code from `code/plot/plot.py`
### Content

The script that performs all computations in this project is `code/main.cpp` which includes:
* `code/src/Particle.cpp`: sets up the class `Particle`, intances of which are initiated with a charge, mass, position and velocity
* `code/src/PenningTrap.cpp`: sets up the class `PenningTrap`, intances of which are initiated with a magnetic field strength, electric potential scale and characteristic dimesion
* `code/src/utils.cpp`: contains functions for writing `.txt`-files, as well as physical constants and other parameters we do not intend to vary in this project

Corresponding **header-files with user-explanations** of member functions etc. are found in `code/include/`. 

## Results

You will find the final results as `.txt`-files in `output/data/RK4` and `output/data/FE`, and presented as figures in `.pdf`-files in `output/plots/`.
