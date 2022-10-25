# FYS4150 Project 3: 

We present the work of Vetle A. Vikenes, Johan M. Kroken and Nanna Bryne in FYS4150.

The answers to the problems are presented in `latex/project3.pdf`.

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
* `code/src/Particle.cpp`: sets up the class `Particle`, intances of which are initiated with a charge, mass, position and velocity. 
* `code/src/PenningTrap.cpp` 
* `code/src/utils.cpp`
* `code/src/algorithms.cpp` ????

Corresponding header-files with explanations of member functions etc. are found in `code/include/`. 