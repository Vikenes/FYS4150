# FYS4150 Project 2: Buckling beam

We present the work of Vetle A. Vikenes, Johan M. Kroken and Nanna Bryne in FYS4150.

The answers to the problems are presented in `latex/project2.pdf`.

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
To genereate plots, write:
```
make plots
```

### Content

The script that performs all computations in this project is `code/main.cpp` which includes:

* `code/src/algorithms.cpp` for creating triagonal matrices and solving eigenvalue problem using Jacobi rotation method, whose header is found in `code/include`
* `code/src/utils.cpp` with functions for writing data to file, whose header is found in `code/include`

The code for plotting is found in `code/plot/plot.py`.