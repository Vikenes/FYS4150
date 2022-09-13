# FYS4150 Project1 - Boundary value problem
The work of Vetle A. Vikenes, Nanna Bryne and Johan Mylius Kroken in FYS4150 - Computational Physics, fall 2022.

The answer to the problem is given in project1_FYS4150.pdf, and is found in the latex directory.

## How to run the code
```
# move into the code directory
g++ main.cpp src/* -I include -o main.exe

```

Inside src/algorithms.cpp you will find the algorithms for general Thomas, special Thomas, as well as functions for computing u(x) and f(x). 

Functions for converting floats to scientific string formats as well as writing to arrays to a txt file is included in src/utils.cpp

The headers for these two scripts is found in the include/ folder. 

## Reproducing the plots

To show the plots from the generated txt files, use run the script code/plot/plot.py. 

At the bottom of the script, there is a function called *show_plots*, which plots all txt files without saving them.  


