//
// A collection of useful functions
//


// First we add an "include guard". It ensures that this header file can 
// only end up being included *once* for the compilation of any given .cpp file
#ifndef __utils_hpp__  
#define __utils_hpp__

// Now we include headers we need
#include <armadillo>
#include <sstream>
#include <string>
#include <iomanip>


// Below we give some function *declarations*.
// The function *definitions* (the actual code) 
// lives in src/utils.cpp

// Return a string with a double in scientific notation
std::string scientific_format(double d, const int width=20, const int prec=10);

// Return a string with an armadillo vector in scientific notation
std::string scientific_format(const arma::vec& v, const int width=20, const int prec=10);

#endif  // end of include guard __utils_hpp__