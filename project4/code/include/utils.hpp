#ifndef __utils_hpp__
#define __utils_hpp__

#include <sstream> 
#include <string>
#include <iomanip>
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <armadillo>
#include <typeinfo>

std::string scientific_format(const double d, const int width=15, const int prec=10);
std::string float_to_string(const double d, const int prec=2);

int write_to_file_scientific(std::vector<double> col1, std::vector<double> col2, std::string filename, int width=15, int prec=10);
int write_to_file_scientific(std::vector<double> col1, std::vector<double> col2, std::string filename, std::string header, int width=15, int prec=10);





#endif