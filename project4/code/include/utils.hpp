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
// std::string scientific_format(double de, const int width, const int prec);


int write_to_file_scientific(std::vector<double> x, std::vector<double> v, std::string filename, int width=15, int prec=10);
int write_arma_to_file_scientific(arma::cube R, std::string filename, int width=15, int prec=10);

int write_to_file(std::vector<double> a, std::vector<double> b, std::string filename, int width=15);
int write_to_file(std::vector<double> a, std::vector<int> b, std::string filename, int width=15);





#endif