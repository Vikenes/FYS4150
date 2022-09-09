#ifndef __utils_hpp__
#define __utils_hpp__

#include <sstream> 
#include <string>
#include <iomanip>
#include <vector>
#include <fstream>
#include <iostream>

std::string scientific_format(const double d, const int width, const int prec);

int write_to_file(std::vector<double> x, 
                    std::vector<double> v, 
                    std::string filename,
                    int width=15, int prec=10);

#endif