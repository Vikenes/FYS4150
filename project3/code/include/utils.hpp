#ifndef __utils_hpp__
#define __utils_hpp__

#include <sstream> 
#include <string>
#include <iomanip>
#include <vector>
#include <fstream>
#include <iostream>

std::string scientific_format(const double d, const int width, const int prec);

int write_to_file_scientific(std::vector<double> x, std::vector<double> v, std::string filename, int width=15, int prec=10);
int write_to_file(std::vector<double> a, std::vector<double> b, std::string filename, int width=15);


//  Constans for this project in units of micrometer, microseconds, atomic mass unit and elementary charge

double k_e; //  Coulomb constant
double T;   //  Tesla
double V;   //  Volt
double B0;  //  Magnetic field configuration in Penning Trap
double V0;  //  Electric potential configuration in Penning Trap
double d;   //  Size measure of Penning Trap
double Vdr;    //  Ratio of V0/d^2

#endif