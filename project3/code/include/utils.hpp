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
int write_arma_to_file_scientific(arma::cube R, arma::vec t, std::string filename, int width=15, int prec=10);

int write_to_file(std::vector<double> a, std::vector<double> b, std::string filename, int width=15);


//  Constans for this project in units of micrometer, microseconds, atomic mass unit and elementary charge

extern double k_e; //  Coulomb constant
extern double T;   //  Tesla
extern double V;   //  Volt
extern double B0;  //  Magnetic field configuration in Penning Trap
extern double V0;  //  Electric potential configuration in Penning Trap
extern double d;   //  Size measure of Penning Trap
extern double Vdr;    //  Ratio of V0/d^2

// extern double k_e = 1.38935333 * std::pow(10,5);
// extern double T = 9.64852558 * 10;
// extern double V = 9.64852558 * std::pow(10,7);
// extern double B0 = 1 * T;
// extern double V0 = 10 * V;
// extern double d = std::pow(10,4);
// extern double Vdr = 9.65;



#endif