#ifndef __algorithms_hpp__
#define __algorithms_hpp__

#include <cmath>
#include <vector>

double u(const double x);

double f(const double x);



std::vector<double> generalThomas(std::vector<double> a, 
                                std::vector<double> b, 
                                std::vector<double> c, 
                                std::vector<double> g,
                                double v_min=0, double v_max=0);

std::vector<double> specialThomas(std::vector<double> g, double v_min=0,
                                double v_max=0, const int b_tilde0=2);

#endif 