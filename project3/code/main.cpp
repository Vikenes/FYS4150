#include <algorithm>
#include "utils.hpp"
#include "algorithms.hpp"
#include "Particle.hpp"
#include "PenningTrap.hpp"

/**
 * PROJECT 3 FYS4150
 *      By: Johan Mylius Kroken, Vetle Vikenes and Nanna Bryne
 */

double k_e = 1.38935333 * std::pow(10,5);
double T = 9.64852558 * 10;
double V = 9.64852558 * std::pow(10,7);
double B0 = 1 * T;
double V0 = 10 * V;
double d = std::pow(10,4);
double Vdr = 9.65;

int main(){
    arma::vec r = arma::vec(3, arma::fill::ones) * 2;
    arma::vec v = arma::vec(3, arma::fill::ones) * 4;
    
    Particle calcium(1, 20, r, v);

    std::cout << "Testing" << std::endl;
    std::cout << "Particle mass: " << calcium.m() << std::endl;
    return 0;
}