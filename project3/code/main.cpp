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
    Particle calcium(0, 20, arma::vec(3).randu(), arma::vec(3).randu());
    Particle calcium_ion(1,20, arma::vec(3).randu(), arma::vec(3).randu());

    std::cout << "Neutral calcium:" << std::endl;
    std::cout << "charge: " << calcium.q() << std::endl;
    std::cout << "mass: " << calcium.m() << std::endl;
    std::cout << "calcium position: " << calcium.r() << std::endl;
    std::cout << "velocity: " << calcium.v() << std::endl;

    std::cout << "Singly ionised calcium:" << std::endl;
    std::cout << "charge: " << calcium_ion.q() << std::endl;
    std::cout << "mass: " << calcium_ion.m() << std::endl;
    std::cout << "calcium position: " << calcium_ion.r() << std::endl;
    std::cout << "velocity: " << calcium_ion.v() << std::endl;
    return 0;
}