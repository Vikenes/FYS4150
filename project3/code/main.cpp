#include <algorithm>
#include "utils.hpp"
#include "algorithms.hpp"
#include "Particle.hpp"
#include "PenningTrap.hpp"

/**
 * PROJECT 3 FYS4150
 *      By: Johan Mylius Kroken, Vetle Vikenes and Nanna Bryne
 */

// double k_e = 1.38935333 * std::pow(10,5);
// double T = 9.64852558 * 10;
// double V = 9.64852558 * std::pow(10,7);
// double B0 = 1 * T;
// double V0 = 10 * V;
// double d = std::pow(10,4);
// double Vdr = 9.65;

void test_single_part(double x0, double y0, double z0, 
                     double vx0, double vy0, double vz0){

    arma::vec r0 = arma::vec({x0, y0, z0});
    arma::vec v0 = arma::vec({vx0, vy0, vz0});

    Particle test = Particle(1, 40, r0, v0);
    PenningTrap Trap = PenningTrap(B0, V0, d);

    double T = 10; // simulation duration 
    double dt = 0.01; // time step 

    Trap.add_particle(test);

    // Check if the particle 'test' is advanced by solver. It is not... 
    // std::cout << test.r() << std::endl;
    // std::cout << Trap.particles[0].r() << std::endl;

    Trap.simulate(T, dt, "FE");

    // std::cout << test.r() << std::endl;
    // std::cout << Trap.particles[0].r() << std::endl;


}

int main(){    
    Particle calcium = Particle(0, 20, arma::vec(3).randu(), arma::vec(3).randu());
    Particle calcium_ion = Particle(1,20, arma::vec(3).randu(), arma::vec(3).randu());



    // std::cout << "Neutral calcium:" << std::endl;
    // std::cout << "charge: " << calcium.q() << std::endl;
    // std::cout << "mass: " << calcium.m() << std::endl;
    // std::cout << "calcium position: " << calcium.r() << std::endl;
    // std::cout << "velocity: " << calcium.v() << std::endl;

    // std::cout << "Singly ionised calcium:" << std::endl;
    // std::cout << "charge: " << calcium_ion.q() << std::endl;
    // std::cout << "mass: " << calcium_ion.m() << std::endl;
    // std::cout << "calcium position: " << calcium_ion.r() << std::endl;
    // std::cout << "velocity: " << calcium_ion.v() << std::endl;

    test_single_part(0.1, 0, 1, 0, 0.1, 0);

    // //  Armadillo vector testing
    // arma::vec A = arma::vec(3).fill(2.);
    // arma::vec B = arma::vec(3).fill(3.);
    // arma::vec C = arma::vec("1.0, 2.0, 3.0");
    // std::cout << "Testing armadillo vectors" << std::endl;
    // std::cout << "A: " << A << std::endl;
    // std::cout << "B: " << B << std::endl;
    // std::cout << "C: " << C << std::endl;
    // std::cout << "A+B: " << A+B << std::endl;
    // std::cout << "C-(A+B)norm: " << arma::norm(C-(A+B)) << std::endl;


    return 0;
}