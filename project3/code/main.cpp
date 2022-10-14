#include <algorithm>
#include "utils.hpp"
#include "algorithms.hpp"
#include "Particle.hpp"
#include "PenningTrap.hpp"

/**
 * PROJECT 3 FYS4150
 *      By: Johan Mylius Kroken, Vetle Vikenes and Nanna Bryne
 */


void test_single_part(double T, double dt, std::string method, bool interactions=true){

    // Use provided initial conditions 
    arma::vec r0 = arma::vec({20, 0, 20});
    arma::vec v0 = arma::vec({0, 25, 0});

    Particle test = Particle(1, 40, r0, v0);
    PenningTrap Trap = PenningTrap(B0, V0, d, interactions);

    Trap.add_particle(test);
    Trap.simulate(T, dt, method);

}

void test_two_particles(double T, double dt, std::string method, bool interactions=true){

    PenningTrap Trap = PenningTrap(B0,V0,d,interactions);

    Particle p1 = Particle(1, 40, arma::vec({20,0,20}), arma::vec({0,25,0}));
    Particle p2 = Particle(1, 40, arma::vec({25,25,0}), arma::vec({0,40,5}));
    Trap.add_particle(p1);
    Trap.add_particle(p2);
    Trap.simulate(T, dt, method);

}


int main(){    
    Particle calcium = Particle(0, 20, arma::vec(3).randu(), arma::vec(3).randu());
    Particle calcium_ion = Particle(1,20, arma::vec(3).randu(), arma::vec(3).randu());

    //test_single_part(50, 0.1, "Euler");
    //test_single_part(50, 0.1, "RK4");
    test_two_particles(50, 0.01, "FE", false);
    test_two_particles(50, 0.01, "RK4", false);



    return 0;
}