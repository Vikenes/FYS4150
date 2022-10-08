#include "PenningTrap.hpp"

//  Definition of constructor
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in){
    B0 = B0_in;
    V0 = V0_in;
    d = d_in;
}

void PenningTrap::add_particle(Particle p_in){
    particles.push_bakc(p_in);
}