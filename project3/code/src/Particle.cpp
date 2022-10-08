#include "Particle.hpp"

//  Definitions of constructor
Particle::Particle(double q, double m, arma::vec r, arma::vec, v){
    q_ = q;
    m_ = m;
    r_ = r;
    v_ = v;
}