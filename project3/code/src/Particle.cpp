#include "Particle.hpp"

//  Definitions of constructor
Particle::Particle(double q, double m, arma::vec r, arma::vec v){
    q_ = q;
    m_ = m;
    r_ = r;
    v_ = v;
}

double Particle::q(){
    return q_;
}

double Particle::m(){
    return m_;
}
