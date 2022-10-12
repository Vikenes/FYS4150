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

arma::vec Particle::r(){
    return r_;
}

arma::vec Particle::v(){
    return v_;
}

void Particle::superpose_position(arma::vec add_r){
    r_ += add_r;
}

void Particle::superpose_velocity(arma::vec add_v){
    v_ += add_v;
}