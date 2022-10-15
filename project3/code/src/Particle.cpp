#include "Particle.hpp"

Particle::Particle(double q, double m, arma::vec r, arma::vec v){
    q_ = q;
    m_ = m;
    r_ = r;
    v_ = v;
}


double Particle::charge(){
    return q_;
}

double Particle::mass(){
    return m_;
}

arma::vec Particle::position(){
    return r_;
}

arma::vec Particle::velocity(){
    return v_;
}

void Particle::new_position(arma::vec r){
    r_ = r;
}

void Particle::new_velocity(arma::vec v){
    v_ += v;
}