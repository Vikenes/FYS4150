#include "PenningTrap.hpp"
#include "utils.hpp"

//  Definition of constructor
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in){
    B0 = B0_in;
    V0 = V0_in;
    d = d_in;
}

std::vector<Particle> particles;

//  Add a particle to the the trap
void PenningTrap::add_particle(Particle p_in){
    particles.push_back(p_in);
}

// External electric field at point r=(x,y,z)
arma::vec external_E_field(arma::vec r){
    //  Should perhaps check if particles is emtpy or not
    arma::vec E_ext = r * Vdr;
    E_ext(2) *= -2;
    return E_ext;

    //  This is probably in the wrong place
    // int n = particles.size();
    // arma::vec E_ext = arma::vec(3).fill(0.);
    // for(int j=0;j<n;j++){
    //     //  Detailed
    //     arma::vec rj = particles[j].r();
    //     double qj = particles[j].q();
    //     double norm = arma::norm(r-rj);
    //     E_ext += qj * (r - rj) / std::pow(norm, 3);
    //     //  Memory efficient
    //     // E_ext += particles[j].q() * (r - particles[j].r()) / std::pow(arma::norm(r-particles[j].r()), 3);
    // }
    // return k_e * E_ext;
}


// External magnetic field at point r=(x,y,z)
arma::vec external_B_field(arma::vec r){
    arma::vec B_ext = arma::vec(3).fill(0.);
    B_ext(2) = B0;
    return B_ext;
  }


// Force on particle_i from particle_j
arma::vec force_particle(int i, int j){
    //  Detailed
    double qi = particles[i].q();
    double qj = particles[j].q();
    arma::vec ri = particles[i].r();
    arma::vec rj = particles[j].r();
    double norm = arma::norm(ri-rj);
    return k_e * qi * qj * (ri-rj) / std::pow(norm, 3);

    //  Memory effecient
    // return k_e * particles[i].q() * particles[j].q() * (particles[i].r() - particles[j].r()) / std::pow(arma::norm(particles[i].r() - particles[j].r()), 3);
  }

// The total force on particle_i from the external fields
arma::vec total_force_external(int i){
        //  F = qE + qv X B
        double qi = particles[i].q();
        arma::vec ri = particles[i].r();
        arma::vec vi = particles[i].v();
        return qi*external_E_field(ri) + qi*arma::cross(vi, external_B_field(ri));
    }

// The total force on particle_i from the other particles
arma::vec total_force_particles(int i){
    // Perhaps include some length tests or something. 
    arma::vec particle_force = arma::vec(3).fill(0.);
    for(int j=0;j<particles.size();j++){
        if(j!=i){
            particle_force += force_particle(i,j);
        }
    }
    return particle_force;
}

// The total force on particle_i from both external fields and other particles
arma::vec total_force(int i){
    return total_force_external(i) + total_force_particles(i);
}