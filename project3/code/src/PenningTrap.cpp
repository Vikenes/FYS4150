#include "PenningTrap.hpp"
#include "Particle.hpp"
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
    N++;
}

// External electric field at point r=(x,y,z)
arma::vec PenningTrap::external_E_field(arma::vec r){
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
arma::vec PenningTrap::external_B_field(arma::vec r){
    arma::vec B_ext = arma::vec(3).fill(0.);
    B_ext(2) = B0;
    return B_ext;
  }


// Force on particle_i from particle_j
arma::vec PenningTrap::force_particle(int i, int j){
    //  Detailed
    double qi = particles[i].q_;
    double qj = particles[j].q_;
    arma::vec ri = particles[i].r_;
    arma::vec rj = particles[j].r_;
    double norm = arma::norm(ri-rj);
    return k_e * qi * qj * (ri-rj) / std::pow(norm, 3);

    //  Memory effecient
    // return k_e * particles[i].q() * particles[j].q() * (particles[i].r() - particles[j].r()) / std::pow(arma::norm(particles[i].r() - particles[j].r()), 3);
  }

// The total force on particle_i from the external fields
arma::vec PenningTrap::total_force_external(int i){
        //  F = qE + qv X B
        double qi = particles[i].q_;
        arma::vec ri = particles[i].r_;
        arma::vec vi = particles[i].v_;

        return qi*external_E_field(ri) + qi*arma::cross(vi, external_B_field(ri));
    }

// The total force on particle_i from the other particles
arma::vec PenningTrap::total_force_particles(int i){
    // Perhaps include some length tests or something. 
    arma::vec particle_force = arma::vec(3).fill(0.);
    for(int j=0; j<N; j++){
        if(j!=i){
            particle_force += force_particle(i,j);
        }
    }
    return particle_force;
}

// The total force on particle_i from both external fields and other particles
arma::vec PenningTrap::total_force(int i){
    return total_force_external(i) + total_force_particles(i);
}


void PenningTrap::simulate(double T, double dt, std::string method){
    // Temporary, performs forward euler on a single particle and writes to file 
    int Nt = int(T/dt) + 1; // Number of time steps 

    // Particle p1 = particles[0]; // Previous method, doesn't work... 

    std::vector<double> t (Nt, 0);

    arma::cube R = arma::cube(3, Nt, N).fill(0.);
    arma::cube U = arma::cube(3, Nt, N).fill(0.);

    for(int i=0; i<N; i++){
        // Initialize all particles
        R.slice(i).rows(0,2).col(0) = particles[i].r_;
        U.slice(i).rows(0,2).col(0) = particles[i].v_;
    }


    
    
    if(method=="Euler"){
        for(int i=1; i<Nt; i++){
    
            // evolve_forward_Euler(dt);
            particles[0].v_ += total_force(0) * dt / particles[0].m_;
            particles[0].r_ += particles[0].v_ * dt;  
            
            R.col(i) = particles[0].r_;
            U.col(i) = particles[0].v_;

            t[i] = t[i-1] + dt;
            }
        }



    if(method=="RK4"){
        arma::cube R_ = R; 
        arma::cube U_ = U;

        arma::mat K1r(3, N);
        arma::mat K2r(3, N); 
        arma::mat K3r(3, N); 
        arma::mat K4r(3, N); 

        arma::mat K1v(3, N); 
        arma::mat K2v(3, N); 
        arma::mat K3v(3, N); 
        arma::mat K4v(3, N); 

        for(int i=0; i<Nt-1; i++){
            K1r = U.col(i) * dt;
            K1v = total_force(0) / particles[0].m_ * dt;

            R_.col(i) += K1r/2;
            U_.col(i) += K1v/2;

            K2r = U_.col(i) * dt;
            K2v = total_force(0) / particles[0].m_ * dt;

            R_.col(i) += K2r/2;
            U_.col(i) += K2v/2;

            K3r = U_.col(i) * dt;
            K3v = total_force(0) / particles[0].m_ * dt;

            R_.col(i) += K3r/2;
            U_.col(i) += K3v/2;

            K4r = U_.col(i) * dt;
            K4v = total_force(0) / particles[0].m_ * dt;

            R_.col(i) += K4r;
            U_.col(i) += K4v;



            particles[0].r_ += (K1r + 2*K2r + 2*K3r + K4r)/6;
            particles[0].v_ += (K1v + 2*K2v + 2*K3v + K4v)/6;

            R.col(i+1) = particles[0].r_;
            U.col(i+1) = particles[0].v_;

            t[i+1] = t[i] + dt;

        }

    }

    std::string fname="test_z" + method;

    write_arma_to_file_scientific(R, t, fname);

}


// Evolve system one time step (dt) using Runge-Kutta 4th order 
void PenningTrap::evolve_RK4(double dt){
    
}

// Evolve system one time step (dt) using Forward Euler order
void PenningTrap::evolve_forward_Euler(double dt){

}