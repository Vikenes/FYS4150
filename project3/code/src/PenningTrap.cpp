#include "PenningTrap.hpp"
#include "Particle.hpp"
#include "utils.hpp"

//  Definition of constructor
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in, bool interactions_in){
    B0 = B0_in;
    V0 = V0_in;
    d = d_in;
    interactions = interactions_in;
}

std::vector<Particle*> particles;

//  Add a particle to the the trap
void PenningTrap::add_particle(Particle &p_in){
    particles.push_back(&p_in);
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
    double qi = particles.at(i)->q();
    double qj = particles.at(j)->q();
    arma::vec ri = particles.at(i)->r();
    arma::vec rj = particles.at(j)->r();
    double norm = arma::norm(ri-rj);
    return k_e * qi * qj * (ri-rj) / std::pow(norm, 3);

    //  Memory effecient
    // return k_e * particles[i].q() * particles[j].q() * (particles[i].r() - particles[j].r()) / std::pow(arma::norm(particles[i].r() - particles[j].r()), 3);
  }

// The total force on particle_i from the external fields
arma::vec PenningTrap::total_force_external(int i, arma::vec ri){
        //  F = qE + qv X B
        double qi = particles.at(i)->q();

        return qi*external_E_field(ri.rows(0,2)) + qi*arma::cross(ri.rows(3,5), external_B_field(ri.rows(0,2)));
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
arma::vec PenningTrap::total_force(int i, arma::vec ri){
    if(interactions){
        return total_force_external(i, ri) + total_force_particles(i);
    }
    else{
        return total_force_external(i, ri);
    }
}


arma::vec PenningTrap::RK4_K_val(int i, arma::vec ri){
    arma::vec K(6);
    K.rows(0,2) = ri.rows(3,5);
    K.rows(3,5) = total_force(i, ri) / particles.at(i)->m();
    return K;
}

void PenningTrap::simulate(double T, double dt, std::string method){
    // Temporary, performs forward euler on a single particle and writes to file 

    int Nt = int(T/dt) + 1; // Number of time points 
    std::vector<double> t (Nt, 0);

    /* R cube:
     * rows (0,2): particle position x,y,z
     * rows (3,5): particle velocities vx,vy,vz
     * Nt columns: phase coordinates at each timestep
     * N slices: One slice corresponds to one particle     
    */ 

    arma::cube R = arma::cube(6, Nt, N).fill(0.);


    for(int i=0; i<N; i++){
        // Initialize all particles
        R.slice(i).rows(0,2).col(0) = particles.at(i)->r();
        R.slice(i).rows(3,5).col(0) = particles.at(i)->v();
    }


    
    
    if(method=="Euler"){

        for(int k=0; k<Nt-1; k++){
            // time loop
            arma::cube dR(6, 1, N);
            for(int i=0; i<N; i++){
                // particle loop
                
                dR.slice(i).rows(3,5) = total_force(i, R.slice(i).col(k)) * dt / particles.at(i) -> m(); 
                particles.at(i) -> superpose_velocity(dR.slice(i).rows(3,5));

                dR.slice(i).rows(0,2) = particles.at(i) -> v() * dt;
                particles.at(i) -> superpose_position(dR.slice(i).rows(0,2));
            }
            R.col(k+1) = R.col(k) + dR;
            
            t[k+1] = t[k] + dt;
            }
        }



    if(method=="RK4"){
        
        for(int k=0; k<Nt-1; k++){
            // time loop
            arma::cube dR(6, 1, N);

            for(int i=0; i<N; i++){
                // find dR for all particles
                arma::vec R_ik = R.slice(i).col(k);

                arma::vec K1(6);
                arma::vec K2(6); 
                arma::vec K3(6); 
                arma::vec K4(6); 

                K1 = RK4_K_val(i, R_ik) * dt;
                K2 = RK4_K_val(i, R_ik + K1/2) * dt;
                K3 = RK4_K_val(i, R_ik + K2/2) * dt;
                K4 = RK4_K_val(i, R_ik + K3) * dt;

                dR.slice(i).col(0) = (K1 + 2*K2 + 2*K3 + K4)/6;
            }
            // update position and velocity vectors
            R.col(k+1) = R.col(k) + dR;

            for(int i=0; i<N; i++){
                // update positions and velocities in particle object
                // important for the interaction force
                particles.at(i) -> superpose_position(dR.slice(i).rows(0,2));
                particles.at(i) -> superpose_velocity(dR.slice(i).rows(3,5));
            }
            t[k+1] = t[k] + dt;

        }

    }

    std::string fname = method + "_N" + std::to_string(N);

    write_arma_to_file_scientific(R, t, fname);
}


// Evolve system one time step (dt) using Runge-Kutta 4th order 
void PenningTrap::evolve_RK4(double dt){
    
}

// Evolve system one time step (dt) using Forward Euler order
void PenningTrap::evolve_forward_Euler(double dt){

}