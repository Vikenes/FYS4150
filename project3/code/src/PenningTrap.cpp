#include "PenningTrap.hpp"
#include "Particle.hpp"
#include "utils.hpp"



PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in, bool interactions_in){
    B0 = B0_in;
    V0 = V0_in;
    d = d_in;
    Vd2r = V0/std::pow(d, 2);  
    interactions = interactions_in;
}

void PenningTrap::add_particle(Particle &p_in){
    particles.push_back(&p_in);
    Np++;
}

arma::mat PenningTrap::compute_external_Efield(arma::mat R){
    E_ext = Vd2r * R; // GENERALISE!!
    E_ext.row(2) *= -2;
    return E_ext;
}

arma::mat PenningTrap::compute_external_Bfield(arma::mat R){
    B_ext(2) = B0;
    return B_ext;
}

arma::mat PenningTrap::compute_interaction_field(arma::mat R){
    arma::mat Rp;
    for(int p=0; p<Np; p++){
        //Rp = arma::ones(N, 3);
        //Rp.each_col() *= R.col(p);
        //arma::mat dist = R - Rp;
        //arma::vec norm = arma::norm(dist);
        //E_int.col(p) = k_e * arma::sum(Q%dist.each_row()/std::pow(norm, 3), 0);
        E_int.col(p) = arma::zeros(3);
    }
    return E_int;
}

arma::mat PenningTrap::external_forces(arma::mat RU){
    compute_external_Efield(RU.rows(0,2));
    compute_external_Bfield(RU.rows(0,2));
    F_ext = Q % (E_ext + arma::cross(RU.rows(3,5), B_ext));
    return F_ext;
}

arma::mat PenningTrap::internal_forces(arma::mat RU){
    compute_interaction_field(RU.rows(0,2));
    F_int = Q*E_int;
    return F_int
}


void PenningTrap::simulate(double T, double dt, std::string scheme, bool point){
    int Nt = int(T/dt) + 1;                 //  number of time points 
    arma::vec t = arma::zeros(Nt);          //  vector for time points


    arma::cube system = arma::cube(6, Np, Nt);
    /**
     * 'system' shall contain information about the position and velocity of each particle at any time
     *      slices      -> time steps
     *      cols        -> particles
     *      rows(0,2)   -> position (x,y,z)
     *      rows(3,5)   -> velocity (vx,vy,vz)
     */

    Q = arma::vec(Np);     //  charges
    M = arma::vec(Np);     //  masses

    //  initialise system
    for(int p=0; p<Np; p++){
        Q(p) = particles.at(p) -> charge();
        M(p) = particles.at(p) -> mass();
        system.slice(0).col(p).rows(0,2) = particles.at(p) -> position();
        system.slice(0).col(p).rows(3,5) = particles.at(p) -> velocity();
    }

    RU = arma::zeros(6, Np);      //   r,  v at a given time for all particles
    dRU = arma::zeros(6, Np);     //  dr, dv at a given time for all particles

    E_ext = arma::zeros(6, Np);
    B_ext = arma::zeros(6, Np);
    E_int = arma::zeros(6, Np);
    F_ext = arma::zeros(6, Np);
    F_int = arma::zeros(6, Np);

    for(int i=0; i<Nt-1; i++){
        RU = system.slice(i);
        if(scheme=="RK4"){
            dRU = evolve_RK4(dt, RU);
        }
        else if(scheme=="FE"){
            dRU = evolve_FE(dt, RU);
        }
        else{
            std::cout << "Arguments FE and RK4 are the only valid ones." << std::endl;
            assert(false);
        }

        system.slice(i) = RU + dRU;

        if(point){
            for(int p=0; p<Np; Np++){
                particles.at(p) -> position(system.slice(i).col(p).rows(0,2));
                particles.at(p) -> velocity(system.slice(i).col(p).rows(3,5));
            }
        }
        t[i+1] = t[i]+dt;
    }

    for(int p=0; p<Np; Np++){
            particles.at(p) -> position(system.slice(Nt-1).col(p).rows(0,2));
            particles.at(p) -> velocity(system.slice(Nt-1).col(p).rows(3,5));
        }

}


arma::mat PenningTrap::evolve_FE(double dt, arma::mat RU){
    external_forces(RU);
    if(interactions){
        internal_forces(RU);
    }
    dRU.rows(3,5) = (F_ext + F_int) / M * dt;
    dRU.rows(0,2) = (RU.rows(3,5) + dRU.rows(3,5)) * dt;
    return dRU;
}

arma::mat PenningTrap::evolve_RK4(double dt, arma::mat RU){
    K1 = K_val(RU) * dt;
    K2 = K_val(RU+K1/2) * dt;
    K3 = K_val(RU+K2/2) * dt;
    K4 = K_val(RU+K3) * dt; 

    dRU = (K1+2*K2+2*K3+K4)/6;
    return dRU;
}

arma::mat PenningTrap::K_val(arma::mat RU){
    arma::mat K = arma::zeros(6, Np);
    K.rows(0,2) = RU.rows(3,5);
    external_forces(RU);
    if(interactions){
        internal_forces(RU);
    }
    K.rows(3,5) = (F_ext + F_int) / M;
    return K;
}

