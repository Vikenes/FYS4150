#include "PenningTrap.hpp"
#include "Particle.hpp"
#include "utils.hpp"



PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in, bool interactions_in){
    B0 = B0_in;
    V0 = V0_in;
    d = d_in;
    Vd2r = V0/std::pow(d, 2);  
    interactions = interactions_in;

    filename = "untitled.txt";
    isready = false;
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
    B_ext.row(2).fill(B0);
    return B_ext;
}

arma::mat PenningTrap::compute_interaction_field(arma::mat R){
    for(int p=0; p<Np; p++){
        norm3.fill(0.);
        dist = R; dist.each_col() -= R.col(p);
        norm3.each_row() += arma::sqrt(arma::sum(arma::square(dist), 0));
        norm3.col(p).fill(1.);                                  //  avoid dividing by zero
        E_int.col(p) = k_e * arma::sum(Q%dist/norm3, 1);        //  sum over all particles
    }
    return E_int;
}

arma::mat PenningTrap::external_forces(arma::mat RU){
    compute_external_Efield(RU.rows(0,2));
    compute_external_Bfield(RU.rows(0,2));
   
    arma::mat U(3,Np); U.rows(0,2) = RU.rows(3,5);
    arma::mat UxB(3,Np);    //  cross products
    UxB.row(0) = U.row(1)%B_ext.row(2) - U.row(2)%B_ext.row(1);
    UxB.row(1) = U.row(2)%B_ext.row(0) - U.row(0)%B_ext.row(2);
    UxB.row(2) = U.row(0)%B_ext.row(1) - U.row(1)%B_ext.row(0);
    F_ext = Q % (E_ext + UxB);
    return F_ext;
}

arma::mat PenningTrap::internal_forces(arma::mat RU){
    compute_interaction_field(RU.rows(0,2));
    F_int = Q % E_int;
    return F_int;
}

void PenningTrap::ready(){

    if(Np==1){interactions = false;}

    Q = arma::mat(3, Np).fill(0.);     //  charges
    M = arma::mat(3, Np).fill(0.);     //  masses

    RU = arma::zeros(6, Np);      //   r,  v at a given time for all particles
    dRU = arma::zeros(6, Np);     //  dr, dv at a given time for all particles

    //  fields and forces:
    E_ext = arma::zeros(3, Np);
    B_ext = arma::zeros(3, Np);
    E_int = arma::zeros(3, Np);
    F_ext = arma::zeros(3, Np);
    F_int = arma::zeros(3, Np);

    //  for interactions:
    dist = arma::mat(3,Np);
    norm3 = arma::mat(3,Np);

    //  misc:
    K = arma::zeros(6, Np);

    isready = true;

}

void PenningTrap::simulate(double T, double dt, std::string scheme, bool point){
    int Nt = int(T/dt) + 1;                 //  number of time points 
    

    //  prepare if not already prepared:
    if(not isready){ready();}
    
    system = arma::cube(7, Np, Nt);
    //  initialise system:
    for(int p=0; p<Np; p++){
        for(int j=0; j<3; j++){
            Q.col(p).row(j) = particles.at(p) -> charge();
            M.col(p).row(j) = particles.at(p) -> mass();
        }
        system.slice(0).col(p).rows(1,3) = particles.at(p) -> position();
        system.slice(0).col(p).rows(4,6) = particles.at(p) -> velocity();
    }

    //  run simulation:
    for(int i=0; i<Nt-1; i++){
        RU = system.slice(i).rows(1,6);
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

        system.slice(i+1).rows(1,6) = RU + dRU;
        system.slice(i+1).row(0).fill(system(0,0,i) + dt);


        if(point){
            for(int p=0; p<Np; Np++){
                particles.at(p) -> new_position(system.slice(i).col(p).rows(1,3));
                particles.at(p) -> new_velocity(system.slice(i).col(p).rows(4,6));
            }
        }
    }


    for(int p=0; p<Np; p++){
            particles.at(p) -> new_position(system.slice(Nt-1).col(p).rows(1,3));
            particles.at(p) -> new_velocity(system.slice(Nt-1).col(p).rows(4,6));
        } 

    write_arma_to_file_scientific(system, filename);  
    std::cout << "\n    Written solution to ../output/data/" << filename << ".\n" << std::endl;
}

void PenningTrap::set_solution_filename(std::string filename_in){
    filename = filename_in;
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
    //arma::mat K = arma::zeros(6, Np);
    K.rows(0,2) = RU.rows(3,5);
    external_forces(RU);
    if(interactions){
        internal_forces(RU);
    }
    K.rows(3,5) = (F_ext + F_int) / M;
    return K;
}

