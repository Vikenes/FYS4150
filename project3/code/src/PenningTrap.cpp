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

void PenningTrap::add_particle(Particle &particle){
    particles.push_back(&particle);
    Np++;
}

void PenningTrap::print_particles(){
    std::cout << "\nParticles at time t = " << t << std::endl;
    std::cout << "_____________________________________________________" << std::endl;
    for(int p=0; p<Np; p++){
        arma::vec r = particles.at(p) -> position();
        arma::vec v = particles.at(p) -> velocity();
        std::cout << "r = (" << r(0) << ", " << r(1) << ", " << r(2) << ")    " << "|r| = " << arma::norm(r) << std::endl;
        std::cout << "v = (" << v(0) << ", " << v(1) << ", " << v(2) << ")    " << "|v| = " << arma::norm(v) << std::endl;
        std::cout << "-----------------------------------------------------" << std::endl;
    }
}

void PenningTrap::generate_random_identical_particles(double charge, double mass, int no_of_particles, int seed){
    arma::arma_rng::set_seed(seed);
    arma::vec rr, vv;

    for(int p=0; p<no_of_particles; p++){
        rr = arma::vec(3).randn() * 0.1 * d;    //  random initial position
        vv = arma::vec(3).randn() * 0.1 * d;    //  random initial velocity 
        generate_particle(particles, rr, vv);
    }
}

void PenningTrap::apply_time_dependence(double amplitude, double frequency){
    f = amplitude;
    omega_V = frequency;
    time_dep = true;
}

void PenningTrap::switch_interactions(){
    if(interactions){interactions = false;}
    else{interactions = true;}
}

void PenningTrap::switch_interactions(std::string switch_in){
    if(switch_in=="off" or switch_in=="OFF"){interactions = false;}
    else if(switch_in=="on" or switch_in=="ON"){interactions = true;}
    else{
        std::cout << "Arguments ON/on and OFF/off are the only valid ones." << std::endl;      
        assert(false);}
}

int PenningTrap::count_particles(){
    int Np_trapped = 0;
    for(int p=0; p<Np; p++){
        arma::vec r = particles.at(p) -> position();
        if(arma::norm(r)<d){Np_trapped++;}
    }
    return Np_trapped;
}

void PenningTrap::simulate(double T, double dt, std::string scheme, bool point){
    int Nt = int(T/dt) + 1;         //  number of time points 
    

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
        t = system(0,0,i);
        if(scheme=="RK4"){ // FIXME (want less if/else)
            if(time_dep){dRU = evolve_RK4(dt, t, RU);}
            else{dRU = evolve_RK4(dt, RU);}
        }
        else if(scheme=="FE"){
            if(time_dep){dRU = evolve_FE(dt, t, RU);}
            else{dRU = evolve_FE(dt, RU);}
        }
        else{
            std::cout << "Arguments FE and RK4 are the only valid ones." << std::endl;
            assert(false);
        }

        system.slice(i+1).rows(1,6) = RU + dRU;
        system.slice(i+1).row(0).fill(t + dt);


        if(point){ // fixme
            for(int p=0; p<Np; Np++){
                particles.at(p) -> new_position(system.slice(i).col(p).rows(1,3));
                particles.at(p) -> new_velocity(system.slice(i).col(p).rows(4,6));
            }
        }
    }
    t = system(0, 0, Nt-1);
    for(int p=0; p<Np; p++){
            particles.at(p) -> new_position(system.slice(Nt-1).col(p).rows(1,3));
            particles.at(p) -> new_velocity(system.slice(Nt-1).col(p).rows(4,6));
        }

}

void PenningTrap::save_solution(std::string filename){
    write_arma_to_file_scientific(system, filename);  
    std::cout << "\n    Written solution to ../output/data/" << filename << ".txt.\n" << std::endl;
}

void PenningTrap::ready(){

    if(Np==1){interactions = false;}

    Q = arma::mat(3, Np).fill(0.); M = arma::mat(3, Np).fill(0.);  
    RU = arma::zeros(6, Np); dRU = arma::zeros(6, Np); 

    E_ext = arma::zeros(3, Np); B_ext = arma::zeros(3, Np); E_int = arma::zeros(3, Np);
    F_ext = arma::zeros(3, Np); F_int = arma::zeros(3, Np);

    dist = arma::mat(3, Np);
    norm3 = arma::mat(3, Np);
    K = arma::zeros(6, Np);
    Rnorm = arma::mat(3, Np);

    isready = true;

}

arma::mat PenningTrap::external_forces(arma::mat RU){
    compute_external_Efield(RU.rows(0,2));
    compute_external_Bfield(RU.rows(0,2));
   
    arma::mat U(3,Np); U.rows(0,2) = RU.rows(3,5);
    arma::mat UxB(3,Np);    //  cross products
    UxB.row(0) = U.row(1)%B_ext.row(2) - U.row(2)%B_ext.row(1);
    UxB.row(1) = U.row(2)%B_ext.row(0) - U.row(0)%B_ext.row(2);
    UxB.row(2) = U.row(0)%B_ext.row(1) - U.row(1)%B_ext.row(0);     //  always zero...
    F_ext = Q % (E_ext + UxB);
    return F_ext;
}

arma::mat PenningTrap::external_forces(double t, arma::mat RU){
    compute_external_Efield(t, RU.rows(0,2));
    compute_external_Bfield(RU.rows(0,2));
   
    arma::mat U(3,Np); U.rows(0,2) = RU.rows(3,5);
    arma::mat UxB(3,Np);    //  cross products
    UxB.row(0) = U.row(1)%B_ext.row(2) - U.row(2)%B_ext.row(1);
    UxB.row(1) = U.row(2)%B_ext.row(0) - U.row(0)%B_ext.row(2);
    UxB.row(2) = U.row(0)%B_ext.row(1) - U.row(1)%B_ext.row(0);     //  always zero...
    F_ext = Q % (E_ext + UxB);
    return F_ext;
}

arma::mat PenningTrap::internal_forces(arma::mat RU){
    compute_interaction_field(RU.rows(0,2));
    F_int = Q % E_int;
    return F_int;
}

arma::mat PenningTrap::evolve_FE(double dt, arma::mat RU){
    external_forces(RU);
    if(interactions){
        internal_forces(RU);
    }
    dRU.rows(3,5) = (F_ext + F_int) / M * dt;
    dRU.rows(0,2) = (RU.rows(3,5)) * dt;
    return dRU;
}

arma::mat PenningTrap::evolve_FE(double dt, double t, arma::mat RU){
    external_forces(t, RU);
    if(interactions){
        internal_forces(RU);
    }
    dRU.rows(3,5) = (F_ext + F_int) / M * dt;
    dRU.rows(0,2) = (RU.rows(3,5)) * dt;
    return dRU;
}

arma::mat PenningTrap::evolve_RK4(double dt, arma::mat RU){
    K1 = K_val(RU) * dt;
    K2 = K_val(RU+K1/2) * dt;
    K3 = K_val(RU+K2/2) * dt;
    K4 = K_val(RU+K3) * dt; 
    // std::cout<<"dt: "<<dt<<std::endl;
    // std::cout<<"K1: "<<K1<<std::endl;
    // std::cout<<"K2: "<<K2<<std::endl;
    // std::cout<<"K3: "<<K3<<std::endl;
    // std::cout<<"K4: "<<K4<<std::endl;
    dRU = (K1+2*K2+2*K3+K4)/6;
    return dRU;
}

arma::mat PenningTrap::evolve_RK4(double dt, double t, arma::mat RU){
    K1 = K_val(t, RU) * dt;
    K2 = K_val(t+dt/2, RU+K1/2) * dt;
    K3 = K_val(t+dt/2, RU+K2/2) * dt;
    K4 = K_val(t+dt, RU+K3) * dt; 

    dRU = (K1+2*K2+2*K3+K4)/6;
    return dRU;
}

arma::mat PenningTrap::compute_external_Efield(arma::mat R){
    E_ext = Vd2r * R;
    E_ext.row(2) *= -2;
    E_ext.elem(arma::find(Pnorm(R) > d)).zeros();   //  zero outside of trap
    return E_ext;
}

arma::mat PenningTrap::compute_external_Efield(double t, arma::mat R){
    E_ext = compute_external_Efield(R);
    E_ext *= (1 + f*cos(omega_V*t));
    return E_ext;
}

arma::mat PenningTrap::compute_external_Bfield(arma::mat R){
    B_ext.row(2).fill(B0);
    B_ext.elem(arma::find(Pnorm(R) > d)).zeros();   //  zero outside of trap 
    return B_ext;
}

arma::mat PenningTrap::compute_interaction_field(arma::mat R){
    for(int p=0; p<Np; p++){
        norm3.fill(0.);
        dist = R; dist.each_col() -= R.col(p);
        norm3 = arma::pow(Pnorm(dist), 3);
        norm3.col(p).fill(1.);                                  //  avoid dividing by zero
        E_int.col(p) = k_e * arma::sum(Q%dist/norm3, 1);        //  sum over all particles
    }
    return E_int;
}

void PenningTrap::generate_particle(std::vector<Particle*> &list, arma::vec r, arma::vec v){
    Particle* new_particle = new Particle(q_Ca, m_Ca, r, v);
    list.push_back(new_particle);  // maybe try calling 'add_particle(new_particle)'? or maybe not ...
    Np++;
}

arma::mat PenningTrap::K_val(arma::mat RU_){
    K.rows(0,2) = RU_.rows(3,5);
    external_forces(RU_);
    if(interactions){
        internal_forces(RU_);
    }
    // std::cout<<RU_<<std::endl;
    K.rows(3,5) = (F_ext + F_int) / M;
    return K;
}

arma::mat PenningTrap::K_val(double t, arma::mat RU_){
    K.rows(0,2) = RU_.rows(3,5);
    external_forces(t, RU_);
    if(interactions){
        internal_forces(RU_);
    }
    K.rows(3,5) = (F_ext + F_int) / M;
    return K;
}

arma::mat PenningTrap::Pnorm(arma::mat R){
    Rnorm.fill(0.);
    Rnorm.each_row() += arma::sqrt(arma::sum(arma::square(R), 0));
    return Rnorm; 
}

