#include <algorithm>
#include "utils.hpp"
#include "algorithms.hpp"
#include "Particle.hpp"
#include "PenningTrap.hpp"

/**
 * PROJECT 3 FYS4150
 *      By: Johan Mylius Kroken, Vetle Vikenes and Nanna Bryne
 */




int test_single_particle(std::string scheme){

    return 0;
}

int test_double_particle(std::string scheme){

    return 0;
}



int run_tests(std::string scheme){
    assert(scheme=="RK4" or scheme=="FE");
    std::string folder = "tests/" + scheme + "/";
    /**
     * Function for running all tests we have created
     */
    PenningTrap Trap = PenningTrap(B0, V0, d);
    //  initialise the two test particles:
    Particle p1 = Particle(q_Ca, m_Ca, arma::vec({20,0,20}), arma::vec({0,25,0}));
    Particle p2 = Particle(q_Ca, m_Ca, arma::vec({25,25,0}), arma::vec({0,40,5}));
    double sim_duration = 50;
    arma::vec n = arma::vec({4000, 8000, 16000, 32000});
    double h; 

    /*  (1) Single particle 
            (a) for n =  4000
            (b) for n =  8000
            (c) for n = 16000
            (d) for n = 32000   */

    Particle p1a = p1;
    PenningTrap Trap1a = Trap;
    Trap1a.add_particle(p1a);
    Trap1a.simulate(sim_duration, sim_duration/n(0), scheme);
    Trap1a.save_solution(folder + "single_n1");
    
    Particle p1b = p1;
    PenningTrap Trap1b = Trap;
    Trap1b.add_particle(p1b);
    Trap1b.simulate(sim_duration, sim_duration/n(1), scheme);
    Trap1b.save_solution(folder + "single_n2");

    Particle p1c = p1;
    PenningTrap Trap1c = Trap;
    Trap1c.add_particle(p1c);
    Trap1c.simulate(sim_duration, sim_duration/n(2), scheme);
    Trap1c.save_solution(folder + "single_n3");

    Particle p1d = p1;
    PenningTrap Trap1d = Trap;
    Trap1d.add_particle(p1d);
    Trap1d.simulate(sim_duration, sim_duration/n(3), scheme);
    Trap1d.save_solution(folder + "single_n4");


    /*  (2) Double particle 
            (a) without interactions
            (b) with interactions   */

    h = sim_duration/n(2);
    
    PenningTrap Trap2a = PenningTrap(B0, V0, d, false);
    Particle p1a2 = p1; Particle p2a2 = p2;
    Trap2a.add_particle(p1a2);
    Trap2a.add_particle(p2a2);
    Trap2a.simulate(sim_duration, h, scheme);
    Trap2a.save_solution(folder + "double_without");

    PenningTrap Trap2b = PenningTrap(B0, V0, d, true);
    Particle p1b2 = p1; Particle p2b2 = p2;
    Trap2b.add_particle(p1b2);
    Trap2b.add_particle(p2b2);
    Trap2b.simulate(sim_duration, h, scheme);
    Trap2b.save_solution(folder + "double_with");

    return 0;
}



int generate_particle(std::vector<Particle*> &list, arma::vec r, arma::vec v){
    Particle* new_particle = new Particle(q_Ca, m_Ca, r, v);
    list.push_back(new_particle);

    return 0;
}


int time_dependent_potential(double amplitude, double frequency, std::string scheme="RK4"){
    // or some other name ....  
    assert(scheme=="RK4" or scheme=="FE");
    std::string folder = scheme + "/";
    
    double sim_duration = 500;
    double h = sim_duration/8000;

    PenningTrap Trap = PenningTrap(B0, V0, d); 
    Trap.switch_interactions("off");
    
    Trap.apply_time_dependence(amplitude, frequency);
    std::vector<Particle*> list;
    arma::arma_rng::set_seed(69); // OBS! want to have this in utils
    arma::vec rr;
    arma::vec vv;

    for(int p=0; p<100; p++){
        rr = arma::vec(3).randn() * 0.1 * d;                //  random initial position
        vv = arma::vec(3).randn() * 0.1 * d;                //  random initial velocity 
        generate_particle(list, rr, vv);
    }

    
    for(int i=0; i<10; i++){
        std::cout << list.at(i) -> position() << std::endl;
    }
    

    //Trap.generate_random_identical_particles(q_Ca, m_Ca, 100);
    Trap.simulate(sim_duration, h, scheme);
    Trap.save_solution(folder + "first");

    return 0;
}


int particles_left(PenningTrap trap, int sim_dur, double amplitude, double frequency, std::string scheme="RK4"){
    PenningTrap Trap = trap; // copy
    //Trap.switch_interactions("off");
    Trap.apply_time_dependence(amplitude, frequency);
    std::cout << Trap.Np << std::endl;
    Trap.simulate(sim_dur, sim_dur/1000, scheme);
    std::cout << Trap.Np << std::endl;
    return Trap.count_particles();
}

int particles_left(double amplitude, arma::vec frequency, std::string scheme="RK4"){
    double sim_duration = 500;
    double h = sim_duration/1000;

    PenningTrap Trap = PenningTrap(B0, V0, d);
    //Trap.generate_random_identical_particles(q_Ca, m_Ca, 10);
    Trap.switch_interactions("off");
    std::cout << Trap.Np << std::endl;


    int Nomega = frequency.size();
    std::vector<int> trapped(Nomega);
    for(int j=0; j<Nomega; j++){
        //trapped[j] = particles_left(Trap, sim_duration, amplitude, frequency(j));
        std::cout << particles_left(Trap, sim_duration, amplitude, frequency(j)) << std::endl;
    }
    
    return 0;
}



int main(){   

    auto start_time = std::chrono::high_resolution_clock::now();


    run_tests("FE");
    run_tests("RK4");



    arma::vec f = arma::vec({0.1,0.4,0.7});
    arma::vec omega_V = arma::linspace(0.2, 2.5, 100); // [ MHz ] 

    //time_dependent_potential(f(2), omega_V(70));
    //particles_left(f(0), omega_V);

    auto stop_time = std::chrono::high_resolution_clock::now();

    double run_time = std::chrono::duration<double>(stop_time-start_time).count();

    std::cout << "Running time: " << run_time << " s" << std::endl;


    return 0;
}