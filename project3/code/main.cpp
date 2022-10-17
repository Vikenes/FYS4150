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
    Trap1a.set_solution_filename(folder + "single_n1");
    Trap1a.add_particle(p1a);
    Trap1a.simulate(sim_duration, sim_duration/n(0), scheme);
    
    Particle p1b = p1;
    PenningTrap Trap1b = Trap;
    Trap1b.set_solution_filename(folder + "single_n2");
    Trap1b.add_particle(p1b);
    Trap1b.simulate(sim_duration, sim_duration/n(1), scheme);

    Particle p1c = p1;
    PenningTrap Trap1c = Trap;
    Trap1c.set_solution_filename(folder + "single_n3");
    Trap1c.add_particle(p1c);
    Trap1c.simulate(sim_duration, sim_duration/n(2), scheme);

    Particle p1d = p1;
    PenningTrap Trap1d = Trap;
    Trap1d.set_solution_filename(folder + "single_n4");
    Trap1d.add_particle(p1d);
    Trap1d.simulate(sim_duration, sim_duration/n(3), scheme);


    /*  (2) Double particle 
            (a) without interactions
            (b) with interactions   */
    h = sim_duration/n(2);
    
    PenningTrap Trap2a = PenningTrap(B0, V0, d, false);
    Particle p1a2 = p1; Particle p2a2 = p2;
    Trap2a.set_solution_filename(folder + "double_without");
    Trap2a.add_particle(p1a2);
    Trap2a.add_particle(p2a2);
    Trap2a.simulate(sim_duration, h, scheme);

    PenningTrap Trap2b = PenningTrap(B0, V0, d, true);
    Particle p1b2 = p1; Particle p2b2 = p2;
    Trap2b.set_solution_filename(folder + "double_with");
    Trap2b.add_particle(p1b2);
    Trap2b.add_particle(p2b2);
    Trap2b.simulate(sim_duration, h, scheme);

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
    Trap.set_solution_filename(folder + "first");
    Trap.apply_time_dependence(amplitude, frequency);
    Trap.generate_random_identical_particles(q_Ca, m_Ca, 10);
    Trap.simulate(sim_duration, h, scheme);

    return 0;
}





int main(){   

    auto start_time = std::chrono::high_resolution_clock::now();


    //run_tests("FE");
    //run_tests("RK4");

    arma::vec f = arma::vec({0.1,0.4,0.7});
    arma::vec omega_V = arma::linspace(0.2, 2.5, 100); // [ MHz ] 

    time_dependent_potential(f(2), omega_V(70));

    auto stop_time = std::chrono::high_resolution_clock::now();

    double run_time = std::chrono::duration<double>(stop_time-start_time).count();

    std::cout << "Running time: " << run_time << " s" << std::endl;


    return 0;
}