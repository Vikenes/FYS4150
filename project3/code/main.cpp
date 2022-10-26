#include <algorithm>
#include "utils.hpp"
#include "algorithms.hpp"
#include "Particle.hpp"
#include "PenningTrap.hpp"

/**
 * PROJECT 3 FYS4150
 *      By: Johan Mylius Kroken, Vetle Vikenes and Nanna Bryne
 */

PenningTrap test_single_particle(std::string scheme, double no_of_timesteps, double simulation_duration=50){
    //  initialise the test particle:
    Particle p1 = Particle(q_Ca, m_Ca, arma::vec({20,0,20}), arma::vec({0,25,0}));
    //  initialise the trap:
    PenningTrap Trap = PenningTrap(B0, V0, d);
    Trap.add_particle(p1);
    //  run simulation:
    double dt = simulation_duration/no_of_timesteps;
    std::cout << "n: "<<no_of_timesteps<<"  "<<"dt: "<<dt<< std::endl;
    Trap.simulate(simulation_duration, dt, scheme);
    return Trap;
}

PenningTrap test_double_particle(bool with_interactions, std::string scheme, double no_of_timesteps, double simulation_duration=50){
    //  initialise the two test particles:
    Particle p1 = Particle(q_Ca, m_Ca, arma::vec({20,0,20}), arma::vec({0,25,0}));
    Particle p2 = Particle(q_Ca, m_Ca, arma::vec({25,25,0}), arma::vec({0,40,5}));
    //  initialise the trap:
    PenningTrap Trap = PenningTrap(B0, V0, d, with_interactions);
    Trap.add_particle(p1);
    Trap.add_particle(p2);
    //  run simulation:
    Trap.simulate(simulation_duration, simulation_duration/no_of_timesteps, scheme);
    return Trap;
}

int run_tests(std::string scheme){
    assert(scheme=="RK4" or scheme=="FE");
    std::string folder = scheme + "/";
    /**
     * Function for running the baby case simulations we have set up
     */
    double sim_duration = 50;
    arma::vec n = arma::vec({4000, 8000, 16000, 32000});
    double h; 

    /*  (1) Single particle 
            (a) for n =  4000
            (b) for n =  8000
            (c) for n = 16000
            (d) for n = 32000   */

    PenningTrap Trap1a = test_single_particle(scheme, n(0));
    Trap1a.save_solution(folder + "single_n1");
    
    PenningTrap Trap1b = test_single_particle(scheme, n(1));
    Trap1b.save_solution(folder + "single_n2");

    PenningTrap Trap1c = test_single_particle(scheme, n(2));
    Trap1c.save_solution(folder + "single_n3");

    PenningTrap Trap1d = test_single_particle(scheme, n(3));
    Trap1d.save_solution(folder + "single_n4");


    /*  (2) Double particle 
            (a) without interactions
            (b) with interactions   */
    
    PenningTrap Trap2a = test_double_particle(false, scheme, n(2));
    Trap2a.save_solution(folder + "double_without");

    PenningTrap Trap2b = test_double_particle(true, scheme, n(2));
    Trap2b.save_solution(folder + "double_with");

    return 0;
}

int particles_left(int sim_dur, double h, double amplitude, double frequency, std::string scheme="RK4", std::string interaction_switch="off"){
    PenningTrap Trap = PenningTrap(B0, V0, d);
    Trap.generate_random_identical_particles(q_Ca, m_Ca, 100);
    Trap.switch_interactions(interaction_switch);
    Trap.apply_time_dependence(amplitude, frequency);
    Trap.simulate(sim_dur, h, scheme);

    return Trap.count_particles();
}

int particles_left(double amplitude, arma::vec frequency, std::string filename, std::string scheme="RK4", std::string interaction_switch="off", int time_steps=80000){
    double sim_duration = 500;
    double h = sim_duration/time_steps;

    int Nomega = frequency.size();
    std::vector<int> trapped(Nomega);
    std::vector<double> omega(Nomega);
    for(int j=0; j<Nomega; j++){
        std::cout << "(" << j+1 << "/" << Nomega << ")" << std::endl;
        omega[j] = frequency(j);
        trapped[j] = particles_left(sim_duration, h, amplitude, frequency(j), scheme, interaction_switch);
        //std::cout << trapped[j] << std::endl;
    }
    std::string folder = scheme + "/";
    write_to_file(omega, trapped, folder+filename);
    
    return 0;
}

int main(){   

    auto start_time = std::chrono::high_resolution_clock::now();

    /* PROBLEM 8 */
    run_tests("FE");
    run_tests("RK4");

    /* PROBLEM 9 */
    double f1=0.1, f2=0.4, f3=0.7; // amplitudes

    //  performoing broad-band scan without interactions:
    arma::vec omega_V = arma::linspace(0.2, 2.5, 300); // [ MHz ] 
    double dOmega_coarse = omega_V(1)-omega_V(0);
    assert(dOmega_coarse<0.02);

    particles_left(f1, omega_V, "trapped_f1_without");
    particles_left(f2, omega_V, "trapped_f2_without");
    particles_left(f3, omega_V, "trapped_f3_without"); 
   
    //  performoing narrow-band scan with and without interactions for one amplitude
    arma::vec omega_V_fine = arma::linspace(1.35, 1.45, 50);
    double dOmega_fine = omega_V_fine(1)-omega_V_fine(0);
    assert(dOmega_fine < 0.5*dOmega_coarse);        //  at least twice as fine-grained
    particles_left(f1, omega_V_fine, "trapped_f1_with_fine", "RK4", "on", 20000);  // Running time: 6623.53 s (110.392 min) (O3)
    particles_left(f1, omega_V_fine, "trapped_f1_without_fine", "RK4", "off", 20000); // Running time: 83.3052 s (1.38842 min) (O3)
    
    auto stop_time = std::chrono::high_resolution_clock::now();

    double run_time = std::chrono::duration<double>(stop_time-start_time).count();

    std::cout << "Running time: " << run_time << " s (" << run_time/60.0 << " min)" << std::endl;

    return 0;
}